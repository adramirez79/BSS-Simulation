#Coded by Hector Gómez Márquez speed_nation@hotmail.com
#Last revision: Nov 17th, 2020

library(Rcpp)
library(tidyverse)
library(stringr)
library(logspline)
library(modelr)
library(rstream)
library(fitdistrplus)

sourceCpp("simulation.cpp")

#anixx - dataframe
#anixx$Estacion_Retiro - column with the stations where withdrawals occurred
#anixx$Estacion_Arribo - column with the stations where returns occurred
#anixx$Tiempo - column with trip times
#anixx$abs_salida - column with the times of arrivals since opening
#dt - logspline demand density models
#anixx$abs_regreso - column with the times of returns since opening
#ot - logspline offer density models
#yes - a 450 x 100 matrix where each row pertains a bike station, and each observation is the value of the density function for the horizon divided by 100 pieces.

lnormales <- anixx %>%
  split(.$Estacion_Retiro) %>%
  map(~fitdist(data=.$Tiempo,distr = "lnorm"))

meanlog <- numeric(450)
for(i in 1:length(lnormales)){
  meanlog[i] <- as.numeric(lnormales[[i]]$estimate[1])
}

sdlog <- numeric(450)
for(i in 1:length(lnormales)){
  sdlog[i] <- as.numeric(lnormales[[i]]$estimate[2])
}

dt <- anixx %>%
  split(.$Estacion_Retiro) %>%
  map(~logspline(x=.$abs_salida,lbound = 0,ubound = 68400,
                 error.action = 0))

ot <- anixx %>%
  split(.$Estacion_Arribo) %>%
  map(~logspline(x=.$abs_regreso,lbound = 0,ubound = 86400,
                 error.action = 2))

retiros_anixx <- anixx %>%
  group_by(Estacion_Retiro) %>%
  count() %>%
  arrange(desc(n))

prom_ret <- anixx %>%
  group_by(Estacion_Retiro) %>%
  count() %>%
  transmute(perday = n/21)  #weekly days sampled, example: 21 days

rates_r <- list(prom_ret$perday) #demand rates

regresos_anixx <- anixx %>%
  group_by(Estacion_Arribo) %>%
  count() %>%
  arrange(desc(n))

prom_col <- anixx %>%
  group_by(Estacion_Arribo) %>%
  count() %>%
  transmute(perday = n/21) #weekly days sampled

priors <- prom_col$perday  #offer rates

# TRIP GENERATION
{
  poissones <- rates_r %>%  #number of trips per station
    map(rpois, n = 450) %>%
    unlist()
  
  
  for (i in 1:450) {   # Trips generated
    viajes[[i]] <- rlogspline(poissones[i],fit = dt[[i]])
  }
  
  for (i in 1:450) { # adding traveling times
    finaliza[[i]] <- rlnorm(meanlog = meanlog[i],sdlog = sdlog[i],
                            n = poissones[i]) + viajes[[i]]
  }
  
  for (i in 1:450) {      #create their destination
    if(length(viajes[[i]]) != 0){
      for (j in 1:length(viajes[[i]])) {
        destinations[[i]][j] <- teleport2(finaliza[[i]][j],priors = priors,
                                          x = x,yes = yes)
      }
    }
    else {
      destinations[[i]] <- finaliza[[i]]
    }
  }
}

bikes <- round(rep(3000,450))  # vector of bikes per station
simular(solution = bikes,jumps = 0,actions = decisions,
             quita = NULL,nquita = est_quitar,programa = F,
             infcap = TRUE,revise = 27,revise_ofer = 27,est_exito = 27)


show_failstart() #return a vector with failed starts per station
show_failend()  #return a vector with failed ends per station


#FUNCTION SIMULAR------------------------------------------------
simular <- function(solution,jumps,actions,quita,nquita,
                         programa,revise,revise_ofer,est_exito,infcap) {
  
  #CREAR SIMULACION------
  destroy_queue()
  restart_fo()
  
  add_to_queue(0, s = 0, num = 0,y = 0, z = 0) # adding start event
  add_to_queue(86400,s = 3, num = 0, y = 0, z = 0) # adding finish event
  
  for (i in 1:450) {
    add_to_queue(x = viajes[[i]], s = 1, num = i, y = finaliza[[i]],
                 z = destinations[[i]])
  }
  
  if (programa == TRUE) { 
    for (i in 1:divisiones) {
      add_to_queue(x = schedule[i],s = 4,num = i,y = 0,z = 0)
    }
  }
  
  if (length(quita) != 0 ) {
    for (i in 1:length(quita)) {
      add_to_queue(x = quita[i],s = 6, num = nquita[i],y=0,z=0)
    }
  }
  
  simulation(bicicletas = solution,priors = priors,yes = yes,
             x = x,means = meanlog, desvs = sdlog, jumps = jumps, actions = actions,
             infcap = infcap,pwrup = TRUE, threshold = 0,
             revise = revise,revise_ofer = revise_ofer,est_exito = est_exito)
  
}

