//Coded by Hector Gómez Márquez speed_nation@hotmail.com
//Last revision: Nov 17th, 2020

//[[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <tuple>
#include <queue>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>


using namespace Rcpp;

int INICIO_SIMULACION = 0;
int NUEVA_LLEGADA = 1;
int DESTINO = 2;
int FIN_SIMULACION = 3;
int ACTION = 4;
int COLOCAR = 5;
int QUITAR = 6;

auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
std::mt19937_64 gen(seed);
std::uniform_real_distribution<double> unif(0.0,1.0);


typedef std::tuple<double, int, int, double, int> triple;

int failstart[450] = {};
int failend[450] = {};
int estacion[450] = {0};
int estacion_out[450] = {0};
double first_time[450] = {0};
double first_out[450] = {0};

std::vector<double> salida(0);
std::vector<double> no_in(0);
std::vector<double> sucess(0);

class my_greater  {
public:
  bool operator() (const triple& arg1, const triple& arg2) const
  {
    return std::get<0>(arg1) >  std::get<0>(arg2);
    return false;
  }
};


// [[Rcpp::export]]
void restart_fo()
{
  for ( int i = 0; i < 450; ++i){
    failstart[i] = 0;
    failend[i] = 0;
    estacion[i] = 0;
    first_time[i] = 0;
    first_out[i] = 0;
    estacion_out[i] = 0;
    salida.clear(); 
    no_in.clear();
    sucess.clear();
  }
}

// [[Rcpp::export]]
NumericVector show_failstart()
{
  NumericVector out(450);
  for ( int i = 0; i < 450; ++i){
    out[i] = failstart[i];
  }
  return out;
}

// [[Rcpp::export]]
NumericVector show_estacion()
{
  NumericVector out(450);
  for ( int i = 0; i < 450; ++i){
    out[i] = estacion[i];
  }
  return out;
}

// [[Rcpp::export]]
NumericVector show_failend()
{
  NumericVector out(450);
  for ( int i = 0; i < 450; ++i){
    out[i] = failend[i];
  }
  
  return out;
}

// [[Rcpp::export]]
NumericVector show_first()
{
  NumericVector out(450);
  for ( int i = 0; i < 450; ++i){
    out[i] = first_time[i];
  }
  return out;
}


// [[Rcpp::export]]
NumericVector show_first_out()
{
  NumericVector out(450);
  for ( int i = 0; i < 450; ++i){
    out[i] = first_out[i];
  }
  return out;
}

// [[Rcpp::export]]
std::vector<double> show_return()
{
  return salida;
}

// [[Rcpp::export]]
std::vector<double> show_noin()
{
  return no_in;
}


// [[Rcpp::export]]
std::vector<double> show_sucess()
{
  return sucess;
}

std::priority_queue< triple, std::vector<triple>, my_greater> eventos;

// [[Rcpp::export]]
void add_to_queue(NumericVector x, int s, int num, NumericVector y, NumericVector z)
{
  for ( int i = 0; i < x.size(); ++i){
    triple suceso(x[i],s,num, y[i], z[i]);
    eventos.push(suceso);
  }
}

// [[Rcpp::export]]
void destroy_queue()
{
  eventos = std::priority_queue< triple, std::vector<triple>, my_greater>();
}

// [[Rcpp::export]]
NumericVector top_queue()
{
  NumericVector out(5);
  triple res = eventos.top();
  out[0] = std::get<0>(res);
  out[1] = std::get<1>(res);
  out[2] = std::get<2>(res);
  out[3] = std::get<3>(res);
  out[4] = std::get<4>(res);
  return out;
}

// [[Rcpp::export]]
void pop_queue()
{
  eventos.pop();
}


// [[Rcpp::export]]
double inter(double tiempo, int estacion, NumericMatrix& yes, NumericVector& x, int t)
{
  double xL = x[t], yL = yes(t,estacion), xR = x[t+1], yR = yes(t+1,estacion); 
  if ( tiempo < xL ) yR = yL; 
  if ( tiempo > xR ) yL = yR;
  
  double dydx = ( yR - yL ) / ( xR - xL );
  
  return yL + dydx * ( tiempo - xL ); 
}



// [[Rcpp::export]]
int teleport2(double tiempo,NumericVector& priors, NumericVector x,
              NumericMatrix& yes)
{
  double f_dens[450];
  
  int t;
  auto low1 = std::lower_bound(x.begin(),x.end(),tiempo);
  t = low1 - x.begin() - 1;
  
  double zeta = 0.0;
  for(int ni = 0; ni < 450; ++ni )
  {
    f_dens[ni] = inter(tiempo, ni,  yes, x, t);
    f_dens[ni] = priors[ni] * f_dens[ni];
    zeta += f_dens[ni];
  }
  
  double suma_parcial = 0.0;
  double alea =  unif(gen);
  
  for(int ni = 0; ni < 450; ++ni)
  {
    suma_parcial += f_dens[ni]/zeta;
    
    if(alea <= suma_parcial) return ni+1;
  }
  
  return 0;
}



// [[Rcpp::export]]
void simulation(NumericVector bicicletas, NumericVector& priors,
                NumericMatrix& yes, NumericVector x,NumericVector means,
                NumericVector desvs, NumericVector jumps,NumericMatrix& actions,
                NumericVector threshold ,
                bool infcap = false,bool pwrup = false,int revise = 1,
                int revise_ofer = 1, int est_exito = 1){
  
  int maxcap[] = {27, 12,36,15,12,15,24,12,24,36,14,14,36,36,24,36,27,36,23,
                  23,25,36,36,26,26,18,36,27,36,36,36,36,27,36,36,23,21,36,21,24,36,24,33,26,
                  36,35,24,12,36,26,27,21,15,35,15,36,12,27,34,36,35,18,36,36,18,26,23,30,36,
                  36,27,36,27,14,21,15,36,27,36,27,11,15,21,14,21,24,15,12,27,20,21,36,23,18,
                  33,27,36,30,24,33,15,27,15,24,25,35,36,29,36,24,23,21,24,27,24,33,21,36,24,
                  18,30,29,36,36,35,30,24,24,20,18,24,35,33,35,27,27,24,24,18,30,30,21,27,27,
                  24,27,24,36,27,27,30,21,33,33,18,30,36,36,30,21,27,36,26,20,24,24,19,21,21,
                  27,18,21,33,23,36,21,30,27,30,24,35,36,36,29,18,20,15,18,24,21,36,36,23,24,
                  27,36,27,23,25,18,21,24,27,20,24,24,27,24,27,36,36,27,24,24,24,36,24,36,24,
                  24,30,30,21,24,26,24,21,21,18,27,30,15,27,27,30,30,30,21,36,36,30,27,17,24,
                  14,30,27,18,27,18,18,21,30,36,25,21,28,24,24,29,33,29,36,36,36,36,35,36,36,
                  36,36,36,36,36,30,27,21,30,36,30,18,21,33,15,21,21,18,18,27,27,36,36,24,36,
                  27,32,21,24,21,23,35,27,33,36,36,27,30,27,15,18,27,24,20,21,27,36,21,27,18,
                  22,27,24,18,24,36,30,27,24,18,27,24,21,35,33,21,21,21,21,30,30,36,18,27,21,
                  24,27,18,21,30,18,30,27,24,20,18,24,24,33,26,27,20,21,27,24,21,23,35,30,21,
                  21,18,30,24,26,34,30,36,26,27,24,15,21,36,36,18,21,36,36,27,21,36,36,30,21,
                  30,26,36,30,33,30,29,30,30,21,21,27,24,24,36,36,36,30,30,30,30,36,29,27,30,
                  21,24,15,36,36,21,27,21,30,21,36,21,26,36,30,21,24,18,36,36,24,30,18,36,36,
                  27,30,36,30,36,18 };
  
  
  if (infcap == true){
    std::fill_n(maxcap,450,10000);
  }
  
  // eventos get0 = time,get1 = event type,get2 = generating station,
  //           get3 = ending,get4 = teleport
  
  for ( ; ; ){
    auto e = eventos.top();
    eventos.pop();
    if(std::get<1>(e) == FIN_SIMULACION){
      //Rcpp::Rcout << "Day finishes in EcoBici!!" << '\n';
      break;
    }
    
    if (std::get<1>(e) == INICIO_SIMULACION){
      //Rcpp::Rcout << "Day starts in EcoBici!!" << '\n';
    }
    
    if (std::get<1>(e) == ACTION){
      for (int i = 0; i < 450; ++i){
        bicicletas[i] =  actions(std::get<2>(e) -1,i);
      }
      
    }
    
    if(std::get<1>(e) == NUEVA_LLEGADA){
      
      if ( bicicletas[std::get<2>(e)-1] > 0){
        triple suceso( std::get<3>(e), DESTINO, 
                       std::get<4>(e),0,0 );
        eventos.push(suceso);
        bicicletas[std::get<2>(e)-1]--;
        if (std::get<2>(e) == est_exito ){
          sucess.push_back(std::get<0>(e));
        }
      }
      else {
        if (std::get<0>(e) > threshold[0] && std::get<0>(e) < threshold[1]){
          failstart[std::get<2>(e)-1]++;
        }
        if (std::get<2>(e) == revise ){
          no_in.push_back(std::get<0>(e));
        }
        if (estacion[std::get<2>(e)-1] == 0){
          first_time[std::get<2>(e)-1] = std::get<0>(e);
          estacion[std::get<2>(e)-1] = 1;
          if (jumps[estacion[std::get<2>(e)-1]] > 0){
            jumps[estacion[std::get<2>(e)-1]]--;
            estacion[std::get<2>(e)-1] = 0;
            first_time[std::get<2>(e)-1] = 0;
          }
        }
        
      }
      
    }
    
    if(std::get<1>(e) == DESTINO){
      
      if (bicicletas[std::get<2>(e)-1] >= maxcap[std::get<2>(e)-1]){
        
        if (std::get<0>(e) > threshold[0] && std::get<0>(e) < threshold[1]  ){
          
          failend[std::get<2>(e)-1]++;
          if (estacion_out[std::get<2>(e)-1] == 0){
            first_out[std::get<2>(e)-1] = std::get<0>(e);
            estacion_out[std::get<2>(e)-1] = 1;
          }
        }
      }
      else {
        bicicletas[std::get<2>(e)-1]++;
        if (std::get<2>(e) == revise_ofer ){
          salida.push_back(std::get<0>(e));
        }
        
      }
    }
    
    if (std::get<1>(e) == COLOCAR){
      if ( bicicletas[std::get<2>(e)-1] < maxcap[std::get<2>(e)-1] ){
        bicicletas[std::get<2>(e)-1] = bicicletas[std::get<2>(e)-1] + 
          1 * std::get<3>(e);
      }
      if ( bicicletas[std::get<2>(e)-1] > maxcap[std::get<2>(e)-1] ){
        bicicletas[std::get<2>(e)-1] = maxcap[std::get<2>(e)-1];
        
      }
    }
    
    if (std::get<1>(e) == QUITAR){
      if ( bicicletas[std::get<2>(e)-1] > 0){
        bicicletas[std::get<2>(e)-1] = 1;
      }
    }
    
  }
  
}
