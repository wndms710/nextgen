
#include <iostream>
#include <cmath>
#include "EquilibriumWM.hpp"
#include <fenv.h>
#include <cstdio>

using namespace EquilibriumWM;
using namespace std;

int main() { 

  //feenableexcept(FE_DIVBYZERO|FE_INVALID) ;


  const double gamma   = 1.4;
  const double p_ref   = 1.0/gamma;
  const double rho_ref = 1.0;
  const double T_ref   = 1.0;
  const double R_gas   = p_ref/rho_ref/T_ref;
  const double mu_ref  =  3e-06;
  
  const double cp = (p_ref/(T_ref*rho_ref)) * gamma / ( gamma - 1.0);
  const double Pr_lam  = 0.7;


  /*
  const double T_wall =  1.5;
  const double u1   = 2.55;
  const double T1   =   1.19155122;
  const double rho1 = p_ref/T1/R_gas;
  const double mu1  =  mu_ref*pow(T1/T_ref,0.76);
  const double y1   =  0.01704540;
  */

  cout << " ================================== " << endl;
  cout << " > solving the isothermal problem ... " << endl;


  const double T_wall = 4.5;
  const double u1     = 6.0;
  const double rho1   = 1.0;
  const double mu1    = 0.00087847999999999984;
  const double y1 = 0.064974367799974134;
  const double T1     =  0.99999999999999989;

  double tau_w, q_w;

  solve_qt(tau_w,q_w,u1,T1,rho1,y1,mu1,Pr_lam,cp,T_wall);

  cout << " > tau_w, q_wall : " << tau_w << "   " << q_w << endl;
  cout << " =================================== " << endl;
  //getchar();

  cout << " > solving the adiabatic problem ... " << endl;

  solve_qt(tau_w,q_w,u1,T1,rho1,y1,mu1,Pr_lam,cp,-1.0);

  cout << " > tau_w, q_wall (adiabatic) : " << tau_w << "    " << q_w << endl;
  cout << " =================================== " << endl;

  return 0;

} 
