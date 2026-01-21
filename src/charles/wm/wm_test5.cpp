

#include "NonParallelWm.hpp"

int main() { 

  const double u_tau = 0.123; 
  const double nu    = 1.0e-05;
  const double y1    = 0.01;

  const double yplus = y1*u_tau/nu;
  cout << " yplus : " << yplus << endl;
  
  const double u1    = u_tau*(1.0/0.41*log(yplus) + 5.3);

  NonParallelWm test;

  cout << " ==================================== " << endl;
  test.solve_tau(u1, 0.0, 0.0, nu, y1);
  cout << " ==================================== " << endl;
  test.solve_tau(u1, 0.4*u1, -1.0e-01*u1*u1, nu, y1);
  cout << " ==================================== " << endl;
  test.solve_tau(u1, -0.1*u1, 1.0e-01*u1*u1, nu, y1);
  cout << " ==================================== " << endl;
  test.solve_tau(u1, -0.4*u1, 1.0e-01*u1*u1, nu, y1);
  cout << " ==================================== " << endl;

  return 0;
}

