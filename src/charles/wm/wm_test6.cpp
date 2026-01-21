

#include "AnsatzWm.hpp"

int main() { 

  /*
  const double u_tau = 0.123; 
  const double nu    = 1.0e-05;
  const double y1    = 0.01;

  const double yplus = y1*u_tau/nu;
  cout << " yplus : " << yplus << endl;
  
  const double u1    = u_tau*(1.0/0.41*log(yplus) + 5.2);

  AnsatzWm test;
  double tau = test.solve_tau(u1, 0.0, 0.0, nu, y1, 0.0, 1.0e16);
  cout << " final u_tau = " << sqrt(tau) << endl;
  cout << " ===================================== " << endl;
  
  tau = test.solve_tau(u1, -0.4*u1, -1.0e-01*u1*u1, nu, y1, 0.0, 1.0e16);
  cout << " fina u_tau = " << sqrt(abs(tau)) << endl;
  cout << " ===================================== " << endl;
  */

  { 
    const double u1   = 69.726;
    const double v1   = 30.0255;
    const double dpdx = -298.347;
    const double y1   = 0.000488717;
    const double nu   = 1.8734e-05;
    
    AnsatzWm test;
    double tau = test.solve_tau(u1,v1,dpdx,nu,y1, nu*u1/y1, 1.0e-6);
  }


  return 0;
}

