
#include "AlgebraicWM.hpp"

using namespace AlgebraicWM;

double calc_Tstar(const double re_delta, const double ystar) { 

  const double C_t = 3.9;
  const double kappa_t = 0.85/0.41;

  double Tstar; 

  if ( re_delta <= ystar) { 
    
    Tstar = 0.5 * re_delta;
    
  } else { 
    
    Tstar  = 0.5*ystar*ystar/re_delta;
    Tstar += C_t*(1.0  - ystar/re_delta);
    
    double tmp = log(re_delta) - 1.0;
    tmp       -= ystar/re_delta*(log(ystar) - 1.0);
    tmp       /= kappa_t;
    
    Tstar += tmp;
    
  }
  
  assert( Tstar >= 0.0);
  
  return Tstar;
}

void checkTstar() { 

  const double ystar = 8.3;

  for (int i = -2; i < 4; ++i) { 
    
    const double re_delta = double(pow(10,i));
    cout << " re_delta, tstar " << re_delta << "    " << calc_Tstar(re_delta,ystar) << endl;

  }

}

int main(const int argc, const char* argv[]) { 

  checkTstar();

  const double u_tau_exact = 0.123;
  const double nu          = 1.0e-05;
  const double rho         = 1.2333;
  const double mu          = rho*nu;

  const double yplus[9]    = {2.0, 15.0, 20.0, 23.32, 40.0, 100.0, 1000.0, 10000.0};
  double max_err           = 0.0;

  for (int i =0; i < 8; ++i) { 
    const double y1 = nu*yplus[i]/u_tau_exact;
    const double u1 = u_tau_exact*uplus_fit(yplus[i]);

    const double tau_wall = solve_tau(u1,y1,rho,mu);
    const double u_tau    = sqrt(tau_wall/rho);
    const double err      = abs(u_tau - u_tau_exact); 
    max_err               = max(max_err,err);

    cout << " yplus = " << yplus[i] << ".....  u_tau (should be " << u_tau_exact << " ) : " << u_tau 
         << " ; error = " << err << endl;
  }

  { // separation point test... 
    const double y1 = 0.005;
    const double u1 = 0.0;
    const double tau_wall = solve_tau(u1,y1,rho,mu);
    const double u_tau    = sqrt(tau_wall/rho);
    const double err      = abs(u_tau);
    max_err = max(max_err,err);
    cout << " sep point check (should be zero): " << u_tau << endl;
  }

  if ( max_err < 1.0e-07) { 
    cout << " looks okay. " << endl;
    return 0;
  } else { 
    cout << " ... something is awry. " << endl;
    return 1;
  }
}
