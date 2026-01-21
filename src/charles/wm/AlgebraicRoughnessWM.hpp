#ifndef ALGEBRAICROUGHNESSWM_HPP
#define ALGEBRAICROUGHNESSWM_HPP

#include <iostream>
#include <cmath>
#include <assert.h>

using namespace std;

/// Simple logarithmic profile for a wall
/// surface roughness height z0
///
/// u(z) = u_star/kappa*ln( (z+z0)/z0 )
///
/// where kappa = 0.41

namespace AlgebraicRoughnessWM { 

  static const double kappa = 0.41;
  
  inline double uplus_fit(const double z, const double z0) { 
    return 1.0/kappa*log((z+z0)/z0);
  }

  inline int solve_ustar(double& u_star, const double u1, const double z1, 
                         const double z0, const double nu, const double u_star_guess = -1.0) { 
    
    assert( (u1 >= 0.0) && (z1 > 0.0));

    // set the initial guess 

    if ( u_star_guess <= 0.0 ) { 
      u_star              = sqrt(nu*u1/z1); 
    } else { 
      u_star             = u_star_guess;
    }

    const double relax   = 0.4;
    int iter             = 0;
    const int max_iter   = 1000;
    const double rel_tol = 1.0e-06;
    
    while ( true ) { 

      ++iter;
      double f,dfdustar;
      f       = u1 - (u_star/kappa)*log((z1+z0)/z0);
      dfdustar = -(1.0/kappa)*log((z1+z0)/z0);

      // u_star = 0 case
      if ( (abs(f) < 1.0e-14) && (abs(dfdustar) < 1.0e-14)) { 
        return 0;
      }
     
      if ( (abs(dfdustar) <= 1.0e-14) && (iter < max_iter)) {  
        cout << " problem converging -- reseting initial guess: " << iter << "   " << f << "   " << dfdustar << "    " << u_star << "   " << u1 << endl;
        u_star              = sqrt(nu*u1/z1);
        continue;
      }

      assert( abs(dfdustar) > 1.0e-14);
      double dustar = relax*f/dfdustar;
      u_star -= dustar;

      // this is a relative tolerance for the residual function
      
      double u_rel = abs(u_star);
      if ( u1 > 0.0) { 
        u_rel = min(u1,u_rel);
      }
      
      if ( abs(f/u_rel) < rel_tol) { 
        break;
      } else if ( iter > max_iter) { 
        return 1;
      } else { 
        //cout << "    >>>  iter, f = " << iter << "   " << f << endl;
      }
    }

    return 0;
  }

  inline double solve_tau(const double u1, const double z1, const double z0, const double rho, 
                          const double mu,const double u_star_guess = -1.0) { 
    double u_star; 
    int ierr = solve_ustar(u_star,u1,z1,z0,mu/rho, u_star_guess);
    assert( ierr == 0); 
    return rho*u_star*u_star;
  }


};//namespace AlgebraicRoughnessWM

#endif
