#ifndef ALGEBRAICWM_HPP
#define ALGEBRAICWM_HPP

#include <iostream>
#include <cmath>
#include <assert.h>

using namespace std;

/// we depart from the algebraic closures made by 
/// taylor and von-karman as they are not C1 continuous
/// and spalding's fit does not explicitly satisfy a 
/// logarithmic velocity profile
///
/// instead we construct the following piecewise 
/// law of the wall..
///
///          = y+ + a1*y+*y+  for y^+ < 23.32
/// u(y+) = 
///          = 1/kappa*ln(y+) + B for y^+ > 23.32
///
/// where kappa = 0.41, B = 5.2, a1 = ((1/kappa/y*) - 1)/2/y*
/// and y* \approx 23.32

namespace AlgebraicWM { 

  static const double ystar = 23.3244;
  static const double kappa = 0.41;
  static const double B     = 5.2;
  static const double a1    = 0.5*(1.0/(kappa*ystar) - 1.0)/ystar;
  
  inline double uplus_fit(const double yplus) { 
    if ( yplus <= ystar) { 
      return yplus + a1*yplus*yplus;
    } else { 
      return 1.0/kappa*log(yplus) + B;
    }
  }

  inline int solve_utau(double& u_tau, const double u1, const double y1, 
                        const double nu, const double u_tau_guess = -1.0) { 
    
    // recall  y^+ = u_tau*y/\nui
    //assert( (u1 >= 0.0) && (y1 > 0.0)); // allow NaN's

    // set the initial guess 

    if ( u_tau_guess <= 0.0 ) { 
      u_tau              = sqrt(nu*u1/y1); 
    } else { 
      u_tau              = u_tau_guess;
    }

    const double relax   = 0.4;
    int iter             = 0;
    const int max_iter   = 1000;
    const double rel_tol = 1.0e-06;
    
    while ( true ) { 

      ++iter;
      double f,dfdutau;
      double yplus     = u_tau*y1/nu; // existing yplus value.. 
      if ( yplus <= ystar) { 
        f       = u1 - u_tau*yplus - u_tau*a1*yplus*yplus;
        dfdutau = -2.0*yplus - 3.0*a1*yplus*yplus;
      } else { 
        f       = u1 - u_tau/kappa*log(yplus) - u_tau*B;
        dfdutau = -1.0/kappa*log(yplus) - 1.0/kappa - B;
      }

      // u_tau = 0 case
      if ( (abs(f) < 1.0e-14) && (abs(dfdutau) < 1.0e-14)) { 
        return 0;
      }
     
      if ( (abs(dfdutau) <= 1.0e-14) && (iter < max_iter)) {  
        cout << " problem converging -- reseting initial guess: " << iter << "   " << f << "   " << dfdutau << "    " << u_tau << "   " << u1 << endl;
        u_tau              = sqrt(nu*u1/y1);
        continue;
      }

      //assert( abs(dfdutau) > 1.0e-14);
      double dutau = relax*f/dfdutau;
      u_tau -= dutau;

      // this is a relative tolerance for the residual function
      
      double u_rel = abs(u_tau);
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

  inline double solve_tau(const double u1, const double y1, const double rho, 
                          const double mu,const double u_tau_guess = -1.0) { 
    double u_tau; 
    int ierr = solve_utau(u_tau,u1,y1,mu/rho, u_tau_guess);
    //assert( ierr == 0); 
    return rho*u_tau*u_tau;
  }
  
  
  inline double compute_q_wall_approx(const double T1, const double Tw, const double tau_w, 
                                      const double rho, const double mu, const double locp, 
                                      const double cp, const double y1) { 

    const double Pr_t     = 0.85;
    const double kappa_t  = kappa / Pr_t;
    const double C_t      = 3.9; // officially a function of Pr -- this is approx for air via Bradshaw
    const double ystar    = 8.28314; 

    // there is an approximate damping function and the damping of the eddy viscosity
    // was ignored in the derivation of this approximate q_wall

    const double u_tau = sqrt(tau_w/ rho);
    const double mu_t  = rho * kappa * u_tau * y1;
    const double l_t   = mu_t * cp / Pr_t;
    const double l_tot = locp*cp + l_t;

    const double rhs   = l_tot*(T1 - Tw)/y1;

    const double re_delta = u_tau * rho * y1 / mu;

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

    //assert( Tstar >= 0.0);

    const double lhs      = 1.0 + kappa_t*Tstar; 
    
    if ( (lhs <= 0.0) || (lhs != lhs)) 
      cout << Tstar << "    " << y1 << "    " << "   " << re_delta << "   " << Tstar << endl;
    //assert( lhs > 0.0);
    return rhs/lhs;
  
  }
  
  inline double compute_q_wall_approx_r(const double T1, const double Tw, const double r, const double tau_w, 
                                        const double rho, const double mu, const double locp, 
                                        const double cp, const double y1,const bool debug) {
    
    // This modified version includes a thermal resistance 
    // r is a thermal resistance like a "dx/k"...
    
    const double Pr_t     = 0.85;
    const double kappa_t  = kappa / Pr_t;
    const double C_t      = 3.9; // officially a function of Pr -- this is approx for air via Bradshaw
    const double ystar    = 8.28314; 

    // there is an approximate damping function and the damping of the eddy viscosity
    // was ignored in the derivation of this approximate q_wall

    const double u_tau = sqrt(tau_w/ rho);
    const double mu_t  = rho * kappa * u_tau * y1;
    const double l_t   = mu_t * cp / Pr_t;
    const double l_tot = locp*cp + l_t;

    const double rhs   = l_tot*(T1 - Tw)/y1;

    const double re_delta = u_tau * rho * y1 / mu;

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

    //assert( Tstar >= 0.0);

    const double lhs      = 1.0 + kappa_t*Tstar; 
    
    if ( (lhs <= 0.0) || (lhs != lhs)) 
      cout << Tstar << "    " << y1 << "    " << "   " << re_delta << "   " << Tstar << endl;
    
    //assert( lhs > 0.0);
    //return rhs/lhs;
    
    // for the case with thermal resistance "r"...
    // introduce the intermediate temperature Ti. Replace Tw with Ti in rhs above,
    // because it will be the interface temperture that the wall model sees, 
    // then add the equation q = (Ti-Tw)/r and solve for q (and Ti)...
    
    // for debugging...
    if (debug) {
      double Ti = (l_tot*r*T1 + Tw*lhs*y1)/(l_tot*r + y1*lhs);
      cout << "WALL MODEL got Ti = " << Ti << " for thermal resistance: " << r << " and q: " << rhs*y1 / (l_tot*r + y1*lhs) << endl;
    }
    
    return rhs*y1 / (l_tot*r + y1*lhs);
    
  }

  inline double compute_R_approx(const double tau_w,  const double rho,  const double mu, 
                                      const double locp,   const double cp,   const double y1) { 
    const double Pr_t     = 0.85;
    const double kappa_t  = kappa / Pr_t;
    const double C_t      = 3.9; // officially a function of Pr -- this is approx for air via Bradshaw
    const double ystar    = 8.28314; 
    
    // there is an approximate damping function and the damping of the eddy viscosity
    // was ignored in the derivation of this approximate q_wall
    const double u_tau = sqrt(tau_w/ rho);
    const double mu_t  = rho * kappa * u_tau * y1;
    const double l_t   = mu_t * cp / Pr_t;
    const double l_tot = locp*cp + l_t;
    
    //const double rhs   = l_tot*(T1 - Tw)/y1;
    const double rhs   = l_tot/y1; //NOTE THAT DELTA_T IS REMOVED HERE SO IS THERMAL RESISTANCE NOT Q!

    const double re_delta = u_tau * rho * y1 / mu;

    double Tstar; 

    if (re_delta <= ystar) { 
      Tstar = 0.5 * re_delta;
    }//(re_delta <= ystar) 

    else { 
      Tstar  = 0.5*ystar*ystar/re_delta;
      Tstar += C_t*(1.0  - ystar/re_delta);

      double tmp = log(re_delta) - 1.0;
      tmp       -= ystar/re_delta*(log(ystar) - 1.0);
      tmp       /= kappa_t;

      Tstar += tmp;
    }//(re_delta > ystar)

    //assert(Tstar >= 0.0);

    const double lhs      = 1.0 + kappa_t*Tstar; 
    
    if ((lhs <= 0.0) || (lhs != lhs)) 
      cout << Tstar << "    " << y1 << "    " << "   " << re_delta << "   " << Tstar << endl;
    //assert(lhs > 0.0);
    return(lhs/rhs);
  }//compute_R_wall_approx()


};

#endif
