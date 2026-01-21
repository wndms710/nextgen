#ifndef EQUILIBRIUMWM_HPP
#define EQUILIBRIUMWM_HPP


#include <iostream>
#include <cmath>
#include <assert.h>

using namespace std;

namespace EquilibriumWM { 


  static const double kappa        = 0.41;
  static const double Pr_t         = 0.9;
  static const double mu_power_law = 0.76;
  //static const double step_ctl_fax = 0.05;
  static const double step_ctl_fax = 0.1;
  //static const double relax        = 0.7;

  inline double damp(const double yplus) { 

    const double inv_Aplus = 1.0/17.0;
    const double fax       = 1.0- exp(-yplus*inv_Aplus);
    return fax*fax;

  }


  // =======================================================================
  // solve for qwall, tau_wall 
  // the following assumes an ideal gas with constant heat capacities and a 
  // no-slip, non-moving wall.
  // =======================================================================

  inline int solve_qt(double& tau_w, double& q_w, const double u1, const double T1, const double rho1,
                      const double y1, const double mu1, const double Pr_lam, const double cp, 
                      double& T_w, const double relax, const double y_plus0, const double tol, const bool b_use_guess=false) { 


    double T_wall;

    if ( T_w < 0.0) { 

      // denotes an adiabatic case s.t. q_wall = 0.0
      // set an initial guess for the adiabatic wall temperature

      q_w = 0.0;

      if ( !b_use_guess) 
        T_wall = T1 + 0.5*u1*u1/cp; 
      else 
        T_wall = -T_w;

    } else { 

      // isothermal condition ... 

      T_wall = T_w;

    }

    // fast return for the case where the matching velocity itself is zero .. 

    if ( abs(u1) < 1.0e-12) { 

      tau_w = 0.0;
      q_w   = cp*mu1*(T1-T_wall)/(Pr_lam*y1);  //laminar clousre for heat flux
      return 0;

    }


    double u_tau, l_visc, mu_w;
    //const double tol = 1.0e-4;

       
    // set the initial guesses for tau_wall, q_wall

    if ( !b_use_guess) { 
      tau_w           = mu1*u1/y1;
      q_w             = mu1*(cp*T1 + 0.5*u1*u1 - cp*T_wall)/y1;
    }

    const int max_iter = 500;
    int iter;
    for (iter = 0; iter < max_iter; ++iter) { 

      // approximate constant pressure in the inner layer to compute a wall density.. 

      const double rho_w = rho1*T1/T_wall;
      mu_w               = mu1*pow(T_wall/T1, mu_power_law);

      u_tau              = sqrt(tau_w/rho_w);
      l_visc             = mu_w/(rho_w*u_tau);
      double y           = 0.0;
      double int1        = 0.0;
      double int2        = 0.0;
      double int3        = 0.0;
      double int4        = 0.0;
      double u           = 0.0;
      double T           = T_wall;
      double dy          = y_plus0*l_visc;
      int n_pts          = 0;
      
      while ( y < y1 ) { 

        // TODO adaptive choose the step size in y .. 

        dy                     = min(dy, y1-y);
        const double y_half    = y + 0.5*dy;
        const double A         = kappa*sqrt(rho_w)*sqrt(T_wall/T)*y_half*damp(y_half/l_visc);
        const double mu_t      = A*sqrt(tau_w);
        //const double mu_t      = kappa*rho_w*u_tau*sqrt(T_wall/T)*y_half*damp(y_half/l_visc);
        const double mu_lam    = mu1*pow(T/T1,mu_power_law);
        const double mu_tot    = mu_lam + mu_t;
        const double k_lam     = cp*mu_lam/Pr_lam;
        const double k_t       = cp*mu_t/Pr_t;
        const double k_tot     = k_lam + k_t;

        //cout << " ZZZZ " << n_pts << "    " << y_half << endl;

        // with the linearized approximations for the viscosities, this is an implicit 
        // midpoint update .. 

        const double unew      = u + dy*tau_w/mu_tot;
        const double Tnew      = T + dy/k_tot*(q_w - 0.5*tau_w*u - 0.5*tau_w*unew);

        // update the integrals .. 
        
        const double mu_prime  = Pr_t*mu_lam/Pr_lam;

        int1 += dy/(mu_prime + mu_t);
        int2 += 0.5*(u + unew)*dy/(mu_prime + mu_t);
        int3 += dy/(mu_lam + mu_t);
        int4 += 0.5*A*dy/(mu_tot*mu_tot);

        // update the solutions and y 

        const double dudy = (unew-u)/dy;
        const double dTdy = (Tnew-T)/dy;

        y    += dy;
        u     = unew;
        T     = Tnew;

#ifdef DEBUG
        cout << " XXXX   iter=" << iter << "   " << y << "    " << u << "   " << T << endl;
#endif

        // update the dy .. 

        dy = max(dy,min(step_ctl_fax*u/abs(dudy),step_ctl_fax*T/abs(dTdy)));
        ++n_pts;
      }


      // make sure we got to the right coordinate .. 

      assert( abs(y-y1)/y1 < 1.0e-12); 

      // check for convergence .. 

      if ( (abs(T-T1)/T1 < tol) && ( (abs(u1) < 1.0e-12) || ( abs(u-u1)/u1 < tol))) 
        break;

      // update the estimates for tau_w, q_wall .. 

      //tau_w              = u1/int3;

      const double J       = int3 - sqrt(tau_w)*int4; 
      const double dtau    = (u1-u)/J;
      tau_w               += relax*dtau;

      if ( T_w < 0.0) { 

        // adiabatic case.. need to update the guess for T_wall ... 

        q_w    = 0.0;
        T_wall = (1.0 - relax)*T_wall + relax* ( T1 + tau_w*Pr_t/cp*int2);


      } else { 

        // isothermal case .. 

        const double dhstar = cp*(T1-T_w);
        q_w = relax* ((dhstar/Pr_t + tau_w*int2)/int1) + (1.0-relax)*q_w;

      }

#ifdef DEBUG
      cout << " iter , tau_w, q_w, T_conv, u_conv, n_pts : " << iter 
           << "    " << tau_w << "    " << q_w << "    " << abs(T-T1)/T1 << "    " << abs(u-u1)/u1 
           << "    " << n_pts << endl;
#endif

    }

    if ( iter == max_iter) 
      return 1;

    if ( T_w < 0.0)
      T_w = -T_wall;

    return 0;
  }



};



#endif
