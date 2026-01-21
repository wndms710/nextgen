#ifndef NONPARALLELWM_HPP
#define NONPARALLELWM_HPP

#include <iostream>
#include <assert.h>
#include <cmath>

using namespace std;

class NonParallelWm { 
public:


  double damp(const double yplus) { 

    const double Aplus = 26.0;
    const double fax   = 1.0 - exp(-yplus/Aplus);
    return fax*fax;

  }

  void integrate(double &u, double& int1, int& n, 
                 const double tau, const double alpha, const double dpdx, 
                 const double nu, const double y1, const bool verbose = false) { 

    const double gf    = 1.02;

    const double u_tau  = sqrt(abs(tau));
    double dudy         = tau/nu; 
    u                   = 0.0; 
    double y            = 0.0;
    double dy         = min(0.5*nu/u_tau, 0.01*y1);
    //double dy           = y1/512.0;
    const double kappa  = 0.41; 
    
    double u_prev, y_prev;
    
    // other necessary integrls.. 
    
    double psi_y  = 0.0;    // int 0^y (u*y) dy
    int1          = 0.0;    // int 0^y dy/nu_tot

    n = 0;

    while ( y < y1) { 
      
      y_prev    = y; 
      u_prev    = u;
      
      // advance the spatial coordinate .. 
      
      y += min(y1-y_prev,dy);
      
      // set the turbulent viscosity .. 
      
      const double y_half = 0.5*(y + y_prev);
      const double ypls   = y_half*u_tau/nu;
      const double nu_t   = kappa*kappa*y_half*y_half*abs(dudy)*damp(ypls);
      const double nu_tot = nu + nu_t;
     
      const double fax = dy/nu_tot;
      const double vv  = alpha*y_half*y_half;
      const double r   = u + fax*(tau + dpdx*y_half - 2.0*alpha*psi_y + 0.5*u*vv);
      u                = r/(1.0 - 0.5*fax*vv);
 
      int1            += dy/nu_tot;
      psi_y           += dy*0.5*(u_prev + u)*y_half;
      dudy             = (u - u_prev)/(y - y_prev);

      dy                  = min(0.01*y1,gf*dy);
      ++n;

      if ( verbose) 
        cout << " XXXX " << y << "    " << u << endl;

    }

    // ensure that we reached the proper termination condition.. 
    
    assert( (y1 - y)/y1 < 1.0e-12); 
    
  }


  double func(const double tau, const double u1, const double alpha, const double dpdx, const double nu, 
              const double y1) { 

    double int1,u_outer;
    int n;

    integrate(u_outer,int1,n,tau,alpha,dpdx,nu,y1,false);
    return u_outer - u1;

  }
  
  double solve_tau(const double u1, const double v1, const double dpdx, const double nu, const double y1) { 


    const double tau_nom = nu*u1/y1;
    const double alpha   = v1/y1/y1;

    double tau_hi = 5.0*tau_nom;
    double tau_lo = -1.0*tau_nom;

    double fhi    = func(tau_hi,u1,alpha,dpdx,nu,y1);
    double flo    = func(tau_lo,u1,alpha,dpdx,nu,y1);

    if ( (fhi > 0.0) && (flo > 0.0)) { 

      int iter = 0;
      while ( flo > 0.0) { 
        ++iter;
        tau_lo *= 2.0;
        flo = func(tau_lo,u1,alpha,dpdx,nu,y1);

        if ( iter > 10)
          assert(0);
      }

    } else if ( (fhi < 0.0) && (flo < 0.0)) { 

      int iter = 0;
      while ( fhi < 0.0) { 
        ++iter;
        tau_hi *= 2.0;
        fhi = func(tau_hi,u1,alpha,dpdx,nu,y1);

        if ( iter > 10) 
          assert(0);
      }

    }


    assert( flo*fhi < 0.0);

    const int max_iter = 500;
    const double tol   = 1.0e-6;
    double tau_mid;

    for (int iter = 0; iter < max_iter; ++iter) { 

      tau_mid              = 0.5*(tau_hi + tau_lo);
      const double fmid    = func(tau_mid,u1,alpha,dpdx,nu,y1);

#ifdef DEBUG
      cout << "iter, res, tau: " << iter << "   " << fmid/u1 << "   " << tau_mid << endl;
#endif

      if ( abs(fmid) < tol*u1) { 

        break;

      } else { 

        if ( fmid*fhi > 0.0) { 

          tau_hi = tau_mid;

        } else { 

          assert( flo*fmid > 0.0);
          tau_lo = tau_mid;

        }
      }

    }

    /*
    double u_outer,int1;
    int n;
    integrate(u_outer,int1,n,tau_mid,alpha,dpdx,nu,y1,true);
    */

    return tau_mid;

  }

  
  double solve_tau_(const double u1, const double v1, const double dpdx, 
                                 const double nu, const double y1) { 

    
    // set the initial guess .. 

    double tau         = nu*u1/y1;
    const double alpha = v1/(y1*y1);

    const double tol      = 1.0e-06;
    const double max_iter = 500;
    const double relax    = 0.8;

    double u_outer, int1;
    int n;


    for (int iter = 0; iter < max_iter; ++iter) { 

      integrate(u_outer,int1,n,tau,alpha,dpdx,nu,y1);

      const double res = u1 - u_outer;
#ifdef DEBUG
      cout << " iter, res, tau, u, u1: " << iter << "   " << res/u1 << "    " << tau << "   " << u_outer << "   " 
           << u1 << "   " << n << endl;
#endif

      if ( abs(res)/u1 < tol) { 
        
        break;

      } else { 


        const double dtau = res/int1;
        tau              += relax*dtau;

      }

    }

      
    //integrate(u_outer,int1,n,tau,alpha,dpdx,nu,y1,true);

    return tau;

  }
};


#endif
