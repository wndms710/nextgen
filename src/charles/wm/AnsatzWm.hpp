#ifndef ANSATZWM_HPP
#define ANSATZWM_HPP

#include <iostream>
#include <assert.h>
#include <cmath>
#include <cstdio>

using namespace std;

inline double SIGN_(double val) { 

  if ( val > 0.0) 
    return 1.0;
  else 
    return -1.0;

}


class AnsatzWm { 
public: 

  double ystar_pls;
  double kappa;
  double B;

  AnsatzWm() { 

    ystar_pls = 11.0;
    kappa     = 0.41;
    B         = ystar_pls - 1.0/kappa*log(ystar_pls);
  }

  double funcU(const double y, const double u_tau, const double tau, const double nu) { 

    const double a    = tau/nu;
    const double ypls = u_tau*y/nu;
    const double st   = SIGN_(tau);

    if ( ypls < ystar_pls) { 
    //if ( true) { 

      return a*y;

    } else { 

      return st*u_tau*(1.0/kappa*log(ypls) + B);

    }

  }

  double integralI(const double y, const double u_tau, const double tau, const double nu) { 

    const double a    = tau/nu;
    const double ypls = u_tau*y/nu;
    const double st   = SIGN_(tau);

    if ( ypls < ystar_pls) { 
    //if ( true) { 

      return 0.5*a*y*y;

    } else { 

      const double ystar = ystar_pls*nu/u_tau;
      double sum         = 0.5*a*ystar*ystar;
      sum               += st*u_tau*B*(y - ystar);
      sum               += st*nu*( ypls     *(log(ypls)      - 1.0) - 
                                   ystar_pls*(log(ystar_pls) - 1.0))/kappa;
      return sum;

    }

  }


  double integralJ(const double y, const double u_tau, const double tau, const double nu) { 

    const double a    = tau/nu;
    const double ypls = u_tau*y/nu;
    const double st   = SIGN_(tau);

    if ( ypls < ystar_pls) { 
    //if ( true) { 

      return a*y*y*y/3.0;

    } else { 

      const double ystar = ystar_pls*nu/u_tau;
      double sum         = a*ystar*ystar*ystar/3.0 ;
      sum               += 0.5*st*u_tau*B*( y*y - ystar*ystar);
      sum               += st*nu*nu*( 0.5*ypls*ypls*log(ypls) - 0.25*ypls*ypls - 
                                      0.5*ystar_pls*ystar_pls*log(ystar_pls) + 0.25*ystar_pls*ystar_pls)/(kappa*u_tau);
      return sum;

    }
    
  }

  double damp(const double ypls) { 

    const double Apls = 17.0;
    const double fax  = 1.0 - exp(-ypls/Apls);
    return fax*fax;

  }

  double func(const double u1, const double v1, const double dpdx, const double y1, const double nu, 
              const double tau, const double tau_prv, const double dt) { 

    const double alpha     = v1/(y1*y1);
    const double u_tau     = sqrt(abs(tau));
    const double u_tau_prv = sqrt(abs(tau_prv));


    // for now, we are going to integrate with a fixed dy 

    //const int n        = 500;
    const int n        = 40;
    const double dy    = y1/double(n);
    double y           = 0.0;

    double uu          = 0.0;
    for (int i =0; i < n; ++i) { 

      y                  += dy;

      const double ypls   = u_tau*y/nu;
      const double nu_t   = kappa*y*u_tau*damp(ypls);
      const double nu_tot = nu + nu_t;
      //const double unst   = (integralI(y,u_tau,tau,nu) - integralI(y,u_tau_prv,tau_prv,nu))/dt;
      const double unst   = 0.0;

      /*
      const double vv     = alpha*y*y*funcU(y,u_tau,tau,nu) - 4.0*alpha*integralJ(y,u_tau,tau,nu);
      const double pp     = dpdx*y;
      */

      const double vv     = alpha*y*y*funcU(y,u_tau,tau,nu) - 2.0*alpha*integralJ(y,u_tau,tau,nu);
      const double pp     = 0.0; 

      const double fax    = dy/nu_tot;
      uu                 += fax*(unst + vv + pp + tau);
      //uu += fax*tau;

    }

    assert( (y - y1)/y1 < 1.0e-12);

    return uu - u1;
  }

  double solve_tau(const double u1, const double v1, const double dpdx, const double nu, const double y1, 
                   const double tau_prv, const double dt) { 


    if ( sqrt(u1*u1 + v1*v1) < 1.0e-14) 
      return 0.0;

    // first attempt to solve via a newton iteration .. 
    
    double tau         = tau_prv; // use the previous iteration as the initial guess here.
    const double relax = 0.4;
    const double tol   = 1.0e-3;

    int iter ;
    const int max_iter = 40;
    for (iter = 0; iter < max_iter; ++iter) { 

      const double eps = 1.0e-3*max(abs(tau), u1*nu/y1);
      double f   = func(u1, v1, dpdx, y1, nu, tau, tau_prv, dt);

      if ( abs(f)/u1 < tol) { 
        break;
      }

      double fp           = func(u1, v1, dpdx, y1, nu, tau+eps,tau_prv,dt);
      const double dfdtau = (fp - f)/eps;

      if ( abs(dfdtau) < 1.0e-14) { 

        // stalled...
        iter = max_iter;
        break;

      } else { 

        const double dtau  = f/dfdtau;
        tau               += relax*dtau;

      }
    }

    if ( iter == max_iter) { 

      // the newton iteration failed to converge.  fall back to bisection.. 

      return solve_tau_bisection(u1,v1,dpdx,nu,y1,tau_prv,dt);

    } else { 

      return tau;

    }




  } 

  double solve_tau_bisection(const double u1, const double v1, const double dpdx, const double nu, const double y1,
                             const double tau_prv, const double dt) { 

    const double tau_nom = max(nu*u1/y1, abs(tau_prv));
    //const double alpha   = v1/y1/y1;
    
    double tau_hi = 5.0*tau_nom;
    double tau_lo = -1.0*tau_nom;
 
    double fhi    = func(u1, v1, dpdx, y1, nu, tau_hi, tau_prv, dt);
    double flo    = func(u1, v1, dpdx, y1, nu, tau_lo, tau_prv, dt);
    
    if ( (fhi > 0.0) && (flo > 0.0)) { 

      int iter = 0;
      while ( flo > 0.0) { 
        ++iter;
        tau_lo *= 2.0;
        flo = func(u1, v1, dpdx, y1, nu, tau_lo, tau_prv, dt);

        if ( iter > 10) { 
          print_failure(u1,v1,dpdx,y1,nu,tau_prv,dt);
          assert(0);
        }
      }

    } else if ( (fhi < 0.0) && (flo < 0.0)) { 

      int iter = 0;
      while ( fhi < 0.0) { 
        ++iter;
        tau_hi *= 2.0;
        fhi     = func(u1, v1, dpdx, y1, nu, tau_hi, tau_prv, dt);

        if ( iter > 10) { 
          print_failure(u1,v1,dpdx,y1,nu,tau_prv,dt);
          assert(0);
        }
      }

    }

#ifdef DEBUG
    cout << " tau_lo, tau_hi : " << tau_lo << "   " << tau_hi << endl;
    cout << " flo, fhi : " << flo << "   " << fhi << endl;
#endif

    assert( flo*fhi < 0.0);

    //const int max_iter = 500;
    const int max_iter = 100;
    //const double tol   = 1.0e-6;
    const double tol   = 1.0e-3;
    double tau_mid;

    int iter;
    for (iter = 0; iter < max_iter; ++iter) { 

      tau_mid              = 0.5*(tau_hi + tau_lo);
      const double fmid    = func(u1, v1, dpdx, y1, nu, tau_mid, tau_prv, dt);

#ifdef DEBUG
      printf("iter,res, tau:  %d  %12.8g  %12.8g\n", iter,fmid/u1,tau_mid);
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

    if ( iter >= max_iter) { 
      print_failure(u1,v1,dpdx,y1,nu,tau_prv,dt);
      assert( iter < max_iter);
    }

    return tau_mid;
  }

  void print_failure(double u1, double v1, double dpdx, double y1, double nu, double tau_prv, double dt) const { 
    cout << " failed: " << u1 << "    " << v1 << "    " << dpdx << "    " << y1 << "    " << nu << "    " << tau_prv << "    " << dt << endl;
  }

};



#endif
