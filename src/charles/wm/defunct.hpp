#ifndef NONPARALLELWM_HPP
#define NONPARALLELWM_HPP

#include <iostream>
#include <cstdio>
#include <cmath>
#include <assert.h>

using namespace std;

class NonParallelWm { 
public: 

  double * y;
  double * vol;
  double * u;
  int n;

  NonParallelWm() { 
    
    n   = 0;
    y   = NULL; 
    vol = NULL; 
    u   = NULL; 

  }

  ~NonParallelWm() { 

    if ( y != NULL) delete[] y;
    if ( y != NULL) delete[] vol;
    if ( y != NULL) delete[] u;

  }


  /*
  void build_grid(const double y1, const double dplus) { 

    // y denotes a cell center of the grid.

    for (int iter = 0; iter < 2; ++iter ) { 

      double yy    = 0.0;
      double dy    = dplus;
      const int gf = 1.2;
      n            = 0;

      while ( yy < y1) { 
        
        yy += dy;
        
        dy = min(dy*gf, 0.1*y1); 
        
        if ( iter == 1) 
          y[n] = yy;

        n++;
        
      }

      if ( iter == 0) { 

        if ( y != NULL) delete[] y;
        y = new double[n];

        if ( vol != NULL) delete[] vol;
        vol = new double[n];

        if ( u != NULL) delete[] u;
        u = new double[n];

      }
    }

    this->build_volume();

  } 
  */

  void build_grid(const double y1, const double dplus) { 

    if ( y == NULL) y = new double[n];
    if ( vol == NULL) vol = new double[n];
   
    // y = zeta^p
    // d+ = y1* (dz/2)^p => log(d+/y1)/log(dz/2) = p

    const double dz = 1.0/double(n);
    //const double p  = log(dplus/y1)/log(0.5*dz);
    const double p = 1.0;
    assert(p > 0.0);

    cout << " power ; " << p << endl;

    for (int i = 0; i < n; ++i)  
      y[i] = y1* pow( (double(i)+0.5)*dz, p);

    this->build_volume(y1);

  }

  void build_volume(const double y1) { 

    vol[0]   = 0.5*(y[0] + y[1]) - 0.0;

    for (int i =1 ; i < n-1; ++i) 
      vol[i] = 0.5*(y[i] + y[i+1]) - 0.5*(y[i] + y[i-1]);

    vol[n-1] = y1 - 0.5*(y[n-1] + y[n-2]);

  }

  void dump_solution() const { 

    FILE * fp = fopen("dump.dat", "w");
    for (int i = 0; i < n; ++i) 
      fprintf(fp, "%12.8g  %12.8g\n", y[i], u[i]);
    fclose(fp);

  }

  double solve(const double y1, const double v1, const double u1, const double nu, const int _n=100) { 

    // set the initial guess .. 

    n = _n;

    assert( u == NULL); u = new double[n];

    double tau_hat = nu*u1/y1;
    double dplus = nu/sqrt(abs(tau_hat));
    this->build_grid(y1,dplus);
    
    for (int i =0; i <n; ++i) 
      u[i] = y[i]/y1*u1;

    const double tol     = 1.0e-04;
    const int max_iter   = 5000;
    const double relax   = 0.3;
    int iter;

    for (iter = 0; iter < max_iter; ++iter) { 

      //dplus = nu/sqrt(abs(tau_hat));
      //this->build_grid(y1,dplus);

      double tau_old = tau_hat;
      tau_hat        = solve_iter(tau_old, y1, v1, u1, nu);

      // relax .. 
      tau_hat        = (1.0-relax)*tau_old + relax*tau_hat;

      /*
      const double u2 = u1*u1 + v1*v1;

      if ( abs(tau_old - tau_hat)/u2 < tol) 
        break;

      */

      if ( abs(tau_old-tau_hat)/tau_hat < tol) 
        break;

      cout << " iter, tau_hat : " << iter << "   " << tau_hat << endl;

    }
   
    dump_solution();

    return tau_hat;
  }

  inline double damp(const double y_plus) const { 

    const double A_plus = 26.0;
    const double fax = 1.0 - exp(-y_plus/A_plus);
    return fax*fax;

  }


  double solve_iter(const double tau_hat, const double y1, const double v1, const double u1, const double nu) { 

    const double alpha  = v1/(y1*y1);
    const double u_tau  = sqrt(abs(tau_hat));
    //const double u_tau = 0.123;
    const double kappa  = 0.41;

    double * diag       = new double[n];
    double * subd       = new double[n];
    double * supd       = new double[n];
    double * rhs        = new double[n];

    for (int i = 0; i < n; ++i) { 

      diag[i] = 0.0;
      subd[i] = 0.0;
      supd[i] = 0.0;
      rhs[i]  = 0.0;

    }

    // interior faces .. 

    for (int i = 0; i < n-1; ++i) { 

      const double yf       = 0.5*(y[i] + y[i+1]);
      const double yf_plus  = yf*u_tau/nu;
      const double dudy     = (u[i+1] - u[i])/(y[i+1] - y[i]);
      const double nu_t     = kappa*kappa*yf*yf*abs(dudy)*damp(yf_plus);
      //const double nu_t     = kappa*u_tau*yf*damp(yf_plus);
      const double nu_tot   = nu + nu_t;

      //
      // flux = nu_tot*(u[i+1] - u[i])/(y[i+1]-y[i]) - alpha*0.5*(u[i+1] + u[i])*yf*yf
      //

      const double fax  = alpha*0.5*yf*yf;
      const double fax2 = nu_tot/(y[i+1] - y[i]);

      diag[i]   += -fax2 - fax; 
      supd[i]   +=  fax2 - fax; 

      diag[i+1] += -fax2 + fax;
      subd[i  ] +=  fax2 + fax;

    }


    // wall boundary , flux = nu* (u[0] - 0.0)/(y[0] - 0.0) = nu* u[0]/y[0];
    
    diag[0] -= nu/y[0];


    // far field boundary .. 
    // flux = nu_tot*(u1 - u[n-1])/(y1 - y[n-1]) - alpha*u1*y1*y1

    { 
      const double y1_plus = y1*u_tau/nu;
      const double dudy    = (u1 - u[n-1])/(y1 - y[n-1]);
      const double nu_t    = kappa*kappa*y1*y1*abs(dudy)*damp(y1_plus);
      const double nu_tot  = nu + nu_t;
      
      diag[n-1]           -= nu_tot/(y1 - y[n-1]);
      rhs[n-1]            -= nu_tot/(y1 - y[n-1])*u1;
      rhs[n-1]            += alpha*u1*y1*y1;

    }


    // and lastly the diagonl terms.. 

    for (int i = 0; i <n ; ++i) 
      diag[i] += 2.0*vol[i]*alpha*y[i];


    // as a placeholder, we'll use a thomas algorithm to solve 
    // the tridigonal system.  note that due to the presence of 
    // the v velocity, we may require implementing pivoting.. 

    { 
      double * tmp = new double[n];
      for (int i = 0; i < n; ++i) { 

        u[i]   = 0.0;
        tmp[i] = 0.0;

      }

      double w = diag[0];
      u[0]     = rhs[0]/w;
      
      for (int i = 1; i < n; ++i) { 

        tmp[i-1] = supd[i-1]/w;
        w        = diag[i] - subd[i-1]*tmp[i-1];
        u[i]     = (rhs[i] - subd[i-1]*u[i-1])/w;

      }

      // back substitution .. 

      for (int i = n-2; i >= 0; --i) { 
        u[i] = u[i] - tmp[i]*u[i+1];
      }


      // lets check the result -- to make sure we actually solved the equation .. 

      double max_diff;
      rhs[0]         -= diag[0]*u[0] + supd[0]*u[1];
      max_diff        = abs(rhs[0]);

      for (int i = 1; i < n-1; ++i) {
        rhs[i] -= subd[i-1]*u[i-1] + diag[i]*u[i] + supd[i]*u[i+1];
        max_diff = max(max_diff,abs(rhs[i]));
      }
        
      rhs[n-1] -= subd[n-2]*u[n-2] + diag[n-1]*u[n-1];
      max_diff  = max(max_diff,abs(rhs[n-1]));


      cout << " > max error ( should be zero) : " << max_diff << endl;

      delete[] tmp;

    }

    delete[] diag;
    delete[] supd;
    delete[] subd;
    delete[] rhs;

    const double tau_hat_new = nu*u[0]/y[0];
    return tau_hat_new;

  }
};


#endif
