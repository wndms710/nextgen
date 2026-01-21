#ifndef COMPRESSIBLESIMILARITY_HPP
#define COMPRESSIBLESIMILARITY_HPP

// solve the compressibile similarity in the limit of a 
// zero pressure gradient flat plate, with a linear viscosity
// temperature relationship (C=1 in white 7-32a)

#include <iostream>
#include <cmath>
#include <assert.h>
#include <vector>

extern "C" { 
  
  void dgtsv_(int*,int*,double*,double*,double*,double*,int*,int*);


  void dgtsvx_(char*, char*, 
               int*, int*, double*, double*, double*,
               double*, double*, double*, double*, int*,
               double*, int*, double*, int*, double*, double*, 
               double*, double*, int*, int*);


};

template<int n>
void solve_ge_inplace(double x[n], double b[n], double A[n*n]) { 
  
#define A(i,j) A[(j)*n+i]
  
  // gaussian elim .. 
  
  for (int k = 0; k < n-1; ++k) { 
    for (int i = k+1; i < n; ++i) { 
      
      A(i,k) = A(i,k)/A(k,k);
      for (int j = k+1; j < n; ++j) 
        A(i,j) -= A(i,k)*A(k,j);
      
    }
  }
  
  // forward elim .. 
  
  for (int k = 0; k < n-1; ++k) { 
    for (int i = k+1; i < n; ++i) { 
      
      b[i] -= A(i,k)*b[k];
      
    }
  }
  
  
  // backward solve .. 
  
  for (int i = n-1; i >= 0; --i) { 
    
    double s = b[i];
    for (int j = i+1; j < n; ++j) { 
      s -= A(i,j)*x[j];
    }
    
    x[i] = s/A(i,i);
    
  }
  
#undef A
}


using namespace std;

class ComSim { 

public: 

  double Ma_infty;
  double gamma;
  double Pr;
  double mu_power;

  vector<double> eta;
  vector<double> f;
  vector<double> fp;
  vector<double> fpp;
  vector<double> g;
  vector<double> gp;

  vector<double> g_int;
  vector<double> psi_int;

  double omega; 

  ComSim(double _Ma_infty, double _gamma, double _Pr, double _mu_power) { 

    Ma_infty = _Ma_infty;
    gamma    = _gamma;
    Pr       = _Pr;
    mu_power = _mu_power;

    omega    = -1.0;

  }

  // f(0) = f'(0) = 0
  // f'(\infty)   = 1

  void solve_f(const bool verbose = true) { 

    double fpp0 = 0.3; // initial guess.. 
    double res  = 1.0;
    const double tol = 1.0e-10;
    const double relax = 0.6;

    const double eps = 1.0e-8;
    int iter = 0;

    while ( abs(res) > tol) { 

      const double res2 = solve_f_(fpp0+eps); // XXX revisit the eps
      res = solve_f_(fpp0);

      const double drdfpp = (res2-res)/eps;
      fpp0 = fpp0 - relax*res/drdfpp;

      if ( verbose) 
        cout << "iter, res, fpp: " << iter << "  " << abs(res) << "   " << fpp0 << endl;
      

      ++iter;
      if ( iter > 100) 
        assert(0);
    }

  }

  void solve_fg_isothermal(const double g_wall, const bool verbose = true) { 

    // we need to shoot to find f'', g' at the wall.

    double fpp0 = 0.3;
    double gp0  = 0.2;

    const double tol = 1.0e-12;
    const double relax_f = 0.2;
    const double relax_g = 0.2;

    double res_f = 1.0;
    double res_g = 1.0;

    int iter = 0;

    while ( max(abs(res_f),abs(res_g)) > tol) { 

      // attempt newtwon .. build the two by two jacobian.

      const double eps_gp = 1.0e-8;
      const double eps_fpp = 1.0e-8;

      double J[4];

      double res_f_feps, res_g_feps;
      solve_fg_(res_f_feps,res_g_feps,fpp0+eps_fpp,g_wall,gp0);

      double res_f_geps, res_g_geps;
      solve_fg_(res_f_geps,res_g_geps,fpp0,g_wall,gp0+eps_gp);

      solve_fg_(res_f,res_g,fpp0,g_wall,gp0);

      J[(2)*0+0] = (res_f_feps - res_f)/eps_fpp;
      J[(2)*0+1] = (res_g_feps - res_g)/eps_fpp;

      J[(2)*1+0] = (res_f_geps - res_f)/eps_gp;
      J[(2)*1+1] = (res_g_geps - res_g)/eps_gp;


      const double detJ = J[3]*J[0] - J[1]*J[2];

      double dfpp0      = -(J[3]*res_f -J[2]*res_g)/detJ;
      double dgp0       = -(-J[1]*res_f + J[0]*res_g)/detJ;

      fpp0             += relax_f*dfpp0;
      gp0              += relax_g*dgp0;

      if ( verbose) 
        cout << "iter, res_f, res_g, fpp, gp: " << iter << "  " << abs(res_f) << "   " 
             << abs(res_g) << "   " << fpp0 << "    " << gp0 << endl;
      
      ++iter;
      if ( iter > 1000) 
        assert(0);

    }

  } 

  void solve_fg_(double &res_f, double& res_g, const double fpp0, const double g0, const double gp0) { 


    // (Cf'')' + f f'' = 0
    //  (Cg')' + Prfg' = -pr*C*(gamma-1)*ma^2*f''^2
    // unwrapped is 
    // 
    // Cf''' + (C' + f)f'' = 0
    // Cg''  + (C' + Pr*f)*g' = -pr*C*(gamma-1)*Ma^2f''^2
    // 
    // dividing by C yields
    // 
    //  f''' + (C' + f)f''/C = 0
    //  g''  + (C' + Pr*f)/C*g' = -pr*(gamma-1)*Ma^2f''^2 

    // 


    eta.clear();
    double eta_ = 0.0;
    eta.push_back(eta_);
    double y[5] = {0.0, 0.0, fpp0, g0, gp0};


    f.clear();
    fp.clear();
    fpp.clear();
    g.clear();
    gp.clear();

    f.push_back(y[0]);
    fp.push_back(y[1]);
    fpp.push_back(y[2]);
    g.push_back(y[3]);
    gp.push_back(y[4]);

    const double eta_infty = 100.0;
    double deta            = 1.0e-7; 

    const double B         = -Pr*(gamma-1.0)*Ma_infty*Ma_infty;

    while( eta_ < eta_infty) { 
      
      eta_ += deta;

#define J(i,j) J[5*(j)+i]
      double J[25];
      
      for (int ii = 0; ii < 25; ++ii) 
        J[ii] = 0.0;

      const double C      = pow(y[3],mu_power-1.0);
      const double Cprime = (mu_power-1.0)*pow(y[3],mu_power-2.0)*y[4];

      J(0,1) = 1.0;
      J(1,2) = 1.0;
      J(2,2) = -(Cprime + y[0])/C;

      J(3,4) = 1.0;
      J(4,2) = B*y[2];
      J(4,4) = -(Cprime + Pr*y[0])/C;


      double r[5];
      for (int i = 0; i < 5; ++i) { 

        r[i] = y[i];
        for (int j = 0; j < 5; ++j) 
          r[i] += 0.5*deta*J(i,j)*y[j];

      }


      // modify the jacobian for the lhs crank op

      for (int i = 0; i < 5; ++i) { 

        for (int j = 0; j < 5; ++j) { 

          J(i,j) *= -0.5*deta;
          if ( i == j) 
            J(i,j) += 1.0;

        }
      }
#undef J

      solve_ge_inplace<5>(y,r,J);
     
      eta.push_back(eta_);
      f.push_back(y[0]);
      fp.push_back(y[1]);
      fpp.push_back(y[2]);
      g.push_back(y[3]);
      gp.push_back(y[4]);

      //deta = min(eta_infty-eta_,min(0.9,deta*1.005));
      deta = min(eta_infty-eta_,min(0.001*abs(y[0]/(y[1] + 1.0e-10)),0.9));

    } 

    res_f = y[1] - 1.0;
    res_g = y[3] - 1.0;
    
  }

  void compute_integrals(const bool verbose=true) { 

    assert( eta.size() == g.size());
    const int n = eta.size();

    g_int.clear();
    double sum = 0.0;
    omega      = 0.0;
    g_int.push_back(sum);

    for (int i = 1; i < n; ++i) { 

      const double deta     = eta[i] - eta[i-1]; 
      sum += deta*g[i];
      g_int.push_back(sum);

      omega += (g[i] - fp[i])*deta;

    }
    
    if ( verbose) 
      cout << " omega : " << omega << endl;

  } 

  void getInterval(int& i_lo, int& i_hi, const vector<double>& arr, const double val) const { 

    const int n    = eta.size();

    if ( val > arr[n-1]) { 

      i_lo = n-1;
      i_hi = n-1;

    } else if ( val < arr[0]) { 

      i_lo = 0;
      i_hi = 0;

    } else { 

      // we are looking for eta s.t. arr(eta) = val
      
      i_lo       = 0;
      i_hi       = n-1;
      
      assert( arr[i_hi] >= val);
      assert( arr[i_lo] <= val);
      
      while ( i_hi-i_lo > 1) { 
        
        const int i_mi    = (i_hi+i_lo)/2;
        
        if ( arr[i_mi] < val)  
          i_lo = i_mi;
        else  
          i_hi = i_mi;
        
      }
      
    }
  }

  void compute_psi_int() { 

    const int n = eta.size();
    psi_int.clear();

    double sum = 0.0;
    psi_int.push_back(sum);

    for (int i = 1; i < n; ++i) { 

      const double H  = fp[i]*g_int[i] - f[i]*g[i];
      const double Hp = fpp[i]*g_int[i] - f[i]*gp[i];

      const double psi = (H*Hp - fp[i]*(2.0*g[i]*H + Hp*g_int[i]))/g[i]/g[i];
      sum             += psi*(eta[i+1]-eta[i]);
      psi_int.push_back(sum);
    }

  } 


  void print_psi() const { 

    const int n = eta.size();

    for (int i = 0; i < n; ++i) { 

      const double H  = fp[i]*g_int[i] - f[i]*g[i];
      const double Hp = fpp[i]*g_int[i] - f[i]*gp[i];

      const double psi = (H*Hp - fp[i]*(2.0*g[i]*H + Hp*g_int[i]))/g[i]/g[i];
      cout << " ZZZZ eta psi H = " << eta[i] << "    " << psi << "   " << H << endl;
    }


  } 

  double getEtaFromYonDeltaStar(const double y_hat) const { 

    const double r = omega*y_hat;
    const int n    = eta.size();
    int i_lo, i_hi;

    getInterval(i_lo,i_hi,g_int,r);

    if ( (i_lo == 0)&& (i_hi == 0)) { 

      return eta[0];

    } else if ( (i_lo == n-1)&&(i_hi == n-1)) { 

      return eta[n-1];

    } else { 

      // interpolate on the interval.. 
    
      return eta[i_lo] + (r-g_int[i_lo])*(eta[i_hi]-eta[i_lo])/(g_int[i_hi]-g_int[i_lo]);

    }
  }

  double getUFromEta(const double eta_) const { 

    const int n = eta.size();

    int i_lo,i_hi;
    getInterval(i_lo,i_hi,eta,eta_);

    if ( (i_lo ==0) && (i_hi == 0)) { 

      return fp[0];

    } else if ( (i_lo == n-1) && (i_hi == n-1)) { 

      return fp[n-1];

    } else { 

      return fp[i_lo] + (fp[i_hi]-fp[i_lo])/(eta[i_hi]-eta[i_lo])*(eta_-eta[i_lo]); 

    } 

  } 

  double getPsiIntFromEta(const double eta_) { 

    const int n = eta.size();
    
    int i_lo,i_hi;
    getInterval(i_lo,i_hi,eta,eta_);
    
    if ( (i_lo ==0) && (i_hi == 0)) { 
      
      return psi_int[0];

    } else if ( (i_lo == n-1) && (i_hi == n-1)) { 

      return psi_int[n-1];

    } else { 

      return psi_int[i_lo] + (psi_int[i_hi]-psi_int[i_lo])/(eta[i_hi]-eta[i_lo])*(eta_-eta[i_lo]); 

    } 

  } 

  double getVFromEta(const double eta_, const double Re_x) { 

    const int n = eta.size();
    
    int i_lo,i_hi;
    getInterval(i_lo,i_hi,eta,eta_);

    const double fax = pow(2.0*Re_x,-0.5);

    if ( (i_lo ==0) && (i_hi == 0)) { 

      return fax*( fp[0]*g_int[0] - f[0]*g[0]);

    } else if ( (i_lo == n-1) && (i_hi == n-1)) { 

      return fax*( fp[n-1]*g_int[n-1] - f[n-1]*g[n-1]); 

    } else { 

      const double v_lo = fax*(fp[i_lo]*g_int[i_lo] - f[i_lo]*g[i_lo]);
      const double v_hi = fax*(fp[i_hi]*g_int[i_hi] - f[i_hi]*g[i_hi]);

      return v_lo + (v_hi - v_lo)/(eta[i_hi]-eta[i_lo])*(eta_-eta[i_lo]); 

    } 

  } 

  double getRhoFromEta(const double eta_) const { 

    const int n = eta.size();
    
    int i_lo,i_hi;
    getInterval(i_lo,i_hi,eta,eta_);

    if ( (i_lo ==0) && (i_hi == 0)) { 

      return 1.0/g[0];

    } else if ( (i_lo == n-1) && (i_hi == n-1)) { 

      return 1.0/g[n-1];

    } else { 

      return 1.0/g[i_lo] + (1.0/g[i_hi] - 1.0/g[i_lo])/(eta[i_hi]-eta[i_lo])*(eta_-eta[i_lo]); 

    } 

  } 

  double solve_f_(const double fpp0) { 

    // 
    // f''' + f f'' = 0
    // unwrapped is 
    // 
    // (f2)'    = -f0*f2
    // (f1)'    =  f2
    // (f0)'    =  f1
    // 
    // subject to initial conditions of f0 = 0, f1 = 0, f2 = fpp0

    eta.clear();
    double eta_ = 0.0;
    eta.push_back(eta_);
    double y[3] = {0.0, 0.0, fpp0};


    f.clear();
    fp.clear();
    fpp.clear();

    f.push_back(y[0]);
    fp.push_back(y[1]);
    fpp.push_back(y[2]);

    const double eta_infty = 100.0;
    double deta            = 1.0e-5; 

    while( eta_ < eta_infty) { 
      
      eta_ += deta;

#define J(i,j) J[3*(j)+i]
      double J[9];
      J(0,0) = 0.0;
      J(0,1) = 1.0;
      J(0,2) = 0.0;

      J(1,0) = 0.0;
      J(1,1) = 0.0;
      J(1,2) = 1.0;

      J(2,0) = 0.0; 
      J(2,1) = 0.0;
      J(2,2) = -y[0];


      double r[3];
      for (int i = 0; i < 3; ++i) { 

        r[i] = y[i];
        for (int j = 0; j < 3; ++j) 
          r[i] += 0.5*deta*J(i,j)*y[j];

      }


      // modify the jacobian for the lhs crank op

      for (int i = 0; i < 3; ++i) { 

        for (int j = 0; j < 3; ++j) { 

          J(i,j) *= -0.5*deta;
          if ( i == j) 
            J(i,j) += 1.0;

        }
      }

      // the resultant J is upper triangular.. back substitute.

      assert( J(2,0) == 0.0);
      assert( J(2,1) == 0.0);
      y[2] = r[2]/J(2,2);

      assert( J(1,0) == 0.0);
      y[1] = (r[1] - J(1,2)*y[2])/J(1,1);

      y[0] = (r[0] - J(0,2)*y[2] - J(0,1)*y[1])/J(0,0);
 
#undef J
      eta.push_back(eta_);
      f.push_back(y[0]);
      fp.push_back(y[1]);
      fpp.push_back(y[2]);

      //deta = min(eta_infty-eta_,min(0.9,deta*1.005));
      deta = min(eta_infty-eta_,min(0.0025*abs(y[0]/(y[1] + 1.0e-10)),0.9));

    } 

    return y[1]-1.0;
    
  }

#ifdef NOPE

  // the following is not advised because it quickly becomes terribly conditioned

  void solve_g_isothermal(const double Tw_hat) { 
    
    // assumes that we have solved for f at this stage -- in the limit of C = 1, 
    // the two equations are uncoupled
    // 
    // Tw_hat = T_w/T_infty

    // solve this on the same grid as the f values .. 

    int n = eta.size();

    // finite difference approximation of the two pt boundary value problem.. 

    const double fax = -Pr*(gamma-1.0)*Ma_infty*Ma_infty;

    double * r       = new double[n];
    double * subd    = new double[n-1];
    double * diag    = new double[n];
    double * supd    = new double[n-1];

    for (int ii = 0; ii < n-1; ++ii) { 

      subd[ii] = 0.0;
      supd[ii] = 0.0;

    } 

    for (int ii = 0; ii < n; ++ii) 
      diag[ii] = 0.0;


    for (int i = 0; i < n; ++i) { 

      if ( i == 0 ) { 

        const double inv_denom = 1.0/pow(eta[i+1]-eta[i],2);
        diag[i] = -1.0*inv_denom;
        r[i]    = -Tw_hat*inv_denom;

      } else if ( i == n-1) { 

        const double inv_denom = 1.0/pow(eta[i]-eta[i-1],2);
        diag[i] = -1.0*inv_denom;
        r[i]    = -1.0*inv_denom;

      } else { 

        r[i] = fax*fpp[i]*fpp[i];
        
        // second derivative .. 

        const double eta_plus = 0.5*(eta[i] + eta[i+1]);
        const double eta_mins = 0.5*(eta[i] + eta[i-1]);
        const double deta     = eta_plus-eta_mins;

        const double aa = 1.0/(eta[i] - eta[i-1])/deta;
        const double bb = 1.0/(eta[i+1]-eta[i]  )/deta;

        subd[i-1]       += aa;
        supd[i]         += bb;
        diag[i]         -= aa+bb;

        // first derivative term.. 
       
        const double tmp = Pr*f[i];
        subd[i-1]       -= tmp/(eta[i+1]-eta[i-1]);
        supd[i]         += tmp/(eta[i+1]-eta[i-1]);

      }

    }

    
    double * r_cpy    = new double[n];
    double * diag_cpy = new double[n];
    double * supd_cpy = new double[n-1];
    double * subd_cpy = new double[n-1];

    for (int i = 0; i < n-1; ++i) { 
      
      subd_cpy[i] = subd[i];
      supd_cpy[i] = supd[i];

    } 

    for (int i = 0; i < n; ++i) { 

      diag_cpy[i] = diag[i];
      r_cpy[i]    = r[i];

    } 


    //thomas_solve(r,subd,diag,supd,n);

    { 

      int ione = 1;
      int info;
      dgtsv_(&n,&ione,subd,diag,supd,r,&n,&info);

      if ( info != 0) 
        cout << " uh-oh : " << info << endl;

    }

    { 

      char fact = 'n';
      char trans = 'n';
      int  ione  = 1;

      double * dlf = new double[n];
      double * df  = new double[n];
      double * duf = new double[n];
      double * du2 = new double[n];

      int* ipiv    = new int[n];
      double* x    = new double[n];

      double rcond; 
      double ferr, berr;
      double* work = new double[5*n];
      int* iwork   = new int[n];
      int info;

      dgtsvx_(&fact, &trans, &n, &ione,
              subd,diag,supd,dlf,df,duf,du2, 
              ipiv,r,&n,x,&n,&rcond,&ferr,&berr,work,
              iwork, &info);

      if ( (info != 0) && (info != n+1)) 
        cout << " info issue : " << info << "   " << n << endl;
      else if ( info == n+1) 
        cout << " rcond suggests matrix is singular to machine prec: " << rcond << endl;
      else 
        cout << " rcond = " << rcond << endl;

      for (int i = 0; i < n; ++i) 
        r[i] = x[i];

      delete[] dlf;
      delete[] df;
      delete[] duf;
      delete[] du2;
      delete[] ipiv;
      delete[] x;
      delete[] work;
      delete[] iwork;

    }




    // check the solution .. 

    double linf = 0.0;

    for (int i = 0; i < n; ++i) { 

      double res;

      if ( i ==0 ) { 

        res = diag_cpy[i]*r[i] + supd_cpy[i]*r[i+1] - r_cpy[i];

      } else if ( i == n-1) { 

        res = subd_cpy[i-1]*r[i-1] + diag_cpy[i]*r[i] - r_cpy[i];

      } else { 

        res = subd_cpy[i-1]*r[i-1] + diag_cpy[i]*r[i] + supd_cpy[i]*r[i+1] - r_cpy[i];

      } 

      cout << i << "  " << res << endl;
      linf = max(fabs(res),linf);

    } 

    cout << " > residual in tridiag solve = " << linf << endl;

    delete[] r_cpy;
    delete[] supd_cpy;
    delete[] diag_cpy;
    delete[] subd_cpy;

    g.clear();

    for (int i = 0; i < n; ++i) 
      g.push_back(r[i]);


    delete[] r;
    delete[] supd;
    delete[] diag;
    delete[] subd;

  } 
#endif

  void thomas_solve(double* r, double* subd, double* diag, double* supd, const int n) { 

    r[0] /= diag[0];
    supd[0] /= diag[0];

    for (int i = 1; i < n; ++i) { 

      const double inv_denom = 1.0/( diag[i] - subd[i-1]*supd[i-1]);
      if ( i < n-1) 
        supd[i]           *= inv_denom;
      r[i]               = (r[i] - subd[i-1]*r[i-1])*inv_denom;

    }

    for (int i = n-2; i >= 0; --i) { 

      r[i] -= supd[i]*r[i+1];

    }
    
  }

  void dump_fp() const { 

    const int n = fp.size();
    assert( eta.size() == n);

    for (int i = 0; i < n; ++i) { 

      cout << " YYYY " << eta[i] << "   " << fp[i] << endl;

    }


  } 

  void dump_fpg() const { 

    const int n = fp.size();
    assert( eta.size() == n);
    assert( g.size() == n);

    for (int i = 0; i < n; ++i) { 

      cout << " YYYY " << eta[i] << "   " << fp[i] << "   " << g[i] << endl;

    }


  } 


};



#endif
