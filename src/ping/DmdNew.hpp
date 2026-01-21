#ifndef DMDNEW_HPP
#define DMDNEW_HPP

#include "Lapack.hpp"
#include <iostream>
#include <assert.h>
#include <cstdio>

using namespace std;

#define corr(i,j) corr[(j)*ldc + i]
#define vr(i,j) vr[(j)*nnz+i]
#define Wk(i,j) Wk[(j)*ns+i]
#define A_tilde(i,j) A_tilde[(j)*nnz+i]
#define tmp2(i,j) tmp2[(j)*ns+i]
#define tmp(i,j) tmp[(j)*ns+i]
#define vr_imag(i,j) vr_imag[(j)*nnz+i]

int compute_dmd(const int n, const double * corr, double* &lambda_r, double* &lambda_i,
                double* &Wk, complex<double> * &vr_imag, double *&amp, const double delta_T) { 

  // n : total number of snapshots .. 
  // corr : corerlation matrix 

  assert ( lambda_r == NULL);
  assert ( lambda_i == NULL);
  assert ( Wk       == NULL);
  assert ( vr_imag  == NULL);

  // first lets extract the prinicipal (n-1) submatrix of the correlation .. 

  int ldc = n;
  int ns = n-1;
  
  double * sigma = new double[ns];
  Wk             = new double[ns*ns];
  double * work  = NULL;

  for (int j =0; j < ns; ++j) { 
    for (int i = 0; i < ns; ++i) 
      Wk(i,j) = corr(i,j);
  }

  // Wk on entry  = X^T X = V \Sigma^2 V^T
  // where X = U \Sigma V^T .. 
  
  char jobz = 'v';
  char uplo = 'u';
  int info;
  int lwork = 5*ns;
  work      = new double[lwork];

  cout << " Computing eigs of correlation matrix ... " ; cout.flush();
  dsyev_(&jobz, &uplo, &ns, Wk, &ns, sigma, work, &lwork, &info);

  if ( info != 0) 
    cout << " failed to compute eigs of corr submatrix : " << info << endl;
  else 
    cout << " ok. " << endl;

  const double eps_sigma = 1.0e-14;
  int nnz = 0 ;
  for (int i =0; i < ns ; ++i) {
    // assuming that the data is scaled appropriate on entry .. 
    if ( (fabs(sigma[i]) > eps_sigma)) { 
      assert(sigma[i] >= 0.0);
      ++nnz;
    }
  }//i
  cout << " Rank = " << nnz << endl ; cout.flush();
  
  if ( nnz < ns ) {
    
    // arrange the eigenvalues and eigenvectors in descending order.
    int shft = ns -nnz ;
    for (int j=ns-nnz; j < ns ; ++j) {
      sigma[j-shft] = sigma[j] ;
      for (int i = 0 ; i < ns ; ++i)
      Wk(i,j-shft) = Wk(i,j);
    }//j
  }

  for (int i = 0; i < nnz; ++i) { 
    assert ( sigma[i] >= eps_sigma);
    sigma[i] = sqrt(sigma[i]);
  }

  // at this stage we have the V (frm the svd of X) stored in Wk
  // and the non-zero singular values in Wk.. rescale the columns 
  // of V with sigma^{-1} --> V_tilde

  for (int j = 0; j < nnz; ++j) { 
    for (int i = 0 ; i < ns; ++i) 
      Wk(i,j) /= sigma[j];
  }

  double * A_tilde = new double[nnz * nnz];
  double * tmp     = new double[ns *ns];

  // A_tilde = Wk^T (X^TY) V.. 
  // first extract X^T Y from the correlation matrix .. 

  for (int j = 0; j < ns; ++j) { 
    for (int i = 0; i < ns; ++i) { 
      tmp(i,j) = corr(i,j+1);
    }
  }

  // then multiply (X^T Y) V_tilde = [ns X ns] X [ns X nnz] = [ns X nnz]
  
  double * tmp2 = new double[ns * nnz];

  for (int j = 0; j < nnz; ++j) { 
    for (int i =0 ; i < ns; ++i) { 
      tmp2(i,j) = 0.0;
      for (int k = 0; k < ns; ++k) 
        tmp2(i,j) += tmp(i,k) * Wk(k,j);
    }
  }

  for (int j = 0; j < nnz; ++j) { 
    for (int i = 0; i < nnz; ++i) { 
      
      A_tilde(i,j) = 0.0;
      for (int k = 0; k < ns; ++k) { 
        A_tilde(i,j) += Wk(k,i) * tmp2(k,j);
      }
    }
  }

  // now we are prepared to compute the dmd eigenvalues .. 

  double * eig_r = new double[nnz];
  double * eig_i = new double[nnz];
  double * vr    = new double[nnz*nnz];
  
  delete[] work;
  char jobvl = 'n';
  char jobvr = 'v';
  lwork = 20*nnz ;
  work  = new double[lwork];
 
  double vl;
  int ldvl = 1;
  int ldvr = nnz;

  cout << " Computing dmd eigs ... " ; cout.flush();
  dgeev_( &jobvl, &jobvr, &nnz, A_tilde, &nnz, eig_r, eig_i, &vl, &ldvl,
         vr, &ldvr, work, &lwork, &info ) ;
  
  if ( info != 0) 
    cout << " failed to compute the dmd eigs : " << info << endl;
  else 
    cout << " ok. " << endl;

  // store the eigenvectors of the A_tilde matrix ... 

  vr_imag = new complex<double>[nnz*nnz];
  for (int j= 0; j < nnz ; ++j) {
    if ( fabs(eig_i[j]) < 1.0e-12) {
      for (int i =0; i < nnz ; ++i) {
        vr_imag(i,j) = complex<double>(vr(i,j), 0.0); // real eigenvector
      }
    } else {
      for (int i =0; i < nnz ; ++i) {
        vr_imag(i,j)   = complex<double>(vr(i,j), vr(i,j+1));
        vr_imag(i,j+1) = complex<double>(vr(i,j),-vr(i,j+1));
      }
      ++j; // complex pair handled.
    }
  }//j


  lambda_r = new double[nnz];
  lambda_i = new double[nnz];

  for (int i =0; i < nnz; ++i) { 

    complex<double> tt(eig_r[i],eig_i[i]);
    const double f      = arg(tt);
    const double sigma  = log(abs(tt));

    lambda_i[i]         = f/ (2.0*M_PI*delta_T);
    lambda_r[i]         = sigma/delta_T;

    printf(" > Freq, gr : %12.8g  %12.8g  %12.8g  %12.8g\n", lambda_i[i], lambda_r[i], eig_r[i], eig_i[i]);
  }

  delete[] A_tilde;
  delete[] tmp;
  delete[] tmp2;
  delete[] work;
  delete[] vr;
  delete[] eig_r;
  delete[] eig_i;
  delete[] sigma;

  return nnz;
} 

#undef corr
#undef vr
#undef Wk
#undef Wk
#undef A_tilde
#undef tmp2
#undef tmp
#undef vr_imag

void compute_dmd_adjoint(const int ns, const double* A_corr,const double * Wk,double* &lambda_r_adj,
                         double* &lambda_i_adj, complex<double>* & vrh_imag,const double dt_samp,const int nnz) { 

  
  // recall V_tilde = [nt X nnz ] 
  int nt = ns - 1;
  
  // T = (Y^* X) V_tilde = [nt x nt] X [nt X nnz ] = [nt * nnz]
  
  double * T = new double[nt * nnz];
  for (int j = 0; j < nnz; ++j) { 
    for (int i = 0; i < nt; ++i) { 
      
      T[(j)*nt+i] = 0.0;
      for (int k =0; k < nt; ++k) {
        // phi_{i+1}^T phi_{k}
        T[(j)*nt+i] += A_corr[(i+1)*ns+k] * Wk[(j)*nt+k];
      }
    }
  }
  
  // R = V_tilde^H * T = [nnz * nt] X [nt X nnz] = [nnz * nnz]
  
  double * R = new double[nnz * nnz];
  for (int j = 0; j < nnz; ++j) { 
    for (int i =0; i < nnz; ++i) { 
      
      R[(j)*nnz+i] = 0.0;
      for (int k = 0; k < nt; ++k) { 
        R[(j)*nnz+i] += Wk[(i)*nt+k] * T[(j)*nt+k];
      }
    }
  }
  
  char jobvr = 'v';
  char jobvl = 'n';
  int lwork = 20*nnz;
  double * work = new double[lwork];
  int info;
  
  double vl; 
  int ldvl = 1;
  
  double * vr_tmp = new double[nnz*nnz];
  assert( lambda_r_adj == NULL); lambda_r_adj = new double[nnz];
  assert( lambda_i_adj == NULL); lambda_i_adj = new double[nnz]; 

  int nnz_tmp = nnz;
  dgeev_(&jobvl,&jobvr,&nnz_tmp,R,&nnz_tmp,lambda_r_adj,lambda_i_adj,
         &vl,&ldvl,vr_tmp, &nnz_tmp, work, &lwork, &info);
  
  if ( info != 0) { 
    cout << " proj adj dgeev error : " << info << endl;
    return;
  }
 

  // repack the vr into the imaginary matrix ... 

  assert( vrh_imag == NULL); vrh_imag = new complex<double>[nnz*nnz];

  for (int j= 0; j < nnz ; ++j) {
    if ( fabs(lambda_i_adj[j]) < 1.0e-12) {
      for (int i =0; i < nnz ; ++i) {
        vrh_imag[(j)*nnz+i] = complex<double>(vr_tmp[(j)*nnz+i],0.0); // real eigenvector
      }
    } else {
      for (int i =0; i < nnz ; ++i) {
        vrh_imag[(j)*nnz+i]   = complex<double>(vr_tmp[(j)*nnz+i], vr_tmp[(j+1)*nnz+i]);
        vrh_imag[(j+1)*nnz+i] = complex<double>(vr_tmp[(j)*nnz+i],-vr_tmp[(j+1)*nnz+i]);
      }
      ++j; // complex pair handled.
    }
  }//j

  delete[] vr_tmp;
  delete[] work;
  delete[] R;
  delete[] T;

}

void compute_dmd_amp(double* &amp, const int ns, const double* A_corr, 
                     const double* Wk, const double* lambda_r, 
                     const double* lambda_i, complex<double>* vr_imag, 
                     const double dt_samp, const int nnz) {

  // recall V_tilde = [nt X nnz]

  int nt = ns - 1;

  // compute What = V_tilde W = [nt X nnz] X [nnz X nnz] 

  cout << " Computing What .. " ;

  complex<double>* What = new complex<double>[nt * nnz];

  for (int j = 0; j < nnz; ++j) { 
    for (int i = 0; i < nt; ++i) { 

      What[(j)*nt+i] = complex<double>(0.0,0.0);

      for (int k = 0; k < nnz; ++k) { 
        What[(j)*nt+i] += complex<double>(Wk[(k)*nt+i])*vr_imag[(j)*nnz+k];
      }
    }
  }

  cout << " ok. " << endl;

  // need to form T = What^* A' What .. = What^* R
  // R = [nt X nt] X [nt X nnz] = [nt X nnz]
  
  cout << " Computing R ... " ;

  complex<double>* R = new complex<double>[nt * nnz];

  for (int j = 0; j < nnz; ++j) { 
    for (int i = 0; i < nt; ++i) { 

      R[(j)*nt+i] = complex<double>(0.0,0.0);

      for (int k = 0; k < nt; ++k) { 
        // recall that A' = Y^*Y
        R[(j)*nt+i] += complex<double>(A_corr[(k+1)*ns+(i+1)])*What[(j)*nt+k];
      }
    }
  }

  cout << " ok." << endl;

  // T = What^* R = [nnz X nt] X [nt X nnz] = [nnz X nnz]

  cout << " Computing T... ";

  complex<double>* T = new complex<double>[nnz * nnz];

  for (int j = 0; j < nnz; ++j) { 
    for (int i = 0; i < nnz; ++i) {

      T[(j)*nnz+i] = complex<double>(0.0,0.0);

      for (int k = 0; k < nt; ++k) { 
        T[(j)*nnz+i] += conj(What[(i)*nt+k]) * R[(j)*nt+k];
      }
    }
  }

  cout << " ok. " << endl;

  delete[] R;

  // now compute the rhs vector ..  r = What^* a_{n-1}

  cout << " Setting up rhs ... ";

  complex<double> * rhs = new complex<double>[nnz];

  for (int i = 0; i < nnz; ++i) { 

    rhs[i] = complex<double>(0.0,0.0);
    
    for (int j = 0; j < nt; ++j) { 
      rhs[i] += conj(What[(i)*nt+j]) * complex<double>(A_corr[(nt)*ns+(j+1)]);
    }
  }

  cout << " ok. " << endl;

  // solve T x = r for the adjoint amplitudes .. 

  cout << " Solving for complex amp ... " ;

  { 
    int nnz_tmp = nnz;
    int ione    = 1;
    int * ipiv  = new int[nnz];
    int info;

    zgesv_(&nnz_tmp,&ione,T,&nnz_tmp,ipiv,rhs,&nnz_tmp,&info);

    delete[] ipiv;

    if ( info != 0 ) 
      cout << " error in the solve for amplitudes " << endl;

  } 

  cout << " ok. " << endl;

  // set the amplitudes .. 

  assert( amp == NULL); amp = new double[nnz];

  for (int i = 0; i < nnz; ++i) { 

    amp[i] = sqrt( rhs[i].real()*rhs[i].real() + rhs[i].imag()*rhs[i].imag());

  }


  delete[] What;
  delete[] T;
  delete[] rhs;
}

void compute_dmd_adjoint_amp(double* &adj_amp, const int ns, const double* A_corr, 
                             const double* Wk, const double* lambda_r_adj, 
                             const double* lambda_i_adj, complex<double>* vrh_imag, 
                             const double dt_samp, const int nnz) {

  // recall V_tilde = [nt X nnz]

  int nt = ns - 1;

  // compute Phat = V_tilde P = [nt X nnz] X [nnz X nnz] 

  cout << " Computing Phat .. " ;

  complex<double>* Phat = new complex<double>[nt * nnz];

  for (int j = 0; j < nnz; ++j) { 
    for (int i = 0; i < nt; ++i) { 

      Phat[(j)*nt+i] = complex<double>(0.0,0.0);

      for (int k = 0; k < nnz; ++k) { 
        Phat[(j)*nt+i] += complex<double>(Wk[(k)*nt+i])*vrh_imag[(j)*nnz+k];
      }
    }
  }

  cout << " ok. " << endl;

  // need to form T = Phat^* A Phat .. = Phat^* R
  // R = [nt X nt] X [nt X nnz] = [nt X nnz]
  
  cout << " Computing R ... " ;

  complex<double>* R = new complex<double>[nt * nnz];

  for (int j = 0; j < nnz; ++j) { 
    for (int i = 0; i < nt; ++i) { 

      R[(j)*nt+i] = complex<double>(0.0,0.0);

      for (int k = 0; k < nt; ++k) { 
        R[(j)*nt+i] += complex<double>(A_corr[(k)*ns+i])*Phat[(j)*nt+k];
      }
    }
  }

  cout << " ok." << endl;

  // T = Phat^* R = [nnz X nt] X [nt X nnz] = [nnz X nnz]

  cout << " Computing T... ";

  complex<double>* T = new complex<double>[nnz * nnz];

  for (int j = 0; j < nnz; ++j) { 
    for (int i = 0; i < nnz; ++i) {

      T[(j)*nnz+i] = complex<double>(0.0,0.0);

      for (int k = 0; k < nt; ++k) { 
        T[(j)*nnz+i] += conj(Phat[(i)*nt+k]) * R[(j)*nt+k];
      }
    }
  }

  cout << " ok. " << endl;

  delete[] R;

  // now compute the rhs vector ..  r = Phat^* a_{n-1}

  cout << " Setting up rhs ... ";

  complex<double> * rhs = new complex<double>[nnz];

  for (int i = 0; i < nnz; ++i) { 

    rhs[i] = complex<double>(0.0,0.0);
    
    for (int j = 0; j < nt; ++j) { 
      rhs[i] += conj(Phat[(i)*nt+j]) * complex<double>(A_corr[(nt-1)*ns+j]);
    }
  }

  cout << " ok. " << endl;

  // solve T x = r for the adjoint amplitudes .. 

  cout << " Solving for complex amp ... " ;

  { 
    int nnz_tmp = nnz;
    int ione    = 1;
    int * ipiv  = new int[nnz];
    int info;

    zgesv_(&nnz_tmp,&ione,T,&nnz_tmp,ipiv,rhs,&nnz_tmp,&info);

    delete[] ipiv;

    if ( info != 0 ) 
      cout << " error in the solve for amplitudes " << endl;

  } 

  cout << " ok. " << endl;

  // set the amplitudes .. 

  assert( adj_amp == NULL); adj_amp = new double[nnz];

  for (int i = 0; i < nnz; ++i) { 

    adj_amp[i] = sqrt( rhs[i].real()*rhs[i].real() + rhs[i].imag()*rhs[i].imag());

  }


  delete[] Phat;
  delete[] T;
  delete[] rhs;
}

#endif
