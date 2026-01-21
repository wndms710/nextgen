#ifndef DMD_HPP
#define DMD_HPP

#define corr(i,j) corr[(j)*ldc +i]
#define cij(i,j) cij[(j)*n_snaps + i]
#define vr(i,j) vr[(j)*nnz+i]
#define vr_imag(i,j) vr_imag[(j)*nnz+i]

#define Wk(i,j) Wk[(j)*n_snaps+i]
#define Hcat(i,j) Hcat[(j)*n_snaps+i]
#define tmp(i,j) tmp[(j)*n_snaps+i]
#define Abar(i,j) Abar[(j)*nnz+i]


#include <iostream>
#include <assert.h>
#include <math.h>
#include <iomanip>
#include <complex>
#include <cstdio>
#include "Lapack.hpp"

using namespace std ;

void dumpMatrix(const int n, double * A) {
  for (int i = 0; i < n ; ++i) {
    for (int j = 0 ; j < n ; ++j)
    cout << A[(j)*n + i] << "   " ;
    cout << endl ;
  }
}//dumpMatrix()

void checkLapackReturn(const int info, const char * msg) {
  if ( info != 0 ) {
    cout << msg << " [failed ]; code = " << info << endl;
    throw(0);
  }
}//checkLapackReturn()

static inline int dmdAnalysis( const int n , const double *corr, double* &lambda_r, double* &lambda_i,
                double *& Wk, complex<double> * &vr_imag, double * &amp, const double delta_T) {
  
  //   n   : total number of snapshots
  //   corr: correlation matrix
  assert ( lambda_r == NULL );
  assert ( lambda_i == NULL );
  assert ( Wk       == NULL );
  assert ( vr_imag  == NULL );
  
  int ldc   = n ;
  int n_snaps    = n-1;
  int nnz   = 0   ;
  char jobz, uplo, jobvr, jobvl;
  int lwork,info;
  double * work = NULL;
  
  // extract n-1 submatrix from the full correlation matrix.
  // compute its eigenvalues, eigenvectors.
  double * sigma = new double[n_snaps];
  Wk             = new double[n_snaps*n_snaps] ;
  
  for (int j =0; j < n_snaps ; ++j) {
    for (int i = 0; i < n_snaps ; ++i)
    Wk(i,j) = corr(i,j); // this is H in the notes...
  }
  
  jobz     = 'v' ;
  uplo     = 'u' ;
  lwork    = 5* n_snaps;
  work     = new double[lwork];
  
  cout << " > Computing eigs of correlation matrix ... " ; cout.flush();
  dsyev_(&jobz, &uplo, &n_snaps, Wk, &n_snaps, sigma, work, &lwork, &info);
  checkLapackReturn(info,"compute eigs of correlation submatrix");
  cout << "OK" << endl;
  
  delete[] work; work = NULL;
  
  // the singular values are in ascending order...
  const double sigma_max      = sigma[n_snaps-1] ; assert(sigma_max > 0.0);
  //const double rank_reduction = 0.0;
  
  // count number of non-zero singular values..
  nnz = 0 ;
  for (int i =0; i < n_snaps ; ++i) {
    if ( (fabs(sigma[i]) > 1.0e-06)) { // assuming data is scaled appropriately on entry...
      assert(sigma[i] >= 0.0);
      //cout << nnz << "   " << sigma[i] << endl;
      ++nnz;
    }
  }//i
  cout << " > Rank = " << nnz << endl ;
  
  if ( nnz < n_snaps ) {
    
    // arrange the eigenvalues and eigenvectors in descending order.
    int shft = n_snaps -nnz ;
    for (int j=n_snaps-nnz; j < n_snaps ; ++j) {
      sigma[j-shft] = sigma[j] ;
      for (int i = 0 ; i < n_snaps ; ++i)
      Wk(i,j-shft) = Wk(i,j);
    }//j
  }
  
  // scale columns of Wk (eigenvectors) by sigma..
  for (int j = 0 ; j < nnz ; ++j) {
    for (int i = 0 ; i < n_snaps ; ++i )
    Wk(i,j) *= 1.0/sqrt(sigma[j]) ; // python code P*Sigg
  }//j
  
  double * Hcat = new double[n_snaps * n_snaps];
  
  for (int j =1; j < n_snaps ; ++j)
  for (int i =0; i < n_snaps ; ++i)
  Hcat(i,j-1) = corr(i,j);
  
  for (int i =0; i < n_snaps ; ++i)
  Hcat(i,n_snaps-1) = corr(i,n-1);
  
  // Abar <--- Wk* * tmp * Wk
  double * tmp   = new double[n_snaps * nnz ] ;
  double * Abar  = new double[nnz* nnz ] ; // Abar is M in the test codes...
  double * eig_r = new double[nnz];
  double * eig_i = new double[nnz];
  double * vr    = new double[nnz * nnz  ];
  amp            = new double[nnz] ;
  
  for (int i = 0 ; i < n_snaps ; ++i ) {
    for (int j = 0 ; j < nnz ; ++j ) {
      tmp(i,j) = 0.0 ;
      for (int k = 0 ; k < n_snaps ; ++k)
      tmp(i,j) += Hcat(i,k) * Wk(k,j) ;
    }
  }
  
  for (int j=0 ; j < nnz ; ++j) {
    for (int i = 0 ; i < nnz ; ++i) {
      Abar(i,j) = 0.0 ;
      for (int k = 0 ; k < n_snaps ; ++k)
      Abar(i,j) += Wk(k,i) * tmp(k,j) ;
    }//i
  }//j
  
  delete[] Hcat;
  
  
  // compute the complex eigenvalues and eigenvectors...
  jobvl = 'n';
  jobvr = 'v';
  lwork = 20*nnz ;
  work  = new double[lwork];
  
  double vl ;
  int ldvl =1, ldvr = nnz ;
  
  cout << " > Computing dmd eigs ... " ; cout.flush();
  dgeev_( &jobvl, &jobvr, &nnz, Abar, &nnz, eig_r, eig_i, &vl, &ldvl,
         vr, &ldvr, work, &lwork, &info ) ;
  
  checkLapackReturn(info, "dmd eigs");
  cout << " OK." << endl;
  
  delete[] work; work = NULL;
  
  // frequencies, growth rates...
  lambda_r    = new double[nnz];
  lambda_i    = new double[nnz]; // log(eig) scaling ..
  
  // computing the amplitudes...
  vr_imag = new complex<double>[nnz*nnz];
  for (int j= 0; j < nnz ; ++j) {
    if ( fabs(eig_i[j]) < 1.0e-06) {
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
  
  // compute the dprime first...
  complex<double>* di = new complex<double>[nnz];
  for (int i =0; i < nnz ; ++i) {
    di[i] = complex<double>(0.0,0.0);
    for (int j=0; j < nnz ; ++j) {
      
      double sum = 0.0;
      for (int m=0; m < n_snaps; ++m) {
        sum += Wk(m,j) * corr(m,0);
      }
      
      di[i]  +=  complex<double>(sum,0.0) * vr_imag(j,i);
    }
    di[i] = conj(di[i]);
  }
  
  complex<double>* VstarV= new complex<double>[nnz*nnz];
  {
    char transa = 'C';
    char transb = 'N';
    complex<double> alfa(1.0,0.0);
    complex<double> beta(0.0,0.0);
    zgemm_(&transa, &transb, &nnz, &nnz, &nnz, &alfa, vr_imag, &nnz,
           vr_imag, &nnz,
           &beta, VstarV , &nnz);
    
    
    
    // solve (V*V) d = dprime...
    int ione = 1;
    int * ipiv = new int[nnz];
    zgesv_(&nnz, &ione, VstarV, &nnz, ipiv, di, &nnz, &info);
    checkLapackReturn(info, "VstarV inversion");
    delete[] ipiv;
  }
  
  cout << " > Amplitude setup ... " ; cout.flush();
  double * Tkm = new double[nnz*n_snaps];
  
#if 0
  for (int k = 0; k < nnz ; ++k) {
    for (int m =0; m < ns ; ++m ) {
      
      double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
      for (int r = 0; r < n_snaps ; ++r) {
        sum += Wk(r,k) * corr(r,m);
      }
      
      Tkm[(k)*n_snaps+m] = sum;
    }
  }
#else
  {
    double * corr_small = new double[n_snaps*n_snaps];
    
    //#pragma omp parallel for
    for (int k = 0; k < n_snaps ; ++k) {
      for (int m =0; m < n_snaps ; ++m) {
        corr_small[m+k*n_snaps] = corr[m+k*n];
      }
    }
    
    char transa = 'N';
    char transb = 'N';
    double alpha=1.0;
    double beta =0.0;
    
    //this one does T = Corr x W
    dgemm_(&transa, &transb,
           &n_snaps, &nnz, &n_snaps,
           &alpha, &corr_small[0], &n_snaps,
           &Wk[0], &n_snaps,
           &beta, &Tkm[0], &n_snaps);
    
    delete[] corr_small;
  }
#endif
  
  cout << " ... " ; cout.flush();
  
  
  double * Qjk = new double[nnz*nnz];
#if 0
  for (int j =0; j < nnz ; ++j) {
    for (int k =0; k < nnz ; ++k) {
      
      double sum = 0.0;
# pragma omp parallel for reduction(+:sum)
      for (int m =0; m < n_snaps ; ++m) {
        sum += Wk(m,j) * Tkm[(k)*n_snaps+m];
      }
      Qjk[(j)*nnz+k] = sum ;
    }
  }
#else
  char transa = 'N';
  char transb = 'N';
  double alpha=1.0;
  double beta =0.0;
  
  double * Wk_t = new double[n_snaps*nnz];
  
  //#pragma omp parallel for
  for(int n = 0; n<nnz*n_snaps; n++) {
    int i = n/nnz;
    int j = n%nnz;
    Wk_t[n] = Wk[i + j*n_snaps];
  }
  
  //this does Q = W^t x T (W is transposed manually above cause blas SUCKS!)
  dgemm_(&transa, &transb,
         &nnz, &nnz, &n_snaps,
         &alpha, &Wk_t[0], &nnz,
         &Tkm[0], &n_snaps,
         &beta, &Qjk[0], &nnz);
  
  delete[] Wk_t;
#endif
  
  cout << "OK" << endl;
  cout.flush();
  delete[] Tkm;
  
  // and lastly .. the amplitude ..
  for (int i = 0 ; i < nnz ; ++i) {
    
    complex<double> sum(0.0,0.0);
    
    for (int j=0; j < nnz ; ++j) {
      for (int k =0; k < nnz ; ++k) {
        complex<double> vr_conj(vr_imag(k,i).real(), -vr_imag(k,i).imag());
        sum += vr_imag(j,i)* vr_conj * complex<double>( Qjk[(j)*nnz+k], 0.0);
      }
    }
    
    //const double norm = sqrt( sum.real()*sum.real() + sum.imag() + sum.imag());
    //assert ( fabs(norm - 1.0) < 1.0e-10);
    
    if ( fabs(sum.imag()) >= 1.0e-09) { 
      cout << " non-real value encountered : " << sum.imag() << endl;
    } else if (sum.real() < 0.0) { 
      cout << " negative sum encountered : " << sum.real() << endl;
    }
    cout.flush();
    
    //assert ( (fabs(sum.imag()) < 1.0e-09) && (sum.real() >= 0.0)) ;
    
    complex<double> di_conj(di[i].real(), -di[i].imag());
    sum *= di[i] * di_conj ;
    amp[i] =  sqrt(sum.real());
  }
  
  delete[] Qjk;
  
  // show me
  for (int i = 0 ; i < nnz ; ++i)  {
    complex<double> tt(eig_r[i], eig_i[i]) ;
    double f     = arg(tt);
    double sigma = log(abs(tt));
    
    lambda_i[i]  = f    / delta_T / 2.0/ M_PI ;
    lambda_r[i]  = sigma/ delta_T ;
    
    printf(" > Freq, gr, amplitude =  %12.8g   %12.8g  %12.8g\n", lambda_i[i],
           lambda_r[i],
           amp[i]);
    
  }//i
  
  delete[] VstarV;
  delete[] di;
  delete[] vr ;
  delete[] Abar ;
  delete[] sigma;
  delete[] tmp ;
  delete[] eig_r ;
  delete[] eig_i ;
  
  return nnz;
}//dmdAnalysis()


// reconstruction coefficient for mode i from snapshot m...
void reconstructCoeff(complex<double>& coeff, const double * Wk, const complex<double>* vr_imag, const int nnz,
                      const int n, const int i, const int m) {
  
  const int n_snaps = n-1;
  
  assert ( (m >= 0) && (m < n_snaps)) ;
  assert ( (i >= 0 ) && (i < nnz )) ;
  
  coeff = complex<double>(0.0,0.0);
  for (int j = 0; j < nnz ; ++j) {
    coeff += complex<double>(Wk(m,j)) * vr_imag(j,i) ;
  }
}//reconstructCoeff()


#endif
