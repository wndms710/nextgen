#ifndef POD_HPP
#define POD_HPP

#ifndef corr
#define corr(i,j) corr[(j)*ldc +i]
#endif

#ifndef cij
#define cij(i,j) cij[(j)*ns + i]
#endif

#ifndef vr
#define vr(i,j) vr[(j)*nnz+i]
#endif

#ifndef Wk
#define Wk(i,j) Wk[(j)*ns+i]
#endif

#ifndef Hcat
#define Hcat(i,j) Hcat[(j)*ns+i]
#endif

#ifndef tmp
#define tmp(i,j) tmp[(j)*ns+i]
#endif

#ifndef Abar
#define Abar(i,j) Abar[(j)*nnz+i]
#endif


#include <iostream>
#include <assert.h>
#include <math.h>
#include <iomanip>
#include <complex>
#include <cstdio>
#include "Lapack.hpp"

using namespace std ;


int podAnalysis(const int n , const double *corr, double* &sigma, double *& Wk) {
  
  //   n   : total number of snapshots
  //   corr: correlation matrix
  assert ( sigma    == NULL ); 
  assert ( Wk       == NULL );
  
  int ldc   = n ;
  int ns    = n-1;
  int nnz   = 0   ;
  char jobz, uplo; 
  int lwork,info;
  double * work = NULL;
  
  // extract n-1 submatrix from the full correlation matrix.
  // compute its eigenvalues, eigenvectors.
  sigma          = new double[ns];
  Wk             = new double[ns*ns] ;
  
  for (int j =0; j < ns ; ++j) {
    for (int i = 0; i < ns ; ++i)
      Wk(i,j) = corr(i,j); // this is H in the notes...
  }
  
  jobz     = 'v' ;
  uplo     = 'u' ;
  lwork    = 5* ns;
  work     = new double[lwork];
  
  cout << " > Computing eigs of correlation matrix ... " ; cout.flush();
  dsyev_(&jobz, &uplo, &ns, Wk, &ns, sigma, work, &lwork, &info);
  if ( info != 0 ) 
    cout << "compute eigs of correlation submatrix failed : " << info << endl; 
  cout << "OK" << endl;
  
  delete[] work; work = NULL;
  
  // count number of non-zero singular values..
  nnz = 0 ;
  for (int i =0; i < ns ; ++i) {
    cout << "sigma[" << i << "] = " << sigma[i] << endl;
    if ((fabs(sigma[i]) > 1.0e-06)) { // assuming data is scaled appropriately on entry...
      //assert(sigma[i] >= 0.0);
      //cout << nnz << "   " << sigma[i] << endl;
      ++nnz;
    }
  }//i
  cout << " > Rank = " << nnz << endl ;
  
  if ( nnz < ns ) {
    
    // arrange the eigenvalues and eigenvectors in descending order.
    int shft = ns -nnz ;
    for (int j=ns-nnz; j < ns ; ++j) {
      sigma[j-shft] = sigma[j] ;
      for (int i = 0 ; i < ns ; ++i)
      Wk(i,j-shft) = Wk(i,j);
    }//j
  }
  
  // scale columns of Wk (eigenvectors) by sigma..
  for (int j = 0 ; j < nnz ; ++j) {
    for (int i = 0 ; i < ns ; ++i )
    Wk(i,j) *= 1.0/sqrt(sigma[j]) ; // python code P*Sigg
  }//j
 

  // show me the cumulative energy from the POD reconstruction...
  double sum = 0.0; 
  for (int j =0; j < nnz ; ++j) 
    sum += sigma[j]; 

  FILE * fp = fopen("pod_amps.dat", "w"); 
  cout << " ---- POD normalized modal energy decomp ---- " << endl;

  fprintf(fp, "#Mode #\tMarginal      \tAggregate\n");
  cout <<     "Mode #\tMarginal\tAggregate" << endl;

  double cumsum = 0.0; 
  for (int j = nnz -1; j >= 0; --j) { 
    cumsum += sigma[j];
    cout << nnz-1-j << "\t" << sigma[j]/sum << "\t" << cumsum/sum << endl;
    fprintf(fp, "%7i\t%12.8e\t%12.8e\n", nnz-1-j, sigma[j]/sum, cumsum/sum);
  }//(j < nnz - 1)

  fclose(fp);

  return nnz;

}//podAnalysis
    
#endif
