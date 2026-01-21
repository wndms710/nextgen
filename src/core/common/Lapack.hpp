#ifndef LAPACK_HPP 
#define LAPACK_HPP 

#include <complex> 

extern "C" void dgeev_(char*, char*, int*, double*, int*, double*, double*, double*, int*, 
		  double*, int*, double*, int*, int*) ; 


extern "C" void zgeev_(char*, char*, int*, std::complex<double>*, int*, std::complex<double>*,
                       std::complex<double>*,int*, std::complex<double>*, int*, std::complex<double>*, 
                       int*, double*, int*);

extern "C" void dsysv_(char*, int*, int*, double*, int*, int*, double*, int*, double*, 
		       int*, int* ) ; 


extern "C" void zsysv_(char*, int*, int*, std::complex<double>*, int*, int*, std::complex<double>*, int*, 
                       std::complex<double>*, int*, int*);

extern "C" void dsyev_(char * , char * , int * , double * , int * , double *, double *, int*, int* ) ; 


extern "C" void zgetrf_(int* , int * , std::complex<double> * , int * , int*, int*) ; 
extern "C" void zgetri_(int* , std::complex<double>* , int*, int*, std::complex<double>*, int* , int * ) ; 

extern "C" void zgesv_(int* , int*, std::complex<double>*, int*, int*, std::complex<double>*, int*, int*); 

extern "C" void zgemm_(char*, char*, int*, int*, int*, std::complex<double>*, std::complex<double>*, int*, 
                       std::complex<double>*, int*, std::complex<double>*, std::complex<double>*, int*); 

extern "C" void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*,
                       double*, int*, double*, double*, int*);

extern "C" void zgelss_(int*, int*, int*, std::complex<double>*, int*, std::complex<double>*, int*, 
                        double*, double*, int*, std::complex<double>*,int*,double*,int*);

extern "C" void dgetrf_(int*, int*, double*, int*, int*, int*);

extern "C" void dpotrf_(char*, int*, double*, int*, int*);

extern "C" void dpotrs_(char*, int*, int*, double*, int*, double*, int*, int*);

extern "C" void dgetri_(int*, double*, int*, int*, double*, int*, int*);

extern "C" void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, 
                        double*, int*, int*);

extern "C" void dgeqp3_(int*,int*,double*,int*,int*,double*,double*,int*,int*);

extern "C" void dormqr_(char*,char*,int*,int*,int*,double*,int*,double*,double*,int*,double*,int*,int*);
   
#endif 
