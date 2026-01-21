// ================================================================
// fft interface.  can be built with fftw or the functions below 
// need to be provided with implementations if using a different fft
// ================================================================

#ifndef FFT_HPP
#define FFT_HPP

#ifdef WITH_FFTW

#include <fftw3.h>

#else
#define FFT_ERR_MSG std::cerr << "Error: compile with -DWITH_FFTW or implement FFT routines" << std::endl; \
  throw(0); 

typedef int fftw_plan;
typedef int fftw_complex;

enum fftw_consts {
  FFTW_ESTIMATE,
  FFTW_MEASURE
};

void fftw_destroy_plan(fftw_plan plan) { FFT_ERR_MSG; } 

fftw_plan fftw_plan_dft_c2r_1d(const int n,fftw_complex * data,
                               double * d2,const int flag) { FFT_ERR_MSG; } 

fftw_plan fftw_plan_dft_c2r_2d(const int n1,const int n2,
			       fftw_complex * data,double * d2,const int flag) { FFT_ERR_MSG; } 

fftw_plan fftw_plan_many_dft_r2c(int rank, const int* n, int howmany, 
				 double *in, const int* inembed,
				 int istride, int idist, 
				 fftw_complex *out, const int* onembed,
				 int ostride, int odist, 
				 unsigned flags) { FFT_ERR_MSG; } 

fftw_plan fftw_plan_dft_r2c_1d(const int n,
			       double * d2,fftw_complex * data,const int flag) { FFT_ERR_MSG;} 

fftw_plan fftw_plan_dft_r2c_2d(const int n1,const int n2,
			       double * d2,fftw_complex * data,const int flag) { FFT_ERR_MSG;}

void fftw_execute(fftw_plan plan) {FFT_ERR_MSG;} 

#endif // #ifdef WITH_FFTW
#endif // #ifndef FFT_HPP

