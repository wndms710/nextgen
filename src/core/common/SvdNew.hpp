#ifndef SVDNEW_HPP
#define SVDNEW_HPP

#include <cmath>
#include <iostream>
#include <assert.h>

#define SIGN(x) ((x > 0.0) ? 1.0 : -1.0)

using namespace std;

/*
double SIGN(double x) { 
  if ( x > 0.0) {
    return 1.0;
  } else if ( x < 0.0) { 
    return -1.0;
  } else { 
    return 1.0;
  }
}
*/

inline int calcSvd33(double * U, double* V, double* sigma, const double* A, double& max_err) { 

  // assumes that the matrices are being provided in column major format ... 
  const double eps_tol = 1.0e-12;

  for (int i =0; i < 9; ++i) 
    U[i] = A[i];

  for (int j = 0; j < 3; ++j) { 
    for (int i =0; i < 3; ++i) { 
      if ( i ==j ) { 
        V[(j)*3+i] = 1.0;
      } else { 
        V[(j)*3+i] = 0.0;
      }
    }
  }

  int done = 0;
  int iter = 0;
  while ( done == 0 ) { 

    done = 1; // assume we're finished...
    ++iter;

    max_err = 0.0;

    for (int j = 0; j < 3; ++j) { 
      for (int i = 0; i < j; ++i) { 

        double a = 0.0;
        double b = 0.0;
        double c = 0.0;

        for (int k =0; k < 3; ++k) { 
          a += U[(i)*3+k]*U[(i)*3+k];
          b += U[(j)*3+k]*U[(j)*3+k];
          c += U[(i)*3+k]*U[(j)*3+k];
        }

        //const double zeta = (b-a)/(2.0*c);  
        double zeta = 1.0e+20; 
        if ( c != 0.0) { 
          zeta = (b-a)/(2.0*c); 
        }
        
        const double t    = SIGN(zeta)/(fabs(zeta) + sqrt(1.0+zeta*zeta));
        const double cs   = 1.0/sqrt(1.0+t*t);
        const double sn   = cs*t;

        for (int k =0; k < 3; ++k) { 
          const double tmp = U[(i)*3+k];
          U[(i)*3+k]       = cs*tmp - sn*U[(j)*3+k];
          U[(j)*3+k]       = sn*tmp + cs*U[(j)*3+k];
        }

        for (int k = 0; k < 3; ++k) { 
          const double tmp = V[(i)*3+k];
          V[(i)*3+k]       = cs*tmp - sn*V[(j)*3+k];
          V[(j)*3+k]       = sn*tmp + cs*V[(j)*3+k];
        }

        if ( (fabs(c) > eps_tol*sqrt(a*b)) && (fabs(c) > eps_tol*eps_tol)) {
          if ( sqrt(a*b) != 0.0) { 
            max_err = fmax(max_err,fabs(c)/sqrt(a*b));
          } else { 
            //cout << " fabs(c) : " << fabs(c) << endl;
            if ( a > 0.0) { 
              max_err = fmax(max_err,fabs(c)/a);
            } else if ( b > 0.0) { 
              max_err = fmax(max_err,fabs(c)/b);
            } else { 
              max_err = fmax(max_err,fabs(c));
            }
          }
              
          done = 0;
        }
      }
    }

    if ( iter > 100)  {
      //cout <<  "Svd33:: max_error, tol: " << max_err << "   " << eps_tol << endl;
      break;
    }

  }

  for (int j = 0; j < 3; ++j) { 
    
    sigma[j] = 0.0;
    for (int i = 0; i < 3; ++i) { 
      sigma[j] += U[(j)*3+i]*U[(j)*3+i];
    }

    sigma[j] = sqrt(sigma[j]);

    if ( sigma[j] > 0.0) { 
      for (int i =0; i < 3; ++i) 
        U[(j)*3+i] /= sigma[j];
    }
  }

  return iter;
} 

inline int calcSvd44(double * U, double* V, double* sigma, const double* A, double& max_err) { 

  // assumes that the matrices are being provided in column major format ... 
  const double eps_tol = 1.0e-12;

  for (int i =0; i < 16; ++i) 
    U[i] = A[i];

  for (int j = 0; j < 4; ++j) { 
    for (int i =0; i < 4; ++i) { 
      if ( i ==j ) { 
        V[(j)*4+i] = 1.0;
      } else { 
        V[(j)*4+i] = 0.0;
      }
    }
  }

  int done = 0;
  int iter = 0;
  while ( done == 0 ) { 

    done = 1; // assume we're finished...
    ++iter;

    max_err = 0.0;

    for (int j = 0; j < 4; ++j) { 
      for (int i = 0; i < j; ++i) { 

        double a = 0.0;
        double b = 0.0;
        double c = 0.0;

        for (int k =0; k < 4; ++k) { 
          a += U[(i)*4+k]*U[(i)*4+k];
          b += U[(j)*4+k]*U[(j)*4+k];
          c += U[(i)*4+k]*U[(j)*4+k];
        }

        //const double zeta = (b-a)/(2.0*c);  
        double zeta = 1.0e+20; 
        if ( c != 0.0) { 
          zeta = (b-a)/(2.0*c); 
        }
        
        const double t    = SIGN(zeta)/(fabs(zeta) + sqrt(1.0+zeta*zeta));
        const double cs   = 1.0/sqrt(1.0+t*t);
        const double sn   = cs*t;

        for (int k =0; k < 4; ++k) { 
          const double tmp = U[(i)*4+k];
          U[(i)*4+k]       = cs*tmp - sn*U[(j)*4+k];
          U[(j)*4+k]       = sn*tmp + cs*U[(j)*4+k];
        }

        for (int k = 0; k < 4; ++k) { 
          const double tmp = V[(i)*4+k];
          V[(i)*4+k]       = cs*tmp - sn*V[(j)*4+k];
          V[(j)*4+k]       = sn*tmp + cs*V[(j)*4+k];
        }

        if ( (fabs(c) > eps_tol*sqrt(a*b)) && (fabs(c) > eps_tol*eps_tol)) {
          if ( sqrt(a*b) != 0.0) { 
            max_err = fmax(max_err,fabs(c)/sqrt(a*b));
          } else { 
            //cout << " fabs(c) : " << fabs(c) << endl;
            if ( a > 0.0) { 
              max_err = fmax(max_err,fabs(c)/a);
            } else if ( b > 0.0) { 
              max_err = fmax(max_err,fabs(c)/b);
            } else { 
              max_err = fmax(max_err,fabs(c));
            }
          }
              
          done = 0;
        }
      }
    }

    if ( iter > 100)  {
      //cout <<  "Svd44:: max_error, tol: " << max_err << "   " << eps_tol << endl;
      break;
    }

  }

  for (int j = 0; j < 4; ++j) { 
    
    sigma[j] = 0.0;
    for (int i = 0; i < 4; ++i) { 
      sigma[j] += U[(j)*4+i]*U[(j)*4+i];
    }

    sigma[j] = sqrt(sigma[j]);

    if ( sigma[j] > 0.0) { 
      for (int i =0; i < 4; ++i) 
        U[(j)*4+i] /= sigma[j];
    }
  }

  return iter;
} 

#endif
