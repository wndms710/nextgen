

#include <iostream>
#include <math.h>
#include <assert.h>
#include <cstdlib>
#include "RkWeights.hpp"


using namespace Ascher343;
using namespace std;

int main(const int argc, const char* argv[]) { 

  const double u0   = 1.0;
  double u          = u0;
  //double dt         = 1.0e-2;

  double dt         = atof(argv[1]);

  // du/dt = lam1*u + lam2*u;

  const double lam1 = -1.0;
  const double lam2 = -10.0;

  double time          = 0.0;
  const double t_final = 1.0;

  double k[3];
  double khat[4];

  while ( time < t_final) { 

    
    //===================
    // rk stage 1
    //===================

    { 

      const double dd = 1.0 - dt*A[0][0]*lam2;
      khat[0]         = lam1*u;

      u              += dt*A_hat[1][0]*khat[0];
      u              /= dd;

      k[0]            = lam2*u;

    }


    //===================
    // rk stage 2
    //===================

    { 

      const double dd = 1.0 - dt*A[1][1]*lam2;
      khat[1]         = lam1*u;

      u              += dt*(A[1][0] - A[0][0])*k[0];
      u              += dt*(A_hat[2][0] - A_hat[1][0])*khat[0] + dt*A_hat[2][1]*khat[1];
      u              /= dd;

      k[1]            = lam2*u;

    }


    //==================
    // rk stage 3 
    //==================

    { 

      const double dd = 1.0 - dt*A[2][2]*lam2;
      khat[2]         = lam1*u;

      u              += dt*(A[2][0] - A[1][0])*k[0] + dt*(A[2][1] - A[1][1])*k[1];
      u              += dt*(A_hat[3][0] - A_hat[2][0])*khat[0] + dt*(A_hat[3][1] - A_hat[2][1])*khat[1] + dt*A_hat[3][2]*khat[2];
      u              /= dd;

      k[2]            = lam2*u;


    }


    //===============
    // rk stage 4 
    //===============

    { 

      khat[3]         = lam1*u;

      for (int j = 0; j < 3; ++j) {  
        u  += dt*(bt[j] - A[2][j])*k[j];            
        u  += dt*(bt_hat[j] - A_hat[3][j])*khat[j];
      }

      u    += dt*bt_hat[3]*khat[3];

    }


    time += dt;

    const double exact_u = u0*exp((lam1+lam2)*time);
    const double err     = fabs(u - exact_u);

    cout << "  XXXX " << time << "    " << u << "    " << err << endl;

  }


  return 0;

}


