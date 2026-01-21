
#include <iostream>
#include <cstdlib>
#include "SvdNew.hpp"
#include <assert.h>
#include <cstdio>

using namespace std;

#ifdef OLD
int main() { 

  /*
  double A[9];
  for (int ii = 0; ii < 9; ++ii) 
    A[ii] = 0.0;

  for (int i =0; i < 3; ++i) 
    A[(i)*3+i] = 1.0;

  A[(2)*3+0] = 8.89845e-18;
  A[(2)*3+1] = -8.31061e-18;

  double U[9];
  double V[9];
  double sigma[3];
  int ii = calcSvd33(U,V,sigma,A);

  for (int ii = 0; ii < 9; ++ii)
    cout << U[ii] << "   " << V[ii] << endl;
  cout << sigma[0] << "   " << sigma[1] << "    " << sigma[2] << endl;
  */

  const int m = 500;
 
  double global_max_err = 0.0;

  for (int iter =0; iter < m; ++iter) { 

    double A[9];
    double U[9];
    double V[9];
    double sigma[3];

    /// populate a O(1) random matrix ... 
    for (int i =0; i < 9; ++i) 
      A[i] = (2.0*double(rand())/double(RAND_MAX)-1.0);

    int ii = calcSvd33(U,V,sigma,A);

    // now compute the error in USV^T - A.. 
    // V <--- V\Sigma.. multiply the columns

    for (int j = 0; j < 3; ++j) { 
      for (int i =0; i < 3; ++i) 
        V[(j)*3+i] *= sigma[j];
    }

    double max_err = 0.0;
    double A_nrm   = 0.0;
    for (int j = 0; j < 3; ++j) { 
      for (int i =0; i < 3; ++i) { 

        double aij = 0.0;
        for (int k =0; k < 3; ++k) 
          aij += U[(k)*3+i]*V[(k)*3+j];

        max_err = max(max_err,abs(aij-A[(j)*3+i]));
        A_nrm   = A[(j)*3+i]*A[(j)*3+i];
      }
    }

    A_nrm = sqrt(A_nrm);
    cout << " iter, max_err, normF(A) = " << iter << "    " << max_err << "   " << A_nrm << "   " << ii << endl;
    cout << " iter, sigma vlue: " << iter << "   " << sigma[0] << "   " << sigma[1] << "    " << sigma[2] << endl;
    assert ( sigma[0] >= 0.0);
    assert ( sigma[1] >= 0.0);
    assert ( sigma[2] >= 0.0);
    global_max_err = max(global_max_err,max_err);
  }

  cout << " ========= summary ========== " << endl;
  cout << " max global error : " << global_max_err << endl;

  return 0;
}
#endif 


int main() { 

  double A[9];
  double U[9];
  double V[9];
  double sigma[3];
  
  //FILE * fp = fopen("broken_svd/A.932676.dat","r");
  /*
  FILE * fp = fopen("broken_svd/A.041578.dat", "r");

  for (int i =0; i < 3; ++i) {  
    fscanf(fp,"%lf  %lf  %lf\n", &A[(0)*3+i], &A[(1)*3+i], &A[(2)*3+i]);
    cout << A[(0)*3+i] << "    " <<  A[(1)*3+i] << "   " <<  A[(2)*3+i] << endl;
  }

  fclose(fp);
  

  double eps_tol;
  int ii = calcSvd33(U,V,sigma,A,eps_tol);
  cout << " Sing value, eps_tol : " << sigma[0] << "   " << sigma[1] << "    " << sigma[2] << "   " << eps_tol << endl;
  */

  FILE * fp2 = fopen("A.000568.binary.dat","rb");
  fread(A,sizeof(double),9,fp2);
  fclose(fp2);

  cout << " ==== A === " << endl;
  for (int i =0; i < 3; ++i) { 
    cout << A[(0)*3+i] << "   " << A[(1)*3+i] << "    " << A[(2)*3+i] << endl;
  }


  double eps_tol;
  int ii = calcSvd33(U,V,sigma,A,eps_tol);

  cout << " Sing value, eps_tol : " << sigma[0] << "   " << sigma[1] << "    " << sigma[2] << "   " << eps_tol << endl;

  cout << " ==== U === " << endl;
  for (int i =0; i < 3; ++i) { 
    cout << U[(0)*3+i] << "   " << U[(1)*3+i] << "    " << U[(2)*3+i] << endl;
  }

  cout << endl;
  cout << endl;
  cout << " ==== V === " << endl;
  for (int i =0; i < 3; ++i) { 
    cout << V[(0)*3+i] << "   " << V[(1)*3+i] << "    " << V[(2)*3+i] << endl;
  }

  double A_new[9]; for (int ii =0; ii < 9; ++ii) A_new[ii] = 0.0;

  for (int i =0; i < 3; ++i) { 
    if ( sigma[i] < 0.1) { 
      sigma[i] = 1.0;
    }
  }

  for (int j = 0; j < 3; ++j) { 
    for (int i =0; i < 3; ++i) { 
      
      double aij = 0.0;
      for (int k =0; k < 3; ++k) 
        aij += U[(k)*3+i]*V[(k)*3+j]*sigma[k];
      
      A_new[(3)*j+i] = aij;
    }
  }

  cout << " ==== A_new === " << endl;
  for (int i =0; i < 3; ++i) { 
    cout << A_new[(0)*3+i] << "   " << A_new[(1)*3+i] << "    " << A_new[(2)*3+i] << endl;
  }

  return 0;
} 

