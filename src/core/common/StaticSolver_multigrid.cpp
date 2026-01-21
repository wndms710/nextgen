
#include "StaticSolver.hpp"

void StaticSolver::smoothCvJacobi(double* phi,double *tmp,const double* inv_diag,const double *A,const double *rhs,const int nsmooth,const double relax) {

  for (int ii = 0; ii < nsmooth; ++ii) {

    calcCvResidual(tmp,phi,A,rhs);

    for (int icv = 0; icv < ncv; ++icv) 
      phi[icv] += relax*tmp[icv]*inv_diag[icv];
    updateCvData(phi);

  }

}

void StaticSolver::smoothCvJacobi(double (*u)[3],double (*tmp)[3],const double* inv_diag,const double *A,const double (*rhs)[3],const int nsmooth,const double relax) {

  for (int ii = 0; ii < nsmooth; ++ii) {

    calcCvResidual(tmp,u,A,rhs);

    for (int icv = 0; icv < ncv; ++icv) 
      FOR_I3 u[icv][i] += relax*tmp[icv][i]*inv_diag[icv];
    updateCvData(u,REPLACE_ROTATE_DATA);

  }

}

void StaticSolver::smoothCvJacobi(double (*u)[3],double (*tmp)[3],const double (*inv_diag)[3],const double *A,const double (*A_diag)[3],const double (*rhs)[3],const int nsmooth,const double relax) {

  for (int ii = 0; ii < nsmooth; ++ii) {

    calcCvResidual(tmp,u,A,A_diag,rhs);

    for (int icv = 0; icv < ncv; ++icv) 
      FOR_I3 u[icv][i] += relax*tmp[icv][i]*inv_diag[icv][i];
    updateCvData(u,REPLACE_ROTATE_DATA);

  }

}

void StaticSolver::smoothCvGs(double* phi,const double* inv_diag,const double *A,const double *rhs,const int nsmooth,const double relax) {

  // note that this is only Gauss-Seidel locally

  for (int ii = 0; ii < nsmooth; ++ii) {

    for (int icv = 0; icv < ncv; ++icv) {
      const double phi0 = phi[icv];
      phi[icv] = rhs[icv];
      for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        phi[icv] -= A[coc]*phi[icv_nbr];
      }
      phi[icv] = phi0 + relax*(phi[icv]*inv_diag[icv]-phi0);
    }
    updateCvData(phi);

  }

}

void StaticSolver::smoothCvSgs(double* phi,const double* inv_diag,const double *A,const double *rhs,const int nsmooth,const double relax) {

  // note that this is only symmetric Gauss-Seidel locally

  for (int ii = 0; ii < nsmooth; ++ii) {

    for (int icv = 0; icv < ncv; ++icv) {
      const double phi0 = phi[icv];
      phi[icv] = rhs[icv];
      for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        phi[icv] -= A[coc]*phi[icv_nbr];
      }
      phi[icv] = phi0 + relax*(phi[icv]*inv_diag[icv]-phi0);
    }
    updateCvData(phi);
    for (int icv = ncv-1; icv >= 0; --icv) {
      const double phi0 = phi[icv];
      phi[icv] = rhs[icv];
      for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        phi[icv] -= A[coc]*phi[icv_nbr];
      }
      phi[icv] = phi0 + relax*(phi[icv]*inv_diag[icv]-phi0);
    }
    updateCvData(phi);

  }

}

void StaticSolver::smoothCvSgs(double (*u)[3],const double* inv_diag,const double *A,const double (*rhs)[3],const int nsmooth,const double relax) {

  // note that this is only symmetric Gauss-Seidel locally

  for (int ii = 0; ii < nsmooth; ++ii) {

    for (int icv = 0; icv < ncv; ++icv) {
      const double u0[3] = {u[icv][0],u[icv][1],u[icv][2]};
      FOR_I3 u[icv][i] = rhs[icv][i];
      for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        FOR_I3 u[icv][i] -= A[coc]*u[icv_nbr][i];
      }
      FOR_I3 u[icv][i] = u0[i] + relax*(u[icv][i]*inv_diag[icv]-u0[i]);
    }
    updateCvData(u,REPLACE_ROTATE_DATA);
    for (int icv = ncv-1; icv >= 0; --icv) {
      const double u0[3] = {u[icv][0],u[icv][1],u[icv][2]};
      FOR_I3 u[icv][i] = rhs[icv][i];
      for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        FOR_I3 u[icv][i] -= A[coc]*u[icv_nbr][i];
      }
      FOR_I3 u[icv][i] = u0[i] + relax*(u[icv][i]*inv_diag[icv]-u0[i]);
    }
    updateCvData(u,REPLACE_ROTATE_DATA);

  }

}

void StaticSolver::smoothCvGs(double (*u)[3],const double* inv_diag,const double *A,const double (*rhs)[3],const int nsmooth,const double relax) {

  for (int ii = 0; ii < nsmooth; ++ii) {

    for (int icv = 0; icv < ncv; ++icv) {
      const double u0[3] = {u[icv][0],u[icv][1],u[icv][2]};
      FOR_I3 u[icv][i] = rhs[icv][i];
      for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        FOR_I3 u[icv][i] -= A[coc]*u[icv_nbr][i];
      }
      FOR_I3 u[icv][i] = u0[i] + relax*(u[icv][i]*inv_diag[icv]-u0[i]);
    }
    updateCvData(u,REPLACE_ROTATE_DATA);

  }

}

void StaticSolver::smoothCvPatr(double* phi,
    double* r,double* p,double* Ar,double* Ap, // work arrays
    const double* inv_diag,const double* A,const double* At,const double* rhs,
    const int nsmooth,const double relax1,const double relax2) {

  // assume we come in with a consistent initial condition...

  // residial: r = b-A*x...
  FOR_ICV {
    r[icv] = A[cvocv_i[icv]]*phi[icv] - rhs[icv];
    for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      r[icv] += A[coc]*phi[icv_nbr];
    }
    r[icv] *= inv_diag[icv];
  }

  for (int ii = 0; ii < nsmooth; ++ii) {

    // need residual in ghosts...
    updateCvData(r); 

    // search direction: p = [A^T]{r}...
    FOR_ICV {
      p[icv] = At[cvocv_i[icv]]*inv_diag[icv]*r[icv];
      for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        p[icv] += At[coc]*inv_diag[icv_nbr]*r[icv_nbr];
      }
    }
    updateCvData(p); 

    // second search direction is r...

    // get work arrays: Ap and Ar...
    FOR_ICV {
      Ar[icv] = A[cvocv_i[icv]]*r[icv];
      Ap[icv] = A[cvocv_i[icv]]*p[icv];
      for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        Ar[icv] += A[coc]*r[icv_nbr];
        Ap[icv] += A[coc]*p[icv_nbr];
      }
      Ar[icv] *= inv_diag[icv];
      Ap[icv] *= inv_diag[icv];
    }

    // calc alpha and beta...
    double my_buf[5] = {0.0,0.0,0.0,0.0,0.0}; 
    FOR_ICV {
      my_buf[0] += Ap[icv]*Ap[icv];
      my_buf[1] += Ap[icv]*Ar[icv];
      my_buf[2] += Ar[icv]*Ar[icv];
      my_buf[3] -=  r[icv]*Ap[icv];
      my_buf[4] -=  r[icv]*Ar[icv];
    }
    double buf[5];
    MPI_Allreduce(my_buf,buf,5,MPI_DOUBLE,MPI_SUM,mpi_comm);
    double alpha,beta; 
    const double inv_det = 1.0/(buf[0]*buf[2] - buf[1]*buf[1] + 1.0E-40);
    assert(inv_det > 0.0);
    alpha = relax1*inv_det*( buf[2]*buf[3] - buf[1]*buf[4]);
    beta  = relax2*inv_det*(-buf[1]*buf[3] + buf[0]*buf[4]);

    // phi += alpha*p + beta*r...
    // r += alpha*Ap + beta*Ar...
    FOR_ICV_G phi[icv] += alpha*p[icv]  + beta*r[icv];
    FOR_ICV r[icv] += alpha*Ap[icv] + beta*Ar[icv];


    /*
       {
       double *res = new double[ncv];
       calcCvResidual(res,phi,A,rhs);
       double my_buf[2] = {0.0,0.0};
       FOR_ICV {
       my_buf[0] += vol_cv[icv]*fabs(res[icv]);
       my_buf[1] += vol_cv[icv];
       }
       double buf[2];
       MPI_Reduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
       if (mpi_rank == 0)
       cout << " res_avg: " << buf[0]/buf[1] << endl;;
       }
       */
  } 

}

void StaticSolver::smoothCvPatr(double (*u)[3],
    double (*r)[3],double (*p)[3],double (*Ar)[3],double (*Ap)[3], // work arrays
    const double* inv_diag,const double* A,const double* At,const double (*rhs)[3],
    const int nsmooth,const double relax1,const double relax2) {

  // assume we come in with a consistent initial condition...

  // residial: r = b-A*x...
  FOR_ICV {
    FOR_I3 r[icv][i] = A[cvocv_i[icv]]*u[icv][i] - rhs[icv][i];
    for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      FOR_I3 r[icv][i] += A[coc]*u[icv_nbr][i];
    }
    FOR_I3 r[icv][i] *= inv_diag[icv];
  }

  for (int ii = 0; ii < nsmooth; ++ii) {

    // need residual in ghosts...
    updateCvData(r,REPLACE_ROTATE_DATA); 

    // search direction: p = [A^T]{r}...
    FOR_ICV {
      FOR_I3 p[icv][i] = At[cvocv_i[icv]]*inv_diag[icv]*r[icv][i];
      for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        FOR_I3 p[icv][i] += At[coc]*inv_diag[icv_nbr]*r[icv_nbr][i];
      }
    }
    updateCvData(p,REPLACE_ROTATE_DATA); 

    // second search direction is r...

    // get work arrays: Ap and Ar...
    FOR_ICV {
      FOR_I3 {
        Ar[icv][i] = A[cvocv_i[icv]]*r[icv][i];
        Ap[icv][i] = A[cvocv_i[icv]]*p[icv][i];
      }
      for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        FOR_I3 {
          Ar[icv][i] += A[coc]*r[icv_nbr][i];
          Ap[icv][i] += A[coc]*p[icv_nbr][i];
        }
      }
      FOR_I3 {
        Ar[icv][i] *= inv_diag[icv];
        Ap[icv][i] *= inv_diag[icv];
      }
    }

    // calc alpha and beta...
    //                      0   1   2
    double my_buf[5][3] = {{0.0,0.0,0.0},  // Ap^T*Ap 
      {0.0,0.0,0.0},  // Ap^T*Ar 
      {0.0,0.0,0.0},  // Ar^T*Ar 
      {0.0,0.0,0.0},  //  r^T*Ap 
      {0.0,0.0,0.0}}; //  r^T*Ar 
    FOR_ICV {
      FOR_I3 {
        my_buf[0][i] += Ap[icv][i]*Ap[icv][i];
        my_buf[1][i] += Ap[icv][i]*Ar[icv][i];
        my_buf[2][i] += Ar[icv][i]*Ar[icv][i];
        my_buf[3][i] -=  r[icv][i]*Ap[icv][i];
        my_buf[4][i] -=  r[icv][i]*Ar[icv][i];
      }
    }
    double buf[5][3];
    MPI_Allreduce((double*)my_buf,(double*)buf,15,MPI_DOUBLE,MPI_SUM,mpi_comm);
    double alpha[3],beta[3]; 
    FOR_I3 {
      const double inv_det = 1.0/(buf[0][i]*buf[2][i] - buf[1][i]*buf[1][i] + 1.0E-40);
      assert(inv_det > 0.0);
      alpha[i] = relax1*inv_det*( buf[2][i]*buf[3][i] - buf[1][i]*buf[4][i]);
      beta[i]  = relax2*inv_det*(-buf[1][i]*buf[3][i] + buf[0][i]*buf[4][i]);
    }

    // u += alpha*p + beta*r...
    // r += alpha*Ap + beta*Ar...
    FOR_ICV_G FOR_I3 u[icv][i] += alpha[i]*p[icv][i]  + beta[i]*r[icv][i];
    FOR_ICV FOR_I3 r[icv][i] += alpha[i]*Ap[icv][i] + beta[i]*Ar[icv][i];

  }

  /*
     {
     double (*res)[3] = new double[ncv][3];
     calcCvResidual(res,u,A,rhs);
     double my_buf[2] = {0.0,0.0};
     FOR_ICV {
     FOR_I3 my_buf[0] += vol_cv[icv]*fabs(res[icv][i]);
     my_buf[1] += vol_cv[icv];
     }
     double buf[2];
     MPI_Reduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
     if (mpi_rank == 0)
     cout << " res_avg: " << buf[0]/buf[1] << endl;;
     delete[] res;
     }
     */

}

void StaticSolver::smoothCvPatr(double (*u)[3],
    double (*r)[3],double (*p)[3],double (*Ar)[3],double (*Ap)[3], // work arrays
    const double (*inv_diag)[3],const double* A,const double* At,const double (*A_diag)[3], const double (*rhs)[3],
    const int nsmooth,const double relax1,const double relax2) {

  // assume we come in with a consistent initial condition...

  // residial: r = b-A*x...
  FOR_ICV {
    FOR_I3 r[icv][i] = (A[cvocv_i[icv]]+A_diag[icv][i])*u[icv][i] - rhs[icv][i];
    for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      FOR_I3 r[icv][i] += A[coc]*u[icv_nbr][i];
    }
    FOR_I3 r[icv][i] *= inv_diag[icv][i];
  }

  for (int ii = 0; ii < nsmooth; ++ii) {

    // need residual in ghosts...
    updateCvData(r,REPLACE_ROTATE_DATA); 

    // search direction: p = [A^T]{r}...
    FOR_ICV {
      FOR_I3 p[icv][i] = (At[cvocv_i[icv]]+A_diag[icv][i])*inv_diag[icv][i]*r[icv][i];
      for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        FOR_I3 p[icv][i] += At[coc]*inv_diag[icv_nbr][i]*r[icv_nbr][i];
      }
    }
    updateCvData(p,REPLACE_ROTATE_DATA); 

    // second search direction is r...

    // get work arrays: Ap and Ar...
    FOR_ICV {
      FOR_I3 {
        Ar[icv][i] = (A[cvocv_i[icv]]+A_diag[icv][i])*r[icv][i];
        Ap[icv][i] = (A[cvocv_i[icv]]+A_diag[icv][i])*p[icv][i];
      }
      for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        FOR_I3 {
          Ar[icv][i] += A[coc]*r[icv_nbr][i];
          Ap[icv][i] += A[coc]*p[icv_nbr][i];
        }
      }
      FOR_I3 {
        Ar[icv][i] *= inv_diag[icv][i];
        Ap[icv][i] *= inv_diag[icv][i];
      }
    }

    // calc alpha and beta...
    //                      0   1   2
    double my_buf[5][3] = {{0.0,0.0,0.0},  // Ap^T*Ap 
      {0.0,0.0,0.0},  // Ap^T*Ar 
      {0.0,0.0,0.0},  // Ar^T*Ar 
      {0.0,0.0,0.0},  //  r^T*Ap 
      {0.0,0.0,0.0}}; //  r^T*Ar 
    FOR_ICV {
      FOR_I3 {
        my_buf[0][i] += Ap[icv][i]*Ap[icv][i];
        my_buf[1][i] += Ap[icv][i]*Ar[icv][i];
        my_buf[2][i] += Ar[icv][i]*Ar[icv][i];
        my_buf[3][i] -=  r[icv][i]*Ap[icv][i];
        my_buf[4][i] -=  r[icv][i]*Ar[icv][i];
      }
    }
    double buf[5][3];
    MPI_Allreduce((double*)my_buf,(double*)buf,15,MPI_DOUBLE,MPI_SUM,mpi_comm);
    double alpha[3],beta[3]; 
    FOR_I3 {
      const double inv_det = 1.0/(buf[0][i]*buf[2][i] - buf[1][i]*buf[1][i] + 1.0E-40);
      assert(inv_det > 0.0);
      alpha[i] = relax1*inv_det*( buf[2][i]*buf[3][i] - buf[1][i]*buf[4][i]);
      beta[i]  = relax2*inv_det*(-buf[1][i]*buf[3][i] + buf[0][i]*buf[4][i]);
    }

    // u += alpha*p + beta*r...
    // r += alpha*Ap + beta*Ar...
    FOR_ICV_G FOR_I3 u[icv][i] += alpha[i]*p[icv][i]  + beta[i]*r[icv][i];
    FOR_ICV FOR_I3 r[icv][i] += alpha[i]*Ap[icv][i] + beta[i]*Ar[icv][i];

  }

}



void StaticSolver::smoothCvTim(double (*&u)[3],
    double (*&w)[3], // work array
    const double *inv_diag,const double* A,const double* invDAs,const double (*rhs)[3],
    const int nsmooth,const double tau_bound,const double relax) {

  // NOTE: invDAs is the skew part of D^(-1)*A

  const double tau = relax*tau_bound;
  for (int ii = 0; ii < nsmooth; ++ii) {

    FOR_ICV {
      FOR_I3 w[icv][i] = u[icv][i] + tau*inv_diag[icv]*(rhs[icv][i] - A[cvocv_i[icv]]*u[icv][i]);
      for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        if (icv_nbr > icv) {
          FOR_I3 w[icv][i] -= tau*inv_diag[icv]*A[coc]*u[icv_nbr][i];
        }
        else {
          assert(icv_nbr < icv);
          FOR_I3 w[icv][i] -= tau*((inv_diag[icv]*A[coc]-2*invDAs[coc])*u[icv_nbr][i] + 2.0*invDAs[coc]*w[icv_nbr][i]);
        }
      }
    }
    updateCvData(w,REPLACE_ROTATE_DATA);

    // just swap u and w
    double (*tmp)[3] = u;
    u = w;
    w = tmp;

  }

}

void StaticSolver::smoothCvCg(double *phi,
    double *res,double* v,double *p, // work arrays
    const double *inv_diag,const double *A,const double *rhs,const int nsmooth) {

  // assume we come in with a consistent initial condition...

  for (int icv = 0; icv < ncv; ++icv)
    p[icv] = 0.0;
  double rho = 1.0;

  // calculate the residual in rhs format...
  calcCvResidual(res,phi,A,rhs);

  // diagonal precon/compute normalized residual...
  for (int icv = 0; icv < ncv; ++icv)
    v[icv] = res[icv]*inv_diag[icv];

  for (int ii = 0; ii < nsmooth; ++ii) {

    double rho_prev = rho;
    if (fabs(rho_prev) < 1.0E-20)
      rho_prev = -1.0E-20; // -1.0E-20? seems to help

    double my_rho = 0.0;
    for (int icv = 0; icv < ncv; ++icv)
      my_rho += res[icv]*v[icv];
    MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    double beta = rho/rho_prev;
    for (int icv = 0; icv < ncv; ++icv)
      p[icv] = v[icv] + beta*p[icv];
    updateCvData(p);

    // v = [Ap]{p}...
    for (int icv = 0; icv < ncv; ++icv) {
      v[icv] = A[cvocv_i[icv]]*p[icv];
      for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        v[icv] += A[coc]*p[icv_nbr];
      }
    }

    double my_gamma = 0.0;
    for (int icv = 0; icv < ncv; ++icv)
      my_gamma += p[icv]*v[icv];
    double gamma;
    MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    if (fabs(gamma) < 1.0E-20)
      gamma = 1.0E-20;

    const double alpha = rho/gamma;

    // update full phi including ghosts...
    for (int icv = 0; icv < ncv_g; ++icv)
      phi[icv] += alpha*p[icv];

    for (int icv = 0; icv < ncv; ++icv) {
      // on the other iterations, use this approximation to update
      // the unreduced residual...
      res[icv] -= alpha*v[icv];
      // still need to compute v, diag precon for next iteration...
      v[icv] = res[icv]*inv_diag[icv];
    }

  }

}

void StaticSolver::smoothCvCg(double (*u)[3],
    double (*res)[3],double (*v)[3],double (*p)[3], // work arrays
    const double *inv_diag,const double *A,const double (*rhs)[3],const int nsmooth) {

  // assume we come in with a consistent initial condition...

  for (int icv = 0; icv < ncv; ++icv)
    FOR_I3 p[icv][i] = 0.0;
  double rho[3] = {1.0,1.0,1.0};

  // calculate the residual in rhs format...
  calcCvResidual(res,u,A,rhs);

  // diagonal precon/compute normalized residual...
  for (int icv = 0; icv < ncv; ++icv)
    FOR_I3 v[icv][i] = res[icv][i]*inv_diag[icv];

  for (int ii = 0; ii < nsmooth; ++ii) {

    double rho_prev[3] = {rho[0],rho[1],rho[2]};
    FOR_I3 {
      if (fabs(rho_prev[i]) < 1.0E-20)
        rho_prev[i] = -1.0E-20; // -1.0E-20? seems to help
    }

    double my_rho[3] = {0.0,0.0,0.0};
    for (int icv = 0; icv < ncv; ++icv)
      FOR_I3 my_rho[i] += res[icv][i]*v[icv][i];
    MPI_Allreduce(my_rho,rho,3,MPI_DOUBLE,MPI_SUM,mpi_comm);

    double beta[3]; FOR_I3 beta[i] = rho[i]/rho_prev[i];
    for (int icv = 0; icv < ncv; ++icv)
      FOR_I3 p[icv][i] = v[icv][i] + beta[i]*p[icv][i];
    updateCvData(p);

    // v = [Ap]{p}...
    for (int icv = 0; icv < ncv; ++icv) {
      FOR_I3 v[icv][i] = A[cvocv_i[icv]]*p[icv][i];
      for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        FOR_I3 v[icv][i] += A[coc]*p[icv_nbr][i];
      }
    }

    double my_gamma[3] = {0.0,0.0,0.0};
    for (int icv = 0; icv < ncv; ++icv)
      FOR_I3 my_gamma[i] += p[icv][i]*v[icv][i];
    double gamma[3];
    MPI_Allreduce(my_gamma,gamma,3,MPI_DOUBLE,MPI_SUM,mpi_comm);
    FOR_I3 {
      if (fabs(gamma[i]) < 1.0E-20)
        gamma[i] = 1.0E-20;
    }

    double alpha[3]; FOR_I3 alpha[i] = rho[i]/gamma[i];

    // update full u including ghosts...
    for (int icv = 0; icv < ncv_g; ++icv)
      FOR_I3 u[icv][i] += alpha[i]*p[icv][i];

    for (int icv = 0; icv < ncv; ++icv) {
      // on the other iterations, use this approximation to update
      // the unreduced residual...
      FOR_I3 res[icv][i] -= alpha[i]*v[icv][i];
      // still need to compute v, diag precon for next iteration...
      FOR_I3 v[icv][i] = res[icv][i]*inv_diag[icv];
    }

  }

}

void StaticSolver::CoarseGrid::initCoarseGrid(const int level,const double agglomeration_factor,const bool split_orphaned_colors) {

  // color base grid based to get new global indices

  ncc_global = int8(double(solver->ncv_global)/pow(agglomeration_factor,level)); // 2^3 = 8.0

  if (mpi_rank == 0)
    cout << "CoarseGrid::initCoarseGrid(), level: " << level << ", requested ncc_global: " << ncc_global << endl;

  int8* icc_global_cv = new int8[solver->ncv_g];

  Multigrid::colorCvsPadt(icc_global_cv,ncc_global,solver->x_cv,solver->icv_global,solver->ncv,solver->ncv_global);
  //FILE * fp = fopen("color_cv2.dat","w");
  //for (int icv = 0; icv < solver->ncv; ++icv) 
  //  fprintf(fp,"%lld %lld %lf %lf %lf\n",icc_global_cv[icv],solver->icv_global[icv],solver->x_cv[icv][0],solver->x_cv[icv][1],solver->x_cv[icv][2]);
  //fclose(fp);
  solver->updateCvData(icc_global_cv);
  if (split_orphaned_colors) {
    Multigrid::splitOrphanedColors(icc_global_cv,ncc_global,solver->cvocv_i,solver->cvocv_v,solver->icv_global,
        solver->rbi_g,NULL,solver->ncv,solver->ncv_g,solver->ncv_g,solver->ncv_global); // NULL because we are using compact grid only
  }

  if (mpi_rank == 0)
    cout << " final ncc_global: " << ncc_global << endl;

  map<const int8,int> globalMap;
  assert(ncc == 0);
  for (int icv = 0; icv < solver->ncv; ++icv) {
    map<const int8,int>::iterator iter = globalMap.find(icc_global_cv[icv]);
    if (iter == globalMap.end()) 
      globalMap[icc_global_cv[icv]] = ncc++;
  }
  assert(globalMap.size() == ncc);

  // figure out which coarse cells (colors) are split amongst ranks

  assert(icc_global == NULL); icc_global = new int8[ncc]; 
  for (map<const int8,int>::iterator iter = globalMap.begin(); iter != globalMap.end(); ++iter) 
    icc_global[iter->second] = iter->first;

  // send to striping so we can rectify shared/ghost icc's

  int8 *ccora_striped = NULL;
  MiscUtils::calcUniformDist(ccora_striped,ncc_global,mpi_size);

  int * send_count = new int[mpi_size];
  FOR_RANK send_count[rank] = 0;
  for (int icc = 0; icc < ncc; ++icc) {
    const int rank = MiscUtils::getRankInXora(icc_global[icc],ccora_striped);
    send_count[rank] += 3; // mpi_rank,icc,icc_striped
  }

  int * send_disp = new int[mpi_size];
  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
  int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

  int* send_buf_int = new int[send_count_sum]; 
  for (int icc = 0; icc < ncc; ++icc) {
    const int rank = MiscUtils::getRankInXora(icc_global[icc],ccora_striped);
    send_buf_int[send_disp[rank]++] = mpi_rank;
    send_buf_int[send_disp[rank]++] = icc;
    send_buf_int[send_disp[rank]++] = icc_global[icc]-ccora_striped[rank];
  }

  // rewind...

  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

  // setup recv side stuff

  int * recv_count = new int[mpi_size];
  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

  int * recv_disp = new int[mpi_size];
  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
  int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

  int * recv_buf_int = new int[recv_count_sum];
  MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
      recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
  delete[] send_buf_int; send_buf_int = NULL;

  const int nccs = ccora_striped[mpi_rank+1]-ccora_striped[mpi_rank];
  int* ccoccs_i = new int[nccs+1];
  for (int iccs = 0; iccs < nccs; ++iccs) 
    ccoccs_i[iccs+1] = 0;
  int8* ccoccs_v_global = NULL;
  for (int iter = 0; iter < 2; ++iter) {
    for (int irecv = 0; irecv < recv_count_sum; irecv += 3) {
      const int iccs = recv_buf_int[irecv+2];
      assert((iccs >= 0)&&(iccs < nccs));
      if (iter == 0) 
        ++ccoccs_i[iccs+1];
      else 
        ccoccs_v_global[ccoccs_i[iccs]++] = BitUtils::packRankBitsIndex(recv_buf_int[irecv+0],0,recv_buf_int[irecv+1]);
    }
    if (iter == 0) {
      ccoccs_i[0] = 0;
      for (int iccs = 0; iccs < nccs; ++iccs) ccoccs_i[iccs+1] += ccoccs_i[iccs];
      ccoccs_v_global = new int8[ccoccs_i[nccs]];
    }
    else {
      // rewind
      for (int iccs = nccs; iccs > 0; --iccs)
        ccoccs_i[iccs] = ccoccs_i[iccs-1];
      ccoccs_i[0] = 0;
    }
  }
  delete[] recv_buf_int; recv_buf_int = NULL;

  FOR_RANK send_count[rank] = 0;
  for (int iter = 0; iter < 2; ++iter) {
    for (int iccs = 0; iccs < nccs; ++iccs) {
      const int nbro = ccoccs_i[iccs+1]-ccoccs_i[iccs]-1; // bros have same icc_global (exclude yourself)
      for (int cocs = ccoccs_i[iccs]; cocs < ccoccs_i[iccs+1]; ++cocs) {
        int rank,bits,index; 
        BitUtils::unpackRankBitsIndex(rank,bits,index,ccoccs_v_global[cocs]);
        assert(bits == 0);
        if (iter == 0) {
          send_count[rank] += 2+2*nbro; // icc,nbro,(rank_bro,icc_bro)
        }
        else {
          send_buf_int[send_disp[rank]++] = index; 
          send_buf_int[send_disp[rank]++] = nbro;
          for (int cocs_ = ccoccs_i[iccs]; cocs_ < ccoccs_i[iccs+1]; ++cocs_) {
            if (cocs_ != cocs) {
              int rank_,bits_,index_; 
              BitUtils::unpackRankBitsIndex(rank_,bits_,index_,ccoccs_v_global[cocs_]);
              assert(bits_ == 0);
              send_buf_int[send_disp[rank]++] = rank_; 
              send_buf_int[send_disp[rank]++] = index_;
            }
          }
        }
      }
    }

    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

    if (iter == 0) {
      send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
      assert(send_buf_int == NULL); 
      send_buf_int = new int[send_count_sum];
    }
  }
  delete[] ccoccs_i;
  delete[] ccoccs_v_global;

  // setup recv side stuff

  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

  recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
  assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum];
  MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
      recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
  delete[] send_buf_int; send_buf_int = NULL;

  int* brocc_i = new int[ncc+1];
  for (int icc = 0; icc < ncc; ++icc) 
    brocc_i[icc+1] = 0;
  int8* brocc_v_global = NULL;
  for (int iter = 0; iter < 2; ++iter) {
    FOR_RANK {
      int irecv = recv_disp[rank];
      while (irecv < recv_disp[rank]+recv_count[rank]) {
        const int icc = recv_buf_int[irecv++];
        const int nbro = recv_buf_int[irecv++];
        if (iter == 0) {
          brocc_i[icc+1] += nbro;
          irecv += 2*nbro;
        }
        else {
          for (int ibro = 0; ibro < nbro; ++ibro) {
            const int rank = recv_buf_int[irecv++];
            const int index = recv_buf_int[irecv++];
            brocc_v_global[brocc_i[icc]++] = BitUtils::packRankBitsIndex(rank,0,index);
          }
        }
      }
    }
    if (iter == 0) {
      brocc_i[0] = 0;
      for (int icc = 0; icc < ncc; ++icc) brocc_i[icc+1] += brocc_i[icc];
      brocc_v_global = new int8[brocc_i[ncc]];
    }
    else {
      // rewind
      for (int icc = ncc; icc > 0; --icc)
        brocc_i[icc] = brocc_i[icc-1];
      brocc_i[0] = 0;
    }
  }
  delete[] recv_buf_int; recv_buf_int = NULL;

  // now figure out which cells are split/active (min rank wins)

  vector<pair<int8,int> > rbiVec(ncc);
  assert(ncc_a == 0);
  assert(ncc_in == 0);
  for (int icc = 0; icc < ncc; ++icc) {
    rbiVec[icc].second = icc;
    int min_rank = mpi_rank;
    if (brocc_i[icc+1]-brocc_i[icc] == 0) {
      rbiVec[icc].first = -2; // just set to -2 if internal
      ++ncc_in; // not shared so this is an internal coarse cell
    }
    else {
      rbiVec[icc].first = -1; // just set to -1 if active on this rank
      for (int boc = brocc_i[icc]; boc < brocc_i[icc+1]; ++boc) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,brocc_v_global[boc]);
        assert(rank != mpi_rank);
        if (rank < min_rank) {
          rbiVec[icc].first = brocc_v_global[boc];
          min_rank = rank;
        }
      }
      if (min_rank == mpi_rank) {
        assert(rbiVec[icc].first == -1);
        ++ncc_a; // min rank wins, so this is the active coarse cell
      }
      else {
        assert(rbiVec[icc].first >= 0);
      }
    }
  }
  ncc_a += ncc_in;
  delete[] brocc_i;
  delete[] brocc_v_global; // do we need to keep/reorder these for active communicator

  // now we need to reorder coarse cells like internal,active,inactive (rbi ordered)
  sort(rbiVec.begin(),rbiVec.end());

  int *reorder_cc = new int[ncc];
  for (int ii = 0; ii < ncc; ++ii)
    reorder_cc[rbiVec[ii].second] = ii;
  const int ncc_check = ncc;
  const int ncc_in_check = ncc_in;
  const int ncc_a_check = ncc_a;
  ncc = ncc_a;
  ncc_a = ncc_in;
  ncc_in = 0;
  for (int icc = 0; icc < ncc_check; ++icc) {
    if (rbiVec[icc].first == -2)
      ++ncc_in;
    else if (rbiVec[icc].first == -1)
      ++ncc_a;
    else 
      ++ncc;
  }
  assert(ncc == ncc_check);
  assert(ncc_in == ncc_in_check);
  assert(ncc_a == ncc_a_check);

  int8 my_ncc = (int8)ncc_a;
  int8 ncc_global_check;
  MPI_Reduce(&my_ncc,&ncc_global_check,1,MPI_INT8,MPI_SUM,0,mpi_comm);
  if (mpi_rank == 0) 
    assert(ncc_global_check == ncc_global);

  assert(rbi_i == NULL); rbi_i = new int8[ncc-ncc_a];
  for (map<const int8,int>::iterator iter = globalMap.begin(); iter != globalMap.end(); ++iter) {
    const int icc_old = iter->second; 
    const int icc = reorder_cc[icc_old];
    icc_global[icc] = iter->first;
    iter->second = icc;
    if (icc >= ncc_a) { 
      assert(rbiVec[icc].first >= 0);
      rbi_i[icc-ncc_a] = rbiVec[icc].first;
    }
    else if (icc >= ncc_in) {
      assert(rbiVec[icc].first == -1);
    }
    else {
      assert(rbiVec[icc].first == -2);
    }
  }
  rbiVec.clear();

  // we need to update rbi_inactive to respect the reorder on nbr rank

  FOR_RANK send_count[rank] = 0;
  int * inactive_index = new int[ncc-ncc_a];
  for (int iter = 0; iter < 2; ++iter) {
    for (int icc = ncc_a; icc < ncc; ++icc) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_i[icc-ncc_a]);
      if (iter == 0) {
        ++send_count[rank];
      }
      else {
        inactive_index[send_disp[rank]] = icc-ncc_a;
        send_buf_int[send_disp[rank]++] = index;
      }
    }

    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

    if (iter == 0) {
      send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
      assert(send_buf_int == NULL); 
      send_buf_int = new int[send_count_sum];
    }

  }

  // setup recv side stuff

  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

  recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
  assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum];
  MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
      recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);

  for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
    const int icc_old = recv_buf_int[irecv]; assert((icc_old >= 0)&&(icc_old < ncc)); 
    const int icc = reorder_cc[icc_old]; assert((icc >= ncc_in)&&(icc < ncc_a)); // should be active 
    recv_buf_int[irecv] = icc;
  }
  delete[] reorder_cc; reorder_cc = NULL;

  // send back

  MPI_Alltoallv(recv_buf_int,recv_count,recv_disp,MPI_INT,
      send_buf_int,send_count,send_disp,MPI_INT,mpi_comm);
  delete[] recv_buf_int; recv_buf_int = NULL;

  // overwrite rbi_i

  FOR_RANK {
    int isend = send_disp[rank];
    while (isend < send_disp[rank]+send_count[rank]) {
      const int8 rbi = BitUtils::packRankBitsIndex(rank,0,send_buf_int[isend]);
      rbi_i[inactive_index[isend]] = rbi;
      ++isend;
    }
  }
  delete[] send_buf_int; send_buf_int = NULL;
  delete[] inactive_index;

  assert(cvocc_i == NULL); cvocc_i = new int[ncc+1];
  for (int icc = 0; icc < ncc; ++icc) cvocc_i[icc+1] = 0;
  for (int iter = 0; iter < 2; ++iter) {
    for (int icv = 0; icv < solver->ncv; ++icv) {
      map<const int8,int>::iterator iter2 = globalMap.find(icc_global_cv[icv]);
      assert(iter2 != globalMap.end());
      const int icc = iter2->second;
      if (iter == 0) 
        ++cvocc_i[icc+1];
      else 
        cvocc_v[cvocc_i[icc]++] = icv;
    }
    if (iter == 0) {
      cvocc_i[0] = 0;
      for (int icc = 0; icc < ncc; ++icc) cvocc_i[icc+1] += cvocc_i[icc];
      assert(solver->ncv == cvocc_i[ncc]);
      assert(cvocc_v == NULL); cvocc_v = new int[cvocc_i[ncc]];
    }
    else {
      // rewind
      for (int icc = ncc; icc > 0; --icc)
        cvocc_i[icc] = cvocc_i[icc-1];
      cvocc_i[0] = 0;
    }
  }

  dumpRange(&ncc_in_check,1,"ncc internal");
  dumpRange(&ncc_a,1,"ncc active");
  dumpRange(&ncc,1,"ncc");

  // now lets get the ghosts

  assert(ccocv == NULL); ccocv = new int[solver->ncv_g]; // fill ghosts once we have them
  for (int icc = 0; icc < ncc; ++icc) {
    for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
      const int icv = cvocc_v[coc];
      ccocv[icv] = icc;
    }
  }
  ncc_g = ncc;
  assert(globalMap.size() == ncc); 
  for (int ifa = solver->nfa_i; ifa < solver->nfa; ++ifa) {
    const int icv0 = solver->cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < solver->ncv));
    const int icv1 = solver->cvofa[ifa][1]; assert((icv1 >= solver->ncv)&&(icv1 < solver->ncv_g));
    assert(icc_global[ccocv[icv0]] == icc_global_cv[icv0]);
    if (icc_global_cv[icv0] != icc_global_cv[icv1]) {
      map<const int8,int>::iterator iter = globalMap.find(icc_global_cv[icv1]);
      if (iter == globalMap.end()) {
        ccocv[icv1] = ncc_g;
        globalMap[icc_global_cv[icv1]] = ncc_g++; // add in the traditional ghosts
      }
      else {
        ccocv[icv1] = iter->second;
      }
    }
    else {
      ccocv[icv1] = ccocv[icv0];
    }
  }
  assert(globalMap.size() == ncc_g);

  // add inactive rbi's 
  assert(rbi_g == NULL); rbi_g = new int8[ncc_g-ncc];
  for (int icc = ncc; icc < ncc_g; ++icc) 
    rbi_g[icc-ncc] = -1;

  int8* rbicc_cv = new int8[solver->ncv_g];
  for (int icc = 0; icc < ncc_a; ++icc) {
    for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
      const int icv = cvocc_v[coc];
      rbicc_cv[icv] = BitUtils::packRankBitsIndex(mpi_rank,0,icc);
    }
  }
  for (int icc = ncc_a; icc < ncc; ++icc) {
    for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
      const int icv = cvocc_v[coc];
      rbicc_cv[icv] = rbi_i[icc-ncc_a];
    }
  }
  for (int icv = 0; icv < solver->ncv; ++icv)
    assert(rbicc_cv[icv] >= 0);
  solver->updateCvData(rbicc_cv);
  for (int icv = solver->ncv; icv < solver->ncv_g; ++icv) {
    const int icc = ccocv[icv];
    if (icc >= ncc) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbicc_cv[icv]);
      if (rank != mpi_rank) {
        if (rbi_g[icc-ncc] == -1) {
          rbi_g[icc-ncc] = rbicc_cv[icv];
        }
        // update with min rank rbi (active)
        else if (rbi_g[icc-ncc] != rbicc_cv[icv]) {
          int rank0,bits0,index0;
          BitUtils::unpackRankBitsIndex(rank0,bits0,index0,rbi_g[icc-ncc]); 
          //assert(rank != rank0); 
          if (rank < rank0)
            rbi_g[icc-ncc] = rbicc_cv[icv];
          else if ((rank == rank0)&&(index < index0))
            rbi_g[icc-ncc] = rbicc_cv[icv];
        }
      }
    }
  }
  delete[] rbicc_cv;

  delete[] icc_global; icc_global = new int8[ncc_g];
  for (map<const int8,int>::iterator iter = globalMap.begin(); iter != globalMap.end(); ++iter)
    icc_global[iter->second] = iter->first;

  map<const pair<int,int>,int> cfMap; // only b/w active/active or active/inactive cc's for now
  assert(ncf == 0);

  // active/active
  for (int ifa = 0; ifa < solver->nfa; ++ifa) {
    const int icv0 = solver->cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < solver->ncv));
    const int icv1 = solver->cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < solver->ncv_g));
    const int icc0 = ccocv[icv0]; assert((icc0 >= 0)&&(icc0 < ncc));
    const int icc1 = ccocv[icv1]; assert((icc1 >= 0)&&(icc1 < ncc_g));
    if ((icc0 != icc1)&&((icc0 < ncc_a)&&(icc1 < ncc_a))) {
      map<const pair<int,int>,int>::iterator iter = cfMap.find(pair<int,int>(min(icc0,icc1),max(icc0,icc1)));
      if (iter == cfMap.end()) 
        cfMap[pair<int,int>(min(icc0,icc1),max(icc0,icc1))] = ncf++;
    }
  }

  // active/inactive(ghost_ added to separate set because it will be updated before adding to cfMap
  set<pair<int,int> > cfSet_ai; 
  for (int ifa = 0; ifa < solver->nfa; ++ifa) {
    const int icv0 = solver->cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < solver->ncv));
    const int icv1 = solver->cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < solver->ncv_g));
    const int icc0 = ccocv[icv0]; assert((icc0 >= 0)&&(icc0 < ncc));
    const int icc1 = ccocv[icv1]; assert((icc1 >= 0)&&(icc1 < ncc_g));
    if ((icc0 != icc1)&&((icc0 < ncc_a)||(icc1 < ncc_a))) {
      map<const pair<int,int>,int>::iterator iter = cfMap.find(pair<int,int>(min(icc0,icc1),max(icc0,icc1)));
      if (iter == cfMap.end()) 
        cfSet_ai.insert(pair<int,int>(min(icc0,icc1),max(icc0,icc1)));
    }
  }

  // inactive/inactive(or ghost) added to separate set
  set<pair<int,int> > cfSet_ii; 
  for (int ifa = 0; ifa < solver->nfa; ++ifa) {
    const int icv0 = solver->cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < solver->ncv));
    const int icv1 = solver->cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < solver->ncv_g));
    const int icc0 = ccocv[icv0]; assert((icc0 >= 0)&&(icc0 < ncc));
    const int icc1 = ccocv[icv1]; assert((icc1 >= 0)&&(icc1 < ncc_g));
    if (icc0 != icc1) {
      map<const pair<int,int>,int>::iterator iter = cfMap.find(pair<int,int>(min(icc0,icc1),max(icc0,icc1)));
      set<pair<int,int> >::iterator iter2 = cfSet_ai.find(pair<int,int>(min(icc0,icc1),max(icc0,icc1)));
      if ((iter == cfMap.end())&&(iter2 == cfSet_ai.end())) {
        assert((icc0 >= ncc_a)&&(icc1 >= ncc_a)); // both inactive or ghost
        cfSet_ii.insert(pair<int,int>(min(icc0,icc1),max(icc0,icc1)));
      }
    }
  }

  int * ccocc0_i = new int[ncc+1]; // will expand once we have ghosts
  for (int icc = 0; icc < ncc; ++icc) 
    ccocc0_i[icc+1] = 1; // diagonal
  for (map<const pair<int,int>,int>::iterator iter = cfMap.begin(); iter != cfMap.end(); ++iter) {
    const int icc0 = iter->first.first; assert((icc0 >= 0)&&(icc0 < ncc_a));
    const int icc1 = iter->first.second; assert((icc1 >= 0)&&(icc1 < ncc));
    ++ccocc0_i[icc0+1];
    if (icc1 < ncc)
      ++ccocc0_i[icc1+1];
  }
  for (set<pair<int,int> >::iterator iter = cfSet_ai.begin(); iter != cfSet_ai.end(); ++iter) {
    const int icc0 = iter->first; assert((icc0 >= 0)&&(icc0 < ncc_a));
    const int icc1 = iter->second; assert((icc1 >= 0)&&(icc1 < ncc_g));
    ++ccocc0_i[icc0+1];
    if (icc1 < ncc)
      ++ccocc0_i[icc1+1];
  }
  for (set<pair<int,int> >::iterator iter = cfSet_ii.begin(); iter != cfSet_ii.end(); ++iter) {
    const int icc0 = iter->first; assert((icc0 >= ncc_a)&&(icc0 < ncc));
    const int icc1 = iter->second; assert((icc1 >= ncc_a)&&(icc1 < ncc_g));
    ++ccocc0_i[icc0+1];
    if (icc1 < ncc)
      ++ccocc0_i[icc1+1];
  }
  //  allocate...
  ccocc0_i[0] = 0;
  for (int icc = 0; icc < ncc; ++icc) ccocc0_i[icc+1] += ccocc0_i[icc];
  int *ccocc0_v = new int[ccocc0_i[ncc]]; // will expand once we have ghosts
  //  and set...
  for (int icc = 0; icc < ncc; ++icc) 
    ccocc0_v[ccocc0_i[icc]++] = icc; // diagonal
  for (map<const pair<int,int>,int>::iterator iter = cfMap.begin(); iter != cfMap.end(); ++iter) {
    const int icc0 = iter->first.first; assert((icc0 >= 0)&&(icc0 < ncc_a));
    const int icc1 = iter->first.second; assert((icc1 >= 0)&&(icc1 < ncc));
    ccocc0_v[ccocc0_i[icc0]++] = icc1;
    if (icc1 < ncc)
      ccocc0_v[ccocc0_i[icc1]++] = icc0;
  }
  for (set<pair<int,int> >::iterator iter = cfSet_ai.begin(); iter != cfSet_ai.end(); ++iter) {
    const int icc0 = iter->first; assert((icc0 >= 0)&&(icc0 < ncc_a));
    const int icc1 = iter->second; assert((icc1 >= 0)&&(icc1 < ncc_g));
    ccocc0_v[ccocc0_i[icc0]++] = icc1;
    if (icc1 < ncc)
      ccocc0_v[ccocc0_i[icc1]++] = icc0;
  }
  for (set<pair<int,int> >::iterator iter = cfSet_ii.begin(); iter != cfSet_ii.end(); ++iter) {
    const int icc0 = iter->first; assert((icc0 >= ncc_a)&&(icc0 < ncc));
    const int icc1 = iter->second; assert((icc1 >= ncc_a)&&(icc1 < ncc_g));
    ccocc0_v[ccocc0_i[icc0]++] = icc1;
    if (icc1 < ncc)
      ccocc0_v[ccocc0_i[icc1]++] = icc0;
  }
  //  rewind...    
  for (int icc = ncc; icc > 0; --icc)
    ccocc0_i[icc] = ccocc0_i[icc-1];
  ccocc0_i[0] = 0;

  // pack inactive icc active nbrs and send them to the associated active icc 

  FOR_RANK send_count[rank] = 0;
  for (int iter = 0; iter < 2; ++iter) {
    for (int icc = ncc_a; icc < ncc; ++icc) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_i[icc-ncc_a]);
      if (iter == 0) {
        send_count[rank] += 2+4*(ccocc0_i[icc+1]-ccocc0_i[icc]-1); // index,ncoc,{nbr,rank_nbr,nbr_striped,nbr_rank_striped}
      }
      else {
        send_buf_int[send_disp[rank]++] = index;
        send_buf_int[send_disp[rank]++] = ccocc0_i[icc+1]-ccocc0_i[icc]-1;
        for (int coc = ccocc0_i[icc]+1; coc < ccocc0_i[icc+1]; ++coc) {
          const int icc_nbr = ccocc0_v[coc]; 
          const int rank_striped = MiscUtils::getRankInXora(icc_global[icc_nbr],ccora_striped);
          if (icc_nbr >= ncc_a) {
            int rank_nbr,bits_nbr,index_nbr;
            if (icc_nbr < ncc) 
              BitUtils::unpackRankBitsIndex(rank_nbr,bits_nbr,index_nbr,rbi_i[icc_nbr-ncc_a]);
            else 
              BitUtils::unpackRankBitsIndex(rank_nbr,bits_nbr,index_nbr,rbi_g[icc_nbr-ncc]);
            send_buf_int[send_disp[rank]++] = index_nbr; 
            send_buf_int[send_disp[rank]++] = rank_nbr;
          }
          else {
            send_buf_int[send_disp[rank]++] = icc_nbr;
            send_buf_int[send_disp[rank]++] = mpi_rank;
          }
          send_buf_int[send_disp[rank]++] = icc_global[icc_nbr]-ccora_striped[rank_striped]; 
          send_buf_int[send_disp[rank]++] = rank_striped;
        }
      }
    }

    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

    if (iter == 0) {
      send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
      assert(send_buf_int == NULL); 
      send_buf_int = new int[send_count_sum];
    }
  }
  delete[] ccocc0_i;
  delete[] ccocc0_v;

  // setup recv side stuff

  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

  recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
  assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum];
  MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
      recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
  delete[] send_buf_int; send_buf_int = NULL;

  // first count unique ghosts,store ghost rbi's and add to set (there could be duplicates)

  const int ncc_g0 = ncc_g;
  set<int>* nbrSet = new set<int>[ncc_a]; 
  for (int ifa = 0; ifa < solver->nfa; ++ifa) {
    const int icv0 = solver->cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < solver->ncv));
    const int icv1 = solver->cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < solver->ncv_g));
    const int icc0 = ccocv[icv0]; assert((icc0 >= 0)&&(icc0 < ncc));
    const int icc1 = ccocv[icv1]; assert((icc1 >= 0)&&(icc1 < ncc_g));
    if (icc0 != icc1) {
      if (icc0 < ncc_a)
        nbrSet[icc0].insert(icc1);
      if (icc1 < ncc_a)
        nbrSet[icc1].insert(icc0);
    }
  }
  for (int iter = 0; iter < 2; ++iter) {
    FOR_RANK {
      int irecv = recv_disp[rank];
      while (irecv < recv_disp[rank]+recv_count[rank]) {
        const int icc = recv_buf_int[irecv++]; assert((icc >= 0)&&(icc < ncc_a));
        const int ncoc = recv_buf_int[irecv++];
        for (int coc = 0; coc < ncoc; ++coc) {
          const int icc_nbr = recv_buf_int[irecv++]; 
          const int rank_nbr = recv_buf_int[irecv++]; 
          const int8 rbi = BitUtils::packRankBitsIndex(rank_nbr,0,icc_nbr);
          const int icc_nbr_striped = recv_buf_int[irecv++];
          const int rank_striped = recv_buf_int[irecv++];
          const int8 icc_nbr_global = icc_nbr_striped+ccora_striped[rank_striped];
          assert(icc_nbr_global != icc_global[icc]);
          map<const int8,int>::iterator iter2 = globalMap.find(icc_nbr_global);
          if (iter == 0) {
            if (iter2 == globalMap.end()) {
              nbrSet[icc].insert(ncc_g);
              globalMap[icc_nbr_global] = ncc_g++;
            }
            else {
              nbrSet[icc].insert(iter2->second);
            }
          }
          else {
            if (iter2->second >= ncc) {
              if (rbi_g[iter2->second-ncc] == -1) {
                rbi_g[iter2->second-ncc] = rbi;
              }
              // update with min rank rbi (active)
              else if (rbi_g[iter2->second-ncc] != rbi) {
                int rank0,bits0,index0;
                BitUtils::unpackRankBitsIndex(rank0,bits0,index0,rbi_g[iter2->second-ncc]); 
                if (rank_nbr < rank0)
                  rbi_g[iter2->second-ncc] = rbi;
              }
            }
          }
        }
      }
    }

    if (iter == 0) {
      // add inactive rbi's 
      int8* rbi_g0 = rbi_g;
      rbi_g = new int8[ncc_g-ncc];
      for (int icc = ncc; icc < ncc_g0; ++icc)
        rbi_g[icc-ncc] = rbi_g0[icc-ncc];
      delete[] rbi_g0;
      for (int icc = ncc_g0; icc < ncc_g; ++icc) 
        rbi_g[icc-ncc] = -1;
    }
  }
  delete[] recv_buf_int; recv_buf_int = NULL;
  for (int icc = ncc; icc < ncc_g; ++icc) 
    assert(rbi_g[icc-ncc] != -1);

  // add in new ghosts to icc_global

  delete[] icc_global; icc_global = new int8[ncc_g];
  for (map<const int8,int>::iterator iter = globalMap.begin(); iter != globalMap.end(); ++iter)
    icc_global[iter->second] = iter->first;

  rbiVec.resize(ncc_g-ncc);
  for (int icc = ncc; icc < ncc_g; ++icc)
    rbiVec[icc-ncc] = pair<int8,int>(rbi_g[icc-ncc],icc);
  sort(rbiVec.begin(),rbiVec.end());

  // reorder ncc:ncc_g data
  assert(reorder_cc == NULL); reorder_cc = new int[ncc_g-ncc];
  int8* icc_global0 = icc_global;
  icc_global = new int8[ncc_g];
  for (int icc = 0; icc < ncc; ++icc) 
    icc_global[icc] = icc_global0[icc];
  for (int ii = 0, lim = ncc_g-ncc; ii < lim; ++ii) {
    const int icc = ii+ncc; 
    const int icc_old = rbiVec[ii].second; assert((icc_old >= ncc)&&(icc_old < ncc_g));
    map<const int8,int>::iterator iter = globalMap.find(icc_global0[icc_old]);
    assert(iter != globalMap.end());
    iter->second = icc;
    rbi_g[ii] = rbiVec[ii].first;
    reorder_cc[icc_old-ncc] = icc;
    icc_global[icc] = iter->first;
  }
  delete[] icc_global0;
  rbiVec.clear();

  // update ccocv and cvocc_i/v
  for (int icv = 0; icv < solver->ncv_g; ++icv) {
    if (ccocv[icv] >= ncc)
      ccocv[icv] = reorder_cc[ccocv[icv]-ncc];
  }
  for (int icc = 0; icc < ncc; ++icc) cvocc_i[icc+1] = 0;
  for (int iter = 0; iter < 2; ++iter) {
    for (int icv = 0; icv < solver->ncv; ++icv) {
      const int icc = ccocv[icv];
      if (iter == 0) 
        ++cvocc_i[icc+1];
      else 
        cvocc_v[cvocc_i[icc]++] = icv;
    }
    if (iter == 0) {
      cvocc_i[0] = 0;
      for (int icc = 0; icc < ncc; ++icc) cvocc_i[icc+1] += cvocc_i[icc];
    }
    else {
      // rewind
      for (int icc = ncc; icc > 0; --icc)
        cvocc_i[icc] = cvocc_i[icc-1];
      cvocc_i[0] = 0;
    }
  }

  // rbi check...
  {
    FOR_RANK send_count[rank] = 0;
    for (int iter = 0; iter < 2; ++iter) {
      for (int icc = ncc_a; icc < ncc_g; ++icc) {
        int rank,bits,index;
        if (icc < ncc)
          BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_i[icc-ncc_a]); 
        else
          BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icc-ncc]); 
        if (iter == 0) {
          ++send_count[rank];
        }
        else {
          send_buf_int[send_disp[rank]++] = index;
        }
      }
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      if (iter == 0) {
        send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
        assert(send_buf_int == NULL); 
        send_buf_int = new int[send_count_sum];
      }
    }

    // setup recv side stuff

    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

    recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
        recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; send_buf_int = NULL;

    for (int irecv = 0; irecv < recv_count_sum; ++irecv)
      assert(recv_buf_int[irecv] < ncc_a);
    delete[] recv_buf_int; recv_buf_int = NULL;
  }

  // build ccocc using nbrSet

  assert(ccocc_i == NULL); ccocc_i = new int[ncc_a+1];
  for (int icc = 0; icc < ncc_a; ++icc) 
    ccocc_i[icc+1] = 1+nbrSet[icc].size();
  ccocc_i[0] = 0;
  for (int icc = 0; icc < ncc_a; ++icc) ccocc_i[icc+1] += ccocc_i[icc];
  assert(ccocc_v == NULL); ccocc_v = new int[ccocc_i[ncc_a]];
  for (int icc = 0; icc < ncc_a; ++icc) {
    ccocc_v[ccocc_i[icc]++] = icc;
    for (set<int>::iterator iter = nbrSet[icc].begin(); iter != nbrSet[icc].end(); ++iter) {
      if (*iter < ncc)
        ccocc_v[ccocc_i[icc]++] = *iter;
      else 
        ccocc_v[ccocc_i[icc]++] = reorder_cc[*iter-ncc];
    }
    nbrSet[icc].clear();
  }
  delete[] nbrSet;
  for (int icc = ncc_a; icc > 0; --icc)
    ccocc_i[icc] = ccocc_i[icc-1];
  ccocc_i[0] = 0;

  // add active/active connections from other ranks
  for (int icc = 0; icc < ncc_a; ++icc) {
    for (int coc = ccocc_i[icc]+1; coc < ccocc_i[icc+1]; ++coc) {
      const int icc_nbr = ccocc_v[coc];
      if (icc_nbr < ncc_a) {
        map<const pair<int,int>,int>::iterator iter = cfMap.find(pair<int,int>(min(icc,icc_nbr),max(icc,icc_nbr)));
        if (iter == cfMap.end()) 
          cfMap[pair<int,int>(min(icc,icc_nbr),max(icc,icc_nbr))] = ncf++;
      }
    }
  }
  assert(ncf_aa == 0); ncf_aa = ncf;

  // add active/inactive(ghost) connections from cfSet_ai
  for (set<pair<int,int> >::iterator iter = cfSet_ai.begin(); iter != cfSet_ai.end(); ++iter) {
    const int icc0 = iter->first; assert((icc0 >= 0)&&(icc0 < ncc));
    int icc1 = iter->second; assert((icc1 >= 0)&&(icc1 < ncc_g));
    if (icc1 >= ncc) 
      icc1 = reorder_cc[icc1-ncc];
    assert(icc1 > icc0);
    map<const pair<int,int>,int>::iterator iter2 = cfMap.find(pair<int,int>(icc0,icc1));
    if (iter2 == cfMap.end()) 
      cfMap[pair<int,int>(icc0,icc1)] = ncf++;
  }
  cfSet_ai.clear();

  // add active/ghost from other ranks
  for (int icc = 0; icc < ncc_a; ++icc) {
    for (int coc = ccocc_i[icc]+1; coc < ccocc_i[icc+1]; ++coc) {
      const int icc_nbr = ccocc_v[coc];
      if (icc_nbr >= ncc_a) {
        map<const pair<int,int>,int>::iterator iter = cfMap.find(pair<int,int>(min(icc,icc_nbr),max(icc,icc_nbr)));
        if (iter == cfMap.end()) 
          cfMap[pair<int,int>(min(icc,icc_nbr),max(icc,icc_nbr))] = ncf++;
      }
    }
  }
  assert(ncf_ai == 0); ncf_ai = ncf;

  // throw the inactive/inactive(or ghost) faces at the end
  for (set<pair<int,int> >::iterator iter = cfSet_ii.begin(); iter != cfSet_ii.end(); ++iter) {
    const int icc0 = iter->first; assert((icc0 >= 0)&&(icc0 < ncc));
    int icc1 = iter->second; assert((icc1 >= 0)&&(icc1 < ncc_g));
    if (icc1 >= ncc)
      icc1 = reorder_cc[icc1-ncc];
    assert(icc1 > icc0);
    map<const pair<int,int>,int>::iterator iter2 = cfMap.find(pair<int,int>(icc0,icc1));
    if (iter2 == cfMap.end()) 
      cfMap[pair<int,int>(icc0,icc1)] = ncf++;
  }
  delete[] reorder_cc;
  cfSet_ii.clear();

  // build ccocf

  assert(ccocf == NULL); ccocf = new int[ncf][2];
  for (map<const pair<int,int>,int>::iterator iter = cfMap.begin(); iter != cfMap.end(); ++iter) {
    const int icc0 = iter->first.first;
    const int icc1 = iter->first.second;
    const int icf = iter->second;
    ccocf[icf][0] = icc0;
    ccocf[icf][1] = icc1;
  }

  // build ccIPrcomm so we can check it...

  buildCcIPrcomm();

  // cc geom data

  assert(x_cc == NULL); x_cc = new double[ncc_g][3];
  assert(vol_cc == NULL); vol_cc = new double[ncc_g];
  for (int icc = 0; icc < ncc; ++icc) {
    FOR_I3 x_cc[icc][i] = 0.0;
    vol_cc[icc] = 0.0;
    for (int coc = cvocc_i[icc]; coc != cvocc_i[icc+1]; ++coc) {
      const int icv = cvocc_v[coc];
      FOR_I3 x_cc[icc][i] += solver->vol_cv[icv]*solver->x_cv[icv][i];
      vol_cc[icc] += solver->vol_cv[icv];
    }
  }

  double * vol_cc_check = new double[ncc];
  double (*x_cc_check)[3] = new double[ncc][3];
  for (int icc = 0; icc < ncc; ++icc) {
    vol_cc_check[icc] = vol_cc[icc];
    FOR_I3 x_cc_check[icc][i] = x_cc[icc][i];
  }
  DistributedDataExchanger * dde_cc_striped = new DistributedDataExchanger(icc_global,ncc,ccora_striped);
  double (*x_ccs)[3] = new double[nccs][3];
  for (int iccs = 0; iccs < nccs; ++iccs)
    FOR_I3 x_ccs[iccs][i] = 0.0;
  dde_cc_striped->push(x_ccs,x_cc_check,ADD_DATA);
  dde_cc_striped->pull(x_cc_check,x_ccs);
  delete[] x_ccs;
  double* vol_ccs = new double[nccs];
  for (int iccs = 0; iccs < nccs; ++iccs)
    vol_ccs[iccs] = 0.0;
  dde_cc_striped->push(vol_ccs,vol_cc_check,ADD_DATA);
  dde_cc_striped->pull(vol_cc_check,vol_ccs);
  delete dde_cc_striped;
  delete[] vol_ccs;
  for (int icc = 0; icc < ncc; ++icc) {
    assert(vol_cc_check[icc] > 0);
    FOR_I3 x_cc_check[icc][i] /= vol_cc_check[icc];
  }

  updateCcIData(vol_cc);
  updateCcIData(x_cc);
  for (int icc = 0; icc < ncc_a; ++icc) {
    assert(vol_cc[icc] > 0.0);
    FOR_I3 x_cc[icc][i] /= vol_cc[icc];
  }

  for (int icc = 0; icc < ncc_a; ++icc) {
    FOR_I3 x_cc_check[icc][i] -= x_cc[icc][i];
    vol_cc_check[icc] -= vol_cc[icc];
  }
  dumpRange(vol_cc_check,ncc_a,"vol_cc check");
  delete[] vol_cc_check;
  dumpRange(x_cc_check,ncc_a,"x_cc check");

  restrictCcData(x_cc_check,solver->x_cv);
  for (int icc = 0; icc < ncc_a; ++icc) 
    FOR_I3 x_cc_check[icc][i] -= x_cc[icc][i];
  dumpRange(x_cc_check,ncc_a,"x_cc check2");
  delete[] x_cc_check;

  double my_buf[8]; FOR_I8 my_buf[i] = 0;
  for (int icv = 0; icv < solver->ncv; ++icv) {
    my_buf[0] += solver->vol_cv[icv];
    FOR_I3 my_buf[1+i] += solver->vol_cv[icv]*solver->x_cv[icv][i];
  }
  for (int icc = 0; icc < ncc_a; ++icc) {
    my_buf[4] += vol_cc[icc];
    FOR_I3 my_buf[5+i] += vol_cc[icc]*x_cc[icc][i];
  }
  double buf[8];
  MPI_Reduce(my_buf,buf,8,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if (mpi_rank == 0) {
    FOR_I3 buf[1+i] /= buf[0];
    FOR_I3 buf[5+i] /= buf[4];
    cout << " > vol diff: " << buf[0]-buf[4] << ", centroid diff: " << buf[5]-buf[1] << " " << buf[6]-buf[2] << " " << buf[7]-buf[3] << endl;
  }

  buildCcPrcomm();
  for (int icc = ncc_a; icc < ncc_g; ++icc) {
    vol_cc[icc] = 1.0E+20;
    FOR_I3 x_cc[icc][i] = 1.0E+20;
  }
  updateCcData(vol_cc);
  updateCcData(x_cc); // no bits for now
  dumpRange(vol_cc,ncc_g,"vol_cc");
  dumpRange(x_cc,ncc_a,"x_cc");

  // cf geom data

  assert(x_cf == NULL); x_cf = new double[ncf][3];
  assert(n_cf == NULL); n_cf = new double[ncf][3]; 
  assert(area_cf == NULL); area_cf = new double[ncf]; 
  for (int icf = 0; icf < ncf; ++icf) {
    FOR_I3 x_cf[icf][i] = 0.0;
    FOR_I3 n_cf[icf][i] = 0.0;
    area_cf[icf] = 0.0;
  }

  // store local geometric info
  assert(faocf_i == NULL); faocf_i = new int[ncf+1]; 
  for (int icf = 0; icf < ncf; ++icf) faocf_i[icf+1] = 0;
  for (int iter = 0; iter < 2; ++iter) {
    for (int ifa = 0; ifa < solver->nfa; ++ifa) {
      // use ifa_global to prevent double counting at proc boundaries
      if (solver->ifa_global[ifa] >= 0) {
        const int icv0 = solver->cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < solver->ncv));
        const int icv1 = solver->cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < solver->ncv_g));
        const int icc0 = ccocv[icv0]; assert((icc0 >= 0)&&(icc0 < ncc));
        const int icc1 = ccocv[icv1]; assert((icc1 >= 0)&&(icc1 < ncc_g));
        if (icc0 != icc1) {
          map<const pair<int,int>,int>::iterator iter2 = cfMap.find(pair<int,int>(min(icc0,icc1),max(icc0,icc1)));
          assert(iter2 != cfMap.end());
          const int icf = iter2->second; assert((icf >= 0)&&(icf < ncf));
          if (iter == 0) {
            ++faocf_i[icf+1];
          }
          else {
            FOR_I3 x_cf[icf][i] += solver->area_fa[ifa]*solver->x_fa[ifa][i];
            double cf_sign;
            if (icc0 > icc1) {
              cf_sign = -1.0;
              faocf_v[faocf_i[icf]++] = -ifa-1;
            }
            else {
              cf_sign = 1.0;
              faocf_v[faocf_i[icf]++] = ifa;
            }
            FOR_I3 n_cf[icf][i] += cf_sign*solver->n_fa[ifa][i];
            area_cf[icf] += solver->area_fa[ifa];
          }
        }
      }
    }
    if (iter == 0) {
      // index from counts
      faocf_i[0] = 0;
      for (int icf = 0; icf < ncf; ++icf) 
        faocf_i[icf+1] += faocf_i[icf];
      assert(solver->nfa >= faocf_i[ncf]); // all the faces are not touched
      assert(faocf_v == NULL); faocf_v = new int[faocf_i[ncf]];
    }
    else {
      // rewind
      for (int icf = ncf; icf > 0; --icf)
        faocf_i[icf] = faocf_i[icf-1];
      faocf_i[0] = 0;
    }
  }

  FOR_RANK send_count[rank] = 0;
  double *send_buf_double = NULL;
  int *cf_send = NULL;
  for (int iter = 0; iter < 2; ++iter) {
    for (int icf = ncf_aa; icf < ncf; ++icf) {
      const int icc0 = ccocf[icf][0]; assert((icc0 >= 0)&&(icc0 < ncc));
      const int icc1 = ccocf[icf][1]; assert((icc1 >= 0)&&(icc1 < ncc_g));
      const int rank0_striped = MiscUtils::getRankInXora(icc_global[icc0],ccora_striped);
      const int icc0_striped = icc_global[icc0]-ccora_striped[rank0_striped];
      const int rank1_striped = MiscUtils::getRankInXora(icc_global[icc1],ccora_striped);
      const int icc1_striped = icc_global[icc1]-ccora_striped[rank1_striped];
      //inactive/inactive and inactive/ghost
      if ((icc0 >= ncc_a)&&(icc1 >= ncc_a)) {
        int rank0,bits0,index0;
        if (icc0 < ncc)
          BitUtils::unpackRankBitsIndex(rank0,bits0,index0,rbi_i[icc0-ncc_a]);
        else 
          BitUtils::unpackRankBitsIndex(rank0,bits0,index0,rbi_g[icc0-ncc]);
        int rank1,bits1,index1;
        if (icc1 < ncc)
          BitUtils::unpackRankBitsIndex(rank1,bits1,index1,rbi_i[icc1-ncc_a]);
        else
          BitUtils::unpackRankBitsIndex(rank1,bits1,index1,rbi_g[icc1-ncc]);
        if (rank0 == rank1) {
          if (iter == 0) {
            ++send_count[rank0];
          }
          else {
            send_buf_int[send_disp[rank0]*3+0] = index0;
            send_buf_int[send_disp[rank0]*3+1] = rank1_striped;
            send_buf_int[send_disp[rank0]*3+2] = icc1_striped;
            FOR_I3 send_buf_double[send_disp[rank0]*7+i] = x_cf[icf][i]; // area weighted
            FOR_I3 send_buf_double[send_disp[rank0]*7+3+i] = n_cf[icf][i];
            send_buf_double[send_disp[rank0]*7+6] = area_cf[icf];
            cf_send[send_disp[rank0]] = icf;
            ++send_disp[rank0];
          }
        }
        else {
          if (iter == 0) {
            ++send_count[rank0];
            ++send_count[rank1];
          }
          else {
            send_buf_int[send_disp[rank0]*3+0] = index0;
            send_buf_int[send_disp[rank0]*3+1] = rank1_striped;
            send_buf_int[send_disp[rank0]*3+2] = icc1_striped;
            FOR_I3 send_buf_double[send_disp[rank0]*7+i] = x_cf[icf][i]; 
            FOR_I3 send_buf_double[send_disp[rank0]*7+3+i] = n_cf[icf][i];
            send_buf_double[send_disp[rank0]*7+6] = area_cf[icf];
            cf_send[send_disp[rank0]] = icf;
            ++send_disp[rank0];
            send_buf_int[send_disp[rank1]*3+0] = index1;
            send_buf_int[send_disp[rank1]*3+1] = rank0_striped;
            send_buf_int[send_disp[rank1]*3+2] = icc0_striped;
            FOR_I3 send_buf_double[send_disp[rank1]*7+i] = x_cf[icf][i];
            FOR_I3 send_buf_double[send_disp[rank1]*7+3+i] = -n_cf[icf][i];
            send_buf_double[send_disp[rank1]*7+6] = area_cf[icf];
            cf_send[send_disp[rank1]] = -icf-1;
            ++send_disp[rank1];
          }
        }
      }
      // icc0 inactive, icc1 active 
      else if (icc0 >= ncc_a) {
        int rank0,bits0,index0;
        if (icc0 < ncc)
          BitUtils::unpackRankBitsIndex(rank0,bits0,index0,rbi_i[icc0-ncc_a]);
        else 
          BitUtils::unpackRankBitsIndex(rank0,bits0,index0,rbi_g[icc0-ncc]);
        if (iter == 0) {
          ++send_count[rank0];
        }
        else {
          send_buf_int[send_disp[rank0]*3+0] = index0;
          send_buf_int[send_disp[rank0]*3+1] = rank1_striped;
          send_buf_int[send_disp[rank0]*3+2] = icc1_striped;
          FOR_I3 send_buf_double[send_disp[rank0]*7+i] = x_cf[icf][i];
          FOR_I3 send_buf_double[send_disp[rank0]*7+3+i] = n_cf[icf][i];
          send_buf_double[send_disp[rank0]*7+6] = area_cf[icf];
          cf_send[send_disp[rank0]] = icf;
          ++send_disp[rank0];
        }
      }
      // icc0 active, icc1 inactive
      else {
        assert(icc1 >= ncc_a);
        int rank1,bits1,index1;
        if (icc1 < ncc)
          BitUtils::unpackRankBitsIndex(rank1,bits1,index1,rbi_i[icc1-ncc_a]);
        else
          BitUtils::unpackRankBitsIndex(rank1,bits1,index1,rbi_g[icc1-ncc]);
        if (iter == 0) {
          ++send_count[rank1];
        }
        else {
          send_buf_int[send_disp[rank1]*3+0] = index1;
          send_buf_int[send_disp[rank1]*3+1] = rank0_striped;
          send_buf_int[send_disp[rank1]*3+2] = icc0_striped;
          FOR_I3 send_buf_double[send_disp[rank1]*7+i] = x_cf[icf][i];
          FOR_I3 send_buf_double[send_disp[rank1]*7+3+i] = -n_cf[icf][i];
          send_buf_double[send_disp[rank1]*7+6] = area_cf[icf];
          cf_send[send_disp[rank1]] = -icf-1;
          ++send_disp[rank1];
        }
      }
    }

    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

    if (iter == 0) {
      send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
      assert(send_buf_int == NULL); 
      send_buf_int = new int[3*send_count_sum];
      cf_send = new int[send_count_sum];
      assert(send_buf_double == NULL); 
      send_buf_double = new double[7*send_count_sum];
    }

  }

  // setup recv side stuff

  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

  recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
  assert(recv_buf_int == NULL); recv_buf_int = new int[3*recv_count_sum];
  FOR_RANK {
    send_count[rank] *= 3;
    send_disp[rank] *= 3;
    recv_count[rank] *= 3;
    recv_disp[rank] *= 3;
  }
  MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
      recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);

  double *recv_buf_double = new double[7*recv_count_sum];
  FOR_RANK {
    send_count[rank] = (send_count[rank]/3)*7;
    send_disp[rank] = (send_disp[rank]/3)*7;
    recv_count[rank] = (recv_count[rank]/3)*7;
    recv_disp[rank] = (recv_disp[rank]/3)*7;
  }
  MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
      recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
  delete[] send_buf_double; send_buf_double = NULL;

  FOR_RANK {
    send_count[rank] /= 7;
    send_disp[rank] /= 7;
    recv_count[rank] /= 7;
    recv_disp[rank] /= 7;
  }

  for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
    const int icc0 = recv_buf_int[irecv*3+0]; assert((icc0 >= 0)&&(icc0 < ncc_a));
    const int rank1_striped = recv_buf_int[irecv*3+1];
    const int icc1_striped = recv_buf_int[irecv*3+2];
    const int8 icc1_global = icc1_striped+ccora_striped[rank1_striped]; assert((icc1_global >= 0)&&(icc1_global < ncc_global));
    map<const int8,int>::iterator iter1 = globalMap.find(icc1_global);
    assert(iter1 != globalMap.end());
    const int icc1 = iter1->second; assert((icc1 >= 0)&&(icc1 < ncc_g));
    map<const pair<int,int>,int>::iterator iter = cfMap.find(pair<int,int>(min(icc0,icc1),max(icc0,icc1)));
    assert(iter != cfMap.end());
    const int icf = iter->second; 
    assert((icf >= 0)&&(icf < ncf_ai));
    FOR_I3 x_cf[icf][i] += recv_buf_double[irecv*7+i];
    double cf_sign;
    if (icc0 > icc1) {
      cf_sign = -1.0;
      recv_buf_int[irecv] = -icf-1;
    }
    else {
      cf_sign = 1.0;
      recv_buf_int[irecv] = icf;
    }
    FOR_I3 n_cf[icf][i] += cf_sign*recv_buf_double[irecv*7+3+i];
    area_cf[icf] += recv_buf_double[irecv*7+6];
  }
  cfMap.clear();
  globalMap.clear();
  delete[] recv_buf_double; recv_buf_double = NULL;
  delete[] icc_global_cv;
  delete[] ccora_striped;

  // normalize x_cf and build delta_cf based on distance b/w centroids in normal direction

  assert(delta_cf == NULL); delta_cf = new double[ncf]; 
  for (int icf = 0; icf < ncf_ai; ++icf) {
    const int icc0 = ccocf[icf][0]; assert((icc0 >= 0)&&(icc0 < ncc));
    const int icc1 = ccocf[icf][1]; assert((icc1 >= 0)&&(icc1 < ncc_g));
    assert(icc0 < icc1);
    assert(area_cf[icf] > 0.0);
    FOR_I3 x_cf[icf][i] /= area_cf[icf];
    delta_cf[icf] = DIST(x_cc[icc1],x_cc[icc0]); 
    //const double proj_area = MAG(n_cf[icf]);
    //if (proj_area > 1.0E-20) {
    //  const double dx[3] = DIFF(x_cc[icc1],x_cc[icc0]);
    // is this correct?
    //  delta_cf[icf] = fabs(DOT_PRODUCT(dx,n_cf[icf])/proj_area); // projected length
    //}
    //else {
    //  delta_cf[icf] = DIST(x_cc[icc1],x_cc[icc0]); // total length (needed?)
    //}
  }
  dumpRange(delta_cf,ncf_ai,"delta_cf");
  dumpRange(area_cf,ncf_ai,"area_cf");
  dumpRange(n_cf,ncf_ai,"n_cf");
  dumpRange(x_cf,ncf_ai,"x_cf");

  // send back the icf's for the fa data we used to build communicator

  MPI_Alltoallv(recv_buf_int,recv_count,recv_disp,MPI_INT,
      send_buf_int,send_count,send_disp,MPI_INT,mpi_comm);
  delete[] recv_buf_int; recv_buf_int = NULL;

  set<int8>* rbiSet_cf_i = new set<int8>[ncf-ncf_aa]; 
  FOR_RANK {
    int isend = send_disp[rank];
    while (isend < send_disp[rank]+send_count[rank]) {
      int icf = cf_send[isend]; 
      int icf_nbr = send_buf_int[isend];
      double cf_sign = 1.0;
      if (icf < 0) {
        icf = -icf-1;
        cf_sign *= -1.0;
      }
      if (icf_nbr < 0) {
        icf_nbr = -icf_nbr-1;
        cf_sign *= -1.0;
      }
      assert((icf >= ncf_aa)&&(icf < ncf));
      assert(icf_nbr >= 0);
      const int8 rbi = BitUtils::packRankBitsIndex(rank,0,icf_nbr);
      if (cf_sign >= 0.0) 
        rbiSet_cf_i[icf-ncf_aa].insert(rbi);
      else
        rbiSet_cf_i[icf-ncf_aa].insert(-rbi-1);
      ++isend;
    }
  }
  delete[] cf_send;
  delete[] send_buf_int; send_buf_int = NULL;

  assert(rbiVec_cf_i == NULL); rbiVec_cf_i = new vector<int8>[ncf-ncf_aa];
  for (int icf = ncf_aa; icf < ncf; ++icf) {
    for (set<int8>::iterator iter = rbiSet_cf_i[icf-ncf_aa].begin(); iter != rbiSet_cf_i[icf-ncf_aa].end(); ++iter)
      rbiVec_cf_i[icf-ncf_aa].push_back(*iter);
  }
  delete[] rbiSet_cf_i;

  buildCfIPrcomm();

  double* area_cf_check = new double[ncf];
  for (int icf = 0; icf < ncf; ++icf) {
    area_cf_check[icf] = 0.0;
    for (int fof = faocf_i[icf]; fof < faocf_i[icf+1]; ++fof) {
      const int ifa = max(faocf_v[fof],-faocf_v[fof]-1);;
      assert(solver->ifa_global[ifa] >= 0);
      area_cf_check[icf] += solver->area_fa[ifa];
    }
  }
  updateCfIData(area_cf_check);
  for (int icf = 0; icf < ncf_ai; ++icf) 
    area_cf_check[icf] -= area_cf[icf];
  dumpRange(area_cf_check,ncf_ai,"area_cf check");
  delete[] area_cf_check;

  double (*x_cf_check)[3] = new double[ncf][3];
  restrictCfData(x_cf_check,solver->x_fa);
  for (int icf = 0; icf < ncf_ai; ++icf) 
    FOR_I3 x_cf_check[icf][i] -= x_cf[icf][i];
  dumpRange(x_cf_check,ncf_ai,"x_cf check");
  delete[] x_cf_check;

  double (*n_cf_check)[3] = new double[ncf][3];
  for (int icf = 0; icf < ncf; ++icf) {
    FOR_I3 n_cf_check[icf][i] = 0.0;
    for (int fof = faocf_i[icf]; fof < faocf_i[icf+1]; ++fof) {
      int ifa = faocf_v[fof];
      double fa_sign = 1.0;
      if (ifa < 0) {
        ifa = -ifa-1;
        fa_sign = -1.0;
      }
      assert((ifa >= 0)&&(ifa < solver->nfa));
      assert(solver->ifa_global[ifa] >= 0);
      FOR_I3 n_cf_check[icf][i] += fa_sign*solver->n_fa[ifa][i];
    }
  }
  updateSignedCfIData(n_cf_check);
  for (int icf = 0; icf < ncf_ai; ++icf) 
    FOR_I3 n_cf_check[icf][i] -= n_cf[icf][i];
  dumpRange(n_cf_check,ncf_ai,"n_cf check");
  delete[] n_cf_check;

  // boundary data...

  map<const pair<int,int>,int> cbMap; // zone,icc -> icb
  assert(ncb == 0);
  for (int ibf = 0; ibf < solver->nbf; ++ibf) {
    const int izone = solver->zone_bf[ibf];
    const int icc = ccocv[solver->cvobf[ibf]];
    if (icc < ncc_in) { // completely internal cells
      map<const pair<int,int>,int>::iterator iter = cbMap.find(pair<int,int>(izone,icc));
      if (iter == cbMap.end())
        cbMap[pair<int,int>(izone,icc)] = ncb++;
    }
  }
  map<const pair<int,int>,int> cbMap_i;
  int ncb_i = 0;
  for (int ibf = 0; ibf < solver->nbf; ++ibf) {
    const int izone = solver->zone_bf[ibf];
    const int icc = ccocv[solver->cvobf[ibf]];
    if (icc >= ncc_in) { // split cells
      if (icc < ncc_a) { // active cells
        map<const pair<int,int>,int>::iterator iter = cbMap.find(pair<int,int>(izone,icc));
        if (iter == cbMap.end())
          cbMap[pair<int,int>(izone,icc)] = ncb++;
      }
      else {
        map<const pair<int,int>,int>::iterator iter = cbMap_i.find(pair<int,int>(izone,icc));
        if (iter == cbMap_i.end())
          cbMap_i[pair<int,int>(izone,icc)] = ncb_i++;
      }
    }
  }

  FOR_RANK send_count[rank] = 0;
  int *cb_send;
  for (int iter = 0; iter < 2; ++iter) {
    for (int ibf = 0; ibf < solver->nbf; ++ibf) {
      const int icc = ccocv[solver->cvobf[ibf]];
      if (icc >= ncc_a) { // inactive/ghost cells
        int rank,bits,index;
        assert(icc < ncc);
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_i[icc-ncc_a]);

        if (iter == 0) {
          ++send_count[rank];
        }
        else {
          FOR_I3 send_buf_double[send_disp[rank]*7+i] = solver->x_bf[ibf][i]*solver->area_bf[ibf];
          FOR_I3 send_buf_double[send_disp[rank]*7+3+i] = solver->n_bf[ibf][i];
          send_buf_double[send_disp[rank]*7+6] = solver->area_bf[ibf];
          send_buf_int[send_disp[rank]*2  ] = index;
          send_buf_int[send_disp[rank]*2+1] = solver->zone_bf[ibf];
          map<const pair<int,int>,int>::iterator iter2 = cbMap_i.find(pair<int,int>(solver->zone_bf[ibf],icc));
          assert(iter2 != cbMap_i.end());
          const int icb_i = iter2->second; 
          cb_send[send_disp[rank]] = icb_i;
          ++send_disp[rank];
        }
      }
    }

    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

    if (iter == 0) {
      send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
      assert(send_buf_double == NULL);
      send_buf_double = new double[7*send_count_sum];
      assert(send_buf_int == NULL);
      send_buf_int = new int[2*send_count_sum];
      cb_send = new int[send_count_sum];
    }

  }

  // setup recv side stuff

  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

  recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
  assert(recv_buf_int == NULL); recv_buf_int = new int[2*recv_count_sum];
  FOR_RANK {
    send_count[rank] *= 2;
    send_disp[rank] *= 2;
    recv_count[rank] *= 2;
    recv_disp[rank] *= 2;
  }
  MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
      recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);

  assert(recv_buf_double == NULL); recv_buf_double = new double[7*recv_count_sum];
  FOR_RANK {
    send_count[rank] = (send_count[rank]/2)*7;
    send_disp[rank] = (send_disp[rank]/2)*7;
    recv_count[rank] = (recv_count[rank]/2)*7;
    recv_disp[rank] = (recv_disp[rank]/2)*7;
  }
  MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
      recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
  FOR_RANK {
    send_count[rank] /= 7;
    send_disp[rank] /= 7;
    recv_count[rank] /= 7;
    recv_disp[rank] /= 7;
  }
  for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
    const int icc = recv_buf_int[irecv*2  ]; assert((icc >= ncc_in)&&(icc < ncc_a)); // active!
    const int izone = recv_buf_int[irecv*2+1]; assert((izone >= 0)&&(izone < int(solver->bfZoneVec.size())));
    map<const pair<int,int>,int>::iterator iter = cbMap.find(pair<int,int>(izone,icc));
    if (iter == cbMap.end())
      cbMap[pair<int,int>(izone,icc)] = ncb++;
  }
  assert(ncb_a == 0); ncb_a = ncb;
  ncb += ncb_i;

  // reorder active based on zone, inactive will be ordered by rbi_cb after
  int *reorder_cb = new int[ncb];
  for (int icb = 0; icb < ncb; ++icb) 
    reorder_cb[icb] = -1;
  ncb = 0;
  for (map<const pair<int,int>,int>::iterator iter = cbMap.begin(); iter != cbMap.end(); ++iter) 
    reorder_cb[iter->second] = ncb++; // zone/icc order
  assert(ncb == ncb_a);
  ncb += ncb_i;

  // should have all the active boundary faces now so allocate the data

  assert(area_cb == NULL); area_cb = new double[ncb];
  assert(n_cb == NULL); n_cb = new double[ncb][3];
  assert(x_cb == NULL); x_cb = new double[ncb][3];
  assert(ccocb == NULL); ccocb = new int[ncb];
  assert(zone_cb == NULL); zone_cb = new int[ncb];
  for (int icb = 0; icb < ncb; ++icb) {
    area_cb[icb] = 0.0;
    FOR_I3 n_cb[icb][i] = 0.0;
    FOR_I3 x_cb[icb][i] = 0.0;
    ccocb[icb] = -1;
    zone_cb[icb] = -1;
  }
  for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
    const int icc = recv_buf_int[irecv*2  ]; assert((icc >= ncc_in)&&(icc < ncc_a));
    const int izone = recv_buf_int[irecv*2+1]; assert((izone >= 0)&&(izone < int(solver->bfZoneVec.size())));
    map<const pair<int,int>,int>::iterator iter = cbMap.find(pair<int,int>(izone,icc));
    assert(iter != cbMap.end());
    const int icb = reorder_cb[iter->second]; assert((icb >= 0)&&(icb < ncb_a));
    FOR_I3 x_cb[icb][i] += recv_buf_double[irecv*7+i];
    FOR_I3 n_cb[icb][i] += recv_buf_double[irecv*7+3+i];
    area_cb[icb] += recv_buf_double[irecv*7+6];
    if (zone_cb[icb] == -1)
      zone_cb[icb] = izone;
    else
      assert(zone_cb[icb] == izone);
    if (ccocb[icb] == -1)
      ccocb[icb] = icc;
    else
      assert(ccocb[icb] == icc);
  }

  // send back x_cb,n_cb,area_cb and rbi_cb

  for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
    const int icc = recv_buf_int[irecv*2  ]; assert((icc >= ncc_in)&&(icc < ncc_a));
    const int izone = recv_buf_int[irecv*2+1]; assert((izone >= 0)&&(izone < int(solver->bfZoneVec.size())));
    map<const pair<int,int>,int>::iterator iter = cbMap.find(pair<int,int>(izone,icc));
    assert(iter != cbMap.end());
    const int icb = reorder_cb[iter->second]; assert((icb >= 0)&&(icb < ncb_a));
    FOR_I3 recv_buf_double[irecv*7+i] = x_cb[icb][i];
    FOR_I3 recv_buf_double[irecv*7+3+i] = n_cb[icb][i];
    recv_buf_double[irecv*7+6] = area_cb[icb];
    recv_buf_int[irecv] = icb;
  }
  MPI_Alltoallv(recv_buf_int,recv_count,recv_disp,MPI_INT,
      send_buf_int,send_count,send_disp,MPI_INT,mpi_comm);
  delete[] recv_buf_int; recv_buf_int = NULL;

  FOR_RANK {
    send_count[rank] *= 7;
    send_disp[rank] *= 7;
    recv_count[rank] *= 7;
    recv_disp[rank] *= 7;
  }
  MPI_Alltoallv(recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,
      send_buf_double,send_count,send_disp,MPI_DOUBLE,mpi_comm);
  delete[] recv_buf_double; recv_buf_double = NULL;

  FOR_RANK {
    send_count[rank] /= 7;
    send_disp[rank] /= 7;
    recv_count[rank] /= 7;
    recv_disp[rank] /= 7;
  }

  vector<pair<int8,int> > rbi_cb_pair_vec(ncb_i);
  for (int ii = 0; ii < ncb_i; ++ii)
    rbi_cb_pair_vec[ii] = pair<int8,int>(-1,ii);
  FOR_RANK {
    int isend = send_disp[rank];
    while (isend < send_disp[rank]+send_count[rank]) {
      const int icb_i = cb_send[isend]; assert((icb_i >= 0)&&(icb_i < ncb_i));
      const int icb_nbr = send_buf_int[isend];
      if (rbi_cb_pair_vec[icb_i].first == -1) 
        rbi_cb_pair_vec[icb_i].first = BitUtils::packRankBitsIndex(rank,0,icb_nbr);
      else 
        assert(rbi_cb_pair_vec[icb_i].first == BitUtils::packRankBitsIndex(rank,0,icb_nbr));
      ++isend;
    }
  }
  sort(rbi_cb_pair_vec.begin(),rbi_cb_pair_vec.end());
  assert(rbi_cb_i == NULL); rbi_cb_i = new int8[ncb_i];
  for (int ii = 0; ii < ncb_i; ++ii) {
    assert(rbi_cb_pair_vec[ii].first != -1);
    rbi_cb_i[ii] = rbi_cb_pair_vec[ii].first;
    assert((rbi_cb_pair_vec[ii].second >= 0)&&(rbi_cb_pair_vec[ii].second < ncb_i));
    const int icb0 = ncb_a+rbi_cb_pair_vec[ii].second;
    reorder_cb[icb0] = ncb_a+ii;
  }

  // reorder_cb should now have zone ordering for active cells and rbi ordering for inactive cells

  assert(bfocb_i == NULL); bfocb_i = new int[ncb+1]; // local bf's only
  for (int icb = 0; icb < ncb; ++icb) 
    bfocb_i[icb+1] = 0;
  for (int iter = 0; iter < 2; ++iter) {
    for (int ibf = 0; ibf < solver->nbf; ++ibf) {
      const int icc = ccocv[solver->cvobf[ibf]];
      if (icc < ncc) {
        map<const pair<int,int>,int>::iterator iter2; 
        int icb;
        if (icc < ncc_a) {
          iter2 = cbMap.find(pair<int,int>(solver->zone_bf[ibf],icc));
          assert(iter2 != cbMap.end());
          assert((iter2->second >= 0)&&(iter2->second < ncb_a));
          icb = reorder_cb[iter2->second]; assert((icb >= 0)&&(icb < ncb_a));
        }
        else {
          iter2 = cbMap_i.find(pair<int,int>(solver->zone_bf[ibf],icc));
          assert(iter2 != cbMap_i.end());
          assert((iter2->second >= 0)&&(iter2->second < ncb_i));
          icb = reorder_cb[iter2->second+ncb_a]; assert((icb >= ncb_a)&&(icb < ncb));
        }
        if (iter == 0) {
          ++bfocb_i[icb+1];
        }
        else {
          bfocb_v[bfocb_i[icb]++] = ibf;
          if (zone_cb[icb] == -1)
            zone_cb[icb] = solver->zone_bf[ibf];
          else
            assert(zone_cb[icb] == solver->zone_bf[ibf]);
          if (ccocb[icb] == -1)
            ccocb[icb] = icc;
          else
            assert(ccocb[icb] == icc);
          if (icc < ncc_a) {
            FOR_I3 x_cb[icb][i] += solver->x_bf[ibf][i]*solver->area_bf[ibf];
            FOR_I3 n_cb[icb][i] += solver->n_bf[ibf][i];
            area_cb[icb] += solver->area_bf[ibf];
          }
          else {
            // inactive coarse boundaries unset for now 
            FOR_I3 x_cb[icb][i] = HUGE_VAL;
            FOR_I3 n_cb[icb][i] = HUGE_VAL;
            area_cb[icb] = HUGE_VAL;
          }
        }
      }
    }

    if (iter == 0) {
      // index from counts
      bfocb_i[0] = 0;
      for (int icb = 0; icb < ncb; ++icb) 
        bfocb_i[icb+1] += bfocb_i[icb];
      assert(solver->nbf == bfocb_i[ncb]);
      assert(bfocb_v == NULL); bfocb_v = new int[bfocb_i[ncb]];
    }
    else {
      // rewind
      for (int icb = ncb; icb > 0; --icb)
        bfocb_i[icb] = bfocb_i[icb-1];
      bfocb_i[0] = 0;
    }
  }
  cbMap.clear();
  cbMap_i.clear();

  FOR_RANK {
    int isend = send_disp[rank];
    while (isend < send_disp[rank]+send_count[rank]) {
      const int icb = reorder_cb[cb_send[isend]+ncb_a]; assert((icb >= ncb_a)&&(icb < ncb));
      //const int icb_nbr = send_buf_int[isend];
      if (area_cb[icb] == HUGE_VAL) {
        FOR_I3 assert(x_cb[icb][i] == HUGE_VAL);
        FOR_I3 assert(n_cb[icb][i] == HUGE_VAL);
        FOR_I3 x_cb[icb][i] = send_buf_double[isend*7+i];
        FOR_I3 n_cb[icb][i] = send_buf_double[isend*7+3+i];
        area_cb[icb] = send_buf_double[isend*7+6];
      }
      else {
        FOR_I3 assert(x_cb[icb][i] == send_buf_double[isend*7+i]);
        FOR_I3 assert(n_cb[icb][i] == send_buf_double[isend*7+3+i]);
        assert(area_cb[icb] == send_buf_double[isend*7+6]);
      }
      ++isend;
    }
  }
  delete[] cb_send;
  delete[] send_buf_double; send_buf_double = NULL;
  delete[] send_buf_int; send_buf_int = NULL;

  // normalize x_cb and build delta_cb based on distance b/w centroids in normal direction

  assert(delta_cb == NULL); delta_cb = new double[ncb];
  for (int icb = 0; icb < ncb; ++icb) {
    assert(area_cb[icb] > 0.0);
    FOR_I3 x_cb[icb][i] /= area_cb[icb];
    const int icc = ccocb[icb]; assert((icc >= 0)&&(icc < ncc));
    delta_cb[icb] = DIST(x_cc[icc],x_cb[icb]); // total length (needed?)
    //const double proj_area = MAG(n_cb[icb]);
    //if (proj_area > 1.0E-20) {
    //  const double dx[3] = DIFF(x_cb[icb],x_cc[icc]);
    // is this correct?
    //  delta_cb[icb] = fabs(DOT_PRODUCT(dx,n_cb[icb])/proj_area); // projected length
    //}
    //else {
    //  delta_cb[icb] = DIST(x_cc[icc],x_cb[icb]); // total length (needed?)
    //}
  }
  dumpRange(delta_cb,ncb,"delta_cb");
  dumpRange(area_cb,ncb,"area_cb");
  dumpRange(n_cb,ncb,"n_cb");
  dumpRange(x_cb,ncb,"x_cb");

  // check geometric face information...

  double (*gcl_check)[3] = new double[ncc_a][3]; 
  for (int icc = 0; icc < ncc_a; ++icc) FOR_I3 gcl_check[icc][i] = 0.0;

  for (int icf = 0; icf < ncf_ai; ++icf) {
    const int icc0 = ccocf[icf][0]; assert((icc0 >= 0)&&(icc0 < ncc_a));
    const int icc1 = ccocf[icf][1]; assert((icc1 >= 0)&&(icc1 < ncc_g));
    FOR_I3 gcl_check[icc0][i] += n_cf[icf][i];
    if (icc1 < ncc_a)
      FOR_I3 gcl_check[icc1][i] -= n_cf[icf][i];
  }

  // add in boundary terms and reduce to active for checking...

  for (int icb = 0; icb < ncb; ++icb) {
    const int icc = ccocb[icb]; assert((icc >= 0)&&(icc < ncc));
    if (icc < ncc_a) 
      FOR_I3 gcl_check[icc][i] += n_cb[icb][i];
  }
  dumpRange(gcl_check,ncc_a,"gcl");
  delete[] gcl_check;

  // should be zone ordered for active cb's, so build disps
  assert(cbozn_i == NULL); cbozn_i = new int[solver->bfZoneVec.size()+1];
  FOR_IZONE(solver->bfZoneVec) 
    cbozn_i[izone+1] = 0;
  for (int icb = 0; icb < ncb_a; ++icb) 
    ++cbozn_i[zone_cb[icb]+1];
  cbozn_i[0] = 0;
  FOR_IZONE(solver->bfZoneVec) 
    cbozn_i[izone+1] += cbozn_i[izone];
  assert(cbozn_i[solver->bfZoneVec.size()] == ncb_a);
  delete[] reorder_cb;

  buildCbIPrcomm();

  double* area_cb_check = new double[ncb];
  for (int icb = 0; icb < ncb; ++icb) {
    area_cb_check[icb] = 0.0;
    for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
      const int ibf = bfocb_v[bob];
      area_cb_check[icb] += solver->area_bf[ibf];
    }
  }
  updateCbIData(area_cb_check);
  for (int icb = 0; icb < ncb_a; ++icb) 
    area_cb_check[icb] -= area_cb[icb];
  dumpRange(area_cb_check,ncb_a,"area_cb check");
  delete[] area_cb_check;

  double (*x_cb_check)[3] = new double[ncb][3];
  restrictCbData(x_cb_check,solver->x_bf);
  for (int icb = 0; icb < ncb_a; ++icb) 
    FOR_I3 x_cb_check[icb][i] -= x_cb[icb][i];
  dumpRange(x_cb_check,ncb_a,"x_cb check");
  delete[] x_cb_check;

  // compare sum cb data to zone data

  for (int izone = 0; izone < solver->bfZoneVec.size(); ++izone) {
    double my_buf[7] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0}; 
    for (int icb = cbozn_i[izone]; icb < cbozn_i[izone+1]; ++icb) {
      my_buf[0] += area_cb[icb];
      FOR_I3 my_buf[1+i] += n_cb[icb][i];
      FOR_I3 my_buf[4+i] += area_cb[icb]*x_cb[icb][i];
    }
    double buf[7];
    MPI_Reduce(my_buf,buf,7,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    if (mpi_rank == 0) {
      FOR_I3 buf[4+i] /= buf[0];
      buf[0] -= solver->bfZoneVec[izone].area_global;
      FOR_I3 buf[1+i] -= solver->bfZoneVec[izone].n_global[i];
      FOR_I3 buf[4+i] -= solver->bfZoneVec[izone].x_global[i];
      cout << " > zone: " << izone << ", area_global check: " << buf[0] << ", n_global check: " << 
        COUT_VEC(&buf[1]) << ", x_global check: " << COUT_VEC(&buf[4]) << endl;
    }
  }

  delete[] send_count;
  delete[] send_disp;
  delete[] recv_count;
  delete[] recv_disp;

}

void StaticSolver::CoarseGrid::buildCcIPrcomm() {

  COUT1("CoarseGrid::buildCcIPrcomm()");

  // ------------------------------------------------
  // build the inactive to active communicator
  // ------------------------------------------------

  assert(ccIPrcommVec.empty());

  // build pack side using rbi_i...
  {
    int rank_current = -1;
    int bits_current = -1;
    Prcomm * prcomm = NULL;
    for (int icc = ncc_a; icc < ncc; ++icc) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_i[icc-ncc_a]);
      assert(bits == 0);
      if (rank > rank_current) {
        // we just switched ranks. If there was a prcomm, then
        // complete the size...
        if (prcomm) {
          prcomm->pack_size = icc-prcomm->pack_offset;
          prcomm->packVec.resize(prcomm->pack_size);
          for (int icc2 = prcomm->pack_offset; icc2 < icc; ++icc2) {
            prcomm->packVec[icc2-prcomm->pack_offset] = icc2;
          }
        }
        rank_current = rank;
        bits_current = -1;
        ccIPrcommVec.push_back(Prcomm());
        prcomm = &ccIPrcommVec.back();
        prcomm->rank = rank;
        prcomm->pack_offset = icc;
      }
      else {
        assert(rank_current == rank);
        assert(prcomm);
      }
      // just tests monotoncity of the bits...
      if (bits > bits_current) {
        bits_current = bits;
      }
      else {
        assert(bits_current == bits);
      }
    }
    // we just finished. If there was a prcomm, then
    // complete the size...
    if (prcomm) {
      prcomm->pack_size = ncc-prcomm->pack_offset;
      prcomm->packVec.resize(prcomm->pack_size);
      for (int icc2 = prcomm->pack_offset; icc2 < ncc; ++icc2) {
        prcomm->packVec[icc2-prcomm->pack_offset] = icc2;
      }
    }
  }

  // we need to get the rbi's for the inactive cells that go with our split active cells...

  // get rbi's for the inactive cells that contribute to our active cells 
  int* send_count = new int[mpi_size];
  int* send_disp = new int[mpi_size];
  int* send_buf_int = new int[2*(ncc-ncc_a)]; // nbr index,index
  FOR_RANK send_count[rank] = 0;
  for (int iter = 0; iter < 2; ++iter) {
    for (int icc = ncc_a; icc < ncc; ++icc) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_i[icc-ncc_a]); 
      assert(rank != mpi_rank);
      assert(bits == 0);
      if (iter == 0) {
        send_count[rank] += 2;
      }
      else {
        send_buf_int[send_disp[rank]++] = index;
        send_buf_int[send_disp[rank]++] = icc;
      }
    }
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    assert((send_count[mpi_size-1] + send_disp[mpi_size-1]) == 2*(ncc-ncc_a));
  }

  // setup recv side stuff

  int* recv_count = new int[mpi_size];
  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
  assert(recv_count[mpi_rank] == 0); 

  int* recv_disp = new int[mpi_size];
  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

  int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
  int* recv_buf_int = new int[recv_count_sum];
  MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
      recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
  delete[] send_buf_int; 
  delete[] send_count;
  delete[] send_disp;

  assert(recv_count_sum%2 == 0);
  vector<pair<int8,int> > rbi_index_pair_vec(recv_count_sum/2);
  FOR_RANK {
    int irecv = recv_disp[rank];
    assert(irecv%2 == 0);
    while (irecv < recv_disp[rank]+recv_count[rank]) {
      assert(recv_buf_int[irecv] < ncc_a); // touches my active cell
      rbi_index_pair_vec[irecv/2] = pair<int8,int>(BitUtils::packRankBitsIndex(rank,0,recv_buf_int[irecv+1]),recv_buf_int[irecv]);
      irecv += 2;
    }
  }
  delete[] recv_buf_int;
  delete[] recv_count;
  delete[] recv_disp;
  sort(rbi_index_pair_vec.begin(),rbi_index_pair_vec.end());

  {
    int rank_current = -1;
    int bits_current = -1;
    Prcomm * prcomm = NULL;
    for (int ii = 0, ii_end=rbi_index_pair_vec.size(); ii < ii_end; ++ii) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_index_pair_vec[ii].first);
      assert(rank != mpi_rank);
      assert(bits == 0);
      if (rank > rank_current) {
        // we just switched ranks. If there was a prcomm, then
        // complete the size...
        if (prcomm) {
          prcomm->unpack_size = ii-prcomm->unpack_offset;
          prcomm->unpackVec.resize(prcomm->unpack_size);
          for (int ii2 = prcomm->unpack_offset; ii2 < ii; ++ii2) {
            prcomm->unpackVec[ii2-prcomm->unpack_offset] = rbi_index_pair_vec[ii2].second;
          }
        }
        rank_current = rank;
        bits_current = -1;
        vector<Prcomm>::iterator it;
        for (it = ccIPrcommVec.begin(); it != ccIPrcommVec.end(); ++it) {
          if (it->rank == rank) {
            prcomm = &(*it);
            break;
          }
        }
        if (it == ccIPrcommVec.end()) {
          ccIPrcommVec.push_back(Prcomm());
          prcomm = &ccIPrcommVec.back();
          prcomm->rank = rank;
        }
        prcomm->unpack_offset = ii;
      }
      else {
        assert(rank_current == rank);
        assert(prcomm);
      }
      if (bits > bits_current) {
        bits_current = bits;
      }
      else {
        assert(bits_current == bits);
      }
    }
    // we just finished. If there was a prcomm, then
    // complete the size...
    if (prcomm) {
      prcomm->unpack_size = rbi_index_pair_vec.size()-prcomm->unpack_offset;
      prcomm->unpackVec.resize(prcomm->unpack_size);
      for (int ii2 = prcomm->unpack_offset, ii2_end=rbi_index_pair_vec.size(); ii2 < ii2_end; ++ii2) {
        prcomm->unpackVec[ii2-prcomm->unpack_offset] = rbi_index_pair_vec[ii2].second;
      }
    }
  }

}

void StaticSolver::CoarseGrid::buildCcPrcomm() {

  COUT1("CoarseGrid::buildCcPrcomm()");

  // ------------------------------------------------
  // build the active to inactive/ghost communicator
  // ------------------------------------------------

  assert(ccPrcommVec.empty());

  // we need to get the rbi's for the ghost cells that go with our active cells...

  int* send_count = new int[mpi_size];
  int* send_disp = new int[mpi_size];
  int* send_buf_int = new int[2*(ncc_g-ncc)]; // nbr index,index
  FOR_RANK send_count[rank] = 0;
  for (int iter = 0; iter < 2; ++iter) {
    for (int icc = ncc; icc < ncc_g; ++icc) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icc-ncc]); 
      assert(rank != mpi_rank); 
      assert(bits == 0);
      if (iter == 0) {
        send_count[rank] += 2;
      }
      else {
        send_buf_int[send_disp[rank]++] = index;
        send_buf_int[send_disp[rank]++] = icc;
      }
    }
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    assert((send_count[mpi_size-1] + send_disp[mpi_size-1]) == 2*(ncc_g-ncc));
  }

  // setup recv side stuff

  int* recv_count = new int[mpi_size];
  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
  assert(recv_count[mpi_rank] == 0); 

  int* recv_disp = new int[mpi_size];
  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

  int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
  int* recv_buf_int = new int[recv_count_sum];
  MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
      recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
  delete[] send_buf_int; 
  delete[] send_count;
  delete[] send_disp;

  assert(recv_count_sum%2 == 0);
  vector<pair<int8,int> > rbi_index_pair_vec(recv_count_sum/2);
  set<int> auxRankSet;
  FOR_RANK {
    int irecv = recv_disp[rank];
    assert(irecv%2 == 0);
    while (irecv < recv_disp[rank]+recv_count[rank]) {
      assert(recv_buf_int[irecv] < ncc_a); // touches my active cell
      rbi_index_pair_vec[irecv/2] = pair<int8,int>(BitUtils::packRankBitsIndex(rank,0,recv_buf_int[irecv+1]),recv_buf_int[irecv]);
      auxRankSet.insert(rank);
      irecv += 2;
    }
  }
  delete[] recv_buf_int;
  delete[] recv_count;
  delete[] recv_disp;
  sort(rbi_index_pair_vec.begin(),rbi_index_pair_vec.end());

  // throw ccIPrcommVec ranks into set 
  for (int ii = 0, lim = ccIPrcommVec.size(); ii < lim; ++ii)
    auxRankSet.insert(ccIPrcommVec[ii].getRank());

  // build unpack side using rbi_g...
  {
    set<int>::iterator iter = auxRankSet.begin();
    int rank_current = -1;
    int bits_current = -1;
    Prcomm * prcomm = NULL;
    for (int icc = ncc; icc < ncc_g; ++icc) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icc-ncc]);
      assert(bits == 0);
      if (rank > rank_current) {
        // we just switched ranks. If there was a prcomm, then
        // complete the size...
        if (prcomm) 
          prcomm->unpack_size = icc-prcomm->unpack_offset;
        while ((iter != auxRankSet.end())&&((*iter) <= rank)) {
          if ((*iter) < rank) {
            // add an empty Prcomm...
            ccPrcommVec.push_back(Prcomm());
            ccPrcommVec.back().rank = *iter;
            ccPrcommVec.back().unpack_offset = icc;
            ccPrcommVec.back().unpack_size = 0;
          }
          ++iter;
        }
        rank_current = rank;
        bits_current = -1;
        ccPrcommVec.push_back(Prcomm());
        prcomm = &ccPrcommVec.back();
        prcomm->rank = rank;
        prcomm->unpack_offset = icc;
      }
      else {
        assert(rank_current == rank);
        assert(prcomm);
      }
      // just tests monotoncity of the bits...
      if (bits > bits_current) {
        bits_current = bits;
      }
      else {
        assert(bits_current == bits);
      }
    }
    // we just finished. If there was a prcomm, then
    // complete the size...
    if (prcomm) 
      prcomm->unpack_size = ncc_g-prcomm->unpack_offset;
    // add any not represented in Prcomm at the end...
    while (iter != auxRankSet.end()) {
      // add an empty Prcomm...
      ccPrcommVec.push_back(Prcomm());
      ccPrcommVec.back().rank = *iter;
      ccPrcommVec.back().unpack_offset = ncc_g;
      ccPrcommVec.back().unpack_size = 0;
      ++iter;
    }
    assert(iter == auxRankSet.end());
  }
  auxRankSet.clear();

  {
    int rank_current = -1;
    int bits_current = -1;
    Prcomm * prcomm = NULL;
    for (int ii = 0, ii_end=rbi_index_pair_vec.size(); ii < ii_end; ++ii) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_index_pair_vec[ii].first);
      assert(rank != mpi_rank); 
      assert(bits == 0);
      if (rank > rank_current) {
        // we just switched ranks. If there was a prcomm, then
        // complete the size...
        if (prcomm) {
          prcomm->pack_size = ii-prcomm->pack_offset;
          prcomm->packVec.resize(prcomm->pack_size);
          for (int ii2 = prcomm->pack_offset; ii2 < ii; ++ii2) {
            prcomm->packVec[ii2-prcomm->pack_offset] = rbi_index_pair_vec[ii2].second;
          }
        }
        rank_current = rank;
        bits_current = -1;
        vector<Prcomm>::iterator it;
        for (it = ccPrcommVec.begin(); it != ccPrcommVec.end(); ++it) {
          if (it->rank == rank) {
            prcomm = &(*it);
            break;
          }
        }
        assert(it != ccPrcommVec.end()); // added above so should NOT add any new ones
        prcomm->pack_offset = ii;
      }
      else {
        assert(rank_current == rank);
        assert(prcomm);
      }
      if (bits > bits_current) {
        bits_current = bits;
      }
      else {
        assert(bits_current == bits);
      }
    }
    // we just finished. If there was a prcomm, then
    // complete the size...
    if (prcomm) {
      prcomm->pack_size = rbi_index_pair_vec.size()-prcomm->pack_offset;
      prcomm->packVec.resize(prcomm->pack_size);
      for (int ii2 = prcomm->pack_offset, ii2_end=rbi_index_pair_vec.size(); ii2 < ii2_end; ++ii2) {
        prcomm->packVec[ii2-prcomm->pack_offset] = rbi_index_pair_vec[ii2].second;
      }
    }
  }

}

void StaticSolver::CoarseGrid::buildCbIPrcomm() {

  COUT1("CoarseGrid::buildCbIPrcomm()");

  // ------------------------------------------------
  // build the inactive to active communicator
  // ------------------------------------------------

  assert(cbIPrcommVec.empty());

  // build pack side using rbi_cb_i...
  {
    int rank_current = -1;
    int bits_current = -1;
    Prcomm * prcomm = NULL;
    for (int icb = ncb_a; icb < ncb; ++icb) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_cb_i[icb-ncb_a]);
      assert(bits == 0);
      if (rank > rank_current) {
        // we just switched ranks. If there was a prcomm, then
        // complete the size...
        if (prcomm) {
          prcomm->pack_size = icb-prcomm->pack_offset;
          prcomm->packVec.resize(prcomm->pack_size);
          for (int icb2 = prcomm->pack_offset; icb2 < icb; ++icb2) {
            prcomm->packVec[icb2-prcomm->pack_offset] = icb2;
          }
        }
        rank_current = rank;
        bits_current = -1;
        cbIPrcommVec.push_back(Prcomm());
        prcomm = &cbIPrcommVec.back();
        prcomm->rank = rank;
        prcomm->pack_offset = icb;
      }
      else {
        assert(rank_current == rank);
        assert(prcomm);
      }
      // just tests monotoncity of the bits...
      if (bits > bits_current) {
        bits_current = bits;
      }
      else {
        assert(bits_current == bits);
      }
    }
    // we just finished. If there was a prcomm, then
    // complete the size...
    if (prcomm) {
      prcomm->pack_size = ncb-prcomm->pack_offset;
      prcomm->packVec.resize(prcomm->pack_size);
      for (int icb2 = prcomm->pack_offset; icb2 < ncb; ++icb2) {
        prcomm->packVec[icb2-prcomm->pack_offset] = icb2;
      }
    }
  }

  // we need to get the rbi's for the inactive cells that go with our split active cells...

  // get rbi's for the inactive cells that contribute to our active cells 
  int* send_count = new int[mpi_size];
  int* send_disp = new int[mpi_size];
  int* send_buf_int = new int[2*(ncb-ncb_a)]; // nbr index,index
  FOR_RANK send_count[rank] = 0;
  for (int iter = 0; iter < 2; ++iter) {
    for (int icb = ncb_a; icb < ncb; ++icb) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_cb_i[icb-ncb_a]); 
      assert(rank != mpi_rank);
      assert(bits == 0);
      if (iter == 0) {
        send_count[rank] += 2;
      }
      else {
        send_buf_int[send_disp[rank]++] = index;
        send_buf_int[send_disp[rank]++] = icb;
      }
    }
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    assert((send_count[mpi_size-1] + send_disp[mpi_size-1]) == 2*(ncb-ncb_a));
  }

  // setup recv side stuff

  int* recv_count = new int[mpi_size];
  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
  assert(recv_count[mpi_rank] == 0); 

  int* recv_disp = new int[mpi_size];
  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

  int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
  int* recv_buf_int = new int[recv_count_sum];
  MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
      recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
  delete[] send_buf_int; 
  delete[] send_count;
  delete[] send_disp;

  assert(recv_count_sum%2 == 0);
  vector<pair<int8,int> > rbi_index_pair_vec(recv_count_sum/2);
  FOR_RANK {
    int irecv = recv_disp[rank];
    assert(irecv%2 == 0);
    while (irecv < recv_disp[rank]+recv_count[rank]) {
      assert(recv_buf_int[irecv] < ncb_a); // touches my active cell
      rbi_index_pair_vec[irecv/2] = pair<int8,int>(BitUtils::packRankBitsIndex(rank,0,recv_buf_int[irecv+1]),recv_buf_int[irecv]);
      irecv += 2;
    }
  }
  delete[] recv_buf_int;
  delete[] recv_count;
  delete[] recv_disp;
  sort(rbi_index_pair_vec.begin(),rbi_index_pair_vec.end());

  {
    int rank_current = -1;
    int bits_current = -1;
    Prcomm * prcomm = NULL;
    for (int ii = 0, ii_end=rbi_index_pair_vec.size(); ii < ii_end; ++ii) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_index_pair_vec[ii].first);
      assert(rank != mpi_rank);
      assert(bits == 0);
      if (rank > rank_current) {
        // we just switched ranks. If there was a prcomm, then
        // complete the size...
        if (prcomm) {
          prcomm->unpack_size = ii-prcomm->unpack_offset;
          prcomm->unpackVec.resize(prcomm->unpack_size);
          for (int ii2 = prcomm->unpack_offset; ii2 < ii; ++ii2) {
            prcomm->unpackVec[ii2-prcomm->unpack_offset] = rbi_index_pair_vec[ii2].second;
          }
        }
        rank_current = rank;
        bits_current = -1;
        vector<Prcomm>::iterator it;
        for (it = cbIPrcommVec.begin(); it != cbIPrcommVec.end(); ++it) {
          if (it->rank == rank) {
            prcomm = &(*it);
            break;
          }
        }
        if (it == cbIPrcommVec.end()) {
          cbIPrcommVec.push_back(Prcomm());
          prcomm = &cbIPrcommVec.back();
          prcomm->rank = rank;
        }
        prcomm->unpack_offset = ii;
      }
      else {
        assert(rank_current == rank);
        assert(prcomm);
      }
      if (bits > bits_current) {
        bits_current = bits;
      }
      else {
        assert(bits_current == bits);
      }
    }
    // we just finished. If there was a prcomm, then
    // complete the size...
    if (prcomm) {
      prcomm->unpack_size = rbi_index_pair_vec.size()-prcomm->unpack_offset;
      prcomm->unpackVec.resize(prcomm->unpack_size);
      for (int ii2 = prcomm->unpack_offset, ii2_end=rbi_index_pair_vec.size(); ii2 < ii2_end; ++ii2) {
        prcomm->unpackVec[ii2-prcomm->unpack_offset] = rbi_index_pair_vec[ii2].second;
      }
    }
  }

}

void StaticSolver::CoarseGrid::buildCfIPrcomm() {

  COUT1("CoarseGrid::buildCfIPrcomm()");

  // ------------------------------------------------
  // build the inactive to active communicator
  // ------------------------------------------------

  assert(cfIPrcommVec.empty());

  // need to locally sort the boundary faces to building the pack
  // shift index over 1 bit to store sign without impacting sort
  vector<pair<int8,int> > rbi_index_pair_vec;
  for (int icf = ncf_aa; icf < ncf; ++icf) {
    for (int ii = 0, lim = rbiVec_cf_i[icf-ncf_aa].size(); ii < lim; ++ii) {
      if (rbiVec_cf_i[icf-ncf_aa][ii] >= 0)
        rbi_index_pair_vec.push_back(pair<int8,int>(rbiVec_cf_i[icf-ncf_aa][ii],icf<<1));
      else
        rbi_index_pair_vec.push_back(pair<int8,int>(-rbiVec_cf_i[icf-ncf_aa][ii]-1,(icf<<1)|1));
    }
  }
  sort(rbi_index_pair_vec.begin(),rbi_index_pair_vec.end());

  // build pack side using rbi_index_pair_vec...
  {
    int rank_current = -1;
    int bits_current = -1;
    Prcomm * prcomm = NULL;
    for (int ii = 0, ii_end=rbi_index_pair_vec.size(); ii < ii_end; ++ii) {
      int rank,bits,index;
      assert(rbi_index_pair_vec[ii].first >= 0);
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_index_pair_vec[ii].first);
      assert(rank != mpi_rank);
      assert(bits == 0);
      if (rank > rank_current) {
        // we just switched ranks. If there was a prcomm, then
        // complete the size...
        if (prcomm) {
          prcomm->pack_size = ii-prcomm->pack_offset;
          prcomm->packVec.resize(prcomm->pack_size);
          for (int ii2 = prcomm->pack_offset; ii2 < ii; ++ii2) {
            const int icf  = (rbi_index_pair_vec[ii2].second>>1); assert((icf >= ncf_aa)&&(icf < ncf));
            const int flip = (rbi_index_pair_vec[ii2].second&1); assert((flip == 0)||(flip == 1));
            if (flip == 0)
              prcomm->packVec[ii2-prcomm->pack_offset] = icf; 
            else 
              prcomm->packVec[ii2-prcomm->pack_offset] = -icf-1; // pack-side vec tells you if you need to flip
            rbi_index_pair_vec[ii2].second = ii2-prcomm->pack_offset; // to sort rbi on other rank
          }
        }
        rank_current = rank;
        bits_current = -1;
        cfIPrcommVec.push_back(Prcomm());
        prcomm = &cfIPrcommVec.back();
        prcomm->rank = rank;
        prcomm->pack_offset = ii;
      }
      else {
        assert(rank_current == rank);
        assert(prcomm);
      }
      // just tests monotoncity of the bits...
      if (bits > bits_current) {
        bits_current = bits;
      }
      else {
        assert(bits_current == bits);
      }
    }
    // we just finished. If there was a prcomm, then
    // complete the size...
    if (prcomm) {
      prcomm->pack_size = rbi_index_pair_vec.size()-prcomm->pack_offset;
      prcomm->packVec.resize(prcomm->pack_size);
      for (int ii2 = prcomm->pack_offset, ii2_end=rbi_index_pair_vec.size(); ii2 < ii2_end; ++ii2) {
        const int icf  = (rbi_index_pair_vec[ii2].second>>1); assert((icf >= ncf_aa)&&(icf < ncf));
        const int flip = (rbi_index_pair_vec[ii2].second&1); assert((flip == 0)||(flip == 1));
        if (flip == 0)
          prcomm->packVec[ii2-prcomm->pack_offset] = icf;
        else 
          prcomm->packVec[ii2-prcomm->pack_offset] = -icf-1; // pack-side vec tells you if you need to flip
        rbi_index_pair_vec[ii2].second = ii2-prcomm->pack_offset;
      }
    }
  }

  // we need to get the rbi's for the inactive cells that go with our split active cells...

  // get rbi's for the inactive cells that contribute to our active cells 
  int* send_count = new int[mpi_size];
  int* send_disp = new int[mpi_size];
  int* send_buf_int = new int[2*int(rbi_index_pair_vec.size())]; // nbr index,index
  FOR_RANK send_count[rank] = 0;
  for (int iter = 0; iter < 2; ++iter) {
    for (int ii = 0, ii_end=rbi_index_pair_vec.size(); ii < ii_end; ++ii) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_index_pair_vec[ii].first);
      assert(rank != mpi_rank);
      assert(bits == 0);
      if (iter == 0) {
        send_count[rank] += 2;
      }
      else {
        send_buf_int[send_disp[rank]++] = index;
        send_buf_int[send_disp[rank]++] = rbi_index_pair_vec[ii].second; // rel index into pack
      }
    }
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    assert((send_count[mpi_size-1] + send_disp[mpi_size-1]) == 2*int(rbi_index_pair_vec.size()));
  }
  rbi_index_pair_vec.clear();

  // setup recv side stuff

  int* recv_count = new int[mpi_size];
  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
  assert(recv_count[mpi_rank] == 0); 

  int* recv_disp = new int[mpi_size];
  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

  int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
  int* recv_buf_int = new int[recv_count_sum];
  MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
      recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
  delete[] send_buf_int; 
  delete[] send_count;
  delete[] send_disp;

  assert(recv_count_sum%2 == 0);
  rbi_index_pair_vec.resize(recv_count_sum/2);
  FOR_RANK {
    int irecv = recv_disp[rank];
    assert(irecv%2 == 0);
    while (irecv < recv_disp[rank]+recv_count[rank]) {
      const int icf = recv_buf_int[irecv]; assert((icf >= 0)&&(icf < ncf_ai));
      rbi_index_pair_vec[irecv/2] = pair<int8,int>(BitUtils::packRankBitsIndex(rank,0,recv_buf_int[irecv+1]),icf);
      irecv += 2;
    }
  }
  delete[] recv_buf_int;
  delete[] recv_count;
  delete[] recv_disp;
  sort(rbi_index_pair_vec.begin(),rbi_index_pair_vec.end());

  {
    int rank_current = -1;
    int bits_current = -1;
    Prcomm * prcomm = NULL;
    for (int ii = 0, ii_end=rbi_index_pair_vec.size(); ii < ii_end; ++ii) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_index_pair_vec[ii].first);
      assert(rank != mpi_rank);
      assert(bits == 0);
      if (rank > rank_current) {
        // we just switched ranks. If there was a prcomm, then
        // complete the size...
        if (prcomm) {
          prcomm->unpack_size = ii-prcomm->unpack_offset;
          prcomm->unpackVec.resize(prcomm->unpack_size);
          for (int ii2 = prcomm->unpack_offset; ii2 < ii; ++ii2) {
            prcomm->unpackVec[ii2-prcomm->unpack_offset] = rbi_index_pair_vec[ii2].second;
          }
        }
        rank_current = rank;
        bits_current = -1;
        vector<Prcomm>::iterator it;
        for (it = cfIPrcommVec.begin(); it != cfIPrcommVec.end(); ++it) {
          if (it->rank == rank) {
            prcomm = &(*it);
            break;
          }
        }
        if (it == cfIPrcommVec.end()) {
          cfIPrcommVec.push_back(Prcomm());
          prcomm = &cfIPrcommVec.back();
          prcomm->rank = rank;
        }
        prcomm->unpack_offset = ii;
      }
      else {
        assert(rank_current == rank);
        assert(prcomm);
      }
      if (bits > bits_current) {
        bits_current = bits;
      }
      else {
        assert(bits_current == bits);
      }
    }
    // we just finished. If there was a prcomm, then
    // complete the size...
    if (prcomm) {
      prcomm->unpack_size = rbi_index_pair_vec.size()-prcomm->unpack_offset;
      prcomm->unpackVec.resize(prcomm->unpack_size);
      for (int ii2 = prcomm->unpack_offset, ii2_end=rbi_index_pair_vec.size(); ii2 < ii2_end; ++ii2) {
        prcomm->unpackVec[ii2-prcomm->unpack_offset] = rbi_index_pair_vec[ii2].second;
      }
    }
  }

}

void StaticSolver::CoarseGrid::prolongCcData(double* phi_cv,const double* phi_cc) {

  // make sure that you call updateCcData to update phi_cc prior to calling this routine to get it in inactive cells. 

  for (int icc = 0; icc < ncc; ++icc) {
    for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
      const int icv = cvocc_v[coc];
      phi_cv[icv] = phi_cc[icc];
    }
  }
  //dumpRange(phi_cc,ncc_a,"phi_cc");
  //dumpRange(phi_cv,solver->ncv,"phi_cv");

}

void StaticSolver::CoarseGrid::prolongCcDataAndUpdateGuess(double* phi_cv,const double* phi_cc) {

  // make sure that you call updateCcData to update phi_cc prior to calling this routine to get it in inactive cells. 

  for (int icc = 0; icc < ncc; ++icc) {
    for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
      const int icv = cvocc_v[coc];
      phi_cv[icv] += phi_cc[icc];
    }
  }
  //dumpRange(phi_cc,ncc_a,"phi_cc");
  //dumpRange(phi_cv,solver->ncv,"phi_cv");

}


void StaticSolver::CoarseGrid::prolongCcData(double (*u_cv)[3],const double (*u_cc)[3]) {

  // make sure that you call updateCcData to update u_cc prior to calling this routine to get it in inactive cells. 

  for (int icc = 0; icc < ncc; ++icc) {
    for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
      const int icv = cvocc_v[coc];
      FOR_I3 u_cv[icv][i] = u_cc[icc][i];
    }
  }

}

void StaticSolver::CoarseGrid::prolongCcDataAndUpdateGuess(double (*u_cv)[3],const double (*u_cc)[3],const double relax) {

  // make sure that you call updateCcData to update u_cc prior to calling this routine to get it in inactive cells. 

  for (int icc = 0; icc < ncc; ++icc) {
    for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
      const int icv = cvocc_v[coc];
      FOR_I3 u_cv[icv][i] += relax*u_cc[icc][i];
    }
  }

}

void StaticSolver::CoarseGrid::restrictCcData(double* phi_cc,const double* phi_cv) {

  // restrict inactive cells first and send to active

  for (int icc = ncc_a; icc < ncc; ++icc) {
    phi_cc[icc] = 0.0;
    for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
      const int icv = cvocc_v[coc];
      phi_cc[icc] += solver->vol_cv[icv]*phi_cv[icv];
    }
  }
  updateCcIDataStart(phi_cc); 

  // restrict the active cells

  for (int icc = 0; icc < ncc_a; ++icc) {
    phi_cc[icc] = 0.0;
    for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
      const int icv = cvocc_v[coc];
      phi_cc[icc] += solver->vol_cv[icv]*phi_cv[icv];
    }
  }
  updateCcIDataFinish(phi_cc);

  // normalize

  for (int icc = 0; icc < ncc_a; ++icc)
    phi_cc[icc] /= vol_cc[icc];

  // you probably will need to call updateCcData to have phi_cc in inactive/ghosts after calling this routine 

  //dumpRange(phi_cv,solver->ncv,"phi_cv");
  //dumpRange(phi_cc,ncc_a,"phi_cc");
}

void StaticSolver::CoarseGrid::restrictExtrinsicCcData(double* phiV_cc,const double* phiV_cv) {

  // restrict inactive cells first and send to active

  for (int icc = ncc_a; icc < ncc; ++icc) {
    phiV_cc[icc] = 0.0;
    for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
      const int icv = cvocc_v[coc];
      phiV_cc[icc] += phiV_cv[icv];
    }
  }
  updateCcIDataStart(phiV_cc); 

  // restrict the active cells

  for (int icc = 0; icc < ncc_a; ++icc) {
    phiV_cc[icc] = 0.0;
    for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
      const int icv = cvocc_v[coc];
      phiV_cc[icc] += phiV_cv[icv];
    }
  }
  updateCcIDataFinish(phiV_cc);

  // you probably will need to call updateCcData to have phiV_cc in inactive/ghosts after calling this routine 

}

void StaticSolver::CoarseGrid::restrictCcData(double (*u_cc)[3],const double (*u_cv)[3]) {

  for (int icc = ncc_a; icc < ncc; ++icc) {
    FOR_I3 u_cc[icc][i] = 0.0;
    for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
      const int icv = cvocc_v[coc];
      FOR_I3 u_cc[icc][i] += solver->vol_cv[icv]*u_cv[icv][i];
    }
  }
  updateCcIDataStart(u_cc); 

  for (int icc = 0; icc < ncc_a; ++icc) {
    FOR_I3 u_cc[icc][i] = 0.0;
    for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
      const int icv = cvocc_v[coc];
      FOR_I3 u_cc[icc][i] += solver->vol_cv[icv]*u_cv[icv][i];
    }
  }
  updateCcIDataFinish(u_cc); 

  for (int icc = 0; icc < ncc_a; ++icc)
    FOR_I3 u_cc[icc][i] /= vol_cc[icc];

  // you probably will need to call updateCcData to have phi_cc in inactive/ghosts after calling this routine 

}

void StaticSolver::CoarseGrid::restrictExtrinsicCcData(double (*uV_cc)[3],const double (*uV_cv)[3]) {

  for (int icc = ncc_a; icc < ncc; ++icc) {
    FOR_I3 uV_cc[icc][i] = 0.0;
    for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
      const int icv = cvocc_v[coc];
      FOR_I3 uV_cc[icc][i] += uV_cv[icv][i];
    }
  }
  updateCcIDataStart(uV_cc); 

  for (int icc = 0; icc < ncc_a; ++icc) {
    FOR_I3 uV_cc[icc][i] = 0.0;
    for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
      const int icv = cvocc_v[coc];
      FOR_I3 uV_cc[icc][i] += uV_cv[icv][i];
    }
  }
  updateCcIDataFinish(uV_cc); 

  // you probably will need to call updateCcData to have phi_cc in inactive/ghosts after calling this routine 

}


void StaticSolver::CoarseGrid::restrictCfData(double* phi_cf,const double* phi_fa) {

  // populate faces w/ at least one inactive cell and send to active rank

  for (int icf = ncf_aa; icf < ncf; ++icf) {
    phi_cf[icf] = 0.0;
    for (int fof = faocf_i[icf]; fof < faocf_i[icf+1]; ++fof) {
      const int ifa = max(faocf_v[fof],-faocf_v[fof]-1);
      phi_cf[icf] += solver->area_fa[ifa]*phi_fa[ifa];
    }
  }
  updateCfIDataStart(phi_cf); 

  // restrict internal faces

  for (int icf = 0; icf < ncf_aa; ++icf) {
    phi_cf[icf] = 0.0;
    for (int fof = faocf_i[icf]; fof < faocf_i[icf+1]; ++fof) {
      const int ifa = max(faocf_v[fof],-faocf_v[fof]-1);
      phi_cf[icf] += solver->area_fa[ifa]*phi_fa[ifa];
    }
  }
  updateCfIDataFinish(phi_cf); 

  // normalize

  for (int icf = 0; icf < ncf_ai; ++icf)
    phi_cf[icf] /= area_cf[icf];

}

void StaticSolver::CoarseGrid::restrictCfData(double (*u_cf)[3],const double (*u_fa)[3]) {

  for (int icf = ncf_aa; icf < ncf; ++icf) {
    FOR_I3 u_cf[icf][i] = 0.0;
    for (int fof = faocf_i[icf]; fof < faocf_i[icf+1]; ++fof) {
      const int ifa = max(faocf_v[fof],-faocf_v[fof]-1);
      FOR_I3 u_cf[icf][i] += solver->area_fa[ifa]*u_fa[ifa][i];
    }
  }
  updateCfIDataStart(u_cf); 

  for (int icf = 0; icf < ncf_aa; ++icf) {
    FOR_I3 u_cf[icf][i] = 0.0;
    for (int fof = faocf_i[icf]; fof < faocf_i[icf+1]; ++fof) {
      const int ifa = max(faocf_v[fof],-faocf_v[fof]-1);
      FOR_I3 u_cf[icf][i] += solver->area_fa[ifa]*u_fa[ifa][i];
    }
  }
  updateCfIDataFinish(u_cf); 

  for (int icf = 0; icf < ncf_ai; ++icf)
    FOR_I3 u_cf[icf][i] /= area_cf[icf];

}

void StaticSolver::CoarseGrid::restrictSignedCfData(double* phi_cf,const double* phi_fa) {

  for (int icf = ncf_aa; icf < ncf; ++icf) {
    phi_cf[icf] = 0.0;
    for (int fof = faocf_i[icf]; fof < faocf_i[icf+1]; ++fof) {
      int ifa = faocf_v[fof];
      double fa_sign = 1.0;
      if (ifa < 0) {
        ifa = -ifa-1;
        fa_sign = -1.0;
      }
      phi_cf[icf] += fa_sign*solver->area_fa[ifa]*phi_fa[ifa];
    }
  }
  updateSignedCfIDataStart(phi_cf); 

  for (int icf = 0; icf < ncf_aa; ++icf) {
    phi_cf[icf] = 0.0;
    for (int fof = faocf_i[icf]; fof < faocf_i[icf+1]; ++fof) {
      int ifa = faocf_v[fof];
      double fa_sign = 1.0;
      if (ifa < 0) {
        ifa = -ifa-1;
        fa_sign = -1.0;
      }
      phi_cf[icf] += fa_sign*solver->area_fa[ifa]*phi_fa[ifa];
    }
  }
  updateCfIDataFinish(phi_cf); // same as for unsigned cf data

  for (int icf = 0; icf < ncf_ai; ++icf)
    phi_cf[icf] /= area_cf[icf];

}

void StaticSolver::CoarseGrid::restrictExtrinsicSignedCfData(double* phiA_cf,const double* phiA_fa) {

  for (int icf = ncf_aa; icf < ncf; ++icf) {
    phiA_cf[icf] = 0.0;
    for (int fof = faocf_i[icf]; fof < faocf_i[icf+1]; ++fof) {
      int ifa = faocf_v[fof];
      double fa_sign = 1.0;
      if (ifa < 0) {
        ifa = -ifa-1;
        fa_sign = -1.0;
      }
      phiA_cf[icf] += fa_sign*phiA_fa[ifa];
    }
  }
  updateSignedCfIDataStart(phiA_cf); 

  for (int icf = 0; icf < ncf_aa; ++icf) {
    phiA_cf[icf] = 0.0;
    for (int fof = faocf_i[icf]; fof < faocf_i[icf+1]; ++fof) {
      int ifa = faocf_v[fof];
      double fa_sign = 1.0;
      if (ifa < 0) {
        ifa = -ifa-1;
        fa_sign = -1.0;
      }
      phiA_cf[icf] += fa_sign*phiA_fa[ifa];
    }
  }
  updateCfIDataFinish(phiA_cf); // same as for unsigned cf data

}

void StaticSolver::CoarseGrid::restrictCbData(double* phi_cb,const double* phi_bf,const int izone) {

  if (izone == -1) {
    // full boundary data...

    // restrict inactive boundaries and send to active

    for (int icb = ncb_a; icb < ncb; ++icb) {
      phi_cb[icb] = 0.0;
      for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
        const int ibf = bfocb_v[bob];
        phi_cb[icb] += solver->area_bf[ibf]*phi_bf[ibf];
      }
    }
    updateCbIDataStart(phi_cb);

    // restrict active boundaries

    for (int icb = 0; icb < ncb_a; ++icb) {
      phi_cb[icb] = 0.0;
      for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
        const int ibf = bfocb_v[bob];
        phi_cb[icb] += solver->area_bf[ibf]*phi_bf[ibf];
      }
    }
    updateCbIDataFinish(phi_cb);

    // normalize

    for (int icb = 0; icb < ncb_a; ++icb)
      phi_cb[icb] /= area_cb[icb];

  }
  else {
    // for now all coarse grid boundary data is allocated everywhere
    // TODO better zone management - will need to reorder boundary data again

    if (solver->cg_level == -1) {
      // when your parent is the fine grid, you restrict using zone local ibf's
      for (int icb = ncb_a; icb < ncb; ++icb) {
        phi_cb[icb] = 0.0;
        if (zone_cb[icb] == izone) {
          for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
            const int ibf = bfocb_v[bob];
            assert(solver->zone_bf[ibf] == izone);
            phi_cb[icb] += solver->area_bf[ibf]*phi_bf[ibf-solver->bfZoneVec[izone].ibf_f];
          }
        }
      }
      updateCbIDataStart(phi_cb);

      for (int icb = 0; icb < ncb_a; ++icb) {
        phi_cb[icb] = 0.0;
        if (zone_cb[icb] == izone) {
          for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
            const int ibf = bfocb_v[bob];
            assert(solver->zone_bf[ibf] == izone);
            phi_cb[icb] += solver->area_bf[ibf]*phi_bf[ibf-solver->bfZoneVec[izone].ibf_f];
          }
        }
      }
      updateCbIDataFinish(phi_cb);

      for (int icb = cbozn_i[izone]; icb < cbozn_i[izone+1]; ++icb) 
        phi_cb[icb] /= area_cb[icb];
    }
    else { 
      // when your parent is a coarse grid, you restrict using global ibf's (proc local)
      for (int icb = ncb_a; icb < ncb; ++icb) {
        phi_cb[icb] = 0.0;
        if (zone_cb[icb] == izone) {
          for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
            const int ibf = bfocb_v[bob];
            assert(solver->zone_bf[ibf] == izone);
            phi_cb[icb] += solver->area_bf[ibf]*phi_bf[ibf];
          }
        }
      }
      updateCbIDataStart(phi_cb);

      for (int icb = 0; icb < ncb_a; ++icb) {
        phi_cb[icb] = 0.0;
        if (zone_cb[icb] == izone) {
          for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
            const int ibf = bfocb_v[bob];
            assert(solver->zone_bf[ibf] == izone);
            phi_cb[icb] += solver->area_bf[ibf]*phi_bf[ibf];
          }
        }
      }
      updateCbIDataFinish(phi_cb);

      for (int icb = cbozn_i[izone]; icb < cbozn_i[izone+1]; ++icb) 
        phi_cb[icb] /= area_cb[icb];
    }

  }

}

void StaticSolver::CoarseGrid::restrictExtrinsicCbData(double* phiA_cb,const double* phiA_bf,const int izone) {

  if (izone == -1) {
    // full boundary data...

    // restrict inactive boundaries and send to active

    for (int icb = ncb_a; icb < ncb; ++icb) {
      phiA_cb[icb] = 0.0;
      for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
        const int ibf = bfocb_v[bob];
        phiA_cb[icb] += phiA_bf[ibf];
      }
    }
    updateCbIDataStart(phiA_cb);

    // restrict active boundaries

    for (int icb = 0; icb < ncb_a; ++icb) {
      phiA_cb[icb] = 0.0;
      for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
        const int ibf = bfocb_v[bob];
        phiA_cb[icb] += phiA_bf[ibf];
      }
    }
    updateCbIDataFinish(phiA_cb);

  }
  else {
    // for now all coarse grid boundary data is allocated everywhere
    // TODO better zone management - will need to reorder boundary data again

    if (solver->cg_level == -1) {
      // when your parent is the fine grid, you restrict using zone local ibf's
      for (int icb = ncb_a; icb < ncb; ++icb) {
        phiA_cb[icb] = 0.0;
        if (zone_cb[icb] == izone) {
          for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
            const int ibf = bfocb_v[bob];
            assert(solver->zone_bf[ibf] == izone);
            phiA_cb[icb] += phiA_bf[ibf-solver->bfZoneVec[izone].ibf_f];
          }
        }
      }
      updateCbIDataStart(phiA_cb);

      for (int icb = 0; icb < ncb_a; ++icb) {
        phiA_cb[icb] = 0.0;
        if (zone_cb[icb] == izone) {
          for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
            const int ibf = bfocb_v[bob];
            assert(solver->zone_bf[ibf] == izone);
            phiA_cb[icb] += phiA_bf[ibf-solver->bfZoneVec[izone].ibf_f];
          }
        }
      }
      updateCbIDataFinish(phiA_cb);

    }
    else { 
      // when your parent is a coarse grid, you restrict using global ibf's (proc local)
      for (int icb = ncb_a; icb < ncb; ++icb) {
        phiA_cb[icb] = 0.0;
        if (zone_cb[icb] == izone) {
          for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
            const int ibf = bfocb_v[bob];
            assert(solver->zone_bf[ibf] == izone);
            phiA_cb[icb] += phiA_bf[ibf];
          }
        }
      }
      updateCbIDataStart(phiA_cb);

      for (int icb = 0; icb < ncb_a; ++icb) {
        phiA_cb[icb] = 0.0;
        if (zone_cb[icb] == izone) {
          for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
            const int ibf = bfocb_v[bob];
            assert(solver->zone_bf[ibf] == izone);
            phiA_cb[icb] += phiA_bf[ibf];
          }
        }
      }
      updateCbIDataFinish(phiA_cb);

    }

  }

}

void StaticSolver::CoarseGrid::restrictCbData(double (*u_cb)[3],const double (*u_bf)[3],const int izone) {

  if (izone == -1) {
    // full boundary data...

    for (int icb = ncb_a; icb < ncb; ++icb) {
      FOR_I3 u_cb[icb][i] = 0.0;
      for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
        const int ibf = bfocb_v[bob];
        FOR_I3 u_cb[icb][i] += solver->area_bf[ibf]*u_bf[ibf][i];
      }
    }
    updateCbIDataStart(u_cb);

    for (int icb = 0; icb < ncb_a; ++icb) {
      FOR_I3 u_cb[icb][i] = 0.0;
      for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
        const int ibf = bfocb_v[bob];
        FOR_I3 u_cb[icb][i] += solver->area_bf[ibf]*u_bf[ibf][i];
      }
    }
    updateCbIDataFinish(u_cb);

    for (int icb = 0; icb < ncb_a; ++icb)
      FOR_I3 u_cb[icb][i] /= area_cb[icb];

  }
  else {
    // for now all coarse grid boundary data is allocated everywhere
    // TODO better zone management - will need to reorder boundary data again

    if (solver->cg_level == -1) {
      for (int icb = ncb_a; icb < ncb; ++icb) {
        FOR_I3 u_cb[icb][i] = 0.0;
        if (zone_cb[icb] == izone) {
          for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
            const int ibf = bfocb_v[bob];
            assert(solver->zone_bf[ibf] == izone);
            FOR_I3 u_cb[icb][i] += solver->area_bf[ibf]*u_bf[ibf-solver->bfZoneVec[izone].ibf_f][i];
          }
        }
      }
      updateCbIDataStart(u_cb);

      for (int icb = 0; icb < ncb_a; ++icb) {
        FOR_I3 u_cb[icb][i] = 0.0;
        if (zone_cb[icb] == izone) {
          for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
            const int ibf = bfocb_v[bob];
            assert(solver->zone_bf[ibf] == izone);
            FOR_I3 u_cb[icb][i] += solver->area_bf[ibf]*u_bf[ibf-solver->bfZoneVec[izone].ibf_f][i];
          }
        }
      }
      updateCbIDataFinish(u_cb);

      for (int icb = cbozn_i[izone]; icb < cbozn_i[izone+1]; ++icb) 
        FOR_I3 u_cb[icb][i] /= area_cb[icb];
    }
    else {
      for (int icb = ncb_a; icb < ncb; ++icb) {
        FOR_I3 u_cb[icb][i] = 0.0;
        if (zone_cb[icb] == izone) {
          for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
            const int ibf = bfocb_v[bob];
            assert(solver->zone_bf[ibf] == izone);
            FOR_I3 u_cb[icb][i] += solver->area_bf[ibf]*u_bf[ibf-solver->bfZoneVec[izone].ibf_f][i];
          }
        }
      }
      updateCbIDataStart(u_cb);

      for (int icb = 0; icb < ncb_a; ++icb) {
        FOR_I3 u_cb[icb][i] = 0.0;
        if (zone_cb[icb] == izone) {
          for (int bob = bfocb_i[icb]; bob < bfocb_i[icb+1]; ++bob) {
            const int ibf = bfocb_v[bob];
            assert(solver->zone_bf[ibf] == izone);
            FOR_I3 u_cb[icb][i] += solver->area_bf[ibf]*u_bf[ibf-solver->bfZoneVec[izone].ibf_f][i];
          }
        }
      }
      updateCbIDataFinish(u_cb);

      for (int icb = cbozn_i[izone]; icb < cbozn_i[izone+1]; ++icb) 
        FOR_I3 u_cb[icb][i] /= area_cb[icb];
    }

  }

}

void StaticSolver::initFromCoarseGrid(StaticSolver::CoarseGrid* cg) {

  // NOTE: you MUST allocate prolongated and/or restricted arrays using 
  // the cg->ncc_g,cg->ncb,cg->ncf counts (e.g., your restricted res or prolongated err) 

  // control volume data...
  assert(ncv_global == 0); ncv_global = cg->ncc_global;
  assert(ncv_i == 0); ncv_i = 0;
  assert(ncv == 0); ncv = cg->ncc_a;
  assert(ncv_g == 0); ncv_g = ncv; 

  // flag cc's that are touched by icf = ncf_aa:ncf_ai

  int* reorder_g = new int[cg->ncc_g-cg->ncc_a];
  for (int icc = cg->ncc_a; icc < cg->ncc_g; ++icc)
    reorder_g[icc-cg->ncc_a] = -1;
  for (int icf = cg->ncf_aa; icf < cg->ncf_ai; ++icf) {
    assert((cg->ccocf[icf][1] >= cg->ncc_a)&&(cg->ccocf[icf][1] < cg->ncc_g));
    reorder_g[cg->ccocf[icf][1]-cg->ncc_a] = -2;
  }
  for (int icc = cg->ncc_a; icc < cg->ncc_g; ++icc) {
    if (reorder_g[icc-cg->ncc_a] == -2)
      ++ncv_g;
  }

  // need to rbi order ncc_a:ncc_g in order to have 
  // a single prcomm (is their a better way?)

  vector<pair<int8,int> > rbi_index_pair_vec;
  for (int icc = cg->ncc_a; icc < cg->ncc; ++icc) {
    if (reorder_g[icc-cg->ncc_a] == -2)
      rbi_index_pair_vec.push_back(pair<int8,int>(cg->rbi_i[icc-cg->ncc_a],icc));
  }
  for (int icc = cg->ncc; icc < cg->ncc_g; ++icc) {
    if (reorder_g[icc-cg->ncc_a] == -2)
      rbi_index_pair_vec.push_back(pair<int8,int>(cg->rbi_g[icc-cg->ncc],icc));
  }
  assert(int(rbi_index_pair_vec.size()) == ncv_g-ncv);
  sort(rbi_index_pair_vec.begin(),rbi_index_pair_vec.end());
  assert(rbi_g == NULL); rbi_g = new uint8[ncv_g-ncv];
  for (int ii = 0, lim = rbi_index_pair_vec.size(); ii < lim; ++ii) {
    rbi_g[ii] = rbi_index_pair_vec[ii].first;
    reorder_g[rbi_index_pair_vec[ii].second-ncv] = ii+ncv;
  }

  assert(icv_global == NULL); icv_global = new int8[ncv_g]; 
  assert(vol_cv == NULL); vol_cv = new double[ncv_g]; 
  assert(inv_vol == NULL); inv_vol = new double[ncv_g];
  assert(x_cv == NULL); x_cv = new double[ncv_g][3];
  FOR_ICV {
    icv_global[icv] = cg->icc_global[icv]; 
    vol_cv[icv] = cg->vol_cc[icv];
    inv_vol[icv] = 1.0/vol_cv[icv];
    FOR_I3 x_cv[icv][i] = cg->x_cc[icv][i];
  }
  for (int icv = ncv; icv < cg->ncc_g; ++icv) {
    if (reorder_g[icv-ncv] != -1) {
      const int icv_new = reorder_g[icv-ncv]; assert((icv_new >= ncv)&&(icv_new < ncv_g));
      icv_global[icv_new] = cg->icc_global[icv]; 
      vol_cv[icv_new] = cg->vol_cc[icv];
      inv_vol[icv_new] = 1.0/vol_cv[icv_new];
      FOR_I3 x_cv[icv_new][i] = cg->x_cc[icv][i];
    }
  }

  // cv communicator

  buildCvPrcomm(); 

  // face data...

  assert(nfa_i == 0); nfa_i = cg->ncf_aa;
  assert(nfa == 0); nfa = cg->ncf_ai;
  assert(cvofa == NULL); cvofa = new int[nfa][2];
  assert(n_fa == NULL); n_fa = new double[nfa][3];
  assert(x_fa == NULL); x_fa = new double[nfa][3];
  assert(area_fa == NULL); area_fa = new double[nfa]; 
  assert(area_over_delta_fa == NULL); area_over_delta_fa = new double[nfa];
  FOR_IFA {
    const int icv0 = cg->ccocf[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
    const int icv1 = cg->ccocf[ifa][1]; assert((icv0 >= 0)&&(icv0 < ncv_g));
    cvofa[ifa][0] = icv0;
    if (icv1 < ncv) {
      assert(ifa < nfa_i);
      cvofa[ifa][1] = icv1;
    }
    else {
      assert(ifa >= nfa_i);
      cvofa[ifa][1] = reorder_g[icv1-ncv];
    }
    FOR_I3 x_fa[ifa][i] = cg->x_cf[ifa][i];
    FOR_I3 n_fa[ifa][i] = cg->n_cf[ifa][i];
    area_fa[ifa] = cg->area_cf[ifa];
    area_over_delta_fa[ifa] = area_fa[ifa]/cg->delta_cf[ifa];
  }
  delete[] reorder_g;
  // for now just store sign (+1/-1)
  assert(ifa_global == NULL); ifa_global = new int8[nfa];
  FOR_INTERNAL_IFA ifa_global[ifa] = 1;
  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
    // already reordered above
    const int icv1 = cvofa[ifa][1]; assert((icv0 >= 0)&&(icv0 < ncv_g));
    const int icv0_global = icv_global[icv0]; 
    const int icv1_global = icv_global[icv1]; 
    if (icv0_global > icv1_global) 
      ifa_global[ifa] = 1;
    else 
      ifa_global[ifa] = -1;
  }

  // cvocv stuff...

  buildCvocv();

  assert(cvocv_grad_coeff == NULL);
  cvocv_grad_coeff = new double[cvocv_i[ncv]][3];
  for (int coc = 0; coc < cvocv_i[ncv]; ++coc)
    FOR_I3 cvocv_grad_coeff[coc][i] = 0.0;
  for (int ifa = 0; ifa < nfa; ++ifa) {
    const int icv0 = cvofa[ifa][0];
    const int coc00 = cvocv_i[icv0];
    const int icv1 = cvofa[ifa][1];
    int coc01 = coc00+1;
    while (cvocv_v[coc01] != icv1) {
      ++coc01;
      assert(coc01 != cvocv_i[icv0+1]);
    }

    FOR_I3 cvocv_grad_coeff[coc00][i] += 0.5*n_fa[ifa][i];
    FOR_I3 cvocv_grad_coeff[coc01][i] -= 0.5*n_fa[ifa][i];

    if (icv1 < ncv) {
      const int coc11 = cvocv_i[icv1];
      int coc10 = coc11+1;
      while (cvocv_v[coc10] != icv0) {
        ++coc10;
        assert(coc10 != cvocv_i[icv1+1]);
      }
      FOR_I3 cvocv_grad_coeff[coc11][i] += 0.5*n_fa[ifa][i];
      FOR_I3 cvocv_grad_coeff[coc10][i] -= 0.5*n_fa[ifa][i];
    }

  }
  FOR_ICV {
    for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc)
      FOR_I3 cvocv_grad_coeff[coc][i] *= inv_vol[icv];
  }

  // boundary data...

  assert(nbf == 0); nbf = cg->ncb_a;
  assert(zone_bf == NULL); zone_bf = new int[nbf];
  assert(cvobf == NULL); cvobf = new int[nbf];
  assert(n_bf == NULL); n_bf = new double[nbf][3];
  assert(x_bf == NULL); x_bf = new double[nbf][3]; 
  assert(area_bf == NULL); area_bf = new double[nbf];
  assert(area_over_delta_bf == NULL); area_over_delta_bf = new double[nbf];
  FOR_IBF {
    zone_bf[ibf] = cg->zone_cb[ibf];
    const int icv =  cg->ccocb[ibf]; assert((icv >= 0)&&(icv < ncv));
    cvobf[ibf] = icv;
    FOR_I3 x_bf[ibf][i] = cg->x_cb[ibf][i];
    FOR_I3 n_bf[ibf][i] = cg->n_cb[ibf][i];
    area_bf[ibf] = cg->area_cb[ibf];
    area_over_delta_bf[ibf] = area_bf[ibf]/cg->delta_cb[ibf];
  }

  assert(bfZoneVec.empty());
  const int nzn = cg->solver->bfZoneVec.size();
  for (int izn = 0; izn < nzn; ++izn) {
    bfZoneVec.push_back(BfZone());
    bfZoneVec.back().setName(cg->solver->bfZoneVec[izn].getName());
    bfZoneVec.back().area_global = cg->solver->bfZoneVec[izn].area_global;
    FOR_I3 bfZoneVec.back().n_global[i] = cg->solver->bfZoneVec[izn].n_global[i];
    FOR_I3 bfZoneVec.back().x_global[i] = cg->solver->bfZoneVec[izn].x_global[i];
    bfZoneVec.back().index = izn;
    bfZoneVec.back().nbf = cg->cbozn_i[izn+1]-cg->cbozn_i[izn];
    bfZoneVec.back().ibf_f = cg->cbozn_i[izn];
    const int ibf_f = bfZoneVec.back().ibf_f;
    bfZoneVec.back().ibf_l = cg->cbozn_i[izn+1]-1;
    bfZoneVec.back().cvobf = cvobf+ibf_f;
    bfZoneVec.back().n_bf = n_bf+ibf_f;
    bfZoneVec.back().x_bf = x_bf+ibf_f;
    bfZoneVec.back().area_bf = area_bf+ibf_f;
    bfZoneVec.back().area_over_delta_bf = area_over_delta_bf+ibf_f;
  }

}
