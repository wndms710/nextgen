#ifndef LAMINARPREMIXED_HPP
#define LAMINARPREMIXED_HPP


class LaminarPremixed : public PremixedSolver {
public:

  void initialHook() {
    
    if (step == 0) {
      
      COUT1(" > setting initial condition"); 
      
      const double x_f   = getDoubleParam("X_FLAME");
      const double Z_g   = getDoubleParam("Z_GLOBAL");
      const double p_ref = this->chemtable->pressure; 

      double sL, lF, C_max, C_min = -1e+20;
      this->chemtable->lookup(&C_min,"prog",&Z_g,&C_min,1);
      this->chemtable->lookupReduced(&sL,"sL",&lF,"lF",&C_max,"upper",&Z_g,1); 

      double R, T, rho0;
      this->chemtable->lookup(&R,"R",&T,"T",&Z_g,&C_min,1);
      rho0 = p_ref/R/T; 

      if (mpi_rank == 0) { 
        cout << "   > Z_g   = " << Z_g    << endl; 
        cout << "   > C_min = " << C_min  << endl; 
        cout << "   > C_max = " << C_max  << endl; 
        cout << "   > rho0  = " << rho0   << endl; 
        cout << "   > sL    = " << sL     << endl; 
        cout << "   > lF    = " << lF     << endl; 
      }
     
      for (int icv = 0; icv < ncv; ++icv) { 
        
        p[icv] = p_ref;
        Z[icv] = Z_g;

        if (x_cv[icv][0] < x_f-lF/2.0) {
          C[icv] = C_min; 
        }
        else if (x_cv[icv][0] > x_f+lF/2.0) {
          C[icv] = C_max; 
        }
        else {
          C[icv] = C_min + (x_cv[icv][0] - (x_f-lF/2.0))*(C_max-C_min)/lF; 
        }

        double e_, R_, T_; 
        this->chemtable->lookup(&e_,"e",&R_,"R",&T_,"T",&Z_g,&C[icv],1); 
        
        rho[icv]  = p[icv]/R_/T_; 
        u[icv][0] = sL*rho0/rho[icv];
        u[icv][1] = 0.0; 
        u[icv][2] = 0.0; 
        rhoE[icv] = rho[icv]*(e_ + 0.5*DOT_PRODUCT(u[icv],u[icv])); 
      }

      updateConservativeAndPrimitiveData();
    }
  }

  /*
  void temporalHook() {
    if (step % check_interval == 0) {

      const int nvar=7;
      int     n_cv_g; 
      int     *rdisp       = NULL;
      int     *rcounts     = NULL;
      double (*buf0)[nvar] = NULL;
      double (*buf )[nvar] = new double[ncv][nvar]; 

      if (mpi_rank == 0) rcounts = new int[mpi_size];
      MPI_Gather(&ncv,1,MPI_INT,rcounts,1,MPI_INT,0,mpi_comm); 

      if (mpi_rank == 0) {
        for (int i=0; i<mpi_size; ++i)
          rcounts[i] *= nvar; 
        rdisp = new int[mpi_size]; 
        rdisp[0] = 0; 
        for (int i=1; i<mpi_size; ++i)
          rdisp[i] = rdisp[i-1] + rcounts[i-1];
        n_cv_g = (rdisp[mpi_size-1] + rcounts[mpi_size-1])/nvar;
        buf0  = new double[n_cv_g][nvar]; 
      }

      FOR_ICV {
        buf[icv][0] = x_cv[icv][0];
        buf[icv][1] =    u[icv][0];
        buf[icv][2] =  rho[icv];
        buf[icv][3] =    p[icv];
        buf[icv][4] =    T[icv];
        buf[icv][5] =    Z[icv];
        buf[icv][6] =    C[icv];
      }
      MPI_Gatherv(buf,nvar*ncv,MPI_DOUBLE,buf0,rcounts,rdisp,MPI_DOUBLE,0,mpi_comm);

      if (mpi_rank == 0) {
        cout << ">>>>> x, u, rho, p, T, Z, C" << endl;
        for (int i=0; i<n_cv_g; ++i) {
          cout << "XXXXX";
          for (int j=0; j<nvar; ++j) cout << " " << buf0[i][j]; 
          cout << endl; 
        }
      }
      delete[] rcounts;
      delete[] rdisp;
      delete[] buf0; 
      delete[] buf; 
    }
  }
  */
};


#endif
