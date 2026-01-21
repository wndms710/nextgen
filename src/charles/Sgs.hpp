#ifndef SGS_HPP
#define SGS_HPP


namespace DynamicSmagorinsky { 

  inline void registerDsmLocalVariables(FlowSolver* solver) { 

    solver->_registerScalar("lijmij",1); // not transported .. 
    solver->_registerScalar("mijmij",1); 
    solver->_registerScalar("c_smag",1);
    solver->_registerScalar("c_clip",1);

  }

  inline void computeLijMij(double * lijmij, double * mijmij, 
                            const double * rho, const double (*u)[3], 
                            double (*dudx)[3][3], FlowSolver* solver) { 
    
    
    const double deltar        = getDoubleParam("DELTA_R_DSM", 2.0);
    const double alpha         = deltar*deltar; 
    const int ncv_g            = solver->ncv_g;   
    const int ncv              = solver->ncv;   
    
    double (*lij_d)[3]  = new double[ncv_g][3];
    double (*lij_od)[3] = new double[ncv_g][3];
    double (*mij_d)[3]  = new double[ncv_g][3];
    double (*mij_od)[3] = new double[ncv_g][3];


    // perform the various filtrations necessary for the dsm.
    
    double *rho_hat    = new double[ncv_g];
    double (*u_hat)[3] = new double[ncv_g][3];
    
    // note that the following proceeds assuming that the variables 
    // are all reynolds averaged (not favre averaged).  this is the 
    // consistent interpretation as we are admitting modeling terms 
    // in the continuity equation.  (although, neglecting other terms
    // throughout the formulation)
    
    solver->filterCvR1(rho_hat,rho);
    solver->filterCvR2(u_hat  , u);
    
    solver->updateCvData(rho_hat);
    solver->updateCvData(u_hat, REPLACE_ROTATE_DATA);

    double (*dudx_hat)[3][3] = new double[ncv][3][3];
    solver->StaticSolver::calcCvGrad(dudx_hat,u_hat);

    /*
    MiscUtils::dumpRange(rho,ncv,"rho");
    MiscUtils::dumpRange(u,ncv,"u");
    MiscUtils::dumpRange(dudx,ncv,"dudx");

    MiscUtils::dumpRange(rho_hat,ncv,"rho_hat");
    MiscUtils::dumpRange(u_hat,ncv,"u_hat");
    MiscUtils::dumpRange(dudx_hat,ncv,"dudx_hat");
    */

    // note that at this juncture you should have access to a 
    // completed gradient calculation at the grid level.

    // we also need two additional buffers .. 

    double (*tmp_d)[3]  = new double[ncv_g][3];
    double (*tmp_od)[3] = new double[ncv_g][3];

    for (int icv = 0; icv < ncv; ++icv) { 

      double smag = 0.0;
      double Sij[3][3];
      for (int i = 0; i < 3; ++i) {  
        for (int j = 0; j < 3; ++j) { 
          Sij[i][j] = 0.5*(dudx[icv][i][j] + dudx[icv][j][i]);
          smag += 2.0*Sij[i][j]*Sij[i][j];
        }
      }
      smag = sqrt(smag);

      // remove the trace from sij
      double Skk = Sij[0][0] + Sij[1][1] + Sij[2][2];
      for (int i = 0; i < 3; ++i)
        Sij[i][i] -= Skk/3.0;

      
      // augment the smag to include a factor of the density .. 

      smag *= rho[icv];

      for(int i = 0; i < 3; ++i) 
        lij_d[icv][i] = smag*Sij[i][i];

      lij_od[icv][0] = smag*Sij[1][2];
      lij_od[icv][1] = smag*Sij[0][2];
      lij_od[icv][2] = smag*Sij[0][1];
      
    }
    
    // the following should be checked for rotational periodicity 
    // of the symmetric tensor update.. 
    
    solver->updateCvData(lij_d);
    solver->updateCvData(lij_od);

    solver->filterCvR2(mij_d , lij_d);
    solver->filterCvR2(mij_od, lij_od);

    // start the calculation of L .. 

    for (int icv = 0; icv < ncv ; ++icv)  { 

      for (int i = 0; i < 3; ++i) 
        tmp_d[icv][i] = rho[icv]*u[icv][i]*u[icv][i];

      tmp_od[icv][0]  = rho[icv]*u[icv][1]*u[icv][2];
      tmp_od[icv][1]  = rho[icv]*u[icv][0]*u[icv][2];
      tmp_od[icv][2]  = rho[icv]*u[icv][0]*u[icv][1];

    }

    solver->updateCvData(tmp_d);
    solver->updateCvData(tmp_od);

    solver->filterCvR2(lij_d, tmp_d);
    solver->filterCvR2(lij_od,tmp_od);

    // now we can complete the calculations of L, M

    for (int icv = 0; icv < ncv; ++icv) { 

      double smag_hat = 0.0;
      double Sij_hat[3][3];
      for (int i = 0; i < 3; ++i) {  
        for (int j = 0; j < 3; ++j) { 
          Sij_hat[i][j] = 0.5*(dudx_hat[icv][i][j] + dudx_hat[icv][j][i]);
          smag_hat += 2.0*Sij_hat[i][j]*Sij_hat[i][j];
        }
      }
      smag_hat = sqrt(smag_hat);

      // remove the trace from S_hat.. 
      double Skk_hat = Sij_hat[0][0] + Sij_hat[1][1] + Sij_hat[2][2];
      for (int i = 0; i < 3; ++i)
        Sij_hat[i][i] -= Skk_hat/3.0;

      
      // augment the smag_hat to include a factor of the density .. 

      smag_hat *= rho_hat[icv];

      for (int i = 0; i < 3; ++i)
        mij_d[icv][i] = alpha*smag_hat*Sij_hat[i][i] - mij_d[icv][i];

      mij_od[icv][0]  = alpha*smag_hat*Sij_hat[1][2] - mij_od[icv][0];
      mij_od[icv][1]  = alpha*smag_hat*Sij_hat[0][2] - mij_od[icv][1];
      mij_od[icv][2]  = alpha*smag_hat*Sij_hat[0][1] - mij_od[icv][2];

      // finish the calculation of L.. 

      for (int i =0; i < 3; ++i)
        lij_d[icv][i] -= rho_hat[icv]*u_hat[icv][i]*u_hat[icv][i];

      lij_od[icv][0]  -= rho_hat[icv]*u_hat[icv][1]*u_hat[icv][2];
      lij_od[icv][1]  -= rho_hat[icv]*u_hat[icv][0]*u_hat[icv][2];
      lij_od[icv][2]  -= rho_hat[icv]*u_hat[icv][0]*u_hat[icv][1];

      // remove the trace of L .. 

      const double Lkk = lij_d[icv][0] + lij_d[icv][1] + lij_d[icv][2];
      for (int i = 0; i < 3; ++i)
        lij_d[icv][i] -= Lkk/3.0;

    }

    // perform some time averaging of lijmij, mijmij... 
  
    if ( solver->dt > 0.0 ) { 

      for (int icv = 0; icv < ncv; ++icv) { 
        
        // use the inv strain rate to set the local relaxation time scale
        
        double smag = 0.0; 
        double Sij[3][3];
        for (int i = 0; i < 3; ++i) {  
          for (int j = 0; j < 3; ++j) { 
            Sij[i][j] = 0.5*(dudx[icv][i][j] + dudx[icv][j][i]);
            smag += 2.0*Sij[i][j]*Sij[i][j];
          }
        }
        smag                 = sqrt(smag) + 1.0e-15;
        const double T_relax = 1.0/smag;
        //const double tt      = max(0.0,-lijmij[icv]);
        //const double T_relax = 2.0*pow(vol_cv[icv],1.0/3.0)/(pow(tt,1.0/4.0)+1.0e-12);
        const double epsln   = solver->dt/T_relax;
        
        // build the local contraction of lijmij, mijmij .. 
        
        double lijmij_ = 0.0;
        double mijmij_ = 1.0e-12;  // set to epsilon to avoid div (zero)
        
        for (int i = 0; i < 3; ++i) { 
          
          lijmij_ +=     lij_d[icv][i]*mij_d[icv][i];
          lijmij_ += 2.0*lij_od[icv][i]*mij_od[icv][i];
          
          mijmij_ +=     mij_d[icv][i]*mij_d[icv][i];
          mijmij_ += 2.0*mij_od[icv][i]*mij_od[icv][i];
          
        }
       
        // these should be non-zero if you wanted to use lagragian averaging .. 

        double dlm     = 0.0; //dt*DOT_PRODUCT(u[icv],tmp_d[icv]);
        double dmm     = 0.0; //dt*DOT_PRODUCT(u[icv],tmp_od[icv]);

        const double rhs_lm  = min(0.0,lijmij[icv] - dlm + epsln*lijmij_);
        const double rhs_mm  = mijmij[icv] - dmm + epsln*mijmij_;

        const double lhs_fax = 1.0/(1.0 + epsln);

        lijmij[icv]          = rhs_lm*lhs_fax;
        mijmij[icv]          = rhs_mm*lhs_fax;

      }

    }
  
    // cleanup ... 

    delete[] tmp_d;
    delete[] tmp_od;
    delete[] rho_hat;
    delete[] u_hat;
    delete[] dudx_hat;
    delete[] lij_d;
    delete[] lij_od;
    delete[] mij_d;
    delete[] mij_od;


  }


}

#endif
