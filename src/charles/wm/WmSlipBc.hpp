
void WmSlipApprox::initData() { 

  assert( cdel_w == NULL); cdel_w = new double[zone_ptr->nbf];
  assert( lijpij == NULL); lijpij = new double[zone_ptr->nbf];
  assert( pijpij == NULL); pijpij = new double[zone_ptr->nbf];
  assert( rhou_s == NULL); rhou_s = new double[zone_ptr->nbf][3];

}

void WmSlipApprox::initialHook() { 

  if ( !zone_ptr->checkDataFlag("cdel_w") ||
       !zone_ptr->checkDataFlag("lijpij") || 
       !zone_ptr->checkDataFlag("pijpij") || 
       !zone_ptr->checkDataFlag("rhou_s") ) { 

    const double epsln = 1.0e-12;

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

      cdel_w[ibf] = 0.0;
      lijpij[ibf] = 0.0;
      pijpij[ibf] = epsln;
      
      for (int i = 0; i < 3; ++i)
        rhou_s[ibf][i] = 0.0;

    }
  }

}

void WmSlipApprox::calcRhs(const double time, const int rk_stage) {

  if ( rk_stage == 1) { 


    // compute the slip length and set the rhou_s at the wall ... 
    // there is likely a more efficient mechanism through which this 
    // can be performed by grabbing some previously computed information
    // from the dynamic procedure .. 

    if ( checkParam("CDEL_W_CONSTANT")) { 

      // the following is here mostly for debugging purposes ... 

      const double cdel_const = getDoubleParam("CDEL_W_CONSTANT");

      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

        const int icv      = zone_ptr->cvobf[ibf];
        const double delta = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        cdel_w[ibf]        = cdel_const*delta;

        for (int i = 0; i < 3; ++i)  
          rhou_s[ibf][i] = cdel_w[ibf] * solver->rho[icv] * solver->u[icv][i] / (delta + cdel_w[ibf]);

      }

    } else { 

      // compute the incompressible approximation of the dynamic slip length 
      
      const double delta_r = getDoubleParam("WM_SLIP_WALL:DELTA_R", 1.5);

      double (*u_hat)[3]  = new double[solver->ncv_g][3];
      double (*lij_d)[3]  = new double[solver->ncv_g][3];
      double (*lij_od)[3] = new double[solver->ncv_g][3]; 
      double (*tmp_d)[3]  = new double[solver->ncv_g][3];
      double (*tmp_od)[3] = new double[solver->ncv_g][3];

      solver->filterCvR2(u_hat,solver->u);
      solver->updateCvData(u_hat, REPLACE_ROTATE_DATA);

      for (int icv = 0; icv < solver->ncv_g; ++icv) { 

        for (int i = 0; i < 3; ++i) 
          tmp_d[icv][i] = solver->u[icv][i]*solver->u[icv][i];

        tmp_od[icv][0]  = solver->u[icv][1]*solver->u[icv][2];
        tmp_od[icv][1]  = solver->u[icv][2]*solver->u[icv][0];
        tmp_od[icv][2]  = solver->u[icv][0]*solver->u[icv][1];

      }

      solver->filterCvR2(lij_d, tmp_d);
      solver->filterCvR2(lij_od, tmp_od);
      
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

        const int icv = zone_ptr->cvobf[ibf];
        
        // complete lij ---

        for (int i = 0; i < 3; ++i) 
          lij_d[icv][i] -= solver->u[icv][i]*solver->u[icv][i];

        lij_od[icv][0]  -= solver->u[icv][1]*solver->u[icv][2];
        lij_od[icv][1]  -= solver->u[icv][2]*solver->u[icv][0];
        lij_od[icv][2]  -= solver->u[icv][0]*solver->u[icv][1];

        
        // now compute the gradient of u_hat in this cell .. 

        double dudx_hat[3][3];

        for (int i = 0; i <3; ++i) 
          for (int j = 0; j < 3; ++j)
            dudx_hat[i][j] = 0.0;

        for (int coc = solver->cvocv_i[icv]; coc < solver->cvocv_i[icv+1]; ++coc) { 

          const int icv_nbr = solver->cvocv_v[coc];
          for (int i = 0; i < 3; ++i) 
            for (int j = 0; j < 3; ++j) 
              dudx_hat[i][j] += solver->cvocv_grad_coeff[coc][j]*u_hat[icv_nbr][i];
        }


        const double mag_n    = MAG(zone_ptr->n_bf[ibf]);
        double unit_n[3];
        
        if ( mag_n > 0.0) { 
          
          for (int i = 0; i < 3; ++i) 
            unit_n[i] = zone_ptr->n_bf[ibf][i]/mag_n;
          
        } else { 
          
          for (int i = 0; i < 3; ++i) 
            unit_n[i] = 0.0;
          
        }
        

        double duidn_hat[3], duidn[3];

        for (int i = 0; i < 3; ++i) { 

          duidn_hat[i] = -DOT_PRODUCT(dudx_hat[i],unit_n);
          duidn[i]     = -DOT_PRODUCT(solver->dudx[icv][i],unit_n);

        }

        double pij_d[3], pij_od[3];

        for (int i = 0; i < 3; ++i) 
          pij_d[i] = delta_r*delta_r*duidn_hat[i]*duidn_hat[i] - duidn[i]*duidn[i];
        pij_od[0]  = delta_r*delta_r*duidn_hat[1]*duidn_hat[2] - duidn[1]*duidn[2];
        pij_od[1]  = delta_r*delta_r*duidn_hat[2]*duidn_hat[0] - duidn[2]*duidn[0];
        pij_od[2]  = delta_r*delta_r*duidn_hat[0]*duidn_hat[1] - duidn[0]*duidn[1];


        double lijpij_ = 0.0;
        double pijpij_ = 0.0;

        for (int i = 0; i < 3; ++i) { 

          lijpij_   += pij_d[i]*lij_d[icv][i];
          lijpij_   += 2.0*pij_od[i]*lij_od[icv][i];

          pijpij_   += pij_d[i]*pij_d[i];
          pijpij_   += 2.0*pij_od[i]*pij_od[i];

        }


        // lastly update the lijpij, pijpij contractions and then solve .. 

        const double T_relax = 0.0015; // hack .. at present.
        const double epsln   = 1.0/(T_relax/solver->dt + 1.0);

        lijpij[ibf] += epsln*( lijpij_ - lijpij[ibf]);
        pijpij[ibf] += epsln*( pijpij_ - pijpij[ibf]);


        const double cdel2 = lijpij[ibf] / (pijpij[ibf] + 1.0e-12);

        if ( cdel2 < 0.0) { 
          
          cdel_w[ibf] = 0.0;
          for (int i = 0; i < 3; ++i) 
            rhou_s[ibf][i] = 0.0;

        } else { 

         const double delta = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
         cdel_w[ibf]       = sqrt(cdel2);
         
         for (int i = 0; i < 3; ++i)  
           rhou_s[ibf][i] = cdel_w[ibf] * solver->rho[icv] * solver->u[icv][i] / (delta + cdel_w[ibf]);
         
        }

      }


      delete[] u_hat;
      delete[] lij_d;
      delete[] lij_od;
      delete[] tmp_d;
      delete[] tmp_od;

    }

  }

}

void WmSlipApprox::addBoundaryFlux(IdealGasRhs* rhs) const { 

  const double gogm1 = solver->gamma / (solver->gamma - 1.0);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

    const int icv       = zone_ptr->cvobf[ibf];
    const double p_wall = solver->p[icv];

    // we are proceeding with a fully incompressible approximation 
    // close to the wall; so we ignore any particular changes in the 
    // thermodynamic state, (p, rho, T) to the wall...

    const double mag_n    = MAG(zone_ptr->n_bf[ibf]);
    double unit_n[3];

    if ( mag_n > 0.0) { 

      for (int i = 0; i < 3; ++i) 
        unit_n[i] = zone_ptr->n_bf[ibf][i]/mag_n;

    } else { 

      for (int i = 0; i < 3; ++i) 
        unit_n[i] = 0.0;

    }

    const double rhoun    = DOT_PRODUCT(rhou_s[ibf], unit_n);
    const double rho_wall = solver->rho[icv];
    
    double Frho, Frhou[3], FrhoE;

    if ( rhoun > 0.0) { 

      const double H0       = gogm1 * solver->p[icv]/solver->rho[icv] + 
                              0.5*(DOT_PRODUCT(solver->u[icv],solver->u[icv]));

      Frho                  = rhoun*mag_n;
      for (int i = 0; i < 3; ++i) 
        Frhou[i]            = (p_wall*unit_n[i] + rhoun*solver->u[icv][i])*mag_n;
      FrhoE                 = rhoun*H0*mag_n;

      
    } else { 

      const double H_wall   = gogm1 * p_wall/rho_wall + 
                              0.5*(DOT_PRODUCT(rhou_s[ibf],rhou_s[ibf]))/(rho_wall*rho_wall);

      
      Frho                  = rhoun*mag_n;
      for (int i = 0; i < 3; ++i)
        Frhou[i]            = (p_wall*unit_n[i] + rhoun*rhou_s[ibf][i]/rho_wall)*mag_n;
      FrhoE                 = rhoun*H_wall*mag_n;

    }

    // add the resolved viscous closure; using only the laminar viscosity bc the dynamic 
    // procedure for the slip length above assumes that the slip velocities carry the full 
    // wall stress at the wall.. 
    
    const double visc_coeff = (solver->mu_lam[icv] )* 
                              zone_ptr->area_over_delta_bf[ibf];

    for (int i = 0; i < 3; ++i) 
      Frhou[i]             -= visc_coeff*(rhou_s[ibf][i]/rho_wall - solver->u[icv][i]); 


    // adiabatic closure -- neglecting the thermal gradient close to the wall.  there is 
    // some heat flux in the above momentum eqn that should be additional closed.  

    rhs[icv].rho -= Frho;
    for (int i = 0; i < 3; ++i) 
      rhs[icv].rhou[i] -= Frhou[i];
    rhs[icv].rhoE -= FrhoE;

  }

}

