#include "NonpremixedSolver.hpp"

// turn on a_sgs limiting by default
/*
void NonpremixedSolver::limitAsgs() { 

  if ( (checkParam("LIMIT_ASGS")) && (dt > 0.0) ) {

    // place a ceiling on the acoustic sgs based s.t. it won't incur a viscous
    // cfl stability limit for the explicitly stepped equations ..

    // dt may be initially zero for the very first time step if you are not restarting,
    // so we will check that condition here ..

    int8 my_buf = 0;

    for (int icv = 0; icv < ncv; ++icv) {

      // note that a_sgs has been divided by sos2.. we have to unwrap that first.. 

      //const double dx2      = r_vv[icv]*r_vv[icv];
      const double cfl_visc = 1.0; // approx >2 for rk3, but we'll see a safety factor here ..
      const double inv_sos2 = 1.0/(sos[icv]*sos[icv]);
      const double usq      = DOT_PRODUCT(u[icv],u[icv]);
      const double fax      = 1.0/(gamma[icv]-1.0) + 0.5*usq*inv_sos2;

      //const double ceil_    = cfl_visc*dx2/(dt*fax*4.0)*inv_sos2;
      const double ceil_    = cfl_visc/(dt*fax*4.0)*inv_sos2;

      if ( a_sgs[icv] > ceil_) { 

        a_sgs[icv] = ceil_;
        ++my_buf;

      }

    }

    int8 buf;
    MPI_Reduce(&my_buf,&buf,1,MPI_INT8,MPI_SUM,0,mpi_comm);

    if ( (mpi_rank == 0) && (step%check_interval ==0))
      cout << " > a_sgs ceiling hit in : " << buf << " ; of " << getNcvGlobal() << " cells [ "
        << 100.0*double(buf)/double(getNcvGlobal()) << "  % ]" << endl;

  }

}
*/

void NonpremixedSolver::computeSgsNone() {
  for (int icv = 0; icv < ncv_g2; ++icv) {
    mu_sgs[icv]  = 0.0;
    loc_sgs[icv] = 0.0;
    a_sgs[icv]   = 0.0;
  }
}

void NonpremixedSolver::computeSgsVreman() {

  const double vreman_coeff   = getDoubleParam("VREMAN_COEFF", 0.07);
  const double a_sgs_coeff    = getDoubleParam("A_SGS_COEFF",0.05); 
  const double a_sgs_coeff2   = getDoubleParam("A_SGS_COEFF2",0.0015);
  const double locp_dil_coeff = getDoubleParam("LAMBDA_DIL_COEFF",0.0);
  const double Pr_t           = getDoubleParam("PR_T", 0.9);
  const double invPr_t        = 1.0/Pr_t;

  // we are also augmenting the subgrid coefficients based on non-smoothness
  // arising from the tabulated lookups.  this can arise due to poor sampling 
  // of the table or gibbs oscillations in the scalars leading to lookups that 
  // cross realizable boundaries in the chemtable.

  const double c_table   = getDoubleParam("C_TABLE", 5.0);
  const double c_table_p = getDoubleParam("C_TABLE_P", 5.0);

  double *dT        = new double[ncv_g2];
  for (int icv = 0; icv < ncv_g2; ++icv) 
    dT[icv] = T[icv] - T_p[icv];
  double (*dTdx)[3] = new double[ncv][3];
  StaticSolver::calcCvGrad(dTdx,dT);
  delete[] dT;

  for (int icv = 0; icv < ncv; ++icv) { 

    const double dx        = pow(vol_cv[icv], 1.0/3.0);
    const double dx2       = dx*dx;

    const double alpha11   = dudx[icv][0][0];
    const double alpha11_2 = alpha11*alpha11; 
    const double alpha22   = dudx[icv][1][1];
    const double alpha22_2 = alpha22*alpha22;
    const double alpha33   = dudx[icv][2][2];
    const double alpha33_2 = alpha33*alpha33;

    const double alpha12   = dudx[icv][0][1];
    const double alpha12_2 = alpha12*alpha12;
    const double alpha13   = dudx[icv][0][2];
    const double alpha13_2 = alpha13*alpha13;
    const double alpha23   = dudx[icv][1][2];
    const double alpha23_2 = alpha23*alpha23;

    const double alpha21   = dudx[icv][1][0];
    const double alpha21_2 = alpha21*alpha21;
    const double alpha31   = dudx[icv][2][0];
    const double alpha31_2 = alpha31*alpha31;
    const double alpha32   = dudx[icv][2][1]; // units: 1/T
    const double alpha32_2 = alpha32*alpha32; // units: 1/T2

    const double beta11    = dx2*(alpha11_2+alpha21_2+alpha31_2);
    const double beta22    = dx2*(alpha12_2+alpha22_2+alpha32_2);
    const double beta33    = dx2*(alpha13_2+alpha23_2+alpha33_2); // units: L2/T2

    const double beta12    = dx2*(alpha11*alpha12+alpha21*alpha22+alpha31*alpha32);
    const double beta13    = dx2*(alpha11*alpha13+alpha21*alpha23+alpha31*alpha33);
    const double beta23    = dx2*(alpha12*alpha13+alpha22*alpha23+alpha32*alpha33);

    double B               = (beta11*beta22-beta12*beta12)+
                             (beta11*beta33-beta13*beta13)+
                             (beta22*beta33-beta23*beta23);
    B                      = (B + abs(B))*0.5; // units: L4/T4

    const double alal      = ((alpha11_2+alpha22_2)+(alpha33_2+alpha12_2))+
                             ((alpha13_2+alpha23_2)+(alpha21_2+alpha31_2))+alpha32_2; // units: 1/T2

    const double s_mag     = sqrt(B/(alal+1.0E-20)); // includes lengthscale squared ... units: sqrt(L4/T4*T2) = L2/T

    // mu_sgs...
    mu_sgs[icv]            = rho[icv]*vreman_coeff*s_mag; // units: M/L3*L2/T = M/(L*T) ... correct units for viscosity: kg/(m*s)

    // acoustic dissipation (to treat unresolved density fluctuations; see notes)
    const double inv_sos2  = 1.0/(sos[icv]*sos[icv]); 
    a_sgs[icv]             = a_sgs_coeff*s_mag/dx2; // using inv time scale from the vreman kernel ... units: 1/T

    const double dil_      = dudx[icv][0][0] + dudx[icv][1][1] + dudx[icv][2][2];
    /*
    // skip this part for now: developed for high-Ma flows in the non-reacting sovler, it
    // should be investigated carefully before incorporating in reacting flows where dil can be
    // present due to heat release...
    const double tmp_      = abs(dil_)*dx/sos[icv];
    if ( tmp_ > 0.01) {
      const double Masq    = (DOT_PRODUCT(u[icv],u[icv]))*inv_sos2;
      a_sgs[icv]          += a_sgs_coeff2*sa_to_vol_ratio[icv]*fmin(Masq,10.0)*sos[icv]; // units: L2/L3*L/T = 1/T
    }
    */

    a_sgs[icv]            *= inv_sos2; // for use with the new flux ... units: 1/T*T2/L2 = T/L2

    // dilatation corrections to the turbulent thermal conductivity .. 

    loc_sgs[icv]           = locp_dil_coeff*rho[icv]*abs(dil_)*dx2;

    const double invgm1    = 1.0/(gamma[icv] - 1.0);
    const double tt        = invgm1*dx*MAG(dTdx[icv])/T[icv]; // note: this is non-dimensional

    loc_sgs[icv]          += mu_sgs[icv] * invPr_t;

    // augment the scalar sgs model based on the non-smoothness of the lookup...

    double dudx_mag2       = 0.0;
    FOR_I3 FOR_J3 dudx_mag2 += 2.0*dudx[icv][i][j]*dudx[icv][i][j]; 
    loc_sgs[icv]          += c_table*dx2*rho[icv]*sqrt(dudx_mag2)*tt;

    // also augment the acoustic subgrid closure ..

    a_sgs[icv]            += c_table_p*tt/sos[icv]/dx; // units: 1/(L/T)/L = T/L2

    // and finally limit a_sgs to control stiffness...

    a_sgs[icv] = min(a_sgs[icv],0.5/sos[icv]/dx); // units: 1/(L/T)/L = T/L2

  }

  // limiting now handled above...
  //limitAsgs();

  updateCv2Class<SgsComm,double>(sgs_comm_container); // packed mu_sgs, a_sgs communication..

  delete[] dTdx;
}
