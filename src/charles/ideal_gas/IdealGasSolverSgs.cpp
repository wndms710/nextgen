
#include "IdealGasSolver.hpp"

void IdealGasSolver::limitAsgs() { 

  if ( (checkParam("LIMIT_ASGS")) && (dt > 0.0) ) {

    // place a ceiling on the acoustic sgs based s.t. it won't incur a viscous
    // cfl stability limit for the explicitly stepped equations ..

    // dt may be initially zero for the very first time step if you are not restarting,
    // so we will check that condition here ..

    const double invgm1 = 1.0/(gamma-1.0);
    int8 my_buf = 0;

    for (int icv = 0; icv < ncv; ++icv) {

      // note that a_sgs has been divided by sos2.. we have to unwrap that first.. 

      //const double dx2      = r_vv[icv]*r_vv[icv];
      const double cfl_visc = 1.0; // approx >2 for rk3, but we'll see a safety factor here ..
      const double inv_sos2 = 1.0/(sos[icv]*sos[icv]);
      const double usq      = DOT_PRODUCT(u[icv],u[icv]);
      const double fax      = invgm1 + 0.5*usq*inv_sos2;

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

void IdealGasSolver::completeSgsStuff() { 

  if ( sgs_model == "DSM_LOCAL") { 

    if ( !checkDataFlag("lijmij") || !checkDataFlag("mijmij")) { 

      if ( mpi_rank == 0)
        cout << " resetting lijmij, mijmij " << endl;

      double *lijmij, *mijmij;
      CtiRegister::CtiData * data = CtiRegister::getRegisteredCtiData("lijmij"); assert( data != NULL);
      lijmij         = data->getDNptr();

      data           = CtiRegister::getRegisteredCtiData("mijmij"); assert( data != NULL);
      mijmij         = data->getDNptr();

      for (int icv = 0; icv < ncv_g2; ++icv) {
        lijmij[icv] = 0.0;
        mijmij[icv] = 1.0e-12;
      }

    }

  }
}

void IdealGasSolver::computeSgsNone() {
  for (int icv = 0; icv < ncv_g2; ++icv) {
    mu_sgs[icv]  = 0.0;
    loc_sgs[icv] = 0.0;
    a_sgs[icv]   = 0.0;
  }
}

void IdealGasSolver::computeSgsDsm() {

  using namespace CtiRegister;

  const double a_sgs_coeff  = getDoubleParam("A_SGS_COEFF",0.05); 
  const double a_sgs_coeff2 = getDoubleParam("A_SGS_COEFF2",0.0015); 
  const double Pr_t          = getDoubleParam("PR_T", 0.9);


  // we are passing this as a FlowSolver object so that the 
  // routine knows how to update and filter cv data ... 

  double *lijmij, *mijmij, *c_smag, *c_clip ;
  {

    CtiData * data = getRegisteredCtiData("lijmij"); assert( data != NULL);
    lijmij         = data->getDNptr();

    data           = getRegisteredCtiData("mijmij"); assert( data != NULL);
    mijmij         = data->getDNptr();

    data           = getRegisteredCtiData("c_smag"); assert( data != NULL);
    c_smag         = data->getDNptr();

    data           = getRegisteredCtiData("c_clip"); assert( data != NULL);
    c_clip         = data->getDNptr();

  }


  DynamicSmagorinsky::computeLijMij(lijmij, mijmij,
      rho,u,dudx,this);

  // and finally complete the calculation ... 

  int8 my_clip_count = 0; 

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

    double cdel2 = -0.5*lijmij[icv]/(mijmij[icv]+1.0e-12);
    c_clip[icv]  = 0.0;

    if ( cdel2 < 0.0) {

      if ( cdel2*smag < mu_ref) {
        ++my_clip_count;
        c_clip[icv] = 1.0;
      }

      cdel2 = 0.0;

    }

    const double dx2   = pow(vol_cv[icv],2.0/3.0);
    c_smag[icv]        = sqrt(cdel2/dx2);

    // ceiling on the smagorinsky coefficient ..

    if ( c_smag[icv] > 0.5) { 

      c_smag[icv] = 0.5;
      cdel2       = c_smag[icv]*dx2;

    }

    mu_sgs[icv]        = 1.0*cdel2*rho[icv]*smag; // recall smag already has the density factor

    // acoustic dissipation (to treat unresolved density fluctuations; see notes)
    // we are going to use the time scale compute from the dynamic coefficient normalized 
    // by the classical smagorinsky constant (0.19 a la lilly)

    const double inv_sos2  = 1.0/(sos[icv]*sos[icv]); 
    a_sgs[icv]             = a_sgs_coeff*cdel2*smag/(0.19*0.19*dx2); // using inv time scale from the dynamic smag 

    const double dil_      = dudx[icv][0][0] + dudx[icv][1][1] + dudx[icv][2][2];
    const double tmp_      = abs(dil_)*r_vv[icv]/sos[icv];
    if ( tmp_ > 0.01) {
      const double Masq    = (DOT_PRODUCT(u[icv],u[icv]))*inv_sos2;
      a_sgs[icv]          += a_sgs_coeff2*sa_to_vol_ratio[icv]*fmin(Masq,10.0)*sos[icv];
    }

    a_sgs[icv]            *= inv_sos2; // for use with the new flux.. 

  }

  limitAsgs();

  updateCv2ClassStart<SgsComm,double>(sgs_comm_container); // packed mu_sgs, a_sgs communication..

  const double invPr_t = 1.0/Pr_t;

  for (int icv = 0; icv < ncv; ++icv)
    loc_sgs[icv] = mu_sgs[icv] * invPr_t;

  // possible to delay the completion of these calls to later
  // in the material property calculation to hide some more latency..

  updateCv2ClassFinish<SgsComm,double>(sgs_comm_container);

  for (int icv = ncv; icv < ncv_g2; ++icv)
    loc_sgs[icv] = mu_sgs[icv] * invPr_t;

}

void IdealGasSolver::computeSgsVreman() {

  const double vreman_coeff   = getDoubleParam("VREMAN_COEFF", 0.07);
  const double a_sgs_coeff    = getDoubleParam("A_SGS_COEFF",0.05); 
  const double a_sgs_coeff2   = getDoubleParam("A_SGS_COEFF2",0.0015);
  const double locp_dil_coeff = getDoubleParam("LAMBDA_DIL_COEFF",0.0);
  const double Pr_t           = getDoubleParam("PR_T", 0.9);
  const double invPr_t        = 1.0/Pr_t;

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
    const double alpha32   = dudx[icv][2][1];
    const double alpha32_2 = alpha32*alpha32;

    const double beta11  = dx2*(alpha11_2+alpha21_2+alpha31_2);
    const double beta22  = dx2*(alpha12_2+alpha22_2+alpha32_2);
    const double beta33  = dx2*(alpha13_2+alpha23_2+alpha33_2);

    const double beta12  = dx2*(alpha11*alpha12+alpha21*alpha22+alpha31*alpha32);
    const double beta13  = dx2*(alpha11*alpha13+alpha21*alpha23+alpha31*alpha33);
    const double beta23  = dx2*(alpha12*alpha13+alpha22*alpha23+alpha32*alpha33);

    double B       = (beta11*beta22-beta12*beta12)+(beta11*beta33-beta13*beta13)+(beta22*beta33-beta23*beta23);
    B              = (B + abs(B))*0.5;

    const double alal    =
      ((alpha11_2+alpha22_2)+(alpha33_2+alpha12_2))+
      ((alpha13_2+alpha23_2)+(alpha21_2+alpha31_2))+alpha32_2;


    const double s_mag = sqrt(B/(alal+1.0E-20)); // includes lengthscale squared too

    // mu_sgs...
    mu_sgs[icv] = rho[icv]*vreman_coeff*s_mag;

    // acoustic dissipation (to treat unresolved density fluctuations; see notes)
    const double inv_sos2  = 1.0/(sos[icv]*sos[icv]); 
    a_sgs[icv]             = a_sgs_coeff*s_mag/dx2; // using inv time scale from the vreman kernel... 

    const double dil_      = dudx[icv][0][0] + dudx[icv][1][1] + dudx[icv][2][2];
    const double tmp_      = abs(dil_)*dx/sos[icv];
    if ( tmp_ > 0.01) {
      const double Masq    = (DOT_PRODUCT(u[icv],u[icv]))*inv_sos2;
      a_sgs[icv]          += a_sgs_coeff2*sa_to_vol_ratio[icv]*fmin(Masq,10.0)*sos[icv];
    }

    a_sgs[icv]            *= inv_sos2; // for use with the new flux.. 

    // dilatation corrections to the turbulent thermal conductivity .. 

    loc_sgs[icv]           = locp_dil_coeff*rho[icv]*abs(dil_)*dx2;
  }

  limitAsgs();

  updateCv2ClassStart<SgsComm,double>(sgs_comm_container); // packed mu_sgs, a_sgs communication..

  for (int icv = 0; icv < ncv; ++icv)
    loc_sgs[icv] += mu_sgs[icv] * invPr_t;

  // possible to delay the completion of these calls to later
  // in the material property calculation to hide some more latency..

  updateCv2ClassFinish<SgsComm,double>(sgs_comm_container);

  for (int icv = ncv; icv < ncv_g2; ++icv)
    loc_sgs[icv] += mu_sgs[icv] * invPr_t;

}

void IdealGasSolver::computeSgsSigma() {

  const double sigma_coeff  = getDoubleParam("SIGMA_COEFF", 1.5);
  const double a_sgs_coeff  = getDoubleParam("A_SGS_COEFF",0.05); 
  const double a_sgs_coeff2 = getDoubleParam("A_SGS_COEFF2",0.0015); 
  const double Pr_t         = getDoubleParam("PR_T", 0.9);
  const double invPr_t      = 1.0/Pr_t;

  for (int icv = 0; icv < ncv; ++icv) {

    const double dx2 = pow( vol_cv[icv], 2.0/3.0 ); // XXX precompute..

    double A_tmp[9];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j)
        A_tmp[(j)*3+i] = dudx[icv][i][j]; // store in col major format..
    }

    double U[9], V[9], sigma[3], max_err;
    calcSvd33(U,V,sigma,A_tmp,max_err);

    //
    // sort the sigma to be descending order ..
    // there's only 3 values, so we can bubble sort pretty easily ..
    //

    if ( sigma[0] < sigma[1])
      swap(sigma[0],sigma[1]);

    if ( sigma[1] < sigma[2])
      swap(sigma[1], sigma[2]);

    if ( sigma[0] < sigma[1])
      swap(sigma[0],sigma[1]);

    assert( (sigma[0] >= sigma[1]) && (sigma[0] >= sigma[2]));
    assert(  sigma[1] >= sigma[2]);

    double D_sigma = 0.0;
    if ( sigma[0] > 1.0e-12) { // protect against zero gradient case ..

      D_sigma = sigma[2]*(sigma[0] - sigma[1])*(sigma[1] - sigma[2])/(sigma[0]*sigma[0]);

    }

    mu_sgs[icv]            = sigma_coeff*sigma_coeff*dx2*rho[icv]*D_sigma;

    // acoustic dissipation (to treat unresolved density fluctuations; see notes)

    const double inv_sos2  = 1.0/(sos[icv]*sos[icv]);
    a_sgs[icv]             = a_sgs_coeff*D_sigma;

    const double dil_      = dudx[icv][0][0] + dudx[icv][1][1] + dudx[icv][2][2];
    const double tmp_      = abs(dil_)*sqrt(dx2)/sos[icv];
    if ( tmp_ > 0.01) {
      const double Masq    = (DOT_PRODUCT(u[icv],u[icv]))*inv_sos2;
      a_sgs[icv]          += a_sgs_coeff2*sa_to_vol_ratio[icv]*fmin(Masq,10.0)*sos[icv];
    }

    a_sgs[icv]            *= inv_sos2; // for use with the new flux..
  }

  limitAsgs();

  updateCv2ClassStart<SgsComm,double>(sgs_comm_container); // packed mu_sgs, a_sgs communication..

  for (int icv = 0; icv < ncv; ++icv)
    loc_sgs[icv] = mu_sgs[icv] * invPr_t;

  // possible to delay the completion of these calls to later
  // in the material property calculation to hide some more latency..

  updateCv2ClassFinish<SgsComm,double>(sgs_comm_container);

  for (int icv = ncv; icv < ncv_g2; ++icv)
    loc_sgs[icv] = mu_sgs[icv] * invPr_t;

}

