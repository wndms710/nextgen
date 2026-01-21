#include "PremixedSolver.hpp"
#include "RkWeights.hpp"
#include "wm/AlgebraicWM.hpp"
#include "SpongeCommon.hpp"
#include "it/InflowTurbulenceP.hpp"

using namespace RkWeights;

void SlipWallPBc::addBoundaryFlux(PremixedRhs* rhs) const {

  assert( solver != NULL);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv = zone_ptr->cvobf[ibf];
    for (int i =0; i < 3; ++i)
      rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];
  }
}

MsTableP::MsTableP(Param *param, BfZone *p, PremixedSolver *s) : PremixedBc(p,s) {
  p->lb_cost = 17; // account for extra work associated with this bc
  mdot_bc = param->getDouble(1);
  Z_bc    = param->getDouble(2);
  C_bc    = param->getDouble(3);
  assert(mf == NULL);
  rho_bf = NULL;
  mflux_bf = NULL;
  if (mpi_rank == 0) {
    cout << "MsTableP() mdot_bc: " << mdot_bc << " Z_bc: " << Z_bc << " C_bc: " << C_bc << endl;
  }
}

void MsTableP::initData() {
  assert(mf == NULL); mf = new double[zone_ptr->nbf];
  assert(rho_bf == NULL); rho_bf = new double[zone_ptr->nbf];
  assert(mflux_bf == NULL); mflux_bf = new double[zone_ptr->nbf];
  solver->chemtable->lookup(&R_bc,"R",&T_bc,"T",&e_bc,"e",&gamma_bc,"gamma",&Z_bc,&C_bc,1);
  p_bc = solver->chemtable->pressure;
  rho_bc = p_bc/(R_bc*T_bc);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    rho_bf[ibf] = rho_bc; // TODO: should read this in from file?
    mflux_bf[ibf] = -mdot_bc/zone_ptr->area_global; // the target mass flux
  }
  if (mpi_rank == 0) {
    double un_bc = mdot_bc/(rho_bc*zone_ptr->area_global);
    cout << "MsTableP::initData() p_bc: " << p_bc << " rho_bc: " << rho_bc << " un_bc: " << un_bc << " gamma_bc: " << gamma_bc << " Ma: " << un_bc/sqrt(gamma_bc*p_bc/rho_bc) << endl;
  }
}

MsTableP::~MsTableP() {
  DELETE(mf);
  DELETE(rho_bf);
  DELETE(mflux_bf);
}

inline double calcRiemannFluxRhoOnly(const double unit_n[3],
				     const PremixedState& cv0, const PremixedState& cv1,
				     const double gammaL, const double gammaR) {

  // This special excerpt from the HLLC Riemann flux solver just computes the rho-flux
  // for the constant mass flux boundary conditions. ..
  
  double unL  = DOT_PRODUCT(cv0.u, unit_n);
  double uLuL = DOT_PRODUCT(cv0.u, cv0.u);
  double gamma_p_over_rhoL = gammaL*cv0.p*cv0.sp_vol;
  double cL   = sqrt(gamma_p_over_rhoL);
  //double HL   = cv0.h + 0.5*uLuL; // capitalization means ke included
  double rhoL = 1.0/cv0.sp_vol;

  double unR  = DOT_PRODUCT(cv1.u, unit_n);
  double uRuR = DOT_PRODUCT(cv1.u, cv1.u);
  double gamma_p_over_rhoR = gammaR*cv1.p*cv1.sp_vol;
  double cR   = sqrt(gamma_p_over_rhoR);
  //double HR   = cv1.h + 0.5*uRuR; // include ke here
  double rhoR = 1.0/cv1.sp_vol;

  // Roe's averaging
  double Rrho = sqrt(rhoR/rhoL);
  double tmp  = 1.0/(1.0+Rrho);
  double uRoe[3];
  FOR_I3 uRoe[i] = tmp*(cv0.u[i] + cv1.u[i]*Rrho);
  double unRoe   = DOT_PRODUCT(uRoe, unit_n);
  //double HRoe    = tmp*(HL + HR*Rrho);

  //  double cRoe  = sqrt((gammaL-1.0)*(hRoe- 0.5*DOT_PRODUCT(uRoe, uRoe)));
  double gamPdivRho = tmp*((gamma_p_over_rhoL+0.5*(gammaL-1.0)*uLuL) + (gamma_p_over_rhoR+0.5*(gammaR-1.0)*uRuR)*Rrho);
  double cRoe       = sqrt(gamPdivRho-(0.5*(gammaL+gammaR)-1.0)*0.5*DOT_PRODUCT(uRoe,uRoe));

  // speed of sound at L and R
  double sL = min(unRoe-cRoe, unL-cL);
  double sR = max(unRoe+cRoe, unR+cR);

  // speed of contact surface
  double sM = (cv0.p-cv1.p-rhoL*unL*(sL-unL)+rhoR*unR*(sR-unR))/(rhoR*(sR-unR)-rhoL*(sL-unL));

  // pressure at right and left (pR=pL) side of contact surface
  //double pStar = rhoR*(unR-sR)*(unR-sM)+cv1.p;

  if (sM >= 0.0) {
    if (sL > 0.0) {
      //cout << "1" << endl;
      return rhoL*unL;
      //FOR_I3 flux.rhou[i] = rhoL*cv0.u[i]*unL + cv0.p*unit_n[i];
      //flux.rhoE = rhoL*HL*unL;
    }
    else {
      //cout << "2" << endl;
      double invSLmSs = 1.0/(sL-sM);
      double sLmuL = sL-unL;
      double rhoSL = rhoL*sLmuL*invSLmSs;
      //double rhouSL[3];
      //FOR_I3 rhouSL[i] = (rhoL*cv0.u[i]*sLmuL+(pStar-cv0.p)*unit_n[i])*invSLmSs;
      //double eSL = (sLmuL*(rhoL*HL-cv0.p)-cv0.p*unL+pStar*sM)*invSLmSs;
      return rhoSL*sM;
      //FOR_I3 flux.rhou[i] = rhouSL[i]*sM + pStar*unit_n[i];
      //flux.rhoE = (eSL+pStar)*sM;
    }
  }
  else {
    if (sR >= 0.0) {
      //cout << "3" << endl;
      double invSRmSs = 1.0/(sR-sM);
      double sRmuR = sR-unR;
      double rhoSR = rhoR*sRmuR*invSRmSs;
      //double rhouSR[3];
      //FOR_I3 rhouSR[i] = (rhoR*cv1.u[i]*sRmuR+(pStar-cv1.p)*unit_n[i])*invSRmSs;
      //double eSR = (sRmuR*(rhoR*HR-cv1.p)-cv1.p*unR+pStar*sM)*invSRmSs;
      return rhoSR*sM;
      //FOR_I3 flux.rhou[i] = rhouSR[i]*sM + pStar*unit_n[i];
      //flux.rhoE = (eSR+pStar)*sM;
    }
    else {
      //cout << "4" << endl;
      return rhoR*unR;
      //FOR_I3 flux.rhou[i] = rhoR*cv1.u[i]*unR + cv1.p*unit_n[i];
      //flux.rhoE = rhoR*HR*unR;
    }
  }

  // leave flux.rho unmultiplied by the area...
  /*
    flux.rho *= area;
    FOR_I3 flux.rhou[i] *= area;
    flux.rhoE *= area;
    flux.rhoZ  = 0.5*( flux.rho*( cv0.Z + cv1.Z ) - abs(flux.rho)*( cv1.Z - cv0.Z) );
    flux.rhoC  = 0.5*( flux.rho*( cv0.C + cv1.C ) - abs(flux.rho)*( cv1.C - cv0.C) );
  */

}

void MsTableP::addBoundaryFlux(PremixedRhs *rhs) const {
  // the idea here: isetropically adjust rho_bc_star and p_bc_star
  // until we get the correct mass flux...
  PremixedState bf_state;
  bf_state.Z = Z_bc;
  bf_state.C = C_bc;
  double unit_n[3];
  PremixedRhs flux;
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const double mf_bc = mflux_bf[ibf]; // this is the target (allow it to vary spatially eventually)
    const int icv = zone_ptr->cvobf[ibf];
    double area = MAG(zone_ptr->n_bf[ibf]);
    if (area > 0.0) {
      FOR_I3 unit_n[i] = zone_ptr->n_bf[ibf][i]/area;
      // use the secant method to converge the mass flux to the target
      // with isentropic perturbations of rho and p. Note to self:
      // I looked at doing this with newton and returning a gradient from
      // the flux function, but it is actually quite a bit of math. This 
      // approach converges in a few iterations, and the RhoOnly version of
      // the flux is pretty light...
      // start with a perterbed rho value...
      double rho1 = rho_bf[ibf]*1.001; // finite, small perturbation
      bf_state.sp_vol = 1.0/rho1;
      FOR_I3 bf_state.u[i] = mf_bc/rho1*unit_n[i];
      bf_state.p      = p_bc*pow(rho1/rho_bc,gamma_bc);
      // note that this partial Riemann solver takes the unit normal and
      // returns the mass flux...
      double mf1 = calcRiemannFluxRhoOnly(unit_n,solver->cv_light[icv],bf_state,solver->gamma[icv],gamma_bc);
      // then the current rho value. In some cases, this will skip the while
      // loop because it is alrady good enough, and we need uX,pX to be 
      // consistent...
      double rho0 = rho_bf[ibf]; // use the stored rho value here
      bf_state.sp_vol = 1.0/rho0;
      FOR_I3 bf_state.u[i] = mf_bc/rho0*unit_n[i];
      bf_state.p      = p_bc*pow(rho0/rho_bc,gamma_bc);
      double mf0 = calcRiemannFluxRhoOnly(unit_n,solver->cv_light[icv],bf_state,solver->gamma[icv],gamma_bc);
      // the mf corresponding to rho_bf[ibf] is mf0...
      int count = 0;
      double mf_test = mf0; 
      while (fabs((mf_test-mf_bc)/mf_bc) >= 1.0E-12) {
	// make sure we don't stall...
	++count;
	assert(count < 50);
	// calc and limit new rho_bf[ibf]...
	rho_bf[ibf] = ((mf_bc-mf0)*rho1 + (mf1-mf_bc)*rho0)/(mf1-mf0);
	rho_bf[ibf] = max(min(rho0,rho1)*0.9,rho_bf[ibf]);
	rho_bf[ibf] = min(max(rho0,rho1)*1.1,rho_bf[ibf]);
	bf_state.sp_vol = 1.0/rho_bf[ibf];
	FOR_I3 bf_state.u[i] = mf_bc/rho_bf[ibf]*unit_n[i];
	bf_state.p = p_bc*pow(rho_bf[ibf]/rho_bc,gamma_bc);
	mf_test = calcRiemannFluxRhoOnly(unit_n,solver->cv_light[icv],bf_state,solver->gamma[icv],gamma_bc);
	// replace the mf0 or mf1 that is further from the target...
	if (fabs(mf1-mf_bc) > fabs(mf0-mf_bc)) {
	  rho1 = rho_bf[ibf];
	  mf1 = mf_test;
	}
	else {
	  rho0 = rho_bf[ibf];
	  mf0 = mf_test;
	}
      }
      // converged! now compute full flux...
      bf_state.h = e_bc + bf_state.p/rho_bf[ibf];
      calcRiemannFlux(flux,zone_ptr->n_bf[ibf],solver->cv_light[icv],bf_state,solver->gamma[icv],gamma_bc);
      // check consistency or RhoOnly flux and full RiemannFlux...
      if (fabs(flux.rho-mf_bc*area) > 1.001E-12*fabs(mf_bc*area)) {
	cout << "WARNING: MS_AT_TABLE bc not converging: " << count << " " << flux.rho << " " << mf_bc << " " << area << " " << flux.rho/area-mf_bc << endl;
      }
      mf[ibf]       = flux.rho;
      rhs[icv].rho -= flux.rho;
      FOR_I3 rhs[icv].rhou[i] -= flux.rhou[i];
      rhs[icv].rhoE -= flux.rhoE;
      rhs[icv].rhoZ -= flux.rhoZ;
      rhs[icv].rhoC -= flux.rhoC;
    }
  }
}

void MsTableP::query(const string& param_str) const {
  double my_buf[3] = { 0.0, 0.0, 0.0 };
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv = zone_ptr->cvobf[ibf];
    my_buf[0] -= mf[ibf]; // inward mf is negative
    double area = MAG(zone_ptr->n_bf[ibf]);
    my_buf[1] += area;
    my_buf[2] += area*solver->p[icv];
  }
  double buf[3];
  MPI_Reduce(my_buf,buf,3,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  if (mpi_rank == 0) {
    cout << " > QUERY_BC " << zone_ptr->getName() << ": time,mdot,avg(p): " << solver->time << " " << buf[0] << " " << buf[2]/buf[1] << endl;
  }
}

void WallAdiabaticP::addBoundaryFlux(PremixedRhs* rhs) const {

  assert( solver != NULL);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv = zone_ptr->cvobf[ibf];
    const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
    for (int i =0; i < 3; ++i)
      rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i] + visc_coeff*(solver->u[icv][i] - u_bc[i]);
  }
}

void WallAdiabaticP::query(const string& param_str) const {

  double my_buf[4] = {0.0,0.0,0.0,0.0};
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

      const double y1         = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
      const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
      double du_mag = 0.0;
      for (int i =0; i < 3; ++i)
        du_mag  += (solver->u[icv][i] - u_bc[i])*(solver->u[icv][i] - u_bc[i]);

      const double tau_wall   = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
      const double y_plus     = y1*sqrt(tau_wall*solver->rho[icv])/solver->mu_lam[icv];

      my_buf[0] += zone_ptr->area_bf[ibf];
      my_buf[1] += zone_ptr->area_bf[ibf]*tau_wall;
      my_buf[2] += zone_ptr->area_bf[ibf]*y_plus;
      my_buf[3] += zone_ptr->area_bf[ibf]*solver->T[icv];

    }
  }

  double buf[4];
  MPI_Reduce(my_buf,buf,4,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if ( mpi_rank == 0 )
    cout << " > " << zone_ptr->getName() << ": time, int(tau_w)dA, avg(y_plus), avg(T): " <<
      solver->time << " " << buf[1] << " " << buf[2]/buf[0] << " " << buf[3]/buf[0] << endl;
}

CtiRegister::CtiDataError WallAdiabaticP::funcEvalCtiData(CtiRegister::CtiData& v,
                                                          const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if ( zone_ptr->isLeadingStrMatched(name) ) {

    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix.  we have to check for the actual data now..

    const string tau_wall_str = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str  = zone_ptr->getName() + ":" + "y_plus";

    // all of the following functions do not take any arguments -- check first

    if ( name == tau_wall_str) {

      double * tmp = zone_ptr->createBfD1Data(v);

      // if stats have been requested on this quantity, this function will be called
      // but its possible that the mu_lam, rho necessary have not been populated yet.
      // check the data flags for these values -- revisit if this impacts performance
      // during the solution run-time.

      if (b_eval_func) {
        if (args.size() == 0) {
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            const int icv           = zone_ptr->cvobf[ibf];
            const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
            double du_mag = 0.0;
            for (int i =0; i < 3; ++i)
              du_mag  += (solver->u[icv][i] - u_bc[i])*(solver->u[icv][i] - u_bc[i]);

            tmp[ibf]                = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
          }
        }
        else if (args.size() == 1) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          double dir[3] = {0.0,0.0,0.0};

          if (arg->getType() == D3_DATA) {
            FOR_I3 dir[i] = arg->d3(i);
            NORMALIZE(dir);
          }
          else if (arg->getType() == D_DATA) {
            const int index = int(arg->d());
            if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;
            dir[index] = 1.0;
          }
          else return CTI_DATA_NOT_VALID;

          // loop faces and compute
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            const int icv           = zone_ptr->cvobf[ibf];
            const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];

            double local_u[3] = DIFF(solver->u[icv],u_bc);
            tmp[ibf]          = visc_coeff*DOT_PRODUCT(dir,local_u)/zone_ptr->area_bf[ibf];
          }
        }
        else return CTI_DATA_ARG_COUNT;
      }

      return CTI_DATA_OK;

    }
    else if ( name == yp_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          const int icv           = zone_ptr->cvobf[ibf];
          const double y1         = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
          const double visc_coeff = solver->mu_lam[icv]/y1;

          double du_mag = 0.0;
          for (int i =0; i < 3; ++i)
            du_mag  += (solver->u[icv][i] - u_bc[i])*(solver->u[icv][i] - u_bc[i]);

          const double tau_wall   = visc_coeff*sqrt(du_mag);
          tmp[ibf]                = y1* sqrt(tau_wall*solver->rho[icv])/solver->mu_lam[icv];

        }

      }

      return CTI_DATA_OK;

    }
  }

  return PremixedBc::funcEvalCtiData(v,name,args,b_eval_func);
}

void WmAlgAdiabaticP::initData() {

  assert( tau_wall == NULL); tau_wall = new double[zone_ptr->nbf];
  assert( u1    == NULL); u1    = new double[zone_ptr->nbf];

  //for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
  //  tau_wall[ibf] = 0.0;

}

void WmAlgAdiabaticP::initialHook() {

  if ( !zone_ptr->checkDataFlag("u1")) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      const int icv = zone_ptr->cvobf[ibf];
      u1[ibf] = MAG(solver->u[icv]);

      // set an initial value for tau_wall, this will constitute the
      // initial guess for the solution when it is finally constructed.

      if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

        const double y1 = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        tau_wall[ibf] = solver->mu_lam[icv] * MAG(solver->u[icv])/y1;

      }
      else {
        tau_wall[ibf] = 0.0;
      }
    }
  }
}

int WmAlgAdiabaticP::storeData(double * buf) const {

  // the wall model needs to store tau_wall because it is used as
  // an initial guess. No need to store u1. It is not used...

  if (buf) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      buf[ibf] = tau_wall[ibf];
    }
  }

  return zone_ptr->nbf;

}

int WmAlgAdiabaticP::restoreData(double * buf) {

  assert(buf);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    tau_wall[ibf] = buf[ibf];
  }

  return zone_ptr->nbf;

}

WmAlgAdiabaticP::~WmAlgAdiabaticP() {

  DELETE(u1);
  DELETE(tau_wall);

}

void WmAlgAdiabaticP::calcRhs(const double time, const int rk_stage) {

  // see notes left in IdealGasSolverBcs.cpp for the time
  // filtering approach

  if ( rk_stage == 1 ) {

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      const int icv = zone_ptr->cvobf[ibf];

      if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

        const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        const double mag_u = MAG(solver->u[icv]);
        //const double eps   = mag_u*solver->dt / (y1 + solver->dt*mag_u);

        // \phi_new = eps*u + (1.0 -eps)*\phi_old
        // \phi_new - \phi_old = eps*(u - \phi_old)

        //u1[ibf]   += eps* ( mag_u - u1[ibf]);
        u1[ibf]   = mag_u; // Re 2000 channel test showed time averaging to have minimal impact, so removed

        // use the existing tau_wall to set a guess for the u_tau

        const double u_tau = sqrt(tau_wall[ibf] / solver->rho[icv]);
        tau_wall[ibf] = AlgebraicWM::solve_tau(u1[ibf], y1, solver->rho[icv],
                                               solver->mu_lam[icv], u_tau);

      } else {

        //assert ( zone_ptr->area_bf[ibf] < 1.0e-16); // this should be a collapsed face.
        tau_wall[ibf] = 0.0;
        u1[ibf]    = MAG(solver->u[icv]);

      }
    }
  }
}

void WmAlgAdiabaticP::addBoundaryFlux(PremixedRhs* rhs) const {

  assert( solver != NULL);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    // tau_wall has already been computed at this point, but we'll assume
    // that the wall stress is instantaneously aligned with the velocity vector

    double u_mag = 0.0;
    for (int i =0; i < 3; ++i)
      u_mag += solver->u[icv][i] * solver->u[icv][i];
    u_mag = sqrt( max(0.0, u_mag));

    if ( u_mag > 0.0) {

      for (int i =0; i < 3;++i) {
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];
        rhs[icv].rhou[i] -= tau_wall[ibf]*solver->u[icv][i]/u_mag*zone_ptr->area_bf[ibf];
      }

    } else {

      for (int i =0; i < 3; ++i)
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];

    }
  }
}

CtiRegister::CtiDataError WmAlgAdiabaticP::funcEvalCtiData(CtiRegister::CtiData& v,
                                                           const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if ( zone_ptr->isLeadingStrMatched(name) ) {

    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix.  we have to check for the actual data now..

    // for consistency, we are supporting a tau_wall(), y_plus()
    // function to have the same syntax as all of the other wall zones..

    const string tau_wall_str = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str  = zone_ptr->getName() + ":" + "y_plus";

    // all of the following functions do not take any arguments -- check first


    if ( name == tau_wall_str) {

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {
        if (args.size() == 0) {
          // use what was computed by wall model
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
            tmp[ibf] = tau_wall[ibf];
          }
        }
        else if (args.size() == 1) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          double dir[3] = {0.0,0.0,0.0};

          if (arg->getType() == D3_DATA) {
            FOR_I3 dir[i] = arg->d3(i);
            NORMALIZE(dir);
          }
          else if (arg->getType() == D_DATA) {
            const int index = int(arg->d());
            if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;
            dir[index] = 1.0;
          }
          else return CTI_DATA_NOT_VALID;

          // loop faces and compute
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
            const int icv      = zone_ptr->cvobf[ibf];
            const double mag_u = MAG(solver->u[icv]);

            double unit_u[3];
            if ( mag_u > 0.0) {
              for (int i = 0; i < 3; ++i) unit_u[i] = solver->u[icv][i]/mag_u;

              tmp[ibf] = tau_wall[ibf]*DOT_PRODUCT(dir,unit_u);
            }
            else {
              tmp[ibf] = 0.0;
            }
          }
        }
        else return CTI_DATA_ARG_COUNT;
      }

      return CTI_DATA_OK;

    } else if ( name == yp_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          const int icv      = zone_ptr->cvobf[ibf];
          const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
          tmp[ibf]           = y1*sqrt(solver->rho[icv]*tau_wall[ibf])/solver->mu_lam[icv];
        }

      }

      return CTI_DATA_OK;

    }
  }

  return PremixedBc::funcEvalCtiData(v,name,args,b_eval_func);
}

void WmAlgAdiabaticP::query(const string& param_str) const {

  // report an estimate of the local y^+ on the boundary and the integrated tau_wall

  double my_buf[4] = {0.0,0.0,0.0,0.0};
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

      const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
      const double u_tau = sqrt( tau_wall[ibf]/solver->rho[icv]);
      const double nu    = solver->mu_lam[icv] / solver->rho[icv];
      const double y_plus = y1*u_tau/nu;

      my_buf[0] += zone_ptr->area_bf[ibf];
      my_buf[1] += zone_ptr->area_bf[ibf]*tau_wall[ibf];
      my_buf[2] += zone_ptr->area_bf[ibf]*y_plus;
      my_buf[3] += zone_ptr->area_bf[ibf]*solver->T[icv];

    }
  }

  double buf[4];
  MPI_Reduce(my_buf,buf,4,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if ( mpi_rank == 0 ) {

    cout << " > " << zone_ptr->getName() << ": time, int(tau_wall)dA, avg(y_plus), avg(T): " <<
      solver->time << " " << buf[1] << " " << buf[2]/buf[0] << " " << buf[3]/buf[0] << endl;

  }
}

void WallIsothermalP::initData() {
  Param* param = getParam(getName());
  T_bc         = param->getDouble(1);
}

void WallIsothermalP::addBoundaryFlux(PremixedRhs* rhs) const {

  assert( solver != NULL);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv = zone_ptr->cvobf[ibf];

    // inherit the R, gamma from the interior cell to compute the cp
    const double cp = solver->R[icv]*solver->gamma[icv]/(solver->gamma[icv]-1.0);
    const double visc_coeff = solver->mu_lam[icv]  * zone_ptr->area_over_delta_bf[ibf];
    const double k_coeff    = solver->loc_lam[icv] * zone_ptr->area_over_delta_bf[ibf]*cp;
    for (int i =0; i < 3; ++i)
      rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i] + visc_coeff*(solver->u[icv][i] - u_bc[i]);
    // XXX if u_bc \neq 0, there should be a u \dot \tau term added here..
    rhs[icv].rhoE -= k_coeff*(solver->T[icv] - T_bc);
  }
}

void WallIsothermalP::query(const string& param_str) const {

  double my_buf[4] = {0.0,0.0,0.0,0.0};

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

      const double y1         = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
      const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
      double du_mag = 0.0;
      for (int i =0; i < 3; ++i)
        du_mag  += (solver->u[icv][i] - u_bc[i])*(solver->u[icv][i] - u_bc[i]);

      const double tau_wall   = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
      const double y_plus      = y1*sqrt(tau_wall*solver->rho[icv])/solver->mu_lam[icv];

      // inherit the R, gamma from the interior cell to compute the cp
      const double cp         = solver->R[icv]*solver->gamma[icv]/(solver->gamma[icv]-1.0);
      const double k_coeff    = solver->loc_lam[icv]*zone_ptr->area_over_delta_bf[ibf]*cp;
      const double q_wall     = k_coeff*abs(solver->T[icv] - T_bc)/zone_ptr->area_bf[ibf];

      my_buf[0]              += zone_ptr->area_bf[ibf];
      my_buf[1]              += zone_ptr->area_bf[ibf]*tau_wall;
      my_buf[2]              += zone_ptr->area_bf[ibf]*y_plus;
      my_buf[3]              += zone_ptr->area_bf[ibf]*q_wall;

    }
  }

  double buf[4];
  MPI_Reduce(my_buf,buf,4,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if ( mpi_rank == 0 ) {

    cout << " > " << zone_ptr->getName() << ": time, int(tau_wall)dA, int(q_wall)dA, avg(y_plus): "
         << solver->time << "     "
         << buf[1] << "     "
         << buf[3] << "     "
         << buf[2]/buf[0] << endl;

  }
}

CtiRegister::CtiDataError WallIsothermalP::funcEvalCtiData(CtiRegister::CtiData& v,
                                                           const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if ( zone_ptr->isLeadingStrMatched(name) ) {

    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix.  we have to check for the actual data now..

    const string tau_wall_str = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str       = zone_ptr->getName() + ":" + "y_plus";
    const string qw_str       = zone_ptr->getName() + ":" + "q_wall";

    // all of the following functions do not take any arguments -- check first

    if ( name == tau_wall_str) {

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {
        if ( args.size() == 0) {
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            const int icv           = zone_ptr->cvobf[ibf];
            const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
            double du_mag = 0.0;
            for (int i =0; i < 3; ++i)
              du_mag  += (solver->u[icv][i] - u_bc[i])*(solver->u[icv][i] - u_bc[i]);

            tmp[ibf]                = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
          }
        }
        else if (args.size() == 1) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          double dir[3] = {0.0,0.0,0.0};

          if (arg->getType() == D3_DATA) {
            FOR_I3 dir[i] = arg->d3(i);
            NORMALIZE(dir);
          }
          else if (arg->getType() == D_DATA) {
            const int index = int(arg->d());
            if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;
            dir[index] = 1.0;
          }
          else return CTI_DATA_NOT_VALID;

          // loop faces and compute
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            const int icv           = zone_ptr->cvobf[ibf];
            const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];

            double local_u[3] = DIFF(solver->u[icv],u_bc);
            tmp[ibf]          = visc_coeff*DOT_PRODUCT(dir,local_u)/zone_ptr->area_bf[ibf];
          }

        }
        else return CTI_DATA_ARG_COUNT;
      }

      return CTI_DATA_OK;

    }

    else if ( name == yp_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          const int icv           = zone_ptr->cvobf[ibf];
          const double y1         = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
          const double visc_coeff = solver->mu_lam[icv]/y1;

          double du_mag = 0.0;
          for (int i =0; i < 3; ++i)
            du_mag  += (solver->u[icv][i] - u_bc[i])*(solver->u[icv][i] - u_bc[i]);

          const double tau_wall   = visc_coeff*sqrt(du_mag);
          tmp[ibf]                = y1* sqrt(tau_wall*solver->rho[icv])/solver->mu_lam[icv];

        }

      }

      return CTI_DATA_OK;

    } else if ( name == qw_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          const int icv        = zone_ptr->cvobf[ibf];

          // use the Cp from the cell attached to the boundary...
          const double cp      = solver->R[icv]*solver->gamma[icv]/(solver->gamma[icv]-1.0);
          const double k_coeff = solver->loc_lam[icv]*zone_ptr->area_over_delta_bf[ibf]*cp;
          tmp[ibf]             = k_coeff*abs(solver->T[icv] - T_bc)/zone_ptr->area_bf[ibf];
        }

      }

      return CTI_DATA_OK;

    }
  }

  return PremixedBc::funcEvalCtiData(v,name,args,b_eval_func);
}

void WmAlgIsothermalP::initData() {

  if (mpi_rank == 0) cout << "WmAlgIsothermalP::initData() for zone " << zone_ptr->getName() << endl;

  assert(tau_wall == NULL); tau_wall = new double[zone_ptr->nbf];
  assert(u1       == NULL); u1       = new double[zone_ptr->nbf];
  assert(q_wall   == NULL); q_wall   = new double[zone_ptr->nbf];
  assert(T_wall   == NULL); T_wall   = new double[zone_ptr->nbf];

  Param * param = getParam(zone_ptr->getName());
  assert(param);

  if (param->size() < 2) {
    CERR("WmAlgIsothermalP::initData() for zone " << zone_ptr->getName() << ": format problems. Valid examples:\n" <<
	 zone_ptr->getName() << " WM_ALG_ISOTHERMAL <temperature>\n" <<
	 zone_ptr->getName() << " WM_ALG_ISOTHERMAL ASCII_FILE <filename> [TSHIFT <double>]");
  }
  else if (param->size() == 2) {
    const double T_bc = param->getDouble(1);
    if (mpi_rank == 0) cout << " > applying constant T " << T_bc << endl;
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
      T_wall[ibf] = T_bc;
  }
  else {
    // must be the form;
    // x0 WM_ALG_ISOTHERMAL ASCII_FILE <filename> [TSHIFT <double>]
    int iarg = 1;
    bool b_ascii_filename = false;
    string filename;
    double shift = 0.0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "ASCII_FILE") {
	filename = param->getString(iarg++);
	b_ascii_filename = true;
	if (mpi_rank == 0) cout << " > ASCII_FILE: " << filename << endl;
      }
      else if (token == "TSHIFT") {
	shift = param->getDouble(iarg++);
	if (mpi_rank == 0) cout << " > TSHIFT: " << shift << endl;
      }
      else {
	CERR("WmAlgIsothermalP::initData() for zone " << zone_ptr->getName() << ": format problems. Valid examples:\n" <<
	     zone_ptr->getName() << " WM_ALG_ISOTHERMAL <temperature>\n" <<
	     zone_ptr->getName() << " WM_ALG_ISOTHERMAL ASCII_FILE <filename> [TSHIFT <double>]");
      }
    }
    if (b_ascii_filename) {
      double * x = NULL;
      double * y = NULL;
      double * z = NULL;
      double * T = NULL;
      int n;
      // read on 1 rank and Bcast...
      if (mpi_rank == 0) {
	if (xcol4(x,y,z,T,n,filename.c_str(),1,2,3,4) != 0) {
	  assert(x == NULL);
	  assert(y == NULL);
	  assert(z == NULL);
	  assert(T == NULL);
	  n = -1;
	}
      }
      MPI_Bcast(&n,1,MPI_INT,0,mpi_comm);
      if (n == -1) {
	CERR("WmAlgIsothermalP::initData() for zone " << zone_ptr->getName() << ": problem with ASCII_FILE " << filename);
      }
      if (mpi_rank != 0) {
	x = new double[n];
	y = new double[n];
	z = new double[n];
	T = new double[n];
      }
      MPI_Bcast(x,n,MPI_DOUBLE,0,mpi_comm);
      MPI_Bcast(y,n,MPI_DOUBLE,0,mpi_comm);
      MPI_Bcast(z,n,MPI_DOUBLE,0,mpi_comm);
      MPI_Bcast(T,n,MPI_DOUBLE,0,mpi_comm);
      // loop through the local points and do an n^2 search to find the closest T...
      double my_minmax[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
	int i_min = -1;
	double d2_min;
	for (int i = 0; i < n; ++i) {
	  const double d2 =
	    pow(zone_ptr->x_bf[ibf][0]-x[i],2) +
	    pow(zone_ptr->x_bf[ibf][1]-y[i],2) +
	    pow(zone_ptr->x_bf[ibf][2]-z[i],2);
	  if ((i_min == -1)||(d2 < d2_min)) {
	    i_min = i;
	    d2_min = d2;
	  }
	}
	assert(i_min != -1);
	T_wall[ibf] = T[i_min] + shift; // add the shift
	my_minmax[0] = min(my_minmax[0],-d2_min);
	my_minmax[1] = min(my_minmax[1],T_wall[ibf]);
	my_minmax[2] = min(my_minmax[2],-T_wall[ibf]);
      }
      delete[] x;
      delete[] y;
      delete[] z;
      delete[] T;
      double minmax[3];
      MPI_Reduce(my_minmax,minmax,3,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
      if (mpi_rank == 0)
	cout << " > x_bf match: " << sqrt(-minmax[0]) << ", T range: " << minmax[1] << " " << -minmax[2] << endl;
    }
    else {
      CERR("WmAlgIsothermalP::initData() for zone " << zone_ptr->getName() << ": format problems. Valid examples:\n" <<
           zone_ptr->getName() << " WM_ALG_ISOTHERMAL <temperature>\n" <<
           zone_ptr->getName() << " WM_ALG_ISOTHERMAL ASCII_FILE <filename> [TSHIFT <double>]");
    }
  }

}

void WmAlgIsothermalP::initialHook() {

  if ( !zone_ptr->checkDataFlag("u1")) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      const int icv = zone_ptr->cvobf[ibf];
      u1[ibf] = MAG(solver->u[icv]);

      // set an initial value for tau_wall, this will constitute the
      // initial guess for the solution when it is finally constructed.

      if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {
        const double y1 = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        tau_wall[ibf]   = solver->mu_lam[icv] * MAG(solver->u[icv])/y1;
        q_wall[ibf]     = 0.0;
      }
      else {
        tau_wall[ibf]   = 0.0;
        q_wall[ibf]     = 0.0;
      }
    }
  }
}

int WmAlgIsothermalP::storeData(double * buf) const {

  // the wall model needs to store tau_wall because it is used as
  // an initial guess. No need to store u1. It is not used...

  if (buf) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      buf[ibf] = tau_wall[ibf];
    }
  }

  return zone_ptr->nbf;

}

int WmAlgIsothermalP::restoreData(double * buf) {

  assert(buf);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    tau_wall[ibf] = buf[ibf];
  }

  return zone_ptr->nbf;

}

WmAlgIsothermalP::~WmAlgIsothermalP() {

  DELETE(u1);
  DELETE(tau_wall);
  DELETE(q_wall);
  DELETE(T_wall);

}

void WmAlgIsothermalP::calcRhs(const double time, const int rk_stage) {

  // see notes left in IdealGasSolverBcs.cpp for the time
  // filtering approach

  if ( rk_stage == 1 ) {

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      const int icv = zone_ptr->cvobf[ibf];

      if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

        const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        const double mag_u = MAG(solver->u[icv]);
        //const double eps   = mag_u*solver->dt / (y1 + solver->dt*mag_u);

        // \phi_new = eps*u + (1.0 -eps)*\phi_old
        // \phi_new - \phi_old = eps*(u - \phi_old)

        //u1[ibf]   += eps* ( mag_u - u1[ibf]);
        u1[ibf]   = mag_u;

        // use the existing tau_wall to set a guess for the u_tau

        const double u_tau = sqrt(tau_wall[ibf] / solver->rho[icv]);
        tau_wall[ibf] = AlgebraicWM::solve_tau(u1[ibf], y1, solver->rho[icv],
                                               solver->mu_lam[icv], u_tau);

        // approximate material properties from the interior -- these are low Ma approximations
        // that ignore property changes in the near wall cell due to compressibility

        const double cp = solver->R[icv]*solver->gamma[icv]/(solver->gamma[icv]-1.0);

        q_wall[ibf] = AlgebraicWM::compute_q_wall_approx(solver->T[icv], T_wall[ibf], tau_wall[ibf],
                                                         solver->rho[icv], solver->mu_lam[icv], solver->loc_lam[icv], cp, y1);

      }
      else {

        //assert ( zone_ptr->area_bf[ibf] < 1.0e-16); // this should be a collapsed face.
        tau_wall[ibf]  = 0.0;
        q_wall[ibf] = 0.0;
        u1[ibf]     = MAG(solver->u[icv]);

      }
    }
  }
}

void WmAlgIsothermalP::addBoundaryFlux(PremixedRhs* rhs) const {

  assert( solver != NULL);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    // tau_wall has already been computed at this point, but we'll assume
    // that the wall stress is instantaneously aligned with the velocity vector

    double u_mag = 0.0;
    for (int i =0; i < 3; ++i)
      u_mag += solver->u[icv][i] * solver->u[icv][i];
    u_mag = sqrt( max(0.0, u_mag));

    if ( u_mag > 0.0) {

      for (int i =0; i < 3;++i) {
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];
        rhs[icv].rhou[i] -= tau_wall[ibf]*solver->u[icv][i]/u_mag*zone_ptr->area_bf[ibf];
      }
      // XXX if u_bc \neq 0, there should be a u \dot \tau term added here..
      rhs[icv].rhoE -= q_wall[ibf]*zone_ptr->area_bf[ibf];

    } else {

      for (int i =0; i < 3; ++i)
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];

      // XXX if u_bc \neq 0, there should be a u \dot \tau term added here..
      rhs[icv].rhoE -= q_wall[ibf]*zone_ptr->area_bf[ibf];

    }
  }
}

CtiRegister::CtiDataError WmAlgIsothermalP::funcEvalCtiData(CtiRegister::CtiData& v,
                                                            const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if ( zone_ptr->isLeadingStrMatched(name) ) {

    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix.  we have to check for the actual data now..

    // for consistency, we are supporting a tau_wall(), y_plus()
    // function to have the same syntax as all of the other wall zones..

    const string tau_wall_str = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str  = zone_ptr->getName() + ":" + "y_plus";
    const string qw_str  = zone_ptr->getName() + ":" + "q_wall";


    if ( name == tau_wall_str) {

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {
        if (args.size() == 0) {
          // use what was computed by wall model
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
            tmp[ibf] = tau_wall[ibf];
          }
        }
        else if (args.size() == 1) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          double dir[3] = {0.0,0.0,0.0};

          if (arg->getType() == D3_DATA) {
            FOR_I3 dir[i] = arg->d3(i);
            NORMALIZE(dir);
          }
          else if (arg->getType() == D_DATA) {
            const int index = int(arg->d());
            if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;
            dir[index] = 1.0;
          }
          else return CTI_DATA_NOT_VALID;

          // loop faces and compute
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
            const int icv      = zone_ptr->cvobf[ibf];
            const double mag_u = MAG(solver->u[icv]);

            double unit_u[3];
            if ( mag_u > 0.0) {
              for (int i = 0; i < 3; ++i) unit_u[i] = solver->u[icv][i]/mag_u;

              tmp[ibf] = tau_wall[ibf]*DOT_PRODUCT(dir,unit_u);
            }
            else {
              tmp[ibf] = 0.0;
            }
          }
        }
        else return CTI_DATA_ARG_COUNT;
      }

      return CTI_DATA_OK;

    } else if ( name == yp_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          const int icv      = zone_ptr->cvobf[ibf];
          const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
          tmp[ibf]           = y1*sqrt(solver->rho[icv]*tau_wall[ibf])/solver->mu_lam[icv];
        }

      }

      return CTI_DATA_OK;

    } else if ( name == qw_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
          tmp[ibf] = q_wall[ibf];

      }

      return CTI_DATA_OK;

    }
  }

  return PremixedBc::funcEvalCtiData(v,name,args,b_eval_func);

}

void WmAlgIsothermalP::query(const string& param_str) const {

  // report an estimate of the local y^+ on the boundary and the integrated tau_wall

  double my_buf[4] = {0.0,0.0,0.0,0.0};
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

      const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
      const double u_tau = sqrt( tau_wall[ibf]/solver->rho[icv]);
      const double nu    = solver->mu_lam[icv] / solver->rho[icv];
      const double y_plus = y1*u_tau/nu;

      my_buf[0] += zone_ptr->area_bf[ibf];
      my_buf[1] += zone_ptr->area_bf[ibf]*tau_wall[ibf];
      my_buf[2] += zone_ptr->area_bf[ibf]*y_plus;
      my_buf[3] += zone_ptr->area_bf[ibf]*q_wall[ibf];

    }
  }

  double buf[4];
  MPI_Reduce(my_buf,buf,4,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if ( mpi_rank == 0 ) {

    cout << " > " << zone_ptr->getName() << ": time, int(tau_wall)dA, int(q_w)dA, avg(y_plus): "
         << solver->time << "    "
         << buf[1] << "    "
         << buf[3] << "    "
         << buf[2]/buf[0] << endl;

  }
}

void WallChtP::initData() {
  if (solver->cht == NULL) CERR("WALL_CHT bcs require CHT");

  int * ssp_flag = new int[solver->subSurface->nsp];
  for (int issp = 0; issp < solver->subSurface->nsp; ++issp)
    ssp_flag[issp] = -1;

  for (int ibf = zone_ptr->ibf_f; ibf <= zone_ptr->ibf_l; ++ibf) {
    for (int sob = solver->sspobf_i[ibf]; sob != solver->sspobf_i[ibf+1]; ++sob) {
      const int issp = solver->sspobf_v[sob];
      ssp_flag[issp] = 0;
    }
  }

  np = 0;
  for (int issp = 0; issp < solver->subSurface->nsp; ++issp)
    if (ssp_flag[issp] == 0)
      ssp_flag[issp] = np++;

  assert(Tp == NULL);
  Tp = new double[np];
  for (int ip = 0; ip < np; ++ip) Tp[ip] = HUGE_VAL; // put garbage here to ensure it set before being used

  assert(qp == NULL);
  qp = new double[np];
  for (int ip = 0; ip < np; ++ip) qp[ip] = HUGE_VAL; // put garbage here to ensure it set before being used

  double (*xp)[3] = new double[np][3];
  int ip = 0;
  for (int issp = 0; issp < solver->subSurface->nsp; ++issp) {
    if (ssp_flag[issp] >= 0) {
      FOR_I3 xp[ip][i] = solver->subSurface->xp[issp][i];
      const int bits = int(solver->subSurface->isp_global_and_bits[issp]>>52);
      if (bits)
        PeriodicData::periodicTranslate(xp[ip],1,BitUtils::flipPeriodicBits(bits));
      ++ip;
    }
  }
  assert(ip == np);

  //solver->cht->registerBoundaryData(xp, Tp, qp, np, zone_ptr->index);
  solver->cht->registerBoundaryData(xp,Tp,qp,ssp_flag,np,solver->subSurface->nsp,zone_ptr->index);
  delete[] xp; // does not need to be persistent

  // build the ipobf_i/v/wgt needed to move information to/from the Tp/qp arrays...
  assert(ipobf_i == NULL);
  assert(zone_ptr->ibf_l-zone_ptr->ibf_f+1 == zone_ptr->nbf);
  ipobf_i = new int[zone_ptr->nbf+1];
  ipobf_i[0] = 0;
  const int ipobf_s = solver->sspobf_i[zone_ptr->ibf_l+1] - solver->sspobf_i[zone_ptr->ibf_f];
  assert(ipobf_v == NULL);
  ipobf_v = new int[ipobf_s];
  assert(ipobf_wgt == NULL);
  ipobf_wgt = new double[ipobf_s];
  int iob = 0;
  for (int ibf = zone_ptr->ibf_f; ibf <= zone_ptr->ibf_l; ++ibf) {
    double wgt_sum = 0.0;
    const int iob0 = iob;
    bool b_unmatched = false;
    for (int sob = solver->sspobf_i[ibf]; sob != solver->sspobf_i[ibf+1]; ++sob) {
      const int issp = solver->sspobf_v[sob];
      //assert((ssp_flag[issp] >= 0)&&(ssp_flag[issp] < np));
      if (ssp_flag[issp] >= 0) {
        assert(ssp_flag[issp] < np);
        ipobf_v[iob] = ssp_flag[issp];
        ipobf_wgt[iob] = solver->sspobf_wgt[sob];
        wgt_sum += ipobf_wgt[iob];
        ++iob;
      }
      else {
        b_unmatched = true;
      }
    }
    //assert(fabs(1.0-wgt_sum) < 1.0E-12);
    ipobf_i[ibf-zone_ptr->ibf_f+1] = iob;
    if (b_unmatched) {
      // redistribute unmatched ip weights...
      assert(wgt_sum > 1.0E-12);
      const double inv_wgt_sum = 1.0/wgt_sum;
      for (int iob_ = iob0; iob_ < iob; ++iob_)
        ipobf_wgt[iob_] *= inv_wgt_sum;
    }
  }
  //assert(iob == ipobf_s);
  assert(iob <= ipobf_s);
  delete[] ssp_flag;
}

void WallChtP::addBoundaryFlux(PremixedRhs* rhs) const {
  assert(solver != NULL);

  //zero the qp so it can be summed below...
  for (int ip = 0; ip < np; ++ip) qp[ip] = 0.0;

  // loop on bfs...
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    // reconstruct the boundary temperature from the Tp data that should
    // be current in the Tp[ip] array. (updated with an earlier call to
    // cht->updateChtTemperatureToFlowBcs())...
    double T_bc = 0.0;
    for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
      const int ip = ipobf_v[iob];
      const double wgt = ipobf_wgt[iob];
      T_bc += wgt*Tp[ip];
    }

    const int icv           = zone_ptr->cvobf[ibf];
    const double cp         = solver->R[icv]*solver->gamma[icv]/(solver->gamma[icv] - 1.0);
    const double visc_coeff = solver->mu_lam[icv] *zone_ptr->area_over_delta_bf[ibf];
    const double k_coeff    = solver->loc_lam[icv]*zone_ptr->area_over_delta_bf[ibf]*cp;

    for (int i = 0; i < 3; ++i)
      rhs[icv].rhou[i]     -= solver->p[icv] *zone_ptr->n_bf[ibf][i] + visc_coeff*(solver->u[icv][i] - 0.0);
    rhs[icv].rhoE          -= k_coeff*(solver->T[icv] - T_bc);

    // scatter/sum q_bf to the qp[ip] data points
    for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
      const int ip = ipobf_v[iob];
      const double wgt = ipobf_wgt[iob];
      qp[ip] += wgt*k_coeff*(solver->T[icv] - Tp[ip]); // should match -q_bf when summed over ip
    }
  }
}

void WallChtP::query(const string& param_str) const {

  double my_buf[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    // compute this...
    double T_bc = 0.0;
    for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
      const int ip = ipobf_v[iob];
      const double wgt = ipobf_wgt[iob];
      T_bc += wgt*Tp[ip];
    }

    if (zone_ptr->area_over_delta_bf[ibf] > 0.0) {

      const double cp = solver->R[icv]*solver->gamma[icv]/(solver->gamma[icv] - 1.0);
      const double y1 = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
      const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
      double du_mag = MAG(solver->u[icv]);

      const double tau_wall = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
      const double y_plus   = y1*sqrt(tau_wall*solver->rho[icv])/solver->mu_lam[icv];

      const double k_coeff  = solver->loc_lam[icv] * zone_ptr->area_bf[ibf]*cp;
      const double q_wall   = k_coeff*abs(solver->T[icv] - T_bc)/zone_ptr->area_bf[ibf];

      my_buf[0]            += zone_ptr->area_bf[ibf];
      my_buf[1]            += zone_ptr->area_bf[ibf]*tau_wall;
      my_buf[2]            += zone_ptr->area_bf[ibf]*y_plus;
      my_buf[3]            += zone_ptr->area_bf[ibf]*q_wall;
      my_buf[4]            += zone_ptr->area_bf[ibf]*T_bc;
    }
  }

  double buf[5];
  MPI_Reduce(my_buf, buf, 5, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);

  if (mpi_rank == 0)
    cout  << zone_ptr->getName() << " time, int(tau_wall)dA, int(q_wall)dA, avg y_plus, avg T_bc = "
          << solver->time   << "\t" << buf[1]         << "\t"
          << buf[3]         << "\t" << buf[2]/buf[0]  << "\t"
          << buf[4]/buf[0] <<  endl;

}

CtiRegister::CtiDataError WallChtP::funcEvalCtiData(CtiRegister::CtiData& v,
                                                    const string& name, list <CtiRegister::CtiData> & args, const bool b_eval_func) {

  using namespace CtiRegister;

  if (zone_ptr->isLeadingStrMatched(name)) {

    // we checked if the requested function is a potential match for the data
    // basedon the zone name prefix. We have to check for the actual data now...

    const string tau_wall_str = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str       = zone_ptr->getName() + ":" + "y_plus";
    const string qw_str       = zone_ptr->getName() + ":" + "q_wall";
    const string t_wall_str = zone_ptr->getName() + ":" + "T_wall";

    if (name == tau_wall_str) {

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {
        if (args.size() == 0) {
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
            const int icv = zone_ptr->cvobf[ibf];
            const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
            double du_mag = MAG(solver->u[icv]);
            tmp[ibf] = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
          }
        }
        else if (args.size() == 1) {
          list <CtiRegister::CtiData>::iterator arg = args.begin();
          double dir[3] = {0.0,0.0,0.0};

          if (arg->getType() == D3_DATA) {
            FOR_I3 dir[i] = arg->d3(i);
            NORMALIZE(dir);
          }
          else if (arg->getType() == D_DATA) {
            const int index = int(arg->d());
            if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;
            dir[index] = 1.0;
          }
          else return CTI_DATA_NOT_VALID;

          // loop faces and compute
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ ibf) {

            const int icv             = zone_ptr->cvobf[ibf];
            const double visc_coeff   = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
            tmp[ibf]                  = visc_coeff*DOT_PRODUCT(dir, solver->u[icv])/zone_ptr->area_bf[ibf];
          }
        }
        else return CTI_DATA_ARG_COUNT;
      }
      return CTI_DATA_OK;
    }
    else if (name == yp_str) {
      if (args.size() != 0) return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          const int icv = zone_ptr->cvobf[ibf];
          const double y1 = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
          const double visc_coeff = solver->mu_lam[icv]/y1;

          double du_mag = MAG(solver->u[icv]);
          const double tau_wall = visc_coeff*sqrt(du_mag);
          tmp[ibf] = y1*sqrt(tau_wall*solver->rho[icv])/solver->mu_lam[icv];
        }
      }

      return CTI_DATA_OK;
    }
    else if (name == qw_str) {
      if (args.size() != 0) return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          double T_bc = 0.0;
          for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
            const int ip = ipobf_v[iob];
            const double wgt = ipobf_wgt[iob];
            T_bc += wgt*Tp[ip];
          }

          const int icv = zone_ptr->cvobf[ibf];
          const double cp = solver->R[icv]*solver->gamma[icv]/(solver->gamma[icv] - 1.0);
          const double k_coeff = solver->loc_lam[icv]*zone_ptr->area_over_delta_bf[ibf]*cp;
          tmp[ibf]            = k_coeff*abs(solver->T[icv] - T_bc)/zone_ptr->area_bf[ibf];
        }
      }
      return CTI_DATA_OK;
    }
    else if (name == t_wall_str) {

      if (args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          // reconstruct wall temperature from Tp...
          double T_wall = 0.0;
          for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
            const int ip = ipobf_v[iob];
            const double wgt = ipobf_wgt[iob];
            T_wall += wgt*Tp[ip];
          }
          tmp[ibf] = T_wall;
        }
      }

      return CTI_DATA_OK;

    }
  }

  return PremixedBc::funcEvalCtiData(v,name,args,b_eval_func);
}

WmAlgChtP::WmAlgChtP(BfZone* p, PremixedSolver* s) : PremixedBc(p,s) {

  if (mpi_rank == 0) cout << "WmAlgChtP() for zone: " << p->getName() << endl;

  u1       = NULL; zone_ptr->registerBfData(u1,"u1",READWRITE_DATA); // saved for filtering (if used)
  tau_wall = NULL; zone_ptr->registerBfData(tau_wall,"tau_wall",READWRITE_DATA); // saved because it is an initial guess
  q_wall   = NULL; zone_ptr->registerBfData(q_wall,"q_wall",READWRITE_DATA); // saved for no reason
  k_eff    = NULL;

  // allow up to 3 layers.
  nlayer = 0;
  FOR_I3 T_layer[i] = NULL; // double * for the layer fluid-side temperatures: T_layer1, T_layer2, T_layer3
  FOR_I3 k_layer[i] = NULL; // SimpleFunc pointers for the k's...

  np = 0;
  Tp = NULL;
  qp = NULL;
  ipobf_i = NULL;
  ipobf_v = NULL;
  ipobf_wgt = NULL;

  // add to funcEvalList
  funcEvalList.push_back(p->getName()+":tau_wall()");
  funcEvalList.push_back(p->getName()+":y_plus()");
  funcEvalList.push_back(p->getName()+":q_wall()");
  funcEvalList.push_back(p->getName()+":T_wall()");

  // set the load balance cost. Default for bc is 1
  // see note above
  /*
    if (mpi_rank == 0) cout << "XXXXX setting lb_cost to 5" << endl;
    assert(zone_ptr->lb_cost == 1);
    zone_ptr->lb_cost = 5;
    assert(p->lb_cost == 5);
  */

  // look for optional thin layer model associated with this bc...

  Param * param = getParam(getName());
  assert(param);

  vector<pair<string,double> > layerVec;
  int iarg = 1;
  while (iarg < param->size()) {
    string token = param->getString(iarg++);
    if (token == "LAYER") {
      token = param->getString(iarg++);
      const double delta = param->getDouble(iarg++);
      if (mpi_rank == 0) cout << " > adding LAYER " << token << " with thickness: " << delta << endl;
      layerVec.push_back(pair<string,double>(token,delta));
    }
    else {
      CERR("unrecognized WM_ALG_CHT param: " << token);
    }
  }

  if (!layerVec.empty()) {
    nlayer = layerVec.size();
    if (nlayer > 3) {
      CERR("number of LAYERs in WM_ALG_CHT limited to 3 for now");
    }
    // register a temperature associated with the fluid side of each layer...
    char vname[32];
    for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
      // register the data with read/write...
      // T_layer1,T_layer2...
      sprintf(vname,"T_layer%d",ilayer+1);
      zone_ptr->registerBfData(T_layer[ilayer],vname,READWRITE_DATA);
      // set the layer delta...
      delta_layer[ilayer] = layerVec[ilayer].second;
      // and k+layer is a SimpleFunc...
      assert(k_layer[ilayer] == NULL);
      Param * k_param = getParam(layerVec[ilayer].first+".K");
      if (k_param == NULL) {
        CERR("LAYER " << layerVec[ilayer].first << " requires param for thermal conductivity " << layerVec[ilayer].first << ".K");
      }
      k_layer[ilayer] = processSimpleFunc(k_param);
      assert(k_layer[ilayer] != NULL);
    }
  }

}

void WmAlgChtP::initData() {

  if (mpi_rank == 0) cout << "WmAlgChtP::initData() for zone " << zone_ptr->getName() << endl;

  // we should have CHT by now...
  if (solver->cht == NULL) {
    CERR("WM_ALG_CHT bcs require CHT");
  }

  assert(u1 == NULL);        u1 = new double[zone_ptr->nbf];
  assert(tau_wall == NULL);  tau_wall = new double[zone_ptr->nbf];
  assert(q_wall == NULL);    q_wall = new double[zone_ptr->nbf];
  assert(k_eff == NULL);     k_eff = new double[zone_ptr->nbf]; // not registered

  // set everything to HUGE_VAL...
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    u1[ibf] = HUGE_VAL;
    tau_wall[ibf] = HUGE_VAL;
    q_wall[ibf] = HUGE_VAL;
    k_eff[ibf] = HUGE_VAL;
  }

  // and allocate and set any layer temperatures to HUGE_VAL...
  assert(nlayer <= 3);
  for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
    assert(T_layer[ilayer] == NULL);
    T_layer[ilayer] = new double[zone_ptr->nbf];
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
      T_layer[ilayer][ibf] = HUGE_VAL;
  }

  int * ssp_flag = new int[solver->subSurface->nsp];
  for (int issp = 0; issp < solver->subSurface->nsp; ++issp)
    ssp_flag[issp] = -1;

  for (int ibf = zone_ptr->ibf_f; ibf <= zone_ptr->ibf_l; ++ibf) {
    for (int sob = solver->sspobf_i[ibf]; sob != solver->sspobf_i[ibf+1]; ++sob) {
      const int issp = solver->sspobf_v[sob];
      ssp_flag[issp] = 0;
    }
  }

  np = 0;
  for (int issp = 0; issp < solver->subSurface->nsp; ++issp)
    if (ssp_flag[issp] == 0)
      ssp_flag[issp] = np++;

  assert(Tp == NULL);
  Tp = new double[np];
  for (int ip = 0; ip < np; ++ip) Tp[ip] = HUGE_VAL; // put garbage here to ensure it is set before being used

  assert(qp == NULL);
  qp = new double[np];
  for (int ip = 0; ip < np; ++ip) qp[ip] = HUGE_VAL; // put garbage here to ensure it is set before being used

  double (*xp)[3] = new double[np][3];
  int ip = 0;
  for (int issp = 0; issp < solver->subSurface->nsp; ++issp) {
    if (ssp_flag[issp] >= 0) {
      FOR_I3 xp[ip][i] = solver->subSurface->xp[issp][i];
      const int bits = int(solver->subSurface->isp_global_and_bits[issp]>>52);
      if (bits)
        PeriodicData::periodicTranslate(xp[ip],1,BitUtils::flipPeriodicBits(bits));
      ++ip;
    }
  }
  assert(ip == np);

  // register this bc with the CHT registration...
  //solver->cht->registerBoundaryData(xp,Tp,qp,np,zone_ptr->index);
  solver->cht->registerBoundaryData(xp,Tp,qp,ssp_flag,np,solver->subSurface->nsp,zone_ptr->index);
  delete[] xp;

  // build the ipobf_i/v/wgt needed to move information to/from the Tp/qp arrays...
  assert(ipobf_i == NULL);
  assert(zone_ptr->ibf_l-zone_ptr->ibf_f+1 == zone_ptr->nbf);
  ipobf_i = new int[zone_ptr->nbf+1];
  ipobf_i[0] = 0;
  const int ipobf_s = solver->sspobf_i[zone_ptr->ibf_l+1] - solver->sspobf_i[zone_ptr->ibf_f];
  assert(ipobf_v == NULL);
  ipobf_v = new int[ipobf_s];
  assert(ipobf_wgt == NULL);
  ipobf_wgt = new double[ipobf_s];
  int iob = 0;
  for (int ibf = zone_ptr->ibf_f; ibf <= zone_ptr->ibf_l; ++ibf) {
    double wgt_sum = 0.0;
    const int iob0 = iob;
    bool b_unmatched = false;
    for (int sob = solver->sspobf_i[ibf]; sob != solver->sspobf_i[ibf+1]; ++sob) {
      const int issp = solver->sspobf_v[sob];
      //assert((ssp_flag[issp] >= 0)&&(ssp_flag[issp] < np));
      if (ssp_flag[issp] >= 0) {
        assert(ssp_flag[issp] < np);
        ipobf_v[iob] = ssp_flag[issp];
        ipobf_wgt[iob] = solver->sspobf_wgt[sob];
        wgt_sum += ipobf_wgt[iob];
        ++iob;
      }
      else {
        b_unmatched = true;
      }
    }
    //assert(fabs(1.0-wgt_sum) < 1.0E-12);
    ipobf_i[ibf-zone_ptr->ibf_f+1] = iob;
    if (b_unmatched) {
      // redistribute unmatched ip weights...
      assert(wgt_sum > 1.0E-12);
      const double inv_wgt_sum = 1.0/wgt_sum;
      for (int iob_ = iob0; iob_ < iob; ++iob_)
        ipobf_wgt[iob_] *= inv_wgt_sum;
    }
  }
  //assert(iob == ipobf_s);
  assert(iob <= ipobf_s);
  delete[] ssp_flag;

}

void WmAlgChtP::initialHook() {

  // This is a better way to set boundary data. If we find any
  // single HUGE_VAL in any of the boundary data that should have been
  // read from the restart, we reset everything...

  bool b_init = false;
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    if (u1[ibf] == HUGE_VAL) {
      b_init = true;
      break;
    }
    if (tau_wall[ibf] == HUGE_VAL) {
      b_init = true;
      break;
    }
    if (q_wall[ibf] == HUGE_VAL) {
      b_init = true;
      break;
    }
    for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
      if (T_layer[ilayer][ibf] == HUGE_VAL) {
        b_init = true;
        break;
      }
    }
    if (b_init)
      break;
  }

  // --------------------------------------------------------------------
  // if even one stored value was HUGE_VAL, initialize the whole bc...
  // --------------------------------------------------------------------

  if (b_init) {

    if (mpi_rank == 0) cout << "WM_ALG_CHT " << zone_ptr->getName() << ": initializing all data" << endl;

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv = zone_ptr->cvobf[ibf];
      u1[ibf] = MAG(solver->u[icv]);
      // set an initial value for tau_wall, this will constitute the
      // initial guess for the solution when it is finally constructed.
      if (zone_ptr->area_over_delta_bf[ibf] > 0.0) {
        const double y1 = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        tau_wall[ibf] = solver->mu_lam[icv]*MAG(solver->u[icv])/y1;
      }
      else {
        tau_wall[ibf] = 0.0;
      }
      q_wall[ibf] = 0.0;
      k_eff[ibf] = 0.0;
      for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
        T_layer[ilayer][ibf] = solver->T[icv]; // initialize any layer temps to the fluid T
      }
    }

  }

}

int WmAlgChtP::storeData(double * buf) const {

  // the wall model needs to store tau_wall because it is used as
  // an initial guess. No need to store u1. It is not used...
  // q_wall is computed. If there are any layer temperatures, these
  // need to be stored too...

  // TODO: if you reintroduce filtering, we need u1 here...

  if (buf) {
    int count = 0;
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      buf[count++] = tau_wall[ibf];
    }
    for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
        buf[count++] = T_layer[ilayer][ibf];
      }
    }
    assert(count == zone_ptr->nbf*(1+nlayer));
  }

  return zone_ptr->nbf*(1+nlayer);

}

int WmAlgChtP::restoreData(double * buf) {

  assert(buf);
  int count = 0;
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    tau_wall[ibf] = buf[count++];
  }
  for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      T_layer[ilayer][ibf] = buf[count++];
    }
  }
  assert(count == zone_ptr->nbf*(1+nlayer));

  return zone_ptr->nbf*(1+nlayer);

}

WmAlgChtP::~WmAlgChtP() {

  DELETE(u1);
  DELETE(tau_wall);
  DELETE(q_wall);
  DELETE(k_eff);

  FOR_I3 DELETE(T_layer[i]);
  FOR_I3 if (k_layer[i]) delete k_layer[i];

  DELETE(Tp);
  DELETE(qp);
  DELETE(ipobf_i);
  DELETE(ipobf_v);
  DELETE(ipobf_wgt);

}

void WmAlgChtP::calcRhs(const double time, const int rk_stage) {

  // see notes left in IdealGasSolverBcs.cpp for the time
  // filtering approach

  assert((rk_stage >= 1)&&(rk_stage <= 3));
  if (rk_stage == 1) {

    // on the first stage, do the wall model calcs...

    double thermal_resistance[3];
    const int nip = 7; // number of interpolation points for T-dependent k in each layer
    double Tip[nip];
    double kip[nip];

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      const int icv = zone_ptr->cvobf[ibf];

      if (zone_ptr->area_over_delta_bf[ibf] > 0.0) {

        const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        const double mag_u = MAG(solver->u[icv]);
        //const double eps   = mag_u*solver->dt / (y1 + solver->dt*mag_u);

        // \phi_new = eps*u + (1.0 -eps)*\phi_old
        // \phi_new - \phi_old = eps*(u - \phi_old)

        //u1[ibf]   += eps* ( mag_u - u1[ibf]);
        u1[ibf]   = mag_u;

        // use the existing tau_wall to set a guess for the u_tau

        const double u_tau = sqrt(tau_wall[ibf] / solver->rho[icv]);
        tau_wall[ibf] = AlgebraicWM::solve_tau(u1[ibf], y1, solver->rho[icv],
					       solver->mu_lam[icv], u_tau);

        // approximate material properties from the interior -- these are low Ma approximations
        // that ignore property changes in the near wall cell due to compressibility

        const double cp = solver->R[icv]*solver->gamma[icv]/(solver->gamma[icv]-1.0);

        // reconstruct the wall temperature from the Tp data that should
        // be current in the Tp[ip] array. (updated with an earlier call to
        // cht->updateChtTemperatureToFlowBcs())...
        double T_wall = 0.0;
        for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
          const int ip = ipobf_v[iob];
          const double wgt = ipobf_wgt[iob];
          T_wall += wgt*Tp[ip];
        }

        // recall that flux and thermal resistance are related by:
        //
        // q = deltaT / sum(R)
        //
        // where sum(R) = dx1/k1 + dx2/k2...

        // determine the thermal resistance of each layer and sum it...
        double thermal_resistance_sum = 0.0;
        double T_layer_prev = T_wall; // start at the cht wall
        for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
          // interpolate the temperatures to the nip points inside the layer...
          for (int ip = 0; ip < nip; ++ip) {
            const double w1 = (double(ip)+0.5)/double(nip);
            Tip[ip] = (1.0-w1)*T_layer_prev + w1*T_layer[ilayer][ibf];
          }
          k_layer[ilayer]->eval(kip,Tip,nip);
          thermal_resistance[ilayer] = 0.0;
          for (int ip = 0; ip < nip; ++ip) {
            thermal_resistance[ilayer] += delta_layer[ilayer]/(kip[ip]*double(nip));
          }
          thermal_resistance_sum += thermal_resistance[ilayer];
          T_layer_prev = T_layer[ilayer][ibf];
        }

        //static bool first = (mpi_rank == 132); // debugging

        // now call the wall model routine that includes a thermal resistance...
        q_wall[ibf] = AlgebraicWM::compute_q_wall_approx_r(solver->T[icv], T_wall, thermal_resistance_sum, tau_wall[ibf],
                                                           solver->rho[icv], solver->mu_lam[icv], solver->loc_lam[icv], cp, y1, false);

        // reset the layer temperatures based on thie new flux and the current T_wall...
        T_layer_prev = T_wall; // start at the cht wall
        for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
          T_layer[ilayer][ibf] = T_layer_prev + q_wall[ibf]*thermal_resistance[ilayer];
          T_layer_prev = T_layer[ilayer][ibf];
        }

        /*
        if (first && nlayer) {
          cout << "just got q_wall[ibf]: " << q_wall[ibf] << " T[icv]: " << solver->T[icv] << " T_wall: " << T_wall << " thermal_resistance_sum: " << thermal_resistance_sum << endl;
          for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
            cout << "layer " << ilayer << " temperature: " << T_layer[ilayer][ibf] << " thermal_resistance[ilayer]: " << thermal_resistance[ilayer] << endl;
          }
          first = false;
        }
        */

        // and compute the effective resistance of the wall model + thermal resistance of the layers...
        // and save as k_eff (1/R) including the bf area...
        const double dT = solver->T[icv] - T_wall;
	if (dT != 0.0) {
          k_eff[ibf] = zone_ptr->area_bf[ibf]*q_wall[ibf]/dT;
        }
	else {
	  // just use the laminar value...
	  k_eff[ibf] = zone_ptr->area_bf[ibf]/(y1/(solver->loc_lam[icv]*cp) + thermal_resistance_sum);
	}

      }
      else {

        // this should be a collapsed face.
        u1[ibf] = MAG(solver->u[icv]);
        tau_wall[ibf] = 0.0;
        q_wall[ibf] = 0.0;
	k_eff[ibf] = 0.0;

      }
    }

  }

}

void WmAlgChtP::addBoundaryFlux(PremixedRhs* rhs) const {

  assert( solver != NULL);

  // zero the qp so it can be summed below...
  for (int ip = 0; ip < np; ++ip) qp[ip] = 0.0;

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    // tau_wall has already been computed at this point, but we'll assume
    // that the wall stress is instantaneously aligned with the velocity vector

    double u_mag = 0.0;
    for (int i =0; i < 3; ++i)
      u_mag += solver->u[icv][i] * solver->u[icv][i];
    u_mag = sqrt(max(0.0,u_mag));

    if ( u_mag > 0.0) {
      for (int i =0; i < 3;++i) {
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];
        rhs[icv].rhou[i] -= tau_wall[ibf]*solver->u[icv][i]/u_mag*zone_ptr->area_bf[ibf];
      }
    }
    else {
      for (int i =0; i < 3; ++i)
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];
    }

    double T_wall = 0.0;
    for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
      const int ip = ipobf_v[iob];
      const double wgt = ipobf_wgt[iob];
      T_wall += wgt*Tp[ip];
    }

    // XXX if u_bc \neq 0, there should be a u \dot \tau term added here..
    //rhs[icv].rhoE -= q_wall[ibf]*zone_ptr->area_bf[ibf];
    rhs[icv].rhoE -= k_eff[ibf]*(solver->T[icv] - T_wall);

    // scatter/sum q_bf to the qp[ip] data points for exchange to the FEM solve this sub-step...
    for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
      const int ip = ipobf_v[iob];
      const double wgt = ipobf_wgt[iob];
      qp[ip] += wgt*k_eff[ibf]*(solver->T[icv] - Tp[ip]); // should match -q_bf when summed over ip
    }

  }

}

CtiRegister::CtiDataError WmAlgChtP::funcEvalCtiData(CtiRegister::CtiData& v,
                                                     const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  // all of the following functions do not take any arguments -- check first

  if (name == zone_ptr->getName()+":tau_wall") {

    double * tmp = zone_ptr->createBfD1Data(v);

    if (b_eval_func) {
      if (args.size() == 0) {
        // use what was computed by wall model
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          tmp[ibf] = tau_wall[ibf];
        }
      }
      else if (args.size() == 1) {
        list<CtiRegister::CtiData>::iterator arg = args.begin();

        double dir[3] = {0.0,0.0,0.0};

        if (arg->getType() == D3_DATA) {
          FOR_I3 dir[i] = arg->d3(i);
          NORMALIZE(dir);
        }
        else if (arg->getType() == D_DATA) {
          const int index = int(arg->d());
          if ((index < 0)||(index > 2)) return CTI_DATA_NOT_VALID;
          dir[index] = 1.0;
        }
        else {
          return CTI_DATA_NOT_VALID;
        }

        // loop faces and compute
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          const int icv      = zone_ptr->cvobf[ibf];
          const double mag_u = MAG(solver->u[icv]);
          double unit_u[3];
          if ( mag_u > 0.0) {
            for (int i = 0; i < 3; ++i) unit_u[i] = solver->u[icv][i]/mag_u;
            tmp[ibf] = tau_wall[ibf]*DOT_PRODUCT(dir,unit_u);
          }
          else {
            tmp[ibf] = 0.0;
          }
        }
      }
      else {
        return CTI_DATA_ARG_COUNT;
      }
    }

    return CTI_DATA_OK;

  }
  else if (name == zone_ptr->getName()+":y_plus") {

    if ( args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double * tmp = zone_ptr->createBfD1Data(v);

    if ( b_eval_func) {
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
        const int icv      = zone_ptr->cvobf[ibf];
        const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        tmp[ibf]           = y1*sqrt(solver->rho[icv]*tau_wall[ibf])/solver->mu_lam[icv];
      }
    }

    return CTI_DATA_OK;

  }
  else if (name == zone_ptr->getName()+":q_wall") {

    if ( args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double * tmp = zone_ptr->createBfD1Data(v);

    if ( b_eval_func) {
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
        if (zone_ptr->area_bf[ibf] > 0.0) {
          double T_wall = 0.0;
          for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
            const int ip = ipobf_v[iob];
            const double wgt = ipobf_wgt[iob];
            T_wall += wgt*Tp[ip];
          }
          const int icv = zone_ptr->cvobf[ibf];
          // recall the definition of k_eff:
          //k_eff[ibf] = q_wall[ibf]*zone_ptr->area_bf[ibf]/dT;
          tmp[ibf] = k_eff[ibf]*(solver->T[icv] - T_wall)/zone_ptr->area_bf[ibf];
        }
        else {
          tmp[ibf] = 0.0;
        }
      }
    }

    return CTI_DATA_OK;

  }
  else if (name == zone_ptr->getName()+":T_wall") {

    if (args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double * tmp = zone_ptr->createBfD1Data(v);

    if (b_eval_func) {
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
        // reconstruct wall temperature from Tp...
        double T_wall = 0.0;
        for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
          const int ip = ipobf_v[iob];
          const double wgt = ipobf_wgt[iob];
          T_wall += wgt*Tp[ip];
        }
        tmp[ibf] = T_wall;
      }
    }

    return CTI_DATA_OK;

  }

  return PremixedBc::funcEvalCtiData(v,name,args,b_eval_func);

}

void WmAlgChtP::query(const string& param_str) const {

  // these statically sized buffers have room for up to 3 layers of temperature
  // data. Expand them if you implement more...

  double my_buf[9] = { 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0 };
  double my_minmax[10] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL,
                           HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL,
                           HUGE_VAL, HUGE_VAL };

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    if (zone_ptr->area_over_delta_bf[ibf] > 0.0) {

      const int icv = zone_ptr->cvobf[ibf];
      const double y1    = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
      const double u_tau = sqrt( tau_wall[ibf]/solver->rho[icv]);
      const double nu    = solver->mu_lam[icv]/solver->rho[icv];
      const double y_plus = y1*u_tau/nu;

      my_buf[0] += zone_ptr->area_bf[ibf];
      my_buf[1] += zone_ptr->area_bf[ibf]*tau_wall[ibf];
      my_buf[2] += zone_ptr->area_bf[ibf]*y_plus;
      my_buf[3] += zone_ptr->area_bf[ibf]*q_wall[ibf];

      double T_wall = 0.0;
      for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
        const int ip = ipobf_v[iob];
        const double wgt = ipobf_wgt[iob];
        T_wall += wgt*Tp[ip];
      }
      my_buf[4] += zone_ptr->area_bf[ibf]*T_wall;
      my_minmax[0] = min(my_minmax[0],T_wall);
      my_minmax[1] = min(my_minmax[1],-T_wall);

      my_buf[5] += zone_ptr->area_bf[ibf]*solver->T[icv];
      my_minmax[2] = min(my_minmax[2],solver->T[icv]);
      my_minmax[3] = min(my_minmax[3],-solver->T[icv]);

      for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
        my_buf[6+ilayer] += zone_ptr->area_bf[ibf]*T_layer[ilayer][ibf];
        my_minmax[4+ilayer*2] = min(my_minmax[4+ilayer*2],T_layer[ilayer][ibf]);
        my_minmax[4+ilayer*2+1] = min(my_minmax[4+ilayer*2+1],-T_layer[ilayer][ibf]);
      }

    }
  }

  double buf[9];
  MPI_Reduce(my_buf,buf,9,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  double minmax[10];
  MPI_Reduce(my_minmax,minmax,10,MPI_DOUBLE,MPI_MIN,0,mpi_comm);

  if ( mpi_rank == 0 ) {

    cout << " > " << zone_ptr->getName() << ": time, int(tau_wall)dA, int(q_w)dA, avg(q_w), avg(y_plus): "
         << solver->time << " "
         << buf[1] << " "
         << buf[3] << " "
         << buf[3]/buf[0] << " "
         << buf[2]/buf[0] << " T_wall(min,avg,max): "
         << minmax[0] << " " << buf[4]/buf[0] << " " << -minmax[1];

    for (int ilayer = 0; ilayer < nlayer; ++ilayer)
      cout << " T_layer" << ilayer+1 << ": " << minmax[4+ilayer*2] << " "
           << buf[6+ilayer]/buf[0] << " " << -minmax[4+ilayer*2+1];

    cout << " T[cvobf]: " << minmax[2] << " "
         << buf[5]/buf[0] << " " << -minmax[3] << endl;

  }

}

void WallConductiveP::initData() {
  Param* param = getParam(getName());

  int i_param = 1;
  T_end_bc = param->getDouble(i_param++);

  if (mpi_rank == 0)
    cout <<
      "WallConductiveP::initData() for bc: " << getName() <<
      "\n > T=" << T_end_bc << endl;

  // layers should be inside to out (typically hot to cold)...

  thermal_resistance_sum = 0.0;
  while (param->getString(i_param) != "END") {
    layer_length_bc.push_back(param->getDouble(i_param++));
    thermal_conductivity_bc.push_back(param->getDouble(i_param++));
    if (mpi_rank == 0)
      cout << " > layer " << layer_length_bc.size() << " thickness=" << layer_length_bc.back() << " k=" << thermal_conductivity_bc.back() << endl;
    // XXX Assuming layers are not parallel.
    thermal_resistance_sum += layer_length_bc.back() / thermal_conductivity_bc.back();
  }
  assert(layer_length_bc.size() == thermal_conductivity_bc.size());

}

void WallConductiveP::query(const string& param_str) const {

  // we want to report average layer temperatures as well...

  const int nlayers = layer_length_bc.size();
  double * my_buf = new double[4+nlayers];
  for (int i = 0; i < 4+nlayers; ++i)
    my_buf[i] = 0.0;

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    // only consider non-zero-area faces...
    if ((zone_ptr->area_over_delta_bf[ibf] > 0.0)&&(zone_ptr->area_bf[ibf] > 0.0)) {

      const int icv             = zone_ptr->cvobf[ibf];

      // here we DO NOT include the sgs model in the closures...
      const double y1           = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
      const double visc_coeff   = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
      const double cp           = solver->R[icv]*solver->gamma[icv]/(solver->gamma[icv]-1.0);
      const double k_cv         = solver->loc_lam[icv]*cp;

      double du_mag             = 0.0;
      FOR_I3 du_mag            += (solver->u[icv][i] - u_bc[i])*(solver->u[icv][i] - u_bc[i]);

      const double tau_wall     = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
      const double y_plus       = y1*sqrt(tau_wall*solver->rho[icv])/solver->mu_lam[icv];
      double thermal_resistance = thermal_resistance_sum + y1 / k_cv;
      const double q_wall       = (solver->T[icv] - T_end_bc) / thermal_resistance;

      my_buf[0]                += zone_ptr->area_bf[ibf];
      my_buf[1]                += zone_ptr->area_bf[ibf]*tau_wall;
      my_buf[2]                += zone_ptr->area_bf[ibf]*y_plus;
      my_buf[3]                += zone_ptr->area_bf[ibf]*q_wall;

      // also, report the layer temperatures...

      thermal_resistance        = y1 / k_cv;
      for (int i = 0; i < nlayers; ++i) {
        const double T_layer    = solver->T[icv] - q_wall*thermal_resistance;
        my_buf[4+i]            += zone_ptr->area_bf[ibf]*T_layer;
        thermal_resistance     += layer_length_bc[i] / thermal_conductivity_bc[i];
      }

    }

  }

  double * buf = NULL;
  if (mpi_rank == 0) buf = new double[4+nlayers];

  MPI_Reduce(my_buf,buf,4+nlayers,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  delete[] my_buf;

  if ( mpi_rank == 0 ) {

    cout << " > " << zone_ptr->getName() << ": time, int(tau_wall)dA, int(q_wall)dA, avg(y_plus): "
         << solver->time << " "
         << buf[1] << " "
         << buf[3] << " "
         << buf[2]/buf[0] << " avg layer T: ";
    for (int i = 0; i < nlayers; ++i)
      cout << buf[4+i]/buf[0] << " ";
    cout << endl;

    delete[] buf;

  }

}

void WallConductiveP::addBoundaryFlux(PremixedRhs* rhs) const {
  assert( solver != NULL);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    // only consider non-zero-area faces...
    if ((zone_ptr->area_over_delta_bf[ibf] > 0.0)&&(zone_ptr->area_bf[ibf] > 0.0)) {

      const int icv = zone_ptr->cvobf[ibf];

      // here we DO NOT include the sgs model in the closures...
      const double y1                 = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
      const double visc_coeff         = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
      const double cp                 = solver->R[icv]*solver->gamma[icv]/(solver->gamma[icv]-1.0);
      const double k_cv               = solver->loc_lam[icv]*cp;

      const double thermal_resistance = thermal_resistance_sum + y1 / k_cv;
      //const double q_wall      = (solver->T[icv] - T_end_bc) / thermal_resistance;

      for (int i = 0; i < 3; ++i)
        rhs[icv].rhou[i] -= solver->p[icv] * zone_ptr->n_bf[ibf][i] + visc_coeff * (solver->u[icv][i] - u_bc[i]);

      // XXX if u_bc \neq 0, there should be a u \dot \tau term added here..
      rhs[icv].rhoE -= (solver->T[icv] - T_end_bc) / thermal_resistance * zone_ptr->area_bf[ibf];
    }

  }

}

CtiRegister::CtiDataError WallConductiveP::funcEvalCtiData(CtiRegister::CtiData& v,
							   const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if (zone_ptr->isLeadingStrMatched(name)) {

    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix. We have to check for the actual data now...

    const string tau_wall_str = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str       = zone_ptr->getName() + ":" + "y_plus";
    const string q_wall_str   = zone_ptr->getName() + ":" + "q_wall";
    //const string T1_str     = zone_ptr->getName() + ":" + "T1";
    //const string T2_str     = zone_ptr->getName() + ":" + "T2";

    // all of the following functions do not take any arguments -- check first

    if (name == tau_wall_str) {

      double * tmp = zone_ptr->createBfD1Data(v);

      // if stats have been requested on this quantity, this function will be called
      // but its possible that the mu_lam, rho necessary have not been populated yet.
      // check the data flags for these values -- revisit if this impacts performance
      // during the solution run-time.

      if (b_eval_func) {
        if (args.size() == 0) {
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            const int icv           = zone_ptr->cvobf[ibf];
            const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
            double du_mag           = 0.0;
            for (int i = 0; i < 3; ++i)
              du_mag               += (solver->u[icv][i] - u_bc[i])*(solver->u[icv][i] - u_bc[i]);

            tmp[ibf]                = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
          }
        }
        else if (args.size() == 1) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          double dir[3] = {0.0,0.0,0.0};

          if (arg->getType() == D3_DATA) {
            FOR_I3 dir[i] = arg->d3(i);
            NORMALIZE(dir);
          }
          else if (arg->getType() == D_DATA) {
            const int index = int(arg->d());
            if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;
            dir[index] = 1.0;
          }
          else return CTI_DATA_NOT_VALID;

          // loop faces and compute
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            const int icv           = zone_ptr->cvobf[ibf];
            const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
            double local_u[3]       = DIFF(solver->u[icv],u_bc);
            tmp[ibf]                = visc_coeff*DOT_PRODUCT(dir,local_u)/zone_ptr->area_bf[ibf];
          }
        }
        else return CTI_DATA_ARG_COUNT;
      }

      return CTI_DATA_OK;

    }
    else if (name == yp_str) {

      if (args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          const int icv           = zone_ptr->cvobf[ibf];
          const double y1         = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
          const double visc_coeff = solver->mu_lam[icv]/y1;

          double du_mag = 0.0;
          for (int i = 0; i < 3; ++i)
            du_mag += (solver->u[icv][i] - u_bc[i])*(solver->u[icv][i] - u_bc[i]);

          const double tau_wall   = visc_coeff*sqrt(du_mag);
          tmp[ibf]                = y1*sqrt(tau_wall*solver->rho[icv])/solver->mu_lam[icv];

        }

      }

      return CTI_DATA_OK;

    }

    else if (name == q_wall_str) {

      if (args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          // only consider non-zero-area faces...
          if ((zone_ptr->area_over_delta_bf[ibf] > 0.0) && (zone_ptr->area_bf[ibf] > 0.0)) {
            const int icv                   = zone_ptr->cvobf[ibf];
            // here we DO NOT include the sgs model in the closures...
            const double y1                 = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
            const double cp                 = solver->R[icv]*solver->gamma[icv]/(solver->gamma[icv] - 1.0);
            const double k_cv               = solver->loc_lam[icv]*cp;
            const double thermal_resistance = thermal_resistance_sum + y1 / k_cv;
            tmp[ibf]                        = (solver->T[icv] - T_end_bc)/thermal_resistance;
          }
          else
            tmp[ibf] = 0.0;
        }
      }
      return CTI_DATA_OK;
    }
  }

  return PremixedBc::funcEvalCtiData(v,name,args,b_eval_func);

}

void WmAlgConductiveP::initData() {
  Param * param = getParam(getName());
  int i_param   = 1;
  T_end_bc      = param->getDouble(i_param++);

  // Read layer lengths and conductivities
  thermal_resistance_external_times_area = 0.0;
  while (param->getString(i_param) != "END") {
    layer_length_bc.push_back(param->getDouble(i_param++)); // [m]
    thermal_conductivity_bc.push_back(param->getDouble(i_param++)); // [W/m*K]
    thermal_resistance_external_times_area += layer_length_bc.back() / thermal_conductivity_bc.back(); // [m] / [W/m*K] = [W*m^2/K]
  }

  // layer length and conductivity vectors are the same size!
  assert(layer_length_bc.size() == thermal_conductivity_bc.size());

  assert(tau_wall == NULL); tau_wall  = new double[zone_ptr->nbf];
  assert(u1       == NULL); u1        = new double[zone_ptr->nbf];
  assert(q_wall   == NULL); q_wall    = new double[zone_ptr->nbf];
  assert(T_bc     == NULL); T_bc      = new double[zone_ptr->nbf];
}

void WmAlgConductiveP::initialHook() {
  if (!zone_ptr->checkDataFlag("u1")) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv = zone_ptr->cvobf[ibf];
      u1[ibf] = MAG(solver->u[icv]);

      // set an initial value for tau_wall, this will constitute the
      // initial guess for the solution when it is finally constructed

      if (zone_ptr->area_over_delta_bf[ibf] > 0.0) {
        const double y1 = zone_ptr->area_bf[ibf] /zone_ptr->area_over_delta_bf[ibf];
        tau_wall[ibf]   = solver->mu_lam[icv] * MAG(solver->u[icv])/y1;
        q_wall[ibf]     = 0.0;
        T_bc[ibf]       = 0.0;
      }
      else {
        tau_wall[ibf]   = 0.0;
        q_wall[ibf]     = 0.0;
        T_bc[ibf]       = 0.0;
      }
    }
  }
}

int WmAlgConductiveP::storeData(double * buf) const {

  // the wall model needs to store tau_wall because it is used as
  // an initial guess. No need to store u1. It is not used...
  // q_wall is computed

  if (buf) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      buf[ibf] = tau_wall[ibf];
    }
  }

  return zone_ptr->nbf;

}

int WmAlgConductiveP::restoreData(double * buf) {

  assert(buf);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    tau_wall[ibf] = buf[ibf];
  }

  return zone_ptr->nbf;

}

WmAlgConductiveP::~WmAlgConductiveP() {
  DELETE(u1);
  DELETE(tau_wall);
  DELETE(q_wall);
  DELETE(T_bc);
}

void WmAlgConductiveP::calcRhs(const double time, const int rk_stage) {
  if (rk_stage == 1) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv = zone_ptr->cvobf[ibf];
      // Ensure that the cell actually has area
      if (zone_ptr->area_over_delta_bf[ibf] > 0.0) {
        const double y1     = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        const double mag_u  = MAG(solver->u[icv]);

        u1[ibf]             = mag_u;

        // use the existing tau_wall to set a guess for the u_tau
        const double u_tau = sqrt(tau_wall[ibf]/solver->rho[icv]);
        tau_wall[ibf] = AlgebraicWM::solve_tau(u1[ibf], y1, solver->rho[icv], solver->mu_lam[icv], u_tau);

        const double cp = solver->R[icv]*solver->gamma[icv]/(solver->gamma[icv]-1.0);
        const double thermal_resistance_cv = AlgebraicWM::compute_R_approx(tau_wall[ibf], solver->rho[icv], solver->mu_lam[icv], solver->loc_lam[icv], cp, y1); // [K/W]
        const double thermal_resistance_overall = thermal_resistance_external_times_area/zone_ptr->area_bf[ibf] + thermal_resistance_cv; // [K/W]

        q_wall[ibf] = (solver->T[icv] - T_end_bc)/thermal_resistance_overall/zone_ptr->area_bf[ibf]; // [K]/[K/W]/[m^2] = [W/m^2]
        T_bc[ibf] = solver->T[icv] - q_wall[ibf]*thermal_resistance_cv*zone_ptr->area_bf[ibf]; // [K] - [W/m^2]*[K/W]*[m^2] = [K] - [K]
      }
      else {
        tau_wall[ibf] = 0.0;
        q_wall[ibf]   = 0.0;
        T_bc[ibf]     = solver->T[icv];
      }
    }
  }
}

void WmAlgConductiveP::addBoundaryFlux(NonpremixedRhs * rhs) const {

  assert(solver != NULL);
  for (int ibf = 0; ibf < zone_ptr->nbf; ibf++) {

    const int icv = zone_ptr->cvobf[ibf];

    // tau_wall has already been computed at this point, but we'll assume
    // that the wall stress is instantaneously aligned with the velocity vector

    double u_mag = MAG(solver->u[icv]);

    if (u_mag > 0.0) {
      for (int i = 0; i < 3; ++i) {
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];
        rhs[icv].rhou[i] -= tau_wall[ibf]*solver->u[icv][i]/u_mag*zone_ptr->area_bf[ibf];
      }
      // XXX if u_bc \neq 0, there should be a u \dot \tau term added here...
      rhs[icv].rhoE -= q_wall[ibf]*zone_ptr->area_bf[ibf];
    }
    else {
      for (int i = 0; i < 3; ++i)
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];

      // XXX if u_bc \neq 0, there should be a u \dot \tau term added here...
      rhs[icv].rhoE -= q_wall[ibf]*zone_ptr->area_bf[ibf];
    }
  }
}

CtiRegister::CtiDataError WmAlgConductiveP::funcEvalCtiData(CtiRegister::CtiData& v,
                                                            const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if (zone_ptr->isLeadingStrMatched(name)) {

    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix. We have to check for the actual data now...

    // for consistency, we are supporting a tau_wall(), y_plus(), q_wall()
    // function to have the same syntax as all the other wall zones...

    const string tau_wall_str = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str       = zone_ptr->getName() + ":" + "y_plus";
    const string qw_str       = zone_ptr->getName() + ":" + "q_wall";


    if (name == tau_wall_str) {
      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {
        if (args.size() == 0)
          // use what was computed by wall model
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
            tmp[ibf] = tau_wall[ibf];
        else if (args.size() == 1) {
          list <CtiRegister::CtiData>::iterator arg = args.begin();

          double dir[3] = {0.0, 0.0, 0.0};

          if (arg->getType() == D3_DATA) {
            FOR_I3 dir[i] = arg->d3(i);
            NORMALIZE(dir);
          }
          else if (arg->getType() == D_DATA) {
            const int index = int(arg->d());
            if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;
            dir[index] = 1.0;
          }
          else return CTI_DATA_NOT_VALID;

          // loop faces and compute
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
            const int icv       = zone_ptr->cvobf[ibf];
            const double mag_u  = MAG(solver->u[icv]);

            double unit_u[3];
            if (mag_u > 0.0) {
              for (int i = 0; i < 3; ++i) unit_u[i] = solver->u[icv][i]/mag_u;
              tmp[ibf] = tau_wall[ibf]*DOT_PRODUCT(dir,unit_u);
            }
            else
              tmp[ibf] = 0.0;
          }
        }
        else return CTI_DATA_ARG_COUNT;
      }

      return CTI_DATA_OK;
    }

    else if (name == yp_str) {

      if (args.size() != 0) return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          const int icv   = zone_ptr->cvobf[ibf];
          const double y1 = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
          tmp[ibf]        = y1*sqrt(solver->rho[icv]*tau_wall[ibf])/solver->mu_lam[icv];
        }
      }

      return CTI_DATA_OK;

    }

    else if (name == qw_str) {

      if (args.size() != 0) return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func)
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
          tmp[ibf] = q_wall[ibf];

      return CTI_DATA_OK;

    }

  }

  return PremixedBc::funcEvalCtiData(v,name,args,b_eval_func);
}

void WmAlgConductiveP::query(const string& param_str) const {

  // report an estimate of the local y^+ on the boundary and the integrated tau_wall and q_wall

  double my_buf[4] = {0.0, 0.0, 0.0, 0.0};

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    if (zone_ptr->area_over_delta_bf[ibf] > 0.0) {

      const double y1 = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
      const double u_tau = sqrt(tau_wall[ibf]/solver->rho[icv]);
      const double nu     = solver->mu_lam[icv]/solver->rho[icv];
      const double y_plus = y1*u_tau/nu;

      my_buf[0] += zone_ptr->area_bf[ibf];
      my_buf[1] += zone_ptr->area_bf[ibf]*tau_wall[ibf];
      my_buf[2] += zone_ptr->area_bf[ibf]*y_plus;
      my_buf[3] += zone_ptr->area_bf[ibf]*q_wall[ibf];
    }
  }

  double buf[4];
  MPI_Reduce(my_buf, buf, 4, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);

  if (mpi_rank == 0)
    cout << " > "  << zone_ptr->getName() << ": time, int(tu_wall)dA, int(q_wall)dA, avg(y_plus):"
         << solver->time << "\t" << buf[1]        << "\t"
         << buf[3]       << "\t" << buf[2]/buf[0] << endl;
}

void CbcP::addBoundaryFlux(PremixedRhs* rhs) const {

  assert( solver != NULL);
  assert ( mf != NULL);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    // TODO: XXX assuming that dgamma/dn is weakly varying here;
    // update with the appropriate bc gamma...
    PremixedRhs flux;
    calcRiemannFlux(flux,zone_ptr->n_bf[ibf],solver->cv_light[icv],
                    bf_state[ibf],solver->gamma[icv], solver->gamma[icv]);

    //  mass flux needs to be stored for potential passive scalars..

    mf[ibf]       = flux.rho;

    rhs[icv].rho -= flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv].rhou[i] -= flux.rhou[i];
    rhs[icv].rhoE -= flux.rhoE;
    rhs[icv].rhoZ -= flux.rhoZ;
    rhs[icv].rhoC -= flux.rhoC;
  }
}

void CbcP::addScalarBoundaryFlux(double * rhs_sc) const {
  BcTemplate<PremixedSolver,PremixedRhs>::_addScalarBoundaryFlux(rhs_sc);
}

void CbcP::query(const string& param_str) const {

  // report integrated values for the p_bc, rho_bc, u_bc, Z_bc, C_bc
  double my_buf[7] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double buf[7];

  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv         = zone_ptr->cvobf[ibf];
    const double area     = MAG(zone_ptr->n_bf[ibf]);
    const double undA_bc  = DOT_PRODUCT(zone_ptr->n_bf[ibf], solver->u[icv]);
    my_buf[0]            += area;
    my_buf[1]            += area*solver->rho[icv];
    my_buf[2]            += solver->rho[icv]*undA_bc;
    my_buf[3]            += area*solver->p[icv];
    my_buf[4]            += area*solver->Z[icv];
    my_buf[5]            += area*solver->C[icv];
    my_buf[6]            += area*solver->T[icv];
  }//FOR_IBF

  MPI_Reduce(my_buf,buf,7,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  if ( mpi_rank == 0 ) {
    cout << " > " << zone_ptr->getName() << ": time, rho_bc, mdot, p_bc, Z_bc, C_bc, T_bc: "
         << solver->time << "   " << buf[1]/buf[0] << "    "
         << buf[2] << "    " << buf[3]/buf[0] << "    "
         << buf[4]/buf[0] << "    " << buf[5]/buf[0] << "    "
         << buf[6]/buf[0] << endl;
  }

}

void computeBfWallDistance(double * wall_dist, const BfZone * zone_ptr, const PremixedSolver* solver) {

  const int ncv = solver->ncv;
  double * cv_flag = new double[ncv];
  FOR_ICV cv_flag[icv] = -2.0;

  const int nbf = zone_ptr->nbf;
  // flag cells adjacent to faces of interest
  FOR_IBF {
    cv_flag[zone_ptr->cvobf[ibf]] = -1.0;
  }

  // We have flagged the cvs adjacent to this zone and we can loop over other boundary conditions to find cells at the corner of this zone and wall BCs. In only these cells we will compute the wall distance
  for(vector<PremixedBc*>::const_iterator it = solver->bcs.begin(); it != solver->bcs.end(); ++it) {
    string s=(*it)->getName();
    Param* p=getParam(s);
    if ((p->getString()=="WALL_ADIABATIC") || (p->getString()=="WALL_ISOTHERMAL") || (p->getString()=="WM_ALG_ADIABATIC") || (p->getString()=="WM_ALG_ISOTHERMAL") || (p->getString()=="WM_ALG_CHT") || (p->getString()=="WM_ALG_CONDUCTIVE") || (p->getString()=="WALL_CHT") || (p->getString()=="WALL_CONDUCTIVE")) {
      // TODO This must recognize all wall bcs. How?
      BfZone* zone_ptr_it=(*it)->zone_ptr;
      for (int ibf=0; ibf < zone_ptr_it->nbf; ++ibf) {
        const int icv = zone_ptr_it->cvobf[ibf];

        if (cv_flag[icv] == -1.0) {
          // first time hit, identify as a corner cell
          cv_flag[icv] = HUGE_VAL;
        }

        if (cv_flag[icv] > 0.0) {
          // we are at a zone corner cv, so compute wall-distance and store minimum in cv_flag
          const double dist = DIST(zone_ptr_it->x_bf[ibf],solver->x_cv[icv]);
          if (dist < cv_flag[icv]) {
            cv_flag[icv] = dist;
          }
        }
      }
    }
  }  // bcs iterator

  solver->computeBfDistanceFromFlaggedCvs(wall_dist,zone_ptr,cv_flag);
  DELETE(cv_flag);
}


void CbcProfileP::initData() {
  Param * param  = getParam(getName());
  assert( param != NULL);

  int type = -1;
  double type_double[2] = {0.0,0.0};

  FluentBPReader profile;
  int iarg = 1;
  while (iarg < param->size()) {
    string token = param->getString(iarg++);
    if (token == "FILE") {
      profile.init(param->getString(iarg++));
    }
    else if (token == "TYPE") {
      const string _type = param->getString(iarg++);
      if (_type == "UPTS") {
        type = 0;
      }
      else if (_type == "RUPS") {
        type = 1;
      }
      else if (_type == "MPTS") {
        type = 2;
        type_double[0] = param->getDouble(iarg++);
      }
      else {
        CERR("unsupported CBC_PROFILE TYPE " << _type << "; please choose from UPTS, RUPS, MPTS");
      }
    }
    else if (token == "FROM_MEAN") {
      profile.setUseMean(true);
    }
  }

  if (!profile.isInitialized()) {
    CERR("boundary profile was not properly initialized");
  }
  // profile.printVarsInFile();

  if (profile.getType() == LINE_WD) {
    // need to pass wall-distance for each ibf
    double * wall_dist = new double[zone_ptr->nbf];

    // compute wall distance
    computeBfWallDistance(wall_dist,zone_ptr,solver);

    profile.setPoints(wall_dist,zone_ptr->nbf);
    DELETE(wall_dist);
  }
  else {
    // pass in actual ibf locations
    double (*x_fa_zone)[3] = new double[zone_ptr->nbf][3];

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      FOR_I3 x_fa_zone[ibf][i] = zone_ptr->x_bf[ibf][i];
    }
    profile.setPoints(x_fa_zone,zone_ptr->x_global,zone_ptr->nbf);
    DELETE(x_fa_zone);
  }

  if (type == -1) {
    CERR("CBC_PROFILE \"TYPE\" must be specified");
  }
  else if (type == 0) {
    profile.ensureVar("x-velocity");
    profile.ensureVar("y-velocity");
    profile.ensureVar("z-velocity");
    profile.ensureVar("absolute-pressure");
    profile.ensureVar("temperature");
    profile.ensureVar("mixture-fraction");
    profile.ensureVar("progress-variable");

    assert( bf_state == NULL);
    bf_state = new PremixedState[zone_ptr->nbf];

    for (int ibf = 0, nbf=zone_ptr->nbf; ibf < nbf; ibf++) {
      bf_state[ibf].Z    = profile.getData(ibf,"mixture-fraction");
      bf_state[ibf].C    = profile.getData(ibf,"progress-variable");
      const double T_bc  = profile.getData(ibf,"temperature");
      bf_state[ibf].p    = profile.getData(ibf,"absolute-pressure");

      double R_bc, T0, e0, gamma0, a_gam;
      solver->chemtable->lookup(&R_bc, "R", &T0, "T", &e0, "e", &gamma0, "gamma",
                                &a_gam, "a_gamma", &bf_state[ibf].Z, &bf_state[ibf].C, 1);
      const double e_bc   = calcInternalEnergyFromTS(T_bc,T0,e0,R_bc,gamma0,a_gam);
      const double rho_bc = bf_state[ibf].p/R_bc/T_bc;

      bf_state[ibf].u[0] = profile.getData(ibf,"x-velocity");
      bf_state[ibf].u[1] = profile.getData(ibf,"y-velocity");
      bf_state[ibf].u[2] = profile.getData(ibf,"z-velocity");
      bf_state[ibf].sp_vol = 1.0/rho_bc;
      bf_state[ibf].h      = e_bc + bf_state[ibf].p*bf_state[ibf].sp_vol;
    }
  }
  else if (type == 1) {
    profile.ensureVar("x-velocity");
    profile.ensureVar("y-velocity");
    profile.ensureVar("z-velocity");
    profile.ensureVar("absolute-pressure");
    profile.ensureVar("density");
    profile.ensureVar("mixture-fraction");
    profile.ensureVar("progress-variable");

    assert( bf_state == NULL);
    bf_state = new NonpremixedState[zone_ptr->nbf];

    for (int ibf = 0, nbf=zone_ptr->nbf; ibf < nbf; ibf++) {
      bf_state[ibf].Z    = profile.getData(ibf,"mixture-fraction");
      bf_state[ibf].C    = profile.getData(ibf,"progress-variable");
      bf_state[ibf].p    = profile.getData(ibf,"absolute-pressure");
      bf_state[ibf].sp_vol = 1.0/profile.getData(ibf,"density");

      // get baseline R, T, e, gamma, a_gam at the bc mixture
      double R_bc, T0, e0, gamma0, a_gam;
      solver->chemtable->lookup(&R_bc, "R", &T0, "T", &e0, "e", &gamma0, "gamma",
                                &a_gam, "a_gamma", &bf_state[ibf].Z, &bf_state[ibf].C, 1);
      const double T_bc = bf_state[ibf].p/R_bc*bf_state[ibf].sp_vol;
      const double e_bc = calcInternalEnergyFromTS(T_bc,T0,e0,R_bc,gamma0,a_gam);

      bf_state[ibf].u[0] = profile.getData(ibf,"x-velocity");
      bf_state[ibf].u[1] = profile.getData(ibf,"y-velocity");
      bf_state[ibf].u[2] = profile.getData(ibf,"z-velocity");
      bf_state[ibf].h      = e_bc + bf_state[ibf].p*bf_state[ibf].sp_vol;
    }
  }
  else if (type == 2) {
    profile.ensureVar("absolute-pressure");
    profile.ensureVar("temperature");
    profile.ensureVar("mixture-fraction");
    profile.ensureVar("progress-variable");
    if (type_double[0] == 0.0) CERR("Mass flux cannot be zero");
    const double mdot_bc = type_double[0];

    assert( bf_state == NULL);
    bf_state = new NonpremixedState[zone_ptr->nbf];

    for (int ibf = 0, nbf=zone_ptr->nbf; ibf < nbf; ibf++) {
      bf_state[ibf].Z    = profile.getData(ibf,"mixture-fraction");
      bf_state[ibf].C    = profile.getData(ibf,"progress-variable");
      const double T_bc  = profile.getData(ibf,"temperature");
      bf_state[ibf].p    = profile.getData(ibf,"absolute-pressure");

      // the bf zone has stored an area_global which has the total zone area..
      // get baseline R, T, e, gamma, a_gam at the bc mixture
      double R_bc, T0, e0, gamma0, a_gam;
      solver->chemtable->lookup(&R_bc, "R", &T0, "T", &e0, "e", &gamma0, "gamma",
                                &a_gam, "a_gamma",&bf_state[ibf].Z,&bf_state[ibf].C,1);
      const double e_bc   = calcInternalEnergyFromTS(T_bc,T0,e0,R_bc,gamma0,a_gam);
      const double rho_bc = bf_state[ibf].p/R_bc/T_bc;
      const double un_bc  = -mdot_bc/(rho_bc*zone_ptr->area_global);

      const double bf_area = MAG(zone_ptr->n_bf[ibf]);
      for (int i =0; i < 3; ++i)
        bf_state[ibf].u[i] = un_bc*zone_ptr->n_bf[ibf][i]/bf_area;
      bf_state[ibf].sp_vol = 1.0/rho_bc;
      bf_state[ibf].h      = e_bc + bf_state[ibf].p/rho_bc;
    }
  }
  else CERR("unrecognized type somehow set...");

  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<PremixedSolver,PremixedRhs>::parseScalarBc(param);

  // profile.dumpTecplot("profile_test_dump.dat");
}

void CbcUptsP::initData() {

  double u_bc[3], p_bc, T_bc, Z_bc, C_bc, rho_bc, e_bc;
  Param * param = getParam(getName());
  u_bc[0] = param->getDouble(1);
  u_bc[1] = param->getDouble(2);
  u_bc[2] = param->getDouble(3);
  p_bc    = param->getDouble(4);
  T_bc    = param->getDouble(5);
  Z_bc    = param->getDouble(6);
  C_bc    = param->getDouble(7);

  // get baseline R, T, e, gamma, a_gam at the bc mixture
  // currently uses the chemtable from the solver

  double R_bc, T0, e0, gamma0, a_gam;
  solver->chemtable->lookup(&R_bc, "R", &T0, "T", &e0, "e", &gamma0, "gamma",
                            &a_gam, "a_gamma",&Z_bc,&C_bc,1);
  e_bc   = calcInternalEnergyFromTS(T_bc,T0,e0,R_bc,gamma0,a_gam);
  rho_bc = p_bc/R_bc/T_bc;

  assert( bf_state == NULL);
  bf_state = new PremixedState[zone_ptr->nbf];
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    for (int i =0; i < 3; ++i)
      bf_state[ibf].u[i] = u_bc[i];
    bf_state[ibf].sp_vol = 1.0/rho_bc;
    bf_state[ibf].p      = p_bc;
    bf_state[ibf].h      = e_bc + p_bc/rho_bc;
    bf_state[ibf].Z      = Z_bc;
    bf_state[ibf].C      = C_bc;
  }

  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<PremixedSolver,PremixedRhs>::parseScalarBc(param);

}

void CbcUptsSpongeP::initData() {

  // sponge parsing...
  // <zone> CBC_UPTS_SPONGE <ux> <uy> <uz> <p> <T> <Z> <C> SPONGE_[XYZ] <x0> <x1> <strength>

  double u_bc[3], p_bc, T_bc, Z_bc, C_bc, rho_bc, e_bc;
  Param * param = getParam(getName());
  u_bc[0] = param->getDouble(1);
  u_bc[1] = param->getDouble(2);
  u_bc[2] = param->getDouble(3);
  p_bc    = param->getDouble(4);
  T_bc    = param->getDouble(5);
  Z_bc    = param->getDouble(6);
  C_bc    = param->getDouble(7);

  string sdir;
  string sponge_dir = param->getString(8);
  if        (sponge_dir == "SPONGE_X") {
    idir = 0;
    sdir = "X";
  } else if (sponge_dir == "SPONGE_Y") {
    idir = 1;
    sdir = "Y";
  } else if (sponge_dir == "SPONGE_Z") {
    idir = 2;
    sdir = "Z";
  } else {
    CERR( " > Invalid SPONGE direction... expecting: SPONGE_[XYZ], not: " << sponge_dir);
  }

  x0 = param->getDouble(9);
  x1 = param->getDouble(10);
  sponge_strength = param->getDouble(11);

  if (mpi_rank == 0) {
   cout << " CbcUptsSpongeP::initData()..." << endl
   << " > zone = " << param->getName() << endl
   << " > u_bc = " << u_bc[0] << ", " << u_bc[1] << ", " << u_bc[2] << endl
   << " > p_bc = " << p_bc << endl
   << " > T_bc = " << T_bc << endl
   << " > Z_bc = " << Z_bc << endl
   << " > C_bc = " << C_bc << endl
   << " > SPONGE_" << sdir << " : " << sdir << "0 = " << x0 << ", " << sdir << "1 = " << x1 << endl
   << " > STRENGTH = " << sponge_strength << endl;
  }

  // get baseline R, T, e, gamma, a_gam at the bc mixture
  // currently uses the chemtable from the solver

  double R_bc, T0, e0, gamma0, a_gam;
  solver->chemtable->lookup(&R_bc, "R", &T0, "T", &e0, "e", &gamma0, "gamma",
                            &a_gam, "a_gamma", &Z_bc, &C_bc, 1);
  e_bc   = calcInternalEnergyFromTS(T_bc,T0,e0,R_bc,gamma0,a_gam);
  rho_bc = p_bc/R_bc/T_bc;

  // set the sponge values...
  rho_sponge = rho_bc;
  FOR_I3 u_sponge[i] = u_bc[i];
  rhoE_sponge = rho_bc*(e_bc + 0.5*DOT_PRODUCT(u_bc,u_bc));
  Z_sponge = Z_bc;
  C_sponge = C_bc;

  FOR_I3 bf_state_sponge.u[i] = u_bc[i];
  bf_state_sponge.sp_vol      = 1.0/rho_bc;
  bf_state_sponge.p           = p_bc;
  bf_state_sponge.h           = e_bc + p_bc/rho_bc;
  bf_state_sponge.Z           = Z_bc;
  bf_state_sponge.C           = C_bc;

  assert(mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<PremixedSolver,PremixedRhs>::parseScalarBc(param);

}

void CbcUptsSpongeP::addBoundaryFlux(PremixedRhs* rhs) const {

  assert(solver != NULL);
  assert(mf != NULL);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    // TODO: XXX assuming that dgamma/dn is weakly varying here;
    // update with the appropriate bc gamma...
    PremixedRhs flux;
    calcRiemannFlux(flux,zone_ptr->n_bf[ibf],solver->cv_light[icv],
                    bf_state_sponge, solver->gamma[icv], solver->gamma[icv]);

    // mass flux needs to be stored for potential passive scalars...

    mf[ibf]        = flux.rho;

    rhs[icv].rho  -= flux.rho;
    for (int i = 0; i < 3; ++i)
      rhs[icv].rhou[i] -= flux.rhou[i];
    rhs[icv].rhoE -= flux.rhoE;
    rhs[icv].rhoZ -= flux.rhoZ;
    rhs[icv].rhoC -= flux.rhoC;

  }

  // and apply the sponge source terms to the volume in the sponge region...

  for (int icv = 0; icv < solver->ncv; ++icv) {
    const double x_sp = (solver->x_cv[icv][idir] - x0)/(x1 - x0);
    if ((x_sp > 0.0) && (x_sp <= 1.0)) {
      // sponge profile: use simple quadratic as in IdealGasSolver...
      // sponge strength should be ~ (5 to 15)*U/L or (5 to 15)*C/L
      const double sponge_coeff = sponge_strength*x_sp*x_sp;
      rhs[icv].rho             -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv] - rho_sponge);
      for (int i = 0; i < 3; ++i)
        rhs[icv].rhou[i]       -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv]*solver->u[icv][i] - rho_sponge*u_sponge[i]);
      rhs[icv].rhoE            -= solver->vol_cv[icv]*sponge_coeff*(solver->rhoE[icv] - rhoE_sponge);
      rhs[icv].rhoZ            -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv]*solver->Z[icv] - rho_sponge*Z_sponge);
      rhs[icv].rhoC            -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv]*solver->C[icv] - rho_sponge*C_sponge);
    }
  }

}

void CbcMptsP::initData() {

  double mdot_bc, p_bc, T_bc, Z_bc, C_bc,un_bc,rho_bc,e_bc;
  Param* param = getParam(getName());
  mdot_bc = param->getDouble(1);
  p_bc    = param->getDouble(2);
  T_bc    = param->getDouble(3);
  Z_bc    = param->getDouble(4);
  C_bc    = param->getDouble(5);

  // the bf zone has stored an area_global which has the total zone area..
  // get baseline R, T, e, gamma, a_gam at the bc mixture

  double R_bc, T0, e0, gamma0, a_gam;
  solver->chemtable->lookup(&R_bc, "R", &T0, "T", &e0, "e", &gamma0, "gamma",
                            &a_gam, "a_gamma",&Z_bc,&C_bc,1);
  e_bc   = calcInternalEnergyFromTS(T_bc,T0,e0,R_bc,gamma0,a_gam);
  rho_bc = p_bc/R_bc/T_bc;
  un_bc  = -mdot_bc/(rho_bc*zone_ptr->area_global);

  assert( bf_state == NULL);
  bf_state = new PremixedState[zone_ptr->nbf];
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const double bf_area = MAG(zone_ptr->n_bf[ibf]);
    for (int i =0; i < 3; ++i)
      bf_state[ibf].u[i] = un_bc*zone_ptr->n_bf[ibf][i]/bf_area;
    bf_state[ibf].sp_vol = 1.0/rho_bc;
    bf_state[ibf].p      = p_bc;
    bf_state[ibf].h      = e_bc + p_bc/rho_bc;
    bf_state[ibf].Z      = Z_bc;
    bf_state[ibf].C      = C_bc;
  }

  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<PremixedSolver,PremixedRhs>::parseScalarBc(param);

}

void CbcMsTableP::initData() {
  
  Param* param = getParam(getName()); assert(param);
  const double mdot_bc = param->getDouble(1);
  const double Z_bc    = param->getDouble(2);
  const double C_bc    = param->getDouble(3);

  // the bf zone has stored an area_global which has the total zone area..
  // get baseline R, T, e, gamma, a_gam at the bc mixture
  double R_bc,T_bc,e_bc,gamma_bc; 
  solver->chemtable->lookup(&R_bc,"R",&T_bc,"T",&e_bc,"e",&gamma_bc,"gamma",&Z_bc,&C_bc,1);
  const double p_bc = solver->chemtable->pressure;
  const double rho_bc = p_bc/(R_bc*T_bc);
  const double un_bc  = -mdot_bc/(rho_bc*zone_ptr->area_global);

  assert( bf_state == NULL);
  bf_state = new PremixedState[zone_ptr->nbf];
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const double bf_area = MAG(zone_ptr->n_bf[ibf]);
    for (int i =0; i < 3; ++i)
      bf_state[ibf].u[i] = un_bc*zone_ptr->n_bf[ibf][i]/bf_area;
    bf_state[ibf].sp_vol = 1.0/rho_bc;
    bf_state[ibf].p      = p_bc;
    bf_state[ibf].h      = e_bc + p_bc/rho_bc;
    bf_state[ibf].Z      = Z_bc;
    bf_state[ibf].C      = C_bc;
  }

  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<PremixedSolver,PremixedRhs>::parseScalarBc(param);

}

void CbcRunpsP::initData() {

  double rho_bc, un_bc, p_bc, Z_bc, C_bc, T_bc, e_bc;
  Param * param = getParam(getName());
  rho_bc  = param->getDouble(1);
  un_bc   = param->getDouble(2);
  p_bc    = param->getDouble(3);
  Z_bc    = param->getDouble(4);
  C_bc    = param->getDouble(5);

  // get baseline R, T, e, gamma, a_gam at the bc mixture

  double R_bc, T0, e0, gamma0, a_gam;
  solver->chemtable->lookup(&R_bc, "R", &T0, "T", &e0, "e", &gamma0, "gamma",
                            &a_gam, "a_gamma",&Z_bc,&C_bc,1);
  T_bc = p_bc/R_bc/rho_bc;
  e_bc = calcInternalEnergyFromTS(T_bc,T0,e0,R_bc,gamma0,a_gam);

  assert( bf_state == NULL);
  bf_state = new PremixedState[zone_ptr->nbf];
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const double bf_area = MAG(zone_ptr->n_bf[ibf]);
    for (int i =0; i < 3; ++i)
      bf_state[ibf].u[i] = -un_bc*zone_ptr->n_bf[ibf][i]/bf_area;
    bf_state[ibf].sp_vol = 1.0/rho_bc;
    bf_state[ibf].p      = p_bc;
    bf_state[ibf].h      = e_bc + p_bc/rho_bc;
    bf_state[ibf].Z      = Z_bc;
    bf_state[ibf].C      = C_bc;
  }

  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<PremixedSolver,PremixedRhs>::parseScalarBc(param);

}

void CbcRupsP::initData() {

  double rho_bc, u_bc[3], p_bc, Z_bc, C_bc, T_bc, e_bc;
  Param * param = getParam(getName());
  rho_bc  = param->getDouble(1);
  u_bc[0] = param->getDouble(2);
  u_bc[1] = param->getDouble(3);
  u_bc[2] = param->getDouble(4);
  p_bc    = param->getDouble(5);
  Z_bc    = param->getDouble(6);
  C_bc    = param->getDouble(7);

  // get baseline R, T, e, gamma, a_gam at the bc mixture
  double R_bc, T0, e0, gamma0, a_gam;
  solver->chemtable->lookup(&R_bc, "R", &T0, "T", &e0, "e", &gamma0, "gamma",
                            &a_gam, "a_gamma",&Z_bc,&C_bc,1);
  T_bc = p_bc/R_bc/rho_bc;
  e_bc = calcInternalEnergyFromTS(T_bc,T0,e0,R_bc,gamma0,a_gam);

  assert( bf_state == NULL);
  bf_state = new PremixedState[zone_ptr->nbf];
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    for (int i =0; i < 3; ++i)
      bf_state[ibf].u[i] = u_bc[i];
    bf_state[ibf].sp_vol = 1.0/rho_bc;
    bf_state[ibf].p      = p_bc;
    bf_state[ibf].h      = e_bc + p_bc/rho_bc;
    bf_state[ibf].Z      = Z_bc;
    bf_state[ibf].C      = C_bc;
  }

  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<PremixedSolver,PremixedRhs>::parseScalarBc(param);
}

void CbcTotalPtsP::initData() {
  Param * param = getParam(getName());
  assert(param != NULL);
  p_tot         = param->getDouble(1);
  T_tot         = param->getDouble(2);
  Z_bc          = param->getDouble(3);
  C_bc          = param->getDouble(4);
  t_relax       = param->getDouble(5);

  // Grab the gas constant for later
  solver->chemtable->lookup(&R_0, "R", &e_0, "e", &gamma_0, "gamma", &a_gam_0, "a_gamma", &T_0, "T", &Z_bc, &C_bc, 1);
  // These quantities are immutable, so set them in initData()
  h_0   = e_0 + R_0*T_0;                                                      // Tabulated enthalpy
  e_tot = e_0 + R_0/a_gam_0*log(1.0 + a_gam_0*(T_tot - T_0)/(gamma_0 - 1.0)); // Total internal energy
  h_tot = e_tot + R_0*T_tot;                                                  // Total enthalpy
  // Compute the pressure at which the injection/stagnation isentrope intersects the table isotherm
  p_star = exp(((gamma_0 - a_gam_0*T_tot)*log(T_tot/T_0) \
                + log((gamma_0 - 1.0)/(gamma_0 + a_gam_0*(T_tot - T_0) - 1.0))) \
               / (gamma_0 - a_gam_0*T_tot - 1.0))*solver->chemtable->pressure;

  // XXX This is really just debugging info...
  if (mpi_rank == 0) cout << "Some #s of interest for zone " << zone_ptr->getName()
                          << ": p_tot:\t"   << p_tot    << ", T_tot:\t"   << T_tot
                          << ", Z_bc:\t"    << Z_bc     << ", C_bc:\t"    << C_bc
                          << ", t_relax:\t" << t_relax  << ", R_0:\t"     << R_0
                          << ", gamma_0:\t" << gamma_0  << ", a_gam_0:\t" << a_gam_0
                          << ", T_0:\t"     << T_0      << ", e_0:\t"     << e_0
                          << ", h_0:\t"     << h_0      << ", p_star:\t"  << p_star
                          << ", e_tot:\t" << e_tot      << ", h_tot:\t"   << h_tot
                          << ", p_0:\t" << solver->chemtable->pressure    << endl;

  assert(bf_state == NULL);
  bf_state  = new PremixedState[zone_ptr->nbf];
  un_bc     = new double[zone_ptr->nbf];

  assert(mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<PremixedSolver, PremixedRhs>::parseScalarBc(param);
}

inline void CbcTotalPtsP::calcBfState(const int ibf, const double unit_n[3]) {
  // NB: 0 == tabulated value, tot = stagnation value, bc == static value
  // Use conservation of enthalpy to obtain current temperature
  const double h_sns  = h_tot - 0.5*un_bc[ibf]*un_bc[ibf];// static enthalpy
  const double T_bc   = T_0 + (gamma_0-1.0)/a_gam_0*(exp((h_sns - h_0)*a_gam_0/R_0) - 1.0); // static T
  // Now we need an isentropic relation to define thermodynamic state
  // g(T)/(g(T) - 1) * dT/T = dp/p
  //
  // integrates to, (per Wolfram)
  //
  // ln(T')*(a_gam_0*T_0 - gamma_0) + ln(a_gam_0*(T_0 - T_bc) - gamma_0 + 1)/(a_gam_0*T_0 - gam_0 + 1)|T_0->T_bc = ln(p_bc/p_0)
  //
  // ls = left side, lb/ub = lower/upper bound
  const double lhs    = ((gamma_0 - a_gam_0*T_bc)*log(T_bc/T_0) \
                         + log((gamma_0 - 1.0)/(gamma_0 + a_gam_0*(T_bc - T_0) - 1.0))) \
    / (gamma_0 - a_gam_0*T_bc - 1.0);
  const double p_bc   = exp(lhs)*solver->chemtable->pressure*p_tot/p_star;
  // And lastly, the ideal gas law will provide us with the density
  const double rho_bc = p_bc/(R_0*T_bc);

  for (int i =0; i < 3; ++i)
    bf_state[ibf].u[i] = un_bc[ibf]*unit_n[i];
  bf_state[ibf].sp_vol = 1.0/rho_bc;
  bf_state[ibf].p      = p_bc;
  bf_state[ibf].h      = h_sns;
  bf_state[ibf].Z      = Z_bc;
  bf_state[ibf].C      = C_bc;
}

void CbcTotalPtsP::initialHook() {
  if (!zone_ptr->checkDataFlag("un_bc")) {
    if (mpi_rank == 0)
      cout << " > zone : " << getName() << ": reinitializing un_bc for cbc total pt ..." << endl;

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv         = zone_ptr->cvobf[ibf];
      const double bf_area  = MAG(zone_ptr->n_bf[ibf]);
      double unit_n[3]      = {zone_ptr->n_bf[ibf][0]/bf_area,
                               zone_ptr->n_bf[ibf][1]/bf_area,
                               zone_ptr->n_bf[ibf][2]/bf_area};
      un_bc[ibf]            = DOT_PRODUCT(unit_n, solver->u[icv]);
      calcBfState(ibf, unit_n);
    }
  }
  else {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const double bf_area = MAG(zone_ptr->n_bf[ibf]);
      double unit_n[3]      = {zone_ptr->n_bf[ibf][0]/bf_area,
                               zone_ptr->n_bf[ibf][1]/bf_area,
                               zone_ptr->n_bf[ibf][2]/bf_area};
      calcBfState(ibf, unit_n);
    }
  }
}

void CbcTotalPtsP::rk3Step(const double dt, const int rk_stage) {
  if (rk_stage == 1) { // updates the value at the beginning of the step...
    const double eps = 1.0 / (t_relax / dt + 1.0);
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv         = zone_ptr->cvobf[ibf];
      const double bf_area  = MAG(zone_ptr->n_bf[ibf]);
      double unit_n[3]      = {zone_ptr->n_bf[ibf][0]/bf_area,
                               zone_ptr->n_bf[ibf][1]/bf_area,
                               zone_ptr->n_bf[ibf][2]/bf_area};
      const double un0      = DOT_PRODUCT(solver->u[icv], unit_n);
      un_bc[ibf]            = (1.0 - eps)*un_bc[ibf] + eps*un0;

      // with this new un_bc, we need to recompute the bf_state
      // This is held frozen for the rest of the step...
      calcBfState(ibf, unit_n);
    }
  }
}


NscbcP::~NscbcP() {
  DELETE(p_bc);
  DELETE(rho_bc);
  DELETE(u_bc);
  DELETE(Z_bc);
  DELETE(C_bc);
  DELETE(rhs_bc);

  DELETE(R_bc);
  DELETE(T0);
  DELETE(e0);
  DELETE(gamma0);
  DELETE(agam);
}

void NscbcP::initData() {
  p_bc    = new double[zone_ptr->nbf];
  rho_bc  = new double[zone_ptr->nbf];
  u_bc    = new double[zone_ptr->nbf][3];
  Z_bc    = new double[zone_ptr->nbf];
  C_bc    = new double[zone_ptr->nbf];
  rhs_bc  = new double[zone_ptr->nbf][3][7];

  // also need some auxilliary variables to look up the state..
  R_bc    = new double[zone_ptr->nbf];
  T0      = new double[zone_ptr->nbf];
  e0      = new double[zone_ptr->nbf];
  gamma0  = new double[zone_ptr->nbf];
  agam    = new double[zone_ptr->nbf];

}

void NscbcP::initialHook() {

  if ( (!zone_ptr->checkDataFlag("p_bc"))   ||
       (!zone_ptr->checkDataFlag("rho_bc")) ||
       (!zone_ptr->checkDataFlag("u_bc"))   ||
       (!zone_ptr->checkDataFlag("Z_bc"))   ||
       (!zone_ptr->checkDataFlag("C_bc"))     ) {

    // if any of the variables are not instantiated, then
    // it is considered that all of the variables are not instantiated..
    if ( mpi_rank == 0 )
      cout << "zone : " << getName() << " : reinitializing nscbc vars .. " << endl;

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv    = zone_ptr->cvobf[ibf];
      rho_bc[ibf]      = solver->rho[icv];
      p_bc[ibf]        = solver->p[icv];
      for (int i =0; i < 3; ++i)
        u_bc[ibf][i]   = solver->u[icv][i];
      Z_bc[ibf]        = solver->Z[icv];
      C_bc[ibf]        = solver->C[icv];
    }
  }
  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    for (int i =0; i < 3; ++i)
      for (int j =0; j < 7; ++j)
        rhs_bc[ibf][i][j] = 0.0;
  }

  updatePrimitiveData();
}

int NscbcP::storeData(double * buf) const {

  if (buf) {
    int count = 0;
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      buf[count++] = rho_bc[ibf];
      buf[count++] = p_bc[ibf];
      FOR_I3 buf[count++] = u_bc[ibf][i];
      buf[count++] = Z_bc[ibf];
      buf[count++] = C_bc[ibf];
    }
    assert(count == zone_ptr->nbf*7);
  }

  return zone_ptr->nbf*7; // rho,p,u(3),Z,C

}

int NscbcP::restoreData(double * buf) {

  assert(buf);
  int count = 0;
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    rho_bc[ibf] = buf[count++];
    p_bc[ibf] = buf[count++];
    FOR_I3 u_bc[ibf][i] = buf[count++];
    Z_bc[ibf] = buf[count++];
    C_bc[ibf] = buf[count++];
  }
  assert(count == zone_ptr->nbf*7);

  updatePrimitiveData();

  return zone_ptr->nbf*7;

}

void NscbcP::query(const string& param_str) const {

  // report integrated values for the p_bc, rho_bc, u_bc, Z_bc, C_bc, T_bc
  double my_buf[7] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double buf[7];

  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    const double area    = MAG(zone_ptr->n_bf[ibf]);
    const double undA_bc = DOT_PRODUCT(zone_ptr->n_bf[ibf],u_bc[ibf]);
    my_buf[0]           += area;                                  // int(dA)            Rey
    my_buf[1]           += area*rho_bc[ibf];;                     // int(rho*dA)        Rey = Favre
    my_buf[2]           +=      rho_bc[ibf]*undA_bc;              // int(rho*dot(u,dA)) Favre
    my_buf[3]           += area            *p_bc[ibf];            // int(p*dA)          Rey
    my_buf[4]           += area*rho_bc[ibf]*Z_bc[ibf];            // int(rho*Z*dA)      Favre
    my_buf[5]           += area*rho_bc[ibf]*C_bc[ibf];            // int(rho*C*dA)      Favre
    my_buf[6]           += area            *p_bc[ibf]/R_bc[ibf];  // int(rho*T*dA)      Favre XXX: p/R = rho*T
  }

  MPI_Reduce(my_buf,buf,7,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if (mpi_rank == 0)
    cout  << " > Zone "   << zone_ptr->getName()
          << ": time = "  << solver->time
          << ", mdot = "  << buf[2]
          << ", rho = "   << buf[1]/buf[0]
          << ", un = "    << buf[2]/buf[1]
          << ", p = "     << buf[3]/buf[0]
          << ", T = "     << buf[6]/buf[1]
          << ", Z = "     << buf[4]/buf[1]
          << ", C = "     << buf[5]/buf[1]
          << endl;
}

void NscbcP::rk3Step(const double dt_, const int rk_stage) {
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const double dt = dt_*solver->dd[zone_ptr->cvobf[ibf]];

    for (int i =0; i < rk_stage; ++i) {
      rho_bc[ibf] += erk_wgt3[rk_stage-1][i]*dt*rhs_bc[ibf][i][0];
      p_bc[ibf]   += erk_wgt3[rk_stage-1][i]*dt*rhs_bc[ibf][i][1];
      for (int j =0; j < 3; ++j)
        u_bc[ibf][j] += erk_wgt3[rk_stage-1][i]*dt*rhs_bc[ibf][i][2+j];
      Z_bc[ibf]   += erk_wgt3[rk_stage-1][i]*dt*rhs_bc[ibf][i][5];
      C_bc[ibf]   += erk_wgt3[rk_stage-1][i]*dt*rhs_bc[ibf][i][6];
    }
  }
}

void NscbcP::addBoundaryFlux(PremixedRhs* rhs) const {

  assert( solver != NULL);
  assert( mf != NULL);

  // assumes that the R_bc,gamma0,agam have been updated properly

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const double T_bc     = p_bc[ibf]/rho_bc[ibf]/R_bc[ibf];
    const double e_bc     = calcInternalEnergyFromTS(T_bc,T0[ibf],e0[ibf],R_bc[ibf],gamma0[ibf],agam[ibf]);
    const double gamma_   = gamma0[ibf] + agam[ibf]*(T_bc - T0[ibf]);
    const int icv         = zone_ptr->cvobf[ibf];

    // pack the state to compute the riemann flux; if this packing becomes a
    // critical barrier to the performance, we can consider an alternative
    // riemann flux signature

    PremixedState bf_state;
    for (int i =0; i < 3; ++i)
      bf_state.u[i] = u_bc[ibf][i];
    bf_state.sp_vol = 1.0/rho_bc[ibf];
    bf_state.p      = p_bc[ibf];
    bf_state.h      = e_bc + p_bc[ibf]/rho_bc[ibf];
    bf_state.Z      = Z_bc[ibf];
    bf_state.C      = C_bc[ibf];

    PremixedRhs flux;
    calcRiemannFlux(flux,zone_ptr->n_bf[ibf],solver->cv_light[icv],bf_state,solver->gamma[icv],gamma_);

    mf[ibf]       = flux.rho;

    rhs[icv].rho -= flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv].rhou[i] -= flux.rhou[i];
    rhs[icv].rhoE -= flux.rhoE;
    rhs[icv].rhoZ -= flux.rhoZ;
    rhs[icv].rhoC -= flux.rhoC;
  }
}

void NscbcP::addScalarBoundaryFlux(double* rhs_sc) const {
  BcTemplate<PremixedSolver,PremixedRhs>::_addScalarBoundaryFlux(rhs_sc);
}

void NscbcMtsP::initData() {

  NscbcP::initData();

  Param* param  = getParam(getName());
  assert(param);

  // look for time-dependent bcs by searching for the string "time" in the
  // bcs. The user may introduce time through another parameter so this should
  // be made more general in the future...

  if (mpi_rank == 0) cout << getName() << " NSCBC_MTS:";

  // mdot_bc...
  if (param->getString(1).find("time") != std::string::npos) {
    mdot_bc_str = param->getString(1);
    eval_mdot_bc = true;
    if (mpi_rank == 0) cout << " " << mdot_bc_str;
  }
  else {
    assert(!eval_mdot_bc);
    const double mdot_bc = param->getDouble(1);
    // the global area has already been computed by the zone_ptr
    // if the mass is not changing as a function of time for this
    // bc, we can set rhoun_bc...
    rhoun_bc = -mdot_bc/zone_ptr->area_global;
    if (mpi_rank == 0) cout << " " << mdot_bc;
  }

  // T_bc...
  if (param->getString(2).find("time") != std::string::npos) {
    T_bc_str = param->getString(2);
    eval_T_bc = true;
    if (mpi_rank == 0) cout << " " << T_bc_str;
  }
  else {
    assert(!eval_T_bc);
    T_bc = param->getDouble(2);
    if (mpi_rank == 0) cout << " " << T_bc;
  }

  // Z_bc...
  if (param->getString(3).find("time") != std::string::npos) {
    Z_bc_str = param->getString(3);
    eval_Z_bc = true;
    if (mpi_rank == 0) cout << " " << Z_bc_str;
  }
  else {
    assert(!eval_Z_bc);
    Z_bc_ = param->getDouble(3);
    if (mpi_rank == 0) cout << " " << Z_bc_;
  }

  // C_bc...
  if (param->getString(4).find("time") != std::string::npos) {
    C_bc_str = param->getString(4);
    eval_C_bc = true;
    if (mpi_rank == 0) cout << " " << C_bc_str;
  }
  else {
    assert(!eval_C_bc);
    C_bc_ = param->getDouble(4);
    if (mpi_rank == 0) cout << " " << C_bc_;
  }

  if (mpi_rank == 0) cout << endl;

  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<PremixedSolver,PremixedRhs>::parseScalarBc(param);
}

void NscbcMtsP::evalZCStrAndUpdate() {

  if (eval_Z_bc) {
    CtiRegister::CtiData * var = CtiRegister::getCtiData(Z_bc_str);
    assert((var != NULL) && (var->getType() == D_DATA));
    Z_bc_ = var->d();
  }

  if (eval_C_bc) {
    CtiRegister::CtiData * var = CtiRegister::getCtiData(C_bc_str);
    assert((var != NULL) && (var->getType() == D_DATA));
    C_bc_ = var->d();
  }

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    Z_bc[ibf] = Z_bc_;
    C_bc[ibf] = C_bc_;
  }

  if (eval_Z_bc || eval_C_bc)
    solver->chemtable->lookup(R_bc, "R", T0, "T", e0, "e", gamma0, "gamma",
                              agam, "a_gamma", Z_bc, C_bc, zone_ptr->nbf);

  if (eval_T_bc) {
    CtiRegister::CtiData * var = CtiRegister::getCtiData(T_bc_str);
    assert((var != NULL) && (var->getType() == D_DATA));
    T_bc = var->d();
  }

  if (eval_mdot_bc) {
    CtiRegister::CtiData * var = CtiRegister::getCtiData(mdot_bc_str);
    assert((var != NULL) && (var->getType() == D_DATA));
    rhoun_bc = -var->d()/zone_ptr->area_global;
  }
}

void NscbcMtsP::initialHook() {

  NscbcP::initialHook();

  // if the Z_bc, C_bc are not changing as a function
  // of time for this boundary condition, we can lookup
  // the gas properties during the initialization

  if ((!eval_Z_bc)&&(!eval_C_bc)) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      Z_bc[ibf] = Z_bc_;
      C_bc[ibf] = C_bc_;
    }
    solver->chemtable->lookup(R_bc,"R", T0, "T", e0, "e", gamma0, "gamma",
			      agam, "a_gamma",Z_bc,C_bc,zone_ptr->nbf);
  }

  // for subsonic MTS bc, only rho is evolved and needs to be initialized...
  if ( (!zone_ptr->checkDataFlag("rho_bc")) ) {
    IF_RANK0 cout << "zone : " << getName() << " : reinitializing rho ... " << endl;

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv = zone_ptr->cvobf[ibf];
      rho_bc[ibf]   = solver->rho[icv];
    }
  }

  evalZCStrAndUpdate();
}

void NscbcMtsP::calcRhs(const double time, const int rk_stage) {

  evalZCStrAndUpdate();

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    // ensure that the the state is consistent with the mdot and T...
    const double bf_area = MAG(zone_ptr->n_bf[ibf]);
    double unit_n[3]     = {zone_ptr->n_bf[ibf][0]/bf_area,
                            zone_ptr->n_bf[ibf][1]/bf_area,
                            zone_ptr->n_bf[ibf][2]/bf_area};

    const double un_bc   = rhoun_bc/rho_bc[ibf];

    for (int i =0; i <3; ++i)
      u_bc[ibf][i] = un_bc*unit_n[i];

    p_bc[ibf]      = rho_bc[ibf]*R_bc[ibf]*T_bc;

    const int icv       = zone_ptr->cvobf[ibf];
    const double delta  = bf_area / zone_ptr->area_over_delta_bf[ibf];
    const double gamma_ = gamma0[ibf] + agam[ibf] * (T_bc - T0[ibf]);

    int ierr = calcNscbcIgSubsonicInletRhs(rhs_bc[ibf][rk_stage-1],unit_n,delta,
					   rho_bc[ibf],u_bc[ibf],p_bc[ibf],
					   solver->rho[icv],solver->u[icv],solver->p[icv],gamma_);
    if (ierr != 0) {
      cout << "calcNscbcIgSubsonicInletRhs failed for bc: " << zone_ptr->getName() << endl;
      assert(0);
    }

  }
}

void NscbcMtsSoftP::initData() {

  NscbcP::initData();

  double mdot_bc;
  Param* param  = getParam(getName());
  mdot_bc       = param->getDouble(1);
  T_bc          = param->getDouble(2);
  Z_bc_         = param->getDouble(3);
  C_bc_         = param->getDouble(4);

  if ( param->getString(5) == "T_RELAX") {
    t_relax = param->getDouble(6);
  } else {
    // XXX need to set some default values so that the t_relax
    // below becomes optional..
    CERR( " > NSCBC_MTS_SOFT <mdot> <T> <Z> <C> T_RELAX <t_relax> " );
  }

  // the global area has already been computed by the zone_ptr

  rhoun_bc = -mdot_bc/zone_ptr->area_global;

  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<PremixedSolver,PremixedRhs>::parseScalarBc(param);
}

void NscbcMtsSoftP::initialHook() {

  NscbcP::initialHook();

  // since the Z_bc, C_bc are not changing as a function
  // of time for this boundary condition, we can lookup
  // the gas properties during the initialization

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    Z_bc[ibf] = Z_bc_;
    C_bc[ibf] = C_bc_;
  }

  solver->chemtable->lookup(R_bc,"R", T0, "T", e0, "e", gamma0, "gamma",
                            agam, "a_gamma",Z_bc,C_bc,zone_ptr->nbf);

  // for soft nscbc , we need to transport rho, p (mdot is enforced),
  if ( (!zone_ptr->checkDataFlag("rho_bc")) ||
       (!zone_ptr->checkDataFlag("p_bc"))      ) {

    if ( mpi_rank == 0 )
      cout << "zone : " << getName() << " : reinitializing rho, p ... " << endl;

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv = zone_ptr->cvobf[ibf];
      rho_bc[ibf]   = solver->rho[icv];
      p_bc[ibf]     = solver->p[icv];
      for (int i =0; i < 3; ++i)
        u_bc[ibf][i] = solver->u[icv][i];
    }
  }
}

void NscbcMtsSoftP::calcRhs(const double time, const int rk_stage) {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    // ensure that the the state is consistent with the mdot
    const double bf_area = MAG(zone_ptr->n_bf[ibf]);
    double unit_n[3]     = {zone_ptr->n_bf[ibf][0]/bf_area,
                            zone_ptr->n_bf[ibf][1]/bf_area,
                            zone_ptr->n_bf[ibf][2]/bf_area};
    const double un_bc   = rhoun_bc/rho_bc[ibf];
    for (int i =0; i <3; ++i)
      u_bc[ibf][i] = un_bc*unit_n[i];

    const double T_bc_  = p_bc[ibf]/(rho_bc[ibf]*R_bc[ibf]);
    const int icv       = zone_ptr->cvobf[ibf];
    const double delta  = bf_area / zone_ptr->area_over_delta_bf[ibf];
    const double gamma_ = gamma0[ibf] + agam[ibf] * (T_bc_ - T0[ibf]);
    calcNscbcIgSubsonicInletSoftRhs(rhs_bc[ibf][rk_stage-1],unit_n,delta,
                                    rho_bc[ibf],u_bc[ibf],p_bc[ibf],
                                    solver->rho[icv],solver->u[icv],solver->p[icv],
                                    gamma_, T_bc, t_relax, R_bc[ibf]);
  }
}

void NscbcOutletPP::initData() {

  NscbcP::initData();

  Param * param = getParam(getName());
  if ( param->getString(0) == "NSCBC_OUTLET_P") {

    p_ref         = param->getDouble(1);
    sigma         = param->getDouble(2);
    L_ref         = param->getDouble(3);

  } else if ( param->getString(0) == "NSCBC_OUTLET_PRESSURE") {

    int iarg = 1;
    sigma    = 0.2;
    p_ref    = -1.0; // set to a negative value for error checking below..
    L_ref    = 1.0;  // XXX estimate a better value from the geometry here..

    while ( iarg < param->size()) {

      string token = param->getString(iarg++);
      if ( token == "SIGMA") {
        sigma = param->getDouble(iarg++);
      } else if ( token == "L") {
        L_ref = param->getDouble(iarg++);
      } else if ( token == "P_REF") {
        p_ref = param->getDouble(iarg++);
      } else if ( token == "CONSTANT_PRESSURE") {
        sigma = -1.0; // negative sigma value used to denote a constant pressure condition...
      } else {
        CWARN( " > skipping unrecognized token in NSCBC_OUTLET_PRESSURE: " << token);
      }
    }

    if ( p_ref < 0.0) {
      CERR( " > did not find P_REF specified for bc: " << zone_ptr->getName() );
    }

  } else {
    // if you get here there was a boundary condition parsing problem
    assert(0);
  }

  assert( mf == NULL); mf = new double[zone_ptr->nbf];

  // if reverse flow is encountered at the outlet, the passive scalars
  // will be reintroduced with 0 value as a default.

  assert( scbc_vals == NULL); scbc_vals = new double[solver->nsc];
  for (int isc = 0; isc < solver->nsc; ++isc)
    scbc_vals[isc] = 0.0;

}

void NscbcOutletPP::updatePrimitiveData() {
  // since Z_bc, C_bc are functions of time, we need to update the R_bc, etc..
  solver->chemtable->lookup(R_bc,"R", T0, "T", e0, "e", gamma0, "gamma",
                            agam, "a_gamma",Z_bc,C_bc,zone_ptr->nbf);
}

void NscbcOutletPP::calcRhs(const double time, const int rk_stage) {

  updatePrimitiveData();

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    // set the state to be consistent with the boundary vars..

    if ( sigma < 0.0) {
      p_bc[ibf] = p_ref; // a negative sigma indicates a constant pressure condition..
    }

    const double T_bc   = p_bc[ibf]/rho_bc[ibf]/R_bc[ibf];

    const int icv          = zone_ptr->cvobf[ibf];
    const double bf_area   = MAG(zone_ptr->n_bf[ibf]);
    const double unit_n[3] = {zone_ptr->n_bf[ibf][0]/bf_area,
                              zone_ptr->n_bf[ibf][1]/bf_area,
                              zone_ptr->n_bf[ibf][2]/bf_area};
    const double delta     = bf_area/zone_ptr->area_over_delta_bf[ibf];

    const double gamma_    = gamma0[ibf] + agam[ibf]*(T_bc - T0[ibf]);
    const double un_bc     = calcNscbcIgOutletRhs(rhs_bc[ibf][rk_stage-1],unit_n,delta,
                                                  rho_bc[ibf],u_bc[ibf],p_bc[ibf],
                                                  solver->rho[icv],solver->u[icv],solver->p[icv],
                                                  gamma_,L_ref,sigma,p_ref);

    const double dZdn      = (Z_bc[ibf] - solver->Z[icv])/delta;
    const double dCdn      = (C_bc[ibf] - solver->C[icv])/delta;

    rhs_bc[ibf][rk_stage-1][5] = -un_bc*dZdn;
    rhs_bc[ibf][rk_stage-1][6] = -un_bc*dCdn;
  }
}

void NscbcOutletPP::addBoundaryFlux(PremixedRhs* rhs) const {

  assert( solver != NULL);
  assert( mf != NULL);

  // assumes that the R_bc,gamma0,agam have been updated properly

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const double T_bc     = p_bc[ibf]/rho_bc[ibf]/R_bc[ibf];
    const double e_bc     = calcInternalEnergyFromTS(T_bc,T0[ibf],e0[ibf],R_bc[ibf],gamma0[ibf],agam[ibf]);
    const double gamma_   = gamma0[ibf] + agam[ibf]*(T_bc - T0[ibf]);
    const int icv         = zone_ptr->cvobf[ibf];

    // pack the state to compute the riemann flux; if this packing becomes a
    // critical barrier to the performance, we can consider an alternative
    // riemann flux signature

    PremixedState bf_state;
    for (int i =0; i < 3; ++i)
      bf_state.u[i] = u_bc[ibf][i];
    bf_state.sp_vol = 1.0/rho_bc[ibf];
    bf_state.p      = p_bc[ibf];
    bf_state.h      = e_bc + p_bc[ibf]/rho_bc[ibf];
    bf_state.Z      = Z_bc[ibf];
    bf_state.C      = C_bc[ibf];

    PremixedRhs flux;
    calcRiemannFlux(flux,zone_ptr->n_bf[ibf],solver->cv_light[icv],bf_state,solver->gamma[icv],gamma_);

    if ( flux.rho >= 0.0)  {

      mf[ibf]       = flux.rho;

      rhs[icv].rho -= flux.rho;
      for (int i =0; i < 3; ++i)
        rhs[icv].rhou[i] -= flux.rhou[i];
      rhs[icv].rhoE -= flux.rhoE;
      rhs[icv].rhoZ -= flux.rhoZ;
      rhs[icv].rhoC -= flux.rhoC;

    } else {

      // do not admit any backflow .. treat this face locally as a slip wall

      mf[ibf]        = 0.0;

      for (int i =0; i < 3; ++i)
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];

    }


  }
}

void NscbcOutletMdotP::initData() {

  NscbcP::initData();
  Param * p = getParam(getName());
  rhoun_bc = p->getDouble(1)/zone_ptr->area_global; // positive value denotes outflow here...

  assert( mf == NULL); mf = new double[zone_ptr->nbf];

  // if reverse flow is encountered at the outlet, the passive scalars
  // will be reintroduced with 0 value as a default.

  assert( scbc_vals == NULL); scbc_vals = new double[solver->nsc];
  for (int isc = 0; isc < solver->nsc; ++isc)
    scbc_vals[isc] = 0.0;
}

void NscbcOutletMdotP::updatePrimitiveData() {
  // since Z_bc, C_bc are functions of time, we need to update the R_bc, etc..
  solver->chemtable->lookup(R_bc,"R", T0, "T", e0, "e", gamma0, "gamma",
                            agam, "a_gamma",Z_bc,C_bc,zone_ptr->nbf);
}

void NscbcOutletMdotP::calcRhs(const double time, const int rk_stage) {

  updatePrimitiveData();

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv          = zone_ptr->cvobf[ibf];
    const double bf_area   = MAG(zone_ptr->n_bf[ibf]);
    const double unit_n[3] = {zone_ptr->n_bf[ibf][0]/bf_area,
                              zone_ptr->n_bf[ibf][1]/bf_area,
                              zone_ptr->n_bf[ibf][2]/bf_area};
    const double delta     = bf_area/zone_ptr->area_over_delta_bf[ibf];

    // the characteristic relations should be consistent with this behavior, but
    // set the velocity to be consistent with the prescribed mdot and overwrite the
    // u field..
    const double un_bc = rhoun_bc/rho_bc[ibf];
    for (int i =0; i < 3; ++i) {
      u_bc[ibf][i] = un_bc*unit_n[i];
    }

    const double T_bc      = p_bc[ibf]/rho_bc[ibf]/R_bc[ibf];
    const double gamma_    = gamma0[ibf] + agam[ibf]*(T_bc - T0[ibf]);
    calcNscbcIgOutletRhsMdot(rhs_bc[ibf][rk_stage-1],unit_n,delta,
                             rho_bc[ibf],u_bc[ibf],p_bc[ibf],solver->rho[icv],
                             solver->u[icv],solver->p[icv],gamma_);

    // update the rhs for the scalars as well..
    const double dZdn      = (Z_bc[ibf] - solver->Z[icv])/delta;
    const double dCdn      = (C_bc[ibf] - solver->C[icv])/delta;
    rhs_bc[ibf][rk_stage-1][5] = -un_bc*dZdn;
    rhs_bc[ibf][rk_stage-1][6] = -un_bc*dCdn;
  }
}

//
// sponge boundary condition implementation below.
// XXX the following implementation assumes that there is only one active sponge.
// otherwise the registration that the bc forces in the interior will run into
// name conflicts.  the bc registration is postponed until the initData
//

void SpongeP::initData() {

  Param* param = getParam(getName());

  // default params..
  eps_p_sponge = 0.01;
  sponge_pressure = solver->chemtable->pressure;

  parseSpongeParam(param,sponge_type,sponge_strength,eps_p_sponge,
                   t_relax,sponge_length,sponge_pressure);

  rho_sponge  = NULL; solver->registerCvData(rho_sponge,"rho_sponge",READWRITE_DATA);
  rhou_sponge = NULL; solver->registerCvData(rhou_sponge,"rhou_sponge",READWRITE_DATA);
  rhoE_sponge = NULL; solver->registerCvData(rhoE_sponge,"rhoE_sponge",READWRITE_DATA);
  rhoZ_sponge = NULL; solver->registerCvData(rhoZ_sponge,"rhoZ_sponge",READWRITE_DATA);
  rhoC_sponge = NULL; solver->registerCvData(rhoC_sponge,"rhoC_sponge",READWRITE_DATA);

  rho_sponge  = new double[solver->ncv];
  rhou_sponge = new double[solver->ncv][3];
  rhoE_sponge = new double[solver->ncv];
  rhoZ_sponge = new double[solver->ncv];
  rhoC_sponge = new double[solver->ncv];


  double my_buf[5] = { -1.0E-20, -1.0E-20, -1.0E-20, -1.0E-20, -1.0E-20 };

  for (int icv = 0; icv < solver->ncv; ++icv) {
    my_buf[0] = max(my_buf[0],solver->x_cv[icv][0]);
    my_buf[1] = max(my_buf[1],solver->x_cv[icv][1]);
    my_buf[2] = max(my_buf[2],solver->x_cv[icv][2]);
    my_buf[3] = max(my_buf[3],sqrt(solver->x_cv[icv][1]*solver->x_cv[icv][1] + solver->x_cv[icv][2]*solver->x_cv[icv][2]));
    my_buf[4] = max(my_buf[4],sqrt(solver->x_cv[icv][0]*solver->x_cv[icv][0] + solver->x_cv[icv][1]*solver->x_cv[icv][1]));
  }
  MPI_Allreduce(my_buf,sponge_data,5,MPI_DOUBLE,MPI_MAX,mpi_comm);

  assert( mf == NULL); mf = new double[zone_ptr->nbf];

  // if reverse flow is encountered at the outlet, the passive scalars
  // will be reintroduced with 0 value as a default.

  assert( scbc_vals == NULL); scbc_vals = new double[solver->nsc];
  for (int isc = 0; isc < solver->nsc; ++isc)
    scbc_vals[isc] = 0.0;

}

SpongeP::~SpongeP() {
  DELETE(rho_sponge);
  DELETE(rhou_sponge);
  DELETE(rhoE_sponge);
  DELETE(rhoZ_sponge);
  DELETE(rhoC_sponge);
}

int SpongeP::storeData(double * buf) const {

  // sponge stores cell values for T and u

  if (buf) {
    int count = 0;
    for (int icv = 0; icv < solver->ncv; ++icv) {
      buf[count++] = rho_sponge[icv];
      FOR_I3 buf[count++] = rhou_sponge[icv][i];
      buf[count++] = rhoE_sponge[icv];
      buf[count++] = rhoZ_sponge[icv];
      buf[count++] = rhoC_sponge[icv];
    }
    assert(count = solver->ncv*7);
  }

  return 7*solver->ncv;

}

int SpongeP::restoreData(double * buf) {

  assert(buf);
  int count = 0;
  for (int icv = 0; icv < solver->ncv; ++icv) {
    rho_sponge[icv] = buf[count++];
    FOR_I3 rhou_sponge[icv][i] = buf[count++];
    rhoE_sponge[icv] = buf[count++];
    rhoZ_sponge[icv] = buf[count++];
    rhoC_sponge[icv] = buf[count++];
  }
  assert(count = solver->ncv*7);

  return 7*solver->ncv;

}

void SpongeP::addBoundaryFlux(PremixedRhs* rhs) const {

  assert( mf != NULL);

  double * T_p   = new double[zone_ptr->nbf];
  double * e0    = new double[zone_ptr->nbf];
  double * gamma = new double[zone_ptr->nbf];
  double * R     = new double[zone_ptr->nbf];
  double * a_gam = new double[zone_ptr->nbf];
  double * Z_bf  = new double[zone_ptr->nbf];
  double * C_bf  = new double[zone_ptr->nbf];

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];
    Z_bf[ibf]     = rhoZ_sponge[icv]/rho_sponge[icv];
    C_bf[ibf]     = rhoC_sponge[icv]/rho_sponge[icv];
  }

  solver->chemtable->lookup(R,"R",T_p,"T",e0,"e",gamma,"gamma",a_gam,"a_gamma",
                            Z_bf,C_bf,zone_ptr->nbf);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];
    const double sp_vol = 1.0/rho_sponge[icv];
    const double usq    = DOT_PRODUCT(rhou_sponge[icv],rhou_sponge[icv])*sp_vol;
    const double e_bc   = rhoE_sponge[icv]*sp_vol - 0.5*usq;
    const double deor   = (e_bc - e0[ibf])/R[ibf];
    const double adeor  = deor*a_gam[ibf];

    // series expansion about a_gam=0 for (exp(a_gam*(e-e0)/R)-1)/a_gam
    // original expr is ( exp( agamma*(e-e0)/R ) - 1.0 ) / agamma * (gamma[i]-1.0);
    const double dT     = (gamma[ibf] - 1.0)*deor*(1.0+adeor*0.5000000000000000*
                                                   (1.0+adeor*0.3333333333333333*
                                                    (1.0+adeor*0.2500000000000000*
                                                     (1.0+adeor*0.2000000000000000*
                                                      (1.0+adeor*0.1666666666666667)))));

    const double T_bc   = T_p[ibf] + dT;
    const double gamma_ = gamma[ibf] + a_gam[ibf]*dT;
    const double RT_    = T_bc*R[ibf];
    const double p_bc   = rho_sponge[icv]*RT_;

    PremixedState bf_state;
    bf_state.sp_vol     = sp_vol;
    bf_state.p          = p_bc;
    bf_state.h          = e_bc + p_bc*sp_vol;
    for (int i =0; i < 3; ++i)
      bf_state.u[i] = rhou_sponge[icv][i]*sp_vol;

    PremixedRhs flux;
    calcRiemannFlux(flux,zone_ptr->n_bf[ibf],solver->cv_light[icv],
                    bf_state,solver->gamma[icv],gamma_);

    mf[ibf]       = flux.rho;

    rhs[icv].rho -= flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv].rhou[i] -= flux.rhou[i];
    rhs[icv].rhoE -= flux.rhoE;
    rhs[icv].rhoZ -= flux.rhoZ;
    rhs[icv].rhoC -= flux.rhoC;
  }

  delete[] T_p;
  delete[] e0;
  delete[] gamma;
  delete[] R;
  delete[] a_gam;
  delete[] Z_bf;
  delete[] C_bf;

  switch( sponge_type) {
  case SPONGE_LX:
  case SPONGE_LY:
  case SPONGE_LZ:
    {

      int idir;
      if (sponge_type == SPONGE_LX) {
        idir = 0;
      } else if ( sponge_type == SPONGE_LY) {
        idir = 1;
      } else {
        assert( sponge_type == SPONGE_LZ);
        idir = 2;
      }

      const double Lmax = sponge_data[idir];

      for (int icv = 0; icv < solver->ncv; ++icv) {
        if ( solver->x_cv[icv][idir] > Lmax-sponge_length) {
          const double x_sp = (solver->x_cv[icv][idir] - (Lmax-sponge_length))/sponge_length;
          // sponge profile = a*x^2 + b*x^8..
          const double sponge_coeff = sponge_strength*(0.068*x_sp*x_sp + 0.845*pow(x_sp,8));
          rhs[icv].rho             -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv] - rho_sponge[icv]);
          for (int i =0; i < 3; ++i)
            rhs[icv].rhou[i]       -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv]*solver->u[icv][i] - rhou_sponge[icv][i]);
          rhs[icv].rhoE            -= solver->vol_cv[icv]*sponge_coeff*(solver->rhoE[icv] - rhoE_sponge[icv]);
          rhs[icv].rhoZ            -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv]*solver->Z[icv] - rhoZ_sponge[icv]);
          rhs[icv].rhoC            -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv]*solver->C[icv] - rhoC_sponge[icv]);
        }
      }
    }
    break;

  case SPONGE_LRX:
    {

      const double Lmax = sponge_data[3];
      for (int icv = 0; icv < solver->ncv; ++icv) {
        const double rx = sqrt( solver->x_cv[icv][1]*solver->x_cv[icv][1] +
                                solver->x_cv[icv][2]*solver->x_cv[icv][2]   );
        if ( rx > Lmax - sponge_length) {
          const double x_sp = (rx - (Lmax-sponge_length))/sponge_length;
          const double sponge_coeff = sponge_strength*(0.068*x_sp*x_sp + 0.845*pow(x_sp,8));
          rhs[icv].rho             -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv] - rho_sponge[icv]);
          for (int i =0; i < 3; ++i)
            rhs[icv].rhou[i]       -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv]*solver->u[icv][i] - rhou_sponge[icv][i]);
          rhs[icv].rhoE            -= solver->vol_cv[icv]*sponge_coeff*(solver->rhoE[icv] - rhoE_sponge[icv]);
          rhs[icv].rhoZ            -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv]*solver->Z[icv] - rhoZ_sponge[icv]);
          rhs[icv].rhoC            -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv]*solver->C[icv] - rhoC_sponge[icv]);
        }
      }
    }
    break;

  case SPONGE_LRZ:
    {

      const double Lmax = sponge_data[4];
      for (int icv = 0; icv < solver->ncv; ++icv) {
        const double rz = sqrt ( solver->x_cv[icv][0]*solver->x_cv[icv][0] +
                                 solver->x_cv[icv][1]*solver->x_cv[icv][1]  );
        if ( rz > Lmax - sponge_length) {
          const double x_sp = (rz-(Lmax-sponge_length))/sponge_length;
          const double sponge_coeff = sponge_strength*(0.068*x_sp*x_sp + 0.845*pow(x_sp,8));
          rhs[icv].rho             -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv] - rho_sponge[icv]);
          for (int i =0; i < 3; ++i)
            rhs[icv].rhou[i]       -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv]*solver->u[icv][i] - rhou_sponge[icv][i]);
          rhs[icv].rhoE            -= solver->vol_cv[icv]*sponge_coeff*(solver->rhoE[icv] - rhoE_sponge[icv]);
          rhs[icv].rhoZ            -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv]*solver->Z[icv] - rhoZ_sponge[icv]);
          rhs[icv].rhoC            -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv]*solver->C[icv] - rhoC_sponge[icv]);
        }
      }
    }
    break;

  default:
    assert(0);  // these errors have been checked earlier..
  }
}

void SpongeP::rk3Step(const double dt, const int rk_stage) {

  if ( rk_stage == 1) {

    // relax the state at the beginning of the step (solver
    // data is consistent at time level n.. )

    const double eps = 1.0/ (t_relax/dt + 1.0);

    for (int icv = 0; icv < solver->ncv; ++icv) {
      rho_sponge[icv]   += eps*(solver->rho[icv] - rho_sponge[icv]);
      for (int i =0; i < 3; ++i)
        rhou_sponge[icv][i] += eps*(solver->rho[icv]*solver->u[icv][i] - rhou_sponge[icv][i]);
      rhoE_sponge[icv]      += eps*(solver->rhoE[icv] - rhoE_sponge[icv]);
      rhoZ_sponge[icv]      += eps*(solver->rho[icv]*solver->Z[icv] - rhoZ_sponge[icv]);
      rhoC_sponge[icv]      += eps*(solver->rho[icv]*solver->C[icv] - rhoC_sponge[icv]);
    }
  }
}


void SpongeP::initialHook() {

  if ( !(solver->checkDataFlag("rho_sponge")) ||
       !(solver->checkDataFlag("rhou_sponge")) ||
       !(solver->checkDataFlag("rhoE_sponge")) ||
       !(solver->checkDataFlag("rhoZ_sponge")) ||
       !(solver->checkDataFlag("rhoC_sponge"))    ) {

    for (int icv = 0; icv < solver->ncv; ++icv) {
      rho_sponge[icv] = solver->rho[icv];
      for (int i =0 ; i < 3; ++i)
        rhou_sponge[icv][i] = solver->rho[icv]*solver->u[icv][i];
      rhoE_sponge[icv] = solver->rhoE[icv];
      rhoZ_sponge[icv] = solver->rho[icv]*solver->Z[icv];
      rhoC_sponge[icv] = solver->rho[icv]*solver->C[icv];
    }
  }
}


void SpongeP::query(const string& param_str) const {

  // report integrated values for the p_bc, rho_bc, u_bc, Z_bc, C_bc
  double my_buf[5] = {0.0,0.0,0.0,0.0,0.0};
  double buf[5];

  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv           = zone_ptr->cvobf[ibf];
    const double area       = MAG(zone_ptr->n_bf[ibf]);
    const double rhoundA_bc = DOT_PRODUCT(zone_ptr->n_bf[ibf],rhou_sponge[icv]);
    const double Z_bc       = rhoZ_sponge[icv]/rho_sponge[icv];
    const double C_bc       = rhoC_sponge[icv]/rho_sponge[icv];

    my_buf[0]        += area;
    my_buf[1]        += area*rho_sponge[icv];
    my_buf[2]        += rhoundA_bc;
    my_buf[3]        += area*Z_bc;
    my_buf[4]        += area*C_bc;
  }

  MPI_Reduce(my_buf,buf,5,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  if ( mpi_rank == 0 ) {
    cout << " > " << zone_ptr->getName() << ": time, rho_bc, mdot, Z_bc, C_bc: "
         << solver->time << "   " << buf[1]/buf[0] << "    "
         << buf[2] << "    " << buf[3]/buf[0] << "    " << buf[4]/buf[0] << endl;
  }
}

void SpongeP::addScalarBoundaryFlux(double* rhs_sc) const {
  BcTemplate<PremixedSolver,PremixedRhs>::_addScalarBoundaryFlux(rhs_sc);
}

void InflowTurbulenceP::initData() {
  um_it         = new double [zone_ptr->nbf][3];
  Rd_it         = new double [zone_ptr->nbf][3];
  Rod_it        = new double [zone_ptr->nbf][3];
  Tm_it         = new double [zone_ptr->nbf];
  pm_it         = new double [zone_ptr->nbf];
  Zm_it         = new double [zone_ptr->nbf];
  Cm_it         = new double [zone_ptr->nbf];
  length_scale  = new double [zone_ptr->nbf];
  u_bc          = new double [zone_ptr->nbf][3];
  T_bc          = new double [zone_ptr->nbf];
  p_bc          = new double [zone_ptr->nbf];

  avg           = new double [zone_ptr->nbf][3];
  rms           = new double [zone_ptr->nbf][3];

  //  This store the state at the boundary
  assert( bf_state == NULL);
  bf_state = new PremixedState[zone_ptr->nbf];

  IturbWgt = 0.0;
  max_length = 1.0E+200;
  in_plane_length = 1.0E+200;
  Un = 1.0E+200;
  inflow_zero = 0.001;
  inflow_maxiter = 30;

  Param* param = getParam(getName());
  assert(mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<PremixedSolver, PremixedRhs>::parseScalarBc(param);
}

void InflowTurbulenceP::initialHook() {
  // Set the statistics of the turbulence... Should this only happen at step 0? Or when the user asks?
  setStatistics();

  // We need the reduced connectivity of the boundary attahed cvs. These are computed here
  buildITBoundaryStructures();

  // Now we compute and store the lengthscale based on the distance to the wall and/or
  // whatever other input the user has specified. Same problem as set statistics
  setLengthscale();

  // Create the operators
  buildFilterOperators();

  // This should only be donw when we are at step 0... or if the user would like to reset the
  // inflow turbulence...
  initiateStats();
}

int InflowTurbulenceP::storeData(double * buf) const {

  // should only need to store avg, rms
  // other values are static or updated from solver state

  if (buf) {
    int count = 0;
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      FOR_I3 buf[count++] = avg[ibf][i];
      FOR_I3 buf[count++] = rms[ibf][i];
    }
    assert(count == zone_ptr->nbf*6);
  }

  return 6*zone_ptr->nbf;

}

int InflowTurbulenceP::restoreData(double * buf) {

  assert(buf);
  int count = 0;
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    FOR_I3 avg[ibf][i] = buf[count++];
    FOR_I3 rms[ibf][i] = buf[count++];
  }
  assert(count == zone_ptr->nbf*6);

  return 6*zone_ptr->nbf;

}

void InflowTurbulenceP::calcRhs(const double time, const int rk_stage) {
  if (rk_stage == 1)
    updateInflowTurbulence(u_bc);
}
