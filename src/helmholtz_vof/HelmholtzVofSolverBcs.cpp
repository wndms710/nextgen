#include "HelmholtzVofSolver.hpp"
#include "wm/AlgebraicWM.hpp"
#include "wm/AlgebraicRoughnessWM.hpp"

void SlipWallVBc::force_bf(double (*f_bf)[9]) const {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    FOR_I3 {
      f_bf[ibf][3*0+i] = 0.0;
      f_bf[ibf][3*1+i] = solver->p[icv0]*zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = 0.0;
    }
  }
}

void WallVBc::initData() {
 
  assert(solver != NULL);
  assert(u_bc == NULL);  u_bc = new double[zone_ptr->nbf][3];

  double uwall[3] = {0.0, 0.0, 0.0};
  double xrot[3]  = {0.0, 0.0, 0.0};
  double omega[3] = {0.0, 0.0, 0.0};
  double urot     = 0.0 ;
  
  bool isRotating = false;
  Param * param = getParam(getName()); 
  if ( (param->size() > 1) && (param->getString(1)=="U_WALL")){
    uwall[0]  = param->getDouble(2);
    uwall[1]  = param->getDouble(3);
    uwall[2]  = param->getDouble(4);
    COUT1(" > Found U_WALL " << u_bc[0] << ", " << u_bc[1] << ", " << u_bc[2]);
  }  
  else if ( (param->size() > 1) && (param->getString(1)=="ROTATING")){

    int ii = 2 ;
    while ( ii < param->size()) {
      FOR_I3 xrot[i] = param->getDouble(ii++);
      FOR_I3 omega[i] = param->getDouble(ii++);
      double mag = sqrt(DOT_PRODUCT(omega,omega));
      assert ( mag > 0.0) ;
      FOR_I3 omega[i] /= mag ;
      urot = param->getDouble(ii++)*M_PI*2.0 ;
      
      isRotating = true;
    }
  }

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    
    if (isRotating){
      double r[3];
      FOR_I3 r[i] = zone_ptr->x_bf[ibf][i] - xrot[i];
      double omega_cross_r[3] = CROSS_PRODUCT(omega,r);
      FOR_I3 u_bc[ibf][i] = urot * omega_cross_r[i];
    }
    else{
      FOR_I3 u_bc[ibf][i] = uwall[i];
    }
  }
}

void WallVBc::addMomentumFlux(double * A, double (*rhs)[3]) const { 

  
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    const int coc00 = solver->cvocv_i[icv0];
    const double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])*zone_ptr->area_over_delta_bf[ibf];
    A[coc00] += mu_coeff;
    FOR_I3 rhs[icv0][i] += mu_coeff*u_bc[ibf][i];
  }

}
 
void WallVBc::force_bf(double (*f_bf)[9]) const {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

    const int icv0 = zone_ptr->cvobf[ibf];
    double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])*zone_ptr->area_over_delta_bf[ibf];

    FOR_I3 { 
      f_bf[ibf][3*0+i] = 0.0;
      f_bf[ibf][3*1+i] = solver->p[icv0]*zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = mu_coeff*(solver->u[icv0][i] - u_bc[ibf][i]);
    }
  }

}


void WallVBc::setBc() {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
    const int icv0 = zone_ptr->cvobf[ibf];
    double mag_n = MAG(zone_ptr->n_bf[ibf]);
    double unit_n[3];
    FOR_I3 unit_n[i] = zone_ptr->n_bf[ibf][i]/mag_n;
    FOR_I3 solver->n[icv0][i] = -unit_n[i];
  }

}

void WallVBc::query(const string& param_str) {

  // report min/max slip length and velocity
  double my_buf[3] = {0.0,0.0,0.0}; // area, int tau_wall dA, y_plus
  double buf[3];

  double (*f_bf)[9] = new double [zone_ptr->nbf][9];
  this->force_bf(f_bf);

  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {

    double mag_n = MAG(zone_ptr->n_bf[ibf]);
    double unit_n[3];
    FOR_I3 unit_n[i] = zone_ptr->n_bf[ibf][i]/mag_n;
    double f_p[3];
    //get wall parallel component of force
    FOR_I3 f_p[i] = f_bf[ibf][3*0+i] + f_bf[ibf][3*2+i]; //ignore (wall normal) pressure term
    FOR_I3 f_p[i] -= DOT_PRODUCT(f_p,unit_n)*unit_n[i];

    const int icv = zone_ptr->cvobf[ibf];
    const double tau_wall = MAG(f_p) / zone_ptr->area_bf[ibf];

    const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
    const double u_tau = sqrt( tau_wall/solver->rho[icv]);
    const double nu    = solver->mu_lam[icv] / solver->rho[icv];
    const double y_plus = y1*u_tau/nu;

    my_buf[0] += zone_ptr->area_bf[ibf];
    my_buf[1] += zone_ptr->area_bf[ibf]*tau_wall;
    my_buf[2] += zone_ptr->area_bf[ibf]*y_plus;
  }
  delete[] f_bf;

  MPI_Reduce(my_buf,buf,3,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if (mpi_rank == 0) {
    *ss << " > " << std::left << std::setw(30) << zone_ptr->getName() << " :  time, int tau_wall dA, avg y_plus = " << std::right
        << std::setw(8) << solver->time << " "
        << std::setw(12) << buf[1] << " "
        << std::setw(12) << buf[2]/buf[0] << endl;
  }
  flush();

}

CtiRegister::CtiDataError WallVBc::funcEvalCtiData(CtiRegister::CtiData& v, 
    const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if ( zone_ptr->isLeadingStrMatched(name) ) { 

    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix.  we have to check for the actual data now..

    const string tau_str      = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str       = zone_ptr->getName() + ":" + "y_plus";
   
    // all of the following functions do not take any arguments -- check first

    if ( args.size() != 0) 
      return CTI_DATA_ARG_COUNT;
    
    if ( name == tau_str) { 

      double * tmp = zone_ptr->createBfD1Data(v); 

      // if stats have been requested on this quantity, this function will be called
      // but its possible that the solver->mu_lam, rho necessary have not been populated yet.
      // check the data flags for these values -- revisit if this impacts performance 
      // during the solution run-time.

      if ( b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          
          const int icv           = zone_ptr->cvobf[ibf];
          const double visc_coeff = (solver->mu_lam[icv] + solver->mu_sgs[icv])*
                                    zone_ptr->area_over_delta_bf[ibf];
          double du_mag = 0.0;
          for (int i =0; i < 3; ++i) 
            du_mag  += (solver->u[icv][i] - u_bc[ibf][i])*(solver->u[icv][i] - u_bc[ibf][i]);
        
          tmp[ibf]                = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
        }

      } 

      return CTI_DATA_OK;
    
    } else if ( name == yp_str) { 

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) { 

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

          const int icv           = zone_ptr->cvobf[ibf];
          const double y1         = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
          const double visc_coeff = (solver->mu_lam[icv] + solver->mu_sgs[icv])/y1;
 
          double du_mag = 0.0;
          for (int i =0; i < 3; ++i) 
            du_mag  += (solver->u[icv][i] - u_bc[ibf][i])*(solver->u[icv][i] - u_bc[ibf][i]);
        
          const double tau_w      = visc_coeff*sqrt(du_mag);
          tmp[ibf]                = y1* sqrt(tau_w*solver->rho[icv])/solver->mu_lam[icv]; 

        }

      }

      return CTI_DATA_OK;

    } else { 

      return CTI_DATA_NOT_FOUND;
    }
  }

  return CTI_DATA_NOT_FOUND;
}

void InletVBc::initData() {
 
  assert(u_bc   == NULL); u_bc   = new double[zone_ptr->nbf][3];
  assert(q_bf   == NULL); q_bf   = new double[zone_ptr->nbf];

  // just constant for now...
  double u_in[3],vof_in;
  double u_normal;
  Param * param = getParam(getName());
  if (param->size() <= 5) { // inlet_normal ... 
    COUT1("INLET_NORMAL:");
    rho0_bc  =  param->getDouble(1);
    rho1_bc  =  param->getDouble(2);
    vof_bc   =  param->getDouble(3);
    u_normal =  param->getDouble(4);
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
      double normal[3];
      const double mag = MAG(zone_ptr->n_bf[ibf]);
      assert(mag > 0.0);
      FOR_I3 normal[i] = zone_ptr->n_bf[ibf][i]/mag;
      //boundary normal vector is out-going...
      FOR_I3 u_bc[ibf][i] = -u_normal*normal[i];
      q_bf[ibf] = DOT_PRODUCT(zone_ptr->n_bf[ibf],u_bc[ibf]);
    }
  }
  else {
    rho0_bc = param->getDouble(1);
    rho1_bc = param->getDouble(2);
    vof_bc  = param->getDouble(3);
    u_in[0] = param->getDouble(4);
    u_in[1] = param->getDouble(5);
    u_in[2] = param->getDouble(6);
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
      FOR_I3 u_bc[ibf][i] = u_in[i];
      q_bf[ibf] = DOT_PRODUCT(zone_ptr->n_bf[ibf],u_bc[ibf]);
    }
  }

  dumpRange(q_bf,zone_ptr->nbf,"q_bf inlet");

}

void InletVBc::setBc() {

  // this routine would be non-empty if there was time variation in the inlet bc... 
 
}

void InletVBc::addFlux(double * rhs) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    rhs[icv0] += q_bf[ibf];
  }

}


void InletVBc::addMassFlux(double *rhs_rhoY0, double *rhs_rhoY1) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    rhs_rhoY0[icv0] -= rho0_bc*vof_bc*q_bf[ibf];
    rhs_rhoY1[icv0] -= rho1_bc*(1.0-vof_bc)*q_bf[ibf];
  }
}


void InletVBc::addMomentumFlux(double * A, double (*rhs)[3]) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    const int coc00 = solver->cvocv_i[icv0];
    const double rho_bc = rho1_bc*(1.0-vof_bc) + rho0_bc*vof_bc;
    const double mu_coeff = ((solver->mu1_ref*(1.0-vof_bc)+solver->mu0_ref*vof_bc) + solver->mu_sgs[icv0])*zone_ptr->area_over_delta_bf[ibf];
    A[coc00] += mu_coeff;
    FOR_I3 rhs[icv0][i] += (-rho_bc*q_bf[ibf] + mu_coeff)*u_bc[ibf][i];

    /*
    if (q_bf[ibf] >= 0.0) {
      // outflow -- use first-order upwind...
      //A[coc00] += rho_bc*q_bf[ibf] + mu_coeff;
      A[coc00] += mu_coeff;
      FOR_I3 rhs[icv0][i] += mu_coeff*u_bc[ibf][i];
    }
    else {
      // inflow...
      A[coc00] += mu_coeff;
      FOR_I3 rhs[icv0][i] += (-rho_bc*q_bf[ibf] + mu_coeff)*u_bc[ibf][i];

      //A[coc00] += 0.5*rho_bc*q_bf[ibf]+mu_coeff;
      //FOR_I3 rhs[icv0][i] += (-0.5*rho_bc*q_bf[ibf] + mu_coeff)*u_bc[ibf][i];
    }
    */
  }

}

void InletVBc::force_bf(double (*f_bf)[9]) const {


  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

    const int icv0 = zone_ptr->cvobf[ibf];
    const double rho_bc = rho1_bc*(1.0-vof_bc) + rho0_bc*vof_bc;
    const double mu_coeff = ((solver->mu1_ref*(1.0-vof_bc)+solver->mu0_ref*vof_bc) + solver->mu_sgs[icv0])*zone_ptr->area_over_delta_bf[ibf];
    FOR_I3 { 
      f_bf[ibf][3*0+i] = rho_bc*q_bf[ibf]*u_bc[ibf][i];
      f_bf[ibf][3*1+i] = solver->p[icv0]*zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = mu_coeff*(solver->u[icv0][i] - u_bc[ibf][i]);
    }
  }

}

void InletVBc::query(const string& param_str) {


  // report min/max
  double my_buf[2] = {HUGE_VAL,HUGE_VAL};
  double buf[2];

  double my_psum = 0.0;

  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    my_buf[0] = min(my_buf[0],solver->p[icv0]);
    my_buf[1] = min(my_buf[1],-solver->p[icv0]);

    my_psum += solver->p[icv0]*zone_ptr->area_bf[ibf] ;
  }

  MPI_Reduce(my_buf,buf,2,MPI_DOUBLE,MPI_MIN,0,mpi_comm);

  double pmean;
  MPI_Reduce(&my_psum,&pmean,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  pmean /= zone_ptr->area_global ;

  if ( mpi_rank == 0 ) {
    *ss << " > " << std::left << std::setw(30) << zone_ptr->getName() << " :  time, mean_p, min_p, max_p = " << std::right
        << std::setw(8) << solver->time << " " << std::setw(12) << pmean << " "
        << std::setw(12) << buf[0] << " " << std::setw(12) << -buf[1] << endl;
  }
  flush();

}

void OutletVVBc::initData() {
  
  assert(q_bf == NULL); q_bf = new double[zone_ptr->nbf];
  
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
    q_bf[ibf] = 0.0;
  }

}

void OutletVVBc::setBc() {

  // do nothing...

}

void OutletVVBc::updateBc() {

  double my_buf  = 0.0; 
  double un_conv = 0.0; 

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf]; 
    my_buf += DOT_PRODUCT(solver->u[icv0], zone_ptr->n_bf[ibf]) ; 
  } 

  MPI_Allreduce(&my_buf,&un_conv,1,MPI_DOUBLE,MPI_SUM,mpi_comm); 
  //un_conv = max(0.0,un_conv/solver->sum_outlet_proj_area); 
  un_conv = max(0.0,un_conv/zone_ptr->area_global) ;

  //  if ( mpi_rank == 0) { 
  //    cout << " outlet un_conv = " << un_conv << endl;
  // }

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
    q_bf[ibf] = un_conv*MAG(zone_ptr->n_bf[ibf]);
  }
  // dumpRange(q_bf,zone_ptr->nbf,"q_bf outlet");

}

void OutletVVBc::addFlux(double * rhs) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    //q_bf[ibf] = DOT_PRODUCT(solver->u[icv0], zone_ptr->n_bf[ibf]);
    rhs[icv0] += q_bf[ibf];
  }

}


void OutletVVBc::addMassFlux(double *rhs_rhoY0, double *rhs_rhoY1) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
   
    double rho0_flux = 0.0;
    double rho1_flux = 0.0;

    const double rho0_f = solver->rho0[icv0];
    const double rho1_f = solver->rho1[icv0];
    const double vof_f  = solver->vof[icv0];

    rhs_rhoY0[icv0] -= solver->rhoY0[icv0]*q_bf[ibf];
    rhs_rhoY1[icv0] -= (solver->rho[icv0] - solver->rhoY0[icv0])*q_bf[ibf];
  }

}



void OutletVVBc::addMomentumFlux(double * A, double (*rhs)[3]) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    const int coc00 = solver->cvocv_i[icv0];
    const double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])* zone_ptr->area_over_delta_bf[ibf];

    // first order upwind the outlet
    assert(q_bf[ibf] >= 0.0) ; //no backflow...
    FOR_I3 rhs[icv0][i] -= solver->rho[icv0]*q_bf[ibf]*solver->u[icv0][i];  

  }

}


void OutletVVBc::force_bf(double (*f_bf)[9]) const {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

    const int icv0 = zone_ptr->cvobf[ibf];
    const double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])* zone_ptr->area_over_delta_bf[ibf];
    FOR_I3 { 
      f_bf[ibf][3*0+i] = solver->rho[icv0]*q_bf[ibf]*solver->u[icv0][i];
      f_bf[ibf][3*1+i] = solver->p[icv0]*zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = mu_coeff*solver->u[icv0][i];
    }
  }

}

inline void reportSpongeErrorHBc(const string msg) {

  CERR( " SPONGE bc syntax error: \n" \
        << msg << "\n" \
        << " SPONGE TYPE  <L_X,L_Y,L_Z,L_FR> LENGTH <l> SPEED <u> STRENGTH <s> OUTLET_PRESSURE <p> OUTLET_STRENGTH <sigma>" );

}

inline void parseSpongeParamHBc(Param* param, SpongeType& sponge_type,
                                double& sponge_length, double& sponge_speed,
                                double& sponge_strength,
                                double& outlet_pressure, double& outlet_strength,
                                double*& xr, double*& nr) {

  int iarg = 1;
  while ( iarg < param->size()) {
    string token = param->getString(iarg++);
    if ( token == "SPEED" ) {
      sponge_speed = param->getDouble(iarg++);
    } else if ( token == "STRENGTH") {
      sponge_strength = param->getDouble(iarg++);
    } else if ( token == "TYPE" ) {
      token = param->getString(iarg++);
      if ( token == "L_X") {
        sponge_type = SPONGE_LX;
      } else if ( token == "L_Y") {
        sponge_type = SPONGE_LY;
      } else if ( token == "L_Z") {
        sponge_type = SPONGE_LZ;
      } else if ( token == "L_TH" ) {
        sponge_type = SPONGE_LTH;
        assert(!xr&&!nr);
        xr = new double[3];
        nr = new double[3];
        FOR_I3 xr[i] = param->getDouble(iarg++);
        FOR_I3 nr[i] = param->getDouble(iarg++);
        double mag_nr = MAG(nr);
        if (mag_nr<=0.0){
          CERR("SPONGE .. TYPE L_TH rotation axis must have non-zero magnitude");
        }
        FOR_I3 nr[i] /= mag_nr;
      } else {
        CERR( " > unrecognized sponge type : " << token);
      }
    } else if ( (token == "LENGTH")||(token == "OUTLET_LENGTH") ) {
      sponge_length = param->getDouble(iarg++);
    } else if ( token == "OUTLET_PRESSURE") {
      outlet_pressure = param->getDouble(iarg++);
    } else if ( token == "OUTLET_STRENGTH") {
      outlet_strength = param->getDouble(iarg++);
    } else {
      CERR( " > unrecognized sponge token : " << token);
    }
  }

  // error check the form of the boundary condition...

  if (sponge_length <= 0.0) {
    reportSpongeErrorHBc(" > LENGTH: nonpositive value specified or unspecified");
  } else if (sponge_strength < 0.0) {
    reportSpongeErrorHBc(" > STRENGTH: negative value specified or unspecified");
  } else if (outlet_strength < 0.0) {
    reportSpongeErrorHBc(" > OUTLET_STRENGTH: negative value specified or unspecified");
  } else if (outlet_pressure < 0.0) {
    reportSpongeErrorHBc(" > OUTLET_PRESSURE: negative value specified or unspecified");
  } else if (sponge_speed < 0.0) {
    reportSpongeErrorHBc(" > SPEED: negative value specified or unspecified");
  } else if (sponge_type == SPONGE_UNDEFINED) {
    reportSpongeErrorHBc(" > TYPE: unspecified");
  }

}

void SpongeVVBc::initData() {

  Param* param = getParam(getName());

  // default params...
  //Note L is used to set both sponge length and outlet length scale,
  //adjust sponge_strength to adjust outlet sigma/L ratio.
  //Also, sponge does not parse a Mach number, relying on OutletVV constructor
  //initialization of Ma=0.
  parseSpongeParamHBc(param,sponge_type,L,sponge_speed,sponge_strength,p_ref,sigma,xr,nr);
  initSpongeTypeTheta();

  // use nodes to get the domain limits...

  double my_buf[5] = { -1.0E-20, -1.0E-20, -1.0E-20, -1.0E-20, -1.0E-20 };
  for (int ino = 0; ino < solver->nno; ++ino) {
    my_buf[0] = max(my_buf[0],solver->x_no[ino][0]);
    my_buf[1] = max(my_buf[1],solver->x_no[ino][1]);
    my_buf[2] = max(my_buf[2],solver->x_no[ino][2]);
    my_buf[3] = max(my_buf[3],sqrt(solver->x_no[ino][1]*solver->x_no[ino][1] + solver->x_no[ino][2]*solver->x_no[ino][2]));
    my_buf[4] = max(my_buf[4],sqrt(solver->x_no[ino][0]*solver->x_no[ino][0] + solver->x_no[ino][1]*solver->x_no[ino][1]));
  }
  MPI_Allreduce(my_buf,sponge_data,5,MPI_DOUBLE,MPI_MAX,mpi_comm);

  assert(q_bf == NULL); q_bf = new double[zone_ptr->nbf];
  
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
    q_bf[ibf] = 0.0;
  }

}

void SpongeVVBc::setBc() {

  // do nothing...

}

void SpongeVVBc::updateBc() {

  double my_buf  = 0.0; 
  double un_conv = 0.0; 

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf]; 
    my_buf += DOT_PRODUCT(solver->u[icv0], zone_ptr->n_bf[ibf]) ; 
  } 

  MPI_Allreduce(&my_buf,&un_conv,1,MPI_DOUBLE,MPI_SUM,mpi_comm); 
  //un_conv = max(0.0,un_conv/solver->sum_outlet_proj_area); 
  un_conv = max(0.0,un_conv/zone_ptr->area_global) ;

  //  if ( mpi_rank == 0) { 
  //    cout << " outlet un_conv = " << un_conv << endl;
  // }

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
    q_bf[ibf] = un_conv*MAG(zone_ptr->n_bf[ibf]);
  }
  // dumpRange(q_bf,zone_ptr->nbf,"q_bf outlet");

}

void SpongeVVBc::addFlux(double * rhs) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    //q_bf[ibf] = DOT_PRODUCT(solver->u[icv0], zone_ptr->n_bf[ibf]);
    rhs[icv0] += q_bf[ibf];
  }

}


void SpongeVVBc::addMassFlux(double *rhs_rhoY0, double *rhs_rhoY1) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
   
    double rho0_flux = 0.0;
    double rho1_flux = 0.0;

    const double rho0_f = solver->rho0[icv0];
    const double rho1_f = solver->rho1[icv0];
    const double vof_f  = solver->vof[icv0];

    rhs_rhoY0[icv0] -= solver->rhoY0[icv0]*q_bf[ibf];
    rhs_rhoY1[icv0] -=(solver->rho[icv0] - solver->rhoY0[icv0])*q_bf[ibf];
  }

}



void SpongeVVBc::addMomentumFlux(double * A, double (*rhs)[3]) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    const int coc00 = solver->cvocv_i[icv0];
    const double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])* zone_ptr->area_over_delta_bf[ibf];

    // first order upwind the outlet
    assert(q_bf[ibf] >= 0.0) ; //no backflow...
    FOR_I3 rhs[icv0][i] -= solver->rho[icv0]*q_bf[ibf]*solver->u[icv0][i];  

  }

  if (sponge_type != SPONGE_LTH) {
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
    
    double sponge_velocity[3] = {0.0,0.0,0.0};
    sponge_velocity[idir] = sponge_speed;

    for (int icv = 0; icv < solver->ncv; ++icv) {
      if ( solver->x_cv[icv][idir] > Lmax-L) {
        const double x_sp = (solver->x_cv[icv][idir] - (Lmax-L))/L;
        // sponge profile = a*x^2 + b*x^8..
        // const double sponge_coeff = sponge_strength*(0.068*x_sp*x_sp + 0.845*pow(x_sp,8));
        // quadratic sponge profile (see JCP 212 681-702 (2006) & JCP 231 704-716 (2012))
        const double sponge_coeff = sponge_strength*x_sp*x_sp;
        const int coc00 = solver->cvocv_i[icv];
        
        A[coc00] += solver->vol_cv[icv]*sponge_coeff*solver->rho[icv];
        for (int i =0; i < 3; ++i)
          rhs[icv][i] += solver->vol_cv[icv]*sponge_coeff*solver->rho[icv]*sponge_velocity[i];
      }
    }
  }
  else if (sponge_type == SPONGE_LTH) {
    for (int icv = 0; icv < solver->ncv; ++icv) {
      double r_cv[3];
      FOR_I3 r_cv[i] = solver->x_cv[icv][i] - xr[i];
      double re0 = DOT_PRODUCT(r_cv,th_e0);
      double re1 = DOT_PRODUCT(r_cv,th_e1);
      double theta = atan2(re1,re0);
      if (theta>0&&theta<L){
        const double x_sp = (L-theta)/L;
        // sponge profile = a*x^2 + b*x^8..
        // const double sponge_coeff = sponge_strength*(0.068*x_sp*x_sp + 0.845*pow(x_sp,8));
        // quadratic sponge profile (see JCP 212 681-702 (2006) & JCP 231 704-716 (2012))
        const double sponge_coeff = sponge_strength*x_sp*x_sp;
        
        double sponge_velocity[3] = CROSS_PRODUCT(nr,r_cv);
        FOR_I3 sponge_velocity[i] *=sponge_speed*2.0*M_PI;
        
        const int coc00 = solver->cvocv_i[icv];
        //A[coc00] += solver->vol_cv[icv]*sponge_coeff*rho;
        for (int i =0; i < 3; ++i)
          rhs[icv][i] += solver->vol_cv[icv]*sponge_coeff*solver->rho[icv]*sponge_velocity[i];
      }
    }
  }
  else {
    assert(0);
  }
  

}


void SpongeVVBc::force_bf(double (*f_bf)[9]) const {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

    const int icv0 = zone_ptr->cvobf[ibf];
    const double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])* zone_ptr->area_over_delta_bf[ibf];
    FOR_I3 { 
      f_bf[ibf][3*0+i] = solver->rho[icv0]*q_bf[ibf]*solver->u[icv0][i];
      f_bf[ibf][3*1+i] = solver->p[icv0]*zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = mu_coeff*solver->u[icv0][i];
    }
  }

}


void SpongeVVBc::initSpongeTypeTheta(){
  if (sponge_type == SPONGE_LTH){
    //User should specify the xr, nr such that the cross product
    //with face based normals is generally positive. I.E.
    // nx X (x_bf-xr) dot (n_bf) > 0.  This is used to regularize
    // the mass flux out of the boundary face.
    //
    //Find the "greatest" rbf, ie the maximum theta assuming theta
    //increases about the negative axis of nr
    double nnr[3] = {-nr[0],-nr[1],-nr[2]};

    double rbf0[3];
    int b_hasFaces = 1;
    FOR_RANK {
      if (mpi_rank==rank){
        if ( zone_ptr->nbf > 0){
          if (b_hasFaces==1){
            FOR_I3 rbf0[i] = zone_ptr->x_bf[0][i] - xr[i];
            b_hasFaces = 0;
          }

          for (int ibf = b_hasFaces; ibf < zone_ptr->nbf; ++ibf) {
            double rbf[3];
            FOR_I3 rbf[i] = zone_ptr->x_bf[ibf][i] - xr[i];
            double cpr[3] = CROSS_PRODUCT(rbf0,rbf);
            if (DOT_PRODUCT(cpr,nnr)<0.0){
              FOR_I3 rbf0[i] = rbf[i];
              b_hasFaces = 0;
            }
          }
        }
      }
      MPI_Bcast(&b_hasFaces, 1, MPI_INT, rank, mpi_comm);
      if (b_hasFaces==0)
        MPI_Bcast(&rbf0, 3, MPI_DOUBLE, rank, mpi_comm);
    }
    double rbf0_mag = MAG(rbf0);

    FOR_I3 th_e0[i] = (rbf0[i] - DOT_PRODUCT(rbf0,nnr)*nnr[i])/rbf0_mag;
    double e1[3] = CROSS_PRODUCT(nnr,th_e0);
    FOR_I3 th_e1[i] = e1[i];
    FOR_I3 rbf0[i] += xr[i];
    if (mpi_rank==0){
      cout << "SPONGE TYPE L_TH Information:" << endl;
      cout << "       origin: " << COUT_VEC(xr) << endl;
      cout << "     cyl_axis: " << COUT_VEC(nnr) << endl;
      cout << " theta00_axis: " << COUT_VEC(th_e0) << endl;
      cout << " theta90_axis: " << COUT_VEC(th_e1) << endl;
      cout << "     bf_ref_x: " << COUT_VEC(rbf0)  << endl;
    }
  }
}



//----------------
void OutletVBc::initData() {
  
  assert(q_bf == NULL); q_bf = new double[zone_ptr->nbf];
  
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
    q_bf[ibf] = 0.0;
  }

}

void OutletVBc::setBc() {

  // do nothing...

}

void OutletVBc::updateBc() {

  // re-compute the outlet q_bf... 

  double my_buf  = 0.0; 
  double un_conv = 0.0; 

  //for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
  //  const int icv0 = zone_ptr->cvobf[ibf]; 
  //  my_buf += DOT_PRODUCT(solver->u[icv0], zone_ptr->n_bf[ibf]) ; 
  //} 
  for (int icv = 0; icv < solver->ncv; ++icv) 
    my_buf -= solver->div[icv];

  MPI_Allreduce(&my_buf,&un_conv,1,MPI_DOUBLE,MPI_SUM,mpi_comm); 
  un_conv = max(0.0,un_conv/solver->sum_outlet_proj_area); 

//  if ( mpi_rank == 0) { 
//    cout << " outlet un_conv = " << un_conv << endl;
// }

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
    q_bf[ibf] = un_conv*MAG(zone_ptr->n_bf[ibf]);
  }
 // dumpRange(q_bf,zone_ptr->nbf,"q_bf outlet");

}

void OutletVBc::addFlux(double * rhs) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    //q_bf[ibf] = DOT_PRODUCT(solver->u[icv0], zone_ptr->n_bf[ibf]);
    rhs[icv0] += q_bf[ibf];
  }

}


void OutletVBc::addMassFlux(double *rhs_rhoY0, double *rhs_rhoY1) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
   
    double rho0_flux = 0.0;
    double rho1_flux = 0.0;

    const double rho0_f = solver->rho0[icv0];
    const double rho1_f = solver->rho1[icv0];
    const double vof_f  = solver->vof[icv0];

    /*
    if (solver->rho0[icv0]*solver->vof[icv0]*solver->vol_cv[icv0] >= solver->rho0[icv0]*q_bf[ibf]*solver->dt ) 
      rho0_flux = solver->rho0[icv0]*solver->vof[icv0]*q_bf[ibf];
    else
      rho0_flux = solver->rho0[icv0]*solver->vof[icv0]*solver->vol_cv[icv0]/solver->dt;

    if (solver->rho1[icv0]*(1.0-solver->vof[icv0])*solver->vol_cv[icv0] >= solver->rho1[icv0]*q_bf[ibf]*solver->dt ) 
      rho1_flux = solver->rho1[icv0]*(1.0-solver->vof[icv0])*q_bf[ibf];
    else
      rho1_flux = solver->rho1[icv0]*(1.0-solver->vof[icv0])*solver->vol_cv[icv0]/solver->dt;
      */
    // just upwind for now
    //rhs_rhoY0[icv0] -= rho0_f*vof_f*q_bf[ibf];
    //rhs_rhoY1[icv0] -= rho1_f*(1.0-vof_f)*q_bf[ibf];
    rhs_rhoY0[icv0] -= solver->rhoY0[icv0]*q_bf[ibf];
    rhs_rhoY1[icv0] -= (solver->rho[icv0] - solver->rhoY0[icv0])*q_bf[ibf];
  }

}



void OutletVBc::addMomentumFlux(double * A, double (*rhs)[3]) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    const int coc00 = solver->cvocv_i[icv0];
    const double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])* zone_ptr->area_over_delta_bf[ibf];

    // first order upwind the outlet
    assert(q_bf[ibf] >= 0.0) ; //no backflow...
    FOR_I3 rhs[icv0][i] -= solver->rho[icv0]*q_bf[ibf]*solver->u[icv0][i];  

/*
    if (solver->vof[icv0]*solver->vol_cv[icv0] >= q_bf[ibf]*solver->dt ) {
      double rho_f = solver->rho_l;
      FOR_I3 rhs[icv0][i] -= solver->rho[icv0]*q_bf[ibf]*solver->u[icv0][i];  
    }
    else {
      double vof_f = solver->vof[icv0]*solver->vol_cv[icv0]/q_bf[ibf]/solver->dt;
      assert(vof_f <= 1.0);
      double rho_f = solver->rho_l*vof_f + solver->rho_g*(1.0-vof_f);
      FOR_I3 rhs[icv0][i] -= rho_f*q_bf[ibf]*solver->u[icv0][i];  
    }
 */  
    //A[coc00] += solver->rho[icv0]*q_bf[ibf];                                

    //FOR_I3 rhs[icv0][i] -= 0.5*solver->rho[icv0]*q_bf[ibf]*solver->u[icv0][i]; 
    //A[coc00] += 0.5*solver->rho[icv0]*q_bf[ibf];
  }

}


void OutletVBc::force_bf(double (*f_bf)[9]) const {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

    const int icv0 = zone_ptr->cvobf[ibf];
    const double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])* zone_ptr->area_over_delta_bf[ibf];
    FOR_I3 { 
      f_bf[ibf][3*0+i] = solver->rho[icv0]*q_bf[ibf]*solver->u[icv0][i];
      f_bf[ibf][3*1+i] = solver->p[icv0]*zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = mu_coeff*solver->u[icv0][i];
    }
  }

}

void AlgebraicWallModelVBc::initData() {

  assert( solver != NULL);
  assert(u_bc == NULL);  u_bc = new double[zone_ptr->nbf][3];
  assert(tau_wall == NULL); tau_wall = new double[zone_ptr->nbf];

  double uwall[3] = {0.0, 0.0, 0.0};

  bool isTranslating = false;
  Param * param = getParam(getName());
  int iarg = 1;
  while(iarg<param->size()){
    string token = param->getString(iarg++);
    if (token=="U_WALL"){
      uwall[0]  = param->getDouble(iarg++);
      uwall[1]  = param->getDouble(iarg++);
      uwall[2]  = param->getDouble(iarg++);
      isTranslating = true;
    }
    else if (token=="ROUGHNESS_HEIGHT"){
      z0 =  param->getDouble(iarg++);
      if (z0<=0){
        CERR("invalid roughness height " << z0 << " for WM_ALG_WALL, must be greater than zero.");
      }
      COUT1("Found roughness height " << z0 << " for " << getName() << " WM_ALG_WALL, using rough wall algebraic closure");
    }
  }

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    FOR_I3 u_bc[ibf][i] = uwall[i];
    tau_wall[ibf] = 0.0; //set to zero for now, requires solver->u be set for a better guess
                         //will to this in initialHook
  }

}



void AlgebraicWallModelVBc::initialHook(){
    //set initial guess for tau_wall, if it wasn't read in
  if (!zone_ptr->checkDataFlag("tau_wall")){
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv0 = zone_ptr->cvobf[ibf];
      double du[3]    = DIFF(solver->u[icv0],u_bc[ibf]); //velocity relative to moving wall
      double mag_nbf  = MAG(zone_ptr->n_bf[ibf]);
      double du_dot_nbf = DOT_PRODUCT(du,zone_ptr->n_bf[ibf])/mag_nbf;
      FOR_I3 du[i]   -= du_dot_nbf*zone_ptr->n_bf[ibf][i]/mag_nbf; //velocity tangent to wall
      double du_mag   = MAG(du);
      double mu = getDoubleParam("MU_REF1");
      tau_wall[ibf] = mu*du_mag*zone_ptr->area_over_delta_bf[ibf]/zone_ptr->area_bf[ibf];
    }
  }
}


void AlgebraicWallModelVBc::addMomentumFlux(double * A,double (*rhs)[3]) const {


  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    const int coc00 = solver->cvocv_i[icv0];

    //double u_bc_mag = MAG(u_bc[ibf]);
    //double du[3]    = DIFF(solver->u[icv],u_bc[ibf]); //velocity relative to moving wall
    //double du_mag   = DOT_PRODUCT(u_diff,u_bc[ibf])/u_bc_mag; //velocity magnitude parallel to wall

    double du[3]    = DIFF(solver->u[icv0],u_bc[ibf]); //velocity relative to moving wall
    double mag_nbf  = MAG(zone_ptr->n_bf[ibf]);
    double du_dot_nbf = DOT_PRODUCT(du,zone_ptr->n_bf[ibf])/mag_nbf;
    FOR_I3 du[i]   -= du_dot_nbf*zone_ptr->n_bf[ibf][i]/mag_nbf; //velocity tangent to wall
    double du_mag   = MAG(du);

    double u_tau    = sqrt(tau_wall[ibf]/solver->rho[icv0]); //initial guess
    if (z0<0){//smooth wall
      tau_wall[ibf] = AlgebraicWM::solve_tau(du_mag, zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf], solver->rho[icv0], solver->mu_lam[icv0], u_tau);
    }
    else{//rough wall
      tau_wall[ibf] = AlgebraicRoughnessWM::solve_tau(du_mag,zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf], z0, solver->rho[icv0], solver->mu_lam[icv0],u_tau);
    }
    //TODO keep mu_sgs??
    const double mu_coeff = (tau_wall[ibf]*zone_ptr->area_bf[ibf]/du_mag);
    A[coc00] += mu_coeff;
    FOR_I3 rhs[icv0][i] += mu_coeff*u_bc[ibf][i];
  }

}



void AlgebraicWallModelVBc::force_bf(double (*f_bf)[9]) const {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];

    double du[3]     = DIFF(solver->u[icv0],u_bc[ibf]); //velocity relative to moving wall
    double du_mag    = MAG(du);
    FOR_I3 du[i]    /= du_mag;
    double mag_nbf   = MAG(zone_ptr->n_bf[ibf]);
    double unit_n[3];
    FOR_I3 unit_n[i] = zone_ptr->n_bf[ibf][i]/mag_nbf;
    double du_dot_n  = DOT_PRODUCT(du,unit_n);
    FOR_I3 du[i]    -= du_dot_n*unit_n[i];  //wall parallel unit vector

    //TODO recompute tau_wall here for current time level force??

    FOR_I3 {
      f_bf[ibf][3*0+i] = 0.0;
      f_bf[ibf][3*1+i] = solver->p[icv0]*zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = tau_wall[ibf]*du[i]*mag_nbf; //shear stress * wall parallel area vector
    }
  }
  //addPressureForceDueToFrameRotation(f_bf);
}

void AlgebraicWallModelVBc::query(const string& param_str) {

  // report min/max slip length and velocity
  double my_buf[3] = {0.0,0.0,0.0}; // area, int tau_wall dA, y_plus
  double buf[3];

  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];


    const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
    const double u_tau = sqrt( tau_wall[ibf]/solver->rho[icv0]);
    const double nu    = solver->mu_lam[icv0] /solver->rho[icv0];
    const double y_plus = y1*u_tau/nu;

    my_buf[0] += zone_ptr->area_bf[ibf];
    my_buf[1] += zone_ptr->area_bf[ibf]*tau_wall[ibf];
    my_buf[2] += zone_ptr->area_bf[ibf]*y_plus;
  }

  MPI_Reduce(my_buf,buf,3,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if (mpi_rank == 0) {
    *ss << " > " << std::left << std::setw(30) << zone_ptr->getName() << " :  time, int tau_wall dA, avg y_plus = " << std::right
        << std::setw(8) << solver->time << " "
        << std::setw(12) << buf[1] << " "
        << std::setw(12) << buf[2]/buf[0] << endl;
  }
  flush();

}

void HookVBc::initData() {
 
  assert(u_bc   == NULL); u_bc   = new double[zone_ptr->nbf][3];
  assert(q_bf   == NULL); q_bf   = new double[zone_ptr->nbf];

  // just constant for now...
  double u_in[3],vof_in;
  double u_normal;
  Param * param = getParam(getName());
  if (param->size() <= 5) { // inlet_normal ... 
    COUT1("INLET_NORMAL:");
    rho0_bc  =  param->getDouble(1);
    rho1_bc  =  param->getDouble(2);
    vof_bc   =  param->getDouble(3);
    u_normal =  param->getDouble(4);
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
      double normal[3];
      const double mag = MAG(zone_ptr->n_bf[ibf]);
      assert(mag > 0.0);
      FOR_I3 normal[i] = zone_ptr->n_bf[ibf][i]/mag;
      //boundary normal vector is out-going...
      FOR_I3 u_bc[ibf][i] = -u_normal*normal[i];
      q_bf[ibf] = DOT_PRODUCT(zone_ptr->n_bf[ibf],u_bc[ibf]);
    }
  }
  else {
    rho0_bc = param->getDouble(1);
    rho1_bc = param->getDouble(2);
    vof_bc  = param->getDouble(3);
    u_in[0] = param->getDouble(4);
    u_in[1] = param->getDouble(5);
    u_in[2] = param->getDouble(6);
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
      FOR_I3 u_bc[ibf][i] = u_in[i];
      q_bf[ibf] = DOT_PRODUCT(zone_ptr->n_bf[ibf],u_bc[ibf]);
    }
  }

  dumpRange(q_bf,zone_ptr->nbf,"q_bf inlet");

}

void HookVBc::setBc() {

  // this routine would be non-empty if there was time variation in the inlet bc... 
  const int _step = solver->step;
  const int _check_interval = solver->check_interval;

  double u_in = 0.00570;
  u_in = u_in + u_in*0.01*(_step - 40390);
  if (u_in > 0.057) u_in = 0.057;
  if (mpi_rank == 0 && _step%_check_interval == 0) cout << "u_in = " << u_in << endl;

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
    u_bc[ibf][0] = u_in;
    u_bc[ibf][1] = 0.0;
    u_bc[ibf][2] = 0.0;
    q_bf[ibf] = DOT_PRODUCT(zone_ptr->n_bf[ibf],u_bc[ibf]);
  }

}

void HookVBc::addFlux(double * rhs) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    rhs[icv0] += q_bf[ibf];
  }

}


void HookVBc::addMassFlux(double *rhs_rhoY0, double *rhs_rhoY1) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    rhs_rhoY0[icv0] -= rho0_bc*vof_bc*q_bf[ibf];
    rhs_rhoY1[icv0] -= rho1_bc*(1.0-vof_bc)*q_bf[ibf];
  }
}


void HookVBc::addMomentumFlux(double * A, double (*rhs)[3]) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    const int coc00 = solver->cvocv_i[icv0];
    const double rho_bc = rho1_bc*(1.0-vof_bc) + rho0_bc*vof_bc;
    const double mu_coeff = ((solver->mu1_ref*(1.0-vof_bc)+solver->mu0_ref*vof_bc) + solver->mu_sgs[icv0])*zone_ptr->area_over_delta_bf[ibf];
    A[coc00] += mu_coeff;
    FOR_I3 rhs[icv0][i] += (-rho_bc*q_bf[ibf] + mu_coeff)*u_bc[ibf][i];

    /*
    if (q_bf[ibf] >= 0.0) {
      // outflow -- use first-order upwind...
      //A[coc00] += rho_bc*q_bf[ibf] + mu_coeff;
      A[coc00] += mu_coeff;
      FOR_I3 rhs[icv0][i] += mu_coeff*u_bc[ibf][i];
    }
    else {
      // inflow...
      A[coc00] += mu_coeff;
      FOR_I3 rhs[icv0][i] += (-rho_bc*q_bf[ibf] + mu_coeff)*u_bc[ibf][i];

      //A[coc00] += 0.5*rho_bc*q_bf[ibf]+mu_coeff;
      //FOR_I3 rhs[icv0][i] += (-0.5*rho_bc*q_bf[ibf] + mu_coeff)*u_bc[ibf][i];
    }
    */
  }

}

void HookVBc::force_bf(double (*f_bf)[9]) const {


  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

    const int icv0 = zone_ptr->cvobf[ibf];
    const double rho_bc = rho1_bc*(1.0-vof_bc) + rho0_bc*vof_bc;
    const double mu_coeff = ((solver->mu1_ref*(1.0-vof_bc)+solver->mu0_ref*vof_bc) + solver->mu_sgs[icv0])*zone_ptr->area_over_delta_bf[ibf];
    FOR_I3 { 
      f_bf[ibf][3*0+i] = rho_bc*q_bf[ibf]*u_bc[ibf][i];
      f_bf[ibf][3*1+i] = solver->p[icv0]*zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = mu_coeff*(solver->u[icv0][i] - u_bc[ibf][i]);
    }
  }

}

void HookVBc::query(const string& param_str) {


  // report min/max
  double my_buf[2] = {HUGE_VAL,HUGE_VAL};
  double buf[2];

  double my_psum = 0.0;

  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    my_buf[0] = min(my_buf[0],solver->p[icv0]);
    my_buf[1] = min(my_buf[1],-solver->p[icv0]);

    my_psum += solver->p[icv0]*zone_ptr->area_bf[ibf] ;
  }

  MPI_Reduce(my_buf,buf,2,MPI_DOUBLE,MPI_MIN,0,mpi_comm);

  double pmean;
  MPI_Reduce(&my_psum,&pmean,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  pmean /= zone_ptr->area_global ;

  if ( mpi_rank == 0 ) {
    *ss << " > " << std::left << std::setw(30) << zone_ptr->getName() << " :  time, mean_p, min_p, max_p = " << std::right
        << std::setw(8) << solver->time << " " << std::setw(12) << pmean << " "
        << std::setw(12) << buf[0] << " " << std::setw(12) << -buf[1] << endl;
  }
  flush();

}

CtiRegister::CtiDataError AlgebraicWallModelVBc::funcEvalCtiData(CtiRegister::CtiData& v,
    const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if ( zone_ptr->isLeadingStrMatched(name) ) {

    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix.  we have to check for the actual data now..

    const string tau_str      = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str       = zone_ptr->getName() + ":" + "y_plus";
    const string bl_delta_str = zone_ptr->getName() + ":" + "bl_delta";

    // all of the following functions do not take any arguments -- check first

    if ( name == tau_str) {

      double * tmp = zone_ptr->createBfD1Data(v);

      // if stats have been requested on this quantity, this function will be called
      // but its possible that the mu_lam, rho necessary have not been populated yet.
      // check the data flags for these values -- revisit if this impacts performance
      // during the solution run-time.

      if ( b_eval_func) {

        if (args.size() == 0) {
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            tmp[ibf]                = tau_wall[ibf];
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
            const int icv0 = zone_ptr->cvobf[ibf];

            double du[3]     = DIFF(solver->u[icv0],u_bc[ibf]); //velocity relative to moving wall
            double du_mag    = MAG(du);
            FOR_I3 du[i]    /= du_mag;
            double mag_nbf   = MAG(zone_ptr->n_bf[ibf]);
            double unit_n[3];
            FOR_I3 unit_n[i] = zone_ptr->n_bf[ibf][i]/mag_nbf;
            double du_dot_n  = DOT_PRODUCT(du,unit_n);
            FOR_I3 du[i]    -= du_dot_n*unit_n[i];  //wall parallel unit vector

            tmp[ibf] = tau_wall[ibf]*DOT_PRODUCT(dir,du);
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
          tmp[ibf]                = y1* sqrt(tau_wall[ibf]*solver->rho[icv])/solver->mu_lam[icv];
        }

      }

      return CTI_DATA_OK;

    }
    else if ( name == bl_delta_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {

        //hard code bl inputs for now...
        int nbl = getIntParam("BL_LAYERS",10);
        double pt_inf = getDoubleParam("BL_PT_INF",0.50);
        double pt_cut = getDoubleParam("BL_PT_CUT",0.95);

        double lbl;
        if (checkParam("BL_DELTA"))
          lbl = getDoubleParam("BL_DELTA");
        else{
          lbl = nbl*zone_ptr->area_global/zone_ptr->area_over_delta_global;
          COUT1(" > BfZone " << getName() << " bl_delta: using mean delta global " << zone_ptr->area_global/zone_ptr->area_over_delta_global);
        }

        COUT1(" > BfZone " << getName() << " bl_delta: Layers " << nbl << ", Thickness " << lbl << ", Pt_threshold " << pt_inf*pt_cut);

        if (blde==NULL)
          blde = new BoundaryLayerDataExchanger<HelmholtzVofSolver>(zone_ptr,solver,nbl,lbl);

        blde->computeBlFromPt(tmp,pt_cut*pt_inf);

      }

      return CTI_DATA_OK;
    }
  }

  return VofBc::funcEvalCtiData(v,name,args,b_eval_func);
}
