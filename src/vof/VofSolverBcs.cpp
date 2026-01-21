#include "VofSolver.hpp"

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
  assert(vof_bc == NULL); vof_bc = new double[zone_ptr->nbf];
  assert(q_bf   == NULL); q_bf   = new double[zone_ptr->nbf];

  // just constant for now...
  double u_in[3],vof_in;
  double u_normal;
  Param * param = getParam(getName());
  if (param->size() <= 3) { // inlet_normal ... 
    COUT1("INLET_NORMAL:");
    u_normal = param->getDouble(1);
    vof_in =   param->getDouble(2);
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
      double normal[3];
      const double mag = MAG(zone_ptr->n_bf[ibf]);
      assert(mag > 0.0);
      FOR_I3 normal[i] = zone_ptr->n_bf[ibf][i]/mag;
      //boundary normal vector is out-going...
      FOR_I3 u_bc[ibf][i] = -u_normal*normal[i];
      vof_bc[ibf] = vof_in;
      q_bf[ibf] = DOT_PRODUCT(zone_ptr->n_bf[ibf],u_bc[ibf]);
    }
  }
  else {
    u_in[0] = param->getDouble(1);
    u_in[1] = param->getDouble(2);
    u_in[2] = param->getDouble(3);
    vof_in  = param->getDouble(4);
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
      FOR_I3 u_bc[ibf][i] = u_in[i];
      vof_bc[ibf] = vof_in;
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

void InletVBc::addVofFlux(double * rhs) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    rhs[icv0] += -vof_bc[ibf]*q_bf[ibf];
  }

}

void InletVBc::addMomentumFlux(double * A, double (*rhs)[3]) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    const int coc00 = solver->cvocv_i[icv0];
    const double rho_bc = solver->rho_g*(1.0-vof_bc[ibf])+solver->rho_l*vof_bc[ibf];
    const double mu_coeff = ((solver->mu_g*(1.0-vof_bc[ibf])+solver->mu_l*vof_bc[ibf]) + solver->mu_sgs[icv0])*zone_ptr->area_over_delta_bf[ibf];
    if (q_bf[ibf] >= 0.0) {
      // outflow -- use first-order upwind...
      A[coc00] += rho_bc*q_bf[ibf] + mu_coeff;
      FOR_I3 rhs[icv0][i] += mu_coeff*u_bc[ibf][i];
    }
    else {
      // inflow...
      A[coc00] += mu_coeff;
      FOR_I3 rhs[icv0][i] += (-rho_bc*q_bf[ibf] + mu_coeff)*u_bc[ibf][i];

      //A[coc00] += 0.5*rho_bc*q_bf[ibf]+mu_coeff;
      //FOR_I3 rhs[icv0][i] += (-0.5*rho_bc*q_bf[ibf] + mu_coeff)*u_bc[ibf][i];
    }
  }

}

void InletVBc::force_bf(double (*f_bf)[9]) const {


  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

    const int icv0 = zone_ptr->cvobf[ibf];
    const double rho_bc = solver->rho_g*(1.0-vof_bc[ibf])+solver->rho_l*vof_bc[ibf];
    const double mu_coeff = ((solver->mu_g*(1.0-vof_bc[ibf])+solver->mu_l*vof_bc[ibf]) + solver->mu_sgs[icv0])*zone_ptr->area_over_delta_bf[ibf];
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

  //if ( mpi_rank == 0) { 
  //  cout << " outlet un_conv = " << un_conv << endl;
 // }

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
    q_bf[ibf] = un_conv*MAG(zone_ptr->n_bf[ibf]);
  }
  //dumpRange(q_bf,zone_ptr->nbf,"q_bf outlet");

}

void OutletVBc::addFlux(double * rhs) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    rhs[icv0] += q_bf[ibf];
  }

}

void OutletVBc::addVofFlux(double * rhs) const { 

  assert( solver != NULL);
 
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    rhs[icv0] -= solver->vof[icv0]*q_bf[ibf];
   
   /* 
    if (solver->vof[icv0]*solver->vol_cv[icv0] >= q_bf[ibf]*solver->dt ) {
      rhs[icv0] -= solver->vof[icv0]*q_bf[ibf];
    }
    else {
      rhs[icv0] -= solver->vof[icv0]*solver->vol_cv[icv0]/solver->dt;
    }
   */ 
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
