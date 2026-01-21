#include "HelmholtzSolver.hpp"
#include "bcprofile/ProfileReader.hpp"
#include "wm/AlgebraicWM.hpp"
#include "wm/AlgebraicRoughnessWM.hpp"
#include "it/InflowTurbulenceH.hpp"

void HelmholtzBc::completePressureGrad(double (*dpdx)[3]) const {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0        = zone_ptr->cvobf[ibf];
    FOR_I3 dpdx[icv0][i] += zone_ptr->n_bf[ibf][i]*solver->p[icv0];
  }


  if (solver->frame_rotation){
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv0     = zone_ptr->cvobf[ibf];
      const double dx[3] = DIFF(zone_ptr->x_bf[ibf],solver->x_cv[icv0]);
      const double mag_n_bf = MAG(zone_ptr->n_bf[ibf]);
      if (mag_n_bf>0.0){
        double unit_n[3];
        FOR_I3 unit_n[i] = zone_ptr->n_bf[ibf][i]/mag_n_bf;

        //compute acceleration due to frame rotation at cv and at bf
        //

        // bf
        // term 1: -omega x omega x r
        // centrifugal force

        double r_bf[3];
        FOR_I3 r_bf[i] = zone_ptr->x_bf[ibf][i] - solver->frame_rotation[i+3];

        double coeff[3];
        coeff[0] = solver->frame_rotation[1]*r_bf[2] - solver->frame_rotation[2]*r_bf[1];
        coeff[1] = solver->frame_rotation[2]*r_bf[0] - solver->frame_rotation[0]*r_bf[2];
        coeff[2] = solver->frame_rotation[0]*r_bf[1] - solver->frame_rotation[1]*r_bf[0];

        double a_bf[3];
        a_bf[0] = -(solver->frame_rotation[1]*coeff[2] - solver->frame_rotation[2]*coeff[1]);
        a_bf[1] = -(solver->frame_rotation[2]*coeff[0] - solver->frame_rotation[0]*coeff[2]);
        a_bf[2] = -(solver->frame_rotation[0]*coeff[1] - solver->frame_rotation[1]*coeff[0]);

        // term 2: -2 omega x v
        // Coriolis force

        //TODO: changing this to -omega x v so when v=-w x r , dp = 0;
        //      why? convective terms balance half the coriolis force...

        //Assumes if u_bc is null then u_bc = u_cv (i.e. slip walls, outlets, ...others?)
        if (u_bc){
          a_bf[0] -= 1.0*(solver->frame_rotation[1]*u_bc[ibf][2] - solver->frame_rotation[2]*u_bc[ibf][1]);
          a_bf[1] -= 1.0*(solver->frame_rotation[2]*u_bc[ibf][0] - solver->frame_rotation[0]*u_bc[ibf][2]);
          a_bf[2] -= 1.0*(solver->frame_rotation[0]*u_bc[ibf][1] - solver->frame_rotation[1]*u_bc[ibf][0]);
        }
        else{
          a_bf[0] -= 1.0*(solver->frame_rotation[1]*solver->u[icv0][2] - solver->frame_rotation[2]*solver->u[icv0][1]);
          a_bf[1] -= 1.0*(solver->frame_rotation[2]*solver->u[icv0][0] - solver->frame_rotation[0]*solver->u[icv0][2]);
          a_bf[2] -= 1.0*(solver->frame_rotation[0]*solver->u[icv0][1] - solver->frame_rotation[1]*solver->u[icv0][0]);
        }

        //// use cv values for trapezoid rule estimation of p_bf.  For a hydrostatic pressure profile with isentropic based
        //// variation in density, the estimate of dpdx is better if just the bf values above are used to estimate dp.  The
        //// accuracy of p_bf at the face is worse but cell based dpdx prediction is better.

        //// cv
        //// term 1: -omega x omega x r
        //// centrifugal force
        //
        //double r_cv[3];
        //FOR_I3 r_cv[i] = solver->x_cv[icv0][i] - solver->frame_rotation[i+3];

        //coeff[0] = solver->frame_rotation[1]*r_cv[2] - solver->frame_rotation[2]*r_cv[1];
        //coeff[1] = solver->frame_rotation[2]*r_cv[0] - solver->frame_rotation[0]*r_cv[2];
        //coeff[2] = solver->frame_rotation[0]*r_cv[1] - solver->frame_rotation[1]*r_cv[0];
        //
        //double a_cv[3];
        //a_cv[0] = -(solver->frame_rotation[1]*coeff[2] - solver->frame_rotation[2]*coeff[1]);
        //a_cv[1] = -(solver->frame_rotation[2]*coeff[0] - solver->frame_rotation[0]*coeff[2]);
        //a_cv[2] = -(solver->frame_rotation[0]*coeff[1] - solver->frame_rotation[1]*coeff[0]);
        //
        //// term 2: -2 omega x v
        //// Coriolis force
        //
        //a_cv[0] -= 2.0*(solver->frame_rotation[1]*solver->u[icv0][2] - solver->frame_rotation[2]*solver->u[icv0][1]);
        //a_cv[1] -= 2.0*(solver->frame_rotation[2]*solver->u[icv0][0] - solver->frame_rotation[0]*solver->u[icv0][2]);
        //a_cv[2] -= 2.0*(solver->frame_rotation[0]*solver->u[icv0][1] - solver->frame_rotation[1]*solver->u[icv0][0]);

        //// add to momentum...
        //const double dp  = 0.5*solver->rho[icv0]*DOT_PRODUCT(dx,unit_n)*(DOT_PRODUCT(a_cv,unit_n)+DOT_PRODUCT(a_bf,unit_n));

        const double dp  = solver->rho[icv0]*DOT_PRODUCT(dx,unit_n)*(DOT_PRODUCT(a_bf,unit_n));
        FOR_I3 dpdx[icv0][i] += zone_ptr->n_bf[ibf][i]*dp;
      }
    }
  }
}

void HelmholtzBc::restrictBc(StaticSolver::CoarseGrid* cg,HelmholtzBc* bc) {

  if (mf) {
    cg->restrictExtrinsicCbData(mf,bc->mf,zone_ptr->index);

    //double my_buf = 0.0;
    //for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
    //  my_buf += mf[ibf+zone_ptr->ibf_f];
    //double buf;
    //MPI_Reduce(&my_buf,&buf,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    //if (mpi_rank == 0)
    //  cout << " mf outlet coarse: " << buf << endl;
  }
}

void HelmholtzBc::query(const string& param_str) const {
  CERR("Boundary query not defined for zone " << getName());
}

CtiRegister::CtiDataError HelmholtzBc::funcEvalCtiData(CtiRegister::CtiData& v,
                                                const string& name, list<CtiRegister::CtiData>& args,
                                                const bool b_eval_func) {

  if ( zone_ptr->isLeadingStrMatched(name)) {

    const string proj_str = this->getName() + ":" + "proj";
    const string xbf_str = this->getName() + ":" + "x_bf";
    const string pbf_str = this->getName() + ":" + "p_bf";

    if ( name == proj_str) {


      if ( args.size() != 1) {
        return CtiRegister::CTI_DATA_ARG_COUNT;
      }

      list<CtiRegister::CtiData>::iterator arg = args.begin();
      if ( arg->getTopology() != CV_DATA)
        return CtiRegister::CTI_DATA_NOT_VALID;

      if ( arg->getType() == DN_DATA) {

        double * tmp = zone_ptr->createBfD1Data(v);

        if ( b_eval_func) {

          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            const int icv = zone_ptr->cvobf[ibf];
            tmp[ibf]      = arg->dn(icv);

          }
        }

        return CtiRegister::CTI_DATA_OK;

      } else if ( arg->getType() == DN3_DATA) {

        double (*tmp)[3] = zone_ptr->createBfD2Data(v);

        if ( b_eval_func) {

          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            const int icv = zone_ptr->cvobf[ibf];
            for (int i = 0; i < 3; ++i)
              tmp[ibf][i] = arg->dn3(icv,i);

          }
        }

        return CtiRegister::CTI_DATA_OK;

      } else {

        return CtiRegister::CTI_DATA_NOT_VALID;

      }

    } else if ( name == xbf_str) {

      if ( args.size() != 0)
        return CtiRegister::CTI_DATA_ARG_COUNT;

      double (*tmp)[3] = zone_ptr->createBfD2Data(v);

      if ( b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          for (int i = 0; i < 3; ++i)
            tmp[ibf][i] = zone_ptr->x_bf[ibf][i];

        }
      }

      return CtiRegister::CTI_DATA_OK;

    } else if ( name == pbf_str) {

      if ( args.size() != 0)
        return CtiRegister::CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          const int icv0          = zone_ptr->cvobf[ibf];
          tmp[ibf]                = solver->p[icv0];
        }
      }
      return CtiRegister::CTI_DATA_OK;
    }

  }

  return CtiRegister::CTI_DATA_NOT_FOUND;
}


void HelmholtzBc::addPressureForceDueToFrameRotation(double (*f_bf)[9]) const {
  if (solver->frame_rotation){
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv0     = zone_ptr->cvobf[ibf];
      const double dx[3] = DIFF(zone_ptr->x_bf[ibf],solver->x_cv[icv0]);
      const double mag_n_bf = MAG(zone_ptr->n_bf[ibf]);
      if (mag_n_bf>0.0){
        double unit_n[3];
        FOR_I3 unit_n[i] = zone_ptr->n_bf[ibf][i]/mag_n_bf;

        //compute acceleration due to frame rotation at cv and at bf
        //
        // bf
        // term 1: -omega x omega x r
        // centrifugal force

        double r_bf[3];
        FOR_I3 r_bf[i] = zone_ptr->x_bf[ibf][i] - solver->frame_rotation[i+3];

        double coeff[3];
        coeff[0] = solver->frame_rotation[1]*r_bf[2] - solver->frame_rotation[2]*r_bf[1];
        coeff[1] = solver->frame_rotation[2]*r_bf[0] - solver->frame_rotation[0]*r_bf[2];
        coeff[2] = solver->frame_rotation[0]*r_bf[1] - solver->frame_rotation[1]*r_bf[0];

        double a_bf[3];
        a_bf[0] = -(solver->frame_rotation[1]*coeff[2] - solver->frame_rotation[2]*coeff[1]);
        a_bf[1] = -(solver->frame_rotation[2]*coeff[0] - solver->frame_rotation[0]*coeff[2]);
        a_bf[2] = -(solver->frame_rotation[0]*coeff[1] - solver->frame_rotation[1]*coeff[0]);

        // term 2: -2 omega x v
        // Coriolis force

        //TODO: changing this to -omega x v so when v=-w x r , dp = 0;
        //      why? convective terms balance half the coriolis force...

        //Assumes if u_bc is null then u_bc = u_cv (i.e. slip walls, outlets, ...others?)
        if (u_bc){
          a_bf[0] -= 1.0*(solver->frame_rotation[1]*u_bc[ibf][2] - solver->frame_rotation[2]*u_bc[ibf][1]);
          a_bf[1] -= 1.0*(solver->frame_rotation[2]*u_bc[ibf][0] - solver->frame_rotation[0]*u_bc[ibf][2]);
          a_bf[2] -= 1.0*(solver->frame_rotation[0]*u_bc[ibf][1] - solver->frame_rotation[1]*u_bc[ibf][0]);
        }
        else{
          a_bf[0] -= 1.0*(solver->frame_rotation[1]*solver->u[icv0][2] - solver->frame_rotation[2]*solver->u[icv0][1]);
          a_bf[1] -= 1.0*(solver->frame_rotation[2]*solver->u[icv0][0] - solver->frame_rotation[0]*solver->u[icv0][2]);
          a_bf[2] -= 1.0*(solver->frame_rotation[0]*solver->u[icv0][1] - solver->frame_rotation[1]*solver->u[icv0][0]);
        }

        // use cv values for trapezoid rule estimation of p_bf.  For a hydrostatic pressure profile with isentropic based
        // variation in density, the estimate of dpdx is better if just the bf values above are used to estimate dp.  The
        // accuracy of p_bf at the face is worse but cell based dpdx prediction is better.

        // cv
        // term 1: -omega x omega x r
        // centrifugal force

        double r_cv[3];
        FOR_I3 r_cv[i] = solver->x_cv[icv0][i] - solver->frame_rotation[i+3];

        coeff[0] = solver->frame_rotation[1]*r_cv[2] - solver->frame_rotation[2]*r_cv[1];
        coeff[1] = solver->frame_rotation[2]*r_cv[0] - solver->frame_rotation[0]*r_cv[2];
        coeff[2] = solver->frame_rotation[0]*r_cv[1] - solver->frame_rotation[1]*r_cv[0];

        double a_cv[3];
        a_cv[0] = -(solver->frame_rotation[1]*coeff[2] - solver->frame_rotation[2]*coeff[1]);
        a_cv[1] = -(solver->frame_rotation[2]*coeff[0] - solver->frame_rotation[0]*coeff[2]);
        a_cv[2] = -(solver->frame_rotation[0]*coeff[1] - solver->frame_rotation[1]*coeff[0]);

        // term 2: -2 omega x v
        // Coriolis force

        //TODO: changing this to -omega x v so when v=-w x r , dp = 0;
        //      why? convective terms balance half the coriolis force...

        a_cv[0] -= 1.0*(solver->frame_rotation[1]*solver->u[icv0][2] - solver->frame_rotation[2]*solver->u[icv0][1]);
        a_cv[1] -= 1.0*(solver->frame_rotation[2]*solver->u[icv0][0] - solver->frame_rotation[0]*solver->u[icv0][2]);
        a_cv[2] -= 1.0*(solver->frame_rotation[0]*solver->u[icv0][1] - solver->frame_rotation[1]*solver->u[icv0][0]);

        //estimate of pressure change from cv to bf due to frame rotation
        const double dp  = 0.5*solver->rho[icv0]*DOT_PRODUCT(dx,unit_n)*(DOT_PRODUCT(a_cv,unit_n)+DOT_PRODUCT(a_bf,unit_n));

        FOR_I3 {
          f_bf[ibf][3*1+i] += dp*zone_ptr->n_bf[ibf][i];
        }
      }
    }
  }
}

void HelmholtzBc::parseScalarBc(Param * param, const bool b_required) {

  assert ( solver != NULL);

  if ( solver->nsc_transport == 0 )
    return;

  int * sc_flag = new int[solver->nsc_transport]; //TODO consider making this global so for a wall bc multiple scalars can zero flux or const val. bcs.
  for (int isc = 0; isc < solver->nsc_transport; ++isc)
    sc_flag[isc] = 0;

  assert ( scbc_vals == NULL);
  int iarg = 0;
  while ( iarg < param->size() ) {

    const string token = param->getString(iarg++);

    const int isc  = solver->getScalarIndex(token);
    if ( isc >= 0) {

      if ( sc_flag[isc] == 1) {

        CERR(" > multiply specified scalar bcs for var: " << token);

      } else {

        if (!scbc_vals)
          scbc_vals = new double[solver->nsc_transport]; //only allocate if needed

        sc_flag[isc]   = 1;
        scbc_vals[isc] = param->getDouble(iarg++);

      }
    }
  }

  // check to make sure that all the scalars were set, if required

  int err_flag = 0;
  for (int isc = 0; isc < solver->nsc_transport; ++isc) {

    if ( sc_flag[isc] == 0) {

      if ( mpi_rank == 0 && (b_required||scbc_vals)  )
        cout << " > failed to provide scalar bc for : " << solver->getScalarName(isc)
          << "   at zone : " << zone_ptr->getName() << endl;

      err_flag = 1;

    } else {
      assert( sc_flag[isc] == 1);
      if (mpi_rank == 0) {
        cout << " > Found a scalar bc value of " << scbc_vals[isc] << " for : " << solver->getScalarName(isc)
          << "   at zone : " << zone_ptr->getName() << endl;
      }
    }
  }

  delete[] sc_flag;

  if ( err_flag == 1 && (b_required||scbc_vals)) {
    CERR( " > please correct boundary condition syntax for failed scalar vars above ");
  }

  if (!b_required&&scbc_vals==NULL&&mpi_rank==0){
    cout << " Using zero flux boundary condition for all scalars on zone " << zone_ptr->getName() << endl;
  }

}

void HelmholtzBc::_addScalarBoundaryFluxConvective(double* A,double* At,double * rhs,const ScalarField *scalar) const {

  assert ( solver != NULL);
  assert ( mf != NULL);
  assert ( scbc_vals != NULL);

  const int isc = scalar->idx;

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv = zone_ptr->cvobf[ibf];
    const int coc00 = solver->cvocv_i[icv];
    if (mf[ibf] >= 0.0) {
      // outflow -- use first-order upwind...
      if (A != NULL) A[coc00] += mf[ibf];
      if (At != NULL) At[coc00] += mf[ibf];
    }
    else {
      // inflow...
      rhs[icv] -= mf[ibf]*scbc_vals[isc];
    }
  }

}

void HelmholtzBc::_addScalarBoundaryFluxWall(double* A,double* At,double * rhs,const ScalarField *scalar) const {

  assert ( solver != NULL);

  const int isc = scalar->idx;
  const double Sc_lam = scalar->Sc_lam;
  const double Sc_t = scalar->Sc_t;

  if (scbc_vals!=NULL) { //this is a const. value wall BC
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      assert(mf==NULL||mf[ibf]==0.0);
      const int icv = zone_ptr->cvobf[ibf];
      const int coc00 = solver->cvocv_i[icv];

      const double D_coeff = (Sc_lam*solver->mu_lam[icv] + Sc_t*solver->mu_sgs[icv])*zone_ptr->area_over_delta_bf[ibf]; 
      //assume asgs_coeff==0, p_wall==p[icv]

      if (A!=NULL) A[coc00] += D_coeff;
      if (At!=NULL) At[coc00] += D_coeff;

      rhs[icv] += D_coeff*scbc_vals[isc];

    }
  }
  //else zero-flux bc

}

void SlipWallHBc::initData() {
  Param * param = getParam(getName());
  assert( param != NULL);
  parseScalarBc(param,false);
}

void SlipWallHBc::force_bf(double (*f_bf)[9]) const {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    FOR_I3 {
      f_bf[ibf][3*0+i] = 0.0;
      f_bf[ibf][3*1+i] = solver->p[icv0]*zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = 0.0;
    }
  }
  addPressureForceDueToFrameRotation(f_bf);
}

void WallHBc::initData() {

  assert( solver != NULL);
  assert(u_bc == NULL);  u_bc = new double[zone_ptr->nbf][3];

  double uwall[3] = {0.0, 0.0, 0.0};
  double xrot[3]  = {0.0, 0.0, 0.0};
  double omega[3] = {0.0, 0.0, 0.0};
  double urot     = 0.0 ;

  bool isRotating = false;
  bool isStationaryInRotatingFrame = false;
  Param * param = getParam(getName());
  if ( param->size() > 1 ){
    string token = param->getString(1);
    if ( token=="U_WALL"){
      uwall[0]  = param->getDouble(2);
      uwall[1]  = param->getDouble(3);
      uwall[2]  = param->getDouble(4);
      //COUT1(" > Found U_WALL " << u_bc[0] << ", " << u_bc[1] << ", " << u_bc[2]);
    }
    else if ( token=="ROTATING"){
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
    else if ( (token == "STATIONARY") || (token == "STATIONARY_FRAME")) {
      if ( solver->frame_rotation != NULL) {
        isStationaryInRotatingFrame = true;
      }
    }
  }

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    if (isRotating){
      double r[3];
      FOR_I3 r[i] = zone_ptr->x_bf[ibf][i] - xrot[i];
      double omega_cross_r[3] = CROSS_PRODUCT(omega,r);
      FOR_I3 u_bc[ibf][i] = urot * omega_cross_r[i];
    }
    else if (isStationaryInRotatingFrame){
      double r[3];
      FOR_I3 r[i] = zone_ptr->x_bf[ibf][i] - solver->frame_rotation[i+3];
      u_bc[ibf][0] = -solver->frame_rotation[1]*r[2] +
                      solver->frame_rotation[2]*r[1];

      u_bc[ibf][1] = -solver->frame_rotation[2]*r[0] +
                      solver->frame_rotation[0]*r[2];

      u_bc[ibf][2] = -solver->frame_rotation[0]*r[1] +
                      solver->frame_rotation[1]*r[0];
    }
    else{
      FOR_I3 u_bc[ibf][i] = uwall[i];
    }
  }

  parseScalarBc(param,false);
}

void WallHBc::addMomentumFlux(double * A,double (*rhs)[3]) const {


  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    const int coc00 = solver->cvocv_i[icv0];
    const double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])*zone_ptr->area_over_delta_bf[ibf];
    A[coc00] += mu_coeff;
    FOR_I3 rhs[icv0][i] += mu_coeff*u_bc[ibf][i];
  }

}

void WallHBc::addMomentumFlux(double * A,double * At,double (*rhs)[3]) const {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    const int coc00 = solver->cvocv_i[icv0];
    const double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])*zone_ptr->area_over_delta_bf[ibf];
    A[coc00] += mu_coeff;
    At[coc00] += mu_coeff;
    if (rhs != NULL)
      FOR_I3 rhs[icv0][i] += mu_coeff*u_bc[ibf][i];
  }

}

void WallHBc::force_bf(double (*f_bf)[9]) const {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])*zone_ptr->area_over_delta_bf[ibf];

    FOR_I3 {
      f_bf[ibf][3*0+i] = 0.0;
      f_bf[ibf][3*1+i] = solver->p[icv0]*zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = mu_coeff*(solver->u[icv0][i] - u_bc[ibf][i]);
    }
  }
  addPressureForceDueToFrameRotation(f_bf);
}

void WallHBc::query(const string& param_str) const {

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

CtiRegister::CtiDataError WallHBc::funcEvalCtiData(CtiRegister::CtiData& v,
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

            const int icv           = zone_ptr->cvobf[ibf];
            const double visc_coeff = (solver->mu_lam[icv] + solver->mu_sgs[icv])*
                                      zone_ptr->area_over_delta_bf[ibf];
            double du_mag = 0.0;
            for (int i =0; i < 3; ++i)
              du_mag  += (solver->u[icv][i] - u_bc[ibf][i])*(solver->u[icv][i] - u_bc[ibf][i]);

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
            const double visc_coeff = (solver->mu_lam[icv] + solver->mu_sgs[icv])*
              zone_ptr->area_over_delta_bf[ibf];

            double local_u[3] = DIFF(solver->u[icv],u_bc[ibf]);
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
          const double visc_coeff = (solver->mu_lam[icv] + solver->mu_sgs[icv])/y1;

          double du_mag = 0.0;
          for (int i =0; i < 3; ++i)
            du_mag  += (solver->u[icv][i] - u_bc[ibf][i])*(solver->u[icv][i] - u_bc[ibf][i]);

          const double tau_w      = visc_coeff*sqrt(du_mag);
          tmp[ibf]                = y1* sqrt(tau_w*solver->rho[icv])/solver->mu_lam[icv];

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
          blde = new BoundaryLayerDataExchanger<HelmholtzSolver>(zone_ptr,solver,nbl,lbl);

        blde->computeBlFromPt(tmp,pt_cut*pt_inf);


      }

      return CTI_DATA_OK;
    }
  }

  return HelmholtzBc::funcEvalCtiData(v,name,args,b_eval_func);
}


void InletHBc::initData() {

  assert(u_bc == NULL);  u_bc = new double[zone_ptr->nbf][3];
  assert(mf   == NULL);  mf   = new double[zone_ptr->nbf];

  Param * param = getParam(getName());
  const double rho_bc = getDoubleParam("RHO");

  double u_in[3];
  bool isStationaryInRotatingFrame = false;
  if (param->getString(1) == "MDOT") {
    double mf = param->getDouble(2);
    const double nmag = MAG(zone_ptr->n_global);
    assert(nmag > 0.0);
    FOR_I3 u_in[i] = -mf*zone_ptr->n_global[i]/(nmag*rho_bc); // negative for inward

    if (param->size()>3){
      if (param->getString(3)=="STATIONARY"&&solver->frame_rotation!=NULL)
        isStationaryInRotatingFrame = true;
    }
  }
  else {
    u_in[0] = param->getDouble(1);
    u_in[1] = param->getDouble(2);
    u_in[2] = param->getDouble(3);

    if (param->size()>4){
      if (param->getString(4)=="STATIONARY"&&solver->frame_rotation!=NULL)
        isStationaryInRotatingFrame = true;
    }
  }

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    FOR_I3 u_bc[ibf][i] = u_in[i];
    if (isStationaryInRotatingFrame){
      double r[3];
      FOR_I3 r[i] = zone_ptr->x_bf[ibf][i] - solver->frame_rotation[i+3];
      u_bc[ibf][0] += -solver->frame_rotation[1]*r[2] +
                       solver->frame_rotation[2]*r[1];

      u_bc[ibf][1] += -solver->frame_rotation[2]*r[0] +
                       solver->frame_rotation[0]*r[2];

      u_bc[ibf][2] += -solver->frame_rotation[0]*r[1] +
                       solver->frame_rotation[1]*r[0];
    }
    mf[ibf] = rho_bc*DOT_PRODUCT(zone_ptr->n_bf[ibf],u_bc[ibf]);
  }

  parseScalarBc(param,true);

}

void InletHBc::initFromCoarseGrid(StaticSolver::CoarseGrid* cg) {

  assert(mf == NULL); mf = new double[cg->ncb]; // TODO limit to zone
  for (int ibf = 0; ibf < cg->ncb; ++ibf)
    mf[ibf] = 0.0;

}

void InletHBc::setBc() {

  // this routine would be non-empty if there was time variation in the inlet bc...

}

void InletHBc::addMomentumFlux(double * A,double (*rhs)[3]) const {

  assert( solver != NULL);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    const int coc00 = solver->cvocv_i[icv0];
    const double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])* zone_ptr->area_over_delta_bf[ibf];
    if (mf[ibf] >= 0.0) {
      // outflow -- use first-order upwind...
      A[coc00] += mf[ibf] + mu_coeff;
      FOR_I3 rhs[icv0][i] += mu_coeff*u_bc[ibf][i];
    }
    else {
      // inflow...
      A[coc00] += mu_coeff;
      FOR_I3 rhs[icv0][i] += (-mf[ibf] + mu_coeff)*u_bc[ibf][i];
    }
  }

}

void InletHBc::addMomentumFlux(double *A,double *At,double (*rhs)[3]) const {

  assert( solver != NULL);

  if (rhs != NULL) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv0 = zone_ptr->cvobf[ibf];
      const int coc00 = solver->cvocv_i[icv0];
      const double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])* zone_ptr->area_over_delta_bf[ibf];
      if (mf[ibf] >= 0.0) {
        // outflow -- use first-order upwind...
        A[coc00] += mf[ibf] + mu_coeff;
        At[coc00] += mf[ibf] + mu_coeff;
        FOR_I3 rhs[icv0][i] += mu_coeff*u_bc[ibf][i];
      }
      else {
        // inflow...
        A[coc00] += mu_coeff;
        At[coc00] += mu_coeff;
        FOR_I3 rhs[icv0][i] += (-mf[ibf] + mu_coeff)*u_bc[ibf][i];
      }
    }
  }
  else {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv0 = zone_ptr->cvobf[ibf];
      const int coc00 = solver->cvocv_i[icv0];
      const double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])* zone_ptr->area_over_delta_bf[ibf];
      if (mf[ibf+zone_ptr->ibf_f] >= 0.0) {
        // outflow -- use first-order upwind...
        A[coc00] += mf[ibf+zone_ptr->ibf_f] + mu_coeff;
        At[coc00] += mf[ibf+zone_ptr->ibf_f] + mu_coeff;
      }
      else {
        // inflow...
        A[coc00] += mu_coeff;
        At[coc00] += mu_coeff;
      }
    }

  }

}

void InletHBc::addMassFlux(double * rhs) const {

  assert( solver != NULL);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    rhs[icv0] += mf[ibf];
  }
}

void InletHBc::force_bf(double (*f_bf)[9]) const {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv0 = zone_ptr->cvobf[ibf];
    double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])*zone_ptr->area_over_delta_bf[ibf];

    FOR_I3 {
      f_bf[ibf][3*0+i] = mf[ibf]*u_bc[ibf][i];
      f_bf[ibf][3*1+i] = solver->p[icv0]*zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = mu_coeff*(solver->u[icv0][i] - u_bc[ibf][i]);
    }
  }
  addPressureForceDueToFrameRotation(f_bf);
}

void InletHBc::query(const string& param_str) const {


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

void computeBfWallDistance(double * wall_dist, const BfZone * zone_ptr, const HelmholtzSolver* solver) {

  const int ncv = solver->ncv;
  double * cv_flag = new double[ncv];
  FOR_ICV cv_flag[icv] = -2.0;

  const int nbf = zone_ptr->nbf;
  // flag cells adjacent to faces of interest
  FOR_IBF {
    cv_flag[zone_ptr->cvobf[ibf]] = -1.0;
  }

  // We have flagged the cvs adjacent to this zone and we can loop over other boundary conditions to find cells at the corner of this zone and wall BCs. In only these cells we will compute the wall distance
  for(vector<HelmholtzBc*>::const_iterator it = solver->bcs.begin(); it != solver->bcs.end(); ++it) {
    string s=(*it)->getName();
    Param* p=getParam(s);
    if ((p->getString()=="WALL") || (p->getString()=="WM_ALG_WALL") || (p->getString()=="WM_SLIP_WALL")) {
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

void InletHBcProfile::initData() {
  assert(u_bc == NULL);  u_bc = new double[zone_ptr->nbf][3];
  assert(mf   == NULL);  mf   = new double[zone_ptr->nbf];

  Param * param = getParam(getName());
  assert( param != NULL);

  int type = -1;
  FluentBPReader profile;

  int iarg = 1;
  while (iarg < param->size()) {
    string token = param->getString(iarg++);
    if (token == "FILE") {
      profile.init(param->getString(iarg++));
    }
    else if (token == "TYPE") {
      const string _type = param->getString(iarg++);
      if (_type == "URHO") {
        type = 0;
      }
      else if (_type == "U") {
        type = 1;
      }
      else {
        CERR("unsupported INLET_PROFILE TYPE " << _type << "; please choose from URHO or U");
      }
    }
    else if (token == "FROM_MEAN") {
      profile.setUseMean(true);
    }
  }

  if (!profile.isInitialized()) {
    CERR("boundary profile was not properly initialized");
  }

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
    CERR("INLET_PROFILE \"TYPE\" must be specified");
  }
  else if (type == 0) {
    profile.ensureVar("x-velocity");
    profile.ensureVar("y-velocity");
    profile.ensureVar("z-velocity");
    profile.ensureVar("density");

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      u_bc[ibf][0] = profile.getData(ibf,"x-velocity");
      u_bc[ibf][1] = profile.getData(ibf,"y-velocity");
      u_bc[ibf][2] = profile.getData(ibf,"z-velocity");
      const double rho_ibf = profile.getData(ibf,"density");
      mf[ibf] = rho_ibf*DOT_PRODUCT(zone_ptr->n_bf[ibf],u_bc[ibf]);
    }
  }
  else if (type == 1) {
    profile.ensureVar("x-velocity");
    profile.ensureVar("y-velocity");
    profile.ensureVar("z-velocity");

    const double rho_ibf = getDoubleParam("RHO");

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      u_bc[ibf][0] = profile.getData(ibf,"x-velocity");
      u_bc[ibf][1] = profile.getData(ibf,"y-velocity");
      u_bc[ibf][2] = profile.getData(ibf,"z-velocity");
      mf[ibf] = rho_ibf*DOT_PRODUCT(zone_ptr->n_bf[ibf],u_bc[ibf]);
    }
  }
  else CERR("unrecognized type somehow set...");

  parseScalarBc(param,true);
}

void InletHBcProfile::initFromCoarseGrid(StaticSolver::CoarseGrid* cg) {

  //assert(u_bc == NULL);  u_bc = new double[zone_ptr->nbf][3];
  assert(mf   == NULL);  mf   = new double[cg->ncb];
  for (int ibf = 0; ibf < cg->ncb; ++ibf)
    mf[ibf]     = 0.0;

}

void InflowTurbulenceHBc::initData() {

  assert(u_bc == NULL);  u_bc = new double[zone_ptr->nbf][3];
  assert(mf   == NULL);  mf   = new double[zone_ptr->nbf];

  //IT Stuff
  assert(up_it == NULL); up_it   = new double[zone_ptr->nbf][3];
  assert(um_it == NULL); um_it   = new double[zone_ptr->nbf][3];
  assert(Rd_it == NULL); Rd_it   = new double[zone_ptr->nbf][3];
  assert(Rod_it == NULL); Rod_it   = new double[zone_ptr->nbf][3];
  assert(length_scale == NULL); length_scale  = new double[zone_ptr->nbf];
  assert(avg == NULL); avg   = new double[zone_ptr->nbf][3];
  assert(rms == NULL); rms   = new double[zone_ptr->nbf][3];

  IturbWgt=0.0;
  max_length=1.0E+200;
  in_plane_length=1.0E+200;
  length_ratio=1.0;
  Un=1.0E+200;
  inflow_zero=0.001;
  inflow_maxiter=30;
  //END IT Stuff

  //zero
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    FOR_I3{
      u_bc[ibf][i] = 0.0;
      up_it[ibf][i] = 0.0;
      um_it[ibf][i] = 0.0;
      Rd_it[ibf][i] = 0.0;
      Rod_it[ibf][i] = 0.0;
      avg[ibf][i] = 0.0;
      rms[ibf][i] = 0.0;
    }
    mf[ibf] = 0.0;
    length_scale[ibf] = 0.0;
  }

  Param * param = getParam(getName());
  assert( param != NULL);
  parseScalarBc(param,true);
}

void InflowTurbulenceHBc::initFromCoarseGrid(StaticSolver::CoarseGrid* cg) {
  assert(mf   == NULL);  mf   = new double[cg->ncb];
  for (int ibf = 0; ibf < cg->ncb; ++ibf)
    mf[ibf]     = 0.0;

}

void InflowTurbulenceHBc::initialHook() {

  //set mean flow profile
  // set the statistics of the turbulence...SHould this only happen at step 0 ? Or when the user asks?
  setStatistics();

  // we need the reduced connectivity of the boundary attached cvs. These are computed here
  buildITBoundaryStructures();

  //now we compute and store the lengthscale based on the distance to the wall and/or whatever other input the user has specified
  // same problem as set statistics
  setLengthscale();

  //create the operators
  buildFilterOperators();

  // This should only be done when we are at step 0....or if the user would like to reset the inflow turbulence...
  initiateStats();

  const double rho_bc = getDoubleParam("RHO");
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    if (b_modify_u){
      FOR_I3 u_bc[ibf][i] = um_it[ibf][i];
    }
    else{
      FOR_I3 u_bc[ibf][i] = up_it[ibf][i]+um_it[ibf][i];
    }
    mf[ibf] = rho_bc*DOT_PRODUCT(zone_ptr->n_bf[ibf],u_bc[ibf]);
  }

}

void InflowTurbulenceHBc::syncPostState(){
  const double rho_bc = getDoubleParam("RHO");
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    if (b_modify_u){
      FOR_I3 u_bc[ibf][i] = um_it[ibf][i];
    }
    else{
      FOR_I3 u_bc[ibf][i] = up_it[ibf][i]+um_it[ibf][i];
    }
    mf[ibf] = rho_bc*DOT_PRODUCT(zone_ptr->n_bf[ibf],u_bc[ibf]);
  }
}

void InflowTurbulenceHBc::resetBc(){
  //This should only be done when we are at step 0....or if the user would like to reset the inflow turbulence...
  initiateStats();
}

void InflowTurbulenceHBc::setBc() { 

  updateInflowTurbulence(up_it);
  const double rho_bc = getDoubleParam("RHO");
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    if (b_modify_u){
      FOR_I3 u_bc[ibf][i] = um_it[ibf][i];
    }
    else{
      FOR_I3 u_bc[ibf][i] = up_it[ibf][i]+um_it[ibf][i];
    }
    mf[ibf] = rho_bc*DOT_PRODUCT(zone_ptr->n_bf[ibf],u_bc[ibf]);
  }

  if (!b_modify_u&&solver->step%solver->check_interval==0){
    if (mpi_rank==0)
      cout << getName() << " INFLOW_TURB setBc()" << endl;

    dumpRange(u_bc,zone_ptr->nbf,getName()+":u_bc");
  }

}

void InflowTurbulenceHBc::modifyU(){
  if (b_modify_u){
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv0 = zone_ptr->cvobf[ibf];  
      FOR_I3 solver->u[icv0][i] = up_it[ibf][i]+um_it[ibf][i];
    }
    if (solver->step%solver->check_interval==0){
      if (mpi_rank==0)
        cout << getName() << " INFLOW_TURB modifyU()" << endl;

      dumpRange(um_it,zone_ptr->nbf,getName()+":um_it");
      dumpRange(up_it,zone_ptr->nbf,getName()+":up_it");
    }
  }
}

void OutletHBc::initData() {

  Param * param = getParam(getName());
  sigma         = param->getDouble(1);
  L             = param->getDouble(2);
  p_ref         = param->getDouble(3);
  Ma            = param->getDouble(4);

  assert(p_bc   == NULL);  p_bc   = new double[zone_ptr->nbf];
  assert(p0_bc  == NULL);  p0_bc  = new double[zone_ptr->nbf];
  assert(p00_bc == NULL);  p00_bc = new double[zone_ptr->nbf];
  assert(mf == NULL);      mf     = new double[zone_ptr->nbf];

  //registered data will be overwritten upon reading
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf){
    mf[ibf]     = 0.0;
    p_bc[ibf]   = 0.0;
    p0_bc[ibf]  = 0.0;
    p00_bc[ibf] = 0.0;
  }

  if (!b_helmholtzCorrection){
    assert(mf0 == NULL);     mf0    = new double[zone_ptr->nbf];
    assert(mf00 == NULL);    mf00   = new double[zone_ptr->nbf];
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf){
      mf0[ibf]    = 0.0;
      mf00[ibf]   = 0.0;
    }
  }


  // if reverse flow is encountered at the outlet, the passive scalars
  // will be reintroduced with 0 value as a default.

  assert( scbc_vals == NULL); scbc_vals = new double[solver->nsc];
  for (int isc = 0; isc < solver->nsc; ++isc)
    scbc_vals[isc] = 0.0;

}

void OutletHBc::initFromCoarseGrid(StaticSolver::CoarseGrid* cg) {

//  Param * param = getParam(getName());
//  sigma         = param->getDouble(1);
//  L             = param->getDouble(2);
//  p_ref         = param->getDouble(3);
//  Ma            = param->getDouble(4);
//
  assert(mf == NULL); mf = new double[cg->ncb];

  for (int ibf = 0; ibf < cg->ncb; ++ibf)
    mf[ibf] = 0.0;

}


void OutletHBc::resetBc() {

  assert( zone_ptr);

  if ( mpi_rank == 0)
    cout << " > resetting the outlet bc ... " ;

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    mf[ibf]    = 0.0;
    p_bc[ibf]    = 0.0;
    p0_bc[ibf]   = 0.0;
    p00_bc[ibf]  = 0.0;
  }

  if (!b_helmholtzCorrection){
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      mf0[ibf]     = 0.0;
      mf00[ibf]    = 0.0;
    }
  }


  if ( mpi_rank == 0)
    cout << " ok. " << endl;

}

void OutletHBc::initialHook() {
  // mf is now registered data
  if (!zone_ptr->checkDataFlag("mf"))
    updateMf();
  
  if (!b_helmholtzCorrection){
    if (!zone_ptr->checkDataFlag("mf0"))
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
        mf0[ibf] = mf[ibf];

    if (!zone_ptr->checkDataFlag("mf00"))
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
        mf00[ibf] = mf0[ibf];
  }

  if (!zone_ptr->checkDataFlag("p_bc"))
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
      p_bc[ibf] = solver->p[zone_ptr->cvobf[ibf]];

  if (!zone_ptr->checkDataFlag("p0_bc"))
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
      p0_bc[ibf] = p_bc[ibf];

  if (!zone_ptr->checkDataFlag("p00_bc"))
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
      p00_bc[ibf] = p0_bc[ibf];
}

void OutletHBc::setBc() {
  //update boundary pressure and mf
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    //tmp  = 3.0*p_bc[ibf] - 3.0*p0_bc[ibf] + p00_bc[ibf];
    p00_bc[ibf] = p0_bc[ibf];
    p0_bc[ibf]  = p_bc[ibf];
    //p_bc[ibf]   = tmp;
  }

  if (!b_helmholtzCorrection){
    if( solver->mf_init>0 
     || checkParam("NO_EXTRAPOLATE_MF") 
     || solver->n_outer_iters!=1) {
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
        mf00[ibf] = mf0[ibf];
        mf0[ibf]  = mf[ibf];
      }
    }
    else {
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
        double tmp  = 3.0*mf[ibf] - 3.0*mf0[ibf] + mf00[ibf];
        mf00[ibf] = mf0[ibf];
        mf0[ibf]  = mf[ibf];
        mf[ibf]   = tmp;
      }
    }
  }

}

void OutletHBc::updateMf() {

  // compute mf

  double my_buf    = 0.0;
  double rhou_conv = 0.0;

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    my_buf += solver->rho[icv0] * DOT_PRODUCT(solver->u[icv0], zone_ptr->n_bf[ibf]) ;
  }

  MPI_Allreduce(&my_buf,&rhou_conv,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
  rhou_conv = max(0.0, rhou_conv/zone_ptr->area_global) ;

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    mf[ibf] = rhou_conv*zone_ptr->area_bf[ibf];
    //const int icv0 = zone_ptr->cvobf[ibf];
    //mf[ibf]   = max(0.0,solver->rho[icv0]*DOT_PRODUCT(solver->u[icv0],zone_ptr->n_bf[ibf]));
  }
}

void OutletHBc::updateP() {

  // compute the outlet pressure ...

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv0  = zone_ptr->cvobf[ibf];
    const double ds = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];

    //const double pp = solver->p[icv0];
    double pp = 0.0;
    for (int coc = solver->cvocv_i[icv0]+1; coc != solver->cvocv_i[icv0+1]; ++coc) {
      const int icv_nbr =  solver->cvocv_v[coc];
      pp += solver->p[icv_nbr];
    }

    pp /= double(solver->cvocv_i[icv0+1]-solver->cvocv_i[icv0]-1);
    pp  = 0.5*pp + 0.5*solver->p[icv0];

    p_bc[ibf] = ( (2.0*p0_bc[ibf] - 0.5*p00_bc[ibf])/(solver->helmholtz_sos*(1.0+Ma)*solver->dt) + pp/ds + p_ref*sigma/L ) /
                  ( 1.5/(solver->helmholtz_sos*(1.0+Ma)*solver->dt) + 1.0/ds + sigma/L ) ;

    //p_bc[ibf] = 0.0;
  }

  //MiscUtils::dumpRange(p_bc,zone_ptr->nbf,"p_bc (outlet)");

}

void OutletHBc::addMomentumFlux(double * A,double (*rhs)[3]) const {

  assert( solver != NULL);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    const int coc00 = solver->cvocv_i[icv0];

    // first order upwind the outlet
    assert(mf[ibf] >= 0.0); // no backflow...
    A[coc00] += mf[ibf];
  }

}

void OutletHBc::addMomentumFlux(double * A,double *At,double (*rhs)[3]) const {

  assert( solver != NULL);

  if (rhs != NULL) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv0 = zone_ptr->cvobf[ibf];
      const int coc00 = solver->cvocv_i[icv0];

      // first order upwind the outlet
      assert(mf[ibf] >= 0.0); // no backflow...
      A[coc00] += mf[ibf];
      At[coc00] += mf[ibf];
    }
  }
  else {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv0 = zone_ptr->cvobf[ibf];
      const int coc00 = solver->cvocv_i[icv0];

      // first order upwind the outlet
      assert(mf[ibf+zone_ptr->ibf_f] >= 0.0); // no backflow...
      A[coc00] += mf[ibf+zone_ptr->ibf_f];
      At[coc00] += mf[ibf+zone_ptr->ibf_f];
    }
  }

}

void OutletHBc::addPressureFlux(double *A,double *rhs) const {
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0  = zone_ptr->cvobf[ibf];
    const int coc00 = solver->cvocv_i[icv0];
    A[coc00] -= zone_ptr->area_over_delta_bf[ibf];
    if (rhs != NULL)
      rhs[icv0] -= p_bc[ibf]*zone_ptr->area_over_delta_bf[ibf];
  }
}

void OutletHBc::correctMf() {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    mf[ibf] -= (2.0/3.0)*solver->dt*zone_ptr->area_over_delta_bf[ibf]*(p_bc[ibf]-solver->p[icv0]);
  }

}

void OutletHBc::addMassFlux(double * rhs) const {

  assert( solver != NULL);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    rhs[icv0] += mf[ibf];
  }

}

void OutletHBc::completePressureGrad(double (*dpdx)[3]) const {

  assert( solver != NULL);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    FOR_I3 dpdx[icv0][i] += zone_ptr->n_bf[ibf][i]*p_bc[ibf];
  }

}

void OutletHBc::force_bf(double (*f_bf)[9]) const {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv0 = zone_ptr->cvobf[ibf];
    FOR_I3 {
      f_bf[ibf][3*0+i] = mf[ibf]*solver->u[icv0][i];
      f_bf[ibf][3*1+i] = p_bc[ibf]*zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = 0.0;
    }
  }

}

void OutletHBc::query(const string& param_str) const {

  // report rhou_conv
  double my_rhou_conv = 0.0;
  double rhou_conv    = 0.0;

  // report min/max pressure
  double my_buf[2] = {HUGE_VAL,HUGE_VAL};
  double buf[2];

  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    my_buf[0] = min(my_buf[0],p_bc[ibf]);
    my_buf[1] = min(my_buf[1],-p_bc[ibf]);

    const int icv0 = zone_ptr->cvobf[ibf];
    my_rhou_conv += solver->rho[icv0] * DOT_PRODUCT(solver->u[icv0], zone_ptr->n_bf[ibf]) ;
  }

  MPI_Reduce(&my_rhou_conv,&rhou_conv,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  rhou_conv = max(0.0, rhou_conv/zone_ptr->area_global) ;

  MPI_Reduce(my_buf,buf,2,MPI_DOUBLE,MPI_MIN,0,mpi_comm);

  if ( mpi_rank == 0 ) {
    *ss << " > " << std::left << std::setw(30) << zone_ptr->getName() << " :  time, rhou_conv, dumpRange(p_bc) = " << std::right
        << std::setw(8) << solver->time << " " << std::setw(12) << rhou_conv << " "
        << std::setw(12) << buf[0] << " " << std::setw(12) << -buf[1] << endl;
  }
  flush();

}

CtiRegister::CtiDataError OutletHBc::funcEvalCtiData(CtiRegister::CtiData& v,
    const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if ( zone_ptr->isLeadingStrMatched(name) ) {

    //p_bf function in HelmholtzBc is overriden here to use p_bc instead of p[icv]
    const string pbf_str = zone_ptr->getName() + ":" + "p_bf";

    if ( name == pbf_str) {

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      if ( b_eval_func) {
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          tmp[ibf] = p_bc[ibf];
        }
      }

      return CTI_DATA_OK;

    }
  }

  return HelmholtzBc::funcEvalCtiData(v,name,args,b_eval_func);
}

void OutletVVHBc::initData() {

  Param * param = getParam(getName());
  sigma         = param->getDouble(1);
  L             = param->getDouble(2);
  p_ref         = param->getDouble(3);
  Ma            = param->getDouble(4);

  if (param->size()>5){
    string token = param->getString(5);
    if (token == "LOCAL_U"){
      b_localU = true;
    }
    else if (token == "X_OMEGA"){
      assert(!xr&&!nr);
      xr = new double[3];
      nr = new double[3];
      FOR_I3 xr[i] = param->getDouble(6+i);
      FOR_I3 nr[i] = param->getDouble(6+3+i);
      double mag_nr = MAG(nr);
      if (mag_nr<=0.0){
        CERR("OUTLET_VV ... X_OMEGA rotation axis must have non-zero magnitude");
      }
      FOR_I3 nr[i] /= mag_nr;
    }
    else {
      CERR("Invalid syntax for OUTLET_VV, expecting \n" <<
           "  OUTLET_VV <sigma> <L> <p_ref> <Ma> LOCAL_U\n" <<
           "  OUTLET_VV <sigma> <L> <p_ref> <Ma> X_OMEGA <x> <y> <z> <omega-x> <omega-y> <omega-z>");
    }
  }

  assert(p_bc   == NULL);   p_bc   = new double[zone_ptr->nbf];
  assert(p0_bc   == NULL);  p0_bc  = new double[zone_ptr->nbf];
  assert(p00_bc   == NULL); p00_bc = new double[zone_ptr->nbf];
  assert(mf == NULL);      mf     = new double[zone_ptr->nbf];

  //registered data will be overwritten upon reading
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    p_bc[ibf]   = 0.0;
    p0_bc[ibf]  = 0.0;
    p00_bc[ibf] = 0.0;
    mf[ibf]     = 0.0;
  }

  if (!b_helmholtzCorrection){
    assert(mf0 == NULL);     mf0    = new double[zone_ptr->nbf];
    assert(mf00 == NULL);    mf00   = new double[zone_ptr->nbf];
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf){
      mf0[ibf]    = 0.0;
      mf00[ibf]   = 0.0;
    }
  }

  // if reverse flow is encountered at the outlet, the passive scalars
  // will be reintroduced with 0 value as a default.

  assert( scbc_vals == NULL); scbc_vals = new double[solver->nsc];
  for (int isc = 0; isc < solver->nsc; ++isc)
    scbc_vals[isc] = 0.0;

}

void OutletVVHBc::initFromCoarseGrid(StaticSolver::CoarseGrid* cg) {

  Param * param = getParam(getName());
  sigma         = param->getDouble(1);
  L             = param->getDouble(2);
  p_ref         = param->getDouble(3);
  Ma            = param->getDouble(4);

 if (param->size()>5){
    string token = param->getString(5);
    if (token == "LOCAL_U"){
      b_localU = true;
    }
    else if (token == "X_OMEGA"){
      assert(!xr&&!nr);
      xr = new double[3];
      nr = new double[3];
      FOR_I3 xr[i] = param->getDouble(6+i);
      FOR_I3 nr[i] = param->getDouble(6+3+i);
      double mag_nr = MAG(nr);
      if (mag_nr<=0.0){
        CERR("OUTLET_VV ... X_OMEGA rotation axis must have non-zero magnitude");
      }
      FOR_I3 nr[i] /= mag_nr;
    }
    else {
      CERR("Invalid syntax for OUTLET_VV, expecting \n" <<
           "  OUTLET_VV <sigma> <L> <p_ref> <Ma> LOCAL_U\n" <<
           "  OUTLET_VV <sigma> <L> <p_ref> <Ma> X_OMEGA <x> <y> <z> <omega-x> <omega-y> <omega-z>");
    }
  }

  assert(mf == NULL); mf = new double[cg->ncb];

  for (int ibf = 0; ibf < cg->ncb; ++ibf)
    mf[ibf] = -1.0; //0.0; // for checking

}

void OutletVVHBc::resetBc() {

  assert( zone_ptr);

  if ( mpi_rank == 0)
    cout << " > resetting the outlet bc ... " ;

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    mf[ibf]      = 0.0;
    p_bc[ibf]    = 0.0;
    p0_bc[ibf]   = 0.0;
    p00_bc[ibf]  = 0.0;
  }

  if (!b_helmholtzCorrection){
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      mf0[ibf]     = 0.0;
      mf00[ibf]    = 0.0;
    }
  }

  if ( mpi_rank == 0)
    cout << " ok. " << endl;

}

void OutletVVHBc::initialHook() {
  // mf is now registered data
  if (!zone_ptr->checkDataFlag("mf"))
    updateMf();
  
  if (!b_helmholtzCorrection){
    if (!zone_ptr->checkDataFlag("mf0"))
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
        mf0[ibf] = mf[ibf];

    if (!zone_ptr->checkDataFlag("mf00"))
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
        mf00[ibf] = mf0[ibf];
  }

  if (!zone_ptr->checkDataFlag("p_bc"))
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
      p_bc[ibf] = solver->p[zone_ptr->cvobf[ibf]];

  if (!zone_ptr->checkDataFlag("p0_bc"))
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
      p0_bc[ibf] = p_bc[ibf];

  if (!zone_ptr->checkDataFlag("p00_bc"))
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
      p00_bc[ibf] = p0_bc[ibf];
}

void OutletVVHBc::setBc() {

  // update mass flux, copy pressure history
  // p_bc is not used until after updateP 
  // (not in add pressure flux, only complete pressure grad)
  // p0_bc, p00_bc store p[icv0] history
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf]; 
    p00_bc[ibf] = p0_bc[ibf];
    p0_bc[ibf]  = solver->p[icv0];
  }

  if (!b_helmholtzCorrection){
    if (solver->mf_init>0 
      || checkParam("NO_EXTRAPOLATE_MF") 
      || solver->n_outer_iters!=1) {
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
        mf00[ibf] = mf0[ibf];
        mf0[ibf]  = mf[ibf];
      }
    }
    else {
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
        double tmp  = 3.0*mf[ibf] - 3.0*mf0[ibf] + mf00[ibf];
        mf00[ibf] = mf0[ibf];
        mf0[ibf]  = mf[ibf];
        mf[ibf]   = max(0.0,tmp);
      }
    }
  }
}

void OutletVVHBc::updateMf() {

  // compute mf

  if (!b_localU&&!xr){

    double my_buf    = 0.0;
    double rhou_conv = 0.0;

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv0 = zone_ptr->cvobf[ibf];
      my_buf += solver->rho[icv0] * DOT_PRODUCT(solver->u[icv0], zone_ptr->n_bf[ibf]) ;
    }

    MPI_Allreduce(&my_buf,&rhou_conv,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
    rhou_conv = max(0.0, rhou_conv/zone_ptr->area_global) ;

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
      mf[ibf] = rhou_conv*zone_ptr->area_bf[ibf];

  }
  else if (xr){
    assert(nr);
    double my_buf[2] = {0.0,0.0};

    double buf[2]    = {0.0,0.0};

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv0 = zone_ptr->cvobf[ibf];
      my_buf[0] += solver->rho[icv0] * DOT_PRODUCT(solver->u[icv0], zone_ptr->n_bf[ibf]);
      double r[3] = DIFF(zone_ptr->x_bf[ibf],xr);
      double omega_x_r[3] = CROSS_PRODUCT(nr,r);
      my_buf[1] += DOT_PRODUCT(omega_x_r,zone_ptr->n_bf[ibf]);
    }

    MPI_Allreduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,mpi_comm);
    double rhou_conv = max(0.0, buf[0]/buf[1]) ;

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf){
      double r[3] = DIFF(zone_ptr->x_bf[ibf],xr);
      double omega_x_r[3] = CROSS_PRODUCT(nr,r);
      mf[ibf] = max(0.0,rhou_conv*DOT_PRODUCT(omega_x_r,zone_ptr->n_bf[ibf]));
    }
  }
  else{

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf){
      const int icv0 = zone_ptr->cvobf[ibf];

      double mf_local_avg = 0.0;
      for (int coc = solver->cvocv_i[icv0]+1; coc != solver->cvocv_i[icv0+1]; ++coc) {
        const int icv_nbr = solver->cvocv_v[coc];
        mf_local_avg += solver->rho[icv_nbr] * DOT_PRODUCT(solver->u[icv_nbr], zone_ptr->n_bf[ibf]);
      }

      mf_local_avg /= double(solver->cvocv_i[icv0+1]-solver->cvocv_i[icv0]-1);
      mf_local_avg  = 0.5*mf_local_avg + 0.5*solver->rho[icv0]*DOT_PRODUCT(solver->u[icv0], zone_ptr->n_bf[ibf]);

      mf[ibf] = max(mf_local_avg,0.0);
    }

  }

/* DAP: this was used in the vida vor outlet, originally called before correctContinuity 
 *      and solveU in setBc when two helmholtz solves were used. this was also used after
 *      setBc in the initial run of the single wheel which produced good results...
 *      might still be needed?
 *      moved to setBc...kind of like position before correctContinuity removal
 *      but not computed where single wheel stats were computed.
 
  // update boundary pressure

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    // use an ab3 extrapolator for setting the outlet pressure? ...
    // p_bc is not used until we update it, so the extrapolation shouldn't matter.

    double tmp  = 3.0*p_bc[ibf] - 3.0*p0_bc[ibf] + p00_bc[ibf];
    p00_bc[ibf] = p0_bc[ibf];
    p0_bc[ibf]  = p_bc[ibf];
    p_bc[ibf]   = tmp;

  }
*/

  if (solver->step%solver->check_interval==0){
    dumpRange(mf,zone_ptr->nbf,"set outlet mf");
    double * tmp = new double[zone_ptr->nbf];
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      tmp[ibf] = mf[ibf]/zone_ptr->area_bf[ibf];
    }
    dumpRange(tmp,zone_ptr->nbf,"set outlet mf/area");
    delete[] tmp;
  }


}

void OutletVVHBc::updateP() {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    const double ds = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
    //p_bc[ibf] = ( (2.0*p0_bc[ibf] - 0.5*p00_bc[ibf])/(solver->helmholtz_sos*(1.0+Ma)*solver->dt) + solver->p[icv0]/ds + p_ref*sigma/L ) /
    //            ( 1.5/(solver->helmholtz_sos*(1.0+Ma)*solver->dt) + 1.0/ds + sigma/L ) ;
    //DP: Replace p_bc with p_cv in relaxation term (sigma/L)
    p_bc[ibf] = ( (2.0*p0_bc[ibf] - 0.5*p00_bc[ibf])/(solver->helmholtz_sos*(1.0+Ma)*solver->dt) + solver->p[icv0]/ds + (p_ref-solver->p[icv0])*sigma/L ) / ( 1.5/(solver->helmholtz_sos*(1.0+Ma)*solver->dt) + 1.0/ds ) ;

  }

}

void OutletVVHBc::addMomentumFlux(double * A,double (*rhs)[3]) const {

  assert( solver != NULL);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    const int coc00 = solver->cvocv_i[icv0];

    // first order upwind the outlet
    assert(mf[ibf] >= 0.0) ; // no backflow...
    A[coc00] += mf[ibf];
  }

}

void OutletVVHBc::addMomentumFlux(double * A,double *At,double (*rhs)[3]) const {

  assert( solver != NULL);

  if (rhs != NULL) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv0 = zone_ptr->cvobf[ibf];
      const int coc00 = solver->cvocv_i[icv0];

      // first order upwind the outlet
      assert(mf[ibf] >= 0.0) ; // no backflow...
      A[coc00] += mf[ibf];
      At[coc00] += mf[ibf];
    }
  }
  else {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv0 = zone_ptr->cvobf[ibf];
      const int coc00 = solver->cvocv_i[icv0];

      // first order upwind the outlet
      assert(mf[ibf+zone_ptr->ibf_f] >= 0.0); // no backflow...
      A[coc00] += mf[ibf+zone_ptr->ibf_f];
      At[coc00] += mf[ibf+zone_ptr->ibf_f];
    }
  }

}

void OutletVVHBc::addPressureFlux(double *A,double *rhs) const {
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0  = zone_ptr->cvobf[ibf];
    const int coc00 = solver->cvocv_i[icv0];
    A[coc00]       -= zone_ptr->area_bf[ibf]*1.5/(solver->helmholtz_sos*(1.0+Ma)*solver->dt) ;
    A[coc00]       -= zone_ptr->area_bf[ibf]*sigma/L ;

    if (rhs != NULL){
      rhs[icv0] -= zone_ptr->area_bf[ibf]/(solver->helmholtz_sos*(1.0+Ma)*solver->dt)*(2.0*p0_bc[ibf] - 0.5*p00_bc[ibf]);
      rhs[icv0] -= zone_ptr->area_bf[ibf]*sigma/L*p_ref;
    }
  }
}

void OutletVVHBc::correctMf() {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    double dpcvdt = (1.5*solver->p[icv0] - 2.0*p0_bc[ibf] + 0.5*p00_bc[ibf])/solver->dt;
    double dpbfdx = -dpcvdt/(solver->helmholtz_sos*(1.0+Ma)) - (solver->p[icv0]-p_ref)*sigma/L;

    mf[ibf] = max(0.0,mf[ibf]-(2.0/3.0)*solver->dt*zone_ptr->area_bf[ibf]*dpbfdx);
 
  }

  if (solver->step%solver->check_interval==0){
    dumpRange(mf,zone_ptr->nbf,"corrected outlet mf");
    double * tmp = new double[zone_ptr->nbf];
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      tmp[ibf] = mf[ibf]/zone_ptr->area_bf[ibf];
    }
    dumpRange(tmp,zone_ptr->nbf,"corrected outlet mf/area");
    delete[] tmp;
  }


}

void OutletVVHBc::addMassFlux(double * rhs) const {

  assert( solver != NULL);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    rhs[icv0] += mf[ibf];
  }

}

void OutletVVHBc::completePressureGrad(double (*dpdx)[3]) const {

  assert( solver != NULL);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    FOR_I3 dpdx[icv0][i] += zone_ptr->n_bf[ibf][i]*p_bc[ibf];
  }

}

void OutletVVHBc::force_bf(double (*f_bf)[9]) const {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    FOR_I3 {
      f_bf[ibf][3*0+i] = mf[ibf]*solver->u[icv0][i];
      f_bf[ibf][3*1+i] = p_bc[ibf]*zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = 0.0;
    }
  }

}

void OutletVVHBc::query(const string& param_str) const {

  // report rhou_conv
  double my_rhou_conv = 0.0;
  double rhou_conv    = 0.0;

  // report min/max pressure
  double my_buf[2] = {HUGE_VAL,HUGE_VAL};
  double buf[2];

  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    my_buf[0] = min(my_buf[0],p_bc[ibf]);
    my_buf[1] = min(my_buf[1],-p_bc[ibf]);

    const int icv0 = zone_ptr->cvobf[ibf];
    my_rhou_conv += solver->rho[icv0] * DOT_PRODUCT(solver->u[icv0], zone_ptr->n_bf[ibf]) ;
  }

  MPI_Reduce(&my_rhou_conv,&rhou_conv,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  rhou_conv = max(0.0, rhou_conv/zone_ptr->area_global) ;

  MPI_Reduce(my_buf,buf,2,MPI_DOUBLE,MPI_MIN,0,mpi_comm);

  if ( mpi_rank == 0 ) {
    *ss << " > " << std::left << std::setw(30) << zone_ptr->getName() << " :  time, rhou_conv, dumpRange(p_bc) = " << std::right
        << std::setw(8) << solver->time << " " << std::setw(12) << rhou_conv << " "
        << std::setw(12) << buf[0] << " " << std::setw(12) << -buf[1] << endl;
  }
  flush();

}

CtiRegister::CtiDataError OutletVVHBc::funcEvalCtiData(CtiRegister::CtiData& v,
    const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if ( zone_ptr->isLeadingStrMatched(name) ) {

    //p_bf function in HelmholtzBc is overriden here to use p_bc instead of p[icv]
    const string pbf_str = zone_ptr->getName() + ":" + "p_bf";

    if ( name == pbf_str) {

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      if ( b_eval_func) {
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          tmp[ibf] = p_bc[ibf];
        }
      }

      return CTI_DATA_OK;

    }
  }

  return HelmholtzBc::funcEvalCtiData(v,name,args,b_eval_func);
}



SpongeHBc::SpongeHBc(BfZone* p, HelmholtzSolver* s) : OutletVVHBc(p,s) {

  FOR_I5 sponge_data[i] = 0.0;
  sponge_strength = -1.0;
  sponge_speed = -1.0;
  sponge_type = SPONGE_UNDEFINED;
  p_ref = -1.0;
  sigma = -1.0;
  FOR_I3 th_e0[i] = 0.0;
  FOR_I3 th_e1[i] = 0.0;

}

SpongeHBc::SpongeHBc(BfZone* p, HelmholtzSolver* s,const int icg) : OutletVVHBc(p,s,icg) {

  FOR_I5 sponge_data[i] = 0.0;
  sponge_strength = -1.0;
  sponge_speed = -1.0;
  sponge_type = SPONGE_UNDEFINED;
  p_ref = -1.0;
  sigma = -1.0;
  FOR_I3 th_e0[i] = 0.0;
  FOR_I3 th_e1[i] = 0.0;
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

void SpongeHBc::initData() {

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

  assert(p_bc   == NULL);   p_bc   = new double[zone_ptr->nbf];
  assert(p0_bc   == NULL);  p0_bc  = new double[zone_ptr->nbf];
  assert(p00_bc   == NULL); p00_bc = new double[zone_ptr->nbf];
  assert(mf == NULL);       mf     = new double[zone_ptr->nbf];
  
  // registered data will be overwritten upon reading
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    p_bc[ibf]   = 0.0;
    p0_bc[ibf]  = 0.0;
    p00_bc[ibf] = 0.0;
    mf[ibf]     = 0.0;
  }

  if (!b_helmholtzCorrection){
    assert(mf0 == NULL);      mf0    = new double[zone_ptr->nbf];
    assert(mf00 == NULL);     mf00   = new double[zone_ptr->nbf];
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      mf0[ibf]    = 0.0;
      mf00[ibf]   = 0.0;
    }
  }


}

void SpongeHBc::initFromCoarseGrid(StaticSolver::CoarseGrid* cg) {

  Param* param = getParam(getName());

  // default params...

  parseSpongeParamHBc(param,sponge_type,L,sponge_speed,sponge_strength,p_ref,sigma,xr,nr);
  initSpongeTypeTheta();

  // use faces to get the domain limits (no nodes on coarse grids)...

  double my_buf[5] = { -1.0E-20, -1.0E-20, -1.0E-20, -1.0E-20, -1.0E-20 };
  for (int ifa = 0; ifa < solver->nfa; ++ifa) {
    my_buf[0] = max(my_buf[0],solver->x_fa[ifa][0]);
    my_buf[1] = max(my_buf[1],solver->x_fa[ifa][1]);
    my_buf[2] = max(my_buf[2],solver->x_fa[ifa][2]);
    my_buf[3] = max(my_buf[3],sqrt(solver->x_fa[ifa][1]*solver->x_fa[ifa][1] + solver->x_fa[ifa][2]*solver->x_fa[ifa][2]));
    my_buf[4] = max(my_buf[4],sqrt(solver->x_fa[ifa][0]*solver->x_fa[ifa][0] + solver->x_fa[ifa][1]*solver->x_fa[ifa][1]));
  }
  MPI_Allreduce(my_buf,sponge_data,5,MPI_DOUBLE,MPI_MAX,mpi_comm);

  assert(mf == NULL); mf = new double[cg->ncb];
  for (int ibf = 0; ibf < cg->ncb; ++ibf)
    mf[ibf] = -1.0; // for checking

}

void SpongeHBc::initSpongeTypeTheta(){
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

void SpongeHBc::addMomentumFlux(double * A, double (*rhs)[3]) const {

  assert(solver != NULL);

  OutletVVHBc::addMomentumFlux(A,rhs);

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
      break;
    }
    case SPONGE_LTH:
    {
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
          A[coc00] += solver->vol_cv[icv]*sponge_coeff*solver->rho[icv];
          for (int i =0; i < 3; ++i)
            rhs[icv][i] += solver->vol_cv[icv]*sponge_coeff*solver->rho[icv]*sponge_velocity[i];
        }
      }
      break;
    }
    default:
      assert(0);
  }

}

void SpongeHBc::addMomentumFlux(double * A,double* At,double (*rhs)[3]) const {

  assert( solver != NULL);

  OutletVVHBc::addMomentumFlux(A,At,rhs);

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
          At[coc00] += solver->vol_cv[icv]*sponge_coeff*solver->rho[icv];
          if (rhs != NULL)
            FOR_I3 rhs[icv][i] += solver->vol_cv[icv]*sponge_coeff*solver->rho[icv]*sponge_velocity[i];
        }
      }

      break;
    }
    case SPONGE_LTH:
    {
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

          const int coc00 = solver->cvocv_i[icv];
          A[coc00] += solver->vol_cv[icv]*sponge_coeff*solver->rho[icv];
          At[coc00] += solver->vol_cv[icv]*sponge_coeff*solver->rho[icv];
          if (rhs!=NULL){
            double sponge_velocity[3] = CROSS_PRODUCT(nr,r_cv);
            FOR_I3 sponge_velocity[i] *=sponge_speed*2.0*M_PI;
            for (int i =0; i < 3; ++i)
              rhs[icv][i] += solver->vol_cv[icv]*sponge_coeff*solver->rho[icv]*sponge_velocity[i];
          }
        }
      }
      break;
    }
    default:
      assert(0);
  }

}

void SpongeHBc::addPressureFlux(double *A,double *rhs) const {
  OutletVVHBc::addPressureFlux(A,rhs);

  switch(sponge_type) {
    case SPONGE_LX:
    case SPONGE_LY:
    case SPONGE_LZ:
    {
      int idir;
      if (sponge_type == SPONGE_LX) {
        idir = 0;
      } else if (sponge_type == SPONGE_LY) {
        idir = 1;
      } else {
        assert(sponge_type == SPONGE_LZ);
        idir = 2;
      }

      const double Lmax = sponge_data[idir];
      for (int icv = 0; icv < solver->ncv; ++icv) {
        if (solver->x_cv[icv][idir] > Lmax-L) {
          const double x_sp = (solver->x_cv[icv][idir] - (Lmax-L))/L;
          // sponge profile = a*x^2 + b*x^8..
          // const double sponge_coeff = sigma*(0.068*x_sp*x_sp + 0.845*pow(x_sp,8));
          // quadratic sponge profile (see JCP 212 681-702 (2006) & JCP 231 704-716 (2012))
          const double sponge_coeff = sponge_strength*x_sp*x_sp;
          const int coc00 = solver->cvocv_i[icv];
          A[coc00] -= 1.5*sponge_coeff*solver->vol_cv[icv]/(solver->helmholtz_sos*solver->helmholtz_sos*solver->dt);
        }
      }

      break;
    }
    case SPONGE_LTH:
    {
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
          const int coc00 = solver->cvocv_i[icv];
          A[coc00] -= 1.5*sponge_coeff*solver->vol_cv[icv]/(solver->helmholtz_sos*solver->helmholtz_sos*solver->dt);
        }
      }
      break;
    }


    default:
      assert(0);
  }

}

void SpongeHBc::addContSource(double *cont) const {

  switch(sponge_type) {
    case SPONGE_LX:
    case SPONGE_LY:
    case SPONGE_LZ:
    {
      int idir;
      if (sponge_type == SPONGE_LX) {
        idir = 0;
      } else if (sponge_type == SPONGE_LY) {
        idir = 1;
      } else {
        assert(sponge_type == SPONGE_LZ);
        idir = 2;
      }

      const double Lmax = sponge_data[idir];
      for (int icv = 0; icv < solver->ncv; ++icv) {
        if (solver->x_cv[icv][idir] > Lmax-L) {
          const double x_sp = (solver->x_cv[icv][idir] - (Lmax-L))/L;
          // sponge profile = a*x^2 + b*x^8..
          // const double sponge_coeff = sigma*(0.068*x_sp*x_sp + 0.845*pow(x_sp,8));
          // quadratic sponge profile (see JCP 212 681-702 (2006) & JCP 231 704-716 (2012))
          const double sponge_coeff = sponge_strength*x_sp*x_sp;
          cont[icv] -= sponge_coeff*solver->vol_cv[icv]*(solver->rho_uncorrected[icv]-solver->rho[icv]);
        }
      }

      break;
    }
    case SPONGE_LTH:
    {
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
          cont[icv] -= sponge_coeff*solver->vol_cv[icv]*(solver->rho_uncorrected[icv]-solver->rho[icv]);
        }
      }
      break;
    }

    default:
      assert(0);
  }

}

void SpongeHBc::computeSigma(double *sigma) const {
   switch(sponge_type) {
    case SPONGE_LX:
    case SPONGE_LY:
    case SPONGE_LZ:
    {
      int idir;
      if (sponge_type == SPONGE_LX) {
        idir = 0;
      } else if (sponge_type == SPONGE_LY) {
        idir = 1;
      } else {
        assert(sponge_type == SPONGE_LZ);
        idir = 2;
      }

      const double Lmax = sponge_data[idir];
      for (int icv = 0; icv < solver->ncv; ++icv) {
        if (solver->x_cv[icv][idir] > Lmax-L) {
          const double x_sp = (solver->x_cv[icv][idir] - (Lmax-L))/L;
          // sponge profile = a*x^2 + b*x^8..
          // const double sponge_coeff = sigma*(0.068*x_sp*x_sp + 0.845*pow(x_sp,8));
          // quadratic sponge profile (see JCP 212 681-702 (2006) & JCP 231 704-716 (2012))
          sigma[icv] += sponge_strength*x_sp*x_sp;
        }
      }
      break;
    }
    case SPONGE_LTH:
    {
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
          sigma[icv] = sponge_strength*x_sp*x_sp;
        }
      }
      break;
    }
    default:
      assert(0);
  }
}
  

bool SlipWallModelHBc::b_computeSlipLength = false;

void SlipWallModelHBc::initData() {

  assert(mf   == NULL);  mf       = new double[zone_ptr->nbf];
  assert(u_bc == NULL);  u_bc     = new double[zone_ptr->nbf][3];
  assert(cdel_w == NULL);  cdel_w = new double[zone_ptr->nbf];
  assert(rhou_s == NULL);  rhou_s = new double[zone_ptr->nbf][3];

  couple_fax = getDoubleParam("WM_SLIP_WALL_COUPLE_FAX", 0.5);

  //registered data will be overwritten upon reading
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    cdel_w[ibf] = 0.0;
    mf[ibf] = 0.0;
    for (int i =0; i < 3; ++i)
      rhou_s[ibf][i] = 0.0;
  }

  double uwall[3] = {0.0, 0.0, 0.0};
  double xrot[3]  = {0.0, 0.0, 0.0};
  double omega[3] = {0.0, 0.0, 0.0};
  double urot     = 0.0 ;

  bool isRotating = false;
  bool isStationaryInRotatingFrame = false;
  Param * param = getParam(getName());
  if ( param->size() > 1 ){
    string token = param->getString(1);
    if ( token=="U_WALL"){
      uwall[0]  = param->getDouble(2);
      uwall[1]  = param->getDouble(3);
      uwall[2]  = param->getDouble(4);
    }
    else if ( token=="ROTATING"){
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
    else if ( (token == "STATIONARY") || (token == "STATIONARY_FRAME")) {
      if ( solver->frame_rotation != NULL) {
        isStationaryInRotatingFrame = true;
      }
    }
  }

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    if (isRotating) {
      double r[3];
      FOR_I3 r[i] = zone_ptr->x_bf[ibf][i] - xrot[i];
      double omega_cross_r[3] = CROSS_PRODUCT(omega,r);
      FOR_I3 u_bc[ibf][i] = urot * omega_cross_r[i];
    }
    else if (isStationaryInRotatingFrame){
      double r[3];
      FOR_I3 r[i] = zone_ptr->x_bf[ibf][i] - solver->frame_rotation[i+3];
      u_bc[ibf][0] = -solver->frame_rotation[1]*r[2] +
                      solver->frame_rotation[2]*r[1];

      u_bc[ibf][1] = -solver->frame_rotation[2]*r[0] +
                      solver->frame_rotation[0]*r[2];

      u_bc[ibf][2] = -solver->frame_rotation[0]*r[1] +
                      solver->frame_rotation[1]*r[0];
    }
    else {
      FOR_I3 u_bc[ibf][i] = uwall[i];
    }
  }

  SlipWallModelHBc::b_computeSlipLength = true;

  //parseScalarBc(param,false); TODO const val scalar BCs not supported yet, mf!=0
}

void SlipWallModelHBc::initFromCoarseGrid(StaticSolver::CoarseGrid* cg) {

  assert(mf   == NULL);  mf       = new double[cg->ncb];
  assert(cdel_w == NULL);  cdel_w = new double[cg->ncb];

  for (int ibf = 0; ibf < cg->ncb; ++ibf) {
    mf[ibf]     = 0.0;
    cdel_w[ibf] = 0.0;
  }

}

void SlipWallModelHBc::resetBc() {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    cdel_w[ibf] = 0.0;
    mf[ibf]     = 0.0;
    for (int i =0; i < 3; ++i)
      rhou_s[ibf][i] = 0.0;
  }
}

void SlipWallModelHBc::syncPostState(){
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    //No wall motion in mf, problematic for MRF regions, DP
    mf[ibf] = DOT_PRODUCT(zone_ptr->n_bf[ibf],rhou_s[ibf]);
  }
}


void SlipWallModelHBc::setBc(){

  // DAP: computing slip length here; 
  //      computeSlipLength should only be run once for all slip wall model zones
  if (SlipWallModelHBc::b_computeSlipLength){
    SlipWallModelHBc::b_computeSlipLength = false;
    solver->computeSlipLength();
  }

  //  update rhou_s and mf once at beginning of timestep using solution from
  //  previous timestep
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv0        = zone_ptr->cvobf[ibf];
    const double area     = zone_ptr->area_bf[ibf];
    //const double cdel_aod = zone_ptr->area_over_delta_bf[ibf] * cdel_w[ibf];

    // use a delta from the cv values ... there appears to be a large variation
    // in the delta implied by area_over_delta_bf --- to be revisited ..

    const double delta_   = 0.5*pow(solver->vol_cv[icv0], 1.0/3.0);
    const double cdel_aod = zone_ptr->area_bf[ibf] * cdel_w[ibf]/delta_;
    assert ( cdel_aod >= 0.0) ;

    const double Astar    = cdel_aod/(area + cdel_aod);

    // in order to prevent the onset of 2\delta waves through the bc, we need to couple
    // the mf of this face with adjacent faces -- currently set through couple_fax.

    double u_avg[3] = {0.0,0.0,0.0};
    for (int coc = solver->cvocv_i[icv0]+1; coc != solver->cvocv_i[icv0+1]; ++coc) {
      const int icv_nbr =  solver->cvocv_v[coc];
      for (int i = 0; i < 3; ++i)
        u_avg[i] += solver->u[icv_nbr][i];
    }

    for (int i = 0; i < 3; ++i) {
      u_avg[i] /= double(solver->cvocv_i[icv0+1]-solver->cvocv_i[icv0]-1);
      u_avg[i]  = couple_fax*u_avg[i] + (1.0-couple_fax)*solver->u[icv0][i];
    }

    // recompute the rhouslip.. and the mf.  (neglecting near wall density grad )
    const double Amod = area/(area + cdel_aod);
    FOR_I3 rhou_s[ibf][i]  = solver->rho[icv0]*(Astar*u_avg[i] + Amod*u_bc[ibf][i]);

    // include wall motion in the mass flux calc
    // this seems beneficial for stability in
    // tire/ground contact patch regions
    //mf[ibf] = DOT_PRODUCT(zone_ptr->n_bf[ibf],rhou_s[ibf]);

    //store rhou_s as relative to moving surface
    FOR_I3 rhou_s[ibf][i] -= solver->rho[icv0]*u_bc[ibf][i];

    //No wall motion in mf, problematic for MRF regions, DP
    mf[ibf] = DOT_PRODUCT(zone_ptr->n_bf[ibf],rhou_s[ibf]);
  }

  //MiscUtils::dumpRange(rhou_s,zone_ptr->nbf,getName() +":rhou_s");
 
}

void SlipWallModelHBc::restrictBc(StaticSolver::CoarseGrid* cg,HelmholtzBc* bc) {

  assert(mf != NULL);
  assert(cdel_w != NULL);
  cg->restrictExtrinsicCbData(mf,bc->mf,zone_ptr->index);
  SlipWallModelHBc * slip_bc = dynamic_cast<SlipWallModelHBc*>(bc);
  cg->restrictCbData(cdel_w,slip_bc->cdel_w,zone_ptr->index);

}



void SlipWallModelHBc::addMassFlux(double * rhs) const {

  assert (solver != NULL);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];
    rhs[icv0]     += mf[ibf];
  }

  //ok to compute slip length again at start of next time step...
  SlipWallModelHBc::b_computeSlipLength = true;

}

void SlipWallModelHBc::addMomentumFlux(double * A,double (*rhs)[3]) const {

  assert( solver != NULL);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0        = zone_ptr->cvobf[ibf];
    const int coc00       = solver->cvocv_i[icv0];
    //const double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])*zone_ptr->area_over_delta_bf[ibf];
    const double mu_coeff = (solver->mu_lam[icv0] )*zone_ptr->area_over_delta_bf[ibf];

    // use a choice of the length scale that is consistent with the delta from before ...
    // in this case that is present the volume based length scale.

    //const double cdel_aod = zone_ptr->area_over_delta_bf[ibf] * cdel_w[ibf];

    const double area     = zone_ptr->area_bf[ibf];
    const double delta_   = 0.5*pow(solver->vol_cv[icv0], 1.0/3.0);
    const double cdel_aod = zone_ptr->area_bf[ibf] * cdel_w[ibf]/delta_;
    assert ( cdel_aod >= 0.0) ;

    const double Astar    = cdel_aod/(area + cdel_aod);
    //A[coc00] += mf[ibf]*Astar ; // neglecting near wall density grad...

    if ( mf[ibf] > 0.0) {
      A[coc00] += mf[ibf]; // neglecting near wall density grad...
    } else {
      A[coc00] += mf[ibf]*Astar;

      FOR_I3 rhs[icv0][i] -= mf[ibf]*(1.0 - Astar)*u_bc[ibf][i];
    }

    A[coc00] += mu_coeff * ( 1.0 - Astar ) ;

    FOR_I3 rhs[icv0][i] += mu_coeff * ( 1.0 - Astar)*u_bc[ibf][i];

  }
}

void SlipWallModelHBc::addMomentumFlux(double *A,double *At,double (*rhs)[3]) const {

  assert( solver != NULL);

  if (rhs != NULL) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv0        = zone_ptr->cvobf[ibf];
      const int coc00       = solver->cvocv_i[icv0];
      //const double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])*zone_ptr->area_over_delta_bf[ibf];
      const double mu_coeff = solver->mu_lam[icv0]*zone_ptr->area_over_delta_bf[ibf];

      // use a choice of the length scale that is consistent with the delta from before ...
      // in this case that is present the volume based length scale.

      const double area     = zone_ptr->area_bf[ibf];
      const double delta_   = 0.5*pow(solver->vol_cv[icv0], 1.0/3.0);
      //const double cdel_aod = zone_ptr->area_over_delta_bf[ibf] * cdel_w[ibf];
      const double cdel_aod = zone_ptr->area_bf[ibf] * cdel_w[ibf]/delta_;
      assert ( cdel_aod >= 0.0) ;
      const double Astar    = cdel_aod/(area + cdel_aod);

      if ( mf[ibf] > 0.0) {
        A[coc00] += mf[ibf]; // neglecting near wall density grad...
        At[coc00] += mf[ibf];
      }
      else {
        A[coc00] += mf[ibf]*Astar;
        At[coc00] += mf[ibf]*Astar;

        FOR_I3 rhs[icv0][i] -= mf[ibf]*(1.0 - Astar)*u_bc[ibf][i];
      }

      A[coc00] += mu_coeff * ( 1.0 - Astar ) ;
      At[coc00] += mu_coeff * ( 1.0 - Astar ) ;

      FOR_I3 rhs[icv0][i] += mu_coeff * ( 1.0 - Astar)*u_bc[ibf][i];
    }
  }
  else {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv0        = zone_ptr->cvobf[ibf];
      const int coc00       = solver->cvocv_i[icv0];
      //const double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])*zone_ptr->area_over_delta_bf[ibf];
      const double mu_coeff = solver->mu_lam[icv0]*zone_ptr->area_over_delta_bf[ibf];

      // use a choice of the length scale that is consistent with the delta from before ...
      // in this case that is present the volume based length scale.

      const double area     = zone_ptr->area_bf[ibf];
      const double delta_   = 0.5*pow(solver->vol_cv[icv0], 1.0/3.0);
      //const double cdel_aod = zone_ptr->area_over_delta_bf[ibf] * cdel_w[ibf+zone_ptr->ibf_f];
      const double cdel_aod = zone_ptr->area_bf[ibf] * cdel_w[ibf+zone_ptr->ibf_f]/delta_;
      assert ( cdel_aod >= 0.0) ;
      const double Astar    = cdel_aod/(area + cdel_aod);

      if ( mf[ibf+zone_ptr->ibf_f] > 0.0) {
        A[coc00] += mf[ibf+zone_ptr->ibf_f]; // neglecting near wall density grad...
        At[coc00] += mf[ibf+zone_ptr->ibf_f];
      }
      else {
        A[coc00] += mf[ibf+zone_ptr->ibf_f]*Astar;
        At[coc00] += mf[ibf+zone_ptr->ibf_f]*Astar;
      }

      A[coc00] += mu_coeff * ( 1.0 - Astar ) ;
      At[coc00] += mu_coeff * ( 1.0 - Astar ) ;
    }

  }

}

void SlipWallModelHBc::force_bf(double (*f_bf)[9]) const {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv0 = zone_ptr->cvobf[ibf];
    //double mu_coeff = (solver->mu_lam[icv0] + solver->mu_sgs[icv0])*zone_ptr->area_over_delta_bf[ibf];
    double mu_coeff = (solver->mu_lam[icv0] )*zone_ptr->area_over_delta_bf[ibf];

    const double area     = zone_ptr->area_bf[ibf];
    //const double cdel_aod = zone_ptr->area_over_delta_bf[ibf] * cdel_w[ibf];
    const double delta_   = 0.5*pow(solver->vol_cv[icv0], 1.0/3.0);
    const double cdel_aod = zone_ptr->area_bf[ibf] * cdel_w[ibf]/delta_;

    assert ( cdel_aod >= 0.0) ;

    const double Astar    = cdel_aod/(area + cdel_aod);

    mu_coeff *= (1.0 - Astar);

    //double mf_coeff = cdel_aod/(area + cdel_aod);
    //double mf_coeff_bc = area/(area + cdel_aod); //moving wall contribution

    FOR_I3 {
      f_bf[ibf][3*0+i] = mf[ibf]*((1.0-Astar)*u_bc[ibf][i]+Astar*solver->u[icv0][i]);
      f_bf[ibf][3*1+i] = solver->p[icv0]*zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = mu_coeff*(solver->u[icv0][i] - u_bc[ibf][i]);
    }
  }
  addPressureForceDueToFrameRotation(f_bf);
}

void SlipWallModelHBc::query(const string& param_str) const {

  // report min/max slip length and velocity
  double my_buf[8] = {HUGE_VAL,HUGE_VAL,HUGE_VAL,HUGE_VAL,HUGE_VAL,HUGE_VAL,HUGE_VAL,HUGE_VAL};
  double buf[8];

  double (*f_bf)[9] = new double [zone_ptr->nbf][9];
  this->force_bf(f_bf);

  double my_buf_2[3] = {0.0, 0.0, 0.0}; // area, int (tau_wall) dA, y+
  double buf_2[3];

  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    my_buf[0] = min(my_buf[0],cdel_w[ibf]);
    my_buf[1] = min(my_buf[1],-cdel_w[ibf]);

    FOR_I3 {
       my_buf[2+i*2] = min(my_buf[2+i*2],rhou_s[ibf][i]);
       my_buf[3+i*2] = min(my_buf[3+i*2],-rhou_s[ibf][i]);
    }

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

    my_buf_2[0] += zone_ptr->area_bf[ibf];
    my_buf_2[1] += zone_ptr->area_bf[ibf]*tau_wall;
    my_buf_2[2] += zone_ptr->area_bf[ibf]*y_plus;
  }
  delete[] f_bf;

  MPI_Reduce(my_buf,buf,8,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
  MPI_Reduce(my_buf_2,buf_2,3,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if (mpi_rank == 0) {
    *ss << " > " << std::left << std::setw(30) << zone_ptr->getName() << " :  time, int tau_wall dA, avg y_plus, dumpRange(cdel_w), dumpRange(rhou_s) = " << std::right
        << std::setw(8) << solver->time << " "
        << std::setw(12) << buf_2[1] << " "
        << std::setw(12) << buf_2[2]/buf_2[0] << " "
        << std::setw(12) << buf[0] << " " << std::setw(12) << -buf[1] << " "
        << std::setw(12) << buf[2] << " " << std::setw(12) << -buf[3] << " "
        << std::setw(12) << buf[4] << " " << std::setw(12) << -buf[5] << " "
        << std::setw(12) << buf[6] << " " << std::setw(12) << -buf[7] << endl;
  }
  flush();

}


CtiRegister::CtiDataError SlipWallModelHBc::funcEvalCtiData(CtiRegister::CtiData& v,
    const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if ( zone_ptr->isLeadingStrMatched(name) ) {

    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix.  we have to check for the actual data now..

    const string tau_str = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str  = zone_ptr->getName() + ":" + "y_plus";
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
          double (*f_bf)[9] = new double[zone_ptr->nbf][9];
          force_bf(f_bf);
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
            double mag_n = MAG(zone_ptr->n_bf[ibf]);
            double unit_n[3];
            FOR_I3 unit_n[i] = zone_ptr->n_bf[ibf][i]/mag_n;
            double f_p[3];
            //get wall parallel component of force
            FOR_I3 f_p[i] = f_bf[ibf][3*0+i] + f_bf[ibf][3*2+i]; //ignore (wall normal) pressure term
            FOR_I3 f_p[i] -= DOT_PRODUCT(f_p,unit_n)*unit_n[i];

            //const int icv           = zone_ptr->cvobf[ibf];

            tmp[ibf] = MAG(f_p)/zone_ptr->area_bf[ibf];
          }
          delete[] f_bf;
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

          double (*f_bf)[9] = new double[zone_ptr->nbf][9];
          force_bf(f_bf);
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
            double mag_n = MAG(zone_ptr->n_bf[ibf]);
            double unit_n[3];
            FOR_I3 unit_n[i] = zone_ptr->n_bf[ibf][i]/mag_n;

            double f_p[3];
            //get wall parallel component of force
            FOR_I3 f_p[i] = f_bf[ibf][3*0+i] + f_bf[ibf][3*2+i]; //ignore (wall normal) pressure term
            FOR_I3 f_p[i] -= DOT_PRODUCT(f_p,unit_n)*unit_n[i];

            // scale appropriately by local velocity
            const int icv      = zone_ptr->cvobf[ibf];
            const double local_u[3] = DIFF(solver->u[icv],u_bc[ibf]);
            const double mag_u = MAG(local_u);

            double unit_u[3];
            if ( mag_u > 0.0) {
              for (int i = 0; i < 3; ++i) unit_u[i] = local_u[i]/mag_u;

              tmp[ibf] = MAG(f_p)*DOT_PRODUCT(dir,unit_u)/zone_ptr->area_bf[ibf];
            }
            else {
              tmp[ibf] = 0.0;
            }
          }
          delete[] f_bf;
        }
        else return CTI_DATA_ARG_COUNT;
      }

      return CTI_DATA_OK;

    } else if ( name == yp_str) {

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      if ( b_eval_func) {

        double (*f_bf)[9] = new double[zone_ptr->nbf][9];
        force_bf(f_bf);
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          const int icv           = zone_ptr->cvobf[ibf];
          const double y1         = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];

          double mag_n = MAG(zone_ptr->n_bf[ibf]);
          double unit_n[3];
          FOR_I3 unit_n[i] = zone_ptr->n_bf[ibf][i]/mag_n;
          double f_p[3];
          //get wall parallel component of force
          FOR_I3 f_p[i] = f_bf[ibf][3*0+i] + f_bf[ibf][3*2+i]; //ignore (wall normal) pressure term
          FOR_I3 f_p[i] -= DOT_PRODUCT(f_p,unit_n)*unit_n[i];

          const double tau_w = MAG(f_p)/zone_ptr->area_bf[ibf];
          tmp[ibf]           = y1* sqrt(tau_w*solver->rho[icv])/solver->mu_lam[icv];
        }
        delete[] f_bf;

      }

      return CTI_DATA_OK;

    } else if ( name == bl_delta_str) {

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
          blde = new BoundaryLayerDataExchanger<HelmholtzSolver>(zone_ptr,solver,nbl,lbl);

        blde->computeBlFromPt(tmp,pt_cut*pt_inf);


      }

      return CTI_DATA_OK;
    }

  }

  return HelmholtzBc::funcEvalCtiData(v,name,args,b_eval_func);
}


void AlgebraicWallModelHBc::initData() {

  assert( solver != NULL);
  assert(u_bc == NULL);  u_bc = new double[zone_ptr->nbf][3];
  assert(tau_wall == NULL); tau_wall = new double[zone_ptr->nbf];

  double uwall[3] = {0.0, 0.0, 0.0};
  double xrot[3]  = {0.0, 0.0, 0.0};
  double omega[3] = {0.0, 0.0, 0.0};
  double urot     = 0.0 ;

  bool isTranslating = false;
  bool isRotating = false;
  bool isStationaryInRotatingFrame = false;
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
    else if (token=="ROTATING"){

      FOR_I3 xrot[i] = param->getDouble(iarg++);
      FOR_I3 omega[i] = param->getDouble(iarg++);
      double mag = sqrt(DOT_PRODUCT(omega,omega));
      assert ( mag > 0.0) ;
      FOR_I3 omega[i] /= mag ;
      urot = param->getDouble(iarg)*M_PI*2.0 ;
      isRotating = true;
    }
    else if ( (token == "STATIONARY") || (token == "STATIONARY_FRAME")) {
      if ( solver->frame_rotation != NULL) {
        isStationaryInRotatingFrame = true;
      }
    }
    else if (token=="ROUGHNESS_HEIGHT"){
      z0 =  param->getDouble(iarg++);
      if (z0<=0){
        CERR("invalid roughness height " << z0 << " for WM_ALG_WALL, must be greater than zero.");
      }
      COUT1("Found roughness height " << z0 << " for " << getName() << " WM_ALG_WALL, using rough wall algebraic closure");
    }
    if (isTranslating && isRotating){
      CERR("U_WALL and ROTATING set for WM_ALG_WALL, please set only one moving wall option.");
    }
  }

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    if (isRotating){
      double r[3];
      FOR_I3 r[i] = zone_ptr->x_bf[ibf][i] - xrot[i];
      double omega_cross_r[3] = CROSS_PRODUCT(omega,r);
      FOR_I3 u_bc[ibf][i] = urot * omega_cross_r[i];
    }
    else if (isStationaryInRotatingFrame){
      double r[3];
      FOR_I3 r[i] = zone_ptr->x_bf[ibf][i] - solver->frame_rotation[i+3];
      u_bc[ibf][0] = -solver->frame_rotation[1]*r[2] +
                      solver->frame_rotation[2]*r[1];

      u_bc[ibf][1] = -solver->frame_rotation[2]*r[0] +
                      solver->frame_rotation[0]*r[2];

      u_bc[ibf][2] = -solver->frame_rotation[0]*r[1] +
                      solver->frame_rotation[1]*r[0];
    }
    else{
      FOR_I3 u_bc[ibf][i] = uwall[i];
    }

    tau_wall[ibf] = 0.0; //set to zero for now, requires solver->u be set for a better guess
                         //will to this in initialHook
  }

  parseScalarBc(param,false);

}

void AlgebraicWallModelHBc::initFromCoarseGrid(StaticSolver::CoarseGrid* cg) {

  assert(u_bc == NULL);  u_bc = new double[cg->ncb][3];
  assert(tau_wall == NULL); tau_wall = new double[cg->ncb];
  for (int ibf = 0; ibf < cg->ncb; ++ibf) {
    FOR_I3 u_bc[ibf][i] = 0.0;
    tau_wall[ibf] = 0.0;
  }

}

void AlgebraicWallModelHBc::initialHook(){
    //set initial guess for tau_wall, if it wasn't read in
  if (!zone_ptr->checkDataFlag("tau_wall")){
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv0 = zone_ptr->cvobf[ibf];
      double du[3]    = DIFF(solver->u[icv0],u_bc[ibf]); //velocity relative to moving wall
      double mag_nbf  = MAG(zone_ptr->n_bf[ibf]);
      double du_dot_nbf = DOT_PRODUCT(du,zone_ptr->n_bf[ibf])/mag_nbf;
      FOR_I3 du[i]   -= du_dot_nbf*zone_ptr->n_bf[ibf][i]/mag_nbf; //velocity tangent to wall
      double du_mag   = MAG(du);
      double mu = getDoubleParam("MU");
      tau_wall[ibf] = mu*du_mag*zone_ptr->area_over_delta_bf[ibf]/zone_ptr->area_bf[ibf];
    }
  }
}

void AlgebraicWallModelHBc::restrictBc(StaticSolver::CoarseGrid* cg,HelmholtzBc* bc) {

  assert(u_bc != NULL);
  assert(tau_wall != NULL);
  AlgebraicWallModelHBc * awm_bc = dynamic_cast<AlgebraicWallModelHBc*>(bc);
  cg->restrictCbData(u_bc,awm_bc->u_bc,zone_ptr->index);
  cg->restrictCbData(tau_wall,awm_bc->tau_wall,zone_ptr->index);

}

void AlgebraicWallModelHBc::addMomentumFlux(double * A,double (*rhs)[3]) const {


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

void AlgebraicWallModelHBc::addMomentumFlux(double * A,double * At,double (*rhs)[3]) const {

  if (rhs != NULL){
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
      At[coc00] += mu_coeff;
      FOR_I3 rhs[icv0][i] += mu_coeff*u_bc[ibf][i];
    }
  }
  else{
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      const int icv0 = zone_ptr->cvobf[ibf];
      const int coc00 = solver->cvocv_i[icv0];

      double du[3]    = DIFF(solver->u[icv0],u_bc[zone_ptr->ibf_f+ibf]); //velocity relative to moving wall
      double mag_nbf  = MAG(zone_ptr->n_bf[ibf]);
      double du_dot_nbf = DOT_PRODUCT(du,zone_ptr->n_bf[ibf])/mag_nbf;
      FOR_I3 du[i]   -= du_dot_nbf*zone_ptr->n_bf[ibf][i]/mag_nbf; //velocity tangent to wall
      double du_mag   = MAG(du);

      double u_tau    = sqrt(tau_wall[ibf+zone_ptr->ibf_f]/solver->rho[icv0]); //initial guess
      if (z0<0){//smooth wall
        tau_wall[ibf+zone_ptr->ibf_f] = AlgebraicWM::solve_tau(du_mag, zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf], solver->rho[icv0], solver->mu_lam[icv0], u_tau);
      }
      else{//rough wall
        tau_wall[ibf+zone_ptr->ibf_f] = AlgebraicRoughnessWM::solve_tau(du_mag,zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf], z0, solver->rho[icv0], solver->mu_lam[icv0],u_tau);
      }


      //TODO keep mu_sgs??
      const double mu_coeff = (tau_wall[ibf+zone_ptr->ibf_f]*zone_ptr->area_bf[ibf]/du_mag);
      A[coc00] += mu_coeff;
      At[coc00] += mu_coeff;
    }

  }


}

void AlgebraicWallModelHBc::force_bf(double (*f_bf)[9]) const {

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
  addPressureForceDueToFrameRotation(f_bf);
}

void AlgebraicWallModelHBc::query(const string& param_str) const {

  // report min/max slip length and velocity
  double my_buf[3] = {0.0,0.0,0.0}; // area, int tau_wall dA, y_plus
  double buf[3];

  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv0 = zone_ptr->cvobf[ibf];

    const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
    const double u_tau = sqrt( tau_wall[ibf]/solver->rho[icv0]);
    const double nu    = solver->mu_lam[icv0] / solver->rho[icv0];
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

CtiRegister::CtiDataError AlgebraicWallModelHBc::funcEvalCtiData(CtiRegister::CtiData& v,
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
          blde = new BoundaryLayerDataExchanger<HelmholtzSolver>(zone_ptr,solver,nbl,lbl);

        blde->computeBlFromPt(tmp,pt_cut*pt_inf);

      }

      return CTI_DATA_OK;
    }
  }

  return HelmholtzBc::funcEvalCtiData(v,name,args,b_eval_func);
}
