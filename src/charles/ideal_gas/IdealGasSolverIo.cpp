
#include "IdealGasSolver.hpp"
#include "bcprofile/ProfileReader.hpp"

// convenience macro to loop over the boundary conditions...
#define FOR_BCZONE for(vector<IdealGasBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)

void IdealGasSolver::initFromParams() {

  // parameter-based initialization of data...

  FlowSolver::initFromParams();

  if (Param * param = getParam("INIT_FROM_PROFILE")) {
    COUT1(" > initializing entire flowfield from INIT_FROM_PROFILE; this will over-write INIT_CV_DATA_IN_GEOM or INIT commands");
    if (step > 0) {
      CWARN("Calling INIT_FROM_PROFILE for step > 0");
    }

    int iarg = 0;
    int type = -1;
    double type_double[2] = {0.0,0.0};
    double x0[3] = {0.0,0.0,0.0};
    double n0[3] = {0.0,0.0,0.0};
    const int imax = param->size();

    FluentBPReader profile;

    while (iarg < imax) {
      const string token = param->getString(iarg++);
      if (token == "FILE") {
        profile.init(param->getString(iarg++));
      }
      else if (token == "TYPE") {
        const string _type = param->getString(iarg++);
        if (_type == "RUP") {
          type = 0;
        }
        else if (_type == "UPT") {
          type = 1;
        }
        else if (_type == "U_CONSTANT_TP") {
          type = 2;
          type_double[0] = param->getDouble(iarg++);
          type_double[1] = param->getDouble(iarg++);
        }
      }
      else if (token == "AXIS") {
        FOR_I3 x0[i] = param->getDouble(iarg++);
        FOR_I3 n0[i] = param->getDouble(iarg++);
      }
      else if (token == "FROM_MEAN") {
        profile.setUseMean(true);
      }
      else {
        CWARN("unrecognized INIT_FROM_PROFILE param: " << token);
      }
    }

    if (!profile.isInitialized()) {
      CERR("boundary profile was not properly initialized");
    }

    if (type == -1) {
      CERR("INIT_FROM_PROFILE \"TYPE\" must be specified");
    }

    COUT1(" > initializing profile data");
    if (profile.getType() == LINE_R) {
      COUT1("    > radial along axis x0: " << COUT_VEC(x0) << " " << COUT_VEC(n0));

      // XXXX this is wrong..
      // make sure an axis has been specified
      const double nsum = n0[0]+n0[1]+n0[2];
      if (nsum == 0.0) {
        CERR("INIT_FROM_PROFILE requires specification of an axis");
      }
      else {
        NORMALIZE(n0);
      }

      // for radial, use distance to specified axis as radius
      double * my_r = new double[ncv];
      FOR_ICV {
        my_r[icv] = MiscUtils::getPointToLineDist(x_cv[icv],x0,n0);
      }
      profile.setPoints(my_r,ncv);
      DELETE(my_r);
    }
    else {
      assert(profile.getType() == POINT_3D);
      COUT1("    > using point specification");
      profile.setPointsAdt(&*x_cv,x0,ncv);
    }
    //TODO add support for wall-distance....

    if (type == 0) {
      // RUP
      profile.checkVar("x-velocity");
      profile.checkVar("y-velocity");
      profile.checkVar("z-velocity");
      profile.checkVar("absolute-pressure");
      profile.checkVar("density");
      FOR_ICV {
        rho[icv] = profile.getData(icv,"density");
        u[icv][0] = profile.getData(icv,"x-velocity");
        u[icv][1] = profile.getData(icv,"y-velocity");
        u[icv][2] = profile.getData(icv,"z-velocity");
        rhoE[icv] = profile.getData(icv,"absolute-pressure")/(gamma-1.0) + 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);
      }

    } else if (type == 1) {
      // UPT
      profile.checkVar("x-velocity");
      profile.checkVar("y-velocity");
      profile.checkVar("z-velocity");
      profile.checkVar("absolute-pressure");
      profile.checkVar("temperature");
      FOR_ICV {
        const double p_icv = profile.getData(icv,"absolute-pressure");
        rho[icv] = p_icv/(R_gas*profile.getData(icv,"temperature"));
        u[icv][0] = profile.getData(icv,"x-velocity");
        u[icv][1] = profile.getData(icv,"y-velocity");
        u[icv][2] = profile.getData(icv,"z-velocity");
        rhoE[icv] = p_icv/(gamma-1.0) + 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);
      }

    } else if (type == 2) {
      // U_CONSTANT_TP
      profile.checkVar("x-velocity");
      profile.checkVar("y-velocity");
      profile.checkVar("z-velocity");
      const double _rho = type_double[1]/(R_gas*type_double[0]);
      FOR_ICV {
        rho[icv] = _rho;
        u[icv][0] = profile.getData(icv,"x-velocity");
        u[icv][1] = profile.getData(icv,"y-velocity");
        u[icv][2] = profile.getData(icv,"z-velocity");
        rhoE[icv] = type_double[1]/(gamma-1.0) + 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);
      }

    } else {
      CERR("unrecognized type somehow set...");
    }

    setDataFlag("rho",1);
    setDataFlag("u",1);
    setDataFlag("rhoE",1);

    // sync everyone in case the bcs need these values.

    updateConservativeAndPrimitiveData();
    calcSgsAndMaterialProperties();

    RestartHashUtilities::slesSetHashFromParam(param);

  } else if (Param * param = getParam("INIT_RUP")) {
    COUT1(" > initializing entire flowfield from INIT_RUP; this will over-write INIT_CV_DATA_IN_GEOM or INIT commands");
    if (step > 0) {
      CWARN("Calling INIT_RUP for step > 0");
    }

    const double init_rho = param->getDouble(0);
    const double init_u   = param->getDouble(1);
    const double init_v   = param->getDouble(2);
    const double init_w   = param->getDouble(3);
    const double init_p   = param->getDouble(4);

    bool b_mrf_transform = false;
    if ( (param->size() == 6) && (param->getString(5) == "STATIONARY_FRAME")) {

      if ( frame_rotation != NULL)
        b_mrf_transform = true;

    }

    COUT1("INIT_RUP: " << init_rho << " " << init_u << " " << init_v << " " << init_w
        << " " << init_p);

    if ( b_mrf_transform) {

      assert( frame_rotation); // checked 10 lines above ..

      for (int icv = 0; icv < ncv; ++icv) {

        u[icv][0] = init_u - frame_rotation[1]*x_cv[icv][2] + frame_rotation[2]*x_cv[icv][1];
        u[icv][1] = init_v - frame_rotation[2]*x_cv[icv][0] + frame_rotation[0]*x_cv[icv][2];
        u[icv][2] = init_w - frame_rotation[0]*x_cv[icv][1] + frame_rotation[1]*x_cv[icv][0];
        rho[icv] = init_rho;
        rhoE[icv] = init_p/(gamma-1.0) + 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);

      }

    } else {

      FOR_ICV {
        rho[icv]  = init_rho;
        u[icv][0] = init_u;
        u[icv][1] = init_v;
        u[icv][2] = init_w;
        rhoE[icv] = init_p/(gamma-1.0) + 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);
      }
    }

    setDataFlag("rho",1);
    setDataFlag("u",1);
    setDataFlag("rhoE",1);

    // sync everyone in case the bcs need these values.

    updateConservativeAndPrimitiveData();
    calcSgsAndMaterialProperties();

    // update slesHash id
    RestartHashUtilities::slesSetHashFromParam(param);

  } else if (Param * param = getParam("INIT_UPT")) {
    COUT1(" > initializing entire flowfield from INIT_UPT; this will over-write INIT_CV_DATA_IN_GEOM or INIT commands");
    if (step > 0) {
      CWARN("Calling INIT_UPT for step > 0");
    }

    const double init_u   = param->getDouble(0);
    const double init_v   = param->getDouble(1);
    const double init_w   = param->getDouble(2);
    const double init_p   = param->getDouble(3);
    const double init_T   = param->getDouble(4);

    bool b_mrf_transform = false;
    if ( (param->size() == 6) && (param->getString(5) == "STATIONARY_FRAME")) {

      if ( frame_rotation != NULL)
        b_mrf_transform = true;

    }

    COUT1("INIT_UPT: " << init_u << " " << init_v << " " << init_w
        << " " << init_p << " " << init_T);

    const double init_rho = init_p/(R_gas*init_T);

    if ( b_mrf_transform) {

      assert( frame_rotation); // checked 10 lines above ..

      for (int icv = 0; icv < ncv; ++icv) {

        u[icv][0] = init_u - frame_rotation[1]*x_cv[icv][2] + frame_rotation[2]*x_cv[icv][1];
        u[icv][1] = init_v - frame_rotation[2]*x_cv[icv][0] + frame_rotation[0]*x_cv[icv][2];
        u[icv][2] = init_w - frame_rotation[0]*x_cv[icv][1] + frame_rotation[1]*x_cv[icv][0];
        rho[icv] = init_rho;
        rhoE[icv] = init_p/(gamma-1.0) + 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);

      }

    } else {

      FOR_ICV {
        rho[icv]  = init_rho;
        u[icv][0] = init_u;
        u[icv][1] = init_v;
        u[icv][2] = init_w;
        rhoE[icv] = init_p/(gamma-1.0) + 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);
      }
    }

    setDataFlag("rho",1);
    setDataFlag("u",1);
    setDataFlag("rhoE",1);

    // update slesHash id
    RestartHashUtilities::slesSetHashFromParam(param);

    // sync everyone in case the bcs need these values.

    updateConservativeAndPrimitiveData();
    calcSgsAndMaterialProperties();

  } else if (checkParam("INIT_FROM_PF")||checkParam("INIT_FROM_POTENTIAL_FLOW")) {
    COUT1(" > initializing entire flowfield from INIT_FROM_PF; this will over-write INIT_CV_DATA_IN_GEOM or INIT commands");
    initFromPotentialFlow();

  } else if ( Param * param = getParam("INTERP_FROM_RESTART") ) {

    //TODO?
    // this may fail if multiple INTERP_FROM_RESTARTs are present, some of which are not
    // from Fluent. I suppose this case doens't need to be handled
    // IF multiple INTERPs are present, this assumes they are from the same simulation
    // i.e., same reference pressure. Probably a good assumption; do we need to say so
    // to the user?

    if (param != NULL) {
      int iarg = 0;
      const string filename = param->getString(0);

      // rescale interpolated data if not coming from charles
      if ((filename.find(".ip") != string::npos)||(filename.find(".cas") != string::npos)) {
        //TODO: assumes all interps are coming from same solver (just check first)

        // fluent  files
        IF_RANK0 cout << "modifying Fluent interpolated data for consistency with charles..." << endl;

        // fluent case might be low-ma with p_exit = 0 (for example)
        double dp = 0.0, pmin = 0.0, my_pmin = 1e+10;
        FOR_ICV if (p[icv] < my_pmin) my_pmin = p[icv];
        MPI_Allreduce(&my_pmin,&pmin,1,MPI_DOUBLE,MPI_MIN,mpi_comm);
        if (pmin <= 0.0) {
          if (p_ref > -pmin) {
            dp = p_ref;
          } else {
            CERR("pressure initialization failed: min p = " << pmin << " vs P_REF = " << p_ref);
          }
        }

        // reconstruct state vars ... u & T should already be set from fluent file
        FOR_ICV {
          p[icv]   += dp;
          rho[icv]  = p[icv]/T[icv]/R_gas;
          rhoE[icv] = p[icv]/(gamma-1.0) + 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);
        }
        IF_RANK0 cout << " > shifted pressure by delta p = " << dp << endl;

        resetScalars();

        setDataFlag("u",1);
        setDataFlag("rho",1);
        setDataFlag("rhoE",1);

        // sync everyone in case the bcs need these values.
        updateConservativeAndPrimitiveData();
        calcSgsAndMaterialProperties();

        dumpRange(rho,ncv,"rho");
        dumpRange(u,ncv,"u");
        dumpRange(p,ncv,"p");
        dumpRange(T,ncv,"T");

        // don't update hash
      }
    }
  }

  if (Param * param = getParam("HELMHOLTZ_TO_IDEAL_GAS")) {
    // not part of if/else loop above because helmholtz solutions may have been
    // provided via RESTART or multiple INTERP_FROM_RESTART on a different grid
    // this check will operate on all interpolated data, assuming they all are
    // coming from Helmholtz solution

    // need to build rhoE from variables read from Helmholtz EOS sles
    COUT1(" > building rhoE state values from Helmholtz EOS data");

    if ((!checkParam("P_REF")) || (!checkParam("RHO_REF"))) {
      CERR("P & RHO reference values must be set to properly define ideal-gas state from Helmholtz state");
    }

    const double p_ref   = getDoubleParam("P_REF");
    const double rho_ref = getDoubleParam("RHO_REF");
    bool p_is_relative   = true;

    int iarg=0;
    const int iarg_max = param->size();
    while (iarg < iarg_max) {
      const string token = param->getString(iarg++);
      if (token == "P_ABS") {
        p_is_relative = false;
      }
      else if (token == "P_REL") {
        p_is_relative = true;
      }
    }

    if (p_is_relative) {
      COUT1(" > pressure from file is assumed a perturbation about P_REF");
    }
    else {
      COUT1(" > pressure from file is absolute valued");
    }

    // rho and u should have been read from file
    //TODO make sure works when specifying restart as INTERP or on RESTART line

    const double d_p = (p_is_relative) ? p_ref:0.0;

    FOR_ICV {

      p[icv]    += d_p;

      // limit p_cv to between 0.9 and 1.1 of p_ref (appropriate for low Ma flows)...
      p[icv] = max(0.9*p_ref,min(1.1*p_ref,p[icv]));

      // use the isentropic relation rho/rho_ref = (p/p_ref)^(1/gamma)
      // overwriting rho from Helmholtz solution
      rho[icv] = rho_ref*pow(p[icv]/p_ref,1.0/gamma);

      rhoE[icv] = p[icv]/(gamma-1.0) + 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);
    }

    setDataFlag("rho",1);
    setDataFlag("u",1);
    setDataFlag("rhoE",1);

    // do not update the hash --

    // sync everyone in case the bcs need these values.

    updateConservativeAndPrimitiveData();
    calcSgsAndMaterialProperties();
  }

}

void IdealGasSolver::initFromPotentialFlow() {

  COUT1("initFromPotentialFlow()");

  // 1. solve a potential flow problem from the boundaries to build a rhou
  // field: rhou = -grad(phi)...

  double * phi = new double[ncv_g];
  for (int icv = 0; icv < ncv_g; ++icv) phi[icv] = 0.0;

  double * mdot_bf = new double[nbf];
  FOR_IBF mdot_bf[ibf] = 0.0;

  double mdot_sum = 0.0;
  double outlet_rhoun = 0.0;

  bool *b_outlet_zone = new bool[bfZoneVec.size()];

  FOR_BCZONE {
    const int izone = (*it)->zone_ptr->index;
    b_outlet_zone[izone] = false;
    const string zone_name = bfZoneVec[izone].getName();
    if (Param * param = getParam(zone_name)) {
      int iarg = 0;
      const string bc_type = param->getString(iarg++);
      if ( bc_type == "NSCBC_MT") {
        const double this_mdot = param->getDouble(iarg++);
        mdot_sum += this_mdot;
        const double this_rhoun = this_mdot/bfZoneVec[izone].area_global;
        // report...
        if (mpi_rank == 0)
          cout << " > INLET \"" << zone_name << "\": MDOT=" << this_mdot << " (area=" << bfZoneVec[izone].area_global << ", rhoun=" << this_rhoun << ")" << endl;
        for (int ibf = bfZoneVec[izone].ibf_f; ibf <= bfZoneVec[izone].ibf_l; ++ibf) {
          mdot_bf[ibf] = -this_rhoun*area_bf[ibf]; // inlet negative by convention
        }
      }
      else if ( (bc_type == "CBC_MPT")||(bc_type == "CBC_RUP")||(bc_type == "CBC_RUNP")||(bc_type == "CBC_UPT") ) {

        double this_mdot;
        if (bc_type == "CBC_MPT") {
          this_mdot = param->getDouble(iarg++);
        }
        else if (bc_type == "CBC_RUP") {
          const double this_rho = param->getDouble(iarg++);
          double this_u[3]; FOR_I3 this_u[i] = param->getDouble(iarg++);
          this_mdot = -this_rho*DOT_PRODUCT(bfZoneVec[izone].n_global,this_u); // flip normal to point in
        }
        else if (bc_type == "CBC_RUNP") {
          const double this_rho = param->getDouble(iarg++);
          const double this_un = param->getDouble(iarg++);
          this_mdot = this_rho*bfZoneVec[izone].area_global*this_un;
        }
        else {
          assert(bc_type == "CBC_UPT");
          double this_u[3]; FOR_I3 this_u[i] = param->getDouble(iarg++);
          const double this_p = param->getDouble(iarg++);
          const double this_T = param->getDouble(iarg++);
          const double this_rho = this_p/(this_T*R_gas);
          this_mdot = -this_rho*DOT_PRODUCT(bfZoneVec[izone].n_global,this_u); // flip normal to point in
        }

        if (this_mdot > 0.0) {
          // looks good...
          mdot_sum += this_mdot;
          const double this_rhoun = this_mdot/bfZoneVec[izone].area_global;
          // report...
          if (mpi_rank == 0)
            cout << " > INLET \"" << zone_name << "\": MDOT=" << this_mdot << " (area=" << bfZoneVec[izone].area_global << ", rhoun=" << this_rhoun << ")" << endl;
          for (int ibf = bfZoneVec[izone].ibf_f; ibf <= bfZoneVec[izone].ibf_l; ++ibf) {
            mdot_bf[ibf] = -this_rhoun*area_bf[ibf]; // inlet negative by convention
          }
        }
        else {
          // this is an outlet...
          outlet_rhoun += bfZoneVec[izone].area_global;
          b_outlet_zone[izone] = true;
        }
      }
      else if ( bc_type == "CBC_PROFILE" ) {
        CbcProfile* bc_ptr = static_cast<CbcProfile*>(*it);
        assert(bc_ptr->bf_state);

        double my_mdot = 0.0;
        for (int ibf = bfZoneVec[izone].ibf_f; ibf <= bfZoneVec[izone].ibf_l; ++ibf) {
          mdot_bf[ibf] = DOT_PRODUCT(bc_ptr->bf_state[ibf-bfZoneVec[izone].ibf_f].u,n_bf[ibf])/bc_ptr->bf_state[ibf-bfZoneVec[izone].ibf_f].sp_vol;
          my_mdot += -mdot_bf[ibf];
        }
        double this_mdot;
        MPI_Allreduce(&my_mdot,&this_mdot,1,MPI_DOUBLE,MPI_SUM,mpi_comm);

        if (this_mdot > 0.0) {
          // looks good...
          mdot_sum += this_mdot;
          const double this_rhoun = this_mdot/bfZoneVec[izone].area_global;
          // report...
          if (mpi_rank == 0)
            cout << " > INLET \"" << zone_name << "\": MDOT=" << this_mdot << " (area=" << bfZoneVec[izone].area_global << ", rhoun=" << this_rhoun << ")" << endl;
        }
        else {
          // this is an outlet...
          outlet_rhoun += bfZoneVec[izone].area_global;
          b_outlet_zone[izone] = true;
          // zero out the mdot_bf...
          for (int ibf = bfZoneVec[izone].ibf_f; ibf <= bfZoneVec[izone].ibf_l; ++ibf)
            mdot_bf[ibf] = 0.0;
        }

      }
      else if ((bc_type == "CBC_TOTAL_PT")||(bc_type == "NSCBC_PROFILE")||(bc_type == "INFLOW_TURBULENCE")) {
        // NscbcProfile and InflowTurbulenceIG initialize their data in initialHook because they READWRITE_DATA
        CERR("INIT_FROM_POTENTIAL_FLOW does not support " << bc_type << ". Please use a different different boundary condition for initialization.");
      }
      else if ((bc_type == "NSCBC_OUTLET_P")||(bc_type == "NSCBC_OUTLET_PRESSURE")||(bc_type == "NSCBS_OUTLET_MDOT")) {
        // example...
        // OUTLET NSCBC_OUTLET_P 3.766e5 .1 1
        // for now, accumulate areas in outlet_rhoun...
        outlet_rhoun += bfZoneVec[izone].area_global;
        b_outlet_zone[izone] = true;
      }
    }
  }

  // also add the volumetric mass source in the cone...

  double * rhs = new double[ncv];
  for (int icv = 0; icv < ncv; ++icv)
    rhs[icv] = 0.0;

  // inlets and outlets must balance in potential flow...
  assert(mdot_sum > 0.0);
  assert(outlet_rhoun > 0.0); // recall this is holding the area right now -- should be positive
  outlet_rhoun = mdot_sum/outlet_rhoun;

  // now go back through the outlet(s) and specify the mdot_bf...
  FOR_IZONE(bfZoneVec) {
    if (b_outlet_zone[izone]) {
      if (mpi_rank == 0)
        cout << " > OUTLET \"" << bfZoneVec[izone].getName() << "\": computed MDOT=" << outlet_rhoun*bfZoneVec[izone].area_global << " (area=" << bfZoneVec[izone].area_global << ", rhoun=" << outlet_rhoun << ")" << endl;
      for (int ibf = bfZoneVec[izone].ibf_f; ibf <= bfZoneVec[izone].ibf_l; ++ibf) {
        mdot_bf[ibf] = outlet_rhoun*area_bf[ibf]; // outlet positive by convention
      }
    }
  }
  delete[] b_outlet_zone;

  // ----------------------------------
  // solver...
  // ----------------------------------

  for (int ibf = 0; ibf < nbf; ++ibf) {
    const int icv = cvobf[ibf];
    rhs[icv] += mdot_bf[ibf];
  }

  // check rhs...

  double my_sum = 0.0;
  for (int icv = 0; icv < ncv; ++icv)
    my_sum += rhs[icv];
  double sum;
  MPI_Reduce(&my_sum,&sum,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  if (mpi_rank == 0)
    cout << " > sum(rhs) (should be zero): " << sum << endl;

  const double pf_zero = getDoubleParam("PF_ZERO",1.0E-8);
  const int pf_maxiter = getIntParam("PF_MAXITER",20000);

  double * A = new double[cvocv_i[ncv]];
  buildCvLaplacian(A);
  solveCvCg(phi,A,rhs,pf_zero,pf_maxiter,true); // verbose

  delete[] A;

  MiscUtils::dumpRange(phi,ncv,"phi");

  // now check divergence with corrections...

  MiscUtils::dumpRange(rhs,ncv,"rhs - before");

  for (int ifa = 0; ifa < nfa; ++ifa) {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    const double drhoun = -(phi[icv1]-phi[icv0])*area_over_delta_fa[ifa];
    rhs[icv0] += drhoun;
    if (icv1 < ncv)
      rhs[icv1] -= drhoun;
  }

  MiscUtils::dumpRange(rhs,ncv,"rhs - after");

  delete[] rhs;

  // now compute rhou = -grad(phi)...
  // use u for now...

  calcCvGradLeastSquares(u,phi,mdot_bf); // this gradient can take the boundary normal derrivative
  FOR_ICV FOR_I3 u[icv][i] = -u[icv][i]; // flip to get u sign correct

  delete[] phi;
  delete[] mdot_bf;

  bool b_u_max = false;
  double pf_u_max;
  if (Param * param = getParam("PF_U_MAX")) {
    b_u_max = true;
    pf_u_max = param->getDouble();
    if (mpi_rank == 0)
      cout << " > velocity magnitude will be limited to PF_U_MAX=" << pf_u_max << endl;
  }
  const double pf_u_ref = getDoubleParam("PF_U_REF",0.0); // assume p_ref = p_stag if not spec.

  FOR_ICV {

    rho[icv] = rho_ref;

    // convert rhou to u...
    FOR_I3 u[icv][i] /= rho[icv];

    // limit u...

    if (b_u_max) {
      const double mag_u = MAG(u[icv]);
      if (mag_u > pf_u_max) {
        FOR_I3 u[icv][i] *= pf_u_max/mag_u;
      }
    }

    // energy...

    //FOR_ICV rhoE[icv] = p_ref/(gamma-1.0) + 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);
    FOR_ICV {
      const double this_p = p_ref+0.5*rho_ref*pf_u_ref*pf_u_ref-0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);
      rhoE[icv] = this_p/(gamma-1.0) + 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);
    }

  }

  // tell solver that we initialized the conserved variables...

  setDataFlag("rho",1);
  setDataFlag("u",1);
  setDataFlag("rhoE",1);

  updateConservativeAndPrimitiveData();

}

void IdealGasSolver::queryBcs() {

  // QUERY_BC <bc-name> [INTERVAL = 1]
  FOR_PARAM_MATCHING("QUERY_BC") {

    // check if the bc is matched against a known query, we will
    // use the entire param string as the key for the bc_query map

    map<string,pair<int,IdealGasBc*> >::iterator bc_it = bc_queries.find(param->str());

    if ( bc_it == bc_queries.end() ) {

      // this is a new query-- that has not been parsed yet.

      int interval         = check_interval; // default interval is check_interval
      const string bc_name = param->getString(0);

      int iarg = 1;
      while ( iarg < param->size()) {
        string token = param->getString(iarg++);
        if ( token == "INTERVAL") {
          interval = param->getInt(iarg++);
          if (interval <= 0) {
            CWARN(" > QUERY_BC INTERVAL expects a positive integer; setting to CHECK_INTERVAL");
            // invalid value entered, so treat as unset (default value)
            interval = check_interval;
          }
        }
      }

      IdealGasBc* bc = getBc(bc_name);

      if ( bc == NULL) {

        CWARN(" > unable to find boundary zone for QUERY: " << bc_name);

      } else {

        pair<map<string,pair<int,IdealGasBc*> >::iterator,bool> ret =
          bc_queries.insert(pair<string,pair<int,IdealGasBc*> >( param->str(),
                pair<int,IdealGasBc*>(interval,bc)));

        assert( ret.second); // ensure that the query was properly inserted.
        bc_it = ret.first;
      }
    }

    if ( bc_it != bc_queries.end() ) {

      const int query_interval = bc_it->second.first;
      IdealGasBc* bc           = bc_it->second.second;

      if ( step%query_interval == 0)
        bc->query(bc_it->first);
    }

  }
}

// custom variable output via cti var evaluation; recall on completion
// if the data is not found, then it must be passed down to the base
// StaticSolver class to see if it can evaluate the data or otherwise error

CtiRegister::CtiDataError IdealGasSolver::funcEvalCtiData(CtiRegister::CtiData& v,const string& name,
                                                          list<CtiRegister::CtiData>& args,
                                                          const bool b_eval_func) {

  using namespace CtiRegister;

  //if ( mpi_rank == 0) {
  //  cout << "IdealGasSolver:: evaluating : " << name << " with n args: " << args.size() << endl;
  //}

  if ( name == "mass_flux" ) {

    // need to evaluate a flux probe of a given variable.  this is an approximation
    // of the inviscid flux (no diffusion and no considerations of convective stab)

    if ( args.size() > 1)
      return CTI_DATA_ARG_COUNT;

    // we need to compute the mass flux at the faces.  we need to check to see
    // if its already available.  we could call the getUnregisteredCtiData but
    // we can circumvent the operator parsing since we know exactly what kind of
    // data we are looking for (this also bypasses all of the error throwing
    // that the getUnregisteredCtiData does)

    //CtiData* mf_data = CtiRegister::getUnregisteredCtiData("mdot_fa");
    //double * mf      = NULL;

    map<const string,CtiData>::iterator iter = currentDataMap.find("mdot_fa");
    CtiData* mf_data                         = NULL;
    double * mf                              = NULL;

    if ( iter != currentDataMap.end())
      mf_data = &(iter->second);

    if ( !mf_data) {

      // the data does not exist yet, so we need to populate it..
      pair<map<const string,CtiData>::iterator,bool> return_pair =
        CtiRegister::currentDataMap.insert(
            pair<const string,CtiData>("mdot_fa",CtiData()));

      assert( return_pair.second);

      mf_data = &(return_pair.first->second);
      mf      = createSignedFaD1Data(*mf_data);
      assert( mf);

      if ( b_eval_func ) {
        for (int ifa = 0; ifa < nfa; ++ifa) {

          const int icv0          = cvofa[ifa][0];
          const int icv1          = cvofa[ifa][1];
          const double inv_sp_vol = 2.0*rho[icv0]*rho[icv1]/(rho[icv0] + rho[icv1]);

          double undA_avg       = 0.0;
          for (int i = 0; i < 3; ++i)
            undA_avg += 0.5*(u[icv0][i] + u[icv1][i])*n_fa[ifa][i];

          mf[ifa] = inv_sp_vol*undA_avg;
        }
      }

    } else {

      assert( mf_data->getType() == DN_DATA);
      assert( mf_data->getTopology() == SIGNED_FA_DATA);
      mf = mf_data->getDNptr();
      assert( mf);
    }

    if ( args.size() == 0 ) {

      // just return the mass flux

      double *v_ptr = createSignedFaD1Data(v);

      if ( b_eval_func) {
        for (int ifa = 0; ifa < nfa; ++ifa)
          v_ptr[ifa] = mf[ifa];
      }

      return CTI_DATA_OK;

    } else {

      list<CtiData>::iterator arg = args.begin();
      const int datatype          = arg->getType();

      if ( datatype == DN_DATA) {

        if ( arg->getTopology() != CV_DATA )
          return CTI_DATA_NOT_VALID;

        double * v_ptr = createSignedFaD1Data(v);

        if ( b_eval_func) {

          // for interprocessor/periodic boundaries we add our half of the flux
          // and start the parallel reduction

          double * arg_ptr = arg->getDNptr();
          for (int ifa = nfa_i; ifa < nfa; ++ifa) {

            const int icv0 = cvofa[ifa][0];
            v_ptr[ifa]     = mf[ifa]* 0.5* arg_ptr[icv0];

          }

          // the normal is of the other sign on the other rank...

          updateFaDataStart( v_ptr, SUBTRACT_DATA);

          // internal faces-- no ghost data required

          for (int ifa = 0; ifa < nfa_i; ++ifa) {

            const int icv0 = cvofa[ifa][0];
            const int icv1 = cvofa[ifa][1];
            v_ptr[ifa]     = mf[ifa]* 0.5*( arg_ptr[icv0] + arg_ptr[icv1]);

          }

          updateFaDataFinish( v_ptr, SUBTRACT_DATA);

        }

        return CTI_DATA_OK;

      } else if ( datatype == DN3_DATA) {

        double (*v_ptr)[3] = createFaD2Data(v);

        if ( b_eval_func) {

          // for interprocessor/periodic boundaries we add our half of the flux
          // and start the parallel reduction

          double (*arg_ptr)[3] = arg->getDN3ptr();
          for (int ifa = nfa_i; ifa < nfa; ++ifa) {

            const int icv0 = cvofa[ifa][0];
            for (int i = 0; i < 3; ++i)
              v_ptr[ifa][i] = mf[ifa]* 0.5* arg_ptr[icv0][i];

          }

          updateFaDataStart( v_ptr, SUBTRACT_ROTATE_DATA);

          // internal faces-- no ghost data required

          for (int ifa = 0; ifa < nfa_i; ++ifa) {

            const int icv0 = cvofa[ifa][0];
            const int icv1 = cvofa[ifa][1];
            for (int i =0; i < 3; ++i)
              v_ptr[ifa][i] = mf[ifa]* 0.5*( arg_ptr[icv0][i] + arg_ptr[icv1][i]);

          }

          updateFaDataFinish( v_ptr, SUBTRACT_ROTATE_DATA);

        }

        return CTI_DATA_OK;

      } else {

        return CTI_DATA_NOT_VALID; // cant evaluate a flux of this..

      }
    }

  } else if ( name == "fgr" ) {

    if (args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double *v_ptr = createCvD1Data(v);

    if ( b_eval_func) {

      for (int icv = 0; icv < ncv; ++icv) {
        v_ptr[icv] = cv_compact[icv].fgr;
      }

    }

    return CTI_DATA_OK;

  }
  else if ( name == "fgr" ) {

    if ( args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double * v_ptr = createCvD1Data(v);

    if ( b_eval_func) {
      FOR_ICV v_ptr[icv] = cv_compact[icv].fgr;
    }

    return CTI_DATA_OK;

  }
  else if ( name == "mach" ) {

    if ( args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double * v_ptr = createCvD1Data(v);

    if ( b_eval_func) {
      FOR_ICV v_ptr[icv] = MAG(u[icv])/sos[icv];
    }

    return CTI_DATA_OK;

  }
  else if ( name == "p_total" ) {

    if ( args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double * v_ptr = createCvD1Data(v);
    if ( b_eval_func) {

      if ( frame_rotation == NULL) {

        FOR_ICV {
          const double Ma2 = DOT_PRODUCT(u[icv],u[icv])/(sos[icv]*sos[icv]);
          v_ptr[icv] = p[icv]*pow(1.0 + 0.5*(gamma-1.0)*Ma2, gamma/(gamma-1.0));
        }

      } else {

        FOR_ICV {

          // transformed velocity back into a stationary reference frame..

          double u_[3];
          u_[0] = u[icv][0] + frame_rotation[1]*x_cv[icv][2] - frame_rotation[2]*x_cv[icv][1];
          u_[1] = u[icv][1] + frame_rotation[2]*x_cv[icv][0] - frame_rotation[0]*x_cv[icv][2];
          u_[2] = u[icv][2] + frame_rotation[0]*x_cv[icv][1] - frame_rotation[1]*x_cv[icv][0];

          const double Ma2 = DOT_PRODUCT(u_,u_)/(sos[icv]*sos[icv]);
          v_ptr[icv] = p[icv]*pow(1.0 + 0.5*(gamma-1.0)*Ma2, gamma/(gamma-1.0));
        }
      }
    }

    return CTI_DATA_OK;

  }
  else if ( name == "T_total" ) {

    if ( args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double * v_ptr = createCvD1Data(v);
    if ( b_eval_func) {

      if ( frame_rotation == NULL) {

        FOR_ICV {
          const double Ma2 = DOT_PRODUCT(u[icv],u[icv])/(sos[icv]*sos[icv]);
          v_ptr[icv] = T[icv]*(1.0 + 0.5*(gamma-1.0)*Ma2);
        }

      } else {

        FOR_ICV {

          // transformed velocity back into a stationary reference frame..

          double u_[3];
          u_[0] = u[icv][0] + frame_rotation[1]*x_cv[icv][2] - frame_rotation[2]*x_cv[icv][1];
          u_[1] = u[icv][1] + frame_rotation[2]*x_cv[icv][0] - frame_rotation[0]*x_cv[icv][2];
          u_[2] = u[icv][2] + frame_rotation[0]*x_cv[icv][1] - frame_rotation[1]*x_cv[icv][0];

          const double Ma2 = DOT_PRODUCT(u_,u_)/(sos[icv]*sos[icv]);
          v_ptr[icv] = T[icv]*(1.0 + 0.5*(gamma-1.0)*Ma2);

        }
      }

    }

    return CTI_DATA_OK;

  }
  else if ( name == "gr") {

    if ( args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double * v_ptr = createCvD1Data(v);

    if ( b_eval_func) {

      for (int icv = 0; icv < ncv; ++icv)
        v_ptr[icv] = 0.0;

      const double cp = R_gas*gamma/(gamma-1.0);

      double * ent    = new double[ncv_g];
      computeIgEntropy(ent);
      updateCvData(ent);


      for (int ifa = 0; ifa < nfa; ++ifa) {

        const int icv0 = cvofa[ifa][0];
        const int icv1 = cvofa[ifa][1];

        const double area     = MAG(n_fa[ifa]);
        const double delta    = area/area_over_delta_fa[ifa];
        const double dp       = p[icv1] - p[icv0];
        const double dh       = cp*(T[icv1] - T[icv0]);
        const double v_avg    = 0.5*(1.0/rho[icv0] + 1.0/rho[icv1]);
        const double beta_avg = 0.5*(1.0/(R_gas*T[icv0]) + 1.0/(R_gas*T[icv1]));
        const double dsor     = ent[icv1] - ent[icv0]; // the R has already been divided ..
        const double gr       = -dsor + beta_avg*(dh - v_avg*dp);

        const double c_avg    = fabs(gr)*area*delta/6.0;
        v_ptr[icv0]          += c_avg;
        if ( icv1 < ncv)
          v_ptr[icv1]        += c_avg;

      }

      for (int icv = 0; icv < ncv; ++icv)
        v_ptr[icv] /= vol_cv[icv];

      dumpRange(v_ptr, ncv, "gr");

      delete[] ent;

    }

    return CTI_DATA_OK;


  }
  else if ( name == "shear_flux") {

    if ( args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double (*v_ptr)[3] = createSignedFaD2Data(v);

    if ( b_eval_func) {

      FOR_IFA FOR_I3 v_ptr[ifa][i] = 0.0;

      FOR_IFA {

        const int icv0 = cvofa[ifa][0];
        const int icv1 = cvofa[ifa][1];

        const double area                   =   MAG(n_fa[ifa]);
        const double aod_half               =   0.5*area_over_delta_fa[ifa];
        const double mu_tmp                 =   0.5*(mu_lam[icv0]+mu_lam[icv1]+mu_sgs[icv0]+mu_sgs[icv1])*area;
        const double mu_total_c             =   (mu_lam[icv0]+mu_lam[icv1]+mu_sgs[icv0]+mu_sgs[icv1])*aod_half;

        double u0_cn = 0.0;
        double u1_cn = 0.0;
        double unit_n[3];
        FOR_I3 unit_n[i] = n_fa[ifa][i]/area;
        FOR_I3 {
          u0_cn += u[icv0][i]*unit_n[i];
          u1_cn += u[icv1][i]*unit_n[i];
        }
        const double one_third_dun = (u1_cn - u0_cn)/3.0;

        FOR_I3 v_ptr[ifa][i] -= mu_total_c*(u[icv1][i] - u[icv0][i] + one_third_dun*unit_n[i]);

        // viscous transpose terms...

        double dundx[3]           = {0.0, 0.0, 0.0};
        double one_third_dundxn   = 0.0;
        double two_third_Skk      = 0.0;
        FOR_K3 {
          dundx[k] = 0.5* ( (dudx[icv0][0][k] + dudx[icv1][0][k])*unit_n[0] +
              (dudx[icv0][1][k] + dudx[icv1][1][k])*unit_n[1] +
              (dudx[icv0][2][k] + dudx[icv1][2][k])*unit_n[2] );

          one_third_dundxn  += dundx[k]*unit_n[k];
          two_third_Skk += dudx[icv0][k][k] + dudx[icv1][k][k];
        }

        two_third_Skk /= 3.0;
        one_third_dundxn /= 3.0;

        FOR_I3 v_ptr[ifa][i] -= (dundx[i] - unit_n[i]*(one_third_dundxn + two_third_Skk))*mu_tmp;
      }

      dumpRange(v_ptr, nfa, "shear_flux");

    }

    return CTI_DATA_OK;

  }
  else if ( (name == "ent") || (name == "entropy")) {

    double * v_ptr = createCvD1Data(v);

    if ( b_eval_func)
      computeIgEntropy(v_ptr);

    return CTI_DATA_OK;

  }
  else if (name == "q_criterion") {

    if ( args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double * v_ptr = createCvD1Data(v);

    if ( b_eval_func) {

      double sij[3][3],  wij[3][3];
      double omega2, strte2;
      FOR_ICV {

        sij[0][0] = dudx[icv][0][0];
        sij[1][1] = dudx[icv][1][1];
        sij[2][2] = dudx[icv][2][2];

        sij[0][1] = 0.5*(dudx[icv][0][1]+dudx[icv][1][0]);
        sij[1][2] = 0.5*(dudx[icv][1][2]+dudx[icv][2][1]);
        sij[0][2] = 0.5*(dudx[icv][0][2]+dudx[icv][2][0]);

        wij[0][0] = 0;
        wij[1][1] = 0;
        wij[2][2] = 0;

        wij[0][1] = 0.5*(dudx[icv][0][1]-dudx[icv][1][0]);
        wij[1][2] = 0.5*(dudx[icv][1][2]-dudx[icv][2][1]);
        wij[0][2] = 0.5*(dudx[icv][0][2]-dudx[icv][2][0]);

        sij[1][0] = sij[0][1];
        sij[2][1] = sij[1][2];
        sij[2][0] = sij[0][2];

        wij[1][0] = -wij[0][1];
        wij[2][1] = -wij[1][2];
        wij[2][0] = -wij[0][2];

        omega2 = 0.5*( wij[0][0]*wij[0][0] + wij[0][1]*wij[0][1] + wij[0][2]*wij[0][2] +
            wij[1][0]*wij[1][0] + wij[1][1]*wij[1][1] + wij[1][2]*wij[1][2] +
            wij[2][0]*wij[2][0] + wij[2][1]*wij[2][1] + wij[2][2]*wij[2][2] );
        strte2 = 0.5*( sij[0][0]*sij[0][0] + sij[0][1]*sij[0][1] + sij[0][2]*sij[0][2] +
            sij[1][0]*sij[1][0] + sij[1][1]*sij[1][1] + sij[1][2]*sij[1][2] +
            sij[2][0]*sij[2][0] + sij[2][1]*sij[2][1] + sij[2][2]*sij[2][2] );

        v_ptr[icv] = omega2 - strte2;
      }

    }

    return CtiRegister::CTI_DATA_OK;

  }
  else if (name == "lambda2") {

    if ( args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double * v_ptr = createCvD1Data(v);

    if ( b_eval_func) {

      double sij[3][3],  wij[3][3];
      double s2w2[3][3];

      FOR_ICV {

        sij[0][0] = dudx[icv][0][0];
        sij[1][1] = dudx[icv][1][1];
        sij[2][2] = dudx[icv][2][2];

        sij[0][1] = 0.5*(dudx[icv][0][1]+dudx[icv][1][0]);
        sij[1][2] = 0.5*(dudx[icv][1][2]+dudx[icv][2][1]);
        sij[0][2] = 0.5*(dudx[icv][0][2]+dudx[icv][2][0]);

        wij[0][0] = 0;
        wij[1][1] = 0;
        wij[2][2] = 0;

        wij[0][1] = 0.5*(dudx[icv][0][1]-dudx[icv][1][0]);
        wij[1][2] = 0.5*(dudx[icv][1][2]-dudx[icv][2][1]);
        wij[0][2] = 0.5*(dudx[icv][0][2]-dudx[icv][2][0]);

        sij[1][0] = sij[0][1];
        sij[2][1] = sij[1][2];
        sij[2][0] = sij[0][2];

        wij[1][0] = -wij[0][1];
        wij[2][1] = -wij[1][2];
        wij[2][0] = -wij[0][2];

        // build symm. tensor Omega2 + S2
        FOR_I3 {
          FOR_J3 {
            s2w2[i][j] = sij[i][j]*sij[i][j] + wij[i][j]*wij[i][j];
          }
        }
        // compute eigenvalues
        double lam[3];
        double eV[3][3];
        MiscUtils::eigenDecomposition(s2w2,eV,lam);

        v_ptr[icv] = lam[1];

      }

    }

    return CtiRegister::CTI_DATA_OK;

  }

  // this solver doesnt know how to evaluate this data -- go to FlowSolver

  return FlowSolver::funcEvalCtiData(v,name,args,b_eval_func);
}

#undef FOR_BCZONE
