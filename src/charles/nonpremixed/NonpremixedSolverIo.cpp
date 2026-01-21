
#include "NonpremixedSolver.hpp"

void NonpremixedSolver::loadBalanceHook(StripedMesh* sm) {

  // for now, make this a param. TODO: get rid of this once optimized.
  // Did an experiment on excalibur, and got 17 as best value.
  const int wm_lb_cost = getIntParam("WM_LB_COST", 17);

  //TODO notes: get ris of this and do bc initialization in a way that
  //allows bcs to set their cost during their constructor...
  COUT1("NonpremixedSolver::loadBalanceHook()");
  FOR_IZONE(bfZoneVec) {
    const string zone_name = bfZoneVec[izone].getName();
    if (Param * p = getParam(zone_name)) {
      const string bc_type = p->getString(0);
      if ((bc_type == "WM_ALG_ADIABATIC") || (bc_type == "WM_ALG_ISOTHERMAL") || (bc_type == "WM_ALG_CONDUCTIVE")) {
        bfZoneVec[izone].lb_cost = wm_lb_cost;
      }
    }
  }
}//loadBalanceHook()

void NonpremixedSolver::initFromPotentialFlow(const bool coldflow = false) {

  COUT1("initFromPotentialFlow()");

  // 1. solve a potential flow problem from the boundaries to build a rhou
  // field: rhou = -grad(phi)...

  double * phi = new double[ncv_g];
  for (int icv = 0; icv < ncv_g; ++icv) phi[icv] = 0.0;

  double * mdot_bf = new double[nbf];
  FOR_IBF mdot_bf[ibf] = 0.0;

  double * Z_bf = new double[nbf];
  FOR_IBF Z_bf[ibf] = 0.0;

  double mdot_sum = 0.0;
  double outlet_rhoun = 0.0;

  FOR_IZONE(bfZoneVec) {
    const string zone_name = bfZoneVec[izone].getName();
    if (Param * param = getParam(zone_name)) {
      int iarg = 0;
      const string bc_type = param->getString(iarg++);
      if ( bc_type == "NSCBC_MTS") {
        // examples...
        // INLET NSCBC_MTS 22.23056205 678.55 0 0 Y_NOX 0
        // PM1_INLET NSCBC_MTS 0.117026831 482.04 1 0 Y_NOX 0
        const double this_mdot = param->getDouble(iarg++);
        ++iarg; // skip T
        const double this_Z = param->getDouble(iarg++);
        // looks good...
        mdot_sum += this_mdot;
        const double this_rhoun = this_mdot/bfZoneVec[izone].area_global;
        // report...
        if (mpi_rank == 0)
          cout << " > INLET (NSCBC_MTS) \"" << zone_name << "\": MDOT=" << this_mdot << " (area=" << bfZoneVec[izone].area_global << ", rhoun=" << this_rhoun << ")" << " Z=" << this_Z << endl;
        for (int ibf = bfZoneVec[izone].ibf_f; ibf <= bfZoneVec[izone].ibf_l; ++ibf) {
          mdot_bf[ibf] = -this_rhoun*area_bf[ibf]; // inlet negative by convention
          Z_bf[ibf] = this_Z;
        }
      }
      else if ( bc_type == "CBC_MPTS") {
        // examples...
        // inlet-fuel CBC_MPTS 0.0005252985 101325 300.0 1.0 0.
        const double this_mdot = param->getDouble(iarg++);
        if (this_mdot > 0.0) {
          ++iarg; // skip p
          ++iarg; // skip T
          const double this_Z = param->getDouble(iarg++);
          // looks good...
          mdot_sum += this_mdot;
          const double this_rhoun = this_mdot/bfZoneVec[izone].area_global;
          // report...
          if (mpi_rank == 0)
            cout << " > INLET (CBC_MPTS) \"" << zone_name << "\": MDOT=" << this_mdot << " (area=" << bfZoneVec[izone].area_global << ", rhoun=" << this_rhoun << ")" << " Z=" << this_Z << endl;
          for (int ibf = bfZoneVec[izone].ibf_f; ibf <= bfZoneVec[izone].ibf_l; ++ibf) {
            mdot_bf[ibf] = -this_rhoun*area_bf[ibf]; // inlet negative by convention
            Z_bf[ibf] = this_Z;
          }
        }
        else {
          // this is an outlet...
          outlet_rhoun += bfZoneVec[izone].area_global;
        }
      }
      else if (bc_type == "NSCBC_OUTLET_P") {
        // example...
        // OUTLET NSCBC_OUTLET_P 3.766e5 .1 1
        // for now, accumulate areas in outlet_rhoun...
        outlet_rhoun += bfZoneVec[izone].area_global;
      }
    }
  }

  // also add the volumetric mass source in the cone...

  double * rhs = new double[ncv];
  for (int icv = 0; icv < ncv; ++icv)
    rhs[icv] = 0.0;

  /*
  double my_mdot_cone = 0.0;
  FOR_ICV if (flag_cone[icv]) {
  const double rhoE_src = spray_alpha*max(0.0,(T[icv]-spray_T_target));
  const double rho_src = rhoE_src*12.0114/-2.53296e+07;
  rhs[icv] -= vol_cv[icv]*rho_src;
  my_mdot_cone += vol_cv[icv]*rho_src;
  }
  double mdot_cone;
  MPI_Allreduce(&my_mdot_cone,&mdot_cone,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
  if (mpi_rank == 0) {
  cout << " > cone model: MDOT=" << mdot_cone << endl;
  }
  mdot_sum += mdot_cone;
  */

  // inlets and outlets must balance in potential flow...
  assert(mdot_sum > 0.0);
  assert(outlet_rhoun > 0.0); // recall this is holding the area right now -- should be positive
  outlet_rhoun = mdot_sum/outlet_rhoun;

  // now go back through the outlet(s) and specify the mdot_bf...
  FOR_IZONE(bfZoneVec) {
    const string zone_name = bfZoneVec[izone].getName();
    if (Param * param = getParam(zone_name)) {
      int iarg = 0;
      const string bc_type = param->getString(iarg++);
      if ( bc_type == "CBC_MPTS") {
        // examples...
        // inlet-fuel CBC_MPTS 0.0005252985 101325 300.0 1.0 0.
        const double this_mdot = param->getDouble(iarg++);
        if (this_mdot < 0.0) {
          if (mpi_rank == 0)
            cout << " > OUTLET \"" << zone_name << "\": computed MDOT=" << outlet_rhoun*bfZoneVec[izone].area_global << " (area=" << bfZoneVec[izone].area_global << ", rhoun=" << outlet_rhoun << ")" << endl;
          for (int ibf = bfZoneVec[izone].ibf_f; ibf <= bfZoneVec[izone].ibf_l; ++ibf) {
            mdot_bf[ibf] = outlet_rhoun*area_bf[ibf]; // outlet positive by convention
          }
        }
      }
      else if (bc_type == "NSCBC_OUTLET_P") {
        // example...
        // OUTLET NSCBC_OUTLET_P 3.766e5 .1 1
        if (mpi_rank == 0)
          cout << " > OUTLET \"" << zone_name << "\": computed MDOT=" << outlet_rhoun*bfZoneVec[izone].area_global << " (area=" << bfZoneVec[izone].area_global << ", rhoun=" << outlet_rhoun << ")" << endl;
        for (int ibf = bfZoneVec[izone].ibf_f; ibf <= bfZoneVec[izone].ibf_l; ++ibf) {
          mdot_bf[ibf] = outlet_rhoun*area_bf[ibf]; // outlet positive by convention
        }
      }
    }
  }

  // ----------------------------------
  // solver...
  // ----------------------------------

  for (int ibf = 0; ibf < nbf; ++ibf) {
    const int icv = cvobf[ibf];
    rhs[icv] += mdot_bf[ibf];
  }

  /*
  for (int ifa = 0; ifa < nfa; ++ifa) {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    const double rhoun0 = rho[icv0]*DOT_PRODUCT(u[icv0],n_fa[ifa]); // n_fa contains area
    const double rhoun1 = rho[icv1]*DOT_PRODUCT(u[icv1],n_fa[ifa]); // n_fa contains area
    rhs[icv0] += 0.5*(rhoun0+rhoun1);
    if (icv1 < ncv)
      rhs[icv1] -= 0.5*(rhoun0+rhoun1);
  }
  */

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

  StaticSolver::calcCvGrad(u,phi);
  //calcCvGradLeastSquares(u,phi,mdot_bf); // this gradient can take the boundary normal derrivative
  FOR_ICV FOR_I3 u[icv][i] = -u[icv][i]; // flip to get u sign correct

  // ====================================================
  // now solve the scalars...
  // ====================================================

  int * sweep_order = new int[ncv];
  {
    vector<pair<double,int> > diPairVec(ncv);
    for (int icv = 0; icv < ncv; ++icv) {
      diPairVec[icv].first = -phi[icv]; // we want largest to smallest in phi, so sort based on -phi
      diPairVec[icv].second = icv;
    }
    sort(diPairVec.begin(),diPairVec.end());
    for (int icv = 0; icv < ncv; ++icv) {
      sweep_order[icv] = diPairVec[icv].second;
    }
  }

  double * mdot_fa = new double[nfa];
  for (int ifa = 0; ifa < nfa; ++ifa) {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    mdot_fa[ifa] = area_over_delta_fa[ifa]*(phi[icv0]-phi[icv1]); // note sign here: decreasing potential = positive mass
  }

  const double pf_scalar_zero = getDoubleParam("PF_SCALAR_ZERO",1.0E-12);
  const int pf_scalar_maxiter = getIntParam("PF_MAXITER",2000);

  // --------------------
  // Z...
  // --------------------

  for (int icv = 0; icv < ncv_g; ++icv) Z[icv] = 0.0;

  int iter = 0;
  int done = 0;
  while (done == 0) {

    ++iter;

    double my_dZ_max = 0.0;
    for (int ii = 0; ii < ncv; ++ii) {

      const int icv = sweep_order[ii];
      double sum_mdot[2] = { 0.0, 0.0 };

      for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
        const int ibf = bfocv_v[boc];
        if (mdot_bf[ibf] <= 0.0) {
          // negative mdot_bf means INFLOW to icv...
          sum_mdot[0] += mdot_bf[ibf]*Z_bf[ibf];
          sum_mdot[1] += mdot_bf[ibf];
        }
      }

      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        if ((mdot_fa[ifa] < 0.0)&&(cvofa[ifa][0] == icv)) {
          // this flux is leaving icv1 and going to icv0...
          sum_mdot[0] += mdot_fa[ifa]*Z[cvofa[ifa][1]];
          sum_mdot[1] += mdot_fa[ifa];
        }
        else if ((mdot_fa[ifa] > 0.0)&&(cvofa[ifa][1] == icv)) {
          // this flux is leaving icv0 and going to icv1...
          sum_mdot[0] -= mdot_fa[ifa]*Z[cvofa[ifa][0]];
          sum_mdot[1] -= mdot_fa[ifa];
        }
      }

      if (sum_mdot[1] != 0.0) {
        const double Z_new = sum_mdot[0]/sum_mdot[1];
        my_dZ_max = max(my_dZ_max,fabs(Z[icv]-Z_new));
        Z[icv] = Z_new;
      }

    }

    updateCvData(Z);

    double dZ_max;
    MPI_Reduce(&my_dZ_max,&dZ_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    if (mpi_rank == 0) {
      cout << " > scalar Z solver iter: " << iter << " dZ_max: " << dZ_max << endl;
      if (dZ_max <= pf_scalar_zero)
        done = 1;
      else if (iter >= pf_scalar_maxiter) {
        if (mpi_rank == 0) {
          cout << " > Warning: scalar Z did not converge. Consider increasing PF_SCALAR_MAXITER." << endl;
        }
        done = 2;
      }
    }
    MPI_Bcast(&done,1,MPI_INT,0,mpi_comm);

  }

  delete[] phi;
  delete[] sweep_order;
  delete[] mdot_fa;
  delete[] mdot_bf;
  delete[] Z_bf;

  // ------------------------------------
  // C...
  // ------------------------------------

  assert( chemtable != NULL);
  if (coldflow) {
    COUT1(" > setting coldflow condition...");
    FOR_ICV C[icv] = 0.0;
    FOR_ICV Zvar[icv] = 0.0;
  } else {
    COUT1(" > setting reacting condition...");
    //TODO: Figure out how to lookup from top
    chemtable->lookupReduced(C, "prog", Z, ncv);
    FOR_ICV Zvar[icv] = 0.0;
  }

  // now set other properties...
  // XXX Compute Zvar ?? or at least set it to zero
  chemtable->lookup(R,"R",T_p,"T",e_p,"e",gamma,"gamma",a_gam,"a_gamma",Z,Zvar,C,ncv);

  // step 1, set the pressure. Need an anchoring strategy but for now use table p
  // as p total and reduce using dynamic pressure...

  bool b_u_max = false;
  double pf_u_max;
  if (Param * param = getParam("PF_U_MAX")) {
    b_u_max = true;
    pf_u_max = param->getDouble();
    if (mpi_rank == 0)
      cout << " > velocity magnitude will be limited to PF_U_MAX=" << pf_u_max << endl;
  }

  int8 my_count = 0;
  bool first = true;
  FOR_ICV {

    const double a = 2.0*R[icv]*T_p[icv];
    const double b = -2.0*chemtable->pressure;
    const double c = DOT_PRODUCT(u[icv],u[icv]); // recall u contains rhou now

    double disc = b*b - 4.0*a*c;
    if (disc < 0.0) {
      if (first) {
        cout << "Error: x_cv: " << COUT_VEC(x_cv[icv]) << " rhou: " << COUT_VEC(u[icv]) <<
          " Z,C: " << Z[icv] << " " << C[icv] <<
          " R: " << R[icv] << " T: " << T_p[icv] << " chemtable->pressure: " << chemtable->pressure << endl;
        first  = false;
      }
      ++my_count;
      disc = 0.0;
    }

    rho[icv] = (-b + sqrt(disc))/(2.0*a);

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
    rhoE[icv] = rho[icv]*e_p[icv] + 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);

    // Z,C already set

  }//FOR_ICV


  int8 count;
  MPI_Allreduce(&my_count,&count,1,MPI_INT8,MPI_SUM,mpi_comm);
  if (count != 0) {
    if (mpi_rank == 0)
      cout << "WARNING: PF state computation failed in " << count << " cvs" << endl;
  }

  dumpRange(rho, ncv, "rho");
  dumpRange(u, ncv, "u");
  dumpRange(rhoE, ncv, "rhoE");
  dumpRange(Z, ncv, "Z");
  dumpRange(Zvar, ncv, "Zvar");
  dumpRange(C, ncv, "C");

  updateConservativeAndPrimitiveData();
  if (mpi_rank == 0) cout << "post init_from_pf!" << endl;
  MPI_Barrier(mpi_comm);

}//initFromPotentialFlow()


void NonpremixedSolver::initFromParams() {

  // parameter-based initialization of data...

  FlowSolver::initFromParams();

  if (Param * param = getParam("INIT_US_AT_TABLE")) {

    // initialize the field based nominal thermodynamic state
    // that is present in the table ...

    const double u_init     = param->getDouble(0);
    const double v_init     = param->getDouble(1);
    const double w_init     = param->getDouble(2);
    const double Z_init     = param->getDouble(3);
    //const double Zvar_init  = param->getDouble(4);
    const double Zvar_init  = 0.0;
    const double C_init     = param->getDouble(4);

    double R_init, T_p_init, e_init;
    assert(chemtable != NULL);
    chemtable->lookup(&R_init,"R",&T_p_init,"T",&e_init, "e", &Z_init, &Zvar_init, &C_init,1);

    const double rho_init = chemtable->pressure/R_init/T_p_init;

    if (mpi_rank == 0) {
      cout << " initialization at the following conditions: "
        << "  (u,v,w) : (" << u_init << ",  " << v_init << ", " << w_init << ")" << endl;
      cout << "  p       : " << chemtable->pressure << endl;
      cout << "  T       : " << T_p_init << endl;
      cout << "  rho     : " << rho_init << endl;
      cout << "  Z       : " << Z_init << endl;
      cout << "  Zvar    : " << Zvar_init << endl;
      cout << "  C       : " << C_init << endl;
    }

    for (int icv = 0; icv < ncv; ++icv) {

      rho[icv] = rho_init;
      u[icv][0] = u_init;
      u[icv][1] = v_init;
      u[icv][2] = w_init;
      rhoE[icv] = rho_init*e_init + 0.5*rho_init*DOT_PRODUCT(u[icv],u[icv]);
      Z[icv]    = Z_init;
      Zvar[icv] = Zvar_init;
      C[icv]    = C_init;

    }

    updateConservativeAndPrimitiveData();
    RestartHashUtilities::slesSetHashFromParam(param);
  }
  else if (Param * param = (checkParam("INIT_FROM_PF") ? getParam("INIT_FROM_PF") : getParam("INIT_FROM_POTENTIAL_FLOW"))) {
    if (param->size() == 0)
      initFromPotentialFlow();
    else
      initFromPotentialFlow((param->getString(0) == "COLD"));
  }//INIT_FROM_PF
  //TODO: Add interp from *.ip, should be same as p
}

void NonpremixedSolver::queryBcs() {

  //====================================================
  // syntax:
  // QUERY_BC <bc-name> [INTERVAL = 1]
  //===================================================

  FOR_PARAM_MATCHING("QUERY_BC") {

    // check if the bc is matched against a known query, we will
    // use the entire param string as the key for the bc_query map

    map<string,pair<int,NonpremixedBc*> >::iterator bc_it = bc_queries.find(param->str());

    if (bc_it == bc_queries.end()) {

      // this is a new query-- that has not been parsed yet.

      int interval         = check_interval; // default interval is check_interval
      const string bc_name = param->getString(0);

      int iarg = 1;
      while (iarg < param->size()) {
        string token = param->getString(iarg++);
        if (token == "INTERVAL")
          interval = param->getInt(iarg++);
          if (interval <= 0) {
            CWARN(" > QUERY_BC INTERVAL expects a positive integer; setting to CHECK_INTERVAL");
            // invalid value entered, so treat as unset (default value)
            interval = check_interval;
          }
      }

      NonpremixedBc* bc = getBc(bc_name);

      if (bc == NULL) {
        CWARN(" > unable to find boundary zone for QUERY: " << bc_name);

      } else {

        pair<map<string,pair<int,NonpremixedBc*> >::iterator,bool> ret =
          bc_queries.insert(pair<string,pair<int,NonpremixedBc*> >(param->str(),
                pair<int,NonpremixedBc*>(interval,bc)));

        assert(ret.second); // ensure that the query was properly inserted.
        bc_it = ret.first;

      }
    }

    else {

      const int query_interval = bc_it->second.first;
      NonpremixedBc* bc        = bc_it->second.second;

      if (step%query_interval == 0)
        bc->query(bc_it->first);

    }

  }
}

void NonpremixedSolver::processIgnition() {

  vector<SimpleSphere> pilots;

  FOR_PARAM_MATCHING("IGNITE") {

    if ((mpi_rank == 0) && (step % check_interval == 0))
      cout << "Found an ignition zone\n";

    double tstart = -1.0;
    double tend   = -1.0;

    // revisit error parsing and the forms of ignition...

    //===================================================================================
    // syntax:
    // IGNITE GEOM SPHERE POINT <x0> <x1> <x2> RADIUS <r0> [TSTART <t1> TEND <t2>]
    // if tstart is unspecified it assumes that the ignition should be applied now
    // if tend is unspecified it assumes that the ignition should continue to be applied...
    //===================================================================================

    double xp[3] = {0.0, 0.0, 0.0};
    double rp    = -1.0;
    int iarg     = 0;

    while (iarg < param->size()) {

      string token = param -> getString(iarg++);

      if (token == "GEOM") {
        token = param->getString(iarg++);
        if (token == "SPHERE") {
          token = param->getString(iarg++);
          assert(token == "POINT");
          FOR_I3 xp[i] = param->getDouble(iarg++);
          token = param->getString(iarg++);
          assert(token == "RADIUS");
          rp = param->getDouble(iarg++);
        }
        else {
          CERR("unrecognized ignition geometry\n");
        }
      }

      else if (token == "TSTART") {
        tstart = param->getDouble(iarg++);
      }

      else if (token == "TEND") {
        tend = param->getDouble(iarg++);
      }

      else {
        CERR(" > unrecognized ignition token: " << token);
      }

    }

    if ((mpi_rank == 0) && (step % check_interval == 0))
      cout << "xp = " << xp[0] << ", " << xp[1] << ", " << xp[2] << ", rp = " << rp << "\n";

    if (rp > 0.0) {

      if ((((tstart >= 0.0) && (time > tstart)) || (tstart < 0.0)) &&
          (((tend >= 0.0) && (time < tend))       || (tend   < 0.0))) {
        pilots.push_back(SimpleSphere(xp,rp));
      }

    }

    else {

      CERR("> error in the specification of the ignition geom; e.g. \n \
          IGNITE GEOM SPHERE POINT <x0> <x1> <x2> RADIUS <r> [TSTART <t1> TEND <t2>]");

    }
  }

  if (pilots.size() > 0) {

    double * Cmax_cv = new double[ncv];
    chemtable->lookupReduced(Cmax_cv, "prog", Z, ncv);

    for (int icv = 0; icv <  ncv; ++icv) {

      // loop through the available pilots, assuming that they are small in number
      // otw, we will need to instantiate a tree of some kind...

      for(vector<SimpleSphere>::iterator it = pilots.begin(); it != pilots.end(); ++it) {

        if ( it->pointIsInside(x_cv[icv])) {
          C[icv] = Cmax_cv[icv];
          break;

        }
      }
    }

    delete[] Cmax_cv;

    updateConservativeAndPrimitiveData();

  }
}

CtiRegister::CtiDataError NonpremixedSolver::funcEvalCtiData(CtiRegister::CtiData& v,const string& name,
    list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if (name == "mass_flux") {

    // need to evaluate a flux probe of a given variable.  this is an approximation
    // of the inviscid flux (no diffusion and no considerations of convective stab)

    if (args.size() > 1)
      return CTI_DATA_ARG_COUNT;

    // we need to compute the mass flux at the faces.  we need to check to see
    // if its already available.  we could call the getUnregisteredCtiData but
    // we can circumvent the operator parsing since we know exactly what kind of
    // data we are looking for (this also bypasses all of the error throwing
    // that the getUnregisteredCtiData does)

    map<const string,CtiData>::iterator iter = currentDataMap.find("mdot_fa");
    CtiData* mf_data                         = NULL;
    double * mf                              = NULL;

    if (iter != currentDataMap.end())
      mf_data = &(iter->second);

    if (!mf_data) {

      // the data does not exist yet, so we need to populate it..
      pair<map<const string,CtiData>::iterator,bool> return_pair =
        CtiRegister::currentDataMap.insert(
            pair<const string,CtiData>("mdot_fa",CtiData()));

      assert(return_pair.second);

      mf_data = &(return_pair.first->second);
      mf      = createSignedFaD1Data(*mf_data);
      assert(mf);

      if (b_eval_func) {

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

    }

    else {

      assert(mf_data->getType() == DN_DATA);
      assert(mf_data->getTopology() == SIGNED_FA_DATA);
      mf = mf_data->getDNptr();
      assert(mf);

    }

    if (args.size() == 0) {

      // just return the mass flux

      double * v_ptr = createSignedFaD1Data(v);

      if (b_eval_func) {
        for (int ifa = 0; ifa < nfa; ++ifa)
          v_ptr[ifa] = mf[ifa];
      }

      return CTI_DATA_OK;

    } else {

      list<CtiData>::iterator arg = args.begin();
      const int datatype          = arg->getType();

      if (datatype == DN_DATA) {
        if (arg->getTopology() != CV_DATA)
          return CTI_DATA_NOT_VALID;

        double *v_ptr = createSignedFaD1Data(v);

        if (b_eval_func) {

          // for interprocessor/periodic boundaries we add our half of the flux
          // and start the parallel reduction

          double * arg_ptr = arg->getDNptr();
          for (int ifa = nfa_i; ifa < nfa; ++ifa) {
            const int icv0 = cvofa[ifa][0];
            v_ptr[ifa]      = mf[ifa]* 0.5* arg_ptr[icv0];
          }

          // note -- we need to subtract the data because the normal
          // of the face is different sign on the other rank

          updateFaDataStart(v_ptr, SUBTRACT_DATA);

          // internal faces-- no ghost data required

          for (int ifa = 0; ifa < nfa_i; ++ifa) {
            const int icv0 = cvofa[ifa][0];
            const int icv1 = cvofa[ifa][1];
            v_ptr[ifa]      = mf[ifa]* 0.5*(arg_ptr[icv0] + arg_ptr[icv1]);

          }

          updateFaDataFinish(v_ptr, SUBTRACT_DATA);

        }

        return CTI_DATA_OK;

      } else if (datatype == DN3_DATA) {

        double (*v_ptr)[3] = createFaD2Data(v);

        if (b_eval_func) {

          // for interprocessor/periodic boundaries we add our half of the flux
          // and start the parallel reduction

          double (*arg_ptr)[3] = arg->getDN3ptr();
          for (int ifa = nfa_i; ifa < nfa; ++ifa) {
            const int icv0 = cvofa[ifa][0];
            for (int i = 0; i < 3; ++i)
              v_ptr[ifa][i] = mf[ifa]* 0.5* arg_ptr[icv0][i];
          }

          updateFaDataStart(v_ptr, SUBTRACT_ROTATE_DATA);

          // internal faces-- no ghost data required

          for (int ifa = 0; ifa < nfa_i; ++ifa) {
            const int icv0 = cvofa[ifa][0];
            const int icv1 = cvofa[ifa][1];
            for (int i =0; i < 3; ++i)
              v_ptr[ifa][i] = mf[ifa]* 0.5*(arg_ptr[icv0][i] + arg_ptr[icv1][i]);
          }

          updateFaDataFinish(v_ptr, SUBTRACT_ROTATE_DATA);

        }

        return CTI_DATA_OK;

      } else {

        return CTI_DATA_NOT_VALID;
      }
    }

  }

  else if (name == "csrc") {

    // csrc is registered natively in the non-premixed solver.  but to
    // maintain syntax consistency with the premixed solver (that required
    // a func eval to reduce the face integrated source terms), we leave a
    // dummy copy here as well.

    if (args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double *v_ptr = createCvD1Data(v);
    if (b_eval_func)
      for (int icv = 0; icv < ncv; ++icv)
        v_ptr[icv] = csrc[icv];

    return CTI_DATA_OK;

  }

  else if (name == "fgr") {

    if (args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double *v_ptr = createCvD1Data(v);

    if (b_eval_func)
      for (int icv = 0; icv < ncv; ++icv)
        v_ptr[icv] = cv_compact[icv].fgr;

    return CTI_DATA_OK;

  }

  else if (name == "de") {

    if (args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double *v_ptr = createCvD1Data(v);
    if (b_eval_func) {

      for (int icv = 0; icv < ncv; ++icv) {
        const double e_cv = (rhoE[icv] - 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]))/rho[icv];
        v_ptr[icv] = e_cv - e_p[icv];
      }
    }

    return CTI_DATA_OK;

  }

  else if (name == "dT") {

    if (args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double *v_ptr = createCvD1Data(v);
    if (b_eval_func)
      for (int icv = 0; icv < ncv; ++icv)
        v_ptr[icv] = T[icv] - T_p[icv];

    return CTI_DATA_OK;

  }

  else if (name == "heat_release") {

    CWARN("There currently isn't a heat release variable stored in the chemtable!\n");
    return CTI_DATA_NOT_VALID;


  }

  else if (name == "mach") {

    if (args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double *v_ptr = createCvD1Data(v);
    if (b_eval_func)
      for (int icv = 0; icv < ncv; ++icv)
        v_ptr[icv] = MAG(u[icv])/sos[icv];

    return CTI_DATA_OK;

  }

  else if (name == "p_total") {

    if (args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double *v_ptr = createCvD1Data(v);
    if (b_eval_func) {
      for (int icv = 0; icv < ncv; ++icv) {

        const double gam = gamma[icv] + a_gam[icv]*(T[icv]-T_p[icv]);
        const double Ma2 = DOT_PRODUCT(u[icv],u[icv])/(sos[icv]*sos[icv]);
        v_ptr[icv]       = p[icv]*pow(1.0 + 0.5*(gam-1.0)*Ma2, gam/(gam-1.0));

      }
    }

    return CTI_DATA_OK;

  } else if ( name == "T_total" ) {

    if ( args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double *v_ptr = createCvD1Data(v);
    if ( b_eval_func) {

      FOR_ICV {
        const double gam = gamma[icv] + a_gam[icv]*(T[icv]-T_p[icv]);
        const double Ma2 = DOT_PRODUCT(u[icv],u[icv])/(sos[icv]*sos[icv]);
        v_ptr[icv] = T[icv]*(1.0 + 0.5*(gam-1.0)*Ma2);
      }

    }

    return CTI_DATA_OK;


  } else {

    // the NonpremixedSolver also can parse "functions" of the form
    // chemtable:<chemtable-var>() -- e.g., chemtable:Y_CO(); split the
    // function string to see if this is what was requested..

    size_t colon_pos = name.find(":");
    if (colon_pos != string::npos) {
      if (name.substr(0,colon_pos) == "chemtable") {

        string var_name  = name.substr(colon_pos+1,name.size()-colon_pos);

        if (args.size() != 0)
          return CTI_DATA_ARG_COUNT;

        // check to ensure that the chemtable possesses the requested var

        if (chemtable->varExists(var_name)) {
          double *v_ptr = createCvD1Data(v);
          if (b_eval_func)
            chemtable->lookup(v_ptr, var_name, Z, Zvar, C, ncv);

          return CTI_DATA_OK;

        } else {

          CWARN(" > unable to find " << var_name << "  in the chemtable ");
          return CTI_DATA_NOT_VALID;

        }

      } else {

        return FlowSolver::funcEvalCtiData(v,name,args,b_eval_func);

      }
    }
  }

  return FlowSolver::funcEvalCtiData(v,name,args,b_eval_func);
}
