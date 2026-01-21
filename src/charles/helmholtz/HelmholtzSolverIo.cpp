
#include "HelmholtzSolver.hpp"

void HelmholtzSolver::initFromParams() {
  // parameter-based initialization of data...
  StaticSolver::initFromParams();

  FOR_PARAM_MATCHING("INIT_STATIONARY") {

    if (frame_rotation){

      if (mpi_rank==0){
        cout << "Found INIT_STATIONARY, frame rotation active..." << endl;
      }

      // INIT_STATIONARY <r3-var-name> <value0> <value1> <value2>
      // INIT_STATIONARY u 0.0 0.0 0.0

      int iarg = 0;
      while (iarg < param->size()) {

        const string var_name = param->getString(iarg++);

        if ( CtiRegister::CtiData* data = CtiRegister::getRegisteredCtiData(var_name) ) {

          const int type = data->getType();

          switch ( type ) {
          case DN3_DATA:
            {
              const int n = data->size();
              double val[3];
              for (int i =0; i < 3; ++i)
                val[i] = param->getDouble(iarg++);

              if (n == ncv || n == ncv_g || n == ncv_g2){
                for (int icv =0; icv < n; ++icv){
                  double r[3];
                  FOR_I3 r[i] = x_cv[icv][i] - frame_rotation[i+3];

                  data->dn3(icv,0) = val[0] - frame_rotation[1]*r[2] +
                                              frame_rotation[2]*r[1];
                  data->dn3(icv,1) = val[1] - frame_rotation[2]*r[0] +
                                              frame_rotation[0]*r[2];
                  data->dn3(icv,2) = val[2] - frame_rotation[0]*r[1] +
                                              frame_rotation[1]*r[0];
                }
              }
              else {
                CERR(" > unable to process INIT_STATIONARY command; var " << var_name << " is not control volume data");
              }
            }
            break;

          case I_DATA:
          case D_DATA:
          case IN_DATA:
          case DN_DATA:
            CERR(" > unable to process INIT_STATIONARY command; var " << var_name << " is not vector data, use INIT")
            break;

          default:

            // if you get here, we forgot to add a case.  this isn't the users
            // fault but an oversight on our part.  boo.

            assert(0);

          }

          // data has been set so update the flag -- is the value important?

          data->setFlag(1);

        } else {

          CERR( " > unable to process INIT command; var " << var_name << " is not registered data.");

        }
      }

    } //if (frame_rotation)
    else {
       COUT1(" > Warning, skipping INIT_STATIONARY, only valid if using ROTATING_REFERENCE_FRAME");
    }

  }

  if (checkParam("INIT_FROM_PF")||checkParam("INIT_FROM_POTENTIAL_FLOW")) {

    initFromPotentialFlow();

  }
  else if ( Param * param = getParam("INTERP_FROM_RESTART")) {

    // rescale interpolated data if not coming from charles
    // assumes all interps are coming from same solver (just check first)
    if (param != NULL) {
      const string filename = param->getString(0);
      if ((filename.find(".ip") != string::npos)||(filename.find(".cas") != string::npos)) {

        // fluent interpolation file
        if (mpi_rank == 0)
          cout << "modifying fluent interpolated data for consistency with charles..." << endl;

        // will recompute these (just taking velocity)...
        setDataFlag("p",0);
        setDataFlag("rho",0);
      }
    }
  }

}

void HelmholtzSolver::initFromPotentialFlow() {

  COUT1("initFromPotentialFlow()");

  // 1. solve a potential flow problem from the boundaries to build a rhou
  // field: rhou = -grad(phi)...

  double * phi = new double[ncv_g];
  for (int icv = 0; icv < ncv_g; ++icv) phi[icv] = 0.0;

  double * mdot_bf = new double[nbf];
  FOR_IBF mdot_bf[ibf] = 0.0;

  // const rho and p, include bad vals for checking...
  double rho_init = getDoubleParam("RHO",-1.0);
  double p_init = HUGE_VAL;

  double mdot_sum = 0.0;
  double outlet_rhoun = 0.0;
  FOR_IZONE(bfZoneVec) {
    const string zone_name = bfZoneVec[izone].getName();
    if (Param * param = getParam(zone_name)) {
      int iarg = 0;
      const string bc_type = param->getString(iarg++);
      if ( bc_type == "INLET") {
        // INLET has ordered params: u v w
        double u_in[3];
        u_in[0] = param->getDouble(1);
        u_in[1] = param->getDouble(2);
        u_in[2] = param->getDouble(3);
        for (int ibf = bfZoneVec[izone].ibf_f; ibf <= bfZoneVec[izone].ibf_l; ++ibf) {
          mdot_bf[ibf] = rho_init*DOT_PRODUCT(u_in,n_bf[ibf]);
        }
        mdot_sum = -rho_init*DOT_PRODUCT(u_in,bfZoneVec[izone].n_global);
        // report...
        if (mpi_rank == 0)
          cout << " > INLET \"" << zone_name << "\": MDOT=" << mdot_sum <<
            " (area=" << bfZoneVec[izone].area_global <<
            ", rhoun=" << mdot_sum/bfZoneVec[izone].area_global << ")" << endl;
      }
      else if ((bc_type == "OUTLET")||(bc_type == "OUTLET_VV")) {
        // for now, accumulate areas in outlet_rhoun...
        // OUTLET and OUTLET_VV have ordered params: sigma L p_ref Ma
        p_init = param->getDouble(3);
        outlet_rhoun += bfZoneVec[izone].area_global;
      }
    }
  }
  assert((rho_init > 0.0)&&(p_init != HUGE_VAL));

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
    const string zone_name = bfZoneVec[izone].getName();
    if (Param * param = getParam(zone_name)) {
      int iarg = 0;
      const string bc_type = param->getString(iarg++);
      if ((bc_type == "OUTLET")||(bc_type == "OUTLET_VV")) {
        if (mpi_rank == 0)
          cout << " > OUTLET \"" << zone_name << "\": computed MDOT=" <<
            outlet_rhoun*bfZoneVec[izone].area_global << " (area=" <<
            bfZoneVec[izone].area_global << ", rhoun=" << outlet_rhoun << ")" << endl;
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
  solveCvCg(phi,A,rhs,pf_zero,pf_maxiter,false);

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
  FOR_ICV FOR_I3 u[icv][i] = -u[icv][i]/rho_init; // div by rho to get u, flip to get u sign correct
  FOR_ICV rho[icv] = rho_init;
  FOR_ICV p[icv] = p_init;

  delete[] phi;
  delete[] mdot_bf;

}

void HelmholtzSolver::queryBcs() {

  // QUERY_BC <bc-name> [INTERVAL = <interval>] [WRITE]
  FOR_PARAM_MATCHING("QUERY_BC") {

    // check if the bc is matched against a known query, we will
    // use the entire param string as the key for the bc_query map

    map<string,pair<int,HelmholtzBc*> >::iterator bc_it = bc_queries.find(param->str());

    bool b_write = false;

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
        else if (token == "WRITE") {
          b_write = true;
        }
        else {
          CERR( " Invalid query_bc syntax; QUERY_BC [INTERVAL <interval>] [WRITE]");
        }
      }

      HelmholtzBc* bc = getBc(bc_name);

      if ( bc == NULL) {

        CWARN(" > unable to find boundary zone for QUERY: " << bc_name);

      }
      else {

        if (bc->ss == NULL) bc->ss = new std::stringstream();
        bc->b_write = b_write;

        pair<map<string,pair<int,HelmholtzBc*> >::iterator,bool> ret =
          bc_queries.insert(pair<string,pair<int,HelmholtzBc*> >( param->str(),
                pair<int,HelmholtzBc*>(interval,bc)));

        assert( ret.second); // ensure that the query was properly inserted.
        bc_it = ret.first;
      }
    }

    if ( bc_it != bc_queries.end() ) {

      const int query_interval = bc_it->second.first;
      HelmholtzBc* bc          = bc_it->second.second;

      if ( step%query_interval == 0) {
        bc->query(bc_it->first);
      }

    }

  }
}

CtiRegister::CtiDataError HelmholtzSolver::funcEvalCtiData(CtiRegister::CtiData& v,const string& name,
    list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if ( name == "p_total" ) {

    if ( args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double * v_ptr = createCvD1Data(v);
    if ( b_eval_func) {

      FOR_ICV {
        const double magu2 = DOT_PRODUCT(u[icv],u[icv]);
        v_ptr[icv] = p[icv]+0.5*rho[icv]*magu2;
      }

    }

    return CTI_DATA_OK;

  }
  else if ( name == "p_dynamic" ) {

    if ( args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double * v_ptr = createCvD1Data(v);
    if ( b_eval_func) {

      FOR_ICV {
        const double magu2 = DOT_PRODUCT(u[icv],u[icv]);
        v_ptr[icv] = 0.5*rho[icv]*magu2;
      }

    }

    return CTI_DATA_OK;

  }
  else if ( name == "vorticity" ) {

    if (args.size() != 0)
      return(CtiRegister::CTI_DATA_ARG_COUNT);

    double (*v_ptr)[3] = createCvD2Data(v);

    if ( b_eval_func ) {

      // alloc some ghosts for communication

      double (*arg_g)[3] = new double[ncv_g-ncv][3];
      updateCvDataSeparateGhosts(u,arg_g);

      // use transpose of cvocv_grad_coeff to compute div at cv

      FOR_ICV {
        const int coc_f = cvocv_i[icv];
        v_ptr[icv][0] = CROSS_PRODUCT_0(cvocv_grad_coeff[coc_f],u[icv]);
        v_ptr[icv][1] = CROSS_PRODUCT_1(cvocv_grad_coeff[coc_f],u[icv]);
        v_ptr[icv][2] = CROSS_PRODUCT_2(cvocv_grad_coeff[coc_f],u[icv]);

        for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          if (icv_nbr < ncv) {
            v_ptr[icv][0] += CROSS_PRODUCT_0(cvocv_grad_coeff[coc],u[icv_nbr]);
            v_ptr[icv][1] += CROSS_PRODUCT_1(cvocv_grad_coeff[coc],u[icv_nbr]);
            v_ptr[icv][2] += CROSS_PRODUCT_2(cvocv_grad_coeff[coc],u[icv_nbr]);
          }
          else {
            v_ptr[icv][0] += CROSS_PRODUCT_0(cvocv_grad_coeff[coc],arg_g[icv_nbr-ncv]);
            v_ptr[icv][1] += CROSS_PRODUCT_1(cvocv_grad_coeff[coc],arg_g[icv_nbr-ncv]);
            v_ptr[icv][2] += CROSS_PRODUCT_2(cvocv_grad_coeff[coc],arg_g[icv_nbr-ncv]);
          }
        }
      }
      delete[] arg_g;
    }

    return CtiRegister::CTI_DATA_OK;

  }
  else if ( name == "helicity" ) {

    if (args.size() != 0)
      return(CtiRegister::CTI_DATA_ARG_COUNT);

    double *v_ptr = createCvD1Data(v);

    if ( b_eval_func ) {

      // alloc some ghosts for communication

      double (*arg_g)[3] = new double[ncv_g-ncv][3];
      updateCvDataSeparateGhosts(u,arg_g);

      // use transpose of cvocv_grad_coeff to compute div at cv

      FOR_ICV {
        const int coc_f = cvocv_i[icv];
        double vorticity[3];
        vorticity[0] = CROSS_PRODUCT_0(cvocv_grad_coeff[coc_f],u[icv]);
        vorticity[1] = CROSS_PRODUCT_1(cvocv_grad_coeff[coc_f],u[icv]);
        vorticity[2] = CROSS_PRODUCT_2(cvocv_grad_coeff[coc_f],u[icv]);

        for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          if (icv_nbr < ncv) {
            vorticity[0] += CROSS_PRODUCT_0(cvocv_grad_coeff[coc],u[icv_nbr]);
            vorticity[1] += CROSS_PRODUCT_1(cvocv_grad_coeff[coc],u[icv_nbr]);
            vorticity[2] += CROSS_PRODUCT_2(cvocv_grad_coeff[coc],u[icv_nbr]);
          }
          else {
            vorticity[0] += CROSS_PRODUCT_0(cvocv_grad_coeff[coc],arg_g[icv_nbr-ncv]);
            vorticity[1] += CROSS_PRODUCT_1(cvocv_grad_coeff[coc],arg_g[icv_nbr-ncv]);
            vorticity[2] += CROSS_PRODUCT_2(cvocv_grad_coeff[coc],arg_g[icv_nbr-ncv]);
          }
        }
        v_ptr[icv] = DOT_PRODUCT(u[icv],vorticity);
      }
      delete[] arg_g;
    }

    return CtiRegister::CTI_DATA_OK;

  }
  else if (name == "q_criterion"){
    if ( args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double * v_ptr = createCvD1Data(v);
    if ( b_eval_func) {
      double sij[3][3],  wij[3][3];
      double omega2, strte2;

      double (*dudx)[3][3] = new double[ncv][3][3];
      //StaticSolver::calcCvGrad(dudx,u);
      calcCvGrad(dudx,u);

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
      delete[] dudx;
    }
    return CtiRegister::CTI_DATA_OK;

  }
  else if (name == "lambda2"){
    if ( args.size() != 0)
      return CTI_DATA_ARG_COUNT;

    double * v_ptr = createCvD1Data(v);
    if ( b_eval_func) {
      double sij[3][3],  wij[3][3];
      double s2w2[3][3];

      double (*dudx)[3][3] = new double[ncv][3][3];
      //StaticSolver::calcCvGrad(dudx,u);
      calcCvGrad(dudx,u);

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
      delete[] dudx;
    }
    return CtiRegister::CTI_DATA_OK;

  }

  return FlowSolver::funcEvalCtiData(v,name,args,b_eval_func);

}
