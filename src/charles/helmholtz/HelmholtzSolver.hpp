#ifndef HELMHOLTZSOLVER_HPP
#define HELMHOLTZSOLVER_HPP

#include "FlowSolver.hpp"
#include "HelmholtzSolverBcs.hpp"
#include "CommContainers.hpp"
#include "it/ZoneFilter.hpp"

// convenience macro to loop over the boundary conditions...
#define FOR_BCZONE for(vector<HelmholtzBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)

class HelmholtzSolver : public FlowSolver {
public:

  bool b_helmholtzCorrection;

  double *rho;
  double *rho0;
  double *rho00;
  double *rho_uncorrected;
  double (*u)[3];
  double (*u0)[3];
  double (*u00)[3];
  double *p;
  double (*dpdx)[3];

  double *mu_lam;
  double *mu_sgs;
  string sgs_model;

  double *a_sgs;
  double *a_sgs_sponge;

  // face-based mass flows...

  double *mf;
  double *mf0;
  double *mf00;
  int mf_init;

  // outer iterations to converge the fully implicit time advancement

  int n_outer_iters;

  // currently, assuming a single global speed of sound

  double helmholtz_sos;

  // boundary conditions

  vector<HelmholtzBc*> bcs;
  map<string,HelmholtzBc*> bc_map;
  map<string,pair<int,HelmholtzBc*> > bc_queries;

  // momentum source container...

  double (*rhou_source)[3];
  double (*A_diag_source)[3];
  bool b_bodyForceHook;

  // multigrid stuff

  int ncg_u; // number of multi-grid levels for solveUMg
  int ncg_p; // number of multi-grid levels for solveHelmholtzMg
  double agglomeration_factor; // nominal number of fine cells in a coarse cell
  bool split_orphaned_colors; // colors can cross boundaries
  vector<CoarseGrid *> cgs; // coarse-grids
  vector<HelmholtzSolver *> css; // coarse-solvers
  vector<int> pressureSmootherVec; // 0 jacobi, 1 symmetric gauss-seidel, 2 patr, 3 conjugate gradient
  vector<int> pressureNsmoothVec;
  bool b_solve_p_on_coarsest;
  vector<int> momentumSmootherVec; // 0 jacobi, 1 gauss-seidel, 2 patr, 3 triangular iterative method
  vector<int> momentumNsmoothVec;

  // vida vor's volume approximation
  // planar face assumption at boundary cells...
  double *planar_vol_ratio_cv;

  // for use with a rotating reference frame ..
  double * frame_rotation;

  // use to write out a serial state...
  int serial_state_interval;
  string serial_state_prefix;

  HelmholtzSolver() {

    b_helmholtzCorrection = checkParam("HELMHOLTZ_CORRECTION");

    // register data with i/o

    rho   = NULL; registerCvData(rho,"rho",READWRITE_DATA);
    rho0  = NULL; registerCvData(rho0,"rho0",READWRITE_DATA);
    rho00 = NULL; registerCvData(rho00,"rho00",READWRITE_DATA);
    u     = NULL; registerCvData(u,"u",READWRITE_DATA);
    u0    = NULL; registerCvData(u0,"u0",READWRITE_DATA);
    u00   = NULL; registerCvData(u00,"u00",READWRITE_DATA);
    p     = NULL; registerCvData(p,"p",READWRITE_DATA);
    dpdx  = NULL; registerCvData(dpdx,"dpdx" ,CAN_WRITE_DATA);

    mf_init = 0;
    mf    = NULL; registerSignedFaData(mf,"mf",READWRITE_DATA);
    mf0   = NULL;
    mf00  = NULL;
    if (!b_helmholtzCorrection){
      registerSignedFaData(mf0,"mf0",READWRITE_DATA);
      registerSignedFaData(mf00,"mf00",READWRITE_DATA);
    }

    // data registration without i/o
    mu_lam  = NULL; registerCvData(mu_lam , "mu_lam" , CAN_WRITE_DATA);
    mu_sgs  = NULL; registerCvData(mu_sgs , "mu_sgs" , CAN_WRITE_DATA);
    a_sgs   = NULL; registerCvData(a_sgs  ,  "a_sgs" , CAN_WRITE_DATA);
    a_sgs_sponge   = NULL; registerCvData(a_sgs_sponge  ,  "a_sgs_sponge" , CAN_WRITE_DATA);
    rho_uncorrected = NULL; registerCvData(rho_uncorrected , "rho_uncorrected" , CAN_WRITE_DATA);

    //rhs = NULL;
    rhou_source = NULL;
    A_diag_source = NULL;
    b_bodyForceHook = false;

    planar_vol_ratio_cv = NULL;

    agglomeration_factor = getDoubleParam("MG_AGGLOMERATION_FACTOR",8.0); // 2^3 = 8
    split_orphaned_colors = getBoolParam("MG_SPLIT_ORPHANED_COLORS",true);

    ncg_p = getIntParam("PRESSURE_MG_NLEVEL",1)-1;
    string tmp = getStringParam("PRESSURE_MG_SMOOTHER","3"); // 0 Jacobi, 1 GS, 2 Patr, 3 CG
    MiscUtils::splitCsv(pressureSmootherVec,tmp);
    while (pressureSmootherVec.size() <= ncg_p)
      pressureSmootherVec.push_back(pressureSmootherVec.back());
    tmp = getStringParam("PRESSURE_MG_NSMOOTH","10");
    MiscUtils::splitCsv(pressureNsmoothVec,tmp);
    while (pressureNsmoothVec.size() <= ncg_p)
      pressureNsmoothVec.push_back(pressureNsmoothVec.back());
    b_solve_p_on_coarsest = checkParam("PRESSURE_MG_SOLVE_ON_COARSEST"); // pressureNsmoothVec[ncg_p] holds maxiter

    ncg_u = getIntParam("MOMENTUM_MG_NLEVEL",1)-1;
    tmp = getStringParam("MOMENTUM_MG_SMOOTHER","2"); // 0 Jacobi, 1 GS, 2 Patr, 3 Tim
    MiscUtils::splitCsv(momentumSmootherVec,tmp);
    while (momentumSmootherVec.size() <= ncg_u)
      momentumSmootherVec.push_back(momentumSmootherVec.back());
    tmp = getStringParam("MOMENTUM_MG_NSMOOTH","10");
    MiscUtils::splitCsv(momentumNsmoothVec,tmp);
    while (momentumNsmoothVec.size() <= ncg_u)
      momentumNsmoothVec.push_back(momentumNsmoothVec.back());

    n_outer_iters = getIntParam("N_OUTER_ITERS",1);
    helmholtz_sos = getDoubleParam("HELMHOLTZ_SOS");
    parseSgsModel(sgs_model);

    frame_rotation = NULL;

    if (Param * param = getParam("ROTATING_REFERENCE_FRAME")) {

      // allocated so logic can use NULL checks below
      if (param->size()<7){
        CERR("Too few parameters for ROTATING_REFERENCE_FRAME, correct syntax\n" <<
             "  ROTATING_REFERENCE_FRAME <x> <y> <z> <w-x> <w-y> <w-z> <rot (rev/s)>");
      }

      frame_rotation = new double[6];
      frame_rotation[3] = param->getDouble(0);
      frame_rotation[4] = param->getDouble(1);
      frame_rotation[5] = param->getDouble(2);
      double omega_hat[3];
      FOR_I3 omega_hat[i] = param->getDouble(3+i);
      double omega_mag =  param->getDouble(6)*MAG(omega_hat);
      FOR_I3 frame_rotation[i] = omega_hat[i]*omega_mag;

      if (mpi_rank == 0){
        cout << " > FRAME_ROTATION active with the following rev/time: " << frame_rotation[0] << " " << frame_rotation[1] << " " << frame_rotation[2] << endl;
        cout << " > FRAME_ROTATION active with center of rotation:     " << frame_rotation[3] << " " << frame_rotation[4] << " " << frame_rotation[5] << endl;
      }

      // convert to rads/time unit...
      FOR_I3 frame_rotation[i] *= 2.0*M_PI;

    }

    serial_state_interval = -1;
    if (Param * param = getParam("WRITE_SERIAL_STATE")) {

      int iarg = 0;
      while (iarg < param->size()) {
        string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (token == "NAME") {
          serial_state_prefix = param->getString(iarg++);
        }
        else if (token == "INTERVAL") {
          serial_state_interval = param->getInt(iarg++);
        }
        else {
          if (mpi_rank == 0) cout << "Warning: skipping unrecognized WRITE_SERIAL_STATE token \"" << token << "\"" << endl;
        }
      }

      if ((serial_state_prefix == "")||(serial_state_interval <= 0)) {
        CERR("incorrect WRITE_SERIAL_STATE syntax, should be of the form\n" <<
             "  WRITE_SERIAL_STATE NAME <string> INTERVAL <int>");
      }

    }

  }

  HelmholtzSolver(const int icg);

  ~HelmholtzSolver();

  void init() {

    wtime = MPI_Wtime();

    StaticSolver::init(INIT_COMPACT_FACES|INIT_CV_GRAD);
    initVVOps();

    if (mpi_rank == 0)
      logger->setKillFilename("killcharles");

  }

  void registerBoundaryConditions() {
    assert( bcs.size() == 0);

    StaticSolver::registerBoundaryConditions(); // registers bf geometric data

    int nerr = 0;
    vector<pair<string,string> > errors;

    FOR_IZONE(bfZoneVec) {
      const string zone_name = bfZoneVec[izone].getName();
      if ( Param* p = getParam(zone_name)) {
        const string bc_type = p->getString(0);
        if ( (bc_type == "SLIP")||(bc_type == "SYMMETRY")) {
          bcs.push_back(new SlipWallHBc(&bfZoneVec[izone],this));
        } else if ( bc_type == "WALL") {
          bcs.push_back(new WallHBc(&bfZoneVec[izone],this));
        } else if ( bc_type == "WM_ALG_WALL" ) {
          bcs.push_back(new AlgebraicWallModelHBc(&bfZoneVec[izone],this));
        } else if ( bc_type == "WM_SLIP_WALL") {
          bcs.push_back(new SlipWallModelHBc(&bfZoneVec[izone],this));
        } else if ( bc_type == "INLET") {
          bcs.push_back(new InletHBc(&bfZoneVec[izone],this));
        } else if ( bc_type == "INLET_PROFILE") {
          bcs.push_back(new InletHBcProfile(&bfZoneVec[izone],this));
        } else if ( bc_type == "INFLOW_TURB") {
          bcs.push_back(new InflowTurbulenceHBc(&bfZoneVec[izone],this));
        } else if ( bc_type == "OUTLET") {
          bcs.push_back(new OutletHBc(&bfZoneVec[izone],this));
        } else if ( bc_type == "OUTLET_VV") {
          bcs.push_back(new OutletVVHBc(&bfZoneVec[izone],this));
        } else if ( bc_type == "SPONGE") {
          bcs.push_back(new SpongeHBc(&bfZoneVec[izone],this));
        } else if ( bc_type == "HOOK") {
          bcs.push_back( initHookBc(&bfZoneVec[izone]));
        } else {
          nerr++;
          errors.push_back(pair<string,string>(zone_name,bc_type));
        }
      } else {
        nerr++;
        errors.push_back(pair<string,string>(zone_name,""));
      }
    }

    // ensure that the error handling is synced..
    // should be un-necessary and we'll check below

    MPI_Bcast(&nerr,1,MPI_INT,0,mpi_comm);
    assert( nerr == int(errors.size()));
    reportBcErrors(errors);
  }

  void initBoundaryConditions() {

    // allow the boundary conditions to allocate their necessary data (not allowed to
    // set anythign here, because it may be coming from the restart data yet to be read)..

    FOR_BCZONE (*it)->initData();

    // also create a map based on the boundary conditions. provides easy access to
    // grab the bc object by name without looping through (O(N_boundary conditions))

    FOR_BCZONE bc_map[(*it)->getName()] = *it;
  }

  virtual HelmholtzBc* initHookBc(BfZone* p) {
    CERR(" > user must supply a hook boundary condition for zone: " << p->getName());
  }

  virtual HelmholtzBc* initHookBcMg(BfZone* p, const int icg) {
    CERR(" > user must supply a hook boundary condition for zone: " << p->getName() <<
         ". Override getHelmholtzSolverMg() to use multigrid solver.");
  }

  virtual HelmholtzSolver* getHelmholtzSolverMg(const int icg);

  HelmholtzBc* getBc(const string& name) const {
    map<string,HelmholtzBc*>::const_iterator it = bc_map.find(name);
    if ( it == bc_map.end()) {
      return NULL;
    } else {
      return it->second;
    }
  }

  void initData();

  void initFromParams();

  void initVVOps();

  void initFromPotentialFlow();

  void initFromCoarseGrid(CoarseGrid* cg);

  void initComplete() {

    // update the transport properties, set the (entropic) reference density
    // DAP: no longer computed each timestep, mu_lam and rho_uncorrected set once here.
    updateBaseState();

    if ( checkParam("RESET_ALL")) {

      setDataFlag("rho",0);
      setDataFlag("rho0",0);
      setDataFlag("rho00",0);
      setDataFlag("u0",0);
      setDataFlag("u00",0);
      setDataFlag("p",0);

    }

    // fill missing time levels if the data is not present or is reset

    if (!checkDataFlag("rho")) {
      FOR_ICV_G rho[icv] = rho_uncorrected[icv];
    }
    else {
      updateCvData(rho);
    }

    if (!checkDataFlag("rho0"))
      FOR_ICV_G rho0[icv] = rho[icv];
    if (!checkDataFlag("rho00"))
      FOR_ICV_G rho00[icv] = rho0[icv];

    updateCvData(u,REPLACE_ROTATE_DATA);
    if (!checkDataFlag("u0"))
      FOR_ICV_G FOR_I3 u0[icv][i] = u[icv][i];
    if (!checkDataFlag("u00"))
      FOR_ICV_G FOR_I3 u00[icv][i] = u[icv][i];

    if ( !checkDataFlag("p") ) {
      FOR_ICV_G p[icv] = 0.0;
    }
    else {
      updateCvData(p);
    }

    initASgs();
    updateSgs();

    calcPressureGrad(dpdx,p); // no need to update ghost
    dumpRange(dpdx,ncv,"DPDX");

    // face-based mass flows...

    //if any mf data is missing, extrapolation will be switched off
    //for two steps to populate the field with a "real" time history
    mf_init = 0;
    if (!checkDataFlag("mf")) {
      for (int ifa=0; ifa<nfa; ++ifa) {
        const int icv0 = cvofa[ifa][0];
        const int icv1 = cvofa[ifa][1];
        mf[ifa] = 0.5*(n_fa[ifa][0]*(rho[icv0]*u[icv0][0]+rho[icv1]*u[icv1][0]) +
                       n_fa[ifa][1]*(rho[icv0]*u[icv0][1]+rho[icv1]*u[icv1][1]) +
                       n_fa[ifa][2]*(rho[icv0]*u[icv0][2]+rho[icv1]*u[icv1][2]) );
      }
      mf_init = 2;
    }
    if (!b_helmholtzCorrection){
      if (!checkDataFlag("mf0")) {
        for (int ifa=0; ifa<nfa; ++ifa) {
          mf0[ifa] = mf[ifa];
        }
        mf_init = 2;
      }
      if (!checkDataFlag("mf00")) {
        for (int ifa=0; ifa<nfa; ++ifa) {
          mf00[ifa] = mf0[ifa];
        }
        mf_init = 2;
      }
    }
  }

  void initialHookBcs() {
    if (checkParam("RESET_BC") || checkParam("RESET_ALL"))
      FOR_BCZONE (*it)->resetBc();

    FOR_BCZONE (*it)->initialHook();

    // update the scalar ghosts...
    updateScalars();

  }

  void syncPostState() {

    if (mpi_rank==0)
      cout << "syncPostState()" << endl;

    updateCvData(u,REPLACE_ROTATE_DATA);
    updateCvData(p);
    calcPressureGrad(dpdx,p); // no need to update ghost
    updateCvData(rho);
    const double mu_const = getDoubleParam("MU");
    FOR_ICV_G mu_lam[icv] = mu_const;
    initASgs();
    updateSgs();
    FOR_BCZONE (*it)->syncPostState();

    queryBcs();
  }

  int advanceSolution() {

    if ( (mpi_rank==0) && (step%check_interval==0) )
      cout << " > advanceSolution()" << endl;

    // need to propogate dt down the coarse solvers...
    for (int icg = 0; icg < max(ncg_u,ncg_p); ++icg)
      css[icg]->dt = dt;

    CtiRegister::clearCurrentData();

    updateTimeLevels();

    time += dt;
    timer.split("reset time levels");

    //update boundaries with time-dependent velocities, extrapolate data where needed...
    //FOR_BCZONE (*it)->temporalHookBc();  //consider changing back to setBc
    FOR_BCZONE (*it)->setBc();

    // outer iteration loop...
    for (int iter = 0; iter < n_outer_iters; ++iter) {

      // ensure that the sgs contributions are consistent with the new state..
      updateSgs();
      timer.split("updateSgs");

      // compute body force source terms

      computeBodyForces();

      correctContinuity();
      timer.split("correctContinuity");

      // predict the intermediate velocity state u^*

      if (ncg_u > 0) {
        if (A_diag_source)
          solveUMgAdiag();
        else
          solveUMg();
      }
      else {
        const string solver = getStringParam("MOMENTUM_SOLVER","JACOBI");
        if (solver == "PATR") {
          solveUPatr();
        }
        else if (solver == "TIM") {
          solveUTim();
        }
        else {
          assert(solver == "JACOBI");
          solveU();
        }
      }
      timer.split("solveU");

      //Allow custom inflow turbulence BCs to modify
      //velocity here
      FOR_BCZONE (*it)->modifyU();

      // solve scalars...

      const string solver = getStringParam("SCALAR_SOLVER","JACOBI");
      if (solver == "PATR") {
        solveScalarsPatr();
      }
      else {
        assert(solver == "JACOBI");
        solveScalars();
      }
      timer.split("solveScalars");

      // solve the resultant helmholtz equation for the pressure and ensure
      // the velocity field satisfies the cont eqn..

      solvePAndCorrectU();
      timer.split("solvePAndCorrectU");

    }

    report();
    timer.split("report");

    return 0;

  }

  void updateBaseState() {

    const double mu_const    = getDoubleParam("MU");
    const double rho_const   = getDoubleParam("RHO");

    if (checkParam("NO_WALL_INTERP")) {
      FOR_ICV_G {
        mu_lam[icv]  = mu_const;
        rho_uncorrected[icv]     = rho_const;
      }
    }
    else {
      // vida vor type approximation (should only differ at boundary)...
      FOR_ICV_G {
        mu_lam[icv]  = mu_const*planar_vol_ratio_cv[icv];
        rho_uncorrected[icv] = rho_const*planar_vol_ratio_cv[icv];
      }
    }

    // HACK! to match vida vor bug...
    //for (int icv = ncv; icv < ncv_g;++icv) {
    //  rho[icv] = 0.0;
    //  mu_lam[icv] = 0.0;
   // }

  }

  void parseSgsModel(string& sgs_model) {

    //using namespace DynamicSmagorinsky;

    const string sgs_type = getStringParam("SGS_MODEL", "NONE");

    if (sgs_type == "NONE") {
      sgs_model = "NONE";
    } else if ( sgs_type == "VREMAN") {
      sgs_model = "VREMAN";
    } else if ( sgs_type == "SIGMA") {
      sgs_model = "SIGMA";
    } else if ( sgs_type == "DSM_LOCAL") {
      sgs_model = "DSM_LOCAL";
      //registerDsmLocalVariables(this);
    } else {
      CERR(" > unrecognized sgs model: " << sgs_type);
    }
  }

  void updateTimeLevels();
  void updateSgs();
  void completeSgsStuff();
  void computeSgsNone();
  void computeSgsDsm();
  void computeSgsVreman();
  void computeSgsSigma();
  void calcDSMStuff(double* s_mag, double* lijmij, double* mijmij,
		    const double* rho, const double (*ui)[3],
		    const double (*duidxj)[3][3]);
  void initASgs();
  void calcASgs();

  void calcContinuity(double * cont,const bool b_include_a_sgs = false);
  void correctContinuity();
  void correctContinuityHelmholtz();
  void correctContinuityAlgebraic();

  void solveHelmholtz(double * p,double * rhs,const bool pressure = true);
  void solveHelmholtzPatr(double * phi,double * rhs,const bool b_pressure = true);
  void buildHelmholtzMg(double * &A,double* rhs,const bool b_pressure = true);
  void solveHelmholtzMg(double * phi,double * rhs,const bool b_pressure = true);

  void solveU();
  void solveUPatr();
  void solveUTim();
  void buildUMg(double * &A,double * &At,double (*rhs)[3]);
  void solveUMg();
  void solveUMgAdiag();

  void solveScalars();
  void solveScalarsPatr();

  void solvePAndCorrectU();

  void calcCvGrad(double (*dudx)[3][3],const double (*u)[3]);
  void calcPressureGrad(double (*dpdx)[3],const double * p);

  virtual void momentumSourceHook(double * A,double * At,double (*rhs)[3]) {}
  virtual void momentumSourceHook(double * A,double (*rhs)[3]) {}
  virtual void massSourceHook(double * rhs) {}

  virtual void inflowTurbulenceStatsHook(double (*um)[3],double (*Rd)[3], double (*Rod)[3], BfZone* zone_ptr) {
    CERR(" > user must supply a hook for the statistics of the inflow turbulence");
  }
  virtual void inflowTurbulenceLengthscaleHook(double * length_scale, BfZone* zone_ptr) {
    CERR(" > user must supply a hook for the length scales of the inflow turbulence");
  }

  void queryBcs();

  void report() {

    if ( step%check_interval == 0) {
      dumpRange(rho,ncv_g,"rho");
      dumpRange(u,ncv_g,"u");
      dumpRange(p,ncv_g,"p");
      dumpRange(mu_lam,ncv_g,"mu_lam");
      dumpRange(mu_sgs,ncv_g,"mu_sgs");
      dumpRange(a_sgs,ncv_g,"a_sgs");
      dumpRange(mf,nfa,"mf");
      if (!b_helmholtzCorrection){
        dumpRange(mf0,nfa,"mf0");
        dumpRange(mf00,nfa,"mf00");
      }
      dumpCfl();
      dumpRangeScalars();

      // global conservation...
      double my_buf[4] = { 0.0, 0.0, 0.0, 0.0 };
      FOR_ICV {
	my_buf[0] += rho[icv]*vol_cv[icv];
	my_buf[1] += rho[icv]*u[icv][0]*vol_cv[icv];
	my_buf[2] += rho[icv]*u[icv][1]*vol_cv[icv];
	my_buf[3] += rho[icv]*u[icv][2]*vol_cv[icv];
      }
      double buf[4];
      MPI_Reduce(my_buf,buf,4,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      if (mpi_rank == 0) {
	cout << " > time: " << time << " total mass, momentum: " <<
	  buf[0] << " " << buf[1] << " " << buf[2] << " " << buf[3] << endl;
      }

    }

    queryBcs();

    if ((serial_state_interval > 0)&&(step%serial_state_interval == 0))
      writeSerialState(serial_state_prefix,step);

  }

  double calcCfl(double* cfl_, const double dt_) const;

  void dumpCfl() const {
    double * cfl = new double[ncv];
    calcCfl(cfl,dt);
    MiscUtils::dumpRange(cfl,ncv,"CFL");
    delete[] cfl;
  }

  virtual CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,const string& name,
                                                    list<CtiRegister::CtiData>& args, const bool b_eval_func);

  void computeForces() {
    assert(forceZoneBuf&&!forceZoneMap.empty());
    for (map<const string,int>::const_iterator iter = forceZoneMap.begin(); iter!=forceZoneMap.end(); ++iter){
      double f_buf[3][3], m_buf[3][3];
      bc_map[iter->first]->force(f_buf,m_buf);
      FOR_I3 FOR_J3 forceZoneBuf[iter->second][3*i+j] = f_buf[i][j];
      FOR_I3 FOR_J3 momentZoneBuf[iter->second][3*i+j] = m_buf[i][j];
    }
  }

  void computeBodyForces();
  virtual void initBodyForceHook(int &my_implicit_init_flag);
  virtual void computeBodyForceHook();

  void addFrameRotationSourceTerms(double (*rhs)[3]);

  void computeSlipLength();
  void filterCvR1_mod(double *phif, const double *phi);
  void filterCvR2_mod(double (*phif)[3], const double (*phi)[3]);
  void filterCvR3_mod(double (*phif)[3][3], const double (*phi)[3][3]);

  void writeSerialState(const string prefix,const int step) {
    static bool b_first = true;

    if (mpi_rank == 0)
      cout << " > writeSerialState()" << endl;

    if (b_first) {
      // index and position...
      char filename[128];
      sprintf(filename,"%s.%s.%06d.state",prefix.c_str(),"id",mpi_rank);
      MiscUtils::mkdir_for_file(filename);
      string tmp_filename = MiscUtils::makeTmpPrefix(filename);
      remove(tmp_filename.c_str());
      FILE* fp = fopen(tmp_filename.c_str(),"wb");
      fwrite(icv_global,sizeof(int8),ncv,fp);
      fclose(fp);
      remove(filename);
      rename(tmp_filename.c_str(),filename);
    }

    // pressure and velocity...
    {
      char filename[128];
      sprintf(filename,"%s.%s.%06d.%08d.state",prefix.c_str(),"pu",mpi_rank,step);
      MiscUtils::mkdir_for_file(filename);
      string tmp_filename = MiscUtils::makeTmpPrefix(filename);
      remove(tmp_filename.c_str());
      FILE* fp = fopen(tmp_filename.c_str(),"wb");
      fwrite(p,sizeof(double),ncv,fp);
      fwrite((double*)u,sizeof(double),3*ncv,fp);
      fclose(fp);
      remove(filename);
      rename(tmp_filename.c_str(),filename);
    }

    b_first = false;
  }

};

#undef FOR_BCZONE
#endif
