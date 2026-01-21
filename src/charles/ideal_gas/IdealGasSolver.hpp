#ifndef IDEALGASSOLVER_HPP
#define IDEALGASSOLVER_HPP

#include "FlowSolver.hpp"
#include "IdealGasSolverBcs.hpp"
#include "IdealGasSolverFlux.hpp"
#include "CommContainers.hpp"
#include "Sgs.hpp"
#include "RkWeights.hpp"

using namespace RkWeights;

// convenience macro to loop over the boundary conditions...
#define FOR_BCZONE for(vector<IdealGasBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)

class IdealGasSolver : public FlowSolver {
public:

  double *rho;          // density
  double (*u)[3];       // velocity
  double *rhoE;         // total energy
  double *p;            // pressure
  double *T;            // temperature..
  double *sos;          // speed of sound.

  double *mu_lam;       // lamainar visc
  double *mu_sgs;       // sgs eddy viscosity
  double *loc_lam;      // lambda on cp (molecular)
  double *loc_sgs;      // .. and sgs contribution
  double *a_sgs;        // acoustic sgs (density perturbation) 

  double (*dudx)[3][3]; // velocity gradient.

  double *vv2;
  double *fgr_sponge;   // user-def. region to full apply fgr 

  double R_gas;
  double gamma;
  double rho_ref;
  double p_ref;
  double T_ref;
  double mu_ref;
  double mu_power_law;
  double Pr_lam;

  string sgs_model;

  // the flux function is going to use some packed
  // states for efficiency which are defined here...

  IdealGasState* cv_light;
  AuxGasProp* cv_compact;

  IgComm* comm_container;
  SgsComm* sgs_comm_container;

  // for an rk3 update, we'll allocate the 3 rhs..

  IdealGasRhs* rhs0;
  IdealGasRhs* rhs1;
  IdealGasRhs* rhs2;

  // boundary conditions

  vector<IdealGasBc*> bcs;
  map<string,IdealGasBc*> bc_map;
  map<string,pair<int,IdealGasBc*> > bc_queries;

  // mass flux stored on extended faces -- used for passive scalars

  double * mf;
  double * rhs0_sc;
  double * rhs1_sc;
  double * rhs2_sc;
  double * rho_old;

  // for use with a rotating reference frame ..

  double * frame_rotation;


  IdealGasSolver() {

    cv_light   = NULL;
    cv_compact = NULL;

    comm_container     = NULL;
    sgs_comm_container = NULL;

    rhs0 = NULL;
    rhs1 = NULL;
    rhs2 = NULL;

    // register data with i/o

    rho  = NULL; registerCvData(rho , "rho" , READWRITE_DATA);
    u    = NULL; registerCvData(u   , "u"   , READWRITE_DATA);
    rhoE = NULL; registerCvData(rhoE, "rhoE", READWRITE_DATA);
    p    = NULL; registerCvData(p   , "p"   , READWRITE_DATA);
    T    = NULL; registerCvData(T   , "T"   , READWRITE_DATA);

    // data registration without i/o

    sos     = NULL; registerCvData(sos    , "sos"    , CAN_WRITE_DATA);
    mu_lam  = NULL; registerCvData(mu_lam , "mu_lam" , CAN_WRITE_DATA);
    mu_sgs  = NULL; registerCvData(mu_sgs , "mu_sgs" , CAN_WRITE_DATA);
    loc_lam = NULL; registerCvData(loc_lam, "loc_lam", CAN_WRITE_DATA);
    loc_sgs = NULL; registerCvData(loc_sgs, "loc_sgs", CAN_WRITE_DATA);
    a_sgs   = NULL; registerCvData(a_sgs  , "a_sgs"  , CAN_WRITE_DATA);
    vv2     = NULL; registerCvData(vv2, "vv2", CAN_WRITE_DATA);
    fgr_sponge = NULL; 

    mf        = NULL;
    rhs0_sc   = NULL;
    rhs1_sc   = NULL;
    rhs2_sc   = NULL;
    rho_old   = NULL;

    dudx      = NULL;

    if ( checkParam("RES_TIME"))
      _registerScalar("RES_TIME");

    parseSgsModel(sgs_model);

    frame_rotation = NULL;

    if (Param * param = getParam("ROTATING_REFERENCE_FRAME")) {

      // allocated so logic can use NULL checks below

      frame_rotation = new double[3];
      frame_rotation[0] = param->getDouble(0);
      frame_rotation[1] = param->getDouble(1);
      frame_rotation[2] = param->getDouble(2);

      if (mpi_rank == 0)
        cout << " > FRAME_ROTATION active with the following rev/time: " << frame_rotation[0] << " " << frame_rotation[1] << " " << frame_rotation[2] << endl;

      // convert to rads/time unit...
      FOR_I3 frame_rotation[i] *= 2.0*M_PI;

    }

    // use this opportunity to parse the variables needed for the eos...

    gamma        = getDoubleParam("GAMMA", 1.4);
    p_ref        = getDoubleParam("P_REF");
    rho_ref      = getDoubleParam("RHO_REF");
    T_ref        = getDoubleParam("T_REF");
    mu_ref       = getDoubleParam("MU_REF");
    mu_power_law = getDoubleParam("MU_POWER_LAW", 3.0/4.0);
    Pr_lam       = getDoubleParam("PR_LAM", 0.7);

    R_gas        = p_ref/rho_ref/T_ref;
  }

  //=====================================================================
  // boundary conditions
  //====================================================================

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
          bcs.push_back(new SlipWallBc(&bfZoneVec[izone],this));
        } else if ( bc_type == "WALL_ADIABATIC") {
          bcs.push_back(new WallAdiabatic(&bfZoneVec[izone],this));
        } else if ( bc_type == "WALL_ISOTHERMAL") {
          bcs.push_back(new WallIsothermal(&bfZoneVec[izone],this));
        } else if ( bc_type == "WALL_CHT") {
          bcs.push_back(new WallCht(&bfZoneVec[izone],this));
        }  else if ( bc_type == "WALL_CONDUCTIVE") {
          bcs.push_back(new WallConductive(&bfZoneVec[izone],this));
        } else if ( bc_type == "CBC_RUP") {
          bcs.push_back(new CbcRup(&bfZoneVec[izone],this));
        } else if ( bc_type == "CBC_RUNP") {
          bcs.push_back(new CbcRunp(&bfZoneVec[izone],this));
        } else if ( bc_type == "CBC_MPT") {
          bcs.push_back(new CbcMpt(&bfZoneVec[izone],this));
        } else if ( bc_type == "CBC_PROFILE") {
          bcs.push_back(new CbcProfile(&bfZoneVec[izone],this));
        } else if ( bc_type == "CBC_UPT") {
          bcs.push_back(new CbcUpt(&bfZoneVec[izone],this));
        } else if ( bc_type == "CBC_EXTRAPOLATE") {
          bcs.push_back(new CbcExtrapolate(&bfZoneVec[izone],this));
        } else if ( bc_type == "NSCBC_MT") {
          bcs.push_back(new NscbcMt(&bfZoneVec[izone],this));
        } else if ( (bc_type == "NSCBC_OUTLET_P") || (bc_type == "NSCBC_OUTLET_PRESSURE")||
                    (bc_type == "NSCBC_OUTLET_MDOT")) {
          bcs.push_back(new NscbcOutletP(&bfZoneVec[izone],this));
        } else if ( bc_type == "CBC_TOTAL_PT") {
          bcs.push_back(new CbcTotalPt(&bfZoneVec[izone],this));
        } else if ( bc_type == "WM_ALG_ADIABATIC") {
          bcs.push_back(new WmAdiabatic(&bfZoneVec[izone],this,WM_ALGEBRAIC_ADIABATIC));
        } else if ( bc_type == "WM_EQ_ADIABATIC") {
          bcs.push_back(new WmAdiabatic(&bfZoneVec[izone],this,WM_EQUILIBRIUM_ADIABATIC));
        } else if ( bc_type == "WM_ALG_ISOTHERMAL") {
          bcs.push_back(new WmIsothermal(&bfZoneVec[izone],this,WM_ALGEBRAIC_ISOTHERMAL));
        } else if ( bc_type == "WM_EQ_ISOTHERMAL") {
          bcs.push_back(new WmIsothermal(&bfZoneVec[izone],this,WM_EQUILIBRIUM_ISOTHERMAL));
        } else if ( bc_type == "WM_ALG_CHT") { 
          bcs.push_back(new WmAlgCht(&bfZoneVec[izone],this));
        } else if ((bc_type == "WM_EXCHANGE_ADIABATIC") || (bc_type == "WM_EXCHANGE_ISOTHERMAL")) {
          bcs.push_back(new WmExchange(&bfZoneVec[izone],this));
        } else if (bc_type == "INFLOW_TURB"){
          bcs.push_back(new InflowTurbulenceIG(&bfZoneVec[izone],this));
        } else if ( bc_type == "SPONGE") {
          bcs.push_back(new Sponge(&bfZoneVec[izone],this));
        } else if ( bc_type == "CBC_PROFILE_SIMPLE") { 
          bcs.push_back(new CbcProfileSimple(&bfZoneVec[izone],this));
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

  virtual IdealGasBc* initHookBc(BfZone* p) {
    CERR(" > user must supply a hook boundary condition for zone: " << p->getName());
  }

  IdealGasBc* getBc(const string& name) const {
    map<string,IdealGasBc*>::const_iterator it = bc_map.find(name);
    if ( it == bc_map.end()) {
      return NULL;
    } else {
      return it->second;
    }
  }


  //=====================================================
  // initialization, initial hooks, etc...
  // implementation found in IdealGasSolver.cpp
  //====================================================

  void initData();

  void initFromParams();
  void initFromPotentialFlow();

  void initComplete() {

    computeFgrSponge();
    completeSgsStuff();

    // need some logic to compute rhoE sometimes (interpolating from fluent)

    if (!checkDataFlag("rhoE")) {
      if (checkDataFlag("u")&&checkDataFlag("rho")&&checkDataFlag("p")) {
        FOR_ICV rhoE[icv] = p[icv]/(gamma-1.0) + 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);
      }
      else {
        CWARN("No IC set for some conserved variables; did you have an INIT_* or initialHook defined?");
      }
    }

    updateConservativeAndPrimitiveData();
    calcSgsAndMaterialProperties();

    // the initial temperature in the fem/cht solver needs to be sent
    // to the flow boundary conditions for the first fluid side flux
    // calc. Here we use the blocking version since it is only called
    // at the start...

    if (cht) cht->updateT();

  }

  void initialHookBcs() {

    // not exactly clear -- there is a bunch of non-virtual initializations for
    // volume and bf data that should be handled without requiring the user
    // to write a unique "initialHook": e.g. take nearby cv values to initialize
    // surface data in NSCBC-type bcs when not set from data file, etc...

    FOR_BCZONE (*it)->initialHook();
    
    // this is the last routine called before the run loop, so use it as an opportunity
    // to store the current conditions. In this way, we can use robust timestepping
    // from the very beginning...
    storeData();
    
  }

  void initMin() {
    FlowSolver::initMin();
  }

  void syncPostState() {

    updateConservativeAndPrimitiveData();
    calcSgsAndMaterialProperties();

    // one bit of reporting that isnt handled in processStep is the
    // boundary condition querying -- so this is called explicitly
    queryBcs();
  }

  void storeData() {

    if (mpi_rank == 0) cout << "IdealGasSolver::storeData(): step: " << step << " time: " << time << endl;

    int count = ncv*5; // rho,u(3),rhoE TODO: add scalars
    FOR_BCZONE count += (*it)->storeData(NULL); // pass NULL to return the count

    // set/check the count: it should not change (unless particles get added eventually)...
    if (store_buf_size == -1) store_buf_size = count;
    assert(store_buf_size == count);

    // if there is already something in "recent", move it to the more distant buf
    // used to restore...
    if (store_buf_recent) {
      assert(store_buf);
      memcpy(store_buf,store_buf_recent,sizeof(double)*store_buf_size);
      store_step = store_step_recent;
      store_time = store_time_recent;
    }
    else {
      // must be the first time...
      store_buf_recent = new double[store_buf_size];
    }

    // and store the data in recent...
    store_step_recent = step;
    store_time_recent = time;
    count = 0;
    FOR_ICV store_buf_recent[count++] = rho[icv];
    FOR_ICV FOR_I3 store_buf_recent[count++] = u[icv][i];
    FOR_ICV store_buf_recent[count++] = rhoE[icv];
    FOR_BCZONE count += (*it)->storeData(store_buf_recent+count);
    assert(count == store_buf_size);

    // if this is the first store, also push it to store_buf...
    if (store_buf == NULL) {
      store_buf = new double[store_buf_size];
      memcpy(store_buf,store_buf_recent,sizeof(double)*store_buf_size);
      store_step = store_step_recent;
      store_time = store_time_recent;
    }

  }

  void restoreData() {

    if (mpi_rank == 0) {
      cout << "WARNING: NaNs detected in step: " << step << " time: " << time << 
        "\n > restoring solution from step: " << store_step << " time: " << store_time <<
	"\n > adjusting TIMESTEP control from: " << dt_data[0] << " to: " << dt_data[0]*0.8 << endl;
    }
    
    step = store_step;
    time = store_time;
    int count = 0;
    assert(store_buf);
    FOR_ICV rho[icv] = store_buf[count++];
    FOR_ICV FOR_I3 u[icv][i] = store_buf[count++];
    FOR_ICV rhoE[icv] = store_buf[count++];
    FOR_BCZONE count += (*it)->restoreData(store_buf+count);
    assert(count == store_buf_size);

    // after restoring data, everything needs to be re-computed...

    updateConservativeAndPrimitiveData();
    calcSgsAndMaterialProperties();

    // and reduce the timestep by some factor...
    // for either CFL or DT mode, the timestep is in FlowSolver::dt_data[0]...

    dt_data[0] *= 0.8;

  }

  ~IdealGasSolver();

  //=======================================================================
  // solver stuff -- advanceSolution, rhs calcs ..
  // implementation for these are available in IdealGasSolver.cpp
  //======================================================================

  int advanceSolution() {

    // if using cht (conjugate heat transfer), we should come into this routine
    // with the cht temperature already updated into the flow side CHT bc's:
    // see initComplete() above.

    // zero all rhs's here to allow the rk3Step routine
    // to get vectorized by the compiler...

    for (int icv = 0; icv < ncv; ++icv) {
      rhs0[icv].zero();
      rhs1[icv].zero();
      rhs2[icv].zero();
    }

    if ( rhs0_sc ) {
      for (int ii = 0; ii < nsc*ncv; ++ii) {
        rhs0_sc[ii] = 0.0;
        rhs1_sc[ii] = 0.0;
        rhs2_sc[ii] = 0.0;
      }
    }

    LpSolidRhs  * lpSolidRhs_arr[3];
    LpSolidRhs  * lpSolidRhs_d = NULL;
    LpLiquidRhs * lpLiquidRhs_arr[3];
    LpLiquidRhs * lpLiquidRhs_d = NULL;
    IdealGasRhs * rhs_arr[3] = {rhs0, rhs1, rhs2};

    preprocessLp(lpSolidRhs_arr, lpLiquidRhs_arr, lpSolidRhs_d, lpLiquidRhs_d);

    LpSolidRhs * lpSolidRhs0  = lpSolidRhs_arr[0];
    LpSolidRhs * lpSolidRhs1  = lpSolidRhs_arr[1];
    LpSolidRhs * lpSolidRhs2  = lpSolidRhs_arr[2];

    LpLiquidRhs * lpLiquidRhs0  = lpLiquidRhs_arr[0];
    LpLiquidRhs * lpLiquidRhs1  = lpLiquidRhs_arr[1];
    LpLiquidRhs * lpLiquidRhs2  = lpLiquidRhs_arr[2];

    // SB/CI why not merge bc->calcRhs and (*it)->addBoundaryFlux(rhs) that
    // appears at the top of solver::calcRhs ???

    // rk3 update...
    FOR_BCZONE (*it)->calcRhs(time,1);
    calcRhs(rhs0,rhs0_sc,time,1);
    if (cht) cht->updateQStart();
    if (lpTracer) lpTracer->rk3Step(dt, 1);
    if ( (lpSolid)&&(lpCoupling==ONE_WAY_COUPLED) ) calcRhsLpSolid_1w(lpSolidRhs0, lpSolidRhs_d);
    if ( (lpSolid)&&(lpCoupling==TWO_WAY_COUPLED) ) calcRhsLpSolid_2w(lpSolidRhs0);
    if ( (lpLiquid)&&(lpCoupling==ONE_WAY_COUPLED) ) calcRhsLpLiquid_1w(lpLiquidRhs0, lpLiquidRhs_d);
    if ( (lpLiquid)&&(lpCoupling==TWO_WAY_COUPLED) ) calcRhsLpLiquid_2w(lpLiquidRhs0);
    timer.split("calc rhs 1");

    FOR_BCZONE (*it)->rk3Step(dt,1);
    if ( (lpSolid)||(lpLiquid) ) {
      if (lpCoupling==TWO_WAY_COUPLED) {
        rk3Step_lp_gas(erk_wgt3[0], rhs_arr, lpSolidRhs_arr, lpLiquidRhs_arr, 1);
      }
      else if (lpCoupling==ONE_WAY_COUPLED) {
        if (lpSolid) rk3Step_lpSolid(erk_wgt3[0], lpSolidRhs_arr, lpSolidRhs_d, 1);
        if (lpLiquid) rk3Step_lpLiquid(erk_wgt3[0], lpLiquidRhs_arr, lpLiquidRhs_d, 1);
        rk3Step(erk_wgt3[0],dt,1);
      }
      updateLpDp();
    }
    else {
      rk3Step(erk_wgt3[0],dt,1);
    }
    if (cht) {
      cht->updateQFinish();
      cht->rk3Step(erk_wgt3[0],dt,1); // assembles and solves the fem system
    }
    timer.split("rk step 1");

    time += dt;
    if (cht) cht->updateTStart();
    updateConservativeAndPrimitiveData();
    CtiRegister::clearCurrentData();
    if (cht) cht->updateTFinish();

    timer.split("update 1");

    // stage 2...
    FOR_BCZONE (*it)->calcRhs(time,2);
    calcRhs(rhs1,rhs1_sc,time,2);
    if (cht) cht->updateQStart();
    if (lpTracer) lpTracer->rk3Step(dt, 2);
    if ( (lpSolid)&&(lpCoupling==ONE_WAY_COUPLED) ) calcRhsLpSolid_1w(lpSolidRhs1, lpSolidRhs_d);
    if ( (lpSolid)&&(lpCoupling==TWO_WAY_COUPLED) ) calcRhsLpSolid_2w(lpSolidRhs1);
    if ( (lpLiquid)&&(lpCoupling==ONE_WAY_COUPLED) ) calcRhsLpLiquid_1w(lpLiquidRhs1, lpLiquidRhs_d);
    if ( (lpLiquid)&&(lpCoupling==TWO_WAY_COUPLED) ) calcRhsLpLiquid_2w(lpLiquidRhs1);
    timer.split("calc rhs 2");

    FOR_BCZONE (*it)->rk3Step(dt,2);
    if ( (lpSolid)||(lpLiquid) ) {
      if (lpCoupling==TWO_WAY_COUPLED) {
        rk3Step_lp_gas(erk_wgt3[1], rhs_arr, lpSolidRhs_arr, lpLiquidRhs_arr, 2);
      }
      else if (lpCoupling==ONE_WAY_COUPLED) {
        if (lpSolid) rk3Step_lpSolid(erk_wgt3[1], lpSolidRhs_arr, lpSolidRhs_d, 2);
        if (lpLiquid) rk3Step_lpLiquid(erk_wgt3[1], lpLiquidRhs_arr, lpLiquidRhs_d, 2);
        rk3Step(erk_wgt3[1],dt,2);
      }
      updateLpDp();
    }
    else {
      rk3Step(erk_wgt3[1],dt,2);
    }
    if (cht) {
      cht->updateQFinish();
      cht->rk3Step(erk_wgt3[1],dt,2);
    }
    timer.split("rk step 2");

    time -= 0.5*dt;
    if (cht) cht->updateTStart();
    updateConservativeAndPrimitiveData();
    CtiRegister::clearCurrentData();
    if (cht) cht->updateTFinish();

    timer.split("update 2");

    // stage 3...
    FOR_BCZONE (*it)->calcRhs(time,3);
    calcRhs(rhs2,rhs2_sc,time,3);
    if (cht) cht->updateQStart();
    if (lpTracer) lpTracer->rk3Step(dt, 3);
    if ( (lpSolid)&&(lpCoupling==ONE_WAY_COUPLED) ) calcRhsLpSolid_1w(lpSolidRhs2, lpSolidRhs_d);
    if ( (lpSolid)&&(lpCoupling==TWO_WAY_COUPLED) ) calcRhsLpSolid_2w(lpSolidRhs2);
    if ( (lpLiquid)&&(lpCoupling==ONE_WAY_COUPLED) ) calcRhsLpLiquid_1w(lpLiquidRhs2, lpLiquidRhs_d);
    if ( (lpLiquid)&&(lpCoupling==TWO_WAY_COUPLED) ) calcRhsLpLiquid_2w(lpLiquidRhs2);
    timer.split("calc rhs 3");


    FOR_BCZONE (*it)->rk3Step(dt,3);
    if ( (lpSolid)||(lpLiquid) ) {
      if (lpCoupling==TWO_WAY_COUPLED) {
        rk3Step_lp_gas(erk_wgt3[2], rhs_arr, lpSolidRhs_arr, lpLiquidRhs_arr, 3);
      }
      else if (lpCoupling==ONE_WAY_COUPLED) {
        if (lpSolid) rk3Step_lpSolid(erk_wgt3[2], lpSolidRhs_arr, lpSolidRhs_d, 3);
        if (lpLiquid) rk3Step_lpLiquid(erk_wgt3[2], lpLiquidRhs_arr, lpLiquidRhs_d, 3);
        rk3Step(erk_wgt3[2],dt,3);
      }
      updateLpDp();
    }
    else {
      rk3Step(erk_wgt3[2],dt,3);
    }
    if (cht) {
      cht->updateQFinish();
      cht->rk3Step(erk_wgt3[2],dt,3);
    }
    timer.split("rk step 3");

    time += 0.5*dt;
    if (cht) cht->updateTStart();
    updateConservativeAndPrimitiveData();
    CtiRegister::clearCurrentData();
    if (cht) cht->updateTFinish();

    timer.split("update 3");

    // compute sgs model and ensurue that material properties are
    // consistent with the new state that's been advanced..

    calcSgsAndMaterialProperties();
    timer.split("calc sgs and properties");

    if (checkForNans()) {
      restoreData();
      return -1;
    }

    report();
    postprocessLp(lpSolidRhs_arr, lpLiquidRhs_arr, lpSolidRhs_d, lpLiquidRhs_d);
    timer.split("reporting");

    return 0;

  }

  bool checkForNans() {
    
    // before running the report, look for nans in rho, p, T...
    int my_got_nans = 0;
    FOR_ICV {
      if ((rho[icv] != rho[icv])||(p[icv] != p[icv])||(T[icv] != T[icv])) {
        my_got_nans = 1;
        break;
      }
    }
    int got_nans;
    MPI_Allreduce(&my_got_nans,&got_nans,1,MPI_INT,MPI_MAX,mpi_comm);
    
    return (got_nans == 1);
    
  }

  void advanceSolutionOrig() {

    // zero all rhs's here to allow the rk3Step routine
    // to get vectorized by the compiler...

    for (int icv = 0; icv < ncv; ++icv) {
      rhs0[icv].zero();
      rhs1[icv].zero();
      rhs2[icv].zero();
    }

    if ( rhs0_sc ) {
      for (int ii = 0; ii < nsc_transport*ncv; ++ii) {
        rhs0_sc[ii] = 0.0;
        rhs1_sc[ii] = 0.0;
        rhs2_sc[ii] = 0.0;
      }
    }

    // rk3 update...

    FOR_BCZONE (*it)->calcRhs(time,1);
    calcRhs(rhs0,rhs0_sc,time,1);
    timer.split("calc rhs 1");

    FOR_BCZONE (*it)->rk3Step(dt,1);
    rk3Step(erk_wgt3[0],dt,1);
    timer.split("rk step 1");

    time += dt;
    updateConservativeAndPrimitiveData();
    CtiRegister::clearCurrentData();

    timer.split("update 1");

    // stage 2...

    FOR_BCZONE (*it)->calcRhs(time,2);
    calcRhs(rhs1,rhs1_sc,time,2);
    timer.split("calc rhs 2");

    FOR_BCZONE (*it)->rk3Step(dt,2);
    rk3Step(erk_wgt3[1],dt,2);
    timer.split("rk step 2");

    time -= 0.5*dt;
    updateConservativeAndPrimitiveData();
    CtiRegister::clearCurrentData();

    timer.split("update 2");

    // stage 3...

    FOR_BCZONE (*it)->calcRhs(time,3);
    calcRhs(rhs2,rhs2_sc,time,3);
    timer.split("calc rhs 3");

    FOR_BCZONE (*it)->rk3Step(dt,3);
    rk3Step(erk_wgt3[2],dt,3);
    timer.split("rk step 3");

    time += 0.5*dt;
    updateConservativeAndPrimitiveData();
    CtiRegister::clearCurrentData();

    timer.split("update 3");

    // compute sgs model and ensurue that material properties are
    // consistent with the new state that's been advanced..

    calcSgsAndMaterialProperties();
    timer.split("calc sgs and properties");

    // dump relevant information to stdout..

    report();
    timer.split("reporting");
  }

  void rk3Step(const double rk_wgt[3], const double dt, const int rk_stage);
  void advanceScalars(const double * rho_old, const double rk_wgt[3],
                      const double time, const int rk_stage);


  void calcRhs(IdealGasRhs *__restrict__ rhs, double *__restrict__ rhs_sc,
               const double time, const int rk_stage) {

    calcRhsBase(rhs,rhs_sc,time,rk_stage);
    calcRhsScalars(rhs_sc, mf, time, rk_stage);
    addSourceHook(rhs,time,rk_stage);
    addRotationSourceTerms(rhs,rhs_sc,time,rk_stage);
    addBodyForceSourceTerms(rhs,rhs_sc,time,rk_stage);
  }


  void calcRhsBase(IdealGasRhs *__restrict__ rhs, double *__restrict__ rhs_sc,
                   const double time, const int rk_stage);
  void calcRhsScalars( double* __restrict__ rhs_sc, const double* __restrict__ mf,
                       const double time, const int rk_stage);

  virtual void addSourceHook(IdealGasRhs* rhs, const double time, const int rk_stage) {}

  void addRotationSourceTerms(IdealGasRhs* __restrict__ rhs, double *__restrict__ rhs_sc,
                              const double time, const int rk_stage) ;

  void addBodyForceSourceTerms(IdealGasRhs* __restrict__ rhs, double *__restrict__ rhs_sc,
                              const double time, const int rk_stage) ;


  void setFgrCvsFull();
  void setFgrCvs();
  void calcSgsAndMaterialProperties();

  void updateExtendedState(IdealGasState* __restrict__ cv_light, AuxGasProp* __restrict__ cv_compact,
                           double *__restrict__ p, double * __restrict__ T, const double* __restrict__ rho,
                           const double (*__restrict__ u)[3], const double* __restrict__ rhoE,
                           const int icv_start, const int icv_end);


  void updateConservativeAndPrimitiveData(const int action = COMPLETE_ALL_UPDATES) {

    // start the communication of the rho,u,rhoE...

    updateCv2ClassStart<IgComm,double>(comm_container);
    updateExtendedState(cv_light,cv_compact,p,T,rho,u,rhoE,0,ncv);
    calcCvGradStart(dudx,u);

    // complete the communication of the ghosts...

    updateCv2ClassFinish<IgComm,double>(comm_container);
    updateExtendedState(cv_light,cv_compact,p,T,rho,u,rhoE,ncv,ncv_g2);
    calcCvGradFinish(dudx,u);

    // update the scalar ghosts -- these are NOT packed updates.

    updateScalars();

    // communicate the grad(u) in the ghosts.  if requested to complete
    // all comm, the routine will block until grad(u) in the ghosts are set;
    // otherwise, it is responsibiltiy of the callign process to complete
    // the cv data update if the grad is postponed...

    updateCvDataStart( dudx, REPLACE_ROTATE_DATA);
    updateCvDataFinish(dudx);
    //if ( action == COMPLETE_ALL_UPDATES)
      //updateCvDataFinish(dudx);
  }

  double calcCfl(double* cfl_, const double dt_) const;
  double calcDfl(double* cfl_, const double dt_) const;

  void computeIgEntropy(double* ent) const;


  //=============================================================
  // sgs routines
  // implementation found in IdealGasSolverSgs.cpp
  //=============================================================

  void parseSgsModel(string& sgs_model) {

    using namespace DynamicSmagorinsky;

    const string sgs_type = getStringParam("SGS_MODEL", "NONE");

    if (sgs_type == "NONE") {
      sgs_model = "NONE";
    } else if ( sgs_type == "VREMAN") {
      sgs_model = "VREMAN";
    } else if ( sgs_type == "SIGMA") {
      sgs_model = "SIGMA";
    } else if ( sgs_type == "DSM_LOCAL") {
      sgs_model = "DSM_LOCAL";
      registerDsmLocalVariables(this);
    } else {
      CERR(" > unrecognized sgs model: " << sgs_type);
    }
  }

  void completeSgsStuff();
  void computeSgsNone();
  void computeSgsDsm();
  void computeSgsVreman();
  void computeSgsSigma();
  void limitAsgs();
  void computeFgrSponge();

  //==============================================================
  // solver diagnostics and io
  // implementation found in IdealGasSolverIo.cpp
  //==============================================================

  void report() {

    if ( step%check_interval == 0) {
      
      // store the data every so often...
      if (step%(check_interval*20) == 0)
        storeData();

      // var ranges
      dumpRange(rho,ncv,"rho");
      dumpRange(u,ncv,"u");
      dumpRange(p,ncv,"p");
      dumpRange(T,ncv, "T");
      dumpRange(mu_lam,ncv,"mu_lam");
      dumpRange(mu_sgs,ncv,"mu_sgs");
      dumpRange(a_sgs,ncv,"a_sgs");
      dumpCfl();
      dumpRangeScalars();

      // global conservation...
      double my_buf[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
      FOR_ICV {
	my_buf[0] += rho[icv]*vol_cv[icv];
	my_buf[1] += rho[icv]*u[icv][0]*vol_cv[icv];
	my_buf[2] += rho[icv]*u[icv][1]*vol_cv[icv];
	my_buf[3] += rho[icv]*u[icv][2]*vol_cv[icv];
	my_buf[4] += rhoE[icv]*vol_cv[icv];
      }
      double buf[5];
      MPI_Reduce(my_buf,buf,5,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      if (mpi_rank == 0) {
	cout << " > time: " << time << " total mass, momentum, energy: " <<
	  buf[0] << " " << buf[1] << " " << buf[2] << " " << buf[3] << " " << buf[4] << endl;
      }
    }

    /*
    if (cht) {
      if (step%100 == 0) cht->writeTecplot(step);
    }
    */

    queryBcs();
  }

  void dumpCfl() const {
    double * cfl = new double[ncv];
    calcDfl(cfl, dt);
    MiscUtils::dumpRange(cfl,ncv,"DFL");
    calcCfl(cfl,dt);
    MiscUtils::dumpRange(cfl,ncv,"CFL");
    delete[] cfl;
  }

  void queryBcs();

  // custom variable output via cti var evaluation; recall on completion
  // if the data is not found, then it must be passed down to the base
  // StaticSolver class to see if it can evaluate the data or otherwise error

  virtual CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,const string& name,
      list<CtiRegister::CtiData>& args, const bool b_eval_func);

  void computeForces() {
    assert(forceZoneBuf&&!forceZoneMap.empty());
    for (map<const string,int>::const_iterator iter = forceZoneMap.begin(); iter!=forceZoneMap.end(); ++iter){
      double f_buf[3][3], m_buf[3][3];
      bc_map[iter->first]->force(f_buf,m_buf,p_ref);
      FOR_I3 FOR_J3 forceZoneBuf[iter->second][3*i+j] = f_buf[i][j];
      FOR_I3 FOR_J3 momentZoneBuf[iter->second][3*i+j] = m_buf[i][j];
    }
  }

  //===================
  // Lp stuff...
  //===================

  void preprocessLp(LpSolidRhs * lpSolidRhs_arr[3], LpLiquidRhs * lpLiquidRhs_arr[3], LpSolidRhs * &lpSolidRhs_d, LpLiquidRhs * &lpLiquidRhs_d);
  void postprocessLp(LpSolidRhs * lpSolidRhs_arr[3], LpLiquidRhs * lpLiquidRhs_arr[3], LpSolidRhs * lpSolidRhs_d, LpLiquidRhs * lpLiquidRhs_d);
  void calcRhsLpSolid_1w(LpSolidRhs * rhs, LpSolidRhs * rhs_d);
  void calcRhsLpSolid_2w(LpSolidRhs * rhs);
  void calcRhsLpLiquid_1w(LpLiquidRhs * rhs, LpLiquidRhs * rhs_d);
  void calcRhsLpLiquid_2w(LpLiquidRhs * rhs);
  void rk3Step_lp_gas(const double* rk_wgt, IdealGasRhs ** rhs_arr, LpSolidRhs ** lpSolidRhs_arr, LpLiquidRhs ** lpLiquidRhs_arr, const int rkstep);
  inline void solveLinearSysGE(double* X[3],const double* const A_d, const double* const v, const double* const q, const double D, double* const rhs[3], const int np) const;
  inline void solveLinearSysGE(double* X, const double* const A_d, const double * const v, const double* const q, const double D,double* const rhs, const int np) const;
  void rk3Step_lpSolid(const double rk_wgt[3], LpSolidRhs ** lpSolidRhs_arr, LpSolidRhs * rhs_d, const int rkstep);
  void rk3Step_lpLiquid(const double rk_wgt[3], LpLiquidRhs ** lpLiquidRhs_arr, LpLiquidRhs * rhs_d, const int rkstep);

};

#undef FOR_BCZONE
#endif
