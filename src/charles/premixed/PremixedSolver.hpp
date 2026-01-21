#ifndef PREMIXEDSOLVER_HPP
#define PREMIXEDSOLVER_HPP

#include "FlowSolver.hpp"
#include "PremixedSolverFlux.hpp"
#include "PremixedSolverBcs.hpp"
#include "RkWeights.hpp"
#include "CommContainers.hpp"
#include "Chemtable.hpp"

using namespace RkWeights;

#define FOR_BCZONE for(vector<PremixedBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)

class PremixedSolver : public FlowSolver {
protected:

  bool b_reset_state;
  
public:

  double * rho;          // density
  double (*u)[3];        // velocity
  double * rhoE;         // total energy
  double * Z;            // mixture fraction
  double * C;            // progress variable
  double * p;            // pressure
  double * T;            // temperature
  double * sos;          // speed of sound..
  double * ent;          // entropy..
  double * h;            // enthalpy

  double (*dudx)[3][3];  // velocity gradient

  double *vv2;

  double * mu_lam;   // laminar visc
  double * mu_sgs;   // sgs eddy visc
  double * loc_lam;  // lambda on cp (molecular)
  double * loc_sgs;  // .. sgs contribution
  double * a_sgs;    // acoustic sgs (density perturbation sgs)

  double * R;        // gas constant
  double * e_p;      // internal energy at the chemtable pressure
  double * gamma;    // specific heat ratio
  double * a_gam;    // sp heat ratio linear expansion coeff
  double * T_p;      // nominal temperature at the chemtable pressure

  // perturbation quantities; defined wrt to the state looked up
  // from the flamelets..

  double * h_prime;
  double * s_prime;

  // for the face based (integrated) premixed combustion model..
  double * int_rho_csrc;

  // sgs & combustion models
  string sgs_model;
  string eff_model; //flame efficiency, defined at the faces...
  double * e_fa;

  // flux functions use special data containers to improve
  // memory locality and cache efficiency...
  PremixedState* cv_light;
  AuxGasProp* cv_compact;

  // packed communicator patterns...
  PremixedComm* comm_container;
  SgsComm* sgs_comm_container;

  // rk3 update uses 3 rhs..
  PremixedRhs* rhs0;
  PremixedRhs* rhs1;
  PremixedRhs* rhs2;

  AbstractChemtable2D* chemtable;

  // boundary conditions
  vector<PremixedBc*> bcs;
  map<string,PremixedBc*> bc_map;
  map<string,pair<int,PremixedBc*> > bc_queries;

  // prefactor coefficient to multiply gibbs remainder
  // based stabilization to the flux (perhaps change name)


  // eff_coeff is part of the Charolote efficiency function...

  double eff_coeff;

  // mass flux stored on extended faces -- used to help
  // transport scalars.

  double * mf;
  double * rhs0_sc;
  double * rhs1_sc;
  double * rhs2_sc;
  double * rho_old;

  // nox model -- the nox is handled through the passive
  // scalar infrastructure, but some special logic is added
  // to amortize some of the chemtable lookups that happen.

  bool b_nox;

  PremixedSolver() {

    b_reset_state = false;
    
    rho     = NULL; registerCvData(rho , "rho" , READWRITE_DATA);
    u       = NULL; registerCvData(u   , "u"   , READWRITE_DATA);
    rhoE    = NULL; registerCvData(rhoE, "rhoE", READWRITE_DATA);
    Z       = NULL; registerCvData(Z   , "Z"   , READWRITE_DATA);
    C       = NULL; registerCvData(C   , "C"   , READWRITE_DATA);
    p       = NULL; registerCvData(p   , "p"   , READWRITE_DATA);
    T       = NULL; registerCvData(T   , "T"   , READWRITE_DATA);
    sos     = NULL; registerCvData(sos , "sos" , READWRITE_DATA);
    ent     = NULL; registerCvData(ent , "ent" , READWRITE_DATA);

    mu_lam  = NULL; registerCvData(mu_lam , "mu_lam" , CAN_WRITE_DATA);
    mu_sgs  = NULL; registerCvData(mu_sgs , "mu_sgs" , CAN_WRITE_DATA);
    loc_lam = NULL; registerCvData(loc_lam, "loc_lam", CAN_WRITE_DATA);
    loc_sgs = NULL; registerCvData(loc_sgs, "loc_sgs", CAN_WRITE_DATA);
    a_sgs   = NULL; registerCvData(a_sgs  , "a_sgs"  , CAN_WRITE_DATA);
    h       = NULL; registerCvData(h   , "h"         , CAN_WRITE_DATA);

    R       = NULL; registerCvData(R , "R"           , CAN_WRITE_DATA);
    e_p     = NULL; registerCvData(e_p, "e_p"        , CAN_WRITE_DATA); // internal energy at table pressure
    gamma   = NULL; registerCvData(gamma , "gamma"   , CAN_WRITE_DATA);
    a_gam   = NULL; registerCvData(a_gam , "a_gam"   , CAN_WRITE_DATA);
    T_p     = NULL; registerCvData(T_p , "T_p"       , CAN_WRITE_DATA); // temperature at table pressure

    h_prime = NULL; registerCvData(h_prime, "h_pr"   , CAN_WRITE_DATA);
    s_prime = NULL; registerCvData(s_prime, "s_pr"   , CAN_WRITE_DATA);

    vv2     = NULL; registerCvData(vv2, "vv2", CAN_WRITE_DATA);

    int_rho_csrc = NULL;
    e_fa         = NULL;

    cv_light   = NULL;
    cv_compact = NULL;

    comm_container     = NULL;
    sgs_comm_container = NULL;

    rhs0 = NULL;
    rhs1 = NULL;
    rhs2 = NULL;

    chemtable = NULL;

    mf        = NULL;
    rhs0_sc   = NULL;
    rhs1_sc   = NULL;
    rhs2_sc   = NULL;
    rho_old   = NULL;

    b_nox     = getBoolParam("NOX_MODEL", false);
    if ( b_nox)
      _registerScalar("Y_NOX");

    dudx      = NULL;

    if ( checkParam("RES_TIME"))
      _registerScalar("RES_TIME");

    parseSgsModel(sgs_model);
    parseEfficiencyModel(eff_model);
    eff_coeff    = getDoubleParam("EFF_COEFF", 2.0);

    if (mpi_rank == 0) {
      cout << " > sgs_model: " << sgs_model << endl;
      cout << " > eff_model: " << eff_model << endl;
      cout << " > eff_coeff: " << eff_coeff << endl;
    }
  }

  ~PremixedSolver();

  //==========================================
  // initialization functions
  //=========================================

  void init() {
    initPremixedChemistry();
    FlowSolver::init();
  }

  void initMin() {
    initPremixedChemistry();
    FlowSolver::initMin();
  }

  void initPremixedChemistry() {
    initChemtable(chemtable,getStringParam("CHEMTABLE"));
    vector<string> strVec;
    strVec.push_back("rho");
    strVec.push_back("T");
    strVec.push_back("R");
    strVec.push_back("e");
    strVec.push_back("prog");
    strVec.push_back("src_prog");
    strVec.push_back("gamma");
    strVec.push_back("a_gamma");
    strVec.push_back("mu");
    strVec.push_back("a_mu");
    strVec.push_back("locp");
    strVec.push_back("a_locp");
    strVec.push_back("sL");
    strVec.push_back("lF");
    strVec.push_back("mw");
    strVec.push_back("int_rho_src");
    strVec.push_back("s");

    if ( b_nox) {
      strVec.push_back("src_nox_therm");
      strVec.push_back("int_src_nox_prompt");
    }

    chemtable->loadVariables(strVec);
  }

  void initFromParams();
  void initFromPotentialFlow(const bool coldflow);

  virtual void initData();
  virtual void loadBalanceHook(StripedMesh* sm);


  void registerBoundaryConditions() {

    assert( bcs.size() == 0);

    StaticSolver::registerBoundaryConditions(); // registers bf geometric data

    int nerr = 0;
    vector<pair<string,string> > errors;

    FOR_IZONE(bfZoneVec) {
      const string zone_name = bfZoneVec[izone].getName();
      if ( Param* p = getParam(zone_name)) {
        const string bc_type = p->getString(0);
        if ((bc_type == "SLIP")||(bc_type == "SYMMETRY"))
          bcs.push_back(new SlipWallPBc(&bfZoneVec[izone],this));
        else if ((bc_type == "MS_TABLE")||(bc_type == "MS_AT_TABLE"))
          bcs.push_back(new MsTableP(p,&bfZoneVec[izone],this));
        else if (bc_type == "WALL_ADIABATIC")
          bcs.push_back(new WallAdiabaticP(&bfZoneVec[izone],this));
        else if (bc_type == "WALL_ISOTHERMAL")
          bcs.push_back(new WallIsothermalP(&bfZoneVec[izone],this));
        else if (bc_type == "WALL_CONDUCTIVE")
          bcs.push_back(new WallConductiveP(&bfZoneVec[izone],this));
        else if (bc_type == "WALL_CHT")
          bcs.push_back(new WallChtP(&bfZoneVec[izone],this));
        else if (bc_type == "CBC_RUPS")
          bcs.push_back(new CbcRupsP(&bfZoneVec[izone],this));
        else if (bc_type == "CBC_RUNPS")
          bcs.push_back(new CbcRunpsP(&bfZoneVec[izone],this));
        else if (bc_type == "CBC_MPTS")
          bcs.push_back(new CbcMptsP(&bfZoneVec[izone],this));
        else if ((bc_type == "CBC_MS_TABLE")||(bc_type == "CBC_MS_AT_TABLE"))
          bcs.push_back(new CbcMsTableP(&bfZoneVec[izone],this));
        else if (bc_type == "CBC_PROFILE")
          bcs.push_back(new CbcProfileP(&bfZoneVec[izone], this));
        else if (bc_type == "CBC_UPTS")
          bcs.push_back(new CbcUptsP(&bfZoneVec[izone],this));
        else if (bc_type == "CBC_UPTS_SPONGE")
          bcs.push_back(new CbcUptsSpongeP(&bfZoneVec[izone],this));
        else if (bc_type == "CBC_TOTAL_PTS")
          bcs.push_back(new CbcTotalPtsP(&bfZoneVec[izone], this));
        else if (bc_type == "NSCBC_MTS")
          bcs.push_back(new NscbcMtsP(&bfZoneVec[izone],this));
        else if ((bc_type == "NSCBC_OUTLET_P")||(bc_type == "NSCBC_OUTLET_PRESSURE"))
          bcs.push_back(new NscbcOutletPP(&bfZoneVec[izone],this));
        else if (bc_type == "NSCBC_OUTLET_MDOT")
          bcs.push_back(new NscbcOutletMdotP(&bfZoneVec[izone],this));
        else if (bc_type == "NSCBC_MTS_SOFT")
          bcs.push_back(new NscbcMtsSoftP(&bfZoneVec[izone],this));
        else if (bc_type == "SPONGE")
          bcs.push_back(new SpongeP(&bfZoneVec[izone],this));
        else if (bc_type == "WM_ALG_ADIABATIC")
          bcs.push_back(new WmAlgAdiabaticP(&bfZoneVec[izone],this));
        else if (bc_type == "WM_ALG_ISOTHERMAL")
          bcs.push_back(new WmAlgIsothermalP(&bfZoneVec[izone],this));
        else if (bc_type == "WM_ALG_CHT")
          bcs.push_back(new WmAlgChtP(&bfZoneVec[izone],this));
        else if (bc_type == "WM_ALG_CONDUCTIVE")
          bcs.push_back(new WmAlgConductiveP(&bfZoneVec[izone], this));
        else if (bc_type == "INFLOW_TURB")
          bcs.push_back(new InflowTurbulenceP(&bfZoneVec[izone], this));
        else if (bc_type == "HOOK")
          bcs.push_back(initHookBc(&bfZoneVec[izone]));
        else {
          nerr++;
          errors.push_back(pair<string,string>(zone_name,bc_type));
        }
      } else {
        nerr++;
        errors.push_back(pair<string,string>(zone_name,""));
      }
    }

    // ensure that the error handling is synced..; should be uncessary and checked below

    MPI_Bcast(&nerr,1,MPI_INT,0,mpi_comm);
    assert( nerr == int(errors.size()));
    reportBcErrors(errors);

  }

  void initBoundaryConditions() {

    // allow the boundary conditions to allocate their necessary data (not allowed to
    // set anythign here, because it may be coming from the restart data yet to be read)..

    FOR_BCZONE (*it)->initData();

    // create the name based map of the bc objects for easy query functions late..

    FOR_BCZONE bc_map[(*it)->getName()] = *it;
  }


  virtual PremixedBc* initHookBc(BfZone* p) {
    CERR(" > user must supply a hook boundary condition for zone: " << p->getName());
  }

  void resetState();
  
  void initComplete() {

    if (b_reset_state) resetState();
    
    updateConservativeAndPrimitiveData();
    calcSgsAndMaterialProperties();

    // note that bc initial hooks are have been parsed yet..

    // the initial temperature in the fem/cht solver needs to be sent
    // to the flow boundary conditions for the first fluid side flux
    // calc. Here we use the blocking version since it is only called
    // at the start...

    if (cht) cht->updateT();
  }

  void initialHookBcs() {
    FOR_BCZONE (*it)->initialHook();

    // this is the last routine called before the run loop, so use it as an opportunity
    // to store the current conditions. In this way, we can use robust timestepping
    // from the very beginning...
    storeData();

  }

  PremixedBc* getBc(const string& name) const {
    map<string,PremixedBc*>::const_iterator it = bc_map.find(name);
    if ( it == bc_map.end()) {
      return NULL;
    } else {
      return it->second;
    }
  }

  //=========================================
  // sgs and combustion models
  //=========================================
  void parseSgsModel(string& sgs_model) {
    sgs_model = "NONE"; // default
    if (Param * param = getParam("SGS_MODEL")) {
      const string token = param->getString(0);
      if (token == "NONE") {
	sgs_model = "NONE";
      }
      else if (token == "VREMAN") {
	sgs_model = "VREMAN";
      }
      else {
	CERR(" > unrecognized sgs model: " << token);
      }
    }
  }

  void computeSgsNone();
  void computeSgsVreman();
  //void limitAsgs(); // limiting always on and part of the a_sgs calc now

  void parseEfficiencyModel(string& eff_model) {
    eff_model = "CHARLETTE"; // default
    if ( Param* param = getParam("EOS")) {
      int i = 0;
      while( i < param->size()) {
        string token = param->getString(i++);
        if ( (token == "EFFICIENCY") && (i < param->size())) {
          token = param->getString(i);
          if ( token == "CHARLETTE")
            eff_model = "CHARLETTE";
          else if ( token == "LAMINAR")
            eff_model = "LAMINAR";
          else if ( token == "OLD")
            eff_model = "OLD"; // for debugging purposes..
          else
            CERR(" > unknown efficiency type: " << token);
        }
      }
    }
  }

  void computeEfficiencyLaminar();
  void computeEfficiencyCharletteOld();
  void computeEfficiencyCharlette();

  void syncPostState() {

    updateConservativeAndPrimitiveData();
    calcSgsAndMaterialProperties();

    // boundary conditions are not queried as part of processStep --
    // so it needs to be explicitly called here

    queryBcs();

  }

  void storeData() {

    if (mpi_rank == 0) cout << "PremixedSolver::storeData(): step: " << step << " time: " << time << endl;

    int count = ncv*7; // rho,u(3),Z,C,rhoE TODO: add scalars
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
    FOR_ICV store_buf_recent[count++] = Z[icv];
    FOR_ICV store_buf_recent[count++] = C[icv];
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

    double my_bbminmax[6] = { 1.0E+20, 1.0E+20, 1.0E+20, 1.0E+20, 1.0E+20, 1.0E+20 }; 
    FOR_ICV {
      if ((rho[icv] != rho[icv])||(p[icv] != p[icv])||(T[icv] != T[icv])) {
	FOR_I3 my_bbminmax[i] = min(my_bbminmax[i],x_cv[icv][i]);
	FOR_I3 my_bbminmax[3+i] = min(my_bbminmax[3+i],-x_cv[icv][i]);
      }
    }
    
    double bbminmax[6];
    MPI_Reduce(my_bbminmax,bbminmax,6,MPI_DOUBLE,MPI_MIN,0,mpi_comm);

    if (mpi_rank == 0) {
      cout << "WARNING: NaNs detected in step: " << step << " time: " << time << " bbmin: " << 
	bbminmax[0] << " " << bbminmax[1] << " " << bbminmax[2] << " bbmax: " << 
	-bbminmax[3] << " " << -bbminmax[4] << " " << -bbminmax[5] << 
        "\n > restoring solution from step: " << store_step << " time: " << store_time <<
	"\n > adjusting TIMESTEP control from: " << dt_data[0] << " to: " << dt_data[0]*0.8 << endl;
    }

    step = store_step;
    time = store_time;
    int count = 0;
    assert(store_buf);
    FOR_ICV rho[icv] = store_buf[count++];
    FOR_ICV FOR_I3 u[icv][i] = store_buf[count++];
    FOR_ICV Z[icv] = store_buf[count++];
    FOR_ICV C[icv] = store_buf[count++];
    FOR_ICV rhoE[icv] = store_buf[count++];
    FOR_BCZONE count += (*it)->restoreData(store_buf+count);
    assert(count == store_buf_size);

    // after restoring data, everything needs to be re-computed...

    updateConservativeAndPrimitiveData();
    processIgnition();
    calcSgsAndMaterialProperties();

    // and reduce the timestep by some factor...
    // for either CFL or DT mode, the timestep is in FlowSolver::dt_data[0]...

    dt_data[0] *= 0.8;

    // don't let this go on forever...
    if (dt_data[0] < 1.0E-30) {
      CERR("TIMESTEP control must have failed. Exiting");
    }
    
  }

  int advanceSolution() {
    
    // if using cht (conjugate heat transfer), we should come into this routine
    // with the cht temperature already updated into the flow side CHT bc's:
    // see initComplete() above.

    for (int icv = 0; icv < ncv; ++icv) {
      rhs0[icv].zero();
      rhs1[icv].zero();
      rhs2[icv].zero();
    }

    if ( rhs0_sc) {
      for (int ii = 0; ii < nsc_transport*ncv; ++ii) {
        rhs0_sc[ii] = 0.0;
        rhs1_sc[ii] = 0.0;
        rhs2_sc[ii] = 0.0;
      }
    }

    // rk3 update...
    FOR_BCZONE (*it)->calcRhs(time,1);
    calcRhs(rhs0,rhs0_sc,time,1);
    if (cht) cht->updateQStart();
    timer.split("calc rhs 1");

    FOR_BCZONE (*it)->rk3Step(dt,1);
    rk3Step(erk_wgt3[0],dt,1);
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
    timer.split("calc rhs 2");

    FOR_BCZONE (*it)->rk3Step(dt,2);
    rk3Step(erk_wgt3[1],dt,2);
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
    timer.split("calc rhs 3");

    FOR_BCZONE (*it)->rk3Step(dt,3);
    rk3Step(erk_wgt3[2],dt,3);
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

    // finite time interval duration ignition kernels...
    processIgnition();
    timer.split("ignite");

    // ensure consistency of the sgs, etc with the new state
    calcSgsAndMaterialProperties();
    timer.split("calc sgs and properties");
    
    if (checkForNans()) {
      restoreData();
      return -1;
    }

    report();
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
  
  void queryBcs();

  void report() {

    if (step%check_interval == 0) {
      
      // store the data every so often...
      if (step%(check_interval*20) == 0) 
        storeData();
      
      // variable ranges...
      dumpRange(rho,ncv,"rho");
      dumpRange(u,ncv,"u");
      dumpRange(p,ncv,"p");
      dumpRange(T,ncv, "T");
      dumpRange(Z,ncv, "Z");
      dumpRange(C,ncv, "C");
      dumpRange(e_fa,nfa, "e_fa");
      dumpRange(mu_sgs,ncv,"mu_sgs");
      dumpRange(a_sgs,ncv,"a_sgs");
      dumpCfl();
      dumpRangeScalars();

      // global conservation...
      double my_buf[7] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
      FOR_ICV {
	my_buf[0] += rho[icv]*vol_cv[icv];
	my_buf[1] += rho[icv]*u[icv][0]*vol_cv[icv];
	my_buf[2] += rho[icv]*u[icv][1]*vol_cv[icv];
	my_buf[3] += rho[icv]*u[icv][2]*vol_cv[icv];
	my_buf[4] += rhoE[icv]*vol_cv[icv];
	my_buf[5] += rho[icv]*Z[icv]*vol_cv[icv];
	my_buf[6] += rho[icv]*C[icv]*vol_cv[icv];
      }
      double buf[7];
      MPI_Reduce(my_buf,buf,7,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      if (mpi_rank == 0) {
	cout << " > time: " << time << " total mass, momentum, energy, scalar(s): " <<
	  buf[0] << " " << buf[1] << " " << buf[2] << " " << buf[3] << " " << buf[4] << " " << buf[5] << " " << buf[6] << endl;
      }

    }

    queryBcs();

  }

  void setFgrCvs();
  void setFgrCvsFull();

  void processIgnition();

  void calcRhs(PremixedRhs *__restrict__ rhs, double *__restrict__ rhs_sc, const double time, const int rk_stage);
  void calcRhsScalars( double* __restrict__ rhs_sc, const double* __restrict__ mf,
                       const double time, const int rk_stage);

  virtual void addSourceHook(PremixedRhs* rhs, const double time, const int rk_stage) {}

  void addBodyForceSourceTerms(PremixedRhs* __restrict__ rhs, double *__restrict__ rhs_sc,
                               const double time, const int rk_stage) ;

  void rk3Step(const double rk_wgt[3], const double dt, const int rk_stage);
  void advanceScalars(const double * rho_old, const double rk_wgt[3], const double time, const int rk_stage);

  void calcSgsAndMaterialProperties();

  void updateExtendedState(PremixedState *__restrict__ cv_light, AuxGasProp *__restrict__ cv_compact,
    double *__restrict__ p, double *__restrict__ T,
    double *__restrict__ h_prime, double *__restrict__ s_prime,
    const double *__restrict__ rho, const double (*__restrict__ u)[3],
    const double *__restrict__ rhoE, const double *__restrict__ Z,
    const double *__restrict__ C, const double* __restrict__ R,
    const double *__restrict__ T_p, const double* __restrict__ e_p,
    const double *__restrict__ gamma, const double * __restrict__ a_gam,
    const double *__restrict__ ent, const int icv_start, const int icv_end);


  void updateConservativeAndPrimitiveData(const int action = COMPLETE_ALL_UPDATES) {

    // start the communication of the rho,u,rhoE,Z,C...

    updateCv2ClassStart<PremixedComm,double>(comm_container, REPLACE_ROTATE_DATA);

    if ( !b_nox ) {

      // lookup the gas properties from the chemtable.

      chemtable->lookup(R,"R",T_p,"T",e_p,"e",gamma,"gamma",a_gam,"a_gamma",
                        ent,"s",Z,C,ncv);

    } else {

      // need to lookup the nox_thermal source term in addition -- note
      // that the nox source term is not needed in the ghosts, so it is
      // only done here and not below.

      chemtable->lookup(R,"R",T_p,"T",e_p,"e",gamma,"gamma",a_gam,"a_gamma",ent,"s",
                        scalar_src_map["nox_thermal_src"],"src_nox_therm",Z,C,ncv);

    }

    // the following routines takes all of the variables as arguments to
    // the function, even though they are members of the class so that
    // the pointers are be "__restrict__"'ed; this is important for vectorization
    // supported in recent compilers...

    updateExtendedState(cv_light,cv_compact,p,T,h_prime,s_prime,
                        rho,u,rhoE,Z,C,R,T_p,e_p,gamma,a_gam,ent,0,ncv);


    // start the gradient calc for purely internal cells...

    calcCvGradStart(dudx,u);

    // lookup the source term at the face.  in the previous version of the
    // solver, Z_fa, C0_fa, C1_fa were populated; but this seems inefficient.
    // instead we will pass along the cvofa for the face-based lookup..

    if ( !b_nox) {

      chemtable->lookupSpecialNew(int_rho_csrc,"int_rho_src",Z,C,cvofa,0,nfa_i);

    } else {

      // batch the lookup of the integrated nox source term with the prog var src

      chemtable->lookupSpecialNew(int_rho_csrc,"int_rho_src",
                                  scalar_src_map["nox_prompt_src"], "int_src_nox_prompt",
                                  Z,C,cvofa,0,nfa_i);

    }

    // complete the communication of the ghosts and repeat the above in the ghosts..

    updateCv2ClassFinish<PremixedComm,double>(comm_container);

    // we need these to set the full thermodynamic state in the full ghosts..

    chemtable->lookup(&R[ncv],"R",&T_p[ncv],"T",&e_p[ncv],"e",&gamma[ncv],"gamma",
                      &a_gam[ncv],"a_gamma",&ent[ncv],"s",&Z[ncv],&C[ncv],ncv_g2-ncv);

    // state calculation in the ghosts..

    updateExtendedState(cv_light,cv_compact,p,T,h_prime,s_prime,
                        rho,u,rhoE,Z,C,R,T_p,e_p,gamma,a_gam,ent,ncv,ncv_g2);

    // complete the gradient calculation in non-purely internal cvs..

    calcCvGradFinish(dudx,u);

    // ghost source term lookups [note that we are not looking up any source terms
    // at the face under the assumption that dZ/dn, dC/dn = 0 at the boundaries

    if ( !b_nox) {

      chemtable->lookupSpecialNew(int_rho_csrc,"int_rho_src",Z,C,cvofa,nfa_i,nfa);

    } else {

      chemtable->lookupSpecialNew(int_rho_csrc,"int_rho_src",
                                  scalar_src_map["nox_prompt_src"], "int_src_nox_prompt",
                                  Z,C,cvofa,nfa_i,nfa);

    }

    // lastly update scalars' ghosts -- if there are many scalars we can consider some latency hiding here.

    updateScalars();


    // communicate the grad(u) in the ghosts.  if requested, we can force all communication
    // to complete here.  we could alternative postpone this to hide a bit of latency while
    // the grad(u) communication completes.  to be handled in the futre

    updateCvDataStart(dudx, REPLACE_ROTATE_DATA);
    updateCvDataFinish(dudx);
    // if ( action == COMPLETE_ALL_UPDATES)
         //updateCvDataFinish(dudx);

    if ( checkParam("NO_SOURCE_TERM") ) {

      for (int ifa = 0; ifa < nfa; ++ifa)
        int_rho_csrc[ifa] = 0.0;

    }
    else if (Param * param = getParam("SCALE_SOURCE_TERM")) {

      const double factor = param->getDouble();
      for (int ifa = 0; ifa < nfa; ++ifa)
        int_rho_csrc[ifa] *= factor;

    }

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
      bc_map[iter->first]->force(f_buf,m_buf,chemtable->pressure);
      FOR_I3 FOR_J3 forceZoneBuf[iter->second][3*i+j] = f_buf[i][j];
      FOR_I3 FOR_J3 momentZoneBuf[iter->second][3*i+j] = m_buf[i][j];
    }
  }

};
#undef FOR_BCZONE
#endif
