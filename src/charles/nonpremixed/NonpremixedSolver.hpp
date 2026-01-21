#ifndef NONPREMIXEDSOLVER_HPP
#define NONPREMIXEDSOLVER_HPP

#include "FlowSolver.hpp"
#include "PremixedSolverFlux.hpp"
#include "NonpremixedSolverBcs.hpp"
#include "RkWeights.hpp"
#include "CommContainers.hpp"
#include "Chemtable.hpp"
#include "CtiLiquid.hpp" // XXX Wut?

using namespace RkWeights;

#define FOR_BCZONE for(vector<NonpremixedBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)

class SprayLp {
public:

  double x[3];
  double u[3];
  double T;
  double mass;
  double rho;

};

class SprayContainer {
public:

  int nsp;
  int capacity;

  SprayLp* lsp_vec;
  CtiLiquid* properties;

  SprayContainer() : nsp(0), capacity(0), lsp_vec(NULL), properties(NULL) {}

  void registerData() {
    assert(0);
  }

};

class NonpremixedSolver : public FlowSolver {
public:

  double * rho;           // density
  double (*u)[3];         // velocity
  double * rhoE;          // total energy
  double * p;             // pressure
  double * T;             // temperature
  double * sos;           // speed of sound..

  // scalars used to describe the thermochemical state...
  double * Z;             // mixture fraction
  double * Zvar;          // mixture fraction variance
  double * C;             // progress variable
  double * csrc;          // prog var source term..

  // off manifold extensions
  bool offManifold;       // Flag to in-/ex- clude compressibility/heat transfer extensions
  double * ratio;         // Source term overall enhancement ratio
  double * Tratio;        // Source term temperature enhancement ratio
  double * pratio;        // Source term pressure enhancement ratio
  double Tratio_max;      // Extension limiting
  double Tratio_min;      // Extension limiting
  double pratio_max;      // Extension limiting
  double pratio_min;      // Extension limiting

  double (*dudx)[3][3];   // velocity gradient

  double * mu_lam;        // laminar visc
  double * mu_sgs;        // sgs eddy visc
  double * loc_lam;       // lambda on cp (molecular)
  double * loc_sgs;       // .. sgs contribution
  double * a_sgs;         // acoustic sgs (density perturbation sgs)

  double * R;             // gas constant
  double * e_p;           // internal energy at the chemtable pressure
  double * gamma;         // specific heat ratio
  double * a_gam;         // sp heat ratio linear expansion coeff
  double * T_p;           // nominal temperature at the chemtable pressure

  // other parameters to parse ..
  string sgs_model;

  // spray particles -- assuming for now that we have a single spray
  // container (i.e., injecting a single fuel) -- could keep a vector
  // of these containers if a fuel blend was used.
  SprayContainer spray;

  // packed arrays for cache efficient flux calculations ..
  NonpremixedState* cv_light;
  AuxGasProp* cv_compact;

  // Packed communicator patterns...
  NonpremixedComm* comm_container;
  SgsComm* sgs_comm_container;

   // Rk3 update uses 3 rhs..
  NonpremixedRhs* rhs0;
  NonpremixedRhs* rhs1;
  NonpremixedRhs* rhs2;

  AbstractChemtable3D* chemtable;

  // Boundary conditions
  vector<NonpremixedBc*> bcs;
  map<string,NonpremixedBc*> bc_map;
  map<string,pair<int,NonpremixedBc*> > bc_queries;

  // mass flux stored on extended faces -- used to help transport scalars..
  double * mf;
  double * rhs0_sc;
  double * rhs1_sc;
  double * rhs2_sc;
  double * rho_old;


  NonpremixedSolver() {

    rho     = NULL; registerCvData(rho , "rho" , READWRITE_DATA);
    u       = NULL; registerCvData(u   , "u"   , READWRITE_DATA);
    rhoE    = NULL; registerCvData(rhoE, "rhoE", READWRITE_DATA);
    Z       = NULL; registerCvData(Z   , "Z"   , READWRITE_DATA);
    Zvar    = NULL; registerCvData(Zvar, "Zvar", READWRITE_DATA);
    C       = NULL; registerCvData(C   , "C"   , READWRITE_DATA);
    p       = NULL; registerCvData(p   , "p"   , READWRITE_DATA);
    T       = NULL; registerCvData(T   , "T"   , READWRITE_DATA);
    sos     = NULL; registerCvData(sos , "sos" , READWRITE_DATA);

    mu_lam  = NULL; registerCvData(mu_lam , "mu_lam" ,  CAN_WRITE_DATA);
    mu_sgs  = NULL; registerCvData(mu_sgs , "mu_sgs" ,  CAN_WRITE_DATA);
    loc_lam = NULL; registerCvData(loc_lam, "loc_lam" , CAN_WRITE_DATA);
    loc_sgs = NULL; registerCvData(loc_sgs, "loc_sgs" , CAN_WRITE_DATA);
    a_sgs   = NULL; registerCvData(a_sgs  , "a_sgs" ,   CAN_WRITE_DATA);
    csrc    = NULL; registerCvData(csrc,    "csrc",     CAN_WRITE_DATA);

    // Off manifold extensions
    ratio = NULL; registerCvData(ratio, "ratio",        CAN_WRITE_DATA);
    Tratio = NULL; registerCvData(Tratio, "Tratio",     CAN_WRITE_DATA);
    pratio = NULL; registerCvData(pratio, "pratio",     CAN_WRITE_DATA);

    R       = NULL; registerCvData(R, "R",              CAN_WRITE_DATA);
    e_p     = NULL; registerCvData(e_p, "e_p",          CAN_WRITE_DATA);
    gamma   = NULL; registerCvData(gamma, "gamma",      CAN_WRITE_DATA);
    a_gam   = NULL; registerCvData(a_gam, "a_gam",      CAN_WRITE_DATA);
    T_p     = NULL; registerCvData(T_p, "T_p",          CAN_WRITE_DATA);

    cv_light           = NULL;
    cv_compact         = NULL;

    comm_container     = NULL;
    sgs_comm_container = NULL;

    rhs0               = NULL;
    rhs1               = NULL;
    rhs2               = NULL;

    chemtable          = NULL;

    mf                 = NULL;
    rhs0_sc            = NULL;
    rhs1_sc            = NULL;
    rhs2_sc            = NULL;
    rho_old            = NULL;

    dudx               = NULL;


    // handle the spray registration if it is present.  you should have
    // access to the params by this stage to decide if there is spray

    if (checkParam("LIQUID_FUEL")) {

      // syntax of the above call to be determined ... XXX

      assert ( spray.lsp_vec == NULL);
      spray.registerData();

    }

    parseSgsModel(sgs_model);

    if (mpi_rank == 0)
      cout << " > sgs model: " << sgs_model << endl;

  }//NonpremixedSolver()

  ~NonpremixedSolver();

  //==========================================
  // initialization functions
  //=========================================

  void init() {
    initNonpremixedChemistry();
    FlowSolver::init();
  }

  void initMin() {
    initNonpremixedChemistry();
    FlowSolver::initMin();
  }

  void initNonpremixedChemistry() {
    initChemtable(chemtable,getStringParam("CHEMTABLE"));

    offManifold = getBoolParam("OFF-MANIFOLD", false);

    if (mpi_rank == 0) {
      if (offManifold)
        cout << "OFF-MANIFOLD flag set to TRUE.\n"
             << "\tScaling law extensions will be applied to the combustion model\n"
             << "\tto account for effects such as heat transfer and compressibility\n"
             << "\tthat could potentially force the fluid state to deviate from the\n"
             << "\tbasline flamelet model\n";
      else
        cout << "OFF-MANIFOLD flag set to FALSE.\n"
             << "\tThis is appropriate for low Mach (Ma < 0.3) reacting flows with\n"
             << "\tminimal boundary heat transfer. If this describes your flow,\n"
             << "\trunning this way will provide a (modest) performance boost.\n"
             << "\tIf this does NOT desccribe your flow, please set\n\n"
             << "\t\tOFF-MANIFOLD = TRUE\n\n"
             << "\tin your input file and restart the job.\n";
    }

    // Get the source term enhancement caps
    if (offManifold) {
      Tratio_max = getDoubleParam("TRATIO_MAX", 10.0);
      Tratio_min = getDoubleParam("TRATIO_MIN", 0.1);
      pratio_max = getDoubleParam("PRATIO_MAX", 10.0);
      pratio_min = getDoubleParam("PRATIO_MIN", 0.1);
      if (mpi_rank == 0)
        cout << "OFF-MANIFOLD: Capping source extensions:\n"
             << "\t\t\tT_enhancement: [" << Tratio_min << ", " << Tratio_max << "]\n"
             << "\t\t\tp_enhancement: [" << pratio_min << ", " << pratio_max << "]\n";
    }

    // Identify the chemtable variables to load
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
    if (offManifold) {
      strVec.push_back("src_pos");
      strVec.push_back("src_neg");
      strVec.push_back("E_pos");
      strVec.push_back("E_neg");
      strVec.push_back("n_pos");
      strVec.push_back("n_neg");
    }

    chemtable->loadVariables(strVec);
  }//initNonpremixedChemistry()

  void initFromParams();
  void initFromPotentialFlow(const bool coldflow);

  virtual void initData();
  virtual void loadBalanceHook(StripedMesh* sm);

  void registerBoundaryConditions() {

    assert(bcs.size() == 0);

    StaticSolver::registerBoundaryConditions(); // registers bf geometric data

    int nerr = 0;
    vector<pair<string,string> > errors;

    FOR_IZONE(bfZoneVec) {

      const string zone_name = bfZoneVec[izone].getName();
      if (Param* p = getParam(zone_name)) {
        const string bc_type = p->getString(0);
        if ((bc_type == "SLIP")||(bc_type == "SYMMETRY"))
          bcs.push_back(new SlipWallNBc(&bfZoneVec[izone],this));
        else if (bc_type == "WALL_ADIABATIC")
          bcs.push_back(new WallAdiabaticN(&bfZoneVec[izone],this));
        else if (bc_type == "WALL_ISOTHERMAL")
          bcs.push_back(new WallIsothermalN(&bfZoneVec[izone],this));
        else if (bc_type == "WALL_CONDUCTIVE")
          bcs.push_back(new WallConductiveN(&bfZoneVec[izone], this));
        else if (bc_type == "WALL_CHT")
          bcs.push_back(new WallChtN(&bfZoneVec[izone], this));
        else if (bc_type == "CBC_RUPS")
          bcs.push_back(new CbcRupsN(&bfZoneVec[izone], this));
        else if (bc_type == "CBC_RUNPS")
          bcs.push_back(new CbcRunpsN(&bfZoneVec[izone], this));
        else if (bc_type == "CBC_MPTS")
          bcs.push_back(new CbcMptsN(&bfZoneVec[izone], this));
        else if (bc_type == "CBC_PROFILE")
          bcs.push_back(new CbcProfileN(&bfZoneVec[izone], this));
        else if (bc_type == "CBC_UPTS")
          bcs.push_back(new CbcUptsN(&bfZoneVec[izone], this));
        else if (bc_type == "CBC_UPTS_SPONGE")
          bcs.push_back(new CbcUptsSpongeN(&bfZoneVec[izone], this));
        else if (bc_type == "CBC_TOTAL_PTS")
          bcs.push_back(new CbcTotalPtsN(&bfZoneVec[izone], this));
        else if (bc_type == "NSCBC_MTS")
          bcs.push_back(new NscbcMtsN(&bfZoneVec[izone], this));
        else if ((bc_type == "NSCBC_OUTLET_P")||(bc_type == "NSCBC_OUTLET_PRESSURE")) //
          bcs.push_back(new NscbcOutletPN(&bfZoneVec[izone], this));
        else if (bc_type == "NSCBC_OUTLET_MDOT")
          bcs.push_back(new NscbcOutletMdotN(&bfZoneVec[izone], this));
        else if (bc_type == "NSCBC_MTS_SOFT")
          bcs.push_back(new NscbcMtsSoftN(&bfZoneVec[izone], this));
        else if (bc_type == "SPONGE")
          bcs.push_back(new SpongeN(&bfZoneVec[izone], this));
        else if (bc_type == "WM_ALG_ADIABATIC")
          bcs.push_back(new WmAlgAdiabaticN(&bfZoneVec[izone], this));
        else if (bc_type == "WM_ALG_ISOTHERMAL")
          bcs.push_back(new WmAlgIsothermalN(&bfZoneVec[izone], this));
        else if (bc_type == "WM_ALG_CHT")
          bcs.push_back(new WmAlgChtN(&bfZoneVec[izone], this));
        else if (bc_type == "WM_ALG_CONDUCTIVE")
          bcs.push_back(new WmAlgConductiveN(&bfZoneVec[izone], this));
        else if (bc_type == "INFLOW_TURB")
          bcs.push_back(new InflowTurbulenceN(&bfZoneVec[izone], this));
        else if (bc_type == "HOOK")
          bcs.push_back(initHookBc(&bfZoneVec[izone]));
        else {
          nerr++;
          errors.push_back(pair<string,string>(zone_name,bc_type));
        }
      }
      else {
        nerr++;
        errors.push_back(pair<string,string>(zone_name,""));
      }
    }

    // ensure that the error handling is synced..; should be uncessary and checked below

    MPI_Bcast(&nerr,1,MPI_INT,0,mpi_comm);
    assert(nerr == int(errors.size()));
    reportBcErrors(errors);

  }//registerBoundaryConditions()

  void initBoundaryConditions() {

    // allow the boundary conditions to allocate their necessary data (not allowed to
    // set anything here, because it may be coming from the restart data yet to be read)..

    FOR_BCZONE (*it)->initData();

    // create the name based map of the bc objects for easy query functions late..

    FOR_BCZONE bc_map[(*it)->getName()] = *it;

  }//initBoundaryConditions()

  virtual NonpremixedBc* initHookBc(BfZone* p) {
    CERR(" > user must supply a hook boundary condition for zone: " << p->getName());
  }

  void initComplete() {
    // initial hooks and init from params have both been called at this juncture.
    updateConservativeAndPrimitiveData();
    calcSgsAndMaterialProperties();

    // the initial temperature in the fem/cht solver needs to be sent
    // to the flow boundary conditions for the first fluid side flux
    // calc. Here we use the blocking version since it is only called
    // at the start...

    if (cht) cht->updateT();
  }//initComplete()

  void initialHookBcs() {
    FOR_BCZONE (*it)->initialHook();

    // this is the last routine called before the run loop, so use it as an opportunity
    // to store the current conditions. In this way, we can use robust timestepping
    // from the very beginning...
    storeData();
  }

  NonpremixedBc* getBc(const string& name) const {
    map<string,NonpremixedBc*>::const_iterator it = bc_map.find(name);
    if (it == bc_map.end())
      return NULL;
    else
      return it->second;
  }

  //=========================================
  // Sgs and combustion models
  //=========================================
  void parseSgsModel(string& sgs_model) {
    const string sgs_type = getStringParam("SGS_MODEL", "NONE");
    if (sgs_type == "NONE")
      sgs_model = "NONE";
    else if (sgs_type == "VREMAN")
      sgs_model = "VREMAN";
    else
      CERR(" > unrecognized sgs model: " << sgs_type);

  }

  void computeSgsNone();
  void computeSgsVreman();
  void limitAsgs();


  void syncPostState() {

    updateConservativeAndPrimitiveData();
    calcSgsAndMaterialProperties();

    // Boundary conditions are not queried as part of processStep -- explicit call here.

    queryBcs();

  }

  void storeData() {

    if (mpi_rank == 0) cout << "NonpremixedSolver::storeData(): step: " << step << " time: " << time << endl;

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

    if (rhs0_sc) {
      for (int ii = 0; ii < nsc_transport*ncv; ++ii) {
        rhs0_sc[ii] = 0.0;
        rhs1_sc[ii] = 0.0;
        rhs2_sc[ii] = 0.0;
      }
    }

    // rk3 stage 1...
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

    // rk3 stage 2..
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

    // rk3 stage 3..
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

      // var ranges
      dumpRange(rho,ncv,"rho");
      dumpRange(u,ncv,"u");
      dumpRange(p,ncv,"p");
      dumpRange(T,ncv, "T");
      dumpRange(Z,ncv, "Z");
      dumpRange(Zvar,ncv, "Zvar");
      dumpRange(C,ncv, "C");
      dumpRange(csrc,ncv, "csrc");
      dumpRange(mu_sgs,ncv,"mu_sgs");
      dumpRange(a_sgs,ncv,"a_sgs");
      dumpCfl();
      dumpRangeScalars();

      // global conservation...
      double my_buf[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

      for (int icv = 0; icv < ncv; ++icv) {

        my_buf[0] += rho[icv]*vol_cv[icv];
      	my_buf[1] += rho[icv]*u[icv][0]*vol_cv[icv];
      	my_buf[2] += rho[icv]*u[icv][1]*vol_cv[icv];
      	my_buf[3] += rho[icv]*u[icv][2]*vol_cv[icv];
      	my_buf[4] += rhoE[icv]*vol_cv[icv];
      	my_buf[5] += rho[icv]*Z[icv]*vol_cv[icv];
      	my_buf[6] += rho[icv]*Zvar[icv]*vol_cv[icv];
      	my_buf[7] += rho[icv]*C[icv]*vol_cv[icv];
        my_buf[8] += p[icv]*vol_cv[icv];
        my_buf[9] += vol_cv[icv];

      }

      double buf[10];
      MPI_Reduce(my_buf,buf,10,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

      if (mpi_rank == 0)
        cout  << " > Global tots:"   << setprecision(8)
              << " time = "          << time
              << ", mass = "         << buf[0]
              << ", x-momentum = "   << buf[1]
              << ", y-momentum = "   << buf[2]
              << ", z-momentum = "   << buf[3]
              << ", energy = "       << buf[4]
              << ", Z = "            << buf[5]
              << ", Zvar = "         << buf[6]
              << ", C = "            << buf[7]
              << endl;
      if (mpi_rank == 0)
        cout  << " > Global avgs:"    << setprecision(8)
              << " time = "           << time
              << ", rho = "           << buf[0]/buf[9]  // Rey
              << ", u-velocity = "    << buf[1]/buf[0]  // Favre
              << ", v-velocity = "    << buf[2]/buf[0]  // Favre
              << ", w-velocity = "    << buf[3]/buf[0]  // Favre
              << ", energy = "        << buf[4]/buf[0]  // Favre
              << ", Z = "             << buf[5]/buf[0]  // Favre
              << ", Zvar = "          << buf[6]/buf[0]  // Favre
              << ", C = "             << buf[7]/buf[0]  // Favre
              << ", p = "             << buf[8]/buf[9]
              << endl;
    }

    queryBcs();

  }

  void setFgrCvs();
  void processIgnition();

  void calcRhs(NonpremixedRhs *__restrict__ rhs, double *__restrict__ rhs_sc, const double time, const int rk_stage);
  void calcRhsScalars(double* __restrict__ rhs_sc, const double* __restrict__ mf,
                      const double time, const int rk_stage);

  virtual void addSourceHook(NonpremixedRhs* rhs, const double time, const int rk_stage) {}

  void rk3Step(const double rk_wgt[3], const double dt, const int rk_stage);
  void advanceScalars(const double * rho_old, const double rk_wgt[3], const double time, const int rk_stage);

  void calcSgsAndMaterialProperties();

  void updateExtendedState(NonpremixedState* __restrict__ cv_light, AuxGasProp* __restrict__ cv_compact,
                           double *__restrict__ p,            double* __restrict__ T,
                           const double *__restrict__ rho,    const double (*__restrict__ u)[3],
                           const double *__restrict__ rhoE,   const double *__restrict__ Z,
                           const double *__restrict__ C,      const double* __restrict__ R,
                           const double *__restrict__ T_p,    const double* __restrict__ e_p,
                           const double *__restrict__ gamma,  const double * __restrict__ a_gam,
                           const int icv_start,  const int icv_end);



  void updateConservativeAndPrimitiveData(const int action = COMPLETE_ALL_UPDATES) {

    // post the send requests for the state and the mixture fraction immediately.

    updateCv2DataStart(Z);
    updateCv2ClassStart<NonpremixedComm,double>(comm_container, REPLACE_ROTATE_DATA);


    double * wgt = new double[ncv];

    for (int icv = 0; icv < ncv; ++icv) {

      wgt[icv]   = 0.0;
      Zvar[icv]  = 0.0;

    }

    for (int ifa = 0; ifa < nfa_i; ++ifa) {

      const int icv0  = cvofa[ifa][0];
      const int icv1  = cvofa[ifa][1];

      const double dZ    = Z[icv1] - Z[icv0];
      const double area  = MAG(n_fa[ifa]);
      const double delta = area/area_over_delta_fa[ifa];

      const double vol_     = delta*area;
      const double vol_dZ2  = vol_*dZ*dZ;
      wgt[icv0]            += vol_;
      wgt[icv1]            += vol_;
      Zvar[icv0]           += vol_dZ2;
      Zvar[icv1]           += vol_dZ2;

    }

    updateCv2DataFinish(Z);

    for (int ifa = nfa_i; ifa < nfa; ++ifa) {

      const int icv0  = cvofa[ifa][0];
      const int icv1  = cvofa[ifa][1];

      const double dZ    = Z[icv1] - Z[icv0];
      const double area  = MAG(n_fa[ifa]);
      const double delta = area/area_over_delta_fa[ifa];

      const double vol_     = delta*area;
      const double vol_dZ2  = vol_*dZ*dZ;
      wgt[icv0]            += vol_;
      Zvar[icv0]           += vol_dZ2;

    }

    // the following factor is derived from the integration
    // of the variance based on a linear field between icv0,icv1
    // contained in the pyramid formed by the two cell points and
    // the face... could be revisited.

    const double var_fax    = 1.0/40.0;

    for (int icv = 0; icv < ncv; ++icv) {

      Zvar[icv] /= wgt[icv];
      Zvar[icv] *= var_fax;

    }

    delete[] wgt;

    if (checkParam("NO_ZVAR")) {
      for (int icv = 0; icv < ncv; ++icv) Zvar[icv] = 0.0;
    }

    updateCv2DataStart(Zvar);


    if (offManifold) {
      double *src_pos = new double[ncv];
      double *src_neg = new double[ncv];
      double *E_pos   = new double[ncv];
      double *E_neg   = new double[ncv];
      double *n_pos   = new double[ncv];
      double *n_neg   = new double[ncv];
      chemtable->lookup(R, "R", T_p, "T", e_p, "e", gamma, "gamma", a_gam, "a_gamma",
                        src_pos, "src_pos", src_neg, "src_neg", E_pos, "E_pos", E_neg, "E_neg",
                        csrc, "src_prog", n_pos, "n_pos", n_neg, "n_neg",
                        Z, Zvar, C, ncv);

      const double eps = 1E-3;//regularization for enhancement ratios

      FOR_ICV {
        // The use of fmin and fmax here should protect against overflow errors
        const double csrc_0     = src_pos[icv] - src_neg[icv];
        const double Tpos_ratio = fmin(fmax(exp(-E_pos[icv]*(T_p[icv]/T[icv] - 1.0)), Tratio_min), Tratio_max);
        const double Tneg_ratio = fmin(fmax(exp(-E_neg[icv]*(T_p[icv]/T[icv] - 1.0)), Tratio_min), Tratio_max);
        const double ppos_ratio = fmin(fmax(pow(p[icv]/chemtable->pressure, n_pos[icv]), Tratio_min), Tratio_max);
        const double pneg_ratio = fmin(fmax(pow(p[icv]/chemtable->pressure, n_neg[icv]), Tratio_min), Tratio_max);
        csrc[icv]               = max(src_pos[icv]*Tpos_ratio*ppos_ratio - src_neg[icv]*Tneg_ratio*pneg_ratio, 0.0);
        Tratio[icv]             = (src_pos[icv]*Tpos_ratio - src_neg[icv]*Tneg_ratio + eps)/(csrc_0 + eps);
        pratio[icv]             = (src_pos[icv]*ppos_ratio - src_neg[icv]*pneg_ratio + eps)/(csrc_0 + eps);
        ratio[icv]              = (csrc[icv] + eps)/(csrc_0 + eps);
      }

      delete [] src_pos;
      delete [] src_neg;
      delete [] E_pos;
      delete [] E_neg;
      delete [] n_pos;
      delete [] n_neg;
    }
    else {
      // lookup gas properties for the processor internal cells ..
      chemtable->lookup(R,"R",T_p,"T",e_p,"e",gamma,"gamma",a_gam,"a_gamma",
                        csrc,"src_prog",Z,Zvar,C,ncv);
    }

    // Now that the internal compact cv data has been updated, the internal extend cv state can be updated
    updateExtendedState(cv_light,cv_compact,p,T,rho,u,rhoE,Z,C,R,T_p,e_p,gamma,a_gam,0,ncv);


    // ensure that the packed state and the variance are completed.
    updateCv2ClassFinish<NonpremixedComm,double>(comm_container);
    updateCv2DataFinish(Zvar);

    // complete the extended state in the ghosts.

    chemtable->lookup(&R[ncv],"R",&T_p[ncv],"T",&e_p[ncv],"e",&gamma[ncv],"gamma",
                      &a_gam[ncv],"a_gamma",&Z[ncv],&Zvar[ncv],&C[ncv],ncv_g2-ncv);

    updateExtendedState(cv_light,cv_compact,p,T,rho,u,rhoE,Z,C,R,T_p,e_p,gamma,a_gam,ncv,ncv_g2);

    // start the gradient calc for purely internal cells...

    calcCvGradStart(dudx,u);

    // complete the gradient calculation in non-purely internal cvs..

    calcCvGradFinish(dudx,u);

    // communicate the grad(u) in the ghosts.  if requested, we can force all communication
    // to complete here.  we could alternative postpone this to hide a bit of latency while
    // the grad(u) communication completes.  to be handled in the futre

    updateCvDataStart(dudx, REPLACE_ROTATE_DATA);

    // lastly update scalars' ghosts -- if there are many scalars we can consider some latency hiding here.

    updateScalars();

    // communicate the grad(u) in the ghosts.  if requested, we can force all communication
    // to complete here.  we could alternative postpone this to hide a bit of latency while
    // the grad(u) communication completes.  to be handled in the futre

    updateCvDataFinish(dudx);

    if (checkParam("NO_SOURCE_TERM")) {

      for (int icv = 0; icv < ncv; ++icv)
        csrc[icv] = 0.0;

    }

  }//updateConservativeAndPrimitiveData()

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
