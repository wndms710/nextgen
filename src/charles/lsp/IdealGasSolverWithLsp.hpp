#ifndef _IDEAL_GAS_SOLVER_WITH_LSP_HPP_
#define _IDEAL_GAS_SOLVER_WITH_LSP_HPP_

// ======================================================================================
// TODO
// - in rebinLsp, modify memory management to use new/delete, and MAX_PPC should be int?
// - lspSrc implementation
// - ensure mu,loc is laminar in alll routines
// - can we compute rho(T,p?) once for each substep, instead of as needed. similar for calcCp, 
//   anything else Tp del
// - can we linearize dxp/dt = u with u implicit?
// - latency hiding, faster bcs, nbrs? other? to reduce lsp relocation time -- will this allow 
//   lp[ip].icv and paocv_i/v update in the substep? 
// ======================================================================================

#include "LspInjector.hpp"
#include "Lsp.hpp"
#include "LspStats.hpp"
#include "Adt.hpp"
#include "CtiLiquid.hpp"
#include "CtiSolid.hpp"
//#include "SingleLspSolver.hpp"
#include "LspBcs.hpp"

// convenience macro to loop over the boundary conditions...
#define FOR_BCZONE for(vector<IdealGasBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)

//================================
// revisit:
// recycleParticles()
//================================

enum collisionModel {
  NO_MODEL,
  STOCHASTIC,
  DETERMINISTIC
};

class planeTrashClass {
public:
  double* xc;
  double* n;
  double r;
  string name;
  planeTrashClass(Param* param) {
    xc = new double [3];
    n = new double [3];
    int iarg = 0;
    name = param->getString(iarg++);
    assert(param->getString(iarg) == "X");
    ++iarg;
    FOR_J3 xc[j] = param->getDouble(iarg++);
    assert(param->getString(iarg) == "N");
    ++iarg;
    FOR_J3 n[j] = param->getDouble(iarg++);
    assert(param->getString(iarg) == "R");
    ++iarg;
    r = param->getDouble(iarg++);
    IF_RANK0 cout << " > planeTrash name: " << name << " xc: " << COUT_VEC(xc) << " n: " << COUT_VEC(n) << " r: " << r << endl;
  }
  ~planeTrashClass() {
    DELETE(xc);
    DELETE(n);
  }

  bool passPlane(double* xp0, double* xp) {
    bool pass = false;
    double dx0[3] = DIFF(xp0,xc);
    double dx[3] = DIFF(xp,xc);
    double d0 = DOT_PRODUCT(dx0,n);
    double d1 = DOT_PRODUCT(dx,n);
    if (d0 <= 0. && d1 > 0.) {
      double dxc[3];
      FOR_I3 dxc[i] = (-d0*dx[i] + d1*dx0[i])/(d1-d0);
      double mag_dxc = MAG(dxc);
      if (mag_dxc <= r) {
        pass = true;
      }
    }
    return pass;
  }

};

class IdealGasSolverWithLsp : public IdealGasSolver {

protected:

  LpClass<LspState> * lpCls;

public:
  // Bcs need access to these, but aren't children of the solver class
  CtiMaterial * material;
  CtiLiquid * fuel;
  CtiSolid  * dust;

protected:
  bool SOL_LP, LIQ_LP;

  //===========
  // injectors
  //============
  list<LspInjector*> injectorList;

  int np_inject;

  //=================
  // particle Cv link
  //=================
  int * paocv_i;        // particle of cv index
  int npaocv_max;       // updated during build of paocv_i/v

  // time integration method EXPLICIT/IMPLICIT
  string time_int;

  // source terms
  IdealGasRhs * lspSrc;

  //==================
  // Breakup model
  //==================
  int BreakupModel;
  double weber_cr;
  double core_br;

  double bu_k1; // k1 in breakup model
  double bu_k2; // k2 in breakup model

  //==================
  // lsp statistics
  //==================
  list<LspStats*> lspStatsList;
  list<LspPDF> lspPDFList;

  //==================
  // lsp agglemoration
  //==================
  bool lspAgglom;
  double MAX_PPC, MAX_NPAR;

  //==================
  // monitoring files
  //==================
  //FILE * fh;
  //ofstream fh;

  // length after which particles trashed
  double plane_trash_bool;
  list<planeTrashClass*> planeTrashList;

  // to keep track of the particle injected mass
  double time0;

  // gravity vector
  double grav[3];

  // npar target multiplier to control the number of droplets from breakup
  double npar_target_multiplier;

  // non ideal bouncing model
  bool b_irregCol;
  double resCoef; // restitution coefficient
  double mu_s;    // static friction coefficient
  double mu_k;    // kinetic friction coefficient

  // boundary conditions
  vector<lpBaseBc*> lpZoneBcVec;
  vector<pair<lpBaseBc*, int> > queryBcVec; // second: interval
  int lsp_stats_interval;

  // collision model
  bool b_colMod;
  collisionModel colMod;
  double resCoef_colMod;

  // concentration of particle (# of particles / gas vol)
  double * CPar;
  double * CPar2;
  double CPar_n_old;
  double vol_total;
  double time_stat;
  vector<string>  LpStatZoneNameVec;
  vector<BfZone*> LpStatZoneVec;
  //double * cv_wall_flag;

  bool show_histogram;
  string hist_file;
  int hist_interval;
  int hist_nbins;
  bool hist_verbose_count;
  bool hist_verbose_mass;

  bool EllipsoidDrag_b;

public:

  IdealGasSolverWithLsp() {

    COUT1("IdealGasSolverWithLsp()");

    // create the lp (not allocating at this stage)
    lpCls   = new LpClass<LspState>; // sets the lp to NULL, np and NP_MAX to 0
    // registrater lp data
    registerLp<LspState>(lpCls->lp,"lsp",lpCls->np);
    registerLpData<LspState>(lpCls->lp,lpCls->lp->icv ,"lsp:icv" ,NO_READWRITE_DATA);
    registerLpData<LspState>(lpCls->lp,lpCls->lp->flag,"lsp:flag",NO_READWRITE_DATA);
    registerLpData<LspState>(lpCls->lp,lpCls->lp->xp  ,"lsp:xp"  ,READWRITE_DATA);
    registerLpData<LspState>(lpCls->lp,lpCls->lp->up  ,"lsp:up"  ,READWRITE_DATA);
    registerLpData<LspState>(lpCls->lp,lpCls->lp->mp  ,"lsp:mp"  ,READWRITE_DATA);
    registerLpData<LspState>(lpCls->lp,lpCls->lp->Tp  ,"lsp:Tp"  ,READWRITE_DATA);
    registerLpData<LspState>(lpCls->lp,lpCls->lp->npar,"lsp:npar",READWRITE_DATA);
    registerLpData<LspState>(lpCls->lp,lpCls->lp->tbu ,"lsp:tbu" ,READWRITE_DATA);
    registerLpData<LspState>(lpCls->lp,lpCls->lp->dp  ,"lsp:dp"  ,READWRITE_DATA);
    registerLpData<LspState>(lpCls->lp,lpCls->lp->e   ,"lsp:e"   ,READWRITE_DATA);
    registerLpData<LspState>(lpCls->lp,lpCls->lp->f   ,"lsp:f"   ,READWRITE_DATA);

    paocv_i = NULL;
    npaocv_max = 0;

    CPar = NULL;
    CPar2 = NULL;
    registerCvData(CPar,"CPar",READWRITE_DATA);
    registerCvData(CPar2,"CPar2",READWRITE_DATA);

    //registerCvData(cv_wall_flag,"cv_wall_flag",READWRITE_DATA);
    time_stat = 0.0; registerData(time_stat,"time_stat",READWRITE_DATA);

    material= NULL;
    fuel    = NULL;
    dust    = NULL;
    lspSrc  = NULL;
    SOL_LP  = LIQ_LP = false;

    // hacking the fuel right now
    Param * param = getParam("LSP.MATERIAL");
    if (param == NULL) {
      CERR("could not find param LSP.MATERIAL. Possible choices are:\n" << 
           "LSP.MATERIAL LIQUID <liquid-tokens>\n" <<
           "LSP.MATERIAL SOLID <solid-tokens>");
    }
    int iarg = 0;
    string token = param->getString(iarg++);
    if (token == "LIQUID") {
      fuel = newCtiLiquid(param,iarg);
      material = fuel;
      LIQ_LP = true;
      IF_RANK0 fuel->info(300,350,10);
      //MPI_Pause("check it out");
    } 
    else if (token == "SOLID") {
      dust = newCtiSolid(param,iarg);
      IF_RANK0 dust->info(300,350,10);
      //MPI_Pause("check it out");
      material = dust;
      SOL_LP = true;
    } 
    else {
      CERR("unrecognized LSP.MATERIAL: " << token << ". Possible choices are:\n" << 
           "LSP.MATERIAL LIQUID <liquid-tokens>\n" <<
           "LSP.MATERIAL SOLID <solid-tokens>");
    }

    // ========================
    // set up the breakup model
    // ========================
    weber_cr = 0.0;
    param = getParam("LSP.BREAKUP_MODEL");
    if (!param) {
      BreakupModel = NO_BREAKUP;
    } else {
      const string token = param->getString();
      if (token == "NO_BREAKUP") {
        BreakupModel = NO_BREAKUP;
        IF_RANK0 cout << "BREAKUP_MODEL: " << param->getString() << endl;
      } else if (token == "SBM") {
        IF_RANK0 cout << "BREAKUP_MODEL: " << param->getString() << endl;
        BreakupModel = SBM_BREAKUP;
        // set core_br parameter
        core_br = 1.0/sqrt(3.0); // this was set to 1.76 in the legacy code, however the value in the literature is 1/sqrt(3)
        //core_br = 1.76;
        //core_br = 1.;
        int iarg = 1;
        while (iarg < param->size()) {
          string sbmtoken = param->getString(iarg++);
          if (sbmtoken == "WEBER_CRITICAL") {
            weber_cr = param->getDouble(iarg++);
          } else {
            IF_RANK0 {
              cerr << "\n\n\n=================================================================\n" <<
                      "Unrecognized SBM option: \"" << sbmtoken << "\" in BREAKUP_MODEL, options available are:\n" <<
                      "  ---  LSP.BREAKUP_MODEL NO_BREAKUP \n" <<
                      "  ---  LSP.BREAKUP_MODEL SBM WEBER_CRITICAL <double> \n" <<
                      "  ---  ...\n" <<
                      "=================================================================\n\n\n";
              assert(0);
            }
          }
        }
        if (weber_cr <= 0.0) {
          IF_RANK0 {
            cerr << "\n\n\n=================================================================\n" <<
                    "WEBER_CRITICAL not assigned, options available are:\n" <<
                    "  ---  LSP.BREAKUP_MODEL NO_BREAKUP \n" <<
                    "  ---  LSP.BREAKUP_MODEL SBM WEBER_CRITICAL <double> \n" <<
                    "  ---  ...\n" <<
                    "=================================================================\n\n\n";
            assert(0);
          }
        }
      } else {
        IF_RANK0 {
          cerr << "\n\n\n=================================================================\n" <<
                  "Unrecognized BREAKUP_MODEL: \"" << token << "\" , options available are:\n" <<
                  "  ---  LSP.BREAKUP_MODEL NO_BREAKUP \n" <<
                  "  ---  LSP.BREAKUP_MODEL SBM WEBER_CRITICAL <double> \n" <<
                  "  ---  ...\n" <<
                  "=================================================================\n\n\n";
          assert(0);
        }
      }
    }

    bu_k1 = 0.0;
    bu_k2 = 0.0;
    core_br = 0.0;
    if (BreakupModel != NO_BREAKUP) {
      Param * param = getParam("LSP.BREAKUP_PARAMS");
      int iarg = 0;
      string token = param->getString(iarg++);
      assert(token == "K1");
      bu_k1 = param->getDouble(iarg++);

      token = param->getString(iarg++);
      assert(token == "K2");
      bu_k2 = param->getDouble(iarg++);

      token = param->getString(iarg++);
      assert(token == "CORE_BR");
      core_br = param->getDouble(iarg++);
    }
    IF_RANK0 cout << "Breakup parameters: " << " bu_k1: " << bu_k1 << " , bu_k2: " << bu_k2 << " , core_br: " << core_br << endl;

    // ===========================
    // set up the time integration
    // ===========================
    param = getParam("LSP.TIME_INT");
    if (!param) {
      CERR("IdealGasSolverWithLsp requires param LSP.TIME_INT <string>. Possible choices include\n" << 
           "LSP.TIME_INT 1WAY_COUPLED\n" <<
           "LSP.TIME_INT 2WAY_COUPLED");
    }
    token = param->getString();
    if ((token == "1WAY_COUPLED") || (token == "2WAY_COUPLED")) {
      time_int = token;
    } 
    else {
      CERR("unrecognized LSP.TIME_INT: " << token << ". Possible choices include\n" << 
           "LSP.TIME_INT 1WAY_COUPLED\n" <<
           "LSP.TIME_INT 2WAY_COUPLED");
    }

    // ===========================
    // lsp agglemoration
    // ===========================
    param = getParam("LSP.AGGLOM");
    if(param != NULL) {
      lspAgglom = true;
      MAX_PPC  = param->getDouble(0);
      MAX_NPAR = param->getDouble(1);
    } else {
      lspAgglom = false;
      IF_RANK0 cout << " > WARNING: running without LSP.AGGLOM" << endl;
    }

    // ===========================
    // lsp plane trash
    // ===========================
    bool has_planeTrash = false;
    FOR_PARAM_MATCHING("LSP.PLANE_TRASH") {
      planeTrashList.push_back(new planeTrashClass(&(*param)));
      has_planeTrash = true;
    }
    plane_trash_bool = (planeTrashList.size() > 0);
    if (has_planeTrash == true) {
      IF_RANK0 cout << " > Total of " << planeTrashList.size() << " plane trashes" << endl;
    }
    else {
      IF_RANK0 cout << " > No plane trash is specified!" << endl;
    }

    // ===========================
    // lsp gravity vector
    // ===========================
    FOR_I3 grav[i] = 0.0;
    param = getParam("LSP.GRAVITY");
    if (param) {
      FOR_I3 grav[i] = param->getDouble(i);
      COUT1(" > Has gravity: " << grav[0] << ", " << grav[1] << ", " << grav[2]);
    } else {
      COUT1(" > No gravity");
    }

    // ===========================
    // npar target multiplier
    // ===========================
    npar_target_multiplier = getDoubleParam("LSP.NPAR_MULTIPLIER_BU", 1.0);
    COUT1(" > npar target multiplier: " << npar_target_multiplier);

    // ========================
    // non ideal bounce
    // ========================
    param = getParam("LSP.IRREGULAR_COLLISION");
    resCoef = 0;
    mu_s = 0;
    mu_k = 0;
    if (param == NULL) {
      b_irregCol = false;
    } else {
      b_irregCol = true;
      bool b_hasRes = false;
      bool b_hasMuS = false;
      bool b_hasMuK = false;
      int iarg = 0;
      while (iarg < param->size()) {
        string token = param->getString(iarg++);
        if (token == "RESTITUTION") {
          resCoef = param->getDouble(iarg++);
          b_hasRes = true;
        } else if (token == "MU_S") {
          mu_s = param->getDouble(iarg++);
          b_hasMuS = true;
        } else if (token == "MU_K") {
          mu_k = param->getDouble(iarg++);
          b_hasMuK = true;
        }
      }
      if ((b_hasRes && b_hasMuS && b_hasMuK) != true) {
        CERR("LSP.IRREGULAR_COLLISION parameters are incorrect!");
      }
    }
    COUT1("LSP.IRREGULAR_COLLISION: boolian: "<<b_irregCol<<" , resCoef: "<<resCoef<<" , mu_s: "<<mu_s<<" , mu_k: "<<mu_k);

    colMod = NO_MODEL;
    resCoef_colMod = 1.; // if not present assume 1
    b_colMod = false;
    param = getParam("LSP.COLLISION_MODEL");
    if (param != NULL) {
      string token = param->getString();
      if (token == "NONE") {
        b_colMod = false;
      } else if (token == "STOCHASTIC") {
        b_colMod = true;
        colMod = STOCHASTIC;
      } else if (token == "DETERMINISTIC") {
        b_colMod = true;
        colMod = DETERMINISTIC;
      } else {
        CERR("Unknown collision model: " << token);
      }

      int iarg = 1;
      while (iarg < param->size()) {
        string token = param->getString(iarg++);
        if (token == "RESTITUTION") {
          resCoef_colMod = param->getDouble(iarg++);
        }
      }
    }
    COUT1("LSP.COLLISION_MODEL: boolian: " <<b_colMod<<" , resCoef: "<<resCoef_colMod);

    param = getParam("LSP.LP_STAT_ZONES");
    if (param) {
      token = param->getString();
      MiscUtils::splitCsv(LpStatZoneNameVec,token);
      IF_RANK0 {
        cout << " > zones for LSP.LP_STAT_ZONES: ";
        for (int i = 0; i < LpStatZoneNameVec.size(); i++) {
          cout << LpStatZoneNameVec[i] << " ";
        }
        cout << endl;
      }
    }

    show_histogram = false;
    hist_file = "";
    hist_interval = -1;
    hist_nbins = 120; // default value for bin counts
    hist_verbose_count = false;
    hist_verbose_mass = false;
    param = getParam("LSP.SHOW_HISTOGRAM");
    if (param) {
      show_histogram = true;
      int iarg = 0;
      while (iarg < param->size()) {
        string token = param->getString(iarg++);
        if (token == "FILE") {
          hist_file = param->getString(iarg++);
        } else if (token == "INTERVAL") {
          hist_interval = param->getInt(iarg++);
          if ((hist_interval>0) and ((hist_interval%check_interval != 0) or (hist_interval<check_interval)) )
            CWARN("interval for particle diameter histogram is not a multiply of checking interval");
        } else if (token == "NBINS") {
          hist_nbins = param->getInt(iarg++);
        } else if (token == "V_COUNT") {
          hist_verbose_count = true;
        } else if (token == "V_MASS") {
          hist_verbose_mass = true;
        } else {
          CERR("LSP.SHOW_HISTOGRAM param not correctly specified");
        }
      }
    }

    np_inject = 0;

    EllipsoidDrag_b = getBoolParam("LSP.ELLIPSOID_DRAG",false);
    if ((mpi_rank==0)&&EllipsoidDrag_b)
      cout << " > ellipsoidal drag correlation is applied to all solid particles." << endl;

  }

  virtual ~IdealGasSolverWithLsp() {

    COUT1("~IdealGasSolverWithLsp()");
    //if (fh.is_open()) fh.close();

    DELETE(lspSrc);
    DELETE(paocv_i);
    if (material != NULL) {material = NULL;} // material is just a pointer into the base class of either liquid or solid
    if (fuel != NULL) {delete fuel; fuel = NULL;}
    if (dust != NULL) {delete dust; dust = NULL;}
    if (lpCls != NULL) {delete lpCls; lpCls = NULL;}

    for (list<LspInjector*>::iterator iter = injectorList.begin(); iter != injectorList.end(); ++iter) {
      delete (*iter);
    }

    for (list<planeTrashClass*>::iterator it = planeTrashList.begin(); it != planeTrashList.end(); ++it) {
      delete (*it);
    }

    for (int izone = 0; izone < lpZoneBcVec.size(); izone++) {
      delete lpZoneBcVec[izone];
    }

    for (list<LspStats*>::iterator iter = lspStatsList.begin(); iter != lspStatsList.end(); ++iter) {
      delete (*iter);
    }
    //DELETE(cv_wall_flag);
    DELETE(CPar);
    DELETE(CPar2);
    
  }

  void init() {

    COUT1("IdealGasSolverWithLsp::init()");

    FlowSolver::init();

    // needed by app to control contex menu
    added_data_fields.insert("particles");

    //================
    // lsp injectors
    //================
    setupInjectors();

    //================
    // lsp statistics
    //================
    initLspStats();

    // ==================
    // CvAdt construction
    // ==================
    if (cvAdt == NULL) buildCvAdt();

    // =======================
    // lsp boundary conditions
    // =======================
    initLspBc();

    for (int iname = 0; iname < LpStatZoneNameVec.size(); iname++) {
      map<const string, int>::iterator it = bfZoneNameMap.find(LpStatZoneNameVec[iname]);
      if (it != bfZoneNameMap.end()) {
        LpStatZoneVec.push_back(&bfZoneVec[it->second]);
      }
    }
    IF_RANK0 {
      cout << " > captured zones for LSP.LP_STAT_ZONES: ";
      for (int i = 0; i < LpStatZoneVec.size(); i++) {
        cout << LpStatZoneVec[i]->getName() << " ";
      }
      cout << endl;
    }
    //cv_wall_flag = new double [ncv];
    //FOR_ICV cv_wall_flag[icv] = -1.0;

    if (lpTracer!=NULL) {
      lpTracer->u = this->u;
      lpTracer->ncv = ncv;
    }
  }

  void initMin() {

    IdealGasSolver::initMin();

    // needed by app to control contex menu
    added_data_fields.insert("particles");
  }

  void initData() {

    COUT1("IdealGasSolverWithLsp::initData()");
    IdealGasSolver::initData();

    //===================
    // should now have the correct np form restart
    // allocate the lp inside lpCls
    //===================
    lpCls->init();

    //===================
    // check if ellipsoidal parameters exist
    //===================
    if ( (!checkDataFlag("lsp:e")) || (!checkDataFlag("lsp:f")) ) {
      for (int ip = 0; ip < lpCls->size(); ip++) {
        lpCls->lp[ip].e = 1.;
        lpCls->lp[ip].f = 1.;
      }
    }

    //=================
    // particle Cv link
    //=================
    paocv_i = new int [ncv+1];

    // only used for mass trasfer
    lspSrc = new IdealGasRhs[ncv];
    FOR_ICV lspSrc[icv].zero();

    // particle concentration
    assert(CPar == NULL);  CPar = new double[ncv];
    assert(CPar2 == NULL); CPar2 = new double[ncv];
    double myvol = 0.0;
    for (int icv = 0; icv < ncv; icv++) {
      myvol += vol_cv[icv];
    }
    MPI_Allreduce(&myvol,&vol_total,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
    COUT1(" > vol_total: "<<vol_total);

    CPar_n_old = double(getLspSize()) / vol_total;

    //cout << "lpCls->size(): " << lpCls->size() << " , lpCls->sizeMax(): " << lpCls->sizeMax() << endl;
    //lpCls->init();
    //cout << "lpCls->size(): " << lpCls->size() << " , lpCls->sizeMax(): " << lpCls->sizeMax() << endl;
    //MPI_Pause("check np");

  }

  virtual void initialHook() {
    COUT1("IdealGasSolverWithLsp initialHook()");

    if (checkDataFlag("lsp:xp")) { // the lsp data comes from restart file
      IF_RANK0 cout << " > lp data is read from the restart" << endl;
      initHookLspBasic();
    }

    if ((!checkDataFlag("CPar"))||(!checkDataFlag("CPar2"))) {
      IF_RANK0 cout << " > initialize CPar and CPar2 to zero" << endl;
      for (int icv = 0; icv < ncv; icv++) {
        CPar[icv] = 0.0;
        CPar2[icv] = 0.0;
      }
    }

    setupStaticInjectors();

  }

  void initFromParams() { 

    IdealGasSolver::initFromParams();

    // set time0
    time0 = time;
    IF_RANK0 cout << " > time0: " << time0 << endl;
  }

  virtual void finalHook() {
    // write stats files when we exit
    if (mpi_rank == 0)
      cout << " > write lsp stats/pdf" << endl;

    for (list<LspStats*>::iterator stats = lspStatsList.begin(); stats != lspStatsList.end(); ++stats) {
      (*stats)->write(time,step,true);
    }

    for (list<LspPDF>::iterator pdfs = lspPDFList.begin(); pdfs != lspPDFList.end(); ++pdfs) {
      pdfs->write(time,step,true);
    }

    if (mpi_rank == 0)
      cout << " > lsp query_bc at final step..." << endl;
 
    queryLpBcs();

  }

protected:

  void setupInjectors() {

    IF_RANK0 cout << "setupInjectors()" << endl;
    // ===================
    // set the injectors
    // ===================
    FOR_PARAM_MATCHING("LSP.INJECTOR") {

      // check that the name is unique. It is used to register the ResidualMass, and
      // could be use in the future to identify which particle came from which
      // injector...

      string name = param->getString();
      for (list<LspInjector*>::iterator iter = injectorList.begin(); iter != injectorList.end(); ++iter) {
        if (name == (*iter)->getName()) {
          IF_RANK0 cerr << "LSP.INJECTOR: " << name << " already exists. Please change the name." << endl;
          throw(0);
        }
      }

      injectorList.push_back(new LspInjector);

      LspInjector * injector = injectorList.back();
      injector->initFromParams(&(*param),material,this,u);

    }

  }

  void initLspBc() { 

    if (mpi_rank == 0) cout << "initLspBc()" << endl;
    FOR_IZONE(bfZoneVec) lpZoneBcVec.push_back(NULL); //XXX Is lpZoneBcVec the right name?

    // Check to insure all of the LSP.BC commands are assigned to actual zones (guard against typos)
    FOR_PARAM_MATCHING("LSP.BC") { 
      bool found_zone = false;
      const string zone_name = param->getString(0);
      FOR_IZONE(bfZoneVec) {
        if (bfZoneVec[izone].getName() == zone_name) {
          found_zone = true;
          break;
        }
      }
      if (!found_zone)
        CERR("LSP.BC assigned to non-existent zone " << zone_name << "!"
             "Check your input file for syntax errors/typos!");
    }

    // Assing lp bcs to all zones in the simulation
    FOR_IZONE(bfZoneVec) { 
      // Assign the lsp bcs explicitly identified in the input file...
      FOR_PARAM_MATCHING("LSP.BC") { 
        if (bfZoneVec[izone].getName() == param->getString(0)) {
          const string bc_type = param->getString(1);
          if      (bc_type == "BOUNCE")
            lpZoneBcVec[izone] = new lpPerfectBounceBc(&bfZoneVec[izone],this);
          else if (bc_type == "IRREGULAR_BOUNCE")
            lpZoneBcVec[izone] = new lpIrregularBounceBc(&bfZoneVec[izone],this);
          else if (bc_type == "INELASTIC_BOUNCE") {
            lpZoneBcVec[izone] = new lpInelasticBounceBc(&bfZoneVec[izone],this);
            // HACK TODO solver templating?
            lpInelasticBounceBc * inelastic_bounce_bc = dynamic_cast<lpInelasticBounceBc*>(lpZoneBcVec[izone]);
            inelastic_bounce_bc->solid = dust; 
          }
          else if (bc_type == "OUTLET")
            lpZoneBcVec[izone] = new lpOutletBc(&bfZoneVec[izone],this);
          else
            CERR("ERROR: Unrecognized LSP.BC type " << bc_type << " assigned to zone " << bfZoneVec[izone].getName());
          break;// The bc was assigned so we can stop looking through the LSP.BC strings
        }
      }
      // This zone didn't receive an LSP.BC command, so we make it BOUNCE by default
      if (lpZoneBcVec[izone] == NULL) {
        CWARN("> WARNING: No LSP.BC assigned to zone: " << bfZoneVec[izone].getName() << " in the input fie."
              << "\nAssigning BOUNCE to this zone by default\n");
        lpZoneBcVec[izone]      = new lpPerfectBounceBc(&bfZoneVec[izone],this);
      }
    }

    // Call initialHook() to populate the required data
    for (vector<lpBaseBc*>::iterator it = lpZoneBcVec.begin(); it != lpZoneBcVec.end(); ++it)
      (*it)->initialHook();

    // Get the LSP.QUERY INTERVALS for each bc while we're here
    FOR_PARAM_MATCHING("LSP.QUERY_BC") { 
      string zone_name  = param->getString(0);
      int interval      = -1; 
      int iarg          = 1;
      while (iarg < param->size()) { 
        string token = param->getString(iarg++);
        if (token == "INTERVAL") interval = param->getInt(iarg++);
        else CWARN(" > Didn't recognize the token: \"" << token << "\" in LSP.QUERY_BC.");
      }

      // If not interval was provided, default to check_interval
      if (interval < 0) interval = check_interval;

      // Add this lsp.query to the the list
      FOR_IZONE(bfZoneVec) { 
        if (bfZoneVec[izone].getName() == zone_name) {
          COUT1(" > LSP.QUERY_BC, zone_name: " << zone_name << ", interval: " << interval);
          queryBcVec.push_back(pair<lpBaseBc*,int>(lpZoneBcVec[izone],interval));
          break;
        }
        // Error handling
        if (izone == bfZoneVec.size()-1)
          CWARN(" > LSP.QUERY_BC, zone_name: " << zone_name << " does not correspond to any zone. Skipping!");
      }
    }

    // Pull the interval at which to reset the lsp bc stats
    lsp_stats_interval = getIntParam("LSP.FLUSH_INTERVAL", -1);
    COUT1("LSP.FLUSH_INTERVAL = " << lsp_stats_interval);

  }

  void registerBoundaryConditions() {
    IdealGasSolver::registerBoundaryConditions();
  }

  void initHookLspBasic() {
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      lpCls->lp[ip].flag = KEEP;
      updateDp(ip);
    }
  }

  virtual void resizeData() {
    StaticSolver::resizeData();
    lpCls->realloc(lpCls->np);
  }

  void setupStaticInjectors() {
    IF_RANK0 cout << "setupStaticInjectors()" << endl;

    vector<LspData> lspDataVec;
    FOR_PARAM_MATCHING("LSP.STATIC_INJECTOR") {
      ensureCvAdt();
      StaticLspInjector sInjector; 
      sInjector.initFromParams(&(*param),material,this,u);
      sInjector.addNewLsp(lspDataVec);
    } 

    const int np_old = lpCls->size();
    lpCls->resize((np_old + lspDataVec.size()));
    for (int ip = 0; ip < lspDataVec.size(); ++ip) {
      lpCls->lp[np_old+ip].icv           = lspDataVec[ip].icv;
      FOR_I3 lpCls->lp[np_old+ip].xp[i]  = lspDataVec[ip].xp[i];
      FOR_I3 lpCls->lp[np_old+ip].up[i]  = lspDataVec[ip].up[i];
      lpCls->lp[np_old+ip].mp            = lspDataVec[ip].mp;
      lpCls->lp[np_old+ip].Tp            = lspDataVec[ip].Tp;
      lpCls->lp[np_old+ip].tbu           = 0.0;  // breakup time
      lpCls->lp[np_old+ip].flag          = KEEP; // recycling var
      lpCls->lp[np_old+ip].npar          = 1.0;  // parcel info
      updateDp(np_old+ip);
      const double e = lspDataVec[ip].e;
      const double f = lspDataVec[ip].f;
      assert((e>0)&&(e<=1));
      assert((f>0)&&(f<=1));
      lpCls->lp[np_old+ip].I             = lpCls->lp[np_old+ip].dp * pow(e/f, 1./3.);
      lpCls->lp[np_old+ip].S             = lpCls->lp[np_old+ip].I * f;
      lpCls->lp[np_old+ip].e             = e;
      lpCls->lp[np_old+ip].f             = f;
    }

    lspDataVec.clear();

    //MPI_Pause("HEREHEREHERE");

  }

  virtual int advanceSolution() {

    //wrapping around rhs's
    IdealGasRhs * rhs_arr[3] = {rhs0,rhs1,rhs2};

    //========================
    // inject lsp particles
    //========================
    injectNewLsp(time,dt,1);
    updateLspPrimitiveData(); // for now, just updates Dp from the mass and rho
    CtiRegister::clearCurrentData();
    timer.split("lsp injection");

    //========================
    // breakup lsp particles
    //========================
    switch(BreakupModel) {
    case NO_BREAKUP:
      break;
    case SBM_BREAKUP:
      breakupLsp();
      break;
    default:
      assert(0);
    }
    timer.split("lsp breakup model");

    //========================
    // lsp agglemoration
    //========================
    // ======================================
    // Bound the number of particles in a cv:
    // checkCvParcel() needs setPaocv()
    // checkCvParcel() generates TRASH particles
    // recycleLsp() is needed after checkCvParcel()
    // setPaocv() is needed for rk3Step_lsp_gas
    // ======================================
    // ======================
    // recycle the particles
    // setPaocv() is needed for rk3Step_lsp_gas
    // all cpus should go through recycleLsp(), so don't use agg logic
    // ======================
    setPaocv();     // we need paocv for rk3Step_lsp_gas
    timer.split("lsp set paocv");
    if (lspAgglom) {
      //reportLspSize();
      const bool agg = rebinLsp(); // the ranks which do agglemorate will need to reset the paocv later
      if (agg) {
        recycleLsp();
        setPaocv();  // if agg then particles are recycled -> setPaocv() is needed
      }
      //reportLspSize();
    }

    setLsp0();           // only at the beginning of time step... should be after agglemoration routine...

    //calcSurfParticleFlux(); // this needs Paocv_i/v

    //=======================================
    // for now new these inside the while loop cause the lspVec size changes (but can be more efficient)
    LspRhs * rhslsp0 = new LspRhs [lpCls->size()];
    LspRhs * rhslsp1 = new LspRhs [lpCls->size()];
    LspRhs * rhslsp2 = new LspRhs [lpCls->size()];
    LspRhs * rhslsp_arr[3] = {rhslsp0, rhslsp1, rhslsp2};

    LspRhs * rhslsp_d = new LspRhs [lpCls->size()]; // each substep sets to zero first (so one variable is enough)


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
    timer.split("lsp allocate rhs vectors");

    if (time_int == "2WAY_COUPLED") {
      // rk3 update...
      FOR_BCZONE (*it)->calcRhs(time,1);
      IdealGasSolver::calcRhs(rhs0,rhs0_sc,time,1);
      timer.split("calc rhs 1");

      calcRhsLsp_2w(rhslsp0,rhslsp_d);
      timer.split("calc rhs lsp 1");
      //addLspSrc(rhs0);

      FOR_BCZONE (*it)->rk3Step(dt,1);
      rk3Step_lsp_gas(erk_wgt3[0], rhs_arr, rhslsp_arr, rhslsp_d, 1);
      //rk3Step(erk_wgt3[0],dt,1);
      timer.split("solve coupled gas/lsp 1");

      time += dt;
      IdealGasSolver::updateConservativeAndPrimitiveData();
      updateLspPrimitiveData();
      CtiRegister::clearCurrentData();
      timer.split("update 1");

      // stage 2...
      FOR_BCZONE (*it)->calcRhs(time,2);
      IdealGasSolver::calcRhs(rhs1,rhs1_sc,time,2);
      timer.split("calc rhs 2");

      calcRhsLsp_2w(rhslsp1,rhslsp_d);
      timer.split("calc rhs lsp 2");
      //addLspSrc(rhs1);

      FOR_BCZONE (*it)->rk3Step(dt,2);
      rk3Step_lsp_gas(erk_wgt3[1], rhs_arr, rhslsp_arr, rhslsp_d, 2);
      //rk3Step(erk_wgt3[1],dt,2);
      timer.split("solve coupled gas/lsp 2");

      time -= 0.5*dt;
      updateConservativeAndPrimitiveData();
      updateLspPrimitiveData();
      CtiRegister::clearCurrentData();
      timer.split("update 2");

      // stage 3...
      FOR_BCZONE (*it)->calcRhs(time,3);
      IdealGasSolver::calcRhs(rhs2,rhs2_sc,time,3);
      timer.split("calc rhs 3");

      calcRhsLsp_2w(rhslsp2,rhslsp_d);
      timer.split("calc rhs lsp 3");
      //addLspSrc(rhs2);

      FOR_BCZONE (*it)->rk3Step(dt,3);
      rk3Step_lsp_gas(erk_wgt3[2], rhs_arr, rhslsp_arr, rhslsp_d, 3);
      //rk3Step(erk_wgt3[2],dt,3);
      timer.split("solve coupled gas/lsp 3");

      time += 0.5*dt;
      updateConservativeAndPrimitiveData();
      updateLspPrimitiveData();
      CtiRegister::clearCurrentData();
      timer.split("update 3");
    }
    else if (time_int == "1WAY_COUPLED") { // XXX
      // rk3 update...
      FOR_BCZONE (*it)->calcRhs(time,1);
      IdealGasSolver::calcRhs(rhs0,rhs0_sc,time,1);
      timer.split("calc rhs 1");

      calcRhsLsp_1w(rhslsp0,rhslsp_d);
      timer.split("calc rhs lsp 1");
      //addLspSrc(rhs0);

      FOR_BCZONE (*it)->rk3Step(dt,1);
      rk3Step(erk_wgt3[0],dt,1);
      rk3StepLsp_im(erk_wgt3[0],rhslsp_arr,rhslsp_d,1);
      timer.split("rk step gas/lsp 1");

      time += dt;
      IdealGasSolver::updateConservativeAndPrimitiveData();
      updateLspPrimitiveData();
      CtiRegister::clearCurrentData();
      timer.split("update 1");

      // stage 2...
      FOR_BCZONE (*it)->calcRhs(time,2);
      IdealGasSolver::calcRhs(rhs1,rhs1_sc,time,2);
      timer.split("calc rhs 2");

      calcRhsLsp_1w(rhslsp1,rhslsp_d);
      timer.split("calc rhs lsp 2");
      //addLspSrc(rhs1);

      FOR_BCZONE (*it)->rk3Step(dt,2);
      rk3Step(erk_wgt3[1],dt,2);
      rk3StepLsp_im(erk_wgt3[1],rhslsp_arr,rhslsp_d,2);
      timer.split("rk step gas/lsp 2");

      time -= 0.5*dt;
      updateConservativeAndPrimitiveData();
      updateLspPrimitiveData();
      CtiRegister::clearCurrentData();
      timer.split("update 2");

      // stage 3...
      FOR_BCZONE (*it)->calcRhs(time,3);
      IdealGasSolver::calcRhs(rhs2,rhs2_sc,time,3);
      timer.split("calc rhs 3");

      calcRhsLsp_1w(rhslsp2,rhslsp_d);
      timer.split("calc rhs lsp 3");
      //addLspSrc(rhs2);

      FOR_BCZONE (*it)->rk3Step(dt,3);
      rk3Step(erk_wgt3[2],dt,3);
      rk3StepLsp_im(erk_wgt3[2],rhslsp_arr,rhslsp_d,3);
      timer.split("rk step gas/lsp 3");

      time += 0.5*dt;
      updateConservativeAndPrimitiveData();
      updateLspPrimitiveData();
      CtiRegister::clearCurrentData();
      timer.split("update 3");
    }
    // compute sgs model and ensurue that material properties are
    // consistent with the new state that's been advanced..

    calcSgsAndMaterialProperties();
    timer.split("calc sgs and properties");

    // dump relevant information to stdout..

    report();
    timer.split("reporting");

    // report evaporated mass before recycling, change of size of lspVec, and relocating
    if (LIQ_LP)
      reportLspEvap();

    // at this point, we have updated x, and we need to relocate the particles in their new cvs...

    relocateLsp<LpClass<LspState>, LspState>(lpCls);
    //relocateLsp_rayTrace<LspState>();
    timer.split("lsp relocation");

    // apply additional plane trashes and recycle
    planeTrashSet();
    timer.split("lsp trash");
    recycleLsp(); // after this Paocv_i/v is not current, if you need them you should call setPaocv()
    timer.split("lsp recycle");

    //========================
    // lsp collision
    //========================
    switch(colMod) {
    case NO_MODEL:
      break;
    case STOCHASTIC:
      // need paocv for collision
      setPaocv();
      collision_stochastic();
      break;
    case DETERMINISTIC:
      // need paocv for collision
      setPaocv();
      collision_deterministic();
      break;
    default:
      assert(0);
    }
    timer.split("lsp collision");

    computeParticleConcentration();
    timer.split("lsp compute concentration");

    if (step%check_interval == 0)
      computeNearWallConcentration();

    //================
    // statistics
    //================
    updateLspStats(lpCls->lp,lpCls->size(),time,dt,step,1);
    timer.split("lsp stats");

    // take a look
    reportLsp();
    timer.split("lsp reporting");

    //checkParticleVF();
    //writeSurfData();
    //checkSanity();

    //monitor_jet();
    //calcInvicsidVortexErr(fh);

    delete[] rhslsp0;
    delete[] rhslsp1;
    delete[] rhslsp2;

    delete[] rhslsp_d;

    return 0;

  }


public:

  //void resetLocation() {
  //  for (int ip = 0; ip < lpCls->size(); ++ip) {
  //    FOR_I3 lpCls->lp[ip].xp[i] = lpCls->lp[ip].xp0[i];
  //  }
  //}

  // ===========================================================================
  // the original finction (reflective BC) is applied in StaticSolverWithLsp.hpp
  // ===========================================================================

  //virtual void impactLpBc(LspState& lsp, const double* xi, const double* st_normal) {

  //  /*
  //  const double rhof = fuel->calcRho(lsp.Tp);
  //  const double sigma = fuel->calcSigma(lsp.Tp);
  //  // we need to add fuel viscosity to the data base in CtiLiquid
  //  const double muf = 2.97e-4;
  //  const double umag = MAG(lsp.up);
  //  const double dp_b = lsp.dp; // droplet diameter before impact
  //  const double K = pow( pow(rhof,3) * pow(dp_b,3) * pow(umag,5) / pow(sigma,2) / muf ,1./4.);

  //  cout << "K: " << K << endl;

  //  if (K <= 57.7) {
  //    // Complete deposition of the droplet

  //  } else {
  //    // Partial deposition and partial splashing
  //    const double dp_a = min(1.0, 8.72*exp(-0.0281*K)) * dp_b;
  //    const int np_a = int(1.676e-5*pow(K,2.539));
  //    const double mdep_ratio = max(0.0, 1 - np_a*pow(dp_a,3)/pow(dp_b,3));
  //    const double d3_b = pow(dp_b,3);
  //    const double d3_a = np_a*pow(dp_a,3);
  //    cout << "d3_b: " << d3_b << " , d3_a: " << d3_a << endl;
  //    cout << "dp_b: " << dp_b << " , dp_a: " << dp_a << " , np_a: " << np_a << " , mdep_ratio: " << mdep_ratio << endl;

  //    for (int ip = 0; ip < np_a; ++ip) {
  //      const double dr = dp_a/dp_b;
  //      const double nmag2 = DOT_PRODUCT(st_normal,st_normal);
  //      const double un_b = DOT_PRODUCT(lsp.up,st_normal); assert(un_b>=0);
  //      double ut_b_v[3];
  //      FOR_I3 ut_b_v[i] = lsp.up[i] - un_b*st_normal[i]/nmag2;
  //      const double ut_b = MAG(ut_b_v);

  //      // non random part
  //      const double Un_a = (1.337 - 1.318*dr + 2.339*dr*dr)*un_b;
  //      const double Ut_a = (-0.249 - 2.959*dr + 7.794*dr*dr)*ut_b;

  //      const double std_un = (0.299*dr + 1.428*dr*dr)*un_b;
  //      const double std_ut = (0.130 + 1.941*dr + 6.539*dr*dr)*ut_b;

  //      // random part
  //      //const double un_a =

  //      //const double ut_b
  //      //const double dr =
  //      //const double un_a = 1.337 - 1.318*x + 2.339*x**2
  //    }

  //  }
  //  MPI_Pause("How is K?");
  //  */

  //  const double dx[3] = DIFF(lsp.xp,xi);
  //  const double dp = DOT_PRODUCT(st_normal,dx); assert(dp > 0.0);
  //  const double nmag2 = DOT_PRODUCT(st_normal,st_normal);
  //  FOR_I3 lsp.xp[i] -= 2.0*dp*st_normal[i]/nmag2;
  //  // and update xp0 to the xi_closest, because the remaining part of the
  //  // trajectory may have additional intersections...
  //  FOR_I3 lsp.xp0[i] = xi[i];
  //  const double un = DOT_PRODUCT(st_normal,lsp.up);
  //  FOR_I3 lsp.up[i] -= 2.0*un*st_normal[i]/nmag2;
  //}

protected:

  void checkParticleVF() {
    double* theta_p = new double [ncv];
    FOR_ICV theta_p[icv] = 0.0;
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      theta_p[lpCls->lp[ip].icv] += lpCls->lp[ip].mp/material->calcRho(lpCls->lp[ip].Tp);
    }
    FOR_ICV theta_p[icv] /= vol_cv[icv];
    dumpRange(theta_p,ncv,"theta_p");
  }

  void updateLspPrimitiveData() {
    // should not set xp0 here. it should be set only in the begining of the
    // time step not after each substep, otherwise the collision detection fails if
    // we set it equal to xp at the end of the time step
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      updateDp(ip);
    }
  }

  void updateDp(const int &ip) {
    // does it make sense to calc rho and store is here? 
    lpCls->lp[ip].dp = pow(lpCls->lp[ip].mp/(M_PI/6.0*lpCls->lp[ip].npar*material->calcRho(lpCls->lp[ip].Tp)) , 1./3.);
  }

  // rk3 explicit time advancement
  void rk3StepLsp(const double * rk_wgt, LspRhs ** rhslsp_arr,const int rkstep) {
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      LspRhs rhs_agg;
      rhs_agg.zero();
      for (int irk = 0; irk < rkstep; ++irk) {
        double wgt = dt*rk_wgt[irk];
        rhs_agg.axpy(wgt,rhslsp_arr[irk][ip]);
      }
      lpCls->lp[ip].add(rhs_agg);
    }
  }

  // rk3 implicit time advancement
  void rk3StepLsp_im(const double * rk_wgt, LspRhs ** rhslsp_arr, LspRhs* rhs_d, const int& rkstep) {
    for (int ip=0; ip < lpCls->size(); ++ip) {
      // build the new state
      LspRhs rhs_agg;
      rhs_agg.zero();
      for (int irk=0; irk < rkstep; ++irk) {
        double wgt = dt*rk_wgt[irk];
        rhs_agg.axpy(wgt,rhslsp_arr[irk][ip]);
      }
      const double temp = dt*rk_wgt[rkstep-1];
      lpCls->lp[ip].add_im(rhs_agg,rhs_d[ip],temp);
      // correct the rhs of previous state
      FOR_I3 rhslsp_arr[rkstep-1][ip].xp[i] += rhs_d[ip].xp[i]*lpCls->lp[ip].xp[i];
      FOR_I3 rhslsp_arr[rkstep-1][ip].up[i] += rhs_d[ip].up[i]*lpCls->lp[ip].up[i];
      rhslsp_arr[rkstep-1][ip].Tp += rhs_d[ip].Tp*lpCls->lp[ip].Tp;
      rhslsp_arr[rkstep-1][ip].mp += rhs_d[ip].mp*lpCls->lp[ip].mp;
    }
  }
  
  void rk3Step_lsp_gas(const double* rk_wgt, IdealGasRhs ** rhs_arr, LspRhs ** rhslsp_arr, const LspRhs * rhslsp_d, const int rkstep) {

    // we already know the maximum number of parcels per cv, npaocv_max, so 
    // allocate memory once...
    
    // below we build the system
    //
    // [ A_d v ]{ up,ug } = {rhs}
    // [ q   D ]

    LspState * lp = lpCls->lp;
    const int lp_size = lpCls->size();
    
    double * A_d = new double[npaocv_max];
    double * v   = new double[npaocv_max];
    double * q   = new double[npaocv_max];
    double * rhs[3];
    double * X[3];
    for (int i = 0; i < 3; ++i) {
      rhs[i] = new double [npaocv_max+1];
      X[i] = new double [npaocv_max+1];
    }
    
    const double cp = R_gas*gamma/(gamma-1.0);
    IdealGasRhs rhs_agg;
    FOR_ICV {
      rhs_agg.zero();
      const double tmp = dt/vol_cv[icv];
      for (int irk = 0; irk < rkstep; ++irk) {
        const double wgt = tmp*rk_wgt[irk];
        rhs_agg.rho += wgt*rhs_arr[irk][icv].rho;
        FOR_I3 rhs_agg.rhou[i] += wgt*rhs_arr[irk][icv].rhou[i];
        rhs_agg.rhoE += wgt*rhs_arr[irk][icv].rhoE;
        //rhs_agg.axpy(wgt,rhs_arr[irk][icv]);
      }
      const double rho_old = rho[icv];
      rho[icv] += rhs_agg.rho;
      // SHOULD ADD LATER XXX
      //cv[icv].Y = (rho_old*cv[icv].Y + rhs_agg.rhoY) / cv[icv].rho;
      
      // now solve the gas and particle velocity together
      const int np = paocv_i[icv+1] - paocv_i[icv];
      if (np == 0) {
        // solve gas velocity/energy with regular rk3
        FOR_I3 u[icv][i] = (u[icv][i]*rho_old + rhs_agg.rhou[i])/rho[icv];
        rhoE[icv] += rhs_agg.rhoE;
      }
      else {
        assert((np > 0)&&(np <= npaocv_max));
        ////////////////////////////////////////////////
        // Solve the coupled particles/gas velocities
        ////////////////////////////////////////////////
        double D = 0.0;
        for (int ip = paocv_i[icv]; ip != paocv_i[icv+1]; ++ip) {
          const int ii = ip - paocv_i[icv];
          assert((ii >= 0)&&(ii < np));
          assert((ip >= 0)&&(ip < lp_size));
          assert(rk_wgt[rkstep-1] > 0.0);
          assert(lp[ip].k >= 0.0);
          const double ki = lp[ip].k * rk_wgt[rkstep-1];
          assert(lp[ip].kb >= 0.0);
          const double kbi = lp[ip].kb * rk_wgt[rkstep-1];
          A_d[ii] = ki + 1.0/dt + kbi;
          v[ii] = -ki;
          q[ii] = -ki*lp[ip].mp;
          FOR_I3 rhs[i][ii] = lp[ip].up[i]/dt;
          for (int irk = 0; irk < rkstep; ++irk) {
            FOR_I3 rhs[i][ii] += rhslsp_arr[irk][ip].up[i]*rk_wgt[irk];
          }
          D += -q[ii];
        }
        D += rho[icv]/tmp;
        FOR_I3 rhs[i][np] = rhs_agg.rhou[i]/tmp + rho_old*u[icv][i]/tmp;
        //cout << "A_d[0]-1.0/dt: " << A_d[0]-1.0/dt << " ,rhslsp_arr[irk=0][0].up[0]: " << rhslsp_arr[0][0].up[0] << " ,D: " << D << " ,rhs[0][1]: " << rhs[0][1] << endl;
        //cout << "v[0]: " << v[0] << " ,q[0]:" << q[0] << endl;
        solveLinearSysGE(X, A_d, v, q, D, rhs, np);
        //checkSolution(X, np);
        //MPI_Pause("check the solution");
        ///////////////////////////////////////
        // update gas and particle velocities
        ///////////////////////////////////////
        FOR_I3 u[icv][i] = X[i][np];
        for (int ip = paocv_i[icv]; ip != paocv_i[icv+1]; ++ip) {
          const int ii = ip - paocv_i[icv];
          assert(ip < lp_size);
          FOR_I3 lp[ip].up[i] = X[i][ii];
          
          // no rk_wgt multiply to these rhs corrections:
          // update the rhs of gas rhou for last iteration
          // NOTE: rk_wgt is not multiplied
          FOR_I3 rhs_arr[rkstep-1][icv].rhou[i] += lp[ip].k*lp[ip].mp*(lp[ip].up[i] - u[icv][i]);
          // update the rhs of particle up for last iteration
          FOR_I3 rhslsp_arr[rkstep-1][ip].up[i] += lp[ip].k*(u[icv][i] - lp[ip].up[i]);
          // add the boundary piece: kb * (0-up) , this is the implicit part,
          FOR_I3 rhslsp_arr[rkstep-1][ip].up[i] += lp[ip].kb * (-lp[ip].up[i]);
        }
        ////////////////////////////////////////////////
        // Solve the coupled particle/gas energies
        ////////////////////////////////////////////////
        D = 0.0;
        rhs[0][np] = 0.0;
        const double ug2mag = DOT_PRODUCT(u[icv], u[icv]);
        for (int ip = paocv_i[icv]; ip != paocv_i[icv+1]; ++ip) {
          const int ii = ip - paocv_i[icv];
          assert(ii >= 0);  assert(ii < np);
          assert(ip < lp_size);
          assert(rk_wgt[rkstep-1] > 0.0);
          assert(lp[ip].kt >= 0.0);
          const double ki = lp[ip].kt * rk_wgt[rkstep-1];
          assert(lp[ip].ktb >= 0.0);
          const double kbi = lp[ip].ktb * rk_wgt[rkstep-1];
          A_d[ii] = ki + 1.0/dt + kbi;
          v[ii] = -ki/(rho[icv]*cp/gamma);
          //q[ii] = -ki*(fuel->calcRho(lspVec[ip].Tp) * fuel->calcCp(lspVec[ip].Tp)) * cv[icv].vol * lspVec[ip].npar;
          q[ii] = -ki*(material->calcRho(lp[ip].Tp) * material->calcCp(lp[ip].Tp)) * (M_PI/6.*pow(lp[ip].dp,3)) * lp[ip].npar;
          rhs[0][ii] = lp[ip].Tp/dt - ki * (0.5 * ug2mag) / (cp/gamma);
          rhs[0][np] += -q[ii] * (0.5 * ug2mag) / (cp/gamma);
          for (int irk = 0; irk < rkstep; ++irk) {
            rhs[0][ii] += rhslsp_arr[irk][ip].Tp*rk_wgt[irk];
          }
          D += -q[ii] / (rho[icv] * cp/gamma);
        }
        D += 1.0/tmp;
        rhs[0][np] += rhs_agg.rhoE/tmp + rhoE[icv]/tmp;
        //cout << "A_d[0]-1.0/dt: " << A_d[0]-1.0/dt << " ,rhslsp_arr[irk=0][0].Tp: " << rhslsp_arr[0][0].Tp << " ,D: " << D << " ,rhs[0][1]: " << rhs[0][1] << endl;
        //cout << "v[0]: " << v[0] << " ,q[0]:" << q[0] << endl;
        solveLinearSysGE(X[0], A_d, v, q, D, rhs[0], np);
        //checkSolution(X[0], np);
        ///////////////////////////////////////
        // update gas and particle temp/energy
        ///////////////////////////////////////
        rhoE[icv] = X[0][np];
        const double Tg = (rhoE[icv] - 0.5*rho[icv]*ug2mag)/(rho[icv]*cp/gamma);
        for (int ip = paocv_i[icv]; ip != paocv_i[icv+1]; ++ip) {
          const int ii = ip - paocv_i[icv];
          assert(ip < lp_size);
          lp[ip].Tp = X[0][ii];
          // update the rhs of gas rhoE for last iteration, NOTE: no rk_wgt should be added here
          rhs_arr[rkstep-1][icv].rhoE +=  (-q[ii]/rk_wgt[rkstep-1]) * (lp[ip].Tp - Tg); // npar is added in q above, devide by rk_wgt to elliminate it
          // update the rhs of particle Tp for last iteration (no need to mult npar cause the equation is normalized with mass)
          rhslsp_arr[rkstep-1][ip].Tp += lp[ip].kt * ( Tg - lp[ip].Tp);
          // add the boundary piece of the rk_wgt * ktb * (-Tp) , this is the implicit part, note that the Ts part is added already in calcRhsLsp
          rhslsp_arr[rkstep-1][ip].Tp += lp[ip].ktb * (-lp[ip].Tp);
        }
      }
    }

    // cleanup matrix memory...
    delete[] A_d;
    delete[] v;
    delete[] q;
    for (int i = 0; i < 3; ++i) {
      delete[] rhs[i];
      delete[] X[i];
    }
    
    //cout << "position: " << i_pos++ << endl;
    // update other particle properties, xp, mp...
    LspRhs lsp_rhs_agg;
    for (int ip = 0; ip < lp_size; ++ip) {
      lsp_rhs_agg.zero();
      for (int irk = 0; irk < rkstep; ++irk) {
        const double wgt = dt*rk_wgt[irk];
        lsp_rhs_agg.axpy(wgt,rhslsp_arr[irk][ip]);
      }
      const double tmp = dt*rk_wgt[rkstep-1];
      // following was lp[ip].add_im_xm(lsp_rhs_agg,rhslsp_d[ip],tmp);
      FOR_I3 lp[ip].xp[i] += lsp_rhs_agg.xp[i];
      lp[ip].mp = (lp[ip].mp + lsp_rhs_agg.mp)/(1.0 - rhslsp_d[ip].mp*tmp);
      // make sure mass is possitive
      lp[ip].mp = max(0.0,lp[ip].mp);
      if (lp[ip].mp == 0.0) lp[ip].flag = TRASH;
      FOR_I3 rhslsp_arr[rkstep-1][ip].xp[i] += rhslsp_d[ip].xp[i]*lp[ip].xp[i];
      rhslsp_arr[rkstep-1][ip].mp += rhslsp_d[ip].mp*lp[ip].mp;
    }
    
  }

  // These solvers implement Gauss Elimination that leverages the
  // diagonal structure of the matrix to minimize flop's...
  inline void solveLinearSysGE(double* X[3],const double* const A_d, const double* const v, const double* const q, const double D, double* const rhs[3], const int np) const {
    double Di = D;
    for (int ip = 0; ip<np; ++ip) {
      const double mult = -q[ip]/A_d[ip];
      Di += mult * v[ip];
      FOR_I3 rhs[i][np] += mult * rhs[i][ip];
    }
    FOR_I3 X[i][np] = rhs[i][np]/Di;
    for (int ip = 0; ip<np; ++ip) {
      FOR_I3 X[i][ip] = (rhs[i][ip] - v[ip]*X[i][np]) / A_d[ip];
    }
  }
  
  inline void solveLinearSysGE(double* X, const double* const A_d, const double * const v, const double* const q, const double D,double* const rhs, const int np) const {
    double Di = D;
    for (int ip = 0; ip<np; ++ip) {
      const double mult = -q[ip]/A_d[ip];
      Di += mult * v[ip];
      rhs[np] += mult * rhs[ip];
    }
    X[np] = rhs[np]/Di;
    for (int ip = 0; ip<np; ++ip) {
      X[ip] = (rhs[ip] - v[ip]*X[np]) / A_d[ip];
    }
  }

  void checkSolution(double** X, const int np) {
    double mymin_x = 1e20;
    double mymax_x = -1e20;
    for (int im = 0; im < np+1 ; ++im) {
      FOR_I3 {
        assert (X[i][im] == X[i][im]);
        if (X[i][im] < mymin_x)
          mymin_x = X[i][im];
        if (X[i][im] > mymax_x)
          mymax_x = X[i][im];
      }
    }
    //double min_x, max_x;
    //MPI_Reduce(&mymin_x,&min_x,1,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
    //MPI_Reduce(&mymax_x,&max_x,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    IF_RANK0
      cout << "Implicit solution range: min:max " << mymin_x << ":" << mymax_x << endl;
  }

  void checkSolution(double* X, const int np) {
    double mymin_x = 1e20;
    double mymax_x = -1e20;
    for (int im = 0; im < np ; ++im) {
      assert (X[im] == X[im]);
      if (X[im] < mymin_x)
        mymin_x = X[im];
      if (X[im] > mymax_x)
        mymax_x = X[im];
    }
    //double min_x, max_x;
    //MPI_Reduce(&mymin_x,&min_x,1,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
    //MPI_Reduce(&mymax_x,&max_x,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    IF_RANK0 {
      cout << "Implicit solution range: min:max " << mymin_x << ":" << mymax_x << endl;
      cout << "rhoE:" << X[np] << endl;
    }
  }

  void setLsp0() {
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      FOR_I3 lpCls->lp[ip].xp0[i] = lpCls->lp[ip].xp[i];
      lpCls->lp[ip].mp0 = lpCls->lp[ip].mp;
    }
  }

  // a wraper for calcRhsLsp, since the SOL_LP, and LIQ_LP are constant throughout the simulation
  // it would not add a significant overhead, however it is not optimal
  void calcRhsLsp_2w(LspRhs* rhs, LspRhs* rhs_d) {
    if (SOL_LP)
      calcRhsLsp_2w_sol(rhs, rhs_d);
    else if (LIQ_LP)
      calcRhsLsp_2w_liq(rhs, rhs_d);
    else
      CERR("LP type is not specified");
  }

  void calcRhsLsp_1w(LspRhs* rhs, LspRhs* rhs_d) {
    if (SOL_LP)
      calcRhsLsp_1w_sol(rhs, rhs_d);
    else if (LIQ_LP)
      calcRhsLsp_1w_liq(rhs, rhs_d);
    else
      CERR("LP type is not specified");
  }

  void calcRhsLsp_1w_sol(LspRhs* rhs, LspRhs* rhs_d) {

    //FOR_ICV lspSrc[icv].zero(); // only used for adding evaporated mass to Y
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      // initialize rhs with zero
      rhs[ip].zero();
      rhs_d[ip].zero();

      //#############################################
      // npar: we compute the Rhs of one droplet,
      //       the effect of npar is taken into acount
      //       in energy transfer and mass transfer to
      //       the gas only
      //#############################################
      // make sure mass is not zero otherwise:
      // rhs and rhs_d would remain zero for ip -> no change in xp, mp
      // set k and kt to zero -> no change in up, Tp (coupled to gas)
      // continue
      //#############################################
      if ((lpCls->lp[ip].flag == TRASH) || (lpCls->lp[ip].mp <= 0.0)) {
        lpCls->lp[ip].flag = TRASH; // make sure it will be trashed
        continue;
      }

      // gas properties
      const int icv     = lpCls->lp[ip].icv;
      const double rhog = rho[icv];
      const double mug  = mu_lam[icv];
      const double cpg  = R_gas*gamma/(gamma-1.0);
      const double kg   = loc_lam[icv] * cpg; // loc: lambda(k) over cp
      const double ug[3]= {u[icv][0], u[icv][1], u[icv][2]};
      const double Tg   = T[icv];

      // compute liquid properties at parcel temperature
      const double rhop = dust->calcRho(lpCls->lp[ip].Tp);
      const double Dp   = lpCls->lp[ip].dp;
      assert(Dp > 0.0);

      double du[3];
      FOR_I3 du[i] = ug[i] - lpCls->lp[ip].up[i];
      double dus = MAG(du);
      // compute non-dimensional and model parameters
      const double tau = rhop*Dp*Dp/(18.*mug);
      const double Rep = rhog*Dp*dus/mug;
      const double Prg = mug*cpg/kg;
      const double Nu  = 2. + 0.552*sqrt(Rep)*pow(Prg,1./3.);

      // drag law from  Crowe et al, Multiphase Flows w/ Droplets & Particles, 1998
      //const double f1 = 1.0+0.15*pow(Rep,0.687) + 0.0175*Rep/(1.0+4.25e4*pow(Rep,-1.16));
      // drag law from White, eqn 3-225
      // CD = 24.0/Rep + 6.0/(1.0+sqrt(Rep)) + 0.4 and f1 = Rep/24.0*CD but for numerical reasons when Rep=0 we add a correction to f1
      const double CD_c = 6.0/(1.0+sqrt(Rep)) + 0.4;
      double f1 =  1.0 + Rep/24.0*CD_c;
      //const double f1 = 1.0; // only for verificatiion step where exact solution needs to be compared

      const double f2 = 1.0;

      if (EllipsoidDrag_b) {
        const double e = lpCls->lp[ip].e;
        const double f = lpCls->lp[ip].f;
        const double Fs = f*pow(e,1.3);
        const double Fs3rd = pow(Fs,1./3.);
        const double Ks = 0.5*(Fs3rd + 1./Fs3rd);
        const double Fn = f*f*e;
        const double alpha2 = 0.45 + 10./(exp(2.5*log10(rhop/rhog))+30.);
        const double beta2  = 1. - 37./(exp(3.*log10(rhop/rhog))+100.);
        const double Kn = pow(10.,alpha2*pow(fabs(-log10(Fn)),beta2));
        f1 = Ks*(1. + 0.125*pow(Rep*Kn/Ks,2./3.)) + 0.46/24.*Kn*Rep*Rep/(Rep+5330.*Ks/Kn);
      }

      ///////////////////////////////
      // setting rhs and rhs_d
      // rhs_total = rhs + rhs_d*var
      ///////////////////////////////
      // set xp rhs
      FOR_I3 rhs[ip].xp[i] = lpCls->lp[ip].up[i];
      FOR_I3 rhs_d[ip].xp[i] = 0.0;

      // set up rhs
      FOR_I3 rhs[ip].up[i] = f1/tau*(ug[i]) + grav[i];
      FOR_I3 rhs_d[ip].up[i] = -f1/tau;

      // energy equation
      //rhs[ip].Tp = f2*Nu*(cpg/dust->calcCp(lpCls->lp[ip].Tp))/(3.0*Prg*tau) * (Tg - lpCls->lp[ip].Tp);
      double const coef = f2*Nu*(cpg/dust->calcCp(lpCls->lp[ip].Tp))/(3.0*Prg*tau);
      rhs[ip].Tp = coef * (Tg);
      rhs_d[ip].Tp = -coef;

      // mass transfer equation
      // mass does not change
      rhs[ip].mp = 0.0;
      rhs_d[ip].mp = 0.0;

      // update lspSrc to ideal gas
      //FOR_I3 lspSrc[icv].rhou[i] += -(lpCls->lp[ip].mp * rhs[ip].up[i]);

      // monitor if needed
      //if (step%check_interval==0) {
      //  cout << " > ip: "    << ip   << endl <<
      //          " > mug = "  << mug  << endl <<
      //          " > rhog = " << rhog << endl <<
      //          " > ug = " << COUT_VEC(ug) << endl;

      //  cout << " > rhop = " << rhop << endl <<
      //          " > Dp = " << Dp << endl;

      //  cout << " > tau = " << tau << endl <<
      //          " > Rep = " << Rep << endl;

      //  cout << " > f1 = " << f1 << endl;

      //}

    }
  }

  void calcRhsLsp_1w_liq(LspRhs* rhs, LspRhs* rhs_d) {

    //FOR_ICV lspSrc[icv].zero(); // only used for adding evaporated mass to Y
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      // initialize rhs with zero
      rhs[ip].zero();
      rhs_d[ip].zero();

      //#############################################
      // npar: we compute the Rhs of one droplet,
      //       the effect of npar is taken into acount
      //       in energy transfer and mass transfer to
      //       the gas only
      //#############################################
      // make sure mass is not zero otherwise:
      // rhs and rhs_d would remain zero for ip -> no change in xp, mp
      // set k and kt to zero -> no change in up, Tp (coupled to gas)
      // continue
      //#############################################
      if ((lpCls->lp[ip].flag == TRASH) || (lpCls->lp[ip].mp <= 0.0)) {
        lpCls->lp[ip].flag = TRASH; // make sure it will be trashed
        continue;
      }

      // gas properties
      const int icv     = lpCls->lp[ip].icv;
      const double Tr   = fuel->Tboil;;
      const double Tg   = T[icv];
      const double Pg   = p[icv];
      const double rhog = rho[icv];
      const double mug  = mu_lam[icv];
      const double cpg  = R_gas*gamma/(gamma-1.0);
      const double kg   = loc_lam[icv] * cpg; // loc: lambda(k) over cp
      //const double Yfg  = max(0.0,cv[icv].Y);
      const double Yfg  = 0.0;
      const double ug[3]= {u[icv][0], u[icv][1], u[icv][2]};
      const double MWg  = 28.97; // kg/m^3 for air

      // compute liquid properties at parcel temperature
      const double rhof = fuel->calcRho(lpCls->lp[ip].Tp);  // liquid density
      const double cpf  = fuel->calcCp(lpCls->lp[ip].Tp);   // liquid heat capacity
      const double Lv   = fuel->calcHvap(lpCls->lp[ip].Tp); // latent heat of evaporation
      const double Dp   = lpCls->lp[ip].dp;
      assert(Dp > 0.0);

      // equilibrium fuel composition at droplet surface
      const double nearOne = 0.9999999999;
      const double Xfseq   = min(nearOne,fuel->calcPv(lpCls->lp[ip].Tp)/Pg);
      const double Yfseq   = Xfseq/(Xfseq+(1.0-Xfseq)*MWg/fuel->MW);
      const double BM      = (Yfseq - Yfg)/(1.0 - Yfseq);            // Spalding mass transfer number

      // compute reference prop
      //const double Yfr = (2.*Yfseq + Yfg)/3.;  // 1/3 rule
      const double Yfr = Yfg;
      const double mur = Yfr*fuel->calcMuv(Tr) + (1.-Yfr)*mug;
      const double rhor= Yfr*fuel->calcRhov(Tr) + (1.-Yfr)*rhog;
      const double kr  = Yfr*fuel->calcKv(Tr) + (1.-Yfr)*kg;
      const double cpr = Yfr*fuel->calcCpv(Tr) + (1.-Yfr)*cpg;

      double du[3];
      FOR_I3 du[i] = ug[i] - lpCls->lp[ip].up[i];
      double dus = MAG(du);
      // compute non-dimensional and model parameters
      const double tau = rhof*Dp*Dp/(18.*mur);
      const double Rep = rhor*Dp*dus/mur;
      const double Prg = mur*cpr/kr;
      const double Scg = mur/(rhor*fuel->calcDv(Tr));
      const double Nu  = 2. + 0.552*sqrt(Rep)*pow(Prg,1./3.);
      const double Sh  = 2. + 0.552*sqrt(Rep)*pow(Scg,1./3.);
      const double Lk  = mur*sqrt(2.0*M_PI*lpCls->lp[ip].Tp*fuel->R_UNIVERSAL/fuel->MW)/(Scg*Pg);

      // Newton-iteration for evaporation rate mpdot
      const double F1 = (Sh/Scg/3.0)*(lpCls->lp[ip].mp/tau); // mdot + F1*ln(1+BM) = 0, see eq (4)
      const double F2 = -1.5*Prg*tau/lpCls->lp[ip].mp;        // beta = F2*mdot, see eq (17)
      double beta, Xfsneq, Yfsneq, BMneq, Fm;
      double mdot  = -1.0E-08;
      double mdot0 =  0.0;
      double Fm0    =  mdot0 + F1*log(1.0+BM); // initialize with equilib value
      double eps    =  1.0E-15;
      double resid  =  eps+1.0;
      int    iter   =  0;
      while (resid > eps) {
        iter += 1;
        beta   = F2*mdot;
        Xfsneq = Xfseq - (Lk/Dp*2.0)*beta;
        Yfsneq = Xfsneq/(Xfsneq+(1.0-Xfsneq)*MWg/fuel->MW);
        BMneq  = max(-nearOne,(Yfsneq-Yfg)/(1.0-Yfsneq));
        Fm     = mdot + F1*log(1.0+BMneq);
        double tmp = mdot;
        if (fabs(Fm-Fm0) < eps ) break;
        mdot = mdot - Fm/(Fm-Fm0)*(mdot-mdot0);
        resid = sqrt(pow(Fm-Fm0,2));
        Fm0 = Fm;
        mdot0 = tmp;
        if (iter > 20) cout << "ERROR: LSP Newton iteration did not converge !!!" << endl;
      }
      //cout << " > iter = " << iter << endl;
      beta = F2*mdot;

      // Dont let possitive mdot
      mdot = min(mdot, 0.0);
      assert(mdot<=0);

      const double mdot_over_mp = mdot/lpCls->lp[ip].mp;
      // drag law from  Crowe et al, Multiphase Flows w/ Droplets & Particles, 1998
      //const double f1 = 1.0+0.15*pow(Rep,0.687) + 0.0175*Rep/(1.0+4.25e4*pow(Rep,-1.16));
      // drag law from White, eqn 3-225
      // CD = 24.0/Rep + 6.0/(1.0+sqrt(Rep)) + 0.4 and f1 = Rep/24.0*CD but for numerical reasons when Rep=0 we add a correction to f1
      const double CD_c = 6.0/(1.0+sqrt(Rep)) + 0.4;
      double f1 =  1.0 + Rep/24.0*CD_c;
      double f2 = beta/(exp(beta)-1.0);
      //double f2 = 1.0; // only for non-evaporating case
      if (beta < 1e-16)
        f2 = 1.0;

      ///////////////////////////////
      // setting rhs and rhs_d
      // rhs_total = rhs + rhs_d*var
      ///////////////////////////////
      // set the xp
      FOR_I3 rhs[ip].xp[i] = lpCls->lp[ip].up[i];
      FOR_I3 rhs_d[ip].xp[i] = 0.0;

      // set up rhs
      FOR_I3 rhs[ip].up[i] = f1/tau*(ug[i]);
      FOR_I3 rhs_d[ip].up[i] = -f1/tau;

      // set the Tp
      const double theta1 = cpr/cpf;
      //rhs[ip].Tp = f2*Nu*theta1/(3.*Prg*tau)*(Tg - lpCls->lp[ip].Tp) + Lv/cpf*mdot_over_mp;

      double const coef = f2*Nu*theta1/(3.0*Prg*tau);
      rhs[ip].Tp = coef * (Tg) + Lv/cpf*mdot_over_mp;
      rhs_d[ip].Tp = -coef;

      // mass transfer equation
      rhs[ip].mp = 0.0;
      rhs_d[ip].mp = mdot_over_mp;

      // Y is the mass fraction of vapor in the gas (vapor mass/ carrier mass)
      // the lspSrc.rhoY is the evaporated vapor mass per time (integrated in vol)
      //
      // YOU SHOULD ADD THIS WHEN YOU HAVE RHOY XXX
      //
      //lspSrc[icv].rhoY += -mdot*lpCls->lp[ip].npar;

      // monitor if needed:
      //if (step%check_interval=0 && monitor_bool) {
      //  cout << " > ip: "    << ip   << endl <<
      //          " > Tp: "    << lpCls->lp[ip].Tp << endl;;
      //  cout << " > Yfg = "  << Yfg  << endl <<
      //          " > mug = "  << mug  << endl <<
      //          " > rhog = " << rhog << endl <<
      //          " > kg = "   << kg   << endl <<
      //          " > Pg = "   << Pg   << endl <<
      //          " > Tg = "   << Tg   << endl <<
      //          " > cpg = "  << cpg  << endl;

      //  cout << " > rhof = " << rhof << endl <<
      //          " > cpf = " << cpf << endl <<
      //          " > Lv = " << Lv << endl <<
      //          " > Dp = " << Dp << endl;

      //  cout << " > Xfseq = " << Xfseq << endl <<
      //          " > Yfseq = " << Yfseq << endl <<
      //          " > BM = " << BM << endl <<
      //          " > Tr " << Tr << endl;

      //  cout << " > Yfr = " << Yfr << endl <<
      //          " > mur = " << mur << endl <<
      //          " > rhor = " << rhor << endl <<
      //          " > kr = " << kr << endl <<
      //          " > cpr = " << cpr << endl;

      //  cout << " > tau = " << tau << endl <<
      //          " > Rep = " << Rep << endl <<
      //          " > Prg = " << Prg << endl <<
      //          " > Scg = " << Scg << endl <<
      //          " > Nu = " << Nu << endl <<
      //          " > Sh = " << Sh << endl <<
      //          " > Lk = " << Lk << endl;

      //  cout << "mdot: " << mdot << endl;

      //  cout << " > f1 = " << f1 << endl <<
      //          " > f2 = " << f2 << endl <<
      //          " > beta = " << beta << endl;
      //}
    }
  }

  // for dust particles
  void calcRhsLsp_2w_sol(LspRhs* rhs, LspRhs* rhs_d) {

    //FOR_ICV lspSrc[icv].zero(); // only used for adding evaporated mass to Y

    FOR_ICV {

      // compute fixed cv-based props here...
      // gas properties
      const double rhog = rho[icv];
      const double mug  = mu_lam[icv];
      const double cpg  = R_gas*gamma/(gamma-1.0);
      const double kg   = loc_lam[icv]*cpg; // loc: lambda(k) over cp
      const double ug[3] = { u[icv][0], u[icv][1], u[icv][2] };
      
      for (int ip = paocv_i[icv]; ip != paocv_i[icv+1]; ++ip) {
        
        // check sort...
        assert(icv == lpCls->lp[ip].icv);
        
        // initialize rhs with zero
        rhs[ip].zero();
        rhs_d[ip].zero();
        
        //#############################################
        // npar: we compute the Rhs of one droplet,
        //       the effect of npar is taken into acount
        //       in energy transfer and mass transfer to
        //       the gas only
        //#############################################
        // make sure mass is not zero otherwise:
        // rhs and rhs_d would remain zero for ip -> no change in xp, mp
        // set k and kt to zero -> no change in up, Tp (coupled to gas)
        // continue
        //#############################################
        if ((lpCls->lp[ip].flag == TRASH) || (lpCls->lp[ip].mp <= 0.0)) {
          lpCls->lp[ip].k = 0.0;
          lpCls->lp[ip].kt = 0.0;
          lpCls->lp[ip].kb = 0.0;
          lpCls->lp[ip].ktb = 0.0;
          lpCls->lp[ip].flag = TRASH; // make sure it will be trashed
          continue;
        }
        
        // compute particle properties at parcel temperature
        const double rhop = dust->calcRho(lpCls->lp[ip].Tp);
        const double Dp   = lpCls->lp[ip].dp;
        assert(Dp > 0.0);
        
        double du[3];
        FOR_I3 du[i] = ug[i] - lpCls->lp[ip].up[i];
        double dus = MAG(du); // slip velocity
        // compute non-dimensional and model parameters
        const double tau = rhop*Dp*Dp/(18.*mug);
        const double Rep = rhog*Dp*dus/mug;
        const double Prg = mug*cpg/kg;
        const double Nu  = 2. + 0.552*sqrt(Rep)*pow(Prg,1./3.);
        
        // drag law from  Crowe et al, Multiphase Flows w/ Droplets & Particles, 1998
        //const double f1 = 1.0+0.15*pow(Rep,0.687) + 0.0175*Rep/(1.0+4.25e4*pow(Rep,-1.16));
        // drag law from White, eqn 3-225
        // CD = 24.0/Rep + 6.0/(1.0+sqrt(Rep)) + 0.4 and f1 = Rep/24.0*CD but for numerical reasons when Rep=0 we add a correction to f1
        const double CD_c = 6.0/(1.0+sqrt(Rep)) + 0.4;
        double f1 =  1.0 + Rep/24.0*CD_c;
        //const double f1 = 1.0; // only for verificatiion step where exact solution needs to be compared
        
        if (EllipsoidDrag_b) {
          const double e = lpCls->lp[ip].e;
          const double f = lpCls->lp[ip].f;
          const double Fs = f*pow(e,1.3);
          const double Fs3rd = pow(Fs,1./3.);
          const double Ks = 0.5*(Fs3rd + 1./Fs3rd);
          const double Fn = f*f*e;
          const double alpha2 = 0.45 + 10./(exp(2.5*log10(rhop/rhog))+30.);
          const double beta2  = 1. - 37./(exp(3.*log10(rhop/rhog))+100.);
          const double Kn = pow(10.,alpha2*pow(fabs(-log10(Fn)),beta2));
          f1 = Ks*(1. + 0.125*pow(Rep*Kn/Ks,2./3.)) + 0.46/24.*Kn*Rep*Rep/(Rep+5330.*Ks/Kn);
        }

        const double f2 = 1.0;
        
        ///////////////////////////////
        // setting rhs and rhs_d
        // rhs_total = rhs + rhs_d*var
        // NOTE: rhs_d is only used for variables like xp and mp which are not coupled to gas
        ///////////////////////////////
        // set the xp
        FOR_I3 rhs[ip].xp[i] = lpCls->lp[ip].up[i];
        FOR_I3 rhs_d[ip].xp[i] = 0.0;

        // update the variable k and kt
        lpCls->lp[ip].k =  f1/tau;
        lpCls->lp[ip].kt = f2*Nu*(cpg/dust->calcCp(lpCls->lp[ip].Tp))/(3.0*Prg*tau);
        lpCls->lp[ip].kb =  0.0;
        lpCls->lp[ip].ktb = 0.0;

        // up would be solved implicitly coupled with ug later
        FOR_I3 rhs[ip].up[i] = grav[i];
        FOR_I3 rhs_d[ip].up[i] = 0.0;
        
        // energy equation
        //#############################################
        // coupling with the gas is treated implicitly and will
        // be solved coupled to gas later
        //#############################################
        rhs[ip].Tp = 0.0;
        rhs_d[ip].Tp = 0.0;
        
        // mass transfer equation, mass is not changed
        rhs[ip].mp = 0.0;
        rhs_d[ip].mp = 0.0;
        
        // monitor if needed
        //if (step%check_interval==0) {
        //  cout << " > ip: "    << ip   << endl <<
        //          " > mug = "  << mug  << endl <<
        //          " > rhog = " << rhog << endl <<
        //          " > ug = " << COUT_VEC(ug) << endl;
        
        //  cout << " > rhop = " << rhop << endl <<
        //          " > Dp = " << Dp << endl;
        
        //  cout << " > tau = " << tau << endl <<
        //          " > Rep = " << Rep << endl;
        
        //  cout << " > f1 = " << f1 << endl;
        
        //}
        //if ((lpCls->lp[ip].k<0)||(lpCls->lp[ip].k!=lpCls->lp[ip].k)) {
        //  cout << " > ip: "    << ip   << endl <<
        //          " > L I S = "  << L  << " " << I << " " << S << endl <<
        //          " > Fs Fs3rd = " << Fs << " " << Fs3rd << endl <<
        //          " > Ks = " << Ks << endl << 
        //          " > alpha2 beta2 = " << alpha2 << " " << beta2 << endl << 
        //          " > Fn = " << Fn << endl << 
        //          " > Kn temp temp1 = " << Kn << " " << temp << " " << temp1 << endl;
        //} 
        assert(lpCls->lp[ip].k>=0);
        assert(lpCls->lp[ip].kt>=0);
      }
    }
  }

  void calcRhsLsp_2w_liq(LspRhs* rhs, LspRhs* rhs_d) {

    FOR_ICV lspSrc[icv].zero(); // only used for adding evaporated mass to Y
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      // initialize rhs with zero
      rhs[ip].zero();
      rhs_d[ip].zero();

      //#############################################
      // npar: we compute the Rhs of one droplet,
      //       the effect of npar is taken into acount
      //       in energy transfer and mass transfer to
      //       the gas only
      //#############################################
      // make sure mass is not zero otherwise:
      // rhs and rhs_d would remain zero for ip -> no change in xp, mp
      // set k and kt to zero -> no change in up, Tp (coupled to gas)
      // continue
      //#############################################
      if ((lpCls->lp[ip].flag == TRASH) || (lpCls->lp[ip].mp <= 0.0)) {
        lpCls->lp[ip].k = 0.0;
        lpCls->lp[ip].kt = 0.0;
        lpCls->lp[ip].kb = 0.0;
        lpCls->lp[ip].ktb = 0.0;
        lpCls->lp[ip].flag = TRASH; // make sure it will be trashed
        continue;
      }

      // gas properties
      const int icv     = lpCls->lp[ip].icv;
      const double Tr   = fuel->Tboil;;
      //const double Tg   = T[icv];
      const double Pg   = p[icv];
      const double rhog = rho[icv];
      const double mug  = mu_lam[icv];
      const double cpg  = R_gas*gamma/(gamma-1.0);
      const double kg   = loc_lam[icv] * cpg; // loc: lambda(k) over cp
      //const double Yfg  = max(0.0,cv[icv].Y);
      const double Yfg  = 0.0;
      const double ug[3]= {u[icv][0], u[icv][1], u[icv][2]};
      const double MWg  = 28.97; // kg/m^3 for air

      // compute liquid properties at parcel temperature
      const double rhof = fuel->calcRho(lpCls->lp[ip].Tp);  // liquid density
      const double cpf  = fuel->calcCp(lpCls->lp[ip].Tp);   // liquid heat capacity
      const double Lv   = fuel->calcHvap(lpCls->lp[ip].Tp); // latent heat of evaporation
      const double Dp   = lpCls->lp[ip].dp;
      assert(Dp > 0.0);

      // equilibrium fuel composition at droplet surface
      const double nearOne = 0.9999999999;
      const double Xfseq   = min(nearOne,fuel->calcPv(lpCls->lp[ip].Tp)/Pg);
      const double Yfseq   = Xfseq/(Xfseq+(1.0-Xfseq)*MWg/fuel->MW);
      const double BM      = (Yfseq - Yfg)/(1.0 - Yfseq);            // Spalding mass transfer number

      // compute reference prop
      //const double Yfr = (2.*Yfseq + Yfg)/3.;  // 1/3 rule
      const double Yfr = Yfg;
      const double mur = Yfr*fuel->calcMuv(Tr) + (1.-Yfr)*mug;
      const double rhor= Yfr*fuel->calcRhov(Tr) + (1.-Yfr)*rhog;
      const double kr  = Yfr*fuel->calcKv(Tr) + (1.-Yfr)*kg;
      const double cpr = Yfr*fuel->calcCpv(Tr) + (1.-Yfr)*cpg;

      double du[3];
      FOR_I3 du[i] = ug[i] - lpCls->lp[ip].up[i];
      double dus = MAG(du);
      // compute non-dimensional and model parameters
      const double tau = rhof*Dp*Dp/(18.*mur);
      const double Rep = rhor*Dp*dus/mur;
      const double Prg = mur*cpr/kr;
      const double Scg = mur/(rhor*fuel->calcDv(Tr));
      const double Nu  = 2. + 0.552*sqrt(Rep)*pow(Prg,1./3.);
      const double Sh  = 2. + 0.552*sqrt(Rep)*pow(Scg,1./3.);
      const double Lk  = mur*sqrt(2.0*M_PI*lpCls->lp[ip].Tp*fuel->R_UNIVERSAL/fuel->MW)/(Scg*Pg);

      // Newton-iteration for evaporation rate mpdot
      const double F1 = (Sh/Scg/3.0)*(lpCls->lp[ip].mp/tau); // mdot + F1*ln(1+BM) = 0, see eq (4)
      const double F2 = -1.5*Prg*tau/lpCls->lp[ip].mp;        // beta = F2*mdot, see eq (17)
      double beta, Xfsneq, Yfsneq, BMneq, Fm;
      double mdot  = -1.0E-08;
      double mdot0 =  0.0;
      double Fm0    =  mdot0 + F1*log(1.0+BM); // initialize with equilib value
      double eps    =  1.0E-15;
      double resid  =  eps+1.0;
      int    iter   =  0;
      while (resid > eps) {
        iter += 1;
        beta   = F2*mdot;
        Xfsneq = Xfseq - (Lk/Dp*2.0)*beta;
        Yfsneq = Xfsneq/(Xfsneq+(1.0-Xfsneq)*MWg/fuel->MW);
        BMneq  = max(-nearOne,(Yfsneq-Yfg)/(1.0-Yfsneq));
        Fm     = mdot + F1*log(1.0+BMneq);
        double tmp = mdot;
        if (fabs(Fm-Fm0) < eps ) break;
        mdot = mdot - Fm/(Fm-Fm0)*(mdot-mdot0);
        resid = sqrt(pow(Fm-Fm0,2));
        Fm0 = Fm;
        mdot0 = tmp;
        if (iter > 20) cout << "ERROR: LSP Newton iteration did not converge !!!" << endl;
      }
      //cout << " > iter = " << iter << endl;
      beta = F2*mdot;

      // Dont let possitive mdot
      mdot = min(mdot, 0.0);
      assert(mdot<=0);

      const double mdot_over_mp = mdot/lpCls->lp[ip].mp;
      // drag law from  Crowe et al, Multiphase Flows w/ Droplets & Particles, 1998
      //const double f1 = 1.0+0.15*pow(Rep,0.687) + 0.0175*Rep/(1.0+4.25e4*pow(Rep,-1.16));
      // drag law from White, eqn 3-225
      // CD = 24.0/Rep + 6.0/(1.0+sqrt(Rep)) + 0.4 and f1 = Rep/24.0*CD but for numerical reasons when Rep=0 we add a correction to f1
      const double CD_c = 6.0/(1.0+sqrt(Rep)) + 0.4;
      double f1 =  1.0 + Rep/24.0*CD_c;
      double f2 = beta/(exp(beta)-1.0);
      //double f2 = 1.0; // only for non-evaporating case
      if (beta < 1e-16)
        f2 = 1.0;

      ///////////////////////////////
      // setting rhs and rhs_d
      // rhs_total = rhs + rhs_d*var
      ///////////////////////////////
      // set the xp
      FOR_I3 rhs[ip].xp[i] = lpCls->lp[ip].up[i];
      FOR_I3 rhs_d[ip].xp[i] = 0.0;

      lpCls->lp[ip].k = f1/tau;
      lpCls->lp[ip].kt = f2*Nu*(cpr/cpf)/(3.0*Prg*tau);
      lpCls->lp[ip].kb =  0.0;
      lpCls->lp[ip].ktb = 0.0;

      assert(lpCls->lp[ip].k>0);
      assert(lpCls->lp[ip].kt>0);
      // up would be solved implicitly coupled with ug later
      FOR_I3 rhs[ip].up[i] = 0.0;
      FOR_I3 rhs_d[ip].up[i] = 0.0;

      // energy equation
      //const double theta1 = cpr/cpf;
      //rhs[ip].Tp = f2*Nu*theta1/(3.*Prg*tau)*(Tg - P[ip].Tp) + Lv/cpf*mdot/P[ip].mp;
      //rhs[ip].Tp = f2*Nu*theta1/(3.*Prg*tau)*(Tg) + Lv/cpf*mdot_over_mp;
      //rhs_d[ip].Tp = - f2*Nu*theta1/(3.*Prg*tau);

      //#############################################
      // only evaporation latent heat is treated explicitly
      // coupling with the gas is treated implicitly and will
      // be solved coupled to gas later
      //#############################################
      //rhs[ip].Tp = 0.0;
      rhs[ip].Tp = Lv/cpf*mdot_over_mp;
      rhs_d[ip].Tp = 0.0;

      //#############################################
      // if you like to test a zero evaporation case, make sure
      // in addition to zeroing mass and temperature source terms
      // you make f2=1. This is the limit for non-evaporative case
      // otherwise your kt coeff is not correct
      //#############################################

      // mass transfer equation
      rhs[ip].mp = 0.0;
      rhs_d[ip].mp = mdot_over_mp;
      //rhs[ip].mp = mdot;
      //rhs_d[ip].mp = 0.0;

      // Y is the mass fraction of vapor in the gas (vapor mass/ carrier mass)
      // the lspSrc.rhoY is the evaporated vapor mass per time (integrated in vol)
      //
      // YOU SHOULD ADD THIS WHEN YOU HAVE RHOY XXX
      //
      //lspSrc[icv].rhoY += -mdot*lpCls->lp[ip].npar;

      // monitor if needed:
      //if (step%check_interval=0 && monitor_bool) {
      //  cout << " > ip: "    << ip   << endl <<
      //          " > Tp: "    << lpCls->lp[ip].Tp << endl;;
      //  cout << " > Yfg = "  << Yfg  << endl <<
      //          " > mug = "  << mug  << endl <<
      //          " > rhog = " << rhog << endl <<
      //          " > kg = "   << kg   << endl <<
      //          " > Pg = "   << Pg   << endl <<
      //          " > Tg = "   << Tg   << endl <<
      //          " > cpg = "  << cpg  << endl;

      //  cout << " > rhof = " << rhof << endl <<
      //          " > cpf = " << cpf << endl <<
      //          " > Lv = " << Lv << endl <<
      //          " > Dp = " << Dp << endl;

      //  cout << " > Xfseq = " << Xfseq << endl <<
      //          " > Yfseq = " << Yfseq << endl <<
      //          " > BM = " << BM << endl <<
      //          " > Tr " << Tr << endl;

      //  cout << " > Yfr = " << Yfr << endl <<
      //          " > mur = " << mur << endl <<
      //          " > rhor = " << rhor << endl <<
      //          " > kr = " << kr << endl <<
      //          " > cpr = " << cpr << endl;

      //  cout << " > tau = " << tau << endl <<
      //          " > Rep = " << Rep << endl <<
      //          " > Prg = " << Prg << endl <<
      //          " > Scg = " << Scg << endl <<
      //          " > Nu = " << Nu << endl <<
      //          " > Sh = " << Sh << endl <<
      //          " > Lk = " << Lk << endl;

      //  cout << "mdot: " << mdot << endl;

      //  cout << " > f1 = " << f1 << endl <<
      //          " > f2 = " << f2 << endl <<
      //          " > beta = " << beta << endl;
      //}
    }
  }

  // only for adding the evaporated mass into the scalar
  void addLspSrc(IdealGasRhs * rhs) {
    /*
    FOR_ICV {
      assert(lspSrc[icv].rho == 0.0);
      assert(lspSrc[icv].rhoE == 0.0);
      assert(lspSrc[icv].rhoY == 0.0);
      rhs[icv].add(lspSrc[icv]);
    }
    */
  }

  void reportLspEvap() {
    // check the evaporated mass in this time step
    if (step%check_interval == 0) {
      double myMp_s = 0.0;
      double myMp0_s = 0.0;
      double Mp_s;
      double Mp0_s;
      for (int ip = 0; ip < lpCls->size(); ++ip) {
        myMp_s += lpCls->lp[ip].mp;
        myMp0_s += lpCls->lp[ip].mp0;
      }
      MPI_Reduce(&myMp_s,&Mp_s,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      MPI_Reduce(&myMp0_s,&Mp0_s,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      IF_RANK0 {
        assert(Mp0_s >= Mp_s);
        cout << "evaporated mass: " << Mp0_s - Mp_s << endl;
      }
    }
  }

  void reportLsp() {

    if (step%check_interval == 0) {

      if ( getLspSize() > 0) {

        cti_real* r1 = new cti_real[lpCls->size()];
        cti_real (*r2)[3] = new cti_real[lpCls->size()][3];
        for (int ip = 0; ip < lpCls->size(); ++ip) {
          FOR_I3 r2[ip][i] = lpCls->lp[ip].xp[i];
        }
        dumpRange(r2,lpCls->size(),"xp");
        for (int ip = 0; ip < lpCls->size(); ++ip) {
        FOR_I3 r2[ip][i] = lpCls->lp[ip].up[i];
        }
        dumpRange(r2,lpCls->size(),"up");
        for (int ip = 0; ip < lpCls->size(); ++ip) {
          r1[ip] = lpCls->lp[ip].Tp;
        }
        dumpRange(r1,lpCls->size(),"Tp");
        for (int ip = 0; ip < lpCls->size(); ++ip) {
          r1[ip] = lpCls->lp[ip].mp;
      }
        dumpRange(r1,lpCls->size(),"mp");
        for (int ip = 0; ip < lpCls->size(); ++ip) {
          r1[ip] = lpCls->lp[ip].dp;
        }
        dumpRange(r1,lpCls->size(),"dp");
        for (int ip = 0; ip < lpCls->size(); ++ip) {
          r1[ip] = lpCls->lp[ip].npar;
        }
        dumpRange(r1,lpCls->size(),"npar");
        for (int ip = 0; ip < lpCls->size(); ++ip) {
          r1[ip] = lpCls->lp[ip].tbu;
        }
        dumpRange(r1,lpCls->size(),"tbu");
        for (int ip = 0; ip < lpCls->size(); ++ip) {
          r1[ip] = lpCls->lp[ip].e;
        }
        dumpRange(r1,lpCls->size(),"e");
        for (int ip = 0; ip < lpCls->size(); ++ip) {
          r1[ip] = lpCls->lp[ip].f;
        }
        dumpRange(r1,lpCls->size(),"f");
        for (int ip = 0; ip < lpCls->size(); ++ip) {
          lpCls->lp[ip].I = lpCls->lp[ip].dp * pow(lpCls->lp[ip].e/lpCls->lp[ip].f,1./3.);
          r1[ip] = lpCls->lp[ip].I;
        }
        dumpRange(r1,lpCls->size(),"I");
        for (int ip = 0; ip < lpCls->size(); ++ip) {
          lpCls->lp[ip].S = lpCls->lp[ip].I * lpCls->lp[ip].f;
          r1[ip] = lpCls->lp[ip].S;
        }
        dumpRange(r1,lpCls->size(),"S");

        dumpRange(CPar,ncv,"Cpar");
        dumpRange(CPar2,ncv,"Cpar2");

        COUT1(" > time_stat: "<<time_stat);

        //for (int ip = 0; ip < lpCls->size(); ++ip) {
        //  double urel[3];
        //  FOR_I3 urel[i] = u[lpCls->lp[ip].icv][i] - lpCls->lp[ip].up[i];
        //  double urel_mag = MAG(urel);
        //  double Weber = (rho[lpCls->lp[ip].icv] * urel_mag * urel_mag * lpCls->lp[ip].dp / 2.)/fuel->calcSigma(lpCls->lp[ip].Tp);
        //  r1[ip] = Weber;
        //}
        //dumpRange(r1,lpCls->size(),"Weber");

        int LspSize = getLspSize();
        IF_RANK0 cout << " > LspSize: " << LspSize << endl;

        double myMp = 0.0;
        double Mp;
        for (int ip = 0; ip < lpCls->size(); ++ip) {
          myMp += lpCls->lp[ip].mp;
        }
        MPI_Reduce(&myMp,&Mp,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
        IF_RANK0 cout << " > Lsp Total Mass: " << Mp << endl;

        delete[] r1;
        delete[] r2;
      }
    }

    queryLpBcs_interval();
    reportVolMassFrac();

    if (show_histogram && (step%check_interval==0)) {

      int LspSize = getLspSize();

      if (LspSize <= 1) return;

      // send/receive 
      int * np_rank = NULL;
      if (mpi_rank == 0) {
        np_rank = new int [mpi_size];
        np_rank[0] = lpCls->size(); // for self
        for (int im = 1; im < mpi_size; im++) {
          // get np_rank[im] from other processors
          MPI_Recv(&np_rank[im],1,MPI_INT,im,100,mpi_comm,MPI_STATUS_IGNORE);
        }
      }
      else {
        // give np to processor 0
        int my_np = lpCls->size();
        MPI_Send(&my_np,1,MPI_INT,0,100,mpi_comm);
      }

      if (mpi_rank == 0) {
        int np_total = 0;
        for (int im = 0; im < mpi_size; im++) np_total += np_rank[im];
        assert(np_total == LspSize);
      }

      // send/receive 
      vector<double> dp_vec;
      if (mpi_rank == 0) {
        for (int ip = 0; ip < np_rank[0]; ip++) dp_vec.push_back(lpCls->lp[ip].dp); // for self
        for (int im = 1; im < mpi_size; im++) {
          double * dp_rec = new double [np_rank[im]];
          MPI_Recv(dp_rec,np_rank[im],MPI_DOUBLE,im,99,mpi_comm,MPI_STATUS_IGNORE);
          for (int ip = 0; ip < np_rank[im]; ip++) dp_vec.push_back(dp_rec[ip]);
          delete[] dp_rec;
        }
      }
      else {
        double * dp_send = new double [lpCls->size()];
        for (int ip = 0; ip < lpCls->size(); ip++) dp_send[ip] = lpCls->lp[ip].dp;
        MPI_Send(dp_send,lpCls->size(),MPI_DOUBLE,0,99,mpi_comm);
        delete[] dp_send;
      }

      if (mpi_rank == 0) {
        assert(dp_vec.size()==LspSize);
        // histogram of size of paritlces in the domain
        vector<double> x;
        vector<double> y_count;
        vector<double> y_mass;
        dumpHist0(&dp_vec[0],LspSize,hist_nbins,"particle diameter count distribution",x,y_count,hist_verbose_count);
        x.clear();
        double this_rho = material->calcRho(300);
        dumpHistMass(&dp_vec[0],LspSize,this_rho,hist_nbins,"particle diameter mass distribution",x,y_mass,hist_verbose_mass);
        if ((hist_file!="") && ((step%hist_interval==0) || (hist_interval==-1))) {
          char filename[128];
          buildIndexedFilename(filename,hist_file.c_str(),step,"dat");
          mkdir_for_file(filename); // allow the user to have specified a subdirectory...
          ofstream file;
          file.open(filename, ios::out);
          file << "# step, time, samples: " << step << " " << time << " " << LspSize << endl;
          file << "# 1: diameter, 2: count dist, 3: count cdf 4: mass dist 5: mass cdf" << endl;
          assert(x.size()==y_count.size());
          assert(x.size()==y_mass.size());
          double integral_count = 0.0;
          double integral_mass = 0.0;
          const double dx = x[1] - x[0];
          for (int ibin = 0; ibin < x.size(); ibin++) {
            integral_count += y_count[ibin]*dx;
            integral_mass += y_mass[ibin]*dx;
            file << x[ibin] << " " << y_count[ibin] << " " << integral_count << " " << y_mass[ibin] << " " << integral_mass << endl;
          }
          file.close();
        }
      }

      if (np_rank!=NULL) delete[] np_rank;

    }

  }

  int searchClosestIpX(double* xp_debug) {

    int ip_return = -1;
    int ip_closest = -1;
    double d2_closest = 1.0E+20;
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      const double d2 = DIST2(lpCls->lp[ip].xp,xp_debug);
      if ((ip_closest == -1)||(d2 < d2_closest)) {
        ip_closest = ip;
        d2_closest = d2;
      }
    }

    DoubleInt myDi,di;
    myDi.this_double = d2_closest;
    myDi.this_int = mpi_rank;

    MPI_Allreduce(&myDi,&di,1,MPI_DOUBLE_INT,MPI_MINLOC,mpi_comm);

    if (di.this_int == mpi_rank) {
      assert(ip_closest >= 0);
      ip_return = ip_closest;
      cout << "XXXXXXXXXXXXX closest point is on rank: " << mpi_rank << " dist: " << sqrt(d2_closest) << " xp: " << COUT_VEC(lpCls->lp[ip_closest].xp) << endl;
    }

    return ip_return;

  }

  void checkSanity() {
    double xp_target[3] = {0.0745, -0.005, 0.0029};
    double d_min = 1e20;
    //int ip_min;
    int mycount_abnormal = 0;
    double margin = 0.0012;
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      double* xp = lpCls->lp[ip].xp;
      double* xcv = x_cv[lpCls->lp[ip].icv];
      double dx[3] = DIFF(xp,xcv);
      double mag_dx = MAG(dx);
      if (mag_dx > margin) {
        cout << " > particle ip: " << ip << " rank: " << mpi_rank << " mag_dx: " << mag_dx << " xcv: " << COUT_VEC(xcv) << " xp: " << COUT_VEC(xp) << " OUT_OF_BOUNCE " << endl;
        mycount_abnormal++;
      }
        if (mpi_rank == 98) {
        double d[3] = DIFF(xp, xp_target);
        double d_mag = MAG(d);
        if (d_min > d_mag) {
          d_min = d_mag;
          //ip_min = ip;
        }
      }
    }
    if (mpi_rank == 98) {
      //cout << "closest particle: " << " ip: " << ip_min << " d_min: " << d_min << endl;
      //cout <<  " ip: " << 4 << " xp " << COUT_VEC(lpCls->lp[4].xp) << endl;
    }
    int count_abnormal;
    MPI_Reduce(&mycount_abnormal,&count_abnormal,1,MPI_INT,MPI_SUM,0,mpi_comm);
    MPI_Bcast(&count_abnormal,1,MPI_INT,0,mpi_comm);
    if (step%check_interval == 0)
      IF_RANK0 cout << " > Total of " << count_abnormal << " out_of_bounce particles" << endl;
    if (count_abnormal > 0) {
      IF_RANK0 cout << " > step: " << step << endl;
      //MPI_Pause("Got an out of bounce particle, check!");
    }
  }

  // set the flag of particles which cross the planeTrashList
  // you should call this at the end of the time step so that the xp0 is set
  void planeTrashSet() {
    if (plane_trash_bool) {
      for (int ip = 0; ip < lpCls->size(); ++ip) {
        bool passPlane = false;
        for (list<planeTrashClass*>::iterator it = planeTrashList.begin(); it != planeTrashList.end(); ++it) {
          passPlane = passPlane || (*it)->passPlane(lpCls->lp[ip].xp0, lpCls->lp[ip].xp);
        }
        if (passPlane)
          lpCls->lp[ip].flag = TRASH;
      }
    }
  }
  // this function takes out the particles with flag=TRASH out. You should set the
  // flag of the particles outside of this function preferencially.
  void recycleLsp(bool verbose=true) {
    // set the flag (option)
    int ntrash = lpCls->recycleLp();
    if ((step%check_interval == 0) && mpi_rank==0 && verbose)
      cout << " > recycleLsp(): # of trashed particles: " << ntrash << endl;
  }

  void setPaocv() {
    //IF_RANK0 cout << "setPaocv()" << endl;
  
    // count in paocv_i, offset by +1...
    assert(paocv_i);
    FOR_ICV paocv_i[icv+1] = 0;
    const int np = lpCls->size();
    for (int ip = 0; ip < np; ++ip) {
      const int icv = lpCls->lp[ip].icv;
      assert((icv >= 0)&&(icv < ncv));
      ++paocv_i[icv+1];
    }

    // turn paocv_i into a CSR structure...
    paocv_i[0] = 0;
    npaocv_max = 0;
    FOR_ICV {
      npaocv_max = max(npaocv_max,paocv_i[icv+1]);
      paocv_i[icv+1] += paocv_i[icv];
    }
    assert(paocv_i[ncv] == np);
    
    // completely reorder the lpCls->lp array...
    
    LspState * lp_new = new LspState[lpCls->size_max()]; 
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      const int icv = lpCls->lp[ip].icv;
      lp_new[paocv_i[icv]] = lpCls->lp[ip];
      ++paocv_i[icv];
    }
    delete[] lpCls->lp;
    lpCls->lp = lp_new;
    
    for (int icv = ncv-1; icv > 0; --icv)
      paocv_i[icv] = paocv_i[icv-1];
    paocv_i[0] = 0;
    
  }
  
  bool rebinLsp() {
    
    IF_RANK0 cout << "  >> rebinLsp()";
    //const int MAX_PPC = 10000;
    const int N_BINS = max(1,int(MAX_PPC/2));
    //const int MAX_NPAR = 1000;

    double binM[N_BINS];
    double binT[N_BINS];
    double binD[N_BINS];
    double binTbu[N_BINS];
    double binU[N_BINS][3];
    double binX[N_BINS][3];

    bool done_agg = false;
    FOR_ICV {
      const int npa = paocv_i[icv+1] - paocv_i[icv];
      if (npa >= MAX_PPC) {
        done_agg = true;
        cout << ".";
        cout.flush();
        for (int ibin = 0; ibin < N_BINS; ++ibin) {
          binM[ibin] = 0.0;
          binT[ibin] = 0.0;
          binD[ibin] = 0.0;
          binTbu[ibin] = 0.0;
          FOR_I3 binU[ibin][i] = 0.0;
          FOR_I3 binX[ibin][i] = 0.0;
        }
        // set the min/max of dp
        double d_min = 1e20;
        double d_max = -1e20;
        const double SMALL = 1e-15;
        for (int ip = paocv_i[icv]; ip != paocv_i[icv+1]; ++ip) {
          d_min = min(d_min, lpCls->lp[ip].dp);
          d_max = max(d_max, lpCls->lp[ip].dp);
        }
        assert(d_max>d_min);
        d_max += SMALL;
        d_min -= SMALL;
        const double dD = (d_max - d_min)/N_BINS;
        for (int ibin = 0; ibin < N_BINS; ++ibin) {
          binD[ibin] = d_min + (double(ibin)+0.5)*dD;
          //cout << "binD[" << ibin << "]: " << binD[ibin] << endl;
        }
        //cout << "d_min: " << d_min << " , d_max: " << d_max << endl;
        //cout << "dD: " << dD << endl;
        for (int ip = paocv_i[icv]; ip != paocv_i[icv+1]; ++ip) {
          const int ibin = (int)floor( (lpCls->lp[ip].dp - d_min)/dD );
          //cout << "ibin: " << ibin << " , dp: " << lspVec[ip].dp << " , binD: " << binD[ibin] << " , dp-binD: " << lspVec[ip].dp-binD[ibin] << endl;
          assert(ibin >= 0); assert(ibin < N_BINS);
          if (lpCls->lp[ip].npar < MAX_NPAR) {
            binM[ibin] += lpCls->lp[ip].mp;
            binT[ibin] += lpCls->lp[ip].Tp * lpCls->lp[ip].mp;
            binTbu[ibin] += lpCls->lp[ip].tbu * lpCls->lp[ip].mp;
            FOR_I3 binU[ibin][i] += lpCls->lp[ip].up[i] * lpCls->lp[ip].mp;
            FOR_I3 binX[ibin][i] += lpCls->lp[ip].xp[i] * lpCls->lp[ip].mp;
            // trash all the current particles, we will add new ones at the end later
            lpCls->lp[ip].flag = TRASH;
          }
        }
        double npar_min = 1e20;
        double npar_max = -1e20;
        for (int ibin = 0; ibin < N_BINS; ++ibin) {
          if (binM[ibin] == 0.0)
            continue;
          LspState lspTmp;
          lspTmp.dp = binD[ibin];
          lspTmp.mp = binM[ibin];
          lspTmp.Tp = binT[ibin]/binM[ibin];
          FOR_I3 lspTmp.up[i] = binU[ibin][i]/binM[ibin];
          FOR_I3 {
            lspTmp.xp[i] = binX[ibin][i]/binM[ibin];
            lspTmp.xp0[i] = lspTmp.xp[i];
          }
          lspTmp.npar = lspTmp.mp/(M_PI/6.0*material->calcRho(lspTmp.Tp)*pow(lspTmp.dp,3));
          lspTmp.icv =  icv;
          lspTmp.flag =  KEEP;
          // for now mass average tbu is assigned to the new parcel
          lspTmp.tbu = binTbu[ibin]/binM[ibin];

          //npar_sum += lspTmp.npar;
          npar_min = min(npar_min, lspTmp.npar);
          npar_max = max(npar_max, lspTmp.npar);
          //cout << "dp: " << lspTmp.dp << " , mp: " << lspTmp.mp << " , Tp: " << lspTmp.Tp << " , up: " << COUT_VEC(lspTmp.up) << " , npar: " << lspTmp.npar << " , tbu: " << lspTmp.tbu << endl;

          lpCls->addLp(lspTmp);
        }
        //cout << "npar_min: " << npar_min << " , npar_max: " << npar_max << endl;
      }
    }
    if (mpi_rank == 0) cout << endl;
    return done_agg;
  }

  void monitor_jet() {
    double my_x_max = 0.;
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      my_x_max = max(my_x_max, lpCls->lp[ip].xp[0]);
    }
    double x_max;
    MPI_Reduce(&my_x_max, &x_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
    ofstream fh;
    fh.open("jet_penetration.txt", ios::out | ios::app);
    IF_RANK0 fh << scientific << step << "   " << time << "   " << x_max << endl;
    fh.close();
  }

  // tracks the location of a single particle
  void monitor() {

    assert(getLspSize() == 1);
    if (lpCls->size() == 1 && step%check_interval==0) {
      //cout << "rank: " << rank << " , icv: " << icv << " , ncv: " << ncv << " , icv-ncv: " << lpCls->lp[0].icv-ncv << endl;
      const double tau = material->calcRho(lpCls->lp[0].Tp) * lpCls->lp[0].dp * lpCls->lp[0].dp / (18.0 * mu_lam[lpCls->lp[0].icv]);
      cout << " > rank: " << mpi_rank << " ,icv: " << lpCls->lp[0].icv << " ,time: " << time << " ,tau: " << tau << " ,xp: " << COUT_VEC(lpCls->lp[0].xp)
           << " ,up: " << COUT_VEC(lpCls->lp[0].up) << endl;
      //fprintf(fh,"%10.8e %10.8e %10.8e %10.8e\n",time, tau, lpCls->lp[0].xp[0], lpCls->lp[0].up[0]);
    //  fh << scientific << time << "   " << tau << "   " << lpCls->lp[0].xp[0] << "   " << lpCls->lp[0].up[0] << endl;
    }
    return;
  }

  void writeLspTecplot(const int step) {

    char filename[128];
    sprintf(filename,"lsp.%08d.dat",step);

    COUT1("writeLspTecplot: " << filename);

    FILE * fp;
    if ( mpi_rank == 0 ) {
      fp = fopen(filename,"w");
      assert(fp != NULL);
      fprintf(fp,"TITLE = \"flagged faces\"\n");
      fprintf(fp,"VARIABLES = \"X\"\n");
      fprintf(fp,"\"Y\"\n");
      fprintf(fp,"\"Z\"\n");
      fprintf(fp,"\"U-X\"\n");
      fprintf(fp,"\"U-Y\"\n");
      fprintf(fp,"\"U-Z\"\n");
      fprintf(fp,"\"D\"\n");
    }
    else {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      fp = fopen(filename,"a");
      assert(fp != NULL);
    }

    //cout << "HACK ******************** writing mpi_rank in delta ****************************** HACK" << endl;

    for (int ip = 0; ip < lpCls->size(); ++ip) {
      fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le\n",
        lpCls->lp[ip].xp[0],
        lpCls->lp[ip].xp[1],
        lpCls->lp[ip].xp[2],
        lpCls->lp[ip].up[0],
        lpCls->lp[ip].up[1],
        lpCls->lp[ip].up[2],
        lpCls->lp[ip].dp);
    }

    fclose(fp);

    if ( mpi_rank < mpi_size-1 ) {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }

    MPI_Barrier(mpi_comm);

  }

  void getNewLspFromInjectors(vector<LspData>& lspDataVec,
            const double time,const double dt,const bool verbose) {

    lspDataVec.clear();

    for (list<LspInjector*>::iterator iter = injectorList.begin(); iter != injectorList.end(); ++iter) {

      (*iter)->addNewLsp(lspDataVec,(time-time0),dt,verbose);

    }

  }

  void injectNewLsp(const double time,const double dt,const bool verbose) {
    vector<LspData> lspDataVec;
    getNewLspFromInjectors(lspDataVec,time,dt,verbose);

    // we have already established in the injector the fact that 
    // all the contents of lspDataVec should be injected in this rank

    int my_npnew = lspDataVec.size();
    int npnew_global;
    MPI_Reduce(&my_npnew, &npnew_global, 1, MPI_INT, MPI_SUM, 0, mpi_comm);
    np_inject += npnew_global;

    if ((mpi_rank == 0) && verbose && (step%check_interval == 0)) {
      cout << " > injecting total of " << npnew_global << " new parcels" << endl;
      cout << " > injected so far: " << np_inject << endl;
    }

    int npold = lpCls->size();
    lpCls->resize(npold+my_npnew);
    for (int ip = 0; ip < my_npnew; ++ip) {
      // local cv index
      lpCls->lp[npold+ip].icv = lspDataVec[ip].icv;
      // data...
      FOR_I3 lpCls->lp[npold+ip].xp[i]  = lspDataVec[ip].xp[i];
      FOR_I3 lpCls->lp[npold+ip].up[i]  = lspDataVec[ip].up[i];
      lpCls->lp[npold+ip].mp            = lspDataVec[ip].mp;
      lpCls->lp[npold+ip].Tp            = lspDataVec[ip].Tp;
      lpCls->lp[npold+ip].tbu           = 0.0;  // breakup time
      lpCls->lp[npold+ip].flag          = KEEP; // recycling var
      lpCls->lp[npold+ip].npar          = 1.0;  // parcel info
      updateDp(npold+ip);
      const double e = lspDataVec[ip].e;
      const double f = lspDataVec[ip].f;
      assert((e>0)&&(e<=1));
      assert((f>0)&&(f<=1));
      lpCls->lp[npold+ip].I             = lpCls->lp[npold+ip].dp * pow(e/f, 1./3.);
      lpCls->lp[npold+ip].S             = lpCls->lp[npold+ip].I * f;
      lpCls->lp[npold+ip].e             = e;
      lpCls->lp[npold+ip].f             = f;
    }


    //int npnew_global = lspDataVec.size();
    //np_inject += npnew_global;
    //if ((mpi_rank == 0) && verbose && (step%check_interval == 0)) {
    //  cout << " > injecting total of " << npnew_global << " new parcels" << endl;
    //  cout << " > injected so far: " << np_inject << endl;
    //}

    //if ( npnew_global > 0 ) {
    //  // figure out which particles are ours...
    //  // We are going to use global cv indexing to do this...

    //  int8 * cvora = NULL;
    //  buildXora(cvora,ncv);
    //  assert( cvora[mpi_rank+1] - cvora[mpi_rank] == ncv );

    //  vector<int> cvList;
    //  int * my_icvp_global = new int[npnew_global];

    //  for (int ip = 0; ip < npnew_global; ++ip) {
    //    my_icvp_global[ip] = -1;
    //    cvAdt->buildListForPoint(cvList, lspDataVec[ip].xp);

    //    const int icv = getClosestICVFromList(cvList,lspDataVec[ip].xp);
    //    if (icv >= 0) my_icvp_global[ip] = icv + cvora[mpi_rank]; // local to global cv offset

    //  }
    //  // at this point, my_icvp_global should contain the global cv index of
    //  // the containing cell, or -1... get the maximum across the
    //  // processors to break ties or check for particles that could not be located...

    //  int * icvp_global = new int[npnew_global];
    //  MPI_Allreduce(my_icvp_global,icvp_global,npnew_global,MPI_INT,MPI_MAX,mpi_comm);

    //  int npnew = 0;
    //  for (int ip = 0; ip < npnew_global; ++ip) {
    //    if (icvp_global[ip] == -1) {
    //      if (mpi_rank == 0)
    //        cout << "Warning: could not locate cv for particle at: " <<
    //          lspDataVec[ip].xp[0] << " " <<
    //          lspDataVec[ip].xp[1] << " " <<
    //          lspDataVec[ip].xp[2] << endl;
    //    }
    //    else if (icvp_global[ip] == my_icvp_global[ip]) {
    //      ++npnew;
    //    }
    //  }

    //  // resize and add particles on processors that need it...

    //  if (npnew > 0) {

    //    int npold = lpCls->size();
    //    lpCls->resize(npold+npnew);

    //    // loop through all the particles again
    //    int npnew_check = 0;
    //    for (int ip = 0; ip < npnew_global; ++ip) {
    //      // if the particles belongs to the processors
    //      if ((icvp_global[ip] != -1)&&(icvp_global[ip] == my_icvp_global[ip])) {
    //        int ipnew = npnew_check++;
    //        // local cv index
    //        lpCls->lp[npold+ipnew].icv = my_icvp_global[ip] - cvora[mpi_rank];
    //        // data...
    //        FOR_I3 lpCls->lp[npold+ipnew].xp[i]  = lspDataVec[ip].xp[i];
    //        FOR_I3 lpCls->lp[npold+ipnew].up[i]  = lspDataVec[ip].up[i];
    //        //FOR_I3 lpCls->lp[npold+ipnew].up[i]  = u[lpCls->lp[npold+ipnew].icv][i]; //XXX this is a hack to set the velocity as the gas velocity, shouldn't forget to change it back XXX TODO XXX TODO XXX
    //        lpCls->lp[npold+ipnew].mp            = lspDataVec[ip].mp;
    //        lpCls->lp[npold+ipnew].Tp            = lspDataVec[ip].Tp;
    //        // set the source terms to zero...
    //        lpCls->lp[npold+ipnew].tbu           = 0.0;  // breakup time
    //        lpCls->lp[npold+ipnew].flag          = KEEP; // recycling var
    //        lpCls->lp[npold+ipnew].npar          = 1.0;  // parcel info
    //        updateDp(npold+ipnew);
    //        //lsp->msg[npold+ipnew]           = 0.0;
    //        //FOR_I3 lsp->fsg[npold+ipnew][i] = 0.0;
    //      }
    //    }
    //    assert( npnew_check == npnew );
    //  }
    //  delete[] cvora;
    //  delete[] my_icvp_global;
    //  delete[] icvp_global;
    //}
  }

  int getClosestICVFromList(vector<int> &cvList, const double* xp) {
    double dist_closest = 1e20;
    int icv_closest = -1;
    for (int i = 0; i < cvList.size(); ++i) {
      const int icv = cvList[i];
      double dist = DIST(xp,x_cv[icv]);
      if (dist < dist_closest) {
        dist_closest = dist;
        icv_closest = icv;
      }
    }
    return(icv_closest);
  }

  // Fokker-Plank stocastic breakup model
  void breakupLsp() {

    assert(fuel != NULL);
    //IF_RANK0
    //  cout << " > breakupLsp()" << endl;
    const double SMALL = 1e-15;
    // save the number of parents
    int np_old = lpCls->size();
    //int np_old_tot = getLspSize();
    //IF_RANK0
    //  cout << " > np_old total: " << np_old_tot << endl;
    // go over all the particles
    int ip = 0;
    while (ip < lpCls->size()) {
      // copy parents properties
      const double mp_parent   = lpCls->lp[ip].mp; assert(mp_parent>0.0);
      const double Tp_parent   = lpCls->lp[ip].Tp;
      const double rhop_parent = fuel->calcRho(Tp_parent);
      const double npar_parent = lpCls->lp[ip].npar;
      //const double rp_parent   = 0.5*pow(mp_parent/(M_PI/6.0*rhop_parent),(1.0/3.0));
      const double rp_parent   = 0.5*pow(mp_parent/(M_PI/6.0*rhop_parent*npar_parent),(1.0/3.0));
      double xp_parent[3], up_parent[3];
      FOR_I3 {
        xp_parent[i] = lpCls->lp[ip].xp[i];
        up_parent[i] = lpCls->lp[ip].up[i];
      }

      const int icv = lpCls->lp[ip].icv;
      const double rhog = rho[icv]; // cv[lpCls->lp[ip].icv].rho;

      // setting velocity at the location of particle
      double dx[3];
      FOR_I3 dx[i] = lpCls->lp[ip].xp[i] - x_cv[icv][i];
      double Delta_u[3];
      FOR_I3 Delta_u[i] = dudx[icv][i][0]*dx[0] + dudx[icv][i][1]*dx[1] + dudx[icv][i][2]*dx[2];
      double ug[3];
      FOR_I3 ug[i] = u[icv][i] + Delta_u[i];
      const double sigma = fuel->calcSigma(lpCls->lp[ip].Tp);
      double du[3];
      FOR_I3 du[i] = lpCls->lp[ip].up[i] - ug[i];
      const double up_m_ug = MAG(du);
      const double urel = up_m_ug + SMALL;
      const double rp =  lpCls->lp[ip].dp/2.;
      // compute the Weber number
      const double weber = rhog*rp*urel*urel/sigma;
      //cout << "> weber: " << weber << " , rhog: " << rhog << " , rp: " << rp << " , urel: " << urel << " , sigma: " << sigma << endl;

// getting rid of very small particles ???
//      if (mp_parent < 1e-015) {
//        P.erase(P.begin()+ip);
//        np -= 1;
//        np_old -= 1;
//        // don't need to increase ip since the next drop is going to be at current ip
//        continue;
//      }

      if (weber>weber_cr) {
        // increase tbu
        // only increase the breakup time associated with the parent particles. The time of the children particles is current and don't need to be updated (they have had their time increased before)
        //IF_RANK0
        //cout << "> weber: " << weber << endl;
        if (ip < np_old) lpCls->lp[ip].tbu += dt;
        const double tbu_parent = lpCls->lp[ip].tbu;

        const double rp_cr = weber_cr*sigma/(rhog*urel*urel);
        const double breakup_time = core_br * sqrt(fuel->calcRho(lpCls->lp[ip].Tp)/rhog)*rp/urel;
        //cout << " > tbu_parent = " << tbu_parent << " , breakup_time = " << breakup_time << endl;
        if (lpCls->lp[ip].tbu >= breakup_time && rp>rp_cr) {
          //cout << " > ----- Breakuptime exceeds, create new particles -----" << endl;
          // create new drops and increase np
          createParticles(ip, Tp_parent, tbu_parent, breakup_time, mp_parent, rhop_parent, rp_parent, rp_cr, weber, du, xp_parent, up_parent);

          // we remove the parent drop,
          lpCls->lp[ip].flag = TRASH;
          //lpCls->lp[ip].icv = -lpCls->lp[ip].icv - 1;
          //np -= 1;
          //np_old -= 1;
        }
      } else {
        lpCls->lp[ip].tbu = 0.0;
      }
      // go to the next particle
      ++ip;
    }

    // eliminates the particles with negative icv
    recycleLsp(false);
    //recycleParticles();

  }

  void recycleParticles() {

    int ip_new = 0;
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      if (lpCls->lp[ip].icv >= 0) {
        if (ip_new != ip) {
          lpCls->lp[ip_new].copy(lpCls->lp[ip]);
        }
        ip_new++;
      }
    }
    lpCls->resize(ip_new);
  }

  void createParticles(const int &ipp, const double &Tp_parent, const double &tbu_parent, const double &breakup_time, const double &mp_parent, const double &rho_parent, const double &rp_parent, const double &rp_critical, const double &weber, const double* up_m_ug, const double* xp_parent, const double* up_parent) {

    // parent droplet properties...
    const double npar_parent = lpCls->lp[ipp].npar;
    const double npar_target = npar_target_multiplier*npar_parent;

    // local parameters:
    const double X_GRID_RATIO = 1.2; // power for CDF distribution
    const int NGRID = 100;
    const int NP_CHILD_MAX = 300;

    // create the pdf of particle radius - some newer SVA settings?...
    //const double K1 = 1.0;
    //const double K2 = 1.0;
    //const double X_si  = K1*log(weber_cr/weber);
    //const double X_si2 = -X_si/(K2*log(rp_parent/rp_critical));
    //double X_si2 = -0.015*log(pow((lspb.weber_critical/weber),0.33));
    //double X_si  = X_si2*X_si_over_X_si2;
    //double X_si  = -0.4;

    // settings from liquid jet simulation (2001)...
    // also had constant values X_si = -0.7, X_si2 = 0.1
    // use these:
    //const double X_si          = -0.1;
    //const double X_si2         = -0.1*log(weber_cr/weber);
    //const double nu_t          = 1.0;
    //const double X_si_nut      = X_si*nu_t;
    //const double X_si2_nut_inv = 1.0/sqrt(2.0*X_si2*nu_t);
    //const double Xd_max        = log(rp_parent);

    //const double K1 = bu_k1;
    //const double K2 = bu_k2;
    //const double X_si  = K1*log(weber_cr/weber);
    //const double X_si2 = -X_si/(K2*log(rp_parent/rp_critical));
    //const double X_si  = -0.4;
    //const double X_si2 = -0.015*log(pow((weber_cr/weber),0.33));
    const double X_si  = bu_k1;
    const double X_si2 = bu_k2*log(pow((weber_cr/weber),0.33));
    const double X_si2_sqrt_inv = 1.0/sqrt(2.0*X_si2);
    const double Xd_max = log(rp_parent);

    // Now construct the CDF over a range from rp_min = rp_critical/2 to rp_parent. This
    // will allow (at some small probability) drops near the critical radius to break
    // into reasonably sized children, rather than leaving the lower threshold as
    // rp_critical.
    double Rd[NGRID], Fp[NGRID];
    for (int i = 0;i < NGRID; ++i) {
      Rd[i]           =  rp_critical/2.0 + (rp_parent - rp_critical/2.0)*pow(((double)(i)/(double)(NGRID-1)),X_GRID_RATIO);
      const double Xd       = log(Rd[i]);
      const double err_func = error_function((Xd - Xd_max - X_si)*X_si2_sqrt_inv);
      Fp[i]           = 0.5*(1.0 + err_func);
    }

    // normalize the CDF...
    const double Fp_min = Fp[0];
    const double Fp_max = Fp[NGRID-1];
    if (Fp_min >= Fp_max) {
      cout << " > Error: CDF problems: weber, rp_crit, rp_parent, Fp_min, Fp_max "
     <<weber<<" "<<rp_critical<<" "<<rp_parent<<" "<<Fp_min<<" "<<Fp_max<<endl;
    }
    const double tmp = 1.0/(Fp_max - Fp_min);
    for (int i = 0;i < NGRID; ++i) {
      Fp[i] = (Fp[i] - Fp_min)*tmp;
    }

    double mass_children = 0.0;
    int np_child = 0;

    //child loop
    double dp_child[NP_CHILD_MAX];
    while (np_child < NP_CHILD_MAX) {

      int ii;
      double rand_val = (double)rand()/(double)RAND_MAX;
      for (int i = 0;i < NGRID; ++i) { //find_bin_parcel
        ii = i;
        if ((Fp[i] >= rand_val)||(i==(NGRID-1))) {
          break;
        }
      }

      const double dp_child_temp = 2.0*Rd[ii];
      const double mp_child_temp = npar_target*rho_parent*M_PI/6.0*pow(dp_child_temp,3);

      if ((np_child == 0)||(mass_children + mp_child_temp < mp_parent)) {
        // add child...
        dp_child[np_child] = dp_child_temp;
        mass_children += mp_child_temp;
        np_child++;
      }
      else {
        // decide if we should add the next child and then normalize...
        if (fabs(mass_children + mp_child_temp - mp_parent) < fabs(mass_children - mp_parent)) {
          // add it...
          dp_child[np_child] = dp_child_temp;
          mass_children += mp_child_temp;
          np_child++;
        }
        // Normalize the childs to match the parent's mass
        const double temp = pow(mp_parent/mass_children , 1./3.);
        for (int ich = 0; ich < np_child; ++ich) {
          dp_child[ich] *= temp;
        }

        //cout << " > Child generation finished with np_child = " << np_child << endl;
        break;
      }
    }

    //cout << "> dp_parent: " << 2.*rp_parent << " , weber: " << weber << " , dp_critical: " << 2.*rp_critical << endl;
    //cout << "> tbu_parent: " << tbu_parent << " , breakup_time: " << breakup_time << endl;
    //cout << "> chile dp: " << endl;
    //for (int i = 0; i < np_child; ++i)
    //  cout << "> " << dp_child[i] << endl;

    //getchar();

    if (np_child <= 0) {
      cout <<" > Error: expect some children." << endl;
      throw(-1);
    }

    int ip_child_f = lpCls->size();
    int ip_child_l = lpCls->size()+np_child;

    //cout << " > ip_child_f = " << ip_child_f << ", ip_child_l = " << ip_child_l << endl;

    // resize the particle vector
    lpCls->resize(ip_child_l);

    // set particle data and compute actual mass...
    mass_children = 0.0;
    for (int ip_child = ip_child_f; ip_child <ip_child_l; ++ip_child) {
      lpCls->lp[ip_child].icv = lpCls->lp[ipp].icv;
      lpCls->lp[ip_child].Tp = Tp_parent;
      // children take their parent's original position so fluxes
      // get properly computed across any stats planes...

      const double dp_child_temp = dp_child[ip_child-ip_child_f];
      lpCls->lp[ip_child].dp = dp_child_temp;
      lpCls->lp[ip_child].mp = npar_target*rho_parent*M_PI/6.0*pow(dp_child_temp,3);
      mass_children += lpCls->lp[ip_child].mp;

      // ===================================
      // setting droplet auxillary variables
      // ===================================
      // for setting the child's breakup time, first
      // check the child's weber number. If the child is
      // above the critical weber number, add the remaining
      // time from the parent to the child...
      // need to modify this...
      //lpCls->lp[ip_child].tbu = 0.0;
      if (weber*dp_child_temp/(2.0*rp_parent) > weber_cr) {
        lpCls->lp[ip_child].tbu = tbu_parent - breakup_time;
      }
      else {
        lpCls->lp[ip_child].tbu = 0.0;
      }

      // set the flag for recycling
      lpCls->lp[ip_child].flag = KEEP;
      // set the parceling param
      lpCls->lp[ip_child].npar = npar_target;

    }

    // check mass conservation...
    if (fabs(mass_children-mp_parent)/mp_parent > 1.0E-8) {
      cout << " > Error: mass_children: " << mass_children << " , mp_parent: " << mp_parent <<
              " , error: " << fabs(mass_children-mp_parent)/mp_parent << " , np_child: " << np_child << endl;
      assert(0);
    }

    // get a direction normal to the relative velocity
    // add random components to new droplet positions and
    // velocity. don't worry about the particle moving out
    // of the current cell because of these additions
    // it will be automatically taken care of in the next
    // iteration...

    double dir_cos[NP_CHILD_MAX][3];
    calcDirectionsNormaltoRelVel(np_child,up_m_ug,dir_cos);

    // set position and momentum for each particle
    double momentum_children[3]; FOR_I3 momentum_children[i] = 0.0;
    double momentum_parent[3];   FOR_I3 momentum_parent[i]   = mp_parent*up_parent[i];
    for (int ip_child = ip_child_f; ip_child < ip_child_l; ++ip_child) {

      // introduce randomness
      //double rand_val = (double)rand()/(double)RAND_MAX;
      double rand_val = 1.0;

      // an appropriate velocity scale is assumed to be the
      // child diameter / the parent breakup time.
      //double vel_norm = dp_child[ip_child-ip_child_f]/ breakup_time;
      double vel_norm = rp_parent / breakup_time * 2.;

      // position
      // leave the x intact because the velocity is perturbed and
      // in the next time step, chrlderen will have different positions
      FOR_I3 lpCls->lp[ip_child].xp[i] = xp_parent[i];
      // =================================================================================================
      // you can add a random number to the position of the particle but it may result in different icv!!!
      // =================================================================================================
      //const double dx = dp_child[ip_child-ip_child_f];
      //FOR_I3 lspVec[ip_child].xp[i] = xp_parent[i] + rand_val*dx*dir_cos[ip_child-ip_child_f][i];
      // set child velocity...
      FOR_I3 {
        lpCls->lp[ip_child].up[i] = up_parent[i] + rand_val*vel_norm*dir_cos[ip_child-ip_child_f][i];
        momentum_children[i] += lpCls->lp[ip_child].mp * lpCls->lp[ip_child].up[i];
      }
    }
    // during breakup, momentum should be preserved
    double momentum_res[3];
    FOR_I3 momentum_res[i] = momentum_parent[i] - momentum_children[i];
    // now correct velocity to have the momentum balance
    double velocity_cor[3]; FOR_I3 velocity_cor[i] = momentum_res[i] / mp_parent;
    for (int ip_child = ip_child_f; ip_child < ip_child_l; ++ip_child) {
      FOR_I3 lpCls->lp[ip_child].up[i] += velocity_cor[i];
    }

    //cout << " > np = " << np << endl;
    //cout << " > ----- Finish create particle with np = " << np << endl;
  }

  double error_function(const double &eta) {

    double z = fabs(eta);
    double t = 1.0/(1.0+0.5*z);
    double erfcc = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+
                t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+
                          t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
    if (eta < 0.0) erfcc = 2.0-erfcc;

    double err_func = 1.0 - erfcc;
    return err_func;

  }

  void calcDirectionsNormaltoRelVel(const int n,const double *vec,double (*dir_cos)[3]) {

    // given a vector, we want to find a basis for the
    // plane normal to the vector

    // copy down the vector
    double V[3]; FOR_I3 V[i] = vec[i];

    // we assume that the vector is non-zero
    {
      double tmp = 0;
      FOR_I3 tmp += V[i]*V[i];
      assert(tmp>0.0);
      double one_over_tmp = 1.0/sqrt(tmp);
      FOR_I3 V[i] *= one_over_tmp;
    }

    // wa want to build the orthonormal basis (V,e1,e2)
    double e1[3], e2[3];
    // start by choosing a vector that is guaranteed not aligned with vec...
    if ( fabs(V[0]) <= min(fabs(V[1]),fabs(V[2])) ) {
      e2[0] = 1.0;
      e2[1] = 0.0;
      e2[2] = 0.0;
    }
    else if ( fabs(V[1]) <= min(fabs(V[0]),fabs(V[2])) ) {
      e2[0] = 0.0;
      e2[1] = 1.0;
      e2[2] = 0.0;
    }
    else {
      e2[0] = 0.0;
      e2[1] = 0.0;
      e2[2] = 1.0;
    }
    // the first vector is the cross of this and vec...
    e1[0] = e2[1]*V[2] - e2[2]*V[1];
    e1[1] = e2[2]*V[0] - e2[0]*V[2];
    e1[2] = e2[0]*V[1] - e2[1]*V[0];
    // normalizing the vector is necessary since V and e2
    // are not necessarily orthonormal
    {
      double tmp = 0.0;
      FOR_I3 tmp += e1[i]*e1[i];
      assert(tmp>0.0);
      double one_over_tmp = 1.0/sqrt(tmp);
      FOR_I3 e1[i] *= one_over_tmp;
    }
    // the second vector is the cross product of e1 and V...
    e2[0] = e1[1]*V[2] - e1[2]*V[1];
    e2[1] = e1[2]*V[0] - e1[0]*V[2];
    e2[2] = e1[0]*V[1] - e1[1]*V[0];
    // normaling the vector is not necessary since V and e2
    // are orthonormal ... however normalize to reduce roundoff
    // errors
    {
      double tmp = 0.0;
      FOR_I3 tmp += e2[i]*e2[i];
      assert(tmp>0.0);
      double one_over_tmp = 1.0/sqrt(tmp);
      FOR_I3 e2[i] *= one_over_tmp;
    }

    // now we can build a set of n vectors in the plane
    // orthogonal to V
    //srand(2);
    for (int j = 0; j < n ; ++j) {
      double theta         = 2.0*M_PI*rand()/(double)RAND_MAX;  // 0<theta<2pi
      double cos_theta     = cos(theta);
      double sin_theta     = sin(theta);
      FOR_I3 dir_cos[j][i] = cos_theta*e1[i] + sin_theta*e2[i];
    }

  }

  int getLspSize() {
    int mynp = lpCls->size();
    int np;
    MPI_Allreduce(&mynp,&np,1,MPI_INT,MPI_SUM,mpi_comm);
    return np;
  }

  void reportLspSize() {
    int np = getLspSize();
    IF_RANK0 cout << " > lsp size: " << np << endl;
  }

  void computePDF(const double xmin, const double xmax) {
    if (mpi_size > 1) {
      return;
    }
    int NGrid = 101;
    double rmin = 1e-7;
    double rmax = 3e-4;
    double* dp_ar = new double[NGrid];

    double* ddp_ar = new double[NGrid-1];
    double* PDF_ar = new double[NGrid-1];
    double mass_sum = 0.0;

    for (int iar = 0; iar < NGrid; ++iar) {
      dp_ar[iar] = rmin + (double)iar/(double)(NGrid-1)*(rmax - rmin);
    }

    for (int iar = 0; iar < NGrid-1; ++iar) {
      ddp_ar[iar] = dp_ar[iar+1] - dp_ar[iar];
      PDF_ar[iar] = 0.0;
    }

    for (int ip = 0; ip < lpCls->size(); ++ip) {
      if ((lpCls->lp[ip].xp[0]>=xmin) && (lpCls->lp[ip].xp[0]<xmax)) {
        double dp_temp = lpCls->lp[ip].dp;
        double mp_temp = lpCls->lp[ip].mp;
        for (int i = 0; i < NGrid-1; ++i) {
          if ((dp_temp >= dp_ar[i]) && (dp_temp < dp_ar[i+1])) {
            PDF_ar[i] += mp_temp;
            mass_sum += mp_temp;
            break;
          }
        }
      }
    }

    for (int iar = 0; iar < NGrid-1; ++iar) {
      //cout << " > mass_sum = " << mass_sum << ", ddp_ar[iar] = " << ddp_ar[iar] << endl;
      PDF_ar[iar] /= (mass_sum*ddp_ar[iar])*1e6;
    }

    FILE* fid_pdf;
    fid_pdf = fopen("Droplet_Diam_PDF.txt","w");
    for (int i = 0; i < NGrid-1; ++i) {
      double dp_dump = dp_ar[i] * 1e6; // in micrometers
      fprintf(fid_pdf,"%10.6e %10.6e\n",dp_dump,PDF_ar[i]);
    }

    delete[] dp_ar;
    delete[] ddp_ar;
    delete[] PDF_ar;

  }

  void initLspStats() {

    if (mpi_rank == 0)
      cout << "initLspStats()" << endl;

    FOR_PARAM_MATCHING("LSP.STATS") {
      // add a new stats...
      lspStatsList.push_back(newLspStats(&(*param)));
    }

    FOR_PARAM_MATCHING("LSP.PDF") {
      // add a new pdf...
      lspPDFList.push_back(LspPDF());
      LspPDF * stats = &(lspPDFList.back());
      stats->init(&(*param));
      stats->report();
    }

  }

  void updateLspStats(const LspState * lsp,const int np,const double time,const double dt,const int step,const bool verbose) {

    if ((mpi_rank == 0) && verbose && (step%check_interval == 0))
      cout << " > updateLspStats()" << endl;

    for (list<LspStats*>::iterator stats = lspStatsList.begin(); stats != lspStatsList.end(); ++stats) {
      (*stats)->update(lsp,np,material,dt,verbose);
      if (step%((*stats)->getInterval()) == 0)
        (*stats)->write(time,step,verbose);
    }

    for (list<LspPDF>::iterator pdfs = lspPDFList.begin(); pdfs != lspPDFList.end(); ++pdfs) {
      if ((step%pdfs->getInterval() == 0) && (step>=pdfs->getStart())) {
        pdfs->update(lsp,np,dt,verbose);
        pdfs->write(time,step,verbose);
      }
    }

  }

  void resetLspStats(const bool verbose) {

    if ((mpi_rank == 0)&&verbose)
      cout << "resetLspStats()" << endl;

    for (list<LspStats*>::iterator stats = lspStatsList.begin(); stats != lspStatsList.end(); ++stats) {
      (*stats)->reset(verbose);
    }
  }

  void reportVolMassFrac() {

    if (step%check_interval == 0) {

      if ( getLspSize() > 0) {

        double max_vol_r = -1e20;
        double min_vol_r = 1e20;
        double max_mas_r = -1e20;
        double min_mas_r = 1e20;
        FOR_ICV {
          double sum_vol = 0.0;
          double sum_mas = 0.0;
          for (int ip = paocv_i[icv]; ip != paocv_i[icv+1]; ++ip) {
            double pa_vol = M_PI/6.0*pow(lpCls->lp[ip].dp,3);
            sum_vol += pa_vol;
            sum_mas += pa_vol*material->calcRho(lpCls->lp[ip].Tp);
          }
          double vol_r = sum_vol/vol_cv[icv];
          double mas_r = sum_mas/(vol_cv[icv]*rho[icv]);
          max_vol_r = max(max_vol_r,vol_r);
          min_vol_r = min(min_vol_r,vol_r);
          max_mas_r = max(max_mas_r,mas_r);
          min_mas_r = min(min_mas_r,mas_r);
        }
        double max_VR, min_VR, max_MR, min_MR;
        MPI_Reduce(&max_vol_r,&max_VR,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
        MPI_Reduce(&min_vol_r,&min_VR,1,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
        MPI_Reduce(&max_mas_r,&max_MR,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
        MPI_Reduce(&min_mas_r,&min_MR,1,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
        IF_RANK0 {
          cout << " > particle volume ratio range " << min_VR <<  " : " << max_VR << endl;
          cout << " > particle mass ratio range " << min_MR <<  " : " << max_MR << endl;
        }
      }
    }
  }

  template<class T,class S> // T: lp container, S: lp state
  void relocateLsp(T * lpContainer) {
    
    if ((step%check_interval==0)&&(mpi_rank == 0)) {
      cout << "relocateLsp()";
      cout.flush();
    }

    bool first = true;
    MPI_Request * sendRequestArray = NULL;
    MPI_Request * recvRequestArray = NULL;
    double * send_buf = NULL;
    int send_buf_size = 0;
    double * recv_buf = NULL;
    int recv_buf_size = 0;
    int * send_count = new int[mpi_size];
    int * send_disp = NULL; // only allocate if/when needed
    int ip_start = 0;
    int np_unpack,recvRequestCount,sendRequestCount;
    double per_t[3],per_R[9];

    int done = 0;
    while (done == 0) {

      if ((step%check_interval==0)&&(mpi_rank == 0)) {
        cout << ".";
        cout.flush();
      }
      
      // assume we are done...
      int my_done = 1;

      // use send_count to record the number of particles that need to be 
      // sent to a particular rank...
      FOR_RANK send_count[rank] = 0;
      
      // loop through particles. everyone with a positive icv has not
      // yet been located...
      for (int ip = ip_start; ip < lpContainer->np; ++ip) {
        while (lpContainer->lp[ip].icv >= 0) {
          const int icv = lpContainer->lp[ip].icv;
          // include a small reduction to preference the particles staying in their current cv and prevent 
          // an iterative exchange that does not terminate due to errors in periodic transforms, for example... 
          double d2_min = DIST2(lpContainer->lp[ip].xp,x_vv[icv])*0.999999;
          int icv_nbr_min = -1;
          for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) { // note we skip self here
            const int icv_nbr = cvocv_v[coc];
            const double d2_nbr = DIST2(lpContainer->lp[ip].xp,x_vv[icv_nbr]);
            if (d2_nbr < d2_min) {
              icv_nbr_min = icv_nbr;
              d2_min = d2_nbr;
            }
          }
          if (icv_nbr_min == -1) {
            // we must be in the right cv. Use -1 indexing to lock this in place...
            lpContainer->lp[ip].icv = -icv-1;
          }
          else if (icv_nbr_min < ncv) {
            // we got transfered to a closer local cv. That means we are not locally done...
            lpContainer->lp[ip].icv = icv_nbr_min;
          }
          else {
            // we got transfered to a ghost cv. Use -1 indexing in the cv to indicate this...
            assert(icv_nbr_min < ncv_g);
            lpContainer->lp[ip].icv = -icv_nbr_min-1;
            //lpContainer->lp[ip].d2 = d2_min; // HACK - store the min distance to check
            // and record the processor that will be receiving this particle...
            int rank,bits,index;
            BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr_min-ncv]);
            ++send_count[rank];
            my_done = 0;
          }
        }
      }
      
      // at this point, we have particles that are properly contained in local 
      //cvs, possibly particles in
      // ghost cvs...
      
      MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);
      if (done == 0) {
        
        // someone has particles in their ghost data, so we all have to exchange... 

        // on the first time, we allocate the send and recv request arrays to the 
        // same size as the cv-communicator. These are the cvs where the ip's have been 
        // located...

        if (first) {
          assert(send_disp == NULL);
          send_disp = new int[mpi_size];
          int requestCount = 0;
          for (int ii = 0; ii < cvPrcommVec.size(); ii++) {
            if ((cvPrcommVec[ii].rank != mpi_rank)) { 
              requestCount++;
            }
          }
          assert(sendRequestArray == NULL); sendRequestArray = new MPI_Request[requestCount];
          assert(recvRequestArray == NULL); recvRequestArray = new MPI_Request[requestCount];
        }
        
        // now post the non-blocking send and recv's for the single int, the count. At the same
        // time prepare the offsets...
        
        int requestCount = 0;
        int np_pack = 0;
        for (int ii = 0; ii < cvPrcommVec.size(); ii++) {
          cvPrcommVec[ii].nunpack_v = -1; // for checking
          if ((cvPrcommVec[ii].rank != mpi_rank)) { 
            // non-local: send number of parcels to receiver...
            MPI_Irecv(&(cvPrcommVec[ii].nunpack_v),1,MPI_INT,cvPrcommVec[ii].rank,
                      1812,mpi_comm,recvRequestArray+requestCount);
            cvPrcommVec[ii].npack_v = send_count[cvPrcommVec[ii].rank];
            MPI_Issend(&(cvPrcommVec[ii].npack_v),1,MPI_INT,cvPrcommVec[ii].rank,
                       1812,mpi_comm,sendRequestArray+requestCount);
            send_disp[cvPrcommVec[ii].rank] = np_pack;
            np_pack += cvPrcommVec[ii].npack_v;
            requestCount++;
          }
          else {
            // local (must be self-periodic)...
            cvPrcommVec[ii].npack_v = cvPrcommVec[ii].nunpack_v = send_count[mpi_rank];
            send_disp[mpi_rank] = np_pack;
            np_pack += cvPrcommVec[ii].npack_v; // include the self-periodic in np_pack, because we need a buffer to transform into
          }
        }

        if (first) {
          assert(send_buf == NULL);
          send_buf_size = np_pack*S::data_size();
          send_buf = new double[send_buf_size];
        }
        else if (np_pack*S::data_size() > send_buf_size) {
          delete[] send_buf;
          send_buf_size = np_pack*S::data_size();
          send_buf = new double[send_buf_size];
        }

        // pack using offset in send_disp...
        for (int ip = ip_start; ip < lpContainer->np; ++ip) {
          // look for particles in the ghost region...
          if (lpContainer->lp[ip].icv <= -ncv-1) {
            const int icv = -lpContainer->lp[ip].icv-1;
            int rank,bits,index;
            BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv-ncv]);
            
            bool b_per_t = false;
            bool b_per_R = false;
            if (bits) {
              b_per_t = PeriodicData::getPeriodicT(per_t,bits);
              b_per_R = PeriodicData::getPeriodicR(per_R,bits);
            }
            
            lpContainer->lp[ip].icv = index; // put the icv on rank in icv for packing...
            if (b_per_R && b_per_t) {
              lpContainer->lp[ip].packRt(send_buf+send_disp[rank]*S::data_size(),per_R,per_t);
            }
            else if (b_per_R) {
              lpContainer->lp[ip].packR(send_buf+send_disp[rank]*S::data_size(),per_R);
            }
            else if (b_per_t) {
              lpContainer->lp[ip].packt(send_buf+send_disp[rank]*S::data_size(),per_t);
            }
            else {
              lpContainer->lp[ip].pack(send_buf+send_disp[rank]*S::data_size());
            }
            lpContainer->lp[ip].icv = -icv-1; // put back icv so we discard below...
            ++send_disp[rank];
          }
        }
        
        // once packed, complete the send/recv of the counts, and post the send/recv of the data...
        
        if (requestCount > 0) {
          MPI_Waitall(requestCount,recvRequestArray,MPI_STATUSES_IGNORE);
          MPI_Waitall(requestCount,sendRequestArray,MPI_STATUSES_IGNORE);
        }
        
        np_unpack = 0;
        for (int ii = 0; ii < cvPrcommVec.size(); ii++) {
          assert(cvPrcommVec[ii].nunpack_v >= 0); // check
          np_unpack += cvPrcommVec[ii].nunpack_v; // include self-periodic here
        }
        
        // and allocate the recv_buf...
        if (first) {
          assert(recv_buf == NULL);
          recv_buf_size = np_unpack*S::data_size();
          recv_buf = new double[recv_buf_size];
          // this is the last first...
          first = false;
        }
        else if (np_unpack*S::data_size() > recv_buf_size) {
          delete[] recv_buf;
          recv_buf_size = np_unpack*S::data_size();
          recv_buf = new double[recv_buf_size];
        }
          
        // now post the sends and recvs, but ONLY when they are finite...

        np_pack = 0;
        np_unpack = 0;
        recvRequestCount = 0;
        sendRequestCount = 0;
        for (int ii = 0; ii < cvPrcommVec.size(); ii++) {
          if ((cvPrcommVec[ii].rank != mpi_rank)) { 
            if (cvPrcommVec[ii].nunpack_v > 0) {
              MPI_Irecv(recv_buf+np_unpack*S::data_size(),cvPrcommVec[ii].nunpack_v*S::data_size(),MPI_DOUBLE,cvPrcommVec[ii].rank,
                        1813,mpi_comm,recvRequestArray+recvRequestCount);
              ++recvRequestCount;
              np_unpack += cvPrcommVec[ii].nunpack_v;
            }
            if (cvPrcommVec[ii].npack_v > 0) {
              MPI_Issend(send_buf+np_pack*S::data_size(),cvPrcommVec[ii].npack_v*S::data_size(),MPI_DOUBLE,cvPrcommVec[ii].rank,
                         1813,mpi_comm,sendRequestArray+sendRequestCount);
              ++sendRequestCount;
              np_pack += cvPrcommVec[ii].npack_v;
            }
          }
          else if (cvPrcommVec[ii].nunpack_v > 0) {
            // self-periodic: local copy from send_buf into recv_buf...
            assert(cvPrcommVec[ii].nunpack_v == cvPrcommVec[ii].npack_v);
            memcpy(recv_buf+np_unpack*S::data_size(),send_buf+np_pack*S::data_size(),sizeof(double)*cvPrcommVec[ii].npack_v*S::data_size());
            np_unpack += cvPrcommVec[ii].nunpack_v;
            np_pack += cvPrcommVec[ii].npack_v;
          }
        }
      
      } // if (done == 0)
      
      // ================================================================
      // with everything posted, take some time to do some work on 
      // any of the particles that have been properly located. This means 
      // bcs and compression...
      // ================================================================
        
      for (int ip = ip_start; ip < lpContainer->np; ++ip) {
        assert(lpContainer->lp[ip].icv < 0);
        // check if we are in a local cv...
        if (lpContainer->lp[ip].icv > -ncv-1) {
          const int icv = -lpContainer->lp[ip].icv-1;
          lpContainer->lp[ip].icv = icv; // flip back to valid
          // check for boundary faces...
          if (bfocv_i[icv+1] == bfocv_i[icv]) {
            // no boundary faces, so keep this one!...
            if (ip_start != ip) lpContainer->lp[ip_start] = lpContainer->lp[ip];
            ++ip_start;
          }
          else {
            // this one has boundary faces. We need a length scale...
            const double tol2 = 1.0E-10*pow(vol_cv[icv],2.0/3.0);
            // and we need the closest ibf,ist_ss...
            int ibf_closest = -1;
            //int ist_ss_closest;
            int zone_closest;
            double d2_closest,dx_closest[3],normal_closest[3],dist_closest;
            for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
              const int ibf = bfocv_v[boc];
              for (int sob = sstobf_i[ibf]; sob != sstobf_i[ibf+1]; ++sob) {
                const int ist_ss = sstobf_v[sob];
                // the corner coords of this sub-surface tri are (recall that the
                // subsurface tri already has any periodic transfrom taken care of)...
                const double * const x0 = subSurface->xp[subSurface->spost[ist_ss][0]];
                const double * const x1 = subSurface->xp[subSurface->spost[ist_ss][1]];
                const double * const x2 = subSurface->xp[subSurface->spost[ist_ss][2]];
                double xp[3]; getClosestPointOnTriRobust(xp,lpContainer->lp[ip].xp,x0,x1,x2);
                const double dx[3] = DIFF(xp,lpContainer->lp[ip].xp);
                const double d2 = DOT_PRODUCT(dx,dx);
                if ((ibf_closest == -1)||(d2 < d2_closest+tol2)) {
                  const double normal[3] = TRI_NORMAL_2(x0,x1,x2);
                  const double normal_mag = MAG(normal);
                  assert(normal_mag != 0); // the subsurface should not have zero-area tris!
                  const double dist = DOT_PRODUCT(dx,normal)/normal_mag;
                  if ((ibf_closest == -1)||(d2 < d2_closest-tol2)||(fabs(dist) > fabs(dist_closest))) {
                    ibf_closest = ibf;
                    //ist_ss_closest = ist_ss;
                    d2_closest = d2;
                    FOR_I3 dx_closest[i] = dx[i];
                    FOR_I3 normal_closest[i] = normal[i];
                    dist_closest = dist;
                    zone_closest = zone_bf[ibf];
                    assert(zone_closest < bfZoneVec.size());
                  }
                }
              }
            }
            assert(ibf_closest >= 0);
            if (dist_closest <= 0.0) {
	      // dist with negative sign means we are outside the fluid volume, so apply bc...
              //const int izn = zone_bf[ibf_closest];
              //assert((izn >= 0)&&(izn < lpContainer->bcVec.size()));
              //const bool keep = lpContainer->bcVec[izn]->applyBc(lpContainer->lp[ip],ibf_closest,ist_ss_closest,dx_closest,normal_closest);
              const bool keep = lpZoneBcVec[zone_closest]->applyBc2(lpContainer->lp[ip], dx_closest, normal_closest);
              if (keep) {
                // bc says keep it, so we keep it...
                if (ip_start != ip) lpContainer->lp[ip_start] = lpContainer->lp[ip];
                ++ip_start;
              }
            }
            else {
              // we are inside the boundary, so keep...
              if (ip_start != ip) lpContainer->lp[ip_start] = lpContainer->lp[ip];
              ++ip_start;
            }
          }
        }
      }

      // prepare to recv the sent ip's into the end of the lp...
      if (done == 0) {
        
        // we are expecting np_unpack particles...
        lpContainer->resize(ip_start+np_unpack);
        
        // wait for send/recv pairs to complete...
        if (recvRequestCount > 0) MPI_Waitall(recvRequestCount,recvRequestArray,MPI_STATUSES_IGNORE);
        if (sendRequestCount > 0) MPI_Waitall(sendRequestCount,recvRequestArray,MPI_STATUSES_IGNORE);

        // and unpack...
        int offset = 0;
        for (int ip = ip_start; ip < lpContainer->np; ++ip) {
          lpContainer->lp[ip].unpack(recv_buf+offset);
          assert((lpContainer->lp[ip].icv >= 0)&&(lpContainer->lp[ip].icv < ncv));
          // HACK: check d2 with our nbr: should match d2_min stored in mp...
          //const double d2 = DIST2(lpContainer->lp[ip].xp,x_vv[lpContainer->lp[ip].icv]);
          //assert(d2 == lpContainer->lp[ip].d2);
          offset += S::data_size();
        }

      }
      else {

        // we are done: resize to handle any compression...
        lpContainer->resize(ip_start);
        
      }
      
    } // while (done == 0)
    
    if ((step%check_interval==0)&&(mpi_rank == 0))
      cout << "OK" << endl;
    
    if (sendRequestArray != NULL) delete[] sendRequestArray; 
    if (recvRequestArray != NULL) delete[] recvRequestArray; 
    if (send_buf != NULL) delete[] send_buf;
    if (recv_buf != NULL) delete[] recv_buf;
    delete[] send_count;
    if (send_disp != NULL) delete[] send_disp;
    
  }

  template<class T>
  void relocateLsp_rayTrace() {

    if ((mpi_rank == 0) && (step%check_interval == 0))
      cout << " > relocateLsp_rayTrace";

    //XXX
    //assert(surface);

    //double tol = 0.0;

    int * send_count = new int[mpi_size];
    int * send_disp = new int[mpi_size];
    int * recv_count = new int[mpi_size];
    int * recv_disp = new int[mpi_size];
    double per_t[3],per_R[9];

    int done = 0;
    while (done == 0) {

      if ((mpi_rank == 0) && (step%check_interval == 0)) {
        cout << ".";
        cout.flush();
      }

      FOR_RANK send_count[rank] = 0; // the number of particles we need to send to each rank

      int my_done = 1;

      // hack to get closest... function should return ip_debug for the right rank and for the
      // other ranks will return -1
      //double xp_debug[3] = { 0.0752, -0.0049, 0.003 };
      //double xp_debug[3] = { 0.0739, -0.005, 0.0038 };
      //int ip_debug = searchClosestIpX(xp_debug);

      for (int ip = 0; ip < lpCls->size(); ++ip) {

        //bool debug = false;
        //if (ip == ip_debug && step >144600 && step <= 144649)
        //if (ip == ip_debug && step >145900 && step <= 146100)
        //  debug = false;

        // this particle is currently in cv "icv"...
        // we use -ve icv indiexing top indicate that a particle has been
        // properly located, and does not need to be revisited...

        const int icv = lpCls->lp[ip].icv;
        if (icv >= 0) {

          //if (icv >= ncv)
          //  cout << "XXXXXX: " << icv << " " << ncv << endl;
          assert(icv < ncv);

          // this particle still needs to be checked...
          // and has followed a trajectory from xp0 to xp...

          double * xp0 = lpCls->lp[ip].xp0;
          double * xp  = lpCls->lp[ip].xp;

          //if (debug) cout << "debug: " << " ,ip: " << ip << " ,rank: " << mpi_rank << " xp0: " << COUT_VEC(xp0) << " xp: " << COUT_VEC(xp) << endl;

          // if this cell has boundaries, then check for crossing...

          double xi_closest[3];
          double st_normal_closest[3];
          int ist_closest = -1;
          double d2_closest = HUGE_VAL;
          int zone_closest;
          for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
            const int ibf = bfocv_v[boc];
            for (int sob = sstobf_i[ibf]; sob != sstobf_i[ibf+1]; ++sob) {
              const int ist = sstobf_v[sob];
              // the corner coords of this surface tri are...
              const double * const x0 = subSurface->xp[subSurface->spost[ist][0]];
              const double * const x1 = subSurface->xp[subSurface->spost[ist][1]];
              const double * const x2 = subSurface->xp[subSurface->spost[ist][2]];
              // did we cross this tri?...
              const double st_normal[3] = TRI_NORMAL_2(x0,x1,x2); // should be an outward-pointing normal (wrt flow)
              // get the side associated with both points...
              const double d0 =
                (xp0[0]-x1[0])*st_normal[0] +
                (xp0[1]-x1[1])*st_normal[1] +
                (xp0[2]-x1[2])*st_normal[2];
              const double d1 =
                (xp[0]-x1[0])*st_normal[0] +
                (xp[1]-x1[1])*st_normal[1] +
                (xp[2]-x1[2])*st_normal[2];
              // NOTE: this requires tolerance
              if ((d0 <= 0)&&(d1 > 0)) {
                // this is a crossing!...
                //if (debug) cout << "GOT POTENTIAL CROSSING! icv: " << icv << " d0: " << d0 << " d1: " << d1 << " st_normal: " << COUT_VEC(st_normal) << " xp: " << COUT_VEC(xp) << " step: " << step << endl;
                // compute xi...
                double xi[3]; FOR_I3 xi[i] = (d1*xp0[i] - d0*xp[i])/(d1-d0);
                const double n01[3] = TRI_NORMAL_2(x0,x1,xi);
                //if (debug) cout << "n01: " << DOT_PRODUCT(st_normal,n01) << " step: " << step << endl;
                if (DOT_PRODUCT(st_normal,n01) >= 0.0) {
                  const double n12[3] = TRI_NORMAL_2(x1,x2,xi);
                  //if (debug) cout << "n12: " << DOT_PRODUCT(st_normal,n12) << " step: " << step << endl;
                  if (DOT_PRODUCT(st_normal,n12) >= 0.0) {
                    const double n20[3] = TRI_NORMAL_2(x2,x0,xi);
                    //if (debug) cout << "n20: " << DOT_PRODUCT(st_normal,n20) << " step: " << step << endl;
                    if (DOT_PRODUCT(st_normal,n20) >= 0.0) {
                      // xi is inside the tri...
                      // we want the point that is closes to xp0 -- i.e. the starting point of the particle.
                      const double d2 = DIST2(xp0,xi);
                      //if (debug) cout << "THIS WAS AN ACTUAL CROSSING with dist: " << sqrt(d2) << endl;
                      if ((ist_closest == -1)||(d2 < d2_closest)) {
                        d2_closest = d2;
                        ist_closest = ist;
                        FOR_I3 xi_closest[i] = xi[i];
                        FOR_I3 st_normal_closest[i] = st_normal[i];
                        zone_closest = zone_bf[ibf];
                        assert(zone_closest < bfZoneVec.size());
                      }
                    }
                  }
                }
                //getchar();
              }
            }
          }
          
          // if we got a closest intersection, check if it is closest to us, or
          // instead closer to one of our neighbors...
          
          if (ist_closest != -1) {
            double d2_min = DIST2(xi_closest,x_vv[icv]);
            int icv_min = icv;
            // check if there is a nbr in which we are closer...
            assert(cvocv_v[cvocv_i[icv]] == icv); // self check
            for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) { // note we skip self here
              const int icv_nbr = cvocv_v[coc];
              const double d2_icv_nbr = DIST2(xi_closest,x_vv[icv_nbr]);
              if (d2_icv_nbr < d2_min) {
                d2_min = d2_icv_nbr;
                icv_min = icv_nbr;
              }
            }
            if (icv_min == icv) {
              // apply the bc...
              //if (debug) cout << "REFLECTION: " << " ip: " << ip << " step: " << step << endl;
              double buf[7];
              lpZoneBcVec[zone_closest]->applyBc(lpCls->lp[ip], xi_closest, st_normal_closest, buf);
              //if (b_irregCol) {
              //  irregularBounceLpBc(lpCls->lp[ip],xi_closest,st_normal_closest);
              //} else {
              //  impactLpBc(lpCls->lp[ip],xi_closest,st_normal_closest,tol);
              //}
            }
            else {
              lpCls->lp[ip].icv = icv_min;
              if (icv_min >= ncv) {
                // icv_min is in the ghost range, so we will be sending this lsp to another rank...
                int rank,bits,index;
                BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_min-ncv]);
                ++send_count[rank];
              }
            }
            // in either case we are not done...
            my_done = 0;
          }
          else {

            // step 1: loop through any surface tris that are part of this
            // cv's boundary and find the closest tri intersection: xi. Note
            // that because the tris are the original surface tris, they may
            // extend beyond the icv's local voronoi diagram.
            // If there is no intersection, then just continue to the nbr check below.
            // If there is an intersection, run the nbr check below on the
            // point xi. If the intersection point is closest to icv, simply apply the bc
            // associated with this surface tri (e.g. update position) and continue
            // to point check. If xi is closer to an icv_nbr, send the particle to that nbr,
            // and releat everything, including the calculation of nearest xi (because it
            // might be different due to additional surface geometry in the nbr)...

            // calc the current distance to this point...
            double d2_min = DIST2(lpCls->lp[ip].xp,x_vv[icv]);
            int icv_min = icv;
            // check if there is a nbr in which we are closer...
            assert(cvocv_v[cvocv_i[icv]] == icv); // self check
            for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) { // note we skip self here
              const int icv_nbr = cvocv_v[coc];
              const double d2_icv_nbr = DIST2(lpCls->lp[ip].xp,x_vv[icv_nbr]);
              if (d2_icv_nbr < d2_min) {
                d2_min = d2_icv_nbr;
                icv_min = icv_nbr;
              }
            }
            if (icv_min != icv) {
              lpCls->lp[ip].icv = icv_min;
              if (icv_min >= ncv) {
                // icv_min is in the ghost range, so we will be sending this lsp to another rank...
                int rank,bits,index;
                BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_min-ncv]);
                ++send_count[rank];
              }
              my_done = 0;
            }
            else {
              // we usea negative icv to indicate that this particle is done. i.e.
              // it is in the right place and doe NOT need to be visited again.
              lpCls->lp[ip].icv = -lpCls->lp[ip].icv-1; // flip icv to negative...
            }
          }
        }
      } // for ip...

      // redistibute the lpVec for elements that have passed into ghost cvs...
      // note: we could also remove elements here in the future if required
      // (e.g. evaporated particles). Not sure if theis is the right place...

      ///////////////////////////

      // prepare sending buffers
      map<int, int*> rank_send_buf_int_map;         // maps the rank to the buffer index
      map<int, double*> rank_send_buf_double_map;   // maps the rank to the buffer doubles
      for (int irank = 0; irank < mpi_size; irank++) send_disp[irank] = 0;
      for (int irank = 0; irank < mpi_size; irank++) {
        if (send_count[irank] > 0) {
          rank_send_buf_int_map[irank] = new int [send_count[irank]];
          rank_send_buf_double_map[irank] = new double [send_count[irank]*T::data_size()];
        }
      }

      int np_new = 0;
      for (int ip = 0; ip < lpCls->size(); ++ip) {
        if (lpCls->lp[ip].icv >= ncv) {
          int rank,bits,index;
          BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[lpCls->lp[ip].icv-ncv]);

          map<int, int*>::iterator iter_int = rank_send_buf_int_map.find(rank);
          assert(iter_int != rank_send_buf_int_map.end());

          map<int, double*>::iterator iter_double = rank_send_buf_double_map.find(rank);
          assert(iter_double != rank_send_buf_double_map.end());

          iter_int->second[send_disp[rank]] = index;

          bool b_per_t = false;
          bool b_per_R = false;
          if (bits) {
            b_per_t = PeriodicData::getPeriodicT(per_t,bits);
            b_per_R = PeriodicData::getPeriodicR(per_R,bits);
          }
          //const int ipR = iRt_g[lpCls->lp[ip].icv-ncv][0];
          //const int ipt = iRt_g[lpCls->lp[ip].icv-ncv][1];

          if (b_per_R && b_per_t) {
            // rotate/translate xp/xp0 and rotate up
            lpCls->lp[ip].packRt(iter_double->second+send_disp[rank]*T::data_size(),per_R,per_t);
          }
          else if (b_per_R) {
            // rotate xp/xp0 and rotate up
            lpCls->lp[ip].packR(iter_double->second+send_disp[rank]*T::data_size(),per_R);
          }
          else if (b_per_t) {
            // translate xp/xp0
            lpCls->lp[ip].packt(iter_double->second+send_disp[rank]*T::data_size(),per_t);
          }
          else {
            lpCls->lp[ip].pack(iter_double->second+send_disp[rank]*T::data_size());
          }
          ++send_disp[rank];
        }
        else {
          if (np_new != ip)
            lpCls->lp[np_new].copy(lpCls->lp[ip]);
          ++np_new;
        }
      }

      for (int irank = 0; irank < mpi_size; irank++) {
        assert(send_disp[irank]==send_count[irank]);
      }
      
      // Now find out how many particles we receive 
      const int UPDATE_TAG1 = 121212;
      int requestCount = 0;
      for (int ivec = 0; ivec < cvPrcommVec.size(); ivec++) {
        if ((cvPrcommVec[ivec].rank != mpi_rank)) { 
          requestCount++;
        }
      }

      MpiRequestStuff * mrs = new MpiRequestStuff; 
      mrs->sendRequestArray = new MPI_Request[requestCount];
      mrs->recvRequestArray = new MPI_Request[requestCount];

      requestCount = 0;
      for (int ivec = 0; ivec < cvPrcommVec.size(); ivec++) {
        if ((cvPrcommVec[ivec].rank != mpi_rank)) {

          MPI_Irecv(&(cvPrcommVec[ivec].nunpack_v),1,MPI_INT,cvPrcommVec[ivec].rank,
                    UPDATE_TAG1,mpi_comm,&(mrs->recvRequestArray[requestCount]));

          MPI_Issend(&(send_count[cvPrcommVec[ivec].rank]),1,MPI_INT,cvPrcommVec[ivec].rank,
                     UPDATE_TAG1,mpi_comm,&(mrs->sendRequestArray[requestCount]));

          requestCount++;
        }
      }
      
      // wait for all messages to be received (and sent)...
      if (requestCount > 0) {
        MPI_Waitall(requestCount,mrs->recvRequestArray,MPI_STATUSES_IGNORE);
        MPI_Waitall(requestCount,mrs->sendRequestArray,MPI_STATUSES_IGNORE);
      }
      
      for (int ivec = 0; ivec < cvPrcommVec.size(); ivec++) {
        assert(cvPrcommVec[ivec].nunpack_v>=0);
      }

      delete [] mrs->sendRequestArray;
      delete [] mrs->recvRequestArray;
      mrs->nsend = 0;
      mrs->nrecv = 0;

      // Now send and recv data only from relevant ranks
      const int UPDATE_TAG2 = 535353;
      const int UPDATE_TAG3 = 646464;
      int nLpRecv = 0;
      for (int ivec = 0; ivec < cvPrcommVec.size(); ivec++) {
        nLpRecv += cvPrcommVec[ivec].nunpack_v;
      }
      mrs->unpack_buf_double = new double [nLpRecv * T::data_size()];
      mrs->unpack_buf_int = new int [nLpRecv];

      // how many send and recv
      for (int ivec = 0; ivec < cvPrcommVec.size(); ivec++) {
        
        if (cvPrcommVec[ivec].nunpack_v > 0) {
          if ((cvPrcommVec[ivec].rank != mpi_rank)) {
            mrs->nrecv += 2;
          }
        }

        if (send_count[cvPrcommVec[ivec].rank] > 0) {
            mrs->nsend += 2;
        }

      }

      mrs->recvRequestArray = new MPI_Request [mrs->nrecv];
      mrs->sendRequestArray = new MPI_Request [mrs->nsend];

      int sendRequestCount = 0;
      int recvRequestCount = 0;
      int offset = 0;
      for (int ivec = 0; ivec < cvPrcommVec.size(); ivec++) {
        
        // post the Irecv...
        if (cvPrcommVec[ivec].nunpack_v > 0) {
          if ((cvPrcommVec[ivec].rank != mpi_rank)) {
            const int nbuf_int = cvPrcommVec[ivec].nunpack_v;
            const int nbuf_double = cvPrcommVec[ivec].nunpack_v * T::data_size();
            MPI_Irecv(mrs->unpack_buf_int+offset,nbuf_int,MPI_INT,cvPrcommVec[ivec].rank,
                      UPDATE_TAG2,mpi_comm,&(mrs->recvRequestArray[recvRequestCount++]));
            MPI_Irecv(mrs->unpack_buf_double+offset*T::data_size(),nbuf_double,MPI_DOUBLE,cvPrcommVec[ivec].rank,
                      UPDATE_TAG3,mpi_comm,&(mrs->recvRequestArray[recvRequestCount++]));
          }
        }
        
        // and send...      
        if (send_count[cvPrcommVec[ivec].rank] > 0) {
          const int nbuf_int = send_count[cvPrcommVec[ivec].rank];
          const int nbuf_double = send_count[cvPrcommVec[ivec].rank] * T::data_size();
          map<int, int*>::iterator iter_int = rank_send_buf_int_map.find(cvPrcommVec[ivec].rank);
          map<int, double*>::iterator iter_double = rank_send_buf_double_map.find(cvPrcommVec[ivec].rank);
          MPI_Issend(iter_int->second,nbuf_int,MPI_INT,cvPrcommVec[ivec].rank,
                     UPDATE_TAG2,mpi_comm,&(mrs->sendRequestArray[sendRequestCount++]));
          MPI_Issend(iter_double->second,nbuf_double,MPI_DOUBLE,cvPrcommVec[ivec].rank,
                     UPDATE_TAG3,mpi_comm,&(mrs->sendRequestArray[sendRequestCount++]));
        }

        offset += cvPrcommVec[ivec].nunpack_v;
      }
      assert(offset == nLpRecv);
      assert(sendRequestCount==mrs->nsend);
      assert(recvRequestCount==mrs->nrecv);
      
      if (requestCount > 0) {
        if (recvRequestCount > 0)
          MPI_Waitall(recvRequestCount,mrs->recvRequestArray,MPI_STATUSES_IGNORE);
        if (sendRequestCount > 0)
          MPI_Waitall(sendRequestCount,mrs->sendRequestArray,MPI_STATUSES_IGNORE);
      }

      lpCls->resize(np_new+nLpRecv);

      // unpack...

      for (int irecv = 0; irecv < nLpRecv; ++irecv) {
        // here you should over write the icv after calling unpack
        // the new unpack routine brings icv over too
        // in this version of relocation icv of the sent particle is not set to the destination icv
        // you can set the sending particle's icv to the destination icv and we don't need to send over the index
        lpCls->lp[np_new].unpack(mrs->unpack_buf_double+irecv*T::data_size());
        lpCls->lp[np_new].icv = mrs->unpack_buf_int[irecv]; assert((lpCls->lp[np_new].icv >= 0)&&(lpCls->lp[np_new].icv < ncv));
        ++np_new;
      }

      for (map<int,int*>::iterator iter = rank_send_buf_int_map.begin(); iter != rank_send_buf_int_map.end(); iter++)
        delete[] iter->second;
      for (map<int,double*>::iterator iter = rank_send_buf_double_map.begin(); iter != rank_send_buf_double_map.end(); iter++)
        delete[] iter->second;
      delete [] mrs->sendRequestArray; mrs->sendRequestArray = NULL;
      delete [] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
      delete [] mrs->unpack_buf_double; mrs->unpack_buf_double = NULL;
      delete [] mrs->unpack_buf_int; mrs->unpack_buf_int = NULL;
      delete mrs;

      ///////////////////////////

/*
      ///////////////////////////////////////////
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
      const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

      int * send_buf_int = new int[send_count_sum]; // for the icv's on the recv'ing rank
      double * send_buf_double = new double[send_count_sum*T::data_size()];

      // now pack the buffer and compress lpVec...

      int np_new = 0;
      //double R[9], t[3];
      for (int ip = 0; ip < lpCls->size(); ++ip) {
        if (lpCls->lp[ip].icv >= ncv) {
          int rank,bits,index;
          BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[lpCls->lp[ip].icv-ncv]);
          send_buf_int[send_disp[rank]] = index; // the icv on rank where we are to unpack
          //lpCls->lp[ip].pack(send_buf_double+send_disp[rank]*T::data_size());

          const int ipR = iRt_g[lpCls->lp[ip].icv-ncv][0];
          const int ipt = iRt_g[lpCls->lp[ip].icv-ncv][1];
          if ((ipR >= 0) && (ipt >= 0)) {
            // rotate/translate xp/xp0 and rotate up
            lpCls->lp[ip].packRt(send_buf_double+send_disp[rank]*T::data_size(),per_R[ipR],per_t[ipt]);
          }
          else if (ipR >= 0) {
            // rotate xp/xp0 and rotate up
            lpCls->lp[ip].packR(send_buf_double+send_disp[rank]*T::data_size(),per_R[ipR]);
          }
          else if (ipt >= 0) {
            // translate xp/xp0
            lpCls->lp[ip].packt(send_buf_double+send_disp[rank]*T::data_size(),per_t[ipt]);
          }
          else {
            lpCls->lp[ip].pack(send_buf_double+send_disp[rank]*T::data_size());
          }
          ++send_disp[rank];
        }
        else {
          if (np_new != ip)
            lpCls->lp[np_new].copy(lpCls->lp[ip]);
          ++np_new;
        }
      }

      // figure out how many we will be receiving and make space for them...

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      int * recv_buf_int = new int[recv_count_sum];
      MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
        recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
      delete[] send_buf_int;

      FOR_RANK {
        send_count[rank] *= T::data_size();
        send_disp[rank] *= T::data_size();
        recv_count[rank] *= T::data_size();
        recv_disp[rank] *= T::data_size();
      }

      double * recv_buf_double = new double[recv_count_sum*T::data_size()];
      MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
        recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
      delete[] send_buf_double;

      // make room for data in lpVec...

      lpCls->resize(np_new+recv_count_sum);

      // unpack...

      for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
        lpCls->lp[np_new].icv = recv_buf_int[irecv]; assert((lpCls->lp[np_new].icv >= 0)&&(lpCls->lp[np_new].icv < ncv));
        lpCls->lp[np_new].unpack(recv_buf_double+irecv*T::data_size());
        // no need to set this because we are unpacking flag too
        //lpVec[np_new].flag = KEEP;
        ++np_new;
      }

      delete[] recv_buf_int;
      delete[] recv_buf_double;
*/
      ///////////////////////////////////////////

      // check if we are done...

      MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);

    } // while loop

    // check that everyone is -ve, and reset to positive...

    for (int ip = 0; ip < lpCls->size(); ++ip) {
      assert(lpCls->lp[ip].icv < 0);
      lpCls->lp[ip].icv = -lpCls->lp[ip].icv-1;
    }

    delete[] send_count;
    delete[] send_disp;
    delete[] recv_count;
    delete[] recv_disp;

    if ((mpi_rank == 0) && (step%check_interval == 0))
      cout << "OK" << endl;

    //MPI_Pause("ALL FINISHED");
  }

  virtual void impactLpBc(LpTracerState& lp, const double* xi, const double* st_normal) {
    const double dx[3] = DIFF(lp.xp,xi);
    const double dp = DOT_PRODUCT(st_normal,dx); assert(dp > 0.0);
    const double nmag2 = DOT_PRODUCT(st_normal,st_normal);
    FOR_I3 lp.xp[i] -= 2.0*dp*st_normal[i]/nmag2;
    // and update xp0 to the xi_closest, because the remaining part of the
    // trajectory may have additional intersections...
    FOR_I3 lp.xp0[i] = xi[i];
    const double un = DOT_PRODUCT(st_normal,lp.up);
    FOR_I3 lp.up[i] -= 2.0*un*st_normal[i]/nmag2;
  }

  virtual void impactLpBc(LspState& lp, const double* xi, const double* st_normal, const double& tol) {
    const double dx[3] = DIFF(lp.xp,xi);
    const double dp = DOT_PRODUCT(st_normal,dx); assert(dp > -tol);
    const double nmag2 = DOT_PRODUCT(st_normal,st_normal);
    FOR_I3 lp.xp[i] -= 2.0*dp*st_normal[i]/nmag2;
    // and update xp0 to the xi_closest, because the remaining part of the
    // trajectory may have additional intersections...
    FOR_I3 lp.xp0[i] = xi[i];
    const double un = DOT_PRODUCT(st_normal,lp.up);
    FOR_I3 lp.up[i] -= 2.0*un*st_normal[i]/nmag2;
  }

  void irregularBounceLpBc(LspState& lp, const double* xi, const double* st_normal) {
    const double dx[3] = DIFF(lp.xp,xi);
    const double dp = DOT_PRODUCT(st_normal,dx); assert(dp > 0.0);
    const double nmag2 = DOT_PRODUCT(st_normal,st_normal);
    FOR_I3 lp.xp[i] -= 2.0*dp*st_normal[i]/nmag2;
    // and update xp0 to the xi_closest, because the remaining part of the
    // trajectory may have additional intersections...
    FOR_I3 lp.xp0[i] = xi[i];

    // perturb the normal vector first
    double u_orth[3];     // the orthagonal to the normal part of the velocity vector
    FOR_I3 u_orth[i] = lp.up[i] - DOT_PRODUCT(st_normal,lp.up)*st_normal[i]/nmag2;
    double mag_u_orth = MAG(u_orth);
    double random = (double(rand())/double(RAND_MAX) - 0.5)*2.0;
    double MAX_GAMMA = 3.0/180.0*M_PI; // maximum change in angle is 3 degres, can be set as a parameter or a measure based on the hight and horizontal lengths of the roughness.
    double tan_Gamma = tan(MAX_GAMMA * random);
    double st_normal_pert[3];     // perturbed normal vector, tan_Gamma = |dn|/|n|
    FOR_I3 st_normal_pert[i] = st_normal[i] + sqrt(nmag2)*tan_Gamma*u_orth[i]/mag_u_orth;
    const double nmag2_pert = DOT_PRODUCT(st_normal_pert,st_normal_pert);
    const double un = DOT_PRODUCT(st_normal_pert,lp.up);
    double un_old[3];
    double ut_old[3];
    FOR_I3 un_old[i] = st_normal_pert[i]*un/nmag2_pert;
    FOR_I3 ut_old[i] = lp.up[i] - un_old[i];
    const double mag_un_old = MAG(un_old);
    const double mag_ut_old = MAG(ut_old);
    double un_new[3];
    double ut_new[3];
    FOR_I3 un_new[i] = -un_old[i]*resCoef;
    double slide = mag_ut_old - 7./2.*mu_s*(1.+resCoef)*mag_un_old; // |u_p| < 7/2 mu_s (1+e) |v_p| -> no sliding
    if (slide < 0) { // no sliding
      //cout << "impact_mode: "<< "NO slide" << " , 7./2.*mu_s*(1.+resCoef): " << 7./2.*mu_s*(1.+resCoef) << " , u/v: " << mag_ut_old/mag_un_old << " , tan_Gamma: " << tan_Gamma << endl;
      FOR_I3 ut_new[i] = 5./7.*ut_old[i];
    } else { // sliding
      //cout << "impact_mode: "<< "slide" << " , 7./2.*mu_s*(1.+resCoef): " << 7./2.*mu_s*(1.+resCoef) << " , u/v: " << mag_ut_old/mag_un_old << " , tan_Gamma: " << tan_Gamma << endl;
      FOR_I3 ut_new[i] = ut_old[i]*(1. - mu_k*(1+resCoef)*mag_un_old/mag_ut_old);
    }
    FOR_I3 lp.up[i] = ut_new[i] + un_new[i];
  }

  //bool sortByMass(const LspState &lhs, const LspState &rhs) {return (lhs.mp<rhs.mp);}
  //struct sortOperation{
  //  bool operator()(const LspState &lhs, const LspState &rhs) {return (lhs.mp<rhs.mp);}
  //} sortByMass;

  double getParticleMdot() {

    double mdot_total = 0.0;

    FOR_PARAM_MATCHING("LSP.INJECTOR") {

      double mdot = 0.0;
      if ((mpi_rank == 0) && (step%check_interval == 0)) cout << " > LSP.INJECTOR: " << param->getString() << " mdot: ";

      int i = 1;
      while (i < param->size()) {
        string token = param->getString(i++);
        if (token == "TYPE") {
          // ==========================================
          // TYPE of injector:
          // ==========================================
          token = param->getString(i++);
          if (token == "COUNT") {
            mdot = 0;
          }
          else if ( token == "MDOT" ) {
            mdot =  param->getDouble(i++);
          }
          else {
            if (mpi_rank == 0)
              cerr << "\n\n\n*********************************************************\n" <<
                "Error: unrecognized TYPE in LSP.INJECTOR: " << token << ". Valid types are: \n" <<
                "TYPE COUNT <n>   -- injects n particles once\n" <<
                "TYPE MDOT <mdot> -- injects at rate mdot\n" <<
                "*********************************************************\n\n\n" << endl;
            throw(0);
          }
        }
      }

      mdot_total += mdot;
      if ((mpi_rank == 0) && (step%check_interval == 0)) cout << mdot << endl;
    }
    if ((mpi_rank == 0) && (step%check_interval == 0)) cout << " > total mdot: " << mdot_total << endl;

    return mdot_total;
  }

  void collision_stochastic() {
    if (step%check_interval == 0) {
      COUT1(" > collision stochastic");
    }
    assert(SOL_LP);
    vector<double> St_monitor;
    vector<double> Pcol_monitor;
    vector<double> duMag_monitor;
    bool any_collision = false;
    int my_ncol = 0;
    const double dist_std = 1.0;

    FOR_ICV {
      // gas properties
      //const double rhog = rho[icv];
      const double mug  = mu_lam[icv];
      //const double cpg  = R_gas*gamma/(gamma-1.0);
      //const double kg   = cv_compact[icv].loc_total * cpg; // loc: lambda(k) over cp
      //const double ug[3]= {u[icv][0], u[icv][1], u[icv][2]};
      //const double Tg   = T[icv];

      const int np = paocv_i[icv+1] - paocv_i[icv];
      double up_avg[3] = {0,0,0};
      double up_rms[3] = {0,0,0};
      // only if there are at least 2 particles in a cell
      if (np <= 1) continue;

      // compute the mean and rms particle velocities
      vector<int> ipVec;
      for (int ip = paocv_i[icv]; ip != paocv_i[icv+1]; ip++) {
        FOR_I3 up_avg[i] += lpCls->lp[ip].up[i];
        FOR_I3 up_rms[i] += lpCls->lp[ip].up[i]*lpCls->lp[ip].up[i];
        ipVec.push_back(ip);
      }
      assert(ipVec.size()>1);

      FOR_I3 up_avg[i] /= np;
      FOR_I3 {
        double up_rms2 = max(0.0,(up_rms[i]/np - up_avg[i]*up_avg[i]));
        assert(up_rms2>=0);
        up_rms[i] = sqrt( up_rms2 );
      }

      // compute tke and dissipation
      double SijSij = 0;
      FOR_I3 {
        SijSij += dudx[icv][i][i]*dudx[icv][i][i];
      }
      SijSij += (dudx[icv][0][1] + dudx[icv][1][0])*(dudx[icv][0][1] + dudx[icv][1][0])/2.;
      SijSij += (dudx[icv][0][2] + dudx[icv][2][0])*(dudx[icv][0][2] + dudx[icv][2][0])/2.;
      SijSij += (dudx[icv][1][2] + dudx[icv][2][1])*(dudx[icv][1][2] + dudx[icv][2][1])/2.;

      const double nu_cv = cv_compact[icv].mu_total/rho[icv];
      const double epsilon = 2.0*nu_cv*SijSij;
      assert(epsilon>0);

      // using sgs tke, Yoshizawa 1985
      const double Cv = 0.05;
      const double delta_cv = r_vv[icv];
      const double tke = pow((nu_cv/(Cv*delta_cv)),2);
      const double flow_ts = 0.16*tke/epsilon;
      assert(flow_ts>0);

      for (int ip = paocv_i[icv]; ip != paocv_i[icv+1]; ip++) {
        assert((ip>=0)&&(ip<lpCls->size()));
        //  compute Stokes
        const double rhop = dust->calcRho(lpCls->lp[ip].Tp);
        const double Dp   = lpCls->lp[ip].dp;
        assert(Dp > 0.0);
        const double tau = rhop*Dp*Dp/(18.*mug);
        assert(tau>0);
        const double St = tau / flow_ts;
        //const double St = 1.0;
        St_monitor.push_back(St);
        const double R = exp(-0.55*pow(St,0.4));
        assert(R<1.0);

        // fictitious particle
        const double up[3] = DIFF(lpCls->lp[ip].up,up_avg);
        double up_f[3];
        FOR_I3 {
          // normal distribution variable
          double xsi = MiscUtils::randn() * dist_std;
          up_f[i] = R*up[i] + up_rms[i]*sqrt(1.0 - R*R)*xsi;
        }

        // sample diameter of the fictitious particle
        // We choose a diameter from the existing particles in the cell
        // This way the probability of the choice is according to the histogram
        const double almost_one = 0.999999;
        const int randInd = (int)floor(double(rand())/double(RAND_MAX)*almost_one*double(ipVec.size()));
        assert((randInd>=0)&&(randInd<ipVec.size()));
        const double Dp_f = lpCls->lp[ipVec[randInd]].dp;

        // collision frequency
        const double du[3] = DIFF(up,up_f);
        const double du_mag = MAG(du);
        const double np_vol = double(np)/vol_cv[icv];
        const double P_col = M_PI/4.0*(Dp+Dp_f)*(Dp+Dp_f)*du_mag*np_vol*dt;
        Pcol_monitor.push_back(P_col);
        duMag_monitor.push_back(du_mag);
        //cout << "rank: " << mpi_rank << " P_col: " << P_col << " du_mag: " << du_mag << endl;
        //assert(P_col <= 1);

        double randVal = double(rand())/double(RAND_MAX);
        if (randVal < P_col) {
          // collision happens
          ++my_ncol;
          any_collision = true;
          const double mp = lpCls->lp[ip].mp;
          const double dp_ratio = Dp_f/Dp;
          const double mp_f = mp*dp_ratio*dp_ratio*dp_ratio;
          //cout << "rank: " << mpi_rank << " Dp_f: " << Dp_f << endl;

          double u_rel[3] = DIFF(up, up_f); // the original coordinate system
          const double u_rel_mag = MAG(u_rel);
          //if ((u_rel_mag/MAG(up))<1e-4) continue; // if the velocity difference is small ignore collision
          assert(u_rel_mag>0);
          double e1[3], e2[3], e3[3];
          FOR_I3 e1[i] = u_rel[i]/u_rel_mag;
          MiscUtils::getBestE1E2FromE0(e2, e3, e1);
          //if ((abs(e1[0]) >= abs(e1[1])) && (abs(e1[0]) >= abs(e1[2]))) {
          //  e2[0] = 0.0;
          //  e2[1] = e1[2];
          //  e2[2] = -e1[1];
          //  NORMALIZE(e2);
          //  const double temp[3] = CROSS_PRODUCT(e1, e2);
          //  FOR_I3 e3[i] = temp[i];
          //} else if ((abs(e1[1]) >= abs(e1[0])) && (abs(e1[1]) >= abs(e1[2]))){
          //  e2[0] = e1[2];
          //  e2[1] = 0.0;
          //  e2[2] = -e1[0];
          //  NORMALIZE(e2);
          //  const double temp[3] = CROSS_PRODUCT(e1, e2);
          //  FOR_I3 e3[i] = temp[i];
          //} else {
          //  e2[0] = e1[1];
          //  e2[1] = -e1[0];
          //  e2[2] = 0.0;
          //  NORMALIZE(e2);
          //  const double temp[3] = CROSS_PRODUCT(e1, e2);
          //  FOR_I3 e3[i] = temp[i];
          //}
          //assert((MAG(e1)-1.0<1e-12));
          //assert((MAG(e2)-1.0<1e-12));
          //assert((MAG(e3)-1.0<1e-12));
          //assert(DOT_PRODUCT(e1, e2) < 1e-12);
          //assert(DOT_PRODUCT(e1, e3) < 1e-12);
          //assert(DOT_PRODUCT(e2, e3) < 1e-12);

          // random variables with uniform distribution
          const double beta = double(rand())/double(RAND_MAX);
          const double alpha = double(rand())/double(RAND_MAX)*2.0*M_PI;   // 0 <= alpha <= 2pi

          const double beta_sqrt = sqrt(beta);
          const double L_pr = beta_sqrt*Dp/2.0;
          //const double theta = asin(beta_sqrt); // 0 <= theta <= pi/2
          const double sin_theta = beta_sqrt;
          const double cos_theta = sqrt(1.0 - sin_theta*sin_theta);

          const double X0[3] = {Dp/2.*cos_theta, L_pr*sin(alpha), L_pr*cos(alpha)};
          const double X0_mag = MAG(X0);
          assert(X0_mag>0);
          double ee1[3], ee2[3], ee3[3];
          FOR_I3 ee1[i] = X0[i]/X0_mag;
          MiscUtils::getBestE1E2FromE0(ee2, ee3, ee1);
          //if ((abs(ee1[0]) >= abs(ee1[1])) && (abs(ee1[0]) >= abs(ee1[2]))) {
          //  ee2[0] = 0.0;
          //  ee2[1] = ee1[2];
          //  ee2[2] = -ee1[1];
          //  NORMALIZE(ee2);
          //  const double temp[3] = CROSS_PRODUCT(ee1, ee2);
          //  FOR_I3 ee3[i] = temp[i];
          //} else if ((abs(ee1[1]) >= abs(ee1[0])) && (abs(ee1[1]) >= abs(ee1[2]))){
          //  ee2[0] = ee1[2];
          //  ee2[1] = 0.0;
          //  ee2[2] = -ee1[0];
          //  NORMALIZE(ee2);
          //  const double temp[3] = CROSS_PRODUCT(ee1, ee2);
          //  FOR_I3 ee3[i] = temp[i];
          //} else {
          //  ee2[0] = ee1[1];
          //  ee2[1] = -ee1[0];
          //  ee2[2] = 0.0;
          //  NORMALIZE(ee2);
          //  const double temp[3] = CROSS_PRODUCT(ee1, ee2);
          //  FOR_I3 ee3[i] = temp[i];
          //}
          //assert((MAG(ee1)-1.0<1e-12));
          //assert((MAG(ee2)-1.0<1e-12));
          //assert((MAG(ee3)-1.0<1e-12));
          //assert(DOT_PRODUCT(ee1, ee2) < 1e-12);
          //assert(DOT_PRODUCT(ee1, ee3) < 1e-12);
          //assert(DOT_PRODUCT(ee2, ee3) < 1e-12);

          double u_rel_p[3] = {u_rel_mag, 0.0, 0.0}; // u_rel in the prime coordinate system
          double u_rel_pp[3] = {DOT_PRODUCT(u_rel_p, ee1), DOT_PRODUCT(u_rel_p, ee2), DOT_PRODUCT(u_rel_p, ee3)}; // u_rel in the double prime coordinate system

          // compute post collision velocity in double prime coordinate system -  only the first component should be modified
          u_rel_pp[0] *= (1.0 - (1.0+resCoef_colMod)/(1.0+mp/mp_f));

          // transfer relative velocity to prime coordinate system
          FOR_I3 u_rel_p[i] = u_rel_pp[0]*ee1[i] + u_rel_pp[1]*ee2[i] + u_rel_pp[2]*ee3[i];

          // transfer relative velocity to the main coordinate system
          FOR_I3 u_rel[i] = u_rel_p[0]*e1[i] + u_rel_p[1]*e2[i] + u_rel_p[2]*e3[i];

          // updating the particle velocity
          FOR_I3 lpCls->lp[ip].up[i] = u_rel[i] + up_f[i] + up_avg[i];

        } // if collision
      } // poc
    } // icv
    int any_collision_i = int(any_collision);
    int any_collision_global;
    MPI_Allreduce(&any_collision_i, &any_collision_global, 1, MPI_INT, MPI_LOR, mpi_comm);
    if ((any_collision_global)&&(step%check_interval==0)) {
      double * r1 = new double[St_monitor.size()];
      for (int i = 0; i < St_monitor.size(); i++) {
        r1[i] = St_monitor[i];
      }
      dumpRange(r1, St_monitor.size(), "Stokes number");
      delete [] r1;

      r1 = new double[Pcol_monitor.size()];
      for (int i = 0; i < Pcol_monitor.size(); i++) {
        r1[i] = Pcol_monitor[i];
      }
      dumpRange(r1, Pcol_monitor.size(), "P_col");
      delete [] r1;

      r1 = new double[duMag_monitor.size()];
      for (int i = 0; i < duMag_monitor.size(); i++) {
        r1[i] = duMag_monitor[i];
      }
      dumpRange(r1, duMag_monitor.size(), "du_mag");
      delete [] r1;

      int ncol;
      MPI_Reduce(&my_ncol,&ncol,1,MPI_INT,MPI_SUM,0,mpi_comm);
      IF_RANK0 cout << " > time: " << time << " ,number-of-collisions: " << ncol << " ,dist_std: " << dist_std << endl;
    }

  }

  void collision_deterministic() {
    if (step%check_interval == 0) {
      COUT1(" > collision deterministic");
    }
    assert(SOL_LP);
    //vector<double> St_monitor;
    //vector<double> Pcol_monitor;
    //vector<double> duMag_monitor;
    bool any_collision = false;
    int my_ncol = 0;

    FOR_ICV {

      const int np = paocv_i[icv+1] - paocv_i[icv];
      // only if there are at least 2 particles in a cell
      if (np <= 1) continue;

      for (int ip0 = paocv_i[icv]; ip0 != paocv_i[icv+1]; ip0++) {
        double closest_t = -1;
        int closest_ip = -1;
        const double * const u0 = lpCls->lp[ip0].up;
        const double * const x0 = lpCls->lp[ip0].xp;
        const double r0 = lpCls->lp[ip0].dp/2.;
        for (int ip1 = paocv_i[icv]; ip1 != paocv_i[icv+1]; ip1++) {
          if (ip1==ip0) continue; // skip self collision

          //cout << "working on ip0: " << ip0 << " ip1: " << ip1 << endl;

          // check if particle ip0 collides with any other particle during the next time step
          const double * const u1 = lpCls->lp[ip1].up;
          const double * const x1 = lpCls->lp[ip1].xp;
          const double r1 = lpCls->lp[ip1].dp/2.;

          const double u_rel[3] = DIFF(u0,u1);
          double xp0[3];
          FOR_I3 xp0[i] = x0[i]+dt*u_rel[i];

          double dx1x0[3] = DIFF(x1,x0);
          double dx1xp0[3] = DIFF(x1,xp0);
          double dxp0x0[3] = DIFF(xp0,x0);

          double dist_min = 0; // minimum distance between x1 and the x0xp0 segment
          const double t1 = DOT_PRODUCT(dx1x0,dxp0x0);
          const double t2 = DOT_PRODUCT(dxp0x0,dxp0x0);

          //cout << "t1: " << t1 << " t2: " << t2 << endl;

          if (t1 <= 0) { // no impact possible
            //cout << "continue " << endl;
            continue;
          } else if (t1 < t2) { // the projection of x1x0 on xp0x0 is within xp0x0
            double xt[3];
            FOR_I3 xt[i] = x0[i] + t1/t2*dxp0x0[i];
            double dxt[3] = DIFF(x1,xt);
            dist_min = MAG(dxt);
            //cout << "d_min inside " << endl;
          }
          else {
            double mag_1 = MAG(dx1x0);
            double mag_2 = MAG(dx1xp0);
            dist_min = min(mag_1, mag_2);
            //cout << "d_min outside " << endl;
          }
          const double t = t1/t2;
          assert(DOT_PRODUCT(dx1x0,dxp0x0)>0); // particles are moving twards each other
          assert(dist_min > 0);
          assert((t>0));

          if (dist_min < (r0+r1)) { // collision with ip1, should check if this is the closest collision
            if ((closest_ip==-1)||(t<closest_t)) {
              closest_t = t;
              closest_ip = ip1;
            }
          }// else {
            //cout << "no collision, dist_min > r0+r1" << endl;
          //}
        }
        if (closest_ip != -1) { // collision for ip0 is found
          //cout << "found collision for ip0: " << ip0 << " closest_ip: " << closest_ip << endl;
          any_collision = true;
          my_ncol++;
          assert((closest_ip>=0)&&(closest_ip<lpCls->size()));
          assert(lpCls->lp[closest_ip].icv == icv);

          const double * const u1 = lpCls->lp[closest_ip].up;
          const double * const x1 = lpCls->lp[closest_ip].xp;
          const double r1 = lpCls->lp[closest_ip].dp/2.;

          const double mp0 = lpCls->lp[ip0].mp;
          const double mp1 = lpCls->lp[closest_ip].mp;

          const double u_rel[3] = DIFF(u0,u1);
          double xp0[3];
          FOR_I3 xp0[i] = x0[i]+dt*u_rel[i];

          //cout << "u_rel: " << COUT_VEC(u_rel) << endl;
          const double dx1x0[3] = DIFF(x1,x0);
          //const double dx1xp0[3] = DIFF(x1,xp0);
          const double dxp0x0[3] = DIFF(xp0,x0);

          // solve an equation of form: a.t^2 + 2b.t + c = 0 to find the intersection of x0xp0 with a sphere with radius (r0+r1) and center x1
          const double a = DOT_PRODUCT(dxp0x0,dxp0x0);
          const double b = -DOT_PRODUCT(dx1x0,dxp0x0);
          const double c = DOT_PRODUCT(dx1x0,dx1x0) - (r0+r1)*(r0+r1);
          assert(b<=0);

          // solve the 2nd order equation
          const double delta = b*b - a*c;
          if (delta < 0) continue; // we should have at least 1 solution, but if for numerical reasons it is a very small negative number skip...
          assert(delta >= 0);

          const double t1 = (-b+sqrt(delta))/a;
          const double t2 = (-b-sqrt(delta))/a;
          const double t_min = min(t1,t2);

          //cout << "t1: " << t1 << " t2: " << t2 << " t_min: " << t_min << endl;
          double xc[3];
          FOR_I3 xc[i] = x0[i] + t_min*dxp0x0[i];

          const double dx1xc[3] = DIFF(x1,xc);
          const double mag2 = DOT_PRODUCT(dx1xc,dx1xc);
          //cout << "sqrt(mag2): " << sqrt(mag2) << " r0+r1: " << r0+r1 << endl;
          assert((mag2-(r0+r1)*(r0+r1))<1e-10);

          double u_rel_n[3];
          double proj = DOT_PRODUCT(u_rel,dx1xc);
          if (proj<0) continue; // at the impact point they should come twards each other
          assert(proj>=0);
          FOR_I3 u_rel_n[i] = proj*dx1xc[i]/mag2;

          //const double u_rel_t[3] = DIFF(u_rel,u_rel_n);

          // update the velocity
          FOR_I3 lpCls->lp[ip0].up[i] -= (1.0+resCoef_colMod)/(1.0+mp0/mp1)*u_rel_n[i];
          FOR_I3 lpCls->lp[closest_ip].up[i] += (1.0+resCoef_colMod)/(1.0+mp1/mp0)*u_rel_n[i];
          //cout << "  after up0: " << COUT_VEC(lpCls->lp[ip0].up) << endl;
          //cout << "  after up_closest: " << COUT_VEC(lpCls->lp[closest_ip].up) << endl;
        }
      }

    }

    // monitoring collision
    int any_collision_global;
    int any_collision_i = int(any_collision);
    MPI_Allreduce(&any_collision_i, &any_collision_global, 1, MPI_INT, MPI_LOR, mpi_comm);
    if ((any_collision_global)&&(step%check_interval==0)) {
      //double * r1 = new double[St_monitor.size()];
      //for (int i = 0; i < St_monitor.size(); i++) {
      //  r1[i] = St_monitor[i];
      //}
      //dumpRange(r1, St_monitor.size(), "Stokes number");
      //delete [] r1;

      //r1 = new double[Pcol_monitor.size()];
      //for (int i = 0; i < Pcol_monitor.size(); i++) {
      //  r1[i] = Pcol_monitor[i];
      //}
      //dumpRange(r1, Pcol_monitor.size(), "P_col");
      //delete [] r1;

      //r1 = new double[duMag_monitor.size()];
      //for (int i = 0; i < duMag_monitor.size(); i++) {
      //  r1[i] = duMag_monitor[i];
      //}
      //dumpRange(r1, duMag_monitor.size(), "du_mag");
      //delete [] r1;

      int ncol;
      MPI_Reduce(&my_ncol,&ncol,1,MPI_INT,MPI_SUM,0,mpi_comm);
      IF_RANK0 cout << " > time: " << time << " ,number-of-collisions: " << ncol << endl;
    }
  }

  void computeParticleConcentration() {
    const double CPar_n = double(getLspSize()) / vol_total;
    if (CPar_n > 0) {
      for (int icv = 0; icv < ncv; icv++) {
        const int np = paocv_i[icv+1]-paocv_i[icv];
        CPar[icv] = (CPar[icv]*time_stat*CPar_n_old + double(np)/vol_cv[icv]*dt)/(time_stat+dt)/CPar_n;
        CPar2[icv] = (double(np)/vol_cv[icv])/CPar_n;
      }
    }
    CPar_n_old = CPar_n;
    time_stat += dt;
  }

  void computeNearWallConcentration() {
    if (LpStatZoneVec.size() > 0) {
      double * delta_of_zone = new double [LpStatZoneVec.size()];
      for (int izone = 0; izone < LpStatZoneVec.size(); izone++) {
        BfZone * zone_ptr = LpStatZoneVec[izone];
        double buf[4];
        FOR_I4 buf[i] = 0.0;
        for (int ibf = 0; ibf < zone_ptr->nbf; ibf++) {
          const int icv = zone_ptr->cvobf[ibf];
          const double visc_coeff = (mu_lam[icv]+mu_sgs[icv])*zone_ptr->area_over_delta_bf[ibf];
          // assuming the boundary zone has zero velocity
          const double du_mag = MAG(u[icv]);
          const double tau_wall = visc_coeff*(du_mag)/area_bf[ibf];
          const double nu = (mu_lam[icv]+mu_sgs[icv])/rho[icv];

          buf[0] += area_bf[ibf];
          buf[1] += area_bf[ibf]*tau_wall;
          buf[2] += area_bf[ibf]*nu;
          buf[3] += area_bf[ibf]*rho[icv];
        }
        double buf_sum[4];
        MPI_Allreduce(buf,buf_sum,4,MPI_DOUBLE,MPI_SUM,mpi_comm);
        const double rho_avg = buf_sum[3]/buf_sum[0];
        const double nu_avg = buf_sum[2]/buf_sum[0];
        const double tau_wall_avg = buf_sum[1]/buf_sum[0];
        const double u_tau_avg = sqrt(tau_wall_avg/rho_avg);
        const double delta_avg = nu_avg/u_tau_avg;
        delta_of_zone[izone] = delta_avg;
        COUT1(" > zone "<<zone_ptr->getName()<<", time, avg tau_wall, avg u_tau, avg delta: "<<
              time<<" "<<tau_wall_avg<<" "<<u_tau_avg<<" "<<delta_avg);
      }

      double * my_NPar_zone = new double [LpStatZoneVec.size()];
      double * my_vol_zone = new double [LpStatZoneVec.size()];
      for (int izone = 0; izone < LpStatZoneVec.size(); izone++) {
        my_NPar_zone[izone] = 0.0;
        my_vol_zone[izone] = 0.0;
      }
      int * flag_cv = new int [ncv];
      FOR_ICV flag_cv[icv] = 0;
      //double my_NPar_corners = 0.0;
      //double my_vol_corners = 0.0;
      double my_NPar_nearWalls = 0.0;
      double my_vol_nearWalls = 0.0;

      const double y_plus_tol = 15.0;
      FOR_ICV {
        for (int izone = 0; izone < LpStatZoneVec.size(); izone++) {
          BfZone * zone_ptr = LpStatZoneVec[izone];
          double dx[3];
          FOR_I3 dx[i] = x_cv[icv][i] - zone_ptr->x_global[i];
          double normal_mag = MAG(zone_ptr->n_global);
          double normal[3];
          FOR_I3 normal[i] = zone_ptr->n_global[i]/normal_mag;
          double dist = abs(DOT_PRODUCT(dx,normal));

          if (dist < delta_of_zone[izone]*y_plus_tol) {
            const int np = paocv_i[icv+1]-paocv_i[icv];
            my_NPar_zone[izone] += np;
            my_vol_zone[izone] += vol_cv[icv];
            //cv_wall_flag[icv] = 1.0;
            if (flag_cv[icv]==0) {
              my_NPar_nearWalls += np;
              my_vol_nearWalls += vol_cv[icv];
              flag_cv[icv] = 1;
            }
          }
        }
      }
      double * NPar_zone = new double [LpStatZoneVec.size()];
      double * vol_zone = new double [LpStatZoneVec.size()];
      double NPar_nearWalls;
      double vol_nearWalls;
      MPI_Reduce(my_NPar_zone,NPar_zone,LpStatZoneVec.size(),MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      MPI_Reduce(my_vol_zone,vol_zone,LpStatZoneVec.size(),MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      MPI_Reduce(&my_NPar_nearWalls,&NPar_nearWalls,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      MPI_Reduce(&my_vol_nearWalls,&vol_nearWalls,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

      const double NP_over_vol_total = getLspSize()/vol_total;
      const int np_total = getLspSize();
      IF_RANK0 {
        for (int izone = 0; izone < LpStatZoneVec.size(); izone++) {
          BfZone * zone_ptr = LpStatZoneVec[izone];
          const double NP_over_vol = NPar_zone[izone]/vol_zone[izone];
          cout << " > zone " << zone_ptr->getName() << ", time, number of particles, concentration: " << time << " "
              << NPar_zone[izone] << " " << NP_over_vol/NP_over_vol_total << endl;
        }
        const double NP_over_vol = NPar_nearWalls/vol_nearWalls;
        cout << " > Near-Walls time, number of particles, Np over Ntot, concentration: " << time << " "
            << NPar_nearWalls << " " << NPar_nearWalls/double(np_total) << " "
            << NP_over_vol/NP_over_vol_total << endl;
      }

      delete [] delta_of_zone;
      delete [] my_NPar_zone;
      delete [] my_vol_zone;
      delete [] NPar_zone;
      delete [] vol_zone;
      delete [] flag_cv;

    }

  }

  void dumpHistMass(const double * dp,const int n,const double rho,const int nbin,const string& message, vector<double>& x, vector<double>& y, bool verbose) {

    assert(mpi_rank == 0);

    // determine the range...

    assert(n > 0);
    double dp_min = dp[0];
    double dp_max = dp[0];
    for (int i = 0; i < n; ++i) {
      assert( dp[i] == dp[i] ); // nan check
      dp_min = min(dp_min,dp[i]);
      dp_max = max(dp_max,dp[i]);
    }

    // expand very slightly...

    //double eps = 1.0E-6*(dp_max-dp_min);
    double eps = 1e-10;
    dp_min -= eps;
    dp_max += eps;

    // build historgram...

    //const int nbin = 120;
    double count[nbin];
    for (int ib = 0; ib < nbin; ++ib)
      count[ib] = 0;

    for (int i = 0; i < n; ++i) {
      int ib = (int)((double)nbin*(dp[i]-dp_min)/(dp_max-dp_min));
      assert((ib >= 0)&&(ib < nbin));
      count[ib] += rho * M_PI/6.*pow(dp[i],3);
    }

    assert(x.size()==0);
    assert(y.size()==0);
    double sum_count = 0;
    for (int ib = 0; ib < nbin; ++ib) {
      x.push_back(dp_min + double(ib)*(dp_max-dp_min)/double(nbin));
      sum_count += count[ib];
    }
    double dx = x[1] - x[0];
    for (int ib = 0; ib < nbin; ++ib) {
      y.push_back(double(count[ib])/double(sum_count)/dx);
    }

    if (not verbose) return;

    double count_max = 0;
    for (int ib = 0; ib < nbin; ++ib)
      count_max = max(count_max,count[ib]);

    // and print it...

    cout << endl;
    for (int ib = 0; ib < nbin; ++ib)
      cout << "-";
    cout << endl;
    cout << "Historgram: " << message << ", range: " << dp_min << " " << dp_max << ", samples: " << n << ", nbins: " << nbin << endl;
    int nrows = 40;
    for (int ir = nrows-1; ir >= 0; --ir) {
      for (int ib = 0; ib < nbin; ++ib) {
	      if (count[ib]*nrows >= count_max*ir)
	        cout << "*";
	      else
	        cout << " ";
      }
      cout << endl;
    }
    for (int ib = 0; ib < nbin; ++ib)
      cout << "-";
    cout << endl;

  }

  void computeForces() {

    assert(0);

  }

  void queryLpBcs_interval() {
    for (int iv = 0; iv < queryBcVec.size(); iv++) {
      if (step%queryBcVec[iv].second == 0) {
        queryBcVec[iv].first->queryBc();
      }
      if ((step % lsp_stats_interval == 0)  && (lsp_stats_interval > 0))
        queryBcVec[iv].first->updateStats();
    }
  }

  void queryLpBcs() {
    for (int iv = 0; iv < queryBcVec.size(); iv++) {
      queryBcVec[iv].first->queryBc();
    }
  }

  void checkLpContain() {

    //MPI_Barrier(mpi_comm);

    int my_count_in = 0;
    int my_count_out = 0;
    for (int ip = 0; ip < lpCls->size(); ip++) {
      const double * const xp = lpCls->lp[ip].xp;
      
      int icv_closest;
      const int isInside = pointIsInside(xp,icv_closest);

      if (isInside) {
        my_count_in++;
      } else {
        my_count_out++;
        //cout << "rank: " << mpi_rank << " xp: " << COUT_VEC(xp) << " ip: " << ip << " icv: " << lpCls->lp[ip].icv << " ncv: " << ncv << endl;
      }
    }

    int count_in, count_out;
    MPI_Reduce(&my_count_in, &count_in, 1, MPI_INT, MPI_SUM, 0, mpi_comm);
    MPI_Reduce(&my_count_out, &count_out, 1, MPI_INT, MPI_SUM, 0, mpi_comm);

    int np_total = getLspSize();
    if ((mpi_rank == 0)&&(step%check_interval==0))
      cout << " > out of " << np_total << " particles " << count_in << " are in, " << count_out << " are out... " << endl;
  }

  int pointIsInside(const double xp[3], int &icv_ret) {
  
    icv_ret = -1;
  
    // returns:
    // 1: is inside   icv_ret: icv_closest
    // 0: is outside  icv_ret: -1
    StaticSolver::ensureCvAdt();
  
    vector<int> cvList;
    cvAdt->buildListForPoint(cvList, xp);
  
    int icv_closest = -1;
    double d2_closest;
    for (int ivec = 0; ivec < cvList.size(); ivec++) {
      const int icv = cvList[ivec];
      const double d2 = DIST2(x_vv[icv], xp);
      if ((icv_closest==-1)||(d2<d2_closest)) {
        d2_closest = d2;
        icv_closest = icv;
      }
    }
    if (icv_closest==-1) {
      // could not find a cell to own the point
      return 0;
    }
  
    // find the closest cell to the point
    assert((icv_closest>=0)&&(icv_closest<ncv));
    int done = 0;
    int counter = 0;
    int counter_max = 1000;
    while (!done) {
      
      done = 1;
      int icv_nbr_closest = -1;
      double d2_nbr_closest;
      for (int coc = cvocv_i[icv_closest]+1; coc != cvocv_i[icv_closest+1]; coc++) {
        const int icv_nbr = cvocv_v[coc];
        assert(icv_nbr!=icv_closest);
        const double d2_nbr = DIST2(x_vv[icv_nbr], xp);
        if ((icv_nbr_closest==-1)||(d2_nbr < d2_nbr_closest)) {
          d2_nbr_closest = d2_nbr;
          icv_nbr_closest = icv_nbr;
        }
      }
  
      if (d2_nbr_closest < d2_closest) {
        // we got a closer neighbor
        done = 0;
        if (icv_nbr_closest >= ncv) {
          // the neighbor is a ghost so we are outside
          return 0;
        } else {
          icv_closest = icv_nbr_closest;
          d2_closest = d2_nbr_closest;
        }
      }
  
      ++counter;
      if (counter>counter_max) cout << "WARNING, reached maximum iteration for rank: " << mpi_rank << " counter: " << counter << " d2_closest: " << d2_closest << " icv_closest: " << icv_closest << " xp: " << COUT_VEC(xp) << endl;
  
    }
  
    // now we should have icv_closest with the point inside of it
    assert((icv_closest>=0)&&(icv_closest<ncv));
  
    // check if the closest cell has a boundary surface and decide if it is inside
    const int my_nbf = bfocv_i[icv_closest+1] - bfocv_i[icv_closest];
    if (my_nbf==0) {
      // no boundary faces so the point is inside
      icv_ret = icv_closest;
      return 1;
    }
    else {
      
      // and we need the closest ibf,ist_ss...
      // We need a length scale...
      const double tol2 = 1.0E-10*pow(vol_cv[icv_closest],2.0/3.0);
      int ibf_closest = -1;
      int ist_ss_closest;
      double d2_closest,dist_closest;
      //double dx_closest[3],normal_closest[3];
      for (int boc = bfocv_i[icv_closest]; boc != bfocv_i[icv_closest+1]; ++boc) {
        const int ibf = bfocv_v[boc];
        for (int sob = sstobf_i[ibf]; sob != sstobf_i[ibf+1]; ++sob) {
          const int ist_ss = sstobf_v[sob];
                
          // the corner coords of this sub-surface tri are (recall that the
          // subsurface tri already has any periodic transfrom taken care of)...
          const double * const x0 = subSurface->xp[subSurface->spost[ist_ss][0]];
          const double * const x1 = subSurface->xp[subSurface->spost[ist_ss][1]];
          const double * const x2 = subSurface->xp[subSurface->spost[ist_ss][2]];
          
          double xp_tri[3]; MiscUtils::getClosestPointOnTriRobust(xp_tri,xp,x0,x1,x2);
          
          const double dx[3] = DIFF(xp_tri, xp);
          const double d2 = DOT_PRODUCT(dx,dx);
          
          if ((ibf_closest == -1)||(d2 < d2_closest+tol2)) {
            const double normal[3] = TRI_NORMAL_2(x0,x1,x2);
            const double normal_mag = MAG(normal);
            assert(normal_mag != 0); // the subsurface should not have zero-area tris!
            const double dist = DOT_PRODUCT(dx,normal)/normal_mag;
            if ((ibf_closest == -1)||(d2 < d2_closest-tol2)||(fabs(dist) > fabs(dist_closest))) {
              ibf_closest = ibf;
              ist_ss_closest = ist_ss;
              d2_closest = d2;
              //FOR_I3 dx_closest[i] = dx[i];
              //FOR_I3 normal_closest[i] = normal[i];
              dist_closest = dist;
            }
          }
        }
      }
      assert(ist_ss_closest >= 0);
      if (dist_closest >= 0.0) {
        // dist with positive sign means we are inside the fluid volume...
        icv_ret = icv_closest;
        return 1;
      } else {
        // dist with negative sign means we are outside the fluid volume...
        return 0;
      }
    }
  
    // we should never reach here
    assert(0);
    return 0;
  }

};

#undef FOR_BCZONE
#endif
