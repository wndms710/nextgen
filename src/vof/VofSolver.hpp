#ifndef VOFSOLVER_HPP
#define VOFSOLVER_HPP

#include "FlowSolver.hpp"
#include "VofSolverBcs.hpp"
#include "DropStats.hpp"
#include "TrilinosInterface.hpp"
#include "AverageOperator.hpp"
#include "LspInjector.hpp"
#include "Lsp.hpp"
#include "DropStats.hpp"
#include "Adt.hpp"
#include "CtiLiquid.hpp"
#include "CtiSolid.hpp"
#include "LspBcs.hpp"
#include "Multigrid.hpp"

enum Solvers {
  AMG_SOLVER,
  BCGSTAB_SOLVER,
  TRILINOS_SOLVER,
  JACOBI_SOLVER
};


enum GeomWindowType {
  FAZONE,
  SPHERE,
  TCONE,
  ELLIPSE_TCONE,
  ANNULAR_TCONE,
  ANNULAR_CONE_X,
  PRISM,
  SBIN,
};

// simple mean...
//#define MU_FA(MU_ICV0,MU_ICV1) (.5*((MU_ICV0)+(MU_ICV1)))
//
#define SP_VOL_FA(RHO_ICV0,RHO_ICV1) (2.0/((RHO_ICV0)+(RHO_ICV1)))

// harmonic mean...
#define MU_FA(MU_ICV0,MU_ICV1) (2.0/(1.0/(MU_ICV0)+1.0/(MU_ICV1)))
//#define SP_VOL_FA(RHO_ICV0,RHO_ICV1) (0.5*(1.0/(RHO_ICV0)+1.0/(RHO_ICV1)))

// convenience macro to loop over the boundary conditions...
#define FOR_BCZONE for(vector<VofBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)

class VofSolver : public FlowSolver {

public:

  double rho_l,rho_g; // liquid and gas densities, repsectively
  double mu_l,mu_g;   // liquid and gas viscosities, repsectively
  double gravity[3];
  double sigma;       // surface tension coefficient
  double p_ref;
  double T_ref;
  double R_gas; // ideal gas
  double gamma;

  double *vof;           // volume fraction
  double *vof0;           // volume fraction
  double (*u)[3];        // velocity
  double (*u0)[3];        // velocity
  double *p;             // pressure
  double *rho;           // density
  double *mu_lam;        // laminar viscosity
  double *mu_sgs;        // sgs eddy viscosity
  double *kappa;         // interface curvature
  double (*n)[3];        // interface normal
  double *g;             // interface distance g : n.(x_interface-x_vv)=g
  double (*sp_dpdx)[3];    // density-weighted cv pressure gradient

  // THINC sharpness parameter ...
  double *beta;
  int rk_step;
  int *iband;

  double *div;          // velocity divergence...
  double (*dudx)[3][3]; // velocity gradient
  double *q_fa;         // conservative face velocity...
  double *q0_fa;        // prev conservative face velocity...


  // 2-level plic area-weighted paraboloid centroid fit...
  list<double> * plicPointList; // 3 x number of edge intersections b/w PLIC and cv (ordered)
  double (*plic_xc)[3]; // plic barycenter
  double *plic_area; // plic area


  int *cv_flag;      // flag interface cells
  int *cv_flag_real;
  int nftb; // number of faces touching boundary
  int *fa_ftb; // face index for a face touching boundary (not in fpnpocv_i/v)

  string sgs_model; // subgrid scale model name
  int    p_maxiter; // max pressure iterations
  double p_zero;    // pressure zero
  int    u_maxiter; // max velocity iterations
  double u_zero;    // velocity zero
  double u_relax;   // relaxation for jacobi solver
  double p_relax;   // relaxation for pressure solver
  double vof_zero;  // min geometrically-supported volume fraction
  bool use_xp;      // moves centroid w/ u for extra accuracy, but may reduce cfl
  int normal_iters; // iterations in least-squares procedure...
  int kappa_iters;  // iterations in kappa smoothing ...
  int g_maxiter;    // max volume enforcement iterations
  int drop_interval; // drop counting interval ...
  bool probe_first; // probing SMD

  int8* icv_global_g; // global indices for ghosts...

  int pressure_solver;
  // Trilinos solver...
  TrilinosSolver * trilinosSolverForPressure;
  int trilinos_idata[3];
  bool rebuildPoissPrec;

  Multigrid::AlgebraicMultigrid* am;

  // boundary conditions...

  vector<VofBc*> bcs;
  map<string,VofBc*> bc_map;
  map<string,pair<int,VofBc*> > bc_queries;
  double sum_outlet_proj_area;

  // lsp stuff ...
  LpClass<LspState> * lpCls;
  double *lspMassSource;
  double (*lspMomentumSource)[3];
  bool lsp;

  CtiMaterial * material;
  CtiLiquid * fuel;
  CtiSolid  * dust;
  bool SOL_LP, LIQ_LP;
  bool evaporation;

  //===========
  // injectors
  //============
  list<LspInjector*> injectorList;
  double time0;

  // particle cv link
  int * paocv_i;        // particle of cv index
  vector<int> paocv_v;  // particle of cv value

  list<LspStats*> lspStatsList; // lsp stats
  list<LspPDF> lspPDFList; // lsp stats
  list<DropPDF> dropPDFList; // lsp stats

  int BreakupModel;
  double weber_cr;
  double core_br;

  double bu_k1; // k1 in breakup model
  double bu_k2; // k2 in breakup model
  double npar_target_multiplier; // npar target in breakup

  //==================
  // lsp agglemoration
  //==================
  bool lspAgglem;
  double MAX_PPC, MAX_NPAR;

  // boundary conditions
  vector<lpBaseBc*> lpZoneBcVec;

  //==================
  // lsp size dist
  //==================
  bool show_histogram;
  string hist_file;
  int hist_interval;
  int hist_nbins;
  bool hist_verbose_count;
  bool hist_verbose_mass;

  VofSolver() {

    COUT1("VofSolver()");

    // field data registration...

    // we only read need u and vof to restart calc...

    vof      = NULL; registerCvData(vof      , "vof"      , READWRITE_DATA);
    vof0     = NULL; registerCvData(vof0     , "vof0"     , CAN_WRITE_DATA);
    u        = NULL; registerCvData(u        , "u"        , READWRITE_DATA);
    u0       = NULL; registerCvData(u0       , "u0"       , CAN_WRITE_DATA);
    p        = NULL; registerCvData(p        , "p"        , CAN_WRITE_DATA);
    rho      = NULL; registerCvData(rho      , "rho"      , CAN_WRITE_DATA);
    mu_lam   = NULL; registerCvData(mu_lam   , "mu_lam"   , CAN_WRITE_DATA);
    mu_sgs   = NULL; registerCvData(mu_sgs   , "mu_sgs"   , CAN_WRITE_DATA);
    kappa    = NULL; registerCvData(kappa    , "kappa"    , CAN_WRITE_DATA);
    n        = NULL; registerCvData(n        , "n"        , CAN_WRITE_DATA);
    g        = NULL; registerCvData(g        , "g"        , CAN_WRITE_DATA);
    q_fa     = NULL; registerSignedFaData(q_fa     , "q_fa"     , READWRITE_DATA);
    q0_fa    = NULL; registerSignedFaData(q0_fa    , "q0_fa"    , READWRITE_DATA);

    dudx = NULL;
    sp_dpdx  = NULL;   registerCvData(sp_dpdx   , "sp_dpdx"    , READWRITE_DATA);
    div = NULL;


    plicPointList = NULL;
    plic_xc = NULL;
    plic_area = NULL;

    cv_flag = NULL;
    cv_flag_real = NULL;
    registerCvData(cv_flag_real, "cv_flag_real", READWRITE_DATA);


    nftb = 0;
    fa_ftb = NULL;

    icv_global_g = NULL;
    trilinosSolverForPressure = NULL;

    //hyperbolic tangent sharpness parameter
    beta = NULL;
    iband = NULL;

    // moved global parameters in constructor so that we had access to them earlier than during init...

    rho_l = getDoubleParam("rho_l");
    rho_g = getDoubleParam("rho_g");

    COUT2(" > rho_l: " << rho_l);
    COUT2(" > rho_g: " << rho_g);

    mu_l  = getDoubleParam("mu_l");
    mu_g  = getDoubleParam("mu_g");

    COUT2(" > mu_l: " << mu_l);
    COUT2(" > mu_g: " << mu_g);

    sigma = getDoubleParam("sigma");

    COUT2(" > sigma: " << sigma);

    p_ref = getDoubleParam("P_REF", 101325.0);
    T_ref = getDoubleParam("T_REF", 290.0);

    R_gas        = p_ref/rho_g/T_ref; // rho_gas is used as a reference
    gamma   = getDoubleParam("GAMMA", 1.4);


    if (Param * param = getParam("gravity")) {
      FOR_I3 gravity[i] = param->getDouble(i);
    }
    else {
      FOR_I3 gravity[i] = 0.0;
    }

    COUT2(" > gravity: " << COUT_VEC(gravity));

    p_zero = getDoubleParam("P_ZERO",1.0E-8);
    p_maxiter = getIntParam("P_MAXITER",2000);
    u_zero = getDoubleParam("U_ZERO",1.0E-8);
    u_relax = getDoubleParam("U_RELAX", 0.7);
    p_relax = getDoubleParam("P_RELAX", 0.7);
    u_maxiter = getIntParam("U_MAXITER",1000);
    vof_zero = getDoubleParam("VOF_ZERO",1.0E-8);
    use_xp = getBoolParam("USE_XP",true);
    g_maxiter = getIntParam("G_MAXITER",100);

    parseSgsModel(sgs_model); // TODO may need to promote for models that need registered data

    // need more than one iterations to get second-order accuracy (1-3 seems to be good)...
    normal_iters = getIntParam("NORMAL_ITERS",2);
    kappa_iters = getIntParam("KAPPA_ITERS", 15);

    rk_step = getIntParam("RK_STEP",2);
    assert(rk_step <= 2);

    drop_interval = getIntParam("DROP_INTERVAL",500);
    probe_first = true;

    // Lsp stuff...
    lsp = false;
    lsp = getBoolParam("LSP",false);
    
    if (lsp) {
      COUT1("LSP is on");
      // create the lp (not allocating at this stage)
      lpCls   = new LpClass<LspState>; // sets the lp to NULL, np and NP_MAX to 0

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

      material= NULL;
      fuel    = NULL;
      dust    = NULL;
      lspMassSource = NULL;
      lspMomentumSource = NULL;
      SOL_LP  = LIQ_LP = false;

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

      //// hacking the fuel right now
      //Param * param = getParam("LSP.MATERIAL");
      //if (param == NULL) CERR("could not find param LSP.MATERIAL");
      //string token = param->getString();
      //if (token == "LIQUID") {
      //  fuel = newCtiLiquid(param->getString(1),101325);
      //  material = fuel;
      //  LIQ_LP = true;
      //  IF_RANK0 fuel->info(300,350,10);
      //} else if (token == "SOLID") {
      //  dust = newCtiSolid(param->getString(1));
      //  IF_RANK0 dust->info(300,350,10);
      //  material = dust;
      //  SOL_LP = true;
      //} else {
      //  CERR("token for LSP.MATERIAL was not recognized!");
      //}
      
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
      // lsp agglemoration
      // ===========================
      param = getParam("LSP.AGGLEM");
      if(param != NULL) {
        lspAgglem = true;
        MAX_PPC  = param->getDouble(0);
        MAX_NPAR = param->getDouble(1);
      } else {
        lspAgglem = false;
        IF_RANK0 cout << " > WARNING: running without lspAgglemoration..." << endl;
      }

      show_histogram = false;
      hist_file = "";
      hist_interval = -1;
      hist_nbins = 120; // default value for bin counts
      hist_verbose_count = false;
      hist_verbose_mass = false;
      param = getParam("LSP.HISTOGRAM");
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
              CWARN("intarval for particle diameter histogram is not a multiply of checking interval");
          } else if (token == "NBINS") {
            hist_nbins = param->getDouble(iarg++);
          } else if (token == "V_COUNT") {
            hist_verbose_count = true;
          } else if (token == "V_MASS") {
            hist_verbose_mass = true;
          } else {
            CERR("LSP.SHOW_HISTOGRAM param not correctly specified");
          }
        }
      }

      evaporation = getBoolParam("EVAPORATION",false);
      npar_target_multiplier = getDoubleParam("LSP.NPAR_MULTIPLIER_BU", 1.0);
      COUT1(" > npar target multiplier: " << npar_target_multiplier);
    }

  }

  virtual ~VofSolver() {

    COUT1("~VofSolver()");

    DELETE(vof);
    DELETE(vof0);
    DELETE(u);
    DELETE(u0);
    DELETE(q_fa);
    DELETE(q0_fa);

    DELETE(p);
    DELETE(rho);
    DELETE(mu_lam);
    DELETE(mu_sgs);
    DELETE(kappa);
    DELETE(n);
    DELETE(g);


    DELETE(plicPointList);
    DELETE(plic_xc);
    DELETE(plic_area);


    DELETE(dudx);
    DELETE(sp_dpdx);
    DELETE(div);
    DELETE(cv_flag);

    DELETE(beta);
    DELETE(iband);

    FOR_BCZONE delete *it;

    DELETE(fa_ftb);

    DELETE(icv_global_g);
    if (trilinosSolverForPressure) {
      delete trilinosSolverForPressure;
      trilinosSolverForPressure = NULL;
    }

    if (am) {
      delete am;
      am = NULL;
    }

    if (lsp) {
      DELETE(lspMassSource);
      DELETE(lspMomentumSource);
      DELETE(paocv_i);
      if (material != NULL) {material = NULL;}
      if (fuel != NULL) {delete fuel; fuel = NULL;}
      if (dust != NULL) {delete dust; dust = NULL;}
      if (lpCls != NULL) {delete lpCls; lpCls = NULL;}

      for (list<LspInjector*>::iterator iter = injectorList.begin(); iter != injectorList.end(); ++iter) {
        delete (*iter);
      }

      for (int izone = 0; izone < lpZoneBcVec.size(); izone++) {
        delete lpZoneBcVec[izone];
      }
    }

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
          bcs.push_back(new SlipWallVBc(&bfZoneVec[izone],this));
        }
        else if ( bc_type == "WALL") {
          bcs.push_back(new WallVBc(&bfZoneVec[izone],this));
        }
        else if ( bc_type == "INLET") {
          bcs.push_back(new InletVBc(&bfZoneVec[izone],this));
        }
        else if ( bc_type == "OUTLET") {
          bcs.push_back(new OutletVBc(&bfZoneVec[izone],this));
        }
        else if ( bc_type == "HOOK") {
         //
        }
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

    // ensure that the error handling is synced..
    // should be un-necessary and we'll check below

    MPI_Bcast(&nerr,1,MPI_INT,0,mpi_comm);
    assert( nerr == int(errors.size()));
    reportBcErrors(errors);
  }

  void initBoundaryConditions() {

    // allow the boundary conditions to allocate their necessary data (not allowed to
    // set anything here, because it may be coming from the restart data yet to be read)..

    FOR_BCZONE (*it)->initData();

    // also create a map based on the boundary conditions. provides easy access to
    // grab the bc object by name without looping through (O(N_boundary conditions))

    FOR_BCZONE bc_map[(*it)->getName()] = *it;
  }

  virtual VofBc* initHookBc(BfZone* p) {
    CERR(" > user must supply a hook boundary condition for zone: " << p->getName());
  }

  VofBc* getBc(const string& name) const {
    map<string,VofBc*>::const_iterator it = bc_map.find(name);
    if ( it == bc_map.end()) {
      return NULL;
    } else {
      return it->second;
    }
  }

  void initData() {

    assert(vof   == NULL); vof   = new double[ncv_g];
    assert(vof0  == NULL); vof0   = new double[ncv_g];
    assert(u     == NULL); u     = new double[ncv_g][3];
    assert(u0    == NULL); u0     = new double[ncv_g][3];
    assert(q_fa  == NULL); q_fa  = new double[nfa];
    assert(q0_fa == NULL); q0_fa = new double[nfa];

    assert(p        == NULL); p        = new double[ncv_g];
    assert(rho      == NULL); rho      = new double[ncv_g];
    assert(mu_lam   == NULL); mu_lam   = new double[ncv_g];
    assert(mu_sgs   == NULL); mu_sgs   = new double[ncv_g];

    assert(kappa    == NULL); kappa    = new double[ncv_g];
    assert(n        == NULL); n        = new double[ncv_g][3];
    assert(g        == NULL); g        = new double[ncv_g];


    assert(plicPointList == NULL); plicPointList = new list<double>[ncv];
    assert(plic_xc == NULL); plic_xc = new double[ncv_g][3];
    assert(plic_area == NULL); plic_area = new double[ncv_g];

    assert(dudx == NULL);        dudx = new double[ncv_g][3][3];
    assert(sp_dpdx == NULL);     sp_dpdx = new double[ncv][3];
    assert(div == NULL);         div = new double[ncv];
    assert(cv_flag == NULL);     cv_flag = new int[ncv_g];
    assert(cv_flag_real == NULL); cv_flag_real = new int[ncv];


    icv_global_g = new int8[ncv_g-ncv];

    assert(beta == NULL); beta = new double[ncv_g];
    assert(iband == NULL); iband = new int[ncv_g];

    am = NULL;

    //lsp stuff

    if (lsp) {
      //===================
      // should now have the correct np form restart
      // allocate the lp inside lsp
      //===================
      lpCls->init();

      //=================
      // particle Cv link
      //=================
      paocv_i = new int [ncv+1];

      //=================
      // lsp injectors
      //=================
      setupInjectors();

      // only used for mass trasfer
      lspMassSource = new double [ncv_g];
      lspMomentumSource = new double [ncv_g][3];
      FOR_ICV_G {
        lspMassSource[icv] = 0.0;
        FOR_I3 lspMomentumSource[icv][i] = 0.0;
      }
    }

    // compute local sharpness parameter for THINC ...
    double h_avg = 0.0;
    double beta_c;
    beta_c = getDoubleParam("BETA",3.0);
    FOR_ICV_G {
      h_avg = pow((6.0*vol_cv[icv]),1.0/3.0);
      beta[icv] = beta_c/h_avg;
    }

  }

  void init() {

    wtime = MPI_Wtime();

    //FlowSolver::init();
    StaticSolver::init(INIT_COMPACT_FACES|INIT_CV_GRAD);

    // needed by app to control contex menu
    added_data_fields.insert("interfaces");

    // needed for curvature...
    updateCvData(x_vv,REPLACE_TRANSLATE_DATA); // TODO just make x_vv updated in ghosts in StaticSolver

    // needed for outlets...
    sum_outlet_proj_area = 0.0;
    FOR_BCZONE {
      Param * param = getParam((*it)->getName());
      if ( param->getString(0) == "OUTLET")
        sum_outlet_proj_area += MAG((*it)->zone_ptr->n_global);
    }

    // needed for trilinos...
    updateCvDataSeparateGhosts(icv_global,icv_global_g);

    trilinosSolverForPressure = NULL;
    setPressureSolver(pressure_solver);

    initDropStats();

    // needed for lsp
    if (lsp) {
      added_data_fields.insert("particles");
      initLspStats();
      if (cvAdt == NULL) buildCvAdt();
      initLspBc();
    }

    // overwrite area over delta with forming point version...
    FOR_IFA {
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      const double area = MAG(n_fa[ifa]);
      const double delta = DIST(x_vv[icv0],x_vv[icv1]);
      area_over_delta_fa[ifa] = area/delta;
    }

    if (mpi_rank == 0)
      logger->setKillFilename("killcharles");
  }

  inline int8 getIcvGlobal(const int icv) const {
    assert((icv >= 0)&&(icv < ncv_g));
    assert(icv_global && icv_global_g);
    if (icv < ncv)
      return icv_global[icv];
    else
      return icv_global_g[icv-ncv];
  }

  void initFromParams() {
    // parameter-based initialization of data...
    StaticSolver::initFromParams();
  }

  void initialHookBcs() {

    // not exactly clear -- there is a bunch of non-virtual initializations for
    // volume and bf data that should be handled without requiring the user
    // to write a unique "initialHook": e.g. take nearby cv values to initialize
    // surface data in NSCBC-type bcs when not set from data file, etc...

    FOR_BCZONE (*it)->initialHook();
  }

  void initComplete() {

    // if we do not have the rho/mu_lam fields, then initialize the rho/mu_lam fields from the vof field...

    COUT1("initComplete()");

    updateCvData(vof); // assume vof has been read or set
    updateCvData(u);
    updateBaseState();

    // get mass conserving velocity..

    FOR_ICV_G p[icv] = 0.0; // for initial guess...
    if (!checkDataFlag("q_fa")) {
      assert(!checkDataFlag("q0_fa"));
      //assert(!checkDataFlag("sp_dpdx"));
      // project out divergence in u...
      FOR_ICV_G kappa[icv] = 0.0;
      double (*u_copy)[3] = new double[ncv][3];
      FOR_ICV FOR_I3 u_copy[icv][i] = u[icv][i];
      const double dt_copy = dt;
      dt = 1.0;
      FOR_ICV_G kappa[icv] = 0.0;
      FOR_ICV_G cv_flag[icv] = 0;
      FOR_ICV FOR_I3 sp_dpdx[icv][i] = 0.0;
      FOR_ICV_G plic_area[icv] = 0.0;
      dumpRange(u,ncv,"u before correction");
      solvePAndCorrectU();
      dt = dt_copy;
      FOR_ICV FOR_I3 u[icv][i] = u_copy[icv][i];
      delete[] u_copy;
      FOR_IFA q0_fa[ifa] = q_fa[ifa]; // first step is forward euler...
    }

    // calculate the interface properties...
    if (!checkDataFlag("n")) {
      assert(!checkDataFlag("g"));
      assert(!checkDataFlag("kappa"));
      updateInterface();
    }

    StaticSolver::calcCvGrad(dudx,u);
    updateCvData(dudx);

    if (lsp) {
      time0 = time;
      for (int ip = 0; ip < lpCls->size(); ++ip) {
        lpCls->lp[ip].flag = KEEP;
        updateDp(ip);
      }
    }
    // lsp data initialization
//    if (lsp) {
//      updateLspPrimitiveData();
//      //setPaoCv();
//      setLspX0();
//      recycleLsp();
//      relocateLsp<LspState>();
//    }

    // subgrid stress...
    calcSgs();

    // update xp...
    //
    buildIband();

    if (lsp) updateLspDrop();

    // report mass and momentum...

    report();

    /*
    if (lsp) {
      time0 = time;
      if (lp_flag) {
        for (int ip = 0; ip < lpCls->size(); ++ip) {
          lpCls->lp[ip].flag = KEEP;
          updateDp(ip);
        }
      }
    }
*/

  }

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

  void setPressureSolver (int& pressure_solver) {

    Param * param = getParam("PRESSURE_SOLVER");

    if (param != NULL) {
      string name = param->getString();
      if (name == "TRILINOS") {
        pressure_solver = TRILINOS_SOLVER;
        trilinos_idata[0] = getIntParam("OPTION", 4); //option
        trilinos_idata[1] = getIntParam("MAX_LEVEL",5); //max ml lvels
        trilinos_idata[2] = getIntParam("MAX_REPART",200); //maximum repart level;
        if (mpi_rank == 0)
          cout << " > PRESSURE_SOLVER = TRILINOS TOL=" << p_zero <<
            " OPTION=" << trilinos_idata[0] <<
            " MAXITER=" << p_maxiter <<
            " ML_MAX_LEVELS = " << trilinos_idata[1] <<
            " REPART_LEVEL = " << trilinos_idata[2] << endl;
      }
      else if (name == "BCGSTAB") {
        pressure_solver = BCGSTAB_SOLVER;
        COUT1(" > PRESSURE_SOLVER = BSGSTAB TOL=" << p_zero <<
            " MAXITER=" << p_maxiter );
      }
      else if (name == "AMG") {
        pressure_solver = AMG_SOLVER;
        COUT1(" > PRESSURE_SOLVER = AMG TOL=" << p_zero <<
            " MAXITER=" << p_maxiter );
      }
      else if (name == "JACOBI") {
        pressure_solver = JACOBI_SOLVER;
        COUT1(" > PRESSURE_SOLVER = JACOBI TOL=" << p_zero << "MAX_ITER= " << p_maxiter );
      }
      else {
        if (mpi_rank == 0 ) {
          cout << name << " is not the supported solver. please use either TRILINOS or BCGSTAB" << endl;
        }
        throw(0);
      }
    }
    else {
      pressure_solver = BCGSTAB_SOLVER;
    }

  }

  int advanceSolution() {

    if ( (mpi_rank==0) && (step%check_interval==0) )
      cout << " > advanceSolution()" << endl;

    CtiRegister::clearCurrentData();

    FOR_ICV_G vof0[icv] = vof[icv];
    FOR_ICV_G FOR_I3 u0[icv][i] = u[icv][i];

    // advance time step...
    time += dt;

    // inject/breakup/advance any Lagrangian spray...
    if (lsp) advanceLsp(time,dt,step,step%check_interval==0);

    // apply flux bcs...
    FOR_BCZONE (*it)->setBc();

    timer.split("set bc's");

    for (int iter = 0; iter < rk_step ; iter++) {
      if ( (mpi_rank==0) && (step%check_interval==0) )
        cout << " > RK sub_step: " << iter+1 << endl;

      calcRhs();
      updateInterface();

      FOR_ICV FOR_I3 u[icv][i] += dt*gravity[i];
      updateCvData(u);

      calcCurvature();
      // pressure correction...
      if (getBoolParam("SOLVER",true)) solvePAndCorrectU();

      // velocity gradient...

      StaticSolver::calcCvGrad(dudx,u);
      updateCvData(dudx);

      calcSgs();
      timer.split("calc dudx");
    } // rk step ends..

    if (rk_step == 2) {
      FOR_ICV_G vof[icv] = 0.5*(vof0[icv]+vof[icv]);
      FOR_ICV_G FOR_I3 u[icv][i] = 0.5*(u0[icv][i]+u[icv][i]);
      updateBaseState();
      updateInterface();
    }

    buildPlic();
    updateDrop();

    doProbes();
    report();

    return 0;
  }

  void report() {

    if (step%check_interval == 0) {
      double my_buf[5] = {0.0,0.0,0.0,0.0,0.0};

      FOR_ICV {
        my_buf[0] += vof[icv]*vol_cv[icv];
        my_buf[1] += rho[icv]*vol_cv[icv];
        my_buf[2] += rho[icv]*u[icv][0]*vol_cv[icv];
        my_buf[3] += rho[icv]*u[icv][1]*vol_cv[icv];
        my_buf[4] += rho[icv]*u[icv][2]*vol_cv[icv];
      }

      double buf[5];
      MPI_Reduce(my_buf,&buf,5,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      if ( mpi_rank==0) {
        cout << " > liquid volume: " << buf[0];
        cout << ", total mass: " << buf[1];
        cout << ", momentum: " << COUT_VEC(&buf[2]) << endl;
      }
      dumpCfl();
      dumpRange(vof,ncv,"vof");
      dumpRange(u,ncv,"u");
      queryVof();
      queryBcs();
    }

    if (lsp) {
      reportLsp();
      int lspwrite_interval = getDoubleParam("LSP_WRITE",1000);
      if (step%lspwrite_interval == 0 ) writeLspTecplot(step);
    }
    
  }

  void updateDrop(){

    bool dropReport = false;
    for (list<DropPDF>::iterator pdfs = dropPDFList.begin(); pdfs != dropPDFList.end(); ++pdfs) {
      if ((step%pdfs->getInterval() == 0) && (step>=pdfs->getStart())) {
        dropReport = true;
        break;
      }
    }

    const int lsp_interval = getIntParam("LSP_TRANSFER_INTERVAL", 100);

    if (dropReport || (step%lsp_interval == 0) ) buildIband();

    if (dropReport) dropIdentification();
    if (step%lsp_interval == 0 && lsp) updateLspDrop();

  }

   
  void queryVof() {
    FOR_PARAM_MATCHING("QUERY_VOF") {
      processQueryVof(&(*param));
    }
  }

  void doProbes() {
    //custom probe here...
  }

  void processQueryVof(Param * param) {

    int interval = 1;

    bool b_name = false;
    string name;
    bool b_xp = false;
    int geom_type = 0;
    double geom_data[9];
    bool kappa_range = false;
    double bound[2];
    double xp[3];
    bool b_np = false;
    double np[3];

    int iarg = 0;
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "INTERVAL") {
        interval = param->getInt(iarg++);
        if (step%interval != 0)
          return;
      }
      else if (token == "NAME") {
        b_name = true;
        name = param->getString(iarg++);
      }
      else if (token == "GEOM") {
        b_xp = true;
        const string geom = param->getString(iarg++);
        if (geom == "BOX") {
          geom_type = BOX;
          for (int i=0; i<6; i++) geom_data[i] = param->getDouble(iarg++);
        }
        else if (token == "ANNULAR_CONE_X") {
          //---------------------------------------------------------------
          // HCP_WINDOW ANNULAR_TCONE <x0> <y0> <z0> <x1> <y1> <z1> <r0> <r1>
          // ---------------------------------------------------------------
          geom_type = ANNULAR_TCONE;
          geom_data[0] = param->getDouble(iarg++);
          geom_data[1] = param->getDouble(iarg++);
          geom_data[2] = param->getDouble(iarg++);
          geom_data[3] = param->getDouble(iarg++);
          geom_data[4] = param->getDouble(iarg++);
          geom_data[5] = param->getDouble(iarg++);
          geom_data[6] = param->getDouble(iarg++);
          geom_data[7] = param->getDouble(iarg++);
        }

        else {
          CERR("unsupported GEOM: " << token);
        }
      }
      else if (token == "RANGE") {
        kappa_range = true;
        bound[0] = param->getDouble(iarg++);
        bound[1] = param->getDouble(iarg++);
      }
      else {
        if (mpi_rank == 0) cout << "Warning: skipping unrecognized QUERY_VOF token: " << token << endl;
      }
    }

    int ierr = 0;
    if (!b_name) {
      if (mpi_rank == 0) cout << "Warning: QUERY_VOF missing NAME <string>" << endl;
      ierr = -1;
    }
    if (!b_xp) {
      if (mpi_rank == 0) cout << "Warning: QUERY_VOF missing GEOM BOX <x0> <x1> <y0> <y1> <z0> <z1> " << endl;
      ierr = -1;
    }
    if (ierr != 0)
      return;

    if (probe_first) {
      if (mpi_rank == 0)
        cout << "[QUERY_VOF:" << name<< "]# 1:step 2:time 3:total volume 4:liquid volume 5:vof_avg 6:mixing_variance 7:liquid surface  8:SMD32 (6volume/area) " << endl;
      probe_first = false;
    }

    // sum the plic area in a given area ...
    double my_stats_sum[3] = {0.0,0.0,0.0};
    FOR_ICV {
      if (geom_type == BOX) {
        if ( x_cv[icv][0] >= geom_data[0] && x_cv[icv][0] <= geom_data[1] ) {
          if ( x_cv[icv][1] >= geom_data[2] && x_cv[icv][1] <= geom_data[3] ) {
            if ( x_cv[icv][2] >= geom_data[4] && x_cv[icv][2] <= geom_data[5] ) {
              my_stats_sum[0] += vol_cv[icv];
              my_stats_sum[1] += vof[icv]*vol_cv[icv];
              if (kappa_range) {
               if (kappa[icv] >= bound[0] && kappa[icv] <= bound[1])  my_stats_sum[2] += plic_area[icv];
              }
              else {
                my_stats_sum[2] += plic_area[icv];
              }
            }
          }
        }
      }
      else if (geom_type == ANNULAR_CONE_X) {
        double radius = sqrt(x_cv[icv][1]*x_cv[icv][1] + x_cv[icv][2]*x_cv[icv][2]);
        if (x_cv[icv][0] >= geom_data[0] && x_cv[icv][0] <= geom_data[3] ) {
          if (radius >= geom_data[6] && radius <= geom_data[7]) {
              my_stats_sum[0] += vol_cv[icv];
              my_stats_sum[1] += vof[icv]*vol_cv[icv];
              if (kappa_range) {
                if (kappa[icv] >= bound[0] && kappa[icv] <= bound[1])  my_stats_sum[2] += plic_area[icv];
              }
              else {
                my_stats_sum[2] += plic_area[icv];
              }
          }
        }
      }
    }

    double stats_sum[3];
    MPI_Allreduce(my_stats_sum,stats_sum,3,MPI_DOUBLE,MPI_SUM,mpi_comm);

    const double vof_avg = stats_sum[1]/stats_sum[0];
    double my_var_sum = 0.0;
    FOR_ICV {
      if (geom_type == BOX) {
        if ( x_cv[icv][0] >= geom_data[0] && x_cv[icv][0] <= geom_data[1] ) {
          if ( x_cv[icv][1] >= geom_data[2] && x_cv[icv][1] <= geom_data[3] ) {
            if ( x_cv[icv][2] >= geom_data[4] && x_cv[icv][2] <= geom_data[5] ) {
              my_var_sum += (vof[icv]-vof_avg)*(vof[icv]-vof_avg)*vol_cv[icv];
            }
          }
        }
      }
      else if (geom_type == ANNULAR_CONE_X) {
        double radius = sqrt(x_cv[icv][1]*x_cv[icv][1] + x_cv[icv][2]*x_cv[icv][2]);
        if (x_cv[icv][0] >= geom_data[0] && x_cv[icv][0] <= geom_data[3] ) {
          if (radius >= geom_data[6] && radius <= geom_data[7]) {
            my_var_sum += (vof[icv]-vof_avg)*(vof[icv]-vof_avg)*vol_cv[icv];
          }
        }
      }
    }


    double var_sum = 0.0;
    MPI_Reduce(&my_var_sum,&var_sum,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);


    if (mpi_rank == 0) {
      cout << "[QUERY_VOF:" << name << "] " << step << " " << time << " " <<
        stats_sum[0] << " " <<
        stats_sum[1] << " " << vof_avg << " " << var_sum/stats_sum[0] << " " << stats_sum[2] << " " <<
        6.0*stats_sum[1]/stats_sum[2]  << " " << endl;
    }

  }

  void queryBcs() {


    // QUERY_BC <bc-name> [INTERVAL = <interval>] [WRITE]
    FOR_PARAM_MATCHING("QUERY_BC") {

      // check if the bc is matched against a known query, we will
      // use the entire param string as the key for the bc_query map

      map<string,pair<int,VofBc*> >::iterator bc_it = bc_queries.find(param->str());

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
          }
          else if (token == "WRITE") {
            b_write = true;
          }
          else {
            CERR( " Invalid query_bc syntax; QUERY_BC [INTERVAL <interval>] [WRITE]");
          }
        }

        VofBc* bc = getBc(bc_name);

        if ( bc == NULL) {

          CWARN(" > unable to find boundary zone for QUERY: " << bc_name);

        }
        else {

          if (bc->ss == NULL) bc->ss = new std::stringstream();
          bc->b_write = b_write;

          pair<map<string,pair<int,VofBc*> >::iterator,bool> ret =
            bc_queries.insert(pair<string,pair<int,VofBc*> >( param->str(),
                  pair<int,VofBc*>(interval,bc)));

          assert( ret.second); // ensure that the query was properly inserted.
          bc_it = ret.first;
        }
      }

      if ( bc_it != bc_queries.end() ) {

        const int query_interval = bc_it->second.first;
        VofBc* bc          = bc_it->second.second;

        if ( step%query_interval == 0) {
          bc->query(bc_it->first);
        }

      }

    }
  }



  void updateLspDrop() {

    lspTransfer();

    // lsp statistics
    updateLspStats(lpCls->lp,lpCls->size(),time,dt,step,1);
  } 

  void updateInterface() {

    if ( step%check_interval==0 )
      COUT1("updateInterface()");

    cleanupVof();
    flagInterfaceCvs();

    calcNormal();
    //build Plic Surface and compute Area in plic_area[] ...
    //buildPlic();
    // normal vector is required to reconstruct hyperbolic tangent...
    calcGfromVof();

  }

  void calcCurvature() {

    //calcCurvatureSignedDistance();
    calcCurvatureSimple();
  }

  void calcCurvatureFromG() {


    double (*dGdx)[3] = new double[ncv][3];
    double (*dGxdx)[3][3] = new double[ncv][3][3];

    StaticSolver::calcCvGrad(dGdx, g);

    // calc Gxx, Gyy, Gzz,Gxy  etc
    StaticSolver::calcCvGrad(dGxdx, dGdx);

    FOR_ICV {
      // we need curvature in the vof nodes and the first band...
      if (cv_flag[icv] >= 1) {
        const double s0 =
          dGxdx[icv][0][0]*(dGdx[icv][1]*dGdx[icv][1]+dGdx[icv][2]*dGdx[icv][2]) +
          dGxdx[icv][1][1]*(dGdx[icv][0]*dGdx[icv][0]+dGdx[icv][2]*dGdx[icv][2]) +
          dGxdx[icv][2][2]*(dGdx[icv][0]*dGdx[icv][0]+dGdx[icv][1]*dGdx[icv][1]);
        const double s1 = sqrt(DOT_PRODUCT(dGdx[icv],dGdx[icv]));
        kappa[icv] = ( s0 - 2.0*(dGxdx[icv][0][1]*dGdx[icv][0]*dGdx[icv][1] +
              dGxdx[icv][0][2]*dGdx[icv][0]*dGdx[icv][2] +
              dGxdx[icv][1][2]*dGdx[icv][1]*dGdx[icv][2]) )/(s1*s1*s1);
      }
      else {
        kappa[icv] = 0.0;
      }
    }
    // bound by cell length scale...
    FOR_ICV {
      if (kappa[icv] > 0.0) {
        kappa[icv] = min(kappa[icv],1.0/r_vv[icv]);
      }
      else if (kappa[icv] < 0.0) {
        kappa[icv] = max(kappa[icv],-1.0/r_vv[icv]);
      }
      if (cv_flag[icv] < 1 ) kappa[icv] = 0.0;
    }


    delete[] dGxdx;
    delete[] dGdx;

  }

  void calcCurvatureDirect() {

    if (cvAdt == NULL) buildCvAdt();

    vector<int> cvList;
    // direct front curvature computation ...
    FOR_ICV {
      if (cv_flag[icv] >= 1.0 ) {
        double xbase[3];
        double sum_weights = 0.0;
        double kappa_d = 0.0;
        double mag = DOT_PRODUCT(n[icv],n[icv]);
        //cout << "icv, n = " << icv << " " << mag << endl;
        FOR_I3  xbase[i] = x_cv[icv][i] + g[icv]*n[icv][i];
        //cout << "x = " << COUT_VEC(x_cv[icv]) << endl;
        //cout << "xbase= " << COUT_VEC(xbase) << endl;

       // cout << "Xbase=? "<< g[icv] << " " << COUT_VEC(n[icv]) << " " << " " << sqrt(x_cv[icv][0]*x_cv[icv][0]+x_cv[icv][1]*x_cv[icv][1]) << " " << sqrt(xbase[0]*xbase[0]+xbase[1]*xbase[1]) << endl;
        // caution ..not parallel for now.
        cvAdt->buildListForPoint(cvList, xbase);
        const int my_icv = getClosestICVFromList(cvList,xbase);
        //if (icv > 0 ) my_icvp[ip] = icv;
       // cout << "my_icv = " << icv << " " << vof[icv] << " "<<   my_icv << " " << vof[my_icv] << " " << COUT_VEC(x_cv[my_icv]) << " " << kappa[my_icv] << endl;
      //  assert(my_icv > 0 && my_icv < ncv); // assume points is in the same proc...not true for sure...
        if (my_icv > 0 ) {
          for (int coc = cvocv_i[my_icv]; coc != cvocv_i[my_icv+1]; ++coc) {
            const int icv_nbr = cvocv_v[coc];
            double d2 = DIST2(x_cv[icv_nbr], xbase);
            double this_weight = 1.0/(sqrt(d2) + 1.0E-15); // add tol to avoid 1/zero
            sum_weights += this_weight;
            kappa_d += this_weight*kappa[icv_nbr];
          }
          // cout << "kappa difference = " << kappa[icv] << " " << kappa_d/sum_weights << " " << kappa[my_icv] << " " << vof[icv] <<  endl;

          kappa[icv]  = kappa_d/sum_weights;
        }
      }
      else kappa[icv] = 0.0;
    }

    updateCvData(kappa);

  }


  void calcCurvatureSimple() {

    // calc kappa = - div.n...
    FOR_ICV kappa[icv] = 0.0;
    FOR_IFA {
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      double n_f[3];
      FOR_I3 n_f[i] = 0.5*(n[icv0][i]+n[icv1][i]);
      const double mag = MAG(n_f);
      if (mag > 1.0E-14) {
        const double flux = DOT_PRODUCT(n_f,n_fa[ifa])/mag;
        kappa[icv0] -= flux;
        if (icv1 < ncv)
          kappa[icv1] += flux;
      }
    }

    FOR_ICV kappa[icv] *= -inv_vol[icv];

    updateCvData(kappa);

    //calcCurvatureFromG(n);

    // iterate for smoothing curvature

    double *wt_sum = new double[ncv_g];
    double *kappa_s = new double[ncv_g];


    for (int iter =0 ; iter < kappa_iters; iter++) {

      FOR_ICV_G {
        wt_sum[icv] = 0.0;
        kappa_s[icv] = 0.0;
      }

      FOR_ICV {
        if (cv_flag[icv] >= 1.0 ) {
          for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) {
            const int icv_nbr = cvocv_v[coc];
            const double wt = pow(vof[icv_nbr]*(1.0-vof[icv_nbr]),1.0);
            kappa_s[icv] += wt*kappa[icv_nbr];
            wt_sum[icv] += wt;
          }
        }
      }

      FOR_ICV {
        if (wt_sum[icv] > 1.0E-10) {
          kappa[icv] = kappa_s[icv]/wt_sum[icv];
        }
        else kappa[icv] = 0.0;
      }

      updateCvData(kappa);

    }

    // bound by cell length scale...
    FOR_ICV {
      if (kappa[icv] > 0.0) {
        kappa[icv] = min(kappa[icv],1.0/r_vv[icv]);
      }
      else if (kappa[icv] < 0.0) {
        kappa[icv] = max(kappa[icv],-1.0/r_vv[icv]);
      }
    //  if (cv_flag[icv] < 1 ) kappa[icv] = 0.0;
      // zero out for the subgrid cell?
     // if (iband[icv] == 10) kappa[icv] = 0.0;
    }

    updateCvData(kappa);
    delete[] wt_sum;
    delete[] kappa_s;

  }

  void buildPlic() {


    if (step%check_interval == 0 ) COUT1(" > buildPlic() ..");

    buildPlicPolys();

  }

  void buildPlicPolys() {

    // buid plicPointList...

    FOR_ICV {
      plicPointList[icv].clear();
      if (cv_flag[icv] >= 1) {
        intersectCvWithPlic(plicPointList[icv],icv);
      }
    }

    // now calculate centroid...
    FOR_ICV {
      if (cv_flag[icv] >= 1) {
        plic_area[icv] = computePointListAreaAndCentroid(plic_xc[icv],plicPointList[icv]);
      }
      else {
        FOR_I3 plic_xc[icv][i] = 0.0;
        plic_area[icv] = 0.0;
      }
    }
    updateCvData(plic_area);
  }

  void intersectCvWithPlic(list<double>& points, const int icv) {

    for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
      const int ifa = faocv_v[foc];

      const int size0 = points.size();
      int ino0 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof1 = noofa_i[ifa]; nof1 != noofa_i[ifa+1]; ++nof1) {
        const int ino1 = noofa_v[nof1];
        const double g0 = g[icv] - ((x_no[ino0][0]-x_vv[icv][0])*n[icv][0] +
                                    (x_no[ino0][1]-x_vv[icv][1])*n[icv][1] +
                                    (x_no[ino0][2]-x_vv[icv][2])*n[icv][2]);
        const double g1 = g[icv] - ((x_no[ino1][0]-x_vv[icv][0])*n[icv][0] +
                                    (x_no[ino1][1]-x_vv[icv][1])*n[icv][1] +
                                    (x_no[ino1][2]-x_vv[icv][2])*n[icv][2]);
        if (g1*g0 < 0.0) {
          // plic intersects edge, so store intersection point...
          const double factor = g0/(g1-g0);
          FOR_I3 points.push_back(x_no[ino0][i]-factor*(x_no[ino1][i]-x_no[ino0][i]));
          if (points.size() == size0+6) {
            break;
          }
        }
        ino0 = ino1;
      }
      if (points.size() == size0+3) {
        FOR_I3 points.pop_back(); // just remove unpaired point
      }
      assert((points.size() == size0)||(points.size() == size0+6)); // added nothing or single edge

    }
    for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
      const int ibf = bfocv_v[boc];

      const int size0 = points.size();
      int ino0 = noobf_v[noobf_i[ibf+1]-1];
      for (int nob1 = noobf_i[ibf]; nob1 != noobf_i[ibf+1]; ++nob1) {
        const int ino1 = noobf_v[nob1];
        const double g0 = g[icv] - ((x_no[ino0][0]-x_vv[icv][0])*n[icv][0] +
                                    (x_no[ino0][1]-x_vv[icv][1])*n[icv][1] +
                                    (x_no[ino0][2]-x_vv[icv][2])*n[icv][2]);
        const double g1 = g[icv] - ((x_no[ino1][0]-x_vv[icv][0])*n[icv][0] +
                                    (x_no[ino1][1]-x_vv[icv][1])*n[icv][1] +
                                    (x_no[ino1][2]-x_vv[icv][2])*n[icv][2]);
        if (g1*g0 < 0.0) {
          // plic intersects edge, so store intersection point...
          const double factor = g0/(g1-g0);
          FOR_I3 points.push_back(x_no[ino0][i]-factor*(x_no[ino1][i]-x_no[ino0][i]));
          if (points.size() == size0+6) {
            break;
          }
        }
        ino0 = ino1;
      }
      if (points.size() == size0+3) {
        FOR_I3 points.pop_back(); // just remove unpaired point
      }
      assert((points.size() == size0)||(points.size() == size0+6)); // added nothing or single edge
    }
    assert(points.size()%6 == 0);
  }


  virtual void buildSolverSurface(vector<SimpleTri>& triVec,const string& surface_name) {

    if (surface_name == "PLIC") {
      FOR_ICV {
        if (cv_flag[icv] >= 1) {
          list<double>::iterator it = plicPointList[icv].begin();
          double x0[3],x1[3];
          while (it != plicPointList[icv].end()) {
            FOR_I3 {
              x0[i] = *it;
              ++it;
            }
            FOR_I3 {
              x1[i] = *it;
              ++it;
            }
            const double nA[3] = TRI_NORMAL_2(x0,x1,plic_xc[icv]);
            if (DOT_PRODUCT(nA,n[icv]) > 0)
              triVec.push_back(SimpleTri(x0,x1,plic_xc[icv]));
            else
              triVec.push_back(SimpleTri(x1,x0,plic_xc[icv]));
          }
        }
      }
    }
    else {
      CWARN(" > unrecognized VofSolver surface name " << surface_name << " in write image.");
    }

  }


  double computePointListAreaAndCentroid(double xc[3],list<double> points) {
    double area = 0.0;
    FOR_I3 xc[i] = 0.0;
    if (points.size() > 0) {
      double xp[3],x0[3],x1[3];
      list<double>::iterator it = points.begin();
      FOR_I3 {
        xp[i] = *it;
        ++it;
      }
      it = points.begin();
      while (it != points.end()) {
        FOR_I3 {
          x0[i] = *it;
          ++it;
        }
        FOR_I3 {
          x1[i] = *it;
          ++it;
        }
        const double nA[3] = TRI_NORMAL_2(x0,x1,xp);
        const double mag = MAG(nA);
        area += mag;
        FOR_I3 xc[i] += mag*(x0[i]+x1[i]+xp[i]);
      }
      FOR_I3 xc[i] /= 3.0*area;
      area *= 0.5;
    }
    return area;
  }

  void buildTransform(double R[3][3], double t[3], const double gcm, const double ncm[3], const double xcm[3]) {

    // translation is just the negative of the interface position

    FOR_I3 t[i] = -(gcm*ncm[i]+xcm[i]);

    // rotation based on coordinate system aligned with n...

    double u[3],v[3],w[3];
    FOR_I3 u[i] = ncm[i];
    const double pu[3] = {fabs(u[0]),fabs(u[1]),fabs(u[2])};

    if ( (pu[2] >= pu[1]) && (pu[2] >= pu[0]) ) {
      // z axis closest to n, use x and y for building v and w

      double y[3] = {0.0,1.0,0.0};
      double proj_yu = u[1]; // (0,1,0).(u0,u1,u2)/|u|^2 = u1
      FOR_I3 v[i] = y[i] - proj_yu*u[i];
      const double mag_v = MAG(v); assert(mag_v > 0.0);
      FOR_I3 v[i] /= mag_v;

      double x[3] = {1.0,0.0,0.0};
      double proj_xu = u[0]; // (1,0,0).(u0,u1,u2)/|u|^2 = u0
      double proj_xv = v[0]; // (1,0,0).(v0,v1,v2)/|v|^2 = v0
      FOR_I3 w[i] = x[i] - proj_xu*u[i] - proj_xv*v[i];
      const double mag_w = MAG(w); assert(mag_w > 0.0);
      FOR_I3 w[i] /= mag_w;

    }
    else if ( (pu[1] >= pu[0]) && (pu[1] >= pu[2]) ) {
      // y axis closest to n, use x and z for building v and w

      double z[3] = {0.0,0.0,1.0};
      double proj_zu = u[2]; // (0,0,1).(u0,u1,u2)/|u|^2 = u2
      FOR_I3 v[i] = z[i] - proj_zu*u[i];
      const double mag_v = MAG(v); assert(mag_v > 0.0);
      FOR_I3 v[i] /= mag_v;

      double x[3] = {1.0,0.0,0.0};
      double proj_xu = u[0]; // (1,0,0).(u0,u1,u2)/|u|^2 = u0
      double proj_xv = v[0]; // (1,0,0).(v0,v1,v2)/|v|^2 = v0
      FOR_I3 w[i] = x[i] - proj_xu*u[i] - proj_xv*v[i];
      const double mag_w = MAG(w); assert(mag_w > 0.0);
      FOR_I3 w[i] /= mag_w;

    }
    else {
      // x axis closest to n, use y and z for building v and w
      assert((pu[0] >= pu[1]) && (pu[0] >= pu[2]));

      double y[3] = {0.0,1.0,0.0};
      double proj_yu = u[1]; // (0,1,0).(u0,u1,u2)/|u|^2 = u1
      FOR_I3 v[i] = y[i] - proj_yu*u[i];
      const double mag_v = MAG(v); assert(mag_v > 0.0);
      FOR_I3 v[i] /= mag_v;

      double z[3] = {0.0,0.0,1.0};
      double proj_zu = u[2]; // (0,0,1).(u0,u1,u2)/|u|^2 = u2
      double proj_zv = v[2]; // (0,0,1).(v0,v1,v2)/|v|^2 = v2
      FOR_I3 w[i] = z[i] - proj_zu*u[i] - proj_zv*v[i];
      const double mag_w = MAG(w); assert(mag_w > 0.0);
      FOR_I3 w[i] /= mag_w;
    }

    FOR_I3 R[0][i] = u[i];
    FOR_I3 R[1][i] = v[i];
    FOR_I3 R[2][i] = w[i];
  }

  void applyTransform(double xp[3], const double R[3][3], const double t[3]) {
    double tmp[3] = {xp[0],xp[1],xp[2]};
    FOR_I3 tmp[i] += t[i];
    FOR_I3 xp[i] = DOT_PRODUCT(R[i],tmp);
  }

  void applyInverseRotation(double xp[3], const double R[3][3]) {
    const double tmp[3] = {xp[0],xp[1],xp[2]};
    FOR_I3 xp[i] = R[0][i]*tmp[0]+R[1][i]*tmp[1]+R[2][i]*tmp[2];
  }

  void applyRotation(double xp[3], const double R[3][3]) {
    double tmp[3] = {xp[0],xp[1],xp[2]};
    FOR_I3 xp[i] = DOT_PRODUCT(R[i],tmp);
  }


  void calcNormal() {

    calcNormalLeastSquares();
    //calcNormalSimple();

  }

  void calcNormalSimple() {

    // define smoothed vof field
    double (*vof_s) = new double[ncv_g];
    const double alpha = 0.1;

    //dumpRange(vof,ncv,"vof check");
    //cout << "pow= " << pow(1.0E-14,0.1) << " " << pow(-1.0E-14,0.1) << endl;

    FOR_ICV_G vof_s[icv ] = 0.0;

    FOR_ICV {
      vof_s[icv] = pow(vof[icv],alpha) / ( pow(vof[icv],alpha) + pow(1.0-vof[icv],alpha) );
      if (vof_s[icv] != vof_s[icv]) cout << "what?" << " " << icv << " "<< pow(vof[icv],alpha)<< " " <<pow(1.0-vof[icv],alpha)  << " " << ( pow(vof[icv],alpha) + pow(1.0-vof[icv],alpha) ) << " " << pow(vof[icv],alpha) << " " <<  vof_s[icv] << endl;
    }

    updateCvData(vof_s);

    dumpRange(vof_s,ncv,"vof_s");
    StaticSolver::calcCvGrad(n,vof_s);

    FOR_ICV {
      const double mag = MAG(n[icv]);
      if (mag > 0.0) {
        const double inv_mag = 1.0/mag;
        FOR_I3 n[icv][i] *= -inv_mag;
      }
    }
    updateCvData(n);

    dumpRange(n,ncv,"simple normal");
    delete[] vof_s;

  }
  void calcNormalLeastSquares() {

    // n = -grad(vof)/|grad(vof)|...

    StaticSolver::calcCvGrad(n,vof);
    updateCvData(n);

    FOR_ICV {
      const double mag = MAG(n[icv]);
      if (mag > 0.0) {
          FOR_I3 n[icv][i] = -n[icv][i]/mag;
      }
    }
    updateCvData(n);


    for (int iter = 0; iter < normal_iters; ++iter) {

      // calc (x_interface-x_cv).n = g...

      calcGfromVof();

      // use the above n/g as initial guess for least squares method...

      double R[3][3];
      double t[3];
      FOR_ICV {
        if (cv_flag[icv] >= 1) {

          // count number of well-defined nbr n's
          int np = 1;
          for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
            const int icv_nbr = cvocv_v[coc];
            if (cv_flag[icv_nbr] >= 1 && (DOT_PRODUCT(n[icv],n[icv_nbr]) > 0.0) ) np++; // fix for thin filaments
          }
          if (np >= 3) {

            // transform coordinate system to be aligned with n[icv] and centered on xg[icv]...

            buildTransform(R,t,g[icv],n[icv],x_cv[icv]);

            double sxx = 0.0;
            double syy = 0.0;
            double sxy = 0.0;
            double sxz = 0.0;
            double syz = 0.0;
            for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) { // everything taken wrt xg...
              const int icv_nbr = cvocv_v[coc];
              if (cv_flag[icv_nbr] >= 1 && (DOT_PRODUCT(n[icv],n[icv_nbr]) > 0.0) ) { // fix for thin filaments

                // get interface point...
                double xp[3]; FOR_I3 xp[i] = g[icv_nbr]*n[icv_nbr][i]+x_cv[icv_nbr][i];

                // transform point...
                applyTransform(xp,R,t);
                FOR_I3 xp[i] /= r_vv[icv];

                // weight data...
                //const double wgt = 1.0/MAG(xp);
                //const double wgt = 1.0/DOT_PRODUCT(xp,xp);
                const double wgt = exp(-DOT_PRODUCT(xp,xp));
                //const double wgt = 1.0;

                // calculate terms for least-squares: zp_i = A*xp_i+B*yp_i...
                // Cramers rule terms...
                sxx += wgt*xp[1]*xp[1];
                syy += wgt*xp[2]*xp[2];
                sxy += wgt*xp[1]*xp[2];
                sxz += wgt*xp[1]*xp[0];
                syz += wgt*xp[2]*xp[0];

              }
            }
            const double den = sxx*syy - sxy*sxy;

            if (den != 0.0) {
              const double A = (sxz*syy-syz*sxy)/den;
              const double B = (sxx*syz-sxz*sxy)/den;

              const double mag = sqrt(1.0+A*A+B*B);
              n[icv][0] = 1.0/mag;
              n[icv][1] = -A/mag;
              n[icv][2] = -B/mag;

              // rotate normal back...
              applyInverseRotation(n[icv],R);

            }

          }
        }
      }
      updateCvData(n);

    }

  }


  void checkVofSanity() {

    FOR_ICV {
      assert(vof[icv] == vof[icv]);
      if (vof[icv] < -vof_zero) cout << "undershoots: icv, cv_flag, vof, vof/vof_zero = " << icv << " " << cv_flag[icv] << " " << vof[icv] << " " <<   vof[icv]/vof_zero << endl;
      if (vof[icv] > 1.0+vof_zero) cout << "overshoots: icv, cv_flag, vof,  (vof-1)/vof_zero = " << icv << " " << cv_flag[icv] << " " << vof[icv] << " " <<   (vof[icv]-1.0)/vof_zero << endl;
    }
  }

  void cleanupVof() {
    // clean up vof residue ...
    FOR_ICV {
      if (vof[icv] < vof_zero) vof[icv] = 0.0;
      if (vof[icv] > 1.0 - vof_zero) vof[icv] = 1.0;
    }
  }


  void limitVof() {
    // clean up vof residue ...
    FOR_ICV {
      if (vof[icv] <= vof_zero) vof[icv] = 0.0;
      if (vof[icv] >= 1.0 - vof_zero) vof[icv] = 1.0;
    }
  }


  void updateBaseState(){

    FOR_ICV_G rho[icv] = rho_g*(1.0-vof[icv]) + rho_l*vof[icv];
    FOR_ICV_G mu_lam[icv] = mu_g*(1.0-vof[icv]) + mu_l*vof[icv];
  }

  void calcRhs() {

    double * rhs_vof = new double[ncv];
    double (*rhs_rhou)[3] = new double[ncv][3];
    double * A = new double[cvocv_i[ncv]];

    // calc mass flux, momentum flux w/ part of viscous and full explicit viscous transpose...
    // rreset LHS A matrix
    for (int coc = 0; coc < cvocv_i[ncv]; ++coc) A[coc] = 0.0;

    FOR_ICV {
      rhs_vof[icv] = 0.0;
      FOR_I3 rhs_rhou[icv][i] = 0.0;
    }

    // compute vof Rhs without/with flux limiter ...
    if (step%check_interval == 0 ) COUT1("solve vof");
    bool flux_limiter = false;
    calcVofandMomRhs(rhs_vof, rhs_rhou, flux_limiter);
    flux_limiter = true;
    calcVofandMomRhs(rhs_vof, rhs_rhou, flux_limiter);
    FOR_ICV {
      vof[icv] = vof[icv] + rhs_vof[icv]*dt*inv_vol[icv];
    }

    //checkVofSanity();
    limitVof();
    updateCvData(vof);

    timer.split("solve vof");

    //add time term...
    FOR_ICV FOR_I3 {
      rhs_rhou[icv][i] += vol_cv[icv]*rho[icv]*u[icv][i]/dt;
    }

    FOR_BCZONE (*it)->addMomentumFlux(A,rhs_rhou);

    //update rho and viscosity
    //updateBaseState();
    // viscousity at t^n+1 or n?
    calcViscRhs(rhs_rhou);

    addLspMomentumSource(rhs_rhou);

    //update rho and viscosity
    updateBaseState();
    buildLhs(A);

    timer.split("build rhou lhs");

    // predict u...
    if (step%check_interval == 0 ) COUT1("solve u");
    if (getBoolParam("SOLVER",true)) solveCvJacobi(u,A,rhs_rhou,u_zero,u_relax,u_maxiter,false);

    delete[] rhs_vof;
    delete[] A;
    delete[] rhs_rhou;

    if (step%check_interval == 0 ) dumpRange(u,ncv,"u before correction");

    timer.split("solve u");

  }

  void calcVofandMomRhs(double * rhs_vof, double (*rhs_rhou)[3], bool flux_limiter) {

    double *cfl = NULL;
    cfl  = new double[ncv_g];
    double *vof_out = NULL;
    vof_out = new double[ncv_g];

    const int subgrid_band = 10;

    FOR_ICV {
      vof_out[icv] = 0.0;
    }

    if (flux_limiter) {
      // compute out CFL number ...is it necessary? ....
      calcCfl(cfl,dt);
      updateCvData(cfl);
      FOR_ICV vof_out[icv] = rhs_vof[icv]*dt/vol_cv[icv];
      updateCvData(vof_out);
      //dumpRange(vof_out,ncv,"vof_out");
      FOR_ICV rhs_vof[icv] = 0.0; // reset rhs_vof and recompute with limiters.
    }

    FOR_BCZONE  (*it)->addVofFlux(rhs_vof);
    if (!flux_limiter) FOR_ICV rhs_vof[icv] = fabs(min(rhs_vof[icv],0.0));

    // calc internal flux ...
    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      const double q_mid = 1.25*q_fa[ifa]-0.25*q0_fa[ifa];

      double vof_flux = 0.0;
      double rhou_flux[3] = {0.0,0.0,0.0};
      double vof_c = 0.5*(vof[icv0]+vof[icv1]);
      double gamma; // blending coeffficient gamma = 0: full upwind, gamma =1 : central

      const double interf_fac =fabs( (vof[icv0] + vof[icv1])*(2.0-vof[icv0] - vof[icv1]));
      const double interf_eps = 1.0E-5;
      //if (cv_flag[icv0] >= 1 || cv_flag[icv1] >= 1) {
      //if (interf_fac  > interf_eps) {
      if ( (vof[icv0] > vof_zero && 1.0-vof[icv0] > vof_zero) || (vof[icv1] > vof_zero && 1.0-vof[icv1] > vof_zero) ) {

        double dx0[3], dx1[3],dx[3];
        FOR_I3 dx0[i] = x_fa[ifa][i] - x_cv[icv0][i];
        FOR_I3 dx1[i] = x_fa[ifa][i] - x_cv[icv1][i];
        FOR_I3 dx[i] = dx0[i] - dx1[i];

        double normal0 = -DOT_PRODUCT(dx0,n[icv0]);
        double normal1 = -DOT_PRODUCT(dx1,n[icv1]);


        const double vof0 =  0.5*(1.0+tanh(beta[icv0]*(normal0+g[icv0])));
        const double vof1 =  0.5*(1.0+tanh(beta[icv1]*(normal1+g[icv1])));
        double vof_f = 0.0;

        //computing blending coefficients..
        double dvofdx0 = 0.5*(1.0-pow(tanh(beta[icv0]*(normal0+g[icv0])),2))*beta[icv0];
        double dvofdx1 = 0.5*(1.0-pow(tanh(beta[icv1]*(normal1+g[icv1])),2))*beta[icv1];
        double gradvof = 0.5*fabs(dvofdx0+dvofdx1);
        gamma = ((gradvof > 1.0E-10) ? fabs(vof[icv1]-vof[icv0])/MAG(dx)/gradvof : 0.0);
        gamma = max(0.0, min(gamma,1.0));
        gamma = pow(1.0 - gamma,4);

        double gamma_mom = 0.0;
        if ( vof[icv0] > 0.1 || vof[icv1] > 0.1) gamma_mom = gamma;
        else gamma_mom = 1.0; // central...

        if (vof[icv0] < 0.01  && vof[icv1] < 0.01) {
            gamma = 0.0; // fully upwind ...
        }

        
        if (q_mid >= 0.0) {
          double vof_f0 = (1.0-gamma)*vof0+gamma*vof_c;
          vof_f = vof_f0;
          if (flux_limiter) {
            if (vof_out[icv0] > 0.0) vof_f = min(vof_f0, vof[icv0]/vof_out[icv0]*vof_f0);
            if (cfl[icv0]-vof_out[icv0]  > 0.0) vof_f = max(vof_f, 1.0-(1.0-vof[icv0])/(cfl[icv0]-vof_out[icv0])*(1.0-vof_f0));
            vof_f =max(0.0,min(vof_f,1.0));
          }
        }
        else {
          double vof_f0 = (1.0-gamma)*vof1+gamma*vof_c;
          vof_f = vof_f0;
          if (flux_limiter) {
            if (vof_out[icv1] > 0.0) vof_f = min(vof_f0, vof[icv1]/vof_out[icv1]*vof_f0);
            if (cfl[icv1]-vof_out[icv1]  > 0.0) vof_f = max(vof_f, 1.0-(1.0-vof[icv1])/(cfl[icv1]-vof_out[icv1])*(1.0-vof_f0));
            vof_f =max(0.0,min(vof_f,1.0));
          }
        }

        vof_flux = vof_f*q_mid;
        // compute momentum flux ...
        double rho_flux = rho_g*q_mid + (rho_l-rho_g)*vof_flux;
        FOR_I3 rhou_flux[i]  = 0.5*(1.0-gamma_mom)*(rho_flux*(u[icv0][i]+u[icv1][i])-fabs(rho_flux)*(u[icv1][i]-u[icv0][i]));  // upwind part of momentum flux ..
        FOR_I3 rhou_flux[i] += 0.5*gamma_mom*rho_flux*(u[icv0][i]+u[icv1][i]); // central part of momentum flux ...
      }
      else {
        double vof_f0 = vof_c; // just for symmetry...
        double vof_f = vof_f0;
        if (q_mid >= 0.0) {
          if (flux_limiter) {
            if (vof_out[icv0] > 0.0) vof_f = min(vof_f0, vof[icv0]/vof_out[icv0]*vof_f0);
            if (cfl[icv0]-vof_out[icv0]  > 0.0) vof_f = max(vof_f, 1.0-(1.0-vof[icv0])/(cfl[icv0]-vof_out[icv0])*(1.0-vof_f0));
            vof_f =max(0.0,min(vof_f,1.0));
          }
        }
        else {
          if (flux_limiter) {
            if (vof_out[icv1] > 0.0) vof_f = min(vof_f0, vof[icv1]/vof_out[icv1]*vof_f0);
            if (cfl[icv1]-vof_out[icv1]  > 0.0) vof_f = max(vof_f, 1.0-(1.0-vof[icv1])/(cfl[icv1]-vof_out[icv1])*(1.0-vof_f0));
            vof_f =max(0.0,min(vof_f,1.0));
          }
        }

        vof_flux = vof_f*q_mid;
        double rho_flux = rho_g*q_mid + (rho_l-rho_g)*vof_flux;
        FOR_I3 rhou_flux[i] = 0.25*rho_flux*(u[icv0][i]+u[icv1][i]); // central part of momentum flux ...
      }


      if (flux_limiter) {
        rhs_vof[icv0]            -= vof_flux;
        FOR_I3 rhs_rhou[icv0][i] -= rhou_flux[i];
        if (icv1 < ncv) {
          rhs_vof[icv1]            += vof_flux;
          FOR_I3 rhs_rhou[icv1][i] += rhou_flux[i];
        }
      }
      else { // only compute out-going flux ...
        rhs_vof[icv0]            += max(vof_flux,0.0);
        if (icv1 < ncv)
          rhs_vof[icv1]          += fabs(min(vof_flux,0.0));
      }
    }


   delete[] cfl;
   delete[] vof_out;

  }


  void calcViscRhs2(double (*rhs_rhou)[3]) {

    // calc transpose of viscous terms...
    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      double u_fa[3]; FOR_I3 u_fa[i] = 0.5*(u[icv0][i] + u[icv1][i]);
      const double area_fa = MAG(n_fa[ifa]); assert(area_fa > 0.0);
      const double inv_area_fa = 1.0/area_fa;
      double unit_n_fa[3]; FOR_I3 unit_n_fa[i] = n_fa[ifa][i]*inv_area_fa;
      const double mu_total_fa = MU_FA(mu_lam[icv0]+mu_sgs[icv0],mu_lam[icv1]+mu_sgs[icv1]);
      const double mu_coeff2 = mu_total_fa*area_fa;
      const double mu_coeff = mu_total_fa*area_over_delta_fa[ifa];


      double v_ptr[3] = {0,0,0};

      const double area                   =   MAG(n_fa[ifa]);
      const double aod                    =   area_over_delta_fa[ifa];
      const double mu_tmp                 =   mu_total_fa*area;
      const double mu_total_c             =   mu_total_fa*aod;

      double u0_cn = 0.0;
      double u1_cn = 0.0;
      double unit_n[3];
      FOR_I3 unit_n[i] = n_fa[ifa][i]/area;
      FOR_I3 {
        u0_cn += u[icv0][i]*unit_n[i];
        u1_cn += u[icv1][i]*unit_n[i];
      }
      const double one_third_dun = (u1_cn - u0_cn)/3.0;

      FOR_I3 v_ptr[i] -= mu_total_c*(one_third_dun*unit_n[i]);

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

      FOR_I3 v_ptr[i] -= (dundx[i] - unit_n[i]*(one_third_dundxn + two_third_Skk))*mu_tmp;

      FOR_I3 rhs_rhou[icv0][i] -= v_ptr[i];
      if (icv1 < ncv)
        FOR_I3 rhs_rhou[icv1][i] += v_ptr[i];
    }

  }

  void calcViscRhs(double (*rhs_rhou)[3]) {

    FOR_IFA {
      // viscous terms...
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      double u_fa[3]; FOR_I3 u_fa[i] = 0.5*(u[icv0][i] + u[icv1][i]);
      const double mu_total_fa = MU_FA(mu_lam[icv0]+mu_sgs[icv0],mu_lam[icv1]+mu_sgs[icv1]);
      const double area_fa = MAG(n_fa[ifa]); assert(area_fa > 0.0);
      const double inv_area_fa = 1.0/area_fa;
      double unit_n_fa[3]; FOR_I3 unit_n_fa[i] = n_fa[ifa][i]*inv_area_fa;
      double rhou_flux[3] = {0.0,0.0,0.0};

      double u0_cn = 0.0;
      double u1_cn = 0.0;
      FOR_I3 {
        u0_cn += u[icv0][i]*unit_n_fa[i];
        u1_cn += u[icv1][i]*unit_n_fa[i];
      }
      const double dun = (u1_cn - u0_cn);

      FOR_I3 rhou_flux[i] -= mu_total_fa*area_over_delta_fa[ifa]*(0.5*(u[icv1][i] - u[icv0][i]) + dun*unit_n_fa[i]); // other half handled implicitly

      double dundx[3]           = {0.0, 0.0, 0.0};
      double dundxn   = 0.0;
      FOR_K3 {
        dundx[k] = 0.5*( (dudx[icv0][0][k] + dudx[icv1][0][k])*unit_n_fa[0] +
                         (dudx[icv0][1][k] + dudx[icv1][1][k])*unit_n_fa[1] +
                         (dudx[icv0][2][k] + dudx[icv1][2][k])*unit_n_fa[2] );

        dundxn += dundx[k]*unit_n_fa[k];
      }

      FOR_I3 rhou_flux[i] -= (dundx[i] - unit_n_fa[i]*dundxn)*mu_total_fa*area_fa;

      FOR_I3 rhs_rhou[icv0][i] -= rhou_flux[i];
      if (icv1 < ncv)
        FOR_I3 rhs_rhou[icv1][i] += rhou_flux[i];
    }
/*
    // calc transpose of viscous terms...
    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      double u_fa[3]; FOR_I3 u_fa[i] = 0.5*(u[icv0][i] + u[icv1][i]);
      const double area_fa = MAG(n_fa[ifa]); assert(area_fa > 0.0);
      const double inv_area_fa = 1.0/area_fa;
      double unit_n_fa[3]; FOR_I3 unit_n_fa[i] = n_fa[ifa][i]*inv_area_fa;
      const double mu_total_fa = MU_FA(mu_lam[icv0]+mu_sgs[icv0],mu_lam[icv1]+mu_sgs[icv1]);
      const double mu_coeff2 = mu_total_fa*area_fa;
      const double mu_coeff = mu_total_fa*area_over_delta_fa[ifa];

      double u0_cn = 0.0;
      double u1_cn = 0.0;

      FOR_I3 {
        u0_cn += u[icv0][i]*unit_n_fa[i];
        u1_cn += u[icv1][i]*unit_n_fa[i];
      }
      const double dun = (u1_cn - u0_cn);

      double visc_tr_flux[3] = {0,0,0};
      FOR_I3 visc_tr_flux[i] -= mu_coeff*(dun*unit_n_fa[i]);

      double dundx[3]           = {0.0, 0.0, 0.0};
      double dundxn   = 0.0;
      FOR_K3 {
        dundx[k] = 0.5*( (dudx[icv0][0][k] + dudx[icv1][0][k])*unit_n_fa[0] +
            (dudx[icv0][1][k] + dudx[icv1][1][k])*unit_n_fa[1] +
            (dudx[icv0][2][k] + dudx[icv1][2][k])*unit_n_fa[2] );

        dundxn += dundx[k]*unit_n_fa[k];
      }

      FOR_I3 visc_tr_flux[i] -= mu_coeff2*(dundx[i] - unit_n_fa[i]*dundxn);

      FOR_I3 rhs_rhou[icv0][i] -= visc_tr_flux[i];
      if (icv1 < ncv)
        FOR_I3 rhs_rhou[icv1][i] += visc_tr_flux[i];
    }
*/
  }

  void addLspMomentumSource(double (*rhs_rhou)[3]) {
    // add momentum source from lsp
    if (lsp) {
      double my_factor_min = 0.0;
      FOR_ICV {
        // potentially limit the lspMomentumSource...
        double lspSourceMag = sqrt(DOT_PRODUCT(lspMomentumSource[icv],lspMomentumSource[icv]));
        // the folowing should have a velocity scale as well, but for now assume 20.0. So this is
        // saying that the momentum source term from the particles will not be able to change the
        // velocity by more than 20 or so...
        double this_factor = min(1.0,20.0*(rho[icv])/(dt*lspSourceMag));
        my_factor_min = min(my_factor_min,this_factor);
        FOR_I3 rhs_rhou[icv][i] += this_factor*lspMomentumSource[icv][i]*vol_cv[icv];
      }

      if (step%check_interval == 0) {
        double factor_min;
        MPI_Reduce(&my_factor_min,&factor_min,1,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
        if (mpi_rank == 0)
          cout << " > minimum particle limiting factor: " << factor_min << endl;
      }
    }
  }

  double calcCfl(double* cfl_, const double dt_) const {

    // functions returns the rank-local cfl maximum..

    bool b_memflag = false;
    double * cfl   = NULL ;
    if ( cfl_ == NULL) {
      cfl       = new double[ncv];
      b_memflag = true;
    } else {
      cfl = cfl_;
    }


    FOR_ICV cfl[icv] = 0.0;

    for (vector<VofBc*>::const_iterator it = bcs.begin(); it != bcs.end(); ++it) {
      if ((*it)->q_bf != NULL) {
        for (int ibf = 0; ibf < (*it)->zone_ptr->nbf; ++ibf) {
          const int icv0 = (*it)->zone_ptr->cvobf[ibf];
          cfl[icv0] += fabs(min((*it)->q_bf[ibf],0.0)); // 0.5*dt/vol below...
        }
      }
    }

    for (int ifa = 0; ifa<nfa; ++ifa){
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      const double qmid = 1.25*q_fa[ifa] - 0.25*q0_fa[ifa];
      cfl[icv0] += max(qmid,0.0);
      if (icv1 < ncv)
        cfl[icv1] += fabs(min(qmid,0.0));
    }

    FOR_ICV cfl[icv] *= dt_/vol_cv[icv];

    double my_cfl_max = 0.0;
    FOR_ICV my_cfl_max = max(my_cfl_max,cfl[icv]);

    //assert(my_cfl_max <= 1.0); // necessary CFL condition for vof advection ....

    if ( b_memflag) delete[] cfl;
    return my_cfl_max;
  }


  void dumpCfl() const {
    double * cfl = new double[ncv];
    calcCfl(cfl,dt);
    MiscUtils::dumpRange(cfl,ncv,"CFL");
    delete[] cfl;
  }

  void linkCvs(const int icv0,const int icv1,int * prev,int * next) {

    assert(icv0 != icv1);

    // check the start of icv0's list for icv1...
    int icv0_ = icv0;
    while (prev[icv0_] != -1) {
      icv0_ = prev[icv0_];
      if (icv0_ == icv1)
        return;
    }
    // then check the end...
    icv0_ = icv0;
    while (next[icv0_] != -1) {
      icv0_ = next[icv0_];
      if (icv0_ == icv1)
        return;
    }
    // then do the opposite for ino1...
    int icv1_ = icv1;
    while (next[icv1_] != -1) {
      icv1_ = next[icv1_];
      if (icv1_ == icv0)
        return;
    }
    // the start..
    icv1_ = icv1;
    while (prev[icv1_] != -1) {
      icv1_ = prev[icv1_];
      if (icv1_ == icv0)
        return;
    }
    // now join the end of ino0_ to the start of ino1_...
    next[icv0_] = icv1_;
    prev[icv1_] = icv0_;

  }


  void solvePAndCorrectU() {

    if (step%check_interval == 0 ) COUT1("solvePAndCorrectU()...");

    // build the div using uncorrected u...

    FOR_ICV div[icv] = 0.0;

    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      const double sp_vol_fa = SP_VOL_FA(rho[icv0],rho[icv1]);
      const double area_fa = MAG(n_fa[ifa]); assert(area_fa > 0.0);
      double kappa_fa = 0.0;
      if (cv_flag[icv0] >= 1 || cv_flag[icv1] >= 1) {
        const double wt0 = pow(vof[icv0]*(1.0-vof[icv0]),1);
        const double wt1 = pow(vof[icv1]*(1.0-vof[icv1]),1);
        const double wt_sum = wt0 + wt1;
        if (wt_sum > 0.0) {
          kappa_fa = (wt0*kappa[icv0]+wt1*kappa[icv1])/wt_sum;
        }
        else {
          kappa_fa = 0.5*(kappa[icv0]+kappa[icv1]);
        }
      }
      const double sp_tension_A_fa = sp_vol_fa*sigma*kappa_fa*(vof[icv1]-vof[icv0])*area_over_delta_fa[ifa];
      //double sp_tension_A_fa;
      //const double vof_norm = fabs(vof[icv1]-vof[icv0]);
      //if ( vof_norm != 0.0 )  sp_tension_A_fa = sp_vol_fa*sigma*kappa_fa*(vof[icv1]-vof[icv0])/vof_norm*area_over_delta_fa[ifa];
      //else sp_tension_A_fa = sp_vol_fa*sigma*kappa_fa*(vof[icv1]-vof[icv0])*area_over_delta_fa[ifa];

      double u_fa[3]; FOR_I3 u_fa[i] = 0.5*(u[icv0][i]+u[icv1][i]);
      const double un_A = DOT_PRODUCT(u_fa,n_fa[ifa]) + dt*sp_tension_A_fa;
      div[icv0] += un_A;
      if (icv1 < ncv)
        div[icv1] -= un_A;
    }

    FOR_BCZONE {
      Param * param = getParam((*it)->getName());
      if ( param->getString(0) != "OUTLET")
        (*it)->addFlux(div);
    }
    FOR_BCZONE (*it)->updateBc();
    FOR_BCZONE {
      Param * param = getParam((*it)->getName());
      if ( param->getString(0) == "OUTLET")
        (*it)->addFlux(div);
    }


    // check rhs...
    if (step%check_interval == 0 ) {
      dumpRange(div,ncv,"div before correction");
      double my_sum = 0.0;
      FOR_ICV my_sum += div[icv];
      double sum;
      MPI_Reduce(&my_sum,&sum,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      if (mpi_rank == 0)
        cout << " > sum(rhs) (should be zero): " << sum << endl;
    }

    //FOR_ICV div[icv] /= dt;
    FOR_ICV div[icv] /= -dt;

    timer.split("calc poisson rhs");

    // build the lhs matrix...

    double * A = new double[cvocv_i[ncv]];
    buildCvDivSpVolGrad(A);

    timer.split("build poisson lhs");

    // solve p...
    switch (pressure_solver) {
      case TRILINOS_SOLVER:
        {

          // jacobi preconditioning: Cl_inv*A*Cr_inv*Cr*x=Cl_inv*b
          double *inv_sqrt_diag = NULL;
          bool jacobi_prec = checkParam("JACOBI_PREC");
          //double p_tol = p_zero;
          if (jacobi_prec) {
            inv_sqrt_diag = new double[ncv_g];
            FOR_ICV {
              assert(A[cvocv_i[icv]] < 0.0);
              inv_sqrt_diag[icv] = 1.0/sqrt(-A[cvocv_i[icv]]);
            }
            updateCvData(inv_sqrt_diag);
            FOR_ICV {
              // A' = Cl_inv*A*Cr_ing...
              A[cvocv_i[icv]] = -1.0;
              for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
                const int icv_nbr = cvocv_v[coc];
                A[coc] *= inv_sqrt_diag[icv]*inv_sqrt_diag[icv_nbr];
              }

              // b' = Cl_inv*b...
              div[icv] *= inv_sqrt_diag[icv];

              // x0' = Cr*x0 (better initial guess)...
              p[icv] /= inv_sqrt_diag[icv];
            }
          }

          // hypre and trilinos requires full reduced matrix and cell information ..
          //double * A_full;
          //int * cvocv_full_i;
          int * cvocv_full_v;
          int * icv_global;
          buildFullCvoCv(icv_global, cvocv_full_v);
          int niters = -1;

          if (trilinosSolverForPressure == NULL) {
            trilinosSolverForPressure = new TrilinosSolver();
            trilinosSolverForPressure->setup(ncv, icv_global, cvocv_i, cvocv_full_v,mpi_comm);

            MPI_Barrier(mpi_comm);
            niters = trilinosSolverForPressure->solve_first(p,div,A, trilinos_idata[0], // OPTION
                                                                     p_maxiter, // MAXITER
                                                                     trilinos_idata[1], // ML_MAX_LEVELS
                                                                     trilinos_idata[2], // REPART_LEVEL
                                                                     p_zero   ); //TOL
          }
          else {
            trilinosSolverForPressure->setMatrixGlobalValues(A) ; // update the entries in the matrix ..
            if ( rebuildPoissPrec){
              if (mpi_rank == 0 && step%check_interval == 0 ) cout << " > " << " " << "update ML Preconditioner "<< endl;
              trilinosSolverForPressure->updateMLPreconditioner() ;
            }

            // and solve...
            niters = trilinosSolverForPressure->solve_again( p, div, step%check_interval==0) ;

            if (niters < 0 || niters == p_maxiter) { // removed !rebuildPoissPrec because sometimes you need to try again even if you rebuilt
              if (mpi_rank == 0 && step%check_interval == 0 ) cout << " > " << " " << "update ML Preconditioner "<< endl;
              trilinosSolverForPressure->updateMLPreconditioner() ;
              niters = trilinosSolverForPressure->solve_again(p,div, step%check_interval==0);
            }

          }
          if ( niters > p_maxiter/2 )
            rebuildPoissPrec = true ;
          else
            rebuildPoissPrec = false ; // always rebuild..

          if ( checkParam("FORCE_ML_REBUILD"))
            rebuildPoissPrec = true ;

          delete[] cvocv_full_v;
          delete[] icv_global;

          if (jacobi_prec) {
            // x = Cr_inv*x'...
            FOR_ICV p[icv] *= inv_sqrt_diag[icv];
            delete[] inv_sqrt_diag;
          }

          break;
        }
      case BCGSTAB_SOLVER:
        {
          solveCvCg(p,A,div,p_zero,p_maxiter,false); // need to replace eventually...
          break;
        }
      case JACOBI_SOLVER:
        {
          solveCvJacobi(p,A,div,p_zero,p_relax,p_maxiter,false);
          break;
        }
      case AMG_SOLVER:
        {
          if (am == NULL) {
            const double agglomeration_factor = getDoubleParam("AGGLOMERATION_FACTOR",8.0); // 2^3 = 8
            const bool split_orphaned_colors = getBoolParam("SPLIT_ORPHANED_COLORS",true); 
            const double amg_coeff = getDoubleParam("AMG_COEFF",pow(agglomeration_factor,-1.0/3.0)); // 1/2
            const int ncg = getIntParam("NCG",0); 
            am = new Multigrid::AlgebraicMultigrid();
            am->init(vol_cv,x_cv,A,cvocv_i,cvocv_v,icv_global,rbi_g,NULL,ncv,ncv_g,ncv_g,ncv_global,ncg,
                agglomeration_factor,split_orphaned_colors,amg_coeff);
          }
          else {
            am->update(A);
          }
          const double relax = getDoubleParam("RELAX",0.7); // jacobi needs something around this
          const int nsmooth = getIntParam("NSMOOTH",10);
          const int dnsmooth = getIntParam("DNSMOOTH",0);
          const string smoother = getStringParam("SMOOTHER","JACOBI");
          const int maxiter_coarsest = getIntParam("MAXITER_COARSEST",-1);
          const string solver_coarsest = getStringParam("SOLVER_COARSEST",smoother);
          const bool verbose = (step%check_interval == 0);
          const int maxcycle = getIntParam("MAXCYCLE",p_maxiter);
          const string cycle = getStringParam("CYCLE","V_CYCLE");
          const bool b_linear_prolongation = checkParam("LINEAR_PROLONGATION");
          am->solve(p,div,p_zero,cycle,maxcycle,b_linear_prolongation,nsmooth,dnsmooth,smoother,relax,maxiter_coarsest,solver_coarsest,verbose);
        }
    }

    delete[] A;
    updateCvData(p); // FH: should be done in solver, so not necessary?

    timer.split("solve poisson");

    // correct u...

    calcSpDpdxandFlux();
    FOR_ICV FOR_I3 u[icv][i] -= dt*sp_dpdx[icv][i];
    updateCvData(u);
    //dumpRange(u,ncv,"u");

    timer.split("correct u");

    {
      FOR_ICV div[icv] = 0.0;
      FOR_IFA {
        const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
        const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
        div[icv0] += q_fa[ifa];
        if (icv1 < ncv)
          div[icv1] -= q_fa[ifa];
      }
      FOR_BCZONE (*it)->addFlux(div);
      if (step%check_interval == 0) dumpRange(div,ncv,"div after correction");
    }
    timer.split("calc div");
  }


  void buildLhs(double * A) {

    FOR_ICV {
      const int coc_f = cvocv_i[icv];
      A[coc_f] += vol_cv[icv]*rho[icv]/dt;
    }

    // internal faces...
    for (int ifa = 0; ifa < nfa; ++ifa){

      const int icv0 = cvofa[ifa][0];
      assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1];
      assert((icv1 >= 0)&&(icv0 < ncv_g));

      const int coc00 = cvocv_i[icv0];
      int coc01 = coc00+1;
      while (cvocv_v[coc01] != icv1) {
        ++coc01;
        assert(coc01 != cvocv_i[icv0+1] );
      }

      const double mu_coeff =0.5*MU_FA(mu_lam[icv0] + mu_sgs[icv0] , mu_lam[icv1] + mu_sgs[icv1])*area_over_delta_fa[ifa];

      double mf_coeff = 0.0;
      if ((cv_flag[icv0] >= 1)||(cv_flag[icv1] >= 1)) {
        // fully explicit at interface to be consistent with vof
        mf_coeff = 0.0;
      }
      else {
        // ck away the interface where
        // density is going to constant b/w icv0 & icv1 during the step.
        // will average just for symmetry sake (shouldn't matter)...
        mf_coeff = 0.125*(rho[icv0]+rho[icv1])*(1.25*q_fa[ifa]-0.25*q0_fa[ifa]); // other half handled explicitly
        //mf_coeff = 0.25*(rho[icv0]+rho[icv1])*(1.25*q_fa[ifa]-0.25*q0_fa0[ifa]); // fully implicit
        //  mf_coeff = 0.0;  // fully explicit
      }

      // no convection in the operator -- just diffusion...
      A[coc00] += mf_coeff + mu_coeff;
      A[coc01] += mf_coeff - mu_coeff;

      if (icv1 < ncv) {

        const int coc11 = cvocv_i[icv1];
        int coc10 = coc11+1;
        while (cvocv_v[coc10] != icv0) {
          ++coc10;
          assert(coc10 != cvocv_i[icv1+1] );
        }

        A[coc11] -= mf_coeff - mu_coeff;
        A[coc10] -= mf_coeff + mu_coeff;

      }

    }

  }

  void buildCvDivSpVolGrad(double * A) {

    for (int coc = 0; coc < cvocv_i[ncv]; ++coc)
      A[coc] = 0.0;

    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      //const double sp_vol_fa = 1.0/rho_f;
      //const double coeff = sp_vol_fa*area_over_delta_fa[ifa];
      //const double coeff = SP_VOL_FA(rho[icv0],rho[icv1])*area_over_delta_fa[ifa];
      const double coeff = -SP_VOL_FA(rho[icv0],rho[icv1])*area_over_delta_fa[ifa];

      {
	int coc = cvocv_i[icv0];
	assert(cvocv_v[coc] == icv0);
	A[coc] -= coeff;
	++coc;
	while (cvocv_v[coc] != icv1)
	  ++coc;
	assert(coc < cvocv_i[icv0+1]);
	A[coc] += coeff;
      }
      if (icv1 < ncv) {
	int coc = cvocv_i[icv1];
	assert(cvocv_v[coc] == icv1);
	A[coc] -= coeff;
	++coc;
	while (cvocv_v[coc] != icv0)
	  ++coc;
	assert(coc < cvocv_i[icv1+1]);
	A[coc] += coeff;
      }
    }

  }

  void buildFullCvoCv(int * &icv_global, int * &cvocv_full_v) {

    // number of neighbors and active cells are same....
    const int cvocv_s = cvocv_i[ncv];
    cvocv_full_v = new int[cvocv_s];
    icv_global = new int[ncv];

    FOR_ICV icv_global[icv] = getIcvGlobal(icv);

    FOR_ICV {
      for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_local = cvocv_v[coc];
        cvocv_full_v[coc] = getIcvGlobal(icv_local);
      }
    }

  }

  void dropIdentification() {

    if (mpi_rank == 0 && step%check_interval == 0 ) cout << " > starting Drop Identification " << endl;

    double *vofs = new double[ncv_g];
    const int subgrid = 2;
    const double vof_tol = getDoubleParam("VOF_TOL",0.1);
    FOR_ICV_G {
      vofs[icv] = vof[icv];
      if (vofs[icv] < vof_tol) vofs[icv] = 0.0;
    }

 //   FOR_ICV_G {
 //     if (iband[icv] <= subgrid) vofs[icv] = vof[icv];
 //     else vofs[icv] = 0.0;
 //   }
    updateCvData(vofs);

    // grouping algorithm uses cv_flag... copy cv_flag and restore back after drop Transfer done ....
    int *cv_flag_copy = new int[ncv_g];
    FOR_ICV_G cv_flag_copy[icv] = cv_flag[icv];

    int * fa_flag = new int[nfa];
    FOR_IFA fa_flag[ifa] = 0;
    // if twp adjacent cvs connected by an face have vof > threshold, they belong to a group.
    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      //if (iband[icv0] <= subgrid && iband[icv1]  <= subgrid) {
      if (vofs[icv0] >= vof_tol && vofs[icv1]  >= vof_tol) {
        fa_flag[ifa] = 1;
      }
    }

    AverageOperator * a = NULL;
    initCvAverageOperator(a, fa_flag);

    delete[] fa_flag;


    // for each group, we need:
    // 1. sum of liquid volumes,
    // 2. center of mass
    // 3. mass-averaged velocity of the liquid...
    // 4. sum of the volume fraction

    double * group_liquid_volume  = new double[a->ngr];
    double (*group_x)[3] = new double[a->ngr][3];
    double (*group_u)[3] = new double[a->ngr][3];
    double * group_max_vof = new double[a->ngr];

    for (int igr = 0; igr < a->ngr; ++igr) {
      group_liquid_volume[igr] = 0.0;
      FOR_I3 group_x[igr][i] = 0.0;
      FOR_I3 group_u[igr][i] = 0.0;
      group_max_vof[igr] = 0.0;
    }

    FOR_ICV {
      const int igr = a->group[icv];
      group_liquid_volume[igr] += vofs[icv]*vol_cv[icv];
      FOR_I3 group_x[igr][i] += vofs[icv]*vol_cv[icv]*x_cv[icv][i];
      FOR_I3 group_u[igr][i] += vofs[icv]*vol_cv[icv]*u[icv][i];
      if (iband[icv] == 1) group_max_vof[igr] += 1.0;
    }

    a->push(group_liquid_volume,  group_liquid_volume+a->ngr_a,  ADD_DATA);
    a->push(group_x,           group_x+a->ngr_a,           ADD_DATA);
    a->push(group_u,           group_u+a->ngr_a,           ADD_DATA);
    a->push(group_max_vof,       group_max_vof+a->ngr_a,       ADD_DATA);

    // go through the active groups (these are the groups we own with reduced
    // values)  and look at the total volume...
    int * group_flag = new int[a->ngr];
    // buffers for gathering drops to root 
    double * bufDropD = new double[a->ngr_a];
    double (*bufDropX)[3] = new double[a->ngr_a][3];
    double (*bufDropU)[3] = new double[a->ngr_a][3];

    int ndrops = 0;
    for (int igr = 0; igr < a->ngr_a; ++igr) {
      group_flag[igr] = 0;
      //if (group_max_vof[igr] >= 1.0) {
        FOR_I3 group_x[igr][i] /= group_liquid_volume[igr];
        FOR_I3 group_u[igr][i] /= group_liquid_volume[igr];
        bufDropD[ndrops] = pow(6.0*group_liquid_volume[igr]/M_PI,1.0/3.0);
        FOR_I3 {
          bufDropX[ndrops][i] = group_x[igr][i];
          bufDropU[ndrops][i] = group_u[igr][i];
        }
        //cout << "mpi_rank, idrop, diameter = " << mpi_rank << " " << ndrops << " " << bufDropD[ndrops] << endl;
        ndrops++;
        group_flag[igr] = 1;
    //  }
    }

    int ndrops_global = 0;
    MPI_Allreduce(&ndrops, &ndrops_global, 1, MPI_INT, MPI_SUM, mpi_comm);
    if (mpi_rank == 0 && step%check_interval == 0 ) {
      cout << "Total number of resolved drops  = " << " " << ndrops_global << endl;
    }


    updateDropStats(bufDropD, bufDropX, bufDropU, ndrops,time,dt,step);

    delete[] group_liquid_volume;
    delete[] group_x;
    delete[] group_u;
    delete[] group_max_vof;
    delete[] vofs;
    delete a;
    delete[] group_flag;

    // now restore back cv_flag
    FOR_ICV_G cv_flag[icv] = cv_flag_copy[icv];

    delete[] cv_flag_copy;
    delete[] bufDropD;
    delete[] bufDropX;
    delete[] bufDropU;
    
  }


 /*
  // calc inviscid specific forces (pressure and surface tension forces per mass) - units of accelaration
  void calcSpTension() {

    // uniformly weighted - definately non-conservative, but maybe more stable?

    double (*ninjdA_diag)[3] = new double[ncv][3];
    double (*ninjdA_offd)[3] = new double[ncv][3];

    FOR_ICV FOR_I3 ninjdA_diag[icv][i] = 0.0;
    FOR_ICV FOR_I3 ninjdA_offd[icv][i] = 0.0;
    FOR_ICV FOR_I3 sp_tension[icv][i] = 0.0; // use for rhs

    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      const double sp_vol_fa = SP_VOL_FA(rho[icv0],rho[icv1]);
      const double area_fa = MAG(n_fa[ifa]); assert(area_fa > 0.0);
      const double inv_area_fa = 1.0/area_fa;
      double unit_n_fa[3]; FOR_I3 unit_n_fa[i] = n_fa[ifa][i]*inv_area_fa;
      double kappa_fa = 0.0;
      if (cv_flag[icv0] >= 1 && cv_flag[icv1] >= 1) {
        const double wt0 = pow(vof[icv0]*(1.0-vof[icv0]),1);
        const double wt1 = pow(vof[icv1]*(1.0-vof[icv1]),1);
        const double wt_sum = wt0 + wt1;
        if (wt_sum > 0.0) {
          kappa_fa = (wt0*kappa[icv0]+wt1*kappa[icv1])/wt_sum;
        }
        else {
          kappa_fa = 0.5*(kappa[icv0]+kappa[icv1]);
        }
      }
      else if (cv_flag[icv0] >= 1) kappa_fa = kappa[icv0];
      else if (cv_flag[icv1] >= 1) kappa_fa = kappa[icv1];
      const double sp_tension_n_fa = sp_vol_fa*sigma*kappa_fa*(vof[icv1]-vof[icv0])*area_over_delta_fa[ifa]*inv_area_fa;

      FOR_I3 ninjdA_diag[icv0][i] += unit_n_fa[i]*n_fa[ifa][i];
      ninjdA_offd[icv0][0] += unit_n_fa[1]*n_fa[ifa][2];
      ninjdA_offd[icv0][1] += unit_n_fa[2]*n_fa[ifa][0];
      ninjdA_offd[icv0][2] += unit_n_fa[0]*n_fa[ifa][1];
      FOR_I3 sp_tension[icv0][i] += sp_tension_n_fa*n_fa[ifa][i];

      if (icv1 < ncv) {
        FOR_I3 ninjdA_diag[icv1][i] += unit_n_fa[i]*n_fa[ifa][i];
        ninjdA_offd[icv1][0] += unit_n_fa[1]*n_fa[ifa][2];
        ninjdA_offd[icv1][1] += unit_n_fa[2]*n_fa[ifa][0];
        ninjdA_offd[icv1][2] += unit_n_fa[0]*n_fa[ifa][1];
        FOR_I3 sp_tension[icv1][i] += sp_tension_n_fa*n_fa[ifa][i];
      }
    }

    FOR_IBF {
      const double area_bf = MAG(n_bf[ibf]);
      if (area_bf > 0.0) {
	const double unit_n[3] = {
	  n_bf[ibf][0]/area_bf,
	  n_bf[ibf][1]/area_bf,
	  n_bf[ibf][2]/area_bf
	};
	const int icv = cvobf[ibf];
	FOR_I3 ninjdA_diag[icv][i] += area_bf*unit_n[i]*unit_n[i];
	FOR_I3 ninjdA_offd[icv][i] += area_bf*unit_n[(i+1)%3]*unit_n[(i+2)%3];
	// assume dpdn == 0, so nothing to add...
      }
    }

    FOR_ICV {

      double inv_denom = 1.0/(     ninjdA_diag[icv][0]*ninjdA_diag[icv][1]*ninjdA_diag[icv][2] +
			       2.0*ninjdA_offd[icv][0]*ninjdA_offd[icv][1]*ninjdA_offd[icv][2] -
			           ninjdA_diag[icv][0]*ninjdA_offd[icv][0]*ninjdA_offd[icv][0] -
			           ninjdA_diag[icv][1]*ninjdA_offd[icv][1]*ninjdA_offd[icv][1] -
			           ninjdA_diag[icv][2]*ninjdA_offd[icv][2]*ninjdA_offd[icv][2] );

      const double rhs[3] = { sp_tension[icv][0], sp_tension[icv][1], sp_tension[icv][2] };

      sp_tension[icv][0] = inv_denom*( (ninjdA_diag[icv][1]*ninjdA_diag[icv][2]-ninjdA_offd[icv][0]*ninjdA_offd[icv][0])*rhs[0] +
                                       (ninjdA_offd[icv][0]*ninjdA_offd[icv][1]-ninjdA_diag[icv][2]*ninjdA_offd[icv][2])*rhs[1] +
                                       (ninjdA_offd[icv][2]*ninjdA_offd[icv][0]-ninjdA_diag[icv][1]*ninjdA_offd[icv][1])*rhs[2] );

      sp_tension[icv][1] = inv_denom*( (ninjdA_offd[icv][0]*ninjdA_offd[icv][1]-ninjdA_diag[icv][2]*ninjdA_offd[icv][2])*rhs[0] +
                                       (ninjdA_diag[icv][2]*ninjdA_diag[icv][0]-ninjdA_offd[icv][1]*ninjdA_offd[icv][1])*rhs[1] +
                                       (ninjdA_offd[icv][1]*ninjdA_offd[icv][2]-ninjdA_diag[icv][0]*ninjdA_offd[icv][0])*rhs[2] );

      sp_tension[icv][2] = inv_denom*( (ninjdA_offd[icv][2]*ninjdA_offd[icv][0]-ninjdA_diag[icv][1]*ninjdA_offd[icv][1])*rhs[0] +
                                       (ninjdA_offd[icv][1]*ninjdA_offd[icv][2]-ninjdA_diag[icv][0]*ninjdA_offd[icv][0])*rhs[1] +
                                       (ninjdA_diag[icv][0]*ninjdA_diag[icv][1]-ninjdA_offd[icv][2]*ninjdA_offd[icv][2])*rhs[2] );

    }

    delete[] ninjdA_diag;
    delete[] ninjdA_offd;

    updateCvData(sp_tension);
    dumpRange(sp_tension,ncv,"sp_tension");
  }
*/
  // calc inviscid specific forces (pressure and surface tension forces per mass) - units of accelaration
  void calcSpDpdxandFlux() {

    // store previous for AB extrapolation...
    FOR_IFA q0_fa[ifa] = q_fa[ifa];

    // uniformly weighted - definately non-conservative, but maybe more stable?

    double (*ninjdA_diag)[3] = new double[ncv][3];
    double (*ninjdA_offd)[3] = new double[ncv][3];

    FOR_ICV FOR_I3 ninjdA_diag[icv][i] = 0.0;
    FOR_ICV FOR_I3 ninjdA_offd[icv][i] = 0.0;
    FOR_ICV FOR_I3 sp_dpdx[icv][i] = 0.0;

    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      //const double sp_vol_fa = 1.0/rho_f;
      const double sp_vol_fa = SP_VOL_FA(rho[icv0],rho[icv1]);
      const double area_fa = MAG(n_fa[ifa]); assert(area_fa > 0.0);
      const double inv_area_fa = 1.0/area_fa;
      double unit_n_fa[3]; FOR_I3 unit_n_fa[i] = n_fa[ifa][i]*inv_area_fa;
      double kappa_fa = 0.0;
      if (cv_flag[icv0] >= 1 || cv_flag[icv1] >= 1) {
        const double wt0 = pow(vof[icv0]*(1.0-vof[icv0]),1);
        const double wt1 = pow(vof[icv1]*(1.0-vof[icv1]),1);
        const double wt_sum = wt0 + wt1;
        if (wt_sum > 0.0) {
          kappa_fa = (wt0*kappa[icv0]+wt1*kappa[icv1])/wt_sum;
        }
        else {
          kappa_fa = 0.5*(kappa[icv0]+kappa[icv1]);
        }
      }

      const double sp_force_inv_A_fa = -sp_vol_fa*(p[icv1]-p[icv0]-sigma*kappa_fa*(vof[icv1]-vof[icv0]))*area_over_delta_fa[ifa];
      //double sp_force_inv_A_fa;
      //const double vof_norm = fabs(vof[icv1]-vof[icv0]);
      //if (vof_norm != 0.0) sp_force_inv_A_fa = -sp_vol_fa*(p[icv1]-p[icv0]-sigma*kappa_fa*(vof[icv1]-vof[icv0])/vof_norm)*area_over_delta_fa[ifa];
      //else sp_force_inv_A_fa = -sp_vol_fa*(p[icv1]-p[icv0]-sigma*kappa_fa*(vof[icv1]-vof[icv0]))*area_over_delta_fa[ifa];
      const double sp_dpdx_n_fa      = -sp_force_inv_A_fa*inv_area_fa;
      double u_fa[3]; FOR_I3 u_fa[i] = 0.5*(u[icv0][i]+u[icv1][i]);

      // conservative volume flux...
      q_fa[ifa] = DOT_PRODUCT(u_fa,n_fa[ifa]) + dt*sp_force_inv_A_fa;

      FOR_I3 ninjdA_diag[icv0][i] += unit_n_fa[i]*n_fa[ifa][i];
      ninjdA_offd[icv0][0] += unit_n_fa[1]*n_fa[ifa][2];
      ninjdA_offd[icv0][1] += unit_n_fa[2]*n_fa[ifa][0];
      ninjdA_offd[icv0][2] += unit_n_fa[0]*n_fa[ifa][1];
      FOR_I3 sp_dpdx[icv0][i] += sp_dpdx_n_fa*n_fa[ifa][i];

      if (icv1 < ncv) {
        FOR_I3 ninjdA_diag[icv1][i] += unit_n_fa[i]*n_fa[ifa][i];
        ninjdA_offd[icv1][0] += unit_n_fa[1]*n_fa[ifa][2];
        ninjdA_offd[icv1][1] += unit_n_fa[2]*n_fa[ifa][0];
        ninjdA_offd[icv1][2] += unit_n_fa[0]*n_fa[ifa][1];
        FOR_I3 sp_dpdx[icv1][i] += sp_dpdx_n_fa*n_fa[ifa][i];
      }
    }

    FOR_IBF {
      const double area_bf = MAG(n_bf[ibf]);
      if (area_bf > 0.0) {
	const double unit_n[3] = {
	  n_bf[ibf][0]/area_bf,
	  n_bf[ibf][1]/area_bf,
	  n_bf[ibf][2]/area_bf
	};
	const int icv = cvobf[ibf];
	FOR_I3 ninjdA_diag[icv][i] += area_bf*unit_n[i]*unit_n[i];
	FOR_I3 ninjdA_offd[icv][i] += area_bf*unit_n[(i+1)%3]*unit_n[(i+2)%3];
	// assume dpdn == 0, so nothing to add...
      }
    }


    FOR_ICV {

      double inv_denom = 1.0/(     ninjdA_diag[icv][0]*ninjdA_diag[icv][1]*ninjdA_diag[icv][2] +
			       2.0*ninjdA_offd[icv][0]*ninjdA_offd[icv][1]*ninjdA_offd[icv][2] -
			           ninjdA_diag[icv][0]*ninjdA_offd[icv][0]*ninjdA_offd[icv][0] -
			           ninjdA_diag[icv][1]*ninjdA_offd[icv][1]*ninjdA_offd[icv][1] -
			           ninjdA_diag[icv][2]*ninjdA_offd[icv][2]*ninjdA_offd[icv][2] );

      const double rhs[3] = { sp_dpdx[icv][0], sp_dpdx[icv][1], sp_dpdx[icv][2] };

      sp_dpdx[icv][0] = inv_denom*( (ninjdA_diag[icv][1]*ninjdA_diag[icv][2]-ninjdA_offd[icv][0]*ninjdA_offd[icv][0])*rhs[0] +
                                    (ninjdA_offd[icv][0]*ninjdA_offd[icv][1]-ninjdA_diag[icv][2]*ninjdA_offd[icv][2])*rhs[1] +
                                    (ninjdA_offd[icv][2]*ninjdA_offd[icv][0]-ninjdA_diag[icv][1]*ninjdA_offd[icv][1])*rhs[2] );

      sp_dpdx[icv][1] = inv_denom*( (ninjdA_offd[icv][0]*ninjdA_offd[icv][1]-ninjdA_diag[icv][2]*ninjdA_offd[icv][2])*rhs[0] +
                                    (ninjdA_diag[icv][2]*ninjdA_diag[icv][0]-ninjdA_offd[icv][1]*ninjdA_offd[icv][1])*rhs[1] +
                                    (ninjdA_offd[icv][1]*ninjdA_offd[icv][2]-ninjdA_diag[icv][0]*ninjdA_offd[icv][0])*rhs[2] );

      sp_dpdx[icv][2] = inv_denom*( (ninjdA_offd[icv][2]*ninjdA_offd[icv][0]-ninjdA_diag[icv][1]*ninjdA_offd[icv][1])*rhs[0] +
                                    (ninjdA_offd[icv][1]*ninjdA_offd[icv][2]-ninjdA_diag[icv][0]*ninjdA_offd[icv][0])*rhs[1] +
                                    (ninjdA_diag[icv][0]*ninjdA_diag[icv][1]-ninjdA_offd[icv][2]*ninjdA_offd[icv][2])*rhs[2] );

    }

    delete[] ninjdA_diag;
    delete[] ninjdA_offd;

    //dumpRange(sp_dpdx,ncv,"sp_dpdx");
    //dumpRange(q_fa,nfa,"q_fa");
  }

  void lspTransfer() {

    if (mpi_rank == 0 && step%check_interval == 0 ) cout << " > staring LSP TRANSFER " << endl;
    dumpRange(cv_flag_real,ncv,"BAND");

    double *vofs = new double[ncv_g];
    const int subgrid = 10;
    const int nband = 20;

    FOR_ICV_G {
      if (iband[icv] > subgrid) vofs[icv] = vof[icv];
      else vofs[icv] = 0.0;
    }
    updateCvData(vofs);

    // grouping algorithm uses cv_flag... copy cv_flag and restore back after drop Transfer done ....
    int *cv_flag_copy = new int[ncv_g];
    FOR_ICV_G cv_flag_copy[icv] = cv_flag[icv];

    int * fa_flag = new int[nfa];
    FOR_IFA fa_flag[ifa] = 0;
    // if twp adjacent cvs connected by an face have vof > threshold, they belong to a group.
    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      //if (vofs[icv0] > vof_zero && vofs[icv1] > vof_zero) {
      if (iband[icv0] > subgrid && iband[icv1]  > subgrid) {
        fa_flag[ifa] = 1;
      }
    }

    AverageOperator * a = NULL;
    initCvAverageOperator(a, fa_flag);

    delete[] fa_flag;


    // for each group, we need:
    // 1. sum of liquid volumes,
    // 2. center of mass
    // 3. mass-averaged velocity of the liquid...
    // 4. sum of the volume fraction

    double * group_liquid_volume  = new double[a->ngr];
    double * group_liquid_surface = new double[a->ngr];
    double (*group_x)[3] = new double[a->ngr][3];
    double (*group_u)[3] = new double[a->ngr][3];
    double * group_vof_sum = new double[a->ngr];
    double * group_max_vof = new double[a->ngr];
    double * group_max_dist = new double[a->ngr];
    double * group_mean_dist = new double[a->ngr];

    for (int igr = 0; igr < a->ngr; ++igr) {
      group_liquid_volume[igr] = 0.0;
      group_liquid_surface[igr] = 0.0;
      FOR_I3 group_x[igr][i] = 0.0;
      FOR_I3 group_u[igr][i] = 0.0;
      group_vof_sum[igr] = 0.0;
      group_max_vof[igr] = 0.0;
      group_max_dist[igr] = 0.0;
      group_mean_dist[igr] = 0.0;
    }

    FOR_ICV {
      const int igr = a->group[icv];
      group_liquid_volume[igr] += vofs[icv]*vol_cv[icv];
      //group_liquid_surface[igr] += 0.65*pow(vof[icv]*vol_cv[icv], 2.0/3.0);
      //group_liquid_surface[igr] += plic_area[icv]; 
      FOR_I3 group_x[igr][i] += vofs[icv]*vol_cv[icv]*x_cv[icv][i];
      FOR_I3 group_u[igr][i] += vofs[icv]*vol_cv[icv]*u[icv][i];
      group_vof_sum[igr] += vofs[icv];
      if (iband[icv] > nband) group_max_vof[igr] += 1.0;
    }

    a->push(group_liquid_volume,  group_liquid_volume+a->ngr_a,  ADD_DATA);
    a->push(group_liquid_surface, group_liquid_surface+a->ngr_a, ADD_DATA);
    a->push(group_x,           group_x+a->ngr_a,           ADD_DATA);
    a->push(group_u,             group_u+a->ngr_a,             ADD_DATA);
    a->push(group_vof_sum,       group_vof_sum+a->ngr_a,       ADD_DATA);
    a->push(group_max_vof,       group_max_vof+a->ngr_a,       ADD_DATA);

    // go through the active groups (these are the groups we own with reduced
    // values)  and look at the total volume...
    int * group_flag = new int[a->ngr];
    // buffers for gathering drops to root for injection
    double * bufDropM = new double[a->ngr_a];
    double * bufDropS = new double[a->ngr_a];
    double (*bufDropX)[3] = new double[a->ngr_a][3];
    double (*bufDropU)[3] = new double[a->ngr_a][3];

    int ntrans = 0;
    int nblobs = 0;
    for (int igr = 0; igr < a->ngr_a; ++igr) {
      group_flag[igr] = 0;
      if (group_vof_sum[igr] > vof_zero) {
	FOR_I3 group_x[igr][i] /= group_liquid_volume[igr];
	FOR_I3 group_u[igr][i] /= group_liquid_volume[igr];
        nblobs++;
        // cout << "Group surface, volume , SMD = " << igr << " " <<  group_liquid_surface[igr] << " " << group_liquid_volume[igr] << " " << 6.0*group_liquid_volume[igr]/group_liquid_surface[igr] << endl;
        if (group_max_vof[igr] >= 1.0) { // transfer totally separated drops
          // transfer this drop... -> add it to the send buffer
          bufDropM[ntrans] = group_liquid_volume[igr]*rho_l;
          bufDropS[ntrans] = group_liquid_surface[igr];
          FOR_I3 {
            bufDropX[ntrans][i] = group_x[igr][i];
            bufDropU[ntrans][i] = group_u[igr][i];
          }
          ntrans++;
          group_flag[igr] = 1;
        }
      }
    }

    /*
    FOR_ICV {
      const int igr = a->group[icv];
      group_max_dist[igr] = max(group_max_dist[igr], DIST(x_cv[icv],group_x[igr])); 
      group_mean_dist[igr] += vol_cv[icv];
    }
    a->push(group_max_dist,       group_max_dist+a->ngr_a,       MAX_DATA);
    a->push(group_mean_dist,       group_mean_dist+a->ngr_a,       ADD_DATA);
    for (int igr = 0; igr < a->ngr_a; ++igr) {
      if (group_vof_sum[igr] > vof_zero) {
        const double radius = pow(3.0*group_mean_dist[igr]/4.0/M_PI,1.0/3.0);
        if ( radius < 10.0*group_max_dist[igr] ) { // transfer totally separated drops
          // transfer this drop... -> add it to the send buffer
          bufDropM[ntrans] = group_liquid_volume[igr]*rho_l;
          bufDropS[ntrans] = group_liquid_surface[igr];
          FOR_I3 {
            bufDropX[ntrans][i] = group_x[igr][i];
            bufDropU[ntrans][i] = group_u[igr][i];
          }
          ntrans++;
          group_flag[igr] = 1;
        }
      }
    }
    */
 
    // nblobs : total number of contiguous liquid structures ..
    double * blob_surf = new double[nblobs];
    double * blob_vol  = new double[nblobs];
    double (*blobX)[3] = new double[nblobs][3];
    int iblob = 0;
    for (int igr = 0; igr < a->ngr_a; ++igr) {
      if (group_vof_sum[igr] > vof_zero) {
        blob_surf[iblob] = group_liquid_surface[igr];
        blob_vol[iblob] = group_liquid_volume[igr];
        FOR_I3 blobX[iblob][i] = group_x[igr][i];
        iblob++;
        //cout << "contiguous liquid structure, iblob, surf, vol, x = " << iblob << " " << blob_surf[iblob] << " " << blob_vol[iblob] << " " << COUT_VEC(blobX[iblob]) << endl;
      }
    }
    assert(iblob==nblobs);
    int nblobs_global = 0;
    MPI_Allreduce(&nblobs, &nblobs_global, 1, MPI_INT, MPI_SUM, mpi_comm);
    if (mpi_rank == 0 && step%check_interval == 0 ) {
      cout << "Total number of contiguous liquidF  = " << " " << nblobs_global << endl;
    }


    //if (lsp) updateDropStats(blob_surf, blob_vol, blobX, nblobs, lpCls->lp,lpCls->size(),time,dt,step,1);
    //else  updateDropStats(blob_surf, blob_vol, blobX, nblobs,time,dt,step,1);


    delete[] blob_surf;
    delete[] blob_vol;
    delete[] blobX;

    delete[] group_liquid_volume;
    delete[] group_liquid_surface;
    delete[] group_x;
    delete[] group_u;
    delete[] group_vof_sum;
    delete[] group_max_vof;
    delete[] group_max_dist;
    delete[] vofs;
    delete[] group_mean_dist;


    // push the group flag out to the others...
    a->pull(group_flag+a->ngr_a,group_flag);
    
    // remove the drop from the VOF field
    FOR_ICV {
      const int igr = a->group[icv];
      if (group_flag[igr] == 1) {
        vof[icv] = 0.0;
      }
    }
    updateCvData(vof);

    // insert LSP drop
    dropTransfer(bufDropX,bufDropU,bufDropM, bufDropS, ntrans);
    
    //now update interface and LSP after transfer
    updateInterface();
    updateBaseState();
    updateLspPrimitiveData();
    //relocateLsp<LspState>();
    relocateLsp_rayTrace<LspState>();
    recycleLsp();

    delete a;
    delete[] group_flag;

    // now restore back cv_flag
    FOR_ICV_G cv_flag[icv] = cv_flag_copy[icv];

    delete[] cv_flag_copy;
    delete[] bufDropM;
    delete[] bufDropX;
    delete[] bufDropU;
    delete[] bufDropS;

  }

  void dropTransfer(double (*bufDropX)[3], double (*bufDropU)[3], double* bufDropM, double*bufDropS, int ntrans) {

    // reduce drops to root for injection into dpm. We can use the group buffers for this, since all information has been
    int npnew_global = 0;
    MPI_Allreduce(&ntrans, &npnew_global, 1, MPI_INT, MPI_SUM, mpi_comm);
    if (mpi_rank == 0 && step%check_interval == 0 ) {
      cout << "Total number of blobs injected to LSP:::" << " " << npnew_global << endl;
    }

    int *bufDropN = new int[ntrans];

    //compute number of drops ...
    for (int ip=0; ip<ntrans; ip++) {
      const double Dmin = 4.0*bufDropM[ip]/rho_l/bufDropS[ip];
      const double L = bufDropS[ip]/M_PI/Dmin;
      bufDropN[ip] = 1;
      //spericity criterion ....
      if (L/Dmin >= 3.0) {
        bufDropN[ip] = floor(L/(0.697*Dmin));
        //cout << "number of subdrops =" << " " << bufDropN[ip] << " " << L << " " << Dmin <<  endl;
      }
    }

    if ( npnew_global > 0 ) {

      int lspsize[mpi_size];
      MPI_Allgather(&ntrans, 1, MPI_INT, lspsize, 1, MPI_INT,  mpi_comm);

      int shift = 0;
      for (int irank =0; irank < mpi_rank; irank++) shift += lspsize[irank];

      double (*DropX)[3] = new double[npnew_global][3];
      double (*DropU)[3] = new double[npnew_global][3];
      double (*DropM) = new double[npnew_global];
      int (*DropN) = new int[npnew_global];
      double (*myDropX)[3] = new double[npnew_global][3];
      double (*myDropU)[3] = new double[npnew_global][3];
      double (*myDropM)    = new double[npnew_global];
      int (*myDropN)    = new int[npnew_global];

      for (int ip=0; ip<npnew_global; ip++) {
        FOR_I3 myDropX[ip][i] = 0.0;
        FOR_I3 myDropU[ip][i] = 0.0;
        myDropM[ip] = 0.0;
        myDropN[ip] = 0;
      }

      for (int ip = 0; ip < ntrans; ip++) {
        FOR_I3 myDropX[ip+shift][i] = bufDropX[ip][i];
        FOR_I3 myDropU[ip+shift][i] = bufDropU[ip][i];
        myDropM[ip+shift] = bufDropM[ip];
        myDropN[ip+shift] = bufDropN[ip];
      }

      MPI_Allreduce(myDropX, DropX, npnew_global*3, MPI_DOUBLE, MPI_SUM,  mpi_comm);
      MPI_Allreduce(myDropU, DropU, npnew_global*3, MPI_DOUBLE, MPI_SUM,  mpi_comm);
      MPI_Allreduce(myDropM, DropM, npnew_global, MPI_DOUBLE, MPI_SUM,  mpi_comm);
      MPI_Allreduce(myDropN, DropN, npnew_global, MPI_INT, MPI_SUM,  mpi_comm);

      // figure out which particles are ours...
      // We are going to use global cv indexing to do this...

      int8 * cvora = NULL;
      buildXora(cvora,ncv);
      assert( cvora[mpi_rank+1] - cvora[mpi_rank] == ncv );

      vector<int> cvList;
      int * my_icvp_global = new int[npnew_global];

      for (int ip = 0; ip < npnew_global; ++ip) {
        my_icvp_global[ip] = -1;
        cvAdt->buildListForPoint(cvList, DropX[ip]);

        const int icv = getClosestICVFromList(cvList,DropX[ip]);
        if (icv >= 0) my_icvp_global[ip] = icv + cvora[mpi_rank]; // local to global cv offset
      }
      // at this point, my_icvp_global should contain the global cv index of
      // the containing cell, or -1... get the maximum across the
      // processors to break ties or check for particles that could not be located...

      int * icvp_global = new int[npnew_global];
      MPI_Allreduce(my_icvp_global,icvp_global,npnew_global,MPI_INT,MPI_MAX,mpi_comm);

      int npnew = 0;
      for (int ip = 0; ip < npnew_global; ++ip) {
        if (icvp_global[ip] == -1) {
          if (mpi_rank == 0)
            cout << "Warning: could not locate cv for particle at: " <<
              COUT_VEC(DropX[ip]) << endl;
        }
        else if (icvp_global[ip] == my_icvp_global[ip]) {
          for (int ichild=0; ichild<DropN[ip]; ichild++) {
            ++npnew;
          }
        }
      }

      // resize and add particles on processors that need it...
      if (npnew > 0) {

        int npold = lpCls->size();
        lpCls->resize(npold+npnew);

        // loop through all the particles again
        int npnew_check = 0;
        for (int ip = 0; ip < npnew_global; ++ip) {
          // if the particles belongs to the processors
          if ((icvp_global[ip] != -1)&&(icvp_global[ip] == my_icvp_global[ip])) {
            for (int ichild=0; ichild<DropN[ip]; ichild++) {
              int ipnew = npnew_check++;
              // local cv index
              lpCls->lp[npold+ipnew].icv = my_icvp_global[ip] - cvora[mpi_rank];
              // data...
              FOR_I3 lpCls->lp[npold+ipnew].xp[i]  = DropX[ip][i];
              FOR_I3 lpCls->lp[npold+ipnew].up[i]  = DropU[ip][i];
              lpCls->lp[npold+ipnew].mp            = DropM[ip];
              lpCls->lp[npold+ipnew].Tp            = T_ref; // hack... fixed to T_ref temperature for now
              // set the source terms to zero...
              lpCls->lp[npold+ipnew].tbu           = 0.0;  // breakup time
              lpCls->lp[npold+ipnew].flag          = KEEP; // recycling var
              lpCls->lp[npold+ipnew].npar          = 1.0;  // parcel info
            }
          }
        }
        assert(npnew_check == npnew);
      }
      delete[] cvora;
      delete[] my_icvp_global;
      delete[] icvp_global;
      delete[] myDropX;
      delete[] myDropU;
      delete[] myDropM;
      delete[] myDropN;
      delete[] DropX;
      delete[] DropU;
      delete[] DropM;
      delete[] DropN;
      delete[] bufDropN;
    }


  }

  void initCvAverageOperator(AverageOperator * &averageOp, int * fa_flag) {

    // this one uses fa_flag == 1 to join nbring cvs into a group...

    if (mpi_rank == 0 && step%check_interval == 0 )
      cout << "initCvAverageOperator()" << endl;

    int * next = new int[ncv];
    int * prev = new int[ncv];

    FOR_ICV next[icv] = prev[icv] = -1;

    // internal faces...
    for (int ifa = 0; ifa < nfa_i; ++ifa) {
      if (fa_flag[ifa] == 1) {
        int icv0 = cvofa[ifa][0];
        assert((icv0 >= 0)&&(icv0 < ncv));
        int icv1 = cvofa[ifa][1];
        assert((icv1 >= 0)&&(icv1 < ncv));
        // these cvs should be linked...
        linkCvs(icv0,icv1,prev,next);
      }
    }

    FOR_ICV cv_flag[icv] = -1;

    int ngr = 0;
    FOR_ICV {
      if (cv_flag[icv] == -1) {
        cv_flag[icv] = ngr;
        int icv_ = icv;
        while (prev[icv_] != -1) {
          icv_ = prev[icv_];
          assert(cv_flag[icv_] == -1);
          cv_flag[icv_] = ngr;
        }
        icv_ = icv;
        while (next[icv_] != -1) {
          icv_ = next[icv_];
          assert(cv_flag[icv_] == -1);
          cv_flag[icv_] = ngr;
        }
        ++ngr;
      }
    }

    delete[] next;
    delete[] prev;

    // cv_flag now contains the unique group number for our local cvs.
    initCvAverageOperatorCommon(averageOp, ngr, fa_flag);

  }

  void initCvAverageOperatorCommon(AverageOperator * &averageOp, int ngr, int * fa_flag) {

    // we now have a unique index in cv_flag that describes each cv's group. We also
    // have a 1 in any fa_flag that connects groups across processors. Now
    // reduce across processors to combine groups.

    int * cvogr_i = new int[ngr+1];
    FOR_IGR cvogr_i[igr+1] = 0;
    FOR_ICV {
      const int igr = cv_flag[icv];
      cvogr_i[igr+1] += 1;
    }
    cvogr_i[0] = 0;
    FOR_IGR cvogr_i[igr+1] += cvogr_i[igr];
    assert(cvogr_i[ngr] == ncv);
    int * cvogr_v = new int[ncv];
    FOR_ICV {
      const int igr = cv_flag[icv];
      cvogr_v[cvogr_i[igr]] = icv;
      cvogr_i[igr] += 1;
    }
    for (int igr = ngr-1; igr > 0; --igr)
      cvogr_i[igr] = cvogr_i[igr-1];
    cvogr_i[0] = 0;

    // check...

    FOR_IGR {
      for (int cog = cvogr_i[igr]; cog != cvogr_i[igr+1]; ++cog) {
        const int icv = cvogr_v[cog];
        assert(cv_flag[icv] == igr);
      }
    }

    // now build global groups...

    int offset;
    MPI_Scan(&ngr,&offset,1,MPI_INT,MPI_SUM,mpi_comm);
    offset -= ngr;

    FOR_ICV cv_flag[icv] += offset;
    const int nfa_bpi = nfa-nfa_i;
    int * fa_flag_copy = new int[nfa];

    int done = 0;
    while (done == 0) {

      int my_done = 1;

      // loop processor boundaries including periodic faces
      for (int ifa = nfa_i; ifa < nfa; ++ifa) {
        const int icv0 = cvofa[ifa][0];
        fa_flag_copy[ifa] = cv_flag[icv0];
      }

      //replace ghost data
      updateFaData(fa_flag_copy,REPLACE_DATA);
      for (int ifa = nfa_i; ifa < nfa; ++ifa) {
        if (fa_flag[ifa] == 1) {
          const int icv0 = cvofa[ifa][0];
          if (fa_flag_copy[ifa] < cv_flag[icv0]) {
            my_done = 0;
            cv_flag[icv0] = fa_flag_copy[ifa];
          }
        }
      }

      FOR_IGR {
        int igr_global_min = igr + offset;
        for (int cvg = cvogr_i[igr]; cvg != cvogr_i[igr+1]; ++cvg) {
          const int icv = cvogr_v[cvg];
          igr_global_min = min(igr_global_min,cv_flag[icv]);
        }
        for (int cvg = cvogr_i[igr]; cvg != cvogr_i[igr+1]; ++cvg) {
          const int icv = cvogr_v[cvg];
          if (cv_flag[icv] !=  igr_global_min) {
            assert(cv_flag[icv] > igr_global_min);
            cv_flag[icv] = igr_global_min;
          }
        }
      }

      MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);

    }

    vector< std::pair<int,int> > icvVec(ncv);
    FOR_ICV {
      icvVec[icv].first = cv_flag[icv];
      icvVec[icv].second = icv;
    }
    std::sort(icvVec.begin(),icvVec.end(),compareFirstSecondIntIntPair);
    int ngr_ag = 0;
    int ngr_a = 0;
    int igr_global = -1;
    FOR_ICV {
      if (icvVec[icv].first != igr_global) {
        igr_global = icvVec[icv].first;
        ++ngr_ag;
        if (igr_global >= offset)
          ++ngr_a;
      }
      if (igr_global >= offset)
        cv_flag[icvVec[icv].second] = ngr_a-1;
      else
        cv_flag[icvVec[icv].second] = -1;
    }

    int * grora = NULL;
    buildXora(grora,ngr_a);

    if (mpi_rank == 0)
      cout << " > global number of groups: " << grora[mpi_size] << endl;

    FOR_ICV {
      if (cv_flag[icv] >= 0)
        cv_flag[icv] += grora[mpi_rank];
      else
        cv_flag[icv] = grora[mpi_size]; // one larger than the maximum possible
    }

    done = 0;
    while (done == 0) {

      int my_done = 1;

      for (int ifa = nfa_i; ifa < nfa; ++ifa) {
        const int icv0 = cvofa[ifa][0];
        fa_flag_copy[ifa] = cv_flag[icv0];
      }
      updateFaData(fa_flag_copy,REPLACE_DATA);
      for (int ifa = nfa_i; ifa < nfa; ++ifa) {
        if (fa_flag[ifa] == 1) {
          const int icv0 = cvofa[ifa][0];
          if (fa_flag_copy[ifa] < cv_flag[icv0]) {
            assert(cv_flag[icv0] == grora[mpi_size]);
            cv_flag[icv0] = fa_flag_copy[ifa];
          }
        }
      }

      FOR_IGR {
        int igr_global_min = grora[mpi_size];
        for (int cvg = cvogr_i[igr]; cvg != cvogr_i[igr+1]; ++cvg) {
          const int icv = cvogr_v[cvg];
          igr_global_min = min(igr_global_min,cv_flag[icv]);
        }
        if (igr_global_min == grora[mpi_size])
          my_done = 0;
        else {
          for (int cvg = cvogr_i[igr]; cvg != cvogr_i[igr+1]; ++cvg) {
            const int icv = cvogr_v[cvg];
            if (cv_flag[icv] != igr_global_min) {
              assert(cv_flag[icv] == grora[mpi_size]);
              cv_flag[icv] = igr_global_min;
            }
          }
        }
      }

      MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);

    }

    delete[] cvogr_i;
    delete[] cvogr_v;
    delete[] fa_flag_copy;

    // and use the icvVec (already sorted - should still be correct) to
    // build the group numbers...

    int * my_igr_global = new int[ngr_ag-ngr_a];
    int ngr_a_check = 0;
    int ngr_ag_check = ngr_a;
    igr_global = -1;
    int igr_global_new = -1;
    FOR_ICV {
      if (icvVec[icv].first != igr_global) {
        assert(icvVec[icv].first > igr_global);
        igr_global = icvVec[icv].first;
        assert(cv_flag[icvVec[icv].second] > igr_global_new);
        igr_global_new = cv_flag[icvVec[icv].second];
        if (igr_global_new >= grora[mpi_rank]) {
          assert(igr_global_new < grora[mpi_rank+1]);
          assert(ngr_ag_check == ngr_ag); // other groups should be done...
          assert(igr_global_new == ngr_a_check + grora[mpi_rank]);
          ++ngr_a_check;
        }
        else {
          my_igr_global[ngr_ag_check-ngr_a] = igr_global_new;
          ++ngr_ag_check;
        }
      }
      assert(cv_flag[icvVec[icv].second] == igr_global_new);
      if (igr_global_new >= grora[mpi_rank])
        cv_flag[icvVec[icv].second] = ngr_a_check-1; // put the groups we own first
      else
        cv_flag[icvVec[icv].second] = ngr_ag_check-1; // then the other groups...
    }
    assert(ngr_a_check == ngr_a);
    assert(ngr_ag_check == ngr_ag);
    icvVec.clear();

    // now build the average stuff...

    assert(averageOp == NULL);
    averageOp = new AverageOperator(my_igr_global,ngr_ag-ngr_a,grora);

    delete[] my_igr_global;
    delete[] grora;

    averageOp->n = ncv;
    averageOp->ngr = ngr_ag; // not passed ngr, this one could be smaller...
    averageOp->ngr_a = ngr_a;

    // the group number is stored in cv_flag...
    averageOp->group = new int[ncv];
    FOR_ICV averageOp->group[icv] = cv_flag[icv];


   // if (cv_volume == NULL)
   //   calcCvGeometry();
    averageOp->weight = vol_cv;

    // and reduce the group weights...

    double * weight_sum = new double[ngr];
    FOR_IGR weight_sum[igr] = 0.0;
    FOR_ICV {
      const int igr = averageOp->group[icv];
      weight_sum[igr] += averageOp->weight[icv];
    }

    averageOp->push(weight_sum,weight_sum+ngr_a,ADD_DATA);
    dumpRange(weight_sum,ngr_a,"WEIGHT_SUM");

    averageOp->inv_weight_sum = new double[ngr_a];
    FOR_ACTIVE_IGR averageOp->inv_weight_sum[igr] = 1.0/weight_sum[igr];

    delete[] weight_sum;


  }



  void calcSgs() {

    if ( sgs_model == "NONE") {
      computeSgsNone();
    } else if ( sgs_model == "VREMAN") {
      computeSgsVreman();
    } else {
      // this error is checked earlier, so we shouldnt get here..
      assert(0);
    }

  }

  //=============================================================
  // sgs routines
  //=============================================================

  void parseSgsModel(string& sgs_model) {
    const string sgs_type = getStringParam("SGS_MODEL", "NONE");
    if (sgs_type == "NONE") {
      sgs_model = "NONE";
    } else if ( sgs_type == "VREMAN") {
      sgs_model = "VREMAN";
    } else {
      CERR(" > unrecognized sgs model: " << sgs_type);
    }
  }

  void computeSgsNone() {
    FOR_ICV_G {
      mu_sgs[icv]  = 0.0;
    }
  }

  void computeSgsVreman() {

    const double vreman_coeff = getDoubleParam("VREMAN_COEFF", 0.07);

    // need to compute a divided difference of |S| so we need to
    // populate the ghosts here...
    double * smag     = new double[ncv_g];
    double * lap_smag = new double[ncv];
    for (int icv = 0; icv < ncv; ++icv) {
      double smag2 = 0.0;
      FOR_I3 FOR_J3 smag2 += 2.0*dudx[icv][i][j]*dudx[icv][i][j];
      smag[icv]     = sqrt(smag2);
      lap_smag[icv] = 0.0;
    }
    updateCvDataStart(smag);

    FOR_INTERNAL_IFA {
      const int icv0  = cvofa[ifa][0];
      const int icv1  = cvofa[ifa][1];
      const double dd = smag[icv1] - smag[icv0];

      // square length scale is built in here...
      lap_smag[icv0] += dd;
      lap_smag[icv1] -= dd;
    }

    updateCvDataFinish(smag);

    FOR_INTERPROC_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g));
      const double dd = smag[icv1] - smag[icv0];
      lap_smag[icv0] += dd;
    }

    for (int icv = 0; icv < ncv; ++icv) {

      const double dx2 = pow( vol_cv[icv], 2.0/3.0 ); // XXX precompute..

      const double alpha11   = dudx[icv][0][0];
      const double alpha11_2 = alpha11*alpha11;
      const double alpha22   = dudx[icv][1][1];
      const double alpha22_2 = alpha22*alpha22;
      const double alpha33   = dudx[icv][2][2];
      const double alpha33_2 = alpha33*alpha33;

      const double alpha12   = dudx[icv][0][1];
      const double alpha12_2 = alpha12*alpha12;
      const double alpha13   = dudx[icv][0][2];
      const double alpha13_2 = alpha13*alpha13;
      const double alpha23   = dudx[icv][1][2];
      const double alpha23_2 = alpha23*alpha23;

      const double alpha21   = dudx[icv][1][0];
      const double alpha21_2 = alpha21*alpha21;
      const double alpha31   = dudx[icv][2][0];
      const double alpha31_2 = alpha31*alpha31;
      const double alpha32   = dudx[icv][2][1];
      const double alpha32_2 = alpha32*alpha32;

      const double beta11  = dx2*(alpha11_2+alpha21_2+alpha31_2);
      const double beta22  = dx2*(alpha12_2+alpha22_2+alpha32_2);
      const double beta33  = dx2*(alpha13_2+alpha23_2+alpha33_2);

      const double beta12  = dx2*(alpha11*alpha12+alpha21*alpha22+alpha31*alpha32);
      const double beta13  = dx2*(alpha11*alpha13+alpha21*alpha23+alpha31*alpha33);
      const double beta23  = dx2*(alpha12*alpha13+alpha22*alpha23+alpha32*alpha33);

      double B       = (beta11*beta22-beta12*beta12)+(beta11*beta33-beta13*beta13)+(beta22*beta33-beta23*beta23);
      B              = (B + abs(B))*0.5;

      const double alal    =
        ((alpha11_2+alpha22_2)+(alpha33_2+alpha12_2))+
        ((alpha13_2+alpha23_2)+(alpha21_2+alpha31_2))+alpha32_2;

      const double s_mag = sqrt(B/(alal+1.0E-20)); // includes lengthscale squared too

      // mu_sgs...
      mu_sgs[icv] = rho[icv]*vreman_coeff*s_mag;

      // clip at the constant smagorinsky level with c2 = 0.5^2...
      mu_sgs[icv] = min(mu_sgs[icv],rho[icv]*0.05*dx2*sqrt(2.0*alal));


    }//for_icv

    updateCvData(mu_sgs);

    delete[] smag;
    delete[] lap_smag;
  }

  void syncPostState() {

    updateCvData(vof); // assume vof has been read or set
    FOR_ICV_G {
      rho[icv] = rho_g*(1.0-vof[icv]) + rho_l*vof[icv];
      mu_lam[icv] = mu_g*(1.0-vof[icv]) + mu_l*vof[icv];
    }
    updateCvData(u);
    FOR_ICV_G p[icv] = 0.0; // for initial guess...

    // calculate the interface properties...
    updateInterface();

    StaticSolver::calcCvGrad(dudx,u);
    updateCvData(dudx);

    // subgrid stress...

    calcSgs();
  };

  virtual void initialHook() {}
  virtual void temporalHook() {}
  virtual void finalHook() {}

  void computeForces() {
    assert(forceZoneBuf&&!forceZoneMap.empty());
    for (map<const string,int>::const_iterator iter = forceZoneMap.begin(); iter!=forceZoneMap.end(); ++iter){
      double f_buf[3][3], m_buf[3][3];
      bc_map[iter->first]->force(f_buf,m_buf);
      FOR_I3 FOR_J3 forceZoneBuf[iter->second][3*i+j] = f_buf[i][j];
      FOR_I3 FOR_J3 momentZoneBuf[iter->second][3*i+j] = m_buf[i][j];
    }
  }

  // vof solver special function evaluations

//  virtual CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,const string& name,
//                                                    list<CtiRegister::CtiData>& args, const bool b_eval_func) {

//    using namespace CtiRegister;

    // this solver doesnt know how to evaluate this data -- go to StaticSolver

//    return FlowSolver::funcEvalCtiData(v,name,args,b_eval_func);
//  }

  //==========================================================================
  // Vof routines...
  //==========================================================================

  void calcGfromVof() {

    FOR_ICV {
      if (vof[icv] >= vof_zero && vof[icv] <= 1.0-vof_zero) {
        g[icv] = atanh(2.0*vof[icv]-1.0)/beta[icv];
      }
      else {
        if (vof[icv] < 0.5)
          g[icv] = -1.0E+20;
        else
          g[icv] = +1.0E+20;
      }
      if (g[icv] != g[icv]) {
        if (vof[icv] < 0.5)
          g[icv] = -1.0E+20;
        else
          g[icv] = +1.0E+20;
      }
    }
    updateCvData(g);

  }

  void calcGasCfl(double *cfl, const double dt_) const {

    // functions returns the

    FOR_ICV_G cfl[icv] = 0.0;

    for (int ibf = 0; ibf < nbf; ++ibf) {
      const int icv = cvobf[ibf];
      double dx[3];
      FOR_I3 dx[i] = x_bf[ibf][i] - x_cv[icv][i];
      double normal = -DOT_PRODUCT(dx,n[icv]);
      const double vof =  0.5*(1.0+tanh(beta[icv]*(normal+g[icv])));
      cfl[icv] += (1.0-vof)*max(0.0, DOT_PRODUCT(u[icv],n_bf[ibf]));
    }

    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      const double q_mid = 1.25*q_fa[ifa]-0.25*q0_fa[ifa];
      double dx0[3], dx1[3],dx[3];
      FOR_I3 dx0[i] = x_fa[ifa][i] - x_cv[icv0][i];
      FOR_I3 dx1[i] = x_fa[ifa][i] - x_cv[icv1][i];

      double normal0 = -DOT_PRODUCT(dx0,n[icv0]);
      double normal1 = -DOT_PRODUCT(dx1,n[icv1]);
      const double vof0 =  0.5*(1.0+tanh(beta[icv0]*(normal0+g[icv0])));
      const double vof1 =  0.5*(1.0+tanh(beta[icv1]*(normal1+g[icv1])));

      double vof_f;

      if (q_mid >= 0.0) {
        double beta;
        if (vof[icv1] != vof[icv0] ) beta = (vof0 - vof[icv0])/(vof[icv1]-vof[icv0]);
        else beta = 1.0;
        beta = max(0.0, min(beta,1.0));
        beta = pow(beta,3.0);
        vof_f = (1.0-beta)*vof0 + beta*vof[icv1];
        cfl[icv0] += (1.0-vof_f)*q_mid;
      }
      else {
        double beta;
        if (vof[icv1] != vof[icv0] ) beta = (vof1 - vof[icv1])/(vof[icv0]-vof[icv1]);
        else beta = 1.0;
        beta = max(0.0, min(beta,1.0));
        beta = pow(beta,3.0);
        vof_f = (1.0-beta)*vof1+beta*vof[icv0];
        cfl[icv1] += (1.0-vof_f)*q_mid;
      }
    }

    FOR_ICV cfl[icv] *= dt_/vol_cv[icv];

    double my_cfl_max = 0.0;
    double cfl_max = 0.0;
    for (int icv = 0; icv < ncv; ++icv)
      my_cfl_max = max(my_cfl_max,cfl[icv]);

    MPI_Allreduce(&my_cfl_max, &cfl_max, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    COUT1(" > Gas CFL max = " << cfl_max);

  }

  void calcGasCfl2(double *cfl, const double dt_) const {

    // functions returns the

    FOR_ICV_G cfl[icv] = 0.0;

    for (int ibf = 0; ibf < nbf; ++ibf) {
      const int icv = cvobf[ibf];
      double dx[3];
      FOR_I3 dx[i] = x_bf[ibf][i] - x_cv[icv][i];
      double normal = -DOT_PRODUCT(dx,n[icv]);
      const double vof =  0.5*(1.0+tanh(beta[icv]*(normal+g[icv])));
      cfl[icv] += (1.0-vof)*max(0.0, DOT_PRODUCT(u[icv],n_bf[ibf]));
    }

    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      const double q_mid = 1.25*q_fa[ifa]-0.25*q0_fa[ifa];
      double dx0[3], dx1[3],dx[3];
      FOR_I3 dx0[i] = x_fa[ifa][i] - x_cv[icv0][i];
      FOR_I3 dx1[i] = x_fa[ifa][i] - x_cv[icv1][i];

      double normal0 = -DOT_PRODUCT(dx0,n[icv0]);
      double normal1 = -DOT_PRODUCT(dx1,n[icv1]);
      const double vof0 =  0.5*(1.0+tanh(beta[icv0]*(normal0+g[icv0])));
      const double vof1 =  0.5*(1.0+tanh(beta[icv1]*(normal1+g[icv1])));

      double dvofdx0 = 0.5*(1.0-pow(tanh(beta[icv0]*(normal0+g[icv0])),2))*beta[icv0];
      double dvofdx1 = 0.5*(1.0-pow(tanh(beta[icv1]*(normal1+g[icv1])),2))*beta[icv1];
      double dvofdx = 0.5*fabs(dvofdx0+dvofdx1);
      double gamma = 0.0;
      if (dvofdx > 1.0E-10) {
        gamma = fabs(vof[icv1]-vof[icv0])/MAG(dx)/dvofdx;
      }
      else {
        gamma = 0.0;
      }

      gamma = max(0.0, min(gamma,1.0));
      gamma = pow(1.0 - gamma,2);

      gamma = 1.0;

      double vof_f;

      if (q_mid >= 0.0) {
        vof_f = (1.0-gamma)*vof0 + 0.5*gamma*(vof[icv0]+vof[icv1]);
        cfl[icv0] += fabs((1.0-vof_f)*q_mid);
      }
      else {
        vof_f = (1.0-gamma)*vof1+0.5*gamma*(vof[icv0]+vof[icv1]);
        cfl[icv1] += fabs((1.0-vof_f)*q_mid);
      }
    }

    FOR_ICV cfl[icv] *= dt_/vol_cv[icv];

    double my_cfl_max = 0.0;
    double cfl_max = 0.0;
    for (int icv = 0; icv < ncv; ++icv)
      my_cfl_max = max(my_cfl_max,cfl[icv]);

    MPI_Allreduce(&my_cfl_max, &cfl_max, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    COUT1(" > Gas outgoing CFL max = " << cfl_max);

  }



  void calcLiquidCfl2(double *cfl, const double dt_) const {

    // functions returns the

    FOR_ICV_G cfl[icv] = 0.0;

    for (int ibf = 0; ibf < nbf; ++ibf) {
      const int icv = cvobf[ibf];
      double dx[3];
      FOR_I3 dx[i] = x_bf[ibf][i] - x_cv[icv][i];
      double normal = -DOT_PRODUCT(dx,n[icv]);
      const double vof =  0.5*(1.0+tanh(beta[icv]*(normal+g[icv])));
      cfl[icv] += vof*max(0.0, DOT_PRODUCT(u[icv],n_bf[ibf]));
    }

    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      const double q_mid = 1.25*q_fa[ifa]-0.25*q0_fa[ifa];
      double dx0[3], dx1[3],dx[3];
      FOR_I3 dx0[i] = x_fa[ifa][i] - x_cv[icv0][i];
      FOR_I3 dx1[i] = x_fa[ifa][i] - x_cv[icv1][i];

      double normal0 = -DOT_PRODUCT(dx0,n[icv0]);
      double normal1 = -DOT_PRODUCT(dx1,n[icv1]);
      const double vof0 =  0.5*(1.0+tanh(beta[icv0]*(normal0+g[icv0])));
      const double vof1 =  0.5*(1.0+tanh(beta[icv1]*(normal1+g[icv1])));


      double dvofdx0 = 0.5*(1.0-pow(tanh(beta[icv0]*(normal0+g[icv0])),2))*beta[icv0];
      double dvofdx1 = 0.5*(1.0-pow(tanh(beta[icv1]*(normal1+g[icv1])),2))*beta[icv1];
      double dvofdx = 0.5*fabs(dvofdx0+dvofdx1);
      double gamma = 0.0;
      if (dvofdx > 1.0E-10) {
        gamma = fabs(vof[icv1]-vof[icv0])/MAG(dx)/dvofdx;
      }
      else {
        gamma = 0.0;
      }

      gamma = max(0.0, min(gamma,1.0));
      gamma = pow(1.0 - gamma,2);

      gamma = 1.0;


      double vof_f;

      if (q_mid >= 0.0) {
        vof_f = (1.0-gamma)*vof0 + 0.5*gamma*(vof[icv0]+vof[icv1]);
        cfl[icv0] += fabs(vof_f*q_mid);
      }
      else {
        vof_f = (1.0-gamma)*vof1+0.5*gamma*(vof[icv0]+vof[icv1]);
        cfl[icv1] += fabs(vof_f*q_mid);
      }
    }

    FOR_ICV cfl[icv] *= dt_/vol_cv[icv];

    double my_cfl_max = 0.0;
    double cfl_max = 0.0;
    for (int icv = 0; icv < ncv; ++icv)
      my_cfl_max = max(my_cfl_max,cfl[icv]);

    MPI_Allreduce(&my_cfl_max, &cfl_max, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    COUT1(" > Liquid out going CFL max = " << cfl_max);

  }




  void calcLiquidCfl(double *cfl, const double dt_) const {

    // functions returns the

    FOR_ICV_G cfl[icv] = 0.0;

    for (int ibf = 0; ibf < nbf; ++ibf) {
      const int icv = cvobf[ibf];
      double dx[3];
      FOR_I3 dx[i] = x_bf[ibf][i] - x_cv[icv][i];
      double normal = -DOT_PRODUCT(dx,n[icv]);
      const double vof =  0.5*(1.0+tanh(beta[icv]*(normal+g[icv])));
      cfl[icv] += vof*max(0.0, DOT_PRODUCT(u[icv],n_bf[ibf]));
    }

    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      const double q_mid = 1.25*q_fa[ifa]-0.25*q0_fa[ifa];
      double dx0[3], dx1[3],dx[3];
      FOR_I3 dx0[i] = x_fa[ifa][i] - x_cv[icv0][i];
      FOR_I3 dx1[i] = x_fa[ifa][i] - x_cv[icv1][i];

      double normal0 = -DOT_PRODUCT(dx0,n[icv0]);
      double normal1 = -DOT_PRODUCT(dx1,n[icv1]);
      const double vof0 =  0.5*(1.0+tanh(beta[icv0]*(normal0+g[icv0])));
      const double vof1 =  0.5*(1.0+tanh(beta[icv1]*(normal1+g[icv1])));

      double vof_f;

      if (q_mid >= 0.0) {
        double beta;
        if (vof[icv1] != vof[icv0] ) beta = 0.5*fabs((vof0 - vof[icv0])/(vof[icv1]-vof[icv0]));
        else beta = 1.0;
        beta = max(0.0, min(beta,1.0));
        beta = pow(beta,3.0);
        vof_f = (1.0-beta)*vof0 + beta*vof[icv1];
        cfl[icv0] += vof_f*q_mid;
      }
      else {
        double beta;
        if (vof[icv1] != vof[icv0] ) beta = 0.5*fabs((vof1 - vof[icv1])/(vof[icv0]-vof[icv1]));
        else beta = 1.0;
        beta = max(0.0, min(beta,1.0));
        beta = pow(beta,3.0);
        vof_f = (1.0-beta)*vof1+beta*vof[icv0];
        cfl[icv1] += vof_f*q_mid;
      }
    }

    FOR_ICV cfl[icv] *= dt_/vol_cv[icv];

    double my_cfl_max = 0.0;
    double cfl_max = 0.0;
    for (int icv = 0; icv < ncv; ++icv)
      my_cfl_max = max(my_cfl_max,cfl[icv]);

    MPI_Allreduce(&my_cfl_max, &cfl_max, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    COUT1(" > Liquid CFL max = " << cfl_max);

  }


  void calcCfl2(double *cfl, const double dt_) const {

    // functions returns the

    FOR_ICV_G cfl[icv] = 0.0;

    for (vector<VofBc*>::const_iterator it = bcs.begin(); it != bcs.end(); ++it) {
      if ((*it)->q_bf != NULL) {
        for (int ibf = 0; ibf < (*it)->zone_ptr->nbf; ++ibf) {
          const int icv0 = (*it)->zone_ptr->cvobf[ibf];
          cfl[icv0] += max((*it)->q_bf[ibf],0.0); // 0.5*dt/vol below...
        }
      }
    }


    for (int ifa = 0; ifa<nfa; ++ifa){
      const int icv0 = cvofa[ifa][0];
      cfl[icv0] += max(q_fa[ifa],0.0);
      const int icv1 = cvofa[ifa][1];
      if (icv1 < ncv)
        cfl[icv1] -= min(q_fa[ifa],0.0);
    }


    FOR_ICV cfl[icv] *= dt_/vol_cv[icv];

    double my_cfl_max = 0.0;
    double cfl_max = 0.0;
    for (int icv = 0; icv < ncv; ++icv)
      my_cfl_max = max(my_cfl_max,cfl[icv]);

    MPI_Allreduce(&my_cfl_max, &cfl_max, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    COUT1(" > outgoing CFL max = " << cfl_max);

  }

  void flagInterfaceCvs() {

    // vof must be updated to at least 1st level ghosts at this point
    // basically flag cvs where 0 < vof < 1. also flag cv w/ vof = 1 next to cv w/ vof = 0...

    // cv_flag = -1 : a nbr has an interface (for cfl calc)
    // cv_flag = 0 : no interface
    // cv_flag = 1 : 0 < vof < 1
    // cv_flag = 2 : vof = 1 and atleast 1 nbr has vof = 0
    // cv_flag = 3 : isolated cell (determined later w/ PLIC n)
    // use vof_zero to control geometric support...

    FOR_ICV_G cv_flag[icv] = 0;

    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));

      if (vof[icv0] > vof_zero && vof[icv0] < (1.0-vof_zero)) {
        cv_flag[icv0] = 1;
      }
      else if (vof[icv0] >= (1.0-vof_zero) && vof[icv1] <= vof_zero) {
        cv_flag[icv0] = 2;
      }

      if (vof[icv1] > vof_zero && vof[icv1] < (1.0-vof_zero)) {
        cv_flag[icv1] = 1;
      }
      else if (vof[icv0] <= vof_zero && vof[icv1] >= (1.0-vof_zero)) {
        cv_flag[icv1] = 2;
      }

    }

    // flag interace nbrs as -1...

    FOR_IFA {
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      if (cv_flag[icv0] >= 1 && cv_flag[icv1] == 0)
        cv_flag[icv1] = -1;
      else if (cv_flag[icv0] == 0 && cv_flag[icv1] >= 1)
        cv_flag[icv0] = -1;
    }

    updateCvData(cv_flag);

  }

  void buildIband() {

    if (mpi_rank == 0 && step%check_interval == 0 )
      cout << " > Detecting subgrid vof cell..." << endl;

    double eps = vof_zero;

    const int isolated = 81;
    const int nband = 20;

    FOR_ICV_G iband[icv] = 0;

    FOR_ICV {
      if (vof[icv] >= 0.15)
        iband[icv] = 1; // well-resolved liquid interfaced
      else if (vof[icv] < vof_zero)
        iband[icv] = 0; // pus gas ....
      else
        iband[icv] = isolated; // pure interface 0<vof<1
    }

    updateCvData(iband);

    // construct iband structure
    // 0 : pure gas
    // 1 : well resolved  liquid
    // 2 ~ nband: band structure from well-resolved liquid cell
    // 81 : out of band .


    for (int index = 1; index < nband; ++index) {
      FOR_IFA {
        const int icv0 = cvofa[ifa][0];
        const int icv1 = cvofa[ifa][1];
        if (iband[icv0] == index && iband[icv1] > index+1 ) {
          iband[icv1] = index+1;
        }
        else if (iband[icv0] > index+1  && iband[icv1] == index) {
          iband[icv0] = index + 1;
        }
      }
      updateCvData(iband);
    }

//    FOR_ICV {
//      if (kappa[icv] > 1.0/r_vv[icv] ) iband[icv] = 20;
//      else if (kappa[icv] < -1.0/r_vv[icv]) iband[icv] = 20;
//    }

    FOR_ICV cv_flag_real[icv] = double(iband[icv]);
    dumpRange(cv_flag_real, ncv, "iband");
  }

/*

  void buildIband2() {

    if (mpi_rank == 0 && step%check_interval == 0 )
      cout << " > Detecting subgrid vof cell..." << endl;

    double eps = vof_zero;

    const int isolated = 81;
    const int subgrid = 10;

    FOR_ICV_G iband[icv] = 0;

    FOR_ICV {
      if (vof[icv] >= 0.5)
        iband[icv] = 1; // well-resolved liquid interfaced
      else if (vof[icv] < vof_zero)
        iband[icv] = 0; // pus gas ....
      else
        iband[icv] = isolated; // pure interface 0<vof<1
    }

    updateCvData(iband);

    // detect non-isolated interface cell
    FOR_ICV {
      if (iband[icv] == isolated ) {
        for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          assert(icv!=icv_nbr);
          if ( icv != icv_nbr) {
            if (iband[icv_nbr] != 0) {
              iband[icv] = subgrid; //non_isolated cell
              break;
            }
          }
        }
      }
    }

    updateCvData(iband);

    // confirm non-isolated interface cell
    FOR_ICV {
      if (iband[icv] == isolated ) {
        for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          assert(iband[icv_nbr] == 0);
        }
      }
    }

    // compute number of isolated cells.
    int my_np = 0;
    int np = 0;
    double my_vofsum = 0.0;
    double vofsum = 0.0;
    FOR_ICV {
      if (iband[icv] == isolated) {
        ++my_np;
        my_vofsum += vof[icv];
        //cout << "isolated drop, icv, vof, x, n = " << icv << " " << vof[icv] << " " << COUT_VEC(x_cv[icv]) << " "<< COUT_VEC(n[icv]) << endl;
        FOR_I3 n[icv][i] = -u[icv][i];
        double mag_n = MAG(n[icv]);
        assert(mag_n > 0.0);
        FOR_I3 n[icv][i] = n[icv][i]/mag_n;
      }
    }

    MPI_Reduce(&my_np,&np, 1, MPI_INT, MPI_SUM,0,mpi_comm);
    MPI_Reduce(&my_vofsum,&vofsum, 1, MPI_DOUBLE, MPI_SUM,0,mpi_comm);

    if (np > 0) {
      COUT1(" > total number of isolated drop = " << np);
      COUT1(" > average isolated vof = " << vofsum << " " <<  vofsum/double(np) );
    }

    // construct iband structure
    // 0 : pure gas
    // 1 : well resolved  liquid
    // 2 : first interface band next to pure liquid
    // 3 : second band next to the first band..
    // 4 : third band next to the second band...
    // 5 : fourth band next to the third band ...
    // 10: subgrid interface
    // 81 : isolated cell.


    for (int index = 1; index <= 8; ++index) {
      FOR_IFA {
        const int icv0 = cvofa[ifa][0];
        const int icv1 = cvofa[ifa][1];
        if (iband[icv0] == index && iband[icv1] > index) {
          iband[icv1] = index+1;
        }
        else if (iband[icv0] > index && iband[icv1] == index) {
          iband[icv0] = index + 1;
        }
      }
      updateCvData(iband);
    }


//    FOR_ICV {
//      if (kappa[icv] > 1.0/r_vv[icv] ) iband[icv] = 20;
//      else if (kappa[icv] < -1.0/r_vv[icv]) iband[icv] = 20;
//    }

    FOR_ICV cv_flag_real[icv] = double(iband[icv]);
    dumpRange(cv_flag_real, ncv, "iband");
  }


  void buildIband() {

    if (step%check_interval == 0 ) COUT1(" > construct band structure and detect subgrid cells");

    double eps = vof_zero;

    //updateGbandandNormal(); // reconstruct Gband and normal ...

    const int max_even_bin = 80;
    int (*iband) = new int[ncv_g];
    // bool *cv_check = new bool[ncv];

    //FOR_ICV cv_check[icv] = false;

    // pure liquid max_even_bin +1
    // pure gas max_even_bin


    FOR_ICV {
      if (vof[icv] > 1.0-vof_zero)
        iband[icv] = max_even_bin+1;
      else if (vof[icv] < vof_zero)
        iband[icv] = max_even_bin;
      else
        iband[icv] = max_even_bin+3; // pure interface 0<vof<1
    }

    updateCvData(iband);
   // confirm non-isolated interface cell
    FOR_ICV {
      if (iband[icv] == max_even_bin + 3 ) {
        for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          if ( icv != icv_nbr) {
            assert(vof[icv_nbr] < vof_zero);
          }
        }
      }
    }


    // detect non-isolated interface cell
    FOR_ICV {
      if (iband[icv] == max_even_bin + 3 ) {
        for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          assert((icv_nbr >= 0)&&(icv_nbr < ncv_g));
          if ( icv != icv_nbr) {
            if (iband[icv_nbr] != max_even_bin) {
              iband[icv] = max_even_bin-1; //non_isolated cell
              break;
            }
          }
        }
      }
    }

    updateCvData(iband);
    // confirm non-isolated interface cell
    FOR_ICV {
      if (iband[icv] == max_even_bin + 3 ) {
        for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          if ( icv != icv_nbr) {
            assert(vof[icv_nbr] < vof_zero);
          }
        }
      }
    }

    // compute number of isolated cells.
    int my_np = 0;
    int np = 0;
    double my_vofsum = 0.0;
    double vofsum = 0.0;
    FOR_ICV {
      if (iband[icv] == max_even_bin+3) {
        ++my_np;
        my_vofsum += vof[icv];
        cout << "isolated drop = " << icv << " " << vof[icv] << " " << COUT_VEC(x_cv[icv]) << " " << cv_flag[icv] <<" " << COUT_VEC(n[icv]) << endl;
        FOR_I3 n[icv][i] = -u[icv][i];
        double mag_n = MAG(n[icv]);
        assert(mag_n > 0.0);
        FOR_I3 n[icv][i] = n[icv][i]/mag_n;
       // vof[icv] = 0.0;
      }
    }

    MPI_Reduce(&my_np,&np, 1, MPI_INT, MPI_SUM,0,mpi_comm);
    MPI_Reduce(&my_vofsum,&vofsum, 1, MPI_DOUBLE, MPI_SUM,0,mpi_comm);

    if (np > 0) {
      COUT1("total number of isolated drop = " << np);
      COUT1("Average Vof sum = " << vofsum << " " <<  vofsum/double(np) );
    }
    FOR_ICV cv_flag_real[icv] = (iband[icv]);
    delete[] iband;

  }

  
  void calcVofSFromG() {

   if (mpi_rank == 0 && step%check_interval == 0 ) cout << " > calcVofSFromG() " << endl;
    FOR_ICV_G {
      const double beta_c = getDoubleParam("BETA_S",10.0);
      const double h_avg = pow((6.0*vol_cv[icv]),1.0/3.0);
      const double beta_s =  beta_c/h_avg;
      vofs[icv] = 0.5*(1.0+tanh(beta_s*(g[icv])));
      vofs[icv] = min(1.0,max(0.0,vofs[icv]));
    }

  }
*/

  void buildVofFromG() {

    FOR_ICV {
      vof[icv] = 0.5*(1.0+tanh(beta[icv]*(g[icv])));
      vof[icv] = min(1.0,max(0.0,vof[icv]));
    }

    limitVof();

  }

  virtual void loadBalanceHook(StripedMesh* sm) {
    COUT1("VofSolver::loadBalanceHook()");
    // no load balance is needed
  /*
    // parse infile...
    FOR_PARAM_MATCHING("LOAD_BALANCE") {
      string token = param->getString(0);
      if ((token == "vof")||(token == "VOF")) {
        const int vof_wgt = param->getInt(1);
        for (list<DnData>::iterator iter = sm->dnList.begin(); iter != sm->dnList.end(); ++iter) {
          if (iter->name == "vof") {
            COUT1(" > load balancing based on vof with weight: " << vof_wgt);
            if (sm->cv_lb_wgt == NULL) {
              sm->cv_lb_wgt = new int8[sm->ncv];
              for (int icv = 0; icv < sm->ncv; ++icv) sm->cv_lb_wgt[icv] = 0;
            }
            for (int icv = 0; icv < sm->ncv; ++icv) {
              if ((iter->data[icv] > vof_zero)&&(iter->data[icv] < 1.0-vof_zero)) {
                sm->cv_lb_wgt[icv] += vof_wgt;
              }
            }
            break;
          }
        }
      }
    }
*/
  }

  // lsp stuff ..
  void updateLspPrimitiveData() {
    // should not set xp0 here. it should be set only in the begining of the
    // time step not after each substep, otherwise the collision detection fails if
    // we set it equal to xp at the end of the time step
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      updateDp(ip);
    }
  }

  void updateDp(const int &ip) {
    lpCls->lp[ip].dp = pow(lpCls->lp[ip].mp/(M_PI/6.0*lpCls->lp[ip].npar*material->calcRho(lpCls->lp[ip].Tp)) , 1./3.);
  }


   void setLspX0() {
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      FOR_I3 lpCls->lp[ip].xp0[i] = lpCls->lp[ip].xp[i];
      lpCls->lp[ip].mp0 = lpCls->lp[ip].mp;
    }
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
    //0.0751, -0.0052, 0.0025};
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

  // this function takes out the particles with flag=TRASH out. You should set the
  // flag of the particles outside of this function preferencially.
  void recycleLsp(bool verbose=true) {
    // set the flag (option)
    int ntrash = lpCls->recycleLp();
    if ((step%check_interval == 0) && mpi_rank==0 && verbose)
      cout << " > recycleLsp(): # of trashed particles: " << ntrash << endl;
  }

  void setPaoCv() {

    //IF_RANK0 cout << "setPaoCv()" << endl;
    //const int MAX_PPC_viz = 300;
    int * paocv_n = new int[ncv];
    FOR_ICV paocv_n[icv] = 0;
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      assert(lpCls->lp[ip].icv < ncv);
      paocv_n[lpCls->lp[ip].icv]++;
    }
    // check the number of particles in a cell
    //FOR_ICV {
    //  if (paocv_n[icv] > MAX_PPC_viz)
    //    cout << "# Particles in this cell: " << paocv_n[icv] << " , icv: " << icv << " , mpi_rank: " << mpi_rank << endl;
    //}
    assert(paocv_i);
    paocv_v.resize(lpCls->size());
    paocv_i[0] = 0;
    for (int icv = 1; icv < ncv+1; ++icv) {
      paocv_i[icv] = paocv_i[icv-1] + paocv_n[icv-1];
      paocv_n[icv-1] = 0;
    }
    //cout << "lspVec.size(): " << lspVec.size() << " , paocv_i[ncv]: " << paocv_i[ncv] << endl;
    assert(paocv_i[ncv] == lpCls->size());
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      const int icv = lpCls->lp[ip].icv;
      paocv_v[paocv_i[icv] + paocv_n[icv]] = ip;
      paocv_n[icv]++;
    }
    delete[] paocv_n;
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
      fprintf(fp,"\"DP\"\n");
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
  bool checkCvParcel() {
    IF_RANK0 cout << "  >> checkCvParcel()";
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
        for (int poc = paocv_i[icv]; poc != paocv_i[icv+1]; ++poc) {
          const int ip = paocv_v[poc];
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
        for (int poc = paocv_i[icv]; poc != paocv_i[icv+1]; ++poc) {
          //const int ii = poc - paocv_i[icv];
          const int ip = paocv_v[poc];
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

    const double K1 = bu_k1;
    const double K2 = bu_k2;
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


  void injectNewLsp(const double time,const double dt,const bool verbose) {
    vector<LspData> lspDataVec;
    getNewLspFromInjectors(lspDataVec,time,dt,verbose);

    // we have already established in the injector the fact that 
    // all the contents of lspDataVec should be injected in this rank

    int my_npnew = lspDataVec.size();
    int npnew_global;
    MPI_Reduce(&my_npnew, &npnew_global, 1, MPI_INT, MPI_SUM, 0, mpi_comm);

    if ((mpi_rank == 0) && verbose && (step%check_interval == 0)) {
      cout << " > injecting total of " << npnew_global << " new parcels" << endl;
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

    /*
    if ( npnew_global > 0 ) {
      // figure out which particles are ours...
      // We are going to use global cv indexing to do this...

      int8 * cvora = NULL;
      buildXora(cvora,ncv);
      assert( cvora[mpi_rank+1] - cvora[mpi_rank] == ncv );

      vector<int> cvList;
      int * my_icvp_global = new int[npnew_global];

      for (int ip = 0; ip < npnew_global; ++ip) {
        my_icvp_global[ip] = -1;
        cvAdt->buildListForPoint(cvList, lspDataVec[ip].xp);

        const int icv = getClosestICVFromList(cvList,lspDataVec[ip].xp);
        if (icv >= 0) my_icvp_global[ip] = icv + cvora[mpi_rank]; // local to global cv offset

      }
      // at this point, my_icvp_global should contain the global cv index of
      // the containing cell, or -1... get the maximum across the
      // processors to break ties or check for particles that could not be located...

      int * icvp_global = new int[npnew_global];
      MPI_Allreduce(my_icvp_global,icvp_global,npnew_global,MPI_INT,MPI_MAX,mpi_comm);

      int npnew = 0;
      for (int ip = 0; ip < npnew_global; ++ip) {
        if (icvp_global[ip] == -1) {
          if (mpi_rank == 0)
            cout << "Warning: could not locate cv for particle at: " <<
              lspDataVec[ip].xp[0] << " " <<
              lspDataVec[ip].xp[1] << " " <<
              lspDataVec[ip].xp[2] << endl;
        }
        else if (icvp_global[ip] == my_icvp_global[ip]) {
          ++npnew;
        }
      }

      // resize and add particles on processors that need it...

      if (npnew > 0) {

        int npold = lpCls->size();
        lpCls->resize(npold+npnew);

        // loop through all the particles again
        int npnew_check = 0;
        for (int ip = 0; ip < npnew_global; ++ip) {
          // if the particles belongs to the processors
          if ((icvp_global[ip] != -1)&&(icvp_global[ip] == my_icvp_global[ip])) {
            int ipnew = npnew_check++;
            // local cv index
            lpCls->lp[npold+ipnew].icv = my_icvp_global[ip] - cvora[mpi_rank];
            // data...
            FOR_I3 lpCls->lp[npold+ipnew].xp[i]  = lspDataVec[ip].xp[i];
            FOR_I3 lpCls->lp[npold+ipnew].up[i]  = lspDataVec[ip].up[i];
            //FOR_I3 lpCls->lp[npold+ipnew].up[i]  = u[lpCls->lp[npold+ipnew].icv][i]; //XXX this is a hack to set the velocity as the gas velocity, shouldn't forget to change it back XXX TODO XXX TODO XXX
            lpCls->lp[npold+ipnew].mp            = lspDataVec[ip].mp;
            lpCls->lp[npold+ipnew].Tp            = lspDataVec[ip].Tp;
            // set the source terms to zero...
            lpCls->lp[npold+ipnew].tbu           = 0.0;  // breakup time
            lpCls->lp[npold+ipnew].flag          = KEEP; // recycling var
            lpCls->lp[npold+ipnew].npar          = 1.0;  // parcel info
          }
        }
        assert( npnew_check == npnew );
      }
      delete[] cvora;
      delete[] my_icvp_global;
      delete[] icvp_global;
    }
    */
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
    if (mpi_rank == 0 && step%check_interval == 0) cout << " > breakupLsp() " << endl;
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

  void getNewLspFromInjectors(vector<LspData>& lspDataVec,
            const double time,const double dt,const bool verbose) {

    lspDataVec.clear();

    for (list<LspInjector*>::iterator iter = injectorList.begin();
        iter != injectorList.end(); ++iter) {
      // no injection after passing inject_time
      if ( ((*iter)->getInjectTime() > 0) && ((time-time0) >= (*iter)->getInjectTime()) ) continue;

      (*iter)->addNewLsp(lspDataVec,time,dt,verbose);
    }

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

  }

  void initDropStats() {

    FOR_PARAM_MATCHING("DROP.PDF") {
      // add a new pdf...
      dropPDFList.push_back(DropPDF());
      DropPDF * stats = &(dropPDFList.back());
      stats->init(&(*param));
      stats->report();
    }



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


  void updateDropStats(const double* dropD, const double (*dropX)[3], const double (*dropU)[3], const int ndrops, const double time,const double dt,const int step) {

    if ((mpi_rank == 0) && (step%check_interval == 0))
      cout << "updateDropStats()" << endl;
    
    for (list<DropPDF>::iterator pdfs = dropPDFList.begin(); pdfs != dropPDFList.end(); ++pdfs) {
      if ((step%pdfs->getInterval() == 0) && (step>=pdfs->getStart())) {
        pdfs->update(dropD, dropX, dropU, ndrops, time, dt);
        pdfs->write(time,step,1);
      }
    }

  }


  void updateDropStats(const double* blob_surf, const double* blob_vol, const double (*blobX)[3], const int nblobs, const double time,const double dt,const int step,const bool verbose) {

    if ((mpi_rank == 0) && verbose && (step%check_interval == 0))
      cout << "updateDropStats()" << endl;

    for (list<DropPDF>::iterator pdfs = dropPDFList.begin(); pdfs != dropPDFList.end(); ++pdfs) {
      if ((step%pdfs->getInterval() == 0) && (step>=pdfs->getStart())) {
        pdfs->update(blob_surf, blob_vol, blobX, nblobs, time, dt,verbose);
        pdfs->write(time,step,verbose);
      }
    }

  }


  void updateDropStats(const double* blob_surf, const double* blob_vol, const double (*blobX)[3], const int nblobs, const LspState * lsp,const int np,const double time,const double dt,const int step,const bool verbose) {

    if ((mpi_rank == 0) && verbose && (step%check_interval == 0))
      cout << "updateDropStats()" << endl;

    for (list<DropPDF>::iterator pdfs = dropPDFList.begin(); pdfs != dropPDFList.end(); ++pdfs) {
      if ((step%pdfs->getInterval() == 0) && (step>=pdfs->getStart())) {
        pdfs->update(blob_surf, blob_vol, blobX, nblobs, lsp,np,material,time,dt,verbose);
        pdfs->write(time,step,verbose);
      }
    }

  }

  void updateLspStats(const LspState * lsp,const int np,const double time,const double dt,const int step,const bool verbose) {

    if ((mpi_rank == 0) && verbose && (step%check_interval == 0))
      cout << "updateLspStats()" << endl;

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

  void calcVolMassFrac() {

    double max_vol_r = -1e20;
    double min_vol_r = 1e20;
    double max_mas_r = -1e20;
    double min_mas_r = 1e20;
    FOR_ICV {
      double sum_vol = 0.0;
      double sum_mas = 0.0;
      for (int ii = paocv_i[icv]; ii != paocv_i[icv+1]; ++ii) {
        int ip = paocv_v[ii];
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
      cout << "particle volume ratio " << min_VR <<  ":" << max_VR << endl;
      cout << "particle mass ratio " << min_MR <<  ":" << max_MR << endl;
    }
  }

  //template<class T>
  //void relocateLsp() {

  //  if ((mpi_rank == 0) && (step%check_interval == 0))
  //    cout << "relocateLsp";

  //  //XXX
  //  //assert(surface);

  //  double tol = 0.0;

  //  int * send_count = new int[mpi_size];
  //  int * send_disp = new int[mpi_size];
  //  int * recv_count = new int[mpi_size];
  //  int * recv_disp = new int[mpi_size];

  //  int done = 0;
  //  while (done == 0) {

  //    if ((mpi_rank == 0) && (step%check_interval == 0)) {
  //      cout << ".";
  //      cout.flush();
  //    }

  //    FOR_RANK send_count[rank] = 0; // the number of particles we need to send to each rank

  //    int my_done = 1;

  //    // hack to get closest... function should return ip_debug for the right rank and for the
  //    // other ranks will return -1
  //    //double xp_debug[3] = { 0.0752, -0.0049, 0.003 };
  //    //double xp_debug[3] = { 0.0739, -0.005, 0.0038 };
  //    //int ip_debug = searchClosestIpX(xp_debug);

  //    for (int ip = 0; ip < lpCls->size(); ++ip) {

  //      //bool debug = false;
  //      //if (ip == ip_debug && step >144600 && step <= 144649)
  //      //if (ip == ip_debug && step >145900 && step <= 146100)
  //      //  debug = false;

  //      // this particle is currently in cv "icv"...
  //      // we use -ve icv indiexing top indicate that a particle has been
  //      // properly located, and does not need to be revisited...

  //      const int icv = lpCls->lp[ip].icv;
  //      if (icv >= 0) {

  //        if (icv >= ncv)
  //          cout << "XXXXXX: " << icv << " " << ncv << endl;
  //        assert(icv < ncv);

  //        // this particle still needs to be checked...
  //        // and has followed a trajectory from xp0 to xp...

  //        double * xp0 = lpCls->lp[ip].xp0;
  //        double * xp  = lpCls->lp[ip].xp;

  //        //if (debug) cout << "debug: " << " ,ip: " << ip << " ,rank: " << mpi_rank << " xp0: " << COUT_VEC(xp0) << " xp: " << COUT_VEC(xp) << endl;

  //        // if this cell has boundaries, then check for crossing...

  //        double xi_closest[3];
  //        double st_normal_closest[3];
  //        int ist_closest = -1;
  //        double d2_closest;
  //        int zone_closest;
  //        for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
  //          const int ibf = bfocv_v[boc];
  //          for (int sob = sstobf_i[ibf]; sob != sstobf_i[ibf+1]; ++sob) {
  //            const int ist = sstobf_v[sob];
  //            // the corner coords of this surface tri are...
  //            const double * const x0 = subSurface->xp[subSurface->spost[ist][0]];
  //            const double * const x1 = subSurface->xp[subSurface->spost[ist][1]];
  //            const double * const x2 = subSurface->xp[subSurface->spost[ist][2]];
  //            // did we cross this tri?...
  //            const double st_normal[3] = TRI_NORMAL_2(x0,x1,x2); // should be an outward-pointing normal (wrt flow)
  //            // get the side associated with both points...
  //            const double d0 =
  //              (xp0[0]-x1[0])*st_normal[0] +
  //              (xp0[1]-x1[1])*st_normal[1] +
  //              (xp0[2]-x1[2])*st_normal[2];
  //            const double d1 =
  //              (xp[0]-x1[0])*st_normal[0] +
  //              (xp[1]-x1[1])*st_normal[1] +
  //              (xp[2]-x1[2])*st_normal[2];
  //            // NOTE: this requires tolerance
  //            if ((d0 <= 0)&&(d1 > 0)) {
  //              // this is a crossing!...
  //              //if (debug) cout << "GOT POTENTIAL CROSSING! icv: " << icv << " d0: " << d0 << " d1: " << d1 << " st_normal: " << COUT_VEC(st_normal) << " xp: " << COUT_VEC(xp) << " step: " << step << endl;
  //              // compute xi...
  //              double xi[3]; FOR_I3 xi[i] = (d1*xp0[i] - d0*xp[i])/(d1-d0);
  //              const double n01[3] = TRI_NORMAL_2(x0,x1,xi);
  //              //if (debug) cout << "n01: " << DOT_PRODUCT(st_normal,n01) << " step: " << step << endl;
  //              if (DOT_PRODUCT(st_normal,n01) >= 0.0) {
  //                const double n12[3] = TRI_NORMAL_2(x1,x2,xi);
  //                //if (debug) cout << "n12: " << DOT_PRODUCT(st_normal,n12) << " step: " << step << endl;
  //                if (DOT_PRODUCT(st_normal,n12) >= 0.0) {
  //                  const double n20[3] = TRI_NORMAL_2(x2,x0,xi);
  //                  //if (debug) cout << "n20: " << DOT_PRODUCT(st_normal,n20) << " step: " << step << endl;
  //                  if (DOT_PRODUCT(st_normal,n20) >= 0.0) {
  //                    // xi is inside the tri...
  //                    // we want the point that is closes to xp0 -- i.e. the starting point of the particle.
  //                    const double d2 = DIST2(xp0,xi);
  //                    //if (debug) cout << "THIS WAS AN ACTUAL CROSSING with dist: " << sqrt(d2) << endl;
  //                    if ((ist_closest == -1)||(d2 < d2_closest)) {
  //                      d2_closest = d2;
  //                      ist_closest = ist;
  //                      FOR_I3 xi_closest[i] = xi[i];
  //                      FOR_I3 st_normal_closest[i] = st_normal[i];
  //                      zone_closest = zone_bf[ibf];
  //                      assert(zone_closest < bfZoneVec.size());
  //                    }
  //                  }
  //                }
  //              }
  //              //getchar();
  //            }
  //          }
  //        }

  //        // if we got a closest intersection, check if it is closest to us, or
  //        // instead closer to one of our neighbors...

  //        if (ist_closest != -1) {
  //          double d2_min = DIST2(xi_closest,x_vv[icv]);
  //          int icv_min = icv;
  //          // check if there is a nbr in which we are closer...
  //          assert(cvocv_v[cvocv_i[icv]] == icv); // self check
  //          for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) { // note we skip self here
  //            const int icv_nbr = cvocv_v[coc];
  //            const double d2_icv_nbr = DIST2(xi_closest,x_vv[icv_nbr]);
  //            if (d2_icv_nbr < d2_min) {
  //              d2_min = d2_icv_nbr;
  //              icv_min = icv_nbr;
  //            }
  //          }
  //          if (icv_min == icv) {
  //            // apply the bc...
  //            //if (debug) cout << "REFLECTION: " << " ip: " << ip << " step: " << step << endl;
  //            double buf[7];
  //            lpZoneBcVec[zone_closest]->applyBc(lpCls->lp[ip], xi_closest, st_normal_closest, buf);
  //          }
  //          else {
  //            lpCls->lp[ip].icv = icv_min;
  //            if (icv_min >= ncv) {
  //              // icv_min is in the ghost range, so we will be sending this lsp to another rank...
  //              int rank,bits,index;
  //              BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_min-ncv]);
  //              ++send_count[rank];
  //            }
  //          }
  //          // in either case we are not done...
  //          my_done = 0;
  //        }
  //        else {

  //          // step 1: loop through any surface tris that are part of this
  //          // cv's boundary and find the closest tri intersection: xi. Note
  //          // that because the tris are the original surface tris, they may
  //          // extend beyond the icv's local voronoi diagram.
  //          // If there is no intersection, then just continue to the nbr check below.
  //          // If there is an intersection, run the nbr check below on the
  //          // point xi. If the intersection point is closest to icv, simply apply the bc
  //          // associated with this surface tri (e.g. update position) and continue
  //          // to point check. If xi is closer to an icv_nbr, send the particle to that nbr,
  //          // and releat everything, including the calculation of nearest xi (because it
  //          // might be different due to additional surface geometry in the nbr)...

  //          // calc the current distance to this point...
  //          double d2_min = DIST2(lpCls->lp[ip].xp,x_vv[icv]);
  //          int icv_min = icv;
  //          // check if there is a nbr in which we are closer...
  //          assert(cvocv_v[cvocv_i[icv]] == icv); // self check
  //          for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) { // note we skip self here
  //            const int icv_nbr = cvocv_v[coc];
  //            const double d2_icv_nbr = DIST2(lpCls->lp[ip].xp,x_vv[icv_nbr]);
  //            if (d2_icv_nbr < d2_min) {
  //              d2_min = d2_icv_nbr;
  //              icv_min = icv_nbr;
  //            }
  //          }
  //          if (icv_min != icv) {
  //            lpCls->lp[ip].icv = icv_min;
  //            if (icv_min >= ncv) {
  //              // icv_min is in the ghost range, so we will be sending this lsp to another rank...
  //              int rank,bits,index;
  //              BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_min-ncv]);
  //              ++send_count[rank];
  //            }
  //            my_done = 0;
  //          }
  //          else {
  //            // we usea negative icv to indicate that this particle is done. i.e.
  //            // it is in the right place and doe NOT need to be visited again.
  //            lpCls->lp[ip].icv = -lpCls->lp[ip].icv-1; // flip icv to negative...
  //          }
  //        }
  //      }
  //    } // for (int ip...

  //    // redistibute the lpVec for elements that have passed into ghost cvs...
  //    // note: we could also remove elements here in the future if required
  //    // (e.g. evaporated particles). Not sure if theis is the right place...

  //    send_disp[0] = 0;
  //    for (int rank = 1; rank < mpi_size; ++rank)
  //      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
  //    const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

  //    int * send_buf_int = new int[send_count_sum]; // for the icv's on the recv'ing rank
  //    double * send_buf_double = new double[send_count_sum*T::data_size()];

  //    // now pack the buffer and compress lpVec...

  //    int np_new = 0;
  //    //double R[9], t[3];
  //    for (int ip = 0; ip < lpCls->size(); ++ip) {
  //      if (lpCls->lp[ip].icv >= ncv) {
  //        int rank,bits,index;
  //        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[lpCls->lp[ip].icv-ncv]);
  //        send_buf_int[send_disp[rank]] = index; // the icv on rank where we are to unpack
  //        //lpCls->lp[ip].pack(send_buf_double+send_disp[rank]*T::data_size());

  //        const int ipR = iRt_g[lpCls->lp[ip].icv-ncv][0];
  //        const int ipt = iRt_g[lpCls->lp[ip].icv-ncv][1];
  //        if ((ipR >= 0) && (ipt >= 0)) {
  //          // rotate/translate xp/xp0 and rotate up
  //          lpCls->lp[ip].packRt(send_buf_double+send_disp[rank]*T::data_size(),per_R[ipR],per_t[ipt]);
  //        }
  //        else if (ipR >= 0) {
  //          // rotate xp/xp0 and rotate up
  //          lpCls->lp[ip].packR(send_buf_double+send_disp[rank]*T::data_size(),per_R[ipR]);
  //        }
  //        else if (ipt >= 0) {
  //          // translate xp/xp0
  //          lpCls->lp[ip].packt(send_buf_double+send_disp[rank]*T::data_size(),per_t[ipt]);
  //        }
  //        else {
  //          lpCls->lp[ip].pack(send_buf_double+send_disp[rank]*T::data_size());
  //        }
  //        ++send_disp[rank];
  //      }
  //      else {
  //        if (np_new != ip)
  //          lpCls->lp[np_new].copy(lpCls->lp[ip]);
  //        ++np_new;
  //      }
  //    }

  //    // figure out how many we will be receiving and make space for them...

  //    send_disp[0] = 0;
  //    for (int rank = 1; rank < mpi_size; ++rank)
  //      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

  //    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

  //    recv_disp[0] = 0;
  //    for (int rank = 1; rank < mpi_size; ++rank)
  //      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
  //    const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

  //    int * recv_buf_int = new int[recv_count_sum];
  //    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
  //      recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
  //    delete[] send_buf_int;

  //    FOR_RANK {
  //      send_count[rank] *= T::data_size();
  //      send_disp[rank] *= T::data_size();
  //      recv_count[rank] *= T::data_size();
  //      recv_disp[rank] *= T::data_size();
  //    }

  //    double * recv_buf_double = new double[recv_count_sum*T::data_size()];
  //    MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
  //      recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
  //    delete[] send_buf_double;

  //    // make room for data in lpVec...

  //    lpCls->resize(np_new+recv_count_sum);

  //    // unpack...

  //    for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
  //      lpCls->lp[np_new].icv = recv_buf_int[irecv]; assert((lpCls->lp[np_new].icv >= 0)&&(lpCls->lp[np_new].icv < ncv));
  //      lpCls->lp[np_new].unpack(recv_buf_double+irecv*T::data_size());
  //      // >>> no need to set this because we are unpacking flag too
  //      //lpVec[np_new].flag = KEEP;
  //      ++np_new;
  //    }

  //    delete[] recv_buf_int;
  //    delete[] recv_buf_double;

  //    // check if we are done...

  //    MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);

  //  } // while loop

  //  // check that everyone is -ve, and reset to positive...

  //  for (int ip = 0; ip < lpCls->size(); ++ip) {
  //    assert(lpCls->lp[ip].icv < 0);
  //    lpCls->lp[ip].icv = -lpCls->lp[ip].icv-1;
  //  }

  //  delete[] send_count;
  //  delete[] send_disp;
  //  delete[] recv_count;
  //  delete[] recv_disp;

  //  if ((mpi_rank == 0) && (step%check_interval == 0))
  //    cout << "OK" << endl;

  //  //MPI_Pause("ALL FINISHED");
  //}

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

  //virtual void impactLpBc(LpState& lp, const double* xi, const double* st_normal) {
  //  const double dx[3] = DIFF(lp.xp,xi);
  //  const double dp = DOT_PRODUCT(st_normal,dx); assert(dp > 0.0);
  //  const double nmag2 = DOT_PRODUCT(st_normal,st_normal);
  //  FOR_I3 lp.xp[i] -= 2.0*dp*st_normal[i]/nmag2;
  //  // and update xp0 to the xi_closest, because the remaining part of the
  //  // trajectory may have additional intersections...
  //  FOR_I3 lp.xp0[i] = xi[i];
  //  const double un = DOT_PRODUCT(st_normal,lp.up);
  //  FOR_I3 lp.up[i] -= 2.0*un*st_normal[i]/nmag2;
  //}

  //virtual void impactLpBc(LspState& lp, const double* xi, const double* st_normal, const double& tol) {
  //  const double dx[3] = DIFF(lp.xp,xi);
  //  const double dp = DOT_PRODUCT(st_normal,dx); assert(dp > -tol);
  //  const double nmag2 = DOT_PRODUCT(st_normal,st_normal);
  //  FOR_I3 lp.xp[i] -= 2.0*dp*st_normal[i]/nmag2;
  //  // and update xp0 to the xi_closest, because the remaining part of the
  //  // trajectory may have additional intersections...
  //  FOR_I3 lp.xp0[i] = xi[i];
  //  const double un = DOT_PRODUCT(st_normal,lp.up);
  //  FOR_I3 lp.up[i] -= 2.0*un*st_normal[i]/nmag2;
  //}

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

  void advanceLsp(const double time,const double dt,const int step,const bool verbose) {

    if ((mpi_rank == 0)&&verbose)
      cout << "advanceLsp()" << endl;

    assert( lpCls != NULL );


    //========================
    // inject lsp particles
    //========================
    injectNewLsp(time,dt,1);
    updateLspPrimitiveData();

    switch(BreakupModel) {
      case NO_BREAKUP:
        break;
      case SBM_BREAKUP:
        breakupLsp();
        break;
      default:
        assert(0);
    }

    setPaoCv();     // we need paocv for rk3Step_lsp_gas
    bool agg = false;
    if (lspAgglem) {
      agg = checkCvParcel(); // the ranks which do agglemorate will need to reset the paocv later
      recycleLsp();
    }

    setLspX0();           // only at the beginning of time step... should be after agglemoration routine...
    if (agg) setPaoCv();  // if agg then particles are recycled -> setPaoCv() is needed

    // now we have and ordered recycled particle list
    // advance particles in time, set the mass to zero for particles that evaporate
    // completely in this time step, and update the particle source terms

    solveODEandUpdateSources(dt,verbose);
    updateLspPrimitiveData();


    // report evaporated mass before recycling, change of size of lspVec, and relocating
    if (LIQ_LP && evaporation)
      reportLspEvap();

    //relocateLsp<LspState>();
    relocateLsp_rayTrace<LspState>();
    recycleLsp(); // after this PaoCv_i/v is not current, if you need them you should call setPaoCv()

    // this virtual routine allows the inheriting class to set the source terms
    // for 2-way coupling. Here we call this BEFORE relocating the particles, so we
    // use the existing, valid cvs...
    //reduceLspSourcesToGasPhase();


  }

  void solveODEandUpdateSources(const double dt,const bool verbose) {

    // evaporation model based on Langmuir-Knudsen model
    // see Miller et al. Int J Multiphase Flow, 24:1025-1055, 1998.

    // NOTE: "XXXX" terms should eventually come from solver/chemtable

    int np = lpCls->size();
    int mynpevap      = 0;
    double mymassevap = 0.0;
    const double nearOne = .9999999999;
    const double SMALL = 1.0E-15;
    double dp_tol;
    dp_tol =getDoubleParam("LSP.DP_TOL",1.0E-8);

    // initialize the source terms
    FOR_ICV_G lspMassSource[icv] = 0.0;
    FOR_ICV_G FOR_I3 lspMomentumSource[icv][i] = 0.0;

    // get the gas-phase properties from the solver
    const double Tr    = fuel->Tboil;;
    const double Tg  = T_ref;
    const double Pg   = p_ref;
    const double rhog = rho_g;
    const double mug  = mu_g;
    const double cpg  = R_gas*gamma/(gamma-1.0);
    const double kg   = 26.2; //[W/K] air thermal conductivity at 20 degree
    //const double Yfg  = max(0.0,cv[icv].Y);
    const double Yfg  = 0.0;
    const double MWg  = 28.97; // kg/m^3 for air // hack only for air???

    // loop over all particles

    for (int ip = 0; ip < np; ip++) {

      assert(lpCls->lp[ip].icv>=0); // assert that there are no recycled parcels
      assert(lpCls->lp[ip].flag == KEEP);       // flag to designate evaporated particles
      assert(lpCls->lp[ip].mp > 0.0); // mass is positive ...

      // pull the velocity (no interpolation here ...)
      const int icv = lpCls->lp[ip].icv;
      double ug[3] = {0.0,0.0,0.0}; // gas velocity

      // interpolate gas velocity ...
      double sum_weights = 0.0;
      for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        double d2 = 0.0; FOR_I3 { double dx = x_cv[icv_nbr][i] - lpCls->lp[ip].xp[i]; d2 += dx*dx; }
        double this_weight = 1.0 / sqrt(d2+SMALL);
        sum_weights += this_weight;
        FOR_I3 ug[i] += this_weight*u[icv_nbr][i];
      }

      FOR_I3 ug[i] /= sum_weights;

      // compute liquid properties at parcel temperature
      const double rhof = fuel->calcRho(lpCls->lp[ip].Tp);  // liquid density
      const double cpf  = fuel->calcCp(lpCls->lp[ip].Tp);   // liquid heat capacity
      const double Hvap = fuel->calcHvap(lpCls->lp[ip].Tp); // latent heat of evaporation
      const double du[3]  = DIFF(lpCls->lp[ip].up, ug);  // relative velocity
      const double du_mag = sqrt(DOT_PRODUCT(du,du));  // relative velocity magnitude
      const double Dp   = lpCls->lp[ip].dp;
      assert(Dp > 0.0);

      // equilibrium fuel composition at droplet surface
      //Yfg =  min(nearOne,Yfg);
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


      // compute non-dimensional and model parameters
      const double tau = rhof*Dp*Dp/(18.*mur);
      const double Rep = rhor*Dp*du_mag/mur;
      const double Prg = mur*cpr/kr;
      const double Scg = mur/(rhor*fuel->calcDv(Tr));
      const double Nug = 2. + 0.552*sqrt(Rep)*pow(Prg,1./3.);
      const double Shg  = 2. + 0.552*sqrt(Rep)*pow(Scg,1./3.);
      const double Lk  = mur*sqrt(2.0*M_PI*lpCls->lp[ip].Tp*fuel->R_UNIVERSAL/fuel->MW)/(Scg*Pg);


      // Newton-iteration for evaporation rate mpdot
      const double F1 = (Shg/Scg/3.0)*(lpCls->lp[ip].mp/tau); // mpdot + F1*ln(1+BM) = 0, see eq (4)
      const double F2 = -1.5*Prg*tau/lpCls->lp[ip].mp;        // beta = F2*mpdot, see eq (17)
      double beta, Xfsneq, Yfsneq, BMneq, Fm;
      double mpdot  = -1.0E-08;
      double mpdot0 =  0.0;
      double Fm0    =  mpdot0 + F1*log(1.0+BM); // initialize with equilib value
      double eps    =  1.0E-15;
      double resid  =  eps+1.0;
      int    iter   =  0;
      if (evaporation) {
        while (resid > eps) {
          iter += 1;
          beta   = F2*mpdot;
          Xfsneq = Xfseq - (Lk/Dp*2.0)*beta;
          Yfsneq = Xfsneq/(Xfsneq+(1.0-Xfsneq)*MWg/fuel->MW);
          BMneq  = max(-nearOne,(Yfsneq-Yfg)/(1.0-Yfsneq));
          Fm     = mpdot + F1*log(1.0+BMneq);
          double tmp = mpdot;
          if (fabs(Fm-Fm0) < eps ) break;
          mpdot = mpdot - Fm/(Fm-Fm0)*(mpdot-mpdot0);
          resid = sqrt(pow(Fm-Fm0,2));
          Fm0 = Fm;
          mpdot0 = tmp;
          if (iter > 20) cout << "ERROR: LSP Newton iteration did not converge !!!" << endl;
        }
        beta = F2*mpdot;
      }else {
        mpdot = 0.0;
        beta = 0.0;
      }

      // integrate droplet equations
      //   dxpdt = up
      //   dupdt = C1*(up-ug) + C2
      //   dmpdt = mpdot
      //   dTpdt = C4*(Tp-Tg) + C5
      double up0[3]; FOR_I3 up0[i] = lpCls->lp[ip].up[i];
      double mp0 = lpCls->lp[ip].mp;
      double Tp0 = lpCls->lp[ip].Tp;

      // drag law from  Crowe et al, Multiphase Flows w/ Droplets & Particles, 1998
      double f1 = 1.0+0.15*pow(Rep,0.687) + 0.0175*Rep/(1.0+4.25e4*pow(Rep,-1.16));

      // solve u first order implicit
      double C1 = -f1/tau;
      double C2[3] = { 0.0, 0.0, 0.0 }; // placeholder for gravity or other body force
      FOR_I3 lpCls->lp[ip].up[i] = ( dt*(C1*ug[i]-C2[i]) - lpCls->lp[ip].up[i] )/( C1*dt-1.0 );

      // solve xp using Crank-Nicolson
      FOR_I3 lpCls->lp[ip].xp[i] += dt*0.5*( up0[i] + lpCls->lp[ip].up[i] );

      // update mp using mpdot from Newton iteration above (don't allow negative mass)
      lpCls->lp[ip].mp = max( 0.0, lpCls->lp[ip].mp+dt*mpdot );

      if (evaporation) {
        // Tp using Crank-Nicolson
        double f2 =  beta/(exp(beta)-1.0);                   // Nusselt correction, see eq (19)
        double C4 = -f2*Nug/Prg/3.0*(cpr/cpf)/tau;           // first term in eq (3)
        double C5 =  Hvap*mpdot/(cpf*0.5*(mp0+lpCls->lp[ip].mp)); // note: using mp0 and mp
        lpCls->lp[ip].Tp = ( dt*(C4*Tg-C5) - lpCls->lp[ip].Tp*(1.0+dt*C4/2.0) )/( C4*dt/2.0-1.0 );

        // limit Tp to Tboil and prevent 2nd law violations
        lpCls->lp[ip].Tp = min( nearOne*fuel->Tboil, lpCls->lp[ip].Tp );
        lpCls->lp[ip].Tp = max(     min(Tp0,Tg), lpCls->lp[ip].Tp );

        // compute new particle size, recycle if dp < threshold
        double dp = pow(lpCls->lp[ip].mp/(M_PI/6.0*lpCls->lp[ip].npar*fuel->calcRho(lpCls->lp[ip].Tp)),1.0/3.0);
        if ( dp<=dp_tol ) {
          lpCls->lp[ip].mp   = 0.0;
          lpCls->lp[ip].flag = TRASH;
          ++mynpevap;
        }
      }

      // update the source terms for each particle
      // Force on gas-phase due to particle = -d/dt(mp*up) = -(mp_t*up_t - mp_0*up_0)/dt
      // Mass evaporation rate = d/dt(mp) = (mp_t - mp_0)/dt
      // Mass source due to evaporation = -d/dt(mp)
      // This should ensure mass/momentum conservation. Note that we do not
      // divide by dt here, just in case the time step is varying (e.g. the case of constant cfl).
      double msg = -(lpCls->lp[ip].mp - mp0);
      mymassevap +=msg;
      double fsg[3];
      FOR_I3 fsg[i] = -(lpCls->lp[ip].mp*lpCls->lp[ip].up[i] - mp0*up0[i]);

      // first compute the total weight
      double totalweight = 0.0;
      for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        // for now, just a simple inverse distance averaging ...
        double d2 = 0.0; FOR_I3 { double dx = x_cv[icv_nbr][i] - lpCls->lp[ip].xp[i]; d2 += dx*dx; }
        double weight = 1.0 / sqrt(d2+SMALL);
        totalweight += weight;
      }

      double inv_totalweight = 1.0/totalweight;
      for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        double d2 = 0.0; FOR_I3 { double dx = x_cv[icv_nbr][i] - lpCls->lp[ip].xp[i]; d2 += dx*dx; }
        double weight = 1.0 / sqrt(d2+SMALL);
        double coeff = weight*inv_totalweight;
        lspMassSource[icv_nbr] += msg*coeff;
        FOR_I3 lspMomentumSource[icv_nbr][i] += fsg[i]*coeff;
      }

    } //FOR_IP loop

    updateCvData(lspMassSource);
    updateCvData(lspMomentumSource);

    FOR_ICV_G lspMassSource[icv] /= dt*vol_cv[icv];
    FOR_ICV_G FOR_I3 lspMomentumSource[icv][i] /= dt*vol_cv[icv];

    if (verbose && evaporation) {

      // number of evaporated particles
      int npevap = 0;
      MPI_Reduce(&mynpevap,&npevap,1,MPI_INT,MPI_SUM,0,mpi_comm);

      // total evaporated mass
      double  massevap = 0.0;
      MPI_Reduce(&mymassevap,&massevap,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

      IF_RANK0 {
        cout << " -------------------------------------------------" << endl;
        cout << " Lsp ODE solver" << endl;
        cout << " -------------------------------------------------" << endl;
        cout << " > total " << npevap << " particles have entirely evaporated" << endl;
        cout << " > total evaporated mass = " << massevap << endl;
      }
    }

  }

  // custom variable output via cti var evaluation; recall on completion
  // if the data is not found, then it must be passed down to the base
  // StaticSolver class to see if it can evaluate the data or otherwise error

  virtual CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,const string& name,
                                                    list<CtiRegister::CtiData>& args, const bool b_eval_func) {

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
    
    // this solver doesnt know how to evaluate this data -- go to FlowSolver
    return FlowSolver::funcEvalCtiData(v,name,args,b_eval_func);
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

    double eps = 1.0E-6*(dp_max-dp_min);
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


};

#undef FOR_BCZONE
#endif
