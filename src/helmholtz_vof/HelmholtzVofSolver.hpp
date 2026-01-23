#ifndef HELMHOLTZVOFSOLVER_HPP
#define HELMHOLTZVOFSOLVER_HPP

#include "FlowSolver.hpp"
#include "HelmholtzVofSolverBcs.hpp"
#include "LiquidFlux.hpp"
#include "TrilinosInterface.hpp"
#include "AverageOperator.hpp"
#include "Adt.hpp"
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
  BOX,
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

class HelmholtzVofSolver : public FlowSolver {

public:

  double gravity[3];
  double sigma;       // surface tension coefficient
  double gamma;

  double *vof;           // volume fraction
  double *vofold;           // volume fraction
  double (*u)[3];        // velocity
  double (*uold)[3];    // velocity
  double *p;             // pressure
  double *pa;            // advection pressure
  double *p_eq;          // equilibrium pressure
  double *rho0;          // liquid phasedensity
  double *rho1;          // gas phase density
  double *rho0old;
  double *rho1old;
  double *rhoY0;
  double *rho;
  double *rhoY0old;
  double *rhoold;
  double *mu_lam;        // laminar viscosity
  double *mu_sgs;        // sgs eddy viscosity
  double *kappa;         // interface curvature
  double (*n)[3];        // interface normal
  double *g;             // interface distance g : n.(x_interface-x_vv)=g
  double (*sp_dpdx)[3];    // density-weighted cv pressure gradient
  double *lsd;         // liquid surface density 

  // varaibles for compressible flow
  double *inv_K; // inverse of bulk modulus (compressiblity)

  // entropy viscous flux coefficents
  double *fgr;
  double *fgr0;
  double *fgr1;

  double *cfl2;

  // level-set stuff...
  double *sd;           // signed distance for curvature
  double d_max;         // max allowable distance (for signed distance initialization)...

  // THINC sharpness parameter ...
  double *beta;
  double beta_c;
  int rk_scheme;
  double (*rk_wgt)[3];
  double *rk_mom;
  int *iband;
  bool incomp;

  double *div;          // velocity divergence...
  double (*dudx)[3][3]; // velocity gradient
  double *q_fa;         // conservative face velocity...
  //double *q0_fa;       //  old conservative face velocity...

  // 2-level plic area-weighted paraboloid centroid fit...
  list<double> * plicPointList; // 3 x number of edge intersections b/w PLIC and cv (ordered)
  double (*plic_xc)[3]; // plic barycenter
  double *plic_area; // plic area


  int *cv_flag;      // flag interface cells
  int nftb; // number of faces touching boundary
  int *fa_ftb; // face index for a face touching boundary (not in fpnpocv_i/v)

  string sgs_model; // subgrid scale model name
  int    p_maxiter; // max pressure iterations
  double p_zero;    // pressure zero
  int    u_maxiter; // max velocity iterations
  double u_zero;    // velocity zero
  double u_relax;   // relaxation for jacobi solver
  double p_relax;   // relaxation for jacobi solver pressure
  double vof_zero;  // min geometrically-supported volume fraction
  double mass_zero; // min mass allowed in a cell
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

  
  // EOS stuff ...
  double mu0_ref, mu1_ref;
  double lambda0, lambda1;
  double rho0_ref, rho1_ref;
  double cp0, cp1;
  double gamma0, gamma1;
  double Pinf0, Pinf1;
  double p_ref;
  double T_ref;
  double q0,q1;
  double e0_ref;
  double e1_ref;
  double s0_ref;
  double s1_ref;
  double R_gas;
  double sos0;
  double sos1;

  double *rhs_rhoY0;
  double *rhs_rhoY1;

  list<LiquidFlux> liquidFluxList; // liquid flux list

  HelmholtzVofSolver() {

    COUT1("HelmholtzVofSolver()");

    // field data registration...

    // we only read need u and vof to restart calc...

    vof      = NULL; registerCvData(vof      , "vof"      , READWRITE_DATA);
    vofold   = NULL; registerCvData(vofold   , "vofold"   , CAN_WRITE_DATA);
    u        = NULL; registerCvData(u        , "u"        , READWRITE_DATA);
    uold     = NULL; registerCvData(uold     , "uold"     , CAN_WRITE_DATA);
    p        = NULL; registerCvData(p        , "p"        , READWRITE_DATA);
    pa       = NULL; registerCvData(pa       , "pa"       , CAN_WRITE_DATA);
    p_eq     = NULL; registerCvData(p_eq     , "p_eq"     , CAN_WRITE_DATA);
    rho0     = NULL; registerCvData(rho0     , "rho0"     , READWRITE_DATA);
    rho1     = NULL; registerCvData(rho1     , "rho1"     , READWRITE_DATA);
    rho0old  = NULL; registerCvData(rho0old  , "rho0old"  , CAN_WRITE_DATA);
    rho1old  = NULL; registerCvData(rho1old  , "rho1old"  , CAN_WRITE_DATA);
    rhoY0    = NULL; registerCvData(rhoY0    , "rhoY0"    , READWRITE_DATA);
    rho      = NULL; registerCvData(rho      , "rho"      , READWRITE_DATA);
    rhoY0old = NULL; registerCvData(rhoY0old , "rhoY0old" , CAN_WRITE_DATA);
    rhoold   = NULL; registerCvData(rhoold   , "rhoold" , CAN_WRITE_DATA);
    mu_lam   = NULL; registerCvData(mu_lam   , "mu_lam"   , CAN_WRITE_DATA);
    mu_sgs   = NULL; registerCvData(mu_sgs   , "mu_sgs"   , CAN_WRITE_DATA);
    kappa    = NULL; registerCvData(kappa    , "kappa"    , CAN_WRITE_DATA);
    n        = NULL; registerCvData(n        , "n"        , CAN_WRITE_DATA);
    g        = NULL; registerCvData(g        , "g"        , CAN_WRITE_DATA);
    q_fa     = NULL; registerSignedFaData(q_fa     , "q_fa"     , READWRITE_DATA);
    //q0_fa    = NULL; registerSignedFaData(q0_fa    , "q0_fa"    , CAN_WRITE_DATA);
    sd       = NULL; registerCvData(sd       , "sd"       , NO_READWRITE_DATA);
    lsd      = NULL; registerCvData(lsd      , "lsd"       , READWRITE_DATA);

    inv_K    =NULL; registerCvData(inv_K      , "inv_K"      , CAN_WRITE_DATA);

    fgr     =NULL; registerCvData(fgr       , "fgr"       , CAN_WRITE_DATA);
    fgr0    =NULL; registerCvData(fgr0      , "fgr0"      , CAN_WRITE_DATA);
    fgr1    =NULL; registerCvData(fgr1      , "fgr1"      , CAN_WRITE_DATA);
    cfl2     =NULL; registerCvData(cfl2      , "cfl2"      , CAN_WRITE_DATA);

    rhs_rhoY0 = NULL; registerCvData(rhs_rhoY0 , "rhs_rhoY0" , CAN_WRITE_DATA);
    rhs_rhoY1 = NULL; registerCvData(rhs_rhoY1 , "rhs_rhoY1" , CAN_WRITE_DATA);
    
    
    cv_flag      = NULL; registerCvData(cv_flag      , "cv_flag"      , READWRITE_DATA);

    dudx = NULL;
    sp_dpdx  = NULL;   registerCvData(sp_dpdx   , "sp_dpdx"    , READWRITE_DATA);
    div = NULL;


    plicPointList = NULL;
    plic_xc = NULL;
    plic_area = NULL;

    //cv_flag = NULL;

    nftb = 0;
    fa_ftb = NULL;

    icv_global_g = NULL;
    trilinosSolverForPressure = NULL;

    //hyperbolic tangent sharpness parameter
    beta = NULL;
    iband = NULL;

    // moved global parameters in constructor so that we had access to them earlier than during init...

    sigma = getDoubleParam("SIGMA");

    COUT2(" > SIGMA: " << sigma);

    T_ref = getDoubleParam("T_REF", 290.0);

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
    u_zero = getDoubleParam("U_ZERO",1.0E-10);
    u_relax = getDoubleParam("U_RELAX", 0.7);
    p_relax = getDoubleParam("P_RELAX", 0.8);
    u_maxiter = getIntParam("U_MAXITER",1000);
    vof_zero = getDoubleParam("VOF_ZERO",1.0E-8);
    mass_zero = getDoubleParam("MASS_ZERO", 1.0E-12);
    use_xp = getBoolParam("USE_XP",true);
    g_maxiter = getIntParam("G_MAXITER",100);

    parseSgsModel(sgs_model); // TODO may need to promote for models that need registered data

    // need more than one iterations to get second-order accuracy (1-3 seems to be good)...
    normal_iters = getIntParam("NORMAL_ITERS",2);
    kappa_iters = getIntParam("KAPPA_ITERS", 2);

    rk_scheme = getIntParam("RK_SCHEME",2);
    assert(rk_scheme <= 3);
    rk_wgt = new double[rk_scheme][3];
    rk_mom = new double[rk_scheme];

    if (rk_scheme == 1 ) {
      rk_wgt[0][0] = 0.0;
      rk_wgt[0][1] = 1.0;
      rk_wgt[0][2] = 1.0;
      rk_mom[0] = 1.0;
    }
    else if (rk_scheme == 2) {
      rk_wgt[0][0] = 0.0;
      rk_wgt[0][1] = 1.0;
      rk_wgt[0][2] = 1.0;
      rk_wgt[1][0] = 0.5;
      rk_wgt[1][1] = 0.5;
      rk_wgt[1][2] = 0.5;

      rk_mom[0] = 0.5;
      rk_mom[1] = 0.5;
    }
    else if (rk_scheme == 3) {
      rk_wgt[0][0] = 0.0;
      rk_wgt[0][1] = 1.0;
      rk_wgt[0][2] = 1.0;
      rk_wgt[1][0] = 0.75;
      rk_wgt[1][1] = 0.25;
      rk_wgt[1][2] = 0.25;
      rk_wgt[2][0] = 1.0/3.0;
      rk_wgt[2][1] = 2.0/3.0;
      rk_wgt[2][2] = 2.0/3.0;

      rk_mom[0] = 1.0/6.0;
      rk_mom[1] = 1.0/6.0;
      rk_mom[2] = 4.0/6.0;
    }


    // EOS stuff
    mu0_ref     = getDoubleParam("MU_REF0");
    mu1_ref     = getDoubleParam("MU_REF1");
    rho0_ref = getDoubleParam("RHO_REF0");
    rho1_ref = getDoubleParam("RHO_REF1");
    p_ref   = getDoubleParam("P_REF",0.0);

    incomp = false;
    incomp = checkParam("INCOMP");
    if (mpi_rank == 0 ) cout << "incomp = " << incomp << endl;
    if (incomp) {
      if (mpi_rank == 0 ) cout << " > sos0 = " << "inf" << " " << ", sos1 = " << "inf" << endl;
    }
    
    sos0 = getDoubleParam("SOS0");
    sos1 = getDoubleParam("SOS1");
    if (mpi_rank == 0 ) cout << " > sos0 = " << sos0 << " " << ", sos1 = " << sos1 << endl;
   
    probe_first = true;

  }

  virtual ~HelmholtzVofSolver() {

    COUT1("~HelmholtzVofSolver()");

    DELETE(vof);
    DELETE(vofold);
    DELETE(u);
    DELETE(uold);
    DELETE(q_fa);
    //DELETE(q0_fa);

    DELETE(p);
    DELETE(pa);
    DELETE(p_eq);
    DELETE(rho0);
    DELETE(rho1);
    DELETE(rho0old);
    DELETE(rho1old);
    DELETE(rhoY0);
    DELETE(rho);
    DELETE(rhoY0old);
    DELETE(rhoold);
    DELETE(mu_lam);
    DELETE(mu_sgs);
    DELETE(kappa);
    DELETE(n);
    DELETE(g);

    DELETE(inv_K);

    DELETE(fgr);
    DELETE(fgr0);
    DELETE(fgr1);
    DELETE(cfl2);

    DELETE(plicPointList);
    DELETE(plic_xc);
    DELETE(plic_area);

    DELETE(sd);
    DELETE(lsd);


    DELETE(dudx);
    DELETE(sp_dpdx);
    DELETE(div);
    //DELETE(cv_flag);

    DELETE(beta);
    DELETE(iband);

    DELETE(rhs_rhoY0);
    DELETE(rhs_rhoY1);

    am = NULL;
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
        else if ( bc_type == "WM_ALG_WALL" ) {
          bcs.push_back(new AlgebraicWallModelVBc(&bfZoneVec[izone],this));
        }
        else if ( bc_type == "INLET") {
          bcs.push_back(new InletVBc(&bfZoneVec[izone],this));
        }
        else if ( bc_type == "OUTLETC") {
          bcs.push_back(new OutletVBc(&bfZoneVec[izone],this));
        }  
        else if ( bc_type == "OUTLET") {
          if (incomp) 
            bcs.push_back(new OutletVBc(&bfZoneVec[izone],this));
          else 
            bcs.push_back(new OutletVVBc(&bfZoneVec[izone],this));
        }
        else if ( bc_type == "SPONGE") {
          bcs.push_back(new SpongeVVBc(&bfZoneVec[izone],this));
        }
        else if ( bc_type == "HOOK") {
          bcs.push_back(new HookVBc(&bfZoneVec[izone],this));
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
    assert(vofold == NULL); vofold   = new double[ncv_g];
    assert(u     == NULL); u     = new double[ncv_g][3];
    assert(uold  == NULL); uold  = new double[ncv_g][3];
    assert(q_fa  == NULL); q_fa  = new double[nfa];
    //assert(q0_fa == NULL); q0_fa = new double[nfa];

    assert(p        == NULL); p        = new double[ncv_g];
    assert(pa       == NULL); pa       = new double[ncv_g];
    assert(p_eq     == NULL); p_eq     = new double[ncv_g];
    assert(rho0     == NULL); rho0     = new double[ncv_g];
    assert(rho1     == NULL); rho1     = new double[ncv_g];
    assert(rho0old  == NULL); rho0old  = new double[ncv_g];
    assert(rho1old  == NULL); rho1old  = new double[ncv_g];
    assert(rhoY0    == NULL); rhoY0    = new double[ncv_g];
    assert(rho      == NULL); rho      = new double[ncv_g];
    assert(rhoY0old == NULL); rhoY0old = new double[ncv_g];
    assert(rhoold   == NULL); rhoold   = new double[ncv_g];
    assert(mu_lam   == NULL); mu_lam   = new double[ncv_g];
    assert(mu_sgs   == NULL); mu_sgs   = new double[ncv_g];

    assert(kappa    == NULL); kappa    = new double[ncv_g];
    assert(n        == NULL); n        = new double[ncv_g][3];
    assert(g        == NULL); g        = new double[ncv_g];

    assert(lsd      == NULL); lsd      = new double[ncv_g];

    assert(inv_K    == NULL); inv_K    = new double[ncv_g];

    assert(fgr     == NULL); fgr     = new double[ncv_g];
    assert(fgr0    == NULL); fgr0    = new double[ncv_g];
    assert(fgr1    == NULL); fgr1    = new double[ncv_g];
    assert(cfl2    == NULL); cfl2    = new double[ncv];

    assert(plicPointList == NULL); plicPointList = new list<double>[ncv];
    //assert(plic_xc == NULL); plic_xc = new double[ncv_g][3];
    //assert(plic_area == NULL); plic_area = new double[ncv_g];
    assert(plic_xc == NULL); plic_xc = new double[ncv_g2][3];
    assert(plic_area == NULL); plic_area = new double[ncv_g2];
    assert(sd       == NULL); sd       = new double[ncv_g2];


    assert(dudx == NULL);        dudx = new double[ncv_g][3][3];
    assert(sp_dpdx == NULL);     sp_dpdx = new double[ncv][3];
    assert(div == NULL);         div = new double[ncv];
    assert(cv_flag == NULL);     cv_flag = new int[ncv_g];

    icv_global_g = new int8[ncv_g-ncv];

    assert(beta == NULL); beta = new double[ncv_g];
    assert(iband == NULL); iband = new int[ncv_g];

    assert(rhs_rhoY0 == NULL); rhs_rhoY0 = new double[ncv];
    assert(rhs_rhoY1 == NULL); rhs_rhoY1 = new double[ncv];

    am = NULL;


    // compute local sharpness parameter for THINC ...
    double h_avg = 0.0;
    beta_c = getDoubleParam("BETA",1.0);
    FOR_ICV_G {
      h_avg = pow((vol_cv[icv]),1.0/3.0);
      beta[icv] = beta_c/h_avg;
    }
  }

  void init() {

    //FlowSolver::init();
    StaticSolver::init(INIT_COMPACT_FACES|INIT_CV_GRAD);

    // needed by app to control contex menu
    added_data_fields.insert("interfaces");

    // needed for curvature...
    updateCvData(x_vv,REPLACE_TRANSLATE_DATA); // TODO just make x_vv updated in ghosts in StaticSolver


    sum_outlet_proj_area = 0.0;
    FOR_BCZONE {
      Param * param = getParam((*it)->getName());
      if ( param->getString(0) == "OUTLETC")
        sum_outlet_proj_area += MAG((*it)->zone_ptr->n_global);
    }

    // get dmax for signed-distance...
    d_max = 0.0;
    FOR_ICV d_max = max(d_max,r_vv[icv]);
    d_max *= 6; // 3 cell diameters



    // needed for trilinos...
    updateCvDataSeparateGhosts(icv_global,icv_global_g);

    trilinosSolverForPressure = NULL;
    setPressureSolver(pressure_solver);

    initLiquidFlux();

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

  bool b_init;

  void initComplete() {

    // if we do not have the rho/mu_lam fields, then initialize the rho/mu_lam fields from the vof field...

    COUT1("initComplete()");

    //support incompressible result.sles : only vof is available 
    //
    if (!checkDataFlag("rho0")) {
      FOR_ICV_G rho0[icv] = rho0_ref;
    }

    if (!checkDataFlag("rho1")) {
      FOR_ICV_G rho1[icv] = rho1_ref;
    }

    // if (!checkDataFlag("vof")) {
    //   CWARN(" >>>> vof variable should be defined first !!!!!! >>>>>>>> ");
    // }

    if (!checkDataFlag("rho")) {
      FOR_ICV_G rho[icv] = rho0[icv]*vof[icv] + rho1[icv]*(1.0-vof[icv]);
    }

    if (!checkDataFlag("rhoY0")) {
      FOR_ICV_G rhoY0[icv] = rho0[icv]*vof[icv]; 
    }

    if (!checkDataFlag("p")) {
      FOR_ICV_G p[icv] = 0.0; // for initial guess...
    }

    updateCvData(vof); // assume vof has been read or set
    updateCvData(u);
    updateCvData(p);
    updateCvData(rho0);
    updateCvData(rho1);
    updateCvData(rhoY0);
    updateCvData(rho);
    updateMu_lam();

    // reset advection pressure and modulus...
    FOR_ICV_G {
      inv_K[icv] = 0.0;
      pa[icv] =0.0;
    }
    updateCvData(pa);
    updateCvData(inv_K);

 
    assert( am == NULL);
    if (!checkDataFlag("q_fa")) {
      //assert(!checkDataFlag("sp_dpdx"));
      // project out divergence in u...
      FOR_ICV_G kappa[icv] = 0.0;
      FOR_IFA {
        const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
        const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
        double u_fa[3]; FOR_I3 u_fa[i] = 0.5*(u[icv0][i]+u[icv1][i]);
        q_fa[ifa] = DOT_PRODUCT(u_fa,n_fa[ifa]);
      }
      double (*u_copy)[3] = new double[ncv_g][3];
      FOR_ICV_G FOR_I3 u_copy[icv][i] = u[icv][i];
      const double dt_copy = dt;
      dt = 1.0;
      FOR_ICV_G kappa[icv] = 0.0;
      FOR_ICV_G cv_flag[icv] = 0;
      FOR_ICV FOR_I3 sp_dpdx[icv][i] = 0.0;
      FOR_ICV_G plic_area[icv] = 0.0;
      FOR_ICV_G sd[icv] = 0.0;
      dumpRange(u,ncv,"u before correction");
      if(getBoolParam("SOLVER",true)) {
        if (incomp) calcVof();
        else solveEqmState(0.0);
        b_init = false;
        solvePAndCorrectU();
        b_init = true;
      }
      dt = dt_copy;
      FOR_ICV_G FOR_I3 u[icv][i] = u_copy[icv][i];
      delete[] u_copy;
    }
    else {
      if (incomp) calcVof();
      else solveEqmState(0.0);
    }

    // calculate the interface properties...
    updateInterface();
    calcCurvature();

    //if (!checkDataFlag("n")) {
    //  assert(!checkDataFlag("g"));
    //  assert(!checkDataFlag("kappa"));
    //  updateInterface();
    //}

    StaticSolver::calcCvGrad(dudx,u);
    updateCvData(dudx);

    setFgrCvsFull();

    /// subgrid stress...
    calcSgs();

    if (incomp) {
      FOR_ICV {
        inv_K[icv] = 0.0;
        pa[icv] =0.0;
      }
    }
    updateCvData(pa);
    updateCvData(inv_K);

    // update xp...
    //computeliquid surface area : lsd
    calcLiqVolandSurf();

    //
    //buildIband();

    // report mass and momentum...
    //processLiquidFlux();

    FOR_ICV {
      rhs_rhoY0[icv] = 0.0;
      rhs_rhoY1[icv] = 0.0;
    }

    report();

  }

  void setPressureSolver (int& pressure_solver) {

    Param * param = getParam("PRESSURE_SOLVER");

    if (param != NULL) {
      string name = param->getString();
      if (name == "TRILINOS") {
        pressure_solver = TRILINOS_SOLVER;
        p_zero = getDoubleParam("P_ZERO",1.0E-8);
        p_maxiter = getIntParam("P_MAXITER",15);
        trilinos_idata[0] = getIntParam("TRILINOS_OPTION",4); //option
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
        p_zero = getDoubleParam("P_ZERO",1.0E-8);
        p_maxiter = getIntParam("P_MAXITER",2000);
        COUT1(" > PRESSURE_SOLVER = BSGSTAB TOL=" << p_zero <<
            " MAXITER=" << p_maxiter );
      }
      else if (name == "AMG") {
        pressure_solver = AMG_SOLVER;
        p_zero = getDoubleParam("P_ZERO",1.0E-8);
        p_maxiter = getIntParam("P_MAXITER",1000);
        COUT1(" > PRESSURE_SOLVER = AMG TOL=" << p_zero <<
            " MAXITER=" << p_maxiter );
      }
      else if (name == "JACOBI") {
        pressure_solver = JACOBI_SOLVER;
        p_zero = getDoubleParam("P_ZERO",1.0E-7);
        p_maxiter = getIntParam("P_MAXITER",2000);
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


  void checkDivandDt() {

    double my_div = 0.0;
    double max_div = 0.0;
    FOR_ICV {
      max_div = max(max_div, fabs(div[icv])/vol_cv[icv]);
    }
    MPI_Allreduce(&my_div, &max_div, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);

    if ( max_div > 1.0) {
      dt = dt/max_div;
      COUT1("Wanrning, dt is adjusted to satisfy divergence condition");
    }

  }




  int advanceSolution() {

    if ( (mpi_rank==0) && (step%check_interval==0) )
      cout << " > advanceSolution()" << endl;

    CtiRegister::clearCurrentData();

    // cv variables ..
    FOR_ICV_G {
      vofold[icv] = vof[icv];
      rho0old[icv] = rho0[icv];
      rho1old[icv] = rho1[icv];
      rhoY0old[icv] = rhoY0[icv];
      rhoold[icv] = rho[icv];
      FOR_I3 uold[icv][i] = u[icv][i];
    }

    // advance time step...
    //
    //checkDivandDt();
    time += dt;

    // set normal for wall film
    FOR_BCZONE (*it)->setBc();

    double max_iter = getIntParam("MAX_ITER",2);
    bool force_iter = checkParam("MAX_ITER");
    if (incomp) max_iter = 1;

    //timer.split("restore old");

    for (int iter = 0; iter < max_iter ; iter++) {

      if ( (mpi_rank==0) && (step%check_interval==0) )
        cout << " > iter: " << iter+1 << endl;

      if (iter > 0 ) {
        // restore states
        FOR_ICV_G {
          vof[icv] = (vofold[icv]);
          rhoY0[icv] =(rhoY0old[icv]);
          rho[icv] =(rhoold[icv]);;
          rho0[icv] = (rho0old[icv]);
          rho1[icv] = (rho1old[icv]);;
        }
        updateInterface();
      }

      //updateInterface();
      //timer.split("update interface1");
      double (*rhs_rhou)[3] = new double[ncv][3];
      double * A = new double[cvocv_i[ncv]];

      // reset RHS and LHS A matrix
      FOR_ICV FOR_I3 rhs_rhou[icv][i] = 0.0;
      for (int coc = 0; coc < cvocv_i[ncv]; ++coc) A[coc] = 0.0;

      solveRho(rhs_rhou, A);

      if (incomp) calcVof();
      //else solveEqmState(dt);

      updateMu_lam();
      calcSgs();
      //timer.split("calcSgs");
      solveU(rhs_rhou, A);

      delete[] rhs_rhou;
      delete[] A;
      //timer.split("solveU");

      updateInterface();
      calcCurvature();
      //timer.split("calcCurvature");
      
      // pressure correction...
      if (getBoolParam("SOLVER",true)) solvePAndCorrectU();
      //timer.split("solveP");

      // velocity gradient...
      setFgrCvsFull();
      //timer.split("Fgrset");
      StaticSolver::calcCvGrad(dudx,u);
      updateCvData(dudx);
      //timer.split("dudx");
      
      if ( (!force_iter) && checkDivandFgr(iter) ) break;
    }
   
    //calcSubIntRes();

    doProbes();
    report();

    //timer.split("report");
    return 0;
  }

  bool checkDivandFgr(const int iter) {

    const double div_ref = getDoubleParam("MAX_DIV",0.1);
    const double fgr_ref = getDoubleParam("MAX_FGR",0.01);

    double my_residue = 0.0;
    const double SMALL = 1.0E-14;
    FOR_ICV {
      my_residue = max(my_residue, fabs((p[icv] - p_eq[icv])/(p[icv]+SMALL)));
    }
    
    double residue = 0.0;
    MPI_Allreduce(&my_residue, &residue, 1, MPI_DOUBLE, MPI_MAX,mpi_comm);

    double my_div = 0.0;
    double max_div = 0.0;
    FOR_ICV {
      my_div = max(my_div, fabs(div[icv])*dt/vol_cv[icv]);
    }
    MPI_Allreduce(&my_div, &max_div, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);

    double my_fgr = 0.0;
    double max_fgr = 0.0;
    FOR_ICV {
      my_fgr = max(my_fgr, fgr[icv]);
    }

    MPI_Allreduce(&my_fgr, &max_fgr, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    if (mpi_rank == 0 && step%check_interval == 0 ) cout << " > check param:  iter, max_prel, max_div, max_fgr = " << iter+1 << " " << residue << " " << max_div << " " <<   max_fgr << endl;

    bool check = false;

    if (max_div < div_ref && max_fgr < fgr_ref) check = true;

    return(check);
  }

  void report() {

    if (step%check_interval == 0) {
     // FOR_ICV {
     //   if (rho1[icv] < 1.0E-5) cout << "rho1 " << rho1[icv] << " " <<  COUT_VEC(x_cv[icv]) << endl;
     // }

      /*
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
      */
      dumpRange(vof,ncv,"vof");
      dumpRange(u,ncv,"u");
      dumpRange(p,ncv,"p");
      dumpRange(rhoY0,ncv,"rhoY0");
      dumpRange(rho,ncv,"rho");
      dumpCfl();
      //dumpRange(kappa,ncv,"kappa");
      //dumpRange(rho0,ncv,"rho0");
      //dumpRange(rho1,ncv,"rho1");
      //dumpRange(inv_K,ncv,"inv_K");
      //dumpRange(fgr,ncv,"fgr");
      //queryVof();
      queryBcs();
    }

    
  }

  void processLiquidPMD(Param * param) {

    //list<pair<string,stringstream*> > cpFileList;

    int interval;
    string name;
    int ny, nz;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double x[3];
    double width=-1.0;
    int nr = 100;
    vector<string> varNameVec;
    
    int iarg = 0;
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "INTERVAL") {
        interval = param->getInt(iarg++);
        if (step%interval != 0)
          return;
      } else if (token == "NAME") {
        name = param->getString(iarg++);
      } else if (token == "X") {
        x[0] = param->getDouble(iarg++);
        x[1] = param->getDouble(iarg++);
        x[2] = param->getDouble(iarg++);
      } else if (token == "WIDTH") {
        width = param->getDouble(iarg++);
      } else if (token == "N") {
        nr = param->getInt(iarg++);
      } else if (token == "VARS") {
        while (iarg < param->size())
          varNameVec.push_back(param->getString(iarg++));
      }
    }
    const double radius = width*0.5;
    
    const int nvar = varNameVec.size();
    if (nvar > 1) {
      CERR("VAR shoulde be one : " <<  nvar);
    }
    assert(radius > 0.0);
    ny = 2*nr;
    nz = 2*nr;
    double xtarget = x[0];
    ymin = x[1] - radius;
    ymax = x[1] + radius;
    zmin = x[2] - radius;
    zmax = x[2] + radius;

    assert(nr>0);
    
    const double dy= (ymax-ymin)/double(ny);
    const double dz= (zmax-zmin)/double(nz);
    
    const int nbin = 4*nr*nr;
    //const int nbin = 0000;
 
    double *my_data = new double[nbin];
    double *global_data = new double[nbin];
    //double my_data[nvar][nbin];
    //double global_data[nvar][nbin];
    
    //double xtarget = (xmin+xmax)*0.5;
    double yp[ny];
    double zp[nz];
    
    for (int ybin=0; ybin< ny; ybin++) {
      yp[ybin] = ymin + 0.5*dy + ybin*dy;
    }
    for (int zbin=0; zbin < nz; zbin++) {
      zp[zbin] = zmax - 0.5*dz - zbin*dz;
    }
    
    for (int ibin=0; ibin<nbin; ibin++) {
      my_data[ibin] = 0.0;
    }
    
    // find out the correct local icv for each bin
    int8 * cvora = NULL;
    buildXora(cvora,ncv);
    assert( cvora[mpi_rank+1] - cvora[mpi_rank] == ncv );
    vector<int> cvList;
    int * my_icvp_global = new int[nbin];
    
    for (int ybin=0; ybin< ny; ybin++) {
      for (int zbin=0; zbin< nz; zbin++) {
        const int ibin = ybin + ny*zbin;
        double xbase[3] = {x[0], yp[ybin], zp[zbin]};
        my_icvp_global[ibin] = -1;
        cvAdt->buildListForPoint(cvList, xbase);
        const int icv = getClosestICVFromList(cvList,xbase);
        if (icv >= 0) my_icvp_global[ibin] = icv + cvora[mpi_rank]; // local to global cv offset
      }
    }
    
    // at this point, my_icvp_global should contain the global cv index of
    // the containing cell, or -1... get the maximum across the
    // processors to break ties or check for xbasa ethat could not be located...
    
    int * icvp_global = new int[nbin];
    MPI_Allreduce(my_icvp_global,icvp_global,nbin,MPI_INT,MPI_MAX,mpi_comm);
    
    int npnew = 0;
    
    for (int ip = 0; ip < nbin; ++ip) {
      if (icvp_global[ip] == -1) {
        if (mpi_rank == 0)
          cout << "Warning: could not locate cv for the PMD grid at " << endl;
      }
      else if (icvp_global[ip] == my_icvp_global[ip]) {
        ++npnew;
      }
    }
    
    int npnew_global =0;
    MPI_Allreduce(&npnew, &npnew_global,1,MPI_INT,MPI_SUM,mpi_comm);
    assert(npnew_global == nbin);
    
    // resize and add particles on processors that need it...
    
    if (npnew > 0) {
      for (int ybin=0; ybin<ny; ybin++) {
        for (int zbin=0; zbin<nz; zbin++) {
          const int ibin = ybin + ny*zbin;
          if ((icvp_global[ibin] != -1)&&(icvp_global[ibin] == my_icvp_global[ibin])) {
            const int ibin = ybin + ny*zbin;
            double xbase[3] = {xtarget, yp[ybin], zp[zbin]};
            const int my_icv = my_icvp_global[ibin] - cvora[mpi_rank]; 
            assert(my_icv >= 0 && my_icv < ncv);
            double sum_weights = 0.0;
            for (int coc = cvocv_i[my_icv]; coc != cvocv_i[my_icv+1]; ++coc) {
              const int icv_nbr = cvocv_v[coc];
              if (icv_nbr < ncv) {
                double d2 = DIST2(x_cv[icv_nbr], xbase);
                double this_weight = 1.0/(sqrt(d2) + 1.0E-15); // add tol to avoid 1/zero
                sum_weights += this_weight;
                CtiRegister::CtiData * data = CtiRegister::getRegisteredCtiData(varNameVec[0]);
                assert(data!= NULL);
                assert(data->getType() == DN_DATA);
                my_data[ibin] += this_weight*data->dn(icv_nbr);
              }
            }
            my_data[ibin] /= sum_weights;
          }
        }
      }
    }
    
    
    MPI_Reduce(my_data,global_data,nvar*nbin,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    
    delete[] icvp_global, my_icvp_global;
    
    if (mpi_rank == 0 ) {
      
      double pmd_y[ny];
      double pmd_z[nz];
      
      for (int ybin=0; ybin<ny; ybin++) pmd_y[ybin] = 0.0;
      for (int zbin=0; zbin<nz; zbin++) pmd_z[zbin] = 0.0;
      
      // perform integration in y and z direction....
      for (int ybin=0; ybin<ny; ybin++) {
        for (int zbin=0; zbin<nz; zbin++) {
          const int ibin = ybin + ny*zbin;
          assert(ibin>= 0); assert(ibin<nbin);
          pmd_y[ybin] += global_data[ibin]*dz;
          pmd_z[zbin] += global_data[ibin]*dy;
        }
      }
      
      ofstream ofile;
      char filename[32];
      //buildIndexedFilename(filename1,name.c_str(),step,"dat");
      sprintf(filename,"%s.%08d.%s",name.c_str(),step,"dat");
      //MiscUtils::mkdir_for_file(filename); // allow the user to have specified a subdirectory...
      
      if (MiscUtils::fileExists(filename)) {
        ofile.open(filename,ios::app);
        if (!ofile.is_open()) {
          cerr << ERRSTART << "Error: could not open WRITE PMD file: " << filename << endl;
        }
      }
      else {
        MiscUtils::mkdir_for_file(filename);
        ofile.open(filename);
        if (!ofile.is_open()) {
          cerr << ERRSTART << "Error: could not open WRITE PMDfile: " << name << ERREND << endl;
        }
        
        else { 
          ofile << setprecision(8);
          ofile << "# step " << step << " time " << time << endl;
          
          ofile << "# 1: R " << " ";
          for (int ii=0; ii< nvar; ii++) ofile << ii+2 << ": " << varNameVec[ii] << " "  ;
          ofile << endl;
          
          for (int zbin=0; zbin<nz; zbin++) {
            ofile << zp[zbin] << " " ;
            ofile << 0.5*(pmd_y[zbin]+pmd_z[zbin]) << " ";
            ofile << endl;
          }
        }
        
      ofile.close();
      
      }
    }
    delete[] my_data;
    delete[] global_data;
  
  }

  void updateLiquidFlux(list<LiquidFlux>::iterator lmf) {

    double np[3] = {1.0,0.0,0.0};
    double xp[3] = {lmf->xp[0], lmf->xp[1], lmf->xp[2]};
    int nbins = lmf->nbins;

    for (int ii = 0; ii < lmf->groupVec.size();  ii++) {
      const int icv = lmf->groupVec[ii][0];
      int ifa = lmf->groupVec[ii][1];
      int fa_sign = 1.0;
      if (ifa < 0) {
        ifa = -ifa -1;
        fa_sign = -1.0;
      }
      const int ybin = lmf->groupVec[ii][2];
      const int zbin = lmf->groupVec[ii][3];
      const int ibin = ybin + zbin*lmf->ny;
      assert(ibin>= 0 && ibin<nbins);
      const double q_flux = max(DOT_PRODUCT(u[icv], np),0.0);
      lmf->liquid_flux[ibin]    += rhoY0[icv]*q_flux*fabs(DOT_PRODUCT(n_fa[ifa],np))*dt;
      lmf->projected_area[ibin] += fabs(DOT_PRODUCT(n_fa[ifa],np))*dt;
    }

//    for (int ibin=0; ibin<nbins; ibin++) {
//      cout << "liquid flux = " << ibin << " " << lmf->liquid_flux[ibin] << " " << lmf->liquid_flux[ibin]/lmf->projected_area[ibin] << endl;
//    }

  }


  void constructLiquidFluxGroup(list<LiquidFlux>::iterator lmf) {

    //cout << "what a change ? " << lmf->ymin << " " << lmf->ymax << " " << lmf->zmin << " " << lmf->xp[0] << endl;
    double np[3] = {1.0,0.0,0.0};
    double xp[3] = {lmf->xp[0], lmf->xp[1], lmf->xp[2]};
    int nbins = lmf->nbins;

    // construct band ..
    //int * cv_flag     = new int[ncv_g];
    for (int icv = 0; icv < ncv; ++icv) {
      const double dn =
        (x_cv[icv][0]-xp[0])*np[0] +
        (x_cv[icv][1]-xp[1])*np[1] +
        (x_cv[icv][2]-xp[2])*np[2];
      if (dn >= 0.0) {
        cv_flag[icv] = 1;
      }
      else {
        cv_flag[icv] = -1;
      }
    }
    
    updateCvData(cv_flag);

    for (int icv = 0; icv < ncv; ++icv) {
      for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        if (cv_flag[icv]*cv_flag[icv_nbr] < 0) {
        // an odd product means this cv is on the boundary...
          cv_flag[icv] *= 2; // makes -1 -> -2
          assert((cv_flag[icv] == 2)||(cv_flag[icv] == -2));
          break;
        }
      }
    }
    
    updateCvData(cv_flag);
    
    for (int icv = 0; icv < ncv; ++icv) {
      if ((cv_flag[icv] == 1)||(cv_flag[icv] == -1)) {
        for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        if ((cv_flag[icv_nbr] == 2)||(cv_flag[icv_nbr] == -2)) {
          cv_flag[icv] *= 4;
          assert((cv_flag[icv] == 4)||(cv_flag[icv] == -4));
          break;
        }
        }
      }
    }
    
    updateCvData(cv_flag);

    double my_area = 0.0;

    int index = 0;
    for (int icv = 0; icv < ncv; ++icv) {
      if (abs(cv_flag[icv]) > 1) {
        for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
          const int ifa = faocv_v[foc];
          const int icv0 = cvofa[ifa][0];
          const int icv1 = cvofa[ifa][1];
          assert((icv0 == icv)||(icv1 == icv));
          if (cv_flag[icv0]*cv_flag[icv1] < 0.0) {
            if (cv_flag[icv] < 0.0) {
              const int ybin = (int)floor( (x_cv[icv][1]-lmf->ymin)/(lmf->ymax - lmf->ymin)*double(lmf->ny) );
              const int zbin = (int)floor( (x_cv[icv][2]-lmf->zmin)/(lmf->zmax - lmf->zmin)*double(lmf->nz) );
              if (( ybin >=0 && ybin < lmf->ny ) && (zbin >=0 && zbin < lmf->nz) ) {
                //set icv, ifa group
                if (icv == icv0) {
                  assert(DOT_PRODUCT(n_fa[ifa],np)>= -1.0E-6);
                  //if (DOT_PRODUCT(n_fa[ifa],np) < 0.0) cout << "what " << DOT_PRODUCT(n_fa[ifa],np) << " " << COUT_VEC(n_fa[ifa]) << endl;
                  int* group;
                  group = new int[4];
                  group[0] = icv;
                  group[1] = ifa;
                  group[2] = ybin;
                  group[3] = zbin;
                  lmf->groupVec.push_back(group);
                  double fa_sign = 1.0;
                  my_area += fabs(DOT_PRODUCT(n_fa[ifa],np));
               //   cout << "index ifa area = " << index << " " << ifa << " " << fa_sign*DOT_PRODUCT(n_fa[ifa],np) << endl;
                }
                else {
                  assert(icv == icv1);
                  assert(DOT_PRODUCT(n_fa[ifa],np)<= 1.0E-6);
                  double fa_sign = -1.0;
                  int* group;
                  group = new int[4];
                  group[0] = icv;
                  group[1] = -ifa-1;
                  group[2] = ybin;
                  group[3] = zbin;
                  lmf->groupVec.push_back(group);
                  my_area += fabs(DOT_PRODUCT(n_fa[ifa],np));
               //   cout << "index ifa area = " << index << " " << ifa << " " << fa_sign*DOT_PRODUCT(n_fa[ifa],np) << endl;
                }
                index++;
              }
            }
          }
        }
      }
    }
   
    double total_flux_area = 0.0;
    MPI_Reduce(&my_area,&total_flux_area, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
    if (mpi_rank == 0 ) cout << "sum flux area = " << total_flux_area << endl;

    //lmf->projected_area = new double[nbins];
    //lmf->liquid_flux = new double[nbins];

    //reset projected area;
    //for (int ibin =0; ibin<nbins; ibin++) {
    //  lmf->projected_area[ibin] = 0.0;
    //  lmf->liquid_flux[ibin] = 0.0;
    // }

 /*
    for (int ii = 0; ii < lmf->groupVec.size();  ii++) {
      int icv = lmf->groupVec[ii][0];
      int ifa = lmf->groupVec[ii][1];
      int fasign = 1.0;
      if (ifa < 0) {
        ifa = -ifa -1;
        fasign = -1.0;
      }
      int ybin = lmf->groupVec[ii][2];
      int zbin = lmf->groupVec[ii][3];
      int ibin = ybin + zbin*lmf->ny;
      lmf->projected_area[ibin] += fasign*DOT_PRODUCT(n_fa[ifa],np);
    }

    
    double *global_projected_area;
    global_projected_area = new double[nbins];
    MPI_Reduce(lmf->projected_area, global_projected_area, nbins, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);

    if (mpi_rank == 0 ) {
      double sum_area = 0.0;
      for (int ibin=0; ibin<nbins; ibin++){
        sum_area += global_projected_area[ibin];
        cout << "projected area = " << ibin << " " << global_projected_area[ibin] << endl;
      }
      cout << "total flux area = " << sum_area << endl;
    }

    delete[] global_projected_area;

*/
  }
   
  void queryVof() {
    FOR_PARAM_MATCHING("QUERY_VOF") {
      processQueryVof(&(*param));
    }
  }

  void updateLiquidPMD() {
    FOR_PARAM_MATCHING("ECN.PMD") {
      processLiquidPMD(&(*param));
    }
  }

  void doProbes() {
    calcLiqVolandSurf();
    processLiquidFlux();
    queryVof();
    updateLiquidPMD();
    //custom probe here...
  }

  void calcLiqVolandSurf() {

    // compute surface area density (minimum surface area density from Chesnel et al 2011. 
    // Large Eddy Simulation of Liquid Jet Atomization, Atomization Sprays, vol. 21, no. 9//
    //FOR_ICV_G lsd[icv]= 2.4*sqrt(vof[icv]*(1.0-vof[icv]))/pow(vol_cv[icv],1.0/3.0);
    FOR_ICV_G lsd[icv]= 2.31*sqrt(vof[icv]*(1.0-vof[icv]))/(5.0*r_vv[icv]);

    updateCvData(lsd);

    // if (step%check_interval == 0 ) {
    // test volume and surface ....
    //   double my_volume = 0.0;
    //   double my_surface = 0.0;
    
    //   FOR_ICV {
    //     my_volume += vof[icv]*vol_cv[icv];
    //     my_surface += lsd[icv]*vol_cv[icv];
    //   }
    
    //   double volume, surface;
    //   MPI_Reduce(&my_volume, &volume, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
    //   MPI_Reduce(&my_surface, &surface ,1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
    
    //   if (mpi_rank == 0 ) {
    //     cout << " > total phase0 volume = " << volume << " " << " total phase0 surface - "  << " " << surface << " SMD = " << 6*volume/surface << endl;
    //   }
    
    //  }  
    
  }

  void processQueryVof(Param * param) {

    //if (mpi_rank == 0 && step%check_interval == 0 ) cout << " >> Query VoF starts" << endl;

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
        else if (geom == "ANNULAR_CONE_X") {
          //---------------------------------------------------------------
          // HCP_WINDOW ANNULAR_TCONE <x0> <y0> <z0> <x1> <y1> <z1> <r0> <r1>
          // ---------------------------------------------------------------
          geom_type = ANNULAR_CONE_X;
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
        cout << "[QUERY_VOF:" << name<< "]# 1:step 2:time 3:total volume 4:phase0 volume 5:vof_avg 6:mixing_variance 7:plic surface 8: vof surface  9:SMD_plic 10:SMD_vof " << endl;
      probe_first = false;
    }

    // compute plic surface
    buildPlic();

    // sum the plic area in a given area ...
    double my_stats_sum[4] = {0.0,0.0,0.0,0.0};
    const double dp_tol = getDoubleParam("DP_TOL",1.0E-6);

    FOR_ICV {
      if (geom_type == BOX) {
        if ( x_cv[icv][0] >= geom_data[0] && x_cv[icv][0] <= geom_data[1] ) {
          if ( x_cv[icv][1] >= geom_data[2] && x_cv[icv][1] <= geom_data[3] ) {
            if ( x_cv[icv][2] >= geom_data[4] && x_cv[icv][2] <= geom_data[5] ) {
              my_stats_sum[0] += vol_cv[icv];
              my_stats_sum[1] += vof[icv]*vol_cv[icv];
              my_stats_sum[2] += plic_area[icv];
              my_stats_sum[3] += 2.31*sqrt(vof[icv]*(1.0-vof[icv]))/(5.0*r_vv[icv])*vol_cv[icv];
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
              my_stats_sum[2] += plic_area[icv];
              my_stats_sum[3] += 2.31*sqrt(vof[icv]*(1.0-vof[icv]))/(5.0*r_vv[icv])*vol_cv[icv];
          }
        }
      }
    }

    double stats_sum[4];
    MPI_Allreduce(my_stats_sum,stats_sum,4,MPI_DOUBLE,MPI_SUM,mpi_comm);

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
        stats_sum[1] << " " << vof_avg << " " << var_sum/stats_sum[0] << " " << stats_sum[2] << " " << stats_sum[3] << " " << 
        6.0*stats_sum[1]/stats_sum[2]  << " " << 6.0*stats_sum[1]/stats_sum[3] <<  endl;
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



  void updateInterface() {

    if ( step%check_interval==0 )
      COUT1(" > updateInterface()");

    //cleanupVof();
    flagInterfaceCvs();
    //timer.split("flaginterface");

    calcNormal();
    //timer.split("calcnormal");
    //build Plic Surface and compute Area in plic_area[] ...
    //buildPlic();
    // normal vector is required to reconstruct hyperbolic tangent...
    calcGfromVof();
    //timer.split("calcGfromVof");

  }

  void calcCurvature() {

    //calcCurvatureSignedDistance();
    calcCurvatureSimple();
    //calcCurvatureDirect();
    
    //calcCurvatureFromG();
    //smoothCurvature();
    //calcCurvatureDirect();
  }

  void buildSignedDistance() {
    
    //if (step%check_interval == 0) COUT1(" > buildSignedDistance()...");

    assert(n_ef); // need extended faces for this approach

    buildPlicPolys();

    // get smallest distance to interface...
    // TODO: investigate weighted combo instead of min. Also look into filtering.
    // init signed distance to max allowable distance...
   
    FOR_ICV_G2 sd[icv] = d_max;
    FOR_ICV {
      if (cv_flag[icv] >= 1) {
        //sd[icv] = computeDistToPlicPoly(x_cv[icv],icv);
        sd[icv] = fabs(g[icv]); // same as |g| in interface cells, will sign again later
      }
    }

    FOR_IEF {
      const int icv0 = cvoef[ief][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvoef[ief][1]; assert((icv1 >= 0)&&(icv1 < ncv_g2));
      if ((cv_flag[icv0] >= 1)&&(cv_flag[icv1] <= 0)) {
        sd[icv1] = min(sd[icv1],computeDistToPlicPoly(x_cv[icv1],icv0));
      }
      else if ((icv1 < ncv)&&(cv_flag[icv1] >= 1)&&(cv_flag[icv0] <= 0)) {
        sd[icv0] = min(sd[icv0],computeDistToPlicPoly(x_cv[icv0],icv1));
      }
    }

    //updateCv2Data(vof);
    updateCv2Data(sd,MIN_DATA);

    // sign the distance based on vof
    FOR_ICV_G {
      if (vof[icv] < 0.5) {
        sd[icv] = -sd[icv];
      }
    }

  }


  void calcCurvatureSignedDistance() {

    // build signed distance from plic and smooth...
    buildSignedDistance();

    //timer.split("build signed distance");

    // compute d(sd)dx /||d(sd)dx||...
    
    double (*n_sd)[3] = new double[ncv_g][3];
    StaticSolver::calcCvGrad(n_sd,sd);
    
    updateCvData(n_sd,REPLACE_ROTATE_DATA); 
    FOR_ICV_G {
      const double mag = MAG(n_sd[icv]);
      if (mag > 1.0E-14) 
        FOR_I3 n_sd[icv][i] /= mag;
    }

    // calc kappa = div.n...
    FOR_ICV kappa[icv] = 0.0;
    FOR_IFA {
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      double n_sd_fa[3]; 
      FOR_I3 n_sd_fa[i] = 0.5*(n_sd[icv0][i]+n_sd[icv1][i]);
      //FOR_I3 n_sd_fa[i] = (abs(sd[icv1])*n_sd[icv0][i]+abs(sd[icv0])*n_sd[icv1][i])/(abs(sd[icv1])+abs(sd[icv0]));
      const double mag = MAG(n_sd_fa);
      if (mag > 1.0E-14) {
        const double flux = DOT_PRODUCT(n_sd_fa,n_fa[ifa])/mag;
        kappa[icv0] -= flux;
        kappa[icv1] += flux;
      }
    }
    delete[] n_sd;
    FOR_ICV kappa[icv] *= inv_vol[icv];

    // bound by cell length scale...
    FOR_ICV {
      if (kappa[icv] > 0.0) {
        kappa[icv] = min(kappa[icv],1.0/r_vv[icv]);
      }
      else if (kappa[icv] < 0.0) {
        kappa[icv] = max(kappa[icv],-1.0/r_vv[icv]);
      }
    }
    updateCvData(kappa);

    //dumpRange(kappa,ncv,"kappa from sd");

  }

  double computeDistToPlicPoly(const double xp[3],const int icv) {

    double dist = d_max;
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
      dist = min(dist,sqrt(getPointToTriDist2(xp,x0,x1,plic_xc[icv])));
    }
    return dist;
  }



  void calcCurvatureFromG() {


    //buildSignedDistance();

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

  void smoothCurvature() {

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


  void buildPlic() {


    //if (step%check_interval == 0 ) COUT1(" > buildPlic() ..");

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

    if (step%check_interval==0) {

      double my_minX = 1e20;
      double my_maxX = -1e20;
      double my_minY = 1e20;
      double my_maxY = -1e20;
      double my_minZ = 1e20;
      double my_maxZ = -1e20;

      FOR_ICV {
        list<double>::iterator it = plicPointList[icv].begin();
        if (plicPointList[icv].size() > 0) {
          while (it != plicPointList[icv].end()) {
              my_minX = min(my_minX, *it); 
              my_maxX = max(my_maxX, *it); ++it;
              my_minY = min(my_minY, *it); 
              my_maxY = max(my_maxX, *it); ++it;
              my_minZ = min(my_minZ, *it); 
              my_maxZ = max(my_maxZ, *it); ++it;
          }
        }
      }

      double minX = 1e20;
      double maxX = -1e20;
      double minY = 1e20;
      double maxY = -1e20;
      double minZ = 1e20;
      double maxZ = -1e20;
      MPI_Allreduce(&my_minX, &minX, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&my_maxX, &maxX, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&my_minY, &minY, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&my_maxY, &maxY, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&my_minZ, &minZ, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&my_maxZ, &maxZ, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      if (mpi_rank==0) cout << "[CHECK] minX = " << minX << ", maxX = " << maxX << endl;
      if (mpi_rank==0) cout << "[CHECK] minY = " << minY << ", maxY = " << maxY << endl;
      if (mpi_rank==0) cout << "[CHECK] minZ = " << minZ << ", maxZ = " << maxZ << endl;


      int my_npts = 0;
      int my_max_npts = 0;

      for (int icv = 0; icv < ncv; ++icv) {
        int size = plicPointList[icv].size();
        my_npts += size;
        my_max_npts = max(my_max_npts, size);
      }
      int npts = 0;
      int max_npts = 0;
      MPI_Allreduce(&my_npts, &npts, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&my_max_npts, &max_npts, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if (mpi_rank==0) cout << "[CHECK] npts = " << npts/3 << ", max_npts = " << max_npts/3 << endl;


      FOR_ICV {
        if (cv_flag[icv] >= 1 && x_cv[icv][0] > 0.2 && x_cv[icv][0] < 0.8) {
          //cout << "icv=" << icv << " n=(" << n[icv][0] << "," << n[icv][1] << "," << n[icv][2] << ") g=" << g[icv] << endl;
        }

      }
    }

    // now calculate centroid...
    FOR_ICV {
      if (cv_flag[icv] >= 1) {
        if (step%check_interval == 0 && icv == 5880) {
          cout << "icv=" << icv << " npts=" << int(plicPointList[icv].size()/3) << endl;
          list<double>::iterator it = plicPointList[icv].begin();
          while (it != plicPointList[icv].end()) {
            //cout << "plicPoints.x=" << *(++it) << endl;
            ++it;
            cout << "plicPoints.y=" << *(++it) << endl;
            ++it;
            //cout << "plicPoints.z=" << *(++it) << endl;
          }
          
        }
        plic_area[icv] = computePointListAreaAndCentroid(plic_xc[icv],plicPointList[icv]);
      }
      else {
        FOR_I3 plic_xc[icv][i] = 0.0;
        plic_area[icv] = 0.0;
      }
    }
    if (step%check_interval==0) {
      double my_area = 0.0;
      for (int icv = 0; icv < ncv; ++icv) {
        my_area += plic_area[icv];
      }
      double area = 0.0;
      MPI_Allreduce(&my_area, &area, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if (mpi_rank==0) cout << "[CHECK] area = " << area << endl;

      dumpRange(plic_area,ncv,"plic_area");

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

    //calcNormalLeastSquares();
    calcNormalSimple();
    // update normal near the wall
    FOR_BCZONE (*it)->setBc();

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

    //if (step%check_interval == 0) dumpRange(vof_s,ncv,"vof_s");
    StaticSolver::calcCvGrad(n,vof_s);

    FOR_ICV {
      const double mag = MAG(n[icv]);
      if (mag > 0.0) {
        const double inv_mag = 1.0/mag;
        FOR_I3 n[icv][i] *= -inv_mag;
      }
    }
    updateCvData(n);

    //if (step%check_interval == 0) dumpRange(n,ncv,"simple normal");
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
 //   FOR_ICV {
 //     if (vof[icv] < vof_zero) vof[icv] = 0.0;
 //     if (vof[icv] > 1.0 - vof_zero) vof[icv] = 1.0;
 //   }
  }


  void limitVof() {
    // clean up vof residue ...
    FOR_ICV {
      vof[icv] = max(0.0,min(vof[icv],1.0));
    }
  }


  void updateMu_lam(){
    FOR_ICV_G mu_lam[icv] = mu1_ref*(1.0-vof[icv]) + mu0_ref*vof[icv];
  }

  void solveRho( double (*rhs_rhou)[3], double* A) {

    // compute vof Rhs without/with flux limiter ...
    if (step%check_interval == 0 ) COUT1(" > solve rho");

    //double *rhs_rhoY0 = new double[ncv];
    //double *rhs_rhoY1 = new double[ncv];

    //reset Rhs
    FOR_ICV {
      rhs_rhoY0[icv]  = 0.0;
      rhs_rhoY1[icv]  = 0.0;
    }

    //add boundary flux 
    FOR_BCZONE {
      (*it)->addMassFlux(rhs_rhoY0, rhs_rhoY1);
    }

    calcInternalFlux(rhs_rhoY0, rhs_rhoY1, rhs_rhou, A);

    FOR_ICV {
      rhoY0[icv]  = rhoY0old[icv]  + rhs_rhoY0[icv]*dt*inv_vol[icv];
      rho[icv]  = rhoold[icv]  + (rhs_rhoY0[icv]+rhs_rhoY1[icv])*dt*inv_vol[icv];
    }

    // vof is not used any more....
    updateCvData(rhoY0);
    updateCvData(rho);

    //delete[] rhs_rhoY0;
    ///delete[] rhs_rhoY1;

  }



  void solveU(double (*rhs_rhou)[3], double* A) {

    // compute vof Rhs without/with flux limiter ...
    if (step%check_interval == 0 ) COUT1(" > solve u");

    //add boundary flux 
    FOR_BCZONE {
      (*it)->addMomentumFlux(A,rhs_rhou);
    }

    FOR_ICV {
      FOR_I3 rhs_rhou[icv][i] += vol_cv[icv]*rhoold[icv]*uold[icv][i]/dt;
    }

    calcViscRhs(rhs_rhou);
    //updateMu_lam();
    buildLhs(A);

    // predict u...
    if (getBoolParam("SOLVER",true)) solveCvJacobi(u,A,rhs_rhou,u_zero,u_relax,u_maxiter,false);

    FOR_ICV FOR_I3 u[icv][i] += dt*gravity[i];
    updateCvData(u);

    if (step%check_interval == 0 ) dumpRange(u,ncv,"u before correction");

  }


  void calcDivergence() {

    FOR_ICV div[icv] = 0.0;
    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; 
      const int icv1 = cvofa[ifa][1]; 
      //const double q_mid = 0.5*(q_fa[ifa]+q0_fa[ifa]);
      const double q_mid = q_fa[ifa];
      div[icv0] += q_mid;
      if (icv1 < ncv)
        div[icv1] -= q_mid;
    }

    FOR_BCZONE (*it)->addFlux(div);

  }

  void checkDivergence(double _dt) {

    FOR_ICV {
      if (vol_cv[icv] < fabs(_dt*div[icv])) 
       CWARN("divergence is too big, please reduce dt");
      //assert( vol_cv[icv] >= fabs(_dt*div[icv]));
    }

  }

  void checkPressure() {

    FOR_ICV {
      const double residue = (p[icv]*inv_K[icv] - pa[icv] + dt*div[icv]*inv_vol[icv])*vol_cv[icv];
      //assert(fabs(residue) < 1.0E-5);
    }

  }

  void calcVof() {

    FOR_ICV {
      vof[icv] = rhoY0[icv]/rho0[icv];
      vof[icv] = max(0.0,min(1.0,vof[icv]));
    }

    updateCvData(vof);

  }

  void solveEqmState(double _dt) {

    // compute equilibrium pressure from advection 
    const double sos0sq = sos0*sos0;
    const double sos1sq = sos1*sos1;

    calcDivergence();
    //checkDivergence(_dt);

    FOR_ICV {
      const double m0 = rhoY0[icv];
      const double m1 = rho[icv] - rhoY0[icv];

      if (m0 > mass_zero && m1 > mass_zero) {

        // solve 2nd order polynomial for p
        const double a = 1.0;
        const double b = rho0_ref*sos0sq + rho1_ref*sos1sq - (m0*sos0sq + m1*sos1sq);
        const double c = (rho0_ref*rho1_ref - (m0*rho1_ref + m1*rho0_ref))*sos0sq*sos1sq;

        assert( (b*b-4.0*c) >= 0.0);
        
        double R1 = (-b + sqrt(b*b-4.0*c))/2.0;
        double R2 = (-b - sqrt(b*b-4.0*c))/2.0;
        //cout << "R1 R2 = " << R1 << " " << R2 << " " << b << " " << c << endl;
       
        p_eq[icv] = R1;
        //if(fabs(R1) > fabs(R2)) {
         // cout << "negative p = " << R1 << " " << R2 << " "<< endl;
          //p_eq[icv] = R2;
        //}
        //assert(fabs(R1) < fabs(R2));
        //cout << "p_eq = ? " << R1 << " " << R2 << endl;
        // check residue ...
        double residue = a*p_eq[icv]*p_eq[icv] + b*p_eq[icv] + c;
       // if (fabs(residue)/rho1_ref/sos1sq > 1.0E-5) {
       //   cout << "Not converged peq = " << p_eq[icv] << " " << residue << endl;
       // }
        
        // compute equilibrium vof ..
        const double vof0 = max(m0,0.0)*sos0*sos0/(p_eq[icv] + rho0_ref*sos0*sos0);
        const double vof1 = max(m1,0.0)*sos1*sos1/(p_eq[icv] + rho1_ref*sos1*sos1);
        // check volume compatibility condition
        assert(p_eq[icv] + rho0_ref*sos0sq > 0.0);
        assert(p_eq[icv] + rho1_ref*sos1sq > 0.0);
        assert(vof0 > -vof_zero);
        assert(vof1 > -vof_zero);
        //assert(fabs(vof0+vof1-1.0) < 1.0E-11);
      
        //hack...
        //vof[icv] = vof0/(vof0+vof1);
        vof[icv] = vof0;

        // compute density from EOS for p_eq = c^2(rho_k-rho_ref)
        rho0[icv] = p_eq[icv]/sos0sq + rho0_ref;
        rho1[icv] = p_eq[icv]/sos1sq + rho1_ref;
        assert(rho0[icv] > 0.0);
        assert(rho1[icv] > 0.0);
        //rho0[icv] = rho*(p_eq[icv]/sos0sq + rho0_ref) / (vof[icv]*(p_eq[icv]/sos0sq + rho0_ref) + (1.0-vof[icv])*(p_eq[icv]/sos1sq + rho1_ref));
        //rho1[icv] = rho*(p_eq[icv]/sos1sq + rho1_ref) / (vof[icv]*(p_eq[icv]/sos0sq + rho0_ref) + (1.0-vof[icv])*(p_eq[icv]/sos1sq + rho1_ref));
      }
      else {
        if (m0 > m1) {
          // single water phase 
          p_eq[icv] = sos0sq*(m0 - rho0_ref);
          vof[icv] = 1.0;
          rho0[icv] = m0;
          rho1[icv] = max(p_eq[icv]/sos1sq + rho1_ref,mass_zero);
        }
        else {
          // single gas phase
          p_eq[icv] = sos1sq*(m1 - rho1_ref);
          vof[icv] = 0.0;
          rho0[icv] =max(p_eq[icv]/sos0sq + rho0_ref,mass_zero);
          rho1[icv] =m1;
        }
      } // single phase ends

      //compute advection pressure ...
      const double p0 = p_eq[icv] + rho0[icv]*sos0sq*_dt*inv_vol[icv]*div[icv];
      const double p1 = p_eq[icv] + rho1[icv]*sos1sq*_dt*inv_vol[icv]*div[icv];
      
      //compute mixture compressibility (inverse of  bulk modulus)
      inv_K[icv] = vof[icv]/rho0[icv]/sos0sq + (1.0-vof[icv])/rho1[icv]/sos1sq;
      assert(inv_K[icv] > 0.0);
      
      // compute mixture advection pressure p^* weighted by bulk modulus
      pa[icv] = (vof[icv]*p0/rho0[icv]/sos0sq + (1.0-vof[icv])*p1/rho1[icv]/sos1sq );
      
      // check pressure compatibility condition
      //const double residue2 = ( p_eq[icv]*inv_K[icv] - pa[icv]  + _dt*div[icv]*inv_vol[icv])*vol_cv[icv];

      // initial condition for the next pressure ...
      p[icv]  = p_eq[icv];
      
      //cout << "icv div pa residue = " << icv << " " << div[icv] << " " <<  pa[icv]/inv_K[icv] << " " << residue2 <<endl;
      //  if (fabs(residue2) > 1.0E-10) 
      //    cout << "residue2 = " << icv << " " << residue2 << endl;
    }
      
   
    updateCvData(rho0);
    updateCvData(rho1);
    updateCvData(pa);
    updateCvData(p_eq);
    updateCvData(p);
    updateCvData(inv_K);
    updateCvData(vof);

    //setPhaseStateBc(); 

  }

  
  void setPhaseStateBc() {


    // this is not necessary, but just in case...
    //
    // calc internal flux ...
    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      
      if ( rhoY0[icv0] <= mass_zero && rhoY0[icv1] > mass_zero) {
        rho0[icv0] = rho0[icv1];
        //sos0[icv0] = sos0[icv1];
      }
      else if (vof[icv0] > vof_zero && vof[icv1] <= vof_zero) {
        rho0[icv1] = rho0[icv0];
        //sos0[icv1] = sos0[icv0];
      }
      
      if ( 1.0-vof[icv0] <= vof_zero && 1.0-vof[icv1] > vof_zero){
        rho1[icv0] = rho1[icv1];
        //sos1[icv0] = sos1[icv1];
      }
      else if (1.0-vof[icv0] > vof_zero && 1.0-vof[icv1] <= vof_zero) {
        rho1[icv1] = rho1[icv0];
        //sos1[icv1] = sos1[icv0];
      }
    }
    
  }


  void calcInternalFlux(double* rhs_rhoY0, double* rhs_rhoY1,  double (*rhs_rhou)[3], double* A) {

    // compute flux limiters to guarantee the positiveness 
    double *rho0_out = NULL;
    double *rho1_out = NULL;

    rho0_out = new double[ncv_g];
    rho1_out = new double[ncv_g];

    //calcVofLimiters(vof_out, vol_out);
    calcMassLimiters(rho0_out, rho1_out);

    //if (mpi_rank == 1 ) cout << " out " << rho0_out[61790] << " " << rhoY0[61790] << " " << vof[61790] << " " << 1.0-vof[61790] <<  endl; 

    const double fgr_max = getDoubleParam("FGR_MAX",0.25);
    const double eps = 1.0E-9;

    // calc internal flux ...
    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));

      const int coc00 = cvocv_i[icv0];
      int coc01 = coc00+1;
      while (cvocv_v[coc01] != icv1) {
        ++coc01;
        assert(coc01 != cvocv_i[icv0+1] );
      }
      
      //const double q_mid = 0.5*( q_fa[ifa]+ q0_fa[ifa]);
      const double q_mid = q_fa[ifa];

      double rhoY0_flux = 0.0;
      double rhoY1_flux = 0.0;
      double rho_flux = 0.0;
      double rhou_flux[3] = {0.0,0.0,0.0};

      double vof_avg = 0.5*(vof[icv0]+vof[icv1]);
      double vof_f = 0.0;

      const double rhoY1_0 = rho[icv0] - rhoY0[icv0];
      const double rhoY1_1 = rho[icv1] - rhoY0[icv1];

      const double rhoY0_avg = 0.5*(rhoY0[icv0] + rhoY0[icv1]);
      const double rhoY1_avg = 0.5*(rhoY1_0 + rhoY1_1);

      double rhoY0_f = 0.0;
      double rhoY1_f = 0.0;

      double fgr_avg = max(fgr[icv0], fgr[icv1]);
      fgr_avg = min(fgr_avg,fgr_max);
      double fgr0_max = max(fgr0[icv1],fgr0[icv0]);
      double fgr1_max = max(fgr1[icv1],fgr1[icv0]);
      double gamma; // blending coeffficient gamma = 0: full upwind, gamma =1 : central

      if ( (vof[icv0] > vof_zero && 1.0-vof[icv0] > vof_zero) || (vof[icv1] > vof_zero && 1.0-vof[icv1] > vof_zero)) {
/*
        double dx0[3], dx1[3],dx[3];
        FOR_I3 dx0[i] = x_fa[ifa][i] - x_cv[icv0][i];
        FOR_I3 dx1[i] = x_fa[ifa][i] - x_cv[icv1][i];
        FOR_I3 dx[i] = dx0[i] - dx1[i];

        double normal0 = -DOT_PRODUCT(dx0,n[icv0]);
        double normal1 = -DOT_PRODUCT(dx1,n[icv1]);


        const double vof0 =  0.5*(1.0+tanh(beta[icv0]*(normal0+g[icv0])));
        const double vof1 =  0.5*(1.0+tanh(beta[icv1]*(normal1+g[icv1])));

        const double vof_mid = 0.5*(vof0+vof1);
        const double norm = pow(0.25,1.0/10.0);
        const double factor = pow(vof_mid*(1.0-vof_mid),1.0/10.0)/norm;
        //computing blending coefficients..
        double dvofdx0 = 0.5*(1.0-pow(tanh(beta[icv0]*(normal0+g[icv0])),2))*beta[icv0];
        double dvofdx1 = 0.5*(1.0-pow(tanh(beta[icv1]*(normal1+g[icv1])),2))*beta[icv1];
        double gradvof = 0.5*fabs(dvofdx0+dvofdx1);
        gamma = ((gradvof > 1.0E-10) ? fabs(vof[icv1]-vof[icv0])/MAG(dx)/gradvof : 0.0);
        gamma = max(0.0, min(gamma,1.0));
        gamma = pow(1.0 - gamma,4);

        if (vof[icv0] < 0.001  && vof[icv1] < 0.001) {
            gamma = 0.0; // fully upwind ...
        }

        if (q_mid >= 0.0) {
          vof_f = (1.0-gamma)*vof0+gamma*vof_avg;
          rhoY0_f = (1.0-gamma)*rho0[icv0]*vof0   +  gamma*rhoY0_avg;
          rhoY1_f = (1.0-gamma)*rho1[icv0]*(1.0-vof0) + gamma*rhoY1_avg;
        //  if (rho0_out[icv0] > 0.0 ) rhoY0_f = min(rhoY0_f, rhoY0[icv0]/rho0_out[icv0]*rhoY0_f);
        //  if (rho1_out[icv0] > 0.0 ) rhoY1_f = min(rhoY1_f, rhoY1[icv0]/rho1_out[icv0]*rhoY1_f);
        }
        else {
          vof_f = (1.0-gamma)*vof1+gamma*vof_avg;
          rhoY0_f = (1.0-gamma)*rho0[icv1]*vof1   +  gamma*rhoY0_avg;
          rhoY1_f = (1.0-gamma)*rho1[icv1]*(1.0-vof1) + gamma*rhoY1_avg;
       //   if (rho0_out[icv1] > 0.0 ) rhoY0_f = min(rhoY0_f, rhoY0[icv1]/rho0_out[icv1]*rhoY0_f);
       //   if (rho1_out[icv1] > 0.0 ) rhoY1_f = min(rhoY1_f, rhoY1[icv1]/rho1_out[icv1]*rhoY1_f);
        }
       */
       if (q_mid >= 0.0)
       {
         rhoY0_f = rho0[icv0] * vof[icv0];
         rhoY1_f = rho1[icv0] * (1.0 - vof[icv0]);
       }
       else
       {
         rhoY0_f = rho0[icv1] * vof[icv1];
         rhoY1_f = rho1[icv1] * (1.0 - vof[icv1]);
       }
       
      }
      //single-phase cell
      else {
        vof_f = vof_avg;
        rhoY0_f = rhoY0_avg;
        rhoY1_f = rhoY1_avg;

    //    if (q_mid >= 0.0) {
    //      if (rho0_out[icv0] > 0.0 ) rhoY0_f = min(rhoY0_f, rhoY0[icv0]/rho0_out[icv0]*rhoY0_f);
    //      if (rho1_out[icv0] > 0.0 ) rhoY1_f = min(rhoY1_f, rhoY1[icv0]/rho1_out[icv0]*rhoY1_f);
    //    }
    //    else {
    //      if (rho0_out[icv1] > 0.0 ) rhoY0_f = min(rhoY0_f, rhoY0[icv1]/rho0_out[icv1]*rhoY0_f);
    //      if (rho1_out[icv1] > 0.0 ) rhoY1_f = min(rhoY1_f, rhoY1[icv1]/rho1_out[icv1]*rhoY1_f);
    //    }
       // flux for non-vof cells 
      }

      if (q_mid >= 0.0) {
        if (rho0_out[icv0] > rhoY0[icv0]) fgr0_max = 0.5;
        if (rho1_out[icv0] > rhoY1_0) fgr1_max = 0.5;
      }
      else {
        if (rho0_out[icv1] > rhoY0[icv1]) fgr0_max = 0.5;
        if (rho1_out[icv1] > rhoY1_1) fgr1_max = 0.5;
      }
      
      //rhoY0_flux = rhoY0_f*q_mid - max(fgr0_max,fgr_avg)*fabs(q_mid)*(rhoY0[icv1]-rhoY0[icv0]);
      rhoY0_flux = rhoY0_f*q_mid;
      //rhoY1_flux = rhoY1_f*q_mid - max(fgr1_max,fgr_avg)*fabs(q_mid)*(rhoY1_1-rhoY1_0);
      rhoY1_flux = rhoY1_f*q_mid;
      rho_flux = rhoY0_flux + rhoY1_flux;
      //FOR_I3 rhou_flux[i]  =0.5*( 0.5*rho_flux*(u[icv0][i]+u[icv1][i]) - fgr_avg*fabs(rho_flux)*(u[icv1][i]-u[icv0][i]) );  
      FOR_I3 rhou_flux[i]  =0.5*( 0.5*rho_flux*(u[icv0][i]+u[icv1][i]) );  
      
      rhs_rhoY0[icv0]           -= rhoY0_flux;
      rhs_rhoY1[icv0]           -= rhoY1_flux;
      FOR_I3 rhs_rhou[icv0][i] -= rhou_flux[i];
      if (icv1 < ncv) {
        rhs_rhoY0[icv1]           += rhoY0_flux;
        rhs_rhoY1[icv1]           += rhoY1_flux;
        FOR_I3 rhs_rhou[icv1][i] += rhou_flux[i];
      }

      const double mf_coeff1 = 0.25*rho_flux;
      const double mf_coeff2 = 0.5*fgr_avg*fabs(rho_flux);

      //build Lhs matrix
      //A[coc00] += mf_coeff1 + mf_coeff2;
      A[coc00] += mf_coeff1;
      //A[coc01] += mf_coeff1 - mf_coeff2;
      A[coc01] += mf_coeff1;

      if (icv1 < ncv) {

        const int coc11 = cvocv_i[icv1];
        int coc10 = coc11+1;
        while (cvocv_v[coc10] != icv0) {
          ++coc10;
          assert(coc10 != cvocv_i[icv1+1] );
        }

        //A[coc11] -= mf_coeff1 - mf_coeff2;
        A[coc11] -= mf_coeff1;
        //A[coc10] -= mf_coeff1 + mf_coeff2;
        A[coc10] -= mf_coeff1;

      }

     
    }

    delete[] rho0_out;
    delete[] rho1_out;

  }

  /*
  void calcVofLimiters(double* vof_out, double* vol_out) {

    // reset
    FOR_ICV_G {
      vof_out[icv] = 0.0;
      vol_out[icv] = 0.0;
    }

    // boundary flux...
    FOR_BCZONE {
      (*it)->addVofFlux(vof_out);
    }
    FOR_ICV_G {
      vof_out[icv]  = fabs(min(vof_out[icv],0.0));
    }

    for (vector<VofBc*>::const_iterator it = bcs.begin(); it != bcs.end(); ++it) {
      if ((*it)->q_bf != NULL) {
        for (int ibf = 0; ibf < (*it)->zone_ptr->nbf; ++ibf) {
          const int icv0 = (*it)->zone_ptr->cvobf[ibf];
          vol_out[icv0] += fabs(min((*it)->q_bf[ibf],0.0)); // 0.5*dt/vol below...
        }
      }
    }

    // calc internal flux ...
    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      const double q_mid = q_fa[ifa];

      double vof_flux = 0.0;

      double vof_f = 0.0;
      double vof_avg = 0.5*(vof[icv0]+vof[icv1]);

      if ( (vof[icv0] > vof_zero && 1.0-vof[icv0] > vof_zero) || (vof[icv1] > vof_zero && 1.0-vof[icv1] > vof_zero)) {

        double dx0[3], dx1[3],dx[3];
        FOR_I3 dx0[i] = x_fa[ifa][i] - x_cv[icv0][i];
        FOR_I3 dx1[i] = x_fa[ifa][i] - x_cv[icv1][i];
        FOR_I3 dx[i] = dx0[i] - dx1[i];

        double normal0 = -DOT_PRODUCT(dx0,n[icv0]);
        double normal1 = -DOT_PRODUCT(dx1,n[icv1]);

        const double vof0 =  0.5*(1.0+tanh(beta[icv0]*(normal0+g[icv0])));
        const double vof1 =  0.5*(1.0+tanh(beta[icv1]*(normal1+g[icv1])));

        //computing blending coefficients..
        double dvofdx0 = 0.5*(1.0-pow(tanh(beta[icv0]*(normal0+g[icv0])),2))*beta[icv0];
        double dvofdx1 = 0.5*(1.0-pow(tanh(beta[icv1]*(normal1+g[icv1])),2))*beta[icv1];
        double gradvof = 0.5*fabs(dvofdx0+dvofdx1);
        double gamma = ((gradvof > 1.0E-10) ? fabs(vof[icv1]-vof[icv0])/MAG(dx)/gradvof : 0.0);
        gamma = max(0.0, min(gamma,1.0));
        gamma = pow(1.0 - gamma,4);

        const bool upwind = checkParam("UPWIND");
        if (upwind) gamma = 0.0;

        if (vof[icv0] < 0.01  && vof[icv1] < 0.01) {
            gamma = 0.0; // fully upwind ...
        }
        
        if (q_mid >= 0.0) {
          const double vof_f0 = (1.0-gamma)*vof0+gamma*vof_avg;
          vof_f = vof_f0;
        }
        else {
          const double vof_f0 = (1.0-gamma)*vof1+gamma*vof_avg;
          vof_f = vof_f0;
        }
        // interface cell 
        vof_flux = vof_f*q_mid;
      }
      else {
        double vof_f0 = vof_avg; // just for symmetry...
        vof_f = vof_f0;
        vof_flux = vof_f*q_mid;
      }

      vof_out[icv0] += max(vof_flux,0.0);
      vol_out[icv0] += max(q_mid,0.0);
      if (icv1 < ncv) {
        vof_out[icv1] += fabs(min(vof_flux, 0.0));
        vol_out[icv1] += fabs(min(q_mid,0.0));
      }
    }

    
    FOR_ICV {
      vof_out[icv] *= dt/vol_cv[icv];
      vol_out[icv] *= dt/vol_cv[icv];
    }
    // compute out CFL number ...is it necessary? ....
    //calcCfl(vol_out,dt);
    updateCvData(vof_out);
    updateCvData(vol_out);

  }
*/

  void calcMassLimiters(double* rho0_out, double* rho1_out) {

    // reset
    FOR_ICV_G {
      rho0_out[icv] = 0.0;
      rho1_out[icv] = 0.0;
    }

    // boundary flux...
    FOR_BCZONE {
      (*it)->addMassFlux(rho0_out, rho1_out);
    }
    FOR_ICV_G {
      rho0_out[icv] = fabs(min(rho0_out[icv],0.0));
      rho1_out[icv] = fabs(min(rho1_out[icv],0.0));
    }

    // calc internal flux ...
    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      //const double q_mid = 1.25*q_fa[ifa]-0.25*q0_fa[ifa];
      //const double q_mid = 0.5*(q_fa[ifa]+q0_fa[ifa]);
      const double q_mid = q_fa[ifa];


      double vof_flux = 0.0;
      double rhoY0_flux = 0.0;
      double rhoY1_flux = 0.0;

      double vof_f = 0.0;
      double vof_avg = 0.5*(vof[icv0]+vof[icv1]);

      const double v0_avg  = 0.5*(1.0/rho0[icv0] + 1.0/rho0[icv1]);
      const double v1_avg  = 0.5*(1.0/rho1[icv0] + 1.0/rho1[icv1]);
      const double rho0_avg = 1.0/v0_avg;
      const double rho1_avg = 1.0/v1_avg;

      const double rhoY1_0 = rho[icv0] - rhoY0[icv0];
      const double rhoY1_1 = rho[icv1] - rhoY0[icv1];
      const double Y0_icv0  = rhoY0[icv0]/rho[icv0];
      const double Y0_icv1  = rhoY0[icv1]/rho[icv1];
      const double Y1_icv0  = 1.0 - Y0_icv0;
      const double Y1_icv1  = 1.0 - Y0_icv1;

      const double v_avg   = 0.5*(1.0/rho[icv0] + 1.0/rho[icv1]);
      double Y0_avg  = 0.5*(Y0_icv0 + Y0_icv1);
      double Y1_avg  = 0.5*(Y1_icv0 + Y1_icv1);

      const double rhoY0_avg = 0.5*(rhoY0[icv0] + rhoY0[icv1]);
      const double rhoY1_avg = 0.5*(rhoY1_0 + rhoY1_1);

      double rhoY0_f = 0.0;
      double rhoY1_f = 0.0;

      const double fgr0_max = max(fgr0[icv0],fgr0[icv1]);
      const double fgr1_max = max(fgr1[icv0],fgr1[icv1]);

      if ( (vof[icv0] > vof_zero && 1.0-vof[icv0] > vof_zero) || (vof[icv1] > vof_zero && 1.0-vof[icv1] > vof_zero)) {

        double dx0[3], dx1[3],dx[3];
        FOR_I3 dx0[i] = x_fa[ifa][i] - x_cv[icv0][i];
        FOR_I3 dx1[i] = x_fa[ifa][i] - x_cv[icv1][i];
        FOR_I3 dx[i] = dx0[i] - dx1[i];

        double normal0 = -DOT_PRODUCT(dx0,n[icv0]);
        double normal1 = -DOT_PRODUCT(dx1,n[icv1]);


        const double vof0 =  0.5*(1.0+tanh(beta[icv0]*(normal0+g[icv0])));
        const double vof1 =  0.5*(1.0+tanh(beta[icv1]*(normal1+g[icv1])));

        //computing blending coefficients..
        double dvofdx0 = 0.5*(1.0-pow(tanh(beta[icv0]*(normal0+g[icv0])),2))*beta[icv0];
        double dvofdx1 = 0.5*(1.0-pow(tanh(beta[icv1]*(normal1+g[icv1])),2))*beta[icv1];
        double gradvof = 0.5*fabs(dvofdx0+dvofdx1);
        double gamma = ((gradvof > 1.0E-10) ? fabs(vof[icv1]-vof[icv0])/MAG(dx)/gradvof : 0.0);
        gamma = max(0.0, min(gamma,1.0));
        gamma = pow(1.0 - gamma,4);

        if (vof[icv0] < 0.001  && vof[icv1] < 0.001) {
            gamma = 0.0; // fully upwind ...
        }
        
        if (q_mid >= 0.0) {
          const double vof_f0 = (1.0-gamma)*vof0+gamma*vof_avg;
          vof_f = vof_f0;
          rhoY0_f = (1.0-gamma)*rho0[icv0]*vof0   +  gamma*rhoY0_avg;
          rhoY1_f = (1.0-gamma)*rho1[icv0]*(1.0-vof0) + gamma*rhoY1_avg;
        }
        else {
          const double vof_f0 = (1.0-gamma)*vof1+gamma*vof_avg;
          vof_f = vof_f0;
          rhoY0_f = (1.0-gamma)*rho0[icv1]*vof1   +  gamma*rhoY0_avg;
          rhoY1_f = (1.0-gamma)*rho1[icv1]*(1.0-vof1) + gamma*rhoY1_avg;
        }
      }
      else {
        rhoY0_f = rhoY0_avg;
        rhoY1_f = rhoY1_avg;
      }

      rhoY0_flux = rhoY0_f*q_mid - fgr0_max*fabs(q_mid)*(rhoY0[icv1]-rhoY0[icv0]);
      rhoY1_flux = rhoY1_f*q_mid - fgr1_max*fabs(q_mid)*(rhoY1_1-rhoY1_0);

      if (rhoY0_f > 0.0 && q_mid > 0.0 ) rho0_out[icv0] += max(rhoY0_flux,0.0);
      if (rhoY1_f > 0.0 && q_mid > 0.0 ) rho1_out[icv0] += max(rhoY1_flux,0.0);
      if (icv1 < ncv) {
        if (rhoY0_f > 0.0 && q_mid < 0.0 ) rho0_out[icv1] += fabs(min(rhoY0_flux,0.0));
        if (rhoY1_f > 0.0 && q_mid < 0.0 ) rho1_out[icv1] += fabs(min(rhoY1_flux,0.0));
      }
    }

    
    FOR_ICV {
      rho0_out[icv] *= dt/vol_cv[icv];
      rho1_out[icv] *= dt/vol_cv[icv];
    }
    updateCvData(rho0_out);
    updateCvData(rho1_out);

  }

  void setFgrCvsFull() {

    // parameter that controls how large of a nondimensional gibbs remainder
    // is considered large.  this value is O(10^-2).  if you decrease the parameter
    // the schemes become less oscillatory as it aims to disallow entropy violations
    // the analysis would suggest that the full discrete gibbs remainder should be
    // used to set the dissipation, but it appears that only the pressure contributions
    // are being used.

    // for Gibbs remainder for the mixture flow.
    // no need to compute the mixture entropy because the gibb free energy changes compensate entropy changes due to mass fraction changes.


    const double fgr_max = getDoubleParam("FGR_MAX",0.5);
    const double gr_limit = getDoubleParam("FGR_LIMIT", 0.001);

    double * vv2 = new double[ncv_g];
    const double eps = 1.0E-16;

    for (int icv = 0; icv < ncv_g; ++icv){
      vv2[icv]    = 0.0;
      fgr0[icv]   = 0.0;
      fgr1[icv]   = 0.0;
      fgr[icv]    = 0.0;
    }
   
   
    FOR_IFA {

      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];

      const double rhoY1_0 = rho[icv0] - rhoY0[icv0];
      const double rhoY1_1 = rho[icv1] - rhoY0[icv1];
      const double rho_avg  = 0.5*(rho[icv0] + rho[icv1]);

      const double p_avg = 0.5*(p[icv1] + p[icv0]);
      const double vof_avg = 0.5*(vof[icv1] + vof[icv0]);

      const double dp = p[icv1] - p[icv0];
      const double drho0 = rho0[icv1] - rho0[icv0];
      const double drho1 = rho1[icv1] - rho1[icv0];

      // isentropic check ...
      const double rho0_avg = 0.5*(rho0[icv0] + rho0[icv1]);
      const double rho1_avg = 0.5*(rho1[icv0] + rho1[icv1]);
      const double  dgr0 = fabs( (drho0 - dp/sos0/sos0)/max(rho0_avg,mass_zero));
      const double  dgr1 = fabs( (drho1 - dp/sos1/sos1)/max(rho1_avg,mass_zero));
      const double factor = 1.0; //#100%
      const double gr = max(factor*(vof_avg*dgr0 + (1.0-vof_avg)*dgr1)-gr_limit, 0.0);

      // icv0, icv1 are both valid in this loop ..
      vv2[icv0] = max(vv2[icv0],fabs(gr));
      vv2[icv1] = max(vv2[icv1],fabs(gr));


      if (rhoY0[icv0] < 0.0 || rhoY0[icv1] < 0.0) {
        fgr0[icv0] = 0.5;
        fgr0[icv1] = 0.5;
      }

      if (rhoY1_0 < 0.0 || rhoY1_1 < 0.0) {
        fgr1[icv0] = 0.5;
        fgr1[icv1] = 0.5;
      }
      
    }
    updateCvData(vv2);
    updateCvData(fgr0);
    updateCvData(fgr1);

//    if (step%check_interval == 0) {
//      dumpRange(vv2,ncv_g,"vv2");
//    }

    FOR_ICV { 
      //fgr[icv] = gr_func(vv2[icv],gr_nrm,gr_width,fgr_max);
      fgr[icv] = min(vv2[icv], fgr_max);
    }

    updateCvData(fgr); 
    
    if ( checkParam("NO_DISS") ) {
      FOR_ICV_G {
        fgr[icv] = 0.0;
      } 
    } else if ( checkParam("ALL_DISS")) {
      FOR_ICV_G { 
        fgr[icv] = fgr_max;
      }
    }
    
    delete[] vv2;
    
  }

  inline double gr_func(const double &x, const double &nrm, const double &l_width, const double &fgr_max) {
    
    //const double nrm4 = nrm*nrm*nrm*nrm;
    //return erf(x*x*x*x/nrm4);
    
    //const double l_width = 5.0e-03; // logistic function width param
    return fgr_max/(1.0+exp(-2.0*(x-nrm)/l_width));
    
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
      const double qmid = q_fa[ifa];
      cfl[icv0] += max(qmid,0.0);
      if (icv1 < ncv)
        cfl[icv1] += fabs(min(qmid,0.0));
    }

    FOR_ICV cfl[icv] *= 2.0*dt_/vol_cv[icv];

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
    FOR_ICV cfl2[icv] = cfl[icv];
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

    if (step%check_interval == 0 ) COUT1(" > solvePAndCorrectU()...");

    // add virtual divergence to u^*...
    if (b_init) {
      FOR_ICV_G {
        u[icv][0] = (2*M_PI*x_cv[icv][0]/1.0)*0.01;
        u[icv][1] = 0.0;
        u[icv][2] = 0.0;
      }
    } 

    updateCvData(u);

    // reset div for n+1
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
      //const double sp_tension_A_fa = sp_vol_fa*sigma*kappa_fa*(vof[icv1]-vof[icv0])*area_over_delta_fa[ifa];
      const double sp_tension_A_fa = sp_vol_fa*sigma*kappa_fa*(pow(vof[icv1],1)-pow(vof[icv0],1))*area_over_delta_fa[ifa];
      //double sp_tension_A_fa;
      //const double vof_norm = fabs(vof[icv1]-vof[icv0]);
      //if ( vof_norm != 0.0 )  sp_tension_A_fa = sp_vol_fa*sigma*kappa_fa*(vof[icv1]-vof[icv0])/vof_norm*area_over_delta_fa[ifa];
      //else sp_tension_A_fa = sp_vol_fa*sigma*kappa_fa*(vof[icv1]-vof[icv0])*area_over_delta_fa[ifa];

      double u_fa[3]; FOR_I3 u_fa[i] = 0.5*(u[icv0][i]+u[icv1][i]);
      //const double un_A = DOT_PRODUCT(u_fa,n_fa[ifa]) + dt*sp_tension_A_fa;
      const double un_A = DOT_PRODUCT(u_fa,n_fa[ifa]);
      div[icv0] += un_A;
      if (icv1 < ncv)
        div[icv1] -= un_A;
    }

    if (!incomp) {
      FOR_BCZONE (*it)->updateBc();
      FOR_BCZONE (*it)->addFlux(div);
    } 
    else {
      FOR_BCZONE {
        Param * param = getParam((*it)->getName());
        if ( param->getString(0) != "OUTLET" && param->getString(0) != "OUTLETC")
          (*it)->addFlux(div);
      }
      FOR_BCZONE (*it)->updateBc();
      FOR_BCZONE {
        Param * param = getParam((*it)->getName());
        if ( param->getString(0) == "OUTLET" || param->getString(0) == "OUTLETC")
          (*it)->addFlux(div);
      }
    }

    if (step%check_interval==0){
      dumpRange(div,ncv,"div before correction");

      double my_sum = 0.0;
      FOR_ICV  my_sum += div[icv]*div[icv];
      double sum;
      MPI_Reduce(&my_sum,&sum,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      sum = sqrt(sum)/ncv_global;
      if (mpi_rank == 0)
        cout << " > div L2 norm before correction: " << sum << endl;
    }

    if (!incomp) {
      FOR_ICV {
        div[icv] -= vol_cv[icv]*pa[icv]/dt;
      }
    }

    FOR_ICV div[icv] /= dt;

    //timer.split("calc poisson rhs");

    //if (step%check_interval == 0 ) dumpRange(div,ncv,"div before correction");
    // build the lhs matrix...

    double * A = new double[cvocv_i[ncv]];
    //buildCvDivSpVolGrad(A);
    buildHelmholtzOp(A);

    //timer.split("build poisson lhs");

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

            //if ( (niters < 0 || niters == p_maxiter) && !rebuildPoissPrec) { // Aztec failed.. rebuild the prec and try again ..
            if (niters < 0 || niters == p_maxiter) { // removed !rebuildPoissPrec because sometimes you need to try again even if you rebuilt
              if (mpi_rank == 0 && step%check_interval == 0 ) cout << " > " << " " << "update ML Preconditioner "<< endl;
              trilinosSolverForPressure->updateMLPreconditioner() ;
              niters = trilinosSolverForPressure->solve_again(p,div, step%check_interval==0);
            }

          }
          //if ( niters > p_maxiter/2 )
          //  rebuildPoissPrec = true ;
          //else
          //  rebuildPoissPrec = false ; // always rebuild..

          //if ( checkParam("FORCE_ML_REBUILD"))
          // always rebuild Poisson
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
          solveCvCg2(p,A,div,p_zero,p_maxiter,true); // need to replace eventually...
          break;
        }
      case JACOBI_SOLVER:
        {
          solveCvJacobi2(p,A,div,p_zero,p_relax,p_maxiter,false);
          break;
        }
      case AMG_SOLVER:
        {
          if (am == NULL) {
            const double agglomeration_factor = getDoubleParam("AGGLOMERATION_FACTOR",8.0); // 2^3 = 8
            const bool split_orphaned_colors = getBoolParam("SPLIT_ORPHANED_COLORS",true);
            const double amg_coeff = getDoubleParam("AMG_COEFF",pow(agglomeration_factor,-1.0/3.0)); // 1/2
            const int ncg = getIntParam("NCG",5);
            am = new Multigrid::AlgebraicMultigrid();
            am->init(vol_cv,x_cv,A,cvocv_i,cvocv_v,icv_global,rbi_g,NULL,ncv,ncv_g,ncv_g,ncv_global,ncg,
                agglomeration_factor,split_orphaned_colors,amg_coeff);
          }
          else {
            am->update(A);
          }
          const double relax = getDoubleParam("RELAX",0.7); // jacobi needs something around this
          const int nsmooth = getIntParam("NSMOOTH",5);
          const int dnsmooth = getIntParam("DNSMOOTH",0);
          const string smoother = getStringParam("SMOOTHER","GMRES");
          const int maxiter_coarsest = getIntParam("MAXITER_COARSEST",15);
          const string solver_coarsest = getStringParam("SOLVER_COARSEST","GMRES");
          const bool verbose = (step%check_interval == 0);
          const int maxcycle = getIntParam("MAXCYCLE",p_maxiter);
          const string cycle = getStringParam("CYCLE","V_CYCLE");
          const bool b_linear_prolongation = checkParam("LINEAR_PROLONGATION");
          am->solve(p,div,p_zero,cycle,maxcycle,b_linear_prolongation,nsmooth,dnsmooth,smoother,relax,maxiter_coarsest,solver_coarsest,verbose);
        }
    }
        
    delete[] A;
    updateCvData(p); // FH: should be done in solver, so not necessary?

    // correct u...
    calcSpDpdxandFlux();

    FOR_ICV FOR_I3 u[icv][i] -= dt*sp_dpdx[icv][i];
    updateCvData(u);
    
    //timer.split("correct u");

    {
      if (step%check_interval == 0 ) {
        FOR_ICV div[icv] = 0.0;
        FOR_IFA {
          const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
          const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
          div[icv0] += q_fa[ifa];
          if (icv1 < ncv)
            div[icv1] -= q_fa[ifa];
        }
        FOR_BCZONE (*it)->addFlux(div);
        //checkPressure();
        dumpRange(div,ncv,"div after correction");
  
        double my_sum = 0.0;
        FOR_ICV  my_sum += div[icv]*div[icv];
        double sum;
        MPI_Reduce(&my_sum,&sum,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
        sum = sqrt(sum)/ncv_global;
        if (mpi_rank == 0)
          cout << " > div L2 norm after correction: " << sum << endl;
      }
    }
    //timer.split("calc div");
  }
 
  void buildHelmholtzOp(double * A) {

    for (int coc = 0; coc < cvocv_i[ncv]; ++coc)
      A[coc] = 0.0;

    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      //const double sp_vol_fa = 1.0/rho_f;
      //const double coeff = sp_vol_fa*area_over_delta_fa[ifa];
      const double coeff = SP_VOL_FA(rho[icv0], rho[icv1])*area_over_delta_fa[ifa];

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

    if (!incomp) {
      // diagonal term for Pressure`
      FOR_ICV {
        const int coc = cvocv_i[icv];
        A[coc] -= vol_cv[icv]/dt/dt*inv_K[icv];
      }// diagonal operator 
    }
    
  }

  int solveCvJacobi2(double * phi,const double *A,const double *rhs,
      const double zero,const double relax,const int maxiter,const bool verbose) {

    double (*res)      = new double[ncv];

    int iter = 0;
    int done = 0;
    while (done == 0) {
      
      iter++;
      
      // calculate the residual...
      FOR_ICV {
        const int coc_f = cvocv_i[icv];
        res[icv] = rhs[icv] - A[coc_f]*phi[icv];
        for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          res[icv] -= A[coc]*phi[icv_nbr];
        }
        res[icv] /= A[coc_f];
      }
      
      // update the active u's...
      FOR_ICV phi[icv]   += relax*res[icv];
      
      // and ghosts...
      updateCvData(phi);
      
      // check...
      double my_res_max = 0.0;
      FOR_ICV {
        my_res_max += res[icv]*res[icv];
      }
      double res_max=0.0;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
      res_max = sqrt(res_max)/ncv_global;
      if (mpi_rank == 0) {
        if ((verbose)||(iter > maxiter/2))
          cout << " > solveCvJacobi iter, res_max: " << iter << " " << res_max <<  endl;
        if (res_max < zero) {
          done = 1;
        }
        else if (iter > maxiter) {
          cout << " > Warning: solveCvJacobi did not converge after " << maxiter << " iters, res_max: " << res_max << " " << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
      
    }
    
    delete[] res;
    
    return( done == 1 );
    
  }

  int solveCvCg2(double * phi,const double * const A,const double * const rhs,const double zero,const int maxiter,const bool verbose) {

    // assume we come in with a consistent initial condition...
    
    // we need the following work arrays...
    
    double * res      = new double[ncv];
    double * v        = new double[ncv];
    double * p        = new double[ncv_g];
    double * inv_diag = new double[ncv];
    
    // initialize...
    for (int icv = 0; icv < ncv; ++icv)
      inv_diag[icv] = 1.0/A[cvocv_i[icv]];
    
    for (int icv = 0; icv < ncv; ++icv)
      p[icv] = 0.0;
    double rho = 1.0;
    
    // calculate the residual in rhs format...
    for (int icv = 0; icv < ncv; ++icv) {
      res[icv] = rhs[icv] - A[cvocv_i[icv]]*phi[icv];
      for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        res[icv] -= A[coc]*phi[icv_nbr];
      }
    }
    
    // diagonal precon/compute normalized residual...
    for (int icv = 0; icv < ncv; ++icv)
      v[icv] = res[icv]*inv_diag[icv];
    
    int iter = 0;
    int done = 0;
    while (done == 0) {
      
      ++iter;
      
      double rho_prev = rho;
      if (fabs(rho_prev) < 1.0E-20)
        rho_prev = -1.0E-20; // -1.0E-20? seems to help
      
      double my_rho = 0.0;
      for (int icv = 0; icv < ncv; ++icv)
        my_rho += res[icv]*v[icv];
      MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

      if (rho == 0.0) {
        if ( (mpi_rank == 0) && verbose && step%check_interval == 0) 
          cout << " > iter, rho = " << iter << "    " << rho << ", escaping..." << endl;
        break;
      }
      
      double beta = rho/rho_prev;
      for (int icv = 0; icv < ncv; ++icv)
        p[icv] = v[icv] + beta*p[icv];
      updateCvData(p);
      
      // v = [Ap]{p}...
      for (int icv = 0; icv < ncv; ++icv) {
        v[icv] = A[cvocv_i[icv]]*p[icv];
        for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          v[icv] += A[coc]*p[icv_nbr];
        }
      }
      
      double my_gamma = 0.0;
      for (int icv = 0; icv < ncv; ++icv)
        my_gamma += p[icv]*v[icv];
      double gamma;
      MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
      if (fabs(gamma) < 1.0E-20)
        gamma = 1.0E-20;
      
      const double alpha = rho/gamma;
      
      // check if we are done...
      if (iter%3 == 0) {
        
        for (int icv = 0; icv < ncv; ++icv)
          phi[icv] += alpha*p[icv];
        updateCvData(phi);
        
        // recompute the residual...
        for (int icv = 0; icv < ncv; ++icv) {
          res[icv] = rhs[icv] - A[cvocv_i[icv]]*phi[icv];
          for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
            const int icv_nbr = cvocv_v[coc];
            res[icv] -= A[coc]*phi[icv_nbr];
          }
        }
        
        for (int icv = 0; icv < ncv; ++icv)
          v[icv] = res[icv]*inv_diag[icv];
        
        // compute the max (L-infinity) normalized residual...
        double  my_res_max = 0.0;
        for (int icv = 0; icv < ncv; ++icv)
          my_res_max += v[icv]*v[icv];
        double res_max=0.0;
        MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
        res_max = sqrt(res_max)/ncv_global;
        if (mpi_rank == 0) {
          // only share the last half of the convergence behaviour...
          if ((verbose || (iter > maxiter/2)) && step%check_interval == 0)
            cout << " > iter, res_max = " << iter << "    " << res_max << endl;
          if (res_max <= zero) {
            //cout << "-> Successfully converged error to " << res_max << endl;
            done = 1;
          }
          else if (iter > maxiter) {
            cout << "Warning: solveCvCg did not converge after " << maxiter <<
              " iters, res_max: " << res_max << endl;
            done = 2;
          }
        }
        MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
        
      }
      else {
        
        // update full phi including ghosts...
        for (int icv = 0; icv < ncv_g; ++icv)
          phi[icv] += alpha*p[icv];
        
        for (int icv = 0; icv < ncv; ++icv) {
          // on the other iterations, use this approximation to update
          // the unreduced residual...
          res[icv] -= alpha*v[icv];
          // still need to compute v, diag precon for next iteration...
          v[icv] = res[icv]*inv_diag[icv];
        }
        
      }
      
    }
    
    delete[] res;
    delete[] v;
    delete[] p;
    delete[] inv_diag;
    
    // let the calling routine know if we were successful...
    return( done == 1 );
    
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
      /*
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
      */
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
      const double coeff = SP_VOL_FA(rho[icv0],rho[icv1])*area_over_delta_fa[ifa];

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

  // calc inviscid specific forces (pressure and surface tension forces per mass) - units of accelaration
  void calcSpDpdxandFlux() {

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

      //const double sp_force_inv_A_fa = -sp_vol_fa*(p[icv1]-p[icv0]-sigma*kappa_fa*(pow(vof[icv1],1)-pow(vof[icv0],1)))*area_over_delta_fa[ifa];
      const double sp_force_inv_A_fa = -sp_vol_fa*(p[icv1]-p[icv0])*area_over_delta_fa[ifa];
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
      //rho[icv] = rho_g*(1.0-vof[icv]) + rho_l*vof[icv];
      mu_lam[icv] = mu1_ref*(1.0-vof[icv]) + mu0_ref*vof[icv];
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
      cout << " > Building interface band..." << endl;

    double eps = vof_zero;

    const int isolated = 81;
    const int nband = 15;

    FOR_ICV_G iband[icv] = 0;

    FOR_ICV {
      if (vof[icv] >= 0.5)
        iband[icv] = 1; // well-resolved liquid interfaced
      else if (vof[icv] < vof_zero)
        iband[icv] = 0; // pus gas ....
      else
        iband[icv] = isolated; // less resolved liquid 
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

    double beta_init = getDoubleParam("BETA_INIT",beta_c);
    FOR_ICV {
      double beta0 = beta_init / pow(vol_cv[icv],1.0/3.0);
      vof[icv] = 0.5*(1.0+tanh(beta0*(g[icv])));
      vof[icv] = min(1.0,max(0.0,vof[icv]));
    }

    updateCvData(vof);
    limitVof();

  }

  virtual void loadBalanceHook(StripedMesh* sm) {
    COUT1("HelmholtzVofSolver::loadBalanceHook()");
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

  void initLiquidFlux() {

    FOR_PARAM_MATCHING("LIQUID.FLUX"){
      // add a new liquid flux ...
      liquidFluxList.push_back(LiquidFlux());
      LiquidFlux *lmf = &(liquidFluxList.back());
      lmf->init(&(*param));
    }

    for (list<LiquidFlux>::iterator lmf = liquidFluxList.begin(); lmf != liquidFluxList.end(); ++lmf) {
      constructLiquidFluxGroup(lmf);
      //lmf->report();
    }

  }


  void processLiquidFlux() {

    // if ((mpi_rank == 0) && (step%check_interval == 0))
    //  cout << "processLiquidFlux()" << endl;

    for (list<LiquidFlux>::iterator lmf = liquidFluxList.begin(); lmf != liquidFluxList.end(); ++lmf) {
      if (step%lmf->getSampleInterval() == 0) {
        updateLiquidFlux(lmf);
      }
      if (step%lmf->getWriteInterval() == 0) {
        lmf->write(time,step);
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
