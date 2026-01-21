#ifndef STATESOLVER_HPP
#define STATESOLVER_HPP

#include "StaticSolver.hpp"
#include "CtiMemory.hpp"
#include "MpiTimer.hpp"

// for conjugate heat transfer
#include "Cht.hpp"
#include "FluentReader.hpp"

// for Lp
#include "LpStuff.hpp"
#include "LpInjector.hpp"
#include "LpBcs.hpp"

// and new Lp...
//#include "Lsp.hpp"

using namespace MiscUtils;

#define scalar_transport_kind 0
#define scalar_transport_src_kind 1
#define scalar_non_transport_kind 2

class ScalarField {
public:

  double * phi;
  string name;
  int idx;

  // transported (0), not transported (1),
  // could be other kinds in the future ..

  int kind;

  // although the compressible solvers assume unity Lewis number,
  // the Helmholtz solvers can have scalar specific diffusion

  double Sc_lam;
  double Sc_t;

  ScalarField(const string& _name) {

    name    = _name;
    phi     = NULL;

    kind    = 0; // default is zero (transported)
    idx     = -1;

    Sc_lam = 0.0;
    Sc_t   = 0.0;

  }

  ScalarField(const string& _name,const int _kind) {

    name    = _name;
    phi     = NULL;

    kind    =  _kind;
    idx     = -1;

    Sc_lam = 0.0;
    Sc_t   = 0.0;

  }

  ScalarField(const string& _name,const int _kind,const double _Sc_lam,const double _Sc_t) {

    name    = _name;
    phi     = NULL;

    kind    =  _kind;
    idx     = -1;

    Sc_lam = _Sc_lam;
    Sc_t   = _Sc_t;

  }

  ~ScalarField() {

    DELETE(phi);

  }

};

//==============================================
// the face hash holds the two cells and
// indices into the compressed extended
// and compact faces.  this presumes that
// there was sufficient repetition in the
// faces such that this indirection is beneficial
// at the the level of the amt of memory loaded
//===============================================

class FaHash {
public:
  int cvofa[2];
  int ief_compressed;
  int icf_compressed;
};

class FaNC {
public:
  double n[3];
  double c[3];
};

class FaCompact {
public:
  double unit_n[3];
  double area_over_delta;
  double area;
};

//====================================================
// for cache ordering, here are some auxilliary
// sorting comparators to aid in the final
// extended/compact face ordering..
//====================================================

class cmpFaHashIcf {
public:
  bool operator()(const FaHash& a, const FaHash& b) const {
    return a.icf_compressed < b.icf_compressed;
  }
};

class cmpFaHashIef {
public:
  bool operator()(const FaHash& a, const FaHash& b) const {
    return a.ief_compressed < b.ief_compressed;
  }
};

// sort based on the icv1 candidate.. this will
// place all of the processor boundary faces
class cmpFaHashCv1 {
public:
  bool operator()(const FaHash& a, const FaHash& b) const {
    return a.cvofa[1] < b.cvofa[1];
  }
};

class cmpMinCv {
public:
  bool operator()(const FaHash& a, const FaHash& b) const {
    int min_icva = min(a.cvofa[0],a.cvofa[1]);
    int min_icvb = min(b.cvofa[0],b.cvofa[1]);
    return min_icva < min_icvb;
  }
};

class ForceIO {
public:

  string name;
  int interval;
  vector<string> bf_zones;
  vector<double> moment_centroids;
  int write_rank;

  ForceIO() {
    name       = "";
    interval   = -1;
    write_rank = 0;
  }

};

// simple container for pmm/mrf stuff...
enum BodeForceType {
  PMM_TYPE, // porous media model
  MRF_TYPE, // moving reference frame
  HOOK_TYPE // user defined
};
class BodyForce {
public:
  string name;        // name
  int type;           // PMM_TYPE/MRF_TYPE
  int interval;       // force/work reporting interval
  vector<double> moment_centroids;
  double volume;      // total volume (local)
  double force[3];    // total force  (local)
  double work;        // total work   (local)
  double moment[3];   // total moment about origin (local)
  int write_rank;     // rank to reduce force/work on for write

  // pmm stuff...
  double C[3][3];     // intertial resistance coeff
  double D[3][3];     // viscous resistance coeff

  // mrf stuff (could just use C and D containers)...
  double axis_rot[3]; // axis of rotation
  double x_rot[3];  // point on axis
  double omega_rot; // speed of rot (rad/s)

  BodyForce() {
    name       = "";
    interval   = -1;
    type       = -1;
    write_rank = 0;

    FOR_I3 FOR_J3 C[i][j] = 0.0;
    FOR_I3 FOR_J3 D[i][j] = 0.0;

    FOR_I3 axis_rot[i]    = 0.0;
    FOR_I3 x_rot[i]       = 0.0;
    omega_rot             = 0.0;

    FOR_I3 force[i]       = 0.0;
    work                  = 0.0;
    FOR_I3 moment[i]      = 0.0;
    volume                = 0.0;
  }

  void zero() {
    FOR_I3 force[i] = 0.0;
    FOR_I3 moment[i] = 0.0;
    work = 0.0;
    volume = 0.0;
  }
};

//======================================================
//======================================================

enum TimestepMode {
  DT_MODE_UNKNOWN,
  DT_MODE_CONSTANT_DT,
  DT_MODE_CONSTANT_CFL,
  DT_MODE_CFL_DT,
  DT_MODE_CFL_T,
};


class FlowSolver : public StaticSolver {
private:

  // stopping criteria handled by FlowSolver...
  int nsteps;
  double runtime_hours;
  double simtime_stop;
  TimestepMode dt_mode;
  bool b_dt_fxn;

protected:

  double dt_data[2];

  // store/restore of solution for NAN robustness...
  // Here we store 2 time steps, a recent one, and a more distant
  // one. We want to use the more distant one to recover the
  // solution because it is further from the failure and may actually
  // be recoverable...
  double * store_buf;
  double * store_buf_recent;
  int store_buf_size;
  int store_step,store_step_recent;
  double store_time,store_time_recent;

public:

  // time stepping and diagnostics

  int step;
  double time, dt;
  int check_interval;
  double wtime;
  double wtime_per_interval;

  // face indexing ranges (for rhs evaluations)...
  int n_e__i; // [0     :n_e__i) non-compact, internal faces, with no c.
  int n_c__i; // [n_e__i:n_c__i) non-compact & compact, internal faces with no c..
  int n_ec_i; // [n_c__i:n_ec_i) non-compact, internal faces, with c..
  int n_cc_i; // [n_ec_i:n_cc_i) non-compact & compact, internal faces, with c..

  int n_e__b; // [n_cc_i:n_e__b) non-compact, periodic/proc boundary faces, with no c..
  int n_c__b; // [n_e__b:n_c__b) non-compact & compact, periodic/proc boundary fa, with no c..
  int n_ec_b; // [n_c__b:n_ec_b) non-compact, per/proc boundary fa, with c...
  int n_cc_b; // [n_ec_b:n_cc_b) non-compact & compact, periodic/proc obundary fa, with c..

  // geometric structures that the specialized solvers can use
  // particularly for rhs calculations

  FaHash*    fa_hashed;
  FaNC*      ef_geom;
  FaCompact* cf_geom;
  int        nef_compressed;
  int        ncf_compressed;

  // convenience map that links the sorted extended faces to their ifa
  // counterpart (-1 if the icv0,icv1 pair does not denote a compact face)

  int* efsorted2fa;

  // measures of the symmetric part of the extended face operator

  double * cmag_cv;
  double * dcmag_hat_cv;
  double * sa_to_vol_ratio;

  MpiTimer timer;

  // scalar infrastructure for the different solvers

  //int nsc;
  //vector<double*> scalar_vec;
  //vector<string> scalar_name;
  //map<string, int> scalar_idx;
  //map<string,double*> scalar_src_map;

  int nsc_transport;      // number of transported scalars
  int nsc;                // number of total scalars..

  map<string,double*> scalar_src_map;
  vector<double*> transport_scalar_vec;
  map<string,ScalarField*> scalar_map;

  // cv-based diagonal dissipation used for
  // local timestepping...

  double * dd;

  // all flow solvers need to comprehend force data requests..

  map<string,ForceIO> fioMap;
  map<const string,int> forceZoneMap;
  double (*forceZoneBuf)[9];
  double (*momentZoneBuf)[9];

  vector<BodyForce> bodyForceVec;
  map<const string,int> bodyForceNameMap;
  int * body_force_flag;

  // conjugate heat transfer stuff...
  Cht * cht;

  // Lp stuff...
  LpTracer* lpTracer;
  LpSolid*  lpSolid;
  LpLiquid* lpLiquid;
  vector<LpInjectorTracer*> tracerInjectorVec;  // tracer injectors
  vector<LpInjector*>       solidInjectorVec;   // solid injectors
  vector<LpInjector*>       liquidInjectorVec;  // liquid injectors
  vector<CtiLiquid*>        fuelVec;            // fuels
  vector<CtiSolid*>         dustVec;            // dusts
  double grav[3];                               // gravity
  int npaocv_max;                               // maximum number of solid/liquid particle in a cv
  int lpCoupling;                               // coupling Lp with gas: NA/ONE_WAY_COUPLED/TWO_WAY_COUPLED

  FlowSolver() {

    dt_mode = DT_MODE_UNKNOWN;
    b_dt_fxn = false;
    dt_data[0] = dt_data[1] = 0.0;

    step            = 0;     registerData(step,"step",READWRITE_DATA);
    time            = 0.0;   registerData(time,"time",READWRITE_DATA);
    dt              = 0.0;   registerData(dt  ,"dt"  ,READWRITE_DATA);
    check_interval  =  1;
    wtime_per_interval = -1;

    fa_hashed       = NULL;
    ef_geom         = NULL;
    cf_geom         = NULL;
    efsorted2fa     = NULL;
    nef_compressed  = -1;
    ncf_compressed  = -1;

    n_e__i          = -1;
    n_c__i          = -1;
    n_ec_i          = -1;
    n_cc_i          = -1;
    n_e__b          = -1;
    n_c__b          = -1;
    n_ec_b          = -1;
    n_cc_b          = -1;

    cmag_cv         = NULL; registerCvData(cmag_cv,"cmag_cv",CAN_WRITE_DATA);
    dcmag_hat_cv    = NULL;
    sa_to_vol_ratio = NULL;

    nsc_transport   = 0;
    nsc             = 0;
    registerScalars();

    store_buf_size = -1;
    store_buf = NULL;
    store_buf_recent = NULL;

    // diagonal dissipation used in DT_LOCAL...
    dd      = NULL; registerCvData(dd, "dd" , CAN_WRITE_DATA);

    forceZoneBuf    = NULL;
    momentZoneBuf   = NULL;

    body_force_flag        = NULL;

    check_interval = getIntParam("CHECK_INTERVAL",1);

    cht = NULL;

    lpTracer = NULL;
    if (getBoolParam("LP_TRACER",false)) {
      lpTracer = new LpTracer();
      registerLp<LpTracerState>(lpTracer->lp,"lpt",lpTracer->np);
      registerLpData<LpTracerState>(lpTracer->lp,lpTracer->lp->icv,"lpt:icv",NO_READWRITE_DATA);
      registerLpData<LpTracerState>(lpTracer->lp,lpTracer->lp->flag,"lpt:flag",NO_READWRITE_DATA);
      registerLpData<LpTracerState>(lpTracer->lp,lpTracer->lp->xp,"lpt:xp",READWRITE_DATA);
      registerLpData<LpTracerState>(lpTracer->lp,lpTracer->lp->dp,"lpt:dp",READWRITE_DATA);
    }

    lpSolid = NULL;
    if (getBoolParam("LP_SOLID",false)) {
      lpSolid = new LpSolid();
      registerLp<LpSolidState>(lpSolid->lp,"lps",lpSolid->np);
      registerLpData<LpSolidState>(lpSolid->lp,lpSolid->lp->icv,"lps:icv",NO_READWRITE_DATA);
      registerLpData<LpSolidState>(lpSolid->lp,lpSolid->lp->flag,"lps:flag",NO_READWRITE_DATA);
      registerLpData<LpSolidState>(lpSolid->lp,lpSolid->lp->xp,"lps:xp",READWRITE_DATA);
      registerLpData<LpSolidState>(lpSolid->lp,lpSolid->lp->up,"lps:up",READWRITE_DATA);
      registerLpData<LpSolidState>(lpSolid->lp,lpSolid->lp->dp,"lps:dp",READWRITE_DATA);
      registerLpData<LpSolidState>(lpSolid->lp,lpSolid->lp->Tp,"lps:Tp",READWRITE_DATA);
      registerLpData<LpSolidState>(lpSolid->lp,lpSolid->lp->mp,"lps:mp",READWRITE_DATA);
      registerLpData<LpSolidState>(lpSolid->lp,lpSolid->lp->npar,"lps:npar",READWRITE_DATA);
    }

    lpLiquid = NULL;
    if (getBoolParam("LP_LIQUID",false)) {
      lpLiquid = new LpLiquid();
      registerLp<LpLiquidState>(lpLiquid->lp,"lpl",lpLiquid->np);
      registerLpData<LpLiquidState>(lpLiquid->lp,lpLiquid->lp->icv,"lpl:icv",NO_READWRITE_DATA);
      registerLpData<LpLiquidState>(lpLiquid->lp,lpLiquid->lp->flag,"lpl:flag",NO_READWRITE_DATA);
      registerLpData<LpLiquidState>(lpLiquid->lp,lpLiquid->lp->xp,"lpl:xp",READWRITE_DATA);
      registerLpData<LpLiquidState>(lpLiquid->lp,lpLiquid->lp->up,"lpl:up",READWRITE_DATA);
      registerLpData<LpLiquidState>(lpLiquid->lp,lpLiquid->lp->dp,"lpl:dp",READWRITE_DATA);
      registerLpData<LpLiquidState>(lpLiquid->lp,lpLiquid->lp->Tp,"lpl:Tp",READWRITE_DATA);
      registerLpData<LpLiquidState>(lpLiquid->lp,lpLiquid->lp->mp,"lpl:mp",READWRITE_DATA);
      registerLpData<LpLiquidState>(lpLiquid->lp,lpLiquid->lp->tbu,"lpl:tbu",READWRITE_DATA);
      registerLpData<LpLiquidState>(lpLiquid->lp,lpLiquid->lp->npar,"lpl:npar",READWRITE_DATA);
    }

    setupLpMaterial();

    FOR_I3 grav[i] = 0.0;
    Param * param = getParam("LP.GRAVITY");
    if (param) {
      FOR_I3 grav[i] = param->getDouble(i);
      COUT1(" > Has gravity: " << grav[0] << ", " << grav[1] << ", " << grav[2]);
    } else {
      COUT1(" > No gravity");
    }

    npaocv_max = -1;
    string token = getStringParam("LP.COUPLING","NA");
    if (token == "1WAY_COUPLED") {
      lpCoupling = ONE_WAY_COUPLED;
      if (mpi_rank==0) cout << " > One-way coupling for particles" << endl;
    }
    else if (token == "2WAY_COUPLED") {
      lpCoupling = TWO_WAY_COUPLED;
      if (mpi_rank==0) cout << " > Two-way coupling for particles" << endl;
    }
    else if (token == "NA") {
      lpCoupling = NOT_COUPLED;
      if ((lpSolid)||(lpLiquid)) {
        CERR( "Error: LP.COUPLING is not provided when particles exist. Options are:\n" <<
                      "     LP.COUPLING 1WAY_COUPLED\n" <<
                      "     LP.COUPLING 2WAY_COUPLED");
      }
    }

  }

  FlowSolver(const int icg) : StaticSolver(icg) {

    step            = 0;
    time            = 0.0;
    dt              = 0.0;
    check_interval  =  1;
    wtime_per_interval = -1;

    fa_hashed       = NULL;
    ef_geom         = NULL;
    cf_geom         = NULL;
    efsorted2fa     = NULL;
    nef_compressed  = -1;
    ncf_compressed  = -1;

    n_e__i          = -1;
    n_c__i          = -1;
    n_ec_i          = -1;
    n_cc_i          = -1;
    n_e__b          = -1;
    n_c__b          = -1;
    n_ec_b          = -1;
    n_cc_b          = -1;

    cmag_cv         = NULL;
    dcmag_hat_cv    = NULL;
    sa_to_vol_ratio = NULL;

    // assume multigrid is not used for scalar transport...
    nsc_transport   = 0;
    nsc             = 0;

    store_buf_size = -1;
    store_buf = NULL;
    store_buf_recent = NULL;

    // diagonal dissipation used in DT_LOCAL...
    dd      = NULL;

    forceZoneBuf    = NULL;
    momentZoneBuf   = NULL;

    body_force_flag        = NULL;

    check_interval = getIntParam("CHECK_INTERVAL",1);

    cht = NULL;
    lpTracer = NULL;
    lpSolid = NULL;
    lpLiquid = NULL;
  }

  void registerScalars() {

    FOR_PARAM_MATCHING("REGISTER_SCALAR") {

      int iarg = 0;
      while ( iarg < param->size() ) {

        // the registration infrastructure is going to check for duplicated
        // names so we are just going to add these scalars to the list..

        const string name = param->getString(iarg++);
        const double Sc_lam = getDoubleParam(name+".SC_LAM",0.0);
        const double Sc_t   = getDoubleParam(name+".SC_T",0.0);
        _registerScalar(name,0,Sc_lam,Sc_t);

      }
    }

  }

  void _registerScalar(const string& name,const int _kind = 0,const double _Sc_lam = 0.0,const double _Sc_t = 0.0) {

    ScalarField * sf = new ScalarField(name,_kind,_Sc_lam,_Sc_t);

    if (mpi_rank == 0)
      cout << " registering scalar with name: " << name << endl;

    registerCvData(sf->phi,name,READWRITE_DATA);

    if ( scalar_map.find(name) == scalar_map.end()) {
      scalar_map[name] = sf;
    } else {
      CERR(" > already have registered a scalar with name: " << name);
    }

    ++nsc;

  }

  virtual void initFromParams() {

    //StaticSolver::initFromParams();  // now called in StaticSolver::init()

    if (checkParam("ZERO_SCALARS"))
     resetScalars();

    // the injectors need to know about the solver's velocity
    // which is set in the solver initData function
    // we should make sure this function is called after the solver's initData call
    setupLpInjectors();

  }

  void resetScalars() {

    // resets all transported and NOT transported scalrs back to zero..

    for (int isc = 0; isc < nsc; ++isc) {
      for (int icv = 0; icv < ncv_g2; ++icv)
        transport_scalar_vec[isc][icv] = 0.0;
    }
  }

  virtual ~FlowSolver() {

    COUT1("~FlowSolver()");

    DELETE(fa_hashed);
    DELETE(ef_geom);
    DELETE(cf_geom);
    DELETE(efsorted2fa);
    DELETE(cmag_cv);
    DELETE(dcmag_hat_cv);
    DELETE(sa_to_vol_ratio);


    for(map<string,ScalarField*>::iterator it = scalar_map.begin(); it != scalar_map.end(); ++it)
      delete it->second;

    scalar_map.clear();

    DELETE(store_buf);
    DELETE(store_buf_recent);

    DELETE(dd);

    DELETE(momentZoneBuf);
    DELETE(forceZoneBuf);

    DELETE(body_force_flag);

    if (cht != NULL) delete cht;

    if (lpTracer != NULL) {
      for (int ivec = 0; ivec < lpTracer->bcVec.size(); ++ivec) {
        delete lpTracer->bcVec[ivec];
      }
      delete lpTracer;
    }
    for (int ivec = 0; ivec < tracerInjectorVec.size(); ++ivec)
      delete tracerInjectorVec[ivec];

    if (lpSolid != NULL) {
      for (int ivec = 0; ivec < lpSolid->bcVec.size(); ++ivec) {
        delete lpSolid->bcVec[ivec];
      }
      delete lpSolid;
    }
    for (int ivec = 0; ivec < solidInjectorVec.size(); ++ivec)
      delete solidInjectorVec[ivec];

    if (lpLiquid != NULL) {
      for (int ivec = 0; ivec < lpLiquid->bcVec.size(); ++ivec) {
        delete lpLiquid->bcVec[ivec];
      }
      delete lpLiquid;
    }
    for (int ivec = 0; ivec < liquidInjectorVec.size(); ++ivec)
      delete liquidInjectorVec[ivec];

    for (int ivec = 0; ivec < dustVec.size(); ++ivec)
      delete dustVec[ivec];

    for (int ivec = 0; ivec < fuelVec.size(); ++ivec)
      delete fuelVec[ivec];

  }

  virtual void init() {

    wtime = MPI_Wtime();

    StaticSolver::init();

    srand(2+mpi_rank);

    if ( checkParam("DISABLE_FA_GROUPING")) {

      for (int ief = 0; ief < nef; ++ief)
        group_ef[ief] = -1;

      for (int ifa = 0; ifa < nfa; ++ifa)
        group_fa[ifa] = -1;

    }

    dumpRange(group_ef,nef,"GROUP_EF");
    dumpRange(group_fa,nfa,"GROUP_FA");

    if ( !checkParam("USE_COMPACT_FACES")) {

      initFaStructures(n_ef,c_ef,group_ef,cvoef,nef,
                       n_fa,area_over_delta_fa,group_fa,cvofa,nfa);

    } else {

      if ( mpi_rank == 0)
        cout << " > Warning: bypassing extended face operators. " << endl;

      double (*c_fa)[3] = new double[nfa][3];
      for (int ifa = 0; ifa < nfa; ++ifa)
        for (int i = 0; i < 3; ++i)
          c_fa[ifa][i] = 0.0;

      initFaStructures(n_fa,c_fa,group_fa,cvofa,nfa,
                       n_fa,area_over_delta_fa,group_fa,cvofa,nfa);

      delete[] c_fa;

    }

    // also build cmag_cv...
    assert(cmag_cv == NULL);
    assert(dcmag_hat_cv == NULL);

    cmag_cv      = new double[ncv_g2];
    dcmag_hat_cv = new double[ncv_g2];
    FOR_ICV {
      cmag_cv[icv]      = 0.0;
      dcmag_hat_cv[icv] = 0.0;
    }

    // as a note on convention -- the stabilization contributions are scaled against
    // the row norm of || D + D^T || (stored in the c-vector), while the compressed
    // fluxes that are used at a later point consider the c-vector = 1/2 ( D+ D^T)
    // that would arise from a decomposition into symmetric and skew-symmetric parts

    FOR_INTERNAL_IEF {
      const double cmag = MAG(c_ef[ief]);
      const int icv0 = cvoef[ief][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvoef[ief][1]; assert((icv1 >= 0)&&(icv1 < ncv));
      cmag_cv[icv0] += cmag;
      cmag_cv[icv1] += cmag;

      const double delta = DIST(x_cv[icv1],x_cv[icv0]);
      dcmag_hat_cv[icv0] += delta*cmag;
      dcmag_hat_cv[icv1] += delta*cmag;
    }

    FOR_INTERPROC_IEF {
      const double cmag = MAG(c_ef[ief]);
      const int icv0 = cvoef[ief][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      cmag_cv[icv0] += cmag;

      const int icv1 = cvoef[ief][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g2));
      const double delta = DIST(x_cv[icv1],x_cv[icv0]);
      dcmag_hat_cv[icv0] += delta*cmag;
    }

    // we need to normalize the dcmag with some measure of the area...
    assert( faocv_i != NULL);
    assert( faocv_v != NULL);
    for (int icv = 0; icv < ncv; ++icv) {
      double area_sum = 0.0;
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        area_sum     += MAG(n_fa[ifa]);
      }
      assert( area_sum > 0.0);
      dcmag_hat_cv[icv] /= area_sum;
    }


    updateCv2Data(cmag_cv);
    updateCv2Data(dcmag_hat_cv);

    MiscUtils::dumpRange(cmag_cv,ncv_g2,"cmag_cv");
    MiscUtils::dumpRange(dcmag_hat_cv,ncv_g2, "dcmag_hat_cv");

    // *****************************************************
    // cmag hack -- use vol_cv differences...
    if (checkParam("DVOL_FORM")) {
      
      double *wgt_sum = new double[ncv]; 
      FOR_ICV {
	dcmag_hat_cv[icv] = 0.0;
	wgt_sum[icv] = 0.0;
      }

      // compact version...
      /*
      FOR_INTERNAL_IFA {
	const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
	const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
	const double area = MAG(n_fa[ifa]);
	const double coeff = area*fabs(vol_cv[icv1]-vol_cv[icv0])/(vol_cv[icv1]+vol_cv[icv0]);
	dcmag_hat_cv[icv0] += coeff;
	wgt_sum[icv0]      += area;
	dcmag_hat_cv[icv1] += coeff;
	wgt_sum[icv1]      += area;
      }

      FOR_INTERPROC_IFA {
	const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
	const int icv1 = cvofa[ifa][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g2));
	const double area = MAG(n_fa[ifa]);
	const double coeff = area*fabs(vol_cv[icv1]-vol_cv[icv0])/(vol_cv[icv1]+vol_cv[icv0]);
	dcmag_hat_cv[icv0] += coeff;
	wgt_sum[icv0]      += area;
      }
      */

      // ef version...
      FOR_INTERNAL_IEF {
	const int icv0 = cvoef[ief][0]; assert((icv0 >= 0)&&(icv0 < ncv));
	const int icv1 = cvoef[ief][1]; assert((icv1 >= 0)&&(icv1 < ncv));
	const double area = MAG(n_ef[ief]);
	const double coeff = area*fabs(vol_cv[icv1]-vol_cv[icv0])/(vol_cv[icv1]+vol_cv[icv0]);
	dcmag_hat_cv[icv0] += coeff;
	wgt_sum[icv0]      += area;
	dcmag_hat_cv[icv1] += coeff;
	wgt_sum[icv1]      += area;
      }

      FOR_INTERPROC_IEF {
	const int icv0 = cvoef[ief][0]; assert((icv0 >= 0)&&(icv0 < ncv));
	const int icv1 = cvoef[ief][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g2));
	const double area = MAG(n_ef[ief]);
	const double coeff = area*fabs(vol_cv[icv1]-vol_cv[icv0])/(vol_cv[icv1]+vol_cv[icv0]);
	dcmag_hat_cv[icv0] += coeff;
	wgt_sum[icv0]      += area;
      }

      FOR_ICV {
	dcmag_hat_cv[icv] = dcmag_hat_cv[icv]/wgt_sum[icv]*pow(vol_cv[icv],1.0/3.0); // works pretty well
	//dcmag_hat_cv[icv] = dcmag_hat_cv[icv]/wgt_sum[icv]; // leave non-dimensional for filtering...
      }

      /*
      updateCv2Data(dcmag_hat_cv);
      
      double *filter = new double[ncv];
      FOR_ICV filter[icv] = 0.0;

      FOR_INTERNAL_IEF {
	const int icv0 = cvoef[ief][0]; assert((icv0 >= 0)&&(icv0 < ncv));
	const int icv1 = cvoef[ief][1]; assert((icv1 >= 0)&&(icv1 < ncv));
	const double coeff = MAG(n_ef[ief])*(dcmag_hat_cv[icv0]+dcmag_hat_cv[icv1]);
	filter[icv0] += coeff;
	filter[icv1] += coeff;
      }

      FOR_INTERPROC_IEF {
	const int icv0 = cvoef[ief][0]; assert((icv0 >= 0)&&(icv0 < ncv));
	const int icv1 = cvoef[ief][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g2));
	const double coeff = MAG(n_ef[ief])*(dcmag_hat_cv[icv0]+dcmag_hat_cv[icv1]);
	filter[icv0] += coeff;
      }
      
      FOR_ICV {
	dcmag_hat_cv[icv] = filter[icv]/wgt_sum[icv]*pow(vol_cv[icv],1.0/3.0); 
      }
      
      delete[] filter;
      */

      delete[] wgt_sum; 

      updateCv2Data(dcmag_hat_cv);
      
      MiscUtils::dumpRange(dcmag_hat_cv,ncv_g2, "dcmag_hat_cv");
      
      MPI_Pause("modified dcmag_hat_cv");
      
    }

    setupFilterMoments();

    assert(sa_to_vol_ratio == NULL);
    sa_to_vol_ratio = new double[ncv];

    for (int icv =0; icv <ncv; ++icv)
        sa_to_vol_ratio[icv] = 0.0;

    for (int ifa =0; ifa < nfa; ++ifa) {

      const int icv0    = cvofa[ifa][0];
      const int icv1    = cvofa[ifa][1];

      const double nmag = MAG(n_fa[ifa]);
      sa_to_vol_ratio[icv0]       += nmag;
      if ( icv1 < ncv)
        sa_to_vol_ratio[icv1]     += nmag;

    }

    for (int icv = 0; icv < ncv; ++icv)
      sa_to_vol_ratio[icv] = sa_to_vol_ratio[icv]/vol_cv[icv];

    // eliminating cht zones...
    //if (cht) cht->setZnost(bfZoneVec.size());

    registerLpBoundaryConditions();

    // needed by app to control contex menu
    if ( (lpTracer)||(lpSolid)||(lpLiquid) )
      added_data_fields.insert("particles");

    if (mpi_rank == 0)
      logger->setKillFilename("killcharles");


  }

  virtual void initData() {

    COUT1("FlowSolver::initData()");

    // dd is always set to 1, unless the user has specified
    // local timstepping...
    assert( dd == NULL); dd = new double[ncv]; FOR_ICV dd[icv] = 1.0;

    // ========================================================
    // CHT stuff
    // ========================================================

    if (Param * param = getParam("CHT")) {

      if (mpi_rank == 0) cout << "CHT..." << endl;

      if (sspobf_i == NULL) {
        assert(sspobf_v == NULL);
        assert(sspobf_wgt == NULL);
        CERR("CHT requires sspobf_i/v/wgt in mles.\nRebuild mles with current stitch.exe");
      }
      assert(sspobf_v);
      assert(sspobf_wgt);

      assert(cht == NULL);
      cht = new Cht();

      // registration...
      CtiRegister::_registerData(cht->T_no,"cht:T_no",NO_DATA,0,cht->nno_g); // READ/WRITE evantually, but 0 for now
      CtiRegister::_registerData(cht->x_no,"cht:x_no",NO_DATA,X_DATA,cht->nno_g);
      CtiRegister::_registerName("cht",cht->nno_g);
      CtiRegister::_registerNoote(cht->noote,cht->nte_ib,cht->nno_g);

      // had to register out here because cht registration lives in CtiRegister and StaticSolver
      //registerCht("cht",cht->noote,cht->noost,cht->znost,cht->nte_ib,cht->nst_ib,cht->nno_g); // need to include all nodes for tets you own
      //registerChtData(cht->x_no,"cht:x_no");
      //registerChtData(cht->T_no,"cht:T_no");

      // as a HACK, look for any cht: in STATS. This should not be neccessary after the
      // registration is moved to the pre-sles read part of the code...
      if (Param * param = getParam("STATS")) {
        for (int iarg = 0; iarg < param->size(); ++iarg) {
          const string varname = param->getString(iarg);
          if (MiscUtils::startsWith(varname,"cht:")) {
            CtiRegister::CtiData * data = CtiRegister::getCtiData(varname,false);
            assert(data != NULL);
            CtiRegister::registerStats(varname,data,false);
          }
        }
      }

      const string mesh_filename = param->getString(0);
      if ((mesh_filename.size() > 5)&&(mesh_filename.compare(mesh_filename.size()-5,5,".tles") == 0)) {
        // NEED to rethink this format and its support for surface faces and their
        // link to boundary faces...
        assert(0);
        // *.tles file coming from stitch
        cht->read_tles(mesh_filename);
        // HACK initial condition
        assert(cht->T_no == NULL); cht->T_no = new double[cht->nno];
        for (int ino = 0; ino < cht->nno; ++ino) {
          cht->T_no[ino] = cht->x_no[ino][0]+cht->x_no[ino][1]+cht->x_no[ino][2];
        }
      }
      else if ((mesh_filename.size() > 8)&&(mesh_filename.compare(mesh_filename.size()-8,8,".tet_bin") == 0)) {
        // simple binary tet file...
        if (mpi_rank == 0) {
          FILE * fp = fopen(mesh_filename.c_str(),"rb");
          if (fp == NULL) {
            cout << "Error: cannot open cht tet_bin file: " << mesh_filename << endl;
            cht->nno = -1;
          }
          else {
            bool b_zones = false;
            fread(&cht->nno,sizeof(int),1,fp);
            if (cht->nno == -1) {
              b_zones = true;
              fread(&cht->nno,sizeof(int),1,fp);
            }
            fread(&cht->nte,sizeof(int),1,fp);
            cout << " > cht tet_bin file has nno: " << cht->nno << ", nte: " << cht->nte << " b_zones: " << b_zones << endl;
            assert(cht->x_no == NULL);
            cht->x_no = new double[cht->nno][3];
            fread(cht->x_no,sizeof(double),cht->nno*3,fp);
            assert(cht->noote == NULL);
            cht->noote = new int[cht->nte][4];
            fread(cht->noote,sizeof(int),cht->nte*4,fp);
            assert(cht->znote == NULL);
            if (b_zones) {
              cht->znote = new int[cht->nte];
              fread(cht->znote,sizeof(int),cht->nte,fp);
              set<int> zoneSet;
              for (int ite = 0; ite < cht->nte; ++ite) {
                const int izn = cht->znote[ite];
                zoneSet.insert(izn);
              }
              cout << " > cht tet_bin file has these volume zones: ";
              for (set<int>::iterator it = zoneSet.begin(); it != zoneSet.end(); ++it) {
                cout << " " << *it;
              }
              cout << endl;
            }
            fclose(fp);
            // finally allocate T...
            assert(cht->T_no == NULL); cht->T_no = new double[cht->nno];
            for (int ino = 0; ino < cht->nno; ++ino) {
              cht->T_no[ino] = cht->x_no[ino][0]+cht->x_no[ino][1]+cht->x_no[ino][2]; // HACK initial condition
            }
          }
        }
        assert(cht->nno_global == 0);
        cht->nno_global = cht->nno;
        MPI_Bcast(&cht->nno_global,1,MPI_INT8,0,mpi_comm);
        if (cht->nno_global == -1) throw(0);
        assert(cht->noora_striped == NULL);
        MiscUtils::calcThresholdDist(cht->noora_striped,cht->nno_global,mpi_size,DIST_THRESHOLD);
      }
      else {
        // assume coming from fluent...
        // read the msh on rank 0 only...
        if (mpi_rank == 0) {
          FluentMsh * msh = new FluentMsh(mesh_filename);
          //msh->writeTecplot("cht.dat");
          // for now, everything on rank 0...
          cht->nno = msh->nno;
          assert(cht->noote == NULL);
          msh->buildTets(cht->nte,cht->noote);
          assert(cht->x_no == NULL); cht->x_no = new double[msh->nno][3];
          assert(cht->T_no == NULL); cht->T_no = new double[msh->nno];
          for (int ino = 0; ino < msh->nno; ++ino) {
            FOR_I3 cht->x_no[ino][i] = msh->x_no[ino][i];
            cht->T_no[ino] = cht->x_no[ino][0]+cht->x_no[ino][1]+cht->x_no[ino][2]; // HACK initial condition
          }
          delete msh;
        }

        assert(cht->nno_global == 0);
        cht->nno_global = cht->nno;
        MPI_Bcast(&cht->nno_global,1,MPI_INT8,0,mpi_comm);
        assert(cht->noora_striped == NULL);
        MiscUtils::calcThresholdDist(cht->noora_striped,cht->nno_global,mpi_size,DIST_THRESHOLD);
      }

      cht->init(); // load balances the above striping

      // if the restart file has a thermal solution, then set it...

      if (param->size() > 1) {
        assert(param->size() == 2);
        const string soln_filename = param->getString(1);
        cht->readData(soln_filename);
      }

      // even if a thermal solution was found, allow the user to overide
      // with CHT.T_INIT...

      if (Param * param = getParam("CHT.T_INIT")) {
        const double T_init = param->getDouble();
        if (mpi_rank == 0) cout << " > applying CHT.T_INIT: " << T_init << endl;
        for (int ino = 0; ino < cht->nno; ++ino)
        cht->T_no[ino] = T_init;
        cht->updateNoData(cht->T_no);
      }

      // alternatively, user can interpolate an initial condition from previous
      // cht/mesh combo (must have both currently b/c cht:x_no not written to solution file)

      if (Param * param = getParam("CHT.INTERP")) {
        COUT1(" > processing CHT.INTERP parameters");

        // build bboxes based on my loaded cht mesh, using edge lengths to define
        // bbox extents

        COUT1("    > determining local max edge lengths");
        // want to build adt bounding boxes based on largest touching edge length;
        // traverse noote and store largest edge for each touching node
        // use ghosts to determine largest touching edge, but only store info for local nodes
        double * max_d2 = new double[cht->nno_g];  // largest edge including ghost nodes
        for (int ino=0; ino<cht->nno_g; ++ino) max_d2[ino] = -HUGE_VAL;

        for (int ite=0; ite<cht->nte; ++ite) {
          const double d2_01 = DIST2(cht->x_no[cht->noote[ite][0]],cht->x_no[cht->noote[ite][1]]);
          const double d2_02 = DIST2(cht->x_no[cht->noote[ite][0]],cht->x_no[cht->noote[ite][2]]);
          const double d2_03 = DIST2(cht->x_no[cht->noote[ite][0]],cht->x_no[cht->noote[ite][3]]);
          const double d2_12 = DIST2(cht->x_no[cht->noote[ite][1]],cht->x_no[cht->noote[ite][2]]);
          const double d2_13 = DIST2(cht->x_no[cht->noote[ite][1]],cht->x_no[cht->noote[ite][3]]);
          const double d2_23 = DIST2(cht->x_no[cht->noote[ite][2]],cht->x_no[cht->noote[ite][3]]);

          double d2_tmp = max(d2_01,max(d2_02,d2_03));
          if (d2_tmp > max_d2[cht->noote[ite][0]]) max_d2[cht->noote[ite][0]] = d2_tmp;

          d2_tmp = max(d2_01,max(d2_12,d2_13));
          if (d2_tmp > max_d2[cht->noote[ite][1]]) max_d2[cht->noote[ite][1]] = d2_tmp;

          d2_tmp = max(d2_02,max(d2_12,d2_23));
          if (d2_tmp > max_d2[cht->noote[ite][2]]) max_d2[cht->noote[ite][2]] = d2_tmp;

          d2_tmp = max(d2_03,max(d2_13,d2_23));
          if (d2_tmp > max_d2[cht->noote[ite][3]]) max_d2[cht->noote[ite][3]] = d2_tmp;
        }

        // only build adt for local nodes (not ghosts)
        double (*xno_min)[3] = new double[cht->nno][3];
        double (*xno_max)[3] = new double[cht->nno][3];
        double my_max = -HUGE_VAL;
        for (int ino=0; ino<cht->nno; ++ino) {
          assert(max_d2[ino] >= 0.0);
          my_max = max(my_max,sqrt(max_d2[ino]));
          const double delta = 1.0*sqrt(max_d2[ino]);  // edge length should overlap enough
          for (int i=0; i<3; ++i) {
            xno_min[ino][i] = cht->x_no[ino][i] - delta;
            xno_max[ino][i] = cht->x_no[ino][i] + delta;
          }
        }
        DELETE(max_d2);
        MiscUtils::dumpRange(&my_max,1,"maximum tet edge length");

        COUT1("    > building cht Adt based on edge lengths");

        Adt<double> * cht_adt = new Adt<double>(cht->nno,xno_min,xno_max);
        DELETE(xno_min);
        DELETE(xno_max);

        // build rank bbox adt from my nodes...
        double my_bbmin[3],my_bbmax[3];
        cht_adt->getBbox(my_bbmin,my_bbmax);
        double (*bbmin)[3] = new double[mpi_size][3];
        double (*bbmax)[3] = new double[mpi_size][3];
        MPI_Allgather(my_bbmin,3,MPI_DOUBLE,bbmin,3,MPI_DOUBLE,mpi_comm);
        MPI_Allgather(my_bbmax,3,MPI_DOUBLE,bbmax,3,MPI_DOUBLE,mpi_comm);
        Adt<double> * cht_bbox_adt = new Adt<double>(mpi_size,bbmin,bbmax);
        DELETE(bbmin);
        DELETE(bbmax);

        int * xno_flag = new int[cht->nno];
        for (int ino=0; ino<cht->nno; ++ino) {
          // flag whether node has been set
          // -2: unset, -1: previously set by interp, 0+ current index
          xno_flag[ino] = -2;
        }

        // rank sized, so can allocate outside
        int * send_count = new int[mpi_size];
        int * send_disp = new int[mpi_size];
        int * recv_count = new int[mpi_size];
        int * recv_disp = new int[mpi_size];

        int nodes_set = 0;  // count how many nodes populated by interp files

        // NOTE: current behavior is that if a node has been set from an interp,
        // subsequent interp will not over-write those values

        FOR_PARAM_MATCHING("CHT.INTERP") {
          // need both mesh/solution for now
          int iarg=0;
          const string meshfile = param->getString(iarg++);
          if (meshfile.compare(meshfile.size()-8,8,".tet_bin") != 0) {
            WUI(WARN,"CHT.INTERP mesh file must be *.tet_bin; skipping mesh: " << meshfile);
            break;
          }

          const string chtfile  = param->getString(iarg++);
          if (mpi_rank == 0) {
            WUI(INFO," > interpolating cht:T_no from: " << meshfile << " " << chtfile);
          }

          // rank 0 reads header and distributes
          int interp_nno_global = -1;
          int mesh_b_zones = 0;

          if (mpi_rank == 0) {
            FILE * fpm = fopen(meshfile.c_str(),"rb");
            if (fpm == NULL) {
              cout << "Error: cannot open cht tet_bin file: " << meshfile << endl;
            }
            else {
              FILE * fpc = fopen(chtfile.c_str(),"rb");
              if (fpc == NULL) {
                cout << "Error: cannot open cht solution file: " << chtfile << endl;
                fclose(fpm);
              }
              else {
                COUT1(" > processing interpolation file headers");
                // read nno header only
                fread(&interp_nno_global,sizeof(int),1,fpm);
                if (interp_nno_global == -1) {
                  mesh_b_zones = 1;
                  fread(&interp_nno_global,sizeof(int),1,fpm);
                }
                fclose(fpm);

                int solution_nno_global = -1;  // should this be int8?

                fread(&solution_nno_global,sizeof(int),1,fpc);  // int8?
                fclose(fpc);

                if (solution_nno_global != interp_nno_global) {
                  CWARN("nno between tet_bin and .cht do not match; skipping");
                  interp_nno_global = -1;  // set to not process these files
                }
                else {
                  COUT1(" > mesh/cht has nno: " << interp_nno_global);
                }
              }
            }
          }  // rank0

          MPI_Bcast(&interp_nno_global,1,MPI_INT,0,mpi_comm);
          if (interp_nno_global == -1) {
            CWARN(" > interpolation files have issues; ignoring");
            break;  // move to next INTERP command
          }
          MPI_Bcast(&mesh_b_zones,1,MPI_INT,0,mpi_comm);

          // determine how many nodes read on each rank
          int8 * noora = NULL;
          MiscUtils::calcUniformDist(noora,interp_nno_global,mpi_size);

          // local nodes (only ones that we read)
          const int interp_nno = noora[mpi_rank+1]-noora[mpi_rank];

          int interp_nte = -1;
          double (*interp_x_no)[3] = new double[interp_nno][3];

          char dummy[128]; assert(meshfile.length() < 128);
          sprintf(dummy,"%s",meshfile.c_str());
          MPI_File fh;
          MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

          // header offsetting
          MPI_Offset offset = int_size;  // nno/zone info
          if (mesh_b_zones) offset += int_size;
          offset += int_size;  // nte

          bool byte_swap = false;  // nothing specified (magic number) so no way to know...

          // x,y,z
          readChunkedData<double>(fh,offset+(noora[mpi_rank]*double_size*3),interp_x_no[0],interp_nno*3,byte_swap,mpi_comm);
          MPI_File_close(&fh);
          MiscUtils::dumpRange(interp_x_no,interp_nno,"interpolated x_no");

          // read T record
          sprintf(dummy,"%s",chtfile.c_str());
          MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

          double * interp_T_no = new double[interp_nno];
          // header offsetting
          offset = int8_size;  // nno
          readChunkedData<double>(fh,offset+noora[mpi_rank]*double_size,interp_T_no,interp_nno,byte_swap,mpi_comm);
          MPI_File_close(&fh);
          MiscUtils::dumpRange(interp_T_no,interp_nno,"interpolated T_no");
          DELETE(noora);

          // find which rank(s) to send points to based on bbox adt
          FOR_RANK send_count[rank] = 0;

          // TODO: right now hard coded variable count based on only x_no, T_no
          // available as CHT vars. Once more vars available need to do more elegantly
          vector<int> intVec;
          for (int iino = 0; iino < interp_nno; ++iino) {
            assert(intVec.empty());
            cht_bbox_adt->buildListForPoint(intVec,interp_x_no[iino]);
            for (int ii = 0, ii_end=intVec.size(); ii < ii_end; ++ii) {
              const int rank = intVec[ii]; assert((rank >= 0)&&(rank < mpi_size));
              send_count[rank] += 4;  // x_no, T_no (4*double)
            }
            intVec.clear();
          }

          // see how much info we will send out
          send_disp[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank) {
            send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
          }
          const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

          MiscUtils::dumpRange(&send_count_sum,1,"packed data size");

          // allocate and pack the buffer with send data
          double * send_buf = new double[send_count_sum];
          for (int iino = 0; iino < interp_nno; ++iino) {
            assert(intVec.empty());
            cht_bbox_adt->buildListForPoint(intVec,interp_x_no[iino]);
            for (int ii = 0, ii_end=intVec.size(); ii < ii_end; ++ii) {
              const int rank = intVec[ii]; assert((rank >= 0)&&(rank < mpi_size));
              send_buf[send_disp[rank]  ] = interp_x_no[iino][0];
              send_buf[send_disp[rank]+1] = interp_x_no[iino][1];
              send_buf[send_disp[rank]+2] = interp_x_no[iino][2];
              send_buf[send_disp[rank]+3] = interp_T_no[iino];
              send_disp[rank] += 4;
            }
            intVec.clear();
          }
          DELETE(interp_x_no);
          DELETE(interp_T_no);

          // rewind...
          send_disp[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank) {
            send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
          }

          // now send...
          MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

          recv_disp[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank) {
            recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
          }
          const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

          MiscUtils::dumpRange(&recv_count_sum,1,"unpacked data size");

          double * recv_buf = new double[recv_count_sum];
          MPI_Alltoallv(send_buf,send_count,send_disp,MPI_DOUBLE,recv_buf,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
          DELETE(send_buf);

          // now unpack and set closest node index...
          COUT1(" > finding closest interpolation nodes");

          double * min_dist2 = new double[cht->nno];
          for (int ino=0; ino<cht->nno; ++ino) min_dist2[ino] = HUGE_VAL;

          int * irono = new int[cht->nno];
          for (int ino=0; ino<cht->nno; ++ino) irono[ino] = -1;
          for (int irecv = 0; irecv < recv_count_sum; irecv += 4) {  // increment by 4 doubles (x_no, T_no)
            double xp[3]; FOR_I3 xp[i] = recv_buf[irecv+i];
            assert(intVec.empty());
            cht_adt->buildListForPoint(intVec,xp);
            for (int ii = 0, ii_end=intVec.size(); ii < ii_end; ++ii) {
              const int ino = intVec[ii]; assert((ino >= 0)&&(ino < cht->nno));
              const double dist2 = DIST2(xp,cht->x_no[ino]);
              if (dist2 < min_dist2[ino]) {
                min_dist2[ino] = dist2;
                irono[ino] = irecv;
              }
            }
            intVec.clear();
          }

          // now set data values...
          COUT1(" > setting data at nodes");
          int8 my_count[2] = { cht->nno, 0 };
          for (int ino=0; ino<cht->nno; ++ino) {
            const int irecv = irono[ino];
            if (irecv >= 0) {
              assert(irecv < recv_count_sum);

              if (xno_flag[ino] != -1) {
                cht->T_no[ino] = recv_buf[irecv+3];
                xno_flag[ino] = -1;  // has been set
              }
            }

            if (xno_flag[ino] == -1) ++my_count[1];  // counts points that have been set from previous interp command
          }
          DELETE(irono);
          DELETE(min_dist2);
          DELETE(recv_buf);

          int8 count[2];
          MPI_Allreduce(my_count,count,2,MPI_INT8,MPI_SUM,mpi_comm);  // everybody needs to know so they can decide to march or not later
          COUT1(" > set " << count[1] << " out of " << count[0] << " cht nodes");
          nodes_set = count[1];
        }
        delete cht_adt;
        delete cht_bbox_adt;

        if (nodes_set > 0) {
          // all interp params done; now perform march of T
          if (mpi_rank == 0) cout << " > marching solution";  // will over-write T_INIT
          assert(cht->noono_i);
          assert(cht->noono_v);
          double * min_dist2 = new double[cht->nno_g];
          for (int ino=0; ino < cht->nno; ++ino) {
            if (xno_flag[ino] == -1) min_dist2[ino] = -1.0;  // data set already
            else min_dist2[ino] = HUGE_VAL;
          }

          // poppulate data at ghosts
          cht->updateNoData(cht->T_no);
          cht->updateNoData(min_dist2);  // populate min_dist2 at ghosts too

          int done = 0;
          while (done == 0) {

            if (mpi_rank == 0) {
              cout << ".";
              cout.flush();
            }

            int my_done = 1;

            for (int ino=0; ino < cht->nno; ++ino) {
              if (min_dist2[ino] == HUGE_VAL) {
                // this node needs to be set
                // loop node neighbors and take value of closest set nbr
                for (int non = cht->noono_i[ino]; non != cht->noono_i[ino+1]; ++non) {
                  const int ino_nbr = cht->noono_v[non];
                  assert((ino_nbr >= 0)&&(ino_nbr < cht->nno_g));

                  // only pull data from set nodes (value = -1.0)
                  if (min_dist2[ino_nbr] < 0.0) {
                    const double dist2 = DIST2(cht->x_no[ino],cht->x_no[ino_nbr]);
                    if (dist2 < min_dist2[ino]) {
                      min_dist2[ino] = dist2;
                      cht->T_no[ino] = cht->T_no[ino_nbr];
                    }
                  }
                }
              }
            }  // ino

            // now count to see who has been set; mark marched values
            for (int ino=0; ino < cht->nno; ++ino) {
              if (min_dist2[ino] == HUGE_VAL) {
                my_done = 0;  // not done yet
              }
              else if (min_dist2[ino] >= 0.0) {
                min_dist2[ino] = -1.0;  // already set criteria
              }
              else {
                assert(min_dist2[ino] == -1.0);
              }
            }
            MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);

            cht->updateNoData(min_dist2);
            cht->updateNoData(cht->T_no);
          }
          if (mpi_rank == 0) cout << "done" << endl;
          DELETE(min_dist2);
        }
        else {
          CWARN(" > no cht nodes set during CHT.INTERP; geoemtries do no overlap");
        }

        DELETE(xno_flag);
        DELETE(send_count);
        DELETE(send_disp);
        DELETE(recv_count);
        DELETE(recv_disp);

        // check that initial condition has been applied somehow for all nodes
        int my_ierr = 0;
        for (int ino = 0; ino < cht->nno; ++ino) {
          if (cht->T_no[ino] == HUGE_VAL) {
            my_ierr = -1;
          }
        }
        int ierr;
        MPI_Allreduce(&my_ierr,&ierr,1,MPI_INT,MPI_MIN,mpi_comm);
        if (ierr != 0) {
          CERR(" > use cht solution file, CHT.T_INIT, or CHT.INTERP to provide initial temperature.");
        }
      }
    }

    // ========================================================
    // Tracer stuff
    // ========================================================

    if (lpTracer)
      lpTracer->resize(lpTracer->np);

  }

  void initScalars() {

    // we have not set the transported scalar count yet... we need
    // to link scalar_vec to allocated scalar fields that are going
    // to be transported.

    if ( mpi_rank == 0)
      cout << "FlowSolver::initScalars() " << endl;

    assert( nsc_transport == 0);

    for ( map<string,ScalarField*>::iterator it = scalar_map.begin();
          it != scalar_map.end(); ++it) {

      // allocate the scalar phi data ..

      ScalarField * sf = it->second;

      const int kind = sf->kind;
      assert ( sf->phi == NULL);

      sf->phi     =  new double[ncv_g2];

      if ( kind == scalar_transport_kind) {

        const int ii         = transport_scalar_vec.size();
        sf->idx              = ii;
        transport_scalar_vec.push_back(sf->phi);

      }
    }

    nsc_transport = transport_scalar_vec.size();

    if ( mpi_rank == 0)
      cout << "> transporting passive scalars : " << nsc_transport << endl;

  }

  void updateScalars() {

    // if the number of scalars became large, we should probably define
    // a communicator structure for the scalars to communicate them into the ghosts.
    // for now we will just issue a series of blocking update calls for the scalar data

    for (int isc =0; isc < nsc_transport; ++isc)
      updateCv2Data( transport_scalar_vec[isc]);

  }

  void dumpRangeScalars() const {

    for (map<string,ScalarField*>::const_iterator it = scalar_map.begin();it != scalar_map.end(); ++it) {

      ScalarField * sf = it->second;
      MiscUtils::dumpRange(sf->phi,ncv,sf->name);

    }

  }


  int getScalarIndex(const string& name) const {

    map<string,ScalarField*>::const_iterator it = scalar_map.find(name);
    if ( it != scalar_map.end())
      return it->second->idx;
    else
      return -1;

  }

  string getScalarName(const int isc) const {

    for (map<string,ScalarField*>::const_iterator it = scalar_map.begin();
         it != scalar_map.end(); ++it) {

      if ( it->second->idx == isc)
        return it->second->name;

    }

    return "unknown_scalar_name";

  }


  void initMin() {

    wtime = MPI_Wtime();

    // the min init inside of the post mode for now will not include
    // any of the extended face infrastructure or the other compressed
    // face structures used to advance the solution.  this is just the
    // bare bones.

    StaticSolver::init((INIT_COMPACT_FACES | INIT_CV_GRAD));

    // also some of the solvers will use dcmag_hat_cv.. need to make sure
    // that this is not-null, but currently empty

    assert( dcmag_hat_cv == NULL);
    dcmag_hat_cv = new double[ncv_g2];
    for (int icv = 0; icv < ncv_g2; ++icv)
      dcmag_hat_cv[icv] = 0.0;

    assert(sa_to_vol_ratio == NULL);
    sa_to_vol_ratio = new double[ncv];

    for (int icv =0; icv <ncv; ++icv)
        sa_to_vol_ratio[icv] = 0.0;

    for (int ifa =0; ifa < nfa; ++ifa) {

      const int icv0    = cvofa[ifa][0];
      const int icv1    = cvofa[ifa][1];

      const double nmag = MAG(n_fa[ifa]);
      sa_to_vol_ratio[icv0]       += nmag;
      if ( icv1 < ncv)
        sa_to_vol_ratio[icv1]     += nmag;

    }

    for (int icv = 0; icv < ncv; ++icv)
      sa_to_vol_ratio[icv] = sa_to_vol_ratio[icv]/vol_cv[icv];

    // eliminating cht zones...
    //if (cht) cht->setZnost(bfZoneVec.size());

    if ( mpi_rank == 0 )
      logger->setKillFilename("killcharles");

  }

  void checkGclExtendedHashed(const int _nef) const {
    double (*gcl)[3] = new double[ncv][3];
    for (int icv = 0; icv < ncv; ++icv) {
      for (int i =0; i < 3; ++i)
        gcl[icv][i] = 0.0;
    }

    for (int ief = 0; ief < _nef; ++ief) {
      const int icv0 = fa_hashed[ief].cvofa[0];
      const int icv1 = fa_hashed[ief].cvofa[1];
      const int ief_compressed = fa_hashed[ief].ief_compressed;

      for (int i =0; i < 3; ++i)
        gcl[icv0][i] += ef_geom[ief_compressed].n[i];

      // not assuming anything about the ordering yet..
      if ( icv1 < ncv)
        for (int i =0; i < 3; ++i)
          gcl[icv1][i] -= ef_geom[ief_compressed].n[i];
    }

    // need to close the boundary terms as well..
    FOR_IZONE(bfZoneVec) {
      const int nbf = bfZoneVec[izone].nbf;
      for (int ibf= 0; ibf < nbf; ++ibf) {
        const int icv = bfZoneVec[izone].cvobf[ibf];
        for (int i =0; i < 3; ++i)
          gcl[icv][i] += bfZoneVec[izone].n_bf[ibf][i];
      }
    }

    dumpRange(gcl,ncv,"gcl (extended n) [should be zero]");
    delete[] gcl;
  }

  void checkGclCompactHashed(const int _nef) const {
    double (*gcl)[3] = new double[ncv][3];
    for (int icv = 0; icv < ncv; ++icv) {
      for (int i =0; i < 3; ++i)
        gcl[icv][i] = 0.0;
    }

    // not assuming anything about the ordering yet...
    for (int ief = 0; ief < _nef; ++ief) {
      const int icf_compressed = fa_hashed[ief].icf_compressed;
      if ( icf_compressed >= 0 ) {
        const int icv0 = fa_hashed[ief].cvofa[0];
        const int icv1 = fa_hashed[ief].cvofa[1];

        for (int i=0; i < 3; ++i)
          gcl[icv0][i] += cf_geom[icf_compressed].unit_n[i]*cf_geom[icf_compressed].area;

        if ( icv1 < ncv)
          for (int i =0; i < 3; ++i)
            gcl[icv1][i] -= cf_geom[icf_compressed].unit_n[i]*cf_geom[icf_compressed].area;
      }
    }

    // need to close the boundary terms as well..
    for (int ibf = 0; ibf < nbf; ++ibf) {
      const int icv = cvobf[ibf];
      for (int i =0; i < 3; ++i)
        gcl[icv][i] += n_bf[ibf][i];
    }

    dumpRange(gcl,ncv,"gcl (compact n) [should be zero]");
    delete[] gcl;
  }

  void checkLinearGrad() {

    const double g[3] = {1.1,2.3,4.5};
    double *phi       = new double[ncv_g2];
    for (int icv = 0; icv < ncv_g2; ++icv)
      phi[icv] = DOT_PRODUCT(x_cv[icv],g);

    double (*gphi)[3] = new double[ncv][3];
    for (int icv = 0; icv < ncv; ++icv)
      for (int i =0; i < 3; ++i)
        gphi[icv][i] = 0.0;

    for (int ief =0; ief < n_c__i; ++ief) { // no c contributions..
      const int icv0           = fa_hashed[ief].cvofa[0];
      const int icv1           = fa_hashed[ief].cvofa[1];
      const int ief_compressed = fa_hashed[ief].ief_compressed;
      const double phi_avg     = 0.5*(phi[icv0]+phi[icv1]);

      for (int i =0; i < 3; ++i)
        gphi[icv0][i] += phi_avg*ef_geom[ief_compressed].n[i];

      for (int i =0; i < 3; ++i)
        gphi[icv1][i] -= phi_avg*ef_geom[ief_compressed].n[i];
    }

    for (int ief = n_c__i; ief < n_cc_i; ++ief) {
      const int icv0           = fa_hashed[ief].cvofa[0];
      const int icv1           = fa_hashed[ief].cvofa[1];
      const int ief_compressed = fa_hashed[ief].ief_compressed;
      const double phi_avg     = 0.5*(phi[icv0]+phi[icv1]);
      const double phi_diff    = phi[icv1]-phi[icv0];

      for (int i =0; i < 3; ++i)
        gphi[icv0][i] += phi_avg*ef_geom[ief_compressed].n[i] + phi_diff*ef_geom[ief_compressed].c[i];

      for (int i =0; i < 3; ++i)
        gphi[icv1][i] -= phi_avg*ef_geom[ief_compressed].n[i] + phi_diff*ef_geom[ief_compressed].c[i];
    }

    for (int ief = n_cc_i; ief < n_c__b ; ++ief) { // no c contributions..
      const int icv0           = fa_hashed[ief].cvofa[0];
      const int icv1           = fa_hashed[ief].cvofa[1];
      const int ief_compressed = fa_hashed[ief].ief_compressed;
      const double phi_avg     = 0.5*(phi[icv0]+phi[icv1]);

      for (int i =0; i < 3; ++i)
        gphi[icv0][i] += phi_avg*ef_geom[ief_compressed].n[i];
    }

    for (int ief = n_c__b; ief < n_cc_b; ++ief) {
      const int icv0           = fa_hashed[ief].cvofa[0];
      const int icv1           = fa_hashed[ief].cvofa[1];
      const int ief_compressed = fa_hashed[ief].ief_compressed;
      const double phi_avg     = 0.5*(phi[icv0]+phi[icv1]);
      const double phi_diff    = phi[icv1]-phi[icv0];

      for (int i =0; i < 3; ++i)
        gphi[icv0][i] += phi_avg*ef_geom[ief_compressed].n[i] + phi_diff*ef_geom[ief_compressed].c[i];
    }

    for (int ibf = 0; ibf < nbf; ++ibf) {
      const int icv = cvobf[ibf];
      const double phi_bf = DOT_PRODUCT(x_bf[ibf],g);
      for (int i =0; i < 3; ++i)
        gphi[icv][i] += phi_bf*n_bf[ibf][i];
    }

    for (int icv = 0; icv < ncv; ++icv) {

      for (int i =0; i < 3; ++i) {
        gphi[icv][i] *= inv_vol[icv];
        gphi[icv][i] /= g[i];
        gphi[icv][i] -= 1.0;
      }
    }

    dumpRange(gphi,ncv,"linear grad phi error (normalized; should order one)");

    delete[] phi;
    delete[] gphi;
  }

  void initFaStructures(const double (*_n_ef)[3], const double (*_c_ef)[3],
                        const int * _group_ef, const int (*_cvoef)[2], const int _nef,
                        const double (*_n_fa)[3], const double * _area_over_delta_fa,
                        const int * _group_fa, const int (*_cvofa)[2], const int _nfa) {

    // at this stage, the data should be distributed in memory..
    // padding the end of the array in case we are issuing software
    // prefetching instructions...
    assert(fa_hashed == NULL);
    fa_hashed = new FaHash[_nef+PFD_DISTANCE];

    // build the compressed extended face array cases...
    map<int,int> group_map;
    map<pair<int,int>,int> cvoef_map;
    nef_compressed = 0;
    int nef_compressed_check = -1;
    int nef_compressed_group =  0;
    for (int iter = 0; iter < 2; ++iter) {
      for (int ief = 0; ief < _nef; ++ief) {
        if ( iter == 0 ) {
          if ( _group_ef[ief] >= 0) {
            group_map.insert(pair<int,int>(_group_ef[ief],ief));
          } else {
            ++nef_compressed;
          }
        } else {
          assert( iter == 1);
          if ( _group_ef[ief] >= 0) {
            // geometry of this face doesn't need to be stored since
            // this is a compressed face..
            map<int,int>::iterator it = group_map.find(_group_ef[ief]);
            assert( it != group_map.end());
            fa_hashed[ief].cvofa[0]       = _cvoef[ief][0];
            fa_hashed[ief].cvofa[1]       = _cvoef[ief][1];
            fa_hashed[ief].ief_compressed = it->second;
            fa_hashed[ief].icf_compressed = -1; // tbd still..
          } else {
            // this is a incompressible face (dont worry, there's
            // no poisson equation--ha!), so it's geometry needs
            // to be stored..
            for (int i =0; i < 3; ++i) {
              ef_geom[nef_compressed].n[i] = _n_ef[ief][i];
              ef_geom[nef_compressed].c[i] = _c_ef[ief][i];
            }
            fa_hashed[ief].cvofa[0]       = _cvoef[ief][0];
            fa_hashed[ief].cvofa[1]       = _cvoef[ief][1];
            fa_hashed[ief].ief_compressed = nef_compressed;
            fa_hashed[ief].icf_compressed = -1; // tbd still..
            ++nef_compressed;
          }

          // store the cvoef pair for later, ensure that there
          // are no more duplicated entries of the cvoef..
          pair<int,int> cvoef_pair(_cvoef[ief][0],_cvoef[ief][1]);
          pair<map<pair<int,int>,int>::iterator,bool> ret;
          ret = cvoef_map.insert(pair<pair<int,int>,int>(cvoef_pair,ief));
          assert( ret.second == true);
        }
      }

      if ( iter == 0 ) {
        nef_compressed += int(group_map.size()); // total number of compressed faces
        ef_geom = new FaNC[nef_compressed];
        nef_compressed_check = nef_compressed;

        // set a local ordering for the group map, based on the original
        // lexical ordering of the groups..
        nef_compressed = 0;
        for(map<int,int>::iterator it = group_map.begin(); it != group_map.end(); ++it) {
          const int ief = it->second;
          for (int i =0; i < 3; ++i) {
            ef_geom[nef_compressed].n[i] = _n_ef[ief][i];
            ef_geom[nef_compressed].c[i] = _c_ef[ief][i];
          }
          // update the group_map to store the new ief_compressed
          it->second = nef_compressed;
          ++nef_compressed;
        }
        nef_compressed_group = int(group_map.size());
        assert( nef_compressed == int(group_map.size()));
      }
    }

    assert( nef_compressed == nef_compressed_check);

    // we now need to build the companion compressed structures
    // for the compact faces.  in addition, the cvofa needs to be
    // linked with the ef structures..
    // unlike the infrastructure, we are currently using "cf" to
    // denote the compact face to avoid confusion with the "ef"
    group_map.clear();
    ncf_compressed = 0;
    int ncf_compressed_check = -1;
    int ncf_compressed_group =  0;
    for (int iter = 0; iter < 2; ++iter) {
      for (int ifa = 0; ifa < _nfa; ++ifa) {
        if ( iter == 0 ) {
          if ( _group_fa[ifa] >= 0) {
            group_map.insert(pair<int,int>(_group_fa[ifa],ifa));
          } else {
            ++ncf_compressed;
          }
        } else {
          assert( iter == 1);

          // in order to coalesce some of the cvo{ef,fa} accesses,
          // we are storing the compact faces inside of the same
          // hashed fa structure.  we need to locate our counterpart
          // in the ef loop here ..
          map<pair<int,int>,int>::iterator it_ef
            = cvoef_map.find(pair<int,int>(_cvofa[ifa][0],_cvofa[ifa][1]));
          int ief, orientation;
          if ( it_ef != cvoef_map.end()) {
            // we found the entry that is oriented in the same direction..
            ief         = it_ef->second;
            orientation = 1;
            assert( (ief >= 0) && (ief < _nef));
            assert( _cvoef[ief][0] == _cvofa[ifa][0]);
            assert( _cvoef[ief][1] == _cvofa[ifa][1]);
          } else {
            // check the opposite orientation of the face..
            it_ef = cvoef_map.find(pair<int,int>(_cvofa[ifa][1],_cvofa[ifa][0]));
            assert( it_ef != cvoef_map.end()); // currently assuming it exists..
            ief = it_ef->second;
            orientation = -1;
            assert( (ief >= 0) && (ief < _nef));
            assert( _cvoef[ief][1] == _cvofa[ifa][0]);
            assert( _cvoef[ief][0] == _cvofa[ifa][1]);
          }

          assert( fa_hashed[ief].icf_compressed == -1);
          assert( (ief >= 0)&& (ief < _nef));

          if ( _group_fa[ifa] >= 0) {
            // geometry of this face doesn't need to be stored since
            // this is a compressed face..
            map<int,int>::iterator it = group_map.find(_group_fa[ifa]);
            assert( it != group_map.end());
            assert( orientation != -1); // both the flipped bits should be here..
            fa_hashed[ief].icf_compressed = it->second;
          } else {
            // this is a incompressible compact face
            const double area     = MAG(_n_fa[ifa]);
            const double inv_area = 1.0/area;
            if ( orientation == 1) {
              for (int i =0; i < 3; ++i)
                cf_geom[ncf_compressed].unit_n[i]     = _n_fa[ifa][i]*inv_area;
            } else {
              assert( orientation == -1);
              // need to flip the face orientation
              for (int i =0; i < 3; ++i)
                cf_geom[ncf_compressed].unit_n[i]     = -_n_fa[ifa][i]*inv_area;
            }
            cf_geom[ncf_compressed].area            = area;
            cf_geom[ncf_compressed].area_over_delta = _area_over_delta_fa[ifa];

            fa_hashed[ief].icf_compressed = ncf_compressed;
            ++ncf_compressed;
          }
        }
      }

      if ( iter == 0 ) {
        ncf_compressed += int(group_map.size()); // total number of compressed faces
        cf_geom = new FaCompact[ncf_compressed];
        ncf_compressed_check = ncf_compressed;

        // set a local ordering for the group map, based on the original
        // lexical ordering of the groups..
        ncf_compressed = 0;
        for(map<int,int>::iterator it = group_map.begin(); it != group_map.end(); ++it) {
          const int ifa         = it->second;
          const double area     = MAG(_n_fa[ifa]);
          const double inv_area = 1.0/area;
          for (int i =0; i < 3; ++i)
            cf_geom[ncf_compressed].unit_n[i]     = _n_fa[ifa][i]*inv_area;
          cf_geom[ncf_compressed].area            = area;
          cf_geom[ncf_compressed].area_over_delta = _area_over_delta_fa[ifa];

          //update the group_map to store the new icf_compressed
          it->second = ncf_compressed;
          ++ncf_compressed;
        }
        ncf_compressed_group = int(group_map.size());
        assert( ncf_compressed == int(group_map.size()));
      }//iter==0
    }//for iter

    assert( ncf_compressed == ncf_compressed_check);

    // these values should never be accessed... but
    // need to be present in case of any sftwre prefetching..
    for (int ief = _nef; ief < _nef + PFD_DISTANCE; ++ief) {
      fa_hashed[ief].cvofa[0] = 0;
      fa_hashed[ief].cvofa[1] = 0;
      fa_hashed[ief].ief_compressed = 0;
      fa_hashed[ief].icf_compressed = 0;
    }

    checkGclExtendedHashed(_nef);
    checkGclCompactHashed(_nef);

    // now we need to sort the faces into the 8 different ranges..
    // first sort the faces into their internal and processor boundary faces..
    // this sorts on the value of icv1 in the cvofa index..
    if ( mpi_rank == 0)
      cout << " ... sorting internal vs processor boundaries ... " << endl;
    sort(fa_hashed,fa_hashed+_nef,cmpFaHashCv1());

    n_cc_i = _nef; // this is the first of the interprocessor/periodic boundaries..
    for (int ief = 0; ief < _nef; ++ief) {
      if ( (n_cc_i == _nef)&&(fa_hashed[ief].cvofa[1] >= ncv)) {
        n_cc_i = ief;
        // could break here but we'll do some checking
      } else if ( n_cc_i != _nef) {
        // these entries should all be ghosts..
        assert((fa_hashed[ief].cvofa[0] >= 0) &&(fa_hashed[ief].cvofa[0] < ncv));
        assert((fa_hashed[ief].cvofa[1] >= 0) &&(fa_hashed[ief].cvofa[1] < ncv_g2));
      } else {
        // these are processor internal faces..
        assert((fa_hashed[ief].cvofa[0] >= 0)&&(fa_hashed[ief].cvofa[0] < ncv));
        assert((fa_hashed[ief].cvofa[1] >= 0)&&(fa_hashed[ief].cvofa[1] < ncv));
      }
    }

    // in the two existing ranges, we need to subdivide the regions into
    // components that have non-zero c vs zero-c.  this is identified by
    // having a group_ef[ief] >= 0, so sort on the ief_compressed index..
    if ( mpi_rank == 0 )
      cout << " ... sorting compressed face regions... " << endl;
    sort(fa_hashed,fa_hashed+n_cc_i,cmpFaHashIef());
    sort(fa_hashed+n_cc_i,fa_hashed+_nef,cmpFaHashIef());

    n_c__i =  n_cc_i;
    for (int ief = 0; ief < n_cc_i; ++ief) {
      if ( (n_c__i == n_cc_i) && (fa_hashed[ief].ief_compressed >= nef_compressed_group)) {
        n_c__i = ief;
        // could break here but we'll do some checking..
      } else if ( n_c__i != n_cc_i) {
        assert( fa_hashed[ief].ief_compressed >= nef_compressed_group);
      } else {
        assert( fa_hashed[ief].ief_compressed >= 0);
        assert( fa_hashed[ief].ief_compressed < nef_compressed_group);
      }
    }

    n_c__b = _nef;
    for (int ief = n_cc_i; ief < _nef; ++ief) {
      if ( (n_c__b == _nef) && (fa_hashed[ief].ief_compressed >= nef_compressed_group)) {
        n_c__b = ief;
      } else if ( n_c__b != _nef) {
        assert( fa_hashed[ief].ief_compressed >= nef_compressed_group);
      } else {
        assert( fa_hashed[ief].ief_compressed >= 0);
        assert( fa_hashed[ief].ief_compressed < nef_compressed_group);
      }
    }

    n_cc_b = _nef;
    // inside of the four ranges, we now need to sort based on whether the cvofa
    // pair points to a purely extended pair or a compact+extended pair; so we'll
    // sort based on the icf_compressed index ...
    if ( mpi_rank == 0 )
      cout << " ... sorting compact face regions ... " << endl;

    sort(fa_hashed       , fa_hashed+n_c__i, cmpFaHashIcf());
    sort(fa_hashed+n_c__i, fa_hashed+n_cc_i, cmpFaHashIcf());
    sort(fa_hashed+n_cc_i, fa_hashed+n_c__b, cmpFaHashIcf());
    sort(fa_hashed+n_c__b, fa_hashed+n_cc_b, cmpFaHashIcf());

    // recall and icf_compressed index of -1 denotes a purely extended face...
    n_e__i = n_c__i;
    for (int ief = 0; ief < n_c__i; ++ief) {
      if ( (n_e__i == n_c__i) && (fa_hashed[ief].icf_compressed >= 0)) {
        n_e__i = ief;
      } else if ( n_e__i != n_c__i) {
        assert( fa_hashed[ief].icf_compressed >= 0);
      } else {
        assert( fa_hashed[ief].icf_compressed == -1);
      }
    }

    n_ec_i = n_cc_i;
    for (int ief = n_c__i; ief < n_cc_i; ++ief) {
      if ( (n_ec_i == n_cc_i) && (fa_hashed[ief].icf_compressed >= 0)) {
        n_ec_i = ief;
      } else if ( n_ec_i != n_cc_i) {
        assert( fa_hashed[ief].icf_compressed >= 0);
      } else {
        assert( fa_hashed[ief].icf_compressed == -1);
      }
    }

    n_e__b = n_c__b;
    for (int ief = n_cc_i; ief < n_c__b; ++ief) {
      if ( (n_e__b == n_c__b) && (fa_hashed[ief].icf_compressed >= 0)) {
        n_e__b = ief;
      } else if ( n_e__b != n_c__b) {
        assert( fa_hashed[ief].icf_compressed >= 0);
      } else {
        assert( fa_hashed[ief].icf_compressed == -1);
      }
    }

    n_ec_b = n_cc_b;
    for (int ief = n_c__b; ief < n_cc_b; ++ief) {
      if ( (n_ec_b == n_cc_b) && (fa_hashed[ief].icf_compressed >= 0)) {
        n_ec_b = ief;
      } else if ( n_ec_b != n_cc_b) {
        assert( fa_hashed[ief].icf_compressed >= 0);
      } else {
        assert( fa_hashed[ief].icf_compressed == -1);
      }
    }

    // congrats.  we have now sent the 8 ranges; we're going to do some additional
    // sorting to help the cache efficiency now.  resort the 8 ranges based on their
    // ief_compressed index (if the sort was stable, this would not be necessary)
    // this will hopefully order both the icf_compressed and ief_compressed as they
    // are linked through the underlying geometry..

    if ( mpi_rank == 0)
      cout << " ... sorting the cv ief_compressed regions " << endl;
    sort(fa_hashed        ,fa_hashed+n_e__i, cmpFaHashIef());
    sort(fa_hashed+n_e__i ,fa_hashed+n_c__i, cmpFaHashIef());
    sort(fa_hashed+n_c__i ,fa_hashed+n_ec_i, cmpFaHashIef());
    sort(fa_hashed+n_ec_i ,fa_hashed+n_cc_i, cmpFaHashIef());
    sort(fa_hashed+n_cc_i ,fa_hashed+n_e__b, cmpFaHashIef());
    sort(fa_hashed+n_e__b, fa_hashed+n_c__b, cmpFaHashIef());
    sort(fa_hashed+n_c__b, fa_hashed+n_ec_b, cmpFaHashIef());
    sort(fa_hashed+n_ec_b, fa_hashed+n_cc_b, cmpFaHashIef());

    // we'll sort to try to order the cells further by sorting on the min cv
    // of the cvofa pair for each group of ief_compressed...
    if ( mpi_rank == 0)
      cout << " ... final sorting sweep " << endl;

    for (int iter =0; iter < 8; ++iter) {
      int ief_begin, ief_end;
      if ( iter == 0) {
        ief_begin = 0;
        ief_end   = n_e__i;
      } else if ( iter == 1) {
        ief_begin = n_e__i;
        ief_end   = n_c__i;
      } else if ( iter == 2) {
        ief_begin = n_c__i;
        ief_end   = n_ec_i;
      } else if ( iter == 3) {
        ief_begin = n_ec_i;
        ief_end   = n_cc_i;
      } else if ( iter == 4) {
        ief_begin = n_cc_i;
        ief_end   = n_e__b;
      } else if ( iter == 5) {
        ief_begin = n_e__b;
        ief_end   = n_c__b;
      } else if ( iter == 6) {
        ief_begin = n_c__b;
        ief_end   = n_ec_b;
      } else {
        ief_begin = n_ec_b;
        ief_end   = n_cc_b;
      }

      int ief_start =  ief_begin;
      int ief_comp  = -1;

      for (int ief = ief_begin; ief < ief_end; ++ief) {
        if ( ief_comp != fa_hashed[ief].ief_compressed) {
          if ( ief - ief_start > 1)
            sort(fa_hashed+ief_start,fa_hashed+ief,cmpMinCv());
          ief_start = ief;
          ief_comp  = fa_hashed[ief].ief_compressed;
        }
      }

      // check the end..
      if ( (ief_end - ief_start > 1) && (ief_comp == fa_hashed[ief_end-1].ief_compressed))
        sort(fa_hashed+ief_start,fa_hashed+ief_end,cmpMinCv());

      // there is some reference code available elsewhere that will do some additional
      // reordering to generate linked lists of the variety
      // (icv0,icv1):(20,30)-->(30,40)-->(40,50), where these lists can be chained to
      // a length of O(20) using HCP packed grids.  this is cache beneficial in some cases,
      // but hurts in others, so we are omitting it here...
    }

    // check that the gcls are still ok after all of the sorting...
    checkGclExtendedHashed(_nef);
    checkGclCompactHashed(_nef);

    { // report some compression statistics ...
      // check the face comrpession count..
      int count = 0;
      for (int ief = 0; ief < _nef; ++ief) {
        if ( _group_ef[ief] >= 0 )
          ++count;
      }
      assert( count == _nef-nef_compressed+nef_compressed_group);

      double my_cr = double(nef_compressed)/double(_nef);
      assert( (my_cr >= 0.0) && (my_cr <= 1.0+1.0e-10));

      double max_cr, min_cr, avg_cr;
      MPI_Reduce(&my_cr,&max_cr,1,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
      MPI_Reduce(&my_cr,&min_cr,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      MPI_Reduce(&my_cr,&avg_cr,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      if ( mpi_rank == 0) {
        max_cr = 1.0 - max_cr;
        min_cr = 1.0 - min_cr;
        avg_cr = 1.0 - avg_cr/double(mpi_size);
        cout << " > face compression ratio stats [min,avg,max]: "
             << min_cr  << " , " << avg_cr << " ,  " << max_cr << endl;
      }
    }

    {
      int count = 0;
      for (int ifa = 0; ifa < _nfa; ++ifa) {
        if ( _group_fa[ifa] >= 0 )
          ++count;
      }
      assert( count == _nfa-ncf_compressed+ncf_compressed_group);

      double my_cr = double(ncf_compressed)/double(_nfa);
      assert( (my_cr >= 0.0) && (my_cr <= 1.0+1.0e-10));

      double max_cr, min_cr, avg_cr;
      MPI_Reduce(&my_cr,&max_cr,1,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
      MPI_Reduce(&my_cr,&min_cr,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      MPI_Reduce(&my_cr,&avg_cr,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      if ( mpi_rank == 0) {
        max_cr = 1.0 - max_cr;
        min_cr = 1.0 - min_cr;
        avg_cr = 1.0 - avg_cr/double(mpi_size);
        cout << " > compact compression ratio stats [min,avg,max]: "
          << min_cr  << " , " << avg_cr << " ,  " << max_cr << endl;
      }
    }

    // lastly, there are some instances where we need to reference
    // unsorted compact faces with their sorted ef counterparts.  so
    // we need to build a efsorted2fa structure..
    {
      assert( efsorted2fa == NULL); efsorted2fa  = new int[_nef];
      for (int ief = 0; ief < _nef; ++ief)
        efsorted2fa[ief] = -1;

      // build the map of the extended faces to their ief indices
      // through the cvofa.. restrict to the 4 ranges that have
      // compact faces...

      map<pair<int,int>,int> cvoef_sorted_map;
      for (int iter = 0; iter < 4; ++iter) {
        int ief_begin, ief_end;
        if ( iter == 0) {
          ief_begin = n_e__i;
          ief_end   = n_c__i;
        } else if ( iter == 1) {
          ief_begin = n_ec_i;
          ief_end   = n_cc_i;
        } else if ( iter == 2) {
          ief_begin = n_e__b;
          ief_end   = n_c__b;
        } else {
          assert( iter == 3);
          ief_begin = n_ec_b;
          ief_end   = n_cc_b;
        }

        for (int ief = ief_begin; ief < ief_end; ++ief) {
          const int icv0 = fa_hashed[ief].cvofa[0];
          const int icv1 = fa_hashed[ief].cvofa[1];
          pair<int,int> cvoef_pair(icv0,icv1);
          cvoef_sorted_map[cvoef_pair] = ief;
        }
      }

      for (int ifa = 0; ifa < _nfa; ++ifa) {
        const int icv0 = _cvofa[ifa][0]; assert( (icv0 >=0) && (icv0 < ncv));
        const int icv1 = _cvofa[ifa][1]; assert( (icv1 >=0) && (icv1 < ncv_g));

        int ief = -1;
        pair<int,int> cvofa_pair(icv0,icv1);
        map<pair<int,int>,int>::iterator it_ef = cvoef_sorted_map.find(pair<int,int>(icv0,icv1));
        if ( it_ef != cvoef_sorted_map.end()) {
          ief = it_ef->second;
          assert( (ief >= 0)&&(ief < nef));
          assert( icv0 == fa_hashed[ief].cvofa[0]);
          assert( icv1 == fa_hashed[ief].cvofa[1]);
        } else {
          // check the opposite orientation..
          it_ef = cvoef_sorted_map.find(pair<int,int>(icv1,icv0));
          ief = it_ef->second;
          assert( (ief >= 0)&&(ief < nef));
          assert( it_ef != cvoef_sorted_map.end());
          assert( icv1 == fa_hashed[ief].cvofa[0]);
          assert( icv0 == fa_hashed[ief].cvofa[1]);
        }

        efsorted2fa[ief] = ifa;
      }
    }

    // lastly, for some computational efficiency reasons, 0.5 of the c-vector is needed,
    // so we will renormalize here based on that....

    for (int ief_compressed = 0; ief_compressed < nef_compressed; ++ief_compressed) {
      for (int i =0; i < 3; ++i)
        ef_geom[ief_compressed].c[i] *= 0.5;
    }

    checkLinearGrad();
  }

  virtual int advanceSolution() = 0;
  virtual void syncPostState() = 0;
  virtual void initComplete() = 0;
  virtual double calcCfl(double* cfl,const double dt_) const = 0;

  virtual void initialHook() {}
  virtual void initialHookBcs() {}
  virtual void temporalHook() {}
  virtual void finalHook() {}

  void setDt(double& dt,const int step) {

    // user can functionalize the time step so we recheck it every check_interval...
    if ((dt_mode == DT_MODE_UNKNOWN)||(b_dt_fxn&&(step%check_interval == 0))) {
      if (Param * param = getParam("TIMESTEP")) {
        processTimestep(param,false,false); // b_help and b_killfile false
        // processTimestep is robust to errors, so check it has been set to something...
        if (dt_mode == DT_MODE_UNKNOWN) {
          // bring the solver down...
          CERR("TIMESTEP not formatted properly. Please specify one of\n" <<
              "TIMESTEP DT <dt-value>\n " <<
              "TIMESTEP CFL <cfl-target>\n" <<
              "TIMESTEP CFL_DT <cfl-target> <dt-max>");
        }
      }
      else {
        // bring the solver down...
        CERR("TIMESTEP not found. Please specify one of\n" <<
            "TIMESTEP DT <dt-value>\n" <<
            "TIMESTEP CFL <cfl-target>\n" <<
            "TIMESTEP CFL_DT <cfl-target> <dt-max>");
      }
    }

    if (dt_mode == DT_MODE_CONSTANT_DT) {
      dt = dt_data[0];
    }
    else if (dt_mode == DT_MODE_CONSTANT_CFL) {
      const double cfl_target = dt_data[0];
      double my_cfl_max       = calcCfl(NULL,1.0);
      double cfl_max;
      MPI_Allreduce(&my_cfl_max,&cfl_max,1,MPI_DOUBLE,MPI_MAX,mpi_comm);
      dt = cfl_target/cfl_max;
    }
    else if (dt_mode == DT_MODE_CFL_DT) {
      const double cfl_target = dt_data[0];
      const double dt_max     = dt_data[1];
      double my_cfl_max       = calcCfl(NULL,1.0);
      double cfl_max;
      MPI_Allreduce(&my_cfl_max,&cfl_max,1,MPI_DOUBLE,MPI_MAX,mpi_comm);
      dt = min(dt_max,cfl_target/max(cfl_max,1.0E-12));
    }
    else if (dt_mode == DT_MODE_CFL_T) {
      const double cfl_target = dt_data[0];
      const double T_end      = dt_data[1];
      double my_cfl_max       = calcCfl(NULL,1.0);
      double cfl_max;
      MPI_Allreduce(&my_cfl_max,&cfl_max,1,MPI_DOUBLE,MPI_MAX,mpi_comm);
      dt = min(T_end-time,cfl_target/max(cfl_max,1.0E-12));
    }
    else {
      assert(0);
    }

  }

  void processTimestep(Param * param,const bool b_help,const bool b_killfile = true) {
    if (b_help) {
      helpTimestep();
    }
    else {
      // for robust parsing, check that there is atleast a string...
      if (param->size() == 0) {
	WUI(WARN,"Problem parsing TIMESTEP: missing mode");
	helpTimestep();
	return;
      }
      // string should be the mode...
      const string token = MiscUtils::toUpperCase(param->getString());
      if (token == "DT") {
	double tmp;
        bool b_tmp;
	try {
	  tmp = param->getDouble(1);
          b_tmp = !param->isConstDouble(1);
	}
	catch (int e) {
	  WUI(WARN,"Problem parsing TIMESTEP DT. Expecting TIMESTEP DT <dt>");
	  return;
	}
	dt_mode = DT_MODE_CONSTANT_DT;
	dt_data[0] = tmp;
        if (b_killfile) {
          WUI(INFO,"TIMESTEP set to constant dt: " << dt_data[0]);
        }
        else {
          b_dt_fxn = b_tmp;
        }
      }
      else if (token == "CFL") {
	double tmp;
        bool b_tmp;
	try {
	  tmp = param->getDouble(1);
          b_tmp = !param->isConstDouble(1);
	}
	catch (int e) {
	  WUI(WARN,"Problem parsing TIMESTEP CFL. Expecting TIMESTEP CFL <cfl>");
	  return;
	}
	dt_mode = DT_MODE_CONSTANT_CFL;
	dt_data[0] = tmp;
        if (b_killfile) {
          WUI(INFO,"TIMESTEP set to constant cfl: " << dt_data[0]);
        }
        else {
          b_dt_fxn = b_tmp;
        }
      }
      else if (token == "CFL_DT") {
	double tmp[2];
        bool b_tmp;
	try {
	  tmp[0] = param->getDouble(1);
          b_tmp = !param->isConstDouble(1);
	  tmp[1] = param->getDouble(2);
          b_tmp = (b_tmp||(!param->isConstDouble(2)));
	}
	catch (int e) {
	  WUI(WARN,"Problem parsing TIMESTEP CFL_DT. Expecting TIMESTEP CFL_DT <cfl-max> <dt-max>");
	  return;
	}
	dt_mode = DT_MODE_CFL_DT;
	dt_data[0] = tmp[0];
	dt_data[1] = tmp[1];
        if (b_killfile) {
          WUI(INFO,"TIMESTEP set to min of cfl-max: " << dt_data[0] << " and dt-max: " << dt_data[1]);
        }
        else {
          b_dt_fxn = b_tmp;
        }
      }
      else if (token == "CFL_T") {
	double tmp[2];
        bool b_tmp;
	try {
	  tmp[0] = param->getDouble(1);
          b_tmp = !param->isConstDouble(1);
	  tmp[1] = param->getDouble(2);
          b_tmp = (b_tmp||(!param->isConstDouble(2)));
	}
	catch (int e) {
	  WUI(WARN,"Problem parsing TIMESTEP CFL_T. Expecting TIMESTEP CFL_T <cfl-target> <T-end>");
	  return;
	}
	dt_mode = DT_MODE_CFL_T;
	dt_data[0] = tmp[0];
	dt_data[1] = tmp[1];
        if (b_killfile) {
          WUI(INFO,"TIMESTEP set to min of cfl-target: " << dt_data[0] << " and T-end: " << dt_data[1]);
        }
        else {
          b_dt_fxn = b_tmp;
        }
      }
      else {
	WUI(WARN,"Unrecognized TIMESTEP mode: " << token);
      }
    }
  }

  void helpTimestep() {
    WUI(INFO,
        "TIMESTEP sets the mode of timestep advancement. Examples:\n" <<
        "TIMESTEP DT <dt>                   # constant dt mode\n" <<
        "TIMESTEP CFL <cfl-max>             # constant cfl mode\n" <<
        "TIMESTEP CFL_DT <cfl-max> <dt-max> # combination cfl and dt mode\n" <<
	"TIMESTEP CFL_T <cfl-max> <T-end>   # combination cfl and T-end mode\n" <<
	"for more detail see XXX");
  }

  void run() {

    // each flow solver can have a variety of ways that it initializes
    // its data: e.g. INIT_RUP, etc...

    initFromParams();

    // note: important to call initialHook first to allow solver
    // to set the volume data BEFORE setting boundary data...
    // Maybe make this one call so probelms like setting bcs from
    // primitive data where other primitive data has been set are
    // more obvious...

    initialHook();

    // make sure that everything is set/exchanged (this includes
    // derived data and ghost data).  this needs to be called
    // before the boundary condition initialHookBcs is called in
    // case they need data from derived quantities that were not
    // present in the restart file

    initComplete();

    // the boundary conditions may need to pull data from the
    // interior of the domain -- or be otherwise set.  the bcs
    // are local to the processor so there isnt an explicit call
    // to sync in parallel

    initialHookBcs();

    // we may have changed the solution in the previous routines, so
    // clear any evaluations made to this point. Note that this should be done
    // EVERY TIME AFTER the solution has changed to avoid the risk of using an
    // incorrect expression evaluation. In the step loop below, it is done inside
    // advanceSolution().

    CtiRegister::clearCurrentData();

    // solver should be ready to go now..

    nsteps        = getIntParam("NSTEPS",-1);
    runtime_hours = getDoubleParam("RUNTIME",-1.0);
    simtime_stop  = getDoubleParam("SIMTIME",-1.0);

    COUT1("=====================================================================\n" <<
          "running with CHECK_INTERVAL=" << check_interval << "...\n" <<
          "=====================================================================");

    // tell the central log a solver is running...
    if (mpi_rank==0) logger->insertSolverAction(SIMULATION_RUN);

    // on the first step, let the user inspect the initial condition...
    if (step == 0) processStep(step,time,dt,step%check_interval==0);

    int done = 0;
    while (done == 0) {

      timer.start("init");

      // NOTE: time is advanced in advanceSolution below to support time-dependent evaluations
      ++step;
      setDt(dt,step);

      if ((step%check_interval == 0)&&(mpi_rank == 0)) {
	cout <<
          "\n----------------------------------------------------------\n" <<
          " starting step: " << step << " time: " << time+dt << " dt: " << dt <<
          "\n----------------------------------------------------------" << endl;
      }

      timer.split("calc dt");

      addNewLp();
      timer.split("inject Lp");

      buildLpocv();
      timer.split("build Lpocv");

      setLp0();
      timer.split("set Lp0");

      const int ierr = advanceSolution(); // time advanced in here with solution to support time-dependent evaluations
      if (ierr != 0)
        continue;

      temporalHook();

      timer.split("temporal hook");

      processLp();

      timer.split("process Lp");

      // process step should normally return 0, but may return -1 or -2, and
      // the solver should take appropriate action:
      //  0 : no action, not done
      // -1 : soft "stop", i.e. stop and write final result, flush any buffers, etc...
      // -2 : hard "stop!", i.e. stop IMMEDIATELY and release resources

      done = processStep(step,time,dt,step%check_interval==0); // verbose

      if (cht) {
        // ----------------------------------------------
        // for WRITE_CHT we only parse on the first time,
        // as if we were reading the params in interpreted order
        // like surfer and stitch. On subsequet times, use the write_cht_interval,
        // if set, to determine if we write...
        // ----------------------------------------------
        if (cht->write_cht_interval == -1) {
          FOR_PARAM_MATCHING("WRITE_CHT") {
            cht->processWriteCht(&*param,step,false); // false: NOT killfile
          }
          // if still -1 at this point, then set to zero meaning no cht...
          if (cht->write_cht_interval == -1) cht->write_cht_interval = 0;
        }
        // no "else if" here because we may need to process the result...
        if ((cht->write_cht_interval > 0)&&(step%cht->write_cht_interval == 0))
          cht->writeCht(step);

        FOR_PARAM_MATCHING("WRITE_CHT_TECPLOT") {
          
          if (param->size() != 1) {
 
            CWARN(" > parsing error ... expecting: WRITE_CHT_TECPLOT <interval>");

          } else {

            int interval = param->getInt(0);
            if (interval <= 0) { 

              CWARN(" > interval must be positive, not " << interval << " ... skipping");

            } else { 

              if (step%interval == 0) cht->writeTecplotByZone(step); 

            }
          }
        }
      }

      timer.split("process step");

      if (done == -2)
	return;

      // report some timing information about the code...
      if (step%check_interval == 0) {

        if (cht) cht->report(step,time);

        if (mpi_rank == 0) {
          double wtime_prev = wtime;
          wtime = MPI_Wtime();
          wtime_per_interval = wtime-wtime_prev;
          const int8 ncv_global = getNcvGlobal();
          cout << " > time since last check: " << wtime_per_interval << " [s] normalized speed: "
               << (wtime_per_interval)/double(check_interval)*1.0E+6/double(ncv_global)*double(mpi_size)
               << " [core-s/Mcv/step] cores: " << mpi_size << " cvs: " << ncv_global
               << " cvs/core: " << double(ncv_global)/double(mpi_size) << endl;
        }

      }

      if ((done == 0)&&(nsteps >= 0)) {
        if (step >= nsteps) {
          done = -1;
          if (mpi_rank == 0) cout << " > reached NSTEPS " << nsteps << ". Stopping" << endl;
        }
        else if (step%check_interval == 0) {
          // Stats resetting reporting if relevant
          if ( reset_stats_time != HUGE_VAL ) {
            if (mpi_rank == 0) cout << " > STATS will be reset in " << reset_stats_time-time << "s" << endl;
          }
          else if (reset_stats_int > 0) {
            if (mpi_rank == 0) cout << " > STATS will be reset in " << reset_stats_int - (step%reset_stats_int) << " steps" << endl;
          }

          if (mpi_rank == 0) cout << " > NSTEPS set to " << nsteps << ". Stopping in " << nsteps-step << " more steps" << endl;
        }
      }

      if ((done == 0)&&(runtime_hours >= 0.0)) {
        if (mpi_rank == 0) {
          const double wtime = MPI_Wtime();
          if (wtime-mpi_wtime_0  > runtime_hours*3600.0) {
            done = -1;
            cout << " > reached RUNTIME " << runtime_hours << ". Stopping" << endl;
          }
          else if (step%check_interval == 0) {
            cout << " > RUNTIME set to " << runtime_hours << " [hrs]. Wall-clock time so far: " <<
              (wtime-mpi_wtime_0)/3600.0 << " [hrs]. Time left: " << runtime_hours-(wtime-mpi_wtime_0)/3600.0 << " [hrs]" << endl;
          }
        }
        MPI_Bcast(&done,1,MPI_INT,0,mpi_comm);
      }

      if ((done == 0)&&(simtime_stop >= 0.0)) {
        if (mpi_rank == 0) {

          if ( time >= simtime_stop) {
            done = -1;
            cout << " > reached SIMTIME " << simtime_stop << ". Stopping" << endl;
          }
          else if ( step%check_interval == 0) {
            cout << " > SIMTIME set to " << simtime_stop << " .  Remaining time: " << simtime_stop - time << endl;
          }

        }
        MPI_Bcast(&done,1,MPI_INT,0,mpi_comm);
      }

      timer.accumulate();

    }

    flushProbes();
    flushImages();
    flushFwhSurfaces();
    finalHook();

    // for the final sles, use the NEW WRITE_SLES behavior if write_sles_interval has
    // been set...

    if (write_sles_interval > 0) {
      writeSles(step);
    }
    else {
      writeResult();
    }
    if (cht) cht->writeCht(step);

    timer.report();
  }

  void runPost() {

    // for the solvers running inside of post mode- solutions are not
    // advanced and instead we either operate on the result file pair
    // that was loaded or process a series of snapshots

    initialHook();

    // tell the central log a solver is running...
    if (mpi_rank==0) logger->insertSolverAction(SIMULATION_RUN);


    if ( Param* param = getParam("SNAPSHOT")) {

      if ( (param->getString(0) != "NAME") || (param->getString(2) != "RANGE")) {
        CERR( " Invalid snapshot syntax; SNAPSHOT NAME <prefix> RANGE <start> <inc> <last>");
      }

      const string prefix         = param->getString(1);
      const int snapshot_first    = param->getInt(3);
      const int snapshot_inc      = param->getInt(4);
      const int snapshot_last     = param->getInt(5);

      // if we are in snapshot mode, we also need to check for the intended behavior
      // for statistical averaging behavior if stats are already present in the snapshot
      // files.  we can either take ensemble avg data over the snapshots (default) or
      // we can allow the user to step through their averages.

      const bool b_overwrite_stats = getBoolParam("OVERWRITE_STATS",true);

      // reset the stats so that we can take ensemble statistics from the available data

      double snapshot_stats_wgt = 1.0;   // ensemble weighting

      if ( b_overwrite_stats) {

        // we need to clear the read bits for any stats data and reset the stats ..
        // so our own ensemble averaging data is not clobbered by the snapshot files

        CtiRegister::resetStats();
        CtiRegister::clearStatsBits(READ_DATA);

      } else {

        snapshot_stats_wgt = 0.0;  // to prevent updating the stats info..

      }

      // loop over the snapshots and give access to the step processing ..

      int done = 0;

      for (int snapshot = snapshot_first; snapshot <= snapshot_last; snapshot += snapshot_inc) {

        CtiRegister::clearCurrentData();

        char filename[128];
        MiscUtils::buildIndexedFilename(filename,prefix.c_str(),snapshot,"sles");
        RestartHashUtilities::clearSlesHash();
        readData(filename);

        // we should have a valid step and time at this stage from the snapshot file
        // we will report to stdout just to confirm

        if ( mpi_rank == 0)
          cout << " > processing snapshot: " << filename << " ; step = "
               << step << " , time = " << time << endl;

        // solver specific sync of the data and population of primitive data fields..

        syncPostState();

        temporalHook();

        // interactive viewing of snapshots

        if ((snapshot > snapshot_first)&&(checkParam("INTERACTIVE")))
          kfr.setHold(true);

        // process step is called to activate the probing, imaging, etc.

        done = processStep(step, time,snapshot_stats_wgt, false);

        if ( done == -2)
          return;
        else if ( done == -1)
          break;

      }

      flushProbes();
      flushImages();
      flushFwhSurfaces();
      finalHook();

      if ( (done >= -1) && (checkParam("RESULT")) ) {

        string result_name = getStringParam("RESULT");
        MiscUtils::mkdir_for_file(result_name.c_str());

        if ( mpi_rank == 0)
          cout << " Writing solution file \"" << result_name << "\"..." << endl;
        writeData(result_name);

      }

    } else {

      // run the post operations on a single piece of data that has been read in...

      CtiRegister::clearCurrentData();

      syncPostState();

      temporalHook();

      // set the dt to 0 in process step so that the stats are not overwritten ..

      processStep(step,time,0.0,false);

      flushProbes();
      flushImages();
      flushFwhSurfaces();
      finalHook();

      if ( checkParam("RESULT")) {

        string result_name = getStringParam("RESULT");
        MiscUtils::mkdir_for_file(result_name.c_str());

        if ( mpi_rank == 0)
          cout << " Writing solution file \"" << result_name << "\"..." << endl;
        writeData(result_name);
      }

    }

  }

  void calcCvGrad(double (*__restrict__ dudx)[3][3], double (*__restrict__ dpdx)[3],
                  const double (*__restrict__ u)[3], const double *__restrict__ p) {

    for (int icv = 0; icv < ncv; ++icv) {
      const int coc_f = cvocv_i[icv];
      for (int i = 0; i < 3; ++i)
        dpdx[icv][i] = cvocv_grad_coeff[coc_f][i] * p[icv];

      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          dudx[icv][i][j] = cvocv_grad_coeff[coc_f][i]* u[icv][j];

      for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        for (int i =0; i < 3; ++i)
          dpdx[icv][i] += cvocv_grad_coeff[coc][i] * p[icv_nbr];

        for (int i =0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            dudx[icv][i][j] += cvocv_grad_coeff[coc][i] * u[icv_nbr][j];
      }//for_coc

      // the cv grad coefficient already has a inv_vol in it...
    }
  }

  void checkCvGradPacked() {

    double (*v)[3] = new double[ncv_g2][3];
    double * phi   = new double[ncv_g2];

    const double b0   = 0.88;
    const double b[3] = {1.1,2.1,3.4};
    for (int icv = 0; icv < ncv; ++icv) {
      phi[icv] = b0 + DOT_PRODUCT(x_cv[icv],b);
      for (int i =0; i < 3; ++i)
        v[icv][i] = phi[icv];
    }

    updateCv2Data(v);
    updateCv2Data(phi);

    double (*dvdx)[3][3] = new double[ncv][3][3];
    double (*dphidx)[3]  = new double[ncv][3];

    calcCvGrad(dvdx,dphidx,v,phi);

    for (int icv = 0; icv < ncv ; ++icv) {
      for (int i =0; i < 3; ++i) {
        dphidx[icv][i] /= b[i];
        dphidx[icv][i] -= 1.0;
      }

      for (int i =0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          dvdx[icv][i][j] /= b[i];
          dvdx[icv][i][j] -= 1.0;
        }
      }
    }

    MiscUtils::dumpRange(dphidx,ncv, "dphidx check ( should be zero)");
    MiscUtils::dumpRange(dvdx  ,ncv, "dvdx check (should be zero)");

    delete[] dvdx;
    delete[] dphidx;
    delete[] v;
    delete[] phi;

  }

  // flow solver implementing specific processParam...

  virtual void processParam(Param * param, const int step, const double time) {

    bool b_help = false;
    string token = MiscUtils::toUpperCase(param->getName());
    if (token == "HELP") {
      b_help = true;
      if (param->size() == 0) {
        helpHelp(true); // true == write header
        return;
      }
      token = MiscUtils::toUpperCase(param->getString(0));  // first token is keyword for help
    }

    bool found = false;
    if (token == "RUNTIME") {
      processRuntime(param,b_help);
      found = true;
    }
    else if (token == "NSTEPS") {
      processNsteps(param,b_help);
      found = true;
    }
    else if (token == "TIMESTEP") {
      processTimestep(param,b_help);
      found = true;
    }
    else if (token == "CHECK_INTERVAL") {
      processCheckInterval(param,b_help);
      found = true;
    }
    else if (token == "WRITE_CHT") {
      if (cht) {
        cht->processWriteCht(param,step,true);
        found = true;
      }
    }

    // we did not recognize the param, so call StaticSolver's implementation...
    if (!found) StaticSolver::processParam(param,step,time);

  }

  virtual void helpHelp(const bool first) {

    if (first) {
      WUI(WARN,"HELP should be followed by the command you want help with: e.g. HELP WRITE_IMAGE\n" <<
	  "HELP is available for:");
    }
    WUI(WARN,"  RUNTIME NSTEPS TIMESTEP CHECK_INTERVAL");

    return StaticSolver::helpHelp(false); // false avoids repeating help header

  }

  void processRuntime(Param * param,const bool b_help) {

    /*
    CTI_MESSAGE(CTI_INFO,"blah");
    CTI_MESSAGE(CTI_INFO_AND_IMAGE);
    _IMAGE,"size scaled successfully");
    CTI_MESSAGE(CTI_WARN,"blah");
    */

    //double INFO = 2.1;

    if (b_help) {
      WUI(INFO,"RUNTIME <hours> sets the total wall-clock time at which the solver will stop, writing a final sles file");
    }
    else {
      double new_runtime_hours;
      try {
	new_runtime_hours = param->getDouble();
      }
      catch (int e) {
	WUI(WARN,"problem parsing RUNTIME token. RUNTIME <double> expects time in hours. Leaving current RUNTIME: " << runtime_hours);
	return;
      }
      runtime_hours = new_runtime_hours;
      WUI(INFO,"RUNTIME set to " << runtime_hours << " hours");
    }

  }

  void processNsteps(Param * param,const bool b_help) {

    if (b_help) {
      WUI(INFO,"NSTEPS <int> sets the final step index at which the solver will stop, writing a final sles file");
    }
    else {
      int new_nsteps;
      try {
	new_nsteps = param->getInt();
      }
      catch (int e) {
	WUI(WARN,"problem parsing NSTEPS token. NSTEPS <int> expects an int. Leaving current NSTEPS: " << nsteps);
	return;
      }
      nsteps = new_nsteps;
      WUI(INFO,"NSTEPS set to " << nsteps);
    }

  }

  void processCheckInterval(Param * param,const bool b_help) {

    if (b_help) {
      WUI(INFO,"CHECK_INTERVAL <int> sets the interval between steps that the solver reports various checks to output.");
    }
    else {
      int new_check_interval;
      try {
	new_check_interval = param->getInt();
      }
      catch (int e) {
	WUI(WARN,"problem parsing CHECK_INTERVAL. CHECK_INTERVAL <int> expects an int. Leaving current CHECK_INTERVAL: " << check_interval);
	return;
      }
      check_interval = new_check_interval;
      WUI(INFO,"CHECK_INTERVAL set to " << check_interval);
    }

  }

  // force data requests to the flow solver ...

  virtual void processSolverSpecificDiagnostics(const int _step, const double _time, const bool _verbose) {

    doForces(_step,_time);

    doBodyForces(_step,_time);
  }

  void reportTiming(Param * param) {
    //write json file with timing information
    // only rank 0 write the json...
    if (mpi_rank == 0) {

      string dt_mode="";
      if (Param * dtModeParam = getParam("TIMESTEP")) {
        dt_mode = dtModeParam->getString();
      }

      string prefix = "timingReport";
      if (param->size()>1 && param->getString(0) == "NAME"){
        prefix = param->getString(1);
      }

      double wtime_elapsed = (MPI_Wtime()-mpi_wtime_0)/3600.0; //hours

      const string filename = prefix+".json";
      MiscUtils::mkdir_for_file(filename);
      const string tmp_filename = MiscUtils::makeTmpPrefix(filename);
      FILE * fp = fopen(tmp_filename.c_str(),"w");
      assert(fp);

      fprintf(fp,"{\n\"dt_mode\":\"%s\"",dt_mode.c_str());
      fprintf(fp,",\n\"check_interval\":%d",check_interval);
      fprintf(fp,",\n\"wtime_per_interval\":%g",wtime_per_interval);
      fprintf(fp,",\n\"stime_per_interval\":%g",dt*check_interval); //Approximate if CFL mode
      fprintf(fp,",\n\"stime_elapsed\":%g",time);
      fprintf(fp,",\n\"wtime_elapsed\":%g",wtime_elapsed);
      fprintf(fp,",\n\"wtime_max\":%g",runtime_hours);
      fprintf(fp,",\n\"steps_elapsed\":%d",step);
      fprintf(fp,",\n\"steps_max\":%d",nsteps);
      fprintf(fp,"\n}\n");
      fclose(fp);

      remove(filename.c_str());
      rename(tmp_filename.c_str(),filename.c_str());

    }

  }

  virtual void computeForces()  = 0;

  void doForces(const int _step = 0,const double _time = 0.0) {

    //populate forceZoneMap with a list of zones on which to compute the force at this step

    static int force_write_rank = 0;

    assert(forceZoneMap.empty()&&!forceZoneBuf);

    FOR_PARAM_MATCHING("WRITE_FORCE") {

      if ( fioMap.find(param->str()) == fioMap.end()) {

        ForceIO force_request;

        int iarg = 0;
        while ( iarg < param->size()) {

          string token = param->getString(iarg++);

          if ( token == "INTERVAL") {

            force_request.interval = param->getInt(iarg++);
            if (force_request.interval <= 0) {
              CWARN(" > FORCE INTERVAL expects a positive integer");
              // invalid value entered, so treat as unset (default value)
              force_request.interval = -1;
            }
          } else if ( (token == "BFZONE") || (token == "FAZONE")) {

            string zone_name = param->getString(iarg++);
            MiscUtils::splitCsv(force_request.bf_zones,zone_name);

          }

          else if ( token == "NAME" ) {

            force_request.name = param->getString(iarg++);

          } else if ( token == "MOMENT") {

            for (int i = 0; i < 3; ++i)
              force_request.moment_centroids.push_back(param->getDouble(iarg++));

          } else {

            CERR( " > unrecognized token in WRITE_FORCE : " << token);

          }
        }

        // do some checking of the arguments now ..

        string err_string = " > WRITE_FORCE request : " + param->str();
        bool err_set      = false;

        if ( force_request.name == "") {

          err_set = true;
          err_string += "\n  >    missing NAME";

        }

        if ( force_request.interval == -1 ) {

          err_set = true;
          err_string += "\n >     missing valid INTERVAL";

        }

        if ( force_request.bf_zones.empty()) {

          err_set = true;
          err_string += "\n >     missing comma separated BFZONE or FAZONE list";

        } else {

          for (int i = 0, lim = force_request.bf_zones.size(); i < lim; ++i) {

            if ( bfZoneNameMap.find(force_request.bf_zones[i]) == bfZoneNameMap.end()) {

              err_set = true;
              err_string += "\n >   invalid bf zone name : " + force_request.bf_zones[i];

            }
          }
        }


        if ( err_set) {

          CERR( err_string);

        } else {

          // this force request has finally succeeded-- hooray.  set the write rank round robin.

          force_request.write_rank = force_write_rank++;
          if ( force_write_rank == mpi_size)
            force_write_rank = 0;

          fioMap[param->str()] = force_request;


          double my_buf[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
          for (int iz = 0, lim = force_request.bf_zones.size(); iz < lim; ++iz) {
            const string zone_name = force_request.bf_zones[iz];
            BfZone * bf_zone       = getBfZone(zone_name);
            for (int ibf = 0; ibf < bf_zone->nbf; ++ibf) {
              FOR_I3 my_buf[i] += bf_zone->n_bf[ibf][i];
              const double tmp[3] = CROSS_PRODUCT(bf_zone->x_bf[ibf],bf_zone->n_bf[ibf]);
              FOR_I3 my_buf[i+3] += tmp[i];
            }
          }
          double buf[6];
          MPI_Reduce(my_buf,buf,6,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
          assert( force_request.moment_centroids.size()%3 == 0);
          const int n_pts = force_request.moment_centroids.size()/3;
          if (mpi_rank == 0) {
            cout << " > " << force_request.name << " int(n_bf)dA: " << COUT_VEC(buf) << " int((x_bf-x_mom) x n_bf)dA: ";
            for (int im = 0; im < n_pts; ++im) {
              const double xp[3] = {force_request.moment_centroids[im*3+0],
                                    force_request.moment_centroids[im*3+1],
                                    force_request.moment_centroids[im*3+2]};
              const double int_mom_x_n_dA[3] = CROSS_PRODUCT(xp,buf);
              const double tmp[3] = DIFF(buf+3,int_mom_x_n_dA);
              cout << COUT_VEC(tmp) << " ";
            }
            cout << endl;
          }

        }

      }
    }//FOR_PARAM_MATCHING("WRITE_FORCE")


    int index = 0;
    for(map<string,ForceIO>::iterator it = fioMap.begin(); it != fioMap.end(); ++it) {

      if ( _step%it->second.interval == 0 ) {

        // this is a force request that we need to process during this step ...

        for (int i = 0, lim = it->second.bf_zones.size(); i < lim; ++i) {

          if ( forceZoneMap.find(it->second.bf_zones[i]) == forceZoneMap.end())
            forceZoneMap[it->second.bf_zones[i]] = index++;

        }
      }
    }

    // possibly we have no force requests at this time step to process ...

    if ( index == 0)
      return;


    assert( index == int(forceZoneMap.size()));
    assert( forceZoneBuf == NULL);
    assert( momentZoneBuf == NULL);

    forceZoneBuf  = new double[index][9];
    momentZoneBuf = new double[index][9];
    computeForces();

    for ( map<string,ForceIO>::iterator it = fioMap.begin(); it != fioMap.end(); ++it) {

      if ( _step%it->second.interval == 0 )
        processWriteForce(_step,_time,it->second,it->first);

    }

    DELETE(forceZoneBuf);
    DELETE(momentZoneBuf);
    forceZoneMap.clear();
  }

  void processWriteForce(const int step_, const double time_, const ForceIO& force_request, const string& param_str) {

    assert(forceZoneBuf&&momentZoneBuf&&!forceZoneMap.empty());

    //sum zone grouping for this WRITE_FORCE param

    double bf_info[4] = {0.0,0.0,0.0,0.0}; //area, n_bf
    double my_f_buf[9], my_m_buf[9]; //conv_force, p_force, visc_force

    for (int i = 0; i < 9; ++i) {

      my_f_buf[i] = 0.0;
      my_m_buf[i] = 0.0;
    }

    for (int iz = 0, lim = force_request.bf_zones.size(); iz < lim; ++iz) {

      const string zone_name = force_request.bf_zones[iz];
      BfZone * bf_zone       = getBfZone(zone_name);
      bf_info[0]            += bf_zone->area_global;
      FOR_I3 bf_info[1+i]   += bf_zone->n_global[i];

      for (int i = 0; i < 9; ++i) {

        const int idx = forceZoneMap[zone_name];
        my_f_buf[i] += forceZoneBuf[idx][i];
        my_m_buf[i] += momentZoneBuf[idx][i];
      }
    }


    double f_buf[9];
    MPI_Reduce(my_f_buf,f_buf,9,MPI_DOUBLE,MPI_SUM,force_request.write_rank,mpi_comm);

    //do file io for force writing...

    writeForces(force_request,step_,time_,bf_info,f_buf,param_str);

    //do file io for moment writiing

    assert( force_request.moment_centroids.size()%3 == 0);
    const int n_pts = force_request.moment_centroids.size()/3;

    for (int im=0; im < n_pts; ++im) {

      double my_mp_buf[9];
      const double xp[3] = {force_request.moment_centroids[im*3+0],
                            force_request.moment_centroids[im*3+1],
                            force_request.moment_centroids[im*3+2]};

      FOR_I3 {
        const double fb[3] = {my_f_buf[i*3+0],my_f_buf[i*3+1],my_f_buf[i*3+2]};
        const double mp[3] = CROSS_PRODUCT(xp,fb);
        FOR_J3 my_mp_buf[i*3+j] = mp[j];
      }

      // could coalsece the following two reduces into 1 since the data is going to
      // the same rank -- SB..

      double m_buf[9], mp_buf[9];
      MPI_Reduce(my_m_buf,m_buf,9,MPI_DOUBLE,MPI_SUM,force_request.write_rank,mpi_comm);
      MPI_Reduce(my_mp_buf,mp_buf,9,MPI_DOUBLE,MPI_SUM,force_request.write_rank,mpi_comm);

      writeMoments(force_request,step_,time_,bf_info,m_buf,mp_buf,xp,param_str);
    }

  }

  void writeForces(const ForceIO& force_request, const int &step, const double &time, const double bf_info[4], const double buf[9], const string& param_str) {
    int ierr = 0;

    if (mpi_rank == force_request.write_rank) {

      ofstream ofile;
      size_t fwidth = 13;

      string name = force_request.name+".dat";

      if (MiscUtils::fileExists(name)) {
        ofile.open(name.c_str(),ios::app);
        if (!ofile.is_open()) {
         cerr << ERRSTART << "Error: could not open WRITE_FORCE file: " << name << ERREND << endl;
         ierr = 1;
        }
      }
      else {
        MiscUtils::mkdir_for_file(name);
        ofile.open(name.c_str());
        if (!ofile.is_open()) {
         cerr << ERRSTART << "Error: could not open WRITE_FORCE file: " << name << ERREND << endl;
         ierr = 1;
        }
        else {

         // on the first time, write the header...

         // the parameter...
         ofile << "# " << param_str << endl;

         // the column index...
         ofile << "# " <<
           setw(fwidth-2) << 1;
         for (int i = 2; i <= 18; ++i)
           ofile << setw(fwidth) << i;
         ofile << endl;

         // the variable name...
         ofile << "# " <<
           setw(fwidth-2) << "step" <<
           setw(fwidth) << "time" <<
           setw(fwidth) << "area";

         ofile <<
           setw(fwidth) << "area-x" << setw(fwidth) << "area-y" << setw(fwidth) << "area-z";

         ofile <<
           setw(fwidth) << "f-conv-x" << setw(fwidth) << "f-conv-y" << setw(fwidth) << "f-conv-z";

         ofile <<
           setw(fwidth) << "f-p-x" << setw(fwidth) << "f-p-y" << setw(fwidth) << "f-p-z";

         ofile <<
           setw(fwidth) << "f-visc-x" << setw(fwidth) << "f-visc-y" << setw(fwidth) << "f-visc-z";

         ofile <<
           setw(fwidth) << "f-total-x" << setw(fwidth) << "f-total-y" << setw(fwidth) << "f-total-z";

         ofile << endl;

        }
      }
      if (ofile.is_open()) {

        // parse the buf data...
        double area                            = bf_info[0];
        double area_vec[3]; FOR_I3 area_vec[i] = bf_info[1+i];
        double f_conv[3]; FOR_I3 f_conv[i] = buf[i];
        double f_p[3];    FOR_I3 f_p[i]    = buf[3+i];
        double f_visc[3]; FOR_I3 f_visc[i] = buf[3+3+i];

        ofile <<
         setw(fwidth) << step <<
         setw(fwidth) << time <<
         setw(fwidth) << area;

        FOR_I3 ofile <<
         setw(fwidth) << area_vec[i];

        FOR_I3 ofile <<
         setw(fwidth) << f_conv[i];

        FOR_I3 ofile <<
         setw(fwidth) << f_p[i];

        FOR_I3 ofile <<
         setw(fwidth) << f_visc[i];

        FOR_I3 ofile <<
         setw(fwidth) << f_conv[i] + f_p[i] + f_visc[i];

        ofile << endl;

        ofile.close();
      }
    }

    MPI_Bcast(&ierr,1,MPI_INT,force_request.write_rank,mpi_comm);

    if ( ierr != 0) {

      CERR( " > unable open the file for the WRITE_FORCE request: " << param_str);

    }

  }

  void writeMoments(const ForceIO& force_request, const int &step, const double &time, const double bf_info[4],const double buf[9], const double p_buf[9],const double xp[3], const string& param_str) {

    int ierr = 0;

    if (mpi_rank == force_request.write_rank) {

      ofstream ofile;
      size_t fwidth = 13;

      stringstream ss;
      ss << force_request.name << "_moment_x" << xp[0] << "y" << xp[1] << "z" << xp[2] << ".dat";
      string name = ss.str();

      if (MiscUtils::fileExists(name)) {
        ofile.open(name.c_str(),ios::app);
        if (!ofile.is_open()) {
         cerr << ERRSTART << "Error: could not open WRITE_FORCE file: " << name << ERREND << endl;
         ierr = 1;
        }
      }
      else {
        MiscUtils::mkdir_for_file(name);
        ofile.open(name.c_str());
        if (!ofile.is_open()) {
         cerr << ERRSTART << "Error: could not open WRITE_FORCE file: " << name << ERREND << endl;
         ierr = 1;
        }
        else {

         // on the first time, write the header...

         // the parameter...
         ofile << "# " << param_str << endl;

         // on the second line write the point about which the moment is reported...
         ofile << "# xp " << xp[0] << " " << xp[1] << " " << xp[2] << endl;

         // the column index...
         ofile << "# " <<
           setw(fwidth-2) << 1;
         for (int i = 2; i <= 18; ++i)
           ofile << setw(fwidth) << i;
         ofile << endl;

         // the variable name...
         ofile << "# " <<
           setw(fwidth-2) << "step" <<
           setw(fwidth) << "time" <<
           setw(fwidth) << "area";

         ofile <<
           setw(fwidth) << "area-x" << setw(fwidth) << "area-y" << setw(fwidth) << "area-z";

         ofile <<
           setw(fwidth) << "m-conv-x" << setw(fwidth) << "m-conv-y" << setw(fwidth) << "m-conv-z";

         ofile <<
           setw(fwidth) << "m-p-x" << setw(fwidth) << "m-p-y" << setw(fwidth) << "m-p-z";

         ofile <<
           setw(fwidth) << "m-visc-x" << setw(fwidth) << "m-visc-y" << setw(fwidth) << "m-visc-z";

         ofile <<
           setw(fwidth) << "m-total-x" << setw(fwidth) << "m-total-y" << setw(fwidth) << "m-total-z";

         ofile << endl;

        }
      }
      if (ofile.is_open()) {

        // parse the buf data...
        double area                            = bf_info[0];
        double area_vec[3]; FOR_I3 area_vec[i] = bf_info[1+i];
        double mom_conv[3]; FOR_I3 mom_conv[i] = buf[i] - p_buf[i];
        double mom_p[3];    FOR_I3 mom_p[i]    = buf[3+i] - p_buf[3+i];
        double mom_visc[3]; FOR_I3 mom_visc[i] = buf[3+3+i] - p_buf[3+3+i];

        ofile <<
         setw(fwidth) << step <<
         setw(fwidth) << time <<
         setw(fwidth) << area;

        FOR_I3 ofile <<
         setw(fwidth) << area_vec[i];

        FOR_I3 ofile <<
         setw(fwidth) << mom_conv[i];

        FOR_I3 ofile <<
         setw(fwidth) << mom_p[i];

        FOR_I3 ofile <<
         setw(fwidth) << mom_visc[i];

        FOR_I3 ofile <<
         setw(fwidth) << mom_conv[i] + mom_p[i] + mom_visc[i];

        ofile << endl;

        ofile.close();
      }
    }

    MPI_Bcast(&ierr,1,MPI_INT,force_request.write_rank,mpi_comm);

    if ( ierr != 0 ) {

      CERR( " > unable to open file for WRITE_FORCE : " << param_str << "\n with moment xp : "
             << COUT_VEC(xp) );

    }
  }

  // --------------------------
  // porous media models...
  // --------------------------

  //
  //
  virtual void computeBodyForces(){
    // Define at the solver level where body forces should be computed in a given
    // time step. This function must first call parseBodyForces()...
    //
    // parseBodyForces()
    //
    // ...solver computation of body forces...
    //
    //
  }

  void parseBodyForces() {
    //populate forceZoneMap with a list of zones on which to compute the force at this step

    static int body_force_write_rank = 0;

    FOR_PARAM_MATCHING("ADD_BODY_FORCE") {

      if ( bodyForceNameMap.find(param->str()) == bodyForceNameMap.end()) {

        double xr[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}}; // principal axis (default xyz)
        double Cr[3];    // inertial coefficient in relative reference frame
        double Dr[3];    // viscous coefficient in relative reference frame
        bool b_dir1 = false, b_dir2 = false, b_Cr = false, b_Dr = false; // some booleans for checking
        bool b_axis_rot = false, b_x_rot = false, b_omega_rot = false;   // some booleans for cheacking
        string sbin_filename = "";
        int geom = -1;
        SimpleGeom* simple_geom = NULL;

        BodyForce body_force;

        int iarg = 0;
        while ( iarg < param->size()) {

          string token = param->getString(iarg++);

          if ( token == "INTERVAL") {

            body_force.interval = param->getInt(iarg++);
            if (body_force.interval <= 0) {
              CWARN(" > ADD_BODY_FORCE INTERVAL expects a positive integer");
              // invalid value entered, so treat as unset (default value)
              body_force.interval = -1;
            }

          }

          else if ( token == "NAME" ) {

            body_force.name = param->getString(iarg++);

          }

          else if ( token == "MOMENT") {

            for (int i = 0; i < 3; ++i)
              body_force.moment_centroids.push_back(param->getDouble(iarg++));

          }

          // TODO allow for simple geoms...
          else if ( token == "GEOM") {

            token = param->getString(iarg++);

            if ( token == "SBIN") {
              geom = FILE_GEOM;
              sbin_filename = param->getString(iarg++);
            }
            else {
              --iarg;
              geom = SIMPLE_GEOM;
              simple_geom = newSimpleGeom(&(*param),iarg);
            }

          }

          else if ( token == "TYPE" ) {

            token = param->getString(iarg++);
            if (token == "PMM")
              body_force.type = PMM_TYPE;
            else if (token == "MRF")
              body_force.type = MRF_TYPE;
            else if (token == "HOOK")
              body_force.type = HOOK_TYPE;

          }

          // PMM specific args...

          else if ( token == "DIR_1") {

            b_dir1 = true;
            FOR_I3 xr[0][i] = param->getDouble(iarg++);

          }

          else if ( token == "DIR_2") {

            b_dir2 = true;
            FOR_I3 xr[1][i] = param->getDouble(iarg++);

          }

          else if ( token == "INERTIAL_COEFF") {

            b_Cr = true;
            FOR_I3 Cr[i] = param->getDouble(iarg++);

          }

          else if ( token == "VISCOUS_COEFF") {

            b_Dr = true;
            FOR_I3 Dr[i] = param->getDouble(iarg++);

          }

          // MRF specific args...

          else if ( token == "AXIS_ROT" ) {

            b_axis_rot = true;
            FOR_I3 body_force.axis_rot[i] = param->getDouble(iarg++);
            double axis_mag = MAG(body_force.axis_rot);
            FOR_I3 body_force.axis_rot[i] /= axis_mag; //normalize

          }

          else if ( token == "X_ROT" ) {

            b_x_rot = true;
            FOR_I3 body_force.x_rot[i] = param->getDouble(iarg++);

          }

          else if ( token == "OMEGA_ROT" ) {

            b_omega_rot = true;
            //assumes rev/s, convert to rad/s
            body_force.omega_rot = M_PI*2.0*param->getDouble(iarg++);

          }

          else {

            CERR( " > unrecognized token in ADD_BODY_FORCE : " << token);

          }
        }

        // do some checking of the arguments now ..

        string err_string = " > ADD_BODY_FORCE request : " + param->str();
        bool err_set      = false;

        if ( body_force.name == "") {

          err_set = true;
          err_string += "\n  >    missing NAME";

        }

        if ( geom == -1 ) {

          err_set = true;
          err_string += "\n >     missing GEOM";

        }
        else if (geom == FILE_GEOM && sbin_filename == "") {

          err_set = true;
          err_string += "\n >     missing SBIN";

        }

        if ( body_force.type == PMM_TYPE) {

          if (b_dir1 && !b_dir2) {

            err_set = true;
            err_string += "\n >     missing DIR_2 to match with DIR_1 to form axis";

          }

          if (!b_dir1 && b_dir2) {

            err_set = true;
            err_string += "\n >     missing DIR_1 to match with DIR_2 to form axis";

          }

          if (!b_Cr) {

            err_set = true;
            err_string += "\n >     missing INERTIAL_COEFF";

          }

          if (!b_Dr) {

            err_set = true;
            err_string += "\n >     missing VISCOUS_COEFF";

          }
        }
        else if ( body_force.type == MRF_TYPE) {

          if (!b_axis_rot) {

            err_set = true;
            err_string += "\n >     missing AXIS_ROT";

          }

          if (!b_x_rot) {

            err_set = true;
            err_string += "\n >     missing X_ROT";

          }

          if (!b_omega_rot) {

            err_set = true;
            err_string += "\n >     missing OMEGA_ROT (revolutions/s)";

          }

        }
        else if ( body_force.type == HOOK_TYPE ){
          //no error checking here
        }
        else {

          err_set = true;
          err_string += "\n  >    missing/invalid TYPE";

        }

        if ( err_set) {

          CERR( err_string);

        }
        else {

          // this body_force request has finally succeeded-- hooray.
          // set the write rank round robin and get coeffs in solver reference frame (if need be)...

          body_force.write_rank = body_force_write_rank++;
          if ( body_force_write_rank == mpi_size)
            body_force_write_rank = 0;

          if (body_force.type == PMM_TYPE) {
            if (b_dir1) {
              assert(b_dir2); // should of been caught above
              xr[2][0] = CROSS_PRODUCT_0(xr[0],xr[1]);
              xr[2][1] = CROSS_PRODUCT_1(xr[0],xr[1]);
              xr[2][2] = CROSS_PRODUCT_2(xr[0],xr[1]);
              FOR_I3 {
                const double mag = MAG(xr[i]);
                assert(mag > 0.0); // if you hit this we need add some more chacking
                FOR_J3 xr[i][j] /= mag;
              }
            }

            // now transform into absolute basis

            const double x[3][3] = { {1.0, 0.0, 0.0},
                                     {0.0, 1.0, 0.0},
                                     {0.0, 0.0, 1.0} }; // absolute basis

            double R[3][3]; // transformation matrix...
            FOR_I3 FOR_J3 R[i][j] = DOT_PRODUCT(x[i],xr[j]);
            FOR_I3 FOR_J3 FOR_M3 FOR_N3 {
              if (m == n) body_force.C[i][j] += R[i][m]*R[j][n]*Cr[m];
            }

            FOR_I3 FOR_J3 FOR_M3 FOR_N3 {
              if (m == n) body_force.D[i][j] += R[i][m]*R[j][n]*Dr[m];
            }
          }

          // add body_force to bodyForceVec and add its index to the bodyForceNameMap
          bodyForceNameMap[param->str()] = bodyForceVec.size();

          // flag cv's with inside of surface...
          flagBodyForce(bodyForceVec.size(),geom,simple_geom,sbin_filename);

          bodyForceVec.push_back(body_force);

        }

      }

    }//FOR_PARAM_MATCHING("ADD_BODY_FORCE")

    for (int ii = 0, limit = bodyForceVec.size(); ii < limit; ++ii) bodyForceVec[ii].zero();
  }

  void doBodyForces(const int _step = 0,const double _time = 0.0) {
    // see if we need to write any body...

    for(map<const string,int>::iterator it = bodyForceNameMap.begin(); it != bodyForceNameMap.end(); ++it) {

      if ((bodyForceVec[it->second].interval != -1)&&(_step%bodyForceVec[it->second].interval == 0)) {

        // this is a body_force force request that we need to process during this step...

        processWriteBodyForce(_step,_time,bodyForceVec[it->second],it->first);

      }

    }

  }

  void flagBodyForce(const int index,const int geom,SimpleGeom* simple_geom,const string& sbin_filename) {

    if (body_force_flag == NULL) {
      body_force_flag = new int[ncv_g];
      FOR_ICV_G body_force_flag[icv] = -1;
    }

    double my_buf[2] = { 0.0, 0.0 };
    assert(vol_cv);

    if (geom == FILE_GEOM) {

      SurfaceShm surface;
      surface.readBinary(sbin_filename);
      surface.calcVolumeAndCentroid();

      // build the bbox of the points to be queried...
      double bbmin[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
      double bbmax[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
      FOR_ICV {
        // skip ones already flagged (first wins)
        if (body_force_flag[icv] == -1) {
          FOR_I3 {
            bbmin[i] = min(bbmin[i],x_cv[icv][i]);
            bbmax[i] = max(bbmax[i],x_cv[icv][i]);
          }
        }
      }

      surface.initPointIsInside(bbmin,bbmax);

      FOR_ICV {
        if (body_force_flag[icv] == -1) {
          if (surface.pointIsInside(x_cv[icv])) {
            body_force_flag[icv] = index;
            // sum count and vol...
            ++my_buf[0];
            my_buf[1] += vol_cv[icv];
          }
        }
      }

    }
    else if (geom == SIMPLE_GEOM) {
      assert(simple_geom);

      FOR_ICV {
        if (body_force_flag[icv] == -1) {
          if (simple_geom->pointIsInside(x_cv[icv])) {
            body_force_flag[icv] = index;
            // sum count and vol...
            ++my_buf[0];
            my_buf[1] += vol_cv[icv];
          }
        }
      }
      delete simple_geom; simple_geom = NULL;

    }

    updateCvData(body_force_flag);

    double buf[2];
    MPI_Reduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    if (mpi_rank == 0)
      cout << " > flagBodyForce with index: " << index << " cv count: " << int8(buf[0]) << " volume: " << buf[1] << endl;

    // write the points...
    /*
      {
      int count = int(my_buf[0]);
      double (*xp)[3] = new double[count][3];
      int ii = 0;
      FOR_ICV {
      if (body_force_flag[icv] == index) {
      FOR_I3 xp[ii][i] = x_cv[icv][i];
      ++ii;
      }
      }
      assert(ii == count);
      char filename[128];
      sprintf(filename,"pts.%04d.dat",index);
      GeomUtils::writePtsTecplot(filename,xp,count);
      }
      MPI_Pause("========== OKOKOK =============");
    */

  }

  void processWriteBodyForce(const int step,const double time,const BodyForce& body_force,const string& param_str) {

    // reduce force for body_force...

    double buf[8];
    double my_buf[8] = {body_force.volume,body_force.force[0],body_force.force[1],body_force.force[2],body_force.work,
                                          body_force.moment[0],body_force.moment[1],body_force.moment[2]};
    MPI_Reduce(my_buf,buf,8,MPI_DOUBLE,MPI_SUM,body_force.write_rank,mpi_comm);

    // do file io for force writing...

    int ierr = 0;

    if (mpi_rank == body_force.write_rank) {

      assert( body_force.moment_centroids.size()%3 == 0);
      int n_pts = body_force.moment_centroids.size()/3;

      //if moments about specific locations were reqested, report them
      //otherwise report moment about origin
      double (*moments)[3];
      if (n_pts>0){
        moments = new double[n_pts][3];
        for (int im=0; im < n_pts; ++im) {
          const double xp[3] = {body_force.moment_centroids[im*3+0],
                                body_force.moment_centroids[im*3+1],
                                body_force.moment_centroids[im*3+2]};

          const double fb[3] = {buf[1],buf[2],buf[3]};
          const double mp[3] = CROSS_PRODUCT(xp,fb);
          FOR_I3 moments[im][i] = buf[5+i] - mp[i];
        }
      }
      else{
        moments = new double[1][3]; //report moment about origin
        FOR_I3 moments[0][i] = buf[5+i];
        n_pts = 1;
      }

      ofstream ofile;
      size_t fwidth = 13;

      string name = body_force.name+".dat";

      if (MiscUtils::fileExists(name)) {
        ofile.open(name.c_str(),ios::app);
        if (!ofile.is_open()) {
         cerr << ERRSTART << "Error: could not open ADD_BODY_FORCE file: " << name << ERREND << endl;
         ierr = 1;
        }
      }
      else {
        MiscUtils::mkdir_for_file(name);
        ofile.open(name.c_str());
        if (!ofile.is_open()) {
         cerr << ERRSTART << "Error: could not open ADD_BODY_FORCE file: " << name << ERREND << endl;
         ierr = 1;
        }
        else {

         // on the first time, write the header...

         // the parameter...
         ofile << "# " << param_str << endl;

         // moment indexing
         if (body_force.moment_centroids.size()>0){
           for (int im=0; im < n_pts; ++im) {
             ofile << "# Moment m" << im << " ("<<body_force.moment_centroids[im*3+0]<<","<<body_force.moment_centroids[im*3+1]<<","<<body_force.moment_centroids[im*3+2] <<")" << endl;
           }
         }
         else{
           ofile << "# Moment m0 (0.0,0.0,0.0)" << endl;
         }

         // the column index...
         ofile << "# " <<
           setw(fwidth-2) << 1;
         for (int i = 2; i <= 7+(n_pts*3); ++i)
           ofile << setw(fwidth) << i;
         ofile << endl;

         // the variable name...
         ofile << "# " <<
           setw(fwidth-2) << "step" <<
           setw(fwidth) << "time" <<
           setw(fwidth) << "volume";

         ofile <<
           setw(fwidth) << "f-total-x" << setw(fwidth) << "f-total-y" << setw(fwidth) << "f-total-z";

         ofile <<
           setw(fwidth) << "work_total";

         for (int im=0; im < n_pts; ++im) {
           ofile <<
             setw(fwidth-9) << "m" << im << "-total-x" << setw(fwidth-9) << "m" << im << "-total-y" << setw(fwidth-9) << "m" << im << "-total-z";
         }

         ofile << endl;

        }
      }
      if (ofile.is_open()) {

        ofile <<
         setw(fwidth) << step <<
         setw(fwidth) << time <<
         setw(fwidth) << buf[0];

        FOR_I3 ofile <<
         setw(fwidth) << buf[1+i];

        ofile <<
         setw(fwidth) << buf[4];

        for (int im=0; im < n_pts; ++im) {
          FOR_I3 ofile <<
           setw(fwidth) << moments[im][i];
        }

        ofile << endl;

        ofile.close();
      }
      delete[] moments;
    }

    MPI_Bcast(&ierr,1,MPI_INT,body_force.write_rank,mpi_comm);

    if ( ierr != 0) {

      CERR( " > unable open the file for the ADD_BODY_FORCE request: " << param_str);

    }

  }

  // common flow solver special function evaluations

  virtual CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,const string& name,
                                                    list<CtiRegister::CtiData>& args, const bool b_eval_func) {

    using namespace CtiRegister;

    if ( name == "sn_prod" ) {

      // return the scalar buffer psi(ifa) = scalar * sign(n_fa[i],normal)

      if ( args.size() != 4)
        return CTI_DATA_ARG_COUNT;

      list<CtiData>::iterator arg = args.begin();

      if ( (arg->getType() != DN_DATA) || (arg->getTopology() != CV_DATA)) {

        return CTI_DATA_NOT_VALID;

      }

      double normal[3]; // parse to find the global normal specified...

      int ii = 0;
      for (list<CtiData>::iterator it = args.begin(); it != args.end(); ++it) {

        // skip the first entry ..
        if ( it == args.begin())
          continue;

        if ( it->getType() != D_DATA)
          return CTI_DATA_NOT_VALID;

        normal[ii++] = it->d();
      }

      // make this the unit normal..

      const double n_mag = MAG(normal);
      for (int i = 0; i < 3; ++i)
        normal[i] /= n_mag;

      double *v_ptr = createSignedFaD1Data(v);

      if ( b_eval_func ) {

        // since we don't know about whether the data has ghost information, we are
        // forced to do an interprocess reduction here ..

        double * arg_ptr = arg->getDNptr();
        for (int ifa = nfa_i; ifa < nfa; ++ifa) {

          const int icv0  = cvofa[ifa][0];
          const double dp = DOT_PRODUCT(n_fa[ifa], normal);
          v_ptr[ifa] = 0.5 * arg_ptr[icv0] * dp;

        }

        updateFaDataStart( v_ptr, SUBTRACT_DATA);

        // complete with the internal faces

        for (int ifa = 0; ifa < nfa_i; ++ifa) {

          const int icv0  = cvofa[ifa][0];
          const int icv1  = cvofa[ifa][1];
          const double dp = DOT_PRODUCT(n_fa[ifa], normal);
          v_ptr[ifa]      = 0.5* ( arg_ptr[icv0] + arg_ptr[icv1] )* dp;

        }

        updateFaDataFinish( v_ptr, SUBTRACT_DATA);
      }

      return CTI_DATA_OK;

    } else if ( name == "cfl") {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double *v_ptr = createCvD1Data(v);
      if ( b_eval_func) {
	calcCfl( v_ptr, dt);
      }

      return CTI_DATA_OK;

    }

    // this solver doesnt know how to evaluate this data -- go to StaticSolver

    return StaticSolver::funcEvalCtiData(v,name,args,b_eval_func);
  }

  // =============================================================================
  // CHT-related functions where the flow solver interacts with the cht solver...
  // =============================================================================





  // =============================
  // LpTracer-related functions...
  // =============================

  void setupLpMaterial() {

    FOR_PARAM_MATCHING("LP.MATERIAL") {

      int iarg = 0;
      string token = param->getString(iarg++);
      if (token == "LIQUID") {
        CtiLiquid * fuel = newCtiLiquid(&(*param),iarg);
        fuelVec.push_back(fuel);
        IF_RANK0 fuel->info(300,350,10);
      }
      else if (token == "SOLID") {
        CtiSolid * dust = newCtiSolid(&(*param),iarg);
        dustVec.push_back(dust);
        IF_RANK0 dust->info(300,350,10);
      }
      else {
        CERR("unrecognized LP.MATERIAL: " << token << ". Possible choices are:\n" <<
             "LP.MATERIAL LIQUID <liquid-tokens>\n" <<
             "LP.MATERIAL SOLID <solid-tokens>");
      }

    }

  }

  void setupLpInjectors() {

    IF_RANK0 cout << "FlowSolver::setupLpInjectors()" << endl;

    FOR_PARAM_MATCHING("LP_TRACER.INJECTOR") {

      string token = param->getString();
      if (token != "NAME") CERR("LP_TRACER.INJECTOR needs NAME. Unrecogniced token "<<token<<"\n");
      string name = param->getString(1);
      for (int ivec = 0; ivec < tracerInjectorVec.size(); ++ivec) {
        if (name == tracerInjectorVec[ivec]->getName())
          CERR("LP_TRACER.INJECTOR: " << name << " already exists. Please change the name.\n");
      }

      tracerInjectorVec.push_back(new LpInjectorTracer);

      LpInjectorTracer * injector = tracerInjectorVec.back();
      int injector_id = tracerInjectorVec.size()-1; // the index of the injector in the injector vector
      injector->initFromParams(&(*param),this,time,injector_id);

    }

    FOR_PARAM_MATCHING("LP_SOLID.INJECTOR") {

      string token = param->getString();
      if (token != "NAME") CERR("LP_SOLID.INJECTOR needs NAME. Unrecogniced token "<<token<<"\n");
      string name = param->getString(1);
      for (int ivec = 0; ivec < solidInjectorVec.size(); ++ivec) {
        if (name == solidInjectorVec[ivec]->getName())
          CERR("LP_SOLID.INJECTOR: " << name << " already exists. Please change the name.\n");
      }

      solidInjectorVec.push_back(new LpInjector);

      LpInjector * injector = solidInjectorVec.back();
      int injector_id = solidInjectorVec.size()-1; // the index of the injector in the injector vector
      injector->initFromParams(&(*param),this,time,fuelVec,dustVec,lpSolid->u,injector_id);

    }

    FOR_PARAM_MATCHING("LP_LIQUID.INJECTOR") {

      string token = param->getString();
      if (token != "NAME") CERR("LP_LIQUID.INJECTOR needs NAME. Unrecogniced token "<<token<<"\n");
      string name = param->getString(1);
      for (int ivec = 0; ivec < liquidInjectorVec.size(); ++ivec) {
        if (name == liquidInjectorVec[ivec]->getName())
          CERR("LP_LIQUID.INJECTOR: " << name << " already exists. Please change the name.\n");
      }

      liquidInjectorVec.push_back(new LpInjector);

      LpInjector * injector = liquidInjectorVec.back();
      int injector_id = liquidInjectorVec.size()-1; // the index of the injector in the injector vector
      injector->initFromParams(&(*param),this,time,fuelVec,dustVec,lpLiquid->u,injector_id);

    }

  }

  void addNewLp() {

    if (lpTracer) {

      vector<LpDataBase> lpDataVec;
      for (int ivec = 0; ivec < tracerInjectorVec.size(); ++ivec) {
        tracerInjectorVec[ivec]->addNewLp(lpDataVec,time,dt,true);
      }

      const int npnew = lpDataVec.size();

      if (npnew > 0) {
        const int npold = lpTracer->np;
        lpTracer->resize(npold+npnew);
        for (int ip = 0; ip < npnew; ++ip) {
          lpTracer->lp[npold+ip].icv          = lpDataVec[ip].icv;
          lpTracer->lp[npold+ip].flag         = lpDataVec[ip].flag;
          lpTracer->lp[npold+ip].dp           = lpDataVec[ip].dp;
          FOR_I3 lpTracer->lp[npold+ip].xp[i] = lpDataVec[ip].xp[i];
          FOR_I3 lpTracer->lp[npold+ip].up[i] = lpTracer->u[lpTracer->lp[npold+ip].icv][i]; // Maybe not needed
        }
      }

    }

    if (lpSolid) {

      vector<LpData> lpDataVec;
      for (int ivec = 0; ivec < solidInjectorVec.size(); ++ivec) {
        solidInjectorVec[ivec]->addNewLp(lpDataVec,time,dt,true);
      }

      const int npnew = lpDataVec.size();

      if (npnew > 0) {
        const int npold = lpSolid->np;
        lpSolid->resize(npold+npnew);
        for (int ip = 0; ip < npnew; ++ip) {
          lpSolid->lp[npold+ip].icv           = lpDataVec[ip].icv;
          lpSolid->lp[npold+ip].flag          = lpDataVec[ip].flag;
          FOR_I3 lpSolid->lp[npold+ip].xp[i]  = lpDataVec[ip].xp[i];
          FOR_I3 lpSolid->lp[npold+ip].up[i]  = lpDataVec[ip].up[i];
          lpSolid->lp[npold+ip].mp            = lpDataVec[ip].mp;
          lpSolid->lp[npold+ip].dp            = lpDataVec[ip].dp;
          lpSolid->lp[npold+ip].Tp            = lpDataVec[ip].Tp;
          lpSolid->lp[npold+ip].npar          = 1.0;
        }
      }

    }

    if (lpLiquid) {

      vector<LpData> lpDataVec;
      for (int ivec = 0; ivec < liquidInjectorVec.size(); ++ivec) {
        liquidInjectorVec[ivec]->addNewLp(lpDataVec,time,dt,true);
      }

      const int npnew = lpDataVec.size();

      if (npnew > 0) {
        const int npold = lpLiquid->np;
        lpLiquid->resize(npold+npnew);
        for (int ip = 0; ip < npnew; ++ip) {
          lpLiquid->lp[npold+ip].icv            = lpDataVec[ip].icv;
          lpLiquid->lp[npold+ip].flag           = lpDataVec[ip].flag;
          FOR_I3 lpLiquid->lp[npold+ip].xp[i]   = lpDataVec[ip].xp[i];
          FOR_I3 lpLiquid->lp[npold+ip].up[i]   = lpDataVec[ip].up[i];
          lpLiquid->lp[npold+ip].mp             = lpDataVec[ip].mp;
          lpLiquid->lp[npold+ip].dp             = lpDataVec[ip].dp;
          lpLiquid->lp[npold+ip].Tp             = lpDataVec[ip].Tp;
          lpLiquid->lp[npold+ip].npar           = 1.0;
          lpLiquid->lp[npold+ip].tbu            = 0.0;
        }
      }

    }

  }

  void processLp() {

    if (lpTracer) {
      if ((step%check_interval==0)&&(mpi_rank == 0)) cout << " > +++++++++++++++++++++" << endl << " > process tracers" << endl << " > +++++++++++++++++++++" << endl;

      // this also implements the bcs...
      relocateLp<LpTracer,LpTracerState>(lpTracer);

      rebuildLpocv<LpTracer,LpTracerState>(lpTracer);

      if (step%check_interval == 0) {
        reportContainmentLp<LpTracer>(lpTracer);
        lpTracer->report();
      }

    }

    if (lpSolid) {
      if ((step%check_interval==0)&&(mpi_rank == 0)) cout << " > +++++++++++++++++++++" << endl << " > process solid lp" << endl << " > +++++++++++++++++++++" << endl;

      // this also implements the bcs...
      relocateLp<LpSolid,LpSolidState>(lpSolid);

      rebuildLpocv<LpSolid,LpSolidState>(lpSolid);

      int ntrash = recycleLp<LpSolidState>(lpSolid);
      if ((step%check_interval == 0)&&(mpi_rank==0)&&(ntrash>0)) cout << " > # of trashed solid particles: " << ntrash << endl;

      if (step%check_interval == 0) {
        reportContainmentLp<LpSolid>(lpSolid);
        lpSolid->report();
      }

    }

    if (lpLiquid) {
      if ((step%check_interval==0)&&(mpi_rank == 0)) cout << " > +++++++++++++++++++++" << endl << " > process liquid lp" << endl << " > +++++++++++++++++++++" << endl;

      // this also implements the bcs...
      relocateLp<LpLiquid,LpLiquidState>(lpLiquid);

      rebuildLpocv<LpLiquid,LpLiquidState>(lpLiquid);

      int ntrash = recycleLp<LpLiquidState>(lpLiquid);
      if ((step%check_interval == 0)&&(mpi_rank==0)&&(ntrash>0)) cout << " > # of trashed liquid particles: " << ntrash << endl;

      if (step%check_interval == 0) {
        reportContainmentLp<LpLiquid>(lpLiquid);
        lpLiquid->report();
      }

    }

  }

  void buildLpocv() {

    if (lpTracer)
      rebuildLpocv<LpTracer,LpTracerState>(lpTracer);

    if (lpSolid)
      rebuildLpocv<LpSolid,LpSolidState>(lpSolid);

    if (lpLiquid)
      rebuildLpocv<LpLiquid,LpLiquidState>(lpLiquid);

  }

  void setLp0() {

    if (lpTracer) {
      for (int ip = 0; ip < lpTracer->size(); ++ip) {
        FOR_I3 lpTracer->lp[ip].xp0[i] = lpTracer->lp[ip].xp[i];
      }
    }

    if (lpSolid) {
      for (int ip = 0; ip < lpSolid->size(); ++ip) {
        FOR_I3 lpSolid->lp[ip].xp0[i] = lpSolid->lp[ip].xp[i];
      }
    }

    if (lpLiquid) {
      for (int ip = 0; ip < lpLiquid->size(); ++ip) {
        FOR_I3 lpLiquid->lp[ip].xp0[i] = lpLiquid->lp[ip].xp[i];
        lpLiquid->lp[ip].mp0 = lpLiquid->lp[ip].mp;
      }
    }

  }

  void updateLpDp() {
    if (lpLiquid) {
      for (int ip = 0; ip < lpLiquid->size(); ++ip) {
        LpLiquidState * lpL = &lpLiquid->lp[ip];
        const int material_id = getLpMaterial(lpL->flag);
        assert(material_id>=0); assert(material_id<fuelVec.size());
        CtiLiquid * fuel = fuelVec[material_id];
        lpL->dp = pow(lpL->mp/(M_PI/6.0*lpL->npar*fuel->calcRho(lpL->Tp)) , 1./3.);
      }
    }
  }

  int getLpMaterial(const int flag) const {

    const char * material_keep_injector = (char *)&(flag);
    const int material = int(material_keep_injector[0]);
    return material;

  }

  int getLpKeep(const int flag) const {

    const char * material_keep_injector = (char *)&(flag);
    const int keep = int(material_keep_injector[1]);
    return keep;

  }

  int getLpInjector(const int flag) const {

    const short int * material_keep_injector = (short int *)&(flag);
    const int injector = int(material_keep_injector[1]);
    return injector;

  }

  void setLpKeep(int &flag, const int keep) {

    char * material_keep_injector = (char *)&(flag);
    material_keep_injector[1] = char(keep);

  }

  template<class T>
  int recycleLp(LpContainer<T> * lpContainer) {
    int ip_new = 0;
    for (int ip = 0; ip < lpContainer->size(); ++ip) {
      const int keep = getLpKeep(lpContainer->lp[ip].flag);
      if (keep == KEEP) {
        if (ip != ip_new) {
          lpContainer->lp[ip_new] = lpContainer->lp[ip];
        }
        ip_new++;
      } else if (keep == TRASH) {
      } else {
        CERR("ERROR:particle flag not recognized" << " , particle: " << ip << " , flag: " << keep << "\n");
      }
    }
    int myNTrash = lpContainer->size() - ip_new;
    int NTrash;
    MPI_Reduce(&myNTrash,&NTrash,1,MPI_INT,MPI_SUM,0,mpi_comm);
    lpContainer->resize(ip_new);
    return NTrash;
  }

  template<class T,class S> // T: lp container, S: lp state
  void relocateLp(T * lpContainer) {

    if ((step%check_interval==0)&&(mpi_rank == 0)) {
      cout << " > relocateLp()";
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
              lpContainer->lp[ip].packRT(send_buf+send_disp[rank]*S::data_size(),per_R,per_t);
            }
            else if (b_per_R) {
              lpContainer->lp[ip].packR(send_buf+send_disp[rank]*S::data_size(),per_R);
            }
            else if (b_per_t) {
              lpContainer->lp[ip].packT(send_buf+send_disp[rank]*S::data_size(),per_t);
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
            int ist_ss_closest;
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
                    ist_ss_closest = ist_ss;
                    d2_closest = d2;
                    FOR_I3 dx_closest[i] = dx[i];
                    FOR_I3 normal_closest[i] = normal[i];
                    dist_closest = dist;
                  }
                }
              }
            }
            assert(ibf_closest >= 0);
            if (dist_closest <= 0.0) {
	      // dist with negative sign means we are outside the fluid volume, so apply bc...
              const int izn = zone_bf[ibf_closest];
              assert((izn >= 0)&&(izn < lpContainer->bcVec.size()));
              const bool keep = lpContainer->bcVec[izn]->applyBc(&lpContainer->lp[ip],ibf_closest,ist_ss_closest,dx_closest,normal_closest);
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

  template<class T,class S>
  void rebuildLpocv(T * lpContainer) {

    // use this as an opportunity to trash particles that have their
    // icv set to -1, and to reorder all particles in lpocv order...

    if ((step%check_interval==0)&&(mpi_rank == 0)) cout << " > rebuildLpocv()" << endl;

    if (lpContainer->lpocv_i == NULL) lpContainer->lpocv_i = new int[ncv+1];

    FOR_ICV lpContainer->lpocv_i[icv+1] = 0;

    for (int ip = 0; ip < lpContainer->np; ++ip) {
      // another routine may have set icv == -1, which means "TRASH"
      // the particle, so skip it in this counting...
      if (lpContainer->lp[ip].icv >= 0) {
        assert(lpContainer->lp[ip].icv < ncv);
        ++lpContainer->lpocv_i[lpContainer->lp[ip].icv+1];
      }
    }

    // turn it into a CSR index for the ip's...
    lpContainer->lpocv_i[0] = 0;
    FOR_ICV lpContainer->lpocv_i[icv+1] += lpContainer->lpocv_i[icv];
    const int np_new = lpContainer->lpocv_i[ncv];

    // this reordering routine works "in-place" in the existing lp memory...
    S tmp1,tmp2;
    for (int ip = 0; ip < lpContainer->np; ++ip) {
      if (lpContainer->lp[ip].icv >= 0) {
        // where does this one go?...
        int ip_target = lpContainer->lpocv_i[lpContainer->lp[ip].icv]++;
        lpContainer->lp[ip].icv = -lpContainer->lp[ip].icv-2; // move to -2 indexing...
        if (ip_target != ip) {
          // this particle is in the wrong place...
          if (lpContainer->lp[ip_target].icv < 0) {
            // the particle at ip_target has been moved or trashed, so just
            // overwrite it...
            lpContainer->lp[ip_target] = lpContainer->lp[ip];
          }
          else {
            tmp1 = lpContainer->lp[ip_target];
            lpContainer->lp[ip_target] = lpContainer->lp[ip];
            // now what to do with tmp1?...
            ip_target = lpContainer->lpocv_i[tmp1.icv]++;
            tmp1.icv = -tmp1.icv-2; // move icv to -2 indexing...
            while (lpContainer->lp[ip_target].icv >= 0) {
              tmp2 = lpContainer->lp[ip_target];
              lpContainer->lp[ip_target] = tmp1;
              tmp1 = tmp2;
              ip_target = lpContainer->lpocv_i[tmp1.icv]++;
              tmp1.icv = -tmp1.icv-2; // move icv to -2 indexing...
            }
            lpContainer->lp[ip_target] = tmp1;
          }
        }
      }
    }

    // return icv to positive index...
    lpContainer->resize(np_new);
    for (int ip = 0; ip < lpContainer->np; ++ip) {
      lpContainer->lp[ip].icv = -lpContainer->lp[ip].icv-2;
      assert((lpContainer->lp[ip].icv >= 0)&&(lpContainer->lp[ip].icv < ncv));
    }

    // and rewind csr...
    for (int icv = ncv-1; icv > 0; --icv)
      lpContainer->lpocv_i[icv] = lpContainer->lpocv_i[icv-1];
    lpContainer->lpocv_i[0] = 0;

    // finally, check!...
    FOR_ICV {
      for (int ip = lpContainer->lpocv_i[icv]; ip != lpContainer->lpocv_i[icv+1]; ++ip) {
        assert(lpContainer->lp[ip].icv == icv);
      }
    }

  }

  template<class T> // T: lp container
  void reportContainmentLp(T * lpContainer) {

    int8 my_count[5] = { lpContainer->np, 0, 0, 0, 0 };

    for (int ip = 0; ip < lpContainer->np; ++ip) {
      const int icv = lpContainer->lp[ip].icv;
      assert((icv >= 0)&&(icv < ncv));
      // check that we are closer to icv than to any neighbor:
      // use the same tolerances as in the relocateLp routine...
      const double d2_icv = DIST2(lpContainer->lp[ip].xp,x_vv[icv])*0.999999;
      bool nbr_closer = false;
      for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) { // note we skip self here
	const int icv_nbr = cvocv_v[coc];
	const double d2_icv_nbr = DIST2(lpContainer->lp[ip].xp,x_vv[icv_nbr]);
	if (d2_icv_nbr < d2_icv) {
	  nbr_closer = true;
	  break;
	}
      }
      if (bfocv_i[icv+1] == bfocv_i[icv]) {
	if (nbr_closer) ++my_count[1];
      }
      else {
	if (pointIsInsideCvBoundary(lpContainer->lp[ip].xp,icv)) {
	  if (nbr_closer) ++my_count[2];
	}
	else {
	  if (nbr_closer) ++my_count[4];
	  else ++my_count[3];
	}
      }
    }

    int8 count[5];
    MPI_Reduce(my_count,count,5,MPI_INT8,MPI_SUM,0,mpi_comm);
    if (mpi_rank == 0) {
      cout << " > lp containment report for " << count[0] << " parcels:" << endl;
      cout << " > ip inside internal cv but closer to nbr : " << count[1] << endl;
      cout << " > ip inside boundary cv but closer to nbr : " << count[2] << endl;
      cout << " > ip outside boundary cv                  : " << count[3] << endl;
      cout << " > ip outside boundary cv and closer to nbr: " << count[4] << endl;
    }

  }

  inline bool pointIsInsideCvBoundary(const double xp[3],const int icv) const {
    if (bfocv_i[icv+1] == bfocv_i[icv]) {
      return true; // no boundary faces, so assume we are inside...
    }
    else {
      // this one has boundary faces. We need a length scale...
      const double tol2 = 1.0E-10*pow(vol_cv[icv],2.0/3.0);
      // and we need the closest ibf,ist_ss...
      int ibf_closest = -1;
      double d2_closest,dist_closest;
      for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
        const int ibf = bfocv_v[boc];
        for (int sob = sstobf_i[ibf]; sob != sstobf_i[ibf+1]; ++sob) {
          const int ist_ss = sstobf_v[sob];
          // the corner coords of this sub-surface tri are (recall that the
          // subsurface tri already has any periodic transfrom taken care of)...
          const double * const x0 = subSurface->xp[subSurface->spost[ist_ss][0]];
          const double * const x1 = subSurface->xp[subSurface->spost[ist_ss][1]];
          const double * const x2 = subSurface->xp[subSurface->spost[ist_ss][2]];
          double xs[3]; getClosestPointOnTriRobust(xs,xp,x0,x1,x2);
          const double dx[3] = DIFF(xs,xp);
          const double d2 = DOT_PRODUCT(dx,dx);
          if ((ibf_closest == -1)||(d2 < d2_closest+tol2)) {
            const double normal[3] = TRI_NORMAL_2(x0,x1,x2);
            const double normal_mag = MAG(normal);
            assert(normal_mag != 0); // the subsurface should not have zero-area tris!
            const double dist = DOT_PRODUCT(dx,normal)/normal_mag;
            if ((ibf_closest == -1)||(d2 < d2_closest-tol2)||(fabs(dist) > fabs(dist_closest))) {
              ibf_closest = ibf;
              d2_closest = d2;
              dist_closest = dist;
            }
          }
        }
      }
      assert(ibf_closest >= 0);
      if (dist_closest <= 0.0) {
        // dist with negative sign means we are outside the fluid volume...
        return false;
      }
      else {
        // we are inside the boundary, so keep...
        return true;
      }
    }
  }

  void registerLpBoundaryConditions() {

    if (lpTracer) {
      assert(lpTracer->bcVec.empty());
      lpTracer->bcVec.resize(bfZoneVec.size());
      for (int izn = 0; izn < lpTracer->bcVec.size(); ++izn)
        lpTracer->bcVec[izn] = NULL;

      FOR_PARAM_MATCHING("LP_TRACER.BC") {
        int iarg = 0;
        const string zone_name = param->getString(iarg++);
        map<const string,int>::iterator iter = bfZoneNameMap.find(zone_name);
        if (iter == bfZoneNameMap.end()) {
          CERR("LP_TRACER.BC: unrecognized zone name: " << zone_name);
        }
        else {
          const int izn = iter->second;
          if (lpTracer->bcVec[izn] != NULL) {
            CERR("LP_TRACER.BC: bc already set for zone_name: " << zone_name);
          }
          else {
            const string bc_type = param->getString(iarg++);
            if (bc_type == "OUTLET") {
              lpTracer->bcVec[izn] = new LpTracerBcOutlet(&bfZoneVec[izn],this);
            }
            else if (bc_type == "BOUNCE") {
              lpTracer->bcVec[izn] = new LpTracerBcBounce(&bfZoneVec[izn],this);
            }
            else {
              CERR("LP_TRACER.BC: unrecognized bc type: " << bc_type);
            }
          }
        }
      }

      // and apply default bc "OUTLET" to remaining...

      const string bc_type_default = getStringParam("LP_TRACER.BC_DEFAULT","OUTLET");
      if (bc_type_default == "OUTLET") {
        for (int izn = 0; izn < lpTracer->bcVec.size(); ++izn)
          if (lpTracer->bcVec[izn] == NULL)
            lpTracer->bcVec[izn] = new LpTracerBcOutlet(&bfZoneVec[izn],this);
      }
      else if (bc_type_default == "BOUNCE") {
        for (int izn = 0; izn < lpTracer->bcVec.size(); ++izn)
          if (lpTracer->bcVec[izn] == NULL)
            lpTracer->bcVec[izn] = new LpTracerBcBounce(&bfZoneVec[izn],this);
      }
      else {
        CERR("LP_TRACER.BC_DEFAULT: unrecognized bc type: " << bc_type_default);
      }
    }

    if (lpSolid) {
      assert(lpSolid->bcVec.empty());
      lpSolid->bcVec.resize(bfZoneVec.size());
      for (int izn = 0; izn < lpSolid->bcVec.size(); ++izn)
        lpSolid->bcVec[izn] = NULL;

      FOR_PARAM_MATCHING("LP_SOLID.BC") {
        int iarg = 0;
        const string zone_name = param->getString(iarg++);
        map<const string,int>::iterator iter = bfZoneNameMap.find(zone_name);
        if (iter == bfZoneNameMap.end()) {
          CERR("LP_SOLID.BC: unrecognized zone name: " << zone_name);
        }
        else {
          const int izn = iter->second;
          if (lpSolid->bcVec[izn] != NULL) {
            CERR("LP_SOLID.BC: bc already set for zone_name: " << zone_name);
          }
          else {
            const string bc_type = param->getString(iarg++);
            if (bc_type == "OUTLET") {
              lpSolid->bcVec[izn] = new LpTracerBcOutlet(&bfZoneVec[izn],this);
            }
            else if (bc_type == "BOUNCE") {
              lpSolid->bcVec[izn] = new LpTracerBcBounce(&bfZoneVec[izn],this);
            }
            else {
              CERR("LP_SOLID.BC: unrecognized bc type: " << bc_type);
            }
          }
        }
      }

      // and apply default bc "OUTLET" to remaining...

      const string bc_type_default = getStringParam("LP_SOLID.BC_DEFAULT","BOUNCE");
      if (bc_type_default == "OUTLET") {
        for (int izn = 0; izn < lpSolid->bcVec.size(); ++izn)
          if (lpSolid->bcVec[izn] == NULL)
            lpSolid->bcVec[izn] = new LpTracerBcOutlet(&bfZoneVec[izn],this);
      }
      else if (bc_type_default == "BOUNCE") {
        for (int izn = 0; izn < lpSolid->bcVec.size(); ++izn)
          if (lpSolid->bcVec[izn] == NULL)
            lpSolid->bcVec[izn] = new LpTracerBcBounce(&bfZoneVec[izn],this);
      }
      else {
        CERR("LP_SOLID.BC_DEFAULT: unrecognized bc type: " << bc_type_default);
      }
    }

    if (lpLiquid) {
      assert(lpLiquid->bcVec.empty());
      lpLiquid->bcVec.resize(bfZoneVec.size());
      for (int izn = 0; izn < lpLiquid->bcVec.size(); ++izn)
        lpLiquid->bcVec[izn] = NULL;

      FOR_PARAM_MATCHING("LP_LIQUID.BC") {
        int iarg = 0;
        const string zone_name = param->getString(iarg++);
        map<const string,int>::iterator iter = bfZoneNameMap.find(zone_name);
        if (iter == bfZoneNameMap.end()) {
          CERR("LP_LIQUID.BC: unrecognized zone name: " << zone_name);
        }
        else {
          const int izn = iter->second;
          if (lpLiquid->bcVec[izn] != NULL) {
            CERR("LP_LIQUID.BC: bc already set for zone_name: " << zone_name);
          }
          else {
            const string bc_type = param->getString(iarg++);
            if (bc_type == "OUTLET") {
              lpLiquid->bcVec[izn] = new LpTracerBcOutlet(&bfZoneVec[izn],this);
            }
            else if (bc_type == "BOUNCE") {
              lpLiquid->bcVec[izn] = new LpTracerBcBounce(&bfZoneVec[izn],this);
            }
            else {
              CERR("LP_LIQUID.BC: unrecognized bc type: " << bc_type);
            }
          }
        }
      }

      // and apply default bc "OUTLET" to remaining...

      const string bc_type_default = getStringParam("LP_LIQUID.BC_DEFAULT","BOUNCE");
      if (bc_type_default == "OUTLET") {
        for (int izn = 0; izn < lpLiquid->bcVec.size(); ++izn)
          if (lpLiquid->bcVec[izn] == NULL)
            lpLiquid->bcVec[izn] = new LpTracerBcOutlet(&bfZoneVec[izn],this);
      }
      else if (bc_type_default == "BOUNCE") {
        for (int izn = 0; izn < lpLiquid->bcVec.size(); ++izn)
          if (lpLiquid->bcVec[izn] == NULL)
            lpLiquid->bcVec[izn] = new LpTracerBcBounce(&bfZoneVec[izn],this);
      }
      else {
        CERR("LP_LIQUID.BC_DEFAULT: unrecognized bc type: " << bc_type_default);
      }
    }

  }

};

inline void reportBcErrors(vector<pair<string,string> >& errors) {
  const int nerr = int(errors.size());
  if ( nerr > 0) {
    // report the errors from rank 0..
    if ( mpi_rank == 0 ) {
      cout << " -------------------------- " << endl;
      cout << " Errors detected in the boundary conditions setup: " << endl;
      cout << endl;
      for (int i =0; i < nerr; ++i) {
        if ( errors[i].second == "") {
          cout << " > zone: " << errors[i].first << " --> MISSING bc specification " << endl;
        } else {
          cout << " > zone: " << errors[i].first << " --> unknown bc type " << errors[i].second << endl;
        }
      }
      cout << " -------------------------- " << endl;
      cout.flush();
    }
    throw(0);
  }
}
#endif
