#ifndef IDEALGASSOLVERBCS_HPP
#define IDEALGASSOLVERBCS_HPP

#include "BcCommon.hpp"
#include "States.hpp"

// need to provide a forward declaration for the
// add boundary flux interface provided below.
// alternatively, we could make the boundary condition
// classes embedded inside of the IdealGasSolver
// note that the defn of the BC interface is in BcCommon.hpp

class IdealGasSolver;
typedef BcTemplate<IdealGasSolver,IdealGasRhs> IdealGasBc;

class DistributedDataExchanger;
class DataExchanger;

//
// the simplest ideal gas bc is one that merely
// requires you to write an addBoundaryFlux routine
// that introduces your information into the domain..
//

class SimpleIdealGasBc : public IdealGasBc {
public:

  SimpleIdealGasBc(BfZone* p, IdealGasSolver* s): IdealGasBc(p,s) {}

  void calcRhs(const double time, const int rk_stage) {}
  void rk3Step(const double dt, const int rk_stage) {}

  void initData() {
    assert( mf == NULL); mf = new double[zone_ptr->nbf];
  }

  void initialHook() {}

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability. Here we assume that SimpleIdealGasBc bcs do not store
  // or restore any data by default...
  virtual int storeData(double * buf) const { return 0; }
  virtual int restoreData(double * buf) { return 0; }

  // note there will be no force information computed frm this bc..
  void force_bf(double (*f_bf)[9]) const {}

  virtual void addBoundaryFlux(IdealGasRhs* rhs) const {
    assert(0); // user must provide an implementation here..
  }
};

class SlipWallBc : public IdealGasBc {
public:

  SlipWallBc(BfZone* p, IdealGasSolver* s): IdealGasBc(p,s) {}

  void calcRhs(const double time, const int rk_stage) {}
  void rk3Step(const double dt, const int rk_stage) {}
  void initData() {}
  void initialHook() {}

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const { return 0; }
  int restoreData(double * buf) { return 0; }

  void force_bf(double (*f_bf)[9]) const {force_bf_slip_wall(f_bf,this);}

  void addBoundaryFlux(IdealGasRhs* rhs) const;
};

class WallAdiabatic : public IdealGasBc {
public:

  double (*u_bc)[3];
  WallAdiabatic(BfZone* p, IdealGasSolver* s) : IdealGasBc(p,s), u_bc(NULL) {

    funcEvalList.push_back(p->getName()+":tau_wall()");
    funcEvalList.push_back(p->getName()+":y_plus()");
    funcEvalList.push_back(p->getName()+":moment()");
    funcEvalList.push_back(p->getName()+":force()");
  }

  ~WallAdiabatic() { DELETE(u_bc); }

  void calcRhs(const double time, const int rk_stage){}
  void rk3Step(const double dt, const int rk_stage) {}
  void initData();
  void initialHook() {}

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const { return 0; }
  int restoreData(double * buf) { return 0; }

  void force_bf(double (*f_bf)[9]) const {force_bf_wall_new(f_bf,this);}

  void addBoundaryFlux(IdealGasRhs* rhs) const;
  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args,const bool b_eval_func);
};

class WmExchange : public IdealGasBc {
public:

  double * u1;
  double * T1;

  double T_wall_;
  double * tau_wall;
  double * q_wall;
  double * T_wall;

  double (*unit_n)[3];
  DataExchanger* exchanger;
  double del_exch_ratio;


  WmExchange(BfZone* p, IdealGasSolver* s) : IdealGasBc(p,s) {

    u1       = NULL; zone_ptr->registerBfData(u1    , "u1", READWRITE_DATA);
    T1       = NULL; zone_ptr->registerBfData(T1    , "T1", READWRITE_DATA);
    tau_wall = NULL; zone_ptr->registerBfData(tau_wall,"tau_wall", CAN_WRITE_DATA);
    q_wall   = NULL; zone_ptr->registerBfData(q_wall  ,"q_wall", CAN_WRITE_DATA);
    T_wall   = NULL; zone_ptr->registerBfData(T_wall  ,"T_wall", CAN_WRITE_DATA);

    unit_n    = NULL;
    exchanger = NULL;

    T_wall_        = 0.0;
    del_exch_ratio = 0.0;

    // add to funcEvalList

    funcEvalList.push_back(p->getName()+":tau_wall()");
    funcEvalList.push_back(p->getName()+":y_plus()");
    funcEvalList.push_back(p->getName()+":q_wall()");
    funcEvalList.push_back(p->getName()+":T_wall()");

  }


  ~WmExchange();

  void calcRhs(const double time, const int rk_stage);
  void rk3Step(const double dt, const int rk_stage) {}
  void initData();
  void initialHook();

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const;
  int restoreData(double * buf);

  void force_bf(double (*f_bf)[9]) const {force_bf_wm(f_bf,this);}

  void addBoundaryFlux(IdealGasRhs* rhs) const;
  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func);
};


class WmAdiabatic : public IdealGasBc {
public:

  double * u1;
  double * tau_wall;
  double * T_wall;
  //double u_bc[3];
  double (*u_bc)[3];

  // for equilibrium wall model...
  double relax,y_plus0,tol;

  WmType wm_type;

  WmAdiabatic(BfZone* p, IdealGasSolver* s, WmType _type) : IdealGasBc(p,s) {

    // XXX note that this bc currently assumes a no-slip (non-moving) wall.
    //for (int i = 0; i < 3; ++i)
    //  u_bc[i] = 0.0;

    // allow the wall to move
    u_bc     = NULL;

    u1       = NULL; zone_ptr->registerBfData(u1   , "u1", READWRITE_DATA);
    tau_wall = NULL; zone_ptr->registerBfData(tau_wall, "tau_wall", READWRITE_DATA);

    wm_type = _type;

    T_wall = NULL;
    if ( wm_type == WM_EQUILIBRIUM_ADIABATIC)
      zone_ptr->registerBfData(T_wall, "T_wall", READWRITE_DATA);

    relax = 0.7;
    y_plus0 = 1.0;
    tol = 1.0E-4;

    // add to funcEvalList
    funcEvalList.push_back(p->getName()+":tau_wall()");
    funcEvalList.push_back(p->getName()+":y_plus()");
  }

  ~WmAdiabatic();

  void calcRhs(const double time, const int rk_stage);
  void rk3Step(const double dt, const int rk_stage) {}
  void initData();
  void initialHook();

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const;
  int restoreData(double * buf);

  void force_bf(double (*f_bf)[9]) const {force_bf_wm(f_bf,this);}

  void addBoundaryFlux(IdealGasRhs* rhs) const;
  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func);

};

class WmAlgConductive : public IdealGasBc {
public:
  double * u1;
  double * tau_wall;
  double * q_wall;
  double * T_bc;
  double u_bc[3];
  double T_end_bc;
  vector<double> layer_length_bc;
  vector<double> thermal_conductivity_bc;
  double thermal_resistance_external_times_area;

  WmAlgConductive(BfZone* p, IdealGasSolver* s) : IdealGasBc(p,s) {
    // XXX note that this bc currently assumes a no-slip (non-moving) wall.
    FOR_I3 u_bc[i] = 0.0;

    u1        = NULL;   zone_ptr->registerBfData(u1,        "u1",         READWRITE_DATA);
    tau_wall  = NULL;   zone_ptr->registerBfData(tau_wall,  "tau_wall",   READWRITE_DATA);
    q_wall    = NULL;   zone_ptr->registerBfData(q_wall,    "q_wall",     READWRITE_DATA);
    T_bc      = NULL;   zone_ptr->registerBfData(T_bc,      "T_bc",       READWRITE_DATA);

    //add to funcEvalList
    funcEvalList.push_back(p->getName()+":tau_wall()");
    funcEvalList.push_back(p->getName()+":y_plus()");
    funcEvalList.push_back(p->getName()+":q_wall()");
    funcEvalList.push_back(p->getName()+":T_bc()");

    // Make sure there's a problem if these values don't get set...
    thermal_resistance_external_times_area = 0.0;
  }

  ~WmAlgConductive();

  void calcRhs(const double time, const int rk_stage);
  void rk3Step(const double dt, const int rk_stage) {}
  void initData();
  void initialHook();

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const;
  int restoreData(double * buf);

  void force_bf(double (*f_bf)[9]) const {force_bf_wm(f_bf, this);}

  void addBoundaryFlux(IdealGasRhs* rhs) const;
  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func);
};

class WallIsothermal : public IdealGasBc {
public:

  double u_bc[3];
  double T_bc;
  WallIsothermal(BfZone* p, IdealGasSolver* s) : IdealGasBc(p,s) {
    T_bc = 0.0;
    for (int i =0; i < 3; ++i) u_bc[i] = 0.0;

    // add to funcEvalList
    funcEvalList.push_back(p->getName()+":tau_wall()");
    funcEvalList.push_back(p->getName()+":y_plus()");
    funcEvalList.push_back(p->getName()+":q_wall()");
  }

  void calcRhs(const double time, const int rk_stage) {}
  void rk3Step(const double dt, const int rk_stage) {}
  void initData();
  void initialHook() {}

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const { return 0; }
  int restoreData(double * buf) { return 0; }

  void force_bf(double (*f_bf)[9]) const {force_bf_wall(f_bf,this);}

  void addBoundaryFlux(IdealGasRhs* rhs) const;
  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func);

};

class WallCht : public IdealGasBc {
private:
  // fem point data for this zone...
  int np;
  double * Tp;
  double * qp; // heat flow (flux times area)...
  int * ipobf_i;
  int * ipobf_v;
  double * ipobf_wgt;
public:
  WallCht(BfZone* p, IdealGasSolver* s) : IdealGasBc(p,s) {
    np = 0;
    Tp = NULL;
    qp = NULL;
    ipobf_i = NULL;
    ipobf_v = NULL;
    ipobf_wgt = NULL;
    // add to funcEvalList
    funcEvalList.push_back(p->getName()+":tau_wall()");
    funcEvalList.push_back(p->getName()+":y_plus()");
    funcEvalList.push_back(p->getName()+":q_wall()");
    funcEvalList.push_back(p->getName()+":T_wall()");
  }
  ~WallCht() {
    DELETE(Tp);
    DELETE(qp);
    DELETE(ipobf_i);
    DELETE(ipobf_v);
    DELETE(ipobf_wgt);
  }
  void calcRhs(const double time, const int rk_stage) {}
  void rk3Step(const double dt, const int rk_stage) {}
  void initData();
  void initialHook() {}

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const { return 0; }
  int restoreData(double * buf) { return 0; }

  void force_bf(double (*f_bf)[9]) const {
    //force_bf_wall(f_bf,this);
    // for now, do nothing...
  }
  void addBoundaryFlux(IdealGasRhs* rhs) const;
  void query(const string& param_str) const;
  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func);
};


class WmIsothermal : public IdealGasBc {
public:

  double * u1;
  double * tau_wall;
  double * q_wall;
  double u_bc[3];
  double T_bc;

  // for equilibrium wall model...
  double relax,y_plus0,tol;

  WmType wm_type;

  WmIsothermal(BfZone* p, IdealGasSolver* s, WmType _type) : IdealGasBc(p,s) {

    // XXX note that this bc currently assumes a no-slip (non-moving wall)
    for (int i =0; i < 3; ++i)
      u_bc[i] = 0.0;

    T_bc  = 0.0;

    u1        = NULL; zone_ptr->registerBfData(u1,  "u1", READWRITE_DATA);
    tau_wall  = NULL; zone_ptr->registerBfData(tau_wall, "tau_wall", READWRITE_DATA);
    q_wall    = NULL; zone_ptr->registerBfData(q_wall,"q_wall", READWRITE_DATA);

    wm_type   = _type;

    relax = 0.7;
    y_plus0 = 1.0;
    tol = 1.0E-4;

    // add to funcEvalList
    funcEvalList.push_back(p->getName()+":tau_wall()");
    funcEvalList.push_back(p->getName()+":y_plus()");
    funcEvalList.push_back(p->getName()+":q_wall()");
  }

  ~WmIsothermal();

  void calcRhs(const double time, const int rk_stage);
  void rk3Step(const double dt, const int rk_stage) {}
  void initData();
  void initialHook();

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const;
  int restoreData(double * buf);

  void force_bf(double (*f_bf)[9]) const {force_bf_wm(f_bf,this);}

  void addBoundaryFlux(IdealGasRhs* rhs) const;
  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func);
};

class WallConductive : public IdealGasBc {
public:

  // LOOK AT IMPLEMENTATION in premixed solver and fix

  double u_bc[3];
  double T_end_bc;
  vector<double> thermal_conductivity_bc;
  vector<double> layer_length_bc;
  double thermal_resistance_external_times_area;

  WallConductive(BfZone* p, IdealGasSolver* s) : IdealGasBc(p,s) {
    T_end_bc = 0.0;
    thermal_resistance_external_times_area = 0.0;
    for (int i =0; i < 3; ++i) u_bc[i] = 0.0;
  }

  void calcRhs(const double time, const int rk_stage) {}
  void rk3Step(const double dt, const int rk_stage) {}
  void initData();
  void initialHook() {}

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const { return 0; }
  int restoreData(double * buf) { return 0; }

  void force_bf(double (*f_bf)[9]) const {force_bf_wall(f_bf,this);}

  void addBoundaryFlux(IdealGasRhs* rhs) const;
  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>&args, const bool b_eval_func);
};

class Cbc : public IdealGasBc {
public:

  // cbc boundary conditions add the boundary flux by
  // computing a riemann flux associated with a constant
  // boundary state and the internal state.  this cbc
  // class is also pure virtual where the bf_state will
  // be populated in the register Hook

  IdealGasState* bf_state;

  Cbc(BfZone * p, IdealGasSolver* s): IdealGasBc(p,s), bf_state(NULL) {}
  ~Cbc() { DELETE(bf_state); }
  void calcRhs(const double time, const int rk_stage) {}
  void rk3Step(const double dt, const int rk_stage) {}
  void initialHook() {}

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  virtual int storeData(double * buf) const { return 0; }
  virtual int restoreData(double * buf) { return 0; }

  void force_bf(double (*f_bf)[9]) const {force_bf_cbc(f_bf,this);}
  void query(const string& param_str) const;

  void addBoundaryFlux(IdealGasRhs* rhs) const;
  void addScalarBoundaryFlux(double* rhs_sc) const;
};

class CbcProfile : public Cbc {
public:
  CbcProfile(BfZone* p, IdealGasSolver* s) : Cbc(p,s) {}
  void initData();
};

class CbcUpt : public Cbc {
public:
  CbcUpt(BfZone* p, IdealGasSolver* s) : Cbc(p,s) {}
  void initData();
};

class CbcMpt : public Cbc {
public:
  CbcMpt(BfZone* p, IdealGasSolver* s) : Cbc(p,s) {}
  void initData();
};

class CbcRunp : public Cbc {
public:
  CbcRunp(BfZone* p, IdealGasSolver* s) : Cbc(p,s) {}
  void initData();
};

class CbcRup : public Cbc {
public:
  bool b_aligned;
  double omega;
  double rot_axis[3];
  double u_bc[3];
  double (*u_fr)[3]; // the velocity field due to frame rotation ..

  CbcRup(BfZone* p, IdealGasSolver* s) : Cbc(p,s) {
    b_aligned = true;
    omega = HUGE_VAL;
    FOR_I3 rot_axis[i] = HUGE_VAL;
    FOR_I3 u_bc[i] = HUGE_VAL;
    u_fr = NULL;
  }
  ~CbcRup() {
    DELETE(u_fr);
  }
  void initData();
  void calcRhs(const double time, const int rk_stage);
  void calcBfState(const int ibf, const double time);
};

class CbcExtrapolate : public Cbc {
public:
  CbcExtrapolate(BfZone* p, IdealGasSolver* s) : Cbc(p,s) {}
  void initialHook();
  void rk3Step(const double dt, const int rk_stage);
  void calcBfState(const int ibf, const double unit_n[3]);
  void initData();
};

class CbcTotalPt : public Cbc {
public:

  double * un_bc;
  double un_bc_avg;
  double total_p;
  double total_t;
  double t_relax;
  double (*u_fr)[3]; // the velocity field due to frame rotation ..

  bool b_use_avg;

  CbcTotalPt(BfZone* p, IdealGasSolver* s) : Cbc(p,s) {
    un_bc = NULL; zone_ptr->registerBfData(un_bc,"un_bc",READWRITE_DATA);
    u_fr  = NULL;
    b_use_avg = false;
  }
  ~CbcTotalPt() {
    DELETE(un_bc);
    DELETE(u_fr);
  }
  void initData();
  void initialHook();
  void rk3Step(const double dt, const int rk_stage);
  void calcBfState(const int ibf, const double unit_n[3]);
};

class Nscbc : public IdealGasBc {
public:
  double * p_bc;
  double * rho_bc;
  double (*u_bc)[3];
  double (*rhs_bc)[3][5];

  // nscbc boundary conditions evolve some set of odes
  // at the boundary surface.  this class implements
  // basic update and registration commands, but a specific
  // nscbc class needs to implement how the rhs of the
  // boundary condition is computed..

  Nscbc(BfZone * p, IdealGasSolver* s): IdealGasBc(p,s) {
    p_bc     = NULL; zone_ptr->registerBfData(p_bc  ,"p_bc"  ,READWRITE_DATA);
    rho_bc   = NULL; zone_ptr->registerBfData(rho_bc,"rho_bc",READWRITE_DATA);
    u_bc     = NULL; zone_ptr->registerBfData(u_bc  ,"u_bc"  ,READWRITE_DATA);
    rhs_bc   = NULL;
  }
  ~Nscbc();

  void initData();
  void initialHook();

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const;
  int restoreData(double * buf);

  void force_bf(double (*f_bf)[9]) const {force_bf_cbc(f_bf,this);}
  void rk3Step(const double dt, const int rk_stage);

  void addBoundaryFlux(IdealGasRhs* rhs) const;
  void query(const string& param_str) const;

  void addScalarBoundaryFlux(double * rhs_sc) const;
};

class NscbcMt : public Nscbc {
public:
  double rhoun_bc;
  double T_bc;
  double (*u_fr)[3];

  NscbcMt(BfZone* p, IdealGasSolver* s) : Nscbc(p,s),rhoun_bc(0.0),T_bc(0.0),u_fr(NULL) {}
  ~NscbcMt() { DELETE(u_fr);}
  void initData();
  void calcRhs(const double time, const int rk_stage);
};

class NscbcOutletP : public Nscbc {
public:

  // outlet pressure vars ..
  double L_ref;
  double sigma;
  double p_ref;

  // outlet mdot vars
  double rhou_target;
  double T_relax;

  enum OutletType {
    OUTLET_PRESSURE,
    OUTLET_MDOT_WEAK,
    OUTLET_UNDEFINED
  };

  OutletType type;

  NscbcOutletP(BfZone* p, IdealGasSolver* s) : Nscbc(p,s), L_ref(0.0), sigma(0.0), p_ref(0.0),
                                               rhou_target(0.0), T_relax(0.0), type(OUTLET_UNDEFINED) {}
  void initData();
  void calcRhs(const double time, const int rk_stage);

  void calcRhsPressure(const double time, const int rk_stage);
  void calcRhsMdot(const double time, const int rk_stage);
};

class Sponge : public IdealGasBc {
public:

  double * T_sponge;
  double (*u_sponge)[3];

  double sponge_data[5];
  double eps_p_sponge; // not used any more
  double t_relax;
  double sponge_strength;
  double sponge_length;
  double sponge_pressure;

  SpongeType sponge_type;

  Sponge(BfZone* p, IdealGasSolver* s);

  ~Sponge();

  void addBoundaryFlux(IdealGasRhs* rhs) const;
  void calcRhs(const double time, const int rk_stage) {}
  void rk3Step(const double dt, const int rk_stage);

  void initData();
  void initialHook();

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const;
  int restoreData(double * buf);

  //void reportError(const string msg); // not used or defined
  void query(const string& param_str) const;

  void addScalarBoundaryFlux(double* rhs_sc) const;
  void force_bf(double (*f_bf)[9]) const {force_bf_cbc(f_bf,this);}

};



class InflowTurbulenceIG : public Cbc {
public:
  // inflow turbulence fields
  double (*um_it)[3];
  double (*Rd_it)[3];
  double (*Rod_it)[3];
  double* Tm_it;
  double* pm_it;

  // turbulence lengthscale
   double* length_scale;
   double max_length; // this is a maximum in plane lengt scale
   double in_plane_length; // this is a constant length used in the whole plane
   double Un; // this is a constant convective velocity
 //TODO add an out of plane lengthscale variable....or not

  double Umean_it, Ma_it; //these are used to scale the turbulent fluctuations
  double (*u_bc)[3];////TODO  This is not needed It only holds the fluctuating velocity field on the boundary....currently used for restarting smoothly and useful for debugging
  double *T_bc;//TODO This is not needed, and should be removed. Useful for debugging
  double *p_bc;//TODO This is not needed, and should be removed. Useful for debugging
  double (*avg)[3];
  double (*rms)[3];
  double IturbWgt;
  bool reset;

  // connectivity and equation solving related variables
  // use an underscore convention for counting the local active
  int* activeFaces;
  int* activeCvs;
  int ncv_,ncv_g_;
  int nfa_i_,nfa_;
  int (*cvofa_)[2];
  int *cvobf_;
  double * A_fa_;
  double * A_diag_cv_;
  double inflow_zero; // Threshold for convergence of inflow turbulence solver
  int inflow_maxiter;

  InflowTurbulenceIG(BfZone* p, IdealGasSolver* s) : Cbc(p,s) {
    um_it         = NULL; zone_ptr->registerBfData(um_it,         "um_it",       READWRITE_DATA);
    Rd_it         = NULL; zone_ptr->registerBfData(Rd_it,         "Rd_it",       READWRITE_DATA);
    Rod_it        = NULL; zone_ptr->registerBfData(Rod_it,        "Rod_it",      READWRITE_DATA);
    Tm_it         = NULL; zone_ptr->registerBfData(Tm_it,         "Tm_it",       READWRITE_DATA);
    pm_it         = NULL; zone_ptr->registerBfData(pm_it,         "pm_it",       READWRITE_DATA);
    length_scale  = NULL; zone_ptr->registerBfData(length_scale, "length_scale", READWRITE_DATA);
    u_bc          = NULL; zone_ptr->registerBfData(u_bc,          "u_bc",        READWRITE_DATA);
    T_bc          = NULL; zone_ptr->registerBfData(T_bc,          "T_bc",        READWRITE_DATA);
    p_bc          = NULL; zone_ptr->registerBfData(p_bc,          "p_bc",        READWRITE_DATA);
    avg           = NULL; zone_ptr->registerBfData(avg,           "avgIT",       READWRITE_DATA);
    rms           = NULL; zone_ptr->registerBfData(rms,           "rmsIT",       READWRITE_DATA);
    IturbWgt      = 0.0;  zone_ptr->registerBfData(IturbWgt,      "IturbWgt",    READWRITE_DATA);
    reset=false;

    activeFaces = NULL;
    activeCvs = NULL;

    //  persistent operators...
    cvofa_ = NULL;
    A_fa_ = NULL;
    A_diag_cv_ = NULL;
    cvobf_ = NULL;
  }



  ~InflowTurbulenceIG() {

    DELETE(Tm_it);
    DELETE(pm_it);
    DELETE(um_it);
    DELETE(Rd_it);
    DELETE(Rod_it);
    DELETE(length_scale);
    DELETE(u_bc);
    DELETE(T_bc);
    DELETE(p_bc);
    DELETE(avg);
    DELETE(rms);

    DELETE(cvofa_);
    DELETE(A_fa_);
    DELETE(A_diag_cv_);
    DELETE(cvobf_);
    DELETE(activeFaces);
    DELETE(activeCvs);

  }

  void rk3Step(const double dt, const int rk_stage){};
  void initData();
  int storeData(double * buf) const;
  int restoreData(double * buf);
  void initialHook();
  void calcRhs(const double time, const int rk_stage);
  //additional functions in InflowTurbulence.cpp
  void updateInflowTurbulence(double (*FacePtr)[3]);
  void setStatistics();
  void buildITBoundaryStructures();
  void setLengthscale();
  void buildFilterOperators();
  void initiateStats();
  void updateRandomField(double (*FacePtr)[3]);
  void removeMeanAndSetVarianceOne(double (*FacePtr)[3]);
  void computeCholesky(double (*FacePtr)[3]);
  void updateBC(double (*FacePtrNew)[3],double (*FacePtrOld)[3]);
  void query(const string& param_str) const;
  void updateCvData_(double (*var_)[3]);
  void updateCvData_(double *var_);
  int solveVecBoundaryCvCg(double (*sol)[3],const double * const A,double (*rhs)[3], double * AddedDiag,double zero, int maxiter );
  void getCvData_(double *FacePtr, double *CvPtr_);
  void getCvVectorData_(double (*FacePtr)[3], double (*CvPtr)[3]);

  virtual void inflowTurbulenceStatsHook(double (*um)[3],double (*Rd)[3], double (*Rod)[3], double *Tm, double *pm  ) {
    CERR(" > user must supply a hook for the statistics of the inflow turbulence on zone: " << getName());
  }
  virtual void inflowTurbulenceLengthscaleHook(double * length_scale) {
    CERR(" > user must supply a hook for the length scales of the inflow turbulence on zone: " << getName());
  }

};

class WmAlgCht : public IdealGasBc {
public:

  double * u1;
  double * tau_wall;
  double * q_wall;
  double * k_eff;

  // fem point data for this zone...
  int np;
  double * Tp;
  double * qp; // heat flow (flux times area)...
  int * ipobf_i;
  int * ipobf_v;
  double * ipobf_wgt;

  // max NLAYER is hard-coded as 3 for now...
  int nlayer;
  double * T_layer[3];
  double delta_layer[3];
  SimpleFunc * k_layer[3];

  WmAlgCht(BfZone* p, IdealGasSolver* s);
  ~WmAlgCht();

  void calcRhs(const double time, const int rk_stage);
  void rk3Step(const double dt, const int rk_stage) {}
  void initData();
  void initialHook();

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const;
  int restoreData(double * buf);

  void force_bf(double (*f_bf)[9]) const {force_bf_wm(f_bf,this);}

  void addBoundaryFlux(IdealGasRhs* rhs) const;
  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func);
};

class CbcProfileSimple : public Cbc { 
public: 
  CbcProfileSimple(BfZone* p, IdealGasSolver* s): Cbc(p,s) {}
  void initData();
}; 



#endif
