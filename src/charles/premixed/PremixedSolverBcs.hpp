#ifndef PREMIXEDSOLVERBCS_HPP
#define PREMIXEDSOLVERBCS_HPP

// ======================================================
// by convention, use undescores to separate
// subscript words, e.g. tau_wall, q_wall, y_plus, etc...
// ======================================================

#include "BcCommon.hpp"
#include "States.hpp"
#include "SimpleFunc.hpp"

class PremixedSolver;
typedef BcTemplate<PremixedSolver,PremixedRhs> PremixedBc;

class SimplePremixedBc : public PremixedBc {
public:

  SimplePremixedBc(BfZone* p, PremixedSolver* s): PremixedBc(p,s) {}

  void calcRhs(const double time, const int rk_stage) {}
  void rk3Step(const double dt, const int rk_stage) {}

  void initData() {
    assert( mf == NULL); mf = new double[zone_ptr->nbf];
  }

  void initialHook() {}

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability. Here we assume that the Simple bc does not store
  // or restore any data by default. The author of a special BC
  // can overide this behavior by writing their own...
  virtual int storeData(double * buf) const { return 0; }
  virtual int restoreData(double * buf) { return 0; }

  // note there will be no force information computed frm this bc..
  void force_bf(double (*f_bf)[9]) const {}

  virtual void addBoundaryFlux(PremixedRhs* rhs) const {
    assert(0); // user must provide an implementation here..
  }
};

class SlipWallPBc : public PremixedBc {
public:
  SlipWallPBc(BfZone* p, PremixedSolver* s) : PremixedBc(p,s) {}

  void calcRhs(const double time, const int rk_stage) {}
  void rk3Step(const double dt, const int rk_stage) {}
  void initData() {}
  void initialHook() {}

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const { return 0; }
  int restoreData(double * buf) { return 0; }

  void force_bf(double (*f_bf)[9]) const {force_bf_slip_wall(f_bf,this);}

  void addBoundaryFlux(PremixedRhs* rhs) const;
};

class MsTableP : public PremixedBc {
private:
  double mdot_bc,Z_bc,C_bc;
  double R_bc,T_bc,e_bc,gamma_bc,p_bc,rho_bc;
  double *rho_bf;
  double *mflux_bf;
public:
  MsTableP(Param *param, BfZone* p, PremixedSolver* s);
  virtual ~MsTableP();
  void calcRhs(const double time, const int rk_stage) {}
  void rk3Step(const double dt, const int rk_stage) {}
  void initData();
  void initialHook() {}

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const { return 0; }
  int restoreData(double * buf) { return 0; }

  void force_bf(double (*f_bf)[9]) const {} //force_bf_slip_wall(f_bf,this);}

  void addBoundaryFlux(PremixedRhs* rhs) const;
  void query(const string& param_str) const;
};

class WallAdiabaticP : public PremixedBc {
public:

  double u_bc[3];
  WallAdiabaticP(BfZone* p, PremixedSolver* s) : PremixedBc(p,s) {
    for (int i = 0; i < 3; ++i) u_bc[i] = 0.0;

    // add to funcEvalList
    funcEvalList.push_back(p->getName()+":tau_wall()");
    funcEvalList.push_back(p->getName()+":y_plus()");
  }

  void calcRhs(const double time, const int rk_stage) {}
  void rk3Step(const double dt, const int rk_stage) {}
  void initData() {}
  void initialHook() {}

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const { return 0; }
  int restoreData(double * buf) { return 0; }

  void force_bf(double (*f_bf)[9]) const {force_bf_wall(f_bf,this);}

  void addBoundaryFlux(PremixedRhs* rhs) const;
  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args,const bool b_eval_func);

};

class WallIsothermalP : public PremixedBc {
public:

  double u_bc[3];
  double T_bc;
  WallIsothermalP(BfZone* p, PremixedSolver* s) : PremixedBc(p,s) {
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

  void addBoundaryFlux(PremixedRhs* rhs) const;
  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func);

};

class WallConductiveP : public PremixedBc {
public:

  double u_bc[3];
  double T_end_bc;
  vector<double> thermal_conductivity_bc;
  vector<double> layer_length_bc;
  double thermal_resistance_sum;

  WallConductiveP(BfZone* p, PremixedSolver* s) : PremixedBc(p,s) {
    T_end_bc = 0.0;
    thermal_resistance_sum = 0.0;
    for (int i =0; i < 3; ++i) u_bc[i] = 0.0;

    funcEvalList.push_back(p->getName()+":tau_wall()");
    funcEvalList.push_back(p->getName()+":y_plus()");
    funcEvalList.push_back(p->getName()+":q_wall()");
    //funcEvalList.push_back(p->getName()+":T1()");
    //funcEvalList.push_back(p->getName()+":T2()");
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

  void addBoundaryFlux(PremixedRhs* rhs) const;
  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args,const bool b_eval_func);
};

class WallChtP : public PremixedBc {
private:
  //fem data point data for this zone...
  int np;
  double * Tp;
  double * qp;
  int * ipobf_i;
  int * ipobf_v;
  double * ipobf_wgt;

public:

  WallChtP(BfZone* p, PremixedSolver * s) : PremixedBc(p,s) {
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
    funcEvalList.push_back(p->getName()+":T_bc()");
    funcEvalList.push_back(p->getName()+":T_wall()");
  }

  ~WallChtP() {
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

  void force_bf(double (*f_bf)[9]) const {}
  void addBoundaryFlux(NonpremixedRhs* rhs) const;
  void query(const string& param_str) const;
  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func);
};


class WmAlgAdiabaticP : public PremixedBc {
public:

  double * u1;
  double * tau_wall;
  double u_bc[3];

  WmAlgAdiabaticP(BfZone* p, PremixedSolver* s) : PremixedBc(p,s) {

    // XXX note that this bc currently assumes a no-slip (non-moving wall)
    for (int i =0; i < 3; ++i)
      u_bc[i] = 0.0;

    // tau_wall should be NOREADWRITE_DATA
    u1    = NULL; zone_ptr->registerBfData(u1,  "u1", READWRITE_DATA);
    tau_wall = NULL; zone_ptr->registerBfData(tau_wall, "tau_wall", READWRITE_DATA);

    // set the load balance cost. Default for bc is 1
    // can't do this here with current solver because this is called
    // AFTER lb routine...
    /*
    if (mpi_rank == 0) cout << "XXXXX setting lb_cost to 5" << endl;
    assert(zone_ptr->lb_cost == 1);
    zone_ptr->lb_cost = 5;
    assert(p->lb_cost == 5);
    */

    // add to funcEvalList
    funcEvalList.push_back(p->getName()+":tau_wall()");
    funcEvalList.push_back(p->getName()+":y_plus()");
  }

  ~WmAlgAdiabaticP();

  void calcRhs(const double time, const int rk_stage);
  void rk3Step(const double dt, const int rk_stage) {}
  void initData();
  void initialHook();

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const;
  int restoreData(double * buf);

  void force_bf(double (*f_bf)[9]) const {force_bf_wm(f_bf,this);}

  void addBoundaryFlux(PremixedRhs* rhs) const;
  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func);
};

class WmAlgIsothermalP : public PremixedBc {
public:

  double * u1;
  double * tau_wall;
  double * q_wall;
  double * T_wall;
  //double u_bc[3];
  //double T_bc; // now T_bc is treated as a variable over bf's even if constant

  WmAlgIsothermalP(BfZone* p, PremixedSolver* s) : PremixedBc(p,s) {

    // XXX note that this bc currently assumes a no-slip (non-moving wall)
    //for (int i =0; i < 3; ++i)
    //  u_bc[i] = 0.0;

    //T_bc  = 0.0;

    // tau_wall should be NOREADWRITE_DATA
    u1     = NULL; zone_ptr->registerBfData(u1,"u1",READWRITE_DATA);
    tau_wall  = NULL; zone_ptr->registerBfData(tau_wall,"tau_wall",READWRITE_DATA);
    q_wall = NULL; zone_ptr->registerBfData(q_wall,"q_wall",READWRITE_DATA);
    T_wall = NULL; zone_ptr->registerBfData(T_wall,"T_wall",READWRITE_DATA);

    // add to funcEvalList
    funcEvalList.push_back(p->getName()+":tau_wall()");
    funcEvalList.push_back(p->getName()+":y_plus()");
    funcEvalList.push_back(p->getName()+":q_wall()");

    // set the load balance cost. Default for bc is 1
    // see note above
    /*
      if (mpi_rank == 0) cout << "XXXXX setting lb_cost to 5" << endl;
      assert(zone_ptr->lb_cost == 1);
      zone_ptr->lb_cost = 5;
      assert(p->lb_cost == 5);
    */

  }

  ~WmAlgIsothermalP();

  void calcRhs(const double time, const int rk_stage);
  void rk3Step(const double dt, const int rk_stage) {}
  void initData();
  void initialHook();

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const;
  int restoreData(double * buf);

  void force_bf(double (*f_bf)[9]) const {force_bf_wm(f_bf,this);}

  void addBoundaryFlux(PremixedRhs* rhs) const;
  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func);
};

class WmAlgChtP : public PremixedBc {
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

  WmAlgChtP(BfZone* p, PremixedSolver* s);
  ~WmAlgChtP();

  void calcRhs(const double time, const int rk_stage);
  void rk3Step(const double dt, const int rk_stage) {}
  void initData();
  void initialHook();

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const;
  int restoreData(double * buf);

  void force_bf(double (*f_bf)[9]) const {force_bf_wm(f_bf,this);}

  void addBoundaryFlux(PremixedRhs* rhs) const;
  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func);
};

class WmAlgConductiveP : public PremixedBc {
public:
  double * u1;
  double * tau_wall;
  double * q_wall;
  double * T_bc;
  double T_end_bc;
  vector<double> layer_length_bc;
  vector<double> thermal_conductivity_bc;
  double thermal_resistance_external_times_area;

  WmAlgConductiveP(BfZone * p, PremixedSolver* s) : PremixedBc(p,s) {
    u1       = NULL; zone_ptr->registerBfData(u1,       "u1",       READWRITE_DATA);
    tau_wall = NULL; zone_ptr->registerBfData(tau_wall, "tau_wall", READWRITE_DATA);
    q_wall   = NULL; zone_ptr->registerBfData(q_wall,   "q_wall",   READWRITE_DATA);
    T_bc     = NULL; zone_ptr->registerBfData(T_bc,     "T_bc",     READWRITE_DATA);

    // add to funcEvalList
    funcEvalList.push_back(p->getName()+":tau_wall()");
    funcEvalList.push_back(p->getName()+":y_plus()");
    funcEvalList.push_back(p->getName()+":q_wall()");
    funcEvalList.push_back(p->getName()+":T_bc()");

    // Make sure there's a problem if these values don't get set...
    thermal_resistance_external_times_area = 0.0;
  }

  ~WmAlgConductiveP();

  void calcRhs(const double time, const int rk_stage);
  void rk3Step(const double dt, const int rk_stage) {}
  void initData();
  void initialHook();

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const;
  int restoreData(double * buf);

  void force_bf(double (*f_bf)[9]) const {force_bf_wm(f_bf, this);}

  void addBoundaryFlux(NonpremixedRhs * rhs) const;
  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func);
};


class CbcP : public PremixedBc {
public:

  // cbc boundary conditions add the boundary flux by
  // computing a riemann flux associated with a constant
  // boundary state and the internal state.  this cbc
  // class is also pure virtual where the bf_state will
  // be populated in the register Hook
  PremixedState* bf_state;

  CbcP(BfZone * p, PremixedSolver* s): PremixedBc(p,s), bf_state(NULL) {}
  ~CbcP() { DELETE(bf_state); }
  void calcRhs(const double time, const int rk_stage) {}
  void rk3Step(const double dt, const int rk_stage) {}
  void initialHook() {}

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  virtual int storeData(double * buf) const { return 0; }
  virtual int restoreData(double * buf) { return 0; }

  void force_bf(double (*f_bf)[9]) const {force_bf_cbc(f_bf,this);}
  void query(const string& param_str) const;

  void addBoundaryFlux(PremixedRhs* rhs) const;
  void addScalarBoundaryFlux(double * rhs_sc) const;
};

class CbcProfileP : public CbcP {
public:
  CbcProfileP(BfZone* p, PremixedSolver* s) : CbcP(p,s) {}
  void initData();
};

class CbcUptsP : public CbcP {
public:
  CbcUptsP(BfZone* p, PremixedSolver* s) : CbcP(p,s) {}
  void initData();
};

class CbcUptsSpongeP : public CbcP {
public:
  int idir;
  double rho_sponge,u_sponge[3],rhoE_sponge,Z_sponge,C_sponge;
  double x0,x1,sponge_strength;
  PremixedState bf_state_sponge;
  CbcUptsSpongeP(BfZone* p, PremixedSolver* s) : CbcP(p,s) {}
  void initData();
  void addBoundaryFlux(PremixedRhs* rhs) const;
};

class CbcMptsP : public CbcP {
public:
  CbcMptsP(BfZone* p, PremixedSolver* s) : CbcP(p,s) {}
  void initData();
};

class CbcMsTableP : public CbcP {
public:
  CbcMsTableP(BfZone* p, PremixedSolver* s) : CbcP(p,s) {}
  void initData();
};

class CbcRunpsP : public CbcP {
public:
  CbcRunpsP(BfZone* p, PremixedSolver* s) : CbcP(p,s) {}
  void initData();
};

class CbcRupsP : public CbcP {
public:
  CbcRupsP(BfZone* p, PremixedSolver* s) : CbcP(p,s) {}
  void initData();
};

class CbcTotalPtsP : public CbcP {
public:
  double * un_bc;
  double p_tot, T_tot, Z_bc, C_bc, t_relax, R_0, gamma_0, a_gam_0, e_0, h_0, e_tot, h_tot, T_0, p_star;

  CbcTotalPtsP(BfZone* p, PremixedSolver *s) : CbcP(p,s) {
    un_bc = NULL; zone_ptr->registerBfData(un_bc, "un_bc", READWRITE_DATA);
  }

  ~CbcTotalPtsP() {
    DELETE(un_bc);
  }

  void initData();
  void initialHook();
  void rk3Step(const double dt, const int rk_stage);
  void calcBfState(const int ibf, const double unit_n[3]);
};

class NscbcP : public PremixedBc {
public:
  double * p_bc;
  double * rho_bc;
  double (*u_bc)[3];
  double * Z_bc;
  double * C_bc;
  double (*rhs_bc)[3][7];

  double * R_bc;
  double * T0;
  double * e0;
  double * gamma0;
  double * agam;

  // nscbc boundary conditions evolve some set of odes
  // at the boundary surface.  this class implements
  // basic update and registration commands, but a specific
  // nscbc class needs to implement how the rhs of the
  // boundary condition is computed..

  NscbcP(BfZone * p, PremixedSolver* s): PremixedBc(p,s) {
    p_bc     = NULL; zone_ptr->registerBfData(p_bc  ,"p_bc"  ,READWRITE_DATA);
    rho_bc   = NULL; zone_ptr->registerBfData(rho_bc,"rho_bc",READWRITE_DATA);
    u_bc     = NULL; zone_ptr->registerBfData(u_bc  ,"u_bc"  ,READWRITE_DATA);
    Z_bc     = NULL; zone_ptr->registerBfData(Z_bc  ,"Z_bc"  ,READWRITE_DATA);
    C_bc     = NULL; zone_ptr->registerBfData(C_bc  ,"C_bc"  ,READWRITE_DATA);
    rhs_bc   = NULL;

    R_bc     = NULL;
    T0       = NULL;
    e0       = NULL;
    gamma0   = NULL;
    agam     = NULL;
  }
  ~NscbcP();

  void initData();
  void initialHook();

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const;
  int restoreData(double * buf);

  void force_bf(double (*f_bf)[9]) const {force_bf_cbc(f_bf,this);}
  void rk3Step(const double dt, const int rk_stage);

  void addBoundaryFlux(PremixedRhs* rhs) const;
  void query(const string& param_str) const;

  // different nscbc types have different requirements for what
  // data (eg gas properties: R, mu, etc) need to be re-computed
  // prior to the boundaryFlux being computed.  currently needs
  // to be called before the addBoundaryFlux-- by default it is empty.
  // called inside of calcRhs

  virtual void updatePrimitiveData() {}

  void addScalarBoundaryFlux(double* rhs_sc) const;
};

class NscbcMtsP : public NscbcP {
private:
  bool eval_mdot_bc;
  string mdot_bc_str;
  bool eval_T_bc;
  string T_bc_str;
  bool eval_Z_bc;
  string Z_bc_str;
  bool eval_C_bc;
  string C_bc_str;
public:
  double rhoun_bc;
  double T_bc;
  double Z_bc_,C_bc_;
  NscbcMtsP(BfZone* p, PremixedSolver* s) :
    NscbcP(p,s),eval_mdot_bc(false),eval_T_bc(false),eval_Z_bc(false),eval_C_bc(false),
    rhoun_bc(0.0),T_bc(0.0),Z_bc_(0.0),C_bc_(0.0) {}
  void initData();
  void initialHook();
  void calcRhs(const double time, const int rk_stage);
  void evalZCStrAndUpdate();
};

class NscbcMtsSoftP : public NscbcP {
public:
  double rhoun_bc;
  double T_bc;
  double Z_bc_,C_bc_;
  double t_relax;

  NscbcMtsSoftP(BfZone* p, PremixedSolver* s) :
     NscbcP(p,s),rhoun_bc(0.0),T_bc(0.0),Z_bc_(0.0),C_bc_(0.0) {}
  void initData();
  void initialHook();
  void calcRhs(const double time, const int rk_stage);
};

class NscbcOutletPP : public NscbcP {
public:
  double L_ref;
  double sigma;
  double p_ref;

  NscbcOutletPP(BfZone* p, PremixedSolver* s) : NscbcP(p,s), L_ref(0.0), sigma(0.0), p_ref(0.0) {}
  void initData();
  void calcRhs(const double time, const int rk_stage);
  void updatePrimitiveData();
  void addBoundaryFlux(PremixedRhs* rhs) const;
};

class NscbcOutletMdotP : public NscbcP {
public:
  double rhoun_bc;

  NscbcOutletMdotP(BfZone* p, PremixedSolver* s) : NscbcP(p,s),rhoun_bc(0.0) {}
  void initData();
  void calcRhs(const double time, const int rk_stage);
  void updatePrimitiveData();
};


class SpongeP : public PremixedBc {
public:

  double * rho_sponge;
  double (*rhou_sponge)[3];
  double *rhoE_sponge;
  double * rhoZ_sponge;
  double * rhoC_sponge;

  double sponge_data[5];
  double eps_p_sponge;
  double t_relax;
  double sponge_strength;
  double sponge_length;
  double sponge_pressure;

  SpongeType sponge_type;

  SpongeP(BfZone* p, PremixedSolver* s) : PremixedBc(p,s) {
    for (int i =0; i < 5; ++i)
      sponge_data[i] = 0.0;
    eps_p_sponge     = 0.0;
    t_relax          = 0.0;
    sponge_strength  = 0.0;
    sponge_length    = 0.0;
    sponge_type      = SPONGE_UNDEFINED;

    rho_sponge  = NULL;
    rhou_sponge = NULL;
    rhoE_sponge = NULL;
    rhoZ_sponge = NULL;
    rhoC_sponge = NULL;
  }

  ~SpongeP();
  void addBoundaryFlux(PremixedRhs* rhs) const;
  void calcRhs(const double time, const int rk_stage) {}
  void rk3Step(const double dt, const int rk_stage);

  void initData();
  void initialHook();

  // store/restore is for robustly addressing NaN's due to explicit
  // numerical stability.
  int storeData(double * buf) const;
  int restoreData(double * buf);

  void force_bf(double (*f_bf)[9]) const {force_bf_cbc(f_bf,this);}

  //void reportError(const string msg); // not used or defined?
  void query(const string& param_str) const;

  void addScalarBoundaryFlux(double* rhs_sc) const;
};

class InflowTurbulenceP : public CbcP {
public:
  // inflow turbulence fields
  double (*um_it)[3];
  double (*Rd_it)[3];
  double (*Rod_it)[3];
  double* Tm_it;
  double* pm_it;
  double* Zm_it;
  double* Cm_it;

  // turbulence lengthscale
  double* length_scale;
  double max_length; // this is a maximum in plane length scale
  double in_plane_length; // this is a constant length used in the whole plane
  double Un; // this is a constant convective velocity
  // TODO: add an out of plane lengthscale variable... or not

  double Umean_it, Ma_it; // these are used to scale the turbulent fluctuations
  double (*u_bc)[3];//TODO: This is not needed.
  //                        It only holds the fluctuating velocity field on the boundary.
  //                        Currently used for restarting smoothly and useful for debugging.
  double *T_bc;//TODO This is not needed and should be removed. Useful for debugging
  double *p_bc;//TODO This is not needed and should be removed. Useful for debugging
  // double *Z_bc;
  // double *C_bc;
  double (*avg)[3];
  double (*rms)[3];
  double IturbWgt;
  bool reset;

  // connectivity and equation solving related variables
  // use an underscore convention for counting the local active
  int* activeFaces;
  int* activeCvs;
  int ncv_,ncv_g_;
# define FOR_ICV_ for (int icv_ = 0; icv_ < ncv_; ++icv_)
# define FOR_IFA_ for (int ifa_ = 0; ifa_ < nfa_; ++ifa_)
  int nfa_i_,nfa_;
  int (*cvofa_)[2];
  int *cvobf_;
  double * A_fa_;
  double *A_diag_cv_;
  double inflow_zero;//Threshold for convergence of inflow turbulence solver
  int inflow_maxiter;


  InflowTurbulenceP(BfZone* p, PremixedSolver* s) : CbcP(p,s) {
    um_it         = NULL; zone_ptr->registerBfData(um_it,         "um_it",       READWRITE_DATA);
    Rd_it         = NULL; zone_ptr->registerBfData(Rd_it,         "Rd_it",       READWRITE_DATA);
    Rod_it        = NULL; zone_ptr->registerBfData(Rod_it,        "Rod_it",      READWRITE_DATA);
    Zm_it         = NULL; zone_ptr->registerBfData(Zm_it,         "Zm_it",       READWRITE_DATA);
    Cm_it         = NULL; zone_ptr->registerBfData(Cm_it,         "Cm_it",       READWRITE_DATA);
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

    // persistent operators...
    cvofa_ = NULL;
    A_fa_ = NULL;
    A_diag_cv_ = NULL;
    cvobf_ = NULL;

  }//InflowTurbulenceP()

  ~InflowTurbulenceP() {
    //  variables
    DELETE(um_it);
    DELETE(Rd_it);
    DELETE(Rod_it);
    DELETE(Tm_it);
    DELETE(pm_it);
    DELETE(Zm_it);
    DELETE(Cm_it);
    DELETE(length_scale);
    DELETE(u_bc);
    DELETE(avg);
    DELETE(rms);
    //  and operators
    DELETE(cvofa_);
    DELETE(A_fa_);
    DELETE(A_diag_cv_);
    DELETE(cvobf_);
  }//~InflowTurbulenceP()

  void initData();
  void initialHook();
  int storeData(double * buf) const;
  int restoreData(double * buf);
  void calcRhs(const double time, const int rk_stage);
  // Additional functions in InflowTurbulence.cpp
  void updateInflowTurbulence(double (*FacePtr)[3]);
  void setStatistics();
  void buildITBoundaryStructures();
  void setLengthscale();
  void buildFilterOperators();
  void initiateStats();
  void updateRandomField(double (*FacePtr)[3]);
  void removeMeanAndSetVarianceOne(double (*FacePtr)[3]);
  void computeCholesky(double (*FacePtr)[3]);
  void updateBC(double (*FacePtrNew)[3], double (*FacePtrOld)[3]);
  void updateCvData_(double (*var_)[3]);
  void updateCvData_(double *var_);
  int solveVecBoundaryCvCg(double (*sol)[3], const double * const A, double (*rhs)[3], double * AddedDiag, double zero, int maxiter);
  void getCvData_(double *FacePtr, double *CvPtr_);
  void getCvVectorData_(double (*FacePtr)[3], double (*CvPtr)[3]);

  virtual void inflowTurbulenceStatsHook(double (*um)[3], double (*Rd)[3], double (*Rod)[3], double *Tm, double *pm) {
    CERR(" > user must supply a hook for the statistics of the inflow turbulence on zone: " << getName());
  }//inflowTurbulenceStatsHook()

  virtual void inflowTurbulenceLengthscaleHook(double * length_scale) {
    CERR(" > user must supply a hook for the length scales of the inflow turbulence on zone: " << getName());
  }
};

// convenience function to help set the boundary state values

inline double calcInternalEnergyFromTS(const double& T_, const double& T0, const double& e0,
                                       const double& R_gas, const double& gamma0, const double& a_gam) {

  const double gamma = fmax(gamma0 + a_gam*(T_-T0), 1.0001);
  return e0 + R_gas/a_gam * log((gamma-1.0)/(gamma0-1.0)); // invert T-T0 for e-e0 (see updatePrimitiveData)

}
#endif
