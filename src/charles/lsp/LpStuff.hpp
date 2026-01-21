#ifndef __LPTRACER_HPP__
#define __LPTRACER_HPP__

#include "CTI.hpp"
using namespace CTI;

#include "MiscUtils.hpp"
using namespace MiscUtils;

#include "LpContainer.hpp"

class LpBcBase;

// ======================================================
// LpTracer: Lagrangian Particle Tracers
// ======================================================

class LpTracerRhs {
public:

  double xp[3];

  inline void zero(){
    FOR_I3 xp[i] = 0.0;
  };

  // rhs += val*rhs2
  inline void axpy(const double& val, const LpTracerRhs& rhs) {
    for (unsigned int i =0; i <3; ++i)
      xp[i] += val*rhs.xp[i];
  }

};

class LpTracerState {
public:
  
  // flag:  contains the information for material_id, keep status and injector id
  //        material_id: the first char (byte) of the int -> range: 0:127 (neglecting negative ones)
  //        keep status: the second char (byte) of the int -> range: 0:127 (onlye two values are needed for now)
  //        injector id: the last two bytes of the int -> range: 0:32767 (neglecting negative ones)

  // primary data...
  int icv,flag;
  double xp[3];
  double xp0[3];
  double up[3];
  double dp;
  
  static int data_size() { return 11; }  // primary data size
  
  inline void pack(double * buf) const { 
    int * icv_flag = (int*)buf; // use the first 8 byte double to store 2 4 byte ints
    icv_flag[0] = icv;
    icv_flag[1] = flag;
    buf[1] = xp[0];
    buf[2] = xp[1];
    buf[3] = xp[2];
    buf[4] = xp0[0];
    buf[5] = xp0[1];
    buf[6] = xp0[2];
    buf[7] = up[0];
    buf[8] = up[1];
    buf[9] = up[2];
    buf[10] = dp;
  }
  
  inline void packT(double * buf,const double t[3]) const {
    int * icv_flag = (int*)buf; // use the first 8 byte double to store 2 4 byte ints
    icv_flag[0] = icv;
    icv_flag[1] = flag;
    // note: subtract t here...
    buf[1] = xp[0]-t[0];
    buf[2] = xp[1]-t[1];
    buf[3] = xp[2]-t[2];
    buf[4] = xp0[0]-t[0];
    buf[5] = xp0[1]-t[1];
    buf[6] = xp0[2]-t[2];
    buf[7] = up[0];
    buf[8] = up[1];
    buf[9] = up[2];
    buf[10] = dp;
  }

  inline void packR(double * buf,const double R[9]) const {
    int * icv_flag = (int*)buf; // use the first 8 byte double to store 2 4 byte ints
    icv_flag[0] = icv;
    icv_flag[1] = flag;
    // note: apply R here...
    MiscUtils::applyInvRotation(buf+1,R,xp);
    MiscUtils::applyInvRotation(buf+4,R,xp0);
    MiscUtils::applyInvRotation(buf+7,R,up);
    buf[10] = dp;
  }

  inline void packRT(double * buf,const double R[9],const double t[3]) const {
    int * icv_flag = (int*)buf; // use the first 8 byte double to store 2 4 byte ints
    icv_flag[0] = icv;
    icv_flag[1] = flag;
    // note: apply R here...
    MiscUtils::applyInvRotation(buf+1,R,xp);
    MiscUtils::applyInvRotation(buf+4,R,xp0);
    MiscUtils::applyInvRotation(buf+7,R,up);
    // note: subtract t here...
    buf[1] -= t[0];
    buf[2] -= t[1];
    buf[3] -= t[2];
    buf[4] -= t[0];
    buf[5] -= t[1];
    buf[6] -= t[2];
    buf[10] = dp;
  }

  inline void unpack(double * buf) {
    int * icv_flag = (int*)buf;
    icv  = icv_flag[0];
    flag = icv_flag[1];
    xp[0] = buf[1];
    xp[1] = buf[2];
    xp[2] = buf[3];
    xp0[0] = buf[4];
    xp0[1] = buf[5];
    xp0[2] = buf[6];
    up[0] = buf[7];
    up[1] = buf[8];
    up[2] = buf[9];
    dp = buf[10];
  }
  
  inline void add(const LpTracerRhs& other) {
    xp[0] += other.xp[0];
    xp[1] += other.xp[1];
    xp[2] += other.xp[2];
  }

};

class LpSolidRhs {
public:

  double xp[3];
  double up[3];
  double Tp;

  inline void zero(){
    FOR_I3 xp[i] = 0.0;
    FOR_I3 up[i] = 0.0;
    Tp = 0.0;
  };

  // rhs += val*rhs2
  inline void axpy(const double& val, const LpSolidRhs& rhs) {
    for (unsigned int i =0; i <3; ++i)
      xp[i] += val*rhs.xp[i];
    for (unsigned int i =0; i <3; ++i)
      up[i] += val*rhs.up[i];
    Tp += val*rhs.Tp;
  }

};

class LpSolidState : public LpTracerState {
public:
  //int icv;
  //int flag
  //double xp[3];
  //double xp0[3];
  //double up[3];
  //double dp;
  double Tp;
  double mp;
  double npar; // number of particles in each parcel

  // auxillary variable for implicit time integration
  // do not need interprocessor trasport or be copied 
  double k;   // equals f1/tau in the momentum eqn
  double kt;
  // boundary coefficients
  double kb;
  double ktb;

  static int data_size() { return 14; }

  inline void pack(double * buf) const {
    LpTracerState::pack(buf);
    buf[11] = Tp;
    buf[12] = mp;
    buf[13] = npar;
  }

  inline void packRT(double * buf,const double R[9], const double t[3]) const {
    LpTracerState::packRT(buf, R, t);
    buf[11] = Tp;
    buf[12] = mp;
    buf[13] = npar;
  }

  inline void packR(double * buf,const double R[9]) const {
    LpTracerState::packR(buf, R);
    buf[11] = Tp;
    buf[12] = mp;
    buf[13] = npar;
  }

  inline void packT(double * buf, const double t[3]) const {
    LpTracerState::packT(buf, t);
    buf[11] = Tp;
    buf[12] = mp;
    buf[13] = npar;
  }

  inline void unpack(double * buf) {
    LpTracerState::unpack(buf);
    Tp = buf[11];
    mp = buf[12];
    npar = buf[13];
  }

  inline void add(const LpSolidRhs& other) {
    xp[0] += other.xp[0];
    xp[1] += other.xp[1];
    xp[2] += other.xp[2];
    up[0] += other.up[0];
    up[1] += other.up[1];
    up[2] += other.up[2];
    Tp    += other.Tp;
  }

  inline void add_im(const LpSolidRhs& rhs, const LpSolidRhs& rhs_d, const double& dt) {
    for (int i =0; i < 3; ++i) {
      xp[i] = (xp[i] + rhs.xp[i])/(1.0 - rhs_d.xp[i]*dt);
      up[i] = (up[i] + rhs.up[i])/(1.0 - rhs_d.up[i]*dt);
    }
    Tp = (Tp + rhs.Tp)/(1.0 - rhs_d.Tp*dt);
  }

  /*
  inline void add_im_xm(const LspRhs& rhs, const LspRhs& rhs_d, const double& dt) {
    for (int i =0; i < 3; ++i) {
      xp[i] = (xp[i] + rhs.xp[i])/(1.0 - rhs_d.xp[i]*dt);
      //up[i] = (up[i] + rhs.up[i])/(1.0 - rhs_d.up[i]*dt);
    }
    mp = (mp + rhs.mp)/(1.0 - rhs_d.mp*dt);
    //Tp = (Tp + rhs.Tp)/(1.0 - rhs_d.Tp*dt);
  }
  */
  
};

class LpLiquidRhs {
public:

  double xp[3];
  double up[3];
  double mp;
  double Tp;

  inline void zero(){
    FOR_I3 xp[i] = 0.0;
    FOR_I3 up[i] = 0.0;
    mp = 0.0;
    Tp = 0.0;
  };

  // rhs += val*rhs2
  inline void axpy(const double& val, const LpLiquidRhs& rhs) {
    for (unsigned int i =0; i <3; ++i)
      xp[i] += val*rhs.xp[i];
    for (unsigned int i =0; i <3; ++i)
      up[i] += val*rhs.up[i];
    mp += val*rhs.mp;
    Tp += val*rhs.Tp;
  }

};

class LpLiquidState : public LpTracerState {
public:
  //int icv;
  //int flag
  //double xp[3];
  //double xp0[3];
  //double up[3];
  //double dp;
  double Tp;
  double mp;
  double npar;  // number of particles in each parcel
  double tbu;   // breakup time

  // just for inspection
  double mp0;

  // auxillary variable for implicit time integration
  // do not need interprocessor trasport or be copied 
  double k;   // equals f1/tau in the momentum eqn
  double kt;
  // boundary coefficients
  double kb;
  double ktb;

  static int data_size() { return 15; }

  inline void pack(double * buf) const {
    LpTracerState::pack(buf);
    buf[11] = Tp;
    buf[12] = mp;
    buf[13] = tbu;
    buf[14] = npar;
  }

  inline void packRT(double * buf,const double R[9], const double t[3]) const {
    LpTracerState::packRT(buf, R, t);
    buf[11] = Tp;
    buf[12] = mp;
    buf[13] = tbu;
    buf[14] = npar;
  }

  inline void packR(double * buf,const double R[9]) const {
    LpTracerState::packR(buf, R);
    buf[11] = Tp;
    buf[12] = mp;
    buf[13] = tbu;
    buf[14] = npar;
  }

  inline void packT(double * buf, const double t[3]) const {
    LpTracerState::packT(buf, t);
    buf[11] = Tp;
    buf[12] = mp;
    buf[13] = tbu;
    buf[14] = npar;
  }

  inline void unpack(double * buf) {
    LpTracerState::unpack(buf);
    Tp = buf[11];
    mp = buf[12];
    tbu = buf[13];
    npar = buf[14];
  }

  inline void add(const LpLiquidRhs& other) {
    xp[0] += other.xp[0];
    xp[1] += other.xp[1];
    xp[2] += other.xp[2];
    up[0] += other.up[0];
    up[1] += other.up[1];
    up[2] += other.up[2];
    Tp    += other.Tp;
    mp    += other.mp;
  }

  inline void add_im(const LpLiquidRhs& rhs, const LpLiquidRhs& rhs_d, const double& dt) {
    for (int i =0; i < 3; ++i) {
      xp[i] = (xp[i] + rhs.xp[i])/(1.0 - rhs_d.xp[i]*dt);
      up[i] = (up[i] + rhs.up[i])/(1.0 - rhs_d.up[i]*dt);
    }
    mp = (mp + rhs.mp)/(1.0 - rhs_d.mp*dt);
    Tp = (Tp + rhs.Tp)/(1.0 - rhs_d.Tp*dt);
  }

  /*
  inline void add_im_xm(const LspRhs& rhs, const LspRhs& rhs_d, const double& dt) {
    for (int i =0; i < 3; ++i) {
      xp[i] = (xp[i] + rhs.xp[i])/(1.0 - rhs_d.xp[i]*dt);
      //up[i] = (up[i] + rhs.up[i])/(1.0 - rhs_d.up[i]*dt);
    }
    mp = (mp + rhs.mp)/(1.0 - rhs_d.mp*dt);
    //Tp = (Tp + rhs.Tp)/(1.0 - rhs_d.Tp*dt);
  }
  */
  
};

class LpTracer : public LpContainer<LpTracerState> {
private:
public:
  double (*u)[3];
  int ncv;
  vector<LpBcBase*> bcVec;

  // consistant RK3 scheme according to Gottleib, Shu, Tadmor equation 4.2
  // with this formulation we don't need to store rhs for previous stages
  double rk_wgt[3][2];  // = {  {0., 1.},
                        //      {3./4., 1./4.},
                        //      {1./3., 2./3.} };

  LpTracer() : LpContainer<LpTracerState>() {
    u = NULL;
    ncv = -1;
    rk_wgt[0][0] = 0.;
    rk_wgt[0][1] = 1.;
    rk_wgt[1][0] = 3./4.;
    rk_wgt[1][1] = 1./4.;
    rk_wgt[2][0] = 1./3.;
    rk_wgt[2][1] = 2./3.;
  }
  
  virtual ~LpTracer() {}

  // should call setLp0() before this in FlowSolver
  void rk3Step(const double dt, const int stage) {
    for (int ip = 0; ip < np; ++ip) {
      FOR_I3 lp[ip].xp[i] = rk_wgt[stage-1][0]*lp[ip].xp0[i] + rk_wgt[stage-1][1]*(lp[ip].xp[i] + dt*u[lp[ip].icv][i]);
    }
  }

  void report() const {

    const int lpSize  = getLpSize();
  
    if (mpi_rank==0) cout << " > report lpTracer:" << endl;

    cti_real* r1 = new cti_real[this->size()];
    cti_real (*r2)[3] = new cti_real[this->size()][3];
    for (int ip = 0; ip < this->size(); ++ip) {
      FOR_I3 r2[ip][i] = this->lp[ip].xp[i];
    }
    dumpRange(r2,this->size(),"xp");
    for (int ip = 0; ip < this->size(); ++ip) {
      r1[ip] = this->lp[ip].dp;
    }
    dumpRange(r1,this->size(),"dp");

    IF_RANK0 cout << " > lpTracer size: " << lpSize << endl;
  
    delete[] r1;
    delete[] r2;

  }

};
  
class LpSolid : public LpContainer<LpSolidState> {
private:
public:
  double (*u)[3];
  int ncv;
  vector<LpBcBase*> bcVec;

  LpSolid() : LpContainer<LpSolidState>() {
    u = NULL;
    ncv = -1;
  }
  
  virtual ~LpSolid() {}

  void report() const {
  
    const int lpSize  = getLpSize();
  
    if (mpi_rank==0) cout << " > report lpSolid:" << endl;
  
    cti_real* r1 = new cti_real[this->size()];
    cti_real (*r2)[3] = new cti_real[this->size()][3];
    for (int ip = 0; ip < this->size(); ++ip) {
      FOR_I3 r2[ip][i] = this->lp[ip].xp[i];
    }
    dumpRange(r2,this->size(),"xp");
    for (int ip = 0; ip < this->size(); ++ip) {
    FOR_I3 r2[ip][i] = this->lp[ip].up[i];
    }
    dumpRange(r2,this->size(),"up");
    for (int ip = 0; ip < this->size(); ++ip) {
      r1[ip] = this->lp[ip].Tp;
    }
    dumpRange(r1,this->size(),"Tp");
    for (int ip = 0; ip < this->size(); ++ip) {
      r1[ip] = this->lp[ip].mp;
    }
    dumpRange(r1,this->size(),"mp");
    for (int ip = 0; ip < this->size(); ++ip) {
      r1[ip] = this->lp[ip].dp;
    }
    dumpRange(r1,this->size(),"dp");
    for (int ip = 0; ip < this->size(); ++ip) {
      r1[ip] = this->lp[ip].npar;
    }
    dumpRange(r1,this->size(),"npar");
  
    IF_RANK0 cout << " > lpSolid size: " << lpSize << endl;
  
    delete[] r1;
    delete[] r2;
  
  }

};
  
class LpLiquid : public LpContainer<LpLiquidState> {
private:
public:
  double (*u)[3];
  int ncv;
  vector<LpBcBase*> bcVec;

  LpLiquid() : LpContainer<LpLiquidState>() {
    u = NULL;
    ncv = -1;
  }
  
  virtual ~LpLiquid() {}

  void report() const {
  
    const int lpSize  = getLpSize();
  
    if (mpi_rank==0) cout << " > report lpLiquid:" << endl;
  
    cti_real* r1 = new cti_real[this->size()];
    cti_real (*r2)[3] = new cti_real[this->size()][3];
    for (int ip = 0; ip < this->size(); ++ip) {
      FOR_I3 r2[ip][i] = this->lp[ip].xp[i];
    }
    dumpRange(r2,this->size(),"xp");
    for (int ip = 0; ip < this->size(); ++ip) {
    FOR_I3 r2[ip][i] = this->lp[ip].up[i];
    }
    dumpRange(r2,this->size(),"up");
    for (int ip = 0; ip < this->size(); ++ip) {
      r1[ip] = this->lp[ip].Tp;
    }
    dumpRange(r1,this->size(),"Tp");
    for (int ip = 0; ip < this->size(); ++ip) {
      r1[ip] = this->lp[ip].mp;
    }
    dumpRange(r1,this->size(),"mp");
    for (int ip = 0; ip < this->size(); ++ip) {
      r1[ip] = this->lp[ip].dp;
    }
    dumpRange(r1,this->size(),"dp");
    for (int ip = 0; ip < this->size(); ++ip) {
      r1[ip] = this->lp[ip].npar;
    }
    dumpRange(r1,this->size(),"npar");
    for (int ip = 0; ip < this->size(); ++ip) {
      r1[ip] = this->lp[ip].tbu;
    }
    dumpRange(r1,this->size(),"tbu");
  
    IF_RANK0 cout << " > lpLiquid size: " << lpSize << endl;
  
    delete[] r1;
    delete[] r2;
  
  }
};
  
template<class T>
class LpClass : public LpContainer<T> {

public:
  LpClass() {}

  virtual ~LpClass() {}

  // intended to reallocate the whole lp
  // this is not to be used in the solver, but a way to read multiple snapshots one at a time
  void realloc(const int NP) {
    LpContainer<T>::np_max = NP;
    allocLpWithOutCopy(LpContainer<T>::np_max);
    LpContainer<T>::np = NP;
  }

private:

  void allocLpWithOutCopy(const int NP) {
    T * lp_new = new T [NP];
    if (LpContainer<T>::lp !=NULL) delete[] LpContainer<T>::lp;
    LpContainer<T>::lp = lp_new;
  }

public:
  void addLp(T &_lp) {
    LpContainer<T>::resize(LpContainer<T>::np+1);
    LpContainer<T>::lp[LpContainer<T>::np-1].copy(_lp); // np is increased so np-1 is the right position for the new _lp
  }

};

#endif
