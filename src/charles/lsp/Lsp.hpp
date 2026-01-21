#ifndef _LSP_HPP_
#define _LSP_HPP_

#include "LpStuff.hpp"
//#include "Transform.hpp"

class LspRhs {
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
  inline void axpy(const double& val, const LspRhs& rhs) {
    for (unsigned int i =0; i <3; ++i)
      xp[i] += val*rhs.xp[i];
    for (unsigned int i =0; i <3; ++i)
      up[i] += val*rhs.up[i];
    mp += val*rhs.mp;
    Tp += val*rhs.Tp;
  }

};

class LspState : public LpTracerState {
public:
  //dd we can carry only one of mp and dp 
  //int icv;
  //int flag
  //double xp[3];
  //double xp0[3];
  //double up[3];
  //double dp;
  double Tp;
  double mp;
  double npar; // number of particles in each parcel

  // just for inspection
  double mp0;

  // auxillary variable for implicit time integration
  // do not need interprocessor trasport or be copied 
  double k;   // equals f1/tau in the momentum eqn
  double kt;
  // boundary coefficients
  double kb;
  double ktb;

  //===============
  // breakup params
  //===============
  double tbu;

  //===============
  //for ellipsoids
  //===============
  double S, I;    // the short and intermediate dimensions of the ellipsoid. The long side is computed from these and dp
  double e, f;

  static int data_size() { return 17; }
  //static string desc() {return "rho:r1;u:r2;rhoE:r1;T:r1;p:r1;";}

  inline void pack(double * buf) const {  // include a transform in the future
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
    buf[10] = Tp;
    buf[11] = mp;
    buf[12] = dp;
    buf[13] = tbu;
    buf[14] = npar;
    buf[15] = e;
    buf[16] = f;
  }

  inline void packRt(double * buf,const double R[9], const double t[3]) const {
    int * icv_flag = (int*)buf; // use the first 8 byte double to store 2 4 byte ints
    icv_flag[0] = icv;
    icv_flag[1] = flag;
    // note: apply R and t here...
    //buf[1] = xp[0];
    //buf[2] = xp[1];
    //buf[3] = xp[2];
    MiscUtils::applyInvRotation(buf+1,R,xp);
    buf[1] -= t[0];
    buf[2] -= t[1];
    buf[3] -= t[2];

    //buf[4] = xp0[0];
    //buf[5] = xp0[1];
    //buf[6] = xp0[2];
    MiscUtils::applyInvRotation(buf+4,R,xp0);
    buf[4] -= t[0];
    buf[5] -= t[1];
    buf[6] -= t[2];

    //buf[7] = up[0];
    //buf[8] = up[1];
    //buf[9] = up[2];
    MiscUtils::applyInvRotation(buf+7,R,up);

    buf[10] = Tp;
    buf[11] = mp;
    buf[12] = dp;
    buf[13] = tbu;
    buf[14] = npar;
    buf[15] = e;
    buf[16] = f;

  }

  inline void packR(double * buf,const double R[9]) const {
    int * icv_flag = (int*)buf; // use the first 8 byte double to store 2 4 byte ints
    icv_flag[0] = icv;
    icv_flag[1] = flag;
    MiscUtils::applyInvRotation(buf+1,R,xp);
    MiscUtils::applyInvRotation(buf+4,R,xp0);
    MiscUtils::applyInvRotation(buf+7,R,up);
    buf[10] = Tp;
    buf[11] = mp;
    buf[12] = dp;
    buf[13] = tbu;
    buf[14] = npar;
    buf[15] = e;
    buf[16] = f;
  }

  inline void packt(double * buf, const double t[3]) const {
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
    buf[10] = Tp;
    buf[11] = mp;
    buf[12] = dp;
    buf[13] = tbu;
    buf[14] = npar;
    buf[15] = e;
    buf[16] = f;
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
    Tp = buf[10];
    mp = buf[11];
    dp = buf[12];
    tbu = buf[13];
    npar = buf[14];
    e = buf[15];
    f = buf[16];
  }

  inline void copy(const LspState& other) {
    icv = other.icv;
    flag   = other.flag;
    xp[0] = other.xp[0];
    xp[1] = other.xp[1];
    xp[2] = other.xp[2];
    xp0[0] = other.xp0[0];
    xp0[1] = other.xp0[1];
    xp0[2] = other.xp0[2];
    up[0] = other.up[0];
    up[1] = other.up[1];
    up[2] = other.up[2];
    Tp    = other.Tp;
    mp    = other.mp;
    dp    = other.dp;
    tbu   = other.tbu;
    npar   = other.npar;
    e     = other.e;
    f     = other.f;
  }

  inline void add(const LspRhs& other) {
    xp[0] += other.xp[0];
    xp[1] += other.xp[1];
    xp[2] += other.xp[2];
    up[0] += other.up[0];
    up[1] += other.up[1];
    up[2] += other.up[2];
    Tp    += other.Tp;
    mp    += other.mp;
  }

  inline void add_im(const LspRhs& rhs, const LspRhs& rhs_d, const double& dt) {
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

#endif
