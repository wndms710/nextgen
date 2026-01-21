#ifndef _SIMPLE_TRI_HPP_
#define _SIMPLE_TRI_HPP_

class SimpleTri {
public:
  double x0[3];
  double x1[3];
  double x2[3];
  SimpleTri() {}
  SimpleTri(const SimpleTri& other) {
    FOR_I3 x0[i] = other.x0[i];
    FOR_I3 x1[i] = other.x1[i];
    FOR_I3 x2[i] = other.x2[i];
  }
  SimpleTri(const double _x0[3],const double _x1[3],const double _x2[3]) {
    FOR_I3 x0[i] = _x0[i];
    FOR_I3 x1[i] = _x1[i];
    FOR_I3 x2[i] = _x2[i];
  }
  SimpleTri(const double _x0[3],const double _x1[3],const double _x2[3],const double _dx[3]) {
    FOR_I3 x0[i] = _x0[i]+_dx[i];
    FOR_I3 x1[i] = _x1[i]+_dx[i];
    FOR_I3 x2[i] = _x2[i]+_dx[i];
  }
};

class SimpleTriWithData: public SimpleTri {
public:
  double d0;
  double d1;
  double d2;
  SimpleTriWithData(const double _x0[3],const double _x1[3],const double _x2[3],const double _d0,const double _d1,const double _d2) : SimpleTri(_x0,_x1,_x2) {
    d0 = _d0;
    d1 = _d1;
    d2 = _d2;
  }
  SimpleTriWithData(const double _x0[3],const double _x1[3],const double _x2[3],const double _dx[3],const double _d0,const double _d1,const double _d2) : SimpleTri(_x0,_x1,_x2,_dx) {
    d0 = _d0;
    d1 = _d1;
    d2 = _d2;
  }
};

class SimpleTriWithWgts: public SimpleTri {
public:
  int i0_l,i0_r;
  int i1_l,i1_r;
  int i2_l,i2_r;
  double wgt0_l; // wgt_r == 1-wgt_l
  double wgt1_l; // wgt_r == 1-wgt_l
  double wgt2_l; // wgt_r == 1-wgt_l
  SimpleTriWithWgts(const double _x0[3],const double _x1[3],const double _x2[3],const int _i0_l,const int _i0_r,const double _wgt0_l,const int _i1_l,const int _i1_r,const double _wgt1_l,const int _i2_l,const int _i2_r,const double _wgt2_l) : SimpleTri(_x0,_x1,_x2) {
    i0_l = _i0_l; i0_r = _i0_r; wgt0_l = _wgt0_l;
    i1_l = _i1_l; i1_r = _i1_r; wgt1_l = _wgt1_l;
    i2_l = _i2_l; i2_r = _i2_r; wgt2_l = _wgt2_l;
  }
};

#endif
