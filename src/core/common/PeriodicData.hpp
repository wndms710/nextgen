#ifndef _PERIODIC_DATA_HPP_
#define _PERIODIC_DATA_HPP_

#include "CTI.hpp"
using namespace CTI;

#include "PeriodicTransform.hpp"

class PeriodicRt {
public:
  bool b_R,b_t;
  double R[9];
  double t[3];
  PeriodicRt() {
    // default is no transform...
    b_R = b_t = false;
  }
};

namespace PeriodicData {

  extern vector<PeriodicTransform> periodicTransformVec;

  // periodicity is common to all parts...

  bool isPeriodic();
  bool checkPeriodicBitPair(const int i);
  bool getPeriodicR(double R[9],const int bits);
  bool getPeriodicT(double t[3],const int bits);
  void periodicTranslate(double (*xp_t)[3],const int n,const int bits);
  void periodicTranslate(double *xp_t,const int n,const int bits);
  void periodicRotate(double (*xp_t)[3],const int n,const int bits);
  void periodicRotate(double *xp_t,const int n,const int bits);
  void scalePeriodicTransformVec(double s[3]);
  void buildPeriodicRtVec(vector<PeriodicRt>& periodicRtVec);
  void buildBitVec(vector<int>& bitsVec);

}

#endif
