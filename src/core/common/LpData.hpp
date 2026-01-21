#ifndef __LPDATA_HPP__
#define __LPDATA_HPP__

// =============================
// flag of Lp for recycling
// =============================
enum {
  TRASH = 10,
  KEEP
};

// =============================
// Solid/Liquid Lp coupling mode
// =============================
enum {
  NOT_COUPLED,
  ONE_WAY_COUPLED,
  TWO_WAY_COUPLED
};

struct LpDataBase {
  int icv;
  int flag;
  double xp[3];
  double dp;
};

struct LpData : public LpDataBase{
  double up[3];
  double mp;
  double Tp;
};

#endif
