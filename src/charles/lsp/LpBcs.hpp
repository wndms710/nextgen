#ifndef __LPBC_HPP__
#define __LPBC_HPP__

#include "StaticSolver.hpp"
#include "LpStuff.hpp"

class LpBcBase {
private:
  BfZone * bf_zone;
  StaticSolver * solver;
public:
  LpBcBase(BfZone * bf_zone,StaticSolver * solver) {
    this->bf_zone = bf_zone;
    this->solver = solver;
  }
  virtual ~LpBcBase() {}
  // return false if particle is removed by bc, true is particle is to remain in the simulation...
  // dx is the vector from the particle to the intersection point on the relevant tri... 
  // st_normal is the outward normal of the relevant surface tri with area magnitude...
  virtual bool applyBc(LpTracerState *lp,const int ibf,const int ist_ss,const double dx[3],const double st_normal[3]) = 0;
};

class LpTracerBcOutlet : public LpBcBase {
private:
  int8 count;
public:
  LpTracerBcOutlet(BfZone * bf_zone,StaticSolver * solver) : LpBcBase(bf_zone,solver) {
    if (mpi_rank == 0) cout << " > LP_TRACER.BC " << bf_zone->getName() << " OUTLET" << endl;
    count = 0;
  }
  bool applyBc(LpTracerState *lp,const int ibf,const int ist_ss,const double dx[3],const double st_normal[3]) {
    ++count;
    return false;
  }
};

class LpTracerBcBounce : public LpBcBase {
private:
  int8 count;
public:
  LpTracerBcBounce(BfZone * bf_zone,StaticSolver * solver) : LpBcBase(bf_zone,solver) {
    if (mpi_rank == 0) cout << " > LP_TRACER.BC " << bf_zone->getName() << " BOUNCE" << endl;
    count = 0;
  }
  bool applyBc(LpTracerState *lp,const int ibf,const int ist_ss,const double dx[3],const double st_normal[3]) {
    const double dp = DOT_PRODUCT(dx,st_normal); // should be negative
    assert(dp <= 0.0);
    const double st_normal_mag2 = DOT_PRODUCT(st_normal,st_normal);
    assert(st_normal_mag2 > 0.0);
    FOR_I3 lp->xp[i] += 2.0*dp*st_normal[i]/st_normal_mag2;
    const double un = DOT_PRODUCT(lp->up,st_normal);
    FOR_I3 lp->up[i] -= 2.0*un*st_normal[i]/st_normal_mag2;
    ++count;
    return true;
  }
};

class LpSolidBcOutlet : public LpBcBase {
private:
  int8 count;
public:
  LpSolidBcOutlet(BfZone * bf_zone,StaticSolver * solver) : LpBcBase(bf_zone,solver) {
    if (mpi_rank == 0) cout << " > LP_SOLID.BC " << bf_zone->getName() << " OUTLET" << endl;
    count = 0;
  }
  bool applyBc(LpTracerState *lp,const int ibf,const int ist_ss,const double dx[3],const double st_normal[3]) {
    ++count;
    return false;
  }
};

class LpSolidBcBounce : public LpBcBase {
private:
  int8 count;
public:
  LpSolidBcBounce(BfZone * bf_zone,StaticSolver * solver) : LpBcBase(bf_zone,solver) {
    if (mpi_rank == 0) cout << " > LP_SOLID.BC " << bf_zone->getName() << " BOUNCE" << endl;
    count = 0;
  }
  bool applyBc(LpTracerState *lp,const int ibf,const int ist_ss,const double dx[3],const double st_normal[3]) {
    // you should use dynamic_cast if you need variables in LpSolidState
    const double dp = DOT_PRODUCT(dx,st_normal); // should be negative
    assert(dp <= 0.0);
    const double st_normal_mag2 = DOT_PRODUCT(st_normal,st_normal);
    assert(st_normal_mag2 > 0.0);
    FOR_I3 lp->xp[i] += 2.0*dp*st_normal[i]/st_normal_mag2;
    const double un = DOT_PRODUCT(lp->up,st_normal);
    FOR_I3 lp->up[i] -= 2.0*un*st_normal[i]/st_normal_mag2;
    ++count;
    return true;
  }
};

class LpLiquidBcOutlet : public LpBcBase {
private:
  int8 count;
public:
  LpLiquidBcOutlet(BfZone * bf_zone,StaticSolver * solver) : LpBcBase(bf_zone,solver) {
    if (mpi_rank == 0) cout << " > LP_LIQUID.BC " << bf_zone->getName() << " OUTLET" << endl;
    count = 0;
  }
  bool applyBc(LpTracerState *lp,const int ibf,const int ist_ss,const double dx[3],const double st_normal[3]) {
    ++count;
    return false;
  }
};

class LpLiquidBcBounce : public LpBcBase {
private:
  int8 count;
public:
  LpLiquidBcBounce(BfZone * bf_zone,StaticSolver * solver) : LpBcBase(bf_zone,solver) {
    if (mpi_rank == 0) cout << " > LP_LIQUID.BC " << bf_zone->getName() << " BOUNCE" << endl;
    count = 0;
  }
  bool applyBc(LpTracerState *lp,const int ibf,const int ist_ss,const double dx[3],const double st_normal[3]) {
    // you should use dynamic_cast if you need variables in LpLiquidState
    const double dp = DOT_PRODUCT(dx,st_normal); // should be negative
    assert(dp <= 0.0);
    const double st_normal_mag2 = DOT_PRODUCT(st_normal,st_normal);
    assert(st_normal_mag2 > 0.0);
    FOR_I3 lp->xp[i] += 2.0*dp*st_normal[i]/st_normal_mag2;
    const double un = DOT_PRODUCT(lp->up,st_normal);
    FOR_I3 lp->up[i] -= 2.0*un*st_normal[i]/st_normal_mag2;
    ++count;
    return true;
  }
};




#endif
