#ifndef __LP_INJECTOR_HPP__
#define __LP_INJECTOR_HPP__

#include "CtiLiquid.hpp"
#include "CtiSolid.hpp"
#include "LpData.hpp"

using namespace std;

namespace LpInjectorNS {

  // different injection kinds
  enum InjectorKind {
    COUNT_INJ,
    COUNTDOT_INJ,
    STRUCTURED_INJ,
    EACHCELL_INJ,
    
    // only used for solid/liquid particle injectors
    MDOT_INJ,
  };
  
  // injector geometry kinds
  enum GeomKind {
    POINT_GEOMETRY,     // with COUNT_INJ:COUNTDOT_INJ
    CIRCLE_GEOMETRY,    // with COUNT_INJ:COUNTDOT_INJ:STRUCTURED_INJ
    BOX_GEOMETRY,       // with COUNT_INJ:COUNTDOT_INJ:STRUCTURED_INJ:EACHCELL_INJ
    ZONE_GEOMETRY,      // with COUNT_INJ:COUNTDOT_INJ
    PLANE_GEOMETRY,     // with COUNT_INJ:COUNTDOT_INJ
  
    // only used for solid/liquid particle injectors
    RING_GEOMETRY,
    CONE_GEOMETRY,
  };
  
  // distribution of the injected particles
  enum DistKind {
    UNIFORM_DIST,
    RR_DIST,
    GAUSSIAN_DIST
  };
  
  enum cdfKind {
    NA,
    CDF_COUNT_DIST,
    CDF_MASS_DIST
  };

}

// this class is used for tracer particle injection
class LpInjectorTracer {
protected:

  string name;
  
  int InjectorKind;     

  double npdot;         // rate of injection for COUNTDOT_INJ mode
  int np;               // number of particles injected for COUNT_INJ/STRUCTURED mode
  int npx, npy, npz;    // number of particles in each dimension for STRUCTURED mode
  int np_in_cell;        // number of particles in each cell for EACHCELL mode

  double dp;            // particle diameter

  int GeomKind;         // injector geometry kind: POINT/...
  double GeomData[11];  // data for injector geometry

  double ResidualCount;

  double time0;         // time of the solver at the initiation
  double injectDt;      // injection time

  // for zone inlets
  int nst;
  int nsp;
  int (* spost)[3];
  double (* xp)[3];
  double * cdf_area;

  // for plane inlets
  vector<pair<SimpleTri,int> > triVec; 
  double * cdf_area_tri;
  double * cdf_area_rank;
  double e1[3];
  double e2[3];

  StaticSolver * solver;

  int injector_id;

public:

  LpInjectorTracer();
  virtual ~LpInjectorTracer();
  string getName() const { return(name); }
  double getInjectDt() {return injectDt;}
  virtual void initFromParams(Param *, StaticSolver *, double time0, int injector_id);
  virtual void addNewLp(vector<LpDataBase>&, const double time, const double dt, const bool verbose);
};

// this class is used for solid/liquid particle injection
class LpInjector {
protected:

  string name;
  
  int InjectorKind;     

  double Tp;            // particle temperature
  double rhop;          // particle density

  double mdot;          // mass flux through injector for MDOT_INJ mode
  int np;               // number of particles injected in the COUNT_INJ mode

  int DistKind;         // injected particle distribution
  double DistData[5];

  int cdfKind;          // cdf type either count or mass distribution
  double * countCdf_x;  // count cdf x axis
  double * countCdf_y;  // count cdf y axis
  int countCdf_n;       // count cdf nbins
  double c0, c1;        // minimum and mximum values of count cdf for dmin and dmax

  int GeomKind;         // injector geometry kind: POINT/...
  double GeomData[11];  // data for injector geometry

  double ResidualMass;

  double time0;         // time of the solver at the initiation
  double injectDt;

  // for zone inlets
  int nst;
  int nsp;
  int (* spost)[3];
  double (* xp)[3];
  double * cdf_area;

  // for cone and plane inlets
  double e1[3];
  double e2[3];

  // for plane inlets
  vector<pair<SimpleTri,int> > triVec; 
  double * cdf_area_tri;
  double * cdf_area_rank;

  StaticSolver * solver;

  // velocity from solver to overwrite particles' velocity
  double (*u)[3];
  
  bool overWriteVelocity_b;

  int injector_id;
  int material_id;
  int fuel_id;
  int dust_id;

public:

  LpInjector();
  virtual ~LpInjector();
  virtual void initFromParams(Param *, StaticSolver *, double time0, vector<CtiLiquid*> fuelVec, vector<CtiSolid*> dustVec, double (*u)[3], int injector_id);
  string getName() const { return(name); }
  virtual void addNewLp(vector<LpData>&, const double, const double, const bool);
  double diam_rr_countCdf(const double &, const double &, const double &, const double &);
  double diam_rr_massCdf(const double &, const double &, const double &, const double &);
  double diam_gaussian(const double &, const double &, const double &, const double &);
  double getInjectDt() {return injectDt;}
  void build_countCdf_stuff();
  void overWriteVelocity(vector<LpData>& ,int, bool verbose=false);
};

#endif
