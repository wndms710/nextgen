#ifndef LSP_INJECT_HPP
#define LSP_INJECT_HPP

#include "StaticSolver.hpp"
#include "CtiLiquid.hpp"
#include "CtiSolid.hpp"
#include "SubSurface.hpp"
#include "GeomUtils.hpp"

using namespace std;

// different injection kinds
enum InjectorType {
  COUNT_TYPE,
  MDOT_TYPE,
  STRUCTURED_TYPE,
  MASS_TYPE,
  EACH_CELL_TYPE
};

// distribution of the injected particles
enum InjectorDist {
  UNIFORM_DIST,
  RR_DIST,
  GAUSSIAN_DIST
};

// injector geometry kinds
enum InjectorGeom {
  POINT_GEOM,
  POINT3_GEOM,
  CIRCLE_GEOM,
  RING_GEOM,
  LPBOX_GEOM,
  ZONE_GEOM,
  CONE_GEOM,
  LPPLANE_GEOM
};

// breakup models
enum BreakupModel {
  NO_BREAKUP,
  SBM_BREAKUP
};

enum cdfDistType {
  NA,
  CDF_COUNT_DIST,
  CDF_MASS_DIST
};

enum velocityMode {
  DEFAULT,
  MULT,
  MAG,
  MAX
};

enum pShape {
  SPHERE_SHAPE,
  ELLIPSOID_SHAPE
};

typedef struct {
  double xp[3];
  double up[3];
  double mp;
  double Tp;
  int icv;
  double e; // elongation for ellipsoid
  double f; // flatness for ellipsoid
} LspData;

class LspInjector {
protected:

  string name;
  
  int InjectorType;     // injector type: MDOT_TYPE/COUNT_TYPE

  double Tp;            // particle temperature
  double rhop;          // particle density

  double mdot;          // mass flux through injector for MDOT_TYPE
  int np;               // number of particles injected in the COUNT_TYPE mode

  int DistType;         // injected particle distribution
  double DistData[5];

  int cdfType;          // cdf type either count or mass distribution
  double * countCdf_x;  // count cdf x axis
  double * countCdf_y;  // count cdf y axis
  int countCdf_n;       // count cdf nbins
  double c0, c1;        // minimum and mximum values of count cdf for dmin and dmax

  int GeomType;         // injector geometry type: POINT/...
  double GeomData[11];  // data for injector geometry

  double ResidualMass;

  double injectTime;

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
  int overWriteVelocity_mode;
  double OWVData[2];

  int pShape;
  double axis_std, axis_e, axis_f;

public:

  LspInjector();
  virtual ~LspInjector();
  virtual void initFromParams(Param *, CtiMaterial *, StaticSolver *, double (*u)[3]);
  string getName() const { return(name); }
  virtual void addNewLsp(vector<LspData>&, const double, const double, const bool);
  double diam_rr_countCdf(const double &, const double &, const double &, const double &);
  double diam_rr_massCdf(const double &, const double &, const double &, const double &);
  double diam_gaussian(const double &, const double &, const double &, const double &);
  double getInjectTime() {return injectTime;}
  void build_countCdf_stuff();
  int initZoneInjector(string);
  int initPlaneInjector();
  void overWriteVelocity(vector<LspData>& ,int, bool verbose=false);
  int pointIsInside(const double xp[3], int &icv_ret);
};

class StaticLspInjector : public LspInjector {
protected:

  int npx, npy, npz;    // structured static injector vars
  
  double total_mass;    // total injected mass of by the injector -- for MASS_TYPE
  
  int np_in_cell;       // number of particles in each cell for EACH_CELL_TYPE

public:

  StaticLspInjector();
  virtual ~StaticLspInjector();
  void initFromParams(Param *, CtiMaterial *, StaticSolver *, double (*u)[3]);
  void addNewLsp(vector<LspData>&,const double time = 0.0,const double dt = 0.0,const bool verbose = true);
  int getClosestICVFromList(vector<int> &, const double*);
};

#endif
