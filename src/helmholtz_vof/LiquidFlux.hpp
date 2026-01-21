#ifndef _LIQUIDFLUX_HPP_
#define _LIQUIDFLUX_HPP_

#include "FlowSolver.hpp"

class LiquidFlux {

private:

  string name;
  int sample_interval;
  int write_interval;
  string varname;
 
public:

  LiquidFlux();
  ~LiquidFlux();

  // geometry variables
  double ymin, ymax, zmin, zmax;
  int ny, nz;
  int nbins;
  double* liquid_flux;
  double* projected_area;
  double xp[3];
  double *yp;
  double *zp;

  vector<pair<int,int> > groupFaceVec;
  vector<int> ibinVec;
  vector<int*> groupVec;

  // function list
  void init(Param * param);
  void report();
  int getSampleInterval() const { return sample_interval; }
  int getWriteInterval() const { return write_interval; }
  void calcProjectedArea(double *local_area);
  void calcLiquidFlux(double *local_flux);
  void write(const double time,const int step);

};

#endif
