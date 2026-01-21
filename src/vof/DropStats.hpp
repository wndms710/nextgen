#ifndef _DROPSTATS_HPP_
#define _DROPSTATS_HPP_

#include "Lsp.hpp"
#include "CtiLiquid.hpp"
#include "LspStats.hpp"

class DropPDF {

private:

  string name;
  int interval;
  int start;
  PdfGeom geom;
  string varname;

  // geometry variables
  double xmin[3], xmax[3];        // BOX geometry
  double xc[3], normal[3], hight, r0, r1; // cone geometry

  // PDF variables
  int nbins;
  double dmin;
  double dmax;
  double* pdf;
  double* num;
  double* diam;

  //double dt_sum;

public:

  DropPDF();
  ~DropPDF();

  int getInterval() const { return interval; }
  int getStart() const { return start; }

  void init(Param * param);
  void update(const double* dropD, const double (*dropX)[3], const double (*dropU)[3], const int ndrops, const double time,const double dt);
  void update(const double* blob_surf, const double* blob_vol, const double (*blobX)[3], const int nblobs, const double time,const double dt,const bool verbose);
  void update(const double* blob_surf, const double* blob_vol, const double (*blobX)[3],  const int nblobs, 
      const LspState * lsp,const int np, CtiMaterial* fuel,const double time, const double dt,const bool verbose);
  void reset(const bool verbose);
  int getBin(const double d);
  double getD(const int ibin);
  void write(const double time,const int step,const bool verbose);
  void report();
  bool inside(const double xp[3]);

};


#endif
