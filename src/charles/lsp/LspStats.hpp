#ifndef _LSPSTATS_HPP_
#define _LSPSTATS_HPP_

#include "Lsp.hpp"
#include "CtiLiquid.hpp"

class LspStats {

protected:

  string name;
  string param_line;
  int interval;

  double xc[3];
  double normal[3];
  double t1[3]; // plane tangent dir 1
  double t2[3]; // plane tangent dir 2

  // stats will be in...
  int nbins;
  double dt_sum;
  double * np_sum;
  double * nm_sum;
  double * mp_sum;
  double * mm_sum;
  double * mup_sum;
  double * mum_sum;
  double * d2p_sum;
  double * d2m_sum;
  double * d3p_sum;
  double * d3m_sum;

  // PDF variables
  bool got_pdf;
  int nbins_pdf;
  double dmin;
  double dmax;
  double* diam;
  double* pdf_mass_p;
  double* pdf_mass_m;
  double* pdf_num_p;
  double* pdf_num_m;
  
  virtual int getBin(const double dxc[]) = 0;
  int         getBin_pdf(const double);
  double      getD_pdf(const int);

  void writePdf(const double time,const int step,const bool verbose);
  string getParamLine(Param * param);

public:

  LspStats(Param * param);
  virtual ~LspStats();

  int getInterval() const { return interval; }
  virtual void init(Param * param) = 0;
  void update(const LspState * lsp,const int np,CtiMaterial* fuel,const double dt,const bool verbose);
  void reset(const bool verbose);
  virtual void write(const double time,const int step,const bool verbose) = 0;
};

class LspStatsDisk : public LspStats {
private:

  double r0,r1;
  double theta0,theta1;
  int nr,ntheta;
  
  int getBin(const double dxc[]);

public:
  LspStatsDisk(Param * param);
  virtual ~LspStatsDisk() {};
  void init(Param * param);
  void write(const double time,const int step,const bool verbose);
};

class LspStatsRect : public LspStats {
private:

  double s1,s2;  // rectangle side lengths
  int nbins_s1, nbins_s2;

  int getBin(const double dxc[]);

public:
  LspStatsRect(Param * param);
  virtual ~LspStatsRect() {};
  void init(Param * param);
  void write(const double time,const int step,const bool verbose);
};

LspStats * newLspStats(Param * param);

enum PdfGeom{
  BOX,
  CONE,
};

class LspPDF {

private:

  string name;
  int interval;
  int start;
  PdfGeom geom;

  // geometry variables
  double xmin[3], xmax[3];        // BOX geometry
  double xc[3], normal[3], hight, r0, r1; // cone geometry

  // PDF variables
  int nbins;
  double dmin;
  double dmax;
  double* pdf;
  double* diam;

public:

  LspPDF();
  ~LspPDF();

  int getInterval() const { return interval; }
  int getStart() const { return start; }

  void init(Param * param);
  void update(const LspState * lsp,const int np,const double dt,const bool verbose);
  void reset(const bool verbose);
  int getBin(const double d);
  double getD(const int ibin);
  void write(const double time,const int step,const bool verbose);
  void report();
  bool inside(const LspState& lsp);

};

#endif
