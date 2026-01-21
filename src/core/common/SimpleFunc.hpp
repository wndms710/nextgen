#ifndef _SIMPLE_FUNC_HPP_
#define _SIMPLE_FUNC_HPP_

#include "CTI.hpp"

class SimpleFunc {
public:
  virtual double eval(const double x) = 0;
  virtual void eval(double * y,const double * const x,const int n) = 0;
  virtual ~SimpleFunc() {};
};

class ConstantFunc : public SimpleFunc {
private:
  double value;
public:
  ConstantFunc(const double val) {
    value = val;
  }
  double eval(const double x) {
    return value;
  }
  void eval(double * y,const double * const x,const int n) {
    for (int i = 0; i < n; ++i)
      y[i] = value;
  }
};

class TableFunc : public SimpleFunc {
private:
  int n_table;
  double * x_table;
  double * y_table; 
public:
  TableFunc(vector<pair<double,double> >& xyVec) {
    if (xyVec.size() < 2) {
      CERR("TABLE functions require atleast two entries. Use CONSTANT for a single value");
    }
    n_table = xyVec.size();
    x_table = new double[n_table];
    y_table = new double[n_table];
    for (int ii = 0; ii < n_table; ++ii) {
      x_table[ii] = xyVec[ii].first;
      y_table[ii] = xyVec[ii].second;
      if (ii > 0) {
        if (x_table[ii] <= x_table[ii-1]) {
          delete[] x_table;
          delete[] y_table;
          CERR("TABLE functions must be monotonically increasing with no repeated values");
        }
      }
    }
  }
  ~TableFunc() {
    delete[] x_table;
    delete[] y_table;
  }
  double eval(const double x) {
    if (x <= x_table[0])
      return y_table[0];
    else if (x >= x_table[n_table-1])
      return y_table[n_table-1];
    else {
      // bisection...
      int i0 = 0;
      int i1 = n_table-1;
      while (i1 > i0+1) {
        const int imid = (i0+i1)/2;
        if (x <= x_table[imid]) {
          i1 = imid;
        }
        else {
          i0 = imid;
        }
      }
      return ((x_table[i1]-x)*y_table[i0] + (x-x_table[i0])*y_table[i1])/(x_table[i1]-x_table[i0]);
    }
  }
  void eval(double * y,const double * const x,const int n) {
    for (int i = 0; i < n; ++i)
      y[i] = eval(x[i]);
  }
};

class ClippedCubicFunc : public SimpleFunc {
private:
  double c0,c1,c2,c3;
  double x0,x1;
public:
  ClippedCubicFunc(const double coeff[4],const double range[2]) {
    c0 = coeff[0];
    c1 = coeff[1];
    c2 = coeff[2];
    c3 = coeff[3];
    assert(range[1] > range[0]);
    x0 = range[0];
    x1 = range[1];
  }
  double eval(const double x) {
    if (x <= x0)
      return c0 + c1*x0 + c2*x0*x0 + c3*x0*x0*x0;
    else if (x >= x1)
      return c0 + c1*x1 + c2*x1*x1 + c3*x1*x1*x1;
    else 
      return c0 + c1*x + c2*x*x + c3*x*x*x;
  }
  void eval(double * y,const double * const x,const int n) {
    for (int i = 0; i < n; ++i)
      y[i] = eval(x[i]);
  }
};

class ClippedQuarticFunc : public SimpleFunc {
private:
  double c0,c1,c2,c3,c4;
  double x0,x1;
public:
  ClippedQuarticFunc(const double coeff[5],const double range[2]) {
    c0 = coeff[0];
    c1 = coeff[1];
    c2 = coeff[2];
    c3 = coeff[3];
    c4 = coeff[4];
    assert(range[1] > range[0]);
    x0 = range[0];
    x1 = range[1];
  }
  double eval(const double x) {
    if (x <= x0)
      return c0 + c1*x0 + c2*x0*x0 + c3*x0*x0*x0 + c4*x0*x0*x0*x0;
    else if (x >= x1)
      return c0 + c1*x1 + c2*x1*x1 + c3*x1*x1*x1 + c4*x1*x1*x1*x1;
    else 
      return c0 + c1*x + c2*x*x + c3*x*x*x + c4*x*x*x*x;
  }
  void eval(double * y,const double * const x,const int n) {
    for (int i = 0; i < n; ++i)
      y[i] = eval(x[i]);
  }
};

SimpleFunc * processSimpleFunc(Param * param);

#endif
