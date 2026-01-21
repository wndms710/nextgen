#ifndef _SPLINE_STUFF_HPP_
#define _SPLINE_STUFF_HPP_

#include "Common.hpp"

namespace SplineStuff {
  
  class CubicSpline {
  
  private:
    int np;
    double (*xp)[3];
    double *sp;
    double (*dp)[3];
    double length;
    double s_gauss5[5];
    double wgt_gauss5[5];
    bool b_periodic;
  
  public:
    CubicSpline() {
      xp = NULL;
      sp = NULL;
      dp = NULL;
      b_periodic = false;
    }
    
    ~CubicSpline() {
      DELETE(xp);
      DELETE(sp);
      DELETE(dp);
    }
    
    // initize a spline with a set of points... 
    void init(const double (*x)[3],const int n,const bool b_periodic = false);

    // forming point access functions...
    int getNp() const { return np; }
    void getXp(double xpi[3],const int i) const; // get forming point xp[i]
    double getSp(const int i) const; // get fractional distance of the forming point sp[i] 
    
    // spline functions...
    double getLength() const { return length; }
    void getX(double x[3],const double s) const; // 0 <= s <= 1
    void getXAndDx(double x[3],double dx[3],const double s) const; // 0 <= s <= 1

  private:
    // spline interval functions...
    void getX(double x[3],const double t,const int ip0) const;
    void getDxdt(double dxdt[3],const double t,const int ip0) const;
  };
  
  class CubicSplineXy {
  private:
    int np;
    double *xp; // x values
    double *yp; // y values
    double *dx; // x[i+1] - x[i]
    double *dp; // second derivative of y
    //double s_gauss5[5];
    //double wgt_gauss5[5];
  public:
    CubicSplineXy() {
      xp = NULL;
      yp = NULL;
      dx = NULL;
      dp = NULL;
    }
    ~CubicSplineXy() {
      DELETE(xp);
      DELETE(yp);
      DELETE(dx);
      DELETE(dp);
    }
    void init(const double *x, const double *y, const int n);
    double getY(const double x) const;
    int getIp(const double x) const;
  };
  
}

#endif
