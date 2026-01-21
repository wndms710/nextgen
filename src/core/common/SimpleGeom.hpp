#ifndef _SIMPLE_GEOM_
#define _SIMPLE_GEOM_

#include "CTI.hpp"

// =======================================
// SimpleGeom is an abstract class that
// implements the routines:
// 
// bool pointIsInside(const double xp[3]);
// double pointSignedDistance(const double xp[3]);
//
// and:
//
// void getBoundingBox(double bbox[6]); 
//
// where:
//   bbox[0] = xmin
//   bbox[1] = xmax
//   bbox[2] = ymin
//   bbox[3] = ymax
//   bbox[4] = zmin
//   bbox[5] = zmax
//
// =======================================

class SimpleGeom {
public:

  virtual bool pointIsInside(const double xp[3]) = 0;
  virtual double pointSignedDistance(const double xp[3]) { // positive outside convention
    CERR("pointSignedDistance needs to be implemented for this SimpleGeom");
  }; 
  virtual void getBoundingBox(double bbox[6]) = 0;

  virtual ~SimpleGeom() {}
};

class SimpleBox : public SimpleGeom {
private:
  
  // the bottom and top corners...
  double bbmin[3],bbmax[3];
  
public:
  
  SimpleBox(const double bbmin_[3],const double bbmax_[3]) {
    FOR_I3 bbmin[i] = bbmin_[i];
    FOR_I3 bbmax[i] = bbmax_[i];
  }
  
  bool pointIsInside(const double xp[3]) {
    return (xp[0] >= bbmin[0])&&(xp[0] <= bbmax[0])&&(xp[1] >= bbmin[1])&&(xp[1] <= bbmax[1])&&(xp[2] >= bbmin[2])&&(xp[2] <= bbmax[2]);
  }

  double pointSignedDistance(const double xp[3]) {

    if (xp[0] < bbmin[0]) {
      return bbmin[0]-xp[0];
    }
    else if (xp[0] > bbmax[0]) {
      return xp[0]-bbmax[0];
    }
    else if (xp[1] < bbmin[1]) {
      return bbmin[1]-xp[1];
    }
    else if (xp[1] > bbmax[1]) {
      return xp[1]-bbmax[1];
    }
    else if (xp[2] < bbmin[2]) {
      return bbmin[2]-xp[2];
    }
    else if (xp[2] > bbmax[2]) {
      return xp[2]-bbmax[2];
    }
    else {
      double my_sds[6];
      my_sds[0] = bbmin[0]-xp[0];
      my_sds[1] = xp[0]-bbmax[0];
      my_sds[2] = bbmin[1]-xp[1];
      my_sds[3] = xp[1]-bbmax[1];
      my_sds[4] = bbmin[2]-xp[2];
      my_sds[5] = xp[2]-bbmax[2];
      int ind = -1;
      double my_min_d = HUGE_VAL;
      FOR_I6 {
        double my_d = fabs(my_sds[i]);
        if (my_d < my_min_d) {
          my_min_d = my_d;
          ind = i;
        }
      }
      assert(ind >= 0);
      return my_sds[ind];
    }

  }
  
  void getBoundingBox(double bbox[6]) {
    // careful of order here...
    bbox[0] = bbmin[0];
    bbox[1] = bbmax[0];
    bbox[2] = bbmin[1];
    bbox[3] = bbmax[1];
    bbox[4] = bbmin[2];
    bbox[5] = bbmax[2];
  }

};

class SimpleSphere : public SimpleGeom {
private:
  
  // the center and radius...
  double xc[3],r;
  
public:
  
  SimpleSphere(const double xc_[3],const double r_) {
    FOR_I3 xc[i] = xc_[i];
    r = r_;
  }
  
  bool pointIsInside(const double xp[3]) {
    return (DIST(xp, xc) <= r);
  }

  double pointSignedDistance(const double xp[3]) {
    return DIST(xp,xc)-r;
  }
  
  void getBoundingBox(double bbox[6]) {
    // careful of order here...
    FOR_I3 bbox[i*2  ] = xc[i] - r;
    FOR_I3 bbox[i*2+1] = xc[i] + r;
  }

};


class SimpleTcone : public SimpleGeom {
private:
  
  // the starting center and radius
  double x0[3],r0;
  // the ending center and radius ...
  double x1[3],r1;
  // store axis and length
  double axis[3];
  double l;
  
public:
  
  SimpleTcone(const double x0_[3],const double r0_,const double x1_[3],const double r1_) {
    FOR_I3 x0[i] = x0_[i];
    r0 = r0_;
    FOR_I3 x1[i] = x1_[i];
    r1 = r1_;
    FOR_I3 axis[i] = x1_[i] - x0_[i];
    l = MAG(axis);
  }
  
  bool pointIsInside(const double xp[3]) {    
    // compute distance from the point xp[3] to the axis formed between x0[3] -- > x1[3]
    // d = mag( (xp-x0)x(xp-x1) ) / mag(x1-x0)
    const double xpx0[3] = DIFF(xp,x0);
    const double xpx1[3] = DIFF(xp,x1);
    const double numer = CROSS_PRODUCT_MAG(xpx0,xpx1);
    const double radius = numer / l;
    
    // determine what the expected radius is at this plane perpendicular to the axis
    double dist_along_axis = DOT_PRODUCT(xpx0,axis);
    dist_along_axis /= (l*l); // normalize

    // const double r_expected = r0 + (r1-r0)*dist_along_axis/denom;
    const double r_expected = r0 + (r1-r0)*dist_along_axis;

    return ((radius <= r_expected) && (dist_along_axis >= 0.0) && (dist_along_axis <= 1.0));
  }

  double pointSignedDistance(const double xp[3]) {
    // compute distance from the point xp[3] to the axis formed between x0[3] -- > x1[3]
    // d = mag( (xp-x0)x(xp-x1) ) / mag(x1-x0)
    const double xpx0[3] = DIFF(xp,x0);
    const double xpx1[3] = DIFF(xp,x1);
    const double numer = CROSS_PRODUCT_MAG(xpx0,xpx1);
    const double radius = numer / l;
    
    // determine what the expected radius is at this plane perpendicular to the axis
    double dist_along_axis = DOT_PRODUCT(xpx0,axis);
    dist_along_axis /= l; // partially normalized: [0,l]
    const double r_expected = r0 + (r1-r0)*dist_along_axis/l;

    const double iso_r = radius-r_expected;
    const double iso_0 = -dist_along_axis;
    const double iso_1 = dist_along_axis-l;
    if (dist_along_axis <= 0.0)
      return iso_0;
    else if (dist_along_axis >= l)
      return iso_1;
    else if (radius >= r_expected)
      return iso_r;
    else {
      assert(iso_r < 0);
      assert(iso_0 < 0);
      assert(iso_1 < 0);
      return MAX3(iso_r,iso_0,iso_1);
    }

  }
  
  void getBoundingBox(double bbox[6]) {
    // careful of order here...
    // "loose" box around the surface, using sphere at both ends
    FOR_I3 bbox[i*2  ] = MIN(x0[i]-r0,x1[i]-r1);
    FOR_I3 bbox[i*2+1] = MAX(x0[i]+r0,x1[i]+r1);
  }

};


class SimpleRoundedRect_X : public SimpleGeom {
private:
  
  // the starting center, half-widths and radius
  double x0[3],dy0,dz0,r0;
  // the ending center and radius ...
  double x1[3],dy1,dz1,r1;
  // store axis and length
  double axis[3];
  double l;
  
public:
  
  SimpleRoundedRect_X(const double x0_[3],const double dy0_,const double dz0_,const double r0_,const double x1_[3],const double dy1_,const double dz1_,const double r1_) {
    FOR_I3 x0[i] = x0_[i];
    dy0 = dy0_;
    dz0 = dz0_;
    r0 = r0_;
    FOR_I3 x1[i] = x1_[i];
    dy1 = dy1_;
    dz1 = dz1_;
    r1 = r1_;
    FOR_I3 axis[i] = x1_[i] - x0_[i];
    l = MAG(axis);
  }
  
  bool pointIsInside(const double xp[3]) {
    bool inside = false;
    
    double dist_along_axis = (xp[0]-x0[0])/(x1[0]-x0[0]);

    if ((dist_along_axis >= 0.0) && (dist_along_axis <= 1.0)) {

      // determine what the expected dimensions are at this plane perpendicular to the axis
      double dy = dy0 + (dy1-dy0)*dist_along_axis;
      double dz = dz0 + (dz1-dz0)*dist_along_axis;
      double r = r0 + (r1-r0)*dist_along_axis;
      
      // corner radius cannot be bigger than either dy or dz
      if (r > MIN(dy,dz)) r = MIN(dy,dz);

      // is the point inside the two central overlapping rectangles?
      inside =            ((fabs(xp[1]-x0[1]) < dy-r) && (fabs(xp[2]-x0[2]) < dz  ));
      inside = (inside || ((fabs(xp[1]-x0[1]) < dy  ) && (fabs(xp[2]-x0[2]) < dz-r)));

      // is the point inside any of the four rounded corner regions?
      double yc;
      double zc;
      // need to check all 4 quandrants: (1,1), (1,-1), (-1,1) and(-1,-1) 
      FOR_I4 {
	yc = pow(-1.0, floor(i/2))*(dy-r);
	zc = pow(-1.0,i)*(dz-r);
	inside = (inside || ( ((xp[1]-x0[1]-yc)*(xp[1]-x0[1]-yc) + (xp[2]-x0[2]-zc)*(xp[2]-x0[2]-zc)) < r*r ));
      }
      /*
      yc = -1*(dy-r);
      zc = dz-r;
      inside = (inside || ( ((xp[1]-x0[1]-yc)*(xp[1]-x0[1]-yc) + (xp[2]-x0[2]-zc)*(xp[2]-x0[2]-zc)) < r*r ));
      yc = dy-r;
      zc = -1*(dz-r);
      inside = (inside || ( ((xp[1]-x0[1]-yc)*(xp[1]-x0[1]-yc) + (xp[2]-x0[2]-zc)*(xp[2]-x0[2]-zc)) < r*r ));
      yc = -1*(dy-r);
      zc = -1*(dz-r);
      inside = (inside || ( ((xp[1]-x0[1]-yc)*(xp[1]-x0[1]-yc) + (xp[2]-x0[2]-zc)*(xp[2]-x0[2]-zc)) < r*r ));
      */      

    }
    return (inside);
  }
  
  void getBoundingBox(double bbox[6]) {
    // careful of order here...
    bbox[0] = min(x0[0],x1[0]);
    bbox[1] = max(x0[0],x1[0]);
    bbox[2] = min(x0[1]-dy0,x1[1]-dy1);
    bbox[3] = max(x0[1]+dy0,x1[1]+dy1);
    bbox[4] = min(x0[2]-dz0,x1[2]-dz1);
    bbox[5] = max(x0[2]+dz0,x1[2]+dz1);
  }

};

class SimpleDisk : public SimpleGeom {
private:
  
  // the bottom and top, radius...
  double x0[3],x1[3],r;
  // store axis and length
  double axis[3];
  double l;
  
public:
  
  SimpleDisk(const double x0_[3],const double x1_[3],const double r_) {
    FOR_I3 x0[i] = x0_[i];
    FOR_I3 x1[i] = x1_[i];
    r = r_;
    FOR_I3 axis[i] = x1_[i] - x0_[i];
    l = MAG(axis);
  }
  
  bool pointIsInside(const double xp[3]) {
    // project onto axis...
    const double dxp[3] = DIFF(xp,x0);
    const double dp = DOT_PRODUCT(dxp,axis);
    if (dp < 0.0) {
      return false;
    }
    else {
      const double l2 = l*l;
      if (dp > l2) {
        return false;
      }
      else {
        const double my_r2 = DOT_PRODUCT(dxp,dxp) - dp*dp/l2;
        if (my_r2 > r*r) 
          return false;
      }
    }
    return true;
  }

  double pointSignedDistance(const double xp[3]) {
    // compute distance from the point xp[3] to the axis formed between x0[3] -- > x1[3]
    // d = mag( (xp-x0)x(xp-x1) ) / mag(x1-x0)
    const double xpx0[3] = DIFF(xp,x0);
    const double xpx1[3] = DIFF(xp,x1);
    const double numer = CROSS_PRODUCT_MAG(xpx0,xpx1);
    const double radius = numer / l;
    double dist_along_axis = DOT_PRODUCT(xpx0,axis);
    dist_along_axis /= l; // partially normalized: [0,l]

    const double iso_r = radius-r;
    const double iso_0 = -dist_along_axis;
    const double iso_1 = dist_along_axis-l;
    if (dist_along_axis <= 0.0)
      return iso_0;
    else if (dist_along_axis >= l)
      return iso_1;
    else if (radius >= r)
      return iso_r;
    else {
      assert(iso_r < 0);
      assert(iso_0 < 0);
      assert(iso_1 < 0);
      return MAX3(iso_r,iso_0,iso_1);
    }
  }

  void getBoundingBox(double bbox[6]) {
    bbox[0] = min(x0[0],x1[0])-r;
    bbox[1] = max(x0[0],x1[0])+r;
    bbox[2] = min(x0[1],x1[1])-r;
    bbox[3] = max(x0[1],x1[1])+r;
    bbox[4] = min(x0[2],x1[2])-r;
    bbox[5] = max(x0[2],x1[2])+r;
  }

};

class SimpleTruncatedSectorY : public SimpleGeom {
private:

  // assume theta is [0,2*pi] and in radians...
  double xo,zo,theta_min,theta_max,y_min,y_max;

public:
  
  SimpleTruncatedSectorY(const double xo_,const double zo_,const double theta_min_,const double theta_max_,const double y_min_,const double y_max_) {
    xo = xo_;
    zo = zo_;
    theta_min = theta_min_; 
    theta_max = theta_max_;
    y_min = y_min_;
    y_max = y_max_;
  }
  
  bool pointIsInside(const double xp[3]) {
    double theta = atan2(xp[0]-xo,xp[2]-zo);
    if (theta < 0.0) theta += 2.0*M_PI; // make range b/w [0,2*PI]
    return ((xp[1] >= y_min)&&(xp[1] <= y_max)&&(theta >= theta_min)&&(theta <= theta_max));
  }
  
  void getBoundingBox(double bbox[6]) {
    // note that sectors are only semi-bounded...
    // careful of order here...
    if ((cos(theta_min) > 0.0)&&(cos(theta_max) > 0.0))
      bbox[0] = xo;
    else 
      bbox[0] = -HUGE_VAL;
    if ((cos(theta_min) < 0.0)&&(cos(theta_max) < 0.0))
      bbox[1] = xo;
    else 
      bbox[1] = HUGE_VAL;
    bbox[2] = y_min;
    bbox[3] = y_max;
    if ((sin(theta_min) > 0.0)&&(sin(theta_max) > 0.0))
      bbox[4] = zo;
    else 
      bbox[4] = -HUGE_VAL;
    if ((sin(theta_min) < 0.0)&&(sin(theta_max) < 0.0))
      bbox[5] = zo;
    else 
      bbox[5] = HUGE_VAL;
  }

};

class SimpleRotatedRectangle : public SimpleGeom {
private:
 
  int inf_dir; // unbounded direction
  double x0[2]; // bottom  
  double x1[2]; // top 
  double hw; // half-width

  // notes (lexical ordering)
  // inf_dir 0: x0[0] = y0, x0[1] = z0, x1[0] = y1, x1[1] = z1
  // inf_dir 1: x0[0] = x0, x0[1] = z0, x1[0] = x1, x1[1] = z1 
  // inf_dir 2: x0[0] = x0, x0[1] = y0, x1[0] = x1, x1[1] = y1
  
public:
  
  SimpleRotatedRectangle(const int inf_dir_,const double x0_[2],const double x1_[2],const double hw_) {
    inf_dir = inf_dir_; 
    FOR_I2 x0[i] = x0_[i];
    FOR_I2 x1[i] = x1_[i];
    hw = hw_;
  }
  
  bool pointIsInside(const double xp[3]) {
    // project onto axis...
    double dxp[2]; 
    if (inf_dir == 0) {
      dxp[0] = xp[0]-x0[0]; 
      dxp[1] = xp[1]-x0[1];
    }
    else if (inf_dir == 1) {
      dxp[0] = xp[0]-x0[0]; 
      dxp[1] = xp[2]-x0[1];
    }
    else {
      assert(inf_dir == 2);
      dxp[0] = xp[0]-x0[0]; 
      dxp[1] = xp[1]-x0[1];
    }
    const double dx1[2] = DIFF_2D(x1,x0);
    const double dp = DOT_PRODUCT_2D(dxp,dx1);
    if (dp < 0.0) {
      return false;
    }
    else {
      const double dx1_mag2 = DOT_PRODUCT_2D(dx1,dx1);
      if (dp > dx1_mag2) {
        return false;
      }
      else {
        const double my_hw2 = DOT_PRODUCT_2D(dxp,dxp) - dp*dp/dx1_mag2;
        if (my_hw2 > hw*hw) 
          return false;
      }
    }
    return true;
  }

  void getBoundingBox(double bbox[6]) {
    // note that rectangles are only semi-bounded...
    if (inf_dir == 0) {
      bbox[0] = -HUGE_VAL;
      bbox[1] = HUGE_VAL;
      bbox[2] = min(x0[0],x1[0])-hw;
      bbox[3] = max(x0[0],x1[0])+hw;
      bbox[4] = min(x0[1],x1[1])-hw;
      bbox[5] = max(x0[1],x1[1])+hw;
    }
    else if (inf_dir == 1) {
      bbox[0] = min(x0[0],x1[0])-hw;
      bbox[1] = max(x0[0],x1[0])+hw;
      bbox[2] = -HUGE_VAL;
      bbox[3] = HUGE_VAL;
      bbox[4] = min(x0[1],x1[1])-hw;
      bbox[5] = max(x0[1],x1[1])+hw;
    }
    else if (inf_dir == 2) {
      bbox[0] = min(x0[0],x1[0])-hw;
      bbox[1] = max(x0[0],x1[0])+hw;
      bbox[2] = min(x0[1],x1[1])-hw;
      bbox[3] = max(x0[1],x1[1])+hw;
      bbox[4] = -HUGE_VAL;
      bbox[5] = HUGE_VAL;
    }
  }

};

class SimpleHalfPlane : public SimpleGeom {
private:
  
  // point on plane and outward normal...
  double x[3],n[3];
  
public:
  
  SimpleHalfPlane(const double x_[3],const double n_[3]) {
    FOR_I3 x[i] = x_[i];
    FOR_I3 n[i] = n_[i];
  }
  
  bool pointIsInside(const double xp[3]) {
    const double dx[3] = DIFF(xp,x);
    return DOT_PRODUCT(dx,n) < 0.0;
  }

  double pointSignedDistance(const double xp[3]) {
    const double dx[3] = DIFF(xp,x);
    return DOT_PRODUCT(dx,n);
  }
  
  void getBoundingBox(double bbox[6]) {
    FOR_I3 {
      if ((n[(i+1)%3] == 0.0)&&(n[(i+2)%3] == 0.0)) {
        if (n[i] < 0.0) {
          bbox[2*i+0] = -HUGE_VAL;
          bbox[2*i+1] = x[i]*n[i];
        }
        else {
          assert(n[i] > 0.0);
          bbox[2*i+0] = x[i]*n[i];
          bbox[2*i+1] = HUGE_VAL;
        }
      }
      else {
        bbox[2*i+0] = -HUGE_VAL;
        bbox[2*i+1] = HUGE_VAL;
      }
    }
  }

};

// see SimpleGeom.cpp

SimpleGeom * newSimpleGeom(Param * param,int& iarg);

#endif
