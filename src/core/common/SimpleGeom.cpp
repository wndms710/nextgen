#include "SimpleGeom.hpp"

SimpleGeom * newSimpleGeom(Param * param,int& iarg) {

  if (iarg == param->size()) {
    CERR("GEOM missing parameters");
  }
  
  const string token = param->getString(iarg++);
  if ((token == "BOX")||(token == "SIMPLE_BOX")||(token == "CART")) {
    // GEOM BOX x0 x1 y0 y1 z0 z1
    double xmin[3],xmax[3];
    FOR_I3 {
      xmin[i] = param->getDouble(iarg++);
      xmax[i] = param->getDouble(iarg++);
    }
    return new SimpleBox(xmin,xmax);
  }
  else if (token == "SPHERE"){
    // GEOM SPHERE xc yc zc r
    double xc[3],r;
    FOR_I3 xc[i] = param->getDouble(iarg++);
    r = param->getDouble(iarg++);
    if ( r <= 0.0) {
      CERR("SPHERE geometry specified with zero or negative radius - SKIPPING:\n > r:"<< r);
    }
    return new SimpleSphere(xc,r);
  }
  else if (token == "TCONE"){
    // GEOM TCONE x0 y0 z0 r0 x1 y1 z1 r1
    double x0[3],r0;
    double x1[3],r1;
    FOR_I3 x0[i] = param->getDouble(iarg++);
    r0 = param->getDouble(iarg++);
    if ( r0 <= 0.0) {
      CERR("Conical geometry specified with zero or negative radius - SKIPPING:\n > r0:"<< r0);
    }
    FOR_I3 x1[i] = param->getDouble(iarg++);
    r1 = param->getDouble(iarg++);
    if ( r1 <= 0.0) {
      CERR("Conical geometry specified with zero or negative radius - SKIPPING:\n > r1:"<< r1);
    }
    const double axis[3] = DIFF(x1,x0);
    const double l = MAG(axis);
    if (l < 1e-12) {
      CERR("Conical geometry specified with collocated end points - SKIPPING:\n > x0,x1:"<< COUT_VEC(x0) << " " << COUT_VEC(x1));
    }
    return new SimpleTcone(x0,r0,x1,r1);
  }
  else if (token == "ROUNDED_RECT_X"){
    // GEOM ROUNDED_RECT_X x0 y0 z0 dy0 dz0 r0 x1 dy1 dz1 r1
    // this geometry is useful for defining FWH surfaces for rectangular jets
    // 0 and 1 refer either end of the (rounded) rectangular prism (oriented in x-dir)
    // the end sections are centered at (x0,y0,z0) and (x1,y0,z0), i.e., x-aligned
    // dy0, dz0, dy1, dz1 are half the lengths of the edges of the rectangle at either ends
    // r0, r1 are the radii of the rounded corners at either ends
    double x0[3],x1[3];
    double dy0,dz0,r0,dy1,dz1,r1;
    // front end
    FOR_I3 x0[i] = param->getDouble(iarg++);
    dy0 = param->getDouble(iarg++);
    if ( dy0 <= 0.0) {
      CERR("ROUNDED_RECT_X geometry specified with zero or negative half-width - SKIPPING:\n > dy0:"<< dy0);
    }
    dz0 = param->getDouble(iarg++);
    if ( dz0 <= 0.0) {
      CERR("ROUNDED_RECT_X geometry specified with zero or negative half-width - SKIPPING:\n > dz0:"<< dz0);
    }
    r0 = param->getDouble(iarg++);
    if ( r0 < 0.0) {
      CERR("ROUNDED_RECT_X geometry specified with negative radius - SKIPPING:\n > r0:"<< r0);
    }
    // back end
    x1[0] = param->getDouble(iarg++);
    x1[1] = x0[1];
    x1[2] = x0[2];
    dy1 = param->getDouble(iarg++);
    if ( dy1 <= 0.0) {
      CERR("ROUNDED_RECT_X geometry specified with zero or negative half-width - SKIPPING:\n > dy1:"<< dy1);
    }
    dz1 = param->getDouble(iarg++);
    if ( dz1 <= 0.0) {
      CERR("ROUNDED_RECT_X geometry specified with zero or negative half-width - SKIPPING:\n > dz1:"<< dz1);
    }
    r1 = param->getDouble(iarg++);
    if ( r1 < 0.0) {
      CERR("ROUNDED_RECT_X geometry specified with negative radius - SKIPPING:\n > r1:"<< r1);
    }
    const double axis[3] = DIFF(x1,x0);
    const double l = MAG(axis);
    if (l < 1e-12) {
      CERR("ROUNDED_RECT_X geometry specified with collocated end points - SKIPPING:\n > x0,x1:"<< COUT_VEC(x0) << " " << COUT_VEC(x1));
    }
    return new SimpleRoundedRect_X(x0,dy0,dz0,r0,x1,dy1,dz1,r1);
  }
  else if ((token == "DISK")||(token == "SIMPLE_DISK")||(token == "CYL")||(token == "CYLINDER")) {
    // GEOM DISK x0 y0 z0 x1 y1 z1 r
    double x0[3],x1[3],r;
    FOR_I3 x0[i] = param->getDouble(iarg++);
    FOR_I3 x1[i] = param->getDouble(iarg++);
    r = param->getDouble(iarg++);
    return new SimpleDisk(x0,x1,r);
  }
  else if ((token == "TSECTOR_Y")||(token == "TRUNCATED_SECTOR_Y")) {
    // GEOM TSECTOR_Y xo zo theta_min theta_max y_min y_max 
    const double xo        = param->getDouble(iarg++);
    const double zo        = param->getDouble(iarg++);
    const double theta_min = param->getDouble(iarg++);
    const double theta_max = param->getDouble(iarg++);
    const double y_min     = param->getDouble(iarg++);
    const double y_max     = param->getDouble(iarg++);
    return new SimpleTruncatedSectorY(xo,zo,theta_min,theta_max,y_min,y_max);
  }
  else if (token == "ROTATED_RECTANGLE") {
    // GEOM ROTATED_RECTANGLE inf_dir x0[0] x0[1] x1[0] x1[1] width
    // GEOM ROTATED_RECTANGLE 0 y0 z0 y1 z1 half-width
    // GEOM ROTATED_RECTANGLE 1 x0 z0 x1 z1 half-width
    // GEOM ROTATED_RECTANGLE 2 x0 y0 x1 y1 half-width
    const int inf_dir = param->getInt(iarg++);
    double x0[2]; FOR_I2 x0[i] = param->getDouble(iarg++);
    double x1[2]; FOR_I2 x1[i] = param->getDouble(iarg++);
    const double hw = param->getDouble(iarg++);
    return new SimpleRotatedRectangle(inf_dir,x0,x1,hw);
  }
  else if (token == "PLANE") {
    // GEOM PLANE x y z nx ny nz
    double x[3],n[3]; 
    FOR_I3 x[i] = param->getDouble(iarg++); 
    FOR_I3 n[i] = param->getDouble(iarg++); 
    return new SimpleHalfPlane(x,n);
  }
  else if (token == "PLANE_X") {
    // GEOM PLANE_X x
    double x[3],n[3];
    x[0] = param->getDouble(iarg++); 
    x[1] = 0;
    x[2] = 0;
    n[0] = 1.0;
    n[1] = 0;
    n[2] = 0;
    return new SimpleHalfPlane(x,n);
  }
  else if (token == "PLANE_Y") {
    // GEOM PLANE_Y y
    double x[3],n[3];
    x[0] = 0;                        
    x[1] = param->getDouble(iarg++); 
    x[2] = 0;                        
    n[0] = 0;                        
    n[1] = 1.0;
    n[2] = 0;                        
    return new SimpleHalfPlane(x,n);
  }
  else if (token == "PLANE_Z") {
    // GEOM PLANE_Z z
    double x[3],n[3];
    x[0] = 0;
    x[1] = 0;
    x[2] = param->getDouble(iarg++);
    n[0] = 0;
    n[1] = 0;
    n[2] = 1.0;
    return new SimpleHalfPlane(x,n);
  }
  else {
    CWARN("unrecognized GEOM: " << token << ". Available options are:\n" <<
          "BOX <x0> <x1> <y0> <y1> <z0> <z1>\n" <<
          "SPHERE <x> <y> <z> <r>\n" <<
          "TCONE <x0> <y0> <z0> <r0> <x1> <y1> <z1> <r1>\n" <<
          "TSECTOR_Y <xo> <zo> <theta_min> <theta_max> <y_min> <y_max>\n" << 
          "DISK <x0> <y0> <z0> <x1> <y1> <z1> <r> (or CYLINDER)\n"
          "ROUNDED_RECT_X <x0> <y0> <z0> <dy0> <dz0> <r0> <x1> <dy1> <dz1> <r1>\n" <<
          "PLANE <x> <y> <z> <nx> <ny> <nz>");
    return NULL;
  }  
}
