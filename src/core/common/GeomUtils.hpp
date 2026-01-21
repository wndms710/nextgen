#ifndef _GEOM_UTILS_HPP_
#define _GEOM_UTILS_HPP_

#include "Common.hpp"
#include "SimpleTri.hpp"
#include "NewTri.hpp"

namespace GeomUtils {
  
  // series of low-level geom utils. Should not be parallel, and should not rely on complex
  // types or classes...

  // returns the distance starting from 0 of grid location i for a grid with first 
  // spacing dn and last spacing dt, n spaces...

  double stretchingFunction(const int i,const int n,const double dn,const double dt);
    
  // move this over eventually...
  //MiscUtils::createCirclePts(xsp,ss_nsp0+1,x0,axis,r0,n);

  void createCirclePts(double (*xp)[3],const int n,const double xc[3],const double axis[3],const double r); // note that axis need not be normalized
  void facetCircleToPoint(int (*spost)[3], int *znost, const int indexPt, const int indexCircle, const int n, const int zoneId, const bool flip);
  void facetCircleToCircle(int (*spost)[3], int *znost, const int indexC0, const int indexC1, const int n, const int zoneId, const bool closed, const bool flip);

  // these functions return the new facet count
  int facetGap(int (*spost)[3],const vector<int> isp0Vec,const vector<int> isp1Vec,const double (*const xsp)[3],const bool loop);
  int facetGap(int (*spost)[3],const int isp0_f,const int isp0_l,const int isp1_f,const int isp1_l,const double (*const xsp)[3]);

  void writeTecplot(const string& filename,const int (*const spost)[3],const int nst,const double (*const xsp)[3]);
  void writeTecplotZone(FILE * fp,const int (* const spost)[3],const int nst,const double (* const xsp)[3],const int * const flag,const int flag_value);
  void writeTecplot(const string& filename,const int (* const spost)[3],const int nst,const double (* const xsp)[3],const set<pair<int,int> >& istBitsSet);
  void writeTecplotEdges(const string& filename,const double (* const xsp)[3],const vector<pair<int,int> >& linkVec);
  void writeSbin(const string& filename,const int (*const spost)[3],const int nst,const double (*const xsp)[3]);

  // writing of spost triangulated surface...
  void writeTecplot(const string& filename,const int (*const spost)[3],const int nst,const double (*const xsp)[3]);
  void writeSbin(const string& filename,const int (*const spost)[3],const int nst,const double (*const xsp)[3]);

  // writing of triVec to tecplot...
  void writeTecplot(const string& filename,vector<pair<SimpleTri,int> >& triVec);

  // NewTri to tecplot file...
  void writeTecplotGeom(const string& filename,vector<NewTri> &newTriVec,const double (*const xsp)[3]);
    
  // collective writing of distributed points using handshake...
  void writePtsTecplot(const string& filename,const double (*const xp)[3],const int np);

  // n is the number of azimuthal segments...
  int getAnnularCylNodeCount(const int nseg);
  int getAnnularCylTriCount(const int nseg);
  void addAnnularCyl(double (*xsp)[3],int (*spost)[3],const double x0[3],const double x1[3],const double r0,const double r1,const int nseg,const bool flip); // flip=false: fluid inside, flip=true: fluid outside

  // reconcile these with the above some day...
  // the fact these return the tri count is useful...
  int facetCircleToPointUseThisOne(int (* spost)[3], const int indexPt, const int indexCircle, const int n, const bool flip);
  int facetGapUseThisOne(int (*spost)[3], const int isp0_start,const int nsp0,const int isp1_start, const int nsp1,const double (* const xsp)[3], const bool loop);

  // hemisphere...
  void getHemisphereNodeAndTriCount(int &nsp,int &nst,const int ntheta);
  void addHemisphere(double (*xsp)[3],int (*spost)[3],const double xp[3],const double np[3],const double rp,const int ntheta,const bool flip);

  // simle ASCII file reading...
  int readXYZ(vector<double>& dVec,const string& filename);
  int readXY(vector<double>& dVec,const string& filename);

  // geometry routines for bloatedTruiIntersections...
  bool getSphereIntersections(double xmp[2],const double yz[2],const double x[3],const double r);
  bool getMinusSphereIntersection(double &xm,const double yz[2],const double x[3],const double r);
  bool getPlusSphereIntersection(double &xp,const double yz[2],const double x[3],const double r);
  bool calcCappedCylinderIntersections(double xmp[2],const double yz[2],const double x0[3],const double x1[3],const double r);
  int getMinusCylinderIntersection(double &xm,const double yz[2],const double x0[3],const double x1[3],const double r);
  int getPlusCylinderIntersection(double &xp,const double yz[2],const double x0[3],const double x1[3],const double r);
  bool calcBloatedTriIntersections(double xmp[2],const double yz[2],const double x0[3],const double x1[3],const double x2[3],const double delta);
  
  // 2d marching front stuff: no new nodes added...
  bool addTrisMarchingFront2d(vector<NewTri> &newTriVec,set<pair<int,int> > &frontSet,const double (* const x)[2],const int n);

  // coordinate transform...
  enum COORD_TRANSFORM_TYPE {
    TRANSFORM_CYL_X = 1,
    TRANSFORM_CYL_Y = 2,
    TRANSFORM_CYL_Z = 3,
    NO_TRANSFORM = 0
  };
  
  void coordinateTransform(double x_t[3],const double x[3],const double r0,const double cos_theta0,const double sin_theta0,
			   const double cyl_center[3],const COORD_TRANSFORM_TYPE transform_type);

}

#endif
