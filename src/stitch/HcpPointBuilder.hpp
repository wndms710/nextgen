#ifndef _HCP_POINT_BUILDER_NEW_HPP_
#define _HCP_POINT_BUILDER_NEW_HPP_

#include "MiscUtils.hpp"
#include "KillfileReader.hpp"
#include "WriteImageData.hpp"
#include "CtiScene.hpp"
#include "WebUI.hpp"
#include "PartData.hpp"
#include "GeomUtils.hpp"

// for more precise timing, add "-std=c++11" to your CXXFLAGS
// and uncomment the chrono stuff
#ifdef WITH_CHRONO
#include <chrono>
#endif

//#define HCP_BIT_MAX 28      // do not go above 30
// this is in Part.hpp now
//#define LEVEL_MASK_BITS 31  // max level < 2^5-1 (reserve 31 for part points)
#define HCP_WINDOW_SHIFT_BITS 5 // window < 2^(16-5-1)=1023
#define MAX_HCP_WINDOWS 1023

#define HCP_PACKING_ONE   0 // most uniform faces (singularity prone)
#define HCP_PACKING_ROOT2 1 // most spherical
#define HCP_PACKING_ROOT3 2 // similar to root2, but stretched in one-direction
#define HCP_PACKING_CART  3 // simple cartesian grid excluding transitions

enum HcpWindowKind {
  UNKNOWN_HCP_WINDOW,
  FAZONE_HCP_WINDOW,
  BOX_HCP_WINDOW,
  TCONE_HCP_WINDOW,
  ELLIPSE_TCONE_HCP_WINDOW,
  ANNULAR_TCONE_HCP_WINDOW,
  ANNULAR_ELLIPSE_TCONE_HCP_WINDOW,
  PRISM_HCP_WINDOW,
  HEMISPHERE_HCP_WINDOW,
  SBIN_HCP_WINDOW,
};

class HcpWindow {
public:
  int iwow; // index into window vec
  string hash; // hash is used for geom windows only
  string param_str; //used by webui to display active windows
  HcpWindowKind kind;
  //string zoneNames;
  vector<int> zoneIdVec;
  bool b_dx;
  double dx;
  int level;
  bool b_d;
  double d;
  bool b_n;
  int n;
  bool b_nlayers;
  string nlayersString;
  bool b_dlayers;
  string dlayersString;
  double geom_data[20];
  string sbin_filename;
  HcpWindow() {
    iwow = -1;
    kind = UNKNOWN_HCP_WINDOW;
    b_dx = false;
    b_d = false;
    b_n = false;
    b_nlayers = false;
    b_dlayers = false;
  }
  string getKindStr() const {
    switch (kind) {
    case UNKNOWN_HCP_WINDOW:
      return "UNKNOWN_HCP_WINDOW";
    case FAZONE_HCP_WINDOW:
      return "FAZONE_HCP_WINDOW";
    case BOX_HCP_WINDOW:
      return "BOX_HCP_WINDOW";
    case TCONE_HCP_WINDOW:
      return "TCONE_HCP_WINDOW";
    case ELLIPSE_TCONE_HCP_WINDOW:
      return "ELLIPSE_TCONE_HCP_WINDOW";
    case ANNULAR_TCONE_HCP_WINDOW:
      return "ANNULAR_TCONE_HCP_WINDOW";
    case ANNULAR_ELLIPSE_TCONE_HCP_WINDOW:
      return "ANNULAR_ELLIPSE_TCONE_HCP_WINDOW";
    case PRISM_HCP_WINDOW:
      return "PRISM_HCP_WINDOW";
    case HEMISPHERE_HCP_WINDOW:
      return "HEMISPHERE_HCP_WINDOW";
    case SBIN_HCP_WINDOW:
      return "SBIN_HCP_WINDOW";
    default:
      return "kind not specified in HcpWindow::getKindStr()";
    }
  }
  string str() const {
    std::stringstream ss;
    switch (kind) {
    case UNKNOWN_HCP_WINDOW:
      ss << "UNKNOWN_HCP_WINDOW";
      break;
    case FAZONE_HCP_WINDOW:
      ss << "HCP_WINDOW FAZONE ";
      for (int ii = 0; ii < zoneIdVec.size(); ++ii)
	ss << zoneIdVec[ii] << ",";
      break;
    default:
      ss << "str() not implemented for this kind: " << kind;
    }
    if (b_dx) {
      ss << " DX " << dx;
    }
    if (b_d) {
      ss << " D " << d;
    }
    if (b_n) {
      ss << " N " << n;
    }
    return ss.str();
  }
};

class HcpVertex {
public:
  double x[3];
  int level;
  int window;
  HcpVertex() {}
  HcpVertex(const double x[3],const int level,const int window) {
    this->x[0]   = x[0];
    this->x[1]   = x[1];
    this->x[2]   = x[2];
    this->level  = level;
    this->window = window;
  }
};

class InOutData {
public:
  double xp;
  uint8 index;
  int level;
  InOutData() {}
  InOutData(const double xp,const uint8 index,const int level) {
    this->xp = xp;
    this->index = index;
    this->level = level;
  }
  bool operator<(const InOutData& rhs) const {
    return (index < rhs.index) || ((index == rhs.index)&&(xp < rhs.xp));
  }
};

class TriStuff {
public:

  class Tri {
  public:

    double x[3][3];

    int window;

    void init(const double x0[3],const double x1[3],const double x2[3],const int window) {
      FOR_I3 this->x[0][i] = x0[i];
      FOR_I3 this->x[1][i] = x1[i];
      FOR_I3 this->x[2][i] = x2[i];
      this->window  = window;
    }

    void setCentroid(double xc[3]) const {
      FOR_I3 xc[i] = (x[0][i] + x[1][i] + x[2][i])/3.0;
    }

  };

  static void initSbin(vector<Tri>& triVec,const string& sbin_filename,const int window) {

    FILE * fp = NULL;
    const int file_err = MiscUtils::openFile(&fp,sbin_filename,"rb");
    if (file_err != 0) CERR(" > HCP_WINDOW SBIN file DNE");

    bool byte_swap = false;
    int version;
    fread(&version,sizeof(int),1,fp);
    if ((version < 1)||(version > 10)) {
      version = ByteSwap::byteSwap(version);
      if ((version < 1)||(version > 10)) {
        CERRT("byteswap does not produce reasonable version: " << version,20);
      }
      byte_swap = true;
    }
    assert(version == 1);

    int count;
    fread(&count,sizeof(int),1,fp);
    if (byte_swap)
      count = ByteSwap::byteSwap(count);

    for (int izone = 0; izone < count; ++izone) {
      int length;
      fread(&length,sizeof(int),1,fp);
      if (byte_swap)
        length = ByteSwap::byteSwap(length);
      char cbuf[length+1];
      fread(cbuf,sizeof(char),length,fp); // no byteswap for char
    }

    int nsp;
    fread(&nsp,sizeof(int),1,fp);
    if (byte_swap)
      nsp = ByteSwap::byteSwap(nsp);
    double (*xp)[3] = new double[nsp][3];
    fread(xp,sizeof(double),nsp*3,fp);
    if (byte_swap)
      ByteSwap::byteSwap(xp,nsp);
    int nst;
    fread(&nst,sizeof(int),1,fp);
    if (byte_swap)
      nst = ByteSwap::byteSwap(nst);
    int (*spost)[3] = new int[nst][3];
    fread(spost,sizeof(int),nst*3,fp);
    if (byte_swap)
      ByteSwap::byteSwap(spost,nst);
    fclose(fp);

    // push into triVec...
    const int nst0 = triVec.size();
    triVec.resize(nst0+nst);
    for (int ist = 0; ist < nst; ++ist) {
      const double * const x0 = xp[spost[ist][0]];
      const double * const x1 = xp[spost[ist][1]];
      const double * const x2 = xp[spost[ist][2]];
      triVec[nst0+ist].init(x0,x1,x2,window);
    }

    // cleanup...
    delete[] xp;
    delete[] spost;

  }

  static void initBox(vector<Tri>& triVec,const double x0,const double x1,
		      const double y0,const double y1,
		      const double z0,const double z1,const int window) {


    //cout << "initBox: " << x0 << " " << x1 << " " << y0 << " " << y1 << " " << z0 << " " << z1 << endl;

    // add 12 tris describing a box...
    const int nst0 = triVec.size();
    triVec.resize(nst0+12);

    const double x000[3] = { x0, y0, z0 };
    const double x001[3] = { x0, y0, z1 };
    const double x010[3] = { x0, y1, z0 };
    const double x011[3] = { x0, y1, z1 };
    const double x100[3] = { x1, y0, z0 };
    const double x101[3] = { x1, y0, z1 };
    const double x110[3] = { x1, y1, z0 };
    const double x111[3] = { x1, y1, z1 };

    int itri = nst0;
    triVec[itri++].init(x000,x001,x011,window);
    triVec[itri++].init(x000,x011,x010,window);
    triVec[itri++].init(x000,x100,x101,window);
    triVec[itri++].init(x000,x101,x001,window);
    triVec[itri++].init(x001,x101,x111,window);
    triVec[itri++].init(x001,x111,x011,window);
    triVec[itri++].init(x011,x111,x110,window);
    triVec[itri++].init(x011,x110,x010,window);
    triVec[itri++].init(x000,x010,x110,window);
    triVec[itri++].init(x000,x110,x100,window);
    triVec[itri++].init(x100,x110,x111,window);
    triVec[itri++].init(x100,x111,x101,window);
    assert(itri == nst0+12);

  }

  static void initPrism(vector<Tri>& triVec,
                        const double x0,const double y0,const double z0,const double rmaj0,const double rmin0,
                        const double x1,const double y1,const double z1,const double rmaj1,const double rmin1,
                        const double m_axis0,const double m_axis1,const double m_axis2,
                        const double nx,const double ny,const double nz,const int window) {


    const double maj_axis[3] = {m_axis0,m_axis1,m_axis2};
    const double xc0[3] = { x0, y0, z0 };
    const double xc1[3] = { x1, y1, z1 };

    // absolute value ensures no folding as well as positive volume
    const double rM0 = fabs(0.5*rmaj0);
    const double rm0 = fabs(0.5*rmin0);
    const double rM1 = fabs(0.5*rmaj1);
    const double rm1 = fabs(0.5*rmin1);

    const double SMALL = 1E-12;
    if (rm0 < SMALL || rm1 < SMALL || rM0 < SMALL || rM1 < SMALL) CERR("cannot create prism with either major/minor axis length = 0");

    double axis[3];
    if (nx+ny+nz != 0.0) {
      axis[0] = nx;
      axis[1] = ny;
      axis[2] = nz;
    }
    else {
      FOR_I3 axis[i] = xc1[i] - xc0[i];
    }

    double radDir1[3];
    double radDir2[3];
    MiscUtils::getOrthogonalVectors(radDir1,radDir2,xc0,axis,maj_axis);

    double x000[3];
    FOR_I3 x000[i] = xc0[i] - rM0*radDir1[i] - rm0*radDir2[i];
    double x001[3];
    FOR_I3 x001[i] = xc0[i] - rM0*radDir1[i] + rm0*radDir2[i];
    double x010[3];
    FOR_I3 x010[i] = xc0[i] + rM0*radDir1[i] - rm0*radDir2[i];
    double x011[3];
    FOR_I3 x011[i] = xc0[i] + rM0*radDir1[i] + rm0*radDir2[i];

    double x100[3];
    FOR_I3 x100[i] = xc1[i] - rM1*radDir1[i] - rm1*radDir2[i];
    double x101[3];
    FOR_I3 x101[i] = xc1[i] - rM1*radDir1[i] + rm1*radDir2[i];
    double x110[3];
    FOR_I3 x110[i] = xc1[i] + rM1*radDir1[i] - rm1*radDir2[i];
    double x111[3];
    FOR_I3 x111[i] = xc1[i] + rM1*radDir1[i] + rm1*radDir2[i];

    // add 12 tris describing a box...
    const int nst0 = triVec.size();
    triVec.resize(nst0+12);

    int itri = nst0;
    triVec[itri++].init(x000,x001,x011,window);
    triVec[itri++].init(x000,x011,x010,window);
    triVec[itri++].init(x000,x100,x101,window);
    triVec[itri++].init(x000,x101,x001,window);
    triVec[itri++].init(x001,x101,x111,window);
    triVec[itri++].init(x001,x111,x011,window);
    triVec[itri++].init(x011,x111,x110,window);
    triVec[itri++].init(x011,x110,x010,window);
    triVec[itri++].init(x000,x010,x110,window);
    triVec[itri++].init(x000,x110,x100,window);
    triVec[itri++].init(x100,x110,x111,window);
    triVec[itri++].init(x100,x111,x101,window);
    assert(itri == nst0+12);

  }

  static void initTcone(vector<Tri>& triVec,
                        const double x0,const double y0,const double z0,const double rad0,
                        const double x1,const double y1,const double z1,const double rad1,
                        const double nx,const double ny,const double nz,const int window,const int n=128) {

    const double xc0[3] = { x0, y0, z0 };
    const double xc1[3] = { x1, y1, z1 };

    // absolute value ensures no singularity in the middle of the tcone as well as positive volume
    const double r0 = fabs(rad0);
    const double r1 = fabs(rad1);

    const double SMALL = 1E-12;
    if (r0 < SMALL || r1 < SMALL) CERR("cannot create a truncated cone with either radius = 0");

    double axis[3];
    if (nx+ny+nz != 0.0) {
      axis[0] = nx;
      axis[1] = ny;
      axis[2] = nz;
    }
    else {
      FOR_I3 axis[i] = xc1[i] - xc0[i];
    }

    const int nsp = 2*n+2;
    double (*xp_tmp)[3] = new double[nsp][3];  // holds tri points
    // first and last node are cap centers
    FOR_I3 {
      xp_tmp[0][i] = xc0[i];
      xp_tmp[nsp - 1][i] = xc1[i];
    }

    MiscUtils::createCirclePts(xp_tmp,1,xc0,axis,r0,n);
    MiscUtils::createCirclePts(xp_tmp,n+1,xc1,axis,r1,n);

    // now create tris
    const int nst0 = triVec.size();
    triVec.resize(nst0+4*n);

    // cap0, cap1
    facetCircleToPoint(triVec,xp_tmp,0,1,nst0,n,window);
    facetCircleToPoint(triVec,xp_tmp,nsp-1,n+1,nst0+3*n,n,window,true);

    // wall
    facetCircleToCircle(triVec,xp_tmp,1,n+1,nst0+n,n,window,true);

    DELETE(xp_tmp);

  }

  static void initAnnularTcone(vector<Tri>& triVec,
			       const double x0,const double y0,const double z0,const double _r00,const double _r01,
			       const double x1,const double y1,const double z1,const double _r10,const double _r11,
                               const double nx,const double ny,const double nz,const int window,const int n=128) {

    const double xc0[3] = { x0, y0, z0 };
    const double xc1[3] = { x1, y1, z1 };

    // absolute value ensures no singularity in the middle of the tcone as well as positive volume
    // ensure radii are ordered min/max
    const double r00 = min(fabs(_r00),fabs(_r01));
    const double r01 = max(fabs(_r00),fabs(_r01));
    const double r10 = min(fabs(_r10),fabs(_r11));
    const double r11 = max(fabs(_r10),fabs(_r11));

    const double SMALL = 1E-12;
    if (r00 < SMALL || r01 < SMALL || r10 < SMALL || r11 < SMALL) 
      CERR("cannot create a annular truncated cone with any radius = 0");

    double axis[3];
    if (nx+ny+nz != 0.0) {
      axis[0] = nx;
      axis[1] = ny;
      axis[2] = nz;
    }
    else {
      FOR_I3 axis[i] = xc1[i] - xc0[i];
    }

    const int nsp = 4*n;
    double (*xp_tmp)[3] = new double[nsp][3];  // holds tri points

    MiscUtils::createCirclePts(xp_tmp,0,xc0,axis,r00,n);
    MiscUtils::createCirclePts(xp_tmp,n,xc0,axis,r01,n);
    MiscUtils::createCirclePts(xp_tmp,2*n,xc1,axis,r10,n);
    MiscUtils::createCirclePts(xp_tmp,3*n,xc1,axis,r11,n);

    // now create tris
    const int nst0 = triVec.size();
    triVec.resize(nst0+8*n);

    facetCircleToCircle(triVec,xp_tmp,0,n,nst0,n,window,true);
    facetCircleToCircle(triVec,xp_tmp,n,3*n,nst0+2*n,n,window,true);
    facetCircleToCircle(triVec,xp_tmp,2*n,0,nst0+4*n,n,window,true);
    facetCircleToCircle(triVec,xp_tmp,3*n,2*n,nst0+6*n,n,window,true);

    DELETE(xp_tmp);

  }

  static void initEllipseTcone(vector<Tri>& triVec,
                               const double x0,const double y0,const double z0,const double _r0M,const double _r0m,
                               const double x1,const double y1,const double z1,const double _r1M,const double _r1m,
                               const double m_axis0,const double m_axis1,const double m_axis2,
                               const double nx,const double ny,const double nz,const int window,const int n=128) {


    const double maj_axis[3] = {m_axis0,m_axis1,m_axis2};

    const double xc0[3] = { x0, y0, z0 };
    const double xc1[3] = { x1, y1, z1 };

    // ensures no infinitely thin ellipse
    const double r0M = fabs(_r0M);
    const double r0m = fabs(_r0m);
    const double r1M = fabs(_r1M);
    const double r1m = fabs(_r1m);

    const double SMALL = 1E-12;
    if (r0M < SMALL || r0m < SMALL || r1M < SMALL || r1m < SMALL) 
      CERR("cannot create an ellipse truncated cone with any radius = 0");

    double axis[3];
    if (nx+ny+nz != 0.0) {
      axis[0] = nx;
      axis[1] = ny;
      axis[2] = nz;
    }
    else {
      FOR_I3 axis[i] = xc1[i] - xc0[i];
    }

    const int nsp = 2*n+2;
    double (*xp_tmp)[3] = new double[nsp][3];  // holds tri points
    // first and last node are cap centers
    FOR_I3 {
      xp_tmp[0][i] = xc0[i];
      xp_tmp[nsp - 1][i] = xc1[i];
    }

    createEllipsePts(xp_tmp,1,xc0,axis,r0M,r0m,maj_axis,n);
    createEllipsePts(xp_tmp,n+1,xc1,axis,r1M,r1m,maj_axis,n);

    // now create tris
    const int nst0 = triVec.size();
    triVec.resize(nst0+4*n);

    // cap0, cap1
    facetCircleToPoint(triVec,xp_tmp,0,1,nst0,n,window);
    facetCircleToPoint(triVec,xp_tmp,nsp-1,n+1,nst0+3*n,n,window,true);

    // wall
    facetCircleToCircle(triVec,xp_tmp,1,n+1,nst0+n,n,window,true);

    DELETE(xp_tmp);

  }

  static void initAnnularEllipseTcone(vector<Tri>& triVec,
                                      const double x0,const double y0,const double z0,
                                      const double _r0M0,const double _r0m0,const double _r0M1,const double _r0m1,
                                      const double x1,const double y1,const double z1,
                                      const double _r1M0,const double _r1m0,const double _r1M1,const double _r1m1,
                                      const double m_axis0,const double m_axis1,const double m_axis2,
                                      const double nx,const double ny,const double nz,const int window,const int n=128) {


    const double maj_axis[3] = {m_axis0,m_axis1,m_axis2};

    const double xc0[3] = { x0, y0, z0 };
    const double xc1[3] = { x1, y1, z1 };

    // ensures no infinitely thin ellipse
    // absolute value ensures no singularity in the middle of the tcone as well as positive volume
    // ensure radii are ordered min/max
    const double r0M0 = min(fabs(_r0M0),fabs(_r0M1));
    const double r0m0 = min(fabs(_r0m0),fabs(_r0m1));
    const double r0M1 = max(fabs(_r0M0),fabs(_r0M1));
    const double r0m1 = max(fabs(_r0m0),fabs(_r0m1));
    const double r1M0 = min(fabs(_r1M0),fabs(_r1M1));
    const double r1m0 = min(fabs(_r1m0),fabs(_r1m1));
    const double r1M1 = max(fabs(_r1M0),fabs(_r1M1));
    const double r1m1 = max(fabs(_r1m0),fabs(_r1m1));

    const double SMALL = 1E-12;
    if (r0M0 < SMALL || r0m0 < SMALL || r1M0 < SMALL || r1m0 < SMALL ||
        r0M1 < SMALL || r0m1 < SMALL || r1M1 < SMALL || r1m1 < SMALL) {
      CERR("cannot create an annular ellipse truncated tcone with any radius = 0");
    }

    double axis[3];
    if (nx+ny+nz != 0.0) {
      axis[0] = nx;
      axis[1] = ny;
      axis[2] = nz;
    }
    else {
      FOR_I3 axis[i] = xc1[i] - xc0[i];
    }

    const int nsp = 4*n;
    double (*xp_tmp)[3] = new double[nsp][3];  // holds tri points

    MiscUtils::createEllipsePts(xp_tmp,0,xc0,axis,r0M0,r0m0,maj_axis,n);
    MiscUtils::createEllipsePts(xp_tmp,n,xc0,axis,r0M1,r0M1,maj_axis,n);
    MiscUtils::createEllipsePts(xp_tmp,2*n,xc1,axis,r1M0,r1m0,maj_axis,n);
    MiscUtils::createEllipsePts(xp_tmp,3*n,xc1,axis,r1M1,r1m1,maj_axis,n);

    // now create tris
    const int nst0 = triVec.size();
    triVec.resize(nst0+8*n);

    facetCircleToCircle(triVec,xp_tmp,0,n,nst0,n,window,true);
    facetCircleToCircle(triVec,xp_tmp,n,3*n,nst0+2*n,n,window,true);
    facetCircleToCircle(triVec,xp_tmp,2*n,0,nst0+4*n,n,window,true);
    facetCircleToCircle(triVec,xp_tmp,3*n,2*n,nst0+6*n,n,window,true);

    DELETE(xp_tmp);

  }

  static void initHemisphere(vector<Tri>& triVec,const double xc,const double yc,const double zc,
                             const double nx,const double ny,const double nz,const double r,const int window) {

    const double x[3] = { xc, yc, zc };
    const double n[3] = { nx, ny, nz };

    const int ntheta = 128; // could pass this in, but 128 probably ok for window

    int nsp,nst;
    GeomUtils::getHemisphereNodeAndTriCount(nsp,nst,ntheta);

    double (*xsp)[3] = new double[nsp][3];
    int (*spost)[3] = new int[nst][3];
    GeomUtils::addHemisphere(xsp,spost,x,n,r,ntheta,false);

    const int ntri = triVec.size();
    triVec.resize(ntri+nst);

    for (int ist = 0; ist < nst; ++ist) {
      triVec[ntri+ist].init(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],window);
    }

    delete[] xsp;
    delete[] spost;

  }

  // TODO: can this stuff be put in GeomUtils?...

  /*
   * create the facets between a circle of points and a single point
   */
  static void facetCircleToPoint(vector<Tri>& triVec, const double (*xp_tmp)[3], const int indexPt, const int indexCircle, const int st0, const int n,
                                 const int window,const bool flip = false) {
    // flip dictates node ordering (normal direction convention)
    int pt1_inc = 1;
    int pt2_inc = 0;
    if (flip) {
      pt1_inc = 0;
      pt2_inc = 1;
    }

    for (int i=0, limit=n; i < limit; ++i) {
      triVec[st0+i].init(xp_tmp[indexPt],
                         xp_tmp[indexCircle + (i+pt1_inc)%n],
                         xp_tmp[indexCircle + (i+pt2_inc)%n],
                         window);
    }
  }

  /*
   * create the facets between two circles of points, assuming they have the same n
   */
  static void facetCircleToCircle(vector<Tri>& triVec, const double (*xp_tmp)[3], const int indexC0, const int indexC1, const int st0, const int n,
                                  const int window,const bool flip = false) {
    // flip dictates node ordering (normal direction convention)
    int pt1_inc = 1;
    int pt2_inc = 0;
    if (flip) {
      pt1_inc = 0;
      pt2_inc = 1;
    }

    for (int i=0, limit=n; i < limit; ++i) {
      triVec[st0+(2*i)].init(xp_tmp[indexC0 + (i+pt1_inc)%n],
                             xp_tmp[indexC0 + (i+pt2_inc)%n],
                             xp_tmp[indexC1 + (i+pt1_inc)%n],
                             window);

      triVec[st0+(2*i + 1)].init(xp_tmp[indexC0 + (i+pt2_inc)%n],
                                 xp_tmp[indexC1 + (i+pt2_inc)%n],
                                 xp_tmp[indexC1 + (i+pt1_inc)%n],
                                 window);
    }
  }

};

class HcpPointBuilder {
private:
  KillfileReader kfr;

public:
  bool b_hcp_delta;

private:
  double hcp_delta;
  double hcp_delta_factor; // scaling that depends on the packing

  bool b_hcp_x0;
  double hcp_x0[3];

  bool b_hcp_dx0;
  double hcp_dx0[3];

  bool b_hcp_x0pdx0; // just used to indicate that it has a value
  double hcp_x0pdx0[3];

  int hcp_packing;
  double hcp_packing_factor;
  double hcp_offset_factor;

  bool b_hcp_e0;
  double hcp_e0[3];
  double hcp_e1[3],hcp_e2[3];

  int nlayers_default;

  int * send_count;
  int * send_disp;
  int * recv_count;
  int * recv_disp;

  double dxost_factor;

  int hcp_bit_max;

public:

  vector<HcpWindow> hcpWindowVec;

  bool b_plane_mesh_points;
  double plane_mesh_points_x[3];
  double plane_mesh_points_n[3];

  bool b_all_mesh_points;

  // current set of points, not necessarily all points (may be just the plane above)...
  vector<HcpVertex> vertexVec;
  int8 np_global;

  bool b_level_max;
  int level_max;
  int8* vv_count;

  // written to json for client
  //vector<View> viewVec;

  HcpPointBuilder() {

    //COUT1("HcpPointBuilder()");

    b_hcp_delta = false;
    b_hcp_x0 = false;
    b_hcp_dx0 = false;
    b_hcp_x0pdx0 = false;

    hcp_packing        = HCP_PACKING_ROOT3;
    hcp_packing_factor = 0.5*sqrt(3.0);
    hcp_delta_factor   = 1.0;
    hcp_offset_factor  = 1.0; // hcp = 1.0, cart = 0.0

    b_hcp_e0 = false;
    hcp_e0[0] = 1.0;
    hcp_e0[1] = 0.0;
    hcp_e0[2] = 0.0;
    hcp_e1[0] = 0.0;
    hcp_e1[1] = 1.0;
    hcp_e1[2] = 0.0;
    hcp_e2[0] = 0.0;
    hcp_e2[1] = 0.0;
    hcp_e2[2] = 1.0;

    nlayers_default = 10;
    dxost_factor = 1.5;

    hcp_bit_max = 28; 

    send_count = NULL;
    send_disp = NULL;
    recv_count = NULL;
    recv_disp = NULL;

    b_level_max = false;
    level_max = 0;

    b_plane_mesh_points = false;
    b_all_mesh_points = false;

    np_global = 0;
    vv_count = NULL;
  }

  ~HcpPointBuilder() {
    //COUT1("~HcpPointBuilder()");

    DELETE(vv_count);
    DELETE(send_count);
    DELETE(recv_count);
    DELETE(send_disp);
    DELETE(recv_disp);

  }

  void clearPointBooleans() {
    b_plane_mesh_points = false; // if any mesh points are built, they need to be rebuilt
    b_all_mesh_points = false;
    b_level_max = false;
  }

  /*
    void init() {

    COUT1("HcpPointBuilder::init()");
    if (mpi_rank==0) logger->setKillFilename("killstitch");

    }

    //generate points
    {
    Param * param = getParam("WRITE_PBIN");
    if (param) processWritePoints(&(*param));
    }

    //Run Interactive Mode
    {
    Param * param = getParam("INTERACTIVE");
    if (param) runInteractive(param);
    }

    //    for (ParamIter param = paramList_begin(); param != paramList_end(); ++param){
    //      param->touch();
    //      processParam(&(*param));
    //    }

    }
  */

  void processHcpPacking(Param * param,const bool b_help) {
    if (b_help) {
      helpHcpPacking();
    }
    else {
      try {
        const string name = MiscUtils::toUpperCase(param->getString());
        if (name == "ROOT3") {
          hcp_packing        = HCP_PACKING_ROOT3;
          hcp_packing_factor = 0.5*sqrt(3.0);
          hcp_delta_factor   = 1.0;
          hcp_offset_factor  = 1.0; // hcp = 1.0, cart = 0.0
        }
        else if (name == "ROOT2") {
          hcp_packing        = HCP_PACKING_ROOT2;
          hcp_packing_factor = 0.5*sqrt(2.0);
          hcp_delta_factor   = pow(3.0/2.0,1.0/3.0); // to get volumetric number-density of ROOT3
          hcp_offset_factor  = 1.0;
        }
        else if ((name == "ONE")||(name == "ROOT1")) { 
          hcp_packing        = HCP_PACKING_ONE;
          hcp_packing_factor = 0.5;
          hcp_delta_factor   = pow(3.0,1.0/3.0);
          hcp_offset_factor  = 1.0;
        }
        else if (name == "CART") {
          hcp_packing        = HCP_PACKING_CART;
          hcp_packing_factor = 1.0;
          hcp_delta_factor   = 1.0; // pow(3.0/4.0,1.0/3.0) ?
          hcp_offset_factor  = 0.0;
        }
        else {
          WUI(WARN," > unrecognized HCP_PACKING " << name << ". Setting to default.");
          helpHcpPacking();
          hcp_packing        = HCP_PACKING_ROOT3;
          hcp_packing_factor = 0.5*sqrt(3.0);
          hcp_delta_factor   = 1.0;
          hcp_offset_factor  = 1.0;
        }
        clearPointBooleans();
        WUI(INFO,"HCP_PACKING set to " << name);
        WebUI::webUIOutput.ensureImage();
      }
      catch(int e) {
        WUI(WARN,"expecting HCP_PACKING <string>");
        helpHcpPacking();
      }
    }
  }
  void helpHcpPacking() const {
    WUI(INFO,
        "HCP_PACKING sets HCP grid geometry. Supported cases:\n" <<
        "HCP_PACKING ROOT3 # (default)\n" <<
        "HCP_PACKING ROOT2\n" <<
        "HCP_PACKING ONE\n" <<
        "HCP_PACKING CART\n" <<
        "for more detail see [$CWIKB:stitch_topology]");
  }

  void processHcpDelta(Param * param,const bool b_help) {
    if (b_help) {
      helpHcpDelta();
    }
    else {
      try {
        const double hcp_delta_new = param->getDouble();
        hcp_delta = hcp_delta_new;
        b_hcp_delta = true;
        clearPointBooleans();
        WUI(INFO,"HCP_DELTA set to " << hcp_delta);
        WebUI::webUIOutput.ensureImage();
      }
      catch(int e) {
        if (b_hcp_delta) WUI(INFO,"HCP_DELTA currently set to " << hcp_delta);
        WUI(WARN,"expecting HCP_DELTA <double>");
        helpHcpDelta();
      }
    }
  }
  void helpHcpDelta() const {
    const double volume = getFluidVolume();
    WUI(INFO,
        "HCP_DELTA sets the background length scale of the coarsest HCP grid.\n" <<
        "  examples: for the current fluid volume of " << volume << ", try:\n" <<
        "    HCP_DELTA " << pow(volume/1000.0,1.0/3.0) << " # for background mesh of approx 1K\n" <<
        "    HCP_DELTA " << pow(volume/10000.0,1.0/3.0) << " # for background mesh of approx 10K\n" <<
        "    HCP_DELTA " << pow(volume/100000.0,1.0/3.0) << " # for background mesh of approx 100K\n" <<
        "    HCP_DELTA " << pow(volume/1000000.0,1.0/3.0) << " # for background mesh of approx 1M\n" <<
        "    HCP_DELTA " << pow(volume/10000000.0,1.0/3.0) << " # for background mesh of approx 10M\n" <<
        "    HCP_DELTA " << pow(volume/100000000.0,1.0/3.0) << " # for background mesh of approx 100M\n" <<
        "  for more detail see [$CWIKB:stitch_topology]");
  }

  void processHcpX0(Param * param,const bool b_help) {
    if (b_help) {
      helpHcpX0();
    }
    else {
      try {
        hcp_x0[0] = param->getDouble(0);
        hcp_x0[1] = param->getDouble(1);
        hcp_x0[2] = param->getDouble(2);
        b_hcp_x0 = true;
        clearPointBooleans();
        WUI(INFO,"HCP_X0 set to " << hcp_x0[0] << " " << hcp_x0[1] << " " << hcp_x0[2]);
        WebUI::webUIOutput.ensureImage();
      }
      catch(int e) {
        WUI(WARN,"expecting HCP_X0 <x> <y> <z>");
        helpHcpX0();
      }
    }
  }
  void helpHcpX0() const {
    if (b_hcp_x0pdx0) {
      WUI(INFO,
          "HCP_X0 sets the origin for the ray-based seeding. Current setting:\n" <<
          "HCP_X0 " << hcp_x0[0] << " " << hcp_x0[1] << " " << hcp_x0[2] << "\n" <<
          "for more detail see [$CWIKB:stitch_topology]");
    }
    else {
      WUI(INFO,
          "HCP_X0 sets the origin for the ray-based seeding. Example:\n" <<
          "HCP_X0 " << 0.1 << " " << 1.2 << " " << 2.0 << "\n" <<
          "for more detail see [$CWIKB:stitch_topology]");
    }
  }

  void processHcpDX0(Param * param,const bool b_help) {
    if (b_help) {
      helpHcpDX0();
    }
    else {
      try {
        hcp_dx0[0] = param->getDouble(0);
        hcp_dx0[1] = param->getDouble(1);
        hcp_dx0[2] = param->getDouble(2);
        b_hcp_dx0 = true;
        clearPointBooleans();
        WUI(INFO,"HCP_DX0 set to " << hcp_dx0[0] << " " << hcp_dx0[1] << " " << hcp_dx0[2]);
        WebUI::webUIOutput.ensureImage();
      }
      catch(int e) {
        WUI(WARN,"expecting HCP_DX0 <dx> <dy> <dz>");
        helpHcpDX0();
      }
    }
  }
  void helpHcpDX0() const {
    if (b_hcp_x0pdx0) {
      WUI(INFO,
          "HCP_DX0 adds a shift to the origin for the ray-based seeding. Current setting:\n" <<
          "HCP_DX0 " << hcp_dx0[0] << " " << hcp_dx0[1] << " " << hcp_dx0[2] << "\n" <<
          "for more detail see [$CWIKB:stitch_topology");
    }
    else {
      WUI(INFO,
          "HCP_DX0 adds a shift to the origin for the ray-based seeding. Example:\n" <<
          "HCP_DX0 " << 0.00 << " " << 0.01 << " " << 0.02 << "\n" <<
          "for more detail see [$CWIKB:stitch_topology]");

    }
  }

  void processDxostFactor(Param * param,const bool b_help) {
    if (b_help) {
      helpDxostFactor();
    }
    else {
      try {
        const double tmp = param->getDouble();
        if (tmp <= 0.0)
          throw(1);

        dxost_factor = tmp;
        clearPointBooleans();
        WUI(INFO,"DXOST_FACTOR set to " << dxost_factor);
        WebUI::webUIOutput.ensureImage();
      }
      catch(int e) {
        WUI(WARN,"expecting DXOST_FACTOR <+double>");
        helpDxostFactor();
      }
    }
  }
  void helpDxostFactor() const {
    WUI(INFO,
        "DXOST_FACTOR sets the scaling factor for the ff_surface length scale for building transition layers. Example:\n" <<
        "DXOST_FACTOR " << 1.5 << " # (default) ");
  }

  void processNlayersDefault(Param * param,const bool b_help) {
    if (b_help) {
      helpNlayersDefault();
    }
    else {
      try {
        const int tmp = param->getInt();
        if (tmp < 0)
          throw(1);

        nlayers_default = tmp;
        clearPointBooleans();
        WUI(INFO,"NLAYERS_DEFAULT set to " << nlayers_default);
        WebUI::webUIOutput.ensureImage();
      }
      catch(int e) {
        WUI(WARN,"expecting NLAYERS_DEFAULT <+int>");
        helpNlayersDefault();
      }
    }
  }
  void helpNlayersDefault() const {
    WUI(INFO,
        "NLAYERS_DEFAULT sets the default number of layers for each transition layer. Example:\n" <<
        "NLAYERS_DEFAULT " << 10 << " # (default) " <<
        "for more detail see [$CWIKB:stitch_refinement]");
  }

  void processHcpBitMax(Param * param,const bool b_help) {
    if (b_help) {
      helpHcpBitMax();
    }
    else {
      try {
        const int tmp = param->getInt();
        if ((tmp <= 0)||(tmp > 30))
          throw(1);

        hcp_bit_max = tmp;
        clearPointBooleans();
        WUI(INFO,"HCP_BIT_MAX set to " << hcp_bit_max);
        WebUI::webUIOutput.ensureImage();
      }
      catch(int e) {
        WUI(WARN,"expecting HCP_BIT_MAX <+int le 30>");
        helpHcpBitMax();
      }
    }
  }
  void helpHcpBitMax() const {
    WUI(INFO,
        "HCP_BIT_MAX controls integer conversion of the surface for in/out checks. Example:\n" <<
        "HCP_BIT_MAX " << 28 << " # (default)" <<
        "for more detail see [$CWIKB:stitch_refinement]");
  }

  /*
  void processHcpWindows(Param * param,const bool b_help) {
    if (b_help) {
      WUI(INFO,
	  "HCP_WINDOWS report the current list of HCP_WINDOWs");
      return;
    }
    else {
      WUI(INFO,"There are currently " << hcpWindowVec.size() << " HCP_WINDOWS");
      for (int ii = 0; ii < hcpWindowVec.size(); ++ii) {
	WUI(INFO,hcpWindowVec[ii].
      
  */

  void getHcpWindowZoneIds(set<int>& zoneIdSet) const {
    assert(zoneIdSet.empty());
    for (int iwow = 0, size = hcpWindowVec.size(); iwow < size; ++iwow) {
      if (hcpWindowVec[iwow].kind == FAZONE_HCP_WINDOW) {
	for (int ii = 0, size2 = hcpWindowVec[iwow].zoneIdVec.size(); ii < size2; ++ii) {
	  zoneIdSet.insert(hcpWindowVec[iwow].zoneIdVec[ii]);
	}
      }
    }
  }
  
  void processHcpWindow(Param * param,const bool b_help) {
    if (b_help) {
      helpHcpWindow();
    }
    else {
      try {
        HcpWindow hcpWindow;
        if (parseHcpWindow(hcpWindow,param)) {
          hcpWindow.param_str = param->str();
          if (hcpWindow.iwow == hcpWindowVec.size()) {
            // add a new window
            hcpWindowVec.push_back(hcpWindow);
            COUT2(" > adding new refinement window");
          }
          else {
            assert(hcpWindow.iwow < int(hcpWindowVec.size()));
            if (hcpWindow.iwow >= 0) {
              // edited an existing window
              hcpWindowVec[hcpWindow.iwow] = hcpWindow;
              COUT2(" > editing existing refinement window");
            }
            else {
              // erase an existing window
              COUT2(" > removing existing refinement window");
              const int window_index = (-hcpWindow.iwow-1);
              if (window_index >=0 && window_index < int(hcpWindowVec.size())) {
                hcpWindowVec.erase(hcpWindowVec.begin()+window_index);
                for (int ii = 0, lim = hcpWindowVec.size(); ii < lim; ++ii) hcpWindowVec[ii].iwow = ii;
              }
              else {
                WUI(WARN,"trying to delete invalid HCP_WINDOW (index: " << window_index << " of "<<hcpWindowVec.size() << ")");
              }
            }
          }
          WUI(INFO,"HCP_WINDOW set successfully");
          clearPointBooleans();
          WebUI::webUIOutput.ensureImage();
        }
        else {
          WUI(WARN,"HCP_WINDOW: Failed. Fix syntax");
        }
      }
      catch(int e) {
        WUI(WARN,"HCP_WINDOW: incorrect or missing token");
        helpHcpWindow();
      }
    }
  }
  void helpHcpWindow() const {
    WUI(INFO,
        "HCP_WINDOW sets HCP grid refinement windows. Some examples:\n" <<
        "HCP_WINDOW FAZONE nozzle-internal-trip LEVEL=8 NLAYERS=5,10\n" <<
        "HCP_WINDOW TCONE 15  0  0  1.5   20.0 0 0  0.5  LEVEL=5 NLAYERS 5\n" <<
        "HCP_WINDOW ANNULAR_TCONE  0.5 0 0 0.42 0.54 1.0 0 0 0.495 0.505 LEVEL=7 NLAYERS=7\n" <<
        "for more detail see [$CWIKB:stitch_refinement]");
  }

private:

  bool parseHcpWindow(HcpWindow &hcpWindow, Param * param) {

    // HCP_WINDOW FAZONE nozzle-external-wall,nozzle-lip N 10 LEVEL 2 NLAYERS 5,10,10

    bool b_level = false;
    double * geom_data = hcpWindow.geom_data;

    int ierr = 0;
    int iarg = 0;
    bool b_delete = false;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if ((token == "FAZONE")||(token == "ZONE")||(token == "ZONES")) {
        // ---------------------------------------------------------------
        // HCP_WINDOW FAZONE z1,z2,z3 LEVEL <int> NLAYERS <int>
        // ---------------------------------------------------------------
        hcpWindow.kind = FAZONE_HCP_WINDOW;
        string zoneNames = param->getString(iarg++);
	assert(hcpWindow.zoneIdVec.empty());
        hcpWindow.iwow = hcpWindowVec.size();
        if (mpi_rank == 0)
          cout << " > HCP_WINDOW " << token << " " << zoneNames << "..." << endl;
        vector<string> zoneNameVec;
        MiscUtils::splitCsv(zoneNameVec,zoneNames);
        for (int ii = 0; ii < zoneNameVec.size(); ++ii) {
          // allow for a wildcard...
          bool found = false;
          for (int izn = 0; izn < PartData::zoneVec.size(); ++izn) {
            if (strcmp_wildcard(PartData::zoneVec[izn].getName(),zoneNameVec[ii])) {
	      hcpWindow.zoneIdVec.push_back(izn);
	      found = true;
            }
          }
          if (!found) {
	    // maybe it was just a zone index...
	    int izn;
	    if (from_string<int>(izn,zoneNameVec[ii],std::dec)&&(izn >= 0)&&(izn < PartData::zoneVec.size())) {
	      hcpWindow.zoneIdVec.push_back(izn);
	    }
	    else {
	      WUI(WARN,"no matching zone found for " << zoneNameVec[ii]);
	      ierr = -1;
	    }
	  }
	}
        if (ierr == -1)
          return false;
        stringstream hash;
        hash << hcpWindow.kind;
        hash << zoneNames;
        hcpWindow.hash = hash.str();
        for (int ii = 0, lim = hcpWindowVec.size(); ii < lim; ++ii) {
          if (hcpWindow.hash == hcpWindowVec[ii].hash) {
            hcpWindow.iwow = hcpWindowVec[ii].iwow;
            break;
          }
        }
        if (hcpWindow.iwow == -1) {
          hcpWindow.iwow = hcpWindowVec.size();
        }
      }
      else if (token == "BOX") {
        // ---------------------------------------------------------------
        // HCP_WINDOW BOX <x0> <x1> <y0> <y1> <z0> <z1> LEVEL <int> [NLAYERS <int>]
        //
        // for the case of a box, the presence of NLAYERS introduces
        // bloated tris accordingly to introduce levels...
        // ---------------------------------------------------------------
        hcpWindow.kind = BOX_HCP_WINDOW;
        // parse x0,x1,y0,y1,z0,z1...
        geom_data[0] = param->getDouble(iarg++);
        geom_data[1] = param->getDouble(iarg++);
        geom_data[2] = param->getDouble(iarg++);
        geom_data[3] = param->getDouble(iarg++);
        geom_data[4] = param->getDouble(iarg++);
        geom_data[5] = param->getDouble(iarg++);
        if ((geom_data[1] <= geom_data[0])||(geom_data[3] <= geom_data[2])||(geom_data[5] <= geom_data[4])) {
          CERR("BOX syntax problem. Check BOX limits:\n\n" <<
               "HCP_WINDOW BOX <x0> <x1> <y0> <y1> <z0> <z1> LEVEL <int> [NLAYERS <int>]");
        }
        if (mpi_rank == 0) {
          cout << " > HCP_WINDOW BOX " <<
            geom_data[0] << " " <<  geom_data[1] << " " <<
            geom_data[2] << " " <<  geom_data[3] << " " <<
            geom_data[4] << " " <<  geom_data[5] << "..." << endl;
        }
        // search for this box within hcpWindowVec to see if we are overwriting it
        stringstream hash;
        hash << hcpWindow.kind;
        FOR_I5 hash << geom_data[i];
        hcpWindow.hash = hash.str();
        for (int ii = 0, lim = hcpWindowVec.size(); ii < lim; ++ii) {
          if (hcpWindow.hash == hcpWindowVec[ii].hash) {
            hcpWindow.iwow = hcpWindowVec[ii].iwow;
            break;
          }
        }
        if (hcpWindow.iwow == -1) {
          hcpWindow.iwow = hcpWindowVec.size();
        }
      }
      else if (token == "TCONE") {
        // ---------------------------------------------------------------
        // HCP_WINDOW TCONE <x0> <y0> <z0> <r0> <x1> <y1> <z1> <r1> [SHEARN <nx> <ny> <nz>] LEVEL <int> [NLAYERS <int>]
        // ---------------------------------------------------------------
        hcpWindow.kind = TCONE_HCP_WINDOW;
        geom_data[0] = param->getDouble(iarg++);
        geom_data[1] = param->getDouble(iarg++);
        geom_data[2] = param->getDouble(iarg++);
        geom_data[3] = param->getDouble(iarg++);
        geom_data[4] = param->getDouble(iarg++);
        geom_data[5] = param->getDouble(iarg++);
        geom_data[6] = param->getDouble(iarg++);
        geom_data[7] = param->getDouble(iarg++);
        COUT1(" > HCP_WINDOW TCONE " <<
              geom_data[0] << " " <<  geom_data[1] << " " <<
              geom_data[2] << " " <<  geom_data[3] << " " <<
              geom_data[4] << " " <<  geom_data[5] << " " <<
              geom_data[6] << " " <<  geom_data[7] );

        geom_data[8] = 0.0;
        geom_data[9] = 0.0;
        geom_data[10] = 0.0;
        if (param->getString(iarg) == "SHEARN") {
          ++iarg;
          geom_data[8] = param->getDouble(iarg++);
          geom_data[9] = param->getDouble(iarg++);
          geom_data[10] = param->getDouble(iarg++);
          COUT1("    > endcaps fixed with normal: " << geom_data[8] << " " <<  geom_data[9] << " " << geom_data[10]);
        }
        stringstream hash;
        hash << hcpWindow.kind;
        for (int i = 0; i < 10; ++i) hash << geom_data[i];
        hcpWindow.hash = hash.str();
        for (int ii = 0, lim = hcpWindowVec.size(); ii < lim; ++ii) {
          if (hcpWindow.hash == hcpWindowVec[ii].hash) {
            hcpWindow.iwow = hcpWindowVec[ii].iwow;
            break;
          }
        }
        if (hcpWindow.iwow == -1) {
          hcpWindow.iwow = hcpWindowVec.size();
        }
      }
      else if (token == "ANNULAR_TCONE") {
        // ---------------------------------------------------------------
        // HCP_WINDOW ANNULAR_TCONE <x0> <y0> <z0> <r0-0> <r0-1> <x1> <y1> <z1> <r1-0> <r1-1> [SHEARN <nx> <ny> <nz>] LEVEL <int> [NLAYERS <int>]
        // ---------------------------------------------------------------
        hcpWindow.kind = ANNULAR_TCONE_HCP_WINDOW;
        geom_data[0] = param->getDouble(iarg++);
        geom_data[1] = param->getDouble(iarg++);
        geom_data[2] = param->getDouble(iarg++);
        geom_data[3] = param->getDouble(iarg++);
        geom_data[4] = param->getDouble(iarg++);
        geom_data[5] = param->getDouble(iarg++);
        geom_data[6] = param->getDouble(iarg++);
        geom_data[7] = param->getDouble(iarg++);
        geom_data[8] = param->getDouble(iarg++);
        geom_data[9] = param->getDouble(iarg++);
        COUT1(" > HCP_WINDOW ANNULAR_TCONE " <<
              geom_data[0] << " " <<  geom_data[1] << " " <<
              geom_data[2] << " " <<  geom_data[3] << " " <<
              geom_data[4] << " " <<  geom_data[5] << " " <<
              geom_data[6] << " " <<  geom_data[7] << " " <<
              geom_data[8] << " " <<  geom_data[9] );
        geom_data[10] = 0.0;
        geom_data[11] = 0.0;
        geom_data[12] = 0.0;
        if (param->getString(iarg) == "SHEARN") {
          ++iarg;
          geom_data[10] = param->getDouble(iarg++);
          geom_data[11] = param->getDouble(iarg++);
          geom_data[12] = param->getDouble(iarg++);
          COUT1("    > endcaps fixed with normal: " << geom_data[10] << " " <<  geom_data[11] << " " << geom_data[12]);
        }
        stringstream hash;
        hash << hcpWindow.kind;
        for (int i = 0; i < 12; ++i) hash << geom_data[i];
        hcpWindow.hash = hash.str();
        for (int ii = 0, lim = hcpWindowVec.size(); ii < lim; ++ii) {
          if (hcpWindow.hash == hcpWindowVec[ii].hash) {
            hcpWindow.iwow = hcpWindowVec[ii].iwow;
            break;
          }
        }
        if (hcpWindow.iwow == -1) {
          hcpWindow.iwow = hcpWindowVec.size();
        }
      }
      else if (token == "ELLIPSE_TCONE") {
        // ---------------------------------------------------------------
        // HCP_WINDOW ELLIPSE_TCONE <x0> <y0> <z0> <r0-M> <r0-m> <x1> <y1> <z1> <r1-M> <r1-m> [ORTHN <nxM> <nyM> <nzM>] [SHEARN <nx> <ny> <nz>]  LEVEL <int> [NLAYERS <int>]
        // ---------------------------------------------------------------
        hcpWindow.kind = ELLIPSE_TCONE_HCP_WINDOW;
        geom_data[0] = param->getDouble(iarg++);
        geom_data[1] = param->getDouble(iarg++);
        geom_data[2] = param->getDouble(iarg++);
        geom_data[3] = param->getDouble(iarg++);
        geom_data[4] = param->getDouble(iarg++);
        geom_data[5] = param->getDouble(iarg++);
        geom_data[6] = param->getDouble(iarg++);
        geom_data[7] = param->getDouble(iarg++);
        geom_data[8] = param->getDouble(iarg++);
        geom_data[9] = param->getDouble(iarg++);
        COUT1(" > HCP_WINDOW ELLIPSE_TCONE " <<
              geom_data[0] << " " <<  geom_data[1] << " " <<
              geom_data[2] << " " <<  geom_data[3] << " " <<
              geom_data[4] << " " <<  geom_data[5] << " " <<
              geom_data[6] << " " <<  geom_data[7] << " " <<
              geom_data[8] << " " <<  geom_data[9] );
        geom_data[10] = 0.0;
        geom_data[11] = 0.0;
        geom_data[12] = 0.0;
        if (param->getString(iarg) == "ORTHN") {
          ++iarg;
          geom_data[10] = param->getDouble(iarg++);
          geom_data[11] = param->getDouble(iarg++);
          geom_data[12] = param->getDouble(iarg++);
          COUT1("    > major axis based on: " << geom_data[10] << " " <<  geom_data[11] << " " << geom_data[12]);
        }
        geom_data[13] = 0.0;
        geom_data[14] = 0.0;
        geom_data[15] = 0.0;
        if (param->getString(iarg) == "SHEARN") {
          ++iarg;
          geom_data[13] = param->getDouble(iarg++);
          geom_data[14] = param->getDouble(iarg++);
          geom_data[15] = param->getDouble(iarg++);
          COUT1("    > endcaps fixed with normal: " << geom_data[13] << " " <<  geom_data[14] << " " << geom_data[15]);
        }
        stringstream hash;
        hash << hcpWindow.kind;
        FOR_I16 hash << geom_data[i];
        hcpWindow.hash = hash.str();
        for (int ii = 0, lim = hcpWindowVec.size(); ii < lim; ++ii) {
          if (hcpWindow.hash == hcpWindowVec[ii].hash) {
            hcpWindow.iwow = hcpWindowVec[ii].iwow;
            break;
          }
        }
        if (hcpWindow.iwow == -1) {
          hcpWindow.iwow = hcpWindowVec.size();
        }
      }
      else if (token == "ANNULAR_ELLIPSE_TCONE") {
        // ---------------------------------------------------------------
        // HCP_WINDOW ANNULAR_ELLIPSE_TCONE <x0> <y0> <z0> <r0-0M> <r0-0m> <r0-1M> <r0-1m> <x1> <y1> <z1> <r1-0M> <r1-0m> <r1-1M> <r1-1m> [ORTHN <nxM> <nyM> <nzM>] [SHEARN <nx> <ny> <nz>]  LEVEL <int> [NLAYERS <int>]
        // ---------------------------------------------------------------
        hcpWindow.kind = ANNULAR_ELLIPSE_TCONE_HCP_WINDOW;
        geom_data[0] = param->getDouble(iarg++);
        geom_data[1] = param->getDouble(iarg++);
        geom_data[2] = param->getDouble(iarg++);
        geom_data[3] = param->getDouble(iarg++);
        geom_data[4] = param->getDouble(iarg++);
        geom_data[5] = param->getDouble(iarg++);
        geom_data[6] = param->getDouble(iarg++);
        geom_data[7] = param->getDouble(iarg++);
        geom_data[8] = param->getDouble(iarg++);
        geom_data[9] = param->getDouble(iarg++);
        geom_data[10] = param->getDouble(iarg++);
        geom_data[11] = param->getDouble(iarg++);
        geom_data[12] = param->getDouble(iarg++);
        geom_data[13] = param->getDouble(iarg++);
        COUT1(" > HCP_WINDOW ANNULAR ELLIPSE_TCONE " <<
              geom_data[0] << " " <<  geom_data[1] << " " <<
              geom_data[2] << " " <<  geom_data[3] << " " <<
              geom_data[4] << " " <<  geom_data[5] << " " <<
              geom_data[6] << " " <<  geom_data[7] << " " <<
              geom_data[8] << " " <<  geom_data[9] << " " <<
              geom_data[10] << " " <<  geom_data[11] << " " <<
              geom_data[12] << " " <<  geom_data[13] );
        geom_data[14] = 0.0;
        geom_data[15] = 0.0;
        geom_data[16] = 0.0;
        if (param->getString(iarg) == "ORTHN") {
          ++iarg;
          geom_data[14] = param->getDouble(iarg++);
          geom_data[15] = param->getDouble(iarg++);
          geom_data[16] = param->getDouble(iarg++);
          COUT1("    > major axis based on: " << geom_data[14] << " " <<  geom_data[15] << " " << geom_data[16]);
        }
        geom_data[17] = 0.0;
        geom_data[18] = 0.0;
        geom_data[19] = 0.0;
        if (param->getString(iarg) == "SHEARN") {
          ++iarg;
          geom_data[17] = param->getDouble(iarg++);
          geom_data[18] = param->getDouble(iarg++);
          geom_data[19] = param->getDouble(iarg++);
          COUT1("    > endcaps fixed with normal: " << geom_data[17] << " " <<  geom_data[18] << " " << geom_data[19]);
        }
        stringstream hash;
        hash << hcpWindow.kind;
        for (int i = 0; i < 20; ++i) hash << geom_data[i];
        hcpWindow.hash = hash.str();
        for (int ii = 0, lim = hcpWindowVec.size(); ii < lim; ++ii) {
          if (hcpWindow.hash == hcpWindowVec[ii].hash) {
            hcpWindow.iwow = hcpWindowVec[ii].iwow;
            break;
          }
        }
        if (hcpWindow.iwow == -1) {
          hcpWindow.iwow = hcpWindowVec.size();
        }
      }
      else if (token == "PRISM") {
        // ---------------------------------------------------------------
        // HCP_WINDOW PRISM <x0> <y0> <z0> <rM0> <rm0> <x1> <y1> <z1> <rM1> <rm1> [ORTHN <nx1> <ny1> <nz1>] [SHEARN <nx> <ny> <nz>] LEVEL <int> [NLAYERS <int>]
        // ---------------------------------------------------------------
        bool orthn_set = false;
        hcpWindow.kind = PRISM_HCP_WINDOW;
        geom_data[0] = param->getDouble(iarg++);
        geom_data[1] = param->getDouble(iarg++);
        geom_data[2] = param->getDouble(iarg++);
        geom_data[3] = param->getDouble(iarg++);
        geom_data[4] = param->getDouble(iarg++);
        geom_data[5] = param->getDouble(iarg++);
        geom_data[6] = param->getDouble(iarg++);
        geom_data[7] = param->getDouble(iarg++);
        geom_data[8] = param->getDouble(iarg++);
        geom_data[9] = param->getDouble(iarg++);
        geom_data[10] = 0.0;
        geom_data[11] = 0.0;
        geom_data[12] = 0.0;
        if (param->getString(iarg) == "ORTHN") {
          ++iarg;
          orthn_set = true;
          geom_data[10] = param->getDouble(iarg++);
          geom_data[11] = param->getDouble(iarg++);
          geom_data[12] = param->getDouble(iarg++);
        }
        COUT1(" > HCP_WINDOW PRISM");
        COUT1("    > x0,l_major,l_minor: " << geom_data[0] << " " <<  geom_data[1] << " " << geom_data[2] << ", " <<  geom_data[3] << ", " << geom_data[4]);
        COUT1("    > x1,l_major,l_minor: " << geom_data[5] << " " <<  geom_data[6] << " " << geom_data[7] << ", " <<  geom_data[8] << ", " << geom_data[9]);
        if (orthn_set) {
          COUT1("    > major axis based on: " << geom_data[10] << " " <<  geom_data[11] << " " << geom_data[12]);
        }
        geom_data[13] = 0.0;
        geom_data[14] = 0.0;
        geom_data[15] = 0.0;
        if (param->getString(iarg) == "SHEARN") {
          ++iarg;
          geom_data[13] = param->getDouble(iarg++);
          geom_data[14] = param->getDouble(iarg++);
          geom_data[15] = param->getDouble(iarg++);
          COUT1("    > endcaps fixed with normal: " << geom_data[13] << " " <<  geom_data[14] << " " << geom_data[15]);
        }
        stringstream hash;
        hash << hcpWindow.kind;
        FOR_I16 hash << geom_data[i];
        hcpWindow.hash = hash.str();
        for (int ii = 0, lim = hcpWindowVec.size(); ii < lim; ++ii) {
          if (hcpWindow.hash == hcpWindowVec[ii].hash) {
            hcpWindow.iwow = hcpWindowVec[ii].iwow;
            break;
          }
        }
        if (hcpWindow.iwow == -1) {
          hcpWindow.iwow = hcpWindowVec.size();
        }
      }
      else if (token == "HEMISPHERE") {
        // HEMISPHERE xc yc zc nx ny nz r...
        hcpWindow.kind = HEMISPHERE_HCP_WINDOW;
        geom_data[0] = param->getDouble(iarg++);
        geom_data[1] = param->getDouble(iarg++);
        geom_data[2] = param->getDouble(iarg++);
        geom_data[3] = param->getDouble(iarg++);
        geom_data[4] = param->getDouble(iarg++);
        geom_data[5] = param->getDouble(iarg++);
        geom_data[6] = param->getDouble(iarg++);
        COUT1(" > HCP_WINDOW HEMISPHERE");
        COUT1("    > xc: " << geom_data[0] << " " << geom_data[1] << " " << geom_data[2] <<
              ", normal: " << geom_data[3] << " " << geom_data[4] << " " << geom_data[5] <<
              ", r: " << geom_data[6]);
        stringstream hash;
        hash << hcpWindow.kind;
        FOR_I7 hash << geom_data[i];
        hcpWindow.hash = hash.str();
        for (int ii = 0, lim = hcpWindowVec.size(); ii < lim; ++ii) {
          if (hcpWindow.hash == hcpWindowVec[ii].hash) {
            hcpWindow.iwow = hcpWindowVec[ii].iwow;
            break;
          }
        }
        if (hcpWindow.iwow == -1) {
          hcpWindow.iwow = hcpWindowVec.size();
        }
      }
      else if (token == "SBIN") {
        // ---------------------------------------------------------------
        // HCP_WINDOW SBIN <filename> LEVEL <int> [NLAYERS <int>]
        // ---------------------------------------------------------------
        hcpWindow.kind = SBIN_HCP_WINDOW;
        hcpWindow.sbin_filename = param->getString(iarg++);
        COUT1(" > HCP_WINDOW SBIN " << hcpWindow.sbin_filename);
        stringstream hash;
        hash << hcpWindow.kind;
        hash << hcpWindow.sbin_filename;
        hcpWindow.hash = hash.str();
        for (int ii = 0, lim = hcpWindowVec.size(); ii < lim; ++ii) {
          if (hcpWindow.hash == hcpWindowVec[ii].hash) {
            hcpWindow.iwow = hcpWindowVec[ii].iwow;
            break;
          }
        }
        if (hcpWindow.iwow == -1) {
          hcpWindow.iwow = hcpWindowVec.size();
        }
      }
      else if (token == "LEVEL") {
        b_level = true;
        hcpWindow.level = param->getInt(iarg++);
        if (hcpWindow.level == 0)
          b_delete = true;
      }
      else if ((token == "DX")||(token == "DELTA")) {
        // DX can be specified instead of LEVEL, and the level will be computed
        // for you as the largest level with length scale at or below dx...
        hcpWindow.b_dx = true;
        hcpWindow.dx = param->getDouble(iarg++);
      }
      else if (token == "D") {
        // specify a distance or thickness: used for the first layer thickness of FAZONE refinement...
        hcpWindow.b_d = true;
        hcpWindow.d = param->getDouble(iarg++);
      }
      else if (token == "N") {
        // specify the count of layers: used for the first layer count in FAZONE refinement...
        hcpWindow.b_n = true;
        hcpWindow.n = param->getInt(iarg++);
      }
      else if (token == "NLAYERS") {
        // NLAYERS accepts a comma-delimited string of int: e.g. 5,10,20
        hcpWindow.b_nlayers = true;
        hcpWindow.nlayersString = param->getString(iarg++);
      }
      else if (token == "DLAYERS") {
        // DLAYERS accepts a comma-delimited string of doubles: e.g. 0.1,0.1,0.4
        hcpWindow.b_dlayers = true;
        hcpWindow.dlayersString = param->getString(iarg++);
      }
      else {
        if (mpi_rank == 0) cout << "unrecognized HCP_WINDOW token: " << token << ". Skipping..." << endl;
      }
    }

    if (hcpWindow.kind == UNKNOWN_HCP_WINDOW) {
      WUI(WARN,"HCP_WINDOW: window kind not specified.\nHELP HCP_WINDOW:");
      helpHcpWindow();
      return false;
    }

    if ((!b_level)&&(!hcpWindow.b_dx)) {
      WUI(WARN,"HCP_WINDOW: missing either LEVEL <int> or DX <double>");
      return false;
    }

    if (b_level&&((hcpWindow.level < 0)||(hcpWindow.level > LEVEL_MASK_BITS))) {
      WUI(WARN,"HCP_WINDOW: LEVEL out of range: 0 <= level <= " << LEVEL_MASK_BITS);
      return false;
    }

    if (b_delete)
      hcpWindow.iwow = -hcpWindow.iwow-1; // -1 index to indicate deletion

    return true; //successfully parsed

  }

  int getLevelMax() {
    if (!b_level_max) calcLevelMax();
    return level_max;
  }

  void calcLevelMax() {
    assert(!b_level_max);

    // -----------------------------------------------------------------------
    // HCP_WINDOW levels...
    // -----------------------------------------------------------------------
    level_max = 0;
    for (vector<HcpWindow>::iterator iter = hcpWindowVec.begin(); iter != hcpWindowVec.end(); ++iter) {
      if (iter->b_dx) {
        assert(b_hcp_delta);
        iter->level = 0;
        while (iter->dx < 0.667*hcp_delta/double(1<<iter->level))
          ++iter->level;
      }
      level_max = max(level_max,iter->level);
    }

    // -----------------------------------------------------------------------
    // part boundary levels...
    // the ff_surfaces can introduce refinement windows. To get the max level
    // associated with these windows, we just need the min dx for all ff_surface tris...
    // -----------------------------------------------------------------------
    //
    // get ff_level_max each part...
    for (int ipart = 0, lim = PartData::partVec.size(); ipart < lim; ++ipart) {
      if (PartData::partVec[ipart]->ff_surface) {
        assert(PartData::partVec[ipart]->ff_surface_dxost);
        double my_dxost_min = HUGE_VAL;
        for (int ist = mpi_rank; ist < PartData::partVec[ipart]->ff_surface->nst; ist += mpi_size) {
          my_dxost_min = min(my_dxost_min,PartData::partVec[ipart]->ff_surface_dxost[ist]);
        }
        double dxost_min;
        MPI_Allreduce(&my_dxost_min,&dxost_min,1,MPI_DOUBLE,MPI_MIN,mpi_comm);
        PartData::partVec[ipart]->ff_level_max = 0;
        if (dxost_min != HUGE_VAL) {
          while (dxost_factor*dxost_min < hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<PartData::partVec[ipart]->ff_level_max))
            ++PartData::partVec[ipart]->ff_level_max;
        }
        level_max = max(level_max,PartData::partVec[ipart]->ff_level_max);
      }
    }

    b_level_max = true;

  }

public:

  void processCountPoints(Param *param,const bool b_help) {

    // get hcp points...
    ensurePoints();

    // start with the hcp points...
    int8 my_buf[3];
    my_buf[0] = vertexVec.size();

    // need to add points added from parts...
    my_buf[1] = 0;
    my_buf[2] = 0;
    PartData::prepareIsInside();
    for (int ipart = 0; ipart < PartData::partVec.size(); ++ipart) {
      if (PartData::partVec[ipart]->pts) {
        for (int ip = 0; ip < PartData::partVec[ipart]->pts->np; ++ip) {
          // we are going to include this point in the build, if and only if it
          // is NOT inside the FF of every part later than us in the list (these parts take priority)
          // AND it is not inside the solid of any part. Since we have checked the FF of the
          // later parts, we can just check the solid of earlier parts...
          bool in = true;
          for (int ipart2 = ipart+1; ipart2 < PartData::partVec.size(); ++ipart2) {
            if (PartData::partVec[ipart2]->hasFF()) {
              if (PartData::partVec[ipart2]->isInsideFF(PartData::partVec[ipart]->pts->xp[ip])) {
                in = false;
                break;
              }
            }
            else {
              if (PartData::partVec[ipart2]->isInsideSolid(PartData::partVec[ipart]->pts->xp[ip],false)) {
                in = false;
                break;
              }
            }
          }
          if (in) {
            for (int ipart2 = 0; ipart2 < ipart; ++ipart2) {
              if (PartData::partVec[ipart2]->isInsideSolid(PartData::partVec[ipart]->pts->xp[ip],ipart2==0)) { // note flag true when ipart2 is the first part
                in = false;
                break;
              }
            }
          }
          if (in) ++my_buf[1];
          else ++my_buf[2];
        }
      }
    }

    int8 buf[3];
    MPI_Reduce(my_buf,buf,3,MPI_INT8,MPI_SUM,0,mpi_comm);

    if (mpi_rank == 0) {

      const int8 count = buf[0]+buf[1];
      cout << " > global point count: " << count << " (hcp: " << buf[0] << ", part: " << buf[1] << ", blanked part: " << buf[2] << ")" << endl;

      UIMessage_ msg(INFO); // TODO: use other format eventually ?
      msg.setClosable(true);
      stringstream tmp,ss;
      tmp << "total point count: " << count;
      msg.addLine(stringstream().flush() << tmp.str()); ss << tmp.str() << endl; tmp.str(string());
      tmp << "hcp point count: " << np_global << " (" << double(np_global)/double(count)*100.0 << "%)";
      msg.addLine(stringstream().flush() << tmp.str()); ss << tmp.str() << endl; tmp.str(string());
      if (buf[0]) {
        for (int level = 0; level <= level_max; ++level) {
          tmp << " > level: " << level << " count: " << vv_count[level] << " (" << double(vv_count[level])/double(np_global)*100.0 << "%)";
          msg.addLine(stringstream().flush() << tmp.str()); ss << tmp.str() << endl; tmp.str(string());
        }
      }
      WebUI::webUIOutput.add_(msg);
      cout << ss.str() << endl; cout.flush();

    }

    WebUI::webUIOutput.message_flag = true; // HACK so all ranks know about message

  }

  void processReportLevels(Param *param,const bool b_help) {

    if (b_hcp_delta) {

      // make sure level max is up to date...
      int level_max = getLevelMax();

      if (mpi_rank == 0) {
        stringstream ss;
        UIMessage_ msg(INFO);
        msg.setClosable(true);
        stringstream tmp;
        tmp << "hcp_delta: " << hcp_delta;
        msg.addLine(stringstream().flush() << tmp.str()); ss << tmp.str() << endl; tmp.str(string());
        for (int level = 0; level <= level_max; ++level) {
          tmp << " > level: " << level << " resolution: " << hcp_delta/pow(2.0,level);
          msg.addLine(stringstream().flush() << tmp.str()); ss << tmp.str() << endl; tmp.str(string());
        }
        WebUI::webUIOutput.add_(msg);
        cout << ss.str() << endl; cout.flush();
      }
    }
    else {

      if (mpi_rank == 0) {
        UIMessage_ msg(INFO);
        msg.setClosable(true);
        msg.addLine("Please set HCP_DELTA before requesting LEVEL report!");
        WebUI::webUIOutput.add_(msg);
        cout << "Please set HCP_DELTA before requesting LEVEL report!" << endl; cout.flush();
      }

    }

    WebUI::webUIOutput.message_flag = true; // HACK so all ranks know about message

  }

  void processWritePts(Param * param) {

    string filename = "hcp_points.pbin";
    int iarg = 0;
    while (iarg < param->size() && iarg < 1) {
      string token = param->getString(iarg++);
      // assume this is a filename...
      filename = token;
    }

    // if the points are not current or are not all there,
    // (re)generate them
    ensurePoints();
    writeBinary(filename);

  }

  void processReadPts(Param * param) {

    string filename = "hcp_points.pbin";
    int iarg = 0;
    while (iarg < param->size() && iarg < 1) {
      string token = param->getString(iarg++);
      // assume this is a filename...
      filename = token;
    }

    readBinary(filename);

  }
  
public:

  void writeHcpPointsToPartData() {

    // when we are done interacting with the points, we write them to the
    // hcpPts in namespace PartData. Note that these points are not necessarily
    // associated with any given part.

    ensurePoints();

    assert(PartData::hcpPts == NULL);
    PartData::hcpPts = new Points();

    // counts...

    PartData::hcpPts->np = vertexVec.size();
    int8 np_int8 = PartData::hcpPts->np;
    MPI_Allreduce(&np_int8,&PartData::hcpPts->np_global,1,MPI_INT8,MPI_SUM,mpi_comm);

    if (mpi_rank == 0) cout << "writeHcpPointsToPartData: hcp point count: " << PartData::hcpPts->np_global << endl;
    assert(PartData::hcpPts->np_global == np_global);

    // xp...

    assert(PartData::hcpPts->xp == NULL); PartData::hcpPts->xp = new double[PartData::hcpPts->np][3];
    for (int ip = 0; ip < PartData::hcpPts->np; ++ip)
      FOR_I3 PartData::hcpPts->xp[ip][i] = vertexVec[ip].x[i];

    // delta...

    assert(PartData::hcpPts->delta == NULL); PartData::hcpPts->delta = new double[PartData::hcpPts->np];
    for (int ip = 0; ip < PartData::hcpPts->np; ++ip)
      PartData::hcpPts->delta[ip] = getPointsDeltaForLevel(vertexVec[ip].level);

    // PartData::hcpPts are NOT from voronoi diagram. they are from hcp point builder...
    PartData::b_vp = false;
  }

private:

  //Zone based bounding box computation...
  /*
    void computeZoneView(View &zoneView, const string &zoneNames) {

    double zone_bbminmax[6] = {HUGE_VAL,HUGE_VAL,HUGE_VAL,HUGE_VAL,HUGE_VAL,HUGE_VAL};
    for (int ipart = 0; ipart < PartData::partVec.size(); ++ipart) {
    if (PartData::partVec[ipart]->surface) {
    double part_bbminmax[6];
    PartData::partVec[ipart]->surface->calcZoneBbminmax(part_bbminmax, zoneNames);;
    FOR_I6 zone_bbminmax[i] = min(zone_bbminmax[i],part_bbminmax[i]);
    }
    }

    // recall zone_bbminmax[3,4,5] contains the max as negative...
    double center[3] = {0.5*(zone_bbminmax[0]-zone_bbminmax[3]),
    0.5*(zone_bbminmax[1]-zone_bbminmax[4]),
    0.5*(zone_bbminmax[2]-zone_bbminmax[5])};

    double dimensions[3] = {-(zone_bbminmax[3]+zone_bbminmax[0]),
    -(zone_bbminmax[4]+zone_bbminmax[1]),
    -(zone_bbminmax[5]+zone_bbminmax[2])};

    //For now just create a single view...
    FOR_I3 zoneView.plane_xp[i] = center[i];
    FOR_I3 zoneView.view_xp[i] = center[i];
    zoneView.plane_np[0] = 0.0;
    zoneView.plane_np[1] = -1.0;
    zoneView.plane_np[2] = 0.0;
    zoneView.view_width = 1.05*sqrt(dimensions[0]*dimensions[0]+
    dimensions[1]*dimensions[1]+
    dimensions[2]*dimensions[2]);
    }
  */

  double getFluidVolume() const {
    // use the surfaces to efficiently manage the volume calculation...
    double volume = 0.0;
    for (int ipart = 0; ipart < PartData::partVec.size(); ++ipart) {
      if (PartData::partVec[ipart]->surface) {
        // note the sign flip here: surfaces return positive volume in
        // solids...
        volume -= PartData::partVec[ipart]->surface->getVolume();
      }
    }
    return volume;
  }

  void setJKFrom2D(int jk[2],const double xp[2],const double hcp_dxp) const {
    // this routine converts a 2d double already in the (e1,e2) plane into (j,k).
    // We make it odd because rays are even and geometry is odd, ensuring we
    // never have ray-vertex collisions (although ray-edge is possible)...
    FOR_I2 {
      jk[i] = (int)floor(hcp_dxp*xp[i]);
      if (jk[i]%2 == 0) jk[i] += 1; // ensure odd...
    }
  }

  void setJKFrom3D(int jk[2],const double x[3],const double hcp_dxp) const {
    // see note above about odd integer space...
    // j is along hcp_e1...
    jk[0] = (int)floor(hcp_dxp*((x[0]-hcp_x0pdx0[0])*hcp_e1[0] +
                                (x[1]-hcp_x0pdx0[1])*hcp_e1[1] +
                                (x[2]-hcp_x0pdx0[2])*hcp_e1[2]));
    if (jk[0]%2 == 0) jk[0] += 1; // ensure odd...
    // and k is along hcp_e2...
    jk[1] = (int)floor(hcp_dxp*((x[0]-hcp_x0pdx0[0])*hcp_e2[0] +
                                (x[1]-hcp_x0pdx0[1])*hcp_e2[1] +
                                (x[2]-hcp_x0pdx0[2])*hcp_e2[2]));
    if (jk[1]%2 == 0) jk[1] += 1; // ensure odd...
  }

  // note these are slightly different than 2*hcp_packing_factor
  double getPointsDeltaForLevel(const int level) const {

    // this returns the initial guess for the points delta -- i.e. the
    // sphere inside which all nbrs will be collected and cut against to
    // build the VD...

    //double factor = getDoubleParam("FACTOR");

    //assert((level >= 0)&&(level <= HCP_LEVEL_MAX));
    //assert(got_hcp_delta);
    //assert(got_hcp_packing);

    // note: these strange packing factors are designed to set the initial
    // sphere for nbr seach large enough to complete the voronoi diagram
    // for a particular packing in one iteration with minimal nbrs...

    switch (hcp_packing) {
    case HCP_PACKING_ONE:
      return hcp_delta_factor*hcp_delta/double(1<<level)*1.215;
    case HCP_PACKING_ROOT2:
      return hcp_delta_factor*hcp_delta/double(1<<level)*1.405;
    case HCP_PACKING_ROOT3:
      return hcp_delta_factor*hcp_delta/double(1<<level)*1.575;
    case HCP_PACKING_CART:
      return hcp_delta_factor*hcp_delta/double(1<<level)*1.000; // TODO determine good val
    default:
      assert(0);
    }

    // should never get here...
    return 0.0;

  }

public:

  // note these are slightly different than 2*hcp_packing_factor
  double getPointsRvvForLevel(const int level) const {

    assert(b_hcp_delta);
    // used when drawing on a plane...
    return hcp_packing_factor*hcp_delta_factor*hcp_delta/double(1<<level);
  }

  void ensurePoints(const double plane_xp[3] = NULL,const double plane_np[3] = NULL) {

    if (PartData::partVec.empty())
      return;

    if (PartData::partVec[0]->pts) // if the first part has pts, no hcp...
      return;

    if (b_all_mesh_points) // if we have all mesh points, just return...
      return;

    if (plane_xp == NULL) { // NULL means build ALL points
      assert(plane_np == NULL);
      b_all_mesh_points = true;
      b_plane_mesh_points = false;
      if (mpi_rank == 0)
        cout << "ensurePointsNearPlane: building ALL mesh points" << endl;
      // hcp points are going to be updated, so clear them here...
      if (PartData::hcpPts) {
        delete PartData::hcpPts;
        PartData::hcpPts = NULL;
      }
    }
    else {
      // if we have been asked to build a plane, then check if it matches
      // the current plane...
      if (b_plane_mesh_points&&
          (plane_mesh_points_x[0] == plane_xp[0])&&
          (plane_mesh_points_x[1] == plane_xp[1])&&
          (plane_mesh_points_x[2] == plane_xp[2])&&
          (plane_mesh_points_n[0] == plane_np[0])&&
          (plane_mesh_points_n[1] == plane_np[1])&&
          (plane_mesh_points_n[2] == plane_np[2]))
        return;
      // store the plane...
      FOR_I3 plane_mesh_points_x[i] = plane_xp[i];
      FOR_I3 plane_mesh_points_n[i] = plane_np[i];
      b_plane_mesh_points = true;
      if (mpi_rank == 0)
        cout << "ensurePointsNearPlane: " << COUT_VEC(plane_xp) << " " << COUT_VEC(plane_np) << endl;
    }

    // insist the the normal is a unit normal. It is used for signed distance calcs...

    double unit_np[3];
    if (plane_xp != NULL) {
      const double np_mag = MAG(plane_np); assert(np_mag > 0.0);
      FOR_I3 unit_np[i] = plane_np[i]/np_mag;
    }

    // clear any existing point memory...
    // x_vv = NULL;
    // delta_vv = NULL;
    // nvv = 0;

    const double wtime0 = MPI_Wtime();
    double wtime_prev = wtime0;

    // if we don't have hcp_x0, set it...
    if (!b_hcp_x0) {
      double bbmin[3],bbmax[3];
      PartData::getBbox(bbmin,bbmax);
      FOR_I3 hcp_x0[i] = 0.5*(bbmin[i]+bbmax[i]);
      COUT1(" > HCP_X0 set to default: " << hcp_x0[0] << " " << hcp_x0[1] << " " << hcp_x0[2]);
    }

    if (!b_hcp_dx0) {
      FOR_I3 hcp_dx0[i] = 1.0E-4*hcp_delta;
      COUT1(" > HCP_DX0 set to default: " << hcp_dx0[0] << " " << hcp_dx0[1] << " " << hcp_dx0[2]);
    }

    // compute hcp_x0-plus-dx0...
    FOR_I3 hcp_x0pdx0[i] = hcp_x0[i] + hcp_dx0[i];
    b_hcp_x0pdx0 = true;

    // put hcp_delta parsing here to see what has to be changed when this
    // is changed by the user...
    // TODO: this should fail to produce the points, but produce the "HELP HCP_DELTA" result
    // back to the client...

    if (!b_hcp_delta) {
      WUI(WARN,"missing HCP_DELTA <double>");
      //hcpDeltaHelp();
      return;
    }

    // modify hcp_delta based on packing factor...

    const double one_o_hcp_packing_factor = 1.0/hcp_packing_factor;

    // for now, write this for one surface of the first part...

    bool no_pts = false;
    assert(PartData::partVec.size() >= 1);
    for (int ipart = 0; ipart < PartData::partVec.size(); ++ipart) {
      //every ff_surface needs a length scale to build potential refinement windows
      if (PartData::partVec[ipart]->ff_surface) {
        assert(PartData::partVec[ipart]->ff_surface_dxost);
      }
      if (PartData::partVec[ipart]->pts == NULL)
        no_pts = true;
    }
    if (!no_pts) {
      if (mpi_rank == 0) cout << " > all parts have pts: skip hcp point build" << endl;
      return;
    }

    // set the surface_st_flag to the level associated with any FAZONE refinement windows...

    // need to get the current max level for allocation (based on window and parts)...
    int level_max = getLevelMax();
    assert(level_max < LEVEL_MASK_BITS);

    // get ff_nlayers for each part...
    for (int ipart = 0, lim = PartData::partVec.size(); ipart < lim; ++ipart) {
      if (PartData::partVec[ipart]->ff_surface) {
        assert(PartData::partVec[ipart]->ff_level_max >= 0); // must be set by call getLevelMax
        vector<int> ff_nlayersVec;
        MiscUtils::splitCsv(ff_nlayersVec,PartData::partVec[ipart]->ff_nlayersString);
        if (ff_nlayersVec.empty()) {
          ff_nlayersVec.resize(PartData::partVec[ipart]->ff_level_max+1);
          for (int il = 0; il <= PartData::partVec[ipart]->ff_level_max; ++il)
            ff_nlayersVec[il] = PartData::partVec[ipart]->ff_nlayers_default;
        }
        else if (ff_nlayersVec.size() < PartData::partVec[ipart]->ff_level_max+1) {
          const int size0 = ff_nlayersVec.size();
          for (int ii = 0; ii < size0; ++ii)
            assert(ff_nlayersVec[ii] > 0);
          ff_nlayersVec.resize(PartData::partVec[ipart]->ff_level_max+1);
          for (int ii = size0; ii <= PartData::partVec[ipart]->ff_level_max; ++ii)
            ff_nlayersVec[ii] = ff_nlayersVec[size0-1];
        }
        if (mpi_rank == 0) {
          cout << "PART " << PartData::partVec[ipart]->name << " PartData::partVec[ipart]->ff_level_max_MAX " << PartData::partVec[ipart]->ff_level_max<< " ff_nlayers ";
          for (int il = 0; il < PartData::partVec[ipart]->ff_level_max; ++il)
            cout << ff_nlayersVec[il] << ",";
          cout << ff_nlayersVec[PartData::partVec[ipart]->ff_level_max] << endl;
        }
        for (int il = 0; il <= PartData::partVec[ipart]->ff_level_max; ++il)
          PartData::partVec[ipart]->ff_nlayers[il] = ff_nlayersVec[PartData::partVec[ipart]->ff_level_max-il]; // flip because users think from finest to coarsest
      }
    }

    // this routine needs periodicRtVec to do quick transforms...
    vector<PeriodicRt> periodicRtVec;
    PeriodicData::buildPeriodicRtVec(periodicRtVec);

    if (mpi_rank == 0) cout << " > applying " << hcpWindowVec.size() << " HCP_WINDOW(s)" << endl;

    int nwi = hcpWindowVec.size();

    assert(nwi < MAX_HCP_WINDOWS);
    int * wiozn = NULL;
    int * leowi = NULL;
    int * nlowi = NULL;
    int nlayers[LEVEL_MASK_BITS]; // TODO: is this the right size?

    if (nwi>0) {
      leowi = new int[nwi];
      nlowi = new int[nwi*(level_max+1)];
    }

    nwi = 0;
    vector<TriStuff::Tri> triVec;
    vector<int> tvogw_i; // trivec index of geometric window (NOT all windows)
    for (vector<HcpWindow>::iterator iter = hcpWindowVec.begin(); iter != hcpWindowVec.end(); ++iter) {
      // ----------------------------------------------------------------------------------------------
      // when we first access a delta-specified window, we need to convert to a level based on the
      // current HCP_DELTA. This should have been done during getLevelMax, but we do it again here
      // to set the levels...
      // ----------------------------------------------------------------------------------------------
      if (iter->b_dx) {
        assert(b_hcp_delta);
        iter->level = 0;
        while (iter->dx < 0.667*hcp_delta/double(1<<iter->level))
          ++iter->level;
      }
      // ----------------------------------------------------------------------------------------------
      // now that we have iter->level, compute the layers from the string...
      // ----------------------------------------------------------------------------------------------
      if (iter->b_nlayers) {
        // NLAYERS wins, if present...
        vector<int> nlayersVec;
        MiscUtils::splitCsv(nlayersVec,iter->nlayersString);
        int il = iter->level;
        nlayers[il] = nlayersVec[0];
        --il;
        for (int ii = 0; ii < nlayersVec.size(); ++ii) {
          if (il == 0)
            break;
          nlayers[il] = nlayersVec[ii];
          --il;
        }
        // and use the last nlayer entry for everything else...
        const int nlayers_last = nlayersVec.back();
        while (il > 0) {
          nlayers[il] = nlayers_last;
          --il;
        }
      }
      else if (iter->b_dlayers) {
        // then DLAYERS...
        // dlayers are specific distances that the user has specified to achieve
        // in the transition layers...
        vector<double> dlayersVec;
        MiscUtils::splitCsv(dlayersVec,iter->dlayersString);
        // put something in nlayers[il] for consistency with the discussion above...
        int il = iter->level;
        nlayers[il] = max(1,int(dlayersVec[0]*double(1<<il)/hcp_delta+0.5));
        --il;
        for (int ii = 0; ii < dlayersVec.size(); ++ii) {
          if (il == 0)
            break;
          nlayers[il] = max(1,int(dlayersVec[ii]*double(1<<il)/hcp_delta+0.5));
          --il;
        }
        // and use the last dlayer entry for everything else...
        const double dlayers_last = dlayersVec.back();
        while (il > 0) {
          nlayers[il] = max(1,int(dlayers_last*double(1<<il)/hcp_delta+0.5));
          --il;
        }
      }
      else {
        // just use the default...
        for (int il = 1; il <= iter->level; ++il)
          nlayers[il] = nlayers_default;
      }
      // nlayers[0] is not used...
      nlayers[0] = 0;
      /*
        cout << "HCP_WINDOW: " << iter->getKindStr() << endl;
        for (int il = iter->level; il >= 0; --il) {
        cout << " > level: " << il << " nlayers[il]: " << nlayers[il] << endl;
        }
        getchar();
      */
      if (iter->kind == FAZONE_HCP_WINDOW) {
        // for the case of FAZONE windowing, the user can (should) specify N, or D to be used
        // for the layer thickness. For backward compatibility, however, take the first value
        // from the NLAYERS parameter if this is missing...
        if (iter->b_n) {
          nlayers[iter->level] = iter->n;
        }
        else if (iter->b_d) {
          nlayers[iter->level] = max(1,int(iter->d*double(1<<iter->level)/hcp_delta+0.5));
        }
        else {
          if (mpi_rank == 0) cout << " > Warning: HCP_WINDOW FAZONE missing N or D. Using N=" << nlayers[iter->level] << " (first transition layer thickness)" << endl;
        }
        if (mpi_rank_shared == 0) {
          if (wiozn == NULL) {
            wiozn = new int[PartData::zoneVec.size()];
            for (int izn = 0; izn < PartData::zoneVec.size(); ++izn) {
              wiozn[izn] = -1;
            }
          }
          // the zone id's have already been parsed...
          for (int ii = 0; ii < iter->zoneIdVec.size(); ++ii) {
            // allow for a wildcard...
            const int izn = iter->zoneIdVec[ii];
	    assert((izn >= 0)&&(izn < PartData::zoneVec.size()));
	    if (wiozn[izn] != -1) {
	      if (mpi_rank == 0) cout << " > Warning: HCP_WINDOW for ZONE " << PartData::zoneVec[izn].getName() << " ZONE_ID " << izn << " has been specified more than once. Using last" << endl;
	    }
	    // overwrite...
	    wiozn[izn] = nwi;
          }
        }
        MPI_Barrier(mpi_comm_shared);
      }
      else if (iter->kind == BOX_HCP_WINDOW) {
        tvogw_i.push_back((int)triVec.size());
        TriStuff::initBox(triVec,
                          iter->geom_data[0],iter->geom_data[1],
                          iter->geom_data[2],iter->geom_data[3],
                          iter->geom_data[4],iter->geom_data[5],nwi);
      }
      else if (iter->kind == TCONE_HCP_WINDOW) {
        tvogw_i.push_back((int)triVec.size());
        TriStuff::initTcone(triVec,
                            iter->geom_data[0],iter->geom_data[1],iter->geom_data[2],iter->geom_data[3],
                            iter->geom_data[4],iter->geom_data[5],iter->geom_data[6],iter->geom_data[7],
                            iter->geom_data[8],iter->geom_data[9],iter->geom_data[10],nwi);
      }
      else if (iter->kind == ANNULAR_TCONE_HCP_WINDOW) {
        tvogw_i.push_back((int)triVec.size());
        TriStuff::initAnnularTcone(triVec,
                                   iter->geom_data[0],iter->geom_data[1],iter->geom_data[2],iter->geom_data[3],
                                   iter->geom_data[4],iter->geom_data[5],iter->geom_data[6],iter->geom_data[7],
                                   iter->geom_data[8],iter->geom_data[9],iter->geom_data[10],iter->geom_data[11],
                                   iter->geom_data[12],nwi);
      }
      else if (iter->kind == ELLIPSE_TCONE_HCP_WINDOW) {
        tvogw_i.push_back((int)triVec.size());
        TriStuff::initEllipseTcone(triVec,
                                   iter->geom_data[0],iter->geom_data[1],iter->geom_data[2],iter->geom_data[3],
                                   iter->geom_data[4],iter->geom_data[5],iter->geom_data[6],iter->geom_data[7],
                                   iter->geom_data[8],iter->geom_data[9],iter->geom_data[10],iter->geom_data[11],
                                   iter->geom_data[12],iter->geom_data[13],iter->geom_data[14],iter->geom_data[15],
                                   nwi);
      }
      else if (iter->kind == ANNULAR_ELLIPSE_TCONE_HCP_WINDOW) {
        tvogw_i.push_back((int)triVec.size());
        TriStuff::initAnnularEllipseTcone(triVec,
                                          iter->geom_data[0],iter->geom_data[1],iter->geom_data[2],iter->geom_data[3],
                                          iter->geom_data[4],iter->geom_data[5],iter->geom_data[6],iter->geom_data[7],
                                          iter->geom_data[8],iter->geom_data[9],iter->geom_data[10],iter->geom_data[11],
                                          iter->geom_data[12],iter->geom_data[13],iter->geom_data[14],iter->geom_data[15],
                                          iter->geom_data[16],iter->geom_data[17],iter->geom_data[18],iter->geom_data[19],
                                          nwi);
      }
      else if (iter->kind == PRISM_HCP_WINDOW) {
        tvogw_i.push_back((int)triVec.size());
        TriStuff::initPrism(triVec,
                            iter->geom_data[0],iter->geom_data[1],iter->geom_data[2],iter->geom_data[3],
                            iter->geom_data[4],iter->geom_data[5],iter->geom_data[6],iter->geom_data[7],
                            iter->geom_data[8],iter->geom_data[9],iter->geom_data[10],iter->geom_data[11],
                            iter->geom_data[12],iter->geom_data[13],iter->geom_data[14],iter->geom_data[15],
                            nwi);
      }
      else if (iter->kind == HEMISPHERE_HCP_WINDOW) {
        tvogw_i.push_back((int)triVec.size());
        TriStuff::initHemisphere(triVec,
                                 iter->geom_data[0],iter->geom_data[1],iter->geom_data[2], // x
                                 iter->geom_data[3],iter->geom_data[4],iter->geom_data[5], // n
                                 iter->geom_data[6],                                       // r
                                 nwi);
      }
      else if (iter->kind == SBIN_HCP_WINDOW) {
        tvogw_i.push_back((int)triVec.size());
        TriStuff::initSbin(triVec,
                           iter->sbin_filename,nwi);
      }
      else {
        assert(0);
      }
      // fast access of level and nlayers from to-be-stored window index...
      leowi[nwi] = iter->level;
      //nlowi[nwi] = iter->nlayers;
      for (int il = 0; il <= iter->level; ++il)
        nlowi[nwi*(level_max+1)+il] = nlayers[il];
      nwi++;
    }
    assert(nwi == hcpWindowVec.size());
    const int ngw = tvogw_i.size();
    tvogw_i.push_back(triVec.size());

    // once all windows are set, apply any FAZONE refinements stored in leozn/nlozn...
    if (mpi_rank_shared == 0) {
      if (mpi_rank == 0) cout << " > setting surface_st_flag to FAZONE window index..." << endl;
      for (int ipart = 0; ipart < PartData::partVec.size(); ++ipart) {
        // every surface needs an st_flag to map to potential windows...
        if (PartData::partVec[ipart]->surface) {
          if (wiozn != NULL) {
            for (int ist = 0; ist < PartData::partVec[ipart]->surface->nst; ++ist) {
              const int isz = PartData::partVec[ipart]->surface->znost[ist];
              const int izn = PartData::partVec[ipart]->znosz[isz];
              assert((izn >= 0)&&(izn < PartData::zoneVec.size()));
              if (wiozn[izn] >= 0)
                PartData::partVec[ipart]->setWindowForSt(wiozn[izn],ist);
            }
          }
        }
      }
      DELETE(wiozn);
    }
    MPI_Barrier(mpi_comm_shared);

    nwi += PartData::partVec.size();
    if (mpi_rank == 0) {
      cout << " > level_max: " << level_max << ", nwi: " << nwi << endl;
      double wtime = MPI_Wtime();
      cout << " setup time: " << wtime-wtime_prev << " [s], total time so far: " << wtime-wtime0 << " [s]" << endl;
      wtime_prev = wtime;
    }

    int * level_count = new int[(level_max+1)*(nwi+1)];

    // --------------------------------------------------------------------
    // parallel build of the 2D bbox limits for each level (still in
    // double coords)...
    // --------------------------------------------------------------------
    double (*my_bbminmax)[4] = new double[level_max+1][4];
    for (int level = 0; level <= level_max; ++level)
      FOR_I4 my_bbminmax[level][i] = HUGE_VAL;

    double (*my_plane_bbminmax)[4] = new double[level_max+1][4];
    for (int level = 0; level <= level_max; ++level)
      FOR_I4 my_plane_bbminmax[level][i] = HUGE_VAL;

    // now go through the tris and use the shm flag to set the level bbox's.
    // we do this in parallel, so no one rank has to loop on all tris...
    for (int ipart = 0; ipart < PartData::partVec.size(); ++ipart) {
      // every surface needs an st_flag to map to potential windows...
      if (PartData::partVec[ipart]->surface) {
        SurfaceShm * surface = PartData::partVec[ipart]->surface;
        int my_nst_avg = surface->nst/mpi_size;
        if (surface->nst%mpi_size) ++my_nst_avg;
        const int ist0 = min(surface->nst,mpi_rank*my_nst_avg);
        const int ist1 = min(surface->nst,(mpi_rank+1)*my_nst_avg);
        assert(ist1-ist0 <= my_nst_avg);
        for (int ist = ist0; ist < ist1; ++ist) {
          // calculate the bbox for this tri in e1,e2 space...
          double this_bbminmax[4] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };
          double xpyp[3][2];
          FOR_I3 {
            const double * const x = surface->xsp[surface->spost[ist][i]];
            const double dx[3] = DIFF(x,hcp_x0pdx0);
            // x-primed is in the hcp_e1 direction...
            xpyp[i][0] = DOT_PRODUCT(dx,hcp_e1);
            this_bbminmax[0] = min(this_bbminmax[0],xpyp[i][0]); // xp-min
            this_bbminmax[1] = min(this_bbminmax[1],-xpyp[i][0]); // -xp-max
            // and y-primed is in the hcp_e2 direction...
            xpyp[i][1] = DOT_PRODUCT(dx,hcp_e2);
            this_bbminmax[2] = min(this_bbminmax[2],xpyp[i][1]); // yp-min
            this_bbminmax[3] = min(this_bbminmax[3],-xpyp[i][1]); // -yp-max
          }
          // the tri bbox unmodified is used to set the level 0 bbminmax...
          FOR_I4 my_bbminmax[0][i] = min(my_bbminmax[0][i],this_bbminmax[i]);
          // depending on level and nlayer, we may need to
          int st_level = 0;
          int iwindow = -1;
          if (PartData::partVec[ipart]->getWindowForSt(iwindow,ist))
            st_level = (leowi != NULL) ? leowi[iwindow]:0;
          if (!((st_level >= 0)&&(st_level <= level_max)))
            cout << "got st_level: " << st_level << " level_max: " << level_max << endl;
          assert((st_level >= 0)&&(st_level <= level_max));
          if (st_level >= 1) {
            assert(iwindow >= 0);
            double rmax = 0.0;
            for (int il = 1; il <= st_level; ++il) {
              const int st_nlayers = nlowi[iwindow*(level_max+1)+il];
              //rmax += double(st_nlayers)*hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<il);
              rmax += double(st_nlayers)*hcp_delta/double(1<<il);
            }
            for (int il = 1; il <= st_level; ++il) {
              FOR_I4 my_bbminmax[il][i] = min(my_bbminmax[il][i],this_bbminmax[i]-rmax);
              const int st_nlayers = nlowi[iwindow*(level_max+1)+il];
              //rmax -= double(st_nlayers)*hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<il);
              rmax -= double(st_nlayers)*hcp_delta/double(1<<il);
            }
          }
          // if we have a plane, constrain our rays to the thickened plane...
          if (plane_xp != NULL) {
            // at any level, add this tri's bbox to my_plane_bbminmax if it is within
            // delta for the level...
            // determine the signed distance from the plane for the corners of this tri...
            double signed_distance[3];
            FOR_I3 {
              const double * const x = surface->xsp[surface->spost[ist][i]];
              const double dx[3] = DIFF(x,plane_xp);
              signed_distance[i] = DOT_PRODUCT(dx,unit_np);
            }
            for (int il = 0; il <= level_max; ++il) {
              const double delta = hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<il);
              if (((signed_distance[0] > delta)&&(signed_distance[1] > delta)&&(signed_distance[2] > delta))||
                  ((signed_distance[0] < -delta)&&(signed_distance[1] < -delta)&&(signed_distance[2] < -delta)))
                break;
              FOR_I3 {
                const int ip1 = (i+1)%3;
                // several conditions possible along this edge...
                if (signed_distance[i] < -delta) {
                  if (signed_distance[ip1] < -delta) {
                    continue; // no action
                  }
                  else if (signed_distance[ip1] <= delta) {
                    // this edge contributes ip1 and an intersection at -delta. intersection first...
                    const double wip1 = (-delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                    double this_xpyp[2];
                    FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                    my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                    my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                    my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                    my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                    // then ip1...
                    my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[ip1][0]); // xp-min
                    my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[ip1][0]); // -xp-max
                    my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[ip1][1]); // yp-min
                    my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[ip1][1]); // -yp-max
                  }
                  else {
                    // the 2-intersection case...
                    double wip1 = (-delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                    double this_xpyp[2];
                    FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                    my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                    my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                    my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                    my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                    wip1 = (delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                    FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                    my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                    my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                    my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                    my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                  }
                }
                else if (signed_distance[i] < delta) {
                  if (signed_distance[ip1] < -delta) {
                    // this edge contributes i and an intersection at -delta. intersection first...
                    const double wip1 = (-delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                    double this_xpyp[2];
                    FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                    my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                    my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                    my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                    my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                    // then i...
                    my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[i][0]); // xp-min
                    my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[i][0]); // -xp-max
                    my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[i][1]); // yp-min
                    my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[i][1]); // -yp-max
                  }
                  else if (signed_distance[ip1] <= delta) {
                    // i and ip1, i first...
                    my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[i][0]); // xp-min
                    my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[i][0]); // -xp-max
                    my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[i][1]); // yp-min
                    my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[i][1]); // -yp-max
                    // then ip1...
                    my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[ip1][0]); // xp-min
                    my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[ip1][0]); // -xp-max
                    my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[ip1][1]); // yp-min
                    my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[ip1][1]); // -yp-max
                  }
                  else {
                    // this edge contributes i and an intersection at -delta. intersection first...
                    const double wip1 = (delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                    double this_xpyp[2];
                    FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                    my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                    my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                    my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                    my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                    // then i...
                    my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[i][0]); // xp-min
                    my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[i][0]); // -xp-max
                    my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[i][1]); // yp-min
                    my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[i][1]); // -yp-max
                  }
                }
                else {
                  if (signed_distance[ip1] < -delta) {
                    // the 2-intersection case...
                    double wip1 = (-delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                    double this_xpyp[2];
                    FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                    my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                    my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                    my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                    my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                    wip1 = (delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                    FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                    my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                    my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                    my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                    my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                  }
                  else if (signed_distance[ip1] <= delta) {
                    // this edge contributes i and an intersection at -delta. intersection first...
                    const double wip1 = (delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                    double this_xpyp[2];
                    FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                    my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                    my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                    my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                    my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                    // then ip1...
                    my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[ip1][0]); // xp-min
                    my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[ip1][0]); // -xp-max
                    my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[ip1][1]); // yp-min
                    my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[ip1][1]); // -yp-max
                  }
                }
              }
            }
          }
        }
      }
      if (PartData::partVec[ipart]->ff_surface) {
        SurfaceShm * ff_surface = PartData::partVec[ipart]->ff_surface;
        int my_nst_avg = ff_surface->nst/mpi_size;
        if (ff_surface->nst%mpi_size) ++my_nst_avg;
        const int ist0 = min(ff_surface->nst,mpi_rank*my_nst_avg);
        const int ist1 = min(ff_surface->nst,(mpi_rank+1)*my_nst_avg);
        assert(ist1-ist0 <= my_nst_avg);
        for (int ist = ist0; ist < ist1; ++ist) {
          // consider all transformed versions of this tri as well...
          for (int iper = 0; iper < periodicRtVec.size(); ++iper) {
            double xsp_t[3][3];
            FOR_I3 FOR_J3 xsp_t[i][j] = ff_surface->xsp[ff_surface->spost[ist][i]][j];
            if (periodicRtVec[iper].b_t) {
              FOR_I3 FOR_J3 xsp_t[i][j] += periodicRtVec[iper].t[j];
            }
            if (periodicRtVec[iper].b_R) {
              //assert(0);
              double tmp[3];
              FOR_I3 {
                FOR_J3 tmp[j] = xsp_t[i][j];
                MiscUtils::applyRotation(xsp_t[i],periodicRtVec[iper].R,tmp);
              }
            }
            // calculate the bbox for this tri in e1,e2 space...
            double this_bbminmax[4] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };
            double xpyp[3][2];
            FOR_I3 {
              //const double * const x = ff_surface->xsp[ff_surface->spost[ist][i]];
              //const double dx[3] = DIFF(x,hcp_x0pdx0);
              const double dx[3] = DIFF(xsp_t[i],hcp_x0pdx0);
              // x-primed is in the hcp_e1 direction...
              xpyp[i][0] = DOT_PRODUCT(dx,hcp_e1);
              this_bbminmax[0] = min(this_bbminmax[0],xpyp[i][0]); // xp-min
              this_bbminmax[1] = min(this_bbminmax[1],-xpyp[i][0]); // -xp-max
              // and y-primed is in the hcp_e2 direction...
              xpyp[i][1] = DOT_PRODUCT(dx,hcp_e2);
              this_bbminmax[2] = min(this_bbminmax[2],xpyp[i][1]); // yp-min
              this_bbminmax[3] = min(this_bbminmax[3],-xpyp[i][1]); // -yp-max
            }
            // these tris cannot extend the level 0 bbox...
            //FOR_I4 my_bbminmax[0][i] = min(my_bbminmax[0][i],this_bbminmax[i]);
            int st_level = 0;
            while (dxost_factor*PartData::partVec[ipart]->ff_surface_dxost[ist] < hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<st_level))
              ++st_level;
            assert(st_level <= level_max);
            if (st_level >= 1) {
              //const int st_nlayers = PartData::partVec[ipart]->ff_nlayers;
              double rmax = 0.0;
              for (int il = 1; il <= st_level; ++il) {
                const int st_nlayers = PartData::partVec[ipart]->ff_nlayers[il];
                rmax += double(st_nlayers)*hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<il);
              }
              for (int il = 1; il <= st_level; ++il) {
                FOR_I4 my_bbminmax[il][i] = min(my_bbminmax[il][i],this_bbminmax[i]-rmax);
                const int st_nlayers = PartData::partVec[ipart]->ff_nlayers[il];
                rmax -= double(st_nlayers)*hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<il);
              }
            }
            // if we have a plane, constrain our rays to the thickened plane...
            if (plane_xp != NULL) {
              // at any level, add this tri's bbox to my_plane_bbminmax if it is within
              // delta for the level...
              // determine the signed distance from the plane for the corners of this tri...
              double signed_distance[3];
              FOR_I3 {
                //const double * const x = ff_surface->xsp[ff_surface->spost[ist][i]];
                //const double dx[3] = DIFF(x,plane_xp);
                const double dx[3] = DIFF(xsp_t[i],plane_xp);
                signed_distance[i] = DOT_PRODUCT(dx,unit_np);
              }
              for (int il = 0; il <= level_max; ++il) {
                const double delta = hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<il);
                if (((signed_distance[0] > delta)&&(signed_distance[1] > delta)&&(signed_distance[2] > delta))||
                    ((signed_distance[0] < -delta)&&(signed_distance[1] < -delta)&&(signed_distance[2] < -delta)))
                  break;
                FOR_I3 {
                  const int ip1 = (i+1)%3;
                  // several conditions possible along this edge...
                  if (signed_distance[i] < -delta) {
                    if (signed_distance[ip1] < -delta) {
                      continue; // no action
                    }
                    else if (signed_distance[ip1] <= delta) {
                      // this edge contributes ip1 and an intersection at -delta. intersection first...
                      const double wip1 = (-delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                      double this_xpyp[2];
                      FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                      my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                      my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                      my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                      my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                      // then ip1...
                      my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[ip1][0]); // xp-min
                      my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[ip1][0]); // -xp-max
                      my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[ip1][1]); // yp-min
                      my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[ip1][1]); // -yp-max
                    }
                    else {
                      // the 2-intersection case...
                      double wip1 = (-delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                      double this_xpyp[2];
                      FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                      my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                      my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                      my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                      my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                      wip1 = (delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                      FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                      my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                      my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                      my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                      my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                    }
                  }
                  else if (signed_distance[i] < delta) {
                    if (signed_distance[ip1] < -delta) {
                      // this edge contributes i and an intersection at -delta. intersection first...
                      const double wip1 = (-delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                      double this_xpyp[2];
                      FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                      my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                      my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                      my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                      my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                      // then i...
                      my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[i][0]); // xp-min
                      my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[i][0]); // -xp-max
                      my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[i][1]); // yp-min
                      my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[i][1]); // -yp-max
                    }
                    else if (signed_distance[ip1] <= delta) {
                      // i and ip1, i first...
                      my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[i][0]); // xp-min
                      my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[i][0]); // -xp-max
                      my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[i][1]); // yp-min
                      my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[i][1]); // -yp-max
                      // then ip1...
                      my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[ip1][0]); // xp-min
                      my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[ip1][0]); // -xp-max
                      my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[ip1][1]); // yp-min
                      my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[ip1][1]); // -yp-max
                    }
                    else {
                      // this edge contributes i and an intersection at -delta. intersection first...
                      const double wip1 = (delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                      double this_xpyp[2];
                      FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                      my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                      my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                      my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                      my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                      // then i...
                      my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[i][0]); // xp-min
                      my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[i][0]); // -xp-max
                      my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[i][1]); // yp-min
                      my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[i][1]); // -yp-max
                    }
                  }
                  else {
                    if (signed_distance[ip1] < -delta) {
                      // the 2-intersection case...
                      double wip1 = (-delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                      double this_xpyp[2];
                      FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                      my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                      my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                      my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                      my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                      wip1 = (delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                      FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                      my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                      my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                      my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                      my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                    }
                    else if (signed_distance[ip1] <= delta) {
                      // this edge contributes i and an intersection at -delta. intersection first...
                      const double wip1 = (delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                      double this_xpyp[2];
                      FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                      my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                      my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                      my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                      my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                      // then ip1...
                      my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[ip1][0]); // xp-min
                      my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[ip1][0]); // -xp-max
                      my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[ip1][1]); // yp-min
                      my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[ip1][1]); // -yp-max
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // we do this in parallel, so no one rank has to loop on all tris...
    {
      const int triVec_size = triVec.size();
      int my_nst_avg = triVec_size/mpi_size;
      if (triVec_size%mpi_size) ++my_nst_avg;
      const int ist0 = min(triVec_size,mpi_rank*my_nst_avg);
      const int ist1 = min(triVec_size,(mpi_rank+1)*my_nst_avg);
      assert(ist1-ist0 <= my_nst_avg);
      for (int ist = ist0; ist < ist1; ++ist) {
        // calculate the bbox for this tri in e1,e2 space...
        double this_bbminmax[4] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };
        double xpyp[3][2];
        FOR_I3 {
          const double * const x = triVec[ist].x[i];
          const double dx[3] = DIFF(x,hcp_x0pdx0);
          // x-primed is in the hcp_e1 direction...
          xpyp[i][0] = DOT_PRODUCT(dx,hcp_e1);
          this_bbminmax[0] = min(this_bbminmax[0],xpyp[i][0]); // xp-min
          this_bbminmax[1] = min(this_bbminmax[1],-xpyp[i][0]); // -xp-max
          // and y-primed is in the hcp_e2 direction...
          xpyp[i][1] = DOT_PRODUCT(dx,hcp_e2);
          this_bbminmax[2] = min(this_bbminmax[2],xpyp[i][1]); // yp-min
          this_bbminmax[3] = min(this_bbminmax[3],-xpyp[i][1]); // -yp-max
        }
        // note: these tris cannot extend the level 0 bbox...
        //FOR_I4 my_bbminmax[0][i] = min(my_bbminmax[0][i],this_bbminmax[i]);
        assert(triVec[ist].window >= 0);
        const int st_level = (leowi!=NULL) ? leowi[triVec[ist].window]:0;
        assert((st_level >= 0)&&(st_level <= level_max));
        if (st_level >= 1) {
          assert(triVec[ist].window >= 0);
          //const int  st_nlayers = nlowi[triVec[ist].window];
          double rmax = 0.0;
          for (int il = 1; il <= st_level; ++il) {
            const int st_nlayers = nlowi[triVec[ist].window*(level_max+1)+il];
            //rmax += double(st_nlayers)*hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<il);
            rmax += double(st_nlayers)*hcp_delta/double(1<<il);
          }
          for (int il = 1; il <= st_level; ++il) {
            FOR_I4 my_bbminmax[il][i] = min(my_bbminmax[il][i],this_bbminmax[i]-rmax);
            const int st_nlayers = nlowi[triVec[ist].window*(level_max+1)+il];
            //rmax -= double(st_nlayers)*hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<il);
            rmax -= double(st_nlayers)*hcp_delta/double(1<<il);
          }
        }
        // if we have a plane constrain our rays to the thickened plane...
        if (plane_xp != NULL) {
          // at any level, add this tri's bbox to my_plane_bbminmax if it is within
          // delta for the level...
          // determine the signed distance from the plane for the corners of this tri...
          double signed_distance[3];
          FOR_I3 {
            const double * const x = triVec[ist].x[i];
            const double dx[3] = DIFF(x,plane_xp);
            signed_distance[i] = DOT_PRODUCT(dx,unit_np);
          }
          for (int il = 0; il <= level_max; ++il) {
            const double delta = hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<il);
            if (((signed_distance[0] > delta)&&(signed_distance[1] > delta)&&(signed_distance[2] > delta))||
                ((signed_distance[0] < -delta)&&(signed_distance[1] < -delta)&&(signed_distance[2] < -delta)))
              break;
            FOR_I3 {
              const int ip1 = (i+1)%3;
              // several conditions possible along this edge...
              if (signed_distance[i] < -delta) {
                if (signed_distance[ip1] < -delta) {
                  continue; // no action
                }
                else if (signed_distance[ip1] <= delta) {
                  // this edge contributes ip1 and an intersection at -delta. intersection first...
                  const double wip1 = (-delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                  double this_xpyp[2];
                  FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                  my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                  my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                  my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                  my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                  // then ip1...
                  my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[ip1][0]); // xp-min
                  my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[ip1][0]); // -xp-max
                  my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[ip1][1]); // yp-min
                  my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[ip1][1]); // -yp-max
                }
                else {
                  // the 2-intersection case...
                  double wip1 = (-delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                  double this_xpyp[2];
                  FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                  my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                  my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                  my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                  my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                  wip1 = (delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                  FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                  my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                  my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                  my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                  my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                }
              }
              else if (signed_distance[i] < delta) {
                if (signed_distance[ip1] < -delta) {
                  // this edge contributes i and an intersection at -delta. intersection first...
                  const double wip1 = (-delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                  double this_xpyp[2];
                  FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                  my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                  my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                  my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                  my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                  // then i...
                  my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[i][0]); // xp-min
                  my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[i][0]); // -xp-max
                  my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[i][1]); // yp-min
                  my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[i][1]); // -yp-max
                }
                else if (signed_distance[ip1] <= delta) {
                  // i and ip1, i first...
                  my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[i][0]); // xp-min
                  my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[i][0]); // -xp-max
                  my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[i][1]); // yp-min
                  my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[i][1]); // -yp-max
                  // then ip1...
                  my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[ip1][0]); // xp-min
                  my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[ip1][0]); // -xp-max
                  my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[ip1][1]); // yp-min
                  my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[ip1][1]); // -yp-max
                }
                else {
                  // this edge contributes i and an intersection at -delta. intersection first...
                  const double wip1 = (delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                  double this_xpyp[2];
                  FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                  my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                  my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                  my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                  my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                  // then i...
                  my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[i][0]); // xp-min
                  my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[i][0]); // -xp-max
                  my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[i][1]); // yp-min
                  my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[i][1]); // -yp-max
                }
              }
              else {
                if (signed_distance[ip1] < -delta) {
                  // the 2-intersection case...
                  double wip1 = (-delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                  double this_xpyp[2];
                  FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                  my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                  my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                  my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                  my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                  wip1 = (delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                  FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                  my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                  my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                  my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                  my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                }
                else if (signed_distance[ip1] <= delta) {
                  // this edge contributes i and an intersection at -delta. intersection first...
                  const double wip1 = (delta-signed_distance[i])/(signed_distance[ip1]-signed_distance[i]);
                  double this_xpyp[2];
                  FOR_J2 this_xpyp[j] = wip1*xpyp[ip1][j] + (1.0-wip1)*xpyp[i][j];
                  my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],this_xpyp[0]); // xp-min
                  my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-this_xpyp[0]); // -xp-max
                  my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],this_xpyp[1]); // yp-min
                  my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-this_xpyp[1]); // -yp-max
                  // then ip1...
                  my_plane_bbminmax[il][0] = min(my_plane_bbminmax[il][0],xpyp[ip1][0]); // xp-min
                  my_plane_bbminmax[il][1] = min(my_plane_bbminmax[il][1],-xpyp[ip1][0]); // -xp-max
                  my_plane_bbminmax[il][2] = min(my_plane_bbminmax[il][2],xpyp[ip1][1]); // yp-min
                  my_plane_bbminmax[il][3] = min(my_plane_bbminmax[il][3],-xpyp[ip1][1]); // -yp-max
                }
              }
            }
          }
        }
      }
    }

    double (*bbminmax)[4] = new double[level_max+1][4];
    MPI_Allreduce((double*)my_bbminmax,(double*)bbminmax,4*(level_max+1),MPI_DOUBLE,MPI_MIN,mpi_comm);
    delete[] my_bbminmax;

    double (*plane_bbminmax)[4] = new double[level_max+1][4];
    MPI_Allreduce((double*)my_plane_bbminmax,(double*)plane_bbminmax,4*(level_max+1),MPI_DOUBLE,MPI_MIN,mpi_comm);
    delete[] my_plane_bbminmax;

    // this is the first time that we have both the e0/e1/e2 space set, and hcp_delta set along with
    // the full surface bbox, use the largest bbox at level 0 to set delta_max, and the conversion from double to int...
    const double delta_max = max(-bbminmax[0][1]-bbminmax[0][0],-bbminmax[0][3]-bbminmax[0][2]);
    const int Nbase = int( log(hcp_delta*hcp_delta_factor*hcp_packing_factor/delta_max)/log(2.0) + hcp_bit_max );
    const double hcp_dxp = double(1<<Nbase)/(hcp_delta*hcp_delta_factor*hcp_packing_factor); // actually dj/dy == dk/dz
    const double one_o_hcp_dxp = 1.0/hcp_dxp;

    // warn users if any of the bounding boxes will blow the integer precision...
    for (int il = 1; il <= level_max; ++il) {
      // some levels could be empty
      if (bbminmax[il][0] == HUGE_VAL) {
        assert((bbminmax[il][1] == HUGE_VAL)&&(bbminmax[il][2] == HUGE_VAL)&&(bbminmax[il][3] == HUGE_VAL));
        continue;
      }
      double xy[2];
      xy[0] = bbminmax[il][0]; // xp-min
      xy[1] = bbminmax[il][2]; // yp-min
      int bbminL[2];
      setJKFrom2D(bbminL,xy,hcp_dxp);
      xy[0] = -bbminmax[il][1]; // xp-max
      xy[1] = -bbminmax[il][3]; // yp-max
      int bbmaxL[2];
      setJKFrom2D(bbmaxL,xy,hcp_dxp);
      // 2^31-1 = -2147483647
      //cout << bbminL[0] << " " << bbmaxL[0] << " " << bbminL[1] << " " << bbmaxL[1] << endl;
      if ((bbminL[0] == -2147483647)||(bbmaxL[0] == -2147483647)||(bbminL[1] == -2147483647)||(bbmaxL[1] == -2147483647)) {
        CWARN(" > level " << il << " bbox blows the integer precision. You may want to adjust your HCP_WINDOWs.");
      }
    }

    // finally, level 1..level_max bboxes should never extend beyond level 0 (the geometry)...
    for (int il = 1; il <= level_max; ++il) {
      FOR_I4 bbminmax[il][i] = max(bbminmax[il][i],bbminmax[0][i]);
    }

    if (mpi_rank == 0) {
      for (int il = 0; il <= level_max; ++il) {
        cout << " > full bbox for level: " << il << " " << bbminmax[il][0] << ":" << -bbminmax[il][1] << " "
             << bbminmax[il][2] << ":" << -bbminmax[il][3] << endl;
      }
      for (int il = 0; il <= level_max; ++il) {
        const double delta = hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<il);
        cout << " > delta for this level: " << delta << endl;
        if (plane_xp != NULL) {
          cout << " > plane bbox for level: " << il << " xp: " <<
            plane_bbminmax[il][0] << ":" << -plane_bbminmax[il][1] << " dxp: " << -plane_bbminmax[il][1]-plane_bbminmax[il][0] <<
            " yp: " << plane_bbminmax[il][2] << ":" << -plane_bbminmax[il][3] << " dyp: " << -plane_bbminmax[il][3]-plane_bbminmax[il][2] << endl;
        }
      }
    }

    // now apply the plane bbox to all full bbox's. These form the rays we need to consider at each level...

    if (plane_xp != NULL) {
      for (int il = 0; il <= level_max; ++il) {
        FOR_I4 bbminmax[il][i] = max(bbminmax[il][i],plane_bbminmax[il][i]);
      }
    }

    delete[] plane_bbminmax;

    // ultimately, the active rays and their in/out points are stored in a vector
    // of InOutData, one for each level...

    if (mpi_rank == 0) {
      double wtime = MPI_Wtime();
      cout << " bbox time: " << wtime-wtime_prev << " [s], total time so far: " << wtime-wtime0 << " [s]" << endl;
      wtime_prev = wtime;
    }

    // TODO: put this in the loop and clear at the end of each iteration

    vector<InOutData> * inoutVec = new vector<InOutData>[level_max+1];

    int8 * my_vv_count = new int8[level_max+1];
    for (int il = 0; il <= level_max; ++il)
      my_vv_count[il] = 0;

    // ==========================================================================
    // now loop down from the finest level and do all the work to introduce
    // points into vertexVec at each level...
    //
    // This is basically 2 steps for each level:
    // 1. identify exact intersections with the surface: called level=0
    // 2. add in/out points associated with refinement regions (e.g. bloated tris)
    //
    // ==========================================================================

    vertexVec.clear();
    for (int level = level_max; level >= 0; --level) {

      if (mpi_rank == 0)
	cout << "working on level " << level << "..." << endl;

      double xy[2];
      xy[0] = bbminmax[level][0]; // xp-min
      xy[1] = bbminmax[level][2]; // yp-min
      int bbminL[2];
      setJKFrom2D(bbminL,xy,hcp_dxp);
      xy[0] = -bbminmax[level][1]; // xp-max
      xy[1] = -bbminmax[level][3]; // yp-max
      int bbmaxL[2];
      setJKFrom2D(bbmaxL,xy,hcp_dxp);

      // now loop rays and start building points...

      const int djk = (1<<(Nbase-level));
      int jminL = ( bbminL[0] < 0 ? bbminL[0]-(bbminL[0]%djk) : bbminL[0]-(bbminL[0]%djk)+djk );  assert(jminL > bbminL[0]);
      int jmaxL = ( bbmaxL[0] < 0 ? bbmaxL[0]-(bbmaxL[0]%djk)-djk : bbmaxL[0]-(bbmaxL[0]%djk) );  assert(jmaxL <= bbmaxL[0]); // added == because zero is possible
      int kminL = ( bbminL[1] < 0 ? bbminL[1]-(bbminL[1]%djk) : bbminL[1]-(bbminL[1]%djk)+djk );  assert(kminL > bbminL[1]);
      int kmaxL = ( bbmaxL[1] < 0 ? bbmaxL[1]-(bbmaxL[1]%djk)-djk : bbmaxL[1]-(bbmaxL[1]%djk) );  assert(kmaxL <= bbmaxL[1]);

      int8 nray_cart = 0;
      int8 nray = 0;
      for (int j = jminL; j <= jmaxL; j += djk) {
        for (int k = kminL; k <= kmaxL; k += djk) {
          ++nray_cart;
          if ((level == 0)||(j%(2*djk) != 0)||(k%(2*djk) != 0)) {
            // j and k are even...
            assert(j%2 == 0);
            assert(k%2 == 0);
            ++nray;
          }
        }
      }
      if (mpi_rank == 0) cout << " > cart ray array size: " << (jmaxL-jminL)/djk+1 << " " << (kmaxL-kminL)/djk+1 << endl;
      if (mpi_rank == 0) cout << " > nray cart: " << nray_cart << " nray: " << nray << " ratio: " << double(nray)/double(nray_cart) << endl;
      if (nray == 0)
        continue;

      // ======================================================================================
      // which rays at this level pass through the surface geometry...
      // ======================================================================================

      int mcount = 0;
      int pcount = 0;
      vector<uint8> indexVec; // for storing the rays that are "on"
      vector<double> xpVec;
      if (send_count == NULL) send_count = new int[mpi_size];
      FOR_RANK send_count[rank] = 0;

      // stride through tris at mpi_size...
      for (int ipart = 0; ipart < PartData::partVec.size(); ++ipart) {
        // -------------------------------
        // every surface needs an st_flag to map to potential windows...
        // -------------------------------
        if (PartData::partVec[ipart]->surface) {
          SurfaceShm * surface = PartData::partVec[ipart]->surface;
          for (int ist = mpi_rank; ist < surface->nst; ist += mpi_size) {
            // ignore ist's blanked by other parts' ff_surfaces...
            int st_type;
            int ipart2;
            if (PartData::partVec[ipart]->getPartForSt(ipart2,ist)) {
              st_type = -2; // tri is shielded by another part. Do not even consider
            }
            else {
              // this tri may be part of a SURFACE_FAZONE_FF...
              if (PartData::partVec[ipart]->surface_zn_flag) {
                assert(PartData::partVec[ipart]->ff_type == PartData::SURFACE_FAZONE_FF);
                const int isz = surface->znost[ist];
                if (PartData::partVec[ipart]->surface_zn_flag[isz] == -1)
                  st_type = -1; // this tri is part of a FF zone
                else {
                  assert(PartData::partVec[ipart]->surface_zn_flag[isz] == 0);
                  st_type = 0; // this tri is a regular boundary tri
                }
              }
              else {
                st_type = 0; // this tri is a regular boundary tri
              }
            }
            // only consder tris that are either ff (st_type = -1), or regular boundary (st_type = 0)...
            if ((st_type == -1)||(st_type == 0)) {
              // we need an integer version of this tri...
              const double * const xsp0 = surface->xsp[surface->spost[ist][0]];
              const double * const xsp1 = surface->xsp[surface->spost[ist][1]];
              const double * const xsp2 = surface->xsp[surface->spost[ist][2]];
              int jk0[2]; setJKFrom3D(jk0,xsp0,hcp_dxp);
              int jk1[2]; setJKFrom3D(jk1,xsp1,hcp_dxp);
              int jk2[2]; setJKFrom3D(jk2,xsp2,hcp_dxp);
              int bbmin[2]; FOR_I2 bbmin[i] = min(jk0[i],min(jk1[i],jk2[i]));
              int bbmax[2]; FOR_I2 bbmax[i] = max(jk0[i],max(jk1[i],jk2[i]));
              int jmin = max( jminL, ( bbmin[0] < 0 ? bbmin[0]-(bbmin[0]%djk)     : bbmin[0]-(bbmin[0]%djk)+djk ));  assert(jmin > bbmin[0]);
              int jmax = min( jmaxL, ( bbmax[0] < 0 ? bbmax[0]-(bbmax[0]%djk)-djk : bbmax[0]-(bbmax[0]%djk) ));      assert(jmax <= bbmax[0]);
              int kmin = max( kminL, ( bbmin[1] < 0 ? bbmin[1]-(bbmin[1]%djk)     : bbmin[1]-(bbmin[1]%djk)+djk ));  assert(kmin > bbmin[1]);
              int kmax = min( kmaxL, ( bbmax[1] < 0 ? bbmax[1]-(bbmax[1]%djk)-djk : bbmax[1]-(bbmax[1]%djk) ));      assert(kmax <= bbmax[1]);
              for (int j = jmin; j <= jmax; j += djk) {
                for (int k = kmin; k <= kmax; k += djk) {
                  if ((level == 0)||(j%(2*djk) != 0)||(k%(2*djk) != 0)) {
                    // j and k are even...
                    assert(j%2 == 0);
                    assert(k%2 == 0);
                    const int8 dx0[2] = { int8( jk0[0] - j ), int8( jk0[1] - k ) };
                    const int8 dx1[2] = { int8( jk1[0] - j ), int8( jk1[1] - k ) };
                    const int8 dx2[2] = { int8( jk2[0] - j ), int8( jk2[1] - k ) };
                    const int8 V0 = dx1[0]*dx2[1] - dx1[1]*dx2[0];
                    const int8 V1 = dx2[0]*dx0[1] - dx2[1]*dx0[0];
                    const int8 V2 = dx0[0]*dx1[1] - dx0[1]*dx1[0];
                    if ((V0 != 0)||(V1 != 0)||(V2 != 0)) {
                      if ((V0 >= 0)&&(V1 >= 0)&&(V2 >= 0)) {
                        // build the intersection "x" along the e0 direction...
                        const double xp = ( double(V0)*((xsp0[0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp0[1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp0[2]-hcp_x0pdx0[2])*hcp_e0[2]) +
                                            double(V1)*((xsp1[0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp1[1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp1[2]-hcp_x0pdx0[2])*hcp_e0[2]) +
                                            double(V2)*((xsp2[0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp2[1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp2[2]-hcp_x0pdx0[2])*hcp_e0[2]) )/double(V0+V1+V2);
                        xpVec.push_back(xp);
                        const uint8 index = uint8((j-jminL)/djk)*uint8((kmaxL-kminL)/djk+1)+uint8((k-kminL)/djk);
                        //assert(index >= 0);
                        //assert(index < TWO_BILLION/2); // need space for 2 bits

                        // TODO check: remove this eventually...
                        {
                          const int j_check = index/((kmaxL-kminL)/djk+1)*djk+jminL;
                          const int k_check = (index-((j-jminL)/djk)*((kmaxL-kminL)/djk+1))*djk+kminL;
                          assert(j_check == j);
                          assert(k_check == k);
                        }

                        const int rank = index%uint8(mpi_size);
                        assert((rank >=0) && (rank < mpi_size));
                        ++send_count[rank];
                        if ((V0 == 0)||(V1 == 0)||(V2 == 0)) {
                          // this is an edge intersection. Because of the staggered spaces, it cannot
                          // be a vertex intersection. Here we add 1 because there should be another.  Note
                          // that we can assert that only one of the volumes is zero -- i.e. we are on an edge...
                          assert( ((V0 > 0)&&(V1 > 0)) || ((V0 > 0)&&(V2 > 0)) || ((V1 > 0)&&(V2 > 0)) );
                          if (st_type == 0) {
                            pcount += 1;
                            indexVec.push_back(index<<2);     // i.e. 00
                            assert((indexVec.back()>>2) == index);
                          }
                          else {
                            assert(st_type == -1);
                            mcount += 1;
                            indexVec.push_back((index<<2)|2ll); // i.e. 10
                            assert((indexVec.back()>>2) == index);
                          }
                        }
                        else {
                          // all positive -- intersection is good and gets a 2...
                          if (st_type == 0) {
                            pcount += 2;
                            indexVec.push_back((index<<2)|1ll); // i.e. 01
                            assert((indexVec.back()>>2) == index);
                          }
                          else {
                            assert(st_type == -1);
                            mcount += 2;
                            indexVec.push_back((index<<2)|3ll); // i.e. 11
                            assert((indexVec.back()>>2) == index);
                          }
                        }
                      }
                      else if ((V0 <= 0)&&(V1 <= 0)&&(V2 <= 0)) {
                        const double xp = ( double(V0)*((xsp0[0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp0[1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp0[2]-hcp_x0pdx0[2])*hcp_e0[2]) +
                                            double(V1)*((xsp1[0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp1[1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp1[2]-hcp_x0pdx0[2])*hcp_e0[2]) +
                                            double(V2)*((xsp2[0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp2[1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp2[2]-hcp_x0pdx0[2])*hcp_e0[2]) )/double(V0+V1+V2);
                        xpVec.push_back(xp);
                        const uint8 index = uint8((j-jminL)/djk)*uint8((kmaxL-kminL)/djk+1)+uint8((k-kminL)/djk);
                        //assert(index >= 0);
                        //assert(index < TWO_BILLION/2); // need space for 2 bits

                        // TODO check: remove this eventually...
                        {
                          const int j_check = index/((kmaxL-kminL)/djk+1)*djk+jminL;
                          const int k_check = (index-((j-jminL)/djk)*((kmaxL-kminL)/djk+1))*djk+kminL;
                          assert(j_check == j);
                          assert(k_check == k);
                        }

                        const int rank = index%uint8(mpi_size);
                        assert((rank >=0) && (rank < mpi_size));
                        ++send_count[rank];
                        if ((V0 == 0)||(V1 == 0)||(V2 == 0)) {
                          // see note above...
                          assert( ((V0 < 0)&&(V1 < 0)) || ((V0 < 0)&&(V2 < 0)) || ((V1 < 0)&&(V2 < 0)) );
                          if (st_type == 0) {
                            mcount += 1;
                            indexVec.push_back((index<<2)|2ll); // i.e. 10
                            assert((indexVec.back()>>2) == index);
                          }
                          else {
                            assert(st_type == -1);
                            pcount += 1;
                            indexVec.push_back((index<<2)); // i.e. 00
                            assert((indexVec.back()>>2) == index);
                          }
                        }
                        else {
                          // all negative...
                          if (st_type == 0) {
                            mcount += 2;
                            indexVec.push_back((index<<2)|3ll); // i.e. 11
                            assert((indexVec.back()>>2) == index);
                          }
                          else {
                            assert(st_type == -1);
                            pcount += 2;
                            indexVec.push_back((index<<2)|1ll); // i.e. 01
                            assert((indexVec.back()>>2) == index);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        // -------------------------------
        // and ff_surface's can also turn rays on/off...
        // note we flip direction compared to above.
        // -------------------------------
        if (PartData::partVec[ipart]->ff_surface) {
          SurfaceShm * ff_surface = PartData::partVec[ipart]->ff_surface;
          for (int ist = mpi_rank; ist < ff_surface->nst; ist += mpi_size) {
            // consider all periodic version of this tri...
            // HACK: consider just the first (untransformed) version of the ff_surface for this part
            //for (int iper = 0; iper < periodicRtVec.size(); ++iper) 
            int iper = 0;
            {
              double xsp_t[3][3];
              FOR_I3 FOR_J3 xsp_t[i][j] = ff_surface->xsp[ff_surface->spost[ist][i]][j];
              if (periodicRtVec[iper].b_t) {
                FOR_I3 FOR_J3 xsp_t[i][j] += periodicRtVec[iper].t[j];
              }
              if (periodicRtVec[iper].b_R) {
                //assert(0);
                double tmp[3];
                FOR_I3 {
                  FOR_J3 tmp[j] = xsp_t[i][j];
                  MiscUtils::applyRotation(xsp_t[i],periodicRtVec[iper].R,tmp);
                }
              }
              // we need an integer version of this tri...
              //const double * const xsp0 = ff_surface->xsp[ff_surface->spost[ist][0]];
              //const double * const xsp1 = ff_surface->xsp[ff_surface->spost[ist][1]];
              //const double * const xsp2 = ff_surface->xsp[ff_surface->spost[ist][2]];
              int jk0[2]; setJKFrom3D(jk0,xsp_t[0],hcp_dxp);
              int jk1[2]; setJKFrom3D(jk1,xsp_t[1],hcp_dxp);
              int jk2[2]; setJKFrom3D(jk2,xsp_t[2],hcp_dxp);
              int bbmin[2]; FOR_I2 bbmin[i] = min(jk0[i],min(jk1[i],jk2[i]));
              int bbmax[2]; FOR_I2 bbmax[i] = max(jk0[i],max(jk1[i],jk2[i]));
              int jmin = max( jminL, ( bbmin[0] < 0 ? bbmin[0]-(bbmin[0]%djk)     : bbmin[0]-(bbmin[0]%djk)+djk ));  assert(jmin > bbmin[0]);
              int jmax = min( jmaxL, ( bbmax[0] < 0 ? bbmax[0]-(bbmax[0]%djk)-djk : bbmax[0]-(bbmax[0]%djk) ));      assert(jmax <= bbmax[0]);
              int kmin = max( kminL, ( bbmin[1] < 0 ? bbmin[1]-(bbmin[1]%djk)     : bbmin[1]-(bbmin[1]%djk)+djk ));  assert(kmin > bbmin[1]);
              int kmax = min( kmaxL, ( bbmax[1] < 0 ? bbmax[1]-(bbmax[1]%djk)-djk : bbmax[1]-(bbmax[1]%djk) ));      assert(kmax <= bbmax[1]);
              for (int j = jmin; j <= jmax; j += djk) {
                for (int k = kmin; k <= kmax; k += djk) {
                  if ((level == 0)||(j%(2*djk) != 0)||(k%(2*djk) != 0)) {
                    // j and k are even...
                    assert(j%2 == 0);
                    assert(k%2 == 0);
                    const int8 dx0[2] = { int8( jk0[0] - j ), int8( jk0[1] - k ) };
                    const int8 dx1[2] = { int8( jk1[0] - j ), int8( jk1[1] - k ) };
                    const int8 dx2[2] = { int8( jk2[0] - j ), int8( jk2[1] - k ) };
                    const int8 V0 = dx1[0]*dx2[1] - dx1[1]*dx2[0];
                    const int8 V1 = dx2[0]*dx0[1] - dx2[1]*dx0[0];
                    const int8 V2 = dx0[0]*dx1[1] - dx0[1]*dx1[0];
                    if ((V0 != 0)||(V1 != 0)||(V2 != 0)) {
                      if ((V0 >= 0)&&(V1 >= 0)&&(V2 >= 0)) {
                        // build the intersection "x" along the e0 direction...
                        const double xp = ( double(V0)*((xsp_t[0][0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp_t[0][1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp_t[0][2]-hcp_x0pdx0[2])*hcp_e0[2]) +
                                            double(V1)*((xsp_t[1][0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp_t[1][1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp_t[1][2]-hcp_x0pdx0[2])*hcp_e0[2]) +
                                            double(V2)*((xsp_t[2][0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp_t[2][1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp_t[2][2]-hcp_x0pdx0[2])*hcp_e0[2]) )/double(V0+V1+V2);
                        xpVec.push_back(xp);
                        const uint8 index = uint8((j-jminL)/djk)*uint8((kmaxL-kminL)/djk+1)+uint8((k-kminL)/djk);
                        //assert(index >= 0);
                        //assert(index < TWO_BILLION/2); // need space for 2 bits

                        // TODO check: remove this eventually...
                        {
                          const int j_check = index/((kmaxL-kminL)/djk+1)*djk+jminL;
                          const int k_check = (index-((j-jminL)/djk)*((kmaxL-kminL)/djk+1))*djk+kminL;
                          assert(j_check == j);
                          assert(k_check == k);
                        }

                        const int rank = index%uint8(mpi_size);
                        assert((rank >=0) && (rank < mpi_size));
                        ++send_count[rank];
                        if ((V0 == 0)||(V1 == 0)||(V2 == 0)) {
                          // this is an edge intersection. Because of the staggered spaces, it cannot
                          // be a vertex intersection. Here we add 1 because there should be another.  Note
                          // that we can assert that only one of the volumes is zero -- i.e. we are on an edge...
                          assert( ((V0 > 0)&&(V1 > 0)) || ((V0 > 0)&&(V2 > 0)) || ((V1 > 0)&&(V2 > 0)) );
                          mcount += 1;
                          indexVec.push_back((index<<2)|2ll); // i.e. 10
                          assert((indexVec.back()>>2) == index);
                        }
                        else {
                          // all positive -- intersection is good and gets a 2...
                          mcount += 2;
                          indexVec.push_back((index<<2)|3ll); // i.e. 11
                          assert((indexVec.back()>>2) == index);
                        }
                      }
                      else if ((V0 <= 0)&&(V1 <= 0)&&(V2 <= 0)) {
                        const double xp = ( double(V0)*((xsp_t[0][0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp_t[0][1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp_t[0][2]-hcp_x0pdx0[2])*hcp_e0[2]) +
                                            double(V1)*((xsp_t[1][0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp_t[1][1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp_t[1][2]-hcp_x0pdx0[2])*hcp_e0[2]) +
                                            double(V2)*((xsp_t[2][0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp_t[2][1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp_t[2][2]-hcp_x0pdx0[2])*hcp_e0[2]) )/double(V0+V1+V2);
                        xpVec.push_back(xp);
                        const uint8 index = uint8((j-jminL)/djk)*uint8((kmaxL-kminL)/djk+1)+uint8((k-kminL)/djk);
                        //assert(index >= 0);
                        //assert(index < TWO_BILLION/2); // need space for 2 bits

                        // TODO check: remove this eventually...
                        {
                          const int j_check = index/((kmaxL-kminL)/djk+1)*djk+jminL;
                          const int k_check = (index-((j-jminL)/djk)*((kmaxL-kminL)/djk+1))*djk+kminL;
                          assert(j_check == j);
                          assert(k_check == k);
                        }

                        const int rank = index%uint8(mpi_size);
                        assert((rank >=0) && (rank < mpi_size));
                        ++send_count[rank];
                        if ((V0 == 0)||(V1 == 0)||(V2 == 0)) {
                          // see note above...
                          assert( ((V0 < 0)&&(V1 < 0)) || ((V0 < 0)&&(V2 < 0)) || ((V1 < 0)&&(V2 < 0)) );
                          pcount += 1;
                          indexVec.push_back((index<<2)); // i.e. 00
                          assert((indexVec.back()>>2) == index);
                        }
                        else {
                          // all negative...
                          pcount += 2;
                          indexVec.push_back((index<<2)|1ll); // i.e. 01
                          assert((indexVec.back()>>2) == index);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        } // ff_surface
      }

      assert(indexVec.size() == xpVec.size());

      if (send_disp == NULL) send_disp = new int[mpi_size];
      send_disp[0] = 0;
      int my_sum = send_count[0];
      for (int rank = 1; rank < mpi_size; ++rank) {
        my_sum += send_count[rank];
	send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
      }
      const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
      assert(send_count_sum == indexVec.size());

      uint8 * send_buf_uint8 = new uint8[send_count_sum];
      double * send_buf_double = new double[send_count_sum];
      for (int isend = 0; isend < send_count_sum; ++isend) {
        // to get the index back, shift down 2 bits...
        const uint8 index = (indexVec[isend]>>2);
        const int rank = index%uint8(mpi_size);
        send_buf_uint8[send_disp[rank]] = indexVec[isend];
        send_buf_double[send_disp[rank]] = xpVec[isend];
        ++send_disp[rank];
      }
      indexVec.clear();
      xpVec.clear();

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
	send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      // set up recv-side stuff...

      if (recv_count == NULL) recv_count = new int[mpi_size];
      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      if (recv_disp == NULL) recv_disp = new int[mpi_size];
      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
	recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      uint8 * recv_buf_uint8 = new uint8[recv_count_sum];
      //for (int irecv = 0; irecv < recv_count_sum; ++irecv)
      //  recv_buf_uint8[irecv] = ULONG_MAX;
      MPI_Alltoallv(send_buf_uint8,send_count,send_disp,MPI_UINT8,
                    recv_buf_uint8,recv_count,recv_disp,MPI_UINT8,mpi_comm);
      delete[] send_buf_uint8; send_buf_uint8 = NULL;
      //for (int irecv = 0; irecv < recv_count_sum; ++irecv)
      //  assert(recv_buf_uint8[irecv] < ULONG_MAX);

      double * recv_buf_double = new double[recv_count_sum];
      //for (int irecv = 0; irecv < recv_count_sum; ++irecv)
      //  recv_buf_double[irecv] = HUGE_VAL;
      MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
		    recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
      delete[] send_buf_double; send_buf_double = NULL;
      //for (int irecv = 0; irecv < recv_count_sum; ++irecv)
      //  assert(recv_buf_double[irecv] < HUGE_VAL);

      // recall we have rays that satisfy index%mpi_size == mpi_rank

      const int my_ray_flag_size = uint8((jmaxL-jminL)/djk+1)*uint8((kmaxL-kminL)/djk+1)/uint8(mpi_size)+1ll;

      //int my_buf[3] = { mcount, pcount, my_nray_final };
      int8 my_buf[3] = { mcount, pcount, my_ray_flag_size };
      int8 buf[3];
      MPI_Reduce(my_buf,buf,3,MPI_INT8,MPI_SUM,0,mpi_comm);

      int my_buf_max[1] = { my_ray_flag_size };
      int buf_max[1];
      MPI_Reduce(my_buf_max,buf_max,1,MPI_INT,MPI_MAX,0,mpi_comm);

      if (mpi_rank == 0) {
        cout << " > level " << level << " nray_final: " << buf[2] << " avg/max per rank: " << int(double(buf[2])/double(mpi_size)) << " " << buf_max[0] << endl;
        //if (buf[0] != buf[1]) cout << buf[0] << " " << buf[1];
        assert(buf[0] == buf[1]); // mcount == pcount
	double wtime = MPI_Wtime();
        cout << " > surface ray intersection time: " << wtime-wtime_prev << " [s], total time so far: " << wtime-wtime0 << " [s]" << endl;
        wtime_prev = wtime;
      }

      // now pack our rays into inoutVec[level]. Note that we use a flag
      // here to figure out where to insert the rays so we only need to
      // sort within the group...

      int * my_ray_flag = new int[my_ray_flag_size];
      for (int ii = 0; ii < my_ray_flag_size; ++ii)
        my_ray_flag[ii] = 0;

      for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
        const uint8 index = (recv_buf_uint8[irecv]>>2);
        assert((index-uint8(mpi_rank))%uint8(mpi_size) == 0ll);
        const int ii = (index-uint8(mpi_rank))/uint8(mpi_size);
        assert((ii >= 0)&&(ii < my_ray_flag_size));
        ++my_ray_flag[ii];
      }

      int count = 0;
      for (int ii = 0; ii < my_ray_flag_size; ++ii) {
        const int tmp = my_ray_flag[ii];
        my_ray_flag[ii] = count;
        count += tmp;
      }
      assert(count == recv_count_sum);

      assert(inoutVec[level].empty());
      inoutVec[level].resize(count);

      for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
        const uint8 index = (recv_buf_uint8[irecv]>>2);
        assert((index-uint8(mpi_rank))%uint8(mpi_size) == 0ll);
        const int ii = (index-uint8(mpi_rank))/uint8(mpi_size);
        assert((ii >= 0)&&(ii < my_ray_flag_size));
        const int iv = my_ray_flag[ii]++;
        inoutVec[level][iv].xp = recv_buf_double[irecv];
        inoutVec[level][iv].index = index;
        const int bits = (recv_buf_uint8[irecv]&3ll);
        if (bits == 0) {
          // one "out" (associated with pcount) -- must have been an edge intersection...
          inoutVec[level][iv].level = -1;
        }
        else if (bits == 1) {
          // two "out" -- regular triangle intersection...
          inoutVec[level][iv].level = -2;
        }
        else if (bits == 2) {
          // one "in" (associated with mcount) -- edge...
          inoutVec[level][iv].level = 1;
        }
        else {
          assert(bits == 3);
          // two "in"
          inoutVec[level][iv].level = 2;
        }
      }

      delete[] recv_buf_uint8; recv_buf_uint8 = NULL;
      delete[] recv_buf_double; recv_buf_double = NULL;

      // and condense the inoutVec: this is a special condense
      // the first time because of the +/- 1,2 stuff associated with these
      // exact intersections. Once condensed, leave the in/out points as simply
      // +1/-1... i.e. use a level+1 strategy

      static int debug_count = 0;
      const int this_debug_rank = getIntParam("DEBUG_RANK",-1);
      const int this_debug_count = getIntParam("DEBUG_COUNT",-1);
      
      vector<InOutData>::iterator iter_end = inoutVec[level].begin();
      while (iter_end != inoutVec[level].end()) {
        vector<InOutData>::iterator iter_begin = iter_end++;
        while ((iter_end != inoutVec[level].end())&&(iter_end->index == iter_begin->index))
          ++iter_end;
        sort(iter_begin,iter_end);

        ++debug_count;

        // get iters within an "0" length region and see if there sum is 0. If so, make all of them 0
        /*
          vector<InOutData>::iterator iter = iter_begin;
          while (iter != iter_end) {
          int sum = iter->level;
          vector<InOutData>::iterator iter_next = iter; iter_next++;
          while ((iter_next != iter_end)&&(fabs(iter->xp-iter_next->xp) < 2.0*one_o_hcp_dxp)) {
          sum += iter_next->level;
          ++iter_next;
          }
          if (sum == 0) {
          vector<InOutData>::iterator iter2 = iter;
          while (iter2 != iter_next) {
          iter2->level = 0;
          ++iter2;
          }
          }
          iter = iter_next;
          }
        */
        
        if ((mpi_rank == this_debug_rank)&&(debug_count == this_debug_count)) {
          cout << "XXXXXXXXXX see file: debug.dat" << endl;
          const int j = iter_begin->index/((kmaxL-kminL)/djk+1)*djk + jminL;
          const int k = (iter_begin->index-((j-jminL)/djk)*((kmaxL-kminL)/djk+1))*djk + kminL;
          const double y = hcp_x0pdx0[1] + double(j)*one_o_hcp_dxp;
          const double z = hcp_x0pdx0[2] + double(k)*one_o_hcp_dxp;
          FILE * fp = fopen("debug.dat","w");
          for (vector<InOutData>::iterator iter = iter_begin; iter != iter_end; ++iter) {
            fprintf(fp,"%18.15le %18.15le %18.15le %llu %d\n",iter->xp+hcp_x0pdx0[0],y,z,iter->index,iter->level);
          }
          fclose(fp);
        }
        
        //for (vector<InOutData>::iterator iter2 = iter_begin; iter2 != iter_end; ++iter2)
        //  cout << iter2->level << " ";
        //cout << endl; cout.flush();

        vector<InOutData>::iterator iter_next = iter_begin;
        //vector<InOutData>::iterator iter_prev = iter_end;
        //vector<InOutData>::iterator iter_prev2 = iter_end;
        int level_sum = 0;
        while (iter_next != iter_end) {
          vector<InOutData>::iterator iter = iter_next; iter_next++;
          level_sum += iter->level;
          //if (iter_next != iter_end) cout << fabs(iter->xp-iter_next->xp) << " ";
          if ((iter_next != iter_end)&&(iter->level*iter_next->level == -1)) {
            //if ((iter_next != iter_end)&&(iter->level*iter_next->level == -1)&&(fabs(iter->xp-iter_next->xp) < 2.0*one_o_hcp_dxp)) { 
            //if ((iter_next != iter_end)&&((iter->level*iter_next->level == -1)||(iter->level*iter_next->level == -4))&&(fabs(iter->xp-iter_next->xp) < 2.0*one_o_hcp_dxp)) { // changed to go coplanar faces
            // add check to make sure we are not canceling a 1 -1 when the 1 was needed to complete a "face" intersection (1 1 -> 2).
            // this occured when the ff was on top of the surface...
            //if ((iter_prev == iter_end)||(iter_prev->level*iter->level != 1)) {
            bool b_zero = false;
            if (abs(level_sum)%2 == 1) {
              level_sum += iter_next->level;
              b_zero = true;
            }
            //if (iter_prev == iter_end) {
            // / assert(iter_prev2 == iter_end);
            //  b_zero = true;
           // }
            //else {
            //  if (iter_prev->level*iter->level != 1) {
            //    b_zero = true;
            //  }
            //  else {
            //    if ((iter_prev2 != iter_end)&&(iter_prev2->level*iter_prev->level == 1)) {
            //      b_zero = true;
            //    }
            //  }
            //}
            if (b_zero) {
              iter->level = 0;
              iter_next->level = 0;
              iter_next++; // skip to next iter,iter_next pair
            }
          }
          //iter_prev2 = iter_prev;
          //iter_prev = iter;
        }
        //cout << endl;

        //for (vector<InOutData>::iterator iter2 = iter_begin; iter2 != iter_end; ++iter2)
        //  cout << iter2->level << " ";
        //cout << endl; cout.flush();

        int level_count0 = 0;
        bool on_flag = false;

        if ((mpi_rank == this_debug_rank)&&(debug_count == this_debug_count)) {
          cout << "XXXXXXXXXX see file: debug2.dat" << endl;
          const int j = iter_begin->index/((kmaxL-kminL)/djk+1)*djk + jminL;
          const int k = (iter_begin->index-((j-jminL)/djk)*((kmaxL-kminL)/djk+1))*djk + kminL;
          const double y = hcp_x0pdx0[1] + double(j)*one_o_hcp_dxp;
          const double z = hcp_x0pdx0[2] + double(k)*one_o_hcp_dxp;
          FILE * fp = fopen("debug2.dat","w");
          for (vector<InOutData>::iterator iter = iter_begin; iter != iter_end; ++iter) {
            fprintf(fp,"%18.15le %18.15le %18.15le %d %d\n",iter->xp+hcp_x0pdx0[0],y,z,iter->index,iter->level);
          }
          fclose(fp);
        }

        for (vector<InOutData>::iterator iter = iter_begin; iter != iter_end; ++iter) {
          //cout << on_flag << " ";
          if (iter->level == 0) continue;
          level_count0 += iter->level;
          if (level_count0 == 2) {
            // keep this one. assume it is the "in" point...
            if (!(iter->level > 0)) cout << "ERROR: mpi_rank: " << mpi_rank << "  debug_count: " << debug_count << endl; // must have been +1 or +2
            assert(iter->level > 0); // must have been +1 or +2
            if (!(on_flag == false)) cout << "ERROR: mpi_rank: " << mpi_rank << "  debug_count: " << debug_count << endl; // must have been +1 or +2
            assert(on_flag == false);
            iter->level = 1; // force 1
            on_flag = true;
          }
          else if (level_count0 == 0) {
            // if we were positive (i.e. on) then this shuts off the ray...
            if (on_flag) {
              assert(iter->level < 0);
              iter->level = -1; // force -1
              on_flag = false;
            }
            else {
              // we returned to level_count0 == 0 from below, so do not turn on...
              //assert(iter->level > 0); // failed for 1 -2 1 -1 2 -1 (never turned on)
              assert((iter->level > 0)||(iter->level == -1));
              iter->level = 0;
            }
          }
          else {
            // anything else we can discard...
            iter->level = 0;
          }
          //cout << endl; cout.flush();
        }
        if (level_count0 != 0) {
          cout << "About to fail in/out balance in surface intersection: --DEBUG_RANK " << mpi_rank << " --DEBUG_COUNT " << debug_count << endl;
          cout.flush();
        }
        assert(level_count0 == 0);
        // just make sure we are off...
        if (on_flag == true) {
          cout << "About to fail in/out balance in surface intersection: --DEBUG_RANK " << mpi_rank << " --DEBUG_COUNT " << debug_count << endl;
          cout.flush();
        }
        assert(on_flag == false);
      }

      // now all the zeros should be removed...

      int new_size = 0;
      for (int ii = 0, limit = inoutVec[level].size(); ii < limit; ++ii) {
        if (inoutVec[level][ii].level != 0) {
          inoutVec[level][new_size] = inoutVec[level][ii];
          ++new_size;
        }
      }
      inoutVec[level].resize(new_size);

      // =======================================================
      // now we have the possible rays at this level. Remaining
      // code is to build all possible transitions in refinement
      // (i.e. bloated tris) as well as other details like
      // other refinement windows and their bloated tris...
      // =======================================================

      int8 my_window_count[3] = { 0, 0, 0 };
      int exchange_count = 0;

      // loop through surface tris round-robin style...
      int threshold_windowVec_size = 2000000; // tried 1M-10M, 2M seems pretty good.
      vector<InOutData> windowVec;

      int round_robin_rank = 0;
      for (int ipart = 0; ipart < PartData::partVec.size(); ++ipart) {

        // --------------------------------------------------------------
        // every surface needs an st_flag to map to potential windows...
        // --------------------------------------------------------------
        if (PartData::partVec[ipart]->surface) {
          SurfaceShm * surface = PartData::partVec[ipart]->surface;
          for (int ist = 0; ist < surface->nst; ++ist) {
            int st_level = 0;
            int iwindow = -1;
            if (PartData::partVec[ipart]->getWindowForSt(iwindow,ist))
              st_level = (leowi!=NULL) ? leowi[iwindow]:0;
            if (st_level >= max(level,1)) {
              assert(iwindow >= 0);
              // this tri can affect this level...
              if (round_robin_rank == mpi_rank) {

                ++my_window_count[1];

                // calculate the maximum bloat...
                double rmax = 0.0;
                for (int il = max(level,1); il <= st_level; ++il) {
                  const int st_nlayers = nlowi[iwindow*(level_max+1)+il];
                  //rmax += double(st_nlayers)*hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<il);
                  rmax += double(st_nlayers)*hcp_delta/double(1<<il);
                }
                const double delta = hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<level);

                // if the tri+bloat is outside the thickened plane (if it exists), we can skip it...
                double signed_distance[3];
                if (plane_xp != NULL) {
                  FOR_I3 {
                    const double * const x = surface->xsp[surface->spost[ist][i]];
                    const double dx[3] = DIFF(x,plane_xp);
                    signed_distance[i] = DOT_PRODUCT(dx,unit_np);
                  }
                }

                if ((plane_xp == NULL)||
                    (!(((signed_distance[0] > (delta+rmax))&&(signed_distance[1] > (delta+rmax))&&(signed_distance[2] > (delta+rmax)))||
                       ((signed_distance[0] < -(delta+rmax))&&(signed_distance[1] < -(delta+rmax))&&(signed_distance[2] < -(delta+rmax)))))) {

                  ++my_window_count[2];

                  // our turn to deal with this tri...
                  double this_bbminmax[4] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };
                  FOR_I3 {
                    const double * const x = surface->xsp[surface->spost[ist][i]];
                    const double dx[3] = DIFF(x,hcp_x0pdx0);
                    // x-primed is in the hcp_e1 direction...
                    const double xp = DOT_PRODUCT(dx,hcp_e1);
                    this_bbminmax[0] = min(this_bbminmax[0],xp); // xp-min
                    this_bbminmax[1] = min(this_bbminmax[1],-xp); // -xp-max
                    // and y-primed is in the hcp_e2 direction...
                    const double yp = DOT_PRODUCT(dx,hcp_e2);
                    this_bbminmax[2] = min(this_bbminmax[2],yp); // yp-min
                    this_bbminmax[3] = min(this_bbminmax[3],-yp); // -yp-max
                  }
                  FOR_I4 this_bbminmax[i] -= rmax;
                  // don't allow this to exceed the level 0 bbox...
                  FOR_I4 this_bbminmax[i] = max(this_bbminmax[i],bbminmax[0][i]);

                  // TODO also eliminate by considering the bbox in the 3rd dimension...

                  // now compute the j,k range of rays...
                  double xy[2];
                  xy[0] = this_bbminmax[0]; // xp-min
                  xy[1] = this_bbminmax[2]; // yp-min
                  int bbmin[2];
                  setJKFrom2D(bbmin,xy,hcp_dxp);
                  xy[0] = -this_bbminmax[1]; // xp-max
                  xy[1] = -this_bbminmax[3]; // yp-max
                  int bbmax[2];
                  setJKFrom2D(bbmax,xy,hcp_dxp);
                  // compute range in EXACTLY the way it is computed above...
                  int jmin = max( jminL, ( bbmin[0] < 0 ? bbmin[0]-(bbmin[0]%djk)     : bbmin[0]-(bbmin[0]%djk)+djk ));  assert(jmin > bbmin[0]);
                  int jmax = min( jmaxL, ( bbmax[0] < 0 ? bbmax[0]-(bbmax[0]%djk)-djk : bbmax[0]-(bbmax[0]%djk) ));      assert(jmax <= bbmax[0]);
                  int kmin = max( kminL, ( bbmin[1] < 0 ? bbmin[1]-(bbmin[1]%djk)     : bbmin[1]-(bbmin[1]%djk)+djk ));  assert(kmin > bbmin[1]);
                  int kmax = min( kmaxL, ( bbmax[1] < 0 ? bbmax[1]-(bbmax[1]%djk)-djk : bbmax[1]-(bbmax[1]%djk) ));      assert(kmax <= bbmax[1]);
                  // and cycle through all these rays...
                  for (int j = jmin; j <= jmax; j += djk) {
                    for (int k = kmin; k <= kmax; k += djk) {
                      if ((level == 0)||(j%(2*djk) != 0)||(k%(2*djk) != 0)) {
                        // j and k are even...
                        assert(j%2 == 0);
                        assert(k%2 == 0);
                        // only consider this ray if it has surface intersections, which we can check
                        // by seeing if the corresponding bit is set in the shm ray_buf...
                        const uint8 index = uint8((j-jminL)/djk)*uint8((kmaxL-kminL)/djk+1)+uint8((k-kminL)/djk);
                        //assert(index >= 0);

                        //TODO remove this eventually
                        {
                          const int j_check = index/((kmaxL-kminL)/djk+1)*djk+jminL;
                          const int k_check = (index-((j-jminL)/djk)*((kmaxL-kminL)/djk+1))*djk+kminL;
                          assert(j_check == j);
                          assert(k_check == k);
                        }

                        // must be a better way to do this, because tri is fixed...
                        double xray[3];
                        FOR_I3 xray[i] = hcp_x0pdx0[i] + (double(j)*hcp_e1[i] + double(k)*hcp_e2[i])*one_o_hcp_dxp;
                        // compute the tri corners in a reference frame where the ray is at (hcp_x0pdx0,0,0)...
                        double dxp[3][3];
                        FOR_I3 {
                          const double * const xsp = surface->xsp[surface->spost[ist][i]];
                          dxp[i][0] = (xsp[0]-xray[0])*hcp_e0[0] + (xsp[1]-xray[1])*hcp_e0[1] + (xsp[2]-xray[2])*hcp_e0[2]; // xray at hcp_x0pdx0 in e0 direction
                          dxp[i][1] = (xsp[0]-xray[0])*hcp_e1[0] + (xsp[1]-xray[1])*hcp_e1[1] + (xsp[2]-xray[2])*hcp_e1[2];
                          dxp[i][2] = (xsp[0]-xray[0])*hcp_e2[0] + (xsp[1]-xray[1])*hcp_e2[1] + (xsp[2]-xray[2])*hcp_e2[2];
                        }
                        // this function returns xp, the location along the ray's hcp_e0...
                        double xp;
                        double yzp[2];
                        const double d2 = getPointToTriDist2Special(xp,yzp,dxp);
                        // this tri may have several bloated tris, each associated with a level
                        // from st_level down to max(level,1) -- this special case when dealing
                        // with level 0 is because there is no refinement windows associated with
                        // level 0.
                        double r = rmax;
                        for (int il = max(level,1); il <= st_level; ++il) {
                          // break out as soon as this tri's d2 no longer intersects the ray
                          const double dr2 = r*r - d2;
                          if (dr2 < 4.0/(hcp_dxp*hcp_dxp))
                            break;
                          // place the in point at xp - sqrt(r*r > d2)...
                          double xp_minus[3] = { xp-sqrt(dr2), 0.0, 0.0 };
                          double xtri[3];
                          MiscUtils::getClosestPointOnTriRobust(xtri,xp_minus,dxp[0],dxp[1],dxp[2]);
                          int iter = 0;
                          while (fabs(DIST(xtri,xp_minus)-r) > 1.0E-6*2.0*one_o_hcp_dxp) {
                            xp_minus[0] = xtri[0]-sqrt(r*r-xtri[1]*xtri[1]-xtri[2]*xtri[2]);
                            MiscUtils::getClosestPointOnTriRobust(xtri,xp_minus,dxp[0],dxp[1],dxp[2]);
                            if (++iter >= 10000) {
                              //cout << " > warning, max iter reached: " << iter << " error: " << fabs(DIST(xtri,xp_minus)-r) << " tol: " << 1.0E-6*2.0*one_o_hcp_dxp << ". breaking... " << endl;
                              //cout.flush();
                              //cout << COUT_VEC(dxp[0]) << " " <<  COUT_VEC(dxp[1]) << " " <<  COUT_VEC(dxp[2]) << endl;
                              //cout.flush();
                              break;
                            }
                          }
                          double xp_plus[3] = { xp+sqrt(dr2), 0.0, 0.0 };
                          MiscUtils::getClosestPointOnTriRobust(xtri,xp_plus,dxp[0],dxp[1],dxp[2]);
                          iter = 0;
                          while (fabs(DIST(xtri,xp_plus)-r) > 1.0E-6*2.0*one_o_hcp_dxp) {
                            xp_plus[0] = xtri[0]+sqrt(r*r-xtri[1]*xtri[1]-xtri[2]*xtri[2]);
                            MiscUtils::getClosestPointOnTriRobust(xtri,xp_plus,dxp[0],dxp[1],dxp[2]);
                            if (++iter >= 10000) {
                              //cout << " > warning, max iter reached: " << iter << " error: " << fabs(DIST(xtri,xp_plus)-r) << " tol: " << 1.0E-6*2.0*one_o_hcp_dxp << ". breaking...  " << endl;
                              //cout.flush();
                              //cout << COUT_VEC(dxp[0]) << " " <<  COUT_VEC(dxp[1]) << " " <<  COUT_VEC(dxp[2]) << endl;
                              //cout.flush();
                              break;
                            }
                          }
                          assert((xp_plus[0] == xp_plus[0])&&(xp_minus[0] == xp_minus[0]));
                          //double xmp[2];
                          //GeomUtils::calcBloatedTriIntersections(xmp,yzp,dxp[0],dxp[1],dxp[2],r);
                          //double xp_minus[3] = { xmp[0], 0.0, 0.0 };
                          //double xp_plus[3] = { xmp[1], 0.0, 0.0 };
                          // this is an in-out point for level "il" on a level "level" ray...
                          if (windowVec.size() >= threshold_windowVec_size) {
                            // here we assume collective synchronization amongst ranks...
#include "compressAndSend.hpp"
                            assert(windowVec.empty());
                            ++exchange_count;
                            // because some ranks may be done their rays before others, this logic
                            // keeps the calls to compressAndSend.hpp synchronized...
                            // my_done should be either 0 or 1...
                            int my_done = 0;
                            // any "0" has to return "0"...
                            int done;
                            MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);
                            assert(done == 0);
                          }
                          assert(il >= 1);
                          assert(iwindow >= 0);
                          const int flag = ((iwindow+1)<<HCP_WINDOW_SHIFT_BITS)|il;
                          windowVec.push_back(InOutData(xp_minus[0],index,flag+1)); // use level+1 convention here: same as inoutVec
                          windowVec.push_back(InOutData(xp_plus[0],index,-flag-1));
                          my_window_count[0] += 2;
                          // and modify r for next level...
                          const int st_nlayers = nlowi[iwindow*(level_max+1)+il];
                          //r -= double(st_nlayers)*hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<il);
                          r -= double(st_nlayers)*hcp_delta/double(1<<il);
                        }
                      }
                    }
                  }
                }
              }
              ++round_robin_rank;
              if (round_robin_rank == mpi_size)
                round_robin_rank = 0;
            }
          }
        } // surface

        // --------------------------------------------------------------
        // ff_surface windows...
        // --------------------------------------------------------------
        if (PartData::partVec[ipart]->ff_surface) {
          SurfaceShm * ff_surface = PartData::partVec[ipart]->ff_surface;
          for (int ist = 0; ist < ff_surface->nst; ++ist) {
            int st_level = 0;
            while (dxost_factor*PartData::partVec[ipart]->ff_surface_dxost[ist] < hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<st_level))
              ++st_level;
            assert(st_level <= level_max);
            if (st_level >= max(level,1)) {
              // this tri can affect this level...
              if (round_robin_rank == mpi_rank) {

                ++my_window_count[1];

                // calculate the maximum bloat...
                double rmax = 0.0;
                for (int il = max(level,1); il <= st_level; ++il) {
                  const int st_nlayers = PartData::partVec[ipart]->ff_nlayers[il];
                  rmax += double(st_nlayers)*hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<il);
                }
                const double delta = hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<level);

                // -------------------------------------------------------------------------
                // introduce periodic looping through all potential periodicities here...
                // -------------------------------------------------------------------------
                for (int iper = 0; iper < periodicRtVec.size(); ++iper) {

                  double xsp_t[3][3];
                  FOR_I3 FOR_J3 xsp_t[i][j] = ff_surface->xsp[ff_surface->spost[ist][i]][j];
                  if (periodicRtVec[iper].b_t) {
                    //cout << " periodicRtVec[iper].t: " << COUT_VEC(periodicRtVec[iper].t) << endl;
                    FOR_I3 FOR_J3 xsp_t[i][j] += periodicRtVec[iper].t[j];
                  }
                  if (periodicRtVec[iper].b_R) {
                    //assert(0);
                    double tmp[3];
                    FOR_I3 {
                      FOR_J3 tmp[j] = xsp_t[i][j];
                      MiscUtils::applyRotation(xsp_t[i],periodicRtVec[iper].R,tmp);
                    }
                  }

                  // if the tri+bloat is outside the thickened plane (if it exists), we can skip it...
                  double signed_distance[3];
                  if (plane_xp != NULL) {
                    FOR_I3 {
                      const double dx[3] = DIFF(xsp_t[i],plane_xp);
                      signed_distance[i] = DOT_PRODUCT(dx,unit_np);
                    }
                  }

                  if ((plane_xp == NULL)||
                      (!(((signed_distance[0] > (delta+rmax))&&(signed_distance[1] > (delta+rmax))&&(signed_distance[2] > (delta+rmax)))||
                         ((signed_distance[0] < -(delta+rmax))&&(signed_distance[1] < -(delta+rmax))&&(signed_distance[2] < -(delta+rmax)))))) {

                    ++my_window_count[2];

                    // our turn to deal with this tri...
                    double this_bbminmax[4] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };
                    FOR_I3 {
                      const double dx[3] = DIFF(xsp_t[i],hcp_x0pdx0);
                      // x-primed is in the hcp_e1 direction...
                      const double xp = DOT_PRODUCT(dx,hcp_e1);
                      this_bbminmax[0] = min(this_bbminmax[0],xp); // xp-min
                      this_bbminmax[1] = min(this_bbminmax[1],-xp); // -xp-max
                      // and y-primed is in the hcp_e2 direction...
                      const double yp = DOT_PRODUCT(dx,hcp_e2);
                      this_bbminmax[2] = min(this_bbminmax[2],yp); // yp-min
                      this_bbminmax[3] = min(this_bbminmax[3],-yp); // -yp-max
                    }
                    FOR_I4 this_bbminmax[i] -= rmax;
                    // don't allow this to exceed the level 0 bbox...
                    FOR_I4 this_bbminmax[i] = max(this_bbminmax[i],bbminmax[0][i]);

                    // TODO also eliminate by considering the bbox in the 3rd dimension...

                    // now compute the j,k range of rays...
                    double xy[2];
                    xy[0] = this_bbminmax[0]; // xp-min
                    xy[1] = this_bbminmax[2]; // yp-min
                    int bbmin[2];
                    setJKFrom2D(bbmin,xy,hcp_dxp);
                    xy[0] = -this_bbminmax[1]; // xp-max
                    xy[1] = -this_bbminmax[3]; // yp-max
                    int bbmax[2];
                    setJKFrom2D(bbmax,xy,hcp_dxp);
                    // compute range in EXACTLY the way it is computed above...
                    int jmin = max( jminL, ( bbmin[0] < 0 ? bbmin[0]-(bbmin[0]%djk)     : bbmin[0]-(bbmin[0]%djk)+djk ));  assert(jmin > bbmin[0]);
                    int jmax = min( jmaxL, ( bbmax[0] < 0 ? bbmax[0]-(bbmax[0]%djk)-djk : bbmax[0]-(bbmax[0]%djk) ));      assert(jmax <= bbmax[0]);
                    int kmin = max( kminL, ( bbmin[1] < 0 ? bbmin[1]-(bbmin[1]%djk)     : bbmin[1]-(bbmin[1]%djk)+djk ));  assert(kmin > bbmin[1]);
                    int kmax = min( kmaxL, ( bbmax[1] < 0 ? bbmax[1]-(bbmax[1]%djk)-djk : bbmax[1]-(bbmax[1]%djk) ));      assert(kmax <= bbmax[1]);
                    // and cycle through all these rays...
                    for (int j = jmin; j <= jmax; j += djk) {
                      for (int k = kmin; k <= kmax; k += djk) {
                        if ((level == 0)||(j%(2*djk) != 0)||(k%(2*djk) != 0)) {
                          // j and k are even...
                          assert(j%2 == 0);
                          assert(k%2 == 0);
                          // only consider this ray if it has surface intersections, which we can check
                          // by seeing if the corresponding bit is set in the shm ray_buf...
                          const uint8 index = uint8((j-jminL)/djk)*uint8((kmaxL-kminL)/djk+1)+uint8((k-kminL)/djk);
                          //assert(index >= 0);

                          //TODO remove this eventually
                          {
                            const int j_check = index/((kmaxL-kminL)/djk+1)*djk+jminL;
                            const int k_check = (index-((j-jminL)/djk)*((kmaxL-kminL)/djk+1))*djk+kminL;
                            assert(j_check == j);
                            assert(k_check == k);
                          }

                          // must be a better way to do this, because tri is fixed...
                          double xray[3];
                          FOR_I3 xray[i] = hcp_x0pdx0[i] + (double(j)*hcp_e1[i] + double(k)*hcp_e2[i])*one_o_hcp_dxp;
                          // compute the tri corners in a reference frame where the ray is at (hcp_x0pdx0,0,0)...
                          double dxp[3][3];
                          FOR_I3 {
                            dxp[i][0] = (xsp_t[i][0]-xray[0])*hcp_e0[0] + (xsp_t[i][1]-xray[1])*hcp_e0[1] + (xsp_t[i][2]-xray[2])*hcp_e0[2]; // xray at hcp_x0pdx0 in e0 direction
                            dxp[i][1] = (xsp_t[i][0]-xray[0])*hcp_e1[0] + (xsp_t[i][1]-xray[1])*hcp_e1[1] + (xsp_t[i][2]-xray[2])*hcp_e1[2];
                            dxp[i][2] = (xsp_t[i][0]-xray[0])*hcp_e2[0] + (xsp_t[i][1]-xray[1])*hcp_e2[1] + (xsp_t[i][2]-xray[2])*hcp_e2[2];
                          }
                          // this function returns xp, the location along the ray's hcp_e0...
                          double xp;
                          double yzp[2];
                          const double d2 = getPointToTriDist2Special(xp,yzp,dxp);
                          // this tri may have several bloated tris, each associated with a level
                          // from st_level down to max(level,1) -- this special case when dealing
                          // with level 0 is because there is no refinement windows associated with
                          // level 0.
                          double r = rmax;
                          for (int il = max(level,1); il <= st_level; ++il) {
                            // break out as soon as this tri's d2 no longer intersects the ray
                            const double dr2 = r*r - d2;
                            if (dr2 < 4.0/(hcp_dxp*hcp_dxp))
                              break;
                            // place the in point at xp - sqrt(r*r > d2)...
                            double xp_minus[3] = { xp-sqrt(dr2), 0.0, 0.0 };
                            double xtri[3];
                            MiscUtils::getClosestPointOnTriRobust(xtri,xp_minus,dxp[0],dxp[1],dxp[2]);
                            int iter = 0;
                            while (fabs(DIST(xtri,xp_minus)-r) > 1.0E-6*2.0*one_o_hcp_dxp) {
                              xp_minus[0] = xtri[0]-sqrt(r*r-xtri[1]*xtri[1]-xtri[2]*xtri[2]);
                              MiscUtils::getClosestPointOnTriRobust(xtri,xp_minus,dxp[0],dxp[1],dxp[2]);
                              if (++iter >= 10000) break;
                            }
                            double xp_plus[3] = { xp+sqrt(dr2), 0.0, 0.0 };
                            MiscUtils::getClosestPointOnTriRobust(xtri,xp_plus,dxp[0],dxp[1],dxp[2]);
                            iter = 0;
                            while (fabs(DIST(xtri,xp_plus)-r) > 1.0E-6*2.0*one_o_hcp_dxp) {
                              xp_plus[0] = xtri[0]+sqrt(r*r-xtri[1]*xtri[1]-xtri[2]*xtri[2]);
                              MiscUtils::getClosestPointOnTriRobust(xtri,xp_plus,dxp[0],dxp[1],dxp[2]);
                              if (++iter >= 10000) break;
                            }
                            //double xmp[2];
                            //GeomUtils::calcBloatedTriIntersections(xmp,yzp,dxp[0],dxp[1],dxp[2],r);
                            //double xp_minus[3] = { xmp[0], 0.0, 0.0 };
                            //double xp_plus[3] = { xmp[1], 0.0, 0.0 };
                            // this is an in-out point for level "il" on a level "level" ray...
                            if (windowVec.size() >= threshold_windowVec_size) {
                              // here we assume collective synchronization amongst ranks...
#include "compressAndSend.hpp"
                              assert(windowVec.empty());
                              ++exchange_count;
                              // because some ranks may be done their rays before others, this logic
                              // keeps the calls to compressAndSend.hpp synchronized...
                              // my_done should be either 0 or 1...
                              int my_done = 0;
                              // any "0" has to return "0"...
                              int done;
                              MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);
                              assert(done == 0);
                            }
                            assert(il >= 1);
                            const int flag = ((hcpWindowVec.size()+ipart+1)<<HCP_WINDOW_SHIFT_BITS)|il; // put ff windows at the end: TODO fix later
                            windowVec.push_back(InOutData(xp_minus[0],index,flag+1)); // use level+1 convention here: same as inoutVec
                            windowVec.push_back(InOutData(xp_plus[0],index,-flag-1));
                            my_window_count[0] += 2;
                            // and modify r for next level...
                            const int st_nlayers = PartData::partVec[ipart]->ff_nlayers[il];
                            r -= double(st_nlayers)*hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<il);
                          }
                        }
                      }
                    }
                  }
                }
              }
              ++round_robin_rank;
              if (round_robin_rank == mpi_size)
                round_robin_rank = 0;
            }
          }
        } // ff_surface

      }

      // if this rank made it out here, we exchange until finished...
      int done = 0;
      while (done != 1) {
#include "compressAndSend.hpp"
	assert(windowVec.empty());
        ++exchange_count;
        // my_done should be either 0 or 1...
        int my_done = 1;
        MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);
      }

      // do mixed tris for each window...
      for (int igw = 0; igw < ngw; ++igw) {

        // ======================================================================================
        // which rays at this level pass through the window geometry...
        // ======================================================================================

        mcount = 0;
        pcount = 0;
        assert(indexVec.empty());
        assert(xpVec.empty());
        FOR_RANK send_count[rank] = 0;

        for (int ist = tvogw_i[igw]+mpi_rank; ist < tvogw_i[igw+1]; ist += mpi_size) {
          // we need an integer version of this tri...
          const double * const xsp0 = triVec[ist].x[0];
          const double * const xsp1 = triVec[ist].x[1];
          const double * const xsp2 = triVec[ist].x[2];
          int jk0[2]; setJKFrom3D(jk0,xsp0,hcp_dxp);
          int jk1[2]; setJKFrom3D(jk1,xsp1,hcp_dxp);
          int jk2[2]; setJKFrom3D(jk2,xsp2,hcp_dxp);
          int bbmin[2]; FOR_I2 bbmin[i] = min(jk0[i],min(jk1[i],jk2[i]));
          int bbmax[2]; FOR_I2 bbmax[i] = max(jk0[i],max(jk1[i],jk2[i]));
          int jmin = max( jminL, ( bbmin[0] < 0 ? bbmin[0]-(bbmin[0]%djk)     : bbmin[0]-(bbmin[0]%djk)+djk ));  assert(jmin > bbmin[0]);
          int jmax = min( jmaxL, ( bbmax[0] < 0 ? bbmax[0]-(bbmax[0]%djk)-djk : bbmax[0]-(bbmax[0]%djk) ));      assert(jmax <= bbmax[0]);
          int kmin = max( kminL, ( bbmin[1] < 0 ? bbmin[1]-(bbmin[1]%djk)     : bbmin[1]-(bbmin[1]%djk)+djk ));  assert(kmin > bbmin[1]);
          int kmax = min( kmaxL, ( bbmax[1] < 0 ? bbmax[1]-(bbmax[1]%djk)-djk : bbmax[1]-(bbmax[1]%djk) ));      assert(kmax <= bbmax[1]);
          for (int j = jmin; j <= jmax; j += djk) {
            for (int k = kmin; k <= kmax; k += djk) {
              if ((level == 0)||(j%(2*djk) != 0)||(k%(2*djk) != 0)) {
                // j and k are even...
                assert(j%2 == 0);
                assert(k%2 == 0);
                const int8 dx0[2] = { int8( jk0[0] - j ), int8( jk0[1] - k ) };
                const int8 dx1[2] = { int8( jk1[0] - j ), int8( jk1[1] - k ) };
                const int8 dx2[2] = { int8( jk2[0] - j ), int8( jk2[1] - k ) };
                const int8 V0 = dx1[0]*dx2[1] - dx1[1]*dx2[0];
                const int8 V1 = dx2[0]*dx0[1] - dx2[1]*dx0[0];
                const int8 V2 = dx0[0]*dx1[1] - dx0[1]*dx1[0];
                if ((V0 != 0)||(V1 != 0)||(V2 != 0)) {
                  if ((V0 >= 0)&&(V1 >= 0)&&(V2 >= 0)) {
                    // build the intersection "x" along the e0 direction...
                    const double xp = ( double(V0)*((xsp0[0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp0[1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp0[2]-hcp_x0pdx0[2])*hcp_e0[2]) +
                                        double(V1)*((xsp1[0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp1[1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp1[2]-hcp_x0pdx0[2])*hcp_e0[2]) +
                                        double(V2)*((xsp2[0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp2[1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp2[2]-hcp_x0pdx0[2])*hcp_e0[2]) )/double(V0+V1+V2);
                    xpVec.push_back(xp);
                    const uint8 index = uint8((j-jminL)/djk)*uint8((kmaxL-kminL)/djk+1)+uint8((k-kminL)/djk);
                    //assert(index >= 0);
                    //assert(index < TWO_BILLION/2); // need space for 2 bits

                    // TODO check: remove this eventually...
                    {
                      const int j_check = index/((kmaxL-kminL)/djk+1)*djk+jminL;
                      const int k_check = (index-((j-jminL)/djk)*((kmaxL-kminL)/djk+1))*djk+kminL;
                      assert(j_check == j);
                      assert(k_check == k);
                    }

                    const int rank = index%uint8(mpi_size);
                    assert((rank >=0) && (rank < mpi_size));
                    ++send_count[rank];
                    if ((V0 == 0)||(V1 == 0)||(V2 == 0)) {
                      // this is an edge intersection. Because of the staggered spaces, it cannot
                      // be a vertex intersection. Here we add 1 because there should be another.  Note
                      // that we can assert that only one of the volumes is zero -- i.e. we are on an edge...
                      assert( ((V0 > 0)&&(V1 > 0)) || ((V0 > 0)&&(V2 > 0)) || ((V1 > 0)&&(V2 > 0)) );
                      pcount += 1;
                      indexVec.push_back(index<<2);     // i.e. 00
                      assert((indexVec.back()>>2) == index);
                    }
                    else {
                      assert((V0 > 0)&&(V1 > 0)&&(V2 > 0));
                      // all positive -- intersection is good and gets a 2...
                      pcount += 2;
                      indexVec.push_back((index<<2)|1ll); // i.e. 01
                      assert((indexVec.back()>>2) == index);
                    }
                  }
                  else if ((V0 <= 0)&&(V1 <= 0)&&(V2 <= 0)) {
                    const double xp = ( double(V0)*((xsp0[0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp0[1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp0[2]-hcp_x0pdx0[2])*hcp_e0[2]) +
                                        double(V1)*((xsp1[0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp1[1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp1[2]-hcp_x0pdx0[2])*hcp_e0[2]) +
                                        double(V2)*((xsp2[0]-hcp_x0pdx0[0])*hcp_e0[0]+(xsp2[1]-hcp_x0pdx0[1])*hcp_e0[1]+(xsp2[2]-hcp_x0pdx0[2])*hcp_e0[2]) )/double(V0+V1+V2);
                    xpVec.push_back(xp);
                    const uint8 index = uint8((j-jminL)/djk)*uint8((kmaxL-kminL)/djk+1)+uint8((k-kminL)/djk);
                    //assert(index >= 0);
                    //assert(index < TWO_BILLION/2); // need space for 2 bits

                    // TODO check: remove this eventually...
                    {
                      const int j_check = index/((kmaxL-kminL)/djk+1)*djk+jminL;
                      const int k_check = (index-((j-jminL)/djk)*((kmaxL-kminL)/djk+1))*djk+kminL;
                      assert(j_check == j);
                      assert(k_check == k);
                    }

                    const int rank = index%uint8(mpi_size);
                    assert((rank >=0) && (rank < mpi_size));
                    ++send_count[rank];
                    if ((V0 == 0)||(V1 == 0)||(V2 == 0)) {
                      // see note above...
                      assert( ((V0 < 0)&&(V1 < 0)) || ((V0 < 0)&&(V2 < 0)) || ((V1 < 0)&&(V2 < 0)) );
                      mcount += 1;
                      indexVec.push_back((index<<2)|2ll); // i.e. 10
                      assert((indexVec.back()>>2) == index);
                    }
                    else {
                      assert((V0 < 0)&&(V1 < 0)&&(V2 < 0));
                      // all negative...
                      mcount += 2;
                      indexVec.push_back((index<<2)|3ll); // i.e. 11
                      assert((indexVec.back()>>2) == index);
                    }
                  }
                }
              }
            }
          }
        }
        assert(indexVec.size() == xpVec.size());

        send_disp[0] = 0;
        int my_sum = send_count[0];
        for (int rank = 1; rank < mpi_size; ++rank) {
          my_sum += send_count[rank];
          send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
        }
        const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
        assert(send_count_sum == indexVec.size());

        assert(send_buf_uint8 == NULL); send_buf_uint8 = new uint8[send_count_sum];
        assert(send_buf_double == NULL); send_buf_double = new double[send_count_sum];
        for (int isend = 0; isend < send_count_sum; ++isend) {
          // to get the index back, shift down 2 bits...
          const uint8 index = (indexVec[isend]>>2);
          const int rank = index%uint8(mpi_size);
          assert((rank >= 0)&&(rank < mpi_size));
          send_buf_uint8[send_disp[rank]] = indexVec[isend];
          send_buf_double[send_disp[rank]] = xpVec[isend];
          ++send_disp[rank];
        }
        indexVec.clear();
        xpVec.clear();

        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

        // set up recv-side stuff...

        MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

        recv_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
        const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

        assert(recv_count_sum >= 0);
        assert(recv_buf_uint8 == NULL); recv_buf_uint8 = new uint8[recv_count_sum];
        //for (int irecv = 0; irecv < recv_count_sum; ++irecv)
        //  recv_buf_uint8[irecv] = ULONG_MAX;
        MPI_Alltoallv(send_buf_uint8,send_count,send_disp,MPI_UINT8,
                      recv_buf_uint8,recv_count,recv_disp,MPI_UINT8,mpi_comm);
        delete[] send_buf_uint8; send_buf_uint8 = NULL;
        //for (int irecv = 0; irecv < recv_count_sum; ++irecv)
        //  assert(recv_buf_uint8[irecv] < ULONG_MAX);

        assert(recv_buf_double == NULL); recv_buf_double = new double[recv_count_sum];
        //for (int irecv = 0; irecv < recv_count_sum; ++irecv)
        //  recv_buf_double[irecv] = HUGE_VAL;
        MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
                      recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
        delete[] send_buf_double; send_buf_double = NULL;
        //for (int irecv = 0; irecv < recv_count_sum; ++irecv)
        //  assert(recv_buf_double[irecv] < HUGE_VAL);

        int my_buf[2] = { mcount, pcount};
        int buf[2];
        MPI_Reduce(my_buf,buf,2,MPI_INT,MPI_SUM,0,mpi_comm);

        if (mpi_rank == 0) {
          cout << " > level " << level << " non-fazone window " << igw << endl;
          //if (buf[0] != buf[1]) 
          //  cout << "XXXXX: " << buf[0] << " " << buf[1] << endl;
          assert(buf[0] == buf[1]); // mcount == pcount
          double wtime = MPI_Wtime();
          cout << " > window ray intersection time: " << wtime-wtime_prev << " [s], total time so far: " << wtime-wtime0 << " [s]" << endl;
          wtime_prev = wtime;
        }

        for (int ii = 0; ii < my_ray_flag_size; ++ii)
          my_ray_flag[ii] = 0;

        for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
          const uint8 index = (recv_buf_uint8[irecv]>>2);
          assert((index-uint8(mpi_rank))%uint8(mpi_size) == 0ll);
          const int ii = (index-uint8(mpi_rank))/uint8(mpi_size);
          assert((ii >= 0)&&(ii < my_ray_flag_size));
          ++my_ray_flag[ii];
        }

        int count = 0;
        for (int ii = 0; ii < my_ray_flag_size; ++ii) {
          const int tmp = my_ray_flag[ii];
          my_ray_flag[ii] = count;
          count += tmp;
        }
        assert(count == recv_count_sum);

        assert(windowVec.empty());
        windowVec.resize(count);

        for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
          const uint8 index = (recv_buf_uint8[irecv]>>2);
          assert((index-uint8(mpi_rank))%uint8(mpi_size) == 0ll);
          const int ii = (index-uint8(mpi_rank))/uint8(mpi_size);
          assert((ii >= 0)&&(ii < my_ray_flag_size));
          const int iv = my_ray_flag[ii]++;
          windowVec[iv].xp = recv_buf_double[irecv];
          windowVec[iv].index = index;
          const int bits = (recv_buf_uint8[irecv]&3ll);
          if (bits == 0) {
            // one "out" (associated with pcount) -- must have been an edge intersection...
            windowVec[iv].level = -1;
          }
          else if (bits == 1) {
            // two "out" -- regular triangle intersection...
            windowVec[iv].level = -2;
          }
          else if (bits == 2) {
            // one "in" (associated with mcount) -- edge...
            windowVec[iv].level = 1;
          }
          else {
            assert(bits == 3);
            // two "in"
            windowVec[iv].level = 2;
          }
        }

        delete[] recv_buf_uint8; recv_buf_uint8 = NULL;
        delete[] recv_buf_double; recv_buf_double = NULL;

        // and condense the windowVec: this is a special condense
        // the first time because of the +/- 1,2 stuff associated with these
        // exact intersections. Once condensed, leave the in/out points as simply
        // +1+il/-1-il... i.e. use a level+1 strategy

        vector<InOutData>::iterator iter_end = windowVec.begin();
        while (iter_end != windowVec.end()) {
          vector<InOutData>::iterator iter_begin = iter_end++;
          while ((iter_end != windowVec.end())&&(iter_end->index == iter_begin->index))
            ++iter_end;
          sort(iter_begin,iter_end);
          vector<InOutData>::iterator iter_next = iter_begin;
          vector<InOutData>::iterator iter_prev = iter_end;
          while (iter_next != iter_end) {
            vector<InOutData>::iterator iter = iter_next; iter_next++;
            if ((iter_next != iter_end)&&(iter->level*iter_next->level == -1)) {
              //if ((iter_next != iter_end)&&(iter->level*iter_next->level == -1)&&(fabs(iter->xp-iter_next->xp) < 2.0*one_o_hcp_dxp)) {
              // add check to make sure we are not canceling a 1 -1 when the 1 was needed to complete a "face" intersection (1 1 -> 2).
              // this occured when the ff was on top of the surface...
              if ((iter_prev == iter_end)||(iter_prev->level*iter->level != 1)) {
                iter->level = 0;
                iter_next->level = 0;
                iter_next++; // skip to next iter,iter_next pair
              }
            }
            iter_prev = iter;
          }
          int level_count0 = 0;
          bool on_flag = false;
          for (vector<InOutData>::iterator iter = iter_begin; iter != iter_end; ++iter) {
            if (iter->level == 0) continue;
            level_count0 += iter->level;
            assert(triVec[tvogw_i[igw]].window >= 0);
            const int st_level = (leowi != NULL) ? leowi[triVec[tvogw_i[igw]].window]:0;
            const int flag = ((triVec[tvogw_i[igw]].window+1)<<HCP_WINDOW_SHIFT_BITS)|st_level;
            if (level_count0 == 2) {
              // keep this one. assume it is the "in" point...
              assert(on_flag == false);
              iter->level = 1+flag; // force 1+il
              on_flag = true;
            }
            else if (level_count0 == 0) {
              // if we were positive (i.e. on) then this shuts off the ray...
              if (on_flag) {
                iter->level = -1-flag; // force -1-il
                on_flag = false;
              }
              else {
                // we returned to level_count0 == 0 from below, so do not turn on...
                iter->level = 0;
              }
            }
            else {
              // anything else we can discard...
              iter->level = 0;
            }
          }
          if (level_count0 != 0) {
            cout << "About to fail in/out balance in window intersection." << endl;
            for (vector<InOutData>::iterator iter = iter_begin; iter != iter_end; ++iter)
              cout << iter->level << " " << endl;
            cout.flush();
          }
          assert(level_count0 == 0);
          // just make sure we are off...
          assert(on_flag == false);
        }

        // now all the zeros should be removed...

        int new_size = 0;
        for (int ii = 0, limit = windowVec.size(); ii < limit; ++ii) {
          if (windowVec[ii].level != 0) {
            windowVec[new_size] = windowVec[ii];
            ++new_size;
          }
        }
        windowVec.resize(new_size);

        for (int ii = 0; ii < my_ray_flag_size; ++ii)
          my_ray_flag[ii] = 0;
        int jj_prev = -1;
        for (int ii = 0, limit = inoutVec[level].size(); ii < limit; ++ii) {
          const uint8 index = inoutVec[level][ii].index;
          assert((index-uint8(mpi_rank))%uint8(mpi_size) == 0ll);
          const int jj = (index-uint8(mpi_rank))/uint8(mpi_size);
          assert((jj >= 0)&&(jj < my_ray_flag_size));
          assert(jj >= jj_prev); // check sort
          jj_prev = jj;
          ++my_ray_flag[jj];
        }
        bool* req_info = new bool[new_size];
        for (int iv = 0; iv < new_size; ++iv) {
          const uint8 index = windowVec[iv].index;
          assert((index-uint8(mpi_rank))%uint8(mpi_size) == 0ll);
          const int ii = (index-uint8(mpi_rank))/uint8(mpi_size);
          assert((ii >= 0)&&(ii < my_ray_flag_size));
          if (my_ray_flag[ii] > 0) {
            // count only those rays that have existing surface intersections
            // and/or other entries...
            ++my_ray_flag[ii];
            req_info[iv] = true;
          }
          else {
            // this is refinement window information that is not
            // required because the ray does not pass through the
            // geometry.
            req_info[iv] = false;
          }
        }

        // my_ray_flag now has the counts associated with each level
        // in it. Turn it into the jj offset in each level...
        for (int ii = 1; ii < my_ray_flag_size; ++ii)
          my_ray_flag[ii] += my_ray_flag[ii-1];
        // now go backwards through the inoutVec[level], making room
        // for the new info...
        const int8 old_size = inoutVec[level].size();
        inoutVec[level].resize(my_ray_flag[my_ray_flag_size-1]);
        for (int ii = old_size-1; ii >= 0; --ii) {
          const uint8 index = inoutVec[level][ii].index;
          const int jj = (index-uint8(mpi_rank))/uint8(mpi_size);
          --my_ray_flag[jj];
          const int ii_new = my_ray_flag[jj];
          assert(ii_new >= ii);
          inoutVec[level][ii_new] = inoutVec[level][ii];
        }
        for (int iv = 0; iv < new_size; ++iv) {
          if (req_info[iv]) {
            const uint8 index = windowVec[iv].index;
            assert((index-uint8(mpi_rank))%uint8(mpi_size) == 0ll);
            const int jj = (index-uint8(mpi_rank))/uint8(mpi_size);
            --my_ray_flag[jj];
            const int ii_new = my_ray_flag[jj];
            inoutVec[level][ii_new].index = windowVec[iv].index;
            inoutVec[level][ii_new].level = windowVec[iv].level;
            inoutVec[level][ii_new].xp    = windowVec[iv].xp;
          }
        }
        delete[] req_info;
        windowVec.clear();

      } // igw loop

      // do bloating for all windows at once...

      round_robin_rank = 0;
      for (int ist = 0, limit = triVec.size(); ist < limit; ++ist) {
        assert(triVec[ist].window >= 0);
	const int st_level = (leowi != NULL) ? leowi[triVec[ist].window]-1:0;
	if (st_level >= max(level,1)) {
	  // this tri can affect this level...
	  if (round_robin_rank == mpi_rank) {

	    ++my_window_count[1];

            // calculate the maximum bloat...
            assert(triVec[ist].window >= 0);
            //const int st_nlayers = nlowi[triVec[ist].window];
	    double rmax = 0.0;
	    for (int il = max(level,1); il <= st_level; ++il) {
              const int st_nlayers = nlowi[triVec[ist].window*(level_max+1)+il];
              //rmax += double(st_nlayers)*hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<il);
              rmax += double(st_nlayers)*hcp_delta/double(1<<il);
	    }
            const double delta = hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<level);

            // if the tri+bloat is outside the thickened plane (if it exists), we can skip it...
            double signed_distance[3];
            if (plane_xp != NULL) {
              FOR_I3 {
                const double * const x = triVec[ist].x[i];
                const double dx[3] = DIFF(x,plane_xp);
                signed_distance[i] = DOT_PRODUCT(dx,unit_np);
              }
            }

            if ((plane_xp == NULL)||
                (!(((signed_distance[0] > (delta+rmax))&&(signed_distance[1] > (delta+rmax))&&(signed_distance[2] > (delta+rmax)))||
                   ((signed_distance[0] < -(delta+rmax))&&(signed_distance[1] < -(delta+rmax))&&(signed_distance[2] < -(delta+rmax)))))) {

	      ++my_window_count[2];

              // our turn to deal with this tri...
              double this_bbminmax[4] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };
              FOR_I3 {
                const double * const x = triVec[ist].x[i];
                const double dx[3] = DIFF(x,hcp_x0pdx0);
                // x-primed is in the hcp_e1 direction...
                const double xp = DOT_PRODUCT(dx,hcp_e1);
                this_bbminmax[0] = min(this_bbminmax[0],xp); // xp-min
                this_bbminmax[1] = min(this_bbminmax[1],-xp); // -xp-max
                // and y-primed is in the hcp_e2 direction...
                const double yp = DOT_PRODUCT(dx,hcp_e2);
                this_bbminmax[2] = min(this_bbminmax[2],yp); // yp-min
                this_bbminmax[3] = min(this_bbminmax[3],-yp); // -yp-max
              }
              FOR_I4 this_bbminmax[i] -= rmax;
              // don't allow this to exceed the level 0 bbox...
              FOR_I4 this_bbminmax[i] = max(this_bbminmax[i],bbminmax[0][i]);

              // TODO also eliminate by considering the bbox in the 3rd dimension...

              // now compute the j,k range of rays...
              double xy[2];
              xy[0] = this_bbminmax[0]; // xp-min
              xy[1] = this_bbminmax[2]; // yp-min
              int bbmin[2];
              setJKFrom2D(bbmin,xy,hcp_dxp);
              xy[0] = -this_bbminmax[1]; // xp-max
              xy[1] = -this_bbminmax[3]; // yp-max
              int bbmax[2];
              setJKFrom2D(bbmax,xy,hcp_dxp);
              // compute range in EXACTLY the way it is computed above...
              int jmin = max( jminL, ( bbmin[0] < 0 ? bbmin[0]-(bbmin[0]%djk)     : bbmin[0]-(bbmin[0]%djk)+djk ));  assert(jmin > bbmin[0]);
              int jmax = min( jmaxL, ( bbmax[0] < 0 ? bbmax[0]-(bbmax[0]%djk)-djk : bbmax[0]-(bbmax[0]%djk) ));      assert(jmax <= bbmax[0]);
              int kmin = max( kminL, ( bbmin[1] < 0 ? bbmin[1]-(bbmin[1]%djk)     : bbmin[1]-(bbmin[1]%djk)+djk ));  assert(kmin > bbmin[1]);
              int kmax = min( kmaxL, ( bbmax[1] < 0 ? bbmax[1]-(bbmax[1]%djk)-djk : bbmax[1]-(bbmax[1]%djk) ));      assert(kmax <= bbmax[1]);
              // and cycle through all these rays...
              for (int j = jmin; j <= jmax; j += djk) {
                for (int k = kmin; k <= kmax; k += djk) {
                  if ((level == 0)||(j%(2*djk) != 0)||(k%(2*djk) != 0)) {
                    // j and k are even...
                    assert(j%2 == 0);
                    assert(k%2 == 0);
                    // only consider this ray if it has surface intersections, which we can check
                    // by seeing if the corresponding bit is set in the shm ray_buf...
                    const uint8 index = uint8((j-jminL)/djk)*uint8((kmaxL-kminL)/djk+1)+uint8((k-kminL)/djk);
                    //assert(index >= 0);

		    //TODO remove this eventually
		    {
		      const int j_check = index/((kmaxL-kminL)/djk+1)*djk+jminL;
		      const int k_check = (index-((j-jminL)/djk)*((kmaxL-kminL)/djk+1))*djk+kminL;
		      assert(j_check == j);
		      assert(k_check == k);
		    }

                    // must be a better way to do this, because tri is fixed...
                    double xray[3];
                    FOR_I3 xray[i] = hcp_x0pdx0[i] + (double(j)*hcp_e1[i] + double(k)*hcp_e2[i])*one_o_hcp_dxp;
                    // compute the tri corners in a reference frame where the ray is at (hcp_x0pdx0,0,0)...
                    double dxp[3][3];
                    FOR_I3 {
                      const double * const xsp = triVec[ist].x[i];
                      dxp[i][0] = (xsp[0]-xray[0])*hcp_e0[0] + (xsp[1]-xray[1])*hcp_e0[1] + (xsp[2]-xray[2])*hcp_e0[2]; // xray at hcp_x0pdx0 in e0 direction
                      dxp[i][1] = (xsp[0]-xray[0])*hcp_e1[0] + (xsp[1]-xray[1])*hcp_e1[1] + (xsp[2]-xray[2])*hcp_e1[2];
                      dxp[i][2] = (xsp[0]-xray[0])*hcp_e2[0] + (xsp[1]-xray[1])*hcp_e2[1] + (xsp[2]-xray[2])*hcp_e2[2];
                    }
                    // this function returns xp, the location along the ray's hcp_e0...
                    double xp;
                    double yzp[2];
                    const double d2 = getPointToTriDist2Special(xp,yzp,dxp);
                    // this tri may have several bloated tris, each associated with a level
                    // from st_level down to max(level,1) -- this special case when dealing
                    // with level 0 is because there is no refinement windows associated with
                    // level 0.
                    double r = rmax;
                    for (int il = max(level,1); il <= st_level; ++il) {
                      // break out as soon as this tri's d2 no longer intersects the ray
                      const double dr2 = r*r - d2;
                      if (dr2 < 4.0/(hcp_dxp*hcp_dxp))
                        break;
                      // place the in point at xp - sqrt(r*r > d2)...
                      double xp_minus[3] = { xp-sqrt(dr2), 0.0, 0.0 };
                      double xtri[3];
                      MiscUtils::getClosestPointOnTriRobust(xtri,xp_minus,dxp[0],dxp[1],dxp[2]);
                      while (fabs(DIST(xtri,xp_minus)-r) > 1.0E-6*2.0*one_o_hcp_dxp) {
                        xp_minus[0] = xtri[0]-sqrt(r*r-xtri[1]*xtri[1]-xtri[2]*xtri[2]);
                        MiscUtils::getClosestPointOnTriRobust(xtri,xp_minus,dxp[0],dxp[1],dxp[2]);
                      }
                      double xp_plus[3] = { xp+sqrt(dr2), 0.0, 0.0 };
                      MiscUtils::getClosestPointOnTriRobust(xtri,xp_plus,dxp[0],dxp[1],dxp[2]);
                      while (fabs(DIST(xtri,xp_plus)-r) > 1.0E-6*2.0*one_o_hcp_dxp) {
                        xp_plus[0] = xtri[0]+sqrt(r*r-xtri[1]*xtri[1]-xtri[2]*xtri[2]);
                        MiscUtils::getClosestPointOnTriRobust(xtri,xp_plus,dxp[0],dxp[1],dxp[2]);
                      }
                      //double xmp[2];
                      //GeomUtils::calcBloatedTriIntersections(xmp,yzp,dxp[0],dxp[1],dxp[2],r);
                      //double xp_minus[3] = { xmp[0], 0.0, 0.0 };
                      //double xp_plus[3] = { xmp[1], 0.0, 0.0 };
                      // this is an in-out point for level "il" on a level "level" ray...
                      if (windowVec.size() >= threshold_windowVec_size) {
                        // here we assume collective synchronization amongst ranks...
#include "compressAndSend.hpp"
                        assert(windowVec.empty());
                        ++exchange_count;
                        // because some ranks may be done their rays before others, this logic
                        // keeps the calls to compressAndSend.hpp synchronized...
                        // my_done should be either 0 or 1...
                        int my_done = 0;
                        // any "0" has to return "0"...
                        int done;
                        MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);
                        assert(done == 0);
                      }
                      assert(il >= 1);
                      const int flag = ((triVec[ist].window+1)<<HCP_WINDOW_SHIFT_BITS)|il;
                      windowVec.push_back(InOutData(xp_minus[0],index,flag+1)); // use level+1 convention here: same as inoutVec
                      windowVec.push_back(InOutData(xp_plus[0],index,-flag-1));
                      my_window_count[0] += 2;
                      // and modify r for next level...
                      const int st_nlayers = nlowi[triVec[ist].window*(level_max+1)+il];
                      //r -= double(st_nlayers)*hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<il);
                      r -= double(st_nlayers)*hcp_delta/double(1<<il);
                    }
                  }
                }
              }
            }
          }
	  ++round_robin_rank;
	  if (round_robin_rank == mpi_size)
	    round_robin_rank = 0;
	}
      }

      // if this rank made it out here, we exchange until finished...
      done = 0;
      while (done != 1) {
#include "compressAndSend.hpp"
	assert(windowVec.empty());
        ++exchange_count;
        // my_done should be either 0 or 1...
        int my_done = 1;
        MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);
      }

      int8 window_count[3];
      MPI_Reduce(my_window_count,window_count,3,MPI_INT8,MPI_SUM,0,mpi_comm);
      int8 window_count_max;
      MPI_Reduce(my_window_count,&window_count_max,1,MPI_INT8,MPI_MAX,0,mpi_comm);

      if (mpi_rank == 0) {
	cout << " > in/out points found: " << window_count[0] <<
	  ", per rank avg/max: " << int8(double(window_count[0])/double(mpi_size)) << " " << window_count_max << endl;
        cout << " > in/out exchanges: " << exchange_count << endl;
        cout << " > tris flagged at or above this level: " << window_count[1] << " tris passing signed distance check: " <<
	  window_count[2] << " fraction: " << double(window_count[2])/double(window_count[1]) << endl;
	double wtime = MPI_Wtime();
        cout << " > bloated tri time: " << wtime-wtime_prev << " [s], total time so far: " << wtime-wtime0 << " [s]" << endl;
        wtime_prev = wtime;
      }

      delete[] my_ray_flag;

      // we can now build the points for this level...

      for (int ilevel = 0; ilevel <= level_max; ++ilevel)
        for (int iwindow = 0; iwindow <= nwi; ++iwindow)
          level_count[(nwi+1)*ilevel+iwindow] = 0; // level takes precendence
      iter_end = inoutVec[level].begin();
      while (iter_end != inoutVec[level].end()) {
	vector<InOutData>::iterator iter_begin = iter_end++;
	while ((iter_end != inoutVec[level].end())&&(iter_end->index == iter_begin->index))
	  ++iter_end;
	// we need the j and k for this index...
	// recall: int index = ((j-jminL)/djk)*((kmaxL-kminL)/djk+1)+((k-kminL)/djk);
	const int j = iter_begin->index/((kmaxL-kminL)/djk+1)*djk + jminL;
	const int k = (iter_begin->index-((j-jminL)/djk)*((kmaxL-kminL)/djk+1))*djk + kminL;
	// should already be sorted and level_count should already be zero...
        for (int ilevel = 0; ilevel <= level_max; ++ilevel)
          for (int iwindow = 0; iwindow <= nwi; ++iwindow)
            assert(level_count[(nwi+1)*ilevel+iwindow] == 0);
	int i1 = 0; // doesn't matter -- not used in below logic until after first time through
	for (vector<InOutData>::iterator iter = iter_begin; iter != iter_end; ++iter) {
	  //cout << "index: " << iter->index << " level: " << iter->level << " xp: " << iter->xp << endl;
	  const double xp_end = iter->xp;
	  const int i0 = i1;
	  i1 = int( (xp_end + hcp_offset_factor*double(j-k)*one_o_hcp_dxp*0.5*one_o_hcp_packing_factor)*(hcp_dxp*hcp_packing_factor) );
	  if (i1%2 == 0) i1 += 1;
	  // add points in this gap as appropriate...
	  // to start, find the highest level where we have a positive count...
          int ilevel, iwindow = nwi; // need to initialize because level_max could be  0
          bool found = false;
          //         for (ilevel = level_max; ilevel > 0; --ilevel) {
          //           for (iwindow = nwi; iwindow > 0; --iwindow) {
          //             if (level_count[(nwi+1)*ilevel+iwindow] > 0) {
          //         cout << level_count[0] << " " << level << " " << ilevel << " " << iwindow << endl;
          //             }
          //           }
          //         }
          for (ilevel = level_max; ilevel > 0; --ilevel) {
            for (iwindow = nwi; iwindow > 0; --iwindow) {
              if (level_count[(nwi+1)*ilevel+iwindow] > 0) {
                found = true;
                break;
              }
            }
            if (found)
              break;
          }
          //cout << "level,ilevel,iwindow: " << level << " " << ilevel << " " << iwindow << " " <<endl;
	  if ((level_count[0] > 0)&&(ilevel >= level)) {
	    // it is possible that with extreme optimization, the fact that xpIntersectionVec[ii].xp are sorted still
	    // doesn't ensure monotonic i0,i1, etc...
	    if (i1 < i0) i1 = i0;
	    // now build the imin/imax values for this level...
	    const int di = (1<<(Nbase-ilevel));
	    const int imin = ( i0 < 0 ? i0-(i0%di) : i0-(i0%di)+di );  assert(imin > i0);
	    const int imax = ( i1 < 0 ? i1-(i1%di)-di : i1-(i1%di) );  assert(imax < i1);
	    assert(imin%di == 0);
	    assert(imax%di == 0);
	    for (int i = imin; i <= imax; i += di) {
	      // TODO: modify for different hcp_e0/e1/e2 eventually...
	      assert(hcp_e0[0] == 1); assert(hcp_e0[1] == 0.0); assert(hcp_e0[2] == 0.0);
	      double x[3];
	      x[0] = hcp_x0pdx0[0] - hcp_offset_factor*double(j-k)*one_o_hcp_dxp*0.5*one_o_hcp_packing_factor + double(i)*one_o_hcp_dxp*one_o_hcp_packing_factor;
	      x[1] = hcp_x0pdx0[1] + double(j)*one_o_hcp_dxp;
	      x[2] = hcp_x0pdx0[2] + double(k)*one_o_hcp_dxp;

              // in an effort to break degeneracies at internal transitions: 
              if (Param * param = getParam("DEGENERACY_HACK")) {
                const double eps = param->getDouble(); 
                FOR_I3 x[i] += eps*(double(rand())/double(RAND_MAX)-0.5)*hcp_delta/double(1<<ilevel);
              }
              
              // only need points that are within thickened plane (if it exists)...
              if (plane_xp == NULL) {
                vertexVec.push_back(HcpVertex(x,ilevel,iwindow));
                my_vv_count[ilevel] += 1; // and level-based count
              }
              else {
                const double dx[3] = DIFF(x,plane_xp);
                if (fabs(DOT_PRODUCT(dx,unit_np)) <= hcp_delta*hcp_delta_factor*hcp_packing_factor/double(1<<ilevel)) {
                  vertexVec.push_back(HcpVertex(x,ilevel,iwindow));
                  my_vv_count[ilevel] += 1; // and level-based count
                }
              }
	    }
	  }
	  // and update the level count...
	  if (iter->level > 0) {
	    //const int il = iter->level-1;
	    //assert((il >= 0)&&(il <= level_max));
            const int il = (nwi+1)*((iter->level-1)&LEVEL_MASK_BITS)+((iter->level-1)>>HCP_WINDOW_SHIFT_BITS);
            //cout << il << " " << ((iter->level-1)&LEVEL_MASK_BITS) << " " << ((iter->level-1)>>HCP_WINDOW_SHIFT_BITS) << endl;
	    assert((il >= 0)&&(il < (nwi+1)*(level_max+1)));
	    ++level_count[il];
            // can go negative for a bit if the ff comes before the surface
	    //assert(level_count[il] == 1);
          }
          else {
	    assert(iter->level < 0);
	    //const int il = -iter->level-1;
	    //assert((il >= 0)&&(il <= level_max));
            const int il = (nwi+1)*((-iter->level-1)&LEVEL_MASK_BITS)+((-iter->level-1)>>HCP_WINDOW_SHIFT_BITS);
            //cout << il << " " << ((-iter->level-1)&LEVEL_MASK_BITS) << " " << ((-iter->level-1)>>HCP_WINDOW_SHIFT_BITS) << endl;
	    assert((il >= 0)&&(il < (nwi+1)*(level_max+1)));
	    --level_count[il];
            // can go negative for a bit if the ff comes before the surface
	    //assert(level_count[il] == 0);
	  }
	}
	// check that everything occured in pairs...
        for (int ilevel = 0; ilevel <= level_max; ++ilevel)
          for (int iwindow = 0; iwindow <= nwi; ++iwindow)
            assert(level_count[(nwi+1)*ilevel+iwindow] == 0);
      }

    }

    if (mpi_rank == 0) {
      DELETE(vv_count);
      vv_count = new int8[level_max+1];
    }
    MPI_Reduce(my_vv_count,vv_count,level_max+1,MPI_INT8,MPI_SUM,0,mpi_comm);
    delete[] my_vv_count;
    np_global = 0;
    if (mpi_rank == 0) {
      for (int level = 0; level <= level_max; ++level)
	np_global += vv_count[level];
      cout << " > hcp point count: " << np_global << endl;
      cout << " > hcp point distribution per level:" << endl;
      for (int level = 0; level <= level_max; ++level) {
        cout << " > level: " << level << " count: " << vv_count[level] << " (" << double(vv_count[level])/double(np_global)*100.0 << "%)" << endl;
      }
    }
    MPI_Bcast(&np_global,1,MPI_INT8,0,mpi_comm);
    
    // cleanup
    delete[] inoutVec;
    delete[] bbminmax;
    delete[] level_count;
    DELETE(nlowi);
    DELETE(leowi);

    //writeVertexVecToTecplot("vv.dat");
    //writeVertexVecToVtk("vv.vtk",np_global);
    
  }
	    
  void queryHcpWindows() const {
    WUI(INFO,"There are currently " << hcpWindowVec.size() << " HCP_WINDOWS:");
    for (int ii = 0; ii < hcpWindowVec.size(); ++ii) {
      WUI(INFO,"  " << ii << " " << hcpWindowVec[ii].str() << " from param: " << hcpWindowVec[ii].param_str);
    }
  }
  
private:

  void writeBinary(const string& filename) const {
    
    COUT1("HcpPointBuilder::writeBinary(): " << filename);
    
    // copy the relevant data into a Points object and write...
    // NOTE: this is a somewhat different version of the pbin binary
    // that supports level instead of delta...

    Points pts;
    pts.np = vertexVec.size();
    pts.xp = new double[pts.np][3];
    pts.setHcpData(hcp_packing,hcp_delta_factor,hcp_delta);
    for (int ip = 0; ip < pts.np; ++ip)
      FOR_I3 pts.xp[ip][i] = vertexVec[ip].x[i];
    pts.level = new int[pts.np];
    for (int ip = 0; ip < pts.np; ++ip)
      pts.level[ip] = vertexVec[ip].level;
    pts.writeBinary(filename); // note: this writes a special version of the points binary file based on level, not delta
    
  }

  void readBinary(const string& filename) {
    
    COUT1("HcpPointBuilder::readBinary(): " << filename);
    
    // when using READ_HCP_DATA, you should turn off ALL hcp-related stuff...

    assert(b_hcp_delta == false);
    assert(vertexVec.empty());
    assert(b_all_mesh_points == false);
    
    Points pts;
    pts.readBinary(filename);

    // global grid size...
    np_global = pts.np_global;

    // copy over the hcp data...
    assert(pts.b_hcp);
    b_hcp_delta = true;
    pts.getHcpData(hcp_packing,hcp_delta_factor,hcp_delta);
    
    // xp and level...
    vertexVec.resize(pts.np);
    for (int ip = 0; ip < pts.np; ++ip) {
      FOR_I3 vertexVec[ip].x[i] = pts.xp[ip][i];
      vertexVec[ip].level = pts.level[ip];
    }

    b_all_mesh_points = true;
    
  }

  void writeVertexVecToVtk(const string& filename,const int8 np_global) const {

    if (mpi_rank == 0) cout << "writeVertexVecToVtk: " << filename << endl;

    // do round-robin writing of points to vtk...

    // position
    FILE * fp;
    if ( mpi_rank == 0 ) {
      fp = fopen(filename.c_str(),"w");
      assert(fp != NULL);
      fprintf(fp,"# vtk DataFile Version 2.0\n");
      fprintf(fp,"3d scatter data\n");
      fprintf(fp,"ASCII\n\n");
      fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
      fprintf(fp,"POINTS %d double\n",(int)np_global);
    }
    else {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
    }
    for (int ivv = 0; ivv < vertexVec.size(); ++ivv) {
      fprintf(fp,"%18.15le %18.15le %18.15le\n",
	      vertexVec[ivv].x[0],
	      vertexVec[ivv].x[1],
	      vertexVec[ivv].x[2]);
    }
    if ( mpi_rank < mpi_size-1 ) {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }
    fclose(fp);
    MPI_Barrier(mpi_comm);

    // cell type
    if ( mpi_rank == 0 ) {
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
      fprintf(fp,"CELL_TYPES %lld\n",np_global);
    }
    else {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
    }
    for (int ivv = 0; ivv < vertexVec.size(); ++ivv) {
      fprintf(fp,"1\n");
    }
    if ( mpi_rank < mpi_size-1 ) {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }
    fclose(fp);
    MPI_Barrier(mpi_comm);

    // delta
    if ( mpi_rank == 0 ) {
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
      fprintf(fp,"POINT_DATA %lld\n",np_global);
      fprintf(fp,"SCALARS delta double\n");
      fprintf(fp,"LOOKUP_TABLE default\n");
    }
    else {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
    }
    for (int ivv = 0; ivv < vertexVec.size(); ++ivv) {
      fprintf(fp,"%18.15le\n",getPointsDeltaForLevel(vertexVec[ivv].level));
    }
    if ( mpi_rank < mpi_size-1 ) {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }
    fclose(fp);
    MPI_Barrier(mpi_comm);

    // level
    if ( mpi_rank == 0 ) {
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
      fprintf(fp,"SCALARS level int\n");
      fprintf(fp,"LOOKUP_TABLE default\n");
    }
    else {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
    }
    for (int ivv = 0; ivv < vertexVec.size(); ++ivv) {
      fprintf(fp,"%d\n",vertexVec[ivv].level);
    }
    if ( mpi_rank < mpi_size-1 ) {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }
    fclose(fp);
    MPI_Barrier(mpi_comm);

    // rank
    if ( mpi_rank == 0 ) {
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
      fprintf(fp,"SCALARS rank int\n");
      fprintf(fp,"LOOKUP_TABLE default\n");
    }
    else {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
    }
    for (int ivv = 0; ivv < vertexVec.size(); ++ivv) {
      fprintf(fp,"%d\n",mpi_rank);
    }
    if ( mpi_rank < mpi_size-1 ) {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }
    fclose(fp);
    MPI_Barrier(mpi_comm);
  }

  void writeVertexVecToTecplot(const string& filename) const {

    if (mpi_rank == 0) cout << "writeVertexVecToTecplot: " << filename << endl;

    // do round-robin writing of points to tecplot...

    FILE * fp;
    if ( mpi_rank == 0 ) {
      fp = fopen(filename.c_str(),"w");
      assert(fp != NULL);
      fprintf(fp,"TITLE = \"vertexVec\"\n");
      fprintf(fp,"VARIABLES = \"X\"\n");
      fprintf(fp,"\"Y\"\n");
      fprintf(fp,"\"Z\"\n");
      fprintf(fp,"\"DELTA\"\n");
      fprintf(fp,"\"LEVEL\"\n");
      fprintf(fp,"\"RANK\"\n");
    }
    else {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
    }

    //cout << "HACK ******************** writing mpi_rank in delta ****************************** HACK" << endl;

    for (int ivv = 0; ivv < vertexVec.size(); ++ivv) {
      fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %d %d\n",
	      vertexVec[ivv].x[0],
	      vertexVec[ivv].x[1],
	      vertexVec[ivv].x[2],
              getPointsDeltaForLevel(vertexVec[ivv].level),
              vertexVec[ivv].level,
              mpi_rank);
    }

    fclose(fp);

    if ( mpi_rank < mpi_size-1 ) {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }

    MPI_Barrier(mpi_comm);

  }

  double getPointToTriDist2Special(double& xp,double yzp[2],double dxp[3][3]) {

    // get the closest intersection point...

    // 2d areas in j,k...
    double A[3];
    double Asum = 0.0;
    FOR_I3 {
      A[i] = dxp[(i+1)%3][1]*dxp[(i+2)%3][2] - dxp[(i+1)%3][2]*dxp[(i+2)%3][1];
      Asum += A[i];
    }

    if (Asum == 0) {
      // degenerate...
      //cout << "A: " << COUT_VEC(A) << endl;
      // zero Asum means the tri is colinear...
      // find the largest area to figure out which 2 points are furthest out...
      if (fabs(A[0]) > max(fabs(A[1]),fabs(A[2]))) {
	// extremum points are x1,x2
        double dxe[2] = { dxp[1][1]-dxp[2][1], dxp[1][2]-dxp[2][2] };
        const double dp = dxe[0]*dxp[1][1] + dxe[1]*dxp[1][2];
        if (dp <= 0.0) {
          // we are closest to corner 1...
          xp = dxp[1][0];
          FOR_J2 yzp[j] = dxp[1][j+1];
          return dxp[1][1]*dxp[1][1] + dxp[1][2]*dxp[1][2];
        }
        else {
          const double dxe_mag2 = dxe[0]*dxe[0] + dxe[1]*dxe[1];
          if (dp > dxe_mag2) {
            // we are closest to corner 2...
            xp = dxp[2][0];
            FOR_J2 yzp[j] = dxp[2][j+1];
            return dxp[2][1]*dxp[2][1] + dxp[2][2]*dxp[2][2];
          }
          else {
            // this is an edge intersection...
            xp = dp/dxe_mag2*dxp[2][0] + (1.0-dp/dxe_mag2)*dxp[1][0];
            FOR_J2 yzp[j] = dp/dxe_mag2*dxp[2][j+1] + (1.0-dp/dxe_mag2)*dxp[1][j+1];
            //cout << "got an EDGE!" << endl;
            return dxp[1][1]*dxp[1][1] + dxp[1][2]*dxp[1][2] - dp*dp/dxe_mag2;
          }
        }
      }
      else if (abs(A[1]) > abs(A[2])) {
	// extremum points are x2,x0
        double dxe[2] = { dxp[2][1]-dxp[0][1], dxp[2][2]-dxp[0][2] };
        const double dp = dxe[0]*dxp[2][1] + dxe[1]*dxp[2][2];
        if (dp <= 0.0) {
          // we are closest to corner 2...
          xp = dxp[2][0];
          FOR_J2 yzp[j] = dxp[2][j+1];
          return dxp[2][1]*dxp[2][1] + dxp[2][2]*dxp[2][2];
        }
        else {
          const double dxe_mag2 = dxe[0]*dxe[0] + dxe[1]*dxe[1];
          if (dp > dxe_mag2) {
            // we are closest to corner 0...
            xp = dxp[0][0];
            FOR_J2 yzp[j] = dxp[0][j+1];
            return dxp[0][1]*dxp[0][1] + dxp[0][2]*dxp[0][2];
          }
          else {
            // this is an edge intersection...
            xp = dp/dxe_mag2*dxp[0][0] + (1.0-dp/dxe_mag2)*dxp[2][0];
            FOR_J2 yzp[j] = dp/dxe_mag2*dxp[0][j+1] + (1.0-dp/dxe_mag2)*dxp[2][j+1];
            //cout << "got an EDGE!" << endl;
            return dxp[2][1]*dxp[2][1] + dxp[2][2]*dxp[2][2] - dp*dp/dxe_mag2;
          }
        }
      }
      else if (abs(A[2]) > 0) {
	// extremum points are x0,x1
        double dxe[2] = { dxp[0][1]-dxp[1][1], dxp[0][2]-dxp[1][2] };
        const double dp = dxe[0]*dxp[0][1] + dxe[1]*dxp[0][2];
        if (dp <= 0.0) {
          // we are closest to corner 0...
          xp = dxp[0][0];
          FOR_J2 yzp[j] = dxp[0][j+1];
          return dxp[0][1]*dxp[0][1] + dxp[0][2]*dxp[0][2];
        }
        else {
          const double dxe_mag2 = dxe[0]*dxe[0] + dxe[1]*dxe[1];
          if (dp > dxe_mag2) {
            // we are closest to corner 1...
            xp = dxp[1][0];
            FOR_J2 yzp[j] = dxp[1][j+1];
            return dxp[1][1]*dxp[1][1] + dxp[1][2]*dxp[1][2];
          }
          else {
            // this is an edge intersection...
            xp = dp/dxe_mag2*dxp[1][0] + (1.0-dp/dxe_mag2)*dxp[0][0];
            FOR_J2 yzp[j] = dp/dxe_mag2*dxp[1][j+1] + (1.0-dp/dxe_mag2)*dxp[0][j+1];
            //cout << "got an EDGE!" << endl;
            return dxp[0][1]*dxp[0][1] + dxp[0][2]*dxp[0][2] - dp*dp/dxe_mag2;
          }
        }
      }
      else {
	// all areas must be zero. Points must be coincident, so return
	// distance to point 0...
	FOR_I3 assert(A[i] == 0);
        xp = dxp[0][0];
        FOR_J2 yzp[j] = dxp[0][j+1];
        return dxp[0][1]*dxp[0][1] + dxp[0][2]*dxp[0][2];
      }
    }
    else {

      FOR_I3 A[i] /= Asum;
      //cout << "A: " << COUT_VEC(A) << endl;

      if ((A[0] >= 0.0)&&(A[1] >= 0.0)&&(A[2] >= 0.0)) {
        // we are inside the tri. return d2 = 0, and use the
        // areas to return the x location...
        xp = A[0]*dxp[0][0] + A[1]*dxp[1][0] + A[2]*dxp[2][0];
        FOR_J2 yzp[j] = A[0]*dxp[0][j+1] + A[1]*dxp[1][j+1] + A[2]*dxp[2][j+1];
        return 0.0;
      }

      bool corner_flag[3] = { false, false, false };
      FOR_I3 {
        if (A[i] < 0.0) {
          // note we flip the sign of both vectors...
          double dxe[2] = { dxp[(i+1)%3][1]-dxp[(i+2)%3][1], dxp[(i+1)%3][2]-dxp[(i+2)%3][2] };
          const double dp = dxe[0]*dxp[(i+1)%3][1] + dxe[1]*dxp[(i+1)%3][2];
          if (dp <= 0.0) {
            // we will need to check corner (i+1)%3
            corner_flag[(i+1)%3] = true;
          }
          else {
            const double dxe_mag2 = dxe[0]*dxe[0] + dxe[1]*dxe[1];
            if (dp > dxe_mag2) {
              corner_flag[(i+2)%3] = true;
            }
            else {
              // this is an edge intersection...
              xp = dp/dxe_mag2*dxp[(i+2)%3][0] + (1.0-dp/dxe_mag2)*dxp[(i+1)%3][0];
              FOR_J2 yzp[j] = dp/dxe_mag2*dxp[(i+2)%3][j+1] + (1.0-dp/dxe_mag2)*dxp[(i+1)%3][j+1];
              //cout << "got an EDGE!" << endl;
              return dxp[(i+1)%3][1]*dxp[(i+1)%3][1] + dxp[(i+1)%3][2]*dxp[(i+1)%3][2] - dp*dp/dxe_mag2;
            }
          }
        }
      }

      // if we got here, we should have one corner set...
      int icorner = -1;
      double d2_min;
      FOR_I3 {
        if (corner_flag[i]) {
          const double d2 = dxp[i][1]*dxp[i][1] + dxp[i][2]*dxp[i][2];
          if ((icorner == -1)||(d2 < d2_min)) {
            icorner = i;
            d2_min = d2;
          }
        }
      }
      assert(icorner != -1);
      xp = dxp[icorner][0];
      FOR_J2 yzp[j] = dxp[icorner][j+1];
      return d2_min;

    }

  }

};

#endif
