#ifndef _CTI_CANVAS_HPP_
#define _CTI_CANVAS_HPP_

#include "CTI.hpp"
#include "PngImage.hpp"
#include "ImageMetadata.hpp"
#include "WriteImageData.hpp"
#include "PlaneData.hpp"
#include "CvImageMap.hpp"
#include "IntFlag.hpp"

enum vecType {
  SURFACE_VECTOR,
  ISOSURFACE_VECTOR,
};

class CtiCanvas {
public:

  int write_rank; // rank to build png on (default to 0, updated by StaticSolver)

private:

  float x0[3]; // top left corner of image in world coords
  float e0[3]; // i-direction (right) vector
  float e1[3]; // j-direction (down) vector
  float e2[3];

  float width; // image width in world coords
  //double x0[3]; // top left corner of image in world coords
  //double e0[3]; // i-direction (right) vector
  //double e1[3]; // j-direction (down) vector
  //double e2[3];

  //double width; // image width in world coords
  int ni,nj; // image pixel size: width,height

  int nI,nJ; // block size (image divided into 8x8 blocks)

  // because only a small part of an image is neccessarily
  // active at any time, the image is allocated in blocks...

  bool b_blank;
  vector< PlaneData<float> > vBlankPlaneData;

  bool b_blank_data_plane;
  PlaneData<float>* dataPlaneData;

  PlaneData<float>* projPlane;

  float range[2]; // volume data
  int n_range;
  bool b_range;
  float range_sVar[2]; // surface data
  int n_range_sVar;
  bool b_range_sVar;
  float range_pVar[2]; // particle data
  bool b_range_pVar;
  float range_iVar[2]; // iso data
  int n_range_iVar;
  bool b_range_iVar;
  bool b_showEdges;

  // buffers for storing PixelBlockData after reduction to a
  // single node (cacheImage routine). flushImage() writes
  // these buffers to png image on disk.  These buffers
  // can persist through multiple flushImage() calls to
  // allow for fast imaging of snapshot sequences (serial surfer)
  int pbuf_size;
  float * pbuf_float;
  uint2 * pbuf_uint2;

  class PixelBlockData {
  public:
    // surface data
    float snormal[64][3];
    float sdepth[64];
    float sdata[64];
    uint2 szone[64]; // unsigned short: 0 to 65535
    // mask surface data
    float mdepth[64];
    // other (plane,particle,etc.) data
    float pnormal[64][3];
    float pdepth[64];
    float pdata[64];
    float paux[64];
    uint2 pzone[64];
    // bits used by surface, mask surface and other data
    uint2 bits[64];
    PixelBlockData() {
      //cout << "PixelBlockData()" << endl;
      for (int ij = 0; ij < 64; ++ij) FOR_K3 snormal[ij][k] = 0.f;
      for (int ij = 0; ij < 64; ++ij) sdepth[ij] = 0.f;
      for (int ij = 0; ij < 64; ++ij) sdata[ij] = 0.f;
      for (int ij = 0; ij < 64; ++ij) szone[ij] = 65535;

      for (int ij = 0; ij < 64; ++ij) mdepth[ij] = 0.f;

      for (int ij = 0; ij < 64; ++ij) FOR_K3 pnormal[ij][k] = 0.f;
      for (int ij = 0; ij < 64; ++ij) pdepth[ij] = 0.f;
      for (int ij = 0; ij < 64; ++ij) pdata[ij] = 0.f;
      for (int ij = 0; ij < 64; ++ij) paux[ij] = 0.f;
      for (int ij = 0; ij < 64; ++ij) pzone[ij] = 65535;

      for (int ij = 0; ij < 64; ++ij) bits[ij] = 0;
    }
  };

  // a local version of the image and depth
  PixelBlockData ** pbd;

  ImageMetadata * imd;

  // volvis stuff...

  bool b_volvis;
  float volvis_aux_data[2]; // min/max or avg/std
  bool b_volvis_aux_data;
  bool b_volvis_surface;
  string volvis_type;

  class VolVisData {
  public:
    int ij,flag;
    float depth,v0,v1;
    VolVisData() {} // needed for resize()
    VolVisData(const int ij,const float depth,const bool show,const bool dir,const float v = 0.f) {
      this->ij = ij;
      this->depth = depth;
      v0 = v1 = v; // could hold some info to aid transparency comp
      // recall that flag 0,1,2 are reserved for internal faces. Here we use
      // flag = 3 as the base, then set show and dir using the next bits of flag...
      flag = 3;
      if (show) flag |= 4;
      if (dir) flag |= 8;
    }
    VolVisData(const int ij,const float depth,const float v0,const float v1,const int flag) {
      this->ij = ij;
      this->depth = depth;
      this->v0 = v0;
      this->v1 = v1;
      assert((flag >= 0)&&(flag <= 2)); // add -1 to indicate v1 is being hijacked by the distance off the line
      this->flag = flag;
    }
    bool operator<(const VolVisData& other) const {
      return (ij < other.ij) || ((ij == other.ij)&&(depth > other.depth));
    }
  };

  vector<VolVisData> volVisDataVec;

  // ultimately volvis stuff is stored in these buffers
  // on write_rank...

  int volvis_size;
  int * volvis_ij;
  float (*volvis_data)[2]; // intensity and depth

  // defines pixel behavior on flush...
  enum PixelFlagCases {
    BACKGROUND_PIXEL,
    SURFACE_PIXEL,
    SURFACE_DATA_PIXEL,
    INTERNAL_SURFACE_PIXEL,
    PARTICLE_DATA_PIXEL,
    ISO_DATA_PIXEL,
    VOLUME_DATA_PIXEL,
    MASKED_VOLUME_PIXEL,
    MASKED_AUX_VOLUME_PIXEL,
    MASKED_VOLUME_DATA_PIXEL,
    MASKED_AUX_VOLUME_DATA_PIXEL,
  };

public:

  CtiCanvas() : imd(NULL) {};

  CtiCanvas(const double target[3],const double camera[3],const double up[3],const double width,const int size[2], ImageMetadata *_imd) : imd(_imd) {

    float target_f[3]; FOR_I3 target_f[i] = target[i];
    float camera_f[3]; FOR_I3 camera_f[i] = camera[i];
    float up_f[3]; FOR_I3 up_f[i] = up[i];
    float width_f = width;
    init(target_f,camera_f,up_f,width_f,size);

  }

  CtiCanvas(const float target[3],const float camera[3],const float up[3],const float width,const int size[2], ImageMetadata *_imd) : imd(_imd) {

    init(target,camera,up,width,size);

  }

  void init(const float target[3],const float camera[3],const float up[3],const float width,const int size[2]) {

    b_blank = false;

    b_range      = false;
    n_range = 0;
    b_range_sVar = false;
    n_range_sVar = 0;
    b_range_pVar = false;
    b_range_iVar = false;
    n_range_iVar = 0;
    b_showEdges  = false;

    // who writes the image...
    write_rank = 0;

    // record width and image size...

    this->width = width;
    ni = size[0]; assert((ni >= 1)&&(ni < 16000));
    nj = size[1]; assert((nj >= 1)&&(nj < 16000));

    // build pixel vectors...

    const float e2_[3] = DIFF(target,camera);
    const float e2_mag = MAG(e2_);
    assert(e2_mag > 0.0); // camera and target cannot be coincident
    FOR_I3 e2[i] = e2_[i]*float(ni)/(width*e2_mag);

    float e0_[3] = CROSS_PRODUCT(e2,up);
    const float e0_mag = MAG(e0_); assert(e0_mag > 0.0); // problem with up

    float e1_[3] = CROSS_PRODUCT(e2,e0_);
    const float e1_mag = MAG(e1_); assert(e1_mag > 0.0);

    // now we can compute the top left corner of the image: this is pixel (i,j) = (0,0)...

    FOR_I3 x0[i] = camera[i] - 0.5*width*float(ni-1)/float(ni)*e0_[i]/e0_mag - 0.5*width*float(nj-1)/float(ni)*e1_[i]/e1_mag;

    // and scale the e0 and e1 vectors to produce the correct image width...

    FOR_I3 e0[i] = e0_[i]*float(ni)/(width*e0_mag);
    FOR_I3 e1[i] = e1_[i]*float(ni)/(width*e1_mag);

    // and allocate the image pixel storage...

    nI = ((ni-1)>>3)+1; assert((nI >= 1)&&(nI < 2000));
    nJ = ((nj-1)>>3)+1; assert((nJ >= 1)&&(nJ < 2000));

    pbd = new PixelBlockData*[nI*nJ];
    for (int ii = 0; ii < nI*nJ; ++ii) pbd[ii] = NULL;

    pbuf_size = 0;
    pbuf_float = NULL;
    pbuf_uint2 = NULL;

    b_volvis          = false;
    b_volvis_surface  = false;
    b_volvis_aux_data = false;
    volvis_type       = "LINEAR";
    volvis_ij         = NULL;
    volvis_data       = NULL;

    b_blank_data_plane = true;
    dataPlaneData      = NULL;
    projPlane          = NULL;
  }

  ~CtiCanvas() {

    if (pbd) {
      for (int ii = 0; ii < nI*nJ; ++ii) {
        if (pbd[ii] != NULL) delete pbd[ii];
      }
      delete[] pbd;
    }
    if (pbuf_float) delete[] pbuf_float;
    if (pbuf_uint2) delete[] pbuf_uint2;

    if (volvis_ij)   delete[] volvis_ij;
    if (volvis_data) delete[] volvis_data;

    if (dataPlaneData) delete dataPlaneData;
    
    if (projPlane) delete projPlane;
  }

  float calcPixelDepth(const float xp[3]) const;
  float calcPixelDepth(const double xp[3]) const;

  // pixel blanking methods
  void addBlankPlane(const float xp[3],const float np[3]);
  void insertBlankPlane(const float xp[3],const float np[3]);
  void deleteFirstBlankPlane();
  bool isBlanked(const float i,const float j,const float d) const;
  bool isBlankedSkipDataPlane(const float i,const float j,const float d,const bool on_data_plane) const;
  bool isBlankedSkipDataPlaneOrClipped(const float i,const float j,const float d,const bool on_data_plane,const float c=0.f,const float * range=NULL) const;
  bool isBlanked(const float i,const float j,const float d,const int iblank) const;
  bool isBehind(const float i,const float j,const float d) const;
  void getBlankPlanes(vector<PlaneData<float> >& blank_planes) const;
  void deleteBlankPlanes();
  void insertBlankPlanes(const vector<PlaneData<float> >& blank_planes);
  void setDataPlane(const float xp[3],const float np[3]);
  void setProjPlane();
  void removeDataPlane();
  void removeProjPlane();

  // New serial surfer imaging routines...
  bool isMaskPixel(const int &image_ij) const;
  bool isMaskPixel(const int IJ, const int ij) const;
  bool isSurfacePixel(const int &image_ij) const;

  void addCimPlaneData(const CvImageMap &cim, double * var);
  void addCimPlaneMesh(const CvImageMap &cim, float * d2_buf);
  void addCimSurfaceData(const CvImageMap &cim, double * var);

  void addCimPlaneLIC(const CvImageMap &cim, float * lic_image);
  void addCimSurfaceLIC(const CvImageMap &cim, float *lic_image);

  //Serial snapshot imaging routines
  void updateCimPlaneData(const CvImageMap &cim, double * var);
  void updateCimSurfaceData(const CvImageMap &cim, double * var);

  //LIC related routines...
  void projectVectorOnPlane(const CvImageMap &cim, double (*vol_vec)[3], const float np[3]);
  void projectVectorOnSurface(const CvImageMap &cim, double (*vol_vec)[3]);
  void projectVectorOnSurfaceAveraged(const CvImageMap &cim, double (*vol_vec)[3]);
  void projectVectorOnSphere(const CvImageMap &cim, double (*vol_vec)[3], double (*vol_xcc)[3]);

  void addParticle(const double xp[3],const double dp,const bool b_skip_blanking = false);
  void addParticle(const double xp[3],const double dp,const double vp,const bool b_skip_blanking = false);

  void addPlaneDataRvvPositive(const double * var,const double (*x_vv)[3],const double * r_vv,const int n,const double factor = 1.0);
  void addSurfaceMeshRvvPositive(const double (*x_vv)[3],const double * r_vv,const int *i_vv,const int *ivv_global,const int n,const double factor = 1.0);
  void addPlaneMeshRvvPositive(const double (*x_vv)[3],const double * r_vv,const int *i_vv,const int *ivv_global,const int n,const double factor = 1.0);
  void addPointsMeshRvvPositive(const double (*x_vv)[3],const double * r_vv,const int * i_vv,const int n,const double factor = 1.0);

  void addSurfaceData(const double * var,const double (*x_vv)[3],const double * r_vv,const int n,const double factor = 1.0);
  void addPartialSurfaceData(const double * var,const double (*x_vv)[3],const double * r_vv,const int * zone_vv,const int n,const double factor = 1.0);
  void rmPbdsWithoutSurfaceData();

  inline void computeBarycentricCoords(float (&lamda)[3],const int i,const int j,const float i2,const float j2,const float j1mj2,const float j0mj2,const float i2mi1,const float i0mi2,const float one_o_det) const;

  inline void computeBarycentricCoords(float& l1,float& l2,float& l3,const int i,const int j,const float i2,const float j2,const float j1mj2,const float j0mj2,const float i2mi1,const float i0mi2,const float one_o_det) const;

  inline void computeBarycentricConsts(float& j1mj2,float& j0mj2,float& i2mi1,float& i0mi2,float& one_o_det,const float i_s[3],const float j_s[3]) const;

  inline void computeBarycentricConsts(float& j1mj2,float& j0mj2,float& i2mi1,float& i0mi2,float& one_o_det,const float i0,const float i1,const float i2,const float j0,const float j1,const float j2) const;

  void addSurfaceTri(const double xt0[3],const double xt1[3],const double xt2[3],const uint2 zone,const bool hide = false,const uchar mesh=0,const bool on_data_plane = false,const double * normal=NULL,const double * data=NULL,const bool cell_flood=false,const double * range=NULL);
  void addSurfaceTriThickened(const double xt0[3],const double xt1[3],const double xt2[3],const uint2 zone,const bool hide = false,const uchar mesh=0,const bool on_data_plane = false,const double * normal=NULL,const double * data=NULL,const bool cell_flood=false,const double * range=NULL);

  void addInternalTri(const double xt0[3],const double xt1[3],const double xt2[3],const uint2 zone,const uchar mesh=0,const bool on_data_plane = false,const double * data=NULL,const bool cell_flood=false,const float * range=NULL);

  // volvis...
  void setVolVis() {
    b_volvis = true;
  }
  void setVolvisAuxData(const float &_f0, const float &_f1) {
    volvis_aux_data[0] = _f0;
    volvis_aux_data[1] = _f1;
    b_volvis_aux_data = true;
  }
  void setVolvisType(const string &_type) {
    volvis_type = _type;
  }
  void setVolvisSurface() {
    b_volvis_surface = true;
  }
  void addSurfaceTriVolVis(const double xt0[3],const double xt1[3],const double xt2[3],const bool show,const double n0=0.0,const double n1=0.0,const double n2=0.0);
  void addVoronoiPointVolVis(const double x_vv[3],const double r_vv,const float v);
  void addEdge(const double xe0[3],const double xe1[3],const double dx[3],const uint2 zone,const bool barbell = false);

  void addVector(const vecType vType,const double xe0[3],const double dx[3],const double factor,const uint2 zone);

  // image writing
  void writeImage(const string& prefix,const int index, const RGB_MODE rgb_mode);
  void writeImage(const string& filename, const RGB_MODE rgb_mode);

  void setBlankDataPlane(const bool blank);
  void setShowEdges(const bool show);
  void setDataRange(const float &_min, const float &_max);
  void setSurfaceDataRange(const float &_min, const float &_max);
  void setParticleDataRange(const float &_min, const float &_max);
  void setIsoDataRange(const float &_min, const float &_max);
  void setDataRangeBins(const int nbin);
  void setSurfaceDataRangeBins(const int nbin);
  void setIsoDataRangeBins(const int nbin);
  void getActiveZones(IntFlag &bfzone_flag,const IntFlag &szozn_i);

  // intensity calc for the volvis...
  void calcVolvisIntensityAndDepth(float &intensity,float &depth,vector<VolVisData>::iterator& iter_begin,vector<VolVisData>::iterator& iter_end) const;

  // low level image writing routines, generally avoid using these
  // unless doing something special like snapshot imaging...
  void cacheImage();
  void flushImage(const string& filename, const RGB_MODE rgb_mode);
  void clearPBuf();

  // getters
  void getX0(float _x0[3]) const;
  void getE0(float _e0[3]) const;
  void getE1(float _e1[3]) const;
  void getE2(float _e2[3]) const;
  int getNi() const;
  int getNj() const;
  float getWidth() const;
  float getDepth(const int &i, const int &j) const;

private:

  void getIJij(int& IJ,int& ij,const int i, const int j);

  bool getPixelPauxBits(float& paux,uint2& bits,const int IJ,const int ij) const;
  bool getPixelPauxPzoneBits(float& paux,uint2& pzone,uint2& bits,const int IJ,const int ij) const;
  bool getPixelPdepthBits(float& pdepth,uint2& bits,const int IJ,const int ij) const;
  bool getPixelSdepthMdepthBits(float& sdepth,float& mdepth,uint2& bits,const int IJ,const int ij) const;
  bool getPixelSdepthBits(float& sdepth,uint2& bits,const int IJ,const int ij) const;
  bool getPixelSdepthPauxBits(float& sdepth,float& paux,uint2& bits,const int IJ,const int ij) const;
  bool getPixelSdepthPauxSzoneBits(float& sdepth,float& paux,uint2& szone,uint2& bits,const int IJ,const int ij) const;

  void setPixelSdataBits(const int IJ,const int ij,const float sdata,const uint2 bits);
  void setPixelSnormalSdepthSzone(const int IJ,const int ij,const float snormal[3],const float sdepth,const uint2 szone);
  void setPixelSnormalSdepthSzoneBits(const int IJ,const int ij,const float snormal[3],const float sdepth,const uint2 szone,const uint2 bits);
  void setPixelSnormalSdepthSdataSzone(const int IJ,const int ij,const float snormal[3],const float sdepth,const float sdata,const uint2 szone);
  void setPixelSnormalSdepthSdataSzoneBits(const int IJ,const int ij,const float snormal[3],const float sdepth,const float sdata,const uint2 szone,const uint2 bits);

  void setPixelMdepth(const int IJ,const int ij,const float mdepth);
  void setPixelMdepthBits(const int IJ,const int ij,const float mdepth,const uint2 bits);

  void setPixelPnormalPdepthPdataPauxPzoneBits(const int IJ,const int ij,const float pnormal[3],const float pdepth,const float pdata,const float paux,const uint2 pzone,const uint2 bits);
  void setPixelPnormalPdepthPdataPauxBits(const int IJ,const int ij,const float pnormal[3],const float pdepth,const float pdata,const float paux,const uint2 bits);
  void setPixelPnormalPdepthPdataPauxPzone(const int IJ,const int ij,const float pnormal[3],const float pdepth,const float pdata,const float paux,const uint2 pzone);
  void setPixelPnormalPdepthPdataPaux(const int IJ,const int ij,const float pnormal[3],const float pdepth,const float pdata,const float paux);
  void setPixelPnormalPdepthPzoneBits(const int IJ,const int ij,const float pnormal[3],const float pdepth,const uint2 pzone,const uint2 bits);
  void setPixelPnormalPdepthPzone(const int IJ,const int ij,const float pnormal[3],const float pdepth,const uint2 pzone);
  void setPixelPnormalPdepthPdataPzone(const int IJ,const int ij,const float pnormal[3],const float pdepth,const float pdata,const uint2 pzone);
  void setPixelPnormalPdepthPdataPzoneBits(const int IJ,const int ij,const float pnormal[3],const float pdepth,const float pdata,const uint2 pzone,const uint2 bits);

  void setPixelSdataPaux(const int IJ,const int ij,const float sdata,const float paux);
  void setPixelSdataPauxBits(const int IJ,const int ij,const float sdata,const float paux,const uint2 bits);

};

#endif
