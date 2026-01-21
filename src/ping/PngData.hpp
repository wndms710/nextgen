#ifndef PNGDATA_HPP
#define PNGDATA_HPP

#include <iostream>
#include <stack>
#include "PngImage.hpp"
#include "ColorMap.hpp"
#include "ImageLegend.hpp"
#include "MinHeap.hpp"
#include <set>

using namespace std;

struct pixel_type {
  int flag;
  string name;
  string varId;
  double rmin;
  double rmax;
} ;

class PngData : public PngImage {

public:

  enum PixelFlagCases {
    FLAG_WHITE   = -10, //keep continous
    FLAG_RED2    = -11, //...
    FLAG_GREEN2  = -12, //...
    FLAG_YELLOW  = -13, //...
    FLAG_BLUE2   = -14, //...
    FLAG_ORANGE  = -15, //...
    FLAG_PURPLE  = -16, //...
    FLAG_CYAN    = -17, //...
    FLAG_MAGENTA = -18, //...
    FLAG_PINK    = -19, //...
    FLAG_TEAL    = -20, //keep continous
    FLAG_BLANK   =  -6,
    FLAG_GRAY    =  -5,
    FLAG_BLUE    =  -4,
    FLAG_GREEN   =  -3,
    FLAG_RED     =  -2,
    BACKGROUND_PIXEL=-1,
    SURFACE_PIXEL,
    VOLUME_DATA_PIXEL,    //keep ordering TODO: are these types written in the image? should be explicitly defined as above
    SURFACE_DATA_PIXEL,   //..
    PARTICLE_DATA_PIXEL,  //..
    ISO_DATA_PIXEL        //keep ordering
  };

  bool b_readImage;

  // for simplicity, we're going to convert the
  // rgb pixels that we have to gray scale
  double * pixel_data;
  int * pixel_flag;
  int id;     // unique identifier
  bool lock;

  PngDataChunk * zoNe;
  PngDataChunk * dpTh;
  PngDataChunk * daTa;
  PngDataChunk * lgFr;

  bool b_data;
  bool b_dataSurface;
  bool b_dataParticles;
  bool b_dataIso;

  int npx;   // number of pixels that have non-trivial info...
  int npxS;  // number of pixels that have non-trivial surface info...
  int npxP;  // number of pixels that have non-trivial particle info...
  int npxI;  // number of pixels that have non-trivial iso info...

  ColorMap theColorMap;
  ColorMap theColorMapSurface;
  ColorMap theColorMapParticles;
  ColorMap theColorMapIso;

  ImageLegend theLegend;
  ImageLegend theLegendSurface;
  ImageLegend theLegendParticles;
  ImageLegend theLegendIso;

  string varId;
  string varIdSurface;
  string varIdParticles;
  string varIdIso;
  
  double range[2];
  double rangeSurface[2];
  double rangeParticles[2];
  double rangeIso[2];
  
  PngData() : PngImage(), b_readImage(true), pixel_data(NULL), pixel_flag(NULL), id(-1), lock(false), zoNe(NULL), dpTh(NULL), daTa(NULL), lgFr(NULL), b_data(false), b_dataSurface(false), b_dataParticles(false), b_dataIso(false), npx(0), npxS(0), npxP(0), npxI(0), theColorMap(), theColorMapSurface(), theColorMapParticles(), theColorMapIso(), theLegend(), theLegendSurface(), theLegendParticles(), theLegendIso(), varId(""), varIdSurface(""), varIdParticles(""), varIdIso("") {

    range[0] = 0.0;
    range[1] = 255.0;
    rangeSurface[0] = 0.0;
    rangeSurface[1] = 255.0;
    rangeParticles[0] = 0.0;
    rangeParticles[1] = 255.0;
    rangeIso[0] = 0.0;
    rangeIso[1] = 255.0;

  }

  PngData(int ni, int nj) : PngImage(ni,nj), b_readImage(false), pixel_data(NULL), pixel_flag(NULL), id(-1), lock(false), zoNe(NULL), dpTh(NULL), daTa(NULL), lgFr(NULL), b_data(false), b_dataSurface(false), b_dataParticles(false), b_dataIso(false), npx(0), npxS(0), npxP(0), npxI(0), theColorMap(), theColorMapSurface(), theColorMapParticles(), theColorMapIso(), theLegend(), theLegendSurface(), theLegendParticles(), theLegendIso(), varId(""), varIdSurface(""), varIdParticles(""), varIdIso("") {

    //not enabling PngImage's memory cleanup option here because it will delete a non-null metadata object
    //we want to be able to construct this outside of this class and attach it prior to writing.
    //PngData will cleanup the buffer's/chunks here after checking b_readImage.

    range[0] = -HUGE_VAL;
    range[1] =  HUGE_VAL;
    rangeSurface[0] = -HUGE_VAL;
    rangeSurface[1] =  HUGE_VAL;
    rangeParticles[0] = -HUGE_VAL;
    rangeParticles[1] =  HUGE_VAL;
    rangeIso[0] = -HUGE_VAL;
    rangeIso[1] =  HUGE_VAL;

    metadata = NULL;
    buffer = new unsigned char[ni*nj][3];

    zoNe = new PngDataChunk("zoNe", ni, nj, 2);
    addPngDataChunk(zoNe);
    dpTh = new PngDataChunk("dpTh", ni, nj, 2);
    addPngDataChunk(dpTh);
    daTa = new PngDataChunk("daTa", ni, nj, 1);
    addPngDataChunk(daTa);
    lgFr = new PngDataChunk("lgFr", ni, nj, 1);
    addPngDataChunk(lgFr);



    pixel_data = new double[nx*ny];
    pixel_flag = new int[nx*ny];

    for (int ipx =0; ipx < nx*ny ; ++ipx) {
      buffer[ipx][0] = (unsigned char) 73;
      buffer[ipx][1] = (unsigned char) 175;
      buffer[ipx][2] = (unsigned char) 205;
      pixel_data[ipx] = 0;
      pixel_flag[ipx] = BACKGROUND_PIXEL;

      unsigned char zoNe_uc[2] = {(unsigned char) 255, (unsigned char) 255};
      zoNe->set(ipx,zoNe_uc);

      unsigned char dpTh_uc[2] = {(unsigned char) 255, (unsigned char) 255};
      dpTh->set(ipx,dpTh_uc);

      unsigned char daTa_uc = (unsigned char) 255;
      daTa->set(ipx,&daTa_uc);
 
      unsigned char lgFr_uc = (unsigned char) 255;
      lgFr->set(ipx,&lgFr_uc);

    }//ipx

  }

  ~PngData() {
    if (pixel_data != NULL) delete[] pixel_data;
    if (pixel_flag != NULL) delete[] pixel_flag;
    if (!b_readImage){
      if (buffer != NULL) delete[] buffer;
      if (zoNe != NULL) delete zoNe;
      if (dpTh != NULL) delete dpTh;
      if (daTa != NULL) delete daTa;
      if (lgFr != NULL) delete lgFr;
    }
  }

  void resize(const int nx_new,const int ny_new) {
    if (nx_new*ny_new < nx*ny) {
      // for smaller total memory, just resize the zoNe, dpTh, daTa, and lgFr...
      if (zoNe != NULL) zoNe->resize(nx_new,ny_new);
      if (dpTh != NULL) dpTh->resize(nx_new,ny_new);
      if (daTa != NULL) daTa->resize(nx_new,ny_new);
      if (lgFr != NULL) lgFr->resize(nx_new,ny_new);
      PngImage::resize(nx_new,ny_new);
    }
    else {
      // need to increase all buffers...
      assert(0);
    }
  }
    
  void crop(const int i0,const int j0,const int i1,const int j1) {
    const int nx_new = i1-i0+1;
    const int ny_new = j1-j0+1;
    assert(nx_new*ny_new < nx*ny); // crop needs to be smaller
    assert(buffer);
    assert(pixel_data);
    assert(pixel_flag);
    assert(zoNe);
    assert(dpTh);
    int ij_new = 0;
    for (int j = j0; j <= j1; ++j) {
      for (int i = i0; i <= i1; ++i) {
        const int ij_old = j*nx+i;
        buffer[ij_new][0] = buffer[ij_old][0];
        buffer[ij_new][1] = buffer[ij_old][1];
        buffer[ij_new][2] = buffer[ij_old][2];
        pixel_data[ij_new] = pixel_data[ij_old];
        pixel_flag[ij_new] = pixel_flag[ij_old];
        unsigned char zoNe_uc[2];
        zoNe->get(ij_old,zoNe_uc);
        zoNe->set(ij_new,zoNe_uc);
        unsigned char dpTh_uc[2];
        dpTh->get(ij_old,dpTh_uc);
        dpTh->set(ij_new,dpTh_uc);
        if (daTa) {
          unsigned char daTa_uc;
          daTa->get(ij_old,&daTa_uc);
          daTa->set(ij_new,&daTa_uc);
        }
        if (lgFr) {
          unsigned char lgFr_uc;
          lgFr->get(ij_old,&lgFr_uc);
          lgFr->set(ij_new,&lgFr_uc);
        }
        ++ij_new;
      }
    }
    // adjust the metadata for the new image offset...
    // recall:
    //  xp[i] = ip*metadata->transformMat[0+i] +
    //          jp*metadata->transformMat[4+i] +
    //          dpTh_us*metadata->transformMat[8+i] +
    //          metadata->transformMat[12+i];
    assert(metadata);
    for (int i = 0; i < 3; ++i) {
      // comprehend the ishift and jshift in the offset...
      metadata->transformMat[12+i] += 
        i0*metadata->transformMat[0+i] +
        j0*metadata->transformMat[4+i];
    }
    // and set the size to the new size...
    resize(nx_new,ny_new);
  }

  void getRGB(unsigned char rgb[3],const int ipx) {

    assert(buffer);
    assert((ipx >= 0)&&(ipx < nx*ny));
    FOR_I3 rgb[i] = buffer[ipx][i];

  }

  void setBackgroundColor(const unsigned char r,const unsigned char g,const unsigned char b) {

    assert(buffer);
    
    if (pixel_flag) {
      for (int ipx = 0; ipx < nx*ny ; ++ipx) {
        if (pixel_flag[ipx] == BACKGROUND_PIXEL) {
          buffer[ipx][0] = r;
          buffer[ipx][1] = g;
          buffer[ipx][2] = b;
        }
      }
    }
    else if (zoNe) {
      for (int ipx = 0; ipx < nx*ny ; ++ipx) {
        unsigned char zoNe_uc[2];
        zoNe->get(ipx,zoNe_uc);
        if ((zoNe_uc[0] == 255)&&(zoNe_uc[1] == 255)) {
          buffer[ipx][0] = r;
          buffer[ipx][1] = g;
          buffer[ipx][2] = b;
        }
      }
    }
    else {
      for (int ipx = 0; ipx < nx*ny ; ++ipx) {
        // look for bahama blue...
        if ((buffer[ipx][0] == 73)&&(buffer[ipx][1] == 175)&&(buffer[ipx][2] == 205)) {
          buffer[ipx][0] = r;
          buffer[ipx][1] = g;
          buffer[ipx][2] = b;
        }
      }
    }
      
  }

  void finalizeRead() {

    //images must have valid chunks now for use with ping_batch...

    assert ( b_readImage );  ///don't call this if you are just creating an empty PngData object
    assert ( pixel_data == NULL );
    assert ( pixel_flag == NULL );
    assert ( zoNe == NULL );
    assert ( dpTh == NULL );
    assert ( daTa == NULL );
    assert ( lgFr == NULL );

    zoNe = getPngDataChunk("zoNe");
    dpTh = getPngDataChunk("dpTh");
    if (!zoNe) {
      CERR("Zone information not found");
    }
    if (!dpTh) {
      CERR("Depth information not found");
    }

    pixel_data = new double[nx*ny];
    pixel_flag = new int[nx*ny];

    daTa = getPngDataChunk("daTa");
    lgFr = getPngDataChunk("lgFr");

    //set data type flags from metadata...
    b_data = b_dataSurface = b_dataParticles = b_dataIso = false;

    if (metadata->hasText("VAR") && metadata->getVarId() != "OFF" &&
        metadata->hasText("RANGE_MAX") && metadata->hasText("RANGE_MIN") &&
        metadata->hasText("COLORMAP")){
      range[0] = metadata->getRangeMin();
      range[1] = metadata->getRangeMax();
      varId    = metadata->getVarId();

      theColorMap.set_type(metadata->getColorMapName());
      //the daTa chunk is always included for other data types, but only for planar
      //data when another data type is present.  If no daTa chunk is present we
      //must get the data value from rgb buffer. TODO: make more maps invertible
      if (theColorMap.set_type(metadata->getColorMapName())){
        if (daTa||theColorMap.get_b_dataFromRgb())
          b_data = true;
        else
          COUT1("Unable to extract data values from volume/planar colormap " << metadata->getColorMapName() << ", ignoring volume/planar data");
      }
      else{
        COUT1("Unrecognized volume/planar data colormap " << metadata->getColorMapName() << ", ignoring volume/planar data");
      }
    }
    // else image has no valid planar data

    if (metadata->hasText("VAR_ON_SURFACE") && metadata->getSurfVarId() != "OFF" &&
        metadata->hasText("RANGE_ON_SURFACE_MAX") && metadata->hasText("RANGE_ON_SURFACE_MIN") &&
        metadata->hasText("COLORMAP_SURFACE")){
      rangeSurface[0] = metadata->getSurfRangeMin();
      rangeSurface[1] = metadata->getSurfRangeMax();
      varIdSurface    = metadata->getSurfVarId();
      if(theColorMapSurface.set_type(metadata->getSurfColorMapName()))
        b_dataSurface = true;
      else
        COUT1("Unrecognized surface data colormap " << metadata->getSurfColorMapName() << ", ignoring surface data");
    }
    // else image has no valid surface data

    if (metadata->hasText("VAR_ON_PARTICLE") && metadata->getPartVarId() != "OFF" &&
        metadata->hasText("RANGE_ON_PARTICLE_MAX") && metadata->hasText("RANGE_ON_PARTICLE_MIN") &&
        metadata->hasText("COLORMAP_PARTICLE")){
      rangeParticles[0] = metadata->getPartRangeMin();
      rangeParticles[1] = metadata->getPartRangeMax();
      varIdParticles    = metadata->getPartVarId();
      if (theColorMapParticles.set_type(metadata->getPartColorMapName()))
        b_dataParticles = true;
      else
        COUT1("Unrecognized particle data colormap " << metadata->getPartColorMapName() << ", ignoring particle data");
    }
    // else image has no valid particle data

    if (metadata->hasText("VAR_ON_ISO") && metadata->getIsoVarId() != "OFF" &&
        metadata->hasText("RANGE_ON_ISO_MAX") && metadata->hasText("RANGE_ON_ISO_MIN") &&
        metadata->hasText("COLORMAP_ISO")){
      rangeIso[0] = metadata->getIsoRangeMin();
      rangeIso[1] = metadata->getIsoRangeMax();
      varIdIso    = metadata->getIsoVarId();
      if (theColorMapIso.set_type(metadata->getIsoColorMapName()))
        b_dataIso = true;
      else
        COUT1("Unrecognized iso data colormap " << metadata->getIsoColorMapName() << ", ignoring iso data");
    }
    // else image has no valid iso data

   //if surface, particle, or iso data exists, then daTa chunk should be present...
    if (b_dataSurface || b_dataParticles || b_dataIso){
      if(!daTa){
        CERR("Image missing daTa record");
      }
      else if (!lgFr){
        COUT1("Warning: missing lighting factor record, lighting will be reconstructed where possible from daTa and RGB");
      }
    }

    //build pixel_flag and pixel_data
    npx = npxP = npxS = npxI = 0;

    for (int ipx =0; ipx < nx*ny ; ++ipx) {

      //zoNe chunk must exist to proceed
      unsigned char zoNe_uc[2];
      zoNe->get(ipx,zoNe_uc);
      uint2 zoNe_us = (zoNe_uc[1] << 8 | zoNe_uc[0]); //little endian

      if (zoNe_us == 65535){
        pixel_flag[ipx] = BACKGROUND_PIXEL;
        pixel_data[ipx] = 0.0;
      }
      else if (zoNe_us < 32768 || (!b_dataSurface && zoNe_us<65532)){
        pixel_flag[ipx] = SURFACE_PIXEL; //could be edges, boundary or internal surface
        pixel_data[ipx] = 0.0;
      }
      else{ //data pixel

        double grayval = 0.0;
        if (daTa){
          unsigned char daTa_uc;
          daTa->get(ipx,&daTa_uc);
          grayval = ((int)daTa_uc);
        }

        if (zoNe_us == 65534){
          if (daTa){
            pixel_flag[ipx] = VOLUME_DATA_PIXEL;
            pixel_data[ipx] = grayval/255.0 * (range[1]-range[0]) + range[0];
            ++npx;
          }
          else if (b_data){
            grayval =  255.0*theColorMap.calcPhiFromRgb(buffer[ipx]);
            pixel_flag[ipx] = VOLUME_DATA_PIXEL;
            pixel_data[ipx] = grayval/255.0 * (range[1]-range[0]) + range[0];
            ++npx;
          }
        }
        else if (zoNe_us == 65533){
          pixel_flag[ipx] = PARTICLE_DATA_PIXEL;
          pixel_data[ipx] = grayval/255.0 * (rangeParticles[1]-rangeParticles[0]) + rangeParticles[0];
          ++npxP;
        }
        else if (zoNe_us == 65532){
          pixel_flag[ipx] = ISO_DATA_PIXEL;
          pixel_data[ipx] = grayval/255.0 * (rangeIso[1]-rangeIso[0]) + rangeIso[0];
          ++npxI;
        }
        else if (b_dataSurface && zoNe_us >= 32768){
          pixel_flag[ipx] = SURFACE_DATA_PIXEL;
          pixel_data[ipx] = grayval/255.0 * (rangeSurface[1]-rangeSurface[0]) + rangeSurface[0];
          ++npxS;
        }
        else
          cerr << "unrecognized data pixel type" << endl;
      }

    }//ipx

  }//finalizeRead

  //return a vector of pixel_type for the current image
  void getPixelTypes(vector<pixel_type> &pixelTypeVec){
    assert(b_readImage); //metadata not present for empty PngData objects

    if (b_data){
      pixelTypeVec.push_back(pixel_type());
      pixelTypeVec.back().flag = VOLUME_DATA_PIXEL;
      pixelTypeVec.back().name = "volume";
      pixelTypeVec.back().varId = metadata->getVarId();
      pixelTypeVec.back().rmin = metadata->getRangeMin();
      pixelTypeVec.back().rmax = metadata->getRangeMax();
    }
    if (b_dataSurface){
      pixelTypeVec.push_back(pixel_type());
      pixelTypeVec.back().flag = SURFACE_DATA_PIXEL;
      pixelTypeVec.back().name = "surface";
      pixelTypeVec.back().varId = metadata->getSurfVarId();
      pixelTypeVec.back().rmin = metadata->getSurfRangeMin();
      pixelTypeVec.back().rmax = metadata->getSurfRangeMax();
    }
    if (b_dataParticles){
      pixelTypeVec.push_back(pixel_type());
      pixelTypeVec.back().flag = PARTICLE_DATA_PIXEL;
      pixelTypeVec.back().name = "particles";
      pixelTypeVec.back().varId = metadata->getPartVarId();
      pixelTypeVec.back().rmin = metadata->getPartRangeMin();
      pixelTypeVec.back().rmax = metadata->getPartRangeMax();
    }
    if (b_dataIso){
      pixelTypeVec.push_back(pixel_type());
      pixelTypeVec.back().flag = ISO_DATA_PIXEL;
      pixelTypeVec.back().name = "iso";
      pixelTypeVec.back().varId = metadata->getIsoVarId();
      pixelTypeVec.back().rmin = metadata->getIsoRangeMin();
      pixelTypeVec.back().rmax = metadata->getIsoRangeMax();
    }
  }

  int getNx() const { return nx ; }
  int getNy() const { return ny ; }

  double getRangeMin() const { return range[0]; }
  double getRangeMax() const { return range[1]; }
  double getRangeSurfaceMin() const { return rangeSurface[0]; }
  double getRangeSurfaceMax() const { return rangeSurface[1]; }
  double getRangeParticlesMin() const { return rangeParticles[0]; }
  double getRangeParticlesMax() const { return rangeParticles[1]; }
  double getRangeIsoMin() const { return rangeIso[0]; }
  double getRangeIsoMax() const { return rangeIso[1]; }

  void setRange(const double _range[2]) {
    range[0] = _range[0];
    range[1] = _range[1];
  }
  void setRangeSurface(const double _range[2]) {
    rangeSurface[0] = _range[0];
    rangeSurface[1] = _range[1];
  }
  void setRangeParticles(const double _range[2]) {
    rangeParticles[0] = _range[0];
    rangeParticles[1] = _range[1];
  }
  void setRangeIso(const double _range[2]) {
    rangeIso[0] = _range[0];
    rangeIso[1] = _range[1];
  }

  double getTime() const {
    assert(b_readImage); //metadata not present for empty PngData objects
    return metadata->getTime();
  }

  uint2 getDepth(int ipx) const {
    assert(ipx<(nx*ny));
    uint2 dpTh_us = 0;
    unsigned char dpTh_uc[2];
    dpTh->get(ipx,dpTh_uc);
    dpTh_us = (dpTh_uc[1] << 8 | dpTh_uc[0]); //little endian
    return dpTh_us;
  }
  void setDepth(int ipx,const uint2 &dpTh_us) {
    assert(ipx<(nx*ny));
    unsigned char dpTh_uc[2];
    dpTh_uc[0] = dpTh_us & 255;
    dpTh_uc[1] = dpTh_us >> 8;
    dpTh->set(ipx,dpTh_uc);
  }

  void getZoneNames(vector<string> &_zoneNames){
    assert(_zoneNames.size()==0);
    for (int iz=0; iz<metadata->zoneNames.size(); ++iz){
      _zoneNames.push_back(metadata->zoneNames[iz]);
    }
  }

  void getZoneIds(vector<uint2> &_zoneIds){
    assert(_zoneIds.size()==0);
    for (int iz=0; iz<metadata->zoneIds.size(); ++iz){
      _zoneIds.push_back(metadata->zoneIds[iz]);
    }
  }

  //metadata object is updated in initializeWrite
  //colormap corresponding to buffer may be needed
  //after output colormap is changed.
  void setColorMap(string _colorMapName){
    theColorMap.set_type(_colorMapName);
  }
  void setColorMapSurface(string _colorMapName){
    theColorMapSurface.set_type(_colorMapName);
  }
  void setColorMapParticles(string _colorMapName){
    theColorMapParticles.set_type(_colorMapName);
  }
  void setColorMapIso(string _colorMapName){
    theColorMapIso.set_type(_colorMapName);
  }

  void setLegend(const double irel, const double jrel){
    theLegend.init(nx, ny, irel, jrel, theColorMap);
  }
  void setLegendSurface(const double irel, const double jrel){
    theLegendSurface.init(nx, ny, irel, jrel, theColorMapSurface);
  }
  void setLegendParticles(const double irel, const double jrel){
    theLegendParticles.init(nx, ny, irel, jrel, theColorMapParticles);
  }
  void setLegendIso(const double irel, const double jrel){
    theLegendIso.init(nx, ny, irel, jrel, theColorMapIso);
  }

  void rescaleRangesToData(){
    //range: plane, surface, particle, iso
    double data_range[8] = {HUGE_VAL,-HUGE_VAL,
                            HUGE_VAL,-HUGE_VAL,
                            HUGE_VAL,-HUGE_VAL,
                            HUGE_VAL,-HUGE_VAL};

    for (int ipx = 0; ipx < nx*ny ; ++ipx) {
      if ( pixel_flag[ipx] > SURFACE_PIXEL ) {
        int imin = -1;
        int imax = -1;
        if (pixel_flag[ipx] == VOLUME_DATA_PIXEL){
          imin=0; imax=1;
        }
        else if (pixel_flag[ipx] == SURFACE_DATA_PIXEL){
          imin=2; imax=3;
        }
        else if (pixel_flag[ipx] == PARTICLE_DATA_PIXEL){
          imin=4; imax=5;
        }
        else if (pixel_flag[ipx] == ISO_DATA_PIXEL){
          imin=6; imax=7;
        }
        assert(imin!=imax);//pixels > SURFACE_PIXEL should be one of the four above
        if(pixel_data[ipx]<data_range[imin]){
          data_range[imin] = pixel_data[ipx];
        }
        if(pixel_data[ipx]>data_range[imax]){
          data_range[imax] = pixel_data[ipx];
        }
      }
    }
    range[0]          = data_range[0];
    range[1]          = data_range[1];
    rangeSurface[0]   = data_range[2];
    rangeSurface[1]   = data_range[3];
    rangeParticles[0] = data_range[4];
    rangeParticles[1] = data_range[5];
    rangeIso[0]       = data_range[6];
    rangeIso[1]       = data_range[7];
  }

  //Update metadata with any changes to range, recoloring is now handled separately through initializeWriteRecolor
  void initializeWrite() {

    assert(b_readImage); //metadata not present for empty PngData objects, use initializeWrite(ImageMetadata &imd)
    assert(buffer);      // this is no longer deleted;
    assert(dpTh&&zoNe);  // these are required during finalizeRead()

    ColorMap colorMapTest;
    double range_scale = 0.0;
    if (b_data) {
      //changes to colormap not supported here, handled through initializeWriteRecolor
      assert(colorMapTest.set_type(metadata->getColorMapName()));
      assert(colorMapTest.getName()==theColorMap.getName());
      metadata->setRangeMin(range[0]);
      metadata->setRangeMax(range[1]);
      metadata->setVarId(varId);
      if (range[1]-range[0]==0.0) range_scale = 0.0;
      else range_scale = 1.0/(double(range[1])-double(range[0]));
    }
    double rangeSurface_scale = 0.0;
    bool b_grayscaleColormapSurface = false;
    if (b_dataSurface) {
      //changes to colormap not supported here, handled through initializeWriteRecolor
      assert(colorMapTest.set_type(metadata->getSurfColorMapName()));
      assert(colorMapTest.getName()==theColorMapSurface.getName());
      if (theColorMapSurface.getName() == "GRAYSCALE_RGB" ||
          theColorMapSurface.getName() == "INVERTED_GRAYSCALE_RGB"){
        b_grayscaleColormapSurface = true;
      }

      metadata->setRangeOnSurfaceMin(rangeSurface[0]);
      metadata->setRangeOnSurfaceMax(rangeSurface[1]);
      metadata->setVarOnSurfaceId(varIdSurface);
      if (rangeSurface[1]-rangeSurface[0]==0.0) rangeSurface_scale = 0.0;
      else rangeSurface_scale = 1.0/(double(rangeSurface[1])-double(rangeSurface[0]));
    }
    double rangeParticles_scale = 0.0;
    bool b_grayscaleColormapParticles = false;
    if (b_dataParticles) {
      //changes to colormap not supported here, handled through initializeWriteRecolor
      assert(colorMapTest.set_type(metadata->getPartColorMapName()));
      assert(colorMapTest.getName()==theColorMapParticles.getName());
      if (theColorMapParticles.getName() == "GRAYSCALE_RGB" ||
          theColorMapParticles.getName() == "INVERTED_GRAYSCALE_RGB"){
        b_grayscaleColormapParticles = true;
      }

      metadata->setRangeOnParticleMin(rangeParticles[0]);
      metadata->setRangeOnParticleMax(rangeParticles[1]);
      metadata->setVarOnParticleId(varIdParticles);
      if (rangeParticles[1]-rangeParticles[0]==0.0) rangeParticles_scale = 0.0;
      else rangeParticles_scale = 1.0/(double(rangeParticles[1])-double(rangeParticles[0]));
    }
    double rangeIso_scale = 0.0;
    bool b_grayscaleColormapIso = false;
    if (b_dataIso) {
      //changes to colormap not supported here, handled through initializeWriteRecolor
      assert(colorMapTest.set_type(metadata->getIsoColorMapName()));
      assert(colorMapTest.getName()==theColorMapIso.getName());
      if (theColorMapIso.getName() == "GRAYSCALE_RGB" ||
          theColorMapIso.getName() == "INVERTED_GRAYSCALE_RGB"){
        b_grayscaleColormapIso = true;
      }

      metadata->setRangeOnIsoMin(rangeIso[0]);
      metadata->setRangeOnIsoMax(rangeIso[1]);
      metadata->setVarOnIsoId(varIdIso);
      if (rangeIso[1]-rangeIso[0]==0.0) rangeIso_scale = 0.0;
      else rangeIso_scale = 1.0/(double(rangeIso[1])-double(rangeIso[0]));
    }

    unsigned int lightFail = 0;
    for (int ipx = 0; ipx < nx*ny ; ++ipx) {
      if (pixel_flag[ipx] >= VOLUME_DATA_PIXEL){

        double phi = 0;
        unsigned char rgb[3] = {0,0,0};
        if (pixel_flag[ipx] == VOLUME_DATA_PIXEL){
          assert(b_data);
          phi = (pixel_data[ipx] - range[0])*range_scale;
          //update daTa chunk if present
          if (daTa){
            unsigned char daTa_uc = (unsigned char) min(max(phi*255.0 + 0.5,0.0),255.0);
            daTa->set(ipx,&daTa_uc);
          }
          //lighting factor set to 1 if present
          if (lgFr){
            unsigned char lgFr_uc = (unsigned char) 255;
            lgFr->set(ipx,&lgFr_uc);
          }
          theColorMap.calcColorMapRGB(rgb, min(max(phi,0.0),1.0));
        }
        else {
          assert(daTa); //must be present for non-planar data

          //attempt to reconstruct the lighting for non-planar data...
          //  notes:
          //1. Lighting is applied in CTI canvas as val_lit = (int) min(max(val*p,0),255)
          //   where p is the lighting, val is unlit 0-255 data value
          //2. Some precision will never be recovered here due to double to int
          //   conversion.
          //3. If val_lit>0 or val_lit<255 for this pixel, then the lighting factor p
          //   can be computed to within the lost precision noted in #2
          //4. If val_lit=255 then the lighting factor cannot be accurately computed;
          //   however, if max range new < max range old than this pixel will still be
          //   flooded.
          //5. Similar logic to 4. if val_lit=0 and min range new > min range old
          //6. If 4 or 5 conditions are not met rescaling may look strange (unresolved?)
          //   Use unlit or approximately lit data when exact lighting cannot be
          //   computed.  This may look better than simply setting phi=0 or phi=255
          //   if the new data range is significantly larger.
          
          unsigned char daTa_uc;
          double lightFactor = 1.0;
          if (lgFr){
            unsigned char lgFr_uc;
            lgFr->get(ipx,&lgFr_uc);
            lightFactor = ((int)lgFr_uc)/255.0;
            //TODO update zone buffer with pixel type?  currently assuming it is unchanged...
            if (pixel_flag[ipx] == SURFACE_DATA_PIXEL && b_dataSurface){
              phi = (pixel_data[ipx] - rangeSurface[0])*rangeSurface_scale;
              theColorMapSurface.calcColorMapRGB(rgb, min(max(phi,0.0),1.0)); //use unlit value to calc rgb
              FOR_I3 rgb[i] = min(max(int(int(rgb[i])*lightFactor),0),255);
            }
            else if (pixel_flag[ipx] == PARTICLE_DATA_PIXEL && b_dataParticles){
              phi = (pixel_data[ipx] - rangeParticles[0])*rangeParticles_scale;
              theColorMapParticles.calcColorMapRGB(rgb, min(max(phi,0.0),1.0)); //use unlit value to calc rgb
              FOR_I3 rgb[i] = min(max(int(int(rgb[i])*lightFactor),0),255);
            }
            else if (pixel_flag[ipx] == ISO_DATA_PIXEL && b_dataIso){
              phi = (pixel_data[ipx] - rangeIso[0])*rangeIso_scale;
              theColorMapIso.calcColorMapRGB(rgb, min(max(phi,0.0),1.0));  //use unlit value to calc rgb
              FOR_I3 rgb[i] = min(max(int(int(rgb[i])*lightFactor),0),255);
            }
          }
          else {
            daTa->get(ipx,&daTa_uc);
            double phi0_unlit = ((int)daTa_uc)/255.0;

            //TODO update zone buffer with pixel type?  currently assuming it is unchanged...
            if (pixel_flag[ipx] == SURFACE_DATA_PIXEL && b_dataSurface){
              phi = (pixel_data[ipx] - rangeSurface[0])*rangeSurface_scale;
              theColorMapSurface.calcColorMapRGB(rgb, min(max(phi,0.0),1.0)); //use unlit value to calc rgb

              unsigned char rgb0_lit[3] = {buffer[ipx][0],buffer[ipx][1],buffer[ipx][2]};
              unsigned char rgb0_unlit[3];
              theColorMapSurface.calcColorMapRGB(rgb0_unlit, min(max(phi0_unlit,0.0),1.0));
              if (!reconstructLightFactorWithRescale(lightFactor,rgb0_unlit,rgb0_lit,phi,phi0_unlit,b_grayscaleColormapSurface)){
                ++lightFail;
              }
              FOR_I3 rgb[i] = min(max(int(int(rgb[i])*lightFactor),0),255);
            }
            else if (pixel_flag[ipx] == PARTICLE_DATA_PIXEL && b_dataParticles){
              phi = (pixel_data[ipx] - rangeParticles[0])*rangeParticles_scale;
              theColorMapParticles.calcColorMapRGB(rgb, min(max(phi,0.0),1.0)); //use unlit value to calc rgb

              unsigned char rgb0_lit[3] = {buffer[ipx][0],buffer[ipx][1],buffer[ipx][2]};
              unsigned char rgb0_unlit[3];
              theColorMapParticles.calcColorMapRGB(rgb0_unlit, min(max(phi0_unlit,0.0),1.0));
              if (!reconstructLightFactorWithRescale(lightFactor,rgb0_unlit,rgb0_lit,phi,phi0_unlit,b_grayscaleColormapParticles)){
                ++lightFail;
              }
              FOR_I3 rgb[i] = min(max(int(int(rgb[i])*lightFactor),0),255);
            }
            else if (pixel_flag[ipx] == ISO_DATA_PIXEL && b_dataIso){
              phi = (pixel_data[ipx] - rangeIso[0])*rangeIso_scale;
              theColorMapIso.calcColorMapRGB(rgb, min(max(phi,0.0),1.0));  //use unlit value to calc rgb

              unsigned char rgb0_lit[3] = {buffer[ipx][0],buffer[ipx][1],buffer[ipx][2]};
              unsigned char rgb0_unlit[3];
              theColorMapIso.calcColorMapRGB(rgb0_unlit, min(max(phi0_unlit,0.0),1.0));
              if (!reconstructLightFactorWithRescale(lightFactor,rgb0_unlit,rgb0_lit,phi,phi0_unlit,b_grayscaleColormapIso)){
                ++lightFail;
              }
              FOR_I3 rgb[i] = min(max(int(int(rgb[i])*lightFactor),0),255);
            }
          }
          //set new unlit value to daTa
          daTa_uc = (unsigned char) min(max(phi*255.0 + 0.5,0.0),255.0);
          daTa->set(ipx,&daTa_uc);
        }
        //apply rgb to image buffer
        for (int k =0; k < 3 ; ++k) buffer[ipx][k] = rgb[k];
      }
      else if (pixel_flag[ipx]==FLAG_RED){
        //used by diff, red
        //base image contains this geometry, this image does not
        buffer[ipx][0] = 255;
        buffer[ipx][1] = 0;
        buffer[ipx][2] = 0;
        // also flip the zone to background to be consistent with the
        // pixel_flag changes...
        unsigned char zoNe_uc[2] = {(unsigned char) 255, (unsigned char) 255};
        zoNe->set(ipx,zoNe_uc);
      }
      else if (pixel_flag[ipx]==FLAG_GREEN){
        //used by diff, green
        //base image does not contain this geometry, this image does
        buffer[ipx][0] = 0;
        buffer[ipx][1] = 255;
        buffer[ipx][2] = 0;
        // also flip the zone to background to be consistent with the
        // pixel_flag changes...
        unsigned char zoNe_uc[2] = {(unsigned char) 255, (unsigned char) 255};
        zoNe->set(ipx,zoNe_uc);
      }
      else if (pixel_flag[ipx]==FLAG_BLUE){
        //used by diff, green
        //base image contains a different data type than this image.
        buffer[ipx][0] = 0;
        buffer[ipx][1] = 0;
        buffer[ipx][2] = 255;
        // also flip the zone to background to be consistent with the
        // pixel_flag changes...
        unsigned char zoNe_uc[2] = {(unsigned char) 255, (unsigned char) 255};
        zoNe->set(ipx,zoNe_uc);
      }
      else if (pixel_flag[ipx]==FLAG_GRAY){
        buffer[ipx][0] = 128;
        buffer[ipx][1] = 128;
        buffer[ipx][2] = 128;
        // also flip the zone to background to be consistent with the
        // pixel_flag changes...
        unsigned char zoNe_uc[2] = {(unsigned char) 255, (unsigned char) 255};
        zoNe->set(ipx,zoNe_uc);
      }
      //=======================================================================================//
      //    BEGIN JOB HACK 4/2/18                                                              //
      //    Hard code a bunch of pixel_flag values to correspond to colors for contour flagging//
      //=======================================================================================//
      else if (pixel_flag[ipx] == FLAG_WHITE) {
        //      Color = White
        buffer[ipx][0] = 255;
        buffer[ipx][1] = 255;
        buffer[ipx][2] = 255;
      }//(pixel_flag[ipx] == -10)
      else if (pixel_flag[ipx] == FLAG_RED2) {
        //      Color = Red
        buffer[ipx][0] = 230;
        buffer[ipx][1] = 25;
        buffer[ipx][2] = 75;
      }//(pixel_flag[ipx] == -11)
      else if (pixel_flag[ipx] == FLAG_GREEN2) {
        //      Color = Green
        buffer[ipx][0] = 60;
        buffer[ipx][1] = 180;
        buffer[ipx][2] = 75;
      }//(pixel_flag[ipx] == -12)
      else if (pixel_flag[ipx] == FLAG_YELLOW) {
        //      Color =  Yellow
        buffer[ipx][0] = 255;
        buffer[ipx][1] = 225;
        buffer[ipx][2] = 25;
      }//(pixel_flag[ipx] == -13)
      else if (pixel_flag[ipx] == FLAG_BLUE2) {
        //      Color = Blue
        buffer[ipx][0] = 0;
        buffer[ipx][1] = 130;
        buffer[ipx][2] = 200;
      }//(pixel_flag[ipx] == -14)
      else if (pixel_flag[ipx] == FLAG_ORANGE) {
        //      Color = Orange
        buffer[ipx][0] = 245;
        buffer[ipx][1] = 130;
        buffer[ipx][2] = 48;
      }//(pixel_flag[ipx] == -15)
      else if (pixel_flag[ipx] == FLAG_PURPLE) {
        //      Color = Purple
        buffer[ipx][0] = 145;
        buffer[ipx][1] = 30;
        buffer[ipx][2] = 180;
      }//(pixel_flag[ipx] == -16)
      else if (pixel_flag[ipx] == FLAG_CYAN) {
        //      Color = Cyan
        buffer[ipx][0] = 70;
        buffer[ipx][1] = 240;
        buffer[ipx][2] = 240;
      }//(pixel_flag[ipx] == -17)
      else if (pixel_flag[ipx] == FLAG_MAGENTA) {
        //      Color = Magenta
        buffer[ipx][0] = 240;
        buffer[ipx][1] = 50;
        buffer[ipx][2] = 230;
      }//(pixel_flag[ipx] == -18)
      else if (pixel_flag[ipx] == FLAG_PINK) {
        //      Color = Pink
        buffer[ipx][0] = 250;
        buffer[ipx][1] = 190;
        buffer[ipx][2] = 190;
      }//(pixel_flag[ipx] == -19)
      else if (pixel_flag[ipx] == FLAG_TEAL) {
        //      Color = Teal
        buffer[ipx][0] = 0;
        buffer[ipx][1] = 128;
        buffer[ipx][2] = 128;
      }//(pixel_flag[ipx] == -20)
      //=======================================================================================//
      //    END JOB HACK 4/2/18                                                                //
      //    Hard code a bunch of pixel_flag values to correspond to colors for contour flagging//
      //=======================================================================================//
      else if (pixel_flag[ipx]==FLAG_BLANK){
        // display based on what was put in buffer by previous method
        //TODO What is this for? -DP
      }
    }

    if (theLegend.isInit()) theLegend.addLegend(buffer,zoNe,range[0],range[1],metadata->getVarId());
    if (theLegendSurface.isInit()) theLegendSurface.addLegend(buffer,zoNe,rangeSurface[0],rangeSurface[1],metadata->getSurfVarId());
    if (theLegendParticles.isInit()) theLegendParticles.addLegend(buffer,zoNe,rangeParticles[0],rangeParticles[1],metadata->getPartVarId());
    if (theLegendIso.isInit()) theLegendIso.addLegend(buffer,zoNe,rangeIso[0],rangeIso[1],metadata->getIsoVarId());

    //compress updated daTa chunk for writing to image.
    if (daTa)
      daTa->compressChunk();
    if (lgFr)
      lgFr->compressChunk();
    if (zoNe)
      zoNe->compressChunk();
    if (dpTh)
      dpTh->compressChunk();

    if (lightFail>0){
      COUT1(" > Warning - writing an image with data where the lighting could not be reconstructed for " << lightFail  << " pixels");
    }
  }

  bool reconstructLightFactor(double &lightFactor, const unsigned char rgb0_unlit[3], const unsigned char rgb0_lit[3]){
    int i_rgb = 0;
    for (i_rgb = 0; i_rgb<3; ++i_rgb){
      int val0_unlit = int(rgb0_unlit[i_rgb]);
      int val0_lit   = int(rgb0_lit[i_rgb]);
      if (val0_unlit<=0){
        lightFactor = 0.0; //unable to reconstruct light factor, make this pixels black
      }
      else{
        lightFactor = ((double) val0_lit)/val0_unlit;
        if (val0_lit>0){
          break; //if the values aren't flooded, stop here, otherwise check the next color.
        }
      }
    }
    return (i_rgb<3); //false if could not accurately reconstruct lighting
  }

  bool reconstructLightFactorWithRescale(double &lightFactor, unsigned char rgb0_unlit[3],unsigned char rgb0_lit[3], const double &phi, const double &phi0_unlit, const bool b_grayscaleColormap){
    
    bool b_lightSuccess = reconstructLightFactor(lightFactor,rgb0_unlit,rgb0_lit);
    if (!b_lightSuccess){
      //if it appears the lighting cannot be reconstructed, for a grayscale colormap this still might
      //be okay if the range is contracting during a rescale
      if (b_grayscaleColormap){
        int val0_unlit = int(rgb0_unlit[2]);
        int val0_lit   = int(rgb0_lit[2]);
        if (val0_lit<=0 || val0_unlit<=0){  //val0_lit < val0_unlit
          if (phi < (phi0_unlit+(0.5/255))){
            b_lightSuccess = true;
          }
        }
        else if (val0_lit>=255){
          if (phi >= (phi0_unlit-(0.5/255))){
            b_lightSuccess = true;
          }

        }
      }
    }
    return b_lightSuccess;

  }

  //When recoloring an image, finalizeRead is not called,
  //therefore initializeWriteRecolor must work only with the
  //underlying PngImage and Metadata values, except for
  //PngData Colormaps which must be set by the calling
  //function. Colormaps not set will default to Grayscale.
  //ColorMap..Buf entries contain the prior Colormap.
  void initializeWriteRecolor(bool b_doVolume, bool b_doSurface, bool b_doParticles, bool b_doIso) {
    assert(b_readImage);
    assert(buffer);

    assert (zoNe == NULL);
    zoNe = getPngDataChunk("zoNe");
    if (!zoNe) {
      CERR("Zone information not found");
    }

    //data chunk will be present whenever lit data included
    assert (daTa == NULL);
    daTa = getPngDataChunk("daTa");
    //lighting factor stored with newer images with lit data
    assert (lgFr == NULL);
    lgFr = getPngDataChunk("lgFr"); 

    ColorMap theColorMapBuf;
    ColorMap theColorMapSurfaceBuf;
    ColorMap theColorMapParticlesBuf;
    ColorMap theColorMapIsoBuf;

    if (b_doVolume){ //volume data is not lit; use daTa if present otherwise try to get value from colormap
      if (!daTa){
        if(!theColorMapBuf.set_type(metadata->getColorMapName())){
          COUT1("Unable to recolor volume/planar data, unrecognized input colormap " << metadata->getColorMapName());
          b_doVolume = false;
        }
        else if (!theColorMapBuf.get_b_dataFromRgb()){
          COUT1("Unable to recolor volume/planar data, data values cannot be recomputed from input colormap " << metadata->getColorMapName());
          b_doVolume = false;
        }
      }
    }

    bool b_doSurfaceLight = true;
    if (b_doSurface){
      if (!daTa){
        COUT1("Warning: missing daTa chunk, unable to recolor surface data");
        b_doSurface = false;
      }
      else{
        if(!lgFr&&!theColorMapSurfaceBuf.set_type(metadata->getSurfColorMapName())){
          COUT1("Unable to recolor surface data with lighting, no lgFr chunk and unrecognized input colormap " << metadata->getSurfColorMapName());
          COUT1("Recoloring from daTa chunk with unlit surface data.");
          b_doSurfaceLight = false;
        }
      }

      //zoneId alone is not enough to flag surface data pixels...
      if (metadata->hasText("VAR_ON_SURFACE") && metadata->getSurfVarId() != "OFF" &&
        metadata->hasText("RANGE_ON_SURFACE_MAX") && metadata->hasText("RANGE_ON_SURFACE_MIN") &&
        metadata->hasText("COLORMAP_SURFACE")){
          b_dataSurface = true;
      }

    }
    bool b_doParticlesLight = true;
    if (b_doParticles){
      if (!daTa){
        COUT1("Warning: missing daTa chunk, unable to recolor particle data");
        b_doParticles = false;
      }
      else{
        if(!lgFr&&!theColorMapParticlesBuf.set_type(metadata->getPartColorMapName())){
          COUT1("Unable to recolor particle data with lighting, no lgFr chunk and unrecognized input colormap " << metadata->getPartColorMapName());
          COUT1("Recoloring from daTa chunk with unlit particle data.");
          b_doParticlesLight = false;
        }
      }
    }
    bool b_doIsoLight = true;
    if (b_doIso){
      if (!daTa){
        COUT1("Warning: missing daTa chunk, unable to recolor iso data");
        b_doIso = false;
      }
      else{
        if(!lgFr&&!theColorMapIsoBuf.set_type(metadata->getIsoColorMapName())){
          COUT1("Unable to recolor iso data with lighting, no lgFr chunk and unrecognized input colormap " << metadata->getIsoColorMapName());
          COUT1("Recoloring from daTa chunk with unlit iso data.");
          b_doIsoLight = false;
        }
      }
    }

    unsigned int lightFail = 0;
    for (int ipx = 0; ipx < nx*ny ; ++ipx) {
      unsigned char zoNe_uc[2];
      zoNe->get(ipx,zoNe_uc);
      uint2 zoNe_us = (zoNe_uc[1] << 8 | zoNe_uc[0]); //little endian

      double phi_unlit = 0.0;
      double phi = 0.0;
      unsigned char rgb[3] = {0,0,0};
      if (b_doVolume && zoNe_us == 65534){
        b_data = true;
        if (daTa){
          unsigned char daTa_uc;
          daTa->get(ipx,&daTa_uc);
          phi = ((int)daTa_uc)/255.0;
        }
        else{
          assert(theColorMapBuf.get_b_dataFromRgb());
          phi = theColorMapBuf.calcPhiFromRgb(buffer[ipx]);
        }
        theColorMap.calcColorMapRGB(rgb, min(max(phi,0.0),1.0));
        for (int k =0; k < 3 ; ++k) buffer[ipx][k] = rgb[k]; //update image
      }
      else if (b_doParticles && zoNe_us == 65533){
        b_dataParticles = true;
        assert(daTa);
        unsigned char daTa_uc;
        daTa->get(ipx,&daTa_uc);
        phi_unlit = ((int)daTa_uc)/255.0;
        theColorMapParticles.calcColorMapRGB(rgb, min(max(phi_unlit,0.0),1.0));
        if (b_doParticlesLight){ //true as long as input colormap is known
          //apply lighting
          double lightFactor = 1.0;
          if (lgFr){
            unsigned char lgFr_uc;
            lgFr->get(ipx,&lgFr_uc);
            lightFactor = ((int)lgFr_uc)/255.0;
          }
          else{
            unsigned char rgb0_unlit[3];
            theColorMapParticlesBuf.calcColorMapRGB(rgb0_unlit,min(max(phi_unlit,0.0),1.0));
            if (!reconstructLightFactor(lightFactor,rgb0_unlit,buffer[ipx]))
              ++lightFail;
          }
          FOR_I3 rgb[i] = min(max(int(int(rgb[i])*lightFactor),0),255);
        }
        for (int k =0; k < 3 ; ++k) buffer[ipx][k] = rgb[k]; //update image
      }
      else if (b_doIso && zoNe_us == 65532){
        b_dataIso = true;
        assert(daTa);
        unsigned char daTa_uc;
        daTa->get(ipx,&daTa_uc);
        phi_unlit = ((int)daTa_uc)/255.0;
        theColorMapIso.calcColorMapRGB(rgb, min(max(phi_unlit,0.0),1.0));
        if (b_doIsoLight){ //true as long as input colormap is known
          //apply lighting
          double lightFactor = 1.0;
          if (lgFr){
            unsigned char lgFr_uc;
            lgFr->get(ipx,&lgFr_uc);
            lightFactor = ((int)lgFr_uc)/255.0;
          }
          else {
            unsigned char rgb0_unlit[3];
            theColorMapIsoBuf.calcColorMapRGB(rgb0_unlit,min(max(phi_unlit,0.0),1.0));
            if (!reconstructLightFactor(lightFactor,rgb0_unlit,buffer[ipx]))
              ++lightFail;
          }
          FOR_I3 rgb[i] = min(max(int(int(rgb[i])*lightFactor),0),255);
        }
        for (int k =0; k < 3 ; ++k) buffer[ipx][k] = rgb[k]; //update image
      }
      else if (b_doSurface && b_dataSurface && zoNe_us >= 32768 && zoNe_us < 65532){
        assert(daTa);//checked earlier
        unsigned char daTa_uc;
        daTa->get(ipx,&daTa_uc);
        phi_unlit = ((int)daTa_uc)/255.0;
        theColorMapSurface.calcColorMapRGB(rgb, min(max(phi_unlit,0.0),1.0));
        if (b_doSurfaceLight){ //true as long as input colormap is known (or lgFr exists)
          //apply lighting; rgb_lit < rgb_unlit
          double lightFactor = 1.0;
          if (lgFr){
            unsigned char lgFr_uc;
            lgFr->get(ipx,&lgFr_uc);
            lightFactor = ((int)lgFr_uc)/255.0;
          }
          else{
            unsigned char rgb0_unlit[3];
            theColorMapSurfaceBuf.calcColorMapRGB(rgb0_unlit,min(max(phi_unlit,0.0),1.0));
            if (!reconstructLightFactor(lightFactor,rgb0_unlit,buffer[ipx]))
              ++lightFail;
          }
          FOR_I3 rgb[i] = min(max(int(int(rgb[i])*lightFactor),0),255);
        }
        for (int k =0; k < 3 ; ++k) buffer[ipx][k] = rgb[k]; //update image
      }
    } //end for all pixels

    if (lightFail>0)
      COUT1(" > Warning - writing an image with data where the lighting could not be reconstructed for " << lightFail  << " pixels");


    //Only set output colormap names for data types in this image...
    if (b_data&&b_doVolume)  metadata->setColorMapName(theColorMap.getName());
    if (b_dataParticles&&b_doParticles) metadata->setColorMapParticlesName(theColorMapParticles.getName());
    if (b_dataIso&&b_doIso) metadata->setColorMapIsoName(theColorMapIso.getName());
    if (b_dataSurface&&b_doSurface) metadata->setColorMapSurfaceName(theColorMapSurface.getName());

    if (theLegend.isInit()) theLegend.addLegend(buffer,zoNe,metadata->getRangeMin(),metadata->getRangeMax(),metadata->getVarId());
    if (theLegendSurface.isInit()) theLegendSurface.addLegend(buffer,zoNe,metadata->getSurfRangeMin(),metadata->getSurfRangeMax(),metadata->getSurfVarId());
    if (theLegendParticles.isInit()) theLegendParticles.addLegend(buffer,zoNe,metadata->getPartRangeMin(),metadata->getPartRangeMax(),metadata->getPartVarId());
    if (theLegendIso.isInit()) theLegendIso.addLegend(buffer,zoNe,metadata->getIsoRangeMin(),metadata->getIsoRangeMax(),metadata->getIsoVarId());

    //must explicitly compress chunks for inclusion in image
    assert(!dpTh);
    dpTh = getPngDataChunk("dpTh");

    if (daTa)
      daTa->compressChunk();
    if (lgFr)
      lgFr->compressChunk();
    if (zoNe)
      zoNe->compressChunk();
    if (dpTh)
      dpTh->compressChunk();
  }

  //Prior to calling this routine, user is responsible for:
  // 1. Setting appropriate data flags (b_data, b_dataSurface, b_dataParticles)
  // 2. Setting appropirate data ranges either manually or via rescaleRangesToData()
  // 3. Setting Colormap/Legend if something other than grayscale/Off is desired
  // 4. Setting the image dpTh buffer
  // 5. Setting other ImageMetadata entries (transformMat, var names, CAM_DEPTH, etc)
  void initializeWrite(ImageMetadata &imd) {

    assert(!b_readImage); // only use this for initialized PngData Objects (not read in)
    assert(buffer);      // this is no longer deleted;
    assert(dpTh&&zoNe);  // these are required during finalizeRead()

    metadata = &imd;

    double range_scale = 0.0;
    if (b_data) {
      metadata->setColorMapName(theColorMap.getName());
      metadata->setRangeMin(range[0]);
      metadata->setRangeMax(range[1]);
      if (range[1]-range[0]==0.0) range_scale = 0.0;
      else range_scale = 1.0/(double(range[1])-double(range[0]));
    }
    double rangeSurface_scale = 0.0;
    ColorMap theColorMapSurfaceBuf;
    if (b_dataSurface) {
      theColorMapSurfaceBuf.set_type(metadata->getSurfColorMapName());
      metadata->setColorMapSurfaceName(theColorMapSurface.getName());
      metadata->setRangeOnSurfaceMin(rangeSurface[0]);
      metadata->setRangeOnSurfaceMax(rangeSurface[1]);
      if (rangeSurface[1]-rangeSurface[0]==0.0) rangeSurface_scale = 0.0;
      else rangeSurface_scale = 1.0/(double(rangeSurface[1])-double(rangeSurface[0]));
    }
    double rangeParticles_scale = 0.0;
    ColorMap theColorMapParticlesBuf;
    if (b_dataParticles) {
      theColorMapParticlesBuf.set_type(metadata->getPartColorMapName());
      metadata->setColorMapParticlesName(theColorMapParticles.getName());
      metadata->setRangeOnParticleMin(rangeParticles[0]);
      metadata->setRangeOnParticleMax(rangeParticles[1]);
      if (rangeParticles[1]-rangeParticles[0]==0.0) rangeParticles_scale = 0.0;
      else rangeParticles_scale = 1.0/(double(rangeParticles[1])-double(rangeParticles[0]));
    }
    double rangeIso_scale = 0.0;
    ColorMap theColorMapIsoBuf;
    if (b_dataIso) {
      theColorMapIsoBuf.set_type(metadata->getIsoColorMapName());
      metadata->setColorMapIsoName(theColorMapIso.getName());
      metadata->setRangeOnIsoMin(rangeIso[0]);
      metadata->setRangeOnIsoMax(rangeIso[1]);
      if (rangeIso[1]-rangeIso[0]==0.0) rangeIso_scale = 0.0;
      else rangeIso_scale = 1.0/(double(rangeIso[1])-double(rangeIso[0]));
    }

    for (int ipx = 0; ipx < nx*ny ; ++ipx) {
      if (pixel_flag[ipx] >= VOLUME_DATA_PIXEL){

        unsigned char zoNe_uc[2] = {(unsigned char) 255, (unsigned char) 255};
        uint2 zoneId = 65534;

        double phi = 0;
        unsigned char rgb[3] = {0,0,0};
        if (pixel_flag[ipx] == VOLUME_DATA_PIXEL){
          assert(b_data);
          phi = (pixel_data[ipx] - range[0])*range_scale;
          //update daTa chunk if present
          if (daTa){
            unsigned char daTa_uc = (unsigned char) min(max(phi*255.0 + 0.5,0.0),255.0);
            daTa->set(ipx,&daTa_uc);
          }
          if (lgFr){
            unsigned char lgFr_uc = (unsigned char) 255;
            lgFr->set(ipx,&lgFr_uc);
          }
          theColorMap.calcColorMapRGB(rgb, min(max(phi,0.0),1.0));
          zoneId = 65534;
        }
        else {

          if (pixel_flag[ipx] == SURFACE_DATA_PIXEL){
            assert(b_dataSurface);
            phi = (pixel_data[ipx] - rangeSurface[0])*rangeSurface_scale;
            theColorMapSurface.calcColorMapRGB(rgb, min(max(phi,0.0),1.0));
            zoneId = 32768; //naive surface with data id specification
          }
          else if (pixel_flag[ipx] == PARTICLE_DATA_PIXEL){
            assert(b_dataParticles);
            phi = (pixel_data[ipx] - rangeParticles[0])*rangeParticles_scale;
            theColorMapParticles.calcColorMapRGB(rgb, min(max(phi,0.0),1.0));
            zoneId = 65533;
          }
          else if (pixel_flag[ipx] == ISO_DATA_PIXEL){
            assert(b_dataIso);
            phi = (pixel_data[ipx] - rangeIso[0])*rangeIso_scale;
            theColorMapIso.calcColorMapRGB(rgb, min(max(phi,0.0),1.0));
            zoneId = 65532;
          }
          unsigned char daTa_uc;
          daTa_uc = (unsigned char) min(max(phi*255.0 + 0.5,0.0),255.0);
          daTa->set(ipx,&daTa_uc);
          //TODO allow for real lighting factor values, apply to rgb values of lit data
          if (lgFr){
            unsigned char lgFr_uc = (unsigned char) 255;
            lgFr->set(ipx,&lgFr_uc);
          }



        }
        //apply rgb to image buffer
        for (int k =0; k < 3 ; ++k) buffer[ipx][k] = rgb[k];
        //set zoneId to chunk
        zoNe_uc[0] = zoneId & 255;
        zoNe_uc[1] = zoneId >> 8;
        zoNe->set(ipx,zoNe_uc);

      }
    }

    if (theLegend.isInit()) theLegend.addLegend(buffer,zoNe,range[0],range[1],metadata->getVarId());
    if (theLegendSurface.isInit()) theLegendSurface.addLegend(buffer,zoNe,rangeSurface[0],rangeSurface[1],metadata->getSurfVarId());
    if (theLegendParticles.isInit()) theLegendParticles.addLegend(buffer,zoNe,rangeParticles[0],rangeParticles[1],metadata->getPartVarId());
    if (theLegendIso.isInit()) theLegendIso.addLegend(buffer,zoNe,rangeIso[0],rangeIso[1],metadata->getIsoVarId());

    //compress chunks for writing to image.
    if (daTa)
      daTa->compressChunk();
    if (lgFr)
      lgFr->compressChunk();
    if (zoNe)
      zoNe->compressChunk();
    if (dpTh)
      dpTh->compressChunk();
    

  }

  //pixel_flag set to FLAG_GREEN for contiguous pixels
  void flagContiguousPixels(const int &pixel){
     //Flag all pixels enclosed by image boundaries and no-data region
     //...of the original pixel's pixel type
     list<int> stack;
     int pixel_type = pixel_flag[pixel];
     if (pixel_type == VOLUME_DATA_PIXEL ||
         pixel_type == SURFACE_DATA_PIXEL ||
         pixel_type == PARTICLE_DATA_PIXEL ||
         pixel_type == ISO_DATA_PIXEL){
       stack.push_back(pixel);
     }

     while (stack.size() > 0){
       int current_pixel = stack.back();
       stack.pop_back();
       if (current_pixel%nx!=0){ //not in left column
         int new_pixel = current_pixel-1;
         if (pixel_flag[new_pixel] == pixel_type){
           pixel_flag[new_pixel] = FLAG_GREEN;
           stack.push_back(new_pixel);
         }
       }
       if ((current_pixel+1)%nx!=0){ //not in right column
         int new_pixel = current_pixel+1;
         if (pixel_flag[new_pixel] == pixel_type){
           pixel_flag[new_pixel] = FLAG_GREEN;
           stack.push_back(new_pixel);
         }
       }
       if (current_pixel >= nx){ //not in top row
         int new_pixel = current_pixel - nx;
         if (pixel_flag[new_pixel] == pixel_type){
           pixel_flag[new_pixel] = FLAG_GREEN;
           stack.push_back(new_pixel);
         }
       }
       if (current_pixel < nx*(ny-1)){ //not in bottom row
         int new_pixel = current_pixel + nx;
         if (pixel_flag[new_pixel] == pixel_type){
           pixel_flag[new_pixel] = FLAG_GREEN;
           stack.push_back(new_pixel);
         }
       }
     }

     for (int ipx = 0; ipx < nx*ny ; ++ipx) {
       if ( pixel_flag[ipx] == pixel_type ){
         pixel_flag[ipx] = FLAG_RED;  //unflagged geometry will be written out red
       }
     }
  }

  //set pixel_flag to FLAG_RED
  void flagYLessThanZero(PixelFlagCases _pflag){
    for (int ipx = 0; ipx < nx*ny ; ++ipx) {
      if (pixel_flag[ipx] == _pflag){
        //DAPDAP
        const int iP = ipx%nx;
        const int jP = (int)ipx/nx;
        double xp[3];
        convertPxtoXp(iP,jP,xp);
        if (xp[1]<=2*getLengthScale()){
          pixel_flag[ipx] = FLAG_RED;
        }
      }
    }
  }

  //pixel_flag set to FLAG_GREEN for contiguous pixels
  void flagContiguousZoneIdPixels(int istart, int jstart, PixelFlagCases _pflag){
    //Flag all continguous pixels of type _pflag
    if (istart>=nx||jstart>=ny){
        CERR("Starting pixel (" << istart << "," << jstart << ") bigger than image size (" << nx << "," << ny << ")");
    }
    if (jstart<0){
      if (istart<0){
        CERR("istart and jstart cannot both be less than zero");
      }
      for (int jp = 0; jp < ny; ++jp) {
        int ipx = jp*nx + istart;
        if (pixel_flag[ipx]==_pflag){
          jstart = jp;
          break;
        }
      }
    }
    else if (istart<0){
      if (jstart<0){
        CERR("istart and jstart cannot both be less than zero");
      }
      for (int ip = 0; ip < nx; ++ip) {
        int ipx = jstart*nx + ip;
        if (pixel_flag[ipx]==_pflag){
          istart = ip;
          break;
        }
      }
    }

    list<int> stack;
    int pixel = jstart*nx + istart;
    int pixel_type = pixel_flag[pixel];
    if (pixel_type == _pflag){
      stack.push_back(pixel);
    }

    while (stack.size() > 0){
      int current_pixel = stack.back();
      stack.pop_back();
      if (current_pixel%nx!=0){ //not in left column
        int new_pixel = current_pixel-1;
        if (pixel_flag[new_pixel] == _pflag){
          pixel_flag[new_pixel] = FLAG_GREEN;
          stack.push_back(new_pixel);
        }
      }
      if ((current_pixel+1)%nx!=0){ //not in right column
        int new_pixel = current_pixel+1;
        if (pixel_flag[new_pixel] == _pflag){
          pixel_flag[new_pixel] = FLAG_GREEN;
          stack.push_back(new_pixel);
        }
      }
      if (current_pixel >= nx){ //not in top row
        int new_pixel = current_pixel - nx;
        if (pixel_flag[new_pixel] == _pflag){
          pixel_flag[new_pixel] = FLAG_GREEN;
          stack.push_back(new_pixel);
        }
      }
      if (current_pixel < nx*(ny-1)){ //not in bottom row
        int new_pixel = current_pixel + nx;
        if (pixel_flag[new_pixel] == _pflag){
          pixel_flag[new_pixel] = FLAG_GREEN;
          stack.push_back(new_pixel);
        }
      }
    }

    for (int ipx = 0; ipx < nx*ny ; ++ipx) {
      if ( pixel_flag[ipx] != FLAG_GREEN ){
        pixel_flag[ipx] = FLAG_RED;  //unflagged geometry will be written out red
        pixel_data[ipx] = 0.0;
      }
      else{
        pixel_data[ipx] = HUGE_VAL;
      }
    }
  }


  //Pixels of type _outflag completely enclose regions of pixel type
  //_inflag; Change all _inflag pixels enclosed by _outflag to _outflag.
  void flagEnclosedPixels(PixelFlagCases _inflag, PixelFlagCases _outflag){

    for (int ip = 0; ip < nx ; ++ip) {
      for (int jp = 0; jp < ny ; ++jp) {
        int ipx = jp*nx + ip;
        if ( pixel_flag[ipx] == _outflag ){ //found an _outflag pixel
        //check if this _outflag pixel borders and inflag pixel
          int ip0 = max(0,ip-1);
          int jp0 = max(0,jp-1);
          int ip1 = min(nx-1,ip+1);
          int jp1 = min(ny-1,jp+1);
        
          int ipx0 = jp*nx+ip0; 
          int ipx1 = jp*nx+ip1;
          int jpx0 = jp0*nx+ip; 
          int jpx1 = jp1*nx+ip; 

          list<int> stack;

          if (pixel_flag[ipx0]==_inflag)
            stack.push_back(ipx0);
          if (pixel_flag[ipx1]==_inflag)
            stack.push_back(ipx1);
          if (pixel_flag[jpx0]==_inflag)
            stack.push_back(jpx0);
          if (pixel_flag[jpx1]==_inflag)
            stack.push_back(jpx1);

          while (stack.size() > 0){
            int current_pixel = stack.back();
            stack.pop_back();
            if (current_pixel%nx!=0){ //not in left column
              int new_pixel = current_pixel-1;
              if (pixel_flag[new_pixel] == _inflag){
                pixel_flag[new_pixel] = _outflag;
                stack.push_back(new_pixel);
              }
            }
            if ((current_pixel+1)%nx!=0){ //not in right column
              int new_pixel = current_pixel+1;
              if (pixel_flag[new_pixel] == _inflag){
                pixel_flag[new_pixel] = _outflag;
                stack.push_back(new_pixel);
              }
            }
            if (current_pixel >= nx){ //not in top row
              int new_pixel = current_pixel - nx;
              if (pixel_flag[new_pixel] == _inflag){
                pixel_flag[new_pixel] = _outflag;
                stack.push_back(new_pixel);
              }
            }
            if (current_pixel < nx*(ny-1)){ //not in bottom row
              int new_pixel = current_pixel + nx;
              if (pixel_flag[new_pixel] == _inflag){
                pixel_flag[new_pixel] = _outflag;
                stack.push_back(new_pixel);
              }
            }
          }
        }
      }
    }
  }


  //Assumes a routine like flagContiguousZoneIdPixels() has been called
  //which sets pixel_flag/_data to FLAG_RED/0.0 or FLAG_GREEN/HUGE_VAL.
  //A fast marching method is used to solve the Eikonal equation for distance
  //from red pixels in all green pixels.
  void computeDistanceFromFlaggedPixels(){
    //BUILD AND INITIALIZE GHOST GRID
    int nx_g = nx+2;
    int ny_g = ny+2;
    double * pixel_data_g = new double[nx_g*ny_g];
    int    * pixel_flag_g = new int[nx_g*ny_g];
    for (int ip = 0; ip < nx; ++ip) {
      for (int jp = 0; jp < ny; ++jp) {
        int ipx = jp*nx + ip;
        int ip_g = ip+1;
        int jp_g = jp+1;
        int ipx_g = jp_g*nx_g + ip_g; 
        pixel_data_g[ipx_g] = pixel_data[ipx];
        pixel_flag_g[ipx_g] = pixel_flag[ipx];
      }
    }
    for (int ip_g = 1; ip_g < nx_g-1; ++ip_g) {
      int jp_g = 0;
      int ipx_g = jp_g*nx_g + ip_g; 
      int ipx_g0 = (jp_g+1)*nx_g + ip_g;
      pixel_data_g[ipx_g] = pixel_data_g[ipx_g0];
      pixel_flag_g[ipx_g] = pixel_flag_g[ipx_g0];
      if (pixel_flag_g[ipx_g]==FLAG_GREEN)
        pixel_flag_g[ipx_g]=VOLUME_DATA_PIXEL;
      //pixel_flag_g[ipx_g] = FLAG_RED;
      //pixel_data_g[ipx_g] = 0.0;

      jp_g = ny_g-1;
      ipx_g = jp_g*nx_g + ip_g; 
      ipx_g0 = (jp_g-1)*nx_g + ip_g;
      pixel_data_g[ipx_g] = pixel_data_g[ipx_g0];
      pixel_flag_g[ipx_g] = pixel_flag_g[ipx_g0];
      if (pixel_flag_g[ipx_g]==FLAG_GREEN)
        pixel_flag_g[ipx_g]=VOLUME_DATA_PIXEL;
      //pixel_flag_g[ipx_g] = FLAG_RED;
      //pixel_data_g[ipx_g] = 0.0;
    }
    for (int jp_g = 0; jp_g < ny_g; ++jp_g) {
      int ip_g = 0;
      int ipx_g = jp_g*nx_g + ip_g; 
      int ipx_g0 = jp_g*nx_g + ip_g+1;
      pixel_data_g[ipx_g] = pixel_data_g[ipx_g0];
      pixel_flag_g[ipx_g] = pixel_flag_g[ipx_g0];
      if (pixel_flag_g[ipx_g]==FLAG_GREEN)
        pixel_flag_g[ipx_g]=VOLUME_DATA_PIXEL;
      //pixel_flag_g[ipx_g] = FLAG_RED;
      //pixel_data_g[ipx_g] = 0.0;

      ip_g = nx_g-1;
      ipx_g = jp_g*nx_g + ip_g; 
      ipx_g0 = jp_g*nx_g + ip_g-1;
      pixel_data_g[ipx_g] = pixel_data_g[ipx_g0];
      pixel_flag_g[ipx_g] = pixel_flag_g[ipx_g0];
      if (pixel_flag_g[ipx_g]==FLAG_GREEN)
        pixel_flag_g[ipx_g]=VOLUME_DATA_PIXEL;
      //pixel_flag_g[ipx_g] = FLAG_RED;
      //pixel_data_g[ipx_g] = 0.0;
    }

    //IF RED PIXEL HAS GREEN NEIGHBOR, FLAG GRAY
    //(Flag interface between alive/far as close)

    //INTERNAL PIXELS, don't flag boundary pixels GRAY
    for (int ip_g = 1; ip_g < nx_g-1 ; ++ip_g) {
      for (int jp_g = 1; jp_g < ny_g-1 ; ++jp_g) {
        int ipx_g = jp_g*nx_g + ip_g;
        if (pixel_flag_g[ipx_g]==FLAG_RED){
          int ipx0_g = jp_g*nx_g + (ip_g-1); 
          int ipx1_g = jp_g*nx_g + (ip_g+1); 
          int ipx2_g = (jp_g-1)*nx_g + ip_g; 
          int ipx3_g = (jp_g+1)*nx_g + ip_g; 
          if (pixel_flag_g[ipx0_g]==FLAG_GREEN&&ip_g>1)
            pixel_flag_g[ipx0_g]=FLAG_GRAY;
          if (pixel_flag_g[ipx1_g]==FLAG_GREEN&&ip_g<(nx_g-1))
            pixel_flag_g[ipx1_g]=FLAG_GRAY;
          if (pixel_flag_g[ipx2_g]==FLAG_GREEN&&jp_g>1)
            pixel_flag_g[ipx2_g]=FLAG_GRAY;
          if (pixel_flag_g[ipx3_g]==FLAG_GREEN&&jp_g<(ny_g-1))
            pixel_flag_g[ipx3_g]=FLAG_GRAY;
        }
      }
    }
    //END
    //IF RED PIXEL HAS GREEN NEIGHBOR, FLAG GRAY
    //(Flag interface between alive/far as close)


    //FLAG_RED Alive 
    //FLAG_GRAY Close
    //FLAG_GREEN Far
    //

    //compute initial value, store in pixel_data_g[ipx_g]
    //for GRAY Close cells
    //const double delta = 1;
    //const double rhs   = 1;
    MinHeap trialHeap(nx_g*ny_g);
    for (int ip_g = 1; ip_g < nx_g-1 ; ++ip_g) {
      for (int jp_g = 1; jp_g < ny_g-1 ; ++jp_g) {
        int ipx_g = jp_g*nx_g + ip_g;
        if (pixel_flag_g[ipx_g]==FLAG_GRAY){
          trialHeap.insert(pair<double,int>(HUGE_VAL,ipx_g));
          updateFmmPt(ipx_g,pixel_flag_g,pixel_data_g,trialHeap);
        }
      }
    }

    //FLAG_RED Alive (Initial condition) VOLUME_DATA_PIXEL
    //FLAG_GRAY Close
    //FLAG_GREEN Far
    //cout << "\nSolving Eikonal Equation...";
    //int iter =0;
    while (trialHeap.size() > 0){
      //++iter;
      //if (iter%10000==0){
      //  cout << ".";
      //}
      //cout << "Trial Heap Size Before " << trialHeap.size();
      // get trial point from close vector...
      pair<double,int> trial = trialHeap.extractMin();
      const int trial_ipx_g = trial.second;
      trialHeap.deleteMin(); // remove from trial
      //cout << " After " << trialHeap.size() << endl;
      if (pixel_flag_g[trial_ipx_g] != FLAG_GRAY){
        cout << "ip_g " << trial_ipx_g%nx_g << " jp_g " << (int)trial_ipx_g/nx_g << endl;
        cout << "pixel_flag_g[trial_ipx_g] = " << pixel_flag_g[trial_ipx_g] << endl;
        //trialHeap.displayHeap();
      }
      assert(pixel_flag_g[trial_ipx_g] == FLAG_GRAY);
      pixel_flag_g[trial_ipx_g] = VOLUME_DATA_PIXEL; // tag as alive
      // mark trial_ipx_g neighbors that aren't alive as close...
      int ip_g = trial_ipx_g%nx_g;
      int jp_g = (int)trial_ipx_g/nx_g;

      assert(ip_g>0&&jp_g>0&&ip_g<(nx_g-1)&&jp_g<(ny_g-1)); //no ghost points in the heap

      int i0px_g = jp_g*nx_g + (ip_g-1); 
      int i1px_g = jp_g*nx_g + (ip_g+1); 
      int j0px_g = (jp_g-1)*nx_g + ip_g; 
      int j1px_g = (jp_g+1)*nx_g + ip_g; 
      if (pixel_flag_g[i0px_g]==FLAG_GREEN){
        pixel_flag_g[i0px_g]=FLAG_GRAY;
        trialHeap.insert(pair<double,int>(pixel_data_g[i0px_g],i0px_g)); // push onto heap
      }
      if (pixel_flag_g[i1px_g]==FLAG_GREEN){
        pixel_flag_g[i1px_g]=FLAG_GRAY;
        trialHeap.insert(pair<double,int>(pixel_data_g[i1px_g],i1px_g)); // push onto heap
      }
      if (pixel_flag_g[j0px_g]==FLAG_GREEN){
        pixel_flag_g[j0px_g]=FLAG_GRAY;
        trialHeap.insert(pair<double,int>(pixel_data_g[j0px_g],j0px_g)); // push onto heap
      }
      if (pixel_flag_g[j1px_g]==FLAG_GREEN){
        pixel_flag_g[j1px_g]=FLAG_GRAY;
        trialHeap.insert(pair<double,int>(pixel_data_g[j1px_g],j1px_g)); // push onto heap
      }

      //(Re)-compute for all "close" (GRAY) neighbors above
      if (pixel_flag_g[i0px_g]==FLAG_GRAY){
        updateFmmPt(i0px_g,pixel_flag_g,pixel_data_g,trialHeap);
      }
      if (pixel_flag_g[i1px_g]==FLAG_GRAY){
        updateFmmPt(i1px_g,pixel_flag_g,pixel_data_g,trialHeap);
      }
      if (pixel_flag_g[j0px_g]==FLAG_GRAY){
        updateFmmPt(j0px_g,pixel_flag_g,pixel_data_g,trialHeap);
      }
      if (pixel_flag_g[j1px_g]==FLAG_GRAY){
        updateFmmPt(j1px_g,pixel_flag_g,pixel_data_g,trialHeap);
      }
       
    }
    //cout << "done" << endl;

    //COPY BACK TO STANDARD GRID, CLEANUP GHOST GRID
    for (int ip = 0; ip < nx; ++ip) {
      for (int jp = 0; jp < ny; ++jp) {
        int ipx = jp*nx + ip;
        int ip_g = ip+1;
        int jp_g = jp+1;
        int ipx_g = jp_g*nx_g + ip_g; 
        pixel_data[ipx] = pixel_data_g[ipx_g];
        pixel_flag[ipx] = pixel_flag_g[ipx_g];
        //TODO this is only valid if the image view and projection plane match, otherwise
        //the depth value will need to be estimate from linear variation in surface pixels
        //that are present.
        //change dpTh of new data pixels to 0 for possible coordinate export
        if (pixel_flag[ipx]==VOLUME_DATA_PIXEL&&getDepth(ipx)==65535){
          setDepth(ipx,0);
        }
      }
    }
    delete[] pixel_data_g;
    delete[] pixel_flag_g;

    b_data = true;
    metadata->setColorMapName("GRAYSCALE_RGB");
    theColorMap.set_type(metadata->getColorMapName());

  }

  void updateFmmPt(const int ipx_g,int * pixel_flag_g,double * pixel_data_g,MinHeap &trialHeap){
    assert(pixel_flag_g[ipx_g]==FLAG_GRAY); //working with a close point

    const double delta = 1;
    const double rhs   = 1;

    int nx_g = nx+2;
    int ip_g = ipx_g%nx_g; 
    int jp_g = (int)ipx_g/nx_g;
    int i0px_g = jp_g*nx_g + (ip_g-1); 
    int i1px_g = jp_g*nx_g + (ip_g+1); 
    int j0px_g = (jp_g-1)*nx_g + ip_g; 
    int j1px_g = (jp_g+1)*nx_g + ip_g; 

    double ival = -HUGE_VAL;
    double jval = -HUGE_VAL;
    int imin=-1;
    int jmin=-1;
    if (pixel_flag_g[i0px_g]==FLAG_RED||pixel_flag_g[i0px_g]==VOLUME_DATA_PIXEL){
      ival = pixel_data_g[i0px_g];
      imin = i0px_g;
    }
    if (pixel_flag_g[i1px_g]==FLAG_RED||pixel_flag_g[i1px_g]==VOLUME_DATA_PIXEL){
      if (imin<0||pixel_data_g[i1px_g]<ival){
        ival = pixel_data_g[i1px_g];
        imin = i1px_g;
      }
    }
    if (pixel_flag_g[j0px_g]==FLAG_RED||pixel_flag_g[j0px_g]==VOLUME_DATA_PIXEL){
      jval = pixel_data_g[j0px_g];
      jmin = j0px_g;
    }
    if (pixel_flag_g[j1px_g]==FLAG_RED||pixel_flag_g[j1px_g]==VOLUME_DATA_PIXEL){
      if (jmin<0||pixel_data_g[j1px_g]<jval){
        jval = pixel_data_g[j1px_g];
        jmin = j1px_g;
      }
    }

    assert(imin>=0||jmin>=0);
    if (imin>=0&&jmin>=0){
      //solve quadratic with two neighbors
      double aa = 2.0/(delta*delta); 
      double bb = (-2.0/(delta*delta))*(ival+jval); 
      double cc = (1.0/(delta*delta))*(ival*ival+jval*jval) - rhs; 
      pixel_data_g[ipx_g] = (-bb + sqrt(bb*bb-4.0*aa*cc) ) / (2.0*aa);
    }
    else if (imin>=0||jmin>=0){
      //solve quadratic with one neighbor
      double data_val;
      if (imin>=0)
        data_val = ival;
      else
        data_val = jval;

      pixel_data_g[ipx_g] = data_val + delta*sqrt(rhs);
    }
    else{
      assert(0);
    }
    trialHeap.updateDouble(trialHeap.index(ipx_g),pixel_data_g[ipx_g]);

  }

  void flagZoneGroup(vector<uint2> &_zoneIds, PixelFlagCases _pflag){

    //build list of zoneIds;
    vector<uint2> selectedZoneIds;
    vector<uint2> rejectedZoneIds;
    vector<string> selectedZoneNames;
    for (int iz=0; iz<metadata->zoneIds.size(); ++iz){
      vector<uint2>::iterator it = find(_zoneIds.begin(), _zoneIds.end(), metadata->zoneIds[iz]);
      if (it != _zoneIds.end()){
        selectedZoneIds.push_back(metadata->zoneIds[iz]);
        selectedZoneNames.push_back(metadata->zoneNames[iz]);
      }
    }
    if (selectedZoneIds.size()<_zoneIds.size()){
      cout << "Warning, unable to match some zone ids" << endl;
      cout << "Zones successfully flagged as single group: " << endl;
      for (int iz=0; iz<selectedZoneIds.size(); ++iz){
        cout << "  " << selectedZoneIds[iz] << " " << selectedZoneNames[iz] << endl;
      }

    }
    
    if (selectedZoneIds.size()>0){
      //cout << "Flagging Zones as single group: " << endl;
      //for (int iz=0; iz<selectedZoneIds.size(); ++iz){
      //  cout << "  " << selectedZoneIds[iz] << " " << selectedZoneNames[iz] << endl;
      //}
      
      for (int ip = 0; ip < nx; ++ip) {
        for (int jp = 0; jp < ny; ++jp) {
          int ipx = jp*nx + ip;
          unsigned char zoNe_uc[2];
          zoNe->get(ipx,zoNe_uc);
          uint2 zoNe_us = (zoNe_uc[1] << 8 | zoNe_uc[0]); //little endian
          if (find(selectedZoneIds.begin(), selectedZoneIds.end(), zoNe_us)!=selectedZoneIds.end()){
            pixel_flag[ipx] = _pflag;
          }
        }
      }
    }
  }


  //Assuming we have three types of pixels:
  //VOLUME_DATA_PIXELS with eikonal distance value
  //FLAG_GREEN and FLAG_BLUE pixels from two zone groups
  //these two zone groups were previously all FLAG_RED during the eikonal solve
  void diffuseZoneGroups(){

    //BUILD AND INITIALIZE GHOST GRID
    int nx_g = nx+2;
    int ny_g = ny+2;
    double * pixel_data_g = new double[nx_g*ny_g];
    int    * pixel_flag_g = new int[nx_g*ny_g];
    for (int ip = 0; ip < nx; ++ip) {
      for (int jp = 0; jp < ny; ++jp) {
        int ipx = jp*nx + ip;
        int ip_g = ip+1;
        int jp_g = jp+1;
        int ipx_g = jp_g*nx_g + ip_g; 
        pixel_data_g[ipx_g] = pixel_data[ipx];
        pixel_flag_g[ipx_g] = pixel_flag[ipx];
      }
    }
    for (int ip_g = 1; ip_g < nx_g-1; ++ip_g) { 
      int jp_g = 0;                             //top row ghosts, no corners
      int ipx_g = jp_g*nx_g + ip_g; 
      int ipx_g0 = (jp_g+1)*nx_g + ip_g;
      pixel_data_g[ipx_g] = pixel_data_g[ipx_g0];
      pixel_flag_g[ipx_g] = pixel_flag_g[ipx_g0];

      jp_g = ny_g-1;                            //bottom row ghosts, no corners
      ipx_g = jp_g*nx_g + ip_g; 
      ipx_g0 = (jp_g-1)*nx_g + ip_g;
      pixel_data_g[ipx_g] = pixel_data_g[ipx_g0];
      pixel_flag_g[ipx_g] = pixel_flag_g[ipx_g0];
    }
    for (int jp_g = 0; jp_g < ny_g; ++jp_g) {
      int ip_g = 0;                              //left column ghosts
      int ipx_g = jp_g*nx_g + ip_g; 
      int ipx_g0 = jp_g*nx_g + ip_g+1;
      pixel_data_g[ipx_g] = pixel_data_g[ipx_g0];
      pixel_flag_g[ipx_g] = pixel_flag_g[ipx_g0];

      ip_g = nx_g-1;                             //right column ghosts
      ipx_g = jp_g*nx_g + ip_g; 
      ipx_g0 = jp_g*nx_g + ip_g-1;
      pixel_data_g[ipx_g] = pixel_data_g[ipx_g0];
      pixel_flag_g[ipx_g] = pixel_flag_g[ipx_g0];
    }

    MinHeap trialHeap(nx_g*ny_g);
    for (int ip_g = 1; ip_g < nx_g-1 ; ++ip_g) {
      for (int jp_g = 1; jp_g < ny_g-1 ; ++jp_g) {
        int ipx_g = jp_g*nx_g + ip_g;
        if (pixel_flag_g[ipx_g]==VOLUME_DATA_PIXEL){
          trialHeap.insert(pair<double,int>(pixel_data_g[ipx_g],ipx_g));
        }
      }
    }

    //cout << "\nDiffusing Zone Groups ...";
    while (trialHeap.size() > 0){
      pair<double,int> trial = trialHeap.extractMin();
      const int trial_ipx_g = trial.second;
      trialHeap.deleteMin(); // remove from trial
      assert(pixel_flag_g[trial_ipx_g]==VOLUME_DATA_PIXEL);

      bool flag_group[2] = {false,false};

      int ip_g = trial_ipx_g%nx_g;
      int jp_g = (int)trial_ipx_g/nx_g;
      assert(ip_g>0&&jp_g>0&&ip_g<(nx_g-1)&&jp_g<(ny_g-1)); //no ghost points in the heap

      int i0px_g = jp_g*nx_g + (ip_g-1); 
      int i1px_g = jp_g*nx_g + (ip_g+1); 
      int j0px_g = (jp_g-1)*nx_g + ip_g; 
      int j1px_g = (jp_g+1)*nx_g + ip_g; 
      if (pixel_flag_g[i0px_g]==FLAG_GREEN||pixel_flag_g[i1px_g]==FLAG_GREEN||pixel_flag_g[j0px_g]==FLAG_GREEN||pixel_flag_g[j1px_g]==FLAG_GREEN){
        flag_group[0] = true;
      }
      if (pixel_flag_g[i0px_g]==FLAG_BLUE||pixel_flag_g[i1px_g]==FLAG_BLUE||pixel_flag_g[j0px_g]==FLAG_BLUE||pixel_flag_g[j1px_g]==FLAG_BLUE){
        flag_group[1] = true;
      }

      if (flag_group[0]&&flag_group[1]){
        //do nothing...this is the medial axis...
      }
      else if (flag_group[0]){
        pixel_flag_g[trial_ipx_g] = FLAG_GREEN;
      }
      else if (flag_group[1]){
        pixel_flag_g[trial_ipx_g] = FLAG_BLUE;
      }
      else{
        //pixel_flag_g[trial_ipx_g] = FLAG_RED;
        //assert(0);
      }
    }
    //cout << "done" << endl;

    //COPY BACK TO STANDARD GRID, CLEANUP GHOST GRID
    for (int ip = 0; ip < nx; ++ip) {
      for (int jp = 0; jp < ny; ++jp) {
        int ipx = jp*nx + ip;
        int ip_g = ip+1;
        int jp_g = jp+1;
        int ipx_g = jp_g*nx_g + ip_g;
        //pixel_flag[ipx] = pixel_flag_g[ipx_g];
        if (pixel_flag_g[ipx_g]==VOLUME_DATA_PIXEL){
          pixel_flag[ipx] = FLAG_YELLOW;
        }
        else if (pixel_flag_g[ipx_g]==FLAG_RED){
          pixel_flag[ipx] = FLAG_RED;
        }
      }
    }
    delete[] pixel_data_g;
    delete[] pixel_flag_g;

    b_data = true;
    metadata->setColorMapName("GRAYSCALE_RGB");
    theColorMap.set_type(metadata->getColorMapName());

  }

  //Walk pixels of type _pFlag, write out information
  void extractLineSegment(vector<int> &pixelVec, PixelFlagCases _pflag,const int istart, const int jstart, const bool b_rewind){
    int il = -1;
    int jl = -1;
    if (jstart<0){
      assert(istart>=0&&istart<nx);
      il = istart;
      for (int jp = 0; jp < ny; ++jp) {
        int ipx = jp*nx + istart;
        if (pixel_flag[ipx]==_pflag){
          jl = jp;
          break;
        }
      }
    }
    else if (istart<0){
      assert(jstart>=0&&jstart<nx);
      jl = jstart;
      for (int ip = 0; ip < nx; ++ip) {
        int ipx = jstart*nx + ip;
        if (pixel_flag[ipx]==_pflag){
          il = ip;
          break;
        }
      }
    }
    else{
      il = istart;
      jl = jstart;
      if (il>=nx||jl>=ny){
        CERR("Starting pixel (" << il << "," << jl << ") bigger than image size (" << nx << "," << ny << ")");
      }
    }

    if (il<0||jl<0){
      CERR("Unable to find a starting pixel for extractLineSegment()");
    }

    bool b_next = true;
    if (b_rewind){ //rewind to beginning of line segment before extracting
      vector<int> pixelVecLocal;
      while (b_next){
        b_next = false;
        const int ipx = jl*nx + il;
        pixel_flag[ipx] = FLAG_ORANGE;
        pixelVecLocal.push_back(ipx);
        int i0 = max(0,il-1);
        int i1 = min(nx-1,il+1);
        int j0 = max(0,jl-1);
        int j1 = min(ny-1,jl+1);

        int ipx_test[8];
        
        ipx_test[0] = j0*nx+i1;
        ipx_test[1] = jl*nx+i1;
        ipx_test[2] = j1*nx+i1;
        ipx_test[3] = j0*nx+il;
        ipx_test[4] = j1*nx+il;
        ipx_test[5] = j0*nx+i0;
        ipx_test[6] = jl*nx+i0;
        ipx_test[7] = j1*nx+i0;

        for (int itest = 0; itest<8; ++itest){
          if (pixel_flag[ipx_test[itest]]==_pflag){
            b_next = true;
            il = ipx_test[itest]%nx;
            jl = ipx_test[itest]/nx;
            break;
          }
        }
      }
      for (int ipV = 0; ipV<pixelVecLocal.size(); ++ipV){
        pixel_flag[pixelVecLocal[ipV]] = _pflag;
      }
      pixelVecLocal.clear();
      b_next = true;
    }

    int count = 0;
    while (b_next){
      ++count;
      b_next = false;
      const int ipx = jl*nx + il;
      pixel_flag[ipx] = FLAG_ORANGE;
      pixelVec.push_back(ipx);
      int i0 = max(0,il-1);
      int i1 = min(nx-1,il+1);
      int j0 = max(0,jl-1);
      int j1 = min(ny-1,jl+1);

      int ipx_test[8];
      
      ipx_test[0] = j0*nx+i1;
      ipx_test[1] = jl*nx+i1;
      ipx_test[2] = j1*nx+i1;
      ipx_test[3] = j0*nx+il;
      ipx_test[4] = j1*nx+il;
      ipx_test[5] = j0*nx+i0;
      ipx_test[6] = jl*nx+i0;
      ipx_test[7] = j1*nx+i0;

      for (int itest = 0; itest<8; ++itest){
        if (pixel_flag[ipx_test[itest]]==_pflag){
          b_next = true;
          il = ipx_test[itest]%nx;
          jl = ipx_test[itest]/nx;
          break;
        }
      }
    }
    cout << "extractLineSegment() found " << count << " pixel(s)" << endl;
  }


  
  int computeStats(double stats[8], int pixel_type) {

    double vmin = HUGE_VAL;
    double vmax = -HUGE_VAL;
    double mean = 0.0;
    double rms = 0.0;
    double x_centroid[3] = { 0.0, 0.0, 0.0 };
    double pixel_area = getLengthScale();
    pixel_area *= pixel_area;
    double area = 0.0;
    int count = 0;
    double xP[3];
    for (int ipx =0; ipx < nx*ny; ipx++) {
      if (pixel_flag[ipx] == pixel_type) {
        const double val = pixel_data[ipx];
        if (val<vmin) vmin=val;
        if (val>vmax) vmax=val;
        mean += val;
        rms += val*val;
        const int iP = ipx%nx;
        const int jP = (int)ipx/nx;
        convertPxtoXp(iP,jP,xP);
        FOR_I3 x_centroid[i] += xP[i];
        ++count;
        area += pixel_area;
      }
    }

    if (count > 0) {
      mean /= double(count);
      rms = sqrt(max(0.0,rms/double(count)-mean*mean));
      FOR_I3 x_centroid[i] /= double(count);
    }

    stats[0] = x_centroid[0];
    stats[1] = x_centroid[1];
    stats[2] = x_centroid[2];
    stats[3] = vmin;
    stats[4] = vmax;
    stats[5] = mean;
    stats[6] = rms;
    stats[7] = area;
    
    return count; // pixel count

  }
  
  int computeStatsMask(double stats[7],int pixel_type,
                       const int mask_i0,const int mask_j0,const int mask_ni,const int mask_nj,const int * const mask) {
   
    // pixels masked with mask[i0*mask_nj+j0] == 1 are included in stats...
 
    double vmin = HUGE_VAL;
    double vmax = -HUGE_VAL;
    double mean = 0.0;
    double rms = 0.0;
    double x_centroid[3] = { 0.0, 0.0, 0.0 };
    double pixel_area = getLengthScale();
    pixel_area *= pixel_area;
    double area = 0.0;
    int count = 0;
    double xP[3];
    for (int ipx =0; ipx < nx*ny; ipx++) {
      if ( pixel_flag[ipx] == pixel_type ) {
        const int iP = ipx%nx;
        const int jP = (int)ipx/nx;
        if ((iP >= mask_i0)&&(iP < mask_i0+mask_ni)&&(jP >= mask_j0)&&(jP < mask_j0+mask_nj)&&(mask[(iP-mask_i0)*mask_nj+(jP-mask_j0)] == 1)) {
          const double val = pixel_data[ipx];
          if (val<vmin) vmin=val;
          if (val>vmax) vmax=val;
          mean += val;
          rms += val*val;
          convertPxtoXp(iP,jP,xP);
          FOR_I3 x_centroid[i] += xP[i];
          ++count;
          area += pixel_area;
        }
      }
    }

    if (count > 0) {
      mean /= double(count);
      rms = sqrt(max(0.0,rms/double(count)-mean*mean));
      FOR_I3 x_centroid[i] /= double(count);
    }

    stats[0] = x_centroid[0];
    stats[1] = x_centroid[1];
    stats[2] = x_centroid[2];
    stats[3] = vmin;
    stats[4] = vmax;
    stats[5] = mean;
    stats[6] = rms;
    stats[7] = area;

    return count; // pixel count

  }
  
  double computeSurfaceArea(bool &b_markSurfaces) {
    const double pixel_area = metadata->getLengthScale()*metadata->getLengthScale();
    double surface_area = 0.0;
    //This doesn't do anything?
    //PngDataChunk * zoNe = getPngDataChunk("zoNe");
    for (int ipx =0; ipx < nx*ny; ipx++) {
      if ( pixel_flag[ipx] == SURFACE_PIXEL ){ //only including non-data surfaces
        surface_area += pixel_area;
        if (b_markSurfaces){  //mark region green for visual confirmation
          pixel_flag[ipx] = FLAG_GREEN;
        }
      }
    }
    return surface_area;
  }//computeSurfaceArea()

//places result of diff in diffImg
  void diff2(PngData * diffImg, double data_range[2]) {
    assert(diffImg->pixel_data!=NULL&&diffImg->pixel_flag!=NULL);
    assert(diffImg->nx==nx && diffImg->ny==ny); //images must be the same size for now
    assert(diffImg->metadata->getVarId()==metadata->getVarId()); //images must display the same variable
    //TODO check transformation matrix to confirm the same simulation space data is being compared?
    //     rounding error could be problematic between geometries?

    for (int ipx =0; ipx < nx*ny; ipx++) {
      if ( !(pixel_flag[ipx] < 0 && diffImg->pixel_flag[ipx] < 0) ){
        if ( pixel_flag[ipx] >= 0 && diffImg->pixel_flag[ipx] >= 0) {
          //both images have data here, compute the diff
          diffImg->pixel_data[ipx] -= pixel_data[ipx];
        }
        else if (diffImg->pixel_flag[ipx] < 0){
          //diffImg has no data here, color red
          diffImg->pixel_flag[ipx] = -2;
        }
        else{
          //this image has no data here, color diffImg green
          diffImg->pixel_flag[ipx] = -3;
        }
      }
    }
  }
  //places result of diff in diffImg
  bool diff(PngData * diffImg, const bool b_signed, string &msg) {
    assert(diffImg->pixel_data!=NULL&&diffImg->pixel_flag!=NULL);
    assert(diffImg->nx==nx && diffImg->ny==ny); //images must be the same size for now
    //TODO check transformation matrix to confirm the same simulation space data is being compared?
    //     rounding error could be problematic between geometries?

    vector<pixel_type> pixelTypeVec;
    getPixelTypes(pixelTypeVec);
    vector<pixel_type> diffPixelTypeVec;
    diffImg->getPixelTypes(diffPixelTypeVec);
    if (pixelTypeVec.size()!=diffPixelTypeVec.size()){
      msg = "diff() - Image contains different data types than reference";
      return false;
    }
    for (int ipt=0; ipt<pixelTypeVec.size(); ++ipt){
      bool b_foundType = false;
      bool b_foundVar = false;
      int jpt;
      for (jpt=0; jpt<diffPixelTypeVec.size(); ++jpt){
        if (pixelTypeVec[ipt].flag==diffPixelTypeVec[jpt].flag){
          b_foundType = true;
          if (pixelTypeVec[ipt].varId==diffPixelTypeVec[jpt].varId)
            b_foundVar = true;
          break;
        }
      }
      stringstream ss;
      if (!b_foundType){
        ss << "diff() - Image missing data of type " << pixelTypeVec[ipt].name;
        msg = ss.str();
        return false;
      }
      else if (!b_foundVar){
        ss << "diff() - Warning, image contains different var (" << diffPixelTypeVec[jpt].varId << ") than reference (" << pixelTypeVec[ipt].varId << ")";
        msg = ss.str();
      }
    }

    double dRangeVolume[2]    = {HUGE_VAL,-HUGE_VAL};
    double dRangeSurface[2]   = {HUGE_VAL,-HUGE_VAL};
    double dRangeParticles[2] = {HUGE_VAL,-HUGE_VAL};
    double dRangeIso[2]       = {HUGE_VAL,-HUGE_VAL};
    int unrecognized_data_type = 0;
    for (int ipx =0; ipx < nx*ny; ipx++) {
      if ( !(pixel_flag[ipx] < VOLUME_DATA_PIXEL && diffImg->pixel_flag[ipx] < VOLUME_DATA_PIXEL) ){
        if ( pixel_flag[ipx] == diffImg->pixel_flag[ipx] ) {
          //both images have data here of the same type, compute the diff
          diffImg->pixel_data[ipx] -= pixel_data[ipx];
          if (!b_signed) //take absolute value of diff
            diffImg->pixel_data[ipx] = fabs(diffImg->pixel_data[ipx]);

          switch (pixel_flag[ipx]){
            case VOLUME_DATA_PIXEL:
              if(diffImg->pixel_data[ipx]<dRangeVolume[0])
                dRangeVolume[0] =  diffImg->pixel_data[ipx];
              if(diffImg->pixel_data[ipx]>dRangeVolume[1])
                dRangeVolume[1] =  diffImg->pixel_data[ipx];
              break;
            case SURFACE_DATA_PIXEL:
              if(diffImg->pixel_data[ipx]<dRangeSurface[0])
                dRangeSurface[0] =  diffImg->pixel_data[ipx];
              if(diffImg->pixel_data[ipx]>dRangeSurface[1])
                dRangeSurface[1] =  diffImg->pixel_data[ipx];
              break;
            case PARTICLE_DATA_PIXEL:
              if(diffImg->pixel_data[ipx]<dRangeParticles[0])
                dRangeParticles[0] =  diffImg->pixel_data[ipx];
              if(diffImg->pixel_data[ipx]>dRangeParticles[1])
                dRangeParticles[1] =  diffImg->pixel_data[ipx];
              break;
            case ISO_DATA_PIXEL:
              if(diffImg->pixel_data[ipx]<dRangeIso[0])
                dRangeIso[0] =  diffImg->pixel_data[ipx];
              if(diffImg->pixel_data[ipx]>dRangeIso[1])
                dRangeIso[1] =  diffImg->pixel_data[ipx];
              break;
            default:
              ++unrecognized_data_type;
              //COUT1("Warning - unrecognized pixel data type " << pixel_data[ipx]);
          }
        }
        else if (diffImg->pixel_flag[ipx] <= VOLUME_DATA_PIXEL){
          //diffImg has no data here, color red
          diffImg->pixel_flag[ipx] = FLAG_RED;
        }
        else if (pixel_flag[ipx] <= VOLUME_DATA_PIXEL){
          //this image has no data here, color diffImg green
          diffImg->pixel_flag[ipx] = FLAG_GREEN;
        }
        else { //data type mismatch
          diffImg->pixel_flag[ipx] = FLAG_BLUE;
        }
      }
    }
    if (diffImg->b_data)
      diffImg->setRange(dRangeVolume);
    if (diffImg->b_dataSurface)
      diffImg->setRangeSurface(dRangeSurface);
    if (diffImg->b_dataParticles)
      diffImg->setRangeParticles(dRangeParticles);
    if (diffImg->b_dataIso)
      diffImg->setRangeIso(dRangeIso);

    stringstream ss;
    if (unrecognized_data_type>0){
      ss << "diff() - " << unrecognized_data_type << " pixel(s) with unrecognized data types";
      msg = ss.str();
      return false;
    }
    else{
      if (msg.length()==0){
        ss << "diff() - image diff successful";
        msg = ss.str();
      }
      //else warnings already in msg...
      return true;
    }

  }

  class color_blend_info {
    public:
      double color[3];
      double factor;

      color_blend_info() {
        for (int i=0; i<3; ++i) color[i] = 0.0;
        factor = 1.0;
      };

      color_blend_info(const double _color[3],const double _factor) {
        for (int i=0; i<3; ++i) color[i] = _color[i];
        factor = _factor;
      }

      color_blend_info(const double r, const double g, const double b,const double _factor) {
        color[0] = r;
        color[1] = g;
        color[2] = b;
        factor = _factor;
      }

      ~color_blend_info() {};
  };

  int diffWithStats(PngData * diffImg,const string metric,const string outputPrefix="image",const double norm_tol=-1.0,const bool force_image=false) {
    // returns error with bits set for different types of errors

    // metadata and image property checks
    bool size_diff = false;
    if (diffImg->nx!=nx) {
      cerr << "image-discrepancy nx " << nx << " " << diffImg->nx << endl;
      size_diff = true;
    }
    if (diffImg->ny!=ny) {
      cerr << "image-discrepancy ny " << ny << " " << diffImg->ny << endl;
      size_diff = true;
    }

    bool meta_diff = false;
    if (diffImg->metadata->getVarId() != metadata->getVarId()) {
      meta_diff = true;
      cerr << "metadata-discrepancy VAR "  << metadata->getVarId() << " " << diffImg->metadata->getVarId() << endl;
    }
    if (diffImg->metadata->getSurfVarId() != metadata->getSurfVarId()) {
      meta_diff = true;
      cerr << "metadata-discrepancy VAR_ON_SURFACE "  << metadata->getSurfVarId() << " " << diffImg->metadata->getSurfVarId() << endl;
    }
    if (diffImg->metadata->getPartVarId() != metadata->getPartVarId()) {
      meta_diff = true;
      cerr << "metadata-discrepancy VAR_ON_PARTICLE "  << metadata->getPartVarId() << " " << diffImg->metadata->getPartVarId() << endl;
    }
    if (diffImg->metadata->getIsoVarId() != metadata->getIsoVarId()) {
      meta_diff = true;
      cerr << "metadata-discrepancy VAR_ON_ISO "  << metadata->getIsoVarId() << " " << diffImg->metadata->getIsoVarId() << endl;
    }
    if (diffImg->metadata->getTime() != metadata->getTime()) {
      meta_diff = true;
      cerr << "metadata-discrepancy TIME "  << metadata->getTime() << " " << diffImg->metadata->getTime() << endl;
    }
    if (diffImg->metadata->getRangeMin() != metadata->getRangeMin()) {
      meta_diff = true;
      cerr << "metadata-discrepancy RANGE_MIN "  << metadata->getRangeMin() << " " << diffImg->metadata->getRangeMin() << endl;
    }
    if (diffImg->metadata->getRangeMax() != metadata->getRangeMax()) {
      meta_diff = true;
      cerr << "metadata-discrepancy RANGE_MAX "  << metadata->getRangeMax() << " " << diffImg->metadata->getRangeMax() << endl;
    }
    if (diffImg->metadata->getRangeMin() != metadata->getRangeMin()) {
      meta_diff = true;
      cerr << "metadata-discrepancy RANGE_MIN "  << metadata->getRangeMin() << " " << diffImg->metadata->getRangeMin() << endl;
    }
    if (diffImg->metadata->getLengthScale() != metadata->getLengthScale()) {
      meta_diff = true;
      cerr << "metadata-discrepancy LENGTH_SCALE "  << metadata->getLengthScale() << " " << diffImg->metadata->getLengthScale() << endl;
    }
    if (diffImg->metadata->getColorMapName() != metadata->getColorMapName()) {
      meta_diff = true;
      cerr << "metadata-discrepancy COLORMAP "  << metadata->getColorMapName() << " " << diffImg->metadata->getColorMapName() << endl;
    }

    if (size_diff) {
      cout << "cannot compute image diff when dimensions are different; skipping";
      return -1;
    }

    //TODO check transformation matrix to confirm the same simulation space data is being compared?
    //     rounding error could be problematic between geometries?

    const double one_o_3 = 1.0/3.0;
    const double one_o_255 = 1.0/255.0;
    const color_blend_info bg_cbi(255.0,255.0,255.0,0.75);
    const color_blend_info prob_cbi(255.0,0.0,0.0,0.25);

    // rgb buffer diff containers (per r,g,b channel)
    int mae[3] = {0,0,0};
    int mse[3] = {0,0,0};
    int n_pixel_diffs = 0;

    double data_range[2] = {HUGE_VAL,-HUGE_VAL};
    int diff[3];
    bool write_image = false;

    for (int ipx =0; ipx < nx*ny; ipx++) {
      bool rgb_same = true;
      for (int i=0; i<3; ++i) {
        diff[i] = abs(int(buffer[ipx][i]) - int(diffImg->buffer[ipx][i]));
      }
      if (diff[0]+diff[1]+diff[2] > 0) rgb_same = false;

      if ( rgb_same ) {
        // blend with white
        for (int i=0; i<3; ++i) {
          diffImg->buffer[ipx][i] = char(floor(sqrt((1.0 - bg_cbi.factor) * double(buffer[ipx][i])*double(buffer[ipx][i]) +  bg_cbi.factor*bg_cbi.color[i]*bg_cbi.color[i])));
        }
      }
      else {
        assert(!rgb_same);
        write_image = true;
        for (int index=0; index<3; ++index) {
          mae[index] += diff[index];
          mse[index] += diff[index]*diff[index];
        }
        ++n_pixel_diffs;

        //TODO could put tolerance on diff
        // magnitude of redness scales with diff
        double diff_mag = 0;
        for (int i=0; i<3; ++i) diff_mag += diff[i];
        diff_mag = diff_mag*one_o_3*one_o_255;

        const double blend_fac = prob_cbi.factor + (1.0-prob_cbi.factor)*diff_mag;
        // color red
        for (int i=0; i<3; ++i) {
          diffImg->buffer[ipx][i] = char(floor(sqrt((1.0 - blend_fac) * double(buffer[ipx][i])*double(buffer[ipx][i]) +  blend_fac*prob_cbi.color[i]*prob_cbi.color[i])));
        }
      }
    }
    // output statistics (rgb based)
    if ((metric != "MAE") && (metric != "MSE") && (metric != "RMSE") && (metric != "AE")) {
      cout << " > unknown metric specified; defaulting to Absolute Error (# of diffing pixels)" << endl;
    }
    else {
      cout << " > diff metric: " << metric << endl;
    }

    const double one_over_npixels = 1.0 / double(nx*ny);
    const double one_over_qr = 1.0/65535.0;
    if (metric == "MAE" || metric == "mae") {
      double d_mae[3];
      FOR_I3 d_mae[i] = double(mae[i])*257.0;
      const double abs_diff = (d_mae[0]+d_mae[1]+d_mae[2])*one_o_3*one_over_npixels;
      const double norm_diff = (d_mae[0]+d_mae[1]+d_mae[2])*one_o_3*one_over_npixels*one_over_qr;
      if (norm_diff < norm_tol) write_image = false;  // within threshold, then don't write out diff image
      cerr << "rgb " << abs_diff << " (" << norm_diff << ")" << endl;
      // cerr << "r " << d_mae[0]*one_over_npixels << " (" << d_mae[0]*one_over_npixels*one_over_qr << ")" << endl;
      // cerr << "g " << d_mae[1]*one_over_npixels << " (" << d_mae[1]*one_over_npixels*one_over_qr << ")" << endl;
      // cerr << "b " << d_mae[2]*one_over_npixels << " (" << d_mae[2]*one_over_npixels*one_over_qr << ")" << endl;
    }
    else if (metric == "MSE" || metric == "mse") {
      double d_mse[3];
      FOR_I3 d_mse[i] = double(mse[i])*257.0*one_o_255;
      const double abs_diff = (d_mse[0]+d_mse[1]+d_mse[2])*one_o_3*one_over_npixels;
      const double norm_diff = (d_mse[0]+d_mse[1]+d_mse[2])*one_o_3*one_over_npixels*one_over_qr;
      if (norm_diff < norm_tol) write_image = false;  // within threshold, then don't write out diff image
      cerr << "rgb " << abs_diff << " (" << norm_diff << ")" << endl;
      // cout << "r " << d_mse[0]*one_over_npixels << " (" << d_mse[0]*one_over_npixels*one_over_qr << ")" << endl;
      // cout << "g " << d_mse[1]*one_over_npixels << " (" << d_mse[1]*one_over_npixels*one_over_qr << ")" << endl;
      // cout << "b " << d_mse[2]*one_over_npixels << " (" << d_mse[2]*one_over_npixels*one_over_qr << ")" << endl;
    }
    else if (metric == "RMSE" || metric == "rmse") {
      double rmse[3];
      FOR_I3 rmse[i] = sqrt(double(mse[i])*66049.0*one_over_npixels);
      const double abs_diff = (rmse[0]+rmse[1]+rmse[2])*one_o_3;
      const double norm_diff = (rmse[0]+rmse[1]+rmse[2])*one_o_3*one_over_qr;
      if (norm_diff < norm_tol) write_image = false;  // within threshold, then don't write out diff image
      cerr << "rgb " << abs_diff << " (" << norm_diff << ")" << endl;
      // cout << "r " << rmse[0] << " (" << rmse[0]*one_over_qr << ")" << endl;
      // cout << "g " << rmse[1] << " (" << rmse[1]*one_over_qr << ")" << endl;
      // cout << "b " << rmse[2] << " (" << rmse[2]*one_over_qr << ")" << endl;
    }
    else {
      // absolute error (n_pixels)
      const double norm_diff = double(n_pixel_diffs)*one_over_npixels;
      if (norm_diff < norm_tol) write_image = false;  // within threshold, then don't write out diff image
      cerr << "rgb " << n_pixel_diffs << "(" << norm_diff << ")" << endl;
    }

    bool chunk_diffs[5] = {false,false,false,false,false};
    if (write_image || force_image) {
      //write rgb diff image with the global range
      chunk_diffs[0] = write_image;
      diffImg->metadata->setVarId("diff-rgb");
      string img_name = outputPrefix+".diff-rgb.png";
      cout << " > writing rgb diff to \"" << img_name << "\"" << endl;

      diffImg->metadata->setRangeMin(data_range[0]);
      diffImg->metadata->setRangeMax(data_range[1]);
      diffImg->write(img_name.c_str());
    }

    //------------- Zone diff
    chunk_diffs[1] = computeChunkDiff(diffImg,outputPrefix,"zoNe",2,metric,norm_tol,force_image);
    chunk_diffs[2] = computeChunkDiff(diffImg,outputPrefix,"dpTh",2,metric,norm_tol,force_image);
    chunk_diffs[3] = computeChunkDiff(diffImg,outputPrefix,"daTa",1,metric,norm_tol,force_image);
    chunk_diffs[4] = computeChunkDiff(diffImg,outputPrefix,"lgFr",1,metric,norm_tol,force_image);

    int err_bits = 0;
    if (size_diff) err_bits |= (1 << 0);
    if (meta_diff) err_bits |= (1 << 1);
    if (chunk_diffs[0]) err_bits |= (1 << 2);
    if (chunk_diffs[1]) err_bits |= (1 << 3);
    if (chunk_diffs[2]) err_bits |= (1 << 4);
    if (chunk_diffs[3]) err_bits |= (1 << 5);
    if (chunk_diffs[4]) err_bits |= (1 << 6);

    return err_bits;
  }

  bool computeChunkDiff(PngData * diffImg,const string outputPrefix,const string chunk,const int n_bytes,const string metric="MAE",const double norm_tol=-1.0,const bool force_image=false) {
    assert((n_bytes > 0) && (n_bytes<3));

    const color_blend_info bg_cbi(255.0,255.0,255.0,0.75);
    const color_blend_info prob_cbi(255.0,0.0,0.0,0.25);
    const double one_o_255 = 1.0/255.0;

    PngDataChunk * chunk0 = getPngDataChunk(chunk);
    PngDataChunk * chunk1 = diffImg->getPngDataChunk(chunk);
    if (chunk0 && chunk1) {
      int mae = 0;
      int mse = 0;
      int count_diff = 0;
      bool write_image = false;
      for (int ipx =0; ipx < nx*ny; ipx++) {
        bool var_same = true;
        double valDiff = 0.0;

        unsigned char val_uc[2];
        if (n_bytes == 2) {
          uint2 val_us0,val_us1;
          chunk0->get(ipx,val_uc);
          val_us0 = (val_uc[1] << 8 | val_uc[0]);
          chunk1->get(ipx,val_uc);
          val_us1 = (val_uc[1] << 8 | val_uc[0]);
          valDiff = double(fabs(int(val_us0) - int(val_us1)));
        }
        else if (n_bytes == 1) {
          chunk0->get(ipx,&val_uc[0]);
          chunk1->get(ipx,&val_uc[1]);
          valDiff = double(fabs(int(val_uc[0]) - int(val_uc[1])));
        }
        else {
          cout << "Warning: cannot compute "<<n_bytes<<"-byte chunk diff; skipping" << endl;
          return false;
        }

        if (valDiff > 0.0) {
          ++count_diff;
          var_same = false;
          mae += int(valDiff);
          mse += int(valDiff)*int(valDiff);
        }

        if ( var_same ) {
          // blend with white
          for (int i=0; i<3; ++i) {
            diffImg->buffer[ipx][i] = char(floor(sqrt((1.0 - bg_cbi.factor) * double(buffer[ipx][i])*double(buffer[ipx][i]) +  bg_cbi.factor*bg_cbi.color[i]*bg_cbi.color[i])));
          }
        }
        else {
          write_image = true;
          for (int i=0; i<n_bytes; ++i) valDiff *= one_o_255;  // normalize by bit-range of zone chunk

          const double blend_fac = prob_cbi.factor + (1.0-prob_cbi.factor)*valDiff;
          // color red
          for (int i=0; i<3; ++i) {
            diffImg->buffer[ipx][i] = char(floor(sqrt((1.0 - blend_fac) * double(buffer[ipx][i])*double(buffer[ipx][i]) +  blend_fac*prob_cbi.color[i]*prob_cbi.color[i])));
          }
        }
      }

      // metric reporting (to stderr to separate from rest of output)
      const double one_o_npixels = 1.0 / double(nx*ny);
      double d_mae = double(mae) * one_o_npixels;
      double norm_mae = d_mae;
      double d_mse = double(mse) * one_o_npixels;
      double norm_mse = d_mse;
      for (int i=0; i<n_bytes; ++i) {
        norm_mae *= one_o_255;
        norm_mse *= one_o_255;
      }
      if (metric == "MAE" || metric == "mae") {
        cerr << chunk << " " << d_mae << " (" << norm_mae << ")" << endl;
        if (norm_mae < norm_tol) write_image = false;  // within threshold, then don't write out diff image
      }
      else if (metric == "MSE" || metric == "mse") {
        cerr << chunk << " " << d_mse << " (" << norm_mse << ")" << endl;
        if (norm_mse < norm_tol) write_image = false;  // within threshold, then don't write out diff image
      }
      else if (metric == "RMSE" || metric == "rmse") {
        cerr << chunk << " " << sqrt(d_mse) << " (" << sqrt(d_mse)/((n_bytes*256.0)-1) << ")" << endl;
        if ((sqrt(d_mse)/((n_bytes*256.0)-1)) < norm_tol) write_image = false;  // within threshold, then don't write out diff image
      }
      else {
        // absolute error (n_pixels)
        const double norm_diff = double(count_diff)*one_o_npixels;
        if (norm_diff < norm_tol) write_image = false;  // within threshold, then don't write out diff image
        cerr << chunk << " " << count_diff << "(" << norm_diff << ")" << endl;
      }

      if (write_image || force_image) {
        string diff_name = "diff-"+chunk;
        string img_name = outputPrefix+"."+diff_name+".png";
        diffImg->metadata->setVarId(diff_name);
        cout << " > writing zone diff to \"" << img_name << "\"" << endl;
        double di_range[2] = {0.0,255.0};
        diffImg->metadata->setRangeMin(di_range[0]);
        diffImg->metadata->setRangeMax(di_range[1]);
        diffImg->write(img_name.c_str());
      }

      return write_image;
    }
    else {
      cout << " > chunk \"" << chunk << "\" did not exist in these images; skipping" << endl; 
    }
    return false;
  }

  int probePoint(const double xin[3], double xout[3],int ij[2], double &val) {

    // find the image pixel ij corresponding to xin by ignored its position along
    // the np (depth) image axis

    // transformation info
    double x0[3], np[3], e0[3], e1[3];
    x0[0] = metadata->transformMat[12];
    x0[1] = metadata->transformMat[13];
    x0[2] = metadata->transformMat[14];
    np[0] = metadata->transformMat[8];
    np[1] = metadata->transformMat[9];
    np[2] = metadata->transformMat[10];
    e0[0] = metadata->transformMat[0];
    e0[1] = metadata->transformMat[1];
    e0[2] = metadata->transformMat[2];
    e1[0] = metadata->transformMat[4];
    e1[1] = metadata->transformMat[5];
    e1[2] = metadata->transformMat[6];

    double np_mag = 0.0;
    for (int i =0; i < 3; ++i)
      np_mag += np[i]*np[i];
    np_mag = sqrt(np_mag);

    for (int i =0; i < 3; ++i)
      np[i] /= np_mag;

    double e0_mag = 0.0;
    for (int i =0; i < 3; ++i)
      e0_mag += e0[i]*e0[i];
    e0_mag = sqrt(e0_mag);

    for (int i =0; i < 3; ++i)
      e0[i] /= e0_mag;

    double e1_mag = 0.0;
    for (int i = 0; i < 3; ++i)
      e1_mag += e1[i]*e1[i];
    e1_mag = sqrt(e1_mag);

    for (int i =0; i < 3; ++i)
      e1[i] /= e1_mag;

    // project xin onto image plane (depth=0)
    double xin_diff[3];
    for (int i =0; i < 3; ++i)
      xin_diff[i] = xin[i] - x0[i];

    double dist = 0.0;
    for (int i = 0; i < 3; ++i)
      dist += xin_diff[i]*np[i];

    // get image i,j of projected point
    double ii = (xin_diff[0]*e0[0] + xin_diff[1]*e0[1] + xin_diff[2]*e0[2])/e0_mag;
    double jj = (xin_diff[0]*e1[0] + xin_diff[1]*e1[1] + xin_diff[2]*e1[2])/e1_mag;

    // compute xout at ij using image depth value
    ij[0] = (int)round(ii);
    ij[1] = (int)round(jj);
    convertPxtoXp(ij[0],ij[1],xout);

    int nn = ij[1]*nx + ij[0];
    if (pixel_flag[nn] >= VOLUME_DATA_PIXEL) {
      val = pixel_data[nn];
      return pixel_flag[nn];
    }
    else {
      // is there something nearby?...
      int i_closest = -1;
      int j_closest = -1;
      double r2_closest = 0.0;
      for (int i = max(0,ij[0]-10); i < min(nx,ij[0]+10); ++i) {
	for (int j = max(0,ij[1]-10); j < min(ny,ij[1]+10); ++j) {
	  if (pixel_flag[j*nx+i] >= VOLUME_DATA_PIXEL) {
	    double r2 = (ii-i)*(ii-i) + (jj-j)*(jj-j);
	    if ((i_closest == -1)||(r2 < r2_closest)) {
	      i_closest = i;
	      j_closest = j;
	      r2_closest = r2;
	    }
	  }
	}
      }
      if (i_closest == -1) {
	cout << "Warning: point " << COUT_VEC(xin) << " closest pixel " << ij[0] << " " << ij[1] <<
	  " has no data and there is no data in nearby pixels. Setting value to 0 (may affect comparison)" << endl;
	val = 0.0;
        return pixel_flag[nn];
      }
      else {
	cout << "Warning: point " << COUT_VEC(xin) << " closest pixel " << ij[0] << " " << ij[1] <<
	  " has no data. Using closest pixel with data at " << i_closest << " " << j_closest << endl;
	val = pixel_data[j_closest*nx+i_closest];
        return pixel_flag[j_closest*nx+i_closest];
      }
    }

  }

  int probeIJ(const int ii,const int jj, double &val) {

    const int nn = jj*nx + ii;
    if (pixel_flag[nn] >= VOLUME_DATA_PIXEL){
      val = pixel_data[nn];
    }
    else{
      val = 0;
      cout << "Warning: no data found at pixel " << ii << ", " << jj << endl;
    }

    return pixel_flag[nn];

  }

  int convertPxtoXp(const int ip,const int jp,double xp[3]){
    assert(dpTh);
    const int nn = jp*nx + ip;
    uint2 dpTh_us = getDepth(nn);
    if (dpTh_us==65535){
      cout << "Warning: requesting simulation coordinates for background pixel " << ip << " " << jp << endl;
    }
    for (int i =0; i < 3; ++i) {
      xp[i] = ip*metadata->transformMat[0+i] +
              jp*metadata->transformMat[4+i] +
              dpTh_us*metadata->transformMat[8+i] +
              metadata->transformMat[12+i];
    }
    return pixel_flag[nn];
  }
  
  //build nodal buffer for a single user specified data type (pixel_flag_val)
  void buildNodalData(double * no_data, const int pixel_flag_val){
    for (int ip = 0; ip<=nx; ++ip){
      for (int jp = 0; jp<=ny; ++jp){
        int ip4[4] = {ip,ip,ip-1,ip-1};
        int jp4[4] = {jp,jp-1,jp-1,jp};
        no_data[jp*(nx+1)+ip] = 0.0;
        double weight = 0.0;
        for (int k = 0; k<4; ++k){
          if (ip4[k]>=0&&jp4[k]>=0&&ip4[k]<nx&&jp4[k]<ny&&pixel_flag[jp4[k]*nx+ip4[k]]==pixel_flag_val){
            no_data[jp*(nx+1)+ip] += pixel_data[jp4[k]*nx+ip4[k]];
            weight += 1.0;
          }
        }
        if (weight>0.0){
          no_data[jp*(nx+1)+ip]/=weight;
        }
      }
    }

  }

  bool interpFromNodalData(double &val, const double &ipf, const double &jpf, const double * no_data,const int pixel_flag_val){
     int ip = int(floor(ipf+0.5));
     int jp = int(floor(jpf+0.5));
     if (ip<0||jp<0||ip>=nx||jp>=nx||pixel_flag[jp*nx+ip]!=pixel_flag_val){
       return false;
     }
     else{
       double ip0 = ip - 0.5f;
       double jp0 = jp - 0.5f;
       double ip1 = ip + 0.5f;
       double jp1 = jp + 0.5f;
       double area[4] = {(ip1-ipf)*(jp1-jpf),(ip1-ipf)*(jpf-jp0),(ipf-ip0)*(jpf-jp0),(ipf-ip0)*(jp1-jpf)};
       val = area[0]*no_data[jp*(nx+1)+ip] +
             area[1]*no_data[(jp+1)*(nx+1)+ip] +
             area[2]*no_data[(jp+1)*(nx+1)+ip+1] +
             area[3]*no_data[jp*(nx+1)+ip+1];

       return true;
     }

  }

  double getLengthScale(){
    return metadata->getLengthScale();
  }

  void dumpInfo(ostream& ofile){
    ofile << "   Reporting image metadata by chunk type...\n" << endl;
    metadata->dumpInfo(ofile);
    for (list<PngDataChunk*>::iterator it = pdcList.begin(); it != pdcList.end(); it++){
      if ((*it)!=NULL) {
        //compute some chunk stats
        double minval = HUGE_VAL;
        double maxval = -HUGE_VAL;
        const int bytesPerPixel = (*it)->getBytesPerPixel();
        unsigned char * chnk_uc = new unsigned char[bytesPerPixel];
        uint2 chnk_us = 0;
        for (int ich =0; ich< ny*nx; ++ich){
          (*it)->get(ich,chnk_uc);
          if (bytesPerPixel==1){
            chnk_us = chnk_uc[0];
          }
          else if (bytesPerPixel==2){
            chnk_us = (chnk_uc[1] << 8 | chnk_uc[0]); //little endian
          }
          if (minval > chnk_us)
            minval = (double) chnk_us;
          if (maxval < chnk_us)
            maxval = (double) chnk_us;
        }
        delete[] chnk_uc;
        ofile << "*****" << (*it)->getName() << " metadata*****" << endl;
        ofile << "compressed size in bytes: " << (*it)->getZdata_size() << endl;
        ofile << "bytes per pixel: " << (*it)->getBytesPerPixel() << endl;
        ofile << "range (min max): " << minval << " " << maxval << endl;
      }
    }
  }

  void blankWithMask(const int * pixel_flag_mask){
    for (int ipx=0; ipx<nx*ny; ++ipx){
      if (pixel_flag[ipx] != pixel_flag_mask[ipx]){
        assert(pixel_flag_mask[ipx] == BACKGROUND_PIXEL);
        buffer[ipx][0]=73;
        buffer[ipx][1]=175;
        buffer[ipx][2]=205;

        pixel_data[ipx] = 0;
        pixel_flag[ipx] = BACKGROUND_PIXEL;
  
        unsigned char zoNe_uc[2] = {(unsigned char) 255, (unsigned char) 255};
        zoNe->set(ipx,zoNe_uc);
  
        unsigned char dpTh_uc[2] = {(unsigned char) 255, (unsigned char) 255};
        dpTh->set(ipx,dpTh_uc);
  
        if (daTa){
          unsigned char daTa_uc = (unsigned char) 0;
          daTa->set(ipx,&daTa_uc);
        }
        if (lgFr){
          unsigned char lgFr_uc = (unsigned char) 0;
          lgFr->set(ipx,&lgFr_uc);
        }
      }
    }
  }

};//class PngData


//dotproducts operate across all data types (combining multiple
//data types into the result if present)
//..this is the behavior pod and dmd routines expect
//TODO consider smarter dot product by type?

inline double dotproduct(const int i, const int j, PngData * g1, PngData * g2) {

  assert ( g1->getNx() == g2->getNx() ) ;
  assert ( g1->getNy() == g2->getNy() ) ;

  int k = i*g1->getNx() + j ;
  if (g1->pixel_flag[k] == g2->pixel_flag[k]) {
    if (g1->pixel_flag[k] >= PngData::VOLUME_DATA_PIXEL) {
      double tmp1 = g1->pixel_data[k];
      double tmp2 = g2->pixel_data[k];
      return tmp1*tmp2;
    }
  }
  return 0.0;
}//dotproduct(i, j, Png, Png)

inline double dotproduct(PngData * g1, PngData * g2 ) {

  // only consider a valid pixel if both images, think the
  // pixel is valid; otherwise, this is an issue...
  assert ( g1->getNx() == g2->getNx() ) ;
  assert ( g1->getNy() == g2->getNy() ) ;

  double dotp = 0.0;
  int n       = 0;
  int nerr    = 0;
  for (int i =0; i < g1->getNx()*g1->getNy() ; ++i) {

    // only contribute if both pixels are valid ...
    if (g1->pixel_flag[i] == g2->pixel_flag[i]) {
      // both pixels should be bahama blue...
      if ( g1->pixel_flag[i] >= PngData::VOLUME_DATA_PIXEL ) {
        n++;
        double tmp1 = g1->pixel_data[i];
        double tmp2 = g2->pixel_data[i];
        dotp       += tmp1*tmp2;
      }
    }
    else {
      //cout << "Warning: pixels are not consistent! " << endl;
      nerr++;
    }
  }//i

  if ( nerr > 0 )
    cout << "Warning: pct of inconsitent px = " << double(nerr)/double(g1->getNx()*g2->getNy()) << endl;

  dotp /= double(n);
  return dotp;
}//dotproduct(Png, Png)

inline double dotproduct(double * g1, PngData * g2 ) {
  double dotp = 0.0;
  int n       = 0;
  for (int i =0; i < g2->getNx()*g2->getNy() ; ++i) {

    if ( g2->pixel_flag[i] >= PngData::VOLUME_DATA_PIXEL ) {
      n++;
      double tmp2 = g2->pixel_data[i];
      dotp       += g1[i]*tmp2;
    }
  }//i

  dotp /= double(n);
  return dotp;
}//dotproduct(double, Png)


// for 2-stacked composite dmd...
class Composite {
public:
  PngData * phi0;
  PngData * phi1;
  bool bMem;

  Composite() : phi0(NULL), phi1(NULL), bMem(false) {}
  Composite(const string& _infile0, const string& _infile1) {

    phi0 = new PngData() ;
    phi0->read(_infile0.c_str());
    phi0->finalizeRead();

    phi1 = new PngData() ;
    phi1->read(_infile1.c_str());
    phi1->finalizeRead();

    bMem = true;
  }

  ~Composite() {
    if ( bMem) {
      if ( phi0 != NULL) delete phi0;
      if ( phi1 != NULL) delete phi1;
    }
  }
};

// need a custom dot product for a composite pair..
inline double dotproduct(Composite* _c0, Composite* _c1) {
  return dotproduct(_c0->phi0, _c1->phi0) + dotproduct(_c0->phi1,_c1->phi1);
}
#endif
