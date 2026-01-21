#ifndef _LES_IMAGE_MAPPER_
#define _LES_IMAGE_MAPPER_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <assert.h>
#include <unistd.h>
#include <sys/stat.h> 
#include <string.h>
#include <ctime>
#include "Defs.hpp"
#include "ByteSwap.hpp"
#include <dirent.h>
using namespace std;

#include "CtiCanvas.hpp"
#include "CvImageMap.hpp"
#include "CvPlaneMap.hpp"
#include "RestartHashUtilities.hpp"
#include "IntFlag.hpp"
#include "SimplePoint.hpp"
#include "SimpleGeom.hpp"

class LesImageMapper{

private:

  FILE * rfp; // restart file pointer
  bool rfp_eof;
  bool rfp_byte_swap;
  int8 rfp_offset;
  int io_version;
  bool got_ncv_global;
  int8 ncv_global;
  bool got_nbf_global;
  int8 nbf_global;
  int8 cvobf_offset;
  string encodedHash;
  bool bSetEncodedHash;
  int step;
  double time;

  // stuff for mles/sles paradigm
  FILE * sfp; // solution file pointer
  bool sfp_eof;
  bool sfp_byte_swap;
  int8 sfp_offset;
  string encodedHash_2;
  bool bSetEncodedHash_2;

  float * d2_buf;
 
public:

  float * mesh_buf;

private:

  map<const string,int8> cvD1OffsetMap;
  map<const string,int8> cvD2OffsetMap;
  map<const string,int> i0ValueMap;
  map<const string,double> d0ValueMap;
  pair<int,int8> vvbboxPair; //mpi_size (nbb) and offset
  vector<pair<int8,int8> > nbfZnOffsetVec; //nbf_zn_global,local offset for each named boundary zone
 
  // stuff for mles/sles paradigm
  map<const string,int8> cvD1OffsetMap_2;
  map<const string,int8> cvD2OffsetMap_2;
  map<const string,int> i0ValueMap_2;
  map<const string,double> d0ValueMap_2;

  // stuff to store bbox record...
  int * ncv_bbox_buf;
  double (*vv_bbox_buf)[6];

public:

  LesImageMapper() : rfp(NULL), rfp_eof(true), rfp_byte_swap(false), rfp_offset(0), got_ncv_global(false), ncv_global(0), got_nbf_global(false), nbf_global(0), cvobf_offset(-1), encodedHash(""), bSetEncodedHash(false), step(0), time(0.0), d2_buf(NULL), mesh_buf(NULL) {
    sfp = NULL;
    sfp_eof = true;
    sfp_byte_swap = false;
    sfp_offset = 0;
    encodedHash_2 = "";
    bSetEncodedHash_2 = false;
    // bbox record stuff...
    ncv_bbox_buf = NULL;
    vv_bbox_buf = NULL;
  }

  LesImageMapper(const string &_filename,const string &_filename_2 = "") : rfp(NULL), rfp_eof(true), rfp_byte_swap(false), rfp_offset(0), got_ncv_global(false), ncv_global(0), got_nbf_global(false), nbf_global(0), cvobf_offset(-1), encodedHash(""), bSetEncodedHash(false), step(0), time(0.0), d2_buf(NULL), mesh_buf(NULL) {
    sfp = NULL;
    sfp_eof = true;
    sfp_byte_swap = false;
    sfp_offset = 0;
    encodedHash_2 = "";
    bSetEncodedHash_2 = false;
    // bbox record stuff...
    ncv_bbox_buf = NULL;
    vv_bbox_buf = NULL;
    init(_filename,_filename_2);
  }

  ~LesImageMapper(){
    // close the restart file pointer and reset the offset maps.
    fclose(rfp); rfp = NULL; rfp_eof = true;
    cvD1OffsetMap.clear();
    cvD2OffsetMap.clear();
    i0ValueMap.clear();
    d0ValueMap.clear();
    if (d2_buf!=NULL){
      delete[] d2_buf;
      d2_buf = NULL;
    }
    if (mesh_buf!=NULL){
      delete[] mesh_buf;
      mesh_buf = NULL;
    }
    DELETE(ncv_bbox_buf);
    DELETE(vv_bbox_buf);
    if (sfp) {
      fclose(sfp);
      sfp = NULL;
      sfp_eof = true;
      cvD1OffsetMap_2.clear();
      cvD2OffsetMap_2.clear();
      i0ValueMap_2.clear();
      d0ValueMap_2.clear();
    }
  }

  double getTime();

  void sparseReadCvR1(double * var,const CvImageMap& cim,const string &img_var);
  void sparseReadCvR1(double * var,const CvImageMap& cim,const string &img_var,const int8 &data_offset);

  void sparseReadCvR2(double (*varR2)[3],const CvImageMap& cim,const string &img_var);
  void sparseReadCvR2(double (*varR2)[3],const CvImageMap& cim,const string &img_var,const int8 &data_offset);
  void sparseReadCvR2Mag(double *varR2Mag,const CvImageMap& cim,const string &img_var,const int8 &data_offset);
  
  void getPointsInsideGeom(vector<SimplePoint>& pointVec,SimpleGeom * geom);
  
  void buildCvPlaneMap(CvPlaneMap& cvPlaneMap, const string &lengthscaleName);
  int buildCvImageMap(CvImageMap& cim, const CtiCanvas * maskCanvas, CvPlaneMap& cvPlaneMap, const string &lengthscaleName="r_vv");
  string buildCvImageMapFilename(const CtiCanvas * canvas, const float xp[3], const float np[3]);

  int buildCvImageMapSurface(CvImageMap& cim, const CtiCanvas * _canvas, const IntFlag &bfzone_flag, const string &lengthscaleName);
  string buildCvImageMapFilenameSurface(const CtiCanvas * canvas, const float bbox[6], const int &nst, const std::vector< PlaneData<float> > &vBlankPlaneData, const IntFlag &show_sz_flag);

  void computeCvPixelD2(const CvImageMap& cim, const CtiCanvas * canvas, const float xp[3], const float np[3]);

  int8 getRfpCvD1Offset(const string& name);
  int8 getRfpCvD2Offset(const string& name);

  // functions added for mles/sles paradigm
  int8 getSfpCvD1Offset(const string& name);
  int8 getSfpCvD2Offset(const string& name);
  void sparseReadCvR1_2(double * var,const CvImageMap& cim,const string &img_var,const int8 &data_offset);
  void sparseReadCvR2_2(double (*varR2)[3],const CvImageMap& cim,const string &img_var,const int8 &data_offset);
  void sparseReadCvR2Mag_2(double *varR2Mag,const CvImageMap& cim,const string &img_var,const int8 &data_offset);

private:
  void init(const string &_filename,const string &_filename_2="");

  void my_fseek(FILE * fp, int8 offset);

  int8 atoi8(const string& var);

  void checkContiguous(FILE * fp,const int n);

  int8 advanceRfp(const int ugp_io_tag,const string& name);
  bool setRfpHashId();
  bool getRfpI0Value(const string& name,int &i0Value);
  bool getRfpD0Value(const string& name,double &d0Value);
  void getNcvGlobalIoV3();
  int8 getRfpVVBBoxOffset(const string& name, int &nbb);
  bool checkVVBBoxIntersection(int &ncv_bbox,const int &ibb,const float xp[3], const float np[3]);
  bool checkVVBBoxBBoxIntersection(int &ncv_bbox,const int &ibb,const double bbox[6]);
  int8 getRfpCvobfOffsetAndNbfGlobal();
  bool getRfpNbfZnOffsetVec();

  // functions added for mles/sles paradigm
  bool setSfpHashId();
  int8 advanceSfp(const int ugp_io_tag,const string& name);
  bool getSfpI0Value(const string& name,int &i0Value);
  bool getSfpD0Value(const string& name,double &d0Value);

};

#endif
