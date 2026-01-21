#ifndef _STITCH_NS_HPP_
#define _STITCH_NS_HPP_

#include "CTI.hpp"
using namespace CTI;

#include "PartData.hpp"
using namespace PartData;

#include "HcpPointBuilder.hpp"

namespace StitchNS {

  extern HcpPointBuilder hcp;

  void finalize();
  void rebuildVd(const int nsmooth);
  void writeMles(const string& filename,const int io_version,const bool write_surface);
  void clearSmoothing();

  void processShow(Param * param,const bool b_help);
  void processWriteJson(Param * param,const bool b_help);
  void helpWriteJson();
  void processWriteImage(Param * param,const bool b_help);
  void helpWriteImage();
  void processForcePlanar(Param * param,const bool b_help);
  void processWriteMles(Param * param,const bool b_help);
  void helpWriteMles();
  void processClearSmoothing(Param *param,const bool b_help);
  void processSmoothLimit(Param *param,const bool b_help);
  void helpSmoothLimit();
  void processSmoothMode(Param *param,const bool b_help);
  void helpSmoothMode();
  void processCreaseAngleDegrees(Param * param,const bool b_help);
  void helpCreaseAngleDegrees();
  void processRcLocalThreshold(Param *param,const bool b_help);
  void helpRcLocalThreshold();
  void processRcLimit(Param *param,const bool b_help);
  void helpRcLimit();
  
} // namespace StitchNS

#endif
