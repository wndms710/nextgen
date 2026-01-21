#ifndef _SIMPLE_RESTART_HPP_
#define _SIMPLE_RESTART_HPP_

#include "SimpleSurface.hpp"

class SimpleRestart {

public:

  string filetype;
  string filename_mesh;
  string filename_data;

  SimpleSurface ss;

  int8 meshMaxOffset;
  map<const string,int8> cvD1OffsetMap;
  map<const string,int8> cvD2OffsetMap;

  SimpleRestart() {
    filetype = "unknown";
    filename_mesh = "";
    filename_data = "";
    ss.clear();
    meshMaxOffset = -1;
    cvD1OffsetMap.clear();
    cvD2OffsetMap.clear();
  }

  ~SimpleRestart() {
  }
  void clear();

  int init(Param * param);
  //helper functions to init
  int initMLES();
  int initMLES_surface();
  int initMLES_data();
  int initSLES();

  bool isInit() const;

  //Used by WRITE_PARAMS to get input file from result.sles
  bool readResultParams(char * &cbuf, const string& filename); //returns true if params found, must delete cbuf


};

#endif
