#ifndef PNG_DATA_CHUNK_HPP
#define PNG_DATA_CHUNK_HPP

// ======================================================================
// PngImage.hpp
//
// This class can write a png image using libpng routines (i.e. png.h).
//
// Sample usage:
//
//
//
// History:
//  author: David Philips - March 2017
//  author: Frank Ham & Phuc Quang - July 2013
// ======================================================================

// ======================================================================
// Copyright (c) 2013 Cascade Technologies Inc.
// All Rights Reserved. www.cascadetechnologies.com
//
// This source is proprietary software and cannot be used, redistributed,
// or modified except under the terms of Cascade's Software License
// Agreement.
//
// THIS SOURCE CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY
// OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
// ======================================================================

#ifdef WITH_PNG

#include <zlib.h>

using namespace std;

class PngDataChunk {

protected:
  string name;
  int nx, ny;
  int bytesPerPixel;
  unsigned char * data;
  unsigned char * zdata;
  uLongf zdata_size;

public:

  static const int MAX_BYTES_PER_PIXEL = 4;

  PngDataChunk() : name(""), nx(0), ny(0), bytesPerPixel(0), data(NULL), zdata(NULL), zdata_size(0) {
  }

  PngDataChunk(const string _name, const int width, const int height, const int bpp) : name(_name), nx(width), ny(height), bytesPerPixel(bpp), data(NULL), zdata(NULL), zdata_size(0){
    assert(_name.length()==4);  //name must be exactly 4 chars like "dpTh"
    data = new unsigned char[width*height*bpp];
  }

  //contstructor used when reading chunks from png files
  PngDataChunk(const string _name, const int width, const int height, unsigned char * _zdata, uLongf _zdata_size) : name(_name), nx(width), ny(height), bytesPerPixel(0), data(NULL), zdata(NULL), zdata_size(_zdata_size) {
    zdata = new unsigned char[width*height*MAX_BYTES_PER_PIXEL]; //allocate max buffer size to inflate chunk data with unknown bytesPerPixel
    uLongf zlib_size = width*height*MAX_BYTES_PER_PIXEL;
    if (uncompress(zdata, &zlib_size, _zdata, _zdata_size) == Z_OK) { //inflate data successful
      assert(zlib_size%(width*height)==0);
      bytesPerPixel = zlib_size/(width*height); //compute bytesPerPixel
      data = new unsigned char[zlib_size];      //allocated exact size for uncompressed data
      memcpy(data,zdata,zlib_size);             //copy from (temp buffer) zdata to data
      //cout << " > Read and inflated chunk data for " << _name << endl;
    }
    else {
      cout << " > Warning, unable to inflate chunk data for " << _name << endl;
    }
    //remove local memory zdata which was used to inflate the data pointed to by _zdata
    //TODO could retain the compressed buffer for faster copying of images??
    delete[] zdata;
    zdata=NULL;
    //We now have uncompressed chunk data with it's name, nx, ny, and bytesPerPixel set.
    //zdata_size is set with the compressed length incase the user wishes to report this length without recomputing
  }

  ~PngDataChunk() {
    if ((data) != NULL) {
      delete[] data;
      data = NULL;
    }
    if ((zdata) != NULL) {
      delete[] zdata;
      zdata = NULL;
    }

  }

  int getNx() {return nx;}
  int getNy() {return ny;}

  void resize(const int nx_new,const int ny_new) {
    if (nx_new*ny_new <= nx*ny) {
      nx = nx_new;
      ny = ny_new;
    }
    else {
      assert(0);
    }
  }

  string getName() {return name;}

  char getName(const int &index){
    assert(index>=0&&index<4);
    return name.at(index);
  }
  int getSize(){
    if (zdata!=NULL) return zdata_size;
    else return nx*ny*bytesPerPixel;
  }
  //being able to specifically request zdata_size is useful
  //for some ping batch operations reporting chunk info.
  int getZdata_size(){
    return zdata_size;
  }
  int getBytesPerPixel(){
    return bytesPerPixel;
  }

  //get value of chunk at a single pixel
  void get(const int index, unsigned char *val){
    assert(index>=0&&index<nx*ny);
    memcpy(val,data+index*bytesPerPixel,bytesPerPixel);
  }
  

  unsigned char* getData(){
    if (zdata!=NULL) return zdata;
    else return data;
  }

  //pass in pointer to unsigned char * of length bytesPerPixel
  //set value of chunk at a single pixel
  void set(const int index, unsigned char *val){
    assert((index>=0)&&(index<nx*ny));
    memcpy(data+index*bytesPerPixel,val,bytesPerPixel);
  }

  //set full data buffer from an input buffer
  //used when reading unknown chunks
  void setData(unsigned char *readBuf){
    memcpy(data,readBuf,bytesPerPixel*ny*nx);
  }

  void compressChunk(){
      //buffer size is one byte less than uncompressed data size
      //if compression fails, reader can treat full length byte
      //stream as uncompressed data
      zdata_size = bytesPerPixel*nx*ny - 1;

      if (zdata==NULL) zdata = new unsigned char[zdata_size];

      if( compress(zdata, &zdata_size, data, bytesPerPixel*nx*ny) != Z_OK ){
        //unable to compress data
        zdata_size = 0;
        delete[] zdata;
        zdata = NULL;
      }
      //compression successful, zdata_size now contains length of compressed
      //data in zdata
  }


};

#else

class PngImage {

public:
  PngImage() {
    cout << "Warning: you must compile -D WITH_PNG and define LIBPNG to use PngImage." << endl;
  }

  PngImage(const int width, const int height, const int _prec_level) {
    cout << "Warning: you must compile -D WITH_PNG and define LIBPNG to use PngImage." << endl;
  }

  PngImage(const int width, const int height, const int _prec_level, unsigned char *inbuffer) {
    cout << "Warning: you must compile -D WITH_PNG and define LIBPNG to use PngImage." << endl;
  }

  void clear() { }

  void set(const int i, const int j, const unsigned char r, const unsigned char g, const unsigned char b) { }

  void write(const char * filename) { }

  void read(const char * filename, const bool withPixelData) { }

  int getNx() {return 0;}
  int getNy() {return 0;}
};

#endif

#endif
