#ifndef PNG_IMAGE_HPP
#define PNG_IMAGE_HPP

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

#include <png.h>
#include <zlib.h>
#include <cstring>
#include <map>
#include <cstdlib>
#include <list>
#include <iomanip>
#include <stdint.h>
#include "Defs.hpp"
#include "PngDataChunk.hpp"
#include "ImageMetadata.hpp"

using namespace std;

class PngImage {

protected:
  int nx, ny;
  ImageMetadata * metadata;
  unsigned char (*buffer)[3]; // 0 <= r, g, b < 255
  list<PngDataChunk*> pdcList;

  bool b_cleanupMem;  //only true when reading png files from disk???

public:

  PngImage() : nx(0), ny(0), metadata(NULL), buffer(NULL), b_cleanupMem(false) {
  }

  PngImage(const int width, const int height) : nx(width), ny(height), metadata(NULL), buffer(NULL), b_cleanupMem(false) {
  }

  PngImage(const int width, const int height, unsigned char (*_rgb)[3]) : nx(width), ny(height), metadata(NULL), buffer(_rgb), b_cleanupMem(false) {
  }

  ~PngImage() {
    //No memory cleanup in PngImage...handle in calling class
    //Exception when using the read() function to build an
    //image from a png file on disk
    if(b_cleanupMem){
      if ((buffer) != NULL) {
        delete[] buffer;
        buffer = NULL;
      }

      if (metadata !=NULL) {
        delete metadata;
        metadata = NULL;
      }

      for (list<PngDataChunk*>::iterator it = pdcList.begin(); it != pdcList.end(); it++){
        if ((*it)!=NULL) {
          //cout << "deleting PngDataChunk " << endl;
          delete (*it);
          (*it) = NULL;
        }
      }
      pdcList.clear();
    }
  }

  int getNx() {return nx;}
  int getNy() {return ny;}

  void resize(const int nx_new,const int ny_new) {
    nx = nx_new;
    ny = ny_new;
  }

  void setMetadata(ImageMetadata* _metadata) {
    metadata = _metadata;
  }

  ImageMetadata* getMetadata() {
    return metadata;
  }
  bool isMetadataSet() {
    return (metadata!=NULL);
  }

  void addPngDataChunk(PngDataChunk * pdc) {
    assert(pdc!=NULL);
    for (list<PngDataChunk*>::iterator it = pdcList.begin(); it != pdcList.end(); it++){
      assert(pdc->getName()!=(*it)->getName());
    }
    pdcList.push_back(pdc);
  }

  bool removePngDataChunk(const string &name){
    for (list<PngDataChunk*>::iterator it = pdcList.begin(); it != pdcList.end(); it++) {
       if ((*it)->getName()==name){
         if (b_cleanupMem){
           delete (*it);
         }
         pdcList.erase(it);
         return true;
       }
    }
    return false;
  }

  PngDataChunk * getPngDataChunk(const string &name){
    for (list<PngDataChunk*>::iterator it = pdcList.begin(); it != pdcList.end(); it++) {
       if ((*it)->getName()==name){
         return (*it);
       }
    }
    return NULL;
  }

  void png_write_cascade(png_structp png_ptr, png_infop info_ptr) {
    //  Prepares image and chunk data in to png file format
    //  and writes data to png_ptr which could be a
    //  file system file or memory buffer.
    png_set_IHDR(png_ptr, info_ptr,
                 nx, ny,             // width, height
                 8,                  // bits per pixel -- 16 does not work with blockbuster
                 PNG_COLOR_TYPE_RGB, // non-alpha options are PNG_COLOR_TYPE_RGB, PNG_COLOR_TYPE_GRAY,
                 PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);

#ifdef PNG_WRITE_UNKNOWN_CHUNKS_SUPPORTED
    //png_unknown_chunkp unk_chunk;
    if (pdcList.size()>0) {
      //unk_chunk = new png_unkown_chunk[pdcList.size()];
      png_unknown_chunk unk_chunk[pdcList.size()];
      int ii=0;
      for (list<PngDataChunk*>::iterator it = pdcList.begin(); it != pdcList.end(); it++) {
        unk_chunk[ii].name[0] = (*it)->getName(0);
        unk_chunk[ii].name[1] = (*it)->getName(1);
        unk_chunk[ii].name[2] = (*it)->getName(2);
        unk_chunk[ii].name[3] = (*it)->getName(3);
        unk_chunk[ii].name[4] = '\0';
        unk_chunk[ii].data = (*it)->getData();
        unk_chunk[ii].size = (*it)->getSize();
#ifndef PNG_AFTER_IDAT
        unk_chunk[ii].location = 0;
        cout << "Warning: cannot write depth png metadata " << endl;
        cout << "  try re-compiling with -DPNG_INTERNAL to "<< endl;
        cout << "  enable additional png.h definitions " << endl;
#else
        unk_chunk[ii].location = PNG_AFTER_IDAT;
#endif
        ++ii;
      }

      png_set_unknown_chunks(png_ptr, info_ptr,
                             unk_chunk, pdcList.size());
      /* Always handle all unknown chunks */
      png_set_keep_unknown_chunks(png_ptr, PNG_HANDLE_CHUNK_ALWAYS, NULL, 0);
# if PNG_LIBPNG_VER < 10600
      /* Deal with unknown chunk location bug in 1.5.x and earlier */
      for (int ichunk=0; ichunk<pdcList.size(); ichunk++){
        png_set_unknown_chunk_location(png_ptr, info_ptr, ichunk, PNG_AFTER_IDAT);
      }
# endif
# if PNG_LIBPNG_VER < 10500
      /* PNG_AFTER_IDAT writes two copies of the chunk prior to libpng-1.5.0,
       * one before IDAT and another after IDAT, so don't use it; only use
       * PNG_HAVE_IHDR location.  This call resets the location previously
       * set by assignment and png_set_unknown_chunk_location().
       */
      for (int ichunk=0; ichunk<pdcList.size(); ichunk++){
        png_set_unknown_chunk_location(png_ptr, info_ptr, ichunk, PNG_HAVE_IHDR);
      }
# endif
    }
#else
    cout << "Warning: Unknown PNG chunks not supported" << endl;
#endif

    //Handles the case where custom chunk can't be written after IDAT
    png_write_info_before_PLTE(png_ptr,info_ptr);

    if(metadata!=NULL) {
      char keybuff[1024];
      //text values may be longer than 1024 so don't use a presized buffer

      //tEXt metadata
      vector<string> textMetadata;
      metadata->buildTextMetadataVector(textMetadata);
      assert(textMetadata.size() % 2 == 0);
      for (unsigned int im = 0; im < textMetadata.size()/2; im++) {
        sprintf(keybuff, "%s", textMetadata.at(im*2).c_str());

        png_text metatext;
        metatext.compression = PNG_TEXT_COMPRESSION_NONE;
        metatext.key = keybuff;
        metatext.text = strdup(textMetadata.at(im*2 + 1).c_str());
        png_set_text(png_ptr, info_ptr, &metatext, 1);
        free(metatext.text);
      }

      //zTXt Metadata (Json Metadata for WebUI)
      textMetadata.clear();
      metadata->buildZTextMetadataVector(textMetadata);
      assert(textMetadata.size() % 2 == 0);
      for (unsigned int im = 0; im < textMetadata.size()/2; im++) {
        sprintf(keybuff, "%s", textMetadata.at(im*2).c_str());

        png_text zmetatext;
        zmetatext.compression = PNG_TEXT_COMPRESSION_zTXt;
        zmetatext.key = keybuff;
        zmetatext.text =  strdup(textMetadata.at(im*2 + 1).c_str());
        png_set_text(png_ptr, info_ptr, &zmetatext, 1);
        free(zmetatext.text);
      }

      // Some bits per pixel notes: 16 does not work with blockbuster, and there are also
      // issues with PNG_COLOR_TYPE_GRAY interpretation, so stick to 8 and PNG_COLOR_TYPE_RGB
      // for now. Note that if you do use 16, pay attention to MSB/LSB order. Endien is
      // flipped on my linux workstation...

      png_write_info(png_ptr, info_ptr);
    }

    // set up row pointers to point into the raw image data in buffer...
    png_byte * row_pointers[ny];
    for (int i = 0; i < ny; ++i) row_pointers[i] = (png_byte*)(buffer + i*nx);

    png_write_image(png_ptr, row_pointers);

    png_write_end(png_ptr, info_ptr);

    png_destroy_write_struct(&png_ptr, &info_ptr);

  }

  void set(const int i, const int j, const unsigned char r, const unsigned char g, const unsigned char b) {
    assert((i>=0)&(i<nx));
    assert((j>=0)&(j<ny));
    assert(buffer);
    buffer[i+j*nx][0] = r;
    buffer[i+j*nx][1] = g;
    buffer[i+j*nx][2] = b;
  }


  void get(unsigned char &r, unsigned char &g, unsigned char &b,const int i,const int j) {
    assert((i>=0)&(i<nx));
    assert((j>=0)&(j<ny));
    assert(buffer);
    r = buffer[i+j*nx][0];
    g = buffer[i+j*nx][1];
    b = buffer[i+j*nx][2];
  }

  //write png data to a file
  void write(const char * filename) {
    // this version writes some metadata into the comment part of the png
    // The metadata char should contain two lines separated by a "\n", the
    // first line will become the key and the second the text of png_text

    // note: we completely skip any error handling treatment here for simplicity.
    FILE * fp = fopen(filename, "wb");
    if (!fp) {
      cout << "Warning: could not open png file: " << filename << endl;
      return;
    }

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, (png_voidp)NULL, NULL, NULL);
    if (!png_ptr) {
      fclose(fp);
      cout << "Warning: could not create png_ptr" << endl;
      return;
    }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
      fclose(fp);
      png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
      cout << "Warning: could not create info_ptr" << endl;
      return;
    }

    png_init_io(png_ptr, fp);

    png_write_cascade(png_ptr,info_ptr);

    fclose(fp);
    // leave the image data in buffer to allow the calling process to change it and/or
    // write another image. It is deleted in the destructor.
  }


  //TODO
  //PngImage.read()
  //An image's metadata and optionally it's pixel data can be read.
  //For now all allowed formats will be converted to RGB.
  //PNG_COLOR_TYPE_GRAY_ALPHA will be written over a bahama blue
  //background
  void read(const char * filename, const bool withPixelData=true) {
    b_cleanupMem = true;

    assert ( nx == 0 );
    assert ( ny == 0 );
    assert ( buffer == NULL);
    assert ( pdcList.size() == 0);
    assert ( metadata == NULL );

    // FILE * fp = NULL;
    // const int file_err = MiscUtils::openFile(&fp,string(filename),"rb");
    // if (file_err != 0) return;

    FILE * fp = fopen(filename, "rb");
    if (!fp) {
      cout << "Warning: could not open png file: " << filename << endl;
      return;
    }

    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, (png_voidp) NULL, NULL, NULL);
    if (!png_ptr) {
      fclose(fp);
      cout << "Warning: could not create png_ptr" << endl;
      return;
    }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
      assert(0);
    }

    png_init_io(png_ptr, fp);

    //read the custom chunks
    png_set_keep_unknown_chunks(png_ptr, PNG_HANDLE_CHUNK_ALWAYS, NULL, 0);

    //  //ignore the custom chunks
    //  png_set_keep_unknown_chunks(png_ptr, PNG_HANDLE_CHUNK_NEVER, NULL, 0);

    // read the png characteristics
    png_byte bit_depth;
    png_read_info(png_ptr, info_ptr);
    nx                      = png_get_image_width(png_ptr, info_ptr);
    ny                      = png_get_image_height(png_ptr, info_ptr);
    png_byte  my_color_prec = png_get_color_type(png_ptr, info_ptr);
    bit_depth               = png_get_bit_depth(png_ptr, info_ptr);

    if (bit_depth == 16) png_set_strip_16(png_ptr);

    if (my_color_prec == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png_ptr);

    if (my_color_prec == PNG_COLOR_TYPE_GRAY_ALPHA){
      png_set_gray_to_rgb(png_ptr);
      //add bahama blue to background of gray-alpha type images
      //gray, red, green, blue
      png_color_16 bahama_blue = {(png_byte)0, (png_uint_16)73, (png_uint_16)175, (png_uint_16)205, (png_uint_16)0};
      png_set_background(png_ptr, &bahama_blue, PNG_BACKGROUND_GAMMA_SCREEN, 0, 1);
    }

    png_read_update_info(png_ptr, info_ptr);

    if (my_color_prec != PNG_COLOR_TYPE_PALETTE &&
        my_color_prec != PNG_COLOR_TYPE_RGB &&
        my_color_prec != PNG_COLOR_TYPE_GRAY_ALPHA) {
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      fclose(fp);
      cout << "Warning: unsupported png color type for file: " << filename << endl;
      return;
    }

    bool verbose = false;
    png_textp text_ptr;
    int num_text  = png_get_text(png_ptr, info_ptr,&text_ptr, NULL);

    if (metadata==NULL){
      metadata = new ImageMetadata();
    }

    for (int itext = 0; itext<num_text; itext++) {
      metadata->setTextFromImage(text_ptr[itext].key,text_ptr[itext].text,verbose);
    }

    //Leave image as type PNG_COLOR_TYPE_RGB if it is used later for writing
    my_color_prec = PNG_COLOR_TYPE_RGB;

    if (withPixelData) { //read the image data, and custom chunks (TBD)
      buffer = new unsigned char[nx*ny][3];

      png_byte * row_pointers[ny];
      for (int i = 0; i < ny; ++i) row_pointers[i] = (png_byte*)(buffer + i*nx);

      png_read_image(png_ptr, row_pointers);
      png_read_end(png_ptr, info_ptr);

      //process cascade custom chunks (depending on version of libpng used for writing,
      //the chunks can be before or after the image data)

      png_unknown_chunkp unknowns;
      int num_unknown_chunks = png_get_unknown_chunks(png_ptr,info_ptr, &unknowns);

      //TODO
      //should set this up to read only compressed data or if requested
      //uncompress it.
      for (int i = 0; i < num_unknown_chunks; i++) {
        string chunkName = reinterpret_cast<char*>(unknowns[i].name);
        //uLongf zlib_size = nx*ny*PngDataChunk::MAX_BYTES_PER_PIXEL;
        PngDataChunk * pdc = new PngDataChunk(chunkName, nx, ny, reinterpret_cast<unsigned char*>(unknowns[i].data), unknowns[i].size);
        if (pdc->getData()!=NULL) pdcList.push_back(pdc);  //if read fails just skip adding this chunk data to pdcList
        else delete pdc;
      }
    }
    else {
      //skip buffer allocation
    }
    //cleanup...
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);

    fclose(fp);
    if ( verbose ) cout << "Read in " << filename << endl;

  } //end PngImage.read()

  //convert RGB to YUV420p (useful for animation encoding)
  void rgbToYuv420p(unsigned char* &destination){
    assert(destination);
    unsigned int upos = nx*ny;
    unsigned int vpos = upos + upos / 4;
    unsigned int i = 0;

    for( unsigned int line = 0, nline=ny; line < nline; ++line ){
      if( !(line % 2) ){
        for( unsigned int x=0, xmax=nx; x < xmax; x += 2 ){
          uint8_t r = buffer[i][0];
          uint8_t g = buffer[i][1];
          uint8_t b = buffer[i][2];

          destination[i++] = ((66*r + 129*g + 25*b) >> 8) + 16;

          destination[upos++] = ((-38*r + -74*g + 112*b) >> 8) + 128;
          destination[vpos++] = ((112*r + -94*g + -18*b) >> 8) + 128;

          r = buffer[i][0];
          g = buffer[i][1];
          b = buffer[i][2];

          destination[i++] = ((66*r + 129*g + 25*b) >> 8) + 16;
        }
      }
      else{
        for( unsigned int x = 0,xmax=nx; x < xmax; x += 1 ){
          uint8_t r = buffer[i][0];
          uint8_t g = buffer[i][1];
          uint8_t b = buffer[i][2];

          destination[i++] = ((66*r + 129*g + 25*b) >> 8) + 16;
        }
      }
    }
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
