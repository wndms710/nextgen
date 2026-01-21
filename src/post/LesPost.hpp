#ifndef _LES_POST_HPP_
#define _LES_POST_HPP_

#include "CTI.hpp"
using namespace CTI;

#include "MiscUtils.hpp"
#include "ByteSwap.hpp"

class SerialLesReader {
public:

  static int readCvDN(double*& v,int& ncv,const string& token,const string& filename) {
    
    FILE * fp = fopen(filename.c_str(),"rb");
    if (fp == NULL) {
      cout << "Error: SerialLesReader: cannot open file: " << filename << endl;
      return -1;
    }

    bool byte_swap = false;
    int itmp[2];
    fread(itmp,sizeof(int),2,fp);
    if ((itmp[0] != UGP_IO_MAGIC_NUMBER)&&(itmp[0] != UGP_IO_MAGIC_NUMBER+1)) {
      ByteSwap::byteSwap(itmp,2);
      if ((itmp[0] != UGP_IO_MAGIC_NUMBER)&&(itmp[0] != UGP_IO_MAGIC_NUMBER+1)) {
        cerr << "Error: file does not start as expected. aborting."<< endl;
	fclose(fp);
	return -1;
      }
      cout << "File requires byte swapping." << endl;
      byte_swap = true;
    }
    assert(!byte_swap); // if you hit this, we need to work byte swapping through the routine more carefully
    
    int8 offset = sizeof(int)*2;
    
    bool found = false;
    Header header;
    while (1) {
      
      fseek (fp,offset,SEEK_SET);
      
      fread(&header,sizeof(Header),1,fp);
      //if (byte_swap) byteSwapHeader(&header,1);
      assert(header.skip > 0);
      
      //cout << " > got header: " << header.name << endl;
      
      if (header.id == UGP_IO_EOF) {
	break;
      }
      else if (header.id == UGP_IO_CV_D1) {
	if (token.compare(header.name) == 0) {
	  int8 ncv_global = ByteSwap::getInt8FromLswMswPair(header.idata);
	  assert(ncv_global < TWO_BILLION);
	  //cout << " > ncv_global: " << ncv_global << endl;
	  if (v == NULL) {
	    ncv = ncv_global;
	    v = new double[ncv_global];
	    if (v == NULL) {
	      cout << "memory allocation failed" << endl;
	      return -1;
	    }
	  }
	  else {
	    assert(ncv == ncv_global);
	  }
	  fread(v,sizeof(double),ncv,fp);
	  found = true;
	  break;
	}
      }

      offset += header.skip;
      
    }
    
    if (!found) {
      cout << "Error: could not find CVDN: " << token << endl;
      offset = sizeof(int)*2;
      while (1) {
	fseek (fp,offset,SEEK_SET);
	fread(&header,sizeof(Header),1,fp);
	//if (byte_swap) byteSwapHeader(&header,1);
	assert(header.skip > 0);
	//cout << " > got header: " << header.name << endl;
	if (header.id == UGP_IO_EOF) {
	  break;
	}
	else if (header.id == UGP_IO_CV_D1) {
	  cout << " > " << header.name << endl;
	}
	offset += header.skip;
      }
      throw(-1);
    }
    

    
    fclose(fp);
    return 0;

  }

  static int readCvDN3(double (*&v)[3],int& ncv,const string& token,const string& filename) {
    
    FILE * fp = fopen(filename.c_str(),"rb");
    if (fp == NULL) {
      cout << "Error: SerialLesReader: cannot open file: " << filename << endl;
      return -1;
    }

    bool byte_swap = false;
    int itmp[2];
    fread(itmp,sizeof(int),2,fp);
    if ((itmp[0] != UGP_IO_MAGIC_NUMBER)&&(itmp[0] != UGP_IO_MAGIC_NUMBER+1)) {
      ByteSwap::byteSwap(itmp,2);
      if ((itmp[0] != UGP_IO_MAGIC_NUMBER)&&(itmp[0] != UGP_IO_MAGIC_NUMBER+1)) {
        cerr << "Error: file does not start as expected. aborting."<< endl;
	fclose(fp);
	return -1;
      }
      cout << "File requires byte swapping." << endl;
      byte_swap = true;
    }
    assert(!byte_swap); // if you hit this, we need to work byte swapping through the routine more carefully
    
    int8 offset = sizeof(int)*2;
    bool found = false;
    Header header;
    while (1) {
      
      fseek (fp,offset,SEEK_SET);
      
      fread(&header,sizeof(Header),1,fp);
      //if (byte_swap) byteSwapHeader(&header,1);
      assert(header.skip > 0);
      
      //cout << " > got header: " << header.name << endl;
      
      if (header.id == UGP_IO_EOF) {
	break;
      }
      else if (header.id == UGP_IO_CV_D2) {
	if (token.compare(header.name) == 0) {
	  int8 ncv_global = ByteSwap::getInt8FromLswMswPair(header.idata);
	  assert(ncv_global < TWO_BILLION);
	  //cout << " > ncv_global: " << ncv_global << endl;
	  if (v == NULL) {
	    ncv = ncv_global;
	    v = new double[ncv_global][3];
	    if (v == NULL) {
	      cout << "memory allocation failed" << endl;
	      return -1;
	    }
	  }
	  else {
	    assert(ncv == ncv_global);
	  }
	  fread(v,sizeof(double),ncv*3,fp);
	  found = true;
	  break;
	}
      }
      
      offset += header.skip;
      
    }

    if (!found) {
      cout << "Error: could not find CVDN3: " << token << endl;
      offset = sizeof(int)*2;
      while (1) {
	fseek (fp,offset,SEEK_SET);
	fread(&header,sizeof(Header),1,fp);
	//if (byte_swap) byteSwapHeader(&header,1);
	assert(header.skip > 0);
	//cout << " > got header: " << header.name << endl;
	if (header.id == UGP_IO_EOF) {
	  break;
	}
	else if (header.id == UGP_IO_CV_D2) {
	  cout << " > " << header.name << endl;
	}
	offset += header.skip;
      }
      throw(0);
    }
    
    fclose(fp);
    return 0;

  }
    
private:

  // Disallow creating an instance of this object
  SerialLesReader() {}

};

class LesPost {
public:

  virtual void initialHook() {}
  virtual void temporalHook(const string& snapshot,const int step) {}
  virtual void finalHook() {}
  
  void init() {
    
    initialHook();
    
  }
  
  void run() {

    // look for the first SNAPSHOT/SNAPSHOTS param...
    Param * param = getParam("SNAPSHOT");
    if (param == NULL)  param = getParam("SNAPSHOTS");
    
    // there are two formats supported: 
    //
    // SNAPSHOT NAME <name> RANGE <start> <inc> <end> 
    // 
    // or
    //
    // SNAPSHOT <snapshot1> [<snapshot2> ...]
    //
    
    if (param) {

      if (param->getString(0) == "NAME") {
	
	const string name = param->getString(1);
	assert(param->getString(2) == "RANGE");
	const int start = param->getInt(3);
	const int inc   = param->getInt(4);
	const int end   = param->getInt(5);
	
	for (int step = start; step <= end; step += inc) {
	  
	  char filename[128];
	  MiscUtils::buildIndexedFilename(filename,name.c_str(),step,"sles");
	  
	  temporalHook(filename,step);

	}

      }
      
      else {

	// we assume these are just a set of snapshot filenames...
	
	for (int iarg = 0; iarg < param->size(); ++iarg) {
	  
	  const string name = param->getString(iarg);
	  const int step = MiscUtils::getIndexFromFilename(name);
	  
	  temporalHook(name,step);
	  
	}

      }

    }

    finalHook();
    
  }
  
};

#endif
