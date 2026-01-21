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

using namespace std;

#include "tomcrypt.hpp"
#include "RestartHashUtilities.hpp"

void my_fseek(FILE * fp, int8 offset) {

  // not very efficient, but works on BG
  int8 cur_offset = offset;
  int8 max_offset = 1000000000ll;

  fseek(fp,0,SEEK_SET);
  while (cur_offset > max_offset) {
    fseek(fp,max_offset,SEEK_CUR);
    cur_offset -= max_offset;
  }
  fseek(fp,cur_offset,SEEK_CUR);
}

string getFaZoneNameForKind(const int kind) {
  switch(kind) {
  case FA_ZONE_BOUNDARY:
    return("BOUNDARY");
  case FA_ZONE_PERIODIC_CART:
    return("PERIODIC_CART");
  case FA_ZONE_PERIODIC_CYL_X:
    return("PERIODIC_CYL_X");
  case FA_ZONE_PERIODIC_CYL_Y:
    return("PERIODIC_CYL_Y");
  case FA_ZONE_PERIODIC_CYL_Z:
    return("PERIODIC_CYL_Z");
  case FA_ZONE_PERIODIC_INTERNAL:
    return("PERIODIC_INTERNAL");
  case FA_ZONE_INTERNAL:
    return("INTERNAL");
  default:
    return("UKNOWN");
  }
}

double atod(const string& var) {
  std::istringstream iss(var);
  double value;
  if ((iss >> value >> std::dec).fail()) {
    //cerr << "\nError: cannot evaluate string to double: \"" << var << "\"" << endl;
    throw(-1);
  }
  return value;
}

int8 atoi8(const string& var) {
  std::istringstream iss(var);
  int8 value;
  if ((iss >> value >> std::dec).fail()) {
    //cerr << "\nError: cannot evaluate string to int8: \"" << var << "\"" << endl;
    throw(-1);
  }
  return value;
}

string getCvZoneNameForKind(const int kind) {
  switch(kind) {
  case CV_ZONE_FLUID:
    return("FLUID");
  default:
    return("UKNOWN");
  }
}

enum LesinfoOptions {
  LESINFO_NO_OPTIONS,
  LESINFO_INPUT_PARAMS,
  LESINFO_JOURNAL,
  LESINFO_MATCH,
  LESINFO_EXTRACT,
  LESINFO_EXTRACT_ALL,
  LESINFO_EXTRACT_TO_FILE,
  LESINFO_CHECK,
  LESINFO_PBIN,
  LESINFO_SBIN,
  LESINFO_ID
};

enum ImageFormat {
  PPM_IMAGE,
  PNG_IMAGE
};

void checkContiguous(FILE * fp,const int n) {
  cout << "# Reading and checking " << n << " contiguous ints";
  const int read_size = 10000000;
  int *buf = new int[min(n,read_size)];
  int8 index0 = 0;
  while (index0 < n) {
    cout << ".";
    cout.flush();
    const int this_read_size = min(read_size,int(n-index0));
    fread(buf,sizeof(int),this_read_size,fp);
    for (int i = 0; i < this_read_size; ++i) {
      if (buf[i] != index0+i) {
	cerr << "\nError: value at position " << index0+i << " does not match: " << buf[i];
	throw(-1);
      }
    }
    index0 += this_read_size;
  }
  delete[] buf;
  cout << "OK" << endl;
}

int main(int argc,char * argv[]) {

  using namespace ByteSwap;

  try {

    if (sizeof(int8) != 8) {
      cerr << "\n\nError: lesinfo is not compiled with the correct int8 size: " << sizeof(int8) << endl;
      cerr << "try compiling with -D INT8_IS_LONG_INT added to your CXXFLAGS and CXXFLAGS_NO_MPI in Makefile.in\n\n" << endl;
      throw(-1);
    }

    assert( sizeof(Header) == 256 );

    if (argc < 2) {
      cerr << "\n\nlesinfo usage: \n\n$ lesinfo [option] <les-filename>" << endl;
      cerr << "\noption is one of:" << endl;
      cerr << " -i                                   prints input file parameters (if present)" << endl;
      cerr << " -j                                   prints journal (if present)" << endl;
      cerr << " -check                               reads and checks contiguous index in no_check,fa_check,cv_check" << endl;
      cerr << " -match <VAR> <X> <Y> <Z>             reports global index of VAR nearest point (X,Y,Z)" << endl;
      cerr << " -extract <VAR> <global_index>        reports value of VAR at global_index" << endl;
      cerr << " -extract-fast <VAR> <global_index>   same as above, but fast version for identical snapshot files" << endl;
      cerr << " -extract-all <NV> <V1> <V2> [...]    extracts all named data to stdout" << endl;
      cerr << " -extract-to-file <V1> [<V2>...]      extracts named data to sles file for INTERP_FROM" << endl;
      cerr << " -pbin <filename>                     extracts coordinate data to filename" << endl;
      cerr << " -sbin <filename>                     extracts surface from mles to filename" << endl;
      cerr << " -id                                  checks for consistency between an mles and one or more sles hash ids" << endl;
      cerr << "\n\n" << endl;
      throw(-1);
    }

    bool verbose = true;

    //for images
    int wireframe_count = 0;
    int imageFormat = PPM_IMAGE;
#ifdef WITH_PNG
    imageFormat = PNG_IMAGE;
#endif

    string name;
    double xp[3];
    int8 global_index_array[20];
    int global_index_array_size = 0;
    int8 fast_offset = -1;
    int fast_type;
    set<string> nameSet;
    int8 offset_surface = -1;

    int8 offset_vol = -1;

    LesinfoOptions option = LESINFO_NO_OPTIONS;
    int iarg0 = 1;
    while ( (iarg0 < argc) && (strncmp(argv[iarg0],"-",1) == 0) ) {
      if (strcmp(argv[iarg0],"-i") == 0) {
	verbose = false;
	if (option != LESINFO_NO_OPTIONS) {
	  cerr << "\n\nError: use only one option at a time.\n\n" << endl;
	  throw(-1);
	}
	option = LESINFO_INPUT_PARAMS;
	++iarg0;
      }
      else if (strcmp(argv[iarg0],"-j") == 0) {
	verbose = false;
	if (option != LESINFO_NO_OPTIONS) {
	  cerr << "\n\nError: use only one option at a time.\n\n" << endl;
	  throw(-1);
	}
	option = LESINFO_JOURNAL;
	++iarg0;
      }
      else if (strcmp(argv[iarg0],"-check") == 0) {
        if (option != LESINFO_NO_OPTIONS) {
          cerr << "\n\nError: use only one option at a time.\n\n" << endl;
          throw(-1);
        }
        option = LESINFO_CHECK;
        ++iarg0;
      }
      else if (strcmp(argv[iarg0],"-match") == 0) {
        verbose = false;
        if (option != LESINFO_NO_OPTIONS) {
          cerr << "\n\nError: use only one option at a time.\n\n" << endl;
          throw(-1);
        }
        option = LESINFO_MATCH;
        ++iarg0;
	try {
	  if (argc < iarg0+4)
	    throw(-1);
	  name = argv[iarg0++];
	  xp[0] = atod(argv[iarg0++]);
	  xp[1] = atod(argv[iarg0++]);
	  xp[2] = atod(argv[iarg0++]);
	}
	catch(...) {
          cerr << "\nError: problem parsing match option.\n\nUsage:\n\n -match <VAR> <X> <Y> <Z> reports global index of VAR nearest point (X,Y,Z)\n" << endl;
          throw(-1);
	}
	cout << "# looking for global index of data " << name << " nearest point " << COUT_VEC(xp) << endl;
      }
      else if ( (strcmp(argv[iarg0],"-extract") == 0) || (strcmp(argv[iarg0],"-extract-fast") == 0) ) {
        verbose = false;
        if (option != LESINFO_NO_OPTIONS) {
          cerr << "\n\nError: use only one option at a time.\n\n" << endl;
          throw(-1);
        }
	option = LESINFO_EXTRACT;
	if (strcmp(argv[iarg0],"-extract-fast") == 0) {
	  assert(fast_offset == -1);
	  fast_offset = -2;
	}
	++iarg0;
	name = argv[iarg0++];
	global_index_array_size = 0;
	try {
	  while (1) {
	    int8 this_global_index = atoi8(argv[iarg0]);
	    assert(global_index_array_size < 20);
	    global_index_array[global_index_array_size] = this_global_index;
	    ++global_index_array_size;
	    ++iarg0;
	  }
	}
	catch(...) {
	}
	cout << "# extracting data " << name << " at global index";
	for (int i = 0; i < global_index_array_size; ++i)
	  cout << " " << global_index_array[i];
	cout << endl;
      }
      else if (strcmp(argv[iarg0],"-extract-all") == 0) {
        verbose = false;
        if (option != LESINFO_NO_OPTIONS) {
          cerr << "\n\nError: use only one option at a time.\n\n" << endl;
          throw(-1);
        }
        option = LESINFO_EXTRACT_ALL;
        ++iarg0;
	const int nvar = atoi(argv[iarg0++]);
	assert((nvar > 0)&&(nvar < 20)); // 20 seems like a lot
	cout << "# extract-all vars:";
	for (int ivar = 0; ivar < nvar; ++ivar) {
	  cout << " " << argv[iarg0];
	  nameSet.insert(argv[iarg0++]);
	}
	cout << endl;
      }
      else if (strcmp(argv[iarg0],"-extract-to-file") == 0) {
        verbose = false;
        if (option != LESINFO_NO_OPTIONS) {
          cerr << "\n\nError: use only one option at a time.\n\n" << endl;
          throw(-1);
        }
        option = LESINFO_EXTRACT_TO_FILE;
        ++iarg0;
	// remaining args should be var names...
	cout << "# extract-to-file:";
	while (iarg0 < argc-1) { // leave 1 arg for the *les file
	  cout << " " << argv[iarg0];
	  nameSet.insert(argv[iarg0++]);
	}
	cout << endl;
      }
      else if (strcmp(argv[iarg0],"-pbin") == 0) {
        //verbose = false;
        if (option != LESINFO_NO_OPTIONS) {
          cerr << "\n\nError: use only one option at a time.\n\n" << endl;
          throw(-1);
        }
        option = LESINFO_PBIN;
        ++iarg0;
	name = argv[iarg0++];
	cout << "# pbin: " << name << endl;
        assert(offset_vol==-1);
      }
      else if (strcmp(argv[iarg0],"-sbin") == 0) {
        //verbose = false;
        if (option != LESINFO_NO_OPTIONS) {
          cerr << "\n\nError: use only one option at a time.\n\n" << endl;
          throw(-1);
        }
        option = LESINFO_SBIN;
        ++iarg0;
	name = argv[iarg0++];
	cout << "# sbin: " << name << endl;
        assert(offset_surface==-1);
      }
      else if (strcmp(argv[iarg0],"-id") == 0) {
        verbose = false;
        if (option != LESINFO_NO_OPTIONS) {
          cerr << "\n\nError: use only one option at a time.\n\n" << endl;
          throw(-1);
        }
        option = LESINFO_ID;
        ++iarg0;
	cout << "# checking for matching mles/sles files..." << endl;
      }
      else {
	cerr << "\n\nError: unrecognized option: " << argv[iarg0] << "\n\n" << endl;
	throw(-1);
      }

    }

    // remaining args should be *.les files...

    if (iarg0 == argc) {
      cerr << "\nError: missing les files.\n" << endl;
      throw(-1);
    }

    for (int iarg = iarg0; iarg < argc; ++iarg) {
      // FILE * fp = NULL;
      // const int file_err = MiscUtils::openFile(&fp,string(argv[iarg]),"rb");
      // if (file_err != 0) throw(file_err);

      FILE *fp_extract = NULL;

      FILE * fp = fopen(argv[iarg],"rb");
      if (fp == NULL) {
	cerr << "\n\nError: could not open file: " << argv[iarg] << "\n\n" << endl;
	throw(-1);
      }

      // the fast_offset mode is triggered by a positive fast_offset...
      // this does not make much difference it seems.

      if (fast_offset >= 0) {
	switch (fast_type) {
	case UGP_IO_CV_F1:
	  cout << "# extract scalar CV_F1: " << name << " from file " << argv[iarg];
	  for (int i = 0; i < global_index_array_size; ++i) {
	    my_fseek(fp,fast_offset+global_index_array[i]*sizeof(float));
	    float var;
	    fread(&var,sizeof(float),1,fp);
	    cout << " " << var;
	  }
	  cout << endl;
	  break;
	default:
	  assert(0);
	}

	fclose(fp);
	continue;
      }


      switch (option) {
      case LESINFO_JOURNAL:
	if (iarg == iarg0)
	  cout << "# journal from les file: " << argv[iarg] << endl;
	else
	  cout << "\n# journal from les file: " << argv[iarg] << endl;
	break;
      case LESINFO_INPUT_PARAMS:
	if (iarg == iarg0)
	  cout << "# input params from les file: " << argv[iarg] << endl;
	else
	  cout << "\n# input params from les file: " << argv[iarg] << endl;
	break;
      case LESINFO_NO_OPTIONS:
	if (iarg == iarg0)
	  cout << "File: " << argv[iarg] << endl;
	else
	  cout << "\nFile: " << argv[iarg] << endl;
	break;
      case LESINFO_MATCH:
      case LESINFO_EXTRACT:
      case LESINFO_EXTRACT_ALL:
      case LESINFO_EXTRACT_TO_FILE:
      case LESINFO_CHECK:
      case LESINFO_PBIN:
      case LESINFO_SBIN:
      case LESINFO_ID:
        break;
      default:
	cerr << "Unrecognized LESINFO option";
      }


      int byte_swap = 0;
      int itmp[2];
      fread(itmp,sizeof(int),2,fp);
      if ((itmp[0] != UGP_IO_MAGIC_NUMBER)&&(itmp[0] != UGP_IO_MAGIC_NUMBER+1)&&(itmp[0] != MOVE_IO_MAGIC_NUMBER)) {
	byteSwap(itmp,2);
	if ((itmp[0] != UGP_IO_MAGIC_NUMBER)&&(itmp[0] != UGP_IO_MAGIC_NUMBER+1)&&(itmp[0] != MOVE_IO_MAGIC_NUMBER)) {
	  cerr << "Error: file does not start as expected. aborting."<< endl;
	  throw(-1);
	}
	if (verbose) cout << "File requires byte swapping." << endl;
	byte_swap = 1;
      }
      
      if (itmp[0] == UGP_IO_MAGIC_NUMBER+1) {
	if (verbose) cout << "looks like a snapshot data file" << endl;
      }
      else if (itmp[0] == MOVE_IO_MAGIC_NUMBER+1) {
	if (verbose) cout << "looks like a moving solver data file" << endl;
      }
      
      const int io_version = itmp[1];
      if (verbose) cout << "io version: " << io_version << endl;

      // surface extraction...
      vector<string> zoneVec;
      bool b_periodic = false;

      // initialize the offset...

      int8 offset = sizeof(int)*2;
      Header header;

      int done = 0;
      while (done != 1) {

	my_fseek(fp,offset);

	fread(&header,sizeof(Header),1,fp);
	if (byte_swap)
	  byteSwapHeader(&header,1);

	if (header.skip <= 0) {
	  cerr << "Error: header has invalid skip: " << header.skip << " at file offset: " << offset << endl;
	  cerr << " > header name =  \"" << header.name << "\"" << endl;
	  cerr << " > header id = " << header.id << endl;
	  cerr << " > header.skip = " << header.skip << endl;
	  cerr << " > header idata: " << header.idata[0] << " " << header.idata[1] << " " << header.idata[2] << "..." << endl;
	  cerr << " > header rdata: " << header.rdata[0] << " " << header.rdata[1] << " " << header.rdata[2] << "..." << endl;
	  throw(-1);
	}

	switch (header.id) {
	case UGP_IO_TIMESTAMP:
	  {
	    struct tm now;
	    now.tm_sec   = header.idata[0]; // seconds of minutes from 0 to 61
	    now.tm_min   = header.idata[1]; // minutes of hour from 0 to 59
	    now.tm_hour  = header.idata[2]; // hours of day from 0 to 24
	    now.tm_mday  = header.idata[3]; // day of month from 1 to 31
	    now.tm_mon   = header.idata[4]; // month of year from 0 to 11
	    now.tm_year  = header.idata[5]; // year since 1900
	    now.tm_wday  = header.idata[6]; // days since sunday
	    now.tm_yday  = header.idata[7]; // days since January 1st
	    now.tm_isdst = header.idata[8]; // hours of daylight savings time
	    char timestamp[24];
	    strftime(timestamp,24,"%y%m%d:%H%M%S",&now);
	    switch (option) {
	    case LESINFO_INPUT_PARAMS:
	    case LESINFO_JOURNAL:
	      cout << "# timestamp: " << timestamp << endl;
	      break;
	    default:
	      if (verbose) cout << "timestamp: " << timestamp << endl;
	    }
	  }
	  break;
	case UGP_IO_PARAMS:
	  if (option == LESINFO_INPUT_PARAMS) {
	    char * cbuf = new char[header.idata[0]+1];
	    fread(cbuf,sizeof(char),header.idata[0],fp);
	    cbuf[header.idata[0]] = '\0';
	    cout << cbuf; cout.flush();
	    delete[] cbuf;
	    done = 1;
	  }
	  else {
	    if (verbose) cout << "params (use option -i for details)" << endl;
	  }
	  break;
	case UGP_IO_JOURNAL:
	  if (option == LESINFO_JOURNAL) {
	    char * cbuf = new char[header.idata[0]+1];
	    fread(cbuf,sizeof(char),header.idata[0],fp);
	    cbuf[header.idata[0]] = '\0';
	    cout << cbuf; cout.flush();
	    delete[] cbuf;
	    done = 1;
	  }
	  else {
	    if (verbose) cout << "journal (use option -j for details)" << endl;
	  }
	  break;
        case UGP_IO_WIREFRAME:
          {
            if ( verbose) cout << header.name << " ; found wireframe" << endl;
            break;
          }
        case UGP_IO_HASHIDS:
          {
            if (verbose){
              RestartHashUtilities::readHashFromRestartHeader_s(fp,header);
              cout << "Restart Hash Id: " << RestartHashUtilities::mlesHash << endl;
              cout << "Parent  Hash Id: " << RestartHashUtilities::myHash << endl;
              RestartHashUtilities::clearHashes();
            }
            else if (option==LESINFO_ID){
              if (itmp[0] == UGP_IO_MAGIC_NUMBER+1){
                RestartHashUtilities::slesReadHashes_s(fp,header);
                if (RestartHashUtilities::mlesHash.getLength()>0){//mles hash has been set
                  if (RestartHashUtilities::slesConsistencyCheck()>0) { //not matching
                    cout << "WARNING: " << argv[iarg] << " does not match the mles file, expecting";
                    cout << " mles hash id " << RestartHashUtilities::myHash << endl;
                    RestartHashUtilities::slesHash.clear();
                    RestartHashUtilities::myHash.clear();
                  }
                  else{ //matching
                    cout << argv[iarg] << " matches mles id" << endl;
                    RestartHashUtilities::slesHash.clear();
                  }
                }
                else{ //no mles hash was set or found
                  cout << "WARNING: no mles consistency check is possible, an mles file was either not specified or contained no id" << endl;
                  iarg = argc; //don't read any more files
                }
              }
              else{
                RestartHashUtilities::clearHashes();
                RestartHashUtilities::mlesReadHashes_s(fp,header);
                cout << argv[iarg] << " has hash id: " << RestartHashUtilities::mlesHash << endl;
              }
              done = 1;
            }
            break;
          }
	case UGP_IO_NO_FA_CV_COUNTS:
	  if (io_version == 1) {
	    if (verbose) cout << "nno: " << header.idata[0] << endl;
	    if (verbose) cout << "nfa: " << header.idata[1] << endl;
	    if (verbose) cout << "ncv: " << header.idata[2] << endl;
	  }
	  else {
	    assert((io_version > 1)&&(io_version <=4));
	    // io_version 2 and higher allows for counts > 2B
	    if (verbose) cout << "nno: " << getInt8FromLswMswPair(header.idata+0) << endl;
	    if (verbose) cout << "nfa: " << getInt8FromLswMswPair(header.idata+2) << endl;
	    if (verbose) cout << "ncv: " << getInt8FromLswMswPair(header.idata+4) << endl;
	  }
	  break;
	case UGP_IO_NO_FA_BF_CV_COUNTS:
          if (verbose) {
            cout << "nno: " << getInt8FromLswMswPair(header.idata+0) << endl;
            cout << "nbf: " << getInt8FromLswMswPair(header.idata+2) << endl;
            cout << "nfa: " << getInt8FromLswMswPair(header.idata+4) << endl;
            cout << "ncv: " << getInt8FromLswMswPair(header.idata+6) << endl;
          }
          {
            int8 nno_p = getInt8FromLswMswPair(header.idata+8);
            if (nno_p > 0)
              b_periodic = true;
	  }
	  break;
	case UGP_IO_MOVING_HEADER:
          if (verbose) {
            cout << "moving solver header:" << endl;
            cout << " > ncv_global: " << header.ui8data[0] << endl;
            cout << " > step:       " << header.idata[0] << endl;
            cout << " > time:       " << header.rdata[0] << endl;
            cout << " > dt:         " << header.rdata[1] << endl;
          }
          break;
	case UGP_IO_NO_CHECK:
	  if (verbose) cout << "no_check nno: " << header.idata[0] << endl;
	  if (option == LESINFO_CHECK) checkContiguous(fp,header.idata[0]);
	  break;
	case UGP_IO_FA_CHECK:
	  if (verbose) cout << "fa_check nfa: " << header.idata[0] << endl;
	  if (option == LESINFO_CHECK) checkContiguous(fp,header.idata[0]);
	  break;
	case UGP_IO_CV_CHECK:
	  if (verbose) cout << "cv_check ncv: " << header.idata[0] << endl;
	  if (option == LESINFO_CHECK) checkContiguous(fp,header.idata[0]);
	  break;
        case UGP_IO_NO_CHECK_INT8:
          if (verbose) cout << "no_check nno: " << getInt8FromLswMswPair(header.idata) << endl;
          break;
        case UGP_IO_FA_CHECK_INT8:
          if (verbose) cout << "fa_check nfa: " << getInt8FromLswMswPair(header.idata) << endl;
          break;
        case UGP_IO_CV_CHECK_INT8:
          if (verbose) cout << "cv_check ncv: " << getInt8FromLswMswPair(header.idata) << endl;
          break;
	case UGP_IO_NOOFA_I_AND_V:
        case UGP_IO_NOOFA_I_AND_V_V3_INT8:
	  if (verbose) cout << "noofa_i/v: " << header.idata[0] << " " << getInt8FromLswMswPair(header.idata+1) << endl;
	  break;
        case UGP_IO_NOOFA_I_AND_V_INT8:
	  if (verbose) cout << "noofa_i/v: " << getInt8FromLswMswPair(header.idata+0) << " " << getInt8FromLswMswPair(header.idata+2) << endl;
	  break;
	case UGP_IO_CVOFA:
        case UGP_IO_CVOFA_V3_INT8:
	  if (verbose) cout << "cvofa: " << header.idata[0] << " " << header.idata[1] << endl;
	  break;
        case UGP_IO_CVOFA_INT8:
	  if (verbose) cout << "cvofa: " << getInt8FromLswMswPair(header.idata+0) << " " << header.idata[2] << endl;
	  break;
	case UGP_IO_DATA:
	  if (verbose) cout << "start of data records" << endl;
	  break;
	case UGP_IO_CV_PART:
	  if (verbose) cout << "cv partition size: " << header.idata[1] << endl;
	  break;
        case UGP_IO_BIT_TRANSFORM_DATA:
          if ( verbose) cout << "found periodic bit : " << header.name << endl;
          break;
	case UGP_IO_FA_ZONE_HEADER:
	  if (io_version == 1) {
	    if (verbose) cout << " > face zone \"" << header.name << "\" " << header.idata[2]-header.idata[1]+1 << " faces" << endl;
	  }
	  else {
	    assert((io_version > 1)&&(io_version <= 4));
	    //if (verbose) cout << "face zone \"" << header.name << "\", kind=" << getFaZoneNameForKind(header.idata[0]) << ", index=" << header.idata[1] << ", nfa=" << header.idata[2] << endl;
	    if (verbose) {
	      if (header.idata[0] == FA_ZONE_BOUNDARY) {
		// these records can have geom...
		const double area_global            = header.rdata[0];
		const double area_over_delta_global = header.rdata[1];
		double n_global[3];
		n_global[0]            = header.rdata[2];
		n_global[1]            = header.rdata[3];
		n_global[2]            = header.rdata[4];
		double x_global[3];
		x_global[0]            = header.rdata[5];
		x_global[1]            = header.rdata[6];
		x_global[2]            = header.rdata[7];
		cout << " > face zone \"" << header.name << "\", kind=" << getFaZoneNameForKind(header.idata[0]) << ", index=" << header.idata[1] <<
		  ", area: " << area_global <<
		  ", n: " << COUT_VEC(n_global) <<
		  ", x: " << COUT_VEC(x_global) <<
		  ", mag(n)/area: " << MAG(n_global)/area_global << endl;
	      }
	      else {
		cout << " > face zone \"" << header.name << "\", kind=" << getFaZoneNameForKind(header.idata[0]) << ", index=" << header.idata[1] << endl;
	      }
	    }
	  }
	  break;
        case UGP_IO_BF_ZONE_HEADER:

          if ( option == LESINFO_SBIN ) {
            zoneVec.push_back(header.name);
          }

          if ( verbose) {
	    // counts...
	    const int8 nbf_global   = ByteSwap::getInt8FromLswMswPair(header.idata+2);
	    const int8 noobf_global = ByteSwap::getInt8FromLswMswPair(header.idata+4);
	    const int8 ibf_f_global = ByteSwap::getInt8FromLswMswPair(header.idata+6);
	    const int8 nob_f_global = ByteSwap::getInt8FromLswMswPair(header.idata+8);
	    const int8 stobf_global = ByteSwap::getInt8FromLswMswPair(header.idata+10);
	    const int8 sob_f_global = ByteSwap::getInt8FromLswMswPair(header.idata+12);
	    // and geom...
	    const double area_global            = header.rdata[0];
	    const double area_over_delta_global = header.rdata[1];
	    double n_global[3];
	    n_global[0]            = header.rdata[2];
	    n_global[1]            = header.rdata[3];
	    n_global[2]            = header.rdata[4];
	    double x_global[3];
	    x_global[0]            = header.rdata[5];
	    x_global[1]            = header.rdata[6];
	    x_global[2]            = header.rdata[7];
	    cout << " > bf zone \"" << header.name << "\", nbf: " << nbf_global <<
	      ", area: " << area_global <<
	      ", n: " << COUT_VEC(n_global) <<
	      ", x: " << COUT_VEC(x_global) <<
	      ", mag(n)/area: " << MAG(n_global)/area_global << endl;
	  }
          break;

        case UGP_IO_SURFACE:
          if (option == LESINFO_SBIN) {
                   
            if (b_periodic) {
              offset_surface = offset;
            }
            else { 
              FILE * fp_sbin = fopen(name.c_str(),"wb");
              assert(fp_sbin != NULL);

              const int version = 1;
              fwrite(&version,sizeof(int),1,fp_sbin);
              const int count = zoneVec.size();
              fwrite(&count,sizeof(int),1,fp_sbin);
              cout << " > writing " << count << " zones" << endl;
              for (int izone = 0, limit = zoneVec.size(); izone < limit; ++izone) {
                const int length = zoneVec[izone].length();
                fwrite(&length,sizeof(int),1,fp_sbin);
                fwrite(zoneVec[izone].c_str(),sizeof(char),length,fp_sbin);
                cout << " > zone " << zoneVec[izone] << endl;
              }

              const int nsp = header.idata[0];
              cout << " > nsp = " << nsp << endl;
              fwrite(&nsp,sizeof(int),1,fp_sbin);
              
              // xsp...
              const int read_size = 1000000;
              double (*dbuf3)[3] = new double[min(nsp,read_size)][3];
              int index0 = 0;
              while (index0 < nsp) {
                const int this_read_size = min(read_size,nsp-index0);
                fread(dbuf3,sizeof(double),this_read_size*3,fp);
                if (byte_swap) ByteSwap::byteSwap(dbuf3,this_read_size);
                fwrite(dbuf3,sizeof(double),this_read_size*3,fp_sbin);
                index0 += this_read_size;
              }
              delete[] dbuf3;

              const int nst = header.idata[1];
              cout << " > nst = " << nst << endl;
              fwrite(&nst,sizeof(int),1,fp_sbin);

              // spost...
              int (*ibuf3)[3] = new int[min(nst,read_size)][3];
              index0 = 0;
              while (index0 < nst) {
                const int this_read_size = min(read_size,nst-index0);
                fread(ibuf3,sizeof(int),this_read_size*3,fp);
                if (byte_swap) ByteSwap::byteSwap(ibuf3,this_read_size);
                fwrite(ibuf3,sizeof(int),this_read_size*3,fp_sbin);
                index0 += this_read_size;
              }
              delete[] ibuf3;

              // znost...
              int *ibuf = new int[min(nst,read_size)];
              index0 = 0;
              while (index0 < nst) {
                const int this_read_size = min(read_size,nst-index0);
                fread(ibuf,sizeof(int),this_read_size,fp);
                if (byte_swap) ByteSwap::byteSwap(ibuf,this_read_size);
                fwrite(ibuf,sizeof(int),this_read_size,fp_sbin);
                index0 += this_read_size;
              }
              delete[] ibuf;
              fclose(fp_sbin);

            }
          }
          break;
        case UGP_IO_SURFACE_PERIODIC_INFO:
          if (option == LESINFO_SBIN) {
            assert(b_periodic);
            assert(offset_surface > -1);
                   
            // add periodic zones to zone vec...
            int8 offset_periodic = offset;
            const int npz  = header.idata[1];
            cout << " > npz: " << npz << endl;
            for (int ipz = 0; ipz < npz; ++ipz) {
              stringstream oss; oss << ipz;
              zoneVec.push_back("periodic_zone_" + oss.str());
            }
            
            // go back to UGP_IO_SURFACE...
            my_fseek(fp,offset_surface);
            fread(&header,sizeof(Header),1,fp);
            if (byte_swap)
              byteSwapHeader(&header,1);

            FILE * fp_sbin = fopen(name.c_str(),"wb");
            assert(fp_sbin != NULL);

            const int version = 1;
            fwrite(&version,sizeof(int),1,fp_sbin);
            const int count = zoneVec.size();
            fwrite(&count,sizeof(int),1,fp_sbin);
            cout << " > writing " << count << " zones" << endl;
            for (int izone = 0, limit = zoneVec.size(); izone < limit; ++izone) {
              const int length = zoneVec[izone].length();
              fwrite(&length,sizeof(int),1,fp_sbin);
              fwrite(zoneVec[izone].c_str(),sizeof(char),length,fp_sbin);
              cout << " > zone " << zoneVec[izone] << endl;
            }

            const int nsp = header.idata[0];
            cout << " > nsp = " << nsp << endl;
            fwrite(&nsp,sizeof(int),1,fp_sbin);
            
            // xsp...
            const int read_size = 1000000;
            double (*dbuf3)[3] = new double[min(nsp,read_size)][3];
            int index0 = 0;
            while (index0 < nsp) {
              const int this_read_size = min(read_size,nsp-index0);
              fread(dbuf3,sizeof(double),this_read_size*3,fp);
              if (byte_swap) ByteSwap::byteSwap(dbuf3,this_read_size);
              fwrite(dbuf3,sizeof(double),this_read_size*3,fp_sbin);
              index0 += this_read_size;
            }
            delete[] dbuf3;

            const int nst = header.idata[1];
            cout << " > nst = " << nst << endl;
            fwrite(&nst,sizeof(int),1,fp_sbin);

            // spost...
            int (*ibuf3)[3] = new int[min(nst,read_size)][3];
            index0 = 0;
            while (index0 < nst) {
              const int this_read_size = min(read_size,nst-index0);
              fread(ibuf3,sizeof(int),this_read_size*3,fp);
              if (byte_swap) ByteSwap::byteSwap(ibuf3,this_read_size);
              fwrite(ibuf3,sizeof(int),this_read_size*3,fp_sbin);
              index0 += this_read_size;
            }
            delete[] ibuf3;

            // znost...
            int *ibuf = new int[min(nst,read_size)];
            index0 = 0;
            while (index0 < nst) {
              const int this_read_size = min(read_size,nst-index0);
              fread(ibuf,sizeof(int),this_read_size,fp);
              if (byte_swap) ByteSwap::byteSwap(ibuf,this_read_size);
              fwrite(ibuf,sizeof(int),this_read_size,fp_sbin);
              index0 += this_read_size;
            }
            delete[] ibuf; ibuf = NULL;

            // go back to UGP_IO_SURFACE_PERIODIC_INFO...
            my_fseek(fp,offset_periodic);
            fread(&header,sizeof(Header),1,fp);
            if (byte_swap)
              byteSwapHeader(&header,1);

            const int npt = header.idata[0];
            cout << " > periodic transforms: " << npt << endl;
            fwrite(&npt,sizeof(int),1,fp_sbin);

            // periodic transforms...
            for (int ipt = 0; ipt < npt; ++ipt) {
              const int kind = header.idata[3+ipt];
              fwrite(&kind,sizeof(int),1,fp_sbin);
              double data[3];
              FOR_I3 data[i] = header.rdata[ipt*3+i];
              fwrite(data,sizeof(double),3,fp_sbin);
            }
            
            assert(npz == header.idata[1]);
            cout << " > periodic zones: " << npz << endl;
            fwrite(&npz,sizeof(int),1,fp_sbin);

            // periodic zone...
            assert(ibuf == NULL); ibuf = new int[npz];
            fread(ibuf,sizeof(int),npz,fp);
            if (byte_swap) ByteSwap::byteSwap(ibuf,npz);
            fwrite(ibuf,sizeof(int),npz,fp_sbin);
            delete[] ibuf; ibuf = NULL;

            // periodic zone bits...
            uint2 * sbuf = new uint2[npz];
            fread(sbuf,sizeof(uint2),npz,fp);
            if (byte_swap) ByteSwap::byteSwap(sbuf,npz);
            fwrite(sbuf,sizeof(uint2),npz,fp_sbin);
            delete[] sbuf; sbuf = NULL;
              
            const int npbi = header.idata[2];
            cout << " > npbi: " << npbi << endl;
            fwrite(&npbi,sizeof(int),1,fp_sbin);

            // periodic sp...
            assert(ibuf == NULL); ibuf = new int[min(npbi,read_size)];
            index0 = 0;
            while (index0 < npbi) {
              const int this_read_size = min(read_size,npbi-index0);
              fread(ibuf,sizeof(int),this_read_size,fp);
              if (byte_swap) ByteSwap::byteSwap(ibuf,this_read_size);
              fwrite(ibuf,sizeof(int),this_read_size,fp_sbin);
              index0 += this_read_size;
            }

            // perioic sp parent ...
            index0 = 0;
            while (index0 < npbi) {
              const int this_read_size = min(read_size,npbi-index0);
              fread(ibuf,sizeof(int),this_read_size,fp);
              if (byte_swap) ByteSwap::byteSwap(ibuf,this_read_size);
              fwrite(ibuf,sizeof(int),this_read_size,fp_sbin);
              index0 += this_read_size;
            }
            delete[] ibuf;

            // periodic sp bits...
            assert(sbuf == NULL); sbuf = new uint2[min(npbi,read_size)];
            index0 = 0;
            while (index0 < npbi) {
              const int this_read_size = min(read_size,npbi-index0);
              fread(sbuf,sizeof(uint2),this_read_size,fp);
              if (byte_swap) ByteSwap::byteSwap(sbuf,this_read_size);
              fwrite(sbuf,sizeof(uint2),this_read_size,fp_sbin);
              index0 += this_read_size;
            }
            delete[] sbuf;

            fclose(fp_sbin);
          }
          break;
	case UGP_IO_FA_ZONE:
	  if (verbose) cout << "FA_I1: fa_zone" << endl;
	  break;
	case UGP_IO_CV_ZONE_HEADER:
	  if (verbose) cout << "cv zone \"" << header.name << "\", kind=" << getCvZoneNameForKind(header.idata[0]) << 
                         ", index=" << header.idata[1] << ", ncv=" << header.idata[2] << endl;
	  break;
	case UGP_IO_CV_ZONE:
	  if (verbose) cout << "CV_I1: cv_zone" << endl;
	  break;
	case UGP_IO_X_NO:
	  if (verbose) cout << "nodal coordinates x_no" << endl;
	  break;

	case UGP_IO_I0:
	  if (verbose) cout << "value I0: " << header.name << "=" << header.idata[0] << endl;
	  break;
	case UGP_IO_D0:
	  if (verbose) cout << "value R0: " << header.name << "=" << header.rdata[0] << endl;
	  break;

	case UGP_IO_NO_I1:
	  if (verbose) cout << "scalar NO_I1: " << header.name << endl;
	  break;
	case UGP_IO_NO_D1:
	  if (verbose) cout << "scalar NO_R1: " << header.name << endl;
	  break;
	case UGP_IO_NO_D2:
	  if (verbose) cout << "vector NO_R2: " << header.name << endl;
	  break;

	case UGP_IO_FA_I1:
	  if (verbose) cout << "scalar FA_I1: " << header.name << endl;
	  break;
	case UGP_IO_FA_D1:
	  if (verbose) cout << "scalar FA_R1: " << header.name << endl;
	  break;
	case UGP_IO_FA_D2:
	  if (verbose) cout << "vector FA_R2: " << header.name << endl;
	  break;

	case UGP_IO_BF_D1:
	  if (verbose) cout << "scalar BF_R1: " << header.name << endl;
	  break;
	case UGP_IO_BF_D2:
	  if (verbose) cout << "vector BF_R2: " << header.name << endl;
	  break;

	case UGP_IO_CV_I1:
	  if (verbose) cout << "scalar CV_I1: " << header.name << endl;
	  break;
	case UGP_IO_CV_I1_NEW:
	  if (verbose) cout << "scalar CV_I1_NEW: " << header.name << " (ncv=" << header.ui8data[0] << ")" << endl;
	  break;
	case UGP_IO_CV_F1:
	  if (verbose) cout << "scalar CV_F1: " << header.name << endl;
	  if ((option == LESINFO_EXTRACT)&&(name == header.name)) {
	    cout << "# extract scalar CV_F1: " << name << " from file " << argv[iarg];
	    // if we are in fast mode, then get the type...
	    if (fast_offset == -2) {
	      fast_offset = offset+256;
	      fast_type = UGP_IO_CV_F1;
	    }
	    else {
	      assert(fast_offset == -1);
	    }
	    for (int i = 0; i < global_index_array_size; ++i) {
	      assert(global_index_array[i] < header.idata[0]);
	      my_fseek(fp,offset+256+global_index_array[i]*sizeof(float)); // 256 is header size
	      float var;
	      fread(&var,sizeof(float),1,fp);
	      cout << " " << var;
	    }
	    cout << endl;
	    done = 1;
	  }
	  break;
	case UGP_IO_CV_D1_NEW:
	  if (verbose) cout << "scalar CV_D1_NEW: " << header.name << " (ncv=" << header.ui8data[0] << ")" << endl;
	  if ((option == LESINFO_EXTRACT_ALL)&&(nameSet.find(header.name) != nameSet.end())) {
	    cout << "# extract scalar CV_D1_NEW " << header.name << " " << header.ui8data[0] << endl;
	    const uint8 read_size = 1000000;
	    double (*buf) = new double[min(header.ui8data[0],read_size)];
	    uint8 index0 = 0;
	    while (index0 < header.ui8data[0]) {
	      const uint8 this_read_size = min(read_size,header.ui8data[0]-index0);
	      fread(buf,sizeof(double),this_read_size,fp);
	      for (int i = 0; i < this_read_size; ++i)
		cout << index0+i << " " << buf[i] << endl;
	      index0 += this_read_size;
	    }
	    delete[] buf;
	  }
          break;
	case UGP_IO_CV_D1:
	  if (verbose) cout << "scalar CV_D1: " << header.name << " (ncv=" << header.idata[0] << ")" << endl;
	  if ((option == LESINFO_EXTRACT)&&(name == header.name)) {
	    cout << "# extract scalar CV_D1: " << name << " from file " << argv[iarg];
	    for (int i = 0; i < global_index_array_size; ++i) {
	      assert(global_index_array[i] < header.idata[0]);
	      my_fseek(fp,offset+256+global_index_array[i]*sizeof(double)); // 256 is header size
	      double var;
	      fread(&var,sizeof(double),1,fp);
	      cout << "  " << var;
	    }
	    cout << endl;
	    done = 1;
	  }
	  else if ((option == LESINFO_EXTRACT_ALL)&&(nameSet.find(header.name) != nameSet.end())) {
	    cout << "CV_R1 " << header.name << " " << header.idata[0] << endl;
	    const int read_size = 1000000;
	    double (*buf) = new double[min(header.idata[0],read_size)];
	    int8 index0 = 0;
	    while (index0 < header.idata[0]) {
	      const int this_read_size = min(read_size,int(header.idata[0]-index0));
	      fread(buf,sizeof(double),this_read_size,fp);
	      for (int i = 0; i < this_read_size; ++i)
		cout << buf[i] << endl;
	      index0 += this_read_size;
	    }
	    delete[] buf;
	  }
	  else if ((option == LESINFO_EXTRACT_TO_FILE)&&(nameSet.find(header.name) != nameSet.end())) {
	    if (fp_extract == NULL) {
	      fp_extract = fopen("extract.sles","wb");
	      int itmp[2] = { UGP_IO_MAGIC_NUMBER+1, 5 }; // v5 io
	      fwrite(itmp,sizeof(int),2,fp_extract);
	    }
	    // write in float...
	    Header header_extract = header;
	    header_extract.id = UGP_IO_CV_F1;
	    assert(header_extract.skip == sizeof(Header) + header_extract.idata[0]*8);
	    header_extract.skip = sizeof(Header) + header_extract.idata[0]*4; // 4 bytes per float
	    fwrite(&header_extract,sizeof(Header),1,fp_extract);
	    cout << "CV_R1 " << header.name << " " << header.idata[0] << endl;
	    const int read_size = 1000000;
	    double (*buf) = new double[min(header.idata[0],read_size)];
	    float (*fbuf) = new float[min(header.idata[0],read_size)];
	    int8 index0 = 0;
	    while (index0 < header.idata[0]) {
	      const int this_read_size = min(read_size,int(header.idata[0]-index0));
	      cout << " > extracting " << index0 << " out of " << header.idata[0] << " doubles..." << endl;
	      fread(buf,sizeof(double),this_read_size,fp);
	      for (int i = 0; i < this_read_size; ++i)
		fbuf[i] = float(buf[i]);
	      fwrite(fbuf,sizeof(float),this_read_size,fp_extract);
	      index0 += this_read_size;
	    }
	    delete[] buf;
	    delete[] fbuf;
	    cout << "OK" << endl;
	  }
          else if ( (option == LESINFO_PBIN) &&
		    ((strcmp(header.name,"vol_vv") == 0)||(strcmp(header.name,"VOL_VV") == 0)||(strcmp(header.name,"vol_cv") == 0)) ) {
            //store offset for use when points header is found,
            //lengthscale must be written to pbin after points
            offset_vol = offset;
            cout << "# Found " << header.name << ", storing offset" << endl;
	  }
          break;
	case UGP_IO_CV_D2_NEW:
	  if (verbose) cout << "vector CV_D2_NEW: " << header.name << " (ncv=" << header.ui8data[0] << ")" << endl;
          break;
	case UGP_IO_CV_D2:
	  if (verbose) cout << "vector CV_D2: " << header.name << " (ncv=" << header.idata[0] << ")" << endl;
	  if ((option == LESINFO_MATCH)&&(name == header.name)) {
	    cout << "# matched vector CV_D2: " << name << ". Reading " << header.idata[0] << "x3 doubles";
	    int8 index_min = -1;
	    double d2_min;
	    const int read_size = 1000000;
	    double (*buf)[3] = new double[min(header.idata[0],read_size)][3];
	    int8 index0 = 0;
	    while (index0 < header.idata[0]) {
	      cout << ".";
	      cout.flush();
	      const int this_read_size = min(read_size,int(header.idata[0]-index0));
	      fread(buf,sizeof(double),this_read_size*3,fp);
	      for (int i = 0; i < this_read_size; ++i) {
		const double this_d2 = DIST2(buf[i],xp);
		if ((index_min == -1)||(this_d2 < d2_min)) {
		  index_min = index0 + i;
		  d2_min = this_d2;
		}
	      }
	      index0 += this_read_size;
	    }
	    delete[] buf;
	    cout << "OK\n# closest global index is " << index_min << ", distance " << sqrt(d2_min) << endl;
	    done = 1;
	  }
	  else if ((option == LESINFO_EXTRACT)&&(name == header.name)) {
	    cout << "# extract vector CV_D2: " << name << " from file " << argv[iarg];
	    for (int i = 0; i < global_index_array_size; ++i) {
	      assert(global_index_array[i] < header.idata[0]);
	      my_fseek(fp,offset+256+global_index_array[i]*sizeof(double)*3); // 256 is header size
	      double var[3];
	      fread(var,sizeof(double),3,fp);
	      cout << "  " << var[0] << " " << var[1] << " " << var[2];
	    }
	    cout << endl;
	    done = 1;
	  }
	  else if ((option == LESINFO_EXTRACT_ALL)&&(nameSet.find(header.name) != nameSet.end())) {
	    cout << "CV_R2 " << header.name << " " << header.idata[0] << " 3" << endl;
	    const int read_size = 1000000;
	    double (*buf)[3] = new double[min(header.idata[0],read_size)][3];
	    int8 index0 = 0;
	    while (index0 < header.idata[0]) {
	      const int this_read_size = min(read_size,int(header.idata[0]-index0));
	      fread(buf,sizeof(double),this_read_size*3,fp);
	      for (int i = 0; i < this_read_size; ++i)
		cout << buf[i][0] << " " << buf[i][1] << " " << buf[i][2] << endl;
	      index0 += this_read_size;
	    }
	    delete[] buf;
	  }
	  else if ((option == LESINFO_EXTRACT_TO_FILE)&&(nameSet.find(header.name) != nameSet.end())) {
	    if (fp_extract == NULL) {
	      fp_extract = fopen("extract.sles","wb");
	      int itmp[2] = { UGP_IO_MAGIC_NUMBER+1, UGP_IO_VERSION };
	      fwrite(itmp,sizeof(int),2,fp_extract);
	    }
	    // write in float...
	    Header header_extract = header;
	    header_extract.id = UGP_IO_CV_F2;
	    assert(header_extract.skip == sizeof(Header) + header_extract.idata[0]*8*3);
	    header_extract.skip = sizeof(Header) + header_extract.idata[0]*4*3; // 4 bytes per float
	    fwrite(&header_extract,sizeof(Header),1,fp_extract);
	    cout << "CV_R2 " << header.name << " " << header.idata[0] << endl;
	    const int read_size = 1000000;
	    double (*buf)[3] = new double[min(header.idata[0],read_size)][3];
	    float (*fbuf)[3] = new float[min(header.idata[0],read_size)][3];
	    int8 index0 = 0;
	    while (index0 < header.idata[0]) {
	      const int this_read_size = min(read_size,int(header.idata[0]-index0));
	      cout << " > extracting " << index0 << " out of " << header.idata[0] << " double[3]s..." << endl;
	      fread(buf,sizeof(double),this_read_size*3,fp);
	      for (int i = 0; i < this_read_size; ++i)
		FOR_J3 fbuf[i][j] = float(buf[i][j]);
	      fwrite(fbuf,sizeof(float),this_read_size*3,fp_extract);
	      index0 += this_read_size;
	    }
	    delete[] buf;
	    delete[] fbuf;
	    cout << "OK" << endl;
	  }
	  else if ( (option == LESINFO_PBIN) &&
		    ((strcmp(header.name,"x_vv") == 0)||(strcmp(header.name,"X_VV") == 0)||(strcmp(header.name,"X_CV") == 0)) ) {
	    FILE * fp_pbin = fopen(name.c_str(),"wb");
	    //const int POINTS_IO_MAGIC_NUMBER = 1235813;
	    //const int POINTS_IO_VERSION = 2;
	    const int ncv_pbin = header.idata[0];
            cout << "# Reading points for " << name << endl;
	    int8 ibuf[4] = { POINTS_IO_MAGIC_NUMBER, POINTS_IO_VERSION, header.idata[0], 1 }; // last is include delta or not maybe?
	    fwrite(ibuf,sizeof(int8),4,fp_pbin);
	    const int read_size = 1000000;
	    double (*buf)[3] = new double[min(header.idata[0],read_size)][3];
	    double bbmin[3] = { 1.0E+20, 1.0E+20, 1.0E+20 };
	    double bbmax[3] = { -1.0E+20, -1.0E+20, -1.0E+20 };
	    int8 index0 = 0;
	    while (index0 < header.idata[0]) {
	      const int this_read_size = min(read_size,int(header.idata[0]-index0));
	      fread(buf,sizeof(double),this_read_size*3,fp);
	      for (int i = 0; i < this_read_size; ++i) {
		FOR_J3 bbmin[j] = min(bbmin[j],buf[i][j]);
		FOR_J3 bbmax[j] = max(bbmax[j],buf[i][j]);
	      }
	      fwrite(buf,sizeof(double),this_read_size*3,fp_pbin);
	      index0 += this_read_size;
	    }
            delete[] buf;

	    double *rvv = new double[min(header.idata[0],read_size)];
            int b_eof = 0;
            //Try to read volume record to compute length scale...
            if (offset_vol>-1) { //record was already found, skip to it and read
              b_eof = 1;
              my_fseek(fp,offset_vol);
              fread(&header,sizeof(Header),1,fp);
              if (byte_swap)
                byteSwapHeader(&header,1);
            }
            else { //record may still be later in the file, continue searching..
              offset += header.skip;
              while (b_eof==0){
                my_fseek(fp,offset);
                fread(&header,sizeof(Header),1,fp);
                if (byte_swap)
                  byteSwapHeader(&header,1);

                if (header.skip <= 0) {
                  b_eof = -1;
                }
                switch (header.id) {
                case UGP_IO_CV_D1:
                  if ((strcmp(header.name,"vol_vv") == 0)||(strcmp(header.name,"VOL_VV") == 0)||(strcmp(header.name,"vol_cv") == 0)){
                    b_eof = 1;
                  }
                  break;
                case UGP_IO_EOF:
                  b_eof = -1;
                  break;
                default:
                  break;
                  //keep searching
                }
                offset += header.skip; 
              }
            }

            double power = 1.0/3.0;
            if (b_eof == 1) { //we found and are ready to read a cv volume record
              cout << "# Estimating delta from cell volume " << endl;
              index0 = 0;
	      while (index0 < header.idata[0]) {
	        const int this_read_size = min(read_size,int(header.idata[0]-index0));
	        fread(rvv,sizeof(double),this_read_size,fp);
	        for (int i = 0; i < this_read_size; ++i) {
                  rvv[i] = pow(rvv[i],power);
	        }
	        fwrite(rvv,sizeof(double),this_read_size,fp_pbin);
	        index0 += this_read_size;
	      }
            }
            else { // no volume record found, estimate a lengthscale
	      // and set the delta to something small...
              // try max dimension of bbox divided by cube root of total point count
	      const double delta = 0.1*max(bbmax[0]-bbmin[0],max(bbmax[1]-bbmin[1],bbmax[2]-bbmin[2])) / pow(ncv_pbin,power);
              cout << "# Estimating delta from domain length scale and point count" << endl;
              cout << "#   delta = " << delta << endl;

	      for (int i = 0; i < min(ncv_pbin,read_size); ++i)
	        rvv[i] = delta;

	      index0 = 0;
	      while (index0 < ncv_pbin) {
	        const int this_read_size = min(read_size,int(ncv_pbin-index0));
	        fwrite(rvv,sizeof(double),this_read_size,fp_pbin);
	        index0 += this_read_size;
	      }
            }
	    fclose(fp_pbin);

            cout << "# " << name << " writing complete."  << endl;
	    delete[] rvv;
            done = 1;
	  }
	  break;
	case UGP_IO_CV_F2:
	  if (verbose) cout << "vector CV_F2: " << header.name << endl;
          break;
        case UGP_IO_FAZONE_NO_I1:
          if ( verbose) cout << "fazone scalar NO_I1: " << header.name << endl ;
          break;
        case UGP_IO_FAZONE_FA_I1:
          if ( verbose) cout << "fazone scalar FA_I1: " << header.name << endl ;
          break;
        case UGP_IO_FAZONE_NO_D1:
          if ( verbose) cout << "fazone scalar NO_R1: " << header.name << endl ;
          break ;
        case UGP_IO_FAZONE_FA_D1:
          if ( verbose) cout << "fazone scalar FA_R1: " << header.name << endl ;
          break;
        case UGP_IO_FAZONE_NO_D2:
          if ( verbose) cout << "fazone vector NO_R2: " << header.name << endl ;
          break;
        case UGP_IO_FAZONE_FA_D2:
          if ( verbose) cout << "fazone vector FA_R2: " << header.name << endl ;
          break;

	case UGP_IO_EOF:
	  if (verbose) cout << "done." << endl;
	  done = 1;
	  break;

	  // chemtable stuff
	case UGP_IO_CT_KD:
	  if (verbose) cout << "Chemical table: " << header.name << endl;
	  if (verbose) cout << " > Coordinate dimension:    " <<  header.idata[0] << endl;
	  if (verbose) cout << " > Data dimension:          " <<  header.idata[1] << endl;
	  if (verbose) cout << " > Number of Verticies:     " <<  header.idata[2] << endl;
	  if (verbose) cout << " > Number of Nodes in Tree: " <<  header.idata[4] << endl;
	  if (verbose) cout << " > Number of Leaves:        " <<  header.idata[3] << endl;
	  if (verbose) cout << " > Input-specified rtol:    " <<  header.rdata[0] << endl;
	  if (strcmp(header.name, "VIDA_NON-PREMIXED_FPV_KDT3D") == 0 ||
	      strcmp(header.name, "CHRIS_NON-PREMIXED_FPV_KDT3D") == 0) {
	    if (verbose) cout << " > Reference Pressure:      " <<  header.rdata[1] << endl;
	    if (verbose) cout << " > Oxidizer Temperature:    " <<  header.rdata[2] << endl;
	    if (verbose) cout << " > Fuel Temperature:        " <<  header.rdata[3] << endl;
	  }
	  break;
	case UGP_IO_CT_CART_1D:
	  if (verbose) cout << "Chemical table: " << header.name << endl;
	  if (verbose) cout << " > Dimensions(n1,nvars): (" <<
                         header.idata[0] << "," << header.idata[1] << ")" << endl;
	  break;
	case UGP_IO_CT_CART_2D:
	  if (verbose) cout << "Chemical table: " << header.name << endl;
	  if (verbose) cout << " > Dimensions(n1,n2,nvars): (" <<
                         header.idata[0] << "," << header.idata[1] <<
                         "," << header.idata[2] << ")" << endl;
	  break;
	case UGP_IO_CT_CART_3D:
	  if (verbose) cout << "Chemical table: " << header.name << endl;
	  if (verbose) cout << " > Dimensions(n1,n2,n3,nvars): (" <<
                         header.idata[0] << "," << header.idata[1] <<
                         "," << header.idata[2] << "," << header.idata[3] << ")" << endl;
	  break;
	case UGP_IO_CT_CART_4D:
	  if (verbose) cout << "Chemical table: " << header.name << endl;
	  if (verbose) cout << " > Dimensions(n1,n2,n3,n4,nvars): (" <<
                         header.idata[0] << "," << header.idata[1] <<
                         "," << header.idata[2] << "," << header.idata[3] <<
                         "," << header.idata[4] << ")" << endl;
	  break;
	case UGP_IO_CT_COOR:
	  if (verbose) cout << " > Chemical table coordinate:  " << header.name << endl;
	  break;
	case UGP_IO_CT_DATA:
	  if (verbose) cout << " > Chemical table data:        " << header.name << endl;
	  break;
	case UGP_IO_CT_KD_HC:
	  // nothing to report here...
	  break;
	case UGP_IO_CT_KD_NODE:
	  // ditto
	  break;
	default:
	  if (verbose) cout << " > got header: " << header.id << " " << header.name << endl;
	  break;
	}

	//getchar();
	offset += header.skip;

      }

      fclose(fp);
      
      if (fp_extract != NULL) {
	Header header_extract;
	header_extract.id = UGP_IO_EOF;
	header_extract.skip = sizeof(Header);
	fwrite(&header_extract,sizeof(Header),1,fp_extract);
	fclose(fp_extract);
      }
      
    }

  }
  catch (int e) {
    return(-1);
  }
  catch(...) {
    return(-1);
  }
  return(0);
}
