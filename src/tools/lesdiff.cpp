
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


int main( const int argc, const char* argv[] ) { 

  using namespace ByteSwap;

  try {

    if (sizeof(int8) != 8) {
      cerr << "\n\nError: lesinfo is not compiled with the correct int8 size: " << sizeof(int8) << endl;
      cerr << "try compiling with -D INT8_IS_LONG_INT added to your CXXFLAGS and CXXFLAGS_NO_MPI in Makefile.in\n\n" << endl;
      throw(-1);
    }
    
    assert( sizeof(Header) == 256 );

    if ( argc < 4) { 
      cerr << "\n\n lesdiff usage: \n\n $ lesdiff <file1> <file2> <record-name>" << endl;
      cerr << "\n\n" << endl;
      throw(-1);
    }

    FILE * fp_arr[2]; 
    Header header_arr[2];
    
    const string var_name(argv[3]);

    for (int ifl = 0; ifl < 2; ++ifl) { 

      fp_arr[ifl] = fopen(argv[ifl+1], "rb");

      if (!fp_arr[ifl] ) { 
        cerr << " > unable to open file : " << argv[ifl+1] << endl;
        throw(-1);
      }

      // initialize the offset -- we're going to skip the magic numbers
      // XXX should check for byte swap
      int8 offset = sizeof(int)*2;

      // move down until we find the record we're interested in ... 

      int done = 0;
      while ( done == 0) { 

        my_fseek(fp_arr[ifl],offset);

        fread(&header_arr[ifl],sizeof(Header),1,fp_arr[ifl]);
        //if ( byte_swap) 
        //  byteSwapHeader(&header,1);


        if ( strcmp(header_arr[ifl].name,var_name.c_str()) == 0 ) { 
          
          // this is the record of interest ; advance to the data record

          offset += sizeof(Header);
          my_fseek(fp_arr[ifl],offset);
          done = 1;

        } else if ( header_arr[ifl].id == UGP_IO_EOF) { 

          done = 2;

        }

        offset += header_arr[ifl].skip;
      }

      // ensure that we found the record that we are interested in.. 
      
      if ( done == 2) { 
        cerr << " > unable to find var " << var_name << " in file : " << argv[ifl+1] << endl;
        throw(-1);
      }
    } // ifl

    // check to ensure that the data is consistent 

    if ( header_arr[0].id != header_arr[1].id ) { 
     
      cerr << " > data is not of consistent type " << endl;

    }

    if ( header_arr[0].idata[0] != header_arr[1].idata[0] ) { 

      cerr << " > data is not of consistent length " << endl; 

    }

    const int8 n     = header_arr[0].idata[0];
    const int nchunk = 1000000;
    double * buf0    = new double[nchunk];
    double * buf1    = new double[nchunk];

    cout << " > checking n records = " << n << endl;

    int8 i          = 0;
    double max_diff = 0.0;
    double norm     = 0.0;
    
    while ( i < n ) { 

      const int nn = min(int(n-i),nchunk);

      // the file pointer was already advanced to the start of the data record.. 

      fread(buf0,sizeof(double),nn,fp_arr[0]);
      fread(buf1,sizeof(double),nn,fp_arr[1]);

      for (int ii = 0; ii < nn; ++ii) {
        const double diff = abs(buf0[ii]-buf1[ii]);
        //if ( diff > 1.0e-12) 
        //  cout <<  " XXXX " << int8(ii)+i << "   " << buf0[ii] << "    " << buf1[ii] << "   " << diff << endl;
        max_diff = max(max_diff,abs(buf0[ii]-buf1[ii]));
        norm    += buf0[ii]*buf0[ii];
      }

      if ( i == 0 ) 
        for (int ii = 0 ; ii < 20; ++ii)
          cout << int8(ii) +i << "    " << buf0[ii] << "    " << buf1[ii] << endl;

      i += int8(nn);

    }

    norm = sqrt(norm/double(n));

    cout << " > max diff in " << var_name << " [abs, rel] = " << max_diff << "   " << max_diff/norm << endl;

    for (int ifl = 0; ifl < 2; ++ifl) 
      fclose(fp_arr[ifl]);
    
    delete[] buf0;
    delete[] buf1;

  } catch(...) { 

    return 1;

  }

  return 0;
}

