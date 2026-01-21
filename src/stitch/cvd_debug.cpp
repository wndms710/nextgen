#include "CTI.hpp"
using namespace CTI;
#include "CuttableVoronoiData.hpp"

int main(int argc, char * argv[]) {

  try {
    
    CTI_Init(argc,argv,"cvd_debug.in");

    {

      int debug_count = 0;
      int debug_failed_count = 0;
      bool b_cvd_debug = false;
      FOR_PARAM_MATCHING("CVD_DEBUG") {
	
	b_cvd_debug = true;
	++debug_count;
	
	cout << "\n================================================================================\n" << 
	  "WORKING ON CVD_DEBUG: " << param->getString() << endl;
	
	assert(mpi_size == 1);
	
	CuttableVoronoiData cvd;
	cvd.readBinary(param->getString());
	
	cvd.writeTecplot(0);
	cout << "take a look at tecplot 0: This is the surface before cube cutting." << endl;
	
	assert(param->size() > 1);
	FILE * fp = fopen(param->getString(1).c_str(),"rb");
	double this_delta_cv;
	fread(&this_delta_cv,sizeof(double),1,fp);
	cout << "this_delta_cv: " << this_delta_cv << endl;
	int recv_buf_int_size;
	fread(&recv_buf_int_size,sizeof(int),1,fp);
	cout << "recv_buf_int_size: " << recv_buf_int_size << endl;
	int * recv_buf_int = new int[recv_buf_int_size];
	fread(recv_buf_int,sizeof(int),recv_buf_int_size,fp);
	
	try {
	  
	  cvd.buildSeed(this_delta_cv,recv_buf_int,true);
	  cout << "SUCCESS. cvd.nno: " << cvd.nno << " cvd.ned: " << cvd.ned << endl;
	  
	}
	catch(...) {
	  
	  cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!! FAILED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	  ++debug_failed_count;
	  
	}
	
	delete[] recv_buf_int;
	
      }
      if (b_cvd_debug) {
	cout << "\n================================================================================" << endl;
	cout << "CVD_DEBUG Summary: tested: " << debug_count << 
	  ", passed: " << debug_count-debug_failed_count << 
	  ", failed: " << debug_failed_count << endl;
	cout << "================================================================================\n" << endl;
      }
      
      
    }
    
    CTI_Finalize();

  }
  catch (int e) {
    // integer exceptions are thrown from some parts of the code.
    // negative exceptions are thrown when only one rank is having a problem,
    // so these require a CTI_Abort to free the resources. Positive exceptions
    // e.g. throw(0), throw(1) are associated with points where everyone is
    // synchronized, and we can call CTI_Finalize and shut down cleanly.
    if (e >= 0) { 
      CTI_Finalize();
    }
    else {
      cout << "catch " << e << " on rank: " << mpi_rank << ": calling abort..." << endl;
      CTI_Abort();
      CTI_Finalize();
    }
  }
  catch(...) {
    // don't know what this could be...
    cout << "catch unknown error on rank: " << mpi_rank << ": calling abort..." << endl;
    CTI_Abort();
    CTI_Finalize();
  }

  return(0);

}
