#include "CTI.hpp"
using namespace CTI;


int main(int argc, char * argv[]) {

  try {
    
    CTI_Init(argc,argv,"testCti.in");
    
    {

      double p1 = getDoubleParam("P1",1.23456);
      int i1 = getIntParam("I1",12345);
      int i2 = getIntParam("I2",12345543);
      int8 i8 = getInt8Param("I8",5732875398);
      string s1 = getStringParam("S1","this-is-string-1");
      bool b = getBoolParam("b",true);
      
      if (mpi_rank == 0) cout << "got params: " << p1 << " " << i1 << " " << i8 << " " << s1 << " " << b << endl;
      
    }
    
    CTI_Finalize();
    
  }
  catch (int e) {
    if (e >= 0)
      CTI_Finalize();
    else
      CTI_Abort();
  }
  catch(...) {
    CTI_Abort();
  }

  return(0);

}

