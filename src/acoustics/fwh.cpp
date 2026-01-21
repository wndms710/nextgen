#include "fwh.hpp"
// ======================================================
// main
// ======================================================

int main(int argc, char * argv[]) {

  try {
       
    // initialize the environment: basically initializes MPI and parameters...
    
    CTI_Init(argc,argv,"fwh.in");
    
    // the solver...
    
    {
      
      FWH solver; 
      
      solver.run();
      
      solver.finalize();
     
    }

    
    // finalize the environment: reports parameter usage, shuts down MPI...
    
    CTI_Finalize();
    
    // rudimentary error management...

  }
  catch (int e) {
    if (e == 0)
      CTI_Finalize();
    else
      CTI_Abort();
  }
  catch(...) {
    CTI_Abort();
  }
    
  return(0);
  
}

