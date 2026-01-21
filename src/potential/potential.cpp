#include "CTI.hpp"
using namespace CTI;

#include "PotentialSolver.hpp"

int main(int argc, char * argv[]) {

  try {
    
    CTI_Init(argc,argv,"potential.in");
    
    {
      
      PotentialSolver solver;
      solver.init();
      solver.run();
      
    }
    
    CTI_Finalize();
    
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

