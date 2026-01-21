#include "CTI.hpp"
using namespace CTI;

#include "StripedMesh.hpp"

int main(int argc, char * argv[]) {

  try {
      
    CTI_Init(argc,argv,"smtest.in");
    
    {
      
      Param * param = getParam("RESTART");
      if (!param) {
	CERR("missing RESTART");
      }
      assert(param->size() >= 1);

      StripedMesh sm;

      // regular read...
      /*
	sm.read_mles(param->getString());
	sm.check();
      */
      
      int ierr = 0;
      if (mpi_rank_shared == 0) {
	ierr = sm.read_mles(param->getString(),mpi_comm_internode);
      }
      MPI_Pause("was that any faster?");

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

