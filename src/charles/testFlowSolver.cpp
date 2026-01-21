
#include "FlowSolver.hpp"

// FlowSolver is an abstract class that requires implementations of the 
// following virtual functions...

class MyFlowSolver : public FlowSolver {
public:

  virtual void advanceSolution() { 
    if (mpi_rank == 0) cout << "MyFlowSolver::advanceSolution()" << endl;
  }
  
  virtual void syncPostState() {
    if (mpi_rank == 0) cout << "MyFlowSolver::syncPostState()" << endl;
  }
  
  virtual void initComplete() {
    if (mpi_rank == 0) cout << "MyFlowSolver::initComplete()" << endl;
  }
  
  virtual double calcCfl(double* cfl,const double dt_) const {
    if (mpi_rank == 0) cout << "MyFlowSolver::calcCfl()" << endl;
  }
  
  virtual void computeForces() {
    if (mpi_rank == 0) cout << "MyFlowSolver::computeForces()" << endl;
  }
  
};

int main(int argc, char* argv[]) {
  
  try {
   
    CTI_Init(argc,argv,"testFlowSolver.in");

    { 

      MyFlowSolver solver;
      solver.init();
      solver.run();
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
  
  return 0;

} 
