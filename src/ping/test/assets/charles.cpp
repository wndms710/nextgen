
#include "FlowSolver.hpp"
#include "IdealGasSolver.hpp"
#include "PremixedSolver.hpp"
#include "IdealGasSolverWithLsp.hpp"
#include "HelmholtzSolver.hpp"
#include "NonpremixedSolver.hpp"
#include "BasicPostpro.hpp"

//==================================================================================
// following solver types are defined below
//      -- IdealGasSolver : non-reacting, compressible solver
//                          (IdealGasSolverWithLsp -- includes lagrangian particles)
//      -- PremixedSolver : reacting, (partially) premixed compressible solver
//                            based on flamelet-prog var approach
//      -- HelmholtzSolver : non-reacting, low-Ma fractional step solver
//      -- BasicPostpro    : eos agnostic postprocessing code
//
// each of these solvers is has a customized skeleton below that inherits from the 
// base flow solver , e.g., 
// 
//    class MyIdealGasSolver : public IdealGasSolver {};
// 
// that exposes the hooks available to customize the flow solver.  Each solver type
// defined above has a similar list of hooks: 
// 
//         initialHook() : for defining initial conditions or setting data before 
//                         the solver runs
// 
//         temporalHook() : access to the data for diagnostics (or manipulation -- 
//                          but be careful) after the conclusion of each time step 
//                          or snapshot processing
//
//         finalHook()    : access to the data for diagnostics (or manipulation) 
//                          just prior to the solver exit
// 
// Flow solvers also have the ability to provide custom boundary conditions defined 
// through the definition of 
// 
//         initHookBc(BfZone* p) { return new BcObject(); }
//
// where an accompanying boundary condition class must also be specified
// 
//         class BcObject : public ... {};
//
// 
// Flow solvers also have hooks for custom forcing that is added to the rhs of the 
// transport equations (for both the transported vars associated with the flow state
// as well as any passive scalars that were requested with the solution)
//
//         addSourceHook(Rhs* rhs, const double time, const int rk_stage) {} 
// 
// Custom data types can be registered in the constructors including setting of 
// requests for the I/O state (data to be read and/or written, or neither) as annotated
// in the class below.  Memory management associated with registered data is left to 
// the user to allocate in initData(), and subsequently de-allocated.
// 
// The IdealGasSolver is annotated with some possible examples below, but the patterns 
// will apply to all of the flow solvers.  Additional examples of custom setups of 
// cases (boundary conditions, hooks, etc), can be found in the src/quiver directory
// 
//==================================================================================

//===============================
// postpro 
//===============================

class MyBasicPostpro : public BasicPostpro { 

  private:

    double * phi;
  
  public: 

  MyBasicPostpro() {
    phi  = NULL; registerCvData(phi , "phi" , CAN_WRITE_DATA);
  }

  void initData() { 

    BasicPostpro::initData();

    assert(phi  == NULL); phi = new double[ncv];
     
  }

  ~MyBasicPostpro() {}

  void initialHook() {
    FOR_ICV {
      const int image_i = int(x_cv[icv][0]);
      const int image_j = 15-int(x_cv[icv][1]);
      phi[icv] = image_j*16 + image_i;
      COUT1("setting phi[" << icv << "] = " << phi[icv] << " ij " << image_i << " " << image_j << " xy " << x_cv[icv][0] << " " << x_cv[icv][1]);
    }
    dumpRange(phi,ncv,"phi");
  }

  void temporalHook() {}
  
  void finalHook() {}

};


int main(int argc, char* argv[]) {
  
  try {
   
    CTI_Init(argc,argv,"charles.in");
    { 
      const bool b_post = checkParam("POST");
      if ( b_post) { 

         // if there is no eos object provided but we are in a post
         // mode then we need to just start a stripped down staticsolver

         MyBasicPostpro solver;
         solver.initMin();
         solver.runPost();
      }
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
