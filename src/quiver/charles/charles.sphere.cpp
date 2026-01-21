
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
// IdealGasSolver
//===============================

class MyIdealGasSolver : public IdealGasSolver { 
public: 

  vector<SimpleSphere> bars;

  MyIdealGasSolver() {
  
    // 
    // parse the bar information ... 
    //


    FOR_PARAM_MATCHING("ADD_SPHERE") { 

      bool b_xp    = false;
      double rp    = -1.0;
      double xp[3] = {0.0,0.0,0.0};

      int iarg     = 0;
      while ( iarg < param->size()) { 

        string token = param->getString(iarg++);
        
        if ( token == "POINT") { 
         
          xp[0] = param->getDouble(iarg++);
          xp[1] = param->getDouble(iarg++);
          xp[2] = param->getDouble(iarg++);
          b_xp  = true;
        
        } else if ( token == "RADIUS") { 

          rp = param->getDouble(iarg++);

        } else { 

          CERR( " > unrecognized add_sphere token: " << token);

        }
      }

      if ( (rp < 0.0) || !b_xp) { 
        CERR( " > bad syntax: ADD_SPHERE POINT <x0> <x1> <x2> RADIUS <r0>");
      }

      bars.push_back(SimpleSphere(xp,rp));


    }
  
  } 

  void initData() { 
    
    IdealGasSolver::initData();

  }

  ~MyIdealGasSolver() {}

  void initialHook() {} 
    
  void temporalHook() {}

  void finalHook() {}

  IdealGasBc* initHookBc(BfZone* p) { 
    CERR(" > user must supply a hook boundary condition for zone: " << p->getName());
  }

  void addBarBodyForce(IdealGasRhs* rhs, vector<SimpleSphere>::iterator& sphere) {

    //
    // perhaps consider transposing the loop arrays for efficiency
    // if this isnt an mrf calculation, we should write a slightly 
    // different loop below .. 
    //

    assert( frame_rotation); 

    const double tau = 20.0*dt;

    for (int icv = 0; icv < ncv; ++icv) { 

      // 
      // need to retrieve the velocity at this x_cv in stationary frame
      // 

      if ( sphere->pointIsInside(x_cv[icv])) { 

        double u_stat[3];

        u_stat[0] = -frame_rotation[1]*x_cv[icv][2] + 
                     frame_rotation[2]*x_cv[icv][1];

        u_stat[1] = -frame_rotation[2]*x_cv[icv][0] + 
                     frame_rotation[0]*x_cv[icv][2];

        u_stat[2] = -frame_rotation[0]*x_cv[icv][1] + 
                     frame_rotation[1]*x_cv[icv][0];

        for (int i =0; i < 3; ++i) { 
          const double fi = vol_cv[icv]*rho[icv]*(u_stat[i] - u[icv][i])/tau;
          rhs[icv].rhou[i]    += fi;
          rhs[icv].rhoE       += fi*u[icv][i];
        }
      }

    }

  }

  void addSourceHook(IdealGasRhs* rhs, const double time, const int rk_stage) {
  
    static bool first = true;
    if ( first && (mpi_rank == 0)) { 
      cout << " > adding source hook for the fake bars " << endl;
    }

    for (vector<SimpleSphere>::iterator it = bars.begin(); it != bars.end(); ++it)  
      addBarBodyForce(rhs,it);

    first = false;
  }

};

//===============================
// PremixedSolver
//===============================

class MyPremixedSolver : public PremixedSolver {
public:

  MyPremixedSolver() {}

  void initData() { 

    PremixedSolver::initData();

  }

  ~MyPremixedSolver() {}

  void initialHook() {}

  void temporalHook() {}

  void finalHook() {}

  PremixedBc* initHookBc(BfZone* p) {
    CERR(" > user must supply a hook boundary condition for zone: " << p->getName());
  }

  void addSourceHook(PremixedRhs* rhs, const double time, const int rk_stage) {}

};

//===============================
// HelmholtzSolver
//===============================

class MyHelmholtzSolver : public HelmholtzSolver { 
public:

  MyHelmholtzSolver() {}

  void initData() { 

    HelmholtzSolver::initData();

  }

  ~MyHelmholtzSolver() {}

  void initialHook() {}

  void temporalHook() {}

  void finalHook() {}

  HelmholtzBc* initHookBc(BfZone* p) {
    CERR(" > user must supply a hook boundary condition for zone: " << p->getName());
  }

  
  // the Helmholtz solver has implicit time advancement in a fractional 
  // step setting; as a result, the hooks for add source hooks are slightly
  // different.

  void momentumSourceHook(double * A,double (*rhs)[3]) {}
  void massSourceHook(double * rhs) {}

};


class MyNonpremixedSolver : public NonpremixedSolver {};

class MyIdealGasSolverWithLsp : public IdealGasSolverWithLsp { 
public: 

  MyIdealGasSolverWithLsp() {} 

  void initData() { 
    
    IdealGasSolverWithLsp::initData();

  }

  ~MyIdealGasSolverWithLsp() {}

  void initialHook() {

    IdealGasSolverWithLsp::initialHook();

  } 
    
  void temporalHook() {}

  void finalHook() {

    IdealGasSolverWithLsp::finalHook();

  }

  //IdealGasBc* initHookBc(BfZone* p) { 
  //  CERR(" > user must supply a hook boundary condition for zone: " << p->getName());
  //}

  void addSourceHook(IdealGasRhs* rhs, const double time, const int rk_stage) {}

};

//===============================
// postpro 
//===============================

class MyBasicPostpro : public BasicPostpro { 
public: 

  MyBasicPostpro() {}

  void initData() { 

    BasicPostpro::initData();

  }

  ~MyBasicPostpro() {}

  void initialHook() {}

  void temporalHook() {}
  
  void finalHook() {}

};


int main(int argc, char* argv[]) {
  
  try {
   
    CTI_Init(argc,argv,"charles.in");

    { 

      const bool b_post = checkParam("POST");

      if (Param * param = getParam("EOS")) { 
        const string eos = param->getString();
        if ( eos == "IDEAL_GAS") {
	
          MyIdealGasSolver solver;

          if (b_post) {

            solver.initMin();
            solver.runPost();

          } else { 

            solver.init();
            solver.run();
          } 
        } 
        else if (eos == "IDEAL_GAS_LSP") {
          
          MyIdealGasSolverWithLsp solver;

          if (b_post) {

            solver.initMin();
            solver.runPost();
  
          } else {
  
            solver.init();
            solver.run();
          }        
        }
        else if ((eos == "PREMIXED_FPV") || (eos == "PREMIXED")) { 
          
          MyPremixedSolver solver;

          if ( b_post) { 

            solver.initMin();
            solver.runPost();

          } else { 

            solver.init();
            solver.run();
          }
        }
        else if ((eos == "NONPREMIXED_FPV") || (eos == "NONPREMIXED") || (eos == "NON-PREMIXED") || (eos == "NON-PREMIXED_FPV")) { 

          MyNonpremixedSolver solver;

          if ( b_post) { 

            solver.initMin();
            solver.runPost();

          } else { 

            solver.init();
            solver.run();

          }
        }
 	else if ( eos == "HELMHOLTZ") { 
          
          MyHelmholtzSolver solver;
          
          if ( b_post) { 

            solver.initMin();
            solver.runPost();

          } else { 

            solver.init();
            solver.run();
          }
        } 
	else { 
          CERR("unrecognized EOS: " << eos << ", possible choices are \n" <<
               "EOS IDEAL_GAS\n" << 
               "EOS PREMIXED_FPV\n" << 
               "EOS NONPREMIXED_FPV\n" << 
               "EOS IDEAL_GAS_LSP");
	}
      }
      else {
        
        if ( b_post) { 

          // if there is no eos object provided but we are in a post
          // mode then we need to just start a stripped down staticsolver

          MyBasicPostpro solver;
          solver.initMin();
          solver.runPost();

        } else { 

          CERR("must specify EOS or POST. Possible choices for EOS are \n" <<
               "EOS IDEAL_GAS\n" << 
               "EOS PREMIXED_FPV\n" << 
               "EOS NONPREMIXED_FPV\n" << 
               "EOS IDEAL_GAS_LSP");
        }
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
