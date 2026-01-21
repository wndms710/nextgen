
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

  /*
  double (* vort)[3];         // = curl(u);
  double * vort_mag;          // = |curl(u)|
  */

  MyIdealGasSolver() {

    /*
    vort     = NULL;  registerCvData( vort, "vort", READWRITE_DATA);
    vort_mag = NULL;  registerCvData( vort_mag, "vort_mag", NO_READWRITE_DATA);
    */
  }

  void initData() {

    IdealGasSolver::initData();

    /*
    assert( vort == NULL);     vort     = new double[ncv_g2][3];  // include vort calculation in ghosts...
    assert( vort_mag == NULL); vort_mag = new double[ncv];        // no ghosts ...
    */

  }

  ~MyIdealGasSolver() {

    /*
    DELETE(vort);
    DELETE(vort_mag);
    */

  }

  void initialHook() {}

  void temporalHook() {

    /*
    if ( step % check_interval == 0) {

      if ( mpi_rank == 0 )
        cout << " > computing vorticity ... " << endl;

      for (int icv = 0; icv < ncv; ++icv) {

        vort[icv][0]  = dudx[icv][2][1] - dudx[icv][1][2];
        vort[icv][1]  = dudx[icv][0][2] - dudx[icv][2][0];
        vort[icv][2]  = dudx[icv][1][0] - dudx[icv][0][1];

        vort_mag[icv] = MAG(vort[icv]);

      }

      updateCv2Data(vort, REPLACE_ROTATE_DATA);

    }
    */

  }

  void finalHook() {}

  IdealGasBc* initHookBc(BfZone* p) {
    CERR(" > user must supply a hook boundary condition for zone: " << p->getName());
  }

  void addSourceHook(IdealGasRhs* rhs, const double time, const int rk_stage) {}

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

  void temporalHook() {

    // if you like to check the containment of particles in the domain
    if (step%check_interval==0) {
      checkLpContain();
    }

  }

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
    if (e >= 0) {
      CTI_Finalize();
    }
    else {
      CTI_Abort();
    }
  }
  catch(...) {
    CTI_Abort();
  }

  return 0;

}
