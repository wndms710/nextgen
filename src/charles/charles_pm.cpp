
#include "FlowSolver.hpp"
#include "PremixedSolver.hpp"

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

  void initialHook() {
    // if you change the state, set b_reset_state to true to recompute the state
    // based on the current scalars and chemtable pressure... 
    /*
    COUT1("WARNING: resetting C for x_vv[0] <= 0.04");
    FOR_ICV {
      if (x_vv[icv][0]<=0.04) {
	C[icv] = 0.0;
      }
    }
    b_reset_state = true;
    */
  }

  void temporalHook() {}

  void finalHook() {}

  PremixedBc* initHookBc(BfZone* p) {
    CERR(" > user must supply a hook boundary condition for zone: " << p->getName());
  }

  void addSourceHook(PremixedRhs* rhs, const double time, const int rk_stage) {}

};

int main(int argc, char* argv[]) {

  try {

    CTI_Init(argc,argv,"charles.in");

    {

      const bool b_post = checkParam("POST");

      if (Param * param = getParam("EOS")) {
        const string eos = param->getString();
        if ((eos == "PREMIXED_FPV") || (eos == "PREMIXED")) {

          MyPremixedSolver solver;

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
               "EOS PREMIXED_FPV");
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
