#ifndef _POSTPRO_SOLVER_HPP_
#define _POSTPRO_SOLVER_HPP_

#include "StaticSolver.hpp"

class PostproSolver : public StaticSolver {
  
private:

public:
  
  PostproSolver() {
    
    COUT1("PostproSolver()");

  }
  
  virtual ~PostproSolver() {
    
    COUT1("~PostproSolver()");

  }
  
  void init() {
    
    COUT1("PostproSolver::init()");
    
    //StaticSolver::init(StaticSolver::INIT_COMPACT_FACES);
   
    Param* restart_param = getParam("RESTART");
    if ( (restart_param ==NULL)||(restart_param->size() < 1)) { 
      CERR("> missing RESTART param.");
    }

    initMesh(restart_param,INIT_COMPACT_FACES);

    initData(); // needs to be populated.. 

    interpFromData();

    CtiRegister::dumpRegisteredData();
  }
  
  //==============================================
  // pure virtual functions from static solver... 
  // these are empty are not used in the init above...
  //==============================================
  void initBoundaryConditions() {}

  void run() {
    
    COUT1("PostproSolver::run()");
    
    temporalHook();

    writeResult();

    doProbes();

    flushProbes();

    processStep();

  }

  virtual void temporalHook() {

    COUT1("PostproSolver::temporalHook()");

  }
  
  
};

#endif
