#include "CTI.hpp"
using namespace CTI;

#include "StaticSolver.hpp"
#include "CtiRegister.hpp"

//
// the requirements of data registration are as follows:
// 
// 1. provide a common name-based registration environment for data
// 2. in some cases, connect data with dde's for input and/or output
// 3. function/expression evaluation involving constants, parameters (consider only REGISTER v 21), etc 
// 4. provide a framework for custom function evaluation in both the solver and solver bcs - 
//    here we try a CtiDataProducer framework. 
// 

class SampleRhs {
public:
  double a;
  double b;
};

class SampleState {
public:
  int it;
  double b;  
  double f[3];
  double a;
  bool b1;
  int g;
  double c;
  double d;
};

class SampleState2 {
public:
  int it;
  double a;
  bool b1;
  bool b2;
  bool b3;
  double gg[23];
  double b;  
  int g;
  double f[3];
  double c;
  double d;
};

class SampleSolverBaseBc {
public:
  BfZone * bfZone; // a pointer to the corresponding bfZone data for this bc
  SampleSolverBaseBc(BfZone * _bfZone) {
    this->bfZone = _bfZone; // just store a ptr
  }
  virtual ~SampleSolverBaseBc() {
  }
  // pull the name from bfZone...
  string getName() const {
    assert(bfZone);
    return bfZone->getName();
  }
  // some virtual methods that each instance must provide...
  virtual void initData() = 0;
  virtual void calcRhs() = 0;
  
};

class SampleSolver : public StaticSolver {
public:
  vector<SampleSolverBaseBc*> bcVec;
  SampleState * state;
  double time,dt;
  int step;
  double (*u)[3];
  double *p;
  SampleSolver() {
    state = NULL;
    u = NULL;
    p = NULL;
  }
  void initData() {
    COUT1("SampleWallBc::initData() " << getName());
  }
};

class SampleInletBc : public SampleSolverBaseBc {
private:
  SampleSolver * solver;
  double (*tau_wall)[3];
  double *p_wall;
public:
  SampleState2 * state;
  SampleInletBc(Param * param,BfZone * _bfZone,SampleSolver * _solver) : SampleSolverBaseBc(_bfZone), solver(_solver) {
    COUT1(" > SampleInletBc: " << getName());
  }
  ~SampleInletBc() {
    DELETE(state);
    DELETE(tau_wall);
    DELETE(p_wall);
  }
  void calcRhs() {
    COUT1("SampleInletBc::calcRhs()");
  }
  void initData() {
    COUT1("SampleInletBc::initData() " << getName());
    // "f"...
    if (!bfZone->checkDataFlag("f")) {
      COUT1(" > setting state[ibf].f to x_bf");
      for (int ibf = 0; ibf < bfZone->nbf; ++ibf) {
	FOR_I3 state[ibf].f[i] = bfZone->x_bf[ibf][i];
      }
    }
    else {
      COUT1(" > checking state[ibf].f with x_bf");
      for (int ibf = 0; ibf < bfZone->nbf; ++ibf) {
	FOR_I3 assert(state[ibf].f[i] == bfZone->x_bf[ibf][i]);
      }
    }
    // "c"...
    if (!bfZone->checkDataFlag("c")) {
      COUT1(" > setting state[ibf].c to 3.2*x_bf[2]");
      for (int ibf = 0; ibf < bfZone->nbf; ++ibf) {
	FOR_I3 state[ibf].c = 3.2*bfZone->x_bf[ibf][2];
      }
    }
    else {
      COUT1(" > checking state[ibf].c with 3.2*x_bf");
      for (int ibf = 0; ibf < bfZone->nbf; ++ibf) {
	FOR_I3 assert(state[ibf].c == 3.2*bfZone->x_bf[ibf][2]);
      }
    }
    // "tau_wall"...
    if (!bfZone->checkDataFlag("tau_wall")) {
      COUT1(" > setting tau_wall to 2.1*x_bf");
      for (int ibf = 0; ibf < bfZone->nbf; ++ibf) {
	FOR_I3 tau_wall[ibf][i] = 2.1*bfZone->x_bf[ibf][i];
      }
    }
    else {
      COUT1(" > checking tau_wall == 2.1*x_bf");
      for (int ibf = 0; ibf < bfZone->nbf; ++ibf) {
	FOR_I3 assert(tau_wall[ibf][i] == 2.1*bfZone->x_bf[ibf][i]);
      }
    }
    // "p_wall"...
    if (!bfZone->checkDataFlag("p_wall")) {
      COUT1(" > setting p_wall to 1.1*x_bf[0]");
      for (int ibf = 0; ibf < bfZone->nbf; ++ibf) {
	FOR_I3 p_wall[ibf] = 1.1*bfZone->x_bf[ibf][0];
      }
    }
    else {
      COUT1(" > checking p_wall == 1.1*x_bf[0]");
      for (int ibf = 0; ibf < bfZone->nbf; ++ibf) {
	FOR_I3 assert(p_wall[ibf] == 1.1*bfZone->x_bf[ibf][0]);
      }
    }
  }


};

void SampleSolver::init() {

  Param * restart_param = getParam("RESTART");
  if ((restart_param == NULL)||(restart_param->size() < 1)) {
    CERR("missing RESTART param.\n\nUsage: on the first run:\n\nRESTART myfile.mles\n\nor, to restart with data:\n\nRESTART myfile.mles snapshot.mles\n");
  }

  // read mles and lsp topology in sles
  
  initMesh(restart_param); 
 
  // ============================================

  // get our solver data registered...
  
  // process some other stuff that might produce registered data -- e.g. stats, probes, etc
  
  Param * stats_param = getParam("STATS");
  if (stats_param != NULL) {
    initStats(stats_param);
  }
  
  // take a look at the registered data...
  
  CtiRegister::dumpRegisteredData();

  // and read if there is a snapshot file to read...
  
  initData(restart_param);
  
  // now data is set -- check if the user wants to reset the stats...
  
  if (checkParam("RESET_STATS")) 
    CtiRegister::resetStats();
  
  // for now, check for the param "INTERACTIVE" -- this will stop the
  // solver in hold mode the first time processStep gets called...

  if (checkParam("INTERACTIVE"))
    kfr.setHold(true);

  // just to see if any init routines have put something in here... 

  CtiRegister::dumpCurrentData();

}

void SampleSolver::run() {
  
  // decide if we are going to set an initial condition...
  
  //initData();

  if (mpi_rank == 0) cout << "top of run loop: time=" << time << " dt=" << dt << " step=" << step << endl;

  int step_final = step + 10;
  int interval = 1;
  while (step < step_final) {
    // should put these in a "processStep" routine down one...
    
    time += dt;
    ++step;

    FOR_ICV FOR_I3 u[icv][i] = -time;
    FOR_ICV p[icv] = time;
    
    if (mpi_rank == 0) cout << 
      "\n==========================================================" <<
      "\n starting step: " << step << " time: " << time << " dt: " << dt << 
      "\n==========================================================" << endl;
    
    processStep(step,time,dt,interval);

  }
  
}

void SampleSolver::initData(Param * restart_param) {
  
  COUT1("SampleSolver::initData()");

  if (!CtiRegister::checkDataFlag("time")) {

    assert(!CtiRegister::checkDataFlag("step"));
    assert(!CtiRegister::checkDataFlag("dt"));
    
    // initialize 

    time = 0.0;
    dt = 1.0E-3;
    step = 0;
    
    COUT1(" > setting state[icv].f to x_vv[icv]");
    COUT1(" > setting state[icv].c to x_cv[icv][2]");

    assert(!CtiRegister::checkDataFlag("f"));
    assert(!CtiRegister::checkDataFlag("c"));
    
    assert(state);
    FOR_ICV {
      state[icv].it = icv;
      state[icv].a = x_cv[icv][0];
      state[icv].b1 = true;
      state[icv].b = x_cv[icv][1];
      state[icv].g = 2*icv;
      FOR_I3 state[icv].f[i] = x_vv[icv][i];
      state[icv].c = x_cv[icv][2];
      state[icv].d = 4.32124;
    }

 //   assert(!CtiRegister::checkDataFlag("u"));
   // COUT1(" > setting u[icv] to vol_cv[icv]*x_cv[icv]");
  //  FOR_ICV FOR_I3 u[icv][i] = x_cv[icv][i]*vol_cv[icv];
    FOR_ICV FOR_I3 u[icv][i] = -time;
    
 //   assert(!CtiRegister::checkDataFlag("p"));
 //   COUT1(" > setting p[icv] to vol_cv[icv]");
 //   FOR_ICV p[icv] = vol_cv[icv];
    FOR_ICV p[icv] = time;
  }
  else {

    assert(CtiRegister::checkDataFlag("step"));
    assert(CtiRegister::checkDataFlag("dt"));
    
    COUT1(" > checking that state[icv].f is set to x_vv[icv]");

 //   assert(CtiRegister::checkDataFlag("f"));
  //  assert(CtiRegister::checkDataFlag("c"));
    
    // check the state...
    FOR_ICV {
      //state[icv].it = icv;
      //state[icv].a = x_cv[icv][0];
      //state[icv].b1 = true;
      //state[icv].b = x_cv[icv][1];
      //state[icv].g = 2*icv;
      FOR_I3 assert(state[icv].f[i] == x_cv[icv][i]);
      assert(state[icv].c == x_cv[icv][2]);
      //state[icv].d = 4.32124;
    }

   // COUT1(" > checking that u[icv] is set to x_cv[icv]*vol_cv[icv]");

    assert(CtiRegister::checkDataFlag("u"));
    
    //FOR_ICV FOR_I3 assert(u[icv][i] == x_cv[icv][i]*vol_cv[icv]);
  //  COUT1(" > checking that p[icv] is set to vol_cv[icv]");

    assert(CtiRegister::checkDataFlag("p"));
    
   // FOR_ICV assert(p[icv] == vol_cv[icv]);
    
  }
  
  FOR_IZONE(bcVec) {
    bcVec[izone]->initData();
  }

}

void * T_ptr;
int offset;
int stride;
int * nptr;

template<class T>
void registerState(T * &state,double &var,const string& name,int &n) {

  T_ptr = (void*)&state;
  double * vptr = &var;
  nptr = &n;

  offset = int(vptr-(*(double**)T_ptr));
  assert(sizeof(T)%8 == 0);
  stride = sizeof(T)/8;

  assert(state != NULL);
  cout << "in registerState: T_ptr: " << T_ptr << " *T_ptr: " << (*(double**)T_ptr) << " vptr: " << vptr << " (vptr-(*(double**)T_ptr)): " << (vptr-(*(double**)T_ptr)) << endl;
  cout << " > offset: " << offset << endl;
  cout << " > stride: " << stride << endl;
  
}

template<class T>
void registerState(T * &state,double (&var)[3],const string& name,int &n) {

  T_ptr = (void*)&state;
  double * vptr = (double*)&var;
  nptr = &n;

  offset = int(vptr-(*(double**)T_ptr));
  assert(sizeof(T)%8 == 0);
  stride = sizeof(T)/8;

  assert(state != NULL);
  cout << "in registerState: T_ptr: " << T_ptr << " *T_ptr: " << (*(double**)T_ptr) << " vptr: " << vptr << " (vptr-(*(double**)T_ptr)): " << (vptr-(*(double**)T_ptr)) << endl;
  cout << " > offset: " << offset << endl;
  cout << " > stride: " << stride << endl;
}

void checkRegisteredStateDN() {

  cout << "checkRegisteredState: T_ptr: " << T_ptr << " *T_ptr: " << (*(double**)T_ptr) << " n: " << (*nptr) << endl;
  
  cout << "first 3 values: " << endl;
  const int n = *nptr;
  for (int i = 0; i < min(n,3); ++i) {
    cout << "i: " << i << " value: " << *(*(double**)T_ptr+offset+stride*i) << endl;
  }
  
  cout << "last 3 values: " << endl;
  for (int i = max(0,n-3); i < n; ++i) {
    cout << "i: " << i << " value: " << *(*(double**)T_ptr+offset+stride*i) << endl;
  }
  
}

void checkRegisteredState() {

  cout << "checkRegisteredState: T_ptr: " << T_ptr << " *T_ptr: " << (*(double**)T_ptr) << " n: " << (*nptr) << endl;
  
  cout << "first 3 values: " << endl;
  const int n = *nptr;
  for (int i = 0; i < min(n,3); ++i) {
    FOR_J3 cout << "i: " << i << " j: " << j << " value: " << *(*(double**)T_ptr+offset+stride*i+j) << endl;
  }
  
  cout << "last 3 values: " << endl;
  for (int i = max(0,n-3); i < n; ++i) {
    FOR_J3 cout << "i: " << i << " j: " << j << " value: " << *(*(double**)T_ptr+offset+stride*i+j) << endl;
  }
  
}

void setRegisteredState() {

  cout << "setRegisteredState: T_ptr: " << T_ptr << " *T_ptr: " << (*(double**)T_ptr) << " n: " << (*nptr) << endl;
  
  const int n = *nptr;
  for (int i = 0; i < n; ++i) {
    *(*(double**)T_ptr+offset+stride*i) = 12345.678+i;
  }

}

  

int main(int argc, char * argv[]) {

  try {
    
    CTI_Init(argc,argv,"testRegister.in");

    {

      /*
	class SampleState {
	public:
	int it;	
	double b;  
	double f[3];
	double a;
	bool b1;
	int g;
	double c;
	double d;
	};
      */
      
      int np = 0;
      SampleState * state = new SampleState[np];
      //registerState<SampleState>(state,state->a,"a",np);
      registerState<SampleState>(state,state->f,"f",np);

      cout << "STEP1: 10 values" << endl;

      for (int i = 0; i < np; ++i) {
	//state[i].a = 101.01+i; cout << "setting state[" << i << "].a to " << 101.01+i << endl;  
	FOR_J3 state[i].f[j] = 101.01+i+j; 
	FOR_J3 cout << "setting state[" << i << "].f[" << j << "] to " << 101.01+i+j << endl;  
      }

      checkRegisteredState();
      setRegisteredState();
      checkRegisteredState();

      cout << "STEP2: nullify state" << endl;

      np = 0;
      SampleState * state2 = state;
      state = NULL;
      
      checkRegisteredState();
      setRegisteredState();
      checkRegisteredState();
      
      cout << "STEP3: 13 values" << endl;
      
      np = 13;
      state = new SampleState[np];
      for (int i = 0; i < np; ++i) {
	//state[i].a = 201.01+i; cout << "setting state[" << i << "].a to " << 201.01+i << endl;  
	FOR_J3 state[i].f[j] = 201.01+i+j; 
	FOR_J3 cout << "setting state[" << i << "].f[" << j << "] to " << 201.01+i+j << endl;  
      }

      checkRegisteredState();
      setRegisteredState();
      checkRegisteredState();


      delete[] state;
      delete[] state2;

      // old check...
      
      /*
	SampleSolver solver;
	solver.init();
	solver.run();
      */

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

