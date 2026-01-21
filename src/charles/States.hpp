#ifndef STATES_HPP
#define STATES_HPP

//=======================================
// state class definitions
//=======================================
#ifndef STATE_ALIGN_SIZE
#define STATE_ALIGN_SIZE 64
#endif

/*
#ifndef NO_STATE_ALIGN
#define cti_state_align __attribute__((aligned(STATE_ALIGN_SIZE)))
#else
#define cti_state_align
#endif
*/

#define cti_state_align

enum SolverCommActions { 
  COMPLETE_ALL_UPDATES,
  POSTPONE_GRAD_GHOSTS
};

class IdealGasState { 
public: 
  double u[3];
  double sp_vol;
  double p;
  double h;
} cti_state_align;

// additional gas properties that are used 
// for the calculation of the viscous fluxes 
// for instance.. 

class AuxGasProp { 
public: 
  double mu_total; 
  double loc_total;
  double half_usq;
  double rho;
  double a_sgs;
  double sos;
  double fgr;
  double symh_fax;
  double symp_fax;
} cti_state_align;

class IdealGasRhs {
public:
  double rho;
  double rhou[3];
  double rhoE;

  void zero() { 
    rho = 0.0;
    for (int i =0; i < 3; ++i) 
      rhou[i] = 0.0;
    rhoE = 0.0;
  }
};

class PremixedState { 
public: 
  double u[3];
  double sp_vol;
  double p;
  double h;
  double Z;
  double C;
} cti_state_align;

class PremixedRhs {
public:
  double rho;
  double rhou[3];
  double rhoE;
  double rhoZ;
  double rhoC;

  void zero() {
    rho = 0.0;
    for (int i =0; i < 3; ++i) 
      rhou[i] = 0.0;
    rhoE = 0.0;
    rhoZ = 0.0;
    rhoC = 0.0;
  }
};


// for the nonpremixed model with algebraic closures, the
// transpored vars are the same as the PremixedState, so 
// we;ll use a typedef here.  if this changes, ie non-alg
// closure of Zvar, fill this in.

typedef PremixedState NonpremixedState;
typedef PremixedRhs NonpremixedRhs;

#endif
