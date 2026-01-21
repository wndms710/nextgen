#include "Lsp.hpp"
#include "CtiSolid.hpp"
#include "CtiLiquid.hpp"

#include "StaticSolver.hpp"
using namespace CTI;

#ifndef __LSPBC__HPP__
#define __LSPBC__HPP__


//class IdealGasSolverWithLsp;
class FlowSolver; // TODO template on solvers? create auxilary class to inherit from that holds lsp stuff?

class lpBaseBc {
public:
  BfZone * bf_zone;
  FlowSolver * solver;

  int nbuf;           // number of doubles passed in applyBc
  int np_hit;         // number of particles hit this bc
  BfZone* zone_ptr;   // Pointer to the zone object
  lpBaseBc(BfZone* _zone_ptr, FlowSolver * _solver) {
    nbuf      = -1;
    np_hit    = 0;
    zone_ptr  = _zone_ptr;
    solver    = _solver;
  }

  virtual ~lpBaseBc();

  virtual void initialHook() {}

  virtual void applyBc(LspState& lp, double* xi, double* st_normal, double* buf) = 0;
  virtual bool applyBc2(LspState& lp, double* dx, double* st_normal) = 0;
  virtual void queryBc();
  virtual void updateStats() { 
    // Reset stats
    //if (solver->step % solver->lsp_stats_interval == 0) 
    np_hit  = 0;
  }
};

class lpOutletBc : public lpBaseBc {
  public:

  lpOutletBc(BfZone* p, FlowSolver *s): lpBaseBc(p,s) {
    COUT1(" > Assigning OUTLET lsp bc to zone " << zone_ptr->getName());
  }
  virtual ~lpOutletBc() {}
  void applyBc(LspState& lp, double* xi, double* st_normal, double* buf);
  bool applyBc2(LspState& lp, double* dx, double* st_normal);

};

class lpPerfectBounceBc : public lpBaseBc {
  public:
  
  lpPerfectBounceBc(BfZone* p, FlowSolver *s): lpBaseBc(p,s) {
    COUT1(" > Assigning BOUNCE lsp bc to zone " << zone_ptr->getName());
  }
  virtual ~lpPerfectBounceBc() {}
  void applyBc(LspState& lp, double* xi, double* st_normal, double* buf);
  bool applyBc2(LspState& lp, double* dx, double* st_normal);

};

class lpIrregularBounceBc : public lpBaseBc {
  public:
  double resCoef;             // restitution coefficient
  double mu_s;                // static friction coef
  double mu_k;                // kinetic friction coef
  double roughness_max_angle; // maximum change in angle due to roughness

  lpIrregularBounceBc(BfZone* p, FlowSolver *s) : lpBaseBc(p,s) { 
    // Assign negative (invalid) values to these parameters so they fail if not properly assigned
    resCoef = -1; mu_s = -1; mu_k = -1; roughness_max_angle = -1;
  }

  virtual ~lpIrregularBounceBc() {}
  void initialHook();
  void applyBc(LspState& lp, double* xi, double* st_normal, double* buf);
  bool applyBc2(LspState& lp, double* dx, double* st_normal);
};

class lpInelasticBounceBc : public lpBaseBc {
  public:
  double mu;                        // friction coef
  double CoR_min, CoR_max, CoR_bar; // CoR stats
  CtiSolid * solid;                 // solid material

  lpInelasticBounceBc(BfZone* p, FlowSolver *s) : lpBaseBc(p,s) { 
    // Assign negative (invalid) values to these parameters so they fail if not properly assigned
    mu = -1.0; solid = NULL;
  }

  virtual ~lpInelasticBounceBc() {}
  void initialHook();
  void applyBc(LspState& lp, double* xi, double* st_normal, double* buf);
  bool applyBc2(LspState& lp, double* dx, double* st_normal);
  void updateStats();
  void queryBc();
};

#endif
