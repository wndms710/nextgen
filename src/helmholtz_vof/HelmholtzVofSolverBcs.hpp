#ifndef HELMHOLTZVOFSOLVERBCS_HPP
#define HELMHOLTZVOFSOLVERBCS_HPP

#include "BcCommon.hpp"
#include "States.hpp"
#include "BoundaryLayerDataExchanger2.hpp"

// need to provide a forward declaration of the solver below...
class HelmholtzVofSolver;

class VofBc : public CtiRegister::CtiDataProducer {
public:
  
  BfZone* zone_ptr;
  HelmholtzVofSolver* solver;
  double* q_bf;
  double (*u_bc)[3];
  std::stringstream *ss;
  bool b_write;

  VofBc(BfZone* p, HelmholtzVofSolver* s): zone_ptr(p), solver(s), q_bf(NULL) {
    addCtiDataProducer(this);
    ss = NULL;
    q_bf = NULL;
    u_bc = NULL;
    b_write = false;
  }

  virtual ~VofBc() {
    DELETE(q_bf);
    DELETE(u_bc);
    if (ss != NULL) delete ss;
  }

  // initialization of boundary condition memory, parsing of bc params...

  virtual void initData() = 0;
  virtual void initialHook() {}

  // perform any pre-time step state copy down and set the (extrapolated) 
  // guess for the boundary conditions at the next time level...

  virtual void setBc() = 0;

  // hooks to apply the boundary conditions (weakly) onto the flow solution

  virtual void addFlux(double * rhs) const = 0;
  virtual void addMassFlux(double *rhs_rhoY0, double *rhs_rhoY1) const = 0;
  virtual void addMomentumFlux(double * A, double (*rhs)[3]) const = 0;

  // some of the boundary conditions (for temporal accuracy/stability) require an update
  // of their mass flux based on their intermediate predicted state in the fractional
  // step algorithm

  virtual void updateBc() {}

  // all boundary conditions require the ability to compute a force...

  virtual void force_bf( double (*f_bf)[9]) const = 0;

  // each boundary condition can also define a query which will report
  // relevant information about the bc (dumpRange, integrated values) depending
  // on what is relevant for its particular case.  the default behavior is
  // that the query is empty..
  
  virtual void query(const string& param_str) {}

  void flush() {
    if (mpi_rank == 0) {
      assert(ss);
      if (b_write) {
        ofstream out_file;
        char filename[128];
        sprintf(filename,"%s.bc",zone_ptr->getName().c_str());
        out_file.open(filename,ofstream::app);
        assert( out_file.is_open());

        out_file << ss->rdbuf();
        out_file.close();
      }
      else {
        // just pipe to cout
        cout << ss->rdbuf();
      }
      ss->str(string()); // clears ss
    }
  }

 
  virtual CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v, 
      const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

    if ( zone_ptr->isLeadingStrMatched(name)) { 

      const string proj_str = zone_ptr->getName() + ":" + "proj";
      //const string xbf_str = zone_ptr->getName() + ":" + "x_bf"; // should work now w/o this - CI
     
      if ( name == proj_str) { 

        if ( args.size() != 1) 
          return CtiRegister::CTI_DATA_ARG_COUNT;

        list<CtiRegister::CtiData>::iterator arg = args.begin();
        if ( arg->getTopology() != CV_DATA) 
          return CtiRegister::CTI_DATA_NOT_VALID;

        if ( arg->getType() == DN_DATA) { 

          double * tmp = zone_ptr->createBfD1Data(v);

          if ( b_eval_func) { 

            for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
              
              const int icv = zone_ptr->cvobf[ibf];
              tmp[ibf]      = arg->dn(icv);
              
            }
          }

          return CtiRegister::CTI_DATA_OK;

        } else if ( arg->getType() == DN3_DATA) { 

          double (*tmp)[3] = zone_ptr->createBfD2Data(v);

          if ( b_eval_func) { 

            for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

              const int icv = zone_ptr->cvobf[ibf];
              for (int i = 0; i < 3; ++i) 
                tmp[ibf][i] = arg->dn3(icv,i);

            }
          }

          return CtiRegister::CTI_DATA_OK;

        } else { 

          return CtiRegister::CTI_DATA_NOT_VALID;

        }

      } 
      /*
      else if ( name == xbf_str) { 

        if ( args.size() != 0) 
          return CtiRegister::CTI_DATA_ARG_COUNT;

        double (*tmp)[3] = zone_ptr->createBfD2Data(v);
        
        if ( b_eval_func) { 
          
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
            
            for (int i = 0; i < 3; ++i) 
              tmp[ibf][i] = zone_ptr->x_bf[ibf][i]; 
            
          }
        }
        
        return CtiRegister::CTI_DATA_OK;

      }
      */
    }

    return CtiRegister::CTI_DATA_NOT_FOUND;

  }

  string getName() const {
    assert(zone_ptr);
    return zone_ptr->getName();
  }

  void force(double (*f_buf)[3],double (*m_buf)[3]) const {
    
    // f_buf: pressure, viscous, convective force components

    FOR_J3 FOR_I3 f_buf[i][j] = 0.0; 

    // m_buf: compute moment about origin

    FOR_J3 FOR_I3 m_buf[i][j] = 0.0;

    double (*f_bf)[9] = new double[zone_ptr->nbf][9];
    force_bf(f_bf);

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      
      FOR_I3 FOR_J3 f_buf[i][j] += f_bf[ibf][3*i+j];
      FOR_I3 {
        const double m_bf[3] = CROSS_PRODUCT(zone_ptr->x_bf[ibf],&f_bf[ibf][3*i]);
        FOR_J3 m_buf[i][j] += m_bf[j];
      }

    }

    delete[] f_bf;

  }

};

class SlipWallVBc : public VofBc {
public:

  SlipWallVBc(BfZone* p, HelmholtzVofSolver* s): VofBc(p,s) {}

  void initData() {}

  // slip walls do not provide any real information; the pressure
  // contribution fro the wall pressure is subsumed through the
  // pressure gradient...

  void setBc() {}
  void addFlux(double *rhs) const {}
  void addMassFlux(double *rhs_rhoY0, double *rhs_rhoY1) const {}
  void addMomentumFlux(double* A, double (*rhs)[3]) const {}
  void force_bf(double (*f_bf)[9]) const;

};

class WallVBc : public VofBc {
public:

  double (*u_bc)[3];

  WallVBc(BfZone* p, HelmholtzVofSolver* s): VofBc(p,s), u_bc(NULL) {
    funcEvalList.push_back(p->getName()+":tau_wall");
    funcEvalList.push_back(p->getName()+":y_plus");
  }

  ~WallVBc() {
    DELETE(u_bc);
  }

  void initData();

  void addFlux(double *rhs) const {}
  void addMassFlux(double *rhs_rhoY0, double *rhs_rhoY1) const {}
  void setBc();
  void addMomentumFlux(double* A, double (*rhs)[3]) const;
  void force_bf(double (*f_bf)[9]) const;
  void query(const string& param_str);

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v, 
    const string& name, list<CtiRegister::CtiData>& args,const bool b_eval_func);

};

class InletVBc : public VofBc {
public:

  double (*u_bc)[3];
  double rho0_bc;
  double rho1_bc;
  double vof_bc;

  InletVBc(BfZone* p, HelmholtzVofSolver* s): VofBc(p,s), u_bc(NULL){}

  ~InletVBc() {
    DELETE(u_bc);
  }
  
  void initData();

  void setBc();
  void addFlux(double *rhs) const;
  void addMassFlux(double *rhs_rhoY0, double *rhs_rhoY1) const;
  void addMomentumFlux(double* A, double (*rhs)[3]) const;
  void force_bf(double (*f_bf)[9]) const;
  void query(const string& param_str);
};

class OutletVVBc : public VofBc {
public:
  
  OutletVVBc(BfZone* p, HelmholtzVofSolver* s): VofBc(p,s) {
    q_bf = NULL; zone_ptr->registerBfData(q_bf,"q_bf",READWRITE_DATA);
  }

  double rho0_bc;
  double rho1_bc;

  void initData();
  
  void setBc();
  void updateBc();
  void addFlux(double * rhs) const;
  void addMassFlux(double *rhs_rhoY0, double *rhs_rhoY1) const;
  void addMomentumFlux(double * A, double (*rhs)[3]) const;
  void force_bf(double (*f_bf)[9]) const;

};

class SpongeVVBc : public VofBc {
public:
  
  SpongeVVBc(BfZone* p, HelmholtzVofSolver* s): VofBc(p,s) {
    q_bf = NULL; zone_ptr->registerBfData(q_bf,"q_bf",READWRITE_DATA);
  }

  double rho0_bc;
  double rho1_bc;

  //sponge
  double (*u_sponge)[3];
  double sponge_data[5];
  double sponge_speed;
  double sponge_strength;
  SpongeType sponge_type;
  
  double th_e0[3];
  double th_e1[3];

  double p_ref;
  double Ma;
  double L;
  double sigma;

  double * xr;
  double * nr;

  void initData();
  
  void setBc();
  void updateBc();
  void addFlux(double * rhs) const;
  void addMassFlux(double *rhs_rhoY0, double *rhs_rhoY1) const;
  void addMomentumFlux(double * A, double (*rhs)[3]) const;
  void force_bf(double (*f_bf)[9]) const;
  void initSpongeTypeTheta();

};



class OutletVBc : public VofBc {
public:
  
  OutletVBc(BfZone* p, HelmholtzVofSolver* s): VofBc(p,s) {
    q_bf = NULL; zone_ptr->registerBfData(q_bf,"q_bf",READWRITE_DATA);
  }

  double rho0_bc;
  double rho1_bc;

  void initData();
  
  void setBc();
  void updateBc();
  void addFlux(double * rhs) const;
  void addMassFlux(double *rhs_rhoY0, double *rhs_rhoY1) const;
  void addMomentumFlux(double * A, double (*rhs)[3]) const;
  void force_bf(double (*f_bf)[9]) const;

};

class SpongeVBc : public VofBc {
public:
  
  SpongeVBc(BfZone* p, HelmholtzVofSolver* s): VofBc(p,s) {
    q_bf = NULL; zone_ptr->registerBfData(q_bf,"q_bf",READWRITE_DATA);
  }

  double rho0_bc;
  double rho1_bc;

  void initData();
  
  void setBc();
  void updateBc();
  void addFlux(double * rhs) const;
  void addMassFlux(double *rhs_rhoY0, double *rhs_rhoY1) const;
  void addMomentumFlux(double * A, double (*rhs)[3]) const;
  void force_bf(double (*f_bf)[9]) const;

};

class AlgebraicWallModelVBc : public VofBc {
public:

  double * tau_wall;
  double z0;
  double (*u_bc)[3];

  BoundaryLayerDataExchanger<HelmholtzVofSolver> * blde;

  AlgebraicWallModelVBc(BfZone* p, HelmholtzVofSolver* s): VofBc(p,s), u_bc(NULL), tau_wall(NULL), z0(-1), blde(NULL) {

    zone_ptr->registerBfData(tau_wall,"tau_wall",READWRITE_DATA);
    funcEvalList.push_back(p->getName()+":tau_wall");
    funcEvalList.push_back(p->getName()+":y_plus");
    funcEvalList.push_back(p->getName()+":bl_delta()");
  }

  ~AlgebraicWallModelVBc() {
    DELETE(u_bc);
    DELETE(tau_wall);
    if (!blde) {
      delete blde;
      blde = NULL;
    }
 }

  void initData();
  void initialHook();

  void setBc() {}
  void addFlux(double *rhs) const {}
  void addMassFlux(double *rhs_rhoY0, double *rhs_rhoY1) const {}
  void addMomentumFlux(double* A, double (*rhs)[3]) const;
  void force_bf(double (*f_bf)[9]) const;
  void query(const string& param_str);

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v, 
    const string& name, list<CtiRegister::CtiData>& args,const bool b_eval_func);

};

class HookVBc : public VofBc {
public:

  double (*u_bc)[3];
  double rho0_bc;
  double rho1_bc;
  double vof_bc;

  HookVBc(BfZone* p, HelmholtzVofSolver* s): VofBc(p,s), u_bc(NULL){}

  ~HookVBc() {
    DELETE(u_bc);
  }
  
  void initData();

  void setBc();
  void addFlux(double *rhs) const;
  void addMassFlux(double *rhs_rhoY0, double *rhs_rhoY1) const;
  void addMomentumFlux(double* A, double (*rhs)[3]) const;
  void force_bf(double (*f_bf)[9]) const;
  void query(const string& param_str);
};


#endif
