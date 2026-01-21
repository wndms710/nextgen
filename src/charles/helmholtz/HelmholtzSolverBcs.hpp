#ifndef HELMHOLTZSOLVERBCS_HPP
#define HELMHOLTZSOLVERBCS_HPP

#include "BcCommon.hpp"
#include "States.hpp"
#include "BoundaryLayerDataExchanger.hpp"

// need to provide a forward declaration of the solver below ..
class HelmholtzSolver;

class HelmholtzBc : public CtiRegister::CtiDataProducer {
public:

  BfZone* zone_ptr;
  HelmholtzSolver* solver;
  double *mf;
  double *scbc_vals;
  double (*u_bc)[3];
  std::stringstream *ss;
  bool b_write;

  bool b_helmholtzCorrection;

  HelmholtzBc(BfZone* p, HelmholtzSolver* s): zone_ptr(p), solver(s) {
    mf = NULL;
    scbc_vals = NULL;
    u_bc = NULL;
    addCtiDataProducer(this);
    ss = NULL;
    b_write = false;

    b_helmholtzCorrection = checkParam("HELMHOLTZ_CORRECTION");

    funcEvalList.push_back(p->getName()+":x_bf()");
    funcEvalList.push_back(p->getName()+":proj()");
    funcEvalList.push_back(p->getName()+":p_bf()");
  }
  HelmholtzBc(BfZone* p, HelmholtzSolver* s,const int icg): zone_ptr(p), solver(s) {
    mf = NULL;
    scbc_vals = NULL;
    u_bc = NULL;
    ss = NULL;
    b_write = false;
    b_helmholtzCorrection = checkParam("HELMHOLTZ_CORRECTION");
  }

  virtual ~HelmholtzBc() {
    DELETE(mf);
    DELETE(scbc_vals);
    DELETE(u_bc);
    if (ss != NULL) delete ss;
  }

  // reset of the bc will re-clear the state information (any time dependent
  // states that are recorded on the boundary)

  virtual void resetBc() = 0;
  virtual void restrictBc(StaticSolver::CoarseGrid* cg,HelmholtzBc* bc); //default defined in *.cpp

  // initialization of boundary condition memory, parsing of bc params...

  virtual void initData() = 0;
  virtual void initFromCoarseGrid(StaticSolver::CoarseGrid* cg) {initData();}
  virtual void initialHook() = 0;

  virtual void modifyU() {} //used for exerimental inflow turbulence

  // used to initialize certain boundary data when in post mode
  // initialHook/temporalHook are not called, so data that was not read in or
  // set in the contructor or initData will need to be set (i.e. mf for inlets, dswm)

  virtual void syncPostState() = 0; 

  // hooks to apply the boundary conditions (weakly) onto the flow solution

  virtual void addMomentumFlux(double *A,double (*rhs)[3]) const = 0;
  virtual void addMomentumFlux(double *A,double *At,double (*rhs)[3]) const = 0;
  virtual void addPressureFlux(double *A,double *rhs) const = 0;
  virtual void addMassFlux(double *rhs) const = 0;

  // the following signature is useful if the pressure gradient calc requires
  // closure of the pressure at the boundary faces ...

  virtual void completePressureGrad(double (*dpdx)[3]) const;  //default defined in *.cpp

  // perform any pre-time step state copy down and set the (extrapolated)
  // guess for the boundary conditions at the next time level...
  virtual void setBc() = 0;

  // some boundary conditions (for temporal accuracy/stability) require an update
  // of their mass flux based on their intermediate predicted (velocity) state in the fractional
  // step algorithm
  virtual void updateMf() = 0;
   
  // boundary conditions who apply a non-zero dpdx during solveHelmholtz must have their
  // mass flux corrected based on that dpdx
  virtual void correctMf() = 0;

  // update the outlet pressure after solving the helmholtz equation

  virtual void updateP() = 0;

  // all boundary conditions require the ability to compute a force ..

  virtual void force_bf( double (*f_bf)[9]) const = 0;
  void addPressureForceDueToFrameRotation(double (*f_bf)[9]) const;

  // each boundary condition can also define a query which will report
  // relevant information about the bc (dumpRange, integrated values) depending
  // on what is relevant for its particular case.  the default behavior is
  // that the query is empty..

  virtual void query(const string& param_str) const; //default defined in *.cpp, throws error

  void flush() const {
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
                                                    const string& name, list<CtiRegister::CtiData>& args,
                                                    const bool b_eval_func); 

  string getName() const {
    assert( zone_ptr) ;
    return zone_ptr->getName();
  }

  void force(double (*f_buf)[3],double (*m_buf)[3]) const {

    //f_buf: pressure, viscous, convective force components

    FOR_J3 FOR_I3 f_buf[i][j] = 0.0;

    //m_buf: compute moment about origin

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

  //-------------------------------------
  // passive scalar transport
  //-------------------------------------

  virtual void addScalarBoundaryFlux(double *A,double *rhs,const ScalarField *scalar) const = 0;
  virtual void addScalarBoundaryFlux(double *A,double *At,double *rhs,const ScalarField *scalar) const = 0;
  void parseScalarBc(Param * param, const bool b_required);
  void _addScalarBoundaryFluxConvective(double* A,double* At,double * rhs,const ScalarField *scalar) const;
  void _addScalarBoundaryFluxWall(double* A,double* At,double * rhs,const ScalarField *scalar) const;

};

class SlipWallHBc : public HelmholtzBc {
public:

  SlipWallHBc(BfZone* p, HelmholtzSolver* s): HelmholtzBc(p,s) {}
  SlipWallHBc(BfZone* p, HelmholtzSolver* s,const int icg): HelmholtzBc(p,s,icg) {}

  //empty overrides of pure virtual functions
  void resetBc() {}
  void syncPostState() {}
  void initialHook() {}
  void addMomentumFlux(double * A, double (*rhs)[3]) const {}
  void addMomentumFlux(double *A,double *At,double (*rhs)[3]) const {}
  void addPressureFlux(double *A,double *rhs) const {}
  void addMassFlux(double * rhs) const {}
  void setBc() {}
  void updateMf() {}
  void correctMf() {}
  void updateP() {}


  void initData();
  void force_bf(double (*f_bf)[9]) const;

  void addScalarBoundaryFlux(double *A,double *rhs,const ScalarField *scalar) const {
    _addScalarBoundaryFluxWall(A,NULL,rhs,scalar);
  }
  void addScalarBoundaryFlux(double *A,double *At,double *rhs,const ScalarField *scalar) const {
    _addScalarBoundaryFluxWall(A,At,rhs,scalar);
  }


};

class WallHBc : public HelmholtzBc {
public:

  BoundaryLayerDataExchanger<HelmholtzSolver> * blde;

  WallHBc(BfZone* p, HelmholtzSolver* s): HelmholtzBc(p,s), blde(NULL) {

    funcEvalList.push_back(p->getName()+":tau_wall()");
    funcEvalList.push_back(p->getName()+":y_plus()");
    funcEvalList.push_back(p->getName()+":bl_delta()");

  }
  WallHBc(BfZone* p, HelmholtzSolver* s,const int icg): HelmholtzBc(p,s,icg), blde(NULL) {}

  ~WallHBc() {

    if (!blde) {
      delete blde;
      blde = NULL;
    }

  }

  //empty overrides of pure virtual functions
  void resetBc() {}
  void initialHook() {}
  void syncPostState() {}
  void addPressureFlux(double *A,double *rhs) const {}
  void addMassFlux(double * rhs) const {}
  void setBc() {}
  void updateMf() {}
  void correctMf() {}
  void updateP() {}

  void initData();
  void addMomentumFlux(double * A, double (*rhs)[3]) const;
  void addMomentumFlux(double * A,double *At,double (*rhs)[3]) const;
  void force_bf(double (*f_bf)[9]) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
    const string& name, list<CtiRegister::CtiData>& args,const bool b_eval_func);

  void query(const string& param_str) const;

  void addScalarBoundaryFlux(double *A,double *rhs,const ScalarField *scalar) const {
    _addScalarBoundaryFluxWall(A,NULL,rhs,scalar);
  }
  void addScalarBoundaryFlux(double *A,double *At,double *rhs,const ScalarField *scalar) const {
    _addScalarBoundaryFluxWall(A,At,rhs,scalar);
  }

};

class InletHBc : public HelmholtzBc {
public:

  InletHBc(BfZone* p, HelmholtzSolver* s): HelmholtzBc(p,s) {}
  InletHBc(BfZone* p, HelmholtzSolver* s,const int icg): HelmholtzBc(p,s,icg) {}

  //empty overrides of pure virtual functions
  virtual void resetBc() {}
  virtual void initialHook() {};
  void addPressureFlux(double *A,double *rhs) const {}
  virtual void updateMf() {}
  virtual void correctMf() {}
  void updateP() {}
  virtual void syncPostState() {};

  virtual void setBc();
  virtual void initData();
  virtual void initFromCoarseGrid(StaticSolver::CoarseGrid* cg);

  void addMomentumFlux(double * A, double (*rhs)[3]) const;
  void addMomentumFlux(double * A,double *At,double (*rhs)[3]) const;
  void addMassFlux(double * rhs) const;

  void addScalarBoundaryFlux(double *A,double *rhs,const ScalarField *scalar) const {
    _addScalarBoundaryFluxConvective(A,NULL,rhs,scalar);
  }
  void addScalarBoundaryFlux(double *A,double *At,double *rhs,const ScalarField *scalar) const {
    _addScalarBoundaryFluxConvective(A,At,rhs,scalar);
  }

  void force_bf(double (*f_bf)[9]) const;

  void query(const string& param_str) const;

};

class InletHBcProfile : public InletHBc {
public:
  InletHBcProfile(BfZone* p, HelmholtzSolver* s) : InletHBc(p,s) {}
  InletHBcProfile(BfZone* p, HelmholtzSolver* s,const int icg) : InletHBc(p,s,icg) {}
  void initData();
  void initFromCoarseGrid(StaticSolver::CoarseGrid* cg);
};

class InflowTurbulenceHBc : public InletHBc {
public:
  
  // inflow turbulence fields
  double (*up_it)[3];
  double (*um_it)[3];
  double (*Rd_it)[3];
  double (*Rod_it)[3];

  // turbulence lengthscale
  double* length_scale;
  double max_length; // this is a maximum in plane lengt scale
  double in_plane_length; // this is a constant length used in the whole plane
  double length_ratio; // ratio of streamwise to in-plane lengthscale
  double Un; // this is a constant convective velocity

  double Umean_it; //these are used to scale the turbulent fluctuations
  double (*avg)[3];
  double (*rms)[3];
  double IturbWgt;
  bool reset;
  bool b_modify_u;

  // connectivity and equation solving related variables
  // use an underscore convention for counting the local active
  int* activeFaces;
  int* activeCvs;
  int ncv_,ncv_g_;
  int nfa_i_,nfa_;
  int (*cvofa_)[2];
  int *cvobf_;
  double * A_fa_;
  double * A_diag_cv_;
  double inflow_zero; // Threshold for convergence of inflow turbulence solver
  int inflow_maxiter;

  InflowTurbulenceHBc(BfZone* p, HelmholtzSolver* s) : InletHBc(p,s) {
    
    up_it         = NULL; zone_ptr->registerBfData(up_it,         "up_it",       READWRITE_DATA); //fluctuating velocity, u'
    um_it         = NULL; zone_ptr->registerBfData(um_it,         "um_it",       READWRITE_DATA); //mean velocity
    Rd_it         = NULL; zone_ptr->registerBfData(Rd_it,         "Rd_it",       READWRITE_DATA);
    Rod_it        = NULL; zone_ptr->registerBfData(Rod_it,        "Rod_it",      READWRITE_DATA);
    length_scale  = NULL; zone_ptr->registerBfData(length_scale, "length_scale", READWRITE_DATA);
    avg           = NULL; zone_ptr->registerBfData(avg,           "avgIT",       READWRITE_DATA);
    rms           = NULL; zone_ptr->registerBfData(rms,           "rmsIT",       READWRITE_DATA);
    IturbWgt      = 0.0;  zone_ptr->registerBfData(IturbWgt,      "IturbWgt",    READWRITE_DATA);
    reset=false;
    b_modify_u=false;

    //  persistent operators...
    activeFaces = NULL;
    activeCvs = NULL;
    cvofa_ = NULL;
    A_fa_ = NULL;
    A_diag_cv_ = NULL;
    cvobf_ = NULL;
  
  }

  InflowTurbulenceHBc(BfZone* p, HelmholtzSolver* s,const int icg) : InletHBc(p,s,icg) {
    
    up_it         = NULL;  
    um_it         = NULL;  
    Rd_it         = NULL; 
    Rod_it        = NULL; 
    length_scale  = NULL; 
    avg           = NULL; 
    rms           = NULL; 
    IturbWgt      = 0.0;  
    reset=false;
    b_modify_u=false;

    //  persistent operators...
    activeFaces = NULL;
    activeCvs = NULL;
    cvofa_ = NULL;
    A_fa_ = NULL;
    A_diag_cv_ = NULL;
    cvobf_ = NULL;
    
  }

  ~InflowTurbulenceHBc(){
    DELETE(up_it);
    DELETE(um_it);
    DELETE(Rd_it);
    DELETE(Rod_it);
    DELETE(length_scale);
    DELETE(avg);
    DELETE(rms);

    DELETE(cvofa_);
    DELETE(A_fa_);
    DELETE(A_diag_cv_);
    DELETE(cvobf_);
    DELETE(activeFaces);
    DELETE(activeCvs);
  }

  //empty overrides of pure virtual functions
  void updateMf() {}
  
  void resetBc();
  void initData();
  void initFromCoarseGrid(StaticSolver::CoarseGrid* cg);
  void syncPostState();

  void initialHook();
  void setBc();
  void modifyU();

  //additional functions in InflowTurbulenceH.hpp
  void setStatistics();
  void updateInflowTurbulence(double (*FacePtr)[3]);
  void buildITBoundaryStructures();
  void setLengthscale();
  void buildFilterOperators();
  void initiateStats();
  void updateRandomField(double (*FacePtr)[3]);
  void removeMeanAndSetVarianceOne(double (*FacePtr)[3]);
  void computeCholesky(double (*FacePtr)[3]);
  void updateBC(double (*FacePtrNew)[3],double (*FacePtrOld)[3]);
  void query(const string& param_str) const;
  void updateCvData_(double (*var_)[3]);
  void updateCvData_(double *var_);
  int solveVecBoundaryCvCg(double (*sol)[3],const double * const A,double (*rhs)[3], double * AddedDiag,double zero, int maxiter );
  void getCvData_(double *FacePtr, double *CvPtr_);
  void getCvVectorData_(double (*FacePtr)[3], double (*CvPtr)[3]);

};

class OutletHBc : public HelmholtzBc {
public:

  double * mf0;
  double * mf00;
  double * p_bc;
  double * p0_bc;
  double * p00_bc;
  double sigma;
  double L;
  double p_ref;
  double Ma;

  OutletHBc(BfZone* p, HelmholtzSolver* s): HelmholtzBc(p,s) {


    mf = NULL; zone_ptr->registerBfData(mf,"mf",READWRITE_DATA);

    mf0 = NULL;
    mf00 = NULL;
    if (!b_helmholtzCorrection){
      zone_ptr->registerBfData(mf0,"mf0",READWRITE_DATA);
      zone_ptr->registerBfData(mf00,"mf00",READWRITE_DATA);
    }

    p_bc   = NULL; zone_ptr->registerBfData(p_bc,"p_bc",READWRITE_DATA);
    p0_bc  = NULL; zone_ptr->registerBfData(p0_bc,"p0_bc",READWRITE_DATA);
    p00_bc = NULL; zone_ptr->registerBfData(p00_bc,"p00_bc",READWRITE_DATA);

    sigma = 0.0;
    L     = 0.0;
    p_ref = 0.0;   // this is the outlet pressure - some global reference pressure.
    Ma    = 0.0;

  }
  OutletHBc(BfZone* p, HelmholtzSolver* s,const int icg): HelmholtzBc(p,s,icg) {
 
    mf     = NULL;
    mf0    = NULL;
    mf00   = NULL;
    p_bc   = NULL; 
    p0_bc  = NULL; 
    p00_bc = NULL; 

    sigma = 0.0;
    L     = 0.0;
    p_ref = 0.0;   // this is the outlet pressure - some global reference pressure.
    Ma    = 0.0;
  }

  ~OutletHBc(){
    //mf deleted by parent
    DELETE(mf0);
    DELETE(mf00);
    DELETE(p_bc);
    DELETE(p0_bc);
    DELETE(p00_bc);
  }

  //empty overrides of pure virtual functions
  void syncPostState() {}

  void resetBc();
  void initData();
  void initFromCoarseGrid(StaticSolver::CoarseGrid* cg);

  void initialHook();
  void setBc();
  void addMomentumFlux(double *A, double (*rhs)[3]) const;
  void addMomentumFlux(double *A,double *At,double (*rhs)[3]) const;
  void addPressureFlux(double *A,double *rhs) const;
  void addMassFlux(double *rhs) const;
  void addScalarBoundaryFlux(double *A,double *rhs,const ScalarField *scalar) const {
    _addScalarBoundaryFluxConvective(A,NULL,rhs,scalar);
  }
  void addScalarBoundaryFlux(double *A,double *At,double *rhs,const ScalarField *scalar) const {
    _addScalarBoundaryFluxConvective(A,At,rhs,scalar);
  }
  void completePressureGrad(double (*dpdx)[3]) const;
  void updateMf();
  void correctMf();
  void updateP();

  void force_bf(double (*f_bf)[9]) const;

  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
    const string& name, list<CtiRegister::CtiData>& args,const bool b_eval_func);
};

// VidaVor outlet

class OutletVVHBc : public HelmholtzBc { 
public:

  double * mf0;
  double * mf00;
  double * p_bc;
  double * p0_bc;
  double * p00_bc;
  double sigma;
  double L;
  double p_ref;
  double Ma;
  bool b_localU;

  double * xr;
  double * nr;

  OutletVVHBc(BfZone* p, HelmholtzSolver* s): HelmholtzBc(p,s) {
 
    mf = NULL; zone_ptr->registerBfData(mf,"mf",READWRITE_DATA);

    mf0 = NULL;
    mf00 = NULL;
    if (!b_helmholtzCorrection){
      zone_ptr->registerBfData(mf0,"mf0",READWRITE_DATA);
      zone_ptr->registerBfData(mf00,"mf00",READWRITE_DATA);
    }

    p_bc   = NULL; zone_ptr->registerBfData(p_bc,"p_bc",READWRITE_DATA);
    p0_bc   = NULL; zone_ptr->registerBfData(p0_bc,"p0_bc",READWRITE_DATA);
    p00_bc   = NULL; zone_ptr->registerBfData(p00_bc,"p00_bc",READWRITE_DATA);
      
    sigma = 0.0;
    L     = 0.0;
    p_ref = 0.0;   // this is the outlet pressure - some global reference pressure.
    Ma    = 0.0;

    b_localU = false;

    xr = NULL;
    nr = NULL;

  }
  OutletVVHBc(BfZone* p, HelmholtzSolver* s,const int icg): HelmholtzBc(p,s,icg) {

    mf = NULL;
    mf0 = NULL; 
    mf00 = NULL;

    p_bc   = NULL; 
    p0_bc   = NULL; 
    p00_bc   = NULL; 
      
    sigma = 0.0;
    L     = 0.0;
    p_ref = 0.0; 
    Ma    = 0.0;

    b_localU = false;

    xr = NULL;
    nr = NULL;

  }

  virtual ~OutletVVHBc() {
    //mf deleted by parent
    DELETE(mf0);
    DELETE(mf00);
    DELETE(p_bc);
    DELETE(p0_bc);
    DELETE(p00_bc);
    DELETE(xr);
    DELETE(nr);
  }
  //empty overrides of pure virtual functions
  void syncPostState() {}

  virtual void initData();
  virtual void initFromCoarseGrid(StaticSolver::CoarseGrid* cg);

  void resetBc();
  void setBc();
  void initialHook();
  virtual void addMomentumFlux(double *A,double (*rhs)[3]) const;
  virtual void addMomentumFlux(double *A,double *At,double (*rhs)[3]) const;
  virtual void addPressureFlux(double *A,double *rhs) const;
  void addMassFlux(double *rhs) const;
  void addScalarBoundaryFlux(double *A,double *rhs,const ScalarField *scalar) const {
    _addScalarBoundaryFluxConvective(A,NULL,rhs,scalar);
  }
  void addScalarBoundaryFlux(double *A,double *At,double *rhs,const ScalarField *scalar) const {
    _addScalarBoundaryFluxConvective(A,At,rhs,scalar);
  }
  void completePressureGrad(double (*dpdx)[3]) const;
  void updateMf();
  void correctMf();
  void updateP();

  void force_bf(double (*f_bf)[9]) const;

  void query(const string& param_str) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
    const string& name, list<CtiRegister::CtiData>& args,const bool b_eval_func);
};

class SpongeHBc : public OutletVVHBc {

public:

  double (*u_sponge)[3];
  double sponge_data[5];
  double sponge_speed;
  double sponge_strength;
  SpongeType sponge_type;

  double th_e0[3];
  double th_e1[3];

  SpongeHBc(BfZone* p, HelmholtzSolver* s);
  SpongeHBc(BfZone* p, HelmholtzSolver* s,const int icg);
  ~SpongeHBc() {}

  void initData();
  void initFromCoarseGrid(StaticSolver::CoarseGrid* cg);
  void initSpongeTypeTheta();
  void addMomentumFlux(double *A, double (*rhs)[3]) const;
  void addMomentumFlux(double *A,double *At,double (*rhs)[3]) const;
  void addPressureFlux(double *A,double *rhs) const;

  void addContSource(double *cont) const;
  void computeSigma(double *sigma) const;

};

class SlipWallModelHBc : public HelmholtzBc {

public:

  static bool b_computeSlipLength;

  double * cdel_w;      // slip length for the wm
  double (*rhou_s)[3];  // slip momentum field
  double couple_fax;    // regularization of the nw gradient to suppress 2 delta..

  BoundaryLayerDataExchanger<HelmholtzSolver> * blde;

  SlipWallModelHBc(BfZone* p, HelmholtzSolver* s): HelmholtzBc(p,s) {

    mf = NULL; zone_ptr->registerBfData(mf,"mf",CAN_WRITE_DATA);
    cdel_w = NULL; zone_ptr->registerBfData(cdel_w,"cdel_w",READWRITE_DATA);
    rhou_s = NULL; zone_ptr->registerBfData(rhou_s,"rhou_s",READWRITE_DATA);
    blde = NULL;

    funcEvalList.push_back(p->getName()+":tau_wall()");
    funcEvalList.push_back(p->getName()+":y_plus()");
    funcEvalList.push_back(p->getName()+":bl_delta()");
  }
  SlipWallModelHBc(BfZone* p, HelmholtzSolver* s,const int icg): HelmholtzBc(p,s,icg) {

    mf = NULL;
    cdel_w = NULL; 
    rhou_s = NULL;

    blde = NULL;
  }

  ~SlipWallModelHBc() {
    DELETE(cdel_w);
    DELETE(rhou_s);

    if (!blde) {
      delete blde;
      blde = NULL;
    }

  }

  //empty overrides of pure virtual functions
  void addPressureFlux(double *A,double *rhs) const {}
  void updateMf() {}  
  void updateP() {}
  void correctMf() {}
  void initialHook() {}

  void initData();
  void initFromCoarseGrid(StaticSolver::CoarseGrid* cg);
  void syncPostState();

  void setBc();
  void resetBc();
  void updateBc();
  void restrictBc(StaticSolver::CoarseGrid* cg,HelmholtzBc* bc);

  void addMomentumFlux(double * A, double (*rhs)[3]) const;
  void addMomentumFlux(double * A,double *At,double (*rhs)[3]) const;
  void addMassFlux(double * rhs) const;
  void force_bf(double (*f_bf)[9]) const;

  void query(const string& param_str) const;

 CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
    const string& name, list<CtiRegister::CtiData>& args,const bool b_eval_func);

  void addScalarBoundaryFlux(double *A,double *rhs,const ScalarField *scalar) const {
    _addScalarBoundaryFluxWall(A,NULL,rhs,scalar);
  }
  void addScalarBoundaryFlux(double *A,double *At,double *rhs,const ScalarField *scalar) const {
    _addScalarBoundaryFluxWall(A,At,rhs,scalar);
  }

};


class AlgebraicWallModelHBc : public HelmholtzBc {
public:

  double * tau_wall;
  double z0;
  BoundaryLayerDataExchanger<HelmholtzSolver> * blde;

  AlgebraicWallModelHBc(BfZone* p, HelmholtzSolver* s): HelmholtzBc(p,s), tau_wall(NULL), z0(-1), blde(NULL) {

    zone_ptr->registerBfData(tau_wall,"tau_wall",READWRITE_DATA);

    funcEvalList.push_back(p->getName()+":tau_wall()");
    funcEvalList.push_back(p->getName()+":y_plus()");
    funcEvalList.push_back(p->getName()+":bl_delta()");

  }
  AlgebraicWallModelHBc(BfZone* p, HelmholtzSolver* s,const int icg): HelmholtzBc(p,s,icg), tau_wall(NULL), z0(-1), blde(NULL) {}

  ~AlgebraicWallModelHBc() {

    DELETE(tau_wall);
    if (!blde) {
      delete blde;
      blde = NULL;
    }

  }

  //empty overrides of pure virtual functions
  void resetBc() {}
  void setBc() {}
  void addMassFlux(double * rhs) const {}
  void addPressureFlux(double *A,double *rhs) const {}
  void updateMf() {}
  void correctMf() {}
  void updateP() {}
  void syncPostState() {}

  void initData();
  void initFromCoarseGrid(StaticSolver::CoarseGrid* cg);
  void initialHook();

  void restrictBc(StaticSolver::CoarseGrid* cg,HelmholtzBc* bc);

  void addMomentumFlux(double * A, double (*rhs)[3]) const;
  void addMomentumFlux(double * A,double *At,double (*rhs)[3]) const;
  void force_bf(double (*f_bf)[9]) const;

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
    const string& name, list<CtiRegister::CtiData>& args,const bool b_eval_func);

  void query(const string& param_str) const;

  void addScalarBoundaryFlux(double *A,double *rhs,const ScalarField *scalar) const {
    _addScalarBoundaryFluxWall(A,NULL,rhs,scalar);
  }
  void addScalarBoundaryFlux(double *A,double *At,double *rhs,const ScalarField *scalar) const {
    _addScalarBoundaryFluxWall(A,At,rhs,scalar);
  }

};


#endif
