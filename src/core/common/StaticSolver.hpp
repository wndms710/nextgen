#ifndef _STATIC_SOLVER_HPP_
#define _STATIC_SOLVER_HPP_

#include "KillfileReader.hpp"
#include "StripedMesh.hpp"
#include "MiscUtils.hpp"
#include "CtiRegister.hpp"
#include "WebUI.hpp"

// for imaging...
#include "CtiScene.hpp"
#include "SimpleTri.hpp"

// for surface stuff...
#include "SubSurface.hpp"
#include "SurfaceShm.hpp"

// probes...
#include "MultiFluxProbe.hpp"
#include "FluxProbe.hpp"
#include "PointProbe.hpp"
#include "ConditionalProbe.hpp"
#include "PointCloudProbe.hpp"
#include "VolumetricProbe.hpp"
#include "FazoneProbe.hpp"
#include "PdfProbe.hpp"
#include "FwhSurface.hpp"

// container for writing data...
#include "DataWriter.hpp"

// snapshots...
#include "SnapshotIO.hpp"

// for interp
#include "DataContainer.hpp"

// prcomm stuff...
#include "Prcomm.hpp"

#include "CtiInterpolateTransform.hpp"

// for Lp
#include "LpData.hpp"

#include "Multigrid.hpp"
// for coordinate transform
#include "GeomUtils.hpp"

class BfZone {
private:
  string name;
public:
  DistributedDataExchanger * dde_striped; // pull and push from the striped dist (the one for reading/writing)
  int8 * bfora_striped;
public:
  int index;
  int8 lb_cost;                  // load balance cost -- set this using loadBalanceHook()
  int8 nbf_global,ibf_f_global;  // global count and first index of boundary faces
  int ibf_f,ibf_l,nbf;           // local inclusive range of boundary faces, local count
  double area_global;            // global area (NOT projected area)
  double area_over_delta_global; // sum of all bf's area_over_delta
  double n_global[3];            // sum of all face normals (area magitude)
  double x_global[3];            // area-weighted center of mass
  int8 * ibf_global;
  double * area_bf;
  double * area_over_delta_bf;
  double (*n_bf)[3];
  double (*x_bf)[3];
  int *cvobf;
  BfZone() {
    dde_striped = NULL;
    bfora_striped = NULL;
    index = -1;
    lb_cost = 1; // default cost to apply this bc - gets modifed by setLoad
    nbf_global = 0;
    ibf_f_global = 0;
    ibf_f = 0;
    ibf_l = -1;
    nbf = 0;
    area_global = 0.0;
    area_over_delta_global = 0.0;
    n_global[0] = 0.0;
    n_global[1] = 0.0;
    n_global[2] = 0.0;
    x_global[0] = 0.0;
    x_global[1] = 0.0;
    x_global[2] = 0.0;
    // simple copies of the solver ptrs offset/indexed to our location [0:nbf)...
    ibf_global = NULL;
    area_bf = NULL;
    area_over_delta_bf = NULL;
    n_bf = NULL;
    x_bf = NULL;
    cvobf = NULL;
  }
  ~BfZone() {
    if (dde_striped) delete dde_striped;
    if (bfora_striped) delete[] bfora_striped;
  }

  void setName(const string& name) { this->name = name; }

  string getName() const { return name; }

  // data registration...
  void registerBfIN(const string& vname,const uint rw_bits);
  void registerBfDN(const string& vname,const uint rw_bits);
  void registerBfDN3(const string& vname,const uint rw_bits);
  void registerBfData(int *&val,const string& vname,const uint rw_bits);
  void registerBfData(double *&val,const string& vname,const uint rw_bits);
  void registerBfData(double (*&val)[3],const string& vname,const uint rw_bits);
  void registerBfData(double& val,const string& vname,const uint rw_bits);

  template<class T>
  void registerBfData(T *& state,double &val,const string& vname,const uint rw_bits) {
    const int topo = ((index<<INDEX_SHIFT)|BF_DATA);
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(state,val,name+":"+vname,topo,rw_bits,nbf,bfora_striped,dde_striped);
    }
    else {
      CtiRegister::_registerData(state,val,name+":"+vname,topo,rw_bits,nbf);
    }
  }
  template<class T>
  void registerBfData(T *& state,int &val,const string& vname,const uint rw_bits) {
    const int topo = ((index<<INDEX_SHIFT)|BF_DATA);
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(state,val,name+":"+vname,topo,rw_bits,nbf,bfora_striped,dde_striped);
    }
    else {
      CtiRegister::_registerData(state,val,name+":"+vname,topo,rw_bits,nbf);
    }
  }
  template<class T>
  void registerBfData(T *& state,double (&val)[3],const string& vname,const uint rw_bits) {
    const int topo = ((index<<INDEX_SHIFT)|BF_DATA);
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(state,val,name+":"+vname,topo,rw_bits,nbf,bfora_striped,dde_striped);
    }
    else {
      CtiRegister::_registerData(state,val,name+":"+vname,topo,rw_bits,nbf);
    }
  }
  bool checkDataFlag(const string& vname) const;
  void setDataFlag(const string& vname,const int val);
  void initDdeStuff();

  // variables and data requests to the bc zones are prefixed with
  // "<zone-name>:", and so provide a fast way to check the leading characters
  // of the boundary condition to see they match

  bool isLeadingStrMatched( const string& str) const {

    // the incoming string must be at least the same length as the zone name

    const size_t name_len  = name.length();
    if ( str.length() > name_len) {
      return ( str.compare( 0, name_len, name) == 0);
    }

    return false;
  }


  double* createBfD1Data(CtiRegister::CtiData& v) {

    v.new_dn(BF_DATA,nbf);
    v.setDdeStuff(bfora_striped,dde_striped);
    v.setBits(CAN_WRITE_DATA);
    v.setTopology(BF_DATA,index);
    return v.getDNptr();

  }

  double (* createBfD2Data(CtiRegister::CtiData& v))[3]  {

    v.new_dn3(BF_DATA,nbf);
    v.setDdeStuff(bfora_striped,dde_striped);
    v.setBits(CAN_WRITE_DATA);
    v.setTopology(BF_DATA,index);
    return v.getDN3ptr();

  }

};

// note that these guys are stored in registration space... so the sizes need to be set to allocate!!!
// -----------------------------------------------------------------------------------------------------
inline void BfZone::registerBfIN(const string& vname,const uint rw_bits) {
  // TODO: we don't have a means to initialize the dde in Stats (Probes) when they create
  // the first read/write data for the bfzone -- FH/CI
  const int topo = ((index<<INDEX_SHIFT)|BF_DATA);
  if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
    CtiRegister::_registerData(vname,topo,IN_DATA,rw_bits,nbf,bfora_striped,dde_striped);
  }
  else {
    CtiRegister::_registerData(vname,topo,IN_DATA,rw_bits,nbf);
  }
}
inline void BfZone::registerBfDN(const string& vname,const uint rw_bits) {
  const int topo = ((index<<INDEX_SHIFT)|BF_DATA);
  if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
    CtiRegister::_registerData(vname,topo,DN_DATA,rw_bits,nbf,bfora_striped,dde_striped);
  }
  else {
    CtiRegister::_registerData(vname,topo,DN_DATA,rw_bits,nbf);
  }
}
inline void BfZone::registerBfDN3(const string& vname,const uint rw_bits) {
  const int topo = ((index<<INDEX_SHIFT)|BF_DATA);
  if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
    CtiRegister::_registerData(vname,topo,DN3_DATA,rw_bits,nbf,bfora_striped,dde_striped);
  }
  else {
    CtiRegister::_registerData(vname,topo,DN3_DATA,rw_bits,nbf);
  }
}

inline void BfZone::registerBfData(int *&val,const string& vname,const uint rw_bits) {
  const int topo = ((index<<INDEX_SHIFT)|BF_DATA);
  if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
    CtiRegister::_registerData(val,name+":"+vname,topo,rw_bits,nbf,bfora_striped,dde_striped);
  }
  else {
    CtiRegister::_registerData(val,name+":"+vname,topo,rw_bits,nbf);
  }
}

inline void BfZone::registerBfData(double *&val,const string& vname,const uint rw_bits) {
  const int topo = ((index<<INDEX_SHIFT)|BF_DATA);
  if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
    CtiRegister::_registerData(val,name+":"+vname,topo,rw_bits,nbf,bfora_striped,dde_striped);
  }
  else {
    CtiRegister::_registerData(val,name+":"+vname,topo,rw_bits,nbf);
  }
}

inline void BfZone::registerBfData(double (*&val)[3],const string& vname,const uint rw_bits) {
  const int topo = ((index<<INDEX_SHIFT)|BF_DATA);
  if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
    CtiRegister::_registerData(val,name+":"+vname,topo,rw_bits,nbf,bfora_striped,dde_striped);
  }
  else {
    CtiRegister::_registerData(val,name+":"+vname,topo,rw_bits,nbf);
  }
}

inline bool BfZone::checkDataFlag(const string& vname) const {
  return CtiRegister::checkDataFlag(name+":"+vname);
}

inline void BfZone::setDataFlag(const string& vname,const int val) {
  CtiRegister::setDataFlag(name+":"+vname,val);
}

inline void BfZone::registerBfData(double& val,const string& vname,const uint rw_bits) {
  CtiRegister::_registerData(val, name+":"+vname, rw_bits);
}

inline void BfZone::initDdeStuff() {
  // to support boundary data io, we need to build a dde for any boundary data that
  // gets registered.
  assert(dde_striped == NULL);
  assert(bfora_striped == NULL);
  assert(ibf_global != NULL);
  int8 * my_ibf_global = new int8[nbf];
  for (int ibf = 0; ibf < nbf; ++ibf) {
    assert((ibf_global[ibf] >= ibf_f_global)&&(ibf_global[ibf] < ibf_f_global+nbf_global));
    my_ibf_global[ibf] = ibf_global[ibf] - ibf_f_global;
  }
  // everyone build bfora based on our common global index...
  MiscUtils::calcThresholdDist(bfora_striped,nbf_global,mpi_size,DIST_THRESHOLD);
  dde_striped = new DistributedDataExchanger(my_ibf_global,nbf,bfora_striped);
  delete[] my_ibf_global;
}
// -----------------------------------------------------------------------------------------------------

class StaticSolver : public CtiRegister::CtiDataProducer {

private:

  bool b_skip_image;
  string skip_image_name;

  // the solver can know about its logfile (where its
  // cout is being collected) using the SET_LOG param...
  bool b_logfile;
  string logfile;

public:

  int8 ncv_global;
  int8 nbf_global;
  int8 nfa_global;
  int8 nno_global;
  int8 nef_global;
  int8 nno_pb_global;

  int8 * icv_global; // always positive
  int8 * ibf_global;
  int8 * ifa_global; // uses -1 indexing for tie-breaking
  int8 * ino_global; // always positive for now
  int cg_level;

protected:

  double reset_stats_time;
  int reset_stats_int;  // interval if specified

  DistributedDataExchanger * dde_cv_striped; // pull and push from the cv striped dist (the one for reading/writing)
  int8 * cvora_striped;
  uint8 * rbi_sm; // this tells us the rbi's our sm needs to send to (persistent for lp)

  DistributedDataExchanger * dde_fa_striped; // pull and push from the fa striped dist (the one for reading/writing)
  int8 * faora_striped;

  // periodicity data and functions...

  //vector<PeriodicTransform> periodicTransformVec; // should this be in a namespace?

  /*
    void periodicTranslate(double (*xp_t)[3],const int n,const int bits) const {
    for (int bit_pair = 0; bit_pair < 3; ++bit_pair) {
    if (bits & (1<<(2*bit_pair))) {
    assert((bits & (1<<(2*bit_pair+1)))==0);
    assert(periodicTransformVec.size() > bit_pair);
    periodicTransformVec[bit_pair].translate(xp_t,n);
    }
    else if (bits & (1<<(2*bit_pair+1))) {
    assert(periodicTransformVec.size() > bit_pair);
    periodicTransformVec[bit_pair].inv_translate(xp_t,n);
    }
    }
    // these ones are the doubles...
    for (int bit_pair = 3; bit_pair < 6; ++bit_pair) {
    if (bits & (1<<(2*bit_pair))) {
    assert((bits & (1<<(2*bit_pair+1)))==0);
    assert(periodicTransformVec.size() > bit_pair-3);
    periodicTransformVec[bit_pair-3].translate(xp_t,n);
    periodicTransformVec[bit_pair-3].translate(xp_t,n);
    }
    else if (bits & (1<<(2*bit_pair+1))) {
    assert(periodicTransformVec.size() > bit_pair-3);
    periodicTransformVec[bit_pair-3].inv_translate(xp_t,n);
    periodicTransformVec[bit_pair-3].inv_translate(xp_t,n);
    }
    }

    }

    void periodicTranslate(double *xp_t,const int n,const int bits) const {
    periodicTranslate((double(*)[3])xp_t,n,bits);
    }

    void periodicRotate(double (*xp_t)[3],const int n,const int bits) const {
    for (int bit_pair = 0; bit_pair < 3; ++bit_pair) {
    if (bits & (1<<(2*bit_pair))) {
    assert((bits & (1<<(2*bit_pair+1)))==0);
    assert(periodicTransformVec.size() > bit_pair);
    periodicTransformVec[bit_pair].rotate(xp_t,n);
    }
    else if (bits & (1<<(2*bit_pair+1))) {
    assert(periodicTransformVec.size() > bit_pair);
    periodicTransformVec[bit_pair].inv_rotate(xp_t,n);
    }
    }
    for (int bit_pair = 3; bit_pair < 6; ++bit_pair) {
    if (bits & (1<<(2*bit_pair))) {
    assert((bits & (1<<(2*bit_pair+1)))==0);
    assert(periodicTransformVec.size() > bit_pair-3);
    periodicTransformVec[bit_pair-3].rotate(xp_t,n);
    periodicTransformVec[bit_pair-3].rotate(xp_t,n);
    }
    else if (bits & (1<<(2*bit_pair+1))) {
    assert(periodicTransformVec.size() > bit_pair-3);
    periodicTransformVec[bit_pair-3].inv_rotate(xp_t,n);
    periodicTransformVec[bit_pair-3].inv_rotate(xp_t,n);
    }
    }
    }

    void periodicRotate(double *xp_t,const int n,const int bits) const {
    periodicRotate((double(*)[3])xp_t,n,bits);
    }
  */

  //bool getPeriodicR(double R[9],const int bits) const;

  //bool getPeriodicT(double t[3],const int bits) const;

protected:

  KillfileReader kfr;

  // costs for load-balancer. These can be adjusted from their defaults
  // by implementing the virtual routine setLoadBalanceCost()

  int8 cv_lb_cost;
  int8 fa_lb_cost;
  int8 ef_lb_cost;

public:

  Adt<double> * cvBboxAdt;
  Adt<double> * cvAdt;

  set<string> added_data_fields; // e.g.: particles,interfaces,tbd

  vector<BfZone> bfZoneVec;
  map<const string,int> bfZoneNameMap;

  int ncv; // active cv count
  int ncv_g; // active + first-layer cv ghosts
  int ncv_g2; // active + first and second-layer ghosts
  int ncv_i; // active and all cv nbrs (first layer) are local
  double (*x_vv)[3];
  double * r_vv;
  double (*x_cv)[3];
  double * vol_cv;
  double * inv_vol;
  bool * b_open_faces_cv;

  int nbf;
  int * zone_bf;
  int * cvobf;
  double (*n_bf)[3];
  double (*x_bf)[3];
  double * area_bf;
  double * area_over_delta_bf;
  double (*Gij_bf)[3][3];
  int * noobf_i;
  int * noobf_v;
  // link between (sub?) surface tris and boundary faces...
  int * sstobf_i;
  int * sstobf_v; // ist | (bits<<52)
  int * sspobf_i;
  int * sspobf_v; // ist | (bits<<52)
  double * sspobf_wgt;

  // surface stuff if present...
  SubSurface * subSurface;

  int nfa_i; // internal face count
  int nfa; // internal + interprocessor faces -- i.e. all faces
  int * group_fa;
  int (*cvofa)[2];
  double *area_fa;
  double (*n_fa)[3];
  double (*x_fa)[3];
  double *area_over_delta_fa;
  int * noofa_i;
  int * noofa_v;

  int nno_b; // first nodes are associated with bf's
  int nno;
  double (*x_no)[3];

  int nef_i;
  int nef;
  int * group_ef;
  int (*cvoef)[2];
  double (*n_ef)[3];
  double (*c_ef)[3];

  // additional operators...
  // cv-nbr-of-cv CSR structure: used for compact gradient and compact implicit solvers...

  int * faocv_i;
  int * faocv_v;

  int * bfocv_i;
  int * bfocv_v;

  int * cvocv_i;
  int * cvocv_v;
  int * cvocv2_i;
  int * cvocv2_v;

  double (*cvocv_grad_coeff)[3]; // gradient for models (same as gradient used to build extended faces)

  double * filter_mom2;

  // -----------------------------------------------
  // the cv-based communicators and ghost data
  // reduction machinery -- now you are down in the
  // guts of the parallel solver ;)
  // -----------------------------------------------

  uint8 * rbi_g;
  uint8 * rbi_g2;

  int (*iRt_g)[2]; // [Rotation][translation] index for each ghost (-1 if DNE for either/both)
  double (*per_R)[9]; // periodic Rotation matrices
  double (*per_t)[3]; // periodic translation vecs

  // stuff for building dual...
  int ncv_d; // active cvs + nbrs within delta part of tets on this rank
  uint8* rbi_d;
  double (*x_cv_d)[3];
  vector<CvPrcomm> cvdPrcommVec; // for ghost nbrs for dual
  class DualTet {
    public:
      int icv[4];
      DualTet() {
        this->icv[0] = 0.0;
        this->icv[1] = 0.0;
        this->icv[2] = 0.0;
        this->icv[3] = 0.0;
      }
      DualTet(const DualTet& other) {
        this->icv[0] = other.icv[0];
        this->icv[1] = other.icv[1];
        this->icv[2] = other.icv[2];
        this->icv[3] = other.icv[3];
      }
      DualTet(const int icv[4]) {
        this->icv[0] = icv[0];
        this->icv[1] = icv[1];
        this->icv[2] = icv[2];
        this->icv[3] = icv[3];
      }
      DualTet(const int icv0,const int icv1,const int icv2,const int icv3) {
        this->icv[0] = icv0;
        this->icv[1] = icv1;
        this->icv[2] = icv2;
        this->icv[3] = icv3;
      }
      bool operator<(const DualTet& other) const {
        return (icv[0] < other.icv[0]) || ((icv[0] == other.icv[0]) &&
              ((icv[1] < other.icv[1]) || ((icv[1] == other.icv[1]) &&
              ((icv[2] < other.icv[2]) || ((icv[2] == other.icv[2]) &&
               (icv[3] < other.icv[3]))))));
      }
      bool operator!=(const DualTet& other) const {
        return (icv[0] != other.icv[0]) || (icv[1] != other.icv[1]) || (icv[2] != other.icv[2]) || (icv[3] != other.icv[3]);
      }
      bool operator==(const DualTet& other) const {
        return (icv[0] == other.icv[0]) && (icv[1] == other.icv[1]) && (icv[2] == other.icv[2]) && (icv[3] == other.icv[3]);
      }
  };
  vector<DualTet> dualTetVec;

public: // made public for namespace copy in StaticSolverNS
  vector<CvPrcomm> cvPrcommVec; // for first layer of ghosts
protected:
  vector<CvPrcomm> cv2PrcommVec; // for second layer of ghosts

  // mpi requests map the void data pointer to a set of
  // non-blocking requests that have yet to be completed...
  // this map relates update*DataStart() and update*DataFinish() non-blocking
  // ghost update routines...
  map<const void*,MpiRequestStuff*> mpiRequestMap;

  // a node-based communicator...

  vector<Prcomm> noPrcommVec;
  vector<Prcomm> faPrcommVec;

  vector<CtiInterpolateTransform*> interp_transforms;

private:

  // geometry...

  double boundingBox[6];
  bool b_boundingBox;

  // coordinate transform...
  double tBoundingBox[6];
  bool b_tBoundingBox;
  int transform_type;  // 0: cartesian; 1: transform_cyl_x; 2: transform_cyl_z
  double r0,cos_theta0,sin_theta0;
  double cyl_center[3];

  // probing...

  map<const string,MultiFluxProbe> mfpMap; // parameter-based multi-flux probe map: why make this a map?
  map<const string,FluxProbe> fpMap;
  map<const string,PointProbe> ppMap;
  map<const string,VolumetricProbe> vpMap;
  map<const string,VolumetricProbe> svpMap;
  map<const string,FazoneProbe> fzpMap;
  map<const string,ConditionalProbe> cpMap;
  map<const string,ConditionalProbe> scpMap;
  map<const string,PointCloudProbe*> pcpMap;
  map<const string,PdfProbe> pdfpMap;
  map<const string,FwhSurface> fwhSurfaceMap;

  // snapshots..
  map<const string,Snapshot> snapMap;

  // tecplot, vtk, ptk ...
  map<const string,DataWriter*> dwMap;

  // images...
  list<CtiScene*> sceneList; // list of buffered images to be flushed

  // write_sles...
  string write_sles_prefix;

protected:

  // make private after adopting WRITE_SLES behavior in FlowSolver...
  int write_sles_interval; // -1 means not parsed, 0 means parsed and not present, 1,2,3... means parsed and prefix set

  // forces..
  map<string,int> forceZoneMap;
  double (*forceZoneBuf)[9];
  double (*momentZoneBuf)[9];

public:

  StaticSolver() {

    // TODO: eventually get rid of lpHelperVec and manage particle io
    // more fundamentally as general data in CtiRegister. For now, we
    // need to prevent this vector from resizing because it wrecks the
    // pointer registration that has been set up...
    lpHelperVec.reserve(10);

    b_skip_image = false;

    reset_stats_time = HUGE_VAL;
    reset_stats_int = -1;

    b_logfile = false;

    cg_level = -1;

    //
    // The following registration allows StaticSolver objects (and any inherited solvers) to implement custom
    // function evaluations to provide data to CtiRegister's evaluation machinery...
    //
    // virtual CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,const string& name,list<CtiRegister::CtiData>& args);
    //

    CtiRegister::addCtiDataProducer(this);

    // ...

    icv_global = NULL;
    ibf_global = NULL;
    ifa_global = NULL;
    ino_global = NULL;

    dde_cv_striped = NULL;
    cvora_striped = NULL;
    rbi_sm = NULL;

    dde_fa_striped = NULL;
    faora_striped = NULL;

    cvBboxAdt = NULL;
    cvAdt = NULL;

    // these values along with the lb_cost in each bfZone get set
    // in the virtual routine setLoadBalanceCost()...

    cv_lb_cost = 10;
    fa_lb_cost = 1;
    ef_lb_cost = 1;

    ncv = 0;
    ncv_i = 0;
    ncv_g = 0;
    ncv_g2 = 0;
    ncv_global = 0;
    x_vv = NULL;
    r_vv = NULL;
    x_cv = NULL;
    vol_cv = NULL;
    inv_vol = NULL;
    b_open_faces_cv = NULL;

    nbf = 0;
    nbf_global = 0;
    zone_bf = NULL;
    cvobf = NULL;
    n_bf = NULL;
    x_bf = NULL;
    area_bf = NULL;
    area_over_delta_bf = NULL;
    Gij_bf = NULL;
    noobf_i = NULL;
    noobf_v = NULL;
    sstobf_i = NULL;
    sstobf_v = NULL;
    sspobf_i = NULL;
    sspobf_v = NULL;
    sspobf_wgt = NULL;

    subSurface = NULL;

    nfa_global = 0;
    nfa_i = nfa = 0;
    group_fa = NULL;
    cvofa = NULL;
    area_fa = NULL;
    n_fa = NULL;
    x_fa = NULL;
    area_over_delta_fa = NULL;
    noofa_i = NULL;
    noofa_v = NULL;

    nno_b = nno = 0;
    nno_global = 0;
    x_no = NULL;

    nef_global = 0;
    nef_i = nef = 0;
    group_ef = NULL;
    cvoef = NULL;
    n_ef = NULL;
    c_ef = NULL;

    rbi_g = NULL;
    rbi_g2 = NULL;

    iRt_g = NULL;
    per_R = NULL;
    per_t = NULL;

    faocv_i = NULL;
    faocv_v = NULL;

    bfocv_i = NULL;
    bfocv_v = NULL;

    cvocv_i = NULL;
    cvocv_v = NULL;
    cvocv2_i = NULL;
    cvocv2_v = NULL;

    cvocv_grad_coeff = NULL;
    filter_mom2      = NULL;

    ncv_d = 0;
    rbi_d = NULL;
    x_cv_d = NULL;

    b_boundingBox = false;
    b_tBoundingBox = false;

    transform_type = -1;
    r0 = HUGE_VAL;
    cos_theta0 = HUGE_VAL;
    sin_theta0 = HUGE_VAL;
    FOR_I3 cyl_center[i] = HUGE_VAL;

    write_sles_interval = -1; // -1 means not parsed, 0 means parsed and not present, 1,2,3... means parsed and prefix set
    write_sles_prefix = "result"; // default prefix

    // data registration now occurs in constructors

    registerData();

  }

  StaticSolver(const int icg) {

    b_skip_image = false;

    reset_stats_time = HUGE_VAL;
    reset_stats_int = -1;

    b_logfile = false;

    cg_level = icg;

    // ...

    icv_global = NULL;
    ibf_global = NULL;
    ifa_global = NULL;
    ino_global = NULL;

    dde_cv_striped = NULL;
    cvora_striped = NULL;
    rbi_sm = NULL;

    dde_fa_striped = NULL;
    faora_striped = NULL;

    cvBboxAdt = NULL;
    cvAdt = NULL;

    // these values along with the lb_cost in each bfZone get set
    // in the virtual routine setLoadBalanceCost()...

    cv_lb_cost = 10;
    fa_lb_cost = 1;
    ef_lb_cost = 1;

    ncv = 0;
    ncv_i = 0;
    ncv_g = 0;
    ncv_global = 0;
    x_vv = NULL;
    r_vv = NULL;
    x_cv = NULL;
    vol_cv = NULL;
    inv_vol = NULL;
    b_open_faces_cv = NULL;

    nbf = 0;
    nbf_global = 0;
    zone_bf = NULL;
    cvobf = NULL;
    n_bf = NULL;
    x_bf = NULL;
    area_bf = NULL;
    area_over_delta_bf = NULL;
    Gij_bf = NULL;
    noobf_i = NULL;
    noobf_v = NULL;
    sstobf_i = NULL;
    sstobf_v = NULL;
    sspobf_i = NULL;
    sspobf_v = NULL;
    sspobf_wgt = NULL;

    subSurface = NULL;

    nfa_global = 0;
    nfa_i = nfa = 0;
    group_fa = NULL;
    cvofa = NULL;
    area_fa = NULL;
    n_fa = NULL;
    x_fa = NULL;
    area_over_delta_fa = NULL;
    noofa_i = NULL;
    noofa_v = NULL;

    nno_b = nno = 0;
    nno_global = 0;
    x_no = NULL;

    nef_global = 0;
    nef_i = nef = 0;
    group_ef = NULL;
    cvoef = NULL;
    n_ef = NULL;
    c_ef = NULL;

    rbi_g = NULL;
    rbi_g2 = NULL;

    iRt_g = NULL;
    per_R = NULL;
    per_t = NULL;

    faocv_i = NULL;
    faocv_v = NULL;

    bfocv_i = NULL;
    bfocv_v = NULL;

    cvocv_i = NULL;
    cvocv_v = NULL;
    cvocv2_i = NULL;
    cvocv2_v = NULL;

    cvocv_grad_coeff = NULL;
    filter_mom2      = NULL;

    ncv_d = 0;
    rbi_d = NULL;
    x_cv_d = NULL;

    b_boundingBox = false;

    write_sles_interval = -1;

  }

  virtual ~StaticSolver() {

    DELETE(icv_global);
    DELETE(ibf_global);
    DELETE(ifa_global);
    DELETE(ino_global);

    if (dde_cv_striped != NULL) delete dde_cv_striped;
    DELETE(cvora_striped);
    DELETE(rbi_sm);

    if (dde_fa_striped != NULL) delete dde_fa_striped;
    DELETE(faora_striped);

    if (cvBboxAdt) delete cvBboxAdt;
    if (cvAdt) delete cvAdt;

    DELETE(x_vv);
    DELETE(r_vv);
    DELETE(x_cv);
    DELETE(vol_cv);
    DELETE(inv_vol);
    DELETE(b_open_faces_cv);

    DELETE(zone_bf);
    DELETE(cvobf);
    DELETE(n_bf);
    DELETE(x_bf);
    DELETE(area_bf);
    DELETE(area_over_delta_bf);
    DELETE(Gij_bf);
    DELETE(noobf_i);
    DELETE(noobf_v);
    DELETE(sstobf_i);
    DELETE(sstobf_v);
    DELETE(sspobf_i);
    DELETE(sspobf_v);
    DELETE(sspobf_wgt);

    if (subSurface) delete subSurface;

    DELETE(group_fa);
    DELETE(cvofa);
    DELETE(area_fa);
    DELETE(n_fa);
    DELETE(x_fa);
    DELETE(area_over_delta_fa);
    DELETE(noofa_i);
    DELETE(noofa_v);

    DELETE(x_no);

    DELETE(group_ef);
    DELETE(cvoef);
    DELETE(n_ef);
    DELETE(c_ef);

    DELETE(rbi_g);
    DELETE(rbi_g2);
    DELETE(iRt_g);
    DELETE(per_R);
    DELETE(per_t);

    DELETE(faocv_i);
    DELETE(faocv_v);

    DELETE(bfocv_i);
    DELETE(bfocv_v);

    DELETE(cvocv_i);
    DELETE(cvocv_v);
    DELETE(cvocv2_i);
    DELETE(cvocv2_v);

    DELETE(cvocv_grad_coeff);
    DELETE(filter_mom2);

    DELETE(rbi_d);
    DELETE(x_cv_d);

    // for now we keep the point cloud and data writer pointers for the entire run...
    // we should do some checking in the future
    for (map<const string,PointCloudProbe*>::iterator iter = pcpMap.begin(); iter != pcpMap.end(); ++iter) {
      delete iter->second;
    }
    for (map<const string,DataWriter*>::iterator iter = dwMap.begin(); iter != dwMap.end(); ++iter) {
      delete iter->second;
    }

    for (int ii = 0, lim = lpHelperVec.size(); ii < lim; ++ii)
      lpHelperVec[ii].clear();


    for(vector<CtiInterpolateTransform*>::iterator it = interp_transforms.begin();
        it != interp_transforms.end(); ++it) {

      delete *it; *it = NULL;

    }
  }

  enum InitBits {
    INIT_COMPACT_FACES  = 1,
    INIT_EXTENDED_FACES = 2,
    INIT_CV_GRAD        = 4,
    INIT_ALL            = 65535, // all bits
  };

  void initBfDdeStuff() {
    FOR_IZONE(bfZoneVec) bfZoneVec[izone].initDdeStuff();
  }

  virtual void resizeData() {

    // allocates any data that hasn't been yet...

    CtiRegister::_initData();

  }

  void readData(const string& filename);

  void readLpocvAndInitDdeStuff(const string filename);

  void redistReorderData(StripedMesh* sm);

  void flipRWSignedFaData();

  void registerFromParams() {

    // data that is registered from the param list will be treated as READWRITE_DATA
    // in the future, you could templatize the following chunks of code to eliminate redundancy.

    // error management for duplicated registration names is currently checked in CtiRegister
    // the current behavior is to throw errors if duplicated data is found; if this behavior
    // is modified, that should happen in CtiRegister

    FOR_ALL_PARAM {

      if ( ( (param->name == "REGISTER_I0") || (param->name == "REGISTER_I") ) && param->touch() ) {

        for (int iarg = 0; iarg < param->size(); ++iarg)
          registerI( param->getString(iarg), READWRITE_DATA);

      } else if ( ( (param->name == "REGISTER_R0") || (param->name == "REGISTER_D") ) && param->touch() ) {

        for (int iarg = 0; iarg < param->size(); ++iarg)
          registerD( param->getString(iarg), READWRITE_DATA);

      } else if ( ( (param->name == "REGISTER_CV_R1") || (param->name == "REGISTER_CV_DN") ) && param->touch() ) {

        for (int iarg = 0; iarg < param->size(); ++iarg)
          registerCvDN( param->getString(iarg), READWRITE_DATA);

      } else if ( ( (param->name == "REGISTER_CV_R2") || (param->name == "REGISTER_CV_DN3") ) && param->touch() ) {

        for (int iarg = 0; iarg < param->size(); ++iarg)
          registerCvDN3( param->getString(iarg), READWRITE_DATA);

      }
      else if ( ( (param->name == "REGISTER_BF_R1") || (param->name == "REGISTER_BF_DN") ) && param->touch() ) {

        for (int iarg = 0; iarg < param->size(); ++iarg) {

          const string name = param->getString(iarg);

          size_t colon_pos = name.find(":");
          if ( colon_pos != string::npos) {

            const string zone_name =  name.substr(0,colon_pos);
            const string var_name  = name.substr(colon_pos+1,name.size()-colon_pos);

            if ( zone_name == "*") {

              // if the wildcard is found, then we will register this data on all zones

              FOR_IZONE(bfZoneVec) {
                const string zone_var_name = bfZoneVec[izone].getName() + ":" + var_name;
                bfZoneVec[izone].registerBfDN(zone_var_name, READWRITE_DATA);
              }

            } else {

              map<const string,int>::const_iterator iter = bfZoneNameMap.find(zone_name);
              if (iter == bfZoneNameMap.end()) {

                CWARN(" > REGISTER_BF_R1:: unable to find zone " << zone_name << ". Skipping...");

              } else {

                const string zone_var_name = bfZoneVec[iter->second].getName() + ":" + var_name;
                bfZoneVec[iter->second].registerBfDN(zone_var_name, READWRITE_DATA);

              }
            }
          } else {

            CWARN(" > REGISTER_BF_R1:: improper zone-indexed variable. Skipping...");

          }

        }
      } else if ( ( (param->name == "REGISTER_BF_R2") || (param->name == "REGISTER_BF_DN3") ) && param->touch() ) {

        for (int iarg = 0; iarg < param->size(); ++iarg) {

          const string name = param->getString(iarg);

          size_t colon_pos = name.find(":");
          if ( colon_pos != string::npos) {

            const string zone_name =  name.substr(0,colon_pos);
            const string var_name  = name.substr(colon_pos+1,name.size()-colon_pos);

            if ( zone_name == "*") {

              // if the wildcard is found, then we will register this data on all zones

              FOR_IZONE(bfZoneVec) {
                const string zone_var_name = bfZoneVec[izone].getName() + ":" + var_name;
                bfZoneVec[izone].registerBfDN3(zone_var_name, READWRITE_DATA);
              }

            } else {

              map<const string,int>::const_iterator iter = bfZoneNameMap.find(zone_name);
              if (iter == bfZoneNameMap.end()) {

                CWARN(" > REGISTER_BF_R2:: unable to find zone " << zone_name << ". Skipping...");

              } else {

                const string zone_var_name = bfZoneVec[iter->second].getName() + ":" + var_name;
                bfZoneVec[iter->second].registerBfDN3(zone_var_name, READWRITE_DATA);

              }
            }
          } else {

            CWARN(" > REGISTER_BF_R2:: improper zone-indexed variable. Skipping...");

          }

        }
      }
      else if (!checkParam("SNAPSHOT")) {

        // particle SNAPSHOT registration is handled in BasicPostpro...

        if ( ( (param->name == "REGISTER_LP_R1") || (param->name == "REGISTER_LP_DN") ) && param->touch() ) {

          // needed by app to control contex menu
          //added_data_fields.insert("particles");

          for (int iarg = 0; iarg < param->size(); ++iarg)
            registerLpDN( param->getString(iarg) );

        } else if ( ( (param->name == "REGISTER_LP_R2") || (param->name == "REGISTER_LP_DN3") ) && param->touch() ) {

          // needed by app to control contex menu
          //added_data_fields.insert("particles");

          for (int iarg = 0; iarg < param->size(); ++iarg)
            registerLpDN3( param->getString(iarg) );

        }

      }

    }

  }

  // these are solver specific...

  virtual void registerBoundaryConditions();

  // didn't make this pure because we could have a solver that doesn't need bc's
  virtual void initBoundaryConditions() {}

  virtual void initData() = 0;

  BfZone * getBfZone(const string& name) {
    map<const string,int>::const_iterator iter = bfZoneNameMap.find(name);
    if (iter == bfZoneNameMap.end())
      return NULL;
    return &bfZoneVec[iter->second];
  }

  void init(const uint init_bits = INIT_ALL) {

    Param * restart_param = getParam("RESTART");
    StripedMesh *sm = NULL;
    string sles_filename = "";
    if ((restart_param == NULL)||(restart_param->size() < 1)) {

      if ( Param * param = getParam("SIMPLE_RESTART")) {

        // if the filename fails to open-- we'll check if this is a simple restart
        // for us to populate a striped mesh, otherwise, this will result in failure.

        sm = new StripedMesh();
        sm->init_simple_restart(param, mpi_comm);

      }
      else {

        CERR("missing RESTART param.\n\nUsage: on the first run:\n\nRESTART restart.mles\n\nor, to restart with data:\n\nRESTART restart.mles result.sles\n");

      }

    }
    else {

      string mles_filename = restart_param->getString(0);
      int ncopy[3] = {1,1,1};
      if (restart_param->size() >= 2) {
        int iarg = 1;
        while (iarg < restart_param->size()) {
          const string token = restart_param->getString(iarg++);
          if ((token == "PATTERN")||(token == "NCOPY")) {
            // the next 3 ints specify the number of unit meshes in each periodic direction. Sign indicates forward/backward copying.
            FOR_I3 {
              ncopy[i] = restart_param->getInt(iarg++);
              // TODO in the future we could allow for ncopy to be signed. i have put in some of the logic for this. but some of the bit logic
              // will have to be expanded.
              if (ncopy[i] <= 0) {
                CERR( " > unable to process PATTERN command; periodic count in direction " << i << " is <= 0.");
              }
            }
          }
          else {
            // assume it is the sles filename...
            sles_filename = token;
          }
        }
      }

      // read mles...

      sm = new StripedMesh();
      sm->read_mles(mles_filename,mpi_comm);

      // perform ncopy operation to build full mesh from periodic component...

      if ((abs(ncopy[0]) > 1)||(abs(ncopy[1]) > 1)||(abs(ncopy[2]) > 1))
        sm->copy_mles(ncopy);

    }

    if (init_bits & INIT_EXTENDED_FACES) sm->buildExtendedFaces();

    // build the bfZoneVec from sm...

    initBfZoneVec(sm);

    // register solver specific bc stuff...

    registerBoundaryConditions();

    // register any user data (needs to come before sles read)...

    registerFromParams();

    // statistics...

    if (Param * param = getParam("STATS")) {
      registerStats(param,false); // false == skip init
    }

    // need to read in registered data into striped, so that solvers
    // can use their data (e.g., vof) to impact the load balance...
    if (sles_filename != "") {

      int ierr = sm->read_sles(sles_filename);
      if ( ierr != 0) {
        CERR( " > unable to read sles with name : " << sles_filename);
      }

      if (sm->lpHelperVec.size() > 0)
        sm->read_sles_lp(restart_param->getString(1));

    }

    // at this point, the solver can do some stuff, like set the load-balance,
    // and initialize their bcs if desired, although no data will be available...

    loadBalanceHook(sm);

    // now we load balance the problem...

    loadBalance(sm);

    // at this point, the sm->x_vv is redistributed and ncv and icv_global are set
    // AND ordered on the StaticSolver's load balanced partition. Now complete
    // the redist and reorder, build communicators, etc...

    redistReorderMesh(sm);

    // print the surface bounding box for user...

    if (!b_boundingBox) calcBoundingBox();

    // init any other required ops -- could control which ones get done
    // from the solver level -- later.

    buildFaocv();

    buildBfocv();

    buildCvocv();

    if (init_bits & INIT_CV_GRAD) calcCvGradCoeff();

    buildNoPrcomm(); // TODO currently does not consider periodicity, but should eventually

    buildFaPrcomm();

    // grouping...

    groupFaces();

    // allocate data solver specific data and set constants

    initData();

    // init solver bc stuff...

    initBfDdeStuff();

    resizeCvDataWithGhosts(icv_global);

    initBoundaryConditions();

    // allocates any data that hasn't been yet...
    CtiRegister::_initData();

    // and read/redistribute if there is a snapshot file to read...

    CtiRegister::clearAllDataFlags();
    if ((restart_param != NULL)&&(restart_param->size() >= 2)) {
      redistReorderData(sm);
    }

    // cleanup StripedMesh...

    delete sm; sm = NULL; // could be moved inside redistReoder

    // interpolate from some data set ...

    interpFromData();

    // take a look at the registered data...
    CtiRegister::dumpRegisteredData();

    // before we potentially reset stats, let's see if we should move any stats
    // vars to their instantaneous counterparts
    initFromParams();

    // now data is set -- check if the user wants to reset the stats...
    if (Param * param = getParam("RESET_STATS")) {
      processResetStats(param,false);  // no help message
    }

    // for now, check for the param "INTERACTIVE" -- this will stop the
    // solver in hold mode the first time processStep gets called...

    if (checkParam("INTERACTIVE"))
      kfr.setHold(true);

    // just to see if any init routines have put something in here...

    CtiRegister::dumpCurrentData();

    // support testing of eval...

    if (checkParam("TEST_EVAL")) {
      while(1) {
        try {
          char c_expression[100];
          int size_expression = 0;
          if (mpi_rank == 0) {
            cout << "\n\nEnter expression to evaluate (exit ends loop): " << endl;
            string expression0;
            cin >> expression0;
            cout << "Evaluating expression \"" << expression0 << "\"..." << endl;
            assert(expression0.size()+1 < 100);
            strcpy(c_expression,expression0.c_str());
            size_expression = expression0.size()+1;
          }
          //         MPI_Barrier(mpi_comm); // need this???
          MPI_Bcast(&size_expression,1,MPI_INT,0,mpi_comm);
          MPI_Bcast(c_expression,size_expression,MPI_CHAR,0,mpi_comm);
          string expression(c_expression);
          if (expression == "exit")
            break;
          CtiRegister::CtiData * data = CtiRegister::getCtiData(expression);
          if (mpi_rank == 0) {
            if (data != NULL) {
              cout << "Got data: " << *data << "\n" << endl;
              CtiRegister::dumpCurrentData();
            }
            else {
              cout << "CtiRegister::getCtiData returned NULL" << endl;
            }
          }
        }
        catch(int ierr) {
          if (mpi_rank == 0) cout << "caught error: " << ierr << endl;
        }
        catch(...) {
          if (mpi_rank == 0) cout << "caught unknown error" << endl;
        }
      }
    }

    // clear the webUIOutput before starting any client interaction.
    // this will prevent the deluge of INFO/WARNING messages that can occur
    // when the session starts...
    WebUI::webUIOutput.clear();

  }

  void initFromParams() {

    if (checkParam("INIT")) COUT1("processing INIT commands:");

    FOR_PARAM_MATCHING("INIT") {

      // INIT <var-name> <value> [<values....>]
      // INIT rho 1.0 u 0.0 0.0 0.0 p 1.0

      int iarg = 0;
      while (iarg < param->size()) {

        const string var_name = param->getString(iarg++);
        CtiRegister::CtiData* data = CtiRegister::getRegisteredCtiData(var_name);
        if ( data != NULL ) {

          const int type = data->getType();

          switch ( type ) {
          case I_DATA:
            {
              const int var_val = param->getInt(iarg++);
              data->i() = var_val;
              COUT1(" > initializing " << var_name << " to " << var_val);
            }
            break;

          case D_DATA:
            {
              const double var_val = param->getDouble(iarg++);
              data->d() = var_val;
              COUT1(" > initializing " << var_name << " to " << var_val);
            }
            break;

          case IN_DATA:
            {
              const int n = data->size();
              const int var_val = param->getInt(iarg++);

              for (int i = 0; i < n; ++i)
                data->in(i) = var_val;

              COUT1(" > initializing " << var_name << " to " << var_val);
            }

            break;

          case DN_DATA:

            {
              const int n = data->size();
              const string from_var = param->getString(iarg++);
              CtiRegister::CtiData * data_from = getCtiData(from_var);  // will perform expression parsing too

              if (data_from != NULL) {
                switch (data_from->getType()) {
                  case DN_DATA:
                    {
                      assert(data_from->size() == n);

                      //TODO: is there a way to determine whether data_from is
                      // a) stats based,
                      // b) whether these stats are populated with real data or just empty (b/c starting now)
                      for (int i =0; i < n; ++i)
                        data->dn(i) = data_from->dn(i);

                      COUT1(" > initializing " << var_name << " to " << from_var);
                    }
                    break;
                  case D_DATA:
                    {
                      for (int i =0; i < n; ++i)
                        data->dn(i) = data_from->d();

                      COUT1(" > initializing " << var_name << " to " << from_var);
                    }
                    break;
                  default:
                    CWARN("cannot set " << data->getTypeAsString() << " from type " << data_from->getTypeAsString() << "; skipping");
                }
              }
            }

            break;

          case DN3_DATA:

            {
              const int n = data->size();

              int count = 0;
              FOR_J3 {
                const string from_var = param->getString(iarg++);
                CtiRegister::CtiData * data_from = getCtiData(from_var);  // will perform expression parsing too

                bool b_break_jloop = false;

                if (data_from != NULL) {
                  switch (data_from->getType()) {
                    case D_DATA:
                      {
                        for (int i =0; i < n; ++i)
                          data->dn3(i,j) = data_from->d();

                        COUT1(" > initializing comp(" << var_name << "," << j << ") to " << from_var);
                        ++count;
                      }
                      break;
                    case DN_DATA:
                      {
                        assert(data_from->size() == n);
                        for (int i =0; i < n; ++i)
                          data->dn3(i,j) = data_from->dn(i);

                        COUT1(" > initializing comp(" << var_name << "," << j << ") to " << from_var);
                        ++count;
                      }
                      break;
                    case DN3_DATA:
                      {
                        assert(data_from->size() == n);
                        assert(j == 0);  // can only be specified on first var

                        for (int i =0; i < n; ++i)
                          for (int k=0; k<3; ++k) data->dn3(i,k) = data_from->dn3(i,k);

                        b_break_jloop = true;
                        COUT1(" > initializing " << var_name << " to " << from_var);
                        count = 3;  // all of vector processed
                      }
                      break;
                    default:
                      CWARN("cannot set component of " << data->getTypeAsString() << " from type " << data_from->getTypeAsString() << "; skipping");
                  }
                }

                if (b_break_jloop) break;  // break for loop as well
              }
              assert(count == 3);  // otherwise not all processed properly
            }

            break;

          default:

            // if you get here, we forgot to add a case.  this isn't the users
            // fault but an oversight on our part.  boo.
            assert(0);
          }

          // data has been set so update the flag -- is the value important?

          data->setFlag(1);

        } else {
          // TODO: just make a warn?
          CERR( " > unable to process INIT command; var " << var_name << " is not registered data.");

        }
      }
    }

    if (checkParam("INIT_CV_DATA_IN_GEOM")) COUT1("processing INIT_CV_DATA_IN_GEOM commands:");
    FOR_PARAM_MATCHING("INIT_CV_DATA_IN_GEOM") {

      // INIT_CV_DATA_IN_GEOM GEOM <geom> INIT <var-name> <value> [INIT <var-name> <value> ...]
      // INIT_CV_DATA_IN_GEOM GEOM BOX 0 1 0 1 0 1 INIT rho 1.0 INIT u 0.0 0.0 0.0 INIT p 1.0

      int geom = -1;
      SimpleGeom* simple_geom = NULL;
      string sbin_filename;
      string iso_var_name;
      double iso_var_value = 0.0;
      vector<Triple<pair<CtiRegister::CtiData*,string>,string,int> > data_vec;  // data/name to set, var specified, int indicates dn3 component
      // vector<pair<CtiRegister::CtiData*,int> > data_vec;
      vector<int> i_vec;
      vector<double> d_vec;

      int iarg = 0;
      while (iarg < param->size()) {
        string token = param->getString(iarg++);
        if (token == "GEOM") {
          token = param->getString(iarg++);
          if (token == "ALL") {
            geom = ALL_GEOM;
          }
          else if (token == "ISO") {
            geom = USERDEF_GEOM;
            iso_var_name = param->getString(iarg++);
            iso_var_value = param->getDouble(iarg++); // iso_var_value
          }
          else if (token == "SBIN") {
            geom = FILE_GEOM;
            sbin_filename = param->getString(iarg++);
          }
          else {
            geom = SIMPLE_GEOM;
            --iarg;
            assert(simple_geom == NULL);
            simple_geom = newSimpleGeom(&(*param),iarg);
            if (simple_geom == NULL) {
              CERR("unrecognized GEOM token in initFromParams: " << token);
            }

          }
        }

        else if (token == "INIT") {
          const string var_name = param->getString(iarg++);

          if ( CtiRegister::CtiData* data = CtiRegister::getRegisteredCtiData(var_name) ) {

            const int type = data->getType();
            if ((data->getUnindexedTopology() != CV_DATA)||((type != IN_DATA)&&(type != DN_DATA)&&(type != DN3_DATA))) {
              CERR("INIT_CV_DATA_IN_GEOM INIT \"" << var_name << "\" does not evaluate to CV IN,DN or DN3 data.")
            }
            else {
              // store pointer to data and another data to hold i,d, or d3 val
              switch ( type ) {
                case IN_DATA:
                  data_vec.push_back(Triple<pair<CtiRegister::CtiData*,string>,string,int>(pair<CtiRegister::CtiData*,string>(data,var_name),param->getString(iarg++),-1));
                  break;

                case DN_DATA:
                  data_vec.push_back(Triple<pair<CtiRegister::CtiData*,string>,string,int>(pair<CtiRegister::CtiData*,string>(data,var_name),param->getString(iarg++),-1));
                  break;

                case DN3_DATA:
                  {
                    FOR_I3 {
                      const string from_name = param->getString(iarg++);
                      CtiRegister::CtiData * from_data = CtiRegister::getCtiData(from_name);

                      if (from_data->getType() == DN3_DATA) {
                        assert(i == 0);
                        data_vec.push_back(Triple<pair<CtiRegister::CtiData*,string>,string,int>(pair<CtiRegister::CtiData*,string>(data,var_name),from_name,-1));
                        break;
                      }
                      else {
                        data_vec.push_back(Triple<pair<CtiRegister::CtiData*,string>,string,int>(pair<CtiRegister::CtiData*,string>(data,var_name),from_name,i));
                      }
                    }
                  }
                  break;

                default:

                  // if you get here, we forgot to add a case.  this isn't the users
                  // fault but an oversight on our part.  boo.

                  assert(0);
              }
            }
          }
          else {
            CERR(" > unable to process INIT_CV_DATA_IN_GEOM command; var " << var_name << " is not registered data.");
          }
        }

        else {
          CWARN(" initFromParams: skipping unrecognized unrecognized INIT_CV_DATA_IN_GEOM token " << token);
        }
      }

      if ( geom == -1 ) {
        CERR(" > INIT_CV_DATA_IN_GEOM missing/unsupported GEOM ");
      }

      // flag cells inside geom...
      // TODO: why is this int8?

      int8 * cv_flag = new int8[ncv];
      if (geom == ALL_GEOM) {
        FOR_ICV cv_flag[icv] = 1;
        COUT1(" > geom: ALL");
      }
      else if (geom == USERDEF_GEOM) {
        double * iso_var_no = new double[nno];
        // have to rebuild every time in case the iso definition is time dependent
        setNoDN(iso_var_no,iso_var_name); // constructs iso from string evaluation
        FOR_ICV cv_flag[icv] = 0;
        FOR_INO iso_var_no[ino] -= iso_var_value;
        flagCvsUnderIso(cv_flag,iso_var_no,1);
        delete[] iso_var_no;
        COUT1(" > geom: ISO");
      }
      else if (geom == FILE_GEOM) {
        SurfaceShm surface;
        surface.readBinary(sbin_filename);
        surface.initPointIsInside(x_cv,ncv);
        FOR_ICV {
          if (surface.pointIsInside(x_cv[icv]))
            cv_flag[icv] = 1;
          else
            cv_flag[icv] = 0;
        }
        COUT1(" > geom: SBIN");
      }
      else {
        assert(geom == SIMPLE_GEOM);
        assert(simple_geom);
        FOR_ICV {
          if (simple_geom->pointIsInside(x_cv[icv]))
            cv_flag[icv] = 1;
          else
            cv_flag[icv] = 0;
        }
        delete simple_geom; simple_geom = NULL;
        COUT1(" > geom: SIMPLE_GEOM");
      }

      // now set variables in flagged regions...

      for (int ii = 0, lim = data_vec.size(); ii < lim; ++ii) {
        assert(data_vec[ii].first.first->getUnindexedTopology() == CV_DATA);
        switch ( data_vec[ii].first.first->getType() ) {
          case IN_DATA:
            {
              CtiRegister::CtiData * set_var = CtiRegister::getCtiData(data_vec[ii].second);
              if (set_var != NULL) {
                switch (set_var->getType()) {
                  case I_DATA:
                    FOR_ICV if (cv_flag[icv] == 1) data_vec[ii].first.first->in(icv) = set_var->i();
                    COUT1("    > initializing " << data_vec[ii].first.second << " to " << data_vec[ii].second);
                    break;
                  case IN_DATA:
                    FOR_ICV if (cv_flag[icv] == 1) data_vec[ii].first.first->in(icv) = set_var->in(icv);
                    COUT1("    > initializing " << data_vec[ii].first.second << " to " << data_vec[ii].second);
                    break;
                  default:
                    CWARN("cannot set " << data_vec[ii].first.second << " from " << data_vec[ii].second << " (type: " << set_var->getTypeAsString() << "); skipping");
                }
                data_vec[ii].first.first->setFlag(1);
              }
              else {
                CWARN("cannot set " << data_vec[ii].first.second << " from " << data_vec[ii].second << " (data not found); skipping");
              }
            }
            break;

          case DN_DATA:
            {
              CtiRegister::CtiData * set_var = CtiRegister::getCtiData(data_vec[ii].second);
              if (set_var != NULL) {
                switch (set_var->getType()) {
                  case I_DATA:
                    {
                      FOR_ICV if (cv_flag[icv] == 1) data_vec[ii].first.first->dn(icv) = (double)set_var->i();
                      COUT1("    > initializing " << data_vec[ii].first.second << " to " << data_vec[ii].second);
                    }
                    break;
                  case IN_DATA:
                    {
                      FOR_ICV if (cv_flag[icv] == 1) data_vec[ii].first.first->dn(icv) = (double)set_var->in(icv);
                      COUT1("    > initializing " << data_vec[ii].first.second << " to " << data_vec[ii].second);
                    }
                    break;
                  case D_DATA:
                    {
                      FOR_ICV if (cv_flag[icv] == 1) data_vec[ii].first.first->dn(icv) = set_var->d();
                      COUT1("    > initializing " << data_vec[ii].first.second << " to " << data_vec[ii].second);
                    }
                    break;
                  case DN_DATA:
                    {
                      FOR_ICV if (cv_flag[icv] == 1) data_vec[ii].first.first->dn(icv) = set_var->dn(icv);
                      COUT1("    > initializing " << data_vec[ii].first.second << " to " << data_vec[ii].second);
                    }
                    break;
                  default:
                    CWARN("cannot set " << data_vec[ii].first.second << " from " << data_vec[ii].second << " (type: " << set_var->getTypeAsString() << "); skipping");
                }
                data_vec[ii].first.first->setFlag(1);
              }
              else {
                CWARN("cannot set " << data_vec[ii].first.second << " from " << data_vec[ii].second << " (data not found); skipping");
              }
            }
            break;

          case DN3_DATA:
            {
              CtiRegister::CtiData * set_var = CtiRegister::getCtiData(data_vec[ii].second);
              if (set_var != NULL) {
                switch (set_var->getType()) {
                  case I_DATA:
                    {
                      assert(data_vec[ii].third != -1);
                      FOR_ICV if (cv_flag[icv] == 1) {
                        data_vec[ii].first.first->dn3(icv,data_vec[ii].third) = set_var->i();
                      }
                      COUT1("    > initializing comp(" << data_vec[ii].first.second << "," << data_vec[ii].third << ") to " << data_vec[ii].second);
                    }
                    break;
                  case IN_DATA:
                    {
                      assert(data_vec[ii].third != -1);
                      FOR_ICV if (cv_flag[icv] == 1) {
                        data_vec[ii].first.first->dn3(icv,data_vec[ii].third) = set_var->in(icv);
                      }
                      COUT1("    > initializing comp(" << data_vec[ii].first.second << "," << data_vec[ii].third << ") to " << data_vec[ii].second);
                    }
                    break;
                  case D_DATA:
                    {
                      assert(data_vec[ii].third != -1);
                      FOR_ICV if (cv_flag[icv] == 1) {
                        data_vec[ii].first.first->dn3(icv,data_vec[ii].third) = set_var->d();
                      }
                      COUT1("    > initializing comp(" << data_vec[ii].first.second << "," << data_vec[ii].third << ") to " << data_vec[ii].second);
                    }
                    break;
                  case DN_DATA:
                    {
                      assert(data_vec[ii].third != -1);
                      FOR_ICV if (cv_flag[icv] == 1) {
                        data_vec[ii].first.first->dn3(icv,data_vec[ii].third) = set_var->dn(icv);
                      }
                      COUT1("    > initializing comp(" << data_vec[ii].first.second << "," << data_vec[ii].third << ") to " << data_vec[ii].second);
                    }
                    break;
                  case DN3_DATA:
                    {
                      assert(data_vec[ii].third == -1);
                      FOR_ICV if (cv_flag[icv] == 1) {
                        FOR_J3 data_vec[ii].first.first->dn3(icv,j) = set_var->dn3(icv,j);
                      }
                      COUT1("    > initializing " << data_vec[ii].first.second << " to " << data_vec[ii].second);
                    }
                    break;
                  default:
                    {
                      if (data_vec[ii].third != -1) {
                        CWARN("cannot set comp(" << data_vec[ii].first.second << "," << data_vec[ii].third << ") from " << data_vec[ii].second << " (type: " << set_var->getTypeAsString() << "); skipping");

                      }
                      else {
                        CWARN("cannot set " << data_vec[ii].first.second << " from " << data_vec[ii].second << " (type: " << set_var->getTypeAsString() << "); skipping");
                      }
                    }
                }
                data_vec[ii].first.first->setFlag(1);
              }
              else {
                CWARN("cannot set " << data_vec[ii].first.second << " from " << data_vec[ii].second << " (data not found); skipping");
              }
            }
            break;

          default:

            // if you get here, we forgot to add a case.  this isn't the users
            // fault but an oversight on our part. boo.
            assert(0);

        }
      }

      delete[] cv_flag;

    }
  }

  virtual void loadBalanceHook(StripedMesh* sm) {
    COUT1("StaticSolver::loadBalanceHook()");
    // hook provided to set the costs/weights in the cv-based load balancer.

    // parse infile...
    FOR_PARAM_MATCHING("LOAD_BALANCE") {
      string token = param->getString(0);
      if (token == "CV") {
        cv_lb_cost = param->getInt(1);
      }
      else if (token == "FA") {
        fa_lb_cost = param->getInt(1);
      }
      else if (token == "EF") {
        ef_lb_cost = param->getInt(1);
      }
      else if ( (token == "BFZONE") || (token == "FAZONE") ) {
        string zone_name = param->getString(1);
        map<const string,int>::iterator iter = bfZoneNameMap.find(zone_name);
        if (iter != bfZoneNameMap.end()) {
          const int izone = iter->second;
          assert(bfZoneVec[izone].getName() == zone_name);
          bfZoneVec[izone].lb_cost = param->getInt(2);
        }
      }
      else if ((token == "LP")||(token == "LSP")) {
        string lp_name = param->getString(1);
        lpHelperNameMap.find(lp_name);
        map<const string,int>::iterator iter = lpHelperNameMap.find(lp_name);
        if (iter != lpHelperNameMap.end()) {
          const int ii = iter->second;
          assert(lpHelperVec[ii].name == lp_name);
          lpHelperVec[ii].lb_cost = param->getInt(2);
        }
      }
    }

  }

  void ensureCvAdt() {

    if (cvAdt == NULL) buildCvAdt();

  }

  void buildCvAdt();

  // ===========================================================================
  // probe infrastructure...
  // ===========================================================================

  void doProbes(const int step = -1,const double time = 0.0) {

    // loop through all of the probe params and init them if necessary by adding
    // them to their respective maps (ppMap, mfpMap, fpMap, etc..)

    FOR_PARAM_MATCHING("PROBE") {
      processPointProbe(&(*param),step,time);
    }

    FOR_PARAM_MATCHING("MULTIFLUX_PROBE") {
      processMultiFluxProbe(&(*param),step,time);
    }

    FOR_PARAM_MATCHING("FLUX_PROBE") {
      processFluxProbe(&(*param),step,time);
    }

    FOR_PARAM_MATCHING("VOLUMETRIC_PROBE") {
      processVolumetricProbe(&(*param),step,time);
    }

    FOR_PARAM_MATCHING("SMD_VOLUMETRIC_PROBE") {
      processSMDVolumetricProbe(&(*param),step,time);
    }

    FOR_PARAM_MATCHING("FAZONE_PROBE") {
      processFazoneProbe(&(*param),step,time);
    }

    FOR_PARAM_MATCHING("PDF_PROBE") {
      processPdfProbe(&(*param),step,time);
    }

    FOR_PARAM_MATCHING("CONDITIONAL_PROBE") {
      processConditionalProbe(&(*param),step,time);
    }

    FOR_PARAM_MATCHING("SMD_CONDITIONAL_PROBE") {
      processSMDConditionalProbe(&(*param),step,time);
    }

    FOR_PARAM_MATCHING("POINTCLOUD_PROBE") {

      // use the full param as the hash for the map..

      if ( pcpMap.find(param->str()) == pcpMap.end()) {

        // couldnt find this probe so we need to insert it to the list

        PointCloudProbe *this_sp = new PointCloudProbe();
        pair<map<const string,PointCloudProbe*>::iterator,bool> ret =
          pcpMap.insert( pair<const string,PointCloudProbe*>(param->str(),this_sp));
        assert( ret.second);

        const int ierr = initPointCloudProbe(ret.first->second,&(*param));

        if ( ierr != 0) {

          // see the note on the error handling logic above for the PointProbe

          if ( ret.first->second->interval > 0) {
            CERR( " > unable to add PointCloudProbe from: \n >   " << param->str() );
          } else {
            CWARN( " > unable to add PointCloudProbe from: \n >  " << param->str() );
            delete[] ret.first->second;
            pcpMap.erase(ret.first);
          }

        } else {

          // if the interval is invalid, then its one-and-done for this probe

          if ( ret.first->second->interval == -1) {

            writePointCloudProbe(ret.first->second,step,time);
            delete[] ret.first->second;
            pcpMap.erase(ret.first);

          }
        }
      }
    }

    // cycle through each map of probes and call doProbe if (step%interval == 0)...
    for (map<const string,PointProbe>::iterator iter = ppMap.begin(); iter != ppMap.end(); ++iter) {
      if ((step == -1)||(step%iter->second.interval == 0))
        iter->second.doProbe(step,time);
    }

    for (map<const string,MultiFluxProbe>::iterator iter = mfpMap.begin(); iter != mfpMap.end(); ++iter) {
      if ((step == -1)||(step%iter->second.interval == 0))
        iter->second.doProbe(step,time);
    }

    for (map<const string,FluxProbe>::iterator iter = fpMap.begin(); iter != fpMap.end(); ++iter) {
      if ((step == -1)||(step%iter->second.interval == 0))
        iter->second.doProbe(step,time);
    }

    for (map<const string,VolumetricProbe>::iterator iter = vpMap.begin(); iter != vpMap.end(); ++iter) {
      if ((step == -1)||(step%iter->second.interval == 0))
        iter->second.doProbe(step,time);
    }

    for (map<const string,VolumetricProbe>::iterator iter = svpMap.begin(); iter != svpMap.end(); ++iter) {
      if ((step == -1)||(step%iter->second.interval == 0))
        iter->second.doSMDProbe(step,time);
    }

    for (map<const string,FazoneProbe>::iterator iter = fzpMap.begin(); iter != fzpMap.end(); ++iter) {
      if ((step == -1)||(step%iter->second.interval == 0))
        iter->second.doProbe(step,time);
    }

    for (map<const string,ConditionalProbe>::iterator iter = cpMap.begin(); iter != cpMap.end(); ++iter) {
      // interval == -1 will force a skip of this probe: it is left in the map, howver so it is not
      // re-initialized above...
      if ((step%iter->second.write_interval != -1) && ((step == -1)||(step%iter->second.sample_interval == 0)) ) // note use of sample interval here
        iter->second.doProbe(step,time);
    }

    for (map<const string,ConditionalProbe>::iterator iter = scpMap.begin(); iter != scpMap.end(); ++iter) {
      // interval == -1 will force a skip of this probe: it is left in the map, howver so it is not
      // re-initialized above...
      if ((step%iter->second.write_interval != -1) && ((step == -1)||(step%iter->second.sample_interval == 0)) ) // note use of sample interval here
        iter->second.doSMDProbe(step,time);
    }

 
    for (map<const string,PdfProbe>::iterator iter = pdfpMap.begin(); iter != pdfpMap.end(); ++iter) {
      // interval == -1 will force a skip of this probe: it is left in the map, howver so it is not
      // re-initialized above...
      if ((step%iter->second.write_interval != -1) && ((step == -1)||(step%iter->second.sample_interval == 0)) ) // note use of sample interval here
        iter->second.doProbe(step,time);
    }

    for (map<const string,PointCloudProbe*>::iterator iter = pcpMap.begin(); iter != pcpMap.end(); ++iter) {
      if ((step == -1)||(step%iter->second->interval == 0))
        writePointCloudProbe(iter->second,step,time);
    }

  }

  void flushProbes();

  void flushFwhSurfaces() {

    // flush everyone who's status is either ACTIVE or CHECK...

    for (map<const string,FwhSurface>::iterator iter = fwhSurfaceMap.begin(); iter != fwhSurfaceMap.end(); ++iter) {
      if (iter->second.status != CTI_STATUS_ERROR) {
        iter->second.flush();
      }
    }

  }

  // ===========================================================================
  // PROBE NAME=probe/Probe_xm2.5 INTERVAL=100 GEOM=CIRCLE_XRN -2.50 0.48 36 VARS=p u T
  // ===========================================================================
  int initPointProbe(PointProbe& pp,Param * param);
  void findClosestBf(double x_bf_seed[3], int& ibf_zone, int &izone, const double xp_seed[3]) const;
  void findClosestBf(double (*x_bf_seed)[3], int * ibf_zone, int * izone, const double (*xp_seed)[3],const int np) const;
  void findClosestBfOnZones(double (*x_bf_seed)[3], int * ibf_zone, int * izone, const double (*xp_seed)[3],const int np,const set<int>& zonesVec) const;

  // ===========================================================================
  // MULTIFLUX_PROBE NAME=three XP -0.045 0 0 NP 1 0 0 DN 0.1 VARS mdot_fa() mdot_fa()*avg_fa(c)
  // ===========================================================================
  int initMultiFluxProbe(MultiFluxProbe& mfp,Param * param);

  int initFluxProbe(FluxProbe& fp, Param* param);
  int initVolumetricProbe(VolumetricProbe& vp, Param* param);
  int initSMDVolumetricProbe(VolumetricProbe& svp, Param* param);
  int initFazoneProbe(FazoneProbe& fzp, Param* param);
  int initPdfProbe(PdfProbe& pdfp,Param* param,const bool killfile);
  int initConditionalProbe(ConditionalProbe& cp,Param* param,const bool killfile);
  int initSMDConditionalProbe(ConditionalProbe& scp,Param* param,const bool killfile);
  int initPointCloudProbe(PointCloudProbe* pcp,Param* param);
  void writePointCloudProbe(PointCloudProbe* pcp,const int step = 0,const double time = 0.0);
  int initFwhSurface(FwhSurface& fs,Param * param,const bool killfile);

  // note we removed default value from step, so that we could use step == 0 to skip restart writing...
  int processStep(const int step,const double time = 0.0,const double dt = 1.0,const bool verbose = false) {

    // processStep should normally return 0. returning -1 or -2 indicates that
    // the calling process should take action:
    // -1 : soft "stop", i.e. stop and write final result, flush any buffers, etc...
    // -2 : hard "stop!", i.e. stop IMMEDIATELY and release resources

    // for the more memory-heavy diagnostics, we cycle ownership of any buffered
    // data through the nodes, then flush when all nodes are full. Recall the
    // number of nodes is in MpiStuff as "mpi_size_internode", and the rank
    // in mpi_comm of the first rank on each node is in "rank_of_rank_internode"...

    static int fwh_write_node = 0;

    // -----------------------------------------------------------------------------------
    // every existing diagnostic should come into this routine with its status in CHECK,
    // or its status in ERROR because it is essentially a placeholder to prevent reparsing
    // and erroneous param...
    // TODO: apply to other diagnostics eventually...
    // -----------------------------------------------------------------------------------

    for (map<const string,FwhSurface>::iterator iter = fwhSurfaceMap.begin(); iter != fwhSurfaceMap.end(); ++iter) {
      assert((iter->second.status == CTI_STATUS_CHECK)||(iter->second.status == CTI_STATUS_ERROR));
    }

    // stats...

    // conditional resetting of statistics
    if ( (reset_stats_time != HUGE_VAL) && (time >= reset_stats_time) && ((time-dt) < reset_stats_time) ) {
      COUT1(" > resetting statistics (TIME threshold: " << reset_stats_time << " hit)");
      CtiRegister::resetStats();
      reset_stats_time = HUGE_VAL;  // don't process again
    }
    if ((reset_stats_int > 0) && (step%reset_stats_int == 0) ) {
      COUT1(" > resetting statistics (INTERVAL criterion: " << reset_stats_int << ")");
      CtiRegister::resetStats();
    }
    CtiRegister::updateStats(dt,verbose);

    // here we loop on all params, which is more efficient that looping over
    // all params each time we are looking for a given param...

    FOR_ALL_PARAM {

      if (param->name == "WRITE_DATA") {
        processWriteData(&*param,step,time);
      }
      else if (param->name == "WRITE_IMAGE") {
        processWriteImage(&*param,step,time);
      }
      else if (param->name == "FWH_SURFACE") {
        // -----------------------------------------------------------------------------------
        // example (note multiple geoms possible to define one FwhSurface):
        // FWH_SURFACE NAME fwh_surf_0 GEOM BOX -0.5 0.5 -0.5 0.5 -0.5 0.5
        //        GEOM BOX 0.5 0.7 -0.25 0.25 -0.25 0.25 INTERVAL 5 ENDCAPS_X 5 0.35 0.7 TEC
        // -----------------------------------------------------------------------------------
        // we use the full param as the hash for the map...
        map<const string,FwhSurface>::iterator iter = fwhSurfaceMap.find(param->str());
        if (iter == fwhSurfaceMap.end()) {
          // couldn't find this one so we need to insert it into the map...
          pair<map<const string,FwhSurface>::iterator,bool> ret = fwhSurfaceMap.insert(pair<const string,FwhSurface>(param->str(),FwhSurface()));
          assert(ret.second); // guaraneed by the find above
          // and try and initialize...
          int ierr = initFwhSurface(ret.first->second,&(*param),false); // killfile == false
          if (ierr != 0) {
            // if this PARAM failed in any way, then LEAVE it in the param map, but set its
            // status to ERROR. This will ensure it not re-parsed in each timestep, but instead
            // skipped efficiently...
            CWARN( " > unable to add FwhSurface from param: \n > " << param->str() );
            ret.first->second.status = CTI_STATUS_ERROR;
          }
          else {
            // created sucessfully, so make ACTIVE...
            ret.first->second.status = CTI_STATUS_ACTIVE;
          }
        }
        else if (iter->second.status == CTI_STATUS_CHECK) {
          // param was found, so switch to ACTIVE. This switch allows
          // us to identify when params associated with diagnostics have
          // been removed...
          iter->second.status = CTI_STATUS_ACTIVE;
        }
        else {
          // only other option is error...
          if (iter->second.status != CTI_STATUS_ERROR) {
            CERR("one or more FWH_SURFACE parameters are identical");
          }
        }
      }

    } // FOR_ALL_PARAM

    // flush images is based on the number of buffered images (not step!)
    // TODO figure out more/better metrics for flushing
    // flush images once there are enough of them...

    const int flush_images = getIntParam("FLUSH_IMAGES",1);
    if (int(sceneList.size()) >= flush_images) {
      flushImages();
    }

    // FWH flushing...
    // use this as a model for all flushing eventually...

    // -----------------------------------------------------------------------------------
    // step 1: any diagnostics still with CHECK status should be flushed and deleted...
    // -----------------------------------------------------------------------------------

    map<const string,FwhSurface>::iterator iter = fwhSurfaceMap.begin();
    while(iter != fwhSurfaceMap.end()) {
      if (iter->second.status == CTI_STATUS_CHECK) {
        iter->second.flush();
        map<const string,FwhSurface>::iterator iter_copy = iter++;
        fwhSurfaceMap.erase(iter_copy);
      }
      else {
        assert((iter->second.status == CTI_STATUS_ACTIVE)||(iter->second.status == CTI_STATUS_ERROR));
        ++iter;
      }
    }

    // -----------------------------------------------------------------------------------
    // step 2: process ACTIVE diagnostics...
    // -----------------------------------------------------------------------------------

    for (map<const string,FwhSurface>::iterator iter = fwhSurfaceMap.begin(); iter != fwhSurfaceMap.end(); ++iter) {
      if (iter->second.status == CTI_STATUS_ACTIVE) {
        if ((step == -1)||(step%iter->second.interval == 0))  {
          iter->second.doProbe(step,time,rank_of_rank_internode[fwh_write_node]);
          ++fwh_write_node;
          if (fwh_write_node == mpi_size_internode) {
            flushFwhSurfaces();
            fwh_write_node = 0;
          }
        }
        iter->second.status = CTI_STATUS_CHECK;
      }
    }

    // if this step is divisible by flush_fwh, then flush...

    const int flush_fwh = getIntParam("FLUSH_FWH",-1); // default of -1 means the fwh surfaces will fill all nodes before flushing
    if ((flush_fwh > 0)&&(step%flush_fwh == 0))
      flushFwhSurfaces();

    // move this up to same loop as fwh above eventually...

    doProbes(step,time);

    // flush probes is based on number of steps
    // for now, flush probes every so often...
    // need to figure out how to

    int flush_probes = getIntParam("FLUSH_PROBES",1000);

    if (step%flush_probes == 0)
      flushProbes();

    processSolverSpecificDiagnostics(step,time,verbose);

    // skip result writing on first step. note, we changed basic postpro
    // default step to 1, so it still can write out a file even if
    // the user does not register a step.
    if (step != 0) {
      FOR_PARAM_MATCHING("WRITE_RESULT") {
        processWriteResult(&*param,step,time);
      }
      // ----------------------------------------------
      // for WRITE_SLES we only parse on the first time,
      // as if we were reading the params in interpreted order
      // like surfer and stitch. On subsequet times, use the write_sles_interval,
      // if set, to determine if we write...
      // ----------------------------------------------
      if (write_sles_interval == -1) {
	FOR_PARAM_MATCHING("WRITE_SLES") {
	  processWriteSles(&*param,step,false); // false: NOT killfile
	}
	// if still -1 at this point, then set to zero meaning no sles...
	if (write_sles_interval == -1) write_sles_interval = 0;
      }
      // no "else if" here because we may need to process the result...
      if ((write_sles_interval > 0)&&(step%write_sles_interval == 0)) {
	writeSles(step);
      }
    }

    FOR_PARAM_MATCHING("WRITE_SNAPSHOT") {
      processWriteSnapshot(&*param,step,time);
    }

    //cout << "HACKHACKHACK" << endl;
    //dumpParams();

    // look for killfile...
    while (Param * param = kfr.getParam("killcharles")) {
      if (mpi_rank == 0) cout << "\n > processing param \"" << param->str() << "\"" << endl;
      const string token = MiscUtils::toUpperCase(param->getName());
      if (token == "STOP") {
        return -1;
      }
      else if (token == "STOP!") {
        return -2;
      }
      else {
        processParam(param,step,time);
      }
    }

    return 0;

  }

  void writeSles(const int step) {
    char filename[128];
    sprintf(filename,"%s.%08d.sles",write_sles_prefix.c_str(),step);
    if (mpi_rank == 0) cout << " Writing sles file \"" << filename << "\"..." << endl;
    writeData(filename);
    if (mpi_rank == 0) {
      unlink("sles");
      symlink(filename,"sles");
    }
    MPI_Barrier(mpi_comm);
  }

  virtual void processSolverSpecificDiagnostics(const int step, const double time, const bool verbose) {}

  virtual void reportTiming(Param * param) {
    if (mpi_rank==0){
      cout << "Timing report not available from this solver" << endl;

      string prefix = "timingReport";
      if (param->size()>1 && param->getString(0) == "NAME"){
        prefix = param->getString(1);
      }
      string filename = prefix+".json";
      MiscUtils::mkdir_for_file(filename);
      const string tmp_filename = MiscUtils::makeTmpPrefix(filename);
      FILE * fp = fopen(tmp_filename.c_str(),"w");
      assert(fp);
      fprintf(fp,"{}");
      fclose(fp);
      remove(filename.c_str());
      rename(tmp_filename.c_str(),filename.c_str());
    }
  }

  virtual void processParam(Param * param, const int step, const double time) {

    bool b_help = false;
    string token = MiscUtils::toUpperCase(param->getName());
    if (token == "HELP") {
      b_help = true;
      if (param->size() == 0) {
        helpHelp(true); // true == write header
        return;
      }
      token = MiscUtils::toUpperCase(param->getString(0));  // first token is keyword for help
    }

    if (token == "STOP") {
      // stop and stop! are now moved to kfr loop directly, however their "help" should
      // be handled here...
      assert(b_help);
      helpStop();
    }
    else if (token == "STOP!") {
      assert(b_help);
      helpStop();
    }
    else if (token == "STATS") {
      processStats(param,b_help);
    }
    else if (token == "RESET_STATS") {
      processResetStats(param,b_help);
    }
    else if (token == "EVAL") {
      processEval(param,b_help);
    }
    else if (token == "WRITE_SLES") {
      processWriteSles(param,step,true);
    }

    else if (token == "SET_LOG") {
      processSetLog(param,b_help);
    }
    else if (token == "WRITE_PLOT") {
      if (b_help) {
        processWritePlotHelp(param);
      }
      else {
        processWritePlot(param,step,time,true);
      }
    }

    // TODO:
    // these ones still need to be formatted properly...

    else if (token == "WRITE_PARAMS") {
      processWriteParams(param);
    }
    else if (token == "WRITE_JSON") {
      processWriteJson(param);
    }
    else if (token == "WRITE_IMAGE") {
      // pass true indicating it is a killfile request
      processWriteImage(param,step,time,true);
    }
    else if (token == "WRITE_DATA") {
      processWriteData(param,step,time,true);
    }
    else if (token == "WRITE_RESULT") {
      processWriteResult(param,step,time,true);
    }
    else if (token == "WRITE_SNAPSHOT") {
      processWriteSnapshot(param,step,time,true);
    }
    else if (token == "PROBE") {
      processPointProbe(param,step,time,true);
    }
    else if (token == "MULTIFLUX_PROBE") {
      processMultiFluxProbe(param,step,time,true);
    }
    else if (token == "FLUX_PROBE") {
      processFluxProbe(param,step,time,true);
    }
    else if (token == "VOLUMETRIC_PROBE") {
      processVolumetricProbe(param,step,time,true);
    }
    else if (token == "SMD_VOLUMETRIC_PROBE") {
      processSMDVolumetricProbe(param,step,time,true);
    }
    else if (token == "FAZONE_PROBE") {
      processFazoneProbe(param,step,time,true);
    }
    else if (token == "PDF_PROBE") {
      processPdfProbe(param,step,time,true);
    }
    else if (token == "CONDITIONAL_PROBE") {
      processConditionalProbe(param,step,time,true);
    }
    else if (token == "SMD_CONDITIONAL_PROBE") {
      processSMDConditionalProbe(param,step,time,true);
    }
    else if (token == "ADD_PARAM") {
      processAddParam(param);
    }
    else if (token == "RM_PARAM") {
      processRmParam(param);
    }
    else if ( token == "FLUSH_ALL") {
      flushProbes();
      flushImages();
      flushFwhSurfaces();
    }
    else if ( token == "FLUSH_PROBES") {
      flushProbes();
    }
    else if ( token == "FLUSH_IMAGES") {
      flushImages();
    }
    else if ( token == "REPORT_TIMING") {
      reportTiming(param);
    }
    else {
      if (b_help) {
        WUI(WARN,"HELP not available for unrecognized param: " << token);
      }
      else {
        WUI(WARN,"skipping unrecognized param: " << token);
      }
    }

  }

  virtual void helpHelp(const bool first) {
    // helpHelp is virtual because
    if (first) {
      WUI(WARN,"HELP should be followed by the command you want help with: e.g. HELP WRITE_IMAGE\n" <<
          "HELP is available for:");
    }
    WUI(WARN,"RESET_STATS stop stop! EVAL...");
  }

  void helpStop() {
    WUI(INFO,"stop and stop! both stop the solver after it completes the current iteration. stop (without\n" <<
        "the exclamation mark) has the same behavior as when the solver finds an emty killfile, writing a\n" <<
        "final solution file (sles) and then stopping. stop! skips this final write and stops immediately.");
  }

  void processStats(Param * param,const bool b_help) {
    if (b_help) {
      WUI(INFO,"STATS <var> [more vars...] adds the requested variables to the solver.");
    }
    else {
      registerStats(param,true); // true == do init
      WUI(INFO,"stats added");
    }
  }

  void processResetStats(Param * param,const bool b_help) {
    if (b_help) {
      WUI(INFO,"RESET_STATS resets all stats declared using the STATS keyword. It can be used in the\n" <<
          "input file to clear stats read in from an SLES at the start of a run, or it can be requested\n"
          "interactively during a run.");
    }
    else {
      // CtiRegister::resetStats();
      // WUI(INFO,"stats are reset");

      try {
        // if the user did not provide a time, then just reset STATS at calc startup
        if (param->size() == 0) {
          WUI(INFO," > resetting statistics...");
          CtiRegister::resetStats();
        }
        else if (param->size() == 1) {
          if (MiscUtils::toUpperCase(param->getString(0))=="OFF") {
            WUI(WARN," > removing automated stats resetting");
            reset_stats_int = -1;  // default
            reset_stats_time = HUGE_VAL;
          }
          else if (param->getBool(0)) {
            // place this second b/c getBool will fail if a string (non numeric) is entered
            // having string parsing first will waterfall naturally here
            WUI(INFO," > resetting statistics...");
            CtiRegister::resetStats();
          }
        }
        else if (param->size() == 2) {
          int iarg=0;
          while(iarg<param->size()) {
            const string token = MiscUtils::toUpperCase(param->getString(iarg++));
            if (token == "TIME") {
              double tmp_time = reset_stats_time;
              reset_stats_time = param->getDouble(iarg++);
              if (reset_stats_time <= 0) {
                WUI(WARN," > TIME expects a positive value; ignoring");
                reset_stats_time = HUGE_VAL;  // default
              }
              else {
                if (tmp_time != HUGE_VAL) {
                  WUI(INFO," > RESET_STATS TIME being changed from " << tmp_time << " to " << reset_stats_time);
                }
                else if (reset_stats_int != -1) {
                  WUI(WARN," > RESET_STATS TIME " << reset_stats_time << " being used instead of INTERVAL " << reset_stats_int << " criterion");
                  reset_stats_int = -1;
                }
                else {
                  WUI(INFO," > RESET_STATS TIME set to " << reset_stats_time);
                }
              }


            }
            else if (token == "INTERVAL") {
              int tmp_int = reset_stats_int;
              reset_stats_int = param->getInt(iarg++);
              if (reset_stats_int <= 0) {
                WUI(WARN," > INTERVAL expects a positive integer; ignoring");
                reset_stats_int = -1;  // default
              }
              else {
                if (reset_stats_time != HUGE_VAL) {
                  WUI(WARN," > RESET_STATS INTERVAL " << reset_stats_int << " being used instead of TIME " << reset_stats_time << " criterion");
                  reset_stats_time = HUGE_VAL;
                }
                else if (tmp_int != -1) {
                  WUI(INFO," > RESET_STATS INTERVAL being changed from " << tmp_int << " to " << reset_stats_int);
                }
                else {
                  WUI(INFO," > RESET_STATS INTERVAL set to " << reset_stats_int);
                }
              }
            }
          }
        }
        else {
          WUI(WARN," > expecting RESET_STATS TIME <time> or RESET_STATS INTERVAL <int>; cannot request both simultaneously");
        }
      }
      catch (int e) {
        WUI(WARN," ! could not parse RESET_STATS properly; ignoring");
      }
    }
  }

  void processEval(Param * param,const bool b_help) {
    if (b_help) {
      WUI(INFO,"EVAL <expression> tries to evaluate the passed expression, and reports the value. It should have\n" <<
          "no effect on the running solver. See CWI EXPRESION LINK.");
    }
    else if (param->size() == 0) {
      WUI(WARN,"EVAL expects an expression to try and evaluate: EVAL <expression>");
    }
    else {
      const string name = param->getString();
      CtiRegister::CtiData * data = CtiRegister::getCtiData(name);
      WUI(INFO,"EVAL " << name << ": " << *data);
    }
  }

  void processSetLog(Param * param,const bool b_help) {
    if (b_help) {
      WUI(INFO,"SET_LOG <filename> sets the log file used to collect the current solver output.");
    }
    else if (param->size() == 0) {
      WUI(WARN,"SET_LOG expects a log filename: e.g. SET_LOG <filename>");
    }
    else {
      b_logfile = true;
      logfile = param->getString();
      WUI(INFO,"logfile successfully set to " << logfile);
    }
  }

  void processWritePlotHelp(Param * param) {

    WUI(INFO,"WRITE_PLOT takes data from this simulation, other simulations, and from data on files to produce 2D plots");

  }

  void processWritePlot(Param * param,const int step = 0,const double time = 0.0,const bool killfile_request = false) {

    // TODO: need wpd eventually, just like the wid, but for now...

    bool b_interval = false;
    int interval;
    bool b_name = false;
    string name;
    bool b_u = false;
    string u;
    bool b_xcol = false;
    int xcol;
    bool b_ycol = false;
    int ycol;
    bool b_filename = false;
    string filename;
    bool b_grep = false;
    string grep;
    bool b_tail = false;
    int tail;

    int iarg = 0;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "INTERVAL") {
        interval = param->getInt(iarg++);
        if (interval <= 0) {
          CWARN(" > INTERVAL expects a positive integer");
          // invalid value entered, so treat as unset (default value)
        }
        else b_interval = true;
      }
      else if (token == "NAME") {
        name = param->getString(iarg++);
        b_name = true;
      }
      else if (token == "U") { // "using"
        u = param->getString(iarg++);
        b_u = true;
      }
      else if (token == "FILE") {
        filename = param->getString(iarg++);
        b_filename = true;
      }
      else if (token == "GREP") {
        grep = param->getString(iarg++);
        b_grep = true;
      }
      else if (token == "TAIL") {
        tail = param->getInt(iarg++);
        b_tail = true;
      }
      else {
        if (mpi_rank == 0) cout << "Error: unrecognized WRITE_PLOT2D token: " << token << endl;
      }
    }

    if (b_u) {
      assert(!b_xcol);
      assert(!b_ycol);
      const size_t pos = u.find(":");
      if (pos == string::npos) {
        // user is requesting just a single column...
        ycol = atoi(u.c_str());
        b_ycol = true;
      }
      else {
        // user is requesting 2 columns in a:b format...
        xcol = atoi(u.substr(0,pos).c_str());
        b_xcol = true;
        ycol = atoi(u.substr(pos+1).c_str());
        b_ycol = true;
      }
    }

    assert(b_ycol);

    // if the filename was not set, use the logfile...
    if (!b_filename) {
      assert(b_grep); // should be grepping the logfile
      assert(b_logfile);
      filename = logfile;
    }

    if (mpi_rank == 0) {

      double * xdata = NULL;
      double * ydata = NULL;
      int n;

      if (b_grep) {
        // with grep...
        if (b_tail) {
          // with tail...
          n = tail;
          if (b_xcol) {
            int ierr = MiscUtils::tailgrepxcol2(xdata,ydata,n,filename.c_str(),grep.c_str(),xcol,ycol);
            assert(ierr == 0);
          }
          else {
            int ierr = MiscUtils::tailgrepxcol(ydata,n,filename.c_str(),grep.c_str(),ycol);
            assert(ierr == 0);
          }
        }
        else {
          // no tail...
          if (b_xcol) {
            int ierr = MiscUtils::grepxcol2(xdata,ydata,n,filename.c_str(),grep.c_str(),xcol,ycol);
            assert(ierr == 0);
          }
          else {
            int ierr = MiscUtils::grepxcol(ydata,n,filename.c_str(),grep.c_str(),ycol);
            assert(ierr == 0);
          }
        }
      }
      else {
        // no grep...
        if (b_tail) {
          // with tail...
          if (b_xcol) {
            int ierr = MiscUtils::tailxcol2(xdata,ydata,n,filename.c_str(),xcol,ycol);
            assert(ierr == 0);
          }
          else {
            int ierr = MiscUtils::tailxcol(ydata,n,filename.c_str(),ycol);
            assert(ierr == 0);
          }
        }
        else {
          // no tail...
          if (b_xcol) {
            int ierr = MiscUtils::xcol2(xdata,ydata,n,filename.c_str(),xcol,ycol);
            assert(ierr == 0);
          }
          else {
            int ierr = MiscUtils::xcol(ydata,n,filename.c_str(),ycol);
            assert(ierr == 0);
          }
        }
      }

      // if no xdata provided, just use a 1-based index...
      if (!b_xcol) {
        assert(xdata == NULL);
        xdata = new double[n];
        for (int i = 0; i < n; ++i)
          xdata[i] = i+1;
      }

      if (killfile_request) {
        char filename[128];
        sprintf(filename,"%s.dat",name.c_str());
        cout << "WRITE_PLOT: writing 2d plot file: " << filename << endl;
        FILE * fp = fopen(filename,"w");
        for (int i = 0; i < n; ++i) {
          fprintf(fp,"%lf %lf\n",xdata[i],ydata[i]);
        }
        fclose(fp);
      }
      else {
        assert(0);
      }

      delete[] xdata;
      delete[] ydata;

    }


  }

  void writeAutoCompleteJson(FILE * fp) {
    fprintf(fp,",\n\"siData\": {");
    fprintf(fp,"\n  \"kbVersion\": \"%s\"",CTI::cti_docs_version);
    const string siData =
#include "static.siData"
      ;
    fprintf(fp,",\n  %s\n",siData.c_str());
    fprintf(fp,"}");
  }

  void processWriteJson(Param * param) {

    // bbox calc is collective...

    //if (!b_boundingBox) calcBoundingBox();

    // if we need to skip the next image, this is coolective....

    int iarg = 0;
    string prefix;
    uint json_bits = 0;

    int transform_type_ = 0;
    double r0_ = HUGE_VAL;
    double cos_theta0_ = HUGE_VAL;
    double sin_theta0_ = HUGE_VAL;
    double cyl_center_[3];
    FOR_I3 cyl_center_[i] = HUGE_VAL;

    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "NAME") {
        prefix = param->getString(iarg++);
      }
      else if (token == "BITS") {
        json_bits = uint(param->getInt(iarg++));
      }
      else if (token == "TRANSFORM_CYL_X") {
	// usage:
	// TRANSFORM_CYL_X <r0> <theta0-radians>
	// projects image and solution data into CYL_X space (x,y,z)->(x,r,r0*theta)
	transform_type_ = 1;
	r0_ = param->getDouble(iarg++);
	const double theta0 = param->getDouble(iarg++);
	cos_theta0_ = cos(theta0);
	sin_theta0_ = sin(theta0);
	// current the App does not have center info, use default
	FOR_I3 cyl_center_[i] = 0.0;
      }
      else if (token == "TRANSFORM_CYL_Z") {
	// usage:
	// TRANSFORM_CYL_Z <r0> <theta0-radians>
	// projects image and solution data into CYL_Z space (x,y,z)->(r,r0*theta,z)
	transform_type_ = 3;
	r0_ = param->getDouble(iarg++);
	const double theta0 = param->getDouble(iarg++);
	cos_theta0_ = cos(theta0);
	sin_theta0_ = sin(theta0);
	// current the App does not have center info, use default
        FOR_I3 cyl_center_[i] = 0.0;
      }
      else {
        if (mpi_rank == 0) cout << "Warning: skipping unrecognized WRITE_JSON token \"" << token << "\"" << endl;
      }

    }

    double bbox[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double xbb[3] = {0.0, 0.0, 0.0};
    double diag   = 0.0;
    if (isNewTransform(transform_type_,r0_,cos_theta0_,sin_theta0_,cyl_center_)) {
      double (*x_no_t)[3] = new double[nno][3];
      FOR_INO {
	GeomUtils::coordinateTransform(x_no_t[ino], x_no[ino], r0_, cos_theta0_, sin_theta0_,
				       cyl_center_, GeomUtils::COORD_TRANSFORM_TYPE(transform_type_));
      }
      calcTBoundingBox(x_no_t);
      getTBoundingBox(bbox);
      getTBoundingBoxCenter(xbb);
      diag = 2.0*getTBoundingBoxRmax();
      delete[] x_no_t;
    }
    else if (transform_type > 0) {
      getTBoundingBox(bbox);
      getTBoundingBoxCenter(xbb);
      diag = 2.0*getTBoundingBoxRmax();
    }
    else {
      getBoundingBox(bbox);
      getBoundingBoxCenter(xbb);
      diag = 2.0*getBoundingBoxRmax();
    }

    // may not need this additional guard: i.e. could use the "newImage" behavior
    // without any messages some day...
    if (WebUI::webUIOutput.message_flag) {
      // the processParam routine that may have pushed messages into this JSON
      // may have set the newImage to false. If this is the case, we can skip
      // the image writing, which should share the same prefix as this JSON...
      if (!WebUI::webUIOutput.newImage) {
        b_skip_image = true;
        skip_image_name = prefix;
      }
    }

    // only rank 0 write the json...

    if (mpi_rank == 0) {

      // write the JSON file...

      const string filename = prefix + ".json";
      MiscUtils::mkdir_for_file(filename);
      const string tmp_filename = MiscUtils::makeTmpPrefix(filename);
      FILE * fp = fopen(tmp_filename.c_str(),"w");
      assert(fp);

      /*
      double bbox[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      double xbb[3] = {0.0, 0.0, 0.0};
      double diag   = 0.0;
      getBoundingBox(bbox);
      getBoundingBoxCenter(xbb);
      diag = 2.0*getBoundingBoxRmax();
      */

      bool b_surface_loaded = true; // allows opening cascade client without any surface -- in this case we have one!
      bool b_volumeData = true;

      fprintf(fp,"{\n\"solver\":\"charles\"");

      // insist on particle interface...
      fprintf(fp,",\n\"b_particles_loaded\":true");

      /*
        it = added_data_fields.find("particles");
        if (it != added_data_fields.end())
        fprintf(fp,",\n\"b_particles_loaded\":true");
        else
        fprintf(fp,",\n\"b_particles_loaded\":false");
      */

      set<string>::iterator it;
      it = added_data_fields.find("interfaces");
      if (it != added_data_fields.end())
        fprintf(fp,",\n\"b_interfaces_loaded\":true");
      else
        fprintf(fp,",\n\"b_interfaces_loaded\":false");
      fprintf(fp,",\n\"b_surface_loaded\":%s",(b_surface_loaded?"true":"false"));
      fprintf(fp,",\n\"BoundingBox\":[%f,%f,%f,%f,%f,%f]",bbox[0],bbox[1],bbox[2],bbox[3],bbox[4],bbox[5]);
      fprintf(fp,",\n\"BoundingBoxCentroid\":[%f,%f,%f]",xbb[0],xbb[1],xbb[2]);
      fprintf(fp,",\n\"BoundingBoxDiagonal\":%f",diag);
      fprintf(fp,",\n\"LightLocation\":[0.039503,0,1.01189,-0.530497,-0.57,0.581893]");

      // bit-based criterion
      if (json_bits & (1<<0)) {
        // siData for ui autocompletion
        writeAutoCompleteJson(fp);
      }

      if (!WebUI::webUIOutput.empty()) {
        assert(WebUI::webUIOutput.message_flag);
        WebUI::webUIOutput.writeJson(fp);
        //WebUI::webUIOutput.clear(); // done below on all ranks...
      }

      // DP/ME: what is filetype
      //fprintf(fp,",\n\"filetype\": \"%s\"", filetype.c_str());

      fprintf(fp,",\n\"zoneAreas\":[\n");
      for (int i = 0, i_end = bfZoneVec.size(); i < i_end; ++i) {
        if (i == 0) fprintf(fp," %f",bfZoneVec[i].area_global);
        else fprintf(fp,",\n %f",bfZoneVec[i].area_global);
      }
      fprintf(fp,"\n]");
      fprintf(fp,",\n\"zoneCentroids\":[\n");
      for (int i = 0, i_end = bfZoneVec.size(); i < i_end; ++i) {
        if (i == 0) fprintf(fp," [%f,%f,%f]",bfZoneVec[i].x_global[0],bfZoneVec[i].x_global[1],bfZoneVec[i].x_global[2]);
        else fprintf(fp,",\n [%f,%f,%f]",bfZoneVec[i].x_global[0],bfZoneVec[i].x_global[1],bfZoneVec[i].x_global[2]);
      }
      fprintf(fp,"\n]");
      fprintf(fp,",\n\"zoneIds\":[\n");
      for (int i = 0, i_end = bfZoneVec.size(); i < i_end; ++i) {
        if (i == 0) fprintf(fp," %d",i);
        else fprintf(fp,",\n %d",i);
      }
      fprintf(fp,"\n]");
      fprintf(fp,",\n\"zoneNames\":[\n");
      for (int i = 0, i_end = bfZoneVec.size(); i < i_end; ++i) {
        if (i == 0) fprintf(fp," \"%s\"",bfZoneVec[i].getName().c_str());
        else fprintf(fp,",\n \"%s\"",bfZoneVec[i].getName().c_str());
      }
      fprintf(fp,"\n]");
      // for subzones, just use the zone index -- i.e. no subzoning
      // DP/ME is this correct?
      fprintf(fp,",\n\"zoneSubzones\":[\n");
      /*
        for (int i = 0, limit = ss.szozn_i.getLength(); i < limit; ++i) {
        if (i == 0) fprintf(fp," \"%d\"",ss.szozn_i[i]);
        else fprintf(fp,",\n \"%d\"",ss.szozn_i[i]);
        }
        fprintf(fp,"\n]");
      */
      for (int i = 0, i_end = bfZoneVec.size(); i <= i_end; ++i) {
        if (i == 0) fprintf(fp," \"%d\"",i);
        else fprintf(fp,",\n \"%d\"",i);
      }
      fprintf(fp,"\n]");
      /*
        if (ss.selectedSubzoneVec.size()) {
        fprintf(fp,",\n\"selectedSubzones\":[\n");
        for (int i = 0, limit = ss.selectedSubzoneVec.size(); i < limit; ++i) {
        if (i == 0) fprintf(fp," \"%d\"",ss.selectedSubzoneVec[i]);
        else fprintf(fp,",\n \"%d\"",ss.selectedSubzoneVec[i]);
        }
        fprintf(fp,"\n]");
        }
      */
      fprintf(fp,",\n\"zonePeriodicPairs\":[]");

      if (b_volumeData) {
        /*
          if (lp_ptr) {
          fprintf(fp,",\n\"particleVars\":{\n");
          fprintf(fp," \"D1\":[\n");
          vector<string> nameVec;
          for (map<const string,pair<int,int> >::iterator it = lp_data_map.begin(); it != lp_data_map.end(); ++it) {
          if (it->second.second == D_DATA) nameVec.push_back(it->first);
          }
          fprintf(fp,"  \"points\"");
          for (int ii = 0, ii_end = nameVec.size(); ii < ii_end; ++ii) {
          fprintf(fp,",\n  \"%s\"",nameVec[ii].c_str());
          }
          fprintf(fp,"\n ],\n");
          fprintf(fp," \"D2\":[\n");
          nameVec.clear();
          for (map<const string,pair<int,int> >::iterator it = lp_data_map.begin(); it != lp_data_map.end(); ++it) {
          if (it->second.second == D3_DATA) nameVec.push_back(it->first);
          }
          bool first = true;
          for (int ii = 0, ii_end = nameVec.size(); ii < ii_end; ++ii) {
          if (first) {
          fprintf(fp,"  \"mag(%s)\"",nameVec[ii].c_str());
          first = false;
          }
          else fprintf(fp,",\n  \"mag(%s)\"",nameVec[ii].c_str());
          }
          fprintf(fp,"\n  ]");
          fprintf(fp,"\n }");
          }
        */

        fprintf(fp,",\n\"volumeVars\":{\n");
        fprintf(fp," \"D1\":[\n");
        fprintf(fp,"  \"mesh\"");
        vector<string> nameVec;
        CtiRegister::setRegisteredCvAndNoDNNames(nameVec);
        for (int ii = 0, ii_end=nameVec.size(); ii < ii_end; ++ii) {
          fprintf(fp,",\n  \"%s\"",nameVec[ii].c_str());
        }
        fprintf(fp,"\n ],\n");
        fprintf(fp," \"D2\":[\n");
        nameVec.clear();
        CtiRegister::setRegisteredCvAndNoDN3Names(nameVec);
        bool first = true;
        for (int ii = 0, ii_end=nameVec.size(); ii < ii_end; ++ii) {
          if (first) {
            fprintf(fp,"  \"mag(%s)\"",nameVec[ii].c_str());
            first = false;
          }
          else fprintf(fp,",\n  \"mag(%s)\"",nameVec[ii].c_str());
        }
        fprintf(fp,"\n  ]");
        fprintf(fp,"\n }");

        fprintf(fp,",\n\"boundaryVarsAndFxns\":{\n");
        fprintf(fp," \"D1\":[\n");
        nameVec.clear();
        set<string> nameSet; // strip out zones and add wildcard
        CtiRegister::setRegisteredBfDNNames(nameVec);
        first = true;
        for (int ii = 0,ii_end=nameVec.size(); ii < ii_end; ++ii) {
          size_t colon_pos = nameVec[ii].find(":");
          assert( colon_pos != string::npos);
          string var_name  = nameVec[ii].substr(colon_pos+1,nameVec[ii].size()-colon_pos);
          nameSet.insert("*:"+var_name);
        }
        for (set<string>::const_iterator cit = nameSet.begin(); cit != nameSet.end(); ++cit) {
          if (first) {
            fprintf(fp,"  \"%s\"",cit->c_str());
            first = false;
          }
          else fprintf(fp,",\n  \"%s\"",cit->c_str());
        }
        nameVec.clear();
        nameSet.clear();
        // for now we only include boundary condition functions that take void and return DN...
        CtiRegister::setRegisteredFuncEvalNames(nameVec);
        for (int ii = 0,ii_end=nameVec.size(); ii < ii_end; ++ii) {
          size_t colon_pos = nameVec[ii].find(":");
          assert( colon_pos != string::npos);
          string var_name  = nameVec[ii].substr(colon_pos+1,nameVec[ii].size()-colon_pos);
          nameSet.insert("*:"+var_name);
        }
        for (set<string>::const_iterator cit = nameSet.begin(); cit != nameSet.end(); ++cit) {
          fprintf(fp,",\n  \"%s\"",cit->c_str());
        }
        fprintf(fp,"\n ],\n");
        fprintf(fp," \"D2\":[\n");
        nameVec.clear();
        nameSet.clear();
        first = true;
        CtiRegister::setRegisteredBfDN3Names(nameVec);
        for (int ii = 0, ii_end=nameVec.size(); ii < ii_end; ++ii) {
          size_t colon_pos = nameVec[ii].find(":");
          assert( colon_pos != string::npos);
          string var_name  = nameVec[ii].substr(colon_pos+1,nameVec[ii].size()-colon_pos);
          nameSet.insert("*:"+var_name);
        }
        for (set<string>::const_iterator cit = nameSet.begin(); cit != nameSet.end(); ++cit) {
          if (first) {
            fprintf(fp,"  \"mag(%s)\"",cit->c_str());
            first = false;
          }
          else fprintf(fp,",\n  \"mag(%s)\"",cit->c_str());
        }
        fprintf(fp,"\n  ]");
        fprintf(fp,"\n }");
      }
      fprintf(fp,"\n}\n");

      fclose(fp);
      remove(filename.c_str());
      rename(tmp_filename.c_str(),filename.c_str());

    }

    // all need to clear...
    WebUI::webUIOutput.clear();

  }

  void processWriteParams(Param * param) {

    // only rank 0 write the json...

    if (mpi_rank == 0) {

      assert(param->getString(0) == "NAME");
      const string prefix = param->getString(1);

      //cout << "JUST GOT WRITE_PARAMS command: " << param->getString() << endl;

      // write the JSON file...
      const string filename = prefix+".in";
      MiscUtils::mkdir_for_file(filename);
      const string tmp_filename = MiscUtils::makeTmpPrefix(filename);
      FILE * fp = fopen(tmp_filename.c_str(),"w");
      assert(fp);

      stringstream ss;
      FOR_ALL_PARAM {
        param->dump(ss);
      }

      fprintf(fp,"%s\n",ss.str().c_str());
      fclose(fp);
      remove(filename.c_str());
      rename(tmp_filename.c_str(),filename.c_str());

    }
  }

  void processAddParam(Param * param) {

    Params::addParamFromAddParam(param);

  }

  void processRmParam(Param * param) {

    Params::rmParam(param);

  }

  void averageBfToNo(double * v_no,const CtiRegister::CtiData * data);
  void averageBfToNo(double * v_no,const vector<CtiRegister::CtiData*> data_vec);
  void averageBfToNo(double (*v_no)[3],const CtiRegister::CtiData * data);
  void averageBfToNo(double (*v_no)[3],const vector<CtiRegister::CtiData*> data_vec);

  void averageCvToNo(double * v_no,const CtiRegister::CtiData * iso);
  void averageCvToNo(double * v_no,const double * v_cv);
  void averageCvToNo(double (*v_no)[3],const CtiRegister::CtiData * iso);
  void averageCvToNo(double (*v_no)[3],const double (*v_cv)[3]);

  //
  // geometry stuff...
  // this is an exact copy of the stuff in SimpleSurface. Maybe there could
  // be a base class called renderable?
  //

  void calcBoundingBox() {

    assert(!b_boundingBox);
    double buf[6] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };

    assert(x_no);
    FOR_INO {
      buf[0] = min(buf[0],x_no[ino][0]);
      buf[1] = min(buf[1],-x_no[ino][0]);
      buf[2] = min(buf[2],x_no[ino][1]);
      buf[3] = min(buf[3],-x_no[ino][1]);
      buf[4] = min(buf[4],x_no[ino][2]);
      buf[5] = min(buf[5],-x_no[ino][2]);
    }

    double buf_global[6];
    MPI_Allreduce(buf,buf_global,6,MPI_DOUBLE,MPI_MIN,mpi_comm);

    boundingBox[0] =  buf_global[0];
    boundingBox[1] =  -buf_global[1];
    boundingBox[2] =  buf_global[2];
    boundingBox[3] =  -buf_global[3];
    boundingBox[4] =  buf_global[4];
    boundingBox[5] =  -buf_global[5];

    if (mpi_rank == 0) {
      cout << "StaticSolver::calcBoundingBox(): X " << boundingBox[0] << ":" << boundingBox[1] <<
        ", Y " << boundingBox[2] << ":" << boundingBox[3] <<
        ", Z " << boundingBox[4] << ":" << boundingBox[5] << endl;
    }

    b_boundingBox = true;

  }

  void getBoundingBox(double _bBox[6]) {
    if (!b_boundingBox) calcBoundingBox();
    assert(b_boundingBox);
    _bBox[0] = boundingBox[0];
    _bBox[1] = boundingBox[1];
    _bBox[2] = boundingBox[2];
    _bBox[3] = boundingBox[3];
    _bBox[4] = boundingBox[4];
    _bBox[5] = boundingBox[5];
  }

  double getBoundingBoxRmax() {
    if (!b_boundingBox) calcBoundingBox();
    assert(b_boundingBox);
    return 0.5*sqrt((boundingBox[0]-boundingBox[1])*(boundingBox[0]-boundingBox[1])+
                    (boundingBox[2]-boundingBox[3])*(boundingBox[2]-boundingBox[3])+
                    (boundingBox[4]-boundingBox[5])*(boundingBox[4]-boundingBox[5]));
  }

  void getBoundingBoxCenter(double _bBoxCenter[3]) {
    if (!b_boundingBox) calcBoundingBox();
    assert(b_boundingBox);
    _bBoxCenter[0] = 0.5*(boundingBox[0]+boundingBox[1]);
    _bBoxCenter[1] = 0.5*(boundingBox[2]+boundingBox[3]);
    _bBoxCenter[2] = 0.5*(boundingBox[4]+boundingBox[5]);
  }

  void calcTBoundingBox(const double (*x_no_t)[3]) {

    assert(!b_tBoundingBox);
    double buf[6] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };

    FOR_INO {
      buf[0] = min(buf[0],x_no_t[ino][0]);
      buf[1] = min(buf[1],-x_no_t[ino][0]);
      buf[2] = min(buf[2],x_no_t[ino][1]);
      buf[3] = min(buf[3],-x_no_t[ino][1]);
      buf[4] = min(buf[4],x_no_t[ino][2]);
      buf[5] = min(buf[5],-x_no_t[ino][2]);
    }

    double buf_global[6];
    MPI_Allreduce(buf,buf_global,6,MPI_DOUBLE,MPI_MIN,mpi_comm);

    tBoundingBox[0] =  buf_global[0];
    tBoundingBox[1] =  -buf_global[1];
    tBoundingBox[2] =  buf_global[2];
    tBoundingBox[3] =  -buf_global[3];
    tBoundingBox[4] =  buf_global[4];
    tBoundingBox[5] =  -buf_global[5];

    if (mpi_rank == 0) {
      cout << "StaticSolver::calcTBoundingBox(): X " << tBoundingBox[0] << ":" << tBoundingBox[1] <<
        ", Y " << tBoundingBox[2] << ":" << tBoundingBox[3] <<
        ", Z " << tBoundingBox[4] << ":" << tBoundingBox[5] << endl;
    }

    b_tBoundingBox = true;

  }

  void getTBoundingBox(double _bBox[6]) {
    //if (!b_tBoundingBox) calcTBoundingBox(x_no_t);
    assert(b_tBoundingBox);
    _bBox[0] = tBoundingBox[0];
    _bBox[1] = tBoundingBox[1];
    _bBox[2] = tBoundingBox[2];
    _bBox[3] = tBoundingBox[3];
    _bBox[4] = tBoundingBox[4];
    _bBox[5] = tBoundingBox[5];
  }

  double getTBoundingBoxRmax() {
    //if (!b_tBoundingBox) calcTBoundingBox(x_no_t);
    assert(b_tBoundingBox);
    return 0.5*sqrt((tBoundingBox[0]-tBoundingBox[1])*(tBoundingBox[0]-tBoundingBox[1])+
                    (tBoundingBox[2]-tBoundingBox[3])*(tBoundingBox[2]-tBoundingBox[3])+
                    (tBoundingBox[4]-tBoundingBox[5])*(tBoundingBox[4]-tBoundingBox[5]));
  }

  void getTBoundingBoxCenter(double _bBoxCenter[3]) {
    //if (!b_tBoundingBox) calcTBoundingBox(x_no_t);
    assert(b_tBoundingBox);
    _bBoxCenter[0] = 0.5*(tBoundingBox[0]+tBoundingBox[1]);
    _bBoxCenter[1] = 0.5*(tBoundingBox[2]+tBoundingBox[3]);
    _bBoxCenter[2] = 0.5*(tBoundingBox[4]+tBoundingBox[5]);
  }

  void processWriteResult(Param* param, const int step = 0, const double time = 0.0, const bool killfile_request = false) {

    string prefix = "result"; // default..
    int interval  = -1;       // the default behavior will not be to write a result..

    // if the user requests a result with "WRITE_RESULT" using the killfile, we will assume they want a one-off
    if ((param->size() == 0)&&(!killfile_request)) {
      CWARN(" >  write_result syntax should be: \n \
            WRITE_RESULT INTERVAL <interval> [NAME <prefix>] or \n \
            WRITE_RESULT <interval>");
      return;
    }
    else if ( param->size() == 1) {
      interval = param->getInt(0);
      if (interval <= 0) {
        CWARN(" > INTERVAL expects a positive integer; ignoring");
        // invalid value entered, so treat as unset (default value)
        interval = -1;
      }

      if (step%interval != 0) return;
    }
    else {
      int iarg = 0;
      while ( iarg < param->size()) {
        string token = param->getString(iarg++);
        if ( token == "NAME") {
          prefix = param->getString(iarg++);
        }
        else if (token == "INTERVAL") {
          interval = param->getInt(iarg++);
          if (interval <= 0) {
            CWARN(" > INTERVAL expects a positive integer; ignoring");
            // invalid value entered, so treat as unset (default value)
            interval = -1;
          }

          if ( step%interval != 0) return;
        }
        else {
          CWARN(" > unrecognized token in WRITE_RESULT : " << token << "; returning...");
          return;
        }
      }
    }

    // interval < 0 means that it wasnt specified-- error here..
    if ((interval < 0)&&(!killfile_request)) {
      CWARN(" >  write_result syntax should be: \n \
            WRITE_RESULT INTERVAL <interval> [NAME <prefix>] or \n      \
            WRITE_RESULT <interval>");
      return;
    }

    char filename[128];
    sprintf(filename,"%s.%08d.sles",prefix.c_str(),step);

    if (mpi_rank == 0) cout << " Writing solution file \"" << filename << "\"..." << endl;
    writeData(filename);

  }

  void processWriteSles(Param* param,const int step,const bool killfile_request) {

    if (param->size() == 0) {
      // ----------------------------------------------------
      // WRITE_SLES
      // ----------------------------------------------------
      if (killfile_request) {
	writeSles(step);
	return;
      }
      else {
	CWARN(" > WRITE_SLES syntax should be:\n" <<
	      "WRITE_SLES INTERVAL <interval> [NAME <prefix>] or \n" <<
	      "WRITE_SLES <interval>");
	if (write_sles_interval == -1) write_sles_interval = 0; // means no WRITE_SLES
	return;
      }
    }
    else if (param->size() == 1) {
      // ----------------------------------------------------
      // WRITE_SLES 1000
      // ----------------------------------------------------
      // if the interval has never been set, then reset it...
      const int interval = param->getInt(0);
      if (interval <= 0) {
	CWARN(" > WRITE_SLES <interval> expects positive integer");
	if (write_sles_interval == -1) write_sles_interval = 0; // means no WRITE_SLES
	return;
      }
      write_sles_interval = interval;
      // if this happens to be from a killfile, then write the current
      // sles filename if the step is divisible by write_sles_interval. Otherwise
      // this just sets the step...
      if (killfile_request && (step%write_sles_interval == 0))
	writeSles(step);
      return;
    }
    else {
      // ----------------------------------------------------
      // WRITE_SLES NAME blah INTERVAL 1000
      // ----------------------------------------------------
      string name;
      bool b_name = false;
      int interval;
      bool b_interval = false;
      int iarg = 0;
      while (iarg < param->size()) {
        const string token = param->getString(iarg++);
        if (token == "NAME") {
          name = param->getString(iarg++);
          b_name = true;
        }
        else if (token == "INTERVAL") {
          interval = param->getInt(iarg++);
          if (interval <= 0) {
            CWARN(" > INTERVAL expects a positive integer; ignoring");
            // invalid value entered, so treat as unset (default value)
          }
          else b_interval = true;

        }
        else {
          CWARN(" > unrecognized token in WRITE_SLES: " << token << "; skipping...");
        }
      }
      if ((!b_interval)||(interval <= 0)) {
        CWARN(" > check WRITE_SLES syntax:\n" <<
        "WRITE_SLES INTERVAL <interval> [NAME <prefix>] or \n" <<
        "WRITE_SLES <interval>");
        if (write_sles_interval == -1) write_sles_interval = 0; // means no WRITE_SLES
        return;
      }
      write_sles_interval = interval;
      if (b_name) write_sles_prefix = name;
      if (killfile_request && (step%write_sles_interval == 0))
      writeSles(step);
      return;
    }

  }

  // methods to initialize the probe if not in map and flush single request if from a killfile
  void processPointProbe(Param* param,const int step=0,const double time = 0.0,const bool killfile_request = false);

  void processMultiFluxProbe(Param* param,const int step=0,const double time = 0.0,const bool killfile_request = false);

  void processFluxProbe(Param* param,const int step=0,const double time = 0.0,const bool killfile_request = false);

  void processVolumetricProbe(Param* param,const int step=0,const double time = 0.0,const bool killfile_request = false);

  void processSMDVolumetricProbe(Param* param,const int step=0,const double time = 0.0,const bool killfile_request = false);

  void processFazoneProbe(Param* param,const int step=0,const double time = 0.0,const bool killfile_request = false);

  void processPdfProbe(Param* param,const int step=0,const double time = 0.0,const bool killfile_request = false);

  void processConditionalProbe(Param* param,const int step=0,const double time = 0.0,const bool killfile_request = false);

  void processSMDConditionalProbe(Param* param,const int step=0,const double time = 0.0,const bool killfile_request = false);


  void processWriteSnapshot(Param* param, const int step = 0, const double time = 0.0, const bool killfile_request = false) {

    if ( param->size() == 0) {

      CWARN(" >  write_snapshot syntax should be: \n \
            WRITE_SNAPSHOT INTERVAL <interval> NAME <prefix> VARS <var-list>");

    } else {

      int iarg           = 0;
      int interval       = -1;
      string prefix_name = "";
      vector<string> var_vec;
      map<const string,Snapshot>::iterator snap_it = snapMap.end();

      while ( iarg < param->size()) {

        string token = param->getString(iarg++);
        if ( token == "NAME" ) {

          prefix_name = param->getString(iarg++);
          // check to see if we already have a snapshot with this name..
          snap_it = snapMap.find(prefix_name);
          if ( snap_it != snapMap.end())
            break;

        } else if ( token == "INTERVAL") {

          interval = param->getInt(iarg++);
          if (interval <= 0) {
            CWARN(" > INTERVAL expects a positive integer; treating as a one-off");
            // invalid value entered, so treat as unset (default value)
            interval = -1;
          }

        } else if ( token == "VARS" ) {

          while ( iarg < param->size() ) {
            token = param->getString(iarg);
            if ( (token == "INTERVAL") || (token == "NAME") ) {
              break;
            } else {

              // want to accomodate the ability to process wildcard for
              // boundary data ..

              size_t found = token.find("*:");
              if ( found != string::npos) {

                for (int izone=0, nzones = bfZoneVec.size(); izone < nzones; ++izone) {
                  string this_var = MiscUtils::replaceAll(token,"*:",bfZoneVec[izone].getName()+":");
                  var_vec.push_back(this_var);
                }

              } else {

                // plain vanilla var...
                var_vec.push_back(token);
              }

              iarg++;
            }
          }

        } else {
          CWARN(" > unrecognized WRITE_SNAPSHOT token: " << token << "; returning... \n WRITE_SNAPSHOT INTERVAL <interval> NAME <prefix> VARS <var-list>");
          return;
        }
      }

      if ( snap_it == snapMap.end() ) {

        // this snapshot has not been parsed before and so it needs to be
        // constructed and the inputs need to be checked for errors..

        if ( prefix_name == "" ) {

          CWARN( " > Error: snapshot NAME is missing; returning... \n WRITE_SNAPSHOT NAME <prefix> INTERVAL <interval> VARS <var-list>");

        } else if ( interval < 0 && !killfile_request) {

          CWARN( " > Error: invalid or missing INTERVAL; returning... \n WRITE_SNAPSHOT NAME <prefix> INTERVAL <interval> VARS <var-list>");

        } else {

          // this is valid snapshot entry and needs to be added to the map
          pair<map<const string,Snapshot>::iterator, bool> ret
            = snapMap.insert(pair<const string,Snapshot>(prefix_name,Snapshot()));

          // we check for membership earlier, so the insertion must have succeeded.
          assert( ret.second == true);

          snap_it = ret.first;
          snap_it->second.interval = interval;
          snap_it->second.name     = prefix_name;

          const int nvar = var_vec.size();
          for (int ivar = 0; ivar < nvar; ++ivar) {
            CtiRegister::CtiData * data = CtiRegister::getRegisteredCtiData(var_vec[ivar],false);

            if ( data != NULL) {

              if ( ((data->getUnindexedTopology() != CV_DATA)&&(data->getUnindexedTopology() != BF_DATA)&&(data->getUnindexedTopology() != LP_DATA)) ||
                   ((data->getType() != DN_DATA)&&(data->getType() != DN3_DATA)) ) {

                if (mpi_rank == 0 )
                  cout << " > Warning : SNAPSHOT var does not evaluate to CV_DN,CV_DN3,BF_DN,BF_DN3,LP_DN,LP_DN3: " << var_vec[ivar] << endl;
                continue;

              } else if ( ( (data->getBits() & WRITE_DATA) == 0) && ( (data->getBits() & CAN_WRITE_DATA) == 0) ) {

                if ( mpi_rank == 0 )
                  cout << " > Warning: SNAPSHOT registered var does not have WRITE_DATA privileges ; skipping " << var_vec[ivar] << endl;
                continue;

              }

              snap_it->second.var_vec.push_back(pair<string,CtiRegister::CtiData*>( var_vec[ivar], data));

            } else {

              // unregistered data so it needs to be evaluated..

              data = CtiRegister::getUnregisteredCtiData(var_vec[ivar]);

              if ( (data == NULL) || ( (data->getUnindexedTopology() != CV_DATA) && (data->getUnindexedTopology() != BF_DATA)) ||
                   ( (data->getType() != DN_DATA) && (data->getType() != DN3_DATA))) {

                if ( data == NULL) {
                  if ( mpi_rank == 0)
                    cout << " Warning: SNAPSHOT unable to find : " << var_vec[ivar] << endl;

                } else {

                  if (mpi_rank == 0 ) {
                    cout << " Warning: SNAPSHOT var does not evaluate to CV_DN, CV_DN3, BF_DN, BF_DN3: " << var_vec[ivar] << endl;
                  }
                }

                continue;

              }

              snap_it->second.var_vec.push_back(pair<string,CtiRegister::CtiData*>( var_vec[ivar],NULL));
            }
          }
        }
      }

      if ( step%snap_it->second.interval == 0) {

        writeSnapshot_(step,time,&(snap_it->second));

      } else if (snap_it->second.interval == -1 && killfile_request) {

        // assume it is a one-off if it is a killfile request without an interval

        writeSnapshot_(step,time,&(snap_it->second));
        snapMap.erase(snap_it);

      }
    }
  }

  void writeSnapshot_(const int step, const double time, const Snapshot* snap) {
    UNUSED(time);
    // should not have inited if none of the data could be found

    const int nvar = snap->var_vec.size();
    if ( nvar == 0 )
      return;

    const double wtime0 = MPI_Wtime();

    char filename[128];
    sprintf(filename,"%s.%08d.sles",snap->name.c_str(),step);
    MiscUtils::mkdir_for_file_collective(filename,0);
    string tmp_filename = MiscUtils::makeTmpPrefix(filename);
    MPI_File_delete(tmp_filename.c_str(),MPI_INFO_NULL);

    MPI_File fh;
    MPI_File_open(mpi_comm,tmp_filename.c_str(),MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    if ( mpi_rank == 0 ) {
      int itmp[2] = {UGP_IO_MAGIC_NUMBER+1, 5}; // convention for data and not mesh data..
      MPI_File_write(fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
    }

    MPI_Offset offset = int_size*2;

    // write out the step and time ...
    // could do this manually, but going through the data map so that
    // the write command is cosistent..
    {
      CtiRegister::CtiData* cti_data = CtiRegister::getRegisteredCtiData("time",false);
      assert( cti_data != NULL); // must have a time registered...
      cti_data->writeData("time",fh,offset);
    }

    {
      CtiRegister::CtiData* cti_data = CtiRegister::getRegisteredCtiData("step",false);
      assert( cti_data != NULL);
      cti_data->writeData("step",fh,offset);
    }

    set<int> lpSet;
    for (int ivar = 0; ivar < nvar; ++ivar) {

      CtiRegister::CtiData* cti_data;
      int write_lpocv = -1;
      bool skip_var = false; // used to skip xp
      if (snap->var_vec[ivar].second != NULL) {

        cti_data = snap->var_vec[ivar].second; // data is registered..

        if (cti_data->getUnindexedTopology() == LP_DATA) {
          string name = snap->var_vec[ivar].first;
          size_t colon_pos = name.find(":");
          assert(colon_pos != string::npos); // names must be prepended with particle class name
          const string lp_name  = name.substr(0,colon_pos);
          const string var_name  = name.substr(colon_pos+1,name.size()-colon_pos);
          if ((var_name == "xp")||(var_name == "dp"))
            skip_var = true;
          map<const string,int>::iterator iter2 = lpHelperNameMap.find(lp_name);
          assert(iter2 != lpHelperNameMap.end());
          const int lp_index = iter2->second;
          if (lpSet.find(lp_index) == lpSet.end()) {
            write_lpocv = lp_index;
            lpSet.insert(lp_index);
          }
        }

      } else {

        // otherwise the data needs to be evaluated..
        cti_data = CtiRegister::getUnregisteredCtiData(snap->var_vec[ivar].first);

        // we also need to prepare the unregistered data for writing..
        // start by checking the topology..

        if ( cti_data->getUnindexedTopology() == CV_DATA)  {

          cti_data->setBits(WRITE_DATA);

          if ( cti_data->xora_ptr == NULL) {

            assert( cti_data->dde_ptr == NULL);
            cti_data->setDdeStuff(cvora_striped,dde_cv_striped);

          }

        } else if ( cti_data->getUnindexedTopology() == BF_DATA) {

          const int izone = cti_data->getIndex();
          cti_data->setBits(WRITE_DATA);

          if ( cti_data->xora_ptr == NULL ) {

            assert( cti_data->dde_ptr ==  NULL);
            cti_data->setDdeStuff(bfZoneVec[izone].bfora_striped,bfZoneVec[izone].dde_striped);

          }

        } else if (cti_data->getUnindexedTopology() == LP_DATA) {

          string name = snap->var_vec[ivar].first;
          size_t colon_pos = name.find(":");
          assert(colon_pos != string::npos); // names must be prepended with particle class name
          const string lp_name  = name.substr(0,colon_pos);
          const string var_name  = name.substr(colon_pos+1,name.size()-colon_pos);
          if ((var_name == "xp")||(var_name == "dp"))
            skip_var = true;
          map<const string,int>::iterator iter2 = lpHelperNameMap.find(lp_name);
          assert(iter2 != lpHelperNameMap.end());
          const int lp_index = iter2->second;
          if (lpSet.find(lp_index) == lpSet.end()) {
            write_lpocv = lp_index;
            lpSet.insert(lp_index);
          }

        }

      }

      if (write_lpocv != -1) {
        assert((write_lpocv >= 0)&&(write_lpocv < lpHelperVec.size()));
        writeLpocvAndInitDdeStuff(fh,offset,write_lpocv);

        // write :xp and :dp by default (needed for imaging - like a minimalist grid)
        CtiRegister::CtiData* cti_data2 = CtiRegister::getRegisteredCtiData(lpHelperVec[write_lpocv].name+":xp",false);
        assert(cti_data2 != NULL);
        cti_data2->writeData(lpHelperVec[write_lpocv].name+":xp",fh,offset);
        cti_data2 = CtiRegister::getRegisteredCtiData(lpHelperVec[write_lpocv].name+":dp",false);
        assert(cti_data2 != NULL);
        cti_data2->writeData(lpHelperVec[write_lpocv].name+":dp",fh,offset);
      }

      assert( cti_data != NULL);
      if (!skip_var)
        cti_data->writeData(snap->var_vec[ivar].first,fh,offset);

    }

    if ( mpi_rank == 0 ) {
      cout << " > EOF" << endl;
      Header header;
      header.id = UGP_IO_EOF;
      sprintf(header.name,"EOF");
      header.skip = header_size;
      MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
    }

    offset += header_size;

    MPI_File_set_size(fh,offset);
    MPI_File_close(&fh);
    MPI_Barrier(mpi_comm);

    if ( mpi_rank == 0) {
      const double seconds = MPI_Wtime() - wtime0;
      cout << " > write file : " << filename << " ; size: " << double(offset)/1.0E+9 << " [GB]; write rate: "
           << double(offset)/(1.0E+9*seconds) << " [GB/sec]" << endl;
    }

    if (mpi_rank == 0) {
      remove(filename);
      rename(tmp_filename.c_str(),filename);
    }

  }

  void flushImages() {

    // flush all images to disk...

    for (list<CtiScene*>::iterator iter = sceneList.begin(); iter != sceneList.end(); ++iter) {
      (*iter)->flushImage();
      delete *iter;
    }
    sceneList.clear();

  }

  //virtual void buildSolverSurface(vector<SimpleTri>& triVec,const string& surface_name) {}

  bool isNewTransform(const int transform_type_,const double r0_,const double cos_theta0_,const double sin_theta0_,const double cyl_center_[3]) {

    if (transform_type_ == 0) {
      // no transform
      transform_type = 0;
      return false;
    }
    else if ((transform_type_ == transform_type) && (fabs(r0_-r0)<1.0e-13) &&
	     (fabs(cos_theta0_-cos_theta0)<1.0e-13) && (fabs(sin_theta0_-sin_theta0)<1.0e-13) &&
	     (DIST(cyl_center_,cyl_center)<1.0e-13)) {
      // the same transform as before
      return false;
    }
    else {
      b_tBoundingBox = false;
      transform_type = transform_type_;
      r0 = r0_;
      cos_theta0 = cos_theta0_;
      sin_theta0 = sin_theta0_;
      FOR_I3 cyl_center[i] = cyl_center_[i];
      return true;
    }
  }

  void doImage(CtiScene* scene,const int step, const double time,const bool killfile_request) {
    if (scene->skip(step) || scene->isCached() ) return;

    // we could have the imaged class pass a "this" into the scene so the
    // scene would call these only when required.

    int transform_type_;
    double r0_, cos_theta0_, sin_theta0_;
    double cyl_center_[3];
    scene->getTransformInfo(transform_type_, r0_, cos_theta0_, sin_theta0_, cyl_center_);

    if (isNewTransform(transform_type_,r0_,cos_theta0_,sin_theta0_,cyl_center_)) {
      double (*x_no_t)[3] = new double[nno][3];
      FOR_INO {
	scene->coordinateTransform(x_no_t[ino], x_no[ino]);
      }
      calcTBoundingBox(x_no_t);
      scene->setRmax(getTBoundingBoxRmax());
      double bbox[6];
      getTBoundingBox(bbox);
      scene->setBoundingBbox(bbox[0],bbox[1],bbox[2],bbox[3],bbox[4],bbox[5]);
      double center[3];
      getTBoundingBoxCenter(center);
      scene->setCenter(center[0],center[1],center[2]);
      delete[] x_no_t;
    }
    else if (transform_type > 0) {
      scene->setRmax(getTBoundingBoxRmax());
      double bbox[6];
      getTBoundingBox(bbox);
      scene->setBoundingBbox(bbox[0],bbox[1],bbox[2],bbox[3],bbox[4],bbox[5]);
      double center[3];
      getTBoundingBoxCenter(center);
      scene->setCenter(center[0],center[1],center[2]);
    }
    else {
      scene->setRmax(getBoundingBoxRmax());
      double bbox[6];
      getBoundingBox(bbox);
      scene->setBoundingBbox(bbox[0],bbox[1],bbox[2],bbox[3],bbox[4],bbox[5]);
      double center[3];
      getBoundingBoxCenter(center);
      scene->setCenter(center[0],center[1],center[2]);
    }

    // this is not optimal -- can the scene know how to call zone info from its "renderable" objects?
    // even the Json file could be produced by the "renderable" class methods

    FOR_IZONE(bfZoneVec) {
      scene->addZoneName(bfZoneVec[izone].getName());
    }
    scene->convertHiddenZoneNamesToIndices();

    scene->initCanvas();

    // need to add data plane before painting so that we can know if were are going to blank...
    if (scene->hasGeomPlane()) {
      scene->setGeomPlaneAsDataPlane();
      if (scene->blankDataPlane()) scene->addGeomPlaneAsBlankPlane();
    }

    // use a static counter to round-robin assignment of image rank0...
    static int image_rank0 = 0;
    scene->setWriteRank(image_rank0++);
    if (image_rank0 == mpi_size) image_rank0 = 0;

    // need step for filename
    scene->setStep(step);

    // needed time in metadata
    scene->setImdTime(time);

    // specify hiding of elements in the scene
    // may be used by both volvis surfaces or surface vis, so do here first
    scene->setHiddenBfsFromHiddenZones(zone_bf,nbf);  // creates a bf_flag based on zones
    // superimpose cv_flag for hiding (used by mesh-plane vis)
    int * cv_flag_geom_plane = NULL;
    int * tet_flag_geom_plane = NULL;
    if (scene->hasGeomPlane() && ((scene->getVar() == "mesh") || (scene->getVar() == "fluid:mesh") || (scene->getVar() == "cht:mesh")) && (!scene->cellFlooded()) ) {
      double geom_plane_xp[3],geom_plane_np[3];
      scene->getGeomPlaneXpAndNp(geom_plane_xp,geom_plane_np);
      cv_flag_geom_plane = new int[ncv_g];  // NOTE! switched flag values
      FOR_ICV {
	double x_vv_t[3];
	scene->coordinateTransform(x_vv_t,x_vv[icv]);
	const double dist =
	  (x_vv_t[0]-geom_plane_xp[0])*geom_plane_np[0] +
	  (x_vv_t[1]-geom_plane_xp[1])*geom_plane_np[1] +
	  (x_vv_t[2]-geom_plane_xp[2])*geom_plane_np[2];
	if (dist < 0.0) {
	  cv_flag_geom_plane[icv] = 0;  // show
	}
	else {
	  cv_flag_geom_plane[icv] = 1;  // hide
	}
      }
      updateCvData(cv_flag_geom_plane); // put into ghosts..

      // overwrite periodic ghosts... why? ignore periodic faces?
      for (int icv = ncv; icv < ncv_g; ++icv) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv-ncv]);
        if (bits) {
          cv_flag_geom_plane[icv] = 2;
        }
      }
      // if cht mesh needs vis
      CtiRegister::CtiData * cht_nodes = getRegisteredCtiData("cht:T_no",false);
      if (cht_nodes != NULL) {
        int (*noote)[4] = NULL;
        int nte;
        const bool b_got_tets = cht_nodes->getTets(noote,nte);
        assert(b_got_tets);
        assert(noote != NULL);
        tet_flag_geom_plane = new int[nte];

        double (*x_no_te)[3] = NULL;
        const bool b_got_x = cht_nodes->getX(x_no_te);
        assert(x_no_te != NULL);
        const int nno_te = cht_nodes->size();
        double * dist_no = new double[nno_te];
        for (int ino = 0; ino < nno_te; ++ino) {
          double x_no_t[3];
	        scene->coordinateTransform(x_no_t,x_no_te[ino]);

          dist_no[ino] =
            (x_no_t[0]-geom_plane_xp[0])*geom_plane_np[0] +
            (x_no_t[1]-geom_plane_xp[1])*geom_plane_np[1] +
            (x_no_t[2]-geom_plane_xp[2])*geom_plane_np[2];
        }
        // find the associated helper...
        for (int ite = 0; ite < nte; ++ite) {
          int n_neg = 0;
          for (int ii=0; ii<4; ++ii) {
            if (dist_no[noote[ite][ii]] < 0.0) ++n_neg;
          }
          if ((n_neg < 4) && (n_neg > 0)) tet_flag_geom_plane[ite] = 0;  // show
          else tet_flag_geom_plane[ite] = 1;
        }
        DELETE(dist_no);
      }

      if (scene->blankDataPlane()) {
        scene->removeGeomPlaneAsBlankPlane();  // crinkle-cut shouldn't use blank plane (cv_Flag takes care of this already)
        scene->setHiddenBfsFromCvFlag(cv_flag_geom_plane,cvobf,nbf);
      }
      // if not blanking we still need the cv_flag for the mesh vis, proper clipping of iso-surfaces

    }
    // could expand to other ways of generating cv_flag...e.g. geometric restriction of vis to local region


    // -----------
    // volume vis
    // -----------
    // some preprocessing of passed parameters so that vis logic is only entered if needed
    bool b_do_volvis = false;
    CtiRegister::CtiData * data_vv = NULL;
    if (scene->hasVolvis()) {
      const string vv_var = scene->getVolvisVar();
      if (vv_var != "none") {
        data_vv = CtiRegister::getCtiData(vv_var);  // get the data we are going to volvis...
        if ((data_vv == NULL)||(data_vv->getUnindexedTopology() != CV_DATA)||(data_vv->getType() != DN_DATA)) {
          WUI(WARN,"VOLVIS requested variable \"" << scene->getVolvisVar() << "\" does not eval to CVDN data; skipping");
          data_vv = NULL;
        }
        else {
          b_do_volvis = true;
        }
      }
      else {
        if ((scene->getVolvisType() != "X-RAY") && (scene->getVolvisType() != "SURFACE")) {
          WUI(WARN,"only VOLVIS_TYPE X-RAY may be specified with VOLVIS_VAR \"none\"; skipping");
        }
        else {
          b_do_volvis = true;  // var is "none" but we are in X-RAY mode
        }
      }
    }

    // Logic portion of volvis image
    if (b_do_volvis) {
      // this transmits the volvis stuff from WriteImageData to Canvas
      scene->prepVolvis();

      // add boundary surfaces...
      // builds the canvas AND sets the in/out points on the volvis pixel lines-of-sight
      scene->addBoundaryFacesVolVis(x_no,x_bf,n_bf,zone_bf,noobf_i,noobf_v,nbf);

      // need to add all periodic internal faces...
      FOR_INTERPROC_IFA {
        const int icv1 = cvofa[ifa][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g));
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv1-ncv]);
        if (bits) scene->addPeriodicFaceVolVis(x_no,x_fa[ifa],noofa_i[ifa+1]-noofa_i[ifa],noofa_v+noofa_i[ifa]);
      }

      // if x-ray we can skip the internal face volvis because it only uses depth b/w surfaces
      if ((scene->getVolvisType() != "X-RAY") && (scene->getVolvisType() != "SURFACE")) {
        FOR_ICV scene->addVoronoiPointVolVis(x_vv[icv],r_vv[icv],data_vv->dn(icv));
      }

      // processing of write/flush is done here so we don't have to skip over all
      // other image types which are over-ridden by volvis anyway
      if ( scene->isClientRequest() ) {
        // don't add step index to filename when interacting with Cascade client...
        scene->writeImage();
      }
      else if ( killfile_request ) {
        // this image request came from a killfile command (outside of app)...
        scene->writeImage(step);
      }
      else {
        // buffer image
        scene->cacheImage();
      }
      return;
    }

    // either volvis failed or wasn't specified; now do other type of image data

    // ------------
    // surface vis
    // ------------
    // preprocessing of passed parameters so that vis logic is only entered if needed
    string surface_var = "";
    CtiRegister::CtiData * data_surf = NULL;
    bool b_has_surf_data = false;
    bool b_surf_tetdata = false;
    bool b_surf_cvdata = false;
    bool b_surf_bfdata = false;
    if (scene->getSurfaceVar(surface_var)) {
      data_surf = CtiRegister::getCtiData(surface_var);
      if (data_surf != NULL) {
        if (data_surf->hasTets()) {
          if (data_surf->getType() != DN_DATA) {
            WUI(WARN,"tet-based surface var: " << surface_var << " does not eval to DN; skipping");
          }
          else {
            b_surf_tetdata = true;
          }
        }
        else if (data_surf->getUnindexedTopology() == CV_DATA) {
          if ((data_surf->getType() == DN_DATA) || (data_surf->getType() == DN3_DATA)) {
            b_surf_cvdata = true;
          }
          else {
            WUI(WARN," > surface var: " << surface_var << " is CV_DATA but not DN_DATA/DN3_DATA; skipping");
          }
        }
        else if (data_surf->getUnindexedTopology() == BF_DATA) {
          if ((data_surf->getType() == DN_DATA) || (data_surf->getType() == DN3_DATA)) {
            b_surf_bfdata = true;
          }
          else {
            WUI(WARN," > surface var: " << surface_var << " is BF_DATA but not DN_DATA/DN3_DATA; skipping");
          }
        }
        else {
          WUI(WARN," > surface var: " << surface_var << " with topo: " << data_surf->getUnindexedTopology() << " not handled; skipping");
        }
      }
      else {
        // look for a wild card...
        size_t colon_pos = surface_var.find("*:");
        if (colon_pos != string::npos) {
          // must be wildcard bf data!
          // leave verification of existence of these vars for later
          b_surf_bfdata = true;
        }
        else {
          WUI(WARN," > surface var: " << surface_var << " not handled; skipping");
        }
      }
      // these values should be exclusive...any need to check?
      b_has_surf_data = (b_surf_tetdata || b_surf_cvdata || b_surf_bfdata);
    }

    // logic for drawing surface(s)
    if (b_has_surf_data) {
      if (b_surf_tetdata) {
        // grab the tets...
        int (*noote)[4] = NULL;
        int nte;
        const bool b_got_tets = data_surf->getTets(noote,nte);
        assert(b_got_tets);
        assert(noote != NULL);
        double (*x_no_te)[3] = NULL;
        const bool b_got_x = data_surf->getX(x_no_te);
        assert(x_no_te != NULL);
        assert(data_surf->getType() == DN_DATA);
        double * data_no = data_surf->getDNptr();
        assert(data_no != NULL);
        scene->addBoundaryTris(x_no_te,data_no,data_surf->size(),noote,nte,bfZoneVec.size()); // HACK: index = try last zone index + 1
      }
      else if (b_surf_cvdata) {
        if (data_surf->getType() == DN_DATA) {
          double * data_no = new double[nno];
          setNoDN(data_no,surface_var,data_surf);
          set<int> zones_with_data;  // empty with populated data means all zones
          scene->addBoundaryFaces(x_no,x_bf,n_bf,data_no,noobf_i,noobf_v,zone_bf,nbf,zones_with_data);
          DELETE(data_no);
        }
        else if (data_surf->getType() == DN3_DATA) {
          // vector vis
          double (*bf_data3)[3] = new double[nbf][3];
          double (*tmp)[3] = NULL;
          setBfDN3(tmp,surface_var,data_surf,bf_data3);
          scene->addBoundaryFacesWithVectors(x_no,x_bf,n_bf,bf_data3,noobf_i,noobf_v,zone_bf,area_bf,nbf);
          DELETE(bf_data3);
        }
        // condition where not DN/DN3 checked previously
      }
      else if (b_surf_bfdata) {
        // containers for various data types
        double * data_no = NULL;        // traditional bf_data (at nodes)
        double * bf_data = NULL;        // cell-flooded
        double (*bf_data3)[3] = NULL;  // vector-data

        set<int> zone_indices; // sorted for addFaces function

        // determine whether wildcard expanded or not
        size_t found = surface_var.find("*:");
        if (found != string::npos) {
          vector<string> zone_surface_vars;
          vector<CtiRegister::CtiData *> zone_data;
          for (int izone = 0, nzones=bfZoneVec.size(); izone < nzones; ++izone) {
            string this_str = MiscUtils::replaceAll(surface_var,"*:",bfZoneVec[izone].getName()+":");
            if ( (data_surf = CtiRegister::getCtiData(this_str)) ) {
              if (data_surf->getType() == DN_DATA) {
                assert(bf_data3 == NULL);
                if (scene->cellFlooded()) {
                  assert(data_no == NULL);
                  if (bf_data == NULL) bf_data = new double[nbf];
                }
                else {
                  assert(bf_data == NULL);
                  if (data_no == NULL) data_no = new double[nno];
                }
              }
              else if (data_surf->getType() == DN3_DATA) {
                assert(data_no == NULL);
                assert(bf_data == NULL);
                if (bf_data3 == NULL) bf_data3 = new double[nbf][3];
              }
              else {
                WUI(WARN," > requested surface var: " << this_str << " is BF_DATA but not DN_DATA/DN3_DATA; skipping\n");
              }
              zone_surface_vars.push_back(this_str);
              zone_data.push_back(data_surf);
              zone_indices.insert(data_surf->getIndex());
            }
          }

          if (data_no) {
            setBfDN(data_no,zone_surface_vars,zone_data);
          }
          else if (bf_data) {
            setBfDN(NULL,zone_surface_vars,zone_data,bf_data);
          }
          else if (bf_data3) {
            double (*tmp)[3] = NULL;
            setBfDN3(tmp,zone_surface_vars,zone_data,bf_data3);
          }
        }
        else {
          // single bf zone
          const int index = data_surf->getIndex();
          zone_indices.insert(index);
          if (data_surf->getType() == DN_DATA) {
            if (scene->cellFlooded()) {
              bf_data = new double[nbf];
              setBfDN(NULL,surface_var,data_surf,bf_data);
            }
            else {
              data_no = new double[nno];
              setBfDN(data_no,surface_var,data_surf);
            }
          }
          else if (data_surf->getType() == DN3_DATA) {
            bf_data3 = new double[nbf][3];
            double (*tmp)[3] = NULL;
            setBfDN3(tmp,surface_var,data_surf,bf_data3);
          }
        }

        // bf data has been allocated & populated
        // now render based on data-type
        if (data_no) {
          // bf_data put on the nodes
          scene->addBoundaryFaces(x_no,x_bf,n_bf,data_no,noobf_i,noobf_v,zone_bf,nbf,zone_indices);
        }
        else if (bf_data) {
          assert(scene->cellFlooded());

          if (subSurface) {
            scene->addSurfaceTris(subSurface->xp,subSurface->spost,subSurface->znost,subSurface->ist_global_and_bits,subSurface->nst);
          }
          else {
            double * data = NULL;
            set<int> empty_zones_set;
            scene->addBoundaryFaces(x_no,x_bf,n_bf,data,noobf_i,noobf_v,zone_bf,nbf,empty_zones_set);  // hiding already applied; empty set so draw all tris...
          }

          int ncv_c = 0;
          double *phi_c = new double[nbf];
          int *zone_vv_c = new int[nbf];
          double (*x_vv_c)[3] = new double[nbf][3];
          double *delta_cv_c  = new double[nbf];

          for (set<int>::iterator it = zone_indices.begin(); it != zone_indices.end(); ++it) {
            const int izone = *it;
            for (int ibf = bfZoneVec[izone].ibf_f; ibf <= bfZoneVec[izone].ibf_l; ++ibf) {
              assert(zone_bf[ibf] == izone);
              const int icv = cvobf[ibf];
              phi_c[ncv_c]            = bf_data[ibf];
              FOR_I3 x_vv_c[ncv_c][i] = x_vv[icv][i];
              delta_cv_c[ncv_c]       = r_vv[icv];
              zone_vv_c[ncv_c]        = izone;
              ++ncv_c;
            }
          }
          scene->addPartialSurfaceData(phi_c,x_vv_c,delta_cv_c,zone_vv_c,ncv_c,1.002);
          delete[] phi_c;
          delete[] zone_vv_c;
          delete[] x_vv_c;
          delete[] delta_cv_c;
        }
        else if (bf_data3) {
          scene->addBoundaryFacesWithPartialVectors(x_no,x_bf,n_bf,bf_data3,noobf_i,noobf_v,zone_bf,area_bf,nbf,zone_indices);
        }
        else {
          WUI(WARN," > unrecognized BF_DATA variable: " << surface_var << "; drawing surface without data\n");
          scene->addBoundaryFaces(x_no,x_bf,n_bf,data_no,noobf_i,noobf_v,zone_bf,nbf,zone_indices);  //data_no should be NULL if we are here
        }

        DELETE(data_no);
        DELETE(bf_data);
        DELETE(bf_data3);
      }
    }
    else {
      // no data on surface, so just draw it with gold/silver
      if (scene->showBfSurface() || (subSurface == NULL)) {
        set<int> zones_with_data;  // empty in this case
        double * data = NULL;      // just surface, so empty
        scene->addBoundaryFaces(x_no,x_bf,n_bf,data,noobf_i,noobf_v,zone_bf,nbf,zones_with_data);
      }
      else {
        // subsurface drawing
        scene->addSurfaceTris(subSurface->xp,subSurface->spost,subSurface->znost,subSurface->ist_global_and_bits,subSurface->nst);
      }
    }

    // ----------
    // plane vis
    // ----------

    // logic preprocessing of data plane
    bool b_has_plane = false;
    bool b_plane_cf = false;        // cell-flooded
    bool b_plane_mesh = false;      // crinkle vis
    bool b_plane_tet_mesh = false;  // crinkle vis
    bool b_plane_tetdata = false;
    CtiRegister::CtiData * data_plane;
    if (scene->hasGeomPlane()) {
      b_has_plane = true;
      if (!scene->hasVar()) {
        WUI(WARN,"missing VAR <variable name> for plane data; skipping");
        b_has_plane = false;
      }
      else if (scene->getVar() == "mesh") {
        b_plane_mesh = true;
        b_plane_tet_mesh = true;
      }
      else if (scene->getVar() == "fluid:mesh") b_plane_mesh = true;
      else if (scene->getVar() == "cht:mesh") b_plane_tet_mesh = true;
      else {
        if ( (data_plane = CtiRegister::getCtiData(scene->getVar())) ) {
          if (data_plane->hasTets()) {
            if (data_plane->getType() != DN_DATA) {
              WUI(WARN,"tet-based plane var: " << scene->getVar() << " does not eval to DN; skipping\n");
            }
            else {
              b_plane_tetdata = true;
            }
          }
          else {
            // not tet data, so better be CV data...
            if ((data_plane == NULL)||(data_plane->getUnindexedTopology() != CV_DATA)) {
              if ((data_plane->getType() != DN_DATA) && (data_plane->getType() != IN_DATA)) {
              WUI(WARN,"plane var: " << scene->getVar() << " does not eval to CVDN or CVIN; skipping\n");
              b_has_plane = false;
              }
            }
          }
        }
        else {
          b_has_plane = false;
        }
      }

      if (scene->cellFlooded()) b_plane_cf = true;
    }

    // logic for data plane vis
    if (b_has_plane) {

      double geom_plane_xp[3],geom_plane_np[3];
      scene->getGeomPlaneXpAndNp(geom_plane_xp,geom_plane_np);
      const double mag = MAG(geom_plane_np); assert(mag > 0.0);
      const double unit_np[3] = {geom_plane_np[0]/mag,geom_plane_np[1]/mag,geom_plane_np[2]/mag};

      if (b_plane_cf) {
        // a cell-flooded planar image
        double (*x_vv_c)[3] = new double[ncv][3];
        double *delta_cv_c  = new double[ncv];
        int ncv_c = 0;

        if (b_plane_mesh) {
          // "mesh" visualization
          int * index_cv_c = new int[ncv];
          FOR_ICV {
            double x_vv_t[3];
            scene->coordinateTransform(x_vv_t,x_vv[icv]);
            const double dx[3] = DIFF(x_vv_t,geom_plane_xp);
            if (fabs(DOT_PRODUCT(dx,unit_np)) <= 2.0*r_vv[icv]) {
              FOR_I3 x_vv_c[ncv_c][i] = x_vv[icv][i];
              delta_cv_c[ncv_c] = r_vv[icv];
              index_cv_c[ncv_c] = 65534-32768;
              ++ncv_c;
            }
          }

          int * ivv_global = new int[ncv_c];
          int8 * vvora = NULL;
          buildXora(vvora,ncv_c);
          for (int iv = 0; iv < ncv_c; ++iv) ivv_global[iv] = (int)vvora[mpi_rank]+iv;
          DELETE(vvora);

          scene->addPlaneMeshRvvPositive(x_vv_c,delta_cv_c,index_cv_c,ivv_global,ncv_c,1.002);

          DELETE(ivv_global);
          DELETE(index_cv_c);
        }
        else {
          // regular CVDN variable
          double * phi_c  = new double[ncv];
          double * data_ptr = NULL;
          int * data_ptr_int = NULL;
          if (data_plane->getType() == DN_DATA) {
            data_ptr = data_plane->getDNptr();
          }
          else if (data_plane->getType() == IN_DATA) {
            data_ptr_int = data_plane->getINptr();
          }
          FOR_ICV {
	    double x_vv_t[3];
	    scene->coordinateTransform(x_vv_t,x_vv[icv]);
	    const double dx[3] = DIFF(x_vv_t,geom_plane_xp);
	    if (fabs(DOT_PRODUCT(dx,unit_np)) <= 2.0*r_vv[icv]) {
	      if (data_ptr != NULL) {
	        phi_c[ncv_c] = data_ptr[icv];
	      }
	      else if (data_ptr_int != NULL) {
	        phi_c[ncv_c] = (double)data_ptr_int[icv];
	      }
	      FOR_I3 x_vv_c[ncv_c][i] = x_vv[icv][i];
	      delta_cv_c[ncv_c]       = r_vv[icv];
	      ++ncv_c;
	    }
          }
          scene->addPlaneDataRvvPositive(phi_c,x_vv_c,delta_cv_c,ncv_c,1.002);
          DELETE(phi_c);
        }

        DELETE(x_vv_c);
        DELETE(delta_cv_c);
      }
      else {
        // not cell-flooded
        if (b_plane_tetdata) {
          // a tet-based data var...
          // grab the tets...
          int (*noote)[4] = NULL;
          int nte;
          const bool b_got_tets = data_plane->getTets(noote,nte);
          assert(b_got_tets);
          assert(noote != NULL);
          double (*x_no_te)[3] = NULL;
          const bool b_got_x = data_plane->getX(x_no_te);
          assert(x_no_te != NULL);
          assert(data_plane->getType() == DN_DATA);
          double * data_var_no = data_plane->getDNptr();
          assert(data_var_no != NULL);
          const int nno_te = data_plane->size();
          double * iso_var_no = new double[nno_te];
          for (int ino = 0; ino < nno_te; ++ino) {
	    double x_no_te_t[3];
	    scene->coordinateTransform(x_no_te_t,x_no_te[ino]);
            iso_var_no[ino] =
              (x_no_te_t[0]-geom_plane_xp[0])*geom_plane_np[0] +
              (x_no_te_t[1]-geom_plane_xp[1])*geom_plane_np[1] +
              (x_no_te_t[2]-geom_plane_xp[2])*geom_plane_np[2];
          }
          vector<pair<SimpleTriWithData,int> > triVec;
          // find the associated helper...
          for (int ite = 0; ite < nte; ++ite) {
            const int ino0 = noote[ite][0];
            const int ino1 = noote[ite][1];
            const int ino2 = noote[ite][2];
            const int ino3 = noote[ite][3];
            tetCutPlaneData(triVec,iso_var_no[ino0],iso_var_no[ino1],iso_var_no[ino2],iso_var_no[ino3],
                            x_no_te[ino0],x_no_te[ino1],x_no_te[ino2],x_no_te[ino3],
                            data_var_no[ino0],data_var_no[ino1],data_var_no[ino2],data_var_no[ino3]);
          }
          scene->addSimpleInternalTrisWithData(triVec,65534);
          delete[] iso_var_no;
        }
        else if (b_plane_mesh || b_plane_tet_mesh) {
          if (b_plane_mesh) {
            assert(cv_flag_geom_plane);
            // crinkle-cut vis
            // earlier we computed the cv_flag, so simply reuse
            FOR_INTERNAL_IFA {
              const int icv0 = cvofa[ifa][0];
              const int icv1 = cvofa[ifa][1];

              if ((cv_flag_geom_plane[icv0] == 1)&&(cv_flag_geom_plane[icv1] == 0)) {
                scene->addInternalFace(x_no,x_fa[ifa],noofa_i[ifa+1]-noofa_i[ifa],noofa_v+noofa_i[ifa]);
              }
              else if ((cv_flag_geom_plane[icv0] == 0)&&(cv_flag_geom_plane[icv1] == 1)) {
                scene->addInternalFaceFlip(x_no,x_fa[ifa],noofa_i[ifa+1]-noofa_i[ifa],noofa_v+noofa_i[ifa]);
              }
            }

            FOR_INTERPROC_IFA {
              const int icv0 = cvofa[ifa][0];
              const int icv1 = cvofa[ifa][1];
              if ((cv_flag_geom_plane[icv0] == 1)&&(cv_flag_geom_plane[icv1] == 0)) {
                scene->addInternalFace(x_no,x_fa[ifa],noofa_i[ifa+1]-noofa_i[ifa],noofa_v+noofa_i[ifa]);
              }
            }
          }

          if (b_plane_tet_mesh) {
            CtiRegister::CtiData * cht_nodes = getRegisteredCtiData("cht:T_no",false);

            if (cht_nodes != NULL) {
              assert(tet_flag_geom_plane);
              int (*noote)[4] = NULL;
              int nte;
              const bool b_got_tets = cht_nodes->getTets(noote,nte);
              assert(b_got_tets);
              assert(noote != NULL);
              double (*x_no_te)[3] = NULL;
              const bool b_got_x = cht_nodes->getX(x_no_te);
              assert(x_no_te != NULL);
              const int nno_te = cht_nodes->size();

              vector<SimpleTri> triVec;
              // find the associated helper...
              for (int ite = 0; ite < nte; ++ite) {
                if (tet_flag_geom_plane[ite] == 0) {
                  const int ino0 = noote[ite][0];
                  const int ino1 = noote[ite][1];
                  const int ino2 = noote[ite][2];
                  const int ino3 = noote[ite][3];

                  triVec.push_back(SimpleTri(x_no_te[ino0],x_no_te[ino1],x_no_te[ino2]) );
                  triVec.push_back(SimpleTri(x_no_te[ino0],x_no_te[ino3],x_no_te[ino1]) );
                  triVec.push_back(SimpleTri(x_no_te[ino3],x_no_te[ino0],x_no_te[ino2]) );
                  triVec.push_back(SimpleTri(x_no_te[ino1],x_no_te[ino3],x_no_te[ino2]) );
                }
              }
              scene->addSimpleInternalTrisWithMesh(triVec,7,false,32766);  // different zone than internal fluid
            }
          }
        }
        else {
          // CV data from intersected cells
          double * iso_var_no = new double[nno];

          FOR_INO {
	    double x_no_t[3];
	    scene->coordinateTransform(x_no_t,x_no[ino]);
	    iso_var_no[ino] =
	      (x_no_t[0]-geom_plane_xp[0])*geom_plane_np[0] +
	      (x_no_t[1]-geom_plane_xp[1])*geom_plane_np[1] +
	      (x_no_t[2]-geom_plane_xp[2])*geom_plane_np[2];
          }

          double * data_plane_no = new double[nno];
          setNoDN(data_plane_no,scene->getVar(),data_plane);

          vector<pair<SimpleTriWithData,int> > triVec;
          buildIsoWithData(triVec,iso_var_no,0.0,data_plane_no);

          delete[] data_plane_no;
          delete[] iso_var_no;

          scene->addSimpleInternalTrisWithData(triVec,65534);
        }
      }
    }

    // ----------
    // iso vis
    // ----------
    // preprocessing of isosurface vis
    bool b_has_iso = false;
    bool b_iso_tetdata = false;
    CtiRegister::CtiData * data_iso = NULL;  // determines whether to render data
    if (scene->hasGeomIso() && scene->getIsoCount()) {
      // evaluate individual isosurface data in loop
      // here just check data var since iso_var may be quite varied
      b_has_iso = true;
      if (scene->hasIsoDataVar()) {
        // data for coloring plot
        if (scene->getIsoDataVar() == "none") {
          // treat as just draw isosurface without data...perhaps give single value instead (for flat coloring)?
        }
        else {
          if ( (data_iso = CtiRegister::getCtiData(scene->getIsoDataVar())) ) {
            if (data_iso->hasTets()) b_iso_tetdata = true; // WHY? -- see below
            else {
              // non-tet-based data variable
              if ((data_iso->getType() != DN3_DATA) && (data_iso->getType() != DN_DATA)) {
                WUI(WARN," > isosurface data var " << scene->getIsoDataVar() << " does not evaluate to DN_DATA/DN3_DATA; skipping data rendering on iso");
                data_iso = NULL;
              }
            }
          }
        }
      }
    }

    // logic portion of iso-surface vis
    if (b_has_iso) {
      if (data_iso == NULL) {
        // draw isos without any data (gold/silver coloring)
        for (int ii = 0; ii < scene->getIsoCount(); ++ii) {
          if (CtiRegister::CtiData * iso_var = CtiRegister::getCtiData(scene->getIsoVar(ii))) {
            if (((iso_var->getUnindexedTopology() == NO_DATA)||(iso_var->getUnindexedTopology() == CV_DATA))&&(iso_var->getType() == DN_DATA)) {
              vector<pair<SimpleTri,int> > triVec;
              if (iso_var->hasTets()) {
                //MPI_Pause("tet iso no data");

                // a cht isosurface variable
                // double * iso_var_no = iso_var->getDNptr();
                buildChtIso(triVec,iso_var,scene->getIsoValue(ii));

              }
              else {
                // non-tet-based variable
                int * hidden_cvs = (scene->blankDataPlane()) ? cv_flag_geom_plane:NULL;
                buildIso(triVec,scene->getIsoVar(ii),iso_var,scene->getIsoValue(ii),hidden_cvs);
              }

              scene->addSimpleInternalTrisWithMeshOnEdge(triVec);
            }
            else {
              WUI(WARN,"isosurface var " << scene->getIsoVar(ii) << " does not evaluate to CVDN/NODN; skipping");
            }
          }
        }
      }
      else {
        // draw isos WITH data
        for (int ii = 0; ii < scene->getIsoCount(); ++ii) {
          if (CtiRegister::CtiData * iso_var = CtiRegister::getCtiData(scene->getIsoVar(ii))) {
            if (((iso_var->getUnindexedTopology() == NO_DATA)||(iso_var->getUnindexedTopology() == CV_DATA))&&(iso_var->getType() == DN_DATA)) {
              if (data_iso->getType() == DN_DATA) {
                vector<pair<SimpleTriWithData,int> > triVec;
                if (iso_var->hasTets()) {
                  if (b_iso_tetdata) { //
                    // a cht isosurface variable
                    // double * iso_var_no = iso_var->getDNptr();
                    // double * data_iso_no = data_iso->getDNptr();
                    //MPI_Pause("tet isodata no data");
                    buildChtIsoWithData(triVec,iso_var,scene->getIsoValue(ii),data_iso);
                  }
                  else {
                    WUI(WARN,"iso var: " << scene->getIsoVar(ii) << " is CHT data but data var: " << scene->getIsoDataVar() << " is not; skipping\n");
                  }
                }
                else {
                  // non-CHT iso variable
                  if (!b_iso_tetdata) {
                    int * hidden_cvs = (scene->blankDataPlane()) ? cv_flag_geom_plane:NULL;
                    buildIsoWithData(triVec,scene->getIsoVar(ii),iso_var,scene->getIsoValue(ii),scene->getIsoDataVar(),data_iso,hidden_cvs);
                  }
                  else {
                    WUI(WARN,"iso var: " << scene->getIsoVar(ii) << " is not CHT data but data var: " << scene->getIsoDataVar() << " is; skipping\n");
                  }
                }

                scene->addSimpleInternalTrisWithData(triVec,65532);
              }
              else if (data_iso->getType() == DN3_DATA) {
                vector<SimpleTriWithWgts> triVec;
                double * iso_var_no = new double[nno];
                setNoDN(iso_var_no,scene->getIsoVar(ii),iso_var);
                double (*dn3_no)[3] = new double[nno][3];
                setNoDN3(dn3_no,scene->getIsoDataVar());
                double (*dn3_tmp)[3] = NULL;
                double *dn_tmp = NULL; // unused
                double *dn_no = NULL; // unused
                int * hidden_cvs = (scene->blankDataPlane()) ? cv_flag_geom_plane:NULL;
                buildIsoWithWgtsAndData(triVec,iso_var_no,scene->getIsoValue(ii),dn_tmp,dn_no,0,dn3_tmp,dn3_no,1,hidden_cvs);
                delete[] iso_var_no;
                delete[] dn_tmp;
                delete[] dn_no;

                int my_count = triVec.size();
                double my_area = 0.0;
                for (int it = 0; it < my_count; ++it) {
                  const double nA[3] = TRI_NORMAL_2(triVec[it].x0,triVec[it].x1,triVec[it].x2);
                  my_area += 0.5*MAG(nA);
                }
                double area;
                MPI_Allreduce(&my_area,&area,1,MPI_DOUBLE,MPI_SUM,mpi_comm);

                double my_mag_max = 0.0;
                FOR_INO my_mag_max = max(my_mag_max,MAG(dn3_no[ino]));
                double mag_max;
                MPI_Allreduce(&my_mag_max,&mag_max,1,MPI_DOUBLE,MPI_MAX,mpi_comm);

                if ((area > 0.0)&&(mag_max > 0.0)) {

                  // user/app controls (defaults are 1.0)...
                  const double factor = 0.2*getBoundingBoxRmax()*scene->getVectorScale()/mag_max;
                  const double target_area = area/(250.0*scene->getVectorDensity());

                  if (factor > 0.0) {

                    vector<double> x0Vec;
                    vector<double> dxVec;
                    double current_area = 0.0;
                    for (int it = 0; it < my_count; ++it) {

                      const double nA[3] = TRI_NORMAL_2(triVec[it].x0,triVec[it].x1,triVec[it].x2);
                      current_area += 0.5*MAG(nA);
                      if (current_area > target_area) {
                        current_area = 0.0; // reset current representative area

                        // d0
                        int i_l = triVec[it].i0_l;
                        int i_r = triVec[it].i0_r;
                        double wgt_l = triVec[it].wgt0_l;
                        double wgt_r = 1.0-wgt_l;
                        double d_l[3], d_r[3];
                        if (i_l < nno)
                          FOR_I3 d_l[i] = wgt_l*dn3_no[i_l][i];
                        else
                          FOR_I3 d_l[i] = wgt_l*dn3_tmp[i_l-nno][i];
                        if (i_r < nno)
                          FOR_I3 d_r[i] = wgt_r*dn3_no[i_r][i];
                        else
                          FOR_I3 d_r[i] = wgt_r*dn3_tmp[i_r-nno][i];
                        double d0[3]; FOR_I3 d0[i] = d_l[i]+d_r[i];
                        // d1
                        i_l = triVec[it].i1_l;
                        i_r = triVec[it].i1_r;
                        wgt_l = triVec[it].wgt1_l;
                        wgt_r = 1.0-wgt_l;
                        if (i_l < nno)
                          FOR_I3 d_l[i] = wgt_l*dn3_no[i_l][i];
                        else
                          FOR_I3 d_l[i] = wgt_l*dn3_tmp[i_l-nno][i];
                        if (i_r < nno)
                          FOR_I3 d_r[i] = wgt_r*dn3_no[i_r][i];
                        else
                          FOR_I3 d_r[i] = wgt_r*dn3_tmp[i_r-nno][i];
                        double d1[3]; FOR_I3 d1[i] = d_l[i]+d_r[i];
                        // d2
                        i_l = triVec[it].i2_l;
                        i_r = triVec[it].i2_r;
                        wgt_l = triVec[it].wgt2_l;
                        wgt_r = 1.0-wgt_l;
                        if (i_l < nno)
                          FOR_I3 d_l[i] = wgt_l*dn3_no[i_l][i];
                        else
                          FOR_I3 d_l[i] = wgt_l*dn3_tmp[i_l-nno][i];
                        if (i_r < nno)
                          FOR_I3 d_r[i] = wgt_r*dn3_no[i_r][i];
                        else
                          FOR_I3 d_r[i] = wgt_r*dn3_tmp[i_r-nno][i];
                        double d2[3]; FOR_I3 d2[i] = d_l[i]+d_r[i];
                        FOR_I3 x0Vec.push_back((triVec[it].x0[i]+triVec[it].x1[i]+triVec[it].x2[i])/3.0);
                        FOR_I3 dxVec.push_back((d0[i]+d1[i]+d2[i])/3.0);

                      }
                    }
                    triVec.clear();
                    delete[] dn3_tmp;
                    delete[] dn3_no;

                    if (x0Vec.size() > 0)
                      scene->addInternalVectors((double(*)[3])(&x0Vec[0]),(double(*)[3]  )(&dxVec[0]),65532,(int)x0Vec.size()/3,factor);
                  }
                }
              }  // end DN3
            }
            else {
              WUI(WARN,"isosurface var " << scene->getIsoVar(ii) << " does not evaluate to CVDN/NODN; skipping");
            }
          }
        }
      }
    }

    // -------------
    // particle vis
    // -------------

    if (scene->hasParticleVar()) {

      double dp_mag = scene->getParticleMagnification();
      string particle_var = scene->getParticleVar();

      if (checkParam("NEW_PARTICLE_VIS")) {

        if (mpi_rank == 0) cout << "NEW_PARTICLE_VIS: particle_var \"" << particle_var << "\", dp_mag: " << dp_mag << endl;

        CtiRegister::CtiData * particle_data = CtiRegister::getCtiData(particle_var);
        if (particle_data != NULL) {
          // see if this data is associated with a position...
          double (*xp)[3] = NULL;
          if (particle_data->getX(xp)) {
            assert(xp != NULL);
            const int np = particle_data->size();
            // is there a length scale?...
            double * delta = NULL;
            if (particle_data->getL(delta)) {
              // we have X and L, so this is enough to render...
              if (particle_data->getType() == DN_DATA) {
                double * data = particle_data->getDNptr();
                for (int ip = 0; ip < np; ++ip) {
                  scene->addParticle(xp[ip],delta[ip]*dp_mag,data[ip]);
                }
              }
              else {
                WUI(WARN,"particle var " << particle_var << " does not eval to scalar: rendering particles without data");
                for (int ip = 0; ip < np; ++ip) {
                  scene->addParticle(xp[ip],delta[ip]*dp_mag);
                }
              }
            }
            else {
              WUI(WARN,"particle var " << particle_var << " is not associated with a length scale; using particle scale: " << dp_mag);
              if (particle_data->getType() == DN_DATA) {
                double * data = particle_data->getDNptr();
                for (int ip = 0; ip < np; ++ip) {
                  scene->addParticle(xp[ip],dp_mag,data[ip]);
                }
              }
              else {
                WUI(WARN,"particle var " << particle_var << " does not eval to scalar: rendering particles without data");
                for (int ip = 0; ip < np; ++ip) {
                  scene->addParticle(xp[ip],dp_mag);
                }
              }
            }
          }
          else {
            WUI(WARN,"particle var " << particle_var << " is not associated with an x position; skipping");
          }
        }
        else {
          // particle_data == NULL
          WUI(WARN,"particle var " << particle_var << " cannot be found; skipping");
        }
      }
      else {

        // old particle viz...
        // remove this as soon as everything works above...

        if (particle_var == "*") {
          // wild card = all particle classes
          for (int lp_index = 0,lim = lpHelperVec.size(); lp_index < lim; ++lp_index) {
            CtiRegister::CtiData * xp_data = CtiRegister::getRegisteredCtiData(lpHelperVec[lp_index].name+":xp",false);
            assert((xp_data)&&(xp_data->getType() == DN3_DATA)&&(xp_data->getUnindexedTopology() == LP_DATA));
            CtiRegister::CtiData * dp_data = CtiRegister::getRegisteredCtiData(lpHelperVec[lp_index].name+":dp",false);
            assert((dp_data)&&(dp_data->getType() == DN_DATA)&&(dp_data->getUnindexedTopology() == LP_DATA));
            const int nlp = xp_data->size();
            for (int ilp = 0; ilp < nlp; ++ilp) {
              const double xp[3] = {xp_data->dn3(ilp,0),xp_data->dn3(ilp,1),xp_data->dn3(ilp,2)};
              scene->addParticle(xp,dp_data->dn(ilp)*dp_mag);
            }
          }
        }
        else {
          map<const string,int>::iterator iter = lpHelperNameMap.find(particle_var);
          if (iter != lpHelperNameMap.end()) {
            // find the associated helper...
            int lp_index = iter->second;
            CtiRegister::CtiData * xp_data = CtiRegister::getRegisteredCtiData(lpHelperVec[lp_index].name+":xp",false);
            assert((xp_data)&&(xp_data->getType() == DN3_DATA)&&(xp_data->getUnindexedTopology() == LP_DATA));
            CtiRegister::CtiData * dp_data = CtiRegister::getRegisteredCtiData(lpHelperVec[lp_index].name+":dp",false);
            assert((dp_data)&&(dp_data->getType() == DN_DATA)&&(dp_data->getUnindexedTopology() == LP_DATA));
            const int nlp = xp_data->size();
            for (int ilp = 0; ilp < nlp; ++ilp) {
              const double xp[3] = {xp_data->dn3(ilp,0),xp_data->dn3(ilp,1),xp_data->dn3(ilp,2)};
              scene->addParticle(xp,dp_data->dn(ilp)*dp_mag);
            }
          }
          else {
            size_t found = particle_var.find("*:");
            // CI: what is this?
            if (found != string::npos) {
              // wild card = all particle classes
              for (int lp_index = 0,lim = lpHelperVec.size(); lp_index < lim; ++lp_index) {
                string this_str = MiscUtils::replaceAll(particle_var,"*:",lpHelperVec[lp_index].name+":");
                if (CtiRegister::CtiData * data = CtiRegister::getCtiData(this_str)) {
                  CtiRegister::CtiData * xp_data = CtiRegister::getRegisteredCtiData(lpHelperVec[lp_index].name+":xp",false);
                  assert((xp_data)&&(xp_data->getType() == DN3_DATA)&&(xp_data->getUnindexedTopology() == LP_DATA));
                  CtiRegister::CtiData * dp_data = CtiRegister::getRegisteredCtiData(lpHelperVec[lp_index].name+":dp",false);
                  assert((dp_data)&&(dp_data->getType() == DN_DATA)&&(dp_data->getUnindexedTopology() == LP_DATA));
                  const int nlp = data->size();
                  for (int ilp = 0; ilp < nlp; ++ilp) {
                    const double xp[3] = {xp_data->dn3(ilp,0),xp_data->dn3(ilp,1),xp_data->dn3(ilp,2)};
                    scene->addParticle(xp,dp_data->dn(ilp)*dp_mag,data->dn(ilp));
                  }
                }
              }
            }
            else if (CtiRegister::CtiData * data = CtiRegister::getCtiData(particle_var)) {

              if ((data->getType() != DN_DATA)||(data->getUnindexedTopology() != LP_DATA)) {
                CWARN(" > particle data variable " << particle_var << " must be DN LP DATA");
              }
              else {

                // to plot as particles, data should have X and possibly L data registered...
                double (*xp)[3] = NULL;
                if (data->getX(xp)) {
                  assert(xp != NULL);
                  // look for the PARTICLE_SIZE override...
                  double dp;
                  if (scene->getParticleSize(dp)) {
                    for (int ilp = 0, nlp = data->size(); ilp < nlp; ++ilp) {
                      scene->addParticle(xp[ilp],dp*dp_mag,data->dn(ilp));
                    }
                  }
                  else {
                    // otherwise there should be a lengthscale for the particles...
                    double *dp = NULL;
                    if (data->getL(dp)) {
                      assert(dp != NULL);
                      for (int ilp = 0, nlp = data->size(); ilp < nlp; ++ilp) {
                        scene->addParticle(xp[ilp],dp[ilp]*dp_mag,data->dn(ilp));
                      }
                    }
                    else {
                      // in this case, we should default to a certain size based on the image resolution: 4 pixels or something
                      CWARN(" > particle data variable " << particle_var << " does not have associated size, or use PARTICLE_SIZE <dp> for constant size");
                    }
                  }
                }
                else {
                  //CWARN(" > particle data variable " << particle_var << " does not have associated position");
                  // find the associated helper...
                  int lp_index = data->getIndex();
                  CtiRegister::CtiData * xp_data = CtiRegister::getRegisteredCtiData(lpHelperVec[lp_index].name+":xp",false);
                  assert((xp_data)&&(xp_data->getType() == DN3_DATA)&&(xp_data->getUnindexedTopology() == LP_DATA));
                  CtiRegister::CtiData * dp_data = CtiRegister::getRegisteredCtiData(lpHelperVec[lp_index].name+":dp",false);
                  assert((dp_data)&&(dp_data->getType() == DN_DATA)&&(dp_data->getUnindexedTopology() == LP_DATA));
                  const int nlp = data->size();
                  for (int ilp = 0; ilp < nlp; ++ilp) {
		    const double xp[3] = {xp_data->dn3(ilp,0),xp_data->dn3(ilp,1),xp_data->dn3(ilp,2)};
		    scene->addParticle(xp,dp_data->dn(ilp)*dp_mag,data->dn(ilp));
                  }
		}
              }
            }
          }
        }

      }
    }

    // -----------
    // solver vis
    // -----------

    // solver defined triangulated surfaces; e.g., PLIC (or FWH, flame-front, etc...)
    /*
      if (scene->hasSolverSurface()) {
      vector<SimpleTri> triVec;
      buildSolverSurface(triVec,scene->getSolverSurface());
      scene->addSimpleInternalTrisWithMeshOnFirstEdge(triVec);
      }
    */

    DELETE(cv_flag_geom_plane);
    DELETE(tet_flag_geom_plane);

    if ( scene->isClientRequest() ) {
      // don't add step index to filename when interacting with Cascade client...
      scene->writeImage();
    }
    else if ( killfile_request ) {
      // this image request came from a killfile command (outside of app)...
      scene->writeImage(step);
    }
    else {
      // buffer image
      scene->cacheImage();
    }

  }

  void processWriteImage(Param * param,const int step = 0,const double time = 0.0,const bool killfile_request = false) {

    // TODO: streamline this a little more later...

    // a request from the killfile that has an INTERVAL specified should be added to the
    // params rather than rendered right now...
    if (killfile_request) {
      int iarg = 0;
      while (iarg < param->size()) {
        string token = param->getString(iarg++);
        if (token == "INTERVAL") {
          // if we got an interval, add this...
          if (mpi_rank == 0) cout << " > adding WRITE_IMAGE to current run..." << endl;
          addParam(param);
          return;
        }
      }
    }

    WriteImageData wid;
    if (wid.init(param,step) != 0) // returns -1 if interval%step != 0, or -2 if parsing error caught
      return;

    // when interacting with the client, sometimes an image is not required,
    // even though one was requested.
    if (b_skip_image && wid.b_name && (skip_image_name == wid.name)) {
      // skip this immage...
      if (mpi_rank == 0) cout << " > skipping this image" << endl;
      b_skip_image = false;
      return;
    }

    // finally, if this is a killfile request, AND the image interval has been

    CtiScene * this_scene = new CtiScene(wid);
    if ( !this_scene->hasName() ) {

      CWARN( " > Warning: image NAME is missing, skipping ... \n \
            WRITE_IMAGE NAME <prefix> INTERVAL <interval> GEOM <geom> VAR <var-list>");
      delete this_scene;

    }
    else if ( (this_scene->getInterval() < 0) && !killfile_request ) {

      // interval required when it isn't a killfile request
      CWARN( " > Warning: invalid or missing INTERVAL, skipping.... \n \
            WRITE_IMAGE NAME <prefix> INTERVAL <interval> GEOM <geom> VAR <var-list>");
      delete this_scene;

    }
    else {
      doImage(this_scene,step,time,killfile_request);

      // if we cached an image store a pointer to the scene in a list so we can flush it later
      if (this_scene->isCached()) {
        sceneList.push_back(this_scene);
      }
      else {
        delete this_scene;
      }
    }
  }

  // TODO put these somewhere else eventually...

#include "tetCutPlane.hpp"
#include "tetCutPlaneData.hpp"
#include "tetCutPlaneWgts.hpp"

  void buildChtIso(vector<pair<SimpleTri,int> >& triVec,const CtiData * iso_var_no,const double iso_var_value);
  void buildChtIsoWithData(vector<pair<SimpleTriWithData,int> >& triVec,const CtiData * iso_var_no,const double iso_var_value,const CtiData * data_var);

  void buildIso(vector<pair<SimpleTri,int> >& triVec,const string& iso_var_name,const double iso_var_value,int * cv_hide_flag = NULL);
  void buildIso(vector<pair<SimpleTri,int> >& triVec,const string& iso_var_name,const CtiRegister::CtiData* iso_var,const double iso_var_value,int * cv_hide_flag = NULL);
  void buildIso(vector<pair<SimpleTri,int> >& triVec,const double * iso_var_no,const double iso_var_value,int * cv_hide_flag = NULL);
  // this version returns the containing cv as the second int...
  void buildIsoWithIcv(vector<pair<SimpleTri,int> >& triVec,const double * iso_var_no,const double iso_var_value);

  void buildIsoWithData(vector<pair<SimpleTriWithData,int> >& triVec,const string& iso_var_name,const double iso_var_value,const string& data_name,int * cv_hide_flag = NULL);

  // pass ctidata's if we already have them
  void buildIsoWithData(vector<pair<SimpleTriWithData,int> >& triVec,const string& iso_var_name,const CtiRegister::CtiData* iso_var,const double iso_var_value,const string& data_name, const CtiRegister::CtiData* data,int * cv_hide_flag = NULL);

  void buildIsoWithData(vector<pair<SimpleTriWithData,int> >& triVec,const double * const iso_var_no,const double iso_var_value,const double * const data_no,int * cv_hide_flag = NULL);

  void buildIsoWithWgtsAndData(vector<SimpleTriWithWgts>& triVec,const double * const iso_var_no,const double iso_var_value,double *&dn_tmp, double *&dn_no, const int nscalars, double (*&dn3_tmp)[3], double (*&dn3_no)[3], const int nvectors,int * cv_hide_flag = NULL);

  void setBfDN(double * no_data,const string& name,double *bf_data = NULL);
  void setBfDN(double * no_data,const string& name,const CtiRegister::CtiData * var,double *bf_data = NULL);

  void setBfDN(double * no_data,const vector<string>& name_vec,double * bf_data = NULL);
  void setBfDN(double * no_data,const vector<string>& name_vec,const vector<CtiRegister::CtiData *>& var_vec,double *bf_data = NULL);

  void setBfDN3(double (*no_data)[3],const string& name,double (*bf_data)[3] = NULL);
  void setBfDN3(double (*no_data)[3],const string& name,const CtiRegister::CtiData * var,double (*bf_data)[3] = NULL);

  void setBfDN3(double (*no_data)[3],const vector<string>& name_vec,double (*bf_data)[3] = NULL);
  void setBfDN3(double (*no_data)[3],const vector<string>& name_vec,const vector<CtiRegister::CtiData *>& var_vec,double (*bf_data)[3] = NULL);

  void setNoDN(double * no_data,const string& name);
  void setNoDN(double * no_data,const string& name,const CtiRegister::CtiData * var);

  void setNoDN3(double (*no_data)[3],const string& name);

  void setCvFlagFromIsoVar(int * cv_flag,const double * iso_var_no,const double iso_var_value);
  int setFlagsToIndexTmpData(int * cv_flag,int * bf_flag,int * fa_flag);
  void flagCvsAdjacentToZones(int8* cv_flag,bool* bfzone_flag,const int index);
  void flagCvsUnderIso(int8* cv_flag, const double *sd_no,const int index);

  // ==========================================================
  // these routines write data in tecplot, vtk, ensight...
  // see StaticSolver_write_data.cpp
  // ==========================================================

  // main routine...
  void buildDelaunayDual();
  void writeData(DataWriter* dw,const int step,const double time);
  void flagCvsWithOpenFaces(const bool (*hnofa_b)[2], const bool *hnobf_b, const int8* cv_flag);
  void flagCvsWithOpenFaces(const int8* cv_flag);
  void writeDataAscii(DataWriter* dw,const int step,const double time);
  int packCellCenteredDbuf(double* dbuf, const double* cv_data, const int8* cv_flag);
  int packCellCenteredDbuf(double* dbuf, const CtiRegister::CtiData* cv_data, const int8* cv_flag);
  int packTecplotDbuf(double* dbuf, const double* no_data, const int8* no_flag, const int8* bf_flag,
                      const int8* fa_flag, const int8* cv_flag);
  void writeFlaggedCvsTecplot(const string& filename,int8* cv_flag,const vector<string>& var_vec,const int nvar); // CI: WHY IS THIS int8?
  void writeFlaggedCvsFluentBoundaryProfile(const string& filename,int8* cv_flag,const vector<string>& var_vec,const int nvar);
  void writeFlaggedBfZonesTecplot(const string& filename,bool* bfzone_flag,const vector<string>& var_vec_cv,const int nvar_cv,
                                  const vector<string>& var_vec_bf,const vector<string>& zones_vec_bf,const int nvar_bf,const int nzones_bf);
  void writeFlaggedFacesASCII(const string& filename,const int * const fa_flag) const;
  // this version is the locally compressed binary
  void writeSimpleTriVecTecplot(const string& filename,const vector<SimpleTriWithWgts>& triVec,
                                const int nst, const int (*spost)[3], const int nsp, const vector<string>& var_vec,const double *data_no,
                                const double *data_tmp, const int nvar, const int nno);
  // this version is the uncompressed ascii
  void writeSimpleTriVecTecplotASCII(const string& filename,const vector<SimpleTriWithWgts>& triVec,
                                     const vector<string>& var_vec,const double *data_no, const double *data_tmp, const int nvar, const int nst);
  void writeSimpleTriVecVTK(const string& filename,const vector<SimpleTriWithWgts>& triVec,
                            const int nst,const int (*spost)[3],const int nsp, const vector<string>& scalar_names,
                            const double *dn_no,const double *dn_tmp,const int nscalars, const vector<string>& vector_names,
                            const double (*dn3_no)[3],const double (*dn3_tmp)[3],const int nvectors,const int nno);
  // multi-file
  void writeFlaggedCvsParallelVTK(const string& master_filename,const string& filename,int8* cv_flag,
                                  const vector<string>& scalar_names,const int nscalars,const vector<string>& vector_names,const int nvectors);
  // single-file
  void writeFlaggedCvsVTK(const string& filename,int8* cv_flag,const vector<string>& scalar_names,
                          const int nscalars, const vector<string>& vector_names,const int nvectors);
  void writeFlaggedBfZonesVTK(const string& filename,bool* bfzone_flag,
                              const vector<string>& scalar_names_cv,const int nscalars_cv,
                              const vector<string>& vector_names_cv,const int nvectors_cv,
                              const vector<string>& scalar_names_bf,const vector<string>& scalar_zones_bf,const int nscalars_bf,
                              const vector<string>& vector_names_bf,const vector<string>& vector_zones_bf,const int nvectors_bf,const int nzones_bf);
  void updateEnsightCaseFile(DataWriter* dw,const int step,const double time);
  void writeFlaggedCvsEnsight(DataWriter* dw,int8* cv_flag,const int step,const double time);
  //void flagHangingNodes(bool * b_hno_flag);
  void flagHangingNodes(bool (*hnofa_b)[2], bool * hnobf_b, const int8 *cv_flag);
  //void buildNodeEdgeCountMap(const int icv, const bool * b_hno_flag, map<const int,int> &nodeEdgeCountMap);
  void buildNodeEdgeCountMap(const int icv, const bool (*hnofa_b)[2], const bool * hnobf_b, map<const int,int> &nodeEdgeCountMap);
  void buildNodeEdgeCountMap(const int icv, map<const int,int> &nodeEdgeCountMap);
  void writeFlaggedCvsDualEnsightMultiPart(DataWriter* dw,const int step,const double time);
  void writeFlaggedCvsEnsightMultiPart(DataWriter* dw,int8* cv_flag,const int step,const double time);
  void writeCvsAndFlaggedBfZonesEnsight(DataWriter* dw,const int step,const double time);
  void writeCvsAndFlaggedBfZonesDualEnsightMultiPart(DataWriter* dw,const int step,const double time);
  void writeCvsAndFlaggedBfZonesEnsightMultiPart(DataWriter* dw,const int step,const double time);
  void writeFlaggedBfZonesDualEnsight(DataWriter* dw,const int step,const double time);
  void writeFlaggedBfZonesEnsight(DataWriter* dw,const int step,const double time);
  void writeSimpleTriVecEnsight(DataWriter* dw,const vector<SimpleTriWithWgts>& triVec,
                                const int nst,const int (*spost)[3],const int nsp,const double *dn_no,
                                const double *dn_tmp,const double (*dn3_no)[3],const double (*dn3_tmp)[3],
                                const int nno,const int step,const double time);
  void writeFlaggedLpsEnsight(DataWriter* dw,const int step,const double time);
  void writeSimpleTriVecStl(DataWriter* dw,const vector<SimpleTriWithWgts>& triVec,const int step,const double time);

  // main calling routine for writing data to tecplot, vtk etc...
  void processWriteData(Param * param,const int step,const double time,const bool killfile_request = false) {
    if (param->size() == 0) {
      CWARN(" > Warning: write_data syntax should be: \n \
            WRITE_DATA INTERVAL <interval> NAME <prefix> FORMAT <format> GEOM <geometry> VARS <var-list>");
    }
    else {
      string name;
      int interval = -1;
      int format = TECPLOT_FORMAT; // default to tecplot
      bool b_dual = false; // default to voronoi poly
      bool b_inside = false; // change ISO_TOPO to CV_TOPO
      int topo = -1;
      int geom = -1;
      SimpleGeom* simple_geom = NULL;
      double iso_var_value = HUGE_VAL; // needs to be set if used
      string iso_var_name; // needed for building user-defined iso-surfaces using registration

      vector<string> bfzone_names; // used for BF_GEOM
      set<pair<int,string> > variable_names; // use a set to protect user against writing the same thing twice
      map<const string,DataWriter*>::iterator dw_it = dwMap.end();
      int nzones_bf = -1;
      string sbin_filename = "";
      int lp_index = -1;

      int iarg = 0;
      while (iarg < param->size()) {
        string token = param->getString(iarg++);
        if (token == "NAME") {
          name = param->getString(iarg++);
          // check to see if we already have a data writer with this name.
          // note that if the user uses the same name for a different
          // WRITE_DATA command, the latter will skipped...
          dw_it = dwMap.find(name);
          if ( dw_it != dwMap.end())
            break;
        }
        else if (token == "FORMAT") {
          token = param->getString(iarg++);
          if (token == "TECPLOT")
            format = TECPLOT_FORMAT;
          else if (token == "VTK")
            format = VTK_FORMAT;
          else if (token == "PVTK")
            format = PVTK_FORMAT;
          else if (token == "ENSIGHT")
            format = ENSIGHT_FORMAT;
          else if (token == "FLUENT_BP")
            format = FLUENT_BP_FORMAT;
          else if (token == "ASCII")
            format = ASCII_FORMAT;
          else if (token == "STL")
            format = STL_FORMAT;
          else {
            CWARN("unrecognized FORMAT token in processWriteData: " << token);
            return;
          }
        }
        else if (token == "INTERVAL") {
          interval = param->getInt(iarg++);
          if (interval <= 0) {
            CWARN(" > INTERVAL expects a positive integer; treating as one-off");
            // invalid value entered, so treat as unset (default value)
            interval = -1;
          }
        }
        else if (token == "DUAL") {
          b_dual = true;
        }
        else if ((token == "IN")||(token == "INSIDE")||(token == "UNDER")) {
          b_inside = true;
        }
        else if (token == "GEOM") {
          token = param->getString(iarg++);
          if (token == "ISO") {
            topo = ISO_TOPO;
            geom = USERDEF_GEOM;
            iso_var_name = param->getString(iarg++);
            iso_var_value = param->getDouble(iarg++);
          }
          else if (token == "LP") {
            topo = LP_TOPO;
            geom = ALL_GEOM; // TODO
          }
          else if ((token == "ZONE")||(token == "ZONES")||(token == "FAZONE")||(token == "BFZONE")) {
            topo = BF_TOPO;
            token = param->getString(iarg);
            while (token != "VARS") {
              ++iarg;
              bfzone_names.push_back(token);

              // get the next zone name...
              if (iarg == param->size()) break;
              token = param->getString(iarg);
            }
            nzones_bf = bfzone_names.size();
          }
          else if (token == "ALL") {
            topo = CV_TOPO;
            geom = ALL_GEOM;
          }
          else if (token == "ALL_WITH_ZONES") {
            // ONLY FOR ENSIGHT!
            topo = CV_TOPO; // it uses cv data to project on boundaries, so it is like a cv topo
            geom = ALL_WITH_ZONES_GEOM;

            token = param->getString(iarg);
            // assume it is all zones if not included (VARS is next token)...
            if (token == "VARS") {
              FOR_IZONE(bfZoneVec)
                bfzone_names.push_back(bfZoneVec[izone].getName());
            }
            else {
              while (token != "VARS") {
                ++iarg;
                bfzone_names.push_back(token);

                // get the next zone name...
                if (iarg == param->size()) break;
                token = param->getString(iarg);
              }
            }
            nzones_bf = bfzone_names.size();
          }
          else if (token == "SBIN") {
            topo = ISO_TOPO;
            geom = FILE_GEOM;
            sbin_filename = param->getString(iarg++);
          }
          else {
            topo = ISO_TOPO;
            geom = SIMPLE_GEOM;
            --iarg;
            assert(simple_geom == NULL);
            simple_geom = newSimpleGeom(param,iarg);
            if (simple_geom == NULL) {
              CWARN("unrecognized GEOM token in processWriteData: " << token);
              return;
            }
          }
        }
        else if (token == "VARS") {
          while (iarg < param->size()) {
            string str = param->getString(iarg++);
            // fazone variables should begin with *:
            if (str.find("*:") != string::npos) {
              if (nzones_bf == -1) {
                CWARN(" > Warning: must specify FAZONE's before VARS, returning.");
                return;
              }
              bool zone_missing_var = false;
              for (int i = 0; i < nzones_bf; ++i) {
                string this_str = MiscUtils::replaceAll(str,"*:",bfzone_names[i]+":");
                CtiRegister::CtiData * data = CtiRegister::getCtiData(this_str);
                if (data == NULL) {
                  zone_missing_var = true;
                  break;
                }
              }
              // Only add variable if it exists across all zones...
              if (!zone_missing_var) {
                for (int i = 0; i < nzones_bf; ++i) {
                  string this_str = MiscUtils::replaceAll(str,"*:",bfzone_names[i]+":");
                  variable_names.insert(pair<int,string>(i,this_str)); // use index into zone names to sort
                }
              }
              else {
                CWARN(" > Warning: variable " << str << " DNE in all zones, skipping...");
              }
            }
            else if ((topo == BF_TOPO)&&(str.find(":") != string::npos)) {
              CWARN(" > Warning: please use wildcard \"*\" instead of zone name before \":\" when specifying FAZONE VARS, skipping...");
            }
            else {
              if (str == "STATS") {
                Param * stats_param = getParam("STATS");
                if (stats_param != NULL) {
                  for (int iarg2 = 0; iarg2 < stats_param->size(); ++iarg2) {
                    const string vname = stats_param->getString(iarg2);

                    vector<string> var_names;
                    if (vname.find("*:") != string::npos) {
                      for (int izn=0, nzn=bfZoneVec.size(); izn<nzn; ++izn) {
                        string data_name = vname;
                        data_name = MiscUtils::replaceAll(data_name,"*:",(bfZoneVec[izn].getName()+":"));
                        var_names.push_back(data_name);
                      }
                    }
                    else {
                      // if not a wild-card string, simply add this variable
                      var_names.push_back(vname);
                    }

                    bool found = false;
                    for (int v=0,vmax=var_names.size(); v<vmax; ++v) {
                      CtiRegister::CtiData * data = CtiRegister::getCtiData(var_names[v]);
                      if (data != NULL) {
                        found = true;
                        if ((data->getUnindexedTopology() == BF_DATA)&&((topo == ISO_TOPO)||(topo == CV_TOPO)))
                          continue; // bf geom supports both cv and bf data; however, cv and iso geoms don't support bf data
                        if (data->getType() == D_DATA || data->getType() == DN_DATA || data->getType() == DN3_DATA) {
                          variable_names.insert(pair<int,string>(-1,vname+"_avg"));
                          variable_names.insert(pair<int,string>(-1,vname+"_rms"));
                          if (data->getType() == DN3_DATA)
                            variable_names.insert(pair<int,string>(-1,vname+"_rey"));
                        }
                      }
                    }

                    if (!found) {
                      CWARN(" > Warning: cannot evaluate \"" << vname << "\" statistics for tecplot data. Skipping.\n");
                    }
                  }
                }
                else {
                  CWARN(" > No STATS command found in the input file, so Tecplot cannot write statistics!");
                }
              }
              else {
                CtiRegister::CtiData * data = CtiRegister::getCtiData(str);
                if (data != NULL) {
                  bool b_insert = true;
                  const int ut = data->getUnindexedTopology();
                  if (ut == LP_DATA) {
                    if (lp_index == -1) {
                      lp_index = data->getIndex();
                    }
                    else if (lp_index != data->getIndex()) {
                      CWARN(" > Skipping mismatched LP class: " << lp_index << " (current) " << data->getIndex() << " (requested), skipping...");
                      b_insert = false;
                    }
                  }
                  else if (ut == BF_DATA) {
                    // if we have bf data and topo is still unset, this is probably fine:
                    if (topo == -1)
                      topo = BF_TOPO_EXPLICIT; // BF_TOPO_EXPLICIT means the user has explicitly specified the boundary data, e.g. x0:area
                  }
                  if (b_insert)
                    variable_names.insert(pair<int,string>(-1,str));
                }
                else {
                  CWARN(" > Warning: unrecognized variable \"" << str << "\", skipping...");
                }
              }
            }
          }
        }
        else {
          CWARN(" > unrecognized token in processWriteData: " << token);
          return;
        }
      }

      if ( dw_it == dwMap.end() ) {

        if (b_inside&&(topo == ISO_TOPO))
          topo = CV_TOPO;

        // this data_writer has not been parsed before and so it needs to be
        // constructed and the inputs need to be checked for errors..

        if (name == "" ) {

          CWARN( " > Warning: write_data NAME is missing; \n \
                WRITE_DATA INTERVAL <interval> NAME <prefix> FORMAT <format> GEOM <geometry> VARS <var-list>");

        }
        else if ( interval < 0 && !killfile_request ) {

          CWARN( " > Warning: invalid or missing INTERVAL; \n \
                WRITE_DATA INTERVAL <interval> NAME <prefix> FORMAT <format> GEOM <geometry> VARS <var-list>");

        }
        else if ( (topo < 0) || ((topo == ISO_TOPO) && (geom < 0)) ) {

          CWARN( " > Warning: invalid or missing GEOM; \n \
                WRITE_DATA INTERVAL <interval> NAME <prefix> FORMAT <format> GEOM <geometry> VARS <var-list>");

        }
        else {

          // expand components of vectors (and eventually tensors)
          // tecplot is written as scalars
          // VTK is written as combination
          vector<string> scalar_names_cv, vector_names_cv;
          vector<string> scalar_names_bf, vector_names_bf;
          vector<string> scalar_zones_bf, vector_zones_bf;
          vector<string> scalar_names_lp, vector_names_lp;
          const int nzones_bf = bfzone_names.size();
          for (set<pair<int,string> >::iterator it = variable_names.begin(); it != variable_names.end(); ++it) {
            CtiRegister::CtiData * data = CtiRegister::getCtiData(it->second);
            assert(data != NULL); // we should of ensure data exists prior to here. if not, logic bug.
            const int data_type = data->getType();
            const int data_topo = data->getUnindexedTopology();
            // COUT1(it->second << " identified as type, topo: " << data->getType() << "," << data->getUnindexedTopology());
            if ((format == TECPLOT_FORMAT) || (format == FLUENT_BP_FORMAT) || (format == ASCII_FORMAT)) {
              if (data_type == DN_DATA || data_type == IN_DATA || data_type == D_DATA || data_type == I_DATA) {
                if (data_topo == CV_DATA) {
                  scalar_names_cv.push_back(it->second);
                }
                else if (data_topo == BF_DATA) {
		  scalar_names_bf.push_back(it->second);
                  scalar_zones_bf.push_back(data->getName());
                }
                else if (data_topo == LP_DATA) {
                  scalar_names_lp.push_back(it->second);
                }
              }
              else if (data_type == DN3_DATA || data_type == D3_DATA) {
                if (data_topo == CV_DATA) {
                  FOR_I3 {
                    stringstream ss;
                    ss << i;
                    scalar_names_cv.push_back("comp("+(it->second)+","+ss.str()+")");
                  }
                }
                else if (data_topo == BF_DATA) {
		  FOR_I3 {
		    stringstream ss;
		    ss << i;
		    scalar_names_bf.push_back("comp("+(it->second)+","+ss.str()+")");
                    scalar_zones_bf.push_back(data->getName());
                  }
                }
                else if (data_topo == LP_DATA) {
                  FOR_I3 {
                    stringstream ss;
                    ss << i;
                    scalar_names_lp.push_back("comp("+(it->second)+","+ss.str()+")");
                  }
                }
              }
            }
            else if (format == VTK_FORMAT || format == PVTK_FORMAT || format == ENSIGHT_FORMAT) {
              if (data_type == DN_DATA || data_type == IN_DATA || data_type == D_DATA || data_type == I_DATA) {
                if (data_topo == CV_DATA) {
                  scalar_names_cv.push_back(it->second);
                }
                else if (data_topo == BF_DATA) {
		  scalar_names_bf.push_back(it->second);
                  scalar_zones_bf.push_back(data->getName());
                }
                else if (data_topo == LP_DATA) {
                  scalar_names_lp.push_back(it->second);
                }
              }
              else if (data_type == DN3_DATA || data_type == D3_DATA) {
                if (data_topo == CV_DATA) {
                  vector_names_cv.push_back(it->second);
                }
                else if (data_topo == BF_DATA) {
		  vector_names_bf.push_back(it->second);
                  vector_zones_bf.push_back(data->getName());
                }
                else if (data_topo == LP_DATA) {
                  vector_names_lp.push_back(it->second);
                }
              }
            }
          }

          // put here to skip when not on interval
          const int nscalars_cv = scalar_names_cv.size();
          const int nvectors_cv = vector_names_cv.size();
          const int nscalars_bf = scalar_names_bf.size();
          const int nvectors_bf = vector_names_bf.size();
          const int nscalars_lp = scalar_names_lp.size();
          const int nvectors_lp = vector_names_lp.size();

          if ((nscalars_cv + nvectors_cv + nscalars_bf + nvectors_bf + nscalars_lp + nvectors_lp) == 0) {
            CWARN(" > Warning: there are no acceptable variables.");
          }

          if (topo == LP_DATA) {
            bool b_err = false;
            if (lp_index == -1) {
              CWARN(" > Warning: there is no indication of particle class from variables, skipping...");
              b_err = true;
            }
            if (format != ENSIGHT_FORMAT) {
              CWARN(" > Warning: ENSIGHT is the only supported format for LP data, skipping...");
              b_err = true;
            }
            if (b_err)
              return;
          }

          // this is valid write data entry and needs to be added to the map
          pair<map<const string,DataWriter*>::iterator, bool> ret
            = dwMap.insert(pair<const string,DataWriter*>(name,new DataWriter()));

          // we check for membership earlier, so the insertion must have succeeded.
          assert( ret.second == true);

          dw_it = ret.first;
          dw_it->second->param_str = param->str();
          dw_it->second->name = name;
          dw_it->second->interval = interval;
          dw_it->second->format = format;
          dw_it->second->topo = topo;
          dw_it->second->geom = geom;
          dw_it->second->simple_geom = simple_geom;
          dw_it->second->b_dual = b_dual;
          dw_it->second->iso_var_name = iso_var_name;
          dw_it->second->iso_var_value = iso_var_value;
          dw_it->second->iso_var = NULL; // allocated during first population
          dw_it->second->cv_flag = NULL;
          dw_it->second->lp_flag = NULL;
          dw_it->second->scalar_names_cv = scalar_names_cv;
          dw_it->second->vector_names_cv = vector_names_cv;
          assert(scalar_names_bf.size() == scalar_zones_bf.size());
          assert(vector_names_bf.size() == vector_zones_bf.size());
          dw_it->second->scalar_names_bf = scalar_names_bf;
          dw_it->second->vector_names_bf = vector_names_bf;
          dw_it->second->scalar_zones_bf = scalar_zones_bf;
          dw_it->second->vector_zones_bf = vector_zones_bf;
          dw_it->second->scalar_names_lp = scalar_names_lp;
          dw_it->second->vector_names_lp = vector_names_lp;
          dw_it->second->nscalars_cv = nscalars_cv;
          dw_it->second->nvectors_cv = nvectors_cv;
          dw_it->second->nscalars_bf = nscalars_bf;
          dw_it->second->nvectors_bf = nvectors_bf;
          dw_it->second->nscalars_lp = nscalars_lp;
          dw_it->second->nvectors_lp = nvectors_lp;
          dw_it->second->nzones_bf = nzones_bf;
          const int nzones = bfZoneVec.size();
          dw_it->second->bfzone_flag = new bool[nzones];
          for (int izone = 0; izone < nzones; ++izone)
            dw_it->second->bfzone_flag[izone] = false;
          for (int i = 0; i < nzones_bf; ++i) {
            map<const string,int>::iterator bf_it = bfZoneNameMap.find(bfzone_names[i]);
            if (bf_it != bfZoneNameMap.end())
              dw_it->second->bfzone_flag[bf_it->second] = true;
          }
          dw_it->second->sbin_filename = sbin_filename;
          if (lp_index != -1) {
            assert((lp_index >= 0)&&(lp_index < int(lpHelperVec.size())));
            dw_it->second->lp_index = lp_index;
          }
        }
      }

      if ( dw_it != dwMap.end() ) {
        if ( dw_it->second->interval == -1 && killfile_request ) {
          // assume it is a one-off if it is a killfile request without an interval
          writeData(dw_it->second,step,time);
          delete dw_it->second;
          dwMap.erase(dw_it);
        }
        else if ( step%dw_it->second->interval == 0) {
          writeData(dw_it->second,step,time);
        }
      }

    }
  }

private:

  void initBfZoneVec(StripedMesh * sm);

  void loadBalance(StripedMesh * sm);

  void redistReorderMesh(StripedMesh * sm);

  void buildPeriodicRt();

  void groupFaces() const {
    groupExtendedFaces();
    groupCompactFaces();
  }
  void groupExtendedFaces() const;
  void groupExtendedFacesOld() const;
  void groupCompactFaces() const;

public:

  void calcCvGrad(double (*__restrict__ dpdx)[3],const double *__restrict__ p);
  void calcCvGrad(double (*__restrict__ dudx)[3][3],const double (*__restrict__ u)[3]);
  void _calcCvGradRange(double (*__restrict__ dudx)[3][3],const double (*__restrict__ u)[3],
                        const int icv_start, const int icv_end);
  void calcCvGradStart(double (*dudx)[3][3], const double (*u)[3]) {
    _calcCvGradRange(dudx,u,0,ncv_i);
  }
  void calcCvGradFinish(double (*dudx)[3][3], const double (*u)[3]) {
    _calcCvGradRange(dudx,u,ncv_i,ncv);
  }
  void checkCvGradCoeff(int * cv_flag = NULL);
  void calcCvGradCoeff();
  void checkGradExtended();

  void buildCvPrcomm();
  void buildCv2Prcomm();
  void buildPrcommSymmetric(vector<CvPrcomm>& prcommVec,uint8 * rbi_g,const int ncv,const int ncv_g);
  void buildPrcomm(vector<CvPrcomm>& prcommVec,uint8 * rbi_g,const set<int>& rankSet,const int ncv,const int ncv_g);

  void buildXfa();

  void buildFaocv();
  void buildBfocv();
  void buildCvocv();
  void buildCvocv2();

  void buildNoPrcomm();
  void buildFaPrcomm();

  void computeBfDistanceFromFlaggedCvs(double * wall_dist, const BfZone * zone_ptr, const double * cv_flag) const;

public:

  bool checkDataFlag(const string& vname) const {
    return CtiRegister::checkDataFlag(vname);
  }

  void setDataFlag(const string& vname,const int val) {
    CtiRegister::setDataFlag(vname,val);
  }

  void registerData() {

    COUT1("StaticSolver::registerData()");

    registerData(mpi_rank,"mpi_rank",0); // rw_bits = 0

    // add so it is present in *.sles
    registerCvData(x_vv,"x_vv",WRITE_DATA|X_DATA);

    // these guys CAN be written to lets say a snapshot (they have a dde/xora)...
    registerCvData(r_vv,"r_vv",CAN_WRITE_DATA);
    registerCvData(x_cv,"x_cv",CAN_WRITE_DATA);
    registerCvData(vol_cv,"vol_cv",CAN_WRITE_DATA);
    registerCvData(inv_vol,"inv_vol",CAN_WRITE_DATA);

    registerFaData(group_fa,"group_fa",CAN_WRITE_DATA);
    //registerSignedFaData(n_fa,"n_fa",CAN_WRITE_DATA);
    registerFaData(n_fa,"n_fa",CAN_WRITE_DATA); // CI - mesh handles flipping, so don't sign face
    registerFaData(area_over_delta_fa,"area_over_delta_fa",CAN_WRITE_DATA);

    // NO_DATA and EF_DATA currently cannot be written...
    CtiRegister::_registerData(x_no,"x_no",NO_DATA,NO_READWRITE_DATA|X_DATA,nno);

    CtiRegister::_registerData(group_ef,"group_ef",EF_DATA,NO_READWRITE_DATA,nef);
    CtiRegister::_registerData(n_ef,"n_ef",EF_DATA,NO_READWRITE_DATA,nef);
    CtiRegister::_registerData(c_ef,"c_ef",EF_DATA,NO_READWRITE_DATA,nef);

    // register user defined transforms based on interpolation..

    FOR_PARAM_MATCHING("REGISTER_INTERPOLATE_TRANSFORM") {

      string filename = "";
      string func_name = "";

      int iarg = 0;
      while ( iarg < param->size()) {

        const string token = param->getString(iarg++);

        if ( token == "FILENAME") {
          filename = param->getString(iarg++);
        } else if ( token == "FUNC_NAME") {
          func_name = param->getString(iarg++);
        } else {
          CERR( " > unknown token : " << token);
        }

      }

      if ( (filename == "") || (func_name == "")) {
        CERR( " > incorrect syntax: REGISTER_INTERPOLATE_TRANSFORM FILENAME <fname> FUNC_NAME <func_name>");
      }

      interp_transforms.push_back(new CtiInterpolateTransform(func_name,filename));

    }


  }

  void registerData(int& val,const string& vname,const uint rw_bits) {
    CtiRegister::_registerData(val,vname,rw_bits);
  }
  void registerData(double& val,const string& vname,const uint rw_bits) {
    CtiRegister::_registerData(val,vname,rw_bits);
  }
  void registerData(double val[3],const string& vname,const uint rw_bits) {
    CtiRegister::_registerData(val,vname,rw_bits);
  }

  void registerI(const string& vname,const uint rw_bits) {
    CtiRegister::_registerData(vname,I_DATA,rw_bits);
  }


  void registerD(const string& vname,const uint rw_bits) {
    CtiRegister::_registerData(vname,D_DATA,rw_bits);
  }

  void registerD3(const string& vname,const uint rw_bits) {
    CtiRegister::_registerData(vname,D3_DATA,rw_bits);
  }

  void registerCvData(int *&val,const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(val,vname,CV_DATA,rw_bits,ncv,cvora_striped,dde_cv_striped);
    }
    else {
      CtiRegister::_registerData(val,vname,CV_DATA,rw_bits,ncv);
    }
  }

  void registerCvData(double *&val,const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(val,vname,CV_DATA,rw_bits,ncv,cvora_striped,dde_cv_striped);
    }
    else {
      CtiRegister::_registerData(val,vname,CV_DATA,rw_bits,ncv);
    }
  }

  void registerCvData(double (*&val)[3],const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(val,vname,CV_DATA,rw_bits,ncv,cvora_striped,dde_cv_striped);
    }
    else {
      CtiRegister::_registerData(val,vname,CV_DATA,rw_bits,ncv);
    }
  }

  template<class T>
  void registerCvData(T *& state,int &val,const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(state,val,vname,CV_DATA,rw_bits,ncv,cvora_striped,dde_cv_striped);
    }
    else {
      CtiRegister::_registerData(state,val,vname,CV_DATA,rw_bits,ncv);
    }
  }

  template<class T>
  void registerCvData(T *& state,double &val,const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(state,val,vname,CV_DATA,rw_bits,ncv,cvora_striped,dde_cv_striped);
    }
    else {
      CtiRegister::_registerData(state,val,vname,CV_DATA,rw_bits,ncv);
    }
  }

  template<class T>
  void registerCvData(T *& state,double (&val)[3],const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(state,val,vname,CV_DATA,rw_bits,ncv,cvora_striped,dde_cv_striped);
    }
    else {
      CtiRegister::_registerData(state,val,vname,CV_DATA,rw_bits,ncv);
    }
  }

  void registerFaData(int *&val,const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(val,vname,FA_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(val,vname,FA_DATA,rw_bits,nfa);
    }
  }

  void registerFaData(double *&val,const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(val,vname,FA_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(val,vname,FA_DATA,rw_bits,nfa);
    }
  }

  void registerFaData(double (*&val)[3],const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(val,vname,FA_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(val,vname,FA_DATA,rw_bits,nfa);
    }
  }

  template<class T>
  void registerFaData(T *& state,int &val,const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(state,val,vname,FA_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(state,val,vname,FA_DATA,rw_bits,nfa);
    }
  }

  template<class T>
  void registerFaData(T *& state,double &val,const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(state,val,vname,FA_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(state,val,vname,FA_DATA,rw_bits,nfa);
    }
  }

  template<class T>
  void registerFaData(T *& state,double (&val)[3],const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(state,val,vname,FA_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(state,val,vname,FA_DATA,rw_bits,nfa);
    }
  }

  void registerSignedFaData(int *&val,const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(val,vname,SIGNED_FA_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(val,vname,SIGNED_FA_DATA,rw_bits,nfa);
    }
  }

  void registerSignedFaData(double *&val,const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(val,vname,SIGNED_FA_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(val,vname,SIGNED_FA_DATA,rw_bits,nfa);
    }
  }

  void registerSignedFaData(double (*&val)[3],const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(val,vname,SIGNED_FA_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(val,vname,SIGNED_FA_DATA,rw_bits,nfa);
    }
  }

  template<class T>
  void registerSignedFaData(T *& state,int &val,const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(state,val,vname,SIGNED_FA_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(state,val,vname,SIGNED_FA_DATA,rw_bits,nfa);
    }
  }

  template<class T>
  void registerSignedFaData(T *& state,double &val,const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(state,val,vname,SIGNED_FA_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(state,val,vname,SIGNED_FA_DATA,rw_bits,nfa);
    }
  }

  template<class T>
  void registerSignedFaData(T *& state,double (&val)[3],const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(state,val,vname,SIGNED_FA_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(state,val,vname,SIGNED_FA_DATA,rw_bits,nfa);
    }
  }

  // -----------------------------------------------------------------------------------------------------
  // note that these guys are stored in registration space...
  // -----------------------------------------------------------------------------------------------------
  void registerCvIN(const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(vname,CV_DATA,IN_DATA,rw_bits,ncv,cvora_striped,dde_cv_striped);
    }
    else {
      CtiRegister::_registerData(vname,CV_DATA,IN_DATA,rw_bits,ncv);
    }
  }
  void registerCvDN(const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(vname,CV_DATA,DN_DATA,rw_bits,ncv,cvora_striped,dde_cv_striped);
    }
    else {
      CtiRegister::_registerData(vname,CV_DATA,DN_DATA,rw_bits,ncv);
    }
  }
  void registerCvDN3(const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(vname,CV_DATA,DN3_DATA,rw_bits,ncv,cvora_striped,dde_cv_striped);
    }
    else {
      CtiRegister::_registerData(vname,CV_DATA,DN3_DATA,rw_bits,ncv);
    }
  }
  void registerFaIN(const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(vname,FA_DATA,IN_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(vname,FA_DATA,IN_DATA,rw_bits,nfa);
    }
  }
  void registerFaDN(const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(vname,FA_DATA,DN_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(vname,FA_DATA,DN_DATA,rw_bits,nfa);
    }
  }
  void registerFaDN3(const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(vname,FA_DATA,DN3_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(vname,FA_DATA,DN3_DATA,rw_bits,nfa);
    }
  }
  void registerSignedFaIN(const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(vname,SIGNED_FA_DATA,IN_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(vname,SIGNED_FA_DATA,IN_DATA,rw_bits,nfa);
    }
  }
  void registerSignedFaDN(const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(vname,SIGNED_FA_DATA,DN_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(vname,SIGNED_FA_DATA,DN_DATA,rw_bits,nfa);
    }
  }
  void registerSignedFaDN3(const string& vname,const uint rw_bits) {
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)) {
      CtiRegister::_registerData(vname,SIGNED_FA_DATA,DN3_DATA,rw_bits,nfa,faora_striped,dde_fa_striped);
    }
    else {
      CtiRegister::_registerData(vname,SIGNED_FA_DATA,DN3_DATA,rw_bits,nfa);
    }
  }

  double* createCvD1Data(CtiRegister::CtiData& v) {

    v.new_dn(CV_DATA,ncv);
    v.setDdeStuff(cvora_striped,dde_cv_striped);
    v.setBits(CAN_WRITE_DATA);
    return v.getDNptr();

  }

  double (* createCvD2Data(CtiRegister::CtiData& v))[3]  {

    v.new_dn3(CV_DATA,ncv);
    v.setDdeStuff(cvora_striped,dde_cv_striped);
    v.setBits(CAN_WRITE_DATA);
    return v.getDN3ptr();

  }

  double (* createCvD3Data(CtiRegister::CtiData& v))[3][3] {

    v.new_dn33(CV_DATA,ncv);
    v.setDdeStuff(cvora_striped,dde_cv_striped);
    v.setBits(CAN_WRITE_DATA);
    return v.getDN33ptr();

  }

  double* createFaD1Data(CtiRegister::CtiData& v) {

    v.new_dn(FA_DATA,nfa);
    v.setDdeStuff(faora_striped,dde_fa_striped);
    v.setBits(CAN_WRITE_DATA);
    return v.getDNptr();

  }

  double (* createFaD2Data(CtiRegister::CtiData& v))[3]  {

    v.new_dn3(FA_DATA,nfa);
    v.setDdeStuff(faora_striped,dde_fa_striped);
    v.setBits(CAN_WRITE_DATA);
    return v.getDN3ptr();

  }

  double* createSignedFaD1Data(CtiRegister::CtiData& v) {

    v.new_dn(SIGNED_FA_DATA,nfa);
    v.setDdeStuff(faora_striped,dde_fa_striped);
    v.setBits(CAN_WRITE_DATA);
    return v.getDNptr();

  }

  double (* createSignedFaD2Data(CtiRegister::CtiData& v))[3]  {

    v.new_dn3(SIGNED_FA_DATA,nfa);
    v.setDdeStuff(faora_striped,dde_fa_striped);
    v.setBits(CAN_WRITE_DATA);
    return v.getDN3ptr();

  }

  // -----------------------------------------------------------------------------------------------------

  void writeLpocvAndInitDdeStuff(MPI_File &fh,MPI_Offset &offset,const int lp_index);
  void writeData(const string& filename);
  void writeResult(const int step) {

    char filename[128];
    sprintf(filename,"result.%08d.sles",step);
    writeData(filename);

  }
  void writeResult() {

    char filename[128];
    sprintf(filename,"result.sles");
    writeData(filename);

  }

public:

  template<class T, typename K>
  void updateCvClassStart(T* s, const int action = REPLACE_ROTATE_DATA) {

    //=======================
    // cv class level 1 only
    //=======================
    // T class must provide its own pack and unpack routines.. .
    const int UPDATE_TAG = 321123;

    map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
    if (iter != mpiRequestMap.end())
      CERR("updateCvClassStart(T*): s is already mapped. updateCvClassFinish required on s.");

    MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
    mpiRequestMap[(void*)s] = mrs;

    mrs->sendRequestArray = new MPI_Request[cvPrcommVec.size()];
    mrs->recvRequestArray = new MPI_Request[cvPrcommVec.size()];

    const int data_size = T::data_size();
    int unpack_size     = 0;
    int pack_size       = 0;
    for (int ii = 0, nn = cvPrcommVec.size(); ii < nn; ++ii) {
      unpack_size += cvPrcommVec[ii].unpack_size*data_size;
      pack_size   += cvPrcommVec[ii].packVec.size()*data_size;
    }

    assert( unpack_size == (ncv_g-ncv)*data_size);
    K* unpack_buf = new K[unpack_size];
    K* pack_buf   = new K[pack_size];
    assert(mrs->unpack_buf_v == NULL); mrs->unpack_buf_v = (void*) unpack_buf;
    assert(mrs->pack_buf_v   == NULL); mrs->pack_buf_v   = (void*) pack_buf;

    MPI_Datatype MPI_TYPE = getMpiDatatype<K>();

    unpack_size = 0;
    pack_size   = 0;
    for (int ii = 0, nn = cvPrcommVec.size(); ii < nn; ++ii) {

      MPI_Irecv(unpack_buf+unpack_size,cvPrcommVec[ii].unpack_size*data_size,
                MPI_TYPE,cvPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
      unpack_size += cvPrcommVec[ii].unpack_size*data_size;

      // pack and send...
      /*
        for (int i =0, n = cvPrcommVec[ii].packVec.size(); i < n; ++i) {
        const int icv = cvPrcommVec[ii].packVec[i];
        assert((icv >= 0)&&(icv < ncv));
        pack_class(pack_buf+pack_size+i*data_size,s,icv,0); // no bits right now..
        }
      */

      switch (action) {
      case REPLACE_ROTATE_DATA: // default transformation
        for (int jj = 0, limit = cvPrcommVec[ii].transformVec.size(); jj < limit; ++jj) {
          if (cvPrcommVec[ii].transformVec[jj].has_R) {
            for (int i = cvPrcommVec[ii].transformVec[jj].start; i != cvPrcommVec[ii].transformVec[jj].end; ++i) {
              const int icv = cvPrcommVec[ii].packVec[i];
              pack_class(pack_buf+pack_size+i*data_size,s,icv,cvPrcommVec[ii].transformVec[jj].R);
            }

          } else {
            for (int i = cvPrcommVec[ii].transformVec[jj].start; i != cvPrcommVec[ii].transformVec[jj].end; ++i) {
              const int icv = cvPrcommVec[ii].packVec[i];
              pack_class(pack_buf+pack_size+i*data_size,s,icv);
            }
          }
        }
        break;

      case REPLACE_DATA:

        // no rotation even if the range exists
        for (int i =0, n = cvPrcommVec[ii].packVec.size(); i < n; ++i) {
          const int icv = cvPrcommVec[ii].packVec[i];
          assert((icv >= 0)&&(icv < ncv));
          pack_class(pack_buf+pack_size+i*data_size,s,icv); // no bits right now..
        }

        break;

      default:
        assert(0);
      }

      MPI_Issend(pack_buf+pack_size,cvPrcommVec[ii].packVec.size()*data_size,
                 MPI_TYPE,cvPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
      pack_size += cvPrcommVec[ii].packVec.size()*data_size;
    }
  }

  template<class T, typename K>
  void updateCvClassFinish(T * s) {

    map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
    if (iter == mpiRequestMap.end())
      CERR("updateCvDataFinish(T*): s is not mapped.");

    MpiRequestStuff * mrs = iter->second;
    mpiRequestMap.erase(iter);

    // now we wait for all messages to be received..
    MPI_Waitall(cvPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);

    // require an unpack for the class strcture..
    K* unpack_buf       = static_cast<K*>(mrs->unpack_buf_v);
    const int data_size = T::data_size();
    for (int icv = ncv; icv < ncv_g; ++icv) {
      //s[icv].unpack(unpack_buf+(icv-ncv)*data_size);
      unpack_class(unpack_buf+(icv-ncv)*data_size,s,icv);
    }

    delete[] unpack_buf;
    mrs->unpack_buf_v = NULL;
    delete[] mrs->recvRequestArray;
    mrs->recvRequestArray = NULL;

    // and finally, wait for the sends..
    MPI_Waitall(cvPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);

    K* pack_buf        = static_cast<K*>(mrs->pack_buf_v);
    delete[] pack_buf;
    mrs->pack_buf_v = NULL;

    delete[] mrs->sendRequestArray;
    mrs->sendRequestArray = NULL;

    delete mrs;
  }

  template<class T, typename K>
  void updateCvClass(T* s, const int action = REPLACE_ROTATE_DATA) {
    updateCvClassStart<T,K>(s, action);
    updateCvClassFinish<T,K>(s);
  }

  //==================================================
  // cv class updates with level-1 and level-2 ghosts..
  //==================================================

  template<class T, typename K>
  void updateCv2ClassStart(T* s, const int action = REPLACE_ROTATE_DATA ) {

    map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
    if (iter != mpiRequestMap.end())
      CERR("updateCv2DataStart(T*): s is already mapped. updateCv2DataFinish required on s.");

    MpiRequestStuff * mrs   = new MpiRequestStuff; // mrs == MpiRequestStuff
    mpiRequestMap[(void*)s] = mrs;

    // count requests and buf sizes. This is a way to
    // step forward simultaneously through the cvPrcommVec and
    // cv2PrcommVec, which should be monotonic in rank by
    // construction...
    const int UPDATE_TAG = 432234;
    const int data_size  = T::data_size();
    int rank_check       = -1; // monotonicity check...
    int unpack_size      = 0;
    int pack_size        = 0;
    int ii               = 0;
    int ii2              = 0;
    int rank_count       = 0;
    CvPrcomm * cvPrcomm;
    CvPrcomm * cv2Prcomm;
    int rank;

    while ((ii != int(cvPrcommVec.size()))||(ii2 != int(cv2PrcommVec.size()))) {

      cvPrcomm  = NULL;
      cv2Prcomm = NULL;
      rank = -1;

      if (ii == int(cvPrcommVec.size())) {
        assert(ii2 != int(cv2PrcommVec.size()));
        cv2Prcomm = &cv2PrcommVec[ii2++];
        rank = cv2Prcomm->getRank();
      }
      else if (ii2 == int(cv2PrcommVec.size())) {
        assert(ii != int(cvPrcommVec.size()));
        cvPrcomm = &cvPrcommVec[ii++];
        rank = cvPrcomm->getRank();
      }
      else if (cvPrcommVec[ii].getRank() < cv2PrcommVec[ii2].getRank()) {
        cvPrcomm = &cvPrcommVec[ii++];
        rank = cvPrcomm->getRank();
      }
      else if (cvPrcommVec[ii].getRank() > cv2PrcommVec[ii2].getRank()) {
        cv2Prcomm = &cv2PrcommVec[ii2++];
        rank = cv2Prcomm->getRank();
      }
      else {
        cvPrcomm = &cvPrcommVec[ii++];
        rank = cvPrcomm->getRank();
        cv2Prcomm = &cv2PrcommVec[ii2++];
        assert(rank == cv2Prcomm->getRank());
      }
      assert(rank != -1);
      assert(rank > rank_check);
      rank_check = rank;
      assert(cvPrcomm || cv2Prcomm); // could be one or both
      assert(cv2Prcomm); // I think we insist this on build -- could simplify the above logic
      if (cvPrcomm) {
        unpack_size += cvPrcomm->unpack_size*data_size;
        pack_size += cvPrcomm->packVec.size()*data_size;
      }
      if (cv2Prcomm) {
        unpack_size += cv2Prcomm->unpack_size*data_size;
        pack_size += cv2Prcomm->packVec.size()*data_size;
      }
      ++rank_count;
    }

    assert(rank_count == cv2PrcommVec.size()); // cv2PrcommVec has at least the ranks from cvPrcommVec
    assert(unpack_size == (ncv_g2-ncv)*data_size);

    // allocate request arrays...
    mrs->sendRequestArray = new MPI_Request[rank_count];
    mrs->recvRequestArray = new MPI_Request[rank_count];

    // and pack and send...
    assert( mrs->unpack_buf_v == NULL);
    K* unpack_buf     = new K[unpack_size];
    mrs->unpack_buf_v = static_cast<void*>(unpack_buf);

    assert( mrs->pack_buf_v   == NULL);
    K* pack_buf       = new K[pack_size];
    mrs->pack_buf_v   = static_cast<void*>(pack_buf);

    MPI_Datatype MPI_TYPE = getMpiDatatype<K>();

    unpack_size = 0;
    pack_size   = 0;
    ii          = 0;
    ii2         = 0;
    rank_count  = 0;
    while ((ii != int(cvPrcommVec.size()))||(ii2 != int(cv2PrcommVec.size()))) {

      cvPrcomm  = NULL;
      cv2Prcomm = NULL;
      rank      = -1;

      if (ii == int(cvPrcommVec.size())) {
        assert(ii2 != int(cv2PrcommVec.size()));
        cv2Prcomm = &cv2PrcommVec[ii2++];
        rank = cv2Prcomm->getRank();
      }
      else if (ii2 == int(cv2PrcommVec.size())) {
        assert(ii != int(cvPrcommVec.size()));
        cvPrcomm = &cvPrcommVec[ii++];
        rank = cvPrcomm->getRank();
      }
      else if (cvPrcommVec[ii].getRank() < cv2PrcommVec[ii2].getRank()) {
        cvPrcomm = &cvPrcommVec[ii++];
        rank = cvPrcomm->getRank();
      }
      else if (cvPrcommVec[ii].getRank() > cv2PrcommVec[ii2].getRank()) {
        cv2Prcomm = &cv2PrcommVec[ii2++];
        rank = cv2Prcomm->getRank();
      }
      else {
        cvPrcomm = &cvPrcommVec[ii++];
        rank = cvPrcomm->getRank();
        cv2Prcomm = &cv2PrcommVec[ii2++];
        assert(rank == cv2Prcomm->getRank());
      }
      assert(rank != -1);
      assert(cvPrcomm || cv2Prcomm); // could be one or both
      assert(cv2Prcomm); // see above

      // post irecv...
      const int unpack_size0 = unpack_size;
      const int pack_size0   = pack_size;

      // pack...
      if (cvPrcomm) {

        /*
          for (int i = 0, n = cvPrcomm->packVec.size(); i < n; ++i) {
          const int icv = cvPrcomm->packVec[i];
          assert((icv >= 0)&&(icv < ncv));
          pack_class(pack_buf+pack_size+i*data_size,s,icv,0); // no bits right now..
          }
        */

        switch (action) {

        case REPLACE_ROTATE_DATA: // default transformation
          for (int jj = 0, limit = cvPrcomm->transformVec.size(); jj < limit; ++jj) {
            if (cvPrcomm->transformVec[jj].has_R) {
              for (int i = cvPrcomm->transformVec[jj].start; i != cvPrcomm->transformVec[jj].end; ++i) {
                const int icv = cvPrcomm->packVec[i];
                pack_class(pack_buf+pack_size+i*data_size,s,icv,cvPrcomm->transformVec[jj].R);
              }

            } else {
              for (int i = cvPrcomm->transformVec[jj].start; i != cvPrcomm->transformVec[jj].end; ++i) {
                const int icv = cvPrcomm->packVec[i];
                pack_class(pack_buf+pack_size+i*data_size,s,icv);
              }
            }
          }
          break;

        case REPLACE_DATA:

          // no rotation even if the range exists
          for (int i =0, n = cvPrcomm->packVec.size(); i < n; ++i) {
            const int icv = cvPrcomm->packVec[i];
            assert((icv >= 0)&&(icv < ncv));
            pack_class(pack_buf+pack_size+i*data_size,s,icv); // no bits right now..
          }

          break;

        default:
          assert(0);
        }

        unpack_size += cvPrcomm->unpack_size*data_size;
        pack_size += cvPrcomm->packVec.size()*data_size;
      }

      if (cv2Prcomm) {

        /*
          for (int i = 0, n = cv2Prcomm->packVec.size(); i < n; ++i) {
          const int icv = cv2Prcomm->packVec[i];
          assert((icv >= 0)&&(icv < ncv));
          pack_class(pack_buf+pack_size+i*data_size,s,icv,0); // no bits right now..
          }
        */

        switch (action) {

        case REPLACE_ROTATE_DATA: // default transformation
          for (int jj = 0, limit = cv2Prcomm->transformVec.size(); jj < limit; ++jj) {
            if (cv2Prcomm->transformVec[jj].has_R) {
              for (int i = cv2Prcomm->transformVec[jj].start; i != cv2Prcomm->transformVec[jj].end; ++i) {
                const int icv = cv2Prcomm->packVec[i];
                pack_class(pack_buf+pack_size+i*data_size,s,icv,cv2Prcomm->transformVec[jj].R);
              }

            } else {
              for (int i = cv2Prcomm->transformVec[jj].start; i != cv2Prcomm->transformVec[jj].end; ++i) {
                const int icv = cv2Prcomm->packVec[i];
                pack_class(pack_buf+pack_size+i*data_size,s,icv);
              }
            }
          }
          break;

        case REPLACE_DATA:

          // no rotation even if the range exists

          for (int i =0, n = cv2Prcomm->packVec.size(); i < n; ++i) {
            const int icv = cv2Prcomm->packVec[i];
            assert((icv >= 0)&&(icv < ncv));
            pack_class(pack_buf+pack_size+i*data_size,s,icv); // no bits right now..
          }

          break;

        default:
          assert(0);
        }

        unpack_size += cv2Prcomm->unpack_size*data_size;
        pack_size += cv2Prcomm->packVec.size()*data_size;
      }

      // send/recv...
      MPI_Irecv(unpack_buf+unpack_size0,unpack_size-unpack_size0,
                MPI_TYPE,rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[rank_count]));

      MPI_Issend(pack_buf+pack_size0,pack_size-pack_size0,
                 MPI_TYPE,rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[rank_count]));

      ++rank_count;
    }

    assert(rank_count == cv2PrcommVec.size());
  }

  template<class T, typename K>
  void updateCv2ClassFinish(T * s) {

    map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
    if (iter == mpiRequestMap.end())
      CERR("updateCv2DataFinish(T*): s is not mapped.");

    MpiRequestStuff * mrs = iter->second;
    mpiRequestMap.erase(iter);

    // wait for the recv's...
    MPI_Waitall(cv2PrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);

    K* unpack_buf = static_cast<K*>(mrs->unpack_buf_v);

    // once all recv's available, we need to unpack...
    const int data_size = T::data_size();
    int unpack_size     = 0;
    int icv_g           = ncv;
    int icv_g2          = ncv_g;
    int ii              = 0;
    int ii2             = 0;
    CvPrcomm * cvPrcomm;
    CvPrcomm * cv2Prcomm;
    int rank;

    while ((ii != int(cvPrcommVec.size()))||(ii2 != int(cv2PrcommVec.size()))) {

      cvPrcomm = NULL;
      cv2Prcomm = NULL;
      rank = -1;

      if (ii == int(cvPrcommVec.size())) {
        assert(ii2 != int(cv2PrcommVec.size()));
        cv2Prcomm = &cv2PrcommVec[ii2++];
        rank = cv2Prcomm->getRank();
      }
      else if (ii2 == int(cv2PrcommVec.size())) {
        assert(ii != int(cvPrcommVec.size()));
        cvPrcomm = &cvPrcommVec[ii++];
        rank = cvPrcomm->getRank();
      }
      else if (cvPrcommVec[ii].getRank() < cv2PrcommVec[ii2].getRank()) {
        cvPrcomm = &cvPrcommVec[ii++];
        rank = cvPrcomm->getRank();
      }
      else if (cvPrcommVec[ii].getRank() > cv2PrcommVec[ii2].getRank()) {
        cv2Prcomm = &cv2PrcommVec[ii2++];
        rank = cv2Prcomm->getRank();
      }
      else {
        cvPrcomm = &cvPrcommVec[ii++];
        rank = cvPrcomm->getRank();
        cv2Prcomm = &cv2PrcommVec[ii2++];
        assert(rank == cv2Prcomm->getRank());
      }
      assert(rank != -1);
      assert(cvPrcomm || cv2Prcomm); // could be one or both

      // unpack...
      if (cvPrcomm) {
        for (int i = 0; i < cvPrcomm->unpack_size; ++i) {
          unpack_class(unpack_buf+unpack_size+i*data_size,s,icv_g);
          icv_g++;
        }
        unpack_size += cvPrcomm->unpack_size*data_size;
      }

      if (cv2Prcomm) {
        for (int i = 0; i < cv2Prcomm->unpack_size; ++i) {
          unpack_class(unpack_buf+unpack_size+i*data_size,s,icv_g2);
          ++icv_g2;
        }
        unpack_size += cv2Prcomm->unpack_size*data_size;
      }
    }

    // make sure we ended up where we expected...
    assert(icv_g == ncv_g);
    assert(icv_g2 == ncv_g2);

    // cleanup...
    delete[] unpack_buf;
    mrs->unpack_buf_v = NULL;
    delete[] mrs->recvRequestArray;
    mrs->recvRequestArray = NULL;

    // now we wait for all messages to be sent and received...
    MPI_Waitall(cv2PrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);

    // as soon as we are all sent, we can clear the buffer...
    K* pack_buf = static_cast<K*>(mrs->pack_buf_v);
    delete[] pack_buf;
    mrs->pack_buf_v = NULL; // nullify for desructor check
    delete[] mrs->sendRequestArray;
    mrs->sendRequestArray = NULL;

    delete mrs;
  }

  // a blocking version of the complete update...
  template<class T, typename K>
  void updateCv2Class(T * s, const int action = REPLACE_ROTATE_DATA) {
    updateCv2ClassStart<T,K>(s,action);
    updateCv2ClassFinish<T,K>(s);
  }

  // =================================================
  // these routines update level-1 ghosts only...
  // =================================================

  // updateCvData(int * s...

#define T int
#define PACK_BUF pack_buf_int
#define MPI_T MPI_INT
#define UPDATE_TAG 12121
#include "updateCvDN.hpp"
#undef T
#undef PACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateCvData(int8 * s...

#define T int8
#define PACK_BUF pack_buf_int8
#define MPI_T MPI_INT8
#define UPDATE_TAG 12126
#include "updateCvDN.hpp"
#undef T
#undef PACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateCvData(double * s and double (*s)[3]...

#define T double
#define PACK_BUF pack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 12122
#include "updateCvDN.hpp"
#include "updateCvDN3.hpp"
#include "updateCvDN33.hpp"
#undef T
#undef PACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateCvDataSeparateGhosts(int * s,int *sg and int (*s)[3],int (*sg)[3]...

#define T int
#define PACK_BUF pack_buf_int
#define MPI_T MPI_INT
#define UPDATE_TAG 12123
#include "updateCvDataSeparateGhosts.hpp"
#undef T
#undef PACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateCvDataSeparateGhosts(double * s,double *sg and double (*s)[3],double (*sg)[3]...

#define T double
#define PACK_BUF pack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 12124
#include "updateCvDataSeparateGhosts.hpp"
#undef T
#undef PACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateCvDataSeparateGhosts(int8 * s,int8 *sg and int8 (*s)[3],int8 (*sg)[3]...

#define T int8
#define PACK_BUF pack_buf_int8
#define MPI_T MPI_INT8
#define UPDATE_TAG 12125
#include "updateCvDataSeparateGhosts.hpp"
#undef T
#undef PACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateCvdData(int * s,int * sg...

#define T int
#define PACK_BUF pack_buf_int
#define MPI_T MPI_INT
#define UPDATE_TAG 12126
#include "updateCvdDN.hpp"
#undef T
#undef PACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateCvdData(int8 * s,int8 * sg...

#define T int8
#define PACK_BUF pack_buf_int8
#define MPI_T MPI_INT8
#define UPDATE_TAG 12127
#include "updateCvdDN.hpp"
#undef T
#undef PACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateCvdData(double * s,double * sg and double (*s)[3], double (*sg)[3]...

#define T double
#define PACK_BUF pack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 12128
#include "updateCvdDN.hpp"
#include "updateCvdDN3.hpp"
#undef T
#undef PACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // =================================================
  // these routines update level-1 AND level-2 ghosts...
  // =================================================

  // updateCv2Data(int * s...

#define T int
#define PACK_BUF pack_buf_int
#define UNPACK_BUF unpack_buf_int
#define MPI_T MPI_INT
#define UPDATE_TAG 21121
#include "updateCv2DN.hpp"
#undef T
#undef PACK_BUF
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateCv2Data(double * s and double (*s)[3]...

#define T double
#define PACK_BUF pack_buf_double
#define UNPACK_BUF unpack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 21122
#include "updateCv2DN.hpp"
#include "updateCv2DN3.hpp"
#include "updateCv2DN33.hpp"
#undef T
#undef PACK_BUF
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateCv2DataReverse(double * s ...

#define T double
#define PACK_BUF pack_buf_double
#define UNPACK_BUF unpack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 21123
#include "updateCv2DNReverse.hpp"
#undef T
#undef PACK_BUF
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // =================================================
  // these routines update nodes with addition...
  // =================================================

  // updateNoData(int * s and int (*s)[3]...

#define T int
#define PACK_BUF pack_buf_int
#define UNPACK_BUF unpack_buf_int
#define MPI_T MPI_INT
#define UPDATE_TAG 32121
#include "updateNoData.hpp"
#undef T
#undef PACK_BUF
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateNoData(double * s and double (*s)[3]...

#define T double
#define PACK_BUF pack_buf_double
#define UNPACK_BUF unpack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 32122
#include "updateNoData.hpp"
#undef T
#undef PACK_BUF
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateNoData(int8 * s and int8 (*s)[3]...

#define T int8
#define PACK_BUF pack_buf_int8
#define UNPACK_BUF unpack_buf_int8
#define MPI_T MPI_INT8
#define UPDATE_TAG 32123
#include "updateNoData.hpp"
#undef T
#undef PACK_BUF
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // =================================================
  // these routines update faces: action must be specified
  // as one of REPLACE_DATA, ADD_DATA, SUBTRACT_DATA
  // =================================================

  // updateFaData(int * s...

#define T int
#define PACK_BUF pack_buf_int
#define UNPACK_BUF unpack_buf_int
#define MPI_T MPI_INT
#define UPDATE_TAG 42121
#include "updateFaDN.hpp"
#undef T
#undef PACK_BUF
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateFaData(double * s and double (*s)[3]...

#define T double
#define PACK_BUF pack_buf_double
#define UNPACK_BUF unpack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 42122
#include "updateFaDN.hpp"
#include "updateFaDN3.hpp"
#undef T
#undef PACK_BUF
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateFaData(int8 * s...

#define T int8
#define PACK_BUF pack_buf_int8
#define UNPACK_BUF unpack_buf_int8
#define MPI_T MPI_INT8
#define UPDATE_TAG 42123
#include "updateFaDN.hpp"
#undef T
#undef PACK_BUF
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

protected:

  void registerStats(Param * param,const bool b_init);

  void interpFromData();
  int prepareInterp(DataContainer& dc,bool first);

public:

  void solvePotentialFlow(double * phi) {

    // -----------------------------------------
    // "PF" == potential flow
    // specify 1 or more PF_INLETS as follows:
    //
    // PF_INLET <zone1> MDOT=1.23
    // PF_INLET <zone2> MDOT=2.1
    //
    // and 1 or more outlets. No MDOT here: it is calculated.
    //
    // PF_OUTLET <zone3>
    //
    // If you have more than one outlet, the solver will distribute the "flow"
    // evenly according to area. If you want another distribution, set
    // one of them to a PF_INLET and give a negative value.
    // -----------------------------------------

    COUT1("solvePotentialFlow()");

    assert(phi);
    for (int icv = 0; icv < ncv_g; ++icv) phi[icv] = 0.0;

    double * rhoun_bf = new double[nbf];
    FOR_IBF rhoun_bf[ibf] = 0.0;

    // ----------------------------------
    // cycle through PF_INLET's first...
    // ----------------------------------
    double inlet_mdot_sum = 0.0;
    int ierr = 0;

    FOR_PARAM_MATCHING("PF_INLET") {
      int this_ierr = 0;
      int iarg = 0;
      const string name = param->getString(iarg++);
      BfZone * bf_zone = getBfZone(name);
      const string token = param->getString(iarg++);
      double mdot;
      bool b_mdot = false;
      if (token == "MDOT") {
        mdot = param->getDouble(iarg++);
        b_mdot = true;
      }
      else {
        if (mpi_rank == 0) cout << " > Unrecognized PF_INLET token: " << token << ". Skipping." << endl;
      }
      // error checking...
      if (bf_zone == NULL) {
        this_ierr = 1;
        if (mpi_rank == 0) cout << " > PF_INLET cannot find zone for name \"" << name << "\"" << endl;
      }
      if (!b_mdot) {
        this_ierr = 1;
        if (mpi_rank == 0) cout << " > PF_INLET missing MDOT." << endl;
      }
      if (this_ierr) {
        ierr = 1;
        continue;
      }

      // looks good...
      inlet_mdot_sum += mdot;
      const double rhoun = mdot/bf_zone->area_global;

      if (mpi_rank == 0)
        cout << " > PF_INLET \"" << name << "\": MDOT=" << mdot << " (rhoun=" << rhoun << ")" << endl;

      for (int ibf = bf_zone->ibf_f; ibf <= bf_zone->ibf_l; ++ibf) {
        rhoun_bf[ibf] = -rhoun*area_bf[ibf]; // inlet negative by convention
      }
    }

    // ----------------------------------
    // then PF_OUTLET's...
    // ----------------------------------
    double outlet_rhoun = 0.0;

    FOR_PARAM_MATCHING("PF_OUTLET") {
      int this_ierr = 0;
      int iarg = 0;
      const string name = param->getString(iarg++);
      BfZone * bf_zone = getBfZone(name);
      // error checking...
      if (bf_zone == NULL) {
        this_ierr = 1;
        if (mpi_rank == 0) cout << " > PF_OUTLET cannot find zone for name \"" << name << "\"" << endl;
      }
      if (this_ierr) {
        ierr = 1;
        continue;
      }
      // looks good...
      // accumulate areas in outlet_rhoun...
      outlet_rhoun += bf_zone->area_global;
    }

    // inlets and outlets must balance in potential flow...
    if (inlet_mdot_sum != 0.0) {
      if (outlet_rhoun == 0.0) { // recall this is holding the area right now -- should be positive
        if (mpi_rank == 0) cout << " > inlets do not balance. You must specify atleast one PF_OUTLET <zone>" << endl;
        ierr = 1;
      }
      else {
        outlet_rhoun = inlet_mdot_sum/outlet_rhoun; // recall outlet_rhoun contained the outlet area
      }
    }

    // return if ierr...
    if (ierr != 0) {
      delete[] rhoun_bf;
      return;
    }

    FOR_PARAM_MATCHING("PF_OUTLET") {
      int iarg = 0;
      const string name = param->getString(iarg++);
      BfZone * bf_zone = getBfZone(name);
      assert(bf_zone != NULL);
      if (mpi_rank == 0)
        cout << " > PF_OUTLET \"" << name << "\": computed MDOT=" << outlet_rhoun*bf_zone->area_global << " (rhoun=" << outlet_rhoun << ")" << endl;
      for (int ibf = bf_zone->ibf_f; ibf <= bf_zone->ibf_l; ++ibf) {
        rhoun_bf[ibf] = outlet_rhoun*area_bf[ibf]; // outlet positive by convention
      }
    }

    // ----------------------------------
    // solver...
    // ----------------------------------

    double * rhs = new double[ncv];
    for (int icv = 0; icv < ncv; ++icv)
      rhs[icv] = 0.0;

    for (int ibf = 0; ibf < nbf; ++ibf) {
      const int icv = cvobf[ibf];
      rhs[icv] += rhoun_bf[ibf];
    }

    // check rhs...

    double my_sum = 0.0;
    for (int icv = 0; icv < ncv; ++icv)
      my_sum += rhs[icv];
    double sum;
    MPI_Reduce(&my_sum,&sum,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    if (mpi_rank == 0)
      cout << " > sum(rhs) (should be zero): " << sum << endl;

    const double pf_zero = getDoubleParam("PF_ZERO",1.0E-6);
    const int pf_maxiter = getIntParam("PF_MAXITER",2000);

    double * A = new double[cvocv_i[ncv]];
    buildCvLaplacian(A);
    solveCvCg(phi,A,rhs,pf_zero,pf_maxiter,true); // verbose

    delete[] A;
    delete[] rhs;

    MiscUtils::dumpRange(phi,ncv,"phi");

    delete[] rhoun_bf;

  }

  void buildCvLaplacian(double * A);
  void buildCv2Laplacian(double * A);

  void calcCvResidual(double *res,const double* phi,const double *A,const double *rhs);
  int solveCvCg(double * phi,const double * const A,const double * const rhs,const double zero,const int maxiter,const bool verbose);
  int solveCvCg(double (*u)[3],const double * const A,const double (*const rhs)[3],const double zero,const int maxiter,const bool verbose);
  int solveCvJacobi(double * phi,const double *A,const double *rhs, const double zero,const double relax,const int maxiter,const bool verbose);
  int solveCvPatr(double *phi,const double * A,const double * At, const double *rhs,
                  const double zero,const double relax1,const double relax2,const int maxiter,const bool verbose);

  void calcCvResidual(double (*res)[3],const double (*u)[3],const double *A,const double (*rhs)[3]);
  int solveCvJacobi(double (*__restrict__ u)[3],const double *__restrict__ A,const double (*__restrict__ rhs)[3],
                    const double zero,const double relax,const int maxiter,const bool verbose);
  int solveCvPatr(double (*u)[3],const double * A,const double * At, const double (*rhs)[3],
                  const double zero,const double relax1,const double relax2,const int maxiter,const bool verbose);
  int solveCvTim(double (*u)[3],const double *A,const double *As,const double (*rhs)[3],
                 const double zero,const double tau_bound,const double relax,const int maxiter,const bool verbose);

  void calcCvResidual(double (*res)[3],const double (*u)[3],const double *A,const double (*A_diag)[3],const double (*rhs)[3]);
  int solveCvJacobi(double (*__restrict__ u)[3],const double *__restrict__ A,const double (*__restrict__ A_diag)[3],const double (*__restrict__ rhs)[3],
                    const double zero,const double relax,const int maxiter,const bool verbose);
  int solveCvPatr(double (*u)[3],const double * A,const double * At, const double (*A_diag)[3],const double (*rhs)[3],
                  const double zero,const double relax1,const double relax2,const int maxiter,const bool verbose);



  void calcCvGradLeastSquares(double (*dphidx)[3],double * phi);
  void calcCvGradLeastSquares(double (*dphidx)[3],double * phi,double * minus_dphidn_dA_bf);

  inline int8 getNcvGlobal() const { return ncv_global; }
  inline int8 getNbfGlobal() const { return nbf_global; }
  inline int8 getNfaGlobal() const { return nfa_global; }
  inline int8 getNnoGlobal() const { return nno_global; }
  inline int8 getNefGlobal() const { return nef_global; }

  inline int8 getIcvGlobal(const int icv) const {
    assert((icv >= 0)&&(icv < ncv));
    assert(icv_global);
    return icv_global[icv];
  }

  inline int8 getIfaGlobal(const int ifa) const {
    assert((ifa >= 0)&&(ifa < nfa));
    assert(ifa_global);
    return ifa_global[ifa];
  }

  inline int8 getAbsIfaGlobal(const int ifa) const {
    assert((ifa >= 0)&&(ifa < nfa));
    assert(ifa_global);
    return max(ifa_global[ifa],-ifa_global[ifa]-1);
  }

  virtual bool getSpecifierByIndex(string& name, const int index) {
    if ((index >= 0) && (index < int(bfZoneVec.size()) )) {
      name = bfZoneVec[index].getName();
      return true;
    }
    return false;  // index out of bounds
  }

  // lagrangian particles...

  class LpHelper {
  public:
    string name;
    int index;
    int8 lb_cost; // load balance cost -- set this using loadBalanceHook()
    int8 nlp_global;
    DistributedDataExchanger * dde_cv_striped; // pull and push from the cv striped dist (the one for reading/writing)
    int8 * lpora_cv_striped;
    int* cvolp;
    int* size_ptr;
    int size;

    LpHelper() {
      name = "";
      index = -1;
      lb_cost = 0; // default is free
      nlp_global = -1;
      dde_cv_striped = NULL;
      lpora_cv_striped = NULL;
      cvolp = NULL;
      size_ptr = NULL;
      size = -1; // used in BasicPostpro
    }

    void clear() {
      if (dde_cv_striped) {
        delete dde_cv_striped;
        dde_cv_striped = NULL;
      }
      DELETE(lpora_cv_striped);
      DELETE(cvolp);
    }

  };
  vector<LpHelper> lpHelperVec;
  map<const string,int> lpHelperNameMap;

  template<class T>
  void registerLp(T * &lp,const string& name,int &np) {

    size_t colon_pos = name.find(":");
    assert(colon_pos == string::npos); // should be name of particle class
    map<const string,int>::iterator iter = lpHelperNameMap.find(name);
    assert(iter == lpHelperNameMap.end()); // no name repeats

    // populate helper vec

    assert(lpHelperVec.size() < 10); // see note on reserve in constructor
    lpHelperVec.push_back(LpHelper());
    const int lp_index = lpHelperVec.size()-1;
    lpHelperNameMap[name]       = lp_index;
    lpHelperVec.back().name     = name;
    lpHelperVec.back().index    = lp_index;
    lpHelperVec.back().size_ptr = &np;

  }

  template<class T>
  void registerLpData(T * &lp,int &data,const string& name,const uint rw_bits) {

    size_t colon_pos = name.find(":");
    assert(colon_pos != string::npos); // names must be prepended with particle class name
    const string lp_name  = name.substr(0,colon_pos);

    // get index into helper vec

    map<const string,int>::iterator iter = lpHelperNameMap.find(lp_name);
    assert(iter != lpHelperNameMap.end());
    const int lp_index = iter->second;
    assert((lp_index >= 0)&&(lp_index < lpHelperVec.size()));
    assert(lpHelperVec[lp_index].name  == lp_name);
    assert(lpHelperVec[lp_index].index == lp_index);

    // now register the data

    const int topo = ((lp_index<<INDEX_SHIFT)|LP_DATA);
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA))
      CtiRegister::_registerData(lp,data,name,topo,rw_bits,*lpHelperVec[lp_index].size_ptr,lpHelperVec[lp_index].lpora_cv_striped,lpHelperVec[lp_index].dde_cv_striped);
    else
      CtiRegister::_registerData(lp,data,name,topo,rw_bits,*lpHelperVec[lp_index].size_ptr);

  }

  template<class T>
  void registerLpData(T * &lp,double &data,const string& name,const uint rw_bits) {

    size_t colon_pos = name.find(":");
    assert(colon_pos != string::npos); // names must be prepended with particle class name
    const string lp_name  = name.substr(0,colon_pos);

    // get index into helper vec

    map<const string,int>::iterator iter = lpHelperNameMap.find(lp_name);
    assert(iter != lpHelperNameMap.end());
    const int lp_index = iter->second;
    assert((lp_index >= 0)&&(lp_index < lpHelperVec.size()));
    assert(lpHelperVec[lp_index].name  == lp_name);
    assert(lpHelperVec[lp_index].index == lp_index);

    // now register the data

    const int topo = ((lp_index<<INDEX_SHIFT)|LP_DATA);
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA))
      CtiRegister::_registerData(lp,data,name,topo,rw_bits,*lpHelperVec[lp_index].size_ptr,lpHelperVec[lp_index].lpora_cv_striped,lpHelperVec[lp_index].dde_cv_striped);
    else
      CtiRegister::_registerData(lp,data,name,topo,rw_bits,*lpHelperVec[lp_index].size_ptr);

  }
  template<class T>
  void registerLpData(T * &lp,double (&data)[3],const string& name,const uint rw_bits) {

    size_t colon_pos = name.find(":");
    assert(colon_pos != string::npos); // names must be prepended with particle class name
    const string lp_name  = name.substr(0,colon_pos);

    // get index into helper vec

    map<const string,int>::iterator iter = lpHelperNameMap.find(lp_name);
    assert(iter != lpHelperNameMap.end());
    const int lp_index = iter->second;
    assert((lp_index >= 0)&&(lp_index < lpHelperVec.size()));
    assert(lpHelperVec[lp_index].name  == lp_name);
    assert(lpHelperVec[lp_index].index == lp_index);

    // now register the data

    const int topo = ((lp_index<<INDEX_SHIFT)|LP_DATA);
    if ((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA))
      CtiRegister::_registerData(lp,data,name,topo,rw_bits,*lpHelperVec[lp_index].size_ptr,lpHelperVec[lp_index].lpora_cv_striped,lpHelperVec[lp_index].dde_cv_striped);
    else
      CtiRegister::_registerData(lp,data,name,topo,rw_bits,*lpHelperVec[lp_index].size_ptr);

  }

  void registerLpDN(const string& name) {


    size_t colon_pos = name.find(":");
    assert(colon_pos != string::npos); // names must be prepended with particle class name
    const string lp_name  = name.substr(0,colon_pos);

    // add/get lp helper

    map<const string,int>::iterator iter = lpHelperNameMap.find(lp_name);
    int lp_index = -1;
    if (iter == lpHelperNameMap.end()) {
      lpHelperVec.push_back(LpHelper());
      lp_index = lpHelperVec.size()-1;
      lpHelperNameMap[lp_name]       = lp_index;
      lpHelperVec.back().name     = lp_name;
      lpHelperVec.back().index    = lp_index;
      lpHelperVec.back().size_ptr = &lpHelperVec.back().size; // internally managed size
      const int topo = ((lp_index<<INDEX_SHIFT)|LP_DATA);
      CtiRegister::_registerData(lp_name+":icv",topo,IN_DATA,READWRITE_DATA,*lpHelperVec[lp_index].size_ptr,lpHelperVec[lp_index].lpora_cv_striped,lpHelperVec[lp_index].dde_cv_striped);
    }
    else {
      lp_index = iter->second;
    }

    // now register the data

    const int topo = ((lp_index<<INDEX_SHIFT)|LP_DATA);
    CtiRegister::_registerData(name,topo,DN_DATA,READWRITE_DATA,*lpHelperVec[lp_index].size_ptr,lpHelperVec[lp_index].lpora_cv_striped,lpHelperVec[lp_index].dde_cv_striped);

  }
  void registerLpDN3(const string& name) {

    size_t colon_pos = name.find(":");
    assert(colon_pos != string::npos); // names must be prepended with particle class name
    const string lp_name  = name.substr(0,colon_pos);

    // add/get lp helper

    map<const string,int>::iterator iter = lpHelperNameMap.find(lp_name);
    int lp_index = -1;
    if (iter == lpHelperNameMap.end()) {
      lpHelperVec.push_back(LpHelper());
      lp_index = lpHelperVec.size()-1;
      lpHelperNameMap[lp_name]       = lp_index;
      lpHelperVec.back().name     = lp_name;
      lpHelperVec.back().index    = lp_index;
      lpHelperVec.back().size_ptr = &lpHelperVec.back().size; // internally managed size
      const int topo = ((lp_index<<INDEX_SHIFT)|LP_DATA);
      CtiRegister::_registerData(lp_name+":icv",topo,IN_DATA,READWRITE_DATA,*lpHelperVec[lp_index].size_ptr,lpHelperVec[lp_index].lpora_cv_striped,lpHelperVec[lp_index].dde_cv_striped);
    }
    else {
      lp_index = iter->second;
    }

    // now register the data

    const int topo = ((lp_index<<INDEX_SHIFT)|LP_DATA);
    CtiRegister::_registerData(name,topo,DN3_DATA,READWRITE_DATA,*lpHelperVec[lp_index].size_ptr,lpHelperVec[lp_index].lpora_cv_striped,lpHelperVec[lp_index].dde_cv_striped);

  }
  void clearLpData() {
    for (map<const string,CtiRegister::CtiData>::iterator it = CtiRegister::registeredDataMap.begin();
         it != CtiRegister::registeredDataMap.end(); ++it) {
      // if registration owns the particle
      if ((it->second.getMemFlag())&&(it->second.getUnindexedTopology() == LP_DATA)&&(it->second.getDataPtr() != NULL)) {
        it->second.clear();

        // reset size_ptr
        string name = it->first;
        size_t colon_pos = name.find(":");
        assert(colon_pos != string::npos);
        const string lp_name  = name.substr(0,colon_pos);
        map<const string,int>::iterator iter = lpHelperNameMap.find(lp_name);
        assert(iter != lpHelperNameMap.end());
        const int lp_index = iter->second;
        it->second.setSizePtr(*lpHelperVec[lp_index].size_ptr);
      }
    }
  }

  virtual CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,const string& name,list<CtiRegister::CtiData>& args,
                                                    const bool b_eval_func) {

    //if (mpi_rank == 0)
    //  cout << "StaticSolver::funcEvalCtiData() being asked to evaluate: \"" << name << "\"" << endl;

    if (name == "diff_fa") {

      if (args.size() != 1)
        return(CtiRegister::CTI_DATA_ARG_COUNT);

      list<CtiRegister::CtiData>::iterator arg = args.begin();

      if ((arg->getTopology() == CV_DATA)&&(arg->getType() == DN_DATA)) {

        double * v_ptr = createSignedFaD1Data(v);

        if ( b_eval_func) {
          double * arg_ptr = arg->getDNptr();

          FOR_INTERPROC_IFA {
            const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
            v_ptr[ifa] = -arg_ptr[icv0];
          }
          updateFaDataStart(v_ptr,SUBTRACT_DATA);
          FOR_INTERNAL_IFA {
            const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
            const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
            v_ptr[ifa] = arg_ptr[icv1] - arg_ptr[icv0];
          }
          updateFaDataFinish(v_ptr,SUBTRACT_DATA);

        }

        return CtiRegister::CTI_DATA_OK;

      }
      else {

        return CtiRegister::CTI_DATA_NOT_VALID;

      }

    }
    else if (name == "avg_fa") {

      if (args.size() != 1)
        return(CtiRegister::CTI_DATA_ARG_COUNT);

      list<CtiRegister::CtiData>::iterator arg = args.begin();

      if (arg->getTopology() == CV_DATA) {
        if (arg->getType() == DN_DATA) {

          double * v_ptr = createFaD1Data(v);

          if ( b_eval_func) {

            double * arg_ptr = arg->getDNptr();
            FOR_INTERPROC_IFA {
              const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
              v_ptr[ifa] = 0.5* arg_ptr[icv0];
            }
            updateFaDataStart(v_ptr,ADD_DATA);
            FOR_INTERNAL_IFA {
              const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
              const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
              v_ptr[ifa] = 0.5* (arg_ptr[icv1] + arg_ptr[icv0]);
            }
            updateFaDataFinish(v_ptr,ADD_DATA);
          }

          return CtiRegister::CTI_DATA_OK;

        }
        else if (arg->getType() == DN3_DATA) {

          double (*v_ptr)[3] = createFaD2Data(v);

          if ( b_eval_func ) {

            double (*arg_ptr)[3] = arg->getDN3ptr();
            FOR_INTERPROC_IFA {
              const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
              FOR_I3 v_ptr[ifa][i] = 0.5*arg_ptr[icv0][i];
            }
            updateFaDataStart(v_ptr,ADD_ROTATE_DATA);
            FOR_INTERNAL_IFA {
              const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
              const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
              FOR_I3 v_ptr[ifa][i] = 0.5*(arg_ptr[icv1][i] + arg_ptr[icv0][i]);
            }
            updateFaDataFinish(v_ptr,ADD_ROTATE_DATA);
          }

          return CtiRegister::CTI_DATA_OK;

        }
        else {

          return CtiRegister::CTI_DATA_NOT_VALID;

        }
      }
      else {

        return CtiRegister::CTI_DATA_NOT_VALID;

      }

    }
    else if (name == "div") {

      if (args.size() != 1)
        return(CtiRegister::CTI_DATA_ARG_COUNT);

      list<CtiRegister::CtiData>::iterator arg = args.begin();

      if (arg->getTopology() == CV_DATA) {
        if (arg->getType() == DN3_DATA) {

          double * v_ptr = createCvD1Data(v);

          if ( b_eval_func) {

            // alloc some ghosts for communication

            double (*arg_ptr)[3] = arg->getDN3ptr();
            double (*arg_g)[3] = new double[ncv_g-ncv][3];
            updateCvDataSeparateGhosts(arg_ptr,arg_g);

            // use transpose of cvocv_grad_coeff to compute div at cv

            FOR_ICV {
              const int coc_f = cvocv_i[icv];
              v_ptr[icv] = 0.0;
              FOR_I3 v_ptr[icv] += cvocv_grad_coeff[coc_f][i] * arg_ptr[icv][i];

              for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
                const int icv_nbr = cvocv_v[coc];
                if (icv_nbr < ncv) {
                  FOR_I3 v_ptr[icv] += cvocv_grad_coeff[coc][i] * arg_ptr[icv_nbr][i];
                }
                else {
                  FOR_I3 v_ptr[icv] += cvocv_grad_coeff[coc][i] * arg_g[icv_nbr-ncv][i];
                }
              }
            }

            delete[] arg_g;
          }

          return CtiRegister::CTI_DATA_OK;

        }
        else if (arg->getType() == DN33_DATA) {

          double (*v_ptr)[3] = createCvD2Data(v);

          if ( b_eval_func) {

            // alloc some ghosts for communication

            double (*arg_ptr)[3][3] = arg->getDN33ptr();
            double (*arg_g)[3][3] = new double[ncv_g-ncv][3][3];
            updateCvDataSeparateGhosts(arg_ptr,arg_g);

            // use transpose of cvocv_grad_coeff to compute div at cv

            FOR_ICV {
              const int coc_f = cvocv_i[icv];
              FOR_I3 v_ptr[icv][i] = 0.0;
              FOR_I3 FOR_J3 v_ptr[icv][i] += cvocv_grad_coeff[coc_f][j] * arg_ptr[icv][j][i];

              for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
                const int icv_nbr = cvocv_v[coc];
                if (icv_nbr < ncv) {
                  FOR_I3 FOR_J3 v_ptr[icv][i] += cvocv_grad_coeff[coc][j] * arg_ptr[icv_nbr][j][i];
                }
                else {
                  FOR_I3 FOR_J3 v_ptr[icv][i] += cvocv_grad_coeff[coc][j] * arg_g[icv_nbr-ncv][j][i];
                }
              }
            }

            delete[] arg_g;
          }

          return CtiRegister::CTI_DATA_OK;

        }

      }

      return CtiRegister::CTI_DATA_NOT_VALID;

    }
    else if (name == "grad") {

      if (args.size() != 1)
        return(CtiRegister::CTI_DATA_ARG_COUNT);

      list<CtiRegister::CtiData>::iterator arg = args.begin();

      // scalars only for now (until we start supporting DN33 io)
      if (arg->getTopology() == CV_DATA) {
        if (arg->getType() == DN_DATA) {

          double (*v_ptr)[3] = createCvD2Data(v);

          if ( b_eval_func) {
            double *arg_ptr = arg->getDNptr();
            double *arg_g = new double[ncv_g-ncv];
            updateCvDataSeparateGhosts(arg_ptr,arg_g);

            // use transpose of cvocv_grad_coeff to compute div at cv

            FOR_ICV {
              const int coc_f = cvocv_i[icv];
              FOR_I3 v_ptr[icv][i] = cvocv_grad_coeff[coc_f][i] * arg_ptr[icv];

              for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
                const int icv_nbr = cvocv_v[coc];
                if (icv_nbr < ncv) {
                  FOR_I3 v_ptr[icv][i] += cvocv_grad_coeff[coc][i] * arg_ptr[icv_nbr];
                }
                else {
                  FOR_I3 v_ptr[icv][i] += cvocv_grad_coeff[coc][i] * arg_g[icv_nbr-ncv];
                }
              }
            }
            delete[] arg_g;
          }

          return CtiRegister::CTI_DATA_OK;
        }
        else if (arg->getType() == DN3_DATA) {

          double (*v_ptr)[3][3] = createCvD3Data(v);

          if ( b_eval_func) {
            double (*arg_ptr)[3] = arg->getDN3ptr();
            double (*arg_g)[3] = new double[ncv_g-ncv][3];
            updateCvDataSeparateGhosts(arg_ptr,arg_g);

            // use transpose of cvocv_grad_coeff to compute div at cv

            FOR_ICV {
              const int coc_f = cvocv_i[icv];
              FOR_I3 FOR_J3 v_ptr[icv][i][j] = cvocv_grad_coeff[coc_f][j] * arg_ptr[icv][i];

              for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
                const int icv_nbr = cvocv_v[coc];
                if (icv_nbr < ncv) {
                  FOR_I3 FOR_J3 v_ptr[icv][i][j] += cvocv_grad_coeff[coc][j] * arg_ptr[icv_nbr][i];
                }
                else {
                  FOR_I3 FOR_J3 v_ptr[icv][i][j] += cvocv_grad_coeff[coc][j] * arg_g[icv_nbr-ncv][i];
                }
              }
            }
            delete[] arg_g;
          }

          return CtiRegister::CTI_DATA_OK;
        }

      }

      return CtiRegister::CTI_DATA_NOT_VALID;

    }
    else if (name == "curl") {

      if (args.size() != 1)
        return(CtiRegister::CTI_DATA_ARG_COUNT);

      list<CtiRegister::CtiData>::iterator arg = args.begin();

      if ((arg->getTopology() == CV_DATA)&&(arg->getType() == DN3_DATA)) {

        double (*v_ptr)[3] = createCvD2Data(v);

        if ( b_eval_func ) {

          // alloc some ghosts for communication

          double (*arg_ptr)[3] = arg->getDN3ptr();
          double (*arg_g)[3] = new double[ncv_g-ncv][3];
          updateCvDataSeparateGhosts(arg_ptr,arg_g);

          // use transpose of cvocv_grad_coeff to compute div at cv

          FOR_ICV {
            const int coc_f = cvocv_i[icv];
            v_ptr[icv][0] = CROSS_PRODUCT_0(cvocv_grad_coeff[coc_f],(double*)arg_ptr+3*icv);
            v_ptr[icv][1] = CROSS_PRODUCT_1(cvocv_grad_coeff[coc_f],(double*)arg_ptr+3*icv);
            v_ptr[icv][2] = CROSS_PRODUCT_2(cvocv_grad_coeff[coc_f],(double*)arg_ptr+3*icv);

            for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
              const int icv_nbr = cvocv_v[coc];
              if (icv_nbr < ncv) {
                v_ptr[icv][0] += CROSS_PRODUCT_0(cvocv_grad_coeff[coc],(double*)arg_ptr+3*icv_nbr);
                v_ptr[icv][1] += CROSS_PRODUCT_1(cvocv_grad_coeff[coc],(double*)arg_ptr+3*icv_nbr);
                v_ptr[icv][2] += CROSS_PRODUCT_2(cvocv_grad_coeff[coc],(double*)arg_ptr+3*icv_nbr);
              }
              else {
                v_ptr[icv][0] += CROSS_PRODUCT_0(cvocv_grad_coeff[coc],arg_g[icv_nbr-ncv]);
                v_ptr[icv][1] += CROSS_PRODUCT_1(cvocv_grad_coeff[coc],arg_g[icv_nbr-ncv]);
                v_ptr[icv][2] += CROSS_PRODUCT_2(cvocv_grad_coeff[coc],arg_g[icv_nbr-ncv]);
              }
            }
          }
          delete[] arg_g;
        }

        return CtiRegister::CTI_DATA_OK;

      }
      else {

        return CtiRegister::CTI_DATA_NOT_VALID;

      }

    }
    else if (name == "lap") {

      if (args.size() != 1)
        return(CtiRegister::CTI_DATA_ARG_COUNT);

      list<CtiRegister::CtiData>::iterator arg = args.begin();

      if (arg->getTopology() == CV_DATA) {
        if (arg->getType() == DN_DATA) {

          double * v_ptr = createCvD1Data(v);

          if ( b_eval_func) {

            // diff data across faces

            double *arg_diff_fa = new double[nfa];
            FOR_INTERPROC_IFA {
              const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
              arg_diff_fa[ifa] = -arg->dn(icv0);
            }
            updateFaDataStart(arg_diff_fa,SUBTRACT_DATA);
            FOR_INTERNAL_IFA {
              const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
              const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
              arg_diff_fa[ifa] = arg->dn(icv1) - arg->dn(icv0);
            }
            updateFaDataFinish(arg_diff_fa,SUBTRACT_DATA);

            // compute lap using area_over_delta

            FOR_ICV {
              v_ptr[icv] = 0.0;
              for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
                const int ifa = faocv_v[foc];
                double fa_sign;
                if ( cvofa[ifa][0] == icv ) {
                  fa_sign = 1.0;
                }
                else {
                  assert( cvofa[ifa][1] == icv );
                  fa_sign = -1.0;
                }
                v_ptr[icv] += arg_diff_fa[ifa]*area_over_delta_fa[ifa]*fa_sign;
              }
              v_ptr[icv] *= inv_vol[icv];
            }
            delete[] arg_diff_fa;

          }

          return CtiRegister::CTI_DATA_OK;

        }
        else if (arg->getType() == DN3_DATA) {

          double (*v_ptr)[3] = createCvD2Data(v);

          if ( b_eval_func ) {

            // diff data across faces

            double (*arg_diff_fa)[3] = new double[nfa][3];
            FOR_INTERPROC_IFA {
              const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
              FOR_I3 arg_diff_fa[ifa][i] = -arg->dn3(icv0,i);
            }
            updateFaDataStart(arg_diff_fa,SUBTRACT_ROTATE_DATA);
            FOR_INTERNAL_IFA {
              const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
              const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
              FOR_I3 arg_diff_fa[ifa][i] = arg->dn3(icv1,i) - arg->dn3(icv0,i);
            }
            updateFaDataFinish(arg_diff_fa,SUBTRACT_ROTATE_DATA);

            // compute lap using area_over_delta

            FOR_ICV {
              FOR_I3 v_ptr[icv][i] = 0.0;
              for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
                const int ifa = faocv_v[foc];
                double fa_sign;
                if ( cvofa[ifa][0] == icv ) {
                  fa_sign = 1.0;
                }
                else {
                  assert( cvofa[ifa][1] == icv );
                  fa_sign = -1.0;
                }
                FOR_I3 v_ptr[icv][i] += arg_diff_fa[ifa][i]*area_over_delta_fa[ifa]*fa_sign;
              }
            }
            FOR_ICV FOR_I3 v_ptr[icv][i] /= vol_cv[icv];
            delete[] arg_diff_fa;

          }

          return CtiRegister::CTI_DATA_OK;

        }
        else {

          return CtiRegister::CTI_DATA_NOT_VALID;

        }
      }
      else {

        return CtiRegister::CTI_DATA_NOT_VALID;

      }

    }
    else if (name == "const_cv") {

      if (args.size() != 1)
        return(CtiRegister::CTI_DATA_ARG_COUNT);

      list<CtiRegister::CtiData>::iterator arg = args.begin();

      if (arg->getType() == D_DATA) {
        double * v_ptr = createCvD1Data(v);
        if ( b_eval_func) {
          const double val = arg->d();
          FOR_ICV v_ptr[icv] = val;
        }
        return CtiRegister::CTI_DATA_OK;
      }
      return CtiRegister::CTI_DATA_NOT_VALID;

    }

    return CtiRegister::CTI_DATA_NOT_FOUND;

  }

  template<typename T>
  void resizeCvDataWithGhosts(T* &arr)  {

    assert( arr != NULL);
    T *tmp = new T[ncv_g];
    for (int icv = 0; icv < ncv; ++icv)
      tmp[icv] = arr[icv];

    updateCvData(tmp);

    delete[] arr; arr = tmp;

  }

  // some prototype filtering routines -- to support dynamic models.
  // these are based on integration of the the pyramids that connect
  // two adjacent voronoi sites that share a face.

  void setupFilterMoments();
  void filterCvR1(double* phif, const double* phi);
  void filterCvR2(double (*phif)[3], const double (*phi)[3]);

  //-------------------------------------------------------------
  // Geometric multigrid stuff -- see StaticSolver_multigrid.cpp
  //-------------------------------------------------------------

  // coloring for agglomeration
  void splitOrphanedColors(int8 * color_cv,int8 &ncolor);

  // smoothers
  void smoothCvJacobi(double* phi,double *tmp,const double* inv_diag,const double *A,const double *rhs,const int nsmooth,const double relax);
  void smoothCvJacobi(double (*u)[3],double (*tmp)[3],const double* inv_diag,const double *A,const double (*rhs)[3],const int nsmooth,const double relax);
  void smoothCvJacobi(double (*u)[3],double (*tmp)[3],const double (*inv_diag)[3],const double *A,const double (*A_diag)[3],const double (*rhs)[3],const int nsmooth,const double relax);
  void smoothCvGs(double* phi,const double* inv_diag,const double *A,const double *rhs,const int nsmooth,const double relax);
  void smoothCvGs(double (*u)[3],const double* inv_diag,const double *A,const double (*rhs)[3],const int nsmooth,const double relax);
  void smoothCvSgs(double* phi,const double* inv_diag,const double *A,const double *rhs,const int nsmooth,const double relax);
  void smoothCvSgs(double (*u)[3],const double* inv_diag,const double *A,const double (*rhs)[3],const int nsmooth,const double relax);
  void smoothCvPatr(double* phi,
                    double* r,double* p,double* Ar,double* Ap, // work arrays
                    const double* inv_diag,const double* A,const double* At,const double* rhs,
                    const int nsmooth,const double relax1,const double relax2);
  void smoothCvPatr(double (*u)[3],
                    double (*r)[3],double (*p)[3],double (*Ar)[3],double (*Ap)[3], // work arrays
                    const double* inv_diag,const double* A,const double* At,const double (*rhs)[3],
                    const int nsmooth,const double relax1,const double relax2);
  void smoothCvPatr(double (*u)[3],
                    double (*r)[3],double (*p)[3],double (*Ar)[3],double (*Ap)[3], // work arrays
                    const double (*inv_diag)[3],const double* A,const double* At,const double (*A_diag)[3], const double (*rhs)[3],
                    const int nsmooth,const double relax1,const double relax2);
  void smoothCvTim(double (*&u)[3],
                   double (*&w)[3], // work array
                   const double *inv_diag,const double* A,const double* invDAs,const double (*rhs)[3],
                   const int nsmooth,const double tau_bound,const double relax);
  void smoothCvCg(double *phi,
                  double *res,double* v,double *p, // work arrays
                  const double *inv_diag,const double *A,const double *rhs,const int nsmooth);
  void smoothCvCg(double (*u)[3],
                  double (*res)[3],double (*v)[3],double (*p)[3], // work arrays
                  const double *inv_diag,const double *A,const double (*rhs)[3],const int nsmooth);


  class CoarseGrid {
  public:

    // note that cvs/fas can be split over processor boundaries, so we sometimes store local properties

    // need pointer to solver to access base grid public data

    StaticSolver* solver;

    // cc data, cc == coarse cell

    int8 ncc_global;
    int ncc_in; // internal
    int ncc_a; // +active cc's
    int ncc;   // +inactive cc's
    int ncc_g; // +ghosts cc's (note that inactive are treated as ghosts too)
    int8* icc_global;
    int8* rbi_i;  // rank-bits-index inactive cc's (for pack)
    int8* rbi_g; // rank-bits-index ghost cc's (for unpack)
    double (*x_cc)[3];
    double * vol_cc;
    double * inv_vol_cc;
    int *ccocc_i;
    int *ccocc_v;
    int *cvocc_i; // coarse to fine
    int *cvocc_v;
    int *ccocv; // fine to coarse

    // cf data, cf == coarse face

    int ncf_aa; // cf's b/w active cc's
    int ncf_ai; // +cf's b/w active/inactive+ghost cc's
    int ncf; // +cf's b/w inactive cc's (thrown after ncf)
    vector<int8>* rbiVec_cf_i; // only for cf's w/ inactive cc (for pack)
    int (*ccocf)[2];
    double (*x_cf)[3];
    double (*n_cf)[3]; // proj area mag
    double *delta_cf;
    double *area_cf;
    int *faocf_i; // coarse to fine
    int *faocf_v; // signed!

    // cb data, cb == coarse boundary face

    int ncb_a; // cb's touching active cc's
    int ncb;   // +cb's touching inactive cc's
    int8 *rbi_cb_i; // only for inactive (for pack)
    double* area_cb;
    double* delta_cb;
    double (*n_cb)[3]; // proj area mag
    double (*x_cb)[3];
    int* ccocb;
    int* zone_cb;
    int* cbozn_i; // only active!
    int* bfocb_i; // coarse to fine
    int* bfocb_v;

    // TODO ef data and level-2 ghosts

    vector<Prcomm> ccIPrcommVec; // reduce inactive cc's onto active cc's
    vector<Prcomm> ccPrcommVec; // replace inactive/ghost cc's with active cc's
    vector<Prcomm> cbIPrcommVec; // reduce inactive cb's onto active cb's
    vector<Prcomm> cfIPrcommVec; // exchange cf's w/ an inactive cc
    map<const void*,MpiRequestStuff*> mpiRequestMap;

    CoarseGrid(StaticSolver* _solver) {
      solver = _solver;
      ncc_global = 0;
      ncc_in = 0;
      ncc_a = 0;
      ncc = 0;
      ncc_g = 0;
      icc_global = NULL;
      rbi_i = NULL;
      rbi_g = NULL;
      x_cc = NULL;
      vol_cc = NULL;
      ccocc_i = NULL;
      ccocc_v = NULL;
      cvocc_i = NULL;
      cvocc_v = NULL;
      ccocv = NULL;

      ncf_aa = 0;
      ncf_ai = 0;
      ncf = 0;
      rbiVec_cf_i = NULL;
      ccocf = NULL;
      x_cf = NULL;
      n_cf = NULL;
      delta_cf = NULL;
      area_cf = NULL;
      faocf_i = NULL;
      faocf_v = NULL;

      ncb_a = 0;
      ncb = 0;
      rbi_cb_i = NULL;
      area_cb = NULL;
      delta_cb = NULL;
      n_cb = NULL;
      x_cb = NULL;
      ccocb = NULL;
      zone_cb = NULL;
      cbozn_i = NULL;
      bfocb_i = NULL;
      bfocb_v = NULL;
    }

    ~CoarseGrid() {
      clear();
    }

    void clear() {
      DELETE(icc_global);
      DELETE(rbi_i);
      DELETE(rbi_g);
      DELETE(x_cc);
      DELETE(vol_cc);
      DELETE(ccocc_i);
      DELETE(ccocc_v);
      DELETE(cvocc_i);
      DELETE(cvocc_v);
      DELETE(ccocv);

      DELETE(rbiVec_cf_i);
      DELETE(ccocf);
      DELETE(x_cf);
      DELETE(n_cf);
      DELETE(delta_cf);
      DELETE(area_cf);
      DELETE(faocf_i);
      DELETE(faocf_v);

      DELETE(rbi_cb_i);
      DELETE(area_cb);
      DELETE(delta_cb);
      DELETE(n_cb);
      DELETE(x_cb);
      DELETE(ccocb);
      DELETE(zone_cb);
      DELETE(cbozn_i);
      DELETE(bfocb_i);
      DELETE(bfocb_v);
    }

    void initCoarseGrid(const int level,const double agglomeration_factor = 8.0,const bool split_orphaned_colors = true);

    void buildCcIPrcomm();

    // updateCcIData(double * s and double (*s)[3]...

#define T double
#define UNPACK_BUF unpack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 52122
#include "updateCcIData.hpp"
#undef T
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

    // updateCcIDataReverse(double * s and double (*s)[3]...

#define T double
#define PACK_BUF pack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 52126
#include "updateCcIDataReverse.hpp"
#undef T
#undef PACK_BUF
#undef MPI_T
#undef UPDATE_TAG

    void buildCcPrcomm();

    // updateCcData(double * s and double (*s)[3]...

#define T double
#define PACK_BUF pack_buf_double
#define UNPACK_BUF unpack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 52123
#include "updateCcData.hpp"
#undef T
#undef PACK_BUF
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

    void buildCbIPrcomm();

    // updateCbIData(double * s and double (*s)[3]...

#define T double
#define UNPACK_BUF unpack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 52124
#include "updateCbIData.hpp"
#undef T
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

    void buildCfIPrcomm();

    // updateCfIData(double * s and double (*s)[3]...

#define T double
#define PACK_BUF pack_buf_double
#define UNPACK_BUF unpack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 52125
#include "updateCfIData.hpp"
#undef T
#undef PACK_BUF
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

    void prolongCcData(double* phi_cv,const double* phi_cc);
    void prolongCcDataAndUpdateGuess(double* phi_cv,const double* phi_cc);
    void prolongCcData(double (*u_cv)[3],const double (*u_cc)[3]);
    void prolongCcDataAndUpdateGuess(double (*u_cv)[3],const double (*u_cc)[3],const double relax = 1.0);


    void restrictCcData(double* phi_cc,const double* phi_cv);
    void restrictExtrinsicCcData(double* phiV_cc,const double* phiV_cv);
    void restrictCcData(double (*u_cc)[3],const double (*u_cv)[3]);
    void restrictExtrinsicCcData(double (*uV_cc)[3],const double (*uV_cv)[3]);

    void restrictCfData(double* phi_cf,const double* phi_fa);
    void restrictCfData(double (*u_cf)[3],const double (*u_fa)[3]);
    void restrictSignedCfData(double* phi_cf,const double* phi_fa);
    void restrictExtrinsicSignedCfData(double* phiA_cf,const double* phiA_fa);

    void restrictCbData(double* phi_cb,const double* phi_bf,const int izone = -1);
    void restrictExtrinsicCbData(double* phiA_cb,const double* phiA_bf,const int izone = -1);
    void restrictCbData(double (*u_cb)[3],const double (*u_bf)[3],const int izone = -1);

  };

  virtual void initFromCoarseGrid(CoarseGrid* cg);

  // =============================
  // LpTracer-related stuff...
  // =============================

  int pointIsInside(const double xp[3], int &icv_ret);
  int initPlaneInjector(double * &cdf_area_rank, double * &cdf_area_tri, vector<pair<SimpleTri,int> >& triVec, const double * const xp_plane, const double * const np_plane);
  int initZoneInjector(double * &cdf_area_tri, int (* &spost_zn)[3] , double (* &xp_zn)[3], int& nst_zn, int& nsp_zn, const string zone_name);
  int addStructuredParticlesInBox(vector<LpDataBase>& lpDataVec, const double xBox[6], const int npx, const int npy, const int npz, const double dp, const string injectorName);
  int addEachCellParticlesInBox(vector<LpDataBase>& lpDataVec, const double xBox[6], const int np_in_cell, const double dp, string injectorName);

};

#endif
