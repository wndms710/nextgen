#ifndef _STRIPED_MESH_HPP_
#define _STRIPED_MESH_HPP_

#include "ByteSwap.hpp"
#include "MiscUtils.hpp"
#include "CommonIo.hpp" // chunked data reading
#include "RestartHashUtilities.hpp"
#include "PeriodicData.hpp"
#include "SvdNew.hpp"
#include "Histogram.hpp"
// added to read registered sles data 
#include "CtiRegister.hpp"
#include "DataContainer.hpp"
#include "Adt.hpp"

class Coeff {
public:
  uint8 rbi,rbi_nbr;
  double coeff[3];
  bool transpose;
  Coeff() {
    transpose = false;
  }
  Coeff(const uint8 rbi,const uint8 rbi_nbr,const double coeff[3]) {
    this->rbi = rbi;
    this->rbi_nbr = rbi_nbr;
    this->coeff[0] = coeff[0];
    this->coeff[1] = coeff[1];
    this->coeff[2] = coeff[2];
    transpose = false;
  }
  void copy(const Coeff& other) {
    rbi       = other.rbi;
    rbi_nbr   = other.rbi_nbr;
    coeff[0]  = other.coeff[0];
    coeff[1]  = other.coeff[1];
    coeff[2]  = other.coeff[2];
    transpose = other.transpose;
  }
  bool operator<(const Coeff& other) const { 
    return (rbi < other.rbi)||((rbi == other.rbi)&&(rbi_nbr < other.rbi_nbr)); 
  }
  // remove this??
  bool operator==(const Coeff& other) const { 
    return (rbi == other.rbi)&&(rbi_nbr == other.rbi_nbr); 
  }
};

struct rbiPairComparator { 
  bool operator()(const pair<uint8,uint8>& s, const pair<uint8,uint8>& other) const {
    return (s.first < other.first)||((s.first == other.first)&&(s.second < other.second)); 
  }
};

class StripedMesh {
  
public:

  int8 nno_global;
  int8 nbf_global;
  int8 nfa_global;
  int8 nfa_p_global;
  int8 ncv_global;
  int8 nno_p_global;
  int8 nno_pb_global;
  int8 nef_global;
  //int nbb_global;

  int8 nfa_zone[27];
  int8 noofa_zone[27];

  class BfZone {
  public:
    string name;
    int8 nbf_global;
    int8 noobf_global;
    int8 ibf_f_global;
    int8 nob_f_global;
    int nbf;
    int8* bfora;
    // surface
    int8 spobf_global;
    int8 spob_f_global;
    int8 sbobf_global;
    int8 stob_f_global; // goes with sb for now: TODO: change to st everywhere
    // geometry
    double area_global;
    double area_over_delta_global;
    double n_global[3];
    double x_global[3];
    BfZone() {
      bfora = NULL;
    }
    void clear() {
      DELETE(bfora);
    }
    // moved from staticsolver...
    void setName(const string& name) {

      this->name = name;

      // XXX this can eventually be replaced with assertions, but
      // we need to replace any "-+/" tokens in the zone name with
      // "_" in order to keep the expression evaluation from misinterpreting
      // the zone name as math.  not terribly efficient, but it will work for now

      std::replace(this->name.begin(), this->name.end(), '-', '_');
      std::replace(this->name.begin(), this->name.end(), '+', '_');
      std::replace(this->name.begin(), this->name.end(), '/', '_');
      std::replace(this->name.begin(), this->name.end(), '*', '_');
      std::replace(this->name.begin(), this->name.end(), '=', '_');
    }
  };
  
  vector<BfZone> bfZoneVec;
  map<const string,int> bfZoneNameMap;

  // put all this stuff somewhere useful eventually. The same thing is in 
  // stitch's Surface class...

  /*
    vector<PeriodicTransform> periodicTransformVec; // should this be in a namespace?
  
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

    bool getPeriodicR(double R[9],const int bits) const {
    
    // set I...
    R[0] = 1.0; R[1] = 0.0; R[2] = 0.0;
    R[3] = 0.0; R[4] = 1.0; R[5] = 0.0;
    R[6] = 0.0; R[7] = 0.0; R[8] = 1.0;

    // loop on bit pairs...
    bool has_R = false;
    FOR_I3 {
    if (bits & (1<<(2*i))) {
    // even bit pair set...
    assert(!(bits & (1<<(2*i+1)))); // only one pair set
    assert(periodicTransformVec.size() > i); // make sure transform vec supports this bit pair
    double this_R[9];
    if (periodicTransformVec[i].getR(this_R)) {
    has_R = true;
    double Rtmp[9] = { 0.0, 0.0, 0.0, 
    0.0, 0.0, 0.0, 
    0.0, 0.0, 0.0 };
    FOR_M3 {
    FOR_J3 {
    FOR_K3 {
    Rtmp[3*m+j] += R[3*m+k]*this_R[3*k+j];
    }
    }
    }
    for (int ii =0; ii < 9; ++ii) 
    R[ii] = Rtmp[ii];
    }
    }
    else if (bits & (1<<(2*i+1))) {
    // odd bit pair set...
    assert(periodicTransformVec.size() > i); // make sure transform vec supports this bit pair
    double this_R[9];
    if (periodicTransformVec[i].getInvR(this_R)) {
    has_R = true;
    double Rtmp[9] = { 0.0, 0.0, 0.0, 
    0.0, 0.0, 0.0, 
    0.0, 0.0, 0.0 };
    FOR_M3 {
    FOR_J3 {
    FOR_K3 {
    Rtmp[3*m+j] += R[3*m+k]*this_R[3*k+j];
    }
    }
    }
    for (int ii =0; ii < 9; ++ii) 
    R[ii] = Rtmp[ii];
    }
    }
    }
    return has_R;
    }

    bool getPeriodicT(double t[3],const int bits) const {

    // zero t...
    t[0] = 0.0; t[1] = 0.0; t[2] = 0.0;

    // loop on bit pairs...
    bool has_t = false;
    FOR_I3 {
    if (bits & (1<<(2*i))) {
    // even bit pair set...
    assert(!(bits & (1<<(2*i+1)))); // only one pair set
    assert(periodicTransformVec.size() > i); // make sure transform vec supports this bit pair
    double this_t[3];
    if (periodicTransformVec[i].getT(this_t)) {
    has_t = true;
    FOR_J3 t[j] += this_t[j];
    }
    }
    else if (bits & (1<<(2*i+1))) {
    // odd bit pair set...
    assert(periodicTransformVec.size() > i); // make sure transform vec supports this bit pair
    double this_t[3];
    if (periodicTransformVec[i].getInvT(this_t)) {
    has_t = true;
    FOR_J3 t[j] += this_t[j];
    }
    }
    }
    return has_t;
    }
  */
  
  int8 * bfora;
  int nbf;
  int8 * cvobf_global;
  int * noobf_i;
  int8 * noobf_v_global;
  int * zone_bf;
  
  int8 * faora;
  int nfa;
  int8 (*cvofa_global)[2];
  int * noofa_i;
  int8 * noofa_v_global;
  
  int8 * noora;
  int nno;
  double (*x_no)[3];

  uint8 * pbi_pno;
  int nno_p;

  int8 * cvora;
  int ncv;
  
  class CvZone {
  public:
    string name;
    int8 ncv_global;
    int icv_f_global;
    CvZone(const string& name) {
      this->name = name;
    }
  };

  vector<CvZone> cvZoneVec;

  double * vol_cv;
  double (*x_cv)[3];
  double * r_vv;
  double (*x_vv)[3];

  double (*n_bf)[3];
  double (*x_bf)[3];
  double * area_bf;
  double * area_over_delta_bf;
  double (*Gij_bf)[3][3];

  double (*n_fa)[3];
  double (*x_fa)[3];
  //double * area_over_delta_fa;
  //int * group_fa;
  
  int8 * efora;
  int nef;
  double (*n_ef)[3];
  double (*c_ef)[3];
  int8 (*cvoef_global)[2];
  //int * group_ef;
  
  //int * bbora;
  //int * ncvobb;
  //double (*vv_bbox)[6];
  //int nbb; 

  // surface stuff 
  
  bool b_surface;
  int surface_nsp_global,surface_nst_global; // int8 in future?
  int surface_nsp,surface_nst;
  int8 * surface_spora; // had to make int8 for dde to work
  int8 * surface_stora; // had to make int8 for dde to work
  double (*surface_xp)[3];
  int (*surface_spost)[3]; // int8 in future?
  int *surface_znost;
  int *surface_spobf_i;
  uint8 *surface_spobf_v_global;
  double *surface_spobf_wgt_global;
  int *surface_sbobf_i;
  uint8 *surface_sbobf_v_global; // ist | (bits<<52)
  
  //int surface_npp;
  //int surface_npp_global;
  //int *surface_isp_pp; // point associated to periodic point 
  //int *surface_isp_p_pp; // point parent of periodic point
  //uint2 *surface_isp_bits_pp; // point bits of periodic point
  //int8 *surface_ppora; // periodic point striping
  
  // lagrangian particles...

  class LpHelper {
  public:
    string name;
    int nlp;
    int8 nlp_global;
    int index;
    int8* lpocv_i_global;

    LpHelper() {
      name = "";
      lpocv_i_global = NULL;
      index = -1;
      nlp = -1;
      nlp_global = -1;
    }

    void clear() {
      DELETE(lpocv_i_global);
    }

  };
  vector<LpHelper> lpHelperVec;
  map<const string,int> lpHelperNameMap;

  // sles stuff...
  
  list<pair<CtiRegister::CtiData*,int> > iList; // int values
  list<pair<CtiRegister::CtiData*,double> > dList; // double values
  list<InData> inList; // int scalar fields
  list<DnData> dnList; // double scalar fields
  list<Dn3Data> dn3List; // double vector fields
  int8* cv_lb_wgt; // soln-based load balance wgt
  
  // io_version 3 stuff...
  int nfabf;
  int8 nfabf_global;
  int8* fabfora;
  int* zone_fabf;
  int8 my_nfa_disp;
  int8 *my_nbf_disp;
  vector<pair<int,int> > zone_ifabf_vec;

  StripedMesh() {
    
    // add others eventually...

    nef = 0;
    nef_global = 0;

    bfora = NULL;
    cvobf_global = NULL;
    noobf_i = NULL;
    noobf_v_global = NULL;
    zone_bf = NULL;
    
    faora = NULL;
    cvofa_global = NULL;
    noofa_i = NULL;
    noofa_v_global = NULL;

    noora = NULL;
    x_no = NULL;
    pbi_pno = NULL;
    
    cvora = NULL;
    vol_cv = NULL;
    x_cv = NULL;
    r_vv = NULL;
    x_vv = NULL;
    
    n_bf = NULL;
    x_bf = NULL;
    area_bf = NULL;
    area_over_delta_bf = NULL;
    Gij_bf = NULL;

    n_fa = NULL;
    x_fa = NULL;
    //area_over_delta_fa = NULL;
    //group_fa = NULL;

    efora = NULL;
    n_ef = NULL;
    c_ef = NULL;
    cvoef_global = NULL;
    //group_ef = NULL;
    
    // bbora = NULL;
    // ncvobb = NULL;
    // vv_bbox = NULL;

    // surface info...
    b_surface = false;
    surface_nsp_global = 0;
    surface_nst_global = 0;
    surface_nsp = 0;
    surface_nst = 0;
    surface_spora = NULL;
    surface_stora = NULL;
    surface_xp = NULL;
    surface_spost = NULL;
    surface_znost = NULL;
    surface_spobf_i = NULL;
    surface_spobf_v_global = NULL;
    surface_spobf_wgt_global = NULL;
    surface_sbobf_i = NULL;
    surface_sbobf_v_global = NULL;
    // periodic point info 
    //surface_npp = 0;
    //surface_npp_global = 0;
    //surface_isp_pp = NULL; 
    //surface_isp_p_pp = NULL; 
    //surface_isp_bits_pp = NULL; 
    //surface_ppora = NULL;; 
    
    cv_lb_wgt = NULL;

    fabfora = NULL;
    zone_fabf = NULL;
    my_nbf_disp = NULL;
    my_nfa_disp = -1;
    zone_ifabf_vec.clear();
  }
  
  ~StripedMesh() {

    DELETE(bfora);
    DELETE(cvobf_global);

    DELETE(noobf_i);
    DELETE(noobf_v_global);
    DELETE(zone_bf);
    
    DELETE(faora);
    DELETE(cvofa_global);

    DELETE(noofa_i);
    DELETE(noofa_v_global);
    
    DELETE(noora);
    DELETE(x_no);
    DELETE(pbi_pno);
    
    DELETE(cvora);
    
    DELETE(vol_cv);
    DELETE(x_cv);
    DELETE(r_vv);
    DELETE(x_vv);
    
    DELETE(n_bf);
    DELETE(x_bf);
    DELETE(area_bf);
    DELETE(area_over_delta_bf);
    DELETE(Gij_bf);

    //DELETE(group_fa);
    
    DELETE(n_fa);
    DELETE(x_fa);
    //DELETE(area_over_delta_fa);
    
    DELETE(efora);
    DELETE(n_ef);
    DELETE(c_ef);
    DELETE(cvoef_global);
    //DELETE(group_ef);

    // DELETE(bbora);
    // DELETE(ncvobb);
    // DELETE(vv_bbox);

    DELETE(surface_spora);
    DELETE(surface_stora);
    DELETE(surface_xp);
    DELETE(surface_spost);
    DELETE(surface_znost);
    DELETE(surface_spobf_i);
    DELETE(surface_spobf_v_global);
    DELETE(surface_spobf_wgt_global);
    DELETE(surface_sbobf_i);
    DELETE(surface_sbobf_v_global);
    //DELETE(surface_isp_pp);
    //DELETE(surface_isp_p_pp);
    //DELETE(surface_isp_bits_pp);
    //DELETE(surface_ppora);
    
    DELETE(cv_lb_wgt);

    FOR_IZONE(bfZoneVec) bfZoneVec[izone].clear();

    DELETE(fabfora);
    DELETE(zone_fabf);
    DELETE(my_nbf_disp);

    for (int ii = 0, lim = lpHelperVec.size(); ii < lim; ++ii)
      lpHelperVec[ii].clear();

  }

  void dumpSurface() {
    
    if (!b_surface) {
      if (mpi_rank == 0)
	cout << "StripedMesh::dumpSurface(): no surface" << endl;
      return;
    }
    
    assert(surface_stora);
    assert(surface_stora[mpi_size] == surface_nst_global);
    assert(surface_stora[mpi_rank+1]-surface_stora[mpi_rank] == surface_nst);
    
    assert(surface_spora);
    assert(surface_spora[mpi_size] == surface_nsp_global);
    assert(surface_spora[mpi_rank+1]-surface_spora[mpi_rank] == surface_nsp);
    
    if (mpi_rank == 0) {
      cout << "StripedMesh::dumpSurface(): surface_nst_global: " << surface_nst_global << " surface_nsp_global: " << surface_nsp_global << endl;
      FOR_RANK {
	cout << " > rank: " << rank << " nst: " << surface_stora[rank+1]-surface_stora[rank] << " nsp: " << surface_spora[rank+1]-surface_spora[rank] << endl;
      }
      cout << " > rank 0 xp[0:20): " << endl;
      for (int isp = 0; isp < min(20,surface_nsp); ++isp)
	cout << " >> isp: " << isp << " xp: " << COUT_VEC(surface_xp[isp]) << endl;
      cout << " > rank 0 spost[0:20), znost[0:20): " << endl;
      for (int ist = 0; ist < min(20,surface_nst); ++ist)
	cout << " >> ist: " << ist << " sposp: " << COUT_VEC(surface_spost[ist]) << " znost: " << surface_znost[ist] << endl;
    }

    // check indexing in surface_sbobf_i/v_global....
    
    for (int ibf = 0; ibf < nbf; ++ibf) {
      for (int stob = surface_sbobf_i[ibf]; stob != surface_sbobf_i[ibf+1]; ++stob) {
	const int ist_global = int(surface_sbobf_v_global[stob]&MASK_52BITS); // can be int8 some day
	assert((ist_global >= 0)&&(ist_global < surface_nst_global));
	const int bits = int(surface_sbobf_v_global[stob]>>52);
	assert((bits >= 0)&&(bits < (1<<6))); // actually (1<<12) is supported by the 52 bit shift (i.e. 64-52=12)
      }
    }

    // check indexing in surface_spobf_i/v/wgt...
    
    for (int ibf = 0; ibf < nbf; ++ibf) {
      double wgt_sum = 0.0;
      for (int spob = surface_spobf_i[ibf]; spob != surface_spobf_i[ibf+1]; ++spob) {
	const int isp_global = int(surface_spobf_v_global[spob]&MASK_52BITS); // can be int8 some day
	assert((isp_global >= 0)&&(isp_global < surface_nsp_global));
	const int bits = int(surface_spobf_v_global[spob]>>52);
	assert(bits == 0.0); // no periodicity for now!
        // wgts should be positive...
        assert(surface_spobf_wgt_global[spob] >= 0.0);
        wgt_sum += surface_spobf_wgt_global[spob];
      }
      // and sum to 1.0...
      assert(fabs(1.0-wgt_sum) < 1.0E-12);
    }
    
  }
  
  void read_mles(const string& filename) {
    
    int ierr = read_mles(filename,mpi_comm);
   
    if ( ierr != 0) { 

      CERR("read_mles returned error: " << ierr);
      
    }
    
  }

  int init_simple_restart(Param * param,MPI_Comm& mpi_comm) { 

    int ierr = 0;

    // conversion between 0-based (I,J,K) and 0-based global indexing...
#define GLOBAL_NODE(I,J,K) ((((int8)(K))*(ny+1)+((int8)(J)))*(nx+1)+((int8)(I)))
#define GLOBAL_CV(I,J,K) ((((int8)(K))*(ny)+((int8)(J)))*(nx)+((int8)(I)))

#define GLOBAL_NODE_X0(i) (x0+double(i)*dx)
#define GLOBAL_NODE_X1(j) (y0+double(j)*dy)
#define GLOBAL_NODE_X2(k) (z0+double(k)*dz)

    if ( param->getString(0) == "CART") { 

      const double x0 = param->getDouble(1);
      const double x1 = param->getDouble(2);
      const double y0 = param->getDouble(3);
      const double y1 = param->getDouble(4);
      const double z0 = param->getDouble(5);
      const double z1 = param->getDouble(6);
      const int8 nx   = param->getInt(7);
      const int8 ny   = param->getInt(8);
      const int8 nz   = param->getInt(9);

      int mpi_rank; MPI_Comm_rank(mpi_comm,&mpi_rank);
      int mpi_size; MPI_Comm_size(mpi_comm,&mpi_size);
 
      ncv_global      = nx*ny*nz;
      nno_global      = (nx+1)*(ny+1)*(nz+1);
      nfa_global      = (nx-1)*ny*nz + (ny-1)*nx*nz + (nz-1)*nx*ny;
      nbf_global      =  2*( ny*nz + nx*ny + nx*nz);

      // no periodicity -- use stitch if you want it.

      //nno_p_global    = 0;
      //nno_pb_global   = 8*1 + 4*((nx-1)+(ny-1)+(nz-1))+2*((nx-1)*(ny-1)+(nx-1)*(nz-1)+(ny-1)*(nz-1));

      nno_p_global    = -1;
      nno_pb_global   = -1;

      nfa_p_global    = 0;

      for (int i = 0; i < 26; ++i) { 
        //nfa_zone[i]   = 0;
        //noofa_zone[i] = 0;
        nfa_zone[i]     = -1;
        noofa_zone[i]   = -1;
      }

      //nfa_zone[26]   = nfa_global;
      //noofa_zone[26] = 4*nfa_global; 
      nfa_zone[26] = -1;
      noofa_zone[26] = -1;

      assert(bfora == NULL);
      MiscUtils::calcUniformDist(bfora,nbf_global,mpi_size); // mpi_size is the size of the passed comm
      assert(bfora[mpi_rank+1]-bfora[mpi_rank] < TWO_BILLION);
      nbf = bfora[mpi_rank+1]-bfora[mpi_rank];
      
      assert(faora == NULL);
      MiscUtils::calcUniformDist(faora,nfa_global,mpi_size);
      assert(faora[mpi_rank+1]-faora[mpi_rank] < TWO_BILLION);
      nfa = faora[mpi_rank+1]-faora[mpi_rank];
      
      assert(noora == NULL);
      MiscUtils::calcUniformDist(noora,nno_global,mpi_size);
      assert(noora[mpi_rank+1]-noora[mpi_rank] < TWO_BILLION);
      nno = noora[mpi_rank+1]-noora[mpi_rank];
      
      assert(cvora == NULL);
      MiscUtils::calcUniformDist(cvora,ncv_global,mpi_size);
      assert(cvora[mpi_rank+1]-cvora[mpi_rank] < TWO_BILLION);
      ncv = cvora[mpi_rank+1]-cvora[mpi_rank];
      
      // extended face build not part of stitch anymore...
      assert(efora == NULL);
     

      const double Ly  = y1-y0;
      const double Lz  = z1-z0;
      const double Lx  = x1-x0;

      // and geom...
      
      int8 nbf_global_check = 0;
 

      for (int i = 0; i < 6; ++i) { 

        double L1, L2, L3;
        double this_n[3] = {0.0,0.0,0.0};
        double this_x[3] = {(x0+x1)/2.0,(y0+y1)/2.0,(z0+z1)/2.0};
        int n1, n2, n3;
        string this_name;

        switch(i) { 
        case 0: 
          { 
            this_name = "x0";
            n1        = nx;
            n2        = ny;
            n3        = nz;
            L1        = Lx;
            L2        = Ly;
            L3        = Lz;

            this_n[0] = -Ly*Lz;
            this_x[0] = x0;
          }
          break;

        case 1:
          { 
            this_name = "x1";
            n1        = nx;
            n2        = ny;
            n3        = nz;
            L1        = Lx;
            L2        = Ly;
            L3        = Lz;

            this_n[0] = Ly*Lz;
            this_x[0] = x1;
          }
          break;

        case 2: 
          {
            this_name = "y0";
            n1        = ny;
            n2        = nx;
            n3        = nz;
            L1        = Ly;
            L2        = Lx;
            L3        = Lz;

            this_n[1] = -Lx*Lz;
            this_x[1] = y0;
          }
          break;

        case 3: 
          {
            this_name = "y1";
            n1        = ny;
            n2        = nx;
            n3        = nz;
            L1        = Ly;
            L2        = Lx;
            L3        = Lz;

            this_n[1] = Lx*Lz;
            this_x[1] = y1;
          }
          break;

        case 4: 
          {
            this_name = "z0";
            n1        = nz;
            n2        = nx;
            n3        = ny;
            L1        = Lz;
            L2        = Lx;
            L3        = Ly;

            this_n[2] = -Lx*Ly;
            this_x[2] = z0;
          }
          break;

        case 5: 
          {
            this_name = "z1";
            n1        = nz;
            n2        = nx;
            n3        = ny;
            L1        = Lz;
            L2        = Lx;
            L3        = Ly;

            this_n[2] = Lx*Ly;
            this_x[2] = z1;
          }
          break;

        default: 
          assert(0);
        }


        bfZoneVec.push_back(BfZone()); bfZoneVec.back().setName(this_name);
        bfZoneVec.back().nbf_global   = n2*n3; 
        bfZoneVec.back().noobf_global = 4*n2*n3; 
        bfZoneVec.back().ibf_f_global = nbf_global_check; 
        bfZoneVec.back().nob_f_global = 4*nbf_global_check;  
        bfZoneVec.back().sbobf_global = 0;
        bfZoneVec.back().stob_f_global = 0;
        
        bfZoneVec.back().area_global            = L2*L3; 
        bfZoneVec.back().area_over_delta_global = L2*L3/double(n2*n3)/(L1/double(n1)/2.0); 
        bfZoneVec.back().n_global[0]            = this_n[0];
        bfZoneVec.back().n_global[1]            = this_n[1]; 
        bfZoneVec.back().n_global[2]            = this_n[2]; 
        bfZoneVec.back().x_global[0]            = this_x[0]; 
        bfZoneVec.back().x_global[1]            = this_x[1];
        bfZoneVec.back().x_global[2]            = this_x[2]; 
        
        assert(bfZoneVec.back().bfora == NULL);
        MiscUtils::calcThresholdDist(bfZoneVec.back().bfora,bfZoneVec.back().nbf_global,mpi_size,DIST_THRESHOLD); 
        assert(bfZoneVec.back().bfora[mpi_rank+1]-bfZoneVec.back().bfora[mpi_rank] < TWO_BILLION);
        bfZoneVec.back().nbf = bfZoneVec.back().bfora[mpi_rank+1]-bfZoneVec.back().bfora[mpi_rank];
        
        if (mpi_rank == 0) {
          cout << " > bf zone \"" << bfZoneVec.back().name << "\", nbf_global: " << bfZoneVec.back().nbf_global << 
            ", area_global: " << bfZoneVec.back().area_global <<
            ", n_global: " << COUT_VEC(bfZoneVec.back().n_global) <<
            ", x_global: " << COUT_VEC(bfZoneVec.back().x_global) << endl;
        }
      

        nbf_global_check += bfZoneVec.back().nbf_global;
        
      } 

      assert( nbf_global_check == nbf_global);


      // the hard work starts now.  all the bf stuff..  
      
      assert( noobf_i        == NULL); noobf_i        = new int[nbf+1];
      assert( noobf_v_global == NULL); noobf_v_global = new int8[4*nbf];
      assert( cvobf_global   == NULL); cvobf_global   = new int8[nbf];
      assert( x_bf           == NULL); x_bf           = new double[nbf][3];
      assert( area_bf        == NULL); area_bf        = new double[nbf];
      assert( n_bf           == NULL); n_bf           = new double[nbf][3];
      assert( Gij_bf         == NULL); Gij_bf         = new double[nbf][3][3];
      assert( area_over_delta_bf == NULL); area_over_delta_bf = new double[nbf];
      assert( zone_bf        == NULL); zone_bf        = new int[nbf];
      
      int8 ibf_global = 0;
      noobf_i[0]      = 0;
      
      const double dx = Lx/double(nx);
      const double dy = Ly/double(ny);
      const double dz = Lz/double(nz);
      
      // x0... 
      
      for (int j = 0; j < ny; ++j) { 
        for (int k = 0; k < nz; ++k) { 
          
          if ( (ibf_global >= bfora[mpi_rank]) && (ibf_global < bfora[mpi_rank+1])) { 
            
            const int ibf  = ibf_global - bfora[mpi_rank];
            zone_bf[ibf] = 0;
            noobf_i[ibf+1] = 4*(ibf+1);
            noobf_v_global[4*ibf+3] = GLOBAL_NODE(0,j,k);
            noobf_v_global[4*ibf+2] = GLOBAL_NODE(0,j+1,k);
            noobf_v_global[4*ibf+1] = GLOBAL_NODE(0,j+1,k+1);
            noobf_v_global[4*ibf+0] = GLOBAL_NODE(0,j,k+1);
            
            cvobf_global[ibf]       = GLOBAL_CV(0,j,k);
            area_bf[ibf]            = dy*dz;
            area_over_delta_bf[ibf] = dy*dz/(dx/2.0);
            n_bf[ibf][0]            = -area_bf[ibf];
            n_bf[ibf][1]            = 0.0;
            n_bf[ibf][2]            = 0.0;
            
            x_bf[ibf][0] = GLOBAL_NODE_X0(0); 
            x_bf[ibf][1] = 0.5*(GLOBAL_NODE_X1(j) + GLOBAL_NODE_X1(j+1)); 
            x_bf[ibf][2] = 0.5*(GLOBAL_NODE_X2(k) + GLOBAL_NODE_X2(k+1));
            
            double dx_[3] = {-dx/2.0,0.0,0.0};
            for (int r = 0; r < 3; ++r) { 
              for (int s = 0; s < 3; ++s) { 
                Gij_bf[ibf][r][s] = -dx_[s]*n_bf[ibf][r];
              }
            }
            
          }
          
          ++ibf_global;
        }
      }
      
      // x1... 
      
      for (int j = 0; j < ny; ++j) { 
        for (int k = 0; k < nz; ++k) { 
          
          if ( (ibf_global >= bfora[mpi_rank]) && (ibf_global < bfora[mpi_rank+1])) { 
            
            const int ibf  = ibf_global - bfora[mpi_rank];
            zone_bf[ibf] = 1;
            noobf_i[ibf+1] = 4*(ibf+1);
            noobf_v_global[4*ibf+0] = GLOBAL_NODE(nx,j,k);
            noobf_v_global[4*ibf+1] = GLOBAL_NODE(nx,j+1,k);
            noobf_v_global[4*ibf+2] = GLOBAL_NODE(nx,j+1,k+1);
            noobf_v_global[4*ibf+3] = GLOBAL_NODE(nx,j,k+1);
            
            cvobf_global[ibf]       = GLOBAL_CV(nx-1,j,k);
            area_bf[ibf]            = dy*dz;
            area_over_delta_bf[ibf] = dy*dz/(dx/2.0);
            n_bf[ibf][0]            = area_bf[ibf];
            n_bf[ibf][1]            = 0.0;
            n_bf[ibf][2]            = 0.0;
            
            x_bf[ibf][0] = GLOBAL_NODE_X0(nx); 
            x_bf[ibf][1] = 0.5*(GLOBAL_NODE_X1(j) + GLOBAL_NODE_X1(j+1)); 
            x_bf[ibf][2] = 0.5*(GLOBAL_NODE_X2(k) + GLOBAL_NODE_X2(k+1));
            
            double dx_[3] = {dx/2.0,0.0,0.0};
            for (int r = 0; r < 3; ++r) { 
              for (int s = 0; s < 3; ++s) { 
                Gij_bf[ibf][r][s] = -dx_[s]*n_bf[ibf][r];
              }
            }
            
          }
          
          ++ibf_global;
        }
      }
      
      
      // y0... 
      
      for (int i = 0; i < nx; ++i) { 
        for (int k = 0; k < nz; ++k) { 
          
          if ( (ibf_global >= bfora[mpi_rank]) && (ibf_global < bfora[mpi_rank+1])) { 
            
            const int ibf  = ibf_global - bfora[mpi_rank];
            zone_bf[ibf] = 2;
            noobf_i[ibf+1] = 4*(ibf+1);
            noobf_v_global[4*ibf+0] = GLOBAL_NODE(i,0,k);
            noobf_v_global[4*ibf+1] = GLOBAL_NODE(i+1,0,k);
            noobf_v_global[4*ibf+2] = GLOBAL_NODE(i+1,0,k+1);
            noobf_v_global[4*ibf+3] = GLOBAL_NODE(i,0,k+1);
            
            cvobf_global[ibf]       = GLOBAL_CV(i,0,k);
            area_bf[ibf]            = dx*dz;
            area_over_delta_bf[ibf] = dx*dz/(dy/2.0);
            n_bf[ibf][1]            = -area_bf[ibf];
            n_bf[ibf][2]            = 0.0;
            n_bf[ibf][0]            = 0.0;
            
            x_bf[ibf][1] = GLOBAL_NODE_X1(0); 
            x_bf[ibf][2] = 0.5*(GLOBAL_NODE_X2(k) + GLOBAL_NODE_X2(k+1));
            x_bf[ibf][0] = 0.5*(GLOBAL_NODE_X0(i) + GLOBAL_NODE_X0(i+1)); 
            
            double dx_[3] = {0.0,-dy/2.0,0.0};
            for (int r = 0; r < 3; ++r) { 
              for (int s = 0; s < 3; ++s) { 
                Gij_bf[ibf][r][s] = -dx_[s]*n_bf[ibf][r];
              }
            }
            
          }
          
          ++ibf_global;
        }
      }
      
      // y1... 
      
      for (int i = 0; i < nx; ++i) { 
        for (int k = 0; k < nz; ++k) { 
          
          if ( (ibf_global >= bfora[mpi_rank]) && (ibf_global < bfora[mpi_rank+1])) { 
            
            const int ibf  = ibf_global - bfora[mpi_rank];
            zone_bf[ibf] = 3;
            noobf_i[ibf+1] = 4*(ibf+1);
            noobf_v_global[4*ibf+3] = GLOBAL_NODE(i,ny,k);
            noobf_v_global[4*ibf+2] = GLOBAL_NODE(i+1,ny,k);
            noobf_v_global[4*ibf+1] = GLOBAL_NODE(i+1,ny,k+1);
            noobf_v_global[4*ibf+0] = GLOBAL_NODE(i,ny,k+1);
            
            cvobf_global[ibf]       = GLOBAL_CV(i,ny-1,k);
            area_bf[ibf]            = dx*dz;
            area_over_delta_bf[ibf] = dx*dz/(dy/2.0);
            n_bf[ibf][1]            = area_bf[ibf];
            n_bf[ibf][2]            = 0.0;
            n_bf[ibf][0]            = 0.0;
            
            x_bf[ibf][1] = GLOBAL_NODE_X1(ny); 
            x_bf[ibf][2] = 0.5*(GLOBAL_NODE_X2(k) + GLOBAL_NODE_X2(k+1));
            x_bf[ibf][0] = 0.5*(GLOBAL_NODE_X0(i) + GLOBAL_NODE_X0(i+1)); 
            
            double dx_[3] = {0.0,dy/2.0,0.0};
            for (int r = 0; r < 3; ++r) { 
              for (int s = 0; s < 3; ++s) { 
                Gij_bf[ibf][r][s] = -dx_[s]*n_bf[ibf][r];
              }
            }
            
          }
          
          ++ibf_global;
        }
      }
      
      // z0... 
      
      for (int i = 0; i < nx; ++i) { 
        for (int j = 0; j < ny; ++j) { 
          
          if ( (ibf_global >= bfora[mpi_rank]) && (ibf_global < bfora[mpi_rank+1])) { 
            
            const int ibf  = ibf_global - bfora[mpi_rank];
            zone_bf[ibf] = 4;
            noobf_i[ibf+1] = 4*(ibf+1);
            noobf_v_global[4*ibf+3] = GLOBAL_NODE(i,j,0);
            noobf_v_global[4*ibf+2] = GLOBAL_NODE(i+1,j,0);
            noobf_v_global[4*ibf+1] = GLOBAL_NODE(i+1,j+1,0);
            noobf_v_global[4*ibf+0] = GLOBAL_NODE(i,j+1,0);
            
            cvobf_global[ibf]       = GLOBAL_CV(i,j,0);
            area_bf[ibf]            = dx*dy;
            area_over_delta_bf[ibf] = dx*dy/(dz/2.0);
            n_bf[ibf][2]            = -area_bf[ibf];
            n_bf[ibf][0]            = 0.0;
            n_bf[ibf][1]            = 0.0;
            
            x_bf[ibf][2] = GLOBAL_NODE_X2(0); 
            x_bf[ibf][0] = 0.5*(GLOBAL_NODE_X0(i) + GLOBAL_NODE_X0(i+1)); 
            x_bf[ibf][1] = 0.5*(GLOBAL_NODE_X1(j) + GLOBAL_NODE_X1(j+1));
            
            double dx_[3] = {0.0,0.0,-dz/2.0};
            for (int r = 0; r < 3; ++r) { 
              for (int s = 0; s < 3; ++s) { 
                Gij_bf[ibf][r][s] = -dx_[s]*n_bf[ibf][r];
              }
            }
            
          }
          
          ++ibf_global;
        }
      }
      
      // z1... 
      
      for (int i = 0; i < nx; ++i) { 
        for (int j = 0; j < ny; ++j) { 
          
          if ( (ibf_global >= bfora[mpi_rank]) && (ibf_global < bfora[mpi_rank+1])) { 
            
            const int ibf  = ibf_global - bfora[mpi_rank];
            zone_bf[ibf] = 5;
            noobf_i[ibf+1] = 4*(ibf+1);
            noobf_v_global[4*ibf+0] = GLOBAL_NODE(i,j,nz);
            noobf_v_global[4*ibf+1] = GLOBAL_NODE(i+1,j,nz);
            noobf_v_global[4*ibf+2] = GLOBAL_NODE(i+1,j+1,nz);
            noobf_v_global[4*ibf+3] = GLOBAL_NODE(i,j+1,nz);
            
            cvobf_global[ibf]       = GLOBAL_CV(i,j,nz-1);
            area_bf[ibf]            = dx*dy;
            area_over_delta_bf[ibf] = dx*dy/(dz/2.0);
            n_bf[ibf][2]            = area_bf[ibf];
            n_bf[ibf][0]            = 0.0;
            n_bf[ibf][1]            = 0.0;
            
            x_bf[ibf][2] = GLOBAL_NODE_X2(nz); 
            x_bf[ibf][0] = 0.5*(GLOBAL_NODE_X0(i) + GLOBAL_NODE_X0(i+1)); 
            x_bf[ibf][1] = 0.5*(GLOBAL_NODE_X1(j) + GLOBAL_NODE_X1(j+1));
            
            double dx_[3] = {0.0,0.0,dz/2.0};
            for (int r = 0; r < 3; ++r) { 
              for (int s = 0; s < 3; ++s) { 
                Gij_bf[ibf][r][s] = -dx_[s]*n_bf[ibf][r];
              }
            }
            
          }
          
          ++ibf_global;
        }
      }
      
      
      assert( ibf_global == nbf_global);
      
      assert( noofa_i        == NULL); noofa_i        = new int[nfa+1];
      assert( noofa_v_global == NULL); noofa_v_global = new int8[4*nfa];
      assert( cvofa_global   == NULL); cvofa_global   = new int8[nfa][2];
      assert( x_fa           == NULL); x_fa           = new double[nfa][3];
      assert( n_fa           == NULL); n_fa           = new double[nfa][3];
      //assert( area_over_delta_fa == NULL); area_over_delta_fa = new double[nfa];
      
      int8 ifa_global = 0;
      noofa_i[0]      = 0;
      
      
      // internal x-faces...
      
      for (int i = 1; i < nx; ++i) { 
        for (int j = 0; j < ny; ++j) { 
          for (int k = 0; k < nz; ++k) {
            
            
            if ((ifa_global >= faora[mpi_rank])&&(ifa_global < faora[mpi_rank+1])) {
              
              const int ifa  = ifa_global - faora[mpi_rank];
              noofa_i[ifa+1] = 4*(ifa+1);
              
              noofa_v_global[4*ifa+0] = GLOBAL_NODE(i,j,k);
              noofa_v_global[4*ifa+1] = GLOBAL_NODE(i,j+1,k);
              noofa_v_global[4*ifa+2] = GLOBAL_NODE(i,j+1,k+1);
              noofa_v_global[4*ifa+3] = GLOBAL_NODE(i,j,k+1);
              
              cvofa_global[ifa][0] = GLOBAL_CV(i-1,j,k);
              cvofa_global[ifa][1] = GLOBAL_CV(i,j,k);
              
              x_fa[ifa][0]         = GLOBAL_NODE_X0(i);
              x_fa[ifa][1]         = 0.5*(GLOBAL_NODE_X1(j) + GLOBAL_NODE_X1(j+1));
              x_fa[ifa][2]         = 0.5*(GLOBAL_NODE_X2(k) + GLOBAL_NODE_X2(k+1));
              
              n_fa[ifa][0]         = dy*dz;
              n_fa[ifa][1]         = 0.0;
              n_fa[ifa][2]         = 0.0;
              
            }
            
            ++ifa_global;
            
          }
        }
      }
      
      
      for (int j = 1; j < ny; ++j) { 
        for (int k = 0; k < nz; ++k) { 
          for (int i = 0; i < nx; ++i) {
            
            if ((ifa_global >= faora[mpi_rank])&&(ifa_global < faora[mpi_rank+1])) {
              
              const int ifa  = ifa_global - faora[mpi_rank];
              noofa_i[ifa+1] = 4*(ifa+1); 
              
              noofa_v_global[4*ifa+0] = GLOBAL_NODE(i,j,k);
              noofa_v_global[4*ifa+1] = GLOBAL_NODE(i,j,k+1);
              noofa_v_global[4*ifa+2] = GLOBAL_NODE(i+1,j,k+1);
              noofa_v_global[4*ifa+3] = GLOBAL_NODE(i+1,j,k);
              
              cvofa_global[ifa][0] = GLOBAL_CV(i,j-1,k);
              cvofa_global[ifa][1] = GLOBAL_CV(i,j,k);
              
              x_fa[ifa][1]         = GLOBAL_NODE_X1(j);
              x_fa[ifa][2]         = 0.5*(GLOBAL_NODE_X2(k) + GLOBAL_NODE_X2(k+1));
              x_fa[ifa][0]         = 0.5*(GLOBAL_NODE_X0(i) + GLOBAL_NODE_X0(i+1));
              
              n_fa[ifa][1]         = dx*dz;
              n_fa[ifa][2]         = 0.0;
              n_fa[ifa][0]         = 0.0;
              
            }
            
            ++ifa_global;
          }
        }
      }
      
      
      // internal z-faces...
      
      for (int k = 1; k < nz; ++k) { 
        for (int i = 0; i < nx; ++i) { 
          for (int j = 0; j < ny; ++j) {
            
            if ((ifa_global >= faora[mpi_rank])&&(ifa_global < faora[mpi_rank+1])) {
              
              const int ifa  = ifa_global - faora[mpi_rank];
              noofa_i[ifa+1] = 4*(ifa+1);
              
              noofa_v_global[4*ifa+0] = GLOBAL_NODE(i,j,k);
              noofa_v_global[4*ifa+1] = GLOBAL_NODE(i+1,j,k);
              noofa_v_global[4*ifa+2] = GLOBAL_NODE(i+1,j+1,k);
              noofa_v_global[4*ifa+3] = GLOBAL_NODE(i,j+1,k);
              
              cvofa_global[ifa][0] = GLOBAL_CV(i,j,k-1);
              cvofa_global[ifa][1] = GLOBAL_CV(i,j,k);
              
              x_fa[ifa][2]         = GLOBAL_NODE_X2(k);
              x_fa[ifa][0]         = 0.5*(GLOBAL_NODE_X0(i) + GLOBAL_NODE_X0(i+1));
              x_fa[ifa][1]         = 0.5*(GLOBAL_NODE_X1(j) + GLOBAL_NODE_X1(j+1));
              
              n_fa[ifa][2]         = dx*dy;
              n_fa[ifa][0]         = 0.0;
              n_fa[ifa][1]         = 0.0;
              
            }
            
            ++ifa_global;
            
          }
        }
      }
      
      
      assert( ifa_global == nfa_global);
      
      
      assert( x_cv == NULL); x_cv = new double[ncv][3];
      assert( x_vv == NULL); x_vv = new double[ncv][3];
      assert( r_vv == NULL); r_vv = new double[ncv];
      assert( vol_cv == NULL); vol_cv = new double[ncv];
      
      int8 icv_global = 0;
      
      for (int k = 0; k < nz; ++k) { 
        for (int j = 0; j < ny; ++j) { 
          for (int i = 0; i < nx; ++i) { 
            
            if ( (icv_global >= cvora[mpi_rank]) && (icv_global < cvora[mpi_rank+1])) { 
             
              const int icv = icv_global - cvora[mpi_rank];

              x_cv[icv][0] = 0.5*(GLOBAL_NODE_X0(i) + GLOBAL_NODE_X0(i+1));
              x_cv[icv][1] = 0.5*(GLOBAL_NODE_X1(j) + GLOBAL_NODE_X1(j+1));
              x_cv[icv][2] = 0.5*(GLOBAL_NODE_X2(k) + GLOBAL_NODE_X2(k+1));
              
              for (int r = 0; r < 3; ++r) 
                x_vv[icv][r] = x_cv[icv][r];
              
              vol_cv[icv]  = dx*dy*dz;
              r_vv[icv]    = sqrt(dx*dx + dy*dy + dz*dz)/2.0;
              
            }
            
            ++icv_global;
          }
        }
      }
      
      assert( icv_global == ncv_global);
      
      assert ( x_no == NULL); x_no = new double[nno][3];
      
      int8 ino_global = 0;
      
      for (int k = 0; k <= nz; ++k) { 
        for (int j = 0; j <= ny; ++j) { 
          for (int i = 0; i <= nx; ++i) { 
            
            if ( (ino_global >= noora[mpi_rank]) && (ino_global < noora[mpi_rank+1])) { 
             
              const int ino = ino_global - noora[mpi_rank];

              x_no[ino][0] = GLOBAL_NODE_X0(i);
              x_no[ino][1] = GLOBAL_NODE_X1(j);
              x_no[ino][2] = GLOBAL_NODE_X2(k);
              
            }
            
            ++ino_global;
          }
        }
      }
      
      assert (ino_global == nno_global);
      
    } else { 
  
      ierr = -2;
      if (mpi_rank == 0) 
        cout << " > unrecognized simple restart type : " << param->getString(2) << endl;
      
    }
#undef GLOBAL_NODE
#undef GLOBAL_CV
    
    { //set a default mles hash if one was not found
      std::stringstream ss;
      ss << "Missing mles Hash";
      RestartHashUtilities::mlesHash.init(ss,RestartHashUtilities::sha1hashlength);
      if (mpi_rank==0) cout << " > setting default mles hash id " <<  RestartHashUtilities::mlesHash << endl;
    }

    return ierr;
  }

  int read_mles(const string& filename, MPI_Comm& mpi_comm) { 

    int ierr = read_mles_(filename,mpi_comm);
   
    if ( ierr != 0) { 

      if ( ierr == -1) { 
        CERR( " read_mles unable to find restart or does not behave expectedly: "  << filename );
      } else { 
        CERR("read_mles returned error: " << ierr);
      }
      
    }

    return ierr; 

  } 

  int read_mles_(const string& filename,MPI_Comm& mpi_comm) {
  
    // here we use integer error handling to allow the calling process to pass
    // a communicator that is different from the solver communicator -- e.g.
    // something more matched to the io node configuration...
  
    int mpi_rank; MPI_Comm_rank(mpi_comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(mpi_comm,&mpi_size);
    
    if (mpi_rank == 0) 
      cout << "StripedMesh::read_mles " << filename << " on " << mpi_size << " ranks..." << endl; 
   
    bool b_foundHash = false;

    MPI_Barrier(mpi_comm);
    double wtime0 = 0.0; 
    if (mpi_rank == 0)
      wtime0 = MPI_Wtime();
    
    MPI_File fh; 
    char dummy[128];
    sprintf(dummy,"%s",filename.c_str());
    int ierr = MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
    
    if (ierr != 0) { 
      //if (mpi_rank == 0) cout << "Error: cannot open restart file: " << filename << endl;
      return -1;
    } 
    
    int itmp[2] = { 0, 0 };
    if (mpi_rank == 0) { 
      // first 2 ints are: 0. magic number, 1. io version
      MPI_File_read(fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
    }
    MPI_Bcast(itmp,2,MPI_INT,0,mpi_comm); 

    bool byte_swap = false;
    if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
      ByteSwap::byteSwap(itmp,2);
      if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
	if (mpi_rank == 0) cout << "Error: restart file does not start as expected." << endl;
	return -1;
      }
      if (mpi_rank == 0) cout << " > file requires byte swapping." << endl;
      byte_swap = true;
    }
    
    int io_version = itmp[1];
    if ((io_version != 3)&&(io_version != 5)) {
      if (mpi_rank == 0) cout << "Error: restart file version not 3 or 5: " << io_version << endl;
      return -1;
    }

    MPI_Offset offset = 8; // 2 ints
    Header header; 
    int done = 0;
    int8 nbf_global_check = 0;
    int8 ncv_global_check = 0;
    set<string> nameSet;
    while (done != 1) {
      
      if (mpi_rank == 0) {
	MPI_File_read_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE); 
	if (byte_swap) ByteSwap::byteSwapHeader(&header,1);
      }
      MPI_Bcast(&header,1,MPI_Header,0,mpi_comm); 
      
      switch (header.id) {
	
      case UGP_IO_NO_FA_CV_COUNTS:

        assert(io_version == 3);

	nno_global = ByteSwap::getInt8FromLswMswPair(header.idata+0);
	nfa_global = ByteSwap::getInt8FromLswMswPair(header.idata+2);
	ncv_global = ByteSwap::getInt8FromLswMswPair(header.idata+4);
        nno_p_global = 0; // no periodic
	nfa_p_global = 0;

	if (mpi_rank == 0) {
	  cout << " > nno_global: " << nno_global << " nfa_global: " << nfa_global << " ncv_global: " << ncv_global << endl;
	  cout << " > nfa_global/ncv_global: " << double(nfa_global)/double(ncv_global) << endl;
 	}
	
        // note that these include boundary faces...
	assert(faora == NULL);
        MiscUtils::calcUniformDist(faora,nfa_global,mpi_size);
        assert(faora[mpi_rank+1]-faora[mpi_rank] < TWO_BILLION);
	nfa = faora[mpi_rank+1]-faora[mpi_rank];
	
	assert(noora == NULL);
        MiscUtils::calcUniformDist(noora,nno_global,mpi_size);
        assert(noora[mpi_rank+1]-noora[mpi_rank] < TWO_BILLION);
	nno = noora[mpi_rank+1]-noora[mpi_rank];
	
	assert(cvora == NULL);
        MiscUtils::calcUniformDist(cvora,ncv_global,mpi_size);
        assert(cvora[mpi_rank+1]-cvora[mpi_rank] < TWO_BILLION);
	ncv = cvora[mpi_rank+1]-cvora[mpi_rank];
        
        break;

      case UGP_IO_FA_CHECK:
        {

          assert(io_version == 3);

          if (mpi_rank == 0) 
            cout << " > face check..." << endl;
          
          assert(nfa_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));
          assert(faora != NULL);
          
          int *fa_flag = new int[nfa];
          readChunkedData<int>(fh,offset+header_size+faora[mpi_rank]*int_size,fa_flag,nfa,byte_swap,mpi_comm);
          
          FOR_IFA assert((ifa + faora[mpi_rank]) == fa_flag[ifa]);

          delete[] fa_flag;
        }

        break;

      case UGP_IO_FA_ZONE_HEADER: 

        assert(io_version == 3);

	if (mpi_rank == 0) 
	  cout << " > boundary and internal face zone header..." << endl;

        if (header.idata[0] == FA_ZONE_BOUNDARY) { 
          assert(header.idata[1] == bfZoneVec.size());

          // ensure the name is unique!...
          {
            string name = header.name;
            while (nameSet.find(name) != nameSet.end())
              name += "1";
            nameSet.insert(name);
            bfZoneVec.push_back(BfZone());
            bfZoneVec.back().setName(name);
          }

          bfZoneVec.back().nbf_global = ByteSwap::getInt8FromLswMswPair(header.idata+2);
          // check...
          nbf_global_check += bfZoneVec.back().nbf_global;
          // striping...
          assert(bfZoneVec.back().bfora == NULL);
          MiscUtils::calcThresholdDist(bfZoneVec.back().bfora,bfZoneVec.back().nbf_global,mpi_size,DIST_THRESHOLD); 
          assert(bfZoneVec.back().bfora[mpi_rank+1]-bfZoneVec.back().bfora[mpi_rank] < TWO_BILLION);
          bfZoneVec.back().nbf = bfZoneVec.back().bfora[mpi_rank+1]-bfZoneVec.back().bfora[mpi_rank];

          if (mpi_rank == 0) 
            cout << " > bf zone \"" << bfZoneVec.back().name << "\", nbf_global: " << bfZoneVec.back().nbf_global << endl;

        }
        else {
          assert(header.idata[0] == FA_ZONE_INTERNAL);
          assert(header.idata[1] >= bfZoneVec.size()); // internal zones came after boundaries (no periodicity for now)
	  assert(strcmp(header.name,"default-internal") == 0);
          const int8 nfa_global_internal = ByteSwap::getInt8FromLswMswPair(header.idata+2);
          
          // set nbf striping...
          
          nbf_global = nfa_global-nfa_global_internal;
          assert(bfora == NULL);
          MiscUtils::calcUniformDist(bfora,nbf_global,mpi_size); // mpi_size is the size of the passed comm
          assert(bfora[mpi_rank+1]-bfora[mpi_rank] < TWO_BILLION);
          nbf = bfora[mpi_rank+1]-bfora[mpi_rank];

          // no periodic...
          for (int i = 0; i < 26; ++i)
            nfa_zone[i] = 0;
          nfa_zone[26] = nfa_global_internal; // nfa_global will be set to this once we read everything in...

        }

        break;

      case UGP_IO_FA_ZONE:

        if (mpi_rank == 0) 
          cout << " > boundary and internal face zone..." << endl;
        
        assert(nfa_global == header.idata[0]);
        assert(zone_bf == NULL); // throw it into zone_bf for now
        assert(faora != NULL);
        
        zone_bf = new int[nfa];
        readChunkedData<int>(fh,offset+header_size+faora[mpi_rank]*int_size,zone_bf,nfa,byte_swap,mpi_comm);
	
	break;

      case UGP_IO_NOOFA_I_AND_V:

        assert(io_version == 3);

	if (mpi_rank == 0) 
	  cout << " > noofa_i and noofa_v_global..." << endl;
	
	assert(nfa_global == header.idata[0]);
	assert(noofa_i == NULL);
	assert(noofa_v_global == NULL);
	assert(faora != NULL);
	
	// recall noofa_i stores the count, so it is a simple int...

	noofa_i = new int[nfa+1];
        readChunkedData<int>(fh,offset+header_size+faora[mpi_rank]*int_size,noofa_i+1,nfa,byte_swap,mpi_comm);
	
	{
	  
	  // it is possible that noofa_i could exceed the int limit when turned into a CSR -- check this...
	  // if this happens, run on more nodes...

	  int8 noofa_s = 0;
	  for (int ifa = 0; ifa < nfa; ++ifa)
	    noofa_s += noofa_i[ifa+1];
	  assert(noofa_s < TWO_BILLION);
	  
	  noofa_i[0] = 0;
	  for (int ifa = 0; ifa < nfa; ++ifa) 
	    noofa_i[ifa+1] += noofa_i[ifa];
	  assert(noofa_s == noofa_i[nfa]);
	  
	  int8 noofa_offset;
	  MPI_Scan(&noofa_s,&noofa_offset,1,MPI_INT8,MPI_SUM,mpi_comm);
	  // the offset on the last rank will be the total -- bcast back to everyone...
	  
	  int8 noofa_s_global = noofa_offset;
	  MPI_Bcast(&noofa_s_global,1,MPI_INT8,mpi_size-1,mpi_comm);
	  noofa_offset -= noofa_s;
	  assert(noofa_s_global == ByteSwap::getInt8FromLswMswPair(header.idata+1));
	  
	  noofa_v_global = new int8[noofa_s];
          int* noofa_v = new int[noofa_s];
          readChunkedData<int>(fh,offset+header_size+faora[mpi_size]*int_size+noofa_offset*int_size,noofa_v,noofa_s,byte_swap,mpi_comm);
          for (int nof = 0; nof < noofa_s; ++nof)
            noofa_v_global[nof] = noofa_v[nof];
          delete[] noofa_v;

	}
	
	break;

      case UGP_IO_CVOFA:
        {

          if (mpi_rank == 0) 
            cout << " > cvofa_global..." << endl;
          
          assert(nfa_global == header.idata[0]);
          assert(header.idata[1] == 2);
          assert(cvofa_global == NULL);
          assert(faora != NULL);
          
          cvofa_global = new int8[nfa][2];
          int (*cvofa)[2] = new int[nfa][2];
          readChunkedData<int>(fh,offset+header_size+faora[mpi_rank]*int_size*2,(int*)cvofa,nfa*2,byte_swap,mpi_comm);
          for (int ifa = 0; ifa < nfa; ++ifa) {
            cvofa_global[ifa][0] = cvofa[ifa][0];
            cvofa_global[ifa][1] = cvofa[ifa][1];
          }
          delete[] cvofa;

        }
	
	break;

      case UGP_IO_NO_FA_BF_CV_COUNTS:
	
	nno_global    = ByteSwap::getInt8FromLswMswPair(header.idata+0);
	nbf_global    = ByteSwap::getInt8FromLswMswPair(header.idata+2);
	nfa_global    = ByteSwap::getInt8FromLswMswPair(header.idata+4);
	ncv_global    = ByteSwap::getInt8FromLswMswPair(header.idata+6);
	nno_p_global  = ByteSwap::getInt8FromLswMswPair(header.idata+8); // first nodes are periodic
	nno_pb_global = ByteSwap::getInt8FromLswMswPair(header.idata+10); // and boundary

	assert(b_surface == false);
	
	// skip surface hack XXXXX
	//if (header.idata[12] == 1) if (mpi_rank == 0) cout << "WARNING: forcing b_surface = false" << endl;
	if (header.idata[12] == 1)
          b_surface = true;

	if (mpi_rank == 0) {
	  MPI_File_read_at(fh,offset+header_size,nfa_zone,27,MPI_INT8,MPI_STATUS_IGNORE); 
	  if (byte_swap) ByteSwap::byteSwap(nfa_zone,27);
	}
	MPI_Bcast(nfa_zone,27,MPI_INT8,0,mpi_comm); 
	
	// recall zones are indexed as follows...
	// 0:25 - periodic of various flavors,
	// 26   - internal 
	// use the nfa_zone to compute nfa_p_global: the periodic face count...
	nfa_p_global = 0;
	for (int izone = 0; izone <= 25; ++izone) 
	  nfa_p_global += nfa_zone[izone];
	assert(nfa_p_global + nfa_zone[26] == nfa_global);
	
	if (mpi_rank == 0) {
	  MPI_File_read_at(fh,offset+header_size+int8_size*27,noofa_zone,27,MPI_INT8,MPI_STATUS_IGNORE); 
	  if (byte_swap) ByteSwap::byteSwap(noofa_zone,27);
	}
	MPI_Bcast(noofa_zone,27,MPI_INT8,0,mpi_comm); 

	if (mpi_rank == 0) {
	  cout << " > nno_global: " << nno_global << 
	    " nbf_global: " << nbf_global << 
	    " nfa_global: " << nfa_global << 
	    " ncv_global: " << ncv_global << endl;
	  cout << " > nfa_global/ncv_global: " << double(nfa_global)/double(ncv_global) << endl;
	  cout << " > b_surface: " << b_surface << endl;
 	}
	
	assert(bfora == NULL);
        MiscUtils::calcUniformDist(bfora,nbf_global,mpi_size); // mpi_size is the size of the passed comm
        assert(bfora[mpi_rank+1]-bfora[mpi_rank] < TWO_BILLION);
	nbf = bfora[mpi_rank+1]-bfora[mpi_rank];
	
	assert(faora == NULL);
        MiscUtils::calcUniformDist(faora,nfa_global,mpi_size);
        assert(faora[mpi_rank+1]-faora[mpi_rank] < TWO_BILLION);
	nfa = faora[mpi_rank+1]-faora[mpi_rank];
	
	assert(noora == NULL);
        MiscUtils::calcUniformDist(noora,nno_global,mpi_size);
        assert(noora[mpi_rank+1]-noora[mpi_rank] < TWO_BILLION);
	nno = noora[mpi_rank+1]-noora[mpi_rank];
	
	assert(cvora == NULL);
        MiscUtils::calcUniformDist(cvora,ncv_global,mpi_size);
        assert(cvora[mpi_rank+1]-cvora[mpi_rank] < TWO_BILLION);
	ncv = cvora[mpi_rank+1]-cvora[mpi_rank];

	// extended face build not part of stitch anymore...
	assert(efora == NULL);
	
	break;

      case UGP_IO_PERIODIC_TRANSFORM:
	
	// actually contains all periodic transforms in the header (1,2, or 3)...
	
	{
	  const int npt = header.idata[0];
	  assert((npt >= 1)&&(npt <= 3)); // if this header is here, should be atleast 1
	  if (mpi_rank == 0) cout << " > periodic transform(s): " << npt << endl;
	  assert(PeriodicData::periodicTransformVec.empty());
          PeriodicData::periodicTransformVec.resize(npt);
	  for (int ipt = 0; ipt < npt; ++ipt) {
            PeriodicData::periodicTransformVec[ipt].setKindAndData(header.idata[1+ipt],header.rdata+3*ipt);
	    if (mpi_rank == 0) {
	      cout << " > bit pair " << ipt << " "; 
              PeriodicData::periodicTransformVec[ipt].dump();
	    }
	  }
          
          int bits = 0;
          for (int ipt = 0; ipt < npt; ++ipt) 
            bits |= (1<<(2*ipt));
          double R[9],t[3];
          const bool has_R = PeriodicData::getPeriodicR(R,bits);
          const bool has_t = PeriodicData::getPeriodicT(t,bits);
          if (has_R&has_t) {
            // check that the rotation is orthogonal to the translation...
            double Rt[3];
	    MiscUtils::matVecMult(Rt,R,t);
            if (mpi_rank == 0)
              cout << " > orthogonality check on CYL/CART periodicity (should be zero): " << DIST2(Rt,t)/DOT_PRODUCT(t,t) << endl;
          }

        }

	break;

      case UGP_IO_BF_ZONE_HEADER:

	assert(header.idata[0] == FA_ZONE_BOUNDARY);
	assert(header.idata[1] == bfZoneVec.size());

	// ensure the name is unique!... TODO do we need this (below it is not used for data lookup)
	{
	  string name = header.name;

          //if (nameSet.find(name) != nameSet.end()) {
          //  CERR("non-unique zone name:" << name);
          //}
	  while (nameSet.find(name) != nameSet.end())
	    name += "1";
	  nameSet.insert(name);
	  bfZoneVec.push_back(BfZone());
          bfZoneVec.back().setName(name);
          // add to map for easy lookup...
          bfZoneNameMap[bfZoneVec.back().name] = bfZoneVec.size()-1;
	}

	bfZoneVec.back().nbf_global   = ByteSwap::getInt8FromLswMswPair(header.idata+2);
	bfZoneVec.back().noobf_global = ByteSwap::getInt8FromLswMswPair(header.idata+4);
	bfZoneVec.back().ibf_f_global = ByteSwap::getInt8FromLswMswPair(header.idata+6);
	bfZoneVec.back().nob_f_global = ByteSwap::getInt8FromLswMswPair(header.idata+8);
	bfZoneVec.back().sbobf_global = ByteSwap::getInt8FromLswMswPair(header.idata+10);
	bfZoneVec.back().stob_f_global = ByteSwap::getInt8FromLswMswPair(header.idata+12);
        // spobf uses the first 2 uint8's in the header...
	bfZoneVec.back().spobf_global = header.ui8data[0];
	bfZoneVec.back().spob_f_global = header.ui8data[1];
	// and geom...
	bfZoneVec.back().area_global            = header.rdata[0];
	bfZoneVec.back().area_over_delta_global = header.rdata[1];
	bfZoneVec.back().n_global[0]            = header.rdata[2];
	bfZoneVec.back().n_global[1]            = header.rdata[3];
	bfZoneVec.back().n_global[2]            = header.rdata[4];
	bfZoneVec.back().x_global[0]            = header.rdata[5];
	bfZoneVec.back().x_global[1]            = header.rdata[6];
	bfZoneVec.back().x_global[2]            = header.rdata[7];
	// check...
	nbf_global_check += bfZoneVec.back().nbf_global;
        // striping...
	assert(bfZoneVec.back().bfora == NULL);
        MiscUtils::calcThresholdDist(bfZoneVec.back().bfora,bfZoneVec.back().nbf_global,mpi_size,DIST_THRESHOLD); 
        assert(bfZoneVec.back().bfora[mpi_rank+1]-bfZoneVec.back().bfora[mpi_rank] < TWO_BILLION);
	bfZoneVec.back().nbf = bfZoneVec.back().bfora[mpi_rank+1]-bfZoneVec.back().bfora[mpi_rank];

	if (mpi_rank == 0) {
	  cout << " > bf zone \"" << bfZoneVec.back().name << "\", nbf_global: " << bfZoneVec.back().nbf_global << 
	    ", area_global: " << bfZoneVec.back().area_global <<
	    ", n_global: " << COUT_VEC(bfZoneVec.back().n_global) <<
	    ", x_global: " << COUT_VEC(bfZoneVec.back().x_global) << endl;
	}
	
	break;
	
      case UGP_IO_CVOBF_INT8:

	if (mpi_rank == 0) 
	  cout << " > cvobf_global..." << endl;
	
	assert(nbf_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));
	assert(cvobf_global == NULL);
	assert(bfora != NULL);
	
	cvobf_global = new int8[nbf];
        readChunkedData<int8>(fh,offset+header_size+bfora[mpi_rank]*int8_size,cvobf_global,nbf,byte_swap,mpi_comm);
	
	break;
	
      case UGP_IO_NOOBF_I_AND_V_INT8:

	if (mpi_rank == 0) 
	  cout << " > noobf_i and noobf_v_global..." << endl;
	
	assert(nbf_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));
	assert(noobf_i == NULL);
	assert(noobf_v_global == NULL);
	assert(bfora != NULL);
	
	// recall noobf_i stores the count, so it is a simple int...

	noobf_i = new int[nbf+1];
        readChunkedData<int>(fh,offset+header_size+bfora[mpi_rank]*int_size,noobf_i+1,nbf,byte_swap,mpi_comm);
	
	{
	  
	  // it is possible that noobf_i could exceed the int limit when turned into a CSR -- check this...
	  // if this happens, run on more nodes...

	  int8 noobf_s = 0;
	  for (int ibf = 0; ibf < nbf; ++ibf)
	    noobf_s += noobf_i[ibf+1];
	  assert(noobf_s < TWO_BILLION);
	  
	  noobf_i[0] = 0;
	  for (int ibf = 0; ibf < nbf; ++ibf) 
	    noobf_i[ibf+1] += noobf_i[ibf];
	  assert(noobf_s == noobf_i[nbf]);
	  
	  int8 noobf_offset;
	  MPI_Scan(&noobf_s,&noobf_offset,1,MPI_INT8,MPI_SUM,mpi_comm);
	  // the offset on the last rank will be the total -- bcast back to everyone...
	  
	  int8 noobf_s_global = noobf_offset;
	  MPI_Bcast(&noobf_s_global,1,MPI_INT8,mpi_size-1,mpi_comm);
	  noobf_offset -= noobf_s;
	  assert(noobf_s_global == ByteSwap::getInt8FromLswMswPair(header.idata+2));
	  
	  noobf_v_global = new int8[noobf_s];
          readChunkedData<int8>(fh,offset+header_size+bfora[mpi_size]*int_size+noobf_offset*int8_size,noobf_v_global,noobf_s,byte_swap,mpi_comm);
	  
	}
	
	break;
	
      case UGP_IO_SPOBF_I_V_WGT:
	
	if (mpi_rank == 0) {
	  cout << " > spobf_i/v/wgt..." << endl;
	}
	
	assert(nbf_global == header.ui8data[0]);
	assert(surface_spobf_i == NULL);
	assert(surface_spobf_v_global == NULL);
	assert(surface_spobf_wgt_global == NULL);
	
	surface_spobf_i = new int[nbf+1];
        readChunkedData<int>(fh,offset+header_size+bfora[mpi_rank]*int_size,surface_spobf_i+1,nbf,byte_swap,mpi_comm);
	
	// recall spobf_i stores the count, so it is a simple int...
        
	{
	  
	  // it is possible that spobf_i could exceed the int limit when turned into a CSR -- check this...
	  // if this happens, run on more nodes...

	  int8 surface_spobf_s = 0;
	  for (int ibf = 0; ibf < nbf; ++ibf)
	    surface_spobf_s += surface_spobf_i[ibf+1];
	  assert(surface_spobf_s < TWO_BILLION);
	  
	  surface_spobf_i[0] = 0;
	  for (int ibf = 0; ibf < nbf; ++ibf) 
	    surface_spobf_i[ibf+1] += surface_spobf_i[ibf];
	  assert(surface_spobf_s == surface_spobf_i[nbf]);
	  
	  int8 surface_spobf_offset;
	  MPI_Scan(&surface_spobf_s,&surface_spobf_offset,1,MPI_INT8,MPI_SUM,mpi_comm);
	  // the offset on the last rank will be the total -- bcast back to everyone...
	  
	  int8 surface_spobf_s_global = surface_spobf_offset;
	  MPI_Bcast(&surface_spobf_s_global,1,MPI_INT8,mpi_size-1,mpi_comm);
	  surface_spobf_offset -= surface_spobf_s;
	  assert(surface_spobf_s_global == header.ui8data[1]);
	  
	  surface_spobf_v_global = new uint8[surface_spobf_s];
          readChunkedData<uint8>(fh,offset+header_size+bfora[mpi_size]*int_size+surface_spobf_offset*uint8_size,
                                 surface_spobf_v_global,surface_spobf_s,byte_swap,mpi_comm);
	  
	  surface_spobf_wgt_global = new double[surface_spobf_s];
          readChunkedData<double>(fh,offset+header_size+bfora[mpi_size]*int_size+surface_spobf_s_global*uint8_size+surface_spobf_offset*double_size,
                                  surface_spobf_wgt_global,surface_spobf_s,byte_swap,mpi_comm);
	
        }
	
	break;

      case UGP_IO_SBOBF_I_AND_V:
	
	if (mpi_rank == 0) 
	  cout << " > sbobf_i and sbobf_v_global..." << endl;
	
	assert(nbf_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));
	assert(surface_sbobf_i == NULL);
	assert(surface_sbobf_v_global == NULL);
	
	surface_sbobf_i = new int[nbf+1];
        readChunkedData<int>(fh,offset+header_size+bfora[mpi_rank]*int_size,surface_sbobf_i+1,nbf,byte_swap,mpi_comm);
	
	// recall sbobf_i stores the count, so it is a simple int...
        
	{
	  
	  // it is possible that sbobf_i could exceed the int limit when turned into a CSR -- check this...
	  // if this happens, run on more nodes...

	  int8 surface_sbobf_s = 0;
	  for (int ibf = 0; ibf < nbf; ++ibf)
	    surface_sbobf_s += surface_sbobf_i[ibf+1];
	  assert(surface_sbobf_s < TWO_BILLION);
	  
	  surface_sbobf_i[0] = 0;
	  for (int ibf = 0; ibf < nbf; ++ibf) 
	    surface_sbobf_i[ibf+1] += surface_sbobf_i[ibf];
	  assert(surface_sbobf_s == surface_sbobf_i[nbf]);
	  
	  int8 surface_sbobf_offset;
	  MPI_Scan(&surface_sbobf_s,&surface_sbobf_offset,1,MPI_INT8,MPI_SUM,mpi_comm);
	  // the offset on the last rank will be the total -- bcast back to everyone...
	  
	  int8 surface_sbobf_s_global = surface_sbobf_offset;
	  MPI_Bcast(&surface_sbobf_s_global,1,MPI_INT8,mpi_size-1,mpi_comm);
	  surface_sbobf_offset -= surface_sbobf_s;
	  assert(surface_sbobf_s_global == ByteSwap::getInt8FromLswMswPair(header.idata+2));
	  
	  surface_sbobf_v_global = new uint8[surface_sbobf_s];
          readChunkedData<uint8>(fh,offset+header_size+bfora[mpi_size]*int_size+surface_sbobf_offset*uint8_size,surface_sbobf_v_global,surface_sbobf_s,byte_swap,mpi_comm);
	  
	}
	
	break;

      case UGP_IO_CVOFA_INT8:

	if (mpi_rank == 0) 
	  cout << " > cvofa_global..." << endl;
	
	assert(nfa_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));
	assert(header.idata[2] == 2);
	assert(cvofa_global == NULL);
	assert(faora != NULL);
	
	cvofa_global = new int8[nfa][2];
        readChunkedData<int8>(fh,offset+header_size+faora[mpi_rank]*int8_size*2,(int8*)cvofa_global,int8(nfa)*2,byte_swap,mpi_comm);
	
	break;

      case UGP_IO_NOOFA_I_AND_V_INT8:
	
	if (mpi_rank == 0) 
	  cout << " > noofa_i and noofa_v_global..." << endl;
	
	assert(nfa_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));
	assert(noofa_i == NULL);
	assert(noofa_v_global == NULL);
	assert(faora != NULL);
	
	// recall noofa_i stores the count, so it is a simple int...

	noofa_i = new int[nfa+1];
        readChunkedData<int>(fh,offset+header_size+faora[mpi_rank]*int_size,noofa_i+1,nfa,byte_swap,mpi_comm);
	
	{
	  
	  // it is possible that noofa_i could exceed the int limit when turned into a CSR -- check this...
	  // if this happens, run on more nodes...

	  int8 noofa_s = 0;
	  for (int ifa = 0; ifa < nfa; ++ifa)
	    noofa_s += noofa_i[ifa+1];
	  assert(noofa_s < TWO_BILLION);
	  
	  noofa_i[0] = 0;
	  for (int ifa = 0; ifa < nfa; ++ifa) 
	    noofa_i[ifa+1] += noofa_i[ifa];
	  assert(noofa_s == noofa_i[nfa]);
	  
	  int8 noofa_offset;
	  MPI_Scan(&noofa_s,&noofa_offset,1,MPI_INT8,MPI_SUM,mpi_comm);
	  // the offset on the last rank will be the total -- bcast back to everyone...
	  
	  int8 noofa_s_global = noofa_offset;
	  MPI_Bcast(&noofa_s_global,1,MPI_INT8,mpi_size-1,mpi_comm);
	  noofa_offset -= noofa_s;
	  assert(noofa_s_global == ByteSwap::getInt8FromLswMswPair(header.idata+2));
	  
	  noofa_v_global = new int8[noofa_s];
          readChunkedData<int8>(fh,offset+header_size+faora[mpi_size]*int_size+noofa_offset*int8_size,noofa_v_global,noofa_s,byte_swap,mpi_comm);


	}
	
	break;

      case UGP_IO_PBI_PNO:
        {
          if (mpi_rank == 0) cout << " > pbi_pno..." << endl;
          
          assert(nno_p_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));
          assert(pbi_pno == NULL);
          assert(noora != NULL);
          
          // read onto nodal striping, figure out which nodes of nno are periodic...
          nno_p = 0;
          if (nno_p_global < noora[mpi_rank]) {
            nno_p = 0;
          }
          else if (nno_p_global >= noora[mpi_rank+1]) {
            nno_p = nno;
          }
          else {
            // TODO bissect...
            FOR_INO {
              if ((ino+noora[mpi_rank]) < nno_p_global) 
                ++nno_p;
            }
          }
          int8 my_nno_p = nno_p;
          int8 my_nno_p_disp;
          MPI_Scan(&my_nno_p,&my_nno_p_disp,1,MPI_INT8,MPI_SUM,mpi_comm);
          my_nno_p_disp -= my_nno_p;
          if (mpi_rank == mpi_size-1) assert(my_nno_p+my_nno_p_disp == nno_p_global);

          pbi_pno = new uint8[nno_p];
          readChunkedData<uint8>(fh,offset+header_size+my_nno_p_disp*int8_size,pbi_pno,nno_p,byte_swap,mpi_comm);
        }
	break;

      case UGP_IO_X_NO:

	if (mpi_rank == 0) cout << " > x_no..." << endl;
	
        if (io_version == 3) {
          assert(nno_global == header.idata[0]);
          assert(header.idata[1] == 3);
        }
        else { 
          assert(nno_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));
	  assert(header.idata[2] == 3);
        }
	assert(x_no == NULL);
	assert(noora != NULL);
	
	x_no = new double[nno][3];
        readChunkedData<double>(fh,offset+header_size+noora[mpi_rank]*double_size*3,(double*)x_no,nno*3,byte_swap,mpi_comm);
	
	break;
	
      case UGP_IO_CV_ZONE_HEADER:

	assert(header.idata[0] == CV_ZONE_FLUID);
	assert(header.idata[1] == cvZoneVec.size());
	cvZoneVec.push_back(CvZone(header.name));
	cvZoneVec.back().ncv_global   = ByteSwap::getInt8FromLswMswPair(header.idata+2);
	cvZoneVec.back().icv_f_global = ByteSwap::getInt8FromLswMswPair(header.idata+4);
	
	ncv_global_check += cvZoneVec.back().ncv_global;
	
	if (mpi_rank == 0) {
	  cout << " > cv zone \"" << cvZoneVec.back().name << "\", ncv_global: " << cvZoneVec.back().ncv_global << endl;
	}
	
	break;
	
      case UGP_IO_CV_D1:

        if (io_version == 3) 
          assert(ncv_global == header.idata[0]);
        else 
          assert(ncv_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));
	assert(cvora != NULL);
	
        if ((strcmp(header.name,"vol_cv") == 0)||(strcmp(header.name,"VOL_VV") == 0)||(strcmp(header.name,"vol_vv") == 0)) {
	  
	  if (mpi_rank == 0) cout << " > vol_cv..." << endl;

	  assert(vol_cv == NULL);
	  vol_cv = new double[ncv];
          readChunkedData<double>(fh,offset+header_size+cvora[mpi_rank]*double_size,vol_cv,ncv,byte_swap,mpi_comm);
	  
	}
        else if (strcmp(header.name,"r_vv") == 0) {
	  
	  if (mpi_rank == 0) cout << " > r_vv..." << endl;

	  assert(r_vv == NULL);
	  r_vv = new double[ncv];
          readChunkedData<double>(fh,offset+header_size+cvora[mpi_rank]*double_size,r_vv,ncv,byte_swap,mpi_comm);
	  
	}
	else {

	  if (mpi_rank == 0) cout << " > unrecognized CV_D1: " << header.name << endl;
	  
	}

	break;
	
      case UGP_IO_CV_D2:

        if (io_version == 3) {
          assert(ncv_global == header.idata[0]);
          assert(header.idata[1] == 3);
        }
        else {
          assert(ncv_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));
          assert(header.idata[2] == 3);
        }
	assert(cvora != NULL);
	
	if ((strcmp(header.name,"x_vv") == 0)||(strcmp(header.name,"X_VV") == 0)) {
	  
	  if (mpi_rank == 0) cout << " > x_vv..." << endl;
	  
	  assert(x_vv == NULL);
	  x_vv = new double[ncv][3];
          readChunkedData<double>(fh,offset+header_size+cvora[mpi_rank]*double_size*3,(double*)x_vv,ncv*3,byte_swap,mpi_comm);
	  
	}
        else if (strcmp(header.name,"x_cv") == 0) {
	  
	  if (mpi_rank == 0) cout << " > x_cv..." << endl;
	  
	  assert(x_cv == NULL);
	  x_cv = new double[ncv][3];
          readChunkedData<double>(fh,offset+header_size+cvora[mpi_rank]*double_size*3,(double*)x_cv,ncv*3,byte_swap,mpi_comm);
	  
	}
	else {

	  if (mpi_rank == 0) cout << " > unrecognized CV_D2: " << header.name << endl;
	  
	}

	break;
	
      case UGP_IO_BF_D1:
	
	assert(nbf_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));
	assert(bfora != NULL);
	
	if (strcmp(header.name,"area_bf") == 0) {
	  
	  if (mpi_rank == 0) cout << " > area_bf..." << endl;
	  
	  assert(area_bf == NULL);
	  area_bf = new double[nbf];
          readChunkedData<double>(fh,offset+header_size+bfora[mpi_rank]*double_size,area_bf,nbf,byte_swap,mpi_comm);
	  
	}
	else if (strcmp(header.name,"area_over_delta_bf") == 0) {
	  
	  if (mpi_rank == 0) cout << " > area_over_delta_bf..." << endl;
	  
	  assert(area_over_delta_bf == NULL);
	  area_over_delta_bf = new double[nbf];
          readChunkedData<double>(fh,offset+header_size+bfora[mpi_rank]*double_size,area_over_delta_bf,nbf,byte_swap,mpi_comm);
	  
	}
	else {

	  if (mpi_rank == 0) cout << " > unrecognized BF_D1: " << header.name << endl;
	  
	}
	
	break;
	
      case UGP_IO_BF_D2:
	
	assert(nbf_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));
	assert(header.idata[2] == 3);
	assert(bfora != NULL);
	
	if (strcmp(header.name,"n_bf") == 0) {
	  
	  if (mpi_rank == 0) cout << " > n_bf..." << endl;
	  
	  assert(n_bf == NULL);
	  n_bf = new double[nbf][3];
          readChunkedData<double>(fh,offset+header_size+bfora[mpi_rank]*double_size*3,(double*)n_bf,nbf*3,byte_swap,mpi_comm);
	  
	}
	else if (strcmp(header.name,"x_bf") == 0) {
	  
	  if (mpi_rank == 0) cout << " > x_bf..." << endl;
	  
	  assert(x_bf == NULL);
	  x_bf = new double[nbf][3];
          readChunkedData<double>(fh,offset+header_size+bfora[mpi_rank]*double_size*3,(double*)x_bf,nbf*3,byte_swap,mpi_comm);
	  
	}
	else {

	  if (mpi_rank == 0) cout << " > unrecognized BF_D2: " << header.name << endl;
	  
	}
	
	break;

      case UGP_IO_BF_DN33:
	
	assert(nbf_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));
	assert(header.idata[2] == 3);
	assert(header.idata[3] == 3);
	assert(bfora != NULL);
	
	if (strcmp(header.name,"Gij_bf") == 0) {
	  
	  if (mpi_rank == 0) cout << " > Gij_bf..." << endl;
	  
	  assert(Gij_bf == NULL);
	  Gij_bf = new double[nbf][3][3];
          readChunkedData<double>(fh,offset+header_size+bfora[mpi_rank]*double_size*9,(double*)Gij_bf,nbf*9,byte_swap,mpi_comm);
	  
	}
	else {

	  if (mpi_rank == 0) cout << " > unrecognized BF_DN33: " << header.name << endl;
	  
	}
	
	break;

      case UGP_IO_FA_D2:
	
	assert(nfa_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));
	assert(header.idata[2] == 3);
	assert(faora != NULL);
	
	if (strcmp(header.name,"n_fa") == 0) {
	  
	  if (mpi_rank == 0) cout << " > n_fa..." << endl;
	  
	  assert(n_fa == NULL);
	  n_fa = new double[nfa][3];
          readChunkedData<double>(fh,offset+header_size+faora[mpi_rank]*double_size*3,(double*)n_fa,nfa*3,byte_swap,mpi_comm);
	  
	}
	else if (strcmp(header.name,"x_fa") == 0) {
	  
	  if (mpi_rank == 0) cout << " > x_fa..." << endl;
	  
	  assert(x_fa == NULL);
	  x_fa = new double[nfa][3];
          readChunkedData<double>(fh,offset+header_size+faora[mpi_rank]*double_size*3,(double*)x_fa,nfa*3,byte_swap,mpi_comm);
	  
	}
	else {

	  if (mpi_rank == 0) cout << " > unrecognized FA_D2: " << header.name << endl;
	  
	}
	
	break;

        /*
          case UGP_IO_VV_BBOX:

          if (strcmp(header.name,"vv_bbox") == 0) {
	  
	  if (mpi_rank == 0) cout << " > ncvobb and vv_bbox..." << endl;
        
          assert(header.idata[0] >= 1); // should have a bounding box
          assert(bbora == NULL);
          assert(ncvobb == NULL);
          assert(vv_bbox == NULL);

          nbb_global = header.idata[0];
          MiscUtils::calcUniformDist(bbora,nbb_global,mpi_size); // should this be thresholded?
          assert(bbora[mpi_rank+1]-bbora[mpi_rank] < TWO_BILLION); // necessary?
          nbb = bbora[mpi_rank+1]-bbora[mpi_rank];

          ncvobb = new int[nbb];
          readChunkedData<int>(fh,offset+header_size+bbora[mpi_rank]*int_size,(int*)ncvobb,nbb,byte_swap,mpi_comm);

          vv_bbox = new double[nbb][6]; // 6 coords of bbox
          readChunkedData<double>(fh,offset+nbb_global*int_size+header_size+bbora[mpi_rank]*double_size*6,(double*)vv_bbox,nbb*6,byte_swap,mpi_comm);

          }

          else {

	  if (mpi_rank == 0) cout << " > unrecognized VV_BBOX: " << header.name << endl;
	  
          }

          break;
        */

      case UGP_IO_SURFACE:

	assert(b_surface);

	if (mpi_rank == 0) cout << " > surface (xp,spost,znost)..." << endl;

	assert(header.idata[0] > 0); // should have at least one tri? what about triple-periodic box?
	assert(header.idata[1] > 0);
        assert(surface_spora == NULL);
        assert(surface_stora == NULL);
        assert(surface_xp == NULL);
        assert(surface_spost == NULL);
        assert(surface_znost == NULL);

	surface_nsp_global = header.idata[0]; 
	surface_nst_global = header.idata[1];

        MiscUtils::calcThresholdDist(surface_spora,surface_nsp_global,mpi_size,DIST_THRESHOLD);
        surface_nsp = surface_spora[mpi_rank+1]-surface_spora[mpi_rank];
	
        MiscUtils::calcThresholdDist(surface_stora,surface_nst_global,mpi_size,DIST_THRESHOLD);
        surface_nst = surface_stora[mpi_rank+1]-surface_stora[mpi_rank];
	
        surface_xp = new double[surface_nsp][3]; 
        readChunkedData<double>(fh,offset+header_size+surface_spora[mpi_rank]*double_size*3,(double*)surface_xp,surface_nsp*3,byte_swap,mpi_comm);

        surface_spost = new int[surface_nst][3]; 
        readChunkedData<int>(fh,offset+header_size+surface_nsp_global*double_size*3+surface_stora[mpi_rank]*int_size*3,(int*)surface_spost,surface_nst*3,byte_swap,mpi_comm);

        surface_znost = new int[surface_nst]; 
        readChunkedData<int>(fh,offset+header_size+surface_nsp_global*double_size*3+surface_nst_global*int_size*3+surface_stora[mpi_rank]*int_size,surface_znost,surface_nst,byte_swap,mpi_comm);

	break;

        /*
          case UGP_IO_SURFACE_PERIODIC_INFO:
          {

          assert(b_surface);

          if (mpi_rank == 0) cout << " > surface periodic info..." << endl;

          const int npt  = header.idata[0]; 
          const int npz  = header.idata[1];
          surface_npp_global = header.idata[2];
          if (mpi_rank == 0) cout << " > " << npt << " periodic transforms, " << npz << " periodic zones and " << surface_npp_global << " periodic points. " << endl;

          MiscUtils::calcThresholdDist(surface_ppora,surface_npp_global,mpi_size,DIST_THRESHOLD);
          surface_npp = surface_ppora[mpi_rank+1]-surface_ppora[mpi_rank];

          assert(surface_isp_pp == NULL); surface_isp_pp = new int[surface_npp];
          assert(surface_isp_p_pp == NULL); surface_isp_p_pp = new int[surface_npp];
          assert(surface_isp_bits_pp == NULL); surface_isp_bits_pp = new uint2[surface_npp];
          const int skip = header_size+npz*(int_size+sizeof(uint2)); // skipped some periodic zone info (needed?)
          readChunkedData<int>(fh,offset+skip+surface_ppora[mpi_rank]*int_size,surface_isp_pp,surface_npp,byte_swap,mpi_comm);
          readChunkedData<int>(fh,offset+skip+surface_ppora[mpi_rank]*int_size+surface_npp_global*int_size,surface_isp_p_pp,surface_npp,byte_swap,mpi_comm);
          readChunkedData<uint2>(fh,offset+skip+surface_ppora[mpi_rank]*sizeof(uint2)+2*surface_npp_global*int_size,surface_isp_bits_pp,surface_npp,byte_swap,mpi_comm);

          }
          break;
        */

      case UGP_IO_HASHIDS:

        b_foundHash = true;
        //read two hash id's from mles file. First id identifies
        //the current file, store in mlesHash.  Ignore the second
        //field for now.  TODO include surface hash as parent
        RestartHashUtilities::mlesReadHashes(fh, offset, header); //will bCast to all ranks
        if (mpi_rank==0){
          cout << " > mles hash: " << RestartHashUtilities::mlesHash << endl;
        }
        break;

      case UGP_IO_EOF:
        
	done = 1;
        break;
	
      default:
	
	if (mpi_rank == 0) cout << " > unknown header: " << header.id << " \"" << header.name << "\", skipping." << endl;
	
      }
      
      offset += header.skip;
    
    }

    assert(nbf_global_check == nbf_global);
    assert(ncv_global_check == ncv_global);
    
    MPI_File_close(&fh);
    
    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0) {
      const double seconds = MPI_Wtime() - wtime0;
      cout << " > read_mles summary: ncv_global: " << ncv_global << 
	" file size: " << double(offset)/1.0E+9 << " [GB] read rate: " << double(offset)/1.0E+9/seconds << " [GB/s]" << endl;
    }
    
    if (!b_foundHash){ //set a default mles hash if one was not found
      std::stringstream ss;
      ss << "Missing mles Hash";
      RestartHashUtilities::mlesHash.init(ss,RestartHashUtilities::sha1hashlength);
      if (mpi_rank==0) cout << " > setting default mles hash id " <<  RestartHashUtilities::mlesHash << endl;
    }

    if (io_version != 3) {

      // we can now set the zone_bf...
      assert(zone_bf == NULL);
      int izone = 0;
      nbf_global_check = 0;
      zone_bf = new int[nbf];
      for (int ibf = 0; ibf < nbf; ++ibf) {
        const int8 ibf_global = bfora[mpi_rank] + ibf;
        assert(izone < bfZoneVec.size());
        while (nbf_global_check + bfZoneVec[izone].nbf_global <= ibf_global) {
          nbf_global_check += bfZoneVec[izone].nbf_global;
          ++izone;
          assert(izone < bfZoneVec.size());
        }
        zone_bf[ibf] = izone;
      }

    }

    // if this is io_version 3 we need to split the read faces into bf and fa data. We also need to
    // populate some missing data...
    if (io_version == 3) {
      
      if (mpi_rank == 0) 
        cout << " > reconstructing mles data from read in v3 les data." << endl;
       
      fabfora = faora;
      int8 (*cvofabf_global)[2] = cvofa_global;
      zone_fabf = zone_bf;
      int* noofabf_i = noofa_i;
      int8* noofabf_v_global = noofa_v_global;

      nfabf_global = nfa_global;
      nfa_global = nfa_zone[26];
      faora = NULL;
      MiscUtils::calcUniformDist(faora,nfa_global,mpi_size);
      nfabf = nfa;
      nfa = faora[mpi_rank+1]-faora[mpi_rank];
      int nbf_tmp = 0;
      int nfa_tmp = 0;
      const int nzone = bfZoneVec.size();
      for (int ifabf = 0; ifabf < nfabf; ++ifabf) {
        if (zone_fabf[ifabf] < nzone) {
          assert(zone_fabf[ifabf] >= 0);
          ++nbf_tmp; 
        }
        else {
          assert(zone_fabf[ifabf] == nzone);
          ++nfa_tmp;
        }
      }

      // ifa...

      int8 my_nfa = (int8)nfa_tmp;
      assert(my_nfa_disp == -1);
      MPI_Scan(&my_nfa,&my_nfa_disp,1,MPI_INT8,MPI_SUM,mpi_comm);
      if (mpi_rank == mpi_size-1)
        assert(my_nfa_disp == nfa_global);
      my_nfa_disp -= my_nfa;;
      int *send_count = new int[mpi_size];
      int *send_disp  = new int[mpi_size];
      FOR_RANK send_count[rank] = 0;
      int8 * send_buf_int8;
      for (int iter = 0; iter < 2; ++iter) {
        int ii = 0;
        for (int ifabf = 0; ifabf < nfabf; ++ifabf) {
          if (zone_fabf[ifabf] == nzone) {
            const int8 ifa_global = ii+my_nfa_disp; ++ii;
            const int rank = MiscUtils::getRankInXora(ifa_global,faora);
            assert((ifa_global >= 0)&&(ifa_global < nfa_global));
            const int nnof = noofabf_i[ifabf+1]-noofabf_i[ifabf];
            if (iter == 0) {
              send_count[rank] += 4 + nnof; // ifa, cvofa[2], nnof and nodes
            }
            else {
              send_buf_int8[send_disp[rank]++] = ifa_global-faora[rank];
              FOR_I2 send_buf_int8[send_disp[rank]++] = cvofabf_global[ifabf][i];
              send_buf_int8[send_disp[rank]++] = nnof;
              for (int nof = noofabf_i[ifabf]; nof != noofabf_i[ifabf+1]; ++nof)
                send_buf_int8[send_disp[rank]++] = noofabf_v_global[nof];
            }
          }
        }

        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

        if (iter == 0) {
          send_buf_int8 = new int8[send_disp[mpi_size-1] + send_count[mpi_size-1]];
        }
      }

      cvofa_global = new int8[nfa][2];
      noofa_i = new int[nfa+1]; noofa_i[0] = 0;
	
      //  exchange...
      
      int *recv_count = new int[mpi_size];
      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
      
      int *recv_disp = new int[mpi_size];
      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
      
      int8* recv_buf_int8 = new int8[recv_count_sum];
      MPI_Alltoallv(send_buf_int8,send_count,send_disp,MPI_INT8,
                    recv_buf_int8,recv_count,recv_disp,MPI_INT8,mpi_comm);
      delete[] send_buf_int8; send_buf_int8 = NULL; 

      for (int iter = 0; iter < 2; ++iter) {
        int irecv = 0;
        while (irecv < recv_count_sum) {
          const int ifa = recv_buf_int8[irecv++];
          if (iter == 0) {
            FOR_I2 cvofa_global[ifa][i] = recv_buf_int8[irecv++];
            noofa_i[ifa+1] = recv_buf_int8[irecv++];
            irecv += noofa_i[ifa+1];
          }
          else {
            irecv += 3; 
            for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) 
              noofa_v_global[nof] = recv_buf_int8[irecv++];
          }
        }

        if (iter == 0) {
          for (int ifa = 0; ifa < nfa; ++ifa) 
            noofa_i[ifa+1] += noofa_i[ifa];
          noofa_v_global = new int8[noofa_i[nfa]];
        }
      }
      delete[] recv_buf_int8; recv_buf_int8 = NULL;

      // ibf...
      
      int8* my_nbf_zone = new int8[nzone];
      assert(zone_ifabf_vec.empty()); zone_ifabf_vec.resize(nbf_tmp);
      for (int izone = 0; izone < nzone; ++izone) my_nbf_zone[izone] = 0;
      {
        int ii = 0;
        for (int ifabf = 0; ifabf < nfabf; ++ifabf) {
          if (zone_fabf[ifabf] < nzone) {
            zone_ifabf_vec[ii].first = zone_fabf[ifabf];
            zone_ifabf_vec[ii].second = ifabf;
            ++ii;
            ++my_nbf_zone[zone_fabf[ifabf]];
          }
        }
        sort(zone_ifabf_vec.begin(),zone_ifabf_vec.end());
      }
      assert(my_nbf_disp == NULL); my_nbf_disp = new int8[nzone];
      MPI_Scan(my_nbf_zone,my_nbf_disp,nzone,MPI_INT8,MPI_SUM,mpi_comm);
      for (int izone = 0; izone < nzone; ++izone)
        my_nbf_disp[izone] -= my_nbf_zone[izone];
      delete[] my_nbf_zone;

      nbf_global_check = 0;
      FOR_IZONE(bfZoneVec) {
	bfZoneVec[izone].ibf_f_global = nbf_global_check;
        nbf_global_check += bfZoneVec[izone].nbf_global;
      }
      assert(nbf_global_check == nbf_global);

      FOR_RANK send_count[rank] = 0;
      assert(send_buf_int8 == NULL);
      for (int iter = 0; iter < 2; ++iter) {
        int disp = 0;
        int izone0 = -1;
        for (int ii = 0; ii < nbf_tmp; ++ii) {
          const int izone = zone_ifabf_vec[ii].first;
          if (izone != izone0) {
            disp = ii;
            izone0 = izone;
          }
          const int ifabf = zone_ifabf_vec[ii].second;
          assert(zone_fabf[ifabf] == izone);
          const int8 ibf_global = ii-disp+(my_nbf_disp[izone]+bfZoneVec[izone].ibf_f_global);
          assert((ibf_global >= 0)&&(ibf_global < nbf_global));
          const int rank = MiscUtils::getRankInXora(ibf_global,bfora);
          const int nnob = noofabf_i[ifabf+1]-noofabf_i[ifabf];
          if (iter == 0) {
            send_count[rank] += 4 + nnob; // ifa, cvobf, zone_bf, nnob and nodes
          }
          else {
            const int ibf = ibf_global-bfora[rank]; 
            send_buf_int8[send_disp[rank]++] = ibf;
            send_buf_int8[send_disp[rank]++] = cvofabf_global[ifabf][0];
            send_buf_int8[send_disp[rank]++] = izone;
            send_buf_int8[send_disp[rank]++] = nnob;
            for (int nob = noofabf_i[ifabf]; nob != noofabf_i[ifabf+1]; ++nob)
              send_buf_int8[send_disp[rank]++] = noofabf_v_global[nob];
          }
        }

        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

        if (iter == 0) {
          send_buf_int8 = new int8[send_disp[mpi_size-1] + send_count[mpi_size-1]];
        }
      }

      cvobf_global = new int8[nbf];
      zone_bf = new int[nbf];
      noobf_i = new int[nbf+1]; noobf_i[0] = 0;
	
      //  exchange...
      
      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
      
      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
      
      assert(recv_buf_int8 == NULL); recv_buf_int8 = new int8[recv_count_sum];
      MPI_Alltoallv(send_buf_int8,send_count,send_disp,MPI_INT8,
                    recv_buf_int8,recv_count,recv_disp,MPI_INT8,mpi_comm);
      delete[] send_buf_int8; send_buf_int8 = NULL; 

      for (int iter = 0; iter < 2; ++iter) {
        int irecv = 0;
        while (irecv < recv_count_sum) {
          const int ibf = recv_buf_int8[irecv++];
          if (iter == 0) {
            cvobf_global[ibf] = recv_buf_int8[irecv++];
            zone_bf[ibf] = recv_buf_int8[irecv++];
            noobf_i[ibf+1] = recv_buf_int8[irecv++];
            irecv += noobf_i[ibf+1];
          }
          else {
            irecv += 3; 
            for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) 
              noobf_v_global[nob] = recv_buf_int8[irecv++];
          }
        }

        if (iter == 0) {
          for (int ibf = 0; ibf < nbf; ++ibf) 
            noobf_i[ibf+1] += noobf_i[ibf];
          noobf_v_global = new int8[noobf_i[nbf]];
        }
      }
      delete[] recv_buf_int8; recv_buf_int8 = NULL;

      // cleanup...
      
      delete[] cvofabf_global;
      delete[] noofabf_i;
      delete[] noofabf_v_global;
      
      // bfZoneVec...

      int8* my_noobf_zone = new int8[nzone];
      for (int izone = 0; izone < nzone; ++izone)
        my_noobf_zone[izone] = 0;
      for (int ibf = 0; ibf < nbf; ++ibf) {
        const int izone = zone_bf[ibf];
        my_noobf_zone[izone] += noobf_i[ibf+1]-noobf_i[ibf];
      }
      int8 * my_noobf_disp = new int8[nzone];
      MPI_Scan(my_noobf_zone,my_noobf_disp,nzone,MPI_INT8,MPI_SUM,mpi_comm);
      int8 * noobf_zone = new int8[nzone];
      if (mpi_rank == mpi_size-1)
        for (int izone = 0; izone < nzone; ++izone)
          noobf_zone[izone] = my_noobf_disp[izone];
      MPI_Bcast(noobf_zone,nzone,MPI_INT8,mpi_size-1,mpi_comm);
      delete[] my_noobf_zone;
      delete[] my_noobf_disp;

      //int nzone_final = 0;
      int8 noobf_end = 0;
      //int8 sbobf_end = 0;
      for (int izone = 0; izone < nzone; ++izone) {
        const int8 noobf_begin = noobf_end; noobf_end += noobf_zone[izone];
        bfZoneVec[izone].nob_f_global = noobf_begin;
        bfZoneVec[izone].noobf_global = noobf_zone[izone];
      }
      delete[] noobf_zone;

      // get nodes on bf striping...
      
      map<const int8,int> globalMap;
      int nno_tmp = 0;
      FOR_IBF {
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int8 ino_global = noobf_v_global[nob];
          assert((ino_global >= 0)&&(ino_global < nno_global));
          map<const int8,int>::iterator iter = globalMap.find(ino_global);
          if (iter == globalMap.end()) 
            globalMap[ino_global] = nno_tmp++;
        }
      }
      int8* ino_global_tmp = new int8[nno_tmp];
      int *noobf_v_tmp = new int[noobf_i[nbf]];
      FOR_IBF {
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int8 ino_global = noobf_v_global[nob];
          assert((ino_global >= 0)&&(ino_global < nno_global));
          map<const int8,int>::iterator iter = globalMap.find(ino_global);
          assert(iter != globalMap.end());
          noobf_v_tmp[nob] = iter->second;
          ino_global_tmp[iter->second] = ino_global;
        }
      }
      globalMap.clear();
      double (*x_no_tmp)[3] = new double[nno_tmp][3];
      {
        DistributedDataExchanger dde(ino_global_tmp,nno_tmp,noora);
        dde.pull(x_no_tmp,x_no);
      }
      delete[] ino_global_tmp; ino_global_tmp = NULL;

      double (*x_vv_tmp)[3] = new double[nbf][3];
      double *vol_cv_tmp = new double[nbf];
      {
        DistributedDataExchanger dde(cvobf_global,nbf,cvora);
        dde.pull(x_vv_tmp,x_vv);
        dde.pull(vol_cv_tmp,vol_cv);
      }
      
      assert(n_bf == NULL); n_bf = new double[nbf][3];
      assert(x_bf == NULL); x_bf = new double[nbf][3];
      assert(area_bf == NULL); area_bf = new double[nbf];
      assert(area_over_delta_bf == NULL); area_over_delta_bf = new double[nbf];
      assert(Gij_bf == NULL); Gij_bf = new double[nbf][3][3];

      // older versions of v3 don't include this record, so lets build it now...
      //bool build_rvv = false;
      double *r_vv_tmp = NULL;
      bool build_r_vv = false;
      if (r_vv == NULL) {
        build_r_vv = true;
        r_vv = new double[ncv]; 
        FOR_ICV r_vv[icv] = 0.0;
        r_vv_tmp = new double[nbf]; 
        FOR_IBF r_vv_tmp[ibf] = 0.0;
      }

      // these are calculated from vida vor's methodology...
      const double area_eps = 1.0E-20;
      int my_count = 0;
      FOR_IBF {
        FOR_I3 x_bf[ibf][i] = 0.0;
        FOR_I3 n_bf[ibf][i] = 0.0;

        // loop to get the normal and approximate x_bf...
        double x_bf_tmp[3] = { 0.0, 0.0, 0.0 };
        const int ino_f = noobf_v_tmp[noobf_i[ibf]];
        double wgt_sum = 0.0;
        int ino1 = noobf_v_tmp[noobf_i[ibf+1]-1];
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int ino0 = ino1;
          ino1 = noobf_v_tmp[nob];
          const double wgt = DIST(x_no_tmp[ino0],x_no_tmp[ino1]);
          FOR_I3 x_bf_tmp[i] += wgt*(x_no_tmp[ino0][i]+x_no_tmp[ino1][i]);
          wgt_sum += wgt;
          const double n[3] = TRI_NORMAL_2(x_no_tmp[ino_f],x_no_tmp[ino0],x_no_tmp[ino1]);
          FOR_I3 n_bf[ibf][i] += 0.5*n[i];
        }
        assert(wgt_sum > 0.0);
        FOR_I3 x_bf_tmp[i] /= wgt_sum*2.0;

        // and get the area of this face...
        area_bf[ibf] = MAG(n_bf[ibf]);

        // get exact x_bf...
        wgt_sum = 0.0;
        ino1 = noobf_v_tmp[noobf_i[ibf+1]-1];
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int ino0 = ino1;
          ino1 = noobf_v_tmp[nob];
          if (build_r_vv) r_vv_tmp[ibf] = max(r_vv_tmp[ibf],DIST(x_vv_tmp[ibf],x_no_tmp[ino1]));
          const double n[3] = TRI_NORMAL_2(x_bf_tmp,x_no_tmp[ino0],x_no_tmp[ino1]);
          const double wgt = DOT_PRODUCT(n_bf[ibf],n); // area squared weighting !?
          FOR_I3 x_bf[ibf][i] += wgt*(x_bf_tmp[i]+x_no_tmp[ino0][i]+x_no_tmp[ino1][i]);
          wgt_sum += wgt;
        }
        if (wgt_sum > 0.0) {
          FOR_I3 x_bf[ibf][i] /= wgt_sum*3.0;
        }
        else {
          // this must be a zero-area face -- so just pick a node...
          FOR_I3 x_bf[ibf][i] = x_no_tmp[ino1][i];
        }
        const double dx[3] = DIFF(x_bf[ibf],x_vv_tmp[ibf]);
        const double unit_n[3] = {n_bf[ibf][0]/max(area_bf[ibf],area_eps),n_bf[ibf][1]/max(area_bf[ibf],area_eps),n_bf[ibf][2]/max(area_bf[ibf],area_eps)};
        const double dp = DOT_PRODUCT(dx,unit_n);
        const double dx_vol = pow(vol_cv_tmp[ibf],1.0/3.0);
        const double dx_vol_area = vol_cv_tmp[ibf]/max(area_bf[ibf],area_eps);
        const double dx_clip = min( 0.25*dx_vol , 0.25*dx_vol_area );
        if (dp < dx_clip) {
          ++my_count;
          area_over_delta_bf[ibf] = area_bf[ibf]/dx_clip;
        }
        else {
          area_over_delta_bf[ibf] = area_bf[ibf]/dp;
        }
        FOR_I3 FOR_J3 Gij_bf[ibf][i][j] = -dx[j]*n_bf[ibf][i];

        //cout << zone_bf[ibf] << " " << COUT_VEC(n_bf[ibf]) << " " << COUT_VEC(x_bf[ibf]) << " " << area_bf[ibf] << " " << area_over_delta_bf[ibf] << endl;
      }
      delete[] noobf_v_tmp;
      delete[] x_no_tmp; x_no_tmp = NULL;
      delete[] x_vv_tmp; x_vv_tmp = NULL;
      delete[] vol_cv_tmp; vol_cv_tmp = NULL;
    
      int count;
      MPI_Reduce(&my_count,&count,1,MPI_INT,MPI_SUM,0,mpi_comm);
      if (mpi_rank == 0)
        cout << " > boundary length scale clipped in " << count << " faces." << endl;

      if (build_r_vv) {
        DistributedDataExchanger dde(cvobf_global,nbf,cvora);
        dde.push(r_vv,r_vv_tmp,MAX_DATA);
        delete[] r_vv_tmp; r_vv_tmp = NULL;
      }
      
      // build bfZoneVec geometry info...

      double (*my_geom_data_zone)[8] = new double[nzone][8];
      for (int izone = 0; izone < nzone; ++izone) {
        FOR_I8 my_geom_data_zone[izone][i] = 0.0;
      }
      FOR_IBF {
        const int izone = zone_bf[ibf];
        my_geom_data_zone[izone][0] += area_bf[ibf];
        my_geom_data_zone[izone][1] += area_over_delta_bf[ibf];
        FOR_I3 my_geom_data_zone[izone][2+i] += n_bf[ibf][i];
        FOR_I3 my_geom_data_zone[izone][5+i] += area_bf[ibf]*x_bf[ibf][i];
      }
      double (*geom_data_zone)[8] = new double[nzone][8];
      MPI_Allreduce((double*)my_geom_data_zone,(double*)geom_data_zone,8*nzone,MPI_DOUBLE,MPI_SUM,mpi_comm);
      delete[] my_geom_data_zone;
      for (int izone = 0; izone < nzone; ++izone) {
        bfZoneVec[izone].area_global = geom_data_zone[izone][0];
        bfZoneVec[izone].area_over_delta_global = geom_data_zone[izone][1];
        FOR_I3 bfZoneVec[izone].n_global[i] = geom_data_zone[izone][2+i];
        FOR_I3 bfZoneVec[izone].x_global[i] = geom_data_zone[izone][5+i]/geom_data_zone[izone][0];
        if (mpi_rank == 0) {
          cout << " > bf zone \"" << bfZoneVec[izone].name << "\"" << 
            ", area_global: " << bfZoneVec[izone].area_global <<
            ", area_over_delta_global: " << bfZoneVec[izone].area_over_delta_global <<
            ", n_global: " << COUT_VEC(bfZoneVec[izone].n_global) <<
            ", x_global: " << COUT_VEC(bfZoneVec[izone].x_global) << endl;
        }
      }
      delete[] geom_data_zone;

      // get nodes on fa striping...
      
      assert(globalMap.empty());
      nno_tmp = 0;
      FOR_IFA {
        for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
          const int8 ino_global = noofa_v_global[nof];
          assert((ino_global >= 0)&&(ino_global < nno_global));
          map<const int8,int>::iterator iter = globalMap.find(ino_global);
          if (iter == globalMap.end()) 
            globalMap[ino_global] = nno_tmp++;
        }
      }
      assert(ino_global_tmp == NULL); ino_global_tmp = new int8[nno_tmp];
      int *noofa_v_tmp = new int[noofa_i[nfa]];
      FOR_IFA {
        for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
          const int8 ino_global = noofa_v_global[nof];
          assert((ino_global >= 0)&&(ino_global < nno_global));
          map<const int8,int>::iterator iter = globalMap.find(ino_global);
          assert(iter != globalMap.end());
          noofa_v_tmp[nof] = iter->second;
          ino_global_tmp[iter->second] = ino_global;
        }
      }
      globalMap.clear();

      assert(x_no_tmp == NULL); x_no_tmp = new double[nno_tmp][3];
      {
        DistributedDataExchanger dde(ino_global_tmp,nno_tmp,noora);
        dde.pull(x_no_tmp,x_no);
      }
      delete[] ino_global_tmp;

      assert(x_vv_tmp == NULL); x_vv_tmp = new double[2*nfa][3];
      assert(vol_cv_tmp == NULL); vol_cv_tmp = new double[2*nfa];
      {
        DistributedDataExchanger dde((int8*)cvofa_global,2*nfa,cvora);
        dde.pull(x_vv_tmp,x_vv);
        dde.pull(vol_cv_tmp,vol_cv);
      }

      if (build_r_vv) {
        assert(r_vv_tmp == NULL); 
        r_vv_tmp = new double[2*nfa]; 
        FOR_IFA {
          r_vv_tmp[2*ifa  ] = 0.0;
          r_vv_tmp[2*ifa+1] = 0.0;
        }
      }
      
      assert(n_fa == NULL); n_fa = new double[nfa][3];
      assert(x_fa == NULL); x_fa = new double[nfa][3];

      FOR_IFA {
        // zero the "face-geometry" elements...
        FOR_I3 x_fa[ifa][i] = 0.0;
        FOR_I3 n_fa[ifa][i] = 0.0;
        // loop to get total face normal and approximate x_fa...
        double x_fa_tmp[3] = { 0.0, 0.0, 0.0 };
        const int ino_f = noofa_v_tmp[noofa_i[ifa]];
        double wgt_sum = 0.0;
        int ino1 = noofa_v_tmp[noofa_i[ifa+1]-1];
        for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
          const int ino0 = ino1;
          ino1 = noofa_v_tmp[nof];
          const double wgt = DIST(x_no_tmp[ino0],x_no_tmp[ino1]);
          FOR_I3 x_fa_tmp[i] += wgt*(x_no_tmp[ino0][i]+x_no_tmp[ino1][i]);
          wgt_sum += wgt;
          const double n[3] = TRI_NORMAL_2(x_no_tmp[ino_f],x_no_tmp[ino0],x_no_tmp[ino1]);
          FOR_I3 n_fa[ifa][i] += 0.5*n[i];
        }
        assert(wgt_sum > 0.0);
        FOR_I3 x_fa_tmp[i] /= wgt_sum*2.0;
        // get exact x_fa...
        wgt_sum = 0.0;
        ino1 = noofa_v_tmp[noofa_i[ifa+1]-1];
        for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
          const int ino0 = ino1;
          ino1 = noofa_v_tmp[nof];
          if (build_r_vv) {
            r_vv_tmp[2*ifa  ] = max(r_vv_tmp[2*ifa  ],DIST(x_vv_tmp[2*ifa  ],x_no_tmp[ino1]));
            r_vv_tmp[2*ifa+1] = max(r_vv_tmp[2*ifa+1],DIST(x_vv_tmp[2*ifa+1],x_no_tmp[ino1]));
          }
          const double n[3] = TRI_NORMAL_2(x_fa_tmp,x_no_tmp[ino0],x_no_tmp[ino1]);
          const double wgt = DOT_PRODUCT(n_fa[ifa],n); // area squared weighting !?
          FOR_I3 x_fa[ifa][i] += wgt*(x_fa_tmp[i]+x_no_tmp[ino0][i]+x_no_tmp[ino1][i]);
          wgt_sum += wgt;
        }
        if (wgt_sum > 0.0) {
          FOR_I3 x_fa[ifa][i] /= wgt_sum*3.0;
        }
        else {
          // this must be a zero-area face -- so just pick a node...
          FOR_I3 x_fa[ifa][i] = x_no_tmp[ino1][i];
        }
      }
      delete[] noofa_v_tmp;
      delete[] x_no_tmp; x_no_tmp = NULL;
      delete[] x_vv_tmp; x_vv_tmp = NULL;
      delete[] vol_cv_tmp; vol_cv_tmp = NULL;

      if (build_r_vv) {
        DistributedDataExchanger dde((int8*)cvofa_global,2*nfa,cvora);
        dde.push(r_vv,r_vv_tmp,MAX_DATA);
        delete[] r_vv_tmp; r_vv_tmp = NULL;
      }

      delete[] send_disp;
      delete[] send_count;
      delete[] recv_disp;
      delete[] recv_count;

      // looks like v3 didn't ever compute an x_cv, so lets just set it x_vv..

      assert(x_cv == NULL); x_cv = new double[ncv][3];
      FOR_ICV FOR_I3 x_cv[icv][i] = x_vv[icv][i];

    }

    return 0;

  }

  inline int getFaceZoneForBits(const int bits) {

    // should be just bit pairs...
    assert(bits < ((1<<5)|(1<<3)|(1<<1)));

    switch (bits) {
    case (0):
      return 26;
    case (1<<0):
      return 0;
    case (1<<1):
      return 1;
    case (1<<2):
      return 2;
    case (1<<3):
      return 3;
    case (1<<2)|(1<<0):
      return 4;
    case (1<<2)|(1<<1):
      return 5;
    case (1<<3)|(1<<0):
      return 6;
    case (1<<3)|(1<<1):
      return 7;
    case (1<<4):
      return 8;
    case (1<<5):
      return 9;
    case (1<<4)|(1<<0):
      return 10;
    case (1<<4)|(1<<1):
      return 11;
    case (1<<5)|(1<<0):
      return 12;
    case (1<<5)|(1<<1):
      return 13;
    case (1<<4)|(1<<2):
      return 14;
    case (1<<4)|(1<<3):
      return 15;
    case (1<<5)|(1<<2):
      return 16;
    case (1<<5)|(1<<3):
      return 17;
    case (1<<4)|(1<<2)|(1<<0):
      return 18;
    case (1<<4)|(1<<2)|(1<<1):
      return 19;
    case (1<<4)|(1<<3)|(1<<0):
      return 20;
    case (1<<4)|(1<<3)|(1<<1):
      return 21;
    case (1<<5)|(1<<2)|(1<<0):
      return 22;
    case (1<<5)|(1<<2)|(1<<1):
      return 23;
    case (1<<5)|(1<<3)|(1<<0):
      return 24;
    case (1<<5)|(1<<3)|(1<<1):
      return 25;
    default:
      cout << "getFaceZoneForBits: cannot find bits: " << bits << " ";
      for (int i = 5; i >= 0; --i)
        if (bits & (1<<i))
          cout << "1";
        else
          cout << "0";
      cout << endl;
      throw(0);

      return -1;
    }
  }

  inline int getBitsForFaceZone(const int zone) {

    switch (zone) {
    case (26):
      return (0);
    case (0):
      return (1<<0);
    case (1):
      return (1<<1);
    case (2):
      return (1<<2);
    case (3):
      return (1<<3);
    case (4):
      return (1<<2)|(1<<0);
    case (5):
      return (1<<2)|(1<<1);
    case (6):
      return (1<<3)|(1<<0);
    case (7):
      return (1<<3)|(1<<1);
    case (8):
      return (1<<4);
    case (9):
      return (1<<5);
    case (10):
      return (1<<4)|(1<<0);
    case (11):
      return (1<<4)|(1<<1);
    case (12):
      return (1<<5)|(1<<0);
    case (13):
      return (1<<5)|(1<<1);
    case (14):
      return (1<<4)|(1<<2);
    case (15):
      return (1<<4)|(1<<3);
    case (16):
      return (1<<5)|(1<<2);
    case (17):
      return (1<<5)|(1<<3);
    case (18):
      return (1<<4)|(1<<2)|(1<<0);
    case (19):
      return (1<<4)|(1<<2)|(1<<1);
    case (20):
      return (1<<4)|(1<<3)|(1<<0);
    case (21):
      return (1<<4)|(1<<3)|(1<<1);
    case (22):
      return (1<<5)|(1<<2)|(1<<0);
    case (23):
      return (1<<5)|(1<<2)|(1<<1);
    case (24):
      return (1<<5)|(1<<3)|(1<<0);
    case (25):
      return (1<<5)|(1<<3)|(1<<1);
    default:
      // should be one of the 27 possible zones (0:26)...
      cout << "getBitsForFaceZone: cannot find zone: " << zone << endl;
      throw(0);

      return -1;
    }
  }

  void copy_mles(const int ncopy[3]) {
    if (mpi_rank == 0) cout << "StripedMesh::copy_mles()" << endl;

    // ncopy can be signed (in the future), so stored abs for convenience...
    int nc[3] = {abs(ncopy[0]),abs(ncopy[1]),abs(ncopy[2])};

    // cylindrical periodicities have maximum number of copies...
    int nc_max[3] = {INT_MAX,INT_MAX,INT_MAX}; 
    for (int ipt = 0; ipt < PeriodicData::periodicTransformVec.size(); ++ipt) {
      const int kind = PeriodicData::periodicTransformVec[ipt].getKind();
      if ( (kind == PERIODIC_TRANSFORM_CYL_X) ||
           (kind == PERIODIC_TRANSFORM_CYL_Y) ||
           (kind == PERIODIC_TRANSFORM_CYL_Z) ) {
        const double degrees = abs(PeriodicData::periodicTransformVec[ipt].getData(2));
        nc_max[ipt] = lrint(360.0/degrees);
        if (nc[ipt] > nc_max[ipt]) {
          CWARN(" > number of copies is greater than max allowable by cylindrical periodicity. taking max...");
          nc[ipt] = nc_max[ipt];
        }
      }
    }
    const int ncopies = nc[0]*nc[1]*nc[2];

    if (mpi_rank == 0) {
      if (pbi_pno != NULL) 
        cout << " > creating " << ncopies << " copies (including original), breakdown: " << COUT_VEC(nc) << endl;
      else
        CERR(" > missing PBI_PNO record... please rebuild base mles file.");
    }

    // indexing
    // a[(((ic*nc[1]+jc)*nc[2]+kc)*ncv+icv] = a[ic][jc][kc][icv]

    // --------------
    // cv data...
    // --------------

    // x_cv .. 

    assert(x_cv != NULL);
    double (*x_cv_base)[3] = x_cv;
    x_cv = new double[ncv*ncopies][3]; // include original
    FOR_ICV FOR_I3 x_cv[icv][i] = x_cv_base[icv][i];
    //delete[] x_cv_base;
    int icopy = 0;
    for (int ic = 0; ic < nc[0]; ++ic) {
      for (int jc = 0; jc < nc[1]; ++jc) {
        for (int kc = 0; kc < nc[2]-1; ++kc) {
          //cout << ((ic*nc[1]+jc)*nc[2]+kc) << " " << icopy << endl; cout.flush();
          assert(((ic*nc[1]+jc)*nc[2]+kc) == icopy);
          // perform k copy/transform...
          FOR_ICV FOR_I3 x_cv[((ic*nc[1]+jc)*nc[2]+(kc+1))*ncv+icv][i] = x_cv[((ic*nc[1]+jc)*nc[2]+kc)*ncv+icv][i];
          if (ncopy[2] > 0) {
            PeriodicData::periodicTransformVec[2].translate(&x_cv[((ic*nc[1]+jc)*nc[2]+(kc+1))*ncv],ncv);
          }
          else {
            assert(ncopy[2] < 0);
            PeriodicData::periodicTransformVec[2].inv_translate(&x_cv[((ic*nc[1]+jc)*nc[2]+(kc+1))*ncv],ncv);
          }
          //cout << "2 dir copy" << endl;
          ++icopy;
        }
        if (jc < nc[1]-1) {
          // perform j copy/transform...
          FOR_ICV FOR_I3 x_cv[(ic*nc[1]+(jc+1))*nc[2]*ncv+icv][i] = x_cv[((ic*nc[1]+jc)*nc[2])*ncv+icv][i];
          if (ncopy[1] > 0) {
            PeriodicData::periodicTransformVec[1].translate(&x_cv[(ic*nc[1]+(jc+1))*nc[2]*ncv],ncv);
          }
          else {
            assert(ncopy[1] < 0);
            PeriodicData::periodicTransformVec[1].inv_translate(&x_cv[(ic*nc[1]+(jc+1))*nc[2]*ncv],ncv);
          }
          //cout << "1 dir copy" << endl;
          ++icopy;
        }
      }
      if (ic < nc[0]-1) {
        // perform i copy/transform...
        FOR_ICV FOR_I3 x_cv[(ic+1)*nc[1]*nc[2]*ncv+icv][i] = x_cv[ic*nc[1]*nc[2]*ncv+icv][i];
        if (ncopy[0] > 0) {
          PeriodicData::periodicTransformVec[0].translate(&x_cv[(ic+1)*nc[1]*nc[2]*ncv],ncv);
        }
        else {
          assert(ncopy[0] < 0);
          PeriodicData::periodicTransformVec[0].inv_translate(&x_cv[(ic+1)*nc[1]*nc[2]*ncv],ncv);
        }
        //cout << "0 dir copy" << endl;
        ++icopy;
      }
    }
    assert(icopy == ncopies-1);

    // x_vv .. 

    assert(x_vv != NULL);
    double (*x_vv_base)[3] = x_vv;
    x_vv = new double[ncv*ncopies][3]; 
    FOR_ICV FOR_I3 x_vv[icv][i] = x_vv_base[icv][i];
    delete[] x_vv_base;
    icopy = 0;
    for (int ic = 0; ic < nc[0]; ++ic) {
      for (int jc = 0; jc < nc[1]; ++jc) {
        for (int kc = 0; kc < nc[2]-1; ++kc) {
          assert(((ic*nc[1]+jc)*nc[2]+kc) == icopy);
          FOR_ICV FOR_I3 x_vv[((ic*nc[1]+jc)*nc[2]+(kc+1))*ncv+icv][i] = x_vv[((ic*nc[1]+jc)*nc[2]+kc)*ncv+icv][i];
          if (ncopy[2] > 0) {
            PeriodicData::periodicTransformVec[2].translate(&x_vv[((ic*nc[1]+jc)*nc[2]+(kc+1))*ncv],ncv);
          }
          else {
            assert(ncopy[2] < 0);
            PeriodicData::periodicTransformVec[2].inv_translate(&x_vv[((ic*nc[1]+jc)*nc[2]+(kc+1))*ncv],ncv);
          }
          ++icopy;
        }
        if (jc < nc[1]-1) {
          FOR_ICV FOR_I3 x_vv[(ic*nc[1]+(jc+1))*nc[2]*ncv+icv][i] = x_vv[((ic*nc[1]+jc)*nc[2])*ncv+icv][i];
          if (ncopy[1] > 0) {
            PeriodicData::periodicTransformVec[1].translate(&x_vv[(ic*nc[1]+(jc+1))*nc[2]*ncv],ncv);
          }
          else {
            assert(ncopy[1] < 0);
            PeriodicData::periodicTransformVec[1].inv_translate(&x_vv[(ic*nc[1]+(jc+1))*nc[2]*ncv],ncv);
          }
          ++icopy;
        }
      }
      if (ic < nc[0]-1) {
        FOR_ICV FOR_I3 x_vv[(ic+1)*nc[1]*nc[2]*ncv+icv][i] = x_vv[ic*nc[1]*nc[2]*ncv+icv][i];
        if (ncopy[0] > 0) {
          PeriodicData::periodicTransformVec[0].translate(&x_vv[(ic+1)*nc[1]*nc[2]*ncv],ncv);
        }
        else {
          assert(ncopy[0] < 0);
          PeriodicData::periodicTransformVec[0].inv_translate(&x_vv[(ic+1)*nc[1]*nc[2]*ncv],ncv);
        }
        ++icopy;
      }
    }
    assert(icopy == ncopies-1);

    // vol_cv .. 

    assert(vol_cv != NULL);
    double *vol_cv_base = vol_cv;
    vol_cv = new double[ncv*ncopies]; 
    FOR_ICV vol_cv[icv] = vol_cv_base[icv];
    delete[] vol_cv_base;
    icopy = 0;
    for (int ic = 0; ic < nc[0]; ++ic) {
      for (int jc = 0; jc < nc[1]; ++jc) {
        for (int kc = 0; kc < nc[2]-1; ++kc) {
          assert(((ic*nc[1]+jc)*nc[2]+kc) == icopy);
          FOR_ICV vol_cv[((ic*nc[1]+jc)*nc[2]+(kc+1))*ncv+icv] = vol_cv[((ic*nc[1]+jc)*nc[2]+kc)*ncv+icv];
          ++icopy;
        }
        if (jc < nc[1]-1) {
          FOR_ICV vol_cv[(ic*nc[1]+(jc+1))*nc[2]*ncv+icv] = vol_cv[((ic*nc[1]+jc)*nc[2])*ncv+icv];
          ++icopy;
        }
      }
      if (ic < nc[0]-1) {
        FOR_ICV vol_cv[(ic+1)*nc[1]*nc[2]*ncv+icv] = vol_cv[ic*nc[1]*nc[2]*ncv+icv];
        ++icopy;
      }
    }
    assert(icopy == ncopies-1);

    // r_vv .. 

    assert(r_vv != NULL);
    double *r_vv_base = r_vv;
    r_vv = new double[ncv*ncopies]; 
    FOR_ICV r_vv[icv] = r_vv_base[icv];
    delete[] r_vv_base;
    icopy = 0;
    for (int ic = 0; ic < nc[0]; ++ic) {
      for (int jc = 0; jc < nc[1]; ++jc) {
        for (int kc = 0; kc < nc[2]-1; ++kc) {
          assert(((ic*nc[1]+jc)*nc[2]+kc) == icopy);
          FOR_ICV r_vv[((ic*nc[1]+jc)*nc[2]+(kc+1))*ncv+icv] = r_vv[((ic*nc[1]+jc)*nc[2]+kc)*ncv+icv];
          ++icopy;
        }
        if (jc < nc[1]-1) {
          FOR_ICV r_vv[(ic*nc[1]+(jc+1))*nc[2]*ncv+icv] = r_vv[((ic*nc[1]+jc)*nc[2])*ncv+icv];
          ++icopy;
        }
      }
      if (ic < nc[0]-1) {
        FOR_ICV r_vv[(ic+1)*nc[1]*nc[2]*ncv+icv] = r_vv[ic*nc[1]*nc[2]*ncv+icv];
        ++icopy;
      }
    }
    assert(icopy == ncopies-1);

    const int ncv0 = ncv;
    ncv *= ncopies;
    const int8 ncv_global0 = ncv_global; 
    ncv_global *= ncopies;
    int8* cvora0 = cvora; cvora = NULL; 
    MiscUtils::buildXora(cvora,ncv); 
    assert(cvora[mpi_rank+1]-cvora[mpi_rank] < TWO_BILLION);
    assert(ncv == cvora[mpi_rank+1]-cvora[mpi_rank]);
    assert(ncv_global == cvora[mpi_size]);
    assert(ncv0 == cvora0[mpi_rank+1]-cvora0[mpi_rank]);
    assert(ncv_global0 == cvora0[mpi_size]);

    // ----------
    // fa data...
    // ----------

    // cvofa_global...

    assert(cvofa_global != NULL);
    int8 (*cvofa_global_base)[2] = cvofa_global;
    vector<pair<int8,int8> > cvofa_global_vec; cvofa_global_vec.reserve(nfa*ncopies);
    vector<pair<int,int> > fazone_ifa0_vec; fazone_ifa0_vec.reserve(nfa*ncopies); 
    vector<pair<int,int> > ifa_icopy_vec; ifa_icopy_vec.reserve(nfa*ncopies);
    int nfa_new = 0;
    icopy = 0;
    int ijkc[3] = {0,0,0};
    int nfa_p = 0;
    for (ijkc[0] = 0; ijkc[0] < nc[0]; ++ijkc[0]) {
      for (ijkc[1] = 0; ijkc[1] < nc[1]; ++ijkc[1]) {
        for (ijkc[2] = 0; ijkc[2] < nc[2]; ++ijkc[2]) {
          assert(((ijkc[0]*nc[1]+ijkc[1])*nc[2]+ijkc[2]) == icopy);
          FOR_IFA {
            // TODO accomodate backwards...
            const int8 icv0_global = cvofa_global_base[ifa][0]; assert((icv0_global >= 0)&&(icv0_global < ncv_global0));
            const int8 icv1_global = (cvofa_global_base[ifa][1]&MASK_55BITS); assert((icv1_global >= 0)&&(icv1_global < ncv_global0));
            const int bits = (cvofa_global_base[ifa][1]>>55);
            const int rank0 = MiscUtils::getRankInXora(icv0_global,cvora0);
            const int rank1 = MiscUtils::getRankInXora(icv1_global,cvora0);
            if (bits) {
              // NOTES: either flip bits or external boundary of full copy
              int ijk_copy[3] = {0,0,0};
              int bits_new = 0;
              bool outer_face = false;
              int count = 0;
              int fcount = 0;
              int bcount = 0;
              for (int bit_pair = 0; bit_pair < 3; ++bit_pair) {
                if (bits & (1<<(2*bit_pair))) {
                  assert((bits & (1<<(2*bit_pair+1))) == 0);
                  ++count;
                  // forward bit...
                  ++fcount;
                  // take periodic copy and copy over bits...
                  if (ijkc[bit_pair] == nc[bit_pair]-1) {
                    ijk_copy[bit_pair] = 0;
                    if (nc[bit_pair] < nc_max[bit_pair]) {
                      outer_face = true;
                      bits_new |= (bits & (1<<(2*bit_pair)));
                    }
                  }
                  else {
                    // take next and clear bits...
                    ijk_copy[bit_pair] = ijkc[bit_pair]+1;
                  }
                }
                else if (bits & (1<<(2*bit_pair+1))) {
                  assert((bits & (1<<(2*bit_pair))) == 0);
                  ++count;
                  // backward bit...
                  ++bcount;
                  if (ijkc[bit_pair] == 0) {
                    // take periodic copy and copy over bits...
                    ijk_copy[bit_pair] = nc[bit_pair]-1;
                    if (nc[bit_pair] < nc_max[bit_pair]) {
                      outer_face = true;
                      bits_new |= (bits & (1<<(2*bit_pair+1)));
                    }
                  }
                  else { 
                    // take previous and clear bits...
                    ijk_copy[bit_pair] = ijkc[bit_pair]-1;
                  }
                }
                else {
                  // no bits in this direction, so just take current copy...
                  ijk_copy[bit_pair] = ijkc[bit_pair];
                }
              }
              assert(outer_face||(bits_new == 0));
              if ((outer_face)||
                  ((rank0 <= rank1)&&((fcount == count)||((fcount > 0)&&(fcount < count)&&(bits > BitUtils::flipPeriodicBits(bits))))) || 
                  ((rank0 <  rank1)&&((bcount == count)||((bcount > 0)&&(bcount < count)&&(bits < BitUtils::flipPeriodicBits(bits)))))) { 
                const int icopy2 = (ijk_copy[0]*nc[1]+ijk_copy[1])*nc[2]+ijk_copy[2];
                const int8 icv0_global_new = icopy*(cvora0[rank0+1]-cvora0[rank0]) + icv0_global - cvora0[rank0] + cvora[rank0];
                const int8 icv1_global_new = icopy2*(cvora0[rank1+1]-cvora0[rank1]) + icv1_global - cvora0[rank1] + cvora[rank1];
                fazone_ifa0_vec.push_back(pair<int,int>(getFaceZoneForBits(bits_new),nfa_new));
                cvofa_global_vec.push_back(pair<int8,int8>(icv0_global_new,icv1_global_new));
                ifa_icopy_vec.push_back(pair<int,int>(ifa,icopy));
                ++nfa_new;
                if (bits_new)
                  ++nfa_p;
              }
            }
            else {
              // just copy over and shift cv global values...
              const int bits_new = 0;
              const int8 icv0_global_new = icopy*(cvora0[rank0+1]-cvora0[rank0]) + icv0_global - cvora0[rank0] + cvora[rank0];
              const int8 icv1_global_new = icopy*(cvora0[rank1+1]-cvora0[rank1]) + icv1_global - cvora0[rank1] + cvora[rank1];
              assert(icv0_global_new < icv1_global_new); 
              fazone_ifa0_vec.push_back(pair<int,int>(getFaceZoneForBits(bits_new),nfa_new));
              cvofa_global_vec.push_back(pair<int8,int8>(icv0_global_new,icv1_global_new));
              ifa_icopy_vec.push_back(pair<int,int>(ifa,icopy));
              ++nfa_new;
            }
          }
          ++icopy;
        }
      }
    }
    delete[] x_cv_base;
    assert(fazone_ifa0_vec.size() == nfa_new);
    assert(cvofa_global_vec.size() == nfa_new);

    // sort faces on face zone...  
    sort(fazone_ifa0_vec.begin(),fazone_ifa0_vec.end());

    // x_fa and n_fa...

    cvofa_global = new int8[nfa_new][2];
    double (*x_fa_base)[3] = x_fa;
    x_fa = new double[nfa_new][3];
    double (*n_fa_base)[3] = n_fa;
    n_fa = new double[nfa_new][3];
    for (int ifa_new = 0; ifa_new < nfa_new; ++ifa_new) {
      const int bits = getBitsForFaceZone(fazone_ifa0_vec[ifa_new].first);
      const int ifa0 = fazone_ifa0_vec[ifa_new].second;
      cvofa_global[ifa_new][0] = cvofa_global_vec[ifa0].first;
      cvofa_global[ifa_new][1] = ( cvofa_global_vec[ifa0].second | (int8(bits)<<55) );
      const int ifa     = ifa_icopy_vec[ifa0].first;
      const int icopy   = ifa_icopy_vec[ifa0].second;
      int ijkc[3];
      ijkc[2] = icopy%nc[2]; assert((ijkc[2] >= 0)&&(ijkc[2] < nc[2]));
      ijkc[1] = (icopy/nc[2])%nc[1]; assert((ijkc[1] >= 0)&&(ijkc[1] < nc[1])); 
      ijkc[0] = icopy/(nc[2]*nc[1]); assert((ijkc[0] >= 0)&&(ijkc[0] < nc[0]));
      assert(icopy == ((ijkc[0]*nc[1]+ijkc[1])*nc[2]+ijkc[2]));
      FOR_I3 x_fa[ifa_new][i] = x_fa_base[ifa][i];
      FOR_I3 {
        for (int ii = 0; ii < ijkc[i]; ++ii) {
          if (ncopy[i] > 0) 
            PeriodicData::periodicTransformVec[i].translate(x_fa[ifa_new]);
          else if (ncopy[i] < 0) 
            PeriodicData::periodicTransformVec[i].inv_translate(x_fa[ifa_new]);
        }
      }
      FOR_I3 n_fa[ifa_new][i] = n_fa_base[ifa][i];
      FOR_I3 {
        for (int ii = 0; ii < ijkc[i]; ++ii) {
          if (ncopy[i] > 0) 
            PeriodicData::periodicTransformVec[i].rotate(n_fa[ifa_new]);
          else if (ncopy[i] < 0) 
            PeriodicData::periodicTransformVec[i].inv_rotate(n_fa[ifa_new]);
        }
      }
    }
    delete[] cvofa_global_base; cvofa_global_base = NULL;
    delete[] x_fa_base; x_fa_base = NULL;
    delete[] n_fa_base; n_fa_base = NULL;
    cvofa_global_vec.clear();

    /*
    //FILE * fp = fopen("xfa.csv","w");
    double min_dist = HUGE_VAL;
    for (int ifa_new = 0; ifa_new < nfa_new; ++ifa_new) {
    const int bits = getBitsForFaceZone(fazone_ifa0_vec[ifa_new].first);
    const int ifa0 = fazone_ifa0_vec[ifa_new].second;
    const int icopy   = ifa_icopy_vec[ifa0].second;
    //fprintf(fp,"%18.15le,%18.15le,%18.15le,%d,%d\n",x_fa[ifa_new][0],x_fa[ifa_new][1],x_fa[ifa_new][2],icopy,bits);
    for (int ii = ifa_new+1; ii < nfa_new; ++ii) {
    const double dist = DIST(x_fa[ifa_new],x_fa[ii]);
    min_dist = min(dist,min_dist);
    if (dist < 1.0E-6) {
    const int bits2 = getBitsForFaceZone(fazone_ifa0_vec[ii].first);
    const int ifa02 = fazone_ifa0_vec[ii].second;
    const int icopy2   = ifa_icopy_vec[ifa02].second;
    cout << "SAME FACE: " << ifa_new << " " << ii << " " << icopy << " " << icopy2 << " " << COUT_VEC(x_fa[ifa_new]) << endl;
    BitUtils::dumpBits(bits,6);
    BitUtils::dumpBits(bits2,6);
    }
    }
    }
    cout << " MIN:  " << min_dist << endl;
    //fclose(fp);
    */

    // noofa_i/v_global, x_no...

    // get bits and indices for nodes from faces...
    
    map<const int8,int> globalMap;
    int nno_tmp = 0;
    int disp_zone[27+1];
    disp_zone[0] = 0;
    for (int izone = 0; izone < 27; ++izone) 
      disp_zone[izone+1] = disp_zone[izone] + nfa_zone[izone];
    assert(disp_zone[27] == nfa_global);
    int izone = 0;
    while (disp_zone[izone+1] < faora[mpi_rank]) 
      ++izone;
    FOR_IFA {
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int8 ino_global = noofa_v_global[nof];
        assert((ino_global >= 0)&&(ino_global < nno_global));
        map<const int8,int>::iterator iter = globalMap.find(ino_global);
        if (iter == globalMap.end()) 
          globalMap[ino_global] = nno_tmp++;
      }
    }
    
    int *no_bits_tmp = new int[nno_tmp];
    int8 *ino_global_tmp = new int8[nno_tmp];
    int *noofa_v_tmp = new int[noofa_i[nfa]];
    for (int ino_tmp = 0; ino_tmp < nno_tmp; ++ino_tmp) no_bits_tmp[ino_tmp] = 0;
    izone = 0;
    while (disp_zone[izone+1] <= faora[mpi_rank]) 
      ++izone;
    FOR_IFA {
      const int8 ifa_global = ifa + faora[mpi_rank];
      while (disp_zone[izone+1] <= ifa_global) 
        ++izone;
      int bits = getBitsForFaceZone(izone);
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int8 ino_global = noofa_v_global[nof];
        map<const int8,int>::iterator iter = globalMap.find(ino_global);
        assert(iter != globalMap.end());
        no_bits_tmp[iter->second] |= bits;
        noofa_v_tmp[nof] = iter->second;
        ino_global_tmp[iter->second] = ino_global;
      }
    }
    delete[] noofa_v_global; noofa_v_global = NULL; // replaced with noofa_v_tmp 
    globalMap.clear();

    //  =================================================
    //  step 1: 
    //  we need to reduce no_bits_tmp onto nodal striping
    //  =================================================

    int * send_count = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;
    for (int ino_tmp = 0; ino_tmp < nno_tmp; ++ino_tmp) {
      const int rank = MiscUtils::getRankInXora(ino_global_tmp[ino_tmp],noora);
      send_count[rank] += 3; // ino_temp, ino (striped), bits
    }

    int * send_disp = new int[mpi_size];
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
    
    int * send_buf_int = new int[send_count_sum];
    for (int ino_tmp = 0; ino_tmp < nno_tmp; ++ino_tmp) {
      const int rank = MiscUtils::getRankInXora(ino_global_tmp[ino_tmp],noora);
      send_buf_int[send_disp[rank]++] = ino_tmp;
      send_buf_int[send_disp[rank]++] = ino_global_tmp[ino_tmp]-noora[rank];
      send_buf_int[send_disp[rank]++] = no_bits_tmp[ino_tmp];
    }
    delete[] no_bits_tmp;
    delete[] ino_global_tmp; ino_global_tmp = NULL; // reused for noobf

    //  rewind
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

    //  exchange...
    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    
    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    
    int * recv_buf_int = new int[recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
		  recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; send_buf_int = NULL;

    int *no_bits = new int[nno];
    FOR_INO no_bits[ino] = 0;
    for (int rank = 0; rank < mpi_size; ++rank) {
      for (int irecv = recv_disp[rank]; irecv < recv_disp[rank]+recv_count[rank]; irecv += 3) {
        //const int ino_tmp = recv_buf_int[irecv  ]; // use later for send back to face striping
        const int ino     = recv_buf_int[irecv+1];
        const int bits    = recv_buf_int[irecv+2];
        no_bits[ino] |= bits;
      }
    }
    /*
    FOR_INO {
      if (ino < nno_p) {
        if (no_bits[ino] == 0) {
          int bits;
          int8 ino_parent_global;
          BitUtils::unpackPbiHash(bits,ino_parent_global,pbi_pno[ino]);
          BitUtils::dumpBits(bits,6);
          cout << COUT_VEC(x_no[ino]) << endl;
        }
      }
      else {
        assert(no_bits[ino] == 0);
      }
    }
    MPI_Pause("check");
    */
    
    // we have pbi_pno for first nno_p nodes...
    for (int ino = 0; ino < nno_p; ++ino) {
      //assert(no_bits[ino]);
      int bits;
      int8 ino_parent_global;
      BitUtils::unpackPbiHash(bits,ino_parent_global,pbi_pno[ino]);
      if (bits == 0) 
        assert(ino_parent_global == ino+noora[mpi_rank]);
      else 
        assert(ino_parent_global < ino+noora[mpi_rank]);
    }
    
    int nparent = 0;
    assert(globalMap.empty());
    vector<int8> parentVec;
    for (int ino = 0; ino < nno_p; ++ino) {
      //assert(no_bits[ino]);
      int8 ino_global = ino + noora[mpi_rank];
      int bits;
      int8 ino_parent_global;
      BitUtils::unpackPbiHash(bits,ino_parent_global,pbi_pno[ino]);
      assert(ino_parent_global <= ino_global);
      map<const int8,int>::iterator iter = globalMap.find(ino_parent_global);
      if (iter == globalMap.end()) {
        parentVec.push_back(ino_parent_global);
        globalMap[ino_parent_global] = nparent++; 
      }
    }
    assert(parentVec.size() == nparent);

    // exchange parents bits to build bits of parent sets for the parents in your striping...
    // we know that ino_parent_global <= ino_global and ino_parent_global == ino_global when its on rank

    int * send_count_3 = new int[mpi_size];
    int * send_disp_3 = new int[mpi_size];
    FOR_RANK send_count_3[rank] = 0;
    int* send_buf_int_3 = NULL;
    int send_count_sum_3 = 0;
    map<const int,int8> * parentBitMap = new map<const int,int8>[nparent]; // bits and ino_global...
    for (int iter = 0; iter < 2; ++iter) {
      for (int ino = 0; ino < nno_p; ++ino) {
        //assert(no_bits[ino]);
        int8 ino_global = ino + noora[mpi_rank];
        int bits;
        int8 ino_parent_global;
        BitUtils::unpackPbiHash(bits,ino_parent_global,pbi_pno[ino]);
        assert(ino_parent_global <= ino_global);
        if (iter == 0) {
          if (ino_parent_global >= noora[mpi_rank]) {
            // parent is on this rank, so add bits to set. could get more from higher ranks...
            map<const int8,int>::iterator iter2 = globalMap.find(ino_parent_global);
            assert(iter2 != globalMap.end());
            parentBitMap[iter2->second].insert(pair<const int,int8>(no_bits[ino],ino_global));
          }
          else {
            // parent is on a lower rank, so prep to send its id and bits to the lower rank...
            const int rank = MiscUtils::getRankInXora(ino_parent_global,noora);
            assert(rank < mpi_rank);
            send_count_3[rank] += 3; // ino_parent on rank, bits and ino on my rank
          }
        }
        else if (ino_parent_global < noora[mpi_rank]) {
          // parent is on a lower rank...
          const int rank = MiscUtils::getRankInXora(ino_parent_global,noora);
          assert(rank < mpi_rank);
          const int ino_parent = ino_parent_global-noora[rank];
          send_buf_int_3[send_disp_3[rank]++] = ino_parent;
          send_buf_int_3[send_disp_3[rank]++] = no_bits[ino];
          send_buf_int_3[send_disp_3[rank]++] = ino;
        }
      }

      // rewind...
      send_disp_3[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp_3[rank] = send_count_3[rank-1] + send_disp_3[rank-1];

      if (iter == 0) {
        send_count_sum_3 = send_disp_3[mpi_size-1] + send_count_3[mpi_size-1];
        assert(send_count_sum_3%3 == 0);
        send_buf_int_3 = new int[send_count_sum_3];
      }
    }

    // exchange...
    int * recv_count_3 = new int[mpi_size];
    MPI_Alltoall(send_count_3,1,MPI_INT,recv_count_3,1,MPI_INT,mpi_comm);

    int * recv_disp_3 = new int[mpi_size];
    recv_disp_3[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp_3[rank] = recv_count_3[rank-1] + recv_disp_3[rank-1];
    int recv_count_sum_3 = recv_disp_3[mpi_size-1] + recv_count_3[mpi_size-1];

    assert(recv_count_sum_3%3 == 0);
    int * recv_buf_int_3 = new int[recv_count_sum_3];
    MPI_Alltoallv(send_buf_int_3,send_count_3,send_disp_3,MPI_INT,
                  recv_buf_int_3,recv_count_3,recv_disp_3,MPI_INT,mpi_comm);
    delete[] send_buf_int_3; send_buf_int_3 = NULL;

    // set and prep for send back...
    set<int> * parentRankSet = new set<int>[nparent]; // should be small (1-8) so just leave as set
    FOR_RANK {
      for (int irecv = recv_disp_3[rank]; irecv < (recv_disp_3[rank]+recv_count_3[rank]); irecv += 3) {
        const int ino_parent = recv_buf_int_3[irecv];
        const int8 ino_parent_global = ino_parent + noora[mpi_rank];
        map<const int8,int>::iterator iter2 = globalMap.find(ino_parent_global);
        assert(iter2 != globalMap.end());
        // then add in bits from recv buf...
        assert(!parentBitMap[iter2->second].empty()); // it should of at least had one element since the parent lives on this rank
        const int bits = recv_buf_int_3[irecv+1];
        const int8 ino_global = recv_buf_int_3[irecv+2]+noora[rank];
        parentRankSet[iter2->second].insert(rank);
        parentBitMap[iter2->second].insert(pair<const int,int8>(bits,ino_global));
      }
    }
    delete[] recv_buf_int_3; recv_buf_int_3 = NULL;

    FOR_RANK send_count_3[rank] = 0;
    int8* send_buf_int8_3 = NULL;
    for (int iter = 0; iter < 2; ++iter) {
      for (int iparent = 0; iparent < nparent; ++iparent) {
        const int ino_parent = parentVec[iparent]-noora[mpi_rank];
        for (set<int>::iterator iter2 = parentRankSet[iparent].begin(); iter2 != parentRankSet[iparent].end(); ++iter2) {
          const int rank = *iter2; assert(rank > mpi_rank);
          if (iter == 0) {
            send_count_3[rank] += 2*parentBitMap[iparent].size(); // iparent and bits. ino_global in int8 buf
          }
          else {
            // you could bring out iparent and include the size, but this doesn't seem worth it when the bits per parent are small.
            for (map<const int,int8>::iterator iter3 = parentBitMap[iparent].begin(); iter3 != parentBitMap[iparent].end(); ++iter3) {
              const int bits = iter3->first;
              send_buf_int_3[send_disp_3[rank]  ] = ino_parent;
              send_buf_int_3[send_disp_3[rank]+1] = bits;
              const int ino_global = iter3->second;
              send_buf_int8_3[send_disp_3[rank]/2] = ino_global;
              send_disp_3[rank] += 2;
            }
          }
        }
      }

      // rewind...
      send_disp_3[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp_3[rank] = send_count_3[rank-1] + send_disp_3[rank-1];

      if (iter == 0) {
        send_count_sum_3 = send_disp_3[mpi_size-1] + send_count_3[mpi_size-1];
        assert(send_count_sum_3%2 == 0);
        send_buf_int_3 = new int[send_count_sum_3];
        send_buf_int8_3 = new int8[send_count_sum_3/2];
      }
    }
    parentVec.clear();
    delete[] parentRankSet;

    // exchange...
    MPI_Alltoall(send_count_3,1,MPI_INT,recv_count_3,1,MPI_INT,mpi_comm);

    recv_disp_3[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp_3[rank] = recv_count_3[rank-1] + recv_disp_3[rank-1];
    recv_count_sum_3 = recv_disp_3[mpi_size-1] + recv_count_3[mpi_size-1];

    assert(recv_count_sum_3%2 == 0);
    assert(recv_buf_int_3 == NULL); recv_buf_int_3 = new int[recv_count_sum_3];
    MPI_Alltoallv(send_buf_int_3,send_count_3,send_disp_3,MPI_INT,
                  recv_buf_int_3,recv_count_3,recv_disp_3,MPI_INT,mpi_comm);
    delete[] send_buf_int_3; send_buf_int_3 = NULL;

    FOR_RANK {
      send_count_3[rank] /= 2;
      send_disp_3[rank] /= 2;
      recv_count_3[rank] /= 2;
      recv_disp_3[rank] /= 2;
    }

    int8* recv_buf_int8_3 = new int8[recv_count_sum_3/2];
    MPI_Alltoallv(send_buf_int8_3,send_count_3,send_disp_3,MPI_INT8,
                  recv_buf_int8_3,recv_count_3,recv_disp_3,MPI_INT8,mpi_comm);
    delete[] send_buf_int8_3; send_buf_int8_3 = NULL;

    FOR_RANK {
      for (int irecv = recv_disp_3[rank]; irecv < (recv_disp_3[rank]+recv_count_3[rank]); ++irecv) {
        const int8 ino_parent_global = recv_buf_int_3[2*irecv] + noora[rank];
        const int bits = recv_buf_int_3[2*irecv+1];
        const int8 ino_global = recv_buf_int8_3[irecv];
        map<const int8,int>::iterator iter = globalMap.find(ino_parent_global);
        assert(iter != globalMap.end());
        // then add in bits amd ino_global from recv bufs...
        parentBitMap[iter->second].insert(pair<const int,int8>(bits,ino_global));
      }
    }
    delete[] recv_buf_int_3; recv_buf_int_3 = NULL;
    delete[] recv_buf_int8_3; recv_buf_int_3 = NULL;

    vector<pair<int,int> > ino_icopy_vec; ino_icopy_vec.reserve(nno*ncopies);
    vector<pair<int,int> > ino_icopy_skip_vec; ino_icopy_skip_vec.reserve(nno*ncopies);
    icopy = 0;
    for (ijkc[0] = 0; ijkc[0] < nc[0]; ++ijkc[0]) {
      for (ijkc[1] = 0; ijkc[1] < nc[1]; ++ijkc[1]) {
        for (ijkc[2] = 0; ijkc[2] < nc[2]; ++ijkc[2]) {
          assert(((ijkc[0]*nc[1]+ijkc[1])*nc[2]+ijkc[2]) == icopy);
          // only certain periodic nodes get copied...
          for (int ino = 0; ino < nno_p; ++ino) {
            //assert(no_bits[ino]);
            int bits;
            int8 ino_parent_global;
            BitUtils::unpackPbiHash(bits,ino_parent_global,pbi_pno[ino]);
            map<const int8,int>::iterator iter = globalMap.find(ino_parent_global);
            assert(iter != globalMap.end());
            // all parents get copied...
            if (bits == 0) { 
              assert(ino_parent_global == ino+noora[mpi_rank]);
              ino_icopy_vec.push_back(pair<int,int>(ino,icopy));
            }
            // children whose appropriately flipped bits are not in its parent set are copied...
            else {
              bool copy = true;
              // 1d...
              for (int bit_pair = 0; bit_pair < 3; ++bit_pair) {
                if ((ijkc[bit_pair] > 0)||(nc[bit_pair] == nc_max[bit_pair])) {
                  if (no_bits[ino] & (1<<(2*bit_pair+1))) {
                    // if the forward flip of this bit pair has the same parent then we will skip...
                    int bits0 = no_bits[ino];
                    bits0 &= ~(1<<(2*bit_pair+1));
                    bits0 |=  (1<<(2*bit_pair  ));
                    if (parentBitMap[iter->second].find(bits0) != parentBitMap[iter->second].end()) {
                      ino_icopy_skip_vec.push_back(pair<int,int>(ino,icopy));
                      copy = false;
                      break;
                    }
                  }
                }
              }
              // 2d...
              if (copy) {
                for (int bit_pair = 0; bit_pair < 3; ++bit_pair) {
                  if ((nc[bit_pair] > 1)&&(nc[(bit_pair+1)%3] > 1)) {
                    const int bit_pair0 = min(bit_pair,(bit_pair+1)%3);
                    const int bit_pair1 = max(bit_pair,(bit_pair+1)%3);
                    if ((ijkc[bit_pair1] > 0)&&(no_bits[ino]&(1<<(2*bit_pair1+1)))) {
                      if ((ijkc[bit_pair0] < nc[bit_pair0]-1)&&(no_bits[ino]&(1<<(2*bit_pair0)))) {
                        int bits0 = no_bits[ino];
                        bits0 &= ~(1<<(2*bit_pair1+1));
                        bits0 |=  (1<<(2*bit_pair1  ));
                        bits0 &= ~(1<<(2*bit_pair0  ));
                        bits0 |=  (1<<(2*bit_pair0+1));
                        if (parentBitMap[iter->second].find(bits0) != parentBitMap[iter->second].end()) {
                          ino_icopy_skip_vec.push_back(pair<int,int>(ino,icopy));
                          copy = false;
                          break;
                        }
                      }
                      if ((ijkc[bit_pair0] > 0)&&(no_bits[ino]&(1<<(2*bit_pair0+1)))) {
                        int bits0 = no_bits[ino];
                        bits0 &= ~(1<<(2*bit_pair1+1));
                        bits0 |=  (1<<(2*bit_pair1  ));
                        bits0 &= ~(1<<(2*bit_pair0+1));
                        bits0 |=  (1<<(2*bit_pair0  ));
                        if (parentBitMap[iter->second].find(bits0) != parentBitMap[iter->second].end()) {
                          ino_icopy_skip_vec.push_back(pair<int,int>(ino,icopy));
                          copy = false;
                          break;
                        }
                      }
                    }
                  }
                }
              }
              // 3d...
              // I don' think this is necessary for our packings, but I am leaving it her for now.
              if (copy) {
                if ((nc[0] > 1)&&(nc[1] > 1)&&(nc[2] > 1)) {
                  if ((ijkc[0] > 0)&&(no_bits[ino]&(1<<(2*0+1)))&&
                      (ijkc[1] > 0)&&(no_bits[ino]&(1<<(2*1+1)))&&
                      (ijkc[2] > 0)&&(no_bits[ino]&(1<<(2*2+1)))) {
                    int bits0 = no_bits[ino];
                    bits0 &= ~(1<<(2*0+1));
                    bits0 |=  (1<<(2*0  ));
                    bits0 &= ~(1<<(2*1+1));
                    bits0 |=  (1<<(2*1  ));
                    bits0 &= ~(1<<(2*2+1));
                    bits0 |=  (1<<(2*2  ));
                    if (parentBitMap[iter->second].find(bits0) != parentBitMap[iter->second].end()) {
                      ino_icopy_skip_vec.push_back(pair<int,int>(ino,icopy));
                      copy = false;
                    }
                  }
                }
              }

              if (copy) 
                ino_icopy_vec.push_back(pair<int,int>(ino,icopy));
            }
          }
          // all internal nodes get copied...
          for (int ino = nno_p; ino < nno; ++ino) {
            assert(no_bits[ino] == 0); 
            ino_icopy_vec.push_back(pair<int,int>(ino,icopy));
          }
          ++icopy;
        }
      }
    }
    assert((ino_icopy_vec.size()+ino_icopy_skip_vec.size())==(nno*ncopies));
    const int nno_new = ino_icopy_vec.size();
    int8 nno_new_int8 = nno_new;
    int8 nno_new_global;
    MPI_Allreduce(&nno_new_int8,&nno_new_global,1,MPI_INT8,MPI_SUM,mpi_comm);

    // we need to store a count for each ino, so lets group the ino_news on ino
 
    sort(ino_icopy_vec.begin(),ino_icopy_vec.end());
    int *nnono_i = new int[nno+1]; // new node of node
    nnono_i[0] = 0;
    FOR_INO nnono_i[ino+1] = 0;
    for (int ino_new = 0; ino_new < nno_new; ++ino_new) {
      const int ino = ino_icopy_vec[ino_new].first;
      ++nnono_i[ino+1];
    }
    FOR_INO nnono_i[ino+1] += nnono_i[ino];

    // reindex ino_global...

    int8* noora_new = NULL;
    MiscUtils::buildXora(noora_new,nno_new); assert(noora_new[mpi_rank+1]-noora_new[mpi_rank] == nno_new);
    assert(nno_new_global = noora_new[mpi_size]);
    //int8 nno_global_new = noora_new[mpi_size];

    // now rebuild x_no's...
    
    double (*x_no_base)[3] = x_no;
    x_no = new double[nno_new][3];
    for (int ino_new = 0; ino_new < nno_new; ++ino_new) {
      const int ino = ino_icopy_vec[ino_new].first;
      const int icopy = ino_icopy_vec[ino_new].second;
      int ijkc[3];
      ijkc[2] = icopy%nc[2]; assert((ijkc[2] >= 0)&&(ijkc[2] < nc[2]));
      ijkc[1] = (icopy/nc[2])%nc[1]; assert((ijkc[1] >= 0)&&(ijkc[1] < nc[1])); 
      ijkc[0] = icopy/(nc[2]*nc[1]); assert((ijkc[0] >= 0)&&(ijkc[0] < nc[0]));
      assert(icopy == ((ijkc[0]*nc[1]+ijkc[1])*nc[2]+ijkc[2]));
      FOR_I3 x_no[ino_new][i] = x_no_base[ino][i];
      FOR_I3 {
        for (int ii = 0; ii < ijkc[i]; ++ii) {
          if (ncopy[i] > 0) 
            PeriodicData::periodicTransformVec[i].translate(x_no[ino_new]);
          else if (ncopy[i] < 0) 
            PeriodicData::periodicTransformVec[i].inv_translate(x_no[ino_new]);
        }
      }
    }

    // get the nodes living on other ranks/copies...

    sort(ino_icopy_skip_vec.begin(),ino_icopy_skip_vec.end());
    const int nno_new2 = ino_icopy_skip_vec.size();
    int *nnono2_i = new int[nno+1]; // new node of node
    nnono2_i[0] = 0;
    FOR_INO nnono2_i[ino+1] = 0;
    for (int ino_new = 0; ino_new < nno_new2; ++ino_new) {
      const int ino = ino_icopy_skip_vec[ino_new].first;
      ++nnono2_i[ino+1];
    }
    FOR_INO nnono2_i[ino+1] += nnono2_i[ino];

    vector<pair<int8,int> > ino_global_icopy2_vec; ino_global_icopy2_vec.resize(nno_new2);
    for (int ii = 0; ii < nno_new2; ++ii) {
      const int ino = ino_icopy_skip_vec[ii].first;
      const int icopy = ino_icopy_skip_vec[ii].second;
      int bits;
      int8 ino_parent_global;
      BitUtils::unpackPbiHash(bits,ino_parent_global,pbi_pno[ino]);
      map<const int8,int>::iterator iter = globalMap.find(ino_parent_global);
      assert(iter != globalMap.end());
      //const int iparent = iter->second;
      int ijkc[3];
      ijkc[2] = icopy%nc[2]; assert((ijkc[2] >= 0)&&(ijkc[2] < nc[2]));
      ijkc[1] = (icopy/nc[2])%nc[1]; assert((ijkc[1] >= 0)&&(ijkc[1] < nc[1])); 
      ijkc[0] = icopy/(nc[2]*nc[1]); assert((ijkc[0] >= 0)&&(ijkc[0] < nc[0]));
      assert(icopy == ((ijkc[0]*nc[1]+ijkc[1])*nc[2]+ijkc[2]));
      map<const int8,int> candidateMap; // ino_global -> icopy2
      // 1d...
      for (int bit_pair = 0; bit_pair < 3; ++bit_pair) {
        if ((ijkc[bit_pair] > 0)||(nc[bit_pair] == nc_max[bit_pair])) {
          if (no_bits[ino] & (1<<(2*bit_pair+1))) {
            // if the forward flip of this bit pair has the same parent then we will skip...
            int bits0 = no_bits[ino];
            bits0 &= ~(1<<(2*bit_pair+1));
            bits0 |=  (1<<(2*bit_pair  ));
            map<const int,int8>::iterator iter2 = parentBitMap[iter->second].find(bits0);
            if (iter2 != parentBitMap[iter->second].end()) {
              const int8 ino_global = iter2->second;
              int ijkc2[3];
              ijkc2[2] = icopy%nc[2]; assert((ijkc[2] >= 0)&&(ijkc[2] < nc[2]));
              ijkc2[1] = (icopy/nc[2])%nc[1]; assert((ijkc[1] >= 0)&&(ijkc[1] < nc[1])); 
              ijkc2[0] = icopy/(nc[2]*nc[1]); assert((ijkc[0] >= 0)&&(ijkc[0] < nc[0]));
              if (ijkc[bit_pair] > 0)
                ijkc2[bit_pair] -= 1; // forward bit has it so the copy is the previous in this direction
              else 
                ijkc2[bit_pair] = nc[bit_pair]-1; // wrap around
              int icopy2 = ((ijkc2[0]*nc[1]+ijkc2[1])*nc[2]+ijkc2[2]);
              assert(candidateMap.find(ino_global) == candidateMap.end());
              candidateMap[ino_global] = icopy2;
            }
          }
        }
      }
      // 2d...
      for (int bit_pair = 0; bit_pair < 3; ++bit_pair) {
        if ((nc[bit_pair] > 1)&&(nc[(bit_pair+1)%3] > 1)) {
          const int bit_pair0 = min(bit_pair,(bit_pair+1)%3);
          const int bit_pair1 = max(bit_pair,(bit_pair+1)%3);
          if ((ijkc[bit_pair1] > 0)&&(no_bits[ino]&(1<<(2*bit_pair1+1)))) {
            if ((ijkc[bit_pair0] < nc[bit_pair0]-1)&&(no_bits[ino]&(1<<(2*bit_pair0)))) {
              int bits0 = no_bits[ino];
              bits0 &= ~(1<<(2*bit_pair1+1));
              bits0 |=  (1<<(2*bit_pair1  ));
              bits0 &= ~(1<<(2*bit_pair0  ));
              bits0 |=  (1<<(2*bit_pair0+1));
              map<const int,int8>::iterator iter2 = parentBitMap[iter->second].find(bits0);
              if (iter2 != parentBitMap[iter->second].end()) {
                const int8 ino_global = iter2->second;
                int ijkc2[3];
                ijkc2[2] = icopy%nc[2]; assert((ijkc[2] >= 0)&&(ijkc[2] < nc[2]));
                ijkc2[1] = (icopy/nc[2])%nc[1]; assert((ijkc[1] >= 0)&&(ijkc[1] < nc[1])); 
                ijkc2[0] = icopy/(nc[2]*nc[1]); assert((ijkc[0] >= 0)&&(ijkc[0] < nc[0]));
                ijkc2[bit_pair1] -= 1; // forward bit has it so the copy is the previous in this direction
                ijkc2[bit_pair0] += 1; // backward bit has it so the copy is the next in this direction
                int icopy2 = ((ijkc2[0]*nc[1]+ijkc2[1])*nc[2]+ijkc2[2]);
                assert(candidateMap.find(ino_global) == candidateMap.end());
                candidateMap[ino_global] = icopy2;
              }
            }
            if ((ijkc[bit_pair0] > 0)&&(no_bits[ino]&(1<<(2*bit_pair0+1)))) {
              int bits0 = no_bits[ino];
              bits0 &= ~(1<<(2*bit_pair1+1));
              bits0 |=  (1<<(2*bit_pair1  ));
              bits0 &= ~(1<<(2*bit_pair0+1));
              bits0 |=  (1<<(2*bit_pair0  ));
              map<const int,int8>::iterator iter2 = parentBitMap[iter->second].find(bits0);
              if (iter2 != parentBitMap[iter->second].end()) {
                const int8 ino_global = iter2->second;
                int ijkc2[3];
                ijkc2[2] = icopy%nc[2]; assert((ijkc[2] >= 0)&&(ijkc[2] < nc[2]));
                ijkc2[1] = (icopy/nc[2])%nc[1]; assert((ijkc[1] >= 0)&&(ijkc[1] < nc[1])); 
                ijkc2[0] = icopy/(nc[2]*nc[1]); assert((ijkc[0] >= 0)&&(ijkc[0] < nc[0]));
                ijkc2[bit_pair1] -= 1; // forward bit has it so the copy is the previous in this direction
                ijkc2[bit_pair0] -= 1; // forward bit has it so the copy is the previous in this direction
                int icopy2 = ((ijkc2[0]*nc[1]+ijkc2[1])*nc[2]+ijkc2[2]);
                assert(candidateMap.find(ino_global) == candidateMap.end());
                candidateMap[ino_global] = icopy2;
              }
            }
          }
        }
      }
      // 3d...
      if ((nc[0] > 1)&&(nc[1] > 1)&&(nc[2] > 1)) {
        if ((ijkc[0] > 0)&&(no_bits[ino]&(1<<(2*0+1)))&&
            (ijkc[1] > 0)&&(no_bits[ino]&(1<<(2*1+1)))&&
            (ijkc[2] > 0)&&(no_bits[ino]&(1<<(2*2+1)))) {
          int bits0 = no_bits[ino];
          bits0 &= ~(1<<(2*0+1));
          bits0 |=  (1<<(2*0  ));
          bits0 &= ~(1<<(2*1+1));
          bits0 |=  (1<<(2*1  ));
          bits0 &= ~(1<<(2*2+1));
          bits0 |=  (1<<(2*2  ));
          map<const int,int8>::iterator iter2 = parentBitMap[iter->second].find(bits0);
          if (iter2 != parentBitMap[iter->second].end()) {
            const int8 ino_global = iter2->second;
            int ijkc2[3];
            ijkc2[2] = icopy%nc[2]; assert((ijkc[2] >= 0)&&(ijkc[2] < nc[2]));
            ijkc2[1] = (icopy/nc[2])%nc[1]; assert((ijkc[1] >= 0)&&(ijkc[1] < nc[1])); 
            ijkc2[0] = icopy/(nc[2]*nc[1]); assert((ijkc[0] >= 0)&&(ijkc[0] < nc[0]));
            ijkc2[0] -= 1; // forward bit has it so the copy is the previous in this direction
            ijkc2[1] -= 1; // forward bit has it so the copy is the previous in this direction
            ijkc2[2] -= 1; // forward bit has it so the copy is the previous in this direction
            int icopy2 = ((ijkc2[0]*nc[1]+ijkc2[1])*nc[2]+ijkc2[2]);
            assert(candidateMap.find(ino_global) == candidateMap.end());
            candidateMap[ino_global] = icopy2;
          }
        }
      }
      assert(!candidateMap.empty());
      ino_global_icopy2_vec[ii] = pair<int8,int>(candidateMap.begin()->first,candidateMap.begin()->second);
      /*
        if (candidateMap.size() > 1) {
        for (map<const int8,int>::iterator iter4 = candidateMap.begin(); iter4 != candidateMap.end(); ++iter4) {
        cout << "?: " << iter4->first << " " << iter4->second << endl;
        }

        double x[3]; FOR_I3 x[i] = x_no_base[ino][i];
        FOR_I3 {
        for (int ii = 0; ii < ijkc[i]; ++ii) {
        if (ncopy[i] > 0) 
        PeriodicData::periodicTransformVec[i].translate(x);
        else if (ncopy[i] < 0) 
        PeriodicData::periodicTransformVec[i].inv_translate(x);
        }
        }

        for (int ino_new2 = 0; ino_new2 < nno_new; ++ino_new2) {
        const int ino2     = ino_icopy_vec[ino_new2].first;
        const int icopy2 = ino_icopy_vec[ino_new2].second;
        int ijkc2[3];
        ijkc2[2] = icopy2%nc[2]; assert((ijkc[2] >= 0)&&(ijkc[2] < nc[2]));
        ijkc2[1] = (icopy2/nc[2])%nc[1]; assert((ijkc[1] >= 0)&&(ijkc[1] < nc[1])); 
        ijkc2[0] = icopy2/(nc[2]*nc[1]); assert((ijkc[0] >= 0)&&(ijkc[0] < nc[0]));
        assert(icopy == ((ijkc[0]*nc[1]+ijkc[1])*nc[2]+ijkc[2]));

        if (DIST(x,x_no[ino_new2]) < 1.0E-6) {
        cout << COUT_VEC(x) << " " << ino_parent_global << endl;
        cout << "1: " << ino << " " << icopy << " " << COUT_VEC(ijkc) << " "; BitUtils::dumpBits(no_bits[ino],6);
        cout << "2: " << ino2 << " " << icopy2 << " " << COUT_VEC(ijkc2) << " "; BitUtils::dumpBits(no_bits[ino2],6);
        for (map<const int,int8>::iterator iter3 = parentBitMap[iparent].begin(); iter3 != parentBitMap[iparent].end(); ++iter3) {
        const int bits = iter3->first;
        const int ino_global = iter3->second;
        cout << " > " <<  ino_global << " "; 
        if (ino_global == ino_parent_global) cout << "! ";
        BitUtils::dumpBits(bits,6);
        }
        break;
        }
        }
        }
      */
    }
    delete[] x_no_base;
    for (int iparent = 0; iparent < nparent; ++iparent)
      parentBitMap[iparent].clear();
    delete[] parentBitMap;
    globalMap.clear();


    // now we need to get the nnono_2_i/v for the nodes that live on other ranks. we know that these are in ino_global_icopy2_vec...

    FOR_RANK send_count_3[rank] = 0;
    int* ino_new2_send;
    for (int iter = 0; iter < 2; ++iter) {
      for (int ino_new2 = 0; ino_new2 < nno_new2; ++ino_new2) {
        const int8 ino_global = ino_global_icopy2_vec[ino_new2].first;
        const int rank = MiscUtils::getRankInXora(ino_global,noora);
        const int ino = ino_global-noora[rank];
        if (iter == 0) {
          send_count_3[rank] += 2;
        }
        else {
          const int icopy2 = ino_global_icopy2_vec[ino_new2].second;
          ino_new2_send[send_disp_3[rank]/2] = ino_new2;
          send_buf_int_3[send_disp_3[rank]  ] = ino;
          send_buf_int_3[send_disp_3[rank]+1] = icopy2;
          send_disp_3[rank] += 2;
        }
      }

      // rewind...
      send_disp_3[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp_3[rank] = send_count_3[rank-1] + send_disp_3[rank-1];

      if (iter == 0) {
        send_count_sum_3 = send_disp_3[mpi_size-1] + send_count_3[mpi_size-1];
        assert(send_count_sum_3%2 == 0);
        send_buf_int_3 = new int[send_count_sum_3];
        ino_new2_send = new int[send_count_sum_3/2];
      }
    }
    
    // exchange...
    MPI_Alltoall(send_count_3,1,MPI_INT,recv_count_3,1,MPI_INT,mpi_comm);

    recv_disp_3[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp_3[rank] = recv_count_3[rank-1] + recv_disp_3[rank-1];
    recv_count_sum_3 = recv_disp_3[mpi_size-1] + recv_count_3[mpi_size-1];

    assert(recv_count_sum_3%2 == 0);
    assert(recv_buf_int_3 == NULL); recv_buf_int_3 = new int[recv_count_sum_3];
    MPI_Alltoallv(send_buf_int_3,send_count_3,send_disp_3,MPI_INT,
                  recv_buf_int_3,recv_count_3,recv_disp_3,MPI_INT,mpi_comm);

    FOR_RANK {
      for (int irecv = recv_disp_3[rank]; irecv < recv_disp_3[rank]+recv_count_3[rank]; irecv += 2) {
        const int ino = recv_buf_int_3[irecv  ];
        const int icopy = recv_buf_int_3[irecv+1];
        //cout << ino << " " << nno_p << " " << icopy << endl; 
        assert((ino >= 0)&&(ino < nno_p)); // should be a periodic one...
        // find the associated new node...
        int non;
        for (non = nnono_i[ino]; non != nnono_i[ino+1]; ++non) {
          assert(ino == ino_icopy_vec[non].first);
          //cout << " > " << non << " " << ino_icopy_vec[non].second << endl; 
          if (ino_icopy_vec[non].second == icopy)
            break;
        }
        assert(non != nnono_i[ino+1]); // found copy
        assert(irecv%2 == 0);
        recv_buf_int_3[irecv/2] = non; // send back ino_new... 
      }
    }

    FOR_RANK {
      send_count_3[rank] /= 2;
      send_disp_3[rank] /= 2;
      recv_count_3[rank] /= 2;
      recv_disp_3[rank] /= 2;
    }


    // send back...
    MPI_Alltoallv(recv_buf_int_3,recv_count_3,recv_disp_3,MPI_INT,
                  send_buf_int_3,send_count_3,send_disp_3,MPI_INT,mpi_comm);
    delete[] recv_buf_int_3; recv_buf_int_3 = NULL;

    int8* ino_global_new2 = new int8[nno_new2];
    {
      //int ino2 = 0;
      FOR_RANK {
        for (int isend = send_disp_3[rank]; isend < send_disp_3[rank]+send_count_3[rank]; ++isend) {
          const int ino_new = send_buf_int_3[isend];
          const int ino_global_new = ino_new + noora_new[rank];
          const int ino_new2 = ino_new2_send[isend];
          ino_global_new2[ino_new2] = ino_global_new;
        }
      }
    }
    delete[] send_buf_int_3; send_buf_int_3 = NULL;
    delete[] ino_new2_send;
    
    // use the recv_counts and nnono_count (from when we got no_bits)...

    FOR_RANK send_count[rank] = 0;
    for (int iter = 0; iter < 2; ++iter) {
      for (int rank = 0; rank < mpi_size; ++rank) {
        for (int irecv = recv_disp[rank]; irecv < recv_disp[rank]+recv_count[rank]; irecv += 3) {
          const int ino_tmp = recv_buf_int[irecv];
          const int ino = recv_buf_int[irecv+1];
          //const int bits = recv_buf_int[irecv+2];
          if (iter == 0) {
            // +3 to store counts and ino_tmp, x3 to store ino_new,rank_new,icopy
            send_count[rank] += 3*(nnono_i[ino+1]-nnono_i[ino]+nnono2_i[ino+1]-nnono2_i[ino])+3;
          }
          else {
            send_buf_int[send_disp[rank]++] = -(nnono_i[ino+1]-nnono_i[ino])-1; // first is -1 indexed on rank count
            send_buf_int[send_disp[rank]++] = nnono2_i[ino+1]-nnono2_i[ino]; // off rank
            send_buf_int[send_disp[rank]++] = ino_tmp;
            for (int non = nnono_i[ino]; non != nnono_i[ino+1]; ++non) {
              assert(ino == ino_icopy_vec[non].first);
              assert((non >= 0)&&(non < nno_new));
              send_buf_int[send_disp[rank]++] = mpi_rank;                  // get rank
              send_buf_int[send_disp[rank]++] = non;                       // get the new nodal index
              send_buf_int[send_disp[rank]++] = ino_icopy_vec[non].second; // get the copy index
            }
            for (int non = nnono2_i[ino]; non != nnono2_i[ino+1]; ++non) {
              const int rank_new = MiscUtils::getRankInXora(ino_global_new2[non],noora_new);
              const int ino_new = ino_global_new2[non]-noora_new[rank_new];
              send_buf_int[send_disp[rank]++] = rank_new;                       // get the rank
              send_buf_int[send_disp[rank]++] = ino_new;                        // get the new nodal index
              send_buf_int[send_disp[rank]++] = ino_icopy_skip_vec[non].second; // get the copy index
            }
          }
        }
      }

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      if (iter == 0) {
        send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
        assert(send_buf_int == NULL); send_buf_int = new int[send_count_sum];
      }
    }
    delete[] no_bits;
    delete[] recv_buf_int; recv_buf_int = NULL; 

    //  exchange...
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    
    assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
		  recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; send_buf_int = NULL; 

    // build nnono_tmp along with the associated icopy indices...
    int *nnono_tmp_i = new int[nno_tmp+1];
    for (int ino_tmp = 0; ino_tmp <= nno_tmp; ++ino_tmp) 
      nnono_tmp_i[ino_tmp] = 0;
    int *cp_nn_tmp = NULL;
    int8 *nnono_tmp_v_global = NULL;
    for (int iter = 0; iter < 2; ++iter) {
      int irecv = 0; //,irecv0 = 0;
      while (irecv < recv_count_sum) {
        // new nodal group...
        assert(recv_buf_int[irecv] < 0);
        const int count_on_rank = -recv_buf_int[irecv++]-1;
        const int count_off_rank = recv_buf_int[irecv++];
        const int count = count_on_rank+count_off_rank;
        const int ino_tmp = recv_buf_int[irecv++];
        if (iter == 0) {
          assert(nnono_tmp_i[ino_tmp+1] == 0); // haven't found this one yet...
          nnono_tmp_i[ino_tmp+1] += count;
          irecv += 3*count; // go ahead and skip the numbers
        }
        else if (iter == 1) {
          for (int ii = 0; ii < count; ++ii) {
            const int rank_new = recv_buf_int[irecv++];
            const int ino_new = recv_buf_int[irecv++];
            nnono_tmp_v_global[nnono_tmp_i[ino_tmp]] = noora_new[rank_new]+ino_new;
            cp_nn_tmp[nnono_tmp_i[ino_tmp]] = recv_buf_int[irecv++];
            ++nnono_tmp_i[ino_tmp];
          }
        }
      }

      if (iter == 0) {
        for (int ino_tmp = 0; ino_tmp < nno_tmp; ++ino_tmp) 
          nnono_tmp_i[ino_tmp+1] += nnono_tmp_i[ino_tmp];
        nnono_tmp_v_global = new int8[nnono_tmp_i[nno_tmp]];
        cp_nn_tmp = new int[nnono_tmp_i[nno_tmp]];
      }
      else {
        // rewind 
        for (int ino_tmp = nno_tmp; ino_tmp > 0; --ino_tmp)
          nnono_tmp_i[ino_tmp] = nnono_tmp_i[ino_tmp-1];
        nnono_tmp_i[0] = 0;
      }
    }
    delete[] recv_buf_int; recv_buf_int = NULL;

    // now that we have the new global nodes and their associated icopy, lets populate the new noofa_i/v_global...

    int *noofa_i_base = noofa_i;
    noofa_i = new int[nfa_new+1];
    noofa_i[0] = 0;
    for (int iter = 0; iter < 2; ++iter) {
      for (int ifa_new = 0; ifa_new < nfa_new; ++ifa_new) {
        const int ifa0 = fazone_ifa0_vec[ifa_new].second;
        const int ifa  = ifa_icopy_vec[ifa0].first;
        if (iter == 0) {
          noofa_i[ifa_new+1] = noofa_i[ifa_new]+(noofa_i_base[ifa+1]-noofa_i_base[ifa]);
        }
        else {
          const int icopy = ifa_icopy_vec[ifa0].second;
          for (int nof = noofa_i_base[ifa]; nof != noofa_i_base[ifa+1]; ++nof) {
            const int ino_tmp = noofa_v_tmp[nof];
            // TODO: bissect?, use set instead?
            int8 ino_global_new = -1;
            for (int non = nnono_tmp_i[ino_tmp]; non != nnono_tmp_i[ino_tmp+1]; ++non) {
              if (cp_nn_tmp[non] == icopy) {
                assert(ino_global_new == -1);
                ino_global_new = nnono_tmp_v_global[non];
                break;
              }
            }
            assert(ino_global_new != -1); // found copy
            noofa_v_global[noofa_i[ifa_new]++] = ino_global_new;
          }
        }
      }
      if (iter == 0) {
        assert(noofa_v_global == NULL); noofa_v_global = new int8[noofa_i[nfa_new]];
      }
      else {
        // rewind 
        for (int ifa_new = nfa_new; ifa_new > 0; --ifa_new)
          noofa_i[ifa_new] = noofa_i[ifa_new-1];
        noofa_i[0] = 0;
      }

    }
    delete[] noofa_v_tmp;
    delete[] noofa_i_base;
    delete[] nnono_tmp_i; nnono_tmp_i = NULL;
    delete[] nnono_tmp_v_global; nnono_tmp_v_global = NULL;
    delete[] cp_nn_tmp; cp_nn_tmp = NULL;
    ifa_icopy_vec.clear();

    // because we killed the bits of the periodic faces that become interior faces, we need to do
    // another shuffle to order them s.t. periodic ones come globally first...

    int8 my_nfa_pi_disp[2];
    int8 my_nfa_pi[2] = {nfa_p,nfa_new-nfa_p};
    MPI_Scan(my_nfa_pi,my_nfa_pi_disp,2,MPI_INT8,MPI_SUM,mpi_comm);
    int8 nfa_pi_global[2] = {my_nfa_pi_disp[0],my_nfa_pi_disp[1]}; 
    MPI_Bcast(nfa_pi_global,2,MPI_INT8,mpi_size-1,mpi_comm);
    FOR_I2 my_nfa_pi_disp[i] -= my_nfa_pi[i];
    
    // get ifa_global_new...

    int8 *ifa_global_new = new int8[nfa_new];
    for (int ifa_new = 0; ifa_new < nfa_p; ++ifa_new) {
      const int bits = getBitsForFaceZone(fazone_ifa0_vec[ifa_new].first);
      assert(bits);
      ifa_global_new[ifa_new] = my_nfa_pi_disp[0]+ifa_new;
    }
    for (int ifa_new = nfa_p; ifa_new < nfa_new; ++ifa_new) {
      const int bits = getBitsForFaceZone(fazone_ifa0_vec[ifa_new].first);
      assert(bits == 0);
      ifa_global_new[ifa_new] = nfa_pi_global[0]+my_nfa_pi_disp[1]+(ifa_new-nfa_p);
    }
    fazone_ifa0_vec.clear();

    delete[] faora; faora = NULL;
    const int8 nfa_global0 = nfa_global;
    nfa_p_global = nfa_pi_global[0];
    nfa_global = nfa_pi_global[0]+nfa_pi_global[1];
    MiscUtils::calcUniformDist(faora,nfa_global,mpi_size);
    assert(faora[mpi_rank+1]-faora[mpi_rank] < TWO_BILLION);
    assert(faora[mpi_size] == nfa_global);
    nfa = faora[mpi_rank+1]-faora[mpi_rank];

    // noofa_i/v_global...

    FOR_RANK {
      send_count[rank] = 0;
      send_count_3[rank] = 0;
    }
    int8 * send_buf_int8;
    double * send_buf_double;
    for (int iter = 0; iter < 2; ++iter) {
      for (int ifa_new = 0; ifa_new < nfa_new; ++ifa_new) {
        const int rank = MiscUtils::getRankInXora(ifa_global_new[ifa_new],faora);
        const int nnof = noofa_i[ifa_new+1]-noofa_i[ifa_new];
        if (iter == 0) {
          send_count[rank] += 4 + nnof; // ifa, cvofa[2], nnof and nodes
          send_count_3[rank] += 6; // x_fa, n_fa
        }
        else {
          send_buf_int8[send_disp[rank]++] = ifa_global_new[ifa_new]-faora[rank];
          FOR_I2 send_buf_int8[send_disp[rank]++] = cvofa_global[ifa_new][i];
          send_buf_int8[send_disp[rank]++] = nnof;
          for (int nof = noofa_i[ifa_new]; nof != noofa_i[ifa_new+1]; ++nof)
            send_buf_int8[send_disp[rank]++] = noofa_v_global[nof];
          FOR_I3 send_buf_double[send_disp_3[rank]++] = x_fa[ifa_new][i];
          FOR_I3 send_buf_double[send_disp_3[rank]++] = n_fa[ifa_new][i];
        }
      }

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      send_disp_3[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp_3[rank] = send_count_3[rank-1] + send_disp_3[rank-1];

      if (iter == 0) {
        send_buf_int8 = new int8[send_disp[mpi_size-1] + send_count[mpi_size-1]];
        send_buf_double = new double[send_disp_3[mpi_size-1] + send_count_3[mpi_size-1]];
      }
    }
    delete[] ifa_global_new;
    delete[] x_fa; x_fa = new double[nfa][3];
    delete[] n_fa; n_fa = new double[nfa][3];
    delete[] cvofa_global; cvofa_global = new int8[nfa][2];
    delete[] noofa_i; noofa_i = new int[nfa+1]; noofa_i[0] = 0;
    delete[] noofa_v_global; // sized below...
	
    //  exchange...
    
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    
    int8* recv_buf_int8 = new int8[recv_count_sum];
    MPI_Alltoallv(send_buf_int8,send_count,send_disp,MPI_INT8,
		  recv_buf_int8,recv_count,recv_disp,MPI_INT8,mpi_comm);
    delete[] send_buf_int8; send_buf_int8 = NULL; 

    MPI_Alltoall(send_count_3,1,MPI_INT,recv_count_3,1,MPI_INT,mpi_comm);

    recv_disp_3[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp_3[rank] = recv_count_3[rank-1] + recv_disp_3[rank-1];
    recv_count_sum_3 = recv_disp_3[mpi_size-1] + recv_count_3[mpi_size-1];

    assert(recv_count_sum_3%6 == 0);
    double * recv_buf_double = new double[recv_count_sum_3];
    MPI_Alltoallv(send_buf_double,send_count_3,send_disp_3,MPI_DOUBLE,
                  recv_buf_double,recv_count_3,recv_disp_3,MPI_DOUBLE,mpi_comm);
    delete[] send_buf_double; send_buf_double = NULL;

    for (int iter = 0; iter < 2; ++iter) {
      int irecv = 0;
      int irecv2 = 0;
      while (irecv < recv_count_sum) {
        const int ifa = recv_buf_int8[irecv++];
        if (iter == 0) {
          FOR_I2 cvofa_global[ifa][i] = recv_buf_int8[irecv++];
          FOR_I3 x_fa[ifa][i] = recv_buf_double[irecv2++];
          FOR_I3 n_fa[ifa][i] = recv_buf_double[irecv2++];
          noofa_i[ifa+1] = recv_buf_int8[irecv++];
          irecv += noofa_i[ifa+1];

          /*
            const int icv1_global = (cvofa_global[ifa][1]&MASK_55BITS);
            const int icv0_global = cvofa_global[ifa][0];
            if (icv1_global < icv0_global) {
            cout << ifa << " " << icv1_global << " " << icv0_global << endl;
            const int bits = (cvofa_global[ifa][1]>>55);
            cvofa_global[ifa][0] = icv1_global;
            cvofa_global[ifa][1] = (icv0_global|(uint8(BitUtils::flipPeriodicBits(bits))<<55));;
            FOR_I3 n_fa[ifa][i] *= -1.0;
            }
          */
        }
        else {
          irecv += 3; 
          for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) 
            noofa_v_global[nof] = recv_buf_int8[irecv++];
        }
      }

      if (iter == 0) {
        delete[] recv_buf_double; recv_buf_double = NULL;
        for (int ifa = 0; ifa < nfa; ++ifa) 
          noofa_i[ifa+1] += noofa_i[ifa];
        noofa_v_global = new int8[noofa_i[nfa]];
      }
    }
    delete[] recv_buf_int8; recv_buf_int8 = NULL;

    delete[] send_disp_3; send_disp_3 = NULL;
    delete[] send_count_3; send_count_3 = NULL;
    delete[] recv_disp_3; recv_disp_3 = NULL;
    delete[] recv_count_3; recv_count_3 = NULL;

    /*
      double max_dist = 0.0;
      for (int ifa_new = 0; ifa_new < nfa_new; ++ifa_new) {
      cout << ifa_new << " " << COUT_VEC(x_fa[ifa_new]) << " ";  
      double x_avg[3] = {0.0,0.0,0.0};
      double area = 0.0;
      int nof0 = noofa_i[ifa_new];
      int nof = noofa_i[ifa_new]+1;
      for (int nof1 = noofa_i[ifa_new]+2; nof1 < noofa_i[ifa_new+1]; ++nof1) {
      const double *x0 = x_no[noofa_v_global[nof0]];
      const double *x  = x_no[noofa_v_global[nof]];
      const double *x1 = x_no[noofa_v_global[nof1]];
      const double n[3] = TRI_NORMAL_2(x0,x,x1);
      const double this_area = MAG(n);
      FOR_I3 x_avg[i] += this_area*(x0[i]+x[i]+x1[i]);
      area += this_area;
      nof++;
      }
      FOR_I3 x_avg[i] /= 3.0*area;
      max_dist = max(max_dist,DIST(x_avg,x_fa[ifa_new]));
      cout << COUT_VEC(x_avg) << " " << DIST(x_avg,x_fa[ifa_new]) << endl;
      }
      cout << "MAX: " << max_dist << endl;
    */
    
    // ----------
    // bf data...
    // ----------

    const int nbf_new = nbf*ncopies;
    vector<pair<int,int> > bfzone_ibf0_vec; bfzone_ibf0_vec.reserve(nbf_new);
    vector<pair<int,int> > ibf_icopy_vec; ibf_icopy_vec.reserve(nbf_new);
    {
      int icopy = 0;
      int ibf0 = 0;
      for (int ic = 0; ic < nc[0]; ++ic) {
        for (int jc = 0; jc < nc[1]; ++jc) {
          for (int kc = 0; kc < nc[2]; ++kc) {
            assert(((ic*nc[1]+jc)*nc[2]+kc) == icopy);
            FOR_IBF {
              bfzone_ibf0_vec.push_back(pair<int,int>(zone_bf[ibf],ibf0++));
              ibf_icopy_vec.push_back(pair<int,int>(ibf,icopy));
            }
            ++icopy;
          }
        }
      }
      assert(ibf0 == nbf_new);
      assert(icopy == ncopies);
    }
    sort(bfzone_ibf0_vec.begin(),bfzone_ibf0_vec.end());

    // area_bf,zone_bf,x_bf,n_bf,Gij_bf,cvobf_global...

    assert(area_bf != NULL);
    double *area_bf_base = area_bf;
    area_bf = new double[nbf_new]; 
    assert(area_over_delta_bf != NULL);
    double *area_over_delta_bf_base = area_over_delta_bf;
    area_over_delta_bf = new double[nbf_new]; 
    assert(zone_bf != NULL);
    int *zone_bf_base = zone_bf;
    zone_bf = new int[nbf_new]; 
    assert(x_bf != NULL);
    double (*x_bf_base)[3] = x_bf;
    x_bf = new double[nbf_new][3]; 
    assert(n_bf != NULL);
    double (*n_bf_base)[3] = n_bf;
    n_bf = new double[nbf_new][3]; 
    assert(Gij_bf != NULL);
    double (*Gij_bf_base)[3][3] = Gij_bf;
    Gij_bf = new double[nbf_new][3][3]; 
    assert(cvobf_global != NULL);
    int8 *cvobf_global_base = cvobf_global;
    cvobf_global = new int8[nbf_new];
    for (int ibf_new = 0; ibf_new < nbf_new; ++ibf_new) {
      const int izone = bfzone_ibf0_vec[ibf_new].first;
      const int ibf0  = bfzone_ibf0_vec[ibf_new].second;
      const int ibf   = ibf_icopy_vec[ibf0].first;
      const int icopy = ibf_icopy_vec[ibf0].second;
      int ijkc[3];
      ijkc[2] = icopy%nc[2]; assert((ijkc[2] >= 0)&&(ijkc[2] < nc[2]));
      ijkc[1] = (icopy/nc[2])%nc[1]; assert((ijkc[1] >= 0)&&(ijkc[1] < nc[1])); 
      ijkc[0] = icopy/(nc[2]*nc[1]); assert((ijkc[0] >= 0)&&(ijkc[0] < nc[0]));
      assert(icopy == ((ijkc[0]*nc[1]+ijkc[1])*nc[2]+ijkc[2]));
      // zone_bf...
      assert(izone == zone_bf_base[ibf]);
      zone_bf[ibf_new] = izone;
      // area_bf...
      area_bf[ibf_new] = area_bf_base[ibf];
      // area_over_delta_bf...
      area_over_delta_bf[ibf_new] = area_over_delta_bf_base[ibf];
      // x_bf...
      FOR_I3 x_bf[ibf_new][i] = x_bf_base[ibf][i];
      FOR_I3 {
        for (int ii = 0; ii < ijkc[i]; ++ii) {
          if (ncopy[i] > 0) 
            PeriodicData::periodicTransformVec[i].translate(x_bf[ibf_new]);
          else if (ncopy[i] < 0) 
            PeriodicData::periodicTransformVec[i].inv_translate(x_bf[ibf_new]);
        }
      }
      // n_bf...
      FOR_I3 n_bf[ibf_new][i] = n_bf_base[ibf][i];
      FOR_I3 {
        for (int ii = 0; ii < ijkc[i]; ++ii) {
          if (ncopy[i] > 0) 
            PeriodicData::periodicTransformVec[i].rotate(n_bf[ibf_new]);
          else if (ncopy[i] < 0) 
            PeriodicData::periodicTransformVec[i].inv_rotate(n_bf[ibf_new]);
        }
      }
      // Gij_bf...
      FOR_I3 FOR_J3 Gij_bf[ibf_new][i][j] = Gij_bf_base[ibf][i][j];
      FOR_I3 {
        for (int ii = 0; ii < ijkc[i]; ++ii) {
          if (ncopy[i] > 0) 
            PeriodicData::periodicTransformVec[i].rotate(Gij_bf[ibf_new]);
          else if (ncopy[i] < 0) 
            PeriodicData::periodicTransformVec[i].inv_rotate(Gij_bf[ibf_new]);
        }
      }
      // cvobf_global...
      //cvobf_global[ibf_new] = int8(icopy)*ncv_global0+cvobf_global_base[ibf];
      const int rank = MiscUtils::getRankInXora(cvobf_global_base[ibf],cvora0);          
      cvobf_global[ibf_new] = int8(icopy)*(cvora0[rank+1]-cvora0[rank]) + cvobf_global_base[ibf] - cvora0[rank] + cvora[rank];
    }
    delete[] cvora0;
    delete[] zone_bf_base;
    delete[] area_bf_base;
    delete[] area_over_delta_bf_base;
    delete[] x_bf_base;
    delete[] n_bf_base;
    delete[] Gij_bf_base;
    delete[] cvobf_global_base;

    if (b_surface) {
      // surface_sbobf_i/v_global...
      assert(surface_sbobf_i != NULL);
      int *surface_sbobf_i_base = surface_sbobf_i;
      surface_sbobf_i = new int[nbf_new+1];
      surface_sbobf_i[0] = 0;
      assert(surface_sbobf_v_global != NULL);
      uint8 *surface_sbobf_v_global_base = surface_sbobf_v_global;
      surface_sbobf_v_global = new uint8[surface_sbobf_i_base[nbf]*ncopies];
      // surface_spobf_i/v/wgt_global...
      assert(surface_spobf_i != NULL);
      int *surface_spobf_i_base = surface_spobf_i;
      surface_spobf_i = new int[nbf_new+1];
      surface_spobf_i[0] = 0;
      assert(surface_spobf_v_global != NULL);
      uint8 *surface_spobf_v_global_base = surface_spobf_v_global;
      surface_spobf_v_global = new uint8[surface_spobf_i_base[nbf]*ncopies];
      assert(surface_spobf_wgt_global != NULL);
      double *surface_spobf_wgt_global_base = surface_spobf_wgt_global;
      surface_spobf_wgt_global = new double[surface_spobf_i_base[nbf]*ncopies];
      for (int ibf_new = 0; ibf_new < nbf_new; ++ibf_new) {
        //const int izone = bfzone_ibf0_vec[ibf_new].first;
        const int ibf0  = bfzone_ibf0_vec[ibf_new].second;
        const int ibf   = ibf_icopy_vec[ibf0].first;
        const int icopy = ibf_icopy_vec[ibf0].second;
        int ijkc[3];
        ijkc[2] = icopy%nc[2]; assert((ijkc[2] >= 0)&&(ijkc[2] < nc[2]));
        ijkc[1] = (icopy/nc[2])%nc[1]; assert((ijkc[1] >= 0)&&(ijkc[1] < nc[1])); 
        ijkc[0] = icopy/(nc[2]*nc[1]); assert((ijkc[0] >= 0)&&(ijkc[0] < nc[0]));
        assert(icopy == ((ijkc[0]*nc[1]+ijkc[1])*nc[2]+ijkc[2]));
        // sbobf_i/v_global...
        surface_sbobf_i[ibf_new+1] = surface_sbobf_i[ibf_new] + (surface_sbobf_i_base[ibf+1]-surface_sbobf_i_base[ibf]);
        for (int stob = surface_sbobf_i_base[ibf]; stob < surface_sbobf_i_base[ibf+1]; ++stob) {
          const int stob_new = surface_sbobf_i[ibf_new] + (stob-surface_sbobf_i_base[ibf]);
          surface_sbobf_v_global[stob_new] = surface_sbobf_v_global_base[stob];
        }
        // spobf_i/v/wgt_global...
        surface_spobf_i[ibf_new+1] = surface_spobf_i[ibf_new] + (surface_spobf_i_base[ibf+1]-surface_spobf_i_base[ibf]);
        for (int spob = surface_spobf_i_base[ibf]; spob < surface_spobf_i_base[ibf+1]; ++spob) {
          const int spob_new = surface_spobf_i[ibf_new] + (spob-surface_spobf_i_base[ibf]);
          surface_spobf_v_global[spob_new] = surface_spobf_v_global_base[spob];
          surface_spobf_wgt_global[spob_new] = surface_spobf_wgt_global_base[spob];
        }
        
      }
      delete[] surface_sbobf_i_base;
      delete[] surface_sbobf_v_global_base;
      delete[] surface_spobf_i_base;
      delete[] surface_spobf_v_global_base;
      delete[] surface_spobf_wgt_global_base;
    }

    // noobf_i/v...

    // get nodes from boundary faces...
    
    assert(globalMap.empty());
    nno_tmp = 0;
    FOR_IBF {
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int8 ino_global = noobf_v_global[nob];
        assert((ino_global >= 0)&&(ino_global < nno_global));
        map<const int8,int>::iterator iter = globalMap.find(ino_global);
        if (iter == globalMap.end()) 
          globalMap[ino_global] = nno_tmp++;
      }
    }
    
    assert(ino_global_tmp == NULL); ino_global_tmp = new int8[nno_tmp];
    int *noobf_v_tmp = new int[noobf_i[nbf]];
    FOR_IBF {
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int8 ino_global = noobf_v_global[nob];
        assert((ino_global >= 0)&&(ino_global < nno_global));
        map<const int8,int>::iterator iter = globalMap.find(ino_global);
        assert(iter != globalMap.end());
        noobf_v_tmp[nob] = iter->second;
        ino_global_tmp[iter->second] = ino_global;
      }
    }
    delete[] noobf_v_global; noobf_v_global = NULL; // replaced with noobf_v_tmp 
    globalMap.clear();

    FOR_RANK send_count[rank] = 0;
    for (int ino_tmp = 0; ino_tmp < nno_tmp; ++ino_tmp) {
      const int rank = MiscUtils::getRankInXora(ino_global_tmp[ino_tmp],noora);
      send_count[rank] += 2; // ino_temp, ino (striped)
    }

    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
    
    assert(send_buf_int == NULL); send_buf_int = new int[send_count_sum];
    for (int ino_tmp = 0; ino_tmp < nno_tmp; ++ino_tmp) {
      const int rank = MiscUtils::getRankInXora(ino_global_tmp[ino_tmp],noora);
      send_buf_int[send_disp[rank]++] = ino_tmp;
      send_buf_int[send_disp[rank]++] = ino_global_tmp[ino_tmp]-noora[rank];
    }
    delete[] ino_global_tmp; ino_global_tmp = NULL; 

    //  rewind
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

    //  exchange...
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    
    assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
		  recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; send_buf_int = NULL;

    FOR_RANK send_count[rank] = 0;
    for (int iter = 0; iter < 2; ++iter) {
      for (int rank = 0; rank < mpi_size; ++rank) {
        for (int irecv = recv_disp[rank]; irecv < recv_disp[rank]+recv_count[rank]; irecv += 2) {
          const int ino_tmp = recv_buf_int[irecv];
          const int ino = recv_buf_int[irecv+1];
          if (iter == 0) {
            // +3 to store counts and ino_tmp, x3 to store ino_new,rank_new,icopy
            send_count[rank] += 3*(nnono_i[ino+1]-nnono_i[ino]+nnono2_i[ino+1]-nnono2_i[ino])+3;
          }
          else {
            send_buf_int[send_disp[rank]++] = -(nnono_i[ino+1]-nnono_i[ino])-1; // first is -1 indexed on rank count
            send_buf_int[send_disp[rank]++] = nnono2_i[ino+1]-nnono2_i[ino]; // off rank
            send_buf_int[send_disp[rank]++] = ino_tmp;
            for (int non = nnono_i[ino]; non != nnono_i[ino+1]; ++non) {
              assert(ino == ino_icopy_vec[non].first);
              assert((non >= 0)&&(non < nno_new));
              send_buf_int[send_disp[rank]++] = mpi_rank;                  // get rank
              send_buf_int[send_disp[rank]++] = non;                       // get the new nodal index
              send_buf_int[send_disp[rank]++] = ino_icopy_vec[non].second; // get the copy index
            }
            for (int non = nnono2_i[ino]; non != nnono2_i[ino+1]; ++non) {
              const int rank_new = MiscUtils::getRankInXora(ino_global_new2[non],noora_new);
              const int ino_new = ino_global_new2[non]-noora_new[rank_new];
              send_buf_int[send_disp[rank]++] = rank_new;                       // get the rank
              send_buf_int[send_disp[rank]++] = ino_new;                        // get the new nodal index
              send_buf_int[send_disp[rank]++] = ino_icopy_skip_vec[non].second; // get the copy index
            }
          }
        }
      }

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      if (iter == 0) {
        send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
        assert(send_buf_int == NULL); send_buf_int = new int[send_count_sum];
      }
    }
    ino_icopy_vec.clear();
    ino_icopy_skip_vec.clear();
    delete[] nnono_i;
    delete[] nnono2_i;
    delete[] ino_global_new2;
    delete[] recv_buf_int; recv_buf_int = NULL;

    //  exchange...
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    
    assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
		  recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; send_buf_int = NULL; 
    delete[] send_count;
    delete[] send_disp;

    // build nnono_tmp along with the associated icopy indices...
    assert(nnono_tmp_i == NULL); nnono_tmp_i = new int[nno_tmp+1];
    for (int ino_tmp = 0; ino_tmp <= nno_tmp; ++ino_tmp) 
      nnono_tmp_i[ino_tmp] = 0;
    assert(cp_nn_tmp == NULL); 
    assert(nnono_tmp_v_global == NULL); 
    for (int iter = 0; iter < 2; ++iter) {
      int irecv = 0; //,irecv0 = 0;
      while (irecv < recv_count_sum) {
        // new nodal group...
        assert(recv_buf_int[irecv] < 0);
        const int count_on_rank = -recv_buf_int[irecv++]-1;
        const int count_off_rank = recv_buf_int[irecv++];
        const int count = count_on_rank+count_off_rank;
        const int ino_tmp = recv_buf_int[irecv++];
        if (iter == 0) {
          assert(nnono_tmp_i[ino_tmp+1] == 0); // haven't found this one yet...
          nnono_tmp_i[ino_tmp+1] += count;
          irecv += 3*count; // go ahead and skip the numbers
        }
        else if (iter == 1) {
          for (int ii = 0; ii < count; ++ii) {
            const int rank_new = recv_buf_int[irecv++];
            const int ino_new = recv_buf_int[irecv++];
            nnono_tmp_v_global[nnono_tmp_i[ino_tmp]] = noora_new[rank_new]+ino_new;
            cp_nn_tmp[nnono_tmp_i[ino_tmp]] = recv_buf_int[irecv++];
            ++nnono_tmp_i[ino_tmp];
          }
        }
      }

      if (iter == 0) {
        for (int ino_tmp = 0; ino_tmp < nno_tmp; ++ino_tmp) 
          nnono_tmp_i[ino_tmp+1] += nnono_tmp_i[ino_tmp];
        nnono_tmp_v_global = new int8[nnono_tmp_i[nno_tmp]];
        cp_nn_tmp = new int[nnono_tmp_i[nno_tmp]];
      }
      else {
        // rewind 
        for (int ino_tmp = nno_tmp; ino_tmp > 0; --ino_tmp)
          nnono_tmp_i[ino_tmp] = nnono_tmp_i[ino_tmp-1];
        nnono_tmp_i[0] = 0;
      }
    }
    delete[] recv_buf_int; recv_buf_int = NULL;
    delete[] recv_count;
    delete[] recv_disp;

    // now that we have the new global nodes and their associated icopy, lets populate the new noobf_i/v_global...

    int *noobf_i_base = noobf_i;
    noobf_i = new int[nbf_new+1];
    noobf_i[0] = 0;
    for (int iter = 0; iter < 2; ++iter) {
      for (int ibf_new = 0; ibf_new < nbf_new; ++ibf_new) {
        const int ibf0  = bfzone_ibf0_vec[ibf_new].second;
        const int ibf   = ibf_icopy_vec[ibf0].first;
        const int icopy = ibf_icopy_vec[ibf0].second;
        int ijkc[3];
        ijkc[2] = icopy%nc[2]; assert((ijkc[2] >= 0)&&(ijkc[2] < nc[2]));
        ijkc[1] = (icopy/nc[2])%nc[1]; assert((ijkc[1] >= 0)&&(ijkc[1] < nc[1])); 
        ijkc[0] = icopy/(nc[2]*nc[1]); assert((ijkc[0] >= 0)&&(ijkc[0] < nc[0]));
        assert(icopy == ((ijkc[0]*nc[1]+ijkc[1])*nc[2]+ijkc[2]));
        if (iter == 0) {
          noobf_i[ibf_new+1] = noobf_i[ibf_new]+(noobf_i_base[ibf+1]-noobf_i_base[ibf]);
        }
        else {
          for (int nob = noobf_i_base[ibf]; nob != noobf_i_base[ibf+1]; ++nob) {
            const int ino_tmp = noobf_v_tmp[nob];
            // TODO: bissect?, use set instead?
            int8 ino_global_new = -1;
            for (int non = nnono_tmp_i[ino_tmp]; non != nnono_tmp_i[ino_tmp+1]; ++non) {
              if (cp_nn_tmp[non] == icopy) {
                assert(ino_global_new == -1);
                ino_global_new = nnono_tmp_v_global[non];
                break;
              }
            }
            assert(ino_global_new != -1); // found copy
            noobf_v_global[noobf_i[ibf_new]++] = ino_global_new;
          }
        }
      }
      if (iter == 0) {
        assert(noobf_v_global == NULL); noobf_v_global = new int8[noobf_i[nbf_new]];
      }
      else {
        // rewind 
        for (int ibf_new = nbf_new; ibf_new > 0; --ibf_new)
          noobf_i[ibf_new] = noobf_i[ibf_new-1];
        noobf_i[0] = 0;
      }

    }
    delete[] noobf_v_tmp;
    delete[] noobf_i_base;
    delete[] nnono_tmp_i;
    delete[] nnono_tmp_v_global;
    delete[] cp_nn_tmp;
    bfzone_ibf0_vec.clear();
    ibf_icopy_vec.clear();

    /*
      max_dist = 0.0;
      for (int ibf_new = 0; ibf_new < nbf_new; ++ibf_new) {
      cout << ibf_new << " " << COUT_VEC(x_bf[ibf_new]) << " ";  
      double x_avg[3] = {0.0,0.0,0.0};
      double area = 0.0;
      int nob0 = noobf_i[ibf_new];
      int nob = noobf_i[ibf_new]+1;
      for (int nob1 = noobf_i[ibf_new]+2; nob1 < noobf_i[ibf_new+1]; ++nob1) {
      const double *x0 = x_no[noobf_v_global[nob0]];
      const double *x  = x_no[noobf_v_global[nob]];
      const double *x1 = x_no[noobf_v_global[nob1]];
      const double n[3] = TRI_NORMAL_2(x0,x,x1);
      const double this_area = MAG(n);
      FOR_I3 x_avg[i] += this_area*(x0[i]+x[i]+x1[i]);
      area += this_area;
      nob++;
      }
      FOR_I3 x_avg[i] /= 3.0*area;
      max_dist = max(max_dist,DIST(x_avg,x_bf[ibf_new]));
      cout << COUT_VEC(x_avg) << " " << DIST(x_avg,x_bf[ibf_new]) << endl;
      }
      cout << "MAX BF: " << max_dist << endl;
    */

    // ------------------
    // BfZoneVec stuff...
    // ------------------
    
    {
      int8 ibf_f_global = 0;
      int8 nob_f_global = 0;
      int8 stob_f_global = 0;
      int8 spob_f_global = 0;
      FOR_IZONE(bfZoneVec) {
        // geometry...
        double x_global[3] = {0.0,0.0,0.0};
        double n_global[3] = {0.0,0.0,0.0};
        for (int icopy = 0; icopy < ncopies; ++icopy) {
          int ijkc[3];
          ijkc[2] = icopy%nc[2]; assert((ijkc[2] >= 0)&&(ijkc[2] < nc[2]));
          ijkc[1] = (icopy/nc[2])%nc[1]; assert((ijkc[1] >= 0)&&(ijkc[1] < nc[1])); 
          ijkc[0] = icopy/(nc[2]*nc[1]); assert((ijkc[0] >= 0)&&(ijkc[0] < nc[0]));
          assert(icopy == ((ijkc[0]*nc[1]+ijkc[1])*nc[2]+ijkc[2]));
          double x[3]; FOR_I3 x[i] = bfZoneVec[izone].x_global[i]; 
          FOR_I3 {
            for (int ii = 0; ii < ijkc[i]; ++ii) {
              if (ncopy[i] > 0) 
                PeriodicData::periodicTransformVec[i].translate(x);
              else if (ncopy[i] < 0) 
                PeriodicData::periodicTransformVec[i].inv_translate(x);
            }
          }
          FOR_I3 x_global[i] += bfZoneVec[izone].area_global*x[i];
          double n[3]; FOR_I3 n[i] = bfZoneVec[izone].n_global[i];
          FOR_I3 {
            for (int ii = 0; ii < ijkc[i]; ++ii) {
              if (ncopy[i] > 0) 
                PeriodicData::periodicTransformVec[i].rotate(n);
              else if (ncopy[i] < 0) 
                PeriodicData::periodicTransformVec[i].inv_rotate(n);
            }
          }
          FOR_I3 n_global[i] += n[i];
        }
        bfZoneVec[izone].area_global *= ncopies;     
        FOR_I3 bfZoneVec[izone].x_global[i] /= bfZoneVec[izone].area_global;
        bfZoneVec[izone].area_over_delta_global *= ncopies;     
        // counts...
        bfZoneVec[izone].nbf *= ncopies;
        bfZoneVec[izone].nbf_global *= ncopies;
        bfZoneVec[izone].ibf_f_global = ibf_f_global; ibf_f_global += bfZoneVec[izone].nbf_global;
        delete[] bfZoneVec[izone].bfora; bfZoneVec[izone].bfora = NULL;
        MiscUtils::buildXora(bfZoneVec[izone].bfora,bfZoneVec[izone].nbf); 
        assert(bfZoneVec[izone].bfora[mpi_rank+1]-bfZoneVec[izone].bfora[mpi_rank] == bfZoneVec[izone].nbf);
        assert(bfZoneVec[izone].bfora[mpi_size] == bfZoneVec[izone].nbf_global);
        bfZoneVec[izone].noobf_global *= ncopies;
        bfZoneVec[izone].nob_f_global = nob_f_global; nob_f_global += bfZoneVec[izone].noobf_global;
        if (b_surface) {
          // st...
          bfZoneVec[izone].sbobf_global *= ncopies;
          bfZoneVec[izone].stob_f_global = stob_f_global; stob_f_global += bfZoneVec[izone].sbobf_global;
          // sp...
          bfZoneVec[izone].spobf_global *= ncopies;
          bfZoneVec[izone].spob_f_global = spob_f_global; spob_f_global += bfZoneVec[izone].spobf_global;
        }
      }
    }

    nbf = nbf_new;
    const int8 nbf_global0 = nbf_global; 
    nbf_global *= ncopies;
    delete[] bfora; bfora = NULL;
    MiscUtils::buildXora(bfora,nbf); 
    assert(bfora[mpi_rank+1]-bfora[mpi_rank] < TWO_BILLION);
    assert(nbf == bfora[mpi_rank+1]-bfora[mpi_rank]);
    assert(nbf_global == bfora[mpi_size]);
	
    delete[] noora; noora = noora_new; 
    nno = nno_new;
    const int8 nno_global0 = nno_global;
    nno_global = nno_new_global;

    if (mpi_rank == 0) {
      cout << "Summary: " << endl;
      cout << " > ncv_global: " << ncv_global0 << " -> " << ncv_global << endl;
      cout << " > nfa_global: " << nfa_global0 << " -> " << nfa_global << endl;
      cout << " > nbf_global: " << nbf_global0 << " -> " << nbf_global << endl;
      cout << " > nno_global: " << nno_global0 << " -> " << nno_global << endl;
    }

    // hacking these because they are unused. if needed, we can reorder the nodes
    // to respect pb coming first. 
    nno_p_global = -1; 
    nno_pb_global = -1; 
    for (int i = 0; i < 27; ++i) {
      nfa_zone[i] = -1;
      noofa_zone[i] = -1;
    }

    // update the transforms!
    for (int ipt = 0, npt = PeriodicData::periodicTransformVec.size(); ipt < npt; ++ipt) {
      double data[3]; PeriodicData::periodicTransformVec[ipt].getData(data);
      const int kind = PeriodicData::periodicTransformVec[ipt].getKind(); 
      if (kind == PERIODIC_TRANSFORM_CART) {
        FOR_I3 data[i] *= nc[ipt];
      }
      else if ( (kind == PERIODIC_TRANSFORM_CYL_X) ||
                (kind == PERIODIC_TRANSFORM_CYL_Y) ||
                (kind == PERIODIC_TRANSFORM_CYL_Z) ) {
        data[2] *= nc[ipt];
        data[0] = cos(data[2]*M_PI/180.0);
        data[1] = sin(data[2]*M_PI/180.0);
      }
      else {
        assert(0);
      }
      PeriodicData::periodicTransformVec[ipt].setData(data);
      if (mpi_rank == 0) {
        cout << " > bit pair " << ipt << " "; 
        PeriodicData::periodicTransformVec[ipt].dump();
      }
    }

  }

  class GroupData {
  public:
    double unit_n[3];
    double dp_tol;
    double area_min;
    double area_max;
    GroupData(const double n[3],const double area) {
      FOR_I3 unit_n[i] = n[i];
      dp_tol = 0.0;
      area_min = area;
      area_max = area;
    }
    
    void add(const double n[3],const double area) {
      const double dp = DOT_PRODUCT(unit_n,n);
      dp_tol = max(dp_tol,1.0-dp);
      area_min = min(area_min,area);
      area_max = max(area_max,area);
    }
  };
  
  void check() {

    if (mpi_size == 1) {

      if (b_surface) {

        /*
          for (int ist = 0; ist < surface_nst; ++ist) {
          cout << ist << ": " << COUT_VEC(surface_spost[ist]) << " " << surface_znost[ist] << endl;
          }
          for (int isp = 0; isp < surface_nsp; ++isp) {
          cout << isp << ": " << COUT_VEC(surface_xp[isp]) << endl;
          }
        */

        // we use a circuitous route on purpose here to visit at elements of sbobf as well as the surface
	assert(surface_nst == surface_nst_global);
	assert(surface_nsp == surface_nsp_global);
        double gcl[3] = {0.0,0.0,0.0};
        int * st_flag = new int[surface_nst];
        for (int ist = 0; ist < surface_nst; ++ist) st_flag[ist] = 0;
        for (int ibf = 0; ibf > nbf; ++ibf) {
          for (int stob = surface_sbobf_i[ibf]; stob != surface_sbobf_i[ibf+1]; ++stob) {
            const int ist = int(surface_sbobf_v_global[stob]&MASK_52BITS);
            assert(zone_bf[ibf] == surface_znost[ist]); // should equate...
            if (st_flag[ist] == 0) {
              st_flag[ist] = 1;
              const double n_st[3] = TRI_NORMAL_2(surface_xp[surface_spost[ist][0]],surface_xp[surface_spost[ist][1]],surface_xp[surface_spost[ist][2]]);
              FOR_I3 gcl[i] += n_st[i]; // TODO what to do with bits here?
            }
          }
        }
        cout << " > GCL 1: surface " << COUT_VEC(gcl) << endl;
        FOR_I3 gcl[i] = 0.0;
        for (int ist = 0; ist < surface_nst; ++ist) st_flag[ist] = 0;
        FOR_IZONE(bfZoneVec) {
          for (int stob = bfZoneVec[izone].stob_f_global; stob < bfZoneVec[izone].stob_f_global+bfZoneVec[izone].sbobf_global; ++stob) {
            const int ist = (surface_sbobf_v_global[stob]&MASK_52BITS)-surface_stora[mpi_rank];
            assert((ist >= 0)&&(ist < surface_nst));
            if (st_flag[ist] == 0) {
              st_flag[ist] = 1;
              const double n_st[3] = TRI_NORMAL_2(surface_xp[surface_spost[ist][0]],surface_xp[surface_spost[ist][1]],surface_xp[surface_spost[ist][2]]);
              FOR_I3 gcl[i] += n_st[i]; // TODO also what to do with bits here?
            }
          }
        }
        for (int ist = 0; ist < surface_nst; ++ist) assert(st_flag[ist] == 1);
	delete[] st_flag;
        cout << " > GCL 2: surface " << COUT_VEC(gcl) << endl;
      }

      double (*gcl)[3] = new double[ncv][3];
      FOR_ICV FOR_I3 gcl[icv][i] = 0.0;

      FOR_IBF {
	const int icv = cvobf_global[ibf]; assert((icv >= 0)&&(icv < ncv));
	FOR_I3 gcl[icv][i] += n_bf[ibf][i];
      }

      FOR_IFA {
	const int icv0 = cvofa_global[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
	const int icv1 = cvofa_global[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
	FOR_I3 gcl[icv0][i] += n_fa[ifa][i];
	FOR_I3 gcl[icv1][i] -= n_fa[ifa][i];
      }

      MiscUtils::dumpRange(gcl,ncv,"GCL: compact faces");

      delete[] gcl;

      MPI_Sync("DO CHECK HERE");
      CWARN("should check the spobf_i/v/wgt reconstruction operators here."); 
      
      /*
	if (ncvobb) {
        for (int ibb = 0; ibb != nbb; ++ibb) {
	cout << "bbox: " << ibb << " ncvobb = " << ncvobb[ibb] << " vv_bbox = " << 
	vv_bbox[ibb][0] << " " << vv_bbox[ibb][1] << " " << 
	vv_bbox[ibb][2] << " " << vv_bbox[ibb][3] << " " << 
	vv_bbox[ibb][4] << " " << vv_bbox[ibb][5] << endl;
        }
	}
      */

      // also compare the n_fa and the node-based normal...
      
      double my_buf[2] = { 0.0, 0.0 };
      for (int ifa = 0; ifa < nfa; ++ifa) {
	double n_fa_check[3] = { 0.0, 0.0, 0.0 };
	const int ino = noofa_v_global[noofa_i[ifa]];
	int ino1 = noofa_v_global[noofa_i[ifa+1]-1];
	for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
	  const int ino0 = ino1;
	  ino1 = noofa_v_global[nof];
	  if ((ino != ino0)&&(ino != ino1)) {
	    const double this_n[3] = TRI_NORMAL_2(x_no[ino],x_no[ino0],x_no[ino1]);
	    FOR_I3 n_fa_check[i] += 0.5*this_n[i];
	  }
	}
	const double n_mag = MAG(n_fa[ifa]);
	const double n_mag_check = MAG(n_fa_check);
	
	//cout << "ifa: " << ifa << " mag: " << n_mag << " " << n_mag_check << " err: " <<  fabs(	1.0 - n_mag_check/n_mag ) << " n_fa: " << COUT_VEC(n_fa[ifa]) << " " << COUT_VEC(n_fa_check) << endl;
	
	my_buf[0] = max( my_buf[0], fabs( n_mag_check - n_mag ) );
	const double dp_check = DOT_PRODUCT(n_fa[ifa],n_fa_check)/(n_mag*n_mag_check);
	my_buf[1] = max( my_buf[1], fabs(	1.0 - dp_check ) );

	//if ( my_buf[0] > 1.0E-12 )
	//MPI_Pause("OKOK");
      }
      
      if (mpi_rank == 0) 
	cout << " > node face normal check (should be small): " << my_buf[0] << " " << my_buf[1] << endl;
      
    }
    
  }

  //MiscUtils

  void buildExtendedFaces() {

    COUT1("StripedMesh::buildExtendedFaces()");
    
    assert(nef == 0);
    assert(nef_global == 0);
    assert(cvoef_global == NULL);
    assert(n_ef == NULL);
    assert(c_ef == NULL);
    assert(efora == NULL);
    if (x_fa == NULL) {
      CERR("rebuild your restart.mles with the current stitch. Missing x_fa needed for new extended operator build.");
    }

    //  ======================================
    //  step 1: 
    //  we need to reduce compact x_fa, n_fa, cvofa to the current striped cv distribution
    //  on sm.
    //  ======================================
    int * send_count = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;
    
    for (int ifa = 0; ifa < nfa; ++ifa) {
      const int8 icv0_global = cvofa_global[ifa][0]; assert((icv0_global >= 0)&&(icv0_global < ncv_global));
      const int rank0 = MiscUtils::getRankInXora(icv0_global,cvora);
      const int8 icv1_global = (cvofa_global[ifa][1]&MASK_55BITS); assert((icv1_global >= 0)&&(icv1_global < ncv_global));
      const int bits = (cvofa_global[ifa][1]>>55);
      const int rank1 = MiscUtils::getRankInXora(icv1_global,cvora);
      if ((rank0 == rank1)||bits) {
        send_count[rank0] += 2;
      }
      else {
      	send_count[rank0] += 2;
      	send_count[rank1] += 2;
      }
    }

    int * send_disp = new int[mpi_size];
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
    
    uint8 * send_buf_uint8 = new uint8[send_count_sum];
    double (*send_buf_double)[3] = new double[send_count_sum][3];
    
    for (int ifa = 0; ifa < nfa; ++ifa) {
      const int8 icv0_global = cvofa_global[ifa][0]; assert((icv0_global >= 0)&&(icv0_global < ncv_global));
      const int rank0 = MiscUtils::getRankInXora(icv0_global,cvora);
      const int8 icv1_global = (cvofa_global[ifa][1]&MASK_55BITS); assert((icv1_global >= 0)&&(icv1_global < ncv_global));
      const int bits = (cvofa_global[ifa][1]>>55);
      const int rank1 = MiscUtils::getRankInXora(icv1_global,cvora);
      //    if there are bits, this must be a face in the periodic range 
      if (bits) {
      	assert(faora[mpi_rank]+ifa < nfa_p_global);
      }
      else {
      	assert(faora[mpi_rank]+ifa >= nfa_p_global);
      }
      if ((rank0 == rank1)||bits) {
        //      faces with bits (i.e. periodic faces) already have their duplicates in the grid, so
        //      do not need to be flipped and sent to rank1...
      	send_buf_uint8[send_disp[rank0]  ] = BitUtils::packRankBitsIndex(rank0,0,icv0_global-cvora[rank0]);
      	send_buf_uint8[send_disp[rank0]+1] = BitUtils::packRankBitsIndex(rank1,bits,icv1_global-cvora[rank1]);
      	FOR_I3 send_buf_double[send_disp[rank0]  ][i] = n_fa[ifa][i];
      	FOR_I3 send_buf_double[send_disp[rank0]+1][i] = x_fa[ifa][i];
      	send_disp[rank0] += 2;
      }
      else {
      	send_buf_uint8[send_disp[rank0]  ] = BitUtils::packRankBitsIndex(rank0,0,icv0_global-cvora[rank0]);
      	send_buf_uint8[send_disp[rank0]+1] = BitUtils::packRankBitsIndex(rank1,0,icv1_global-cvora[rank1]);
      	FOR_I3 send_buf_double[send_disp[rank0]  ][i] = n_fa[ifa][i];
      	FOR_I3 send_buf_double[send_disp[rank0]+1][i] = x_fa[ifa][i];
      	send_disp[rank0] += 2;
        //      flip the orientation when sending to rank1...
      	send_buf_uint8[send_disp[rank1]  ] = BitUtils::packRankBitsIndex(rank1,0,icv1_global-cvora[rank1]);
      	send_buf_uint8[send_disp[rank1]+1] = BitUtils::packRankBitsIndex(rank0,0,icv0_global-cvora[rank0]);
      	FOR_I3 send_buf_double[send_disp[rank1]  ][i] = -n_fa[ifa][i];
      	FOR_I3 send_buf_double[send_disp[rank1]+1][i] = x_fa[ifa][i];
      	send_disp[rank1] += 2;
      }
    }
    
    //  rewind
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

    //  exchange...
    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    
    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    
    uint8 * recv_buf_uint8 = new uint8[recv_count_sum];
    MPI_Alltoallv(send_buf_uint8,send_count,send_disp,MPI_UINT8,
		  recv_buf_uint8,recv_count,recv_disp,MPI_UINT8,mpi_comm);
    delete[] send_buf_uint8; send_buf_uint8 = NULL;

    FOR_RANK {
      send_count[rank] *= 3;
      send_disp[rank] *= 3;
      recv_count[rank] *= 3;
      recv_disp[rank] *= 3;
    }
    
    double (*recv_buf_double)[3] = new double[recv_count_sum][3];
    MPI_Alltoallv((double*)send_buf_double,send_count,send_disp,MPI_DOUBLE,
		  (double*)recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
    delete[] send_buf_double; send_buf_double = NULL;

    assert(recv_count_sum%2 == 0);
    const int nfa_tmp = recv_count_sum/2;
    double (*n_fa_tmp)[3] = new double[nfa_tmp][3];
    double (*x_fa_tmp)[3] = new double[nfa_tmp][3];
    int (*cvofa_tmp)[2] = new int[nfa_tmp][2];
    
    map<const uint8,int> rbiMap_tmp;
    const int ncv_tmp = ncv;
    assert(ncv_tmp == cvora[mpi_rank+1]-cvora[mpi_rank]);
    int ncv_g_tmp = ncv;
    
    for (int ifa = 0; ifa < nfa_tmp; ++ifa) {
      int rank,bits,icv0;
      BitUtils::unpackRankBitsIndex(rank,bits,icv0,recv_buf_uint8[ifa*2  ]);
      assert(rank == mpi_rank);
      assert(bits == 0);
      assert((icv0 >= 0)&&(icv0 < ncv));
      cvofa_tmp[ifa][0] = icv0;
      int icv1;
      BitUtils::unpackRankBitsIndex(rank,bits,icv1,recv_buf_uint8[ifa*2+1]);
      if ((rank == mpi_rank)&&(bits == 0)) {
      	assert((icv1 >= 0)&&(icv1 < ncv));
      	cvofa_tmp[ifa][1] = icv1;
      }
      else {
        map<const uint8,int>::iterator iter = rbiMap_tmp.find(recv_buf_uint8[ifa*2+1]);
      	if (iter == rbiMap_tmp.end()) {
      	  rbiMap_tmp[recv_buf_uint8[ifa*2+1]] = cvofa_tmp[ifa][1] = ncv_g_tmp++;
      	}
      	else {
      	  cvofa_tmp[ifa][1] = iter->second;
      	}
      }
      //    and the geometry...
      FOR_I3 n_fa_tmp[ifa][i] = recv_buf_double[ifa*2  ][i];
      FOR_I3 x_fa_tmp[ifa][i] = recv_buf_double[ifa*2+1][i];
    }
    //cout << nfa_tmp << " " << ncv_g_tmp << endl;

    delete[] recv_buf_uint8;  recv_buf_uint8 = NULL;
    delete[] recv_buf_double; recv_buf_double = NULL;

    //  ======================================
    //  step 2 fill in the coordinates...
    //  copy x_cv for convience...
    //  ======================================
    double (*x_cv_tmp)[3] = new double[ncv_g_tmp][3];
    for (int icv = 0; icv < ncv_tmp; ++icv)
      FOR_I3 x_cv_tmp[icv][i] = x_cv[icv][i];
    
    const int ncv_g_only = ncv_g_tmp-ncv_tmp;
    int8 * icv_global = new int8[ncv_g_only];
    uint8 * rbi_g_tmp = new uint8[ncv_g_only];
    for (map<const uint8,int>::iterator iter = rbiMap_tmp.begin(); iter != rbiMap_tmp.end(); ++iter) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,iter->first);
      assert((rank != mpi_rank)||bits);
      rbi_g_tmp[iter->second-ncv_tmp] = iter->first;
      icv_global[iter->second-ncv_tmp] = cvora[rank] + index;
    }
    rbiMap_tmp.clear();
    
    //  now pull from the striped distribution...
    DistributedDataExchanger * dde = new DistributedDataExchanger(icv_global,ncv_g_only,cvora);
    dde->pull(x_cv_tmp+ncv_tmp,x_cv);
    delete dde;
    delete[] icv_global;
    
    //  apply periodic transforms where required...
    for (int icv_g = 0; icv_g < ncv_g_only; ++icv_g) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g_tmp[icv_g]);
      if (bits) {
        PeriodicData::periodicTranslate(x_cv_tmp+ncv_tmp+icv_g,1,bits);
      }
    }
    
    //  now check that the normals are correct wrt the faces...
    //  we can't check anything more strict because these are centroids x_cv, not
    //  Voronoi points x_vv... 
    for (int ifa = 0; ifa < nfa_tmp; ++ifa) {
      const int icv0 = cvofa_tmp[ifa][0];  
      const double dx_0f[3] = DIFF(x_fa_tmp[ifa],x_cv_tmp[icv0]);
      const int icv1 = cvofa_tmp[ifa][1];
      const double dx_f1[3] = DIFF(x_cv_tmp[icv1],x_fa_tmp[ifa]);
      if (!(DOT_PRODUCT(dx_0f,n_fa_tmp[ifa]) > 0.0) || !(DOT_PRODUCT(dx_f1,n_fa_tmp[ifa]) > 0.0)) {
        cout << " > warning, centroid is on other side of voronoi face. " << endl;
        //cout << mpi_rank << " dp about to fail. dx_0f: " << COUT_VEC(dx_0f) << " mag: " << MAG(dx_0f) << 
        //  " dx_f1: " << COUT_VEC(dx_f1) << " mag: " << MAG(dx_f1) << " n_fa_tmp[ifa]: " << COUT_VEC(n_fa_tmp[ifa]) << endl;
        //cout << mpi_rank << " x_cv_tmp[icv0]: " << COUT_VEC(x_cv_tmp[icv0]) << " x_cv_tmp[icv1]: " << COUT_VEC(x_cv_tmp[icv1]) << " x_fa_tmp[ifa]: " << COUT_VEC(x_fa_tmp[ifa]) << endl;
      }
      //assert(DOT_PRODUCT(dx_0f,n_fa_tmp[ifa]) > 0.0);
      //assert(DOT_PRODUCT(dx_f1,n_fa_tmp[ifa]) > 0.0);
    }

    //  finally, we need faocv_i_tmp, faocv_v_tmp...
    //  and cvocv_i_tmp, cvocv_v_tmp...
    int * faocv_i_tmp = new int[ncv_tmp+1];
    int * cvocv_i_tmp = new int[ncv_tmp+1];
    for (int icv = 0; icv < ncv_tmp; ++icv) {
      faocv_i_tmp[icv+1] = 0;
      cvocv_i_tmp[icv+1] = 1; // diagonal
    }
    for (int ifa = 0; ifa < nfa_tmp; ++ifa) {
      const int icv0 = cvofa_tmp[ifa][0];
      ++faocv_i_tmp[icv0+1];
      ++cvocv_i_tmp[icv0+1];
      const int icv1 = cvofa_tmp[ifa][1];
      if (icv1 < ncv_tmp) {
      	++faocv_i_tmp[icv1+1];
      	++cvocv_i_tmp[icv1+1];
      }
    }
    //  allocate...
    faocv_i_tmp[0] = 0;
    for (int icv = 0; icv < ncv_tmp; ++icv) faocv_i_tmp[icv+1] += faocv_i_tmp[icv];
    const int faocv_s_tmp = faocv_i_tmp[ncv_tmp];
    int * faocv_v_tmp = new int[faocv_s_tmp];
    cvocv_i_tmp[0] = 0;
    for (int icv = 0; icv < ncv_tmp; ++icv) cvocv_i_tmp[icv+1] += cvocv_i_tmp[icv];
    const int cvocv_s_tmp = cvocv_i_tmp[ncv_tmp];
    int * cvocv_v_tmp = new int[cvocv_s_tmp];
    //  and set...
    for (int icv = 0; icv < ncv_tmp; ++icv) {
      cvocv_v_tmp[cvocv_i_tmp[icv]++] = icv; // diagonal
    }
    for (int ifa = 0; ifa < nfa_tmp; ++ifa) {
      const int icv0 = cvofa_tmp[ifa][0];
      const int icv1 = cvofa_tmp[ifa][1];
      faocv_v_tmp[faocv_i_tmp[icv0]++] = ifa;
      cvocv_v_tmp[cvocv_i_tmp[icv0]++] = icv1;
      if (icv1 < ncv_tmp) {
        faocv_v_tmp[faocv_i_tmp[icv1]++] = ifa;
      	cvocv_v_tmp[cvocv_i_tmp[icv1]++] = icv0;
      }
    }
    //  rewind...    
    for (int icv = ncv_tmp; icv > 0; --icv)
      faocv_i_tmp[icv] = faocv_i_tmp[icv-1];
    faocv_i_tmp[0] = 0;
    for (int icv = ncv_tmp; icv > 0; --icv)
      cvocv_i_tmp[icv] = cvocv_i_tmp[icv-1];
    cvocv_i_tmp[0] = 0;
    /*
      for (int icv = 0; icv < ncv_tmp; ++icv) {
      cout << "coc: " << endl;
      for (int coc = cvocv_i_tmp[icv]; coc != cvocv_i_tmp[icv+1]; ++coc) 
      cout << cvocv_v_tmp[coc] << endl;
      cout << "foc: " << endl;
      for (int foc = faocv_i_tmp[icv]; foc != faocv_i_tmp[icv+1]; ++foc) 
      cout << faocv_v_tmp[foc] << endl;
      }
    */
    
    //  ============================================
    //  step 3
    //  we need to reduce cvobf_global, Gij_bf, to the current striped cv distribution on sm.
    //  ============================================
    FOR_RANK send_count[rank] = 0;

    for (int ibf = 0; ibf < nbf; ++ibf) {
      const int rank = MiscUtils::getRankInXora(cvobf_global[ibf],cvora);
      send_count[rank] += 1;
    }

    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
    
    int * send_buf_int = new int[send_count_sum];
    double * send_buf_double_2 = new double[send_count_sum*9];
    
    for (int ibf = 0; ibf < nbf; ++ibf) {
      const int rank = MiscUtils::getRankInXora(cvobf_global[ibf],cvora);
      send_buf_int[send_disp[rank]] = cvobf_global[ibf] - cvora[rank];
      FOR_I3 FOR_J3 send_buf_double_2[send_disp[rank]*9+(3*i+j)] = Gij_bf[ibf][i][j];
      send_disp[rank] += 1;
    }
    
    //  rewind
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

    //  exchange...
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    const int nbf_tmp = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    //  exchange directly into the nbf_tmp data...
    int *cvobf_tmp = new int[nbf_tmp];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
		  cvobf_tmp,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; send_buf_int = NULL;
  
    FOR_RANK {
      send_count[rank] *= 9;
      send_disp[rank] *= 9;
      recv_count[rank] *= 9;
      recv_disp[rank] *= 9;
    }
    
    double (*Gij_bf_tmp)[3][3] = new double[nbf_tmp][3][3];
    MPI_Alltoallv(send_buf_double_2,send_count,send_disp,MPI_DOUBLE,
		  (double*)Gij_bf_tmp,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
    delete[] send_buf_double_2; send_buf_double_2 = NULL;

    //  bfocv_i_tmp/bfocv_v_tmp...
    int * bfocv_i_tmp = new int[ncv_tmp+1];
    for (int icv = 0; icv < ncv_tmp; ++icv) {
      bfocv_i_tmp[icv+1] = 0;
    }
    for (int ibf = 0; ibf < nbf_tmp; ++ibf) {
      const int icv = cvobf_tmp[ibf];
      ++bfocv_i_tmp[icv+1];
    }
    //  allocate...
    bfocv_i_tmp[0] = 0;
    for (int icv = 0; icv < ncv_tmp; ++icv) bfocv_i_tmp[icv+1] += bfocv_i_tmp[icv];
    const int bfocv_s_tmp = bfocv_i_tmp[ncv_tmp];
    int * bfocv_v_tmp = new int[bfocv_s_tmp];
    //  and set...
    for (int ibf = 0; ibf < nbf_tmp; ++ibf) {
      const int icv = cvobf_tmp[ibf];
      bfocv_v_tmp[bfocv_i_tmp[icv]++] = ibf;
    }
    //  rewind...    
    for (int icv = ncv_tmp; icv > 0; --icv)
      bfocv_i_tmp[icv] = bfocv_i_tmp[icv-1];
    bfocv_i_tmp[0] = 0;

    //  some cleanup...
  
    delete[] cvobf_tmp;
  
    //  ------------------------------------
    //  step 4: build the grad operator...
    //  ------------------------------------
    double (*cvocv_grad_coeff_tmp)[3] = new double[cvocv_i_tmp[ncv_tmp]][3];
  
    int8 nthreshold = 0;
    double my_sigma[2] = {1.0e+20,0.0};
    const double sigma_epsilon = getDoubleParam("SIGMA_THRESH", 4.5e-1);
    vector<double> sigma_thresh;

    int * cv_flag = new int[ncv_tmp];

    double my_svd_eps_tol = 0.0;

    for (int icv = 0; icv < ncv_tmp; ++icv) {

      const int coc_f = cvocv_i_tmp[icv];
      const int coc_l = cvocv_i_tmp[icv+1]-1;
      for (int coc = coc_f; coc <= coc_l; ++coc)
      	FOR_I3 cvocv_grad_coeff_tmp[coc][i] = 0.0;
     
      double A[9]; 
      for (int ii = 0; ii < 9; ++ii) 
        A[ii] = 0.0;

      //    start with volume on the diagonal
      for (int i =0; i < 3; ++i) 
        A[(i)*3+i] = vol_cv[icv];

      //cout << vol_cv[icv] << endl;
      for (int boc = bfocv_i_tmp[icv]; boc != bfocv_i_tmp[icv+1]; ++boc) { 
        const int ibf          = bfocv_v_tmp[boc];
        //cout << " boc: " << boc  << " ibf: " << ibf << endl;
        for (int j =0; j < 3; ++j) { 
          //cout << COUT_VEC(Gij_bf_tmp[ibf][j]) << endl;
          for (int i =0; i < 3; ++i) { 
            A[(j)*3+i] += Gij_bf_tmp[ibf][i][j];
          }
        }
      }
   
      for (int foc = faocv_i_tmp[icv]; foc != faocv_i_tmp[icv+1]; ++foc) {
        const int ifa = faocv_v_tmp[foc];
      	double fa_sign;
      	int icv_nbr;
      	if ( cvofa_tmp[ifa][0] == icv ) {
      	  fa_sign = 1.0;
      	  icv_nbr = cvofa_tmp[ifa][1];
      	}
      	else {
      	  assert( cvofa_tmp[ifa][1] == icv );
      	  fa_sign = -1.0;
      	  icv_nbr = cvofa_tmp[ifa][0];
      	}
        // also get the nbr associated with this icv_nbr...
      	const int coc = coc_f+foc-faocv_i_tmp[icv]+1;
      	assert(cvocv_v_tmp[coc] == icv_nbr);
      	double r[3]; FOR_I3 r[i] = x_fa_tmp[ifa][i] - 0.5*(x_cv_tmp[icv][i] + x_cv_tmp[icv_nbr][i]);

        for (int j =0; j < 3; ++j) { 
          for (int i =0; i < 3; ++i) { 
            A[(j)*3+i] -= r[j]*n_fa_tmp[ifa][i]*fa_sign;
          }
        }

        FOR_I3 cvocv_grad_coeff_tmp[coc][i]   += 0.5*n_fa_tmp[ifa][i]*fa_sign;
      	FOR_I3 cvocv_grad_coeff_tmp[coc_f][i] -= 0.5*n_fa_tmp[ifa][i]*fa_sign;
      }

      // rescale the matrices with the inv_vol to help the conditioning of the system...
      const double inv_vol_cv = 1.0/vol_cv[icv];
      for (int j = 0; j < 3; ++j) { 
        for (int i =0; i < 3; ++i) { 
          A[(j)*3+i] *= inv_vol_cv; 
        }
      }

      for (int coc = coc_f; coc <= coc_l; ++coc) {
        for (int i =0; i < 3; ++i) { 
          cvocv_grad_coeff_tmp[coc][i] *= inv_vol_cv;
        }
      }

      double U[9], V[9], sigma[3];
      double eps_tol;
      calcSvd33(U,V,sigma,A,eps_tol);

      my_svd_eps_tol = max(eps_tol,my_svd_eps_tol);

      if ( eps_tol > 0.1) { 
        char filename[128]; sprintf(filename,"A.%06d.dat", icv);
        FILE * fp = fopen(filename,"w");
        for (int i =0; i < 3; ++i) 
          fprintf(fp,"%12.15g  %12.15g  %12.15g\n", A[(0)*3+i], A[(1)*3+i], A[(2)*3+i]);
        fclose(fp);


        sprintf(filename, "A.%06d.binary.dat",icv);
        FILE * fp2 = fopen(filename,"w");
        fwrite(A,sizeof(double),9,fp2);
        fclose(fp2);
      }

      cv_flag[icv] = 0; // assume that this cell wont be thresholded.. 

      //    the singular values are not necessarily sorted from the prior routine..
      //    so we need to loop and build the thresholding here... 
      
      double sigma_inv[3];
      for (int i =0; i < 3; ++i) {
        /*
          if ( (sigma[i] < 0.0) || (sigma[i] != sigma[i])) { 
          cout << "uh oh: " << sigma[i] << endl;
          FOR_I3 { 
          FOR_J3 cout << A[(j)*3+i] << "   ";
          cout << endl;
          }

          cout << "sigma: " << COUT_VEC(sigma) << endl;
          }
        */
        assert( sigma[i] >= 0.0);
        my_sigma[0] = min(sigma[i], my_sigma[0]);
        if ( sigma[i] < sigma_epsilon) { 
          sigma_inv[i] = 0.0; // cant trust the gradient in this dir.. 
          cv_flag[icv] = 1;
          my_sigma[1]  = min(-sigma[i],my_sigma[1]);
        } 
        else { 
          sigma_inv[i] = 1.0/sigma[i];
        }
      }

      if ( cv_flag[icv] == 1) { 
        ++nthreshold;
      } 
      else { 
        assert( cv_flag[icv] == 0);
        // this is a somewhat costly check, but for now, if we didnt threshold
        // the singular values, check to ensure that we have decomposed to matrix 
        // properly.  this check can be removed at a later date ... 
        
        double max_err = 0.0;
        for (int j = 0; j < 3; ++j) { 
          for (int i =0; i < 3; ++i) { 

            double aij = 0.0;
            for (int k =0; k < 3; ++k) 
              aij += U[(k)*3+i]*V[(k)*3+j]*sigma[k];
            
            max_err = max(max_err,abs(aij-A[(j)*3+i]));
          }
        }

        // since the matrix is now O(1), this tolerance should be properly scaled.
        assert( max_err < max(1.0e-12,eps_tol));
      }

      // now supply the regularized inverse .. let V <--- V\Sigma^+
      for (int j = 0; j < 3; ++j) { 
        for (int i =0; i < 3; ++i) { 
          V[(j)*3+i] *= sigma_inv[j];
        }
      }

      // finally supply the inverse for the gradient coefficients ..
      for (int coc = coc_f; coc <= coc_l; ++coc) { 
        double b_tmp[3];
        for (int i =0; i < 3; ++i) { 
          b_tmp[i] = 0.0;
          for (int j =0; j < 3; ++j) { 
            b_tmp[i] += A[(j)*3+i]*cvocv_grad_coeff_tmp[coc][j];
          }
        }
	
        double tmp[3];
        //  V*U^T[r]
        for (int i =0; i <3 ; ++i)  
          tmp[i] = 
      	    U[(i)*3+0]*(cvocv_grad_coeff_tmp[coc][0] - b_tmp[0]) + 
      	    U[(i)*3+1]*(cvocv_grad_coeff_tmp[coc][1] - b_tmp[1]) +
      	    U[(i)*3+2]*(cvocv_grad_coeff_tmp[coc][2] - b_tmp[2]); 

        for (int i =0; i < 3; ++i) 
          cvocv_grad_coeff_tmp[coc][i] += V[0*3+i]*tmp[0] + 
            V[1*3+i]*tmp[1] + 
            V[2*3+i]*tmp[2];

      }

      bool debug = false;
      for (int coc = coc_f; coc <= coc_l; ++coc) { 
        for (int i =0; i < 3; ++i) { 
          if ( cvocv_grad_coeff_tmp[coc][i] != cvocv_grad_coeff_tmp[coc][i] ) { 
            debug = true;
          }
        }
      }

      if ( debug) { 
        char filename[128]; sprintf(filename,"A.%06d.dat", icv);
        FILE * fp = fopen(filename,"w");
        for (int i =0; i < 3; ++i) 
          fprintf(fp,"%12.15g  %12.15g  %12.15g\n", A[(0)*3+i], A[(1)*3+i], A[(2)*3+i]);
        fclose(fp);


        sprintf(filename, "A.%06d.binary.dat",icv);
        FILE * fp2 = fopen(filename,"w");
        fwrite(A,sizeof(double),9,fp2);
        fclose(fp2);
      }

    }

    delete[] Gij_bf_tmp;
    
    /*    
          if (mpi_size == 1) {

          const double exact_grad[3] = { 1.123, 2.134, -1.3423 };
    
          double * phi = new double[ncv_tmp];
          for (int icv = 0; icv < ncv_tmp; ++icv) 
          phi[icv] = DOT_PRODUCT(x_cv_tmp[icv],exact_grad);
    
          double (*grad_phi)[3] = new double[ncv_tmp][3];
          for (int icv = 0; icv < ncv_tmp; ++icv) { 
          FOR_I3 grad_phi[icv][i] = 0.0;
          for (int coc = cvocv_i_tmp[icv]; coc != cvocv_i_tmp[icv+1]; ++coc) { 
          const int icv_nbr = cvocv_v_tmp[coc];
          FOR_I3 grad_phi[icv][i] += vol_cv[icv]*cvocv_grad_coeff_tmp[coc][i] * (phi[icv_nbr] - phi[icv]);
          }
          }
    
          for (int ibf = 0; ibf < nbf; ++ibf) {
          const double phi_bf = DOT_PRODUCT(exact_grad,x_bf[ibf]);
          const int icv = cvobf_global[ibf];
          FOR_I3 grad_phi[icv][i] += (phi_bf - phi[icv])*n_bf[ibf][i];
          }
    
          for (int icv = 0; icv < ncv_tmp; ++icv) {
          FOR_I3 grad_phi[icv][i] = grad_phi[icv][i]/(vol_cv[icv]*exact_grad[i]) - 1.0; // to produce 0
          }
    
          FILE * fp = fopen("grad.dat","w");
          for (int icv = 0; icv < ncv_tmp; ++icv) {
          fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %18.15le %18.15le\n",
          x_cv_tmp[icv][0],x_cv_tmp[icv][1],x_cv_tmp[icv][2],
          grad_phi[icv][0],grad_phi[icv][1],grad_phi[icv][2]);
          }
          fclose(fp);

          // internal only...

          for (int ibf = 0; ibf < nbf; ++ibf) {
          const int icv = cvobf_global[ibf];
          FOR_I3 grad_phi[icv][i] = 0.0;
          }
    
          MiscUtils::dumpRange(grad_phi,ncv_tmp,"extended grad - internal only");
    
          delete[] grad_phi;
          delete[] phi;
    
          }
    */

    // reporting...

    { 

      int8 ng_threshold;
      MPI_Allreduce(&nthreshold,&ng_threshold,1,MPI_INT8,MPI_SUM,mpi_comm);

      double sigma_buf[2];
      MPI_Reduce(my_sigma,sigma_buf,2,MPI_DOUBLE,MPI_MIN,0,mpi_comm);

      double svd_eps_tol;
      MPI_Reduce(&my_svd_eps_tol,&svd_eps_tol,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);

      int8 ng_threshold_internal;
      int8 tmp =0;
      for (int icv = 0; icv < ncv_tmp; ++icv) { 
        const int nboc = bfocv_i_tmp[icv+1] - bfocv_i_tmp[icv];
        if ( (nboc == 0) && (cv_flag[icv] == 1)) 
          ++tmp;
      }

      MPI_Allreduce(&tmp,&ng_threshold_internal,1,MPI_INT8,MPI_SUM,mpi_comm);

      if ( mpi_rank == 0 ) { 

        cout << " > threshold grad calc in (total, internal): ( " 
             << ng_threshold << "    " << ng_threshold_internal << " ) cells of : " << ncv_global << endl;
        cout << " > min normalized sigma   : " << sigma_buf[0] << endl;
        cout << " > max thresholded sigma  : " << -sigma_buf[1] << endl;
        cout << " > svd eps tol (should be small): " << svd_eps_tol << endl;
      }
    }

    delete[] bfocv_i_tmp;
    delete[] bfocv_v_tmp;

    //  check compact gradient
    {
      
      const double grad_phi_exact[3] = { 1.1234, -1.3243, 1.5321}; 
      double my_err_max[3] = { 0.0, 0.0, 0.0 };
      double my_err_max2[3] = {0.0, 0.0, 0.0 };

      double (*grad_phi_cv)[3] = new double[ncv_tmp][3];

      for (int icv = 0; icv < ncv_tmp; ++icv) {
      	const double phi_icv = DOT_PRODUCT(x_cv_tmp[icv],grad_phi_exact);
      	double grad_phi[3] = { 0.0, 0.0, 0.0 };
      	for (int coc = cvocv_i_tmp[icv]; coc != cvocv_i_tmp[icv+1]; ++coc) {
      	  const int icv_nbr = cvocv_v_tmp[coc];
      	  const double phi_icv_nbr = DOT_PRODUCT(x_cv_tmp[icv_nbr],grad_phi_exact);
      	  FOR_I3 grad_phi[i] += cvocv_grad_coeff_tmp[coc][i]*(phi_icv_nbr-phi_icv);
      	}
      	FOR_I3 {
      	  grad_phi[i] -= grad_phi_exact[i];
      	  grad_phi[i] /= grad_phi_exact[i];
          grad_phi_cv[icv][i] = grad_phi[i];
      	  const double err = fabs(grad_phi[i]);
      	  my_err_max[i] = max(my_err_max[i],err);
        }

        if ( cv_flag[icv] == 0) { 
          for (int i =0; i < 3; ++i) { 
            my_err_max2[i] = max(my_err_max2[i],fabs(grad_phi[i]));
          }
        } else { 
          FOR_I3 grad_phi_cv[icv][i] = 0.0;
        }
    
      
      }

      MiscUtils::dumpRange(grad_phi_cv, ncv_tmp, "compact gradient error (range)");
      delete[] grad_phi_cv;

      double err_max[3];
      MPI_Reduce(my_err_max,err_max,3,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      if (mpi_rank == 0)
      	cout << " > compact gradient error -- including thresholded cvs (should be order one): " << COUT_VEC(err_max) << endl;

      MPI_Reduce(my_err_max2,err_max,3,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      if ( mpi_rank == 0)
        cout << " > compact gradient error -- non-thresholded cvs (should be small): " << COUT_VEC(err_max) << endl;

    }

    // report compact nbr complexity before deleting cv_flag...
    
    {
      for (int icv = 0; icv < ncv_tmp; ++icv) { 
        cv_flag[icv] = faocv_i_tmp[icv+1]-faocv_i_tmp[icv];
      }
      MiscUtils::dumpBins(cv_flag,ncv_tmp,"Compact nbr count (not including self)");
    }
    
    delete[] cv_flag;

    // ---------------------------------------------
    // now try to build the operators...
    // EF_OPS param controls method:
    // EF_OPS DEFAULT
    // EF_OPS ALIGNED
    // EF_OPS COMPACT_SKEW
    // EF_OPS COMPACT_LINEAR
    //
    // for the first 2 cases, the follow param is
    // also available:
    // EF_OPS.INTERP_FACTOR <double> with default 1.0/3.0
    // 
    // from the euler tests, a lower value improves
    // the error, e.g.
    // EF_OPS.INTERP_FACTOR 0.273
    // ---------------------------------------------
    
    vector<Coeff> coeffVec;
    
    {

      // lower memory build of the coeffVec
      map<pair<uint8,uint8>,int,rbiPairComparator> cvec_idx;
      cvec_idx.clear();

      const string ef_ops = getStringParam("EF_OPS","DEFAULT");
      if (ef_ops == "DEFAULT") {
      
        const double interp_factor = getDoubleParam("EF_OPS.INTERP_FACTOR",1.0/3.0);

        // ===================================================
        // this is the original operator reconstruction, except
        // the interp_factor controls the relative weighting between
        // the compact centered term and the gradient-correction
        // term. 
        // ===================================================

        for (int iter = 0; iter < 2; ++iter) { 
          for (int icv = 0; icv < ncv_tmp; ++icv) { 
            const uint8 rbi = BitUtils::packRankBitsIndex(mpi_rank,0,icv);
            for (int foc = faocv_i_tmp[icv]; foc != faocv_i_tmp[icv+1]; ++foc) {
              const int ifa = faocv_v_tmp[foc];
              double fa_sign;
              int icv_nbr;
              uint8 rbi_nbr;
              if (icv == cvofa_tmp[ifa][0]) {
                fa_sign = 1.0;
                icv_nbr = cvofa_tmp[ifa][1];
                // nbr could be ghost...
                if (icv_nbr >= ncv_tmp) 
                  rbi_nbr = rbi_g_tmp[icv_nbr-ncv_tmp];
                else
                  rbi_nbr = BitUtils::packRankBitsIndex(mpi_rank,0,icv_nbr);
              }
              else {
                assert(icv == cvofa_tmp[ifa][1]);
                fa_sign = -1.0;
                icv_nbr = cvofa_tmp[ifa][0];
                // nbr is always local in this case...
                rbi_nbr = BitUtils::packRankBitsIndex(mpi_rank,0,icv_nbr);
              }
              // central part first - here do the full central part because we can...
              { 
                pair<uint8,uint8> ij_pair(rbi,rbi_nbr);
                map<pair<uint8,uint8>,int>::iterator it = cvec_idx.find(ij_pair);
                if ( iter == 0 ) { 
                  if ( it == cvec_idx.end())  
                    cvec_idx[ij_pair] = -1;
                }
                else { 
                  assert ( it != cvec_idx.end()); // coeff location must exist
                  assert ( it->second >= 0);
                  FOR_I3 coeffVec[it->second].coeff[i] += 0.5*fa_sign*n_fa_tmp[ifa][i];
                }
              }
              // now the gradient part...
              // note that interp_factor == 1/3 recovers the original default...
              double dx[3];
              //FOR_I3 dx[i] = 0.5*x_fa_tmp[ifa][i] - (5.0*x_cv_tmp[icv][i] + x_cv_tmp[icv_nbr][i])/12.0;
              FOR_I3 dx[i] = 0.5*x_fa_tmp[ifa][i] - (0.5-0.25*interp_factor)*x_cv_tmp[icv][i] - 0.25*interp_factor*x_cv_tmp[icv_nbr][i];
              for (int coc = cvocv_i_tmp[icv]; coc != cvocv_i_tmp[icv+1]; ++coc) {
                const int icv2 = cvocv_v_tmp[coc]; assert((icv2 >= 0)&&(icv2 < ncv_g_tmp));
                uint8 rbi2;
                if (icv2 < ncv_tmp)
                  rbi2 = BitUtils::packRankBitsIndex(mpi_rank,0,icv2);
                else
                  rbi2 = rbi_g_tmp[icv2-ncv_tmp];
                const double dp = DOT_PRODUCT(cvocv_grad_coeff_tmp[coc],dx);
              
                // (rbi,rbi2) ..  
                if ( rbi != rbi2) { 
                  pair<uint8,uint8> ij_pair(rbi,rbi2);
                  map<pair<uint8,uint8>,int>::iterator it = cvec_idx.find(ij_pair);
                  if ( iter == 0 ) { 
                    if ( it == cvec_idx.end())  
                      cvec_idx[ij_pair] = -1;
                  } 
                  else { 
                    assert( it != cvec_idx.end()); // (i,j) must exist.. 
                    assert( it->second >= 0);
                    FOR_I3 coeffVec[it->second].coeff[i] += dp*fa_sign*n_fa_tmp[ifa][i];
                  }
                }
              
                // (rbi_nbr,rbi2)
                if ( rbi_nbr != rbi2) { 
                  pair<uint8,uint8> ij_pair(rbi_nbr,rbi2);
                  map<pair<uint8,uint8>,int>::iterator it = cvec_idx.find(ij_pair);
                  if ( iter == 0 ) { 
                    if ( it == cvec_idx.end())  
                      cvec_idx[ij_pair] = -1;
                  } 
                  else { 
                    assert( it != cvec_idx.end()); // (i,j) must exist...
                    assert( it->second >= 0);
                    FOR_I3 coeffVec[it->second].coeff[i] -= dp*fa_sign*n_fa_tmp[ifa][i];
                  }
                } 
              }//coc
            }//foc
          }//icv
        
          if ( iter == 0 ) { 
          
            int coeff_size = cvec_idx.size();
            coeffVec.resize(coeff_size);
          
            int max_coeff_size;
            MPI_Reduce(&coeff_size,&max_coeff_size,1,MPI_INT,MPI_MAX,0,mpi_comm);
            if ( mpi_rank == 0 ) 
              cout << " > max new memory alloc size [GB]: " << double(max_coeff_size)*24.0/(1024.0*1024.0*1024.0) << endl;
          
            int ii = 0;
            for(map<pair<uint8,uint8>,int>::iterator it = cvec_idx.begin();
                it != cvec_idx.end(); ++it) { 
              it->second = ii;
              coeffVec[ii].rbi = it->first.first;
              coeffVec[ii].rbi_nbr = it->first.second;
              coeffVec[ii].transpose = false;
              FOR_I3 coeffVec[ii].coeff[i] = 0.0;
              ++ii;
            }
          }
        }//iter

      }
      else if (ef_ops == "ALIGNED") {
      
        const double interp_factor = getDoubleParam("EF_OPS.INTERP_FACTOR",1.0/3.0);

        // ===================================================
        // this is like the DEFAULT reconstruction, but it tries to
        // find the opposite nbr and create an essentially 1D stencil
        // when possible...
        // ===================================================
      
        for (int iter = 0; iter < 2; ++iter) { 
          for (int icv = 0; icv < ncv_tmp; ++icv) { 
            const uint8 rbi = BitUtils::packRankBitsIndex(mpi_rank,0,icv);
            for (int foc = faocv_i_tmp[icv]; foc != faocv_i_tmp[icv+1]; ++foc) {
              const int ifa = faocv_v_tmp[foc];
              double fa_sign;
              int icv_nbr;
              uint8 rbi_nbr;
              if (icv == cvofa_tmp[ifa][0]) {
                fa_sign = 1.0;
                icv_nbr = cvofa_tmp[ifa][1];
                // nbr could be ghost...
                if (icv_nbr >= ncv_tmp) 
                  rbi_nbr = rbi_g_tmp[icv_nbr-ncv_tmp];
                else
                  rbi_nbr = BitUtils::packRankBitsIndex(mpi_rank,0,icv_nbr);
              }
              else {
                assert(icv == cvofa_tmp[ifa][1]);
                fa_sign = -1.0;
                icv_nbr = cvofa_tmp[ifa][0];
                // nbr is always local in this case...
                rbi_nbr = BitUtils::packRankBitsIndex(mpi_rank,0,icv_nbr);
              }
              // central part first - here do the full central part because we can...
              { 
                pair<uint8,uint8> ij_pair(rbi,rbi_nbr);
                map<pair<uint8,uint8>,int>::iterator it = cvec_idx.find(ij_pair);
                if ( iter == 0 ) { 
                  if ( it == cvec_idx.end())  
                    cvec_idx[ij_pair] = -1;
                }
                else { 
                  assert ( it != cvec_idx.end()); // coeff location must exist
                  assert ( it->second >= 0);
                  FOR_I3 coeffVec[it->second].coeff[i] += 0.5*fa_sign*n_fa_tmp[ifa][i];
                }
              }

              // also try to find an opposite nbr...
              int icv_opp = -1;
              double dp_min = -0.75;
              const double dx_nbr[3] = DIFF(x_cv_tmp[icv_nbr],x_cv_tmp[icv]);
              const double dx_nbr_mag = MAG(dx_nbr); assert(dx_nbr_mag > 0.0);
              for (int foc2 = faocv_i_tmp[icv]; foc2 != faocv_i_tmp[icv+1]; ++foc2) {
                if (foc2 != foc) {
                  const int ifa2 = faocv_v_tmp[foc2];
                  int icv2;
                  if (icv == cvofa_tmp[ifa2][0]) {
                    icv2 = cvofa_tmp[ifa2][1];
                  }
                  else {
                    assert(icv == cvofa_tmp[ifa2][1]);
                    icv2 = cvofa_tmp[ifa2][0];
                  }
                  const double dx_opp[3] = DIFF(x_cv_tmp[icv2],x_cv_tmp[icv]);
                  const double dx_opp_mag = MAG(dx_opp); assert(dx_opp_mag > 0.0);
                  const double dp = DOT_PRODUCT(dx_nbr,dx_opp)/(dx_nbr_mag*dx_opp_mag);
                  if (dp < dp_min) {
                    icv_opp = icv2;
                    dp_min = dp;
                  }
                }
              }

              // now the gradient part...
              // note that interp_factor == 1/3 recovers the original default...
              double dx[3];
              //FOR_I3 dx[i] = 0.5*x_fa_tmp[ifa][i] - (5.0*x_cv_tmp[icv][i] + x_cv_tmp[icv_nbr][i])/12.0;
              FOR_I3 dx[i] = 0.5*x_fa_tmp[ifa][i] - (0.5-0.25*interp_factor)*x_cv_tmp[icv][i] - 0.25*interp_factor*x_cv_tmp[icv_nbr][i];
            
              if (icv_opp >= 0) {
              
                // we got an icv_opp. So use the simple approximation (icv_nbr-icv_opp)/delta for the 
                // gradient in that direction...
              
                uint8 rbi_opp;
                if (icv_opp >= ncv_tmp) 
                  rbi_opp = rbi_g_tmp[icv_opp-ncv_tmp];
                else
                  rbi_opp = BitUtils::packRankBitsIndex(mpi_rank,0,icv_opp);
              
                // the dx_aligned is nbr minus opp...
                const double dx_aligned[3] = DIFF(x_cv_tmp[icv_nbr],x_cv_tmp[icv_opp]);
                const double dx_aligned_mag = MAG(dx_aligned); assert(dx_aligned_mag > 0.0);
              
                // take the component of dx on dx_aligned and use the simple approximation
                // to estimate the gradient operator...
              
                const double dp = DOT_PRODUCT(dx,dx_aligned)/dx_aligned_mag;
              
                // in rbi's matrix row, we need to add a term at rbi_nbr and rbi_opp
              
                if (rbi != rbi_nbr) { 
                  pair<uint8,uint8> ij_pair(rbi,rbi_nbr);
                  map<pair<uint8,uint8>,int>::iterator it = cvec_idx.find(ij_pair);
                  if ( iter == 0 ) { 
                    if ( it == cvec_idx.end())  
                      cvec_idx[ij_pair] = -1;
                  } 
                  else { 
                    assert( it != cvec_idx.end()); // (i,j) must exist.. 
                    assert( it->second >= 0);
                    FOR_I3 coeffVec[it->second].coeff[i] += dp/dx_aligned_mag*fa_sign*n_fa_tmp[ifa][i];
                  }
                }
                if (rbi != rbi_opp) { 
                  pair<uint8,uint8> ij_pair(rbi,rbi_opp);
                  map<pair<uint8,uint8>,int>::iterator it = cvec_idx.find(ij_pair);
                  if ( iter == 0 ) { 
                    if ( it == cvec_idx.end())  
                      cvec_idx[ij_pair] = -1;
                  } 
                  else { 
                    assert( it != cvec_idx.end()); // (i,j) must exist.. 
                    assert( it->second >= 0);
                    FOR_I3 coeffVec[it->second].coeff[i] -= dp/dx_aligned_mag*fa_sign*n_fa_tmp[ifa][i];
                  }
                }
              
                // in rbi_nbr's row, we need just the rbi_opp term (because the rbi_nbr
                // term is a diagonal term, obviously)...
                if (rbi_nbr != rbi_opp) { 
                  pair<uint8,uint8> ij_pair(rbi_nbr,rbi_opp);
                  map<pair<uint8,uint8>,int>::iterator it = cvec_idx.find(ij_pair);
                  if ( iter == 0 ) { 
                    if ( it == cvec_idx.end())  
                      cvec_idx[ij_pair] = -1;
                  } 
                  else { 
                    assert( it != cvec_idx.end()); // (i,j) must exist.. 
                    assert( it->second >= 0);
                    FOR_I3 coeffVec[it->second].coeff[i] += dp/dx_aligned_mag*fa_sign*n_fa_tmp[ifa][i];
                  }
                }
              
                // finally, removed the aligned component from dx so the below full gradient part can
                // act on the remainder...
              
                FOR_I3 dx[i] -= dp*dx_aligned[i]/dx_aligned_mag;
              
                // if the remaining dx is tiny wrt dx_aligned_mag, then formally skip the full gradient 
                // part. This will leave the coeff matrix more sparse...
              
                if (MAG(dx) < 1.0E-8*dx_aligned_mag) 
                  continue;
              
              }
            
              // any "dx" left treat with the full gradient...
            
              for (int coc = cvocv_i_tmp[icv]; coc != cvocv_i_tmp[icv+1]; ++coc) {
                const int icv2 = cvocv_v_tmp[coc]; assert((icv2 >= 0)&&(icv2 < ncv_g_tmp));
                uint8 rbi2;
                if (icv2 < ncv_tmp)
                  rbi2 = BitUtils::packRankBitsIndex(mpi_rank,0,icv2);
                else
                  rbi2 = rbi_g_tmp[icv2-ncv_tmp];
                const double dp = DOT_PRODUCT(cvocv_grad_coeff_tmp[coc],dx);
              
                // (rbi,rbi2) ..  
                if ( rbi != rbi2) { 
                  pair<uint8,uint8> ij_pair(rbi,rbi2);
                  map<pair<uint8,uint8>,int>::iterator it = cvec_idx.find(ij_pair);
                  if ( iter == 0 ) { 
                    if ( it == cvec_idx.end())  
                      cvec_idx[ij_pair] = -1;
                  } 
                  else { 
                    assert( it != cvec_idx.end()); // (i,j) must exist.. 
                    assert( it->second >= 0);
                    FOR_I3 coeffVec[it->second].coeff[i] += dp*fa_sign*n_fa_tmp[ifa][i];
                  }
                }
              
                // (rbi_nbr,rbi2)
                if ( rbi_nbr != rbi2) { 
                  pair<uint8,uint8> ij_pair(rbi_nbr,rbi2);
                  map<pair<uint8,uint8>,int>::iterator it = cvec_idx.find(ij_pair);
                  if ( iter == 0 ) { 
                    if ( it == cvec_idx.end())  
                      cvec_idx[ij_pair] = -1;
                  } 
                  else { 
                    assert( it != cvec_idx.end()); // (i,j) must exist...
                    assert( it->second >= 0);
                    FOR_I3 coeffVec[it->second].coeff[i] -= dp*fa_sign*n_fa_tmp[ifa][i];
                  }
                } 
              }//coc
            }//foc
          }//icv
        
          if ( iter == 0 ) { 
          
            int coeff_size = cvec_idx.size();
            coeffVec.resize(coeff_size);
          
            int max_coeff_size;
            MPI_Reduce(&coeff_size,&max_coeff_size,1,MPI_INT,MPI_MAX,0,mpi_comm);
            if ( mpi_rank == 0 ) 
              cout << " > max new memory alloc size [GB]: " << double(max_coeff_size)*24.0/(1024.0*1024.0*1024.0) << endl;
          
            int ii = 0;
            for(map<pair<uint8,uint8>,int>::iterator it = cvec_idx.begin();
                it != cvec_idx.end(); ++it) { 
              it->second = ii;
              coeffVec[ii].rbi = it->first.first;
              coeffVec[ii].rbi_nbr = it->first.second;
              coeffVec[ii].transpose = false;
              FOR_I3 coeffVec[ii].coeff[i] = 0.0;
              ++ii;
            }
          }
        }//iter

      }
      else if (ef_ops == "COMPACT_SKEW") {
      
        // ===================================================
        // 1/2*(L+R)...
        // ===================================================

        for (int iter = 0; iter < 2; ++iter) { 
          for (int icv = 0; icv < ncv_tmp; ++icv) { 
            const uint8 rbi = BitUtils::packRankBitsIndex(mpi_rank,0,icv);
            for (int foc = faocv_i_tmp[icv]; foc != faocv_i_tmp[icv+1]; ++foc) {
              const int ifa = faocv_v_tmp[foc];
              double fa_sign;
              int icv_nbr;
              uint8 rbi_nbr;
              if (icv == cvofa_tmp[ifa][0]) {
                fa_sign = 1.0;
                icv_nbr = cvofa_tmp[ifa][1];
                // nbr could be ghost...
                if (icv_nbr >= ncv_tmp) 
                  rbi_nbr = rbi_g_tmp[icv_nbr-ncv_tmp];
                else
                  rbi_nbr = BitUtils::packRankBitsIndex(mpi_rank,0,icv_nbr);
              }
              else {
                assert(icv == cvofa_tmp[ifa][1]);
                fa_sign = -1.0;
                icv_nbr = cvofa_tmp[ifa][0];
                // nbr is always local in this case...
                rbi_nbr = BitUtils::packRankBitsIndex(mpi_rank,0,icv_nbr);
              }
              // just the central part...
              { 
                pair<uint8,uint8> ij_pair(rbi,rbi_nbr);
                map<pair<uint8,uint8>,int>::iterator it = cvec_idx.find(ij_pair);
                if ( iter == 0 ) { 
                  if ( it == cvec_idx.end())  
                    cvec_idx[ij_pair] = -1;
                }
                else { 
                  assert ( it != cvec_idx.end()); // coeff location must exist
                  assert ( it->second >= 0);
                  FOR_I3 coeffVec[it->second].coeff[i] += 0.5*fa_sign*n_fa_tmp[ifa][i];
                }
              }
            }
          }
          if ( iter == 0 ) { 
            int coeff_size = cvec_idx.size();
            coeffVec.resize(coeff_size);
            int max_coeff_size;
            MPI_Reduce(&coeff_size,&max_coeff_size,1,MPI_INT,MPI_MAX,0,mpi_comm);
            if ( mpi_rank == 0 ) 
              cout << " > max new memory alloc size [GB]: " << double(max_coeff_size)*24.0/(1024.0*1024.0*1024.0) << endl;
            int ii = 0;
            for(map<pair<uint8,uint8>,int>::iterator it = cvec_idx.begin();
                it != cvec_idx.end(); ++it) { 
              it->second = ii;
              coeffVec[ii].rbi = it->first.first;
              coeffVec[ii].rbi_nbr = it->first.second;
              coeffVec[ii].transpose = false;
              FOR_I3 coeffVec[ii].coeff[i] = 0.0;
              ++ii;
            }
          }
        }//iter

      }
      else if (ef_ops == "COMPACT_LINEAR") {

        // ===================================================
        // (L,R) based in inverse distance from the face...
        // - minimally better than COMPACT_SKEW. 
        // ===================================================

        for (int iter = 0; iter < 2; ++iter) { 
          for (int icv = 0; icv < ncv_tmp; ++icv) { 
            const uint8 rbi = BitUtils::packRankBitsIndex(mpi_rank,0,icv);
            for (int foc = faocv_i_tmp[icv]; foc != faocv_i_tmp[icv+1]; ++foc) {
              const int ifa = faocv_v_tmp[foc];
              double fa_sign;
              int icv_nbr;
              uint8 rbi_nbr;
              if (icv == cvofa_tmp[ifa][0]) {
                fa_sign = 1.0;
                icv_nbr = cvofa_tmp[ifa][1];
                // nbr could be ghost...
                if (icv_nbr >= ncv_tmp) 
                  rbi_nbr = rbi_g_tmp[icv_nbr-ncv_tmp];
                else
                  rbi_nbr = BitUtils::packRankBitsIndex(mpi_rank,0,icv_nbr);
              }
              else {
                assert(icv == cvofa_tmp[ifa][1]);
                fa_sign = -1.0;
                icv_nbr = cvofa_tmp[ifa][0];
                // nbr is always local in this case...
                rbi_nbr = BitUtils::packRankBitsIndex(mpi_rank,0,icv_nbr);
              }
              // just the central part...
              { 
                // use dist to weight the values. When aligned this is a true linear
                // approimation. When there is some amount of transverse displacement,
                // this could be better conditioned. Also the dist never goes negative
                const double dist = DIST(x_cv_tmp[icv],x_fa_tmp[ifa]); 
                const double dist_nbr = DIST(x_cv_tmp[icv_nbr],x_fa_tmp[ifa]); 
                pair<uint8,uint8> ij_pair(rbi,rbi_nbr);
                map<pair<uint8,uint8>,int>::iterator it = cvec_idx.find(ij_pair);
                if ( iter == 0 ) { 
                  if ( it == cvec_idx.end())  
                    cvec_idx[ij_pair] = -1;
                }
                else { 
                  assert ( it != cvec_idx.end()); // coeff location must exist
                  assert ( it->second >= 0);
                  FOR_I3 coeffVec[it->second].coeff[i] += dist/(dist+dist_nbr)*fa_sign*n_fa_tmp[ifa][i];
                }
              }
            }
          }
          if ( iter == 0 ) { 
            int coeff_size = cvec_idx.size();
            coeffVec.resize(coeff_size);
            int max_coeff_size;
            MPI_Reduce(&coeff_size,&max_coeff_size,1,MPI_INT,MPI_MAX,0,mpi_comm);
            if ( mpi_rank == 0 ) 
              cout << " > max new memory alloc size [GB]: " << double(max_coeff_size)*24.0/(1024.0*1024.0*1024.0) << endl;
            int ii = 0;
            for(map<pair<uint8,uint8>,int>::iterator it = cvec_idx.begin();
                it != cvec_idx.end(); ++it) { 
              it->second = ii;
              coeffVec[ii].rbi = it->first.first;
              coeffVec[ii].rbi_nbr = it->first.second;
              coeffVec[ii].transpose = false;
              FOR_I3 coeffVec[ii].coeff[i] = 0.0;
              ++ii;
            }
          }
        }//iter

      }
      else {

        CERR("unrecognized EF_OPS: " << ef_ops);

      }

    }
          
    int nii_new  = coeffVec.size();

    delete[] cvocv_grad_coeff_tmp;
    delete[] cvocv_i_tmp;
    delete[] cvocv_v_tmp;

    delete[] faocv_i_tmp;
    delete[] faocv_v_tmp;

    delete[] n_fa_tmp;
    delete[] x_fa_tmp;
    delete[] cvofa_tmp;
  
    delete[] rbi_g_tmp; 

    // ==============================================
    // hack -- check serial gradient...
    // this code has NOT been modified to support periodicity,
    // but the final serial check below has, so skip for now.
    // ==============================================
    
    /*
      if (mpi_size == 1) {

      const double exact_grad[3] = { 1.123, 2.134, -1.3423 };
    
      double * phi = new double[ncv_tmp];
      for (int icv = 0; icv < ncv_tmp; ++icv) 
      phi[icv] = DOT_PRODUCT(x_cv_tmp[icv],exact_grad);
    
      double (*grad_phi)[3] = new double[ncv_tmp][3];
      for (int icv = 0; icv < ncv_tmp; ++icv) FOR_I3 grad_phi[icv][i] = 0.0;
    
      for (int ii = 0; ii < nii_new; ++ii) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,coeffVec[ii].rbi);
      //assert(bits == 0);
      assert(rank == mpi_rank);
      int rank_nbr,bits_nbr,index_nbr;
      BitUtils::unpackRankBitsIndex(rank_nbr,bits_nbr,index_nbr,coeffVec[ii].rbi_nbr);
      //assert(bits_nbr == 0);
      assert(rank_nbr == mpi_rank);
      FOR_I3 grad_phi[index][i] += coeffVec[ii].coeff[i]*(phi[index_nbr] - phi[index]);
      }

      for (int ibf = 0; ibf < nbf; ++ibf) {
      const double phi_bf = DOT_PRODUCT(exact_grad,x_bf[ibf]);
      const int icv = cvobf_global[ibf];
      FOR_I3 grad_phi[icv][i] += (phi_bf - phi[icv])*n_bf[ibf][i];
      }
    
      for (int icv = 0; icv < ncv_tmp; ++icv) {
      FOR_I3 grad_phi[icv][i] = grad_phi[icv][i]/(vol_cv[icv]*exact_grad[i]) - 1.0; // to produce 0
      }
    
      MiscUtils::dumpRange(grad_phi,ncv_tmp,"extended grad");
    
      FILE * fp = fopen("grad.dat","w");
      for (int icv = 0; icv < ncv_tmp; ++icv) {
      fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %18.15le %18.15le\n",
      x_cv_tmp[icv][0],x_cv_tmp[icv][1],x_cv_tmp[icv][2],
      grad_phi[icv][0],grad_phi[icv][1],grad_phi[icv][2]);
      }
      fclose(fp);

      // internal only...

      for (int ibf = 0; ibf < nbf; ++ibf) {
      const int icv = cvobf_global[ibf];
      FOR_I3 grad_phi[icv][i] = 0.0;
      }
    
      MiscUtils::dumpRange(grad_phi,ncv_tmp,"extended grad - internal only");
    
      delete[] grad_phi;
      delete[] phi;
    
      }
    */

    // transpose elements in the lower-triangular part of D...

    for (int ii=0; ii < nii_new;  ++ii) {
      int rank,bits,icv;
      BitUtils::unpackRankBitsIndex(rank,bits,icv,coeffVec[ii].rbi);
      int rank_nbr,bits_nbr,icv_nbr;
      BitUtils::unpackRankBitsIndex(rank_nbr,bits_nbr,icv_nbr,coeffVec[ii].rbi_nbr);
      assert((rank != rank_nbr)||(icv != icv_nbr)); // you will hit this for thin periodic domains
      if ((rank > rank_nbr)||((rank == rank_nbr)&&(icv > icv_nbr))) {
	const uint8 rbi = coeffVec[ii].rbi;
	coeffVec[ii].rbi = coeffVec[ii].rbi_nbr;
	coeffVec[ii].rbi_nbr = rbi;
	coeffVec[ii].transpose = true;
      } 
    }

    // now everything is upper triangle. Exchange parallel parts...
  
    FOR_RANK send_count[rank] = 0;
  
    for (int ii=0; ii < nii_new;  ++ii) {
      int rank,bits,icv;
      BitUtils::unpackRankBitsIndex(rank,bits,icv,coeffVec[ii].rbi);
      if (rank != mpi_rank)
	++send_count[rank];
    }
    
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
    send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
    
    assert(send_buf_uint8 == NULL); send_buf_uint8 = new uint8[send_count_sum*2];
    assert(send_buf_double == NULL); send_buf_double = new double[send_count_sum][3];
    
    int ii_new = 0;
    for (int ii=0; ii < nii_new;  ++ii) {
      int rank,bits,icv;
      BitUtils::unpackRankBitsIndex(rank,bits,icv,coeffVec[ii].rbi);
      if (rank != mpi_rank) {
	// here it is possible that the row's bits are non-zero, so everthing needs to
	// be transformed so the row's bits are zero. This allows
	// us to use them to store the transpose boolean...
	if (coeffVec[ii].transpose) 
	  send_buf_uint8[send_disp[rank]*2  ] = BitUtils::packRankBitsIndex(rank,1,icv);
	else
	  send_buf_uint8[send_disp[rank]*2  ] = BitUtils::packRankBitsIndex(rank,0,icv); // can't just use coeffVec[ii].rbi because it may have bits
	// and the column...
	if (bits == 0) {
	  send_buf_uint8[send_disp[rank]*2+1] = coeffVec[ii].rbi_nbr;
	  // and the coeff...
	  FOR_I3 send_buf_double[send_disp[rank]][i] = coeffVec[ii].coeff[i];
	}
	else {
	  int rank_nbr,bits_nbr,icv_nbr;
	  BitUtils::unpackRankBitsIndex(rank_nbr,bits_nbr,icv_nbr,coeffVec[ii].rbi_nbr);
	  const int inv_bits = BitUtils::flipPeriodicBits(bits);
	  send_buf_uint8[send_disp[rank]*2+1] = BitUtils::packRankBitsIndex(rank_nbr,BitUtils::addPeriodicBits(bits_nbr,inv_bits),icv_nbr);
          PeriodicData::periodicRotate(coeffVec[ii].coeff,1,inv_bits);
	  FOR_I3 send_buf_double[send_disp[rank]][i] = coeffVec[ii].coeff[i];
	}
	++send_disp[rank];
      }
      else {
	// already on the right rank, but we need the row's bits to be zero...
	if (bits != 0) {
	  int rank_nbr,bits_nbr,icv_nbr;
	  BitUtils::unpackRankBitsIndex(rank_nbr,bits_nbr,icv_nbr,coeffVec[ii].rbi_nbr);
	  const int inv_bits = BitUtils::flipPeriodicBits(bits);
	  coeffVec[ii].rbi = BitUtils::packRankBitsIndex(rank,0,icv);
	  coeffVec[ii].rbi_nbr = BitUtils::packRankBitsIndex(rank_nbr,BitUtils::addPeriodicBits(bits_nbr,inv_bits),icv_nbr);
          PeriodicData::periodicRotate(coeffVec[ii].coeff,1,inv_bits);
	}
	// copy down to compress if required...
	if (ii_new != ii) coeffVec[ii_new].copy(coeffVec[ii]);
	++ii_new;
      }
    }

    // rewind send_disp...
  
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
    
    // setup recv side...
  
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    
    coeffVec.resize(ii_new+recv_count_sum);
    
    FOR_RANK {
      send_count[rank] *= 2;
      send_disp[rank] *= 2;
      recv_count[rank] *= 2;
      recv_disp[rank] *= 2;
    }
    
    assert(recv_buf_uint8 == NULL); recv_buf_uint8 = new uint8[recv_count_sum*2];
    MPI_Alltoallv(send_buf_uint8,send_count,send_disp,MPI_UINT8,
		  recv_buf_uint8,recv_count,recv_disp,MPI_UINT8,mpi_comm);
    delete[] send_buf_uint8; 
  
    FOR_RANK {
      send_count[rank] = (send_count[rank]/2)*3;
      send_disp[rank]  = (send_disp[rank]/2)*3;
      recv_count[rank] = (recv_count[rank]/2)*3;
      recv_disp[rank]  = (recv_disp[rank]/2)*3;
    }
  
    assert(recv_buf_double == NULL); recv_buf_double = new double[recv_count_sum][3];
    MPI_Alltoallv((double*)send_buf_double,send_count,send_disp,MPI_DOUBLE,
		  (double*)recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
    delete[] send_buf_double;
  
    // unpack...
    
    for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
      int rank,bits,icv;
      BitUtils::unpackRankBitsIndex(rank,bits,icv,recv_buf_uint8[irecv*2]);
      assert(rank == mpi_rank);
      if (bits == 1) {
	coeffVec[ii_new].rbi = BitUtils::packRankBitsIndex(rank,0,icv); 
	coeffVec[ii_new].transpose = true;
      }
      else {
	assert(bits == 0);
	coeffVec[ii_new].rbi = recv_buf_uint8[irecv*2];
	coeffVec[ii_new].transpose = false;
      }
      // and the nbr...
      coeffVec[ii_new].rbi_nbr = recv_buf_uint8[irecv*2+1];
      // and the coeff...
      FOR_I3 coeffVec[ii_new].coeff[i] = recv_buf_double[irecv][i];
      ++ii_new;
    }

    assert(ii_new == coeffVec.size());
    delete[] recv_buf_double;
    delete[] recv_buf_uint8;
  
    delete[] recv_disp;
    delete[] recv_count;
    delete[] send_disp;
    delete[] send_count;
  
    // final sort...
  
    sort(coeffVec.begin(),coeffVec.end());
  
    // now build HO faces [0:nef)...
    // first count...

    int nef_i = 0;
    int nef_ip = 0;    
    uint8 rbi_current = BitUtils::packRankBitsIndex(mpi_rank,0,ncv_tmp); // does not exist
    uint8 rbi_nbr_current = rbi_current;
    int icv = -1;
    for (int ii=0,nii=coeffVec.size(); ii < nii;  ++ii) {
      if (coeffVec[ii].rbi != rbi_current) {
	rbi_current = coeffVec[ii].rbi;
	int rank,bits;
	BitUtils::unpackRankBitsIndex(rank,bits,icv,coeffVec[ii].rbi);
	assert(rank == mpi_rank);
	assert(bits == 0);
	rbi_nbr_current = BitUtils::packRankBitsIndex(mpi_rank,0,ncv_tmp); // does not exist
      }
      if (coeffVec[ii].rbi != coeffVec[ii].rbi_nbr) {
	if (coeffVec[ii].rbi_nbr != rbi_nbr_current) {
	  rbi_nbr_current = coeffVec[ii].rbi_nbr;
	  int rank_nbr,bits_nbr,icv_nbr;
	  BitUtils::unpackRankBitsIndex(rank_nbr,bits_nbr,icv_nbr,coeffVec[ii].rbi_nbr);
	  if (rank_nbr == mpi_rank) {
	    // local...
	    // only count a new internal face when its index is larger than ours...
	    assert(icv_nbr != icv);
	    if (icv_nbr > icv) {
	      ++nef_i;
	    }
	    else {
	      assert(0); // if you get here, there is a logic problem
	    }
	  }
	  else if (mpi_rank < rank_nbr) {
	    ++nef_ip;
	  }
	  else {
	    assert(0); // see comment above
	  }
	}
      }
    }
  
    assert(nef == 0);
    nef = nef_i+nef_ip;
    
    assert(cvoef_global == NULL); cvoef_global = new int8[nef][2];
    assert(n_ef == NULL); n_ef = new double[nef][3];
    assert(c_ef == NULL); c_ef = new double[nef][3];

    for (int ief = 0; ief < nef; ++ief) {
      FOR_I2 cvoef_global[ief][i] = -1;
      FOR_I3 n_ef[ief][i] = 0.0;
      FOR_I3 c_ef[ief][i] = 0.0;
    }

    nef_ip = nef_i;
    nef_i = 0;
    int ief = -1;
    rbi_current = BitUtils::packRankBitsIndex(mpi_rank,0,ncv_tmp); // does not exist
    rbi_nbr_current = rbi_current;
    for (int ii=0,nii=coeffVec.size(); ii < nii;  ++ii) {
      if (coeffVec[ii].rbi != rbi_current) {
	rbi_current = coeffVec[ii].rbi;
	int rank,bits;
	BitUtils::unpackRankBitsIndex(rank,bits,icv,coeffVec[ii].rbi);
	assert(rank == mpi_rank);
	assert(bits == 0);
	rbi_nbr_current = BitUtils::packRankBitsIndex(mpi_rank,0,ncv_tmp); // does not exist
      }
      if (coeffVec[ii].rbi != coeffVec[ii].rbi_nbr) {
	// off-diagonal coeff...
	if (coeffVec[ii].rbi_nbr != rbi_nbr_current) {
	  rbi_nbr_current = coeffVec[ii].rbi_nbr;
	  int rank_nbr,bits_nbr,icv_nbr;
	  BitUtils::unpackRankBitsIndex(rank_nbr,bits_nbr,icv_nbr,coeffVec[ii].rbi_nbr);
	  if (rank_nbr == mpi_rank) {
	    // only count a new internal face when its index is larger than ours...
	    assert(icv_nbr != icv);
	    if (icv_nbr > icv) {
	      // this is a new face...
	      ief = nef_i++;
	      assert(cvoef_global[ief][0] == -1);
	      assert(cvoef_global[ief][1] == -1);
	      cvoef_global[ief][0] = cvora[mpi_rank] + icv;
	      cvoef_global[ief][1] = (cvora[mpi_rank] + icv_nbr) | (int8(bits_nbr)<<52);
	    }
	    else {
	      assert(0);
	      ief = -1;
	    }
	  }
	  else if (mpi_rank < rank_nbr) {
	    // this is always a new outward-pointing face...
	    ief = nef_ip++;
	    cvoef_global[ief][0] = cvora[mpi_rank] + icv;
	    cvoef_global[ief][1] = (cvora[rank_nbr] + icv_nbr) | (int8(bits_nbr)<<52);
	  }
	  else {
	    assert(0);
	    ief = -1;
	  }
	}
      
	// ---------------------------------------------------------------
	// we should have ief and fa2_sign set...
	// 0.5*(D-DT) goes in unit_n (area-weighted normal for now)...
	// ---------------------------------------------------------------
	if (ief >= 0) {
	  if (coeffVec[ii].transpose) {
	    FOR_I3 n_ef[ief][i] -= coeffVec[ii].coeff[i];
	  }
	  else {
	    FOR_I3 n_ef[ief][i] += coeffVec[ii].coeff[i];
	  }
	  // 0.5*(D+DT) goes in c[i] 
	  FOR_I3 c_ef[ief][i] += coeffVec[ii].coeff[i];
	}
      
      }
    
    }

    if (checkParam("ZERO_C_EF")) {
      if (mpi_rank == 0) cout << "got ZERO_C_EF: setting c_ef to zero!" << endl;
      FOR_IEF FOR_I3 c_ef[ief][i] = 0.0;
    }

    assert(nef_ip == nef);
  
    // extended face grad check...
  
    if (mpi_size == 1) {

      const double exact_grad[3] = { 1.123, 2.134, -1.3423 };
      double (*grad_phi)[3] = new double[ncv_tmp][3];
      for (int icv = 0; icv < ncv_tmp; ++icv) FOR_I3 grad_phi[icv][i] = 0.0;
    
      //int icv_debug = getIntParam("ICV_DEBUG");

      assert(nef_i == nef);
      FOR_INTERNAL_IEF {

	//if ((fabs(n_ef[ief][2]) > 1.0E-15)||(fabs(c_ef[ief][2]) > 1.0E-15)) {
	//  cout << "Z part of faces: " << n_ef[ief][2] << " " << c_ef[ief][2] << endl;
	//  getchar();
	//}
	
	const int icv0 = cvoef_global[ief][0]; assert((icv0 >= 0)&&(icv0 < ncv_tmp));
	const int icv1 = (cvoef_global[ief][1]&MASK_52BITS); assert((icv1 >= 0)&&(icv1 < ncv_tmp));
	const int bits = (cvoef_global[ief][1]>>52);
	//double flux[3];	FOR_I3 flux[i] = 0.5*n_ef[ief][i]*(phi[icv1]-phi[icv0]) + 0.5*c_ef[ief][i]*(phi[icv1]-phi[icv0]);
	// write this flux in difference-form by using: phi[icv]*(sum(n_ef)+sum(n_bf)) == 0
	// icv0...
	{
	  const double phi_icv0 = DOT_PRODUCT(exact_grad,x_cv[icv0]);
	  double x_cv_icv1_t[3]; FOR_I3 x_cv_icv1_t[i] = x_cv[icv1][i]; PeriodicData::periodicTranslate(x_cv_icv1_t,1,bits);
	  const double phi_icv1_t = DOT_PRODUCT(exact_grad,x_cv_icv1_t);
	  FOR_I3 grad_phi[icv0][i] += 0.5*n_ef[ief][i]*(phi_icv1_t-phi_icv0) + 0.5*c_ef[ief][i]*(phi_icv1_t-phi_icv0);
	}
	// icv1...
	{
	  double x_cv_icv0_t[3]; FOR_I3 x_cv_icv0_t[i] = x_cv[icv0][i]; PeriodicData::periodicTranslate(x_cv_icv0_t,1,BitUtils::flipPeriodicBits(bits));
	  const double phi_icv0_t = DOT_PRODUCT(exact_grad,x_cv_icv0_t);
	  const double phi_icv1 = DOT_PRODUCT(exact_grad,x_cv[icv1]);
	  FOR_I3 grad_phi[icv1][i] -= 0.5*n_ef[ief][i]*(phi_icv0_t-phi_icv1) + 0.5*c_ef[ief][i]*(phi_icv1-phi_icv0_t);
	}
      }
      for (int ibf = 0; ibf < nbf; ++ibf) {
	const double phi_bf = DOT_PRODUCT(exact_grad,x_bf[ibf]);
	const int icv = cvobf_global[ibf]; assert((icv >= 0)&&(icv < ncv_tmp));


	const double phi_icv = DOT_PRODUCT(exact_grad,x_cv[icv]);
	FOR_I3 grad_phi[icv][i] += (phi_bf - phi_icv)*n_bf[ibf][i]; // see note above

	/*
          if (icv == icv_debug) {
	  cout << "GOT icv_debug: " << icv << " x_cv: " << COUT_VEC(x_cv[icv]) << " x_bf: " << COUT_VEC(x_bf[ibf]) << " n_bf: " << COUT_VEC(n_bf[ibf]) << endl;
          }
	*/
	
      }
    
      for (int icv = 0; icv < ncv_tmp; ++icv) {
	FOR_I3 grad_phi[icv][i] = grad_phi[icv][i]/(vol_cv[icv]*exact_grad[i]) - 1.0; // to produce 0
        //if (fabs(grad_phi[icv][2]) > 1.0E-12) cout << "Z grad error: icv: " << icv << " ncv: " << ncv << " err: " << grad_phi[icv][2] << endl;
      }
    
      MiscUtils::dumpRange(grad_phi,ncv_tmp,"extended grad");

      /*
	FILE * fp = fopen("grad.dat","w");
	for (int icv = 0; icv < ncv_tmp; ++icv) {
	fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %18.15le %18.15le\n",
	x_cv_tmp[icv][0],x_cv_tmp[icv][1],x_cv_tmp[icv][2],
	grad_phi[icv][0],grad_phi[icv][1],grad_phi[icv][2]);
	}
	fclose(fp);
      */

      // internal only...

      for (int ibf = 0; ibf < nbf; ++ibf) {
	const int icv = cvobf_global[ibf];
	FOR_I3 grad_phi[icv][i] = 0.0;
      }
    
      MiscUtils::dumpRange(grad_phi,ncv_tmp,"extended grad - internal only");
    
      delete[] grad_phi;
      
    }

    delete[] x_cv_tmp;
  
    assert(efora == NULL);
    MiscUtils::buildXora(efora,nef); assert(efora[mpi_rank+1]-efora[mpi_rank] == nef);
    nef_global = efora[mpi_size];
    
  }

  int read_sles(const string& sles_filename) {
    
    COUT1("StripedMesh::read_sles(\"" << sles_filename << "\")...");

    bool b_foundHash = false;

    // registered sles data...

    MPI_File sles_fh; 
    char dummy[128];
    sprintf(dummy,"%s",sles_filename.c_str());
    int ierr = MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&sles_fh);
    
    if (ierr != 0) { 
      if (mpi_rank == 0) cout << "Error: cannot open sles file: " << sles_filename << endl;
      return -1;
    } 

    int itmp[2] = { 0, 0 };
    if (mpi_rank == 0) { 
      // first 2 ints are: 0. magic number, 1. io version
      MPI_File_read(sles_fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
    }
    MPI_Bcast(itmp,2,MPI_INT,0,mpi_comm); 
    
    bool byte_swap = false;
    //cout << UGP_IO_MAGIC_NUMBER << " " << itmp[0] << endl;
    int itmp_swap[2] = {itmp[0],itmp[1]};
    ByteSwap::byteSwap(itmp_swap,2);
    if ((itmp[1] == 3)||(itmp_swap[1] == 3)) {
      // had to do this so vida_vor solutions could be read in...
      if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
        ByteSwap::byteSwap(itmp,2);
        if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
          if (mpi_rank == 0) cout << "Error: les file does not start as expected." << endl;
          return -1;
        }
        if (mpi_rank == 0) cout << " > file requires byte swapping." << endl;
        byte_swap = true;
      }
    }
    else {
      if (itmp[0] != UGP_IO_MAGIC_NUMBER+1) {
        ByteSwap::byteSwap(itmp,2);
        if (itmp[0] != UGP_IO_MAGIC_NUMBER+1) {
          if (mpi_rank == 0) cout << "Error: sles file does not start as expected." << endl;
          return -1;
        }
        if (mpi_rank == 0) cout << " > file requires byte swapping." << endl;
        byte_swap = true;
      }
    }
    
    const int io_version = itmp[1];
    if (!((io_version == 3)||(io_version >= 5))) {
      if (mpi_rank == 0) cout << "Warning: sles file version not = 3 or >= 5: " << io_version << endl;
    }

    Header header;
    MPI_Offset offset = 8; // 2 ints
    int done = 0;
    while (done != 1) {
      
      if (mpi_rank == 0) {
	MPI_File_read_at(sles_fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE); 
	if (byte_swap) ByteSwap::byteSwapHeader(&header,1);
      }
      MPI_Bcast(&header,1,MPI_Header,0,mpi_comm); 
      
      switch (header.id) {

      case UGP_IO_HASHIDS:

        b_foundHash = true;
        int slesNonMatchFlag;
        //read two hash id's from sles file. First id identifies
        //the current file, store in slesHash.  Second  id
        //identifies the "parent" mles file.  Store in mlesHash
        if (mpi_rank==0)  cout << "RestartHashUtilities::slesReadHashes()" << endl;
        RestartHashUtilities::slesReadHashes(sles_fh, offset, header); //will bCast to all ranks
        if (mpi_rank==0){
          if (io_version == 3) {
            if (!(RestartHashUtilities::mlesHash == RestartHashUtilities::slesHash)) {
              slesNonMatchFlag = 1;
            }
            else {
              cout << " > Found v3 les data with hash: " << RestartHashUtilities::slesHash << endl;
              RestartHashUtilities::myHash.clear();
              slesNonMatchFlag = 0;
            }
          }
          else {
            slesNonMatchFlag = RestartHashUtilities::slesConsistencyCheck(); //-1 not enough info, 0 match, 1 no match
            cout << " > Found sles hash: " << RestartHashUtilities::slesHash << endl;
            if (slesNonMatchFlag==0){
              cout << " >  with mles hash: " << RestartHashUtilities::mlesHash << endl;
            }
            else{
              cout << " > sles expects mles hash: " << RestartHashUtilities::myHash << endl;
              cout << " >      current mles hash: " << RestartHashUtilities::mlesHash << endl;
            }
          }
        }
        MPI_Bcast(&slesNonMatchFlag,1,MPI_INT,0,mpi_comm);
        if (slesNonMatchFlag>0){
          if (io_version == 3) {
            CERR("v3 support in charles requires the same les file to be used for mles (mesh) and sles (data).");
          }
          else {
            CERR("sles data file is not a match for the existing mesh.\n" <<
                 "Consider using INTERP_FROM to initialize with this data.");
          }
        }
        break;
	
      case UGP_IO_I0:

        {
          map<const string,CtiRegister::CtiData>::iterator iter;
          if (io_version == 3) 
            iter = CtiRegister::find_ignore_case(CtiRegister::registeredDataMap,header.name);
          else
            iter = CtiRegister::registeredDataMap.find(header.name);

          if (iter != CtiRegister::registeredDataMap.end()) {
            
            if (mpi_rank == 0) cout << " > reading  I0    \"" << iter->first << "\" value: " << header.idata[0] << endl;
            
            iList.push_back(pair<CtiRegister::CtiData*,int>(&(iter->second),header.idata[0]));
            
          }
          else {

            if (mpi_rank == 0) cout << " > skipping I0    \"" << header.name << "\" value: " << header.idata[0] << endl;
            
          }
        }
	break;
	
      case UGP_IO_D0:
	
        {
          map<const string,CtiRegister::CtiData>::iterator iter;
          if (io_version == 3) 
            iter = CtiRegister::find_ignore_case(CtiRegister::registeredDataMap,header.name);
          else
            iter = CtiRegister::registeredDataMap.find(header.name);

          if (iter != CtiRegister::registeredDataMap.end()) {
            
            if (mpi_rank == 0) cout << " > reading  D0    \"" << iter->first << "\" value: " << header.rdata[0] << endl;
            
            dList.push_back(pair<CtiRegister::CtiData*,double>(&(iter->second),header.rdata[0]));
            
          }
          else {
            
            if (mpi_rank == 0) cout << " > skipping D0    \"" << header.name << "\" value: " << header.rdata[0] << endl;
            
          }
        }
	break;
	
      case UGP_IO_CV_I1:
	
        {
          map<const string,CtiRegister::CtiData>::iterator iter = CtiRegister::registeredDataMap.find(header.name);

          if (iter != CtiRegister::registeredDataMap.end()) {
            
            assert(cvora != NULL);
            
            if (mpi_rank == 0) cout << " > reading  CV_I1 \"" << header.name << "\"";
            
            inList.push_back(InData(header.name));
            assert(inList.back().data == NULL);
            inList.back().data = new int[ncv];
            inList.back().ctiData = &(iter->second);
            
            readChunkedData<int>(sles_fh,offset+header_size+cvora[mpi_rank]*sizeof(int),inList.back().data,ncv,byte_swap,mpi_comm);

            MiscUtils::dumpRange(inList.back().data,ncv,header.name);

          }
          else {

            if (mpi_rank == 0) cout << " > skipping CV_I1 \"" << header.name << "\"" << endl;
            
          }
        }
	break;

      case UGP_IO_FA_I1:
	
        {
          map<const string,CtiRegister::CtiData>::iterator iter = CtiRegister::registeredDataMap.find(header.name);
          if (iter != CtiRegister::registeredDataMap.end()) {
            
            assert(faora != NULL);
            
            if (mpi_rank == 0) cout << " > reading  FA_I1 \"" << header.name << "\"";
            
            inList.push_back(InData(header.name));
            assert(inList.back().data == NULL);
            inList.back().data = new int[nfa];
            inList.back().ctiData = &(iter->second);
            
            readChunkedData<int>(sles_fh,offset+header_size+cvora[mpi_rank]*sizeof(int),inList.back().data,nfa,byte_swap,mpi_comm);

            MiscUtils::dumpRange(inList.back().data,nfa,header.name);

          }
          else {

            if (mpi_rank == 0) cout << " > skipping FA_I1 \"" << header.name << "\"" << endl;
            
          }
        }
	break;

      case UGP_IO_BF_I1:
	
        {
          map<const string,CtiRegister::CtiData>::iterator iter = CtiRegister::registeredDataMap.find(header.name);
          if (iter != CtiRegister::registeredDataMap.end()) {
            
            assert(bfora != NULL);
            
            if (mpi_rank == 0) cout << " > reading  BF_I1 \"" << header.name << "\"";
            
            inList.push_back(InData(header.name));
            assert(inList.back().data == NULL);
            inList.back().data = new int[nbf];
            inList.back().ctiData = &(iter->second);
            
            readChunkedData<int>(sles_fh,offset+header_size+cvora[mpi_rank]*sizeof(int),inList.back().data,nbf,byte_swap,mpi_comm);

            MiscUtils::dumpRange(inList.back().data,nbf,header.name);

          }
          else {

            if (mpi_rank == 0) cout << " > skipping BF_I1 \"" << header.name << "\"" << endl;
            
          }
        }
	break;

      case UGP_IO_CV_D1:
	
	if (strcmp(header.name,"r_vv") == 0) {

	  if (mpi_rank == 0) cout << " > skipping r_vv because it was read in from mles" << endl;

        }
        else {

          map<const string,CtiRegister::CtiData>::iterator iter;
          if (io_version == 3) 
            iter = CtiRegister::find_ignore_case(CtiRegister::registeredDataMap,header.name);
          else
            iter = CtiRegister::registeredDataMap.find(header.name);

          if (iter != CtiRegister::registeredDataMap.end()) {
            
            assert(cvora != NULL);
            
            if (mpi_rank == 0) cout << " > reading CV_D1 \"" << iter->first << "\"";
            
            dnList.push_back(DnData(iter->first));
            assert(dnList.back().data == NULL);
            dnList.back().data = new double[ncv];
            dnList.back().ctiData = &(iter->second);
            
            readChunkedData<double>(sles_fh,offset+header_size+cvora[mpi_rank]*double_size,dnList.back().data,ncv,byte_swap,mpi_comm);

            MiscUtils::dumpRange(dnList.back().data,ncv,iter->first);

          }
          else {

            if (mpi_rank == 0) cout << " > skipping CV_D1 \"" << header.name << "\"" << endl;
            
          }
        }
	break;

      case UGP_IO_FA_D1:
	
        {
          map<const string,CtiRegister::CtiData>::iterator iter;
          if (io_version == 3) 
            iter = CtiRegister::find_ignore_case(CtiRegister::registeredDataMap,header.name);
          else
            iter = CtiRegister::registeredDataMap.find(header.name);

          if (iter != CtiRegister::registeredDataMap.end()) {
            
            assert(faora != NULL);
            
            if (mpi_rank == 0) cout << " > reading  FA_D1 \"" << iter->first << "\"";
            
            dnList.push_back(DnData(iter->first));
            assert(dnList.back().data == NULL);
            dnList.back().data = new double[nfa];
            dnList.back().ctiData = &(iter->second);

            if (io_version == 3) {

              assert(fabfora);
              assert(zone_fabf);
              assert(my_nfa_disp >= 0);
              assert(my_nbf_disp);

              double *data_fabf = new double[nfabf];
              readChunkedData<double>(sles_fh,offset+header_size+fabfora[mpi_rank]*double_size,
                                      data_fabf,nfabf,byte_swap,mpi_comm);
              
              int *send_count = new int[mpi_size];
              int *send_disp  = new int[mpi_size];
              FOR_RANK send_count[rank] = 0;
              const int nbf_tmp = zone_ifabf_vec.size();
              const int nfa_tmp = nfabf-nbf_tmp;
              int * send_buf_int = new int[nfa_tmp];
              double * send_buf_double = new double[nfa_tmp];
              const int nzone = bfZoneVec.size();
              for (int iter2 = 0; iter2 < 2; ++iter2) {
                int ii = 0;
                for (int ifabf = 0; ifabf < nfabf; ++ifabf) {
                  if (zone_fabf[ifabf] == nzone) {
                    const int8 ifa_global = ii+my_nfa_disp; ++ii;
                    const int rank = MiscUtils::getRankInXora(ifa_global,faora);
                    assert((ifa_global >= 0)&&(ifa_global < nfa_global));
                    if (iter2 == 0) {
                      ++send_count[rank];
                    }
                    else {
                      send_buf_int[send_disp[rank]] = ifa_global-faora[rank]; // ifa
                      send_buf_double[send_disp[rank]] = data_fabf[ifabf];
                      ++send_disp[rank];
                    }
                  }
                }

                send_disp[0] = 0;
                for (int rank = 1; rank < mpi_size; ++rank)
                  send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
                assert((send_count[mpi_size-1] + send_disp[mpi_size-1]) == nfa_tmp);
              }

              //  exchange...
              
              int *recv_count = new int[mpi_size];
              MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
              
              int *recv_disp = new int[mpi_size];
              recv_disp[0] = 0;
              for (int rank = 1; rank < mpi_size; ++rank)
                recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
              assert((recv_disp[mpi_size-1] + recv_count[mpi_size-1]) == nfa);
              
              int* recv_buf_int = new int[nfa];
              MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
                            recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
              delete[] send_buf_int; send_buf_int = NULL; 
              double* recv_buf_double = new double[nfa];
              MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
                            recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
              delete[] send_buf_double; send_buf_double = NULL; 

              for (int irecv = 0; irecv < nfa; ++irecv) {
                const int ifa = recv_buf_int[irecv];
                dnList.back().data[ifa] = recv_buf_double[irecv];
              }
              delete[] recv_buf_int; recv_buf_int = NULL;
              delete[] recv_buf_double; recv_buf_double = NULL;

              MiscUtils::dumpRange(dnList.back().data,nfa,iter->first);

              // ibf...

              // find if any zones have this fa record...
              bool found = false;
              int *dnList_index = new int[nzone];
              for (int izone = 0; izone < nzone; ++izone) dnList_index[izone] = -1;
              FOR_IZONE(bfZoneVec) {
                map<const string,CtiRegister::CtiData>::iterator iter2 = CtiRegister::registeredDataMap.find(bfZoneVec[izone].name+":"+iter->first);
                if (iter2 != CtiRegister::registeredDataMap.end()) {
                  found = true;
                  dnList_index[izone] = dnList.size();

                  dnList.push_back(DnData(iter2->first));
                  assert(dnList.back().data == NULL);
              
                  assert(bfZoneVec[izone].bfora != NULL);
                  dnList.back().data = new double[bfZoneVec[izone].nbf];
                  dnList.back().ctiData = &(iter2->second);

                }
              }
              if (found) {
                
                FOR_RANK send_count[rank] = 0;
                assert(send_buf_int == NULL); send_buf_int = new int[nbf_tmp];
                assert(send_buf_double == NULL); send_buf_double = new double[nbf_tmp];
                for (int iter2 = 0; iter2 < 2; ++iter2) {
                  int disp = 0;
                  int izone0 = -1;
                  for (int ii = 0; ii < nbf_tmp; ++ii) {
                    const int izone = zone_ifabf_vec[ii].first;
                    if (dnList_index[izone] != -1) { 
                      if (izone != izone0) {
                        disp = ii;
                        izone0 = izone;
                      }
                      const int ifabf = zone_ifabf_vec[ii].second;
                      assert(zone_fabf[ifabf] == izone);
                      //const int8 ibf_global = ii-disp+(my_nbf_disp[izone]+bfZoneVec[izone].ibf_f_global);
                      //assert((ibf_global >= 0)&&(ibf_global < nbf_global));
                      //const int rank = MiscUtils::getRankInXora(ibf_global,bfora);
                      const int8 ibf_global = ii-disp+my_nbf_disp[izone];
                      assert((ibf_global >= 0)&&(ibf_global < bfZoneVec[izone].nbf_global));
                      const int rank = MiscUtils::getRankInXora(ibf_global,bfZoneVec[izone].bfora);
                      if (iter2 == 0) {
                        ++send_count[rank];
                      }
                      else {
                        send_buf_int[send_disp[rank]*2  ] = izone;
                        send_buf_int[send_disp[rank]*2+1] = ibf_global-bfZoneVec[izone].bfora[rank];
                        send_buf_double[send_disp[rank]] = data_fabf[ifabf];
                        ++send_disp[rank];
                      }
                    }
                  }

                  send_disp[0] = 0;
                  for (int rank = 1; rank < mpi_size; ++rank)
                    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
                  assert((send_count[mpi_size-1] + send_disp[mpi_size-1]) <= nbf_tmp);
                }
                  
                //  exchange...
                
                MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
                
                recv_disp[0] = 0;
                for (int rank = 1; rank < mpi_size; ++rank)
                  recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
                const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
                
                assert(recv_buf_double == NULL); recv_buf_double = new double[recv_count_sum];
                MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
                              recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
                delete[] send_buf_double; 

                FOR_RANK {
                  send_count[rank] *= 2;
                  send_disp[rank] *= 2;
                  recv_count[rank] *= 2;
                  recv_disp[rank] *= 2;
                }

                assert(recv_buf_int == NULL); recv_buf_int = new int[2*recv_count_sum];
                MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
                              recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
                delete[] send_buf_int; 

                int izone0 = -1;
                list<DnData>::iterator dnList_iter = dnList.begin();
                for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
                  const int izone = recv_buf_int[2*irecv];
                  if (izone != izone0) {
                    izone0 = izone;
                    dnList_iter = dnList.begin();
                    assert(dnList_index[izone] != -1);
                    advance(dnList_iter,dnList_index[izone]);
                  }
                  const int ibf = recv_buf_int[2*irecv+1];
                  assert((ibf >= 0)&&(ibf < bfZoneVec[izone].nbf));
                  dnList_iter->data[ibf] = recv_buf_double[irecv];
                }
                delete[] recv_buf_int; 
                delete[] recv_buf_double; 

                FOR_IZONE(bfZoneVec) {
                  if (dnList_index[izone] != -1) {
                    dnList_iter = dnList.begin();
                    advance(dnList_iter,dnList_index[izone]);
                    MiscUtils::dumpRange(dnList_iter->data,bfZoneVec[izone].nbf,bfZoneVec[izone].name+":"+iter->first);
                  }
                }

              }
              delete[] dnList_index;
              delete[] data_fabf;
              delete[] send_disp;
              delete[] send_count;
              delete[] recv_disp;
              delete[] recv_count;

            }
            else {
            
              readChunkedData<double>(sles_fh,offset+header_size+faora[mpi_rank]*double_size,dnList.back().data,nfa,byte_swap,mpi_comm);
              
              MiscUtils::dumpRange(dnList.back().data,nfa,iter->first);
            }

          }
          else {

            if (mpi_rank == 0) cout << " > skipping FA_D1 \"" << header.name << "\"" << endl;
            
          }
        }
	break;

      case UGP_IO_BF_D1:
      case UGP_IO_FAZONE_FA_D1:
	
        {

          map<const string,CtiRegister::CtiData>::iterator iter;
          if (io_version == 3) {
            string my_string = header.name;
            my_string = MiscUtils::replaceAll(my_string,"CDELW","cdel_w");
            iter = CtiRegister::find_ignore_case(CtiRegister::registeredDataMap,my_string);
          }
          else {
            iter = CtiRegister::registeredDataMap.find(header.name);
          }

          if (iter != CtiRegister::registeredDataMap.end()) {
            
            if (mpi_rank == 0) cout << " > reading BF_D1 \"" << iter->first << "\"";
            
            dnList.push_back(DnData(iter->first));
            assert(dnList.back().data == NULL);
            dnList.back().ctiData = &(iter->second);

            const string zone_name = CtiRegister::getFirstSpecifier(header.name);
            map<const string,int>::iterator iter2 = bfZoneNameMap.find(zone_name);
            if (iter2 != bfZoneNameMap.end()) {
              const int izone = iter2->second;
              assert(bfZoneVec[izone].bfora != NULL);
              dnList.back().data = new double[bfZoneVec[izone].nbf];
              
              readChunkedData<double>(sles_fh,offset+header_size+bfZoneVec[izone].bfora[mpi_rank]*double_size,dnList.back().data,bfZoneVec[izone].nbf,byte_swap,mpi_comm);

              MiscUtils::dumpRange(dnList.back().data,bfZoneVec[izone].nbf,header.name);
            }
            else {
              // no matching zone name, so assume it encompasses all bf data...
              assert(bfora != NULL);
              dnList.back().data = new double[nbf];
              
              readChunkedData<double>(sles_fh,offset+header_size+bfora[mpi_rank]*double_size,dnList.back().data,nbf,byte_swap,mpi_comm);

              MiscUtils::dumpRange(dnList.back().data,nbf,iter->first);
            }

          }
          else {

            if (mpi_rank == 0) cout << " > skipping BF_D1 \"" << header.name << "\"" << endl;
            
          }
        }
	break;
	
      case UGP_IO_CV_D2:

        if ((strcmp(header.name,"x_vv") == 0)||(strcmp(header.name,"X_VV") == 0)) {

          if (mpi_rank == 0) cout << " > skipping x_vv because it was read in from mles" << endl;

        }
        else {

          map<const string,CtiRegister::CtiData>::iterator iter;
          if (io_version == 3) 
            iter = CtiRegister::find_ignore_case(CtiRegister::registeredDataMap,header.name);
          else
            iter = CtiRegister::registeredDataMap.find(header.name);

          if (iter != CtiRegister::registeredDataMap.end()) {

            assert(cvora != NULL);
            
            if (mpi_rank == 0) cout << " > reading  CV_D2 \"" << iter->first << "\"";
            
            dn3List.push_back(Dn3Data(iter->first));
            assert(dn3List.back().data == NULL);
            dn3List.back().data = new double[ncv][3];
            dn3List.back().ctiData = &(iter->second);

            readChunkedData<double>(sles_fh,offset+header_size+cvora[mpi_rank]*double_size*3,(double*)dn3List.back().data,ncv*3,byte_swap,mpi_comm);

            MiscUtils::dumpRange(dn3List.back().data,ncv,iter->first);
            
          }
          else {

            if (mpi_rank == 0) cout << " > skipping CV_D2 \"" << header.name << "\"" << endl;
            
          }
        }

	break;

      case UGP_IO_FA_D2:

        {
          map<const string,CtiRegister::CtiData>::iterator iter = CtiRegister::registeredDataMap.find(header.name);
          if (iter != CtiRegister::registeredDataMap.end()) {

            assert(faora != NULL);
            
            if (mpi_rank == 0) cout << " > reading FA_D2 \"" << header.name << "\"";
            
            dn3List.push_back(Dn3Data(header.name));
            assert(dn3List.back().data == NULL);
            dn3List.back().data = new double[nfa][3];
            dn3List.back().ctiData = &(iter->second);

            readChunkedData<double>(sles_fh,offset+header_size+cvora[mpi_rank]*double_size*3,(double*)dn3List.back().data,nfa*3,byte_swap,mpi_comm);

            MiscUtils::dumpRange(dn3List.back().data,nfa,header.name);
            
          }
          else {

            if (mpi_rank == 0) cout << " > skipping FA_D2 \"" << header.name << "\"" << endl;
            
          }
        }
	break;

      case UGP_IO_BF_D2:
      case UGP_IO_FAZONE_FA_D2:

        {

          map<const string,CtiRegister::CtiData>::iterator iter;
          if (io_version == 3) {
            string my_string = header.name;
            my_string = MiscUtils::replaceAll(my_string,"RHOUSLIP","rhou_s");
            iter = CtiRegister::find_ignore_case(CtiRegister::registeredDataMap,my_string);
          }
          else
            iter = CtiRegister::registeredDataMap.find(header.name);

          if (iter != CtiRegister::registeredDataMap.end()) {
            
            if (mpi_rank == 0) cout << " > reading  BF_D2 \"" << iter->first << "\"";
            
            dn3List.push_back(Dn3Data(iter->first));
            assert(dn3List.back().data == NULL);
            dn3List.back().ctiData = &(iter->second);

            const string zone_name = CtiRegister::getFirstSpecifier(header.name);
            map<const string,int>::iterator iter2 = bfZoneNameMap.find(zone_name);
            if (iter2 != bfZoneNameMap.end()) {
              const int izone = iter2->second;

              assert(bfZoneVec[izone].bfora != NULL);
              dn3List.back().data = new double[bfZoneVec[izone].nbf][3];
              
              readChunkedData<double>(sles_fh,offset+header_size+bfZoneVec[izone].bfora[mpi_rank]*3*double_size,(double*)dn3List.back().data,3*bfZoneVec[izone].nbf,byte_swap,mpi_comm);
              MiscUtils::dumpRange(dn3List.back().data,bfZoneVec[izone].nbf,header.name);
              break;
            }
            else {
              // no matching zone name, so assume it encompasses all bf data...
              assert(bfora != NULL);
              dn3List.back().data = new double[nbf][3];
              
              readChunkedData<double>(sles_fh,offset+header_size+bfora[mpi_rank]*double_size*3,(double*)dn3List.back().data,3*nbf,byte_swap,mpi_comm);
              
              MiscUtils::dumpRange(dn3List.back().data,nbf,iter->first);
            }

          }
          else {

            if (mpi_rank == 0) cout << " > skipping BF_D2 \"" << header.name << "\"" << endl;
            
          }
        }
	break;

      case UGP_IO_LPOCV_I:

        {
          if (mpi_rank == 0) {
            cout << " > LPOCV_I \"" << header.name << "\"" << endl;
            cout.flush();
          }

          assert(ncv_global+1 == ByteSwap::getInt8FromLswMswPair(header.idata+0));

          const string name = header.name;
          size_t colon_pos = name.find(":");
          assert(colon_pos != string::npos); // names must be prepended with particle class name
          const string lp_name  = name.substr(0,colon_pos);

          // add lp helper vec entry...
         
          map<const string,int>::iterator iter2 = lpHelperNameMap.find(lp_name);
          assert(iter2 == lpHelperNameMap.end());
          lpHelperVec.push_back(LpHelper());
          lpHelperVec.back().index = lpHelperVec.size()-1;
          lpHelperVec.back().name  = lp_name;
          lpHelperNameMap[lp_name] = lpHelperVec.back().index;
          lpHelperVec.back().nlp_global = ByteSwap::getInt8FromLswMswPair(header.idata+2);
          lpHelperVec.back().lpocv_i_global = new int8[ncv+1];
          readChunkedData<int8>(sles_fh,offset+header_size+cvora[mpi_rank]*int8_size,lpHelperVec.back().lpocv_i_global,ncv+1,byte_swap,mpi_comm);
          lpHelperVec.back().nlp = lpHelperVec.back().lpocv_i_global[ncv]-lpHelperVec.back().lpocv_i_global[0];
        }
        break;

      case UGP_IO_EOF:
        
	done = 1;
        break;
	
	/*
	  default:
	  if (mpi_rank == 0) cout << " > skipping header: " << header.id << " \"" << header.name << "\"." << endl;
	*/

      }
      
      offset += header.skip;
      
    }
    
    MPI_File_close(&sles_fh);
    
    if (!b_foundHash){ //set a default sles hash if one was not found
      std::stringstream ss;
      ss << "Missing sles Hash";
      RestartHashUtilities::slesHash.init(ss,RestartHashUtilities::sha1hashlength);
      if (mpi_rank==0) cout << " > setting default sles hash id " <<  RestartHashUtilities::slesHash << endl;
    }

    // cleanup version 3 stuff...

    DELETE(fabfora);
    DELETE(zone_fabf);
    DELETE(my_nbf_disp);
    
    return 0;
    
  }

  int read_sles_lp(const string& sles_filename) {
    
    COUT1("StripedMesh::read_sles_lp(\"" << sles_filename << "\")...");

    // registered sles data...

    MPI_File sles_fh; 
    char dummy[128];
    sprintf(dummy,"%s",sles_filename.c_str());
    int ierr = MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&sles_fh);
    
    if (ierr != 0) { 
      if (mpi_rank == 0) cout << "Error: cannot open sles file: " << sles_filename << endl;
      return -1;
    } 

    int itmp[2] = { 0, 0 };
    if (mpi_rank == 0) { 
      // first 2 ints are: 0. magic number, 1. io version
      MPI_File_read(sles_fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
    }
    MPI_Bcast(itmp,2,MPI_INT,0,mpi_comm); 
    
    bool byte_swap = false;
    //cout << UGP_IO_MAGIC_NUMBER << " " << itmp[0] << endl;
    if (itmp[0] != UGP_IO_MAGIC_NUMBER+1) {
      ByteSwap::byteSwap(itmp,2);
      if (itmp[0] != UGP_IO_MAGIC_NUMBER+1) {
	if (mpi_rank == 0) cout << "Error: sles file does not start as expected." << endl;
	return -1;
      }
      if (mpi_rank == 0) cout << " > file requires byte swapping." << endl;
      byte_swap = true;
    }
    
    const int io_version = itmp[1];
    if (!(io_version >= 5)) {
      if (mpi_rank == 0) cout << "Warning: sles file version not >= 5: " << io_version << endl;
    }

    Header header;
    MPI_Offset offset = 8; // 2 ints
    int done = 0;
    while (done != 1) {
      
      if (mpi_rank == 0) {
	MPI_File_read_at(sles_fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE); 
	if (byte_swap) ByteSwap::byteSwapHeader(&header,1);
      }
      MPI_Bcast(&header,1,MPI_Header,0,mpi_comm); 
      
      switch (header.id) {

      case UGP_IO_LP_XP:
      case UGP_IO_LP_D2: // effectively the same code

        {
          map<const string,CtiRegister::CtiData>::iterator iter = CtiRegister::registeredDataMap.find(header.name);
          if (iter != CtiRegister::registeredDataMap.end()) {

            if (mpi_rank == 0) {
              if (header.id == UGP_IO_LP_XP) cout << " > LP_XP \"" << header.name << "\"";
              else cout << " > LP_D2 \"" << header.name << "\"";
            }
            
            const string name = header.name;
            size_t colon_pos = name.find(":");
            assert(colon_pos != string::npos); // names must be prepended with particle class name
            const string lp_name  = name.substr(0,colon_pos);

            // find/build lp helper vec entry...
           
            map<const string,int>::iterator iter2 = lpHelperNameMap.find(lp_name);
            assert(iter2 != lpHelperNameMap.end());
            const int lp_index = iter2->second;
            assert((lp_index >= 0)&&(lp_index < lpHelperVec.size()));
            assert(lpHelperVec[lp_index].name  == lp_name);
            assert(lpHelperVec[lp_index].index == lp_index);
            assert(lpHelperVec[lp_index].nlp_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));

            assert(header.idata[2] == 3);

            // use lp helper to read in particle data...

            dn3List.push_back(Dn3Data(header.name));
            assert(dn3List.back().data == NULL);
            dn3List.back().data = new double[lpHelperVec[lp_index].nlp][3];
            dn3List.back().ctiData = &(iter->second);

            readChunkedData<double>(sles_fh,offset+header_size+lpHelperVec[lp_index].lpocv_i_global[0]*double_size*3,(double*)dn3List.back().data,lpHelperVec[lp_index].nlp*3,byte_swap,mpi_comm);
            MiscUtils::dumpRange(dn3List.back().data,lpHelperVec[lp_index].nlp,header.name);
            
          }
          else {

            if (mpi_rank == 0) cout << " > skipping LP_D2 \"" << header.name << "\"" << endl;
            
          }

        }
        break;

      case UGP_IO_LP_D1: 

        {
          map<const string,CtiRegister::CtiData>::iterator iter = CtiRegister::registeredDataMap.find(header.name);
          if (iter != CtiRegister::registeredDataMap.end()) {

            if (mpi_rank == 0) cout << " > LP_D1 \"" << header.name << "\"";
            
            const string name = header.name;
            size_t colon_pos = name.find(":");
            assert(colon_pos != string::npos); // names must be prepended with particle class name
            const string lp_name  = name.substr(0,colon_pos);

            // find lp helper vec entry...
           
            map<const string,int>::iterator iter2 = lpHelperNameMap.find(lp_name);
            assert(iter2 != lpHelperNameMap.end());
            const int lp_index = iter2->second;
            assert((lp_index >= 0)&&(lp_index < lpHelperVec.size()));
            assert(lpHelperVec[lp_index].name  == lp_name);
            assert(lpHelperVec[lp_index].index == lp_index);
            assert(lpHelperVec[lp_index].nlp_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));

            // use lp helper to read in particle data...

            dnList.push_back(DnData(header.name));
            assert(dnList.back().data == NULL);
            dnList.back().data = new double[lpHelperVec[lp_index].nlp];
            dnList.back().ctiData = &(iter->second);

            readChunkedData<double>(sles_fh,offset+header_size+lpHelperVec[lp_index].lpocv_i_global[0]*double_size,dnList.back().data,lpHelperVec[lp_index].nlp,byte_swap,mpi_comm);
            MiscUtils::dumpRange(dnList.back().data,lpHelperVec[lp_index].nlp,header.name);
            
          }
          else {

            if (mpi_rank == 0) cout << " > skipping LP_D1 \"" << header.name << "\"" << endl;
            
          }

        }
        break;

      case UGP_IO_LP_I1: 

        {
          map<const string,CtiRegister::CtiData>::iterator iter = CtiRegister::registeredDataMap.find(header.name);
          if (iter != CtiRegister::registeredDataMap.end()) {

            if (mpi_rank == 0) cout << " > LP_I1 \"" << header.name << "\"";
            
            const string name = header.name;
            size_t colon_pos = name.find(":");
            assert(colon_pos != string::npos); // names must be prepended with particle class name
            const string lp_name  = name.substr(0,colon_pos);

            // find lp helper vec entry...
           
            map<const string,int>::iterator iter2 = lpHelperNameMap.find(lp_name);
            assert(iter2 != lpHelperNameMap.end());
            const int lp_index = iter2->second;
            assert((lp_index >= 0)&&(lp_index < lpHelperVec.size()));
            assert(lpHelperVec[lp_index].name  == lp_name);
            assert(lpHelperVec[lp_index].index == lp_index);
            assert(lpHelperVec[lp_index].nlp_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));

            // use lp helper to read in particle data...

            inList.push_back(InData(header.name));
            assert(inList.back().data == NULL);
            inList.back().data = new int[lpHelperVec[lp_index].nlp];
            inList.back().ctiData = &(iter->second);

            readChunkedData<int>(sles_fh,offset+header_size+lpHelperVec[lp_index].lpocv_i_global[0]*int_size,inList.back().data,lpHelperVec[lp_index].nlp,byte_swap,mpi_comm);
            MiscUtils::dumpRange(inList.back().data,lpHelperVec[lp_index].nlp,header.name);
            
          }
          else {

            if (mpi_rank == 0) cout << " > skipping LP_I1 \"" << header.name << "\"" << endl;
            
          }

        }
        break;

      case UGP_IO_EOF:
        
	done = 1;
        break;
	
	/*
	  default:
	  if (mpi_rank == 0) cout << " > skipping header: " << header.id << " \"" << header.name << "\"." << endl;
	*/

      }
      
      offset += header.skip;
      
    }
    
    MPI_File_close(&sles_fh);
    
    return 0;
    
  }

};

#endif
