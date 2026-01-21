#ifndef _CUTTABLE_VORONOI_DATA_
#define _CUTTABLE_VORONOI_DATA_

#include "CTI.hpp"
using namespace CTI;

#include "MiscUtils.hpp" // required in the cut routines (cut.hpp)
#include <set>
#include <stack>
#include <queue>

class FaceGeometryData {
public:
  double * x0; // pointer to the first node found for this face
  double xc[3];
  double normal[3];
  double area;
  double r2_max;
  FaceGeometryData(double * x0) {
    this->x0 = x0;
    xc[0] = 0.0;
    xc[1] = 0.0;
    xc[2] = 0.0;
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;
    area = 0.0;
    r2_max = 0.0;
  }
  FaceGeometryData() {
    this->x0 = NULL;
  }
  void zero() {
    this->x0 = NULL;
    xc[0] = 0.0;
    xc[1] = 0.0;
    xc[2] = 0.0;
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;
    area = 0.0;
  }
  void complete() {
    // once all the tris have been added in, this turns the data
    // into the correct values...
    FOR_I3 xc[i] /= 3.0*area;
    FOR_I3 normal[i] /= 2.0;
    area /= 2.0;
  }
};

class BfGeometryData {
public:
  double * x0; // pointer to the first node found for this face
  double xc[3];
  double normal[3];
  double area;
  double area_over_delta;
  double Gij[3][3];
  int index1; // can be used any way you want!
  int index2; // can be used any way you want!
  BfGeometryData(double * x0) {
    this->x0 = x0;
    xc[0] = 0.0;
    xc[1] = 0.0;
    xc[2] = 0.0;
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;
    area = 0.0;
    area_over_delta = 0.0;
    FOR_I3 FOR_J3 Gij[i][j] = 0.0;
  }
  BfGeometryData(double * x0,const int index1) {
    this->x0 = x0;
    this->index1 = index1;
    xc[0] = 0.0;
    xc[1] = 0.0;
    xc[2] = 0.0;
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;
    area = 0.0;
    area_over_delta = 0.0;
    FOR_I3 FOR_J3 Gij[i][j] = 0.0;
  }
  BfGeometryData(double * x0,const int index1,const int index2) {
    this->x0 = x0;
    this->index1 = index1;
    this->index2 = index2;
    xc[0] = 0.0;
    xc[1] = 0.0;
    xc[2] = 0.0;
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;
    area = 0.0;
    area_over_delta = 0.0;
    FOR_I3 FOR_J3 Gij[i][j] = 0.0;
  }
  BfGeometryData() {
    this->x0 = NULL;
    xc[0] = 0.0;
    xc[1] = 0.0;
    xc[2] = 0.0;
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;
    area = 0.0;
    area_over_delta = 0.0;
    FOR_I3 FOR_J3 Gij[i][j] = 0.0;
  }
  void zero() {
    this->x0 = NULL;
    xc[0] = 0.0;
    xc[1] = 0.0;
    xc[2] = 0.0;
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;
    area = 0.0;
    area_over_delta = 0.0;
    FOR_I3 FOR_J3 Gij[i][j] = 0.0;
    index1 = -1;
    index2 = -1;
  }
  void complete() {
    // once all the tris have been added in, this turns the data
    // into the correct values...
    FOR_I3 xc[i] /= 3.0*area;
    FOR_I3 normal[i] /= 2.0;
    area /= 2.0;
  }
};

class CutCubeData {
public:
  int id; // cube edge [0:12)
  int sign; // pointing plus or minus
  int ino;
  double d; // distance
  CutCubeData(const int id,const int sign,const int ino,const double d) {
    this->id = id;
    this->sign = sign;
    this->ino = ino;
    this->d = d;
  }
  // for sort...
  bool operator<(const CutCubeData& rhs) const { return (id < rhs.id) || ((id == rhs.id)&&(d < rhs.d)); }
};

#include "IntersectionStuff.hpp"

class CuttableVoronoiData {
private:

  /*
  enum IntersectionKind {
    NODE_INTERSECTION,
    EDGE_INTERSECTION,
    FACE_INTERSECTION,
    NODE_NODE_INTERSECTION,
    NODE_EDGE_INTERSECTION,
    NODE_FACE_INTERSECTION,
    EDGE_EDGE_INTERSECTION,
    EDGE_FACE_INTERSECTION,
  };

  class IntersectionData {
  public:
    IntersectionKind kind;
    int idata[2];
    double ddata[2];
    int ino;
    IntersectionData() { }
    IntersectionData(const IntersectionKind kind,const int idata0,const int idata1) {
      this->kind = kind;
      if (kind == NODE_NODE_INTERSECTION) {
        assert(idata0 != idata1);
        idata[0] = min(idata0,idata1);
        idata[1] = max(idata0,idata1);
      }
      else {
        assert(kind == NODE_FACE_INTERSECTION);
        idata[0] = idata0;
        idata[1] = idata1;
      }
    }
    IntersectionData(const IntersectionKind kind,const int idata0,const int idata1,const double ddata0) {
      this->kind = kind;
      if (kind == NODE_EDGE_INTERSECTION) {
        idata[0] = idata0;
        idata[1] = idata1;
        // the passed weight is for the edge...
        ddata[1] = ddata0;
      }
      else {
        assert(kind == EDGE_FACE_INTERSECTION);
        idata[0] = idata0;
        idata[1] = idata1;
        // the passed weight is for the edge...
        ddata[0] = ddata0;
      }
    }
    IntersectionData(const IntersectionKind kind,const int idata0,const int idata1,const double ddata0,const double ddata1) {
      this->kind = kind;
      assert(kind == EDGE_EDGE_INTERSECTION);
      assert(idata0 != idata1);
      if (idata0 < idata1) {
        idata[0] = idata0;
        idata[1] = idata1;
        ddata[0] = ddata0;
        ddata[1] = ddata1;
      }
      else {
        idata[0] = idata1;
        idata[1] = idata0;
        ddata[0] = ddata1;
        ddata[1] = ddata0;
      }
    }
    void calcXi(double xi[3],const double (* const x_no)[3],const int (* const nooed)[2]) const {
      switch (kind) {
      case NODE_NODE_INTERSECTION:
        FOR_I3 xi[i] = 0.5*(x_no[idata[0]][i]+x_no[idata[1]][i]);
        break;
      case NODE_EDGE_INTERSECTION:
        FOR_I3 xi[i] = 0.5*(x_no[idata[0]][i] + ddata[1]*x_no[nooed[idata[1]][1]][i] + (1.0-ddata[1])*x_no[nooed[idata[1]][0]][i]);
        break;
      case NODE_FACE_INTERSECTION:
        FOR_I3 xi[i] = x_no[idata[0]][i];
        break;
      case EDGE_EDGE_INTERSECTION:
        FOR_I3 xi[i] = 0.5*(ddata[0]*x_no[nooed[idata[0]][1]][i] + (1.0-ddata[0])*x_no[nooed[idata[0]][0]][i] +
                            ddata[1]*x_no[nooed[idata[1]][1]][i] + (1.0-ddata[1])*x_no[nooed[idata[1]][0]][i]);
        break;
      case EDGE_FACE_INTERSECTION:
        FOR_I3 xi[i] = ddata[0]*x_no[nooed[idata[0]][1]][i] + (1.0-ddata[0])*x_no[nooed[idata[0]][0]][i];
        break;
      default:
        assert(0);
      }
    }
    IntersectionData(const IntersectionKind kind0,const int idata0,const double ddata0,
                     const IntersectionKind kind1,const int idata1,const double ddata1) {
      assert(0);
      if ((kind0 == NODE_INTERSECTION)&&(kind1 == NODE_INTERSECTION)) {
        kind = NODE_NODE_INTERSECTION;
        // sort nodes by index...
        if (!(idata0 != idata1)) {
          cout << "IntersectionData() NODE_NODE_INTERSECTION has the same nodes: " << idata0 << endl;
          assert(0);
        }
        assert(idata0 != idata1); // different parts
        if (idata0 < idata1) {
          idata[0] = idata0;
          idata[1] = idata1;
        }
        else {
          idata[0] = idata1;
          idata[1] = idata0;
        }
      }
      else if ((kind0 == NODE_INTERSECTION)&&(kind1 == EDGE_INTERSECTION)) {
        kind = NODE_EDGE_INTERSECTION;
        idata[0] = idata0; // node index
        idata[1] = idata1; // edge index
        ddata[1] = ddata1; // distance along edge
      }
      else if ((kind0 == EDGE_INTERSECTION)&&(kind1 == NODE_INTERSECTION)) {
        kind = NODE_EDGE_INTERSECTION;
        idata[0] = idata1; // node index
        idata[1] = idata0; // edge index
        ddata[1] = ddata0; // distance along edge
      }
      else if ((kind0 == NODE_INTERSECTION)&&(kind1 == FACE_INTERSECTION)) {
        kind = NODE_FACE_INTERSECTION;
        idata[0] = idata0; // node index
        idata[1] = idata1; // face index
      } 
      else if ((kind0 == EDGE_INTERSECTION)&&(kind1 == EDGE_INTERSECTION)) {
        kind = EDGE_EDGE_INTERSECTION;
        // sort edges by index...
        assert(idata0 != idata1); // different parts
        if (idata0 < idata1) {
          idata[0] = idata0;
          ddata[0] = ddata0;
          idata[1] = idata1;
          ddata[1] = ddata1;
        }
        else {
          idata[0] = idata1;
          ddata[0] = ddata1;
          idata[1] = idata0;
          ddata[1] = ddata0;
        }
      }
      else if ((kind0 == EDGE_INTERSECTION)&&(kind1 == FACE_INTERSECTION)) {
        kind = EDGE_FACE_INTERSECTION;
        idata[0] = idata0; // edge index
        ddata[0] = ddata0; // distance along edge
        idata[1] = idata1; // face index
      }
      else {
        assert(0);
      }
    }
    void dump() {
      switch (kind) {
      case NODE_NODE_INTERSECTION:
        cout << "NODE_NODE " << idata[0] << " " << idata[1] << endl;
        break;
      case NODE_EDGE_INTERSECTION:
        cout << "NODE_EDGE " << idata[0] << " " << idata[1] << ":" << ddata[1] << endl;
        break;
      case NODE_FACE_INTERSECTION:
        cout << "NODE_FACE " << idata[0] << " " << idata[1] << endl;
        break;
      case EDGE_EDGE_INTERSECTION:
        cout << "EDGE_EDGE " << idata[0] << ":" << ddata[0] << " " << idata[1] << ":" << ddata[1] << endl;
        break;
      case EDGE_FACE_INTERSECTION:
        cout << "EDGE_FACE " << idata[0] << ":" << ddata[0] << " " << idata[1] << endl;
        break;
      default:
        assert(0);
      }
    }
    inline bool operator<(const IntersectionData& rhs) const { 
      return (kind < rhs.kind) || ( (kind == rhs.kind) && ((idata[0] < rhs.idata[0]) || ((idata[0] == rhs.idata[0]) && (idata[1] < rhs.idata[1]))));
    }
    inline bool operator==(const IntersectionData& rhs) const { 
      return (kind == rhs.kind) && (idata[0] == rhs.idata[0]) && (idata[1] == rhs.idata[1]);
    }
    inline bool operator!=(const IntersectionData& rhs) const { 
      return (kind != rhs.kind) || (idata[0] != rhs.idata[0]) || (idata[1] != rhs.idata[1]);
    }
  };
  */

  int ned_max,nno_max;
  
  int faoce[12][2];
  int cdoce[12];
  int nooce[12][2];

public:

  int ned,nno;

  int (*nooed)[2];
  int (*faoed)[2];
  double (*x_no)[3];

  double d2_max;
  double Lmin[3];
  double Lmax[3];
  
#define X0_FACE -2
#define X1_FACE -3
#define Y0_FACE -4
#define Y1_FACE -5
#define Z0_FACE -6
#define Z1_FACE -7

  CuttableVoronoiData() {
    
    // int faoce[12][2]
    faoce[0][0] = X0_FACE; faoce[0][1] = Y0_FACE;
    faoce[1][0] = Y1_FACE; faoce[1][1] = X0_FACE;
    faoce[2][0] = Z0_FACE; faoce[2][1] = X0_FACE;
    faoce[3][0] = X0_FACE; faoce[3][1] = Z1_FACE;
    faoce[4][0] = Y0_FACE; faoce[4][1] = X1_FACE;
    faoce[5][0] = X1_FACE; faoce[5][1] = Y1_FACE;
    faoce[6][0] = X1_FACE; faoce[6][1] = Z0_FACE;
    faoce[7][0] = Z1_FACE; faoce[7][1] = X1_FACE;
    faoce[8][0] = Y0_FACE; faoce[8][1] = Z0_FACE;
    faoce[9][0] = Z1_FACE; faoce[9][1] = Y0_FACE;
    faoce[10][0] = Z0_FACE; faoce[10][1] = Y1_FACE;
    faoce[11][0] = Y1_FACE; faoce[11][1] = Z1_FACE;
    
    // each cube edge is associated with a Cartesian direction
    // Cartesian-direction of cube-edge...
    // int cdoce[12]
    cdoce[0] = 2;
    cdoce[1] = 2;
    cdoce[2] = 1;
    cdoce[3] = 1;
    cdoce[4] = 2;
    cdoce[5] = 2;
    cdoce[6] = 1;
    cdoce[7] = 1;
    cdoce[8] = 0;
    cdoce[9] = 0;
    cdoce[10] = 0;
    cdoce[11] = 0;
    
    // and corner node-of cube-edge
    // int nooce[12][2]
    nooce[0][0] = 0; nooce[0][1] = 4;
    nooce[1][0] = 2; nooce[1][1] = 6;
    nooce[2][0] = 0; nooce[2][1] = 2;
    nooce[3][0] = 4; nooce[3][1] = 6;
    nooce[4][0] = 1; nooce[4][1] = 5;
    nooce[5][0] = 3; nooce[5][1] = 7;
    nooce[6][0] = 1; nooce[6][1] = 3;
    nooce[7][0] = 5; nooce[7][1] = 7;
    nooce[8][0] = 0; nooce[8][1] = 1;
    nooce[9][0] = 4; nooce[9][1] = 5;
    nooce[10][0] = 2; nooce[10][1] = 3;
    nooce[11][0] = 6; nooce[11][1] = 7;
    
    ned_max = 32;
    nno_max = 32;
    ned = nno = 0;

    nooed = new int[ned_max][2];
    faoed = new int[ned_max][2];
    x_no = new double[nno_max][3];

    d2_max = 0.0;
    FOR_I3 Lmin[i] = 0.0;
    FOR_I3 Lmax[i] = 0.0;

  }

  CuttableVoronoiData(const CuttableVoronoiData& other) {

    assert(0);

    // the copy construtor is called when vector is resized...
    //cout << "CuttableVoronoiData() copy" << endl;

    ned_max = other.ned_max;
    nno_max = other.nno_max;
    ned     = other.ned;
    nno     = other.nno;

    nooed = new int[ned_max][2];
    faoed = new int[ned_max][2];
    x_no = new double[nno_max][3];

    assert(ned == 0);
    assert(nno == 0);

    d2_max = 0.0;
    FOR_I3 Lmin[i] = 0.0;
    FOR_I3 Lmax[i] = 0.0;

  }

  ~CuttableVoronoiData() {

    //cout << "~CuttableVoronoiData()" << endl;

    delete[] nooed;
    delete[] faoed;
    delete[] x_no;

  }

  void writeBinary(const string& filename) const {

    FILE * fp = fopen(filename.c_str(),"wb");
    assert(fp);

    fwrite(&ned,sizeof(int),1,fp);
    fwrite(&nno,sizeof(int),1,fp);
    fwrite(nooed,sizeof(int),ned*2,fp);
    fwrite(faoed,sizeof(int),ned*2,fp);
    fwrite(x_no,sizeof(double),nno*3,fp);

    fclose(fp);

  }

  void readBinary(const string& filename) {

    FILE * fp = fopen(filename.c_str(),"rb");
    if (fp == NULL) return;

    int ned_,nno_;
    fread(&ned_,sizeof(int),1,fp);
    fread(&nno_,sizeof(int),1,fp);

    resize_ned(ned_);
    fread(nooed,sizeof(int),ned_*2,fp);
    fread(faoed,sizeof(int),ned_*2,fp);

    resize_nno(nno_);
    fread(x_no,sizeof(double),nno_*3,fp);

    fclose(fp);

  }

  void addEdgeToFp(FILE * fp,const int ied) const {
    // used for debugging...
    assert((ied >= 0)&&(ied < ned));
    for (int i = 1; i < 10; ++i) {
      // plot 10 points along the edge, not including the end points, and biased slightly
      // towards node0...
      const double wgt = pow(double(i)/10.0,1.2);
      fprintf(fp,"%18.15le %18.15le %18.15le\n",
              (1.0-wgt)*x_no[nooed[ied][0]][0]+wgt*x_no[nooed[ied][1]][0],
              (1.0-wgt)*x_no[nooed[ied][0]][1]+wgt*x_no[nooed[ied][1]][1],
              (1.0-wgt)*x_no[nooed[ied][0]][2]+wgt*x_no[nooed[ied][1]][2]);
    }
  }

  // =======================================================================
  // routines for processing the cvd that may involve intersecting 
  // multiple parts: used in moving solver...
  void debug();
  double getBestXInRange(const double xmin,const double xmax) const;
  double getBestYInRange(const double ymin,const double ymax) const;
  double getBestZInRange(const double zmin,const double zmax) const;
  void doit(const double L,const int * const ist_and_ipart_buf,const bool debug);
  void doit2(const double Lxm,const double Lxp,const double Lym,const double Lyp,const double Lzm,const double Lzp,
             const int * const ist_and_ipart_buf,const bool debug);
  void buildSeed(const double this_delta_cv,const int * const ist_and_ipart_buf,const bool debug);
  int triangulateFaceForDoit(const int ifa,vector<pair<int,int> >& faceEdgeVec,const vector<int>& edofa_i);
  void packTriEdgeVecForDoit(vector<int>& teVec,const int ifa,const vector<pair<int,int> >& faceEdgeVec,const vector<int>& edofa_i,const int ied_f);
  void addEdgeTriIntersectionForDoit(vector<IntersectionData>& localIntersectionVec,
                                     const int ied0,const int nootr1[3],const int edotr1[3],const int faoed1,
                                     const int8 tet_vol0,const int8 tet_vol1,
                                     const int8 tet_vol_no0,const int8 tet_vol_no1,const int8 tet_vol_no2);
  int processCornerStack(queue<pair<int,int> >& cornerStack,int corner_flag[8][3],vector<pair<double,pair<int,int> > > cutCubeNodeVec[12],const int nno_orig);
  void walkAndFlagEdges(int * ed_flag,const int ned_cube,const bool debug) const;
  int getNextEdge(const int ied,const int ifa,const int ino,const multimap<const pair<int,int>,int>& edgeMultimap) const;
  // end of moving solver routines
  // =======================================================================
  
  void clear() {

    ned = nno = 0;
    d2_max = 0.0;
    FOR_I3 Lmin[i] = 0.0;
    FOR_I3 Lmax[i] = 0.0;

  }

  bool empty() const {
    if (nno == 0) {
      assert(ned == 0);
      return true;
    }
    assert(nno > 0);
    assert(ned > 0);
    return false;
  }

  void init(const CuttableVoronoiData& other) {

    resize_ned(other.ned);
    for (int ied = 0; ied < ned; ++ied) {
      FOR_I2 nooed[ied][i] = other.nooed[ied][i];
      FOR_I2 faoed[ied][i] = other.faoed[ied][i];
    }
    resize_nno(other.nno);
    for (int ino = 0; ino < nno; ++ino) {
      FOR_I3 x_no[ino][i] = other.x_no[ino][i];
    }

  }

  int getNed() const { return ned; }

  int getNno() const { return nno; }

  inline int new_node() {

    if ( nno == nno_max ) {
      nno_max *= 2;
      double (*x_no_new)[3] = new double[nno_max][3];
      FOR_INO FOR_I3 x_no_new[ino][i] = x_no[ino][i];
      delete[] x_no;
      x_no = x_no_new;
    }
    return(nno++);

  }

  inline void resize_nno(const int nno_new) {

    if ( nno_new > nno_max ) {
      nno_max = nno_new;
      double (*x_no_new)[3] = new double[nno_max][3];
      FOR_INO FOR_I3 x_no_new[ino][i] = x_no[ino][i];
      delete[] x_no;
      x_no = x_no_new;
    }
    nno = nno_new;

  }

  void translate(const double dx[3]) {
    for (int ino = 0; ino < nno; ++ino) {
      FOR_I3 x_no[ino][i] += dx[i];
    }
  }

  inline int new_edge() {

    if ( ned == ned_max ) {
      ned_max *= 2;
      int (*nooed_new)[2] = new int[ned_max][2];
      int (*faoed_new)[2] = new int[ned_max][2];
      FOR_IED {
	nooed_new[ied][0] = nooed[ied][0];
	nooed_new[ied][1] = nooed[ied][1];
	faoed_new[ied][0] = faoed[ied][0];
	faoed_new[ied][1] = faoed[ied][1];
      }
      delete[] nooed; nooed = nooed_new;
      delete[] faoed; faoed = faoed_new;
    }
    return(ned++);

  }

  inline void resize_ned(const int ned_new) {

    if ( ned_new > ned_max ) {
      ned_max = ned_new;
      int (*nooed_new)[2] = new int[ned_max][2];
      int (*faoed_new)[2] = new int[ned_max][2];
      FOR_IED {
	nooed_new[ied][0] = nooed[ied][0];
	nooed_new[ied][1] = nooed[ied][1];
	faoed_new[ied][0] = faoed[ied][0];
	faoed_new[ied][1] = faoed[ied][1];
      }
      delete[] nooed; nooed = nooed_new;
      delete[] faoed; faoed = faoed_new;
    }
    ned = ned_new;

  }

  bool hasSeedBoundary() const {

    // return true if the volume has any seed boundary...

    FOR_IED {
      FOR_I2 {
	if ((faoed[ied][i] >= -7)&&(faoed[ied][i] <= -2))
	  return(true);
      }
    }
    return(false);

  }

  bool hasSurfaceBoundary() const {

    // return true if the volume has any surface boundary...

    FOR_IED {
      FOR_I2 {
	if (faoed[ied][i] >= 0)
	  return(true);
      }
    }
    return(false);

  }

  void setD2Max() {
    d2_max = 0.0;
    FOR_I3 Lmax[i] = 0.0;
    FOR_I3 Lmin[i] = 0.0;
    FOR_INO {
      const double d2 = DOT_PRODUCT(x_no[ino],x_no[ino]);
      d2_max = max(d2_max,d2);
      FOR_I3 {
	Lmax[i] = max(Lmax[i],sqrt(d2)+x_no[ino][i]);
	Lmin[i] = min(Lmin[i],-sqrt(d2)+x_no[ino][i]);
      }
    }
  }

  double calcVolume() const {
    double volume = 0.0;
    map<const int,const double *> firstNodeMap;
    FOR_IED {
      FOR_I2 {
        const int ifa = faoed[ied][i];
        map<const int,const double *>::iterator iter = firstNodeMap.find(ifa);
        if (iter == firstNodeMap.end()) {
          firstNodeMap[ifa] = x_no[nooed[ied][i]];
        }
        else {
          const double * const x1 = x_no[nooed[ied][i]];
          const double * const x2 = x_no[nooed[ied][1-i]];
          if ((iter->second != x1)&&(iter->second != x2)) {
            const double this_vol = CROSS_DOT(iter->second,x1,x2);
            volume += this_vol;
          }
        }
      }
    }
    return volume/6.0;
  }

  void calcVolumeAndCentroid(double &volume,double centroid[3]) const {
    volume = 0.0;
    FOR_I3 centroid[i] = 0.0;
    map<const int,const double *> firstNodeMap;
    FOR_IED {
      FOR_I2 {
        const int ifa = faoed[ied][i];
        map<const int,const double *>::iterator iter = firstNodeMap.find(ifa);
        if (iter == firstNodeMap.end()) {
          firstNodeMap[ifa] = x_no[nooed[ied][i]];
        }
        else {
          const double * const x1 = x_no[nooed[ied][i]];
          const double * const x2 = x_no[nooed[ied][1-i]];
          if ((iter->second != x1)&&(iter->second != x2)) {
            const double this_vol = CROSS_DOT(iter->second,x1,x2);
            volume += this_vol;
            FOR_I3 centroid[i] += this_vol*(iter->second[i]+x1[i]+x2[i]); // average 4 tet corners: note one corner is (0,0,0)
          }
        }
      }
    }
    if (volume != 0.0) {
      FOR_I3 centroid[i] /= volume*4.0;
      volume /= 6.0;
    }
    else {
      // for zero-volume, just return zero centroid.
      // centroid not defined for zero volume...
      FOR_I3 centroid[i] = 0.0;
    }
  }

  void calcXcvAndStVec(double x_cv[3],vector<int>& stVec) {

    FOR_I3 x_cv[i] = 0.0;
    double wgt = 0.0;

    map<const int,int> faceMap;
    FOR_IED {
      // a positive faoed is an ist reference in the surface...
      if (faoed[ied][0] >= 0) {
	map<const int,int>::iterator iter = faceMap.find(faoed[ied][0]);
	if (iter == faceMap.end()) {
	  faceMap[faoed[ied][0]] = nooed[ied][0];
	  stVec.push_back(faoed[ied][0]);
	}
	else {
	  const int ino = iter->second;
	  if ((nooed[ied][0] != ino)&&(nooed[ied][1] != ino)) {
	    const double n[3] = TRI_NORMAL_2(x_no[ino],x_no[nooed[ied][0]],x_no[nooed[ied][1]]);
	    const double this_wgt = MAG(n);
	    FOR_I3 x_cv[i] += this_wgt*(x_no[ino][i] + x_no[nooed[ied][0]][i] + x_no[nooed[ied][1]][i]);
	    wgt += this_wgt;
	  }
	}
      }
      if (faoed[ied][1] >= 0) {
	map<const int,int>::iterator iter = faceMap.find(faoed[ied][1]);
	if (iter == faceMap.end()) {
	  faceMap[faoed[ied][1]] = nooed[ied][1];
	  stVec.push_back(faoed[ied][1]);
	}
	else {
	  const int ino = iter->second;
	  if ((nooed[ied][0] != ino)&&(nooed[ied][1] != ino)) {
	    const double n[3] = TRI_NORMAL_2(x_no[ino],x_no[nooed[ied][0]],x_no[nooed[ied][1]]);
	    const double this_wgt = MAG(n);
	    FOR_I3 x_cv[i] += this_wgt*(x_no[ino][i] + x_no[nooed[ied][0]][i] + x_no[nooed[ied][1]][i]);
	    wgt += this_wgt;
	  }
	}
      }
    }

    if (wgt != 0.0) {
      FOR_I3 x_cv[i] /= 3.0*wgt;
    }

  }

  void dump() {

    set<int> faSet;
    FOR_IED {
      faSet.insert(faoed[ied][0]);
      faSet.insert(faoed[ied][1]);
    }

    cout << "Cuttable dump: ned: " << ned << " nno: " << nno << " nfa: " << faSet.size() << endl;
    FOR_IED {
      const double length = DIST(x_no[nooed[ied][0]],x_no[nooed[ied][1]]);
      cout << "ied: " << ied << " length: " << length << endl;
    }

  }

  void writeTecplot(const int iv) const {

    double x0[3] = { 0.0, 0.0, 0.0 };
    writeTecplot(iv,x0);

  }

  void writeTecplot(const int iv,const double x0[3]) const {

    cout << "writeTeclot: " << iv << " x0: " << COUT_VEC(x0) << endl;
    
    if (nno == 0)
      return;
    
    char filename[128];
    sprintf(filename,"cvd.%04d.%06d.dat",mpi_rank,iv);

    writeTecplot(filename,x0);

  }

  void writeTecplot(const string& filename,const double x0[3]) const {

    FILE * fp = fopen(filename.c_str(),"w");
    if (fp == NULL) {
      cout << "Error: cannot open filename: " << filename << endl;
      return;
    }
      
    fprintf(fp,"TITLE = \"%s\"\n","debug.dat");
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");
    fprintf(fp,"ZONE T=\"%s\"\n",filename.c_str());
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nno,ned);

    FOR_INO {
      fprintf(fp,"%18.15e %18.15e %18.15e\n",x_no[ino][0]+x0[0],x_no[ino][1]+x0[1],x_no[ino][2]+x0[2]);
    }

    FOR_IED {
      //cout << "faoed[ied]: " << faoed[ied][0] << " " << faoed[ied][1] << endl;
      fprintf(fp,"%d %d %d\n",nooed[ied][0]+1,nooed[ied][1]+1,nooed[ied][1]+1);
    }

    fclose(fp);

  }

  void writeSplot(const string& filename) const {

    double zero[3] = { 0.0, 0.0, 0.0 };
    writeSplot(filename,zero);
  
  }

  void writeSplot(const string& filename,const double x0[3]) const {
    
    FILE * fp = fopen(filename.c_str(),"w");
    if (fp == NULL) {
      cout << "Error: cannot open filename: " << filename << endl;
      return;
    }

    for (int ied = 0; ied < ned; ++ied) {
      fprintf(fp,"%18.15le %18.15le %18.15le\n",x_no[nooed[ied][0]][0]+x0[0],x_no[nooed[ied][0]][1]+x0[1],x_no[nooed[ied][0]][2]+x0[2]);
      fprintf(fp,"%18.15le %18.15le %18.15le\n",x_no[nooed[ied][1]][0]+x0[0],x_no[nooed[ied][1]][1]+x0[1],x_no[nooed[ied][1]][2]+x0[2]);
      fprintf(fp,"\n\n");
    }

    fclose(fp);
    
  }
  
#define IF_FAOED_0 if (faoed[ied][0] >= 0)
#define IF_FAOED_1 if (faoed[ied][1] >= 0)
#define CUT cut_surf
#include "cut.hpp"
#undef CUT
#undef IF_FAOED_0
#undef IF_FAOED_1

#define IF_FAOED_0
#define IF_FAOED_1
#define CUT cut
#include "cut.hpp"
#undef CUT
#undef IF_FAOED_0
#undef IF_FAOED_1

  void addRhombicDodecahedron(const double L) {

    // allow us to add the rhombic dodecahedron to whatever is there...

    const int nno_orig = nno;
    resize_nno(nno_orig+14);

    const double Lo2 = L/2.0;
    const double Lo4 = L/4.0;
    x_no[nno_orig+ 0][0] = -Lo4; x_no[nno_orig+ 0][1] = -Lo4; x_no[nno_orig+ 0][2] =  Lo4;
    x_no[nno_orig+ 1][0] =    0; x_no[nno_orig+ 1][1] = -Lo2; x_no[nno_orig+ 1][2] =    0;
    x_no[nno_orig+ 2][0] = -Lo4; x_no[nno_orig+ 2][1] = -Lo4; x_no[nno_orig+ 2][2] = -Lo4;
    x_no[nno_orig+ 3][0] = -Lo2; x_no[nno_orig+ 3][1] =    0; x_no[nno_orig+ 3][2] =    0;
    x_no[nno_orig+ 4][0] =    0; x_no[nno_orig+ 4][1] =    0; x_no[nno_orig+ 4][2] =  Lo2;
    x_no[nno_orig+ 5][0] =    0; x_no[nno_orig+ 5][1] =    0; x_no[nno_orig+ 5][2] = -Lo2;
    x_no[nno_orig+ 6][0] = -Lo4; x_no[nno_orig+ 6][1] =  Lo4; x_no[nno_orig+ 6][2] =  Lo4;
    x_no[nno_orig+ 7][0] =    0; x_no[nno_orig+ 7][1] =  Lo2; x_no[nno_orig+ 7][2] =    0;
    x_no[nno_orig+ 8][0] = -Lo4; x_no[nno_orig+ 8][1] =  Lo4; x_no[nno_orig+ 8][2] = -Lo4;
    x_no[nno_orig+ 9][0] =  Lo4; x_no[nno_orig+ 9][1] = -Lo4; x_no[nno_orig+ 9][2] =  Lo4;
    x_no[nno_orig+10][0] =  Lo2; x_no[nno_orig+10][1] =    0; x_no[nno_orig+10][2] =    0;
    x_no[nno_orig+11][0] =  Lo4; x_no[nno_orig+11][1] = -Lo4; x_no[nno_orig+11][2] = -Lo4;
    x_no[nno_orig+12][0] =  Lo4; x_no[nno_orig+12][1] =  Lo4; x_no[nno_orig+12][2] =  Lo4;
    x_no[nno_orig+13][0] =  Lo4; x_no[nno_orig+13][1] =  Lo4; x_no[nno_orig+13][2] = -Lo4;

    const int ned_orig = ned;
    resize_ned(ned_orig+24);

    // note that we are in -8 indexing
    nooed[ned_orig+ 0][0] =  4+nno_orig; nooed[ned_orig+ 0][1] =  9+nno_orig; faoed[ned_orig+ 0][0] =  -8; faoed[ned_orig+ 0][1] = -11;
    nooed[ned_orig+ 1][0] =  9+nno_orig; nooed[ned_orig+ 1][1] = 10+nno_orig; faoed[ned_orig+ 1][0] =  -8; faoed[ned_orig+ 1][1] = -17;
    nooed[ned_orig+ 2][0] = 10+nno_orig; nooed[ned_orig+ 2][1] = 12+nno_orig; faoed[ned_orig+ 2][0] =  -8; faoed[ned_orig+ 2][1] = -16;
    nooed[ned_orig+ 3][0] = 12+nno_orig; nooed[ned_orig+ 3][1] =  4+nno_orig; faoed[ned_orig+ 3][0] =  -8; faoed[ned_orig+ 3][1] =  -9;
    nooed[ned_orig+ 4][0] = 12+nno_orig; nooed[ned_orig+ 4][1] =  7+nno_orig; faoed[ned_orig+ 4][0] =  -9; faoed[ned_orig+ 4][1] = -16;
    nooed[ned_orig+ 5][0] =  7+nno_orig; nooed[ned_orig+ 5][1] =  6+nno_orig; faoed[ned_orig+ 5][0] =  -9; faoed[ned_orig+ 5][1] = -18;
    nooed[ned_orig+ 6][0] =  6+nno_orig; nooed[ned_orig+ 6][1] =  4+nno_orig; faoed[ned_orig+ 6][0] =  -9; faoed[ned_orig+ 6][1] = -10;
    nooed[ned_orig+ 7][0] =  6+nno_orig; nooed[ned_orig+ 7][1] =  3+nno_orig; faoed[ned_orig+ 7][0] = -10; faoed[ned_orig+ 7][1] = -18;
    nooed[ned_orig+ 8][0] =  3+nno_orig; nooed[ned_orig+ 8][1] =  0+nno_orig; faoed[ned_orig+ 8][0] = -10; faoed[ned_orig+ 8][1] = -19;
    nooed[ned_orig+ 9][0] =  0+nno_orig; nooed[ned_orig+ 9][1] =  4+nno_orig; faoed[ned_orig+ 9][0] = -10; faoed[ned_orig+ 9][1] = -11;
    nooed[ned_orig+10][0] =  0+nno_orig; nooed[ned_orig+10][1] =  1+nno_orig; faoed[ned_orig+10][0] = -11; faoed[ned_orig+10][1] = -19;
    nooed[ned_orig+11][0] =  1+nno_orig; nooed[ned_orig+11][1] =  9+nno_orig; faoed[ned_orig+11][0] = -11; faoed[ned_orig+11][1] = -17;
    nooed[ned_orig+12][0] =  5+nno_orig; nooed[ned_orig+12][1] =  2+nno_orig; faoed[ned_orig+12][0] = -12; faoed[ned_orig+12][1] = -15;
    nooed[ned_orig+13][0] =  2+nno_orig; nooed[ned_orig+13][1] =  3+nno_orig; faoed[ned_orig+13][0] = -12; faoed[ned_orig+13][1] = -19;
    nooed[ned_orig+14][0] =  3+nno_orig; nooed[ned_orig+14][1] =  8+nno_orig; faoed[ned_orig+14][0] = -12; faoed[ned_orig+14][1] = -18;
    nooed[ned_orig+15][0] =  8+nno_orig; nooed[ned_orig+15][1] =  5+nno_orig; faoed[ned_orig+15][0] = -12; faoed[ned_orig+15][1] = -13;
    nooed[ned_orig+16][0] =  8+nno_orig; nooed[ned_orig+16][1] =  7+nno_orig; faoed[ned_orig+16][0] = -13; faoed[ned_orig+16][1] = -18;
    nooed[ned_orig+17][0] =  7+nno_orig; nooed[ned_orig+17][1] = 13+nno_orig; faoed[ned_orig+17][0] = -13; faoed[ned_orig+17][1] = -16;
    nooed[ned_orig+18][0] = 13+nno_orig; nooed[ned_orig+18][1] =  5+nno_orig; faoed[ned_orig+18][0] = -13; faoed[ned_orig+18][1] = -14;
    nooed[ned_orig+19][0] = 13+nno_orig; nooed[ned_orig+19][1] = 10+nno_orig; faoed[ned_orig+19][0] = -14; faoed[ned_orig+19][1] = -16;
    nooed[ned_orig+20][0] = 10+nno_orig; nooed[ned_orig+20][1] = 11+nno_orig; faoed[ned_orig+20][0] = -14; faoed[ned_orig+20][1] = -17;
    nooed[ned_orig+21][0] = 11+nno_orig; nooed[ned_orig+21][1] =  5+nno_orig; faoed[ned_orig+21][0] = -14; faoed[ned_orig+21][1] = -15;
    nooed[ned_orig+22][0] = 11+nno_orig; nooed[ned_orig+22][1] =  1+nno_orig; faoed[ned_orig+22][0] = -15; faoed[ned_orig+22][1] = -17;
    nooed[ned_orig+23][0] =  1+nno_orig; nooed[ned_orig+23][1] =  2+nno_orig; faoed[ned_orig+23][0] = -15; faoed[ned_orig+23][1] = -19;

    /*
      double xc[3] = { 0.0, 0.0, 0.0 };
      writeTecplot(0,xc);
      getchar();
    */
  }
  void addTruncatedOctahedronRoot2(const double L) {

    // allow us to add the truncated octahedron to whatever is there...

    const int nno_orig = nno;
    resize_nno(nno_orig+24);

    const double Lo2 = L/2.0;
    const double Lo4 = L/4.0;
    const double Lr2o4 = L*sqrt(2.0)/4.0;
    const double Lr2o8 = L*sqrt(2.0)/8.0;
    const double L3r2o8 = L*3.0*sqrt(2.0)/8.0;
    x_no[nno_orig+ 0][0] =  Lo2; x_no[nno_orig+ 0][1] =  -Lr2o8; x_no[nno_orig+ 0][2] =   Lr2o8;
    x_no[nno_orig+ 1][0] =  Lo2; x_no[nno_orig+ 1][1] =   Lr2o8; x_no[nno_orig+ 1][2] =   Lr2o8;
    x_no[nno_orig+ 2][0] =  Lo2; x_no[nno_orig+ 2][1] =  -Lr2o8; x_no[nno_orig+ 2][2] =  -Lr2o8;
    x_no[nno_orig+ 3][0] =  Lo2; x_no[nno_orig+ 3][1] =   Lr2o8; x_no[nno_orig+ 3][2] =  -Lr2o8;
    x_no[nno_orig+ 4][0] =    0; x_no[nno_orig+ 4][1] =  L3r2o8; x_no[nno_orig+ 4][2] =   Lr2o8;
    x_no[nno_orig+ 5][0] = -Lo4; x_no[nno_orig+ 5][1] =   Lr2o4; x_no[nno_orig+ 5][2] =   Lr2o4;
    x_no[nno_orig+ 6][0] =    0; x_no[nno_orig+ 6][1] =   Lr2o8; x_no[nno_orig+ 6][2] =  L3r2o8;
    x_no[nno_orig+ 7][0] =  Lo4; x_no[nno_orig+ 7][1] =   Lr2o4; x_no[nno_orig+ 7][2] =   Lr2o4;
    x_no[nno_orig+ 8][0] =    0; x_no[nno_orig+ 8][1] = -L3r2o8; x_no[nno_orig+ 8][2] =   Lr2o8;
    x_no[nno_orig+ 9][0] = -Lo4; x_no[nno_orig+ 9][1] =  -Lr2o4; x_no[nno_orig+ 9][2] =   Lr2o4;
    x_no[nno_orig+10][0] =    0; x_no[nno_orig+10][1] =  -Lr2o8; x_no[nno_orig+10][2] =  L3r2o8;
    x_no[nno_orig+11][0] =  Lo4; x_no[nno_orig+11][1] =  -Lr2o4; x_no[nno_orig+11][2] =   Lr2o4;
    x_no[nno_orig+12][0] = -Lo2; x_no[nno_orig+12][1] =  -Lr2o8; x_no[nno_orig+12][2] =   Lr2o8;
    x_no[nno_orig+13][0] = -Lo2; x_no[nno_orig+13][1] =   Lr2o8; x_no[nno_orig+13][2] =   Lr2o8;
    x_no[nno_orig+14][0] = -Lo2; x_no[nno_orig+14][1] =  -Lr2o8; x_no[nno_orig+14][2] =  -Lr2o8;
    x_no[nno_orig+15][0] = -Lo2; x_no[nno_orig+15][1] =   Lr2o8; x_no[nno_orig+15][2] =  -Lr2o8;
    x_no[nno_orig+16][0] =    0; x_no[nno_orig+16][1] =  L3r2o8; x_no[nno_orig+16][2] =  -Lr2o8;
    x_no[nno_orig+17][0] =  Lo4; x_no[nno_orig+17][1] =   Lr2o4; x_no[nno_orig+17][2] =  -Lr2o4;
    x_no[nno_orig+18][0] =    0; x_no[nno_orig+18][1] =   Lr2o8; x_no[nno_orig+18][2] = -L3r2o8;
    x_no[nno_orig+19][0] = -Lo4; x_no[nno_orig+19][1] =   Lr2o4; x_no[nno_orig+19][2] =  -Lr2o4;
    x_no[nno_orig+20][0] =    0; x_no[nno_orig+20][1] = -L3r2o8; x_no[nno_orig+20][2] =  -Lr2o8;
    x_no[nno_orig+21][0] =  Lo4; x_no[nno_orig+21][1] =  -Lr2o4; x_no[nno_orig+21][2] =  -Lr2o4;
    x_no[nno_orig+22][0] =    0; x_no[nno_orig+22][1] =  -Lr2o8; x_no[nno_orig+22][2] = -L3r2o8;
    x_no[nno_orig+23][0] = -Lo4; x_no[nno_orig+23][1] =  -Lr2o4; x_no[nno_orig+23][2] =  -Lr2o4;

    const int ned_orig = ned;
    resize_ned(ned_orig+36);

    // note that we are in -8 indexing
    nooed[ned_orig+ 0][0] =  2+nno_orig; nooed[ned_orig+ 0][1] =  3+nno_orig; faoed[ned_orig+ 0][0] = -16; faoed[ned_orig+ 0][1] = -12;
    nooed[ned_orig+ 1][0] =  3+nno_orig; nooed[ned_orig+ 1][1] =  1+nno_orig; faoed[ned_orig+ 1][0] = -16; faoed[ned_orig+ 1][1] =  -9;
    nooed[ned_orig+ 2][0] =  1+nno_orig; nooed[ned_orig+ 2][1] =  0+nno_orig; faoed[ned_orig+ 2][0] = -16; faoed[ned_orig+ 2][1] = -10;
    nooed[ned_orig+ 3][0] =  0+nno_orig; nooed[ned_orig+ 3][1] =  2+nno_orig; faoed[ned_orig+ 3][0] = -16; faoed[ned_orig+ 3][1] = -11;
    nooed[ned_orig+ 4][0] =  3+nno_orig; nooed[ned_orig+ 4][1] = 17+nno_orig; faoed[ned_orig+ 4][0] =  -9; faoed[ned_orig+ 4][1] = -12;
    nooed[ned_orig+ 5][0] = 17+nno_orig; nooed[ned_orig+ 5][1] = 16+nno_orig; faoed[ned_orig+ 5][0] =  -9; faoed[ned_orig+ 5][1] = -20;
    nooed[ned_orig+ 6][0] = 16+nno_orig; nooed[ned_orig+ 6][1] =  4+nno_orig; faoed[ned_orig+ 6][0] =  -9; faoed[ned_orig+ 6][1] = -14;
    nooed[ned_orig+ 7][0] =  4+nno_orig; nooed[ned_orig+ 7][1] =  7+nno_orig; faoed[ned_orig+ 7][0] =  -9; faoed[ned_orig+ 7][1] = -18;
    nooed[ned_orig+ 8][0] =  7+nno_orig; nooed[ned_orig+ 8][1] =  1+nno_orig; faoed[ned_orig+ 8][0] =  -9; faoed[ned_orig+ 8][1] = -10;
    nooed[ned_orig+ 9][0] =  7+nno_orig; nooed[ned_orig+ 9][1] =  6+nno_orig; faoed[ned_orig+ 9][0] = -10; faoed[ned_orig+ 9][1] = -18;
    nooed[ned_orig+10][0] =  6+nno_orig; nooed[ned_orig+10][1] = 10+nno_orig; faoed[ned_orig+10][0] = -10; faoed[ned_orig+10][1] = -13;
    nooed[ned_orig+11][0] = 10+nno_orig; nooed[ned_orig+11][1] = 11+nno_orig; faoed[ned_orig+11][0] = -10; faoed[ned_orig+11][1] = -19;
    nooed[ned_orig+12][0] = 11+nno_orig; nooed[ned_orig+12][1] =  0+nno_orig; faoed[ned_orig+12][0] = -10; faoed[ned_orig+12][1] = -11;
    nooed[ned_orig+13][0] = 11+nno_orig; nooed[ned_orig+13][1] =  8+nno_orig; faoed[ned_orig+13][0] = -11; faoed[ned_orig+13][1] = -19;
    nooed[ned_orig+14][0] =  8+nno_orig; nooed[ned_orig+14][1] = 20+nno_orig; faoed[ned_orig+14][0] = -11; faoed[ned_orig+14][1] =  -8;
    nooed[ned_orig+15][0] = 20+nno_orig; nooed[ned_orig+15][1] = 21+nno_orig; faoed[ned_orig+15][0] = -11; faoed[ned_orig+15][1] = -21;
    nooed[ned_orig+16][0] = 21+nno_orig; nooed[ned_orig+16][1] =  2+nno_orig; faoed[ned_orig+16][0] = -11; faoed[ned_orig+16][1] = -12;
    nooed[ned_orig+17][0] = 21+nno_orig; nooed[ned_orig+17][1] = 22+nno_orig; faoed[ned_orig+17][0] = -12; faoed[ned_orig+17][1] = -21;
    nooed[ned_orig+18][0] = 22+nno_orig; nooed[ned_orig+18][1] = 18+nno_orig; faoed[ned_orig+18][0] = -12; faoed[ned_orig+18][1] = -15;
    nooed[ned_orig+19][0] = 18+nno_orig; nooed[ned_orig+19][1] = 17+nno_orig; faoed[ned_orig+19][0] = -12; faoed[ned_orig+19][1] = -20;
    nooed[ned_orig+20][0] = 12+nno_orig; nooed[ned_orig+20][1] = 13+nno_orig; faoed[ned_orig+20][0] = -17; faoed[ned_orig+20][1] = -13;
    nooed[ned_orig+21][0] = 13+nno_orig; nooed[ned_orig+21][1] = 15+nno_orig; faoed[ned_orig+21][0] = -17; faoed[ned_orig+21][1] = -14;
    nooed[ned_orig+22][0] = 15+nno_orig; nooed[ned_orig+22][1] = 14+nno_orig; faoed[ned_orig+22][0] = -17; faoed[ned_orig+22][1] = -15;
    nooed[ned_orig+23][0] = 14+nno_orig; nooed[ned_orig+23][1] = 12+nno_orig; faoed[ned_orig+23][0] = -17; faoed[ned_orig+23][1] =  -8;
    nooed[ned_orig+24][0] = 13+nno_orig; nooed[ned_orig+24][1] =  5+nno_orig; faoed[ned_orig+24][0] = -14; faoed[ned_orig+24][1] = -13;
    nooed[ned_orig+25][0] =  5+nno_orig; nooed[ned_orig+25][1] =  4+nno_orig; faoed[ned_orig+25][0] = -14; faoed[ned_orig+25][1] = -18;
    nooed[ned_orig+26][0] = 16+nno_orig; nooed[ned_orig+26][1] = 19+nno_orig; faoed[ned_orig+26][0] = -14; faoed[ned_orig+26][1] = -20;
    nooed[ned_orig+27][0] = 19+nno_orig; nooed[ned_orig+27][1] = 15+nno_orig; faoed[ned_orig+27][0] = -14; faoed[ned_orig+27][1] = -15;
    nooed[ned_orig+28][0] = 19+nno_orig; nooed[ned_orig+28][1] = 18+nno_orig; faoed[ned_orig+28][0] = -15; faoed[ned_orig+28][1] = -20;
    nooed[ned_orig+29][0] = 22+nno_orig; nooed[ned_orig+29][1] = 23+nno_orig; faoed[ned_orig+29][0] = -15; faoed[ned_orig+29][1] = -21;
    nooed[ned_orig+30][0] = 23+nno_orig; nooed[ned_orig+30][1] = 14+nno_orig; faoed[ned_orig+30][0] = -15; faoed[ned_orig+30][1] =  -8;
    nooed[ned_orig+31][0] = 23+nno_orig; nooed[ned_orig+31][1] = 20+nno_orig; faoed[ned_orig+31][0] =  -8; faoed[ned_orig+31][1] = -21;
    nooed[ned_orig+32][0] =  8+nno_orig; nooed[ned_orig+32][1] =  9+nno_orig; faoed[ned_orig+32][0] =  -8; faoed[ned_orig+32][1] = -19;
    nooed[ned_orig+33][0] =  9+nno_orig; nooed[ned_orig+33][1] = 12+nno_orig; faoed[ned_orig+33][0] =  -8; faoed[ned_orig+33][1] = -13;
    nooed[ned_orig+34][0] =  9+nno_orig; nooed[ned_orig+34][1] = 10+nno_orig; faoed[ned_orig+34][0] = -13; faoed[ned_orig+34][1] = -19;
    nooed[ned_orig+35][0] =  6+nno_orig; nooed[ned_orig+35][1] =  5+nno_orig; faoed[ned_orig+35][0] = -13; faoed[ned_orig+35][1] = -18;

    /*
      double xc[3] = { 0.0, 0.0, 0.0 };
      writeTecplot(0,xc);
      getchar();
    */
  }

  void addTruncatedOctahedronRoot3(const double L) {

    // allow us to add the truncated octahedron to whatever is there...

    const int nno_orig = nno;
    resize_nno(nno_orig+24);

    const double Lo2 = L/2.0;
    const double Lo4 = L/4.0;
    const double Lr3o3 = L*sqrt(3.0)/3.0;
    const double Lr3o6 = L*sqrt(3.0)/6.0;
    const double Lr3o4 = L*sqrt(3.0)/4.0;
    x_no[nno_orig+ 0][0] =  Lo2; x_no[nno_orig+ 0][1] =  Lr3o6; x_no[nno_orig+ 0][2] =  Lr3o6;
    x_no[nno_orig+ 1][0] =  Lo2; x_no[nno_orig+ 1][1] = -Lr3o6; x_no[nno_orig+ 1][2] =  Lr3o6;
    x_no[nno_orig+ 2][0] =  Lo2; x_no[nno_orig+ 2][1] = -Lr3o6; x_no[nno_orig+ 2][2] = -Lr3o6;
    x_no[nno_orig+ 3][0] =  Lo2; x_no[nno_orig+ 3][1] =  Lr3o6; x_no[nno_orig+ 3][2] = -Lr3o6;
    x_no[nno_orig+ 4][0] = -Lo2; x_no[nno_orig+ 4][1] =  Lr3o6; x_no[nno_orig+ 4][2] =  Lr3o6;
    x_no[nno_orig+ 5][0] = -Lo2; x_no[nno_orig+ 5][1] = -Lr3o6; x_no[nno_orig+ 5][2] =  Lr3o6;
    x_no[nno_orig+ 6][0] = -Lo2; x_no[nno_orig+ 6][1] = -Lr3o6; x_no[nno_orig+ 6][2] = -Lr3o6;
    x_no[nno_orig+ 7][0] = -Lo2; x_no[nno_orig+ 7][1] =  Lr3o6; x_no[nno_orig+ 7][2] = -Lr3o6;
    x_no[nno_orig+ 8][0] =    0; x_no[nno_orig+ 8][1] =  Lr3o6; x_no[nno_orig+ 8][2] =  Lr3o3;
    x_no[nno_orig+ 9][0] =  Lo4; x_no[nno_orig+ 9][1] =  Lr3o4; x_no[nno_orig+ 9][2] =  Lr3o4;
    x_no[nno_orig+10][0] =    0; x_no[nno_orig+10][1] =  Lr3o3; x_no[nno_orig+10][2] =  Lr3o6;
    x_no[nno_orig+11][0] = -Lo4; x_no[nno_orig+11][1] =  Lr3o4; x_no[nno_orig+11][2] =  Lr3o4;
    x_no[nno_orig+12][0] =    0; x_no[nno_orig+12][1] = -Lr3o6; x_no[nno_orig+12][2] =  Lr3o3;
    x_no[nno_orig+13][0] =  Lo4; x_no[nno_orig+13][1] = -Lr3o4; x_no[nno_orig+13][2] =  Lr3o4;
    x_no[nno_orig+14][0] =    0; x_no[nno_orig+14][1] = -Lr3o3; x_no[nno_orig+14][2] =  Lr3o6;
    x_no[nno_orig+15][0] = -Lo4; x_no[nno_orig+15][1] = -Lr3o4; x_no[nno_orig+15][2] =  Lr3o4;
    x_no[nno_orig+16][0] =    0; x_no[nno_orig+16][1] =  Lr3o3; x_no[nno_orig+16][2] = -Lr3o6;
    x_no[nno_orig+17][0] =  Lo4; x_no[nno_orig+17][1] =  Lr3o4; x_no[nno_orig+17][2] = -Lr3o4;
    x_no[nno_orig+18][0] =    0; x_no[nno_orig+18][1] =  Lr3o6; x_no[nno_orig+18][2] = -Lr3o3;
    x_no[nno_orig+19][0] = -Lo4; x_no[nno_orig+19][1] =  Lr3o4; x_no[nno_orig+19][2] = -Lr3o4;
    x_no[nno_orig+20][0] =    0; x_no[nno_orig+20][1] = -Lr3o3; x_no[nno_orig+20][2] = -Lr3o6;
    x_no[nno_orig+21][0] =  Lo4; x_no[nno_orig+21][1] = -Lr3o4; x_no[nno_orig+21][2] = -Lr3o4;
    x_no[nno_orig+22][0] =    0; x_no[nno_orig+22][1] = -Lr3o6; x_no[nno_orig+22][2] = -Lr3o3;
    x_no[nno_orig+23][0] = -Lo4; x_no[nno_orig+23][1] = -Lr3o4; x_no[nno_orig+23][2] = -Lr3o4;

    const int ned_orig = ned;
    resize_ned(ned_orig+36);

    // note that we are in -8 indexing
    nooed[ned_orig+ 0][0] =  2+nno_orig; nooed[ned_orig+ 0][1] =  3+nno_orig; faoed[ned_orig+ 0][0] =  -8; faoed[ned_orig+ 0][1] = -12;
    nooed[ned_orig+ 1][0] =  3+nno_orig; nooed[ned_orig+ 1][1] =  0+nno_orig; faoed[ned_orig+ 1][0] =  -8; faoed[ned_orig+ 1][1] =  -9;
    nooed[ned_orig+ 2][0] =  0+nno_orig; nooed[ned_orig+ 2][1] =  1+nno_orig; faoed[ned_orig+ 2][0] =  -8; faoed[ned_orig+ 2][1] = -10;
    nooed[ned_orig+ 3][0] =  1+nno_orig; nooed[ned_orig+ 3][1] =  2+nno_orig; faoed[ned_orig+ 3][0] =  -8; faoed[ned_orig+ 3][1] = -11;
    nooed[ned_orig+ 4][0] =  3+nno_orig; nooed[ned_orig+ 4][1] = 17+nno_orig; faoed[ned_orig+ 4][0] =  -9; faoed[ned_orig+ 4][1] = -12;
    nooed[ned_orig+ 5][0] = 17+nno_orig; nooed[ned_orig+ 5][1] = 16+nno_orig; faoed[ned_orig+ 5][0] =  -9; faoed[ned_orig+ 5][1] = -20;
    nooed[ned_orig+ 6][0] = 16+nno_orig; nooed[ned_orig+ 6][1] = 10+nno_orig; faoed[ned_orig+ 6][0] =  -9; faoed[ned_orig+ 6][1] = -14;
    nooed[ned_orig+ 7][0] = 10+nno_orig; nooed[ned_orig+ 7][1] =  9+nno_orig; faoed[ned_orig+ 7][0] =  -9; faoed[ned_orig+ 7][1] = -18;
    nooed[ned_orig+ 8][0] =  9+nno_orig; nooed[ned_orig+ 8][1] =  0+nno_orig; faoed[ned_orig+ 8][0] =  -9; faoed[ned_orig+ 8][1] = -10;
    nooed[ned_orig+ 9][0] =  9+nno_orig; nooed[ned_orig+ 9][1] =  8+nno_orig; faoed[ned_orig+ 9][0] = -10; faoed[ned_orig+ 9][1] = -18;
    nooed[ned_orig+10][0] =  8+nno_orig; nooed[ned_orig+10][1] = 12+nno_orig; faoed[ned_orig+10][0] = -10; faoed[ned_orig+10][1] = -17;
    nooed[ned_orig+11][0] = 12+nno_orig; nooed[ned_orig+11][1] = 13+nno_orig; faoed[ned_orig+11][0] = -10; faoed[ned_orig+11][1] = -19;
    nooed[ned_orig+12][0] = 13+nno_orig; nooed[ned_orig+12][1] =  1+nno_orig; faoed[ned_orig+12][0] = -10; faoed[ned_orig+12][1] = -11;
    nooed[ned_orig+13][0] = 13+nno_orig; nooed[ned_orig+13][1] = 14+nno_orig; faoed[ned_orig+13][0] = -11; faoed[ned_orig+13][1] = -19;
    nooed[ned_orig+14][0] = 14+nno_orig; nooed[ned_orig+14][1] = 20+nno_orig; faoed[ned_orig+14][0] = -11; faoed[ned_orig+14][1] = -16;
    nooed[ned_orig+15][0] = 20+nno_orig; nooed[ned_orig+15][1] = 21+nno_orig; faoed[ned_orig+15][0] = -11; faoed[ned_orig+15][1] = -21;
    nooed[ned_orig+16][0] = 21+nno_orig; nooed[ned_orig+16][1] =  2+nno_orig; faoed[ned_orig+16][0] = -11; faoed[ned_orig+16][1] = -12;
    nooed[ned_orig+17][0] = 21+nno_orig; nooed[ned_orig+17][1] = 22+nno_orig; faoed[ned_orig+17][0] = -12; faoed[ned_orig+17][1] = -21;
    nooed[ned_orig+18][0] = 22+nno_orig; nooed[ned_orig+18][1] = 18+nno_orig; faoed[ned_orig+18][0] = -12; faoed[ned_orig+18][1] = -15;
    nooed[ned_orig+19][0] = 18+nno_orig; nooed[ned_orig+19][1] = 17+nno_orig; faoed[ned_orig+19][0] = -12; faoed[ned_orig+19][1] = -20;
    nooed[ned_orig+20][0] =  4+nno_orig; nooed[ned_orig+20][1] =  7+nno_orig; faoed[ned_orig+20][0] = -13; faoed[ned_orig+20][1] = -14;
    nooed[ned_orig+21][0] =  7+nno_orig; nooed[ned_orig+21][1] =  6+nno_orig; faoed[ned_orig+21][0] = -13; faoed[ned_orig+21][1] = -15;
    nooed[ned_orig+22][0] =  6+nno_orig; nooed[ned_orig+22][1] =  5+nno_orig; faoed[ned_orig+22][0] = -13; faoed[ned_orig+22][1] = -16;
    nooed[ned_orig+23][0] =  5+nno_orig; nooed[ned_orig+23][1] =  4+nno_orig; faoed[ned_orig+23][0] = -13; faoed[ned_orig+23][1] = -17;
    nooed[ned_orig+24][0] =  4+nno_orig; nooed[ned_orig+24][1] = 11+nno_orig; faoed[ned_orig+24][0] = -14; faoed[ned_orig+24][1] = -17;
    nooed[ned_orig+25][0] = 11+nno_orig; nooed[ned_orig+25][1] = 10+nno_orig; faoed[ned_orig+25][0] = -14; faoed[ned_orig+25][1] = -18;
    nooed[ned_orig+26][0] = 16+nno_orig; nooed[ned_orig+26][1] = 19+nno_orig; faoed[ned_orig+26][0] = -14; faoed[ned_orig+26][1] = -20;
    nooed[ned_orig+27][0] = 19+nno_orig; nooed[ned_orig+27][1] =  7+nno_orig; faoed[ned_orig+27][0] = -14; faoed[ned_orig+27][1] = -15;
    nooed[ned_orig+28][0] = 19+nno_orig; nooed[ned_orig+28][1] = 18+nno_orig; faoed[ned_orig+28][0] = -15; faoed[ned_orig+28][1] = -20;
    nooed[ned_orig+29][0] = 22+nno_orig; nooed[ned_orig+29][1] = 23+nno_orig; faoed[ned_orig+29][0] = -15; faoed[ned_orig+29][1] = -21;
    nooed[ned_orig+30][0] = 23+nno_orig; nooed[ned_orig+30][1] =  6+nno_orig; faoed[ned_orig+30][0] = -15; faoed[ned_orig+30][1] = -16;
    nooed[ned_orig+31][0] = 23+nno_orig; nooed[ned_orig+31][1] = 20+nno_orig; faoed[ned_orig+31][0] = -16; faoed[ned_orig+31][1] = -21;
    nooed[ned_orig+32][0] = 14+nno_orig; nooed[ned_orig+32][1] = 15+nno_orig; faoed[ned_orig+32][0] = -16; faoed[ned_orig+32][1] = -19;
    nooed[ned_orig+33][0] = 15+nno_orig; nooed[ned_orig+33][1] =  5+nno_orig; faoed[ned_orig+33][0] = -16; faoed[ned_orig+33][1] = -17;
    nooed[ned_orig+34][0] = 15+nno_orig; nooed[ned_orig+34][1] = 12+nno_orig; faoed[ned_orig+34][0] = -17; faoed[ned_orig+34][1] = -19;
    nooed[ned_orig+35][0] =  8+nno_orig; nooed[ned_orig+35][1] = 11+nno_orig; faoed[ned_orig+35][0] = -17; faoed[ned_orig+35][1] = -18;

    /*
      double xc[3] = { 0.0, 0.0, 0.0 };
      writeTecplot(0,xc);
      getchar();
    */
  }

  void addCube(const double L) {

    addCube(L,L,L);

  }

  void addCube(const double Lx,const double Ly,const double Lz) {
    
    // allow us to add the cube to whatever is there...

    const int nno_orig = nno;
    resize_nno(nno_orig+8);
    
    x_no[nno_orig+0][0] = -Lx; x_no[nno_orig+0][1] = -Ly; x_no[nno_orig+0][2] = -Lz;
    x_no[nno_orig+1][0] =  Lx; x_no[nno_orig+1][1] = -Ly; x_no[nno_orig+1][2] = -Lz;
    x_no[nno_orig+2][0] =  Lx; x_no[nno_orig+2][1] =  Ly; x_no[nno_orig+2][2] = -Lz;
    x_no[nno_orig+3][0] = -Lx; x_no[nno_orig+3][1] =  Ly; x_no[nno_orig+3][2] = -Lz;
    x_no[nno_orig+4][0] = -Lx; x_no[nno_orig+4][1] = -Ly; x_no[nno_orig+4][2] =  Lz;
    x_no[nno_orig+5][0] =  Lx; x_no[nno_orig+5][1] = -Ly; x_no[nno_orig+5][2] =  Lz;
    x_no[nno_orig+6][0] =  Lx; x_no[nno_orig+6][1] =  Ly; x_no[nno_orig+6][2] =  Lz;
    x_no[nno_orig+7][0] = -Lx; x_no[nno_orig+7][1] =  Ly; x_no[nno_orig+7][2] =  Lz;

    const int ned_orig = ned;
    resize_ned(ned_orig+12);

    // z-edges...
    nooed[ned_orig+0][0] = 0+nno_orig;  nooed[ned_orig+0][1] = 4+nno_orig;  faoed[ned_orig+0][0] = -2;  faoed[ned_orig+0][1] = -4;
    nooed[ned_orig+1][0] = 1+nno_orig;  nooed[ned_orig+1][1] = 5+nno_orig;  faoed[ned_orig+1][0] = -4;  faoed[ned_orig+1][1] = -3;
    nooed[ned_orig+2][0] = 2+nno_orig;  nooed[ned_orig+2][1] = 6+nno_orig;  faoed[ned_orig+2][0] = -3;  faoed[ned_orig+2][1] = -5;
    nooed[ned_orig+3][0] = 3+nno_orig;  nooed[ned_orig+3][1] = 7+nno_orig;  faoed[ned_orig+3][0] = -5;  faoed[ned_orig+3][1] = -2;
    // x-edges...
    nooed[ned_orig+4][0] = 0+nno_orig;  nooed[ned_orig+4][1] = 1+nno_orig;  faoed[ned_orig+4][0] = -4;  faoed[ned_orig+4][1] = -6;
    nooed[ned_orig+5][0] = 3+nno_orig;  nooed[ned_orig+5][1] = 2+nno_orig;  faoed[ned_orig+5][0] = -6;  faoed[ned_orig+5][1] = -5;
    nooed[ned_orig+6][0] = 7+nno_orig;  nooed[ned_orig+6][1] = 6+nno_orig;  faoed[ned_orig+6][0] = -5;  faoed[ned_orig+6][1] = -7;
    nooed[ned_orig+7][0] = 4+nno_orig;  nooed[ned_orig+7][1] = 5+nno_orig;  faoed[ned_orig+7][0] = -7;  faoed[ned_orig+7][1] = -4;
    // y-edges...
    nooed[ned_orig+8][0] = 0+nno_orig;  nooed[ned_orig+8][1] = 3+nno_orig;  faoed[ned_orig+8][0] = -6;  faoed[ned_orig+8][1] = -2;
    nooed[ned_orig+9][0] = 4+nno_orig;  nooed[ned_orig+9][1] = 7+nno_orig;  faoed[ned_orig+9][0] = -2;  faoed[ned_orig+9][1] = -7;
    nooed[ned_orig+10][0] = 5+nno_orig; nooed[ned_orig+10][1] = 6+nno_orig; faoed[ned_orig+10][0] = -7; faoed[ned_orig+10][1] = -3;
    nooed[ned_orig+11][0] = 1+nno_orig; nooed[ned_orig+11][1] = 2+nno_orig; faoed[ned_orig+11][0] = -3; faoed[ned_orig+11][1] = -6;

  }

  int completeCube(const double L) {

    static int debug_count_cc = 0;
    ++debug_count_cc;

    const bool debug = ( (getIntParam("DEBUG_RANK",-1)==mpi_rank) && (getIntParam("DEBUG_COUNT_CC",-1)==debug_count_cc) );

    if (ned == 0) {
      assert(nno == 0);
      addCube(L);
      return 1;
    }

    // step 1: there should be a cut surface at this point, with nooed referring
    // to the local nodes, and faoed as follows:
    // faoed[ied][0,1] >= 0: the index of the associated surface tri
    // faoed[ied][0] == -2 : x0 face
    // faoed[ied][0] == -3 : x1 face
    // faoed[ied][0] == -4 : y0 face
    // faoed[ied][0] == -5 : y1 face
    // faoed[ied][0] == -6 : z0 face
    // faoed[ied][0] == -7 : z1 face
    // faoed[ied][0] == -1 : original edge of the surface tri -- should all be cut away

    // note that it is possible that this could fail for certain peroidic boundaries because we are not
    // requesting enough data. Try cutting further to handle most of these cases.
    // everything. Another reason for failure is that the hull is not water tight.

    if (debug) {
      double xc[3] = { 0.0, 0.0, 0.0 };
      writeTecplot(0,xc);
    }

    {
      map<const int,pair<int,int> > faceMap;
      int nfa_ = 0;
      FOR_IED {
	if (faoed[ied][1] < 0) {
	  // we have found an edge that should be cut off. the
	  // other side will need a normal...
	  assert(faoed[ied][0] >= 0);
	  if (faceMap.find(faoed[ied][0]) == faceMap.end())
	    faceMap[faoed[ied][0]] = pair<int,int>(nooed[ied][0],nfa_++);
	}
      }

      if (!faceMap.empty()) {
	bool cut_flag[6] = { false, false, false, false, false, false };
	double cut_position[6];
	double (*fa_normal_)[3] = new double[nfa_][3];
	for (int ifa_ = 0; ifa_ < nfa_; ++ifa_)
	  FOR_I3 fa_normal_[ifa_][i] = 0.0;
	FOR_IED {
	  if (faoed[ied][1] >= 0) {
	    if (faoed[ied][0] >= 0) {
	      map<const int,pair<int,int> >::iterator iter = faceMap.find(faoed[ied][0]);
	      if (iter != faceMap.end()) {
		const int ino  = iter->second.first;
		const int ifa_ = iter->second.second;
		if ((nooed[ied][0] != ino)&&(nooed[ied][1] != ino)) {
		  const double this_n[3] = TRI_NORMAL_2(x_no[ino],x_no[nooed[ied][0]],x_no[nooed[ied][1]]);
		  FOR_I3 fa_normal_[ifa_][i] += this_n[i];
		}
	      }
	    }
	    map<const int,pair<int,int> >::iterator iter = faceMap.find(faoed[ied][1]);
	    if (iter != faceMap.end()) {
	      const int ino  = iter->second.first;
	      const int ifa_ = iter->second.second;
	      if ((nooed[ied][0] != ino)&&(nooed[ied][1] != ino)) {
		const double this_n[3] = TRI_NORMAL_2(x_no[ino],x_no[nooed[ied][1]],x_no[nooed[ied][0]]);
		FOR_I3 fa_normal_[ifa_][i] += this_n[i];
	      }
	    }
	  }
	}
	// now go back through the edges and set the cut positions...
	FOR_IED {
	  if (faoed[ied][1] < 0) {
	    // we have found an edge that should be cut off. the
	    // other side will need a normal...
	    assert(faoed[ied][0] >= 0);
	    map<const int,pair<int,int> >::iterator iter = faceMap.find(faoed[ied][0]);
	    assert(iter != faceMap.end());
	    const int ifa_ = iter->second.second;
	    assert(DOT_PRODUCT(fa_normal_[ifa_],fa_normal_[ifa_]) > 0.0);
	    const double dx_ed[3] = DIFF(x_no[nooed[ied][1]],x_no[nooed[ied][0]]);
	    const double cp[3] = CROSS_PRODUCT(dx_ed,fa_normal_[ifa_]);
	    if (fabs(cp[0]) >= max(fabs(cp[1]),fabs(cp[2]))) {
	      if (cp[0] < 0.0) {
		if (!cut_flag[0]) {
		  cut_flag[0] = true;
		  cut_position[0] = max(x_no[nooed[ied][0]][0],x_no[nooed[ied][1]][0]);
		}
		else {
		  cut_position[0] = max(cut_position[0],max(x_no[nooed[ied][0]][0],x_no[nooed[ied][1]][0]));
		}
	      }
	      else {
		if (!cut_flag[1]) {
		  cut_flag[1] = true;
		  cut_position[1] = min(x_no[nooed[ied][0]][0],x_no[nooed[ied][1]][0]);
		}
		else {
		  cut_position[1] = min(cut_position[1],max(x_no[nooed[ied][0]][0],x_no[nooed[ied][1]][0]));
		}
	      }
	    }
	    else if (fabs(cp[1]) >= max(fabs(cp[0]),fabs(cp[2]))) {
	      if (cp[1] < 0.0) {
		if (!cut_flag[2]) {
		  cut_flag[2] = true;
		  cut_position[2] = max(x_no[nooed[ied][0]][1],x_no[nooed[ied][1]][1]);
		}
		else {
		  cut_position[2] = max(cut_position[2],max(x_no[nooed[ied][0]][1],x_no[nooed[ied][1]][1]));
		}
	      }
	      else {
		if (!cut_flag[3]) {
		  cut_flag[3] = true;
		  cut_position[3] = min(x_no[nooed[ied][0]][1],x_no[nooed[ied][1]][1]);
		}
		else {
		  cut_position[3] = min(cut_position[3],max(x_no[nooed[ied][0]][1],x_no[nooed[ied][1]][1]));
		}
	      }
	    }
	    else {
	      if (cp[2] < 0.0) {
		if (!cut_flag[4]) {
		  cut_flag[4] = true;
		  cut_position[4] = max(x_no[nooed[ied][0]][2],x_no[nooed[ied][1]][2]);
		}
		else {
		  cut_position[4] = max(cut_position[4],max(x_no[nooed[ied][0]][2],x_no[nooed[ied][1]][2]));
		}
	      }
	      else {
		if (!cut_flag[5]) {
		  cut_flag[5] = true;
		  cut_position[5] = min(x_no[nooed[ied][0]][2],x_no[nooed[ied][1]][2]);
		}
		else {
		  cut_position[5] = min(cut_position[5],max(x_no[nooed[ied][0]][2],x_no[nooed[ied][1]][2]));
		}
	      }
	    }
	  }
	}
	delete[] fa_normal_;
	if (cut_flag[0]) {
	  const double x[3] = { cut_position[0]+1.0E-6*L, 0.0, 0.0 };
	  const double n[3] = { -0.5*L, 0.0, 0.0 }; // normal magnitude related to local length scale
	  cut_surf(x,n,-2);
	}
	if (cut_flag[1]) {
	  const double x[3] = { cut_position[1]-1.0E-6*L, 0.0, 0.0 };
	  const double n[3] = { 0.5*L, 0.0, 0.0 };
	  cut_surf(x,n,-3);
	}
	if (cut_flag[2]) {
	  const double x[3] = { 0.0, cut_position[2]+1.0E-6*L, 0.0 };
	  const double n[3] = { 0.0, -0.5*L, 0.0 };
	  cut_surf(x,n,-4);
	}
	if (cut_flag[3]) {
	  const double x[3] = { 0.0, cut_position[3]-1.0E-6*L, 0.0 };
	  const double n[3] = { 0.0, 0.5*L, 0.0 };
	  cut_surf(x,n,-5);
	}
	if (cut_flag[4]) {
	  const double x[3] = { 0.0, 0.0, cut_position[4]+1.0E-6*L };
	  const double n[3] = { 0.0, 0.0, -0.5*L };
	  cut_surf(x,n,-6);
	}
	if (cut_flag[5]) {
	  const double x[3] = { 0.0, 0.0, cut_position[5]-1.0E-6*L };
	  const double n[3] = { 0.0, 0.0, 0.5*L };
	  cut_surf(x,n,-7);
	}
      }
    }

    if (debug) {
      double xc[3] = { 0.0, 0.0, 0.0 };
      writeTecplot(1,xc);
    }

    int no_flag[nno];
    FOR_INO no_flag[ino] = -1;

    {
      FILE * fp = NULL;
      FOR_IED {
	if (faoed[ied][1] < 0) {
	  // if this is still a problem, then write a file...
	  if (fp == NULL) {
	    char filename[128];
	    sprintf(filename,"debug_cuttable.%06d.dat",mpi_rank);
	    fp = fopen(filename,"w");
	    const int ino0 = nooed[ied][0];
	    const int ino1 = nooed[ied][1];
	    fprintf(fp,"%18.15e %18.15e %18.15e\n",x_no[ino0][0],x_no[ino0][1],x_no[ino0][2]);
	    fprintf(fp,"%18.15e %18.15e %18.15e\n",x_no[ino1][0],x_no[ino1][1],x_no[ino1][2]);
	  }
	}
	assert((faoed[ied][0] >= -7)&&(faoed[ied][0] != -1)); // should all be gone -- i.e. cut away.
	if (faoed[ied][0] < -1) {
	  const int ic = -faoed[ied][0]-2;
	  assert((ic >= 0)&&(ic < 6)); // -x,+x,-y,+y,-z,+z
	  const int ino0 = nooed[ied][0]; assert(ino0 >= 0);
	  assert(no_flag[ino0] == -1);
	  no_flag[ino0] = ic;
	}
      }
      if (fp != NULL) {
	cout << "Error: cuttable problem not solved. see file debug_cuttable.dat. --DEBUG_RANK=" << mpi_rank << " --DEBUG_COUNT_CC=" << debug_count_cc << endl;
	fclose(fp);
	throw(-1);
      }
    }

    // now there are several possibilities for cutting cube...
    vector<CutCubeData> cutCubeDataVec;

    FOR_IED {
      if (faoed[ied][0] < -1) {
	const int ic0 = -faoed[ied][0]-2;
	assert((ic0 >= 0)&&(ic0 < 6));
	const int ino1 = nooed[ied][1]; assert(ino1 >= 0);
	const int ic1 = no_flag[ino1];
	assert((ic1 >= 0)&&(ic1 < 6));
	no_flag[ino1] = -1;
	if (ic0 != ic1) {

	  //cout << " ic0,ic1: " << ic0 << " " << ic1 << COUT_VEC(nodeVec[ino1].x) << endl;

	  /*
	    FILE * fp = fopen("p.dat","w");
	    fprintf(fp,"%18.15e %18.15e %18.15e\n",nodeVec[ino1].x[0]+xv[iv][0],nodeVec[ino1].x[1]+xv[iv][1],nodeVec[ino1].x[2]+xv[iv][2]);
	    fclose(fp);
	  */

	  // 4 z-edges x 2 signs...
	  // 0,2,1,3

	  if ((ic0 == 0)&&(ic1 == 2)) {
	    cutCubeDataVec.push_back(CutCubeData(0,1,ino1,x_no[ino1][2]));
	  }
	  else if ((ic0 == 2)&&(ic1 == 0)) {
	    cutCubeDataVec.push_back(CutCubeData(0,-1,ino1,x_no[ino1][2]));
	  }
	  else if ((ic0 == 2)&&(ic1 == 1)) {
	    cutCubeDataVec.push_back(CutCubeData(1,1,ino1,x_no[ino1][2]));
	  }
	  else if ((ic0 == 1)&&(ic1 == 2)) {
	    cutCubeDataVec.push_back(CutCubeData(1,-1,ino1,x_no[ino1][2]));
	  }
	  else if ((ic0 == 1)&&(ic1 == 3)) {
	    cutCubeDataVec.push_back(CutCubeData(2,1,ino1,x_no[ino1][2]));
	  }
	  else if ((ic0 == 3)&&(ic1 == 1)) {
	    cutCubeDataVec.push_back(CutCubeData(2,-1,ino1,x_no[ino1][2]));
	  }
	  else if ((ic0 == 3)&&(ic1 == 0)) {
	    cutCubeDataVec.push_back(CutCubeData(3,1,ino1,x_no[ino1][2]));
	  }
	  else if ((ic0 == 0)&&(ic1 == 3)) {
	    cutCubeDataVec.push_back(CutCubeData(3,-1,ino1,x_no[ino1][2]));
	  }

	  // 4 x-edges x 2 signs...
	  // 2,4,3,5

	  else if ((ic0 == 2)&&(ic1 == 4)) {
	    cutCubeDataVec.push_back(CutCubeData(4,1,ino1,x_no[ino1][0]));
	  }
	  else if ((ic0 == 4)&&(ic1 == 2)) {
	    cutCubeDataVec.push_back(CutCubeData(4,-1,ino1,x_no[ino1][0]));
	  }
	  else if ((ic0 == 4)&&(ic1 == 3)) {
	    cutCubeDataVec.push_back(CutCubeData(5,1,ino1,x_no[ino1][0]));
	  }
	  else if ((ic0 == 3)&&(ic1 == 4)) {
	    cutCubeDataVec.push_back(CutCubeData(5,-1,ino1,x_no[ino1][0]));
	  }
	  else if ((ic0 == 3)&&(ic1 == 5)) {
	    cutCubeDataVec.push_back(CutCubeData(6,1,ino1,x_no[ino1][0]));
	  }
	  else if ((ic0 == 5)&&(ic1 == 3)) {
	    cutCubeDataVec.push_back(CutCubeData(6,-1,ino1,x_no[ino1][0]));
	  }
	  else if ((ic0 == 5)&&(ic1 == 2)) {
	    cutCubeDataVec.push_back(CutCubeData(7,1,ino1,x_no[ino1][0]));
	  }
	  else if ((ic0 == 2)&&(ic1 == 5)) {
	    cutCubeDataVec.push_back(CutCubeData(7,-1,ino1,x_no[ino1][0]));
	  }

	  // 4 y-edges x 2 signs...
	  // 4,0,5,1

	  else if ((ic0 == 4)&&(ic1 == 0)) {
	    cutCubeDataVec.push_back(CutCubeData(8,1,ino1,x_no[ino1][1]));
	  }
	  else if ((ic0 == 0)&&(ic1 == 4)) {
	    cutCubeDataVec.push_back(CutCubeData(8,-1,ino1,x_no[ino1][1]));
	  }
	  else if ((ic0 == 0)&&(ic1 == 5)) {
	    cutCubeDataVec.push_back(CutCubeData(9,1,ino1,x_no[ino1][1]));
	  }
	  else if ((ic0 == 5)&&(ic1 == 0)) {
	    cutCubeDataVec.push_back(CutCubeData(9,-1,ino1,x_no[ino1][1]));
	  }
	  else if ((ic0 == 5)&&(ic1 == 1)) {
	    cutCubeDataVec.push_back(CutCubeData(10,1,ino1,x_no[ino1][1]));
	  }
	  else if ((ic0 == 1)&&(ic1 == 5)) {
	    cutCubeDataVec.push_back(CutCubeData(10,-1,ino1,x_no[ino1][1]));
	  }
	  else if ((ic0 == 1)&&(ic1 == 4)) {
	    cutCubeDataVec.push_back(CutCubeData(11,1,ino1,x_no[ino1][1]));
	  }
	  else if ((ic0 == 4)&&(ic1 == 1)) {
	    cutCubeDataVec.push_back(CutCubeData(11,-1,ino1,x_no[ino1][1]));
	  }

	  // any other edge is passing through a corner and can be avoided
	  // by slightly changing delta...

	  else {
	    cout << "Problem. Do not know how to handle ic0,ic1: " << ic0 << " " << ic1 << endl;
	    assert(0);
	  }

	}
      }
    }

    // ====================================
    // ====================================
    // ====================================

    if (cutCubeDataVec.empty()) {

      // none of the cube edges were intersected by the surface, so we need to
      // check if a region of the surface is cut by a single face...

      // try to figure out if the cut surfaces associated with any/all cartesian
      // faces form positive or negative areas...

      double face_area[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
      bool got_edge = false;

      FOR_IED {
	if (faoed[ied][0] < -1) {
	  // this edge will be enough for us to determine the surface orientation...
	  got_edge = true;
	  const int ic = -faoed[ied][0]-2;
	  assert((ic >= 0)&&(ic < 6)); // -x,+x,-y,+y,-z,+z
	  const int ino0 = nooed[ied][0]; assert(ino0 >= 0);
	  const int ino1 = nooed[ied][1]; assert(ino1 >= 0);
	  switch (ic) {
	  case 0:
	    // surface on -x...
	    assert( fabs(x_no[ino0][0] + L) < 1.0E-10 );
	    assert( fabs(x_no[ino1][0] + L) < 1.0E-10 );
	    face_area[0] += x_no[ino0][1]*x_no[ino1][2] - x_no[ino0][2]*x_no[ino1][1];
	    break;
	  case 1:
	    // surface on +x...
	    assert( fabs(x_no[ino0][0] - L) < 1.0E-10 );
	    assert( fabs(x_no[ino1][0] - L) < 1.0E-10 );
	    face_area[1] -= x_no[ino0][1]*x_no[ino1][2] - x_no[ino0][2]*x_no[ino1][1];
	    break;
	  case 2:
	    // surface on -y...
	    assert( fabs(x_no[ino0][1] + L) < 1.0E-10 );
	    assert( fabs(x_no[ino1][1] + L) < 1.0E-10 );
	    face_area[2] += x_no[ino0][2]*x_no[ino1][0] - x_no[ino0][0]*x_no[ino1][2];
	    break;
	  case 3:
	    // surface on +y...
	    assert( fabs(x_no[ino0][1] - L) < 1.0E-10 );
	    assert( fabs(x_no[ino1][1] - L) < 1.0E-10 );
	    face_area[3] -= x_no[ino0][2]*x_no[ino1][0] - x_no[ino0][0]*x_no[ino1][2];
	    break;
	  case 4:
	    // surface on -z...
	    assert( fabs(x_no[ino0][2] + L) < 1.0E-10 );
	    assert( fabs(x_no[ino1][2] + L) < 1.0E-10 );
	    face_area[4] += x_no[ino0][0]*x_no[ino1][1] - x_no[ino0][1]*x_no[ino1][0];
	    break;
	  case 5:
	    // surface on +z...
	    assert( fabs(x_no[ino0][2] - L) < 1.0E-10 );
	    assert( fabs(x_no[ino1][2] - L) < 1.0E-10 );
	    face_area[5] -= x_no[ino0][0]*x_no[ino1][1] - x_no[ino0][1]*x_no[ino1][0];
	    break;
	  default:
	    assert(0);
	  }
	}
      }

      if (!got_edge) {

	// let this go...
	//cout << "Warning: surface is entirely inside starting volume (floating part of the mesh? or perhaps a very coarse 1-cell grid)" << endl;
	
	addCube(L);

      }
      else if ( (face_area[0] >= 0.0) && (face_area[1] >= 0.0) && (face_area[2] >= 0.0) &&
		(face_area[3] >= 0.0) && (face_area[4] >= 0.0) && (face_area[5] >= 0.0) ) {

	// the surface is penetrating into the cube. We need to add the full outer cube...

	addCube(L);

      }
      else {

	// make sure all areas have the expected sign...

	if (!( (face_area[0] <= 0.0) && (face_area[1] <= 0.0) && (face_area[2] <= 0.0) &&
	       (face_area[3] <= 0.0) && (face_area[4] <= 0.0) && (face_area[5] <= 0.0) ) ) {
	  cout << "Error: CuttableVoronoiData: failed face_area check: " << endl;
	  throw(-1);
	}

	//assert( (face_area[0] <= 0.0) && (face_area[1] <= 0.0) && (face_area[2] <= 0.0) &&
	//(face_area[3] <= 0.0) && (face_area[4] <= 0.0) && (face_area[5] <= 0.0) );

      }

    }
    else {

      // !CutCubeDataVec.empty()
      // use the intersections along edges to draw the corect part of the cube...

      vector< pair<int,int> > corner_flag;
      corner_flag.resize(8);
      for (int i=0; i<8; ++i) {
        corner_flag[i].first = -1;  // flag value
        corner_flag[i].second = 0;  // count of edges touching this corner (checks consistency of edges at cube corners)
      }
      bool edge_flag[12] = { false, false, false, false,
			     false, false, false, false,
			     false, false, false, false };

      std::sort(cutCubeDataVec.begin(),cutCubeDataVec.end());

      int current_id = -1;
      int current_sign = 0; // will be 1 or -1 during edge loop...
      int current_ino = -1;
      for (int ii = 0,ii_end=cutCubeDataVec.size(); ii < ii_end; ++ii) {
	if (cutCubeDataVec[ii].id != current_id) {
	  // we are moving to an new edge. if the last intersection of the previous edge had a positive sign, then
	  // add the completing edge from that node to the relevant corner...
	  if (current_sign == 1) {
#include "completeEdge.hpp"
	  }
	  current_id = cutCubeDataVec[ii].id;
	  current_sign = cutCubeDataVec[ii].sign;
	  current_ino = cutCubeDataVec[ii].ino;
	  if (current_sign == -1) {
#include "startEdge.hpp"
	  }
	  // record that we worked on this edge...
	  edge_flag[current_id] = true;
	}
	else {
	  // this is another node along the same edge. If the current sign
	  // is 1, then add an edge from current_ino to our ino...
	  if (!(cutCubeDataVec[ii].sign == -current_sign)) {
	    // cout << "Error: !(cutCubeDataVec[ii].sign == -current_sign)" << cutCubeDataVec[ii].sign << " " << current_sign << endl;
	    return -1; // should fail check...
	  }
	  assert(cutCubeDataVec[ii].sign == -current_sign);
	  if (current_sign == 1) {
	    const int next_ino = cutCubeDataVec[ii].ino;
#include "addEdge.hpp"
	  }
	  current_sign = cutCubeDataVec[ii].sign;
	  current_ino = cutCubeDataVec[ii].ino;
	}
      }
      if (current_sign == 1) {
#include "completeEdge.hpp"
      }

      // now complete the hull by ensuring all active corners have their
      // edges...

      bool done = false;
      while (!done) {
	done = true;

	// ===============================
	// z-edges
	// ===============================

#define IED  0
#define IFA0 -2
#define IFA1 -4
#define INO0 0
#define INO1 4
#define XNO0 (-L)
#define YNO0 (-L)
#define ZNO0 (-L)
#define XNO1 (-L)
#define YNO1 (-L)
#define ZNO1 (L)
#include "checkOrSetEdge.hpp"
#undef IED
#undef IFA0
#undef IFA1
#undef INO0
#undef INO1
#undef XNO0
#undef YNO0
#undef ZNO0
#undef XNO1
#undef YNO1
#undef ZNO1

#define IED  1
#define IFA0 -4
#define IFA1 -3
#define INO0 1
#define INO1 5
#define XNO0 (L)
#define YNO0 (-L)
#define ZNO0 (-L)
#define XNO1 (L)
#define YNO1 (-L)
#define ZNO1 (L)
#include "checkOrSetEdge.hpp"
#undef IED
#undef IFA0
#undef IFA1
#undef INO0
#undef INO1
#undef XNO0
#undef YNO0
#undef ZNO0
#undef XNO1
#undef YNO1
#undef ZNO1

#define IED  2
#define IFA0 -3
#define IFA1 -5
#define INO0 2
#define INO1 6
#define XNO0 (L)
#define YNO0 (L)
#define ZNO0 (-L)
#define XNO1 (L)
#define YNO1 (L)
#define ZNO1 (L)
#include "checkOrSetEdge.hpp"
#undef IED
#undef IFA0
#undef IFA1
#undef INO0
#undef INO1
#undef XNO0
#undef YNO0
#undef ZNO0
#undef XNO1
#undef YNO1
#undef ZNO1

#define IED  3
#define IFA0 -5
#define IFA1 -2
#define INO0 3
#define INO1 7
#define XNO0 (-L)
#define YNO0 (L)
#define ZNO0 (-L)
#define XNO1 (-L)
#define YNO1 (L)
#define ZNO1 (L)
#include "checkOrSetEdge.hpp"
#undef IED
#undef IFA0
#undef IFA1
#undef INO0
#undef INO1
#undef XNO0
#undef YNO0
#undef ZNO0
#undef XNO1
#undef YNO1
#undef ZNO1

	// ===============================
	// x-edges
	// ===============================

#define IED  4
#define IFA0 -4
#define IFA1 -6
#define INO0 0
#define INO1 1
#define XNO0 (-L)
#define YNO0 (-L)
#define ZNO0 (-L)
#define XNO1 (L)
#define YNO1 (-L)
#define ZNO1 (-L)
#include "checkOrSetEdge.hpp"
#undef IED
#undef IFA0
#undef IFA1
#undef INO0
#undef INO1
#undef XNO0
#undef YNO0
#undef ZNO0
#undef XNO1
#undef YNO1
#undef ZNO1

#define IED  5
#define IFA0 -6
#define IFA1 -5
#define INO0 3
#define INO1 2
#define XNO0 (-L)
#define YNO0 (L)
#define ZNO0 (-L)
#define XNO1 (L)
#define YNO1 (L)
#define ZNO1 (-L)
#include "checkOrSetEdge.hpp"
#undef IED
#undef IFA0
#undef IFA1
#undef INO0
#undef INO1
#undef XNO0
#undef YNO0
#undef ZNO0
#undef XNO1
#undef YNO1
#undef ZNO1

#define IED  6
#define IFA0 -5
#define IFA1 -7
#define INO0 7
#define INO1 6
#define XNO0 (-L)
#define YNO0 (L)
#define ZNO0 (L)
#define XNO1 (L)
#define YNO1 (L)
#define ZNO1 (L)
#include "checkOrSetEdge.hpp"
#undef IED
#undef IFA0
#undef IFA1
#undef INO0
#undef INO1
#undef XNO0
#undef YNO0
#undef ZNO0
#undef XNO1
#undef YNO1
#undef ZNO1

#define IED  7
#define IFA0 -7
#define IFA1 -4
#define INO0 4
#define INO1 5
#define XNO0 (-L)
#define YNO0 (-L)
#define ZNO0 (L)
#define XNO1 (L)
#define YNO1 (-L)
#define ZNO1 (L)
#include "checkOrSetEdge.hpp"
#undef IED
#undef IFA0
#undef IFA1
#undef INO0
#undef INO1
#undef XNO0
#undef YNO0
#undef ZNO0
#undef XNO1
#undef YNO1
#undef ZNO1

	// ===============================
	// y-edges
	// ===============================

#define IED  8
#define IFA0 -6
#define IFA1 -2
#define INO0 0
#define INO1 3
#define XNO0 (-L)
#define YNO0 (-L)
#define ZNO0 (-L)
#define XNO1 (-L)
#define YNO1 (L)
#define ZNO1 (-L)
#include "checkOrSetEdge.hpp"
#undef IED
#undef IFA0
#undef IFA1
#undef INO0
#undef INO1
#undef XNO0
#undef YNO0
#undef ZNO0
#undef XNO1
#undef YNO1
#undef ZNO1

#define IED  9
#define IFA0 -2
#define IFA1 -7
#define INO0 4
#define INO1 7
#define XNO0 (-L)
#define YNO0 (-L)
#define ZNO0 (L)
#define XNO1 (-L)
#define YNO1 (L)
#define ZNO1 (L)
#include "checkOrSetEdge.hpp"
#undef IED
#undef IFA0
#undef IFA1
#undef INO0
#undef INO1
#undef XNO0
#undef YNO0
#undef ZNO0
#undef XNO1
#undef YNO1
#undef ZNO1

#define IED  10
#define IFA0 -7
#define IFA1 -3
#define INO0 5
#define INO1 6
#define XNO0 (L)
#define YNO0 (-L)
#define ZNO0 (L)
#define XNO1 (L)
#define YNO1 (L)
#define ZNO1 (L)
#include "checkOrSetEdge.hpp"
#undef IED
#undef IFA0
#undef IFA1
#undef INO0
#undef INO1
#undef XNO0
#undef YNO0
#undef ZNO0
#undef XNO1
#undef YNO1
#undef ZNO1

#define IED  11
#define IFA0 -3
#define IFA1 -6
#define INO0 1
#define INO1 2
#define XNO0 (L)
#define YNO0 (-L)
#define ZNO0 (-L)
#define XNO1 (L)
#define YNO1 (L)
#define ZNO1 (-L)
#include "checkOrSetEdge.hpp"
#undef IED
#undef IFA0
#undef IFA1
#undef INO0
#undef INO1
#undef XNO0
#undef YNO0
#undef ZNO0
#undef XNO1
#undef YNO1
#undef ZNO1

      }

      // check consistency at cube edges
      bool valid = true;
      for (int i=0; i<8; ++i) {
        // cout << "corner[" << i << "] count: " << corner_flag[i].second << endl;
        if ((corner_flag[i].second != 0) && (corner_flag[i].second != 3)) {
          valid = false;
        }
      }
      if (!valid) return 5;
    }

    return 0;
  }

  int completeCube2(const double Lx,const double Ly,const double Lz,const bool debug);
  
  int check() const {

    // HACK - skip check...
    //return;

    if ((ned <= 0)||(nno <= 0)) {
      cout << "FAILED CuttableVoronoiData::check(): ned,nno problem: " << ned << " " << nno << endl;
      throw(-1);
    }

    // look for nans in node locations...
    FOR_INO {
      FOR_I3
      if (x_no[ino][i] != x_no[ino][i]) {
	cout << "FAILED CuttableVoronoiData::check(): x_no nan problem" << endl;
	throw(-1);
      }
    }

    int * no_flag = new int[nno];
    FOR_INO no_flag[ino] = 0;
    FOR_IED {
      if (faoed[ied][0] == -1) {
	cout << "FAILED CuttableVoronoiData::check(): faoed[ied][0] == -1" << endl;
	throw(-1);
      }
      FOR_I2 {
	const int ino = nooed[ied][i];
	assert((ino >= 0)&&(ino < nno));
	no_flag[ino] = 1;
      }
    }
    FOR_INO {
      if (no_flag[ino] != 1) {
	cout << "FAILED CuttableVoronoiData::check(): some nodes present that do not touch edges" << endl;
	throw(-1);
      }
    }
    delete[] no_flag;

    set<int> faceSet;
    FOR_IED {
      FOR_I2 {
	// faoed can be just about anything:
	// +ve value here represent surface tris, -ve values represet either
	// original hull faces (-2..-7 inclusive), or a voronoi nbr, -8,-9,etc...
	const int ifa_value = faoed[ied][i];
	faceSet.insert(ifa_value);
      }
    }
    const int nfa = faceSet.size();
    double (*face_normal)[3] = new double[nfa][3];
    FOR_IFA FOR_I3 face_normal[ifa][i] = 0.0;
    double (*face_dx)[3] = new double[nfa][3];
    FOR_IFA FOR_I3 face_dx[ifa][i] = 0.0;
    double (*x_fa)[3] = new double[nfa][3];
    FOR_IFA FOR_I3 x_fa[ifa][i] = 0.0;
    int * fa_flag = new int[nfa];
    FOR_IFA fa_flag[ifa] = 0;
    map<const int,int> faceMap;
    {
      int ifa = 0;
      for (set<int>::iterator iter = faceSet.begin(); iter != faceSet.end(); ++iter)
	faceMap[*iter] = ifa++;
      assert(ifa == nfa);
    }

    double d2_max = 0.0;
    FOR_IED {
      assert(faoed[ied][0] != faoed[ied][1]);
      const int ino0 = nooed[ied][0];
      const int ino1 = nooed[ied][1];
      const double n[3] = CROSS_PRODUCT(x_no[ino0],x_no[ino1]);
      const double dx[3] = DIFF(x_no[ino1],x_no[ino0]);
      d2_max = max(d2_max,DOT_PRODUCT(dx,dx));
      const double x_mid[3] = {
	0.5*(x_no[ino1][0]+x_no[ino0][0]),
	0.5*(x_no[ino1][1]+x_no[ino0][1]),
	0.5*(x_no[ino1][2]+x_no[ino0][2]) };
      // ifa0...
      {
	map<const int,int>::iterator iter = faceMap.find(faoed[ied][0]);
	assert(iter != faceMap.end());
	const int ifa = iter->second;
	FOR_I3 face_normal[ifa][i] += n[i];
	FOR_I3 face_dx[ifa][i] += dx[i];
	FOR_I3 x_fa[ifa][i] += x_mid[i];
	fa_flag[ifa] += 1;
      }
      // ifa1...
      {
	map<const int,int>::iterator iter = faceMap.find(faoed[ied][1]);
	assert(iter != faceMap.end());
	const int ifa = iter->second;
	FOR_I3 face_normal[ifa][i] -= n[i];
	FOR_I3 face_dx[ifa][i] -= dx[i];
	FOR_I3 x_fa[ifa][i] += x_mid[i];
	fa_flag[ifa] += 1;
      }
      // notice that by construction, this volume will satisfy gcl...
    }
    
    int ierr = 0;
    FOR_IFA {
      //cout << "DOT_PRODUCT(face_dx[ifa],face_dx[ifa]): " << DOT_PRODUCT(face_dx[ifa],face_dx[ifa]) << endl;
      if (!(DOT_PRODUCT(face_dx[ifa],face_dx[ifa]) < 1.0E-18*d2_max)) {
	cout << "ifa: " << ifa << " DOT_PRODUCT(face_dx[ifa],face_dx[ifa])/d2_max: " << DOT_PRODUCT(face_dx[ifa],face_dx[ifa])/d2_max << endl;
        bool report = true;
	// report associated nodes...
        FOR_IED {
          const int ino0 = nooed[ied][0];
          const int ino1 = nooed[ied][1];
          // ifa0...
          {
            map<const int,int>::iterator iter = faceMap.find(faoed[ied][0]);
            assert(iter != faceMap.end());
            if (ifa == iter->second) {
	      if (report) {
		cout << "faoed: " << faoed[ied][0] << endl;
		report = false;
	      }
              // this is a participating edge...
              cout << " EDGE: " << x_no[ino0][0] << " " << x_no[ino0][1] << " " << x_no[ino0][2] << endl;
              cout << " MID: " << 
                0.5*(x_no[ino0][0]+x_no[ino1][0]) << " " <<
                0.5*(x_no[ino0][1]+x_no[ino1][1]) << " " <<
                0.5*(x_no[ino0][2]+x_no[ino1][2]) << endl;
              cout << " EDGE: " << x_no[ino1][0] << " " << x_no[ino1][1] << " " << x_no[ino1][2] << endl;
            }
          }
          // ifa1...
          {
            map<const int,int>::iterator iter = faceMap.find(faoed[ied][1]);
            assert(iter != faceMap.end());
            if (ifa == iter->second) {
	      if (report) {
		cout << "faoed: " << faoed[ied][1] << endl;
		report = false;
	      }
              // this is a participating edge...
              cout << " EDGE: " << x_no[ino1][0] << " " << x_no[ino1][1] << " " << x_no[ino1][2] << endl;
              cout << " MID: " << 
                0.5*(x_no[ino0][0]+x_no[ino1][0]) << " " <<
                0.5*(x_no[ino0][1]+x_no[ino1][1]) << " " <<
                0.5*(x_no[ino0][2]+x_no[ino1][2]) << endl;
              cout << " EDGE: " << x_no[ino0][0] << " " << x_no[ino0][1] << " " << x_no[ino0][2] << endl;
            }
          }
        }
	//writeTecplot(0,x0);
	ierr = -1;
        //throw(-1);
      }
      FOR_I3 x_fa[ifa][i] /= double(fa_flag[ifa]);
      // actually, this dp can be quite negative...
      //cout << "dp: " << DOT_PRODUCT(x_fa[ifa],face_normal[ifa]) << endl;
      //if (DOT_PRODUCT(x_fa[ifa],face_normal[ifa]) < 0.0)
      //  throw(0);
    }

    if (ierr != 0) {
      cout << "FAILED CuttableVoronoiData::check(): non-zero faces found" << endl;
      throw(-1);
    }

    delete[] face_normal;
    delete[] face_dx;
    delete[] x_fa;
    delete[] fa_flag;

    //assert(ierr == 0);
    return ierr;

  }

  void check(set<int>& stSet,const double x0[3]) const {

    static int debug_index = 0;

    // HACK - skip check...
    //return;

    if ((ned <= 0)||(nno <= 0)) {
      cout << "ned,nno problem: " << ned << " " << nno << endl;
      throw(-1);
    }

    // look for nans in node locations...
    FOR_INO
      FOR_I3
      if (x_no[ino][i] != x_no[ino][i]) {
	cout << "x_no nan problem" << endl;
	throw(-1);
      }

    int * no_flag = new int[nno];
    FOR_INO no_flag[ino] = 0;
    FOR_IED {
      if ((faoed[ied][0] == -1)||(faoed[ied][1] == -1)) {
	cout << "Warning: faoed -1: " << faoed[ied][0] << " " << faoed[ied][1] << endl;
	//throw(-1);
      }
      FOR_I2 {
	const int ino = nooed[ied][i];
	no_flag[ino] = 1;
      }
    }
    FOR_INO assert(no_flag[ino] == 1);
    delete[] no_flag;

    set<int> faceSet;
    FOR_IED {
      FOR_I2 {
	// faoed can be just about anything:
	// +ve value here represent surface tris, -ve values represet either
	// original hull faces (-2..-7 inclusive), or a voronoi nbr, -8,-9,etc...
	const int ifa_value = faoed[ied][i];
	faceSet.insert(ifa_value);
      }
    }
    const int nfa = faceSet.size();
    double (*face_normal)[3] = new double[nfa][3];
    FOR_IFA FOR_I3 face_normal[ifa][i] = 0.0;
    double (*face_dx)[3] = new double[nfa][3];
    FOR_IFA FOR_I3 face_dx[ifa][i] = 0.0;
    double (*x_fa)[3] = new double[nfa][3];
    FOR_IFA FOR_I3 x_fa[ifa][i] = 0.0;
    int * fa_flag = new int[nfa];
    FOR_IFA fa_flag[ifa] = 0;

    map<const int,int> faceMap;
    {
      int ifa = 0;
      for (set<int>::iterator iter = faceSet.begin(); iter != faceSet.end(); ++iter)
	faceMap[*iter] = ifa++;
      assert(ifa == nfa);
    }

    double d2_max = 0.0;
    FOR_IED {
      assert(faoed[ied][0] != faoed[ied][1]);
      const int ino0 = nooed[ied][0];
      const int ino1 = nooed[ied][1];
      const double n[3] = CROSS_PRODUCT(x_no[ino0],x_no[ino1]);
      const double dx[3] = DIFF(x_no[ino1],x_no[ino0]);
      d2_max = max(d2_max,DOT_PRODUCT(dx,dx));
      const double x_mid[3] = {
	0.5*(x_no[ino1][0]+x_no[ino0][0]),
	0.5*(x_no[ino1][1]+x_no[ino0][1]),
	0.5*(x_no[ino1][2]+x_no[ino0][2]) };
      // ifa0...
      {
	map<const int,int>::iterator iter = faceMap.find(faoed[ied][0]);
	assert(iter != faceMap.end());
	const int ifa = iter->second;
	FOR_I3 face_normal[ifa][i] += n[i];
	FOR_I3 face_dx[ifa][i] += dx[i];
	FOR_I3 x_fa[ifa][i] += x_mid[i];
	fa_flag[ifa] += 1;
      }
      // ifa1...
      {
	map<const int,int>::iterator iter = faceMap.find(faoed[ied][1]);
	assert(iter != faceMap.end());
	const int ifa = iter->second;
	FOR_I3 face_normal[ifa][i] -= n[i];
	FOR_I3 face_dx[ifa][i] -= dx[i];
	FOR_I3 x_fa[ifa][i] += x_mid[i];
	fa_flag[ifa] += 1;
      }
      // notice that by construction, this volume will satisfy gcl...
    }

    FOR_IFA {
      //cout << "DOT_PRODUCT(face_dx[ifa],face_dx[ifa]): " << DOT_PRODUCT(face_dx[ifa],face_dx[ifa]) << endl;
      if (!(DOT_PRODUCT(face_dx[ifa],face_dx[ifa]) < 1.0E-18*d2_max)) {
	cout << "DOT_PRODUCT(face_dx[ifa],face_dx[ifa])/d2_max: " << DOT_PRODUCT(face_dx[ifa],face_dx[ifa])/d2_max << endl;
	FOR_IED FOR_I2 if (faoed[ied][i] >= 0) stSet.insert(faoed[ied][i]);
	writeTecplot(debug_index++,x0);
	//throw(-1);
      }
      FOR_I3 x_fa[ifa][i] /= double(fa_flag[ifa]);
      // actually, this dp can be quite negative...
      //cout << "dp: " << DOT_PRODUCT(x_fa[ifa],face_normal[ifa]) << endl;
      //if (DOT_PRODUCT(x_fa[ifa],face_normal[ifa]) < 0.0)
      //  throw(0);
    }
    //getchar();

    delete[] face_normal;
    delete[] face_dx;
    delete[] x_fa;
    delete[] fa_flag;

  }
  
};

#endif
