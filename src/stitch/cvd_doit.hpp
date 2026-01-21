
// there should be a cut surface at this point, with nooed referring
// to the local nodes, and faoed as follows:
// faoed[ied][0,1] >= 0: the index of the associated surface tri
// faoed[ied][0] == -2 : x0 face
// faoed[ied][0] == -3 : x1 face
// faoed[ied][0] == -4 : y0 face
// faoed[ied][0] == -5 : y1 face
// faoed[ied][0] == -6 : z0 face
// faoed[ied][0] == -7 : z1 face
// faoed[ied][0] == -1 : original edge of the surface tri -- should all be cut away

//#define X0_FACE -2
//#define X1_FACE -3
//#define Y0_FACE -4
//#define Y1_FACE -5
//#define Z0_FACE -6
//#define Z1_FACE -7

// 6 times the tet volume of 3 integer points, cast to int8...
#define SIGNED_TET_VOLUME_6_INT8(A,B,C,D) (int8((B)[0]-(A)[0])*(int8((C)[1]-(A)[1])*int8((D)[2]-(A)[2])-int8((C)[2]-(A)[2])*int8((D)[1]-(A)[1])) + \
                                           int8((B)[1]-(A)[1])*(int8((C)[2]-(A)[2])*int8((D)[0]-(A)[0])-int8((C)[0]-(A)[0])*int8((D)[2]-(A)[2])) + \
                                           int8((B)[2]-(A)[2])*(int8((C)[0]-(A)[0])*int8((D)[1]-(A)[1])-int8((C)[1]-(A)[1])*int8((D)[0]-(A)[0])))

// twice the tri normal of 3 integer points, cast to int8...
#define TRI_NORMAL_2_INT8(A,B,C) { int8((B)[1]-(A)[1])*int8((C)[2]-(A)[2])-int8((B)[2]-(A)[2])*int8((C)[1]-(A)[1]), \
                                   int8((B)[2]-(A)[2])*int8((C)[0]-(A)[0])-int8((B)[0]-(A)[0])*int8((C)[2]-(A)[2]), \
                                   int8((B)[0]-(A)[0])*int8((C)[1]-(A)[1])-int8((B)[1]-(A)[1])*int8((C)[0]-(A)[0]) }

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

/*
  #define STACK_BIT 1
  #define IN_BIT 2
  #define OUT_BIT 4
  #define LINK_BIT 8
*/

// ============================================================================================
// start of new routines
// ============================================================================================

class Bbox {
public:

  double xmin[3];
  double xmax[3];
  
  Bbox() {
    FOR_I3 xmin[i] = HUGE_VAL;
    FOR_I3 xmax[i] = -HUGE_VAL;
  }

  Bbox(const double x0[3],const double x1[3]) {
    FOR_I3 xmin[i] = min(x0[i],x1[i]);
    FOR_I3 xmax[i] = max(x0[i],x1[i]);
  }
  
  void expand(const double x[3]) {
    FOR_I3 xmin[i] = min(xmin[i],x[i]);
    FOR_I3 xmax[i] = max(xmax[i],x[i]);
  }
  
  void getBbox(double bbmin[3],double bbmax[3]) {
    FOR_I3 bbmin[i] = xmin[i];
    FOR_I3 bbmax[i] = xmax[i];
  }
  
};

/*
void dump_edbits(const int ed_flag) {

  const int in_bit        = (1<<0);
  const int out_bit       = (1<<1);
  const int interface_bit = (1<<2);
  const int current_group_bit = (1<<3);
  const int grouped_bit      = (1<<4);

  cout << "[";
  bool first = true;
  if (ed_flag&in_bit) {
    assert(first);
    cout << "in";
    first = false;
  }
  if (ed_flag&out_bit) {
    assert(first);
    cout << "out";
    first = false;
  }
  if (ed_flag&interface_bit) {
    assert(!first);
    cout << ":interface";
  }
  if (ed_flag&current_group_bit) {
    if (first) {
      cout << "current_group";
      first = false;
    }
    else cout << ":current_group";
  }
  if (ed_flag&grouped_bit) {
    if (first) {
      cout << "grouped";
      first = false;
    }
    else cout << ":grouped";
  }  
  if (first) cout << "none";
  cout << "]";
  
}
*/

int CuttableVoronoiData::completeCube2(const double Lx,const double Ly,const double Lz,const bool debug) {
  
  // if empty, just add a cube of size L...
  
  if (empty()) {
    addCube(Lx,Ly,Lz);
    return 0;
  }
  
  const int nno_orig = nno;
  const int ned_orig = ned;
  
  int * no_flag = new int[nno+8]; // add 8 for use at corner nodes...
  for (int ino = 0; ino < nno+8; ++ino)
    no_flag[ino] = 0;

  FOR_IED {
    assert(faoed[ied][1] >= -1); 
    //if ((faoed[ied][0] >= -7)&&(faoed[ied][0] <= -2)) {
    if (faoed[ied][0] <= -2) {
      assert(faoed[ied][0] >= -7);
      // this is an edge on one of the cube faces. Set the corresponding bit in both its nodes:
      // here we add 6 to one of the nodes to store orientation information implicitly. This 
      // prevents us from having to compute normals along the cube corners to determine
      // the surface orientation...
      no_flag[nooed[ied][0]] |= (1<<(-faoed[ied][0]-2+6));
      no_flag[nooed[ied][1]] |= (1<<(-faoed[ied][0]-2));
    }
  }
  
  // check that all nodes are at expected locations. The combination of 
  // cutting against the normal with a single component (i.e x,y,or z) and
  // the choice of frac that avoids cutting through existing nodes in the surface
  // should force this to be exactly correct...
  FOR_INO {
    // x...
    if (no_flag[ino]&((1<<(-X0_FACE-2))|(1<<(-X0_FACE-2+6)))) {
      assert(x_no[ino][0] == -Lx);
    }
    else if (no_flag[ino]&((1<<(-X1_FACE-2))|(1<<(-X1_FACE-2+6)))) {
      assert(x_no[ino][0] == Lx);
    }
    else {
      assert((x_no[ino][0] > -Lx)&&(x_no[ino][0] < Lx)); 
    }
    // y...
    if (no_flag[ino]&((1<<(-Y0_FACE-2))|(1<<(-Y0_FACE-2+6)))) {
      assert(x_no[ino][1] == -Ly);
    }
    else if (no_flag[ino]&((1<<(-Y1_FACE-2))|(1<<(-Y1_FACE-2+6)))) {
      assert(x_no[ino][1] == Ly);
    }
    else {
      assert((x_no[ino][1] > -Ly)&&(x_no[ino][1] < Ly));
    }
    // z...
    if (no_flag[ino]&((1<<(-Z0_FACE-2))|(1<<(-Z0_FACE-2+6)))) {
      assert(x_no[ino][2] == -Lz);
    }
    else if (no_flag[ino]&((1<<(-Z1_FACE-2))|(1<<(-Z1_FACE-2+6)))) {
      assert(x_no[ino][2] == Lz);
    }
    else {
      assert((x_no[ino][2] > -Lz)&&(x_no[ino][2] < Lz));
    }
  }
    
  // store surface-cube edge intersections in a vec, one for each of the 
  // 12 cube edges. We use a vec here so we can sort along the edge if there
  // are multiple intersections...
  // edge 0 faces x0,y0 nodes 0,4
  // edge 1 faces y1,x0 nodes 2,6
  // edge 2 faces z0,x0 nodes 0,2
  // edge 3 faces x0,z1 nodes 4,6
  // edge 4 faces y0,x1 nodes 1,5
  // edge 5 faces x1,y1 nodes 3,7
  // edge 6 faces x1,z0 nodes 1,3
  // edge 7 faces z1,x1 nodes 5,7
  // edge 8 faces y0,z0 nodes 0,1
  // edge 9 faces z1,y0 nodes 4,5
  // edge 10 faces z0,y1 nodes 2,3
  // edge 11 faces y1,z1 nodes 6,7

  vector<pair<double,int> > cutCubeNodeVec[12];
  bool got_ccn = false;
    
  for (int ino = 0; ino < nno; ++ino) {
    if (no_flag[ino] == ((1<<(-X0_FACE-2))|(1<<(-Y0_FACE-2+6)))) {
      cutCubeNodeVec[0].push_back(pair<double,int>(x_no[ino][2],ino));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-X0_FACE-2+6))|(1<<(-Y0_FACE-2)))) {
      cutCubeNodeVec[0].push_back(pair<double,int>(x_no[ino][2],-ino-1));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-Y1_FACE-2))|(1<<(-X0_FACE-2+6)))) {
      cutCubeNodeVec[1].push_back(pair<double,int>(x_no[ino][2],ino));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-Y1_FACE-2+6))|(1<<(-X0_FACE-2)))) {
      cutCubeNodeVec[1].push_back(pair<double,int>(x_no[ino][2],-ino-1));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-Z0_FACE-2))|(1<<(-X0_FACE-2+6)))) {
      cutCubeNodeVec[2].push_back(pair<double,int>(x_no[ino][1],ino));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-Z0_FACE-2+6))|(1<<(-X0_FACE-2)))) {
      cutCubeNodeVec[2].push_back(pair<double,int>(x_no[ino][1],-ino-1));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-X0_FACE-2))|(1<<(-Z1_FACE-2+6)))) {
      cutCubeNodeVec[3].push_back(pair<double,int>(x_no[ino][1],ino));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-X0_FACE-2+6))|(1<<(-Z1_FACE-2)))) {
      cutCubeNodeVec[3].push_back(pair<double,int>(x_no[ino][1],-ino-1));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-Y0_FACE-2))|(1<<(-X1_FACE-2+6)))) {
      cutCubeNodeVec[4].push_back(pair<double,int>(x_no[ino][2],ino));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-Y0_FACE-2+6))|(1<<(-X1_FACE-2)))) {
      cutCubeNodeVec[4].push_back(pair<double,int>(x_no[ino][2],-ino-1));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-X1_FACE-2))|(1<<(-Y1_FACE-2+6)))) {
      cutCubeNodeVec[5].push_back(pair<double,int>(x_no[ino][2],ino));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-X1_FACE-2+6))|(1<<(-Y1_FACE-2)))) {
      cutCubeNodeVec[5].push_back(pair<double,int>(x_no[ino][2],-ino-1));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-X1_FACE-2))|(1<<(-Z0_FACE-2+6)))) {
      cutCubeNodeVec[6].push_back(pair<double,int>(x_no[ino][1],ino));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-X1_FACE-2+6))|(1<<(-Z0_FACE-2)))) {
      cutCubeNodeVec[6].push_back(pair<double,int>(x_no[ino][1],-ino-1));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-Z1_FACE-2))|(1<<(-X1_FACE-2+6)))) {
      cutCubeNodeVec[7].push_back(pair<double,int>(x_no[ino][1],ino));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-Z1_FACE-2+6))|(1<<(-X1_FACE-2)))) {
      cutCubeNodeVec[7].push_back(pair<double,int>(x_no[ino][1],-ino-1));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-Y0_FACE-2))|(1<<(-Z0_FACE-2+6)))) {
      cutCubeNodeVec[8].push_back(pair<double,int>(x_no[ino][0],ino));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-Y0_FACE-2+6))|(1<<(-Z0_FACE-2)))) {
      cutCubeNodeVec[8].push_back(pair<double,int>(x_no[ino][0],-ino-1));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-Z1_FACE-2))|(1<<(-Y0_FACE-2+6)))) {
      cutCubeNodeVec[9].push_back(pair<double,int>(x_no[ino][0],ino));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-Z1_FACE-2+6))|(1<<(-Y0_FACE-2)))) {
      cutCubeNodeVec[9].push_back(pair<double,int>(x_no[ino][0],-ino-1));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-Z0_FACE-2))|(1<<(-Y1_FACE-2+6)))) {
      cutCubeNodeVec[10].push_back(pair<double,int>(x_no[ino][0],ino));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-Z0_FACE-2+6))|(1<<(-Y1_FACE-2)))) {
      cutCubeNodeVec[10].push_back(pair<double,int>(x_no[ino][0],-ino-1));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-Y1_FACE-2))|(1<<(-Z1_FACE-2+6)))) {
      cutCubeNodeVec[11].push_back(pair<double,int>(x_no[ino][0],ino));
      got_ccn = true;
    }
    else if (no_flag[ino] == ((1<<(-Y1_FACE-2+6))|(1<<(-Z1_FACE-2)))) {
      cutCubeNodeVec[11].push_back(pair<double,int>(x_no[ino][0],-ino-1));
      got_ccn = true;
    }
  }
  
  if (debug) {
    cout << "Summary of nodes along the cube edges: got_ccn: " << got_ccn << endl;
    FILE * fp = fopen("edge_intersection.dat","w");
    for (int i = 0; i < 12; ++i) {
      cout << "cutCubeNodeVec[" << i << "].size(): " << cutCubeNodeVec[i].size() << endl;
      for (int j = 0; j < cutCubeNodeVec[i].size(); ++j) {
        const double L = cutCubeNodeVec[i][j].first;
        int ino = cutCubeNodeVec[i][j].second;
        cout << " > L,ino: " << L << " " << ino << endl;
        if (ino < 0)
          ino = -ino-1;
        fprintf(fp,"%18.15e %18.15e %18.15e\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
      }
    }
    fclose(fp);
    cout << "checkout edge_intersection.dat" << endl; 
  }

  if (!got_ccn) {
    
    // there were no intersections of the cut surface with any of the 12 edges...
      
    // if we didn't find any cube intersections, then it may still be 
    // necessary to include the cube...

    double face_area[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    bool got_edge = false;
      
    FOR_IED {
      if (faoed[ied][0] < -1) {
	// even one edge will be enough for us to determine the surface orientation...
	got_edge = true;
	const int ic = -faoed[ied][0]-2;
	assert((ic >= 0)&&(ic < 6)); // -x,+x,-y,+y,-z,+z
	const int ino0 = nooed[ied][0]; assert(ino0 >= 0);
	const int ino1 = nooed[ied][1]; assert(ino1 >= 0);
	switch (ic) {
	case 0:
	  // surface on -x...
	  assert(x_no[ino0][0] == -Lx);
	  assert(x_no[ino1][0] == -Lx);
	  face_area[0] += x_no[ino0][1]*x_no[ino1][2] - x_no[ino0][2]*x_no[ino1][1];
	  break;
	case 1:
	  // surface on +x...
	  assert(x_no[ino0][0] == Lx);
	  assert(x_no[ino1][0] == Lx);
	  face_area[1] -= x_no[ino0][1]*x_no[ino1][2] - x_no[ino0][2]*x_no[ino1][1];
	  break;
	case 2:
	  // surface on -y...
	  assert(x_no[ino0][1] == -Ly);
	  assert(x_no[ino1][1] == -Ly);
	  face_area[2] += x_no[ino0][2]*x_no[ino1][0] - x_no[ino0][0]*x_no[ino1][2];
	  break;
	case 3:
	  // surface on +y...
	  assert(x_no[ino0][1] == Ly);
	  assert(x_no[ino1][1] == Ly);
	  face_area[3] -= x_no[ino0][2]*x_no[ino1][0] - x_no[ino0][0]*x_no[ino1][2];
	  break;
	case 4:
	  // surface on -z...
	  assert(x_no[ino0][2] == -Lz);
	  assert(x_no[ino1][2] == -Lz);
	  face_area[4] += x_no[ino0][0]*x_no[ino1][1] - x_no[ino0][1]*x_no[ino1][0];
	  break;
	case 5:
	  // surface on +z...
	  assert(x_no[ino0][2] == Lz);
	  assert(x_no[ino1][2] == Lz);
	  face_area[5] -= x_no[ino0][0]*x_no[ino1][1] - x_no[ino0][1]*x_no[ino1][0];
	  break;
	default:
	  assert(0);
	}
      }
    }
      
    if (!got_edge) {

      double V_sum = 0.0;
      map<const int,double*> xMap;
      FOR_IED {
	map<const int,double*>::iterator iter = xMap.find(faoed[ied][0]);
	if (iter == xMap.end()) {
	  xMap.insert(pair<int,double*>(faoed[ied][0],x_no[nooed[ied][0]]));
	}
	else if (x_no[nooed[ied][1]] != iter->second) {
	  V_sum += CROSS_DOT(iter->second,x_no[nooed[ied][0]],x_no[nooed[ied][1]]);
	}
	iter = xMap.find(faoed[ied][1]);
	if (iter == xMap.end()) {
	  xMap.insert(pair<int,double*>(faoed[ied][1],x_no[nooed[ied][1]]));
	}
	else if (x_no[nooed[ied][0]] != iter->second) {
	  V_sum += CROSS_DOT(iter->second,x_no[nooed[ied][1]],x_no[nooed[ied][0]]);
	}
      }

      if (V_sum < 0.0) {
        
	// here we have grabbed an isolated void. Add the cube around it and continue...

	int ino;
	ino = new_node(); x_no[ino][0] = -Lx; x_no[ino][1] = -Ly; x_no[ino][2] = -Lz;
	ino = new_node(); x_no[ino][0] =  Lx; x_no[ino][1] = -Ly; x_no[ino][2] = -Lz;
	ino = new_node(); x_no[ino][0] = -Lx; x_no[ino][1] =  Ly; x_no[ino][2] = -Lz; 
	ino = new_node(); x_no[ino][0] =  Lx; x_no[ino][1] =  Ly; x_no[ino][2] = -Lz; 
	ino = new_node(); x_no[ino][0] = -Lx; x_no[ino][1] = -Ly; x_no[ino][2] =  Lz;
	ino = new_node(); x_no[ino][0] =  Lx; x_no[ino][1] = -Ly; x_no[ino][2] =  Lz;
	ino = new_node(); x_no[ino][0] = -Lx; x_no[ino][1] =  Ly; x_no[ino][2] =  Lz; 
	ino = new_node(); x_no[ino][0] =  Lx; x_no[ino][1] =  Ly; x_no[ino][2] =  Lz;
	assert(nno == nno_orig+8);

	int ied;
	// x-edges...
	ied = new_edge(); nooed[ied][0] = nno_orig+0; nooed[ied][1] = nno_orig+1; faoed[ied][0] = Y0_FACE; faoed[ied][1] = Z0_FACE;
	ied = new_edge(); nooed[ied][0] = nno_orig+2; nooed[ied][1] = nno_orig+3; faoed[ied][0] = Z0_FACE; faoed[ied][1] = Y1_FACE;
	ied = new_edge(); nooed[ied][0] = nno_orig+4; nooed[ied][1] = nno_orig+5; faoed[ied][0] = Z1_FACE; faoed[ied][1] = Y0_FACE;
	ied = new_edge(); nooed[ied][0] = nno_orig+6; nooed[ied][1] = nno_orig+7; faoed[ied][0] = Y1_FACE; faoed[ied][1] = Z1_FACE;
	// y-edges...
	ied = new_edge(); nooed[ied][0] = nno_orig+0; nooed[ied][1] = nno_orig+2; faoed[ied][0] = Z0_FACE; faoed[ied][1] = X0_FACE;
	ied = new_edge(); nooed[ied][0] = nno_orig+1; nooed[ied][1] = nno_orig+3; faoed[ied][0] = X1_FACE; faoed[ied][1] = Z0_FACE;
	ied = new_edge(); nooed[ied][0] = nno_orig+4; nooed[ied][1] = nno_orig+6; faoed[ied][0] = X0_FACE; faoed[ied][1] = Z1_FACE;
	ied = new_edge(); nooed[ied][0] = nno_orig+5; nooed[ied][1] = nno_orig+7; faoed[ied][0] = Z1_FACE; faoed[ied][1] = X1_FACE;
	// z-edges...
	ied = new_edge(); nooed[ied][0] = nno_orig+0; nooed[ied][1] = nno_orig+4; faoed[ied][0] = X0_FACE; faoed[ied][1] = Y0_FACE;
	ied = new_edge(); nooed[ied][0] = nno_orig+1; nooed[ied][1] = nno_orig+5; faoed[ied][0] = Y0_FACE; faoed[ied][1] = X1_FACE;
	ied = new_edge(); nooed[ied][0] = nno_orig+2; nooed[ied][1] = nno_orig+6; faoed[ied][0] = Y1_FACE; faoed[ied][1] = X0_FACE;
	ied = new_edge(); nooed[ied][0] = nno_orig+3; nooed[ied][1] = nno_orig+7; faoed[ied][0] = X1_FACE; faoed[ied][1] = Y1_FACE;
        
      }
        
      // no need to add the outer cube: a positive V_sum is simply a very coarse grid...
        
    }
    else if ( (face_area[0] >= 0.0) && (face_area[1] >= 0.0) && (face_area[2] >= 0.0) &&
	      (face_area[3] >= 0.0) && (face_area[4] >= 0.0) && (face_area[5] >= 0.0) ) {

      // the surface is penetrating into the cube. We need to add the full outer cube...
      //cout << "add all cube nodes and edges" << endl;

      assert(nno == nno_orig);
      int ino;
      ino = new_node(); x_no[ino][0] = -Lx; x_no[ino][1] = -Ly; x_no[ino][2] = -Lz;
      ino = new_node(); x_no[ino][0] =  Lx; x_no[ino][1] = -Ly; x_no[ino][2] = -Lz;
      ino = new_node(); x_no[ino][0] = -Lx; x_no[ino][1] =  Ly; x_no[ino][2] = -Lz; 
      ino = new_node(); x_no[ino][0] =  Lx; x_no[ino][1] =  Ly; x_no[ino][2] = -Lz; 
      ino = new_node(); x_no[ino][0] = -Lx; x_no[ino][1] = -Ly; x_no[ino][2] =  Lz;
      ino = new_node(); x_no[ino][0] =  Lx; x_no[ino][1] = -Ly; x_no[ino][2] =  Lz;
      ino = new_node(); x_no[ino][0] = -Lx; x_no[ino][1] =  Ly; x_no[ino][2] =  Lz; 
      ino = new_node(); x_no[ino][0] =  Lx; x_no[ino][1] =  Ly; x_no[ino][2] =  Lz;
        
      int ied;
      // x-edges...
      ied = new_edge(); nooed[ied][0] = nno_orig+0; nooed[ied][1] = nno_orig+1; faoed[ied][0] = Y0_FACE; faoed[ied][1] = Z0_FACE;
      ied = new_edge(); nooed[ied][0] = nno_orig+2; nooed[ied][1] = nno_orig+3; faoed[ied][0] = Z0_FACE; faoed[ied][1] = Y1_FACE;
      ied = new_edge(); nooed[ied][0] = nno_orig+4; nooed[ied][1] = nno_orig+5; faoed[ied][0] = Z1_FACE; faoed[ied][1] = Y0_FACE;
      ied = new_edge(); nooed[ied][0] = nno_orig+6; nooed[ied][1] = nno_orig+7; faoed[ied][0] = Y1_FACE; faoed[ied][1] = Z1_FACE;
      // y-edges...
      ied = new_edge(); nooed[ied][0] = nno_orig+0; nooed[ied][1] = nno_orig+2; faoed[ied][0] = Z0_FACE; faoed[ied][1] = X0_FACE;
      ied = new_edge(); nooed[ied][0] = nno_orig+1; nooed[ied][1] = nno_orig+3; faoed[ied][0] = X1_FACE; faoed[ied][1] = Z0_FACE;
      ied = new_edge(); nooed[ied][0] = nno_orig+4; nooed[ied][1] = nno_orig+6; faoed[ied][0] = X0_FACE; faoed[ied][1] = Z1_FACE;
      ied = new_edge(); nooed[ied][0] = nno_orig+5; nooed[ied][1] = nno_orig+7; faoed[ied][0] = Z1_FACE; faoed[ied][1] = X1_FACE;
      // z-edges...
      ied = new_edge(); nooed[ied][0] = nno_orig+0; nooed[ied][1] = nno_orig+4; faoed[ied][0] = X0_FACE; faoed[ied][1] = Y0_FACE;
      ied = new_edge(); nooed[ied][0] = nno_orig+1; nooed[ied][1] = nno_orig+5; faoed[ied][0] = Y0_FACE; faoed[ied][1] = X1_FACE;
      ied = new_edge(); nooed[ied][0] = nno_orig+2; nooed[ied][1] = nno_orig+6; faoed[ied][0] = Y1_FACE; faoed[ied][1] = X0_FACE;
      ied = new_edge(); nooed[ied][0] = nno_orig+3; nooed[ied][1] = nno_orig+7; faoed[ied][0] = X1_FACE; faoed[ied][1] = Y1_FACE;
      
    }
    else {
        
      // no need to add the outer cube, but 
      // make sure all areas have the expected sign...

      if (!( (face_area[0] <= 0.0) && (face_area[1] <= 0.0) && (face_area[2] <= 0.0) &&
	     (face_area[3] <= 0.0) && (face_area[4] <= 0.0) && (face_area[5] <= 0.0) ) ) {
	cout << "Error: CuttableVoronoiData: failed face_area check: " << 
	  face_area[0] << " " << 
	  face_area[1] << " " << 
	  face_area[2] << " " << 
	  face_area[3] << " " << 
	  face_area[4] << " " << 
	  face_area[5] << endl;
	throw(-1);
      }
        
    }

  }
  else {
    
    // set a corner_flag from 1 to 0 (in) or -1 (out) based on the orientation of
    // the edge intersections...
    
    int corner_flag[8];
    FOR_I8 corner_flag[i] = 1;
    
    // edge 0 faces x0,y0 nodes 0,4
#define IED 0
#define INO0 0
#define INO1 4
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
    // edge 1 faces y1,x0 nodes 2,6
#define IED 1
#define INO0 2
#define INO1 6
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
    // edge 2 faces z0,x0 nodes 0,2
#define IED 2
#define INO0 0
#define INO1 2
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
    // edge 3 faces x0,z1 nodes 4,6
#define IED 3
#define INO0 4
#define INO1 6
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
    // edge 4 faces y0,x1 nodes 1,5
#define IED 4
#define INO0 1
#define INO1 5
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
    // edge 5 faces x1,y1 nodes 3,7
#define IED 5
#define INO0 3
#define INO1 7
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
    // edge 6 faces x1,z0 nodes 1,3
#define IED 6
#define INO0 1
#define INO1 3
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
    // edge 7 faces z1,x1 nodes 5,7
#define IED 7
#define INO0 5
#define INO1 7
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
    // edge 8 faces y0,z0 nodes 0,1
#define IED 8
#define INO0 0
#define INO1 1
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
    // edge 9 faces z1,y0 nodes 4,5
#define IED 9
#define INO0 4
#define INO1 5
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
    // edge 10 faces z0,y1 nodes 2,3
#define IED 10
#define INO0 2
#define INO1 3
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
    // edge 11 faces y1,z1 nodes 6,7
#define IED 11
#define INO0 6
#define INO1 7
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1

    // and add the corner nodes (we may need some of these)...
    assert(nno == nno_orig);
    resize_nno(nno_orig+8);
    assert(nno == nno_orig+8);
    x_no[nno_orig  ][0] = -Lx; x_no[nno_orig  ][1] = -Ly; x_no[nno_orig  ][2] = -Lz;
    x_no[nno_orig+1][0] =  Lx; x_no[nno_orig+1][1] = -Ly; x_no[nno_orig+1][2] = -Lz;
    x_no[nno_orig+2][0] = -Lx; x_no[nno_orig+2][1] =  Ly; x_no[nno_orig+2][2] = -Lz; 
    x_no[nno_orig+3][0] =  Lx; x_no[nno_orig+3][1] =  Ly; x_no[nno_orig+3][2] = -Lz; 
    x_no[nno_orig+4][0] = -Lx; x_no[nno_orig+4][1] = -Ly; x_no[nno_orig+4][2] =  Lz;
    x_no[nno_orig+5][0] =  Lx; x_no[nno_orig+5][1] = -Ly; x_no[nno_orig+5][2] =  Lz;
    x_no[nno_orig+6][0] = -Lx; x_no[nno_orig+6][1] =  Ly; x_no[nno_orig+6][2] =  Lz; 
    x_no[nno_orig+7][0] =  Lx; x_no[nno_orig+7][1] =  Ly; x_no[nno_orig+7][2] =  Lz;
      
    // and then loop edges to set all nodes...
    int done = 0;
    while (done == 0) {
      done = 1; // assume done
      // edge 0 faces x0,y0 nodes 0,4
#define IED 0
#define INO0 0
#define INO1 4
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 1 faces y1,x0 nodes 2,6
#define IED 1
#define INO0 2
#define INO1 6
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 2 faces z0,x0 nodes 0,2
#define IED 2
#define INO0 0
#define INO1 2
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 3 faces x0,z1 nodes 4,6
#define IED 3
#define INO0 4
#define INO1 6
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 4 faces y0,x1 nodes 1,5
#define IED 4
#define INO0 1
#define INO1 5
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 5 faces x1,y1 nodes 3,7
#define IED 5
#define INO0 3
#define INO1 7
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 6 faces x1,z0 nodes 1,3
#define IED 6
#define INO0 1
#define INO1 3
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 7 faces z1,x1 nodes 5,7
#define IED 7
#define INO0 5
#define INO1 7
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 8 faces y0,z0 nodes 0,1
#define IED 8
#define INO0 0
#define INO1 1
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 9 faces z1,y0 nodes 4,5
#define IED 9
#define INO0 4
#define INO1 5
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 10 faces z0,y1 nodes 2,3
#define IED 10
#define INO0 2
#define INO1 3
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 11 faces y1,z1 nodes 6,7
#define IED 11
#define INO0 6
#define INO1 7
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
    }
    
    // at this point, the corners that will remain are flagged with 0. Use this to add 
    // all required edges along the sides of the cube...
    
    // also use no_flag here to indicate who is out/in...
    for (int ino = 0; ino < nno; ++ino) no_flag[ino] = 0; // 0 = unknown, 1 = in, -1 = out
    
    // the code below in corner_flag3.hpp will set any "in" nodes back to 0. For now, we assume all nodes 
    // along the edges are out, including the 8 corners. Note that just because a node is used in 
    // producing an edge along the cut cube boundary does NOT make it automatically in, because
    // there are certain pathological cases of curved surfaces crossing an edge in what seems like 
    // a correct orientation, however they are already entirely on the wrong side of another surface
    // and will be cut away in the walking done later. So all we can say is that any edge node that does
    // NOT involve a cube edge is DEFINATELY out...
    FOR_I8 no_flag[nno_orig+i] = -1; // mark as out -- this gets turned back to 0 in corner_flag3.hpp below when connected to an edge
    for (int i = 0; i < 12; ++i) {
      for (int j = 0; j < cutCubeNodeVec[i].size(); ++j) {
	const int ino = max(cutCubeNodeVec[i][j].second,-cutCubeNodeVec[i][j].second-1);
	assert((ino >= 0)&&(ino < nno_orig));
	no_flag[ino] = -1; // mark as out -- this gets turned back to 0 in corner_flag3.hpp below when connected to an edge
      }
    }

    // edge 0 faces x0,y0 nodes 0,4
#define IED 0
#define INO0 0
#define INO1 4
#define IFA0 X0_FACE
#define IFA1 Y0_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
    // edge 1 faces y1,x0 nodes 2,6
#define IED 1
#define INO0 2
#define INO1 6
#define IFA0 Y1_FACE
#define IFA1 X0_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
    // edge 2 faces z0,x0 nodes 0,2
#define IED 2
#define INO0 0
#define INO1 2
#define IFA0 Z0_FACE
#define IFA1 X0_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
    // edge 3 faces x0,z1 nodes 4,6
#define IED 3
#define INO0 4
#define INO1 6
#define IFA0 X0_FACE
#define IFA1 Z1_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
    // edge 4 faces y0,x1 nodes 1,5
#define IED 4
#define INO0 1
#define INO1 5
#define IFA0 Y0_FACE
#define IFA1 X1_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
    // edge 5 faces x1,y1 nodes 3,7
#define IED 5
#define INO0 3
#define INO1 7
#define IFA0 X1_FACE
#define IFA1 Y1_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
    // edge 6 faces x1,z0 nodes 1,3
#define IED 6
#define INO0 1
#define INO1 3
#define IFA0 X1_FACE
#define IFA1 Z0_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
    // edge 7 faces z1,x1 nodes 5,7
#define IED 7
#define INO0 5
#define INO1 7
#define IFA0 Z1_FACE
#define IFA1 X1_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
    // edge 8 faces y0,z0 nodes 0,1
#define IED 8
#define INO0 0
#define INO1 1
#define IFA0 Y0_FACE
#define IFA1 Z0_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
    // edge 9 faces z1,y0 nodes 4,5
#define IED 9
#define INO0 4
#define INO1 5
#define IFA0 Z1_FACE
#define IFA1 Y0_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
    // edge 10 faces z0,y1 nodes 2,3
#define IED 10
#define INO0 2
#define INO1 3
#define IFA0 Z0_FACE
#define IFA1 Y1_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
    // edge 11 faces y1,z1 nodes 6,7
#define IED 11
#define INO0 6
#define INO1 7
#define IFA0 Y1_FACE
#define IFA1 Z1_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1

    // nodes are now flagged as out-for-sure (-1), unknown (0). None are considered in.
    // and we should have added edges...

    assert(ned > ned_orig);
    
    if (debug) {
      writeTecplot(3);
      cout << "tecplot 3: This is the surface after adding cube edges" << endl;
      FILE * fp = fopen("nodes0.dat","w");
      FOR_INO {
	fprintf(fp,"%18.16e %18.16e %18.16e %d\n",x_no[ino][0],x_no[ino][1],x_no[ino][2],no_flag[ino]);
      }
      fclose(fp);
      cout << "nodes0.dat contains nodes as flagged after corner_flag3.hpp code (expect 0,-1)" << endl;
    }
    
    // at this point, shift the meaning of no_flag to store bits...
    
    const int in_bit       = (1<<0);
    const int out_bit      = (1<<1);
    const int in_link_bit  = (1<<2);
    const int out_link_bit = (1<<3);
    //const int check_bit    = (1<<4);
    
    FOR_INO {
      if (no_flag[ino] == -1) {
        no_flag[ino] = out_bit;
      }
      else {
        assert(no_flag[ino] == 0);
      }
    }

    // now we need to flag edges as in (+1) or out (-1)...

    // intersection edges ("links") are pairs of edges that share the same nodes, and are 
    // ordered such that the first is "in" and the second is "out". No other edges should have 
    // this property at this point, so use this to identify them and flag as in/out...
    // TODO: can we retain knowledge of these intersection edges as they pass through the 
    // cutting routine above?: for example, can the cvd.cut_surf() routine carry an edge flag, 
    // or use the sign of nooed, for example...
    // One thing that should be true of these edges is that they are beside
    // each other: so loop and compare ied-1 and ied...
    
    int * ed_flag = new int[ned];
    FOR_IED ed_flag[ied] = 0;

    vector<int> inLinkVec;
    vector<int> outLinkVec;
    for (int ied = 1; ied < ned; ++ied) {
      // intersection edge section. Here we need to skip the first because we
      // work with an edge pair...
      // look for links: same nodes...
      if ((nooed[ied][0] == nooed[ied-1][0])&&(nooed[ied][1] == nooed[ied-1][1])) {
	// these 2 edges point to exactly the same node...
	assert(ed_flag[ied-1] == 0);
	assert(ed_flag[ied] == 0);
	ed_flag[ied-1] = 2; // first is in -- note use 2/-2 for these "link" edges...
	ed_flag[ied] = -2; // second is out
	inLinkVec.push_back(ied-1);
	outLinkVec.push_back(ied);
	// and mark these nodes as touching an inlink and an out link....
	// (same nodes in this case)...
	no_flag[nooed[ied-1][0]] |= in_link_bit;
	no_flag[nooed[ied-1][1]] |= in_link_bit;
	no_flag[nooed[ied][0]] |= out_link_bit;
	no_flag[nooed[ied][1]] |= out_link_bit;
      }
      else if ((x_no[nooed[ied][0]][0] == x_no[nooed[ied-1][0]][0])&&
	       (x_no[nooed[ied][0]][1] == x_no[nooed[ied-1][0]][1])&&
	       (x_no[nooed[ied][0]][2] == x_no[nooed[ied-1][0]][2])&&
	       (x_no[nooed[ied][1]][0] == x_no[nooed[ied-1][1]][0])&&
	       (x_no[nooed[ied][1]][1] == x_no[nooed[ied-1][1]][1])&&
	       (x_no[nooed[ied][1]][2] == x_no[nooed[ied-1][1]][2])) {
	// on the ends of the set of linked edges where it comes to one of the 
	// cube faces, there will be two nodes (because of the edge cutting). These
	// are still link edges, however. The nodes should be EXACTLY identical. 
	// (TODO: see note above about remembering which edges are link edges)
	assert(ed_flag[ied-1] == 0);
	assert(ed_flag[ied] == 0);
	ed_flag[ied-1] = 2; // first is in -- note use 2/-2 for these "link" edges...
	ed_flag[ied] = -2; // second is out
	inLinkVec.push_back(ied-1);
	outLinkVec.push_back(ied);
	// nodes may be different...
	no_flag[nooed[ied-1][0]] |= in_link_bit;
	no_flag[nooed[ied-1][1]] |= in_link_bit;
	no_flag[nooed[ied][0]] |= out_link_bit;
	no_flag[nooed[ied][1]] |= out_link_bit;
      }
      // as a check, make sure nothing matches in reverse...
      if ((x_no[nooed[ied][0]][0] == x_no[nooed[ied-1][1]][0])&&
	  (x_no[nooed[ied][0]][1] == x_no[nooed[ied-1][1]][1])&&
	  (x_no[nooed[ied][0]][2] == x_no[nooed[ied-1][1]][2])&&
	  (x_no[nooed[ied][1]][0] == x_no[nooed[ied-1][0]][0])&&
	  (x_no[nooed[ied][1]][1] == x_no[nooed[ied-1][0]][1])&&
	  (x_no[nooed[ied][1]][2] == x_no[nooed[ied-1][0]][2])) {
	cout << "Error: reverse edge match" << endl;
	throw(-1);
      }
    }
    
    if (debug) {
      FILE * fp = fopen("links.dat","w");
      FOR_IED {
        if ((ed_flag[ied] == 2)||(ed_flag[ied] == -2)) {
          writeEdge(fp,x_no[nooed[ied][0]],x_no[nooed[ied][1]]);
        }
      }
      cout << "links.dat contains intersection edges" << endl;
    }
    
    // now use the links to flag the nodes immediately around the links as either
    // in or out using the link orientation. Note that this also double-checks
    // that the link orientation in any group (in or out) is self-consistent...
    
    set<pair<int,int> > forwardFaNoInSet;
    set<pair<int,int> > backwardFaNoInSet;
    for (int ii = 0; ii < inLinkVec.size(); ++ii) {
      const int ied = inLinkVec[ii];
      assert(ed_flag[ied] = 2); // all of these are in
      bool f0 = (forwardFaNoInSet.find(pair<int,int>(faoed[ied][0],nooed[ied][0])) == forwardFaNoInSet.end());
      bool f1 = (forwardFaNoInSet.find(pair<int,int>(faoed[ied][1],nooed[ied][1])) == forwardFaNoInSet.end());
      bool b0 = (backwardFaNoInSet.find(pair<int,int>(faoed[ied][0],nooed[ied][1])) == backwardFaNoInSet.end());
      bool b1 = (backwardFaNoInSet.find(pair<int,int>(faoed[ied][1],nooed[ied][0])) == backwardFaNoInSet.end());
      if (f0 && f1 && b0 && b1) {
	// looks good...
	forwardFaNoInSet.insert(pair<int,int>(faoed[ied][0],nooed[ied][0]));
	forwardFaNoInSet.insert(pair<int,int>(faoed[ied][1],nooed[ied][1]));
	backwardFaNoInSet.insert(pair<int,int>(faoed[ied][0],nooed[ied][1]));
	backwardFaNoInSet.insert(pair<int,int>(faoed[ied][1],nooed[ied][0]));
      }
      else {
        cout << "In link orientation problem: " << f0 << " " << f1 << " " << b0 << " " << b1 << endl;
        throw(-1);
      }
    }
    
    set<pair<int,int> > forwardFaNoOutSet;
    set<pair<int,int> > backwardFaNoOutSet;
    for (int ii = 0; ii < outLinkVec.size(); ++ii) {
      const int ied = outLinkVec[ii];
      assert(ed_flag[ied] = -2); // all of these are in
      bool f0 = (forwardFaNoOutSet.find(pair<int,int>(faoed[ied][0],nooed[ied][0])) == forwardFaNoOutSet.end());
      bool f1 = (forwardFaNoOutSet.find(pair<int,int>(faoed[ied][1],nooed[ied][1])) == forwardFaNoOutSet.end());
      bool b0 = (backwardFaNoOutSet.find(pair<int,int>(faoed[ied][0],nooed[ied][1])) == backwardFaNoOutSet.end());
      bool b1 = (backwardFaNoOutSet.find(pair<int,int>(faoed[ied][1],nooed[ied][0])) == backwardFaNoOutSet.end());
      if (f0 && f1 && b0 && b1) {
	// looks good...
	forwardFaNoOutSet.insert(pair<int,int>(faoed[ied][0],nooed[ied][0]));
	forwardFaNoOutSet.insert(pair<int,int>(faoed[ied][1],nooed[ied][1]));
	backwardFaNoOutSet.insert(pair<int,int>(faoed[ied][0],nooed[ied][1]));
	backwardFaNoOutSet.insert(pair<int,int>(faoed[ied][1],nooed[ied][0]));
      }
      else {
        cout << "Out link orientation problem: " << f0 << " " << f1 << " " << b0 << " " << b1 << endl;
        throw(-1);
      }
    }

    // final loop through edges to prepare the walk...
    
    FOR_IED {
      if ((faoed[ied][0] == -1)||(faoed[ied][1] == -1)) {
	// an edge along an uncut boundary is out (e.g. the exposed edge of a piston
	// surface passing through the cylinder)...
	ed_flag[ied] = -1;
	no_flag[nooed[ied][0]] |= out_bit;	
	no_flag[nooed[ied][1]] |= out_bit;
      }
      else if (ed_flag[ied] == 0) {
	bool b_link = false;
	if (no_flag[nooed[ied][0]]&in_link_bit) {
	  b_link = true;
	  if ( (forwardFaNoInSet.find(pair<int,int>(faoed[ied][1],nooed[ied][0])) != forwardFaNoInSet.end()) ||
	       (backwardFaNoInSet.find(pair<int,int>(faoed[ied][0],nooed[ied][0])) != backwardFaNoInSet.end()) ) {
	    ed_flag[ied] = 1;
	    no_flag[nooed[ied][1]] |= in_bit;
	  }
	}
	if (no_flag[nooed[ied][0]]&out_link_bit) {
	  b_link = true;
	  if ( (forwardFaNoOutSet.find(pair<int,int>(faoed[ied][1],nooed[ied][0])) != forwardFaNoOutSet.end()) ||
	       (backwardFaNoOutSet.find(pair<int,int>(faoed[ied][0],nooed[ied][0])) != backwardFaNoOutSet.end()) ) {
	    assert(ed_flag[ied] != 1);
	    ed_flag[ied] = -1;
	    no_flag[nooed[ied][1]] |= out_bit;
	  }
	}
	if (no_flag[nooed[ied][1]]&in_link_bit) {
          b_link = true;
	  if ( (forwardFaNoInSet.find(pair<int,int>(faoed[ied][0],nooed[ied][1])) != forwardFaNoInSet.end()) ||
	       (backwardFaNoInSet.find(pair<int,int>(faoed[ied][1],nooed[ied][1])) != backwardFaNoInSet.end()) ) {
            assert(ed_flag[ied] != -1);
            ed_flag[ied] = 1;
            no_flag[nooed[ied][0]] |= in_bit;
	  }
	}
	if (no_flag[nooed[ied][1]]&out_link_bit) {
          b_link = true;
	  if ( (forwardFaNoOutSet.find(pair<int,int>(faoed[ied][0],nooed[ied][1])) != forwardFaNoOutSet.end()) ||
	       (backwardFaNoOutSet.find(pair<int,int>(faoed[ied][1],nooed[ied][1])) != backwardFaNoOutSet.end()) ) {
            assert(ed_flag[ied] != 1);
	    ed_flag[ied] = -1;
            no_flag[nooed[ied][0]] |= out_bit;
	  }
	}
        // if we touch a link, make sure we got set...
        if (b_link) {
          if (ed_flag[ied] == 0) {
            cout << "This edge was not set!" << endl;
            throw(-1);
          }
        }
      }
    }
    
    if (debug) {
      FILE * fp = fopen("nodes.dat","w");
      FOR_INO {
        if ((no_flag[ino]&in_link_bit)||(no_flag[ino]&out_link_bit)) {
          fprintf(fp,"%18.16e %18.16e %18.16e 0\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
        }
        else if (no_flag[ino]&in_bit) {
          assert(!(no_flag[ino]&out_bit));
          fprintf(fp,"%18.16e %18.16e %18.16e 1\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
        }
        else if (no_flag[ino]&out_bit) {
          assert(!(no_flag[ino]&in_bit));
          fprintf(fp,"%18.16e %18.16e %18.16e -1\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
        }
        else {
          assert(no_flag[ino] == 0);
        }
      }
      fclose(fp);
      cout << "nodes.dat contains nodes initial node grouping: intersection (0), in (+1), and out (-1)" << endl;
    }

    // now loop until there are no changes...
    done = 0;
    while (done == 0) {
      done = 1;
      FOR_IED {
        if (ed_flag[ied] == 0) {
          // none of these remaining "0" edges should touch any link nodes...
          assert(!(no_flag[nooed[ied][0]]&in_link_bit));
          assert(!(no_flag[nooed[ied][0]]&out_link_bit));
          const bool b_in = ((no_flag[nooed[ied][0]]&in_bit)||(no_flag[nooed[ied][1]]&in_bit));
          const bool b_out = ((no_flag[nooed[ied][0]]&out_bit)||(no_flag[nooed[ied][1]]&out_bit));
          if (b_in) {
            // cannot be both...
            if (b_out) {
              cout << "in/out walking violation." << endl;
              throw(-1);
            }
            ed_flag[ied] = 1;
            no_flag[nooed[ied][0]] |= in_bit;
            no_flag[nooed[ied][1]] |= in_bit;
            done = 0;
          }
          else if (b_out) {
            ed_flag[ied] = -1;
            no_flag[nooed[ied][0]] |= out_bit;
            no_flag[nooed[ied][1]] |= out_bit;
            done = 0;
          }
        }
      }
    }
    
    /*
      {
      FILE * fp = fopen("edges.dat","w");
      FOR_IED {
      if (ed_flag[ied] == 0) {
      writeEdge(fp,x_no[nooed[ied][0]],x_no[nooed[ied][1]]);
      }
      }
      fclose(fp);
      cout << "checkout zero edges: done: " << done << endl;
      getchar();
      }
    */
    
    // at this point, we may still have unflagged edges? For these edges,
    // assume they are in, althrough we really should do some kind of positive 
    // volume check for the truly pathological cases of intersecting complex surfaces...
    
    FOR_IED {
      if (ed_flag[ied] == 0) { 
        FOR_I2 {
          const int ino = nooed[ied][i];
          no_flag[ino] |= in_bit;
        }
      }
    }

    int nno_new = 0;
    FOR_INO {
      if (no_flag[ino]&(in_bit|in_link_bit)) {
        FOR_I3 x_no[nno_new][i] = x_no[ino][i];
        no_flag[ino] = nno_new++;
      }
      else {
        no_flag[ino] = -1;
      }
    }
    
    int ned_new = 0;
    FOR_IED {
      if (ed_flag[ied] >= 0) { // i.e. 0, 1, or 2...
        FOR_I2 faoed[ned_new][i] = faoed[ied][i];
        FOR_I2 nooed[ned_new][i] = no_flag[nooed[ied][i]];
        FOR_I2 assert((nooed[ned_new][i] >= 0)&&(nooed[ned_new][i] < nno_new));
        ++ned_new;
      }
    }
    
    // and take the new counts...
    nno = nno_new;
    ned = ned_new;
    
    delete[] ed_flag;
    
  }
  
  delete[] no_flag; 
  
  return 0;
  
}

#define FOR_ITR for (int itr = 0; itr < ntr; ++itr)

void CuttableVoronoiData::buildSeed(const double this_delta_cv,const int * const ist_and_ipart_buf,const bool debug) {
 
  // if the cvd is empty, just add a cube of 0.505*this_delta_cv...

  if (empty()) {
   
    addCube(0.505*this_delta_cv);
    return;
  
  }

  if (debug) {
    writeTecplot(0);
    cout << "tecplot 0: This is the surface as passed to CuttableVoronoiData::buildSeed" << endl;
  }
  
  // count unique surfaces...

  set<int> partSet;
  for (int ied = 0; ied < ned; ++ied) {
    if (faoed[ied][1] >= 0) {
      const int ipart = (ist_and_ipart_buf[faoed[ied][1]+1]>>6); // recall first entry is ist_global, second is ipart/bits
      partSet.insert(ipart);
    }
    if (faoed[ied][0] >= 0) {
      const int ipart = (ist_and_ipart_buf[faoed[ied][0]+1]>>6); // recall first entry is ist_global, second is ipart/bits
      partSet.insert(ipart);
    }
  }
  
  if (debug) cout << "got partSet.size(): " << partSet.size() << endl;

  if (partSet.size() >= 2) {
    
    for (set<int>::iterator it = partSet.begin(); it != partSet.end(); ++it) {
      
      // skip the first...
      if (it == partSet.begin())
        continue;

      // set ipart1 from the set...
      const int ipart1 = *it;

      // here we are going to intersect the parts successively:
      // 0+1 => (0+1)
      // (0+1)+2 => ((0+1)+2)
      // etc...
    
      if (debug) cout << " > intersecting ipart1: " << ipart1 << " with everything lower..." << endl;

      map<const int,Bbox> bbox0Map;
      map<const int,Bbox> bbox1Map;
      for (int ied = 0; ied < ned; ++ied) {
	FOR_I2 {
	  if (faoed[ied][i] >= 0) {
	    const int ipart = (ist_and_ipart_buf[faoed[ied][i]+1]>>6); // recall first entry is ist_global, second is ipart/bits
	    if (ipart < ipart1) {
	      // this belongs in part0...
	      map<const int,Bbox>::iterator it = bbox0Map.find(faoed[ied][i]);
	      if (it == bbox0Map.end()) bbox0Map[faoed[ied][i]] = Bbox(x_no[nooed[ied][0]],x_no[nooed[ied][1]]);
	      else {
		it->second.expand(x_no[nooed[ied][0]]);
		it->second.expand(x_no[nooed[ied][1]]);
	      }
	    }
	    else if (ipart == ipart1) {
	      // this belongs in part1...
	      map<const int,Bbox>::iterator it = bbox1Map.find(faoed[ied][i]);
	      if (it == bbox1Map.end()) bbox1Map[faoed[ied][i]] = Bbox(x_no[nooed[ied][0]],x_no[nooed[ied][1]]);
	      else {
		it->second.expand(x_no[nooed[ied][0]]);
		it->second.expand(x_no[nooed[ied][1]]);
	      }
	    }
	  }
	}
      }

      // TODO: only bbox's that overlaped bbox need be considered. 
      
      int nbb0 = bbox0Map.size();
      int * faobb0 = new int[nbb0];
      double (*bbmin0)[3] = new double[nbb0][3];
      double (*bbmax0)[3] = new double[nbb0][3];
      
      int ibb0 = 0;
      for (map<const int,Bbox>::iterator it = bbox0Map.begin(); it != bbox0Map.end(); ++it) {
	faobb0[ibb0] = it->first;
	it->second.getBbox(bbmin0[ibb0],bbmax0[ibb0]);
	++ibb0;
      }
      assert(ibb0 == nbb0);
      
      Adt<double> * adt0 = new Adt<double>(nbb0,bbmin0,bbmax0);
      
      delete[] bbmin0;
      delete[] bbmax0;

      vector<pair<int,int> > facePairVec;
      set<int> faceSet;
      vector<int> intVec;
      for (map<const int,Bbox>::iterator it = bbox1Map.begin(); it != bbox1Map.end(); ++it) {
	double bbmin1[3];
	double bbmax1[3];
	it->second.getBbox(bbmin1,bbmax1);
	assert(intVec.empty());
	adt0->buildListForBBox(intVec,bbmin1,bbmax1);
	for (int ii = 0; ii < intVec.size(); ++ii) {
	  facePairVec.push_back(pair<int,int>(faobb0[intVec[ii]],it->first));
	  faceSet.insert(faobb0[intVec[ii]]);
	  faceSet.insert(it->first);
	}
	intVec.clear();
      }
      
      delete adt0;
      delete[] faobb0;

      // the face set contains the faces involved in the potential intersections. Now build the 
      // tris and references into the cvd edge structure...
      
      int ntr = 0;
      map<const int,int> faceMap;
      for (set<int>::iterator it = faceSet.begin(); it != faceSet.end(); ++it) 
	faceMap[*it] = ntr++;
      faceSet.clear();
      
      // here we assume they are tris. If not, we cannot handle this intersection right now...
      
      int (*edotr)[3] = new int[ntr][3];
      FOR_ITR FOR_I3 edotr[itr][i] = -1;
      
      for (int ied = 0; ied < ned; ++ied) {
	if (faoed[ied][0] >= 0) {
	  map<const int,int>::iterator it = faceMap.find(faoed[ied][0]);
	  if (it != faceMap.end()) {
	    // this face is participating on the current intersection calc...
	    const int ifa = it->second;
	    int j;
	    for (j = 0; j < 3; ++j) {
	      if (edotr[ifa][j] == -1)
		break;
	    }
	    assert(j < 3); // you will hit this if ned per face is more than 3 -- e.g. polyhedral faces...
	    edotr[ifa][j] = ied;
	  }
	}
	if (faoed[ied][1] >= 0) {
	  map<const int,int>::iterator it = faceMap.find(faoed[ied][1]);
	  if (it != faceMap.end()) {
	    // this face is participating on the current intersection calc...
	    const int ifa = it->second;
	    int j;
	    for (j = 0; j < 3; ++j) {
	      if (edotr[ifa][j] == -1)
		break;
	    }
	    assert(j < 3); // you will hit this if ned per face is more than 3 -- e.g. polyhedral faces...
	    edotr[ifa][j] = -ied-2; // -2 indexing
	  }
	}
      }
      
      // then populate nodes, and ensure/reorder edges into rh-rule order...

      int (*nootr)[3] = new int[ntr][3];
      
      FOR_ITR {
	// edge 0 sets nodes 0 and 1...
	if (edotr[itr][0] >= 0) {
	  nootr[itr][0] = nooed[edotr[itr][0]][0];
	  nootr[itr][1] = nooed[edotr[itr][0]][1];
	}
	else {
	  assert(edotr[itr][0] < -1);
	  nootr[itr][0] = nooed[-edotr[itr][0]-2][1];
	  nootr[itr][1] = nooed[-edotr[itr][0]-2][0];
	}
	// edge 1 matches node 1 and sets node 2...
	if (edotr[itr][1] >= 0) {
	  if (nooed[edotr[itr][1]][0] == nootr[itr][1]) {
	    nootr[itr][2] = nooed[edotr[itr][1]][1];
	  }
	  else {
	    const int tmp = edotr[itr][1];
	    edotr[itr][1] = edotr[itr][2];
	    edotr[itr][2] = tmp;
	    if (edotr[itr][1] >= 0) {
	      assert(nooed[edotr[itr][1]][0] == nootr[itr][1]);
	      nootr[itr][2] = nooed[edotr[itr][1]][1];
	    }
	    else {
	      assert(edotr[itr][1] < -1);
	      assert(nooed[-edotr[itr][1]-2][1] == nootr[itr][1]);
	      nootr[itr][2] = nooed[-edotr[itr][1]-2][0];
	    }
	  }
	}
	else {
	  assert(edotr[itr][1] < -1);
	  if (nooed[-edotr[itr][1]-2][1] == nootr[itr][1]) {
	    nootr[itr][2] = nooed[-edotr[itr][1]-2][0];
	  }
	  else {
	    const int tmp = edotr[itr][1];
	    edotr[itr][1] = edotr[itr][2];
	    edotr[itr][2] = tmp;
	    if (edotr[itr][1] >= 0) {
	      assert(nooed[edotr[itr][1]][0] == nootr[itr][1]);
	      nootr[itr][2] = nooed[edotr[itr][1]][1];
	    }
	    else {
	      assert(edotr[itr][1] < -1);
	      assert(nooed[-edotr[itr][1]-2][1] == nootr[itr][1]);
	      nootr[itr][2] = nooed[-edotr[itr][1]-2][0];
	    }
	  }
	}
	// edge 2 matches nodes 2 and 0...
	if (edotr[itr][2] >= 0) {
	  assert(nootr[itr][2] == nooed[edotr[itr][2]][0]);
	  assert(nootr[itr][0] == nooed[edotr[itr][2]][1]);
	}
	else {
	  assert(edotr[itr][2] < -1);
	  assert(nootr[itr][2] == nooed[-edotr[itr][2]-2][1]);
	  assert(nootr[itr][0] == nooed[-edotr[itr][2]-2][0]);
	}
      }

      //const int ned_old = ned;
      vector<IntersectionData> intersectionVec;
      // data required by TriTri routine...
      int idata[6];
      double ddata[4];
      double * x0[3];
      double * x1[3];
      for (int ii = 0; ii < facePairVec.size(); ++ii) {
	const int ifa0 = facePairVec[ii].first;
	const int ifa1 = facePairVec[ii].second;
	assert(ifa0 != ifa1);
	map<const int,int>::iterator it = faceMap.find(ifa0);
	assert(it != faceMap.end());
	const int itr0 = it->second;
	it = faceMap.find(ifa1);
	assert(it != faceMap.end());
	const int itr1 = it->second;
	x0[0] = x_no[nootr[itr0][0]];
	x0[1] = x_no[nootr[itr0][1]];
	x0[2] = x_no[nootr[itr0][2]];
	x1[0] = x_no[nootr[itr1][0]];
	x1[1] = x_no[nootr[itr1][1]];
	x1[2] = x_no[nootr[itr1][2]];	
	if (const int n = calcTriTriIntersection(idata,ddata,x0,x1)) {
	  if (n == -2) {
            cout << "Error: calcTriTriIntersection failed. Use CVD_DEBUG on cvd_debug-*.bin and cvd_debug_support-*.bin files to generate code" << endl;
            // in debug mode, we run the permutations...
            if (debug) {
              FOR_I3 FOR_J3 {
                x0[0] = x_no[nootr[itr0][i]];
                x0[1] = x_no[nootr[itr0][(i+1)%3]];
                x0[2] = x_no[nootr[itr0][(i+2)%3]];
                x1[0] = x_no[nootr[itr1][j]];
                x1[1] = x_no[nootr[itr1][(j+1)%3]];
                x1[2] = x_no[nootr[itr1][(j+2)%3]];	
                calcTriTriIntersection(idata,ddata,x0,x1);
                calcTriTriIntersection(idata,ddata,x1,x0);
                x0[0] = x_no[nootr[itr0][(i+2)%3]];
                x0[1] = x_no[nootr[itr0][(i+1)%3]];
                x0[2] = x_no[nootr[itr0][i]];
                x1[0] = x_no[nootr[itr1][j]];
                x1[1] = x_no[nootr[itr1][(j+1)%3]];
                x1[2] = x_no[nootr[itr1][(j+2)%3]];	
                calcTriTriIntersection(idata,ddata,x0,x1);
                calcTriTriIntersection(idata,ddata,x1,x0);
                x0[0] = x_no[nootr[itr0][i]];
                x0[1] = x_no[nootr[itr0][(i+1)%3]];
                x0[2] = x_no[nootr[itr0][(i+2)%3]];
                x1[0] = x_no[nootr[itr1][(j+2)%3]];
                x1[1] = x_no[nootr[itr1][(j+1)%3]];
                x1[2] = x_no[nootr[itr1][j]];	
                calcTriTriIntersection(idata,ddata,x0,x1);
                calcTriTriIntersection(idata,ddata,x1,x0);
                x0[0] = x_no[nootr[itr0][(i+2)%3]];
                x0[1] = x_no[nootr[itr0][(i+1)%3]];
                x0[2] = x_no[nootr[itr0][i]];
                x1[0] = x_no[nootr[itr1][(j+2)%3]];
                x1[1] = x_no[nootr[itr1][(j+1)%3]];
                x1[2] = x_no[nootr[itr1][j]];	
                calcTriTriIntersection(idata,ddata,x0,x1);
                calcTriTriIntersection(idata,ddata,x1,x0);
              }
            }
	    throw(-1);
	  }
          else if (n == 1) {
	    switch (idata[0]) {
            case EDGE_EDGE_INT:
	      if (edotr[itr0][idata[1]] >= 0) {
                if (edotr[itr1][idata[2]] >= 0) {
                  intersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,edotr[itr0][idata[1]],edotr[itr1][idata[2]],ddata[0],ddata[1]));
                }
                else {
                  assert(edotr[itr1][idata[2]] < -1);
                  intersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,edotr[itr0][idata[1]],-edotr[itr1][idata[2]]-2,ddata[0],1.0-ddata[1]));
                }
              }
              else {
		assert(edotr[itr0][idata[1]] < -1);
                if (edotr[itr1][idata[2]] >= 0) {
                  intersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,-edotr[itr0][idata[1]]-2,edotr[itr1][idata[2]],1.0-ddata[0],ddata[1]));
                }
                else {
                  assert(edotr[itr1][idata[2]] < -1);
                  intersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,-edotr[itr0][idata[1]]-2,-edotr[itr1][idata[2]]-2,1.0-ddata[0],1.0-ddata[1]));
                }
              }
              break;
            default:
              cout << "unsupported n == 1 intersection: " << idata[0] << endl;
              throw(-1);
            }
          }
          else if (n == 2) {
	    // add 2 new edges that join this pair of intersections...
	    int ied = new_edge();
	    faoed[ied][0] = ifa0; // needs to be the actual face int, not the local tri
	    faoed[ied][1] = ifa1;
	    nooed[ied][0] = nno+intersectionVec.size();
	    nooed[ied][1] = nno+intersectionVec.size()+1;
	    // for the second edge, flip the faces and keep the nodes the same. This ensures that
	    // even if they are cut, we can identify these edges later???... 
	    ied = new_edge();
	    faoed[ied][0] = ifa1;
	    faoed[ied][1] = ifa0;
	    nooed[ied][0] = nno+intersectionVec.size();
	    nooed[ied][1] = nno+intersectionVec.size()+1;
	    // ============================
	    // first intersection...
	    // ============================
	    switch (idata[0]) {
	    case EDGE_TRI_INT:
	      // idata[1] contains the edge index on the first tri...
	      if (edotr[itr0][idata[1]] >= 0) {
		intersectionVec.push_back(IntersectionData(EDGE_FACE_INTERSECTION,edotr[itr0][idata[1]],ifa1,ddata[0]));
	      }
	      else {
		assert(edotr[itr0][idata[1]] < -1);
		intersectionVec.push_back(IntersectionData(EDGE_FACE_INTERSECTION,-edotr[itr0][idata[1]]-2,ifa1,1.0-ddata[0]));
	      }
	      break;
	    case TRI_EDGE_INT:
	      // idata[1] contains the edge index on the second tri...
	      if (edotr[itr1][idata[1]] >= 0) {
		intersectionVec.push_back(IntersectionData(EDGE_FACE_INTERSECTION,edotr[itr1][idata[1]],ifa0,ddata[0]));
	      }
	      else {
		assert(edotr[itr1][idata[1]] < -1);
		intersectionVec.push_back(IntersectionData(EDGE_FACE_INTERSECTION,-edotr[itr1][idata[1]]-2,ifa0,1.0-ddata[0]));
	      }
	      break;
            case EDGE_EDGE_INT:
	      if (edotr[itr0][idata[1]] >= 0) {
                if (edotr[itr1][idata[2]] >= 0) {
                  intersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,edotr[itr0][idata[1]],edotr[itr1][idata[2]],ddata[0],ddata[1]));
                }
                else {
                  assert(edotr[itr1][idata[2]] < -1);
                  intersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,edotr[itr0][idata[1]],-edotr[itr1][idata[2]]-2,ddata[0],1.0-ddata[1]));
                }
              }
              else {
		assert(edotr[itr0][idata[1]] < -1);
                if (edotr[itr1][idata[2]] >= 0) {
                  intersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,-edotr[itr0][idata[1]]-2,edotr[itr1][idata[2]],1.0-ddata[0],ddata[1]));
                }
                else {
                  assert(edotr[itr1][idata[2]] < -1);
                  intersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,-edotr[itr0][idata[1]]-2,-edotr[itr1][idata[2]]-2,1.0-ddata[0],1.0-ddata[1]));
                }
              }
              break;
            default:
              cout << "unsupported n == 2 intersection: " << idata[0] << endl;
              throw(-1);
	    }
	    // ============================
	    // second intersection...
	    // ============================
	    switch (idata[3]) {
	    case EDGE_TRI_INT:
	      // idata[4] contains the edge index on the first tri...
	      if (edotr[itr0][idata[4]] >= 0) {
		intersectionVec.push_back(IntersectionData(EDGE_FACE_INTERSECTION,edotr[itr0][idata[4]],ifa1,ddata[2]));
	      }
	      else {
		assert(edotr[itr0][idata[4]] < -1);
		intersectionVec.push_back(IntersectionData(EDGE_FACE_INTERSECTION,-edotr[itr0][idata[4]]-2,ifa1,1.0-ddata[2]));
	      }
	      break;
	    case TRI_EDGE_INT:
	      // idata[4] contains the edge index on the second tri...
	      if (edotr[itr1][idata[4]] >= 0) {
		intersectionVec.push_back(IntersectionData(EDGE_FACE_INTERSECTION,edotr[itr1][idata[4]],ifa0,ddata[2]));
	      }
	      else {
		assert(edotr[itr1][idata[4]] < -1);
		intersectionVec.push_back(IntersectionData(EDGE_FACE_INTERSECTION,-edotr[itr1][idata[4]]-2,ifa0,1.0-ddata[2]));
	      }
	      break;
            case EDGE_EDGE_INT:
	      if (edotr[itr0][idata[4]] >= 0) {
                if (edotr[itr1][idata[5]] >= 0) {
                  intersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,edotr[itr0][idata[4]],edotr[itr1][idata[5]],ddata[2],ddata[3]));
                }
                else {
                  assert(edotr[itr1][idata[5]] < -1);
                  intersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,edotr[itr0][idata[4]],-edotr[itr1][idata[5]]-2,ddata[2],1.0-ddata[3]));
                }
              }
              else {
		assert(edotr[itr0][idata[4]] < -1);
                if (edotr[itr1][idata[5]] >= 0) {
                  intersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,-edotr[itr0][idata[4]]-2,edotr[itr1][idata[5]],1.0-ddata[2],ddata[3]));
                }
                else {
                  assert(edotr[itr1][idata[5]] < -1);
                  intersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,-edotr[itr0][idata[4]]-2,-edotr[itr1][idata[5]]-2,1.0-ddata[2],1.0-ddata[3]));
                }
              }
              break;
	    default:
              cout << "unsupported n == 2 intersection: " << idata[3] << endl;
              throw(-1);
	    }
	  }
	  else {
	    cout << "ERROR: unsupported n: " << n << endl;
	    throw(-1);
	  }
	}
	
      }

      delete[] edotr;
      delete[] nootr;
      
      // at this point, the link edges have been added between intersectionVec entries

      /*
      if (debug) {
        cout << "intersectionVec.size(): " << intersectionVec.size() << endl;
        FILE * fp = fopen("link.dat","w");
        for (int ied = ned_old; ied < ned; ++ied) {
          const int iv0 = nooed[ied][0] - nno; assert((iv0 >= 0)&&(iv0 < intersectionVec.size()));
          const int iv1 = nooed[ied][1] - nno; assert((iv0 >= 0)&&(iv0 < intersectionVec.size()));
          double xi0[3]; intersectionVec[iv0].calcXi(xi0,x_no,nooed);
          double xi1[3]; intersectionVec[iv1].calcXi(xi1,x_no,nooed);
          writeEdge(fp,xi0,xi1);
        }
        fclose(fp);
        cout << "checkout link.dat" << endl;
      }
      */

      if (!intersectionVec.empty()) {
	  
	// assume each intersection is associated with a new node...
	for (int ii = 0; ii < intersectionVec.size(); ++ii) {
	  intersectionVec[ii].ino = nno+ii;
	}
	sort(intersectionVec.begin(),intersectionVec.end());

	// allocate no_flag...
	const int nno_max = nno+intersectionVec.size(); 
	int * no_flag = new int[nno_max];

	for (int ino = 0; ino < nno; ++ino)
	  no_flag[ino] = ino;
	for (int ino = nno; ino < nno+intersectionVec.size(); ++ino)
	  no_flag[ino] = -1;
    
	// add the nodes...
	vector<pair<pair<int,double>,int> > cutEdgeDataVec;
	vector<int> intersectionNodeNode;

	assert(cutEdgeDataVec.empty());
	assert(intersectionNodeNode.empty());
	int ii_prev = -1;
	for (int ii = 0; ii < intersectionVec.size(); ++ii) {
	  if ((ii_prev == -1)||
	      (intersectionVec[ii].kind != intersectionVec[ii_prev].kind)||
	      (intersectionVec[ii].idata[0] != intersectionVec[ii_prev].idata[0])||
	      (intersectionVec[ii].idata[1] != intersectionVec[ii_prev].idata[1])) {
	    ii_prev = ii;
	    switch (intersectionVec[ii].kind) {
	    case NODE_NODE_INTERSECTION:
	      {
		// combine nodes in idata[0] and idata[1]...
		int ino0 = no_flag[intersectionVec[ii].idata[0]];
		while (ino0 != no_flag[ino0])
		  ino0 = no_flag[ino0];
		int ino1 = no_flag[intersectionVec[ii].idata[1]];
		while (ino1 != no_flag[ino1])
		  ino1 = no_flag[ino1];
		no_flag[intersectionVec[ii].ino] = no_flag[ino0] = no_flag[ino1] = min(ino0,ino1);
		// and add this to a special vec of node-node intersections...
		intersectionNodeNode.push_back(min(ino0,ino1));
	      }
	      break;
	    case NODE_EDGE_INTERSECTION:
	      {
		// use the node to define the x...
		no_flag[intersectionVec[ii].ino] = intersectionVec[ii].idata[0];
		const int ied = intersectionVec[ii].idata[1];
		const double wgt = intersectionVec[ii].ddata[1];
		cutEdgeDataVec.push_back(pair<pair<int,double>,int>(pair<int,double>(ied,wgt),intersectionVec[ii].idata[0]));
	      }
	      break;
	    case NODE_FACE_INTERSECTION:
	      {
		// use the node to define the x...
		no_flag[intersectionVec[ii].ino] = intersectionVec[ii].idata[0];
	      }
	      break;
	    case EDGE_EDGE_INTERSECTION:
	      {
		// this requires a new node. Use the simple average of the two edge intersections...
		const int ied0 = intersectionVec[ii].idata[0];
		const double wgt0 = intersectionVec[ii].ddata[0];
		const int ied1 = intersectionVec[ii].idata[1];
		const double wgt1 = intersectionVec[ii].ddata[1];
		const int ino = new_node();
		FOR_I3 x_no[ino][i] = 0.5*( wgt0*x_no[nooed[ied0][1]][i] + (1.0-wgt0)*x_no[nooed[ied0][0]][i] +
					    wgt1*x_no[nooed[ied1][1]][i] + (1.0-wgt1)*x_no[nooed[ied1][0]][i] );
		no_flag[intersectionVec[ii].ino] = ino;
		cutEdgeDataVec.push_back(pair<pair<int,double>,int>(pair<int,double>(ied0,wgt0),ino));
		cutEdgeDataVec.push_back(pair<pair<int,double>,int>(pair<int,double>(ied1,wgt1),ino));
	      }
	      break;
	    case EDGE_FACE_INTERSECTION:
	      { 
		const int ied = intersectionVec[ii].idata[0];
		const double wgt = intersectionVec[ii].ddata[0];
		const int ino = new_node();
		FOR_I3 x_no[ino][i] = wgt*x_no[nooed[ied][1]][i] + (1.0-wgt)*x_no[nooed[ied][0]][i];
		no_flag[intersectionVec[ii].ino] = ino;
		cutEdgeDataVec.push_back(pair<pair<int,double>,int>(pair<int,double>(ied,wgt),ino));
	      }
	      break;
	    default:
	      assert(0);
	    }
	  }
	  else {
	    assert(ii_prev != -1);
	    no_flag[intersectionVec[ii].ino] = no_flag[intersectionVec[ii_prev].ino];
	  }
	}
	
	intersectionVec.clear();

	// now update nooed for ALL edges to the no_flag value...
	for (int ied = 0; ied < ned; ++ied) {
	  FOR_I2 {
	    const int ino = nooed[ied][i];
	    assert((no_flag[ino] >= 0)&&(no_flag[ino] < nno));
	    nooed[ied][i] = no_flag[ino];
	  }
	}

	delete[] no_flag; no_flag = NULL;
	
	// and cut the edges that have intersections...
	sort(cutEdgeDataVec.begin(),cutEdgeDataVec.end());
	int ied_current = -1;
	int ino_current = -1;
	for (int ii = 0; ii < cutEdgeDataVec.size(); ++ii) {
	  const int ied = cutEdgeDataVec[ii].first.first;
	  if (ied != ied_current) {
	    // this is the first of a new edge. Complete the old edge first...
	    if (ied_current >= 0) {
	      assert(ino_current >= 0);
	      nooed[ied_current][0] = ino_current;
	    }
	    ied_current = ied;
	    ino_current = nooed[ied][0];
	  }
	  //const double s = cutEdgeDataVec[ii].first.second;
	  const int ino_new = cutEdgeDataVec[ii].second;
	  assert(ino_new != ino_current);
	  const int ied_new = new_edge();
	  nooed[ied_new][0] = ino_current;
	  nooed[ied_new][1] = ino_new;
	  faoed[ied_new][0] = faoed[ied_current][0];
	  faoed[ied_new][1] = faoed[ied_current][1];
	  ino_current = ino_new;
	}
	cutEdgeDataVec.clear();
	// complete last edge (the original edge)...
	if (ied_current >= 0) {
	  assert(ino_current >= 0);
	  nooed[ied_current][0] = ino_current;
	}
	
      }
      
    }
    
  }

  if (debug) {
    writeTecplot(1);
    cout << "tecplot 1: This is the surface after intersection, but before cube cutting" << endl;
  }
  
  check();
  
  // now cut against cube boundaries. One important robustness feature here is to 
  // avoid cutting through or near any existing node, so we allow the cube boundary to be 
  // selected somwhat optimally in the range L=0.501*this_delta_cv to 0.509*this_delta_cv.
  // Note that any surface that was grabbed used 0.51*this_delta_cv, and ultimately 
  // the nbrs will use a sphere of r=this_delta_cv (exactly), so this range gives
  // plenty of room to avoid slicing nodes. Note that the creating of new nodes is 
  // also considered by selecting 3 frac's, one for each Cartesian direction...
  
  double frac[3] = { 0.505, 0.505, 0.505 };
  if (ned > 0) {
    // compute frac for x-plane cuts...
    vector<double> fracVec;
    fracVec.push_back(0.501);
    fracVec.push_back(0.509);
    for (int ino = 0; ino < nno; ++ino) {
      if ((x_no[ino][0] > 0.501*this_delta_cv)&&(x_no[ino][0] < 0.509*this_delta_cv)) fracVec.push_back(x_no[ino][0]/this_delta_cv);
      if ((-x_no[ino][0] > 0.501*this_delta_cv)&&(-x_no[ino][0] < 0.509*this_delta_cv)) fracVec.push_back(-x_no[ino][0]/this_delta_cv);
    }
    // then sort...
    if (fracVec.size() > 2) {
      sort(fracVec.begin(),fracVec.end());
      // now find the largest gap in fracVec...
      int imax = 1;
      for (int i = 2; i < fracVec.size(); ++i) {
        if ((fracVec[i]-fracVec[i-1]) > (fracVec[imax]-fracVec[imax-1])) {
          imax = i;
        }
      }
      frac[0] = 0.5*(fracVec[imax]+fracVec[imax-1]);
    }
    // -2 == x0
    double n[3] = { -(frac[0]*this_delta_cv), 0.0, 0.0 };
    cut_surf(n,-2);
    if (ned > 0) {
      // -3 == x1
      n[0] = frac[0]*this_delta_cv;
      cut_surf(n,-3);
      if (ned > 0) {
        // the cuts above have potentially added new nodes, so compute a new frac for y-plane cuts...
        fracVec.clear();
        fracVec.push_back(0.501);
        fracVec.push_back(0.509);
        for (int ino = 0; ino < nno; ++ino) {
          if ((x_no[ino][1] > 0.501*this_delta_cv)&&(x_no[ino][1] < 0.509*this_delta_cv)) fracVec.push_back(x_no[ino][1]/this_delta_cv);
          if ((-x_no[ino][1] > 0.501*this_delta_cv)&&(-x_no[ino][1] < 0.509*this_delta_cv)) fracVec.push_back(-x_no[ino][1]/this_delta_cv);
        }
        // then sort...
        if (fracVec.size() > 2) {
          sort(fracVec.begin(),fracVec.end());
          // now find the largest gap in fracVec...
          int imax = 1;
          for (int i = 2; i < fracVec.size(); ++i) {
            if ((fracVec[i]-fracVec[i-1]) > (fracVec[imax]-fracVec[imax-1])) {
              imax = i;
            }
          }
          frac[1] = 0.5*(fracVec[imax]+fracVec[imax-1]);
        }
        // -4 == y0
        n[0] = 0.0;
        n[1] = -(frac[1]*this_delta_cv); // n[2] still zero
        cut_surf(n,-4);
        if (ned > 0) {
          // -5 == y1
          n[1] = frac[1]*this_delta_cv;
          cut_surf(n,-5);
          if (ned > 0) {
            // and for z-plane cuts...
            fracVec.clear();
            fracVec.push_back(0.501);
            fracVec.push_back(0.509);
            for (int ino = 0; ino < nno; ++ino) {
              if ((x_no[ino][2] > 0.501*this_delta_cv)&&(x_no[ino][2] < 0.509*this_delta_cv)) fracVec.push_back(x_no[ino][2]/this_delta_cv);
              if ((-x_no[ino][2] > 0.501*this_delta_cv)&&(-x_no[ino][2] < 0.509*this_delta_cv)) fracVec.push_back(-x_no[ino][2]/this_delta_cv);
            }
            // then sort...
            if (fracVec.size() > 2) {
              sort(fracVec.begin(),fracVec.end());
              // now find the largest gap in fracVec...
              int imax = 1;
              for (int i = 2; i < fracVec.size(); ++i) {
                if ((fracVec[i]-fracVec[i-1]) > (fracVec[imax]-fracVec[imax-1])) {
                  imax = i;
                }
              }
              frac[2] = 0.5*(fracVec[imax]+fracVec[imax-1]);
            }
            // -6 == z0
            n[1] = 0.0; // n[0] already zero
            n[2] = -(frac[2]*this_delta_cv);
            cut_surf(n,-6);
            if (ned > 0) {
              // -7 == z1
              n[2] = frac[2]*this_delta_cv;
              cut_surf(n,-7);
            }
          }
        }
      }
    }
  }
  
  if (debug) {
    writeTecplot(2);
    cout << "tecplot 2: This is the surface after cube cutting with frac = " << COUT_VEC(frac) << endl;
  }

  // sometimes all of the surface is cut away, so just return a cube as the seed...
  
  completeCube2(frac[0]*this_delta_cv,frac[1]*this_delta_cv,frac[2]*this_delta_cv,debug);

  if (debug) {
    writeTecplot(4);
    cout << "tecplot 4: This is the surface after completeCube2" << endl;
  }

  check();

}

#undef FOR_ITR

// ============================================================================================
// end of new routines
// ============================================================================================

void CuttableVoronoiData::debug() {

  // It is necessary to do a bi-directional walk here...
  
  bool debug = true;
  const int ned_cube = ned; // here we do not rely on the cube edges as indicators
  
  int * ed_flag = new int[ned];
  walkAndFlagEdges(ed_flag,ned_cube,debug);
  
  // at this point, all edges are either +2 (in) or -2 (out),
  // so clean up and return...
  
  int * no_flag = new int[nno];
  FOR_INO no_flag[ino] = -1;
  
  FOR_IED {
    if (ed_flag[ied] == 2) {
      no_flag[nooed[ied][0]] = 0;
      no_flag[nooed[ied][1]] = 0;
    }
  }
  
  const int nno_old = nno;
  nno = 0;
  for (int ino = 0; ino < nno_old; ++ino) {
    if (no_flag[ino] == 0) {
      const int ino_new = nno++;
      FOR_I3 x_no[ino_new][i] = x_no[ino][i];
      no_flag[ino] = ino_new;
    }
    else {
      assert(no_flag[ino] == -1);
    }
  }
  
  const int ned_old = ned;
  ned = 0;
  for (int ied = 0; ied < ned_old; ++ied) {
    if (ed_flag[ied] == 2) {
      const int ied_new = ned++;
      nooed[ied_new][0] = no_flag[nooed[ied][0]]; assert(nooed[ied_new][0] >= 0);
      nooed[ied_new][1] = no_flag[nooed[ied][1]]; assert(nooed[ied_new][1] >= 0);
      faoed[ied_new][0] = faoed[ied][0];
      faoed[ied_new][1] = faoed[ied][1];
    }
  }
  
  delete[] ed_flag;
  delete[] no_flag;
  
  cout << "cvd::debug() done successfully" << endl;
  
}

#include "doit.hpp"

#include "doit2.hpp"

// this routine hijacked for now and put at the top of charles_m.cpp...
/*
void CuttableVoronoiData::doit2(const double L,const int * const ist_and_ipart_buf,const bool debug) {

  doit2(-L,L,-L,L,-L,L,ist_and_ipart_buf,debug);
  
}
*/

int CuttableVoronoiData::triangulateFaceForDoit(const int ifa,vector<pair<int,int> >& faceEdgeVec,const vector<int>& edofa_i) {
  
  // returns the index into the new edges...
	
  const int fe_f = edofa_i[ifa];
  const int ned_loop = edofa_i[ifa+1]-fe_f;
  if (ned_loop == 3) {
    // already a tri. All we need to do (potentially) is to flip the first and second 
    // edge if they do not match. This should put the edges in loop-order. This is 
    // checked later...
    int ino1;
    if (faceEdgeVec[fe_f].second >= 0) {
      ino1 = nooed[faceEdgeVec[fe_f].second][1];
    }
    else {
      ino1 = nooed[-faceEdgeVec[fe_f].second-1][0];
    }
    if (faceEdgeVec[fe_f+1].second >= 0) {
      if (ino1 != nooed[faceEdgeVec[fe_f+1].second][0]) {
	// swap first and second...
	const int tmp = faceEdgeVec[fe_f].second;
	faceEdgeVec[fe_f].second = faceEdgeVec[fe_f+1].second;
	faceEdgeVec[fe_f+1].second = tmp;
      }
    }
    else {
      if (ino1 != nooed[-faceEdgeVec[fe_f+1].second-1][1]) {
	// swap first and second...
	const int tmp = faceEdgeVec[fe_f].second;
	faceEdgeVec[fe_f].second = faceEdgeVec[fe_f+1].second;
	faceEdgeVec[fe_f+1].second = tmp;
      }
    }
    return -1;
  }

  assert(ned_loop > 3);
  
  vector<int> nodeLoop;
  vector<int> edgeLoop;
  int this_faoed = -1;
  int done = 0;
  while (done == 0) {
    for (int fe = edofa_i[ifa]; fe != edofa_i[ifa+1]; ++fe) {
      if (faceEdgeVec[fe].second >= 0) {
	if (this_faoed == -1) {
	  this_faoed = faoed[faceEdgeVec[fe].second][0];
	}
	else {
	  assert(this_faoed == faoed[faceEdgeVec[fe].second][0]);
	}
	if (nodeLoop.empty()) {
	  assert(edgeLoop.empty());
	  nodeLoop.push_back(nooed[faceEdgeVec[fe].second][0]);
	  nodeLoop.push_back(nooed[faceEdgeVec[fe].second][1]);
	  edgeLoop.push_back(faceEdgeVec[fe].second);
	}
	else if (nooed[faceEdgeVec[fe].second][0] == nodeLoop.back()) {
	  edgeLoop.push_back(faceEdgeVec[fe].second);
	  if (nooed[faceEdgeVec[fe].second][1] == nodeLoop.front()) {
	    done = 1;
	    break;
	  }
	  nodeLoop.push_back(nooed[faceEdgeVec[fe].second][1]);
	}
      }
      else {
	if (this_faoed == -1) {
	  this_faoed = faoed[-faceEdgeVec[fe].second-1][1];
	}
	else {
	  assert(this_faoed == faoed[-faceEdgeVec[fe].second-1][1]);
	}
	if (nodeLoop.empty()) {
	  nodeLoop.push_back(nooed[-faceEdgeVec[fe].second-1][1]);
	  nodeLoop.push_back(nooed[-faceEdgeVec[fe].second-1][0]);
	  edgeLoop.push_back(faceEdgeVec[fe].second);
	}
	else if (nooed[-faceEdgeVec[fe].second-1][1] == nodeLoop.back()) {
	  edgeLoop.push_back(faceEdgeVec[fe].second);
	  if (nooed[-faceEdgeVec[fe].second-1][0] == nodeLoop.front()) {
	    done = 1;
	    break;
	  }
	  nodeLoop.push_back(nooed[-faceEdgeVec[fe].second-1][0]);
	}
      }
    }
  }
  // for now, we work with single loops. It might be possible in the future
  // to store the tris more carefully and handle cut parts of cut parts... 
  assert(this_faoed >= 0); // assert(this_faoed < (1<<28)) if using bits to distinguish tri faces from polygon?
  assert(nodeLoop.size() == ned_loop);
  assert(edgeLoop.size() == ned_loop);
  // we are going to fan out from one node. Decide where to add the new 
  // such that the minimum tri size produced is maximized.
  int i0_best = -1;
  double a2_best;
  for (int i0 = 0; i0 < ned_loop; ++i0) {
    // we are going to fan out from the first node...
    double a2_min = HUGE_VAL;
    const int ino0 = nodeLoop[i0];
    for (int i1 = i0+1; i1 < i0+ned_loop-1; ++i1) {
      const int ino1 = nodeLoop[i1%ned_loop];
      const int ino2 = nodeLoop[(i1+1)%ned_loop];
      assert(ino1 != ino0);
      assert(ino2 != ino0);
      const double normal[3] = TRI_NORMAL_2(x_no[ino0],x_no[ino1],x_no[ino2]);
      const double this_a2 = DOT_PRODUCT(normal,normal);
      // we need to find the minimum area. We do not want to fan out
      // from the node that produces a small area...
      if (this_a2 < a2_min)
	a2_min = this_a2;
    }
    if ((i0_best == -1)||(a2_min > a2_best)) {
      i0_best = i0;
      a2_best = a2_min;
    }
  }
  // now repopulate the edges into the faceEdgeVec structure and add new edges. We
  // are also going to flag the new edges with bits to separate their tris
  // and allow cutting. So long as the calling process knows about the bits, they can
  // be linked back to the original ist later...
  faceEdgeVec[fe_f].second   = edgeLoop[i0_best];
  faceEdgeVec[fe_f+1].second = edgeLoop[(i0_best+1)%ned_loop];
  const int ied_f = new_edge();
  faoed[ied_f][0] = -1;
  faoed[ied_f][1] = this_faoed;
  nooed[ied_f][0] = nodeLoop[i0_best];
  nooed[ied_f][1] = nodeLoop[(i0_best+2)%ned_loop];
  // intermediate tris...
  int itri;
  for (itri = 1; itri < ned_loop-3; ++itri) {
    assert(itri < 8); // not necessary if still in a loop, but should look at this...
    assert(faoed[ied_f+itri-1][0] == -1);
    faoed[ied_f+itri-1][0] = this_faoed; // does this face need to be different?: (this_faoed|(itri<<28));
    faceEdgeVec[fe_f+1+itri].second = edgeLoop[(i0_best+1+itri)%ned_loop];
    const int ied = new_edge();
    assert(ied == ied_f+itri);
    faoed[ied_f+itri][0] = -1;
    faoed[ied_f+itri][1] = this_faoed;
    nooed[ied_f+itri][0] = nodeLoop[i0_best];
    nooed[ied_f+itri][1] = nodeLoop[(i0_best+2+itri)%ned_loop];
  }
  // final tri...
  assert(itri < 8);
  assert(faoed[ied_f+itri-1][0] == -1);
  faoed[ied_f+itri-1][0] = this_faoed; // was (this_faoed|(itri<<28));
  faceEdgeVec[fe_f+itri+1].second = edgeLoop[(i0_best+itri+1)%ned_loop];
  faceEdgeVec[fe_f+itri+2].second = edgeLoop[(i0_best+itri+2)%ned_loop];
  assert(fe_f+itri+3 == edofa_i[ifa+1]);
  
  return ied_f;
  
}

void CuttableVoronoiData::packTriEdgeVecForDoit(vector<int>& teVec,const int ifa,const vector<pair<int,int> >& faceEdgeVec,const vector<int>& edofa_i,const int ied_f) {
   
  assert(teVec.empty());
  const int fe_f = edofa_i[ifa];
  const int ned_loop = edofa_i[ifa+1] - fe_f;
  if (ned_loop == 3) {
    assert(ied_f == -1);
    teVec.push_back(faceEdgeVec[fe_f].second);
    teVec.push_back(faceEdgeVec[fe_f+1].second);
    teVec.push_back(faceEdgeVec[fe_f+2].second);
  }
  else {
    assert(ned_loop > 3);
    assert(ied_f >= 0);
    teVec.push_back(faceEdgeVec[fe_f].second);
    teVec.push_back(faceEdgeVec[fe_f+1].second);
    teVec.push_back(-ied_f-1);
    int itri;
    for (itri = 1; itri < ned_loop-3; ++itri) {
      teVec.push_back(ied_f+itri-1);
      teVec.push_back(faceEdgeVec[fe_f+1+itri].second);
      teVec.push_back(-(ied_f+itri)-1);
    }
    teVec.push_back(ied_f+itri-1);
    teVec.push_back(faceEdgeVec[fe_f+itri+1].second);
    teVec.push_back(faceEdgeVec[fe_f+itri+2].second);
    assert(fe_f+itri+3 == edofa_i[ifa+1]);
  }
  assert(teVec.size() == (ned_loop-2)*3);
  
}
 
void CuttableVoronoiData::addEdgeTriIntersectionForDoit(vector<IntersectionData>& localIntersectionVec,
                                                        const int ied0,const int nootr1[3],const int edotr1[3],const int faoed1,
                                                        const int8 tet_vol0,const int8 tet_vol1,
                                                        const int8 tet_vol_no0,const int8 tet_vol_no1,const int8 tet_vol_no2) {
   
  // a few possibilities here now...
  if (tet_vol_no0 == 0) {
    if (tet_vol_no1 == 0) {
      // ied0 is intersecting no2...
      if (tet_vol0 == 0) {
        // we are at ino0 of ied0...
        localIntersectionVec.push_back(IntersectionData(NODE_NODE_INTERSECTION,nooed[ied0][0],nootr1[2]));
      }
      else if (tet_vol1 == 0) {
        // we are at ino1 of ied0...
        localIntersectionVec.push_back(IntersectionData(NODE_NODE_INTERSECTION,nooed[ied0][1],nootr1[2]));
      }
      else {
        localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nootr1[2],ied0,
                                                        double(tet_vol0)/double(tet_vol0+tet_vol1)));
      }
    }
    else if (tet_vol_no2 == 0) {
      // ied0 is intersecting no1...
      if (tet_vol0 == 0) {
        // we are at ino0 of ied0...
        localIntersectionVec.push_back(IntersectionData(NODE_NODE_INTERSECTION,nooed[ied0][0],nootr1[1]));
      }
      else if (tet_vol1 == 0) {
        // we are at ino1 of ied0...
        localIntersectionVec.push_back(IntersectionData(NODE_NODE_INTERSECTION,nooed[ied0][1],nootr1[1]));
      }
      else {
        localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nootr1[1],ied0,
                                                        double(tet_vol0)/double(tet_vol0+tet_vol1)));
      }
    }
    else {
      // ied0 is intersecting edge 1 across from no0...
      if (tet_vol0 == 0) {
        // we are at ino0 of ied0...
        if (edotr1[1] >= 0) {
          // edge1 of tri1 is positively oriented...
          localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nooed[ied0][0],edotr1[1],
                                                          double(tet_vol_no2)/double(tet_vol_no1+tet_vol_no2)));
        }
        else {
          // edge1 of tri1 is negatively oriented...
          localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nooed[ied0][0],-edotr1[1]-1,
                                                          double(tet_vol_no1)/double(tet_vol_no1+tet_vol_no2)));
        }
      }
      else if (tet_vol1 == 0) {
        // we are at ino1 of ied0...
        if (edotr1[1] >= 0) {
          // edge1 of tri1 is positively oriented...
          localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nooed[ied0][1],edotr1[1],
                                                          double(tet_vol_no2)/double(tet_vol_no1+tet_vol_no2)));
        }
        else {
          // edge1 of tri1 is negatively oriented...
          localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nooed[ied0][1],-edotr1[1]-1,
                                                          double(tet_vol_no1)/double(tet_vol_no1+tet_vol_no2)));
        }
      }
      else {
        if (edotr1[1] >= 0) {
          localIntersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,ied0,edotr1[1],
                                                          double(tet_vol0)/double(tet_vol0+tet_vol1),
                                                          double(tet_vol_no2)/double(tet_vol_no1+tet_vol_no2)));
        }
        else {
          localIntersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,ied0,-edotr1[1]-1,
                                                          double(tet_vol0)/double(tet_vol0+tet_vol1),
                                                          double(tet_vol_no1)/double(tet_vol_no1+tet_vol_no2)));
        }
      }
    }
  }
  else if (tet_vol_no1 == 0) {
    // the case of tet_vol_no0 == 0 is captured above
    if (tet_vol_no2 == 0) {
      // ied0 is intersecting no0...
      if (tet_vol0 == 0) {
        // we are at ino0 of ied0...
        localIntersectionVec.push_back(IntersectionData(NODE_NODE_INTERSECTION,nooed[ied0][0],nootr1[0]));
      }
      else if (tet_vol1 == 0) {
        // we are at ino1 of ied0...
        localIntersectionVec.push_back(IntersectionData(NODE_NODE_INTERSECTION,nooed[ied0][1],nootr1[0]));
      }
      else {
        localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nootr1[0],ied0,
                                                        double(tet_vol0)/double(tet_vol0+tet_vol1)));
      }
    }
    else {
      // ied0 is intersecting edge across from no1 (edge 2)...
      if (tet_vol0 == 0) {
        // we are at ino0 of ied0...
        if (edotr1[2] >= 0) {
          // edge2 of tri1 is positively oriented...
          localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nooed[ied0][0],edotr1[2],
                                                          double(tet_vol_no0)/double(tet_vol_no2+tet_vol_no0)));
        }
        else {
          // edge2 of tri1 is negatively oriented...
          localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nooed[ied0][0],-edotr1[2]-1,
                                                          double(tet_vol_no2)/double(tet_vol_no2+tet_vol_no0)));
        }
      }
      else if (tet_vol1 == 0) {
        // we are at ino1 of ied0...
        if (edotr1[2] >= 0) {
          // edge2 of tri1 is positively oriented...
          localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nooed[ied0][1],edotr1[2],
                                                          double(tet_vol_no0)/double(tet_vol_no2+tet_vol_no0)));
        }
        else {
          // edge2 of tri1 is negatively oriented...
          localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nooed[ied0][1],-edotr1[2]-1,
                                                          double(tet_vol_no2)/double(tet_vol_no2+tet_vol_no0)));
        }
      }
      else {
        if (edotr1[2] >= 0) {
          localIntersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,ied0,edotr1[2],
                                                          double(tet_vol0)/double(tet_vol0+tet_vol1),
                                                          double(tet_vol_no0)/double(tet_vol_no2+tet_vol_no0)));
        }
        else {
          localIntersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,ied0,-edotr1[2]-1,
                                                          double(tet_vol0)/double(tet_vol0+tet_vol1),
                                                          double(tet_vol_no2)/double(tet_vol_no2+tet_vol_no0)));
        }
      }
    }
  }
  else if (tet_vol_no2 == 0) {
    // the cases of tet_vol_no0 == 0, and/or tet_vol_no1 == 0 are captured above
    // ied0 is intersecting edge across from no2...
    if (tet_vol0 == 0) {
      // we are at ino0 of ied0...
      if (edotr1[0] >= 0) {
        // edge0 of tri1 is positively oriented...
        localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nooed[ied0][0],edotr1[0],
                                                        double(tet_vol_no1)/double(tet_vol_no0+tet_vol_no1)));
      }
      else {
        // edge0 of tri1 is negatively oriented...
        localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nooed[ied0][0],-edotr1[0]-1,
                                                        double(tet_vol_no0)/double(tet_vol_no0+tet_vol_no1)));
      }
    }
    else if (tet_vol1 == 0) {
      // we are at ino1 of ied0...
      if (edotr1[0] >= 0) {
        // edge0 of tri1 is positively oriented...
        localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nooed[ied0][1],edotr1[0],
                                                        double(tet_vol_no1)/double(tet_vol_no0+tet_vol_no1)));
      }
      else {
        // edge0 of tri1 is negatively oriented...
        localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nooed[ied0][1],-edotr1[0]-1,
                                                        double(tet_vol_no0)/double(tet_vol_no0+tet_vol_no1)));
      }
    }
    else {
      if (edotr1[0] >= 0) {
        localIntersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,ied0,edotr1[0],
                                                        double(tet_vol0)/double(tet_vol0+tet_vol1),
                                                        double(tet_vol_no1)/double(tet_vol_no0+tet_vol_no1)));
      }
      else {
        localIntersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,ied0,-edotr1[0]-1,
                                                        double(tet_vol0)/double(tet_vol0+tet_vol1),
                                                        double(tet_vol_no0)/double(tet_vol_no0+tet_vol_no1)));
      }
    }
  }
  else {
    // all 3 tet volumes are non-zero. edge is intersecting tri...
    if (tet_vol0 == 0) {
      localIntersectionVec.push_back(IntersectionData(NODE_FACE_INTERSECTION,nooed[ied0][0],faoed1));
    }
    else if (tet_vol1 == 0) {
      localIntersectionVec.push_back(IntersectionData(NODE_FACE_INTERSECTION,nooed[ied0][1],faoed1));
    }
    else {
      localIntersectionVec.push_back(IntersectionData(EDGE_FACE_INTERSECTION,ied0,faoed1,double(tet_vol0)/double(tet_vol0+tet_vol1)));
    }
  }
}

int CuttableVoronoiData::processCornerStack(queue<pair<int,int> >& cornerStack,int corner_flag[8][3],vector<pair<double,pair<int,int> > > cutCubeNodeVec[12],const int nno_cube) {
  
  int ierr = 0;
  
  while (!cornerStack.empty()) {

    const pair<int,int> p = cornerStack.front(); cornerStack.pop();
    const int icorner = p.first; assert((icorner >= 0)&&(icorner < 8));
    const int icd = p.second; assert((icd >= 0)&&(icd < 3));
    
    //cout << "PPPPPPPPPPPPPPP popped icorner: " << icorner << " icd: " << icd << endl;
    //writeTecplot(2);

    switch (icorner) {

    case 0: // ================================================
      switch (icd) {
      case 0:
#define ICE 8
#include "ice_forward.hpp"
#undef ICE
	break;
      case 1:
#define ICE 2
#include "ice_forward.hpp"
#undef ICE
	break;
      case 2:
#define ICE 0
#include "ice_forward.hpp"
#undef ICE
	break;
      }
      break;
      
    case 1: // ================================================
      switch (icd) {
      case 0:
#define ICE 8
#include "ice_backward.hpp"
#undef ICE
	break;
      case 1:
#define ICE 6
#include "ice_forward.hpp"
#undef ICE
	break;
      case 2:
#define ICE 4
#include "ice_forward.hpp"
#undef ICE
	break;
      }
      break;

    case 2: // ================================================
      switch (icd) {
      case 0:
#define ICE 10
#include "ice_forward.hpp"
#undef ICE
	break;
      case 1:
#define ICE 2
#include "ice_backward.hpp"
#undef ICE
	break;
      case 2:
#define ICE 1
#include "ice_forward.hpp"
#undef ICE
	break;
      }
      break;

    case 3: // ================================================
      switch (icd) {
      case 0:
#define ICE 10
#include "ice_backward.hpp"
#undef ICE
	break;
      case 1:
#define ICE 6
#include "ice_backward.hpp"
#undef ICE
	break;
      case 2:
#define ICE 5
#include "ice_forward.hpp"
#undef ICE
	break;
      }
      break;
      
    case 4: // ================================================
      switch (icd) {
      case 0:
#define ICE 9
#include "ice_forward.hpp"
#undef ICE
	break;
      case 1:
#define ICE 3
#include "ice_forward.hpp"
#undef ICE
	break;
      case 2:
#define ICE 0
#include "ice_backward.hpp"
#undef ICE
	break;
      default:
	assert(0);
      }
      break;

    case 5: // ================================================
      switch (icd) {
      case 0:
#define ICE 9
#include "ice_backward.hpp"
#undef ICE        
	break;
      case 1:
#define ICE 7
#include "ice_forward.hpp"
#undef ICE        
	break;
      case 2:
#define ICE 4
#include "ice_backward.hpp"
#undef ICE        
	break;
      default:
	assert(0);
      }
      break;

    case 6: // ================================================
      switch (icd) {
      case 0:
#define ICE 11
#include "ice_forward.hpp"
#undef ICE
	break;
      case 1:
#define ICE 3
#include "ice_backward.hpp"
#undef ICE
	break;
      case 2:
#define ICE 1
#include "ice_backward.hpp"
#undef ICE
	break;
      default:
	assert(0);
      }
      break;

    case 7: // ================================================
      switch (icd) {
      case 0:
#define ICE 11
#include "ice_backward.hpp"
#undef ICE
	break;
      case 1:
#define ICE 7
#include "ice_backward.hpp"
#undef ICE
	break;
      case 2:
#define ICE 5
#include "ice_backward.hpp"
#undef ICE
	break;
      default:
	assert(0);
      }
      break;
    default:
      assert(0);
    }
  }
  
  return ierr;
  
}

void CuttableVoronoiData::walkAndFlagEdges(int * ed_flag,const int ned_cube,const bool debug) const {
  
  FOR_IED ed_flag[ied] = 0;
    
  // we need multimaps in both the forward and backward directions. There are
  // cases that require this...
    
  multimap<const pair<int,int>,int> forwardEdgeMultimap;
  multimap<const pair<int,int>,int> backwardEdgeMultimap;
  FOR_IED {
    assert(faoed[ied][0] != faoed[ied][1]);
    assert(nooed[ied][0] != nooed[ied][1]);
    forwardEdgeMultimap.insert(pair<const pair<int,int>,int>(pair<int,int>(faoed[ied][0],nooed[ied][0]),ied));
    forwardEdgeMultimap.insert(pair<const pair<int,int>,int>(pair<int,int>(faoed[ied][1],nooed[ied][1]),ied));
    backwardEdgeMultimap.insert(pair<const pair<int,int>,int>(pair<int,int>(faoed[ied][0],nooed[ied][1]),ied));
    backwardEdgeMultimap.insert(pair<const pair<int,int>,int>(pair<int,int>(faoed[ied][1],nooed[ied][0]),ied));
  }
    
  int loop_count = 0;
  stack<int> edgeStack;
  while (1) {
      
    ++loop_count;
      
    if (debug) cout << "working on loop: " << loop_count << endl;
    
    int ied_next = -1;
    FOR_IED { // note now just ned here
      // skip edges already flagged...
      if (ed_flag[ied] != 0)
	continue;
      ied_next = ied;
      break;
    }
    if (ied_next == -1)
      break;
      
    // assume we are "in" unless we learn otherwise...
    bool in_flag = true;
    assert(ed_flag[ied_next] == 0);
    ed_flag[ied_next] = 1; // provisionally include this edge
      
    assert(edgeStack.empty());
    edgeStack.push(ied_next);
    
    if (debug) {
      FILE * fp = fopen("edges.dat","w");
      fclose(fp);
    }
    
    while (!edgeStack.empty()) {
      
      int ied = edgeStack.top(); edgeStack.pop();
      assert(ed_flag[ied] == 1); // should already be flagged
        
      if (debug) {
	cout << "popping edge " << ied << " with faoed: " << faoed[ied][0] << " " << faoed[ied][1] << endl;
	FILE * fp = fopen("edges.dat","a");
	fprintf(fp,"%18.15le %18.15le %18.15le\n",
		0.5*(x_no[nooed[ied][0]][0]+x_no[nooed[ied][1]][0]),
		0.5*(x_no[nooed[ied][0]][1]+x_no[nooed[ied][1]][1]),
		0.5*(x_no[nooed[ied][0]][2]+x_no[nooed[ied][1]][2]));
	fclose(fp);
      }
    
      // --------------------------------------
      // first face of this edge...
      // --------------------------------------
      int ied_next = getNextEdge(ied,faoed[ied][0],nooed[ied][1],forwardEdgeMultimap);
      if (ied_next >= 0) {
	const int ied_prev = getNextEdge(ied_next,faoed[ied][0],nooed[ied][1],backwardEdgeMultimap);
	if (ied_prev == ied) {
	  // this is a proper, bi-directional match. 
	  if (ed_flag[ied_next] == 0) {
	    ed_flag[ied_next] = 1;
	    edgeStack.push(ied_next);
	  }
	  else {
	    assert(ed_flag[ied_next] == 1);
	  }
	}
	else {
	  // must not have been able to decide...
	  if (ied_prev != -2) {
	    cout << "UNEXPECTED ied_prev: " << ied_prev << ", expecting -2" << endl;
	    //assert(ied_prev == -2);
	  }
	}
      }
      else if (ied_next == -1) {
	// could not find the next edge: must be an edge on the cube face
	assert((faoed[ied][0] >= -7)&&(faoed[ied][0] <= -2));
	in_flag = false;
      }
      else {
	// must not have been able to decide...
	assert(ied_next == -2);
      }
        
      // --------------------------------------
      // second face of this edge...
      // --------------------------------------
      ied_next = getNextEdge(ied,faoed[ied][1],nooed[ied][0],forwardEdgeMultimap);
      if (ied_next >= 0) {
	const int ied_prev = getNextEdge(ied_next,faoed[ied][1],nooed[ied][0],backwardEdgeMultimap);
	if (ied_prev == ied) {
	  // this is a proper, bi-directional match. 
	  if (ed_flag[ied_next] == 0) {
	    ed_flag[ied_next] = 1;
	    edgeStack.push(ied_next);
	  }
	  else {
	    assert(ed_flag[ied_next] == 1);
	  }
	}
	else {
	  // must not have been able to decide...
	  assert(ied_prev == -2);
	}
      }
      else if (ied_next == -1) {
	// could not find the next edge: must be an edge on the cube face
	assert(faoed[ied][1] == -1);
	in_flag = false;
      }
      else {
	// must not have been able to decide...
	assert(ied_next == -2);
      }
      
    } //while (!edgeStack.empty())
      
      // we've done the whole edgeStack. If we made it here, convert all
      // ed_flag == 1 (provisional) to ed_flag = 2 (in) or ed_flag = -2 (out)...
      
      // if we think we are in because the whole surface is walkable, we
      // should still make a couple of checks...
      
    if (in_flag) {
        
      // if any edges selected are cube edges, then we stay in. Recall ned_cube
      // is the start of the cube edges, and they are the last... (I think this 
      // is equivalent to a quick volume check)...
      int ied;
      for (ied = ned_cube; ied < ned; ++ied) {
	if (ed_flag[ied] == 1) {
	  break;
	}
      }
      if (ied == ned) {
	// we had no cube edges, then check if we have links. If we crossed links to 
	// get into this part of the geometry, then we can check the sign of the volume
	// it represents. If we did not cross links, then we need to keep the geometry
	// in general, independent of its volume sign...
	for (ied = 0; ied < ned; ++ied) {
	  if (ed_flag[ied] == 1) {
	    const bool ied_is_link = (((ied > 0)&&(faoed[ied-1][0] == faoed[ied][0])&&(faoed[ied-1][1] == faoed[ied][1])&&(nooed[ied-1][0] == nooed[ied][1])&&(nooed[ied-1][1] == nooed[ied][0]))||
				      ((ied < ned-1)&&(faoed[ied+1][0] == faoed[ied][0])&&(faoed[ied+1][1] == faoed[ied][1])&&(nooed[ied+1][0] == nooed[ied][1])&&(nooed[ied+1][1] == nooed[ied][0])));
	    if (ied_is_link)
	      break;
	  }
	}
	if (ied < ned) {
	  // we have links: check the sign of the volume...
	  double V_sum = 0.0;
	  map<const int,double*> xMap;
	  FOR_IED {
	    if (ed_flag[ied] == 1) {
	      map<const int,double*>::iterator iter = xMap.find(faoed[ied][0]);
	      if (iter == xMap.end()) {
		xMap.insert(pair<int,double*>(faoed[ied][0],x_no[nooed[ied][0]]));
	      }
	      else if (x_no[nooed[ied][1]] != iter->second) {
		V_sum += CROSS_DOT(iter->second,x_no[nooed[ied][0]],x_no[nooed[ied][1]]);
	      }
	      iter = xMap.find(faoed[ied][1]);
	      if (iter == xMap.end()) {
		xMap.insert(pair<int,double*>(faoed[ied][1],x_no[nooed[ied][1]]));
	      }
	      else if (x_no[nooed[ied][0]] != iter->second) {
		V_sum += CROSS_DOT(iter->second,x_no[nooed[ied][1]],x_no[nooed[ied][0]]);
	      }
	    }
	  }
	  if (V_sum < 0.0)
	    in_flag = false;
	}
      }
    }

    if (debug) {
      cout << "just finished loop: " << loop_count << " in_flag: " << in_flag << endl;
      writeTecplot(loop_count);
    }
      
    if (in_flag) {
      FOR_IED {
	if (ed_flag[ied] == 1)
	  ed_flag[ied] = 2;
	else {
	  assert(ed_flag[ied] != -1);
	}
      }
    }
    else {
      // eliminate all flagged edges, AND rewind ned...
      FOR_IED {
	if (ed_flag[ied] == 1)
	  ed_flag[ied] = -2;
	else {
	  assert(ed_flag[ied] != -1);
	}
      }
    }
      
    if (debug) {
      FILE * fp = fopen("edges-p2.dat","w");
      for (int ied = 0; ied < ned; ++ied) {
	if (ed_flag[ied] == 2) {
	  fprintf(fp,"%18.15le %18.15le %18.15le\n",
		  0.5*(x_no[nooed[ied][0]][0]+x_no[nooed[ied][1]][0]),
		  0.5*(x_no[nooed[ied][0]][1]+x_no[nooed[ied][1]][1]),
		  0.5*(x_no[nooed[ied][0]][2]+x_no[nooed[ied][1]][2]));
	}
      }
      fclose(fp);
      fp = fopen("edges-m2.dat","w");
      for (int ied = 0; ied < ned; ++ied) {
	if (ed_flag[ied] == -2) {
	  fprintf(fp,"%18.15le %18.15le %18.15le\n",
		  0.5*(x_no[nooed[ied][0]][0]+x_no[nooed[ied][1]][0]),
		  0.5*(x_no[nooed[ied][0]][1]+x_no[nooed[ied][1]][1]),
		  0.5*(x_no[nooed[ied][0]][2]+x_no[nooed[ied][1]][2]));
	}
      }
      fclose(fp);
      cout << "TAKE A LOOK" << endl;
      if (mpi_size == 1) getchar();
    }
      
  } // while (1)...
    
    // at this point, every edge should be either -2 (not included) or a 2 (included)...
    
  FOR_IED assert((ed_flag[ied] == -2)||(ed_flag[ied] == 2));
    
}

int CuttableVoronoiData::getNextEdge(const int ied,const int ifa,const int ino,const multimap<const pair<int,int>,int>& edgeMultimap) const {
  
  multimap<const pair<int,int>,int>::const_iterator iter = edgeMultimap.find(pair<int,int>(ifa,ino));
  if (iter == edgeMultimap.end()) {
    // did not find edge -- must be unclosed...
    assert((ifa >= -7)&&(ifa <= -1)); // must be an edge on the cube face
    return -1;
  }
  else {
    // we found atleast one edge that fits.
    // there may be more than one edge associated with this face and node...
    multimap<const pair<int,int>,int>::const_iterator iter2 = iter;
    ++iter2;
    if ((iter2 != edgeMultimap.end())&&(iter2->first == pair<int,int>(ifa,ino))) {
      // this is a 2-edge node.
      multimap<const pair<int,int>,int>::const_iterator iter_check = iter2;
      ++iter_check;
      assert((iter_check == edgeMultimap.end())||(iter_check->first != pair<int,int>(ifa,ino)));
      const int ied1 = iter->second;
      const int ied2 = iter2->second;
      assert(ied1 != ied2);
      if ((faoed[ied1][0] == faoed[ied][0])&&(faoed[ied1][1] == faoed[ied][1])) {
	assert(!((faoed[ied2][0] == faoed[ied][0])&&(faoed[ied2][1] == faoed[ied][1])));
	return ied2;
      }
      else if ((faoed[ied2][0] == faoed[ied][0])&&(faoed[ied2][1] == faoed[ied][1])) {
	return ied1;
      }
      else {
	// Note that the above check takes the link when appropriate, and avoids back-tracing
	// along the link while ied is the link. If the above does not work and we get here, then 
	// I think we need to choose the link...
	const bool ied1_is_link = (((ied1 > 0)&&(faoed[ied1-1][0] == faoed[ied1][0])&&(faoed[ied1-1][1] == faoed[ied1][1])&&(nooed[ied1-1][0] == nooed[ied1][1])&&(nooed[ied1-1][1] == nooed[ied1][0]))||
				   ((ied1 < ned-1)&&(faoed[ied1+1][0] == faoed[ied1][0])&&(faoed[ied1+1][1] == faoed[ied1][1])&&(nooed[ied1+1][0] == nooed[ied1][1])&&(nooed[ied1+1][1] == nooed[ied1][0])));
	const bool ied2_is_link = (((ied2 > 0)&&(faoed[ied2-1][0] == faoed[ied2][0])&&(faoed[ied2-1][1] == faoed[ied2][1])&&(nooed[ied2-1][0] == nooed[ied2][1])&&(nooed[ied2-1][1] == nooed[ied2][0]))||
				   ((ied2 < ned-1)&&(faoed[ied2+1][0] == faoed[ied2][0])&&(faoed[ied2+1][1] == faoed[ied2][1])&&(nooed[ied2+1][0] == nooed[ied2][1])&&(nooed[ied2+1][1] == nooed[ied2][0])));
	if (ied1_is_link && (!ied2_is_link)) {
	  return ied1;
	}
	else if ((!ied1_is_link) && ied2_is_link) {
	  return ied2;
	}
	else {
	  // I don't think we can make this decision. So don't, and hope
	  // the part of the geometry is traversed by other walks...
	  return -2;
	}
      }
    }
    else {
      // just found one linked edge, so return it...
      return iter->second;
    }
  }

}

