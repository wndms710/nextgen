
void tetCutPlane(vector<SimpleTri>& triVec,
		 const double phi0,const double phi1,const double phi2,const double phi3,
		 const double x0[3],const double x1[3],const double x2[3],const double x3[3]) {

  if (phi0 >= 0.0) {
    if (phi1 >= 0.0) {
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
	  //cout << " case a " << endl;
	  return;
	}
	else {
	  //int debug_int = 1;
  	  //cout << " case b " << endl;
	  // phi3 is the only point on the minus side...
#define PHIM phi3
#define PHIP0 phi0
#define PHIP1 phi2
#define PHIP2 phi1
#define XM x3
#define XP0 x0
#define XP1 x2
#define XP2 x1
#include "iso_one_point_minus.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
	}
      }
      else {
	// phi2 is negative...
	if (phi3 >= 0.0) {
	  //int debug_int = 2;
	  //cout << " case c " << endl;
	  // phi2 is the only minus...
#define PHIM phi2
#define PHIP0 phi0
#define PHIP1 phi1
#define PHIP2 phi3
#define XM x2
#define XP0 x0
#define XP1 x1
#define XP2 x3
#include "iso_one_point_minus.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
	}
	else {
	  //int debug_int = 3;
	  //cout << " case d " << endl;
	  // both phi2 and phi3 are negative...
#define PHIM0 phi2
#define PHIM1 phi3
#define PHIP0 phi0
#define PHIP1 phi1
#define XM0 x2
#define XM1 x3
#define XP0 x0
#define XP1 x1
#include "iso_two_point_minus.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
	}
      }
    }
    else {
      // phi1 is minus...
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
	  //int debug_int = 4;
	  //cout << " case e " << endl;
	  // phi1 is the only minus...
#define PHIM phi1
#define PHIP0 phi0
#define PHIP1 phi3
#define PHIP2 phi2
#define XM x1
#define XP0 x0
#define XP1 x3
#define XP2 x2
#include "iso_one_point_minus.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
	}
	else {
	  //int debug_int = 5;
	  //cout << " case f " << endl;
	  // phi1 and phi3 are both minus...
#define PHIM0 phi1
#define PHIM1 phi3
#define PHIP0 phi2
#define PHIP1 phi0
#define XM0 x1
#define XM1 x3
#define XP0 x2
#define XP1 x0
#include "iso_two_point_minus.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
	  //int debug_int = 6;
	  //cout << " case g " << endl;
	  // phi1 and phi2 are both minus...
#define PHIM0 phi1
#define PHIM1 phi2
#define PHIP0 phi0
#define PHIP1 phi3
#define XM0 x1
#define XM1 x2
#define XP0 x0
#define XP1 x3
#include "iso_two_point_minus.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
	}
	else {
	  //int debug_int = 7;
	  //cout << " case h " << endl;
	  // phi1, phi2, and phi3 are minus...
#define PHIP phi0
#define PHIM0 phi1
#define PHIM1 phi3
#define PHIM2 phi2
#define XP x0
#define XM0 x1
#define XM1 x3
#define XM2 x2
#include "iso_three_point_minus.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
	}
      }
    }
  }
  else {
    // phi0 is negative...
    if (phi1 >= 0.0) {
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
	  //int debug_int = 8;
	  //cout << " case i " << endl;
	  // only phi0 is negative...
#define PHIM phi0
#define PHIP0 phi1
#define PHIP1 phi2
#define PHIP2 phi3
#define XM x0
#define XP0 x1
#define XP1 x2
#define XP2 x3
#include "iso_one_point_minus.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
	}
	else {
	  //int debug_int = 9;
	  //cout << " case j " << endl;
	  // both phi0 and phi3 are minus...
#define PHIM0 phi0
#define PHIM1 phi3
#define PHIP0 phi1
#define PHIP1 phi2
#define XM0 x0
#define XM1 x3
#define XP0 x1
#define XP1 x2
#include "iso_two_point_minus.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
	  //int debug_int = 10;
	  //cout << " case k " << endl;
	  // both phi0 and phi2 are minus...
#define PHIM0 phi0
#define PHIM1 phi2
#define PHIP0 phi3
#define PHIP1 phi1
#define XM0 x0
#define XM1 x2
#define XP0 x3
#define XP1 x1
#include "iso_two_point_minus.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
	}
	else {
	  //int debug_int = 11;
	  //cout << " case l " << endl;
	  // phi0, phi2, phi3 all minus...
#define PHIP phi1
#define PHIM0 phi0
#define PHIM1 phi2
#define PHIM2 phi3
#define XP x1
#define XM0 x0
#define XM1 x2
#define XM2 x3
#include "iso_three_point_minus.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
	}
      }
    }
    else {
      // phi1 is minus...
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
	  //int debug_int = 12;
	  //cout << " case m " << endl;
	  // phi0 and phi1 are minus...
#define PHIM0 phi0
#define PHIM1 phi1
#define PHIP0 phi2
#define PHIP1 phi3
#define XM0 x0
#define XM1 x1
#define XP0 x2
#define XP1 x3
#include "iso_two_point_minus.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
	}
	else {
	  //int debug_int = 13;
	  //cout << " case n " << endl;
	  // phi0 and phi1 and phi3 are minus...
#define PHIP phi2
#define PHIM0 phi0
#define PHIM1 phi3
#define PHIM2 phi1
#define XP x2
#define XM0 x0
#define XM1 x3
#define XM2 x1
#include "iso_three_point_minus.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
	  //int debug_int = 14;
	  //cout << " case o " << endl;
	  // phi0,phi1,and phi2 are all minus...
#define PHIP phi3
#define PHIM0 phi0
#define PHIM1 phi1
#define PHIM2 phi2
#define XP x3
#define XM0 x0
#define XM1 x1
#define XM2 x2
#include "iso_three_point_minus.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
	}
	else {
	  //cout << " case p " << endl;
	  // everybody minus...
	  assert(phi0 < 0.0);
	  assert(phi1 < 0.0);
	  assert(phi2 < 0.0);
	  assert(phi3 < 0.0);
	  return;
	}
      }
    }
  }
}


void tetCutPlane(vector<pair<SimpleTri,int> >& triVec,
		 const double phi0,const double phi1,const double phi2,const double phi3,
		 const double x0[3],const double x1[3],const double x2[3],const double x3[3]) {

  if (phi0 >= 0.0) {
    if (phi1 >= 0.0) {
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
	  //cout << " case a " << endl;
	  return;
	}
	else {
	  //int debug_int = 1;
  	  //cout << " case b " << endl;
	  // phi3 is the only point on the minus side...
#define PHIM phi3
#define PHIP0 phi0
#define PHIP1 phi2
#define PHIP2 phi1
#define XM x3
#define XP0 x0
#define XP1 x2
#define XP2 x1
#define MESH_EDGES 1<<1
#include "iso_one_point_minus_mesh.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef MESH_EDGES
	}
      }
      else {
	// phi2 is negative...
	if (phi3 >= 0.0) {
	  //int debug_int = 2;
	  //cout << " case c " << endl;
	  // phi2 is the only minus...
#define PHIM phi2
#define PHIP0 phi0
#define PHIP1 phi1
#define PHIP2 phi3
#define XM x2
#define XP0 x0
#define XP1 x1
#define XP2 x3
#define MESH_EDGES 1<<1
#include "iso_one_point_minus_mesh.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef MESH_EDGES
	}
	else {
	  //int debug_int = 3;
	  //cout << " case d " << endl;
	  // both phi2 and phi3 are negative...
#define PHIM0 phi2
#define PHIM1 phi3
#define PHIP0 phi0
#define PHIP1 phi1
#define XM0 x2
#define XM1 x3
#define XP0 x0
#define XP1 x1
#define MESH_EDGES0 -1
#define MESH_EDGES1 1<<0
#include "iso_two_point_minus_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
      }
    }
    else {
      // phi1 is minus...
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
	  //int debug_int = 4;
	  //cout << " case e " << endl;
	  // phi1 is the only minus...
#define PHIM phi1
#define PHIP0 phi0
#define PHIP1 phi3
#define PHIP2 phi2
#define XM x1
#define XP0 x0
#define XP1 x3
#define XP2 x2
#define MESH_EDGES 1<<1
#include "iso_one_point_minus_mesh.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef MESH_EDGES
	}
	else {
	  //int debug_int = 5;
	  //cout << " case f " << endl;
	  // phi1 and phi3 are both minus...
#define PHIM0 phi1
#define PHIM1 phi3
#define PHIP0 phi2
#define PHIP1 phi0
#define XM0 x1
#define XM1 x3
#define XP0 x2
#define XP1 x0
#define MESH_EDGES0 1<<0
#define MESH_EDGES1 -1
#include "iso_two_point_minus_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
	  //int debug_int = 6;
	  //cout << " case g " << endl;
	  // phi1 and phi2 are both minus...
#define PHIM0 phi1
#define PHIM1 phi2
#define PHIP0 phi0
#define PHIP1 phi3
#define XM0 x1
#define XM1 x2
#define XP0 x0
#define XP1 x3
#define MESH_EDGES0 -1
#define MESH_EDGES1 1<<0
#include "iso_two_point_minus_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
	else {
	  //int debug_int = 7;
	  //cout << " case h " << endl;
	  // phi1, phi2, and phi3 are minus...
#define PHIP phi0
#define PHIM0 phi1
#define PHIM1 phi3
#define PHIM2 phi2
#define XP x0
#define XM0 x1
#define XM1 x3
#define XM2 x2
#define MESH_EDGES -1
#include "iso_three_point_minus_mesh.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef MESH_EDGES
	}
      }
    }
  }
  else {
    // phi0 is negative...
    if (phi1 >= 0.0) {
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
	  //int debug_int = 8;
	  //cout << " case i " << endl;
	  // only phi0 is negative...
#define PHIM phi0
#define PHIP0 phi1
#define PHIP1 phi2
#define PHIP2 phi3
#define XM x0
#define XP0 x1
#define XP1 x2
#define XP2 x3
#define MESH_EDGES -1
#include "iso_one_point_minus_mesh.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef MESH_EDGES
	}
	else {
	  //int debug_int = 9;
	  //cout << " case j " << endl;
	  // both phi0 and phi3 are minus...
#define PHIM0 phi0
#define PHIM1 phi3
#define PHIP0 phi1
#define PHIP1 phi2
#define XM0 x0
#define XM1 x3
#define XP0 x1
#define XP1 x2
#define MESH_EDGES0 1<<2
#define MESH_EDGES1 -1
#include "iso_two_point_minus_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
	  //int debug_int = 10;
	  //cout << " case k " << endl;
	  // both phi0 and phi2 are minus...
#define PHIM0 phi0
#define PHIM1 phi2
#define PHIP0 phi3
#define PHIP1 phi1
#define XM0 x0
#define XM1 x2
#define XP0 x3
#define XP1 x1
#define MESH_EDGES0 1<<2
#define MESH_EDGES1 -1
#include "iso_two_point_minus_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
	else {
	  //int debug_int = 11;
	  //cout << " case l " << endl;
	  // phi0, phi2, phi3 all minus...
#define PHIP phi1
#define PHIM0 phi0
#define PHIM1 phi2
#define PHIM2 phi3
#define XP x1
#define XM0 x0
#define XM1 x2
#define XM2 x3
#define MESH_EDGES 1<<1
#include "iso_three_point_minus_mesh.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef MESH_EDGES
	}
      }
    }
    else {
      // phi1 is minus...
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
	  //int debug_int = 12;
	  //cout << " case m " << endl;
	  // phi0 and phi1 are minus...
#define PHIM0 phi0
#define PHIM1 phi1
#define PHIP0 phi2
#define PHIP1 phi3
#define XM0 x0
#define XM1 x1
#define XP0 x2
#define XP1 x3
#define MESH_EDGES0 1<<2
#define MESH_EDGES1 -1
#include "iso_two_point_minus_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
	else {
	  //int debug_int = 13;
	  //cout << " case n " << endl;
	  // phi0 and phi1 and phi3 are minus...
#define PHIP phi2
#define PHIM0 phi0
#define PHIM1 phi3
#define PHIM2 phi1
#define XP x2
#define XM0 x0
#define XM1 x3
#define XM2 x1
#define MESH_EDGES 1<<1
#include "iso_three_point_minus_mesh.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef MESH_EDGES
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
	  //int debug_int = 14;
	  //cout << " case o " << endl;
	  // phi0,phi1,and phi2 are all minus...
#define PHIP phi3
#define PHIM0 phi0
#define PHIM1 phi1
#define PHIM2 phi2
#define XP x3
#define XM0 x0
#define XM1 x1
#define XM2 x2
#define MESH_EDGES 1<<1
#include "iso_three_point_minus_mesh.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef MESH_EDGES
	}
	else {
	  //cout << " case p " << endl;
	  // everybody minus...
	  assert(phi0 < 0.0);
	  assert(phi1 < 0.0);
	  assert(phi2 < 0.0);
	  assert(phi3 < 0.0);
	  return;
	}
      }
    }
  }
}

void chtTetCutPlane(vector<pair<SimpleTri,int> >& triVec,
		 const double phi0,const double phi1,const double phi2,const double phi3,
		 const double x0[3],const double x1[3],const double x2[3],const double x3[3]) {

  if (phi0 >= 0.0) {
    if (phi1 >= 0.0) {
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
	  //cout << " case a " << endl;
	  return;
	}
	else {
	  //int debug_int = 1;
  	  //cout << " case b " << endl;
	  // phi3 is the only point on the minus side...
#define PHIM phi3
#define PHIP0 phi0
#define PHIP1 phi2
#define PHIP2 phi1
#define XM x3
#define XP0 x0
#define XP1 x2
#define XP2 x1
#define MESH_EDGES 7
#include "iso_one_point_minus_mesh.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef MESH_EDGES
	}
      }
      else {
	// phi2 is negative...
	if (phi3 >= 0.0) {
	  //int debug_int = 2;
	  //cout << " case c " << endl;
	  // phi2 is the only minus...
#define PHIM phi2
#define PHIP0 phi0
#define PHIP1 phi1
#define PHIP2 phi3
#define XM x2
#define XP0 x0
#define XP1 x1
#define XP2 x3
#define MESH_EDGES 7
#include "iso_one_point_minus_mesh.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef MESH_EDGES
	}
	else {
	  //int debug_int = 3;
	  //cout << " case d " << endl;
	  // both phi2 and phi3 are negative...
#define PHIM0 phi2
#define PHIM1 phi3
#define PHIP0 phi0
#define PHIP1 phi1
#define XM0 x2
#define XM1 x3
#define XP0 x0
#define XP1 x1
#define MESH_EDGES0 5
#define MESH_EDGES1 5
#include "iso_two_point_minus_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
      }
    }
    else {
      // phi1 is minus...
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
	  //int debug_int = 4;
	  //cout << " case e " << endl;
	  // phi1 is the only minus...
#define PHIM phi1
#define PHIP0 phi0
#define PHIP1 phi3
#define PHIP2 phi2
#define XM x1
#define XP0 x0
#define XP1 x3
#define XP2 x2
#define MESH_EDGES 7
#include "iso_one_point_minus_mesh.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef MESH_EDGES
	}
	else {
	  //int debug_int = 5;
	  //cout << " case f " << endl;
	  // phi1 and phi3 are both minus...
#define PHIM0 phi1
#define PHIM1 phi3
#define PHIP0 phi2
#define PHIP1 phi0
#define XM0 x1
#define XM1 x3
#define XP0 x2
#define XP1 x0
#define MESH_EDGES0 5
#define MESH_EDGES1 5
#include "iso_two_point_minus_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
	  //int debug_int = 6;
	  //cout << " case g " << endl;
	  // phi1 and phi2 are both minus...
#define PHIM0 phi1
#define PHIM1 phi2
#define PHIP0 phi0
#define PHIP1 phi3
#define XM0 x1
#define XM1 x2
#define XP0 x0
#define XP1 x3
#define MESH_EDGES0 5
#define MESH_EDGES1 5
#include "iso_two_point_minus_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
	else {
	  //int debug_int = 7;
	  //cout << " case h " << endl;
	  // phi1, phi2, and phi3 are minus...
#define PHIP phi0
#define PHIM0 phi1
#define PHIM1 phi3
#define PHIM2 phi2
#define XP x0
#define XM0 x1
#define XM1 x3
#define XM2 x2
#define MESH_EDGES 7
#include "iso_three_point_minus_mesh.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef MESH_EDGES
	}
      }
    }
  }
  else {
    // phi0 is negative...
    if (phi1 >= 0.0) {
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
	  //int debug_int = 8;
	  //cout << " case i " << endl;
	  // only phi0 is negative...
#define PHIM phi0
#define PHIP0 phi1
#define PHIP1 phi2
#define PHIP2 phi3
#define XM x0
#define XP0 x1
#define XP1 x2
#define XP2 x3
#define MESH_EDGES 7
#include "iso_one_point_minus_mesh.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef MESH_EDGES
	}
	else {
	  //int debug_int = 9;
	  //cout << " case j " << endl;
	  // both phi0 and phi3 are minus...
#define PHIM0 phi0
#define PHIM1 phi3
#define PHIP0 phi1
#define PHIP1 phi2
#define XM0 x0
#define XM1 x3
#define XP0 x1
#define XP1 x2
#define MESH_EDGES0 5
#define MESH_EDGES1 5
#include "iso_two_point_minus_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
	  //int debug_int = 10;
	  //cout << " case k " << endl;
	  // both phi0 and phi2 are minus...
#define PHIM0 phi0
#define PHIM1 phi2
#define PHIP0 phi3
#define PHIP1 phi1
#define XM0 x0
#define XM1 x2
#define XP0 x3
#define XP1 x1
#define MESH_EDGES0 5
#define MESH_EDGES1 5
#include "iso_two_point_minus_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
	else {
	  //int debug_int = 11;
	  //cout << " case l " << endl;
	  // phi0, phi2, phi3 all minus...
#define PHIP phi1
#define PHIM0 phi0
#define PHIM1 phi2
#define PHIM2 phi3
#define XP x1
#define XM0 x0
#define XM1 x2
#define XM2 x3
#define MESH_EDGES 7
#include "iso_three_point_minus_mesh.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef MESH_EDGES
	}
      }
    }
    else {
      // phi1 is minus...
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
	  //int debug_int = 12;
	  //cout << " case m " << endl;
	  // phi0 and phi1 are minus...
#define PHIM0 phi0
#define PHIM1 phi1
#define PHIP0 phi2
#define PHIP1 phi3
#define XM0 x0
#define XM1 x1
#define XP0 x2
#define XP1 x3
#define MESH_EDGES0 5
#define MESH_EDGES1 5
#include "iso_two_point_minus_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
	else {
	  //int debug_int = 13;
	  //cout << " case n " << endl;
	  // phi0 and phi1 and phi3 are minus...
#define PHIP phi2
#define PHIM0 phi0
#define PHIM1 phi3
#define PHIM2 phi1
#define XP x2
#define XM0 x0
#define XM1 x3
#define XM2 x1
#define MESH_EDGES 7
#include "iso_three_point_minus_mesh.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef MESH_EDGES
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
	  //int debug_int = 14;
	  //cout << " case o " << endl;
	  // phi0,phi1,and phi2 are all minus...
#define PHIP phi3
#define PHIM0 phi0
#define PHIM1 phi1
#define PHIM2 phi2
#define XP x3
#define XM0 x0
#define XM1 x1
#define XM2 x2
#define MESH_EDGES 7
#include "iso_three_point_minus_mesh.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef MESH_EDGES
	}
	else {
	  //cout << " case p " << endl;
	  // everybody minus...
	  assert(phi0 < 0.0);
	  assert(phi1 < 0.0);
	  assert(phi2 < 0.0);
	  assert(phi3 < 0.0);
	  return;
	}
      }
    }
  }
}
