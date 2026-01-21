
void tetCutPlaneData(vector<SimpleTriWithData>& triVec,
		     const double phi0,const double phi1,const double phi2,const double phi3,
		     const double x0[3],const double x1[3],const double x2[3],const double x3[3],
		     const double d0,const double d1,const double d2,const double d3) {

  if (phi0 >= 0.0) {
    if (phi1 >= 0.0) {
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
	  //cout << " case a " << endl;
	  return;
	}
	else {
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
#define DM d3
#define DP0 d0
#define DP1 d2
#define DP2 d1
#include "iso_one_point_minus_data.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef DM
#undef DP0
#undef DP1
#undef DP2
	}
      }
      else {
	// phi2 is negative...
	if (phi3 >= 0.0) {
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
#define DM d2
#define DP0 d0
#define DP1 d1
#define DP2 d3
#include "iso_one_point_minus_data.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef DM
#undef DP0
#undef DP1
#undef DP2
	}
	else {
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
#define DM0 d2
#define DM1 d3
#define DP0 d0
#define DP1 d1
#include "iso_two_point_minus_data.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
	}
      }
    }
    else {
      // phi1 is minus...
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
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
#define DM d1
#define DP0 d0
#define DP1 d3
#define DP2 d2
#include "iso_one_point_minus_data.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef DM
#undef DP0
#undef DP1
#undef DP2
	}
	else {
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
#define DM0 d1
#define DM1 d3
#define DP0 d2
#define DP1 d0
#include "iso_two_point_minus_data.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
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
#define DM0 d1
#define DM1 d2
#define DP0 d0
#define DP1 d3
#include "iso_two_point_minus_data.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
	}
	else {
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
#define DP d0
#define DM0 d1
#define DM1 d3
#define DM2 d2
#include "iso_three_point_minus_data.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef DP
#undef DM0
#undef DM1
#undef DM2
	}
      }
    }
  }
  else {
    // phi0 is negative...
    if (phi1 >= 0.0) {
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
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
#define DM d0
#define DP0 d1
#define DP1 d2
#define DP2 d3
#include "iso_one_point_minus_data.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef DM
#undef DP0
#undef DP1
#undef DP2
	}
	else {
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
#define DM0 d0
#define DM1 d3
#define DP0 d1
#define DP1 d2
#include "iso_two_point_minus_data.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
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
#define DM0 d0
#define DM1 d2
#define DP0 d3
#define DP1 d1
#include "iso_two_point_minus_data.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
	}
	else {
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
#define DP d1
#define DM0 d0
#define DM1 d2
#define DM2 d3
#include "iso_three_point_minus_data.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef DP
#undef DM0
#undef DM1
#undef DM2
	}
      }
    }
    else {
      // phi1 is minus...
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
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
#define DM0 d0
#define DM1 d1
#define DP0 d2
#define DP1 d3
#include "iso_two_point_minus_data.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
	}
	else {
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
#define DP d2
#define DM0 d0
#define DM1 d3
#define DM2 d1
#include "iso_three_point_minus_data.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef DP
#undef DM0
#undef DM1
#undef DM2
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
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
#define DP d3
#define DM0 d0
#define DM1 d1
#define DM2 d2
#include "iso_three_point_minus_data.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef DP
#undef DM0
#undef DM1
#undef DM2
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

// int paired with SimpleTri determines which edge of this tri to draw as a valid mesh edge
void tetCutPlaneData(vector<pair<SimpleTriWithData,int> >& triVec,
		     const double phi0,const double phi1,const double phi2,const double phi3,
		     const double x0[3],const double x1[3],const double x2[3],const double x3[3],
		     const double d0,const double d1,const double d2,const double d3) {

  if (phi0 >= 0.0) {
    if (phi1 >= 0.0) {
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
	  //cout << " case a " << endl;
	  return;
	}
	else {
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
#define DM d3
#define DP0 d0
#define DP1 d2
#define DP2 d1
#define MESH_EDGES 1<<1
#include "iso_one_point_minus_data_mesh.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef DM
#undef DP0
#undef DP1
#undef DP2
#undef MESH_EDGES
	}
      }
      else {
	// phi2 is negative...
	if (phi3 >= 0.0) {
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
#define DM d2
#define DP0 d0
#define DP1 d1
#define DP2 d3
#define MESH_EDGES 1<<1
#include "iso_one_point_minus_data_mesh.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef DM
#undef DP0
#undef DP1
#undef DP2
#undef MESH_EDGES
	}
	else {
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
#define DM0 d2
#define DM1 d3
#define DP0 d0
#define DP1 d1
#define MESH_EDGES0 -1
#define MESH_EDGES1 1<<0
#include "iso_two_point_minus_data_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
      }
    }
    else {
      // phi1 is minus...
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
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
#define DM d1
#define DP0 d0
#define DP1 d3
#define DP2 d2
#define MESH_EDGES 1<<1
#include "iso_one_point_minus_data_mesh.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef DM
#undef DP0
#undef DP1
#undef DP2
#undef MESH_EDGES
	}
	else {
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
#define DM0 d1
#define DM1 d3
#define DP0 d2
#define DP1 d0
#define MESH_EDGES0 1<<0
#define MESH_EDGES1 -1
#include "iso_two_point_minus_data_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
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
#define DM0 d1
#define DM1 d2
#define DP0 d0
#define DP1 d3
#define MESH_EDGES0 -1
#define MESH_EDGES1 1<<0
#include "iso_two_point_minus_data_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
	else {
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
#define DP d0
#define DM0 d1
#define DM1 d3
#define DM2 d2
#define MESH_EDGES -1
#include "iso_three_point_minus_data_mesh.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef DP
#undef DM0
#undef DM1
#undef DM2
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
#define DM d0
#define DP0 d1
#define DP1 d2
#define DP2 d3
#define MESH_EDGES -1
#include "iso_one_point_minus_data_mesh.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef DM
#undef DP0
#undef DP1
#undef DP2
#undef MESH_EDGES
	}
	else {
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
#define DM0 d0
#define DM1 d3
#define DP0 d1
#define DP1 d2
#define MESH_EDGES0 1<<2
#define MESH_EDGES1 -1
#include "iso_two_point_minus_data_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
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
#define DM0 d0
#define DM1 d2
#define DP0 d3
#define DP1 d1
#define MESH_EDGES0 1<<2
#define MESH_EDGES1 -1
#include "iso_two_point_minus_data_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
	else {
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
#define DP d1
#define DM0 d0
#define DM1 d2
#define DM2 d3
#define MESH_EDGES 1<<1
#include "iso_three_point_minus_data_mesh.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef DP
#undef DM0
#undef DM1
#undef DM2
#undef MESH_EDGES
	}
      }
    }
    else {
      // phi1 is minus...
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
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
#define DM0 d0
#define DM1 d1
#define DP0 d2
#define DP1 d3
#define MESH_EDGES0 1<<2
#define MESH_EDGES1 -1
#include "iso_two_point_minus_data_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
	else {
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
#define DP d2
#define DM0 d0
#define DM1 d3
#define DM2 d1
#define MESH_EDGES 1<<1
#include "iso_three_point_minus_data_mesh.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef DP
#undef DM0
#undef DM1
#undef DM2
#undef MESH_EDGES
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
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
#define DP d3
#define DM0 d0
#define DM1 d1
#define DM2 d2
#define MESH_EDGES 1<<1
#include "iso_three_point_minus_data_mesh.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef DP
#undef DM0
#undef DM1
#undef DM2
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

// int paired with SimpleTri determines which edge of this tri to draw as a valid mesh edge
void chtTetCutPlaneData(vector<pair<SimpleTriWithData,int> >& triVec,
		     const double phi0,const double phi1,const double phi2,const double phi3,
		     const double x0[3],const double x1[3],const double x2[3],const double x3[3],
		     const double d0,const double d1,const double d2,const double d3) {

  if (phi0 >= 0.0) {
    if (phi1 >= 0.0) {
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
	  //cout << " case a " << endl;
	  return;
	}
	else {
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
#define DM d3
#define DP0 d0
#define DP1 d2
#define DP2 d1
#define MESH_EDGES 7
#include "iso_one_point_minus_data_mesh.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef DM
#undef DP0
#undef DP1
#undef DP2
#undef MESH_EDGES
	}
      }
      else {
	// phi2 is negative...
	if (phi3 >= 0.0) {
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
#define DM d2
#define DP0 d0
#define DP1 d1
#define DP2 d3
#define MESH_EDGES 7
#include "iso_one_point_minus_data_mesh.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef DM
#undef DP0
#undef DP1
#undef DP2
#undef MESH_EDGES
	}
	else {
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
#define DM0 d2
#define DM1 d3
#define DP0 d0
#define DP1 d1
#define MESH_EDGES0 5
#define MESH_EDGES1 5
#include "iso_two_point_minus_data_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
      }
    }
    else {
      // phi1 is minus...
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
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
#define DM d1
#define DP0 d0
#define DP1 d3
#define DP2 d2
#define MESH_EDGES 7
#include "iso_one_point_minus_data_mesh.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef DM
#undef DP0
#undef DP1
#undef DP2
#undef MESH_EDGES
	}
	else {
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
#define DM0 d1
#define DM1 d3
#define DP0 d2
#define DP1 d0
#define MESH_EDGES0 5
#define MESH_EDGES1 5
#include "iso_two_point_minus_data_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
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
#define DM0 d1
#define DM1 d2
#define DP0 d0
#define DP1 d3
#define MESH_EDGES0 5
#define MESH_EDGES1 5
#include "iso_two_point_minus_data_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
	else {
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
#define DP d0
#define DM0 d1
#define DM1 d3
#define DM2 d2
#define MESH_EDGES 7
#include "iso_three_point_minus_data_mesh.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef DP
#undef DM0
#undef DM1
#undef DM2
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
#define DM d0
#define DP0 d1
#define DP1 d2
#define DP2 d3
#define MESH_EDGES 7
#include "iso_one_point_minus_data_mesh.hpp"
#undef PHIM
#undef PHIP0
#undef PHIP1
#undef PHIP2
#undef XM
#undef XP0
#undef XP1
#undef XP2
#undef DM
#undef DP0
#undef DP1
#undef DP2
#undef MESH_EDGES
	}
	else {
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
#define DM0 d0
#define DM1 d3
#define DP0 d1
#define DP1 d2
#define MESH_EDGES0 5
#define MESH_EDGES1 5
#include "iso_two_point_minus_data_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
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
#define DM0 d0
#define DM1 d2
#define DP0 d3
#define DP1 d1
#define MESH_EDGES0 5
#define MESH_EDGES1 5
#include "iso_two_point_minus_data_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
	else {
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
#define DP d1
#define DM0 d0
#define DM1 d2
#define DM2 d3
#define MESH_EDGES 7
#include "iso_three_point_minus_data_mesh.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef DP
#undef DM0
#undef DM1
#undef DM2
#undef MESH_EDGES
	}
      }
    }
    else {
      // phi1 is minus...
      if (phi2 >= 0.0) {
	if (phi3 >= 0.0) {
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
#define DM0 d0
#define DM1 d1
#define DP0 d2
#define DP1 d3
#define MESH_EDGES0 5
#define MESH_EDGES1 5
#include "iso_two_point_minus_data_mesh.hpp"
#undef PHIM0
#undef PHIM1
#undef PHIP0
#undef PHIP1
#undef XM0
#undef XM1
#undef XP0
#undef XP1
#undef DM0
#undef DM1
#undef DP0
#undef DP1
#undef MESH_EDGES0
#undef MESH_EDGES1
	}
	else {
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
#define DP d2
#define DM0 d0
#define DM1 d3
#define DM2 d1
#define MESH_EDGES 7
#include "iso_three_point_minus_data_mesh.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef DP
#undef DM0
#undef DM1
#undef DM2
#undef MESH_EDGES
	}
      }
      else {
	// phi2 is minus...
	if (phi3 >= 0.0) {
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
#define DP d3
#define DM0 d0
#define DM1 d1
#define DM2 d2
#define MESH_EDGES 7
#include "iso_three_point_minus_data_mesh.hpp"
#undef PHIP
#undef PHIM0
#undef PHIM1
#undef PHIM2
#undef XP
#undef XM0
#undef XM1
#undef XM2
#undef DP
#undef DM0
#undef DM1
#undef DM2
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
