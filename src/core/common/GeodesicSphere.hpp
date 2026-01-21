#ifndef _GEODESIC_SPHERE_HPP_
#define _GEODESIC_SPHERE_HPP_

#include "Common.hpp"

namespace GeodesicSphere {
  
  // functions related to the icosahedral sphere...
  
  int getSphereNodeCount(const int n); // n is the number of segments each edge is divided into... 
  int getSphereTriCount(const int n); // n is the number of segments each edge is diviedd into... 
  void addSphere(double (*xsp)[3],int (*spost)[3],const double xc[3],const double r,const int n,const bool flip); // flip=false: fluid inside, flip=true: fluid outside  
  void addSphere(double (*xsp)[3],const int n); // same as above, but ject sphere points as unit normals about (0,0,0)
  
}

#endif
