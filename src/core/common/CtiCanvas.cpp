#include "CtiCanvas.hpp"

#define CALC_PIXEL_I(X) (((X)[0]-x0[0])*e0[0]+((X)[1]-x0[1])*e0[1]+((X)[2]-x0[2])*e0[2])
#define CALC_PIXEL_J(X) (((X)[0]-x0[0])*e1[0]+((X)[1]-x0[1])*e1[1]+((X)[2]-x0[2])*e1[2])
#define CALC_PIXEL_DEPTH(X) (((X)[0]-x0[0])*e2[0]+((X)[1]-x0[1])*e2[1]+((X)[2]-x0[2])*e2[2])
// define some convient accessors into the pixel buffers...

#define PBUF_SNORMAL0  pbuf_float[ip*12   ]
#define PBUF_SNORMAL1  pbuf_float[ip*12+ 1]
#define PBUF_SNORMAL2  pbuf_float[ip*12+ 2]
#define PBUF_SDEPTH    pbuf_float[ip*12+ 3]
#define PBUF_SDATA     pbuf_float[ip*12+ 4]
#define PBUF_MDEPTH    pbuf_float[ip*12+ 5]
#define PBUF_PNORMAL0  pbuf_float[ip*12+ 6]
#define PBUF_PNORMAL1  pbuf_float[ip*12+ 7]
#define PBUF_PNORMAL2  pbuf_float[ip*12+ 8]
#define PBUF_PDEPTH    pbuf_float[ip*12+ 9]
#define PBUF_PDATA     pbuf_float[ip*12+10]
#define PBUF_PAUX      pbuf_float[ip*12+11]

#define PBUF_SZONE   pbuf_uint2[ip*5  ]
#define PBUF_PZONE   pbuf_uint2[ip*5+1]
#define PBUF_BITS    pbuf_uint2[ip*5+2]
#define PBUF_I       pbuf_uint2[ip*5+3]
#define PBUF_J       pbuf_uint2[ip*5+4]

// define some bits used in pbd bits...

#define INTERNAL_DATA_BIT 1   // volume/iso/particle data
#define ISURFACE_BIT      2   // internal surface (volume/iso/particle w/o data)
#define INTERNAL_BITS     3   // if ANY internal bit is set
#define SURFACE_DATA_BIT  4
#define SURFACE_BIT       8
#define SURFACE_BITS      12 // (0100|1000) if ANY surface bits is set
#define MSURFACE_BIT      16  // masking surface
#define MASK_BIT          32  // if using the mdepth buffer for masking
#define AUX_BIT           64  // if using the paux buffer for reduction (e.g. euclidean distance)
#define MESH_BIT          128
#define DELAY_BIT         256 // will construct data at flush, uses pdata for interim data
#define EDGE_BIT          512

// TODO we want to generalize the blanking to use other Geoms as well; below is just the
// implementation to register a blanking plane
void CtiCanvas::addBlankPlane(const float xp[3],const float np[3]) {

  // convert xp,np to x,n in ijk space...
  float ijk_blank[3];
  float nijk_blank[3];

  ijk_blank[0] = CALC_PIXEL_I(xp);
  ijk_blank[1] = CALC_PIXEL_J(xp);
  ijk_blank[2] = CALC_PIXEL_DEPTH(xp);

  nijk_blank[0] = np[0]*e0[0] + np[1]*e0[1] + np[2]*e0[2];
  nijk_blank[1] = np[0]*e1[0] + np[1]*e1[1] + np[2]*e1[2];
  nijk_blank[2] = np[0]*e2[0] + np[1]*e2[1] + np[2]*e2[2];

  vBlankPlaneData.push_back( PlaneData<float>(ijk_blank,nijk_blank) );

  b_blank = true;

}

void CtiCanvas::setDataPlane(const float xp[3],const float np[3]) {

  assert(dataPlaneData == NULL); dataPlaneData = new PlaneData<float>();

  dataPlaneData->center[0] = CALC_PIXEL_I(xp);
  dataPlaneData->center[1] = CALC_PIXEL_J(xp);
  dataPlaneData->center[2] = CALC_PIXEL_DEPTH(xp);

  dataPlaneData->normal[0] = np[0]*e0[0] + np[1]*e0[1] + np[2]*e0[2];
  dataPlaneData->normal[1] = np[0]*e1[0] + np[1]*e1[1] + np[2]*e1[2];
  dataPlaneData->normal[2] = np[0]*e2[0] + np[1]*e2[1] + np[2]*e2[2];

}

void CtiCanvas::removeDataPlane() {
  assert(dataPlaneData);
  delete dataPlaneData; dataPlaneData = NULL;
}

void CtiCanvas::setProjPlane() {
  assert(projPlane == NULL); projPlane = new PlaneData<float>();
  const float xp[3] = {0.0,0.0,0.0};
  const float np[3] = {1.0,0.0,0.0};

  projPlane->center[0] = CALC_PIXEL_I(xp);
  projPlane->center[1] = CALC_PIXEL_J(xp);
  projPlane->center[2] = CALC_PIXEL_DEPTH(xp);

  COUT1("projPlane->center " << COUT_VEC(projPlane->center));

  projPlane->normal[0] = np[0]*e0[0] + np[1]*e0[1] + np[2]*e0[2];
  projPlane->normal[1] = np[0]*e1[0] + np[1]*e1[1] + np[2]*e1[2];
  projPlane->normal[2] = np[0]*e2[0] + np[1]*e2[1] + np[2]*e2[2];

  COUT1("projPlane->normal " << COUT_VEC(projPlane->normal));

  double normal_mag = MAG(projPlane->normal);
  FOR_I3 projPlane->normal[i] /= (float) normal_mag;

  COUT1("projPlane->normal " << COUT_VEC(projPlane->normal));
}

void CtiCanvas::removeProjPlane() {
  assert(projPlane);
  delete projPlane; projPlane = NULL;
}


void CtiCanvas::insertBlankPlane(const float xp[3],const float np[3]) {

  // convert xp,np to x,n in ijk space...
  float ijk_blank[3];
  float nijk_blank[3];

  ijk_blank[0] = CALC_PIXEL_I(xp);
  ijk_blank[1] = CALC_PIXEL_J(xp);
  ijk_blank[2] = CALC_PIXEL_DEPTH(xp);

  nijk_blank[0] = np[0]*e0[0] + np[1]*e0[1] + np[2]*e0[2];
  nijk_blank[1] = np[0]*e1[0] + np[1]*e1[1] + np[2]*e1[2];
  nijk_blank[2] = np[0]*e2[0] + np[1]*e2[1] + np[2]*e2[2];

  // ME: what is going on here: why the strange insertion? so you can use the first for something?
  // just use *.back() ?
  if (vBlankPlaneData.size()>0)
    vBlankPlaneData.insert(vBlankPlaneData.begin(), PlaneData<float>(ijk_blank,nijk_blank) );
  else
    vBlankPlaneData.push_back( PlaneData<float>(ijk_blank,nijk_blank) );

  b_blank = true;

}

void CtiCanvas::deleteFirstBlankPlane() {

  // possible to request removal of blanking when no blanking plane was added (e.g.,stitch 3d_mesh viz)... 
  if (!b_blank)
    return;
  //assert(b_blank);
  assert(vBlankPlaneData.size() > 0);

  vBlankPlaneData.erase(vBlankPlaneData.begin());

  if (vBlankPlaneData.empty())
    b_blank = false;
}

void CtiCanvas::getBlankPlanes(vector<PlaneData<float> >& blank_planes) const {
  for (int iblank=0, limit=vBlankPlaneData.size(); iblank < limit; ++iblank) {
    blank_planes.push_back(vBlankPlaneData[iblank]);
  }
}

void CtiCanvas::deleteBlankPlanes() {
  vBlankPlaneData.clear();
  b_blank = false;
}

void CtiCanvas::insertBlankPlanes(const vector<PlaneData<float> >& blank_planes) {
  for (int iblank=0, limit=blank_planes.size(); iblank < limit; ++iblank) {
    vBlankPlaneData.push_back(blank_planes[iblank]);
  }
  b_blank = true;
}

bool CtiCanvas::isBlanked(const float i,const float j,const float d) const {

  if (!b_blank) return false;
  // else return ((i-ijk_blank[0])*nijk_blank[0]+(j-ijk_blank[1])*nijk_blank[1]+(d-ijk_blank[2])*nijk_blank[2] >= 0.0);
  else {
    for (int iblank=0, limit=vBlankPlaneData.size(); iblank < limit; ++iblank) {
      if (isBlanked(i,j,d,iblank)) return true;
    }

    // no blank encountered if we get here
    return false;
  }

}


bool CtiCanvas::isBlankedSkipDataPlane(const float i,const float j,const float d,const bool on_data_plane) const {

  if (!b_blank) return false;
  // else return ((i-ijk_blank[0])*nijk_blank[0]+(j-ijk_blank[1])*nijk_blank[1]+(d-ijk_blank[2])*nijk_blank[2] >= 0.0);
  else {
    int iblank0 = 0;
    if ((dataPlaneData != NULL) && b_blank_data_plane && on_data_plane) iblank0 = 1;
    for (int iblank=iblank0, limit=vBlankPlaneData.size(); iblank < limit; ++iblank) {
      if (isBlanked(i,j,d,iblank)) return true;
    }

    // no blank encountered if we get here
    return false;
  }

}

bool CtiCanvas::isBlankedSkipDataPlaneOrClipped(const float i,const float j,const float d,const bool on_data_plane,const float c,const float * my_range) const {

  if (my_range) {
    // if clipping range is specified
    if ( (c < my_range[0]) || (c > my_range[1]) ) return true;
  }

  return isBlankedSkipDataPlane(i,j,d,on_data_plane);

}

bool CtiCanvas::isBlanked(const float i,const float j,const float d,const int iblank) const {
  const float ijk[3] = {i,j,d};

  float value = 0.f;
  FOR_K3 value += (ijk[k]-vBlankPlaneData[iblank].center[k])*vBlankPlaneData[iblank].normal[k];

  return (value >= 0.0);
}

bool CtiCanvas::isBehind(const float i,const float j,const float d) const {
  if (dataPlaneData == NULL) {
    return true;
  }
  else  {
    const float ijk[3] = {i,j,d};

    float value = 0.f;
    FOR_K3 value += (ijk[k]-dataPlaneData->center[k])*dataPlaneData->normal[k];

    return (value*dataPlaneData->normal[2] >= 0.0);
  }
}

void CtiCanvas::computeBarycentricConsts(float& j1mj2,float& j0mj2,float& i2mi1,float& i0mi2,float& one_o_det,const float i_s[3],const float j_s[3]) const {
  computeBarycentricConsts(j1mj2,j0mj2,i2mi1,i0mi2,one_o_det,i_s[0],i_s[1],i_s[2],j_s[0],j_s[1],j_s[2]);
}

void CtiCanvas::computeBarycentricConsts(float& j1mj2,float& j0mj2,float& i2mi1,float& i0mi2,float& one_o_det,const float i0,const float i1,const float i2,const float j0,const float j1,const float j2) const {
  j1mj2 = (j1-j2);
  j0mj2 = (j0-j2);
  i2mi1 = (i2-i1);
  i0mi2 = (i0-i2);
  one_o_det = 1.f / (j1mj2*i0mi2 + i2mi1*j0mj2);
}

void CtiCanvas::computeBarycentricCoords(float (&lamda)[3],const int i,const int j,const float i2,const float j2,const float j1mj2,const float j0mj2,const float i2mi1,const float i0mi2,const float one_o_det) const {
  computeBarycentricCoords(lamda[0],lamda[1],lamda[2],i,j,i2,j2,j1mj2,j0mj2,i2mi1,i0mi2,one_o_det);
}

void CtiCanvas::computeBarycentricCoords(float& l1,float& l2,float& l3,const int i,const int j,const float i2,const float j2,const float j1mj2,const float j0mj2,const float i2mi1,const float i0mi2,const float one_o_det) const {
  const float arg1 = float(i)-i2;
  const float arg2 = float(j)-j2;
  l1 = max(0.f,(j1mj2*arg1 + i2mi1*arg2)*one_o_det);
  l2 = max(0.f,(-j0mj2*arg1 + i0mi2*arg2)*one_o_det);
  l3 = 1.0f-l1-l2;

  if (l3 < 0.f) {
    const float sum = l1+l2;
    l1 /= sum;
    l2 /= sum;
    l3 = 0.f;
  }
}

void CtiCanvas::addSurfaceTri(const double xt0[3],const double xt1[3],const double xt2[3],const uint2 zone,const bool hide,const uchar mesh,const bool on_data_plane,const double * _normal,const double * data,const bool cell_flood,const double * range) {
  // convert simulation coordinates to image coordinates (i,j,depth)
  const float d0 = CALC_PIXEL_DEPTH(xt0);
  const float d1 = CALC_PIXEL_DEPTH(xt1);
  const float d2 = CALC_PIXEL_DEPTH(xt2);

  if (max(d0,max(d1,d2)) <= 0.f) return;

  const float i0 = CALC_PIXEL_I(xt0);
  const float i1 = CALC_PIXEL_I(xt1);
  const float i2 = CALC_PIXEL_I(xt2);

  if ( (max(i0,max(i1,i2)) < 0.f) || (min(i0,min(i1,i2)) > float(ni-1)) ) return;

  const float j0 = CALC_PIXEL_J(xt0);
  const float j1 = CALC_PIXEL_J(xt1);
  const float j2 = CALC_PIXEL_J(xt2);

  if ( (max(j0,max(j1,j2)) < 0.f) || (min(j0,min(j1,j2)) > float(nj-1)) ) return;

  // determine tri normal in image coordinates
  float normal[3];
  if (_normal == NULL) {
    // not provided, so generate from tri nodes
    normal[0] = (j1-j0)*(d2-d0)-(d1-d0)*(j2-j0);
    normal[1] = (d1-d0)*(i2-i0)-(i1-i0)*(d2-d0);
    normal[2] = (i1-i0)*(j2-j0)-(j1-j0)*(i2-i0);
  }
  else {
    // provided in simulation coordinates, so convert
    normal[0] = float(_normal[0]*e0[0] + _normal[1]*e0[1] + _normal[2]*e0[2]);
    normal[1] = float(_normal[0]*e1[0] + _normal[1]*e1[1] + _normal[2]*e1[2]);
    normal[2] = float(_normal[0]*e2[0] + _normal[1]*e2[1] + _normal[2]*e2[2]);
  }
  const double mag = sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );
  if (mag == 0.f)
    return;

  FOR_I3 normal[i] /= mag;

  // sort the node info by y...
  float i_s[3],j_s[3],depth_s[3],data_s[3];
  int orig_index[3] = {0,1,2};  // map from sorted index to original tri index
  if (j0 <= j1) {
    if (j1 <= j2) {
      // do nothing; maintains tri indexed order
    }
    else if (j0 <= j2) {
      orig_index[1]=2;
      orig_index[2]=1;
    }
    else {
      orig_index[0]=1;
      orig_index[1]=2;
      orig_index[2]=0;
    }
  }
  else if (j0 <= j2) {
    orig_index[0]=1;
    orig_index[1]=0;
  }
  else if (j1 <= j2) {
    orig_index[0]=2;
    orig_index[1]=0;
    orig_index[2]=1;
  }
  else {
    orig_index[0]=2;
    orig_index[2]=0;
  }

  // store in sorted order
  i_s[orig_index[0]] = i0;
  i_s[orig_index[1]] = i1;
  i_s[orig_index[2]] = i2;

  j_s[orig_index[0]] = j0;
  j_s[orig_index[1]] = j1;
  j_s[orig_index[2]] = j2;

  depth_s[orig_index[0]] = d0;
  depth_s[orig_index[1]] = d1;
  depth_s[orig_index[2]] = d2;

  // check sort was successful
  assert(j_s[0] <= j_s[1]);
  assert(j_s[1] <= j_s[2]);

  if (data != NULL) {
    if (cell_flood) {
      // single data value passed for all nodes of this tri
      FOR_I3 data_s[orig_index[i]] = float( data[0] );
    }
    else {
      // data at triangle nodes
      FOR_I3 data_s[orig_index[i]] = float( data[i] );
    }
  }

  const float pixel_eps = 1.0E-6f;
  const float grid_delta = 2.f;
  const float mesh_fraction = 0.5f;

  // the sorting has divided the tri into 2 vertical pieces
  // we process columns of j, and for each j determine the rows i to update

  // barycentric interpolation constants for this triangle
  float j1mj2,j0mj2,i2mi1,i0mi2,one_o_det;
  computeBarycentricConsts(j1mj2,j0mj2,i2mi1,i0mi2,one_o_det,i_s,j_s);

  float lamda_grid[3] = {-1.f,-1.f,-1.f};  // shouldn't be possible to be less than 0, so this default means we ignore using if not positive
  if (mesh) {
    // indicates that a mesh edge will be highlighted
    // we can determine if distance is within grid_delta using barycentric coords
    // we just need a map  relating lamda_i to pixel distance
    // we compute this by finding the distance in pixels to the opposite edge midpoint
    FOR_K3 {
      if (mesh & (1 << k)) {
        // bit k is the edge we want to draw, so determine the node
        // opposite it; this is the lamda_i index we need to check for
        // proximity to the edge
        const int node_opp_edge = orig_index[(k+2)%3];
        const int n_adj0 = orig_index[k];
        const int n_adj1 = orig_index[(k+1)%3];

        // determine minimum number of pixels from this node to opposite edge
        const float d_pixels = sqrt( pow((i_s[n_adj0]-i_s[node_opp_edge])*(j_s[n_adj1]-j_s[node_opp_edge])-(i_s[n_adj1]-i_s[node_opp_edge])*(j_s[n_adj0]-j_s[node_opp_edge]),2)/( pow(i_s[n_adj0]-i_s[n_adj1],2) + pow(j_s[n_adj0]-j_s[n_adj1],2) ));

        // lamda below which we are within grid_delta of the edge
        lamda_grid[node_opp_edge] = min(1.f,grid_delta/d_pixels);
      }
    }
  }

  if (j_s[0] <= j_s[2]) {
    const float one_o_j1mj0 = 1.0f/(j_s[1] - j_s[0]);
    const float one_o_j2mj0 = 1.0f/(-j0mj2);
    const float one_o_j2mj1 = 1.0f/(-j1mj2);

    const int j0_plus = max(0,(int)ceil(j_s[0]-pixel_eps));
    const int j2_minus = min(nj-1,(int)floor(j_s[2]+pixel_eps));
    for (int j = j0_plus; j <= j2_minus; ++j) {

      int i_f,i_l;  // min/max i values for this j
      {
        // first do the bottom edge (n0 -> n2)
        float lamda_edge[3];
        lamda_edge[2] = (float(j)-j_s[0]) * one_o_j2mj0;
        lamda_edge[1] = 0.f;
        lamda_edge[0] = 1.f - lamda_edge[2];
        const float i_bottom = DOT_PRODUCT(lamda_edge,i_s);

        // now do top edge, which depends on which half of tri we are processing
        if (float(j) <= j_s[1]) {
          lamda_edge[1] = (float(j)-j_s[0]) * one_o_j1mj0;
          lamda_edge[0] = 1.f - lamda_edge[1];
          lamda_edge[2] = 0.f;
        }
        else {
          lamda_edge[0] = 0.f;
          lamda_edge[2] = (float(j)-j_s[1]) * one_o_j2mj1;
          lamda_edge[1] = 1.f - lamda_edge[2];
        }
        const float i_top = DOT_PRODUCT(lamda_edge,i_s);

        // can be in either order, so sort appropriately
        i_f = max(0,(int)ceil(min(i_bottom,i_top)-pixel_eps));
        i_l = min(ni-1,(int)floor(max(i_bottom,i_top)+pixel_eps));
      }

      if (i_f <= i_l) {
        for (int i = i_f; i <= i_l; ++i) {
          // interpolate depth and data from barycentric coordinates based in (i,j)
          float lamda[3];
          computeBarycentricCoords(lamda,i,j,i_s[2],j_s[2],j1mj2,j0mj2,i2mi1,i0mi2,one_o_det);
          const float d = DOT_PRODUCT(depth_s,lamda);

          float c = 0.0;
          if (data != NULL) c = DOT_PRODUCT(data_s,lamda);

          if (d > 0.f) {
            const bool visible = !(isBlankedSkipDataPlane(float(i),float(j),d,on_data_plane)||hide);
            const bool behind  = isBehind(float(i),float(j),d);

            // normal to pass into draw
            // may be adjusted if edge is drawn on this pixel
            float scaled_normal[3] = {normal[0],normal[1],normal[2]};

            if (mesh) {
              float factor = 1.f;  // based on distance to edge it ranges from [0.5*mesh_frac + 0.5, 1.0]

              FOR_K3 {
                if (mesh & (1 << k)) {
                  // the k-th edge should be drawn
                  const uchar node_opp_edge = orig_index[(k+2)%3];
                  const float frac = lamda[node_opp_edge]/lamda_grid[node_opp_edge];
                  if ((frac > 0.f) && (frac < 1.f)) factor = min(factor,float((mesh_fraction*0.5f+0.5f) + (mesh_fraction*0.5f-0.5f)*cos(frac*M_PI)));
                }
              }

              FOR_K3 scaled_normal[k] *= factor;
            }

            int IJ,ij;
            getIJij(IJ,ij,i,j);
            float sdepth0,mdepth0;
            uint2 bits0;
            if (getPixelSdepthMdepthBits(sdepth0,mdepth0,bits0,IJ,ij)) {
              uint2 bits = bits0;
              if (visible) {
                if (bits0&SURFACE_BITS) {
                  if (d < sdepth0) {
                    // closer than original surface
                    if (data != NULL) {
                      bits &= ~SURFACE_BIT;     // clear surface bit
                      bits |= SURFACE_DATA_BIT; // set data bit
                      setPixelSnormalSdepthSdataSzoneBits(IJ,ij,scaled_normal,d,c,zone,bits);
                    }
                    else {
                      bits |= SURFACE_BIT;       // set surface bit
                      bits &= ~SURFACE_DATA_BIT; // clear data bit
                      setPixelSnormalSdepthSzoneBits(IJ,ij,scaled_normal,d,zone,bits);
                    }
                  }
                }
                else {
                  // originally background
                  if (data != NULL) {
                    bits |= SURFACE_DATA_BIT;
                    setPixelSnormalSdepthSdataSzoneBits(IJ,ij,scaled_normal,d,c,zone,bits);
                  }
                  else {
                    bits |= SURFACE_BIT;
                    setPixelSnormalSdepthSzoneBits(IJ,ij,scaled_normal,d,zone,bits);
                  }
                }
              }
              if (behind) {
                if (bits0&MSURFACE_BIT) {
                  if (d < fabs(mdepth0)) {
                    // closer than original masking surface
                    if (normal[2] >= 0.0f) setPixelMdepthBits(IJ,ij,-d,bits);
                    else setPixelMdepth(IJ,ij,d);
                  }
                }
                else {
                  if (normal[2] >= 0.0f) setPixelMdepthBits(IJ,ij,-d,bits|MSURFACE_BIT);
                  else setPixelMdepthBits(IJ,ij,d,bits|MSURFACE_BIT);
                }
              }
            }
            else {
              pbd[IJ] = new PixelBlockData();
              uint2 bits = 0;
              // treat like originally background
              if (visible) {
                if (data != NULL) {
                  bits |= SURFACE_DATA_BIT;
                  setPixelSnormalSdepthSdataSzoneBits(IJ,ij,scaled_normal,d,c,zone,bits);
                }
                else {
                  bits |= SURFACE_BIT;
                  setPixelSnormalSdepthSzoneBits(IJ,ij,scaled_normal,d,zone,bits);
                }
              }
              if (behind) {
                if (normal[2] >= 0.0f) setPixelMdepthBits(IJ,ij,-d,bits|MSURFACE_BIT);
                else setPixelMdepthBits(IJ,ij,d,bits|MSURFACE_BIT);
              }
            }
          }
        }
      }
    }
  }
}

void CtiCanvas::addSurfaceTriThickened(const double xt0[3],const double xt1[3],const double xt2[3],const uint2 zone,const bool hide,const uchar mesh,const bool on_data_plane,const double * _normal,const double * data,const bool cell_flood,const double * range) {

  // convert simulation coordinates to image coordinates (i,j,depth)
  const float d0 = CALC_PIXEL_DEPTH(xt0);
  const float d1 = CALC_PIXEL_DEPTH(xt1);
  const float d2 = CALC_PIXEL_DEPTH(xt2);

  if (max(d0,max(d1,d2)) <= 0.f) return;

  const float i0 = CALC_PIXEL_I(xt0);
  const float i1 = CALC_PIXEL_I(xt1);
  const float i2 = CALC_PIXEL_I(xt2);

  if ( (max(i0,max(i1,i2)) < 0.f) || (min(i0,min(i1,i2)) > float(ni-1)) ) return;

  const float j0 = CALC_PIXEL_J(xt0);
  const float j1 = CALC_PIXEL_J(xt1);
  const float j2 = CALC_PIXEL_J(xt2);

  if ( (max(j0,max(j1,j2)) < 0.f) || (min(j0,min(j1,j2)) > float(nj-1)) ) return;

  // determine tri normal in image coordinates
  float normal[3];
  if (_normal == NULL) {
    // not provided, so generate from tri nodes
    normal[0] = (j1-j0)*(d2-d0)-(d1-d0)*(j2-j0);
    normal[1] = (d1-d0)*(i2-i0)-(i1-i0)*(d2-d0);
    normal[2] = (i1-i0)*(j2-j0)-(j1-j0)*(i2-i0);
  }
  else {
    // provided in simulation coordinates, so convert
    normal[0] = float(_normal[0]*e0[0] + _normal[1]*e0[1] + _normal[2]*e0[2]);
    normal[1] = float(_normal[0]*e1[0] + _normal[1]*e1[1] + _normal[2]*e1[2]);
    normal[2] = float(_normal[0]*e2[0] + _normal[1]*e2[1] + _normal[2]*e2[2]);
  }
  const double mag = sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );
  if (mag == 0.f)
    return;

  FOR_I3 normal[i] /= mag;

  // sort the node info by y...
  float i_s[3],j_s[3],depth_s[3],data_s[3];
  int orig_index[3] = {0,1,2};  // map from sorted index to original tri index
  if (j0 <= j1) {
    if (j1 <= j2) {
      // do nothing; maintains tri indexed order
    }
    else if (j0 <= j2) {
      orig_index[1]=2;
      orig_index[2]=1;
    }
    else {
      orig_index[0]=1;
      orig_index[1]=2;
      orig_index[2]=0;
    }
  }
  else if (j0 <= j2) {
    orig_index[0]=1;
    orig_index[1]=0;
  }
  else if (j1 <= j2) {
    orig_index[0]=2;
    orig_index[1]=0;
    orig_index[2]=1;
  }
  else {
    orig_index[0]=2;
    orig_index[2]=0;
  }

  // store in sorted order
  i_s[orig_index[0]] = i0;
  i_s[orig_index[1]] = i1;
  i_s[orig_index[2]] = i2;

  j_s[orig_index[0]] = j0;
  j_s[orig_index[1]] = j1;
  j_s[orig_index[2]] = j2;

  depth_s[orig_index[0]] = d0;
  depth_s[orig_index[1]] = d1;
  depth_s[orig_index[2]] = d2;

  // check sort was successful
  assert(j_s[0] <= j_s[1]);
  assert(j_s[1] <= j_s[2]);

  if (data != NULL) {
    if (cell_flood) {
      // single data value passed for all nodes of this tri
      FOR_I3 data_s[orig_index[i]] = float( data[0] );
    }
    else {
      // data at triangle nodes
      FOR_I3 data_s[orig_index[i]] = float( data[i] );
    }
  }

  const float pixel_eps = 1.0E-6f;
  const float grid_delta = 2.f;
  const float mesh_fraction = 0.5f;

  // the sorting has divided the tri into 2 vertical pieces
  // we process columns of j, and for each j determine the rows i to update

  // barycentric interpolation constants for this triangle
  float j1mj2,j0mj2,i2mi1,i0mi2,one_o_det;
  computeBarycentricConsts(j1mj2,j0mj2,i2mi1,i0mi2,one_o_det,i_s,j_s);

  float lamda_grid[3] = {-1.f,-1.f,-1.f};  // shouldn't be possible to be less than 0, so this default means we ignore using if not positive
  if (mesh) {
    // indicates that a mesh edge will be highlighted
    // we can determine if distance is within grid_delta using barycentric coords
    // we just need a map  relating lamda_i to pixel distance
    // we compute this by finding the distance in pixels to the opposite edge midpoint
    FOR_K3 {
      if (mesh & (1 << k)) {
        // bit k is the edge we want to draw, so determine the node
        // opposite it; this is the lamda_i index we need to check for
        // proximity to the edge
        const int node_opp_edge = orig_index[(k+2)%3];
        const int n_adj0 = orig_index[k];
        const int n_adj1 = orig_index[(k+1)%3];

        // determine minimum number of pixels from this node to opposite edge
        const float d_pixels = sqrt( pow((i_s[n_adj0]-i_s[node_opp_edge])*(j_s[n_adj1]-j_s[node_opp_edge])-(i_s[n_adj1]-i_s[node_opp_edge])*(j_s[n_adj0]-j_s[node_opp_edge]),2)/( pow(i_s[n_adj0]-i_s[n_adj1],2) + pow(j_s[n_adj0]-j_s[n_adj1],2) ));

        // lamda below which we are within grid_delta of the edge
        lamda_grid[node_opp_edge] = min(1.f,grid_delta/d_pixels);
      }
    }
  }

  if (j_s[0] <= j_s[2]) {
    const float one_o_j1mj0 = 1.0f/(j_s[1] - j_s[0]);
    const float one_o_j2mj0 = 1.0f/(-j0mj2);
    const float one_o_j2mj1 = 1.0f/(-j1mj2);

    const int j0_plus = max(0,(int)ceil(j_s[0]-pixel_eps));
    const int j2_minus = min(nj-1,(int)floor(j_s[2]+pixel_eps));

    //Force single pixel resolution if j_s[0] and j_s[2] fall in same pixel
    int j0_plus_adj  = j0_plus;
    int j2_minus_adj = j2_minus;
    if (j0_plus==j2_minus+1&&j2_minus>=0){
      j0_plus_adj = j2_minus;
    }

    for (int j = j0_plus_adj; j <= j2_minus_adj; ++j) {

      int i_f,i_l;  // min/max i values for this j
      if (j0_plus>=j2_minus){
        i_f = max(0,(int)ceil(min(i_s[0],min(i_s[1],i_s[2]))-pixel_eps));
        i_l = min(ni-1,(int)floor(max(i_s[0],max(i_s[1],i_s[2]))+pixel_eps));
        if (i_f == i_l+1 && i_l>=0){//Force single pixel resolution if i_f[0] and i_l[2] fall in same pixel
          i_f = i_l;
        }
      }
      else{
        // first do the bottom edge (n0 -> n2)
        float lamda_edge[3];
        lamda_edge[2] = (float(j)-j_s[0]) * one_o_j2mj0;
        lamda_edge[1] = 0.f;
        lamda_edge[0] = 1.f - lamda_edge[2];
        const float i_bottom = DOT_PRODUCT(lamda_edge,i_s);

        // now do top edge, which depends on which half of tri we are processing
        if (float(j) <= j_s[1]) {
          lamda_edge[1] = (float(j)-j_s[0]) * one_o_j1mj0;
          lamda_edge[0] = 1.f - lamda_edge[1];
          lamda_edge[2] = 0.f;
        }
        else {
          lamda_edge[0] = 0.f;
          lamda_edge[2] = (float(j)-j_s[1]) * one_o_j2mj1;
          lamda_edge[1] = 1.f - lamda_edge[2];
        }
        const float i_top = DOT_PRODUCT(lamda_edge,i_s);

        // can be in either order, so sort appropriately
        i_f = max(0,(int)ceil(min(i_bottom,i_top)-pixel_eps));
        i_l = min(ni-1,(int)floor(max(i_bottom,i_top)+pixel_eps));
        if (i_f == i_l+1 && i_l>=0){ //Force single pixel resolution if i_f[0] and i_l[2] fall in same pixel
          i_f = i_l;
        }
      }

      if (i_f <= i_l) {
        for (int i = i_f; i <= i_l; ++i) {
          // interpolate depth and data from barycentric coordinates based in (i,j)
          float lamda[3];
          computeBarycentricCoords(lamda,i,j,i_s[2],j_s[2],j1mj2,j0mj2,i2mi1,i0mi2,one_o_det);
          const float d = DOT_PRODUCT(depth_s,lamda);

          float c = 0.0;
          if (data != NULL) c = DOT_PRODUCT(data_s,lamda);

          if (d > 0.f) {
            const bool visible = !(isBlankedSkipDataPlane(float(i),float(j),d,on_data_plane)||hide);
            const bool behind  = isBehind(float(i),float(j),d);

            // normal to pass into draw
            // may be adjusted if edge is drawn on this pixel
            float scaled_normal[3] = {normal[0],normal[1],normal[2]};

            if (mesh) {
              float factor = 1.f;  // based on distance to edge it ranges from [0.5*mesh_frac + 0.5, 1.0]

              FOR_K3 {
                if (mesh & (1 << k)) {
                  // the k-th edge should be drawn
                  const uchar node_opp_edge = orig_index[(k+2)%3];
                  const float frac = lamda[node_opp_edge]/lamda_grid[node_opp_edge];
                  if ((frac > 0.f) && (frac < 1.f)) factor = min(factor,float((mesh_fraction*0.5f+0.5f) + (mesh_fraction*0.5f-0.5f)*cos(frac*M_PI)));
                }
              }

              FOR_K3 scaled_normal[k] *= factor;
            }

            int IJ,ij;
            getIJij(IJ,ij,i,j);
            float sdepth0,mdepth0;
            uint2 bits0;
            if (getPixelSdepthMdepthBits(sdepth0,mdepth0,bits0,IJ,ij)) {
              uint2 bits = bits0;
              if (visible) {
                if (bits0&SURFACE_BITS) {
                  if (d < sdepth0) {
                    // closer than original surface
                    if (data != NULL) {
                      bits &= ~SURFACE_BIT;     // clear surface bit
                      bits |= SURFACE_DATA_BIT; // set data bit
                      setPixelSnormalSdepthSdataSzoneBits(IJ,ij,scaled_normal,d,c,zone,bits);
                    }
                    else {
                      bits |= SURFACE_BIT;       // set surface bit
                      bits &= ~SURFACE_DATA_BIT; // clear data bit
                      setPixelSnormalSdepthSzoneBits(IJ,ij,scaled_normal,d,zone,bits);
                    }
                  }
                }
                else {
                  // originally background
                  if (data != NULL) {
                    bits |= SURFACE_DATA_BIT;
                    setPixelSnormalSdepthSdataSzoneBits(IJ,ij,scaled_normal,d,c,zone,bits);
                  }
                  else {
                    bits |= SURFACE_BIT;
                    setPixelSnormalSdepthSzoneBits(IJ,ij,scaled_normal,d,zone,bits);
                  }
                }
              }
              if (behind) {
                if (bits0&MSURFACE_BIT) {
                  if (d < fabs(mdepth0)) {
                    // closer than original masking surface
                    if (normal[2] >= 0.0f) setPixelMdepthBits(IJ,ij,-d,bits);
                    else setPixelMdepth(IJ,ij,d);
                  }
                }
                else {
                  if (normal[2] >= 0.0f) setPixelMdepthBits(IJ,ij,-d,bits|MSURFACE_BIT);
                  else setPixelMdepthBits(IJ,ij,d,bits|MSURFACE_BIT);
                }
              }
            }
            else {
              pbd[IJ] = new PixelBlockData();
              uint2 bits = 0;
              // treat like originally background
              if (visible) {
                if (data != NULL) {
                  bits |= SURFACE_DATA_BIT;
                  setPixelSnormalSdepthSdataSzoneBits(IJ,ij,scaled_normal,d,c,zone,bits);
                }
                else {
                  bits |= SURFACE_BIT;
                  setPixelSnormalSdepthSzoneBits(IJ,ij,scaled_normal,d,zone,bits);
                }
              }
              if (behind) {
                if (normal[2] >= 0.0f) setPixelMdepthBits(IJ,ij,-d,bits|MSURFACE_BIT);
                else setPixelMdepthBits(IJ,ij,d,bits|MSURFACE_BIT);
              }
            }
          }
        }
      }
    }
  }
}



void CtiCanvas::addInternalTri(const double xt0[3],const double xt1[3],const double xt2[3],const uint2 zone,const uchar mesh,const bool on_data_plane,const double * data,const bool cell_flood,const float * range) {

  // if data range was specified, we can determine whether to fully skip this tri or not here
  if (data && range) {
    // conditions for all nodes being out of range
    if ( (data[0] < double(range[0])) && (data[1] < double(range[0])) && (data[2] < double(range[0])) ) return;

    if ( (data[0] > double(range[1])) && (data[1] > double(range[1])) && (data[2] > double(range[1])) ) return;
  }

  // convert simulation coordinates to image coordinates (i,j,depth)
  const float d0 = CALC_PIXEL_DEPTH(xt0);
  const float d1 = CALC_PIXEL_DEPTH(xt1);
  const float d2 = CALC_PIXEL_DEPTH(xt2);

  if (max(d0,max(d1,d2)) <= 0.f) return;

  const float i0 = CALC_PIXEL_I(xt0);
  const float i1 = CALC_PIXEL_I(xt1);
  const float i2 = CALC_PIXEL_I(xt2);

  if ( (max(i0,max(i1,i2)) < 0.f) || (min(i0,min(i1,i2)) > float(ni-1)) ) return;

  const float j0 = CALC_PIXEL_J(xt0);
  const float j1 = CALC_PIXEL_J(xt1);
  const float j2 = CALC_PIXEL_J(xt2);

  if ( (max(j0,max(j1,j2)) < 0.f) || (min(j0,min(j1,j2)) > float(nj-1)) ) return;

  // determine tri normal in image coordinates
  float normal[3];
  // not provided, so generate from tri nodes
  normal[0] = (j1-j0)*(d2-d0)-(d1-d0)*(j2-j0);
  normal[1] = (d1-d0)*(i2-i0)-(i1-i0)*(d2-d0);
  normal[2] = (i1-i0)*(j2-j0)-(j1-j0)*(i2-i0);

  const double mag = sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );
  if (mag == 0.f) return;

  FOR_I3 normal[i] /= mag;

  // sort the node info by y...
  float i_s[3],j_s[3],depth_s[3],data_s[3];
  int orig_index[3] = {0,1,2};  // map from sorted index to original tri index
  if (j0 <= j1) {
    if (j1 <= j2) {
      // do nothing; maintains tri indexed order
    }
    else if (j0 <= j2) {
      orig_index[1]=2;
      orig_index[2]=1;
    }
    else {
      orig_index[0]=1;
      orig_index[1]=2;
      orig_index[2]=0;
    }
  }
  else if (j0 <= j2) {
    orig_index[0]=1;
    orig_index[1]=0;
  }
  else if (j1 <= j2) {
    orig_index[0]=2;
    orig_index[1]=0;
    orig_index[2]=1;
  }
  else {
    orig_index[0]=2;
    orig_index[2]=0;
  }

  // store in sorted order
  i_s[orig_index[0]] = i0;
  i_s[orig_index[1]] = i1;
  i_s[orig_index[2]] = i2;

  j_s[orig_index[0]] = j0;
  j_s[orig_index[1]] = j1;
  j_s[orig_index[2]] = j2;

  depth_s[orig_index[0]] = d0;
  depth_s[orig_index[1]] = d1;
  depth_s[orig_index[2]] = d2;

  // check sort was successful
  assert(j_s[0] <= j_s[1]);
  assert(j_s[1] <= j_s[2]);

  if (data != NULL) {
    if (cell_flood) {
      // single data value passed for all nodes of this tri
      FOR_I3 data_s[orig_index[i]] = float( data[0] );
    }
    else {
      // data at triangle nodes
      FOR_I3 data_s[orig_index[i]] = float( data[i] );
    }
  }

  const float pixel_eps = 1.0E-6f;
  const float grid_delta = 2.f;
  const float mesh_fraction = 0.5f;

  // the sorting has divided the tri into 2 vertical pieces
  // we process columns of j, and for each j determine the rows i to update

  // barycentric interpolation constants for this triangle
  float j1mj2,j0mj2,i2mi1,i0mi2,one_o_det;
  computeBarycentricConsts(j1mj2,j0mj2,i2mi1,i0mi2,one_o_det,i_s,j_s);

  float lamda_grid[3] = {-1.f,-1.f,-1.f};  // shouldn't be possible to be less than 0, so this default means we ignore using if not positive
  if (mesh) {
    // indicates that a mesh edge will be highlighted
    // we can determine if distance is within grid_delta using barycentric coords
    // we just need a map  relating lamda_i to pixel distance
    // we compute this by finding the distance in pixels to the opposite edge midpoint
    FOR_K3 {
      if (mesh & (1 << k)) {
        // bit k is the edge we want to draw, so determine the node
        // opposite it; this is the lamda_i index we need to check for
        // proximity to the edge
        const int node_opp_edge = orig_index[(k+2)%3];
        const int n_adj0 = orig_index[k];
        const int n_adj1 = orig_index[(k+1)%3];

        // determine minimum number of pixels from this node to opposite edge
        const float d_pixels = sqrt( pow((i_s[n_adj0]-i_s[node_opp_edge])*(j_s[n_adj1]-j_s[node_opp_edge])-(i_s[n_adj1]-i_s[node_opp_edge])*(j_s[n_adj0]-j_s[node_opp_edge]),2)/( pow(i_s[n_adj0]-i_s[n_adj1],2) + pow(j_s[n_adj0]-j_s[n_adj1],2) ));

        // lamda below which we are within grid_delta of the edge
        lamda_grid[node_opp_edge] = min(1.f,grid_delta/d_pixels);
      }
    }
  }

  if (j_s[0] <= j_s[2]) {
    const float one_o_j1mj0 = 1.0f/(j_s[1] - j_s[0]);
    const float one_o_j2mj0 = 1.0f/(-j0mj2);
    const float one_o_j2mj1 = 1.0f/(-j1mj2);

    // range of possible j-values to render
    const int j0_plus = max(0,(int)ceil(j_s[0]-pixel_eps));
    const int j2_minus = min(nj-1,(int)floor(j_s[2]+pixel_eps));
    for (int j = j0_plus; j <= j2_minus; ++j) {
      // for each j we need to determine the i-range to render
      // we can use barycentric coords along the edges to determine this

      int i_f,i_l;  // min/max i values for this j
      {
        // first do the bottom edge (n0 -> n2)
        float lamda_edge[3];
        lamda_edge[2] = (float(j)-j_s[0]) * one_o_j2mj0;
        lamda_edge[1] = 0.f;
        lamda_edge[0] = 1.f - lamda_edge[2];
        const float i_bottom = DOT_PRODUCT(lamda_edge,i_s);

        // now do top edge, which depends on which half of tri we are processing
        if (float(j) <= j_s[1]) {
          lamda_edge[1] = (float(j)-j_s[0]) * one_o_j1mj0;
          lamda_edge[0] = 1.f - lamda_edge[1];
          lamda_edge[2] = 0.f;
        }
        else {
          lamda_edge[0] = 0.f;
          lamda_edge[2] = (float(j)-j_s[1]) * one_o_j2mj1;
          lamda_edge[1] = 1.f - lamda_edge[2];
        }
        const float i_top = DOT_PRODUCT(lamda_edge,i_s);

        // can be in either order, so sort appropriately
        i_f = max(0,(int)ceil(min(i_bottom,i_top)-pixel_eps));
        i_l = min(ni-1,(int)floor(max(i_bottom,i_top)+pixel_eps));
      }

      if (i_f <= i_l) {
        for (int i = i_f; i <= i_l; ++i) {
          // interpolate depth and data from barycentric coordinates based in (i,j)
          float lamda[3];
          computeBarycentricCoords(lamda,i,j,i_s[2],j_s[2],j1mj2,j0mj2,i2mi1,i0mi2,one_o_det);
          const float d = DOT_PRODUCT(depth_s,lamda);

          float c = 0.0;
          if (data != NULL) c = DOT_PRODUCT(data_s,lamda);

          const bool drawPixel = (d > 0.f)&&(!isBlankedSkipDataPlaneOrClipped(float(i),float(j),d,on_data_plane,c,range));  // only uses range if range != NULL

          if (drawPixel) {

            // normal to pass into draw
            // may be adjusted if edge is drawn on this pixel
            float scaled_normal[3] = {normal[0],normal[1],normal[2]};

            if (mesh) {
              float factor = 1.f;  // based on distance to edge it ranges from [0.5*mesh_frac + 0.5, 1.0]

              FOR_K3 {
                if (mesh & (1 << k)) {
                  // the k-th edge should be drawn
                  const uchar node_opp_edge = orig_index[(k+2)%3];
                  const float frac = lamda[node_opp_edge]/lamda_grid[node_opp_edge];
                  if ((frac > 0.f) && (frac < 1.f)) factor = min(factor,float((mesh_fraction*0.5f+0.5f) + (mesh_fraction*0.5f-0.5f)*cos(frac*M_PI)));
                }
              }

              FOR_K3 scaled_normal[k] *= factor;
            }

            int IJ,ij;
            getIJij(IJ,ij,i,j);
            float pdepth;
            uint2 bits;
            if (getPixelPdepthBits(pdepth,bits,IJ,ij)) {
              if (bits&INTERNAL_BITS) {
                if (d < pdepth) {
                  // closer than current internal data
                  if (data != NULL) {
                    bits |= INTERNAL_DATA_BIT;  // set data bit
                    bits &= ~ISURFACE_BIT; // clear surface bit
                    setPixelPnormalPdepthPdataPzoneBits(IJ,ij,scaled_normal,d,c,zone,bits);
                  }
                  else {
                    bits |= ISURFACE_BIT;       // set surface bit
                    bits &= ~INTERNAL_DATA_BIT; // clear data bit
                    setPixelPnormalPdepthPzoneBits(IJ,ij,scaled_normal,d,zone,bits);
                  }
                }
              }
              else {
                // first internal data
                if (data != NULL) setPixelPnormalPdepthPdataPzoneBits(IJ,ij,scaled_normal,d,c,zone,bits|INTERNAL_DATA_BIT);
                else setPixelPnormalPdepthPzoneBits(IJ,ij,scaled_normal,d,zone,bits|ISURFACE_BIT);
              }
            }
            else {
              pbd[IJ] = new PixelBlockData();
              // treat like first internal data
              if (data != NULL) setPixelPnormalPdepthPdataPzoneBits(IJ,ij,scaled_normal,d,c,zone,INTERNAL_DATA_BIT);
              else setPixelPnormalPdepthPzoneBits(IJ,ij,scaled_normal,d,zone,ISURFACE_BIT);
            }
          }
        }
      }
    }
  }
}

void CtiCanvas::addSurfaceMeshRvvPositive(const double (*x_vv)[3],const double * r_vv,const int *i_vv,const int *ivv_global,const int n,const double factor) {

  // use r_vv as a flag...
  for (int ivv = 0; ivv < n; ++ivv) {

    double i = CALC_PIXEL_I(x_vv[ivv]);
    double j = CALC_PIXEL_J(x_vv[ivv]);
    double k = CALC_PIXEL_DEPTH(x_vv[ivv]);

    // the radius of possible influence of the point is...
    assert(r_vv[ivv] > 0.0);
    const double rijk = r_vv[ivv]*factor*double(ni)/width; // add pixel eps

    assert(rijk > 0.0);

    for (int ip = max((int)floor(i-rijk),0), ilimit = min((int)ceil(i+rijk),ni-1); ip <= ilimit; ++ip) {
      for (int jp = max((int)floor(j-rijk),0), jlimit = min((int)ceil(j+rijk),nj-1); jp <= jlimit; ++jp) {
        // only consider pixels that can possibly touch a surface...
        const double d2_plane = (double(ip)-i)*(double(ip)-i)+(double(jp)-j)*(double(jp)-j);
        if (d2_plane < rijk*rijk) {
          // get the value, if any, at this pixel...
          int IJ,ij;
          getIJij(IJ,ij,ip,jp);
          float sdepth,paux;
          uint2 bits,szone;
          if (getPixelSdepthPauxSzoneBits(sdepth,paux,szone,bits,IJ,ij)) {
            if (bits&SURFACE_BITS) {
              if (i_vv[ivv] == szone) {
                // this pixel has either surface mask or surface data (already)...
                const double d2 = d2_plane + (k-double(sdepth))*(k-double(sdepth));
                assert(d2 == d2);
                if (d2 < rijk*rijk) {
                  if (bits & SURFACE_DATA_BIT) {
                    // assert that this routine was the one setting this surface data. We
                    // can check this by seeing if the aux bit is also set...
                    assert(bits&AUX_BIT);
                    assert(bits&MESH_BIT);
                    assert(bits&DELAY_BIT);
                    // there is already surface data set at this pixel, so only use our surface data
                    // if we are closer...
                    if (d2 < paux) {
                      setPixelSdataPaux(IJ,ij,ivv_global[ivv],d2);
                    }
                  }
                  else {
                    // this was just a surface without data...
                    bits &= ~SURFACE_BIT;
                    bits |= (SURFACE_DATA_BIT|AUX_BIT|MESH_BIT|DELAY_BIT);
                    setPixelSdataPauxBits(IJ,ij,ivv_global[ivv],d2,bits);
                  }
                }
              }
            }
          }
        }
      }
    }

  }

}

void CtiCanvas::addPlaneMeshRvvPositive(const double (*x_vv)[3],const double * r_vv,const int *i_vv,const int *ivv_global,const int n,const double factor) {

  double ijk[3];
  double nijk[3];
  FOR_I3 ijk[i] = dataPlaneData->center[i];
  FOR_I3 nijk[i] = dataPlaneData->normal[i];
  if (nijk[2] == 0.f) return;

  const double nijk_mag = MAG(nijk); assert(nijk_mag > 0.0);
  const double unit_nijk[3] = {nijk[0]/nijk_mag,nijk[1]/nijk_mag,nijk[2]/nijk_mag};
  const float normal[3] = {(float)unit_nijk[0],(float)unit_nijk[1],(float)unit_nijk[2]};

  // recall r_vv was flipped to negative when points do not need
  // to be considered...
  for (int ivv = 0; ivv < n; ++ivv) if (r_vv[ivv] > 0.0) {

      // compute formulas below with i,j,k in delta mode relative to the plane ijk.
      // this must be corrected when
      double i = CALC_PIXEL_I(x_vv[ivv]);
      double j = CALC_PIXEL_J(x_vv[ivv]);
      double k = CALC_PIXEL_DEPTH(x_vv[ivv]);

      // compute distance from the plane...

      const double dp = (i-ijk[0])*unit_nijk[0] + (j-ijk[1])*unit_nijk[1] + (k-ijk[2])*unit_nijk[2];
      // the radius of the point on the plane...
      const double r2 = r_vv[ivv]*r_vv[ivv]*factor*factor*double(ni*ni)/(width*width) - dp*dp;
      if (r2 > 0.0) {
        // move the i,j,k point onto the plane...
        i -= dp*unit_nijk[0];
        j -= dp*unit_nijk[1];
        k -= dp*unit_nijk[2];
        // now i,j,k is in the imaging data plane. Figure out the limits of the projected sphere r2 in ip,jp...
        // for these details see maple worksheet voronoi plane.mw...
        const double tmp1 = (nijk[0]*nijk[0]+nijk[2]*nijk[2]); assert(tmp1 > 0.0);
        const double djp = sqrt(r2*tmp1)/nijk_mag;
        const int jp_min = max((int)ceil(j-djp),0);
        const int jp_max = min((int)floor(j+djp),nj-1);
        for (int jp = jp_min; jp <= jp_max; ++jp) {
          const double tmp2 = sqrt(r2*tmp1 - (double(jp)-j)*(double(jp)-j)*nijk_mag*nijk_mag)*nijk[2]*(nijk[0] >= 0.0 ? 1.0 : -1.0);
          const double ip0 = -(tmp2 + (double(jp)-j)*nijk[0]*nijk[1])/tmp1;
          const double ip1 = (tmp2 - (double(jp)-j)*nijk[0]*nijk[1])/tmp1;
          const int ip_min = max((int)ceil(i+min(ip0,ip1)),0);
          const int ip_max = min((int)floor(i+max(ip0,ip1)),ni-1);
          for (int ip = ip_min; ip <= ip_max; ++ip) {
            const float depth = dataPlaneData->center[2] - (dataPlaneData->normal[0]*(ip-dataPlaneData->center[0]) +
                                                            dataPlaneData->normal[1]*(jp-dataPlaneData->center[1]))/dataPlaneData->normal[2];
            if (!isBlankedSkipDataPlane(ip,jp,depth,true)) {
              int IJ,ij;
              getIJij(IJ,ij,ip,jp);
              float paux;
              uint2 bits;
              if (getPixelPauxBits(paux,bits,IJ,ij)) {
                if (bits&INTERNAL_DATA_BIT) {
                  assert(bits&MESH_BIT);
                  assert(bits&MASK_BIT);
                  assert(bits&DELAY_BIT);
                  assert(bits&AUX_BIT); // should not have other volume data concurrently
                  const double dkp = ((double(ip)-i)*nijk[0]+(double(jp)-j)*nijk[1])/nijk[2];
                  const double my_d2 = dp*dp + (double(ip)-i)*(double(ip)-i) + (double(jp)-j)*(double(jp)-j) + dkp*dkp;
                  if (my_d2 < paux)
                    setPixelPnormalPdepthPdataPauxPzone(IJ,ij,normal,depth,(float)ivv_global[ivv],my_d2,i_vv[ivv]); // leave bits unchanged
                }
                else {
                  // use the paux to store the SQUARE of the Euclidean distance of the ORIGINAL point to the rendered pixel...
                  const double dkp = ((double(ip)-i)*nijk[0]+(double(jp)-j)*nijk[1])/nijk[2];
                  paux = dp*dp + (double(ip)-i)*(double(ip)-i) + (double(jp)-j)*(double(jp)-j) + dkp*dkp;
                  setPixelPnormalPdepthPdataPauxPzoneBits(IJ,ij,normal,depth,(float)ivv_global[ivv],paux,i_vv[ivv],bits|INTERNAL_DATA_BIT|MESH_BIT|AUX_BIT|MASK_BIT|DELAY_BIT);
                }
              }
              else {
                pbd[IJ] = new PixelBlockData();
                bits = 0;
                // use the paux to store the SQUARE of the Euclidean distance of the ORIGINAL point to the rendered pixel...
                const double dkp = ((double(ip)-i)*nijk[0]+(double(jp)-j)*nijk[1])/nijk[2];
                paux = dp*dp + (double(ip)-i)*(double(ip)-i) + (double(jp)-j)*(double(jp)-j) + dkp*dkp;
                setPixelPnormalPdepthPdataPauxPzoneBits(IJ,ij,normal,depth,(float)ivv_global[ivv],paux,i_vv[ivv],bits|INTERNAL_DATA_BIT|MESH_BIT|AUX_BIT|MASK_BIT|DELAY_BIT);
              }
            }
          }
        }
      }
    }

}


void CtiCanvas::addPointsMeshRvvPositive(const double (*x_vv)[3],const double * r_vv,const int *i_vv,const int n,const double factor) {

  double ijk[3];
  double nijk[3];
  FOR_I3 ijk[i] = dataPlaneData->center[i];
  FOR_I3 nijk[i] = dataPlaneData->normal[i];
  if (nijk[2] == 0.f) return;

  // use r_vv as a flag...
  for (int ivv = 0; ivv < n; ++ivv) if (r_vv[ivv] > 0.0) {

      const float rp = r_vv[ivv]*factor*float(ni)/width;
      if (rp == 0.f) continue;

      const float dc = CALC_PIXEL_DEPTH(x_vv[ivv]);
      if (dc <= 0.f) continue;

      const float ic = CALC_PIXEL_I(x_vv[ivv]);
      const float jc = CALC_PIXEL_J(x_vv[ivv]);

      const int i0 = max(0,(int)ceil(ic-rp));
      const int i1 = min(ni-1,(int)floor(ic+rp));
      const int j0 = max(0,(int)ceil(jc-rp));
      const int j1 = min(nj-1,(int)floor(jc+rp));
      for (int j = j0; j <= j1; ++j) {
        const float dj = float(j)-jc;
        for (int i = i0; i <= i1; ++i) {
          if (!isBlankedSkipDataPlane(i,j,ijk[2],true)) {
            const float di = float(i)-ic;
            const float dz2 = rp*rp - di*di - dj*dj;
            if (dz2 > 0.f) {
              // the normal...
              const float dz = sqrt(dz2);
              const float inv_rp = 1.f/rp;
              const float normal[3] = { di*inv_rp, dj*inv_rp, dz*inv_rp };
              const float depth = dc-dz;
              int IJ,ij;
              getIJij(IJ,ij,i,j);
              float pdepth;
              uint2 bits;
              if (getPixelPdepthBits(pdepth,bits,IJ,ij)) {
                if (bits&ISURFACE_BIT) {
                  assert(bits&MESH_BIT);
                  assert(bits&MASK_BIT);
                  if (depth < pdepth)
                    setPixelPnormalPdepthPzone(IJ,ij,normal,depth,i_vv[ivv]);
                }
                else {
                  setPixelPnormalPdepthPzoneBits(IJ,ij,normal,depth,i_vv[ivv],bits|ISURFACE_BIT|MESH_BIT|MASK_BIT);
                }
              }
              else {
                pbd[IJ] = new PixelBlockData();
                bits = 0;
                setPixelPnormalPdepthPzoneBits(IJ,ij,normal,depth,i_vv[ivv],bits|ISURFACE_BIT|MESH_BIT|MASK_BIT);
              }
            }
          }
        }
      }
    }
}

void CtiCanvas::addPlaneDataRvvPositive(const double * var,const double (*x_vv)[3],const double * r_vv,const int n,const double factor) {

  double ijk[3];
  double nijk[3];
  FOR_I3 ijk[i] = dataPlaneData->center[i];
  FOR_I3 nijk[i] = dataPlaneData->normal[i];
  if (nijk[2] == 0.0) return;

  const double nijk_mag = MAG(nijk); assert(nijk_mag > 0.0);
  const double unit_nijk[3] = {nijk[0]/nijk_mag,nijk[1]/nijk_mag,nijk[2]/nijk_mag};
  const float normal[3] = {(float)unit_nijk[0],(float)unit_nijk[1],(float)unit_nijk[2]};

  // use +ve r_vv as a flag (for example, to skip data that is "inactive" in moving solver)...
  for (int ivv = 0; ivv < n; ++ivv) if (r_vv[ivv] > 0.0) {

      // compute formulas below with i,j,k in delta mode relative to the plane ijk.
      // this must be corrected when
      double i = CALC_PIXEL_I(x_vv[ivv]);
      double j = CALC_PIXEL_J(x_vv[ivv]);
      double k = CALC_PIXEL_DEPTH(x_vv[ivv]);

      // compute distance from the plane...

      const double dp = (i-ijk[0])*unit_nijk[0] + (j-ijk[1])*unit_nijk[1] + (k-ijk[2])*unit_nijk[2];
      // the radius of the point on the plane...
      const double r2 = r_vv[ivv]*r_vv[ivv]*factor*factor*double(ni*ni)/(width*width) - dp*dp;
      if (r2 > 0.0) {
        // move the i,j,k point onto the plane...
        i -= dp*unit_nijk[0];
        j -= dp*unit_nijk[1];
        k -= dp*unit_nijk[2];
        // now i,j,k is in the imaging data plane. Figure out the limits of the projected sphere r2 in ip,jp...
        // for these details see maple worksheet voronoi plane.mw...
        const double tmp1 = (nijk[0]*nijk[0]+nijk[2]*nijk[2]); assert(tmp1 > 0.0);
        const double djp = sqrt(r2*tmp1)/nijk_mag;
        const int jp_min = max((int)ceil(j-djp),0);
        const int jp_max = min((int)floor(j+djp),nj-1);
        for (int jp = jp_min; jp <= jp_max; ++jp) {
          const double tmp2 = sqrt(r2*tmp1 - (double(jp)-j)*(double(jp)-j)*nijk_mag*nijk_mag)*nijk[2]*(nijk[0] >= 0.0 ? 1.0 : -1.0);
          const double ip0 = -(tmp2 + (double(jp)-j)*nijk[0]*nijk[1])/tmp1;
          const double ip1 = (tmp2 - (double(jp)-j)*nijk[0]*nijk[1])/tmp1;
          const int ip_min = max((int)ceil(i+min(ip0,ip1)),0);
          const int ip_max = min((int)floor(i+max(ip0,ip1)),ni-1);
          for (int ip = ip_min; ip <= ip_max; ++ip) {
            // TODO: consider converting to float...
            const float depth = dataPlaneData->center[2] - (dataPlaneData->normal[0]*(ip-dataPlaneData->center[0]) +
                                                            dataPlaneData->normal[1]*(jp-dataPlaneData->center[1]))/dataPlaneData->normal[2];
            if (!isBlankedSkipDataPlane(ip,jp,depth,true)) {
              int IJ,ij;
              getIJij(IJ,ij,ip,jp);
              float paux;
              uint2 pzone,bits;
              if (getPixelPauxPzoneBits(paux,pzone,bits,IJ,ij)) {
                if (bits&INTERNAL_DATA_BIT) {
                  assert(pzone == 65534);
                  assert(bits&MASK_BIT);
                  assert(bits&AUX_BIT); // should not have other volume data concurrently
                  const double dkp = ((double(ip)-i)*nijk[0]+(double(jp)-j)*nijk[1])/nijk[2];
                  const double my_d2 = dp*dp + (double(ip)-i)*(double(ip)-i) + (double(jp)-j)*(double(jp)-j) + dkp*dkp;
                  if (my_d2 < paux)
                    setPixelPnormalPdepthPdataPaux(IJ,ij,normal,depth,var[ivv],my_d2); // leave bits/zone unchanged
                }
                else {
                  // use the paux to store the SQUARE of the Euclidean distance of the ORIGINAL point to the rendered pixel...
                  const double dkp = ((double(ip)-i)*nijk[0]+(double(jp)-j)*nijk[1])/nijk[2];
                  paux = dp*dp + (double(ip)-i)*(double(ip)-i) + (double(jp)-j)*(double(jp)-j) + dkp*dkp;
                  // TODO throw unique cell qualifier into data and do reduction on flush
                  setPixelPnormalPdepthPdataPauxPzoneBits(IJ,ij,normal,depth,var[ivv],paux,65534,bits|INTERNAL_DATA_BIT|AUX_BIT|MASK_BIT);
                }
              }
              else {
                pbd[IJ] = new PixelBlockData();
                bits = 0;
                // use the paux to store the SQUARE of the Euclidean distance of the ORIGINAL point to the rendered pixel...
                const double dkp = ((double(ip)-i)*nijk[0]+(double(jp)-j)*nijk[1])/nijk[2];
                paux = dp*dp + (double(ip)-i)*(double(ip)-i) + (double(jp)-j)*(double(jp)-j) + dkp*dkp;
                setPixelPnormalPdepthPdataPauxPzoneBits(IJ,ij,normal,depth,var[ivv],paux,65534,bits|INTERNAL_DATA_BIT|AUX_BIT|MASK_BIT);
              }
            }
          }
        }
      }
    }

}

void CtiCanvas::addSurfaceData(const double * var,const double (*x_vv)[3],const double * r_vv,const int n,const double factor) {

  // use r_vv as a flag...
  for (int ivv = 0; ivv < n; ++ivv) {

    double i = CALC_PIXEL_I(x_vv[ivv]);
    double j = CALC_PIXEL_J(x_vv[ivv]);
    double k = CALC_PIXEL_DEPTH(x_vv[ivv]);

    // the radius of possible influence of the point is...
    assert(r_vv[ivv] > 0.0);
    const double rijk = r_vv[ivv]*factor*double(ni)/width; // add pixel eps

    assert(rijk > 0.0);

    for (int ip = max((int)floor(i-rijk),0), ilimit = min((int)ceil(i+rijk),ni-1); ip <= ilimit; ++ip) {
      for (int jp = max((int)floor(j-rijk),0), jlimit = min((int)ceil(j+rijk),nj-1); jp <= jlimit; ++jp) {
        // only consider pixels that can possibly touch a surface...
        const double d2_plane = (double(ip)-i)*(double(ip)-i)+(double(jp)-j)*(double(jp)-j);
        if (d2_plane < rijk*rijk) {
          // get the value, if any, at this pixel...
          int IJ,ij;
          getIJij(IJ,ij,ip,jp);
          float sdepth,paux;
          uint2 bits;
          if (getPixelSdepthPauxBits(sdepth,paux,bits,IJ,ij)) {
            if (bits&SURFACE_BITS) {
              // this pixel has either surface mask or surface data (already)...
              const double d2 = d2_plane + (k-double(sdepth))*(k-double(sdepth));
              assert(d2 == d2);
              if (d2 < rijk*rijk) {
                if (bits & SURFACE_DATA_BIT) {
                  // assert that this routine was the one setting this surface data. We
                  // can check this by seeing if the aux bit is also set...
                  assert(bits&AUX_BIT);
                  // there is already surface data set at this pixel, so only use our surface data
                  // if we are closer...
                  if (d2 < paux) {
                    setPixelSdataPaux(IJ,ij,var[ivv],d2);
                  }
                }
                else {
                  // this was just a surface without data...
                  bits &= ~SURFACE_BIT;
                  bits |= (SURFACE_DATA_BIT|AUX_BIT);
                  setPixelSdataPauxBits(IJ,ij,var[ivv],d2,bits);
                }
              }
            }
          }
        }
      }
    }

  }

}

void CtiCanvas::addPartialSurfaceData(const double * var,const double (*x_vv)[3],const double * r_vv,const int * zone_vv,const int n,const double factor) {

  // use r_vv as a flag...
  for (int ivv = 0; ivv < n; ++ivv) {

    double i = CALC_PIXEL_I(x_vv[ivv]);
    double j = CALC_PIXEL_J(x_vv[ivv]);
    double k = CALC_PIXEL_DEPTH(x_vv[ivv]);

    // the radius of possible influence of the point is...
    assert(r_vv[ivv] > 0.0);
    const double rijk = r_vv[ivv]*factor*double(ni)/width; // add pixel eps

    assert(rijk > 0.0);

    for (int ip = max((int)floor(i-rijk),0), ilimit = min((int)ceil(i+rijk),ni-1); ip <= ilimit; ++ip) {
      for (int jp = max((int)floor(j-rijk),0), jlimit = min((int)ceil(j+rijk),nj-1); jp <= jlimit; ++jp) {
        // only consider pixels that can possibly touch a surface...
        const double d2_plane = (double(ip)-i)*(double(ip)-i)+(double(jp)-j)*(double(jp)-j);
        if (d2_plane < rijk*rijk) {
          // get the value, if any, at this pixel...
          int IJ,ij;
          getIJij(IJ,ij,ip,jp);
          float sdepth,paux;
          uint2 bits,szone;
          if (getPixelSdepthPauxSzoneBits(sdepth,paux,szone,bits,IJ,ij)) {
            if (bits&SURFACE_BITS) {
              if (zone_vv[ivv] == szone) {
                // this pixel has either surface mask or surface data (already)...
                const double d2 = d2_plane + (k-double(sdepth))*(k-double(sdepth));
                assert(d2 == d2);
                if (d2 < rijk*rijk) {
                  if (bits & SURFACE_DATA_BIT) {
                    // assert that this routine was the one setting this surface data. We
                    // can check this by seeing if the aux bit is also set...
                    assert(bits&AUX_BIT);
                    // there is already surface data set at this pixel, so only use our surface data
                    // if we are closer...
                    if (d2 < paux) {
                      setPixelSdataPaux(IJ,ij,var[ivv],d2);
                    }
                  }
                  else {
                    // this was just a surface without data...
                    bits &= ~SURFACE_BIT;
                    bits |= (SURFACE_DATA_BIT|AUX_BIT);
                    setPixelSdataPauxBits(IJ,ij,var[ivv],d2,bits);
                  }
                }
              }
            }
          }
        }
      }
    }

  }

}

void CtiCanvas::addParticle(const double xp[3],const double dp,const bool b_skip_blanking) {

  const float rp = 0.5f*dp*float(ni)/width;
  if (rp == 0.f) return;

  const float dc = CALC_PIXEL_DEPTH(xp);
  if (dc <= 0.f) return;

  const float ic = CALC_PIXEL_I(xp);
  if ( (ic+rp < 0.f) || (ic-rp > float(ni-1)) ) return;

  const float jc = CALC_PIXEL_J(xp);
  if ( (jc+rp < 0.f) || (jc-rp > float(nj-1)) ) return;

  if ((dc - rp < 0.f)||((!b_skip_blanking)&&isBlanked(float(ic),float(jc),dc-rp)))
    return;

  const int i0 = max(0,(int)ceil(ic-rp));
  const int i1 = min(ni-1,(int)floor(ic+rp));
  const int j0 = max(0,(int)ceil(jc-rp));
  const int j1 = min(nj-1,(int)floor(jc+rp));
  for (int j = j0; j <= j1; ++j) {
    const float dj = float(j)-jc;
    for (int i = i0; i <= i1; ++i) {
      // the pbd check seemed to miss data in parrallel...
      const float di = float(i)-ic;
      const float dz2 = rp*rp - di*di - dj*dj;
      if (dz2 > 0.f) {
        const float dz = sqrt(dz2);
        const float inv_rp = 1.f/rp;
        const float normal[3] = { di*inv_rp, dj*inv_rp, dz*inv_rp };
        const float depth = dc-dz;
        int IJ,ij;
        getIJij(IJ,ij,i,j);
        float pdepth;
        uint2 bits;
        if (getPixelPdepthBits(pdepth,bits,IJ,ij)) {
          if (bits&INTERNAL_BITS) {
            if (depth < pdepth) {
              // closer than current internal data
              bits &= ~INTERNAL_DATA_BIT;  // clear data bit
              bits |= ISURFACE_BIT; // set surface bit
              setPixelPnormalPdepthPzoneBits(IJ,ij,normal,depth,32767,bits);
            }
          }
          else {
            // first internal data
            setPixelPnormalPdepthPzoneBits(IJ,ij,normal,depth,32767,bits|ISURFACE_BIT);
          }
        }
        else {
          pbd[IJ] = new PixelBlockData();
          // treat like first internal data
          setPixelPnormalPdepthPzoneBits(IJ,ij,normal,depth,32767,ISURFACE_BIT);
        }
      }
    }
  }
}

void CtiCanvas::addParticle(const double xp[3],const double dp, const double vp,const bool b_skip_blanking) {

  const float rp = 0.5f*dp*float(ni)/width;
  if (rp == 0.f) return;

  const float dc = CALC_PIXEL_DEPTH(xp);
  if (dc <= 0.f) return;

  const float ic = CALC_PIXEL_I(xp);
  if ( (ic+rp < 0.f) || (ic-rp > float(ni-1)) ) return;

  const float jc = CALC_PIXEL_J(xp);
  if ( (jc+rp < 0.f) || (jc-rp > float(nj-1)) ) return;

  if ((dc - rp < 0.f)||((!b_skip_blanking)&&isBlanked(float(ic),float(jc),dc-rp)))
    return;

  const int i0 = max(0,(int)ceil(ic-rp));
  const int i1 = min(ni-1,(int)floor(ic+rp));
  const int j0 = max(0,(int)ceil(jc-rp));
  const int j1 = min(nj-1,(int)floor(jc+rp));
  for (int j = j0; j <= j1; ++j) {
    const float dj = float(j)-jc;
    for (int i = i0; i <= i1; ++i) {
      const float di = float(i)-ic;
      const float dz2 = rp*rp - di*di - dj*dj;
      if (dz2 > 0.f) {
        const float dz = sqrt(dz2);
        const float inv_rp = 1.f/rp;
        const float normal[3] = { di*inv_rp, dj*inv_rp, dz*inv_rp };
        const float depth = dc-dz;
        int IJ,ij;
        getIJij(IJ,ij,i,j);
        float pdepth;
        uint2 bits;
        if (getPixelPdepthBits(pdepth,bits,IJ,ij)) {
          if (bits&INTERNAL_BITS) {
            if (depth < pdepth) {
              // closer than current internal data
              bits |= INTERNAL_DATA_BIT;  // set data bit
              bits &= ~ISURFACE_BIT; // clear surface bit
              setPixelPnormalPdepthPdataPzoneBits(IJ,ij,normal,depth,vp,65533,bits);
            }
          }
          else {
            // first internal data
            setPixelPnormalPdepthPdataPzoneBits(IJ,ij,normal,depth,vp,65533,bits|INTERNAL_DATA_BIT);
          }
        }
        else {
          pbd[IJ] = new PixelBlockData();
          // treat like first internal data
          setPixelPnormalPdepthPdataPzoneBits(IJ,ij,normal,depth,vp,65533,INTERNAL_DATA_BIT);
        }
      }
    }
  }
}

typedef struct {
  float intensity;
  float depth;
} IntensityDepth;

void CtiCanvas::writeImage(const string& filename, const RGB_MODE rgb_mode) {
  double wtime0 = MPI_Wtime();
  cacheImage();
  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0) {
    double wtime = MPI_Wtime();
    cout << " cache time: " << wtime-wtime0 << " [s] " << endl;
    wtime0 = wtime;
  }
  flushImage(filename,rgb_mode);
  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0) {
    double wtime = MPI_Wtime();
    cout << " flush time: " << wtime-wtime0 << " [s] " << endl;
    wtime0 = wtime;
  }
  clearPBuf();
}

void CtiCanvas::rmPbdsWithoutSurfaceData() {

  for (int J = 0; J < nJ; ++J) {
    for (int I = 0; I < nI; ++I) {
      const int IJ = J*nI+I;
      bool b_hasSurfaceData = false;
      if (pbd[IJ]) {
        for (int j = 0; j < 8; ++j) {
          for (int i = 0; i < 8; ++i) {
            const int ij = j*8+i;
            if (pbd[IJ]->bits[ij]&SURFACE_DATA_BIT) {
              b_hasSurfaceData = true;
              break;
            }
          }
          if (b_hasSurfaceData)
            break;
        }
        if (!b_hasSurfaceData) {
          delete pbd[IJ];
          pbd[IJ] = NULL;
        }
      }
    }
  }

}

void CtiCanvas::cacheImage() {

  // reduce/collect main image...

  assert(pbd);

  // count first...

  int my_pbuf_size = 0;
  for (int J = 0; J < nJ; ++J) {
    for (int I = 0; I < nI; ++I) {
      const int IJ = J*nI+I;
      if (pbd[IJ]) {
        for (int j = 0; j < 8; ++j) {
          for (int i = 0; i < 8; ++i) {
            const int ij = j*8+i;
            if (pbd[IJ]->bits[ij])
              ++my_pbuf_size;
          }
        }
      }
    }
  }

  // pack and reduce the pixel buffers...
  // recall the data in the pixel block...

  float * my_pbuf_float = new float[my_pbuf_size*12]; // snormal[3],sdepth,sdata,mdepth,pnormal[3],pdepth,pdata,paux
  uint2 * my_pbuf_uint2 = new uint2[my_pbuf_size*5];  // szone,pzone,bits,i,j

  my_pbuf_size = 0;
  for (int J = 0; J < nJ; ++J) {
    for (int I = 0; I < nI; ++I) {
      const int IJ = J*nI+I;
      if (pbd[IJ]) {
        for (int j = 0; j < 8; ++j) {
          for (int i = 0; i < 8; ++i) {
            const int ij = j*8+i;
            if (pbd[IJ]->bits[ij]) {

              // HACK -- skip SURFACE_BIT only...
              //if (pbd[IJ]->bits[ij]&SURFACE_BIT)
              //  continue;

              // float data...
              my_pbuf_float[my_pbuf_size*12   ] = pbd[IJ]->snormal[ij][0];
              my_pbuf_float[my_pbuf_size*12+ 1] = pbd[IJ]->snormal[ij][1];
              my_pbuf_float[my_pbuf_size*12+ 2] = pbd[IJ]->snormal[ij][2];
              my_pbuf_float[my_pbuf_size*12+ 3] = pbd[IJ]->sdepth[ij];
              my_pbuf_float[my_pbuf_size*12+ 4] = pbd[IJ]->sdata[ij];
              my_pbuf_float[my_pbuf_size*12+ 5] = pbd[IJ]->mdepth[ij];
              my_pbuf_float[my_pbuf_size*12+ 6] = pbd[IJ]->pnormal[ij][0];
              my_pbuf_float[my_pbuf_size*12+ 7] = pbd[IJ]->pnormal[ij][1];
              my_pbuf_float[my_pbuf_size*12+ 8] = pbd[IJ]->pnormal[ij][2];
              my_pbuf_float[my_pbuf_size*12+ 9] = pbd[IJ]->pdepth[ij];
              my_pbuf_float[my_pbuf_size*12+10] = pbd[IJ]->pdata[ij];
              my_pbuf_float[my_pbuf_size*12+11] = pbd[IJ]->paux[ij];
              // uint2 data...
              my_pbuf_uint2[my_pbuf_size*5  ] = pbd[IJ]->szone[ij];
              my_pbuf_uint2[my_pbuf_size*5+1] = pbd[IJ]->pzone[ij];
              my_pbuf_uint2[my_pbuf_size*5+2] = pbd[IJ]->bits[ij];
              my_pbuf_uint2[my_pbuf_size*5+3] = (I<<3)+i;
              my_pbuf_uint2[my_pbuf_size*5+4] = (J<<3)+j;
              // increment...
              ++my_pbuf_size;
            }
          }
        }
        delete pbd[IJ];
      }
    }
  }

  delete[] pbd; pbd = NULL;

  // now gather at the write process (not necessarily rank==0, but set this for now)...

  //int write_rank = 0;

  clearPBuf(); // TODO: why is this called here? why not the asserts below?

  //assert (pbuf_float==NULL);
  //assert (pbuf_uint2==NULL);
  //assert (pbuf_size==0);

  if (mpi_size == 1) {

    // in serial mode, just take copies of the buffers...
    assert(write_rank == 0);
    pbuf_size = my_pbuf_size;
    pbuf_float = my_pbuf_float;
    pbuf_uint2 = my_pbuf_uint2;

  }
  else {

    // in parallel, just gather everyone at write_rank...
    assert((write_rank >= 0)&&(write_rank < mpi_size));

    int * pixel_count = NULL;
    if (mpi_rank == write_rank) pixel_count = new int[mpi_size];
    MPI_Gather(&my_pbuf_size,1,MPI_INT,pixel_count,1,MPI_INT,write_rank,mpi_comm);

    int * pixel_disp = NULL;
    if (mpi_rank == write_rank) {
      pixel_disp = new int[mpi_size];
      pixel_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        pixel_disp[rank] = pixel_disp[rank-1] + pixel_count[rank-1];
      pbuf_size = pixel_disp[mpi_size-1] + pixel_count[mpi_size-1];
      pbuf_float = new float[pbuf_size*12];
      //for (int ip = 0; ip < pbuf_size*12; ++ip)
      //  pbuf_float[ip] = HUGE_VALF;
      pbuf_uint2 = new uint2[pbuf_size*5];
      //for (int ip = 0; ip < pbuf_size*5; ++ip)
      //  pbuf_uint2[ip] = 65535;
      // adjust count and disp for the floats...
      FOR_RANK {
        //cout << pixel_count[rank] << endl;
        pixel_count[rank] *= 12;
        pixel_disp[rank] *= 12;
      }
    }

    MPI_Gatherv(my_pbuf_float,my_pbuf_size*12,MPI_FLOAT,pbuf_float,pixel_count,pixel_disp,MPI_FLOAT,write_rank,mpi_comm);
    delete[] my_pbuf_float;

    if (mpi_rank == write_rank) {
      FOR_RANK {
        pixel_count[rank] = (pixel_count[rank]/12)*5;
        pixel_disp[rank]  = (pixel_disp[rank] /12)*5;
      }
    }

    MPI_Gatherv(my_pbuf_uint2,my_pbuf_size*5,MPI_UINT2,pbuf_uint2,pixel_count,pixel_disp,MPI_UINT2,write_rank,mpi_comm);
    delete[] my_pbuf_uint2;

    if (mpi_rank == write_rank) {
      //for (int ip = 0; ip < pbuf_size*12; ++ip)
      //  assert(pbuf_float[ip] < HUGE_VALF); // check
      delete[] pixel_count;
      delete[] pixel_disp;
    }

  }

  // =======================================================
  // =======================================================
  // =======================================================
  // if there is some volvis, we do that here too...
  // =======================================================
  // =======================================================
  // =======================================================

  if (b_volvis) {

    // we are going to round-robin the data associated with each (i,j), but first we need to
    // have all ij's of the same value on the same rank. Use this as a chance to load balance as well...

    if (mpi_size > 1) {

      int * send_count = new int[mpi_size];
      FOR_RANK send_count[rank] = 0;

      for (int ii = 0,ii_max=volVisDataVec.size(); ii < ii_max; ++ii) {
        const int rank = volVisDataVec[ii].ij%mpi_size;
        send_count[rank] += 1; // 2 ints and 3 floats eventually
      }

      int * send_disp = new int[mpi_size];
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
      const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

      int * send_buf_int     = new int[send_count_sum*2]; // 2 ints, 3 floats
      float * send_buf_float = new float[send_count_sum*3];

      for (int ii = 0,ii_max=volVisDataVec.size(); ii < ii_max; ++ii) {
        const int rank = volVisDataVec[ii].ij%mpi_size;
        // 2 ints...
        send_buf_int[send_disp[rank]*2  ] = volVisDataVec[ii].ij;
        send_buf_int[send_disp[rank]*2+1] = volVisDataVec[ii].flag;
        // 3 floats...
        send_buf_float[send_disp[rank]*3  ] = volVisDataVec[ii].depth;
        send_buf_float[send_disp[rank]*3+1] = volVisDataVec[ii].v0;
        send_buf_float[send_disp[rank]*3+2] = volVisDataVec[ii].v1;
        // and inc...
        ++send_disp[rank];
      }

      // rewind send_disp...
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      // setup recv side stuff...
      int * recv_count = new int[mpi_size];
      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      int * recv_disp = new int[mpi_size];
      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      FOR_RANK {
        send_count[rank] *= 2;
        send_disp[rank] *= 2;
        recv_count[rank] *= 2;
        recv_disp[rank] *= 2;
      }

      int * recv_buf_int = new int[recv_count_sum*2];
      MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
      delete[] send_buf_int;

      FOR_RANK {
        send_count[rank] = (send_count[rank]/2)*3;
        send_disp[rank] = (send_disp[rank]/2)*3;
        recv_count[rank] = (recv_count[rank]/2)*3;
        recv_disp[rank] = (recv_disp[rank]/2)*3;
      }

      float * recv_buf_float = new float[recv_count_sum*3];
      MPI_Alltoallv(send_buf_float,send_count,send_disp,MPI_FLOAT,recv_buf_float,recv_count,recv_disp,MPI_FLOAT,mpi_comm);
      delete[] send_buf_float;

      // now set he data back into the volVisDataVec. This allows us to sort...
      volVisDataVec.resize(recv_count_sum);
      for (int irecv= 0; irecv < recv_count_sum; ++irecv) {
        // 2 ints...
        volVisDataVec[irecv].ij    = recv_buf_int[irecv*2  ];
        volVisDataVec[irecv].flag  = recv_buf_int[irecv*2+1];
        // 3 floats...
        volVisDataVec[irecv].depth = recv_buf_float[irecv*3  ];
        volVisDataVec[irecv].v0    = recv_buf_float[irecv*3+1];
        volVisDataVec[irecv].v1    = recv_buf_float[irecv*3+2];
      }
      delete[] recv_buf_int;
      delete[] recv_buf_float;

      delete[] recv_count;
      delete[] recv_disp;
      delete[] send_count;
      delete[] send_disp;

    }

    // need to compute if not provided by user...
    if (!b_volvis_aux_data) {
      const int n = volVisDataVec.size();
      if (volvis_type == "X-RAY") {
        volvis_aux_data[0] = 0.0;
        volvis_aux_data[1] = 255.0;
      }
      else if (volvis_type == "SURFACE") {
        volvis_aux_data[0] = 0.5;  // this is used as alpha channel; default is opaque
        volvis_aux_data[1] = 255.0;
      }
      else if (volvis_type == "LINEAR") {
        // get min and max
        float my_buf[2] = {HUGE_VALF,HUGE_VALF};
        for (int i = 0; i < n; ++i) {
          if (volVisDataVec[i].flag <= 2) {
            if ((volVisDataVec[i].flag == -1)||(volVisDataVec[i].flag == 0)) { // TODO: no more -1!
              // only meaningful data in v0...
              my_buf[0] = min(my_buf[0],volVisDataVec[i].v0);
              my_buf[1] = min(my_buf[1],-volVisDataVec[i].v0);
            }
            else if (volVisDataVec[i].flag == 1) {
              my_buf[0] = min(my_buf[0],volVisDataVec[i].v1);
              my_buf[1] = min(my_buf[1],-volVisDataVec[i].v1);
            }
            else if (volVisDataVec[i].flag == 2) {
              my_buf[0] = min(my_buf[0],volVisDataVec[i].v0);
              my_buf[1] = min(my_buf[1],-volVisDataVec[i].v0);
              my_buf[0] = min(my_buf[0],volVisDataVec[i].v1);
              my_buf[1] = min(my_buf[1],-volVisDataVec[i].v1);
            }
            else {
              assert(0);
            }
          }
        }
        float buf[2];
        MPI_Allreduce(my_buf,buf,2,MPI_FLOAT,MPI_MIN,mpi_comm);
        volvis_aux_data[0] = buf[0];
        volvis_aux_data[1] = -buf[1];
      }
      else {
        assert(volvis_type == "GAUSSIAN");
        float my_buf[3] = {0.f,0.f,0.f}; // avg,std,cnt
        for (int i = 0; i < n; ++i) {
          if (volVisDataVec[i].flag <= 2) {
            if (volVisDataVec[i].flag == 0) {
              my_buf[0] += volVisDataVec[i].v0;
              my_buf[1] += volVisDataVec[i].v0*volVisDataVec[i].v0;
              my_buf[2] += 1.f;
            }
            else if (volVisDataVec[i].flag == 1) {
              my_buf[0] += volVisDataVec[i].v1;
              my_buf[1] += volVisDataVec[i].v1*volVisDataVec[i].v1;
              my_buf[2] += 1.f;
            }
            else if (volVisDataVec[i].flag == 2) {
              my_buf[0] += volVisDataVec[i].v0;
              my_buf[1] += volVisDataVec[i].v0*volVisDataVec[i].v0;
              my_buf[0] += volVisDataVec[i].v1;
              my_buf[1] += volVisDataVec[i].v1*volVisDataVec[i].v1;
              my_buf[2] += 2.f;
            }
            else {
              // figure this out later!...
              assert(0);
            }
          }
        }
        float buf[3];
        MPI_Allreduce(my_buf,buf,3,MPI_FLOAT,MPI_SUM,mpi_comm);
        volvis_aux_data[0] = buf[0]/buf[2];
        volvis_aux_data[1] = sqrt(buf[1]/buf[2]-volvis_aux_data[0]*volvis_aux_data[0]);
      }
    }
    else {
      assert(b_volvis_aux_data);
      if (volvis_type == "SURFACE") {
        // ensure alpha is [0,1]
        if (volvis_aux_data[0] < 0.0) volvis_aux_data[0] = 0.0;
        else if (volvis_aux_data[0] > 1.0) volvis_aux_data[0] = 1.0;
      }
    }

    // we should now have a set of complete ij's on each rank, with
    // the data associated with them approximately load balanced...

    sort(volVisDataVec.begin(),volVisDataVec.end());

    // put the result of the local pixel line-of-sight integration
    // for each complete ij into these vecs. We use vecs so we can
    // gather in mpi without repacking...
    vector<int> my_ij_vec;
    vector<IntensityDepth> my_data_vec;

    // now the vector should be sorted in terms of line-of-sight rays...
    vector<VolVisData>::iterator iter_end = volVisDataVec.begin();
    while (iter_end != volVisDataVec.end()) {
      vector<VolVisData>::iterator iter_begin = iter_end;
      do {
	++iter_end;
      }
      while ((iter_end != volVisDataVec.end())&&(iter_begin->ij == iter_end->ij));
      // we should now have a range of iterators associated with a particular ij...
      my_ij_vec.push_back(iter_begin->ij);
      my_data_vec.push_back(IntensityDepth());
      calcVolvisIntensityAndDepth(my_data_vec.back().intensity,my_data_vec.back().depth,iter_begin,iter_end);
    }
    assert(my_ij_vec.size() == my_data_vec.size());
    // integration is done -- get rid of the VolVisData...

    volVisDataVec.clear();

    // reduce this to a the volvis_* buffers on the write_rank...

    if (mpi_size == 1) {

      // must be us...
      assert(write_rank == 0);

      // in serial mode, just take copies of the buffers...
      volvis_size = my_ij_vec.size();
      assert(volvis_ij   == NULL); volvis_ij   = new int[volvis_size];
      assert(volvis_data == NULL); volvis_data = new float[volvis_size][2];
      for (int ii = 0; ii < volvis_size; ++ii) {
        volvis_ij[ii]        = my_ij_vec[ii];
        volvis_data[ii][0] = my_data_vec[ii].intensity;
        volvis_data[ii][1] = my_data_vec[ii].depth;
      }

    }
    else {

      int * count = NULL;
      if (mpi_rank == write_rank) count = new int[mpi_size];
      int my_volvis_size = my_ij_vec.size();
      MPI_Gather(&my_volvis_size,1,MPI_INT,count,1,MPI_INT,write_rank,mpi_comm);

      int * disp = NULL;
      if (mpi_rank == write_rank) {
        disp = new int[mpi_size];
        disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank) disp[rank] = disp[rank-1] + count[rank-1];
        volvis_size = disp[mpi_size-1] + count[mpi_size-1];
        assert(volvis_ij   == NULL); volvis_ij   = new int[volvis_size];
        assert(volvis_data == NULL); volvis_data = new float[volvis_size][2];
      }

      MPI_Gatherv(&my_ij_vec.front(),my_volvis_size,MPI_INT,volvis_ij,count,disp,MPI_INT,write_rank,mpi_comm);

      if (mpi_rank == write_rank) {
        FOR_RANK {
          count[rank] *= 2;
          disp[rank] *= 2;
        }
      }

      MPI_Gatherv(&my_data_vec.front().intensity,my_volvis_size*2,MPI_FLOAT,(float*)volvis_data,count,disp,MPI_FLOAT,write_rank,mpi_comm);

      if (mpi_rank == write_rank) {
        delete[] count;
        delete[] disp;
      }

    }

    my_ij_vec.clear();
    my_data_vec.clear();

  } // if (b_volvis)...

  // ==============================================================================
  // this is a logical point to break image rendering and call a second routine
  // later if desired...
  // ==============================================================================

}

void CtiCanvas::flushImage(const string& filename, const RGB_MODE rgb_mode) {

  if (mpi_rank == 0)
    cout << "CtiCanvas::writeImage() " << filename << endl;

  if (mpi_rank == write_rank) {

    // serial render of the image...

    unsigned char (*rgb)[3] = new unsigned char[ni*nj][3];

    PngDataChunk *dpTh = new PngDataChunk("dpTh", ni, nj, 2);
    unsigned char dpTh_uc[2];

    PngDataChunk *zoNe = new PngDataChunk("zoNe", ni, nj, 2);
    unsigned char zoNe_uc[2];

    bool b_volumeData = false;
    bool b_surfaceData = false;
    bool b_particleData = false;
    bool b_isoData = false;
    // for euclidean distance reduction for moving solver and mesh point viz
    bool b_auxData = false;
    // for masking
    bool b_maskData = false;
    bool b_internalSurface = false;
    bool b_delayData = false;
    bool b_SauxData = false; // PAUX used for surface data tie breaking

    if (!b_volvis) {
      // set famous blue in parts of the image where no information exists...
      dpTh_uc[0] = (unsigned char)255;
      dpTh_uc[1] = (unsigned char)255;
      zoNe_uc[0] = (unsigned char)255;
      zoNe_uc[1] = (unsigned char)255;
      for (int ij = 0; ij < ni*nj; ++ij) {
        rgb[ij][0] = 73;
        rgb[ij][1] = 175;
        rgb[ij][2] = 205;
        dpTh->set(ij,dpTh_uc);
        zoNe->set(ij,zoNe_uc);
      }
    }
    else {
      // set black to emphasize volvis...
      dpTh_uc[0] = (uint2)65534 & 255;  // default to background
      dpTh_uc[1] = (uint2)65534 >> 8;
      zoNe_uc[0] = (uint2)65535 & 255;  // background zone
      zoNe_uc[1] = (uint2)65535 >> 8;
      for (int ij = 0; ij < ni*nj; ++ij) {
        FOR_I3 rgb[ij][i] = (unsigned char) 0;
        dpTh->set(ij,dpTh_uc);
        zoNe->set(ij,zoNe_uc);
      }
    }

    // lighting settings...
    const float light1[3] = { 0.0f, 0.0f, 1.0f };
    const float light2[3] = { -0.57f, -0.57f, 0.57f };
    const float base3[3] = { 0.99f, 0.96f, 0.89f }; // * color
    const float base2[3] = { 0.92f, 0.91f, 0.83f };
    const float base00[3] = { 0.40f, 0.48f, 0.51f };

    // allocated always...
    uint2 * pixel_flag       = new uint2[ni*nj];
    float (*pixel_normal)[3] = new float[ni*nj][3]; // used for surface or internal data
    uint2 * pixel_zone       = new uint2[ni*nj]; // used for surface or internal data
 
    // allocated when needed...
    float * pixel_data       = NULL; // used for surface or internal data
    float * pixel_sdepth     = NULL;
    float * pixel_mdepth     = NULL;
    float * pixel_pdepth     = NULL;
    float * pixel_paux       = NULL;

    if (!projPlane){// default to background
      for (int ij = 0; ij < ni*nj; ++ij) {
        pixel_flag[ij] = BACKGROUND_PIXEL; // default to background
        FOR_K3 pixel_normal[ij][k] = HUGE_VALF;
        pixel_zone[ij] = 65535;
      }
    }
    else{ //default to iso-plane, whe projPlane is not null
      if (pixel_pdepth==NULL)
        pixel_pdepth = new float[ni*nj];
      if (pixel_data==NULL)
        pixel_data = new float[ni*nj];

      for (int i=0; i<ni; ++i){
        for (int j=0; j<nj; ++j){
          int ij = j*ni + i;//DAPDAP
          FOR_K3 pixel_normal[ij][k] = projPlane->normal[k];
          pixel_pdepth[ij] = projPlane->center[2] - (projPlane->normal[0]*(i-projPlane->center[0]) +
                                                     projPlane->normal[1]*(j-projPlane->center[1]))/projPlane->normal[2];
          b_isoData = true;
          pixel_flag[ij] = ISO_DATA_PIXEL;
          pixel_zone[ij] = 65532;
          pixel_data[ij] = 0.0f;
        }
      }
    }

    // fist big loop through the pbuf to:
    // 1. figure out which of the above arrays are needed
    // 2. perform surface tie breaking
    // 3.

    float min_depth = HUGE_VALF;
    for (int ip = 0; ip < pbuf_size; ++ip) {
      const int ij = PBUF_J*ni + PBUF_I;
      // every pixel sent from the ranks should have some meaning!...
      assert(PBUF_BITS);
      // surface stuff...
      if (PBUF_BITS&SURFACE_BITS) {
        // can't be internal as well!...
        //assert(!(PBUF_BITS&INTERNAL_BITS)); // yes it can!
        if (pixel_sdepth == NULL) {
          pixel_sdepth = new float[ni*nj];
          for (int ij_ = 0; ij_ < ni*nj; ++ij_) pixel_sdepth[ij_] = HUGE_VALF;
        }
        if (PBUF_BITS&SURFACE_DATA_BIT) {
          // we have surface data...
          b_surfaceData = true;
          if (pixel_data == NULL) {
            pixel_data = new float[ni*nj];
            for (int ij_ = 0; ij_ < ni*nj; ++ij_) pixel_data[ij_] = HUGE_VALF;
          }
          // surface data may be associated with the paux containing an
          // additional tie-breaking distance (squared)...
          if (PBUF_BITS&AUX_BIT) {
            // assert(!b_auxData);  //HACK can we simply ignore?
            b_SauxData = true; // we can make all decisions in this loop...
            if (pixel_paux == NULL) {
              pixel_paux = new float[ni*nj];
              for (int ij_ = 0; ij_ < ni*nj; ++ij_) pixel_paux[ij_] = HUGE_VALF;
            }
            if ( (PBUF_SDEPTH < pixel_sdepth[ij]) || ((PBUF_SDEPTH == pixel_sdepth[ij])&&(PBUF_PAUX < pixel_paux[ij])) ) {
            //if ( (PBUF_SDEPTH <= pixel_sdepth[ij]) ) {
              pixel_sdepth[ij] = PBUF_SDEPTH;
              pixel_paux[ij] = PBUF_PAUX;
              pixel_normal[ij][0] = PBUF_SNORMAL0;
              pixel_normal[ij][1] = PBUF_SNORMAL1;
              pixel_normal[ij][2] = PBUF_SNORMAL2;
              min_depth = min(min_depth,pixel_sdepth[ij]);
              pixel_flag[ij] = SURFACE_DATA_PIXEL;
              pixel_data[ij] = PBUF_SDATA;
              pixel_zone[ij] = PBUF_SZONE+32768; // could be either and edge index (>= 32768) or subzone
            }
          }
          else if (PBUF_SDEPTH < pixel_sdepth[ij]) {
            // surface data without the paux, and distance has won...
            pixel_sdepth[ij] = PBUF_SDEPTH;
            pixel_normal[ij][0] = PBUF_SNORMAL0;
            pixel_normal[ij][1] = PBUF_SNORMAL1;
            pixel_normal[ij][2] = PBUF_SNORMAL2;
            min_depth = min(min_depth,pixel_sdepth[ij]);
            pixel_flag[ij] = SURFACE_DATA_PIXEL;
            pixel_data[ij] = PBUF_SDATA;
            pixel_zone[ij] = PBUF_SZONE+32768; // could be either and edge index (>= 32768) or subzone
            // assert(pixel_paux == NULL); // if you hit this, multiple ways of laying down surface data?
          }
        }
        else {
          assert(PBUF_BITS&SURFACE_BIT);
          // no surface data, so just use depth to decide...
          if (PBUF_SDEPTH < pixel_sdepth[ij]) {
            pixel_sdepth[ij] = PBUF_SDEPTH;
            pixel_normal[ij][0] = PBUF_SNORMAL0;
            pixel_normal[ij][1] = PBUF_SNORMAL1;
            pixel_normal[ij][2] = PBUF_SNORMAL2;
            min_depth = min(min_depth,pixel_sdepth[ij]);
            pixel_flag[ij] = SURFACE_PIXEL;
            pixel_zone[ij] = PBUF_SZONE; // could be either and edge index (>= 32768) or subzone
            // surface without data has just won the depth check. Reset
            // saux data (which uses pixel_paux) to store the data distance to the
            // surface pixel...
            if ((pixel_paux)&&(b_SauxData)) pixel_paux[ij] = HUGE_VALF;
          }
        }
        if (PBUF_BITS&DELAY_BIT)
          b_delayData = true;
      }

      if (PBUF_BITS&INTERNAL_BITS) {
        // internal stuff...
        if (pixel_pdepth == NULL) {
          pixel_pdepth = new float[ni*nj];
          for (int ij_ = 0; ij_ < ni*nj; ++ij_) pixel_pdepth[ij_] = HUGE_VALF;
        }
        pixel_pdepth[ij] = min(pixel_pdepth[ij],PBUF_PDEPTH); // smallest plane/point/iso depth wins
        min_depth = min(min_depth,pixel_pdepth[ij]);

        if (PBUF_BITS&INTERNAL_DATA_BIT) {
          if (pixel_data == NULL) {
            pixel_data = new float[ni*nj];
            for (int ij_ = 0; ij_ < ni*nj; ++ij_) pixel_data[ij_] = HUGE_VALF;
          }
          if (PBUF_PZONE == 65532)
            b_isoData = true;
          else if (PBUF_PZONE == 65533)
            b_particleData = true;
          else { // could be 65534 or something less than 32768 (ilevel/iwindow)
            assert((PBUF_PZONE == 65534)||(PBUF_PZONE < 32767));
            b_volumeData = true;
            if (PBUF_BITS&MASK_BIT) {
              b_maskData = true;
              if (pixel_mdepth == NULL) {
                pixel_mdepth = new float[ni*nj];
                for (int ij_ = 0; ij_ < ni*nj; ++ij_) pixel_mdepth[ij_] = HUGE_VALF;
              }
              // have to set mdepth in a second pass
            }
            if (PBUF_BITS&AUX_BIT) {
              // assert(!b_SauxData);  //HACK: ok to remove?
              b_auxData = true;
              if (pixel_paux == NULL) {
                pixel_paux = new float[ni*nj];
                for (int ij_ = 0; ij_ < ni*nj; ++ij_) pixel_paux[ij_] = HUGE_VALF;
              }
              pixel_paux[ij] = min(pixel_paux[ij],PBUF_PAUX); // get closest to plane
            }
            if (PBUF_BITS&DELAY_BIT)
              b_delayData = true;
          }
        }
        else {
          assert(PBUF_BITS&ISURFACE_BIT);
          b_internalSurface = true;
          assert(PBUF_PZONE <= 32767);
          if (PBUF_BITS&MASK_BIT) {
            //assert(PBUF_BITS&MESH_BIT); // if you hit this that means there is a new case
            // zone stores ilevel/iwindow
            b_maskData = true;
            if (pixel_mdepth == NULL) {
              pixel_mdepth = new float[ni*nj];
              for (int ij_ = 0; ij_ < ni*nj; ++ij_) pixel_mdepth[ij_] = HUGE_VALF;
            }
          }
        }
      }
    }

    // if we used the Saux, then we have to think more carefully about the
    // downstream use of paux_data...
    if (b_SauxData) {
      assert(b_surfaceData);
      // assert(!b_auxData);
      // may not be neccesary...
      // assert(!b_volumeData);
      // assert(!b_particleData);
      // assert(!b_maskData);
      // assert(!b_isoData);
    }

    // include a daTa chunk if using rgb channel for NORMALS (images for cascade app client lighting)
    // or if any lit data is included (anything but planar data)
    // **if lit data is included, also set daTa for planar data

    PngDataChunk *daTa = NULL;
    PngDataChunk *lgFr = NULL;

    if (rgb_mode==NORMALS||b_surfaceData||b_particleData||b_isoData) {
      daTa = new PngDataChunk("daTa",ni,nj,1);
      unsigned char data_uc = (unsigned char) 255;
      for (int ij = 0; ij < ni*nj; ++ij) daTa->set(ij,&data_uc);
      lgFr = new PngDataChunk("lgFr",ni,nj,1);
      for (int ij = 0; ij < ni*nj; ++ij) lgFr->set(ij,&data_uc); 
    }

    if (!b_volvis) {

      // ====================================================
      // regular visualization (i.e. NOT volvis)...
      // ====================================================

      // masking data has to be done after, because of coordination b/w MASK_BIT and MSURFACE_BIT
      if (b_maskData) {
        assert(pixel_mdepth);
        for (int ip = 0; ip < pbuf_size; ++ip) {
          const int ij = PBUF_J*ni + PBUF_I;
          if ((PBUF_BITS&MSURFACE_BIT)&&(fabs(PBUF_MDEPTH) < fabs(pixel_mdepth[ij])))
            pixel_mdepth[ij] = PBUF_MDEPTH;
        }
      }

      // currently data is cell flooded with a global index. take grad of this to find edges (set to 1 or zero).
      // then do laplacian smoothing to thicken edge and to dealias.
      const int niter = 1;
      const float grid_delta = (float)niter;
      const float mesh_fraction = 0.5f;
      if (b_delayData) {
        assert(pixel_data);
        assert(pixel_paux);
        float* pixel_pdata = new float[ni*nj];;
        float* pixel_pdata2 = new float[ni*nj];
        // copy data over to matrix, will post process down below...
        for (int ij = 0; ij < ni*nj; ++ij) {
          pixel_pdata[ij] = -1.f;
          pixel_pdata2[ij] = -1.f;
        }
        if (!b_SauxData) {
          assert(pixel_mdepth);
          assert(pixel_pdepth);
          for (int ip = 0; ip < pbuf_size; ++ip) {
            const int ij = PBUF_J*ni + PBUF_I;
            if ((PBUF_PDEPTH == pixel_pdepth[ij])&&(PBUF_PAUX == pixel_paux[ij])&&(pixel_mdepth[ij] < 0.f)&&((pixel_sdepth == NULL)||(pixel_pdepth[ij] < pixel_sdepth[ij]))) {
              pixel_pdata[ij] = PBUF_PDATA; // should be the unique id
              assert(pixel_pdata[ij] >= 0.f);
            }
          }
        }
        else {
          assert(pixel_sdepth);
          // TODO should I add a SURFACE_MESH_PIXEL so that I do not copy all data?
          for (int ij = 0; ij < ni*nj; ++ij)
            pixel_pdata[ij] = pixel_data[ij];
        }
        // find edges by finding non-zero data gradient (use backward grad)...
        for (int i = 0; i < ni; ++i) {
          for (int j = 0; j < nj; ++j) {
            const int ij = j*ni+i;
            if (pixel_pdata[ij] != -1.f) {
              float dpix = 0.f;
              if ((i > 0)&&(pixel_pdata[j*ni+(i-1)] > pixel_pdata[ij]))
                dpix += abs(pixel_pdata[ij]-pixel_pdata[j*ni+(i-1)]);
              if ((i < ni-1)&&(pixel_pdata[j*ni+(i+1)] > pixel_pdata[ij]))
                dpix += abs(pixel_pdata[ij]-pixel_pdata[j*ni+(i+1)]);
              if ((j > 0)&&(pixel_pdata[(j-1)*ni+i] > pixel_pdata[ij]))
                dpix += abs(pixel_pdata[ij]-pixel_pdata[(j-1)*ni+i]);
              if ((j < nj-1)&&(pixel_pdata[(j+1)*ni+i] > pixel_pdata[ij]))
                dpix += abs(pixel_pdata[ij]-pixel_pdata[(j+1)*ni+i]);
              if (dpix > 0.f)
                pixel_pdata2[ij] = 0.f; // edges black
              else
                pixel_pdata2[ij] = 1.f; // interior white
            }
          }
        }
        // laplacian smoothing (change number of iterations to change thickness)...
        float *original = pixel_pdata;
        float *filtered = pixel_pdata2;
        for (int iter = 0; iter < niter; ++iter) {
          // left/right buffers flip every iteration
          if (iter%2 == 0) {
            original = pixel_pdata2;
            filtered = pixel_pdata;
          }
          else {
            original = pixel_pdata;
            filtered = pixel_pdata2;
          }
          for (int i = 0; i < ni; ++i) {
            for (int j = 0; j < nj; ++j) {
              const int ij = j*ni+i;
              if (original[ij] > 0.f) {
                int cnt = 0;
                float sum = 0.f;
                if ((i > 0)&&(original[j*ni+(i-1)] != -1.f)) {
                  sum += original[j*ni+(i-1)];
                  ++cnt;
                }
                if ((i < ni-1)&&(original[j*ni+(i+1)] != -1.f)) {
                  sum += original[j*ni+(i+1)];
                  ++cnt;
                }
                if ((j > 0)&&(original[(j-1)*ni+i] != -1.f)) {
                  sum += original[(j-1)*ni+i];
                  ++cnt;
                }
                if ((j < nj-1)&&(original[(j+1)*ni+i] != -1.f)) {
                  sum += original[(j+1)*ni+i];
                  ++cnt;
                }
                if (cnt > 0)
                  filtered[ij] = sum/(float)cnt;
              }
              else
                filtered[ij] = original[ij];
            }
          }
        }
        // use sqrt(r) kernal for filter...
        for (int ij = 0; ij < ni*nj; ++ij) {
          if (filtered[ij] >= 0.f) {
            if (!b_SauxData)
              pixel_flag[ij] = MASKED_AUX_VOLUME_PIXEL; // switch data to normal below
            pixel_data[ij] = filtered[ij];
          }
        }
        delete[] pixel_pdata;
        delete[] pixel_pdata2;
      }

      // now set pixel flag/data based on internal data

      for (int ip = 0; ip < pbuf_size; ++ip) {
        const int ij = PBUF_J*ni + PBUF_I;
        if ((PBUF_BITS&INTERNAL_BITS)&&(PBUF_PDEPTH == pixel_pdepth[ij])&&((pixel_sdepth == NULL)||(pixel_pdepth[ij] < pixel_sdepth[ij]))) {
          if (PBUF_BITS&INTERNAL_DATA_BIT) {
            if (b_volumeData) {
              if (b_maskData) {
                if (pixel_mdepth[ij] < 0.f) {
                  if (b_auxData) {
                    if (PBUF_PAUX == pixel_paux[ij]) {
                      pixel_normal[ij][0] = PBUF_PNORMAL0;
                      pixel_normal[ij][1] = PBUF_PNORMAL1;
                      pixel_normal[ij][2] = PBUF_PNORMAL2;

                      if (!b_delayData) {
                        pixel_data[ij] = PBUF_PDATA;
                        pixel_flag[ij] = MASKED_AUX_VOLUME_DATA_PIXEL;
                      }

                      if (PBUF_BITS&MESH_BIT)
                        pixel_zone[ij] = PBUF_PZONE+32768;
                      else
                        pixel_zone[ij] = PBUF_PZONE;
                    }
                  }
                  else {
                    pixel_flag[ij] = MASKED_VOLUME_DATA_PIXEL;
                    pixel_data[ij] = PBUF_PDATA;
                    pixel_normal[ij][0] = PBUF_PNORMAL0;
                    pixel_normal[ij][1] = PBUF_PNORMAL1;
                    pixel_normal[ij][2] = PBUF_PNORMAL2;
                    pixel_zone[ij]      = PBUF_PZONE;
                  }
                }
              }
              else if (PBUF_PZONE == 65534) {
                pixel_flag[ij] = VOLUME_DATA_PIXEL;
                pixel_data[ij] = PBUF_PDATA;
                pixel_normal[ij][0] = PBUF_PNORMAL0;
                pixel_normal[ij][1] = PBUF_PNORMAL1;
                pixel_normal[ij][2] = PBUF_PNORMAL2;
                pixel_zone[ij]      = PBUF_PZONE;
              }
            }
            if ((b_isoData)&&(PBUF_PZONE == 65532)) {
              pixel_flag[ij] = ISO_DATA_PIXEL;
              pixel_data[ij] = PBUF_PDATA;
              pixel_normal[ij][0] = PBUF_PNORMAL0;
              pixel_normal[ij][1] = PBUF_PNORMAL1;
              pixel_normal[ij][2] = PBUF_PNORMAL2;
              pixel_zone[ij]      = PBUF_PZONE;
            }
            if ((b_particleData)&&(PBUF_PZONE == 65533)) {
              pixel_flag[ij] = PARTICLE_DATA_PIXEL;
              pixel_data[ij] = PBUF_PDATA;
              pixel_normal[ij][0] = PBUF_PNORMAL0;
              pixel_normal[ij][1] = PBUF_PNORMAL1;
              pixel_normal[ij][2] = PBUF_PNORMAL2;
              pixel_zone[ij]      = PBUF_PZONE;
            }
          }
          else {
            assert(b_internalSurface);
            assert(PBUF_BITS&ISURFACE_BIT);
            if (b_maskData&&(pixel_mdepth[ij] < 0.f)&&(PBUF_BITS&MESH_BIT)) {
              pixel_flag[ij] = MASKED_VOLUME_PIXEL;
              pixel_zone[ij] = 32768+PBUF_PZONE;
              pixel_normal[ij][0] = PBUF_PNORMAL0;
              pixel_normal[ij][1] = PBUF_PNORMAL1;
              pixel_normal[ij][2] = PBUF_PNORMAL2;
            }
            else if ((PBUF_PZONE == 32767) || (PBUF_PZONE == 32766)) {
              // zones distinguish between fluid/cht internal surface
              pixel_flag[ij] = INTERNAL_SURFACE_PIXEL;
              pixel_zone[ij] = PBUF_PZONE;
              pixel_normal[ij][0] = PBUF_PNORMAL0;
              pixel_normal[ij][1] = PBUF_PNORMAL1;
              pixel_normal[ij][2] = PBUF_PNORMAL2;
            }
          }
        }
      }

      if (b_delayData) {
        for (int ij = 0; ij < ni*nj; ++ij) {
          if ((pixel_flag[ij] == MASKED_AUX_VOLUME_PIXEL)||((b_SauxData)&&(pixel_flag[ij] == SURFACE_DATA_PIXEL))) {
            //if (!((pixel_data[ij] >= 0.f)&&(pixel_data[ij] <= 1.f)))
            //  cout << "ERROR: pixel_data[ij]: " << pixel_data[ij] << endl; // TODO fix this
            const float dn = pixel_data[ij];
            const bool is_edge = (dn < grid_delta);
            if (is_edge) {
              const float factor = mesh_fraction*0.5f+0.5f + (mesh_fraction*0.5f-0.5f)*cos(dn/grid_delta*M_PI);
              FOR_I3 pixel_normal[ij][i] *= factor;
            }
            pixel_data[ij] = 1.f;
          }
        }
        if (!b_SauxData) {
          range[0] = range[1] = 1.f;
          b_range = true;
        }
        else {
          range_sVar[0] = range_sVar[1] = 1.f;
          b_range_sVar = true;
        }
      }

      //compute variable range based only on pixels that win depth check after reduction...

      if (!b_range){
        range[0]    =  HUGE_VALF;
        range[1]    = -HUGE_VALF;
      }
      if (!b_range_sVar) {
        range_sVar[0] =  HUGE_VALF;
        range_sVar[1] = -HUGE_VALF;
      }
      if (!b_range_pVar) {
        range_pVar[0] =  HUGE_VALF;
        range_pVar[1] = -HUGE_VALF;
      }
      if (!b_range_iVar) {
        range_iVar[0] =  HUGE_VALF;
        range_iVar[1] = -HUGE_VALF;
      }

      //compute variable range based only on pixels that win depth check after reduction...

      for (int ij = 0; ij < ni*nj; ++ij) {
        if (pixel_flag[ij] == SURFACE_DATA_PIXEL && !b_range_sVar) {
          range_sVar[0] = min(range_sVar[0],pixel_data[ij]);
          range_sVar[1] = max(range_sVar[1],pixel_data[ij]);
        }
        else if (pixel_flag[ij] == PARTICLE_DATA_PIXEL && !b_range_pVar) {
          range_pVar[0] = min(range_pVar[0],pixel_data[ij]);
          range_pVar[1] = max(range_pVar[1],pixel_data[ij]);
        }
        else if (pixel_flag[ij] == ISO_DATA_PIXEL && !b_range_iVar) {
          range_iVar[0] = min(range_iVar[0],pixel_data[ij]);
          range_iVar[1] = max(range_iVar[1],pixel_data[ij]);
        }
        else if (!b_range && ((pixel_flag[ij] == VOLUME_DATA_PIXEL) ||
                              (pixel_flag[ij] == MASKED_VOLUME_DATA_PIXEL) ||
                              (pixel_flag[ij] == MASKED_AUX_VOLUME_DATA_PIXEL))) {
          range[0] = min(range[0],pixel_data[ij]);
          range[1] = max(range[1],pixel_data[ij]);
        }
      }

      // precompute denominator constant for ranges
      float denom = 0.f;
      if (range[1] != range[0]) denom = 1.f/(range[1]-range[0]);
      float denom_p = 0.f;
      if (range_pVar[1] != range_pVar[0]) denom_p = 1.f/(range_pVar[1]-range_pVar[0]);
      float denom_s = 0.f;
      if (range_sVar[1] != range_sVar[0]) denom_s = 1.f/(range_sVar[1]-range_sVar[0]);
      float denom_i = 0.f;
      if (range_iVar[1] != range_iVar[0]) denom_i = 1.f/(range_iVar[1]-range_iVar[0]);

      // bins for contour computing
      float f_bin[2] = {1.f,1.f};
      float f_sbin[2] = {1.f,1.f};
      float f_ibin[2] = {1.f,1.f};
      if (n_range) {
        f_bin[0] = 255.f/n_range;
        f_bin[1] = n_range/255.f;
      }
      if (n_range_sVar) {
        f_sbin[0] = 255.f/n_range_sVar;
        f_sbin[1] = n_range_sVar/255.f;
      }
      if (n_range_iVar) {
        f_ibin[0] = 255.f/n_range_iVar;
        f_ibin[1] = n_range_iVar/255.f;
      }

      for (int ij = 0; ij < ni*nj; ++ij) {
        if (pixel_flag[ij] == BACKGROUND_PIXEL) {
          continue;
        }
        else {
          // set the zone...
          uint2 zoNe_us = pixel_zone[ij];
          zoNe_uc[0] = zoNe_us & 255;
          zoNe_uc[1] = zoNe_us >> 8;
          zoNe->set(ij,zoNe_uc);

          // set the depth...
          switch (pixel_flag[ij]) {
          case MASKED_VOLUME_PIXEL:
          case MASKED_AUX_VOLUME_PIXEL:
          case PARTICLE_DATA_PIXEL:
          case ISO_DATA_PIXEL:
          case VOLUME_DATA_PIXEL:
          case MASKED_VOLUME_DATA_PIXEL:
          case MASKED_AUX_VOLUME_DATA_PIXEL:
          case INTERNAL_SURFACE_PIXEL:
            {
              uint2 dpTh_us = (uint2)max(0.f,min(65534.f,pixel_pdepth[ij]-min_depth));
              dpTh_uc[0] = dpTh_us & 255;
              dpTh_uc[1] = dpTh_us >> 8;
              dpTh->set(ij,dpTh_uc);
              break;
            }
          case SURFACE_PIXEL:
          case SURFACE_DATA_PIXEL:
            {
              uint2 dpTh_us = (uint2) max(0.f,min(65534.f,pixel_sdepth[ij]-min_depth));
              dpTh_uc[0] = dpTh_us & 255;
              dpTh_uc[1] = dpTh_us >> 8;
              dpTh->set(ij,dpTh_uc);
              break;
            }
          default:
            assert(0);
          }

          // set rgb to normal in NORMALS rgb_mode...
          if (rgb_mode == NORMALS) {
            rgb[ij][0] = (unsigned char) max(0,min(255,int((pixel_normal[ij][0]+1.0f)*128.0f)));
            rgb[ij][1] = (unsigned char) max(0,min(255,int((pixel_normal[ij][1]+1.0f)*128.0f)));
            rgb[ij][2] = (unsigned char) max(0,min(255,int((pixel_normal[ij][2]+1.0f)*128.0f)));
          }

          // set data
          switch (pixel_flag[ij]) {
          case PARTICLE_DATA_PIXEL:
            {
              // clean handling of the case when range_pVar[0] == range_pVar[1]...
              unsigned char daTa_uc;
              if (denom_p == 0.f) {  // range_pVar[0] == range_pVar[1]
                if      (pixel_data[ij] < range_pVar[0]) daTa_uc = (unsigned char)0;
                else if (pixel_data[ij] > range_pVar[0]) daTa_uc = (unsigned char)255;
                else                                     daTa_uc = (unsigned char)127;
              }
              else {
                //Expected by image diagnostics: normalized data 0-1 uses 255 bins of precision, endpoints centered at bins 0, 255
                daTa_uc = (unsigned char)max(0,min(255,int((pixel_data[ij]-range_pVar[0])*denom_p*255.0f+0.5f)));
              }
              assert(daTa);
              daTa->set(ij,&daTa_uc);
              if (rgb_mode == DATA) {
                const float sign = (pixel_normal[ij][2] >= 0.f) ? 1.0f:-1.0f;
                const float a = sign*(pixel_normal[ij][0]*light1[0] + pixel_normal[ij][1]*light1[1] + pixel_normal[ij][2]*light1[2]);
                const float b = sign*(pixel_normal[ij][0]*light2[0] + pixel_normal[ij][1]*light2[1] + pixel_normal[ij][2]*light2[2]);
                const float p = (a*base2[0] + (1.0-a)*base00[0])*0.5 + (b*base3[0] + (1.0-b)*base00[0])*0.5;
                unsigned char lgFr_uc = (unsigned char)max(0,min(255,int(p*255.0f+0.5f)));
                lgFr->set(ij,&lgFr_uc);

                //rgb[ij][0] = daTa_uc;
                //apply lighting to data in 256 bins...
                rgb[ij][0] = max(0,min(255,(int)(p*daTa_uc)));
                // for grayscale, r = g = b...
                rgb[ij][1] = rgb[ij][0];
                rgb[ij][2] = rgb[ij][0];
              }
              //else if (rgb_mode == NORMALS) {
              //  assert(daTa);
              //  daTa->set(ij,&daTa_uc);
              //}
              //else {
              //  assert(0);
              //}
              break;
            }
          case ISO_DATA_PIXEL:
            {
              unsigned char daTa_uc;
              if (denom_i == 0.f) {  // range_iVar[0] == range_iVar[1]
                if      (pixel_data[ij] < range_iVar[0]) daTa_uc = (unsigned char)0;
                else if (pixel_data[ij] > range_iVar[0]) daTa_uc = (unsigned char)255;
                else                                     daTa_uc = (unsigned char)127;
              }
              else {
                //Expected by image diagnostics: normalized data 0-1 uses 255 bins of precision, endpoints centered at bins 0, 255
                float d_tmp = max(0,min(255,int((pixel_data[ij]-range_iVar[0])*denom_i*255.0f+0.5f)));

                if (n_range_iVar) {
                  d_tmp = float(int(d_tmp*f_ibin[1]))*f_ibin[0];
                }
                daTa_uc = (unsigned char)d_tmp;
              }
              assert(daTa);
              daTa->set(ij,&daTa_uc);
              if (rgb_mode == DATA) {
                const float sign = (pixel_normal[ij][2] >= 0.f) ? 1.0f:-1.0f;
                const float a = sign*(pixel_normal[ij][0]*light1[0] + pixel_normal[ij][1]*light1[1] + pixel_normal[ij][2]*light1[2]);
                const float b = sign*(pixel_normal[ij][0]*light2[0] + pixel_normal[ij][1]*light2[1] + pixel_normal[ij][2]*light2[2]);
                const float p = (a*base2[0] + (1.0-a)*base00[0])*0.5 + (b*base3[0] + (1.0-b)*base00[0])*0.5;
                unsigned char lgFr_uc = (unsigned char)max(0,min(255,int(p*255.0f+0.5f)));
                lgFr->set(ij,&lgFr_uc);
                //rgb[ij][0] = daTa_uc;
                //apply lighting to data in 256 bins...
                rgb[ij][0] = max(0,min(255,(int)(p*daTa_uc)));
                rgb[ij][1] = rgb[ij][0];
                rgb[ij][2] = rgb[ij][0];
              }
              //else if (rgb_mode == NORMALS) {
              //  assert(daTa);
              //  daTa->set(ij,&daTa_uc);
              //}
              //else {
              //  assert(0);
              //}
              break;
            }
          case VOLUME_DATA_PIXEL:
          case MASKED_VOLUME_DATA_PIXEL:
          case MASKED_AUX_VOLUME_DATA_PIXEL:
            {
              unsigned char daTa_uc;
              if (denom == 0.f) {  // range[0] == range[1]
                if      (pixel_data[ij] < range[0]) daTa_uc = (unsigned char)0;
                else if (pixel_data[ij] > range[0]) daTa_uc = (unsigned char)255;
                else                                daTa_uc = (unsigned char)127;
              }
              else {
                //Expected by image diagnostics: normalized data 0-1 uses 255 bins of precision, endpoints centered at bins 0, 255
                // daTa_uc = (unsigned char)max(0,min(255,int((pixel_data[ij]-range[0])/(range[1]-range[0])*255.0f+0.5f)));
                float d_tmp = max(0,min(255,int((pixel_data[ij]-range[0])*denom*255.0f+0.5f)));

                if (n_range) {
                  d_tmp = float(int(d_tmp*f_bin[1]))*f_bin[0];
                }
                daTa_uc = (unsigned char)d_tmp;
              }
              if (rgb_mode == DATA) {
                rgb[ij][0] = daTa_uc;
                rgb[ij][1] = rgb[ij][0];
                rgb[ij][2] = rgb[ij][0];
                if (daTa) {
                  daTa->set(ij,&daTa_uc); //if daTa is not null due to other lit data, set it here as well.
                  unsigned char lgFr_uc = (unsigned char)255;
                  lgFr->set(ij,&lgFr_uc);
                }
              }
              else if (rgb_mode == NORMALS) {
                assert(daTa);
                daTa->set(ij,&daTa_uc);
              }
              else {
                assert(0);
              }
              break;
            }
          case MASKED_VOLUME_PIXEL:
          case MASKED_AUX_VOLUME_PIXEL:
            {
              if (rgb_mode == DATA) {
                // should we light surface/iso/particle data with normal (or does that screw up recoloring DAP)?
                const float sign = (pixel_normal[ij][2] >= 0.f) ? 1.0f:-1.0f;
                const float a = sign*(pixel_normal[ij][0]*light1[0] + pixel_normal[ij][1]*light1[1] + pixel_normal[ij][2]*light1[2]);
                const float b = sign*(pixel_normal[ij][0]*light2[0] + pixel_normal[ij][1]*light2[1] + pixel_normal[ij][2]*light2[2]);
                const float p = (a*base2[0] + (1.0-a)*base00[0])*0.5 + (b*base3[0] + (1.0-b)*base00[0])*0.5;
                rgb[ij][0] = max(0,min(255,(int)(p*256.0)));
                rgb[ij][1] = rgb[ij][0];
                rgb[ij][2] = rgb[ij][0];
              }
              else if (rgb_mode == NORMALS) {
                ;
              }
              else {
                assert(0);
              }
              break;
            }
          case INTERNAL_SURFACE_PIXEL:
          case SURFACE_PIXEL:
            {
              if (rgb_mode == DATA) {
                const float sign = ((pixel_normal[ij][2] >= 0.f) || (pixel_zone[ij] >= 32768)) ? 1.0f:-1.0f;  // edge rendering doesn't care about normal direction

                const float a = sign*(pixel_normal[ij][0]*light1[0] + pixel_normal[ij][1]*light1[1] + pixel_normal[ij][2]*light1[2]);
                const float b = sign*(pixel_normal[ij][0]*light2[0] + pixel_normal[ij][1]*light2[1] + pixel_normal[ij][2]*light2[2]);
                const float p = (a*base2[0] + (1.0-a)*base00[0])*0.5 + (b*base3[0] + (1.0-b)*base00[0])*0.5;

                if (pixel_zone[ij] < 32768) {
                  // regular surface zone
                  if (sign >= 0.f) {
                    // silver for wetted surface...
                    rgb[ij][0] = max(0,min(255,(int)(p*256.0)));
                    rgb[ij][1] = rgb[ij][0];
                    rgb[ij][2] = rgb[ij][0];
                  }
                  else {
                    // gold for unwetted surface...
                    rgb[ij][0] = max(0,min(255,(int)(p*256.0)));
                    rgb[ij][1] = max(0,min(255,(int)(0.93*p*256.0)));
                    rgb[ij][2] = max(0,min(255,(int)(0.72*p*256.0)));
                  }
                }
                else if (pixel_zone[ij] == 32768) {
                  // turquoise for misaligned edge...
                  rgb[ij][0] = max(0,min(255,(int)(p*120.0)));
                  rgb[ij][1] = max(0,min(255,(int)(p*211.0)));
                  rgb[ij][2] = max(0,min(255,(int)(p*120.0)));
                }
                else if (pixel_zone[ij] < 49152) {
                  // blue for open edges...
                  rgb[ij][0] = 0;
                  rgb[ij][1] = 0;
                  rgb[ij][2] = max(0,min(255,(int)(p*256.0)));
                }
                else {
                  assert(pixel_zone[ij] < 65535);
                  // purple for multiedges...
                  rgb[ij][0] = max(0,min(255,(int)(p*256.0)));
                  rgb[ij][1] = 0;
                  rgb[ij][2] = max(0,min(255,(int)(p*256.0)));
                }
              }
              else if (rgb_mode == NORMALS) {
                ;
              }
              else {
                assert(0);
              }
              break;
            }
          case SURFACE_DATA_PIXEL:
            {
              unsigned char daTa_uc;
              if (denom_s == 0.f) {  // range_sVar[0] == range_sVar[1]
                if      (pixel_data[ij] < range_sVar[0]) daTa_uc = (unsigned char)0;
                else if (pixel_data[ij] > range_sVar[0]) daTa_uc = (unsigned char)255;
                else                                     daTa_uc = (unsigned char)127;
              }
              else {
                //Expected by image diagnostics: normalized data 0-1 uses 255 bins of precision, endpoints centered at bins 0, 255
                // daTa_uc = (unsigned char)max(0,min(255,int((pixel_data[ij]-range_sVar[0])/(range_sVar[1]-range_sVar[0])*255.0f+0.5f)));
                float d_tmp = max(0,min(255,int((pixel_data[ij]-range_sVar[0])*denom_s*255.0f+0.5f)));

                if (n_range_sVar) {
                  d_tmp = float(int(d_tmp*f_sbin[1]))*f_sbin[0];
                }
                daTa_uc = (unsigned char)d_tmp;
              }
              assert(daTa);
              daTa->set(ij,&daTa_uc);
              if (rgb_mode == DATA) {
                float sign = -1.0f;
                if (pixel_normal[ij][2] >= 0.f) sign = 1.0f;
                const float a = sign*(pixel_normal[ij][0]*light1[0] + pixel_normal[ij][1]*light1[1] + pixel_normal[ij][2]*light1[2]);
                const float b = sign*(pixel_normal[ij][0]*light2[0] + pixel_normal[ij][1]*light2[1] + pixel_normal[ij][2]*light2[2]);
                const float p = (a*base2[0] + (1.0-a)*base00[0])*0.5 + (b*base3[0] + (1.0-b)*base00[0])*0.5;
                unsigned char lgFr_uc = (unsigned char)max(0,min(255,int(p*255.0f+0.5f)));
                lgFr->set(ij,&lgFr_uc);
                //rgb[ij][0] = daTa_uc;
                //apply lighting to data in 256 bins...
                rgb[ij][0] = max(0,min(255,(int)(p*daTa_uc)));
                rgb[ij][1] = rgb[ij][0];
                rgb[ij][2] = rgb[ij][0];
              }
              //else if (rgb_mode == NORMALS) {
              //  assert(daTa);
              //  daTa->set(ij,&daTa_uc);
              //}
              //else {
              //  assert(0);
              //}
              break;
            }
          default:
            assert(0);
          }
        }
      }
    }
    else {

      // ==================================================
      // b_volvis...
      // ==================================================

      for (int ii = 0; ii < volvis_size; ++ii) {
        if (volvis_data[ii][1] != 0.f) // TODO look into this depth = 0.f business
          min_depth = min(min_depth,volvis_data[ii][1]);
      }

      // add the surface if requested...
      if (b_volvis_surface) {
        for (int ij = 0; ij < ni*nj; ++ij) {
          if (pixel_flag[ij] == SURFACE_PIXEL) {
            // set the zone...
            uint2 zoNe_us = pixel_zone[ij];
            zoNe_uc[0] = zoNe_us & 255;
            zoNe_uc[1] = zoNe_us >> 8;
            zoNe->set(ij,zoNe_uc);

            // set the depth...
            uint2 dpTh_us = (uint2) max(0.f,min(65534.f,pixel_sdepth[ij]-min_depth));
            dpTh_uc[0] = dpTh_us & 255;
            dpTh_uc[1] = dpTh_us >> 8;
            dpTh->set(ij,dpTh_uc);

            // set rgb to normal in NORMALS rgb_mode...
            if (rgb_mode == NORMALS) {
              rgb[ij][0] = (unsigned char) max(0,min(255,int((pixel_normal[ij][0]+1.0f)*128.0f)));
              rgb[ij][1] = (unsigned char) max(0,min(255,int((pixel_normal[ij][1]+1.0f)*128.0f)));
              rgb[ij][2] = (unsigned char) max(0,min(255,int((pixel_normal[ij][2]+1.0f)*128.0f)));
            }
            else if (rgb_mode == DATA) {
              // silver...
              if (pixel_normal[ij][2] >= 0.f) {
                const float a = pixel_normal[ij][0]*light1[0] + pixel_normal[ij][1]*light1[1] + pixel_normal[ij][2]*light1[2];
                const float b = pixel_normal[ij][0]*light2[0] + pixel_normal[ij][1]*light2[1] + pixel_normal[ij][2]*light2[2];
                const float p = (a*base2[0] + (1.0-a)*base00[0])*0.5 + (b*base3[0] + (1.0-b)*base00[0])*0.5;
                rgb[ij][0] = max(0,min(255,(int)(p*256.0)));
                rgb[ij][1] = rgb[ij][0];
                rgb[ij][2] = rgb[ij][0];
              }
              // gold...
              else {
                const float a = - pixel_normal[ij][0]*light1[0] - pixel_normal[ij][1]*light1[1] - pixel_normal[ij][2]*light1[2];
                const float b = - pixel_normal[ij][0]*light2[0] - pixel_normal[ij][1]*light2[1] - pixel_normal[ij][2]*light2[2];
                const float p = (a*base2[0] + (1.0-a)*base00[0])*0.5 + (b*base3[0] + (1.0-b)*base00[0])*0.5;
                rgb[ij][0] = max(0,min(255,(int)(p*256.0)));
                rgb[ij][1] = max(0,min(255,(int)(0.93*p*256.0)));
                rgb[ij][2] = max(0,min(255,(int)(0.72*p*256.0)));
              }
            }
            else {
              assert(0);
            }
          }
        }
      }

      if (!b_range) {
        range[0]    =  HUGE_VALF;
        range[1]    = -HUGE_VALF;
        for (int ii = 0; ii < volvis_size; ++ii) {
          if (volvis_data[ii][1] != 0.f) { // what is this?
            range[0] = min(range[0],volvis_data[ii][0]);
            range[1] = max(range[1],volvis_data[ii][0]);
          }
        }
      }

      for (int ii = 0; ii < volvis_size; ++ii) {
        const int ij = volvis_ij[ii];
        if (volvis_data[ii][1] != 0.f) { // and this? (see above)

          // depth
          const uint2 my_depth = (uint2)max(0.f,min(65534.f,volvis_data[ii][1]-min_depth));
          dpTh_uc[0] = (uint2)my_depth & 255;
          dpTh_uc[1] = (uint2)my_depth >> 8;
          dpTh->set(ij,dpTh_uc);

          // zone information
          zoNe_uc[0] = (uint2)65534 & 255;  // HACK write as plane data zone (so no normal coloring)
          zoNe_uc[1] = (uint2)65534 >> 8;
          zoNe->set(ij,zoNe_uc);

          // handle the case of range[0] == range[1] robustly...
          unsigned char daTa_uc;
          if (range[0] == range[1]) {
            if      (volvis_data[ii][0] < range[0]) daTa_uc = (unsigned char)0;
            else if (volvis_data[ii][0] > range[0]) daTa_uc = (unsigned char)255;
            else                                    daTa_uc = (unsigned char)127;
          }
          else {
            //Expected by image diagnostics: normalized data 0-1 uses 255 bins of precision, endpoints centered at bins 0, 255
            daTa_uc = (unsigned char)max(0,min(255,int((volvis_data[ii][0]-range[0])/(range[1]-range[0])*255.0f+0.5f)));
          }

          //TODO any distinction between RGB modes? currently don't setup nrMl
          FOR_I3 rgb[ij][i] = daTa_uc;
          if (daTa) {
            daTa->set(ij,&daTa_uc);
            unsigned char lgFr_uc = (unsigned char)255;
            lgFr->set(ij,&lgFr_uc);
          }

        }
      }
    }

    delete[] pixel_flag;
    DELETE(pixel_data);
    DELETE(pixel_normal);
    DELETE(pixel_zone);
    DELETE(pixel_sdepth);
    DELETE(pixel_mdepth);
    DELETE(pixel_pdepth);
    DELETE(pixel_paux);

    PngImage png(ni,nj,rgb);

    // ===========================
    // metadata...
    // ===========================

    if (imd != NULL) {

      if (b_volumeData||b_volvis){
        imd->setColorMapName("GRAYSCALE_RGB");
        imd->setRangeMin(range[0]);
        imd->setRangeMax(range[1]);
      }
      if (b_surfaceData){
        imd->setColorMapSurfaceName("GRAYSCALE_RGB");
        imd->setRangeOnSurfaceMin(range_sVar[0]);
        imd->setRangeOnSurfaceMax(range_sVar[1]);
      }
      if (b_particleData){
        imd->setColorMapParticlesName("GRAYSCALE_RGB");
        imd->setRangeOnParticleMin(range_pVar[0]);
        imd->setRangeOnParticleMax(range_pVar[1]);
      }
      if (b_isoData){
        imd->setColorMapIsoName("GRAYSCALE_RGB");
        imd->setRangeOnIsoMin(range_iVar[0]);
        imd->setRangeOnIsoMax(range_iVar[1]);
      }

      // light planar mesh data in app during mesh viz
      if (b_delayData)
        imd->setLightData(true);

      const double width_ni_ratio = width/(double)ni;

      imd->setLengthScale(width_ni_ratio);

      imd->setRgbMode(rgb_mode);

      //Map pixel coords to simulation coords
      imd->setCamDepth(min_depth);

      float origin[3] = {0.0,0.0,0.0};
      float p0[3];  //simulation origin in pixel coords
      p0[0] = CALC_PIXEL_I(origin);
      p0[1] = CALC_PIXEL_J(origin);
      p0[2] = CALC_PIXEL_DEPTH(origin) - min_depth; //Adjust for depth buffer size limit

      //simulation axes in pixel coords
      float eX[3] = {e0[0],e1[0],e2[0]};
      float eX_mag = MAG(eX);
      float eY[3] = {e0[1],e1[1],e2[1]};
      float eY_mag = MAG(eY);
      float eZ[3] = {e0[2],e1[2],e2[2]};
      float eZ_mag = MAG(eZ);
      FOR_K3{
        eX[k] /= eX_mag;
        eY[k] /= eY_mag;
        eZ[k] /= eZ_mag;
      }

      imd->transformMat[0]  = eX[0]*width_ni_ratio;
      imd->transformMat[1]  = eY[0]*width_ni_ratio;
      imd->transformMat[2]  = eZ[0]*width_ni_ratio;
      imd->transformMat[3]  = 0.0;
      imd->transformMat[4]  = eX[1]*width_ni_ratio;
      imd->transformMat[5]  = eY[1]*width_ni_ratio;
      imd->transformMat[6]  = eZ[1]*width_ni_ratio;
      imd->transformMat[7]  = 0.0;
      imd->transformMat[8]  = eX[2]*width_ni_ratio;
      imd->transformMat[9]  = eY[2]*width_ni_ratio;
      imd->transformMat[10] = eZ[2]*width_ni_ratio;
      imd->transformMat[11] = 0.0;
      imd->transformMat[12] = (-p0[0]*eX[0]-p0[1]*eX[1]-p0[2]*eX[2])*width_ni_ratio;
      imd->transformMat[13] = (-p0[0]*eY[0]-p0[1]*eY[1]-p0[2]*eY[2])*width_ni_ratio;
      imd->transformMat[14] = (-p0[0]*eZ[0]-p0[1]*eZ[1]-p0[2]*eZ[2])*width_ni_ratio;
      imd->transformMat[15] = 1.0;

    }

    png.setMetadata(imd); // DP: why not have this inside the above NULL check?

    dpTh->compressChunk();
    png.addPngDataChunk(dpTh);

    zoNe->compressChunk();
    png.addPngDataChunk(zoNe);

    if (daTa){
      daTa->compressChunk();
      png.addPngDataChunk(daTa);
    }

    if (lgFr){
      lgFr->compressChunk();
      png.addPngDataChunk(lgFr);
    }

    MiscUtils::mkdir_for_file(filename);
    string tmp_filename = MiscUtils::makeTmpPrefix(filename);
    png.write(tmp_filename.c_str());
    remove(filename.c_str());
    rename(tmp_filename.c_str(),filename.c_str());

    //cleanup buffers, no longer handled by PngImage...
    delete[] rgb; rgb = NULL;
    delete dpTh; dpTh = NULL;
    delete zoNe; zoNe = NULL;
    if (daTa != NULL) {
      delete daTa;
      daTa = NULL;
    }
    if (lgFr != NULL) {
      delete lgFr;
      lgFr = NULL;
    }
      
  } // if (mpi_rank == write_rank)...

}

void CtiCanvas::clearPBuf(){
  pbuf_size=0;
  if (pbuf_float){
    delete[] pbuf_float;
    pbuf_float = NULL;
  }
  if (pbuf_uint2){
    delete[] pbuf_uint2;
    pbuf_uint2 = NULL;
  }
}

// get/set pixel data routines...

void CtiCanvas::getIJij(int& IJ,int& ij,const int i, const int j) {

  assert((i >= 0)&&(i < ni));
  assert((j >= 0)&&(j < nj));
  IJ = (j>>3)*nI + (i>>3);
  assert((IJ >= 0)&&(IJ < nI*nJ));
  ij = (j&7)*8+(i&7);
  assert((ij >= 0)&&(ij < 64));

}

bool CtiCanvas::getPixelPdepthBits(float& pdepth,uint2& bits,const int IJ,const int ij) const {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  if (pbd[IJ] == NULL) return false;
  assert((ij >= 0)&&(ij < 64));

  pdepth = pbd[IJ]->pdepth[ij];
  bits = pbd[IJ]->bits[ij];
  return true;

}

bool CtiCanvas::getPixelPauxBits(float& paux,uint2& bits,const int IJ,const int ij) const {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  if (pbd[IJ] == NULL) return false;
  assert((ij >= 0)&&(ij < 64));

  paux = pbd[IJ]->paux[ij];
  bits = pbd[IJ]->bits[ij];
  return true;

}

bool CtiCanvas::getPixelPauxPzoneBits(float& paux,uint2& pzone,uint2& bits,const int IJ,const int ij) const {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  if (pbd[IJ] == NULL) return false;
  assert((ij >= 0)&&(ij < 64));

  paux = pbd[IJ]->paux[ij];
  pzone = pbd[IJ]->pzone[ij];
  bits = pbd[IJ]->bits[ij];
  return true;

}

bool CtiCanvas::getPixelSdepthBits(float& sdepth,uint2& bits,const int IJ,const int ij) const {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  if (pbd[IJ] == NULL) return false;
  assert((ij >= 0)&&(ij < 64));

  sdepth = pbd[IJ]->sdepth[ij];
  bits = pbd[IJ]->bits[ij];
  return true;

}

bool CtiCanvas::getPixelSdepthPauxBits(float& sdepth,float& paux,uint2& bits,const int IJ,const int ij) const {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  if (pbd[IJ] == NULL) return false;
  assert((ij >= 0)&&(ij < 64));

  sdepth = pbd[IJ]->sdepth[ij];
  paux = pbd[IJ]->paux[ij];
  bits = pbd[IJ]->bits[ij];
  return true;

}

bool CtiCanvas::getPixelSdepthPauxSzoneBits(float& sdepth,float& paux,uint2& szone,uint2& bits,const int IJ,const int ij) const {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  if (pbd[IJ] == NULL) return false;
  assert((ij >= 0)&&(ij < 64));

  sdepth = pbd[IJ]->sdepth[ij];
  paux = pbd[IJ]->paux[ij];
  szone = pbd[IJ]->szone[ij];
  bits = pbd[IJ]->bits[ij];
  return true;

}

bool CtiCanvas::getPixelSdepthMdepthBits(float& sdepth,float& mdepth,uint2& bits,const int IJ,const int ij) const {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  if (pbd[IJ] == NULL) return false;
  assert((ij >= 0)&&(ij < 64));

  sdepth = pbd[IJ]->sdepth[ij];
  mdepth = pbd[IJ]->mdepth[ij];
  bits = pbd[IJ]->bits[ij];
  return true;

}

void CtiCanvas::setPixelSnormalSdepthSzone(const int IJ,const int ij,const float snormal[3],const float sdepth,const uint2 szone) {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  FOR_K3 pbd[IJ]->snormal[ij][k] = snormal[k];
  pbd[IJ]->sdepth[ij] = sdepth;
  pbd[IJ]->szone[ij] = szone;

}

void CtiCanvas::setPixelSdataBits(const int IJ,const int ij,const float sdata,const uint2 bits) {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  pbd[IJ]->sdata[ij] = sdata;
  pbd[IJ]->bits[ij] = bits;

}

void CtiCanvas::setPixelSnormalSdepthSdataSzone(const int IJ,const int ij,const float snormal[3],const float sdepth,const float sdata,const uint2 szone) {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  FOR_K3 pbd[IJ]->snormal[ij][k] = snormal[k];
  pbd[IJ]->sdepth[ij] = sdepth;
  pbd[IJ]->sdata[ij] = sdata;
  pbd[IJ]->szone[ij] = szone;

}

void CtiCanvas::setPixelSnormalSdepthSdataSzoneBits(const int IJ,const int ij,const float snormal[3],const float sdepth,const float sdata,const uint2 szone,const uint2 bits) {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  FOR_K3 pbd[IJ]->snormal[ij][k] = snormal[k];
  pbd[IJ]->sdepth[ij] = sdepth;
  pbd[IJ]->sdata[ij] = sdata;
  pbd[IJ]->szone[ij] = szone;
  pbd[IJ]->bits[ij] = bits;

}

void CtiCanvas::setPixelSnormalSdepthSzoneBits(const int IJ,const int ij,const float snormal[3],const float sdepth,const uint2 szone,const uint2 bits) {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  FOR_K3 pbd[IJ]->snormal[ij][k] = snormal[k];
  pbd[IJ]->sdepth[ij] = sdepth;
  pbd[IJ]->szone[ij] = szone;
  pbd[IJ]->bits[ij] = bits;

}

void CtiCanvas::setPixelPnormalPdepthPdataPaux(const int IJ,const int ij,const float pnormal[3],const float pdepth,const float pdata,const float paux) {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  FOR_K3 pbd[IJ]->pnormal[ij][k] = pnormal[k];
  pbd[IJ]->pdepth[ij] = pdepth;
  pbd[IJ]->pdata[ij] = pdata;
  pbd[IJ]->paux[ij] = paux;

}

void CtiCanvas::setPixelPnormalPdepthPdataPauxBits(const int IJ,const int ij,const float pnormal[3],const float pdepth,const float pdata,const float paux,const uint2 bits) {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  FOR_K3 pbd[IJ]->pnormal[ij][k] = pnormal[k];
  pbd[IJ]->pdepth[ij] = pdepth;
  pbd[IJ]->pdata[ij] = pdata;
  pbd[IJ]->paux[ij] = paux;
  pbd[IJ]->bits[ij] = bits;

}

void CtiCanvas::setPixelPnormalPdepthPdataPauxPzone(const int IJ,const int ij,const float pnormal[3],const float pdepth,const float pdata,const float paux,const uint2 pzone) {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  FOR_K3 pbd[IJ]->pnormal[ij][k] = pnormal[k];
  pbd[IJ]->pdepth[ij] = pdepth;
  pbd[IJ]->pdata[ij] = pdata;
  pbd[IJ]->paux[ij] = paux;
  pbd[IJ]->pzone[ij] = pzone;

}

void CtiCanvas::setPixelPnormalPdepthPdataPauxPzoneBits(const int IJ,const int ij,const float pnormal[3],const float pdepth,const float pdata,const float paux,const uint2 pzone,const uint2 bits) {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  FOR_K3 pbd[IJ]->pnormal[ij][k] = pnormal[k];
  pbd[IJ]->pdepth[ij] = pdepth;
  pbd[IJ]->pdata[ij] = pdata;
  pbd[IJ]->paux[ij] = paux;
  pbd[IJ]->pzone[ij] = pzone;
  pbd[IJ]->bits[ij] = bits;

}

void CtiCanvas::setPixelPnormalPdepthPzone(const int IJ,const int ij,const float pnormal[3],const float pdepth,const uint2 pzone) {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  FOR_K3 pbd[IJ]->pnormal[ij][k] = pnormal[k];
  pbd[IJ]->pdepth[ij] = pdepth;
  pbd[IJ]->pzone[ij] = pzone;

}

void CtiCanvas::setPixelPnormalPdepthPzoneBits(const int IJ,const int ij,const float pnormal[3],const float pdepth,const uint2 pzone,const uint2 bits) {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  FOR_K3 pbd[IJ]->pnormal[ij][k] = pnormal[k];
  pbd[IJ]->pdepth[ij] = pdepth;
  pbd[IJ]->pzone[ij] = pzone;
  pbd[IJ]->bits[ij] = bits;

}

void CtiCanvas::setPixelPnormalPdepthPdataPzone(const int IJ,const int ij,const float pnormal[3],const float pdepth,const float pdata,const uint2 pzone) {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  FOR_K3 pbd[IJ]->pnormal[ij][k] = pnormal[k];
  pbd[IJ]->pdepth[ij] = pdepth;
  pbd[IJ]->pdata[ij] = pdata;
  pbd[IJ]->pzone[ij] = pzone;

}

void CtiCanvas::setPixelPnormalPdepthPdataPzoneBits(const int IJ,const int ij,const float pnormal[3],const float pdepth,const float pdata,const uint2 pzone,const uint2 bits) {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  FOR_K3 pbd[IJ]->pnormal[ij][k] = pnormal[k];
  pbd[IJ]->pdepth[ij] = pdepth;
  pbd[IJ]->pdata[ij] = pdata;
  pbd[IJ]->pzone[ij] = pzone;
  pbd[IJ]->bits[ij] = bits;

}

void CtiCanvas::setPixelMdepth(const int IJ,const int ij,const float mdepth) {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  pbd[IJ]->mdepth[ij] = mdepth;

}

void CtiCanvas::setPixelMdepthBits(const int IJ,const int ij,const float mdepth,const uint2 bits) {

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  pbd[IJ]->mdepth[ij] = mdepth;
  pbd[IJ]->bits[ij] = bits;

}

void CtiCanvas::setPixelSdataPaux(const int IJ,const int ij,const float sdata,const float paux) {

  assert(paux == paux);

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  pbd[IJ]->sdata[ij] = sdata;
  pbd[IJ]->paux[ij] = paux;

}

void CtiCanvas::setPixelSdataPauxBits(const int IJ,const int ij,const float sdata,const float paux,const uint2 bits) {

  assert(paux == paux);

  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);
  assert((ij >= 0)&&(ij < 64));

  pbd[IJ]->sdata[ij] = sdata;
  pbd[IJ]->paux[ij] = paux;
  pbd[IJ]->bits[ij] = bits;

}

//Return true if mdepth is less than zero
bool CtiCanvas::isMaskPixel(const int &image_ij) const {
  //   cout << "image_ij " << image_ij << " " << ni*nj;
  assert(image_ij<ni*nj);

  //image coordinates
  const int image_j = (int) image_ij/ni;
  const int image_i = image_ij%ni;

  // determine the index of the PixelBlock...
  const int IJ = (image_j>>3)*nI + (image_i>>3);
  assert((IJ >= 0)&&(IJ < nI*nJ));

  if (pbd[IJ] == NULL)
    return false;

  const int ij = (image_j&7)*8+(image_i&7);
  assert((ij >= 0)&&(ij < 64));
  return (pbd[IJ]->mdepth[ij] < 0.0f);
}

bool CtiCanvas::isMaskPixel(const int IJ, const int ij) const {

  if (pbd[IJ] == NULL)
    return false;
  else
    return (pbd[IJ]->mdepth[ij] < 0.0f);
}

bool CtiCanvas::isSurfacePixel(const int &image_ij) const {
  //  cout << "image_ij " << image_ij;
  assert(image_ij<ni*nj);

  //image coordinates
  const int image_j = (int) image_ij/ni;
  const int image_i = image_ij%ni;

  // determine the index of the PixelBlock...
  const int IJ = (image_j>>3)*nI + (image_i>>3);
  assert((IJ >= 0)&&(IJ < nI*nJ));

  if (pbd[IJ] == NULL)
    return false;

  const int ij = (image_j&7)*8+(image_i&7);
  assert((ij >= 0)&&(ij < 64));
  //return (pbd[IJ]->szone[ij] < 65534);
  return (pbd[IJ]->bits[ij]&SURFACE_BITS);

}

void CtiCanvas::addCimPlaneData(const CvImageMap &cim, double * var) {

  //int count=0;
  for (int icv = 0; icv < cim.ncv_unique; ++icv) {
    for (int poc = cim.poucv_i[icv]; poc != cim.poucv_i[icv+1]; ++poc) {
      unsigned int index = cim.getIndex(poc);
      unsigned int repeat_count = cim.getRepeatCount(poc);
      assert( (index | (repeat_count<<(32-CIM_REPEAT_BITS))) == cim.poucv_v[poc]);
      for (unsigned int offset = 0; offset <= repeat_count; ++offset) {
        const int image_ij = index+offset;
        const int image_j = (int) image_ij/ni;
        const int image_i = image_ij%ni;
        float planedepth = dataPlaneData->center[2] - (dataPlaneData->normal[0]*(image_i-dataPlaneData->center[0]) +
						       dataPlaneData->normal[1]*(image_j-dataPlaneData->center[1]) )
          / dataPlaneData->normal[2]; //compute depth of data plane

        if (!isBlankedSkipDataPlane(image_i,image_j,planedepth,true)) {
          const int IJ = (image_j>>3)*nI + (image_i>>3);
          const int ij = (image_j&7)*8+(image_i&7);

          setPixelPnormalPdepthPdataPzoneBits(IJ,ij,dataPlaneData->normal,planedepth,var[icv],65534,pbd[IJ]->bits[ij]|INTERNAL_DATA_BIT|MASK_BIT);
        }
        //++count;
      }
    }
  }
  //cout << "CtiCanvas::fillDataPlane() count=" << dap_count << endl;
}

//used for snapshot imaging, called after CacheImage(), pdb==NULL
void CtiCanvas::updateCimPlaneData(const CvImageMap &cim, double * var) {

  //cim maps cvs to image pixels, need map from pixels to PBUF ip...
  //this will only work in serial surfer
  int * ipoij = new int[ni*nj];
  for (int ij = 0; ij < ni*nj; ++ij)
    ipoij[ij] = -1;
  for (int ip = 0; ip < pbuf_size; ++ip) {
    const int ij = PBUF_J*ni + PBUF_I;
    assert(ipoij[ij]==-1); //this only works in serial when ip to ij relation is 1 to 1
    ipoij[ij] = ip;
  }

  for (int icv = 0; icv < cim.ncv_unique; ++icv) {
    for (int poc = cim.poucv_i[icv]; poc != cim.poucv_i[icv+1]; ++poc) {
      unsigned int index = cim.getIndex(poc);
      unsigned int repeat_count = cim.getRepeatCount(poc);
      assert( (index | (repeat_count<<(32-CIM_REPEAT_BITS))) == cim.poucv_v[poc]);
      for (unsigned int offset = 0; offset <= repeat_count; ++offset) {
        const int image_ij = index+offset;

        int ip = ipoij[image_ij];
        if (ip>-1){
          const int image_j = (int) image_ij/ni;
          const int image_i = image_ij%ni;
          float planedepth = dataPlaneData->center[2] - (dataPlaneData->normal[0]*(image_i-dataPlaneData->center[0]) +
                                                         dataPlaneData->normal[1]*(image_j-dataPlaneData->center[1]) )
            / dataPlaneData->normal[2]; //compute depth of data plane

          if (!isBlankedSkipDataPlane(image_i,image_j,planedepth,true)) {
            //pnormal,pdepth,pzone and bits were set in addCimPlaneData, here only pdata has changed and must
            //be updated...
            PBUF_PDATA = var[icv];
          }
        }
      }
    }
  }
  delete[] ipoij;

}

void CtiCanvas::addCimPlaneMesh(const CvImageMap &cim, float * d2_buf) {

  if (cim.ncv_unique > 0) {

    double nijk[3]; FOR_I3 nijk[i] = dataPlaneData->normal[i];
    const double nijk_mag = MAG(nijk); assert(nijk_mag > 0.0);
    const double unit_nijk[3] = {nijk[0]/nijk_mag,nijk[1]/nijk_mag,nijk[2]/nijk_mag};
    const float normal[3] = {(float)unit_nijk[0],(float)unit_nijk[1],(float)unit_nijk[2]};

    for (int icv = 0; icv < cim.ncv_unique; ++icv) {
      for (int poc = cim.poucv_i[icv]; poc != cim.poucv_i[icv+1]; ++poc) {
        unsigned int index = cim.getIndex(poc);
        unsigned int repeat_count = cim.getRepeatCount(poc);
        assert( (index | (repeat_count<<(32-CIM_REPEAT_BITS))) == cim.poucv_v[poc]);
        for (unsigned int offset = 0; offset <= repeat_count; ++offset) {
          const int image_ij = index+offset;
          const int image_j = (int) image_ij/ni;
          const int image_i = image_ij%ni;

          float planedepth = dataPlaneData->center[2] - (dataPlaneData->normal[0]*(image_i-dataPlaneData->center[0]) +
              dataPlaneData->normal[1]*(image_j-dataPlaneData->center[1]) )/dataPlaneData->normal[2]; //compute depth of data plane

          if (!isBlankedSkipDataPlane(image_i,image_j,planedepth,true)) {
            const int IJ = (image_j>>3)*nI + (image_i>>3);
            const int ij = (image_j&7)*8+(image_i&7);
            setPixelPnormalPdepthPdataPzoneBits(IJ,ij,normal,planedepth,(float)icv,65534-32768,pbd[IJ]->bits[ij]|INTERNAL_DATA_BIT|MESH_BIT|AUX_BIT|MASK_BIT|DELAY_BIT);
          }
        }
      }
    }

    /*
    //normalize the range of d2_buf for each unique icv
    float * d2_min = new float[cim.ncv_unique];
    float * d2_max = new float[cim.ncv_unique];

    for (int icv = 0; icv < cim.ncv_unique; ++icv) {
      d2_min[icv] =  HUGE_VALF;
      d2_max[icv] = -HUGE_VALF;
      //const int8 icv_global = cim.ucv_v[icv];
      for (int poc = cim.poucv_i[icv]; poc != cim.poucv_i[icv+1]; ++poc) {
        unsigned int index = cim.poucv_v[poc] & ((1<<(32-CIM_REPEAT_BITS))-1);
        unsigned int repeat_count = cim.poucv_v[poc]>>(32-CIM_REPEAT_BITS);
        assert( (index | (repeat_count<<(32-CIM_REPEAT_BITS))) == cim.poucv_v[poc]);
        for (unsigned int offset = 0; offset <= repeat_count; ++offset) {
          if (d2_buf[index+offset]<d2_min[icv])
            d2_min[icv] = d2_buf[index+offset];
          if (d2_buf[index+offset]>d2_max[icv])
            d2_max[icv] = d2_buf[index+offset];
        }
      }
    }
    //MiscUtils::dumpRange(d2_min,cim.ncv_unique,"d2_min");
    //MiscUtils::dumpRange(d2_max,cim.ncv_unique,"d2_max");
    for (int icv = 0; icv < cim.ncv_unique; ++icv) {
      for (int poc = cim.poucv_i[icv]; poc != cim.poucv_i[icv+1]; ++poc) {
        unsigned int index = cim.getIndex(poc);
        unsigned int repeat_count = cim.getRepeatCount(poc);
        assert( (index | (repeat_count<<(32-CIM_REPEAT_BITS))) == cim.poucv_v[poc]);
        for (unsigned int offset = 0; offset <= repeat_count; ++offset) {
          const int image_ij = index+offset;
          const int image_j = (int) image_ij/ni;
          const int image_i = image_ij%ni;

          float planedepth = dataPlaneData->center[2] - (dataPlaneData->normal[0]*(image_i-dataPlaneData->center[0]) +
                                                         dataPlaneData->normal[1]*(image_j-dataPlaneData->center[1]) )
            /dataPlaneData->normal[2]; //compute depth of data plane

          if (!isBlankedSkipDataPlane(image_i,image_j,planedepth,true)) {
            float d2_norm=1.0f;
            if (d2_max[icv]>d2_min[icv]){
              d2_norm = 1.0-(d2_buf[index+offset]-d2_min[icv])/(d2_max[icv]-d2_min[icv]);
            }

            const int IJ = (image_j>>3)*nI + (image_i>>3);
            const int ij = (image_j&7)*8+(image_i&7);
            setPixelPnormalPdepthPdataPzoneBits(IJ,ij,dataPlaneData->normal,planedepth,d2_norm,65534,pbd[IJ]->bits[ij]|INTERNAL_DATA_BIT|MASK_BIT);
          }
        }
      }
    }
    delete[] d2_min;
    delete[] d2_max;
    */
  }

}

void CtiCanvas::addCimSurfaceData(const CvImageMap &cim, double * var) {

  for (int icv = 0; icv < cim.ncv_unique; ++icv) {
    for (int poc = cim.poucv_i[icv]; poc != cim.poucv_i[icv+1]; ++poc) {
      unsigned int index = cim.getIndex(poc);
      unsigned int repeat_count = cim.getRepeatCount(poc);
      assert( (index | (repeat_count<<(32-CIM_REPEAT_BITS))) == cim.poucv_v[poc]);
      for (unsigned int offset = 0; offset <= repeat_count; ++offset) {
        const int image_ij = index+offset;
        const int image_j = (int) image_ij/ni;
        const int image_i = image_ij%ni;

        assert((image_i >= 0)&&(image_i < ni));
        assert((image_j >= 0)&&(image_j < nj));
        // determine the index of the PixelBlock...
        const int IJ = (image_j>>3)*nI + (image_i>>3);
        const int ij = (image_j&7)*8+(image_i&7);

        setPixelSdataBits(IJ,ij,var[icv],(pbd[IJ]->bits[ij]&~SURFACE_BIT)|SURFACE_DATA_BIT);
      }
    }
  }

  // TODO: confirm no longer needed with new imaging structure...
  // added this for periodicity (sometimes there was no data available)
  //setDataFreePartsOfSurfaceToBackgroundPbd();
}

//used for snapshot imaging, called after cacheImage(), pdb==NULL
void CtiCanvas::updateCimSurfaceData(const CvImageMap &cim, double * var) {

  //cim maps cvs to image pixels, need map from pixels to PBUF ip...
  //this will only work in serial surfer
  int * ipoij = new int[ni*nj];
  for (int ij = 0; ij < ni*nj; ++ij)
    ipoij[ij] = -1;
  for (int ip = 0; ip < pbuf_size; ++ip) {
    const int ij = PBUF_J*ni + PBUF_I;
    assert(ipoij[ij]==-1); //this only works in serial when ip to ij relation is 1 to 1
    ipoij[ij] = ip;
  }

  for (int icv = 0; icv < cim.ncv_unique; ++icv) {
    for (int poc = cim.poucv_i[icv]; poc != cim.poucv_i[icv+1]; ++poc) {
      unsigned int index = cim.getIndex(poc);
      unsigned int repeat_count = cim.getRepeatCount(poc);
      assert( (index | (repeat_count<<(32-CIM_REPEAT_BITS))) == cim.poucv_v[poc]);
      for (unsigned int offset = 0; offset <= repeat_count; ++offset) {
        const int image_ij = index+offset;

        int ip = ipoij[image_ij];

        //the surface cim should only map to valid indices...a new surface cim
        //is generated when the zone hiding or blanking configuration changes...
        assert(ip>-1);

        //surface data bit set in addCimSurfaceData, here we
        //just update the data
        PBUF_SDATA = var[icv];
      }
    }
  }

  delete[] ipoij;

  // TODO: confirm no longer needed with new imaging structure...
  // added this for periodicity (sometimes there was no data available)
  //setDataFreePartsOfSurfaceToBackgroundPbuf();
}

//speedup serial surface data visualization by reading only zones
//present on the canvas
void CtiCanvas::getActiveZones(IntFlag &bfzone_flag,const IntFlag &szozn_i){

  std::set<int> activeZones;
  for (int image_j = 0; image_j < nj; ++image_j) {
    for (int image_i = 0; image_i < ni; ++image_i) {
      // determine the index of the PixelBlock...
      const int IJ = (image_j>>3)*nI + (image_i>>3);
      assert((IJ >= 0)&&(IJ < nI*nJ));
      if (pbd[IJ] != NULL){
        const int ij = (image_j&7)*8+(image_i&7);
        assert((ij >= 0)&&(ij < 64));
        if (pbd[IJ]->szone[ij] < 32768)
          activeZones.insert(pbd[IJ]->szone[ij]);
      }
    }
  }
  //cout << "Found subzones: " << endl;
  //for (std::set<int>::iterator it = activeZones.begin(); it!=activeZones.end(); ++it){
  //   cout << "  subzone id: " << *it << endl;
  //}

  if (activeZones.size()==0)
    return;

  //flag zones with at least one subzone on the canvas...
  for (int izn = 0; izn < szozn_i.getLength() - 1; ++izn) {
    //cout << " Zone " << izn << " contains subzones: " << endl;
    for (int iszn = szozn_i[izn]; iszn < szozn_i[izn+1]; ++iszn) {
      //cout << iszn << endl;
      if ( activeZones.find(iszn) != activeZones.end() ){
        //cout << " Found Zone in image with id: " << izn << endl;
        bfzone_flag[izn] = 1;
        break; //include a zone if a single subzone is present
      }
    }
  }

}


// GI LIC -------------------------
void CtiCanvas::addCimPlaneLIC(const CvImageMap &cim, float *lic_image) {

  for (int icv = 0; icv < cim.ncv_unique; ++icv) {
    for (int poc = cim.poucv_i[icv]; poc != cim.poucv_i[icv+1]; ++poc) {
      unsigned int index = cim.getIndex(poc);
      unsigned int repeat_count = cim.getRepeatCount(poc);
      assert( (index | (repeat_count<<(32-CIM_REPEAT_BITS))) == cim.poucv_v[poc]);
      for (unsigned int offset = 0; offset <= repeat_count; ++offset) {
        const int image_ij = index+offset;
        const int image_j = (int) image_ij/ni;
        const int image_i = image_ij%ni;
        float planedepth = dataPlaneData->center[2] - (dataPlaneData->normal[0]*(image_i-dataPlaneData->center[0]) +
                                                       dataPlaneData->normal[1]*(image_j-dataPlaneData->center[1]) )
          / dataPlaneData->normal[2]; //compute depth of data plane

        if (!isBlankedSkipDataPlane(image_i,image_j,planedepth,true)) {
          const int IJ = (image_j>>3)*nI + (image_i>>3);
          const int ij = (image_j&7)*8+(image_i&7);
          setPixelPnormalPdepthPdataPzoneBits(IJ,ij,dataPlaneData->normal,planedepth,lic_image[image_j*ni+image_i],65534,pbd[IJ]->bits[ij]|INTERNAL_DATA_BIT|MASK_BIT);
        }
      }
    }
  }
}


void CtiCanvas::addCimSurfaceLIC(const CvImageMap &cim, float *lic_image) {
  for (int icv = 0; icv < cim.ncv_unique; ++icv) {
    for (int poc = cim.poucv_i[icv]; poc != cim.poucv_i[icv+1]; ++poc) {
      unsigned int index = cim.getIndex(poc);
      unsigned int repeat_count = cim.getRepeatCount(poc);
      assert( (index | (repeat_count<<(32-CIM_REPEAT_BITS))) == cim.poucv_v[poc]);
      for (unsigned int offset = 0; offset <= repeat_count; ++offset) {
        const int image_ij = index+offset;
        const int image_j = (int) image_ij/ni;
        const int image_i = image_ij%ni;

        assert((image_i >= 0)&&(image_i < ni));
        assert((image_j >= 0)&&(image_j < nj));
        // determine the index of the PixelBlock...
        const int IJ = (image_j>>3)*nI + (image_i>>3);
        const int ij = (image_j&7)*8+(image_i&7);
        setPixelSdataBits(IJ,ij,lic_image[image_j*ni+image_i],(pbd[IJ]->bits[ij]&~SURFACE_BIT)|SURFACE_DATA_BIT);
        //++count;
      }
    }
  }

  // TODO: confirm no longer needed with new imaging structure...
  // added this for periodicity (sometimes there was no data available)
  //setDataFreePartsOfSurfaceToBackgroundPbd();
}

// eliminate the normal component from a vector
void CtiCanvas::projectVectorOnPlane(const CvImageMap &cim, double (*vol_vec)[3],const float np[3]) {

  // GI: get some info about the canvas
  float x0[3];  //simulation coordinates, center of the view
  float e0[3];  //NOTE, not unit vectors, have length ni/width
  float e1[3];
  float e2[3];
  getX0(x0);
  getE0(e0);
  getE1(e1);
  getE2(e2);

  for (int icv = 0; icv < cim.ncv_unique; ++icv) {
    for (int poc = cim.poucv_i[icv]; poc != cim.poucv_i[icv+1]; ++poc) {
      unsigned int index = cim.getIndex(poc);
      unsigned int repeat_count = cim.getRepeatCount(poc);
      assert( (index | (repeat_count<<(32-CIM_REPEAT_BITS))) == cim.poucv_v[poc]);
      for (unsigned int offset = 0; offset <= repeat_count; ++offset) {
        //const int image_ij = index+offset;
        //const int image_j = (int) image_ij/ni;
        //const int image_i = image_ij%ni;

        double ut[3],myn[3];
        FOR_I3 myn[i] = np[i];
        double vn = DOT_PRODUCT(vol_vec[icv],myn);
        FOR_I3 ut[i] = vol_vec[icv][i]-vn*myn[i];
        FOR_I3 vol_vec[icv][i] = ut[i];
      }
    }
  }

}

void CtiCanvas::projectVectorOnSphere(const CvImageMap &cim, double (*vol_vec)[3], double (*vol_xcc)[3]) {
  for (int icv = 0; icv < cim.ncv_unique; ++icv) {
    float myn[3],vn,ut[3],rad;
    rad = 0.0;
    FOR_I3 rad += vol_xcc[icv][i]*vol_xcc[icv][i];
    rad = sqrt(rad);
    FOR_I3 myn[i] = vol_xcc[icv][i]/rad;
    vn = DOT_PRODUCT(vol_vec[icv],myn);
    FOR_I3 ut[i] = vol_vec[icv][i]-vn*myn[i];
    FOR_I3 vol_vec[icv][i] = ut[i];
  }
}

void CtiCanvas::projectVectorOnSurfaceAveraged(const CvImageMap &cim, double (*vol_vec)[3]) {
  for (int icv = 0; icv < cim.ncv_unique; ++icv) {
    int icv_found=0;
    float this_myn[3];
    FOR_I3 this_myn[i] = 0.0;
    for (int poc = cim.poucv_i[icv]; poc != cim.poucv_i[icv+1]; ++poc) {
      unsigned int index = cim.getIndex(poc);
      unsigned int repeat_count = cim.getRepeatCount(poc);
      assert( (index | (repeat_count<<(32-CIM_REPEAT_BITS))) == cim.poucv_v[poc]);
      for (unsigned int offset = 0; offset <= repeat_count; ++offset) {
        const int image_ij = index+offset;
        const int image_j = (int) image_ij/ni;
        const int image_i = image_ij%ni;
        assert((image_i >= 0)&&(image_i < ni));
        assert((image_j >= 0)&&(image_j < nj));
        // determine the index of the PixelBlock...
        const int IJ = (image_j>>3)*nI + (image_i>>3);
        assert((IJ >= 0)&&(IJ < nI*nJ));
        //if blanking is enabled, not all pixels in cim will
        //have been set, so must check pdb is not null and
        //that the zone id is < 65534
        if (pbd[IJ] != NULL){
          const int ij = (image_j&7)*8+(image_i&7);
          assert((ij >= 0)&&(ij < 64));
          if (pbd[IJ]->szone[ij] < 65534) {
            if (icv_found==0) {
              FOR_I3 this_myn[i] += pbd[IJ]->snormal[ij][i];
              //icv_found=1;
            }
          }
        }
      }
    }
    float this_myn_mag,vn,ut[3];
    this_myn_mag = 0.0;
    FOR_I3 this_myn_mag += this_myn[i]*this_myn[i];
    FOR_I3 this_myn[i] = this_myn[i]/this_myn_mag; // this is NOT right....
    vn = DOT_PRODUCT(vol_vec[icv],this_myn);
    FOR_I3 ut[i] = vol_vec[icv][i]-vn*this_myn[i];
    FOR_I3 vol_vec[icv][i] = ut[i];
  }
}

void CtiCanvas::projectVectorOnSurface(const CvImageMap &cim, double (*vol_vec)[3]) {
  for (int icv = 0; icv < cim.ncv_unique; ++icv) {
    for (int poc = cim.poucv_i[icv]; poc != cim.poucv_i[icv+1]; ++poc) {
      unsigned int index = cim.getIndex(poc);
      unsigned int repeat_count = cim.getRepeatCount(poc);
      assert( (index | (repeat_count<<(32-CIM_REPEAT_BITS))) == cim.poucv_v[poc]);
      for (unsigned int offset = 0; offset <= repeat_count; ++offset) {
        const int image_ij = index+offset;
        const int image_j = (int) image_ij/ni;
        const int image_i = image_ij%ni;

        assert((image_i >= 0)&&(image_i < ni));
        assert((image_j >= 0)&&(image_j < nj));
        // determine the index of the PixelBlock...
        const int IJ = (image_j>>3)*nI + (image_i>>3);
        assert((IJ >= 0)&&(IJ < nI*nJ));
        //if blanking is enabled, not all pixels in cim will
        //have been set, so must check pdb is not null and
        //that the zone id is < 65534
        if (pbd[IJ] != NULL){
          const int ij = (image_j&7)*8+(image_i&7);
          assert((ij >= 0)&&(ij < 64));
          if (pbd[IJ]->szone[ij] < 65534) {
            float myn[3],vn,ut[3];
            FOR_I3 myn[i] = pbd[IJ]->snormal[ij][i];
            vn = DOT_PRODUCT(vol_vec[icv],myn);
            FOR_I3 ut[i] = vol_vec[icv][i]-vn*myn[i];
            FOR_I3 vol_vec[icv][i] = ut[i];
          }
        }
      }
    }
  }
}

void CtiCanvas::addSurfaceTriVolVis(const double xt0[3],const double xt1[3],const double xt2[3],const bool show,const double n0,const double n1,const double n2) {

  // this routine adds start/stop points to the volVisDataVec for use
  // in volume-visualization...

  // here show is a flag indicating whether this surface is on in the current image...

  const float d0 = CALC_PIXEL_DEPTH(xt0);
  const float d1 = CALC_PIXEL_DEPTH(xt1);
  const float d2 = CALC_PIXEL_DEPTH(xt2);

  //if (max(d0,max(d1,d2)) <= 0.f) return; // allow negative depth here for now

  const float i0 = CALC_PIXEL_I(xt0);
  const float i1 = CALC_PIXEL_I(xt1);
  const float i2 = CALC_PIXEL_I(xt2);

  if ( (max(i0,max(i1,i2)) < 0.f) || (min(i0,min(i1,i2)) > float(ni-1)) ) return;

  const float j0 = CALC_PIXEL_J(xt0);
  const float j1 = CALC_PIXEL_J(xt1);
  const float j2 = CALC_PIXEL_J(xt2);

  if ( (max(j0,max(j1,j2)) < 0.f) || (min(j0,min(j1,j2)) > float(nj-1)) ) return;

  // blanked applied later...

  float normal[3];
  if (n0+n1+n2 == 0.0) {
    normal[0] = (j1-j0)*(d2-d0)-(d1-d0)*(j2-j0);
    normal[1] = (d1-d0)*(i2-i0)-(i1-i0)*(d2-d0);
    normal[2] = (i1-i0)*(j2-j0)-(j1-j0)*(i2-i0);
  }
  else {
    normal[0] = float(n0*e0[0] + n1*e0[1] + n2*e0[2]);
    normal[1] = float(n0*e1[0] + n1*e1[1] + n2*e1[2]);
    normal[2] = float(n0*e2[0] + n1*e2[1] + n2*e2[2]);
  }
  const float mag = sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );
  if (mag == 0.f)
    return;

  // we don't need the normal after this, but do need the dir, so compute without normalizing...

  const bool dir = (normal[2] > 0.f);

  // sort the data by y...
  float x0_s[3],x1_s[3],x2_s[3];
  if (j0 <= j1) {
    if (j1 <= j2) {
      x0_s[0] = i0; x0_s[1] = j0; x0_s[2] = d0;
      x1_s[0] = i1; x1_s[1] = j1; x1_s[2] = d1;
      x2_s[0] = i2; x2_s[1] = j2; x2_s[2] = d2;
    }
    else if (j0 <= j2) {
      x0_s[0] = i0; x0_s[1] = j0; x0_s[2] = d0;
      x1_s[0] = i2; x1_s[1] = j2; x1_s[2] = d2;
      x2_s[0] = i1; x2_s[1] = j1; x2_s[2] = d1;
    }
    else {
      x0_s[0] = i2; x0_s[1] = j2; x0_s[2] = d2;
      x1_s[0] = i0; x1_s[1] = j0; x1_s[2] = d0;
      x2_s[0] = i1; x2_s[1] = j1; x2_s[2] = d1;
    }
  }
  else if (j0 <= j2) {
    x0_s[0] = i1; x0_s[1] = j1; x0_s[2] = d1;
    x1_s[0] = i0; x1_s[1] = j0; x1_s[2] = d0;
    x2_s[0] = i2; x2_s[1] = j2; x2_s[2] = d2;
  }
  else if (j1 <= j2) {
    x0_s[0] = i1; x0_s[1] = j1; x0_s[2] = d1;
    x1_s[0] = i2; x1_s[1] = j2; x1_s[2] = d2;
    x2_s[0] = i0; x2_s[1] = j0; x2_s[2] = d0;
  }
  else {
    x0_s[0] = i2; x0_s[1] = j2; x0_s[2] = d2;
    x1_s[0] = i1; x1_s[1] = j1; x1_s[2] = d1;
    x2_s[0] = i0; x2_s[1] = j0; x2_s[2] = d0;
  }

  // check y-sort...

  assert(x0_s[1] <= x1_s[1]);
  assert(x1_s[1] <= x2_s[1]);

  // the sorting has divided the tri into 2 vertical pieces...

  const float pixel_eps = 1.0E-6f; // magic?

  // these can be equal, but we skip this case...

  const float light1[3] = { 0.0f, 0.0f, 1.0f };
  const float light2[3] = { -0.57f, -0.57f, 0.57f };
  const float base3[3] = { 0.99f, 0.96f, 0.89f }; // * color
  const float base2[3] = { 0.92f, 0.91f, 0.83f };
  const float base00[3] = { 0.40f, 0.48f, 0.51f };

  if (x0_s[1] < x1_s[1]) {
    const int j0_plus = max(0,(int)ceil(x0_s[1]-pixel_eps));
    const int j1_minus = min(nj-1,(int)floor(x1_s[1]+pixel_eps));
    for (int j = j0_plus; j <= j1_minus; ++j) {
      // get the x01 as a linear weight on the x0_s->x1_s line...
      const float w01 = max(0.f,min(1.f,(x1_s[1]-float(j))/(x1_s[1] - x0_s[1])));
      const float x01 = w01*x0_s[0] + (1.0f-w01)*x1_s[0];
      const float w02 = max(0.f,min(1.f,(x2_s[1]-float(j))/(x2_s[1] - x0_s[1])));
      const float x02 = w02*x0_s[0] + (1.0f-w02)*x2_s[0];
      // the x values can be in either order...
      if (x01 < x02) {
        const int i_f = max(0,(int)ceil(x01-pixel_eps));
        const int i_l = min(ni-1,(int)floor(x02+pixel_eps));
        if (i_f <= i_l) {
          const float d01 = w01*x0_s[2] + (1.0f-w01)*x1_s[2];
          const float d02 = w02*x0_s[2] + (1.0f-w02)*x2_s[2];
          for (int i = i_f; i <= i_l; ++i) {
            const float w = max(0.f,min(1.f,(x02 - float(i))/(x02 - x01)));
            const float d = w*d01 + (1.0f-w)*d02;
            float sign = -1.0f;
            if (normal[2] >= 0.f) sign = 1.0f;
            const float a = sign*(normal[0]*light1[0] + normal[1]*light1[1] + normal[2]*light1[2]);
            const float b = sign*(normal[0]*light2[0] + normal[1]*light2[1] + normal[2]*light2[2]);
            const float p = ((a*base2[0] + (1.0-a)*base00[0])*0.5 + (b*base3[0] + (1.0-b)*base00[0])*0.5)/mag;
	    volVisDataVec.push_back(VolVisData(j*ni+i,d,show,dir,p));
          }
        }
      }
      else if (x02 < x01) {
        const int i_f = max(0,(int)ceil(x02-pixel_eps));
        const int i_l = min(ni-1,(int)floor(x01+pixel_eps));
        if (i_f <= i_l) {
          const float d01 = w01*x0_s[2] + (1.0f-w01)*x1_s[2];
          const float d02 = w02*x0_s[2] + (1.0f-w02)*x2_s[2];
          for (int i = i_f; i <= i_l; ++i) {
            const float w = max(0.f,min(1.f,(x01 - float(i))/(x01 - x02)));
            const float d = w*d02 + (1.0f-w)*d01;
            float sign = -1.0f;
            if (normal[2] >= 0.f) sign = 1.0f;
            const float a = sign*(normal[0]*light1[0] + normal[1]*light1[1] + normal[2]*light1[2]);
            const float b = sign*(normal[0]*light2[0] + normal[1]*light2[1] + normal[2]*light2[2]);
            const float p = ((a*base2[0] + (1.0-a)*base00[0])*0.5 + (b*base3[0] + (1.0-b)*base00[0])*0.5)/mag;
	    volVisDataVec.push_back(VolVisData(j*ni+i,d,show,dir,p));
          }
        }
      }
    }
  }

  if (x1_s[1] < x2_s[1]) {
    const int j1_plus = max(0,(int)ceil(x1_s[1]-pixel_eps));
    const int j2_minus = min(nj-1,(int)floor(x2_s[1]+pixel_eps));
    for (int j = j1_plus; j <= j2_minus; ++j) {
      // get the x12 as a linear weight on the x1_s->x2_s line...
      const float w12 = max(0.f,min(1.f,(x2_s[1]-float(j))/(x2_s[1] - x1_s[1])));
      const float x12 = w12*x1_s[0] + (1.0f-w12)*x2_s[0];
      const float w02 = max(0.f,min(1.f,(x2_s[1]-float(j))/(x2_s[1] - x0_s[1])));
      const float x02 = w02*x0_s[0] + (1.0f-w02)*x2_s[0];
      // the x values can be in either order...
      if (x12 < x02) {
        const int i_f = max(0,(int)ceil(x12-pixel_eps));
        const int i_l = min(ni-1,(int)floor(x02+pixel_eps));
        if (i_f <= i_l) {
          const float d12 = w12*x1_s[2] + (1.0f-w12)*x2_s[2];
          const float d02 = w02*x0_s[2] + (1.0f-w02)*x2_s[2];
          for (int i = i_f; i <= i_l; ++i) {
            const float w = max(0.f,min(1.f,(x02 - float(i))/(x02 - x12)));
            const float d = w*d12 + (1.0f-w)*d02;
            float sign = -1.0f;
            if (normal[2] >= 0.f) sign = 1.0f;
            const float a = sign*(normal[0]*light1[0] + normal[1]*light1[1] + normal[2]*light1[2]);
            const float b = sign*(normal[0]*light2[0] + normal[1]*light2[1] + normal[2]*light2[2]);
            const float p = ((a*base2[0] + (1.0-a)*base00[0])*0.5 + (b*base3[0] + (1.0-b)*base00[0])*0.5)/mag;
	    volVisDataVec.push_back(VolVisData(j*ni+i,d,show,dir,p));
          }
        }
      }
      else if (x02 < x12) {
        const int i_f = max(0,(int)ceil(x02-pixel_eps));
        const int i_l = min(ni-1,(int)floor(x12+pixel_eps));
        if (i_f <= i_l) {
          const float d12 = w12*x1_s[2] + (1.0f-w12)*x2_s[2];
          const float d02 = w02*x0_s[2] + (1.0f-w02)*x2_s[2];
          for (int i = i_f; i <= i_l; ++i) {
            const float w = max(0.f,min(1.f,(x12 - float(i))/(x12 - x02)));
            const float d = w*d02 + (1.0f-w)*d12;
            float sign = -1.0f;
            if (normal[2] >= 0.f) sign = 1.0f;
            const float a = sign*(normal[0]*light1[0] + normal[1]*light1[1] + normal[2]*light1[2]);
            const float b = sign*(normal[0]*light2[0] + normal[1]*light2[1] + normal[2]*light2[2]);
            const float p = ((a*base2[0] + (1.0-a)*base00[0])*0.5 + (b*base3[0] + (1.0-b)*base00[0])*0.5)/mag;
	    volVisDataVec.push_back(VolVisData(j*ni+i,d,show,dir,p));
	  }
        }
      }
    }
  }
}

void CtiCanvas::addVoronoiPointVolVis(const double x_vv[3],const double r_vv,const float v) {

  // for every i,j ray that could possibly pass through the sphere defined by x_vv, with radius r_vv

  const float rp = r_vv*float(ni)/width;
  assert(rp > 0.f);

  const float ic = CALC_PIXEL_I(x_vv);
  if ( (ic+rp < 0.f) || (ic-rp > float(ni-1)) ) return;

  const float jc = CALC_PIXEL_J(x_vv);
  if ( (jc+rp < 0.f) || (jc-rp > float(nj-1)) ) return;

  // we can now do the face geometry...
  const float dc = CALC_PIXEL_DEPTH(x_vv); // depth of the x_vv

  // now loop on rays....
  const int i0 = max(0,(int)ceil(ic-rp));
  const int i1 = min(ni-1,(int)floor(ic+rp));
  const int j0 = max(0,(int)ceil(jc-rp));
  const int j1 = min(nj-1,(int)floor(jc+rp));
  for (int i = i0; i <= i1; ++i) {
    for (int j = j0; j <= j1; ++j) {
      const double rij2 = (float(i)-ic)*(float(i)-ic)+(float(j)-jc)*(float(j)-jc);
      if (rij2 < rp*rp) {
        // store the location of our intersection, d, AND our square distance from the ray...
        volVisDataVec.push_back(VolVisData(j*ni+i,dc,v,rij2,0)); // here we use the v1 spot in the structure to store our distances (squared), and use flag=-1
      }
    }
  }

}

void CtiCanvas::addEdge(const double xe0[3],const double xe1[3],const double dx[3],const uint2 zone,const bool barbell) {

  const float r = 3.f; // cylinder pixel radius
  const float l_width_ni_ratio = 2.f*width/(float)ni;

  const float dd = l_width_ni_ratio*(DOT_PRODUCT(dx,e2));
  const float d0 = CALC_PIXEL_DEPTH(xe0)+dd;
  const float d1 = CALC_PIXEL_DEPTH(xe1)+dd;

  if (max(d0,d1) <= 0.f) return;

  const float di = l_width_ni_ratio*(DOT_PRODUCT(dx,e0));
  const float i0 = CALC_PIXEL_I(xe0)+di;
  const float i1 = CALC_PIXEL_I(xe1)+di;

  if ( (max(i0,i1) < 0.f) || (min(i0,i1) > float(ni-1)) ) return;

  const float dj = l_width_ni_ratio*(DOT_PRODUCT(dx,e1));
  const float j0 = CALC_PIXEL_J(xe0)+dj;
  const float j1 = CALC_PIXEL_J(xe1)+dj;

  if ( (max(j0,j1) < 0.f) || (min(j0,j1) > float(nj-1)) ) return;

  // render...

  float normal[3];

  const float denom = (i1-i0)*(i1-i0) + (j1-j0)*(j1-j0);
  const float r_inv = 1.f/r;

  if (barbell) {

    const float rb = 2.f*r;
    const float rb_over_sqrt_denom = rb/sqrt(denom);
    const float rb_inv = 1.f/rb;

    if (denom == 0.f) {

      // render point at i0,j0...
      const int i_f = max(0,(int)ceil(i0-rb));
      const int i_l = min(ni-1,(int)floor(i0+rb));
      const int j_f = max(0,(int)ceil(j0-rb));
      const int j_l = min(nj-1,(int)floor(j0+rb));
      for (int i = i_f; i <= i_l; ++i) {
        for (int j = j_f; j <= j_l; ++j) {
          normal[0] = float(i)-i0;
          normal[1] = float(j)-j0;
          normal[2] = rb*rb - normal[0]*normal[0] - normal[1]*normal[1];
          if (normal[2] > 0.f) {
            normal[2] = sqrt(normal[2])*rb_inv;
            const float d = d0-normal[2]*rb;
            if ((d > 0.f)&&(!isBlanked(float(i),float(j),d))) {
              // finish normalization...
              normal[0] *= rb_inv;
              normal[1] *= rb_inv;
              int IJ,ij;
              getIJij(IJ,ij,i,j);
              float sdepth;
              uint2 bits;
              if (getPixelSdepthBits(sdepth,bits,IJ,ij)) {
                if (bits&SURFACE_BITS) {
                  if (d < sdepth) {
                    // closer than original surface
                    bits &= ~SURFACE_DATA_BIT; // clear data bit
                    bits |= SURFACE_BIT;  // set surface bit
                    setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,bits);
                  }
                }
                else {
                  // originally background
                  bits |= SURFACE_BIT;
                  setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,bits);
                }
              }
              else {
                pbd[IJ] = new PixelBlockData();
                // treat like originally background
                setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,SURFACE_BIT);
              }
            }
          }
        }
      }

    }
    else {

      // render edge...
      const int i_f = max(0,(int)ceil(min(i0,i1)-rb));
      const int i_l = min(ni-1,(int)floor(max(i0,i1)+rb));
      const int j_f = max(0,(int)ceil(min(j0,j1)-rb));
      const int j_l = min(nj-1,(int)floor(max(j0,j1)+rb));
      for (int i = i_f; i <= i_l; ++i) {
        for (int j = j_f; j <= j_l; ++j) {
          const float t = ((i1-i0)*(float(i)-i0)+(j1-j0)*(float(j)-j0))/denom;
          if (t <= rb_over_sqrt_denom) {
            // sphere at i0,j0...
            normal[0] = float(i)-i0;
            normal[1] = float(j)-j0;
            normal[2] = rb*rb - normal[0]*normal[0] - normal[1]*normal[1];
            if (normal[2] > 0.f) {
              normal[2] = sqrt(normal[2])*rb_inv;
              const float d = d0-normal[2]*rb;
              if ((d > 0.f)&&(!isBlanked(float(i),float(j),d))) {
                normal[0] *= rb_inv;
                normal[1] *= rb_inv;
                int IJ,ij;
                getIJij(IJ,ij,i,j);
                float sdepth;
                uint2 bits;
                if (getPixelSdepthBits(sdepth,bits,IJ,ij)) {
                  if (bits&SURFACE_BITS) {
                    if (d < sdepth) {
                      // closer than original surface
                      bits &= ~SURFACE_DATA_BIT; // clear data bit
                      bits |= SURFACE_BIT;  // set surface bit
                      setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,bits);
                    }
                  }
                  else {
                    // originally background
                    bits |= SURFACE_BIT;
                    setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,bits);
                  }
                }
                else {
                  pbd[IJ] = new PixelBlockData();
                  // treat like originally background
                  setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,SURFACE_BIT);
                }
              }
            }
          }
          else if (t >= (1.f-rb_over_sqrt_denom)) {
            // sphere at i1,j1...
            normal[0] = float(i)-i1;
            normal[1] = float(j)-j1;
            normal[2] = rb*rb - normal[0]*normal[0] - normal[1]*normal[1];
            if (normal[2] > 0.f) {
              normal[2] = sqrt(normal[2])*rb_inv;
              const float d = d1-normal[2]*rb;
              if ((d > 0.f)&&(!isBlanked(float(i),float(j),d))) {
                normal[0] *= rb_inv;
                normal[1] *= rb_inv;
                int IJ,ij;
                getIJij(IJ,ij,i,j);
                float sdepth;
                uint2 bits;
                if (getPixelSdepthBits(sdepth,bits,IJ,ij)) {
                  if (bits&SURFACE_BITS) {
                    if (d < sdepth) {
                      // closer than original surface
                      bits &= ~SURFACE_DATA_BIT; // clear data bit
                      bits |= SURFACE_BIT;  // set surface bit
                      setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,bits);
                    }
                  }
                  else {
                    // originally background
                    bits |= SURFACE_BIT;
                    setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,bits);
                  }
                }
                else {
                  pbd[IJ] = new PixelBlockData();
                  // treat like originally background
                  setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,SURFACE_BIT);
                }
              }
            }
          }
          // overlaps a bit (relation between spherical and cylindrical coordinates for r/rb = 0.5)
          if ((t >= 0.5f*rb_over_sqrt_denom)&&(t <= (1.f-0.5f*rb_over_sqrt_denom))) {
            // along cylinder...
            normal[0] = float(i)-((1.f-t)*i0 + t*i1);
            normal[1] = float(j)-((1.f-t)*j0 + t*j1);
            normal[2] = r*r - normal[0]*normal[0] - normal[1]*normal[1];
            if (normal[2] > 0.f) {
              normal[2] = sqrt(normal[2])*r_inv;
              const float d = (1.f-t)*d0+t*d1-normal[2]*r;
              if ((d > 0.f)&&(!isBlanked(float(i),float(j),d))) {
                normal[0] *= r_inv;
                normal[1] *= r_inv;
                int IJ,ij;
                getIJij(IJ,ij,i,j);
                float sdepth;
                uint2 bits;
                if (getPixelSdepthBits(sdepth,bits,IJ,ij)) {
                  if (bits&SURFACE_BITS) {
                    if (d < sdepth) {
                      // closer than original surface
                      bits &= ~SURFACE_DATA_BIT; // clear data bit
                      bits |= SURFACE_BIT;  // set surface bit
                      setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,bits);
                    }
                  }
                  else {
                    // originally background
                    bits |= SURFACE_BIT;
                    setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,bits);
                  }
                }
                else {
                  pbd[IJ] = new PixelBlockData();
                  // treat like originally background
                  setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,SURFACE_BIT);
                }
              }
            }
          }
        }
      }
    }
  }
  else {

    if (denom == 0.f) {

      // render point at i0,j0...
      const int i_f = max(0,(int)ceil(i0-r));
      const int i_l = min(ni-1,(int)floor(i0+r));
      const int j_f = max(0,(int)ceil(j0-r));
      const int j_l = min(nj-1,(int)floor(j0+r));
      for (int i = i_f; i <= i_l; ++i) {
        for (int j = j_f; j <= j_l; ++j) {
          normal[0] = float(i)-i0;
          normal[1] = float(j)-j0;
          normal[2] = r*r - normal[0]*normal[0] - normal[1]*normal[1];
          if (normal[2] > 0.f) {
            normal[2] = sqrt(normal[2])*r_inv;
            const float d = d0-normal[2]*r;
            if ((d > 0.f)&&(!isBlanked(float(i),float(j),d))) {
              // finish normalization...
              normal[0] *= r_inv;
              normal[1] *= r_inv;
              int IJ,ij;
              getIJij(IJ,ij,i,j);
              float sdepth;
              uint2 bits;
              if (getPixelSdepthBits(sdepth,bits,IJ,ij)) {
                if (bits&SURFACE_BITS) {
                  if (d < sdepth) {
                    // closer than original surface
                    bits &= ~SURFACE_DATA_BIT; // clear data bit
                    bits |= SURFACE_BIT;  // set surface bit
                    setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,bits);
                  }
                }
                else {
                  // originally background
                  bits |= SURFACE_BIT;
                  setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,bits);
                }
              }
              else {
                pbd[IJ] = new PixelBlockData();
                // treat like originally background
                setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,SURFACE_BIT);
              }
            }
          }
        }
      }

    }
    else {

      // render edge...
      // NOT quite right yet...
      const int i_f = max(0,(int)ceil(min(i0,i1)-r));
      const int i_l = min(ni-1,(int)floor(max(i0,i1)+r));
      const int j_f = max(0,(int)ceil(min(j0,j1)-r));
      const int j_l = min(nj-1,(int)floor(max(j0,j1)+r));
      for (int i = i_f; i <= i_l; ++i) {
        for (int j = j_f; j <= j_l; ++j) {
          const float t = ((i1-i0)*(float(i)-i0)+(j1-j0)*(float(j)-j0))/denom;
          if (t <= 0.f) {
            // sphere at i0,j0...
            normal[0] = float(i)-i0;
            normal[1] = float(j)-j0;
            normal[2] = r*r - normal[0]*normal[0] - normal[1]*normal[1];
            if (normal[2] > 0.f) {
              normal[2] = sqrt(normal[2])*r_inv;
              const float d = d0-normal[2]*r;
              if ((d > 0.f)&&(!isBlanked(float(i),float(j),d))) {
                normal[0] *= r_inv;
                normal[1] *= r_inv;
                int IJ,ij;
                getIJij(IJ,ij,i,j);
                float sdepth;
                uint2 bits;
                if (getPixelSdepthBits(sdepth,bits,IJ,ij)) {
                  if (bits&SURFACE_BITS) {
                    if (d < sdepth) {
                      // closer than original surface
                      bits &= ~SURFACE_DATA_BIT; // clear data bit
                      bits |= SURFACE_BIT;  // set surface bit
                      setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,bits);
                    }
                  }
                  else {
                    // originally background
                    bits |= SURFACE_BIT;
                    setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,bits);
                  }
                }
                else {
                  pbd[IJ] = new PixelBlockData();
                  // treat like originally background
                  setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,SURFACE_BIT);
                }
              }
            }
          }
          else if (t >= 1.f) {
            // sphere at i1,j1...
            normal[0] = float(i)-i1;
            normal[1] = float(j)-j1;
            normal[2] = r*r - normal[0]*normal[0] - normal[1]*normal[1];
            if (normal[2] > 0.f) {
              normal[2] = sqrt(normal[2])*r_inv;
              const float d = d1-normal[2]*r;
              if ((d > 0.f)&&(!isBlanked(float(i),float(j),d))) {
                normal[0] *= r_inv;
                normal[1] *= r_inv;
                int IJ,ij;
                getIJij(IJ,ij,i,j);
                float sdepth;
                uint2 bits;
                if (getPixelSdepthBits(sdepth,bits,IJ,ij)) {
                  if (bits&SURFACE_BITS) {
                    if (d < sdepth) {
                      // closer than original surface
                      bits &= ~SURFACE_DATA_BIT; // clear data bit
                      bits |= SURFACE_BIT;  // set surface bit
                      setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,bits);
                    }
                  }
                  else {
                    // originally background
                    bits |= SURFACE_BIT;
                    setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,bits);
                  }
                }
                else {
                  pbd[IJ] = new PixelBlockData();
                  // treat like originally background
                  setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,SURFACE_BIT);
                }
              }
            }
          }
          else {
            // along cylinder...
            normal[0] = float(i)-((1.f-t)*i0 + t*i1);
            normal[1] = float(j)-((1.f-t)*j0 + t*j1);
            normal[2] = r*r - normal[0]*normal[0] - normal[1]*normal[1];
            if (normal[2] > 0.f) {
              normal[2] = sqrt(normal[2])*r_inv;
              const float d = (1.f-t)*d0+t*d1-normal[2]*r;
              if ((d > 0.f)&&(!isBlanked(float(i),float(j),d))) {
                normal[0] *= r_inv;
                normal[1] *= r_inv;
                int IJ,ij;
                getIJij(IJ,ij,i,j);
                float sdepth;
                uint2 bits;
                if (getPixelSdepthBits(sdepth,bits,IJ,ij)) {
                  if (bits&SURFACE_BITS) {
                    if (d < sdepth) {
                      // closer than original surface
                      bits &= ~SURFACE_DATA_BIT; // clear data bit
                      bits |= SURFACE_BIT;  // set surface bit
                      setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,bits);
                    }
                  }
                  else {
                    // originally background
                    bits |= SURFACE_BIT;
                    setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,bits);
                  }
                }
                else {
                  pbd[IJ] = new PixelBlockData();
                  // treat like originally background
                  setPixelSnormalSdepthSzoneBits(IJ,ij,normal,d,zone,SURFACE_BIT);
                }
              }
            }
          }
        }
      }
    }
  }
}

void CtiCanvas::addVector(const vecType vType,const double xe0[3],const double dx[3],const double factor,const uint2 zone) {
  assert(factor > 0.0);
  const float mag = MAG(dx);
  double xe1[3]; FOR_I3 xe1[i] = xe0[i] + factor*dx[i]; // head

  const float r_eps = 1.0E-6f; // smallest r
  const float tfrac = 0.66f; // fraction of arrow that is cylindrical
  const float r = 0.02f*factor*mag*float(ni)/width; // cylinder radius
  const float rb = 2.f*r; // max radius of arrow head

  const float d0 = CALC_PIXEL_DEPTH(xe0);
  const float d1 = CALC_PIXEL_DEPTH(xe1);

  if (max(d0,d1) <= 0.f) return;

  const float i0 = CALC_PIXEL_I(xe0);
  const float i1 = CALC_PIXEL_I(xe1);

  if ( (max(i0,i1) < 0.f) || (min(i0,i1) > float(ni-1)) ) return;

  const float j0 = CALC_PIXEL_J(xe0);
  const float j1 = CALC_PIXEL_J(xe1);

  if ( (max(j0,j1) < 0.f) || (min(j0,j1) > float(nj-1)) ) return;

  bool (CtiCanvas::*getPixelDepth_ptr)(float& ,uint2& ,const int ,const int) const;
  void (CtiCanvas::*setPixel_ptr)(const int ,const int ,const float * ,const float ,const float ,const uint2 ,const uint2 );
  uint2 bit_type_check,bit_data_type,bit_surface_clear;
  if (vType == SURFACE_VECTOR) {
    bit_type_check = SURFACE_BITS;
    bit_data_type = SURFACE_DATA_BIT;
    bit_surface_clear = SURFACE_BIT;
    getPixelDepth_ptr = &CtiCanvas::getPixelSdepthBits;
    setPixel_ptr = &CtiCanvas::setPixelSnormalSdepthSdataSzoneBits;
  }
  else {
    assert(vType == ISOSURFACE_VECTOR);  // for now only two types...
    bit_type_check = INTERNAL_BITS;
    bit_data_type = INTERNAL_DATA_BIT;
    bit_surface_clear = ISURFACE_BIT;
    getPixelDepth_ptr = &CtiCanvas::getPixelPdepthBits;
    setPixel_ptr = &CtiCanvas::setPixelPnormalPdepthPdataPzoneBits;
  }
  // render...

  const float denom = (i1-i0)*(i1-i0) + (j1-j0)*(j1-j0);
  if (denom != 0.f) {
    float normal[3];
    const float r_inv = 1.f/r;
    const float rb_inv = 1.f/rb;

    // render edge...
    const int i_f = max(0,(int)ceil(min(i0,i1)-rb));
    const int i_l = min(ni-1,(int)floor(max(i0,i1)+rb));
    const int j_f = max(0,(int)ceil(min(j0,j1)-rb));
    const int j_l = min(nj-1,(int)floor(max(j0,j1)+rb));
    for (int i = i_f; i <= i_l; ++i) {
      for (int j = j_f; j <= j_l; ++j) {
        const float t = ((i1-i0)*(float(i)-i0)+(j1-j0)*(float(j)-j0))/denom;
        if (t <= 0.f) {
          // sphere at i0,j0...
          normal[0] = float(i)-i0;
          normal[1] = float(j)-j0;
          normal[2] = r*r - normal[0]*normal[0] - normal[1]*normal[1];
          if (normal[2] > 0.f) {
            normal[2] = sqrt(normal[2])*r_inv;
            const float d = d0-normal[2]*r;
            if ((d > 0.f)&&(!isBlanked(float(i),float(j),d))) {
              normal[0] *= r_inv;
              normal[1] *= r_inv;
              int IJ,ij;
              getIJij(IJ,ij,i,j);
              float depth;
              uint2 bits;
              if ((this->*getPixelDepth_ptr)(depth,bits,IJ,ij)) {
                if (bits&bit_type_check) {
                  if (d < depth) {
                    // closer than current internal data
                    bits |= bit_data_type;  // set data bit
                    bits &= ~bit_surface_clear; // clear surface bit
                    (this->*setPixel_ptr)(IJ,ij,normal,d,mag,zone,bits);
                  }
                }
                else {
                  // first internal data
                  (this->*setPixel_ptr)(IJ,ij,normal,d,mag,zone,bits|bit_data_type);
                }
              }
              else {
                pbd[IJ] = new PixelBlockData();
                // treat like first internal data
                (this->*setPixel_ptr)(IJ,ij,normal,d,mag,zone,bit_data_type);
              }
            }
          }
        }
        else if ( (t > 0.f)&&(t <= tfrac) ) {
          // along cylinder...
          normal[0] = float(i)-((1.f-t)*i0 + t*i1);
          normal[1] = float(j)-((1.f-t)*j0 + t*j1);
          normal[2] = r*r - normal[0]*normal[0] - normal[1]*normal[1];
          if (normal[2] > 0.f) {
            normal[2] = sqrt(normal[2])*r_inv;
            const float d = (1.f-t)*d0+t*d1-normal[2]*r;
            if ((d > 0.f)&&(!isBlanked(float(i),float(j),d))) {
              normal[0] *= r_inv;
              normal[1] *= r_inv;
              int IJ,ij;
              getIJij(IJ,ij,i,j);
              float depth;
              uint2 bits;
              if ((this->*getPixelDepth_ptr)(depth,bits,IJ,ij)) {
                if (bits&bit_type_check) {
                  if (d < depth) {
                    // closer than current internal data
                    bits |= bit_data_type;  // set data bit
                    bits &= ~bit_surface_clear; // clear surface bit
                    (this->*setPixel_ptr)(IJ,ij,normal,d,mag,zone,bits);
                  }
                }
                else {
                  // first internal data
                  (this->*setPixel_ptr)(IJ,ij,normal,d,mag,zone,bits|bit_data_type);
                }
              }
              else {
                pbd[IJ] = new PixelBlockData();
                // treat like first internal data
                (this->*setPixel_ptr)(IJ,ij,normal,d,mag,zone,bit_data_type);
              }
            }
          }
        }
        else if ( (t > tfrac)&&(t <= 1.f) ) {
          const float this_r = max(r_eps,rb*(1.f-(t-tfrac)/(1.f-tfrac)));
          const float this_r_inv = 1.f/rb_inv;

          // along cone..
          normal[0] = float(i)-((1.f-t)*i0 + t*i1);
          normal[1] = float(j)-((1.f-t)*j0 + t*j1);
          normal[2] = this_r*this_r - normal[0]*normal[0] - normal[1]*normal[1];
          if (normal[2] > 0.f) {
            normal[2] = sqrt(normal[2])*this_r_inv;
            const float d = (1.f-t)*d0+t*d1-normal[2]*this_r;
            if ((d > 0.f)&&(!isBlanked(float(i),float(j),d))) {
              normal[0] *= this_r_inv;
              normal[1] *= this_r_inv;
              int IJ,ij;
              getIJij(IJ,ij,i,j);
              float depth;
              uint2 bits;
              if ((this->*getPixelDepth_ptr)(depth,bits,IJ,ij)) {
                if (bits&bit_type_check) {
                  if (d < depth) {
                    // closer than current internal data
                    bits |= bit_data_type;  // set data bit
                    bits &= ~bit_surface_clear; // clear surface bit
                    (this->*setPixel_ptr)(IJ,ij,normal,d,mag,zone,bits);
                  }
                }
                else {
                  // first internal data
                  (this->*setPixel_ptr)(IJ,ij,normal,d,mag,zone,bits|bit_data_type);
                }
              }
              else {
                pbd[IJ] = new PixelBlockData();
                // treat like first internal data
                (this->*setPixel_ptr)(IJ,ij,normal,d,mag,zone,bit_data_type);
              }
            }
          }
        }
      }
    }
  }
}

void CtiCanvas::calcVolvisIntensityAndDepth(float &intensity,float &depth,vector<VolVisData>::iterator& iter_begin,vector<VolVisData>::iterator& iter_end) const {

  // recall flag can be 0,1,2 for internal face, or
  // flag = 3 +
  // if (show) flag |= 4;
  // if (dir) flag |= 8;

  const double intensity_eps = 1.0E-32;

  // integrate in double, then return float
  double intensity_dble = 0.0;
  double fabs_intensity_dble = 0.0;
  double depth_dble = 0.0;
  double intensity_current = 0.0;
  float depth_current = iter_begin->depth;
  float delta_sum = 0.f;

  if (volvis_type == "X-RAY") {

    //if (debug) cout << "starting depth: " << depth_current << endl;
    bool b_closed = true;
    for (vector<VolVisData>::iterator iter = iter_begin; iter != iter_end; ++iter) {

      // surface
      if ((iter->flag & 3) == 3) {
        // the 8 bit determines the direction...
        if (iter->flag & 8) {
          // this is a start.
          if (!b_closed)
            delta_sum += depth_current-iter->depth; // opposite sign due to direction of depth
          depth_current = iter->depth;
          b_closed = false; // we use this just in case pixel round off switches order b/w open and close
        }
        else {
          // this is an end.
          delta_sum += depth_current-iter->depth; // opposite sign due to direction of depth
          intensity_dble += delta_sum;
          depth_dble += 0.5*(depth_current+iter->depth)*(intensity_dble-intensity_current);
          intensity_current = intensity_dble;
          depth_current = iter->depth;
          delta_sum = 0.f;
          b_closed = true;
        }
        // the 4 bit determines show or no-show...
        if (iter->flag & 4) {
          // if the surface is "show", then we clear out the integrated intensity and depth
          intensity_dble = 0.0;
          depth_dble = 0.0;
        }
        else {
          //if (debug) cout << "This is no-show" << endl;
        }
      }
      // internal/interproc face (assumed empty!)
      else {
        assert(0);
      }

    }

  }
  else if (volvis_type == "LINEAR") {

    // start by counting the number of "intersections" and also ensure there are some
    // in/out points. It is very possible to have rays that only contain internal proximity
    // points, and these cannot contribute to the intesity...
    int np_internal = 0;
    int np_surface_start = 0;
    int np_surface_end = 0;
    for (vector<VolVisData>::iterator iter = iter_begin; iter != iter_end; ++iter) {
      //cout << "iter->depth: " << iter->depth << " iter->flag: " << iter->flag << endl;
      if (iter->flag == 0) {
        ++np_internal;
      }
      else if ((iter->flag & 3) == 3) {
        // the 8 bit determines the direction...
        if (iter->flag & 8) {
          ++np_surface_start;
        }
        else {
          ++np_surface_end;
        }
      }
      else {
        assert(0);
      }
    }

    if ((np_surface_start > 0)&&(np_surface_end > 0)&&(np_internal > 0)) {

      // we can compute the intensity along this ray...

      // start by transfering everything to arrays...
      const int np = np_surface_start + np_surface_end + np_internal;
      int * flag = new int[np];
      float * depth = new float[np];
      float * v = new float[np];
      float * d2 = new float[np];
      int ip = 0;
      for (vector<VolVisData>::iterator iter = iter_begin; iter != iter_end; ++iter) {
        // all have flag and depth...
        flag[ip] = iter->flag;
        depth[ip] = iter->depth;
        if (iter->flag == 0) {
          // for internal, also store the data value and the 2d off the ray...
          v[ip] = iter->v0;
          d2[ip] = iter->v1;
        }
        ++ip;
      }
      assert(ip == np);

      // now set the index of the closest internal point at each point. We can do this with a 2-sweep
      // algo: once down, once back...

      int * index = new int[np];
      for (int ip = 0; ip < np; ++ip) {
        if (flag[ip] == 0)
          index[ip] = ip;
        else
          index[ip] = -1; // start the surface points without any index...
      }

      // forward sweep...
      for (int ip = 0; ip < np-1; ++ip) {
        const int ipp1 = ip+1;
        if ((index[ip] == -1)&&(index[ipp1] != -1)) {
          index[ip] = index[ipp1];
        }
        else if ((index[ip] != -1)&&(index[ipp1] == -1)) {
          index[ipp1] = index[ip];
        }
        else if ((index[ip] != -1)&&(index[ipp1] != -1)) {
          assert(index[ipp1] == ipp1); // has not been touched
          // both points refer to a valid internal point. compute the
          // various distances...
          const float d2_ip_ipp1 = d2[index[ip]] + (depth[index[ip]]-depth[ipp1])*(depth[index[ip]]-depth[ipp1]);
          const float d2_ipp1_ipp1 = d2[ipp1];
          if (d2_ip_ipp1 < d2_ipp1_ipp1) {
            index[ipp1] = index[ip];
          }
        }
      }

      // backward sweep...
      for (int ip = np-1; ip > 0; --ip) {
        const int ipm1 = ip-1;
        if ((index[ip] == -1)&&(index[ipm1] != -1)) {
          assert(0);
        }
        else if ((index[ip] != -1)&&(index[ipm1] == -1)) {
          index[ipm1] = index[ip];
        }
        else if ((index[ip] != -1)&&(index[ipm1] != -1)) {
          // both points refer to a valid internal point. compute the
          // various distances...
          const float d2_ip_ipm1 = d2[index[ip]] + (depth[index[ip]]-depth[ipm1])*(depth[index[ip]]-depth[ipm1]);
          const float d2_ipm1_ipm1 = d2[index[ipm1]] + (depth[index[ipm1]]-depth[ipm1])*(depth[index[ipm1]]-depth[ipm1]);
          if (d2_ip_ipm1 < d2_ipm1_ipm1) {
            index[ipm1] = index[ip];
          }
        }
      }

      // now perform the integration. Instead of using inside/outside, we count the
      // start and end's positively and negatively respectively. In this way, when the
      // count is greater than 1, we have more starts than ends, and we can add to the
      // integration...
      bool inside = false;
      for (int ip = 0; ip < np; ++ip) {
	// this is a surface intersection. See if it is start or end...
	if ((flag[ip] & 3) == 3) {
	  // the 8 bit determines the direction...
	  if (flag[ip] & 8) {
	    // start...
	    inside = true;
	  }
	  else {
	    // end...
	    inside = false;
	  }
	  if (flag[ip] & 4) {
	    //if (debug) cout << "This is show" << endl;
	    // if the surface is "show", then we clear out the integrated intensity and depth
	    intensity_dble = 0.0;
	    depth_dble = 0.0;
	    fabs_intensity_dble = 0.0;
	  }
	}
	if (ip < np-1) {
	  const int ipp1 = ip+1;
	  if (inside) {
	    // we are inside the geometry, so include this interval in the integration...
	    if (index[ip] == index[ipp1]) {
	      // same value at both ends of the interval...
	      const double this_intensity = max(min(v[index[ip]],volvis_aux_data[1]),volvis_aux_data[0])*(depth[ip]-depth[ipp1]);
	      intensity_dble += this_intensity;
	      const double this_fabs_intensity = fabs(this_intensity)+intensity_eps*(depth[ip]-depth[ipp1]);
	      fabs_intensity_dble += this_fabs_intensity;
	      depth_dble += this_fabs_intensity*0.5*(depth[ip]+depth[ipp1]);
	    }
	    else {
	      // value changes across this interval: somewhere in this interval there is a depth_mid
	      // that is equidistant to both the index points: see Maple worksheet depth_mid.mw...
	      const float denom = depth[index[ip]]-depth[index[ipp1]];
	      if (denom > 0.f) {
		const float depth_mid = 0.5*(depth[index[ip]]+depth[index[ipp1]]+(d2[index[ip]]-d2[index[ipp1]])/denom);
		// first part of interval associated with ip...
		{
		  const double this_intensity = max(min(v[index[ip]],volvis_aux_data[1]),volvis_aux_data[0])*(depth[ip]-depth_mid);
		  intensity_dble += this_intensity;
		  const double this_fabs_intensity = fabs(this_intensity)+intensity_eps*fabs(depth[ip]-depth_mid); // can get flipped due to machine precision
		  fabs_intensity_dble += this_fabs_intensity;
		  depth_dble += this_fabs_intensity*0.5*(depth[ip]+depth_mid);
		}
		// second part of interval associated with ipp1...
		{
		  const double this_intensity = max(min(v[index[ipp1]],volvis_aux_data[1]),volvis_aux_data[0])*(depth_mid-depth[ipp1]);
		  intensity_dble += this_intensity;
		  const double this_fabs_intensity = fabs(this_intensity)+intensity_eps*fabs(depth_mid-depth[ipp1]); // can get flipped due to machine precision
		  fabs_intensity_dble += this_fabs_intensity;
		  depth_dble += this_fabs_intensity*0.5*(depth_mid+depth[ipp1]);
		}
	      }
	    }
	  }
	}
      }

      delete[] flag;
      delete[] depth;
      delete[] v;
      delete[] d2;
      delete[] index;

    }

  }
  else if (volvis_type == "GAUSSIAN") {

    // start by counting the number of "intersections" and also ensure there are some
    // in/out points. It is very possible to have rays that only contain internal proximity
    // points, and these cannot contribute to the intesity...
    int np_internal = 0;
    int np_surface_start = 0;
    int np_surface_end = 0;
    for (vector<VolVisData>::iterator iter = iter_begin; iter != iter_end; ++iter) {
      //cout << "iter->depth: " << iter->depth << " iter->flag: " << iter->flag << endl;
      if (iter->flag == 0) {
        ++np_internal;
      }
      else if ((iter->flag & 3) == 3) {
        // the 8 bit determines the direction...
        if (iter->flag & 8) {
          ++np_surface_start;
        }
        else {
          ++np_surface_end;
        }
      }
      else {
        assert(0);
      }
    }

    if ((np_surface_start > 0)&&(np_surface_end > 0)&&(np_internal > 0)) {

      // we can compute the intensity along this ray...

      // start by transfering everything to arrays...
      const int np = np_surface_start + np_surface_end + np_internal;
      int * flag = new int[np];
      float * depth = new float[np];
      float * v = new float[np];
      float * d2 = new float[np];
      int ip = 0;
      for (vector<VolVisData>::iterator iter = iter_begin; iter != iter_end; ++iter) {
        // all have flag and depth...
        flag[ip] = iter->flag;
        depth[ip] = iter->depth;
        if (iter->flag == 0) {
          // for internal, also store the data value and the 2d off the ray...
          v[ip] = iter->v0;
          d2[ip] = iter->v1;
        }
        ++ip;
      }
      assert(ip == np);

      // now set the index of the closest internal point at each point. We can do this with a 2-sweep
      // algo: once down, once back...

      int * index = new int[np];
      for (int ip = 0; ip < np; ++ip) {
        if (flag[ip] == 0)
          index[ip] = ip;
        else
          index[ip] = -1; // start the surface points without any index...
      }

      // forward sweep...
      for (int ip = 0; ip < np-1; ++ip) {
        const int ipp1 = ip+1;
        if ((index[ip] == -1)&&(index[ipp1] != -1)) {
          index[ip] = index[ipp1];
        }
        else if ((index[ip] != -1)&&(index[ipp1] == -1)) {
          index[ipp1] = index[ip];
        }
        else if ((index[ip] != -1)&&(index[ipp1] != -1)) {
          assert(index[ipp1] == ipp1); // has not been touched
          // both points refer to a valid internal point. compute the
          // various distances...
          const float d2_ip_ipp1 = d2[index[ip]] + (depth[index[ip]]-depth[ipp1])*(depth[index[ip]]-depth[ipp1]);
          const float d2_ipp1_ipp1 = d2[ipp1];
          if (d2_ip_ipp1 < d2_ipp1_ipp1) {
            index[ipp1] = index[ip];
          }
        }
      }

      // backward sweep...
      for (int ip = np-1; ip > 0; --ip) {
        const int ipm1 = ip-1;
        if ((index[ip] == -1)&&(index[ipm1] != -1)) {
          assert(0);
        }
        else if ((index[ip] != -1)&&(index[ipm1] == -1)) {
          index[ipm1] = index[ip];
        }
        else if ((index[ip] != -1)&&(index[ipm1] != -1)) {
          // both points refer to a valid internal point. compute the
          // various distances...
          const float d2_ip_ipm1 = d2[index[ip]] + (depth[index[ip]]-depth[ipm1])*(depth[index[ip]]-depth[ipm1]);
          const float d2_ipm1_ipm1 = d2[index[ipm1]] + (depth[index[ipm1]]-depth[ipm1])*(depth[index[ipm1]]-depth[ipm1]);
          if (d2_ip_ipm1 < d2_ipm1_ipm1) {
            index[ipm1] = index[ip];
          }
        }
      }

      // now perform the integration. Instead of using inside/outside, we count the
      // start and end's positively and negatively respectively. In this way, when the
      // count is greater than 1, we have more starts than ends, and we can add to the
      // integration...
      bool inside = false;
      for (int ip = 0; ip < np; ++ip) {
	// this is a surface intersection. See if it is start or end...
	if ((flag[ip] & 3) == 3) {
	  // the 8 bit determines the direction...
	  if (flag[ip] & 8) {
	    // start...
	    inside = true;
	  }
	  else {
	    // end...
	    inside = false;
	  }
	  if (flag[ip] & 4) {
	    //if (debug) cout << "This is show" << endl;
	    // if the surface is "show", then we clear out the integrated intensity and depth
	    intensity_dble = 0.0;
	    depth_dble = 0.0;
	    fabs_intensity_dble = 0.0;
	  }
	}
	if (ip < np-1) {
	  const int ipp1 = ip+1;
	  if (inside) {
	    // we are inside the geometry, so include this interval in the integration...
	    if (index[ip] == index[ipp1]) {
	      // same value at both ends of the interval...
	      //const double this_intensity = max(min(v[index[ip]],volvis_aux_data[1]),volvis_aux_data[0])*(depth[ip]-depth[ipp1]);
              const double this_intensity = exp(-0.5*(v[index[ip]]-volvis_aux_data[0])*(v[index[ip]]-volvis_aux_data[0])/
                                                (volvis_aux_data[1]*volvis_aux_data[1]))/(volvis_aux_data[1]*sqrt(2.0*M_PI))*(depth[ip]-depth[ipp1]);
	      intensity_dble += this_intensity;
	      const double this_fabs_intensity = fabs(this_intensity)+intensity_eps*(depth[ip]-depth[ipp1]);
	      fabs_intensity_dble += this_fabs_intensity;
	      depth_dble += this_fabs_intensity*0.5*(depth[ip]+depth[ipp1]);
	    }
	    else {
	      // value changes across this interval: somewhere in this interval there is a depth_mid
	      // that is equidistant to both the index points: see Maple worksheet depth_mid.mw...
	      const float denom = depth[index[ip]]-depth[index[ipp1]];
	      if (denom > 0.f) {
		const float depth_mid = 0.5*(depth[index[ip]]+depth[index[ipp1]]+(d2[index[ip]]-d2[index[ipp1]])/denom);
		// first part of interval associated with ip...
		{
		  //const double this_intensity = max(min(v[index[ip]],volvis_aux_data[1]),volvis_aux_data[0])*(depth[ip]-depth_mid);
                  const double this_intensity = exp(-0.5*(v[index[ip]]-volvis_aux_data[0])*(v[index[ip]]-volvis_aux_data[0])/
                                                    (volvis_aux_data[1]*volvis_aux_data[1]))/(volvis_aux_data[1]*sqrt(2.0*M_PI))*(depth[ip]-depth_mid);
		  intensity_dble += this_intensity;
		  const double this_fabs_intensity = fabs(this_intensity)+intensity_eps*fabs(depth[ip]-depth_mid); // can get flipped due to machine precision
		  fabs_intensity_dble += this_fabs_intensity;
		  depth_dble += this_fabs_intensity*0.5*(depth[ip]+depth_mid);
		}
		// second part of interval associated with ipp1...
		{
		  //const double this_intensity = max(min(v[index[ipp1]],volvis_aux_data[1]),volvis_aux_data[0])*(depth_mid-depth[ipp1]);
                  const double this_intensity = exp(-0.5*(v[index[ipp1]]-volvis_aux_data[0])*(v[index[ipp1]]-volvis_aux_data[0])/
                                                    (volvis_aux_data[1]*volvis_aux_data[1]))/(volvis_aux_data[1]*sqrt(2.0*M_PI))*(depth_mid-depth[ipp1]);
		  intensity_dble += this_intensity;
		  const double this_fabs_intensity = fabs(this_intensity)+intensity_eps*fabs(depth_mid-depth[ipp1]); // can get flipped due to machine precision
		  fabs_intensity_dble += this_fabs_intensity;
		  depth_dble += this_fabs_intensity*0.5*(depth_mid+depth[ipp1]);
		}
	      }
	    }
	  }
	}
      }

      delete[] flag;
      delete[] depth;
      delete[] v;
      delete[] d2;
      delete[] index;

    }

  }
  else {
    assert(volvis_type == "SURFACE");

    //if (debug) cout << "starting depth: " << depth_current << endl;

    for (vector<VolVisData>::iterator iter = iter_begin; iter != iter_end; ++iter) {

      // surface
      if ((iter->flag & 3) == 3) {
        // the 8 bit determines the direction...
        if ((depth_current-iter->depth > 1.0E-6)||iter == iter_begin) {
          intensity_dble = intensity_dble*volvis_aux_data[0]+(1.0-volvis_aux_data[0])*iter->v0; // fresnel effect
          depth_dble += 0.5*(depth_current+iter->depth)*(intensity_dble-intensity_current);
          intensity_current = intensity_dble;
          depth_current = iter->depth;
          //depth_dble = depth_current;
        }
        if (iter->flag & 8) {
          // this is a start.
        }
        else {
          // this is an end.
        }
        // the 4 bit determines show or no-show...
        if (iter->flag & 4) {
          // if the surface is "show", then we clear out the integrated intensity and depth
          intensity_dble = 0.0;
          depth_dble = 0.0;
        }
        else {
          //if (debug) cout << "This is no-show" << endl;
        }
      }
      // internal/interproc face (assumed empty!)
      else {
        assert(0);
      }

    }
  }

  if (intensity_dble != 0.0) {
    intensity = float(intensity_dble);
    // TODO: hack: if fabs_intensity_dble not worked through above, set it to intensity_dble here...
    // at present, it is only worked through LINEAR...
    if (fabs_intensity_dble == 0.0)
      fabs_intensity_dble = intensity_dble;
  }
  else {
    intensity = 0.f;
  }

  if (fabs_intensity_dble != 0.0) {
    depth = float(depth_dble/fabs_intensity_dble); // center of intensity based depth
  }
  else {
    depth = 0.f;
  }

  return;
}

float CtiCanvas::calcPixelDepth(const float xp[3]) const {
  return CALC_PIXEL_DEPTH(xp);
}

float CtiCanvas::calcPixelDepth(const double xp[3]) const {
  return CALC_PIXEL_DEPTH(xp);
}

void CtiCanvas::setBlankDataPlane(const bool blank) {
  b_blank_data_plane = blank;
}

void CtiCanvas::setShowEdges(const bool show) {
  b_showEdges = show;
}

void CtiCanvas::setDataRange(const float &_min, const float &_max){
  range[0] = _min;
  range[1] = _max;
  b_range = true;
}
void CtiCanvas::setSurfaceDataRange(const float &_min, const float &_max){
  range_sVar[0] = _min;
  range_sVar[1] = _max;
  b_range_sVar = true;
}
void CtiCanvas::setParticleDataRange(const float &_min, const float &_max){
  range_pVar[0] = _min;
  range_pVar[1] = _max;
  b_range_pVar = true;
}
void CtiCanvas::setIsoDataRange(const float &_min, const float &_max){
  range_iVar[0] = _min;
  range_iVar[1] = _max;
  b_range_iVar = true;
}
void CtiCanvas::setDataRangeBins(const int nbin) {
  n_range = nbin;
}
void CtiCanvas::setSurfaceDataRangeBins(const int nbin) {
  n_range_sVar = nbin;
}
void CtiCanvas::setIsoDataRangeBins(const int nbin) {
  n_range_iVar = nbin;
}

// getters
void CtiCanvas::getX0(float _x0[3]) const{
  FOR_I3 _x0[i] = x0[i];
}
void CtiCanvas::getE0(float _e0[3]) const{
  FOR_I3 _e0[i] = e0[i];
}
void CtiCanvas::getE1(float _e1[3]) const{
  FOR_I3 _e1[i] = e1[i];
}
void CtiCanvas::getE2(float _e2[3]) const{
  FOR_I3 _e2[i] = e2[i];
}
int CtiCanvas::getNi() const{
  return ni;
}
int CtiCanvas::getNj() const{
  return nj;
}
float CtiCanvas::getWidth() const{
  return width;
}
//TODO DAP remove checks, assume valid i, j?
float CtiCanvas::getDepth(const int &i, const int &j) const{
  assert((i >= 0)&&(i < ni));
  assert((j >= 0)&&(j < nj));
  const int IJ = (j>>3)*nI + (i>>3);
  assert((IJ >= 0)&&(IJ < nI*nJ));
  assert(pbd[IJ]);

  const int ij = (j&7)*8+(i&7);
  assert((ij >= 0)&&(ij < 64));
  assert(pbd[IJ]->szone[ij] != 65535);
  return pbd[IJ]->sdepth[ij];
}

#undef CALC_PIXEL_I
#undef CALC_PIXEL_J
#undef CALC_PIXEL_DEPTH

#undef PBUF_SNORMAL0
#undef PBUF_SNORMAL1
#undef PBUF_SNORMAL2
#undef PBUF_SDEPTH
#undef PBUF_SDATA
#undef PBUF_MDEPTH
#undef PBUF_PNORMAL0
#undef PBUF_PNORMAL1
#undef PBUF_PNORMAL2
#undef PBUF_PDEPTH
#undef PBUF_PDATA
#undef PBUF_PAUX

#undef PBUF_SZONE
#undef PBUF_PZONE
#undef PBUF_BITS
#undef PBUF_I
#undef PBUF_J

#undef INTERNAL_DATA_BIT
#undef ISURFACE_BIT
#undef INTERNAL_BITS
#undef SURFACE_DATA_BIT
#undef SURFACE_BIT
#undef SURFACE_BITS
#undef MSURFACE_BIT
#undef MASK_BIT
#undef AUX_BIT
#undef MESH_BIT
#undef DELAY_BIT
