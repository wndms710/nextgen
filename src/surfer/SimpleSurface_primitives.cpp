#include "SimpleSurface.hpp"
#include "GeomUtils.hpp"
#include "GeodesicSphere.hpp"
#include "SplineStuff.hpp"

int SimpleSurface::addGridForFlaggedTris(const double diam,const double spacing,const double xc[3],const double normal[3],const string& zone_name) {

  // adds a turbulence grid of spheres within (inside and touching) the flagged tris...

  sp_flag.setLength(nsp);

  const double nmag = MAG(normal); assert(nmag > 0.0);
  const double unit_normal[3] = { normal[0]/nmag, normal[1]/nmag, normal[2]/nmag };
  for (int isp = 0; isp < nsp; ++isp) {
    const double dx[3] = DIFF(xsp[isp],xc);
    const double dist = DOT_PRODUCT(dx,unit_normal);
    if (dist > 0.5*diam) {
      sp_flag[isp] = 1;
    }
    else if (dist < -0.5*diam) {
      sp_flag[isp] = -1;
    }
    else {
      sp_flag[isp] = 0;
    }
  }

  int nst_active = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist]) {
      const int sp_sum = sp_flag[spost[ist][0]] + sp_flag[spost[ist][1]] + sp_flag[spost[ist][2]];
      if ((sp_sum != 3)&&(sp_sum != -3)) {
        ++nst_active;
      }
    }
  }

  if (nst_active == 0) {
    cout << "Error: no tris bounding requested grid plane" << endl;
    return -1;
  }

  // use the smallest value in normal to determine a direction not aligned with normal...
  double e1[3];
  if (fabs(normal[0]) < max(fabs(normal[1]),fabs(normal[2]))) {
    // use (1,0,0) cross normal for e1...
    e1[0] = 0.0;
    e1[1] = -normal[2];
    e1[2] = normal[1];
  }
  else if (fabs(normal[1]) < fabs(normal[2])) {
    // use (0,1,0) cross normal for e1...
    e1[0] = normal[2];
    e1[1] = 0.0;
    e1[2] = -normal[0];
  }
  else {
    // use (0,0,1) cross normal for e1...
    e1[0] = -normal[1];
    e1[1] = normal[0];
    e1[2] = 0.0;
  }
  const double e1_mag = MAG(e1);
  FOR_I3 e1[i] /= e1_mag;
  const double e0[3] = CROSS_PRODUCT(unit_normal,e1);

  double bb_min[2] = { 1.0E+20, 1.0E+20 };
  double bb_max[2] = { -1.0E+20, -1.0E+20 };
  int * active_tris = new int[nst_active];
  nst_active = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist]) {
      const int sp_sum = sp_flag[spost[ist][0]] + sp_flag[spost[ist][1]] + sp_flag[spost[ist][2]];
      if ((sp_sum != 3)&&(sp_sum != -3)) {
        active_tris[nst_active++] = ist;
        // also use this active tri to enlarge the bounding box...
        double dist[3];
        FOR_I3 {
          const int isp = spost[ist][i];
          const double dx[3] = DIFF(xsp[isp],xc);
          dist[i] = DOT_PRODUCT(dx,unit_normal);
        }
        FOR_I3 {
          const int isp0 = spost[ist][i];
          const int isp1 = spost[ist][(i+1)%2];
          if (dist[i]*dist[(i+1)%3] < 0.0) {
            double xp[3];
            FOR_J3 xp[j] = (dist[(i+1)%3]*xsp[isp0][j]-dist[i]*xsp[isp1][j])/(dist[(i+1)%3]-dist[i]);
            double x_primed = (xp[0]-xc[0])*e0[0] + (xp[1]-xc[1])*e0[1] + (xp[2]-xc[2])*e0[2];
            double y_primed = (xp[0]-xc[0])*e1[0] + (xp[1]-xc[1])*e1[1] + (xp[2]-xc[2])*e1[2];
            bb_min[0] = min(bb_min[0],x_primed);
            bb_max[0] = max(bb_max[0],x_primed);
            bb_min[1] = min(bb_min[1],y_primed);
            bb_max[1] = max(bb_max[1],y_primed);
          }
          else {
            if (dist[i] == 0.0) {
              double xp[3];
              FOR_J3 xp[j] = xsp[isp0][j];
              double x_primed = (xp[0]-xc[0])*e0[0] + (xp[1]-xc[1])*e0[1] + (xp[2]-xc[2])*e0[2];
              double y_primed = (xp[0]-xc[0])*e1[0] + (xp[1]-xc[1])*e1[1] + (xp[2]-xc[2])*e1[2];
              bb_min[0] = min(bb_min[0],x_primed);
              bb_max[0] = max(bb_max[0],x_primed);
              bb_min[1] = min(bb_min[1],y_primed);
              bb_max[1] = max(bb_max[1],y_primed);
            }
            if (dist[(i+1)%3] == 0.0) {
              double xp[3];
              FOR_J3 xp[j] = xsp[isp1][j];
              double x_primed = (xp[0]-xc[0])*e0[0] + (xp[1]-xc[1])*e0[1] + (xp[2]-xc[2])*e0[2];
              double y_primed = (xp[0]-xc[0])*e1[0] + (xp[1]-xc[1])*e1[1] + (xp[2]-xc[2])*e1[2];
              bb_min[0] = min(bb_min[0],x_primed);
              bb_max[0] = max(bb_max[0],x_primed);
              bb_min[1] = min(bb_min[1],y_primed);
              bb_max[1] = max(bb_max[1],y_primed);
            }
          }
        }
      }
    }
  }

  // here we could introduce an adt to store the tris and find the nearest, but
  // for now, do it brute force...

  // now figure out the grid of points...

  const int i0 = (int)ceil(bb_min[0]/spacing);
  const int i1 = (int)floor(bb_max[0]/spacing);
  const int j0 = (int)ceil(bb_min[1]/spacing);
  const int j1 = (int)floor(bb_max[1]/spacing);

  // check the n*m algorithm...
  if ((i1-i0+1)*(j1-j0+1)*nst_active > 100000) {
    cout << "Warning: this GRID seems very fine: " << (i1-i0+1)*(j1-j0+1) << " spheres, " << nst_active << " active tris (may take a while)." << endl;
  }

  //const double mpoint_factor = 1.0E-3;
  //vector<pair<int,int> > ijVec;
  vector<double> xNoIntersectVec;
  vector<double> xIntersectVec;
  const double diam_tol = 0.05;
  for (int i = i0; i <= i1; ++i) {
    for (int j = j0; j <= j1; ++j) {
      double xp[3];
      FOR_K3 xp[k] = xc[k] + double(i)*spacing*e0[k] + double(j)*spacing*e1[k];
      // for each point, find its closest tri and then determine if the closest
      // tri normal is pointing away from the line of sight. If it is, then
      // we are inside the requested tris...
      int ist_closest = -1;
      double d2_min;
      for (int ist_active = 0; ist_active < nst_active; ++ist_active) {
        const int ist = active_tris[ist_active];
        const double this_d2 = MiscUtils::getPointToTriDist2(xp,xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
        if ((ist_closest == -1)||(this_d2 < d2_min)) {
          ist_closest = ist;
          d2_min = this_d2;
        }
      }
      assert(ist_closest != -1);
      const double tri_normal[3] = TRI_NORMAL_2(xsp[spost[ist_closest][0]],xsp[spost[ist_closest][1]],xsp[spost[ist_closest][2]]);
      const double dx[3] = DIFF(xsp[spost[ist_closest][0]],xp);
      const double dp = DOT_PRODUCT(tri_normal,dx);
      if (dp > 0.0) {
        // dp > 0 means the center of the sphere is inside the geom. It may still intersect the geometry,
        // so treat these cases here...
        if (d2_min > diam*diam/4.0*(1.0+diam_tol)) {
          // this point is completely inside, so some of these points may intersect...
          xNoIntersectVec.push_back(xp[0]);
          xNoIntersectVec.push_back(xp[1]);
          xNoIntersectVec.push_back(xp[2]);
        }
        else if (d2_min > diam*diam/4.0) { // radius stupid!
          // move this point inside so it does not intersect the surface...
          const double delta = sqrt(diam*diam/4.0*(1.0+diam_tol)) - sqrt(d2_min);
          const double mag_n = MAG(tri_normal);
          xNoIntersectVec.push_back(xp[0]-delta*tri_normal[0]/mag_n);
          xNoIntersectVec.push_back(xp[1]-delta*tri_normal[1]/mag_n);
          xNoIntersectVec.push_back(xp[2]-delta*tri_normal[2]/mag_n);
        }
        else {
          // this sphere is touching, but make sure it is touching enough...
          const double delta = max(0.0,sqrt(d2_min) - sqrt(diam*diam/4.0*(1.0-diam_tol)));
          const double mag_n = MAG(tri_normal);
          xIntersectVec.push_back(xp[0]+delta*tri_normal[0]/mag_n);
          xIntersectVec.push_back(xp[1]+delta*tri_normal[1]/mag_n);
          xIntersectVec.push_back(xp[2]+delta*tri_normal[2]/mag_n);
        }
      }
      else if (d2_min < diam*diam/4.0*(1.0-diam_tol)) {
        xIntersectVec.push_back(xp[0]);
        xIntersectVec.push_back(xp[1]);
        xIntersectVec.push_back(xp[2]);
      }
    }
  }

  delete[] active_tris;

  const int n_sphere_edge = 4; // the decimation of the tri...
  int nsp0 = nsp;
  int nst0 = nst;

  const int nsp_inc = GeodesicSphere::getSphereNodeCount(n_sphere_edge);
  const int nst_inc = GeodesicSphere::getSphereTriCount(n_sphere_edge);
  nsp += (xNoIntersectVec.size()+xIntersectVec.size())/3*nsp_inc;
  nst += (xNoIntersectVec.size()+xIntersectVec.size())/3*nst_inc;
  growNspData(nsp,nsp0);
  growNstData(nst,nst0);

  // need a new (or existing) zone for the internal spheres...

  int izone;
  int nzn0 = zoneVec.size();
  for (izone = 0; izone < nzn0; ++izone) {
    if (zoneVec[izone].getName() == zone_name)
      break;
  }
  if (izone == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name));
    ++nsz;
  }

  // and the bool ones -- put these in the same zone as above...

  /*
    const string zone_name_bool = zone_name+"_bool";
    int izone_bool;
    nzn0 = zoneVec.size();
    for (izone_bool = 0; izone_bool < nzn0; ++izone_bool) {
    if (zoneVec[izone_bool].getName() == zone_name_bool)
    break;
    }
    if (izone_bool == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name_bool));
    ++nsz;
    }
  */
  int izone_bool = izone;

  // NOTE: addSphere assumes memory is available and returns incremented nsp0,nst0...
  for (int ii = 0,limit=xNoIntersectVec.size(); ii < limit; ii += 3) {
    double xp[3];
    xp[0] = xNoIntersectVec[ii];
    xp[1] = xNoIntersectVec[ii+1];
    xp[2] = xNoIntersectVec[ii+2];
    GeodesicSphere::addSphere(xsp+nsp0,spost+nst0,xp,0.5*diam,n_sphere_edge,true);
    for (int ist = nst0; ist < nst0+nst_inc; ++ist) {
      znost[ist] = izone;
      FOR_I3 spost[ist][i] += nsp0;
    }
    nsp0 += nsp_inc;
    nst0 += nst_inc;
  }

  for (int ii = 0,limit=xIntersectVec.size(); ii < limit; ii += 3) {
    double xp[3];
    xp[0] = xIntersectVec[ii];
    xp[1] = xIntersectVec[ii+1];
    xp[2] = xIntersectVec[ii+2];
    GeodesicSphere::addSphere(xsp+nsp0,spost+nst0,xp,0.5*diam,n_sphere_edge,true);
    for (int ist = nst0; ist < nst0+nst_inc; ++ist) {
      znost[ist] = izone_bool;
      FOR_I3 spost[ist][i] += nsp0;
    }
    nsp0 += nsp_inc;
    nst0 += nst_inc;
  }

  assert(nsp0 == nsp);
  assert(nst0 == nst);

  cout << " > added " << xNoIntersectVec.size()+xIntersectVec.size() << " turbulence grid spheres to zone \"" << zone_name << "\"" << endl;

  return 0;

}

int SimpleSurface::addSphere(const double xc0[3],const double r,const int n,const bool flip) {

  // current surface values
  int nsp0 = nsp;
  int nst0 = nst;

  // recall n contains the edge split size: e.g. 1,2,3 or some relatively low int
  nsp += GeodesicSphere::getSphereNodeCount(n);
  nst += GeodesicSphere::getSphereTriCount(n);
  growNspData(nsp,nsp0);
  growNstData(nst,nst0);

  GeodesicSphere::addSphere(xsp+nsp0,spost+nst0,xc0,r,n,flip);

  const string zone_name = "SIMPLE_SPHERE";
  int izone;
  int nzn0 = zoneVec.size();
  for (izone = 0; izone < nzn0; ++izone) {
    if (zoneVec[izone].getName() == zone_name)
      break;
  }
  if (izone == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name));
    ++nsz;
  }

  for (int ist = nst0; ist < nst; ++ist) {
    znost[ist] = izone;
    FOR_I3 spost[ist][i] += nsp0;
  }

  return 0;

}

int SimpleSurface::addHemisphere(const double xp[3],const double np[3],const double rp,const int ntheta,const bool flip) {

  int nsp0 = nsp;
  int nst0 = nst;

  int nsp_inc,nst_inc;
  GeomUtils::getHemisphereNodeAndTriCount(nsp_inc,nst_inc,ntheta);

  nsp += nsp_inc;
  nst += nst_inc;

  growNspData(nsp,nsp0);
  growNstData(nst,nst0);

  GeomUtils::addHemisphere(xsp+nsp0,spost+nst0,xp,np,rp,ntheta,flip);

  // finally the zone and also offset the spost if we weren't
  // the first tris...
  string zone_name = "HEMISPHERE";
  int izone,nzn;
  for (izone = 0,nzn=zoneVec.size(); izone < nzn; ++izone) {
    if (zoneVec[izone].getName() == zone_name)
      break;
  }
  if (izone == int(zoneVec.size())) {
    zoneVec.push_back(SurfaceZone(zone_name));
  }

  for (int ist = nst0; ist < nst; ++ist) {
    znost[ist] = izone;
    FOR_I3 spost[ist][i] += nsp0;
  }

  return 0;

}

int SimpleSurface::addPistonFF(const double dr0,const double dz0,const double dr1,const double dz1) {

  assert(0); // march 2019
  // probably use lifted surface and cylinder trimming...

  ensureTeost();

  // put the edges in a vec...

  int * sp_flag = new int[nsp];
  FOR_ISP sp_flag[isp] = -1;

  vector<pair<int,int> > edgeVec;
  double xcc[3] = { 0.0, 0.0, 0.0 };
  double normal[3] = { 0.0, 0.0, 0.0 };
  double d2min[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
  double d2max[3] = { 0.0, 0.0, 0.0 };
  double wgt_sum = 0.0;
  FOR_IST {
    FOR_I3 {
      if (isEdgeOpen(ist,i)) {
        edgeVec.push_back(pair<int,int>(ist,i));
        const int isp0 = spost[ist][i];
        const int isp1 = spost[ist][(i+1)%3];
        sp_flag[isp0] = 0;
        sp_flag[isp1] = 0;
        const double wgt = DIST(xsp[isp0],xsp[isp1]);
        FOR_I3 xcc[i] += wgt*(xsp[isp0][i]+xsp[isp1][i]);
        wgt_sum += wgt;
        const double cp[3] = CROSS_PRODUCT(xsp[isp0],xsp[isp1]);
        FOR_I3 normal[i] += cp[i];
        // also, assuming the piston is centered in either x,y,z, compute the min
        // and max radius^2...
        double xmid[3] = { 0.0, 0.0, 0.0 };
        FOR_I3 {
          xmid[i] += 0.5*(xsp[isp0][i]+xsp[isp1][i]);
          d2max[i] = max(d2max[i],xsp[isp0][(i+1)%3]*xsp[isp0][(i+1)%3]+xsp[isp0][(i+2)%3]*xsp[isp0][(i+2)%3]);
          d2max[i] = max(d2max[i],xsp[isp1][(i+1)%3]*xsp[isp1][(i+1)%3]+xsp[isp1][(i+2)%3]*xsp[isp1][(i+2)%3]);
        }
        FOR_I3 {
          d2min[i] = min(d2min[i],xmid[(i+1)%3]*xmid[(i+1)%3]+xmid[(i+2)%3]*xmid[(i+2)%3]);
        }
      }
    }
  }
  FOR_I3 xcc[i] /= 2.0*wgt_sum;
  double mag = MAG(normal);
  FOR_I3 normal[i] /= mag;

  cout << " > xcc, normal: " << COUT_VEC(xcc) << " " << COUT_VEC(normal) << endl;

  int nsp_edge = 0;
  FOR_ISP {
    if (sp_flag[isp] == 0)
      sp_flag[isp] = nsp_edge++;
    else
      assert(sp_flag[isp] == -1);
  }

  // if this is a loop, the nsp_edge size should be the same as the edgeVec.size()...
  assert(nsp_edge == int(edgeVec.size()));

  int id = -1;
  if (fabs(normal[0]) > 10.0*max(fabs(normal[1]),fabs(normal[2]))) {
    // this is an x-aligned piston...
    assert(0);
  }
  else if (fabs(normal[1]) >10.0*max(fabs(normal[0]),fabs(normal[2]))) {
    // this is an y-aligned piston...
    assert(0);
  }
  else if (fabs(normal[2]) >10.0*max(fabs(normal[0]),fabs(normal[1]))) {
    // this is an z-aligned piston...
    cout << " > assuming z-alignment" << endl;
    id = 2;
  }
  else {
    // no obvious alignment direction...
    assert(0);
  }

  cout << " > min/max radius of open edges: " << sqrt(d2min[id]) << " " << sqrt(d2max[id]) << endl;

  // the number of new tris will be edgeVec.size*5, and new nodes is 2*nsp_edge+1 (for the center)...
  const int nst0 = nst;
  nst += edgeVec.size()*5;
  growNstData(nst,nst0);

  const int nsp0 = nsp;
  nsp += nsp_edge*2+1;
  growNspData(nsp,nsp0);

  const string zone_name = "PISTON_FF";
  int izone;
  int nzn0 = zoneVec.size();
  for (izone = 0; izone < nzn0; ++izone) {
    if (zoneVec[izone].getName() == zone_name)
      break;
  }
  if (izone == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name));
    ++nsz;
  }

  // z-only for now...
  assert(id == 2);

  // nodes first...
  for (int isp = 0; isp < nsp0; ++isp) {
    if (sp_flag[isp] >= 0) {
      // figure out our "r"...
      const double rmag = sqrt(xsp[isp][0]*xsp[isp][0] + xsp[isp][1]*xsp[isp][1]);
      // first ring of points...
      xsp[nsp0+sp_flag[isp]][0] = xsp[isp][0] + dr0*xsp[isp][0]/rmag;
      xsp[nsp0+sp_flag[isp]][1] = xsp[isp][1] + dr0*xsp[isp][1]/rmag;
      xsp[nsp0+sp_flag[isp]][2] = xsp[isp][2] + dz0;
      // second ring of points...
      xsp[nsp0+nsp_edge+sp_flag[isp]][0] = xsp[isp][0] + dr1*xsp[isp][0]/rmag;
      xsp[nsp0+nsp_edge+sp_flag[isp]][1] = xsp[isp][1] + dr1*xsp[isp][1]/rmag;
      xsp[nsp0+nsp_edge+sp_flag[isp]][2] = xsp[isp][2] + dz1;
    }
  }
  assert(nsp0+nsp_edge*2 == nsp-1);
  xsp[nsp0+nsp_edge*2][0] = xcc[0];
  xsp[nsp0+nsp_edge*2][1] = xcc[1];
  xsp[nsp0+nsp_edge*2][2] = xcc[2]+dz1;

  // then new tris...
  for (int ied = 0,ned=edgeVec.size(); ied < ned; ++ied) {
    const int ist = edgeVec[ied].first;
    const int i = edgeVec[ied].second;
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    assert(sp_flag[isp0] >= 0);
    assert(sp_flag[isp1] >= 0);
    // first pair of tris...
    spost[nst0+ied*5  ][0] = isp1;
    spost[nst0+ied*5  ][1] = isp0;
    spost[nst0+ied*5  ][2] = nsp0+sp_flag[isp1];
    spost[nst0+ied*5+1][0] = isp0;
    spost[nst0+ied*5+1][1] = nsp0+sp_flag[isp0];
    spost[nst0+ied*5+1][2] = nsp0+sp_flag[isp1];
    // second pair of tris...
    spost[nst0+ied*5+2][0] = nsp0+sp_flag[isp1];
    spost[nst0+ied*5+2][1] = nsp0+sp_flag[isp0];
    spost[nst0+ied*5+2][2] = nsp0+nsp_edge+sp_flag[isp1];
    spost[nst0+ied*5+3][0] = nsp0+sp_flag[isp0];
    spost[nst0+ied*5+3][1] = nsp0+nsp_edge+sp_flag[isp0];
    spost[nst0+ied*5+3][2] = nsp0+nsp_edge+sp_flag[isp1];
    // third tri...
    spost[nst0+ied*5+4][0] = nsp0+nsp_edge+sp_flag[isp1];
    spost[nst0+ied*5+4][1] = nsp0+nsp_edge+sp_flag[isp0];
    spost[nst0+ied*5+4][2] = nsp0+nsp_edge*2;
    // znost...
    znost[nst0+ied*5  ] = izone;
    znost[nst0+ied*5+1] = izone;
    znost[nst0+ied*5+2] = izone;
    znost[nst0+ied*5+3] = izone;
    znost[nst0+ied*5+4] = izone;
  }

  delete[] sp_flag;

  return 0;

}

int SimpleSurface::addOffset(const double dn) {

  ensureTeost();

  // put the edges in a vec...

  vector<pair<int,int> > edgeVec;
  double xcc[3] = { 0.0, 0.0, 0.0 };
  double normal[3] = { 0.0, 0.0, 0.0 };
  double d2min[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
  double d2max[3] = { 0.0, 0.0, 0.0 };
  double wgt_sum = 0.0;
  FOR_IST {
    FOR_I3 {
      if (isEdgeOpen(ist,i)) {
        edgeVec.push_back(pair<int,int>(ist,i));
        const int isp0 = spost[ist][i];
        const int isp1 = spost[ist][(i+1)%3];
        const double wgt = DIST(xsp[isp0],xsp[isp1]);
        FOR_I3 xcc[i] += wgt*(xsp[isp0][i]+xsp[isp1][i]);
        wgt_sum += wgt;
        const double cp[3] = CROSS_PRODUCT(xsp[isp0],xsp[isp1]);
        FOR_I3 normal[i] += cp[i];
        // also, assuming the piston is centered in either x,y,z, compute the min
        // and max radius^2...
        double xmid[3] = { 0.0, 0.0, 0.0 };
        FOR_I3 {
          xmid[i] += 0.5*(xsp[isp0][i]+xsp[isp1][i]);
          d2max[i] = max(d2max[i],xsp[isp0][(i+1)%3]*xsp[isp0][(i+1)%3]+xsp[isp0][(i+2)%3]*xsp[isp0][(i+2)%3]);
          d2max[i] = max(d2max[i],xsp[isp1][(i+1)%3]*xsp[isp1][(i+1)%3]+xsp[isp1][(i+2)%3]*xsp[isp1][(i+2)%3]);
        }
        FOR_I3 {
          d2min[i] = min(d2min[i],xmid[(i+1)%3]*xmid[(i+1)%3]+xmid[(i+2)%3]*xmid[(i+2)%3]);
        }
      }
    }
  }
  FOR_I3 xcc[i] /= 2.0*wgt_sum;
  double mag = MAG(normal);
  FOR_I3 normal[i] /= mag;

  cout << " > xcc, normal: " << COUT_VEC(xcc) << " " << COUT_VEC(normal) << endl;

  // if this is a loop, the nsp_edge size should be the same as the edgeVec.size()...
  //assert(nsp_edge == edgeVec.size());

  // the number of new tris will be edgeVec.size*5, and new nodes is 2*nsp_edge+1 (for the center)...
  const int nst0 = nst;
  nst += nst0 + edgeVec.size()*2;
  growNstData(nst,nst0);

  const int nsp0 = nsp;
  nsp += nsp0;
  growNspData(nsp,nsp0);

  const string zone_name = "OFFSET";
  int izone;
  int nzn0 = zoneVec.size();
  for (izone = 0; izone < nzn0; ++izone) {
    if (zoneVec[izone].getName() == zone_name)
      break;
  }
  if (izone == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name));
    ++nsz;
  }

  // build the normals at the old nodes...
  double (*n_sp)[3] = new double[nsp0][3];
  for (int isp = 0; isp < nsp0; ++isp)
    FOR_I3 n_sp[isp][i] = 0.0;

  for (int ist = 0; ist < nst0; ++ist) {
    const double n[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    FOR_I3 {
      const int isp = spost[ist][i];
      FOR_J3 n_sp[isp][j] += n[j];
    }
  }

  // normalize...
  for (int isp = 0; isp < nsp0; ++isp) {
    const double mag = MAG(n_sp[isp]);
    assert(mag > 0.0);
    FOR_I3 n_sp[isp][i] /= mag;
  }

  // place new nodes at old + dn*n_sp...
  for (int isp = 0; isp < nsp0; ++isp) {
    const int isp_new = isp+nsp0;
    FOR_I3 xsp[isp_new][i] = xsp[isp][i] + n_sp[isp][i]*dn;
  }

  // then new tris: tri copies first: note orientation must be flipped...
  for (int ist = 0; ist < nst0; ++ist) {
    const int ist_new = ist+nst0;
    spost[ist_new][0] = spost[ist][0]+nsp0;
    spost[ist_new][1] = spost[ist][2]+nsp0;
    spost[ist_new][2] = spost[ist][1]+nsp0;
    znost[ist_new] = izone;
  }

  // then any tris from open edges...
  for (int ied = 0,ned=edgeVec.size(); ied < ned; ++ied) {
    const int ist = edgeVec[ied].first;
    const int i = edgeVec[ied].second;
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    // pair of tris...
    // 1...
    spost[nst0*2+ied*2][0] = isp0;
    spost[nst0*2+ied*2][1] = isp0+nsp0;
    spost[nst0*2+ied*2][2] = isp1;
    znost[nst0*2+ied*2]    = izone;
    // 2...
    spost[nst0*2+ied*2+1][0] = isp1;
    spost[nst0*2+ied*2+1][1] = isp0+nsp0;
    spost[nst0*2+ied*2+1][2] = isp1+nsp0;
    znost[nst0*2+ied*2+1]    = izone;
  }

  delete[] n_sp;
  return 0;

}

int SimpleSurface::constructPlane(int& iarg,Param * param) {
  int ierr = -1;
  const double pt[3] = {param->getDouble(iarg),param->getDouble(iarg+1),param->getDouble(iarg+2)};
  iarg += 3;
  const double norm[3] = {param->getDouble(iarg),param->getDouble(iarg+1),param->getDouble(iarg+2)};
  iarg += 3;

  int nx = 1;
  int ny = 1;
  double width,height = -1.0;

  while (iarg < param->size()) {
    const string tok = param->getString(iarg++);
    if (tok == "WIDTH") {
      width = param->getDouble(iarg++);
    }
    else if (tok == "HEIGHT") {
      height = param->getDouble(iarg++);
    }
    else if (tok == "NX") {
      nx = param->getInt(iarg++);
    }
    else if (tok == "NY") {
      ny = param->getInt(iarg++);
    }
    else {
      CWARN("unrecognized parameter \"" << tok << "\"; skipping");
    }
  }

  if (width<=0.0 || height <= 0.0) {
    CWARN("a positive WIDTH and HEIGHT must be specified");
    ierr = -1;
  }
  else ierr = addPlane(pt,norm,width,height,nx,ny);

  return ierr;
}

// initialize a plane surface by specifying the center point, normal axis, width, and height
int SimpleSurface::addPlane(const double xp[3],const double np[3],const double width,const double height,const int _nx, const int _ny) {

  // make sure enough points to discretize edges
  const int nx = max(_nx,1);
  const int ny = max(_ny,1);

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  nsp += (nx+1)*(ny+1);  // number of nodes is +1 number of facets in each dimension
  growNspData(nsp,ss_nsp0);

  double up[3] = {0.0,0.0,0.0};
  if (fabs(np[0]) < min(fabs(np[1]),fabs(np[2]))  || (np[0] == 0.0 && np[1] == 0.0)) up[0] = 1.0;
  else if (fabs(np[1]) < min(fabs(np[0]),fabs(np[2]))) up[1] = 1.0;
  else up[2] = 1.0;

  double e0[3] = CROSS_PRODUCT(up,np);
  NORMALIZE(e0);
  double e1[3] = CROSS_PRODUCT(np,e0);
  NORMALIZE(e1);

  const double delta_w = fabs(width)/nx;
  const double delta_h = fabs(height)/ny;

  double top_left[3];
  FOR_I3 top_left[i] = xp[i] - 0.5*fabs(height)*e0[i] - 0.5*fabs(width)*e1[i];

  int count = 0;
  for (int row=0, nrow=ny+1; row<nrow; ++row) {
    for (int col=0,ncol=nx+1; col<ncol; ++col) {
      ++count;
      FOR_I3 xsp[ss_nsp0 + (row*ncol) + col][i] = top_left[i] + row*delta_h*e0[i] + col*delta_w*e1[i];
    }
  }
  assert(count == (nsp-ss_nsp0));

  const string zone_name = "SIMPLE_PLANE";
  int new_zone;
  int nzn0 = zoneVec.size();
  for (new_zone = 0; new_zone < nzn0; ++new_zone) {
    if (zoneVec[new_zone].getName() == zone_name)
      break;
  }
  if (new_zone == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name));
    ++nsz;
  }

  nst += 2*nx*ny;  // number of new tris
  growNstData(nst,ss_nst0);

  count = 0;
  for (int row=0; row<ny; ++row) {
    for (int col=0; col<nx; ++col) {
      // below indices are relative to new nodes; need to index appropriately
      const int ul = row*(nx+1) + col;
      const int ur = ul+1;
      const int ll = ul+nx+1;
      const int lr = ll+1;

      // add 2 tris for this square
      spost[ss_nst0 + 2*(row*nx + col)][0] = ss_nsp0 + ul;
      spost[ss_nst0 + 2*(row*nx + col)][1] = ss_nsp0 + ll;
      spost[ss_nst0 + 2*(row*nx + col)][2] = ss_nsp0 + lr;

      spost[ss_nst0 + 2*(row*nx + col)+1][0] = ss_nsp0 + ul;
      spost[ss_nst0 + 2*(row*nx + col)+1][1] = ss_nsp0 + lr;
      spost[ss_nst0 + 2*(row*nx + col)+1][2] = ss_nsp0 + ur;

      count+=2;
    }
  }
  assert(count == nst-ss_nst0);

  for (int ist=ss_nst0; ist<nst; ++ist) {
    znost[ist] = new_zone;
  }

  // ====================================================
  // make sure all zones were set...
  for (int ist = 0; ist < nst; ++ist) {
    assert((znost[ist] >= 0)&&(znost[ist] < int(zoneVec.size())));
  }

  return 0;
}

int SimpleSurface::addCircle(const double xp[3],const double np[3],const double r,const int n) {
  assert(n>=3);
  assert(MAG(np)>0.0);

  double e1[3];
  double e2[3];
  MiscUtils::getBestE1E2FromE0(e1,e2,np);

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  nsp += n+1;  // number of nodes is n+1
  growNspData(nsp,ss_nsp0);

  FOR_I3 xsp[ss_nsp0][i] = xp[i];
  MiscUtils::createCirclePts(xsp,ss_nsp0+1,xp,np,r,n);

  // populate tri information
  const int new_zone = zoneVec.size();
  zoneVec.push_back(SurfaceZone("circle"));

  nst += n;
  growNstData(nst,ss_nst0);

  facetCircleToPoint(spost,znost,ss_nsp0+0,ss_nsp0+1,ss_nst0+0,new_zone,n,false);

  return 0;
}

int SimpleSurface::addAnnulus(const double xp[3],const double np[3],const double r0,const double r1,const int n) {
  assert(n>=3);
  assert(MAG(np)>0.0);

  double e1[3];
  double e2[3];
  MiscUtils::getBestE1E2FromE0(e1,e2,np);

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  nsp += 2*n;
  growNspData(nsp,ss_nsp0);

  MiscUtils::createCirclePts(xsp,ss_nsp0,xp,np,r0,n);
  MiscUtils::createCirclePts(xsp,ss_nsp0+n,xp,np,r1,n);

  // populate tri information
  const int new_zone = zoneVec.size();
  zoneVec.push_back(SurfaceZone("annulus"));

  nst += 2*n;
  growNstData(nst,ss_nst0);

  facetCircleToCircle(spost,znost,ss_nsp0+0,ss_nsp0+n,ss_nst0+0,new_zone,n,true,true);

  return 0;
}

// initialize a box surface from 2 corner points...
int SimpleSurface::addBox(const double x0[3],const double x1[3],const bool flip) {

  // make a box consisting of 6 sides, and 12 tris...
  // hex: standard node numbering is...
  //
  //        7------6
  //       /      /|
  //      4------5 |
  //      | 3    | 2
  //      |      |/
  //      0------1
  //

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  nsp += 8;
  growNspData(nsp,ss_nsp0);

  xsp[ss_nsp0+0][0] = x0[0]; xsp[ss_nsp0+0][1] = x0[1]; xsp[ss_nsp0+0][2] = x0[2];
  xsp[ss_nsp0+1][0] = x1[0]; xsp[ss_nsp0+1][1] = x0[1]; xsp[ss_nsp0+1][2] = x0[2];
  xsp[ss_nsp0+2][0] = x1[0]; xsp[ss_nsp0+2][1] = x1[1]; xsp[ss_nsp0+2][2] = x0[2];
  xsp[ss_nsp0+3][0] = x0[0]; xsp[ss_nsp0+3][1] = x1[1]; xsp[ss_nsp0+3][2] = x0[2];
  xsp[ss_nsp0+4][0] = x0[0]; xsp[ss_nsp0+4][1] = x0[1]; xsp[ss_nsp0+4][2] = x1[2];
  xsp[ss_nsp0+5][0] = x1[0]; xsp[ss_nsp0+5][1] = x0[1]; xsp[ss_nsp0+5][2] = x1[2];
  xsp[ss_nsp0+6][0] = x1[0]; xsp[ss_nsp0+6][1] = x1[1]; xsp[ss_nsp0+6][2] = x1[2];
  xsp[ss_nsp0+7][0] = x0[0]; xsp[ss_nsp0+7][1] = x1[1]; xsp[ss_nsp0+7][2] = x1[2];

  // do this a little more robustly in the future in case there are matching zone
  // names already there...
  const int zone_x0 = zoneVec.size();
  const int zone_x1 = zoneVec.size()+1;
  const int zone_y0 = zoneVec.size()+2;
  const int zone_y1 = zoneVec.size()+3;
  const int zone_z0 = zoneVec.size()+4;
  const int zone_z1 = zoneVec.size()+5;
  zoneVec.push_back(SurfaceZone("x0"));
  zoneVec.push_back(SurfaceZone("x1"));
  zoneVec.push_back(SurfaceZone("y0"));
  zoneVec.push_back(SurfaceZone("y1"));
  zoneVec.push_back(SurfaceZone("z0"));
  zoneVec.push_back(SurfaceZone("z1"));
  nsz += 6;

  nst += 12;
  growNstData(nst,ss_nst0);

  // x0...
  spost[ss_nst0+0][0] = 0; spost[ss_nst0+0][1] = 4; spost[ss_nst0+0][2] = 7; znost[ss_nst0+0] = zone_x0;
  spost[ss_nst0+1][0] = 0; spost[ss_nst0+1][1] = 7; spost[ss_nst0+1][2] = 3; znost[ss_nst0+1] = zone_x0;

  // x1...
  spost[ss_nst0+2][0] = 1; spost[ss_nst0+2][1] = 6; spost[ss_nst0+2][2] = 5; znost[ss_nst0+2] = zone_x1;
  spost[ss_nst0+3][0] = 1; spost[ss_nst0+3][1] = 2; spost[ss_nst0+3][2] = 6; znost[ss_nst0+3] = zone_x1;

  // y0...
  spost[ss_nst0+4][0] = 0; spost[ss_nst0+4][1] = 1; spost[ss_nst0+4][2] = 5; znost[ss_nst0+4] = zone_y0;
  spost[ss_nst0+5][0] = 0; spost[ss_nst0+5][1] = 5; spost[ss_nst0+5][2] = 4; znost[ss_nst0+5] = zone_y0;

  // y1...
  spost[ss_nst0+6][0] = 2; spost[ss_nst0+6][1] = 3; spost[ss_nst0+6][2] = 6; znost[ss_nst0+6] = zone_y1;
  spost[ss_nst0+7][0] = 3; spost[ss_nst0+7][1] = 7; spost[ss_nst0+7][2] = 6; znost[ss_nst0+7] = zone_y1;

  // z0...
  spost[ss_nst0+8][0] = 0; spost[ss_nst0+8][1] = 3; spost[ss_nst0+8][2] = 2; znost[ss_nst0+8] = zone_z0;
  spost[ss_nst0+9][0] = 0; spost[ss_nst0+9][1] = 2; spost[ss_nst0+9][2] = 1; znost[ss_nst0+9] = zone_z0;

  // z1...
  spost[ss_nst0+10][0] = 4; spost[ss_nst0+10][1] = 5; spost[ss_nst0+10][2] = 6; znost[ss_nst0+10] = zone_z1;
  spost[ss_nst0+11][0] = 4; spost[ss_nst0+11][1] = 6; spost[ss_nst0+11][2] = 7; znost[ss_nst0+11] = zone_z1;

  if (ss_nsp0 > 0) {
    // if adding surface, properly offset node indices, zones
    for (int ist = ss_nst0; ist < nst; ++ist) {
      FOR_I3 spost[ist][i] += ss_nsp0;
    }
  }

  // if flip requested...
  if (flip) {
    for (int ist = ss_nst0; ist < nst; ++ist) {
      const int tmp = spost[ist][1];
      spost[ist][1] = spost[ist][2];
      spost[ist][2] = tmp;
    }
  }

  return 0;

}

int SimpleSurface::addDisk(const double x0[3],const double x1[3],const double r,const int n,const bool b_flip) {
  const int ss_nst0 = nst;
  const int ierr = addTcone(x0,x1,r,r,n);
  if ((ierr == 0) && b_flip) {
    for (int ist = ss_nst0; ist < nst; ++ist) {
      const int tmp = spost[ist][1];
      spost[ist][1] = spost[ist][2];
      spost[ist][2] = tmp;
    }
  }
  return ierr;
}

int SimpleSurface::addTcone(const double x0[3],const double x1[3],const double rad0,const double rad1,const int n) {

  // n points around base...
  const double dx[3] = DIFF(x1,x0);
  const double dx_mag = MAG(dx);
  if (dx_mag <= 0.0) {
    CERR("cannot create TCONE with overlapping endpoints");
  }

  // absolute value ensures no singularity in the middle of the tcone as well as positive volume
  const double r0 = fabs(rad0);
  const double r1 = fabs(rad1);

  const double SMALL = 1E-12;
  if (r0 < SMALL || r1 < SMALL) {
    CERR("cannot create a TCONE with either radius = 0");
  }

  double axis[3];
  FOR_I3 axis[i] = x1[i] - x0[i];
  COUT2(" TCONE axis = " << COUT_VEC( axis ) );

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  nsp += 2 + 2*n;
  growNspData(nsp,ss_nsp0);

  // first and last node are cap centers
  FOR_I3 {
    xsp[ss_nsp0+0][i] = x0[i];
    xsp[nsp - 1][i] = x1[i];
  }

  // create nodes around each cap, cap0 first [1:n],
  // then cap1 [n+1:2n]
  // TODO: modify to GeomUtils...
  MiscUtils::createCirclePts(xsp,ss_nsp0+1,x0,axis,r0,n);
  MiscUtils::createCirclePts(xsp,ss_nsp0+n+1,x1,axis,r1,n);

  const string zone_name = "tcone";
  int new_zone;
  int nzn0 = zoneVec.size();
  for (new_zone = 0; new_zone < nzn0; ++new_zone) {
    if (zoneVec[new_zone].getName() == zone_name)
      break;
  }
  if (new_zone == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name));
    ++nsz;
  }

  nst += n*4;
  growNstData(nst,ss_nst0);

  // cap0, cap1
  // TODO: and here...
  facetCircleToPoint(spost,znost,ss_nsp0+0,ss_nsp0+1,ss_nst0+0,new_zone,n,false);
  facetCircleToPoint(spost,znost,nsp-1,ss_nsp0+n+1,ss_nst0+3*n,new_zone,n,true);

  // wall
  facetCircleToCircle(spost,znost,ss_nsp0+1,ss_nsp0+n+1,ss_nst0+n,new_zone,n,true,true);

  return 0;
}

int SimpleSurface::addAnnularTcone(const double x0[3],const double x1[3],const double rad00,const double rad01,const double rad10,const double rad11,const int n) {

  // n points around base...
  const double dx[3] = DIFF(x1,x0);
  const double dx_mag = MAG(dx);
  if (dx_mag <= 0.0) {
    CERR("cannot create ANNULAR_TCONE with overlapping endpoints");
  }

  // absolute value ensures no singularity in the middle of the tcone as well as positive volume
  const double r00 = fabs(rad00);
  const double r01 = fabs(rad01);
  const double r10 = fabs(rad10);
  const double r11 = fabs(rad11);

  const double SMALL = 1E-12;
  if (r00 < SMALL || r01 < SMALL || r10 < SMALL || r11 < SMALL) {
    CERR("cannot create an ANNUULAR_TCONE with any radius = 0");
  }

  double axis[3];
  FOR_I3 axis[i] = x1[i] - x0[i];
  COUT2(" ANNULAR_TCONE axis = " << COUT_VEC( axis ) );

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  nsp += 4*n;
  growNspData(nsp,ss_nsp0);


  // create nodes around each cap
  MiscUtils::createCirclePts(xsp,ss_nsp0+0*n,x0,axis,r00,n);
  MiscUtils::createCirclePts(xsp,ss_nsp0+1*n,x0,axis,r01,n);

  MiscUtils::createCirclePts(xsp,ss_nsp0+2*n,x1,axis,r10,n);
  MiscUtils::createCirclePts(xsp,ss_nsp0+3*n,x1,axis,r11,n);

  const string zone_name = "annular_tcone";
  int new_zone;
  int nzn0 = zoneVec.size();
  for (new_zone = 0; new_zone < nzn0; ++new_zone) {
    if (zoneVec[new_zone].getName() == zone_name)
      break;
  }
  if (new_zone == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name));
    ++nsz;
  }

  nst += n*8;
  growNstData(nst,ss_nst0);

  // cap0
  facetCircleToCircle(spost,znost,ss_nsp0+0*n,ss_nsp0+1*n,ss_nst0+0*n,new_zone,n,true,true);
  // r1
  facetCircleToCircle(spost,znost,ss_nsp0+1*n,ss_nsp0+3*n,ss_nst0+2*n,new_zone,n,true,true);
  // r0
  facetCircleToCircle(spost,znost,ss_nsp0+2*n,ss_nsp0+0*n,ss_nst0+4*n,new_zone,n,true,true);
  // cap1
  facetCircleToCircle(spost,znost,ss_nsp0+3*n,ss_nsp0+2*n,ss_nst0+6*n,new_zone,n,true,true);

  return 0;
}

int SimpleSurface::addRevolve(const double rAxis[3],const double x0[3],double (*ir_vals)[2],const int n_vals,const int _n_theta, bool capSurf, const double maxAR) {
  COUT1("SimpleSurface::initRevolve()");

  if (n_vals < 2) {
    CWARN("not enough points were specified in the profile; skipping");
    return -1;
  }
  if (_n_theta < 3) {
    CWARN("N_THETA of \"" << _n_theta << "\" is too small; defaulting to 3");
  }
  const int n_theta = max(3,_n_theta);

  double axis[3];
  FOR_I3 {
    axis[i] = rAxis[i];
  }
  NORMALIZE(axis); // make a unit vector
  COUT2(" > revolve axis = " << COUT_VEC( axis ) );

  vector< pair<double,double> > irVec;
  for (int v=0; v < n_vals; ++v) {
    irVec.push_back(pair<double,double> (ir_vals[v][0],ir_vals[v][1]));
  }

  for (int ixy = 1, limit=irVec.size()-1; ixy < limit; ++ixy) {
    if (irVec[ixy].second <= 0.0) {
      CWARN("cannot support profile points on the axis itself unless the starting or ending point; skipping");
      return -1;
    }
  }

  bool specifiedCap0 = false;
  bool specifiedCap1 = false;
  if (irVec[0].second == 0.0) {
    specifiedCap0 = true;
  } else {
    assert(irVec[0].second > 0.0);
  }

  if (irVec[irVec.size()-1].second == 0.0) {
    specifiedCap1 = true;
  } else {
    assert(irVec[irVec.size()-1].second > 0.0);
  }

  // if specified, add points within irVec to maintain certain aspect ratio
  // of facets. Since we use a fixed "n", we only look at adjustments introduced
  // by adding points in (I,R) space.
  if (maxAR != -1.0) {
    assert(maxAR >= 1.0);
    COUT2(" > aspect ratio limit: " << maxAR);

    if (capSurf) {
      // add point(s) to close the surface
      if (!specifiedCap0) {
        irVec.insert( irVec.begin(), pair<double,double>(irVec.front().first,0.0) );
        COUT2(" > adding starting point to surface at: (" << irVec.front().first << "," << irVec.front().second << ")");
      }

      if (!specifiedCap1) {
        irVec.push_back( pair<double,double>(irVec.back().first,0.0) );
        COUT2(" > adding closing point to surface at: (" << irVec.back().first << "," << irVec.back().second << ")");
      }
    }

    vector< pair<double,double> > irVec_new;
    irVec_new.push_back(irVec[0]);  //start from first point in IR profile

    int pos = 0;
    for (int i=1,limit=irVec.size(); i < limit; i++) {
      // COUT2("segment[" << i << "]");
      const double dI = irVec[i].first - irVec_new[pos].first;
      const double dR = irVec[i].second - irVec_new[pos].second;
      const double L0 = sqrt(dI*dI + dR*dR);
      double L1;
      if (i != limit-1) {
        L1 = 2*irVec[i].second*sin(M_PI/double(n_theta));
      }
      else {
        L1 = 2*irVec_new[pos].second*sin(M_PI/double(n_theta));
      }
      const double L2 = sqrt(L0*L0 + L1*L1);

      // compute aspect ratio based on ratio of circumradius/(2*inradius)
      // if ratio below criteria, add point(s) to profile
      const double s = 0.5*(L0+L1+L2);
      const double AR = L0*L1*L2 / (8.0*(s-L0)*(s-L1)*(s-L2));

      // COUT2(" > lengths L0,L1,L2,s: " << L0 << " " << L1 << " " << L2 << " " << s);
      // COUT2(" > computed AR: " << AR);

      if ((AR > maxAR) && (max(L0,L2) > L1)) {
        const int newPoints = (int)floor(AR/maxAR);
        // COUT2("    > inserting " << newPoints << " points");
        const double deltaI = dI / double(newPoints+1);
        const double deltaR = dR / double(newPoints+1);
        for (int p=0; p < newPoints; p++) {
          const double newI = irVec_new[pos].first + (p+1)*deltaI;
          const double newR = irVec_new[pos].second + (p+1)*deltaR;
          irVec_new.push_back(pair<double,double>(newI,newR));
        }
        pos += newPoints;
      }
      // add specified point
      irVec_new.push_back(irVec[i]);
      ++pos;
    }
    irVec.swap(irVec_new);

    if (capSurf) {
      // remove endpoints from irVec - used only for AR limiting
      if (!specifiedCap0) irVec.erase(irVec.begin());
      if (!specifiedCap1) irVec.pop_back();
    }
  }

  double x1[3];
  if (capSurf) {
    FOR_I3 x1[i] = x0[i] + irVec.back().first*axis[i];
  }

  if (specifiedCap0) {
    // remove endpoints from irVecs
    irVec.erase(irVec.begin());
  }

  if (specifiedCap1) {
    FOR_I3 x1[i] = x0[i] + irVec.back().first*axis[i];
    irVec.pop_back();
  }
  COUT2("x0:" << COUT_VEC(x0));
  COUT2("x1:" << COUT_VEC(x1));

  for (int i=0,limit=irVec.size(); i < limit; i++) {
    COUT2("[" << i << "]: " << irVec[i].first << " " << irVec[i].second);
  }
  COUT2(" > number of IR points prescribed: " << irVec.size());


  // Add new points
  const int ss_nsp0 = nsp;

  nsp += irVec.size()*n_theta;
  if (capSurf || specifiedCap0) nsp += 1;  // account for starting point
  if (capSurf || specifiedCap1) nsp += 1;  // account for ending point

  growNspData(nsp,ss_nsp0);

  int isp = ss_nsp0;
  if (capSurf || specifiedCap0) {
    FOR_I3 xsp[isp][i] = x0[i];
    isp++;  // cap0 endpoint
  }

  // place all irVec points
  for (int ixy=0, limit=irVec.size(); ixy < limit; ++ixy) {
    double center[3];
    FOR_I3 {
      center[i] = x0[i] + irVec[ixy].first*axis[i];
    }
    MiscUtils::createCirclePts(xsp,isp,center,axis,irVec[ixy].second,n_theta,false); // do NOT stagger the points
    isp += n_theta;
  }

  if (capSurf || specifiedCap1) {
    FOR_I3 xsp[isp][i] = x1[i];
    isp++;
  }

  assert(isp == nsp);

  // Add new tris
  const int ss_nst0 = nst;

  nst += 2*n_theta*(irVec.size()-1);
  if (capSurf || specifiedCap0) nst += n_theta;  // accounts for starting cap facets
  if (capSurf || specifiedCap1) nst += n_theta;  // accounts for ending cap facets

  growNstData(nst,ss_nst0);
  const int ss_nzn0 = zoneVec.size();

  COUT1(" > zones:");
  zoneVec.push_back(SurfaceZone("wall")); COUT1("    > \"wall\"");
  int cap0_idx=0;  // default write to wall zone
  int cap1_idx=0;  // default write to wall zone
  if (capSurf || specifiedCap0) {
    zoneVec.push_back(SurfaceZone("cap0")); COUT1("    > \"cap0\"");
    cap0_idx++;
    cap1_idx++;
  }
  if (capSurf || specifiedCap1) {
    zoneVec.push_back(SurfaceZone("cap1")); COUT1("    > \"cap1\"");
    cap1_idx++;
  }

  // place all surface tris
  int ist = ss_nst0;

  if (capSurf || specifiedCap0) {
    // cap0
    facetCircleToPoint(spost,znost,ss_nsp0,ss_nsp0+1,ist,cap0_idx,n_theta,false);
    ist += n_theta;
  }

  for (int ixy=0,limit=irVec.size()-1; ixy < limit; ++ixy) {
    // revolve
    int c0  = ss_nsp0 + ixy*n_theta;
    int c1  = ss_nsp0 + ixy*n_theta + n_theta;
    if (capSurf || specifiedCap0) {
      c0 += 1;
      c1 += 1;
    }
    facetCircleToCircle(spost,znost,c0,c1,ist,0,n_theta,true,true);
    ist += 2*n_theta;
  }

  if (capSurf || specifiedCap1) {
    // cap1
    facetCircleToPoint(spost,znost,nsp-1,nsp-n_theta-1,ist,cap1_idx,n_theta,true);
    ist += n_theta;
  }

  assert(ist == nst);

  if (ss_nst0 > 0) {
    // if adding surface, properly offset zone indices
    for (int ist = ss_nst0; ist < nst; ++ist) {
      // FOR_I3 spost[ist][i] += ss_nsp0;
      znost[ist] += ss_nzn0;
    }
  }

  COUT1(" > addRevolve done");
  // ====================================================
  // make sure all zones were set...
  for (int ist = 0; ist < nst; ++ist) {
    assert((znost[ist] >= 0)&&(znost[ist] < int(zoneVec.size())));
  }

  return 0;  // success
}

int SimpleSurface::addExtrudedProfile(const double x0[3], const double x1[3], double (*xyz_vals)[3], const int n_vals, const int _n_span, const bool closed) {

  COUT1("SimpleSurface::addExtrudedProfile()");
  COUT1(" > assumes profile is closed loop; if not this will produce errors currently");
  if (n_vals < 2) {
    CWARN("not enough points were specified in the profile; skipping");
    return -1;
  }

  int n_span = (_n_span+1);  // number of profile stations total, including start and end
  if (_n_span < 1) {
    CWARN("invalid specification of N_SPAN; defaulting to N_SPAN = 1");
    n_span = 2;
  }

  const double axis[3] = DIFF(x1,x0);
  if (MAG(axis) <= 0.0) {
    CWARN("extrude direction is not properly defined; skipping");
    return -1;
  }
  COUT2(" > extrusion direction: " << COUT_VEC(axis));

  // Add new points
  const int ss_nsp0 = nsp;
  nsp += n_vals*n_span;

  growNspData(nsp,ss_nsp0);

  // fill in points
  const double denom = 1.0 / (n_span-1);
  for (int ispan=0; ispan < n_span; ++ispan) {
    const double frac = ispan*denom;
    // double xc[3];
    // FOR_I3 xc[i] = x0[i] + frac*axis[i];
    for (int ii=0; ii<n_vals; ++ii) {
      FOR_I3 xsp[ss_nsp0+(ispan*n_vals)+ii][i] = x0[i] + xyz_vals[ii][i] + frac*axis[i];
    }
    // MiscUtils::create3DFrom2DPtAndNormal(xsp,ispan*n_vals,xc,axis,xy_vals,n_vals);
  }

  // create all tris
  const int ss_nst0 = nst;

  int _n_vals = n_vals;
  if (!closed) --_n_vals;  // last panel (between end and start node) doesn't get generated
  nst += 2*(n_span-1)*_n_vals;

  growNstData(nst,ss_nst0);
  const int ss_nzn0 = zoneVec.size();
  COUT1(" > zones:");
  zoneVec.push_back(SurfaceZone("extrusion")); COUT1("    > \"extrusion\"");

  int ist = ss_nst0;
  for (int ispan=0, limit = n_span-1; ispan < limit; ++ispan) {
    // start indices for circles-of-points nodes
    const int c0 = ss_nsp0 + ispan*n_vals;
    const int c1 = ss_nsp0 + (ispan+1)*n_vals;
    facetCircleToCircle(spost,znost,c0,c1,ist,ss_nzn0,n_vals,closed,true);
    ist += 2*_n_vals;
  }

  assert(ist == nst);

  COUT1(" > addExtrude done");
  // ====================================================
  // make sure all zones were set...
  for (int ist = 0; ist < nst; ++ist) {
    assert((znost[ist] >= 0)&&(znost[ist] < int(zoneVec.size())));
  }

  return 0;  // success

}

// create the facets between a circle of points and a single point

void SimpleSurface::facetCircleToPoint(int (* spost)[3], int * znost, const int indexPt, const int indexCircle, const int st0, const int zoneId, const int n, const bool flip) {
  // flip dictates node ordering (normal direction convention)
  int pt1_inc = 1;
  int pt2_inc = 0;
  if (flip) {
    pt1_inc = 0;
    pt2_inc = 1;
  }

  for (int i=0, limit=n; i < limit; ++i) {
    spost[st0+i][0] = indexPt;
    spost[st0+i][1] = indexCircle + (i+pt1_inc)%n;
    spost[st0+i][2] = indexCircle + (i+pt2_inc)%n;
    znost[st0+i] = zoneId;
  }
}


// create the facets between two circles of points, assuming they have the same n

void SimpleSurface::facetCircleToCircle(int (* spost)[3], int * znost, const int indexC0, const int indexC1, const int st0, const int zoneId, const int n, const bool closed, const bool flip) {
  // flip dictates node ordering (normal direction convention)
  int pt1_inc = 1;
  int pt2_inc = 0;
  if (flip) {
    pt1_inc = 0;
    pt2_inc = 1;
  }

  int limit = n;
  if (!closed) --limit;
  for (int i=0; i < limit; ++i) {
    spost[st0+(2*i)][0] = indexC0 + (i+pt1_inc)%n;
    spost[st0+(2*i)][1] = indexC0 + (i+pt2_inc)%n;
    spost[st0+(2*i)][2] = indexC1 + (i+pt1_inc)%n;
    znost[st0+(2*i)] = zoneId;

    spost[st0+(2*i + 1)][0] = indexC0 + (i+pt2_inc)%n;;
    spost[st0+(2*i + 1)][1] = indexC1 + (i+pt2_inc)%n;
    spost[st0+(2*i + 1)][2] = indexC1 + (i+pt1_inc)%n;
    znost[st0+(2*i + 1)] = zoneId;
  }
}

int SimpleSurface::addLofted(const vector<string>& profileVec,const int nr) {

  using namespace SplineStuff;

  const int npr = profileVec.size();
  assert(npr >= 2);
  vector<double> * dVec = new vector<double>[npr];
  int np = -1; // number of points in all profiles (should be the same)...
  for (int ipr = 0; ipr < npr; ++ipr) {
    //cout << " > reading BLADEGEN profile: " << profileVec[ipr] << endl;
    GeomUtils::readXYZ(dVec[ipr],profileVec[ipr]);
    assert(dVec[ipr].size()%3 == 0);
    const int this_np = dVec[ipr].size()/3;
    //cout << " > profile has np: " << this_np << endl;
    if (np == -1) {
      np = this_np;
    }
    else if (np != this_np) {
      cout << "Error: profiles with different point counts not supported" << endl;
      return -1;
    }
  }

  CubicSpline * cspline = new CubicSpline[np];
  double (*xpr)[3] = new double[npr][3];
  for (int ip = 0; ip < np; ++ip) {
    for (int ipr = 0; ipr < npr; ++ipr) {
      FOR_I3 xpr[ipr][i] = dVec[ipr][ip*3+i];
    }
    cspline[ip].init(xpr,npr);
  }
  delete[] xpr;

  // the number of tris is np*nr*2, and the number of nodes is np*(nr+1)...

  const int nsp = np*(nr+1);
  double (*xsp)[3] = new double[nsp][3];
  int isp = 0;
  for (int ip = 0; ip < np; ++ip) {
    for (int ir = 0; ir <= nr; ++ir) {
      cspline[ip].getX(xsp[isp],double(ir)/double(nr));
      ++isp;
    }
  }
  assert(isp == nsp);

  const int nst = np*nr*2;
  int (*spost)[3] = new int[nst][3];
  int ist = 0;
  for (int ip0 = 0; ip0 < np; ++ip0) {
    const int ip1 = (ip0+1)%np;
    for (int ir = 0; ir < nr; ++ir) {
      spost[ist][0] = ip0*(nr+1)+ir;
      spost[ist][1] = ip1*(nr+1)+ir;
      spost[ist][2] = ip0*(nr+1)+ir+1;
      spost[ist+1][0] = ip1*(nr+1)+ir;
      spost[ist+1][1] = ip1*(nr+1)+ir+1;
      spost[ist+1][2] = ip0*(nr+1)+ir+1;
      ist += 2;
    }
  }
  assert(ist == nst);

  // for now, just output an sbin...

  GeomUtils::writeTecplot("nasa37.dat",spost,nst,xsp);
  GeomUtils::writeSbin("nasa37.sbin",spost,nst,xsp);

  cout << " > looks good: about to throw" << endl;
  getchar();
  throw(0);

}

int SimpleSurface::addExtrudeSpline(vector<double>& dVec,const double x0[3],const double x1[3],const int nspan) {

  using namespace SplineStuff;

  // the dVec contains x0 and x1 locations...

  assert(dVec.size()%2 == 0);
  const int np = dVec.size()/2;

  // splines are 3D, so use z = 0...

  double (*xp)[3] = new double[np][3];
  for (int ip = 0; ip < np; ++ip) {
    xp[ip][0] = dVec[ip*2  ];
    xp[ip][1] = dVec[ip*2+1];
    xp[ip][2] = 0.0; // 2d spline would be nice
  }
  CubicSpline cspline;
  cspline.init(xp,np,true);
  delete[] xp;
  const double length = cspline.getLength();

  cout << " > length: " << length << endl;

  // take a look...
  {
    FILE * fp = fopen("spline.dat","w");
    const int n = 2000;
    for (int i = 0; i <= n; ++i) {
      double xr[3];
      cspline.getX(xr,double(i)/double(n));
      fprintf(fp,"%18.15le %18.15le\n",xr[0],xr[1]);
    }
    fclose(fp);
  }

  // to determine nx, we use aspect ratio considerations...

  const double dx[3] = DIFF(x1,x0);
  const double dx_mag = MAG(dx);
  assert(dx_mag > 0.0);
  double e1[3],e2[3];
  MiscUtils::getBestE1E2FromE0(e1,e2,dx);

  const int nx = max(1,int(length/dx_mag*double(nspan)));
  cout << " > nspan: " << nspan << ", nx (computed from isotropic considerations): " << nx << endl;

  int nsp = nx*(nspan+1);
  int nst = 2*nspan*nx;

  cout << " > nsp: " << nsp << ", nst: " << nst << endl;

  double (*xsp)[3] = new double[nsp][3];
  int isp = 0;
  for (int i = 0; i < nx; ++i) {
    double x[3];
    cspline.getX(x,double(i)/double(nx));
    for (int j = 0; j <= nspan; ++j) {
      FOR_K3 xsp[isp][k] = x0[k] + double(j)/double(nspan)*dx[k] + x[0]*e1[k] + x[1]*e2[k];
      ++isp;
    }
  }
  assert(isp == nsp);

  int ist = 0;
  int (*spost)[3] = new int[nst][3];
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nspan; ++j) {
      // this patch has 4 points:
      const int isp0 = i*(nspan+1)+j;
      const int isp1 = ((i+1)%nx)*(nspan+1)+j;
      const int isp2 = ((i+1)%nx)*(nspan+1)+(j+1);
      const int isp3 = i*(nspan+1)+(j+1);
      spost[ist][0] = isp0;
      spost[ist][1] = isp1;
      spost[ist][2] = isp2;
      spost[ist+1][0] = isp0;
      spost[ist+1][1] = isp2;
      spost[ist+1][2] = isp3;
      ist += 2;
    }
  }
  assert(ist == nst);

  //GeomUtils::writePtsTecplot("pts.dat",xsp,nsp);
  GeomUtils::writeSbin("extrude_spline.sbin",spost,nst,xsp);
  GeomUtils::writeTecplot("extrude_spline.dat",spost,nst,xsp);

  delete[] spost;
  delete[] xsp;

  cout << "EXTRUDE_SPLINE just wrote extrude_spline.sbin -- finish this off some day" << endl;

  return -1;

}

int SimpleSurface::addRevolveSpline(vector<double>& dVec,const double axis[3],const double origin[3],const int ntheta) {

  using namespace SplineStuff;

  // the dVec contains axial distance, radial distance...

  assert(dVec.size()%2 == 0);
  const int np = dVec.size()/2;

  // splines are 3D, so use z = 0...

  double r_avg = 0.0;
  double (*xp)[3] = new double[np][3];
  for (int ip = 0; ip < np; ++ip) {
    xp[ip][0] = dVec[ip*2  ];
    xp[ip][1] = dVec[ip*2+1];
    r_avg += xp[ip][1];
    xp[ip][2] = 0.0; // 2d spline would be nice
  }
  CubicSpline cspline;
  cspline.init(xp,np);
  delete[] xp;
  r_avg /= double(np);
  const double length = cspline.getLength();

  cout << " > r_avg: " << r_avg << " length: " << length << endl;

  //// take a look...
  //{
  //  FILE * fp = fopen("spline.dat","w");
  //  const int n = 2000;
  //  for (int i = 0; i <= n; ++i) {
  //    double xr[3];
  //    cspline.getX(xr,double(i)/double(n));
  //    fprintf(fp,"%18.15le %18.15le\n",xr[0],xr[1]);
  //  }
  //  fclose(fp);
  //}

  const double axis_mag = MAG(axis);
  assert(axis_mag > 0.0);
  double e1[3],e2[3];
  MiscUtils::getBestE1E2FromE0(e1,e2,axis);

  // to determine nx, we use aspect ratio considerations:

  const int nx = max(1,int(length/(2.0*M_PI*r_avg/double(ntheta))));
  cout << " > ntheta: " << ntheta << ", nx (computed from isotropic considerations): " << nx << endl;

  // current surface values
  int nsp0 = nsp;
  int nst0 = nst;

  nsp = (nx+1)*ntheta;
  nst = 2*ntheta*nx;

  cout << " > nsp: " << nsp << ", nst: " << nst << endl;

  growNspData(nsp,nsp0);
  growNstData(nst,nst0);

  //double (*xsp)[3] = new double[nsp][3];
  int isp = nsp0;
  for (int i = 0; i <= nx; ++i) {
    double xr[3];
    cspline.getX(xr,double(i)/double(nx));
    for (int j = 0; j < ntheta; ++j) {
      const double theta = M_PI*2.0*double(j)/double(ntheta);
      FOR_K3 xsp[isp][k] = origin[k] + xr[0]*axis[k] + xr[1]*(cos(theta)*e1[k] + sin(theta)*e2[k]);
      ++isp;
    }
  }
  assert(isp == nsp);

  int ist = nst0;
  //int (*spost)[3] = new int[nst][3];
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ntheta; ++j) {
      // this patch has 4 points:
      const int isp0 = i*ntheta+j;
      const int isp1 = (i+1)*ntheta+j;
      const int isp2 = (i+1)*ntheta+((j+1)%ntheta);
      const int isp3 = i*ntheta+((j+1)%ntheta);
      spost[ist][0] = isp0+nsp0;
      spost[ist][1] = isp1+nsp0;
      spost[ist][2] = isp2+nsp0;
      spost[ist+1][0] = isp0+nsp0;
      spost[ist+1][1] = isp2+nsp0;
      spost[ist+1][2] = isp3+nsp0;
      ist += 2;
    }
  }
  assert(ist == nst);

  //GeomUtils::writePtsTecplot("pts.dat",xsp,nsp);
  //GeomUtils::writeSbin("revolve_spline.sbin",spost,nst,xsp);
  //GeomUtils::writeTecplot("revolve_spline.dat",spost,nst,xsp);

  //delete[] spost;
  //delete[] xsp;

  //cout << "REVOLVE_SPLINE just wrote revolve_spline.sbin -- finish this off some day" << endl;
  //return -1;

  const string zone_name = "REVOLVE_SPLINE";
  int izone;
  int nzn0 = zoneVec.size();
  for (izone = 0; izone < nzn0; ++izone) {
    if (zoneVec[izone].getName() == zone_name)
      break;
  }
  if (izone == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name));
    ++nsz;
  }

  for (int ist = nst0; ist < nst; ++ist) {
    znost[ist] = izone;
    //FOR_I3 spost[ist][i] += nsp0; //offset applied above
  }

  // ====================================================
  // make sure all zones were set...
  for (int ist = 0; ist < nst; ++ist) {
    assert((znost[ist] >= 0)&&(znost[ist] < int(zoneVec.size())));
  }

}

void SimpleSurface::interfaceFlaggedZones(const double xc[3],const double axis[3]) {
  
  cout << "SimpleSurface::interfaceFlaggedZones(): xc: " << COUT_VEC(xc) << " axis: " << COUT_VEC(axis) << endl;

  const double axis_mag = MAG(axis);
  double unit_axis[3] = { axis[0]/axis_mag, axis[1]/axis_mag, axis[2]/axis_mag };
  
  ensureTeost();

  int * sp_flag = new int[nsp];
  double (*dxsp)[3] = new double[nsp][3];
  for (int isp = 0; isp < nsp; ++isp) {
    sp_flag[isp] = -1;
    FOR_I3 dxsp[isp][i] = 0.0;
  }
  for (int ist = 0; ist < nst; ++ist) {
    if (zoneVec[znost[ist]].flag == 1) {
      FOR_I3 {
	int ist_nbr;
	bool b_nbr = getAlignedTriNbr(ist_nbr,ist,i);
	assert(b_nbr);
	assert(ist_nbr != ist);
	if (zoneVec[znost[ist_nbr]].flag == 0) {
	  // this is an edge of the flagged region: make sure it is as we expect...
	  sp_flag[spost[ist][i]] = spost[ist][(i+1)%3];
	}
      }
    }
  }
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] >= 0) {
      // this is the start of a loop...
      double xsum = 0.0;
      double rsum = 0.0;
      double wsum = 0.0;
      double xcsum[3] = { 0.0, 0.0, 0.0 };
      double nsum[3] = { 0.0, 0.0, 0.0 };
      const int isp_start = isp;
      double dx[3] = DIFF(xsp[isp_start],xc);
      double x_start = DOT_PRODUCT(dx,unit_axis);
      double r_start = sqrt(max(0.0,DOT_PRODUCT(dx,dx)-x_start*x_start));
      int isp_next = isp_start;
      double x_next = x_start;
      double r_next = r_start;
      int count = 0;
      while (1) {
	++count;
	const int isp_prev = isp_next;
	const double x_prev = x_next;
	const double r_prev = r_next;
	isp_next = sp_flag[isp_prev];
	assert(isp_next >= 0);
	double dx[3] = DIFF(xsp[isp_next],xc);
	x_next = DOT_PRODUCT(dx,unit_axis);
	r_next = sqrt(max(0.0,DOT_PRODUCT(dx,dx)-x_next*x_next));
	const double wgt = DIST(xsp[isp_prev],xsp[isp_next]);
	xsum += wgt*(x_prev+x_next);
	rsum += wgt*(r_prev+r_next);
	wsum += wgt;
	// and the center and normal...
	FOR_I3 xcsum[i] += wgt*(xsp[isp_prev][i]+xsp[isp_next][i]);
	double n2[3] = TRI_NORMAL_2(xsp[isp_start],xsp[isp_prev],xsp[isp_next]);
	FOR_I3 nsum[i] += n2[i];
	if (isp_next == isp_start)
	  break;
      }
      xsum /= 2.0*wsum;
      rsum /= 2.0*wsum;
      FOR_I3 xcsum[i] /= 2.0*wsum;
      double nmag = MAG(nsum);
      FOR_I3 nsum[i] /= nmag;
      // loop again, setting dxsp to snap to exact x,r...
      double d2_max = 0.0;
      double x_max[3];
      isp_next = isp_start;
      while (1) {
	const int isp_prev = isp_next;
	isp_next = sp_flag[isp_prev];
	assert(isp_next >= 0);
	sp_flag[isp_prev] = -2;
	double dx[3] = DIFF(xsp[isp_next],xc);
	x_next = DOT_PRODUCT(dx,unit_axis);
	r_next = sqrt(max(0.0,DOT_PRODUCT(dx,dx)-x_next*x_next));
	double dx_next = x_next-xsum;
	double dr_next = r_next-rsum;
	// the unit_axis gives the dx correction...
	FOR_I3 dxsp[isp_next][i] = -dx_next*unit_axis[i];
	// the r is...
	//double dr[3]; FOR_I3 dr[i] = dx[i] - x_next*unit_axis[i];
	FOR_I3 dxsp[isp_next][i] -= dr_next*(dx[i] - x_next*unit_axis[i])/r_next;
	const double d2 = DOT_PRODUCT(dxsp[isp],dxsp[isp]);
	if (d2 > d2_max) {
	  d2_max = d2;
	  FOR_I3 x_max[i] = xsp[isp][i];
	}
	if (isp_next == isp_start)
	  break;
      }
      cout << " > open edge loop with " << count << 
	" edges: avg r: " << rsum << " n: " << COUT_VEC(nsum) << " xc: " << COUT_VEC(xcsum) << 
	" max dx: " << sqrt(d2_max) << " at x: " << COUT_VEC(x_max) << endl;
    }
  }
  
  // now we have dxsp where sp_flag is -2, zero where sp_flag is -1
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] == -1) {
      assert(DOT_PRODUCT(dxsp[isp],dxsp[isp]) == 0.0);
    }
    else {
      assert(sp_flag[isp] == -2);
      assert(DOT_PRODUCT(dxsp[isp],dxsp[isp]) >= 0.0); // could still be zero
    }
  }

  // do some jacobi smoothing iterations...
  // TODO: must be a way to accelerate this...
  // e.g. flag all the tris that have nodes that can possibly 
  // be changed and just work on those...
  double rate = 0.2;
  double *wgt = new double[nsp];
  double (*sum)[3] = new double[nsp][3];
  int niter = 15;
  cout << " > " << niter << " smoothing iterations";
  for (int iter = 0; iter < niter; ++iter) {
    cout << ".";
    cout.flush();
    for (int isp = 0; isp < nsp; ++isp) {
      if (sp_flag[isp] == -1) {
	wgt[isp] = 0.0;
	FOR_I3 sum[isp][i] = 0.0;
      }
    }
    double mag[3];
    double length[3];
    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
	mag[i] = MAG(dxsp[spost[ist][i]]);
	length[i] = DIST(xsp[spost[ist][i]],xsp[spost[ist][(i+1)%3]]);
      }
      FOR_I3 {
	const int isp = spost[ist][i];
	if (sp_flag[isp] == -1) {
	  // we linearly reduce the nbr magnitude based on distance and a rate. When the
	  // rate is 1, the magnitude reduces to zero after length == mag. For 
	  // lower rates (recommended), the magnitude reduces to zero after length*rate == mag.
	  if (length[i]*rate < mag[(i+1)%3]) {
	    double new_mag = mag[(i+1)%3]-length[i]*rate + 1.0E-10; assert(new_mag > 0.0);
	    wgt[isp] += new_mag;
	    FOR_J3 sum[isp][j] += dxsp[spost[ist][(i+1)%3]][j]*new_mag*new_mag/mag[(i+1)%3];
	  }
	  if (length[(i+2)%3]*rate < mag[(i+2)%3]) {
	    double new_mag = mag[(i+2)%3]-length[(i+2)%3]*rate + 1.0E-10; assert(new_mag > 0.0);
	    wgt[isp] += new_mag;
	    FOR_J3 sum[isp][j] += dxsp[spost[ist][(i+2)%3]][j]*new_mag*new_mag/mag[(i+2)%3];
	  }
	}
      }
    }
    for (int isp = 0; isp < nsp; ++isp) {
      if ((sp_flag[isp] == -1)&&(wgt[isp] > 0.0)) {
	FOR_I3 dxsp[isp][i] = sum[isp][i]/wgt[isp];
      }
    }
  }
  delete[] wgt;
  delete[] sum;
  
  cout << "done" << endl;
  
  // finally, add the perturbation...
  for (int isp = 0; isp < nsp; ++isp) {
    FOR_I3 xsp[isp][i] += dxsp[isp][i];
  }
  
  delete[] dxsp;
  delete[] sp_flag;
  
}

class SegmentPoint {
public:
  double xp[2];
  int prev;
  int next;
};

void getPointEdgeDist2(double &d2,double &dd2dx,double &dd2dy,
		       const double xp,const double yp,const double xe0,const double ye0,const double xe1,const double ye1) {
  double dxp = xp-xe0;
  double dyp = yp-ye0;
  const double dxe = xe1-xe0;
  const double dye = ye1-ye0;
  const double dp = dxp*dxe+dyp*dye;
  if (dp <= 0.0) {
    d2 = dxp*dxp + dyp*dyp;
    dd2dx = 2.0*dxp;
    dd2dy = 2.0*dyp;
  }
  else {
    const double dxe2 = dxe*dxe+dye*dye;
    if (dp >= dxe2) {
      // we are closest to the second point v1...
      dxp = xp-xe1;
      dyp = yp-ye1;
      d2 = dxp*dxp + dyp*dyp;
      dd2dx = 2.0*dxp;
      dd2dy = 2.0*dyp;
    }
    else {
      // d2 can get slightly negative due to machine precision errors...
      d2 = max(0.0, dxp*dxp + dyp*dyp - dp*dp/dxe2);
      dd2dx = 2.0*dxp - 2.0*dp*dxe/dxe2;
      dd2dy = 2.0*dyp - 2.0*dp*dye/dxe2;
    }
  }
}

// this fun

#define ANGLEFUNC(I) ((1-(((I)&1)<<1))*((I)>>1))

// i: 0 ANGLEFUNC(i): 0
// i: 1 ANGLEFUNC(i): 0
// i: 2 ANGLEFUNC(i): 1
// i: 3 ANGLEFUNC(i): -1
// i: 4 ANGLEFUNC(i): 2
// i: 5 ANGLEFUNC(i): -2
// i: 6 ANGLEFUNC(i): 3
// i: 7 ANGLEFUNC(i): -3
// i: 8 ANGLEFUNC(i): 4
// i: 9 ANGLEFUNC(i): -4
// i: 10 ANGLEFUNC(i): 5
// i: 11 ANGLEFUNC(i): -5
// i: 12 ANGLEFUNC(i): 6
// i: 13 ANGLEFUNC(i): -6
// i: 14 ANGLEFUNC(i): 7
// i: 15 ANGLEFUNC(i): -7
// i: 16 ANGLEFUNC(i): 8
// i: 17 ANGLEFUNC(i): -8
// i: 18 ANGLEFUNC(i): 9
// i: 19 ANGLEFUNC(i): -9

void SimpleSurface::makeFF(const double xc[3],const double axis[3],const double dn,const double ds) {

  cout << "SimpleSurface::makeFF(): xc: " << COUT_VEC(xc) << " axis: " << COUT_VEC(axis) << endl;
  
  const double axis_mag = MAG(axis);
  double unit_axis[3] = { axis[0]/axis_mag, axis[1]/axis_mag, axis[2]/axis_mag };
  
  ensureTeost();
  
  int *sp_flag = new int[nsp];
  for (int isp = 0; isp < nsp; ++isp) {
    sp_flag[isp] = -1;
  }

  int *st_flag = new int[nst];
  for (int ist = 0; ist < nst; ++ist) {
    st_flag[ist] = -1;
  }
  
  // the 2D projected coordinates are...
  double (*xp_sp)[2] = new double[nsp][2]; 
  for (int isp = 0; isp < nsp; ++isp) {
    const double dx[3] = DIFF(xsp[isp],xc);
    xp_sp[isp][1] = DOT_PRODUCT(dx,unit_axis);
    xp_sp[isp][0] = sqrt(max(0.0,DOT_PRODUCT(dx,dx)-xp_sp[isp][1]*xp_sp[isp][1]));
  }
  
  // the 2D projected grid (xp,yp) is based on yp = dx.dot.unit_axis, xp = r 
  double xp_max = 0.0;
  double yp_min = 1.0E+20;
  double yp_max = -1.0E+20;
  
  bool b_axis0 = false;
  double yp_axis0 = 0.0;
  
  bool b_axis1 = false;
  double yp_axis1 = 0.0;

  for (int ist = 0; ist < nst; ++ist) {
    if (zoneVec[znost[ist]].flag == 1) {
      assert(st_flag[ist] == -1);
      st_flag[ist] = 1;
      double A[3];
      FOR_I3 {
	// start by expanding the 2D bounding box based on this tri's nodes...
	const int isp = spost[ist][i];
	xp_max = max(xp_max,xp_sp[isp][0]);
	yp_min = min(yp_min,xp_sp[isp][1]);
	yp_max = max(yp_max,xp_sp[isp][1]);
	// next, we are going to check if this tri intersects the axis...
	const double dx[3] = DIFF(xsp[isp],xc);
	const int isp_next = spost[ist][(i+1)%3];
	const double dx_next[3] = DIFF(xsp[isp_next],xc);
	const double cp[3] = CROSS_PRODUCT(dx,dx_next);
	A[(i+2)%3] = DOT_PRODUCT(cp,unit_axis); // populate the barycentric weight of the opposite node
	// next flag any nodes that touch open edges: set their sp_flag to the NEXT
	// node along the edge... 
	int ist_nbr;
	bool b_nbr = getAlignedTriNbr(ist_nbr,ist,i);
	assert(b_nbr);
	assert(ist_nbr != ist);
	if (zoneVec[znost[ist_nbr]].flag == 0) {
	  // this is an edge of the flagged region: make sure it is as we expect...
	  sp_flag[isp] = ist; // isp_next; use the tri instead of the isp_next, which we can get from the tri
	}
      }
      if ((A[0] >= 0.0)&&(A[1] >= 0.0)&&(A[2] >= 0.0)&&(A[0]+A[1]+A[2] > 0.0)) {
	cout << "XXXXXX got positive projected area!" << endl;
	assert(0);
      }
      else if ((A[0] <= 0.0)&&(A[1] <= 0.0)&&(A[2] <= 0.0)&&(A[0]+A[1]+A[2] < 0.0)) {
	// turn off this tri for later consideration...
	st_flag[ist] = 0;
	// use the A's as baricentric weights to store the yp intersection...
	double yp_avg = 0.0;
	FOR_I3 {
	  // start by expanding the 2D bounding box based on this tri's nodes...
	  const int isp = spost[ist][i];
	  yp_avg += A[i]*xp_sp[isp][1];
	}
	yp_avg /= A[0]+A[1]+A[2];
	if (!b_axis1) {
	  yp_axis1 = yp_avg;
	  b_axis1 = true;

	  cout << "AXIS1 " << 0.0 << " " << yp_avg << endl;
	  cout << "AXIS1 " << 0.0 << " " << yp_avg+dn << endl;
	  
	}
	else {
	  cout << " got multiple axis crossings: " << yp_axis1 << " " << yp_avg << endl;
	  assert(0); // figure this out later
  	}
      }
    }
  }
  
  // report what we found...
  if (b_axis0) {
    cout << " > got yp_axis0: " << yp_axis0 << endl;
    yp_min = min(yp_min,yp_axis0);
  }
  if (b_axis1) {
    cout << " > got yp_axis1: " << yp_axis1 << endl;
    yp_max = min(yp_max,yp_axis1);
  }
  cout << " > projected grid size: xp_max: " << xp_max << " yp_min, yp_max: " << yp_min << " " << yp_max << endl; 
  
  bool b_loop0 = false;
  double xp_loop0 = 0.0;
  double yp_loop0 = 0.0;

  bool b_loop1 = false;
  double xp_loop1 = 0.0;
  double yp_loop1 = 0.0;
  
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] >= 0) {
      // this is the start of a loop...
      double xp_sum = 0.0;
      double yp_sum = 0.0;
      double w_sum = 0.0;
      double xp_min = 1.0E+20;
      double xp_max = 0.0;
      double xc_sum[3] = { 0.0, 0.0, 0.0 };
      double normal_sum[3] = { 0.0, 0.0, 0.0 };
      double dxp_m_sum = 0.0;
      double dyp_m_sum = 0.0;
      const int isp_start = isp;
      int isp_next = isp_start;
      int count = 0;
      while (1) {
	++count;
	const int isp_prev = isp_next;
	//isp_next = sp_flag[isp_prev];
	const int ist = sp_flag[isp_prev];
	assert(zoneVec[znost[ist]].flag == 1); // this is the moving side
	int i;
	for (i = 0; i < 3; ++i) {
	  if (spost[ist][i] == isp_prev)
	    break;
	}
	assert(i < 3);
	isp_next = spost[ist][(i+1)%3];
	assert(isp_next >= 0);
	int ist_nbr;
	bool b_nbr = getAlignedTriNbr(ist_nbr,ist,i);
	assert(b_nbr);
	assert(ist_nbr != ist);
	assert(zoneVec[znost[ist_nbr]].flag == 0); // this is the stationary side
	sp_flag[isp_prev] = -2;
	const double wgt = DIST(xsp[isp_prev],xsp[isp_next]);
	xp_sum += wgt*(xp_sp[isp_prev][0]+xp_sp[isp_next][0]);
	yp_sum += wgt*(xp_sp[isp_prev][1]+xp_sp[isp_next][1]);
	xp_min = min(xp_min,xp_sp[isp_prev][0]);
	xp_max = max(xp_max,xp_sp[isp_prev][0]);
	w_sum += wgt;
	// and the center and normal...
	FOR_I3 xc_sum[i] += wgt*(xsp[isp_prev][i]+xsp[isp_next][i]);
	double n2[3] = TRI_NORMAL_2(xsp[isp_start],xsp[isp_prev],xsp[isp_next]);
	FOR_I3 normal_sum[i] += n2[i];
	// also the vector out of the shared point...
	const int isp_m = spost[ist][(i+2)%3];
	const double dxp_m = xp_sp[isp_m][0] - 0.5*(xp_sp[isp_prev][0]+xp_sp[isp_next][0]);
	const double dyp_m = xp_sp[isp_m][1] - 0.5*(xp_sp[isp_prev][1]+xp_sp[isp_next][1]);
	const double dxp_m_mag = sqrt(dxp_m*dxp_m+dyp_m*dyp_m);
	dxp_m_sum += wgt*dxp_m/dxp_m_mag;
	dyp_m_sum += wgt*dyp_m/dxp_m_mag;
	if (isp_next == isp_start)
	  break;
      }
      yp_sum /= 2.0*w_sum;
      xp_sum /= 2.0*w_sum;
      FOR_I3 xc_sum[i] /= 2.0*w_sum;
      double nmag = MAG(normal_sum);
      FOR_I3 normal_sum[i] /= nmag;
      cout << " > open edge loop with " << count << 
	" edges: avg xp: " << xp_sum << " dxp: " << xp_min-xp_sum << " " << xp_max-xp_sum << " normal: " << COUT_VEC(normal_sum) << " xc: " << COUT_VEC(xc_sum) << endl;
      
      const double dxp_m_sum_mag = sqrt(dxp_m_sum*dxp_m_sum+dyp_m_sum*dyp_m_sum);

      // the orientation of the normal determines which loop this is...
      const double dp = DOT_PRODUCT(normal_sum,unit_axis);
      if (fabs(1.0-dp) < 1.0E-10) {
	cout << "pos dp!" << endl;
	//assert(0);
	b_loop1 = true;
	xp_loop1 = xp_sum;
	yp_loop1 = yp_sum;
	cout << "LOOP1 " << xp_loop1 << " " << yp_loop1 << endl;
	cout << "LOOP1 " << xp_loop1+0.005*dxp_m_sum/dxp_m_sum_mag << " " << yp_loop1+0.005*dyp_m_sum/dxp_m_sum_mag << endl;
	cout << "LOOP1 " << xp_loop1+0.010*dxp_m_sum/dxp_m_sum_mag << " " << yp_loop1+0.010*dyp_m_sum/dxp_m_sum_mag << endl;
      }
      else if (fabs(dp+1.0) < 1.0E-10) {
	cout << "minus dp!" << endl;
	//assert(!b_loop0);
	b_loop0 = true;
	xp_loop0 = xp_sum;
	yp_loop0 = yp_sum;
	cout << "LOOP0 " << xp_loop0 << " " << yp_loop0 << endl;
	cout << "LOOP0 " << xp_loop0+0.005*dxp_m_sum/dxp_m_sum_mag << " " << yp_loop0+0.005*dyp_m_sum/dxp_m_sum_mag << endl;
	cout << "LOOP0 " << xp_loop0+0.010*dxp_m_sum/dxp_m_sum_mag << " " << yp_loop0+0.010*dyp_m_sum/dxp_m_sum_mag << endl;
      }
      else {
	cout << "misaligned. check normal" << endl;
	assert(0);
      }
      
    }
  }

  delete[] sp_flag; sp_flag = NULL; 

  // now put the tris (edges actually) into a 2D grid for fast lookup...
  // recall that at this point, the st_flag is set as follows:
  // -1: un-checked tri: part of the stationary domain
  //  0: moving tri that touches the axis. depending on orientation, has its axis intersection in yp_axis0/1
  //  1: moving tri that maps to the projected space without crossing the axis

  // extend the domain by 2x the requested dn in all directions...
  yp_min -= 2.0*dn;
  yp_max += 2.0*dn;
  xp_max += 2.0*dn;

  vector<uint> mTeostVec; // vec of moving teost
  vector<uint> sTeostVec; // vec of stationary teost
  
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] == -1) {
      // this tri is part of the stationary domain. First, check if its coordinate 
      // range eliminates it from consideration...
      const int isp0 = spost[ist][0];
      const int isp1 = spost[ist][1];
      const int isp2 = spost[ist][2];
      // note that a tri that intersects the axis can be very large, so we cannot use an
      // r-check (xp-check) to eliminate these checks. We can, hovever, use a yp check...
      const double yp_tri_min = min(xp_sp[isp0][1],min(xp_sp[isp1][1],xp_sp[isp2][1]));
      const double yp_tri_max = max(xp_sp[isp0][1],max(xp_sp[isp1][1],xp_sp[isp2][1]));
      if ((yp_tri_min < yp_max)&&(yp_tri_max > yp_min)) {
	// this is a stationary tri that passes the yp check: check if it intersects the axis...
	double A[3];
	FOR_I3 {
	  // start by expanding the 2D bounding box based on this tri's nodes...
	  const int isp = spost[ist][i];
	  const double dx[3] = DIFF(xsp[isp],xc);
	  // next, we are going to check if this tri intersects the axis...
	  const int isp_next = spost[ist][(i+1)%3];
	  const double dx_next[3] = DIFF(xsp[isp_next],xc);
	  const double cp[3] = CROSS_PRODUCT(dx,dx_next);
	  A[(i+2)%3] = DOT_PRODUCT(cp,unit_axis); // populate the barycentric weight of the opposite node
	}
	if ((A[0] >= 0.0)&&(A[1] >= 0.0)&&(A[2] >= 0.0)&&(A[0]+A[1]+A[2] > 0.0)) {
	  cout << "XXXXXX got positive projected area for stationary tri" << endl;
	  // turn off this tri for later consideration...
	  st_flag[ist] = 0;
	  // use the A's as baricentric weights to store the yp intersection...
	  double yp_avg = 0.0;
	  FOR_I3 {
	    // start by expanding the 2D bounding box based on this tri's nodes...
	    const int isp = spost[ist][i];
	    yp_avg += A[i]*xp_sp[isp][1];
	  }
	  yp_avg /= A[0]+A[1]+A[2];
	  cout << "CCCCCC " << 0.0 << " " << yp_avg << endl;
	  //MPI_Pause("take a look");
	}
	else if ((A[0] <= 0.0)&&(A[1] <= 0.0)&&(A[2] <= 0.0)&&(A[0]+A[1]+A[2] < 0.0)) {
	  cout << "XXXXXX got negative projected area for stationary tri" << endl;
	  assert(0);
	}
	else {
	  FOR_I3 {
	    // consider 1 edge at a time...
	    const int isp = spost[ist][i];
	    const int isp_next = spost[ist][(i+1)%3];
	    // consider only 1 orientation of each edge...
	    if ((xp_sp[isp_next][0] > xp_sp[isp][0])||
		((xp_sp[isp_next][0] == xp_sp[isp][0])&&
		 (xp_sp[isp_next][1] > xp_sp[isp][1]))) {
	      // the orientation is good: check that its r-min is in the region of interest...
	      // note that because of the orientation check above, r_min is isp...
	      if (xp_sp[isp][0] < xp_max) {
		sTeostVec.push_back((uint(ist)<<2)|i);
	      }
	    }
	  }
	}
      }
    }
    else if (st_flag[ist] == 1) {
      // this is a moving tri that needs to be considered...
      FOR_I3 {
	// consider 1 edge at a time...
	const int isp = spost[ist][i];
	const int isp_next = spost[ist][(i+1)%3];
	// consider only 1 orientation of each edge...
	if ((xp_sp[isp_next][0] > xp_sp[isp][0])||
	    ((xp_sp[isp_next][0] == xp_sp[isp][0])&&
	     (xp_sp[isp_next][1] > xp_sp[isp][1]))) {
	  // the orientation is good: check that its r-min is in the region of interest...
	  // note that because of the orientation check above, r_min is isp...
	  if (xp_sp[isp][0] < xp_max) {
	    mTeostVec.push_back((uint(ist)<<2)|i);
	  }
	}
      }
    }
  }
  
  delete[] st_flag; st_flag = NULL;
  
  const int sTeostVec_size = sTeostVec.size();
  const int mTeostVec_size = mTeostVec.size();

  // write out the moving and stationary projected points...
  {
    // stationary...
    FILE * fp = fopen("s.dat","w");
    for (int ii = 0; ii < sTeostVec_size; ++ii) {
      const int ist = sTeostVec[ii]>>2;
      const int i = sTeostVec[ii]&3;
      const int isp = spost[ist][i];
      const int isp_next = spost[ist][(i+1)%3];
      fprintf(fp,"%18.15e %18.15e\n",xp_sp[isp][0],xp_sp[isp][1]);
      fprintf(fp,"%18.15e %18.15e\n\n",xp_sp[isp_next][0],xp_sp[isp_next][1]);
    }
    fclose(fp);
    // moving...
    fp = fopen("m.dat","w");
    for (int ii = 0; ii < mTeostVec_size; ++ii) {
      const int ist = mTeostVec[ii]>>2;
      const int i = mTeostVec[ii]&3;
      const int isp = spost[ist][i];
      const int isp_next = spost[ist][(i+1)%3];
      fprintf(fp,"%18.15e %18.15e\n",xp_sp[isp][0],xp_sp[isp][1]);
      fprintf(fp,"%18.15e %18.15e\n\n",xp_sp[isp_next][0],xp_sp[isp_next][1]);
    }
    fclose(fp);
  }

  double (*bbmin)[2] = new double[max(sTeostVec_size,mTeostVec_size)][2];
  double (*bbmax)[2] = new double[max(sTeostVec_size,mTeostVec_size)][2];

  for (int ii = 0; ii < sTeostVec_size; ++ii) {
    const int ist = sTeostVec[ii]>>2;
    const int i = sTeostVec[ii]&3;
    const int isp = spost[ist][i];
    const int isp_next = spost[ist][(i+1)%3];
    FOR_J2 bbmin[ii][j] = min(xp_sp[isp][j],xp_sp[isp_next][j]);
    FOR_J2 bbmax[ii][j] = max(xp_sp[isp][j],xp_sp[isp_next][j]);
  }
  
  Adt2d<double> sAdt(sTeostVec_size,bbmin,bbmax);

  for (int ii = 0; ii < mTeostVec_size; ++ii) {
    const int ist = mTeostVec[ii]>>2;
    const int i = mTeostVec[ii]&3;
    const int isp = spost[ist][i];
    const int isp_next = spost[ist][(i+1)%3];
    FOR_J2 bbmin[ii][j] = min(xp_sp[isp][j],xp_sp[isp_next][j]);
    FOR_J2 bbmax[ii][j] = max(xp_sp[isp][j],xp_sp[isp_next][j]);
  }
  
  Adt2d<double> mAdt(mTeostVec_size,bbmin,bbmax);

  delete[] bbmin; bbmin = NULL;
  delete[] bbmax; bbmax = NULL;
  
  // now do the walking...
  
  double xp_prev,yp_prev;
  double xp_end,yp_end;
  
  if (b_loop0 && b_axis1 && !b_loop1 && !b_axis0) {
    
    // loop to top axis...
    xp_prev = xp_loop0;
    yp_prev = yp_loop0;
    xp_end = 0.0;
    yp_end = yp_axis1+dn;
    
    // axis to bottom loop...
    //xp_end = xp_loop0;
    //yp_end = yp_loop0;
    //xp_prev = 0.0;
    //yp_prev = yp_axis1+dn;
    
  }
  else {

    assert(0);
    
  }
  
  cout << "CCCCCC " << xp_prev << " " << yp_prev << endl;
  cout << "CCCCCC " << xp_end << " " << yp_end << endl;

  //MPI_Pause("about to start...");
  
  // how to start: should be left to the user as potential input...
  double dxp = ds;
  double dyp = 0.0;
  //int count = 0;

  vector<double> segVec;
  segVec.push_back(xp_prev);
  segVec.push_back(yp_prev);
  cout << "BBBBBB " << xp_prev << " " << yp_prev << endl;

  vector<int> intVec;
  
  while (1) {

    // check distance to end point, and finish if we are close...
    double delta_end = sqrt((xp_end-xp_prev)*(xp_end-xp_prev) + (yp_end-yp_prev)*(yp_end-yp_prev));
    if (delta_end < 1.5*ds) { 
      // we are done: add the end and break out...
      segVec.push_back(xp_end);
      segVec.push_back(yp_end);
      cout << "BBBBBB " << xp_end << " " << yp_end << endl;
      break;
    }

    // now solve using a sampling then bisection approach. The iangle is an index that 
    // is used to access an oscillating angle function... 
    int iangle = 1;
    bool b_below = false;
    bool b_above = false;
    double angle0,angle1;
    double xp_next,yp_next;
    
    while (1) {
      
      double angle;
      if (iangle >= 1) {
	assert(iangle < 40); // fail for large iangle
	angle = M_PI*40.0/180*ANGLEFUNC(iangle);
      }
      else {
	// bisection mode...
	angle = 0.5*(angle0+angle1);
      }
      
      xp_next = xp_prev + dxp*cos(angle) - dyp*sin(angle); 
      yp_next = yp_prev + dyp*cos(angle) + dxp*sin(angle);
      double xpyp_next[2] = { xp_next, yp_next };

      assert(intVec.empty());
      mAdt.buildListForClosestPoint(intVec,xpyp_next);
      
      double d2_m = 1.0E+20;
      double dd2dx_m,dd2dy_m;
      for (int ii = 0, nii = intVec.size(); ii < nii; ++ii) {
	const int ist = mTeostVec[intVec[ii]]>>2;
	const int i = mTeostVec[intVec[ii]]&3;
	const int isp = spost[ist][i];
	const int isp_next = spost[ist][(i+1)%3];
	double d2,dd2dx,dd2dy;
	getPointEdgeDist2(d2,dd2dx,dd2dy,xp_next,yp_next,xp_sp[isp][0],xp_sp[isp][1],xp_sp[isp_next][0],xp_sp[isp_next][1]);
	if (d2 < d2_m) {
	  d2_m = d2;
	  dd2dx_m = dd2dx;
	  dd2dy_m = dd2dy;
	}
      }
      
      intVec.clear();
      sAdt.buildListForClosestPoint(intVec,xpyp_next,1.01*d2_m);
      
      double d2_s = 1.0E+20;
      double dd2dx_s,dd2dy_s;
      for (int ii = 0, nii = intVec.size(); ii < nii; ++ii) {
	const int ist = sTeostVec[intVec[ii]]>>2;
	const int i = sTeostVec[intVec[ii]]&3;
	const int isp = spost[ist][i];
	const int isp_next = spost[ist][(i+1)%3];
	double d2,dd2dx,dd2dy;
	getPointEdgeDist2(d2,dd2dx,dd2dy,xp_next,yp_next,xp_sp[isp][0],xp_sp[isp][1],xp_sp[isp_next][0],xp_sp[isp_next][1]);
	if (d2 < d2_s) {
	  d2_s = d2;
	  dd2dx_s = dd2dx;
	  dd2dy_s = dd2dy;
	}
      }
      
      intVec.clear();

      cout << "angle: " << angle*180.0/M_PI << " d_m: " << sqrt(d2_m) << " d_s: " << sqrt(d2_s) << " dn: " << dn << endl;
	
      if (d2_s < dn*dn) {
	if (fabs(d2_m-d2_s) < 1.0E-12*(d2_m+d2_s))
	  break;
      }
      else {
	const double eps2 = d2_m-dn*dn;
	if (fabs(eps2) < 1.0E-12*dn*dn) {
	  break;
	}
      }
      
      if (d2_m > min(d2_s,dn*dn)) {
	// this is "ABOVE": i.e. too far from moving surface...
	if (iangle >= 1) {
	  b_above = true;
	  if (b_below) {
	    // the sampling has discovered angle limits, so start bisection...
	    angle1 = angle;
	    angle0 = M_PI*40.0/180*ANGLEFUNC(iangle-2); // make sure this is the same as the use above
	    iangle = -1; // trigger bisection mode
	  }
	  else {
	    // advance to the next sampling angle...
	    ++iangle;
	  }
	}
	else {
	  angle1 = angle;
	}
      }
      else {
	// this is "BELOW": i.e. too close to moving surface...
	if (iangle >= 1) {
	  b_below = true;
	  if (b_above) {
	    // the sampling has discovered angle limits, so start bisection...
	    angle0 = angle;
	    angle1 = M_PI*40.0/180*ANGLEFUNC(iangle-2); // see above
	    iangle = -1; // trigger bisection mode
	  }
	  else {
	    // advance to the next sampling angle...
	    ++iangle;
	  }
	}
	else {
	  angle0 = angle;
	}
      }
      
      //MPI_Pause("ok");
      
    }

    // set dxp,dyp...
    segVec.push_back(xp_next);
    segVec.push_back(yp_next);
    cout << "BBBBBB " << xp_next << " " << yp_next << endl;
    
    dxp = xp_next-xp_prev;
    dyp = yp_next-yp_prev;
    xp_prev = xp_next;
    yp_prev = yp_next;
    
    //MPI_Pause("done this point");
    
  }

  delete[] xp_sp;

  //cout << "segVec.size(): " << segVec.size() << endl;

  addRevolve(segVec,xc,axis,128);

  //MPI_Pause("done! take a look");
  //assert(0);

}

// newton method based on gradients...
/*



-       if (fabs(d2_m-d2_s) < 1.0E-6*(d2_m+d2_s))
+      if (d2_s < dn*dn) {
+       if (fabs(d2_m-d2_s) < 1.0E-12*(d2_m+d2_s))
          break;
-         
-       double t1 = d2_m-d2_s;
-       double t2 = dd2dx_m-dd2dx_s;
-       double t3 = dd2dy_m-dd2dy_s;
-       // HACK: adding a factor of 2 here makes convergence faster...
-       xp_next -= t1*t2/(t2*t2+t3*t3);
-       yp_next -= t1*t3/(t2*t2+t3*t3);
-         
       }
       else {
-
-       cout << "d2_s GT: dd2dx_m: " << dd2dx_m << " dd2dy_m: " << dd2dy_m << endl;
-         
-       const double eps2 = d2_m-standoff*standoff;
-       if (fabs(eps2) < 1.0E-6*standoff*standoff) {
+       const double eps2 = d2_m-dn*dn;
+       if (fabs(eps2) < 1.0E-12*dn*dn) {
          break;
        }
-
-       xp_next -= eps2*dd2dx_m/(dd2dx_m*dd2dx_m+dd2dy_m*dd2dy_m);
-       yp_next -= eps2*dd2dy_m/(dd2dx_m*dd2dx_m+dd2dy_m*dd2dy_m);
-
-       cout << "AAAAAA " << xp_next << " " << yp_next << endl;
-
       }
*/

int SimpleSurface::addRevolve(vector<double>& segVec,const double xc[3],const double axis[3],const int ntheta) {
  
  // step 1: decide if the start and end of the segment is on the axis. It affects the new 
  // count of tris. This is pretty easy because segVec constains [r0,z0,r1,z1...], so we just need to look
  // for a r == 0...

  // define R,Z accessors into segVec so we can easily change this some day...
#define RSEG(I) segVec[2*(I)]
#define ZSEG(I) segVec[2*(I)+1]
  
  // define the number of segments as 1 less than the number of point pairs in segVec...
  assert(segVec.size()%2 == 0);
  const int ns = segVec.size()/2-1;
  assert(ns >= 2);

  // the start and end may be on the axis: this is easily detected by
  // R == 0...
  bool b_start = false;
  if (RSEG(0) == 0.0) b_start = true;

  bool b_end = false;
  if (RSEG(ns) == 0.0) b_end = true;
  
  cout << "b_start: " << b_start << " b_end: " << b_end << endl;
  
  int nsp_new = 0;
  int nst_new = 0;
  // start...
  if (b_start) {
    nsp_new += 1 + ntheta;
    nst_new += ntheta;
  }
  else {
    nsp_new += 2*ntheta;
    nst_new += 2*ntheta;
  }
  // middle...
  nsp_new += (ns-2)*ntheta;
  nst_new += 2*(ns-2)*ntheta;
  // end...
  if (b_end) {
    nsp_new += 1;
    nst_new += ntheta;
  }
  else {
    nsp_new += ntheta;
    nst_new += 2*ntheta;
  }

  // resize...
  const int nsp0 = nsp;
  nsp += nsp_new;
  growNspData(nsp,nsp0);
  
  const int nst0 = nst;
  nst += nst_new;
  growNstData(nst,nst0);

  // build the basis vectors...
  const double axis_mag = MAG(axis);
  assert(axis_mag > 0.0);
  const double e0[3] = { axis[0]/axis_mag, axis[1]/axis_mag, axis[2]/axis_mag };
  double e1[3],e2[3];
  MiscUtils::getBestE1E2FromE0(e1,e2,e0);
  
  // add the geometry...
  int isp = nsp0;
  int ist = nst0;
  if (b_start) {
    // first point on axis...
    assert(RSEG(0) == 0.0);
    FOR_J3 xsp[isp][j] = xc[j] + ZSEG(0)*e0[j];
    for (int i = 0; i < ntheta; ++i) {
      const double theta = 2.0*M_PI*double(i)/double(ntheta); 
      const double cos_theta = cos(theta);
      const double sin_theta = sin(theta);
      FOR_J3 xsp[isp+1+i][j] = xc[j] + ZSEG(1)*e0[j] + RSEG(1)*(cos_theta*e1[j] + sin_theta*e2[j]);
      spost[ist+i][0] = isp;
      spost[ist+i][1] = isp+1+i;
      spost[ist+i][2] = isp+1+(i+1)%ntheta;
    }
    isp += 1 + ntheta;
    ist += ntheta;
  }
  else {
    for (int i = 0; i < ntheta; ++i) {
      const double theta = 2.0*M_PI*double(i)/double(ntheta); 
      const double cos_theta = cos(theta);
      const double sin_theta = sin(theta);
      FOR_J3 xsp[isp+i][j]        = xc[j] + ZSEG(0)*e0[j] + RSEG(0)*(cos_theta*e1[j] + sin_theta*e2[j]);
      FOR_J3 xsp[isp+ntheta+i][j] = xc[j] + ZSEG(1)*e0[j] + RSEG(1)*(cos_theta*e1[j] + sin_theta*e2[j]);
      spost[ist+2*i][0]   = isp+i;
      spost[ist+2*i][1]   = isp+ntheta+i;
      spost[ist+2*i][2]   = isp+ntheta+(i+1)%ntheta;
      spost[ist+2*i+1][0] = isp+(i+1)%ntheta;
      spost[ist+2*i+1][1] = isp+i;
      spost[ist+2*i+1][2] = isp+ntheta+(i+1)%ntheta;
    }
    isp += 2*ntheta;
    ist += 2*ntheta;
  }
  // middle...
  for (int is = 1; is < ns-1; ++is) {
    for (int i = 0; i < ntheta; ++i) {
      const double theta = 2.0*M_PI*double(i)/double(ntheta); 
      const double cos_theta = cos(theta);
      const double sin_theta = sin(theta);
      FOR_J3 xsp[isp+i][j] = xc[j] + ZSEG(is+1)*e0[j] + RSEG(is+1)*(cos_theta*e1[j] + sin_theta*e2[j]);
      spost[ist+2*i][0]   = isp-ntheta+i;
      spost[ist+2*i][1]   = isp+i;
      spost[ist+2*i][2]   = isp+(i+1)%ntheta;
      spost[ist+2*i+1][0] = isp-ntheta+(i+1)%ntheta;
      spost[ist+2*i+1][1] = isp-ntheta+i;
      spost[ist+2*i+1][2] = isp+(i+1)%ntheta;
    }
    isp += ntheta;
    ist += 2*ntheta;
  }
  // end...
  if (b_end) {
    FOR_J3 xsp[isp][j] = xc[j] + ZSEG(ns)*e0[j];
    for (int i = 0; i < ntheta; ++i) {
      const double theta = 2.0*M_PI*double(i)/double(ntheta); 
      const double cos_theta = cos(theta);
      const double sin_theta = sin(theta);
      spost[ist+i][0] = isp-ntheta+(i+1)%ntheta;
      spost[ist+i][1] = isp-ntheta+i;
      spost[ist+i][2] = isp;
    }
    isp += 1;
    ist += ntheta;
  }
  else {
    assert(0);
    //nsp_new += n;
    //nst_new += 2*n;
  }
  assert(isp == nsp);
  assert(ist == nst);

  // finally the zone...
  const string zone_name = "FARFIELD";
  int izone;
  int nzn0 = zoneVec.size();
  for (izone = 0; izone < nzn0; ++izone) {
    if (zoneVec[izone].getName() == zone_name)
      break;
  }
  if (izone == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name));
    ++nsz;
  }
  for (int ist = nst0; ist < nst; ++ist) {
    znost[ist] = izone;
  }
  
  //MPI_Pause("in add revolve -- LOOKS GOOD");

  clearTeost();
  
#undef RSEG
#undef ZSEG
  
}
  
  
  
#ifdef fsdhjfhsdksfd
  
  
  int isp = ss_nsp0;
  if (capSurf || specifiedCap0) {
    FOR_I3 xsp[isp][i] = x0[i];
    isp++;  // cap0 endpoint
  }

  // place all irVec points
  for (int ixy=0, limit=irVec.size(); ixy < limit; ++ixy) {
    double center[3];
    FOR_I3 {
      center[i] = x0[i] + irVec[ixy].first*axis[i];
    }
    MiscUtils::createCirclePts(xsp,isp,center,axis,irVec[ixy].second,n_theta,false); // do NOT stagger the points
    isp += n_theta;
  }

  if (capSurf || specifiedCap1) {
    FOR_I3 xsp[isp][i] = x1[i];
    isp++;
  }

  assert(isp == nsp);

  // Add new tris
  const int ss_nst0 = nst;

  nst += 2*n_theta*(irVec.size()-1);
  if (capSurf || specifiedCap0) nst += n_theta;  // accounts for starting cap facets
  if (capSurf || specifiedCap1) nst += n_theta;  // accounts for ending cap facets

  
  


  
  



const double rAxis[3],const double x0[3],double (*ir_vals)[2],const int n_vals,const int _n_theta, bool capSurf, const double maxAR) {
  COUT1("SimpleSurface::initRevolve()");

  if (n_vals < 2) {
    CWARN("not enough points were specified in the profile; skipping");
    return -1;
  }
  if (_n_theta < 3) {
    CWARN("N_THETA of \"" << _n_theta << "\" is too small; defaulting to 3");
  }
  const int n_theta = max(3,_n_theta);

  double axis[3];
  FOR_I3 {
    axis[i] = rAxis[i];
  }
  NORMALIZE(axis); // make a unit vector
  COUT2(" > revolve axis = " << COUT_VEC( axis ) );

  vector< pair<double,double> > irVec;
  for (int v=0; v < n_vals; ++v) {
    irVec.push_back(pair<double,double> (ir_vals[v][0],ir_vals[v][1]));
  }

  for (int ixy = 1, limit=irVec.size()-1; ixy < limit; ++ixy) {
    if (irVec[ixy].second <= 0.0) {
      CWARN("cannot support profile points on the axis itself unless the starting or ending point; skipping");
      return -1;
    }
  }

  bool specifiedCap0 = false;
  bool specifiedCap1 = false;
  if (irVec[0].second == 0.0) {
    specifiedCap0 = true;
  } else {
    assert(irVec[0].second > 0.0);
  }

  if (irVec[irVec.size()-1].second == 0.0) {
    specifiedCap1 = true;
  } else {
    assert(irVec[irVec.size()-1].second > 0.0);
  }

  // if specified, add points within irVec to maintain certain aspect ratio
  // of facets. Since we use a fixed "n", we only look at adjustments introduced
  // by adding points in (I,R) space.
  if (maxAR != -1.0) {
    assert(maxAR >= 1.0);
    COUT2(" > aspect ratio limit: " << maxAR);

    if (capSurf) {
      // add point(s) to close the surface
      if (!specifiedCap0) {
        irVec.insert( irVec.begin(), pair<double,double>(irVec.front().first,0.0) );
        COUT2(" > adding starting point to surface at: (" << irVec.front().first << "," << irVec.front().second << ")");
      }

      if (!specifiedCap1) {
        irVec.push_back( pair<double,double>(irVec.back().first,0.0) );
        COUT2(" > adding closing point to surface at: (" << irVec.back().first << "," << irVec.back().second << ")");
      }
    }

    vector< pair<double,double> > irVec_new;
    irVec_new.push_back(irVec[0]);  //start from first point in IR profile

    int pos = 0;
    for (int i=1,limit=irVec.size(); i < limit; i++) {
      // COUT2("segment[" << i << "]");
      const double dI = irVec[i].first - irVec_new[pos].first;
      const double dR = irVec[i].second - irVec_new[pos].second;
      const double L0 = sqrt(dI*dI + dR*dR);
      double L1;
      if (i != limit-1) {
        L1 = 2*irVec[i].second*sin(M_PI/double(n_theta));
      }
      else {
        L1 = 2*irVec_new[pos].second*sin(M_PI/double(n_theta));
      }
      const double L2 = sqrt(L0*L0 + L1*L1);

      // compute aspect ratio based on ratio of circumradius/(2*inradius)
      // if ratio below criteria, add point(s) to profile
      const double s = 0.5*(L0+L1+L2);
      const double AR = L0*L1*L2 / (8.0*(s-L0)*(s-L1)*(s-L2));

      // COUT2(" > lengths L0,L1,L2,s: " << L0 << " " << L1 << " " << L2 << " " << s);
      // COUT2(" > computed AR: " << AR);

      if ((AR > maxAR) && (max(L0,L2) > L1)) {
        const int newPoints = (int)floor(AR/maxAR);
        // COUT2("    > inserting " << newPoints << " points");
        const double deltaI = dI / double(newPoints+1);
        const double deltaR = dR / double(newPoints+1);
        for (int p=0; p < newPoints; p++) {
          const double newI = irVec_new[pos].first + (p+1)*deltaI;
          const double newR = irVec_new[pos].second + (p+1)*deltaR;
          irVec_new.push_back(pair<double,double>(newI,newR));
        }
        pos += newPoints;
      }
      // add specified point
      irVec_new.push_back(irVec[i]);
      ++pos;
    }
    irVec.swap(irVec_new);

    if (capSurf) {
      // remove endpoints from irVec - used only for AR limiting
      if (!specifiedCap0) irVec.erase(irVec.begin());
      if (!specifiedCap1) irVec.pop_back();
    }
  }

  double x1[3];
  if (capSurf) {
    FOR_I3 x1[i] = x0[i] + irVec.back().first*axis[i];
  }

  if (specifiedCap0) {
    // remove endpoints from irVecs
    irVec.erase(irVec.begin());
  }

  if (specifiedCap1) {
    FOR_I3 x1[i] = x0[i] + irVec.back().first*axis[i];
    irVec.pop_back();
  }
  COUT2("x0:" << COUT_VEC(x0));
  COUT2("x1:" << COUT_VEC(x1));

  for (int i=0,limit=irVec.size(); i < limit; i++) {
    COUT2("[" << i << "]: " << irVec[i].first << " " << irVec[i].second);
  }
  COUT2(" > number of IR points prescribed: " << irVec.size());


  // Add new points
  const int ss_nsp0 = nsp;

  nsp += irVec.size()*n_theta;
  if (capSurf || specifiedCap0) nsp += 1;  // account for starting point
  if (capSurf || specifiedCap1) nsp += 1;  // account for ending point

  growNspData(nsp,ss_nsp0);

  int isp = ss_nsp0;
  if (capSurf || specifiedCap0) {
    FOR_I3 xsp[isp][i] = x0[i];
    isp++;  // cap0 endpoint
  }

  // place all irVec points
  for (int ixy=0, limit=irVec.size(); ixy < limit; ++ixy) {
    double center[3];
    FOR_I3 {
      center[i] = x0[i] + irVec[ixy].first*axis[i];
    }
    MiscUtils::createCirclePts(xsp,isp,center,axis,irVec[ixy].second,n_theta,false); // do NOT stagger the points
    isp += n_theta;
  }

  if (capSurf || specifiedCap1) {
    FOR_I3 xsp[isp][i] = x1[i];
    isp++;
  }

  assert(isp == nsp);

  // Add new tris
  const int ss_nst0 = nst;

  nst += 2*n_theta*(irVec.size()-1);
  if (capSurf || specifiedCap0) nst += n_theta;  // accounts for starting cap facets
  if (capSurf || specifiedCap1) nst += n_theta;  // accounts for ending cap facets

  growNstData(nst,ss_nst0);
  const int ss_nzn0 = zoneVec.size();

  COUT1(" > zones:");
  zoneVec.push_back(SurfaceZone("wall")); COUT1("    > \"wall\"");
  int cap0_idx=0;  // default write to wall zone
  int cap1_idx=0;  // default write to wall zone
  if (capSurf || specifiedCap0) {
    zoneVec.push_back(SurfaceZone("cap0")); COUT1("    > \"cap0\"");
    cap0_idx++;
    cap1_idx++;
  }
  if (capSurf || specifiedCap1) {
    zoneVec.push_back(SurfaceZone("cap1")); COUT1("    > \"cap1\"");
    cap1_idx++;
  }

  // place all surface tris
  int ist = ss_nst0;

  if (capSurf || specifiedCap0) {
    // cap0
    facetCircleToPoint(spost,znost,ss_nsp0,ss_nsp0+1,ist,cap0_idx,n_theta,false);
    ist += n_theta;
  }

  for (int ixy=0,limit=irVec.size()-1; ixy < limit; ++ixy) {
    // revolve
    int c0  = ss_nsp0 + ixy*n_theta;
    int c1  = ss_nsp0 + ixy*n_theta + n_theta;
    if (capSurf || specifiedCap0) {
      c0 += 1;
      c1 += 1;
    }
    facetCircleToCircle(spost,znost,c0,c1,ist,0,n_theta,true,true);
    ist += 2*n_theta;
  }

  if (capSurf || specifiedCap1) {
    // cap1
    facetCircleToPoint(spost,znost,nsp-1,nsp-n_theta-1,ist,cap1_idx,n_theta,true);
    ist += n_theta;
  }

  assert(ist == nst);

  if (ss_nst0 > 0) {
    // if adding surface, properly offset zone indices
    for (int ist = ss_nst0; ist < nst; ++ist) {
      // FOR_I3 spost[ist][i] += ss_nsp0;
      znost[ist] += ss_nzn0;
    }
  }

  COUT1(" > addRevolve done");
  // ====================================================
  // make sure all zones were set...
  for (int ist = 0; ist < nst; ++ist) {
    assert((znost[ist] >= 0)&&(znost[ist] < int(zoneVec.size())));
  }

  return 0;  // success
}

#endif
