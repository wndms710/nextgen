#include "SimpleSurface.hpp"
#include "WebUI.hpp"

void SimpleSurface::setFeatureAngle(const double degrees) {
  if (degrees == getFeatureAngle()) {
    COUT1(" > feature angle already set to: " << degrees << " degrees");
    return;
  }
  COUT1(" > feature angle set to: " << degrees << " degrees");
  feature_cos = cos((180.0-degrees)*M_PI/180.0);

  if (eoi_type == FEATURE) clearDynamicEdgeGroups();  // if new feature angle we need to rebuild edge groups accordingly
  // COUT2("verified: " << getFeatureAngle());
}

double SimpleSurface::getFeatureAngle() const {
  return -1.0*((acos(feature_cos)*180.0/M_PI)-180.0);
}

void SimpleSurface::calcGcl(double gcl[3], const bool include_open_edge_groups) {

  // surface unit fluxes
  FOR_I3 gcl[i] = 0.0;
  for (int ist = 0; ist < nst; ++ist) {
    const double n_st[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    FOR_I3 gcl[i] += n_st[i];
  }

  if (include_open_edge_groups) {

    // need teost with open edge group indexing
    ensureOpenEdgeGroups();

    // add hole unit fluxes
    int *isp0_igr = new int[n_open_edge_groups];
    for (int igr = 0; igr < n_open_edge_groups; ++igr) isp0_igr[igr] = -1;
    FOR_IST {
      FOR_I3 {
        int igr;
        if (getOpenEdgeGroup(igr,ist,i)) {
          const int isp0 = spost[ist][i];
          if (isp0_igr[igr] < 0) {
            isp0_igr[igr] = isp0; // store first point
          }
          else {
            const int isp1 = spost[ist][(i+1)%3];
            const double n_st[3] = TRI_NORMAL_2(xsp[isp0_igr[igr]],xsp[isp1],xsp[isp0]);
            FOR_I3 gcl[i] += n_st[i];
          }
        }
      }
    }
    delete[] isp0_igr;

    cout << "GCL (+ open edge group unit fluxes): " << COUT_VEC(gcl) << endl;

  }
  else {

    cout << "GCL: " << COUT_VEC(gcl) << endl;

  }

}

void SimpleSurface::ensureSposp() {
  if (!b_sposp) {
    buildSposp();
    assert(b_sposp);
  }
}

void SimpleSurface::buildSposp() {

  assert(!b_sposp);
  ensureTeost();

  assert(sposp_i == NULL);
  sposp_i = new int[nsp+1];
  FOR_ISP sposp_i[isp+1] = 1; // include yourself

  // first pass to count...
  FOR_IST {
    FOR_I3 {
      int ist_nbr, i_nbr, orient_nbr;
      const bool has_single_edge_nbr = getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i);
      if (has_single_edge_nbr) {
        if (ist_nbr > ist) {
          // only increment for lower tri index
          ++sposp_i[spost[ist][i]+1];
          ++sposp_i[spost[ist][(i+1)%3]+1];
        }
      }
      else if (isEdgeOpen(ist,i)) {
        // this is the only link between these nodes, so increment both...
        ++sposp_i[spost[ist][i]+1];
        ++sposp_i[spost[ist][(i+1)%3]+1];
      }
    }
  }

  // convert to range...
  sposp_i[0] = 0;
  FOR_ISP sposp_i[isp+1] += sposp_i[isp];
  const int sposp_s = sposp_i[nsp];

  assert(sposp_v == NULL);
  sposp_v = new int[sposp_s];

  // populate...
  FOR_ISP sposp_v[sposp_i[isp]++] = isp; // include yourself
  FOR_IST {
    FOR_I3 {
      int ist_nbr, i_nbr, orient_nbr;
      const bool has_single_edge_nbr = getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i);
      if (has_single_edge_nbr) {
        if (ist_nbr > ist) {
          // only increment for lower tri index
          sposp_v[sposp_i[spost[ist][i]]++] = spost[ist][(i+1)%3];
          sposp_v[sposp_i[spost[ist][(i+1)%3]]++] = spost[ist][i];
        }
      }
      else if (isEdgeOpen(ist,i)) {
        // this is the only link between these nodes, so increment both...
        sposp_v[sposp_i[spost[ist][i]]++] = spost[ist][(i+1)%3];
        sposp_v[sposp_i[spost[ist][(i+1)%3]]++] = spost[ist][i];
      }
    }
  }

  // and rewind...
  for (int isp = nsp-1; isp > 0; --isp) sposp_i[isp] = sposp_i[isp-1];
  sposp_i[0] = 0;

  b_sposp = true;

}

void SimpleSurface::clearSposp() {

  b_sposp = false;
  DELETE(sposp_i);
  DELETE(sposp_v);

}

// this routine determines if a point is in a triangle by using barycentric coordinates
bool SimpleSurface::pointInTriangle(const double a[3], const double b[3], const double c[3], const double p[3], const double eps) {

  // determine if p projected onto the plane parallel to tri(a,b,c) is in tri(a,b,c)
  const double u[3] = DIFF(c,a);
  const double v[3] = DIFF(b,a);
  const double w[3] = DIFF(p,a);
  const double n2A[3] = CROSS_PRODUCT(u,v);
  const double mag_n2A_sq = DOT_PRODUCT(n2A,n2A);
  const double gamma = CROSS_DOT(u,w,n2A);
  const double beta = CROSS_DOT(w,v,n2A);
  const double alpha = mag_n2A_sq-gamma-beta;
  bool point_in_tri;
  // eps is relative to squared area
  if (gamma < -eps*mag_n2A_sq || gamma > (1.0+eps)*mag_n2A_sq ||
      beta < -eps*mag_n2A_sq || beta > (1.0+eps)*mag_n2A_sq ||
      alpha < -eps*mag_n2A_sq || alpha > (1.0+eps)*mag_n2A_sq ||
      // check if p is in plane
      fabs(SIGNED_TET_VOLUME_6(a,b,c,p))/pow((1.0+eps)*mag_n2A_sq,0.75) > 1.0E-6) {
    //cout << alpha << " " << beta << " " << gamma << " " << fabs(SIGNED_TET_VOLUME_6(a,b,c,p))/pow(mag_n2A_sq,0.75) << endl;
    point_in_tri = false;
  }
  else {
    //cout << "flipped tri near" << COUT_VEC(p) << endl;
    //cout << alpha << " " << beta << " " << gamma << " " << fabs(SIGNED_TET_VOLUME_6(a,b,c,p))/pow(mag_n2A_sq,0.75) << endl;
    point_in_tri = true;
  }
  return point_in_tri;
}

// this routine determines if a point is in a triangle by using barycentric coordinates
int SimpleSurface::howIsPointOnTriangle(const double a[3], const double b[3], const double c[3], const double p[3], const double dist_tol) {
  // categorize as either:
  // 0: not on triangle
  // 1: on triangle interior
  // -1,-2,-3: at a node
  // 2,3,4: at an edge

  //determine how close point is to edges
  double d2[3];
  d2[0] = MiscUtils::getPointToEdgeDist2(p,a,b);
  d2[1] = MiscUtils::getPointToEdgeDist2(p,b,c);
  d2[2] = MiscUtils::getPointToEdgeDist2(p,c,a);

  const double d2_tol = dist_tol*dist_tol;

  // categorize if close to any edge
  if (d2[2] <= d2_tol && d2[0] <= d2_tol) return -1;
  if (d2[0] <= d2_tol && d2[1] <= d2_tol) return -2;
  if (d2[2] <= d2_tol && d2[1] <= d2_tol) return -3;
  if (d2[0] <= d2_tol) return 2;
  if (d2[1] <= d2_tol) return 3;
  if (d2[2] <= d2_tol) return 4;

  // otherwise it is in or out of tri
  const bool inTri = pointInTriangle(a,b,c,p);

  return ((inTri) ? 1:0);
}

// these routines are used to triangulate a 2d polygon using ear clipping

// compute the signed area of triangle(a,b,c) (ccw ordering)
// maybe replace this with a macro
double SimpleSurface::triangleSignedArea(const double xa, const double ya, const double xb, const double yb,
    const double xc, const double yc) {

  return 0.5*( (xb - xa)*(yc - ya) - (xc - xa)*(yb - ya) );

}

// Is c between a and b?
bool SimpleSurface::between(const double xa, const double ya, const double xb, const double yb,const double xc, const double yc) {

  bool value;
  if ( !collinear(xa, ya, xb, yb, xc, yc) ) {
    value = false;
  }
  else if ( fabs(ya - yb) < fabs(xa - xb) ) {
    double xmax = max(xa, xb);
    double xmin = min(xa, xb);
    value = (xmin <= xc && xc <= xmax);
  }
  else {
    double ymax = max(ya, yb);
    double ymin = min(ya, yb);
    value = (ymin <= yc && yc <= ymax);
  }
  return value;
}

// are a, b and c collinear to numerical precision?
bool SimpleSurface::collinear(const double xa, const double ya, const double xb, const double yb,const double xc, const double yc) {

  double area = triangleSignedArea(xa, ya, xb, yb, xc, yc);

  double side_ab_sq = (xa - xb)*(xa - xb) + (ya - yb)*(ya - yb);
  double side_bc_sq = (xb - xc)*(xb - xc) + (yb - yc)*(yb - yc);
  double side_ca_sq = (xc - xa)*(xc - xa) + (yc - ya)*(yc - ya);

  double side_max_sq = max(side_ab_sq, max(side_bc_sq, side_ca_sq));

  bool value;
  const double EPS = 2.220446049250313E-16;
  if (side_max_sq <= EPS) {
    value = true;
  }
  else if (2.0*fabs(area) <= EPS*side_max_sq) {
    value = true;
  }
  else {
    value = false;
  }

  return value;
}

// is im1 -> ip1 a proper internal diagonal?
bool SimpleSurface::diagonal(const int im1, const int ip1, const int n, const int* prev_node, const int* next_node,const double* x, const double* y) {

  bool value1 = inCone(im1, ip1, n, prev_node, next_node, x, y);
  bool value2 = inCone(ip1, im1, n, prev_node, next_node, x, y);
  bool value3 = diagonalie(im1, ip1, n, next_node, x, y);

  return ( value1 && value2 && value3 );
}

// is im1->ip1 a proper diagonal?
bool SimpleSurface::diagonalie(const int im1, const int ip1, const int n, const int* next_node, const double* x, const double* y) {
  UNUSED(n);
  int first = im1;
  int j = first;
  int jp1 = next_node[first];
  bool value = true;

  // for all edges
  while ( 1 ) {

    // skip edges that include im1 or ip1
    if ( j == im1 || j == ip1 || jp1 == im1 || jp1 == ip1 ) {
      ;
    }
    else {
      bool value2 = intersect ( x[im1], y[im1], x[ip1], y[ip1], x[j], y[j], x[jp1], y[jp1] );

      if ( value2 ) {
        value = false;
        break;
      }
    }
    j = jp1;
    jp1 = next_node[j];

    if ( j == first ) {
      break;
    }
  }

  return value;
}

// is im1->ip1 an internal diagonal?
bool SimpleSurface::inCone (const int im1, const int ip1, const int n, const int* prev_node, const int* next_node,const double* x, const double* y) {
  UNUSED(n);
  int im2 = prev_node[im1];
  int i = next_node[im1];
  double t1 = triangleSignedArea(x[im1], y[im1], x[i], y[i], x[im2], y[im2] );
  bool value;
  if ( 0.0 <= t1 ) {
    double t2 = triangleSignedArea(x[im1], y[im1], x[ip1], y[ip1], x[im2], y[im2] );
    double t3 = triangleSignedArea(x[ip1], y[ip1], x[im1], y[im1], x[i], y[i] );
    value = ( (0.0 < t2) && (0.0 < t3) );
  }
  else {
    double t4 = triangleSignedArea(x[im1], y[im1], x[ip1], y[ip1], x[i], y[i] );
    double t5 = triangleSignedArea(x[ip1], y[ip1], x[im1], y[im1], x[im2], y[im2] );
    value = !( (0.0 <= t4) && (0.0 <= t5) );
  }
  return value;
}

// does a,b intersect with b,c?
bool SimpleSurface::intersect (const double xa, const double ya, const double xb, const double yb,const double xc, const double yc, const double xd, const double yd) {

  bool value;
  if ( intersectProp(xa, ya, xb, yb, xc, yc, xd, yd) ) {
    value = true;
  }
  else if ( between(xa, ya, xb, yb, xc, yc) ) {
    value = true;
  }
  else if ( between(xa, ya, xb, yb, xd, yd) ) {
    value = true;
  }
  else if ( between(xc, yc, xd, yd, xa, ya) ) {
    value = true;
  }
  else if ( between(xc, yc, xd, yd, xb, yb) ) {
    value = true;
  }
  else {
    value = false;
  }
  return value;
}

// do a,b and c,d have a proper intersection?
bool SimpleSurface::intersectProp (const double xa, const double ya, const double xb, const double yb,const double xc, const double yc, const double xd, const double yd) {

  bool value;
  if ( collinear(xa, ya, xb, yb, xc, yc) ) {
    value = false;
  }
  else if ( collinear(xa, ya, xb, yb, xd, yd) ) {
    value = false;
  }
  else if ( collinear(xc, yc, xd, yd, xa, ya) ) {
    value = false;
  }
  else if ( collinear(xc, yc, xd, yd, xb, yb) ) {
    value = false;
  }
  else {
    double t1 = triangleSignedArea ( xa, ya, xb, yb, xc, yc );
    double t2 = triangleSignedArea ( xa, ya, xb, yb, xd, yd );
    double t3 = triangleSignedArea ( xc, yc, xd, yd, xa, ya );
    double t4 = triangleSignedArea ( xc, yc, xd, yd, xb, yb );

    bool value1 = ( 0.0 < t1 );
    bool value2 = ( 0.0 < t2 );
    bool value3 = ( 0.0 < t3 );
    bool value4 = ( 0.0 < t4 );

    value = ( value1 != value2 ) && ( value3 != value4 );
  }
  return value;
}

// the signed area of a polygon
double SimpleSurface::polygonSignedArea (const int n, const double* x, const double* y) {

  // find mean
  double xm = 0.0;
  double ym = 0.0;
  for (int i = 0; i < n; i++) {
    xm += x[i];
    ym += y[i];
  }
  xm /= n;
  ym /= n;

  // calc area
  double area = 0.0;
  /*
  for (int i = 0, im1 = n-1; i < n; i++) {
    //area = area + x[im1] * y[i] - x[i] * y[im1];
    im1 = i;
  }
  area = 0.5 * area;
  */
  for (int i = 0, im1 = n-1; i < n; i++) {
    area = area + 2.0*triangleSignedArea(xm,ym,x[im1],y[im1],x[i],y[i]); // more accurate, little slower
    im1 = i;
  }

  return area;
}

// triangulate polygon by ear clipping
// return false if some assertion failed
bool SimpleSurface::polygonTriangulate(int* triangles, const int n, const double* x, const double* y) {

  // polygon must atleast have 3 points
  //assert( n >= 3 );
  if (n < 3) return false;

  // neighboring points can't be collocated
  for (int node = 0, node_m1 = n-1; node < n; node++) {
    //assert( x[node_m1] != x[node] && y[node_m1] != y[node] );
    if ( x[node_m1] == x[node] || y[node_m1] == y[node] ) return false;
    node_m1 = node;
  }

  // min angle that we can support?
  for (int node2 = 0, node1 = n-1; node2 < n; node2++ ) {
    int node3 = (node2 + 1)%n;
    double dx = (x[node3] - x[node2])*(x[node1] - x[node2]) + (y[node3] - y[node2])*(y[node1] - y[node2]);
    double dy = (x[node3] - x[node2])*(y[node1] - y[node2]) - (y[node3] - y[node2])*(x[node1] - x[node2]);
    double angle;
    if ( dx == 0.0 && dy == 0.0 ) {
      angle = 0.0;
    }
    else {
      angle  = atan2(dy, dx);
      if ( angle < 0.0 ) {
        angle = angle+2.0*M_PI;
      }
    }
    // this assertion seems too strict, because the algo worked when it failed
    //assert( fabs(angle) > M_PI/180.0);
    node1 = node2;
  }

  // area must be positive.
  double area = polygonSignedArea(n, x, y);
  //assert( area > 0.0 );
  if ( area <= 0.0 ) return false;

  // previous and next nodes
  int *prev_node = new int[n];
  int *next_node = new int[n];
  prev_node[0] = n-1;
  next_node[0] = 1;
  for (int i = 1; i < n - 1; i++) {
    prev_node[i] = i - 1;
    next_node[i] = i + 1;
  }
  prev_node[n-1] = n-2;
  next_node[n-1] = 0;

  // can the ear between a node and its neighbors be sliced off
  bool *ear = new bool[n];
  for (int i = 0; i < n; i++) {
    ear[i] = diagonal(prev_node[i], next_node[i], n, prev_node, next_node, x, y);
  }

  // continue until we fully triangulated the polygon
  int triangle_num = 0;
  int i2 = 0;
  while ( triangle_num < (n - 3) ) {

    // if i2 is an ear
    if ( ear[i2] ) {
      int i3 = next_node[i2];
      int i4 = next_node[i3];
      int i1 = prev_node[i2];
      int i0 = prev_node[i1];

      // delete i2
      next_node[i1] = i3;
      prev_node[i3] = i1;

      // update the neighboring info
      ear[i1] = diagonal(i0, i3, n, prev_node, next_node, x, y);
      ear[i3] = diagonal(i1, i4, n, prev_node, next_node, x, y);

      // add diagonal to list of tris
      triangles[0+triangle_num*3] = i3;
      triangles[1+triangle_num*3] = i1;
      triangles[2+triangle_num*3] = i2;
      triangle_num++;
    }

    // go to next vertex
    i2 = next_node[i2];
  }

  // add remaining vertices to the list of triangles
  int i3 = next_node[i2];
  int i1 = prev_node[i2];
  triangles[0+triangle_num*3] = i3;
  triangles[1+triangle_num*3] = i1;
  triangles[2+triangle_num*3] = i2;
  triangle_num++;

  // clean
  delete[] ear;
  delete[] next_node;
  delete[] prev_node;

  return true;
}

bool SimpleSurface::pointIsInside(const double xp[3]) {

  // pointIsInside uses an integer adt to count intersections
  // from x = -infinity to the point's x. If the count is ODD, the point is inside.

  if (stAdt2d == NULL) {

    ensureCentroid();
    ensureBoundingBox();

    // consider only (y,z)...
    const double delta_max = max(boundingBox[3]-boundingBox[2],boundingBox[5]-boundingBox[4]);
    const int BIT_MAX = 28; // do not go above 30
    d2i_factor = double(1<<BIT_MAX)/delta_max;

    // now build the 2D integer Adt based on y,z bounding boxes...

    int (*bbmin)[2] = new int[nst][2];
    int (*bbmax)[2] = new int[nst][2];

    for (int ist = 0; ist < nst; ++ist) {
      const int isp0 = spost[ist][0];
      const int isp1 = spost[ist][1];
      const int isp2 = spost[ist][2];
      const double ymin = min(xsp[isp0][1],min(xsp[isp1][1],xsp[isp2][1]));
      const double ymax = max(xsp[isp0][1],max(xsp[isp1][1],xsp[isp2][1]));
      const double zmin = min(xsp[isp0][2],min(xsp[isp1][2],xsp[isp2][2]));
      const double zmax = max(xsp[isp0][2],max(xsp[isp1][2],xsp[isp2][2]));
      // convert to odd integers...
      bbmin[ist][0] = int( (ymin - centroid[1])*d2i_factor ); if (bbmin[ist][0]%2 == 0) bbmin[ist][0] += 1;
      bbmin[ist][1] = int( (zmin - centroid[2])*d2i_factor ); if (bbmin[ist][1]%2 == 0) bbmin[ist][1] += 1;
      bbmax[ist][0] = int( (ymax - centroid[1])*d2i_factor ); if (bbmax[ist][0]%2 == 0) bbmax[ist][0] += 1;
      bbmax[ist][1] = int( (zmax - centroid[2])*d2i_factor ); if (bbmax[ist][1]%2 == 0) bbmax[ist][1] += 1;
    }

    stAdt2d = new Adt2d<int>(nst,bbmin,bbmax);

    delete[] bbmin;
    delete[] bbmax;

  }

  // convert the y,z of the passed point to EVEN integers...

  int xp_int[2];
  xp_int[0] = int( (xp[1] - centroid[1])*d2i_factor ); if (xp_int[0]%2 != 0) xp_int[0] += 1;
  xp_int[1] = int( (xp[2] - centroid[2])*d2i_factor ); if (xp_int[1]%2 != 0) xp_int[1] += 1;

  // grab everything from the adt that matches...

  vector<int> candidateVec;
  stAdt2d->buildListForPoint(candidateVec,xp_int);

  int count = 0;
  for (int ii = 0, limit = candidateVec.size(); ii < limit; ++ii) {
    const int ist = candidateVec[ii];
    // since we are only counting intersections on one side, do not even check tris that
    // cannot support an appropriate intersection...
    if (min(xsp[spost[ist][0]][0],min(xsp[spost[ist][1]][0],xsp[spost[ist][2]][0])) < xp[0]) {
      int8 dx[3][2];
      FOR_I3 {
        const int isp = spost[ist][i];
        int this_xsp_int = int( (xsp[isp][1] - centroid[1])*d2i_factor );
        if (this_xsp_int%2 == 0) {this_xsp_int += 1;}
        int this_ysp_int = int( (xsp[isp][2] - centroid[2])*d2i_factor );
        if (this_ysp_int%2 == 0) {this_ysp_int += 1;}
        dx[i][0] = this_xsp_int - xp_int[0];
        dx[i][1] = this_ysp_int - xp_int[1];
      }
      const int8 A0 = dx[1][0]*dx[2][1] - dx[1][1]*dx[2][0];
      const int8 A1 = dx[2][0]*dx[0][1] - dx[2][1]*dx[0][0];
      const int8 A2 = dx[0][0]*dx[1][1] - dx[0][1]*dx[1][0];
      if ((A0 == 0)&&(A1 == 0)&&(A2 == 0)) continue;
      else if ((A0 >= 0)&&(A1 >= 0)&&(A2 >= 0)) {
        // use the A's to check that the x is less than the xp...
        if ((xsp[spost[ist][0]][0]*double(A0)+xsp[spost[ist][1]][0]*double(A1)+xsp[spost[ist][2]][0]*double(A2))/double(A0+A1+A2) < xp[0]) {
          // this is positive tri or edge intersection...
          if ((A0 == 0)||(A1 == 0)||(A2 == 0)) {
            // this is an edge intersection...
            // the only possibility here is that ONLY ONE A is zero. The staggered space makes
            // vertex intersections impossible...
            assert( (((A0 == 0)&&(A1 > 0))||(A2 > 0)) || (((A0 > 0)&&(A1 == 0))||(A2 > 0)) || (((A0 > 0)&&(A1 > 0))||(A2 == 0)) );
            // this counts as a half positive intersection...
            count += 1;
          }
          else {
            // full positive intersection...
            count += 2;
          }
        }
      }
      else if ((A0 <= 0)&&(A1 <= 0)&&(A2 <= 0)) {
        // use the A's to check that the x is less than the xp...
        if ((xsp[spost[ist][0]][0]*double(A0)+xsp[spost[ist][1]][0]*double(A1)+xsp[spost[ist][2]][0]*double(A2))/double(A0+A1+A2) < xp[0]) {
          // this is negative tri or edge intersection...
          if ((A0 == 0)||(A1 == 0)||(A2 == 0)) {
            // this is an edge intersection...
            // the only possibility here is that ONLY ONE A is zero. The staggered space makes
            // vertex intersections impossible...
            assert( (((A0 == 0)&&(A1 < 0))||(A2 < 0)) || (((A0 < 0)&&(A1 == 0))||(A2 < 0)) || (((A0 < 0)&&(A1 < 0))||(A2 == 0)) );
            // this counts as a half negative intersection...
            count -= 1;
          }
          else {
            // full negative intersection...
            count -= 2;
          }
        }
      }
    }
  }
  return (count != 0);

}
