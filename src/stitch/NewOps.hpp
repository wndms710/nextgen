// ======================================================================
// revisited these stabilized reconstructions in March/April 2019
//
// To run: stitch.exe on one core with input stitch.newd.in
// try: switch between regular and cart packing. Cart is helpful
// to confirm implementations based on maple worksheets mentioned below.   
//
// Conclusions: initial integration of (p(x,y,z)-P)^2 over the whole
// voronoi cvs were not k-exact, even in the case of uniform packings.
// The squaring of the distance from the presumed piecewise constant 
// solution sets up some subtle competition that pulls the solution 
// away from the k-exact linear solution.
// 
// I then tried a slice-based integration, where the volume made up
// of the face polygon and its 2 forming points (P and NBR for example)
// is integrated for each internal face. This seemed to work in 2d, recovering exact 
// linear gradients in simple cases, but introduced a factor of 25/26 (cartesian), or 15/14 
// (hcp) wrt exact gradients in perfect cases. This was reproduced (for Cartesian)
// in the family of maple worksheets new_polynomial_3d_fa_lambda*mw.
//
// At this point, I think the way forward is to use the constraints to force
// the polynomial to formally respect the solution in the neighborhood of
// the integration region, and then use the (p(x,y,z)-P)^2 integration
// to regularize unclosed coefficients. Unfortunately, this approach will not
// work in the neighborhood of the boundary? 
//
// F. Ham April 2, 2019
// ======================================================================

#define SR 0 // holds the length scale that normalizes all x's

#define SW 1 // volume (*6?)

#define SWX 2
#define SWY 3
#define SWZ 4

#define SWXX 5
#define SWXY 6
#define SWXZ 7
#define SWYY 8
#define SWYZ 9
#define SWZZ 10

#define SWXXX 11
#define SWXXY 12
#define SWXXZ 13
#define SWXYY 14
#define SWXYZ 15
#define SWXZZ 16
#define SWYYY 17
#define SWYYZ 18
#define SWYZZ 19
#define SWZZZ 20

#define SWXXXX 21
#define SWXXXY 22
#define SWXXXZ 23
#define SWXXYY 24
#define SWXXYZ 25
#define SWXXZZ 26
#define SWXYYY 27
#define SWXYYZ 28
#define SWXYZZ 29
#define SWXZZZ 30
#define SWYYYY 31
#define SWYYYZ 32
#define SWYYZZ 33
#define SWYZZZ 34
#define SWZZZZ 35

inline void addTetQuadrature(double sums[36],const double x0[3],const double x1[3],const double x2[3],const double x3[3],const double delta) {

  using GaussQuadrature::tet4;
  using GaussQuadrature::tet10;
  using GaussQuadrature::tet20;

  const double this_vol = SIGNED_TET_VOLUME_6(x0,x1,x2,x3);
  sums[SW] += this_vol;

  // the linear terms can be computed at the mid-point...
  sums[SWX] += this_vol*0.25*(x0[0]+x1[0]+x2[0]+x3[0])/delta;
  sums[SWY] += this_vol*0.25*(x0[1]+x1[1]+x2[1]+x3[1])/delta;
  sums[SWZ] += this_vol*0.25*(x0[2]+x1[2]+x2[2]+x3[2])/delta;
  
  // the quadratic and above terms require quadrature...
  double xp[3]; 
  for (int ip = 0; ip < 4; ++ip) {
    FOR_J3 xp[j] = (tet4[ip][0]*x0[j] + tet4[ip][1]*x1[j] + tet4[ip][2]*x2[j] + tet4[ip][3]*x3[j])/delta;
    // tet4[ip][4] is the weight... 
    // diagonal...
    sums[SWXX] += this_vol*tet4[ip][4]*xp[0]*xp[0];
    sums[SWYY] += this_vol*tet4[ip][4]*xp[1]*xp[1];
    sums[SWZZ] += this_vol*tet4[ip][4]*xp[2]*xp[2];
    // offd...
    sums[SWYZ] += this_vol*tet4[ip][4]*xp[1]*xp[2];
    sums[SWXZ] += this_vol*tet4[ip][4]*xp[0]*xp[2];
    sums[SWXY] += this_vol*tet4[ip][4]*xp[0]*xp[1];
  }
  
  for (int ip = 0; ip < 10; ++ip) {
    FOR_J3 xp[j] = (tet10[ip][0]*x0[j] + tet10[ip][1]*x1[j] + tet10[ip][2]*x2[j] + tet10[ip][3]*x3[j])/delta;
    sums[SWXXX] += this_vol*tet10[ip][4]*xp[0]*xp[0]*xp[0];
    sums[SWXXY] += this_vol*tet10[ip][4]*xp[0]*xp[0]*xp[1];
    sums[SWXXZ] += this_vol*tet10[ip][4]*xp[0]*xp[0]*xp[2];
    sums[SWXYY] += this_vol*tet10[ip][4]*xp[0]*xp[1]*xp[1];
    sums[SWXYZ] += this_vol*tet10[ip][4]*xp[0]*xp[1]*xp[2];
    sums[SWXZZ] += this_vol*tet10[ip][4]*xp[0]*xp[2]*xp[2];
    sums[SWYYY] += this_vol*tet10[ip][4]*xp[1]*xp[1]*xp[1];
    sums[SWYYZ] += this_vol*tet10[ip][4]*xp[1]*xp[1]*xp[2];
    sums[SWYZZ] += this_vol*tet10[ip][4]*xp[1]*xp[2]*xp[2];
    sums[SWZZZ] += this_vol*tet10[ip][4]*xp[2]*xp[2]*xp[2];
  }
  
  for (int ip = 0; ip < 20; ++ip) {
    FOR_J3 xp[j] = (tet20[ip][0]*x0[j] + tet20[ip][1]*x1[j] + tet20[ip][2]*x2[j] + tet20[ip][3]*x3[j])/delta;
    sums[SWXXXX] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[0]*xp[0];
    sums[SWXXXY] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[0]*xp[1];
    sums[SWXXXZ] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[0]*xp[2];
    sums[SWXXYY] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[1]*xp[1];
    sums[SWXXYZ] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[1]*xp[2];
    sums[SWXXZZ] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[2]*xp[2];
    sums[SWXYYY] += this_vol*tet20[ip][4]*xp[0]*xp[1]*xp[1]*xp[1];
    sums[SWXYYZ] += this_vol*tet20[ip][4]*xp[0]*xp[1]*xp[1]*xp[2];
    sums[SWXYZZ] += this_vol*tet20[ip][4]*xp[0]*xp[1]*xp[2]*xp[2];
    sums[SWXZZZ] += this_vol*tet20[ip][4]*xp[0]*xp[2]*xp[2]*xp[2];
    sums[SWYYYY] += this_vol*tet20[ip][4]*xp[1]*xp[1]*xp[1]*xp[1];
    sums[SWYYYZ] += this_vol*tet20[ip][4]*xp[1]*xp[1]*xp[1]*xp[2];
    sums[SWYYZZ] += this_vol*tet20[ip][4]*xp[1]*xp[1]*xp[2]*xp[2];
    sums[SWYZZZ] += this_vol*tet20[ip][4]*xp[1]*xp[2]*xp[2]*xp[2];
    sums[SWZZZZ] += this_vol*tet20[ip][4]*xp[2]*xp[2]*xp[2]*xp[2];
  }
}

inline void transformQuadratureTerms(double sums[36],const double sums0[36],const double dx_actual[3],const double r_vv) {

  // sums0[36] contains the quadrature terms computed about the voronoi seed point of the 
  // cv in question. These will be normalized wrt their own r_vv stored in sums0[0].
  // this routine returns the quadrature terms computed as if there were an actual displacement
  // of dx_actual of the voronoi seed point (e.g. the nbr distance) and normalization was
  // wrt r_vv, which will not neccessarily equal the current normalization in sums0[0].
  
  sums[SR] = r_vv; // sums0[SR]; // length scale
  
  sums[SW] = sums0[SW]; // volume

  double dx[3]; FOR_I3 dx[i] = dx_actual[i]/sums0[0];
  const double ratio = sums0[0]/r_vv;
  
  sums[SWX] = ( dx[0]*sums0[SW] +
                sums0[SWX] )*ratio;
      
  sums[SWY] = ( dx[1]*sums0[SW] +
                sums0[SWY] )*ratio;
      
  sums[SWZ] = ( dx[2]*sums0[SW] +
                sums0[SWZ] )*ratio;

  const double ratio2 = ratio*ratio;
  
  sums[SWXX] = ( dx[0]*dx[0]*sums0[SW] +
                 dx[0]*sums0[SWX]*2.0 +
                 sums0[SWXX] )*ratio2;
    
  sums[SWXY] = ( dx[0]*dx[1]*sums0[SW] +
                 dx[0]*sums0[SWY] +
                 dx[1]*sums0[SWX] +
                 sums0[SWXY] )*ratio2;
  
  sums[SWXZ] = ( dx[0]*dx[2]*sums0[SW] +
                 dx[0]*sums0[SWZ] +
                 dx[2]*sums0[SWX] +
                 sums0[SWXZ] )*ratio2;

  sums[SWYY] = ( dx[1]*dx[1]*sums0[SW] +
                 dx[1]*sums0[SWY]*2.0 +
                 sums0[SWYY] )*ratio2;

  sums[SWYZ] = ( dx[1]*dx[2]*sums0[SW] +
                 dx[1]*sums0[SWZ] +
                 dx[2]*sums0[SWY] +
                 sums0[SWYZ] )*ratio2;

  sums[SWZZ] = ( dx[2]*dx[2]*sums0[SW] +
                 dx[2]*sums0[SWZ]*2.0 +
                 sums0[SWZZ] )*ratio2;

  const double ratio3 = ratio2*ratio;
    
  sums[SWXXX] = ( dx[0]*dx[0]*dx[0]*sums0[SW] +
                  dx[0]*dx[0]*sums0[SWX]*3.0 +
                  dx[0]*sums0[SWXX]*3.0 +
                  sums0[SWXXX] )*ratio3;

  sums[SWXXY] = ( dx[0]*dx[0]*dx[1]*sums0[SW] +
                  dx[0]*dx[0]*sums0[SWY] +
                  dx[0]*dx[1]*sums0[SWX]*2.0 +
                  dx[0]*sums0[SWXY]*2.0 +
                  dx[1]*sums0[SWXX] +
                  sums0[SWXXY] )*ratio3;
    
  sums[SWXXZ] = ( dx[0]*dx[0]*dx[2]*sums0[SW] +
                  dx[0]*dx[0]*sums0[SWZ] +
                  dx[0]*dx[2]*sums0[SWX]*2.0 +
                  dx[0]*sums0[SWXZ]*2.0 +
                  dx[2]*sums0[SWXX] +
                  sums0[SWXXZ] )*ratio3;
    
  sums[SWXYY] = ( dx[0]*dx[1]*dx[1]*sums0[SW] +
                  dx[0]*dx[1]*sums0[SWY]*2.0 +
                  dx[1]*dx[1]*sums0[SWX] +
                  dx[0]*sums0[SWYY] +
                  dx[1]*sums0[SWXY]*2.0 +
                  sums0[SWXYY] )*ratio3;
    
  sums[SWXYZ] = ( dx[0]*dx[1]*dx[2]*sums0[SW] +
                  dx[0]*dx[1]*sums0[SWZ] +
                  dx[0]*dx[2]*sums0[SWY] +
                  dx[1]*dx[2]*sums0[SWX] +
                  dx[0]*sums0[SWYZ] +
                  dx[1]*sums0[SWXZ] +
                  dx[2]*sums0[SWXY] +
                  sums0[SWXYZ] )*ratio3;
  
  sums[SWXZZ] = ( dx[0]*dx[2]*dx[2]*sums0[SW] +
                  dx[0]*dx[2]*sums0[SWZ]*2.0 +
                  dx[2]*dx[2]*sums0[SWX] +
                  dx[0]*sums0[SWZZ] +
                  dx[2]*sums0[SWXZ]*2.0 +
                  sums0[SWXZZ] )*ratio3;
    
  sums[SWYYY] = ( dx[1]*dx[1]*dx[1]*sums0[SW] +
                  dx[1]*dx[1]*sums0[SWY]*3.0 +
                  dx[1]*sums0[SWYY]*3.0 +
                  sums0[SWYYY] )*ratio3;
  
  sums[SWYYZ] = ( dx[1]*dx[1]*dx[2]*sums0[SW] +
                  dx[1]*dx[1]*sums0[SWZ] +
                  dx[1]*dx[2]*sums0[SWY]*2.0 +
                  dx[1]*sums0[SWYZ]*2.0 +
                  dx[2]*sums0[SWYY] +
                  sums0[SWYYZ] )*ratio3;
  
  sums[SWYZZ] = ( dx[1]*dx[2]*dx[2]*sums0[SW] +
                  dx[1]*dx[2]*sums0[SWZ]*2.0 +
                  dx[2]*dx[2]*sums0[SWY] +
                  dx[1]*sums0[SWZZ] +
                  dx[2]*sums0[SWYZ]*2.0 +
                  sums0[SWYZZ] )*ratio3;
  
  sums[SWZZZ] = ( dx[2]*dx[2]*dx[2]*sums0[SW] +
                  dx[2]*dx[2]*sums0[SWZ]*3.0 +
                  dx[2]*sums0[SWZZ]*3.0 +
                  sums0[SWZZZ] )*ratio3;
  
  const double ratio4 = ratio2*ratio2;
  
  sums[SWXXXX] = ( dx[0]*dx[0]*dx[0]*dx[0]*sums0[SW] +
                   dx[0]*dx[0]*dx[0]*sums0[SWX]*4.0 +
                   dx[0]*dx[0]*sums0[SWXX]*6.0 +
                   dx[0]*sums0[SWXXX]*4.0 +
                   sums0[SWXXXX] )*ratio4;
  
  sums[SWXXXY] = ( dx[0]*dx[0]*dx[0]*dx[1]*sums0[SW] +
                   dx[0]*dx[0]*dx[0]*sums0[SWY] +
                   dx[0]*dx[0]*dx[1]*sums0[SWX]*3.0 +
                   dx[0]*dx[0]*sums0[SWXY]*3.0 +
                   dx[0]*dx[1]*sums0[SWXX]*3.0 +
                   dx[0]*sums0[SWXXY]*3.0 +
                   dx[1]*sums0[SWXXX] +
                   sums0[SWXXXY] )*ratio4;
  
  sums[SWXXXZ] = ( dx[0]*dx[0]*dx[0]*dx[2]*sums0[SW] +
                   dx[0]*dx[0]*dx[0]*sums0[SWZ] +
                   dx[0]*dx[0]*dx[2]*sums0[SWX]*3.0 +
                   dx[0]*dx[0]*sums0[SWXZ]*3.0 +
                   dx[0]*dx[2]*sums0[SWXX]*3.0 +
                   dx[0]*sums0[SWXXZ]*3.0 +
                   dx[2]*sums0[SWXXX] +
                   sums0[SWXXXZ] )*ratio4;
  
  sums[SWXXYY] = ( dx[0]*dx[0]*dx[1]*dx[1]*sums0[SW] +
                   dx[0]*dx[0]*dx[1]*sums0[SWY]*2.0 +
                   dx[0]*dx[1]*dx[1]*sums0[SWX]*2.0 +
                   dx[0]*dx[0]*sums0[SWYY] +
                   dx[0]*dx[1]*sums0[SWXY]*4.0 +
                   dx[1]*dx[1]*sums0[SWXX] +
                   dx[0]*sums0[SWXYY]*2.0 +
                   dx[1]*sums0[SWXXY]*2.0 +
                   sums0[SWXXYY] )*ratio4;
  
  sums[SWXXYZ] = ( dx[0]*dx[0]*dx[1]*dx[2]*sums0[SW] +
                   dx[0]*dx[0]*dx[1]*sums0[SWZ] + 
                   dx[0]*dx[0]*dx[2]*sums0[SWY] + 
                   dx[0]*dx[1]*dx[2]*sums0[SWX]*2.0 +
                   dx[0]*dx[0]*sums0[SWYZ] +
                   dx[0]*dx[1]*sums0[SWXZ]*2.0 +
                   dx[0]*dx[2]*sums0[SWXY]*2.0 +
                   dx[1]*dx[2]*sums0[SWXX] +
                   dx[0]*sums0[SWXYZ]*2.0 +
                   dx[1]*sums0[SWXXZ] +
                   dx[2]*sums0[SWXXY] +
                   sums0[SWXXYZ] )*ratio4;
  
  sums[SWXXZZ] = ( dx[0]*dx[0]*dx[2]*dx[2]*sums0[SW] +
                   dx[0]*dx[0]*dx[2]*sums0[SWZ]*2.0 +
                   dx[0]*dx[2]*dx[2]*sums0[SWX]*2.0 +
                   dx[0]*dx[0]*sums0[SWZZ] +
                   dx[0]*dx[2]*sums0[SWXZ]*4.0 +
                   dx[2]*dx[2]*sums0[SWXX] +
                   dx[0]*sums0[SWXZZ]*2.0 +
                   dx[2]*sums0[SWXXZ]*2.0 +
                   sums0[SWXXZZ] )*ratio4;
  
  sums[SWXYYY] = ( dx[0]*dx[1]*dx[1]*dx[1]*sums0[SW] +
                   dx[0]*dx[1]*dx[1]*sums0[SWY]*3.0 +
                   dx[1]*dx[1]*dx[1]*sums0[SWX] +
                   dx[1]*dx[1]*sums0[SWXY]*3.0 +
                   dx[0]*dx[1]*sums0[SWYY]*3.0 +
                   dx[1]*sums0[SWXYY]*3.0 +
                   dx[0]*sums0[SWYYY] +
                   sums0[SWXYYY] )*ratio4;
  
  sums[SWXYYZ] = ( dx[0]*dx[1]*dx[1]*dx[2]*sums0[SW] +
                   dx[0]*dx[1]*dx[1]*sums0[SWZ] + 
                   dx[0]*dx[1]*dx[2]*sums0[SWY]*2.0 + 
                   dx[1]*dx[1]*dx[2]*sums0[SWX] +
                   dx[0]*dx[1]*sums0[SWYZ]*2.0 +
                   dx[1]*dx[1]*sums0[SWXZ] +
                   dx[1]*dx[2]*sums0[SWXY]*2.0 +
                   dx[0]*dx[2]*sums0[SWYY] +
                   dx[0]*sums0[SWYYZ] +
                   dx[1]*sums0[SWXYZ]*2.0 +
                   dx[2]*sums0[SWXYY] +
                   sums0[SWXYYZ] )*ratio4;
  
  sums[SWXYZZ] = ( dx[0]*dx[1]*dx[2]*dx[2]*sums0[SW] +
                   dx[0]*dx[1]*dx[2]*sums0[SWZ]*2.0 + 
                   dx[0]*dx[2]*dx[2]*sums0[SWY] + 
                   dx[1]*dx[2]*dx[2]*sums0[SWX] +
                   dx[0]*dx[1]*sums0[SWZZ] +
                   dx[0]*dx[2]*sums0[SWYZ]*2.0 +
                   dx[1]*dx[2]*sums0[SWXZ]*2.0 +
                   dx[2]*dx[2]*sums0[SWXY] +
                   dx[0]*sums0[SWYZZ] +
                   dx[1]*sums0[SWXZZ] +
                   dx[2]*sums0[SWXYZ]*2.0 +
                   sums0[SWXYZZ] )*ratio4;
  
  sums[SWXZZZ] = ( dx[0]*dx[2]*dx[2]*dx[2]*sums0[SW] +
                   dx[0]*dx[2]*dx[2]*sums0[SWZ]*3.0 +
                   dx[2]*dx[2]*dx[2]*sums0[SWX] +
                   dx[2]*dx[2]*sums0[SWXZ]*3.0 +
                   dx[0]*dx[2]*sums0[SWZZ]*3.0 +
                   dx[0]*sums0[SWZZZ] +
                   dx[2]*sums0[SWXZZ]*3.0 +
                   sums0[SWXZZZ] )*ratio4;
  
  sums[SWYYYY] = ( dx[1]*dx[1]*dx[1]*dx[1]*sums0[SW] +
                   dx[1]*dx[1]*dx[1]*sums0[SWY]*4.0 +
                   dx[1]*dx[1]*sums0[SWYY]*6.0 +
                   dx[1]*sums0[SWYYY]*4.0 +
                   sums0[SWYYYY] )*ratio4;
  
  sums[SWYYYZ] = ( dx[1]*dx[1]*dx[1]*dx[2]*sums0[SW] +
                   dx[1]*dx[1]*dx[1]*sums0[SWZ] +
                   dx[1]*dx[1]*dx[2]*sums0[SWY]*3.0 +
                   dx[1]*dx[1]*sums0[SWYZ]*3.0 +
                   dx[1]*dx[2]*sums0[SWYY]*3.0 +
                   dx[1]*sums0[SWYYZ]*3.0 +
                   dx[2]*sums0[SWYYY] +
                   sums0[SWYYYZ] )*ratio4;
  
  sums[SWYYZZ] = ( dx[1]*dx[1]*dx[2]*dx[2]*sums0[SW] +
                   dx[1]*dx[1]*dx[2]*sums0[SWZ]*2.0 +
                   dx[1]*dx[2]*dx[2]*sums0[SWY]*2.0 +
                   dx[1]*dx[1]*sums0[SWZZ] +
                   dx[1]*dx[2]*sums0[SWYZ]*4.0 +
                   dx[2]*dx[2]*sums0[SWYY] +
                   dx[1]*sums0[SWYZZ]*2.0 +
                   dx[2]*sums0[SWYYZ]*2.0 +
                   sums0[SWYYZZ] )*ratio4;
  
  sums[SWYZZZ] = ( dx[1]*dx[2]*dx[2]*dx[2]*sums0[SW] +
                   dx[1]*dx[2]*dx[2]*sums0[SWZ]*3.0 +
                   dx[2]*dx[2]*dx[2]*sums0[SWY] +
                   dx[1]*dx[2]*sums0[SWZZ]*3.0 +
                   dx[2]*dx[2]*sums0[SWYZ]*3.0 +
                   dx[1]*sums0[SWZZZ] +
                   dx[2]*sums0[SWYZZ]*3.0 +
                   sums0[SWYZZZ] )*ratio4;
  
  sums[SWZZZZ] = ( dx[2]*dx[2]*dx[2]*dx[2]*sums0[SW] +
                   dx[2]*dx[2]*dx[2]*sums0[SWZ]*4.0 +
                   dx[2]*dx[2]*sums0[SWZZ]*6.0 +
                   dx[2]*sums0[SWZZZ]*4.0 +
                   sums0[SWZZZZ] )*ratio4;
  
}

#define L(i,j) L[(j)*10 + i]
#define A(i,j) A[(j)*10 + i]

inline void addSumsToAAndRhs(double * A,double (*rhs)[10],const double sums[36],const int inbr) {
  // c0...
  rhs[inbr][0] += sums[SW];
  A(0,0) += sums[SW];
  A(0,1) += sums[SWX];
  A(0,2) += sums[SWY];
  A(0,3) += sums[SWZ];
  A(0,4) += sums[SWXX];
  A(0,5) += sums[SWXY];
  A(0,6) += sums[SWXZ];
  A(0,7) += sums[SWYY];
  A(0,8) += sums[SWYZ];
  A(0,9) += sums[SWZZ];
  // c1...
  rhs[inbr][1] += sums[SWX];
  A(1,0) += sums[SWX];
  A(1,1) += sums[SWXX];
  A(1,2) += sums[SWXY];
  A(1,3) += sums[SWXZ];
  A(1,4) += sums[SWXXX];
  A(1,5) += sums[SWXXY];
  A(1,6) += sums[SWXXZ];
  A(1,7) += sums[SWXYY];
  A(1,8) += sums[SWXYZ];
  A(1,9) += sums[SWXZZ];
  // c2...
  rhs[inbr][2] += sums[SWY];
  A(2,0) += sums[SWY];
  A(2,1) += sums[SWXY];
  A(2,2) += sums[SWYY];
  A(2,3) += sums[SWYZ];
  A(2,4) += sums[SWXXY];
  A(2,5) += sums[SWXYY];
  A(2,6) += sums[SWXYZ];
  A(2,7) += sums[SWYYY];
  A(2,8) += sums[SWYYZ];
  A(2,9) += sums[SWYZZ];
  // c3...
  rhs[inbr][3] += sums[SWZ];
  A(3,0) += sums[SWZ];
  A(3,1) += sums[SWXZ];
  A(3,2) += sums[SWYZ];
  A(3,3) += sums[SWZZ];
  A(3,4) += sums[SWXXZ];
  A(3,5) += sums[SWXYZ];
  A(3,6) += sums[SWXZZ];
  A(3,7) += sums[SWYYZ];
  A(3,8) += sums[SWYZZ];
  A(3,9) += sums[SWZZZ];
  // c4...
  rhs[inbr][4] += sums[SWXX];
  A(4,0) += sums[SWXX];
  A(4,1) += sums[SWXXX];
  A(4,2) += sums[SWXXY];
  A(4,3) += sums[SWXXZ];
  A(4,4) += sums[SWXXXX];
  A(4,5) += sums[SWXXXY];
  A(4,6) += sums[SWXXXZ];
  A(4,7) += sums[SWXXYY];
  A(4,8) += sums[SWXXYZ];
  A(4,9) += sums[SWXXZZ];
  // c5...
  rhs[inbr][5] += sums[SWXY];
  A(5,0) += sums[SWXY];
  A(5,1) += sums[SWXXY];
  A(5,2) += sums[SWXYY];
  A(5,3) += sums[SWXYZ];
  A(5,4) += sums[SWXXXY];
  A(5,5) += sums[SWXXYY];
  A(5,6) += sums[SWXXYZ];
  A(5,7) += sums[SWXYYY];
  A(5,8) += sums[SWXYYZ];
  A(5,9) += sums[SWXYZZ];
  // c6...
  rhs[inbr][6] += sums[SWXZ];
  A(6,0) += sums[SWXZ];
  A(6,1) += sums[SWXXZ];
  A(6,2) += sums[SWXYZ];
  A(6,3) += sums[SWXZZ];
  A(6,4) += sums[SWXXXZ];
  A(6,5) += sums[SWXXYZ];
  A(6,6) += sums[SWXXZZ];
  A(6,7) += sums[SWXYYZ];
  A(6,8) += sums[SWXYZZ];
  A(6,9) += sums[SWXZZZ];
  // c7...
  rhs[inbr][7] += sums[SWYY];
  A(7,0) += sums[SWYY];
  A(7,1) += sums[SWXYY];
  A(7,2) += sums[SWYYY];
  A(7,3) += sums[SWYYZ];
  A(7,4) += sums[SWXXYY];
  A(7,5) += sums[SWXYYY];
  A(7,6) += sums[SWXYYZ];
  A(7,7) += sums[SWYYYY];
  A(7,8) += sums[SWYYYZ];
  A(7,9) += sums[SWYYZZ];
  // c8...
  rhs[inbr][8] += sums[SWYZ];
  A(8,0) += sums[SWYZ];
  A(8,1) += sums[SWXYZ];
  A(8,2) += sums[SWYYZ];
  A(8,3) += sums[SWYZZ];
  A(8,4) += sums[SWXXYZ];
  A(8,5) += sums[SWXYYZ];
  A(8,6) += sums[SWXYZZ];
  A(8,7) += sums[SWYYYZ];
  A(8,8) += sums[SWYYZZ];
  A(8,9) += sums[SWYZZZ];
  // c9...
  rhs[inbr][9] += sums[SWZZ];
  A(9,0) += sums[SWZZ];
  A(9,1) += sums[SWXZZ];
  A(9,2) += sums[SWYZZ];
  A(9,3) += sums[SWZZZ];
  A(9,4) += sums[SWXXZZ];
  A(9,5) += sums[SWXYZZ];
  A(9,6) += sums[SWXZZZ];
  A(9,7) += sums[SWYYZZ];
  A(9,8) += sums[SWYZZZ];
  A(9,9) += sums[SWZZZZ];
}
    
inline void cholesky_fac_n10(double* L, const double* A) { 

  // for now , we're not going to overwrite A
  // but we can always change that later..
  for (int i =0; i < 100; ++i) 
    L[i] = A[i];

  for (int j = 0; j < 10; ++j) { 
    for (int k = 0; k < j; ++k) 
      for (int i = 0; i < 10; ++i) 
        L(i,j) -= L(i,k)*L(j,k);
    
    L(j,j) = sqrt(L(j,j));
    for (int k = j+1; k < 10; ++k) 
      L(k,j) /= L(j,j); 
  }
    
}

inline void cholesky_solve_n10(double* r, const double *L) { 
  // solve LL^T[x] = r

  // solve Ly = r
  for(int i = 0; i < 10; ++i) {
    for(int j = 0; j < i; ++j) {                                      
      r[i] -= L(i,j)*r[j]; 
    }
    r[i] = r[i]/L(i,i); 
  }

  // solve L^Tx = y, overwrite r..
  for (int i = 9; i >= 0; --i) { 
    for (int j = 9; j > i; --j) 
      r[i] -= L(j,i)*r[j];
    r[i] = r[i]/L(i,i);
  }

}

void VoronoiBuilder::processComputeNewOps(Param * param,const bool b_help) {

  if (b_help) {
    WUI(INFO,"COMPUTE_NEW_OPS computes accurate operators base on finite volume integration");
    return;
  }
  
  // need to build voronoi first...
  // TODO: just NSMOOTH 0 to 
  if (!PartData::b_vp) {
    WUI(WARN,"expecting NSMOOTH <int> prior to COMPUTE_NEW_OPS");
    return;
  }

  const double x_debug[3] = { 0.15, 0, -0.519615 }; // for regular packing
  //const double x_debug[3] = { 0.15, 0, -0.525 }; // for cart packing
  const double phi0_exact = 0.456;
  //const double phi_grad_exact[3] = { 1.123, 2.765, -1.987 };
  const double phi_grad_exact[3] = { 1.123, 1.123, 1.123 };

  double (*firstSums)[5] = new double[points->np][5];
  int nnbr_max = 0;
  double dx[3] = { 0.0, 0.0, 0.0 }; // HACK -- get rid of this once everything has been checked
  for (int ip = 0; ip < points->np; ++ip) {
    calcQuadratureLengthScale(firstSums[ip][0],ip);
    const int nnbr = calcFirstQuadratureTerms(firstSums[ip],ip,dx,firstSums[ip][0]); // delete dx, r_vv eventually
    nnbr_max = max(nnbr_max,nnbr);
  }
  
  cout << " > got nnbr_max: " << nnbr_max << endl;
  
  // the centroids are going to be helpful for checking the operators. Produce them using the
  // firstSums...
  
  double (*x_cv)[3] = new double[points->np][3];
  double *phi = new double[points->np];
  for (int ip = 0; ip < points->np; ++ip) {
    x_cv[ip][0] = points->xp[ip][0] + firstSums[ip][SWX]/firstSums[ip][SW]*firstSums[ip][0];
    x_cv[ip][1] = points->xp[ip][1] + firstSums[ip][SWY]/firstSums[ip][SW]*firstSums[ip][0];
    x_cv[ip][2] = points->xp[ip][2] + firstSums[ip][SWZ]/firstSums[ip][SW]*firstSums[ip][0];
    phi[ip] = phi0_exact + DOT_PRODUCT(x_cv[ip],phi_grad_exact);
  }
  
  GeomUtils::writePtsTecplot("x_cv.dat",x_cv,points->np);
  cout << "take a look at x_cv.dat" << endl;
  
  vector<VdIndexAndTri> vdFaTriVec;
  vector<VdTri> vdBfTriVec;
  vector<uint8> rbiNbrVec;
  
  double sums0[36],sums1[36];
  double A[100],L[100];
  double r[10],r_check[10];
  int * iponb = new int[nnbr_max];
  double (*rhs)[10] = new double[nnbr_max+1][10]; // need space for ip
  const double zero[3] = { 0.0, 0.0, 0.0 }; 
  
  for (int ip = 0; ip < points->np; ++ip) {

    // select debug...
    bool debug = false;
    if (DIST(points->xp[ip],x_debug) < 1.0E-5) {
      const double dx_cv[3] = DIFF(x_cv[ip],points->xp[ip]);
      cout << "got debug: ip: " << ip << " dx_cv: " << COUT_VEC(dx_cv) << endl;
      debug = true;
    }
    
    // this routine returns the vdFaTriVec sorted by associated nbr rbi...
    buildTriVecsAndNbrs(vdFaTriVec,vdBfTriVec,rbiNbrVec,ip);
    const int nnbr = rbiNbrVec.size(); assert(nnbr <= nnbr_max);
    
    // the normalization length scale has already been computed... 
    double r_vv = firstSums[ip][0];
    
    // zero everything...
    for (int inbr = 0; inbr <= nnbr; ++inbr) 
      for (int i = 0; i < 10; ++i)
        rhs[inbr][i] = 0.0;
    for (int i = 0; i < 100; ++i)
      A[i] = 0.0;
    
    int inbr = 0;
    for (int i = 0; i < 36; ++i) sums0[i] = sums1[i] = 0.0;
    for (int ii = 0; ii < vdFaTriVec.size(); ++ii) {
      while (vdFaTriVec[ii].index != rbiNbrVec[inbr]) {
        cout << "adding sums to nbr: " << inbr << " SW: " << sums0[SW] << " " << sums1[SW] << endl;
        addSumsToAAndRhs(A,rhs,sums0,0);
        addSumsToAAndRhs(A,rhs,sums1,inbr+1); // pass inbr+1 here so 0 is "P"...
        for (int i = 0; i < 36; ++i) sums0[i] = sums1[i] = 0.0;
        ++inbr;
      }
      addTetQuadrature(sums0,zero,vdFaTriVec[ii].x[0],vdFaTriVec[ii].x[1],vdFaTriVec[ii].x[2],r_vv);
      int rank,bits,ip_nbr;
      BitUtils::unpackRankBitsIndex(rank,bits,ip_nbr,vdFaTriVec[ii].index);
      assert(rank == 0);
      assert(bits == 0);
      assert((ip_nbr >= 0)&&(ip_nbr < points->np)); // local for now
      assert(ip_nbr != ip);
      // store this for use below...
      iponb[inbr] = ip_nbr;
      // one of the vertices of the tet is the nbr's forming point. Note we permute
      // the order of the integration points to get the sign of the tet correct...
      const double dx[3] = DIFF(points->xp[ip_nbr],points->xp[ip]);
      addTetQuadrature(sums1,vdFaTriVec[ii].x[0],vdFaTriVec[ii].x[1],vdFaTriVec[ii].x[2],dx,r_vv);
    }
    // and add the final sums...
    cout << "adding sums to nbr: " << inbr << " SW: " << sums0[SW] << " " << sums1[SW] << endl;
    addSumsToAAndRhs(A,rhs,sums0,0);
    addSumsToAAndRhs(A,rhs,sums1,inbr+1);
    
    // symmetry check...
    
    cout << "symmetry check...";
    for (int i = 0; i < 10; ++i) {
      for (int j = i+1; j < 10; ++j) {
        assert(A(i,j) == A(j,i));
      }
    }
    cout << "OK" << endl;

    // cholesky factor A...

    cholesky_fac_n10(L,A);
    
    // the full rhs for this polynomial is then...
    
    for (int i = 0; i < 10; ++i) {
      r[i] = rhs[0][i]*phi[ip];
      for (int inbr = 0; inbr < nnbr; ++inbr) {
        const int ip_nbr = iponb[inbr];
        assert(ip_nbr != ip);
        r[i] += rhs[inbr+1][i]*phi[ip_nbr];
      }
    }
    
    for (int i = 0; i < 10; ++i) 
      r_check[i] = r[i];

    cholesky_solve_n10(r,L);

    cout << "cholesky_solve_n10: " << endl;
    for (int i = 0; i < 10; ++i) {
      cout << "r[" << i << "] = " << r[i] << endl;
    }
    cout << " r[1]/r_vv: " << r[1]/r_vv << " vs grad: " << r[1]/r_vv/phi_grad_exact[0] << endl;
    cout << " r[2]/r_vv: " << r[2]/r_vv << " vs grad: " << r[2]/r_vv/phi_grad_exact[1] << endl;
    cout << " r[3]/r_vv: " << r[3]/r_vv << " vs grad: " << r[3]/r_vv/phi_grad_exact[2] << endl;
    
    cout << "phi at points->xp: " << phi0_exact + DOT_PRODUCT(points->xp[ip],phi_grad_exact) << " phi(x_cv): " << phi[ip] << " nnbr: " << nnbr << endl;

    // check...
    
    cout << "check: " << endl;
    for (int i = 0; i < 10; ++i) {
      double this_r = 0.0;
      for (int j = 0; j < 10; ++j) {
        this_r += A(i,j)*r[j];
      }
      cout << " i: " << i << " this_r: " << this_r << " r_check[i]: " << r_check[i] << " diff: " << this_r-r_check[i] << endl;
      assert(fabs(this_r-r_check[i]) < 1.0E-15);
    }
    
    cout << "done check" << endl;
    
    if (debug) 
      getchar();

    vdFaTriVec.clear();
    vdBfTriVec.clear();
    rbiNbrVec.clear();
    
  }

  MPI_Pause("SOOOOOOO FAR SOOOOOOOOO GOOOOOOOOOOD");
  
}

#undef L
#undef A

void VoronoiBuilder::calcQuadratureLengthScale(double& r_vv,const int ip) {

  // the quadrature length scale for a given vd is the maximum node distance
  // of the vd (including orphans) from the forming point. put this in sums[0]...

  // these integrations must be normalized by a length scale. Use r_vv...
  r_vv = 0.0;
  for (int ied = 0; ied < vdArray[ip].getNedMain(); ++ied) {
    FOR_I2 {
      double * x = vdArray[ip].getXnoPtr(vdArray[ip].getNooed(ied,i));
      r_vv = max(r_vv,DOT_PRODUCT(x,x));
    }
  }
  int next = vdArray[ip].orphan_chunk_data; // -1 or the index in the ocdVec of our first orphan chunk
  while (next != -1) {
    // edge loop on ocdVec[next]...
    for (int ied = 0; ied < ocdVec[next].getNed(); ++ied) {
      FOR_I2 {
        double * x = ocdVec[next].getXnoPtr(ocdVec[next].getNooed(ied,i));
        r_vv = max(r_vv,DOT_PRODUCT(x,x));
      }
    }
    next = ocdVec[next].next;
  }
  // and take the sqrt...
  r_vv = sqrt(r_vv);
  
}

void VoronoiBuilder::calcQuadratureTermsAndNbrs(double sums[36],vector<uint8>& rbiNbrVec,const int ip,const double dx[3],const double r_vv) {

  // coming into this routine, sums[0] contains ip's r_vv. However the passed r_vv is
  // considered separate to allow us to check...

  map<const uint8,double*> internalFaceGeometryMap; // map of face geometry with rbi key
  map<const pair<pair<int,int>,int>,double*> boundaryFaceGeometryMap; // map of face geometry with ((ipart,ist),bits) key
    
  double gcl[3] = { 0.0, 0.0, 0.0 };
    
  using GaussQuadrature::tet4;
  using GaussQuadrature::tet10;
  using GaussQuadrature::tet20;
  
  for (int i = 1; i < 36; ++i) sums[i] = 0.0;
    
  // loop through main edges (i.e. no orphans)...
  for (int ied = 0; ied < vdArray[ip].getNedMain(); ++ied) {
    // each edge has 2 faces...
    FOR_I2 {
      const int ibf_or_ifa = vdArray[ip].getFaoed(ied,i);
      // a negative ibf_or_ifa indicates an internal nbr (positive is a boundary)...
      if (ibf_or_ifa < 0) {
        // recall vdArray[ip] uses -1 indexing for faces...
        assert((-ibf_or_ifa-1 >= 0)&&(-ibf_or_ifa-1 < vdArray[ip].getNfa()));
        // confirm that this is a group0 face. It has to be because
        // we are looping on group0 edges...
        // recall that faces have been paired in Step2, and active faces are paired faces where
        // neither face is a "zero" face (i.e. its area is below a certain tolerance)...
        // Recall "group" refers to the orphan indexing, with the main part always
        // being group == 0. Orphans are then group=1 (or group=2, etc if more than 1 orphan associated
        // with the cell), in no particular order.
        if (vdArray[ip].faceIsActive(-ibf_or_ifa-1)) {
          int group,rank,bits,index;
          vdArray[ip].getGrbiForFace(group,rank,bits,index,-ibf_or_ifa-1);
          assert(group == 0); // must be group 0 -- this is a check.
          int group_nbr,ifa_nbr;
          vdArray[ip].getGiForFaceNbr(group_nbr,ifa_nbr,-ibf_or_ifa-1);
          if (group_nbr == 0) {
            // this is a boundary face that should probably be included, unless it shares
            // its edge with a non-group 0 face nbr...
            if (vdArray[ip].getFaoed(ied,1-i) < 0) {
              const int ifa_opp = -vdArray[ip].getFaoed(ied,1-i)-1;
              assert((ifa_opp >= 0)&&(ifa_opp < vdArray[ip].getNfa()));
              if (vdArray[ip].faceIsActive(ifa_opp)) {
                int group_nbr,ifa_nbr;
                vdArray[ip].getGiForFaceNbr(group_nbr,ifa_nbr,ifa_opp);
                if (group_nbr != 0)
                  continue;
              }
            }
            // this is an active face edge - add it to the faceMap...
            const uint8 rbiHash = BitUtils::packRankBitsIndex(rank,bits,index); // rbi uniquely represents the nbr we have been cut against
            // recall i = 0,1, so edge goes form node 0->1 or 1->0, depending on the face
            double * x0 = vdArray[ip].getXnoPtr(vdArray[ip].getNooed(ied,i)); // the coordinates of this node would be x0[0],x0[1],x0[2].
            double * x1 = vdArray[ip].getXnoPtr(vdArray[ip].getNooed(ied,1-i));
            //assert(faMap.find(pair<uint8,double*>(rbiHash,x0)) == faMap.end());
            //faMap[pair<uint8,double*>(rbiHash,x0)] = x1;
            // also, build the geometric stuff...
            map<const uint8,double*>::iterator iter = internalFaceGeometryMap.find(rbiHash);
            if (iter == internalFaceGeometryMap.end()) {
              internalFaceGeometryMap[rbiHash] = x0;
            }
            else if ((x0 != iter->second)&&(x1 != iter->second)) {
              assert(iter->second != NULL);
              const double this_n[3] = TRI_NORMAL_2(iter->second,x0,x1);
              FOR_J3 gcl[j] += this_n[j];
              const double this_vol = DOT_PRODUCT(iter->second,this_n);
              sums[SW] += this_vol;
              // the linear terms can be computed at the mid-point...
              sums[SWX] += this_vol*(dx[0]+0.25*(iter->second[0]+x0[0]+x1[0]))/r_vv;
              sums[SWY] += this_vol*(dx[1]+0.25*(iter->second[1]+x0[1]+x1[1]))/r_vv;
              sums[SWZ] += this_vol*(dx[2]+0.25*(iter->second[2]+x0[2]+x1[2]))/r_vv;
              // the quadratic and above terms require quadrature...
              double xp[3]; 
              for (int ip = 0; ip < 4; ++ip) {
                // note that the first point is the forming point, zero by definition
                FOR_J3 xp[j] = (dx[j] + tet4[ip][1]*iter->second[j] + tet4[ip][2]*x0[j] + tet4[ip][3]*x1[j])/r_vv;
                // tet4[ip][4] is the weight... 
                // diagonal...
                sums[SWXX] += this_vol*tet4[ip][4]*xp[0]*xp[0];
                sums[SWYY] += this_vol*tet4[ip][4]*xp[1]*xp[1];
                sums[SWZZ] += this_vol*tet4[ip][4]*xp[2]*xp[2];
                // offd...
                sums[SWYZ] += this_vol*tet4[ip][4]*xp[1]*xp[2];
                sums[SWXZ] += this_vol*tet4[ip][4]*xp[0]*xp[2];
                sums[SWXY] += this_vol*tet4[ip][4]*xp[0]*xp[1];
              }
              for (int ip = 0; ip < 10; ++ip) {
                FOR_J3 xp[j] = (dx[j] + tet10[ip][1]*iter->second[j] + tet10[ip][2]*x0[j] + tet10[ip][3]*x1[j])/r_vv;
                sums[SWXXX] += this_vol*tet10[ip][4]*xp[0]*xp[0]*xp[0];
                sums[SWXXY] += this_vol*tet10[ip][4]*xp[0]*xp[0]*xp[1];
                sums[SWXXZ] += this_vol*tet10[ip][4]*xp[0]*xp[0]*xp[2];
                sums[SWXYY] += this_vol*tet10[ip][4]*xp[0]*xp[1]*xp[1];
                sums[SWXYZ] += this_vol*tet10[ip][4]*xp[0]*xp[1]*xp[2];
                sums[SWXZZ] += this_vol*tet10[ip][4]*xp[0]*xp[2]*xp[2];
                sums[SWYYY] += this_vol*tet10[ip][4]*xp[1]*xp[1]*xp[1];
                sums[SWYYZ] += this_vol*tet10[ip][4]*xp[1]*xp[1]*xp[2];
                sums[SWYZZ] += this_vol*tet10[ip][4]*xp[1]*xp[2]*xp[2];
                sums[SWZZZ] += this_vol*tet10[ip][4]*xp[2]*xp[2]*xp[2];
              }
              for (int ip = 0; ip < 20; ++ip) {
                FOR_J3 xp[j] = (dx[j] + tet20[ip][1]*iter->second[j] + tet20[ip][2]*x0[j] + tet20[ip][3]*x1[j])/r_vv;
                sums[SWXXXX] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[0]*xp[0];
                sums[SWXXXY] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[0]*xp[1];
                sums[SWXXXZ] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[0]*xp[2];
                sums[SWXXYY] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[1]*xp[1];
                sums[SWXXYZ] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[1]*xp[2];
                sums[SWXXZZ] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[2]*xp[2];
                sums[SWXYYY] += this_vol*tet20[ip][4]*xp[0]*xp[1]*xp[1]*xp[1];
                sums[SWXYYZ] += this_vol*tet20[ip][4]*xp[0]*xp[1]*xp[1]*xp[2];
                sums[SWXYZZ] += this_vol*tet20[ip][4]*xp[0]*xp[1]*xp[2]*xp[2];
                sums[SWXZZZ] += this_vol*tet20[ip][4]*xp[0]*xp[2]*xp[2]*xp[2];
                sums[SWYYYY] += this_vol*tet20[ip][4]*xp[1]*xp[1]*xp[1]*xp[1];
                sums[SWYYYZ] += this_vol*tet20[ip][4]*xp[1]*xp[1]*xp[1]*xp[2];
                sums[SWYYZZ] += this_vol*tet20[ip][4]*xp[1]*xp[1]*xp[2]*xp[2];
                sums[SWYZZZ] += this_vol*tet20[ip][4]*xp[1]*xp[2]*xp[2]*xp[2];
                sums[SWZZZZ] += this_vol*tet20[ip][4]*xp[2]*xp[2]*xp[2]*xp[2];
              }
              /*
                vdFaTriVec.push_back(VdTri(iter->second,x0,x1));
                iter->second.area += this_vol;
                FOR_J3 {
                const double tmp = this_vol*(iter->second[j]+x0[j]+x1[j]);
                iter->second.xc[j] += tmp;
                x_cv[ip][j] += tmp;
                iter->second.normal[j] += this_n[j];
                }
              */
            }
            else {
              // face loop should be 1-directional...
              assert(iter->second == x1);
            }
          }
        }
      }
      else {

        // a positive ibf_or_ifa is the ibf index, meaning it is an edge of a boundary face.
        // If the other face is an internal face, and
        // it has a connection with an orphan (i.e. group > 0), we skip it...

        assert((ibf_or_ifa >= 0)&&(ibf_or_ifa < vdArray[ip].getNbf()));
        int ipart,ist,bits;
        vdArray[ip].getPartStAndBitsOfBf(ipart,ist,bits,ibf_or_ifa);
          
        //map<const pair<pair<int,int>,int>,int>::iterator iter0 = colorMap.end();
        if (vdArray[ip].getFaoed(ied,1-i) < 0) {
          const int ifa_opp = -vdArray[ip].getFaoed(ied,1-i)-1;
          assert((ifa_opp >= 0)&&(ifa_opp < vdArray[ip].getNfa()));
          if (vdArray[ip].faceIsActive(ifa_opp)) {
            int group_nbr,ifa_nbr;
            vdArray[ip].getGiForFaceNbr(group_nbr,ifa_nbr,ifa_opp);
            if (group_nbr != 0) {
              continue;
            }
          }
        }
        // recall i = 0,1, so edge goes form node 0->1 or 1->0, depending on the face
        double * x0 = vdArray[ip].getXnoPtr(vdArray[ip].getNooed(ied,i)); // the coordinates of this node would be x0[0],x0[1],x0[2].
        double * x1 = vdArray[ip].getXnoPtr(vdArray[ip].getNooed(ied,1-i));
        // geometry...
        map<const pair<pair<int,int>,int>,double*>::iterator iter =
          boundaryFaceGeometryMap.find(pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits));
        if (iter == boundaryFaceGeometryMap.end()) {
          boundaryFaceGeometryMap[pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits)] = x0;
        }
        else if ((x0 != iter->second)&&(x1 != iter->second)) {
          assert(iter->second != NULL);
          const double this_n[3] = TRI_NORMAL_2(iter->second,x0,x1);
          FOR_J3 gcl[j] += this_n[j];
          const double this_vol = DOT_PRODUCT(iter->second,this_n);
          sums[SW] += this_vol;
          // the linear terms sums[1,2,3] can be computed at the mid-point...
          sums[SWX] += this_vol*(dx[0]+0.25*(iter->second[0]+x0[0]+x1[0]))/r_vv;
          sums[SWY] += this_vol*(dx[1]+0.25*(iter->second[1]+x0[1]+x1[1]))/r_vv;
          sums[SWZ] += this_vol*(dx[2]+0.25*(iter->second[2]+x0[2]+x1[2]))/r_vv;
          // the quadratic and above terms require quadrature...
          double xp[3]; 
          for (int ip = 0; ip < 4; ++ip) {
            // note that the first point is the forming point, zero by definition
            FOR_J3 xp[j] = (dx[j] + tet4[ip][1]*iter->second[j] + tet4[ip][2]*x0[j] + tet4[ip][3]*x1[j])/r_vv;
            // tet4[ip][4] is the weight... 
            // diagonal...
            sums[SWXX] += this_vol*tet4[ip][4]*xp[0]*xp[0];
            sums[SWYY] += this_vol*tet4[ip][4]*xp[1]*xp[1];
            sums[SWZZ] += this_vol*tet4[ip][4]*xp[2]*xp[2];
            // offd...
            sums[SWYZ] += this_vol*tet4[ip][4]*xp[1]*xp[2];
            sums[SWXZ] += this_vol*tet4[ip][4]*xp[0]*xp[2];
            sums[SWXY] += this_vol*tet4[ip][4]*xp[0]*xp[1];
          }
          for (int ip = 0; ip < 10; ++ip) {
            FOR_J3 xp[j] = (dx[j] + tet10[ip][1]*iter->second[j] + tet10[ip][2]*x0[j] + tet10[ip][3]*x1[j])/r_vv;
            sums[SWXXX] += this_vol*tet10[ip][4]*xp[0]*xp[0]*xp[0];
            sums[SWXXY] += this_vol*tet10[ip][4]*xp[0]*xp[0]*xp[1];
            sums[SWXXZ] += this_vol*tet10[ip][4]*xp[0]*xp[0]*xp[2];
            sums[SWXYY] += this_vol*tet10[ip][4]*xp[0]*xp[1]*xp[1];
            sums[SWXYZ] += this_vol*tet10[ip][4]*xp[0]*xp[1]*xp[2];
            sums[SWXZZ] += this_vol*tet10[ip][4]*xp[0]*xp[2]*xp[2];
            sums[SWYYY] += this_vol*tet10[ip][4]*xp[1]*xp[1]*xp[1];
            sums[SWYYZ] += this_vol*tet10[ip][4]*xp[1]*xp[1]*xp[2];
            sums[SWYZZ] += this_vol*tet10[ip][4]*xp[1]*xp[2]*xp[2];
            sums[SWZZZ] += this_vol*tet10[ip][4]*xp[2]*xp[2]*xp[2];
          }
          for (int ip = 0; ip < 20; ++ip) {
            FOR_J3 xp[j] = (dx[j] + tet20[ip][1]*iter->second[j] + tet20[ip][2]*x0[j] + tet20[ip][3]*x1[j])/r_vv;
            sums[SWXXXX] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[0]*xp[0];
            sums[SWXXXY] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[0]*xp[1];
            sums[SWXXXZ] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[0]*xp[2];
            sums[SWXXYY] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[1]*xp[1];
            sums[SWXXYZ] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[1]*xp[2];
            sums[SWXXZZ] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[2]*xp[2];
            sums[SWXYYY] += this_vol*tet20[ip][4]*xp[0]*xp[1]*xp[1]*xp[1];
            sums[SWXYYZ] += this_vol*tet20[ip][4]*xp[0]*xp[1]*xp[1]*xp[2];
            sums[SWXYZZ] += this_vol*tet20[ip][4]*xp[0]*xp[1]*xp[2]*xp[2];
            sums[SWXZZZ] += this_vol*tet20[ip][4]*xp[0]*xp[2]*xp[2]*xp[2];
            sums[SWYYYY] += this_vol*tet20[ip][4]*xp[1]*xp[1]*xp[1]*xp[1];
            sums[SWYYYZ] += this_vol*tet20[ip][4]*xp[1]*xp[1]*xp[1]*xp[2];
            sums[SWYYZZ] += this_vol*tet20[ip][4]*xp[1]*xp[1]*xp[2]*xp[2];
            sums[SWYZZZ] += this_vol*tet20[ip][4]*xp[1]*xp[2]*xp[2]*xp[2];
            sums[SWZZZZ] += this_vol*tet20[ip][4]*xp[2]*xp[2]*xp[2]*xp[2];
          }
          /*
            vdBfTriVec.push_back(VdTri(iter->second,x0,x1));
            vdBfColorVec.push_back(bf_color[iter0->second]);
            const double this_area = MAG(this_n); // because these cut surface tris can only be cut by voronoi faces, they are convex
            iter->second.area += this_area;
            FOR_J3 {
            iter->second.xc[j]     += this_area*(iter->second[j]+x0[j]+x1[j]);
            x_cv[ip][j]            += this_vol*(iter->second[j]+x0[j]+x1[j]);
            iter->second.normal[j] += this_n[j];
            }
          */
        }
      }

    }
  }

  // =============================================================
  // any orphan chunks?...
  // =============================================================

  int next = vdArray[ip].orphan_chunk_data; // -1 or the index in the ocdVec of our first orphan chunk
  while (next != -1) {

    // edge loop on ocdVec[next]...
    for (int ied = 0; ied < ocdVec[next].getNed(); ++ied) {
      // loop through faces...
      FOR_I2 {
        const int ibf_or_ifa = ocdVec[next].getFaoed(ied,i);
        if (ibf_or_ifa <= -ORPHAN_FACE_OFFSET) {
          // note that we skip open faces just as we skipped group_nbr != 0 faces in the main above...
          assert((-ibf_or_ifa-ORPHAN_FACE_OFFSET >= 0)&&(-ibf_or_ifa-ORPHAN_FACE_OFFSET < ocdVec[next].getNfa()));
          // if this edge shares its face with an open face, then DON'T include...
          if ((ocdVec[next].getFaoed(ied,1-i) <= -1)&&(ocdVec[next].getFaoed(ied,1-i) > -ORPHAN_FACE_OFFSET))
            continue;
          uint8 rbiHash = ocdVec[next].getRankBitsIndexForFace(-ibf_or_ifa-ORPHAN_FACE_OFFSET);
          double * x0 = ocdVec[next].getXnoPtr(ocdVec[next].getNooed(ied,i));
          double * x1 = ocdVec[next].getXnoPtr(ocdVec[next].getNooed(ied,1-i));
          //assert(faMap.find(pair<uint8,double*>(rbiHash,x0)) == faMap.end());
          //faMap[pair<uint8,double*>(rbiHash,x0)] = x1;
          // also, build the geometric stuff...
          map<const uint8,double*>::iterator iter = internalFaceGeometryMap.find(rbiHash);
          if (iter == internalFaceGeometryMap.end()) {
            internalFaceGeometryMap[rbiHash] = x0;
          }
          else if ((x0 != iter->second)&&(x1 != iter->second)) {
            assert(iter->second != NULL);
            const double this_n[3] = TRI_NORMAL_2(iter->second,x0,x1);
            FOR_J3 gcl[j] += this_n[j];
            const double this_vol = DOT_PRODUCT(iter->second,this_n);
            sums[SW] += this_vol;
            // the linear terms sums[1,2,3] can be computed at the mid-point...
            sums[SWX] += this_vol*(dx[0]+0.25*(iter->second[0]+x0[0]+x1[0]))/r_vv;
            sums[SWY] += this_vol*(dx[1]+0.25*(iter->second[1]+x0[1]+x1[1]))/r_vv;
            sums[SWZ] += this_vol*(dx[2]+0.25*(iter->second[2]+x0[2]+x1[2]))/r_vv;
            // the quadratic and above terms require quadrature...
            double xp[3]; 
            for (int ip = 0; ip < 4; ++ip) {
              // note that the first point is the forming point, zero by definition
              FOR_J3 xp[j] = (dx[j] + tet4[ip][1]*iter->second[j] + tet4[ip][2]*x0[j] + tet4[ip][3]*x1[j])/r_vv;
              // tet4[ip][4] is the weight... 
              // diagonal...
              sums[SWXX] += this_vol*tet4[ip][4]*xp[0]*xp[0];
              sums[SWYY] += this_vol*tet4[ip][4]*xp[1]*xp[1];
              sums[SWZZ] += this_vol*tet4[ip][4]*xp[2]*xp[2];
              // offd...
              sums[SWYZ] += this_vol*tet4[ip][4]*xp[1]*xp[2];
              sums[SWXZ] += this_vol*tet4[ip][4]*xp[0]*xp[2];
              sums[SWXY] += this_vol*tet4[ip][4]*xp[0]*xp[1];
            }
            for (int ip = 0; ip < 10; ++ip) {
              FOR_J3 xp[j] = (dx[j] + tet10[ip][1]*iter->second[j] + tet10[ip][2]*x0[j] + tet10[ip][3]*x1[j])/r_vv;
              sums[SWXXX] += this_vol*tet10[ip][4]*xp[0]*xp[0]*xp[0];
              sums[SWXXY] += this_vol*tet10[ip][4]*xp[0]*xp[0]*xp[1];
              sums[SWXXZ] += this_vol*tet10[ip][4]*xp[0]*xp[0]*xp[2];
              sums[SWXYY] += this_vol*tet10[ip][4]*xp[0]*xp[1]*xp[1];
              sums[SWXYZ] += this_vol*tet10[ip][4]*xp[0]*xp[1]*xp[2];
              sums[SWXZZ] += this_vol*tet10[ip][4]*xp[0]*xp[2]*xp[2];
              sums[SWYYY] += this_vol*tet10[ip][4]*xp[1]*xp[1]*xp[1];
              sums[SWYYZ] += this_vol*tet10[ip][4]*xp[1]*xp[1]*xp[2];
              sums[SWYZZ] += this_vol*tet10[ip][4]*xp[1]*xp[2]*xp[2];
              sums[SWZZZ] += this_vol*tet10[ip][4]*xp[2]*xp[2]*xp[2];
            }
            for (int ip = 0; ip < 20; ++ip) {
              FOR_J3 xp[j] = (dx[j] + tet20[ip][1]*iter->second[j] + tet20[ip][2]*x0[j] + tet20[ip][3]*x1[j])/r_vv;
              sums[SWXXXX] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[0]*xp[0];
              sums[SWXXXY] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[0]*xp[1];
              sums[SWXXXZ] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[0]*xp[2];
              sums[SWXXYY] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[1]*xp[1];
              sums[SWXXYZ] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[1]*xp[2];
              sums[SWXXZZ] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[2]*xp[2];
              sums[SWXYYY] += this_vol*tet20[ip][4]*xp[0]*xp[1]*xp[1]*xp[1];
              sums[SWXYYZ] += this_vol*tet20[ip][4]*xp[0]*xp[1]*xp[1]*xp[2];
              sums[SWXYZZ] += this_vol*tet20[ip][4]*xp[0]*xp[1]*xp[2]*xp[2];
              sums[SWXZZZ] += this_vol*tet20[ip][4]*xp[0]*xp[2]*xp[2]*xp[2];
              sums[SWYYYY] += this_vol*tet20[ip][4]*xp[1]*xp[1]*xp[1]*xp[1];
              sums[SWYYYZ] += this_vol*tet20[ip][4]*xp[1]*xp[1]*xp[1]*xp[2];
              sums[SWYYZZ] += this_vol*tet20[ip][4]*xp[1]*xp[1]*xp[2]*xp[2];
              sums[SWYZZZ] += this_vol*tet20[ip][4]*xp[1]*xp[2]*xp[2]*xp[2];
              sums[SWZZZZ] += this_vol*tet20[ip][4]*xp[2]*xp[2]*xp[2]*xp[2];
            }
            /*
              vdFaTriVec.push_back(VdTri(iter->second,x0,x1));
              iter->second.area += this_vol;
              FOR_J3 {
              const double tmp = this_vol*(iter->second[j]+x0[j]+x1[j]);
              iter->second.xc[j] += tmp;
              x_cv[ip][j]        += tmp;
              iter->second.normal[j] += this_n[j];
              }
            */
          }
          else {
            // face loop should be 1-directional...
            assert(iter->second == x1);
          }
        }
        else if (ibf_or_ifa >= 0) {
            
          assert((ibf_or_ifa >= 0)&&(ibf_or_ifa < ocdVec[next].getNbf()));

          // a positive ibf_or_ifa is a surface tri index -- same as above...
          // if this edge shares its face with an open face, then DON'T include...
          if ((ocdVec[next].getFaoed(ied,1-i) <= -1)&&(ocdVec[next].getFaoed(ied,1-i) > -ORPHAN_FACE_OFFSET))
            continue;
            
          int ipart,ist,bits;
          ocdVec[next].getPartStAndBitsOfBf(ipart,ist,bits,ibf_or_ifa);
	    
          double * x0 = ocdVec[next].getXnoPtr(ocdVec[next].getNooed(ied,i));
          double * x1 = ocdVec[next].getXnoPtr(ocdVec[next].getNooed(ied,1-i));
          // and the geometry...
          map<const pair<pair<int,int>,int>,double*>::iterator iter = 
            boundaryFaceGeometryMap.find(pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits)); // only ifa's are guaranteed planar
          if (iter == boundaryFaceGeometryMap.end()) {
            boundaryFaceGeometryMap[pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits)] = x0;
          }
          else if ((x0 != iter->second)&&(x1 != iter->second)) {
            assert(iter->second != NULL);
            const double this_n[3] = TRI_NORMAL_2(iter->second,x0,x1);
            FOR_J3 gcl[j] += this_n[j];
            const double this_vol = DOT_PRODUCT(iter->second,this_n);
            sums[SW] += this_vol;
            // the linear terms sums[1,2,3] can be computed at the mid-point...
            sums[SWX] += this_vol*(dx[0]+0.25*(iter->second[0]+x0[0]+x1[0]))/r_vv;
            sums[SWY] += this_vol*(dx[1]+0.25*(iter->second[1]+x0[1]+x1[1]))/r_vv;
            sums[SWZ] += this_vol*(dx[2]+0.25*(iter->second[2]+x0[2]+x1[2]))/r_vv;
            // the quadratic and above terms require quadrature...
            double xp[3]; 
            for (int ip = 0; ip < 4; ++ip) {
              // note that the first point is the forming point, zero by definition
              FOR_J3 xp[j] = (dx[j] + tet4[ip][1]*iter->second[j] + tet4[ip][2]*x0[j] + tet4[ip][3]*x1[j])/r_vv;
              // tet4[ip][4] is the weight... 
              // diagonal...
              sums[SWXX] += this_vol*tet4[ip][4]*xp[0]*xp[0];
              sums[SWYY] += this_vol*tet4[ip][4]*xp[1]*xp[1];
              sums[SWZZ] += this_vol*tet4[ip][4]*xp[2]*xp[2];
              // offd...
              sums[SWYZ] += this_vol*tet4[ip][4]*xp[1]*xp[2];
              sums[SWXZ] += this_vol*tet4[ip][4]*xp[0]*xp[2];
              sums[SWXY] += this_vol*tet4[ip][4]*xp[0]*xp[1];
            }
            for (int ip = 0; ip < 10; ++ip) {
              FOR_J3 xp[j] = (dx[j] + tet10[ip][1]*iter->second[j] + tet10[ip][2]*x0[j] + tet10[ip][3]*x1[j])/r_vv;
              sums[SWXXX] += this_vol*tet10[ip][4]*xp[0]*xp[0]*xp[0];
              sums[SWXXY] += this_vol*tet10[ip][4]*xp[0]*xp[0]*xp[1];
              sums[SWXXZ] += this_vol*tet10[ip][4]*xp[0]*xp[0]*xp[2];
              sums[SWXYY] += this_vol*tet10[ip][4]*xp[0]*xp[1]*xp[1];
              sums[SWXYZ] += this_vol*tet10[ip][4]*xp[0]*xp[1]*xp[2];
              sums[SWXZZ] += this_vol*tet10[ip][4]*xp[0]*xp[2]*xp[2];
              sums[SWYYY] += this_vol*tet10[ip][4]*xp[1]*xp[1]*xp[1];
              sums[SWYYZ] += this_vol*tet10[ip][4]*xp[1]*xp[1]*xp[2];
              sums[SWYZZ] += this_vol*tet10[ip][4]*xp[1]*xp[2]*xp[2];
              sums[SWZZZ] += this_vol*tet10[ip][4]*xp[2]*xp[2]*xp[2];
            }
            for (int ip = 0; ip < 20; ++ip) {
              FOR_J3 xp[j] = (dx[j] + tet20[ip][1]*iter->second[j] + tet20[ip][2]*x0[j] + tet20[ip][3]*x1[j])/r_vv;
              sums[SWXXXX] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[0]*xp[0];
              sums[SWXXXY] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[0]*xp[1];
              sums[SWXXXZ] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[0]*xp[2];
              sums[SWXXYY] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[1]*xp[1];
              sums[SWXXYZ] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[1]*xp[2];
              sums[SWXXZZ] += this_vol*tet20[ip][4]*xp[0]*xp[0]*xp[2]*xp[2];
              sums[SWXYYY] += this_vol*tet20[ip][4]*xp[0]*xp[1]*xp[1]*xp[1];
              sums[SWXYYZ] += this_vol*tet20[ip][4]*xp[0]*xp[1]*xp[1]*xp[2];
              sums[SWXYZZ] += this_vol*tet20[ip][4]*xp[0]*xp[1]*xp[2]*xp[2];
              sums[SWXZZZ] += this_vol*tet20[ip][4]*xp[0]*xp[2]*xp[2]*xp[2];
              sums[SWYYYY] += this_vol*tet20[ip][4]*xp[1]*xp[1]*xp[1]*xp[1];
              sums[SWYYYZ] += this_vol*tet20[ip][4]*xp[1]*xp[1]*xp[1]*xp[2];
              sums[SWYYZZ] += this_vol*tet20[ip][4]*xp[1]*xp[1]*xp[2]*xp[2];
              sums[SWYZZZ] += this_vol*tet20[ip][4]*xp[1]*xp[2]*xp[2]*xp[2];
              sums[SWZZZZ] += this_vol*tet20[ip][4]*xp[2]*xp[2]*xp[2]*xp[2];
            }
            /*
              vdBfTriVec.push_back(VdTri(iter->second,x0,x1));
              vdBfColorVec.push_back(bf_color[iter0->second]);
              const double this_area = MAG(this_n); // because individual cut surface tris can only be cut by voronoi faces, they are convex
              iter->second.area += this_area;
              FOR_J3 {
              iter->second.xc[j]     += this_area*(iter->second[j]+x0[j]+x1[j]);
              x_cv[ip][j]            += this_vol*(iter->second[j]+x0[j]+x1[j]);
              iter->second.normal[j] += this_n[j];
              }
            */
          }
        }
      }
    }
    next = ocdVec[next].next;
  }
  //if (ip == 2) writeTriColorOFF(vdBfTriVec,vdBfColorVec,0,ip,points->xp[ip]);

  // normalize the x_cv and vol_cv for this ip...
  //vol_cv /= 6.0;

  //cout << "got vol_cv: " << sums[0] << " vol^(1/3): " << pow(sums[0],1.0/3.0) << endl;
  //cout << "got GCL: " << COUT_VEC(gcl) << endl;
    
  FOR_J3 assert(fabs(gcl[j]) < 1.0E-15);
  
  for (map<const uint8,double*>::iterator iter = internalFaceGeometryMap.begin(); iter != internalFaceGeometryMap.end(); ++iter)
    rbiNbrVec.push_back(iter->first);
  
}

int VoronoiBuilder::calcFirstQuadratureTerms(double sums[5],const int ip,const double dx[3],const double r_vv) {
  
  // just computes SW,SWX,SWY,SWZ
  // returns the number of internal nbrs
  
  // coming into this routine, sums[0] contains ip's r_vv. However the passed r_vv is
  // considered separate to allow us to check transformations of these terms...
  
  map<const uint8,double*> internalFaceGeometryMap; // map of face geometry with rbi key
  map<const pair<pair<int,int>,int>,double*> boundaryFaceGeometryMap; // map of face geometry with ((ipart,ist),bits) key
  
  double gcl[3] = { 0.0, 0.0, 0.0 };
  for (int i = 1; i < 5; ++i) sums[i] = 0.0;
    
  // loop through main edges (i.e. no orphans)...
  for (int ied = 0; ied < vdArray[ip].getNedMain(); ++ied) {
    // each edge has 2 faces...
    FOR_I2 {
      const int ibf_or_ifa = vdArray[ip].getFaoed(ied,i);
      // a negative ibf_or_ifa indicates an internal nbr (positive is a boundary)...
      if (ibf_or_ifa < 0) {
        // recall vdArray[ip] uses -1 indexing for faces...
        assert((-ibf_or_ifa-1 >= 0)&&(-ibf_or_ifa-1 < vdArray[ip].getNfa()));
        // confirm that this is a group0 face. It has to be because
        // we are looping on group0 edges...
        // recall that faces have been paired in Step2, and active faces are paired faces where
        // neither face is a "zero" face (i.e. its area is below a certain tolerance)...
        // Recall "group" refers to the orphan indexing, with the main part always
        // being group == 0. Orphans are then group=1 (or group=2, etc if more than 1 orphan associated
        // with the cell), in no particular order.
        if (vdArray[ip].faceIsActive(-ibf_or_ifa-1)) {
          int group,rank,bits,index;
          vdArray[ip].getGrbiForFace(group,rank,bits,index,-ibf_or_ifa-1);
          assert(group == 0); // must be group 0 -- this is a check.
          int group_nbr,ifa_nbr;
          vdArray[ip].getGiForFaceNbr(group_nbr,ifa_nbr,-ibf_or_ifa-1);
          if (group_nbr == 0) {
            // this is a boundary face that should probably be included, unless it shares
            // its edge with a non-group 0 face nbr...
            if (vdArray[ip].getFaoed(ied,1-i) < 0) {
              const int ifa_opp = -vdArray[ip].getFaoed(ied,1-i)-1;
              assert((ifa_opp >= 0)&&(ifa_opp < vdArray[ip].getNfa()));
              if (vdArray[ip].faceIsActive(ifa_opp)) {
                int group_nbr,ifa_nbr;
                vdArray[ip].getGiForFaceNbr(group_nbr,ifa_nbr,ifa_opp);
                if (group_nbr != 0)
                  continue;
              }
            }
            // this is an active face edge - add it to the faceMap...
            const uint8 rbiHash = BitUtils::packRankBitsIndex(rank,bits,index); // rbi uniquely represents the nbr we have been cut against
            // recall i = 0,1, so edge goes form node 0->1 or 1->0, depending on the face
            double * x0 = vdArray[ip].getXnoPtr(vdArray[ip].getNooed(ied,i)); // the coordinates of this node would be x0[0],x0[1],x0[2].
            double * x1 = vdArray[ip].getXnoPtr(vdArray[ip].getNooed(ied,1-i));
            //assert(faMap.find(pair<uint8,double*>(rbiHash,x0)) == faMap.end());
            //faMap[pair<uint8,double*>(rbiHash,x0)] = x1;
            // also, build the geometric stuff...
            map<const uint8,double*>::iterator iter = internalFaceGeometryMap.find(rbiHash);
            if (iter == internalFaceGeometryMap.end()) {
              internalFaceGeometryMap[rbiHash] = x0;
            }
            else if ((x0 != iter->second)&&(x1 != iter->second)) {
              assert(iter->second != NULL);
              const double this_n[3] = TRI_NORMAL_2(iter->second,x0,x1);
              FOR_J3 gcl[j] += this_n[j];
              const double this_vol = DOT_PRODUCT(iter->second,this_n);
              sums[SW] += this_vol;
              // the linear terms can be computed at the mid-point...
              sums[SWX] += this_vol*(dx[0]+0.25*(iter->second[0]+x0[0]+x1[0]))/r_vv;
              sums[SWY] += this_vol*(dx[1]+0.25*(iter->second[1]+x0[1]+x1[1]))/r_vv;
              sums[SWZ] += this_vol*(dx[2]+0.25*(iter->second[2]+x0[2]+x1[2]))/r_vv;
            }
            else {
              // face loop should be 1-directional...
              assert(iter->second == x1);
            }
          }
        }
      }
      else {
        // a positive ibf_or_ifa is the ibf index, meaning it is an edge of a boundary face.
        // If the other face is an internal face, and
        // it has a connection with an orphan (i.e. group > 0), we skip it...
        assert((ibf_or_ifa >= 0)&&(ibf_or_ifa < vdArray[ip].getNbf()));
        int ipart,ist,bits;
        vdArray[ip].getPartStAndBitsOfBf(ipart,ist,bits,ibf_or_ifa);
        //map<const pair<pair<int,int>,int>,int>::iterator iter0 = colorMap.end();
        if (vdArray[ip].getFaoed(ied,1-i) < 0) {
          const int ifa_opp = -vdArray[ip].getFaoed(ied,1-i)-1;
          assert((ifa_opp >= 0)&&(ifa_opp < vdArray[ip].getNfa()));
          if (vdArray[ip].faceIsActive(ifa_opp)) {
            int group_nbr,ifa_nbr;
            vdArray[ip].getGiForFaceNbr(group_nbr,ifa_nbr,ifa_opp);
            if (group_nbr != 0) {
              continue;
            }
          }
        }
        // recall i = 0,1, so edge goes form node 0->1 or 1->0, depending on the face
        double * x0 = vdArray[ip].getXnoPtr(vdArray[ip].getNooed(ied,i)); // the coordinates of this node would be x0[0],x0[1],x0[2].
        double * x1 = vdArray[ip].getXnoPtr(vdArray[ip].getNooed(ied,1-i));
        // geometry...
        map<const pair<pair<int,int>,int>,double*>::iterator iter =
          boundaryFaceGeometryMap.find(pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits));
        if (iter == boundaryFaceGeometryMap.end()) {
          boundaryFaceGeometryMap[pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits)] = x0;
        }
        else if ((x0 != iter->second)&&(x1 != iter->second)) {
          assert(iter->second != NULL);
          const double this_n[3] = TRI_NORMAL_2(iter->second,x0,x1);
          FOR_J3 gcl[j] += this_n[j];
          const double this_vol = DOT_PRODUCT(iter->second,this_n);
          sums[SW] += this_vol;
          // the linear terms sums[1,2,3] can be computed at the mid-point...
          sums[SWX] += this_vol*(dx[0]+0.25*(iter->second[0]+x0[0]+x1[0]))/r_vv;
          sums[SWY] += this_vol*(dx[1]+0.25*(iter->second[1]+x0[1]+x1[1]))/r_vv;
          sums[SWZ] += this_vol*(dx[2]+0.25*(iter->second[2]+x0[2]+x1[2]))/r_vv;
        }
      }

    }
  }

  // =============================================================
  // any orphan chunks?...
  // =============================================================

  int next = vdArray[ip].orphan_chunk_data; // -1 or the index in the ocdVec of our first orphan chunk
  while (next != -1) {
    // edge loop on ocdVec[next]...
    for (int ied = 0; ied < ocdVec[next].getNed(); ++ied) {
      // loop through faces...
      FOR_I2 {
        const int ibf_or_ifa = ocdVec[next].getFaoed(ied,i);
        if (ibf_or_ifa <= -ORPHAN_FACE_OFFSET) {
          // note that we skip open faces just as we skipped group_nbr != 0 faces in the main above...
          assert((-ibf_or_ifa-ORPHAN_FACE_OFFSET >= 0)&&(-ibf_or_ifa-ORPHAN_FACE_OFFSET < ocdVec[next].getNfa()));
          // if this edge shares its face with an open face, then DON'T include...
          if ((ocdVec[next].getFaoed(ied,1-i) <= -1)&&(ocdVec[next].getFaoed(ied,1-i) > -ORPHAN_FACE_OFFSET))
            continue;
          uint8 rbiHash = ocdVec[next].getRankBitsIndexForFace(-ibf_or_ifa-ORPHAN_FACE_OFFSET);
          double * x0 = ocdVec[next].getXnoPtr(ocdVec[next].getNooed(ied,i));
          double * x1 = ocdVec[next].getXnoPtr(ocdVec[next].getNooed(ied,1-i));
          //assert(faMap.find(pair<uint8,double*>(rbiHash,x0)) == faMap.end());
          //faMap[pair<uint8,double*>(rbiHash,x0)] = x1;
          // also, build the geometric stuff...
          map<const uint8,double*>::iterator iter = internalFaceGeometryMap.find(rbiHash);
          if (iter == internalFaceGeometryMap.end()) {
            internalFaceGeometryMap[rbiHash] = x0;
          }
          else if ((x0 != iter->second)&&(x1 != iter->second)) {
            assert(iter->second != NULL);
            const double this_n[3] = TRI_NORMAL_2(iter->second,x0,x1);
            FOR_J3 gcl[j] += this_n[j];
            const double this_vol = DOT_PRODUCT(iter->second,this_n);
            sums[SW] += this_vol;
            // the linear terms sums[1,2,3] can be computed at the mid-point...
            sums[SWX] += this_vol*(dx[0]+0.25*(iter->second[0]+x0[0]+x1[0]))/r_vv;
            sums[SWY] += this_vol*(dx[1]+0.25*(iter->second[1]+x0[1]+x1[1]))/r_vv;
            sums[SWZ] += this_vol*(dx[2]+0.25*(iter->second[2]+x0[2]+x1[2]))/r_vv;
          }
          else {
            // face loop should be 1-directional...
            assert(iter->second == x1);
          }
        }
        else if (ibf_or_ifa >= 0) {
          assert((ibf_or_ifa >= 0)&&(ibf_or_ifa < ocdVec[next].getNbf()));
          // a positive ibf_or_ifa is a surface tri index -- same as above...
          // if this edge shares its face with an open face, then DON'T include...
          if ((ocdVec[next].getFaoed(ied,1-i) <= -1)&&(ocdVec[next].getFaoed(ied,1-i) > -ORPHAN_FACE_OFFSET))
            continue;
          int ipart,ist,bits;
          ocdVec[next].getPartStAndBitsOfBf(ipart,ist,bits,ibf_or_ifa);
          double * x0 = ocdVec[next].getXnoPtr(ocdVec[next].getNooed(ied,i));
          double * x1 = ocdVec[next].getXnoPtr(ocdVec[next].getNooed(ied,1-i));
          // and the geometry...
          map<const pair<pair<int,int>,int>,double*>::iterator iter = 
            boundaryFaceGeometryMap.find(pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits)); // only ifa's are guaranteed planar
          if (iter == boundaryFaceGeometryMap.end()) {
            boundaryFaceGeometryMap[pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits)] = x0;
          }
          else if ((x0 != iter->second)&&(x1 != iter->second)) {
            assert(iter->second != NULL);
            const double this_n[3] = TRI_NORMAL_2(iter->second,x0,x1);
            FOR_J3 gcl[j] += this_n[j];
            const double this_vol = DOT_PRODUCT(iter->second,this_n);
            sums[SW] += this_vol;
            // the linear terms sums[1,2,3] can be computed at the mid-point...
            sums[SWX] += this_vol*(dx[0]+0.25*(iter->second[0]+x0[0]+x1[0]))/r_vv;
            sums[SWY] += this_vol*(dx[1]+0.25*(iter->second[1]+x0[1]+x1[1]))/r_vv;
            sums[SWZ] += this_vol*(dx[2]+0.25*(iter->second[2]+x0[2]+x1[2]))/r_vv;
          }
        }
      }
    }
    next = ocdVec[next].next;
  }
  FOR_J3 assert(fabs(gcl[j]) < 1.0E-15);
  
  return internalFaceGeometryMap.size();
  
}

void VoronoiBuilder::buildTriVecsAndNbrs(vector<VdIndexAndTri>& vdFaTriVec,vector<VdTri>& vdBfTriVec,vector<uint8>& rbiNbrVec,const int ip) {
  
  // and build the tri vec -- basically the tris that, when integrated with the
  // center point -- which is zero -- could be used to compute the volume of
  // the Voronoi diagram.
  
  assert(vdFaTriVec.empty());
  assert(vdBfTriVec.empty());
  //assert(vdBfColorVec.empty());

  map<const uint8,double*> internalFaceGeometryMap; // map of face geometry with rbi key
  map<const pair<pair<int,int>,int>,double*> boundaryFaceGeometryMap; // map of face geometry with ((ipart,ist),bits) key
    
  //double gcl[3] = { 0.0, 0.0, 0.0 };
    
  // loop through main edges (i.e. no orphans)...
  for (int ied = 0; ied < vdArray[ip].getNedMain(); ++ied) {
    // each edge has 2 faces...
    FOR_I2 {
      const int ibf_or_ifa = vdArray[ip].getFaoed(ied,i);
      // a negative ibf_or_ifa indicates an internal nbr (positive is a boundary)...
      if (ibf_or_ifa < 0) {
        // recall vdArray[ip] uses -1 indexing for faces...
        assert((-ibf_or_ifa-1 >= 0)&&(-ibf_or_ifa-1 < vdArray[ip].getNfa()));
        // confirm that this is a group0 face. It has to be because
        // we are looping on group0 edges...
        // recall that faces have been paired in Step2, and active faces are paired faces where
        // neither face is a "zero" face (i.e. its area is below a certain tolerance)...
        // Recall "group" refers to the orphan indexing, with the main part always
        // being group == 0. Orphans are then group=1 (or group=2, etc if more than 1 orphan associated
        // with the cell), in no particular order.
        if (vdArray[ip].faceIsActive(-ibf_or_ifa-1)) {
          int group,rank,bits,index;
          vdArray[ip].getGrbiForFace(group,rank,bits,index,-ibf_or_ifa-1);
          assert(group == 0); // must be group 0 -- this is a check.
          int group_nbr,ifa_nbr;
          vdArray[ip].getGiForFaceNbr(group_nbr,ifa_nbr,-ibf_or_ifa-1);
          if (group_nbr == 0) {
            // this is a boundary face that should probably be included, unless it shares
            // its edge with a non-group 0 face nbr...
            if (vdArray[ip].getFaoed(ied,1-i) < 0) {
              const int ifa_opp = -vdArray[ip].getFaoed(ied,1-i)-1;
              assert((ifa_opp >= 0)&&(ifa_opp < vdArray[ip].getNfa()));
              if (vdArray[ip].faceIsActive(ifa_opp)) {
                int group_nbr,ifa_nbr;
                vdArray[ip].getGiForFaceNbr(group_nbr,ifa_nbr,ifa_opp);
                if (group_nbr != 0)
                  continue;
              }
            }
            // this is an active face edge - add it to the faceMap...
            const uint8 rbiHash = BitUtils::packRankBitsIndex(rank,bits,index); // rbi uniquely represents the nbr we have been cut against
            // recall i = 0,1, so edge goes form node 0->1 or 1->0, depending on the face
            double * x0 = vdArray[ip].getXnoPtr(vdArray[ip].getNooed(ied,i)); // the coordinates of this node would be x0[0],x0[1],x0[2].
            double * x1 = vdArray[ip].getXnoPtr(vdArray[ip].getNooed(ied,1-i));
            //assert(faMap.find(pair<uint8,double*>(rbiHash,x0)) == faMap.end());
            //faMap[pair<uint8,double*>(rbiHash,x0)] = x1;
            // also, build the geometric stuff...
            map<const uint8,double*>::iterator iter = internalFaceGeometryMap.find(rbiHash);
            if (iter == internalFaceGeometryMap.end()) {
              internalFaceGeometryMap[rbiHash] = x0;
            }
            else if ((x0 != iter->second)&&(x1 != iter->second)) {
              assert(iter->second != NULL);
              vdFaTriVec.push_back(VdIndexAndTri(rbiHash,iter->second,x0,x1));
            }
            else {
              // face loop should be 1-directional...
              assert(iter->second == x1);
            }
          }
        }
      }
      else {
        // a positive ibf_or_ifa is the ibf index, meaning it is an edge of a boundary face.
        // If the other face is an internal face, and
        // it has a connection with an orphan (i.e. group > 0), we skip it...
        assert((ibf_or_ifa >= 0)&&(ibf_or_ifa < vdArray[ip].getNbf()));
        int ipart,ist,bits;
        vdArray[ip].getPartStAndBitsOfBf(ipart,ist,bits,ibf_or_ifa);
        //map<const pair<pair<int,int>,int>,int>::iterator iter0 = colorMap.end();
        if (vdArray[ip].getFaoed(ied,1-i) < 0) {
          const int ifa_opp = -vdArray[ip].getFaoed(ied,1-i)-1;
          assert((ifa_opp >= 0)&&(ifa_opp < vdArray[ip].getNfa()));
          if (vdArray[ip].faceIsActive(ifa_opp)) {
            int group_nbr,ifa_nbr;
            vdArray[ip].getGiForFaceNbr(group_nbr,ifa_nbr,ifa_opp);
            if (group_nbr != 0) {
              continue;
            }
          }
        }
        // recall i = 0,1, so edge goes form node 0->1 or 1->0, depending on the face
        double * x0 = vdArray[ip].getXnoPtr(vdArray[ip].getNooed(ied,i)); // the coordinates of this node would be x0[0],x0[1],x0[2].
        double * x1 = vdArray[ip].getXnoPtr(vdArray[ip].getNooed(ied,1-i));
        // geometry...
        map<const pair<pair<int,int>,int>,double*>::iterator iter =
          boundaryFaceGeometryMap.find(pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits));
        if (iter == boundaryFaceGeometryMap.end()) {
          boundaryFaceGeometryMap[pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits)] = x0;
        }
        else if ((x0 != iter->second)&&(x1 != iter->second)) {
          assert(iter->second != NULL);
          vdBfTriVec.push_back(VdTri(iter->second,x0,x1)); // TODO: eventually we need to put the above bf identifier (and coloring?) into a uint8 -- why not use this from the start?
        }
      }
      
    }
  }

  // =============================================================
  // any orphan chunks?...
  // =============================================================

  int next = vdArray[ip].orphan_chunk_data; // -1 or the index in the ocdVec of our first orphan chunk
  while (next != -1) {
    // edge loop on ocdVec[next]...
    for (int ied = 0; ied < ocdVec[next].getNed(); ++ied) {
      // loop through faces...
      FOR_I2 {
        const int ibf_or_ifa = ocdVec[next].getFaoed(ied,i);
        if (ibf_or_ifa <= -ORPHAN_FACE_OFFSET) {
          // note that we skip open faces just as we skipped group_nbr != 0 faces in the main above...
          assert((-ibf_or_ifa-ORPHAN_FACE_OFFSET >= 0)&&(-ibf_or_ifa-ORPHAN_FACE_OFFSET < ocdVec[next].getNfa()));
          // if this edge shares its face with an open face, then DON'T include...
          if ((ocdVec[next].getFaoed(ied,1-i) <= -1)&&(ocdVec[next].getFaoed(ied,1-i) > -ORPHAN_FACE_OFFSET))
            continue;
          uint8 rbiHash = ocdVec[next].getRankBitsIndexForFace(-ibf_or_ifa-ORPHAN_FACE_OFFSET);
          double * x0 = ocdVec[next].getXnoPtr(ocdVec[next].getNooed(ied,i));
          double * x1 = ocdVec[next].getXnoPtr(ocdVec[next].getNooed(ied,1-i));
          //assert(faMap.find(pair<uint8,double*>(rbiHash,x0)) == faMap.end());
          //faMap[pair<uint8,double*>(rbiHash,x0)] = x1;
          // also, build the geometric stuff...
          map<const uint8,double*>::iterator iter = internalFaceGeometryMap.find(rbiHash);
          if (iter == internalFaceGeometryMap.end()) {
            internalFaceGeometryMap[rbiHash] = x0;
          }
          else if ((x0 != iter->second)&&(x1 != iter->second)) {
            assert(iter->second != NULL);
            vdFaTriVec.push_back(VdIndexAndTri(rbiHash,iter->second,x0,x1));
          }
          else {
            // face loop should be 1-directional...
            assert(iter->second == x1);
          }
        }
        else if (ibf_or_ifa >= 0) {
          assert((ibf_or_ifa >= 0)&&(ibf_or_ifa < ocdVec[next].getNbf()));
          // a positive ibf_or_ifa is a surface tri index -- same as above...
          // if this edge shares its face with an open face, then DON'T include...
          if ((ocdVec[next].getFaoed(ied,1-i) <= -1)&&(ocdVec[next].getFaoed(ied,1-i) > -ORPHAN_FACE_OFFSET))
            continue;
          int ipart,ist,bits;
          ocdVec[next].getPartStAndBitsOfBf(ipart,ist,bits,ibf_or_ifa);
          double * x0 = ocdVec[next].getXnoPtr(ocdVec[next].getNooed(ied,i));
          double * x1 = ocdVec[next].getXnoPtr(ocdVec[next].getNooed(ied,1-i));
          // and the geometry...
          map<const pair<pair<int,int>,int>,double*>::iterator iter = 
            boundaryFaceGeometryMap.find(pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits)); // only ifa's are guaranteed planar
          if (iter == boundaryFaceGeometryMap.end()) {
            boundaryFaceGeometryMap[pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits)] = x0;
          }
          else if ((x0 != iter->second)&&(x1 != iter->second)) {
            assert(iter->second != NULL);
            vdBfTriVec.push_back(VdTri(iter->second,x0,x1)); // TODO: eventually we need to put the above bf identifier (and coloring?) into a uint8 -- why not use this from the start?
          }
        }
      }
    }
    next = ocdVec[next].next;
  }
  
  sort(vdFaTriVec.begin(),vdFaTriVec.end());
       
  for (map<const uint8,double*>::iterator iter = internalFaceGeometryMap.begin(); iter != internalFaceGeometryMap.end(); ++iter)
    rbiNbrVec.push_back(iter->first);
  
}

#undef SR // holds the length scale that normalizes all x's

#undef SW // volume (*6?)

#undef SWX 
#undef SWY 
#undef SWZ 

#undef SWXX 
#undef SWXY 
#undef SWXZ 
#undef SWYY 
#undef SWYZ 
#undef SWZZ 

#undef SWXXX 
#undef SWXXY 
#undef SWXXZ 
#undef SWXYY 
#undef SWXYZ 
#undef SWXZZ 
#undef SWYYY 
#undef SWYYZ 
#undef SWYZZ 
#undef SWZZZ 

#undef SWXXXX 
#undef SWXXXY 
#undef SWXXXZ 
#undef SWXXYY 
#undef SWXXYZ 
#undef SWXXZZ 
#undef SWXYYY 
#undef SWXYYZ 
#undef SWXYZZ 
#undef SWXZZZ 
#undef SWYYYY 
#undef SWYYYZ 
#undef SWYYZZ 
#undef SWYZZZ 
#undef SWZZZZ 

