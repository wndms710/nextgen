#include "SimpleSurface.hpp"

// generates 4-digit naca airfoils

int SimpleSurface::addNaca00Airfoil(const int xx,const int _n_panels,const double chord,const int _n_span,const double span,const bool b_sharp) {


  // make sure enough points to discretize profile
  int n_panels = _n_panels;
  if (n_panels < 10) {
    CWARN("using insufficient umber of panels to properly characterize airfoil; setting N = 10");
    n_panels = 10;
  }

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  // compute number of new nodes, tris
  const double d_chord = chord/n_panels;

  int n_spanwise_panels = int(ceil(span/d_chord));
  assert(n_spanwise_panels>0);

  if (_n_span > 0) n_spanwise_panels = _n_span;  // user can override default spanwise panelling

  const double d_span = span/n_spanwise_panels;
  const double thickness = double(xx)*.01*chord;

  // because these are symmetric airfoils the leading edge and trailing edge nodes are the same.
  // total number of nodes in a single profile is: 2*n_panels
  // for all profiles multiply by (n_spanwise_panels+1)
  nsp += 2*n_panels*(n_spanwise_panels+1);
  growNspData(nsp,ss_nsp0);

  // assume initial reference point is 0,0,0; move in +X and +Y direction for chord/span
  const double x_ref[3] = {0.0,0.0,0.0};

  // initialize all points to leading edge at first profile; simply update span and chord as necessary
  for (int isp=ss_nsp0; isp<nsp; ++isp) {
    FOR_I3 xsp[isp][i] = x_ref[i];
  }

  double scaling;
  double coeffs[5];
  if (!b_sharp) {
    scaling = 1.0;
    coeffs[0] = 0.2969;
    coeffs[1] = -0.126;
    coeffs[2] = -0.3516;
    coeffs[3] = 0.2843;
    coeffs[4] = -0.1015;
  }
  else {
    // corresponds to the sharp & scaled trailing edge defined by NASAs Turb model validation fro NACA0012
    scaling = 1.0/1.008930411365;
    coeffs[0] = 0.298222773;
    coeffs[1] = -0.127125232;
    coeffs[2] = -0.357907906;
    coeffs[3] = 0.291984971;
    coeffs[4] = -0.105174606;
  }

  for (int i_span=0, n_span = (n_spanwise_panels+1); i_span < n_span; ++i_span) {
    // do leading and trailing edge points first
    const int le_offset = ss_nsp0 + i_span*2*n_panels;  // leading edge node offset
    const int le_next_offset = ss_nsp0 + (i_span+1)*2*n_panels;  // leading edge node offset of next profile

    // update spanwise position of everybody on this profile
    for (int isp=le_offset; isp<le_next_offset; ++isp) xsp[isp][1] += i_span*d_span;

    // update trailing edge node
    xsp[le_offset + n_panels][0] += chord;

    // points between endpoints
    for (int ic=1, nc=n_panels; ic<nc; ++ic) {
      // uniform sampling
      const double x_uniform = double(ic)/double(n_panels);  // percentage in [0,1]
      const double x = 0.5*(-1.0*cos(M_PI*x_uniform) + 1);  // smaller sampling near endpoints

      // update both top and bottom
      const double x_ic = x*chord;
      xsp[le_offset+ic][0] += x_ic;
      xsp[le_next_offset-ic][0] += x_ic;

      double dz = coeffs[0]*sqrt(x) + coeffs[1]*x + coeffs[2]*x*x + coeffs[3]*x*x*x + coeffs[4]*x*x*x*x;
      dz *= 5*thickness*scaling;  // thickness includes scaling by chord

      xsp[le_offset+ic][2] += dz;
      xsp[le_next_offset-ic][2] -= dz;
    }
  }
  COUT2(" > computed location for " << (nsp-ss_nsp0) << " new nodes");

  const int zone_top = zoneVec.size();
  zoneVec.push_back(SurfaceZone("blade_top"));
  const int zone_bottom = zoneVec.size();
  zoneVec.push_back(SurfaceZone("blade_bottom"));
  nsz += 2;

  nst += 2*2*n_panels*n_spanwise_panels;  // number of new tris
  growNstData(nst,ss_nst0);

  int count = 0;
  for (int i_span=0; i_span<n_spanwise_panels; ++i_span) {
    const int le_offset = ss_nsp0 + i_span*2*n_panels;  // leading edge node offset
    const int le_next_offset = ss_nsp0 + (i_span+1)*2*n_panels;  // leading edge node offset of next profile

    // first add top panels
    for (int i_chord=0; i_chord<n_panels; ++i_chord) {
      const int ist0 = ss_nst0 + count;
      const int ist1 = ss_nst0 + count + 1;

      if (ist0 >= nst || ist0 < 0) COUT2("whoah ist0: " << ist0 << " i_chord: " << i_chord);
      if (ist1 >= nst || ist1 < 0) COUT2("whoah ist1: " << ist1 << " i_chord: " << i_chord);

      spost[ist0][0] = le_offset + i_chord + 0;
      spost[ist0][1] = le_offset + i_chord + 1;
      spost[ist0][2] = le_next_offset + i_chord + 0;
      znost[ist0] = zone_top;

      spost[ist1][0] = le_offset + i_chord + 1;
      spost[ist1][1] = le_next_offset + i_chord + 1;
      spost[ist1][2] = le_next_offset + i_chord + 0;
      znost[ist1] = zone_top;
      szost[ist1] = nsz-2;

      count +=2;
    }

    // do bottom tris, peeling last loop becausenodes come around
    for (int i_chord=0,end=(n_panels-1); i_chord<end; ++i_chord) {
      const int ist0 = ss_nst0 + count + 0;
      const int ist1 = ss_nst0 + count + 1;
      spost[ist0][0] = le_offset + n_panels + i_chord + 0;
      spost[ist0][1] = le_offset + n_panels + i_chord + 1;
      spost[ist0][2] = le_next_offset + n_panels + i_chord + 0;
      znost[ist0] = zone_bottom;
      szost[ist0] = nsz-1;

      spost[ist1][0] = le_offset + n_panels + i_chord + 1;
      spost[ist1][1] = le_next_offset + n_panels + i_chord + 1;
      spost[ist1][2] = le_next_offset + n_panels + i_chord + 0;
      znost[ist1] = zone_bottom;
      szost[ist1] = nsz-1;

      count+=2;
    }

    // last two tris
    const int ist0 = ss_nst0 + count + 0;
    const int ist1 = ss_nst0 + count + 1;
    spost[ist0][0] = le_offset + 0;
    spost[ist0][1] = le_next_offset + 0;
    spost[ist0][2] = le_next_offset + 2*n_panels - 1;
    znost[ist0] = zone_bottom;
    szost[ist0] = nsz-1;

    spost[ist1][0] = le_offset + 0;
    spost[ist1][1] = le_next_offset + 2*n_panels - 1;
    spost[ist1][2] = le_next_offset - 1;
    znost[ist1] = zone_bottom;
    szost[ist1] = nsz-1;
    count += 2;
  }
  assert(count == nst-ss_nst0);
  COUT2(" > added " << (nsp-ss_nsp0) << " new tris");

  // ====================================================
  // make sure all zones were set...
  for (int ist = 0; ist < nst; ++ist) {
    assert((znost[ist] >= 0)&&(znost[ist] < int(zoneVec.size())));
  }

  return 0;
}

int SimpleSurface::addNaca4DigitAirfoil(const int _m, const int _p,const int xx,const int _n_panels,const double chord,const int _n_span,const double span) {


  // make sure enough points to discretize profile
  int n_panels = _n_panels;
  if (n_panels < 10) {
    CWARN("using insufficient umber of panels to properly characterize airfoil; setting N = 10");
    n_panels = 10;
  }

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  // compute number of new nodes, tris
  const double d_chord = chord/n_panels;

  int n_spanwise_panels = int(ceil(span/d_chord));
  assert(n_spanwise_panels>0);

  if (_n_span > 0) n_spanwise_panels = _n_span;  // user can override default spanwise panelling

  const double d_span = span/n_spanwise_panels;
  const double thickness = double(xx)*.01;
  const double m = double(_m)*.01;
  const double p = double(_p)*.1;

  COUT2(" > max-thickness (chord normalized): " << thickness);
  COUT2(" > max-camber: " << m);
  COUT2(" > max-camber loction (chord normalized): " << p);

  // because these are symmetric airfoils the leading edge and trailing edge nodes are the same.
  // total number of nodes in a single profile is: 2*n_panels
  // for all profiles multiply by (n_spanwise_panels+1)
  nsp += 2*n_panels*(n_spanwise_panels+1);
  growNspData(nsp,ss_nsp0);

  // assume initial reference point is 0,0,0; move in +X and +Y direction for chord/span
  const double x_ref[3] = {0.0,0.0,0.0};

  // initialize all points to leading edge at first profile; simply update span and chord as necessary
  for (int isp=ss_nsp0; isp<nsp; ++isp) {
    FOR_I3 xsp[isp][i] = x_ref[i];
  }

  for (int i_span=0, n_span = (n_spanwise_panels+1); i_span < n_span; ++i_span) {
    // do leading and trailing edge points first
    const int le_offset = ss_nsp0 + i_span*2*n_panels;  // leading edge node offset
    const int le_next_offset = ss_nsp0 + (i_span+1)*2*n_panels;  // leading edge node offset of next profile

    // update spanwise position of everybody on this profile
    for (int isp=le_offset; isp<le_next_offset; ++isp) xsp[isp][1] += i_span*d_span;

    // update trailing edge node
    xsp[le_offset + n_panels][0] += chord;

    // points between endpoints
    for (int ic=1, nc=n_panels; ic<nc; ++ic) {
      // uniform sampling
      const double x_uniform = double(ic)/double(n_panels);  // percentage in [0,1]
      const double x = 0.5*(-1.0*cos(M_PI*x_uniform) + 1);  // smaller sampling near endpoints

      double dz = 0.2969*sqrt(x) - 0.126*x - 0.3516*x*x + 0.2843*x*x*x - 0.1015*x*x*x*x;  // symmetric component about camber
      double dz_c = m*(2*p*x - x*x);
      double grad_dz_c = 2*m*(p-x);
      if (x <= p) {
        dz_c /= (p*p);
        grad_dz_c /= (p*p);
      }
      else {
        dz_c += m*(1-2*p);
        dz_c /= (1-p)*(1-p);
        grad_dz_c /= (1-p)*(1-p);
      }
      const double theta = atan(grad_dz_c);
      dz *= 5*thickness*chord;  // thickness does not include chord
      dz_c *= chord;

      // update both top and bottom nodes
      // chord-dimension
      xsp[le_offset+ic][0] += chord*x - dz*sin(theta);
      xsp[le_next_offset-ic][0] += chord*x + dz*sin(theta);

      // chord normal position
      xsp[le_offset+ic][2] += dz_c + dz*cos(theta);
      xsp[le_next_offset-ic][2] += dz_c - dz*cos(theta);
    }
  }
  COUT2(" > computed location for " << (nsp-ss_nsp0) << " new nodes");

  const int zone_top = zoneVec.size();
  zoneVec.push_back(SurfaceZone("blade_top"));
  const int zone_bottom = zoneVec.size();
  zoneVec.push_back(SurfaceZone("blade_bottom"));
  nsz += 2;

  nst += 2*2*n_panels*n_spanwise_panels;  // number of new tris
  growNstData(nst,ss_nst0);

  int count = 0;
  for (int i_span=0; i_span<n_spanwise_panels; ++i_span) {
    const int le_offset = ss_nsp0 + i_span*2*n_panels;  // leading edge node offset
    const int le_next_offset = ss_nsp0 + (i_span+1)*2*n_panels;  // leading edge node offset of next profile

    // first add top panels
    for (int i_chord=0; i_chord<n_panels; ++i_chord) {
      const int ist0 = ss_nst0 + count;
      const int ist1 = ss_nst0 + count + 1;

      if (ist0 >= nst || ist0 < 0) COUT2("whoah ist0: " << ist0 << " i_chord: " << i_chord);
      if (ist1 >= nst || ist1 < 0) COUT2("whoah ist1: " << ist1 << " i_chord: " << i_chord);

      spost[ist0][0] = le_offset + i_chord + 0;
      spost[ist0][1] = le_offset + i_chord + 1;
      spost[ist0][2] = le_next_offset + i_chord + 0;
      znost[ist0] = zone_top;

      spost[ist1][0] = le_offset + i_chord + 1;
      spost[ist1][1] = le_next_offset + i_chord + 1;
      spost[ist1][2] = le_next_offset + i_chord + 0;
      znost[ist1] = zone_top;
      szost[ist1] = nsz-2;

      count +=2;
    }

    // do bottom tris, peeling last loop becausenodes come around
    for (int i_chord=0,end=(n_panels-1); i_chord<end; ++i_chord) {
      const int ist0 = ss_nst0 + count + 0;
      const int ist1 = ss_nst0 + count + 1;
      spost[ist0][0] = le_offset + n_panels + i_chord + 0;
      spost[ist0][1] = le_offset + n_panels + i_chord + 1;
      spost[ist0][2] = le_next_offset + n_panels + i_chord + 0;
      znost[ist0] = zone_bottom;
      szost[ist0] = nsz-1;

      spost[ist1][0] = le_offset + n_panels + i_chord + 1;
      spost[ist1][1] = le_next_offset + n_panels + i_chord + 1;
      spost[ist1][2] = le_next_offset + n_panels + i_chord + 0;
      znost[ist1] = zone_bottom;
      szost[ist1] = nsz-1;

      count+=2;
    }

    // last two tris
    const int ist0 = ss_nst0 + count + 0;
    const int ist1 = ss_nst0 + count + 1;
    spost[ist0][0] = le_offset + 0;
    spost[ist0][1] = le_next_offset + 0;
    spost[ist0][2] = le_next_offset + 2*n_panels - 1;
    znost[ist0] = zone_bottom;
    szost[ist0] = nsz-1;

    spost[ist1][0] = le_offset + 0;
    spost[ist1][1] = le_next_offset + 2*n_panels - 1;
    spost[ist1][2] = le_next_offset - 1;
    znost[ist1] = zone_bottom;
    szost[ist1] = nsz-1;
    count += 2;
  }
  assert(count == nst-ss_nst0);
  COUT2(" > added " << (nsp-ss_nsp0) << " new tris");

  // ====================================================
  // make sure all zones were set...
  for (int ist = 0; ist < nst; ++ist) {
    assert((znost[ist] >= 0)&&(znost[ist] < int(zoneVec.size())));
  }

  return 0;
}
