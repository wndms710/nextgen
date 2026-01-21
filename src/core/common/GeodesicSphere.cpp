#include "GeodesicSphere.hpp"

namespace GeodesicSphere {
  
  // functions realted to the icosahedral sphere...
  
  int getSphereNodeCount(const int n) {
    
    // n is the number of segments each edge is divied into... 
  
    assert(n>=1);
    return 12 + 30*(n-1) + 20*(n-1)*(n-2)/2;
    
  }
  
  int getSphereTriCount(const int n) {

    // n is the number of segments each edge is divied into... 
    
    assert(n>=1);
    return 20*((n+1)*n + (n)*(n-1))/2;
    
  }

  int getNodeIJ(const int i,const int j,const int itr,const int ip0,const int ip1,const int ip2,const int ied0,const int ied1,const int ied2,const int ied0_sign,const int ied1_sign,const int ied2_sign,const int n) {

    // ======================================================
    // node indexing function for the sphere geodesic...
    // ======================================================

    //cout << " >> working on i,j: " << i << " " << j << endl;

    // recall surface->nsp = 12 + 30*(n-1) + 20*(n-1)*(n-2)/2;

    if (i == 0) {
      // bottom edge ied0....
      if (j == 0)
        return ip0;
      else if (j == n)
        return ip1;
      else {
        assert((j > 0)&&(j < n));
        if (ied0_sign == 1) {
          // return the j-1 point in ied0...
          return 12 + ied0*(n-1) + (j-1);
        }
        else {
          assert(ied0_sign == -1);
          return 12 + ied0*(n-1) + (n-j-1);
        }
      }
    }
    else if (i == n) {
      assert(j == 0);
      return ip2;
    }
    else {
      assert((i > 0)&&(i < n));
      if (j == 0) {
        // left edge ied1...
        if (ied1_sign == 1) {
          // return the i-1 point in ied1...
          return 12 + ied1*(n-1) + (i-1);
        }
        else {
          assert(ied1_sign == -1);
          return 12 + ied1*(n-1) + (n-i-1);
        }
      }
      else if (j == n-i) {
        // right edge ied2...
        if (ied2_sign == 1) {
          // return the i-1 point in ied2...
          return 12 + ied2*(n-1) + (i-1);
        }
        else {
          assert(ied2_sign == -1);
          return 12 + ied2*(n-1) + (n-i-1);
        }
      }
      else {
        assert((j > 0)&&(j < (n-i)));
        return 12 + 30*(n-1) + itr*(n-1)*(n-2)/2 + (n-1)*(n-2)/2 - (n-i)*(n-i-1)/2 + (j-1);
      }
    }
  }


  void addSphere(double (*xsp)[3],int (*spost)[3],const double xc[3],const double r,const int n,const bool flip) {
    
    // default for this sphere is as an outer fluid boundary -- i.e. flow on the inside, normals
    // pointing out. If you want the normal pointing in (i.e. flow over the outside of the sphere),
    // simply pass flip = true.

    // note: this routine assumes you have made space in xsp,spost,znost based on the count routines above...
    
    assert(n >= 1);

    // 12 unit vectors of the regular icosahedron...
    // for unit r, the edge length is...
    const double theta = 2.0*asin(2.0/sqrt(10.0+2*sqrt(5.0)));
    //cout << "theta: " << theta << endl;
    // polar angle is...
    const double cos_theta = cos(theta);
    const double sin_theta = sin(theta);
    const double xv[12][3] = {
      { 0.0, 0.0, 1.0 },
      { sin_theta, 0.0, cos_theta },
      { sin_theta*cos(2.0*M_PI/5.0),sin_theta*sin(2.0*M_PI/5.0),cos_theta },
      { sin_theta*cos(4.0*M_PI/5.0),sin_theta*sin(4.0*M_PI/5.0),cos_theta },
      { sin_theta*cos(6.0*M_PI/5.0),sin_theta*sin(6.0*M_PI/5.0),cos_theta },
      { sin_theta*cos(8.0*M_PI/5.0),sin_theta*sin(8.0*M_PI/5.0),cos_theta },
      { sin_theta*cos(1.0*M_PI/5.0),sin_theta*sin(1.0*M_PI/5.0),-cos_theta },
      { sin_theta*cos(3.0*M_PI/5.0),sin_theta*sin(3.0*M_PI/5.0),-cos_theta },
      { sin_theta*cos(5.0*M_PI/5.0),sin_theta*sin(5.0*M_PI/5.0),-cos_theta },
      { sin_theta*cos(7.0*M_PI/5.0),sin_theta*sin(7.0*M_PI/5.0),-cos_theta },
      { sin_theta*cos(9.0*M_PI/5.0),sin_theta*sin(9.0*M_PI/5.0),-cos_theta },
      { 0.0, 0.0, -1.0 }
    };
    
    const int nooed[30][2] = {
      { 0, 1 },
      { 1, 2 },
      { 0, 2 },
      { 2, 3 },
      { 0, 3 },
      { 3, 4 },
      { 0, 4 },
      { 4, 5 },
      { 0, 5 },
      { 1, 5 },
      { 1, 6 },
      { 2, 6 },
      { 6, 7 },
      { 2, 7 },
      { 3, 7 },
      { 7, 8 },
      { 3, 8 },
      { 4, 8 },
      { 8, 9 },
      { 4, 9 },
      { 5, 9 },
      { 9, 10 },
      { 5, 10 },
      { 1, 10 },
      { 6, 10 },
      { 6, 11 },
      { 7, 11 },
      { 8, 11 },
      { 9, 11 },
      { 10, 11 }
    };
    
    const int nootr[20][3] = {
      { 0, 1, 2 },
      { 0, 2, 3 },
      { 0, 3, 4 },
      { 0, 4, 5 },
      { 0, 5, 1 },
      { 1, 6, 2 },
      { 2, 6, 7 },
      { 2, 7, 3 },
      { 3, 7, 8 },
      { 3, 8, 4 },
      { 4, 8, 9 },
      { 4, 9, 5 },
      { 5, 9, 10 },
      { 5, 10, 1 },
      { 1, 10, 6 },
      { 6, 11, 7 },
      { 7, 11, 8 },
      { 8, 11, 9 },
      { 9, 11, 10 },
      { 10, 11, 6 }
    };
    
    const int edotr[20][3] = {
      { 0, 2, 1 },
      { 2, 4, 3 },
      { 4, 6, 5 },
      { 6, 8, 7 },
      { 8, 0, 9 },
      { 10, 1, 11 },
      { 11, 13, 12 },
      { 13, 3, 14 },
      { 14, 16, 15 },
      { 16, 5, 17 },
      { 17, 19, 18 },
      { 19, 7, 20 },
      { 20, 22, 21 },
      { 22, 9, 23 },
      { 23, 10, 24 },
      { 25, 12, 26 },
      { 26, 15, 27 },
      { 27, 18, 28 },
      { 28, 21, 29 },
      { 29, 24, 25 }
    };
    
    const int edotr_sign[20][3] = {
      { 1, 1, 1 },
      { 1, 1, 1 },
      { 1, 1, 1 },
      { 1, 1, 1 },
      { 1, 1, -1 },
      { 1, 1, -1 },
      { 1, 1, 1 },
      { 1, 1, -1 },
      { 1, 1, 1 },
      { 1, 1, -1 },
      { 1, 1, 1 },
      { 1, 1, -1 },
      { 1, 1, 1 },
      { 1, -1, -1 },
      { 1, 1, -1 },
      { 1, 1, -1 },
      { 1, 1, -1 },
      { 1, 1, -1 },
      { 1, 1, -1 },
      { 1, -1, -1 }
    };
    
    // the points on the sphere are ordered as follows:
    // nodes: first 12
    // edges: next 30*(N-1)
    // tris: next (N-1)*(N-2)/2
    
    int nsp = 0;

    // main nodes...
    for (int ip = 0; ip < 12; ++ip) {
      FOR_I3 xsp[nsp][i] = xc[i] + xv[ip][i]*r;
      ++nsp;
    }
    
    // edges...
    for (int ied = 0; ied < 30; ++ied) {
      const int ip0 = nooed[ied][0];
      const int ip1 = nooed[ied][1];
      // recall xv are already unit vectors...
      // here we use the so-called SLERP (spherical linear interpolation)...
      //const double theta = acos(DOT_PRODUCT(xv[ip0],xv[ip1])); // all 63.4349 degrees
      //cout << "theta: " << theta << endl;
      //const double inv_sin_theta = 1.0/sin(theta);
      for (int ii = 1; ii < n; ++ii) {
        const double w1 = double(ii)/double(n);
        double xv01[3]; FOR_I3 xv01[i] = (sin(w1*theta)*xv[ip1][i] + sin((1.0-w1)*theta)*xv[ip0][i])/sin_theta; // unit vector
        FOR_I3 xsp[nsp][i] = xc[i] + xv01[i]*r;
        ++nsp;
      }
    }
    assert(nsp == 12 + 30*(n-1));
    
    // internal nodes to each tri...
    for (int itr = 0; itr < 20; ++itr) {
      const int ip0 = nootr[itr][0];
      const int ip1 = nootr[itr][1];
      const int ip2 = nootr[itr][2];
      //const double theta = acos(DOT_PRODUCT(xv[ip0],xv[ip2])); // 63.4349 degrees
      for (int i = 1; i < n; ++i) {
        const double wgt2 = double(i)/double(n);
        for (int j = 1; j < n-i; ++j) {
          const double wgt1 = double(j)/double(n);
          double xp[3]; FOR_K3 xp[k] = sin((1.0-wgt1-wgt2)*theta)*xv[ip0][k] + sin(wgt1*theta)*xv[ip1][k] + sin(wgt2*theta)*xv[ip2][k];
          const double mag = MAG(xp);
          FOR_K3 xsp[nsp][k] = xc[k] + xp[k]*r/mag;
          ++nsp;
        }
      }
    }
    assert(nsp == 12 + 30*(n-1) + 20*(n-1)*(n-2)/2);
    assert(nsp == getSphereNodeCount(n));

    // =========================================================
    // now the tris...
    // =========================================================
    
    int nst = 0;
    
    for (int itr = 0; itr < 20; ++itr) {
      
      // stuff about this tri...
      const int ip0 = nootr[itr][0];
      const int ip1 = nootr[itr][1];
      const int ip2 = nootr[itr][2];
      const int ied0 = edotr[itr][0];
      const int ied1 = edotr[itr][1];
      const int ied2 = edotr[itr][2];
      const int ied0_sign = edotr_sign[itr][0];
      const int ied1_sign = edotr_sign[itr][1];
      const int ied2_sign = edotr_sign[itr][2];
      
      // upright tris...
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n-i; ++j) {
          //cout << "first i,j: " << i << " " << j << endl;
          spost[nst][0] = getNodeIJ(i,j,itr,ip0,ip1,ip2,ied0,ied1,ied2,ied0_sign,ied1_sign,ied2_sign,n);
          if (flip) {
            spost[nst][2] = getNodeIJ(i,j+1,itr,ip0,ip1,ip2,ied0,ied1,ied2,ied0_sign,ied1_sign,ied2_sign,n);
            spost[nst][1] = getNodeIJ(i+1,j,itr,ip0,ip1,ip2,ied0,ied1,ied2,ied0_sign,ied1_sign,ied2_sign,n);
          }
          else {
            spost[nst][1] = getNodeIJ(i,j+1,itr,ip0,ip1,ip2,ied0,ied1,ied2,ied0_sign,ied1_sign,ied2_sign,n);
            spost[nst][2] = getNodeIJ(i+1,j,itr,ip0,ip1,ip2,ied0,ied1,ied2,ied0_sign,ied1_sign,ied2_sign,n);
          }
          ++nst;
        }
      }

      for (int i = 1; i < n; ++i) {
        for (int j = 0; j < n-i; ++j) {
          //cout << "second i,j: " << i << " " << j << endl;
          spost[nst][0] = getNodeIJ(i,j,itr,ip0,ip1,ip2,ied0,ied1,ied2,ied0_sign,ied1_sign,ied2_sign,n);
          if (flip) {
            spost[nst][2] = getNodeIJ(i-1,j+1,itr,ip0,ip1,ip2,ied0,ied1,ied2,ied0_sign,ied1_sign,ied2_sign,n);
            spost[nst][1] = getNodeIJ(i,j+1,itr,ip0,ip1,ip2,ied0,ied1,ied2,ied0_sign,ied1_sign,ied2_sign,n);
          }
          else {
            spost[nst][1] = getNodeIJ(i-1,j+1,itr,ip0,ip1,ip2,ied0,ied1,ied2,ied0_sign,ied1_sign,ied2_sign,n);
            spost[nst][2] = getNodeIJ(i,j+1,itr,ip0,ip1,ip2,ied0,ied1,ied2,ied0_sign,ied1_sign,ied2_sign,n);
          }
          ++nst;
        }
      }

    }
    
    assert(nst == 20*((n+1)*n + (n)*(n-1))/2);
    assert(nst == getSphereTriCount(n));
    
    // finally, report the average spacing...

    double dmin = 1.0E+20;
    double dmax = 0.0;
    double davg = 0.0;
    double wgt = 0.0;
    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
        const double d = DIST(xsp[spost[ist][i]],xsp[spost[ist][(i+1)%3]]);
        davg += d;
        wgt += 1.0;
        dmax = max(dmax,d);
        dmin = min(dmin,d);
      }
    }
    davg /= wgt;
    cout << " > edge n: " << n << " dmin davg dmax: " << dmin << " " << davg << " " << dmax << ", davg/r: " << davg/r << endl;
    
  }

  void addSphere(double (*xsp)[3],const int n) {
    
    // just add sphere points about (0,0,0)...
    
    assert(n >= 1);

    // 12 unit vectors of the regular icosahedron...
    // for unit r, the edge length is...
    const double theta = 2.0*asin(2.0/sqrt(10.0+2*sqrt(5.0)));
    //cout << "theta: " << theta << endl;

    // polar angle is...
    const double cos_theta = cos(theta);
    const double sin_theta = sin(theta);
    const double xv[12][3] = {
      { 0.0, 0.0, 1.0 },
      { sin_theta, 0.0, cos_theta },
      { sin_theta*cos(2.0*M_PI/5.0),sin_theta*sin(2.0*M_PI/5.0),cos_theta },
      { sin_theta*cos(4.0*M_PI/5.0),sin_theta*sin(4.0*M_PI/5.0),cos_theta },
      { sin_theta*cos(6.0*M_PI/5.0),sin_theta*sin(6.0*M_PI/5.0),cos_theta },
      { sin_theta*cos(8.0*M_PI/5.0),sin_theta*sin(8.0*M_PI/5.0),cos_theta },
      { sin_theta*cos(1.0*M_PI/5.0),sin_theta*sin(1.0*M_PI/5.0),-cos_theta },
      { sin_theta*cos(3.0*M_PI/5.0),sin_theta*sin(3.0*M_PI/5.0),-cos_theta },
      { sin_theta*cos(5.0*M_PI/5.0),sin_theta*sin(5.0*M_PI/5.0),-cos_theta },
      { sin_theta*cos(7.0*M_PI/5.0),sin_theta*sin(7.0*M_PI/5.0),-cos_theta },
      { sin_theta*cos(9.0*M_PI/5.0),sin_theta*sin(9.0*M_PI/5.0),-cos_theta },
      { 0.0, 0.0, -1.0 }
    };
    
    const int nooed[30][2] = {
      { 0, 1 },
      { 1, 2 },
      { 0, 2 },
      { 2, 3 },
      { 0, 3 },
      { 3, 4 },
      { 0, 4 },
      { 4, 5 },
      { 0, 5 },
      { 1, 5 },
      { 1, 6 },
      { 2, 6 },
      { 6, 7 },
      { 2, 7 },
      { 3, 7 },
      { 7, 8 },
      { 3, 8 },
      { 4, 8 },
      { 8, 9 },
      { 4, 9 },
      { 5, 9 },
      { 9, 10 },
      { 5, 10 },
      { 1, 10 },
      { 6, 10 },
      { 6, 11 },
      { 7, 11 },
      { 8, 11 },
      { 9, 11 },
      { 10, 11 }
    };
    
    const int nootr[20][3] = {
      { 0, 1, 2 },
      { 0, 2, 3 },
      { 0, 3, 4 },
      { 0, 4, 5 },
      { 0, 5, 1 },
      { 1, 6, 2 },
      { 2, 6, 7 },
      { 2, 7, 3 },
      { 3, 7, 8 },
      { 3, 8, 4 },
      { 4, 8, 9 },
      { 4, 9, 5 },
      { 5, 9, 10 },
      { 5, 10, 1 },
      { 1, 10, 6 },
      { 6, 11, 7 },
      { 7, 11, 8 },
      { 8, 11, 9 },
      { 9, 11, 10 },
      { 10, 11, 6 }
    };
    
    // the points on the sphere are ordered as follows:
    // nodes: first 12
    // edges: next 30*(N-1)
    // tris: next (N-1)*(N-2)/2
    
    int nsp = 0;

    // main nodes...
    for (int ip = 0; ip < 12; ++ip) {
      FOR_I3 xsp[nsp][i] = xv[ip][i];
      ++nsp;
    }
    
    // edges...
    for (int ied = 0; ied < 30; ++ied) {
      const int ip0 = nooed[ied][0];
      const int ip1 = nooed[ied][1];
      // recall xv are already unit vectors...
      // here we use the so-called SLERP (spherical linear interpolation)...
      //const double theta = acos(DOT_PRODUCT(xv[ip0],xv[ip1])); // all 63.4349 degrees
      //cout << "theta: " << theta << endl;
      //const double inv_sin_theta = 1.0/sin(theta);
      for (int ii = 1; ii < n; ++ii) {
        const double w1 = double(ii)/double(n);
        double xv01[3]; FOR_I3 xv01[i] = (sin(w1*theta)*xv[ip1][i] + sin((1.0-w1)*theta)*xv[ip0][i])/sin_theta; // unit vector
        FOR_I3 xsp[nsp][i] = xv01[i];
        ++nsp;
      }
    }
    assert(nsp == 12 + 30*(n-1));
    
    // internal nodes to each tri...
    for (int itr = 0; itr < 20; ++itr) {
      const int ip0 = nootr[itr][0];
      const int ip1 = nootr[itr][1];
      const int ip2 = nootr[itr][2];
      //const double theta = acos(DOT_PRODUCT(xv[ip0],xv[ip2])); // 63.4349 degrees
      for (int i = 1; i < n; ++i) {
        const double wgt2 = double(i)/double(n);
        for (int j = 1; j < n-i; ++j) {
          const double wgt1 = double(j)/double(n);
          double xp[3]; FOR_K3 xp[k] = sin((1.0-wgt1-wgt2)*theta)*xv[ip0][k] + sin(wgt1*theta)*xv[ip1][k] + sin(wgt2*theta)*xv[ip2][k];
          const double mag = MAG(xp);
          FOR_K3 xsp[nsp][k] = xp[k]/mag;
          ++nsp;
        }
      }
    }
    assert(nsp == 12 + 30*(n-1) + 20*(n-1)*(n-2)/2);
    assert(nsp == getSphereNodeCount(n));

  }
  
}
