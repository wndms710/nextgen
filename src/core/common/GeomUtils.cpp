#include "GeomUtils.hpp"
#include "MiscUtils.hpp"
#include "WebUI.hpp"

namespace GeomUtils {

  double stretchingFunction(const int i,const int n,const double dn,const double dt) {
    assert((i >= 0)&&(i <= n));
    assert(dt >= dn);
    if (i == 0) return 0.0;
    if (dt == dn) return double(i)*dn;
    const double sf = pow(dt/dn,1.0/(double(n)-1.0)); assert(sf > 1.0);
    return dn*(pow(sf,i)-1.0)/(sf-1.0);
  }
  
  void createCirclePts(double (*xp)[3],const int n,const double xc[3],const double axis[3],const double r) {

    // note that axis need not be normalized

    double e1[3],e2[3];
    MiscUtils::getBestE1E2FromE0(e1,e2,axis);

    // populate circle of points into xp at specified index (clockwise ordering)
    for (int i = 0; i < n; ++i) {
      const double theta = double(i)/double(n)*2.0*M_PI; // note: put a point at top-dead-center
      FOR_J3 xp[i][j] = xc[j] + r*e1[j]*cos(theta) + r*e2[j]*sin(theta);
    }

  }

  void facetCircleToPoint(int (* spost)[3], int * znost, const int indexPt, const int indexCircle, const int n, const int zoneId, const bool flip) {

    // XXXXX: see version in surfer/SimpleSurface_primitives when condensing to GeomUtils in core

    // create the facets between a circle of points and a single point

    // flip dictates node ordering (normal direction convention)
    int pt1_inc = 1;
    int pt2_inc = 0;
    if (flip) {
      pt1_inc = 0;
      pt2_inc = 1;
    }

    for (int i=0, limit=n; i < limit; ++i) {
      spost[i][0] = indexPt;
      spost[i][1] = indexCircle + (i+pt1_inc)%n;
      spost[i][2] = indexCircle + (i+pt2_inc)%n;
      znost[i] = zoneId;
    }
  }


  void facetCircleToCircle(int (* spost)[3], int * znost, const int indexC0, const int indexC1, const int n, const int zoneId, const bool closed, const bool flip) {

    // create the facets between two circles of points, assuming they have the same n

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
      spost[2*i][0] = indexC0 + (i+pt1_inc)%n;
      spost[2*i][1] = indexC0 + (i+pt2_inc)%n;
      spost[2*i][2] = indexC1 + (i+pt1_inc)%n;
      znost[2*i] = zoneId;

      spost[2*i+1][0] = indexC0 + (i+pt2_inc)%n;;
      spost[2*i+1][1] = indexC1 + (i+pt2_inc)%n;
      spost[2*i+1][2] = indexC1 + (i+pt1_inc)%n;
      znost[2*i+1] = zoneId;
    }

  }


  int facetGap(int (*spost)[3],const vector<int> isp0Vec,const vector<int> isp1Vec,const double (* const xsp)[3],const bool loop) {

    // adds new tris to spost, and returns the new facet count

    assert(!isp0Vec.empty());
    assert(!isp1Vec.empty());

    if (loop) {

      int i0_start = -1;
      int i1_start = -1;
      double d2_start = HUGE_VAL;
      for (int i0 = 0, lim0 = isp0Vec.size(); i0 < lim0; i0++) { 
        const int isp0 = isp0Vec[i0];
        for (int i1 = 0, lim1 = isp1Vec.size(); i1 < lim1; i1++) { 
          const int isp1 = isp1Vec[i1];
          const double d2 = DIST2(xsp[isp0],xsp[isp1]);
          if (d2 < d2_start) { 
            i0_start = i0;
            i1_start = i1;
            d2_start = d2;
          }
        }//i1 
      }//i0

      assert(i0_start != -1);
      assert(i1_start != -1);

      int i0 = i0_start;
      int i0_next = i0+1;
      if (i0_next == int(isp0Vec.size()))
        i0_next = 0;
      assert(i0_next != i0);
      int i1 = i1_start;
      int i1_next = i1-1;
      if (i1_next == -1)
        i1_next = isp1Vec.size()-1;
      assert(i1_next != i1);

      int ist = 0;
      while ((i0 != i0_next)||(i1 != i1_next)) {

        if (i0 == i0_next) {

          assert(i0 == i0_start);
          spost[ist][0] = isp0Vec[i0];
          spost[ist][2] = isp1Vec[i1_next];
          spost[ist][1] = isp1Vec[i1];
          ++ist;
          i1 = i1_next;
          if (i1_next != i1_start) {
            --i1_next;
            if (i1_next == -1)
              i1_next = isp1Vec.size()-1;
          }

        }
        else if (i1 == i1_next) {

          assert(i1 == i1_start);
          spost[ist][0] = isp0Vec[i0];
          spost[ist][2] = isp0Vec[i0_next];
          spost[ist][1] = isp1Vec[i1];
          ++ist;
          i0 = i0_next;
          if (i0_next != i0_start) {
            ++i0_next;
            if (i0_next == int(isp0Vec.size()))
              i0_next = 0;
          }

        }
        else {

          assert(i0_next != i0);
          assert(i1_next != i1);

          const int isp0 = isp0Vec[i0];
          const int isp0_next = isp0Vec[i0_next];
          const int isp1 = isp1Vec[i1];
          const int isp1_next = isp1Vec[i1_next];

          // try the smallest next edge...
          
          const double d2_01n = DIST2(xsp[isp0],xsp[isp1_next]);
          const double d2_10n = DIST2(xsp[isp1],xsp[isp0_next]);
          if (d2_10n < d2_01n) {

            spost[ist][0] = isp0;
            spost[ist][2] = isp0_next;
            spost[ist][1] = isp1;
            ++ist;
            i0 = i0_next;
            if (i0_next != i0_start) {
              ++i0_next;
              if (i0_next == int(isp0Vec.size()))
                i0_next = 0;
            }
            
          }
          else {

            spost[ist][0] = isp0;
            spost[ist][2] = isp1_next;
            spost[ist][1] = isp1;
            ++ist;
            i1 = i1_next;
            if (i1_next != i1_start) {
              --i1_next;
              if (i1_next == -1)
                i1_next = isp1Vec.size()-1;
            }
            
          }

        }

      }

      return ist; // return the count of new facets

    }
    else {

      // find the closest...
      int i0_start = 0;
      int i1_start = 0;
      int i0_end = isp0Vec.size()-1;
      int i1_end = isp1Vec.size()-1;

      int i0 = i0_start;
      int i0_next = i0+1;
      int i1 = i1_start;
      int i1_next = i1+1;

      int ist = 0;
      while ((i0 != i0_end)||(i1 != i1_end)) {

        if (i0 == i0_end) {

          spost[ist][0] = isp0Vec[i0];
          spost[ist][2] = isp1Vec[i1_next];
          spost[ist][1] = isp1Vec[i1];
          ++ist;
          i1 = i1_next;
          ++i1_next;

        }
        else if (i1 == i1_end) {

          spost[ist][0] = isp0Vec[i0];
          spost[ist][2] = isp0Vec[i0_next];
          spost[ist][1] = isp1Vec[i1];
          ++ist;
          i0 = i0_next;
          ++i0_next;

        }
        else {

          assert(i0_next != i0);
          assert(i1_next != i1);

          const int isp0 = isp0Vec[i0];
          const int isp0_next = isp0Vec[i0_next];
          const int isp1 = isp1Vec[i1];
          const int isp1_next = isp1Vec[i1_next];

          // try the smallest next edge...

          const double d2_01n = DIST2(xsp[isp0],xsp[isp1_next]);
          const double d2_10n = DIST2(xsp[isp1],xsp[isp0_next]);
          if (d2_10n < d2_01n) {

            spost[ist][0] = isp0;
            spost[ist][2] = isp0_next;
            spost[ist][1] = isp1;
            ++ist;
            i0 = i0_next;
            ++i0_next;

          }
          else {

            spost[ist][0] = isp0;
            spost[ist][2] = isp1_next;
            spost[ist][1] = isp1;
            ++ist;
            i1 = i1_next;
            ++i1_next;

          }

        }

      }

      return ist;

    }

  }


  int facetGap(int (* spost)[3],const int isp0_f,const int isp0_l,const int isp1_f,const int isp1_l,const double (* const xsp)[3]) {

    // returns the new facet count
    // TODO: probably want this routine to operate on 2 consecutive vector<int> to
    // determine the isp's, and consider the case of a closed loop/annular region for example.

    int isp0 = isp0_f;
    int isp1 = isp1_f;
    int ist = 0;

    while ((isp0 != isp0_l)||(isp1 != isp1_l)) {

      if (isp0 == isp0_l) {

        assert(isp1 < isp1_l);
        spost[ist][0] = isp0;
        spost[ist][1] = isp1+1;
        spost[ist][2] = isp1;
        ++ist;
        ++isp1;

      }
      else if (isp1 == isp1_l) {

        assert(isp0 < isp0_l);
        spost[ist][0] = isp0;
        spost[ist][1] = isp0+1;
        spost[ist][2] = isp1;
        ++ist;
        ++isp0;

      }
      else {

        assert(isp0 < isp0_l);
        assert(isp1 < isp1_l);

        // check if isp1+1 is in the circle of isp0->isp0+1->isp1...
        double e0[3] = DIFF(xsp[isp0],xsp[isp1]);
        double mag_e0 = MAG(e0); assert(mag_e0 > 0.0);
        FOR_I3 e0[i] /= mag_e0;

        double det0;
        {

          double e1[3] = DIFF(xsp[isp0+1],xsp[isp0]);
          const double dp = DOT_PRODUCT(e1,e0);
          FOR_I3 e1[i] -= dp*e0[i];
          double mag_e1 = MAG(e1); assert(mag_e1 > 0.0);
          FOR_I3 e1[i] /= mag_e1;

          const double dx00 =
            (xsp[isp0][0]-xsp[isp1+1][0])*e0[0] +
            (xsp[isp0][1]-xsp[isp1+1][1])*e0[1] +
            (xsp[isp0][2]-xsp[isp1+1][2])*e0[2];
          const double dx01 =
            (xsp[isp0][0]-xsp[isp1+1][0])*e1[0] +
            (xsp[isp0][1]-xsp[isp1+1][1])*e1[1] +
            (xsp[isp0][2]-xsp[isp1+1][2])*e1[2];
          const double dx10 =
            (xsp[isp0+1][0]-xsp[isp1+1][0])*e0[0] +
            (xsp[isp0+1][1]-xsp[isp1+1][1])*e0[1] +
            (xsp[isp0+1][2]-xsp[isp1+1][2])*e0[2];
          const double dx11 =
            (xsp[isp0+1][0]-xsp[isp1+1][0])*e1[0] +
            (xsp[isp0+1][1]-xsp[isp1+1][1])*e1[1] +
            (xsp[isp0+1][2]-xsp[isp1+1][2])*e1[2];
          const double dx20 =
            (xsp[isp1][0]-xsp[isp1+1][0])*e0[0] +
            (xsp[isp1][1]-xsp[isp1+1][1])*e0[1] +
            (xsp[isp1][2]-xsp[isp1+1][2])*e0[2];
          const double dx21 =
            (xsp[isp1][0]-xsp[isp1+1][0])*e1[0] +
            (xsp[isp1][1]-xsp[isp1+1][1])*e1[1] +
            (xsp[isp1][2]-xsp[isp1+1][2])*e1[2];
          det0  =
            dx00* ( dx11*(dx20*dx20 + dx21*dx21) - dx21*( dx10*dx10 + dx11*dx11)) -
            dx01* ( dx10*(dx20*dx20 + dx21*dx21) - dx20*( dx10*dx10 + dx11*dx11)) +
            (dx00*dx00 + dx01*dx01)*( dx10*dx21 - dx11*dx20) ;
        }

        double det1;
        {

          double e1[3] = DIFF(xsp[isp1+1],xsp[isp1]);
          const double dp = DOT_PRODUCT(e1,e0);
          FOR_I3 e1[i] -= dp*e0[i];
          double mag_e1 = MAG(e1); assert(mag_e1 > 0.0);
          FOR_I3 e1[i] /= mag_e1;

          const double dx00 =
            (xsp[isp0][0]-xsp[isp0+1][0])*e0[0] +
            (xsp[isp0][1]-xsp[isp0+1][1])*e0[1] +
            (xsp[isp0][2]-xsp[isp0+1][2])*e0[2];
          const double dx01 =
            (xsp[isp0][0]-xsp[isp0+1][0])*e1[0] +
            (xsp[isp0][1]-xsp[isp0+1][1])*e1[1] +
            (xsp[isp0][2]-xsp[isp0+1][2])*e1[2];
          const double dx10 =
            (xsp[isp1+1][0]-xsp[isp0+1][0])*e0[0] +
            (xsp[isp1+1][1]-xsp[isp0+1][1])*e0[1] +
            (xsp[isp1+1][2]-xsp[isp0+1][2])*e0[2];
          const double dx11 =
            (xsp[isp1+1][0]-xsp[isp0+1][0])*e1[0] +
            (xsp[isp1+1][1]-xsp[isp0+1][1])*e1[1] +
            (xsp[isp1+1][2]-xsp[isp0+1][2])*e1[2];
          const double dx20 =
            (xsp[isp1][0]-xsp[isp0+1][0])*e0[0] +
            (xsp[isp1][1]-xsp[isp0+1][1])*e0[1] +
            (xsp[isp1][2]-xsp[isp0+1][2])*e0[2];
          const double dx21 =
            (xsp[isp1][0]-xsp[isp0+1][0])*e1[0] +
            (xsp[isp1][1]-xsp[isp0+1][1])*e1[1] +
            (xsp[isp1][2]-xsp[isp0+1][2])*e1[2];
          det1  =
            dx00* ( dx11*(dx20*dx20 + dx21*dx21) - dx21*( dx10*dx10 + dx11*dx11)) -
            dx01* ( dx10*(dx20*dx20 + dx21*dx21) - dx20*( dx10*dx10 + dx11*dx11)) +
            (dx00*dx00 + dx01*dx01)*( dx10*dx21 - dx11*dx20) ;
        }

        if (det1 > det0) {

          spost[ist][0] = isp0;
          spost[ist][1] = isp0+1;
          spost[ist][2] = isp1;
          ++ist;
          ++isp0;

        }
        else {

          spost[ist][0] = isp0;
          spost[ist][1] = isp1+1;
          spost[ist][2] = isp1;
          ++ist;
          ++isp1;

        }

      }

    }

    assert(ist == isp0_l-isp0_f+isp1_l-isp1_f);

    return ist; // return the count of new facets

  }


  void writeTecplot(const string& filename,const int (* const spost)[3],const int nst,const double (* const xsp)[3]) {

    cout << "writeTecplot: nst: " << nst << endl;

    // 1. figure out the range of isp's in spost...

    int isp0 = TWO_BILLION;
    int isp1 = -1;
    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
        const int isp = spost[ist][i];
        isp0 = min(isp0,isp);
        isp1 = max(isp1,isp);
      }
    }
    assert(isp0 < TWO_BILLION);
    assert(isp1 > -1);

    int nsp_max = isp1-isp0+1;
    int * sp_flag = new int[nsp_max];
    for (int isp = 0; isp < nsp_max; ++isp)
      sp_flag[isp] = -1;

    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
        const int isp = spost[ist][i];
        sp_flag[isp-isp0] = 0;
      }
    }

    int nsp = 0;
    for (int isp = 0; isp < nsp_max; ++isp)
      if (sp_flag[isp] == 0)
        sp_flag[isp] = nsp++;

    // 2. now write...

    FILE * fp = fopen(filename.c_str(),"w");
    fprintf(fp,"TITLE = \"%s\"\n",filename.c_str());
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");

    fprintf(fp,"ZONE T=\"%s\"\n",filename.c_str());
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nsp,nst);

    for (int isp = 0; isp < nsp_max; ++isp) {
      if (sp_flag[isp] >= 0) {
        fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp+isp0][0],xsp[isp+isp0][1],xsp[isp+isp0][2]);
      }
    }

    for (int ist = 0; ist < nst; ++ist) {
      fprintf(fp,"%d %d %d\n",
              sp_flag[spost[ist][0]-isp0]+1,
              sp_flag[spost[ist][1]-isp0]+1,
              sp_flag[spost[ist][2]-isp0]+1);
    }

    fclose(fp);

    delete[] sp_flag;

  }

  void writeTecplotZone(FILE * fp,const int (* const spost)[3],const int nst,const double (* const xsp)[3],const int * const flag,const int flag_value) {

    // 1. figure out the range of isp's in spost...

    int isp0 = TWO_BILLION;
    int isp1 = -1;
    int nst_flagged = 0;
    for (int ist = 0; ist < nst; ++ist) {
      if (flag[ist] == flag_value) {
        ++nst_flagged;
        FOR_I3 {
          const int isp = spost[ist][i];
          isp0 = min(isp0,isp);
          isp1 = max(isp1,isp);
        }
      }
    }
    assert(isp0 < TWO_BILLION);
    assert(isp1 > -1);

    int nsp_max = isp1-isp0+1;
    int * sp_flag = new int[nsp_max];
    for (int isp = 0; isp < nsp_max; ++isp) {
      sp_flag[isp] = -1;
    }

    for (int ist = 0; ist < nst; ++ist) {
      if (flag[ist] == flag_value) {
        FOR_I3 {
          const int isp = spost[ist][i];
          sp_flag[isp-isp0] = 0;
        }
      }
    }

    int nsp = 0;
    for (int isp = 0; isp < nsp_max; ++isp) {
      if (sp_flag[isp] == 0) {
        sp_flag[isp] = nsp++;
      }
    }

    // 2. now write...
    fprintf(fp,"ZONE T=\"Group %d\"\n",flag_value);
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nsp,nst_flagged);

    for (int isp = 0; isp < nsp_max; ++isp) {
      if (sp_flag[isp] >= 0) {
        fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp+isp0][0],xsp[isp+isp0][1],xsp[isp+isp0][2]);
      }
    }

    for (int ist = 0; ist < nst; ++ist) {
      if (flag[ist] == flag_value) {
        fprintf(fp,"%d %d %d\n",
                sp_flag[spost[ist][0]-isp0]+1,
                sp_flag[spost[ist][1]-isp0]+1,
                sp_flag[spost[ist][2]-isp0]+1);
      }
    }

    delete[] sp_flag;

  }
  
  void writeTecplot(const string& filename,const int (* const spost)[3],const int nst,const double (* const xsp)[3],const set<pair<int,int> >& istBitsSet) {

    FILE * fp = fopen(filename.c_str(),"w");
    fprintf(fp,"TITLE = \"%s\"\n",filename.c_str());
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");
    
    // 1. figure out the range of isp's in spost...

    int isp0 = TWO_BILLION;
    int isp1 = -1;
    for (set<pair<int,int> >::iterator it = istBitsSet.begin(); it != istBitsSet.end(); ++it) {
      const int ist = it->first;
      FOR_I3 {
        const int isp = spost[ist][i];
        isp0 = min(isp0,isp);
        isp1 = max(isp1,isp);
      }
    }
    assert(isp0 < TWO_BILLION);
    assert(isp1 > -1);
    
    int nsp_max = isp1-isp0+1;
    int * sp_flag = new int[nsp_max];
    for (int isp = 0; isp < nsp_max; ++isp) {
      sp_flag[isp] = -1;
    }
    
    for (set<pair<int,int> >::iterator it = istBitsSet.begin(); it != istBitsSet.end(); ++it) {
      const int ist = it->first;
      FOR_I3 {
        const int isp = spost[ist][i];
        sp_flag[isp-isp0] = 0;
      }
    }

    int nsp = 0;
    for (int isp = 0; isp < nsp_max; ++isp) {
      if (sp_flag[isp] == 0) {
        sp_flag[isp] = nsp++;
      }
    }
    
    // 2. now write...
    fprintf(fp,"ZONE T=\"istBitsSet\"\n");
    
    int nst_ = istBitsSet.size();
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nsp,nst_);
    
    for (int isp = 0; isp < nsp_max; ++isp) {
      if (sp_flag[isp] >= 0) {
        fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp+isp0][0],xsp[isp+isp0][1],xsp[isp+isp0][2]);
      }
    }
    
    for (set<pair<int,int> >::iterator it = istBitsSet.begin(); it != istBitsSet.end(); ++it) {
      const int ist = it->first;
      fprintf(fp,"%d %d %d\n",
              sp_flag[spost[ist][0]-isp0]+1,
              sp_flag[spost[ist][1]-isp0]+1,
              sp_flag[spost[ist][2]-isp0]+1);
    }

    fclose(fp);
    
    delete[] sp_flag;
    
  }

  void writeTecplotEdges(const string& filename,const double (* const xsp)[3],const vector<pair<int,int> >& linkVec) {
    
    FILE * fp = fopen(filename.c_str(),"w");
    fprintf(fp,"TITLE = \"%s\"\n",filename.c_str());
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");

    // 1. figure out the range of isp's in spost...

    int isp0 = TWO_BILLION;
    int isp1 = -1;
    for (int ii = 0; ii < linkVec.size(); ++ii) {
      int isp = linkVec[ii].first;
      isp0 = min(isp0,isp);
      isp1 = max(isp1,isp);
      isp = linkVec[ii].second;
    }
    assert(isp0 < TWO_BILLION);
    assert(isp1 > -1);
    
    int nsp_max = isp1-isp0+1;
    int * sp_flag = new int[nsp_max];
    for (int isp = 0; isp < nsp_max; ++isp) {
      sp_flag[isp] = -1;
    }
    
    for (int ii = 0; ii < linkVec.size(); ++ii) {
      int isp = linkVec[ii].first;
      sp_flag[isp-isp0] = 0;
      isp = linkVec[ii].second;
      sp_flag[isp-isp0] = 0;
    }

    int nsp = 0;
    for (int isp = 0; isp < nsp_max; ++isp) {
      if (sp_flag[isp] == 0) {
        sp_flag[isp] = nsp++;
      }
    }
    
    // 2. now write...
    fprintf(fp,"ZONE T=\"istBitsSet\"\n");
    
    int nst_ = linkVec.size();
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nsp,nst_);
    
    for (int isp = 0; isp < nsp_max; ++isp) {
      if (sp_flag[isp] >= 0) {
        fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp+isp0][0],xsp[isp+isp0][1],xsp[isp+isp0][2]);
      }
    }
    
    for (int ii = 0; ii < linkVec.size(); ++ii) {
      const int isp_first = linkVec[ii].first;
      const int isp_second = linkVec[ii].second;
      fprintf(fp,"%d %d %d\n",
              sp_flag[isp_first-isp0]+1,
              sp_flag[isp_second-isp0]+1,
              sp_flag[isp_second-isp0]+1);
    }
    
    fclose(fp);
    
    delete[] sp_flag;

  }
  
  void writeSbin(const string& filename,const int (* const spost)[3],const int nst,const double (* const xsp)[3]) {

    cout << "writeSbin: nst: " << nst << endl;

    // 1. figure out the range of isp's in spost. This
    // allows the sp allocations to be as small as possible...

    int isp0 = TWO_BILLION;
    int isp1 = -1;
    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
        const int isp = spost[ist][i];
        isp0 = min(isp0,isp);
        isp1 = max(isp1,isp);
      }
    }
    assert(isp0 < TWO_BILLION);
    assert(isp1 > -1);

    int nsp_max = isp1-isp0+1;
    int * sp_flag = new int[nsp_max];
    for (int isp = 0; isp < nsp_max; ++isp)
      sp_flag[isp] = -1;

    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
        const int isp = spost[ist][i];
        sp_flag[isp-isp0] = 0;
      }
    }

    int nsp = 0;
    for (int isp = 0; isp < nsp_max; ++isp)
      if (sp_flag[isp] == 0)
        sp_flag[isp] = nsp++;

    // 2. now write...

    FILE * fp = fopen(filename.c_str(),"wb");
    const int version = 1;
    fwrite(&version,sizeof(int),1,fp);
    const int count = 1;
    fwrite(&count,sizeof(int),1,fp);
    const string zoneName = "spost";
    const int length = zoneName.length();
    fwrite(&length,sizeof(int),1,fp);
    fwrite(zoneName.c_str(),sizeof(char),length,fp);
    cout << " > zone " << zoneName << endl;

    cout << " > nsp = " << nsp << endl;
    fwrite(&nsp,sizeof(int),1,fp);
    for (int isp = 0; isp < nsp_max; ++isp) {
      if (sp_flag[isp] >= 0) {
        fwrite(xsp+isp0+isp,sizeof(double),3,fp);
      }
    }

    cout << " > nst = " << nst << endl;
    fwrite(&nst,sizeof(int),1,fp);
    for (int ist = 0; ist < nst; ++ist) {
      int spost_local[3]; FOR_I3 spost_local[i] = sp_flag[spost[ist][i]-isp0];
      fwrite(spost_local,sizeof(int),3,fp);
    }

    const int izone = 0;
    for (int ist = 0; ist < nst; ++ist) {
      fwrite(&izone,sizeof(int),1,fp);
    }

    fclose(fp);

    delete[] sp_flag;

  }

  void writeTecplot(const string& filename,vector<pair<SimpleTri,int> >& triVec) {

    int my_ntri = triVec.size();
    int ntri;
    MPI_Reduce(&my_ntri,&ntri,1,MPI_INT,MPI_SUM,0,mpi_comm);
    
    FILE * fp;
    if ( mpi_rank == 0 ) {
      fp = fopen(filename.c_str(),"w");
      assert(fp != NULL);
      fprintf(fp,"TITLE = \"%s\"\n",filename.c_str());
      fprintf(fp,"VARIABLES = \"X\"\n");
      fprintf(fp,"\"Y\"\n");
      fprintf(fp,"\"Z\"\n");
      fprintf(fp,"ZONE T=\"%s\"\n",filename.c_str());
      fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",ntri*3,ntri);
    }
    else {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
    }
    
    for (int itri = 0; itri < my_ntri; ++itri) {
      fprintf(fp,"%18.15le %18.15le %18.15le\n",triVec[itri].first.x0[0],triVec[itri].first.x0[1],triVec[itri].first.x0[2]);
      fprintf(fp,"%18.15le %18.15le %18.15le\n",triVec[itri].first.x1[0],triVec[itri].first.x1[1],triVec[itri].first.x1[2]);
      fprintf(fp,"%18.15le %18.15le %18.15le\n",triVec[itri].first.x2[0],triVec[itri].first.x2[1],triVec[itri].first.x2[2]);
    }
    
    fclose(fp);

    if ( mpi_rank < mpi_size-1 ) {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }
    
    MPI_Barrier(mpi_comm);

    int offset = 0;
    if ( mpi_rank == 0 ) {
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
    }
    else {
      MPI_Status status;
      MPI_Recv(&offset,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
    }
    
    for (int itri = 0; itri < my_ntri; ++itri) {
      fprintf(fp,"%d %d %d\n",offset+1,offset+2,offset+3);
      offset += 3;
    }
    
    fclose(fp);
    
    if ( mpi_rank < mpi_size-1 ) {
      MPI_Send(&offset,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }
    
    MPI_Barrier(mpi_comm);
    
  }

  void writeTecplotGeom(const string& filename,vector<NewTri> &newTriVec,const double (*const xsp)[3]) {

    const int nst = newTriVec.size();
    cout << "writeTecplotGeom: nst: " << nst << endl;

    // 1. figure out the range of isp's in spost...

    int isp0 = TWO_BILLION;
    int isp1 = -1;
    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
        const int isp = newTriVec[ist].spost[i];
        isp0 = min(isp0,isp);
        isp1 = max(isp1,isp);
      }
    }
    assert(isp0 < TWO_BILLION);
    assert(isp1 > -1);

    int nsp_max = isp1-isp0+1;
    int * sp_flag = new int[nsp_max];
    for (int isp = 0; isp < nsp_max; ++isp)
      sp_flag[isp] = -1;

    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
        const int isp = newTriVec[ist].spost[i];
        sp_flag[isp-isp0] = 0;
      }
    }

    int nsp = 0;
    for (int isp = 0; isp < nsp_max; ++isp)
      if (sp_flag[isp] == 0)
        sp_flag[isp] = nsp++;

    // 2. now write...

    FILE * fp = fopen(filename.c_str(),"w");
    fprintf(fp,"TITLE = \"%s\"\n",filename.c_str());
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");

    fprintf(fp,"ZONE T=\"%s\"\n",filename.c_str());
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nsp,nst);

    for (int isp = 0; isp < nsp_max; ++isp) {
      if (sp_flag[isp] >= 0) {
        fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp+isp0][0],xsp[isp+isp0][1],xsp[isp+isp0][2]);
      }
    }

    for (int ist = 0; ist < nst; ++ist) {
      fprintf(fp,"%d %d %d\n",
              sp_flag[newTriVec[ist].spost[0]-isp0]+1,
              sp_flag[newTriVec[ist].spost[1]-isp0]+1,
              sp_flag[newTriVec[ist].spost[2]-isp0]+1);
    }

    fclose(fp);

    delete[] sp_flag;

  }

  void writePtsTecplot(const string& filename,const double (*const xp)[3],const int np) {

    // a handshake ascii point writer for distributed points...

    FILE * fp;
    if ( mpi_rank == 0 ) {
      fp = fopen(filename.c_str(),"w");
      assert(fp != NULL);
      fprintf(fp,"TITLE = \"flagged faces\"\n");
      fprintf(fp,"VARIABLES = \"X\"\n");
      fprintf(fp,"\"Y\"\n");
      fprintf(fp,"\"Z\"\n");
      //fprintf(fp,"\"DELTA\"\n");
      //fprintf(fp,"\"RANK\"\n");
    }
    else {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
    }

    for (int ip = 0; ip < np; ++ip) {
      fprintf(fp,"%18.15le %18.15le %18.15le\n", // %18.15le %d\n",
	      xp[ip][0],
	      xp[ip][1],
	      xp[ip][2]);
    }

    fclose(fp);

    if ( mpi_rank < mpi_size-1 ) {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }

    MPI_Barrier(mpi_comm);

  }

  int getAnnularCylNodeCount(const int nseg) {

    return nseg*4;

  }

  int getAnnularCylTriCount(const int nseg) {

    return nseg*8;

  }

  void addAnnularCyl(double (*xsp)[3],int (*spost)[3],const double x0[3],const double x1[3],const double r0,const double r1,const int nseg,const bool flip) {

    // flip=false: fluid inside, flip=true: fluid outside
    assert(r1 > r0);

    const double dx[3] = DIFF(x1,x0);
    double e1[3],e2[3];
    MiscUtils::getBestE1E2FromE0(e1,e2,dx);

    // populate circle of points into xp at specified index (clockwise ordering)
    for (int i = 0; i < nseg; ++i) {
      const double theta = double(i)/double(nseg)*2.0*M_PI; // note: put a point at top-dead-center
      const double cos_theta = cos(theta);
      const double sin_theta = sin(theta);
      FOR_J3 xsp[i][j]        = x0[j] + r0*e1[j]*cos_theta + r0*e2[j]*sin_theta;
      FOR_J3 xsp[i+nseg][j]   = x0[j] + r1*e1[j]*cos_theta + r1*e2[j]*sin_theta;
      FOR_J3 xsp[i+nseg*2][j] = x1[j] + r1*e1[j]*cos_theta + r1*e2[j]*sin_theta;
      FOR_J3 xsp[i+nseg*3][j] = x1[j] + r0*e1[j]*cos_theta + r0*e2[j]*sin_theta;
    }

    for (int i = 0; i < nseg; ++i) {
      const int ip1 = (i+1)%nseg;
      // inlet...
      spost[i*2][0] = i;
      spost[i*2][1] = ip1;
      spost[i*2][2] = ip1+nseg;
      spost[i*2+1][0] = i;
      spost[i*2+1][1] = ip1+nseg;
      spost[i*2+1][2] = i+nseg;
      // outside...
      spost[(i+nseg)*2][0] = i+nseg;
      spost[(i+nseg)*2][1] = ip1+nseg;
      spost[(i+nseg)*2][2] = ip1+2*nseg;
      spost[(i+nseg)*2+1][0] = i+nseg;
      spost[(i+nseg)*2+1][1] = ip1+2*nseg;
      spost[(i+nseg)*2+1][2] = i+2*nseg;
      // outlet...
      spost[(i+2*nseg)*2][0] = i+2*nseg;
      spost[(i+2*nseg)*2][1] = ip1+2*nseg;
      spost[(i+2*nseg)*2][2] = ip1+3*nseg;
      spost[(i+2*nseg)*2+1][0] = i+2*nseg;
      spost[(i+2*nseg)*2+1][1] = ip1+3*nseg;
      spost[(i+2*nseg)*2+1][2] = i+3*nseg;
      // inside...
      spost[(i+3*nseg)*2][0] = i+3*nseg;
      spost[(i+3*nseg)*2][1] = ip1+3*nseg;
      spost[(i+3*nseg)*2][2] = ip1;
      spost[(i+3*nseg)*2+1][0] = i+3*nseg;
      spost[(i+3*nseg)*2+1][1] = ip1;
      spost[(i+3*nseg)*2+1][2] = i;
    }

  }


  // TODO: condense to core along with stuff in stitch...

  int facetCircleToPointUseThisOne(int (* spost)[3], const int indexPt, const int indexCircle, const int n, const bool flip) {

    // create the facets between a circle of points and a single point

    // flip dictates node ordering (normal direction convention)
    int pt1_inc = 1;
    int pt2_inc = 0;
    if (flip) {
      pt1_inc = 0;
      pt2_inc = 1;
    }

    for (int i=0, limit=n; i < limit; ++i) {
      spost[i][0] = indexPt;
      spost[i][1] = indexCircle + (i+pt1_inc)%n;
      spost[i][2] = indexCircle + (i+pt2_inc)%n;
    }

    return n;

  }

  int facetGapUseThisOne(int (*spost)[3], const int isp0_start,const int nsp0,const int isp1_start, const int nsp1,const double (* const xsp)[3], const bool loop) {

    // adds new tris to spost, and returns the count
    // this one assumes the start is known (isp0_start, and isp1_start), AND
    // the nodes loop in the SAME direction. It is normally used in building
    // primitive geometries, not zipping up gaps (see vector-based version).

    if (loop) {

      int isp0 = isp0_start;
      int isp0_next = isp0+1;
      if (isp0_next == isp0_start+nsp0)
        isp0_next = isp0_start;
      int isp1 = isp1_start;
      int isp1_next = isp1+1;
      if (isp1_next == isp1_start+nsp1)
        isp1_next = isp1_start;

      int ist = 0;
      while ((isp0 != isp0_next)||(isp1 != isp1_next)) {

        if (isp0 == isp0_next) {

          assert(isp0 == isp0_start);
          spost[ist][0] = isp0;
          spost[ist][2] = isp1_next;
          spost[ist][1] = isp1;
          ++ist;
          isp1 = isp1_next;
          if (isp1_next != isp1_start) {
            ++isp1_next;
            if (isp1_next == isp1_start+nsp1)
              isp1_next = isp1_start;
          }
        }
        else if (isp1 == isp1_next) {

          assert(isp1 == isp1_start);
          spost[ist][0] = isp0;
          spost[ist][2] = isp0_next;
          spost[ist][1] = isp1;
          ++ist;
          isp0 = isp0_next;
          if (isp0_next != isp0_start) {
            ++isp0_next;
            if (isp0_next == isp0_start+nsp0)
              isp0_next = isp0_start;
          }
        }
        else {

          assert(isp0_next != isp0);
          assert(isp1_next != isp1);

          // try the smallest next edge...

          const double d2_01n = DIST2(xsp[isp0],xsp[isp1_next]);
          const double d2_10n = DIST2(xsp[isp1],xsp[isp0_next]);
          if (d2_10n < d2_01n) {

            spost[ist][0] = isp0;
            spost[ist][2] = isp0_next;
            spost[ist][1] = isp1;
            ++ist;
            isp0 = isp0_next;
            if (isp0_next != isp0_start) {
              ++isp0_next;
              if (isp0_next == isp0_start+nsp0)
                isp0_next = isp0_start;
            }

          }
          else {

            spost[ist][0] = isp0;
            spost[ist][2] = isp1_next;
            spost[ist][1] = isp1;
            ++ist;
            isp1 = isp1_next;
            if (isp1_next != isp1_start) {
              ++isp1_next;
              if (isp1_next == isp1_start+nsp1)
                isp1_next = isp1_start;
            }
          }

        }

      }

      return ist; // return the count of new facets

    }
    else {

      assert(0);

    }

    return -1;

  }

  void getHemisphereNodeAndTriCount(int &nsp,int &nst,const int ntheta) {

    // only requires ntheta...

    nsp = 0;
    nst = 0;
    nsp += 2; // bottom (center of sphere) and top (north pole) points
    const int nphi = max(1,ntheta/4);
    for (int iphi = 0; iphi < nphi; ++iphi) {
      const double phi = 0.5*M_PI*double(nphi-iphi)/double(nphi);
      int this_ntheta;
      if (iphi == 0) this_ntheta = max(3,ntheta);
      else this_ntheta = max(3,(int)ceil(sin(phi)*ntheta));
      nsp += this_ntheta;
      nst += 2*this_ntheta;
    }

  }

  void addHemisphere(double (*xsp)[3],int (*spost)[3],const double xp[3],const double np[3],const double rp,const int ntheta,const bool flip) {

    assert(flip == false); // for now

    // if you want to separate out the base, the first "ntheta" tris are the base...

    // two unit vectors in the equitorial plane
    double e1[3],e2[3];
    MiscUtils::getBestE1E2FromE0(e1,e2,np);
    const double np_mag = MAG(np); // may not be unit

    // add points and tris...
    int isp = 0;
    int ist = 0;
    // first 2 points are center and pole...
    FOR_I3 xsp[isp][i] = xp[i];
    FOR_I3 xsp[isp+1][i] = xp[i] + rp*np[i]/np_mag;
    isp += 2;
    const int nphi = max(1,ntheta/4);
    int isp_prev = -1;
    int ntheta_prev = -1;
    for (int iphi = 0; iphi < nphi; ++iphi) {
      const double phi = 0.5*M_PI*double(nphi-iphi)/double(nphi);
      const double h = rp*cos(phi);
      const double r = sqrt(rp*rp-h*h);
      int this_isp = isp;
      int this_ntheta;
      if (iphi == 0) this_ntheta = max(3,ntheta);
      else this_ntheta = max(3,(int)ceil(sin(phi)*ntheta));
      for (int itheta = 0; itheta < this_ntheta; ++itheta) {
        const double theta = 2.0*M_PI*double(itheta)/double(this_ntheta);
        FOR_I3 xsp[isp][i] = xp[i] + h*np[i]/np_mag + r*(e1[i]*cos(theta) + e2[i]*sin(theta));
        ++isp;
      }
      if (iphi == 0) {
        // this is the equatorial line of points, so connect to the center...
        ist += facetCircleToPointUseThisOne(spost+ist,0,this_isp,this_ntheta,false); // HACK: condense with GeomTools
      }
      else {
        ist += facetGapUseThisOne(spost+ist,this_isp,this_ntheta,isp_prev,ntheta_prev,xsp,true); // HACK: condense with GeomTools  periodic
      }
      // copy this_isp and this_ntheta for use the next time through...
      isp_prev = this_isp;
      ntheta_prev = this_ntheta;
    }

    // close the pole...
    ist += facetCircleToPointUseThisOne(spost+ist,1,isp_prev,ntheta_prev,true); // HACK: condense with GeomTools

    // check...
    int nsp_check,nst_check;
    getHemisphereNodeAndTriCount(nsp_check,nst_check,ntheta);
    assert(isp == nsp_check);
    assert(ist == nst_check);

  }

  int readXYZ(vector<double>& dVec,const string& filename) {

    // read a 3-column ascii file into a double vector...

    FILE * fp = fopen(filename.c_str(),"r");
    if (fp == NULL) {
      WUI(WARN,"3d points file cannot be found or opened: " << filename);
      return -1;
    }

    double bbmin[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
    double bbmax[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
    char line[256];
    int count = 0;
    while (fgets(line,256,fp) != NULL) {
      // skip possible comment lines...
      if (strncmp(line,"#",1) == 0)
        continue;
      double x,y,z;
      sscanf(line,"%lf %lf %lf\n",&x,&y,&z);
      ++count;
      bbmin[0] = min(bbmin[0],x);
      bbmax[0] = max(bbmax[0],x);
      bbmin[1] = min(bbmin[1],y);
      bbmax[1] = max(bbmax[1],y);
      bbmin[2] = min(bbmin[2],z);
      bbmax[2] = max(bbmax[2],z);
      dVec.push_back(x);
      dVec.push_back(y);
      dVec.push_back(z);
    }
    fclose(fp);

    if (mpi_rank == 0) {
      cout << " > readXYZ: " << filename << " pts: " << count <<
        " bbox 0: " << bbmin[0] << " " << bbmax[0] <<
        " 1: " << bbmin[1] << " " << bbmax[1] <<
        " 2: " << bbmin[2] << " " << bbmax[2] << endl;
    }

    return 0;

  }

  int readXY(vector<double>& dVec,const string& filename) {

    // read a 2-column ascii file into a double vector...

    FILE * fp = fopen(filename.c_str(),"r");
    if (fp == NULL) {
      WUI(WARN,"2d points file cannot be found or opened: " << filename);
      return -1;
    }

    double bbmin[2] = { HUGE_VAL, HUGE_VAL };
    double bbmax[2] = { -HUGE_VAL, -HUGE_VAL };
    char line[256];
    int count = 0;
    while (fgets(line,256,fp) != NULL) {
      // skip possible comment lines...
      if (strncmp(line,"#",1) == 0)
        continue;
      double x,y;
      sscanf(line,"%lf %lf\n",&x,&y);
      ++count;
      bbmin[0] = min(bbmin[0],x);
      bbmax[0] = max(bbmax[0],x);
      bbmin[1] = min(bbmin[1],y);
      bbmax[1] = max(bbmax[1],y);
      dVec.push_back(x);
      dVec.push_back(y);
    }
    fclose(fp);

    if (mpi_rank == 0) {
      cout << " > readXY: " << filename << " pts: " << count <<
        " bbox 0: " << bbmin[0] << " " << bbmax[0] <<
        " 1: " << bbmin[1] << " " << bbmax[1] << endl;
    }

    return 0;

  }

  bool getSphereIntersections(double xmp[2],const double yz[2],const double x[3],const double r) {

    const double d2 = r*r*0.999999999 - (x[1]-yz[0])*(x[1]-yz[0]) - (x[2]-yz[1])*(x[2]-yz[1]);
    if ( d2 > 0 ) {
      xmp[0] = x[0] - sqrt(d2);
      xmp[1] = x[0] + sqrt(d2);
      return true;
    }

    return false;

  }

  bool getMinusSphereIntersection(double &xm,const double yz[2],const double x[3],const double r) {

    const double d2 = r*r*0.999999999 - (x[1]-yz[0])*(x[1]-yz[0]) - (x[2]-yz[1])*(x[2]-yz[1]);
    if ( d2 > 0 ) {
      xm = x[0] - sqrt(d2);
      return true;
    }

    return false;

  }

  bool getPlusSphereIntersection(double &xp,const double yz[2],const double x[3],const double r) {

    // note one more 9 here to bias finding plus intersection...

    const double d2 = r*r*0.9999999999 - (x[1]-yz[0])*(x[1]-yz[0]) - (x[2]-yz[1])*(x[2]-yz[1]);
    if ( d2 > 0 ) {
      xp = x[0] + sqrt(d2);
      return true;
    }

    return false;

  }

  bool calcCappedCylinderIntersections(double xmp[2],const double yz[2],const double x0[3],const double x1[3],const double r) {

    // returns true if the passed cylinder capped with spheres of the same radius has in and out intersections...

    // in the y-z plane, we can check if an intersection exists...
    const double dx01[3] = DIFF(x1,x0);

    // when cylinder is aligned with x-ray, use limiting spheres...
    const double a = dx01[1]*dx01[1] + dx01[2]*dx01[2];
    if (a == 0.0) {
      // calculate the intersections based on one of the spheres...
      if (getSphereIntersections(xmp,yz,x0,r)) {
        // then adjust one of the intersections depending on whether x1 is less or greater in x...
        if (x1[0] < x0[0]) {
          // the minus intersection needs adjustment...
          xmp[0] -= (x0[0]-x1[0]);
        }
        else {
          // the plus intersection needs adjustment...
          xmp[1] += (x1[0]-x0[0]);
        }
        return true;
      }
      else {
        return false;
      }
    }

    const double dp = (yz[0]-x0[1])*dx01[1] + (yz[1]-x0[2])*dx01[2];
    const double b = (yz[0]-x0[1])*(yz[0]-x0[1]) + (yz[1]-x0[2])*(yz[1]-x0[2]);
    if (b - dp*dp/a < 0.9999999999*r*r) {

      // there is a cylinder intersection...

      const double dx01_2 = dx01[0]*dx01[0];
      const double A = dx01_2 + a;

      // quadratic coeffs...
      const double aq = 1.0 - dx01_2/A;
      if (aq < 1.0E-12) {
        // assuming O(1) normalization here...
        return false;
      }
      const double bq = -2.0*( aq*x0[0] + dp*dx01[0]/A );
      const double cq = aq*x0[0]*x0[0] + (2.0*dx01[0]*x0[0]-dp)*dp/A + b - r*r;

      const double t = bq*bq - 4.0*aq*cq;
      if (t <= 0.0) {
        //cout << "Surface::calcCappedCylinderIntersections: unexpected t. Is this near zero?: " << t << endl;
        return false;
      }

      const double sqrt_t = sqrt(t);
      xmp[0] = 0.5*(-bq - sqrt_t)/aq;
      xmp[1] = 0.5*(-bq + sqrt_t)/aq;

      const double sm = ((xmp[0]-x0[0])*dx01[0] + dp)/A;
      if (sm < 0.0) {
        // the minus intersection is on the < 0.0 side...
        if (!getMinusSphereIntersection(xmp[0],yz,x0,r))
          return false;
      }
      else if (sm > 1.0) {
        // the minus intersection is on the > 1.0 side...
        if (!getMinusSphereIntersection(xmp[0],yz,x1,r))
          return false;
      }

      const double sp = ((xmp[1]-x0[0])*dx01[0] + dp)/A;
      if (sp < 0.0) {
        // the plus intersection is on the < 0.0 side...
        if (!getPlusSphereIntersection(xmp[1],yz,x0,r)) {
          assert(0); // if we get a minus, we should get a plus: see tolerance settings in getPlusSphereIntersection...
          return false;
        }
      }
      else if (sp > 1.0) {
        // the plus intersection is on the > 1.0 side...
        if (!getPlusSphereIntersection(xmp[1],yz,x1,r)) {
          assert(0); // see comment above
          return false;
        }
      }

      return true;

    }

    return false;

  }

  int getMinusCylinderIntersection(double &xm,const double yz[2],const double x0[3],const double x1[3],const double r) {

    // returns
    // -1: no intersection
    // 0: intersection on part of cylinder below x0
    // 1: intersection on part of cylinder above x1
    // 2: intersection between x0 and x1 -- proper intersection

    // in the y-z plane, we can check if an intersection exists...
    const double dx01[3] = DIFF(x1,x0);

    // no intersection when cylinder is aligned with x-ray...
    if (dx01[0] == 0.0)
      return -1;

    const double a = dx01[1]*dx01[1] + dx01[2]*dx01[2];
    if (a == 0.0)
      return -1;

    const double dp = (yz[0]-x0[1])*dx01[1] + (yz[1]-x0[2])*dx01[2];
    const double b = (yz[0]-x0[1])*(yz[0]-x0[1]) + (yz[1]-x0[2])*(yz[1]-x0[2]);
    if (b - dp*dp/a < 0.9999999999*r*r) {

      // there is a cylinder intersection...

      const double dx01_2 = dx01[0]*dx01[0];
      const double A = dx01_2 + a;

      // quadratic coeffs...
      const double aq = 1.0 - dx01_2/A; assert(aq > 0.0);
      const double bq = -2.0*( aq*x0[0] + dp*dx01[0]/A );
      const double cq = aq*x0[0]*x0[0] + (2.0*dx01[0]*x0[0]-dp)*dp/A + b - r*r;

      const double t = bq*bq - 4.0*aq*cq; assert(t > 0.0);
      const double sqrt_t = sqrt(t);
      xm = 0.5*(-bq - sqrt_t)/aq;

      const double s = ((xm-x0[0])*dx01[0] + dp)/A;
      if (s < 0.0)
        return 0;
      else if (s > 1.0)
        return 1;
      else
        return 2;

    }
    return -1;

  }

  int getPlusCylinderIntersection(double &xp,const double yz[2],const double x0[3],const double x1[3],const double r) {

    // returns
    // -1: no intersection
    // 0: intersection on part of cylinder below x0
    // 1: intersection on part of cylinder above x1
    // 2: intersection between x0 and x1 -- proper intersection

    // in the y-z plane, we can check if an intersection exists...
    const double dx01[3] = DIFF(x1,x0);

    // no intersection when cylinder is aligned with x-ray...
    if (dx01[0] == 0.0)
      return -1;

    const double a = dx01[1]*dx01[1] + dx01[2]*dx01[2];
    if (a == 0.0)
      return -1;

    const double dp = (yz[0]-x0[1])*dx01[1] + (yz[1]-x0[2])*dx01[2];
    const double b = (yz[0]-x0[1])*(yz[0]-x0[1]) + (yz[1]-x0[2])*(yz[1]-x0[2]);
    if (b - dp*dp/a < 0.9999999999*r*r) {

      // there is a cylinder intersection...

      const double dx01_2 = dx01[0]*dx01[0];
      const double A = dx01_2 + a;

      // quadratic coeffs...
      const double aq = 1.0 - dx01_2/A; assert(aq > 0.0);
      const double bq = -2.0*( aq*x0[0] + dp*dx01[0]/A );
      const double cq = aq*x0[0]*x0[0] + (2.0*dx01[0]*x0[0]-dp)*dp/A + b - r*r;

      const double t = bq*bq - 4.0*aq*cq; assert(t > 0.0);
      const double sqrt_t = sqrt(t);
      xp = 0.5*(-bq + sqrt_t)/aq;

      const double s = ((xp-x0[0])*dx01[0] + dp)/A;
      if (s < 0.0)
        return 0;
      else if (s > 1.0)
        return 1;
      else
        return 2;

    }
    return -1;

  }

  bool calcBloatedTriIntersections(double xmp[2],const double yz[2],const double x0[3],const double x1[3],const double x2[3],const double delta) {

    // the above routine was an effort to manage the amount of intersections being calculated. It
    // was not very robust. Here we do all cylinders first, then the triangular caps...

    bool found = false;
    double this_xmp[2];
    if (calcCappedCylinderIntersections(this_xmp,yz,x0,x1,delta)) {
      if (!found) {
        found = true;
        xmp[0] = this_xmp[0];
        xmp[1] = this_xmp[1];
      }
      else {
        xmp[0] = min(xmp[0],this_xmp[0]);
        xmp[1] = max(xmp[1],this_xmp[1]);
      }
    }
    if (calcCappedCylinderIntersections(this_xmp,yz,x1,x2,delta)) {
      if (!found) {
        found = true;
        xmp[0] = this_xmp[0];
        xmp[1] = this_xmp[1];
      }
      else {
        xmp[0] = min(xmp[0],this_xmp[0]);
        xmp[1] = max(xmp[1],this_xmp[1]);
      }
    }
    if (calcCappedCylinderIntersections(this_xmp,yz,x2,x0,delta)) {
      if (!found) {
        found = true;
        xmp[0] = this_xmp[0];
        xmp[1] = this_xmp[1];
      }
      else {
        xmp[0] = min(xmp[0],this_xmp[0]);
        xmp[1] = max(xmp[1],this_xmp[1]);
      }
    }

    // finally consider the triangle caps...

    bool found_minus = false;
    bool found_plus = false;

    const double Ax = ((x1)[1]-(x0)[1])*((x2)[2]-(x0)[2])-((x1)[2]-(x0)[2])*((x2)[1]-(x0)[1]);
    if (Ax > 1.0E-8*delta*delta) {
      // the projected area is positive and large enough to test the intersection of
      // the lifted area...
      const double Ay = ((x1)[2]-(x0)[2])*((x2)[0]-(x0)[0])-((x1)[0]-(x0)[0])*((x2)[2]-(x0)[2]);
      const double Az = ((x1)[0]-(x0)[0])*((x2)[1]-(x0)[1])-((x1)[1]-(x0)[1])*((x2)[0]-(x0)[0]);
      const double factor = delta/sqrt(Ax*Ax+Ay*Ay+Az*Az);
      // for the minus intersections, shift by -n/nmag*delta...
      {
        const double x0_[3] = { x0[0] - Ax*factor, x0[1] - Ay*factor, x0[2] - Az*factor };
        const double x1_[3] = { x1[0] - Ax*factor, x1[1] - Ay*factor, x1[2] - Az*factor };
        const double x2_[3] = { x2[0] - Ax*factor, x2[1] - Ay*factor, x2[2] - Az*factor };
        const double A1x = (yz[0]-(x0_)[1])*((x2_)[2]-(x0_)[2])-(yz[1]-(x0_)[2])*((x2_)[1]-(x0_)[1]);
        const double A2x = ((x1_)[1]-(x0_)[1])*(yz[1]-(x0_)[2])-((x1_)[2]-(x0_)[2])*(yz[0]-(x0_)[1]);
        const double A0x = Ax - A1x - A2x;
        if ((A0x >= -1.0E-6*Ax)&&(A1x >= -1.0E-6*Ax)&&(A2x >= -1.0E-6*Ax)) {
          // the minus intersection is in the offset tri...
          xmp[0] = (A0x*x0_[0] + A1x*x1_[0] + A2x*x2_[0])/Ax;
          found_minus = true;
        }
      }
      // for the plus...
      {
        const double x0_[3] = { x0[0] + Ax*factor, x0[1] + Ay*factor, x0[2] + Az*factor };
        const double x1_[3] = { x1[0] + Ax*factor, x1[1] + Ay*factor, x1[2] + Az*factor };
        const double x2_[3] = { x2[0] + Ax*factor, x2[1] + Ay*factor, x2[2] + Az*factor };
        const double A1x = (yz[0]-(x0_)[1])*((x2_)[2]-(x0_)[2])-(yz[1]-(x0_)[2])*((x2_)[1]-(x0_)[1]);
        const double A2x = ((x1_)[1]-(x0_)[1])*(yz[1]-(x0_)[2])-((x1_)[2]-(x0_)[2])*(yz[0]-(x0_)[1]);
        const double A0x = Ax - A1x - A2x;
        if ((A0x >= -1.0E-6*Ax)&&(A1x >= -1.0E-6*Ax)&&(A2x >= -1.0E-6*Ax)) {
          // the plus intersection is in the offset tri...
          xmp[1] = (A0x*x0_[0] + A1x*x1_[0] + A2x*x2_[0])/Ax;
          found_plus = true;
        }
      }
    }
    else if (Ax < -1.0E-8*delta*delta) {
      // the projected area is positive and large enough to test the intersection of
      // the lifted area...
      const double Ay = ((x1)[2]-(x0)[2])*((x2)[0]-(x0)[0])-((x1)[0]-(x0)[0])*((x2)[2]-(x0)[2]);
      const double Az = ((x1)[0]-(x0)[0])*((x2)[1]-(x0)[1])-((x1)[1]-(x0)[1])*((x2)[0]-(x0)[0]);
      const double factor = delta/sqrt(Ax*Ax+Ay*Ay+Az*Az);
      // for the minus intersections, shift by -n/nmag*delta...
      {
        const double x0_[3] = { x0[0] + Ax*factor, x0[1] + Ay*factor, x0[2] + Az*factor };
        const double x1_[3] = { x1[0] + Ax*factor, x1[1] + Ay*factor, x1[2] + Az*factor };
        const double x2_[3] = { x2[0] + Ax*factor, x2[1] + Ay*factor, x2[2] + Az*factor };
        const double A1x = (yz[0]-(x0_)[1])*((x2_)[2]-(x0_)[2])-(yz[1]-(x0_)[2])*((x2_)[1]-(x0_)[1]);
        const double A2x = ((x1_)[1]-(x0_)[1])*(yz[1]-(x0_)[2])-((x1_)[2]-(x0_)[2])*(yz[0]-(x0_)[1]);
        const double A0x = Ax - A1x - A2x;
        if ((A0x <= 1.0E-6*Ax)&&(A1x <= 1.0E-6*Ax)&&(A2x <= 1.0E-6*Ax)) {
          // the minus intersection is in the offset tri...
          xmp[0] = (A0x*x0_[0] + A1x*x1_[0] + A2x*x2_[0])/Ax;
          found_minus = true;
        }
      }
      // for the plus...
      {
        const double x0_[3] = { x0[0] - Ax*factor, x0[1] - Ay*factor, x0[2] - Az*factor };
        const double x1_[3] = { x1[0] - Ax*factor, x1[1] - Ay*factor, x1[2] - Az*factor };
        const double x2_[3] = { x2[0] - Ax*factor, x2[1] - Ay*factor, x2[2] - Az*factor };
        const double A1x = (yz[0]-(x0_)[1])*((x2_)[2]-(x0_)[2])-(yz[1]-(x0_)[2])*((x2_)[1]-(x0_)[1]);
        const double A2x = ((x1_)[1]-(x0_)[1])*(yz[1]-(x0_)[2])-((x1_)[2]-(x0_)[2])*(yz[0]-(x0_)[1]);
        const double A0x = Ax - A1x - A2x;
        if ((A0x <= 1.0E-6*Ax)&&(A1x <= 1.0E-6*Ax)&&(A2x <= 1.0E-6*Ax)) {
          // the plus intersection is in the offset tri...
          xmp[1] = (A0x*x0_[0] + A1x*x1_[0] + A2x*x2_[0])/Ax;
          found_plus = true;
        }
      }
    }

    if (found_plus && found_minus) {
      return true;
    }
    else if (found_plus) {
      assert(found);
      return true;
    }
    else if (found_minus) {
      assert(found);
      return true;
    }
    else {
      return found;
    }

  }


  // also: this should probably not use NewTri, but int3 instead or something like that...

  bool addTrisMarchingFront2d(vector<NewTri> &newTriVec,set<pair<int,int> > &frontSet,const double (* const x)[2],const int n) {

    // expected number of new tris is nEdges-2 if all are used
    //const int n_expected = frontSet.size()-2;
    // int count[n];
    // for (int i = 0; i < n; ++i)
    // count[i] = 0;

    int * count = new int[n];
    for (int i = 0; i < n; ++i) count[i] = 0;

    // int failed[n];
    // for (int i = 0; i < n; ++i)
    // failed[i] = 0;

    int * failed = new int[n];
    for (int i = 0; i < n; ++i) failed[i] = 0;

    int index = 0;

    for (set< std::pair<int,int> >::iterator iter = frontSet.begin(); iter != frontSet.end(); ++iter) {
      int i0 = iter->first; assert((i0 >= 0)&&(i0 < n));
      int i1 = iter->second; assert((i1 >= 0)&&(i1 < n));

      count[i0] += 1;
      count[i1] += 1;
    }
    for (int i = 0; i < n; ++i) {
      if (count[i] != 2) {
        // cout << "count[i] < 2: " << count[i] << " debug_count=" << debug_count << endl;
        // assert(0);
        CWARN("some nodes on the passed loops have a valence != 2. Cannot cap these edges; skipping.");
        return -1;
      }
    }

    while (!frontSet.empty()) {

      // take the first element from the front...
      set< std::pair<int,int> >::iterator iter = frontSet.begin();
      const int i0 = iter->first; assert((i0 >= 0)&&(i0 < n)); assert(count[i0] > 0);
      const int i1 = iter->second; assert((i1 >= 0)&&(i1 < n)); assert(count[i1] > 0);

      // and remove the element...

      frontSet.erase(iter);
      --count[i0];
      --count[i1];

      ++index;

      // look for a matching edge pair...

      set< std::pair<int,int> >::iterator iter0 = frontSet.find( std::pair<int,int>(i1,i0) );
      if ( (iter0 != frontSet.end()) && ((count[i0] == 1)||(count[i1] == 1)) ) {

        frontSet.erase(iter0);
        --count[i0];
        --count[i1];

        // why do this?
        assert(0);
        /*
          newTrisVec.push_back(NewTri(i0,i1,i0,izone,nsz));
          CWARN("Adding a colapsed tri: " << i0+1 << " " << i1+1 << " " << i0+1);
        */
      }
      else {

        int i2;
        for (i2 = 0; i2 < n; ++i2) if ((i2 != i0)&&(i2 != i1)&&(count[i2] > 0)&&(failed[i2] < index)) {

            // check orientation of node...
            const double dx02[2] = { x[i2][0]-x[i0][0], x[i2][1]-x[i0][1] };
            const double dx21[2] = { x[i1][0]-x[i2][0], x[i1][1]-x[i2][1] };
            if ((dx02[0]*dx21[1] - dx02[1]*dx21[0]) < 0.0) {

              // check that no other valid node is in-circle...
              int i;
              for (i = 0; i < n; ++i) if ((i != i0)&&(i != i1)&&(i != i2)&&(count[i] > 0)&&(failed[i] < index)) {
                  const double dx00 = x[i0][0] - x[i][0];
                  const double dx01 = x[i0][1] - x[i][1];
                  const double dx10 = x[i1][0] - x[i][0];
                  const double dx11 = x[i1][1] - x[i][1];
                  const double dx20 = x[i2][0] - x[i][0];
                  const double dx21 = x[i2][1] - x[i][1];
                  const double det  =
                    dx00* ( dx11*(dx20*dx20 + dx21*dx21) - dx21*( dx10*dx10 + dx11*dx11)) -
                    dx01* ( dx10*(dx20*dx20 + dx21*dx21) - dx20*( dx10*dx10 + dx11*dx11)) +
                    (dx00*dx00 + dx01*dx01)*( dx10*dx21 - dx11*dx20) ;

                  if (det > 0.0) {

                    // cout << " > > potentially failed in-circle test, i: " << i+1 << " det: " << det << " i0: " << i0+1 << " i1: " << i1+1 << endl;

                    // point i is in-circle, but if tri i0->i1->i is invalid, then this is not
                    // a problem...

                    const double dx0i[2] = { x[i][0]-x[i0][0], x[i][1]-x[i0][1] };
                    const double dxi1[2] = { x[i1][0]-x[i][0], x[i1][1]-x[i][1] };
                    if ((dx0i[0]*dxi1[1] - dx0i[1]*dxi1[0]) >= 0.0) {
                      // cout << " > > potential check: failed orientation test, i: " << i+1 << endl;
                      failed[i] = index;
                      continue;
                    }

                    // cout << "1 looking for front: " << i+1 << " " << i0+1 << endl;

                    // new edge i0->i first...
                    set< std::pair<int,int> >::iterator iter0 = frontSet.find( std::pair<int,int>(i,i0) );
                    if (iter0 == frontSet.end()) {

                      // cout << "1 did not find: " << i+1 << " " << i0+1 << endl;

                      // we did not find the front edge i->i0, so we need to ensure that new
                      // front edge i0->i does NOT intersect any other front edges...
                      set< std::pair<int,int> >::iterator iter_check;
                      for (iter_check = frontSet.begin(); iter_check != frontSet.end(); ++iter_check) {
                        const int i0_check = iter_check->first; assert((i0_check >= 0)&&(i0_check < n)); assert(count[i0_check] > 0);
                        const int i1_check = iter_check->second; assert((i1_check >= 0)&&(i1_check < n)); assert(count[i1_check] > 0);

                        // cout << "1 doing edge check: " << i0 << " " << i << " against " << i0_check << " " << i1_check << endl;

                        if ((i0_check != i0)&&(i0_check != i)&&(i1_check != i0)&&(i1_check != i)) {
                          const double dx0_check[2] = { x[i0_check][0]-x[i0][0], x[i0_check][1]-x[i0][1] };
                          const double dx01_check[2] = { x[i1_check][0]-x[i0_check][0], x[i1_check][1]-x[i0_check][1] };
                          const double denom = dx01_check[0]*dx0i[1] - dx01_check[1]*dx0i[0];
                          if (denom != 0.0) {
                            const double s = (dx01_check[0]*dx0_check[1] - dx01_check[1]*dx0_check[0])/denom;
                            if ((s >= 0.0)&&(s <= 1.0)) {
                              // the first line intersected. Now check the second...
                              const double t = (dx0i[0]*dx0_check[1] - dx0i[1]*dx0_check[0])/denom;
                              if ((t >= 0.0)&&(t <= 1.0)) {
                                // cout << " > > potential check: got intersection with edge: " << i0_check+1 << " " << i1_check+1 << endl;
                                failed[i] = index;
                                break;
                              }
                            }
                          }
                        }
                      }
                      if (iter_check != frontSet.end()) {
                        assert(failed[i] == index);
                        // cout << " > > failed in-circle ignored, i: " << i+1 << endl;
                        continue;
                      }

                    }

                    // cout << "2 looking for front: " << i1+1 << " " << i+1 << endl;

                    // new edge i->i1 next...
                    set< std::pair<int,int> >::iterator iter1 = frontSet.find( std::pair<int,int>(i1,i) );
                    if (iter1 == frontSet.end()) {

                      // cout << "2 did not find: " << i1+1 << " " << i+1 << endl;

                      // we did not find the front edge i1->i, so we need to ensure that new
                      // front edge i->i1 does NOT intersect any other front edges...
                      set< std::pair<int,int> >::iterator iter_check;
                      for (iter_check = frontSet.begin(); iter_check != frontSet.end(); ++iter_check) {
                        const int i0_check = iter_check->first; assert((i0_check >= 0)&&(i0_check < n)); assert(count[i0_check] > 0);
                        const int i1_check = iter_check->second; assert((i1_check >= 0)&&(i1_check < n)); assert(count[i1_check] > 0);

                        // cout << "2 doing edge check: " << i << " " << i1 << " against " << i0_check << " " << i1_check << endl;

                        if ((i0_check != i1)&&(i0_check != i)&&(i1_check != i1)&&(i1_check != i)) {
                          const double dx0_check[2] = { x[i0_check][0]-x[i][0], x[i0_check][1]-x[i][1] };
                          const double dx01_check[2] = { x[i1_check][0]-x[i0_check][0], x[i1_check][1]-x[i0_check][1] };
                          const double denom = dx01_check[0]*dxi1[1] - dx01_check[1]*dxi1[0];
                          if (denom != 0.0) {
                            const double s = (dx01_check[0]*dx0_check[1] - dx01_check[1]*dx0_check[0])/denom;
                            if ((s >= 0.0)&&(s <= 1.0)) {
                              // the first line intersected. Now check the second...
                              const double t = (dxi1[0]*dx0_check[1] - dxi1[1]*dx0_check[0])/denom;
                              if ((t >= 0.0)&&(t <= 1.0)) {
                                // if (debug) cout << " > > potential check: got intersection with edge: " << i0_check+1 << " " << i1_check+1 << endl;
                                failed[i] = index;
                                break;
                              }
                            }
                          }
                        }
                      }
                      if (iter_check != frontSet.end()) {
                        assert(failed[i] == index);
                        // if (debug) cout << " > > failed in-circle ignored, i: " << i+1 << endl;
                        continue;
                      }
                    }

                    break;

                  }
                }

              // if we did not get through all the points, then i2 will not work, so continue...
              if (i < n) {
                // cout << " > > failed in-circle test, i: " << i+1 << endl;
                failed[i2] = index;
                continue;
              }

              // we are considering forming the tri i0->i1->i2...

              // new edge i0->i2 first...
              set< std::pair<int,int> >::iterator iter0 = frontSet.find( std::pair<int,int>(i2,i0) );
              if (iter0 == frontSet.end()) {
                // we did not find the front edge i2->i0, so we need to ensure that new
                // front edge i0->i2 does NOT intersect any other front edges...
                // we already have dx02...
                set< std::pair<int,int> >::iterator iter_check;
                for (iter_check = frontSet.begin(); iter_check != frontSet.end(); ++iter_check) {
                  const int i0_check = iter_check->first; assert((i0_check >= 0)&&(i0_check < n)); assert(count[i0_check] > 0);
                  const int i1_check = iter_check->second; assert((i1_check >= 0)&&(i1_check < n)); assert(count[i1_check] > 0);
                  if ((i0_check != i0)&&(i0_check != i2)&&(i1_check != i0)&&(i1_check != i2)) {
                    const double dx0_check[2] = { x[i0_check][0]-x[i0][0], x[i0_check][1]-x[i0][1] };
                    const double dx01_check[2] = { x[i1_check][0]-x[i0_check][0], x[i1_check][1]-x[i0_check][1] };
                    const double denom = dx01_check[0]*dx02[1] - dx01_check[1]*dx02[0];
                    if (denom != 0.0) {
                      const double s = (dx01_check[0]*dx0_check[1] - dx01_check[1]*dx0_check[0])/denom;
                      if ((s >= 0.0)&&(s <= 1.0)) {
                        // the first line intersected. Now check the second...
                        const double t = (dx02[0]*dx0_check[1] - dx02[1]*dx0_check[0])/denom;
                        if ((t >= 0.0)&&(t <= 1.0)) {
                          // cout << " > > got intersection with edge: " << i0_check+1 << " " << i1_check+1 << endl;
                          failed[i2] = index;
                          break;
                        }
                      }
                    }
                  }
                }
                if (iter_check != frontSet.end())
                  continue;
              }

              // new edge i2->i1 next...
              set< std::pair<int,int> >::iterator iter1 = frontSet.find( std::pair<int,int>(i1,i2) );
              if (iter1 == frontSet.end()) {
                // we did not find the front edge i1->i2, so we need to ensure that new
                // front edge i2->i1 does NOT intersect any other front edges...
                set< std::pair<int,int> >::iterator iter_check;
                for (iter_check = frontSet.begin(); iter_check != frontSet.end(); ++iter_check) {
                  const int i0_check = iter_check->first; assert((i0_check >= 0)&&(i0_check < n)); assert(count[i0_check] > 0);
                  const int i1_check = iter_check->second; assert((i1_check >= 0)&&(i1_check < n)); assert(count[i1_check] > 0);
                  if ((i0_check != i1)&&(i0_check != i2)&&(i1_check != i1)&&(i1_check != i2)) {
                    const double dx0_check[2] = { x[i0_check][0]-x[i2][0], x[i0_check][1]-x[i2][1] };
                    const double dx01_check[2] = { x[i1_check][0]-x[i0_check][0], x[i1_check][1]-x[i0_check][1] };
                    const double denom = dx01_check[0]*dx21[1] - dx01_check[1]*dx21[0];
                    if (denom != 0.0) {
                      const double s = (dx01_check[0]*dx0_check[1] - dx01_check[1]*dx0_check[0])/denom;
                      if ((s >= 0.0)&&(s <= 1.0)) {
                        // the first line intersected. Now check the second...
                        const double t = (dx21[0]*dx0_check[1] - dx21[1]*dx0_check[0])/denom;
                        if ((t >= 0.0)&&(t <= 1.0)) {
                          // cout << " > > got intersection with edge: " << i0_check+1 << " " << i1_check+1 << endl;
                          failed[i2] = index;
                          break;
                        }
                      }
                    }
                  }
                }
                if (iter_check != frontSet.end())
                  continue;

              }

              // if we made it here, then i2 is our guy!...
              newTriVec.push_back(NewTri(i0,i1,i2));

              // delete the current front element from the front...
              if (iter0 == frontSet.end()) {
                frontSet.insert( std::pair<int,int>(i0,i2) );
                ++count[i0];
                ++count[i2];
              }
              else {
                frontSet.erase(iter0);
                --count[i0];
                --count[i2];
              }

              if (iter1 == frontSet.end()) {
                frontSet.insert( std::pair<int,int>(i2,i1) );
                ++count[i1];
                ++count[i2];
              }
              else {
                frontSet.erase(iter1);
                --count[i1];
                --count[i2];
              }

              break;

            }
          }

        if (i2 == n) {
          cout << endl << "WARNING triangulate marching front: no valid i2" << endl;
          cout.flush();
          // assert(0);
          return false;
        }
      }
    }

    delete[] failed;
    
    // int ierr = 0;
    for (int i = 0; i < n; ++i) {
      if (count[i] != 0) {
        CWARN("front is finished but some nodes were not processed");
        return -1;
      }
    }

    delete[] count;
    return 0;

  }

  void coordinateTransform(double x_t[3],const double x[3],const double r0,const double cos_theta0,const double sin_theta0,
			   const double cyl_center[3],const COORD_TRANSFORM_TYPE transform_type) {

    switch (transform_type) {
    case TRANSFORM_CYL_X:
      // (x,y,z) -> (x,r,r0*theta)
      x_t[0] = (x[0]-cyl_center[0]);
      x_t[1] = sqrt((x[1]-cyl_center[1])*(x[1]-cyl_center[1])+(x[2]-cyl_center[2])*(x[2]-cyl_center[2]));
      x_t[2] = atan2((x[2]-cyl_center[2])*cos_theta0-(x[1]-cyl_center[1])*sin_theta0,
		     (x[1]-cyl_center[1])*cos_theta0+(x[2]-cyl_center[2])*sin_theta0)*r0;
      break;

    case TRANSFORM_CYL_Y:
      CWARN("TRANSFORM_CYL_Y not supported yet");
      break;

    case TRANSFORM_CYL_Z:
      // (x,y,z) -> (r,r0*theta,z)
      x_t[0] = sqrt((x[0]-cyl_center[0])*(x[0]-cyl_center[0])+(x[1]-cyl_center[1])*(x[1]-cyl_center[1]));
      x_t[1] = atan2((x[1]-cyl_center[1])*cos_theta0-(x[0]-cyl_center[0])*sin_theta0,
		     (x[0]-cyl_center[0])*cos_theta0+(x[1]-cyl_center[1])*sin_theta0)*r0;
      x_t[2] = (x[2]-cyl_center[2]);
      break;

    default:
      // just copy
      FOR_I3 x_t[i] = x[i];
    }
  }

} // namespace
