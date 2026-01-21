#include "SimpleSurface.hpp"
#include "WebUI.hpp"

#include "GeomUtils.hpp"
#include "../stitch/CuttableVoronoiData.hpp"

/*
  #ifdef WITH_CHRONO
  #include <chrono>
  #endif
*/

void lloydIterate2d(vector<NewTri>& triVec,double (*xp)[3],double * deltap,int * stoxp,const int np,
                    const int (* const spost)[3],const int (* const stost)[3],const int nst,const double (* const xsp)[3],const int nsp) {

  /*
    #ifdef WITH_CHRONO
    // timing...
    double time_buf[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    #endif
  */

  int maxiter = 10;
  double zero = 0.0;

  int * st_flag = new int[nst];
  double (*dxp)[3] = new double[np][3]; // lloyd motion
  int *stoxp_new = new int[np];
  int * xpost_i = new int[nst+1];
  int * xpost_v = new int[np];

  CuttableVoronoiData cvd;
  stack<int> stack;
  vector<pair<double,int> > nbrVec;
  map<const pair<int,int>,int> edgeMap;
  map<const int,int> nodeMap;
  int noost[3],edost[3];
  map<const int,double*> x0Map;
  set<int> nbrSet;

  // also the last iteration will put the nbrs in these...
  int * nboxp_i = new int[np+1];
  nboxp_i[0] = 0;
  vector<int> nboxp_v;

  int debug_iter = 0;
  int iter = 0;
  int done = 0;
  while (done != 2) {

    /*
      #ifdef WITH_CHRONO
      auto t0 = std::chrono::high_resolution_clock::now();
      #endif
    */

    ++iter;

    // build reverse-lookup xpost_i/v...
    for (int ist = 0; ist < nst; ++ist)
      xpost_i[ist+1] = 0;
    for (int ip = 0; ip < np; ++ip) {
      const int ist = stoxp[ip];
      assert((ist >= 0)&&(ist < nst));
      ++xpost_i[ist+1];
    }
    xpost_i[0] = 0;
    for (int ist = 0; ist < nst; ++ist) xpost_i[ist+1] += xpost_i[ist];
    assert(xpost_i[nst] == np);
    for (int ip = 0; ip < np; ++ip) {
      const int ist = stoxp[ip];
      xpost_v[xpost_i[ist]++] = ip;
    }
    for (int ist = nst-1; ist > 0; --ist)
      xpost_i[ist] = xpost_i[ist-1];
    xpost_i[0] = 0;

    /*
    // check...
    for (int ist = 0; ist < nst; ++ist) {
    for (int xos = xpost_i[ist]; xos != xpost_i[ist+1]; ++xos) {
    const int ip = xpost_v[xos];
    assert(stoxp[ip] == ist);
    }
    }
    */

    // st_flag is used to build the surface patches. Start it off
    // with -1...

    for (int ist = 0; ist < nst; ++ist) st_flag[ist] = -1;

    /*
      #ifdef WITH_CHRONO
      auto t1 = std::chrono::high_resolution_clock::now();
      #endif
    */

    // loop through our ip's...
    int my_count[3] = { 0, 0, 0 };
    for (int ip = 0; ip < np; ++ip) {

      //cout << "working on ip: " << ip << " stoxp: " << stoxp[ip] << endl;

      /*
        #ifdef WITH_CHRONO
        auto t20 = std::chrono::high_resolution_clock::now();
        #endif
      */

      ++debug_iter;

      assert(cvd.nno == 0);
      assert(cvd.ned == 0);
      assert(nbrVec.empty());
      assert(edgeMap.empty());
      assert(nodeMap.empty());
      assert(stack.empty());

      assert(st_flag[stoxp[ip]] != ip);
      st_flag[stoxp[ip]] = ip;
      stack.push(-stoxp[ip]-1); // push as negative because BOTH tris and nbrs for the tri we live on...

      while (!stack.empty()) {
        int ist = stack.top(); stack.pop();
        if (ist < 0) {
          // -----------------------------------------------------------------
          // a negative ist means both tris AND nbrs, so add the tris here...
          // -----------------------------------------------------------------
          ist = -ist-1;
          // use edge-matching to set as much of the noost and edost that we can...
          FOR_I3 noost[i] = -1;
          FOR_I3 {
            edost[i] = -1;
            // here the edgeMap is used to find edges in terms of surface node pairs...
            map<const pair<int,int>,int>::iterator it = edgeMap.find(pair<int,int>(spost[ist][(i+1)%3],spost[ist][i]));
            if (it != edgeMap.end()) {
              edost[i] = it->second;
              edgeMap.erase(it);
              if (noost[i] == -1) noost[i] = cvd.nooed[edost[i]][1];
              else assert(noost[i] == cvd.nooed[edost[i]][1]);
              if (noost[(i+1)%3] == -1) noost[(i+1)%3] = cvd.nooed[edost[i]][0];
              else assert(noost[(i+1)%3] == cvd.nooed[edost[i]][0]);
            }
          }
          // step 3: now noost and edost are set from edges where possible. At this point, complete any nodes not set...
          FOR_I3 {
            if (noost[i] == -1) {
              map<const int,int>::iterator it = nodeMap.find(spost[ist][i]);
              if (it == nodeMap.end()) {
                const int ino = cvd.new_node();
                noost[i] = ino;
                nodeMap[spost[ist][i]] = ino;
                FOR_J3 cvd.x_no[ino][j] = xsp[spost[ist][i]][j] - xp[ip][j];
              }
              else {
                noost[i] = it->second;
              }
            }
          }
          // step 4: finally complete the edges...
          FOR_I3 {
            if (edost[i] == -1) {
              const int ied = cvd.new_edge();
              edgeMap[pair<int,int>(spost[ist][i],spost[ist][(i+1)%3])] = ied;
              cvd.nooed[ied][0] = noost[i];
              cvd.nooed[ied][1] = noost[(i+1)%3];
              cvd.faoed[ied][0] = ist;
              const int ist_nbr = stost[ist][i];
              if (ist_nbr >= 0) cvd.faoed[ied][1] = -1;
              else cvd.faoed[ied][1] = -2; // used to mark external boundary
            }
            else {
              const int ied = edost[i];
              assert(cvd.nooed[ied][1] == noost[i]);
              assert(cvd.nooed[ied][0] == noost[(i+1)%3]);
              assert(cvd.faoed[ied][0] != ist); // no edge can have the same tri on both sides
              assert(cvd.faoed[ied][1] == -1);
              cvd.faoed[ied][1] = ist;
            }
          }
        }
        // -----------------------------------------------------------------
        // always add nbrs...
        // -----------------------------------------------------------------
        assert(st_flag[ist] == ip);
        for (int xos = xpost_i[ist]; xos != xpost_i[ist+1]; ++xos) {
          const int ip_nbr = xpost_v[xos];
          if (ip_nbr != ip) {
            const double d2 = DIST2(xp[ip_nbr],xp[ip]);
            assert(d2 > 0.0); // possible to seed 2 points at the same location?
            if (d2 < deltap[ip]*deltap[ip]) nbrVec.push_back(pair<double,int>(d2,ip_nbr));
          }
        }
        // if the neighboring tri is closer than 0.5*deltap[ip], then we need to add it AND
        // its nbrs. If it is closer than deltap[ip], then we just need to add the nbrs.
        FOR_I3 {
          const int ist_nbr = stost[ist][i];
          assert(ist_nbr != ist);
          if ((ist_nbr >= 0)&&(st_flag[ist_nbr] != ip)) {
            st_flag[ist_nbr] = ip; // do not ever revisit this tri
            // the 3 nodes of ist_nbr...
            const int isp0 = spost[ist_nbr][0];
            const int isp1 = spost[ist_nbr][1];
            const int isp2 = spost[ist_nbr][2];
            const double d2 = MiscUtils::getPointToTriDist2(xp[ip],xsp[isp0],xsp[isp1],xsp[isp2]);
            // if this tri is closer than 0.5*deltap[ip], we add it to the surface, and if less than
            // deltap[ip], we need to potentially add the nbrs (but not the surface)...
            if (d2 < 0.2501*deltap[ip]*deltap[ip])  {
              // both tri and nbrs...
              stack.push(-ist_nbr-1);
            }
            else if (d2 < deltap[ip]*deltap[ip])  {
              // just nbrs...
              stack.push(ist_nbr);
            }
          }
        }
      }
      nodeMap.clear();
      edgeMap.clear();

      /*
        #ifdef WITH_CHRONO
        auto t21 = std::chrono::high_resolution_clock::now();
        #endif
      */

      cvd.setD2Max();
      sort(nbrVec.begin(),nbrVec.end());
      // we now have a sorted list of nbrs. Cut the cvd against the nbrs until it cannot possibly
      // be cut anymore...
      for (int in = 0,nnbr=nbrVec.size(); in < nnbr; ++in) {
        // recall that the .first contains the d2 nbr distance. and cvd.d2_max contains
        // the maximum node radius. as soon as d2 >= 4*cvd.d2_max, there cannot be any
        // more nbrs that cut this cv...
        if (nbrVec[in].first > 4.0*cvd.d2_max) break; // done: all remaining nbrs are too far away to possibly cut the current vd...
        const int ip_nbr = nbrVec[in].second;
        double dn[3]; FOR_I3 dn[i] = 0.5*(xp[ip_nbr][i] - xp[ip][i]);
        // only cut if we are inside the 6 paraboloids...
        if ((2.0*dn[0] > cvd.Lmax[0])||
            (2.0*dn[1] > cvd.Lmax[1])||
            (2.0*dn[2] > cvd.Lmax[2])||
            (2.0*dn[0] < cvd.Lmin[0])||
            (2.0*dn[1] < cvd.Lmin[1])||
            (2.0*dn[2] < cvd.Lmin[2]) ) {
          continue;
        }
        // if we got here, we are cutting...
        cvd.cut_surf(dn,-8-ip_nbr); // note that local nbrs use a -8 convention. cvd.d2_max updated automatically
      }
      nbrVec.clear();

      // compute the centroid of the remaining parts...
      bool has_seed_boundary = false;
      //bool has_boundary = false;
      double dxc[3] = { 0.0, 0.0, 0.0 };
      double area2_sum = 0.0;
      assert(x0Map.empty());
      for (int ied = 0; ied < cvd.ned; ++ied) {
        if (cvd.faoed[ied][0] >= 0) {
          map<const int,double*>::const_iterator it = x0Map.find(cvd.faoed[ied][0]);
          if (it == x0Map.end()) {
            // on the first time we visit any face, store a node to use in all
            // the tris with other edges (the face is planar and convex, so this is fine)...
            x0Map[cvd.faoed[ied][0]] = cvd.x_no[cvd.nooed[ied][0]];
          }
          else if (cvd.x_no[cvd.nooed[ied][1]] != it->second) {
            const double normal2[3] = TRI_NORMAL_2(it->second,cvd.x_no[cvd.nooed[ied][0]],cvd.x_no[cvd.nooed[ied][1]]);
            const double area2 = MAG(normal2);
            FOR_I3 dxc[i] += area2*(it->second[i]+cvd.x_no[cvd.nooed[ied][0]][i]+cvd.x_no[cvd.nooed[ied][1]][i]);
            area2_sum += area2;
          }
        }
        if (cvd.faoed[ied][1] >= 0) {
          map<const int,double*>::const_iterator it = x0Map.find(cvd.faoed[ied][1]);
          if (it == x0Map.end()) {
            x0Map[cvd.faoed[ied][1]] = cvd.x_no[cvd.nooed[ied][1]];
          }
          else if (cvd.x_no[cvd.nooed[ied][0]] != it->second) {
            const double normal2[3] = TRI_NORMAL_2(it->second,cvd.x_no[cvd.nooed[ied][1]],cvd.x_no[cvd.nooed[ied][0]]);
            const double area2 = MAG(normal2);
            FOR_I3 dxc[i] += area2*(it->second[i]+cvd.x_no[cvd.nooed[ied][0]][i]+cvd.x_no[cvd.nooed[ied][1]][i]);
            area2_sum += area2;
          }
        }
        else if (cvd.faoed[ied][1] == -1) {
          has_seed_boundary = true;
          break;
        }
        else {
          // must be a -2: this is ok (flagged region boundary)...
          assert(cvd.faoed[ied][1] == -2);
          //has_boundary = true;
        }
      }

      if (has_seed_boundary) {
        // we are too small. Increase deltap[ip], and do NOT apply any motion to the point. The remaining
        // patch is arbitrarily large...
        ++my_count[0];
        assert(done == 0);
        deltap[ip] *= 1.5;
        stoxp_new[ip] = stoxp[ip];
        FOR_I3 dxp[ip][i] = 0.0;
      }
      else {
        // Even though we have cut away all the original seed boundary of this patch, there still
        // may be nbrs out there that can cut some of our corners off if deltap was not large
        // enough...
        const double delta_half = sqrt(cvd.d2_max); // cvd.d2_max stores the furthest distance of any node
        if (deltap[ip] <= 2.0*delta_half) {
          // other neighbors might exist that could still cut this volume. Find them in
          // the next iteration...
          ++my_count[1];
          assert(done == 0);
          deltap[ip] = 2.5*delta_half; // could be just 2.0 here? - use a little larger so face checks pass for sure
          stoxp_new[ip] = stoxp[ip];
          FOR_I3 dxp[ip][i] = 0.0;
        }
        else {
          ++my_count[2];
          // this is good. Do some more work to decide where to move this point...
          if (done == 0) {
            // regular lloyd iteration pass
            assert(area2_sum > 0.0);
            FOR_I3 dxc[i] /= area2_sum*3.0;
            deltap[ip] = 2.5*delta_half;
            // in addition to lloyd, add random perturbation...
            //FOR_I3 dxc[i] += (double(rand())/double(RAND_MAX)-0.5)*deltap[ip]*factor;
            // now we need to figure out which patch we are closest to...
            int ist_closest = -1;
            double d2_closest;
            for (int ied = 0; ied < cvd.ned; ++ied) {
              if (cvd.faoed[ied][0] >= 0) {
                map<const int,double*>::const_iterator it = x0Map.find(cvd.faoed[ied][0]);
                assert(it != x0Map.end());
                if ((cvd.x_no[cvd.nooed[ied][0]] != it->second)&&(cvd.x_no[cvd.nooed[ied][1]] != it->second)) {
                  double this_dxp[3];
                  MiscUtils::getClosestPointOnTriRobust(this_dxp,dxc,it->second,cvd.x_no[cvd.nooed[ied][0]],cvd.x_no[cvd.nooed[ied][1]]);
                  const double this_d2 = DIST2(this_dxp,dxc);
                  if ((ist_closest == -1)||(this_d2 < d2_closest)) {
                    FOR_I3 dxp[ip][i] = this_dxp[i];
                    d2_closest = this_d2;
                    ist_closest = cvd.faoed[ied][0];
                  }
                }
              }
              if (cvd.faoed[ied][1] >= 0) {
                map<const int,double*>::const_iterator it = x0Map.find(cvd.faoed[ied][1]);
                assert(it != x0Map.end());
                if ((cvd.x_no[cvd.nooed[ied][0]] != it->second)&&(cvd.x_no[cvd.nooed[ied][1]] != it->second)) {
                  double this_dxp[3];
                  MiscUtils::getClosestPointOnTriRobust(this_dxp,dxc,it->second,cvd.x_no[cvd.nooed[ied][0]],cvd.x_no[cvd.nooed[ied][1]]);
                  const double this_d2 = DIST2(this_dxp,dxc);
                  if ((ist_closest == -1)||(this_d2 < d2_closest)) {
                    FOR_I3 dxp[ip][i] = this_dxp[i];
                    d2_closest = this_d2;
                    ist_closest = cvd.faoed[ied][1];
                  }
                }
              }
            }
            // make sure we found one!...
            assert(ist_closest != -1);
            stoxp_new[ip] = ist_closest;
          }
          else {
            assert(done == 1);
            // this is a recompute of the last time...
            assert(nbrSet.empty());
            for (int ied = 0; ied < cvd.ned; ++ied) {
              if (cvd.faoed[ied][0] < 0) {
                assert(cvd.faoed[ied][0] <= -8);
                nbrSet.insert(-cvd.faoed[ied][0]-8);
              }
            }
            assert(!nbrSet.empty());
            for (set<int>::const_iterator it = nbrSet.begin(); it != nbrSet.end(); ++it) nboxp_v.push_back(*it);
            nboxp_i[ip+1] = nboxp_v.size();
            nbrSet.clear();
          }
        }
      }
      // clear cvd for next time...
      x0Map.clear();
      cvd.clear();

      /*
        #ifdef WITH_CHRONO
        auto t22 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed20_21 = t21-t20;
        std::chrono::duration<double> elapsed21_22 = t22-t21;
        time_buf[0] += elapsed20_21.count();
        time_buf[1] += elapsed21_22.count();
        #endif
      */

    } // ip loop

    /*
      #ifdef WITH_CHRONO
      auto t2 = std::chrono::high_resolution_clock::now();
      #endif
    */

    double dx2_max = 0.0;
    double dx2_sum = 0.0;
    for (int ip = 0; ip < np; ++ip) {
      const double dx2 = DOT_PRODUCT(dxp[ip],dxp[ip])/(deltap[ip]*deltap[ip]);
      dx2_max = max(dx2_max,dx2);
      dx2_sum += dx2;
    }
    cout << "iter: " << iter << " counts: " << my_count[0] << " " << my_count[1] << " " << my_count[2] << " l2, linf: " << sqrt(dx2_sum/double(np)) << " " << sqrt(dx2_max) << endl;

    if ((done == 0)&&(my_count[0]+my_count[1]==0)&&((iter >= maxiter)||(sqrt(dx2_sum/double(np)) < zero))) {
      // we need to go through again and compute the nbrs...
      done = 1;
    }
    else if (done == 1) {
      // we are done...
      done = 2;
    }

    // only modify xp by dxp when done == 0. If done has been switched
    // to 1, we leave xp unmodified and recompute, storing nbr csr structures...

    if (done == 0) {
      for (int ip = 0; ip < np; ++ip)  {
        stoxp[ip] = stoxp_new[ip];
        FOR_I3 xp[ip][i] += dxp[ip][i];
      }
    }

    /*
      #ifdef WITH_CHRONO
      auto t3 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed0_1 = t1-t0;
      std::chrono::duration<double> elapsed1_2 = t2-t1;
      std::chrono::duration<double> elapsed2_3 = t3-t2;
      time_buf[2] += elapsed0_1.count();
      time_buf[3] += elapsed1_2.count();
      time_buf[4] += elapsed2_3.count();
      cout << "timing summary" << endl;
      cout << " t01: " << time_buf[2] << endl;
      cout << " t12: " << time_buf[3] << endl;
      cout << "  > t12 walking: " << time_buf[0] << endl;
      cout << "  > t12 cutting: " << time_buf[1] << endl;
      cout << " t23: " << time_buf[4] << endl;
      #endif
    */

  }

  delete[] st_flag;
  delete[] dxp;
  delete[] stoxp_new;
  delete[] xpost_i;
  delete[] xpost_v;

  // build the triangulation...

  for (int ip0 = 0; ip0 < np; ++ip0) {
    for (int nox0 = nboxp_i[ip0]; nox0 != nboxp_i[ip0+1]; ++nox0) {
      const int ip1 = nboxp_v[nox0];
      if (ip1 > ip0) {
        int nox1;
        for (nox1 = nboxp_i[ip1]; nox1 != nboxp_i[ip1+1]; ++nox1) {
          if (nboxp_v[nox1] == ip0)
            break;
        }
        if (nox1 != nboxp_i[ip1+1]) {
          for (nox1 = nboxp_i[ip1]; nox1 != nboxp_i[ip1+1]; ++nox1) {
            const int ip2 = nboxp_v[nox1];
            if (ip2 > ip1) {
              // ip0 has ip1,
              // ip1 has ip0 and ip2
              // first, confirm that ip0 has ip2...
              int nox0_;
              for (nox0_ = nboxp_i[ip0]; nox0_ != nboxp_i[ip0+1]; ++nox0_) {
                if (nboxp_v[nox0_] == ip2)
                  break;
              }
              if (nox0_ != nboxp_i[ip0+1]) {
                // finally, check for ip0 and ip1 in ip2's nbrs...
                int count = 0;
                for (int nox2 = nboxp_i[ip2]; nox2 != nboxp_i[ip2+1]; ++nox2) {
                  if ((nboxp_v[nox2] == ip0)||(nboxp_v[nox2] == ip1)) {
                    ++count;
                    if (count == 2) {
                      const double this_normal2[3] = TRI_NORMAL_2(xp[ip0],xp[ip1],xp[ip2]);
                      const int ist0 = stoxp[ip0];
                      const double st0_normal2[3] = TRI_NORMAL_2(xsp[spost[ist0][0]],xsp[spost[ist0][1]],xsp[spost[ist0][2]]);
                      const double st0_normal2_mag = MAG(st0_normal2); assert(st0_normal2_mag > 0.0);
                      const int ist1 = stoxp[ip1];
                      const double st1_normal2[3] = TRI_NORMAL_2(xsp[spost[ist1][0]],xsp[spost[ist1][1]],xsp[spost[ist1][2]]);
                      const double st1_normal2_mag = MAG(st1_normal2); assert(st1_normal2_mag > 0.0);
                      const int ist2 = stoxp[ip2];
                      const double st2_normal2[3] = TRI_NORMAL_2(xsp[spost[ist2][0]],xsp[spost[ist2][1]],xsp[spost[ist2][2]]);
                      const double st2_normal2_mag = MAG(st2_normal2); assert(st2_normal2_mag > 0.0);
                      // just sum the dot products here (should we normalize? sts can be any size)...
                      const double dp =
                        DOT_PRODUCT(this_normal2,st0_normal2)/st0_normal2_mag +
                        DOT_PRODUCT(this_normal2,st1_normal2)/st1_normal2_mag +
                        DOT_PRODUCT(this_normal2,st2_normal2)/st2_normal2_mag;
                      if (dp > 0.0) {
                        triVec.push_back(NewTri(ip0,ip1,ip2));
                      }
                      else {
                        triVec.push_back(NewTri(ip0,ip2,ip1));
                      }
                      break;
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

  delete[] nboxp_i;

}

void SimpleSurface::retri() {

  // hard-coded hack for a particular case for now...

  const int isz = 57;
  const double delta = 0.0075;
  const double area_coeff = 0.75;

  // st_flag holds the index in the local list of tris, or -1...
  int nst_local = 0;
  FOR_IST {
    if (szost[ist] == isz) {
      ++nst_local;
    }
  }

  int * stost_local = new int[nst_local];
  st_flag.resize(nst);
  st_flag.setAll(-1);
  nst_local = 0;
  FOR_IST {
    if (szost[ist] == isz) {
      stost_local[nst_local] = ist;
      st_flag[ist] = nst_local++;
    }
  }

  // this is the base condition for retriangulation. st_flag holds the index
  // into the nst_local tris or -1, and stost_local is the reverse...

  // ==========================================================

  // build a local version of the current patch...

  int (*spost_local)[3] = new int[nst_local][3];
  sp_flag.resize(nsp);
  sp_flag.setAll(-1);
  int nsp_local = 0;
  for (int ist_local = 0; ist_local < nst_local; ++ist_local) {
    const int ist = stost_local[ist_local];
    assert(st_flag[ist] == ist_local);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == -1)
        sp_flag[isp] = nsp_local++;
      spost_local[ist_local][i] = sp_flag[isp];
    }
  }
  double (*xsp_local)[3] = new double[nsp_local][3];
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] >= 0) {
      const int isp_local = sp_flag[isp];
      FOR_I3 xsp_local[isp_local][i] = xsp[isp][i];
    }
  }

  /*
    GeomUtils::writeTecplot("patch.dat",spost_local,nst_local,xsp_local);
    cout << "take a look" << endl;
    getchar();
  */

  // use the area to produce the new xsp counts and what tris they sit on...

  int * st_flag_local = new int[nst_local];
  int nsp_new = 0;
  double area_local = 0.0;
  for (int ist_local = 0; ist_local < nst_local; ++ist_local) {
    const double normal2[3] = TRI_NORMAL_2(xsp_local[spost_local[ist_local][0]],
                                           xsp_local[spost_local[ist_local][1]],
                                           xsp_local[spost_local[ist_local][2]]);
    const double area = 0.5*MAG(normal2);
    area_local += area;
    const int nsp_new_next = int(area_local/(area_coeff*delta*delta));
    st_flag_local[ist_local] = nsp_new_next - nsp_new;
    nsp_new = nsp_new_next;
  }

  cout << "area_local: " << area_local << " nsp_new: " << nsp_new << endl;

  // and now seed...
  double (*xsp_new)[3] = new double[nsp_new][3];
  double *delta_sp_new = new double[nsp_new];
  int * stosp_new = new int[nsp_new]; // holds the ist_local associated with the new point
  int isp_new = 0;
  for (int ist_local = 0; ist_local < nst_local; ++ist_local) {
    if (st_flag_local[ist_local] > 0) {
      for (int ii = 0; ii < st_flag_local[ist_local]; ++ii) {
        // uniform distribution on a tri...
        const double r0 = 0.999*double(rand())/double(RAND_MAX)+0.0005;
        const double r1 = 0.999*double(rand())/double(RAND_MAX)+0.0005;
        FOR_I3 xsp_new[isp_new][i] =
          (1.0-sqrt(r0))*xsp_local[spost_local[ist_local][0]][i] +
          sqrt(r0)*(1.0-r1)*xsp_local[spost_local[ist_local][1]][i] +
          sqrt(r0)*r1*xsp_local[spost_local[ist_local][2]][i];
        delta_sp_new[isp_new] = delta*1.25;
        stosp_new[isp_new] = ist_local;
        ++isp_new;
      }
    }
  }
  assert(isp_new == nsp_new);

  // and build a local nbost teost-like structure for our 3 tri nbrs...

  int (*nbost_local)[3] = new int[nst_local][3];
  map<const pair<int,int>,uint> edgeMap;
  for (int ist_local = 0; ist_local < nst_local; ++ist_local) {
    // loop on edges...
    FOR_I3 {
      const int isp0_local = spost_local[ist_local][i];
      const int isp1_local = spost_local[ist_local][(i+1)%3];
      map<const pair<int,int>,uint>::iterator iter = edgeMap.find(pair<int,int>(isp1_local,isp0_local));
      if (iter == edgeMap.end()) {
        // no edge found. Set stost to -1...
        nbost_local[ist_local][i] = -1;
        edgeMap[pair<int,int>(isp0_local,isp1_local)] = (uint(ist_local)<<2)|uint(i); // mix the tri index and edge index like teost
      }
      else {
        const int i_nbr = iter->second&3;
        const int ist_local_nbr = iter->second>>2;
        assert(nbost_local[ist_local_nbr][i_nbr] == -1);
        nbost_local[ist_local_nbr][i_nbr] = ist_local;
        nbost_local[ist_local][i] = ist_local_nbr;
        edgeMap.erase(iter);
      }
    }
  }

  // lloyd iterate the points...

  vector<NewTri> triVec;
  lloydIterate2d(triVec,xsp_new,delta_sp_new,stosp_new,nsp_new,spost_local,nbost_local,nst_local,xsp_local,nsp_local);

  int nst_new = triVec.size();
  int (*spost_new)[3] = new int[nst_new][3];
  for (int ist_new = 0; ist_new < nst_new; ++ist_new) {
    FOR_I3 spost_new[ist_new][i] = triVec[ist_new].spost[i];
  }
  triVec.clear();

  //GeomUtils::writePtsTecplot("xsp_new.dat",xsp_new,nsp_new);
  GeomUtils::writeTecplot("new_tris.dat",spost_new,nst_new,xsp_new);

  cout << "check out tecplot file new_tri.dat" << endl;


  delete[] spost_new;
  delete[] nbost_local;
  delete[] xsp_new;
  delete[] delta_sp_new;
  delete[] stosp_new;

  delete[] stost_local;
  delete[] spost_local;
  delete[] xsp_local;
  delete[] st_flag_local;

}

void SimpleSurface::lloydIterateInternalVorPoints(double (*xp)[3],int * stoxp, double * deltap,const int np_fixed,const int np_local,int * nboxp_i,vector<int>& nboxp_v,const int (*stost)[3],const int maxiter,const double delta,const double zero,const int igr,const bool b_power_diagram,const double growth_factor,const double growth_power) const {
  int * my_st_flag = new int[nst];

  double (*dxp)[3] = new double[np_local][3]; // lloyd motion
  int *stoxp_new = new int[np_local];
  int * xpost_i = new int[nst+1];
  int * xpost_v = new int[np_local];
  int * xp_flag = new int[np_local];

  CuttableVoronoiData cvd;
  stack<int> stack;
  vector<pair<double,int> > nbrVec;
  map<const pair<int,int>,int> edgeMap;
  map<const int,int> nodeMap;
  int noost[3],edost[3];
  map<const int,double*> x0Map;
  set<int> nbrSet;

  double *wgt = new double[np_local];
  for (int ip = 0; ip < np_local; ++ip)
    wgt[ip] = 1.0;

  // sp_buf is holding distance away from zone boundary...
  if (b_power_diagram) {
    assert(sp_buf != NULL);
    double sp_max = 0.0; // sp_min is by definition 0
    for (int ip = 0; ip < np_local; ++ip) {
      const int ist = stoxp[ip];
      FOR_I3 {
        const int isp = spost[ist][i];
        if (sp_buf[isp] > 0.0) {
          sp_max = max(sp_max,sp_buf[isp]);
          sp_buf[isp] *= -1.0; // flag as visited
        }
      }
    }
    for (int ip = 0; ip < np_local; ++ip) {
      const int ist = stoxp[ip];
      FOR_I3 {
        const int isp = spost[ist][i];
        if (sp_buf[isp] < 0.0) {
          sp_buf[isp] = -sp_buf[isp]/sp_max; // 0 -> 1
        }
        else if (sp_buf[isp] == 0.0) {
          sp_buf[isp] = 2.0*deltap[ip]/sp_max; // push points away from boundary
        }
      }
    }
  }

  int debug_iter = 0;
  int iter = 0;
  int done = 0;

  double * area2_ip = new double[np_local];
  for (int ip = 0; ip < np_local; ++ip) area2_ip[ip] = 0.0;

  while (done != 2) {
    ++iter;

    if (b_power_diagram) {
      for (int ip = 0; ip < np_local; ++ip) {
        const int ist = stoxp[ip];
        double xp_proj[3];
        MiscUtils::getClosestPointOnTriRobust(xp_proj,xp[ip],xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
        double n_sum = 0.0;
        wgt[ip] = 0.0;
        FOR_I3 {
          const double n[3] = TRI_NORMAL_2(xp_proj,xsp[spost[ist][(i+1)%3]],xsp[spost[ist][(i+2)%3]]);
          const double n_mag = MAG(n);
          wgt[ip] += n_mag*sp_buf[spost[ist][i]];
          n_sum += n_mag;
        }
        if (n_sum > 0.0)
          wgt[ip] /= n_sum;
      }
    }

    // build reverse-lookup xpost_i/v...
    for (int ist = 0; ist < nst; ++ist) xpost_i[ist+1] = 0;
    for (int ip = 0; ip < np_local; ++ip) {
      const int ist = stoxp[ip];
      ++xpost_i[ist+1];
    }
    xpost_i[0] = 0;
    for (int ist = 0; ist < nst; ++ist) xpost_i[ist+1] += xpost_i[ist];
    assert(xpost_i[nst] == np_local);
    for (int ip = 0; ip < np_local; ++ip) {
      const int ist = stoxp[ip];
      xpost_v[xpost_i[ist]++] = ip;
    }
    for (int ist = nst-1; ist > 0; --ist) xpost_i[ist] = xpost_i[ist-1];
    xpost_i[0] = 0;

    // check reverse look up...
    for (int ist = 0; ist < nst; ++ist) {
      for (int xos = xpost_i[ist]; xos != xpost_i[ist+1]; ++xos) {
        const int ip = xpost_v[xos];
        assert(stoxp[ip] == ist);
      }
    }

    // my_st_flag is used to build the surface patches. Start it off
    // with -1...
    for (int ist = 0; ist < nst; ++ist) my_st_flag[ist] = -1;

    // loop through our ip's...
    int my_count[3] = { 0, 0, 0 };
    for (int ip = 0; ip < np_local; ++ip) {
      // cout << "working on ip: " << ip << " stoxp: " << stoxp[ip] << endl;
      ++debug_iter;

      assert(cvd.nno == 0);
      assert(cvd.ned == 0);
      assert(nbrVec.empty());
      assert(edgeMap.empty());
      assert(nodeMap.empty());
      assert(my_st_flag[stoxp[ip]] != ip);
      my_st_flag[stoxp[ip]] = ip;
      assert(stack.empty());
      stack.push(-stoxp[ip]-1);  // -1 indexed
      while (!stack.empty()) {
        int ist = stack.top(); stack.pop();
        if (ist < 0) {
          // -----------------------------------------------------------------
          // a negative ist means both tris AND nbrs, so add the tris here...
          // -----------------------------------------------------------------
          ist = -ist-1;
          // use edge-matching to set as much of the noost and edost that we can...
          FOR_I3 noost[i] = -1;
          FOR_I3 {
            edost[i] = -1;
            // here the edgeMap is used to find edges in terms of surface node pairs...
            map<const pair<int,int>,int>::iterator it = edgeMap.find(pair<int,int>(spost[ist][(i+1)%3],spost[ist][i]));
            if (it != edgeMap.end()) {
              edost[i] = it->second;
              edgeMap.erase(it);
              if (noost[i] == -1) noost[i] = cvd.nooed[edost[i]][1];
              else assert(noost[i] == cvd.nooed[edost[i]][1]);

              if (noost[(i+1)%3] == -1) noost[(i+1)%3] = cvd.nooed[edost[i]][0];
              else assert(noost[(i+1)%3] == cvd.nooed[edost[i]][0]);
            }
          }
          // step 3: now noost and edost are set from edges where possible. At this point, complete any nodes not set...
          FOR_I3 {
            if (noost[i] == -1) {
              map<const int,int>::iterator it = nodeMap.find(spost[ist][i]);
              if (it == nodeMap.end()) {
                const int ino = cvd.new_node();
                noost[i] = ino;
                nodeMap[spost[ist][i]] = ino;
                FOR_J3 cvd.x_no[ino][j] = xsp[spost[ist][i]][j] - xp[ip][j];
              }
              else {
                noost[i] = it->second;
              }
            }
          }
          // step 4: finally complete the edges...
          FOR_I3 {
            if (edost[i] == -1) {
              const int ied = cvd.new_edge();
              edgeMap[pair<int,int>(spost[ist][i],spost[ist][(i+1)%3])] = ied;
              cvd.nooed[ied][0] = noost[i];
              cvd.nooed[ied][1] = noost[(i+1)%3];
              cvd.faoed[ied][0] = ist;
              const int ist_nbr = stost[ist][i];
              if ((ist_nbr >= 0)&&(st_flag[ist_nbr] == igr)) cvd.faoed[ied][1] = -1;  // check to see if valid edge nbr based on current group
              else cvd.faoed[ied][1] = -2; // used to mark external boundary
            }
            else {
              const int ied = edost[i];
              assert(cvd.nooed[ied][1] == noost[i]);
              assert(cvd.nooed[ied][0] == noost[(i+1)%3]);
              assert(cvd.faoed[ied][0] != ist); // no edge can have the same tri on both sides
              assert(cvd.faoed[ied][1] == -1);
              cvd.faoed[ied][1] = ist;
            }
          }
        }
        // -----------------------------------------------------------------
        // always add nbrs...
        // -----------------------------------------------------------------
        assert(my_st_flag[ist] == ip);
        //cout << "adding nbrs from tri: " << ist << endl;
        for (int xos = xpost_i[ist]; xos != xpost_i[ist+1]; ++xos) {
          const int ip_nbr = xpost_v[xos];
          if (ip_nbr != ip) {
            const double d2 = DIST2(xp[ip_nbr],xp[ip]);
            if (d2 == 0.0) {
              cout << "got zero: " << ip << " " << ip_nbr << endl;
              FILE * fp = fopen("xp.dat","w");
              fprintf(fp,"VARIABLES=\"X\" \"Y\" \"Z\"\n");
              fprintf(fp,"%18.15e %18.15e %18.15e\n",xp[ip][0],xp[ip][1],xp[ip][2]);
              fclose(fp);
            }
            assert(d2 > 0.0); // possible to seed 2 points at the same location?
            if (d2 < deltap[ip]*deltap[ip]) nbrVec.push_back(pair<double,int>(d2,ip_nbr));
          }
        }
        // if the neighboring tri is closer than 0.5*deltap[ip], then we need to add it AND
        // its nbrs. If it is closer than deltap[ip], then we just need to add the nbrs.
        FOR_I3 {
          const int ist_nbr = stost[ist][i];
          assert(ist_nbr != ist);
          if ((ist_nbr >= 0)&&(st_flag[ist_nbr] == igr)&&(my_st_flag[ist_nbr] != ip)) {
            my_st_flag[ist_nbr] = ip; // do not ever revisit this tri
            // the 3 nodes of ist_nbr...
            const int isp0 = spost[ist_nbr][0];
            const int isp1 = spost[ist_nbr][1];
            const int isp2 = spost[ist_nbr][2];
            const double d2 = MiscUtils::getPointToTriDist2(xp[ip],xsp[isp0],xsp[isp1],xsp[isp2]);
            // if this tri is closer than 0.5*deltap[ip], we add it to the surface, and if less than
            // deltap[ip], we need to potentially add the nbrs (but not the surface)...
            if (d2 < 0.2501*deltap[ip]*deltap[ip])  {
              // both tri and nbrs...
              stack.push(-ist_nbr-1);
            }
            else if (d2 < deltap[ip]*deltap[ip])  {
              // just nbrs...
              stack.push(ist_nbr);
            }
          }
        }
      }
      nodeMap.clear();
      edgeMap.clear();

      cvd.setD2Max();
      sort(nbrVec.begin(),nbrVec.end());
      // we now have a sorted list of nbrs. Cut the cvd against the nbrs until it cannot possibly
      // be cut anymore...
      for (int in = 0,nnbr=nbrVec.size(); in < nnbr; ++in) {
        // recall that the .first contains the d2 nbr distance. and cvd.d2_max contains
        // the maximum node radius. as soon as d2 >= 4*cvd.d2_max, there cannot be any
        // more nbrs that cut this cv...
        if (nbrVec[in].first > 4.0*cvd.d2_max) break; // done: all remaining nbrs are too far away to possibly cut the current vd...
        const int ip_nbr = nbrVec[in].second;
        double dn[3]; FOR_I3 dn[i] = 0.5*(xp[ip_nbr][i] - xp[ip][i]);
        if ((2.0*dn[0] > cvd.Lmax[0])||
            (2.0*dn[1] > cvd.Lmax[1])||
            (2.0*dn[2] > cvd.Lmax[2])||
            (2.0*dn[0] < cvd.Lmin[0])||
            (2.0*dn[1] < cvd.Lmin[1])||
            (2.0*dn[2] < cvd.Lmin[2]) ) {
          continue;
        }
        // only cut if we are inside the 6 paraboloids...
        // if we got here, we are cutting...
        cvd.cut_surf(dn,-8-ip_nbr); // note that local nbrs use a -8 convention. cvd.d2_max updated automatically
      }
      nbrVec.clear();

      // compute the centroid of the remaining parts...
      bool has_seed_boundary = false;
      bool has_boundary = false;
      double dxc[3] = { 0.0, 0.0, 0.0 };
      double area2_sum = 0.0;
      assert(x0Map.empty());
      for (int ied = 0; ied < cvd.ned; ++ied) {
        if (cvd.faoed[ied][0] >= 0) {
          map<const int,double*>::const_iterator it = x0Map.find(cvd.faoed[ied][0]);
          if (it == x0Map.end()) {
            // on the first time we visit any face, store a node to use in all
            // the tris with other edges (the face is planar and convex, so this is fine)...
            x0Map[cvd.faoed[ied][0]] = cvd.x_no[cvd.nooed[ied][0]];
          }
          else if (cvd.x_no[cvd.nooed[ied][1]] != it->second) {
            const double normal2[3] = TRI_NORMAL_2(it->second,cvd.x_no[cvd.nooed[ied][0]],cvd.x_no[cvd.nooed[ied][1]]);
            double area2 = MAG(normal2);
            if (cvd.faoed[ied][1] <= -8) {
              const int ip_nbr = -cvd.faoed[ied][1]-8;
              area2 *= growth_factor*pow(wgt[ip]/wgt[ip_nbr],growth_power);
            }
            FOR_I3 dxc[i] += area2*(it->second[i]+cvd.x_no[cvd.nooed[ied][0]][i]+cvd.x_no[cvd.nooed[ied][1]][i]);
            area2_sum += area2;
          }
        }
        if (cvd.faoed[ied][1] >= 0) {
          map<const int,double*>::const_iterator it = x0Map.find(cvd.faoed[ied][1]);
          if (it == x0Map.end()) {
            x0Map[cvd.faoed[ied][1]] = cvd.x_no[cvd.nooed[ied][1]];
          }
          else if (cvd.x_no[cvd.nooed[ied][0]] != it->second) {
            const double normal2[3] = TRI_NORMAL_2(it->second,cvd.x_no[cvd.nooed[ied][1]],cvd.x_no[cvd.nooed[ied][0]]);
            double area2 = MAG(normal2);
            if (cvd.faoed[ied][0] <= -8) {
              const int ip_nbr = -cvd.faoed[ied][0]-8;
              area2 *= growth_factor*pow(wgt[ip]/wgt[ip_nbr],growth_power);
            }
            FOR_I3 dxc[i] += area2*(it->second[i]+cvd.x_no[cvd.nooed[ied][0]][i]+cvd.x_no[cvd.nooed[ied][1]][i]);
            area2_sum += area2;
          }
        }
        else if (cvd.faoed[ied][1] == -1) {
          has_seed_boundary = true;
          break;
        }
        else {
          // must be a -2: this is ok (flagged region boundary)...
          assert(cvd.faoed[ied][1] == -2);
          has_boundary = true;
        }
      }

      area2_ip[ip] = area2_sum;

      if (has_seed_boundary) {
        // we are too small. Increase deltap[ip], and do NOT apply any motion to the point. The remaining
        // patch is arbitrarily large...
        ++my_count[0]; xp_flag[ip] = 0;
        assert(done == 0);
        deltap[ip] *= 1.5;
        stoxp_new[ip] = stoxp[ip];
        FOR_I3 dxp[ip][i] = 0.0;
      }
      else {
        // Even though we have cut away all the original seed boundary of this patch, there still
        // may be nbrs out there that can cut some of our corners off if deltap was not large
        // enough...
        const double delta_half = sqrt(cvd.d2_max); // cvd.d2_max stores the furthest distance of any node
        if (deltap[ip] <= 2.0*delta_half) {
          // other neighbors might exist that could still cut this volume. Find them in
          // the next iteration...
          ++my_count[1]; xp_flag[ip] = 1;
          assert(done == 0);
          deltap[ip] = 2.5*delta_half; // could be just 2.0 here? - use a little larger so face checks pass for sure
          stoxp_new[ip] = stoxp[ip];
          FOR_I3 dxp[ip][i] = 0.0;
        }
        else {
          ++my_count[2]; xp_flag[ip] = 2;
          // this is good. Do some more work to decide where to move this point...

          if (done == 0) {
            // regular lloyd iteration pass
            assert(area2_sum > 0.0);
            FOR_I3 dxc[i] /= area2_sum*3.0;
            deltap[ip] = 2.5*delta_half;
            // in addition to lloyd, add random perturbation...
            //FOR_I3 dxc[i] += (double(rand())/double(RAND_MAX)-0.5)*deltap[ip]*factor;
            // now we need to figure out which patch we are closest to...
            int ist_closest = -1;
            double d2_closest;
            for (int ied = 0; ied < cvd.ned; ++ied) {
              if (cvd.faoed[ied][0] >= 0) {
                map<const int,double*>::const_iterator it = x0Map.find(cvd.faoed[ied][0]);
                assert(it != x0Map.end());
                if ((cvd.x_no[cvd.nooed[ied][0]] != it->second)&&(cvd.x_no[cvd.nooed[ied][1]] != it->second)) {
                  double this_dxp[3];
                  MiscUtils::getClosestPointOnTriRobust(this_dxp,dxc,it->second,cvd.x_no[cvd.nooed[ied][0]],cvd.x_no[cvd.nooed[ied][1]]);
                  const double this_d2 = DIST2(this_dxp,dxc);
                  if ((ist_closest == -1)||(this_d2 < d2_closest)) {
                    FOR_I3 dxp[ip][i] = this_dxp[i];
                    d2_closest = this_d2;
                    ist_closest = cvd.faoed[ied][0];
                  }
                }
              }
              if (cvd.faoed[ied][1] >= 0) {
                map<const int,double*>::const_iterator it = x0Map.find(cvd.faoed[ied][1]);
                assert(it != x0Map.end());
                if ((cvd.x_no[cvd.nooed[ied][0]] != it->second)&&(cvd.x_no[cvd.nooed[ied][1]] != it->second)) {
                  double this_dxp[3];
                  MiscUtils::getClosestPointOnTriRobust(this_dxp,dxc,it->second,cvd.x_no[cvd.nooed[ied][0]],cvd.x_no[cvd.nooed[ied][1]]);
                  const double this_d2 = DIST2(this_dxp,dxc);
                  if ((ist_closest == -1)||(this_d2 < d2_closest)) {
                    FOR_I3 dxp[ip][i] = this_dxp[i];
                    d2_closest = this_d2;
                    ist_closest = cvd.faoed[ied][1];
                  }
                }
              }
            }
            // make sure we found one!...
            assert(ist_closest != -1);
            stoxp_new[ip] = ist_closest;
          }
          else {
            assert(done == 1);
            // this is a recompute of the last time...
            if (has_boundary) xp_flag[ip] = 1;
            assert(nbrSet.empty());
            for (int ied = 0; ied < cvd.ned; ++ied) {
              if (cvd.faoed[ied][0] < 0) {
                assert(cvd.faoed[ied][0] <= -8);
                nbrSet.insert(-cvd.faoed[ied][0]-8);
              }
            }
            assert(!nbrSet.empty());
            for (set<int>::const_iterator it = nbrSet.begin(); it != nbrSet.end(); ++it) nboxp_v.push_back(*it);
            nboxp_i[ip+1] = nboxp_v.size();
            nbrSet.clear();
          }
        }
      }
      // clear cvd for next time...
      x0Map.clear();
      cvd.clear();
    } // ip loop

    double dx2_max = 0.0;
    double dx2_sum = 0.0;
    for (int ip = np_fixed; ip < np_local; ++ip) {
      const double dx2 = DOT_PRODUCT(dxp[ip],dxp[ip]);
      dx2_max = max(dx2_max,dx2);
      dx2_sum += dx2;
    }

    if (false) {
      FILE * fp;
      char filename[128];
      sprintf(filename,"xp.%03d.dat",igr);
      if (iter == 1) {
        fp = fopen(filename,"w");
        fprintf(fp,"VARIABLES=\"X\" \"Y\" \"Z\" \"area\"\n");
      }
      else {
        fp = fopen(filename,"a");
      }
      fprintf(fp,"ZONE T=\"iteration_%d\"\n",iter-1);
      double area_range[3] = {HUGE_VAL,-HUGE_VAL,0.0};
      for (int ip = np_fixed; ip < np_local; ++ip) {
        area_range[0] = min(area_range[0],area2_ip[ip]);
        area_range[1] = max(area_range[1],area2_ip[ip]);
        area_range[2] += area2_ip[ip];
      }
      const double area_mean = area_range[2]*0.5/(np_local-np_fixed);
      // const double one_o_delta_area = 1.0 / (M_PI*delta*delta*0.25);  // circle
      const double one_o_delta_area = 1.0 / (6.0*tan(M_PI/6.0)*delta*delta*0.25);  // hexagon
      COUT2(" > iter[" << iter << "] area min,mean,max: " << area_range[0]*0.5 << " " << area_mean << " " << 0.5*area_range[1]);

      for (int ip = np_fixed; ip < np_local; ++ip) {
        fprintf(fp,"%18.15e %18.15e %18.15e %18.15e\n",xp[ip][0],xp[ip][1],xp[ip][2],0.5*area2_ip[ip]*one_o_delta_area);
      }
      fclose(fp);
    }

    cout << "iter: " << iter << " counts: " << my_count[0] << " (has seed), " << my_count[1] << " (grow delta), " << my_count[2] << " (passed); l2,linf: " << sqrt(dx2_sum/double(np_local))/delta << " " << sqrt(dx2_max)/delta << endl;

    if ((done == 0)&&(my_count[0]+my_count[1]==0)&&((iter >= maxiter)||(sqrt(dx2_sum/double(np_local))/delta < zero))) {
      // we need to go through again and compute the nbrs...
      done = 1;
    }
    else if (done == 1) {
      // we are done...
      done = 2;
    }

    // only modify xp by dxp when done == 0. If done has been switched
    // to 1, we leave xp unmodified and recompute, storing nbr csr structures...

    if (done == 0) {
      for (int ip = np_fixed; ip < np_local; ++ip)  {  // only update non-fixed nodes
        stoxp[ip] = stoxp_new[ip];
        FOR_I3 xp[ip][i] += dxp[ip][i];
      }
    }
  }

  DELETE(area2_ip);

  DELETE(xp_flag);
  DELETE(my_st_flag);
  DELETE(dxp);
  DELETE(stoxp_new);
  DELETE(xpost_i);
  DELETE(xpost_v);
  DELETE(wgt);
}
