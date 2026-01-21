#include "SimpleSurface.hpp"
#include "Histogram.hpp"
#define MESH_INSPECTOR_IST 100100100
#include "../stitch/CuttableVoronoiData.hpp"

void SimpleSurface::testSeed(Param * param,const bool help) {

  cout << "SimpleSurface::testSeed()" << endl;

  ensureTeost();
  Histogram h(-1.0,1.0,MPI_COMM_NULL);

  st_flag.setLength(nst);
  st_flag.setAll(0);

  for (int ist = 0; ist < nst; ++ist) {
    const double n_st[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    const double mag_n_st = MAG(n_st);
    if (mag_n_st > 0.0) {
      FOR_I3 {
        int ist_nbr;
        if (getAlignedTriNbr(ist_nbr,ist,i)) {
          const double n_st_nbr[3] = TRI_NORMAL_2(xsp[spost[ist_nbr][0]],xsp[spost[ist_nbr][1]],xsp[spost[ist_nbr][2]]);
          const double mag_n_st_nbr = MAG(n_st_nbr);
          if (mag_n_st_nbr > 0.0) {
            const double dp = DOT_PRODUCT(n_st,n_st_nbr)/mag_n_st/mag_n_st_nbr;
            h.add(&dp,1);
            if (dp < -0.75) {
              st_flag[ist] = 1;
              st_flag[ist_nbr] = 1;
            }
          }
        }
      }
    }
  }

  h.write("hist.dat");

  // now build a table of these tris...

  int nst_problem = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] == 1)
    ++nst_problem;
  }
  cout << " > nst_problem: " << nst_problem << endl;

  int * st_problem = new int[nst_problem];
  int ist_problem = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] == 1)
    st_problem[ist_problem++] = ist;
  }
  assert(ist_problem == nst_problem);

  // set default delta values: can override...

  double delta_min = 1.0E-6*getRmax();
  double delta_max = 1.0E-3*getRmax();

  cout << " > default delta_min: " << delta_min << " delta_max: " << delta_max << endl;

  int iarg = 0;
  while (iarg < param->size()) {
    const string token = param->getString(iarg++);
    if (token == "DELTA_MIN") {
      delta_min = param->getDouble(iarg++);
      cout << " > delta_min: " << delta_min << endl;
    }
    else if (token == "DELTA_MAX") {
      delta_max = param->getDouble(iarg++);
      cout << " > delta_max: " << delta_max << endl;
    }
    else {
      CERR("unrecognized TEST_SEED token: " << token);
    }
  }

  // step 1: build stAdt...

  cout << " > building stAdt..." << endl;

  double (*bbmin)[3] = new double[nst][3];
  double (*bbmax)[3] = new double[nst][3];
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 bbmin[ist][i] = HUGE_VAL;
    FOR_I3 bbmax[ist][i] = -HUGE_VAL;
    FOR_I3 {
      const int isp = spost[ist][i];
      // orig
      FOR_J3 bbmin[ist][j] = min(bbmin[ist][j],xsp[isp][j]);
      FOR_J3 bbmax[ist][j] = max(bbmax[ist][j],xsp[isp][j]);
    }
  }

  // build adt...
  Adt<double> * stAdt = new Adt<double>(nst,bbmin,bbmax);
  delete[] bbmin;
  delete[] bbmax;

  // now loop forever...

  cout << " > starting loop..." << endl;

  CuttableVoronoiData cvd;
  map<const int,int> nodeMap;
  map<const pair<int,int>,int> edgeMap;
  int noost[3];
  vector<int> stVec; // surface-tri-vector

  int iter = 0;
  while (1) {

    ++iter;
    if (iter%1000 == 0)
    cout << " > iter: " << iter << endl;

    cvd.clear();
    nodeMap.clear();
    edgeMap.clear();
    assert(cvd.nno == 0);
    assert(cvd.ned == 0);

    // decide where to build the seed. Randomly choose a surface tri, and
    // build inside the tri at some random location...

    const int seed = rand()%nst_problem;
    const int ist_seed = st_problem[seed];

    const double n_st[3] = TRI_NORMAL_2(xsp[spost[ist_seed][0]],xsp[spost[ist_seed][1]],xsp[spost[ist_seed][2]]);
    const double n_mag = MAG(n_st);
    if (n_mag > 0.0) {

      // choose a delta...
      const double delta = delta_min + double(rand())/double(RAND_MAX)*(delta_max-delta_min);

      // choose 2 baricentric randoms...
      const double w0 = double(rand())/double(RAND_MAX);
      const double w1 = double(rand())/double(RAND_MAX);
      const double w2 = double(rand())/double(RAND_MAX);
      if (w0+w1+w2 > 0.0) {

        double xp[3];
        FOR_I3 xp[i] = ( w0*xsp[spost[ist_seed][0]][i] + w1*xsp[spost[ist_seed][1]][i] + w2*xsp[spost[ist_seed][2]][i] )/(w0+w1+w2);

        // recall that the seed is built with an L = 0.505*delta (don't ask why)...

        // perturb the location of the seed off the boundary in the inward direction by some fraction
        // of 0.505*delta...

        const double r = double(rand())/double(RAND_MAX);
        FOR_I3  xp[i] -= r*n_st[i]/n_mag*0.5*delta;

        // now grab the surfaces that intersect the sphere...
        // note: this borrows from the code in VoronoiPart::addNearbySurfaceTrisAndTransforms...

        const double bbmin[3] = {
          xp[0]-0.51*delta,
          xp[1]-0.51*delta,
          xp[2]-0.51*delta
        };
        const double bbmax[3] = {
          xp[0]+0.51*delta,
          xp[1]+0.51*delta,
          xp[2]+0.51*delta
        };

        // grab the nearby tris from the subSurface. Recall subSurface is the part of the
        // surface that may intersect the local points being rebuilt. It should be guaranteed to
        // service this request...

        assert(stVec.empty());
        stAdt->buildListForBBox(stVec,bbmin,bbmax);

        // add to cvd...

        for (int ii = 0, ii_end=stVec.size(); ii < ii_end; ++ii) {
          // the tri is uniquely associated with its combination of its "ist" and its bits if any...
          const int ist  = stVec[ii];
          // everyone has access to the full surface, so grab and assemble tris directly...
          // nodes first...
          FOR_I3 {
            const int isp = spost[ist][i];
            map<const int,int>::iterator iter = nodeMap.find(isp);
            if (iter == nodeMap.end()) {
              const int ino = cvd.new_node();
              nodeMap[isp] = ino;
              FOR_J3 cvd.x_no[ino][j] = xsp[isp][j] - xp[j];
              noost[i] = ino; // store the nodes in noost[3]. They are req'd to connect edges
            }
            else {
              noost[i] = iter->second;
            }
          }
          // edges...
          FOR_I3 {
            const int ino0 = noost[i];
            const int ino1 = noost[(i+1)%3];
            map<const pair<int,int>,int>::iterator iter = edgeMap.find(pair<int,int>(ino1,ino0));
            if (iter == edgeMap.end()) {
              const int ied = cvd.new_edge();
              edgeMap[pair<int,int>(ino0,ino1)] = ied;
              cvd.nooed[ied][0] = ino0;
              cvd.nooed[ied][1] = ino1;
              cvd.faoed[ied][0] = ist;
              cvd.faoed[ied][1] = -1;
            }
            else {
              const int ied = iter->second;
              edgeMap.erase(iter);
              assert(cvd.faoed[ied][0] != ist); // no edge can have the same tri on both sides
              assert(cvd.faoed[ied][1] == -1);
              cvd.faoed[ied][1] = ist;
            }
          }
        }

        // cut surface with the box...

        // ================================================
        // now cut surface patch wrt corners...
        // this introduces new edges around the surface boundary
        // that have their faoed[ied][1] linked to a surface face,
        // and their faoed[ied][0] linked to a cut surface indicated
        // by -2,-3,-4,-5,-6,-7 for -x,+x,-y,+y,etc...
        // ================================================
        if (cvd.ned > 0) {
          // -2 == x0
          const double n[3] = { -0.505*delta, 0.0, 0.0 };
          cvd.cut_surf(n,-2);
          if (cvd.ned > 0) {
            // -3 == x1
            const double n[3] = { 0.505*delta, 0.0, 0.0 };
            cvd.cut_surf(n,-3);
            if (cvd.ned > 0) {
              // -4 == y0
              const double n[3] = { 0.0, -0.505*delta, 0.0 };
              cvd.cut_surf(n,-4);
              if (cvd.ned > 0) {
                // -5 == y1
                const double n[3] = { 0.0, 0.505*delta, 0.0 };
                cvd.cut_surf(n,-5);
                if (cvd.ned > 0) {
                  // -6 == z0
                  const double n[3] = { 0.0, 0.0, -0.505*delta };
                  cvd.cut_surf(n,-6);
                  if (cvd.ned > 0) {
                    // -7 == z1
                    const double n[3] = { 0.0, 0.0, 0.505*delta };
                    cvd.cut_surf(n,-7);
                    if (cvd.ned > 0) {
                      // check: if we still have stuff left, all cvd.faoed[ied][1]
                      // that were -1 should have been cut off...
                      for (int ied = 0; ied < cvd.ned; ++ied) {
                        if (!(cvd.faoed[ied][1] != -1)) {
                          cvd.writeTecplot(0,xp);
                        }
                        assert(cvd.faoed[ied][1] != -1);
                      }
                    }
                  }
                }
              }
            }
          }
        }

        // if any surface was left,
        assert(cvd.ned > 0);

        // build the starting region considering the remaining surface patch...
        int ierr = cvd.completeCube(0.505*delta);
        if (ierr != 0) {
          cvd.writeTecplot(0,xp);

          IntFlag my_st_flag(nst); my_st_flag.setAll(0);
          for (int ii = 0,ii_end=stVec.size(); ii < ii_end; ++ii) {
            my_st_flag[stVec[ii]] = 1;
          }

          writeSelectedFacesByZoneToTecplot("seed_surface_patch.dat",my_st_flag);

          cout << "press enter to continue..." << endl;
          getchar();
        }

        stVec.clear();
      }

    }

  }

}

void SimpleSurface::testNonmanifoldSeeds(Adt<double> * stAdt,const int n_samples_max) {

  COUT2("SimpleSurface::testNonmanifoldSeeds()");

  ensureTeost();
  ensureNonManifoldData();

  stack<int> stStack;
  for (int ist = 0; ist < nst; ++ist) {
    nmTriData.unsetVisited(ist);
    if (nmTriData.isOverlapping(ist) || nmTriData.isImprinting(ist) || nmTriData.isIntersecting(ist)) stStack.push(ist);
  }
  int n_problem_tris = stStack.size();

  COUT1("    > additional filtering of " << n_problem_tris << " flagged tris:");
  const double start_time = MPI_Wtime();
  COUT2("       > max samples per candidate: " << n_samples_max);
  cout << "       > using Stitch cell build checks" << endl;
  cout << "         0%" << std::flush;
  const int interval = int(floor(0.1*(double)n_problem_tris));
  n_problem_tris = 0;
  int percent_done = 10;
  vector<int> neighborhood_tris;
  // only look at diagnostics returned candidates

  while (!stStack.empty()) {
    const int ist = stStack.top(); stStack.pop();
    assert(nmTriData.isOverlapping(ist) || nmTriData.isImprinting(ist) || nmTriData.isIntersecting(ist));
    if ((interval > 0) && (n_problem_tris > 0)) {
      if (n_problem_tris%interval == 0) {  // separate modulo b/c seems to produce problem if both are 0
        cout << ".." << percent_done << "\%" << std::flush;
        percent_done += 10;
      }
    }

    ++n_problem_tris;
    neighborhood_tris.clear();
    if ( isTriValidBasedOnSeeds(neighborhood_tris,ist,stAdt,n_samples_max) ) nmTriData.setClean(ist);

    if (!neighborhood_tris.empty()) {
      for (int ii=0,end=neighborhood_tris.size(); ii<end; ++ii) {
        const int iist = neighborhood_tris[ii];
        nmTriData.setIntersecting(iist);
      }
    }
  }

  // stack-based search over problem tri regions...?
  // vector<int> neighborhood_tris;
  // while (!stStack.empty()) {
  //   const int ist = stStack.top(); stStack.pop();
  //
  //   if (!nmTriData.isVisited(ist)) {
  //     neighborhood_tris.clear();
  //     if ( isTriValidBasedOnSeeds(neighborhood_tris,ist,stAdt,n_samples) ) {
  //       nmTriData.setClean(ist);
  //       nmTriData.setVisited(ist);
  //     }
  //     else {
  //       // invalid seed, insert all seed-cut tris to stack for vicinity checking
  //       for (int ii=0,end=neighborhood_tris.size(); ii<end; ++ ii) {
  //         const int iist = neighborhood_tris[ii];
  //         if (!nmTriData.isVisited(iist)) {
  //           stStack.push(iist);
  //           nmTriData.setIntersecting(iist);
  //         }
  //       }
  //     }
  //   }
  // }
  cout << "..100%" << endl;
  const double processing_time = MPI_Wtime()-start_time;
  COUT2("    > finished checking in: " << processing_time << "s");
}

bool SimpleSurface::isTriValidBasedOnSeeds(vector<int>& neighborhood_tris,const int ist, Adt<double> * stAdt, const int n_samples_max) {
  UNUSED(neighborhood_tris);
  // tri to evaulate is passed in

  Histogram h(-1.0,1.0,MPI_COMM_NULL);

  // now loop forever...
  CuttableVoronoiData cvd;
  map<const int,int> nodeMap;
  map<const pair<int,int>,int> edgeMap;
  int noost[3];

  const double n_st[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
  const double n_mag = MAG(n_st);
  if (n_mag <= 0.0) return false;

  // determine min/max lengthscales to test (min altitude and max edge length to start...)
  double delta_min = HUGE_VAL;
  double delta_max = -HUGE_VAL;

  double edge_length[3];
  edge_length[0] = DIST(xsp[spost[ist][0]],xsp[spost[ist][1]]);
  edge_length[1] = DIST(xsp[spost[ist][1]],xsp[spost[ist][2]]);
  edge_length[2] = DIST(xsp[spost[ist][2]],xsp[spost[ist][0]]);
  // FOR_I3 delta_max = max(delta_max,edge_length[i]);

  const double s = 0.5*(edge_length[0]+edge_length[1]+edge_length[2]);
  // const double altitude_numerator = 2*sqrt( s*(s-edge_length[0])*(s-edge_length[1])*(s-edge_length[2]) );
  // FOR_I3 delta_min = min(delta_min, (altitude_numerator/edge_length[i]) );

  // make sure much smaller than tri itself...use inradius as starting point
  const double n[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
  const double area = MAG(n);
  const double inradius = area/s;

  delta_max = 0.25*inradius;
  delta_min = 0.01*inradius;

  int n_samples = int(ceil(1.5*area/(delta_min*delta_min)));  // conservative? dunno
  n_samples = min(n_samples,n_samples_max);

  vector<int> stVec; // surface-tri-vector

  int iter = 0;
  // MPI_Pause("pausing");
  while (iter < n_samples) {

    ++iter;

    cvd.clear();
    nodeMap.clear();
    edgeMap.clear();
    assert(cvd.nno == 0);
    assert(cvd.ned == 0);

    // decide where to build the seed somewhat randomly
    // choose a delta...
    const double delta = delta_min + double(rand())/double(RAND_MAX)*(delta_max-delta_min);

    // choose 2 baricentric randoms...
    const double r1 = double(rand())/double(RAND_MAX);
    const double r2 = double(rand())/double(RAND_MAX);
    const double sqrt_r1 = sqrt(r1);
    if (true) {

      double xp[3];
      FOR_I3 xp[i] = ( (1.0-sqrt_r1)*xsp[spost[ist][0]][i] + sqrt_r1*(1.0-r2)*xsp[spost[ist][1]][i] + sqrt_r1*r2*xsp[spost[ist][2]][i] );

      // recall that the seed is built with an L = 0.505*delta (don't ask why)...

      // perturb the location of the seed off the boundary in the inward direction by some fraction
      // of 0.505*delta...

      double r = double(rand())/double(RAND_MAX);
      FOR_I3  xp[i] -= r*n_st[i]/n_mag*0.5*delta;

      // now grab the surfaces that intersect the sphere...
      // note: this borrows from the code in VoronoiPart::addNearbySurfaceTrisAndTransforms...

      const double bbmin[3] = {
        xp[0]-0.51*delta,
        xp[1]-0.51*delta,
        xp[2]-0.51*delta
      };
      const double bbmax[3] = {
        xp[0]+0.51*delta,
        xp[1]+0.51*delta,
        xp[2]+0.51*delta
      };

      // grab the nearby tris from the subSurface. Recall subSurface is the part of the
      // surface that may intersect the local points being rebuilt. It should be guaranteed to
      // service this request...

      assert(stVec.empty());
      stAdt->buildListForBBox(stVec,bbmin,bbmax);

      // add to cvd...

      for (int ii = 0,ii_end=stVec.size(); ii < ii_end; ++ii) {
        // the tri is uniquely associated with its combination of its "ist" and its bits if any...
        const int ist  = stVec[ii];
        // everyone has access to the full surface, so grab and assemble tris directly...
        // nodes first...
        FOR_I3 {
          const int isp = spost[ist][i];
          map<const int,int>::iterator iter = nodeMap.find(isp);
          if (iter == nodeMap.end()) {
            const int ino = cvd.new_node();
            nodeMap[isp] = ino;
            FOR_J3 cvd.x_no[ino][j] = xsp[isp][j] - xp[j];
            noost[i] = ino; // store the nodes in noost[3]. They are req'd to connect edges
          }
          else {
            noost[i] = iter->second;
          }
        }
        // edges...
        FOR_I3 {
          const int ino0 = noost[i];
          const int ino1 = noost[(i+1)%3];
          map<const pair<int,int>,int>::iterator iter = edgeMap.find(pair<int,int>(ino1,ino0));
          if (iter == edgeMap.end()) {
            const int ied = cvd.new_edge();
            edgeMap[pair<int,int>(ino0,ino1)] = ied;
            cvd.nooed[ied][0] = ino0;
            cvd.nooed[ied][1] = ino1;
            cvd.faoed[ied][0] = ist;
            cvd.faoed[ied][1] = -1;
          }
          else {
            const int ied = iter->second;
            edgeMap.erase(iter);
            assert(cvd.faoed[ied][0] != ist); // no edge can have the same tri on both sides
            assert(cvd.faoed[ied][1] == -1);
            cvd.faoed[ied][1] = ist;
          }
        }
      }

      // cut surface with the box...

      // ================================================
      // now cut surface patch wrt corners...
      // this introduces new edges around the surface boundary
      // that have their faoed[ied][1] linked to a surface face,
      // and their faoed[ied][0] linked to a cut surface indicated
      // by -2,-3,-4,-5,-6,-7 for -x,+x,-y,+y,etc...
      // ================================================
      if (cvd.ned > 0) {
        // -2 == x0
        const double n[3] = { -0.505*delta, 0.0, 0.0 };
        cvd.cut_surf(n,-2);
        if (cvd.ned > 0) {
          // -3 == x1
          const double n[3] = { 0.505*delta, 0.0, 0.0 };
          cvd.cut_surf(n,-3);
          if (cvd.ned > 0) {
            // -4 == y0
            const double n[3] = { 0.0, -0.505*delta, 0.0 };
            cvd.cut_surf(n,-4);
            if (cvd.ned > 0) {
              // -5 == y1
              const double n[3] = { 0.0, 0.505*delta, 0.0 };
              cvd.cut_surf(n,-5);
              if (cvd.ned > 0) {
                // -6 == z0
                const double n[3] = { 0.0, 0.0, -0.505*delta };
                cvd.cut_surf(n,-6);
                if (cvd.ned > 0) {
                  // -7 == z1
                  const double n[3] = { 0.0, 0.0, 0.505*delta };
                  cvd.cut_surf(n,-7);
                  if (cvd.ned > 0) {
                    // check: if we still have stuff left, all cvd.faoed[ied][1]
                    // that were -1 should have been cut off...
                    for (int ied = 0; ied < cvd.ned; ++ied) {
                      if (!(cvd.faoed[ied][1] != -1)) {
                        cvd.writeTecplot(0,xp);
                      }
                      // do not assert here; user may have purposefully down selected surfaces
                      // to reduce checking time and we want to count this simply as failing
                      // assert(cvd.faoed[ied][1] != -1);
                    }
                  }
                }
              }
            }
          }
        }
      }

      // if any surface was left,
      assert(cvd.ned > 0);

      // build the starting region considering the remaining surface patch...
      int ierr = cvd.completeCube(0.505*delta);

      if (ierr != 0) return false;
      stVec.clear();
    }
  }
  return true;
}
