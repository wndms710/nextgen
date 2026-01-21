

void testSeed(Param * param) {

  cout << "testSeed()" << endl;

  ss.ensureTeost();
  Histogram h(-1.0,1.0,MPI_COMM_NULL);

  ss.st_flag.setAll(0);

  for (int ist = 0; ist < ss.nst; ++ist) {
    const double n_st[3] = TRI_NORMAL_2(ss.xsp[ss.spost[ist][0]],ss.xsp[ss.spost[ist][1]],ss.xsp[ss.spost[ist][2]]);
    const double mag_n_st = MAG(n_st);
    if (mag_n_st > 0.0) {
      FOR_I3 {
        int ist_nbr;
        if (ss.getAlignedTriNbr(ist_nbr,ist,i)) {
          const double n_st_nbr[3] = TRI_NORMAL_2(ss.xsp[ss.spost[ist_nbr][0]],ss.xsp[ss.spost[ist_nbr][1]],ss.xsp[ss.spost[ist_nbr][2]]);
          const double mag_n_st_nbr = MAG(n_st_nbr);
          if (mag_n_st_nbr > 0.0) {
            const double dp = DOT_PRODUCT(n_st,n_st_nbr)/mag_n_st/mag_n_st_nbr;
            h.add(&dp,1);
            if (dp < -0.75) {
              ss.st_flag[ist] = 1;
              ss.st_flag[ist_nbr] = 1;
            }
          }
        }
      }
    }
  }

  h.write("hist.dat");

  // now build a table of these tris...

  int nst_problem = 0;
  for (int ist = 0; ist < ss.nst; ++ist) {
    if (ss.st_flag[ist] == 1)
    ++nst_problem;
  }
  cout << " > nst_problem: " << nst_problem << endl;

  int * st_problem = new int[nst_problem];
  int ist_problem = 0;
  for (int ist = 0; ist < ss.nst; ++ist) {
    if (ss.st_flag[ist] == 1)
    st_problem[ist_problem++] = ist;
  }
  assert(ist_problem == nst_problem);

  // set default delta values: can override...

  double delta_min = 1.0E-6*ss.getRmax();
  double delta_max = 1.0E-3*ss.getRmax();

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

  double (*bbmin)[3] = new double[ss.nst][3];
  double (*bbmax)[3] = new double[ss.nst][3];
  for (int ist = 0; ist < ss.nst; ++ist) {
    FOR_I3 bbmin[ist][i] = HUGE_VAL;
    FOR_I3 bbmax[ist][i] = -HUGE_VAL;
    FOR_I3 {
      const int isp = ss.spost[ist][i];
      // orig
      FOR_J3 bbmin[ist][j] = min(bbmin[ist][j],ss.xsp[isp][j]);
      FOR_J3 bbmax[ist][j] = max(bbmax[ist][j],ss.xsp[isp][j]);
    }
  }

  // build adt...
  Adt<double> * stAdt = new Adt<double>(ss.nst,bbmin,bbmax);
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

    const double n_st[3] = TRI_NORMAL_2(ss.xsp[ss.spost[ist_seed][0]],ss.xsp[ss.spost[ist_seed][1]],ss.xsp[ss.spost[ist_seed][2]]);
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
        FOR_I3 xp[i] = ( w0*ss.xsp[ss.spost[ist_seed][0]][i] + w1*ss.xsp[ss.spost[ist_seed][1]][i] + w2*ss.xsp[ss.spost[ist_seed][2]][i] )/(w0+w1+w2);

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

        for (int ii = 0; ii < stVec.size(); ++ii) {
          // the tri is uniquely associated with its combination of its "ist" and its bits if any...
          const int ist  = stVec[ii];
          // everyone has access to the full surface, so grab and assemble tris directly...
          // nodes first...
          FOR_I3 {
            const int isp = ss.spost[ist][i];
            map<const int,int>::iterator iter = nodeMap.find(isp);
            if (iter == nodeMap.end()) {
              const int ino = cvd.new_node();
              nodeMap[isp] = ino;
              FOR_J3 cvd.x_no[ino][j] = ss.xsp[isp][j] - xp[j];
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

          IntFlag my_st_flag(ss.nst); my_st_flag.setAll(0);
          for (int ii = 0; ii < stVec.size(); ++ii) {
            my_st_flag[stVec[ii]] = 1;
          }

          ss.writeSelectedFacesByZoneToTecplot("seed_surface_patch.dat",my_st_flag);

          cout << "press enter to continue..." << endl;
          getchar();
        }

        stVec.clear();
      }

    }

  }

}
