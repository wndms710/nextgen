
#include "StaticSolver.hpp"

void StaticSolver::flushProbes() {

  // flush all probes to disk...

  for (map<const string,PointProbe>::iterator iter = ppMap.begin(); iter != ppMap.end(); ++iter) {
    iter->second.flush();
  }

  for (map<const string,MultiFluxProbe>::iterator iter = mfpMap.begin(); iter != mfpMap.end(); ++iter) {
    iter->second.flush();
  }

  for (map<const string,FluxProbe>::iterator iter = fpMap.begin(); iter != fpMap.end(); ++iter) {
    iter->second.flush();
  }

  for (map<const string,VolumetricProbe>::iterator iter = vpMap.begin(); iter != vpMap.end(); ++iter) {
    iter->second.flush();
  }

  for (map<const string,VolumetricProbe>::iterator iter = svpMap.begin(); iter != svpMap.end(); ++iter) {
    iter->second.flush();
  }

  for (map<const string,FazoneProbe>::iterator iter = fzpMap.begin(); iter != fzpMap.end(); ++iter) {
    iter->second.flush();
  }

  for (map<const string,ConditionalProbe>::iterator iter = cpMap.begin(); iter != cpMap.end(); ++iter) {
    iter->second.flush();
  }

  for (map<const string,ConditionalProbe>::iterator iter = scpMap.begin(); iter != scpMap.end(); ++iter) {
    iter->second.flush();
  }


  for (map<const string,PdfProbe>::iterator iter = pdfpMap.begin(); iter != pdfpMap.end(); ++iter) {
    iter->second.flush();
  }

  // point cloud probes not cached because they could be big...

}

void StaticSolver::processPointProbe(Param* param,const int step,const double time,const bool killfile_request) {
  // this method only handles the initialization of probe
  // if a killfile request, the probe is sampled; if not it is simply registered
  // to the probeMap and sampled appropriately during doProbes()
  // use the full param as the hash for the map..

  // if doesn't exist in map yet, add it
  if ( ppMap.find(param->str()) == ppMap.end()) {

    // couldnt find this probe, so it needs to inserted
    pair<map<const string,PointProbe>::iterator,bool> ret =
      ppMap.insert( pair<const string,PointProbe>(param->str(),PointProbe(this)) );
    assert( ret.second);

    const int ierr = initPointProbe(ret.first->second,param);

    if ( ierr != 0) {

      // if you supply an interval for a persistent probe, then failure to
      // complete the point probe will result in termination of the solver
      // but a diagnostic probe will only yield a warning..
      if ( !killfile_request ) {
        CERR( " > unable to add PointProbe from: \n >    " << param->str() );
      } else {
        CWARN( " > unable to add PointProbe from: \n >   " << param->str() );
        ppMap.erase(ret.first);
      }
    } else {
      // if the interval is invalid, then its one-and-done for this probe
      // otherwise register for continuous processing
      if ( ret.first->second.interval == -1 ) {
        ret.first->second.doProbe(step,time);
        ret.first->second.flush();
        ppMap.erase(ret.first);
      }
    }
  }
}

int StaticSolver::initPointProbe(PointProbe& pp,Param * param) {

  COUT1("initPointProbe: \"" << param->str() << "\"...");

  bool got_name = false;
  double (*xp)[3] = NULL;
  vector<string> varNameVec;
  set<int> snapBfZoneSet;

  int iarg = 0;
  while (iarg < param->size()) {
    string token = MiscUtils::toUpperCase(param->getString(iarg++));
    if (token == "NAME") {
      pp.name = param->getString(iarg++);
      got_name = true;
    }
    else if (token == "INTERVAL" || token == "WRITE_INTERVAL") {
      pp.interval = param->getInt(iarg++); // default is -1
      if (pp.interval <= 0) {
        CWARN(" > INTERVAL expects a positive integer; treating as one-off request");
        // invalid value entered, so treat as unset (default value)
        pp.interval = -1;
      }
    }
    else if (token == "GEOM") {
      token = param->getString(iarg++);
      if (token == "CIRCLE_XRN") {
        const double x = param->getDouble(iarg++);
        const double r = param->getDouble(iarg++);
        pp.np = param->getInt(iarg++);
        assert(xp == NULL);
        xp = new double[pp.np][3];
        for (int ip = 0; ip < pp.np; ++ip) {
          const double theta = (double)(2*ip)/(double)pp.np*M_PI;
          xp[ip][0] = x;
          xp[ip][1] = r*cos(theta);
          xp[ip][2] = r*sin(theta);
        }

      }
      else if ( token == "LINE" ) {

        const double x0 = param->getDouble(iarg++);
        const double y0 = param->getDouble(iarg++);
        const double z0 = param->getDouble(iarg++);
        const double x1 = param->getDouble(iarg++);
        const double y1 = param->getDouble(iarg++);
        const double z1 = param->getDouble(iarg++);
        pp.np           = param->getInt(iarg++);

        assert( pp.np >= 2);
        assert( xp == NULL);
        xp = new double[pp.np][3];
        for (int ip = 0; ip < pp.np; ++ip) {
          xp[ip][0] = x0 + (double)ip/(double)(pp.np-1)*(x1-x0);
          xp[ip][1] = y0 + (double)ip/(double)(pp.np-1)*(y1-y0);
          xp[ip][2] = z0 + (double)ip/(double)(pp.np-1)*(z1-z0);
        }

      }
      else if ( token == "POINT") {

        pp.np     = 1;
        assert( xp == NULL);
        xp        = new double[pp.np][3];
        xp[0][0]  = param->getDouble(iarg++);
        xp[0][1]  = param->getDouble(iarg++);
        xp[0][2]  = param->getDouble(iarg++);
      }
      else if ( token == "FILE") {
        // grab xyz locations from file
        const string filename = param->getString(iarg++);
        double (*x_vec)[3];
        pp.np = MiscUtils::read3DAsciiTable(x_vec,filename);
        assert(pp.np > 0);
        assert( xp == NULL);
        xp = new double[pp.np][3];
        for (int p=0; p<pp.np; ++p) {
          FOR_I3 xp[p][i]  = x_vec[p][i];
        }
        DELETE(x_vec);
      }
      else if ( token == "SURFACE_POINT") {
        double xp_req[3];
        xp_req[0] = param->getDouble(iarg++);
        xp_req[1] = param->getDouble(iarg++);
        xp_req[2] = param->getDouble(iarg++);

        // find the nearest x_bf to the point requested..
        int ibf_zone,izone;
        double xbf_req[3];
        findClosestBf(xbf_req,ibf_zone,izone,xp_req);

        // set point to cv location...
        pp.np     = 1;
        assert( xp == NULL);
        xp        = new double[pp.np][3];

        int current_rank;
        if ( ibf_zone >= 0) {
          // you own the seed location ..
          const int icv = bfZoneVec[izone].cvobf[ibf_zone];
          FOR_I3 xp[0][i] = x_cv[icv][i];
          current_rank = mpi_rank;
        }
        else {
          current_rank = -ibf_zone-1; // owner_rank -1 indexed from findClosestBf
          assert( current_rank != mpi_rank);
          assert( (current_rank >= 0) && ( current_rank < mpi_size));
        }

        MPI_Bcast(xp[0],3,MPI_DOUBLE,current_rank,mpi_comm);

        if ( mpi_rank == 0 )
          cout << " > surface point, x_req: " << COUT_VEC(xp_req) << ", x_bf: " << COUT_VEC(xbf_req) << ", x_cv: " << COUT_VEC(xp[0]) << ", zone: " << bfZoneVec[izone].getName() << endl;
      }
      else if (token == "XDISK")  {
        const double x = param->getDouble(iarg++);
        const double r = param->getDouble(iarg++);
        const int Nradial = param->getInt(iarg++);
        const int Ntheta = param->getInt(iarg++);
        if (mpi_rank == 0) {
          cout << "Building an X-DISK PROBE IN PLANE X = " << x << " with max diameter " << r;
          cout << ", Nr = " << Nradial << " and N_t = " << Ntheta << "\n";
        }//(mpi_rank == 0)
        pp.np = Nradial*Ntheta;
        assert(xp == NULL);
        xp = new double[pp.np][3];
        int ip = 0;
        for (int i_r = 0; i_r < Nradial; ++i_r) {
          double my_r = 0 + (double)i_r/(double)(Nradial)*(r);
          for (int i_t = 0; i_t < Ntheta; ++i_t) {
            double my_t = 0 + (double)i_t/(double)(Ntheta)*2*3.14159265359;
            xp[ip][0] = x;
            xp[ip][1] = my_r*cos(my_t);
            xp[ip][2] = my_r*sin(my_t);
            ip++;
          }//(i_t < Ntheta)
        }//(i_r < Nradial)
      }//"XDISK"
      else if ( token == "SURFACE_SLICE") {
        const double px = param->getDouble(iarg++);
        const double py = param->getDouble(iarg++);
        const double pz = param->getDouble(iarg++);
        const double nx = param->getDouble(iarg++);
        const double ny = param->getDouble(iarg++);
        const double nz = param->getDouble(iarg++);
        const int nzn = bfZoneVec.size();
        bool * bfzone_flag = new bool[nzn];
        for (int izn = 0; izn < nzn; ++izn) bfzone_flag[izn] = false;

        token = param->getString(iarg);
        if (token == "FAZONE" || token == "BFZONE") {
          ++iarg;
          token = param->getString(iarg);
          while (token != "VARS") {
            ++iarg;
            map<const string,int>::iterator bf_it = bfZoneNameMap.find(token);
            if (bf_it != bfZoneNameMap.end()) {
              bfzone_flag[bf_it->second] = true;
              snapBfZoneSet.insert(bf_it->second);
            }

            // get the next zone name...
            if (iarg == param->size()) break;
            token = param->getString(iarg);
          }
        }
        else {
          // assume it is the entire domain...
          for (int izn = 0; izn < nzn; ++izn) bfzone_flag[izn] = true;
        }

        const double iso_var_value = nx*px+ny*py+nz*pz;
        vector<double> my_xps;
        FOR_IZONE(bfZoneVec) {
          if (bfzone_flag[izone]) {
            for (int ibf = bfZoneVec[izone].ibf_f; ibf <= bfZoneVec[izone].ibf_l; ++ibf) {
              const int icv = cvobf[ibf];
              int cv_flag = 0;
              for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
                const int ino = noobf_v[nob];
                const double iso_var_no = nx*x_no[ino][0]+ny*x_no[ino][1]+nz*x_no[ino][2];
                if (iso_var_no < iso_var_value) {
                  cv_flag |= 1;
                }
                else if (iso_var_no > iso_var_value) {
                  cv_flag |= 2;
                }
                else {
                  // equals zero!
                  cv_flag = 3;
                  break;
                }
                if (cv_flag == 3) {
                  break;
                }
              }

              // slice plane goes through this boundary cell...
              if (cv_flag == 3) {
                FOR_I3 my_xps.push_back(x_cv[icv][i]);
              }
            }
          }
        }
        delete[] bfzone_flag;
        assert(my_xps.size()%3 == 0);
        int my_np = my_xps.size()/3;

        // gather counts and points to every rank...

        int *recv_count = new int[mpi_size];
        MPI_Allgather(&my_np,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

        pp.np = 0;
        FOR_RANK pp.np += recv_count[rank];
        assert(xp == NULL); xp = new double[pp.np][3];

        FOR_RANK recv_count[rank] *= 3; // 3 doubles per point
        int * recv_disp = new int[mpi_size];
        recv_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

        MPI_Allgatherv(&(my_xps[0]),(int)my_xps.size(),MPI_DOUBLE,(double*)xp,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
        delete[] recv_disp;
        delete[] recv_count;

        // sort points so that probe output/order is rank independent
        vector<Triple<double,double,double> > xp_vec;
        for (int ip=0; ip<pp.np; ++ip) xp_vec.push_back(Triple<double,double,double>(xp[ip][0],xp[ip][1],xp[ip][2]));
        std::sort(xp_vec.begin(),xp_vec.end());
        for (int ip=0; ip<pp.np; ++ip) {
          xp[ip][0] = xp_vec[ip].first;
          xp[ip][1] = xp_vec[ip].second;
          xp[ip][2] = xp_vec[ip].third;
        }
      }
      else if ( token == "SURFACE_WALK") {

        double xp_seed[3];
        double direction[3];
        double max_dist = -1.0;

        token = param->getString(iarg++);

        if ( token != "X_SEED") {
          CERR( " > syntax error; should be PROBE GEOM SURFACE_WALK X_SEED <x> <y> <z> DIRECTION <nx> <ny> <nz> [MAX_DIST dist]");
        }

        xp_seed[0] = param->getDouble(iarg++);
        xp_seed[1] = param->getDouble(iarg++);
        xp_seed[2] = param->getDouble(iarg++);

        token = param->getString(iarg++);

        if ( token != "DIRECTION") {
          CERR( " > syntax error; should be PROBE GEOM SURFACE_WALK X_SEED <x> <y> <z> DIRECTION <nx> <ny> <nz> [MAX_DIST dist]");
        }

        direction[0] = param->getDouble(iarg++);
        direction[1] = param->getDouble(iarg++);
        direction[2] = param->getDouble(iarg++);

        const double tmp = MAG(direction);
        for (int i = 0; i < 3; ++i)
          direction[i] /= tmp;

        if ( param->getString(iarg) == "MAX_DIST") {

          iarg++;
          max_dist = param->getDouble(iarg++);

        }

        // now we need to perform the surface walk ... first find the x_bf closest to the x_seed point

        double x_bf_seed[3],x_cv_seed[3];
        int ibf_zone,izone;

        findClosestBf(x_bf_seed,ibf_zone,izone,xp_seed);

        if ( mpi_rank == 0 )
          cout << " > surface seed point : " << COUT_VEC(x_bf_seed) << endl;

        // we are going to walk participating cvs because there is no connectivity structure
        // that we can walk on the bfs.

        int * cv_flag = new int[ncv_g];
        for (int icv = 0; icv < ncv; ++icv)
          cv_flag[icv] = -2; // cant participate in the walk..

        /*
        for (int ibf = 0; ibf < bfZoneVec[izone].nbf; ++ibf) {

          const int icv = bfZoneVec[izone].cvobf[ibf];
          cv_flag[icv]  = -1; // potentially can participate, but are not touched now

        }
        */


        // allow all boundary surfaces to participate ...
        for (int ibf = 0; ibf < nbf; ++ibf) {
          cv_flag[cvobf[ibf]] = -1;
        }

        updateCvData(cv_flag);

        int icv_seed;
        int current_rank;

        if ( ibf_zone >= 0) {

          // you own the seed location ..
          icv_seed = bfZoneVec[izone].cvobf[ibf_zone];
          current_rank = mpi_rank;
          for (int i = 0; i < 3; ++i)
            x_cv_seed[i] = x_cv[icv_seed][i];

          // set the flag for the initial seed...

          cv_flag[icv_seed] = 1;

        } else {

          icv_seed = -1;
          current_rank = -ibf_zone-1;
          assert( current_rank != mpi_rank);
          assert( (current_rank >= 0) && ( current_rank < mpi_size));

        }

        MPI_Bcast(x_cv_seed,3,MPI_DOUBLE,current_rank,mpi_comm);

        int done = 0;
        double tot_dist = 0.0;
        while ( done == 0) {

          // march via cvocv.. assume that we are done ..

          done = 1;

          uint8 rbi_nbr = 0;

          while ( icv_seed >= 0) {

            double dp_max = 0.0;
            int ii        = -1;
            for (int coc = cvocv_i[icv_seed]; coc != cvocv_i[icv_seed+1]; ++coc) {

              const int icv_nbr    = cvocv_v[coc];
              if ( cv_flag[icv_nbr] == -1) {

                double dx[3]         = DIFF(x_cv[icv_nbr],x_cv[icv_seed]);
                const double tmp     = MAG(dx);
                const double this_dp = DOT_PRODUCT(dx,direction)/tmp; // direction is unit vector arleady.

                if ( this_dp > dp_max ) {
                  dp_max    = this_dp;
                  ii        = icv_nbr;
                }
              }
            }

            if ( ii >= 0 && ( (tot_dist <= max_dist) || (max_dist < 0.0)) ) {

              tot_dist += DIST(x_cv[ii],x_cv[icv_seed]);

              if ( ii < ncv) {

                icv_seed = ii;
                cv_flag[icv_seed] = 1;

              } else {

                assert( ii < ncv_g);
                icv_seed = -1;
                rbi_nbr  = rbi_g[ii-ncv];
                done     = 0;

              }

            } else {
              icv_seed = -1;
              done     =  1;
            }
          }

          MPI_Bcast(&done,1,MPI_INT,current_rank,mpi_comm);

          if ( done == 0) {

            updateCvData(cv_flag);

            // we could make done an int8 and bcast all together

            MPI_Bcast(&rbi_nbr,1,MPI_INT8,current_rank,mpi_comm);
            MPI_Bcast(&tot_dist,1,MPI_DOUBLE,current_rank,mpi_comm);

            int bits,index;
            BitUtils::unpackRankBitsIndex(current_rank,bits,index,rbi_nbr);

            if ( mpi_rank == current_rank ) {
              icv_seed = index;
              cv_flag[icv_seed] = 1;
            } else {
              icv_seed = -1;
            }

          } else {

            assert( done == 1);

          }

        }

        int np_local = 0;
        FOR_ICV if (cv_flag[icv] == 1) np_local++;
        double (*xp_local)[3] = new double[np_local][3];
        np_local = 0;
        FOR_ICV if (cv_flag[icv] == 1) {
          FOR_I3 xp_local[np_local][i] = x_cv[icv][i];
          np_local++;
        }
        delete[] cv_flag;

        int *recv_count = new int[mpi_size];
        MPI_Allgather(&np_local,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

        pp.np = 0;
        FOR_RANK pp.np += recv_count[rank];
        xp = new double[pp.np][3];

        FOR_RANK recv_count[rank] *= 3; // 3 doubles per point
        int * recv_disp = new int[mpi_size];
        recv_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

        MPI_Allgatherv((double*)xp_local,np_local*3,MPI_DOUBLE,(double*)xp,recv_count,
                       recv_disp,MPI_DOUBLE,mpi_comm);


        vector<pair<double,int> > dist_vec;
        for (int ip = 0; ip < pp.np; ++ip) {

          const double this_dist = DIST(xp[ip],x_cv_seed);
          dist_vec.push_back(pair<double,int>(this_dist,ip));

        }

        sort(dist_vec.begin(), dist_vec.end());

        double (*xp_reorder)[3] = new double[pp.np][3];
        for (int ip = 0; ip < pp.np; ++ip) {

          const int ip_old = dist_vec[ip].second;
          assert( (ip_old >= 0) && (ip_old < pp.np));

          for (int i = 0; i < 3; ++i)
            xp_reorder[ip][i] = xp[ip_old][i];

        }

        delete[] xp; xp = xp_reorder;

        delete[] xp_local;
        delete[] recv_count;
        delete[] recv_disp;
      }
      else {
        CWARN("unrecognized PROBE GEOM: " << token << ". Skipping.");
      }
    }
    else if ((token == "VAR")||(token == "VARS")) {
      // for vars, don't push directly into pp's varNameVec, because
      // we want to make sure they eval to CV_DN or CV_DN3...
      while (iarg < param->size())
      varNameVec.push_back(param->getString(iarg++));
    }
    else if (token == "SNAP_TO_SURFACE") {
      pp.b_surface = true;
    }
    else if (token == "SNAP_TO_CHT") {
      pp.b_cht = true;
    }
    else if (token == "SNAP_FAZONES") {
      const string bfzonesCsv = param->getString(iarg++);
      vector<string> bfzonesVec;
      MiscUtils::splitCsv(bfzonesVec,bfzonesCsv);
      for (vector<string>::const_iterator it=bfzonesVec.begin(); it!=bfzonesVec.end(); ++it) {
        map<const string,int>::iterator bf_it = bfZoneNameMap.find(*it);
        if (bf_it != bfZoneNameMap.end()) {
          snapBfZoneSet.insert(bf_it->second);
        }
        else {
          CWARN("did not find SNAP_FAZONES zone: " << *it << "; skipping");
        }
      }
    }
    else {
      COUT1("processPointProbe: skipping unrecognized PROBE token \"" << token << "\"");
    }
  }

  int ierr = 0;
  if (!got_name) {
    if (mpi_rank == 0)
    cout << " > Warning: PROBE missing NAME." << endl;
    ierr = -1;
  }
  if (xp == NULL) {
    if (mpi_rank == 0)
    cout << " > Warning: PROBE missing GEOM." << endl;
    ierr = -1;
  }
  if (pp.b_surface && pp.b_cht) {
    if (mpi_rank == 0)
    cout << " > Warning: Cannot use both SNAP_TO_SURFACE and SNAP_TO_CHT simultaneously." << endl;
    ierr = -1;
  }
  if (ierr != 0) {
    if (mpi_rank == 0)
    cout << " > Warning: Skipping this PROBE. Correct syntax warnings." << endl;
    DELETE(xp);
    return -1;
  }

  // if snapped-to-surface (accesses bf data), need to store which ibf
  int * ibf_zone = NULL;
  int * izone = NULL;
  IntFlag active_zones(bfZoneVec.size());  // for surface snapped, keeps track of participating zones
  active_zones.setAll(0);

  double (*x_snap)[3] = NULL;
  if (pp.b_surface) {
    COUT1(" > probe is being snapped to surface; it can only access BF_DATA");
    // update probe locations based on nearest bf
    // find the nearest x_bf to the point requested..
    assert(ibf_zone == NULL);
    assert(izone == NULL);
    ibf_zone = new int[pp.np];
    izone = new int[pp.np];
    x_snap = new double[pp.np][3];

    if (snapBfZoneSet.empty()) {
      findClosestBf(x_snap,ibf_zone,izone,xp,pp.np);  // each point probe in array should now know its ibf/izone in case bf_data requested
    }
    else {
      findClosestBfOnZones(x_snap,ibf_zone,izone,xp,pp.np,snapBfZoneSet);
    }

    // set point to cv location...
    for (int ip=0; ip<pp.np; ++ip) {
      if ( ibf_zone[ip] >= 0) {
        // you own the seed location ..

        const int icv = bfZoneVec[izone[ip]].cvobf[ibf_zone[ip]];
        FOR_I3 xp[ip][i] = x_cv[icv][i];  // need to store x_cv for now b/c used for load balancing

        active_zones[izone[ip]] |= 2;  // izone was distributed in findClosestBf, this count is rank_local
      }
      else {
        FOR_I3 xp[ip][i] = 0.0;  // reduction via summation
        const int current_rank = -ibf_zone[ip]-1; // owner_rank -1 indexed from findClosestBf
        assert( current_rank != mpi_rank);
        assert( (current_rank >= 0) && ( current_rank < mpi_size));
      }
      active_zones[izone[ip]] |= 1;  // izone was distributed in findClosestBf
    }
    MPI_Allreduce(MPI_IN_PLACE,(double*)xp,pp.np*3,MPI_DOUBLE,MPI_SUM,mpi_comm);

    if (mpi_rank == 0) {
      cout << "    > number of paticipating boundary zones: " << active_zones.countPositive() << endl;
    }
  }
  else if (pp.b_cht) {
    COUT1(" > probe is being snapped to CHT; it can only access cht nodal data");
    // update probe locations based on nearest cht node
    // find the nearest cht:x_no to the point requested
    // will reuse x_snap above
    assert(izone == NULL);  // eventually may need cht-zone specific snapping...
    //izone = new int[pp.np];
    x_snap = new double[pp.np][3];
    for (int ip=0; ip < pp.np; ++ip) {
      // izone[ip] = -1;
      FOR_J3 x_snap[ip][j] = 0.0;
    }
    // during load balancing we need to determine which ranks own the nearest point
    // so defer the determination of the nearest point to then as well
  }

  // figure out the variables using CtiRegister::getCtiData...
  assert(pp.varVec.empty());
  vector<int> ctiDataType;
  for (int ii = 0,ii_end = varNameVec.size(); ii < ii_end; ++ii) {
    // registered data is persistent/immutable (i.e. mapped ptrs do not change), so we store
    // the data ptr when the data is registered to save time on future evaluations...

    vector<CtiRegister::CtiData*> zoneVarVec;  // holds registered data pointers per-active-zone

    if (varNameVec[ii].find("*:") != string::npos) {
      // requires zonename expansion, assumption is that it must access BF_DATA...
      // test first to ensure this is surface-snapped
      if ((ibf_zone == NULL)||(izone == NULL)||(!pp.b_surface)) {
        if (mpi_rank == 0) cout << " > Warning: PointProbe data \"" << varNameVec[ii] << "\" requires zonename expansion, but probe is not snapped-to-surface. Skipping." << endl;
        continue;
      }

      // ensure variable is available on all participating zones
      //TODO possibly we may want to store a vector<vector> of the registered data pointers; currently we re-fetch during each doProbe
      int d_type = 0;  // data_type (dn vs. dn3)
      bool b_problem = false;
      for (int izn=0, nzn=bfZoneVec.size(); izn<nzn; ++izn) {
        if (active_zones[izn]) {
          string data_name = varNameVec[ii];
          data_name = MiscUtils::replaceAll(data_name,"*:",(bfZoneVec[izn].getName()+":"));

          bool b_registered = true;
          CtiRegister::CtiData * data = CtiRegister::getRegisteredCtiData(data_name,false);
          if (data == NULL) {
            data = CtiRegister::getUnregisteredCtiData(data_name);
            b_registered = false;
          }

          if ((data == NULL)||(data->getUnindexedTopology() != BF_DATA)) {
            if (mpi_rank == 0) cout << " > Warning: PointProbe is surface-snapped but data \"" << varNameVec[ii] << "\" does not evaluate to BF_DATA. Skipping." << endl;
            b_problem = true;
            break;
          }

          if (data->getType() == DN_DATA) {
            d_type |= 1;
          }
          else if (data->getType() == DN3_DATA) {
            d_type |= 2;
          }
          else {
            if (mpi_rank == 0) cout << " > Warning: PointProbe data \"" << varNameVec[ii] << "\" does not evaluate to DN or DN3 on BF_ZONE " << bfZoneVec[izn].getName() << ". Skipping." << endl;
            b_problem = true;
            break;
          }

          // add to list of rank_local active-zone-data-pointers
          if (active_zones[izn] == 3) {
            if (b_registered) zoneVarVec.push_back(data);
            else zoneVarVec.push_back(NULL);
          }
        }
      }

      if ((d_type != 1)&&(d_type !=2)) {
        if (mpi_rank == 0) cout << " > Warning: PointProbe data \"" << varNameVec[ii] << "\" doesn't have consistent data type across zones. Skipping." << endl;
        b_problem = true;
      }

      if (b_problem) continue;  // skip registration of this variable

      // seems legit; register with probe and output info
      if (d_type == 1) {
        COUT1(" > adding scalar \"" << varNameVec[ii] << "\"");
        pp.ns += 1;
      }
      else if (d_type == 2) {
        COUT1(" > adding vector \"" << varNameVec[ii] << "\"");
        pp.ns += 3;
      }
      else {
        assert(0);
      }
      // register with NULL - zone expansion will be caught by name
      pp.varVec.push_back(pair<string,vector<CtiRegister::CtiData*> >(varNameVec[ii],zoneVarVec));
      ctiDataType.push_back(d_type); // just used to determine whether dn or dn3 for header writing
    }
    else {
      // this should be an unexpanded variable, i.e., a named var or a named-zone:var pair
      bool b_registered = true;
      CtiRegister::CtiData * data = CtiRegister::getRegisteredCtiData(varNameVec[ii],false); // verbose == false (no warning if NULL)
      if (data == NULL) {
        b_registered = false;
        data = CtiRegister::getUnregisteredCtiData(varNameVec[ii]);
      }

      if ((data == NULL)||((data->getType() != DN_DATA)&&(data->getType() != DN3_DATA))) {
        COUT1(" > Warning: PointProbe data \"" << varNameVec[ii] << "\" does not evaluate to DN or DN3. Skipping.");
        continue;
      }

      if (data->getUnindexedTopology() == BF_DATA) {
        // if it is BF_DATA make sure this is a surface-snapped probe
        if ((ibf_zone == NULL)||(izone == NULL)||(!pp.b_surface)) {
          COUT1(" > Warning: PointProbe data \"" << varNameVec[ii] << "\" evaluates to BF_DATA, but SNAP_TO_SURFACE flag was not specified. Skipping.");
          continue;
        }
        // also need to ensure this data exists on all probe ibfs
        //TODO could reduce this requirement and just return "0" during doProbe if not available...
        bool b_exit = false;
        for (int ip=0; ip<pp.np; ++ip) {
          if (izone[ip] != data->getIndex()) {
            COUT1(" > Warning: PointProbe data \"" << varNameVec[ii] << "\" doesn't exist at all probe locations. If probe points span multiple zones, consider using the wildcard notation \"*:<var>\". Skipping.");
            b_exit = true;
            break;
          }
        }
        if (b_exit) continue;
      }
      else if ((data->getName() == "cht")&&(!pp.b_cht)) {
        COUT1(" > Warning: PointProbe variable \"" << varNameVec[ii] << "\" is CHT data but SNAP_TO_CHT flag was not specified. Skipping.");
        continue;
      }
      else if ((data->getUnindexedTopology() == CV_DATA)&&pp.b_surface) {
        COUT1(" > Warning: PointProbe is surface-snapped but \"" << varNameVec[ii] << "\" evaluates to CV_DATA. If this is intended please use the proj(CV_DATA) operator. Skipping.");
        continue;
      }
      else if ((data->getUnindexedTopology() == CV_DATA)&&pp.b_cht) {
        COUT1(" > Warning: PointProbe is specified as CHT data but \"" << varNameVec[ii] << "\" evaluates to CV_DATA. Skipping.");
        continue;
      }
      else if ((data->getUnindexedTopology() != CV_DATA)&&(!pp.b_cht)) {
        // otherwise make sure it is CV_*...
        COUT1(" > Warning: PointProbe data \"" << varNameVec[ii] << "\" does not evaluate to BF_DATA or CV_DATA. Skipping.");
        continue;
      }

      if (b_registered) {
        zoneVarVec.push_back(data);
        pp.varVec.push_back(pair<string,vector<CtiRegister::CtiData*> >(varNameVec[ii],zoneVarVec));
      }
      else {
        zoneVarVec.push_back(NULL);
        pp.varVec.push_back(pair<string,vector<CtiRegister::CtiData*> >(varNameVec[ii],zoneVarVec));
      }
      ctiDataType.push_back((data->getType() == DN_DATA) ? 1:2);

      if (data->getType() == DN_DATA) {
        if (mpi_rank == 0) cout << " > adding scalar \"" << varNameVec[ii] << "\"" << endl;
        pp.ns += 1;
      }
      else if (data->getType() == DN3_DATA) {
        if (mpi_rank == 0) cout << " > adding vector \"" << varNameVec[ii] << "\"" << endl;
        pp.ns += 3;
      }
      else {
        assert(0);
      }
    }
  }

  if (pp.varVec.empty()) {
    CWARN("PointProbe has no valid data. Skipping.");
    delete[] xp;
    return -1;
  }

  // -----------------------------------
  // now decide who owns each probe...
  // -----------------------------------

  if (mpi_rank == 0)
    cout << " > requesting: " << pp.np << " point probes...";

  DoubleInt * my_di = new DoubleInt[pp.np];
  int * icv_closest = new int[pp.np];  // icv or cht:ino, depending

  if (pp.b_cht) {

    // for faster point lookup, build the xnoAdt...but what to use as bbox delta?
    // no access to FlowSolver cht, but maybe could build the noote from registered data
    // and do bboxes on based on tets...

    CtiRegister::CtiData * x_no_data = CtiRegister::getRegisteredCtiData("cht:x_no");
    assert(x_no_data);
    assert(x_no_data->getType() == DN3_DATA);

    const int nno = x_no_data->size();
    double (*x_data)[3] = x_no_data->getDN3ptr();
    // here we use the MPI_MINLOC capability to break any containment ties for
    // probe points...

    // brute force b/c not using adt
    for (int ip = 0; ip < pp.np; ++ip) {
      icv_closest[ip]       = -1;
      my_di[ip].this_double = HUGE_VAL;

      for (int ino=0; ino < nno; ++ino) {
        const double d2 = DIST2(x_data[ino],xp[ip]);
        if ((icv_closest[ip] == -1)||(d2 < my_di[ip].this_double)) {
          icv_closest[ip]      = ino;
          my_di[ip].this_double = d2;
          FOR_J3 x_snap[ip][j] = x_data[ino][j];  // store local min location for each probe
        }
      }

      my_di[ip].this_int = -1; // set the int to -1 unless we found one...
      if (icv_closest[ip] != -1) my_di[ip].this_int = mpi_rank;
    }

    // eventually use bbox based approach similar to cv below

  }
  else {
    // cv or bf-based data
    // even for BF_DATA probes still use this routine to see which rank own which cv, since ibf is tied to icv anyway

    // for fast point lookup, build the cvAdt...
    if (cvAdt == NULL) buildCvAdt();

    // here we use the MPI_MINLOC capability to break any containment ties for
    // probe points...

    vector<int> bboxVec;
    for (int ip = 0; ip < pp.np; ++ip) {
      assert(bboxVec.empty());
      cvAdt->buildListForPoint(bboxVec,xp[ip]);
      icv_closest[ip]       = -1;
      my_di[ip].this_double = 1.0E+20; // big number
      for (int ii = 0,ii_end = bboxVec.size(); ii < ii_end; ++ii) {
        const int icv = bboxVec[ii];
        const double d2 = DIST2(x_vv[icv],xp[ip]);
        if ((icv_closest[ip] == -1)||(d2 < my_di[ip].this_double)) {
          icv_closest[ip]      = icv;
          my_di[ip].this_double = d2;
        }
      }
      bboxVec.clear();
      my_di[ip].this_int = -1; // set the int to -1 unless we found one...
      if (icv_closest[ip] != -1) my_di[ip].this_int = mpi_rank;
    }
  }

  // TODO: Thunder,Mira,Honda SGI (old) failed on large? np's using the following...
  // Thunder fixed with compiler params. Mira was failing on 6M points. Honda?
  DoubleInt * di = new DoubleInt[pp.np];
  MPI_Allreduce(my_di,di,pp.np,MPI_DOUBLE_INT,MPI_MINLOC,mpi_comm);
  delete[] my_di;

  if (pp.b_surface) {
    // for bf snapped probes swap the xp and x_snap locations
    // so probes know their x_bf locations and not the cvobf(ibf) location (which was used for cv-based load balancing)
    double max_d2 = -HUGE_VAL;
    for (int ip=0; ip<pp.np; ++ip) {
      if (di[ip].this_double > max_d2) max_d2 = di[ip].this_double;
      if ( ibf_zone[ip] >= 0) {
        // you own the seed location ..
        // Voronoi diagram should ensure bf is closest to cv anyway, butin case of orphans or something
        // we use x_cv in the alg. above
        FOR_I3 xp[ip][i] = x_snap[ip][i];  // snapped probes store their surface locations now
      }
    }
    COUT1(" > maximum snapped distance: " << sqrt(max_d2));
  }
  else if (pp.b_cht) {
    // make sure all probes know actual probe locations (nearest cht nodes)
    double max_d2 = -HUGE_VAL;
    for (int ip=0; ip<pp.np; ++ip) {
      if (di[ip].this_double > max_d2) max_d2 = di[ip].this_double;

      if (di[ip].this_int == mpi_rank) {
        FOR_J3 xp[ip][j] = x_snap[ip][j];
      }
      else {
        FOR_J3 xp[ip][j] = 0.0;  // will use MPI_SUM
      }
    }
    MPI_Allreduce(MPI_IN_PLACE,(double*)xp,pp.np*3,MPI_DOUBLE,MPI_SUM,mpi_comm);
    COUT1(" > maximum snapped distance: " << sqrt(max_d2));
  }

  DELETE(x_snap);

  // put the probe file i/o on the rank with the largest number of matched points: i.e.
  // minimize message size associated with probing...

  // set rank0...

  {

    pp.rank0 = 0;
    pp.active_rank_count = 0;
    map<const int,int> rankCountMap;
    for (int ip = 0; ip < pp.np; ++ip) {
      map<const int,int>::iterator iter = rankCountMap.find(di[ip].this_int);
      if (iter == rankCountMap.end()) {
        rankCountMap[di[ip].this_int] = 1;
        // when we first encounter a non-(-1) rank, increment the count...
        if (di[ip].this_int != -1) {
          assert((di[ip].this_int >= 0)&&(di[ip].this_int < mpi_size));
          ++pp.active_rank_count;
        }
      }
      else {
        iter->second += 1;
      }
    }

    int probe_count_max = 0;
    for (map<const int,int>::iterator iter = rankCountMap.begin(); iter != rankCountMap.end(); ++iter) {
      if (iter->first == -1) {
        if (mpi_rank == 0)
        cout << " > WARNING: " << iter->second << " of " << pp.np << " probe points could not be located. These will appear as zeros in the probe file." << endl;
      }
      else if (iter->second > probe_count_max) {
        pp.rank0 = iter->first;
        probe_count_max = iter->second;
      }
    }

    if (mpi_rank == 0)
    cout << " > rank0 set to " << pp.rank0 << " because it owns " << probe_count_max << " of " << pp.np << " probe points" << endl;

    // on the active rank, build the rank-of-active-rank...

    assert(pp.raoar == NULL);
    if (mpi_rank == pp.rank0) {

      pp.raoar = new int[pp.active_rank_count];

      int irank = 0;
      for (map<const int,int>::iterator iter = rankCountMap.begin(); iter != rankCountMap.end(); ++iter) {
        if (iter->first != -1) {
          assert((iter->first >= 0)&&(iter->first < mpi_size));
          pp.raoar[irank] = iter->first;
          ++irank;
        }
      }
      assert(irank == pp.active_rank_count);

    }

  }

  // now figure out pack and unpack buffers...

  // pack first...

  pp.pack_buf_size = 0;
  for (int ip = 0; ip < pp.np; ++ip) {
    if (di[ip].this_int == mpi_rank) {
      pp.pack_buf_size += 1; // eventually this can be larger to handle weighted data...
    }
  }

  int * pack_buf_int = NULL; // to store the ip...

  assert(pp.pack_buf_index == NULL);
  if (pp.b_surface) {
    assert(pp.pack_buf_zone == NULL);
    assert((ibf_zone != NULL)&&(izone != NULL));
  }

  if (pp.pack_buf_size > 0) {
    pp.pack_buf_index = new int[pp.pack_buf_size];  // icv or ibf
    if (pp.b_surface) pp.pack_buf_zone = new int[pp.pack_buf_size];
    pp.pack_buf_wgt   = new double[pp.pack_buf_size];
    pack_buf_int = new int[pp.pack_buf_size];
    int ipack = 0;
    for (int ip = 0; ip < pp.np; ++ip) {
      if (di[ip].this_int == mpi_rank) {
        assert(icv_closest[ip] >= 0);
        if (pp.b_surface) {
          pp.pack_buf_index[ipack] = ibf_zone[ip];
          pp.pack_buf_zone[ipack] = izone[ip];
        }
        else {
          pp.pack_buf_index[ipack] = icv_closest[ip];
        }
        pp.pack_buf_wgt[ipack]   = 1.0; // unit wgt for now
        pack_buf_int[ipack]      = ip; // the associated probe
        ++ipack;
      }
    }
    assert(ipack == pp.pack_buf_size);

    // specifically for surface-snapped probes, we want to resort the buf_index by zone
    // to leverage the zone-based zoneVarVec
    if (pp.b_surface) {
      vector<pair<int,int> > v_zone_ipack;  // zone, ipack
      int * _pack_buf_index  = new int[pp.pack_buf_size];
      int * _pack_buf_zone   = new int[pp.pack_buf_size];
      double * _pack_buf_wgt = new double[pp.pack_buf_size];
      int * _pack_buf_int    = new int[pp.pack_buf_size];

      for (ipack=0; ipack < pp.pack_buf_size; ++ipack) {
        v_zone_ipack.push_back(pair<int,int> (pp.pack_buf_zone[ipack],ipack));
        _pack_buf_index[ipack] = pp.pack_buf_index[ipack];
        _pack_buf_zone[ipack] = pp.pack_buf_zone[ipack];
        _pack_buf_wgt[ipack] = pp.pack_buf_wgt[ipack];
        _pack_buf_int[ipack] = pack_buf_int[ipack];
      }
      // sort based on zone
      std::sort(v_zone_ipack.begin(),v_zone_ipack.end());

      for (ipack=0; ipack < pp.pack_buf_size; ++ipack) {
        const int sorted_index = v_zone_ipack[ipack].second;
        pp.pack_buf_index[ipack] = _pack_buf_index[sorted_index];
        pp.pack_buf_zone[ipack] = _pack_buf_zone[sorted_index];
        pp.pack_buf_wgt[ipack] = _pack_buf_wgt[sorted_index];
        pack_buf_int[ipack] = _pack_buf_int[sorted_index];
      }

      delete[] _pack_buf_index;
      delete[] _pack_buf_zone;
      delete[] _pack_buf_wgt;
      delete[] _pack_buf_int;
    }
  }

  delete[] di;
  delete[] icv_closest;
  DELETE(ibf_zone);
  DELETE(izone);

  // and exchange the counts...
  // even though the rank0 process could determine the counts for this case, we
  // prefer this more general approach so that the ranks can send back a more
  // complex pattern. To be honest, the number of participating ranks could grow
  // too, so changing from closest point to some sort of reconstruction will
  // be quite a bit of work...

  assert(pp.ioar_i == NULL);
  assert(pp.ioar_v == NULL);

  if (mpi_rank == pp.rank0) {

    pp.ioar_i = new int[pp.active_rank_count+1];

    // pull counts from active ranks...

    for (int iar = 0; iar < pp.active_rank_count; ++iar) {
      const int rank = pp.raoar[iar]; assert((rank >= 0)&&(rank < mpi_size));
      if (rank == mpi_rank) {
        pp.ioar_i[iar+1] = pp.pack_buf_size;
      }
      else {
        MPI_Recv(pp.ioar_i+iar+1,1,MPI_INT,rank,52526,mpi_comm,MPI_STATUS_IGNORE);
      }
    }

    pp.ioar_i[0] = 0;
    for (int iar = 0; iar < pp.active_rank_count; ++iar) {
      //cout << " RECIEVING " << pp.ioar_i[iar+1] << " from rank: " << pp.raoar[iar] << endl;
      pp.ioar_i[iar+1] += pp.ioar_i[iar];
    }

    pp.unpack_buf_size = pp.ioar_i[pp.active_rank_count];
    pp.ioar_v = new int[pp.unpack_buf_size];
    for (int ioa = 0; ioa < pp.unpack_buf_size; ++ioa) pp.ioar_v[ioa] = -1;

    // pull index data from active ranks...

    for (int iar = 0; iar < pp.active_rank_count; ++iar) {
      const int rank = pp.raoar[iar]; assert((rank >= 0)&&(rank < mpi_size));
      if (rank == mpi_rank) {
        for (int ioa = pp.ioar_i[iar]; ioa != pp.ioar_i[iar+1]; ++ioa) {
          pp.ioar_v[ioa] = pack_buf_int[ioa-pp.ioar_i[iar]];
        }
      }
      else {
        MPI_Recv(pp.ioar_v+pp.ioar_i[iar],pp.ioar_i[iar+1]-pp.ioar_i[iar],MPI_INT,rank,52527,mpi_comm,MPI_STATUS_IGNORE);
      }
    }

    // check index corresponds to one of the probe points...

    for (int ioa = 0; ioa < pp.unpack_buf_size; ++ioa)
    assert((pp.ioar_v[ioa] >= 0)&&(pp.ioar_v[ioa] < pp.np));

  }
  else if (pp.pack_buf_size > 0) {

    MPI_Ssend(&pp.pack_buf_size,1,MPI_INT,pp.rank0,52526,mpi_comm);

    MPI_Ssend(pack_buf_int,pp.pack_buf_size,MPI_INT,pp.rank0,52527,mpi_comm);

  }

  DELETE(pack_buf_int);

  // ---------------------------
  // final prep...
  // ---------------------------

  assert(pp.ssVec.empty());
  if (mpi_rank == pp.rank0) {

    // build and initiate ssVec on rank0...
    // there is one entry in ssVec for each probe file.

    // assert(ctiDataVec.size() == pp.varVec.size());
    assert(ctiDataType.size() == pp.varVec.size());
    for (int ii = 0,ii_end = ctiDataType.size(); ii < ii_end; ++ii) {
      if (ctiDataType[ii] == 1) {
        // DN_DATA
        pp.ssVec.push_back(pair<string,std::stringstream*>(pp.name+"."+pp.varVec[ii].first,new std::stringstream));
        *pp.ssVec.back().second << setprecision(8);
        *pp.ssVec.back().second << "# " << param->str() << endl;
        *pp.ssVec.back().second << "# mles hash id " << RestartHashUtilities::mlesHash << endl;
        *pp.ssVec.back().second << "# sles hash id " << RestartHashUtilities::slesHash << endl;
        *pp.ssVec.back().second << "# VAR " << pp.varVec[ii].first << endl;
        *pp.ssVec.back().second << "# 1:step 2:time 3:count 4-" << pp.np+3 << ":data" << endl;
      }
      else if (ctiDataType[ii] == 2) {
        // DN3_DATA
        // x...
        pp.ssVec.push_back(pair<string,std::stringstream*>(pp.name+"."+pp.varVec[ii].first+"-x",new std::stringstream));
        *pp.ssVec.back().second << setprecision(8);
        *pp.ssVec.back().second << "# " << param->str() << endl;
        *pp.ssVec.back().second << "# mles hash id " << RestartHashUtilities::mlesHash << endl;
        *pp.ssVec.back().second << "# sles hash id " << RestartHashUtilities::slesHash << endl;
        *pp.ssVec.back().second << "# VAR " << pp.varVec[ii].first << "-x" << endl;
        *pp.ssVec.back().second << "# 1:step 2:time 3:count 4-" << pp.np+3 << ":data" << endl;
        // y...
        pp.ssVec.push_back(pair<string,std::stringstream*>(pp.name+"."+pp.varVec[ii].first+"-y",new std::stringstream));
        *pp.ssVec.back().second << setprecision(8);
        *pp.ssVec.back().second << "# " << param->str() << endl;
        *pp.ssVec.back().second << "# mles hash id " << RestartHashUtilities::mlesHash << endl;
        *pp.ssVec.back().second << "# sles hash id " << RestartHashUtilities::slesHash << endl;
        *pp.ssVec.back().second << "# VAR " << pp.varVec[ii].first << "-y" << endl;
        *pp.ssVec.back().second << "# 1:step 2:time 3:count 4-" << pp.np+3 << ":data" << endl;
        // z...
        pp.ssVec.push_back(pair<string,std::stringstream*>(pp.name+"."+pp.varVec[ii].first+"-z",new std::stringstream));
        *pp.ssVec.back().second << setprecision(8);
        *pp.ssVec.back().second << "# " << param->str() << endl;
        *pp.ssVec.back().second << "# mles hash id " << RestartHashUtilities::mlesHash << endl;
        *pp.ssVec.back().second << "# sles hash id " << RestartHashUtilities::slesHash << endl;
        *pp.ssVec.back().second << "# VAR " << pp.varVec[ii].first << "-z" << endl;
        *pp.ssVec.back().second << "# 1:step 2:time 3:count 4-" << pp.np+3 << ":data" << endl;
      }
      else {
        assert(0);
      }
    }

    // make sure the file can exist (i.e. build sub-directories if requested in name)...

    MiscUtils::mkdir_for_file(pp.name);

    // write the README file containing the coordinates...

    const string filename = pp.name+".README";
    const string tmp_filename = MiscUtils::makeTmpPrefix(filename);

    ofstream outFile;
    if (pp.interval == -1)
      outFile.open(tmp_filename.c_str(),ofstream::trunc);
    else
      outFile.open(filename.c_str(),ofstream::app);
    assert(outFile.is_open());
    outFile << setprecision(8);
    outFile << "# " << param->str() << endl;
    outFile << "# mles hash id " << RestartHashUtilities::mlesHash << endl;
    outFile << "# sles hash id " << RestartHashUtilities::slesHash << endl;
    outFile << "# 1:probe 2:x 3:y 4:z" << endl;
    for (int ip = 0; ip < pp.np; ++ip) outFile << ip << " " << xp[ip][0] << " " << xp[ip][1] << " " << xp[ip][2] << endl;
    outFile.close();
    if (pp.interval == -1) {
      remove(filename.c_str());
      rename(tmp_filename.c_str(),filename.c_str());
    }

    // for the case of single-file-mode, put the header in pp.ss right now. Otherwise,
    // we store the header to dump at the top of each file...

    /*
      if (mfp.format == SINGLE_FILE) {
      *mfp.ss << hss.rdbuf();
      }
      else {
      assert(mfp.format == MULTI_FILE);
      mfp.header = hss.str();
      }
    */

  }

  delete[] xp;

  return 0; // OK

}

void StaticSolver::findClosestBf(double x_bf_seed[3], int& ibf_zone, int &izone, const double xp_seed[3]) const {

  ibf_zone = -1;
  izone    = -1;

  DoubleInt d2_rank;
  d2_rank.this_int = mpi_rank;
  d2_rank.this_double = 1.0e+20;

  for (int ibf = 0; ibf < nbf; ++ibf) {

    double this_d2 = DIST2(xp_seed,x_bf[ibf]);

    if ( this_d2 < d2_rank.this_double) {
      d2_rank.this_double = this_d2;
      ibf_zone            = ibf; // rank local ibf is stored right now...
    }
  }

  DoubleInt d2_rank_global;
  MPI_Allreduce(&d2_rank,&d2_rank_global,1,MPI_DOUBLE_INT,MPI_MINLOC,mpi_comm);

  if ( mpi_rank == d2_rank_global.this_int) {

    for (int i =0; i < 3; ++i)
      x_bf_seed[i] = x_bf[ibf_zone][i];

    izone = zone_bf[ibf_zone];

    // return the ibf_zone to a zone local indexing

    ibf_zone -= bfZoneVec[izone].ibf_f;
    assert( (ibf_zone >= 0) && (ibf_zone < bfZoneVec[izone].nbf));

  } else {

    ibf_zone = -d2_rank_global.this_int-1; // you dont own the seed, get -1 indexing of rank that does
    izone    = -1;

  }

  MPI_Bcast(x_bf_seed,3,MPI_DOUBLE,d2_rank_global.this_int,mpi_comm);
  MPI_Bcast(&izone, 1, MPI_INT, d2_rank_global.this_int, mpi_comm);

  assert( (izone >= 0) && (izone < int(bfZoneVec.size())));
}

void StaticSolver::findClosestBf(double (*x_bf_seed)[3], int * ibf_zone, int * izone, const double (*xp_seed)[3],const int np) const {

  DoubleInt * d2_rank = new DoubleInt[np];

  for (int ip=0; ip<np; ++ip) {
    ibf_zone[ip] = -1;
    izone[ip] = -1;
    d2_rank[ip].this_int = mpi_rank;
    d2_rank[ip].this_double = 1.0e+20;
  }

  for (int ibf = 0; ibf < nbf; ++ibf) {
    for (int ip=0; ip<np; ++ip) {
      double this_d2 = DIST2(xp_seed[ip],x_bf[ibf]);

      if ( this_d2 < d2_rank[ip].this_double) {
        d2_rank[ip].this_double = this_d2;
        ibf_zone[ip]            = ibf; // rank local ibf is stored right now...
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE,d2_rank,np,MPI_DOUBLE_INT,MPI_MINLOC,mpi_comm);

  for (int ip=0; ip<np; ++ip) {
    if ( mpi_rank == d2_rank[ip].this_int) {

      FOR_I3 x_bf_seed[ip][i] = x_bf[ibf_zone[ip]][i];

      izone[ip] = zone_bf[ibf_zone[ip]];

      // return the ibf_zone to a zone local indexing
      ibf_zone[ip] -= bfZoneVec[izone[ip]].ibf_f;
      assert( (ibf_zone[ip] >= 0) && (ibf_zone[ip] < bfZoneVec[izone[ip]].nbf));

    }
    else {
      FOR_I3 x_bf_seed[ip][i] = 0.0;  // can reduce via summation
      ibf_zone[ip] = -d2_rank[ip].this_int-1; // you dont own the seed, -1 index of which rank does
      izone[ip] = -1;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,(double*)x_bf_seed,3*np,MPI_DOUBLE,MPI_SUM,mpi_comm);
  MPI_Allreduce(MPI_IN_PLACE,izone,np,MPI_INT,MPI_MAX,mpi_comm);  // non-owners are -1, so max should be only owner
  // check validitiy
  for (int ip=0, nzn=bfZoneVec.size(); ip<np; ++ip) {
    assert( (izone[ip] >= 0) && (izone[ip] < nzn));
  }
  DELETE(d2_rank);
}

void StaticSolver::findClosestBfOnZones(double (*x_bf_seed)[3], int * ibf_zone, int * izone, const double (*xp_seed)[3],const int np,const set<int>& zonesVec) const {

  DoubleInt * d2_rank = new DoubleInt[np];

  bool * bfz_flag = new bool[bfZoneVec.size()];
  for (int izn=0,nzn=bfZoneVec.size(); izn<nzn; ++izn) bfz_flag[izn] = 0;
  for (set<int>::const_iterator it=zonesVec.begin(); it!=zonesVec.end(); ++it) {
    bfz_flag[*it] = true;
  }

  for (int ip=0; ip<np; ++ip) {
    ibf_zone[ip] = -1;
    izone[ip] = -1;
    d2_rank[ip].this_int = mpi_rank;
    d2_rank[ip].this_double = 1.0e+20;
  }

  for (int ibf = 0; ibf < nbf; ++ibf) {
    if (bfz_flag[zone_bf[ibf]]) {
      for (int ip=0; ip<np; ++ip) {
        double this_d2 = DIST2(xp_seed[ip],x_bf[ibf]);

        if ( this_d2 < d2_rank[ip].this_double) {
          d2_rank[ip].this_double = this_d2;
          ibf_zone[ip]            = ibf; // rank local ibf is stored right now...
        }
      }
    }
  }
  DELETE(bfz_flag);
  MPI_Allreduce(MPI_IN_PLACE,d2_rank,np,MPI_DOUBLE_INT,MPI_MINLOC,mpi_comm);

  for (int ip=0; ip<np; ++ip) {
    if ( mpi_rank == d2_rank[ip].this_int) {

      FOR_I3 x_bf_seed[ip][i] = x_bf[ibf_zone[ip]][i];

      izone[ip] = zone_bf[ibf_zone[ip]];

      // return the ibf_zone to a zone local indexing
      ibf_zone[ip] -= bfZoneVec[izone[ip]].ibf_f;
      assert( (ibf_zone[ip] >= 0) && (ibf_zone[ip] < bfZoneVec[izone[ip]].nbf));

    }
    else {
      FOR_I3 x_bf_seed[ip][i] = 0.0;  // can reduce via summation
      ibf_zone[ip] = -d2_rank[ip].this_int-1; // you dont own the seed, -1 index of which rank does
      izone[ip] = -1;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,(double*)x_bf_seed,3*np,MPI_DOUBLE,MPI_SUM,mpi_comm);
  MPI_Allreduce(MPI_IN_PLACE,izone,np,MPI_INT,MPI_MAX,mpi_comm);  // non-owners are -1, so max should be only owner
  // check validitiy
  for (int ip=0, nzn=bfZoneVec.size(); ip<np; ++ip) {
    assert( (izone[ip] >= 0) && (izone[ip] < nzn));
  }
  DELETE(d2_rank);
}

void StaticSolver::processMultiFluxProbe(Param* param,const int step,const double time,const bool killfile_request) {
  // use the full param as the hash for the map..
  if ( mfpMap.find(param->str()) == mfpMap.end()) {

    // couldnt find this probe so we need to insert it to the list
    pair<map<const string,MultiFluxProbe>::iterator,bool> ret = mfpMap.insert(pair<const string,MultiFluxProbe>(param->str(), MultiFluxProbe()));
    assert( ret.second);

    const int ierr = initMultiFluxProbe(ret.first->second,param);

    if ( ierr != 0) {

      // see the note on the error handling logic above for the PointProbe
      if ( !killfile_request ) {
        CERR( " > unable to add MultiFluxProbe from: \n >   " << param->str() );
      } else {
        CWARN( " > unable to add MultiFluxProbe from: \n >  " << param->str() );
        mfpMap.erase(ret.first);
      }

    } else {
      // if the interval is invalid, then its one-and-done for this probe
      // otherwise register for continuous processing
      if ( ret.first->second.interval == -1) {
        ret.first->second.doProbe(step,time);
        ret.first->second.flush();
        mfpMap.erase(ret.first);
      }
    }
  }
}

int StaticSolver::initMultiFluxProbe(MultiFluxProbe& mfp,Param * param) {

  // use a static counter to round-robin assignment of probe rank0...

  static int mfp_rank0 = 0;

  mfp.probe_rank0 = mfp_rank0++;
  if (mfp_rank0 == mpi_size) mfp_rank0 = 0;

  std::stringstream hss; // header-stringstream

  // start the header with the param line and the hash...

  if (mpi_rank == mfp.probe_rank0) {
    hss << setprecision(8);
    hss << "# " << param->str() << endl;
    hss << "# mles hash id " << RestartHashUtilities::mlesHash << endl;
    hss << "# sles hash id " << RestartHashUtilities::slesHash << endl;
  }

  bool got_name = false;
  double xp[3],np[3],dn = HUGE_VAL;
  bool got_xp = false;
  bool got_np = false;
  bool got_dn = false;
  //bool got_write_images = false;
  //int image_nx = 800;
  //int image_ny = 800;
  vector<string> varNameVec;

  int iarg = 0;
  while (iarg < param->size()) {
    string token = MiscUtils::toUpperCase(param->getString(iarg++));
    if (token == "NAME") {
      mfp.name = param->getString(iarg++);
      got_name = true;
    }
    else if (token == "INTERVAL" || token == "WRITE_INTERVAL") {
      mfp.interval = param->getInt(iarg++); // default is -1
      if (mfp.interval <= 0) {
        CWARN(" > INTERVAL expects a positive integer; treating as one-off request");
        // invalid value entered, so treat as unset (default value)
        mfp.interval = -1;
      }
    }
    else if (token == "FORMAT") {
      token = param->getString(iarg++);
      if (token == "SINGLE_FILE") {
        mfp.format = SINGLE_FILE;
      }
      else if (token == "MULTI_FILE") {
        mfp.format = MULTI_FILE;
      }
      else {
        CWARN("unrecognized MULTIFLUX_PROBE FORMAT: " << token << ". Skipping.");
      }
    }
    else if (token == "XP") {
      FOR_I3 xp[i] = param->getDouble(iarg++);
      got_xp = true;
    }
    else if (token == "NP") {
      FOR_I3 np[i] = param->getDouble(iarg++);
      // turn np into a unit vector...
      const double mag = MAG(np);
      assert(mag > 0.0);
      FOR_I3 np[i] /= mag;
      got_np = true;
    }
    else if (token == "DN") {
      dn = param->getDouble(iarg++);
      got_dn = true;
    }
    // TODO used to be able to request image of the flux probe -- add this later...
    /*
      else if (token == "WRITE_IMAGES") {
      got_write_images = true;
      }
      else if ((token == "SIZE")||(token == "IMAGE_SIZE")) {
      image_nx = param->getInt(iarg++);
      image_ny = param->getInt(iarg++);
      }
    */
    else if (token == "VARS") {
      // for vars, don't push directly into mfp's varNameVec, because
      // we want to make sure they eval to FA_DN...
      while (iarg < param->size())
        varNameVec.push_back(param->getString(iarg++));
    }
    else {
      if (mpi_rank == 0)
        cout << "processFluxProbe: skipping unrecognized MULTIFLUX_PROBE token \"" << token << "\"" << endl;
    }
  }

  int ierr = 0;
  if (!got_name) {
    if (mpi_rank == 0)
      cout << "Warning: MULTIFLUX_PROBE missing NAME." << endl;
    ierr = -1;
  }
  if (!got_xp) {
    if (mpi_rank == 0)
      cout << "Warning: MULTIFLUX_PROBE missing XP." << endl;
    ierr = -1;
  }
  if (!got_np) {
    if (mpi_rank == 0)
      cout << "Warning: MULTIFLUX_PROBE missing NP." << endl;
    ierr = -1;
  }
  if (ierr != 0) {
    if (mpi_rank == 0)
      cout << "Warning: Skipping this MULTIFLUX_PROBE. Correct syntax warnings." << endl;
    return -1;
  }

  if (got_dn) {
    // adjust xp by dn in the unit normal direction...
    FOR_I3 xp[i] += dn*np[i];
  }

  if (mpi_rank == 0)
    cout << "Building MULTIFLUX_PROBE \"" << mfp.name << "\" xp=" << COUT_VEC(xp) << " np=" << COUT_VEC(np) << "..." << endl;

  // check that the data evaluates to FA_DN. If is doesn't, discard and let user know...

  assert(mfp.varNameVec.empty());
  for (int ii = 0,ii_end = varNameVec.size(); ii < ii_end; ++ii) {
    CtiRegister::CtiData * data = CtiRegister::getCtiData(varNameVec[ii]);
    if ((data == NULL)||(data->getTopology() != SIGNED_FA_DATA)||(data->getType() != DN_DATA)) {
      CWARN("MultiFluxProbe data \"" << varNameVec[ii] << "\" does not evaluate to SIGNED_FA_DN. Skipping.");
    }
    else {
      mfp.varNameVec.push_back(varNameVec[ii]);
    }
  }

  if (mfp.varNameVec.empty()) {
    CWARN("MultiFluxProbe has no valid data. Skipping.");
    return -1;
  }

  int * cv_flag     = new int[ncv_g];
  for (int icv = 0; icv < ncv; ++icv) {
    const double dn =
      (x_cv[icv][0]-xp[0])*np[0] +
      (x_cv[icv][1]-xp[1])*np[1] +
      (x_cv[icv][2]-xp[2])*np[2];
    if (dn >= 0.0) {
      cv_flag[icv] = 1;
    }
    else {
      cv_flag[icv] = -1;
    }
  }

  updateCvData(cv_flag);

  for (int icv = 0; icv < ncv; ++icv) {
    for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      if (cv_flag[icv]*cv_flag[icv_nbr] < 0) {
        // an odd product means this cv is on the boundary...
        cv_flag[icv] *= 2; // makes -1 -> -2
        assert((cv_flag[icv] == 2)||(cv_flag[icv] == -2));
        break;
      }
    }
  }

  updateCvData(cv_flag);

  for (int icv = 0; icv < ncv; ++icv) {
    if ((cv_flag[icv] == 1)||(cv_flag[icv] == -1)) {
      for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        if ((cv_flag[icv_nbr] == 2)||(cv_flag[icv_nbr] == -2)) {
          cv_flag[icv] *= 4;
          assert((cv_flag[icv] == 4)||(cv_flag[icv] == -4));
          break;
        }
      }
    }
  }

  updateCvData(cv_flag);

  // now we have the zone of cvs around flux faces 2-deep all flagged with
  // +/- 2 or +/- 4. Hopefully 2-deep is sufficient to keep cv groups
  // together. Now re-index these cvs...

  int my_ncv_active = 0;
  for (int icv = 0; icv < ncv; ++icv) {
    if ((cv_flag[icv] > 1)||(cv_flag[icv] < -1)) {
      cv_flag[icv] = my_ncv_active++;
    }
    else {
      cv_flag[icv] = -1;
    }
  }

  int * acora = NULL;
  MiscUtils::buildXora(acora,my_ncv_active);
  assert(acora[mpi_size] < TWO_BILLION);

  for (int icv = 0; icv < ncv; ++icv) {
    if (cv_flag[icv] >= 0) {
      cv_flag[icv] += acora[mpi_rank];
    }
  }

  if (mpi_rank == 0)
    cout << " > group";

  int done = 0;
  while (done == 0) {

    if (mpi_rank == 0) {
      cout << ".";
      cout.flush();
    }

    int my_done = 1;

    updateCvData(cv_flag);

    for (int icv = 0; icv < ncv; ++icv) {
      if (cv_flag[icv] >= 0) {
        for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          if ((cv_flag[icv_nbr] >= 0)&&(cv_flag[icv_nbr] < cv_flag[icv])) {
            cv_flag[icv] = cv_flag[icv_nbr];
            my_done = 0;
          }
        }
      }
    }

    MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);

  }

  if (mpi_rank == 0)
    cout << "OK" << endl;

  // now we need to compress...

  int * active_flag = new int[my_ncv_active];
  for (int icv_active = 0; icv_active < my_ncv_active; ++icv_active)
    active_flag[icv_active] = 0;

  for (int icv = 0; icv < ncv; ++icv) {
    if ((cv_flag[icv] >= acora[mpi_rank])&&(cv_flag[icv] < acora[mpi_rank+1])) {
      const int icv_active = cv_flag[icv] - acora[mpi_rank];
      active_flag[icv_active] = 1;
    }
  }

  int my_ngr = 0;
  for (int icv_active = 0; icv_active < my_ncv_active; ++icv_active) {
    if (active_flag[icv_active] == 1)
      active_flag[icv_active] = my_ngr++;
    else
      active_flag[icv_active] = -1;
  }

  int * grora = NULL;
  MiscUtils::buildXora(grora,my_ngr);

  for (int icv = 0; icv < ncv; ++icv) {
    if ((cv_flag[icv] >= acora[mpi_rank])&&(cv_flag[icv] < acora[mpi_rank+1])) {
      const int icv_active = cv_flag[icv] - acora[mpi_rank];
      assert(active_flag[icv_active] >= 0);
      cv_flag[icv] = active_flag[icv_active] + grora[mpi_rank];
    }
  }
  delete[] active_flag;
  delete[] acora;

  // now do reduction again to get final coloring...

  if (mpi_rank == 0)
    cout << " > reduce";

  done = 0;
  while (done == 0) {

    if (mpi_rank == 0) {
      cout << ".";
      cout.flush();
    }

    int my_done = 1;

    updateCvData(cv_flag);

    for (int icv = 0; icv < ncv; ++icv) {
      if (cv_flag[icv] >= 0) {
        for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          if ((cv_flag[icv_nbr] >= 0)&&(cv_flag[icv_nbr] < cv_flag[icv])) {
            cv_flag[icv] = cv_flag[icv_nbr];
            my_done = 0;
          }
        }
      }
    }

    MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);

  }

  if (mpi_rank == 0)
    cout << "OK" << endl;

  // although it is not strictly necessary, if the number of flux probes
  // exceeds some large number, like 1000s, we might consider an even more
  // distributed way of managing them. For now, they all are gathered/reduced
  // to a single rank...

  mfp.probe_ngr = grora[mpi_size];
  delete[] grora;

  if (mpi_rank == 0)
    cout << " > flux probe count: " << mfp.probe_ngr << endl;

  // switch the group indexing to -2 indexing on the opposite side...

  for (int icv = 0; icv < ncv; ++icv) {
    if (cv_flag[icv] >= 0) {
      const double dn =
        (x_cv[icv][0]-xp[0])*np[0] +
        (x_cv[icv][1]-xp[1])*np[1] +
        (x_cv[icv][2]-xp[2])*np[2];
      if (dn < 0.0) {
        cv_flag[icv] = -cv_flag[icv]-2; // use -2 indexing...
      }
    }
    else {
      assert(cv_flag[icv] == -1);
    }
  }

  updateCvData(cv_flag);

  // now build groups of (group,face) pairs and (group,izone) set...

  set<int> groupSet;
  vector<pair<int,int> > groupFaceVec;
  set<pair<int,int> > groupZoneSet;

  for (int icv = 0; icv < ncv; ++icv) {
    if (cv_flag[icv] != -1) {
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        const int icv0 = cvofa[ifa][0];
        const int icv1 = cvofa[ifa][1];
        assert((icv0 == icv)||(icv1 == icv));
        if (cv_flag[icv0] == -cv_flag[icv1]-2) {
          const int igr = max(cv_flag[icv],-cv_flag[icv]-2);
          assert((igr >= 0)&&(igr < mfp.probe_ngr));
          // this is a flux face. Depending on its orientation, it
          // may or may not be included in the flux calculation. Any
          // boundaries this flux face touches, however, should be
          // recorded...
          if (cv_flag[icv] < 0) {
            assert(igr == -cv_flag[icv]-2);
            groupSet.insert(igr);
            // by convention, we are going to add the flux face from the
            // cv on the negative side. That is this cv. The direction
            // of the face, however, could be -ve or positive depending
            // on its orientation wrt this cv...
            if (icv == icv0) {
              // face is outward pointing...
              groupFaceVec.push_back(pair<int,int>(igr,ifa));
            }
            else {
              assert(icv == icv1);
              // face is inward pointing...
              groupFaceVec.push_back(pair<int,int>(igr,-ifa-1));
            }
          }
          for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
            const int ibf = bfocv_v[boc];
            const int izone = zone_bf[ibf];
            groupZoneSet.insert(pair<int,int>(igr,izone));
          }
        }
      }
    }
  }

  delete[] cv_flag;

  // reduce groupZoneSet info to probe rank0...

  vector<string> mfp_zoneVec;

  int my_count = groupZoneSet.size();
  int * count = NULL;
  if (mpi_rank == mfp.probe_rank0)
    count = new int[mpi_size];

  MPI_Gather(&my_count,1,MPI_INT,count,1,MPI_INT,mfp.probe_rank0,mpi_comm);

  // add everyone's set to our set...
  if (mpi_rank == mfp.probe_rank0) {
    int ibuf_size = 0;
    FOR_RANK {
      if (rank != mpi_rank) {
        ibuf_size = max(ibuf_size,count[rank]);
      }
    }
    int * ibuf = new int[ibuf_size*2];
    FOR_RANK {
      if ((rank != mpi_rank)&&(count[rank] > 0)) {
        MPI_Recv(ibuf,count[rank]*2,MPI_INT,rank,42424,mpi_comm,MPI_STATUS_IGNORE);
        for (int i = 0; i < count[rank]; ++i) {
          const int igr = ibuf[2*i]; assert((igr >= 0)&&(igr < mfp.probe_ngr));
          const int izone = ibuf[2*i+1]; assert((izone >= 0)&&(izone < int(bfZoneVec.size())));
          groupZoneSet.insert(pair<int,int>(igr,izone));
        }
      }
    }
    delete[] ibuf;

    // and build the probe's zoneVec...

    mfp_zoneVec.resize(mfp.probe_ngr);

    int igr_current = -1;
    for (set<pair<int,int> >::const_iterator iter = groupZoneSet.begin(); iter != groupZoneSet.end(); ++iter) {
      if (iter->first != igr_current) {
        assert(iter->first > igr_current); // sets are sorted
        igr_current = iter->first;
        mfp_zoneVec[igr_current] = bfZoneVec[iter->second].getName();
      }
      else {
        mfp_zoneVec[igr_current].append(":"+bfZoneVec[iter->second].getName());
      }
    }
    groupZoneSet.clear();

  }
  else if (my_count > 0) {

    // non-rank0 probe with zone info...

    int * ibuf = new int[my_count*2];
    int i = 0;
    for (set<pair<int,int> >::const_iterator iter = groupZoneSet.begin(); iter != groupZoneSet.end(); ++iter) {
      ibuf[i++] = iter->first;
      ibuf[i++] = iter->second;
    }
    assert(i == my_count*2);
    MPI_Ssend(ibuf,my_count*2,MPI_INT,mfp.probe_rank0,42424,mpi_comm);
    delete[] ibuf;
    groupZoneSet.clear();

  }

  // ====================================================
  // now the internal face pack/unpack info...
  // ====================================================

  mfp.probe_pack_buf_size = groupSet.size();
  int * ibuf = new int[mfp.probe_pack_buf_size];

  mfp.probe_pack_index_size =  groupFaceVec.size();
  assert(mfp.probe_pack_index == NULL);
  mfp.probe_pack_index = new int[mfp.probe_pack_index_size][2]; // 2 ints: pack index, and -1-indexed flux face...

  sort(groupFaceVec.begin(),groupFaceVec.end());
  int ibuf_pos = -1;
  int igr_current = -1;
  for (int ii = 0,ii_end = groupFaceVec.size(); ii < ii_end; ++ii) {
    const int igr = groupFaceVec[ii].first;
    if (igr != igr_current) {
      assert(igr > igr_current);
      igr_current = igr;
      ++ibuf_pos;
      ibuf[ibuf_pos] = igr_current;
    }
    mfp.probe_pack_index[ii][0] = ibuf_pos;
    mfp.probe_pack_index[ii][1] = groupFaceVec[ii].second; // the face
  }
  assert(ibuf_pos+1 == mfp.probe_pack_buf_size);
  groupFaceVec.clear();

  MPI_Gather(&mfp.probe_pack_buf_size,1,MPI_INT,count,1,MPI_INT,mfp.probe_rank0,mpi_comm);

  if (mpi_rank == mfp.probe_rank0) {

    // note we include our own count in this buf because it makes the
    // reduction much cleaner...
    int ioar_s = 0; // index-of-rank-size
    FOR_RANK {
      mfp.probe_unpack_buf_size = max(mfp.probe_unpack_buf_size,count[rank]);
      if (count[rank] > 0) {
        ++mfp.probe_active_rank_count;
        ioar_s += count[rank];
      }
    }

    mfp.probe_raoar = new int[mfp.probe_active_rank_count]; // rank-of-active-rank
    mfp.probe_ioar_i = new int[mfp.probe_active_rank_count+1];
    mfp.probe_ioar_i[0] = 0;
    mfp.probe_ioar_v = new int[ioar_s];

    int iar = 0;
    int ioa = 0;
    FOR_RANK {
      if (count[rank] > 0) {
        mfp.probe_raoar[iar] = rank;
        mfp.probe_ioar_i[iar+1] = mfp.probe_ioar_i[iar] + count[rank];
        ++iar;
        if (rank == mpi_rank) {
          assert(count[rank] == mfp.probe_pack_buf_size);
          for (int i = 0; i < mfp.probe_pack_buf_size; ++i)
            mfp.probe_ioar_v[ioa++] = ibuf[i];
        }
        else {
          MPI_Recv(mfp.probe_ioar_v+ioa,count[rank],MPI_INT,rank,42425,mpi_comm,MPI_STATUS_IGNORE);
          ioa += count[rank];
        }
      }
    }
    assert(iar == mfp.probe_active_rank_count);
    assert(ioa == ioar_s);

  }
  else if (mfp.probe_pack_buf_size > 0) {

    MPI_Ssend(ibuf,mfp.probe_pack_buf_size,MPI_INT,mfp.probe_rank0,42425,mpi_comm);

  }

  delete[] ibuf;
  MPI_Barrier(mpi_comm);

  // cleanup...

  if (mpi_rank == mfp.probe_rank0)
    delete[] count;

  // ---------------------------------------------------
  // now all the structures are built to support exchanges.
  // start with reducing the projected area and the centroid
  // of the flux zones...
  // ---------------------------------------------------

  assert(mfp.probe_pack_buf_size <= mfp.probe_pack_index_size);

  double * pack_buf = new double[mfp.probe_pack_buf_size*4];
  for (int ipack = 0; ipack < mfp.probe_pack_buf_size*4; ++ipack)
    pack_buf[ipack] = 0.0;

  double my_proj_area_min = 0.0;
  for (int ii = 0; ii < mfp.probe_pack_index_size; ++ii) {
    const int ipack = mfp.probe_pack_index[ii][0];
    assert((ipack >= 0)&&(ipack < mfp.probe_pack_buf_size));
    int ifa = mfp.probe_pack_index[ii][1];
    double fa_sign = 1.0;
    if (ifa < 0) {
      ifa = -ifa-1;
      fa_sign = -1.0;
    }
    assert((ifa >= 0)&&(ifa < nfa));

    const double proj_area = fa_sign*DOT_PRODUCT(n_fa[ifa],np);

    // the projectd area should not be less than 0.0
    my_proj_area_min = min(my_proj_area_min,proj_area);

    pack_buf[ipack*4  ] += proj_area*x_fa[ifa][0];
    pack_buf[ipack*4+1] += proj_area*x_fa[ifa][1];
    pack_buf[ipack*4+2] += proj_area*x_fa[ifa][2];
    pack_buf[ipack*4+3] += proj_area;
  }

  /*
    for (int ipack = 0; ipack < mfp.probe_pack_buf_size; ++ipack) {
    cout << " ipack: " << ipack << " proj_area: " << pack_buf[ipack*4+3] << " centroid: " <<
    pack_buf[ipack*4  ]/pack_buf[ipack*4+3] << " " <<
    pack_buf[ipack*4+1]/pack_buf[ipack*4+3] << " " <<
    pack_buf[ipack*4+2]/pack_buf[ipack*4+3] << endl;
    }
  */

  double * unpack_buf = NULL;
  if (mpi_rank == mfp.probe_rank0)
    unpack_buf = new double[mfp.probe_ngr*4];

  mfp.probe_reduce(pack_buf,unpack_buf,4);

  delete[] pack_buf;

  // this is an un-normalized check -- may fail...

  //assert(my_proj_area_min > -1.0E-12);

  // now report...

  if (mpi_rank == mfp.probe_rank0) {

    {

      list<pair<string,pair<double,int> > > sortList; // sort by name, d2, original order

      for (int igr = 0; igr < mfp.probe_ngr; ++igr) {
        const double proj_area = unpack_buf[igr*4+3];
        assert(proj_area > 0.0);
        const double xc[3] = {
          unpack_buf[igr*4  ]/proj_area,
          unpack_buf[igr*4+1]/proj_area,
          unpack_buf[igr*4+2]/proj_area
        };
        const double d2 = DIST2(xc,xp);

        /*
          cout << " probe: " << igr << " zones: \"" << mfp_zoneVec[igr] << "\", projected area: " <<
          proj_area << ", centroid: " << COUT_VEC(xc) << ", dist: " << sqrt(d2) << endl;
        */

        // set the sortVec elements...

        // NOTE: if you want to change the way that probes are sorted, change this double d2...

        sortList.push_back(pair<string,pair<double,int> >(mfp_zoneVec[igr],pair<double,int>(d2,igr)));

      }

      sortList.sort();

      int * new_order = new int[mfp.probe_ngr];
      int ngr_new = 0;
      int supergroup_ngr = 0;
      assert(mfp.probe_supergroup_last == NULL); mfp.probe_supergroup_last = new int[mfp.probe_ngr];
      assert(mfp.probe_proj_area == NULL); mfp.probe_proj_area = new double[mfp.probe_ngr];
      assert(mfp.probe_x == NULL); mfp.probe_x = new double[mfp.probe_ngr][3];

      list<pair<string,pair<double,int> > >::const_iterator iter_prev;
      for (list<pair<string,pair<double,int> > >::const_iterator iter = sortList.begin(); iter != sortList.end(); ++iter) {
        if (ngr_new > 0) {
          if (iter->first != iter_prev->first) {
            mfp.probe_supergroup_last[ngr_new-1] = supergroup_ngr;
            supergroup_ngr = 0;
          }
          else {
            mfp.probe_supergroup_last[ngr_new-1] = 0;
          }
        }
        // this is needed for reordering the unpack stuff...
        const int igr = iter->second.second;
        new_order[igr] = ngr_new;
        // set the name, proj area and x...
        mfp_zoneVec[ngr_new] = iter->first;
        mfp.probe_proj_area[ngr_new] = unpack_buf[igr*4+3];
        mfp.probe_x[ngr_new][0] = unpack_buf[igr*4  ]/mfp.probe_proj_area[ngr_new];
        mfp.probe_x[ngr_new][1] = unpack_buf[igr*4+1]/mfp.probe_proj_area[ngr_new];
        mfp.probe_x[ngr_new][2] = unpack_buf[igr*4+2]/mfp.probe_proj_area[ngr_new];
        // increment the stuff...
        ++supergroup_ngr;
        ++ngr_new;
        iter_prev = iter;
      }
      // complete mfp.probe_supergroup_last to support reporting...
      assert(ngr_new == mfp.probe_ngr);
      if (ngr_new > 0) {
        mfp.probe_supergroup_last[ngr_new-1] = supergroup_ngr;
      }

      delete[] unpack_buf;

      // all we need to update is the local index in mfp.probe_ioar_v into which we put the reduced
      // flux results.

      const int mfp_ioar_s = mfp.probe_ioar_i[mfp.probe_active_rank_count];
      for (int ioa = 0; ioa < mfp_ioar_s; ++ioa) {
        const int igr = mfp.probe_ioar_v[ioa]; assert((igr >= 0)&&(igr < mfp.probe_ngr));
        mfp.probe_ioar_v[ioa] = new_order[igr];
      }

      delete[] new_order;

    }

    // ==========================================
    // report probes...
    // ==========================================

    int supergroup_index = 0;
    int supergroup_ngr = 0;
    double total_buf[4] = { 0.0, 0.0, 0.0, 0.0 };
    double supergroup_buf[4] = { 0.0, 0.0, 0.0, 0.0 }; // sum(xc*proj_area), sum(proj_area)
    for (int igr = 0; igr < mfp.probe_ngr; ++igr) {
      // add supergroup and total bufs...
      supergroup_buf[0] += mfp.probe_x[igr][0]*mfp.probe_proj_area[igr];
      supergroup_buf[1] += mfp.probe_x[igr][1]*mfp.probe_proj_area[igr];
      supergroup_buf[2] += mfp.probe_x[igr][2]*mfp.probe_proj_area[igr];
      supergroup_buf[3] += mfp.probe_proj_area[igr];
      total_buf[0] += mfp.probe_x[igr][0]*mfp.probe_proj_area[igr];
      total_buf[1] += mfp.probe_x[igr][1]*mfp.probe_proj_area[igr];
      total_buf[2] += mfp.probe_x[igr][2]*mfp.probe_proj_area[igr];
      total_buf[3] += mfp.probe_proj_area[igr];
      // and report...
      if (mfp.probe_ngr == 1) {
        hss << "# probe: count: 1 zones: \"" << mfp_zoneVec[igr] << "\"" << endl;
      }
      else if (mfp.probe_supergroup_last[igr] == 1) {
        assert(supergroup_ngr == 0);
        // this is a probe that is the ONLY member of the current supergroup. report it separately with single indexing...
        // unless it is the only group. In this case, we will only report the total...
        //if (mfp.probe_ngr > 1) {
        /*
          hss << "# probe-" << (supergroup_index+1) << ": zones: \"" << mfp_zoneVec[igr] << "\", xc: " <<
          mfp.probe_x[igr][0] << " " << mfp.probe_x[igr][1] << " " << mfp.probe_x[igr][2] << ", proj_area: " << mfp.probe_proj_area[igr] << endl;
          hss << "# probe-sum-" << (supergroup_index+1) << ": count: 1 zones: \"" << mfp_zoneVec[igr] << "\", xc: " <<
          mfp.probe_x[igr][0] << " " << mfp.probe_x[igr][1] << " " << mfp.probe_x[igr][2] << ", proj_area: " << mfp.probe_proj_area[igr] << endl;
        */
        hss << "# probe-sum-" << (supergroup_index+1) << ": count: 1 zones: \"" << mfp_zoneVec[igr] << "\"" << endl;
        //}
        ++supergroup_index;
        // leave supergroup_ngr = 0
        for (int i = 0; i < 4; ++i) supergroup_buf[i] = 0.0;
      }
      else {
        // this is a probe that is part of the current supergroup. report it separately with double indexing...
        //hss << "# probe-" << (supergroup_index+1) << "-" << (supergroup_ngr+1) << ": zones: \"" << mfp_zoneVec[igr] << "\", xc: " <<
        //  mfp.probe_x[igr][0] << " " << mfp.probe_x[igr][1] << " " << mfp.probe_x[igr][2] << ", proj_area: " << mfp.probe_proj_area[igr] << endl;
        // if this is the last one, then report the supergroup as well...
        if (mfp.probe_supergroup_last[igr] > 1) {
          //hss << "# probe-sum-" << (supergroup_index+1) << ": count: " << mfp.probe_supergroup_last[igr] << " zones: \"" << mfp_zoneVec[igr] << "\", xc: " <<
          //  supergroup_buf[0]/supergroup_buf[3] << " " << supergroup_buf[1]/supergroup_buf[3] << " " << supergroup_buf[2]/supergroup_buf[3] <<
          //  ", proj_area: " << supergroup_buf[3] << endl;
          hss << "# probe-sum-" << (supergroup_index+1) << ": count: " << mfp.probe_supergroup_last[igr] << " zones: \"" << mfp_zoneVec[igr] << "\"" << endl;
          ++supergroup_index;
          supergroup_ngr = 0;
          for (int i = 0; i < 4; ++i) supergroup_buf[i] = 0.0;
        }
        else {
          ++supergroup_ngr;
        }
      }
    }
    // and finally report the total if there is more than 1 probe...
    if (mfp.probe_ngr > 1)
      hss << "# probe-total: count: " << mfp.probe_ngr << endl;
    //" xc: " << total_buf[0]/total_buf[3] << " " << total_buf[1]/total_buf[3] << " " << total_buf[2]/total_buf[3] <<
    //", proj_area: " << total_buf[3] << endl;
  }

  // =========================================================
  // finally image writing
  // =========================================================

#ifdef SKIP_IMAGE_WRITING

  if (got_write_images) {

    if (mpi_rank == 0)
      cout << " > WRITE_IMAGES..." << endl;

    // do a reverse send of the sorted group information to the probe on the pack side...

    int * pack_ibuf = new int[mfp.probe_pack_buf_size];

    if (mpi_rank == mfp.probe_rank0) {

      int * ibuf = new int[mfp.probe_unpack_buf_size];
      for (int iar = 0; iar < mfp.probe_active_rank_count; ++iar) {
        const int rank = mfp.probe_raoar[iar]; assert((rank >= 0)&&(rank < mpi_size));
        if (rank == mpi_rank) {
          for (int ioa = mfp.probe_ioar_i[iar]; ioa != mfp.probe_ioar_i[iar+1]; ++ioa) {
            const int igr = mfp.probe_ioar_v[ioa]; assert((igr >= 0)&&(igr < mfp.probe_ngr));
            pack_ibuf[ioa-mfp.probe_ioar_i[iar]] = igr; // note that this is the new group index...
          }
        }
        else {
          for (int ioa = mfp.probe_ioar_i[iar]; ioa != mfp.probe_ioar_i[iar+1]; ++ioa) {
            const int igr = mfp.probe_ioar_v[ioa]; assert((igr >= 0)&&(igr < mfp.probe_ngr));
            ibuf[ioa-mfp.probe_ioar_i[iar]] = igr;
          }
          MPI_Ssend(ibuf,mfp.probe_ioar_i[iar+1]-mfp.probe_ioar_i[iar],MPI_INT,rank,42427,mpi_comm);
        }
      }
      delete[] ibuf;

    }
    else if (mfp.probe_pack_buf_size > 0) {

      MPI_Recv(pack_ibuf,mfp.probe_pack_buf_size,MPI_INT,mfp.probe_rank0,42427,mpi_comm,MPI_STATUS_IGNORE);

    }

    // now the pack_ibuf has the group index associated with each face...

    // start by getting the width...

    double my_d2_max = 0.0;
    for (int ii = 0; ii < mfp.probe_pack_index_size; ++ii) {
      int ifa = mfp.probe_pack_index[ii][1];
      if (ifa < 0) {
        ifa = -ifa-1;
      }
      assert((ifa >= 0)&&(ifa < nfa));
      const double d2 = DIST2(x_fa[ifa],xp);
      my_d2_max = max(my_d2_max,d2);
    }
    double d2_max;
    MPI_Allreduce(&my_d2_max,&d2_max,1,MPI_DOUBLE,MPI_MAX,mpi_comm);

    // Note: do this more efficiently in the future with a single reduction of the
    // surface image, then tiny reductions of the red regions. But for now do it the
    // hard way...

    int * stora = NULL;
    buildUniformXora(stora,surface->nst);

    // also put r_cv, the cv radius of the voronoi cv, into the ghost data...
    double * r_cv_with_ghost = new double[ncv_g];
    for (int icv = 0; icv < ncv; ++icv)
      r_cv_with_ghost[icv] = r_cv[icv];
    updateCvData(r_cv_with_ghost);

    for (int igr = 0; igr < mfp.probe_ngr; ++igr) {

      if (mpi_rank == 0)
        cout << " > Imaging fluxprobe group: " << igr << endl;

      CtiCanvas canvas(xp,np,2.1*sqrt(d2_max),image_nx,image_ny);

      canvas.addSurfaceTris(stora[mpi_rank+1]-stora[mpi_rank],surface->spost+stora[mpi_rank],surface->xp);

      // now the faces...

      canvas.activateRgbFlood(); // this must be done collective

      for (int ii = 0; ii < mfp.probe_pack_index_size; ++ii) {
        const int ipack = mfp.probe_pack_index[ii][0];
        if (pack_ibuf[ipack] == igr) {
          const int ifa = max(mfp.probe_pack_index[ii][1],-mfp.probe_pack_index[ii][1]-1);
          assert((ifa >= 0)&&(ifa < nfa));
          // the distance (squared) between cvs...
          const double dist2 = DIST2(x_cv[cvofa[ifa][0]],x_cv[cvofa[ifa][1]]);
          // the minimum radius of the 2 cvs...
          const double r_min = min(r_cv_with_ghost[cvofa[ifa][0]],r_cv_with_ghost[cvofa[ifa][1]]);
          // so the face radius must be less than this...
          const double r2_fa = r_min*r_min - 0.25*dist2; assert(r2_fa > 0.0);
          double r_fa = sqrt(r2_fa);
          // we need to also use the midpoint between cvs...
          double x_mid[3]; FOR_I3 x_mid[i] = 0.5*(x_cv[cvofa[ifa][0]][i] + x_cv[cvofa[ifa][1]][i]);
          canvas.addRgbFlood(220,0,0,&x_mid,&r_fa,1);
          canvas.forceRgbPixels(255,0,0,&fa[ifa].x,1);
        }
      }

      char filename[256];
      sprintf(filename,"%s-%04d.png",name.c_str(),igr);
      canvas.writeImage(filename);

    }

    delete[] stora;
    delete[] r_cv_with_ghost;
    delete[] pack_ibuf;

  }

#endif

  // final prep...

  if (mpi_rank == mfp.probe_rank0) {

    // finish the header with the variable line...

    hss << "# 1:name 2:step 3:time 4:xp 5:yp 6:zp 7:proj_area";
    for (int i = 0, i_end = mfp.varNameVec.size(); i < i_end; ++i)
      hss << " " << 8+i << ":" << mfp.varNameVec[i];
    hss << endl;

    // allocate ss (ss has a private copy constructor, so must use a ptr)...

    assert(mfp.ss == NULL);
    mfp.ss = new std::stringstream();
    *mfp.ss << setprecision(8);

    // make sure the file can exist (i.e. build sub-directories if requested in name)...

    MiscUtils::mkdir_for_file(mfp.name);

    // for the case of single-file-mode, put the header in mfp.ss right now. Otherwise,
    // we store the header to dump at the top of each file...

    if (mfp.format == SINGLE_FILE) {
      *mfp.ss << hss.rdbuf();
    }
    else {
      assert(mfp.format == MULTI_FILE);
      mfp.header = hss.str();
    }

    cout << "============= HEADER HEADER ===============\n" << hss.str() << "=============== HEADER HEADER ==================" << endl;

  }

  return 0; // OK

}

void StaticSolver::processFluxProbe(Param* param,const int step,const double time,const bool killfile_request) {
  // use the full param as the hash for the map..
  if ( fpMap.find(param->str()) == fpMap.end()) {

    // couldnt find this probe so we need to insert it to the list
    pair<map<const string,FluxProbe>::iterator,bool> ret = fpMap.insert(pair<const string,FluxProbe>(param->str(), FluxProbe()));
    assert( ret.second);

    const int ierr = initFluxProbe(ret.first->second,&(*param));

    if ( ierr != 0) {
      // see the note on the error handling logic above for the PointProbe
      if ( !killfile_request ) {
        CERR( " > unable to add FluxProbe from: \n >   " << param->str() );
      } else {
        CWARN( " > unable to add FluxProbe from: \n >  " << param->str() );
        fpMap.erase(ret.first);
      }

    } else {
      // if the interval is invalid, then its one-and-done for this probe
      // otherwise register for continuous processing
      if ( ret.first->second.interval == -1) {
        ret.first->second.doProbe(step,time);
        ret.first->second.flush();
        fpMap.erase(ret.first);
      }
    }
  }
}

int StaticSolver::initFluxProbe(FluxProbe& fp, Param* param) {

  // static counter to round robin the assignment of the rank0..

  static int fp_rank0 = 0;
  fp.probe_rank0      = fp_rank0++;
  if ( fp_rank0 == mpi_size)
    fp_rank0 = 0;

  stringstream hss;

  if ( mpi_rank == fp.probe_rank0) {
    hss << setprecision(8);
    hss << "# " << param->str() << endl;
    hss << "# mles hash id " << RestartHashUtilities::mlesHash << endl;
    hss << "# sles hash id " << RestartHashUtilities::slesHash << endl;
  }

  double xp[3], np[3];
  bool got_xp = false;
  bool got_np = false;
  vector<string> var_vec;

  double dn   = 0.0; // default normal displacement..

  int iarg = 0;
  while ( iarg < param->size()) {

    string token = MiscUtils::toUpperCase(param->getString(iarg++));

    if ( token == "NAME" ) {
      fp.name = param->getString(iarg++);

    } else if (token == "INTERVAL" || token == "WRITE_INTERVAL") {
      fp.interval = param->getInt(iarg++); // default is -1
      if (fp.interval <= 0) {
        CWARN(" > INTERVAL expects a positive integer; treating as one-off request");
        // invalid value entered, so treat as unset (default value)
        fp.interval = -1;
      }

    } else if ( token == "XP") {

      for (int i =0; i < 3; ++i)
        xp[i] = param->getDouble(iarg++);
      got_xp = true;

    } else if ( token == "NP") {

      for (int i =0; i < 3; ++i)
        np[i] = param->getDouble(iarg++);

      const double mag = MAG(np);
      assert( mag > 0.0);

      for (int i = 0; i < 3; ++i)
        np[i] /= mag;
      got_np = true;

    } else if ( token == "DN") {
      dn = param->getDouble(iarg++);


    } else if ( token == "VARS") {

      while (iarg < param->size()) {
        var_vec.push_back(param->getString(iarg++));
      }

    } else {

      if ( mpi_rank == 0)
        cout << " processFluxProbe: skipping unrecognized FLUX_PROBE token " << token << endl;

    }
  }

  int ierr = 0;
  if ( fp.name == "") {

    if ( mpi_rank == 0 )
      cout << "Warning: FLUX_PROBE missing NAME for probe \n " << param->str() << endl;

    ierr = -1;
  }

  if ( !got_xp) {

    if ( mpi_rank == 0 )
      cout << "Warning: FLUX_PROBE missing XP for probe \n " << param->str() << endl;

    ierr = -1;
  }

  if ( !got_np) {

    if ( mpi_rank == 0 )
      cout << "Warning: FLUX_PROBE missing NP for probe \n " << param->str() << endl;

    ierr = -1;
  }

  if ( ierr != 0) {
    if ( mpi_rank == 0 )
      cout << " Warning: skipping this FLUX_PROBE. Correct syntax. " << endl;
    return -1;
  }

  // adjust the xp by the dn in the unit direction

  for (int i =0; i < 3; ++i)
    xp[i] += dn*np[i];

  // var type check ..

  assert( fp.var_vec.empty());
  for (int ii = 0; ii < int(var_vec.size()); ++ii) {

    CtiRegister::CtiData * data = CtiRegister::getCtiData(var_vec[ii]);

    if ( (data== NULL)||(data->getTopology() != SIGNED_FA_DATA) || (data->getType() != DN_DATA)) {
      CWARN(" > flux probe data : " << var_vec[ii] << " does not evaluate to SIGNED_FA_DN.  skipping.");
    } else {
      fp.var_vec.push_back(var_vec[ii]);
    }
  }

  if ( fp.var_vec.empty()) {
    CWARN(" flux probe has no valid data. ");
    return -1;
  }

  int * cv_flag     = new int[ncv_g];
  for (int icv = 0; icv < ncv; ++icv) {
    const double dn =
      (x_cv[icv][0]-xp[0])*np[0] +
      (x_cv[icv][1]-xp[1])*np[1] +
      (x_cv[icv][2]-xp[2])*np[2];
    if (dn >= 0.0) {
      cv_flag[icv] = 1;
    }
    else {
      cv_flag[icv] = -1;
    }
  }

  updateCvData(cv_flag);

  for (int icv = 0; icv < ncv; ++icv) {
    for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      if (cv_flag[icv]*cv_flag[icv_nbr] < 0) {
        // an odd product means this cv is on the boundary...
        cv_flag[icv] *= 2; // makes -1 -> -2
        assert((cv_flag[icv] == 2)||(cv_flag[icv] == -2));
        break;
      }
    }
  }

  updateCvData(cv_flag);

  for (int icv = 0; icv < ncv; ++icv) {
    if ((cv_flag[icv] == 1)||(cv_flag[icv] == -1)) {
      for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        if ((cv_flag[icv_nbr] == 2)||(cv_flag[icv_nbr] == -2)) {
          cv_flag[icv] *= 4;
          assert((cv_flag[icv] == 4)||(cv_flag[icv] == -4));
          break;
        }
      }
    }
  }

  updateCvData(cv_flag);

  // now we have the zone of cvs around flux faces 2-deep all flagged with
  // +/- 2 or +/- 4. Hopefully 2-deep is sufficient to keep cv groups
  // together. Now re-index these cvs...

  int my_ncv_active = 0;
  for (int icv = 0; icv < ncv; ++icv) {
    if ((cv_flag[icv] > 1)||(cv_flag[icv] < -1)) {
      cv_flag[icv] = my_ncv_active++;
    }
    else {
      cv_flag[icv] = -1;
    }
  }

  int * acora = NULL;
  MiscUtils::buildXora(acora,my_ncv_active);
  assert(acora[mpi_size] < TWO_BILLION);

  for (int icv = 0; icv < ncv; ++icv) {
    if (cv_flag[icv] >= 0) {
      cv_flag[icv] += acora[mpi_rank];
    }
  }

  // group contiguous chunks of flagged cells together -- this is similar
  // to the operations done by the multifluxprobe.  eventually we are going
  // to operate on the single group we are interested in.

  if (mpi_rank == 0)
    cout << " > group";

  int done = 0;
  while (done == 0) {

    if (mpi_rank == 0) {
      cout << ".";
      cout.flush();
    }

    int my_done = 1;

    updateCvData(cv_flag);

    for (int icv = 0; icv < ncv; ++icv) {
      if (cv_flag[icv] >= 0) {
        for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          if ((cv_flag[icv_nbr] >= 0)&&(cv_flag[icv_nbr] < cv_flag[icv])) {
            cv_flag[icv] = cv_flag[icv_nbr];
            my_done = 0;
          }
        }
      }
    }

    MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);

  }

  if (mpi_rank == 0)
    cout << "OK" << endl;

  // the multifluxprobe further compressed to get a final set of groups [0,ngr]
  // but we only need to identify the unindexed "color" associated with the
  // point specified by the flux probe ..

  // reporting the total number of groups is actually not necessary -- but helps debug for now

  set<int> groups;
  for (int icv = 0; icv < ncv; ++icv) {
    if ( (cv_flag[icv] >= acora[mpi_rank]) && (cv_flag[icv] < acora[mpi_rank+1]))
      groups.insert(cv_flag[icv]);
  }

  int my_ngr = groups.size();
  int ngr_global;
  MPI_Allreduce(&my_ngr,&ngr_global,1,MPI_INT,MPI_SUM,mpi_comm);

  if ( mpi_rank == 0 )
    cout << " > total number of flagged grps : " << ngr_global << endl;

  if (ngr_global == 0) {
    delete[] cv_flag;
    delete[] acora;
    if (mpi_rank == 0) cout << "Warning: Skipping this FLUX_PROBE.\n"
                            << "No cvs were tagged to participate in the flux probe."
                            << "Check to confirm your probe location is inside your computational domain."
                            << endl;
    return(-1);
  }

  double d2_min = HUGE_VAL;
  int icv_min   = -1;

  for (int icv = 0; icv < ncv; ++icv) {

    if ( cv_flag[icv] >= 0 ) {

      const double this_d2 = DIST2(xp,x_cv[icv]);
      if ( this_d2 < d2_min ) {
        d2_min  = this_d2;
        icv_min = icv;
      }
    }
  }

  DoubleInt my_di, di;
  my_di.this_double = d2_min;
  my_di.this_int    = mpi_rank;
  MPI_Allreduce(&my_di,&di,1,MPI_DOUBLE_INT,MPI_MINLOC, mpi_comm);

  // broadcast the final color back to all of the ranks..

  int icv_color;
  if ( di.this_int == mpi_rank) {

    assert( (icv_min >= 0) && (icv_min < ncv));
    icv_color = cv_flag[icv_min];
    assert ( icv_color >= 0);
    MPI_Bcast(&icv_color,1,MPI_INT,mpi_rank,mpi_comm);

  } else {

    MPI_Bcast(&icv_color,1,MPI_INT,di.this_int,mpi_comm);

  }

  if ( mpi_rank == 0 )
    cout << " > group color associated with the fp point : " << icv_color << endl;

  // flag all the cells that retain this single color

  for (int icv = 0; icv < ncv; ++icv) {

    if ( cv_flag[icv] == icv_color ) {

      // swap the indexing to be -2 on the opposite side
      const double dn =
        (x_cv[icv][0]-xp[0])*np[0] +
        (x_cv[icv][1]-xp[1])*np[1] +
        (x_cv[icv][2]-xp[2])*np[2];
      if (dn < 0.0) {
        cv_flag[icv] = -cv_flag[icv]-2; // use -2 indexing...
      }

    } else if ( cv_flag[icv] >= 0 ) {

      cv_flag[icv] = -1; // not associated with a group we are interested in..

    } else {

      assert( cv_flag[icv] == -1); // should have already been -1

    }
  }

  updateCvData(cv_flag);

  vector<int> face_vec;
  set<int> zone_set;

  for (int icv = 0; icv < ncv; ++icv) {
    if ( cv_flag[icv] != -1) {
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa  = faocv_v[foc];
        const int icv0 = cvofa[ifa][0];
        const int icv1 = cvofa[ifa][1];
        assert( (icv0 == icv)||(icv1 == icv));

        if ( cv_flag[icv0] == -cv_flag[icv1]-2) {

          // this is a flux face--- add it into the face_vec; the
          // sign of the face will denote its orientation

          if ( cv_flag[icv] < 0 ) {

            if ( icv == icv0) {
              face_vec.push_back(ifa); // outward pointing
            } else {
              assert( icv == icv1);
              face_vec.push_back(-ifa-1);
            }
          }

          for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
            const int ibf = bfocv_v[boc];
            zone_set.insert(zone_bf[ibf]);
          }

        }
      }
    }
  }

  delete[] cv_flag;
  delete[] acora;

  // reduce the zone_set info to the probe rank0..

  string fp_zone;

  int my_count = zone_set.size();
  int * count  = NULL;
  if ( mpi_rank == fp.probe_rank0)
    count = new int[mpi_size];

  MPI_Gather(&my_count,1,MPI_INT,count,1,MPI_INT,fp.probe_rank0,mpi_comm);

  if ( mpi_rank == fp.probe_rank0) {

    int ibuf_size = 0;
    for (int rank = 0; rank < mpi_size; ++rank) {
      if ( rank != mpi_rank)
        ibuf_size = max(ibuf_size,count[rank]);
    }

    int * ibuf = new int[ibuf_size];

    for (int rank = 0; rank < mpi_size; ++rank) {

      if ( (rank != mpi_rank)&&(count[rank] > 0)) {

        MPI_Recv(ibuf,count[rank],MPI_INT,rank,42424,mpi_comm,MPI_STATUS_IGNORE);
        for (int i =0; i < count[rank]; ++i) {
          const int izone = ibuf[i];
          assert( (izone >= 0) && ( izone < int(bfZoneVec.size())));
          zone_set.insert(izone);
        }
      }
    }

    delete[] ibuf;

    for ( set<int>::iterator it = zone_set.begin(); it != zone_set.end(); ++it) {

      if ( it == zone_set.begin())
        fp_zone = bfZoneVec[*it].getName();
      else
        fp_zone.append(":"+bfZoneVec[*it].getName());

    }

    zone_set.clear();

  } else if ( my_count > 0 ) {

    int * ibuf = new int[my_count];
    int i = 0;
    for (set<int>::iterator it = zone_set.begin(); it != zone_set.end(); ++it)
      ibuf[i++] = *it;

    assert ( i == my_count);
    MPI_Ssend(ibuf,my_count,MPI_INT,fp.probe_rank0,42424,mpi_comm);

    delete[] ibuf;

    zone_set.clear();
  }


  if ( mpi_rank == fp.probe_rank0)
    delete[] count;

  // since there is only one set of faces to account for (no grouping)
  // we will directly store the face vec here ..

  sort( face_vec.begin(), face_vec.end());
  fp.face_vec.swap(face_vec);
  face_vec.clear();


  // reduce the projected area and centroid of the flux probe

  double pack_buf[4] = {0.0,0.0,0.0,0.0};

  double my_proj_area_min = 0.0;
  for (int ii = 0; ii < int(fp.face_vec.size()); ++ii) {
    int ifa        = fp.face_vec[ii];
    double fa_sign = 1.0;
    if ( ifa < 0) {
      ifa = -ifa -1;
      fa_sign = -1.0;
    }
    assert ((ifa >=0) &&(ifa < nfa));

    const double proj_area = fa_sign * DOT_PRODUCT(n_fa[ifa],np);
    my_proj_area_min       = min(my_proj_area_min,proj_area);

    pack_buf[0]           += proj_area*x_fa[ifa][0];
    pack_buf[1]           += proj_area*x_fa[ifa][1];
    pack_buf[2]           += proj_area*x_fa[ifa][2];
    pack_buf[3]           += proj_area;
  }

  double unpack_buf[4];
  fp.probe_reduce(pack_buf,unpack_buf,4);

  // report ..
  if ( mpi_rank == fp.probe_rank0 ) {

    for (int i =0; i < 3; ++i)
      fp.probe_x[i] = unpack_buf[i] / unpack_buf[3];

    fp.projected_area = unpack_buf[3];

    cout << COUT_VEC(fp.probe_x) << "  " << fp.projected_area << endl;

    hss << "# probe: zones: \" " << fp_zone << "\"" << endl;
    hss << "# 1:name 2:step 3:time 4:xp 5:yp 6:zp 7:proj_area";
    for (int i =0; i < int(fp.var_vec.size()); ++i)
      hss << " " << 8+i << ":" << fp.var_vec[i];
    hss << endl;

    assert( fp.ss == NULL);
    fp.ss = new stringstream();
    *fp.ss << setprecision(8);

    MiscUtils::mkdir_for_file(fp.name);

    *fp.ss << hss.rdbuf();
    cout << " ===== HEADER =====\n" << hss.str() << " ==== HEADER ==== " << endl;
  }

  return 0;
}

void StaticSolver::processVolumetricProbe(Param* param,const int step,const double time,const bool killfile_request) {
  // use the full param as the hash for the map..
  if ( vpMap.find(param->str()) == vpMap.end()) {
    // couldnt find this probe so we need to insert it to the list
    pair<map<const string,VolumetricProbe>::iterator,bool> ret = vpMap.insert(pair<const string,VolumetricProbe>(param->str(), VolumetricProbe()));
    assert( ret.second);

    const int ierr = initVolumetricProbe(ret.first->second,&(*param));

    if ( ierr != 0) {

      // see the note on the error handling logic above for the PointProbe
      if ( !killfile_request ) {
        CERR( " > unable to add VolumetricProbe from: \n >   " << param->str() );
      } else {
        CWARN( " > unable to add VolumetricProbe from: \n >  " << param->str() );
        vpMap.erase(ret.first);
      }

    } else {
      // if the interval is invalid, then its one-and-done for this probe
      // otherwise register for continuous processing
      if ( ret.first->second.interval == -1 ) {
        ret.first->second.doProbe(step,time);
        ret.first->second.flush();
        vpMap.erase(ret.first);
      }
    }
  }
}

int StaticSolver::initVolumetricProbe(VolumetricProbe& vp, Param* param) {

  // static counter to round robin the assignment of the rank0..

  static int vp_rank0 = 0;
  vp.probe_rank0 = vp_rank0++;
  if (vp_rank0 == mpi_size)
    vp_rank0 = 0;

  stringstream hss;

  if ( mpi_rank == vp.probe_rank0) {
    hss << setprecision(8);
    hss << "# " << param->str() << endl;
    hss << "# mles hash id " << RestartHashUtilities::mlesHash << endl;
    hss << "# sles hash id " << RestartHashUtilities::slesHash << endl;
  }

  SimpleGeom* geom = NULL;
  vector<string> var_vec;

  int iarg = 0;
  bool b_surface = false;
  while ( iarg < param->size()) {

    string token = MiscUtils::toUpperCase(param->getString(iarg++));

    if ( token == "NAME" ) {
      vp.name = param->getString(iarg++);

    }
    else if (token == "INTERVAL" || token == "WRITE_INTERVAL") {
      vp.interval = param->getInt(iarg++); // default is -1
      if (vp.interval <= 0) {
        CWARN(" > INTERVAL expects a positive integer; treating as one-off request");
        // invalid value entered, so treat as unset (default value)
        vp.interval = -1;
      }

    }
    else if (token == "GEOM") {

      geom = newSimpleGeom(param,iarg);

    }
    else if ( token == "SURFACE_ONLY" ) {

      b_surface = true;

    }
    else if ( token == "VARS") {

      while (iarg < param->size()) {
        var_vec.push_back(param->getString(iarg++));
      }

    }
    else {

      if ( mpi_rank == 0)
        cout << " processVolumetricProbe: skipping unrecognized VOLUMETRIC_PROBE token " << token << endl;

    }
  }

  int ierr = 0;
  if ( vp.name == "") {

    if ( mpi_rank == 0 )
      cout << "Warning: VOLUMETRIC_PROBE missing NAME for probe \n " << param->str() << endl;

    ierr = -1;
  }

  if ( geom == NULL) {

    if ( mpi_rank == 0 )
      cout << "Warning: VOLUMETRIC_PROBE missing/unsupported GEOM for probe \n " << param->str() << endl;

    ierr = -1;
  }

  if ( ierr != 0) {
    if ( mpi_rank == 0 )
      cout << " Warning: skipping this VOLUMETRIC_PROBE. Correct syntax. " << endl;
    return -1;
  }

  // var type check ..

  assert( vp.var_vec.empty());
  for (int ii = 0; ii < int(var_vec.size()); ++ii) {

    CtiRegister::CtiData * data = CtiRegister::getCtiData(var_vec[ii]);

    if ( (data== NULL)||(data->getUnindexedTopology() != CV_DATA) || (data->getType() != DN_DATA)) {
      CWARN(" > volumetric probe data : " << var_vec[ii] << " does not evaluate to CV_DN.  skipping.");
    }
    else {
      vp.var_vec.push_back(var_vec[ii]);
    }
  }

  if ( vp.var_vec.empty()) {
    CWARN(" volumetric probe has no valid data. ");
    return -1;
  }

  vector<int> cv_vec;
  if (b_surface) {
    FOR_IBF {
      const int icv = cvobf[ibf];
      if (geom->pointIsInside(x_cv[icv]))
        cv_vec.push_back(icv);
    }
  }
  else {
    FOR_ICV {
      if (geom->pointIsInside(x_cv[icv]))
        cv_vec.push_back(icv);
    }
  }
  delete geom; geom = NULL;

  // already sorted by icv...
  vp.cv_vec.swap(cv_vec);
  cv_vec.clear();

  // reduce the volume and centroid of the volumetric probe

  double pack_buf[4] = {0.0,0.0,0.0,0.0};

  for (int ii = 0; ii < int(vp.cv_vec.size()); ++ii) {
    int icv = vp.cv_vec[ii];
    assert ((icv >=0) &&(icv < ncv));
    pack_buf[0] += vol_cv[icv]*x_cv[icv][0];
    pack_buf[1] += vol_cv[icv]*x_cv[icv][1];
    pack_buf[2] += vol_cv[icv]*x_cv[icv][2];
    pack_buf[3] += vol_cv[icv];
  }

  double unpack_buf[4];
  vp.probe_reduce(pack_buf,unpack_buf,4);

  int8 my_count = vp.cv_vec.size();
  int8 count;
  MPI_Reduce(&my_count,&count,1,MPI_INT8,MPI_SUM,vp.probe_rank0,mpi_comm);

  // report ..
  if ( mpi_rank == vp.probe_rank0 ) {

    for (int i =0; i < 3; ++i)
      vp.probe_x[i] = unpack_buf[i] / unpack_buf[3];

    vp.probe_vol = unpack_buf[3];

    cout << COUT_VEC(vp.probe_x) << "  " << vp.probe_vol << endl;

    hss << "# name: " << vp.name << " location: " << vp.probe_x[0] << " " << vp.probe_x[1] << " " << vp.probe_x[2] << " volume: " << vp.probe_vol << " count: " << count << endl;
    hss << "# 1:step 2:time ";
    for (int i =0; i < int(vp.var_vec.size()); ++i)
      hss << " " << 3+i << ":" << vp.var_vec[i];
    hss << endl;

    assert( vp.ss == NULL);
    vp.ss = new stringstream();
    *vp.ss << setprecision(8);

    MiscUtils::mkdir_for_file(vp.name);

    *vp.ss << hss.rdbuf();
    cout << " ===== HEADER =====\n" << hss.str() << " ==== HEADER ==== " << endl;
  }

  return 0;
}

void StaticSolver::processSMDVolumetricProbe(Param* param,const int step,const double time,const bool killfile_request) {

  // use the full param as the hash for the map..
  if ( svpMap.find(param->str()) == svpMap.end()) {
    // couldnt find this probe so we need to insert it to the list
    pair<map<const string,VolumetricProbe>::iterator,bool> ret = svpMap.insert(pair<const string,VolumetricProbe>(param->str(), VolumetricProbe()));
    assert( ret.second);

    const int ierr = initSMDVolumetricProbe(ret.first->second,&(*param));

    if ( ierr != 0) {

      // see the note on the error handling logic above for the PointProbe
      if ( !killfile_request ) {
        CERR( " > unable to add SMDVolumetricProbe from: \n >   " << param->str() );
      } else {
        CWARN( " > unable to add SMDVolumetricProbe from: \n >  " << param->str() );
        svpMap.erase(ret.first);
      }

    } else {
      // if the interval is invalid, then its one-and-done for this probe
      // otherwise register for continuous processing
      if ( ret.first->second.interval == -1 ) {
        ret.first->second.doSMDProbe(step,time);
        ret.first->second.flush();
        svpMap.erase(ret.first);
      }
    }
  }
}

int StaticSolver::initSMDVolumetricProbe(VolumetricProbe& svp, Param* param) {

  // static counter to round robin the assignment of the rank0..

  static int svp_rank0 = 0;
  svp.probe_rank0 = svp_rank0++;
  if (svp_rank0 == mpi_size)
    svp_rank0 = 0;

  stringstream hss;

  if ( mpi_rank == svp.probe_rank0) {
    hss << setprecision(8);
    hss << "# " << param->str() << endl;
    hss << "# mles hash id " << RestartHashUtilities::mlesHash << endl;
    hss << "# sles hash id " << RestartHashUtilities::slesHash << endl;
  }

  SimpleGeom* geom = NULL;
  vector<string> var_vec;

  int iarg = 0;
  bool b_surface = false;
  while ( iarg < param->size()) {

    string token = MiscUtils::toUpperCase(param->getString(iarg++));

    if ( token == "NAME" ) {
      svp.name = param->getString(iarg++);

    }
    else if (token == "INTERVAL" || token == "WRITE_INTERVAL") {
      svp.interval = param->getInt(iarg++); // default is -1
      if (svp.interval <= 0) {
        CWARN(" > INTERVAL expects a positive integer; treating as one-off request");
        // invalid value entered, so treat as unset (default value)
        svp.interval = -1;
      }

    }
    else if (token == "GEOM") {

      geom = newSimpleGeom(param,iarg);

    }
    else if ( token == "SURFACE_ONLY" ) {

      b_surface = true;

    }
    else if ( token == "VOL") {
      var_vec.push_back(param->getString(iarg++));
    }
    else if ( token == "SUR") {
      var_vec.push_back(param->getString(iarg++));
    }
    else {

      if ( mpi_rank == 0)
        cout << " processSMDVolumetricProbe: skipping unrecognized SMD_VOLUMETRIC_PROBE token " << token << endl;

    }
  }

  // enforcing variables to compute SMD

  int ierr = 0;
  if ( svp.name == "") {

    if ( mpi_rank == 0 )
      cout << "Warning: SMD_VOLUMETRIC_PROBE missing NAME for probe \n " << param->str() << endl;

    ierr = -1;
  }

  if ( geom == NULL) {

    if ( mpi_rank == 0 )
      cout << "Warning: SMD_VOLUMETRIC_PROBE missing/unsupported GEOM for probe \n " << param->str() << endl;

    ierr = -1;
  }

  if ( ierr != 0) {
    if ( mpi_rank == 0 )
      cout << " Warning: skipping this SMD_VOLUMETRIC_PROBE. Correct syntax. " << endl;
    return -1;
  }

  // var type check ..

  assert( svp.var_vec.empty());
  for (int ii = 0; ii < int(var_vec.size()); ++ii) {

    CtiRegister::CtiData * data = CtiRegister::getCtiData(var_vec[ii]);

    if ( (data== NULL)||(data->getUnindexedTopology() != CV_DATA) || (data->getType() != DN_DATA)) {
      CWARN(" > volumetric probe data : " << var_vec[ii] << " does not evaluate to CV_DN.  skipping.");
    }
    else {
      svp.var_vec.push_back(var_vec[ii]);
    }
  }

  if ( svp.var_vec.empty()) {
    CWARN(" smd volumetric probe has no valid data. ");
    return -1;
  }

  vector<int> cv_vec;
  if (b_surface) {
    FOR_IBF {
      const int icv = cvobf[ibf];
      if (geom->pointIsInside(x_cv[icv]))
        cv_vec.push_back(icv);
    }
  }
  else {
    FOR_ICV {
      if (geom->pointIsInside(x_cv[icv]))
        cv_vec.push_back(icv);
    }
  }
  delete geom; geom = NULL;

  // already sorted by icv...
  svp.cv_vec.swap(cv_vec);
  cv_vec.clear();

  // reduce the volume and centroid of the volumetric probe

  double pack_buf[4] = {0.0,0.0,0.0,0.0};

  for (int ii = 0; ii < int(svp.cv_vec.size()); ++ii) {
    int icv = svp.cv_vec[ii];
    assert ((icv >=0) &&(icv < ncv));
    pack_buf[0] += vol_cv[icv]*x_cv[icv][0];
    pack_buf[1] += vol_cv[icv]*x_cv[icv][1];
    pack_buf[2] += vol_cv[icv]*x_cv[icv][2];
    pack_buf[3] += vol_cv[icv];
  }

  double unpack_buf[4];
  svp.probe_reduce(pack_buf,unpack_buf,4);

  int8 my_count = svp.cv_vec.size();
  int8 count;
  MPI_Reduce(&my_count,&count,1,MPI_INT8,MPI_SUM,svp.probe_rank0,mpi_comm);

  // report ..
  if ( mpi_rank == svp.probe_rank0 ) {

    for (int i =0; i < 3; ++i)
      svp.probe_x[i] = unpack_buf[i] / unpack_buf[3];

   svp.probe_vol = unpack_buf[3];

    cout << COUT_VEC(svp.probe_x) << "  " << svp.probe_vol << endl;

    hss << "# name: " << svp.name << " location: " << svp.probe_x[0] << " " << svp.probe_x[1] << " " << svp.probe_x[2] << " volume: " << svp.probe_vol << " count: " << count << endl;
    hss << "# 1:step 2:time ";
    for (int i =0; i < int(svp.var_vec.size()); ++i)
      hss << " " << 3+i << ":" << svp.var_vec[i];
    hss <<  " " << 3+int(svp.var_vec.size()) << ":" << "SMD" << endl;

    assert( svp.ss == NULL);
    svp.ss = new stringstream();
    *svp.ss << setprecision(8);

    MiscUtils::mkdir_for_file(svp.name);

    *svp.ss << hss.rdbuf();
    cout << " ===== HEADER =====\n" << hss.str() << " ==== HEADER ==== " << endl;
  }

  return 0;
}

void StaticSolver::processFazoneProbe(Param* param,const int step,const double time,const bool killfile_request) {
  // use the full param as the hash for the map..
  if ( fzpMap.find(param->str()) == fzpMap.end()) {
    // couldnt find this probe so we need to insert it to the list
    pair<map<const string,FazoneProbe>::iterator,bool> ret = fzpMap.insert(pair<const string,FazoneProbe>(param->str(), FazoneProbe()));
    assert( ret.second);

    const int ierr = initFazoneProbe(ret.first->second,&(*param));

    if ( ierr != 0) {

      // see the note on the error handling logic above for the PointProbe
      if ( !killfile_request ) {
        CERR( " > unable to add FazoneProbe from: \n >   " << param->str() );
      } else {
        CWARN( " > unable to add FazoneProbe from: \n >  " << param->str() );
        fzpMap.erase(ret.first);
      }

    } else {
      // if the interval is invalid, then its one-and-done for this probe
      // otherwise register for continuous processing
      if ( ret.first->second.interval == -1 ) {
        ret.first->second.doProbe(step,time);
        ret.first->second.flush();
        fzpMap.erase(ret.first);
      }
    }
  }
}

int StaticSolver::initFazoneProbe(FazoneProbe& fzp, Param* param) {

  // static counter to round robin the assignment of the rank0..

  static int fzp_rank0 = 0;
  fzp.probe_rank0 = fzp_rank0++;
  if (fzp_rank0 == mpi_size)
    fzp_rank0 = 0;

  stringstream hss;

  if ( mpi_rank == fzp.probe_rank0) {
    hss << setprecision(8);
    hss << "# " << param->str() << endl;
    hss << "# mles hash id " << RestartHashUtilities::mlesHash << endl;
    hss << "# sles hash id " << RestartHashUtilities::slesHash << endl;
  }

  SimpleGeom* geom = NULL;
  vector<string> var_vec;

  int iarg = 0;
  while ( iarg < param->size()) {

    string token = MiscUtils::toUpperCase(param->getString(iarg++));

    if ( token == "NAME" ) {
      fzp.name = param->getString(iarg++);

    }
    else if ( (token == "ZONE") || (token == "FAZONE") || (token == "BFZONE") ) {
      fzp.setZone(param->getString(iarg++));

    }
    else if (token == "INTERVAL" || token == "WRITE_INTERVAL") {
      fzp.interval = param->getInt(iarg++); // default is -1
      if (fzp.interval <= 0) {
        CWARN(" > INTERVAL expects a positive integer; treating as one-off request");
        // invalid value entered, so treat as unset (default value)
        fzp.interval = -1;
      }

    }
    else if (token == "GEOM") {

      geom = newSimpleGeom(param,iarg);

    }
    else if ( token == "VARS") {

      while (iarg < param->size()) {
        var_vec.push_back(param->getString(iarg++));
      }

    }
    else {

      if ( mpi_rank == 0)
        cout << " processFazoneProbe: skipping unrecognized FAZONE_PROBE token " << token << endl;

    }
  }

  int ierr = 0;
  if ( fzp.name == "") {

    if ( mpi_rank == 0 )
      cout << "Warning: FAZONE_PROBE missing NAME for probe \n " << param->str() << endl;

    ierr = -1;
  }

  BfZone * bf_zone = NULL;
  bf_zone = getBfZone(fzp.zone);

  if ((fzp.zone == "") || (bf_zone == NULL)) {

    if ( mpi_rank == 0 )
      cout << "Warning: FAZONE_PROBE missing valid ZONE for probe \n " << param->str() << endl;

    ierr = -1;
  }

  // if ( geom == NULL) {
  //
  //   if ( mpi_rank == 0 )
  //     cout << "Warning: VOLUMETRIC_PROBE missing/unsupported GEOM for probe \n " << param->str() << endl;
  //
  //   ierr = -1;
  // }

  if ( ierr != 0) {
    if ( mpi_rank == 0 )
      cout << " Warning: skipping this FAZONE_PROBE. Correct syntax. " << endl;
    return -1;
  }

  // var type check ..

  assert( fzp.var_vec.empty());
  for (int ii = 0, iend=var_vec.size(); ii < iend; ++ii) {
    // update var name based on zone:variable combination
    string full_var_name = var_vec[ii];
    if (full_var_name.find(":") == string::npos) full_var_name = fzp.zone + ":" + var_vec[ii];
    CtiRegister::CtiData * data = CtiRegister::getCtiData(full_var_name);

    if ( (data == NULL)||(data->getUnindexedTopology() != BF_DATA) || (data->getType() != DN_DATA)) {
      CWARN(" > fazone probe variable : " << var_vec[ii] << " does not evaluate to FA_DN.  skipping.");
    }
    else {
      fzp.var_vec.push_back(full_var_name);
    }
  }

  if ( fzp.var_vec.empty()) {
    CWARN(" fazone probe has no valid data. ");
    return -1;
  }

  vector<int> bf_vec;
  if (geom) {
    for (int ibf = 0; ibf < bf_zone->nbf; ++ibf) {
      if (geom->pointIsInside(bf_zone->x_bf[ibf]))
        bf_vec.push_back(ibf);
    }
  }
  else {
    for (int ibf = 0; ibf < bf_zone->nbf; ++ibf)
      bf_vec.push_back(ibf);
  }
  delete geom; geom = NULL;

  // already sorted by icv...
  fzp.bf_vec.swap(bf_vec);
  bf_vec.clear();

  // reduce the area and centroid of the fazone probe

  double pack_buf[4] = {0.0,0.0,0.0,0.0};

  for (int ii = 0, iend=fzp.bf_vec.size(); ii < iend; ++ii) {
    int ibf = fzp.bf_vec[ii];
    assert ((ibf >=0) && (ibf < bf_zone->nbf));
    pack_buf[0] += bf_zone->area_bf[ibf]*bf_zone->x_bf[ibf][0];
    pack_buf[1] += bf_zone->area_bf[ibf]*bf_zone->x_bf[ibf][1];
    pack_buf[2] += bf_zone->area_bf[ibf]*bf_zone->x_bf[ibf][2];
    pack_buf[3] += bf_zone->area_bf[ibf];
  }

  double unpack_buf[4];
  fzp.probe_reduce(pack_buf,unpack_buf,4);

  int8 my_count = fzp.bf_vec.size();
  int8 count;
  MPI_Reduce(&my_count,&count,1,MPI_INT8,MPI_SUM,fzp.probe_rank0,mpi_comm);

  // report ..
  if ( mpi_rank == fzp.probe_rank0 ) {

    for (int i =0; i < 3; ++i)
      fzp.probe_x[i] = unpack_buf[i] / unpack_buf[3];

    fzp.probe_area = unpack_buf[3];

    cout << COUT_VEC(fzp.probe_x) << "  " << fzp.probe_area << endl;

    hss << "# name: " << fzp.name << ", zone: " << fzp.zone << ", centroid: " << fzp.probe_x[0] << " " << fzp.probe_x[1] << " " << fzp.probe_x[2] << ", area: " << fzp.probe_area << ", count: " << count << endl;
    hss << "# 1:step 2:time ";
    for (int ii=0, iend=fzp.var_vec.size(); ii < iend; ++ii)
      hss << " " << 3+(2*ii) << ":int(" << fzp.var_vec[ii] << ")dA " <<  4+(2*ii) << ":avg(" << fzp.var_vec[ii] << ")";
    hss << endl;

    assert( fzp.ss == NULL);
    fzp.ss = new stringstream();
    *fzp.ss << setprecision(8);

    MiscUtils::mkdir_for_file(fzp.name);

    *fzp.ss << hss.rdbuf();
    cout << " ===== HEADER =====\n" << hss.str() << " ==== HEADER ==== " << endl;
  }

  return 0;
}

void StaticSolver::processPdfProbe(Param* param,const int step,const double time,const bool killfile_request) {
  // use the full param as the hash for the map..

  if (pdfpMap.find(param->str()) == pdfpMap.end()) {

    // couldn't find this probe so we need to insert it to the list
    pair<map<const string,PdfProbe>::iterator,bool> ret = pdfpMap.insert(pair<const string,PdfProbe>(param->str(),PdfProbe()));
    assert( ret.second); // guaraneed by the find above

    const int ierr = initPdfProbe(ret.first->second,&(*param),killfile_request);
    if (ierr != 0) {
      if ( !killfile_request ) {
        CERR( " > unable to add PdfProbe from: \n >  " << param->str() );
      } else {
        CWARN( " > unable to add PdfProbe from: \n >  " << param->str() );
        pdfpMap.erase(ret.first);
      }
    }
    else {
      if ( killfile_request && (ret.first->second.write_interval == -1) && (ret.first->second.sample_interval == -1) ) {
        // killfile request with no sample or write interval; implies this should
        // be a one-off, so OK to remove it from map
        ret.first->second.write_interval = ret.first->second.sample_interval = 1;
        ret.first->second.doProbe(step,time);
        // return to no sample or write interval to tell flush its a one and done
        ret.first->second.write_interval = ret.first->second.sample_interval = -1;
        ret.first->second.flush();
        pdfpMap.erase(ret.first);
      }
    }

  }
}

int StaticSolver::initPdfProbe(PdfProbe& pdfp,Param* param,const bool killfile) {

  // if killfile == true, then this param is being parsed from a killfile. If
  // false, then this param is one of the input file parameters. This affects the behavior
  // of certain options...

  if (mpi_rank == 0)
    cout << "initPdfProbe() \"" << param->str() << "\"..." << endl;

  // parse args...

  string var_name = "";
  string wvar_name = ""; // for weighted pdf probe (1,vol_cv,rho_cv,etc...)

  int iarg = 0;
  while ( iarg < param->size()) {

    string token = MiscUtils::toUpperCase(param->getString(iarg++));

    if ( token == "NAME" ) {
      pdfp.name = param->getString(iarg++);
    }
    else if ((token == "INTERVAL")||(token == "WRITE_INTERVAL")) {
      pdfp.write_interval = param->getInt(iarg++);
      if (pdfp.write_interval <= 0) {
        CWARN(" > INTERVAL expects a positive integer; treating as one-off request");
        // invalid value entered, so treat as unset (default value)
        pdfp.write_interval = -1;
      }
    }
    else if (token == "SAMPLE_INTERVAL" ) {
      pdfp.sample_interval = param->getInt(iarg++);
      if (pdfp.sample_interval <= 0) {
        CWARN(" > SAMPLE_INTERVAL expects a positive integer; treating as one-off request");
        // invalid value entered, so treat as unset (default value)
        pdfp.sample_interval = -1;
      }
    }
    else if ( token == "VAR") {
      var_name = param->getString(iarg++);
    }
    else if ((token == "NBIN")||(token == "NBINS")||(token == "SIZE")||(token == "N")) {
      pdfp.nbin = param->getInt(iarg++);
    }
    else if (token == "RANGE") {
      // user can fix the range. If not done, the range will be set from the first set of conditional data...
      pdfp.range[0] = param->getDouble(iarg++);
      pdfp.range[1] = param->getDouble(iarg++);
      pdfp.b_range = true;
    }
    else if ( token == "WVAR") {
      wvar_name = param->getString(iarg++);
    }
    else if (token == "RESET_AFTER_WRITE") {
      pdfp.reset_after_write = true;
    }
    else if (token == "GEOM") {
      assert(pdfp.geom == NULL);
      pdfp.geom = newSimpleGeom(param,iarg);
    }
    else {
      if ( mpi_rank == 0)
        cout << "Warning: PDF_PROBE skipping unrecognized token " << token << endl;
    }
  }

  int ierr = 0;

  // NAME...

  if ( pdfp.name == "") {

    if ( mpi_rank == 0 )
      cout << "Warning: PDF_PROBE missing NAME." << endl;

    ierr = -1;
  }

  // whether this is being registered from killfile or not
  const bool one_off = (killfile && (pdfp.write_interval == -1) && (pdfp.sample_interval == -1)) ? true:false;

  // INTERVAL/WRITE_INTERVAL...
  if (!killfile || !one_off) {
    // check values are valid if from input OR
    // being registered from a killfile
    if (pdfp.write_interval <= 0) {

      if ( mpi_rank == 0 )
        cout << "Warning: PDF_PROBE WRITE_INTERVAL wasn't properly specified" << endl;

      ierr = -1;
    }

    if (pdfp.sample_interval <= 0) {

      if ( mpi_rank == 0 )
        cout << "Warning: PDF_PROBE SAMPLE_INTERVAL wasn't properly specified" << endl;

      ierr = -1;
    }

    if (pdfp.write_interval%pdfp.sample_interval != 0) {

      if ( mpi_rank == 0 )
        cout << "Warning: PDF_PROBE WRITE_INTERVAL not divisible by SAMPLE_INTERVAL" << endl;

      ierr = -1;
    }
  }

  // VAR...

  if ( var_name == "") {

    if ( mpi_rank == 0 )
      cout << "Warning: PDF_PROBE missing VAR" << endl;

    ierr = -1;
  }

  // VAR type check...

  CtiRegister::CtiData * var_data = CtiRegister::getCtiData(var_name);
  if ((var_data == NULL)||(var_data->getType() != DN_DATA)) {
    if ( mpi_rank == 0 )
      cout << "Warning: PDF_PROBE var " << var_name << " must evaluate to DN_DATA." << endl;
    ierr = -1;
  }
  else {
    // if this is a killfile version of this init, leave the pdfp.var.second connected to the CtiData, even though it is
    // not registered data. It is just going to be "one-and-done" anyways. For the case of registered data,
    // the pointer is immutable so this saves time...
    if (one_off||(CtiRegister::isRegisteredData(var_data))) {
      pdfp.var.first = var_name;
      pdfp.var.second = var_data;
    }
    else {
      pdfp.var.first = var_name;
      pdfp.var.second = NULL; // use NULL to indicate this in not registered and needs eval...
    }

  }

  // WVAR type check...

  if (wvar_name != "") {
    CtiRegister::CtiData * wvar_data = CtiRegister::getCtiData(wvar_name);
    if ((wvar_data == NULL)||(wvar_data->getType() != DN_DATA)||(wvar_data->getUnindexedTopology() != var_data->getUnindexedTopology())) {
      if ( mpi_rank == 0 )
        cout << "Warning: PDF_PROBE WVAR " << wvar_name << " must evaluate to same topology/type as VAR." << endl;
      ierr = -1;
    }
    else {
      if (one_off||(CtiRegister::isRegisteredData(wvar_data))) {
        pdfp.wvar.first = wvar_name;
        pdfp.wvar.second = wvar_data;
      }
      else {
        pdfp.wvar.first = wvar_name;
        pdfp.wvar.second = NULL; // use NULL to indicate this in not registered and needs eval...
      }
    }
  }

  if (ierr != 0)
    return ierr;

  // complete initialization...

  MiscUtils::mkdir_for_file_collective(pdfp.name);

  return 0;

}

void StaticSolver::processConditionalProbe(Param* param,const int step,const double time,const bool killfile_request) {
  // use the full param as the hash for the map..

  if (cpMap.find(param->str()) == cpMap.end()) {

    // couldn't find this probe so we need to insert it to the list
    pair<map<const string,ConditionalProbe>::iterator,bool> ret = cpMap.insert(pair<const string,ConditionalProbe>(param->str(),ConditionalProbe()));
    assert( ret.second); // guaraneed by the find above

    const int ierr = initConditionalProbe(ret.first->second,&(*param),killfile_request);

    if (ierr != 0) {
      if ( !killfile_request ) {
        CERR( " > unable to add ConditionalProbe from: \n >  " << param->str() );
      } else {
        CWARN( " > unable to add ConditionalProbe from: \n >  " << param->str() );
        cpMap.erase(ret.first);
      }
    }
    else {
      if ( killfile_request && (ret.first->second.write_interval == -1) && (ret.first->second.sample_interval == -1) ) {
        // killfile request with no sample or write interval; implies this should
        // be a one-off, so OK to remove it from map
        ret.first->second.write_interval = ret.first->second.sample_interval = 1;
        ret.first->second.doProbe(step,time);
        ret.first->second.write_interval = ret.first->second.sample_interval = -1;
        ret.first->second.flush();
        cpMap.erase(ret.first);
      }
    }

  }
}

int StaticSolver::initConditionalProbe(ConditionalProbe& cp,Param* param,const bool killfile) {

  // if killfile == true, then this param is being parsed from a killfile. If
  // false, then this param is one of the input file parameters. This affects the behavior
  // of certain options...

  if (mpi_rank == 0)
    cout << "initConditionalProbe() \"" << param->str() << "\"..." << endl;

  // cp header: everyone gets this, because everyone can be a writer...

  /*
  assert(cp.hss == NULL);
  cp.hss = new stringstream();
  *cp.hss << setprecision(8);
  *cp.hss << "# " << param->str() << endl;
  *cp.hss << "# mles hash id " << RestartHashUtilities::mlesHash << endl;
  *cp.hss << "# sles hash id " << RestartHashUtilities::slesHash << endl;
  */

  // parse args...

  string cvar_name;
  bool got_cvar = false;
  vector<string> var_vec;
  string cvar_2_name = ""; // for 2d binning
  string wvar_name = ""; // for weighted conditional probe (1,vol_cv,rho_cv,etc...)

  int iarg = 0;
  while ( iarg < param->size()) {

    string token = MiscUtils::toUpperCase(param->getString(iarg++));

    if ( token == "NAME" ) {
      cp.name = param->getString(iarg++);
    }
    else if ((token == "INTERVAL")||(token == "WRITE_INTERVAL")) {
            cp.write_interval = param->getInt(iarg++);
            if (cp.write_interval <= 0) {
              CWARN(" > INTERVAL expects a positive integer; treating as one-off request");
              // invalid value entered, so treat as unset (default value)
              cp.write_interval = -1;
            }
    }
    else if (token == "SAMPLE_INTERVAL" ) {
            cp.sample_interval = param->getInt(iarg++);
            if (cp.sample_interval <= 0) {
              CWARN(" > SAMPLE_INTERVAL expects a positive integer; treating as one-off request");
              // invalid value entered, so treat as unset (default value)
              cp.sample_interval = -1;
            }
    }
    else if ( token == "CVAR") {
            cvar_name = param->getString(iarg++);
            got_cvar = true;
    }
    else if ((token == "NBIN")||(token == "NBINS")||(token == "SIZE")||(token == "N")) {
            cp.nbin = param->getInt(iarg++);
    }
    else if (token == "RANGE") {
            // user can fix the range. If not done, the range will be set from the first set of conditional data...
            cp.range[0] = param->getDouble(iarg++);
            cp.range[1] = param->getDouble(iarg++);
            cp.b_range = true;
    }
    else if ( token == "CVAR_2") {
            cvar_2_name = param->getString(iarg++);
    }
    else if ((token == "NBIN_2")||(token == "NBINS_2")||(token == "SIZE_2")||(token == "N_2")) {
            cp.nbin_2 = param->getInt(iarg++);
    }
    else if (token == "RANGE_2") {
            // user can fix the range. If not done, the range will be set from the first set of conditional data...
            cp.range_2[0] = param->getDouble(iarg++);
            cp.range_2[1] = param->getDouble(iarg++);
            cp.b_range_2 = true;
    }
    else if (token == "RESET_AFTER_WRITE") {
            cp.reset_after_write = true;
    }
    else if (token == "GEOM") {
            assert(cp.geom == NULL);
            cp.geom = newSimpleGeom(param,iarg);
    }
    else if ( token == "WVAR") {
      wvar_name = param->getString(iarg++);
    }
    else if ((token == "VAR") || ( token == "VARS")) {
            // the rest are vars...
      while (iarg < param->size()) {
        var_vec.push_back(param->getString(iarg++));
      }
    }
    else {
      if ( mpi_rank == 0)
        cout << "Warning: CONDITIONAL_PROBE skipping unrecognized token " << token << endl;
    }
  }

  int ierr = 0;

  // NAME...
  if ( cp.name == "") {
    if ( mpi_rank == 0 )
      cout << "Warning: CONDITIONAL_PROBE missing NAME." << endl;
    ierr = -1;
  }

  // whether this is being registered from killfile or not
  const bool one_off = (killfile && (cp.write_interval == -1) && (cp.sample_interval == -1)) ? true:false;

  // INTERVAL/WRITE_INTERVAL...
  if (!killfile || !one_off) {
    // check values are valid if from input OR
    // if from a killfile with either WRITE or SAMPLE interval specified
    if (cp.write_interval <= 0) {

      if ( mpi_rank == 0 )
        cout << "Warning: CONDITIONAL_PROBE WRITE_INTERVAL wasn't properly specified" << endl;

      ierr = -1;
    }

    if (cp.sample_interval <= 0) {

      if ( mpi_rank == 0 )
        cout << "Warning: CONDITIONAL_PROBE SAMPLE_INTERVAL wasn't properly specified" << endl;

      ierr = -1;
    }

    if (cp.write_interval%cp.sample_interval != 0) {

      if ( mpi_rank == 0 )
        cout << "Warning: CONDITIONAL_PROBE WRITE_INTERVAL not divisible by SAMPLE_INTERVAL" << endl;

      if ((cp.sample_interval > 0)&&(cp.write_interval > 0)) {
        cp.write_interval = max(cp.sample_interval,cp.write_interval - cp.write_interval%cp.sample_interval);
        if ( mpi_rank == 0 )
          cout << " > resetting WRITE_INTERVAL to: " << cp.write_interval << endl;
      }
      else {
        ierr = -1;
      }
    }
  }

  // CVAR...
  if ( !got_cvar) {
    if ( mpi_rank == 0 )
      cout << "Warning: CONDITIONAL_PROBE missing CVAR" << endl;
    ierr = -1;
  }

  // CVAR type check...
  CtiRegister::CtiData * cvar_data = CtiRegister::getCtiData(cvar_name);
  if ((cvar_data == NULL) || (cvar_data->getType() != DN_DATA)) {
    if ( mpi_rank == 0 )
      cout << "Warning: CONDITIONAL_PROBE CVAR " << cvar_name << " must evaluate to DN_DATA." << endl;
    ierr = -1;
  }
  else {
    // if this is a killfile version of this init, leave the cp.cvar.second connected to the CtiData, even though it is
    // not registered data. It is just going to be "one-and-done" anyways. For the case of registered data,
    // the pointer is immutable so this saves time...
    if (one_off||(CtiRegister::isRegisteredData(cvar_data))) {
            cp.cvar.first = cvar_name;
            cp.cvar.second = cvar_data;
    }
    else {
            cp.cvar.first = cvar_name;
            cp.cvar.second = NULL; // use NULL to indicate this in not registered and needs eval...
    }

    // CVAR_2 type check...
    if (cvar_2_name != "") {
      CtiRegister::CtiData * cvar_2_data = CtiRegister::getCtiData(cvar_2_name);
      if ((cvar_2_data == NULL)||(cvar_2_data->getType() != DN_DATA)||(cvar_2_data->getUnindexedTopology() != cvar_data->getUnindexedTopology())) {
        if ( mpi_rank == 0 )
          cout << "Warning: CONDITIONAL_PROBE CVAR_2 " << cvar_2_name << " must evaluate to same topology/type as CVAR." << endl;
        ierr = -1;
      }
      else {
        if (one_off||(CtiRegister::isRegisteredData(cvar_2_data))) {
          cp.cvar_2.first = cvar_2_name;
          cp.cvar_2.second = cvar_2_data;
        }
        else {
          cp.cvar_2.first = cvar_2_name;
          cp.cvar_2.second = NULL; // use NULL to indicate this in not registered and needs eval...
        }
      }
    }

    // WVAR type check...
    if (wvar_name != "") {
      CtiRegister::CtiData * wvar_data = CtiRegister::getCtiData(wvar_name);
      if ((wvar_data == NULL)||(wvar_data->getType() != DN_DATA)||(wvar_data->getUnindexedTopology() != cvar_data->getUnindexedTopology())) {
        if ( mpi_rank == 0 )
          cout << "Warning: CONDITIONAL_PROBE WVAR " << wvar_name << " must evaluate to same topology/type as CVAR." << endl;
        ierr = -1;
      }
      else {
        if (one_off||(CtiRegister::isRegisteredData(wvar_data))) {
          cp.wvar.first = wvar_name;
          cp.wvar.second = wvar_data;
        }
        else {
          cp.wvar.first = wvar_name;
          cp.wvar.second = NULL; // use NULL to indicate this in not registered and needs eval...
        }
      }
    }

    if (ierr != -1) {
      // VARS type check...
      assert( cp.var_vec.empty());
      for (int ii = 0; ii < int(var_vec.size()); ++ii) {
        CtiRegister::CtiData * data = CtiRegister::getCtiData(var_vec[ii]);

        if ((data == NULL)||(data->getUnindexedTopology() != cvar_data->getUnindexedTopology())||(data->getType() != DN_DATA)) {
          if ( mpi_rank == 0 )
            cout << "Warning: CONDITIONAL_PROBE VAR \"" << var_vec[ii] << "\" does not evaluate to same topology/type as CVAR. Skipping." << endl;
        }

        else {
          // see note above on killfile...
          if (one_off||(CtiRegister::isRegisteredData(data))) {
            cp.var_vec.push_back(pair<string,CtiRegister::CtiData*>(var_vec[ii],data));
          }
          else {
            cp.var_vec.push_back(pair<string,CtiRegister::CtiData*>(var_vec[ii],NULL));
          }
        }
      }

      if ( cp.var_vec.empty()) {
        if ( mpi_rank == 0 )
          cout << "Warning: CONDITIONAL_PROBE has no valid data." << endl;
        ierr = -1;
      }
    }

  }

  if (ierr != 0)
    return ierr;

  // complete initialization...

  /*
  *cp.hss << "# 1:bin 2:count 3:" << cvar_name << "_avg 4:" << cvar_name << "_rms";
  for (int i = 0; i < int(cp.var_vec.size()); ++i)
    *cp.hss << " " << 5+i*4 << ":" << cp.var_vec[i].first << "_avg " <<
      6+i*4 << ":" << cp.var_vec[i].first << "_rms " <<
      7+i*4 << ":" << cp.var_vec[i].first << "_min " <<
      8+i*4 << ":" << cp.var_vec[i].first << "_max";
  *cp.hss << endl;
  */

  MiscUtils::mkdir_for_file_collective(cp.name);

  return 0;

}

void StaticSolver::processSMDConditionalProbe(Param* param,const int step,const double time,const bool killfile_request) {
  // use the full param as the hash for the map..

  if (scpMap.find(param->str()) == scpMap.end()) {

    // couldn't find this probe so we need to insert it to the list
    pair<map<const string,ConditionalProbe>::iterator,bool> ret = scpMap.insert(pair<const string,ConditionalProbe>(param->str(),ConditionalProbe()));
    assert( ret.second); // guaraneed by the find above

    const int ierr = initSMDConditionalProbe(ret.first->second,&(*param),killfile_request);

    if (ierr != 0) {
      if ( !killfile_request ) {
        CERR( " > unable to add SMDConditionalProbe from: \n >  " << param->str() );
      } else {
        CWARN( " > unable to add SMDConditionalProbe from: \n >  " << param->str() );
        scpMap.erase(ret.first);
      }
    }
    else {
      if ( killfile_request && (ret.first->second.write_interval == -1) && (ret.first->second.sample_interval == -1) ) {
        // killfile request with no sample or write interval; implies this should
        // be a one-off, so OK to remove it from map
        ret.first->second.write_interval = ret.first->second.sample_interval = 1;
        ret.first->second.doProbe(step,time);
        ret.first->second.write_interval = ret.first->second.sample_interval = -1;
        ret.first->second.flush();
        cpMap.erase(ret.first);
      }
    }

  }
}

int StaticSolver::initSMDConditionalProbe(ConditionalProbe& scp,Param* param,const bool killfile) {

  // if killfile == true, then this param is being parsed from a killfile. If
  // false, then this param is one of the input file parameters. This affects the behavior
  // of certain options...

  if (mpi_rank == 0)
    cout << "initSMDConditionalProbe() \"" << param->str() << "\"..." << endl;

  string cvar_name;
  bool got_cvar = false;
  vector<string> var_vec;
  string cvar_2_name = ""; // for 2d binning
  string wvar_name = ""; // for weighted conditional probe (1,vol_cv,rho_cv,etc...)

  int iarg = 0;
  while ( iarg < param->size()) {

    string token = MiscUtils::toUpperCase(param->getString(iarg++));

    if ( token == "NAME" ) {
      scp.name = param->getString(iarg++);
    }
    else if ((token == "INTERVAL")||(token == "WRITE_INTERVAL")) {
            scp.write_interval = param->getInt(iarg++);
            if (scp.write_interval <= 0) {
              CWARN(" > INTERVAL expects a positive integer; treating as one-off request");
              // invalid value entered, so treat as unset (default value)
              scp.write_interval = -1;
            }
    }
    else if (token == "SAMPLE_INTERVAL" ) {
            scp.sample_interval = param->getInt(iarg++);
            if (scp.sample_interval <= 0) {
              CWARN(" > SAMPLE_INTERVAL expects a positive integer; treating as one-off request");
              // invalid value entered, so treat as unset (default value)
              scp.sample_interval = -1;
            }
    }
    else if ( token == "CVAR") {
            cvar_name = param->getString(iarg++);
            got_cvar = true;
    }
    else if ((token == "NBIN")||(token == "NBINS")||(token == "SIZE")||(token == "N")) {
            scp.nbin = param->getInt(iarg++);
    }
    else if (token == "RANGE") {
            // user can fix the range. If not done, the range will be set from the first set of conditional data...
            scp.range[0] = param->getDouble(iarg++);
            scp.range[1] = param->getDouble(iarg++);
            scp.b_range = true;
    }
    else if ( token == "CVAR_2") {
            cvar_2_name = param->getString(iarg++);
    }
    else if ((token == "NBIN_2")||(token == "NBINS_2")||(token == "SIZE_2")||(token == "N_2")) {
            scp.nbin_2 = param->getInt(iarg++);
    }
    else if (token == "RANGE_2") {
            // user can fix the range. If not done, the range will be set from the first set of conditional data...
            scp.range_2[0] = param->getDouble(iarg++);
            scp.range_2[1] = param->getDouble(iarg++);
            scp.b_range_2 = true;
    }
    else if (token == "RESET_AFTER_WRITE") {
            scp.reset_after_write = true;
    }
    else if (token == "GEOM") {
            assert(scp.geom == NULL);
            scp.geom = newSimpleGeom(param,iarg);
    }
    else if ( token == "WVAR") {
      wvar_name = param->getString(iarg++);
    }
    else if ( token == "VOL") {
      var_vec.push_back(param->getString(iarg++));
    }
    else if ( token == "SUR") { 
      var_vec.push_back(param->getString(iarg++));
    }
    else {
      if ( mpi_rank == 0)
        cout << "Warning: SMD CONDITIONAL_PROBE skipping unrecognized token " << token << endl;
    }
  }

  int ierr = 0;
  if (wvar_name != "vol_cv") 
    wvar_name = "vol_cv";

  // NAME...
  if ( scp.name == "") {
    if ( mpi_rank == 0 )
      cout << "Warning: SMD_CONDITIONAL_PROBE missing NAME." << endl;
    ierr = -1;
  }

  // whether this is being registered from killfile or not
  const bool one_off = (killfile && (scp.write_interval == -1) && (scp.sample_interval == -1)) ? true:false;

  // INTERVAL/WRITE_INTERVAL...
  if (!killfile || !one_off) {
    // check values are valid if from input OR
    // if from a killfile with either WRITE or SAMPLE interval specified
    if (scp.write_interval <= 0) {

      if ( mpi_rank == 0 )
        cout << "Warning: SMD_CONDITIONAL_PROBE WRITE_INTERVAL wasn't properly specified" << endl;

      ierr = -1;
    }

    if (scp.sample_interval <= 0) {

      if ( mpi_rank == 0 )
        cout << "Warning: SMD_CONDITIONAL_PROBE SAMPLE_INTERVAL wasn't properly specified" << endl;

      ierr = -1;
    }

    if (scp.write_interval%scp.sample_interval != 0) {

      if ( mpi_rank == 0 )
        cout << "Warning: SMD_CONDITIONAL_PROBE WRITE_INTERVAL not divisible by SAMPLE_INTERVAL" << endl;

      if ((scp.sample_interval > 0)&&(scp.write_interval > 0)) {
        scp.write_interval = max(scp.sample_interval,scp.write_interval - scp.write_interval%scp.sample_interval);
        if ( mpi_rank == 0 )
          cout << " > resetting WRITE_INTERVAL to: " << scp.write_interval << endl;
      }
      else {
        ierr = -1;
      }
    }
  }

  // CVAR...
  if ( !got_cvar) {
    if ( mpi_rank == 0 )
      cout << "Warning: SMD_CONDITIONAL_PROBE missing CVAR" << endl;
    ierr = -1;
  }

  // CVAR type check...
  CtiRegister::CtiData * cvar_data = CtiRegister::getCtiData(cvar_name);
  if ((cvar_data == NULL) || (cvar_data->getType() != DN_DATA)) {
    if ( mpi_rank == 0 )
      cout << "Warning: SMD_CONDITIONAL_PROBE CVAR " << cvar_name << " must evaluate to DN_DATA." << endl;
    ierr = -1;
  }
  else {
    // if this is a killfile version of this init, leave the scp.cvar.second connected to the CtiData, even though it is
    // not registered data. It is just going to be "one-and-done" anyways. For the case of registered data,
    // the pointer is immutable so this saves time...
    if (one_off||(CtiRegister::isRegisteredData(cvar_data))) {
            scp.cvar.first = cvar_name;
            scp.cvar.second = cvar_data;
    }
    else {
            scp.cvar.first = cvar_name;
            scp.cvar.second = NULL; // use NULL to indicate this in not registered and needs eval...
    }

    // CVAR_2 type check...
    if (cvar_2_name != "") {
      CtiRegister::CtiData * cvar_2_data = CtiRegister::getCtiData(cvar_2_name);
      if ((cvar_2_data == NULL)||(cvar_2_data->getType() != DN_DATA)||(cvar_2_data->getUnindexedTopology() != cvar_data->getUnindexedTopology())) {
        if ( mpi_rank == 0 )
          cout << "Warning: SMD_CONDITIONAL_PROBE CVAR_2 " << cvar_2_name << " must evaluate to same topology/type as CVAR." << endl;
        ierr = -1;
      }
      else {
        if (one_off||(CtiRegister::isRegisteredData(cvar_2_data))) {
          scp.cvar_2.first = cvar_2_name;
          scp.cvar_2.second = cvar_2_data;
        }
        else {
          scp.cvar_2.first = cvar_2_name;
          scp.cvar_2.second = NULL; // use NULL to indicate this in not registered and needs eval...
        }
      }
    }

    // WVAR type check...
    if (wvar_name != "") {
      CtiRegister::CtiData * wvar_data = CtiRegister::getCtiData(wvar_name);
      if ((wvar_data == NULL)||(wvar_data->getType() != DN_DATA)||(wvar_data->getUnindexedTopology() != cvar_data->getUnindexedTopology())) {
        if ( mpi_rank == 0 )
          cout << "Warning: SMD_CONDITIONAL_PROBE WVAR " << wvar_name << " must evaluate to same topology/type as CVAR." << endl;
        ierr = -1;
      }
      else {
        if (one_off||(CtiRegister::isRegisteredData(wvar_data))) {
          scp.wvar.first = wvar_name;
          scp.wvar.second = wvar_data;
        }
        else {
          scp.wvar.first = wvar_name;
          scp.wvar.second = NULL; // use NULL to indicate this in not registered and needs eval...
        }
      }
    }

    if (ierr != -1) {
      // VARS type check...
      assert( scp.var_vec.empty());
      for (int ii = 0; ii < int(var_vec.size()); ++ii) {
        CtiRegister::CtiData * data = CtiRegister::getCtiData(var_vec[ii]);

        if ((data == NULL)||(data->getUnindexedTopology() != cvar_data->getUnindexedTopology())||(data->getType() != DN_DATA)) {
          if ( mpi_rank == 0 )
            cout << "Warning: SMD_CONDITIONAL_PROBE VAR \"" << var_vec[ii] << "\" does not evaluate to same topology/type as CVAR. Skipping." << endl;
        }

        else {
          // see note above on killfile...
          if (one_off||(CtiRegister::isRegisteredData(data))) {
            scp.var_vec.push_back(pair<string,CtiRegister::CtiData*>(var_vec[ii],data));
          }
          else {
            scp.var_vec.push_back(pair<string,CtiRegister::CtiData*>(var_vec[ii],NULL));
          }
        }
      }

      if ( scp.var_vec.empty()) {
        if ( mpi_rank == 0 )
          cout << "Warning: SMD_CONDITIONAL_PROBE has no valid data." << endl;
        ierr = -1;
      }
    }

  }

  if (ierr != 0)
    return ierr;

  // complete initialization...

  /*
  *cp.hss << "# 1:bin 2:count 3:" << cvar_name << "_avg 4:" << cvar_name << "_rms";
  for (int i = 0; i < int(cp.var_vec.size()); ++i)
    *cp.hss << " " << 5+i*4 << ":" << cp.var_vec[i].first << "_avg " <<
      6+i*4 << ":" << cp.var_vec[i].first << "_rms " <<
      7+i*4 << ":" << cp.var_vec[i].first << "_min " <<
      8+i*4 << ":" << cp.var_vec[i].first << "_max";
  *cp.hss << endl;
  */

  MiscUtils::mkdir_for_file_collective(scp.name);

  return 0;

}


int StaticSolver::initPointCloudProbe(PointCloudProbe* pcp,Param* param) {

  if (mpi_rank == 0)
    cout << "initPointCloudProbe() \"" << param->str() << "\"..." << endl;

  // everyone gets the point cloud probe because everyone particpates in the striped write. the probe data
  // are parsed on the processor which has the cv they sit in. then the data is sent to striped distribution.

  int geom = -1;
  string filename = "";
  int iarg = 0;
  set<string> varSet;
  bool write_readme = true;
  while ( iarg < param->size()) {

    string token = param->getString(iarg++);

    if ( token == "NAME" ) {
      pcp->name = param->getString(iarg++);
    }
    else if (token == "INTERVAL") {
      pcp->interval = param->getInt(iarg++);
      if (pcp->interval <= 0) {
        CWARN(" > INTERVAL expects a positive integer; treating as one-off request");
        // invalid value entered, so treat as unset (default value)
        pcp->interval = -1;
      }
    }
    else if ( token == "GEOM" ) {
      token = param->getString(iarg++);
      if (token == "FILE") {
        geom = FILE_GEOM;
        filename = param->getString(iarg++);
      }
      // should add BOX, CYLINDER, TCONE, etc...
      else {
        CWARN("unrecognized PROBE GEOM: " << token << ". Skipping.");
      }
    }
    else if (token == "PRECISION") {
      token = param->getString(iarg++);
      if ((token == "FLOAT")||(token == "FLOAT_PRECISION")||(token == "SINGLE")||(token == "SINGLE_PRECISION")) {
        pcp->precision = FLOAT_PRECISION;
      }
      else if (token == "DOUBLE" || token == "DOUBLE_PRECISION") {
        pcp->precision = DOUBLE_PRECISION;
      }
    }
    else if ( token == "FORMAT" ) {
      token = param->getString(iarg++);
      if (token == "ASCII")
        pcp->b_ascii = true;
      else if (token == "BINARY")
        pcp->b_ascii = false;
    }
    else if ( token == "VARS") {
      // the rest are vars...
      while (iarg < param->size()) {
        varSet.insert(param->getString(iarg++));
      }
    }
    else if ( token == "NO_README") {
      write_readme = false;
    }
    else {
      if ( mpi_rank == 0)
        cout << "Warning: POINTCLOUD_PROBE skipping unrecognized token " << token << endl;
    }
  }

  int ierr = 0;

  // NAME...

  if ( pcp->name == "") {

    if ( mpi_rank == 0 )
      cout << "Warning: POINTCLOUD_PROBE missing NAME." << endl;

    ierr = -1;
  }

  // GEOM...

  if (geom <= 0) {

    if ( mpi_rank == 0 )
      cout << "Warning: POINTCLOUD_PROBE GEOM invalid: " << geom << endl;

    ierr = -1;
  }
  else if (geom == FILE_GEOM && filename == "") {

    if ( mpi_rank == 0 )
      cout << "Warning: POINTCLOUD_PROBE GEOM FILE name invalid: " << filename << endl;

    ierr = -1;
  }

  // VARS...

  if (varSet.size() == 0) {
    if ( mpi_rank == 0 )
      cout << "Warning: POINTCLOUD_PROBE VARS empty." << endl;

    ierr = -1;
  }

  if ( ierr != 0) {
    if ( mpi_rank == 0 )
      cout << " Warning: skipping this POINTCLOUD_PROBE. Correct syntax. " << endl;
    return -1;
  }

  // get cti data ptrs to requested vars...

  for (set<string>::iterator it = varSet.begin(); it != varSet.end(); ++it) {
    CtiRegister::CtiData * data = CtiRegister::getRegisteredCtiData(*it,false);
    if ( data != NULL) {

      if ( ((data->getUnindexedTopology() != CV_DATA)&&(data->getUnindexedTopology() != NO_DATA)) ||
           (data->getType() != DN_DATA)) {
        if (mpi_rank == 0 )
          cout << " > Warning : POINTCLOUD_PROBE var does not evaluate to CV_DN,NO_DN: " << *it << endl;
        continue;
      }

      pcp->var_vec.push_back(pair<string,CtiRegister::CtiData*>(*it, data));

    } else {

      // unregistered data so it needs to be evaluated..
      data = CtiRegister::getUnregisteredCtiData(*it);
      if ( (data == NULL) ||
           ((data->getUnindexedTopology() != CV_DATA)&&(data->getUnindexedTopology() != NO_DATA)) ||
           (data->getType() != DN_DATA) ) {
        if (mpi_rank == 0 ) {
          cout << " Warning: POINTCLOUD_PROBE var does not evaluate to CV_DN,NO_DN: " << *it << endl;
        }
        continue;
      }

      pcp->var_vec.push_back(pair<string,CtiRegister::CtiData*>(*it,NULL));
    }
  }

  // get points...

  vector<double> x_vec;
  pcp->readPointsAscii(x_vec,filename);

  assert(x_vec.size()%3 == 0);
  const int npt0 = x_vec.size()/3;
  int8 my_npt0 = (int8)npt0;
  int8 ipt0_offset;
  MPI_Scan(&my_npt0,&ipt0_offset,1,MPI_INT8,MPI_SUM,mpi_comm);
  //if (mpi_rank == mpi_size-1) cout << " > number of read points: " << ipt0_offset << endl;
  ipt0_offset -= my_npt0;

  // ------------------------------------------------------------------------
  // now send the probe points to the ranks owning the cv's that contain them
  // ------------------------------------------------------------------------

  // build the bbox...

  // make sure we have the stuff we need...

  if (cvAdt == NULL) buildCvAdt();

  // get all periodicity combinations...
  vector<int> bitVec;
  PeriodicData::buildBitVec(bitVec);

  int * send_count = new int[mpi_size];
  FOR_RANK send_count[rank] = 0;

  vector<int> intVec;
  for (int ipt0 = 0; ipt0 < npt0; ++ipt0) {
    for (int ibv = 0; ibv < bitVec.size(); ++ibv) {
      const int bits = bitVec[ibv];
      double xp_t[3]; FOR_I3 xp_t[i] = x_vec[ipt0*3+i];
      if (bits != 0) PeriodicData::periodicTranslate(xp_t,1,bits);
      assert(intVec.empty());
      cvBboxAdt->buildListForPoint(intVec,xp_t);
      for (int ii = 0,ii_end = intVec.size(); ii < ii_end; ++ii) {
        const int rank = intVec[ii]; assert((rank >= 0)&&(rank < mpi_size));
        send_count[rank] += 3;
      }
      intVec.clear();
    }
  }

  int * send_disp = new int[mpi_size];
  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
  const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

  double * send_buf = new double[send_count_sum];
  assert(send_count_sum%3 == 0);
  int8 *send_buf_int8 = new int8[send_count_sum/3];
  for (int ipt0 = 0; ipt0 < npt0; ++ipt0) {
    for (int ibv = 0; ibv < bitVec.size(); ++ibv) {
      const int bits = bitVec[ibv];
      double xp_t[3]; FOR_I3 xp_t[i] = x_vec[ipt0*3+i];
      if (bits != 0) PeriodicData::periodicTranslate(xp_t,1,bits);
      assert(intVec.empty());
      cvBboxAdt->buildListForPoint(intVec,xp_t);
      for (int ii = 0, ii_end=intVec.size(); ii < ii_end; ++ii) {
        const int rank = intVec[ii]; assert((rank >= 0)&&(rank < mpi_size));
        FOR_I3 send_buf[send_disp[rank]+i] = xp_t[i];
        send_buf_int8[send_disp[rank]/3] = ipt0+ipt0_offset;
        send_disp[rank] += 3;
      }
      intVec.clear();
    }
  }
  x_vec.clear();

  // rewind...

  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

  // now send...

  int * recv_count = new int[mpi_size];
  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

  int * recv_disp = new int[mpi_size];
  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
  const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

  double * recv_buf = new double[recv_count_sum];
  MPI_Alltoallv(send_buf,send_count,send_disp,MPI_DOUBLE,
                recv_buf,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);

  FOR_RANK {
    send_count[rank] /= 3;
    send_disp[rank] /= 3;
    recv_count[rank] /= 3;
    recv_disp[rank] /= 3;
  }

  assert(recv_count_sum%3 == 0);
  int8 * recv_buf_int8 = new int8[recv_count_sum/3];
  MPI_Alltoallv(send_buf_int8,send_count,send_disp,MPI_INT8,
                recv_buf_int8,recv_count,recv_disp,MPI_INT8,mpi_comm);

  // now unpack and set...

  double* recv_buf2 = new double[recv_count_sum/3]; // holds dist2
  int * flag = new int[recv_count_sum/3]; // -1
  for (int ipt = 0; ipt < recv_count_sum/3; ++ipt) flag[ipt] = -1; // when flag[ipt] >= 0, x_pt[ipt] lies on this rank
  for (int irecv = 0; irecv < recv_count_sum; irecv += 3) {
    const int ipt = irecv/3;
    double xp[3]; FOR_I3 xp[i] = recv_buf[irecv+i];
    assert(intVec.empty());
    cvAdt->buildListForPoint(intVec,xp);
    double dist_sq = HUGE_VAL;
    for (int ii = 0, ii_end =intVec.size(); ii < ii_end; ++ii) {
      const int this_icv = intVec[ii]; assert((this_icv >= 0)&&(this_icv < ncv));
      const double this_dist_sq = DIST2(xp,x_vv[this_icv]);
      if (this_dist_sq <= dist_sq && this_dist_sq <= (1.0+1.0E-6)*r_vv[this_icv]*r_vv[this_icv]) {
        flag[ipt] = this_icv;
        dist_sq = this_dist_sq;
      }
    }
    intVec.clear();
    recv_buf2[irecv/3] = dist_sq; // for tie breaking...
  }

  // tie breaking and reindexing..

  MPI_Alltoallv(recv_buf2,recv_count,recv_disp,MPI_DOUBLE,
                send_buf,send_count,send_disp,MPI_DOUBLE,mpi_comm);

  int *ipt0_new = new int[npt0];
  for (int ipt0 = 0; ipt0 < npt0; ++ipt0) {
    ipt0_new[ipt0] = -1;
  }

  // tie breaking... 
  
  for (int isend = 0; isend < send_count_sum/3; ++isend) {
    const int ipt0 = send_buf_int8[isend]-ipt0_offset;
    if (send_buf[isend] < HUGE_VAL) {
      if (ipt0_new[ipt0] == -1) {
        ipt0_new[ipt0] = isend;
      }
      else {
        const int isend0 = ipt0_new[ipt0];
        assert((isend0 >= 0)&&(isend0 < send_count_sum/3));
        if (send_buf[isend] < send_buf[isend0]) {
          send_buf[isend0] = HUGE_VAL; 
          ipt0_new[ipt0] = isend;
        }
        else {
          send_buf[isend] = HUGE_VAL; 
        }
      }
    }
  }

  // reindex ipt0 excluding points outside domain...
  
  int8 my_npt0_new = 0;
  for (int ipt0 = 0; ipt0 < npt0; ++ipt0) {
    if (ipt0_new[ipt0] >= 0)
      ipt0_new[ipt0] = my_npt0_new++;
  }
  int8 ipt0_offset_new;
  MPI_Scan(&my_npt0_new,&ipt0_offset_new,1,MPI_INT8,MPI_SUM,mpi_comm);
  ipt0_offset_new -= my_npt0_new;
  for (int isend = 0; isend < send_count_sum/3; ++isend) {
    if (send_buf[isend] < HUGE_VAL) {
      const int ipt0 = send_buf_int8[isend]-ipt0_offset;
      assert((ipt0 >= 0)&&(ipt0 < npt0));
      assert(ipt0_new[ipt0] >= 0);
      send_buf_int8[isend] = ipt0_new[ipt0]+ipt0_offset_new;
    }
  }
  delete[] ipt0_new;

  MPI_Alltoallv(send_buf,send_count,send_disp,MPI_DOUBLE,
                recv_buf2,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
  delete[] send_buf;

  // send updated index back to ranks containing them...
  MPI_Alltoallv(send_buf_int8,send_count,send_disp,MPI_INT8,
                recv_buf_int8,recv_count,recv_disp,MPI_INT8,mpi_comm);
  delete[] send_buf_int8;
  delete[] send_disp;
  delete[] send_count;
  delete[] recv_disp;
  delete[] recv_count;

  pcp->npt = 0;
  for (int irecv = 0; irecv < recv_count_sum/3; ++irecv) {
    if (recv_buf2[irecv] < HUGE_VAL) {
      assert(flag[irecv] >= 0);
      ++pcp->npt;
    }
  }
  pcp->x_pt = new double[pcp->npt][3];
  pcp->icv_pt = new int[pcp->npt];
  pcp->ipt0_pt = new int8[pcp->npt];
  int ipt = 0;
  for (int irecv = 0; irecv < recv_count_sum; irecv += 3) {
    const int ind = irecv/3;
    if (recv_buf2[ind] < HUGE_VAL) {
      FOR_I3 pcp->x_pt[ipt][i] = recv_buf[ind*3+i];
      pcp->icv_pt[ipt] = flag[ind];
      pcp->ipt0_pt[ipt] = recv_buf_int8[ind];
      ipt++;
    }
  }
  assert(ipt <= recv_count_sum/3);
  assert(ipt == pcp->npt);

  delete[] flag;
  delete[] recv_buf;
  delete[] recv_buf2;
  delete[] recv_buf_int8;

  // build striping for the write...

  pcp->initDdeStuff();

  if (mpi_rank == 0) cout << " > point cloud probe npt_global: " << pcp->npt_global << endl;
  MiscUtils::dumpRange(pcp->x_pt,pcp->npt,"point cloud probe x_pt");

  // write points...

  if (pcp->b_ascii) {
    pcp->writePointsAscii();
  }
  else {
    pcp->writePointsBinary();

    // write read me...

    if (write_readme)
      pcp->writeReadMe(param,filename);
  }
  delete[] pcp->ipt0_pt; pcp->ipt0_pt = NULL; // dont need it anymore (for now)...



  // build node weights...

  assert(pcp->wtopt_i == NULL);
  assert(pcp->wtopt_v == NULL);

  pcp->wtopt_i = new int[pcp->npt+1];
  pcp->wtopt_i[0] = 0;
  for (int ipt = 0; ipt < pcp->npt; ++ipt) {
    const int icv = pcp->icv_pt[ipt];
    set<int> no_list;
    for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
      const int ifa = faocv_v[foc];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        no_list.insert(noofa_v[nof]);
      }
    }
    for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
      const int ibf = bfocv_v[boc];
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        no_list.insert(noobf_v[nob]);
      }
    }
    pcp->wtopt_i[ipt+1] = pcp->wtopt_i[ipt]+no_list.size();
  }
  int nwt = pcp->wtopt_i[pcp->npt];

  pcp->wtopt_v = new pair<int,double>[nwt];
  nwt = 0;
  double my_max_dlambda = 0.0;
  for (int ipt = 0; ipt < pcp->npt; ++ipt) {
    const int icv = pcp->icv_pt[ipt];
    map<const int,int> nodeMap;
    for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
      const int ifa = faocv_v[foc];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino = noofa_v[nof];
        if (nodeMap.find(ino) == nodeMap.end()) {
          nodeMap[ino] = nwt;
          pcp->wtopt_v[nwt].first = ino;
          //const double dist = DIST(x_no[ino],pcp->x_pt[ipt]);
          //assert(dist > 0.0);
          //pcp->wtopt_v[nwt].second = 1.0/dist; // inverse distance weighted (normalized below)
          pcp->wtopt_v[nwt].second = 0.0;
          nwt++;
        }
      }
    }
    for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
      const int ibf = bfocv_v[boc];
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino = noobf_v[nob];
        if (nodeMap.find(ino) == nodeMap.end()) {
          nodeMap[ino] = nwt;
          pcp->wtopt_v[nwt].first = ino;
          //const double dist = DIST(x_no[ino],pcp->x_pt[ipt]);
          //assert(dist > 0.0);
          //pcp->wtopt_v[nwt].second = 1.0/dist; // inverse distance weighted (normalized below)
          pcp->wtopt_v[nwt].second = 0.0;
          nwt++;
        }
      }
    }

    // find the sub-tet where the point lies. this should be where all lambdas are [0,1]...

    int tet[4] = {-1,-1,-1,-1}; // icv, ifa or -ibf-1, ino0, ino1
    double lambda[4]; // lambda_cv,lambda_fa,lambda_ino0,lambda_ino1
    // some times we need to choose the "best" candidate which does not strictly satisfy the inside check
    double dlambda = HUGE_VAL; // how bad our best candidate is
    bool b_tet = false;
    for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
      const int ifa = faocv_v[foc];
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        double this_lambda[4];
        if (MiscUtils::getBarycentricCoordinates(this_lambda,pcp->x_pt[ipt],x_cv[icv],x_fa[ifa],x_no[ino0],x_no[ino1])) {
          dlambda = 0.0;
          FOR_I4 lambda[i] = this_lambda[i];
          tet[0] = icv; tet[1] = ifa; tet[2] = ino0; tet[3] = ino1;
          b_tet = true;
          break;
        }
        else {
          double this_dlambda = 0.0;
          FOR_I4 this_dlambda = max(this_dlambda,max(-this_lambda[i],this_lambda[i]-1.0));
          if (this_dlambda < dlambda) {
            dlambda = this_dlambda;
            FOR_I4 lambda[i] = this_lambda[i];
            tet[0] = icv; tet[1] = ifa; tet[2] = ino0; tet[3] = ino1;
            //b_tet = true;
            //break;
          }
        }
      }
      if (b_tet) break;
    }
    if (!b_tet) {
      for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
        const int ibf = bfocv_v[boc];
        int ino1 = noobf_v[noobf_i[ibf+1]-1];
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int ino0 = ino1;
          ino1 = noobf_v[nob];
          double this_lambda[4];
          if (MiscUtils::getBarycentricCoordinates(this_lambda,pcp->x_pt[ipt],x_cv[icv],x_bf[ibf],x_no[ino0],x_no[ino1])) {
            dlambda = 0.0;
            FOR_I4 lambda[i] = this_lambda[i];
            tet[0] = icv; tet[1] = -ibf-1; tet[2] = ino0; tet[3] = ino1;
            b_tet = true;
            break;
          }
          else {
            double this_dlambda = 0.0;
            FOR_I4 this_dlambda = max(this_dlambda,max(-this_lambda[i],this_lambda[i]-1.0));
            if (this_dlambda < dlambda) {
              dlambda = this_dlambda;
              FOR_I4 lambda[i] = this_lambda[i];
              tet[0] = icv; tet[1] = -ibf-1; tet[2] = ino0; tet[3] = ino1;
              //b_tet = true;
              //break;
            }
          }
        }
        if (b_tet) break;
      }
    }
    //assert(b_tet); // must have some tet...
    assert(tet[0] >= 0); // must have some tet
    my_max_dlambda = max(my_max_dlambda,dlambda);

    // get cv and fa coeff...

    double cv_den = 0.0;
    double fa_coeff = 0.0;
    for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
      const int ifa = faocv_v[foc];
      const double area = MAG(n_fa[ifa]);
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      double sum_dist = 0.0;
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        sum_dist += DIST(x_no[ino0],x_no[ino1]);
      }
      cv_den += area*sum_dist;
      if (tet[1] == ifa) {
        assert(sum_dist > 0.0);
        fa_coeff = lambda[1]/sum_dist;
      }
    }
    for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
      const int ibf = bfocv_v[boc];
      const double area = MAG(n_bf[ibf]);
      int ino1 = noobf_v[noobf_i[ibf+1]-1];
      double sum_dist = 0.0;
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino0 = ino1;
        ino1 = noobf_v[nob];
        sum_dist += DIST(x_no[ino0],x_no[ino1]);
      }
      cv_den += area*sum_dist;
      if (tet[1] == (-ibf-1)) {
        assert(sum_dist > 0.0);
        fa_coeff = lambda[1]/sum_dist;
      }
    }
    assert(cv_den > 0.0);
    const double cv_coeff = lambda[0]/cv_den;

    // now populate nodal weights...

    for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
      const int ifa = faocv_v[foc];
      const double area = MAG(n_fa[ifa]);
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      map<const int,int>::iterator it1 = nodeMap.find(ino1);
      assert(it1 != nodeMap.end());
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        map<const int,int>::iterator it0 = it1;
        ino1 = noofa_v[nof];
        const double half_dist = 0.5*DIST(x_no[ino0],x_no[ino1]); // edge is split b/w ino0 and ino1
        it1 = nodeMap.find(ino1);
        assert(it1 != nodeMap.end());
        // contribution of phi_cv...
        pcp->wtopt_v[it0->second].second += area*(half_dist)*cv_coeff;
        pcp->wtopt_v[it1->second].second += area*(half_dist)*cv_coeff;
        if (tet[1] == ifa) {
          // contribution of phi_fa...
          pcp->wtopt_v[it0->second].second += half_dist*fa_coeff;
          pcp->wtopt_v[it1->second].second += half_dist*fa_coeff;
          if (tet[2] == ino0) {
            assert(tet[3] == ino1);
            // contribution of phi_ino0 and phi_ino1...
            pcp->wtopt_v[it0->second].second += lambda[2];
            pcp->wtopt_v[it1->second].second += lambda[3];
          }
        }
      }
    }
    for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
      const int ibf = bfocv_v[boc];
      const double area = MAG(n_bf[ibf]);
      int ino1 = noobf_v[noobf_i[ibf+1]-1];
      map<const int,int>::iterator it1 = nodeMap.find(ino1);
      assert(it1 != nodeMap.end());
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino0 = ino1;
        map<const int,int>::iterator it0 = it1;
        ino1 = noobf_v[nob];
        const double half_dist = 0.5*DIST(x_no[ino0],x_no[ino1]); // edge is split b/w ino0 and ino1
        it1 = nodeMap.find(ino1);
        assert(it1 != nodeMap.end());
        // contribution of phi_cv...
        pcp->wtopt_v[it0->second].second += area*(half_dist)*cv_coeff;
        pcp->wtopt_v[it1->second].second += area*(half_dist)*cv_coeff;
        if (tet[1] == (-ibf-1)) {
          // contribution of phi_fa...
          pcp->wtopt_v[it0->second].second += half_dist*fa_coeff;
          pcp->wtopt_v[it1->second].second += half_dist*fa_coeff;
          if (tet[2] == ino0) {
            assert(tet[3] == ino1);
            // contribution of phi_ino0 and phi_ino1...
            pcp->wtopt_v[it0->second].second += lambda[2];
            pcp->wtopt_v[it1->second].second += lambda[3];
          }
        }
      }
    }
  }
  assert(nwt == pcp->wtopt_i[pcp->npt]);

  // report max barycentric distance outside of a chosen tet...
  double max_dlambda;
  MPI_Reduce(&my_max_dlambda,&max_dlambda,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
  if (mpi_rank == 0) {
    cout << " > max barycentric distance outside of tet (non-convexity): " << max_dlambda << endl;
  }

  // normalize inverse dist weighting...
  //for (int ipt = 0; ipt < pcp->npt; ++ipt) {
  //  double sum_wt = 0.0;
  //  for (int wop = pcp->wtopt_i[ipt]; wop != pcp->wtopt_i[ipt+1]; ++wop)
  //    sum_wt += pcp->wtopt_v[wop].second;
  //  const double inv_sum_wt = 1.0/sum_wt;
  //  for (int wop = pcp->wtopt_i[ipt]; wop != pcp->wtopt_i[ipt+1]; ++wop)
  //    pcp->wtopt_v[wop].second *= inv_sum_wt;
  //}

  return ierr;
}

void StaticSolver::writePointCloudProbe(PointCloudProbe* pcp,const int step,const double time) {
  UNUSED(time);
  if (mpi_rank == 0)
    cout << "writePointCloudProbe(): " << pcp->name << endl;

  assert(pcp->var_vec.size() > 0); // shouldn't be here if we don't have vars...

  char dummy[128]; assert(pcp->name.length() < 128);
  sprintf(dummy,"%s.%08d.pcd",pcp->name.c_str(),step);
  MiscUtils::mkdir_for_file_collective(dummy,0);
  string tmp_filename = MiscUtils::makeTmpPrefix(dummy);
  MPI_File_delete(tmp_filename.c_str(),MPI_INFO_NULL);

  MPI_File fh;
  MPI_File_open(mpi_comm,tmp_filename.c_str(),MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

  const int npt_sm = pcp->ptora_striped[mpi_rank+1]-pcp->ptora_striped[mpi_rank];

  // header...

  MPI_Offset offset = 0;
  char tmp[23+1];
  char* carray = NULL;
  if (pcp->b_ascii) {
    if (pcp->precision == DOUBLE_PRECISION) {
      if (mpi_rank == 0) {
        carray = new char[23*pcp->var_vec.size()*max(npt_sm,1)+1];
        for (int ivar = 0,lim = pcp->var_vec.size(); ivar < lim; ++ivar) {
          snprintf(tmp,23,"%22s",pcp->var_vec[ivar].first.c_str());
          if (ivar < lim-1)
            sprintf(tmp+22,"%s"," ");
          else
            sprintf(tmp+22,"%s","\n");
          for (int i = 0; i < 23; ++i) carray[23*ivar+i] = tmp[i];
        }
        MPI_File_write_at(fh,0,carray,pcp->var_vec.size()*23,MPI_CHAR,MPI_STATUS_IGNORE);
      }
      else {
        carray = new char[23*pcp->var_vec.size()*npt_sm+1];
      }
      offset += pcp->var_vec.size()*23;
    }
    else {
      if (mpi_rank == 0) {
        carray = new char[16*pcp->var_vec.size()*max(npt_sm,1)+1];
        for (int ivar = 0,lim = pcp->var_vec.size(); ivar < lim; ++ivar) {
          snprintf(tmp,16,"%15s",pcp->var_vec[ivar].first.c_str());
          if (ivar < lim-1)
            sprintf(tmp+15,"%s"," ");
          else
            sprintf(tmp+15,"%s","\n");
          for (int i = 0; i < 16; ++i) carray[16*ivar+i] = tmp[i];
        }
        MPI_File_write_at(fh,0,carray,pcp->var_vec.size()*16,MPI_CHAR,MPI_STATUS_IGNORE);
      }
      else {
        carray = new char[16*pcp->var_vec.size()*npt_sm+1];
      }
      offset += pcp->var_vec.size()*16;
    }
  }
  else {
    if (mpi_rank == 0) {
      int8 ibuf[5] = { POINTS_IO_MAGIC_NUMBER, POINTS_IO_VERSION, pcp->npt_global, int8(pcp->var_vec.size()), pcp->precision};
      MPI_File_write_at(fh,0,ibuf,5,MPI_INT8,MPI_STATUS_IGNORE);
    }
    offset += 5*int8_size;
  }

  // scalar data (user responsible for breaking data into components)...

  double *dn_no = new double[nno];
  double *dn_pt = new double[pcp->npt];
  double *dn_pt_sm = new double[npt_sm];
  float *fn_pt_sm = NULL;
  if ((pcp->precision == FLOAT_PRECISION)&&(!pcp->b_ascii)) fn_pt_sm = new float[npt_sm];
  for (int ivar = 0,lim = pcp->var_vec.size(); ivar < lim; ++ivar) {

    // average the requested data to the nodes...

    if (pcp->var_vec[ivar].second != NULL) {
      // previously evaluated to registered data
      setNoDN(dn_no,pcp->var_vec[ivar].first,pcp->var_vec[ivar].second);
    }
    else {
      // needs to be evaluated
      setNoDN(dn_no,pcp->var_vec[ivar].first);
    }

    // now use inverse distance weighting to get values at the points...

    for (int ipt = 0; ipt < pcp->npt; ++ipt) {
      dn_pt[ipt] = 0.0;
      for (int wop = pcp->wtopt_i[ipt]; wop != pcp->wtopt_i[ipt+1]; ++wop) {
        const int ino = pcp->wtopt_v[wop].first;
        assert((ino >= 0)&&(ino < nno));
        const double wgt = pcp->wtopt_v[wop].second;
        //assert((wgt > 0.0)&&(wgt < 1.0)); // not convex anymore
        dn_pt[ipt] += wgt*dn_no[ino];
      }
    }

    // push to striped...
    pcp->dde_striped->push(dn_pt_sm,dn_pt);

    // write...

    if (pcp->precision == DOUBLE_PRECISION) {
      if (pcp->b_ascii) {
        for (int ipt = 0; ipt < npt_sm; ++ipt) {
          if (ivar < lim-1)
            sprintf(tmp,"% 18.15e ",dn_pt_sm[ipt]);
          else
            sprintf(tmp,"% 18.15e\n",dn_pt_sm[ipt]);
          for (int i = 0; i < 23; ++i) carray[23*(pcp->var_vec.size()*ipt+ivar)+i] = tmp[i];
        }
      }
      else {
        writeChunkedData<double>(fh,offset+pcp->ptora_striped[mpi_rank]*double_size,dn_pt_sm,npt_sm,mpi_comm);
        offset += pcp->npt_global*double_size;
      }
    }
    else {
      assert(pcp->precision == FLOAT_PRECISION);
      if (pcp->b_ascii) {
        for (int ipt = 0; ipt < npt_sm; ++ipt) {
          if (ivar < lim-1)
            sprintf(tmp,"% 12.8e ",(float)dn_pt_sm[ipt]);
          else
            sprintf(tmp,"% 12.8e\n",(float)dn_pt_sm[ipt]);
          for (int i = 0; i < 16; ++i) carray[16*(pcp->var_vec.size()*ipt+ivar)+i] = tmp[i];
        }
      }
      else {
        for (int ipt = 0; ipt < npt_sm; ++ipt) fn_pt_sm[ipt] = (float)dn_pt_sm[ipt];
        writeChunkedData<float>(fh,offset+pcp->ptora_striped[mpi_rank]*sizeof(float),fn_pt_sm,npt_sm,mpi_comm);
        offset += pcp->npt_global*sizeof(float);
      }
    }

  }
  if (pcp->b_ascii) {
    if (pcp->precision == DOUBLE_PRECISION)
      writeChunkedData<char>(fh,offset+pcp->ptora_striped[mpi_rank]*pcp->var_vec.size()*23,carray,pcp->var_vec.size()*23*npt_sm,mpi_comm);
    else
      writeChunkedData<char>(fh,offset+pcp->ptora_striped[mpi_rank]*pcp->var_vec.size()*16,carray,pcp->var_vec.size()*16*npt_sm,mpi_comm);
    delete[] carray;
  }
  delete[] dn_no;
  delete[] dn_pt;
  delete[] dn_pt_sm;
  DELETE(fn_pt_sm);

  // and close...

  MPI_File_close(&fh);
  if (mpi_rank == 0) {
    remove(dummy);
    rename(tmp_filename.c_str(),dummy);
  }

}

int StaticSolver::initFwhSurface(FwhSurface& fs,Param * param,const bool killfile) {

  // if killfile == true, then this param is being parsed from a killfile. If
  // false, then this param is one of the input file parameters. This affects the behavior
  // of certain options...

  if (mpi_rank == 0)
    cout << "[FWH_SURFACE] initFwhSurface() \"" << param->str() << "\"..." << endl;

  bool b_name = false;
  vector<SimpleGeom*> simpleGeomVec;
  // endcaps...
  double x_ec[3],dx_ec[3];
  int n_ec = 1;
  // tecplot output of geom...
  bool b_tecplot = false;
  bool var_mode = false;
  vector<string> var_names;
  vector<int> zone_indices;

  int ierr = 0;
  int iarg = 0;
  while (iarg < param->size()) {
    string token = param->getString(iarg++);
    if (token == "NAME") {
      fs.name = param->getString(iarg++);
      // make sure the name is unique...
      int match_count;
      int name_index = 0;
      const string base_name = fs.name;
      do {
        match_count = 0;
        for (map<const string,FwhSurface>::iterator iter = fwhSurfaceMap.begin(); iter != fwhSurfaceMap.end(); ++iter) {
          if ((iter->second.status != CTI_STATUS_ERROR)&&(iter->second.name == fs.name))
            ++match_count;
        }
        if (match_count != 1) {
          ++name_index;
          // non const lvalue error 
          //fs.name = base_name+"_"+SSTR(name_index);
          stringstream tmp_ss;
          tmp_ss << name_index;
          fs.name = base_name+"_"+tmp_ss.str();
        }
      }
      while (match_count != 1);
      if (fs.name != base_name) {
        if (mpi_rank == 0) cout << "[FWH_SURFACE] Warning: Duplicate name found. Name modified to: " << fs.name << endl;
      }
      b_name = true;
      var_mode = false;
    }
    else if (token == "GEOM") {
      if (iarg == param->size()) {
        CERR("GEOM missing parameters");
      }
      else {
        token = MiscUtils::toUpperCase(param->getString(iarg));
        if ((token == "FAZONE")||(token == "BFZONE")||(token == "ZONE")) {
          ++iarg;
          const string zonesCsv = param->getString(iarg++);
          vector<string> zonesVec;
          MiscUtils::splitCsv(zonesVec,zonesCsv);
          for (vector<string>::iterator it=zonesVec.begin(); it!=zonesVec.end(); ++it) {
            map<const string,int>::iterator iter = bfZoneNameMap.find(*it);
            if (iter != bfZoneNameMap.end()) {
              const int izone = iter->second;
              zone_indices.push_back(izone);
            }
            else {
              CWARN(" > skipping invalid zone name: " << *it);
            }
          }
        }
        else {
          SimpleGeom * sg = newSimpleGeom(&(*param),iarg);
          if (sg == NULL) ierr = -1;
          simpleGeomVec.push_back(sg);
        }
        var_mode = false;
      }
    }
    else if ((token == "TEC")||(token == "TECPLOT")) {
      b_tecplot = true;
      var_mode = false;
    }
    else if (token == "INTERVAL") {
      fs.interval = param->getInt(iarg++);
      if (fs.interval < 1) {
        if (mpi_rank == 0) cout << "[FWH_SURFACE] Warning: INTERVAL must be >= 1" << endl;
        ierr = -1;
      }
      var_mode = false;
    }
    else if (token == "ENDCAPS_X") {
      // ENDCAP_X n x0 x1
      assert(n_ec == 1);
      n_ec = param->getInt(iarg++);
      assert(n_ec > 1);
      const double x0 = param->getDouble(iarg++);
      const double x1 = param->getDouble(iarg++);
      // set the general endcap data...
      x_ec[0] = x0;
      x_ec[1] = 0.0;
      x_ec[2] = 0.0;
      dx_ec[0] = (x1-x0)/double(n_ec-1);
      dx_ec[1] = 0.0;
      dx_ec[2] = 0.0;
      var_mode = false;
    }
    else if (token == "ENDCAPS_Y") {
      // ENDCAP_Y n y0 y1
      assert(n_ec == 1);
      n_ec = param->getInt(iarg++);
      assert(n_ec > 1);
      const double y0 = param->getDouble(iarg++);
      const double y1 = param->getDouble(iarg++);
      // set the general endcap data...
      x_ec[0] = 0.0;
      x_ec[1] = y0;
      x_ec[2] = 0.0;
      dx_ec[0] = 0.0;
      dx_ec[1] = (y1-y0)/double(n_ec-1);
      dx_ec[2] = 0.0;
    }
    else if (token == "ENDCAPS_Z") {
      // ENDCAP_X n z0 z1
      assert(n_ec == 1);
      n_ec = param->getInt(iarg++);
      assert(n_ec > 1);
      const double z0 = param->getDouble(iarg++);
      const double z1 = param->getDouble(iarg++);
      // set the general endcap data...
      x_ec[0] = 0.0;
      x_ec[1] = 0.0;
      x_ec[2] = z0;
      dx_ec[0] = 0.0;
      dx_ec[1] = 0.0;
      dx_ec[2] = (z1-z0)/double(n_ec-1);

    }
    else if (token == "ENDCAPS") {
      // Note: this arbitrary encap can be difficult to specify because
      // it requires x0,y0,z0 to be computed along the axis of the surface.
      // There are probably better ways to do this -- revisit when
      // required.
      // ENDCAP n x0 y0 z0 x1 y1 z1
      assert(n_ec == 1);
      n_ec = param->getInt(iarg++);
      assert(n_ec > 1);
      const double x0 = param->getDouble(iarg++);
      const double y0 = param->getDouble(iarg++);
      const double z0 = param->getDouble(iarg++);
      const double x1 = param->getDouble(iarg++);
      const double y1 = param->getDouble(iarg++);
      const double z1 = param->getDouble(iarg++);
      // set the general endcap data...
      x_ec[0] = x0;
      x_ec[1] = y0;
      x_ec[2] = z0;
      dx_ec[0] = (x1-x0)/double(n_ec-1);
      dx_ec[1] = (y1-y0)/double(n_ec-1);
      dx_ec[2] = (z1-z0)/double(n_ec-1);
    }
    else if ((token == "VARS")||(token == "VAR")) {
      var_mode = true;
    }
    else if (var_mode) {
      // var_mode means we have just read the "VAR" param, and we do not
      // recognize the next param as a keyword, so it must be a var...
      var_names.push_back(token);
    }
    else {
      if (mpi_rank == 0) cout << "[FWH_SURFACE] Warning: unrecognized token: " << token << endl;
    }
  }

  if (!b_name) {
    if (mpi_rank == 0) cout << "[FWH_SURFACE] Warning: FWH_SURFACE missing NAME" << endl;
    ierr = -1;
  }

  if (simpleGeomVec.empty()&&zone_indices.empty()) {
    if (mpi_rank == 0) cout << "Warning: FWH_SURFACE missing GEOM" << endl;
    ierr = -1;
  }

  if (fs.interval == -1) {
    if (mpi_rank == 0) cout << "Warning: FWH_SURFACE missing INTERVAL" << endl;
    ierr = -1;
  }

  if (var_names.empty()) {
    // for now, if no vars are specified, we introduce a default behavior here. If
    // GEOM = FAZONE <name1>,<name2>,<name3> then we just add "p". Otherwise we expect
    // one of more simpleGeom's in simpleGeomVec, and we add p and u...
    // TODO: in the future, this should call a virtual function?...
    var_names.push_back("p");
    // when we have zones specified we are solid and do not need u
    if (zone_indices.empty())
      var_names.push_back("u");
  }

  assert(fs.var_vec.empty());
  for (int ii = 0; ii < int(var_names.size()); ++ii) {
    CtiRegister::CtiData * data = CtiRegister::getCtiData(var_names[ii]);
    if ( (data == NULL) || (data->getUnindexedTopology() != CV_DATA) || !((data->getType() == DN_DATA)||(data->getType() == DN3_DATA)) ) {
      if (mpi_rank == 0)
        cout << "Warning: FWH_SURFACE VAR \"" << var_names[ii] << "\" does not evaluate to CV DN or CV DN3 data." << endl;
      ierr = -1;
    }
    else {
      // if the killfile is present, we are going to be using the data right away,
      // so keep the reference. Otherwise, the reference is no good, and we
      // keep NULL and lookup every time...
      if (killfile||(CtiRegister::isRegisteredData(data))) {
        fs.var_vec.push_back(pair<string,CtiRegister::CtiData*>(var_names[ii],data));
      }
      else {
        fs.var_vec.push_back(pair<string,CtiRegister::CtiData*>(var_names[ii],NULL));
      }
    }
  }

  if (fs.var_vec.empty()) {
    if (mpi_rank == 0)
      cout << "Warning: FWH_SURFACE has no valid data." << endl;
    ierr = -1;
  }

  if (ierr != 0)
    return -1;

  int nfa_fwh = 0;
  int * count = NULL;
  int * disp = NULL;
  int8 * my_ifa_fwh_global = NULL;
  int * cv_flag = NULL;
  vector<pair<int,double> > ifaVec;
  if (!zone_indices.empty()) {
    assert(simpleGeomVec.empty());

    // loop through zones and set to 1 when next to boundary...
    for (int ii = 0, lim = zone_indices.size(); ii < lim; ++ii) {
      const int izone = zone_indices[ii];
      nfa_fwh += bfZoneVec[izone].nbf;
    }

    // build the stuff needed to gather the faces at one location...

    if (mpi_rank == 0) count = new int[mpi_size];
    MPI_Gather(&nfa_fwh,1,MPI_INT,count,1,MPI_INT,0,mpi_comm);

    fs.nfa = -1;
    if (mpi_rank == 0) {
      disp = new int[mpi_size];
      disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        disp[rank] = disp[rank-1] + count[rank-1];
      fs.nfa = disp[mpi_size-1] + count[mpi_size-1];
    }

    // now gather the global face indices to rank 0...

    my_ifa_fwh_global = new int8[nfa_fwh];
    int ind = 0;
    for (int ii = 0, lim = zone_indices.size(); ii < lim; ++ii) {
      const int izone = zone_indices[ii];
      for (int ibf = 0; ibf < bfZoneVec[izone].nbf; ++ibf)
        my_ifa_fwh_global[ind++] = bfZoneVec[izone].ibf_global[ibf];
    }
    assert(ind == nfa_fwh);

  }
  else {
    assert(zone_indices.empty());

    cv_flag = new int[ncv_g];
    FOR_ICV_G cv_flag[icv] = 0;

    // loop through GEOM's and set to 1 when inside...
    for (int ii = 0; ii < simpleGeomVec.size(); ++ii) {
      FOR_ICV if (cv_flag[icv] == 0) {
        if (simpleGeomVec[ii]->pointIsInside(x_cv[icv])) {
          cv_flag[icv] = 1;
        }
      }
      delete simpleGeomVec[ii];
    }
    simpleGeomVec.clear();

    // now the cv's inside the outer FWH surface are flagged. The following code
    // changes some of these 1's to 2,3,4... depending on how many nested surfaces
    // a given point is inside. At present, this implementation assumes arbitrary planes
    // (presumably perpendicular to the main wake/flow direction), but could be easily
    // extended to handle side wakes, etc.

    if (n_ec > 1) {

      // user has requested endcaps...

      for (int i_ec = 0; i_ec < n_ec-1; ++i_ec) {
        double this_x_ec[3]; FOR_I3 this_x_ec[i] = x_ec[i] + double(i_ec)*dx_ec[i];
        FOR_ICV {
          if (cv_flag[icv] != 0) {
            const double dx[3] = DIFF(x_cv[icv],this_x_ec);
            if (DOT_PRODUCT(dx,dx_ec) < 0.0) {
              // this is a point behind this endcap, so add 1...
              cv_flag[icv] += 1;
            }
          }
        }
      }

    }

    updateCvData(cv_flag);

    for (int ifa = 0; ifa < nfa; ++ifa) {
      if (getIfaGlobal(ifa) >= 0) {
        const int delta = cv_flag[cvofa[ifa][0]] - cv_flag[cvofa[ifa][1]];
        if (delta != 0) {
          const double this_wgt = double(delta)/double(n_ec);
          ifaVec.push_back(pair<int,double>(ifa,this_wgt));
        }
      }
    }

    nfa_fwh = ifaVec.size();

    // build the stuff needed to gather the faces at one location...

    if (mpi_rank == 0) count = new int[mpi_size];
    MPI_Gather(&nfa_fwh,1,MPI_INT,count,1,MPI_INT,0,mpi_comm);

    fs.nfa = -1;
    if (mpi_rank == 0) {
      disp = new int[mpi_size];
      disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        disp[rank] = disp[rank-1] + count[rank-1];
      fs.nfa = disp[mpi_size-1] + count[mpi_size-1];
    }

    // now gather the global face indices to rank 0...

    my_ifa_fwh_global = new int8[nfa_fwh];
    for (int ii = 0; ii < nfa_fwh; ++ii) {
      const int ifa = ifaVec[ii].first;
      my_ifa_fwh_global[ii] = getIfaGlobal(ifa);
    }
  }

  int8 * ifa_fwh_global = NULL;
  if (mpi_rank == 0) ifa_fwh_global = new int8[fs.nfa];

  MPI_Gatherv(my_ifa_fwh_global,nfa_fwh,MPI_INT8,ifa_fwh_global,count,disp,MPI_INT8,0,mpi_comm);
  delete[] my_ifa_fwh_global;

  if (mpi_rank == 0) {

    // we need to sort these faces so they are always written in the
    // same order...

    vector<pair<int8,int> > sortVec(fs.nfa);
    for (int ifa_fwh = 0; ifa_fwh < fs.nfa; ++ifa_fwh) {
      sortVec[ifa_fwh].first = ifa_fwh_global[ifa_fwh];
      sortVec[ifa_fwh].second = ifa_fwh;
    }

    sort(sortVec.begin(),sortVec.end());

    // put the sorted order into ifa_fwh_global...

    for (int ifa_fwh = 0; ifa_fwh < fs.nfa; ++ifa_fwh) {
      ifa_fwh_global[sortVec[ifa_fwh].second] = ifa_fwh;
    }

  }

  // now gather the face coordinate and weighted face normal data from all ranks
  // to write the geom file...

  double my_buf[4] = { 0.0, 0.0, 0.0, 0.0 };
  double * my_buf_double = new double[nfa_fwh*6]; // x,y,z,nx,ny,nz
  if (!zone_indices.empty()) {
    int ind = 0;
    for (int ii = 0, lim = zone_indices.size(); ii < lim; ++ii) {
      const int izone = zone_indices[ii];
      for (int ibf = 0; ibf < bfZoneVec[izone].nbf; ++ibf) {
        FOR_I3 my_buf_double[ind*6+i] = bfZoneVec[izone].x_bf[ibf][i];
        FOR_I3 my_buf_double[ind*6+3+i] = -bfZoneVec[izone].n_bf[ibf][i]; // flip normal to be outward pointing
        ++ind;
        // also sum for checking...
        FOR_I3 my_buf[i] += bfZoneVec[izone].n_bf[ibf][i];
        my_buf[3] += MAG(bfZoneVec[izone].n_bf[ibf]);
      }
    }
    assert(ind == nfa_fwh);
  }
  else {
    for (int ii = 0; ii < nfa_fwh; ++ii) {
      const int ifa = ifaVec[ii].first;
      FOR_I3 my_buf_double[ii*6+i] = x_fa[ifa][i];
      const double wgt = ifaVec[ii].second;
      double this_n[3]; FOR_I3 this_n[i] = wgt*n_fa[ifa][i];
      FOR_I3 my_buf_double[ii*6+3+i] = this_n[i];
      // also sum for checking...
      const double area = MAG(this_n);
      FOR_I3 my_buf[i] += this_n[i];
      my_buf[3] += area;
    }
  }

  // report counts/sums...

  double buf[4];
  MPI_Reduce(my_buf,buf,4,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  if (mpi_rank == 0)
    cout << " > fwh face count: " << fs.nfa << " normal sum: " << COUT_VEC(buf) << " weighted area: " << buf[3] << endl;

  // gather on rank0...

  double * buf_double = NULL;
  if (mpi_rank == 0) {
    buf_double = new double[fs.nfa*6];
    FOR_RANK {
      count[rank] *= 6;
      disp[rank] *= 6;
    }
  }
  MPI_Gatherv(my_buf_double,nfa_fwh*6,MPI_DOUBLE,buf_double,count,disp,MPI_DOUBLE,0,mpi_comm);
  delete[] my_buf_double;

  if (mpi_rank == 0) {

    char filename[256];
    sprintf(filename,"%s.geom",fs.name.c_str());
    string tmp_filename = MiscUtils::makeTmpPrefix(filename);
    MiscUtils::mkdir_for_file(filename);
    FILE * fp = fopen(tmp_filename.c_str(),"wb");

    int ibuf[4];
    ibuf[0] = UGP_IO_MAGIC_NUMBER;
    ibuf[1] = UGP_IO_VERSION;
    ibuf[2] = FFW_DATA_D2D2;
    ibuf[3] = fs.nfa;
    fwrite(ibuf, sizeof(int),4,fp);

    // we need a double[3] array large enough to hold the sorted x's, then
    // the sorted n's...

    double (*buf_sort)[3] = new double[fs.nfa][3];

    // pack x's first in sorted order...
    for (int ifa_fwh = 0; ifa_fwh < fs.nfa; ++ifa_fwh) {
      const int ifa_fwh_sort = ifa_fwh_global[ifa_fwh];
      FOR_I3 buf_sort[ifa_fwh_sort][i] = buf_double[ifa_fwh*6+i];
    }

    // write x's...
    fwrite(buf_sort, sizeof(double),fs.nfa*3,fp);

    // pack n's second in sorted order...
    for (int ifa_fwh = 0; ifa_fwh < fs.nfa; ++ifa_fwh) {
      const int ifa_fwh_sort = ifa_fwh_global[ifa_fwh];
      FOR_I3 buf_sort[ifa_fwh_sort][i] = buf_double[ifa_fwh*6+3+i];
    }
    delete[] ifa_fwh_global;

    delete[] buf_double;

    // write n's...
    fwrite(buf_sort, sizeof(double),fs.nfa*3,fp);
    delete[] buf_sort;

    // and close...
    fclose(fp);
    remove(filename);
    rename(tmp_filename.c_str(),filename);

  }

  // if we have been asked to write the surface, then we ALL participate...

  if (b_tecplot) {

    int * fa_flag = new int[nfa];
    FOR_IFA fa_flag[ifa] = 0;
    for (int ii = 0; ii < nfa_fwh; ++ii) {
      const int ifa = ifaVec[ii].first;
      fa_flag[ifa] = 1;
    }

    char filename[256];
    sprintf(filename,"%s.dat",fs.name.c_str());
    writeFlaggedFacesASCII(filename,fa_flag);
    delete[] fa_flag;

  }

  // build structures to write the data...

  vector<pair<int,pair<int8,double> > > dpwVec; // data-index/probe-index/wgt: icv, ifa_global, wgt

  if (!zone_indices.empty()) {
    int ind = 0;
    dpwVec.resize(nfa_fwh);
    for (int ii = 0, lim = zone_indices.size(); ii < lim; ++ii) {
      const int izone = zone_indices[ii];
      for (int ibf = 0; ibf < bfZoneVec[izone].nbf; ++ibf) {
        const int icv0 = bfZoneVec[izone].cvobf[ibf];
        dpwVec[ind++] = pair<int,pair<int8,double> >(icv0,pair<int8,double>(bfZoneVec[izone].ibf_global[ibf],1.0));
      }
    }
    assert(ind == nfa_fwh);
  }
  else {
    for (int ifa = 0; ifa < nfa; ++ifa) {
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      const int delta = cv_flag[icv0] - cv_flag[icv1];
      if (delta != 0) {
        // this face's value will be 0.5*(icv0+icv1)...
        const int8 ifa_global = getAbsIfaGlobal(ifa); // return the positive version here -- no tie breaking required
        dpwVec.push_back(pair<int,pair<int8,double> >(icv0,pair<int8,double>(ifa_global,0.5)));
        if (icv1 < ncv) {
          dpwVec.push_back(pair<int,pair<int8,double> >(icv1,pair<int8,double>(ifa_global,0.5)));
        }
      }
    }

    delete[] cv_flag;
  }

  // we now have a vector of all cv's we need to send. Sort this so we only pack and send each piece of data once...

  sort(dpwVec.begin(),dpwVec.end()); // cv's determine sort here

  // now count how many unique pieces of cv data are involved, and need to be
  // sent for each field requested...

  vector<int> dVec; // a compressed vector containing just the unique cv indices
  vector<pair<int,pair<int8,double> > >::iterator iter_end = dpwVec.begin();
  while (iter_end != dpwVec.end()) {
    vector<pair<int,pair<int8,double> > >::iterator iter_begin = iter_end;
    do {
      ++iter_end;
    }
    while ((iter_end != dpwVec.end())&&(iter_end->first == iter_begin->first));
    // store the unique cv numbers in the compressed cv vec, and replace the
    // icv int in the dpwVec with an incrementing index. This gets adjusted
    // below with an offset to account for process striping...
    const int index = dVec.size();
    dVec.push_back(iter_begin->first);
    for (vector<pair<int,pair<int8,double> > >::iterator iter = iter_begin; iter != iter_end; ++iter)
      iter->first = index;
  }

  // before we reduce, adjust the offset of the icv index...

  int8 cv_count = dVec.size();
  int8 cv_offset;
  MPI_Scan(&cv_count,&cv_offset,1,MPI_INT8,MPI_SUM,mpi_comm);
  assert(cv_offset < TWO_BILLION); // make sure global cvs are below 2B so we can use ints
  cv_offset -= cv_count;
  for (vector<pair<int,pair<int8,double> > >::iterator iter = dpwVec.begin(); iter != dpwVec.end(); ++iter)
    iter->first += cv_offset;

  // now reduce the full arrays (icv-as-index, ifa_global, wgt) to rank 0 (for now)...

  int my_ndpw = dpwVec.size();
  MPI_Gather(&my_ndpw,1,MPI_INT,count,1,MPI_INT,0,mpi_comm);

  fs.ndpw = -1;
  int8 * probe_index_global = NULL;
  if (mpi_rank == 0) {
    disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      disp[rank] = disp[rank-1] + count[rank-1];
    fs.ndpw = disp[mpi_size-1] + count[mpi_size-1];
    cout << " > fs.ndpw: " << fs.ndpw << endl;
    if (zone_indices.empty()) {
      // for fwh data that is based on the simple average of cvs, this
      // should be 2x the face count...
      assert(fs.ndpw == 2*fs.nfa);
    }
    else {
      assert(fs.ndpw == fs.nfa);
    }
    // for gathering at the write process (rank == 0 for now)...
    assert(fs.cv_unpack_of_dpw == NULL); fs.cv_unpack_of_dpw = new int[fs.ndpw];
    probe_index_global                                       = new int8[fs.ndpw];
    assert(fs.wgt_of_dpw == NULL); fs.wgt_of_dpw             = new double[fs.ndpw];
  }

  // repack local sorted dpwVec to send to write process...
  int * data_index   = new int[my_ndpw];
  int8 * probe_index = new int8[my_ndpw];
  double * wgt       = new double[my_ndpw];

  for (int ii = 0; ii < my_ndpw; ++ii) {
    data_index[ii]  = dpwVec[ii].first;
    probe_index[ii] = dpwVec[ii].second.first;
    wgt[ii]         = dpwVec[ii].second.second;
  }
  dpwVec.clear();

  MPI_Gatherv(data_index,my_ndpw,MPI_INT,fs.cv_unpack_of_dpw,count,disp,MPI_INT,0,mpi_comm);
  MPI_Gatherv(probe_index,my_ndpw,MPI_INT8,probe_index_global,count,disp,MPI_INT8,0,mpi_comm);
  MPI_Gatherv(wgt,my_ndpw,MPI_DOUBLE,fs.wgt_of_dpw,count,disp,MPI_DOUBLE,0,mpi_comm);

  delete[] data_index;
  delete[] probe_index;
  delete[] wgt;

  if (mpi_rank == 0) {

    vector<pair<int8,int> > sortVec(fs.ndpw);
    for (int ii = 0; ii < fs.ndpw; ++ii) {
      sortVec[ii].first = probe_index_global[ii];
      sortVec[ii].second = ii;
    }
    delete[] probe_index_global;

    sort(sortVec.begin(),sortVec.end());

    // because we have 2 cvs associated with every data point,
    // this should sort into each global face represented twice...

    assert(fs.fa_fwh_of_dpw == NULL); fs.fa_fwh_of_dpw = new int[fs.ndpw];

    int nfa_check = 0;
    int8 ifa_global_current = -1;
    for (int ii = 0; ii < fs.ndpw; ++ii) {
      if (sortVec[ii].first != ifa_global_current) {
        // this is a new (or the first) ifa_fwh...
        ++nfa_check;
        ifa_global_current = sortVec[ii].first;
      }
      // every time, we do this...
      fs.fa_fwh_of_dpw[sortVec[ii].second] = nfa_check-1; // 0-indexed fwh face
    }
    assert(nfa_check == fs.nfa);

  }

  // to gather the cv data, everyone needs to know who to pack and in what order. This
  // is currently in dVec...

  fs.ncv_pack = dVec.size();
  assert(fs.cv_pack_array == NULL); fs.cv_pack_array = new int[fs.ncv_pack];
  for (int ii = 0; ii < dVec.size(); ++ii) {
    fs.cv_pack_array[ii] = dVec[ii];
    assert((fs.cv_pack_array[ii] >= 0)&&(fs.cv_pack_array[ii] < ncv));
  }
  dVec.clear();

  // on the unpack side (the writer) we need the count and disp for the gatherv routine...
  // We put this in local count/disp for now...

  MPI_Gather(&fs.ncv_pack,1,MPI_INT,count,1,MPI_INT,0,mpi_comm);
  if (mpi_rank == 0) {
    disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      disp[rank] = disp[rank-1] + count[rank-1];
    fs.ncv_unpack = disp[mpi_size-1] + count[mpi_size-1];
    cout << " > fs.ncv_unpack: " << fs.ncv_unpack << endl;
  }

  // just take count,disp into the fs and
  // return without delete[]...
  assert(fs.count == NULL); fs.count = count;
  assert(fs.disp == NULL);  fs.disp = disp;

  // for the case of a killfile fwh request, everything will go through
  // rank 0, so no need to do these bcast's. Otherwise, we setup one rank
  // on each node with the necessary buffers to handle the gather and
  // pack the data properly...

  if (!killfile) {

    // and populate the other nodes that may be called
    // on to write...

    if (mpi_rank_shared == 0) {

      // these are the nodes that can have write status...
      // make sure of this...
      int rank_internode;
      for (rank_internode = 0; rank_internode < mpi_size_internode; ++rank_internode)
        if (mpi_rank == rank_of_rank_internode[rank_internode])
          break;
      // ensure we exited by the break...
      assert(rank_internode < mpi_size_internode);

      int buf[3] = { fs.nfa, fs.ncv_unpack, fs.ndpw };
      MPI_Bcast(buf,3,MPI_INT,0,mpi_comm_internode);

      if (mpi_rank != 0) {
        fs.nfa        = buf[0];
        fs.ncv_unpack = buf[1];
        fs.ndpw       = buf[2];
        assert(fs.count == NULL); fs.count = new int[mpi_size];
        // fs.disp done below...
        assert(fs.cv_unpack_of_dpw == NULL); fs.cv_unpack_of_dpw = new int[fs.ndpw];
        assert(fs.fa_fwh_of_dpw == NULL);    fs.fa_fwh_of_dpw = new int[fs.ndpw];
        assert(fs.wgt_of_dpw == NULL);       fs.wgt_of_dpw = new double[fs.ndpw];
      }

      MPI_Bcast(fs.count,           mpi_size, MPI_INT,    0, mpi_comm_internode);
      MPI_Bcast(fs.cv_unpack_of_dpw,fs.ndpw,  MPI_INT,    0, mpi_comm_internode);
      MPI_Bcast(fs.fa_fwh_of_dpw,   fs.ndpw,  MPI_INT,    0, mpi_comm_internode);
      MPI_Bcast(fs.wgt_of_dpw,      fs.ndpw,  MPI_DOUBLE, 0, mpi_comm_internode);

      if (mpi_rank != 0) {
        assert(fs.disp == NULL);
        fs.disp = new int[mpi_size];
        fs.disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          fs.disp[rank] = fs.disp[rank-1] + fs.count[rank-1];
        assert(fs.ncv_unpack == fs.disp[mpi_size-1] + fs.count[mpi_size-1]);
      }

    }

  }

  return 0;

}
