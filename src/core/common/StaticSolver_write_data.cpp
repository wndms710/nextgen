#include "StaticSolver.hpp"
#include "../../stitch/CuttableVoronoiData.hpp"

void StaticSolver::buildDelaunayDual() {

  // perturbation with global affine transform...
  double I_plus_eps_R[3][3];
  const double dual_eps = getDoubleParam("DUAL_EPS",1.0E-3);
  if (mpi_rank == 0) {
    FOR_I3 {
      FOR_J3 I_plus_eps_R[i][j] = dual_eps*(double(rand())/double(RAND_MAX)-0.5);
      I_plus_eps_R[i][i] += 1.0;
    }
  }
  MPI_Bcast((double*)I_plus_eps_R,9,MPI_DOUBLE,0,mpi_comm);

  double (*x_vv_tmp)[3] = new double[ncv][3];
  FOR_ICV FOR_I3 x_vv_tmp[icv][i] = DOT_PRODUCT(I_plus_eps_R[i],x_vv[icv]);
  double *delta_cv = new double[ncv];
  FOR_ICV delta_cv[icv] = 2.0*r_vv[icv] + DIST(x_vv_tmp[icv],x_vv[icv]); // I think this is guaranteed to be large enough

  // to start the "even faster" rebuild, we need nbocv_i/v based on the
  // current delta_cv...
  int *nbocv_i = NULL;
  int *nbocv_v = NULL;
  map<uint8,int> rbiMap;

  // get all periodicity combinations...
  vector<int> bitVec;
  PeriodicData::buildBitVec(bitVec);

  // ALL cvs get adt'd...
  Adt<double> * x_vv_adt = new Adt<double>(ncv,x_vv_tmp,x_vv_tmp);
  double my_bbmin[3],my_bbmax[3];
  x_vv_adt->getBbox(my_bbmin,my_bbmax);
  double (*bbmin)[3] = new double[mpi_size][3];
  double (*bbmax)[3] = new double[mpi_size][3];
  MPI_Allgather(my_bbmin,3,MPI_DOUBLE,bbmin,3,MPI_DOUBLE,mpi_comm);
  MPI_Allgather(my_bbmax,3,MPI_DOUBLE,bbmax,3,MPI_DOUBLE,mpi_comm);
  Adt<double> * x_vv_bbox_adt = new Adt<double>(mpi_size,bbmin,bbmax);
  delete[] bbmin;
  delete[] bbmax;

  // now collect nbrs...
  int * send_count = new int[mpi_size];
  FOR_RANK send_count[rank] = 0;

  vector<int> intVec;
  FOR_ICV {
    if (bfocv_i[icv] == bfocv_i[icv+1]) {
      for (int ibv = 0; ibv < bitVec.size(); ++ibv) {
        assert(intVec.empty());
        const int bits = bitVec[ibv];
        double xp_t[3]; FOR_I3 xp_t[i] = x_vv_tmp[icv][i];
        if (bits != 0) PeriodicData::periodicTranslate(xp_t,1,bits);
        // get the ranks where we have to send the (periodically transformed) points...
        x_vv_bbox_adt->buildListForSphere(intVec,xp_t,delta_cv[icv]);
        for (int ii = 0; ii < intVec.size(); ++ii) {
          const int rank = intVec[ii]; assert((rank >= 0)&&(rank < mpi_size));
          if (bits||(rank != mpi_rank)) {
            ++send_count[rank];
          }
        }
        intVec.clear();
      }
    }
  }
  int * send_disp = new int[mpi_size];
  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
  int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

  // send 2 ints and 4 doubles...
  int * send_buf_int = new int[send_count_sum*2];
  double * send_buf_double = new double[send_count_sum*4];

  // now pack...
  FOR_ICV {
    if (bfocv_i[icv] == bfocv_i[icv+1]) {
      for (int ibv = 0; ibv < bitVec.size(); ++ibv) {
        // get the ranks where we have to send the (periodically transformed) points...
        assert(intVec.empty());
        const int bits = bitVec[ibv];
        double xp_t[3]; FOR_I3 xp_t[i] = x_vv_tmp[icv][i];
        if (bits != 0) PeriodicData::periodicTranslate(xp_t,1,bits);
        x_vv_bbox_adt->buildListForSphere(intVec,xp_t,delta_cv[icv]);
        const int inv_bits = BitUtils::flipPeriodicBits(bits);
        for (int ii = 0; ii < intVec.size(); ++ii) {
          const int rank = intVec[ii]; assert((rank >= 0)&&(rank < mpi_size));
          if (bits||(rank != mpi_rank)) {
            send_buf_int[send_disp[rank]*2] = BitUtils::packRankBits(mpi_rank,inv_bits);
            send_buf_int[send_disp[rank]*2+1] = icv;
            FOR_I3 send_buf_double[send_disp[rank]*4+i] = xp_t[i];
            send_buf_double[send_disp[rank]*4+3] = delta_cv[icv];
            ++send_disp[rank];
          }
        }
        intVec.clear();
      }
    }
  }

  delete x_vv_bbox_adt;

  // rewind...
  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

  // recv-side stuff...
  int * recv_count = new int[mpi_size];
  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

  int * recv_disp = new int[mpi_size];
  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
  int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

  int * recv_buf_int = new int[recv_count_sum*2];
  FOR_RANK {
    send_count[rank] *= 2;
    send_disp[rank] *= 2;
    recv_count[rank] *= 2;
    recv_disp[rank] *= 2;
  }
  MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
      recv_buf_int,recv_count,recv_disp,MPI_INT,
      mpi_comm);
  delete[] send_buf_int; send_buf_int = NULL;

  double * recv_buf_double = new double[recv_count_sum*4];
  FOR_RANK {
    send_count[rank] *= 2;
    send_disp[rank] *= 2;
    recv_count[rank] *= 2;
    recv_disp[rank] *= 2;
  }
  MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
      recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,
      mpi_comm);
  delete[] send_buf_double; send_buf_double = NULL;

  FOR_RANK send_count[rank] = 0;
  for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
    assert(intVec.empty());
    x_vv_adt->buildListForSphere(intVec,recv_buf_double+4*irecv,recv_buf_double[4*irecv+3]);
    if (!intVec.empty()) {
      int rank,bits; BitUtils::unpackRankBits(rank,bits,recv_buf_int[irecv*2]); assert((rank >= 0)&&(rank < mpi_size));
      send_count[rank] += 2+2*intVec.size(); // send back icv, n, then (rank/bits,index) for n nbrs
      intVec.clear();
    }
  }

  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
  send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

  assert(send_buf_int == NULL);
  send_buf_int = new int[send_count_sum];
  set<int> rankSet; // store ranks we need to send to for prcomm build
  for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
    assert(intVec.empty());
    x_vv_adt->buildListForSphere(intVec,recv_buf_double+4*irecv,recv_buf_double[4*irecv+3]);
    if (!intVec.empty()) {
      int rank,bits; BitUtils::unpackRankBits(rank,bits,recv_buf_int[irecv*2]); assert((rank >= 0)&&(rank < mpi_size));
      rankSet.insert(rank);
      const int icv = recv_buf_int[irecv*2+1]; // icv on calling process
      send_buf_int[send_disp[rank]++] = icv; // send back icv, n, then (rank/bits,index) for n nbrs
      send_buf_int[send_disp[rank]++] = intVec.size();
      for (int ii = 0; ii < intVec.size(); ++ii) {
        const int icv_nbr = intVec[ii]; assert((icv_nbr >= 0)&&(icv_nbr < ncv));
        send_buf_int[send_disp[rank]++] = BitUtils::packRankBits(mpi_rank,bits);
        if (bfocv_i[icv_nbr+1] == bfocv_i[icv_nbr])
          send_buf_int[send_disp[rank]++] = icv_nbr;
        else
          send_buf_int[send_disp[rank]++] = -icv_nbr-1; // -1 idx to indicate boundary
      }
      intVec.clear();
    }
  }
  delete[] recv_buf_int; recv_buf_int = NULL;
  delete[] recv_buf_double; recv_buf_double = NULL;

  // rewind...
  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
  recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

  assert(recv_buf_int == NULL);
  recv_buf_int = new int[recv_count_sum];
  MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
      recv_buf_int,recv_count,recv_disp,MPI_INT,
      mpi_comm);
  delete[] send_buf_int; send_buf_int = NULL;
  delete[] send_count;
  delete[] send_disp;
  delete[] recv_count;
  delete[] recv_disp;

  // count the local nbrs...

  assert(nbocv_i == NULL);
  nbocv_i = new int[ncv+1];
  FOR_ICV {
    if (bfocv_i[icv] == bfocv_i[icv+1]) {
      assert(intVec.empty());
      x_vv_adt->buildListForSphere(intVec,x_vv_tmp[icv],delta_cv[icv]);
      nbocv_i[icv+1] = intVec.size()-1; // should include self
      intVec.clear();
    }
    else {
      nbocv_i[icv+1] = 0;
    }
  }

  int irecv = 0;
  while (irecv < recv_count_sum) {
    const int icv = recv_buf_int[irecv++]; assert((icv >= 0)&&(icv < ncv));
    const int n = recv_buf_int[irecv++];
    nbocv_i[icv+1] += n;
    irecv += n*2;
  }
  assert(irecv == recv_count_sum);

  // now record nbrs...

  nbocv_i[0] = 0;
  FOR_ICV nbocv_i[icv+1] += nbocv_i[icv];
  assert(nbocv_v == NULL);
  nbocv_v = new int[nbocv_i[ncv]];

  FOR_ICV {
    if (bfocv_i[icv] == bfocv_i[icv+1]) {
      assert(intVec.empty());
      x_vv_adt->buildListForSphere(intVec,x_vv_tmp[icv],delta_cv[icv]);
      bool b_found = false;
      for (int ii = 0; ii < intVec.size(); ++ii) {
        const int icv_nbr = intVec[ii];
        if (icv_nbr == icv) {
          assert(!b_found);
          b_found = true;
        }
        else {
          nbocv_v[nbocv_i[icv]++] = icv_nbr;
        }
      }
      assert(b_found);
      intVec.clear();
    }
  }

  delete x_vv_adt;

  ncv_d = ncv; // start at ghost
  irecv = 0;
  while (irecv < recv_count_sum) {
    const int icv = recv_buf_int[irecv++]; assert((icv >= 0)&&(icv < ncv));
    const int n = recv_buf_int[irecv++];
    for (int ii = 0; ii < n; ++ii) {
      int rank,bits; BitUtils::unpackRankBits(rank,bits,recv_buf_int[irecv++]); assert((rank >= 0)&&(rank < mpi_size));
      assert(bits||(rank != mpi_rank));
      const int index = recv_buf_int[irecv++];
      uint8 rbi;
      if (index < 0) 
        rbi = BitUtils::packRankBitsIndex(rank,bits,-index-1);
      else
        rbi = BitUtils::packRankBitsIndex(rank,bits,index);
      // get an index for this guy...
      map<uint8,int>::const_iterator it = rbiMap.find(rbi);
      if (it == rbiMap.end()) {
        if (index < 0)
          rbiMap[rbi] = -ncv_d-1; // -1 idx indicates boundary
        else
          rbiMap[rbi] = ncv_d;
        nbocv_v[nbocv_i[icv]++] = ncv_d++;
      }
      else {
        if (it->second < 0)
          nbocv_v[nbocv_i[icv]++] = -it->second-1;
        else
          nbocv_v[nbocv_i[icv]++] = it->second;
      }
    }
  }
  assert(irecv == recv_count_sum);
  delete[] recv_buf_int; recv_buf_int = NULL;

  // rebuild csr...
  for (int icv = ncv-1; icv > 0; --icv)
    nbocv_i[icv] = nbocv_i[icv-1];
  nbocv_i[0] = 0;

  assert(rbiMap.size() == ncv_d-ncv);
  int *cvocv0 = new int[ncv_d-ncv];
  rbi_d = new uint8[ncv_d-ncv];
  bool *b_boundary = new bool[ncv_d-ncv];
  int ncv_d_check = ncv;
  for (map<uint8,int>::const_iterator it = rbiMap.begin(); it != rbiMap.end(); ++it) {
    int rank,bits,index;
    BitUtils::unpackRankBitsIndex(rank,bits,index,it->first);
    assert((bits)||(rank != mpi_rank));
    rbi_d[ncv_d_check-ncv] = it->first;
    if (it->second < 0) {
      b_boundary[ncv_d_check-ncv] = true;
      cvocv0[(-it->second-1)-ncv] = ncv_d_check++;
      const int icv = -it->second-1;
      assert((icv >= ncv)&&(icv < ncv_d));
    }
    else {
      b_boundary[ncv_d_check-ncv] = false;
      cvocv0[it->second-ncv] = ncv_d_check++;
      assert((it->second >= ncv)&&(it->second < ncv_d));
    }
  }
  assert(ncv_d_check == ncv_d);
  rbiMap.clear();

  FOR_ICV {
    for (int noc = nbocv_i[icv]; noc < nbocv_i[icv+1]; ++noc) {
      const int icv_nbr0 = nbocv_v[noc];
      assert((icv_nbr0 >= 0)&&(icv_nbr0 < ncv_d));
      assert(icv_nbr0 != icv);
      if (icv_nbr0 >= ncv)
        nbocv_v[noc] = cvocv0[icv_nbr0-ncv];
    }
  }
  delete[] cvocv0;

  buildPrcomm(cvdPrcommVec,rbi_d,rankSet,ncv,ncv_d);

  // get the ghost x_vv_tmp...
  double (*x_vv_tmp_g)[3] = new double[ncv_d-ncv][3];
  updateCvdDataSeparateGhosts(x_vv_tmp,x_vv_tmp_g,REPLACE_TRANSLATE_DATA);

  CuttableVoronoiData cvd;
  vector<pair<double,int> > nbrVec;
  FOR_ICV {
    if (bfocv_i[icv+1] == bfocv_i[icv]) {
      assert(nbrVec.empty());

      for (int noc = nbocv_i[icv]; noc != nbocv_i[icv+1]; ++noc) {
        const int icv_nbr = nbocv_v[noc];
        assert((icv_nbr >= 0)&&(icv_nbr < ncv_d)&&(icv_nbr != icv));
        if (icv_nbr < ncv)
          nbrVec.push_back(pair<double,int>(DIST2(x_vv_tmp[icv],x_vv_tmp[icv_nbr]),icv_nbr));
        else
          nbrVec.push_back(pair<double,int>(DIST2(x_vv_tmp[icv],x_vv_tmp_g[icv_nbr-ncv]),icv_nbr));
      }

      // sort nbrs...
      sort(nbrVec.begin(),nbrVec.end());

      assert(cvd.nno == 0);
      cvd.addCube(0.505*delta_cv[icv]);
      cvd.setD2Max();

      // and cut...
      // we now have a sorted list of nbrs. Cut the cvd against the nbrs until it cannot possibly
      // be cut anymore...
      for (int in = 0; in < nbrVec.size(); ++in) {
        // recall that the .first contains the d2 nbr distance. and cvd.d2_max contains
        // the maximum node radius. as soon as d2 >= 4*cvd.d2_max, there cannot be any
        // more nbrs that cut this cv...
        if (nbrVec[in].first > 4.0*cvd.d2_max)
          break; // done: all remaining nbrs are too far away to possibly cut the current vd...
        const int icv_nbr = nbrVec[in].second;
        double dn[3];
        if (icv_nbr < ncv) {
          FOR_I3 dn[i] = 0.5*(x_vv_tmp[icv_nbr][i]-x_vv_tmp[icv][i]);
        }
        else {
          FOR_I3 dn[i] = 0.5*(x_vv_tmp_g[icv_nbr-ncv][i]-x_vv_tmp[icv][i]);
        }
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
        cvd.cut(dn,-8-nbrVec[in].second); // note that nbrs use a -8 convention. cvd.d2_max updated automatically
      }
      nbrVec.clear();

      // skip vd's that are not completely internal...
      if (cvd.hasSeedBoundary()) {
        cvd.clear();
        continue;
      }

      try {
        cvd.check();
      }
      catch(...) {
        // just skips problem cells, so that the io doesn't bring down house, report errors though...
        char filename[128];
        sprintf(filename,"cvd_debug.%06d.bin",mpi_rank);
        cvd.writeBinary(filename);
        double zero[3] = { 0.0, 0.0, 0.0 };
        cvd.writeTecplot(mpi_rank,zero);
        cout << "Error: failed check: rank: " << mpi_rank << ", icv: " << icv << ", x_vv: " << COUT_VEC(x_vv[icv]) << endl;
        //throw(0);
        cvd.clear();
        continue;
      }

      // extract tets...
      int (*cvono)[3] = new int[cvd.nno][3];
      for (int ino = 0; ino < cvd.nno; ++ino)
        FOR_I3 cvono[ino][i] = -1;
      for (int ied = 0; ied < cvd.ned; ++ied) {
        FOR_I2 {
          const int ino = cvd.nooed[ied][i];
          FOR_J2 {
            const int icv_nbr = -cvd.faoed[ied][j]-8;
            FOR_L3 {
              if (cvono[ino][l] == -1) {
                cvono[ino][l] = icv_nbr;
                break;
              }
            }
          }
        }
      }
      const uint8 rbi = BitUtils::packRankBitsIndex(mpi_rank,0,icv);
      double xtet[4][3]; FOR_I3 xtet[0][i] = x_vv_tmp[icv][i];
      for (int ino = 0; ino < cvd.nno; ++ino) {
        // see if icv owns tet (min rbi)...
        bool b_add = true;
        FOR_I3 {
          const int icv_nbr = cvono[ino][i];
          assert((icv_nbr >= 0)&&(icv_nbr < ncv_d));
          if (icv_nbr < ncv) {
            FOR_J3 xtet[1+i][j] = x_vv_tmp[icv_nbr][j];
            if ((icv_nbr < icv)&&(bfocv_i[icv_nbr+1]==bfocv_i[icv_nbr])) {
              b_add = false;
              break;
            }
          }
          else if (icv_nbr >= ncv) {
            FOR_J3 xtet[1+i][j] = x_vv_tmp_g[icv_nbr-ncv][j];
            if ((rbi_d[icv_nbr-ncv] < rbi)&&(!b_boundary[icv_nbr-ncv])) {
              b_add = false;
              break;
            }
          }
        }
        if (b_add) {
          // add dual tet such that it has positive signed volume...
          if (SIGNED_TET_VOLUME_6(xtet[0],xtet[1],xtet[2],xtet[3]) > 0.0) {
            dualTetVec.push_back(DualTet(icv,cvono[ino][0],cvono[ino][1],cvono[ino][2]));
          }
          else {
            assert(SIGNED_TET_VOLUME_6(xtet[0],xtet[1],xtet[3],xtet[2]) > 0.0);
            dualTetVec.push_back(DualTet(icv,cvono[ino][0],cvono[ino][2],cvono[ino][1]));
          }
        }

      }

      delete[] cvono;
      cvd.clear();
    }
    else {
      assert(nbocv_i[icv+1] == nbocv_i[icv]);
    }
  }
  delete[] b_boundary;
  delete[] nbocv_i;
  delete[] nbocv_v;
  delete[] x_vv_tmp_g;
  delete[] x_vv_tmp;
  delete[] delta_cv;

  assert(x_cv_d == NULL);
  x_cv_d = new double[ncv_d-ncv][3];
  updateCvdDataSeparateGhosts(x_cv,x_cv_d,REPLACE_TRANSLATE_DATA);

  //dumpRange(x_cv,ncv,"x_cv");
  //dumpRange(x_cv_d,ncv_d-ncv,"x_cv_d");

  uint8 my_ntet = dualTetVec.size();
  uint8 ntet;
  MPI_Allreduce(&my_ntet,&ntet,1,MPI_UINT8,MPI_SUM,mpi_comm);
  if (mpi_rank == 0)
    cout << " > created " << ntet << " dual tetrahedra." << endl;

}

// convenience macro to set hanging node cv-face index
#define IHN_VAL (cvofa[ifa][0]==icv) ? 0 : 1;

// main routine...
void StaticSolver::writeData(DataWriter* dw,const int step,const double time) {

  if (dw->b_dual&&(x_cv_d == NULL)) {
    buildDelaunayDual();
    subSurface->buildPrcomms();
  }

  if (dw->topo == ISO_TOPO) {
    if (dw->format == FLUENT_BP_FORMAT) {
      WUI(WARN," ! the FLUENT_PROFILE format is not supported for this GEOM; only FAZONE or PLANE are allowed");
      return;
    }

    assert((dw->nscalars_bf == 0)&&(dw->nvectors_bf == 0));
    assert((dw->nscalars_lp == 0)&&(dw->nvectors_lp == 0));
    double *dn_no = NULL;
    double (*dn3_no)[3] = NULL;
    double * dn_tmp = NULL; // allocated in buildIsoWithWgtsAndData
    double (*dn3_tmp)[3] = NULL; // allocated in buildIsoWithWgtsAndData
    vector<SimpleTriWithWgts> triVec;

    if (dw->b_dual) {

      if (dw->geom == USERDEF_GEOM) {
        if (dw->iso_var == NULL) dw->iso_var = new double[ncv_d];
        // have to rebuild every time in case the iso definition is time dependent
        CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->iso_var_name);
        assert((var != NULL)&&(var->getType() == DN_DATA)&&(var->getUnindexedTopology() == CV_DATA));
        double* ptr = var->getDNptr();
        FOR_ICV dw->iso_var[icv] = ptr[icv];
      }
      else if (dw->geom == FILE_GEOM) {
        if (dw->iso_var == NULL) {
          dw->iso_var = new double[ncv_d];
          SurfaceShm surface;
          surface.readBinary(dw->sbin_filename);
          surface.initPointIsInside(x_cv,ncv);
          if (cvAdt == NULL) buildCvAdt();
          FOR_ICV dw->iso_var[icv] = HUGE_VAL;
          vector<int> intVec;
          for (int ist = 0; ist < surface.nst; ++ist) {
            const double x_st[3] = {(surface.xsp[surface.spost[ist][0]][0]+
                surface.xsp[surface.spost[ist][1]][0]+
                surface.xsp[surface.spost[ist][2]][0])/3.0,
                  (surface.xsp[surface.spost[ist][0]][1]+
                   surface.xsp[surface.spost[ist][1]][1]+
                   surface.xsp[surface.spost[ist][2]][1])/3.0,
                  (surface.xsp[surface.spost[ist][0]][2]+
                   surface.xsp[surface.spost[ist][1]][2]+
                   surface.xsp[surface.spost[ist][2]][2])/3.0};

            assert(intVec.empty());
            cvAdt->buildListForPoint(intVec,x_st);
            for (int ii = 0, ii_end=intVec.size(); ii < ii_end; ++ii) {
              const int icv = intVec[ii]; assert((icv >= 0)&&(icv < ncv));
              dw->iso_var[icv] = min(dw->iso_var[icv],DIST(x_st,x_cv[icv]));
            }
            intVec.clear();
          }
          FOR_ICV {
            if (surface.pointIsInside(x_cv[icv]))
              dw->iso_var[icv] *= -1.0;
          }
          dw->iso_var_value = 0.0;
        }
      }
      else {
        assert(dw->geom == SIMPLE_GEOM);
        // only need to build once
        if (dw->iso_var == NULL) {
          dw->iso_var_value = 0.0;
          dw->iso_var = new double[ncv_d];
          FOR_ICV dw->iso_var[icv] = dw->simple_geom->pointSignedDistance(x_cv[icv]);
        }
      }
      updateCvdData(dw->iso_var);

      if (dw->nscalars_cv > 0) {
        dn_no = new double[ncv_d*dw->nscalars_cv];
        for (int ivar = 0; ivar < dw->nscalars_cv; ++ivar) {
          CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->scalar_names_cv[ivar]);
          assert((var != NULL)&&(var->getType() == DN_DATA)&&(var->getUnindexedTopology() == CV_DATA));
          double* ptr = var->getDNptr();
          FOR_ICV dn_no[icv+ncv_d*ivar] = ptr[icv];
          updateCvdData(dn_no+ncv_d*ivar);
        }
      }
      if (dw->nvectors_cv > 0) {
        dn3_no = new double[ncv_d*dw->nvectors_cv][3]; // reallocated/shuffled in buildIsoWithWgtsAndData
        for (int ivar = 0; ivar < dw->nvectors_cv; ++ivar) {
          CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->vector_names_cv[ivar]);
          assert((var != NULL)&&(var->getType() == DN3_DATA)&&(var->getUnindexedTopology() == CV_DATA));
          double (*ptr)[3] = var->getDN3ptr();
          FOR_ICV FOR_I3 dn3_no[icv+ncv_d*ivar][i] = ptr[icv][i];
          updateCvdData(dn3_no+ncv_d*ivar,REPLACE_ROTATE_DATA);
        }
      }

      // shuffle to group collocated data
      if (dw->nscalars_cv > 0) {
        assert(dn_no != NULL);
        double *_dn_no = dn_no;
        dn_no = new double[ncv_d*dw->nscalars_cv];
        for (int icv = 0; icv < ncv_d; ++icv)
          for (int ivar = 0; ivar < dw->nscalars_cv; ++ivar)
            dn_no[icv*dw->nscalars_cv+ivar] = _dn_no[ncv_d*ivar+icv];
        delete[] _dn_no;
      }
      if (dw->nvectors_cv > 0) {
        assert(dn3_no != NULL);
        double (*_dn3_no)[3] = dn3_no;
        dn3_no = new double[ncv_d*dw->nvectors_cv][3];
        for (int icv = 0; icv < ncv_d; ++icv)
          for (int ivar = 0; ivar < dw->nvectors_cv; ++ivar)
            FOR_I3 dn3_no[icv*dw->nvectors_cv+ivar][i] = _dn3_no[ncv_d*ivar+icv][i];
        delete[] _dn3_no;
      }

      const int ndt = dualTetVec.size();
      for (int idt = 0; idt < ndt; ++idt) {
        double xtet[4][3];
        FOR_I4 {
          const int icv = dualTetVec[idt].icv[i];
          if (icv < ncv)
            FOR_J3 xtet[i][j] = x_cv[icv][j];
          else
            FOR_J3 xtet[i][j] = x_cv_d[icv-ncv][j];
        }
        tetCutPlaneWgts(triVec,
            dw->iso_var[dualTetVec[idt].icv[0]]-dw->iso_var_value,dw->iso_var[dualTetVec[idt].icv[1]]-dw->iso_var_value,
            dw->iso_var[dualTetVec[idt].icv[2]]-dw->iso_var_value,dw->iso_var[dualTetVec[idt].icv[3]]-dw->iso_var_value,
            xtet[0],xtet[1],xtet[2],xtet[3],
            dualTetVec[idt].icv[0],dualTetVec[idt].icv[1],dualTetVec[idt].icv[2],dualTetVec[idt].icv[3]);
      }

    }
    else {

      if (dw->geom == USERDEF_GEOM) {
        if (dw->iso_var == NULL) dw->iso_var = new double[nno];
        // have to rebuild every time in case the iso definition is time dependent
        setNoDN(dw->iso_var,dw->iso_var_name); // constructs iso from string evaluation
      }
      else if (dw->geom == FILE_GEOM) {
        if (dw->iso_var == NULL) {
          dw->iso_var = new double[nno];
          double *sgn_no = new double[nno];
          SurfaceShm surface;
          surface.readBinary(dw->sbin_filename);
          surface.initPointIsInside(x_no,nno);
          if (cvAdt == NULL) buildCvAdt();
          FOR_INO dw->iso_var[ino] = HUGE_VAL;
          FOR_INO {
            if (surface.pointIsInside(x_no[ino]))
              sgn_no[ino] = -1.0;
            else
              sgn_no[ino] = +1.0;
          }
          vector<int> intVec;
          for (int ist = 0; ist < surface.nst; ++ist) {
            const double x_st[3] = {(surface.xsp[surface.spost[ist][0]][0]+
                surface.xsp[surface.spost[ist][1]][0]+
                surface.xsp[surface.spost[ist][2]][0])/3.0,
                  (surface.xsp[surface.spost[ist][0]][1]+
                   surface.xsp[surface.spost[ist][1]][1]+
                   surface.xsp[surface.spost[ist][2]][1])/3.0,
                  (surface.xsp[surface.spost[ist][0]][2]+
                   surface.xsp[surface.spost[ist][1]][2]+
                   surface.xsp[surface.spost[ist][2]][2])/3.0};

            assert(intVec.empty());
            cvAdt->buildListForPoint(intVec,x_st);
            for (int ii = 0, ii_end=intVec.size(); ii < ii_end; ++ii) {
              const int icv = intVec[ii]; assert((icv >= 0)&&(icv < ncv));
              for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
                const int ifa = faocv_v[foc];
                for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
                  const int ino = noofa_v[nof]; assert((ino >= 0)&&(ino < nno));
                  dw->iso_var[ino] = min(dw->iso_var[ino],DIST(x_st,x_no[ino]));
                }
              }
            }
            intVec.clear();
          }
          updateNoData(dw->iso_var,MIN_NO_PERIODIC_DATA);
          FOR_INO dw->iso_var[ino] *= sgn_no[ino];
          dw->iso_var_value = 0.0;
          delete[] sgn_no;
        }
      }
      else {
        assert(dw->geom == SIMPLE_GEOM);
        // only need to build once
        if (dw->iso_var == NULL) {
          dw->iso_var_value = 0.0;
          dw->iso_var = new double[nno];
          FOR_INO dw->iso_var[ino] = dw->simple_geom->pointSignedDistance(x_no[ino]);
        }
      }
      //if (mpi_rank == 0) cout << " > building iso-surface for " <<  iso_var_name << " = " << iso_var_value << endl;

      if (dw->nscalars_cv > 0) {
        dn_no = new double[nno*dw->nscalars_cv];      // reallocated/shuffled in buildIsoWithWgtsAndData
        for (int ivar = 0; ivar < dw->nscalars_cv; ++ivar) {
          setNoDN(dn_no+nno*ivar,dw->scalar_names_cv[ivar]);
        }
      }
      if (dw->nvectors_cv > 0) {
        dn3_no = new double[nno*dw->nvectors_cv][3]; // reallocated/shuffled in buildIsoWithWgtsAndData
        for (int ivar = 0; ivar < dw->nvectors_cv; ++ivar) {
          setNoDN3(dn3_no+nno*ivar,dw->vector_names_cv[ivar]);
        }
      }

      buildIsoWithWgtsAndData(triVec,dw->iso_var,dw->iso_var_value,dn_tmp,dn_no,dw->nscalars_cv,dn3_tmp,dn3_no,dw->nvectors_cv);
    }

    // reduce using map of points bounding intersection point (LOCALLY)
    int nst = triVec.size();
    map<const pair<int,int>, int> edgeMap;
    map<const pair<int,int>, int>::iterator iter;
    int (*spost)[3] = new int[nst][3];
    int nsp = 0;
    for (int itri = 0; itri < nst; ++itri) {
      pair<int,int> pair0 = pair<int,int>(min(triVec[itri].i0_l,triVec[itri].i0_r),max(triVec[itri].i0_l,triVec[itri].i0_r));
      iter = edgeMap.find(pair0);
      if (iter == edgeMap.end()) {
        edgeMap[pair0] = nsp;
        spost[itri][0] = nsp++;
      }
      else {
        spost[itri][0] = iter->second;
      }
      pair<int,int> pair1 = pair<int,int>(min(triVec[itri].i1_l,triVec[itri].i1_r),max(triVec[itri].i1_l,triVec[itri].i1_r));
      iter = edgeMap.find(pair1);
      if (iter == edgeMap.end()) {
        edgeMap[pair1] = nsp;
        spost[itri][1] = nsp++;
      }
      else {
        spost[itri][1] = iter->second;
      }
      pair<int,int> pair2 = pair<int,int>(min(triVec[itri].i2_l,triVec[itri].i2_r),max(triVec[itri].i2_l,triVec[itri].i2_r));
      iter = edgeMap.find(pair2);
      if (iter == edgeMap.end()) {
        edgeMap[pair2] = nsp;
        spost[itri][2] = nsp++;
      }
      else {
        spost[itri][2] = iter->second;
      }
    }
    edgeMap.clear();

    // write to tecplot file...
    char filename[128];
    if (dw->format == TECPLOT_FORMAT) {
      assert(dw->nvectors_cv == 0); // tecplot only writes scalar data
      sprintf(filename,"%s.%08d.plt",dw->name.c_str(),step);
      if (dw->b_dual)
        writeSimpleTriVecTecplot(filename,triVec,nst,spost,nsp,dw->scalar_names_cv,dn_no,dn_tmp,dw->nscalars_cv,ncv_d);
      else
        writeSimpleTriVecTecplot(filename,triVec,nst,spost,nsp,dw->scalar_names_cv,dn_no,dn_tmp,dw->nscalars_cv,nno);
    }
    else if ((dw->format == VTK_FORMAT)||(dw->format == PVTK_FORMAT)) {
      if (dw->format == PVTK_FORMAT) {
        CWARN("Using single-file VTK_FORMAT for writing iso-surface data.");
      }
      sprintf(filename,"%s.%08d.vtu",dw->name.c_str(),step);
      if (dw->b_dual)
        writeSimpleTriVecVTK(filename,triVec,nst,spost,nsp,dw->scalar_names_cv,dn_no,dn_tmp,dw->nscalars_cv,dw->vector_names_cv,dn3_no,dn3_tmp,dw->nvectors_cv,ncv_d);
      else
        writeSimpleTriVecVTK(filename,triVec,nst,spost,nsp,dw->scalar_names_cv,dn_no,dn_tmp,dw->nscalars_cv,dw->vector_names_cv,dn3_no,dn3_tmp,dw->nvectors_cv,nno);
    }
    else if (dw->format == ENSIGHT_FORMAT) {
      if (dw->b_dual)
        writeSimpleTriVecEnsight(dw,triVec,nst,spost,nsp,dn_no,dn_tmp,dn3_no,dn3_tmp,ncv_d,step,time);
      else
        writeSimpleTriVecEnsight(dw,triVec,nst,spost,nsp,dn_no,dn_tmp,dn3_no,dn3_tmp,nno,step,time);
    }
    else if (dw->format == STL_FORMAT) {
      WUI(WARN,"STL format ignores VARS");
      writeSimpleTriVecStl(dw,triVec,step,time);
    }
    else {
      assert(0); // should of been handled earlier
    }

    delete[] spost;
    DELETE(dn_no);
    DELETE(dn_tmp);
    DELETE(dn3_no);
    DELETE(dn3_tmp);

  }
  else if ((dw->topo == BF_TOPO)||(dw->topo == BF_TOPO_EXPLICIT)) {
    assert((dw->nscalars_lp == 0)&&(dw->nvectors_lp == 0));

    char filename[128];
    if (dw->format == TECPLOT_FORMAT) {
      assert((dw->nvectors_cv == 0)&&(dw->nvectors_bf == 0)); // tecplot only writes scalar data
      sprintf(filename,"%s.%08d.plt",dw->name.c_str(),step);
      writeFlaggedBfZonesTecplot(filename,dw->bfzone_flag,dw->scalar_names_cv,dw->nscalars_cv,
                                 dw->scalar_names_bf,dw->scalar_zones_bf,dw->nscalars_bf,dw->nzones_bf);
    }
    else if (dw->format == VTK_FORMAT) {
      sprintf(filename,"%s.%08d.vtu",dw->name.c_str(),step);
      writeFlaggedBfZonesVTK(filename,dw->bfzone_flag,dw->scalar_names_cv,dw->nscalars_cv,dw->vector_names_cv,dw->nvectors_cv,
                             dw->scalar_names_bf,dw->scalar_zones_bf,dw->nscalars_bf,
                             dw->vector_names_bf,dw->vector_zones_bf,dw->nvectors_bf,dw->nzones_bf);
    }
    else if (dw->format == PVTK_FORMAT) {
      CWARN("Using single-file VTK_FORMAT for writing boundary data.");
      sprintf(filename,"%s.%08d.vtu",dw->name.c_str(),step);
      writeFlaggedBfZonesVTK(filename,dw->bfzone_flag,dw->scalar_names_cv,dw->nscalars_cv,dw->vector_names_cv,dw->nvectors_cv,
                             dw->scalar_names_bf,dw->scalar_zones_bf,dw->nscalars_bf,
                             dw->vector_names_bf,dw->vector_zones_bf,dw->nvectors_bf,dw->nzones_bf);
    }
    else if (dw->format == ENSIGHT_FORMAT) {
      if (dw->b_dual)
        writeFlaggedBfZonesDualEnsight(dw,step,time);
      else
        writeFlaggedBfZonesEnsight(dw,step,time);
    }
    else if (dw->format == FLUENT_BP_FORMAT) {
      if (dw->cv_flag == NULL) {
        dw->cv_flag = new int8[ncv];
        FOR_ICV dw->cv_flag[icv] = 0;
        flagCvsAdjacentToZones(dw->cv_flag,dw->bfzone_flag,1);
        sprintf(filename,"%s.%08d.prof",dw->name.c_str(),step);
        writeFlaggedCvsFluentBoundaryProfile(filename,dw->cv_flag,dw->scalar_names_cv,dw->nscalars_cv);
      }
    }
    else if (dw->format == ASCII_FORMAT) {
      writeDataAscii(dw,step,time);
    }
    else if (dw->format == STL_FORMAT) {
      WUI(WARN,"STL format is unavailable for boundary data; skipping");
      return;
    }
    else {
      assert(0);
    }

  }
  else if (dw->topo == CV_TOPO) {
    assert((dw->nscalars_bf == 0)&&(dw->nvectors_bf == 0));
    assert((dw->nscalars_lp == 0)&&(dw->nvectors_lp == 0));

    if (dw->geom == ALL_GEOM) {
      if (dw->cv_flag == NULL) {
        dw->cv_flag = new int8[ncv];
        FOR_ICV dw->cv_flag[icv] = 1;
      }
    }
    else if (dw->geom == USERDEF_GEOM) {
      if (dw->iso_var == NULL) {
        assert(dw->cv_flag == NULL);
        dw->iso_var = new double[nno];
        dw->cv_flag = new int8[ncv];
      }
      // have to rebuild every time in case the iso definition is time dependent
      setNoDN(dw->iso_var,dw->iso_var_name); // constructs iso from string evaluation
      FOR_ICV dw->cv_flag[icv] = 0;
      FOR_INO dw->iso_var[ino] -= dw->iso_var_value;
      flagCvsUnderIso(dw->cv_flag,dw->iso_var,1);
    }
    else if (dw->geom == FILE_GEOM) {
      if (dw->cv_flag == NULL) {
        SurfaceShm surface;
        surface.readBinary(dw->sbin_filename);
        surface.initPointIsInside(x_cv,ncv);
        dw->cv_flag = new int8[ncv];
        FOR_ICV {
          if (surface.pointIsInside(x_cv[icv]))
            dw->cv_flag[icv] = 1;
          else
            dw->cv_flag[icv] = 0;
        }
      }
    }
    else if (dw->geom == ALL_WITH_ZONES_GEOM) {
      //no action, handled later.
    }
    else {
      assert(dw->geom == SIMPLE_GEOM);
      if (dw->cv_flag == NULL) {
        dw->cv_flag = new int8[ncv];
        FOR_ICV {
          if (dw->simple_geom->pointIsInside(x_cv[icv]))
            dw->cv_flag[icv] = 1;
          else
            dw->cv_flag[icv] = 0;
        }
      }
    }

    char filename[128];
    if (dw->format == TECPLOT_FORMAT) {
      assert(dw->nvectors_cv == 0); // tecplot only writes scalar data
      sprintf(filename,"%s.%08d.plt",dw->name.c_str(),step);
      writeFlaggedCvsTecplot(filename,dw->cv_flag,dw->scalar_names_cv,dw->nscalars_cv);
    }
    else if (dw->format == FLUENT_BP_FORMAT) {
      assert(dw->nvectors_cv == 0); // only write scalar data
      sprintf(filename,"%s.%08d.prof",dw->name.c_str(),step);
      writeFlaggedCvsFluentBoundaryProfile(filename,dw->cv_flag,dw->scalar_names_cv,dw->nscalars_cv);
    }
    else if (dw->format == VTK_FORMAT) {
      sprintf(filename,"%s.%08d.vtu",dw->name.c_str(),step);
      writeFlaggedCvsVTK(filename,dw->cv_flag,dw->scalar_names_cv,dw->nscalars_cv,dw->vector_names_cv,dw->nvectors_cv);
    }
    else if (dw->format == PVTK_FORMAT) {
      char master_filename[128];
      sprintf(master_filename,"%s.%08d.pvtu",dw->name.c_str(),step);
      sprintf(filename,"%s.%08d/p.%07d.vtu",dw->name.c_str(),step,mpi_rank);
      writeFlaggedCvsParallelVTK(master_filename,filename,dw->cv_flag,dw->scalar_names_cv,dw->nscalars_cv,dw->vector_names_cv,dw->nvectors_cv);
    }
    else if (dw->format == ENSIGHT_FORMAT) {
      // ensight case format has a transient multiple-file format, so we only need to create the
      // file once, then just edit to it later. The geo file is written once per run (handled internally),
      // because the mesh is stationary and each field scalar/vector is given their own file...
      if (dw->geom == ALL_WITH_ZONES_GEOM) {
        // special case where we write out the full cv and boundary data (requested by honda)...
        //writeCvsAndFlaggedBfZonesEnsight(dw,step,time);
        if (dw->b_dual) {
          assert(x_cv_d);
          writeCvsAndFlaggedBfZonesDualEnsightMultiPart(dw,step,time);
        }
        else {
          writeCvsAndFlaggedBfZonesEnsightMultiPart(dw,step,time);
        }
      }
      else {
        //writeFlaggedCvsEnsight(dw,dw->cv_flag,step,time);
        if (dw->b_dual) {
          assert(x_cv_d);
          writeFlaggedCvsDualEnsightMultiPart(dw,step,time);
        }
        else {
          writeFlaggedCvsEnsightMultiPart(dw,dw->cv_flag,step,time);
        }
      }
    }
    else if (dw->format == STL_FORMAT) {
      WUI(WARN,"STL format is unavailable for volume data; skipping");
      return;
    }
    else {
      assert(0);
    }

  }
  else if (dw->topo == LP_TOPO) {
    if (dw->format != ENSIGHT_FORMAT) {
      WUI(WARN," ! only ENSIGHT format is available for particle data; skipping");
      return;
    }

    assert((dw->nscalars_bf == 0)&&(dw->nvectors_bf == 0));
    assert((dw->nscalars_cv == 0)&&(dw->nvectors_cv == 0));

    DELETE(dw->lp_flag);
    const int nlp = *lpHelperVec[dw->lp_index].size_ptr;
    dw->lp_flag = new int[nlp];
    if (dw->geom == ALL_GEOM) {
      for (int ip = 0; ip < nlp; ++ip)
        dw->lp_flag[ip] = 1;
    }
    else {
      assert(0);
    }

    if (dw->format == ENSIGHT_FORMAT) {
      writeFlaggedLpsEnsight(dw,step,time);
    }
    else {
      assert(0);
    }
  }
  else {
    assert(0);
  }
}


void StaticSolver::flagHangingNodes(bool (*hnofa_b)[2], bool *hnobf_b, const int8* cv_flag){

  assert(hnofa_b);
  assert(hnobf_b);
  for (int nof = 0; nof < noofa_i[nfa]; ++nof){
    hnofa_b[nof][0] = false;
    hnofa_b[nof][1] = false;
  }
  for (int nob = 0; nob < noobf_i[nbf]; ++nob){
    hnobf_b[nob] = false;
  }

  if (checkParam("WRITE_DATA_REMOVE_HANGING_NODES")){
    bool * b_hcv_flag = new bool[ncv];
    FOR_ICV b_hcv_flag[icv] = true;

    int my_hanging_node_cv_count,hanging_node_cv_count;
    int nloop =0;
    COUT1("Flagging hanging nodes...");
    do {
      my_hanging_node_cv_count = 0;
      FOR_ICV {
        // here we assume that flagged cvs are identified by (cv_flag[icv] != 0)
        // or if cv_flag is null then go through ALL cvs
        if ((cv_flag==NULL||cv_flag[icv]!=0)&&b_hcv_flag[icv]){
          b_hcv_flag[icv] = false;
          map<const int,int> nodeEdgeCountMap;
          map<const int,int>::iterator nodeEdgeCountMapIter;
          buildNodeEdgeCountMap(icv,hnofa_b,hnobf_b,nodeEdgeCountMap);

          for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
            const int ifa = faocv_v[foc];
            int ihn = IHN_VAL;
            for (int nof = noofa_i[ifa]; nof < noofa_i[ifa+1]; ++nof){
              if (!hnofa_b[nof][ihn]){
                nodeEdgeCountMapIter = nodeEdgeCountMap.find(noofa_v[nof]);
                assert(nodeEdgeCountMapIter!=nodeEdgeCountMap.end());
                int edgeCount = nodeEdgeCountMapIter->second;
                if (edgeCount==2){
                  hnofa_b[nof][ihn] = true;
                  b_hcv_flag[icv] = true;
                }
              }
            }
          }
          for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
            const int ibf = bfocv_v[boc];
            for (int nob = noobf_i[ibf] ; nob < noobf_i[ibf+1]; ++nob){
              if (!hnobf_b[nob]){
                nodeEdgeCountMapIter = nodeEdgeCountMap.find(noobf_v[nob]);
                assert(nodeEdgeCountMapIter!=nodeEdgeCountMap.end());
                int edgeCount = nodeEdgeCountMapIter->second;
                if (edgeCount==2){
                  hnobf_b[nob] = true;
                  b_hcv_flag[icv] = true;
                }
              }
            }
          }
          if (b_hcv_flag[icv])
            ++my_hanging_node_cv_count;
        }//if (b_hcv_flag[icv])
      }//FOR_ICV
      MPI_Allreduce(&my_hanging_node_cv_count,&hanging_node_cv_count,1,MPI_INT,MPI_SUM,mpi_comm);
      COUT1("  iter " << nloop << " : found " << hanging_node_cv_count << " cvs with hanging nodes");
      ++nloop;
    } while (hanging_node_cv_count > 0);
    delete[] b_hcv_flag;
  }
}

//skips any flagged hanging nodes
void StaticSolver::buildNodeEdgeCountMap(const int icv, const bool (*hnofa_b)[2], const bool * hnobf_b, map<const int,int> &nodeEdgeCountMap){

  assert(hnofa_b);
  assert(hnobf_b);
  nodeEdgeCountMap.clear();

  map<pair<const int,const int>,int> edgeMap;
  map<pair<const int,const int>,int>::iterator iter;

  // look for edge pairs...
  for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
    const int ifa = faocv_v[foc];
    int ihn = IHN_VAL;
    int nof_l = noofa_i[ifa+1]-1;
    int nof_f = noofa_i[ifa];
    int nnof = 0;
    bool b_fa_has_hanging_node = false;
    for (int nof = nof_f; nof < noofa_i[ifa+1]; ++nof) {
      if (hnofa_b[nof][ihn]){
        b_fa_has_hanging_node = true;
        if (nof==nof_f)
          ++nof_f;
      }
      else{
        ++nnof;
        nof_l = nof;
      }
    }
    if (nnof<3){
      assert(b_fa_has_hanging_node); //skip faces that have become degenerate due to hanging nodes...
    }
    else{
      assert(nof_f<nof_l);
      int ino1 = noofa_v[nof_l];
      for (int nof = nof_f; nof <= nof_l; ++nof) {
        int ino0 = ino1;
        ino1     = noofa_v[nof];
        while (hnofa_b[nof][ihn]){           //skip known hanging nodes when building edges
          ino1 = noofa_v[++nof];
          assert(nof<=nof_l);
        }
        double fa_sign;
        if ( cvofa[ifa][0] == icv ) {
          fa_sign = 1.0;
          iter = edgeMap.find(pair<const int,const int>(ino0,ino1));
          if (iter == edgeMap.end()) {
            //not found, so add the reverse edge to edgeMap...
            edgeMap[pair<const int,const int>(ino1,ino0)] = ifa;
          }
        }
        else {
          assert( cvofa[ifa][1] == icv );
          fa_sign = -1.0;
          iter = edgeMap.find(pair<const int,const int>(ino1,ino0));
          if (iter == edgeMap.end()) {
            // not found, so add the reverse edge to edgeMap...
            edgeMap[pair<const int,const int>(ino0,ino1)] = ifa;
          }
        }
      }
    }
  }
  for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
    const int ibf = bfocv_v[boc];
    int nob_l = noobf_i[ibf+1]-1;
    int nob_f = noobf_i[ibf];
    int nnob = 0;
    bool b_bf_has_hanging_node = false;
    for (int nob = nob_f; nob < noobf_i[ibf+1]; ++nob) {
      if (hnobf_b[nob]){
        b_bf_has_hanging_node = true;
        if (nob==nob_f)
          ++nob_f;
      }
      else{
        ++nnob;
        nob_l = nob;
      }
    }
    if (nnob<3){
      assert(b_bf_has_hanging_node); //skip faces that have become degenerate due to hanging nodes...
    }
    else{
      assert(nob_f<nob_l);
      int ino1 = noobf_v[nob_l];
      for (int nob = nob_f; nob <= nob_l; ++nob) {
        int ino0 = ino1;
        ino1     = noobf_v[nob];
        while (hnobf_b[nob]){           //skip known hanging nodes when building edges
          ino1 = noobf_v[++nob];
          assert(nob<=nob_l);
        }
        iter = edgeMap.find(pair<const int,const int>(ino0,ino1));
        if (iter == edgeMap.end()) {
          // not found, so add the reverse edge to edgeMap...
          edgeMap[pair<const int,const int>(ino1,ino0)] = -ibf-1;
        }
      }
    }
  }
  //Now we have a list of all edges for the cv, count the node frequency
  assert(!edgeMap.empty());
  for (iter = edgeMap.begin(); iter!=edgeMap.end(); ++iter){
    pair<const int,const int> nodePair = iter->first;
    int ino0 = nodePair.first;
    int ino1 = nodePair.second;
    if (nodeEdgeCountMap.find(ino0)==nodeEdgeCountMap.end())
      nodeEdgeCountMap[ino0] = 1;
    else
      nodeEdgeCountMap[ino0]++;

    if (nodeEdgeCountMap.find(ino1)==nodeEdgeCountMap.end())
      nodeEdgeCountMap[ino1] = 1;
    else
      nodeEdgeCountMap[ino1]++;
  }
  edgeMap.clear();

}

void StaticSolver::buildNodeEdgeCountMap(const int icv, map<const int,int> &nodeEdgeCountMap){

  nodeEdgeCountMap.clear();

  map<pair<const int,const int>,int> edgeMap;
  map<pair<const int,const int>,int>::iterator iter;

  // look for edge pairs...
  for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
    const int ifa = faocv_v[foc];
    int ino1 = noofa_v[noofa_i[ifa+1]-1];
    for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      int ino0 = ino1;
      ino1     = noofa_v[nof];
      double fa_sign;
      if ( cvofa[ifa][0] == icv ) {
        fa_sign = 1.0;
        iter = edgeMap.find(pair<const int,const int>(ino0,ino1));
        if (iter == edgeMap.end()) {
          //not found, so add the reverse edge to edgeMap...
          edgeMap[pair<const int,const int>(ino1,ino0)] = ifa;
        }
      }
      else {
        assert( cvofa[ifa][1] == icv );
        fa_sign = -1.0;
        iter = edgeMap.find(pair<const int,const int>(ino1,ino0));
        if (iter == edgeMap.end()) {
          // not found, so add the reverse edge to edgeMap...
          edgeMap[pair<const int,const int>(ino0,ino1)] = ifa;
        }
      }
    }
  }
  for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
    const int ibf = bfocv_v[boc];
    int ino1 = noobf_v[noobf_i[ibf+1]-1];
    for (int nof = noobf_i[ibf]; nof != noobf_i[ibf+1]; ++nof) {
      int ino0 = ino1;
      ino1     = noobf_v[nof];
      iter = edgeMap.find(pair<const int,const int>(ino0,ino1));
      if (iter == edgeMap.end()) {
        // not found, so add the reverse edge to edgeMap...
        edgeMap[pair<const int,const int>(ino1,ino0)] = -ibf-1;
      }
    }
  }
  //Now we have a list of all edges for the cv, count the node frequency
  assert(!edgeMap.empty());
  for (iter = edgeMap.begin(); iter!=edgeMap.end(); ++iter){
    pair<const int,const int> nodePair = iter->first;
    int ino0 = nodePair.first;
    int ino1 = nodePair.second;
    if (nodeEdgeCountMap.find(ino0)==nodeEdgeCountMap.end())
      nodeEdgeCountMap[ino0] = 1;
    else
      nodeEdgeCountMap[ino0]++;

    if (nodeEdgeCountMap.find(ino1)==nodeEdgeCountMap.end())
      nodeEdgeCountMap[ino1] = 1;
    else
      nodeEdgeCountMap[ino1]++;
  }
  edgeMap.clear();

}

//skips flagged hanging nodes
void StaticSolver::flagCvsWithOpenFaces(const bool (*hnofa_b)[2], const bool *hnobf_b, const int8* cv_flag) {

  if (b_open_faces_cv)
    return;

  b_open_faces_cv = new bool[ncv];
  map<pair<const int,const int>,int> edgeMap;
  map<pair<const int,const int>,int>::iterator iter;

  FOR_ICV {
    b_open_faces_cv[icv] = false;
    if ((cv_flag == NULL)||(cv_flag[icv] > 0)) {
      assert(edgeMap.empty());
      bool b_has_neg_tri = false;
      // look for edge pairs...
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        int ihn = IHN_VAL;
        int nof_l = noofa_i[ifa+1]-1;
        int nof_f = noofa_i[ifa];
        int nnof = 0;
        bool b_fa_has_hanging_node = false;
        for (int nof = nof_f; nof < noofa_i[ifa+1]; ++nof) {
          if (hnofa_b[nof][ihn]){
            b_fa_has_hanging_node = true;
            if (nof==nof_f)
              ++nof_f;
          }
          else{
            ++nnof;
            nof_l = nof;
          }
        }
        if (nnof<3){
          assert(b_fa_has_hanging_node); //skip faces that have become degenerate due to hanging nodes...
        }
        else{
          assert(nof_f<nof_l);
          int ino1 = noofa_v[nof_l];
          for (int nof = nof_f; nof <= nof_l; ++nof) {
            int ino0 = ino1;
            ino1     = noofa_v[nof];
            while (hnofa_b[nof][ihn]){           //skip known hanging nodes when building edges
              ino1 = noofa_v[++nof];
              assert(nof<=nof_l);
            }
            double fa_sign;
            if ( cvofa[ifa][0] == icv ) {
              fa_sign = 1.0;
              iter = edgeMap.find(pair<const int,const int>(ino0,ino1));
              if (iter != edgeMap.end()) {
                // found a match...
                edgeMap.erase(iter);
              }
              else {
                // not found, so add the reverse edge to edgeMap...
                edgeMap[pair<const int,const int>(ino1,ino0)] = ifa;
              }
            }
            else {
              assert( cvofa[ifa][1] == icv );
              fa_sign = -1.0;
              iter = edgeMap.find(pair<const int,const int>(ino1,ino0));
              if (iter != edgeMap.end()) {
                // found a match...
                edgeMap.erase(iter);
              }
              else {
                // not found, so add the reverse edge to edgeMap...
                edgeMap[pair<const int,const int>(ino0,ino1)] = ifa;
              }
            }
          }
        }
      }
      for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
        const int ibf = bfocv_v[boc];
        int nob_l = noobf_i[ibf+1]-1;
        int nob_f = noobf_i[ibf];
        int nnob = 0;
        bool b_bf_has_hanging_node = false;
        for (int nob = nob_f; nob < noobf_i[ibf+1]; ++nob) {
          if (hnobf_b[nob]){
            b_bf_has_hanging_node = true;
            if (nob==nob_f)
              ++nob_f;
          }
          else{
            ++nnob;
            nob_l = nob;
          }
        }
        if (nnob<3){
          assert(b_bf_has_hanging_node); //skip faces that have become degenerate due to hanging nodes...
        }
        else{
          assert(nob_f<nob_l);
          int ino1 = noobf_v[nob_l];
          for (int nob = nob_f; nob <= nob_l; ++nob) {
            int ino0 = ino1;
            ino1     = noobf_v[nob];
            while (hnobf_b[nob]){           //skip known hanging nodes when building edges
              ino1 = noobf_v[++nob];
              assert(nob<=nob_l);
            }
            iter = edgeMap.find(pair<const int,const int>(ino0,ino1));
            if (iter != edgeMap.end()) {
              // found a match...
              edgeMap.erase(iter);
            }
            else {
              // not found, so add the reverse edge to edgeMap...
              edgeMap[pair<const int,const int>(ino1,ino0)] = -ibf-1;
            }
          }
        }
      }
      if (!edgeMap.empty()) {
        b_open_faces_cv[icv] = true;
        //for (iter = edgeMap.begin(); iter != edgeMap.end(); ++iter)
        //  cout << icv << " " << iter->second << " " << iter->first.first << " " << iter->first.second << endl;
      }
      else {
        b_open_faces_cv[icv] = false;
      }
      edgeMap.clear();
    }
  }
}

void StaticSolver::flagCvsWithOpenFaces(const int8* cv_flag) {

  if (b_open_faces_cv)
    return;

  b_open_faces_cv = new bool[ncv];
  map<pair<const int,const int>,int> edgeMap;
  map<pair<const int,const int>,int>::iterator iter;

  FOR_ICV {
    b_open_faces_cv[icv] = false;
    if ((cv_flag == NULL)||(cv_flag[icv] > 0)) {
      assert(edgeMap.empty());
      bool b_has_neg_tri = false;
      // look for edge pairs...
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        int ino1 = noofa_v[noofa_i[ifa+1]-1];
        for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
          int ino0 = ino1;
          ino1     = noofa_v[nof];
          double fa_sign;
          if ( cvofa[ifa][0] == icv ) {
            fa_sign = 1.0;
            iter = edgeMap.find(pair<const int,const int>(ino0,ino1));
            if (iter != edgeMap.end()) {
              // found a match...
              edgeMap.erase(iter);
            }
            else {
              // not found, so add the reverse edge to edgeMap...
              edgeMap[pair<const int,const int>(ino1,ino0)] = ifa;
            }
          }
          else {
            assert( cvofa[ifa][1] == icv );
            fa_sign = -1.0;
            iter = edgeMap.find(pair<const int,const int>(ino1,ino0));
            if (iter != edgeMap.end()) {
              // found a match...
              edgeMap.erase(iter);
            }
            else {
              // not found, so add the reverse edge to edgeMap...
              edgeMap[pair<const int,const int>(ino0,ino1)] = ifa;
            }
          }
          // check edge ordering
          //const double n[3] = TRI_NORMAL_2(x_fa[ifa],x_no[ino0],x_no[ino1]);
          //if (DOT_PRODUCT(n,n_fa[ifa]) < 0.0)
          //  b_has_neg_tri = true;
        }
      }
      //if (b_has_neg_tri) {
      //  cout << x_cv[icv][0] << "," << x_cv[icv][1] << "," << x_cv[icv][2] << endl;
      //  cout.flush();
      //}
      for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
        const int ibf = bfocv_v[boc];
        int ino1 = noobf_v[noobf_i[ibf+1]-1];
        for (int nof = noobf_i[ibf]; nof != noobf_i[ibf+1]; ++nof) {
          int ino0 = ino1;
          ino1     = noobf_v[nof];
          iter = edgeMap.find(pair<const int,const int>(ino0,ino1));
          if (iter != edgeMap.end()) {
            // found a match...
            edgeMap.erase(iter);
          }
          else {
            // not found, so add the reverse edge to edgeMap...
            edgeMap[pair<const int,const int>(ino1,ino0)] = -ibf-1;
          }
        }
      }

      if (!edgeMap.empty()) {
        b_open_faces_cv[icv] = true;
        //for (iter = edgeMap.begin(); iter != edgeMap.end(); ++iter)
        //  cout << icv << " " << iter->second << " " << iter->first.first << " " << iter->first.second << endl;
      }
      else {
        b_open_faces_cv[icv] = false;
      }
      edgeMap.clear();
    }
  }
}

void StaticSolver::writeDataAscii(DataWriter* dw,const int step,const double time) {

  // TODO: this is a minimal data writer for the ASCII format. It can only write
  // bf scalars from a single bf_zone for now, requested explicitly using the zone
  // naming ":" format. Some notes on this ASCII writer and data writing in general:
  //
  // 1. Make a WriteDataData container class like the "wid" and use that for parsing
  //    WRITE_DATA params once and managing the subsequent data writing. Probably
  //    remove the DataWriter class as a result.
  // 2. Add store and flush concepts to data writing when only one rank is
  //    involved, with round-robin write_rank: e.g. this routine
  // 3. robustly handle the mix of FAZONE zone1 zone2... VARS *:tau_wall and
  //    VARS zone1:tau_wall zone2:tau_wall. As a consequence, eliminate the
  //    BF_TOPO_EXPLICIT topology which is just a workaround for now.
  // 4. write data in the order it was requested as much as possible. E.g. intersperse
  //    vector components and scalars in the VAR order requested.
  // 5. convert FAZONE to comma-delimited list of zones. It will never be
  //    possible to convert VARS to such a list because the variables themselves
  //    may legitimately have commas, e.g. comp(y1:x_bf,1)
  //
  // To make this easier, the data registration process in CtiRegister now supports the
  // NData class that can store helpful pointers to information common to a class. This currently
  // includes an "x" (position), and "l" (lengthscale) and an "index". "x" and "index" are
  // used in this routine to ensure consistency and order of the writen data. We could add
  // csr-like structures like noobf_i/v to make the writing of other output formats (e.g.
  // tecplot, ensight, vtk) even more straight-forward, with fallbacks for writing data
  // as point-like data when these structures are not present.
  //
  // F. Ham, Feb 2020

  char filename[128];
  sprintf(filename,"%s.%08d.dat",dw->name.c_str(),step);

  if (mpi_rank == 0) cout << "StaticSolver::writeDataAscii(): " << filename << endl;

  // the idea here is to write bf data on one rank using gather/gatherv...
  // for now this has a hard assert enforcement that the data is the same
  // topology...

  vector<pair<string,double*> > dataVec;
  int * n_ptr = NULL;
  double (*x)[3] = NULL;
  int8 *index = NULL;
  for (int ii = 0; ii < dw->scalar_names_bf.size(); ++ii) {
    //if (mpi_rank == 0) cout << " > scalar: " << dw->scalar_names_bf[ii] << endl;
    CtiRegister::CtiData * data = CtiRegister::getCtiData(dw->scalar_names_bf[ii]);
    assert(data != NULL);
    assert(data->getType() == DN_DATA);
    if (n_ptr == NULL) {
      n_ptr = data->sizePtr();
      assert(*n_ptr >= 0);
    }
    else if (n_ptr != data->sizePtr()) {
      // for now, the topology is determined from the first data encountered (n_ptr)
      if (mpi_rank == 0) cout << " > WARNING: skipping VAR \"" << dw->scalar_names_bf[ii] << "\": ASCII write requires all data to have the same topology" << endl;
      continue;
    }
    double * ptr = data->getDNptr();
    assert(ptr != NULL);
    dataVec.push_back(pair<string,double*>(dw->scalar_names_bf[ii],ptr));
    if (x == NULL) {
      const bool b_got_x = data->getX(x);
      assert(b_got_x);
      assert(x != NULL);
    }
    else {
      double (*x_check)[3] = NULL;
      const bool b_got_x = data->getX(x_check);
      assert(b_got_x);
      assert(x_check != NULL);
      assert(x_check == x);
    }
    if (index == NULL) {
      const bool b_got_index = data->getGlobalIndex(index);
      assert(b_got_index);
      assert(index != NULL);
    }
    else {
      int8 *index_check = NULL;
      const bool b_got_index = data->getGlobalIndex(index_check);
      assert(b_got_index);
      assert(index_check != NULL);
      assert(index_check == index);
    }
  }
  assert(!dataVec.empty());

  // send/recv the global index first to prepare the writing order...

  int send_count = *n_ptr;
  int * recv_count = NULL;
  int write_rank = 0;
  if (mpi_rank == write_rank) {
    recv_count = new int[mpi_size];
  }

  MPI_Gather(&send_count,1,MPI_INT,recv_count,1,MPI_INT,write_rank,mpi_comm);

  int * recv_disp = NULL;
  int8 * recv_buf_int8 = NULL;
  int recv_count_sum;
  if (mpi_rank == write_rank) {
    recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_disp[rank-1] + recv_count[rank-1];
    recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    recv_buf_int8 = new int8[recv_count_sum];
  }

  MPI_Gatherv(index,send_count,MPI_INT8,recv_buf_int8,recv_count,recv_disp,MPI_INT8,write_rank,mpi_comm);

  // now the double gather will have n = (3+dataVec.size())*n, i.e. x,y,z,data1,data2...
  const int nvar = dataVec.size();
  int * order = NULL; // build an unpack order based on the global index...
  double * recv_buf_double = NULL;
  if (mpi_rank == write_rank) {
    int8 index_min = TWO_BILLION;
    int8 index_max = 0;
    for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
      index_min = min(index_min,recv_buf_int8[irecv]);
      index_max = max(index_max,recv_buf_int8[irecv]);
    }
    const int n_max = index_max-index_min+1;
    if (n_max <= recv_count_sum*2) {
      order = new int[n_max];
      for (int i = 0; i < n_max; ++i)
        order[i] = -1;
      for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
        const int i = recv_buf_int8[irecv]-index_min;
        assert((i >= 0)&&(i < n_max));
        assert(order[i] == -1); // no repeated indices
        order[i] = irecv;
      }
      // and compress...
      int ii = 0;
      for (int i = 0; i < n_max; ++i) {
        if (order[i] >= 0) {
          if (ii < i) order[ii] = order[i];
          ++ii;
        }
      }
      assert(ii == recv_count_sum);
    }
    else {
      // TODO: the sparsity of the global index is signifcant (>50%):
      // probably use a sort to produce the order array...
      assert(0);
    }

    delete[] recv_buf_int8;

    FOR_RANK {
      recv_count[rank] *= (3+nvar);
      recv_disp[rank] *= (3+nvar);
    }

    recv_buf_double = new double[recv_count_sum*(3+nvar)];

  }

  // everyone packs their stuff...

  send_count *= (3+nvar);
  double * send_buf_double = new double[send_count];
  int isend = 0;
  for (int i = 0; i < *n_ptr; ++i) {
    send_buf_double[isend++] = x[i][0];
    send_buf_double[isend++] = x[i][1];
    send_buf_double[isend++] = x[i][2];
    for (int ii = 0; ii < dataVec.size(); ++ii)
      send_buf_double[isend++] = dataVec[ii].second[i];
  }
  assert(isend == send_count);

  MPI_Gatherv(send_buf_double,send_count,MPI_DOUBLE,
              recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,
              write_rank,mpi_comm);

  delete[] send_buf_double;

  if (mpi_rank == write_rank) {

    delete[] recv_count;
    delete[] recv_disp;

    MiscUtils::mkdir_for_file(filename);
    FILE * fp = fopen(filename,"w");
    fprintf(fp,"# %s\n",dw->param_str.c_str());
    fprintf(fp,"# step %d time %18.15le\n",step,time);
    fprintf(fp,"# 1:x 2:y 3:z");
    for (int ii = 0; ii < dataVec.size(); ++ii)
      fprintf(fp," %d:%s",4+ii,dataVec[ii].first.c_str());
    fprintf(fp,"\n");
    for (int i = 0; i < recv_count_sum; ++i) {
      const int irecv = order[i]*(3+nvar);
      // write x,y,z...
      fprintf(fp,"%18.15le %18.15le %18.15le",
              recv_buf_double[irecv],
              recv_buf_double[irecv+1],
              recv_buf_double[irecv+2]);
      // and then the vars...
      for (int ivar = 0; ivar < nvar; ++ivar)
        fprintf(fp," %18.15le",recv_buf_double[irecv+3+ivar]);
      fprintf(fp,"\n");
    }
    fclose(fp);

    delete[] recv_buf_double;
    delete[] order;

  }

}

// ===========================================================================
// tecplot functions...
// ===========================================================================

int StaticSolver::packTecplotDbuf(double* dbuf, const double* no_data, const int8* no_flag,
    const int8* bf_flag, const int8* fa_flag, const int8* cv_flag) {

  int index     = 0;

  FOR_INO if (no_flag[ino] >= 0) {
    dbuf[index++] = no_data[ino];
  }

  FOR_IBF if (bf_flag[ibf] >= 0) {
    double x_bf = 0.0;
    for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
      const int ino = noobf_v[nob];
      x_bf         += no_data[ino];
    }
    x_bf /= double(noobf_i[ibf+1]-noobf_i[ibf]);
    dbuf[index++] = x_bf;
  }

  FOR_IFA if (fa_flag[ifa] >= 0) {
    double x_fa = 0.0;
    for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      const int ino = noofa_v[nof];
      x_fa         += no_data[ino];
    }
    x_fa /= double(noofa_i[ifa+1]-noofa_i[ifa]);
    dbuf[index++] = x_fa;
  }

  double* cv_wgt = new double[ncv];
  double* cv_val = new double[ncv];

  FOR_ICV if (cv_flag[icv] >= 0) {
    cv_wgt[icv] = 0.0;
    cv_val[icv] = 0.0;
    set<int> no_list;
    for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
      const int ifa = faocv_v[foc];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino = noofa_v[nof];
        if (no_list.find(ino) == no_list.end()) {
          no_list.insert(ino);
          cv_val[icv] += no_data[ino];
          cv_wgt[icv] += 1.0;
        }
      }
    }
    for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
      const int ibf = bfocv_v[boc];
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino = noobf_v[nob];
        if (no_list.find(ino) == no_list.end()) {
          no_list.insert(ino);
          cv_val[icv] += no_data[ino];
          cv_wgt[icv] += 1.0;
        }
      }
    }
  }

  FOR_ICV if(cv_flag[icv] >= 0) {
    assert(cv_wgt[icv] > 0.0);
    dbuf[index++] = cv_val[icv]/cv_wgt[icv];
  }

  delete[] cv_val;
  delete[] cv_wgt;
  return index;
}

int StaticSolver::packCellCenteredDbuf(double* dbuf, const CtiRegister::CtiData* cv_data, const int8* cv_flag) {
  int index = 0;
  FOR_ICV if (cv_flag[icv] >= 0) {
    dbuf[index++] = cv_data->dn(icv);
  }
  return index;
}

int StaticSolver::packCellCenteredDbuf(double* dbuf, const double* cv_data, const int8* cv_flag) {
  int index = 0;
  FOR_ICV if (cv_flag[icv] >= 0) {
    dbuf[index++] = cv_data[icv];
  }
  return index;
}

void StaticSolver::writeFlaggedCvsFluentBoundaryProfile(const string& filename,int8* cv_flag,const vector<string>& var_vec,const int nvar) {

  COUT1("StaticSolver::writeFlaggedCvsFluentBoundaryProfile(" << filename << ")");

  int8 my_count = 0;

  FOR_ICV {
    if (cv_flag[icv] != 0) {
      cv_flag[icv] = my_count++;
    }
    else {
      cv_flag[icv] = -1;
    }
  }

  int8 my_disp; // cv displacement..
  MPI_Scan(&my_count,&my_disp,1,MPI_INT8,MPI_SUM,mpi_comm);
  my_disp -= my_count;

  // put our global numbers into the cv_flag
  FOR_ICV if (cv_flag[icv] >= 0)
    cv_flag[icv] += my_disp;

  int8 cv_count; // everyone gets a copy of the cv counts
  MPI_Allreduce(&my_count,&cv_count,1,MPI_INT8,MPI_SUM,mpi_comm);

  if (mpi_rank == 0) cout << "cell count: " << int(cv_count) << endl;

  if (cv_count > TWO_BILLION) {
    CERR(" > fluent boundary profile io cannot support >2B cell count currently");
  }
  else if (cv_count == 0) {
    WUI(WARN," ! no cells properly flagged for i/o; skipping");
    return;
  }

  // =============================================
  // ready to write!
  // =============================================
  MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...

  MPI_File fh;
  if (MPI_File_open(mpi_comm,(char*)filename.c_str(),
        MPI_MODE_WRONLY | MPI_MODE_CREATE,
        MPI_INFO_NULL,&fh) != 0) {
    cerr << "Error: could not open : " << filename << endl;
    throw(0);
  }

  // for the header, we construct the file on the rank 0 processor...
  const int cbuf_size = 2048; // should never be bigger than this
  char * cbuf         = new char[cbuf_size];

  int n = sprintf(cbuf,"((test point %lld)\n",cv_count);

  // only rank 0 writes the header...
  MPI_Offset offset = 0;
  if (mpi_rank == 0) MPI_File_write_at(fh,offset,cbuf,n,MPI_CHAR,MPI_STATUS_IGNORE);
  offset += n;

  // we all write data at the cell centers...
  double* dbuf    = new double[cv_count];
  double* cv_data = new double[ncv];

  // start with x,y,z...
  FOR_I3 {
    FOR_ICV cv_data[icv] = x_cv[icv][i];
    int index = packCellCenteredDbuf(dbuf,cv_data,cv_flag);
    assert(index == my_count);

    // variable header
    // only rank=0 writes, but all need to know about offset
    switch (i) {
      case 0:
        n = sprintf(cbuf,"(x\n");
        break;
      case 1:
        n = sprintf(cbuf,"(y\n");
        break;
      case 2:
        n = sprintf(cbuf,"(z\n");
        break;
      default:
        assert(0);
    }
    if (mpi_rank == 0) MPI_File_write_at(fh,offset,cbuf,n,MPI_CHAR,MPI_STATUS_IGNORE);
    offset += n;

    writeChunkedDataAscii<double>(fh,offset,dbuf,my_count,mpi_comm,"%18.15f ");  // offset globally updated within chunked write

    n = sprintf(cbuf,")\n");
    if (mpi_rank == 0) MPI_File_write_at(fh,offset,cbuf,n,MPI_CHAR,MPI_STATUS_IGNORE);
    offset += n;
  }

  // requested variables are now written ..
  for (int ivar = 0; ivar < nvar; ++ivar) {
    // variable name remapping
    string var_write_name = var_vec[ivar];  // header string
    if (var_write_name == "p") var_write_name = "absolute-pressure";
    else if (var_write_name == "p_avg") var_write_name = "mean-absolute-pressure";
    else if (var_write_name == "T") var_write_name = "temperature";
    else if (var_write_name == "T_avg") var_write_name = "mean-temperature";
    else if (var_write_name == "rho") var_write_name = "density";
    else if (var_write_name == "rho_avg") var_write_name = "mean-density";
    else if (var_write_name == "comp(u,0)") var_write_name = "x-velocity";
    else if (var_write_name == "comp(u,1)") var_write_name = "y-velocity";
    else if (var_write_name == "comp(u,2)") var_write_name = "z-velocity";
    else if (var_write_name == "comp(u_avg,0)") var_write_name = "mean-x-velocity";
    else if (var_write_name == "comp(u_avg,1)") var_write_name = "mean-y-velocity";
    else if (var_write_name == "comp(u_avg,2)") var_write_name = "mean-z-velocity";
    else if ((var_write_name == "comp(u_rms,0)*comp(u_rms,0)") || (var_write_name == "pow(comp(u_rms,0),2)")) var_write_name = "uu-reynolds-stress";
    else if ((var_write_name == "comp(u_rms,1)*comp(u_rms,1)") || (var_write_name == "pow(comp(u_rms,1),2)")) var_write_name = "vv-reynolds-stress";
    else if ((var_write_name == "comp(u_rms,2)*comp(u_rms,2)") || (var_write_name == "pow(comp(u_rms,2),2)")) var_write_name = "ww-reynolds-stress";
    else if (var_write_name == "comp(u_rey,0)") var_write_name = "vw-reynolds-stress";
    else if (var_write_name == "comp(u_rey,1)") var_write_name = "uw-reynolds-stress";
    else if (var_write_name == "comp(u_rey,2)") var_write_name = "uv-reynolds-stress";
    else if (var_write_name == "Z") var_write_name = "mixture-fraction";
    else if (var_write_name == "Z_avg") var_write_name = "mean-mixture-fraction";
    else if (var_write_name == "C") var_write_name = "progress-variable";
    else if (var_write_name == "C_avg") var_write_name = "mean-progress-variable";

    CtiRegister::CtiData * cv_data = CtiRegister::getUnregisteredCtiData(var_vec[ivar]);
    n = sprintf(cbuf,"(%s\n",var_write_name.c_str());
    if (mpi_rank == 0) MPI_File_write_at(fh,offset,cbuf,n,MPI_CHAR,MPI_STATUS_IGNORE);
    offset += n;

    int index = packCellCenteredDbuf(dbuf,cv_data,cv_flag);
    assert(index == my_count);

    writeChunkedDataAscii<double>(fh,offset,dbuf,my_count,mpi_comm,"%18.15f ");  // offset incremented inside here

    n = sprintf(cbuf,")\n");
    if (mpi_rank == 0) MPI_File_write_at(fh,offset,cbuf,n,MPI_CHAR,MPI_STATUS_IGNORE);
    offset += n;
  }

  // close file
  n = sprintf(cbuf,")\n");
  if (mpi_rank == 0) MPI_File_write_at(fh,offset,cbuf,n,MPI_CHAR,MPI_STATUS_IGNORE);
  offset += n;
  delete[] cbuf;

  delete[] cv_data;
  delete[] dbuf;

  MPI_File_set_size(fh,offset);
  MPI_File_close(&fh);

}

// write flagged cv's to tecplot
void StaticSolver::writeFlaggedCvsTecplot(const string& filename,int8* cv_flag,const vector<string>& var_vec,const int nvar) {

  COUT1("StaticSolver::writeFlaggedCvsTecplot() " + filename);

  int8* no_flag = new int8[nno];
  int8* fa_flag = new int8[nfa];
  int8* bf_flag = new int8[nbf];

  FOR_INO no_flag[ino] = mpi_size;
  FOR_IFA fa_flag[ifa] = mpi_size;
  FOR_IBF bf_flag[ibf] = mpi_size;

  FOR_ICV {
    if (cv_flag[icv] != 0) { // this cv gets dumped

      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        fa_flag[ifa]  = mpi_rank-mpi_size;

        for(int nof = noofa_i[ifa]; nof != noofa_i[ifa+1];++nof) {
          const int ino = noofa_v[nof];
          no_flag[ino]  = mpi_rank;
        }
      }
    }
  }

  FOR_IBF {
    const int icv0 = cvobf[ibf];
    if (cv_flag[icv0] != 0) {
      bf_flag[ibf] = mpi_rank-mpi_size;
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino = noobf_v[nob];
        no_flag[ino]  = mpi_rank;
      }
    }
  }

  int8 my_count[2] = {0,0};

  updateNoData(no_flag,MIN_NO_PERIODIC_DATA);
  FOR_INO {
    if (no_flag[ino] == mpi_rank)
      no_flag[ino] = my_count[0]++;     // this node is getting dumped
    else if (no_flag[ino] < mpi_size)
      no_flag[ino] = -1;                // this node gets dumped, but not by us
    else {
      assert(no_flag[ino] == mpi_size);
      no_flag[ino] = -2;                // nobody owns this node
    }
  }

  // need to tiebreak faces, but we need to go through the cells
  // since we don't have a face communicator at this stage.
  updateFaData(fa_flag,MIN_NO_PERIODIC_DATA);

  FOR_IBF {
    if (bf_flag[ibf] == mpi_rank-mpi_size) {
      bf_flag[ibf] = my_count[0]++; // we own this face and it gets dumped
    } else {
      assert(bf_flag[ibf] == mpi_size); // we own this face, but its not getting dumped
      bf_flag[ibf] = -2;
    }
  }

  FOR_IFA {
    if (fa_flag[ifa] == mpi_rank-mpi_size)
      fa_flag[ifa] = my_count[0]++; // we own this face, and it does get dumped.
    else if (fa_flag[ifa] < 0)
      fa_flag[ifa] = -1; // somebody owns this face and it is getting dumped
    else if (fa_flag[ifa] == mpi_rank)
      fa_flag[ifa] = -2; // we own this face, but it is not getting dumped.
    else if (fa_flag[ifa] < mpi_size)
      fa_flag[ifa] = -2; // someone else owns this face, but it is not getting dumped.
    else {
      assert(fa_flag[ifa] == mpi_size);
      fa_flag[ifa] = -2;
    }
  }

  FOR_ICV {
    if (cv_flag[icv] == 0) {
      cv_flag[icv] = -2; // skip completely
    } else {
      int ne = 0;
      // tesselate the cell and dump a point at the centroid..
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa  = faocv_v[foc];
        const int nnof = noofa_i[ifa+1] - noofa_i[ifa];
        ne            += nnof;
      }
      cv_flag[icv] = my_count[0]++;
      my_count[1] += ne;
    }
  }

  // count the additional elements associated with the boundary faces..
  FOR_IBF {
    const int icv = cvobf[ibf];
    if (cv_flag[icv] != -2) {
      assert(cv_flag[icv] >= 0);
      my_count[1] += noobf_i[ibf+1] - noobf_i[ibf];
    }
  }

  int8 my_disp[2]; // nodal, element displacements..
  MPI_Scan(my_count,my_disp,2,MPI_INT8,MPI_SUM,mpi_comm);
  FOR_I2 my_disp[i] -= my_count[i];

  // put our global numbers into the no_flag, fa_flag, cv_flag
  FOR_INO if (no_flag[ino] >= 0)
    no_flag[ino] += my_disp[0];

  FOR_IBF if (bf_flag[ibf] >= 0)
    bf_flag[ibf] += my_disp[0];

  FOR_IFA if (fa_flag[ifa] >= 0)
    fa_flag[ifa] += my_disp[0];

  FOR_ICV if (cv_flag[icv] >= 0)
    cv_flag[icv] += my_disp[0];

  int8 count[2]; // everyone gets a copy of the node,element counts
  MPI_Allreduce(my_count,count,2,MPI_INT8,MPI_SUM,mpi_comm);

  if (mpi_rank == 0)
    cout << "Node count: " << int(count[0]) << " Element count: " << int(count[1]) << endl;

  if ((count[0] > TWO_BILLION)||(count[1] > TWO_BILLION)) {
    CERR(" > tecplot io cannot support >2B node/element count currently");
  }

  // for the header, we construct the file on the rank 0 processor...
  const int cbuf_size = 2048; // should never be bigger than this
  char * cbuf         = new char[cbuf_size];

  // 8 chars...
  int n = sprintf(cbuf,"#!TDV75 ");

  // 1 int...
  int ibuf[5]; // used below
  ibuf[0] = 1;
  n += TecplotIO::writeIntsTecplot(cbuf+n,ibuf,1);

  // title is the filename...
  n += TecplotIO::writeCharsTecplot(cbuf+n,filename.c_str());

  // variable count...
  int num_var = 3+nvar; // for X,Y,Z
  n += TecplotIO::writeIntsTecplot(cbuf+n,&num_var,1);

  // variable names...
  n += TecplotIO::writeCharsTecplot(cbuf+n,"x");
  n += TecplotIO::writeCharsTecplot(cbuf+n,"y");
  n += TecplotIO::writeCharsTecplot(cbuf+n,"z");

  string my_string;
  for (int ivar = 0; ivar < nvar; ++ivar) {
    my_string = var_vec[ivar];
    if (my_string == "x" || my_string == "y" || my_string == "z" || my_string == "X" || my_string == "Y" || my_string == "Z") {
      my_string += "-scalar";
    }
    n += TecplotIO::writeCharsTecplot(cbuf+n,my_string.c_str());
  }

  // zone marker...
  float flt = 299.0;
  n += TecplotIO::writeFloatsTecplot(cbuf+n,&flt,1);

  // for the zone name, use the filename...
  n += TecplotIO::writeCharsTecplot(cbuf+n,filename.c_str());

  // block data...
  ibuf[0] = 2;      // format = FEBLOCK
  ibuf[1] = -1;     // Zone color - tecplot choose
  ibuf[2] = int(count[0]); // global number of "points"
  ibuf[3] = int(count[1]); // global number of "elements"
  ibuf[4] = 3;      // Bricks
  n      += TecplotIO::writeIntsTecplot(cbuf+n,ibuf,5);

  // demarkation between header and data...
  flt = 357.0;
  n += TecplotIO::writeFloatsTecplot(cbuf+n,&flt,1);

  // data...
  flt = 299.0;
  n += TecplotIO::writeFloatsTecplot(cbuf+n,&flt,1);

  ibuf[0] = 0; // number of "repeat variables?"
  n += TecplotIO::writeIntsTecplot(cbuf+n,ibuf,1);

  // the precision for each variable...
  int precision = 2; // 1 - single precision, 2 - double precision
  for (int i = 0; i < num_var; ++i)
    n += TecplotIO::writeIntsTecplot(cbuf+n,&precision,1);
  assert(n <= cbuf_size);

  // =============================================
  // ready to write!
  // =============================================
  MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...
  MPI_File_delete((char*)filename.c_str(),MPI_INFO_NULL);

  MPI_File fh;
  if (MPI_File_open(mpi_comm,(char*)filename.c_str(),
        MPI_MODE_WRONLY | MPI_MODE_CREATE,
        MPI_INFO_NULL,&fh) != 0) {
    cerr << "Error: could not open : " << filename << endl;
    throw(0);
  }

  MPI_Offset offset = 0;
  offset += n;

  // only rank 0 writes the header...
  if (mpi_rank == 0) MPI_File_write(fh,cbuf,n,MPI_CHAR,MPI_STATUS_IGNORE);

  delete[] cbuf;

  // we all write data at the nodes...
  double* dbuf    = new double[my_count[0]];
  double* no_data = new double[nno];

  // start with x,y,z...
  FOR_I3 {
    FOR_INO no_data[ino] = x_no[ino][i];
    int index = packTecplotDbuf(dbuf,no_data,no_flag,bf_flag,fa_flag,cv_flag);
    assert(index == my_count[0]);

    writeChunkedData<double>(fh,offset+my_disp[0]*double_size,dbuf,my_count[0],mpi_comm);
    offset += count[0]*double_size;
  }

  // requested variables are now written ..
  for (int ivar = 0; ivar < nvar; ++ivar) {
    setNoDN(no_data,var_vec[ivar]);
    int index = packTecplotDbuf(dbuf,no_data,no_flag,bf_flag,fa_flag,cv_flag);
    assert(index == my_count[0]);

    writeChunkedData<double>(fh,offset+my_disp[0]*double_size,dbuf,my_count[0],mpi_comm);
    offset += count[0]*double_size;
  }

  delete[] no_data;
  delete[] dbuf;

  // ======================================
  // connectivity...
  // ======================================
  // each element gets 8 ints + the first rank writes one additional int (a zero)
  // to mark the start of the connectivity...

  my_count[1] *= 8;
  my_disp[1]  *= 8;
  if (mpi_rank == 0)
    my_count[1] += 1;
  else
    my_disp[1] += 1;

  // also update the global node and face counts to
  // overwrite the -1 values on processors that did not
  // write data...

  updateNoData(no_flag,MAX_NO_PERIODIC_DATA);
  updateFaData(fa_flag,MAX_NO_PERIODIC_DATA);

  int * ibuffer = new int[my_count[1]];
  int index     = 0;

  // rank 0 writes a zero...
  if (mpi_rank == 0) ibuffer[index++] = 0;

  FOR_ICV {
    if (cv_flag[icv] >= 0) {
      // this cv is getting dumped by tesselating to each face...
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; foc++) {
        const int ifa   = faocv_v[foc];
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;

        { // tesselating ...
          assert(fa_flag[ifa] >= 0);
          int ino1 = noofa_v[nof_l];
          for (int nof = nof_f; nof <= nof_l; ++nof) {
            int ino0 = ino1;
            ino1     = noofa_v[nof];
            // joint to the fa and cv forming a tet...
            if (cvofa[ifa][0] == icv) {
              ibuffer[index++] = no_flag[ino1]+1; // 1-indexed
              ibuffer[index++] = no_flag[ino0]+1;
              ibuffer[index++] = fa_flag[ifa]+1;
              ibuffer[index++] = fa_flag[ifa]+1;
              ibuffer[index++] = cv_flag[icv]+1;
              ibuffer[index++] = cv_flag[icv]+1;
              ibuffer[index++] = cv_flag[icv]+1;
              ibuffer[index++] = cv_flag[icv]+1;
            } else {
              assert(cvofa[ifa][1] == icv);
              ibuffer[index++] = no_flag[ino0]+1; // 1-indexed
              ibuffer[index++] = no_flag[ino1]+1;
              ibuffer[index++] = fa_flag[ifa]+1;
              ibuffer[index++] = fa_flag[ifa]+1;
              ibuffer[index++] = cv_flag[icv]+1;
              ibuffer[index++] = cv_flag[icv]+1;
              ibuffer[index++] = cv_flag[icv]+1;
              ibuffer[index++] = cv_flag[icv]+1;
            }
          }
        }
      }
    }
  }
  // count the additional elements associated with the boundary faces..
  FOR_IBF {
    const int icv = cvobf[ibf];
    if (cv_flag[icv] >= 0) {
      assert(bf_flag[ibf] >=0);
      const int nof_f = noobf_i[ibf];
      const int nof_l = noobf_i[ibf+1]-1;
      int ino1 = noobf_v[nof_l];
      for (int nof = nof_f; nof <= nof_l; ++nof) {
        int ino0 = ino1;
        ino1     = noobf_v[nof];
        // joint to the fa and cv forming a tet...
        ibuffer[index++] = no_flag[ino1]+1; // 1-indexed
        ibuffer[index++] = no_flag[ino0]+1;
        ibuffer[index++] = bf_flag[ibf]+1;
        ibuffer[index++] = bf_flag[ibf]+1;
        ibuffer[index++] = cv_flag[icv]+1;
        ibuffer[index++] = cv_flag[icv]+1;
        ibuffer[index++] = cv_flag[icv]+1;
        ibuffer[index++] = cv_flag[icv]+1;
      }
    }
  }

  assert(index == my_count[1]);

  // and write...
  writeChunkedData<int>(fh,offset+my_disp[1]*sizeof(int),ibuffer,my_count[1],mpi_comm);

  offset += (count[1]*8+1)*sizeof(int); // recall count[1] contains the element count, so x8 and add the one additional int...

  delete[] ibuffer;

  MPI_File_set_size(fh,offset);
  MPI_File_close(&fh);

  delete[] no_flag;
  delete[] fa_flag;
  delete[] bf_flag;

}

// write bf tris in flagged bf zones to tecplot
void StaticSolver::writeFlaggedBfZonesTecplot(const string& filename,bool* bfzone_flag,const vector<string>& var_vec_cv,const int nvar_cv,
                                              const vector<string>& var_vec_bf,const vector<string>& zones_vec_bf,const int nvar_bf,const int nzones_bf) {

  COUT1("StaticSolver::writeFlaggedBfZonesTecplot() " + filename);

  int8* no_flag = new int8[nno];
  int8* bf_flag = new int8[nbf];

  FOR_INO no_flag[ino] = mpi_size;
  FOR_IBF bf_flag[ibf] = mpi_size;

  for (map<const string,int>::const_iterator it = bfZoneNameMap.begin(); it != bfZoneNameMap.end(); ++it) {
    const int iz = it->second; // zone index
    if (bfzone_flag[iz]) {
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int ino = noobf_v[nob];
          assert((ino >= 0)&&(ino < nno));
          no_flag[ino] = mpi_rank;
        }
        int nnof = noobf_i[ibf+1]-noobf_i[ibf];
        assert( nnof >= 3 ); //check later ...
        if ( (nnof == 3) || (nnof == 4) )
          bf_flag[ibf] = mpi_rank;
        else
          bf_flag[ibf] = mpi_rank-mpi_size;
      }
    }
  }

  int8 my_count[2] = {0, 0}; // local ncount and ecount

  updateNoData(no_flag,MIN_NO_PERIODIC_DATA); // smallest rank wins

  FOR_INO {
    if (no_flag[ino] == mpi_rank)
      no_flag[ino] = my_count[0]++;     // this node is getting dumped
    else if (no_flag[ino] < mpi_size)
      no_flag[ino] = -1;                // this node gets dumped, but not by us
    else {
      assert(no_flag[ino] == mpi_size);
      no_flag[ino] = -2;                // nobody owns this node
    }
  }

  // need to tiebreak faces, but we need to go through the cells
  // since we don't have a face communicator at this stage.
  //updateBfData(bf_flag,MIN_NO_PERIODIC_DATA);

  FOR_IBF {
    if (bf_flag[ibf] == mpi_rank-mpi_size) {
      // this is a face that must be triangulated...
      bf_flag[ibf] = my_count[0]++;
      my_count[1] +=  noobf_i[ibf+1]-noobf_i[ibf];
    }
    else if (bf_flag[ibf] == mpi_rank) {
      bf_flag[ibf] = -1; // this face does not need to be tessellated, but needs to be dumped ...
      my_count[1] += 1;
    }
    else {
      assert(bf_flag[ibf] == mpi_size); // it will fail if smaller rank wins ...
      bf_flag[ibf] = -2; // this face is not dumped
    }
  }

  // figure out our nodal and element displacement...

  int8 my_disp[2];
  MPI_Scan(my_count,my_disp,2,MPI_INT8,MPI_SUM,mpi_comm);
  my_disp[0] -= my_count[0];
  my_disp[1] -= my_count[1];

  for (int ino = 0; ino < nno; ++ino)
    if (no_flag[ino] >= 0)
      no_flag[ino] += my_disp[0];

  FOR_IBF {
    if (bf_flag[ibf] >= 0)
      bf_flag[ibf] += my_disp[0];
  }

  // everyone gets the counts...

  int8 count[2]={0,0};
  MPI_Allreduce(my_count,count,2,MPI_INT8,MPI_SUM,mpi_comm);

  if ((count[0] > TWO_BILLION)||(count[1] > TWO_BILLION)) {
    CERR(" > tecplot io cannot support >2B node/element count currently");
  }

  // for the header, we construct the file on the rank 0 processor...

  const int cbuf_size = 4096;
  char * cbuf = new char[cbuf_size];

  // 8 chars...

  int n = sprintf(cbuf,"#!TDV75 ");

  // 1 int...

  int ibuf[5]; // used below
  ibuf[0] = 1;
  n += TecplotIO::writeIntsTecplot(cbuf+n,ibuf,1);

  // title is the filename...
  n += TecplotIO::writeCharsTecplot(cbuf+n,filename.c_str());

  // variable count...
  int num_var = 3+nvar_cv; // for X,Y,Z
  if (nvar_bf > 0) {
    assert(nzones_bf > 0);
    assert(nvar_bf%nzones_bf == 0);
    num_var += nvar_bf/nzones_bf; // global boundary vars
  }
  n += TecplotIO::writeIntsTecplot(cbuf+n,&num_var,1);

  // variable names...
  n += TecplotIO::writeCharsTecplot(cbuf+n,"x");
  n += TecplotIO::writeCharsTecplot(cbuf+n,"y");
  n += TecplotIO::writeCharsTecplot(cbuf+n,"z");

  string my_string;
  for (int ivar = 0; ivar < nvar_cv; ++ivar) {
    my_string = var_vec_cv[ivar];
    if (my_string == "x" || my_string == "y" || my_string == "z" || my_string == "X" || my_string == "Y" || my_string == "Z") {
      my_string += "-scalar";
    }
    n += TecplotIO::writeCharsTecplot(cbuf+n,my_string.c_str());
  }
  if (nvar_bf > 0) {
    assert(nzones_bf > 0);
    assert(nvar_bf%nzones_bf == 0);
    // only do each var once (not for every zone)
    // sorted by zone, so just write out entire
    for (int ivar = 0; ivar < (nvar_bf/nzones_bf); ++ivar) {
      my_string = var_vec_bf[ivar];
      size_t found = my_string.find(zones_vec_bf[ivar]);
      assert(found != std::string::npos);
      my_string = MiscUtils::replaceAll(my_string,zones_vec_bf[ivar]+":","");
      if (my_string == "x" || my_string == "y" || my_string == "z" || my_string == "X" || my_string == "Y" || my_string == "Z") {
        my_string += "-scalar";
      }
      n += TecplotIO::writeCharsTecplot(cbuf+n,my_string.c_str());
    }
  }

  // zone marker...
  float flt = 299.0;
  n += TecplotIO::writeFloatsTecplot(cbuf+n,&flt,1);

  // for the zone name, use the filename...
  n += TecplotIO::writeCharsTecplot(cbuf+n,filename.c_str());

  // block data...
  ibuf[0] = 2;      // format = FEBLOCK
  ibuf[1] = -1;     // Zone color - tecplot choose
  ibuf[2] = int(count[0]); // global number of "points"
  ibuf[3] = int(count[1]); // global number of "elements"
  ibuf[4] = 1;      // 1=Quald 3= Bricks
  n      += TecplotIO::writeIntsTecplot(cbuf+n,ibuf,5);

  // demarkation between header and data...
  flt = 357.0;
  n += TecplotIO::writeFloatsTecplot(cbuf+n,&flt,1);

  // data...
  flt = 299.0;
  n += TecplotIO::writeFloatsTecplot(cbuf+n,&flt,1);

  ibuf[0] = 0; // number of "repeat variables?"
  n += TecplotIO::writeIntsTecplot(cbuf+n,ibuf,1);

  // the precision for each variable...
  int precision = 2; // 1 - single precision, 2 - double precision
  for (int i = 0; i < num_var; ++i)
    n += TecplotIO::writeIntsTecplot(cbuf+n,&precision,1);
  assert(n <= cbuf_size);

  // =============================================
  // ready to write!
  // =============================================
  MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...
  MPI_File_delete((char*)filename.c_str(),MPI_INFO_NULL);

  MPI_File fh;
  if (MPI_File_open(mpi_comm,(char*)filename.c_str(),
        MPI_MODE_WRONLY | MPI_MODE_CREATE,
        MPI_INFO_NULL,&fh) != 0) {
    cerr << "Error: could not open : " << filename.c_str() << endl;
    throw(0);
  }

  MPI_Offset offset = 0;
  offset += n;

  // only rank 0 writes the header...
  if (mpi_rank == 0) MPI_File_write(fh,cbuf,n,MPI_CHAR,MPI_STATUS_IGNORE);

  delete[] cbuf;

  // we all write data at the nodes...
  double* dbuf = new double[my_count[0]];

  // start with x,y,z...
  FOR_I3 {
    int index = 0;
    for (int ino = 0; ino < nno; ++ino) {
      if (no_flag[ino] >= 0) {
        dbuf[index++] = x_no[ino][i];
      }
    }
    FOR_IBF {
      if(bf_flag[ibf] >=0) {
        double x_fa = 0.0;
        int nob_f = noobf_i[ibf];
        int nob_l = noobf_i[ibf+1]-1;
        for (int nob = nob_f; nob <= nob_l; nob++) {
          const int ino = noobf_v[nob];
          assert((ino >= 0)&&(ino < nno));
          x_fa += x_no[ino][i]; // ith component only
        }
        x_fa /= (double)(nob_l-nob_f+1);
        dbuf[index++] = x_fa;
      }
    }
    assert( index == my_count[0] );

    writeChunkedData<double>(fh,offset+my_disp[0]*double_size,dbuf,my_count[0],mpi_comm);
    offset += count[0]*double_size;

  }

  // now the var data...
  if (nvar_cv > 0) {
    double* no_data = new double[nno];
    for (int ivar = 0; ivar < nvar_cv; ++ivar) {
      setNoDN(no_data,var_vec_cv[ivar]);

      int index = 0;
      for (int ino = 0; ino < nno; ++ino) {
        if (no_flag[ino] >= 0) {
          dbuf[index++] = no_data[ino];
        }
      }
      for (int ibf = 0; ibf < nbf; ++ibf) {
        if (bf_flag[ibf] >= 0) {
          double phi_fa = 0.0;
          int nob_f = noobf_i[ibf];
          int nob_l = noobf_i[ibf+1]-1;
          for (int nob = nob_f; nob <= nob_l; nob++) {
            int ino = noobf_v[nob];
            phi_fa += no_data[ino];
          }
          phi_fa /= (double)(nob_l-nob_f+1);
          dbuf[index++] = phi_fa;
        }
      }
      assert(index == my_count[0]);

      writeChunkedData<double>(fh,offset+my_disp[0]*double_size,dbuf,my_count[0],mpi_comm);
      offset += count[0]*8;

    }
    delete[] no_data;
  }
  if (nvar_bf > 0) {
    double* bf_data = new double[nbf];
    double* no_data = new double[nno];
    vector<string> var_subvec_bf(nzones_bf);
    for (int ivar = 0; ivar < nvar_bf/nzones_bf; ++ivar) {

      for (int izone = 0; izone < nzones_bf; ++izone) {
        var_subvec_bf[izone] = var_vec_bf[izone*(nvar_bf/nzones_bf)+ivar];
      }
      setBfDN(no_data,var_subvec_bf,bf_data);

      int index = 0;
      for (int ino = 0; ino < nno; ++ino) {
        if (no_flag[ino] >= 0) {
          dbuf[index++] = no_data[ino];
        }
      }
      for (int ibf = 0; ibf < nbf; ++ibf) {
        if (bf_flag[ibf] >= 0) {
          dbuf[index++] = bf_data[ibf];
        }
      }
      assert(index == my_count[0]);

      writeChunkedData<double>(fh,offset+my_disp[0]*double_size,dbuf,my_count[0],mpi_comm);
      offset += count[0]*8;

    }
    delete[] no_data;
    delete[] bf_data;
  }
  delete[] dbuf;

  // ======================================
  // connectivity...
  // ======================================

  // each element gets 4 ints + the first rank writes one additional int (a zero)
  // to mark the start of the connectivity...

  my_count[1] *= 4;
  my_disp[1]  *= 4;
  if (mpi_rank == 0)
    my_count[1] += 1;
  else
    my_disp[1] += 1;

  // also update the global node and face counts to
  // overwrite the -1 values on processors that did not
  // write data...

  updateNoData(no_flag,MAX_NO_PERIODIC_DATA);

  int * ibuffer = new int[my_count[1]];
  int index = 0;

  // rank 0 writes a zero...
  if (mpi_rank == 0) ibuffer[index++] = 0;

  FOR_IBF {
    if (bf_flag[ibf] == -1) {
      int nob_f = noobf_i[ibf];
      int nob_l = noobf_i[ibf+1]-1;
      int nnof = nob_l - nob_f + 1;
      switch (nnof) {
        case 4:
          // quad...
          ibuffer[index++] = no_flag[noobf_v[nob_f+0]]+1; // 1-indexed
          ibuffer[index++] = no_flag[noobf_v[nob_f+1]]+1;
          ibuffer[index++] = no_flag[noobf_v[nob_f+2]]+1;
          ibuffer[index++] = no_flag[noobf_v[nob_f+3]]+1;
          break;
        case 3:
          // tri...
          ibuffer[index++] = no_flag[noobf_v[nob_f+0]]+1; // 1-indexed
          ibuffer[index++] = no_flag[noobf_v[nob_f+1]]+1;
          ibuffer[index++] = no_flag[noobf_v[nob_f+2]]+1;
          ibuffer[index++] = no_flag[noobf_v[nob_f+2]]+1;
          break;
        default:
          cerr << "Error: expecting a standard primative: " << nnof << endl;
          throw(-1);
      }
    }
    else if (bf_flag[ibf] >=0) {
      // this face is getting dumped by tessellating to each edge...
      int nob_f = noobf_i[ibf];
      int nob_l = noobf_i[ibf+1]-1;
      int ino1 = noobf_v[nob_l];
      for (int nob = nob_f; nob <= nob_l; ++nob) {
        int ino0 = ino1;
        ino1 = noobf_v[nob];
        // joint to the fa and edge forming a tri...
        ibuffer[index++] = no_flag[ino0]+1; // 1-indexed
        ibuffer[index++] = no_flag[ino1]+1;
        ibuffer[index++] = bf_flag[ibf]+1;
        ibuffer[index++] = bf_flag[ibf]+1;
      }
    }
  }

  if ( index != my_count[1] )
    cout << mpi_rank << " index: " << index << " my_count: " << my_count[1] << endl;
  assert( index == my_count[1] );

  // and write...

  writeChunkedData<int>(fh,offset+my_disp[1]*sizeof(int),ibuffer,my_count[1],mpi_comm);
  offset += (count[1]*4+1)*sizeof(int); // recall count[1] contains the element count, so x4 and add the one additional int...

  // and free...
  delete[] ibuffer;

  // and trim the file to the final size...

  MPI_File_set_size(fh,offset);
  MPI_File_close(&fh);

  delete[] no_flag;
  delete[] bf_flag;
}

void StaticSolver::writeFlaggedFacesASCII(const string& filename,const int * const fa_flag) const {

  // write each face as a set of tri's, triangulated to its face centroid...

  int * no_flag = new int[nno];
  FOR_INO no_flag[ino] = 0;

  // loop on the flagged faces, and prepare local counts...
  int my_ntri = 0;
  int my_nno = 0;
  FOR_IFA {
    if (fa_flag[ifa] != 0) {
      // there is always a new node at teh face centroid...
      ++my_nno;
      // for each flagged face, the number of new tris is the same as the
      // number of edges/nodes around the noofa_i/v loop...
      my_ntri += noofa_i[ifa+1] - noofa_i[ifa];
      // also flag the nodes...
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino = noofa_v[nof];
        no_flag[ino] = 1;
      }
    }
  }

  // index the remining nodes...
  FOR_INO {
    if (no_flag[ino] == 1) {
      no_flag[ino] = my_nno++;
    }
    else {
      no_flag[ino] = -1;
    }
  }

  // use round-robin writing to write the ASCII file...

  int8 my_count[2] = { my_ntri, my_nno };
  int8 count[2];
  MPI_Reduce(my_count,count,2,MPI_INT8,MPI_SUM,0,mpi_comm);

  FILE * fp;
  if ( mpi_rank == 0 ) {
    fp = fopen(filename.c_str(),"w");
    assert(fp != NULL);
    fprintf(fp,"VARIABLES = \"x\"\n");
    fprintf(fp,"\"y\"\n");
    fprintf(fp,"\"z\"\n");
    fprintf(fp,"ZONE T=\"iso\"\n");

    // the tecplot indexing is not going to comprehend extremely large
    // meshes -- so just check that we havent exceeded the 32 bit size

    assert( count[1] < TWO_BILLION);
    assert( count[0] < TWO_BILLION);

    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",int(count[1]),int(count[0]));
  }
  else {
    int dummy;
    MPI_Status status;
    MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
    fp = fopen(filename.c_str(),"a");
    assert(fp != NULL);
  }

  // recall face centroids are indexed first on any rank...

  FOR_IFA {
    if (fa_flag[ifa] != 0) {
      // the local index of this face is...
      fprintf(fp,"%18.15le %18.15le %18.15le\n",x_fa[ifa][0],x_fa[ifa][1],x_fa[ifa][2]);
    }
  }

  // then our nodes...

  FOR_INO {
    if (no_flag[ino] >= 0) {
      fprintf(fp,"%18.15le %18.15le %18.15le\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
    }
  }

  fclose(fp);

  if ( mpi_rank < mpi_size-1 ) {
    int dummy = 1;
    MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
  }
  else if (mpi_size > 1) { // guard against serial dead-lock
    assert(mpi_rank == mpi_size-1);
    int dummy = 1;
    MPI_Send(&dummy,1,MPI_INT,0,1235,mpi_comm); // switch the tag
  }

  // now the tri connectivity...

  int ino_offset;
  if ( mpi_rank == 0 ) {
    // wait for other ranks to finish their points when not serial...
    if (mpi_size > 1) {
      MPI_Status status;
      MPI_Recv(&ino_offset,1,MPI_INT,mpi_size-1,1235,mpi_comm,&status);
    }
    ino_offset = 0; // rank 0 starts ino_offset at 0
    fp = fopen(filename.c_str(),"a");
  }
  else {
    MPI_Status status;
    MPI_Recv(&ino_offset,1,MPI_INT,mpi_rank-1,1235,mpi_comm,&status);
    fp = fopen(filename.c_str(),"a");
  }
  assert(fp != NULL);

  // output face connectivity, offset appropriately...

  // loop on the flagged faces, and prepare local counts...
  int nno_fa = 0;
  FOR_IFA {
    if (fa_flag[ifa] != 0) {
      // there is always a new node at the face centroid...
      int ino_fa = nno_fa++;
      // loop around the outside edges...
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        // the tri is ino_fa,ino0,ino1 in the offset coordinates...
        fprintf(fp,"%d %d %d\n",ino_offset+ino_fa+1,ino_offset+no_flag[ino0]+1,ino_offset+no_flag[ino1]+1);
      }
    }
  }
  delete[] no_flag;

  fclose(fp);

  ino_offset += my_nno;

  if ( mpi_rank < mpi_size-1 ) {
    MPI_Send(&ino_offset,1,MPI_INT,mpi_rank+1,1235,mpi_comm);
  }

  MPI_Barrier(mpi_comm);

}

// this version is the locally compressed binary
void StaticSolver::writeSimpleTriVecTecplot(const string& filename,const vector<SimpleTriWithWgts>& triVec,
    const int nst, const int (*spost)[3], const int nsp, const vector<string>& var_vec,const double *data_no,
    const double *data_tmp, const int nvar,const int nno) {

  COUT1("StaticSolver::writeSimpleTriVecTecplot() " + filename);

  int8 my_count[2] = {nsp,nst}; // my nodes/elements

  int8 my_disp[2]; // nodal, element displacements..
  MPI_Scan(my_count,my_disp,2,MPI_INT8,MPI_SUM,mpi_comm);
  FOR_I2 my_disp[i] -= my_count[i];

  int8 count[2]; // everyone gets a copy of the node,element counts
  MPI_Allreduce(my_count,count,2,MPI_INT8,MPI_SUM,mpi_comm);

  if ((count[0] > TWO_BILLION)||(count[1] > TWO_BILLION)) {
    CERR(" > tecplot io cannot support >2B node/element count currently");
  }

  // for the header, we construct the file on the rank 0 processor...
  const int cbuf_size = 2048; // should never be bigger than this
  char * cbuf         = new char[cbuf_size];

  // 8 chars...
  int n = sprintf(cbuf,"#!TDV75 ");

  // 1 int...
  int ibuf[5]; // used below
  ibuf[0] = 1;
  n += TecplotIO::writeIntsTecplot(cbuf+n,ibuf,1);

  // title is the filename...
  n += TecplotIO::writeCharsTecplot(cbuf+n,filename.c_str());

  // variables count...
  int num_var = nvar+3;
  n += TecplotIO::writeIntsTecplot(cbuf+n,&num_var,1);

  // variable names...
  n += TecplotIO::writeCharsTecplot(cbuf+n,"x");
  n += TecplotIO::writeCharsTecplot(cbuf+n,"y");
  n += TecplotIO::writeCharsTecplot(cbuf+n,"z");

  string my_string;
  for (int ivar = 0; ivar < nvar; ++ivar) {
    my_string = var_vec[ivar];
    if (my_string == "x" || my_string == "y" || my_string == "z" || my_string == "X" || my_string == "Y" || my_string == "Z") {
      my_string += "-scalar";
    }
    n += TecplotIO::writeCharsTecplot(cbuf+n,my_string.c_str());
  }

  // zone marker...
  float flt = 299.0;
  n += TecplotIO::writeFloatsTecplot(cbuf+n,&flt,1);

  // for the zone name, use the filename...
  n += TecplotIO::writeCharsTecplot(cbuf+n,filename.c_str());

  // block data...
  ibuf[0] = 2;      // format = FEBLOCK
  ibuf[1] = -1;     // Zone color - tecplot choose
  ibuf[2] = int(count[0]); // global number of "points"
  ibuf[3] = int(count[1]); // global number of "elements"
  ibuf[4] = 0;      // Triangles
  n      += TecplotIO::writeIntsTecplot(cbuf+n,ibuf,5);

  // demarkation between header and data...
  flt = 357.0;
  n += TecplotIO::writeFloatsTecplot(cbuf+n,&flt,1);

  // data...
  flt = 299.0;
  n += TecplotIO::writeFloatsTecplot(cbuf+n,&flt,1);

  ibuf[0] = 0; // number of "repeat variables?"
  n += TecplotIO::writeIntsTecplot(cbuf+n,ibuf,1);

  // the precision for each variable...
  int precision = 2; // 1 - single precision, 2 - double precision
  for (int i = 0; i < num_var; ++i)
    n += TecplotIO::writeIntsTecplot(cbuf+n,&precision,1);
  assert(n <= cbuf_size);

  // =============================================
  // ready to write!
  // =============================================
  MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...
  MPI_File_delete((char*)filename.c_str(),MPI_INFO_NULL);

  MPI_File fh;
  if (MPI_File_open(mpi_comm,(char*)filename.c_str(),
        MPI_MODE_WRONLY | MPI_MODE_CREATE,
        MPI_INFO_NULL,&fh) != 0) {
    cerr << "Error: could not open : " << filename << endl;
    throw(0);
  }

  MPI_Offset offset = 0;

  // only rank 0 writes the header...
  if (mpi_rank == 0) MPI_File_write(fh,cbuf,n,MPI_CHAR,MPI_STATUS_IGNORE);
  offset += n;

  delete[] cbuf;

  // we all write data at the nodes...
  double* dbuf = new double[my_count[0]];

  // start with x,y,z...
  FOR_I3 {
    for (int itri = 0; itri < nst; ++itri) {
      // just overwrite the same data (quicker than flagging?)
      dbuf[spost[itri][0]] = triVec[itri].x0[i];
      dbuf[spost[itri][1]] = triVec[itri].x1[i];
      dbuf[spost[itri][2]] = triVec[itri].x2[i];
    }

    // write
    writeChunkedData<double>(fh,offset+my_disp[0]*double_size,dbuf,my_count[0],mpi_comm);
    offset += count[0]*double_size;
  }

  // do the vars...
  for (int ivar = 0; ivar < nvar; ++ivar) {
    for (int itri = 0; itri < nst; ++itri) {
      // d0
      int i_l = triVec[itri].i0_l;
      int i_r = triVec[itri].i0_r;
      double wgt_l = triVec[itri].wgt0_l;
      double wgt_r = 1.0-wgt_l;
      double d_l, d_r;
      if (i_l < nno)
        d_l = wgt_l*data_no[i_l*nvar+ivar];
      else
        d_l = wgt_l*data_tmp[(i_l-nno)*nvar+ivar];
      if (i_r < nno)
        d_r = wgt_r*data_no[i_r*nvar+ivar];
      else
        d_r = wgt_r*data_tmp[(i_r-nno)*nvar+ivar];
      dbuf[spost[itri][0]] = d_l+d_r;
      // d1
      i_l = triVec[itri].i1_l;
      i_r = triVec[itri].i1_r;
      wgt_l = triVec[itri].wgt1_l;
      wgt_r = 1.0-wgt_l;
      if (i_l < nno)
        d_l = wgt_l*data_no[i_l*nvar+ivar];
      else
        d_l = wgt_l*data_tmp[(i_l-nno)*nvar+ivar];
      if (i_r < nno)
        d_r = wgt_r*data_no[i_r*nvar+ivar];
      else
        d_r = wgt_r*data_tmp[(i_r-nno)*nvar+ivar];
      dbuf[spost[itri][1]] = d_l+d_r;
      // d2
      i_l = triVec[itri].i2_l;
      i_r = triVec[itri].i2_r;
      wgt_l = triVec[itri].wgt2_l;
      wgt_r = 1.0-wgt_l;
      if (i_l < nno)
        d_l = wgt_l*data_no[i_l*nvar+ivar];
      else
        d_l = wgt_l*data_tmp[(i_l-nno)*nvar+ivar];
      if (i_r < nno)
        d_r = wgt_r*data_no[i_r*nvar+ivar];
      else
        d_r = wgt_r*data_tmp[(i_r-nno)*nvar+ivar];
      dbuf[spost[itri][2]] = d_l+d_r;
    }

    writeChunkedData<double>(fh,offset+my_disp[0]*double_size,dbuf,my_count[0],mpi_comm);
    offset += count[0]*double_size;
  }

  delete[] dbuf;

  // ======================================
  // connectivity...
  // ======================================
  // each element gets 3 ints + the first rank writes one additional int (a zero)
  // to mark the start of the connectivity...

  my_count[1] *= 3;
  my_disp[1]  *= 3;
  if (mpi_rank == 0)
    my_count[1] += 1;
  else
    my_disp[1] += 1;

  int * ibuffer = new int[my_count[1]];
  int index     = 0;

  // rank 0 writes a zero...
  if (mpi_rank == 0) ibuffer[index++] = 0;

  // everyone else writes the tris pts...
  for (int itri = 0; itri < nst; ++itri) FOR_I3 ibuffer[index++] = spost[itri][i]+1+my_disp[0]; //1-indexed

  // and write...
  writeChunkedData<int>(fh,offset+my_disp[1]*sizeof(int),ibuffer,my_count[1],mpi_comm);
  offset += (count[1]*3+1)*sizeof(int); // recall count[1] contains the element count, so x3 and add the one additional int...

  delete[] ibuffer;

  MPI_File_set_size(fh,offset);
  MPI_File_close(&fh);

}


// this version is the uncompressed ascii
void StaticSolver::writeSimpleTriVecTecplotASCII(const string& filename,const vector<SimpleTriWithWgts>& triVec,
    const vector<string>& var_vec,const double *data_no, const double *data_tmp, const int nvar, const int nst) {

  COUT1("StaticSolver::writeSimpleTriVecTecplotASCII() " + filename);

  // a simple handshake write of the tris in tecplot format...
  int nst_global;
  int my_nst = nst;
  MPI_Reduce(&my_nst,&nst_global,1,MPI_INT,MPI_SUM,0,mpi_comm);

  FILE * fp;
  if ( mpi_rank == 0 ) {
    fp = fopen(filename.c_str(),"w");
    assert(fp != NULL);
    fprintf(fp,"VARIABLES = \"x\"\n");
    fprintf(fp,"\"y\"\n");
    fprintf(fp,"\"z\"\n");
    string str = "comp(";
    string str0 = ",0)";
    string str1 = ",1)";
    string str2 = ",2)";
    string my_string;
    for (int ivar = 0; ivar < nvar; ++ivar) {
      my_string = var_vec[ivar];
      if (my_string.find(str) != std::string::npos) {
        my_string.replace(my_string.find(str),str.length(),"");
        if (my_string.find(str0) != std::string::npos) {
          my_string.replace(my_string.find(str0),str0.length(),"-x");
        }
        else if (my_string.find(str1) != std::string::npos) {
          my_string.replace(my_string.find(str1),str1.length(),"-y");
        }
        else {
          assert(my_string.find(str2) != std::string::npos);
          my_string.replace(my_string.find(str2),str2.length(),"-z");
        }
      }
      fprintf(fp,"\"%s\"\n",my_string.c_str());
    }
    fprintf(fp,"ZONE T=\"iso\"\n");
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nst_global*3,nst_global);
  }
  else {
    int dummy;
    MPI_Status status;
    MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
    fp = fopen(filename.c_str(),"a");
    assert(fp != NULL);
  }

  for (int it = 0; it < triVec.size(); ++it) {
    fprintf(fp,"%18.15le %18.15le %18.15le",triVec[it].x0[0],triVec[it].x0[1],triVec[it].x0[2]);
    for (int ivar = 0; ivar < nvar; ++ivar) {
      const int i_l = triVec[it].i0_l;
      const int i_r = triVec[it].i0_r;
      const double wgt_l = triVec[it].wgt0_l;
      const double wgt_r = 1.0-wgt_l;
      double d_l, d_r;
      if (i_l < nno)
        d_l = wgt_l*data_no[i_l*nvar+ivar];
      else
        d_l = wgt_l*data_tmp[(i_l-nno)*nvar+ivar];
      if (i_r < nno)
        d_r = wgt_r*data_no[i_r*nvar+ivar];
      else
        d_r = wgt_r*data_tmp[(i_r-nno)*nvar+ivar];
      fprintf(fp," %18.15le",d_l+d_r);
    }
    fprintf(fp,"\n");
    fprintf(fp,"%18.15le %18.15le %18.15le",triVec[it].x1[0],triVec[it].x1[1],triVec[it].x1[2]);
    for (int ivar = 0; ivar < nvar; ++ivar) {
      const int i_l = triVec[it].i1_l;
      const int i_r = triVec[it].i1_r;
      const double wgt_l = triVec[it].wgt1_l;
      const double wgt_r = 1.0-wgt_l;
      double d_l, d_r;
      if (i_l < nno)
        d_l = wgt_l*data_no[i_l*nvar+ivar];
      else
        d_l = wgt_l*data_tmp[(i_l-nno)*nvar+ivar];
      if (i_r < nno)
        d_r = wgt_r*data_no[i_r*nvar+ivar];
      else
        d_r = wgt_r*data_tmp[(i_r-nno)*nvar+ivar];
      fprintf(fp," %18.15le",d_l+d_r);
    }
    fprintf(fp,"\n");
    fprintf(fp,"%18.15le %18.15le %18.15le",triVec[it].x2[0],triVec[it].x2[1],triVec[it].x2[2]);
    for (int ivar = 0; ivar < nvar; ++ivar) {
      const int i_l = triVec[it].i2_l;
      const int i_r = triVec[it].i2_r;
      const double wgt_l = triVec[it].wgt2_l;
      const double wgt_r = 1.0-wgt_l;
      double d_l, d_r;
      if (i_l < nno)
        d_l = wgt_l*data_no[i_l*nvar+ivar];
      else
        d_l = wgt_l*data_tmp[(i_l-nno)*nvar+ivar];
      if (i_r < nno)
        d_r = wgt_r*data_no[i_r*nvar+ivar];
      else
        d_r = wgt_r*data_tmp[(i_r-nno)*nvar+ivar];
      fprintf(fp," %18.15le",d_l+d_r);
    }
    fprintf(fp,"\n");
  }

  fclose(fp);

  if ( mpi_rank < mpi_size-1 ) {
    int dummy = 1;
    MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
  }

  MPI_Barrier(mpi_comm);

  // now the tri connectivity...

  int tri_offset = 0;
  if ( mpi_rank == 0 ) {
    fp = fopen(filename.c_str(),"a");
  }
  else {
    MPI_Status status;
    MPI_Recv(&tri_offset,1,MPI_INT,mpi_rank-1,1235,mpi_comm,&status);
    fp = fopen(filename.c_str(),"a");
  }
  assert(fp != NULL);

  for (int it = 0; it < triVec.size(); ++it) {
    fprintf(fp,"%d %d %d\n",(tri_offset+it)*3+1,(tri_offset+it)*3+2,(tri_offset+it)*3+3);
  }

  fclose(fp);

  tri_offset += triVec.size();

  if ( mpi_rank < mpi_size-1 ) {
    MPI_Send(&tri_offset,1,MPI_INT,mpi_rank+1,1235,mpi_comm);
  }

  MPI_Barrier(mpi_comm);

}

// ===========================================================================
// vtk functions...
// ===========================================================================

// first define the face types
#define VTK_TRIANGLE 5
#define VTK_POLYGON 7
#define VTK_QUAD 9
#define VTK_TET 10
#define VTK_HEX 12
#define VTK_PRISM 13
#define VTK_PYRAMID 14
#define VTK_POLY 42

void StaticSolver::writeSimpleTriVecVTK(const string& filename,const vector<SimpleTriWithWgts>& triVec,
    const int nst,const int (*spost)[3],const int nsp, const vector<string>& scalar_names,
    const double *dn_no,const double *dn_tmp,const int nscalars, const vector<string>& vector_names,
    const double (*dn3_no)[3],const double (*dn3_tmp)[3],const int nvectors,const int nno) {

  COUT1("StaticSolver::writeSimpleTriVecVTK() " + filename);

  int my_count[2] = {nsp,nst}; // my nodes/elements

  int my_disp[2]; // nodal, element displacements..
  MPI_Scan(my_count,my_disp,2,MPI_INT,MPI_SUM,mpi_comm);
  FOR_I2 my_disp[i] -= my_count[i];

  int count[2]; // everyone gets a copy of the node,element counts
  MPI_Allreduce(my_count,count,2,MPI_INT,MPI_SUM,mpi_comm);

  int myNumNodes = nsp;
  int myNumCells = nst;
  int myNumCellNodes = 3*nst;
  int myNodeDisp = my_disp[0];
  int myCellDisp = my_disp[1];
  int myCellNodeDisp = 3*my_disp[1];
  int NumNodes = count[0];
  int NumCells = count[1];
  int NumCellNodes = 3*count[1];

  // open the file...
  MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...
  MPI_File_delete((char*)filename.c_str(),MPI_INFO_NULL);

  MPI_File fh;
  if (MPI_File_open(mpi_comm,(char*)filename.c_str(),
        MPI_MODE_WRONLY | MPI_MODE_CREATE,
        MPI_INFO_NULL,&fh) != 0) {
    cerr << "Error: could not open : " << filename << endl;
    throw(0);
  }

  // the char buffer should be big enough for the largest part...
  MPI_Offset cbuf_size = 10000 + max(4*3*myNumNodes,max(4*3*myNumCells,4*myNumCellNodes+9*myNumCells));
  char * cbuf = new char[cbuf_size];

  // Checking Endian-ness of the machine. The answer is in Endian[endian]
  const char *Endian[] = { "BigEndian", "LittleEndian" };
  const int endian = MiscUtils::getEndianness();

  // everyone prepare the header - it helps us all know the offset...
  MPI_Offset offset = 0;

  int n = sprintf(cbuf, "<?xml version=\"1.0\"?>\n");
  n += sprintf(cbuf+n, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n",Endian[endian]);
  n += sprintf(cbuf+n, "  <UnstructuredGrid>\n");
  n += sprintf(cbuf+n, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", NumNodes, NumCells);

  n += sprintf(cbuf+n, "      <Points>\n");
  n += sprintf(cbuf+n, "        <DataArray type=\"Float32\" Name=\"coordinates\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += sizeof(int) + 3 * NumNodes * sizeof(float);

  n += sprintf(cbuf+n, "      </Points>\n");

  n += sprintf(cbuf+n, "      <Cells>\n");
  n += sprintf(cbuf+n, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += sizeof(int) + NumCellNodes * sizeof(int);

  n += sprintf(cbuf+n, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += sizeof(int) + NumCells * sizeof(int);

  n += sprintf(cbuf+n, "        <DataArray type=\"Int8\" Name=\"types\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += sizeof(int) + NumCells * sizeof(char);

  n += sprintf(cbuf+n, "      </Cells>\n");

  n += sprintf(cbuf+n, "      <PointData>\n");

  string my_string;
  for (int ivar = 0; ivar < nscalars; ++ivar) {
    my_string = scalar_names[ivar];
    if (my_string == "x" || my_string == "y" || my_string == "z" || my_string == "X" || my_string == "Y" || my_string == "Z") {
      my_string += "-scalar";
    }
    my_string = MiscUtils::replaceAll(my_string,"<~","_le_");
    my_string = MiscUtils::replaceAll(my_string,"<","_lt_");
    n += sprintf(cbuf+n, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"appended\" offset=\"%lld\" />\n",
        my_string.c_str(), (long long int) offset);
    offset += sizeof(int) + NumNodes * sizeof(float);
  }
  for (int ivar = 0; ivar < nvectors; ++ivar) {
    my_string = vector_names[ivar];
    my_string = MiscUtils::replaceAll(my_string,"<~","_le_");
    my_string = MiscUtils::replaceAll(my_string,"<","_lt_");
    n += sprintf(cbuf+n, "        <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lld\" />\n",
        my_string.c_str(), (long long int) offset);
    offset += sizeof(int) + 3 * NumNodes * sizeof(float);
  }

  n += sprintf(cbuf+n, "      </PointData>\n");

  n += sprintf(cbuf+n, "    </Piece>\n");

  n += sprintf(cbuf+n, "  </UnstructuredGrid>\n");

  n += sprintf(cbuf+n, "<AppendedData encoding=\"raw\">\n_");

  // the header...
  //if (mpi_rank == 0)
  //  cout << cbuf << endl;

  int n_header = n;

  // rank 0 writes the header and the first bit of coord data. Ranks 1,2,etc
  // reset n and just write their coord data...

  if (mpi_rank != 0)
    n = 0;

  // ==============================
  // x_no...
  // ==============================

  if (mpi_rank == 0) {
    int bytes = sizeof(int) + 3 * NumNodes * sizeof(float);
    memcpy(cbuf+n,&bytes,sizeof(int));
    n += sizeof(int);
  }

  double (*xsp)[3] = new double[nsp][3];
  for (int ist = 0; ist < nst; ++ist) {
    // just overwrite the same data
    FOR_I3 xsp[spost[ist][0]][i] = triVec[ist].x0[i];
    FOR_I3 xsp[spost[ist][1]][i] = triVec[ist].x1[i];
    FOR_I3 xsp[spost[ist][2]][i] = triVec[ist].x2[i];
  }
  for (int isp = 0; isp < nsp; ++isp) {
    float xyz[3];
    FOR_I3 xyz[i] = float(xsp[isp][i]);
    memcpy(cbuf+n,xyz,sizeof(float)*3);
    n += sizeof(float)*3;
  }
  delete[] xsp;

  // should be large enough...
  assert( n <= cbuf_size );

  // now write header + x_no...
  offset = 0;

  MPI_Aint dispA;
  if (mpi_rank == 0) {
    dispA = 0;
    assert(n == n_header + sizeof(int) + sizeof(float)*myNumNodes*3);
  }
  else {
    dispA = n_header + sizeof(int) + sizeof(float)*myNodeDisp*3;
    assert(n == sizeof(float)*myNumNodes*3);
  }
  writeChunkedData<char>(fh,dispA,cbuf,n,mpi_comm);
  offset += n_header + sizeof(int) + sizeof(float)*3*NumNodes;

  // ==============================
  // connectivity...
  // this is in 3 parts: global indices, the
  // offsets, and the types. It can all be done at once...
  // ==============================

  {

    n = 0;
    int n1 = sizeof(int)*myNumCellNodes;
    if (mpi_rank == 0) n1 += sizeof(int);
    int n2 = n1 + sizeof(int)*myNumCells;
    if (mpi_rank == 0) n2 += sizeof(int);
    int nf3 = n2 + sizeof(char)*myNumCells;
    if (mpi_rank == 0) nf3 += sizeof(int);

    if (mpi_rank == 0) {

      // connectivity part...
      int bytes = sizeof(int) + NumCellNodes * sizeof(int);
      memcpy(cbuf+n,&bytes,sizeof(int));
      n += sizeof(int);

      // offset part...
      bytes = sizeof(int) + NumCells * sizeof(int);
      memcpy(cbuf+n1,&bytes,sizeof(int));
      n1 += sizeof(int);

      // type part...
      bytes = sizeof(int) + NumCells * sizeof(char);
      memcpy(cbuf+n2,&bytes,sizeof(int));
      n2 += sizeof(int);

    }

    // start the offset at our node-of-cell disp...
    int disp = myCellNodeDisp;
    for (int ist = 0; ist < nst; ++ist) {

      // connectivity...
      int ibuf[3];
      FOR_I3 ibuf[i] = spost[ist][i]+myNodeDisp;
      memcpy(cbuf+n,ibuf,sizeof(int)*3);
      n += sizeof(int)*3;

      // offset...
      disp += 3;
      memcpy(cbuf+n1,&disp,sizeof(int));
      n1 += sizeof(int);

      // type...
      cbuf[n2++] = VTK_TRIANGLE;
    }

    assert(n2 <= cbuf_size);

    int n3[3];
    MPI_Aint disp3[3];
    n3[0] = sizeof(int)*myNumCellNodes;
    if (mpi_rank == 0) {
      n3[0] += sizeof(int);
      assert(myCellNodeDisp == 0);
      disp3[0] = 0;
    }
    else {
      disp3[0] = sizeof(int)*(1+myCellNodeDisp);
    }
    n3[1] = sizeof(int)*myNumCells;
    if (mpi_rank == 0) {
      n3[1] += sizeof(int);
      assert(myCellDisp == 0);
      disp3[1] = sizeof(int)*(1+NumCellNodes);
    }
    else {
      disp3[1] = sizeof(int)*(1+NumCellNodes+1+myCellDisp);
    }
    n3[2] = sizeof(char)*myNumCells;
    if (mpi_rank == 0) {
      n3[2] += sizeof(int);
      assert(myCellDisp == 0);
      disp3[2] = sizeof(int)*(1+NumCellNodes+1+NumCells);
    }
    else {
      disp3[2] = sizeof(int)*(1+NumCellNodes+1+NumCells+1) + sizeof(char)*myCellDisp;
    }

    writeChunkedData<char>(fh,offset+disp3[0],cbuf,n3[0],mpi_comm);
    writeChunkedData<char>(fh,offset+disp3[1],cbuf+n3[0],n3[1],mpi_comm);
    writeChunkedData<char>(fh,offset+disp3[2],cbuf+n3[0]+n3[1],n3[2],mpi_comm);
  }

  offset += sizeof(int)*(1+NumCellNodes+1+NumCells+1) + sizeof(char)*NumCells;

  // we all write data at the nodes...

  // scalars...
  if (nscalars > 0) {
    if (mpi_rank == 0) {
      assert(myNodeDisp == 0);
      dispA = 0;
    }
    else {
      dispA = sizeof(int) + sizeof(float)*myNodeDisp;
    }
    double* dbuf = new double[myNumNodes];
    for (int ivar = 0; ivar < nscalars; ++ivar) {
      for (int itri = 0; itri < nst; ++itri) {
        // d0
        int i_l = triVec[itri].i0_l;
        int i_r = triVec[itri].i0_r;
        double wgt_l = triVec[itri].wgt0_l;
        double wgt_r = 1.0-wgt_l;
        double d_l, d_r;
        if (i_l < nno)
          d_l = wgt_l*dn_no[i_l*nscalars+ivar];
        else
          d_l = wgt_l*dn_tmp[(i_l-nno)*nscalars+ivar];
        if (i_r < nno)
          d_r = wgt_r*dn_no[i_r*nscalars+ivar];
        else
          d_r = wgt_r*dn_tmp[(i_r-nno)*nscalars+ivar];
        dbuf[spost[itri][0]] = d_l+d_r;
        // d1
        i_l = triVec[itri].i1_l;
        i_r = triVec[itri].i1_r;
        wgt_l = triVec[itri].wgt1_l;
        wgt_r = 1.0-wgt_l;
        if (i_l < nno)
          d_l = wgt_l*dn_no[i_l*nscalars+ivar];
        else
          d_l = wgt_l*dn_tmp[(i_l-nno)*nscalars+ivar];
        if (i_r < nno)
          d_r = wgt_r*dn_no[i_r*nscalars+ivar];
        else
          d_r = wgt_r*dn_tmp[(i_r-nno)*nscalars+ivar];
        dbuf[spost[itri][1]] = d_l+d_r;
        // d2
        i_l = triVec[itri].i2_l;
        i_r = triVec[itri].i2_r;
        wgt_l = triVec[itri].wgt2_l;
        wgt_r = 1.0-wgt_l;
        if (i_l < nno)
          d_l = wgt_l*dn_no[i_l*nscalars+ivar];
        else
          d_l = wgt_l*dn_tmp[(i_l-nno)*nscalars+ivar];
        if (i_r < nno)
          d_r = wgt_r*dn_no[i_r*nscalars+ivar];
        else
          d_r = wgt_r*dn_tmp[(i_r-nno)*nscalars+ivar];
        dbuf[spost[itri][2]] = d_l+d_r;
      }

      n = 0;
      if (mpi_rank == 0) {
        int bytes = sizeof(int) + NumNodes * sizeof(float);
        memcpy(cbuf+n,&bytes,sizeof(int));
        n += sizeof(int);
      }

      for (int isp = 0; isp < nsp; ++isp) {
        float phi = float(dbuf[isp]);
        memcpy(cbuf+n,&phi,sizeof(float));
        n += sizeof(float);
      }

      writeChunkedData<char>(fh,offset+dispA,cbuf,n,mpi_comm);
      offset += sizeof(int)+sizeof(float)*NumNodes;

    }
    delete[] dbuf;
  }

  // vectors...
  if (nvectors > 0) {
    if (mpi_rank == 0) {
      assert(myNodeDisp == 0);
      dispA = 0;
    }
    else {
      dispA = sizeof(int) + sizeof(float)*myNodeDisp*3;
    }
    double (*dbuf)[3] = new double[myNumNodes][3];
    for (int ivar = 0; ivar < nvectors; ++ivar) {
      for (int itri = 0; itri < nst; ++itri) {
        // d0
        int i_l = triVec[itri].i0_l;
        int i_r = triVec[itri].i0_r;
        double wgt_l = triVec[itri].wgt0_l;
        double wgt_r = 1.0-wgt_l;
        double d_l[3], d_r[3];
        if (i_l < nno)
          FOR_I3 d_l[i] = wgt_l*dn3_no[i_l*nvectors+ivar][i];
        else
          FOR_I3 d_l[i] = wgt_l*dn3_tmp[(i_l-nno)*nvectors+ivar][i];
        if (i_r < nno)
          FOR_I3 d_r[i] = wgt_r*dn3_no[i_r*nvectors+ivar][i];
        else
          FOR_I3 d_r[i] = wgt_r*dn3_tmp[(i_r-nno)*nvectors+ivar][i];
        FOR_I3 dbuf[spost[itri][0]][i] = d_l[i]+d_r[i];
        // d1
        i_l = triVec[itri].i1_l;
        i_r = triVec[itri].i1_r;
        wgt_l = triVec[itri].wgt1_l;
        wgt_r = 1.0-wgt_l;
        if (i_l < nno)
          FOR_I3 d_l[i] = wgt_l*dn3_no[i_l*nvectors+ivar][i];
        else
          FOR_I3 d_l[i] = wgt_l*dn3_tmp[(i_l-nno)*nvectors+ivar][i];
        if (i_r < nno)
          FOR_I3 d_r[i] = wgt_r*dn3_no[i_r*nvectors+ivar][i];
        else
          FOR_I3 d_r[i] = wgt_r*dn3_tmp[(i_r-nno)*nvectors+ivar][i];
        FOR_I3 dbuf[spost[itri][1]][i] = d_l[i]+d_r[i];
        // d2
        i_l = triVec[itri].i2_l;
        i_r = triVec[itri].i2_r;
        wgt_l = triVec[itri].wgt2_l;
        wgt_r = 1.0-wgt_l;
        if (i_l < nno)
          FOR_I3 d_l[i] = wgt_l*dn3_no[i_l*nvectors+ivar][i];
        else
          FOR_I3 d_l[i] = wgt_l*dn3_tmp[(i_l-nno)*nvectors+ivar][i];
        if (i_r < nno)
          FOR_I3 d_r[i] = wgt_r*dn3_no[i_r*nvectors+ivar][i];
        else
          FOR_I3 d_r[i] = wgt_r*dn3_tmp[(i_r-nno)*nvectors+ivar][i];
        FOR_I3 dbuf[spost[itri][2]][i] = d_l[i]+d_r[i];
      }

      n = 0;
      if (mpi_rank == 0) {
        int bytes = sizeof(int) + NumNodes * sizeof(float) * 3;
        memcpy(cbuf+n,&bytes,sizeof(int));
        n += sizeof(int);
      }

      for (int isp = 0; isp < nsp; ++isp) {
        float vec[3]; FOR_I3 vec[i] = float(dbuf[isp][i]);
        memcpy(cbuf+n,vec,sizeof(float)*3);
        n += sizeof(float)*3;
      }

      writeChunkedData<char>(fh,offset+dispA,cbuf,n,mpi_comm);
      offset += sizeof(int)+sizeof(float)*NumNodes*3;

    }
    delete[] dbuf;
  }

  // ==========================================================================
  // end of file: just rank0 writes...
  // ==========================================================================

  n = sprintf(cbuf, "\n  </AppendedData>\n</VTKFile>\n");

  if (mpi_rank == 0) {
    MPI_File_write_at(fh,offset,cbuf,n,MPI_CHAR,MPI_STATUS_IGNORE);
  }

  offset += n;

  MPI_File_set_size(fh,offset);
  MPI_File_close(&fh);

  delete[] cbuf;
}

// multi-file
void StaticSolver::writeFlaggedCvsParallelVTK(const string& master_filename,const string& filename,int8* cv_flag,
    const vector<string>& scalar_names,const int nscalars,const vector<string>& vector_names,const int nvectors) {

  COUT1("StaticSolver::writeFlaggedCvsParallelVTK() " + master_filename);

  int* no_flag = new int[nno];
  FOR_INO no_flag[ino] = mpi_size;
  FOR_ICV {
    if (cv_flag[icv] != 0) { // this cv gets dumped
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        for(int nof = noofa_i[ifa]; nof != noofa_i[ifa+1];++nof) {
          const int ino = noofa_v[nof];
          no_flag[ino]  = mpi_rank;
        }
      }
    }
  }

  FOR_IBF {
    const int icv0 = cvobf[ibf];
    if (cv_flag[icv0] != 0) {
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino = noobf_v[nob];
        no_flag[ino] = mpi_rank;
      }
    }
  }

  int my_count[4] = {0,0,0,0};

  FOR_INO {
    if (no_flag[ino] == mpi_rank)
      no_flag[ino] = my_count[0]++;     // this node is getting dumped
    else {
      assert(no_flag[ino] == mpi_size);
      no_flag[ino] = -2;                // nobody owns this node
    }
  }

  FOR_ICV {
    if (cv_flag[icv] != 0) {
      //int my_nnofoc = 0; // total number of nodes for each face of this cell
      set<int> no_list; // noocv is not available...
      my_count[3] += 1 + faocv_i[icv+1]-faocv_i[icv] + bfocv_i[icv+1]-bfocv_i[icv];
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        const int nnof = nof_l - nof_f + 1;
        my_count[3] += nnof; // add face nodes
        // all nodes have to be dumped...
        for (int nof = nof_f; nof <= nof_l; ++nof) {
          const int ino = noofa_v[nof];
          no_list.insert(ino);
        }
      }
      for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
        const int ibf = bfocv_v[boc];
        const int nob_f = noobf_i[ibf];
        const int nob_l = noobf_i[ibf+1]-1;
        const int nnob = nob_l - nob_f + 1;
        my_count[3] += nnob; // add boundary face nodes
        // all nodes have to be dumped...
        for (int nob = nob_f; nob <= nob_l; ++nob) {
          const int ino = noobf_v[nob];
          no_list.insert(ino);
        }
      }
      const int nnoc = no_list.size();
      my_count[1]++; // this cell is getting dumped
      my_count[2] += nnoc; // add cell nodes
    }
  }

  // here we assume that flagged cvs are identified by (cv_flag[icv] != 0)...
  int myNumNodes = my_count[0];
  int myNumCells = my_count[1];
  int myNumCellNodes = my_count[2];
  int myNumCellFaNo = my_count[3];

  // open the file...
  MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...

  // write control file
  if (mpi_rank == 0) {
    //Setup piece filenames, the master filename itself may be in a
    //subdirectory but this should not be written to the pvtu file
    //which lists the RELATIVE location of the vtu files.
    size_t last_slash = master_filename.find_last_of('/');
    string fn_clean;
    if (last_slash != string::npos) {
      // we now have the last position of "/"...
      int sublength = master_filename.length()-last_slash-1;
      fn_clean = master_filename.substr(last_slash+1,sublength-5);
    }
    else {
      fn_clean = master_filename.substr(0,master_filename.length()-5);
    }

    // open the file...
    FILE * fh_m;
    if ( (fh_m=fopen(master_filename.c_str(),"w"))==NULL ) {
      cerr << "Error: cannot open file " << master_filename << endl;
      throw(-1);
    }

    char * cbuf_m = new char[10000+mpi_size*100];

    // Checking Endian-ness of the machine. The answer is in Endian[endian]
    const char *Endian[] = { "BigEndian", "LittleEndian" };
    const int endian = MiscUtils::getEndianness();

    int n=sprintf(cbuf_m, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"%s\" header_type=\"UInt64\">\n",Endian[endian]);
    n += sprintf(cbuf_m+n,  "  <PUnstructuredGrid GhostLevel=\"0\">\n");

    n += sprintf(cbuf_m+n, "    <PPoints>\n");
    n += sprintf(cbuf_m+n, "      <PDataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"3\"/>\n");
    n += sprintf(cbuf_m+n, "    </PPoints>\n");

    n += sprintf(cbuf_m+n, "    <PPointData>\n");

    string my_string;
    for (int ivar = 0; ivar < nscalars; ++ivar) {
      my_string = scalar_names[ivar];
      if (my_string == "x" || my_string == "y" || my_string == "z" || my_string == "X" || my_string == "Y" || my_string == "Z") {
        my_string += "-scalar";
      }
      my_string = MiscUtils::replaceAll(my_string,"<~","_le_");
      my_string = MiscUtils::replaceAll(my_string,"<","_lt_");
      n += sprintf(cbuf_m+n, "        <PDataArray type=\"Float64\" Name=\"%s\"/>\n",
          my_string.c_str());
    }
    for (int ivar = 0; ivar < nvectors; ++ivar) {
      my_string = vector_names[ivar];
      my_string = MiscUtils::replaceAll(my_string,"<~","_le_");
      my_string = MiscUtils::replaceAll(my_string,"<","_lt_");
      n += sprintf(cbuf_m+n, "        <PDataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\"/>\n",
          my_string.c_str());
    }
    n += sprintf(cbuf_m+n, "    </PPointData>\n");


    assert(mpi_size < 10000000);
    for (int i=0;i<mpi_size;i++){
      n += sprintf(cbuf_m+n, "    <Piece Source=\"%s/p.%07d.vtu\"/>\n",fn_clean.c_str(),i);
    }

    n += sprintf(cbuf_m+n, "  </PUnstructuredGrid>\n");
    n += sprintf(cbuf_m+n, "</VTKFile>");

    //Write and close;
    fwrite(cbuf_m,1,n,fh_m);
    fclose(fh_m);

    delete[] cbuf_m;
  }

  // open the file...
  FILE * fh;

  if ( (fh=fopen(filename.c_str(),"w"))==NULL ) {
    cerr << "Error: cannot open file " << filename << endl;
    throw(-1);
  }

  // the char buffer should be big enough for the largest part...
  MPI_Offset cbuf_size = 10000; //+max(3*4*myNumNodes,8*myNumCellNodes+17*myNumCells+8*myNumCellFaNo);
  char * cbuf = new char[cbuf_size];

  // Checking Endian-ness of the machine. The answer is in Endian[endian]
  const char *Endian[] = { "BigEndian", "LittleEndian" };
  const int endian = MiscUtils::getEndianness();

  // everyone prepare the header - it helps us all know the offset...
  MPI_Offset offset = 0;

  int n = sprintf(cbuf, "<?xml version=\"1.0\"?>\n");
  n += sprintf(cbuf+n, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n",Endian[endian]);
  n += sprintf(cbuf+n, "  <UnstructuredGrid>\n");
  n += sprintf(cbuf+n, "    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n", (long long)(myNumNodes), (long long) myNumCells); // stb casting for format compliance..

  n += sprintf(cbuf+n, "      <Points>\n");
  n += sprintf(cbuf+n, "        <DataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += int_size + 3 * myNumNodes * double_size;

  n += sprintf(cbuf+n, "      </Points>\n");

  n += sprintf(cbuf+n, "      <Cells>\n");
  n += sprintf(cbuf+n, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += int_size + myNumCellNodes * int_size;

  n += sprintf(cbuf+n, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += int_size + myNumCells * int_size;

  n += sprintf(cbuf+n, "        <DataArray type=\"Int8\" Name=\"types\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += int_size + myNumCells * sizeof(char);

  // face connectivity info --> this should be of size sum_poly(1 + poly_nFace + poly_nNodes)
  n += sprintf(cbuf+n, "        <DataArray type=\"Int32\" Name=\"faces\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += int_size + myNumCellFaNo * int_size;

  // storage for the offsets of face connectivity array
  n += sprintf(cbuf+n, "        <DataArray type=\"Int32\" Name=\"faceoffsets\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += int_size + myNumCells * int_size;

  n += sprintf(cbuf+n, "      </Cells>\n");

  n += sprintf(cbuf+n, "      <PointData>\n");

  string my_string;
  for (int ivar = 0; ivar < nscalars; ++ivar) {
    my_string = scalar_names[ivar];
    if (my_string == "x" || my_string == "y" || my_string == "z" || my_string == "X" || my_string == "Y" || my_string == "Z") {
      my_string += "-scalar";
    }
    my_string = MiscUtils::replaceAll(my_string,"<~","_le_");
    my_string = MiscUtils::replaceAll(my_string,"<","_lt_");
    n += sprintf(cbuf+n, "        <DataArray type=\"Float64\" Name=\"%s\" format=\"appended\" offset=\"%lld\" />\n",
        my_string.c_str(), (long long int) offset);
    offset += int_size + myNumNodes * double_size;
  }
  for (int ivar = 0; ivar < nvectors; ++ivar) {
    my_string = vector_names[ivar];
    my_string = MiscUtils::replaceAll(my_string,"<~","_le_");
    my_string = MiscUtils::replaceAll(my_string,"<","_lt_");
    n += sprintf(cbuf+n, "        <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lld\" />\n",
        my_string.c_str(), (long long int) offset);
    offset += int_size + 3 * myNumNodes * double_size;
  }

  n += sprintf(cbuf+n, "      </PointData>\n");

  n += sprintf(cbuf+n, "    </Piece>\n");

  n += sprintf(cbuf+n, "  </UnstructuredGrid>\n");

  n += sprintf(cbuf+n, "<AppendedData encoding=\"raw\">\n_");

  // the header...
  //if (mpi_rank == 0) {
  //  cout << cbuf << endl;
  //  cout << "end - header: " << offset << endl;
  //  cout << "n " << n << endl;
  //}

  int n_header = n;

  // rank 0 writes the header and the first bit of coord data. Ranks 1,2,etc
  // reset n and just write their coord data...

  offset = 0;

  // write out header
  fwrite(cbuf,1,n,fh);
  delete[] cbuf; cbuf = NULL;

  // ==============================
  // x_no...
  // ==============================

  // now write header + x_no...
  offset += n_header; assert(sizeof(char)==1);

  double *dbuf3 = new double[myNumNodes*3];
  {
    int ind = 0;
    FOR_INO {
      if (no_flag[ino] >= 0) {
        FOR_I3 dbuf3[3*ind+i] = x_no[ino][i];
        ind++;
      }
    }
    assert(ind == myNumNodes);
  }

  int nbytes = 1*int_size+myNumNodes*3*double_size;
  fwrite(&nbytes,int_size,1,fh);
  offset += int_size*1;
  fwrite(dbuf3,double_size,myNumNodes*3,fh);
  offset += double_size*3*myNumNodes;
  if (nvectors == 0) delete[] dbuf3;

  //if (mpi_rank == 0)
  //  cout << "end of nodes: " << offset << ", rel to end of header: " << offset-n_header << endl;

  // ==============================
  // connectivity...
  // this is in 5 parts: global indices, the
  // offsets, types, faces and face offets. It can all be done at once...
  // ==============================

  // start the offset at our node-of-cell disp...
  int disp = 0;
  int fdisp = 0;
  int *int_buf = new int[myNumCellNodes];
  int *int_buf_2 = new int[myNumCells];
  assert(cbuf == NULL); cbuf = new char[myNumCells];
  int *int_buf_3 = new int[myNumCellFaNo];
  int *int_buf_4 = new int[myNumCells];
  int cv_index = 0, index = 0, index_3 = 0;
  FOR_ICV {
    if (cv_flag[icv] != 0) {
      set<int> no_list; // noocv is not available...
      int foc_f = faocv_i[icv];
      int foc_l = faocv_i[icv+1]-1;
      int nfoc = foc_l - foc_f + 1;
      for (int foc = foc_f; foc <= foc_l; ++foc) {
        const int ifa = faocv_v[foc];
        for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
          const int ino = noofa_v[nof];
          no_list.insert(ino);
        }
      }
      int boc_f = bfocv_i[icv];
      int boc_l = bfocv_i[icv+1]-1;
      int nboc = boc_l - boc_f + 1;
      for (int boc = boc_f; boc <= boc_l; ++boc) {
        const int ibf = bfocv_v[boc];
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int ino = noobf_v[nob];
          no_list.insert(ino);
        }
      }
      int nnoc = no_list.size();

      // node numbers...
      for (set<int>::iterator it = no_list.begin(); it != no_list.end(); ++it) {
        int_buf[index++] = no_flag[*it];
      }

      // offset...
      disp += nnoc;
      int_buf_2[cv_index] = disp;

      // type...
      cbuf[cv_index] = VTK_POLY;

      // face connectivity...
      // write number of faces for this cell
      int_buf_3[index_3++] = nfoc+nboc;
      fdisp += 1;

      // loop over faces of cell
      for (int foc = foc_f; foc <= foc_l; foc++) {
        int ifa   = faocv_v[foc];
        int nof_f = noofa_i[ifa];
        int nof_l = noofa_i[ifa+1]-1;
        int nnof  = nof_l - nof_f + 1;

        // number of nodes of this face
        int_buf_3[index_3++] = nnof;

        // loop over nodes of face
        for (int nof = nof_f ; nof <= nof_l; nof++) {
          int_buf_3[index_3++] = no_flag[noofa_v[nof]];
        }
        fdisp += 1+nnof;
      }
      for (int boc = boc_f; boc <= boc_l; ++boc) {
        const int ibf = bfocv_v[boc];
        int nob_f = noobf_i[ibf];
        int nob_l = noobf_i[ibf+1]-1;
        int nnob  = nob_l - nob_f + 1;

        // number of nodes of this face
        int_buf_3[index_3++] = nnob;

        // loop over nodes of face
        for (int nob = nob_f; nob <= nob_l; nob++) {
          int_buf_3[index_3++] = no_flag[noobf_v[nob]];
        }
        fdisp += 1+nnob;
      }

      // face offset...
      int_buf_4[cv_index] = fdisp;
      cv_index++;
    }
  }
  assert(cv_index == myNumCells);

  // connectivity part...
  nbytes = (1+myNumCellNodes)*int_size;
  fwrite(&nbytes,int_size,1,fh);
  offset += int_size;
  fwrite(int_buf,int_size,myNumCellNodes,fh);
  offset += int_size*myNumCellNodes;
  delete[] int_buf;

  // offset part...
  nbytes = (1+myNumCells)*int_size;
  fwrite(&nbytes,1,int_size,fh);
  offset += int_size;
  fwrite(int_buf_2,myNumCells,int_size,fh);
  offset += int_size*myNumCells;
  delete[] int_buf_2;

  // type part...
  nbytes = 1*int_size+myNumCells*1;
  fwrite(&nbytes,int_size,1,fh);
  offset += int_size;
  fwrite(cbuf,1,myNumCells,fh);
  offset += myNumCells;
  delete[] cbuf; cbuf = NULL;

  // face part...
  nbytes = (1+myNumCellFaNo)*int_size;
  fwrite(&nbytes,int_size,1,fh);
  offset += int_size;
  fwrite(int_buf_3,int_size,myNumCellFaNo,fh);
  offset += int_size*myNumCellFaNo;
  delete[] int_buf_3;

  // face offset part...
  nbytes = (1+myNumCells)*int_size;
  fwrite(&nbytes,int_size,1,fh);
  offset += int_size;
  fwrite(int_buf_4,int_size,myNumCells,fh);
  offset += int_size*myNumCells;
  delete[] int_buf_4;

  //if (mpi_rank == 0)
  //  cout << "after connectivity: offset - n_header: " << offset-n_header << endl;

  // scalars...
  if (nscalars > 0) {
    double *no_data = new double[nno];
    double *dbuf = new double[myNumNodes];
    for (int ivar = 0; ivar < nscalars; ++ivar) {
      setNoDN(no_data,scalar_names[ivar]);

      int ind = 0;
      FOR_INO {
        if (no_flag[ino] >= 0) {
          dbuf[ind++] = no_data[ino];
        }
      }
      assert(ind == myNumNodes);

      nbytes = 1*int_size+myNumNodes*double_size;
      fwrite(&nbytes,int_size,1,fh);
      offset += int_size;
      fwrite(dbuf,double_size,myNumNodes,fh);
      offset += double_size*myNumNodes;

    }
    delete[] no_data;
    delete[] dbuf;
  }

  // vectors...
  if (nvectors > 0) {
    double (*no_data)[3] = new double[nno][3];
    for (int ivar = 0; ivar < nvectors; ++ivar) {
      setNoDN3(no_data,vector_names[ivar]);

      int ind = 0;
      FOR_INO {
        if (no_flag[ino] >= 0) {
          FOR_I3 dbuf3[ind*3+i] = no_data[ino][i];
          ind++;
        }
      }
      assert(ind == myNumNodes);

      nbytes = 1*int_size+myNumNodes*3*double_size;
      fwrite(&nbytes,int_size,1,fh);
      offset += int_size;
      fwrite(dbuf3,double_size,myNumNodes*3,fh);
      offset += double_size*3*myNumNodes;
      //if (mpi_rank == 0)
      //  cout << "after no vector: " << offset << ", rel to end of header: " << offset-n_header << endl;
    }
    delete[] no_data;
    delete[] dbuf3;
  }
  delete[] no_flag;

  assert(cbuf == NULL); cbuf = new char[256];
  n = sprintf(cbuf, "\n  </AppendedData>\n</VTKFile>\n");
  delete[] cbuf;
  offset += n;
  //if (mpi_rank == 0)
  //  cout << "end size - header:" << offset-n_header << endl;

  fclose(fh);

}

// single-file
void StaticSolver::writeFlaggedCvsVTK(const string& filename,int8* cv_flag,const vector<string>& scalar_names,const int nscalars,
    const vector<string>& vector_names,const int nvectors) {

  COUT1("StaticSolver::writeFlaggedCvsVTK() " + filename);

  int8* no_flag = new int8[nno];
  FOR_INO no_flag[ino] = mpi_size;
  FOR_ICV {
    if (cv_flag[icv] != 0) { // this cv gets dumped
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        for(int nof = noofa_i[ifa]; nof != noofa_i[ifa+1];++nof) {
          const int ino = noofa_v[nof];
          no_flag[ino]  = mpi_rank;
        }
      }
    }
  }

  FOR_IBF {
    const int icv0 = cvobf[ibf];
    if (cv_flag[icv0] != 0) {
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino = noobf_v[nob];
        no_flag[ino] = mpi_rank;
      }
    }
  }

  int8 my_count[4] = {0,0,0,0};

  updateNoData(no_flag,MIN_NO_PERIODIC_DATA);
  FOR_INO {
    if (no_flag[ino] == mpi_rank)
      no_flag[ino] = my_count[0]++;     // this node is getting dumped
    else if (no_flag[ino] < mpi_size)
      no_flag[ino] = -1;                // this node gets dumped, but not by us
    else {
      assert(no_flag[ino] == mpi_size);
      no_flag[ino] = -2;                // nobody owns this node
    }
  }

  FOR_ICV {
    if (cv_flag[icv] != 0) {
      //int my_nnofoc = 0; // total number of nodes for each face of this cell
      set<int> no_list; // noocv is not available...
      my_count[3] += 1 + faocv_i[icv+1]-faocv_i[icv] + bfocv_i[icv+1]-bfocv_i[icv];
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        const int nnof = nof_l - nof_f + 1;
        my_count[3] += nnof; // add face nodes
        // all nodes have to be dumped...
        for (int nof = nof_f; nof <= nof_l; ++nof) {
          const int ino = noofa_v[nof];
          no_list.insert(ino);
        }
      }
      for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
        const int ibf = bfocv_v[boc];
        const int nob_f = noobf_i[ibf];
        const int nob_l = noobf_i[ibf+1]-1;
        const int nnob = nob_l - nob_f + 1;
        my_count[3] += nnob; // add boundary face nodes
        // all nodes have to be dumped...
        for (int nob = nob_f; nob <= nob_l; ++nob) {
          const int ino = noobf_v[nob];
          no_list.insert(ino);
        }
      }
      const int nnoc = no_list.size();
      my_count[1]++; // this cell is getting dumped
      my_count[2] += nnoc; // add cell nodes
    }
  }

  int8 my_disp[4]; // nodal, element displacements..
  MPI_Scan(my_count,my_disp,4,MPI_INT8,MPI_SUM,mpi_comm);
  FOR_I4 my_disp[i] -= my_count[i];

  // put our global numbers into the no_flag
  FOR_INO if (no_flag[ino] >= 0)
    no_flag[ino] += my_disp[0];

  int8 count[4]; // everyone gets a copy of the node,element counts
  MPI_Allreduce(my_count,count,4,MPI_INT8,MPI_SUM,mpi_comm);

  // here we assume that flagged cvs are identified by (cv_flag[icv] != 0)...
  int myNumNodes = my_count[0];
  int myNumCells = my_count[1];
  int myNumCellNodes = my_count[2];
  int myNumCellFaNo = my_count[3];
  int8 myNodeDisp = my_disp[0];
  int8 myCellDisp = my_disp[1];
  int8 myCellNodeDisp = my_disp[2];
  int8 myCellFaDisp = my_disp[3];
  int8 NumNodes = count[0];
  int8 NumCells = count[1];
  int8 NumCellNodes = count[2];
  int8 NumCellFaNo = count[3];

  if (mpi_rank == 0)
    cout << " > global counts: " << NumNodes << " " << NumCells << " " << NumCellNodes
      << " " << double(NumCellNodes)/double(NumCells) << " " << NumCellFaNo << endl;

  // open the file...
  MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...
  MPI_File_delete((char*)filename.c_str(),MPI_INFO_NULL);

  MPI_File fh;
  int ierr = MPI_File_open(mpi_comm,(char*)filename.c_str(),MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);
  MPI_Check(ierr,"writeFlaggedCvsVTK::MPI_File_open");

  // the char buffer should be big enough for the largest part...
  MPI_Offset cbuf_size = 10000; //+max(3*4*myNumNodes,8*myNumCellNodes+17*myNumCells+8*myNumCellFaNo);
  char * cbuf = new char[cbuf_size];

  // Checking Endian-ness of the machine. The answer is in Endian[endian]
  const char *Endian[] = { "BigEndian", "LittleEndian" };
  const int endian = MiscUtils::getEndianness();

  // everyone prepare the header - it helps us all know the offset...
  MPI_Offset offset = 0;

  int n = sprintf(cbuf, "<?xml version=\"1.0\"?>\n");
  n += sprintf(cbuf+n, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"%s\" header_type=\"UInt64\">\n",Endian[endian]);
  n += sprintf(cbuf+n, "  <UnstructuredGrid>\n");
  n += sprintf(cbuf+n, "    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n", NumNodes, NumCells);

  n += sprintf(cbuf+n, "      <Points>\n");
  n += sprintf(cbuf+n, "        <DataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += int8_size + 3 * NumNodes * double_size;

  n += sprintf(cbuf+n, "      </Points>\n");

  n += sprintf(cbuf+n, "      <Cells>\n");
  n += sprintf(cbuf+n, "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += int8_size + NumCellNodes * int8_size;

  n += sprintf(cbuf+n, "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += int8_size + NumCells * int8_size;

  n += sprintf(cbuf+n, "        <DataArray type=\"Int8\" Name=\"types\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += int8_size + NumCells * sizeof(char);

  // face connectivity info --> this should be of size sum_poly(1 + poly_nFace + poly_nNodes)
  n += sprintf(cbuf+n, "        <DataArray type=\"Int64\" Name=\"faces\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += int8_size + NumCellFaNo * int8_size;

  // storage for the offsets of face connectivity array
  n += sprintf(cbuf+n, "        <DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += int8_size + NumCells * int8_size;

  n += sprintf(cbuf+n, "      </Cells>\n");

  n += sprintf(cbuf+n, "      <PointData>\n");

  string my_string;
  for (int ivar = 0; ivar < nscalars; ++ivar) {
    my_string = scalar_names[ivar];
    if (my_string == "x" || my_string == "y" || my_string == "z" || my_string == "X" || my_string == "Y" || my_string == "Z") {
      my_string += "-scalar";
    }
    // for some reason paraview doesn't like "<"
    my_string = MiscUtils::replaceAll(my_string,"<~","_le_");
    my_string = MiscUtils::replaceAll(my_string,"<","_lt_");
    n += sprintf(cbuf+n, "        <DataArray type=\"Float64\" Name=\"%s\" format=\"appended\" offset=\"%lld\" />\n",
        my_string.c_str(), (long long int) offset);
    offset += int8_size + NumNodes * double_size;
  }
  for (int ivar = 0; ivar < nvectors; ++ivar) {
    my_string = vector_names[ivar];
    my_string = MiscUtils::replaceAll(my_string,"<~","_le_");
    my_string = MiscUtils::replaceAll(my_string,"<","_lt_");
    n += sprintf(cbuf+n, "        <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lld\" />\n",
        my_string.c_str(), (long long int) offset);
    offset += int8_size + 3 * NumNodes * double_size;
  }

  n += sprintf(cbuf+n, "      </PointData>\n");

  n += sprintf(cbuf+n, "    </Piece>\n");

  n += sprintf(cbuf+n, "  </UnstructuredGrid>\n");

  n += sprintf(cbuf+n, "<AppendedData encoding=\"raw\">\n_");

  // the header...
  //if (mpi_rank == 0) {
  //  cout << cbuf << endl;
  //  cout << "end - header: " << offset << endl;
  //  cout << "n " << n << endl;
  //}

  int n_header = n;

  // rank 0 writes the header and the first bit of coord data. Ranks 1,2,etc
  // reset n and just write their coord data...

  offset = 0;

  if (mpi_rank == 0) {
    // write out header
    MPI_File_write(fh,cbuf,n_header,MPI_CHAR,MPI_STATUS_IGNORE);
  }
  delete[] cbuf; cbuf = NULL;

  // ==============================
  // x_no...
  // ==============================

  // now write header + x_no...
  offset += n_header; assert(sizeof(char)==1);

  double *dbuf3 = new double[myNumNodes*3];
  {
    int ind = 0;
    FOR_INO {
      if (no_flag[ino] >= 0) {
        FOR_I3 dbuf3[3*ind+i] = x_no[ino][i];
        ind++;
      }
    }
    assert(ind == myNumNodes);
  }

  int8 nbytes = 1*int8_size+NumNodes*3*double_size;
  if (mpi_rank == 0) MPI_File_write_at(fh,offset,&nbytes,1,MPI_INT8,MPI_STATUS_IGNORE);
  offset += int8_size*1;
  writeChunkedData(fh,offset+myNodeDisp*3*double_size,(double*)dbuf3,myNumNodes*3,mpi_comm);
  offset += double_size*3*NumNodes;
  if (nvectors == 0) delete[] dbuf3;

  // ==============================
  // connectivity...
  // this is in 5 parts: global indices, the
  // offsets, types, faces and face offets. It can all be done at once...
  // ==============================

  // recall that the 0-indexed ino_global is in no_flag in only the
  // nodes that actually OWN the nodes. otherwise, there is a -1. Same for
  // faces, so we need to have temporary

  // also update the global node counts to
  // overwrite the -1 values on processors that did not
  // write data...

  int no_minus1_check = 0;
  int no_minus2_check = 0;
  FOR_INO {
    if (no_flag[ino] >= 0) {
      assert((no_flag[ino] >= myNodeDisp)&&(no_flag[ino] < myNodeDisp+myNumNodes));
    }
    else if (no_flag[ino] == -1) {
      ++no_minus1_check;
    }
    else {
      assert(no_flag[ino] == -2);
      ++no_minus2_check;
    }
  }

  updateNoData(no_flag,MAX_NO_PERIODIC_DATA);

  // should have gotten rid of all -1's...
  FOR_INO assert(no_flag[ino] != -1);

  // start the offset at our node-of-cell disp...
  int8 disp = myCellNodeDisp;
  int8 fdisp = myCellFaDisp;
  int8 *int8_buf = new int8[myNumCellNodes];
  int8 *int8_buf_2 = new int8[myNumCells];
  assert(cbuf == NULL); cbuf = new char[myNumCells];
  int8 *int8_buf_3 = new int8[myNumCellFaNo];
  int8 *int8_buf_4 = new int8[myNumCells];
  int cv_index = 0, index = 0, index_3 = 0;
  FOR_ICV {
    if (cv_flag[icv] != 0) {
      set<int> no_list; // noocv is not available...
      int foc_f = faocv_i[icv];
      int foc_l = faocv_i[icv+1]-1;
      int nfoc = foc_l - foc_f + 1;
      for (int foc = foc_f; foc <= foc_l; ++foc) {
        const int ifa = faocv_v[foc];
        for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
          const int ino = noofa_v[nof];
          no_list.insert(ino);
        }
      }
      int boc_f = bfocv_i[icv];
      int boc_l = bfocv_i[icv+1]-1;
      int nboc = boc_l - boc_f + 1;
      for (int boc = boc_f; boc <= boc_l; ++boc) {
        const int ibf = bfocv_v[boc];
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int ino = noobf_v[nob];
          no_list.insert(ino);
        }
      }
      int nnoc = no_list.size();

      // node numbers...
      for (set<int>::iterator it = no_list.begin(); it != no_list.end(); ++it) {
        int8_buf[index++] = no_flag[*it];
      }

      // offset...
      disp += nnoc;
      int8_buf_2[cv_index] = disp;

      // type...
      cbuf[cv_index] = VTK_POLY;

      // face connectivity...
      // write number of faces for this cell
      int8_buf_3[index_3++] = nfoc+nboc; // int8 to keep it simple
      fdisp += 1;

      // loop over faces of cell
      for (int foc = foc_f; foc <= foc_l; foc++) {
        int ifa   = faocv_v[foc];
        int nof_f = noofa_i[ifa];
        int nof_l = noofa_i[ifa+1]-1;
        int nnof  = nof_l - nof_f + 1;

        // number of nodes of this face
        int8_buf_3[index_3++] = nnof;

        // loop over nodes of face
        for (int nof = nof_f ; nof <= nof_l; nof++) {
          int8_buf_3[index_3++] = no_flag[noofa_v[nof]];
        }
        fdisp += 1+nnof;
      }
      for (int boc = boc_f; boc <= boc_l; ++boc) {
        const int ibf = bfocv_v[boc];
        int nob_f = noobf_i[ibf];
        int nob_l = noobf_i[ibf+1]-1;
        int nnob  = nob_l - nob_f + 1;

        // number of nodes of this face
        int8_buf_3[index_3++] = nnob;

        // loop over nodes of face
        for (int nob = nob_f; nob <= nob_l; nob++) {
          int8_buf_3[index_3++] = no_flag[noobf_v[nob]];
        }
        fdisp += 1+nnob;
      }

      // face offset...
      int8_buf_4[cv_index] = fdisp;
      cv_index++;
    }
  }
  assert(cv_index == myNumCells);

  // connectivity part...
  nbytes = (1+NumCellNodes)*int8_size;
  if (mpi_rank == 0) MPI_File_write_at(fh,offset,&nbytes,1,MPI_INT8,MPI_STATUS_IGNORE);
  offset += int8_size;
  writeChunkedData<int8>(fh,offset+myCellNodeDisp*int8_size,int8_buf,myNumCellNodes,mpi_comm);
  offset += int8_size*NumCellNodes;
  delete[] int8_buf;

  // offset part...
  nbytes = (1+NumCells)*int8_size;
  if (mpi_rank == 0) MPI_File_write_at(fh,offset,&nbytes,1,MPI_INT8,MPI_STATUS_IGNORE);
  offset += int8_size;
  writeChunkedData<int8>(fh,offset+myCellDisp*int8_size,int8_buf_2,myNumCells,mpi_comm);
  offset += int8_size*NumCells;
  delete[] int8_buf_2;

  // type part...
  nbytes = 1*int8_size+NumCells*1;
  if (mpi_rank == 0) MPI_File_write_at(fh,offset,&nbytes,1,MPI_INT8,MPI_STATUS_IGNORE);
  offset += int8_size;
  writeChunkedData<char>(fh,offset+myCellDisp,cbuf,myNumCells,mpi_comm);
  offset += NumCells;
  delete[] cbuf; cbuf = NULL;

  // face part...
  nbytes = (1+NumCellFaNo)*int8_size;
  if (mpi_rank == 0) MPI_File_write_at(fh,offset,&nbytes,1,MPI_INT8,MPI_STATUS_IGNORE);
  offset += int8_size;
  writeChunkedData<int8>(fh,offset+myCellFaDisp*int8_size,int8_buf_3,myNumCellFaNo,mpi_comm);
  offset += int8_size*NumCellFaNo;
  delete[] int8_buf_3;

  // face offset part...
  nbytes = (1+NumCells)*int8_size;
  if (mpi_rank == 0) MPI_File_write_at(fh,offset,&nbytes,1,MPI_INT8,MPI_STATUS_IGNORE);
  offset += int8_size;
  writeChunkedData<int8>(fh,offset+myCellDisp*int8_size,int8_buf_4,myNumCells,mpi_comm);
  offset += int8_size*NumCells;
  delete[] int8_buf_4;

  // now return no_flag to only flagged nodes...

  int no_minus1_check2 = 0;
  int no_minus2_check2 = 0;
  FOR_INO {
    if (no_flag[ino] >= 0) {
      if ((no_flag[ino] < myNodeDisp)||(no_flag[ino] >= myNodeDisp+myNumNodes)) {
        no_flag[ino] = -1;
        ++no_minus1_check2;
      }
    }
    else {
      assert(no_flag[ino] == -2);
      ++no_minus2_check2;
    }
  }
  assert(no_minus1_check2 == no_minus1_check);
  assert(no_minus2_check2 == no_minus2_check);

  //if (mpi_rank == 0)
  //  cout << "after connectivity: offset - n_header: " << offset-n_header << endl;

  // scalars...
  if (nscalars > 0) {
    double *no_data = new double[nno];
    double *dbuf = new double[myNumNodes];
    for (int ivar = 0; ivar < nscalars; ++ivar) {
      setNoDN(no_data,scalar_names[ivar]);

      int ind = 0;
      FOR_INO {
        if (no_flag[ino] >= 0) {
          dbuf[ind++] = no_data[ino];
        }
      }
      assert(ind == myNumNodes);

      nbytes = 1*int8_size+NumNodes*double_size;
      if (mpi_rank == 0) MPI_File_write_at(fh,offset,&nbytes,1,MPI_INT8,MPI_STATUS_IGNORE);
      offset += int8_size;
      MPI_File_write_at(fh,offset+myNodeDisp*double_size,dbuf,myNumNodes,MPI_DOUBLE,MPI_STATUS_IGNORE);
      offset += double_size*NumNodes;

    }
    delete[] no_data;
    delete[] dbuf;
  }

  // vectors...
  if (nvectors > 0) {
    double (*no_data)[3] = new double[nno][3];
    for (int ivar = 0; ivar < nvectors; ++ivar) {
      setNoDN3(no_data,vector_names[ivar]);

      int ind = 0;
      FOR_INO {
        if (no_flag[ino] >= 0) {
          FOR_I3 dbuf3[ind*3+i] = no_data[ino][i];
          ind++;
        }
      }
      assert(ind == myNumNodes);

      nbytes = 1*int8_size+NumNodes*3*double_size;
      if (mpi_rank == 0) MPI_File_write_at(fh,offset,&nbytes,1,MPI_INT8,MPI_STATUS_IGNORE);
      offset += int8_size;
      MPI_File_write_at(fh,offset+myNodeDisp*3*double_size,dbuf3,myNumNodes*3,MPI_DOUBLE,MPI_STATUS_IGNORE);
      offset += double_size*3*NumNodes;
      //if (mpi_rank == 0)
      //  cout << "after no vector: " << offset << ", rel to end of header: " << offset-n_header << endl;
    }
    delete[] no_data;
    delete[] dbuf3;
  }
  delete[] no_flag;

  // ==========================================================================
  // end of file: just rank0 writes...
  // ==========================================================================

  assert(cbuf == NULL); cbuf = new char[256];
  n = sprintf(cbuf, "\n  </AppendedData>\n</VTKFile>\n");
  if (mpi_rank == 0) MPI_File_write_at(fh,offset,cbuf,n,MPI_CHAR,MPI_STATUS_IGNORE);
  delete[] cbuf;

  offset += n;
  //if (mpi_rank == 0)
  //  cout << "end size - header:" << offset-n_header << endl;

  MPI_File_set_size(fh,offset);
  MPI_File_close(&fh);

}

void StaticSolver::writeFlaggedBfZonesVTK(const string& filename,bool* bfzone_flag,
					  const vector<string>& scalar_names_cv,const int nscalars_cv,
					  const vector<string>& vector_names_cv,const int nvectors_cv,
					  const vector<string>& scalar_names_bf,const vector<string>& scalar_zones_bf,const int nscalars_bf,
					  const vector<string>& vector_names_bf,const vector<string>& vector_zones_bf,const int nvectors_bf,const int nzones_bf) {

  COUT1("StaticSolver::writeFlaggedBfZonesVTK() " + filename);

  int* no_flag = new int[nno];
  int* bf_flag = new int[nbf];

  FOR_INO no_flag[ino] = mpi_size;
  FOR_IBF bf_flag[ibf] = mpi_size;

  for (map<const string,int>::const_iterator it = bfZoneNameMap.begin(); it != bfZoneNameMap.end(); ++it) {
    const int iz = it->second; // zone index
    if (bfzone_flag[iz]) {
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int ino = noobf_v[nob];
          assert((ino >= 0)&&(ino < nno));
          no_flag[ino] = mpi_rank;
        }
        int nnof = noobf_i[ibf+1]-noobf_i[ibf];
        assert( nnof >= 3 ); //check later ...
        bf_flag[ibf] = mpi_rank;
      }
    }
  }

  int my_count[3] = {0,0,0}; // local ncount,ecount and necount

  updateNoData(no_flag,MIN_NO_PERIODIC_DATA); // smallest rank wins

  FOR_INO {
    if (no_flag[ino] == mpi_rank)
      no_flag[ino] = my_count[0]++;     // this node is getting dumped
    else if (no_flag[ino] < mpi_size)
      no_flag[ino] = -1;                // this node gets dumped, but not by us
    else {
      assert(no_flag[ino] == mpi_size);
      no_flag[ino] = -2;                // nobody owns this node
    }
  }

  // need to tiebreak faces, but we need to go through the cells
  // since we don't have a face communicator at this stage.
  //updateBfData(bf_flag,MIN_NO_PERIODIC_DATA);

  int my_max_nnof = 0;
  FOR_IBF {
    if (bf_flag[ibf] == mpi_rank) {
      bf_flag[ibf] = -1; // this but needs to be dumped ...
      my_count[1] += 1;
      const int nnof = noobf_i[ibf+1]-noobf_i[ibf];
      my_count[2] += nnof;
      my_max_nnof = max(my_max_nnof,nnof);
    }
    else {
      assert(bf_flag[ibf] == mpi_size); // it will fail if smaller rank wins ...
      bf_flag[ibf] = -2; // this face is not dumped
    }
  }

  // figure out our nodal and element displacement...

  int my_disp[3];
  MPI_Scan(my_count,my_disp,3,MPI_INT,MPI_SUM,mpi_comm);
  my_disp[0] -= my_count[0];
  my_disp[1] -= my_count[1];
  my_disp[2] -= my_count[2];

  for (int ino = 0; ino < nno; ++ino)
    if (no_flag[ino] >= 0)
      no_flag[ino] += my_disp[0];

  // everyone gets the counts...

  int count[3]={0,0,0};
  MPI_Allreduce(my_count,count,3,MPI_INT,MPI_SUM,mpi_comm);

  int myNumNodes = my_count[0];
  int myNumCells = my_count[1];
  int myNumCellNodes = my_count[2];
  int myNodeDisp = my_disp[0];
  int myCellDisp = my_disp[1];
  int myCellNodeDisp = my_disp[2];
  int NumNodes = count[0];
  int NumCells = count[1];
  int NumCellNodes = count[2];

  // open the file...
  MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...
  MPI_File_delete((char*)filename.c_str(),MPI_INFO_NULL);

  MPI_File fh;
  int ierr = MPI_File_open(mpi_comm,(char*)filename.c_str(),MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);
  MPI_Check(ierr,"writeFlaggedBfZonesVTK::MPI_File_open");

  // the char buffer should be big enough for the largest part...
  MPI_Offset cbuf_size = 10000 + max(4*3*myNumNodes,max(4*3*myNumCells,4*myNumCellNodes+9*myNumCells));
  char * cbuf = new char[cbuf_size];

  // Checking Endian-ness of the machine. The answer is in Endian[endian]
  const char *Endian[] = { "BigEndian", "LittleEndian" };
  const int endian = MiscUtils::getEndianness();

  // everyone prepare the header - it helps us all know the offset...
  MPI_Offset offset = 0;

  int n = sprintf(cbuf, "<?xml version=\"1.0\"?>\n");
  n += sprintf(cbuf+n, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n",Endian[endian]);
  n += sprintf(cbuf+n, "  <UnstructuredGrid>\n");
  n += sprintf(cbuf+n, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", NumNodes, NumCells);

  n += sprintf(cbuf+n, "      <Points>\n");
  n += sprintf(cbuf+n, "        <DataArray type=\"Float32\" Name=\"coordinates\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += sizeof(int) + 3 * NumNodes * sizeof(float);

  n += sprintf(cbuf+n, "      </Points>\n");

  n += sprintf(cbuf+n, "      <Cells>\n");
  n += sprintf(cbuf+n, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += sizeof(int) + NumCellNodes * sizeof(int);

  n += sprintf(cbuf+n, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += sizeof(int) + NumCells * sizeof(int);

  n += sprintf(cbuf+n, "        <DataArray type=\"Int8\" Name=\"types\" format=\"appended\" offset=\"%lld\" />\n", (long long int) offset);

  offset += sizeof(int) + NumCells * sizeof(char);

  n += sprintf(cbuf+n, "      </Cells>\n");

  n += sprintf(cbuf+n, "      <PointData>\n");

  string my_string;
  for (int ivar = 0; ivar < nscalars_cv; ++ivar) {
    my_string = scalar_names_cv[ivar];
    if (my_string == "x" || my_string == "y" || my_string == "z" || my_string == "X" || my_string == "Y" || my_string == "Z") {
      my_string += "-scalar";
    }
    my_string = MiscUtils::replaceAll(my_string,"<~","_le_");
    my_string = MiscUtils::replaceAll(my_string,"<","_lt_");
    n += sprintf(cbuf+n, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"appended\" offset=\"%lld\" />\n",
        my_string.c_str(), (long long int) offset);
    offset += sizeof(int) + NumNodes * sizeof(float);
  }
  if (nscalars_bf > 0) {
    assert(nzones_bf > 0);
    assert(nscalars_bf%nzones_bf == 0);
    for (int ivar = 0; ivar < nscalars_bf/nzones_bf; ++ivar) {
      my_string = scalar_names_bf[ivar];
      size_t found = my_string.find(scalar_zones_bf[ivar]);
      assert(found != std::string::npos);
      my_string = MiscUtils::replaceAll(my_string,scalar_zones_bf[ivar]+":","");
      if (my_string == "x" || my_string == "y" || my_string == "z" || my_string == "X" || my_string == "Y" || my_string == "Z") {
        my_string += "-scalar";
      }
      my_string = MiscUtils::replaceAll(my_string,"<~","_le_");
      my_string = MiscUtils::replaceAll(my_string,"<","_lt_");
      n += sprintf(cbuf+n, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"appended\" offset=\"%lld\" />\n",
          my_string.c_str(), (long long int) offset);
      offset += sizeof(int) + NumNodes * sizeof(float);
    }
  }
  for (int ivar = 0; ivar < nvectors_cv; ++ivar) {
    my_string = vector_names_cv[ivar];
    my_string = MiscUtils::replaceAll(my_string,"<~","_le_");
    my_string = MiscUtils::replaceAll(my_string,"<","_lt_");
    n += sprintf(cbuf+n, "        <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lld\" />\n",
        my_string.c_str(), (long long int) offset);
    offset += sizeof(int) + 3 * NumNodes * sizeof(float);
  }
  if (nvectors_bf > 0) {
    assert(nzones_bf > 0);
    assert(nvectors_bf%nzones_bf == 0);
    for (int ivar = 0; ivar < nvectors_bf/nzones_bf; ++ivar) {
      my_string = vector_names_bf[ivar];
      size_t found = my_string.find(vector_zones_bf[ivar]);
      assert(found != std::string::npos);
      my_string = MiscUtils::replaceAll(my_string,vector_zones_bf[ivar]+":","");
      my_string = MiscUtils::replaceAll(my_string,"<~","_le_");
      my_string = MiscUtils::replaceAll(my_string,"<","_lt_");
      n += sprintf(cbuf+n, "        <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lld\" />\n",
          my_string.c_str(), (long long int) offset);
      offset += sizeof(int) + 3 * NumNodes * sizeof(float);
    }
  }

  n += sprintf(cbuf+n, "      </PointData>\n");

  n += sprintf(cbuf+n, "    </Piece>\n");

  n += sprintf(cbuf+n, "  </UnstructuredGrid>\n");

  n += sprintf(cbuf+n, "<AppendedData encoding=\"raw\">\n_");

  // the header...
  //if (mpi_rank == 0)
  //  cout << cbuf << endl;

  int n_header = n;

  // rank 0 writes the header and the first bit of coord data. Ranks 1,2,etc
  // reset n and just write their coord data...

  if (mpi_rank != 0)
    n = 0;

  // ==============================
  // x_no...
  // ==============================

  if (mpi_rank == 0) {
    int bytes = sizeof(int) + 3 * NumNodes * sizeof(float);
    memcpy(cbuf+n,&bytes,sizeof(int));
    n += sizeof(int);
  }

  for (int ino = 0; ino < nno; ++ino) {
    if (no_flag[ino] >= 0) {
      float xyz[3];
      FOR_I3 xyz[i] = float(x_no[ino][i]);
      memcpy(cbuf+n,xyz,sizeof(float)*3);
      n += sizeof(float)*3;
    }
  }

  // should be large enough...
  assert( n <= cbuf_size );

  // now write header + x_no...
  offset = 0;

  MPI_Aint dispA;
  if (mpi_rank == 0) {
    dispA = 0;
    assert(n == n_header + sizeof(int) + sizeof(float)*myNumNodes*3);
  }
  else {
    dispA = n_header + sizeof(int) + sizeof(float)*myNodeDisp*3;
    assert(n == sizeof(float)*myNumNodes*3);
  }

  writeChunkedData<char>(fh,dispA,cbuf,n,mpi_comm);
  offset += n_header + sizeof(int) + sizeof(float)*3*NumNodes;

  // ==============================
  // connectivity...
  // this is in 3 parts: global indices, the
  // offsets, and the types. It can all be done at once...
  // ==============================

  {

    n = 0;
    int n1 = sizeof(int)*myNumCellNodes;
    if (mpi_rank == 0) n1 += sizeof(int);
    int n2 = n1 + sizeof(int)*myNumCells;
    if (mpi_rank == 0) n2 += sizeof(int);
    int nf3 = n2 + sizeof(char)*myNumCells;
    if (mpi_rank == 0) nf3 += sizeof(int);

    if (mpi_rank == 0) {

      // connectivity part...
      int bytes = sizeof(int) + NumCellNodes * sizeof(int);
      memcpy(cbuf+n,&bytes,sizeof(int));
      n += sizeof(int);

      // offset part...
      bytes = sizeof(int) + NumCells * sizeof(int);
      memcpy(cbuf+n1,&bytes,sizeof(int));
      n1 += sizeof(int);

      // type part...
      bytes = sizeof(int) + NumCells * sizeof(char);
      memcpy(cbuf+n2,&bytes,sizeof(int));
      n2 += sizeof(int);

    }

    // also update the global node and face counts to
    // overwrite the -1 values on processors that did not
    // write data...

    int no_minus1_check = 0;
    int no_minus2_check = 0;
    FOR_INO {
      if (no_flag[ino] >= 0) {
        assert((no_flag[ino] >= myNodeDisp)&&(no_flag[ino] < myNodeDisp+myNumNodes));
      }
      else if (no_flag[ino] == -1) {
        ++no_minus1_check;
      }
      else {
        assert(no_flag[ino] == -2);
        ++no_minus2_check;
      }
    }

    updateNoData(no_flag,MAX_NO_PERIODIC_DATA);

    // start the offset at our node-of-cell disp...
    int disp = myCellNodeDisp;
    int * ibuf = new int[my_max_nnof];
    FOR_IBF {
      if (bf_flag[ibf] == -1) {

        int nnof = noobf_i[ibf+1]-noobf_i[ibf];

        // connectivity...
        for (int nof = noobf_i[ibf]; nof != noobf_i[ibf+1]; ++nof) {
          ibuf[nof-noobf_i[ibf]] = no_flag[noobf_v[nof]];
        }
        memcpy(cbuf+n,ibuf,sizeof(int)*nnof);
        n += sizeof(int)*nnof;

        // offset...
        disp += nnof;
        memcpy(cbuf+n1,&disp,sizeof(int));
        n1 += sizeof(int);

        // type...
        switch (nnof) {
          case 3:
            cbuf[n2++] = VTK_TRIANGLE;
            break;
          case 4:
            cbuf[n2++] = VTK_QUAD;
            break;
          default:
            cbuf[n2++] = VTK_POLYGON;
        }
      }
    }
    delete[] ibuf;

    assert(n2 <= cbuf_size);

    int n3[3];
    MPI_Aint disp3[3];
    n3[0] = sizeof(int)*myNumCellNodes;
    if (mpi_rank == 0) {
      n3[0] += sizeof(int);
      assert(myCellNodeDisp == 0);
      disp3[0] = 0;
    }
    else {
      disp3[0] = sizeof(int)*(1+myCellNodeDisp);
    }
    n3[1] = sizeof(int)*myNumCells;
    if (mpi_rank == 0) {
      n3[1] += sizeof(int);
      assert(myCellDisp == 0);
      disp3[1] = sizeof(int)*(1+NumCellNodes);
    }
    else {
      disp3[1] = sizeof(int)*(1+NumCellNodes+1+myCellDisp);
    }
    n3[2] = sizeof(char)*myNumCells;
    if (mpi_rank == 0) {
      n3[2] += sizeof(int);
      assert(myCellDisp == 0);
      disp3[2] = sizeof(int)*(1+NumCellNodes+1+NumCells);
    }
    else {
      disp3[2] = sizeof(int)*(1+NumCellNodes+1+NumCells+1) + sizeof(char)*myCellDisp;
    }

    writeChunkedData<char>(fh,offset+disp3[0],cbuf,n3[0],mpi_comm);
    writeChunkedData<char>(fh,offset+disp3[1],cbuf+n3[0],n3[1],mpi_comm);
    writeChunkedData<char>(fh,offset+disp3[2],cbuf+n3[0]+n3[1],n3[2],mpi_comm);

    // now return no_flag to only flagged nodes...

    int no_minus1_check2 = 0;
    int no_minus2_check2 = 0;
    FOR_INO {
      if (no_flag[ino] >= 0) {
        if ((no_flag[ino] < myNodeDisp)||(no_flag[ino] >= myNodeDisp+myNumNodes)) {
          no_flag[ino] = -1;
          ++no_minus1_check2;
        }
      }
      else {
        assert(no_flag[ino] == -2);
        ++no_minus2_check2;
      }
    }
    assert(no_minus1_check2 == no_minus1_check);
    assert(no_minus2_check2 == no_minus2_check);
  }

  offset += sizeof(int)*(1+NumCellNodes+1+NumCells+1) + sizeof(char)*NumCells;

  // we all write data at the nodes...

  // scalars...
  if (nscalars_cv > 0) {
    if (mpi_rank == 0) {
      assert(myNodeDisp == 0);
      dispA = 0;
    }
    else {
      dispA = sizeof(int) + sizeof(float)*myNodeDisp;
    }
    double *no_data = new double[nno];
    for (int ivar = 0; ivar < nscalars_cv; ++ivar) {
      setNoDN(no_data,scalar_names_cv[ivar]);

      n = 0;
      if (mpi_rank == 0) {
        int bytes = sizeof(int) + NumNodes * sizeof(float);
        memcpy(cbuf+n,&bytes,sizeof(int));
        n += sizeof(int);
      }

      int index = 0;
      for (int ino = 0; ino < nno; ++ino) {
        if (no_flag[ino] >= 0) {
          float phi = float(no_data[ino]);
          memcpy(cbuf+n,&phi,sizeof(float));
          n += sizeof(float);
          index++;
        }
      }
      assert(index == myNumNodes);

      writeChunkedData<char>(fh,offset+dispA,cbuf,n,mpi_comm);
      offset += sizeof(int)+sizeof(float)*NumNodes;

    }
    delete[] no_data;
  }

  if (nscalars_bf > 0) {
    if (mpi_rank == 0) {
      assert(myNodeDisp == 0);
      dispA = 0;
    }
    else {
      dispA = sizeof(int) + sizeof(float)*myNodeDisp;
    }
    double* no_data = new double[nno];
    vector<string> var_subvec_bf(nzones_bf);
    for (int ivar = 0; ivar < nscalars_bf/nzones_bf; ++ivar) {
      for (int izone = 0; izone < nzones_bf; ++izone) {
        var_subvec_bf[izone] = scalar_names_bf[izone*(nscalars_bf/nzones_bf)+ivar];
      }
      setBfDN(no_data,var_subvec_bf);

      n = 0;
      if (mpi_rank == 0) {
        int bytes = sizeof(int) + NumNodes * sizeof(float);
        memcpy(cbuf+n,&bytes,sizeof(int));
        n += sizeof(int);
      }

      int index = 0;
      for (int ino = 0; ino < nno; ++ino) {
        if (no_flag[ino] >= 0) {
          float phi = float(no_data[ino]);
          memcpy(cbuf+n,&phi,sizeof(float));
          n += sizeof(float);
          index++;
        }
      }
      assert(index == myNumNodes);

      writeChunkedData<char>(fh,offset+dispA,cbuf,n,mpi_comm);
      offset += sizeof(int)+sizeof(float)*NumNodes;

    }
    delete[] no_data;
  }

  // vectors...
  if (nvectors_cv > 0) {
    if (mpi_rank == 0) {
      assert(myNodeDisp == 0);
      dispA = 0;
    }
    else {
      dispA = sizeof(int) + sizeof(float)*myNodeDisp*3;
    }
    double (*no_data)[3] = new double[nno][3];
    for (int ivar = 0; ivar < nvectors_cv; ++ivar) {
      setNoDN3(no_data,vector_names_cv[ivar]);

      n = 0;
      if (mpi_rank == 0) {
        int bytes = sizeof(int) + NumNodes * sizeof(float) * 3;
        memcpy(cbuf+n,&bytes,sizeof(int));
        n += sizeof(int);
      }

      int index = 0;
      for (int ino = 0; ino < nno; ++ino) {
        if (no_flag[ino] >= 0) {
          float vec[3]; FOR_I3 vec[i] = float(no_data[ino][i]);
          memcpy(cbuf+n,vec,sizeof(float)*3);
          n += sizeof(float)*3;
          index++;
        }
      }
      assert(index == myNumNodes);

      writeChunkedData<char>(fh,offset+dispA,cbuf,n,mpi_comm);
      offset += sizeof(int)+sizeof(float)*NumNodes*3;

    }
    delete[] no_data;
  }

  if (nvectors_bf > 0) {
    if (mpi_rank == 0) {
      assert(myNodeDisp == 0);
      dispA = 0;
    }
    else {
      dispA = sizeof(int) + sizeof(float)*myNodeDisp*3;
    }
    double (*no_data)[3] = new double[nno][3];
    vector<string> var_subvec_bf(nzones_bf);
    for (int ivar = 0; ivar < nvectors_bf/nzones_bf; ++ivar) {
      for (int izone = 0; izone < nzones_bf; ++izone) {
        var_subvec_bf[izone] = vector_names_bf[izone*(nvectors_bf/nzones_bf)+ivar];
      }
      setBfDN3(no_data,var_subvec_bf);

      n = 0;
      if (mpi_rank == 0) {
        int bytes = sizeof(int) + NumNodes * sizeof(float) * 3;
        memcpy(cbuf+n,&bytes,sizeof(int));
        n += sizeof(int);
      }

      int index = 0;
      for (int ino = 0; ino < nno; ++ino) {
        if (no_flag[ino] >= 0) {
          float vec[3]; FOR_I3 vec[i] = float(no_data[ino][i]);
          memcpy(cbuf+n,vec,sizeof(float)*3);
          n += sizeof(float)*3;
          index++;
        }
      }
      assert(index == myNumNodes);

      writeChunkedData<char>(fh,offset+dispA,cbuf,n,mpi_comm);
      offset += sizeof(int)+sizeof(float)*NumNodes*3;

    }
    delete[] no_data;
  }
  delete[] no_flag;
  delete[] bf_flag;

  // ==========================================================================
  // end of file: just rank0 writes...
  // ==========================================================================

  n = sprintf(cbuf, "\n  </AppendedData>\n</VTKFile>\n");

  if (mpi_rank == 0) {
    MPI_File_write_at(fh,offset,cbuf,n,MPI_CHAR,MPI_STATUS_IGNORE);
  }

  offset += n;

  MPI_File_set_size(fh,offset);
  MPI_File_close(&fh);

  delete[] cbuf;
}

#undef VTK_TRIANGLE
#undef VTK_QUAD
#undef VTK_POLYGON
#undef VTK_TET
#undef VTK_HEX
#undef VTK_PRISM
#undef VTK_PYRAMID
#undef VTK_POLY

void StaticSolver::updateEnsightCaseFile(DataWriter* dw,const int step,const double time) {
  COUT1("StaticSolver::updateEnsightCaseFile()");

  if (mpi_rank == 0) {
    string filename = dw->name + ".case";
    MiscUtils::mkdir_for_file(filename.c_str());

    //if the name contains a file path,
    //strip out the prefix for the case control file

    size_t pos_path = dw->name.find_last_of("/\\");
    if (pos_path == string::npos) pos_path = 0;
    else ++pos_path;
    string prefix = dw->name.substr(pos_path);

    FILE * pFile;
    if (dw->ensight_casefile_loc == -1) {
      pFile = fopen ( filename.c_str() , "w" );

      // create FORMAT section...
      fprintf(pFile,"FORMAT\n");
      fprintf(pFile,"type: ensight gold\n\n");

      // create GEOMETRY section...
      fprintf(pFile,"GEOMETRY\n");
      if (dw->lp_index == -1)
        fprintf(pFile,"model: %s.geo\n\n",prefix.c_str()); // model stationary
      else
        fprintf(pFile,"model: 1 %s.********.geo\n\n",prefix.c_str()); // model transient

      // create VARIABLE section... nodal variables for now (need to look into cell vars to see
      fprintf(pFile,"VARIABLE\n");
      // can also store constants (that can chagned in time)...
      // scalars...
      for (int i = 0; i < dw->nscalars_cv; ++i) {
        // get rid of reserved chars (TODO make faster)...
        string name = MiscUtils::replaceAll(dw->scalar_names_cv[i],"(","<");
        name = MiscUtils::replaceAll(name,")",">");
        name = MiscUtils::replaceAll(name,"+","_p_");
        name = MiscUtils::replaceAll(name,"-","_m_");
        name = MiscUtils::replaceAll(name,"*","_t_");
        name = MiscUtils::replaceAll(name,"/","_d_");
        name = MiscUtils::replaceAll(name," ","");
        name = name.substr(0,19); // limited to 19 characters
        // the stars are used to hold the increment (time step for now)
        // the latter two digits are for storing the scalar index
        fprintf(pFile,"scalar per node: 1 %s %s.********.sca%02d\n",name.c_str(),prefix.c_str(),i);
        dw->scalar_descs_cv.push_back(name); // store so we don't need to do this again
      }
      if (dw->nscalars_bf > 0) {
        for (int i = 0; i < dw->nscalars_bf/dw->nzones_bf; ++i) {
          // get rid of reserved chars (TODO make faster)...
          string name = MiscUtils::replaceAll(dw->scalar_names_bf[i],"(","<");
          name = MiscUtils::replaceAll(name,dw->scalar_zones_bf[i]+":","");
          name = MiscUtils::replaceAll(name,")",">");
          name = MiscUtils::replaceAll(name,"+","_p_");
          name = MiscUtils::replaceAll(name,"-","_m_");
          name = MiscUtils::replaceAll(name,"*","_t_");
          name = MiscUtils::replaceAll(name,"/","_d_");
          name = MiscUtils::replaceAll(name," ","");
          name = name.substr(0,19); // limited to 19 characters
          // the stars are used to hold the increment (time step for now)
          // the latter two digits are for storing the scalar index
          fprintf(pFile,"scalar per node: 1 %s %s.********.sca%02d\n",name.c_str(),prefix.c_str(),i+dw->nscalars_cv);
          dw->scalar_descs_bf.push_back(name); // store so we don't need to do this again
        }
      }
      for (int i = 0; i < dw->nscalars_lp; ++i) {
        // get rid of reserved chars (TODO make faster)...
        string name = MiscUtils::replaceAll(dw->scalar_names_lp[i],"(","<");
        name = MiscUtils::replaceAll(name,")",">");
        name = MiscUtils::replaceAll(name,"+","_p_");
        name = MiscUtils::replaceAll(name,"-","_m_");
        name = MiscUtils::replaceAll(name,"*","_t_");
        name = MiscUtils::replaceAll(name,"/","_d_");
        name = MiscUtils::replaceAll(name," ","");
        name = name.substr(0,19); // limited to 19 characters
        // the stars are used to hold the increment (time step for now)
        // the latter two digits are for storing the scalar index
        fprintf(pFile,"scalar per node: 1 %s %s.********.sca%02d\n",name.c_str(),prefix.c_str(),i);
        dw->scalar_descs_lp.push_back(name); // store so we don't need to do this again
      }
      // vectors...
      for (int i = 0; i < dw->nvectors_cv; ++i) {
        // get rid of reserved chars (TODO make faster)...
        string name = MiscUtils::replaceAll(dw->vector_names_cv[i],"(","<");
        name = MiscUtils::replaceAll(name,")",">");
        name = MiscUtils::replaceAll(name,"+","_p_");
        name = MiscUtils::replaceAll(name,"-","_m_");
        name = MiscUtils::replaceAll(name,"*","_t_");
        name = MiscUtils::replaceAll(name,"/","_d_");
        name = MiscUtils::replaceAll(name," ","");
        // the stars are used to hold the increment (time step for now)
        // the latter two digits are for storing the vector index
        fprintf(pFile,"vector per node: 1 %s %s.********.vec%02d\n",name.c_str(),prefix.c_str(),i);
        dw->vector_descs_cv.push_back(name); // store so we don't need to do this again
      }
      if (dw->nvectors_bf > 0) {
        for (int i = 0; i < dw->nvectors_bf/dw->nzones_bf; ++i) {
          // get rid of reserved chars (TODO make faster)...
          string name = MiscUtils::replaceAll(dw->vector_names_bf[i],"(","<");
          name = MiscUtils::replaceAll(name,dw->vector_zones_bf[i]+":","");
          name = MiscUtils::replaceAll(name,")",">");
          name = MiscUtils::replaceAll(name,"+","_p_");
          name = MiscUtils::replaceAll(name,"-","_m_");
          name = MiscUtils::replaceAll(name,"*","_t_");
          name = MiscUtils::replaceAll(name,"/","_d_");
          name = MiscUtils::replaceAll(name," ","");
          // the stars are used to hold the increment (time step for now)
          // the latter two digits are for storing the vector index
          fprintf(pFile,"vector per node: 1 %s %s.********.vec%02d\n",name.c_str(),prefix.c_str(),i+dw->nvectors_cv);
          dw->vector_descs_bf.push_back(name); // store so we don't need to do this again
        }
      }
      for (int i = 0; i < dw->nvectors_lp; ++i) {
        // get rid of reserved chars (TODO make faster)...
        string name = MiscUtils::replaceAll(dw->vector_names_lp[i],"(","<");
        name = MiscUtils::replaceAll(name,")",">");
        name = MiscUtils::replaceAll(name,"+","_p_");
        name = MiscUtils::replaceAll(name,"-","_m_");
        name = MiscUtils::replaceAll(name,"*","_t_");
        name = MiscUtils::replaceAll(name,"/","_d_");
        name = MiscUtils::replaceAll(name," ","");
        // the stars are used to hold the increment (time step for now)
        // the latter two digits are for storing the vector index
        fprintf(pFile,"vector per node: 1 %s %s.********.vec%02d\n",name.c_str(),prefix.c_str(),i);
        dw->vector_descs_lp.push_back(name); // store so we don't need to do this again
      }
      fprintf(pFile,"\n");

      // create TIME section...
      fprintf(pFile,"TIME\n");
      fprintf(pFile,"time set: 1\n"); // only one time series for now
      fprintf(pFile,"number of steps: ");
      dw->ensight_casefile_loc = ftell(pFile);
      fprintf(pFile,"%8d",++dw->ensight_casefile_nsteps);
      fprintf(pFile,"\nfilename start number: %8d\n",step);
      fprintf(pFile,"filename increment: %8d\n",dw->interval);
      fprintf(pFile,"time values: %18.15le\n",time);

      fclose(pFile);
    }
    else {
      pFile = fopen ( filename.c_str() , "r+" );

      // update TIME section...
      fseek(pFile,dw->ensight_casefile_loc,SEEK_SET);
      fprintf(pFile,"%8d",++dw->ensight_casefile_nsteps);
      fseek(pFile,0,SEEK_END);
      fprintf(pFile,"            %18.15e\n",time);

      fclose(pFile);
    }
  }
  else {
    ++dw->ensight_casefile_nsteps; // used to determine first vs latter writes
  }
}

void StaticSolver::writeCvsAndFlaggedBfZonesDualEnsightMultiPart(DataWriter* dw,const int step,const double time) {

  // update case file with info about this step...
  updateEnsightCaseFile(dw,step,time);

  // here we assume that flagged cvs are identified by (cv_flag[icv] != 0)...

  int8* cv_flag = new int8[ncv];
  FOR_ICV cv_flag[icv] = 1; // dw->cv_flag[icv];
  FOR_IBF cv_flag[cvobf[ibf]] = 0; // internal only
  int8* cv_flag_d = new int8[ncv_d-ncv];
  updateCvdDataSeparateGhosts(cv_flag,cv_flag_d);

  int8 my_count[2] = {0,0}; // flagged cvs (including dual ghosts), and flagged tets
  FOR_ICV {
    if (cv_flag[icv] != 0)
      cv_flag[icv] = my_count[0]++;
    else
      cv_flag[icv] = -1;
  }
  for (int icv = ncv; icv < ncv_d; ++icv) {
    if (cv_flag_d[icv-ncv] != 0)
      cv_flag_d[icv-ncv] = my_count[0]++;
    else
      cv_flag_d[icv-ncv] = -1;
  }

  const int ndt = dualTetVec.size();
  int* dt_flag = new int[ndt];
  for (int idt = 0; idt < ndt; ++idt) {
    dt_flag[idt] = 0;
    FOR_I4 {
      const int icv = dualTetVec[idt].icv[i];
      assert((icv >= 0)&&(icv < ncv_d));
      if (icv < ncv) {
        if (cv_flag[icv] == -1) {
          dt_flag[idt] = -1;
          break;
        }
      }
      else {
        if (cv_flag_d[icv-ncv] == -1) {
          dt_flag[idt] = -1;
          break;
        }
      }
    }
    if (dt_flag[idt] == 0)
      dt_flag[idt] = my_count[1]++;
  }

  int my_ncv = my_count[0];
  int my_ndt = my_count[1];

  int8 my_disp[2]; // nodal, element displacements..
  MPI_Scan(my_count,my_disp,2,MPI_INT8,MPI_SUM,mpi_comm);
  FOR_I2 my_disp[i] -= my_count[i];

  int8 count[2]; // everyone gets a copy of the node,element counts
  MPI_Allreduce(my_count,count,2,MPI_INT8,MPI_SUM,mpi_comm);

  if (mpi_rank == 0)
    cout << " > global cells and tets: " << count[0] << " " << count[1] << endl;

  // geom and scalar files need these array, so just allocate here...
  char * cbuf = new char[256];
  float * fbuf = new float[my_ncv];

  // for static solver the mesh is static, so we only need to build it once.
  // This assumes updateEnsightCasFile is called first!
  vector<MPI_Offset> offsets;
  if (dw->ensight_casefile_nsteps == 1) {

    // lets build the geo file used for all transient data associated to this ensight case...
    COUT1("StaticSolver::writeFlaggedCvsDualEnsightGeoFile()");

    // open the file...

    string filename = dw->name + ".geo";
    MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...
    MPI_File_delete((char*)filename.c_str(),MPI_INFO_NULL);
    MPI_File fh;
    MPI_File_open(mpi_comm,(char*)filename.c_str(),MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    // ==============================
    // header...
    // ==============================

    MPI_Offset offset = 0;
    if (mpi_rank == 0) {
      // write out header

      sprintf(cbuf,"%-80s","C Binary");
      MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Geometry file"); // should probably include hashes here
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Produced by Cascade Technologies, Inc.");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","node id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","element id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
    }
    offset = 5*80*sizeof(char);

    // all ranks write there flagged data...

    offset += mpi_rank*(int_size+2*80*sizeof(char)); // part headers
    offset += mpi_rank*(80*sizeof(char)+int_size*1); // coordinates, nno
    offset += sizeof(float)*3*my_disp[0]; // x_no,y_no,z_no
    offset += mpi_rank*(80*sizeof(char)+int_size*1); // element type, ncv
    offset += int_size*4*my_disp[1];

    // ==============================
    // part header...
    // ==============================

    sprintf(cbuf,"%-80s","part");
    writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
    int part = 1+mpi_rank;
    writeChunkedData<int>(fh,offset+80*sizeof(char),&part,1,mpi_comm);
    sprintf(cbuf,"%-80d",mpi_rank);
    writeChunkedData<char>(fh,offset+int_size+80*sizeof(char),cbuf,80,mpi_comm);
    offset += int_size+2*80*sizeof(char);

    // ==============================
    // coordinates...
    // ==============================

    sprintf(cbuf,"%-80s","coordinates"); // convex polyhedron format...
    writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
    writeChunkedData<int>(fh,offset+80*sizeof(char),&my_ncv,1,mpi_comm);
    offset += 80*sizeof(char)+int_size*1;

    FOR_I3 {
      int ind = 0;
      FOR_ICV
        if (cv_flag[icv] >= 0)
          fbuf[ind++] = (float)x_cv[icv][i];
      for (int icv = ncv; icv < ncv_d; ++icv)
        if (cv_flag_d[icv-ncv] >= 0)
          fbuf[ind++] = (float)x_cv_d[icv-ncv][i];
      assert(ind == my_ncv);
      writeChunkedData<float>(fh,offset,fbuf,my_ncv,mpi_comm);
      offset += sizeof(float)*my_ncv;
    }

    // ==============================
    // connectivity...
    // ==============================

    sprintf(cbuf,"%-80s","tetra4"); // convex polyhedron format...
    // write out element block header...
    writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
    writeChunkedData<int>(fh,offset+80*sizeof(char),&my_ndt,1,mpi_comm);
    offset += 80*sizeof(char)+int_size*1;

    int* ibuf = new int[4*my_ndt];
    {
      int ind = 0;
      for (int idt = 0; idt < ndt; ++idt) {
        if (dt_flag[idt] >= 0) {
          //ibuf[ind++] = dt_flag[idt]+1;
          FOR_I4 {
            const int icv = dualTetVec[idt].icv[i];
            if (icv < ncv) {
              assert(cv_flag[icv] >= 0);
              ibuf[ind++] = cv_flag[icv]+1;
            }
            else {
              assert(cv_flag_d[icv-ncv] >= 0);
              ibuf[ind++] = cv_flag_d[icv-ncv]+1;
            }
          }
        }
      }
      assert(ind == 4*my_ndt);
    }

    writeChunkedData<int>(fh,offset,ibuf,4*my_ndt,mpi_comm);
    delete[] ibuf;
    offset += int_size*4*my_ndt;

    // give everyone the same offset in end to set size

    offset = 5*80*sizeof(char);
    offset += mpi_size*(int_size+2*80*sizeof(char)); // part headers
    offset += mpi_size*(80*sizeof(char)+int_size); // coordinates, ncv
    offset += sizeof(float)*3*count[0]; // x_cv,y_cv,z_cv
    offset += mpi_size*(80*sizeof(char)+int_size); // element type, ntet
    offset += int_size*4*count[1];

    MPI_File_set_size(fh,offset);
    MPI_File_close(&fh);

    // add to offset vec...
    offsets.push_back(offset);
  }

  // now lets write out the scalar and vector files for this step...
  COUT1("StaticSolver::writeFlaggedCvsDualEnsightVarFiles()");

  if (dw->nscalars_cv > 0) {
    double *data_cv_d = new double[ncv_d-ncv];
    for (int isc = 0; isc < dw->nscalars_cv; ++isc) {

      // open file...
      sprintf(cbuf,"%s.%08d.sca%02d",dw->name.c_str(),step,isc); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);
      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->scalar_descs_cv[isc].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      }
      offset = 80*sizeof(char);
      offset += mpi_rank*(int_size+2*80*sizeof(char)); // part header
      offset += sizeof(float)*my_disp[0]; // part scalar nodal record

      sprintf(cbuf,"%-80s","part");
      writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
      offset += 80*sizeof(char);
      int part = 1+mpi_rank;
      writeChunkedData<int>(fh,offset,&part,1,mpi_comm);
      offset += int_size;
      sprintf(cbuf,"%-80s","coordinates");
      writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
      offset += 80*sizeof(char);

      // write data...
      CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->scalar_names_cv[isc]);
      assert((var->getType() == DN_DATA)&&(var->getUnindexedTopology() == CV_DATA));
      double* data_cv = var->getDNptr();
      updateCvdDataSeparateGhosts(data_cv,data_cv_d);
      {
        int ind = 0;
        FOR_ICV
          if (cv_flag[icv] >= 0)
            fbuf[ind++] = (float)data_cv[icv];
        for (int icv = ncv; icv < ncv_d; ++icv)
          if (cv_flag_d[icv-ncv] >= 0)
            fbuf[ind++] = (float)data_cv_d[icv-ncv];
        assert(ind == my_ncv);
      }
      writeChunkedData<float>(fh,offset,fbuf,my_ncv,mpi_comm);
      offset += sizeof(float)*my_ncv;

      // give everyone the same offset to size size

      offset = 80*sizeof(char);
      offset += mpi_size*(int_size+2*80*sizeof(char)); // part header
      offset += sizeof(float)*count[0]; // part scalar nodal record

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);

      // add to offset vec...
      offsets.push_back(offset);
    }
    delete[] data_cv_d;
  }
  if (dw->nvectors_cv > 0) {
    double (*data_cv_d)[3] = new double[ncv_d-ncv][3];
    for (int ivec = 0; ivec < dw->nvectors_cv; ++ivec) {

      // open file...
      sprintf(cbuf,"%s.%08d.vec%02d",dw->name.c_str(),step,ivec); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);
      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->vector_descs_cv[ivec].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      }
      offset += 80*sizeof(char);

      offset += mpi_rank*(int_size+2*80*sizeof(char)); // part header
      offset += 3*sizeof(float)*my_disp[0]; // part vector nodal record

      sprintf(cbuf,"%-80s","part");
      writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
      offset += 80*sizeof(char);
      int part = 1+mpi_rank;
      writeChunkedData<int>(fh,offset,&part,1,mpi_comm);
      offset += int_size;
      sprintf(cbuf,"%-80s","coordinates");
      writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
      offset += 80*sizeof(char);

      // write data...
      CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->vector_names_cv[ivec]);
      assert((var->getType() == DN3_DATA)&&(var->getUnindexedTopology() == CV_DATA));
      double (*data_cv)[3] = var->getDN3ptr();
      updateCvdDataSeparateGhosts(data_cv,data_cv_d);
      FOR_I3 {
        int ind = 0;
        FOR_ICV
          if (cv_flag[icv] >= 0)
            fbuf[ind++] = (float)data_cv[icv][i];
        for (int icv = ncv; icv < ncv_d; ++icv)
          if (cv_flag_d[icv-ncv] >= 0)
            fbuf[ind++] = (float)data_cv_d[icv-ncv][i];
        assert(ind == my_ncv);
        writeChunkedData<float>(fh,offset,fbuf,my_ncv,mpi_comm);
        offset += sizeof(float)*my_ncv;
      }

      // give everyone the same offset to size size

      offset = 80*sizeof(char);
      offset += mpi_size*(int_size+2*80*sizeof(char)); // part header
      offset += 3*sizeof(float)*count[0]; // part scalar nodal record

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);

      // add to offset vec...
      offsets.push_back(offset);
    }
    delete[] data_cv_d;
  }

  // cleanup...
  delete[] fbuf; fbuf = NULL;
  delete[] dt_flag;
  delete[] cv_flag;
  delete[] cv_flag_d;

  // boundary part of the write...
  {
    int ioffset = 0;

    vector<int> znofz;
    for (map<const string,int>::const_iterator it = bfZoneNameMap.begin(); it != bfZoneNameMap.end(); ++it) {
      const int iz = it->second; // zone index
      if (dw->bfzone_flag[iz]) 
        znofz.push_back(iz);
    }
    const int nfz = znofz.size();

    int* ssp_flag_fz = new int[subSurface->nsp*nfz];
    for (int ifz = 0; ifz < nfz; ++ifz) {
      int* const ssp_flag = &ssp_flag_fz[ifz*subSurface->nsp];
      for (int issp = 0; issp < subSurface->nsp; ++issp)
        ssp_flag[issp] = mpi_size;
      const int iz = znofz[ifz];
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
          const int issp = sspobf_v[pob];
          assert((issp >= 0)&&(issp < subSurface->nsp));
          ssp_flag[issp] = mpi_rank;
        }
      }
      subSurface->updateSpData(ssp_flag,MIN_DATA); // smallest rank wins
    }

    // points b/w zones are multiply written if the proc own it

    int *my_sp_count = new int[nfz];
    for (int ifz = 0; ifz < nfz; ++ifz) {
      const int iz = znofz[ifz];
      my_sp_count[ifz] = 0;
      int* const ssp_flag = &ssp_flag_fz[ifz*subSurface->nsp];
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
          const int issp = sspobf_v[pob];
          assert((issp >= 0)&&(issp < subSurface->nsp));
          if (ssp_flag[issp] == mpi_rank) {
            ++my_sp_count[ifz]; // surface point dumped by us
            ssp_flag[issp] = -1;
          }
        }
      }
      // reset...
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
          const int issp = sspobf_v[pob];
          assert((issp >= 0)&&(issp < subSurface->nsp));
          if (ssp_flag[issp] == -1)
            ssp_flag[issp] = mpi_rank;
          assert((ssp_flag[issp] >= 0)&&(ssp_flag[issp] <= mpi_size));
        }
      }
    }

    // figure out our nodal and element displacement...

    int *my_sp_disp = new int[nfz];
    MPI_Scan(my_sp_count,my_sp_disp,nfz,MPI_INT,MPI_SUM,mpi_comm);
    for (int ifz = 0; ifz < nfz; ++ifz)
      my_sp_disp[ifz] -= my_sp_count[ifz];

    // everyone gets the counts...

    int *sp_count = new int[nfz];
    MPI_Allreduce(my_sp_count,sp_count,nfz,MPI_INT,MPI_SUM,mpi_comm);

    int max_my_sp_count = 0;
    for (int ifz = 0; ifz < nfz; ++ifz)
      max_my_sp_count = max(max_my_sp_count,my_sp_count[ifz]);

    // geom and scalar files need these array, so just allocate here...
    assert(fbuf == NULL); fbuf = new float[max_my_sp_count];

    // for static solver the mesh is static, so we only need to build it once.
    // This assumes updateEnsightCasFile is called first!
    if (dw->ensight_casefile_nsteps == 1) {

      int* sst_flag_fz = new int[subSurface->nst*nfz];
      for (int ifz = 0; ifz < nfz; ++ifz) {
        int* const sst_flag = &sst_flag_fz[ifz*subSurface->nst];
        for (int isst = 0; isst < subSurface->nst; ++isst)
          sst_flag[isst] = mpi_size;
        const int iz = znofz[ifz];
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);
          for (int tob = sstobf_i[ibf]; tob != sstobf_i[ibf+1]; ++tob) {
            const int isst = sstobf_v[tob];
            assert((isst >= 0)&&(isst < subSurface->nst));
            sst_flag[isst] = mpi_rank;
          }
        }
        subSurface->updateStData(sst_flag,MIN_DATA); // smallest rank wins
      }

      // points b/w zones are multiply written if the proc own it

      int *my_st_count = new int[nfz]; 
      for (int ifz = 0; ifz < nfz; ++ifz) {
        const int iz = znofz[ifz];
        int* const sst_flag = &sst_flag_fz[ifz*subSurface->nst];
        my_st_count[ifz] = 0;
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);
          for (int tob = sstobf_i[ibf]; tob != sstobf_i[ibf+1]; ++tob) {
            const int isst = sstobf_v[tob];
            assert((isst >= 0)&&(isst < subSurface->nst));
            assert(dw->bfzone_flag[subSurface->znost[isst]]);
            if (sst_flag[isst] == mpi_rank) {
              ++my_st_count[ifz]; // surface tri dumped by us
              sst_flag[isst] = -1;
            }
          }
        }
        // reset...
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);
          for (int tob = sstobf_i[ibf]; tob != sstobf_i[ibf+1]; ++tob) {
            const int isst = sstobf_v[tob];
            assert((isst >= 0)&&(isst < subSurface->nst));
            if (sst_flag[isst] == -1)
              sst_flag[isst] = mpi_rank;
            assert((sst_flag[isst] >= 0)&&(sst_flag[isst] <= mpi_size));
          }
        }
      }

      // figure out our nodal and element displacement...

      int *my_st_disp = new int[nfz];
      MPI_Scan(my_st_count,my_st_disp,nfz,MPI_INT,MPI_SUM,mpi_comm);
      for (int ifz = 0; ifz < nfz; ++ifz)
        my_st_disp[ifz] -= my_st_count[ifz];

      // everyone gets the counts...

      int *st_count = new int[nfz];
      MPI_Allreduce(my_st_count,st_count,nfz,MPI_INT,MPI_SUM,mpi_comm);

      int max_my_st_count = 0;
      for (int ifz = 0; ifz < nfz; ++ifz)
        max_my_st_count = max(max_my_st_count,my_st_count[ifz]);

      // lets build the geo file used for all transient data associated to this ensight case...
      if (mpi_rank == 0)
        cout << "StaticSolver::writeFlaggedBfZonesDualEnsightGeoFile()" << endl;

      // open the file...

      string filename = dw->name + ".geo";
      // we are appending boundary data to the geo file, do not delete here...
      MPI_File fh;
      MPI_File_open(mpi_comm,(char*)filename.c_str(),MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
      // note file pointer set at end of geo file...
      MPI_Offset offset = offsets[ioffset];

      // ==============================
      // header...
      // ==============================

      // go through all flagged bf zones...

      int *ibuf = new int[3*max_my_st_count]; // global point indices for each tri
      int* ssp_flag2 = new int[subSurface->nsp];
      for (int ifz = 0; ifz < nfz; ++ifz) {
        const int iz = znofz[ifz];
        int* const sst_flag = &sst_flag_fz[ifz*subSurface->nst];
        int* const ssp_flag = &ssp_flag_fz[ifz*subSurface->nsp];

        // ==============================
        // part header...
        // ==============================

        if (mpi_rank == 0) {
          sprintf(cbuf,"%-80s","part");
          MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          int part = 1+mpi_size+ifz;
          MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
          sprintf(cbuf,"%-80s",bfZoneVec[iz].getName().c_str());
          MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        }
        offset += int_size+2*80*sizeof(char);

        // ==============================
        // part coordinates...
        // ==============================

        if (mpi_rank == 0) {
          sprintf(cbuf,"%-80s","coordinates"); // convex polyhedron format...
          MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          MPI_File_write_at(fh,offset+80*sizeof(char),&sp_count[ifz],1,MPI_INT,MPI_STATUS_IGNORE);
        }
        offset += 80*sizeof(char)+int_size*1;

        // write out x[3]...
        FOR_I3 {
          int ind = 0;
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
              const int issp = sspobf_v[pob];
              assert((issp >= 0)&&(issp < subSurface->nsp));
              if (ssp_flag[issp] == mpi_rank) {
                ssp_flag[issp] = -1;
                fbuf[ind++] = (float)subSurface->xp[issp][i]; // x[N],y[N],z[N]
              }
            }
          }
          // reset...
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
              const int issp = sspobf_v[pob];
              assert((issp >= 0)&&(issp < subSurface->nsp));
              if (ssp_flag[issp] == -1)
                ssp_flag[issp] = mpi_rank;
              assert((ssp_flag[issp] >= 0)&&(ssp_flag[issp] <= mpi_size));
            }
          }
          assert(ind == my_sp_count[ifz]);

          writeChunkedData<float>(fh,offset+(i*sp_count[ifz]+my_sp_disp[ifz])*sizeof(float),fbuf,my_sp_count[ifz],mpi_comm);
        }
        offset += sizeof(float)*3*sp_count[ifz];

        // ==============================
        // connectivity...
        // ==============================

        sprintf(cbuf,"%-80s","tria3"); // triangles
        if (mpi_rank == 0) {
          // write out element block header...
          MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          MPI_File_write_at(fh,offset+80*sizeof(char),&st_count[ifz],1,MPI_INT,MPI_STATUS_IGNORE);
        }
        offset += 80*sizeof(char)+int_size*1;

        // global index...
        int my_count0 = 0;
        for (int issp = 0; issp < subSurface->nsp; ++issp)
          ssp_flag2[issp] = -2;
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);
          for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
            const int issp = sspobf_v[pob];
            assert((issp >= 0)&&(issp < subSurface->nsp));
            // first time this node was visited by this zone
            if (ssp_flag[issp] >= 0) {
              assert(ssp_flag[issp] <= mpi_size);
              if (ssp_flag[issp] == mpi_rank) {
                ssp_flag2[issp] = my_sp_disp[ifz]+my_count0++; // this node is dumped by us
              }
              else if (ssp_flag[issp] < mpi_size) {
                ssp_flag2[issp] = -1; // this node gets dumped, but not by us
              }
              else {
                assert(ssp_flag[issp] == mpi_size);
                ssp_flag2[issp] = -2; // nobody owns this node
              }
              ssp_flag[issp] = -ssp_flag[issp]-1;
            }
          }
        }
        // reset...
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);
          for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
            const int issp = sspobf_v[pob];
            assert((issp >= 0)&&(issp < subSurface->nsp));
            if (ssp_flag[issp] < 0)
              ssp_flag[issp] = -ssp_flag[issp]-1;
            assert((ssp_flag[issp] >= 0)&&(ssp_flag[issp] <= mpi_size));
          }
        }
        assert(my_count0 == my_sp_count[ifz]);
        subSurface->updateSpData(ssp_flag2,MAX_DATA);

        int ind = 0;
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);
          for (int tob = sstobf_i[ibf]; tob != sstobf_i[ibf+1]; ++tob) {
            const int isst = sstobf_v[tob];
            assert((isst >= 0)&&(isst < subSurface->nst));
            if (sst_flag[isst] == mpi_rank) {
              sst_flag[isst] = -1;
              FOR_I3 {
                const int issp = subSurface->spost[isst][i];
                assert((issp >= 0)&&(issp < subSurface->nsp));
                ibuf[ind++] = ssp_flag2[issp]+1; // 1-indexed
              }
            }
          }
        }
        assert(ind == 3*my_st_count[ifz]);
        // reset...
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);
          for (int tob = sstobf_i[ibf]; tob != sstobf_i[ibf+1]; ++tob) {
            const int isst = sstobf_v[tob];
            assert((isst >= 0)&&(isst < subSurface->nst));
            if (sst_flag[isst] == -1)
              sst_flag[isst] = mpi_rank;
            assert((sst_flag[isst] >= 0)&&(sst_flag[isst] <= mpi_size));
          }
        }

        writeChunkedData<int>(fh,offset+my_st_disp[ifz]*3*int_size,ibuf,3*my_st_count[ifz],mpi_comm);
        offset += int_size*3*st_count[ifz];
      }
      delete[] ibuf;
      delete[] ssp_flag2;
      delete[] sst_flag_fz;
      delete[] st_count;
      delete[] my_st_count;
      delete[] my_st_disp;

      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);

      offsets[ioffset] = offset;
      ioffset++;
    }

    // now lets write out the scalar and vector files for this step...
    if (mpi_rank == 0)
      cout << "StaticSolver::writeFlaggedBfZonesDualEnsightVarFiles()" << endl;

    if ((dw->nscalars_cv + dw->nvectors_cv) > 0) {

      double *inv_ssp_wgt_sum = new double[subSurface->nsp];
      for (int issp = 0; issp < subSurface->nsp; ++issp)
        inv_ssp_wgt_sum[issp] = 0.0;
      for (int ifz = 0; ifz < nfz; ++ifz) {
        const int iz = znofz[ifz];
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
            const int issp = sspobf_v[pob];
            inv_ssp_wgt_sum[issp] += sspobf_wgt[pob];
          }
        }
      }
      subSurface->updateSpData(inv_ssp_wgt_sum,ADD_DATA);
      for (int issp = 0; issp < subSurface->nsp; ++issp) {
        if (inv_ssp_wgt_sum[issp] > 0.0)
          inv_ssp_wgt_sum[issp] = 1.0/inv_ssp_wgt_sum[issp];
      }
      dumpRange(inv_ssp_wgt_sum,subSurface->nsp,"inv_ssp_wgt_sum");

      if (dw->nscalars_cv > 0) {
        double *ssp_data = new double[subSurface->nsp];
        for (int isc = 0; isc < dw->nscalars_cv; ++isc) {

          // open file...
          sprintf(cbuf,"%s.%08d.sca%02d",dw->name.c_str(),step,isc); // see case file
          // We are appending data to the file, do not delete here...

          MPI_File fh;
          MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
          // note file pointer set at end...
          MPI_Offset offset = offsets[ioffset];

          // get data on subsurface points...
          for (int issp = 0; issp < subSurface->nsp; ++issp)
            ssp_data[issp] = 0.0;
          CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->scalar_names_cv[isc]);
          assert((var->getType() == DN_DATA)&&(var->getUnindexedTopology() == CV_DATA));
          double* var_d = var->getDNptr();
          for (int ifz = 0; ifz < nfz; ++ifz) {
            const int iz = znofz[ifz];
            for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
              const int icv = cvobf[ibf];
              for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
                const int issp = sspobf_v[pob];
                ssp_data[issp] += sspobf_wgt[pob]*var_d[icv];
              }
            }
          }
          subSurface->updateSpData(ssp_data,ADD_DATA);
          for (int issp = 0; issp < subSurface->nsp; ++issp)
            ssp_data[issp] *= inv_ssp_wgt_sum[issp];
          for (int ifz = 0; ifz < nfz; ++ifz) {
            const int iz = znofz[ifz];
            int* const ssp_flag = &ssp_flag_fz[ifz*subSurface->nsp];

            if (mpi_rank == 0) {
              sprintf(cbuf,"%-80s","part");
              MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
              int part = 1+mpi_size+ifz;
              MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
              sprintf(cbuf,"%-80s","coordinates");
              MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
            }
            offset += int_size+2*80*sizeof(char);

            // write data...
            int ind = 0;
            for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
              assert(zone_bf[ibf] == iz);
              for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
                const int issp = sspobf_v[pob];
                assert((issp >= 0)&&(issp < subSurface->nsp));
                if (ssp_flag[issp] == mpi_rank) {
                  ssp_flag[issp] = -1;
                  fbuf[ind++] = (float)ssp_data[issp];
                }
              }
            }
            // reset...
            for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
              assert(zone_bf[ibf] == iz);
              for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
                const int issp = sspobf_v[pob];
                assert((issp >= 0)&&(issp < subSurface->nsp));
                if (ssp_flag[issp] == -1)
                  ssp_flag[issp] = mpi_rank;
                assert((ssp_flag[issp] >= 0)&&(ssp_flag[issp] <= mpi_size));
              }
            }
            assert(ind == my_sp_count[ifz]);

            writeChunkedData<float>(fh,offset+my_sp_disp[ifz]*sizeof(float),fbuf,my_sp_count[ifz],mpi_comm);
            offset += sizeof(float)*sp_count[ifz];
          }

          // close file...
          MPI_File_set_size(fh,offset);
          MPI_File_close(&fh);

          offsets[ioffset] = offset;
          ioffset++;
        }
        delete[] ssp_data;
      }
      if (dw->nvectors_cv > 0) {
        double (*ssp_data)[3] = new double[subSurface->nsp][3];
        for (int ivec = 0; ivec < dw->nvectors_cv; ++ivec) {

          // open file...
          sprintf(cbuf,"%s.%08d.vec%02d",dw->name.c_str(),step,ivec); // see case file
          // we are appending data to the file, do not delete here...
          MPI_File fh;
          MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
          // note file pointer set at end...
          MPI_Offset offset = offsets[ioffset];

          // get data on subsurface points...
          for (int issp = 0; issp < subSurface->nsp; ++issp)
            FOR_I3 ssp_data[issp][i] = 0.0;
          CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->vector_names_cv[ivec]);
          assert((var->getType() == DN3_DATA)&&(var->getUnindexedTopology() == CV_DATA));
          double (*var_d3)[3] = var->getDN3ptr();
          for (int ifz = 0; ifz < nfz; ++ifz) {
            const int iz = znofz[ifz];
            for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
              const int icv = cvobf[ibf];
              for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
                const int issp = sspobf_v[pob];
                FOR_I3 ssp_data[issp][i] += sspobf_wgt[pob]*var_d3[icv][i];
              }
            }
          }
          subSurface->updateSpData(ssp_data,ADD_DATA);
          for (int issp = 0; issp < subSurface->nsp; ++issp)
            FOR_I3 ssp_data[issp][i] *= inv_ssp_wgt_sum[issp];

          for (int ifz = 0; ifz < nfz; ++ifz) {
            const int iz = znofz[ifz];
            int* const ssp_flag = &ssp_flag_fz[ifz*subSurface->nsp];

            if (mpi_rank == 0) {
              sprintf(cbuf,"%-80s","part");
              MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
              int part = 1+mpi_size+ifz;
              MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
              sprintf(cbuf,"%-80s","coordinates");
              MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
            }
            offset += int_size+2*80*sizeof(char);

            // write data...
            FOR_I3 {
              int ind = 0;
              for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
                assert(zone_bf[ibf] == iz);
                for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
                  const int issp = sspobf_v[pob];
                  assert((issp >= 0)&&(issp < subSurface->nsp));
                  if (ssp_flag[issp] == mpi_rank) {
                    ssp_flag[issp] = -1;
                    fbuf[ind++] = (float)ssp_data[issp][i];
                  }
                }
              }
              // reset...
              for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
                assert(zone_bf[ibf] == iz);
                for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
                  const int issp = sspobf_v[pob];
                  assert((issp >= 0)&&(issp < subSurface->nsp));
                  if (ssp_flag[issp] == -1)
                    ssp_flag[issp] = mpi_rank;
                  assert((ssp_flag[issp] >= 0)&&(ssp_flag[issp] <= mpi_size));
                }
              }
              assert(ind == my_sp_count[ifz]);

              writeChunkedData<float>(fh,offset+(i*sp_count[ifz]+my_sp_disp[ifz])*sizeof(float),fbuf,my_sp_count[ifz],mpi_comm);
            }
            offset += sizeof(float)*3*sp_count[ifz];
          }

          // close file...
          MPI_File_set_size(fh,offset);
          MPI_File_close(&fh);

          offsets[ioffset] = offset;
          ioffset++;
        }
        delete[] ssp_data;
      }

      delete[] inv_ssp_wgt_sum;
    }

    // cleanup...
    delete[] my_sp_count;
    delete[] my_sp_disp;
    delete[] sp_count;
    delete[] ssp_flag_fz;
    delete[] fbuf;
  }

  delete[] cbuf;

}

void StaticSolver::writeCvsAndFlaggedBfZonesEnsightMultiPart(DataWriter* dw,const int step,const double time) {

  bool (*hnofa_b)[2] = new bool[noofa_i[nfa]][2];
  bool * hnobf_b = new bool[noobf_i[nbf]];
  flagHangingNodes(hnofa_b,hnobf_b,NULL);

  // need to know which cv's have open faces so we can prism them up...
  if (checkParam("FLAG_OPEN_FACES_HACK")) {
    if (b_open_faces_cv == NULL) 
      b_open_faces_cv = new bool[ncv];
    FOR_ICV b_open_faces_cv[icv] = true; 
  }
  else {
    flagCvsWithOpenFaces(hnofa_b,hnobf_b,NULL);
  }

  // update case file with info about this step...
  updateEnsightCaseFile(dw,step,time);

  COUT1("Checking Face Orientation...");

  int* no_flag = new int[nno];
  FOR_INO no_flag[ino] = mpi_size;

  int* fa_sign = new int[nfa];
  FOR_IFA fa_sign[ifa] = 0;
 
  int my_face_sum = 0;

  FOR_ICV {
    for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
      const int ifa = faocv_v[foc];
      int ihn = IHN_VAL;
      int nof_l = noofa_i[ifa+1]-1;
      while (hnofa_b[nof_l][ihn]&&nof_l>noofa_i[ifa]){
        --nof_l;
      }
      double fa_n[3] = {0.0,0.0,0.0};
      int ino1 = noofa_v[nof_l];
      for(int nof = noofa_i[ifa]; nof<=nof_l; ++nof) {
        if (!hnofa_b[nof][ihn]){
          const int ino0 = ino1;
          ino1 = noofa_v[nof];
          no_flag[ino1] = mpi_rank;
          double tri_n[3] = TRI_NORMAL_2(x_no[ino0],x_no[ino1],x_fa[ifa]);
          FOR_I3 fa_n[i] += tri_n[i];
        }
      }
      //assert(DOT_PRODUCT(fa_n,n_fa[ifa])>0);
      if(DOT_PRODUCT(fa_n,n_fa[ifa])<0){
        fa_sign[ifa] = 1;
      }
    }
  }

  int* bf_sign = new int[nbf];
  FOR_IBF {

    bf_sign[ibf] = 0;
    //const int icv0 = cvobf[ibf];

    int nob_l = noobf_i[ibf+1]-1;
    while (hnobf_b[nob_l]&&nob_l>noobf_i[ibf]){
      --nob_l;
    }
    double bf_n[3] = {0.0,0.0,0.0};
    int ino1 = noobf_v[nob_l];
    for (int nob = noobf_i[ibf]; nob<=nob_l; ++nob) {
      if (!hnobf_b[nob]){
        const int ino0 = ino1;
        ino1 = noobf_v[nob];
        no_flag[ino1] = mpi_rank;
        double tri_n[3] = TRI_NORMAL_2(x_no[ino0],x_no[ino1],x_bf[ibf]);
        FOR_I3 bf_n[i] += tri_n[i];
      }
    }
    //assert(DOT_PRODUCT(bf_n,n_bf[ibf])>0);
    if(DOT_PRODUCT(bf_n,n_bf[ibf])<0){
      bf_sign[ibf] = 1;
    }
    my_face_sum += bf_sign[ibf];
  }

  FOR_IFA my_face_sum += fa_sign[ifa];

  int face_sum;
  MPI_Reduce(&my_face_sum,&face_sum,1,MPI_INT,MPI_SUM,0,mpi_comm);
  if (mpi_rank==0){
    if (face_sum>0)
      cout << "  Warning: face normal orientation may be flipped in " << face_sum << " cell(s) due to node ordering" << endl;
    else
      cout << "  Node ordering consistent with face normal orientation" << endl;
  }
  delete[] fa_sign;
  delete[] bf_sign;

  int8 my_count[4] = {0,0,0,0};

  FOR_INO {
    if (no_flag[ino] == mpi_rank)
      no_flag[ino] = my_count[0]++;     // this node is getting dumped
    else {
      assert(no_flag[ino] == mpi_size);
      no_flag[ino] = -2;                // nobody owns this node
    }
  }

  int8* cv_flag = new int8[ncv];
  FOR_ICV {
    // for flagged cells...
    if (b_open_faces_cv[icv]) {
      cv_flag[icv] = my_count[0]++;
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        int ihn = IHN_VAL;
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        int nnof = 0;
        for (int nof = nof_f ; nof <= nof_l; ++nof){
          if (!hnofa_b[nof][ihn]) //skip hanging nodes
            ++nnof;
        }
        if (nnof>=3){
          my_count[1]++; // prism for each face
          my_count[2] += nnof+1; // face for each edge + 1
          my_count[3] += nnof*4; // add cell face nodes
        }
      }
      for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
        const int ibf = bfocv_v[boc];
        const int nob_f = noobf_i[ibf];
        const int nob_l = noobf_i[ibf+1]-1;
        int nnob = 0;
        for (int nob = nob_f ; nob <= nob_l; ++nob){
          if (!hnobf_b[nob]) //skip hanging nodes
            ++nnob;
        }
        if (nnob>=3){
          my_count[1]++;
          my_count[2] += nnob+1;
          my_count[3] += nnob*4; 
        }
      }
    }
    else {
      cv_flag[icv] = -1;
      ++my_count[1]; // this cell is getting dumped
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        int ihn = IHN_VAL;
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        int nnof = 0;
        for (int nof = nof_f ; nof <= nof_l; ++nof){
          if (!hnofa_b[nof][ihn]) //skip hanging nodes
            ++nnof;
        }
        if (nnof>=3){
          ++my_count[2];       // add cell face
          my_count[3] += nnof; // add cell face nodes
        }
      }
      for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
        const int ibf = bfocv_v[boc];
        const int nob_f = noobf_i[ibf];
        const int nob_l = noobf_i[ibf+1]-1;
        int nnob = 0;
        for (int nob = nob_f ; nob <= nob_l; ++nob){
          if (!hnobf_b[nob]) //skip hanging nodes
            ++nnob;
        }
        if (nnob>=3){
          ++my_count[2];       // add cell boundary face
          my_count[3] += nnob; // add cell boundary face nodes
        }
      }
    }
  }
  int myNumNodes = my_count[0];
  int myNumCells = my_count[1];
  int myNumCellFaces = my_count[2];
  int myNumCellFaceNodes = my_count[3];

  int8 my_disp[4]; // nodal, element displacements..
  MPI_Scan(my_count,my_disp,4,MPI_INT8,MPI_SUM,mpi_comm);
  FOR_I4 my_disp[i] -= my_count[i];

  int8 count[4]; // everyone gets a copy of the node,element counts
  MPI_Allreduce(my_count,count,4,MPI_INT8,MPI_SUM,mpi_comm);

  if (mpi_rank == 0)
    cout << " > global nodes, cells, cell-faces, cell-face-nodes: " << count[0] << " " << count[1] << " " << count[2] << " " << count[3] << endl;

  // geom and scalar files need these array, so just allocate here...
  char * cbuf = new char[256];
  float * fbuf = new float[myNumNodes];

  // for static solver the mesh is static, so we only need to build it once.
  // This assumes updateEnsightCasFile is called first!
  vector<MPI_Offset> offsets;
  if (dw->ensight_casefile_nsteps == 1) {

    // lets build the geo file used for all transient data associated to this ensight case...
    COUT1("StaticSolver::writeFlaggedCvsEnsightGeoFile()");

    // open the file...

    string filename = dw->name + ".geo";
    MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...
    MPI_File_delete((char*)filename.c_str(),MPI_INFO_NULL);
    MPI_File fh;
    MPI_File_open(mpi_comm,(char*)filename.c_str(),MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    // ==============================
    // header...
    // ==============================

    MPI_Offset offset = 0;
    if (mpi_rank == 0) {
      // write out header

      sprintf(cbuf,"%-80s","C Binary");
      MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Geometry file"); // should probably include hashes here
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Produced by Cascade Technologies, Inc.");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","node id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","element id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
    }
    offset = 5*80*sizeof(char);

    // all ranks write there flagged data...

    offset += mpi_rank*(int_size+2*80*sizeof(char)); // part headers
    offset += mpi_rank*(80*sizeof(char)+int_size*1); // coordinates, nno
    offset += sizeof(float)*3*my_disp[0]; // x_no,y_no,z_no
    offset += mpi_rank*(80*sizeof(char)+int_size*1); // element type, ncv
    offset += int_size*my_disp[1];
    offset += int_size*my_disp[2];
    offset += int_size*my_disp[3];

    // ==============================
    // part header...
    // ==============================

    sprintf(cbuf,"%-80s","part");
    writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
    int part = 1+mpi_rank;
    writeChunkedData<int>(fh,offset+80*sizeof(char),&part,1,mpi_comm);
    sprintf(cbuf,"%-80d",mpi_rank);
    writeChunkedData<char>(fh,offset+int_size+80*sizeof(char),cbuf,80,mpi_comm);
    offset += int_size+2*80*sizeof(char);

    // ==============================
    // coordinates...
    // ==============================

    sprintf(cbuf,"%-80s","coordinates"); // convex polyhedron format...
    writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
    writeChunkedData<int>(fh,offset+80*sizeof(char),&myNumNodes,1,mpi_comm);
    offset += 80*sizeof(char)+int_size*1;

    FOR_I3 {
      int ind = 0;
      FOR_INO {
        if (no_flag[ino] >= 0) {
          fbuf[ind++] = (float)x_no[ino][i]; // x[N],y[N],z[N]
        }
      }
      FOR_ICV {
        if (cv_flag[icv] >= 0) {
          fbuf[ind++] = (float)x_cv[icv][i]; // x[N],y[N],z[N]
        }
      }
      assert(ind == myNumNodes);
      writeChunkedData<float>(fh,offset,fbuf,myNumNodes,mpi_comm);
      offset += sizeof(float)*myNumNodes;
    }

    // ==============================
    // connectivity...
    // ==============================

    sprintf(cbuf,"%-80s","nfaced"); // convex polyhedron format...
    // write out element block header...
    writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
    writeChunkedData<int>(fh,offset+80*sizeof(char),&myNumCells,1,mpi_comm);
    offset += 80*sizeof(char)+int_size*1;

    int *int_buf = new int[myNumCells]; // num faces for cell
    int *int_buf_2 = new int[myNumCellFaces]; // num nodes for face
    int *int_buf_3 = new int[myNumCellFaceNodes]; // node index for face node
    {
      int ind = 0, ind_2 = 0, ind_3 = 0;
      FOR_ICV {
        if (!b_open_faces_cv[icv]) {
          int_buf[ind] = 0;
          // loop on internal faces...
          for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
            const int ifa = faocv_v[foc];
            int ihn = IHN_VAL;
            const int nof_f = noofa_i[ifa];
            const int nof_l = noofa_i[ifa+1]-1;
            int nnof = 0;
            vector<int> ino_vec;
            // loop over nodes of face....
            for (int nof = nof_f ; nof <= nof_l; nof++){
              if (!hnofa_b[nof][ihn]){ //skip hanging nodes
                ++nnof;
                ino_vec.push_back(noofa_v[nof]);
                //int_buf_3[ind_3++] = (int)no_flag[noofa_v[nof]]+1; // node index
              }
            }
            if (nnof>=3) {
              ++int_buf[ind];            // cell face
              int_buf_2[ind_2++] = nnof; // cell face nodes
              if (cvofa[ifa][0] == icv){ // flip node ordering for outward normal 
                for (int iiv = 0; iiv<ino_vec.size(); ++iiv)
                  int_buf_3[ind_3++] = (int)no_flag[ino_vec[iiv]]+1; //node index
              }
              else{
                for (int iiv = ino_vec.size()-1; iiv>=0; --iiv)
                  int_buf_3[ind_3++] = (int)no_flag[ino_vec[iiv]]+1; //node index
              }
            }
//            else{
//              assert(0);
//            }
          }
          // loop on boundary faces...
          for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
            const int ibf = bfocv_v[boc];
            const int nob_f = noobf_i[ibf];
            const int nob_l = noobf_i[ibf+1]-1;
            int nnob = 0;
            vector<int> ino_vec;
            // loop over nodes of face...
            for (int nob = nob_f; nob <= nob_l; nob++){
              if (!hnobf_b[nob]){ //skip hanging nodes
                ++nnob;
                ino_vec.push_back(noobf_v[nob]);
                //int_buf_3[ind_3++] = (int)no_flag[noobf_v[nob]]+1; // node index
              }
            }
            if (nnob>=3){
              ++int_buf[ind];            // cell boundary face
              int_buf_2[ind_2++] = nnob; // cell boundary face nodes
              for (int iiv = 0; iiv<ino_vec.size(); ++iiv)
                int_buf_3[ind_3++] = (int)no_flag[ino_vec[iiv]]+1; //node index
            }
//            else{
//              assert(0);
//            }
          }
          ++ind;
        }
        else {
          for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
            const int ifa = faocv_v[foc];
            int ihn = IHN_VAL;
            int nof_f = noofa_i[ifa];
            int nof_l = noofa_i[ifa+1]-1;
            int nnof = 0;
            for (int nof = nof_f ; nof < noofa_i[ifa+1]; ++nof){
              if (!hnofa_b[nof][ihn]){ //skip hanging nodes
                ++nnof;
                nof_l = nof;
              }
              else if (nof_f==nof){ //ensure first node is not a hanging node
                ++nof_f;
              }
            }
            if (nnof>=3){
              assert(nof_l>nof_f); //nof_l and nof_f are both set to non-hanging nodes
              int_buf[ind++] = nnof+1; 
              int ino1 = noofa_v[nof_l];
              for (int nof = nof_f ; nof <= nof_l; nof++) {
                int_buf_2[ind_2++] = 3; 
                int ino0 = ino1;
                while (hnofa_b[nof][ihn])
                  ++nof;
                assert(nof<=nof_l);
                ino1 = noofa_v[nof];
                
                if (cvofa[ifa][0] == icv){ // flip node ordering for outward normal 
                  int_buf_3[ind_3++] = (int)cv_flag[icv]+1;
                  int_buf_3[ind_3++] = (int)no_flag[ino1]+1; 
                  int_buf_3[ind_3++] = (int)no_flag[ino0]+1; 
                }
                else{
                  int_buf_3[ind_3++] = (int)cv_flag[icv]+1;
                  int_buf_3[ind_3++] = (int)no_flag[ino0]+1; 
                  int_buf_3[ind_3++] = (int)no_flag[ino1]+1; 
                }
              }
              int_buf_2[ind_2++] = nnof; 

              if (cvofa[ifa][0] == icv){ // flip node ordering for outward normal 
                for (int nof = nof_f ; nof <= nof_l; nof++){
                  if (!hnofa_b[nof][ihn])
                    int_buf_3[ind_3++] = (int)no_flag[noofa_v[nof]]+1; 
                }
              }
              else {
                for (int nof = nof_l ; nof >= nof_f; nof--){
                  if (!hnofa_b[nof][ihn])
                    int_buf_3[ind_3++] = (int)no_flag[noofa_v[nof]]+1; 
                }
              }
            }
          }
          for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
            const int ibf = bfocv_v[boc];
            int nob_f = noobf_i[ibf];
            int nob_l = noobf_i[ibf+1]-1;
            int nnob = 0;
            for (int nob = nob_f ; nob < noobf_i[ibf+1]; ++nob){
              if (!hnobf_b[nob]){ //skip hanging nodes
                ++nnob;
                nob_l = nob;
              }
              else if (nob_f==nob){ //ensure first node is not a hanging node
                ++nob_f;
              }
            }
            if (nnob>=3){
              assert(nob_l>nob_f); //nob_l and nob_f are both set to non-hanging nodes
              int_buf[ind++] = nnob+1; 
              int ino1 = noobf_v[nob_l];
              for (int nob = nob_f ; nob <= nob_l; nob++) {
                int_buf_2[ind_2++] = 3; 
                int ino0 = ino1;
                while (hnobf_b[nob])
                  ++nob;
                assert(nob<=nob_l);
                ino1 = noobf_v[nob];
                int_buf_3[ind_3++] = (int)cv_flag[icv]+1;
                int_buf_3[ind_3++] = (int)no_flag[ino0]+1; 
                int_buf_3[ind_3++] = (int)no_flag[ino1]+1; 
              }
              int_buf_2[ind_2++] = nnob; 
              for (int nob = nob_f; nob <= nob_l; nob++){
                if (!hnobf_b[nob]) //skip hanging nodes
                  int_buf_3[ind_3++] = (int)no_flag[noobf_v[nob]]+1; 
              }
            }
          }
        }
      }
      assert(ind == myNumCells);
      assert(ind_2 == myNumCellFaces);
      assert(ind_3 == myNumCellFaceNodes);
    }

    writeChunkedData<int>(fh,offset,int_buf,myNumCells,mpi_comm);
    delete[] int_buf;
    offset += int_size*my_count[1];
    writeChunkedData<int>(fh,offset,int_buf_2,myNumCellFaces,mpi_comm);
    delete[] int_buf_2;
    offset += int_size*my_count[2];
    writeChunkedData<int>(fh,offset,int_buf_3,myNumCellFaceNodes,mpi_comm);
    delete[] int_buf_3;
    offset += int_size*my_count[3];

    // give everyone the same offset in end to set size

    offset = 5*80*sizeof(char);
    offset += mpi_size*(int_size+2*80*sizeof(char)); // part headers
    offset += mpi_size*(80*sizeof(char)+int_size); // coordinates, nno
    offset += sizeof(float)*3*count[0]; // x_no,y_no,z_no
    offset += mpi_size*(80*sizeof(char)+int_size); // element type, ncv
    offset += int_size*count[1];
    offset += int_size*count[2];
    offset += int_size*count[3];

    MPI_File_set_size(fh,offset);
    MPI_File_close(&fh);

    // add to offset vec...
    offsets.push_back(offset);
  }

  // now lets write out the scalar and vector files for this step...
  if (mpi_rank == 0)
    cout << "StaticSolver::writeFlaggedCvsEnsightVarFiles()" << endl;

  if (dw->nscalars_cv > 0) {
    double *no_data = new double[nno];
    for (int isc = 0; isc < dw->nscalars_cv; ++isc) {

      // open file...
      sprintf(cbuf,"%s.%08d.sca%02d",dw->name.c_str(),step,isc); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);
      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->scalar_descs_cv[isc].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      }
      offset = 80*sizeof(char);
      offset += mpi_rank*(int_size+2*80*sizeof(char)); // part header
      offset += sizeof(float)*my_disp[0]; // part scalar nodal record

      sprintf(cbuf,"%-80s","part");
      writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
      offset += 80*sizeof(char);
      int part = 1+mpi_rank;
      writeChunkedData<int>(fh,offset,&part,1,mpi_comm);
      offset += int_size;
      sprintf(cbuf,"%-80s","coordinates");
      writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
      offset += 80*sizeof(char);

      // write data...
      setNoDN(no_data,dw->scalar_names_cv[isc]);
      CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->scalar_names_cv[isc]);
      assert((var->getType() == DN_DATA)&&(var->getUnindexedTopology() == CV_DATA));
      double* var_d = var->getDNptr();
      {
        int ind = 0;
        FOR_INO {
          if (no_flag[ino] >= 0) {
            fbuf[ind++] = (float)no_data[ino];
          }
        }
        FOR_ICV {
          if (b_open_faces_cv[icv]) {
            fbuf[ind++] = (float)var_d[icv];
          }
        }
        assert(ind == myNumNodes);
      }
      writeChunkedData<float>(fh,offset,fbuf,myNumNodes,mpi_comm);
      offset += sizeof(float)*myNumNodes;

      // give everyone the same offset to size size

      offset = 80*sizeof(char);
      offset += mpi_size*(int_size+2*80*sizeof(char)); // part header
      offset += sizeof(float)*count[0]; // part scalar nodal record

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);

      // add to offset vec...
      offsets.push_back(offset);
    }
    delete[] no_data;
  }
  if (dw->nvectors_cv > 0) {
    double (*no_data)[3] = new double[nno][3];
    for (int ivec = 0; ivec < dw->nvectors_cv; ++ivec) {

      // open file...
      sprintf(cbuf,"%s.%08d.vec%02d",dw->name.c_str(),step,ivec); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);
      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->vector_descs_cv[ivec].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      }
      offset += 80*sizeof(char);

      offset += mpi_rank*(int_size+2*80*sizeof(char)); // part header
      offset += 3*sizeof(float)*my_disp[0]; // part vector nodal record

      sprintf(cbuf,"%-80s","part");
      writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
      offset += 80*sizeof(char);
      int part = 1+mpi_rank;
      writeChunkedData<int>(fh,offset,&part,1,mpi_comm);
      offset += int_size;
      sprintf(cbuf,"%-80s","coordinates");
      writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
      offset += 80*sizeof(char);

      // write data...
      setNoDN3(no_data,dw->vector_names_cv[ivec]);
      CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->vector_names_cv[ivec]);
      assert((var->getType() == DN3_DATA)&&(var->getUnindexedTopology() == CV_DATA));
      double (*var_d3)[3] = var->getDN3ptr();
      FOR_I3 {
        int ind = 0;
        FOR_INO {
          if (no_flag[ino] >= 0) {
            fbuf[ind++] = (float)no_data[ino][i];
          }
        }
        FOR_ICV {
          if (b_open_faces_cv[icv]) {
            fbuf[ind++] = (float)var_d3[icv][i];
          }
        }
        assert(ind == myNumNodes);
        writeChunkedData<float>(fh,offset,fbuf,myNumNodes,mpi_comm);
        offset += sizeof(float)*myNumNodes;
      }

      // give everyone the same offset to size size

      offset = 80*sizeof(char);
      offset += mpi_size*(int_size+2*80*sizeof(char)); // part header
      offset += 3*sizeof(float)*count[0]; // part scalar nodal record

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);

      // add to offset vec...
      offsets.push_back(offset);
    }
    delete[] no_data;
  }

  // cleanup...
  delete[] fbuf; fbuf = NULL; // reused below
  delete[] cv_flag;
  delete[] hnofa_b;
  delete[] hnobf_b;

  // boundary part of the write...
  {
    int ioffset = 0;
    FOR_INO no_flag[ino] = mpi_size;

    int nfz = 0;
    for (map<const string,int>::const_iterator it = bfZoneNameMap.begin(); it != bfZoneNameMap.end(); ++it) {
      const int iz = it->second; // zone index
      if (dw->bfzone_flag[iz]) {
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);
          for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
            const int ino = noobf_v[nob];
            assert((ino >= 0)&&(ino < nno));
            no_flag[ino] = mpi_rank;
          }
          int nnof = noobf_i[ibf+1]-noobf_i[ibf];
          assert( nnof >= 3 ); //check later ...
        }
        nfz += 1;
      }
    }

    int *znofz = new int[nfz]; // zn of flagged zn
    nfz = 0;
    for (map<const string,int>::const_iterator it = bfZoneNameMap.begin(); it != bfZoneNameMap.end(); ++it) {
      const int iz = it->second; // zone index
      if (dw->bfzone_flag[iz]) {
        znofz[nfz++] = iz;
      }
    }

    updateNoData(no_flag,MIN_NO_PERIODIC_DATA); // smallest rank wins

    // nodes b/w zones are multiply written if the proc own it

    bool *visited_no = new bool[nno]; FOR_INO visited_no[ino] = false;
    int (*my_count)[3] = new int[nfz][3]; // local ncount,ecount and necount
    for (int ifz = 0; ifz < nfz; ++ifz) {
      const int iz = znofz[ifz];
      // get zone node,element,node-element counts...
      my_count[ifz][0] = 0;
      my_count[ifz][1] = bfZoneVec[iz].nbf;
      my_count[ifz][2] = 0;
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        my_count[ifz][2] += noobf_i[ibf+1]-noobf_i[ibf];
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int ino = noobf_v[nob];
          assert((ino >= 0)&&(ino < nno));
          // first time this node was visited by this zone
          if (!visited_no[ino]) {
            visited_no[ino] = true;
            if (no_flag[ino] == mpi_rank) {
              my_count[ifz][0] += 1; // this node is dumped by us
            }
          }
        }
      }
      // reset visited_no...
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int ino = noobf_v[nob];
          assert((ino >= 0)&&(ino < nno));
          visited_no[ino] = false;
        }
      }
    }

    // figure out our nodal and element displacement...

    int (*my_disp)[3] = new int[nfz][3];
    MPI_Scan((int*)my_count,(int*)my_disp,3*nfz,MPI_INT,MPI_SUM,mpi_comm);
    for (int ifz = 0; ifz < nfz; ++ifz) {
      my_disp[ifz][0] -= my_count[ifz][0];
      my_disp[ifz][1] -= my_count[ifz][1];
      my_disp[ifz][2] -= my_count[ifz][2];
    }

    // everyone gets the counts...

    int (*count)[3] = new int[nfz][3];
    MPI_Allreduce((int*)my_count,(int*)count,3*nfz,MPI_INT,MPI_SUM,mpi_comm);

    int max_my_count[3] = {0,0,0};
    for (int ifz = 0; ifz < nfz; ++ifz) {
      FOR_I3 max_my_count[i] = max(max_my_count[i],my_count[ifz][i]);
    }

    // geom and scalar files need these array, so just allocate here...
    assert(fbuf == NULL); fbuf = new float[max_my_count[0]];

    // for static solver the mesh is static, so we only need to build it once.
    // This assumes updateEnsightCasFile is called first!
    if (dw->ensight_casefile_nsteps == 1) {

      // lets build the geo file used for all transient data associated to this ensight case...
      COUT1("StaticSolver::writeFlaggedBfZonesEnsightGeoFile()");

      // open the file...

      string filename = dw->name + ".geo";
      // we are appending boundary data to the geo file, do not delete here...
      MPI_File fh;
      MPI_File_open(mpi_comm,(char*)filename.c_str(),MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
      // note file pointer set at end of geo file...
      MPI_Offset offset = offsets[ioffset];

      // go through all flagged bf zones...

      int *ibuf = new int[max_my_count[1]];
      int *ibuf_2 = new int[max_my_count[2]];
      int* no_flag2 = new int[nno];
      for (int ifz = 0; ifz < nfz; ++ifz) {
        const int iz = znofz[ifz];

        // ==============================
        // part header...
        // ==============================

        if (mpi_rank == 0) {
          sprintf(cbuf,"%-80s","part");
          MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          int part = 1+mpi_size+ifz;
          MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
          sprintf(cbuf,"%-80s",bfZoneVec[iz].getName().c_str());
          MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        }
        offset += int_size+2*80*sizeof(char);

        // ==============================
        // part coordinates...
        // ==============================

        if (mpi_rank == 0) {
          sprintf(cbuf,"%-80s","coordinates"); // convex polyhedron format...
          MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          MPI_File_write_at(fh,offset+80*sizeof(char),&count[ifz][0],1,MPI_INT,MPI_STATUS_IGNORE);
        }
        offset += 80*sizeof(char)+int_size*1;

        // write out x[3]...
        FOR_I3 {
          int ind = 0;
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
              const int ino = noobf_v[nob];
              assert((ino >= 0)&&(ino < nno));
              // first time this node was visited by this zone
              if (!visited_no[ino]) {
                visited_no[ino] = true;
                if (no_flag[ino] == mpi_rank) {
                  fbuf[ind++] = (float)x_no[ino][i]; // x[N],y[N],z[N]
                }
              }
            }
          }
          // reset visited_no...
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
              const int ino = noobf_v[nob];
              assert((ino >= 0)&&(ino < nno));
              visited_no[ino] = false;
            }
          }
          assert(ind == my_count[ifz][0]);

          writeChunkedData<float>(fh,offset+(i*count[ifz][0]+my_disp[ifz][0])*sizeof(float),fbuf,my_count[ifz][0],mpi_comm);
        }
        offset += sizeof(float)*3*count[ifz][0];

        // ==============================
        // connectivity...
        // ==============================

        sprintf(cbuf,"%-80s","nsided"); // polygon...
        if (mpi_rank == 0) {
          // write out element block header...
          MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          MPI_File_write_at(fh,offset+80*sizeof(char),&count[ifz][1],1,MPI_INT,MPI_STATUS_IGNORE);
        }
        offset += 80*sizeof(char)+int_size*1;

        // put global index in no_flag2...
        int my_count0 = 0;
        FOR_INO no_flag2[ino] = -2;
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);
          for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
            const int ino = noobf_v[nob];
            assert((ino >= 0)&&(ino < nno));
            // first time this node was visited by this zone
            if (!visited_no[ino]) {
              visited_no[ino] = true;
              if (no_flag[ino] == mpi_rank) {
                no_flag2[ino] = my_disp[ifz][0]+my_count0++; // this node is dumped by us
              }
              else if (no_flag[ino] < mpi_size) {
                no_flag2[ino] = -1; // this node gets dumped, but not by us
              }
              else {
                assert(no_flag[ino] == mpi_size);
                no_flag2[ino] = -2; // nobody owns this node
              }
            }
          }
        }
        // reset visited_no...
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);
          for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
            const int ino = noobf_v[nob];
            assert((ino >= 0)&&(ino < nno));
            visited_no[ino] = false;
          }
        }
        assert(my_count0 == my_count[ifz][0]);
        updateNoData(no_flag2,MAX_NO_PERIODIC_DATA);

        int ind = 0;
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);

          // counts...
          ibuf[ibf-bfZoneVec[iz].ibf_f] = noobf_i[ibf+1]-noobf_i[ibf];

          // connectivity...
          for (int nof = noobf_i[ibf]; nof != noobf_i[ibf+1]; ++nof) {
            ibuf_2[ind++] = 1+no_flag2[noobf_v[nof]]; // 1-indexed
          }
        }
        assert(ind == my_count[ifz][2]);

        writeChunkedData<int>(fh,offset+my_disp[ifz][1]*int_size,ibuf,my_count[ifz][1],mpi_comm);
        offset += int_size*count[ifz][1];
        writeChunkedData<int>(fh,offset+my_disp[ifz][2]*int_size,ibuf_2,my_count[ifz][2],mpi_comm);
        offset += int_size*count[ifz][2];
      }
      delete[] ibuf;
      delete[] ibuf_2;
      delete[] no_flag2;

      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);

      offsets[ioffset] = offset;
      ioffset++;
    }

    // now lets write out the scalar and vector files for this step...
    COUT1("StaticSolver::writeFlaggedBfZonesEnsightVarFiles()");

    if (dw->nscalars_cv > 0) {
      double *no_data = new double[nno];
      for (int isc = 0; isc < dw->nscalars_cv; ++isc) {

        // open file...
        sprintf(cbuf,"%s.%08d.sca%02d",dw->name.c_str(),step,isc); // see case file
        // we are appending data to the file, do not delete here...
        MPI_File fh;
        MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
        // note file pointer set at end...
        MPI_Offset offset = offsets[ioffset];

        // get data on nodes...
        setNoDN(no_data,dw->scalar_names_cv[isc]);
        for (int ifz = 0; ifz < nfz; ++ifz) {
          const int iz = znofz[ifz];

          if (mpi_rank == 0) {
            sprintf(cbuf,"%-80s","part");
            MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
            int part = 1+mpi_size+ifz;
            MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
            sprintf(cbuf,"%-80s","coordinates");
            MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          }
          offset += int_size+2*80*sizeof(char);

          // write data...

          int ind = 0;
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
              const int ino = noobf_v[nob];
              assert((ino >= 0)&&(ino < nno));
              // first time this node was visited by this zone
              if (!visited_no[ino]) {
                visited_no[ino] = true;
                if (no_flag[ino] == mpi_rank) {
                  fbuf[ind++] = (float)no_data[ino]; // x[N],y[N],z[N]
                }
              }
            }
          }
          // reset visited_no...
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
              const int ino = noobf_v[nob];
              assert((ino >= 0)&&(ino < nno));
              visited_no[ino] = false;
            }
          }
          assert(ind == my_count[ifz][0]);

          writeChunkedData<float>(fh,offset+my_disp[ifz][0]*sizeof(float),fbuf,my_count[ifz][0],mpi_comm);

          offset += sizeof(float)*count[ifz][0];
        }

        // close file...
        MPI_File_set_size(fh,offset);
        MPI_File_close(&fh);

        offsets[ioffset] = offset;
        ioffset++;
      }
      delete[] no_data;
    }
    if (dw->nvectors_cv > 0) {
      double (*no_data)[3] = new double[nno][3];
      for (int ivec = 0; ivec < dw->nvectors_cv; ++ivec) {

        // open file...
        sprintf(cbuf,"%s.%08d.vec%02d",dw->name.c_str(),step,ivec); // see case file
        // we are appending data to the file, do not delete here...
        MPI_File fh;
        MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
        // note file pointer set at end...
        MPI_Offset offset = offsets[ioffset];

        setNoDN3(no_data,dw->vector_names_cv[ivec]);
        for (int ifz = 0; ifz < nfz; ++ifz) {
          const int iz = znofz[ifz];

          if (mpi_rank == 0) {
            sprintf(cbuf,"%-80s","part");
            MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
            int part = 1+mpi_size+ifz;
            MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
            sprintf(cbuf,"%-80s","coordinates");
            MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          }
          offset += int_size+2*80*sizeof(char);

          // write data...
          FOR_I3 {
            int ind = 0;
            for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
              assert(zone_bf[ibf] == iz);
              for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
                const int ino = noobf_v[nob];
                assert((ino >= 0)&&(ino < nno));
                // first time this node was visited by this zone
                if (!visited_no[ino]) {
                  visited_no[ino] = true;
                  if (no_flag[ino] == mpi_rank) {
                    fbuf[ind++] = (float)no_data[ino][i];
                  }
                }
              }
            }
            // reset visited_no...
            for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
              assert(zone_bf[ibf] == iz);
              for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
                const int ino = noobf_v[nob];
                assert((ino >= 0)&&(ino < nno));
                visited_no[ino] = false;
              }
            }
            assert(ind == my_count[ifz][0]);

            writeChunkedData<float>(fh,offset+(i*count[ifz][0]+my_disp[ifz][0])*sizeof(float),fbuf,my_count[ifz][0],mpi_comm);
          }
          offset += sizeof(float)*3*count[ifz][0];
        }

        // close file...
        MPI_File_set_size(fh,offset);
        MPI_File_close(&fh);

        offsets[ioffset] = offset;
        ioffset++;
      }
      delete[] no_data;
    }

    // cleanup
    delete[] my_count;
    delete[] my_disp;
    delete[] count;
    delete[] visited_no;
    delete[] znofz;
  }

  // cleanup...
  delete[] cbuf;
  delete[] fbuf;
  delete[] no_flag;
}

void StaticSolver::writeCvsAndFlaggedBfZonesEnsight(DataWriter* dw,const int step,const double time) {

  // need to know which cv's have open faces so we can prism them up...
  if (checkParam("FLAG_OPEN_FACES_HACK")) {
    if (b_open_faces_cv == NULL)
      b_open_faces_cv = new bool[ncv];
    FOR_ICV b_open_faces_cv[icv] = true;
  }
  else {
    flagCvsWithOpenFaces(NULL);
  }

  // update case file with info about this step...
  updateEnsightCaseFile(dw,step,time);

  // lets build the geo file used for all transient data associated to this ensight case...
  COUT1("StaticSolver::writeCvsAndFlaggedBfZonesEnsight()");

  // get counts and no flag...
  int8* no_flag = new int8[nno];
  FOR_INO no_flag[ino] = mpi_size;
  FOR_ICV {
    for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
      const int ifa = faocv_v[foc];
      for(int nof = noofa_i[ifa]; nof != noofa_i[ifa+1];++nof) {
        const int ino = noofa_v[nof];
        no_flag[ino]  = mpi_rank;
      }
    }
  }

  FOR_IBF {
    //const int icv0 = cvobf[ibf];
    for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
      const int ino = noobf_v[nob];
      no_flag[ino] = mpi_rank;
    }
  }

  int8 my_count[4] = {0,0,0,0}; // number of nodes, number of cells, number of cell faces, number of cell face nodes

  updateNoData(no_flag,MIN_NO_PERIODIC_DATA);
  FOR_INO {
    if (no_flag[ino] == mpi_rank)
      no_flag[ino] = my_count[0]++;     // this node is getting dumped
    else if (no_flag[ino] < mpi_size)
      no_flag[ino] = -1;                // this node gets dumped, but not by us
    else {
      assert(no_flag[ino] == mpi_size);
      no_flag[ino] = -2;                // nobody owns this node
    }
  }

  int8* cv_flag = new int8[ncv];
  FOR_ICV {
    if (b_open_faces_cv[icv]) {
      cv_flag[icv] = my_count[0]++;
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        my_count[1]++; // prism for each face
        const int ifa = faocv_v[foc];
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        const int nnof = nof_l - nof_f + 1;
        my_count[2] += nnof+1; // face for each edge + 1
        my_count[3] += nnof*4; // add cell face nodes
      }
      for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
        my_count[1]++;
        const int ibf = bfocv_v[boc];
        const int nob_f = noobf_i[ibf];
        const int nob_l = noobf_i[ibf+1]-1;
        const int nnob = nob_l - nob_f + 1;
        my_count[2] += nnob+1;
        my_count[3] += nnob*4;
      }
    }
    else {
      cv_flag[icv] = -1;
      my_count[1]++; // this cell is getting dumped
      my_count[2] += faocv_i[icv+1]-faocv_i[icv] + bfocv_i[icv+1]-bfocv_i[icv]; // cell faces
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        const int nnof = nof_l - nof_f + 1;
        my_count[3] += nnof; // add cell face nodes
      }
      for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
        const int ibf = bfocv_v[boc];
        const int nob_f = noobf_i[ibf];
        const int nob_l = noobf_i[ibf+1]-1;
        const int nnob = nob_l - nob_f + 1;
        my_count[3] += nnob; // add cell boundary face nodes
      }
    }
  }

  int8 my_disp[4]; // nodal, element displacements..
  MPI_Scan(my_count,my_disp,4,MPI_INT8,MPI_SUM,mpi_comm);
  FOR_I4 my_disp[i] -= my_count[i];

  // put our global numbers into the no_flag
  FOR_INO if (no_flag[ino] >= 0)
    no_flag[ino] += my_disp[0];
  FOR_ICV if (b_open_faces_cv[icv])
    cv_flag[icv] += my_disp[0];

  int8 count[4]; // everyone gets a copy of the node,element counts
  MPI_Allreduce(my_count,count,4,MPI_INT8,MPI_SUM,mpi_comm);

  if ((count[0] > TWO_BILLION)||(count[1] > TWO_BILLION)) {
    CERR(" > ensight io cannot support >2B node/element count currently");
  }

  int myNumNodes = my_count[0];
  int myNumCells = my_count[1];
  int myNumCellFaces = my_count[2];
  int myNumCellFaceNodes = my_count[3];
  int8 myNodeDisp = my_disp[0];
  int8 myCellDisp = my_disp[1];
  int8 myCellFaceDisp = my_disp[2];
  int8 myCellFaceNodeDisp = my_disp[3];
  int8 NumNodes = count[0];
  int8 NumCells = count[1];
  int8 NumCellFaces = count[2];
  int8 NumCellFaceNodes = count[3];

  if (mpi_rank == 0)
    cout << " > global counts: " << NumNodes << " " << NumCells << " " << NumCellFaces << " " << (double)NumCellFaceNodes/(double)NumCellFaces << endl;

  // geom and scalar files need these array, so just allocate here...
  char * cbuf = new char[256];
  float * fbuf = new float[myNumNodes];

  // for static solver the mesh is static, so we only need to build it once.
  //   This assumes updateEnsightCasFile is called first!
  vector<MPI_Offset> offsets;
  if (dw->ensight_casefile_nsteps == 1) {

    // open the file...

    string filename = dw->name + ".geo";
    MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...
    MPI_File_delete((char*)filename.c_str(),MPI_INFO_NULL);

    MPI_File fh;
    MPI_File_open(mpi_comm,(char*)filename.c_str(),MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    // ==============================
    // header...
    // ==============================

    MPI_Offset offset = 0;
    if (mpi_rank == 0) {
      // write out header

      sprintf(cbuf,"%-80s","C Binary");
      MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Geometry file"); // should probably include hashes here
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Produced by Cascade Technologies, Inc.");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","node id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","element id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","part");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      int part = 1;
      MPI_File_write_at(fh,offset,&part,1,MPI_INT,MPI_STATUS_IGNORE);
      offset += int_size;
      sprintf(cbuf,"%-80s","all");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
    }
    offset = int_size+7*80*sizeof(char);

    // ==============================
    // coordinates...
    // ==============================

    int nn = (int)NumNodes;
    sprintf(cbuf,"%-80s","coordinates"); // convex polyhedron format...
    if (mpi_rank == 0) {
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      MPI_File_write_at(fh,offset+80*sizeof(char),&nn,1,MPI_INT,MPI_STATUS_IGNORE);
    }
    offset += 80*sizeof(char)+int_size*1;

    FOR_I3 {
      int ind = 0;
      FOR_INO {
        if (no_flag[ino] >= 0) {
          fbuf[ind++] = (float)x_no[ino][i]; // x[N],y[N],z[N]
        }
      }
      FOR_ICV {
        if (b_open_faces_cv[icv]) {
          fbuf[ind++] = (float)x_cv[icv][i]; // x[N],y[N],z[N]
        }
      }
      assert(ind == myNumNodes);

      writeChunkedData<float>(fh,offset+(i*NumNodes+myNodeDisp)*sizeof(float),fbuf,myNumNodes,mpi_comm);
    }
    offset += sizeof(float)*3*NumNodes;

    // ==============================
    // connectivity...
    // ==============================

    // recall that the 0-indexed ino_global is in no_flag in only the
    // nodes that actually OWN the nodes. otherwise, there is a -1. Same for
    // faces, so we need to have temporary

    // also update the global node counts to
    // overwrite the -1 values on processors that did not
    // write data...

    int no_minus1_check = 0;
    int no_minus2_check = 0;
    FOR_INO {
      if (no_flag[ino] >= 0) {
        assert((no_flag[ino] >= myNodeDisp)&&(no_flag[ino] < myNodeDisp+myNumNodes));
      }
      else if (no_flag[ino] == -1) {
        ++no_minus1_check;
      }
      else {
        assert(no_flag[ino] == -2);
        ++no_minus2_check;
      }
    }

    updateNoData(no_flag,MAX_NO_PERIODIC_DATA);

    // should have gotten rid of all -1's...
    FOR_INO assert(no_flag[ino] != -1);

    sprintf(cbuf,"%-80s","nfaced"); // convex polyhedron format...
    int ne = (int)NumCells;
    if (mpi_rank == 0) {
      // write out element block header...
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      MPI_File_write_at(fh,offset+80*sizeof(char),&ne,1,MPI_INT,MPI_STATUS_IGNORE);
    }
    offset += 80*sizeof(char)+int_size*1;

    int *int_buf = new int[myNumCells]; // num faces for cell
    int *int_buf_2 = new int[myNumCellFaces]; // num nodes for face
    int *int_buf_3 = new int[myNumCellFaceNodes]; // node index for face node
    {
      int ind = 0, ind_2 = 0, ind_3 = 0;
      FOR_ICV {
        if (!b_open_faces_cv[icv]) {
          int_buf[ind++] = faocv_i[icv+1]-faocv_i[icv] + bfocv_i[icv+1]-bfocv_i[icv]; // cell faces
          // loop on internal faces...
          for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
            const int ifa = faocv_v[foc];
            const int nof_f = noofa_i[ifa];
            const int nof_l = noofa_i[ifa+1]-1;
            const int nnof = nof_l - nof_f + 1;
            int_buf_2[ind_2++] = nnof; // cell face nodes
            // loop over nodes of face....
            for (int nof = nof_f ; nof <= nof_l; nof++)
              int_buf_3[ind_3++] = (int)no_flag[noofa_v[nof]]+1; // node index
          }
          // loop on boundary faces...
          for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
            const int ibf = bfocv_v[boc];
            const int nob_f = noobf_i[ibf];
            const int nob_l = noobf_i[ibf+1]-1;
            const int nnob = nob_l - nob_f + 1;
            int_buf_2[ind_2++] = nnob; // cell boundary face nodes
            // loop over nodes of face...
            for (int nob = nob_f; nob <= nob_l; nob++)
              int_buf_3[ind_3++] = (int)no_flag[noobf_v[nob]]+1; // node index
          }
        }
        else {
          for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
            const int ifa = faocv_v[foc];
            const int nof_f = noofa_i[ifa];
            const int nof_l = noofa_i[ifa+1]-1;
            const int nnof = nof_l - nof_f + 1;
            int_buf[ind++] = nnof+1;
            int ino1 = noofa_v[nof_l];
            for (int nof = nof_f ; nof <= nof_l; nof++) {
              int_buf_2[ind_2++] = 3;
              int ino0 = ino1;
              ino1 = noofa_v[nof];
              int_buf_3[ind_3++] = (int)cv_flag[icv]+1;
              int_buf_3[ind_3++] = (int)no_flag[ino0]+1;
              int_buf_3[ind_3++] = (int)no_flag[ino1]+1;
            }
            int_buf_2[ind_2++] = nnof;
            for (int nof = nof_f ; nof <= nof_l; nof++)
              int_buf_3[ind_3++] = (int)no_flag[noofa_v[nof]]+1;
          }
          for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
            const int ibf = bfocv_v[boc];
            const int nob_f = noobf_i[ibf];
            const int nob_l = noobf_i[ibf+1]-1;
            const int nnob = nob_l - nob_f + 1;
            int_buf[ind++] = nnob+1;
            int ino1 = noobf_v[nob_l];
            for (int nob = nob_f ; nob <= nob_l; nob++) {
              int_buf_2[ind_2++] = 3;
              int ino0 = ino1;
              ino1 = noobf_v[nob];
              int_buf_3[ind_3++] = (int)cv_flag[icv]+1;
              int_buf_3[ind_3++] = (int)no_flag[ino0]+1;
              int_buf_3[ind_3++] = (int)no_flag[ino1]+1;
            }
            int_buf_2[ind_2++] = nnob;
            for (int nob = nob_f; nob <= nob_l; nob++)
              int_buf_3[ind_3++] = (int)no_flag[noobf_v[nob]]+1;
          }
        }
      }
      assert(ind == myNumCells);
      assert(ind_2 == myNumCellFaces);
      assert(ind_3 == myNumCellFaceNodes);
    }

    writeChunkedData<int>(fh,offset+myCellDisp*int_size,int_buf,myNumCells,mpi_comm);
    delete[] int_buf;
    offset += int_size*NumCells;
    writeChunkedData<int>(fh,offset+myCellFaceDisp*int_size,int_buf_2,myNumCellFaces,mpi_comm);
    delete[] int_buf_2;
    offset += int_size*NumCellFaces;
    writeChunkedData<int>(fh,offset+myCellFaceNodeDisp*int_size,int_buf_3,myNumCellFaceNodes,mpi_comm);
    delete[] int_buf_3;
    offset += int_size*NumCellFaceNodes;

    // now return no_flag to only flagged nodes...

    int no_minus1_check2 = 0;
    int no_minus2_check2 = 0;
    FOR_INO {
      if (no_flag[ino] >= 0) {
        if ((no_flag[ino] < myNodeDisp)||(no_flag[ino] >= myNodeDisp+myNumNodes)) {
          no_flag[ino] = -1;
          ++no_minus1_check2;
        }
      }
      else {
        assert(no_flag[ino] == -2);
        ++no_minus2_check2;
      }
    }
    assert(no_minus1_check2 == no_minus1_check);
    assert(no_minus2_check2 == no_minus2_check);

    // get out of here...
    MPI_File_set_size(fh,offset);
    MPI_File_close(&fh);

    // add to offset vec...
    offsets.push_back(offset);
  }

  if (dw->nscalars_cv > 0) {
    double *no_data = new double[nno];
    for (int isc = 0; isc < dw->nscalars_cv; ++isc) {

      // open file...
      sprintf(cbuf,"%s.%08d.sca%02d",dw->name.c_str(),step,isc); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);

      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->scalar_descs_cv[isc].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
        sprintf(cbuf,"%-80s","part");
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
        int part = 1; // 1 part for now
        MPI_File_write_at(fh,offset,&part,1,MPI_INT,MPI_STATUS_IGNORE);
        offset += int_size;
        sprintf(cbuf,"%-80s","coordinates");
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
      }
      offset = 3*80*sizeof(char)+int_size;

      // write data...
      setNoDN(no_data,dw->scalar_names_cv[isc]);
      CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->scalar_names_cv[isc]);
      assert((var->getType() == DN_DATA)&&(var->getUnindexedTopology() == CV_DATA));
      double* var_d = var->getDNptr();
      {
        int ind = 0;
        FOR_INO {
          if (no_flag[ino] >= 0) {
            fbuf[ind++] = (float)no_data[ino];
          }
        }
        FOR_ICV {
          if (b_open_faces_cv[icv]) {
            fbuf[ind++] = (float)var_d[icv];
          }
        }
        assert(ind == myNumNodes);
      }

      writeChunkedData<float>(fh,offset+myNodeDisp*sizeof(float),fbuf,myNumNodes,mpi_comm);

      offset += sizeof(float)*NumNodes;

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);

      // add to offset vec...
      offsets.push_back(offset);
    }
    delete[] no_data;
  }
  if (dw->nvectors_cv > 0) {
    double (*no_data)[3] = new double[nno][3];
    for (int ivec = 0; ivec < dw->nvectors_cv; ++ivec) {

      // open file...
      sprintf(cbuf,"%s.%08d.vec%02d",dw->name.c_str(),step,ivec); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);

      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->vector_descs_cv[ivec].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
        sprintf(cbuf,"%-80s","part");
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
        int part = 1; // 1 part for now
        MPI_File_write_at(fh,offset,&part,1,MPI_INT,MPI_STATUS_IGNORE);
        offset += int_size;
        sprintf(cbuf,"%-80s","coordinates");
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
      }
      offset = 3*80*sizeof(char)+int_size;

      // write data...
      setNoDN3(no_data,dw->vector_names_cv[ivec]);
      CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->vector_names_cv[ivec]);
      assert((var->getType() == DN3_DATA)&&(var->getUnindexedTopology() == CV_DATA));
      double (*var_d3)[3] = var->getDN3ptr();
      FOR_I3 {
        int ind = 0;
        FOR_INO {
          if (no_flag[ino] >= 0) {
            fbuf[ind++] = (float)no_data[ino][i];
          }
        }
        FOR_ICV {
          if (b_open_faces_cv[icv]) {
            fbuf[ind++] = (float)var_d3[icv][i];
          }
        }
        assert(ind == myNumNodes);
        writeChunkedData<float>(fh,offset+(i*NumNodes+myNodeDisp)*sizeof(float),fbuf,myNumNodes,mpi_comm);
      }
      offset += sizeof(float)*3*NumNodes;

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);

      // add to offset vec...
      offsets.push_back(offset);
    }
    delete[] no_data;
  }

  // cleanup...
  delete[] fbuf; fbuf = NULL; // reused below

  // now go through all FLAGGED boundary zones...
  int* bf_flag = new int[nbf];
  FOR_IZONE(bfZoneVec) if (dw->bfzone_flag[izone]) {

    // index into offset vec...
    int ioffset = 0;

    FOR_INO no_flag[ino] = mpi_size;
    FOR_IBF bf_flag[ibf] = mpi_size;

    for (int ibf = bfZoneVec[izone].ibf_f; ibf <= bfZoneVec[izone].ibf_l; ++ibf) {
      assert(zone_bf[ibf] == izone);
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino = noobf_v[nob];
        assert((ino >= 0)&&(ino < nno));
        no_flag[ino] = mpi_rank;
      }
      int nnof = noobf_i[ibf+1]-noobf_i[ibf];
      assert( nnof >= 3 ); //check later ...
      bf_flag[ibf] = mpi_rank;
    }

    int my_count[3] = {0,0,0}; // local ncount,ecount and necount

    updateNoData(no_flag,MIN_NO_PERIODIC_DATA); // smallest rank wins

    FOR_INO {
      if (no_flag[ino] == mpi_rank)
        no_flag[ino] = my_count[0]++;     // this node is getting dumped
      else if (no_flag[ino] < mpi_size)
        no_flag[ino] = -1;                // this node gets dumped, but not by us
      else {
        assert(no_flag[ino] == mpi_size);
        no_flag[ino] = -2;                // nobody owns this node
      }
    }

    // need to tiebreak faces, but we need to go through the cells
    // since we don't have a face communicator at this stage.
    //updateBfData(bf_flag,MIN_NO_PERIODIC_DATA);

    FOR_IBF {
      if (bf_flag[ibf] == mpi_rank) {
        bf_flag[ibf] = -1; // this but needs to be dumped ...
        my_count[1] += 1;
        const int nnof = noobf_i[ibf+1]-noobf_i[ibf];
        my_count[2] += nnof;
      }
      else {
        assert(bf_flag[ibf] == mpi_size); // it will fail if smaller rank wins ...
        bf_flag[ibf] = -2; // this face is not dumped
      }
    }

    // figure out our nodal and element displacement...

    int my_disp[3];
    MPI_Scan(my_count,my_disp,3,MPI_INT,MPI_SUM,mpi_comm);
    my_disp[0] -= my_count[0];
    my_disp[1] -= my_count[1];
    my_disp[2] -= my_count[2];

    for (int ino = 0; ino < nno; ++ino)
      if (no_flag[ino] >= 0)
        no_flag[ino] += my_disp[0];

    // everyone gets the counts...

    int count[3]={0,0,0};
    MPI_Allreduce(my_count,count,3,MPI_INT,MPI_SUM,mpi_comm);

    int myNumNodes = my_count[0];
    int myNumCells = my_count[1];
    int myNumCellNodes = my_count[2];
    int myNodeDisp = my_disp[0];
    int myCellDisp = my_disp[1];
    int myCellNodeDisp = my_disp[2];
    int NumNodes = count[0];
    int NumCells = count[1];
    int NumCellNodes = count[2];

    if (mpi_rank == 0)
      cout << " > global counts for bf zone " << izone << ": " << NumNodes << " " << NumCells << " " << (double)NumCellNodes/(double)NumCells << endl;

    // geom and scalar files need these array, so just allocate here...
    assert(fbuf == NULL); float * fbuf = new float[myNumNodes];

    // for static solver the mesh is static, so we only need to build it once.
    // This assumes updateEnsightCasFile is called first!
    if (dw->ensight_casefile_nsteps == 1) {

      // open the file...

      string filename = dw->name + ".geo";
      //We are appending boundary data to the geo file, do not delete here...

      MPI_File fh;
      MPI_File_open(mpi_comm,(char*)filename.c_str(),MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
      // note file pointer set at end...
      MPI_Offset offset = offsets[ioffset];

      // ==============================
      // header...
      // ==============================

      // write out header

      if (mpi_rank == 0) {

        sprintf(cbuf,"%-80s","part");
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
        int part = 2+izone;
        MPI_File_write_at(fh,offset,&part,1,MPI_INT,MPI_STATUS_IGNORE);
        offset += int_size;
        sprintf(cbuf,"%-80s",bfZoneVec[izone].getName().c_str());
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
      }
      offset = offsets[ioffset]+2*80*sizeof(char)+int_size;

      // ==============================
      // coordinates...
      // ==============================

      sprintf(cbuf,"%-80s","coordinates"); // convex polyhedron format...
      if (mpi_rank == 0) {
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        MPI_File_write_at(fh,offset+80*sizeof(char),&NumNodes,1,MPI_INT,MPI_STATUS_IGNORE);
      }
      offset += 80*sizeof(char)+int_size*1;

      FOR_I3 {
        int ind = 0;
        FOR_INO {
          if (no_flag[ino] >= 0) {
            fbuf[ind++] = (float)x_no[ino][i]; // x[N],y[N],z[N]
          }
        }
        assert(ind == myNumNodes);

        writeChunkedData<float>(fh,offset+(i*NumNodes+myNodeDisp)*sizeof(float),fbuf,myNumNodes,mpi_comm);
      }
      offset += sizeof(float)*3*NumNodes;

      // ==============================
      // connectivity...
      // ==============================

      sprintf(cbuf,"%-80s","nsided"); // polygon...
      if (mpi_rank == 0) {
        // write out element block header...
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        MPI_File_write_at(fh,offset+80*sizeof(char),&NumCells,1,MPI_INT,MPI_STATUS_IGNORE);
      }
      offset += 80*sizeof(char)+int_size*1;

      // also update the global node and face counts to
      // overwrite the -1 values on processors that did not
      // write data...

      int no_minus1_check = 0;
      int no_minus2_check = 0;
      FOR_INO {
        if (no_flag[ino] >= 0) {
          assert((no_flag[ino] >= myNodeDisp)&&(no_flag[ino] < myNodeDisp+myNumNodes));
        }
        else if (no_flag[ino] == -1) {
          ++no_minus1_check;
        }
        else {
          assert(no_flag[ino] == -2);
          ++no_minus2_check;
        }
      }

      updateNoData(no_flag,MAX_NO_PERIODIC_DATA);

      int *ibuf = new int[myNumCells];
      int *ibuf_2 = new int[myNumCellNodes];
      int ind = 0, ind_2 = 0;
      FOR_IBF {
        if (bf_flag[ibf] == -1) {

          // counts...
          ibuf[ind++] = noobf_i[ibf+1]-noobf_i[ibf];

          // connectivity...
          for (int nof = noobf_i[ibf]; nof != noobf_i[ibf+1]; ++nof) {
            ibuf_2[ind_2++] = 1+no_flag[noobf_v[nof]]; // 1-indexed
          }
        }
      }
      assert(ind == myNumCells);
      assert(ind_2 == myNumCellNodes);
      writeChunkedData<int>(fh,offset+myCellDisp*int_size,ibuf,myNumCells,mpi_comm);
      delete[] ibuf;
      offset += int_size*NumCells;
      writeChunkedData<int>(fh,offset+myCellNodeDisp*int_size,ibuf_2,myNumCellNodes,mpi_comm);
      delete[] ibuf_2;
      offset += int_size*NumCellNodes;

      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);

      offsets[ioffset] = offset;
      ioffset++;

      // now return no_flag to only flagged nodes...

      int no_minus1_check2 = 0;
      int no_minus2_check2 = 0;
      FOR_INO {
        if (no_flag[ino] >= 0) {
          if ((no_flag[ino] < myNodeDisp)||(no_flag[ino] >= myNodeDisp+myNumNodes)) {
            no_flag[ino] = -1;
            ++no_minus1_check2;
          }
        }
        else {
          assert(no_flag[ino] == -2);
          ++no_minus2_check2;
        }
      }
      assert(no_minus1_check2 == no_minus1_check);
      assert(no_minus2_check2 == no_minus2_check);
    }

    if (dw->nscalars_cv > 0) {
      double *no_data = new double[nno];
      for (int isc = 0; isc < dw->nscalars_cv; ++isc) {

        // open file...
        sprintf(cbuf,"%s.%08d.sca%02d",dw->name.c_str(),step,isc); // see case file
        //We are appending data to the file, do not delete here...

        MPI_File fh;
        MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
        // note file pointer set at end...
        MPI_Offset offset = offsets[ioffset];

        // write header...
        if (mpi_rank == 0) {
          sprintf(cbuf,"%-80s","part");
          MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          offset += 80*sizeof(char);
          int part = 2+izone;
          MPI_File_write_at(fh,offset,&part,1,MPI_INT,MPI_STATUS_IGNORE);
          offset += int_size;
          sprintf(cbuf,"%-80s","coordinates");
          MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          offset += 80*sizeof(char);
        }
        offset = offsets[ioffset]+2*80*sizeof(char)+int_size;

        // write data...
        setNoDN(no_data,dw->scalar_names_cv[isc]);
        {
          int ind = 0;
          FOR_INO {
            if (no_flag[ino] >= 0) {
              fbuf[ind++] = (float)no_data[ino];
            }
          }
          assert(ind == myNumNodes);
        }
        writeChunkedData<float>(fh,offset+myNodeDisp*sizeof(float),fbuf,myNumNodes,mpi_comm);
        offset += sizeof(float)*NumNodes;

        // close file...
        MPI_File_set_size(fh,offset);
        MPI_File_close(&fh);

        offsets[ioffset] = offset;
        ioffset++;
      }
      delete[] no_data;
    }
    if (dw->nvectors_cv > 0) {
      double (*no_data)[3] = new double[nno][3];
      for (int ivec = 0; ivec < dw->nvectors_cv; ++ivec) {

        // open file...
        sprintf(cbuf,"%s.%08d.vec%02d",dw->name.c_str(),step,ivec); // see case file
        //We are appending data to the file, do not delete here...

        MPI_File fh;
        MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
        // note file pointer set at end...
        MPI_Offset offset = offsets[ioffset];

        // write header...
        if (mpi_rank == 0) {
          sprintf(cbuf,"%-80s","part");
          MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          offset += 80*sizeof(char);
          int part = 2+izone;
          MPI_File_write_at(fh,offset,&part,1,MPI_INT,MPI_STATUS_IGNORE);
          offset += int_size;
          sprintf(cbuf,"%-80s","coordinates");
          MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          offset += 80*sizeof(char);
        }
        offset = offsets[ioffset]+2*80*sizeof(char)+int_size;

        // write data...
        setNoDN3(no_data,dw->vector_names_cv[ivec]);
        FOR_I3 {
          int ind = 0;
          FOR_INO {
            if (no_flag[ino] >= 0) {
              fbuf[ind++] = (float)no_data[ino][i];
            }
          }
          assert(ind == myNumNodes);
          writeChunkedData(fh,offset+(i*NumNodes+myNodeDisp)*sizeof(float),fbuf,myNumNodes,mpi_comm);
        }
        offset += sizeof(float)*3*NumNodes;

        // close file...
        MPI_File_set_size(fh,offset);
        MPI_File_close(&fh);


        offsets[ioffset] = offset;
        ioffset++;
      }
      delete[] no_data;
    }
    delete[] fbuf; fbuf = NULL;
  }

  delete[] cv_flag;
  delete[] no_flag;
  delete[] bf_flag;
  delete[] cbuf;
}

void StaticSolver::writeFlaggedCvsDualEnsightMultiPart(DataWriter* dw,const int step,const double time) {

  // update case file with info about this step...
  updateEnsightCaseFile(dw,step,time);

  // here we assume that flagged cvs are identified by (cv_flag[icv] != 0)...

  int8* cv_flag = new int8[ncv];
  FOR_ICV cv_flag[icv] = dw->cv_flag[icv];
  FOR_IBF cv_flag[cvobf[ibf]] = 0; // internal only
  int8* cv_flag_d = new int8[ncv_d-ncv];
  updateCvdDataSeparateGhosts(cv_flag,cv_flag_d);

  int8 my_count[2] = {0,0}; // flagged cvs (including dual ghosts), and flagged tets
  FOR_ICV {
    if (cv_flag[icv] != 0)
      cv_flag[icv] = my_count[0]++;
    else
      cv_flag[icv] = -1;
  }
  for (int icv = ncv; icv < ncv_d; ++icv) {
    if (cv_flag_d[icv-ncv] != 0)
      cv_flag_d[icv-ncv] = my_count[0]++;
    else
      cv_flag_d[icv-ncv] = -1;
  }

  const int ndt = dualTetVec.size();
  int* dt_flag = new int[ndt];
  for (int idt = 0; idt < ndt; ++idt) {
    dt_flag[idt] = 0;
    FOR_I4 {
      const int icv = dualTetVec[idt].icv[i];
      assert((icv >= 0)&&(icv < ncv_d));
      if (icv < ncv) {
        if (cv_flag[icv] == -1) {
          dt_flag[idt] = -1;
          break;
        }
      }
      else {
        if (cv_flag_d[icv-ncv] == -1) {
          dt_flag[idt] = -1;
          break;
        }
      }
    }
    if (dt_flag[idt] == 0)
      dt_flag[idt] = my_count[1]++;
  }

  int my_ncv = my_count[0];
  int my_ndt = my_count[1];

  int8 my_disp[2]; // nodal, element displacements..
  MPI_Scan(my_count,my_disp,2,MPI_INT8,MPI_SUM,mpi_comm);
  FOR_I2 my_disp[i] -= my_count[i];

  int8 count[2]; // everyone gets a copy of the node,element counts
  MPI_Allreduce(my_count,count,2,MPI_INT8,MPI_SUM,mpi_comm);

  if (mpi_rank == 0)
    cout << " > global cells and tets: " << count[0] << " " << count[1] << endl;

  // geom and scalar files need these array, so just allocate here...
  char * cbuf = new char[256];
  float * fbuf = new float[my_ncv];

  // for static solver the mesh is static, so we only need to build it once.
  // This assumes updateEnsightCasFile is called first!
  if (dw->ensight_casefile_nsteps == 1) {

    // lets build the geo file used for all transient data associated to this ensight case...
    COUT1("StaticSolver::writeFlaggedCvsDualEnsightGeoFile()");

    // open the file...

    string filename = dw->name + ".geo";
    MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...
    MPI_File_delete((char*)filename.c_str(),MPI_INFO_NULL);
    MPI_File fh;
    MPI_File_open(mpi_comm,(char*)filename.c_str(),MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    // ==============================
    // header...
    // ==============================

    MPI_Offset offset = 0;
    if (mpi_rank == 0) {
      // write out header

      sprintf(cbuf,"%-80s","C Binary");
      MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Geometry file"); // should probably include hashes here
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Produced by Cascade Technologies, Inc.");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","node id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","element id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
    }
    offset = 5*80*sizeof(char);

    // all ranks write there flagged data...

    offset += mpi_rank*(int_size+2*80*sizeof(char)); // part headers
    offset += mpi_rank*(80*sizeof(char)+int_size*1); // coordinates, nno
    offset += sizeof(float)*3*my_disp[0]; // x_no,y_no,z_no
    offset += mpi_rank*(80*sizeof(char)+int_size*1); // element type, ncv
    offset += int_size*4*my_disp[1];

    // ==============================
    // part header...
    // ==============================

    sprintf(cbuf,"%-80s","part");
    writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
    int part = 1+mpi_rank;
    writeChunkedData<int>(fh,offset+80*sizeof(char),&part,1,mpi_comm);
    sprintf(cbuf,"%-80d",mpi_rank);
    writeChunkedData<char>(fh,offset+int_size+80*sizeof(char),cbuf,80,mpi_comm);
    offset += int_size+2*80*sizeof(char);

    // ==============================
    // coordinates...
    // ==============================

    sprintf(cbuf,"%-80s","coordinates"); // convex polyhedron format...
    writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
    writeChunkedData<int>(fh,offset+80*sizeof(char),&my_ncv,1,mpi_comm);
    offset += 80*sizeof(char)+int_size*1;

    FOR_I3 {
      int ind = 0;
      FOR_ICV
        if (cv_flag[icv] >= 0)
          fbuf[ind++] = (float)x_cv[icv][i];
      for (int icv = ncv; icv < ncv_d; ++icv)
        if (cv_flag_d[icv-ncv] >= 0)
          fbuf[ind++] = (float)x_cv_d[icv-ncv][i];
      assert(ind == my_ncv);
      writeChunkedData<float>(fh,offset,fbuf,my_ncv,mpi_comm);
      offset += sizeof(float)*my_ncv;
    }

    // ==============================
    // connectivity...
    // ==============================

    sprintf(cbuf,"%-80s","tetra4"); // convex polyhedron format...
    // write out element block header...
    writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
    writeChunkedData<int>(fh,offset+80*sizeof(char),&my_ndt,1,mpi_comm);
    offset += 80*sizeof(char)+int_size*1;

    int* ibuf = new int[4*my_ndt];
    {
      int ind = 0;
      for (int idt = 0; idt < ndt; ++idt) {
        if (dt_flag[idt] >= 0) {
          //ibuf[ind++] = dt_flag[idt]+1;
          FOR_I4 {
            const int icv = dualTetVec[idt].icv[i];
            if (icv < ncv) {
              assert(cv_flag[icv] >= 0);
              ibuf[ind++] = cv_flag[icv]+1;
            }
            else {
              assert(cv_flag_d[icv-ncv] >= 0);
              ibuf[ind++] = cv_flag_d[icv-ncv]+1;
            }
          }
        }
      }
      assert(ind == 4*my_ndt);
    }

    writeChunkedData<int>(fh,offset,ibuf,4*my_ndt,mpi_comm);
    delete[] ibuf;
    offset += int_size*4*my_ndt;

    // give everyone the same offset in end to set size

    offset = 5*80*sizeof(char);
    offset += mpi_size*(int_size+2*80*sizeof(char)); // part headers
    offset += mpi_size*(80*sizeof(char)+int_size); // coordinates, ncv
    offset += sizeof(float)*3*count[0]; // x_cv,y_cv,z_cv
    offset += mpi_size*(80*sizeof(char)+int_size); // element type, ntet
    offset += int_size*4*count[1];

    MPI_File_set_size(fh,offset);
    MPI_File_close(&fh);
  }

  // now lets write out the scalar and vector files for this step...
  COUT1("StaticSolver::writeFlaggedCvsDualEnsightVarFiles()");

  if (dw->nscalars_cv > 0) {
    double *data_cv_d = new double[ncv_d-ncv];
    for (int isc = 0; isc < dw->nscalars_cv; ++isc) {

      // open file...
      sprintf(cbuf,"%s.%08d.sca%02d",dw->name.c_str(),step,isc); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);
      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->scalar_descs_cv[isc].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      }
      offset = 80*sizeof(char);
      offset += mpi_rank*(int_size+2*80*sizeof(char)); // part header
      offset += sizeof(float)*my_disp[0]; // part scalar nodal record

      sprintf(cbuf,"%-80s","part");
      writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
      offset += 80*sizeof(char);
      int part = 1+mpi_rank;
      writeChunkedData<int>(fh,offset,&part,1,mpi_comm);
      offset += int_size;
      sprintf(cbuf,"%-80s","coordinates");
      writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
      offset += 80*sizeof(char);

      // write data...
      CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->scalar_names_cv[isc]);
      assert((var->getType() == DN_DATA)&&(var->getUnindexedTopology() == CV_DATA));
      double* data_cv = var->getDNptr();
      updateCvdDataSeparateGhosts(data_cv,data_cv_d);
      {
        int ind = 0;
        FOR_ICV
          if (cv_flag[icv] >= 0)
            fbuf[ind++] = (float)data_cv[icv];
        for (int icv = ncv; icv < ncv_d; ++icv)
          if (cv_flag_d[icv-ncv] >= 0)
            fbuf[ind++] = (float)data_cv_d[icv-ncv];
        assert(ind == my_ncv);
      }
      writeChunkedData<float>(fh,offset,fbuf,my_ncv,mpi_comm);
      offset += sizeof(float)*my_ncv;

      // give everyone the same offset to size size

      offset = 80*sizeof(char);
      offset += mpi_size*(int_size+2*80*sizeof(char)); // part header
      offset += sizeof(float)*count[0]; // part scalar nodal record

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);

    }
    delete[] data_cv_d;
  }
  if (dw->nvectors_cv > 0) {
    double (*data_cv_d)[3] = new double[ncv_d-ncv][3];
    for (int ivec = 0; ivec < dw->nvectors_cv; ++ivec) {

      // open file...
      sprintf(cbuf,"%s.%08d.vec%02d",dw->name.c_str(),step,ivec); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);
      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->vector_descs_cv[ivec].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      }
      offset += 80*sizeof(char);

      offset += mpi_rank*(int_size+2*80*sizeof(char)); // part header
      offset += 3*sizeof(float)*my_disp[0]; // part vector nodal record

      sprintf(cbuf,"%-80s","part");
      writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
      offset += 80*sizeof(char);
      int part = 1+mpi_rank;
      writeChunkedData<int>(fh,offset,&part,1,mpi_comm);
      offset += int_size;
      sprintf(cbuf,"%-80s","coordinates");
      writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
      offset += 80*sizeof(char);

      // write data...
      CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->vector_names_cv[ivec]);
      assert((var->getType() == DN3_DATA)&&(var->getUnindexedTopology() == CV_DATA));
      double (*data_cv)[3] = var->getDN3ptr();
      updateCvdDataSeparateGhosts(data_cv,data_cv_d);
      FOR_I3 {
        int ind = 0;
        FOR_ICV
          if (cv_flag[icv] >= 0)
            fbuf[ind++] = (float)data_cv[icv][i];
        for (int icv = ncv; icv < ncv_d; ++icv)
          if (cv_flag_d[icv-ncv] >= 0)
            fbuf[ind++] = (float)data_cv_d[icv-ncv][i];
        assert(ind == my_ncv);
        writeChunkedData<float>(fh,offset,fbuf,my_ncv,mpi_comm);
        offset += sizeof(float)*my_ncv;
      }

      // give everyone the same offset to size size

      offset = 80*sizeof(char);
      offset += mpi_size*(int_size+2*80*sizeof(char)); // part header
      offset += 3*sizeof(float)*count[0]; // part scalar nodal record

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);
    }
    delete[] data_cv_d;
  }

  // cleanup...
  delete[] cbuf;
  delete[] fbuf;
  delete[] dt_flag;
  delete[] cv_flag;
  delete[] cv_flag_d;
}

void StaticSolver::writeFlaggedCvsEnsightMultiPart(DataWriter* dw,int8* cv_flag,const int step,const double time) {

  bool (*hnofa_b)[2] = new bool[noofa_i[nfa]][2];
  bool * hnobf_b = new bool[noobf_i[nbf]];
  flagHangingNodes(hnofa_b,hnobf_b,cv_flag);

  // need to know which cv's have open faces so we can prism them up...
  if (checkParam("FLAG_OPEN_FACES_HACK")) {
    if (b_open_faces_cv == NULL) 
      b_open_faces_cv = new bool[ncv];
    FOR_ICV b_open_faces_cv[icv] = true; 
  }
  else {
    flagCvsWithOpenFaces(hnofa_b,hnobf_b,cv_flag);
  }

  // here we assume that flagged cvs are identified by (cv_flag[icv] != 0)...

  // update case file with info about this step...
  updateEnsightCaseFile(dw,step,time);


  COUT1("Checking Face Orientation...");
  int* fa_sign = new int[nfa];
  FOR_IFA fa_sign[ifa] = 0;
  int my_face_sum = 0;

  int* no_flag = new int[nno];
  FOR_INO no_flag[ino] = mpi_size;
  FOR_ICV {
    if (cv_flag[icv] != 0) { // this cv gets dumped
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        int ihn = IHN_VAL;
        int nof_l = noofa_i[ifa+1]-1;
        while (hnofa_b[nof_l][ihn]&&nof_l>noofa_i[ifa]){
          --nof_l;
        }
        double fa_n[3] = {0.0,0.0,0.0};
        int ino1 = noofa_v[nof_l];
        for(int nof = noofa_i[ifa]; nof<=nof_l; ++nof) {
          if (!hnofa_b[nof][ihn]){
            const int ino0 = ino1;
            ino1 = noofa_v[nof];
            no_flag[ino1] = mpi_rank;
            double tri_n[3] = TRI_NORMAL_2(x_no[ino0],x_no[ino1],x_fa[ifa]);
            FOR_I3 fa_n[i] += tri_n[i];
          }
        }
        //assert(DOT_PRODUCT(fa_n,n_fa[ifa])>0);
        if(DOT_PRODUCT(fa_n,n_fa[ifa])<0){
          fa_sign[ifa] = 1;
        }
      }
    }
  }

  int* bf_sign = new int[nbf];
  FOR_IBF {
    bf_sign[ibf] = 0;
    const int icv0 = cvobf[ibf];
    if (cv_flag[icv0] != 0) {
      int nob_l = noobf_i[ibf+1]-1;
      while (hnobf_b[nob_l]&&nob_l>noobf_i[ibf]){
        --nob_l;
      }
      double bf_n[3] = {0.0,0.0,0.0};
      int ino1 = noobf_v[nob_l];
      for (int nob = noobf_i[ibf]; nob<=nob_l; ++nob) {
        if (!hnobf_b[nob]){
          const int ino0 = ino1;
          ino1 = noobf_v[nob];
          no_flag[ino1] = mpi_rank;
          double tri_n[3] = TRI_NORMAL_2(x_no[ino0],x_no[ino1],x_bf[ibf]);
          FOR_I3 bf_n[i] += tri_n[i];
        }
      }
      //assert(DOT_PRODUCT(bf_n,n_bf[ibf])>0);
      if(DOT_PRODUCT(bf_n,n_bf[ibf])<0){
        bf_sign[ibf] = 1;
      }
      my_face_sum += bf_sign[ibf];
    }
  }

  FOR_IFA my_face_sum += fa_sign[ifa];

  int face_sum;
  MPI_Reduce(&my_face_sum,&face_sum,1,MPI_INT,MPI_SUM,0,mpi_comm);
  if (mpi_rank==0){
    if (face_sum>0)
      cout << "  Warning: face normal orientation may be flipped in " << face_sum << " cell(s) due to node ordering" << endl;
    else
      cout << "  Node ordering consistent with face normal orientation" << endl;
  }
  delete[] fa_sign;
  delete[] bf_sign;

  int8 my_count[4] = {0,0,0,0};

  FOR_INO {
    if (no_flag[ino] == mpi_rank)
      no_flag[ino] = my_count[0]++;     // this node is getting dumped
    else {
      assert(no_flag[ino] == mpi_size);
      no_flag[ino] = -2;                // nobody owns this node
    }
  }

  FOR_ICV {
    if (cv_flag[icv] == 0)
      cv_flag[icv] = -2;
    else 
      cv_flag[icv] = -1;
  }


  FOR_ICV {
    // for flagged cells...
    if (cv_flag[icv] == -1) {

      if (b_open_faces_cv[icv]) {
        cv_flag[icv] = my_count[0]++;
        for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
          const int ifa = faocv_v[foc];
          int ihn = IHN_VAL;
          const int nof_f = noofa_i[ifa];
          const int nof_l = noofa_i[ifa+1]-1;
          int nnof = 0;
          for (int nof = nof_f ; nof <= nof_l; ++nof){
            if (!hnofa_b[nof][ihn]) //skip hanging nodes
              ++nnof;
          }
          if (nnof>=3){
            my_count[1]++; // prism for each face
            my_count[2] += nnof+1; // face for each edge + 1
            my_count[3] += nnof*4; // add cell face nodes
          }
        }
        for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
          const int ibf = bfocv_v[boc];
          const int nob_f = noobf_i[ibf];
          const int nob_l = noobf_i[ibf+1]-1;
          int nnob = 0;
          for (int nob = nob_f ; nob <= nob_l; ++nob){
            if (!hnobf_b[nob]) //skip hanging nodes
              ++nnob;
          }
          if (nnob>=3){
            my_count[1]++;
            my_count[2] += nnob+1;
            my_count[3] += nnob*4; 
          }
        }
      }
      else {
        ++my_count[1]; // this cell is getting dumped
        for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
          const int ifa = faocv_v[foc];
          int ihn = IHN_VAL;
          const int nof_f = noofa_i[ifa];
          const int nof_l = noofa_i[ifa+1]-1;
          int nnof = 0;
          for (int nof = nof_f ; nof <= nof_l; ++nof){
            if (!hnofa_b[nof][ihn]) //skip hanging nodes
              ++nnof;
          }
          if (nnof>=3){
            ++my_count[2];       // add cell face
            my_count[3] += nnof; // add cell face nodes
          }
//          else{
//            assert(0);
//          }
        }
        for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
          const int ibf = bfocv_v[boc];
          const int nob_f = noobf_i[ibf];
          const int nob_l = noobf_i[ibf+1]-1;
          int nnob = 0;
          for (int nob = nob_f ; nob <= nob_l; ++nob){
            if (!hnobf_b[nob]) //skip hanging nodes
              ++nnob;
          }
          if (nnob>=3){
            ++my_count[2];       // add cell boundary face
            my_count[3] += nnob; // add cell boundary face nodes
          }
//          else{
//            assert(0);
//          }
        }
      }
    }
  }

  int myNumNodes = my_count[0];
  int myNumCells = my_count[1];
  int myNumCellFaces = my_count[2];
  int myNumCellFaceNodes = my_count[3];

  int8 my_disp[4]; // nodal, element displacements..
  MPI_Scan(my_count,my_disp,4,MPI_INT8,MPI_SUM,mpi_comm);
  FOR_I4 my_disp[i] -= my_count[i];

  int8 count[4]; // everyone gets a copy of the node,element counts
  MPI_Allreduce(my_count,count,4,MPI_INT8,MPI_SUM,mpi_comm);

  if (mpi_rank == 0)
    cout << " > global nodes, cells, cell-faces, cell-face-nodes: " << count[0] << " " << count[1] << " " << count[2] << " " << count[3] << endl;

  // geom and scalar files need these array, so just allocate here...
  char * cbuf = new char[256];
  float * fbuf = new float[myNumNodes];

  // for static solver the mesh is static, so we only need to build it once.
  // This assumes updateEnsightCasFile is called first!
  if (dw->ensight_casefile_nsteps == 1) {

    // lets build the geo file used for all transient data associated to this ensight case...
    COUT1("StaticSolver::writeFlaggedCvsEnsightGeoFile()");

    // open the file...

    string filename = dw->name + ".geo";
    MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...
    MPI_File_delete((char*)filename.c_str(),MPI_INFO_NULL);
    MPI_File fh;
    MPI_File_open(mpi_comm,(char*)filename.c_str(),MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    // ==============================
    // header...
    // ==============================

    MPI_Offset offset = 0;
    if (mpi_rank == 0) {
      // write out header

      sprintf(cbuf,"%-80s","C Binary");
      MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Geometry file"); // should probably include hashes here
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Produced by Cascade Technologies, Inc.");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","node id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","element id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
    }
    offset = 5*80*sizeof(char);

    // all ranks write there flagged data...

    offset += mpi_rank*(int_size+2*80*sizeof(char)); // part headers
    offset += mpi_rank*(80*sizeof(char)+int_size*1); // coordinates, nno
    offset += sizeof(float)*3*my_disp[0]; // x_no,y_no,z_no
    offset += mpi_rank*(80*sizeof(char)+int_size*1); // element type, ncv
    offset += int_size*my_disp[1];
    offset += int_size*my_disp[2];
    offset += int_size*my_disp[3];

    // ==============================
    // part header...
    // ==============================

    sprintf(cbuf,"%-80s","part");
    writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
    int part = 1+mpi_rank;
    writeChunkedData<int>(fh,offset+80*sizeof(char),&part,1,mpi_comm);
    sprintf(cbuf,"%-80d",mpi_rank);
    writeChunkedData<char>(fh,offset+int_size+80*sizeof(char),cbuf,80,mpi_comm);
    offset += int_size+2*80*sizeof(char);

    // ==============================
    // coordinates...
    // ==============================

    sprintf(cbuf,"%-80s","coordinates"); // convex polyhedron format...
    writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
    writeChunkedData<int>(fh,offset+80*sizeof(char),&myNumNodes,1,mpi_comm);
    offset += 80*sizeof(char)+int_size*1;

    FOR_I3 {
      int ind = 0;
      FOR_INO {
        if (no_flag[ino] >= 0) {
          fbuf[ind++] = (float)x_no[ino][i]; // x[N],y[N],z[N]
        }
      }
      FOR_ICV {
        if (cv_flag[icv] >= 0) {
          fbuf[ind++] = (float)x_cv[icv][i]; // x[N],y[N],z[N]
        }
      }
      assert(ind == myNumNodes);
      writeChunkedData<float>(fh,offset,fbuf,myNumNodes,mpi_comm);
      offset += sizeof(float)*myNumNodes;
    }

    // ==============================
    // connectivity...
    // ==============================

    sprintf(cbuf,"%-80s","nfaced"); // convex polyhedron format...
    // write out element block header...
    writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
    writeChunkedData<int>(fh,offset+80*sizeof(char),&myNumCells,1,mpi_comm);
    offset += 80*sizeof(char)+int_size*1;

    int *int_buf = new int[myNumCells]; // num faces for cell
    int *int_buf_2 = new int[myNumCellFaces]; // num nodes for face
    int *int_buf_3 = new int[myNumCellFaceNodes]; // node index for face node
    {
      int ind = 0, ind_2 = 0, ind_3 = 0;
      FOR_ICV {
        // for flagged cells...
        if (cv_flag[icv] > -2) {
          if (!b_open_faces_cv[icv]) {
//            int_buf[ind++] = faocv_i[icv+1]-faocv_i[icv] + bfocv_i[icv+1]-bfocv_i[icv]; // cell faces
            int_buf[ind] = 0;
            // loop on internal faces...
            for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
              const int ifa = faocv_v[foc];
              int ihn = IHN_VAL;
              const int nof_f = noofa_i[ifa];
              const int nof_l = noofa_i[ifa+1]-1;
              int nnof = 0;
              vector<int> ino_vec;
              // loop over nodes of face....
              for (int nof = nof_f ; nof <= nof_l; nof++){
                if (!hnofa_b[nof][ihn]){ //skip hanging nodes
                  ++nnof;
                  ino_vec.push_back(noofa_v[nof]);
                  //int_buf_3[ind_3++] = (int)no_flag[noofa_v[nof]]+1; // node index
                }
              }
              if (nnof>=3) {
                ++int_buf[ind];            // cell face
                int_buf_2[ind_2++] = nnof; // cell face nodes
                if (cvofa[ifa][0] == icv){ // flip node ordering for outward normal 
                  for (int iiv = 0; iiv<ino_vec.size(); ++iiv)
                    int_buf_3[ind_3++] = (int)no_flag[ino_vec[iiv]]+1; //node index
                }
                else{
                  for (int iiv = ino_vec.size()-1; iiv>=0; --iiv)
                    int_buf_3[ind_3++] = (int)no_flag[ino_vec[iiv]]+1; //node index
                }
              }
//              else{
//                assert(0);
//              }
            }
            // loop on boundary faces...
            for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
              const int ibf = bfocv_v[boc];
              const int nob_f = noobf_i[ibf];
              const int nob_l = noobf_i[ibf+1]-1;
              int nnob = 0;
              vector<int> ino_vec;
              // loop over nodes of face...
              for (int nob = nob_f; nob <= nob_l; nob++){
                if (!hnobf_b[nob]){ //skip hanging nodes
                  ++nnob;
                  ino_vec.push_back(noobf_v[nob]);
                  //int_buf_3[ind_3++] = (int)no_flag[noobf_v[nob]]+1; // node index
                }
              }
              if (nnob>=3){
                ++int_buf[ind];            // cell boundary face
                int_buf_2[ind_2++] = nnob; // cell boundary face nodes
                for (int iiv = 0; iiv<ino_vec.size(); ++iiv)
                  int_buf_3[ind_3++] = (int)no_flag[ino_vec[iiv]]+1; //node index
              }
//              else{
//                assert(0);
//              }
            }
            ++ind;
          }
          else {
            for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
              const int ifa = faocv_v[foc];
              int ihn = IHN_VAL;
              int nof_f = noofa_i[ifa];
              int nof_l = noofa_i[ifa+1]-1;
              int nnof = 0;
              for (int nof = nof_f ; nof < noofa_i[ifa+1]; ++nof){
                if (!hnofa_b[nof][ihn]){ //skip hanging nodes
                  ++nnof;
                  nof_l = nof;
                }
                else if (nof_f==nof){ //ensure first node is not a hanging node
                  ++nof_f;
                }
              }
              if (nnof>=3){
                assert(nof_l>nof_f); //nof_l and nof_f are both set to non-hanging nodes
                int_buf[ind++] = nnof+1; 
                int ino1 = noofa_v[nof_l];
                for (int nof = nof_f ; nof <= nof_l; nof++) {
                  int_buf_2[ind_2++] = 3; 
                  int ino0 = ino1;
                  while (hnofa_b[nof][ihn])
                    ++nof;
                  assert(nof<=nof_l);
                  ino1 = noofa_v[nof];
                  if (cvofa[ifa][0] == icv){ // flip node ordering for outward normal 
                    int_buf_3[ind_3++] = (int)cv_flag[icv]+1;
                    int_buf_3[ind_3++] = (int)no_flag[ino1]+1; 
                    int_buf_3[ind_3++] = (int)no_flag[ino0]+1; 
                  }
                  else{
                    int_buf_3[ind_3++] = (int)cv_flag[icv]+1;
                    int_buf_3[ind_3++] = (int)no_flag[ino0]+1; 
                    int_buf_3[ind_3++] = (int)no_flag[ino1]+1; 
                  }
                }
                int_buf_2[ind_2++] = nnof; 
                if (cvofa[ifa][0] == icv){ // flip node ordering for outward normal 
                  for (int nof = nof_f ; nof <= nof_l; nof++){
                    if (!hnofa_b[nof][ihn])
                      int_buf_3[ind_3++] = (int)no_flag[noofa_v[nof]]+1; 
                  }
                }
                else {
                  for (int nof = nof_l ; nof >= nof_f; nof--){
                    if (!hnofa_b[nof][ihn])
                      int_buf_3[ind_3++] = (int)no_flag[noofa_v[nof]]+1; 
                  }
                }
              }
            }
            for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
              const int ibf = bfocv_v[boc];
              int nob_f = noobf_i[ibf];
              int nob_l = noobf_i[ibf+1]-1;
              int nnob = 0;
              for (int nob = nob_f ; nob < noobf_i[ibf+1]; ++nob){
                if (!hnobf_b[nob]){ //skip hanging nodes
                  ++nnob;
                  nob_l = nob;
                }
                else if (nob_f==nob){ //ensure first node is not a hanging node
                  ++nob_f;
                }
              }
              if (nnob>=3){
                assert(nob_l>nob_f); //nob_l and nob_f are both set to non-hanging nodes
                int_buf[ind++] = nnob+1; 
                int ino1 = noobf_v[nob_l];
                for (int nob = nob_f ; nob <= nob_l; nob++) {
                  int_buf_2[ind_2++] = 3; 
                  int ino0 = ino1;
                  while (hnobf_b[nob])
                    ++nob;
                  assert(nob<=nob_l);
                  ino1 = noobf_v[nob];
                  int_buf_3[ind_3++] = (int)cv_flag[icv]+1;
                  int_buf_3[ind_3++] = (int)no_flag[ino0]+1; 
                  int_buf_3[ind_3++] = (int)no_flag[ino1]+1; 
                }
                int_buf_2[ind_2++] = nnob; 
                for (int nob = nob_f; nob <= nob_l; nob++){
                  if (!hnobf_b[nob]) //skip hanging nodes
                    int_buf_3[ind_3++] = (int)no_flag[noobf_v[nob]]+1; 
                }
              }
            }
          }
        }
      }
//      FOR_RANK{
//        if (mpi_rank == rank)
//          cout << "Rank " << rank << " ind " << ind << " myNumCells " << myNumCells << endl;
//        MPI_Barrier(mpi_comm);
 //     }
      assert(ind == myNumCells);
      assert(ind_2 == myNumCellFaces);
      assert(ind_3 == myNumCellFaceNodes);
    }

    writeChunkedData<int>(fh,offset,int_buf,myNumCells,mpi_comm);
    delete[] int_buf;
    offset += int_size*my_count[1];
    writeChunkedData<int>(fh,offset,int_buf_2,myNumCellFaces,mpi_comm);
    delete[] int_buf_2;
    offset += int_size*my_count[2];
    writeChunkedData<int>(fh,offset,int_buf_3,myNumCellFaceNodes,mpi_comm);
    delete[] int_buf_3;
    offset += int_size*my_count[3];

    // give everyone the same offset in end to set size

    offset = 5*80*sizeof(char);
    offset += mpi_size*(int_size+2*80*sizeof(char)); // part headers
    offset += mpi_size*(80*sizeof(char)+int_size); // coordinates, nno
    offset += sizeof(float)*3*count[0]; // x_no,y_no,z_no
    offset += mpi_size*(80*sizeof(char)+int_size); // element type, ncv
    offset += int_size*count[1];
    offset += int_size*count[2];
    offset += int_size*count[3];

    MPI_File_set_size(fh,offset);
    MPI_File_close(&fh);
  }

  // now lets write out the scalar and vector files for this step...
  COUT1("StaticSolver::writeFlaggedCvsEnsightVarFiles()");

  if (dw->nscalars_cv > 0) {
    double *no_data = new double[nno];
    for (int isc = 0; isc < dw->nscalars_cv; ++isc) {

      // open file...
      sprintf(cbuf,"%s.%08d.sca%02d",dw->name.c_str(),step,isc); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);
      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->scalar_descs_cv[isc].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      }
      offset = 80*sizeof(char);
      offset += mpi_rank*(int_size+2*80*sizeof(char)); // part header
      offset += sizeof(float)*my_disp[0]; // part scalar nodal record

      sprintf(cbuf,"%-80s","part");
      writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
      offset += 80*sizeof(char);
      int part = 1+mpi_rank;
      writeChunkedData<int>(fh,offset,&part,1,mpi_comm);
      offset += int_size;
      sprintf(cbuf,"%-80s","coordinates");
      writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
      offset += 80*sizeof(char);

      // write data...
      setNoDN(no_data,dw->scalar_names_cv[isc]);
      CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->scalar_names_cv[isc]);
      assert((var->getType() == DN_DATA)&&(var->getUnindexedTopology() == CV_DATA));
      double* var_d = var->getDNptr();
      {
        int ind = 0;
        FOR_INO {
          if (no_flag[ino] >= 0) {
            fbuf[ind++] = (float)no_data[ino];
          }
        }
        FOR_ICV {
          if (cv_flag[icv] >= 0) {
            fbuf[ind++] = (float)var_d[icv];
          }
        }
        assert(ind == myNumNodes);
      }
      writeChunkedData<float>(fh,offset,fbuf,myNumNodes,mpi_comm);
      offset += sizeof(float)*myNumNodes;

      // give everyone the same offset to size size

      offset = 80*sizeof(char);
      offset += mpi_size*(int_size+2*80*sizeof(char)); // part header
      offset += sizeof(float)*count[0]; // part scalar nodal record

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);

    }
    delete[] no_data;
  }
  if (dw->nvectors_cv > 0) {
    double (*no_data)[3] = new double[nno][3];
    for (int ivec = 0; ivec < dw->nvectors_cv; ++ivec) {

      // open file...
      sprintf(cbuf,"%s.%08d.vec%02d",dw->name.c_str(),step,ivec); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);
      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->vector_descs_cv[ivec].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      }
      offset += 80*sizeof(char);

      offset += mpi_rank*(int_size+2*80*sizeof(char)); // part header
      offset += 3*sizeof(float)*my_disp[0]; // part vector nodal record

      sprintf(cbuf,"%-80s","part");
      writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
      offset += 80*sizeof(char);
      int part = 1+mpi_rank;
      writeChunkedData<int>(fh,offset,&part,1,mpi_comm);
      offset += int_size;
      sprintf(cbuf,"%-80s","coordinates");
      writeChunkedData<char>(fh,offset,cbuf,80,mpi_comm);
      offset += 80*sizeof(char);

      // write data...
      setNoDN3(no_data,dw->vector_names_cv[ivec]);
      CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->vector_names_cv[ivec]);
      assert((var->getType() == DN3_DATA)&&(var->getUnindexedTopology() == CV_DATA));
      double (*var_d3)[3] = var->getDN3ptr();
      FOR_I3 {
        int ind = 0;
        FOR_INO {
          if (no_flag[ino] >= 0) {
            fbuf[ind++] = (float)no_data[ino][i];
          }
        }
        FOR_ICV {
          if (cv_flag[icv] >= 0) {
            fbuf[ind++] = (float)var_d3[icv][i];
          }
        }
        assert(ind == myNumNodes);
        writeChunkedData<float>(fh,offset,fbuf,myNumNodes,mpi_comm);
        offset += sizeof(float)*myNumNodes;
      }

      // give everyone the same offset to size size

      offset = 80*sizeof(char);
      offset += mpi_size*(int_size+2*80*sizeof(char)); // part header
      offset += 3*sizeof(float)*count[0]; // part scalar nodal record

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);
    }
    delete[] no_data;
  }

  // cleanup...
  delete[] cbuf;
  delete[] fbuf;
  delete[] no_flag;

  delete[] hnofa_b;
  delete[] hnobf_b;
}

void StaticSolver::writeFlaggedCvsEnsight(DataWriter* dw,int8* cv_flag,const int step,const double time) {

  // update case file with info about this step...
  updateEnsightCaseFile(dw,step,time);

  // get counts and no flag...
  int8* no_flag = new int8[nno];
  FOR_INO no_flag[ino] = mpi_size;
  FOR_ICV {
    if (cv_flag[icv] != 0) { // this cv gets dumped
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        for(int nof = noofa_i[ifa]; nof != noofa_i[ifa+1];++nof) {
          const int ino = noofa_v[nof];
          no_flag[ino]  = mpi_rank;
        }
      }
    }
  }

  FOR_IBF {
    const int icv0 = cvobf[ibf];
    if (cv_flag[icv0] != 0) {
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino = noobf_v[nob];
        no_flag[ino] = mpi_rank;
      }
    }
  }

  int8 my_count[4] = {0,0,0,0}; // number of nodes, number of cells, number of cell faces, number of cell face nodes

  updateNoData(no_flag,MIN_NO_PERIODIC_DATA);
  FOR_INO {
    if (no_flag[ino] == mpi_rank)
      no_flag[ino] = my_count[0]++;     // this node is getting dumped
    else if (no_flag[ino] < mpi_size)
      no_flag[ino] = -1;                // this node gets dumped, but not by us
    else {
      assert(no_flag[ino] == mpi_size);
      no_flag[ino] = -2;                // nobody owns this node
    }
  }

  FOR_ICV {
    // for flagged cells...
    if (cv_flag[icv] != 0) {
      my_count[1]++; // this cell is getting dumped
      my_count[2] += faocv_i[icv+1]-faocv_i[icv] + bfocv_i[icv+1]-bfocv_i[icv]; // add cell faces
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        const int nof_f = noofa_i[ifa];
        const int nof_l = noofa_i[ifa+1]-1;
        const int nnof = nof_l - nof_f + 1;
        my_count[3] += nnof; // add cell face nodes
      }
      for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
        const int ibf = bfocv_v[boc];
        const int nob_f = noobf_i[ibf];
        const int nob_l = noobf_i[ibf+1]-1;
        const int nnob = nob_l - nob_f + 1;
        my_count[3] += nnob; // add cell boundary face nodes
      }
    }
  }

  int8 my_disp[4]; // nodal, element displacements..
  MPI_Scan(my_count,my_disp,4,MPI_INT8,MPI_SUM,mpi_comm);
  FOR_I4 my_disp[i] -= my_count[i];

  // put our global numbers into the no_flag
  FOR_INO if (no_flag[ino] >= 0)
    no_flag[ino] += my_disp[0];

  int8 count[4]; // everyone gets a copy of the node,element counts
  MPI_Allreduce(my_count,count,4,MPI_INT8,MPI_SUM,mpi_comm);

  if ((count[0] > TWO_BILLION)||(count[1] > TWO_BILLION)) {
    CERR(" > ensight io cannot support >2B node/element count currently");
  }

  // here we assume that flagged cvs are identified by (cv_flag[icv] != 0)...
  int myNumNodes = my_count[0];
  int myNumCells = my_count[1];
  int myNumCellFaces = my_count[2];
  int myNumCellFaceNodes = my_count[3];
  int8 myNodeDisp = my_disp[0];
  int8 myCellDisp = my_disp[1];
  int8 myCellFaceDisp = my_disp[2];
  int8 myCellFaceNodeDisp = my_disp[3];
  int8 NumNodes = count[0];
  int8 NumCells = count[1];
  int8 NumCellFaces = count[2];
  int8 NumCellFaceNodes = count[3];

  if (mpi_rank == 0)
    cout << " > global counts: " << NumNodes << " " << NumCells << " " << NumCellFaces << " " << (double)NumCellFaceNodes/(double)NumCellFaces << endl;

  // geom and scalar files need these array, so just allocate here...
  char * cbuf = new char[256];
  float * fbuf = new float[myNumNodes];

  // for static solver the mesh is static, so we only need to build it once.
  // This assumes updateEnsightCasFile is called first!
  if (dw->ensight_casefile_nsteps == 1) {

    // lets build the geo file used for all transient data associated to this ensight case...
    COUT1("StaticSolver::writeFlaggedCvsEnsightGeoFile()");

    // open the file...

    string filename = dw->name + ".geo";
    MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...
    MPI_File_delete((char*)filename.c_str(),MPI_INFO_NULL);
    MPI_File fh;
    MPI_File_open(mpi_comm,(char*)filename.c_str(),MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    // ==============================
    // header...
    // ==============================

    MPI_Offset offset = 0;
    if (mpi_rank == 0) {
      // write out header

      sprintf(cbuf,"%-80s","C Binary");
      MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Geometry file"); // should probably include hashes here
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Produced by Cascade Technologies, Inc.");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","node id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","element id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","part");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      int part = 1; // 1 part for now
      MPI_File_write_at(fh,offset,&part,1,MPI_INT,MPI_STATUS_IGNORE);
      offset += int_size;
      sprintf(cbuf,"%-80s","Flagged cells");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
    }
    offset = int_size+7*80*sizeof(char);

    // ==============================
    // coordinates...
    // ==============================

    int nn = (int)NumNodes;
    sprintf(cbuf,"%-80s","coordinates"); // convex polyhedron format...
    if (mpi_rank == 0) {
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      MPI_File_write_at(fh,offset+80*sizeof(char),&nn,1,MPI_INT,MPI_STATUS_IGNORE);
    }
    offset += 80*sizeof(char)+int_size*1;

    FOR_I3 {
      int ind = 0;
      FOR_INO {
        if (no_flag[ino] >= 0) {
          fbuf[ind++] = (float)x_no[ino][i]; // x[N],y[N],z[N]
        }
      }
      assert(ind == myNumNodes);

      writeChunkedData<float>(fh,offset+(i*NumNodes+myNodeDisp)*sizeof(float),fbuf,myNumNodes,mpi_comm);
    }
    offset += sizeof(float)*3*NumNodes;

    // ==============================
    // connectivity...
    // ==============================

    // recall that the 0-indexed ino_global is in no_flag in only the
    // nodes that actually OWN the nodes. otherwise, there is a -1. Same for
    // faces, so we need to have temporary

    // also update the global node counts to
    // overwrite the -1 values on processors that did not
    // write data...

    int no_minus1_check = 0;
    int no_minus2_check = 0;
    FOR_INO {
      if (no_flag[ino] >= 0) {
        assert((no_flag[ino] >= myNodeDisp)&&(no_flag[ino] < myNodeDisp+myNumNodes));
      }
      else if (no_flag[ino] == -1) {
        ++no_minus1_check;
      }
      else {
        assert(no_flag[ino] == -2);
        ++no_minus2_check;
      }
    }

    updateNoData(no_flag,MAX_NO_PERIODIC_DATA);

    // should have gotten rid of all -1's...
    FOR_INO assert(no_flag[ino] != -1);

    sprintf(cbuf,"%-80s","nfaced"); // convex polyhedron format...
    int ne = (int)NumCells;
    if (mpi_rank == 0) {
      // write out element block header...
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      MPI_File_write_at(fh,offset+80*sizeof(char),&ne,1,MPI_INT,MPI_STATUS_IGNORE);
    }
    offset += 80*sizeof(char)+int_size*1;

    int *int_buf = new int[myNumCells]; // num faces for cell
    int *int_buf_2 = new int[myNumCellFaces]; // num nodes for face
    int *int_buf_3 = new int[myNumCellFaceNodes]; // node index for face node
    {
      int ind = 0, ind_2 = 0, ind_3 = 0;
      FOR_ICV {
        // for flagged cells...
        if (cv_flag[icv] != 0) {
          int_buf[ind++] = faocv_i[icv+1]-faocv_i[icv] + bfocv_i[icv+1]-bfocv_i[icv]; // cell faces
          // loop on internal faces...
          for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
            const int ifa = faocv_v[foc];
            const int nof_f = noofa_i[ifa];
            const int nof_l = noofa_i[ifa+1]-1;
            const int nnof = nof_l - nof_f + 1;
            int_buf_2[ind_2++] = nnof; // cell face nodes
            // loop over nodes of face....
            for (int nof = nof_f ; nof <= nof_l; nof++)
              int_buf_3[ind_3++] = (int)no_flag[noofa_v[nof]]+1; // node index
          }
          // loop on boundary faces...
          for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
            const int ibf = bfocv_v[boc];
            const int nob_f = noobf_i[ibf];
            const int nob_l = noobf_i[ibf+1]-1;
            const int nnob = nob_l - nob_f + 1;
            int_buf_2[ind_2++] = nnob; // cell boundary face nodes
            // loop over nodes of face...
            for (int nob = nob_f; nob <= nob_l; nob++)
              int_buf_3[ind_3++] = (int)no_flag[noobf_v[nob]]+1; // node index
          }
        }
      }
      assert(ind == myNumCells);
      assert(ind_2 == myNumCellFaces);
      assert(ind_3 == myNumCellFaceNodes);
    }

    writeChunkedData<int>(fh,offset+myCellDisp*int_size,int_buf,myNumCells,mpi_comm);
    delete[] int_buf;
    offset += int_size*NumCells;
    writeChunkedData<int>(fh,offset+myCellFaceDisp*int_size,int_buf_2,myNumCellFaces,mpi_comm);
    delete[] int_buf_2;
    offset += int_size*NumCellFaces;
    writeChunkedData<int>(fh,offset+myCellFaceNodeDisp*int_size,int_buf_3,myNumCellFaceNodes,mpi_comm);
    delete[] int_buf_3;
    offset += int_size*NumCellFaceNodes;

    // now return no_flag to only flagged nodes...

    int no_minus1_check2 = 0;
    int no_minus2_check2 = 0;
    FOR_INO {
      if (no_flag[ino] >= 0) {
        if ((no_flag[ino] < myNodeDisp)||(no_flag[ino] >= myNodeDisp+myNumNodes)) {
          no_flag[ino] = -1;
          ++no_minus1_check2;
        }
      }
      else {
        assert(no_flag[ino] == -2);
        ++no_minus2_check2;
      }
    }
    assert(no_minus1_check2 == no_minus1_check);
    assert(no_minus2_check2 == no_minus2_check);

    MPI_File_set_size(fh,offset);
    MPI_File_close(&fh);
  }

  // now lets write out the scalar and vector files for this step...
  COUT1("StaticSolver::writeFlaggedCvsEnsightVarFiles()");

  if (dw->nscalars_cv > 0) {
    double *no_data = new double[nno];
    for (int isc = 0; isc < dw->nscalars_cv; ++isc) {

      // open file...
      sprintf(cbuf,"%s.%08d.sca%02d",dw->name.c_str(),step,isc); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);
      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->scalar_descs_cv[isc].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
        sprintf(cbuf,"%-80s","part");
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
        int part = 1; // 1 part for now
        MPI_File_write_at(fh,offset,&part,1,MPI_INT,MPI_STATUS_IGNORE);
        offset += int_size;
        sprintf(cbuf,"%-80s","coordinates");
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
      }
      offset = int_size+3*80*sizeof(char);

      // write data...
      setNoDN(no_data,dw->scalar_names_cv[isc]);
      {
        int ind = 0;
        FOR_INO {
          if (no_flag[ino] >= 0) {
            fbuf[ind++] = (float)no_data[ino];
          }
        }
        assert(ind == myNumNodes);
      }
      writeChunkedData<float>(fh,offset+myNodeDisp*sizeof(float),fbuf,myNumNodes,mpi_comm);
      offset += sizeof(float)*NumNodes;

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);
    }
    delete[] no_data;
  }
  if (dw->nvectors_cv > 0) {
    double (*no_data)[3] = new double[nno][3];
    for (int ivec = 0; ivec < dw->nvectors_cv; ++ivec) {

      // open file...
      sprintf(cbuf,"%s.%08d.vec%02d",dw->name.c_str(),step,ivec); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);
      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->vector_descs_cv[ivec].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
        sprintf(cbuf,"%-80s","part");
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
        int part = 1; // 1 part for now
        MPI_File_write_at(fh,offset,&part,1,MPI_INT,MPI_STATUS_IGNORE);
        offset += int_size;
        sprintf(cbuf,"%-80s","coordinates");
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
      }
      offset = int_size+3*80*sizeof(char);

      // write data...
      setNoDN3(no_data,dw->vector_names_cv[ivec]);
      FOR_I3 {
        int ind = 0;
        FOR_INO {
          if (no_flag[ino] >= 0) {
            fbuf[ind++] = (float)no_data[ino][i];
          }
        }
        assert(ind == myNumNodes);
        writeChunkedData<float>(fh,offset+(i*NumNodes+myNodeDisp)*sizeof(float),fbuf,myNumNodes,mpi_comm);
      }
      offset += sizeof(float)*3*NumNodes;

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);
    }
    delete[] no_data;
  }

  // cleanup...
  delete[] cbuf;
  delete[] fbuf;
  delete[] no_flag;
}

void StaticSolver::writeSimpleTriVecEnsight(DataWriter* dw,const vector<SimpleTriWithWgts>& triVec,
    const int nst,const int (*spost)[3],const int nsp,const double *dn_no,
    const double *dn_tmp,const double (*dn3_no)[3],const double (*dn3_tmp)[3],const int nno,
    const int step,const double time) {

  // update case file with info about this step...
  updateEnsightCaseFile(dw,step,time);

  // get counts...
  int my_count[2] = {nsp,nst}; // my nodes/elements
  int my_disp[2]; // nodal, element displacements..
  MPI_Scan(my_count,my_disp,2,MPI_INT,MPI_SUM,mpi_comm);
  FOR_I2 my_disp[i] -= my_count[i];
  int count[2]; // everyone gets a copy of the node,element counts
  MPI_Allreduce(my_count,count,2,MPI_INT,MPI_SUM,mpi_comm);

  if ((count[0] > TWO_BILLION)||(count[1] > TWO_BILLION)) {
    CERR(" > ensight io cannot support >2B node/element count currently");
  }

  int myNumNodes = nsp;
  int myNumCells = nst;
  //int myNumCellNodes = 3*nst;
  int myNodeDisp = my_disp[0];
  int myCellDisp = my_disp[1];
  //int myCellNodeDisp = 3*my_disp[1];
  int NumNodes = count[0];
  int NumCells = count[1];
  //int NumCellNodes = 3*count[1];

  if (mpi_rank == 0)
    cout << " > global counts: " << NumNodes << " " << NumCells << endl;

  // geom and scalar files need these array, so just allocate here...
  char * cbuf = new char[256];
  float * fbuf = new float[myNumNodes];

  // for static solver the mesh is static, so we only need to build it once.
  // This assumes updateEnsightCasFile is called first!
  if (dw->ensight_casefile_nsteps == 1) {

    // lets build the geo file used for all transient data associated to this ensight case...
    COUT1("StaticSolver::writeSimpleTriVecEnsightGeoFile()");

    // open the file...

    string filename = dw->name + ".geo";
    MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...
    MPI_File_delete((char*)filename.c_str(),MPI_INFO_NULL);
    MPI_File fh;
    MPI_File_open(mpi_comm,(char*)filename.c_str(),MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    // ==============================
    // header...
    // ==============================

    MPI_Offset offset = 0;
    if (mpi_rank == 0) {
      // write out header

      sprintf(cbuf,"%-80s","C Binary");
      MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Geometry file"); // should probably include hashes here
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Produced by Cascade Technologies, Inc.");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","node id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","element id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","part");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      int part = 1; // 1 part for now
      MPI_File_write_at(fh,offset,&part,1,MPI_INT,MPI_STATUS_IGNORE);
      offset += int_size;
      sprintf(cbuf,"%-80s","Triangulated surface");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
    }
    offset = int_size+7*80*sizeof(char);

    // ==============================
    // coordinates...
    // ==============================

    sprintf(cbuf,"%-80s","coordinates"); // convex polyhedron format...
    if (mpi_rank == 0) {
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      MPI_File_write_at(fh,offset+80*sizeof(char),&NumNodes,1,MPI_INT,MPI_STATUS_IGNORE);
    }
    offset += 80*sizeof(char)+int_size*1;

    FOR_I3 {
      for (int ist = 0; ist < nst; ++ist) {
        // just overwrite the same data
        fbuf[spost[ist][0]] = (float)triVec[ist].x0[i];
        fbuf[spost[ist][1]] = (float)triVec[ist].x1[i];
        fbuf[spost[ist][2]] = (float)triVec[ist].x2[i];
      }
      writeChunkedData<float>(fh,offset+(i*NumNodes+myNodeDisp)*sizeof(float),fbuf,myNumNodes,mpi_comm);
    }
    offset += sizeof(float)*3*NumNodes;

    // ==============================
    // connectivity...
    // ==============================

    sprintf(cbuf,"%-80s","tria3"); // 3-node triangle...
    if (mpi_rank == 0) {
      // write out element block header...
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      MPI_File_write_at(fh,offset+80*sizeof(char),&NumCells,1,MPI_INT,MPI_STATUS_IGNORE);
    }
    offset += 80*sizeof(char)+int_size*1;

    int *ibuf = new int[3*nst];
    for (int ist = 0; ist < nst; ++ist) FOR_I3 ibuf[ist*3+i] = 1+spost[ist][i]+myNodeDisp; // 1-indexed
    writeChunkedData<int>(fh,offset+3*myCellDisp*int_size,ibuf,3*myNumCells,mpi_comm);
    delete[] ibuf;
    offset += int_size*NumCells*3;

    MPI_File_set_size(fh,offset);
    MPI_File_close(&fh);
  }

  // now lets write out the scalar and vector files for this step...
  COUT1("StaticSolver::writeSimpleTriVecEnsightVarFiles()");

  if (dw->nscalars_cv > 0) {
    for (int isc = 0; isc < dw->nscalars_cv; ++isc) {

      // open file...
      sprintf(cbuf,"%s.%08d.sca%02d",dw->name.c_str(),step,isc); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);
      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->scalar_descs_cv[isc].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
        sprintf(cbuf,"%-80s","part");
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
        int part = 1; // 1 part for now
        MPI_File_write_at(fh,offset,&part,1,MPI_INT,MPI_STATUS_IGNORE);
        offset += int_size;
        sprintf(cbuf,"%-80s","coordinates");
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
      }
      offset = int_size+3*80*sizeof(char);

      // write data...
      for (int itri = 0; itri < nst; ++itri) {
        // d0
        int i_l = triVec[itri].i0_l;
        int i_r = triVec[itri].i0_r;
        double wgt_l = triVec[itri].wgt0_l;
        double wgt_r = 1.0-wgt_l;
        double d_l, d_r;
        if (i_l < nno)
          d_l = wgt_l*dn_no[i_l*dw->nscalars_cv+isc];
        else
          d_l = wgt_l*dn_tmp[(i_l-nno)*dw->nscalars_cv+isc];
        if (i_r < nno)
          d_r = wgt_r*dn_no[i_r*dw->nscalars_cv+isc];
        else
          d_r = wgt_r*dn_tmp[(i_r-nno)*dw->nscalars_cv+isc];
        fbuf[spost[itri][0]] = (float)(d_l+d_r);
        // d1
        i_l = triVec[itri].i1_l;
        i_r = triVec[itri].i1_r;
        wgt_l = triVec[itri].wgt1_l;
        wgt_r = 1.0-wgt_l;
        if (i_l < nno)
          d_l = wgt_l*dn_no[i_l*dw->nscalars_cv+isc];
        else
          d_l = wgt_l*dn_tmp[(i_l-nno)*dw->nscalars_cv+isc];
        if (i_r < nno)
          d_r = wgt_r*dn_no[i_r*dw->nscalars_cv+isc];
        else
          d_r = wgt_r*dn_tmp[(i_r-nno)*dw->nscalars_cv+isc];
        fbuf[spost[itri][1]] = (float)(d_l+d_r);
        // d2
        i_l = triVec[itri].i2_l;
        i_r = triVec[itri].i2_r;
        wgt_l = triVec[itri].wgt2_l;
        wgt_r = 1.0-wgt_l;
        if (i_l < nno)
          d_l = wgt_l*dn_no[i_l*dw->nscalars_cv+isc];
        else
          d_l = wgt_l*dn_tmp[(i_l-nno)*dw->nscalars_cv+isc];
        if (i_r < nno)
          d_r = wgt_r*dn_no[i_r*dw->nscalars_cv+isc];
        else
          d_r = wgt_r*dn_tmp[(i_r-nno)*dw->nscalars_cv+isc];
        fbuf[spost[itri][2]] = (float)(d_l+d_r);
      }
      writeChunkedData<float>(fh,offset+myNodeDisp*sizeof(float),fbuf,myNumNodes,mpi_comm);
      offset += sizeof(float)*NumNodes;

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);
    }
  }
  if (dw->nvectors_cv > 0) {
    for (int ivec = 0; ivec < dw->nvectors_cv; ++ivec) {

      // open file...
      sprintf(cbuf,"%s.%08d.vec%02d",dw->name.c_str(),step,ivec); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);
      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->vector_descs_cv[ivec].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
        sprintf(cbuf,"%-80s","part");
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
        int part = 1; // 1 part for now
        MPI_File_write_at(fh,offset,&part,1,MPI_INT,MPI_STATUS_IGNORE);
        offset += int_size;
        sprintf(cbuf,"%-80s","coordinates");
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        offset += 80*sizeof(char);
      }
      offset = int_size+3*80*sizeof(char);

      // write data...
      FOR_I3 {
        for (int itri = 0; itri < nst; ++itri) {
          // d0
          int i_l = triVec[itri].i0_l;
          int i_r = triVec[itri].i0_r;
          double wgt_l = triVec[itri].wgt0_l;
          double wgt_r = 1.0-wgt_l;
          double d_l, d_r;
          if (i_l < nno)
            d_l = wgt_l*dn3_no[i_l*dw->nvectors_cv+ivec][i];
          else
            d_l = wgt_l*dn3_tmp[(i_l-nno)*dw->nvectors_cv+ivec][i];
          if (i_r < nno)
            d_r = wgt_r*dn3_no[i_r*dw->nvectors_cv+ivec][i];
          else
            d_r = wgt_r*dn3_tmp[(i_r-nno)*dw->nvectors_cv+ivec][i];
          fbuf[spost[itri][0]] = (float)(d_l+d_r);
          // d1
          i_l = triVec[itri].i1_l;
          i_r = triVec[itri].i1_r;
          wgt_l = triVec[itri].wgt1_l;
          wgt_r = 1.0-wgt_l;
          if (i_l < nno)
            d_l = wgt_l*dn3_no[i_l*dw->nvectors_cv+ivec][i];
          else
            d_l = wgt_l*dn3_tmp[(i_l-nno)*dw->nvectors_cv+ivec][i];
          if (i_r < nno)
            d_r = wgt_r*dn3_no[i_r*dw->nvectors_cv+ivec][i];
          else
            d_r = wgt_r*dn3_tmp[(i_r-nno)*dw->nvectors_cv+ivec][i];
          fbuf[spost[itri][1]] = (float)(d_l+d_r);
          // d2
          i_l = triVec[itri].i2_l;
          i_r = triVec[itri].i2_r;
          wgt_l = triVec[itri].wgt2_l;
          wgt_r = 1.0-wgt_l;
          if (i_l < nno)
            d_l = wgt_l*dn3_no[i_l*dw->nvectors_cv+ivec][i];
          else
            d_l = wgt_l*dn3_tmp[(i_l-nno)*dw->nvectors_cv+ivec][i];
          if (i_r < nno)
            d_r = wgt_r*dn3_no[i_r*dw->nvectors_cv+ivec][i];
          else
            d_r = wgt_r*dn3_tmp[(i_r-nno)*dw->nvectors_cv+ivec][i];
          fbuf[spost[itri][2]] = (float)(d_l+d_r);
        }
        writeChunkedData<float>(fh,offset+(i*NumNodes+myNodeDisp)*sizeof(float),fbuf,myNumNodes,mpi_comm);
      }
      offset += sizeof(float)*3*NumNodes;

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);
    }
  }

  // cleanup...
  delete[] cbuf;
  delete[] fbuf;
}

void StaticSolver::writeFlaggedBfZonesDualEnsight(DataWriter* dw,const int step,const double time) {

  // update case file with info about this step...
  updateEnsightCaseFile(dw,step,time);

  vector<int> znofz;
  for (map<const string,int>::const_iterator it = bfZoneNameMap.begin(); it != bfZoneNameMap.end(); ++it) {
    const int iz = it->second; // zone index
    if (dw->bfzone_flag[iz]) 
      znofz.push_back(iz);
  }
  const int nfz = znofz.size();

  int* ssp_flag_fz = new int[subSurface->nsp*nfz];
  for (int ifz = 0; ifz < nfz; ++ifz) {
    int* const ssp_flag = &ssp_flag_fz[ifz*subSurface->nsp];
    for (int issp = 0; issp < subSurface->nsp; ++issp)
      ssp_flag[issp] = mpi_size;
    const int iz = znofz[ifz];
    for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
      assert(zone_bf[ibf] == iz);
      for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
        const int issp = sspobf_v[pob];
        assert((issp >= 0)&&(issp < subSurface->nsp));
        ssp_flag[issp] = mpi_rank;
      }
    }
    subSurface->updateSpData(ssp_flag,MIN_DATA); // smallest rank wins
  }

  // points b/w zones are multiply written if the proc own it

  int *my_sp_count = new int[nfz];
  for (int ifz = 0; ifz < nfz; ++ifz) {
    const int iz = znofz[ifz];
    my_sp_count[ifz] = 0;
    int* const ssp_flag = &ssp_flag_fz[ifz*subSurface->nsp];
    for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
      assert(zone_bf[ibf] == iz);
      for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
        const int issp = sspobf_v[pob];
        assert((issp >= 0)&&(issp < subSurface->nsp));
        if (ssp_flag[issp] == mpi_rank) {
          ++my_sp_count[ifz]; // surface point dumped by us
          ssp_flag[issp] = -1;
        }
      }
    }
    // reset...
    for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
      assert(zone_bf[ibf] == iz);
      for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
        const int issp = sspobf_v[pob];
        assert((issp >= 0)&&(issp < subSurface->nsp));
        if (ssp_flag[issp] == -1)
          ssp_flag[issp] = mpi_rank;
        assert((ssp_flag[issp] >= 0)&&(ssp_flag[issp] <= mpi_size));
      }
    }
  }

  // figure out our nodal and element displacement...

  int *my_sp_disp = new int[nfz];
  MPI_Scan(my_sp_count,my_sp_disp,nfz,MPI_INT,MPI_SUM,mpi_comm);
  for (int ifz = 0; ifz < nfz; ++ifz)
    my_sp_disp[ifz] -= my_sp_count[ifz];

  // everyone gets the counts...

  int *sp_count = new int[nfz];
  MPI_Allreduce(my_sp_count,sp_count,nfz,MPI_INT,MPI_SUM,mpi_comm);

  int max_my_sp_count = 0;
  for (int ifz = 0; ifz < nfz; ++ifz)
    max_my_sp_count = max(max_my_sp_count,my_sp_count[ifz]);

  // geom and scalar files need these array, so just allocate here...
  char * cbuf = new char[256];
  float * fbuf = new float[max_my_sp_count];

  // for static solver the mesh is static, so we only need to build it once.
  // This assumes updateEnsightCasFile is called first!
  if (dw->ensight_casefile_nsteps == 1) {

    int* sst_flag_fz = new int[subSurface->nst*nfz];
    for (int ifz = 0; ifz < nfz; ++ifz) {
      int* const sst_flag = &sst_flag_fz[ifz*subSurface->nst];
      for (int isst = 0; isst < subSurface->nst; ++isst)
        sst_flag[isst] = mpi_size;
      const int iz = znofz[ifz];
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        for (int tob = sstobf_i[ibf]; tob != sstobf_i[ibf+1]; ++tob) {
          const int isst = sstobf_v[tob];
          assert((isst >= 0)&&(isst < subSurface->nst));
          sst_flag[isst] = mpi_rank;
        }
      }
      subSurface->updateStData(sst_flag,MIN_DATA); // smallest rank wins
    }

    // points b/w zones are multiply written if the proc own it

    int *my_st_count = new int[nfz]; 
    for (int ifz = 0; ifz < nfz; ++ifz) {
      const int iz = znofz[ifz];
      int* const sst_flag = &sst_flag_fz[ifz*subSurface->nst];
      my_st_count[ifz] = 0;
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        for (int tob = sstobf_i[ibf]; tob != sstobf_i[ibf+1]; ++tob) {
          const int isst = sstobf_v[tob];
          assert((isst >= 0)&&(isst < subSurface->nst));
          assert(dw->bfzone_flag[subSurface->znost[isst]]);
          if (sst_flag[isst] == mpi_rank) {
            ++my_st_count[ifz]; // surface tri dumped by us
            sst_flag[isst] = -1;
          }
        }
      }
      // reset...
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        for (int tob = sstobf_i[ibf]; tob != sstobf_i[ibf+1]; ++tob) {
          const int isst = sstobf_v[tob];
          assert((isst >= 0)&&(isst < subSurface->nst));
          if (sst_flag[isst] == -1)
            sst_flag[isst] = mpi_rank;
          assert((sst_flag[isst] >= 0)&&(sst_flag[isst] <= mpi_size));
        }
      }
    }

    // figure out our nodal and element displacement...

    int *my_st_disp = new int[nfz];
    MPI_Scan(my_st_count,my_st_disp,nfz,MPI_INT,MPI_SUM,mpi_comm);
    for (int ifz = 0; ifz < nfz; ++ifz)
      my_st_disp[ifz] -= my_st_count[ifz];

    // everyone gets the counts...

    int *st_count = new int[nfz];
    MPI_Allreduce(my_st_count,st_count,nfz,MPI_INT,MPI_SUM,mpi_comm);

    int max_my_st_count = 0;
    for (int ifz = 0; ifz < nfz; ++ifz)
      max_my_st_count = max(max_my_st_count,my_st_count[ifz]);

    // lets build the geo file used for all transient data associated to this ensight case...
    if (mpi_rank == 0)
      cout << "StaticSolver::writeFlaggedBfZonesDualEnsightGeoFile()" << endl;

    // open the file...

    string filename = dw->name + ".geo";
    MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...
    MPI_File_delete((char*)filename.c_str(),MPI_INFO_NULL);
    MPI_File fh;
    MPI_File_open(mpi_comm,(char*)filename.c_str(),MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    // ==============================
    // header...
    // ==============================

    MPI_Offset offset = 0;
    if (mpi_rank == 0) {
      // write out header
      sprintf(cbuf,"%-80s","C Binary");
      MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Geometry file"); // should probably include hashes here
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Produced by Cascade Technologies, Inc.");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","node id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","element id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
    }
    offset = 5*80*sizeof(char);

    // go through all flagged bf zones...

    int *ibuf = new int[3*max_my_st_count]; // global point indices for each tri
    int* ssp_flag2 = new int[subSurface->nsp];
    for (int ifz = 0; ifz < nfz; ++ifz) {
      const int iz = znofz[ifz];
      int* const sst_flag = &sst_flag_fz[ifz*subSurface->nst];
      int* const ssp_flag = &ssp_flag_fz[ifz*subSurface->nsp];

      // ==============================
      // part header...
      // ==============================

      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s","part");
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        int part = 1+ifz;
        MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
        sprintf(cbuf,"%-80s",bfZoneVec[iz].getName().c_str());
        MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      }
      offset += int_size+2*80*sizeof(char);

      // ==============================
      // part coordinates...
      // ==============================

      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s","coordinates"); // convex polyhedron format...
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        MPI_File_write_at(fh,offset+80*sizeof(char),&sp_count[ifz],1,MPI_INT,MPI_STATUS_IGNORE);
      }
      offset += 80*sizeof(char)+int_size*1;

      // write out x[3]...
      FOR_I3 {
        int ind = 0;
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);
          for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
            const int issp = sspobf_v[pob];
            assert((issp >= 0)&&(issp < subSurface->nsp));
            if (ssp_flag[issp] == mpi_rank) {
              ssp_flag[issp] = -1;
              fbuf[ind++] = (float)subSurface->xp[issp][i]; // x[N],y[N],z[N]
            }
          }
        }
        // reset...
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);
          for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
            const int issp = sspobf_v[pob];
            assert((issp >= 0)&&(issp < subSurface->nsp));
            if (ssp_flag[issp] == -1)
              ssp_flag[issp] = mpi_rank;
            assert((ssp_flag[issp] >= 0)&&(ssp_flag[issp] <= mpi_size));
          }
        }
        assert(ind == my_sp_count[ifz]);

        writeChunkedData<float>(fh,offset+(i*sp_count[ifz]+my_sp_disp[ifz])*sizeof(float),fbuf,my_sp_count[ifz],mpi_comm);
      }
      offset += sizeof(float)*3*sp_count[ifz];

      // ==============================
      // connectivity...
      // ==============================

      sprintf(cbuf,"%-80s","tria3"); // triangles
      if (mpi_rank == 0) {
        // write out element block header...
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        MPI_File_write_at(fh,offset+80*sizeof(char),&st_count[ifz],1,MPI_INT,MPI_STATUS_IGNORE);
      }
      offset += 80*sizeof(char)+int_size*1;

      // global index...
      int my_count0 = 0;
      for (int issp = 0; issp < subSurface->nsp; ++issp)
        ssp_flag2[issp] = -2;
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
          const int issp = sspobf_v[pob];
          assert((issp >= 0)&&(issp < subSurface->nsp));
          // first time this node was visited by this zone
          if (ssp_flag[issp] >= 0) {
            assert(ssp_flag[issp] <= mpi_size);
            if (ssp_flag[issp] == mpi_rank) {
              ssp_flag2[issp] = my_sp_disp[ifz]+my_count0++; // this node is dumped by us
            }
            else if (ssp_flag[issp] < mpi_size) {
              ssp_flag2[issp] = -1; // this node gets dumped, but not by us
            }
            else {
              assert(ssp_flag[issp] == mpi_size);
              ssp_flag2[issp] = -2; // nobody owns this node
            }
            ssp_flag[issp] = -ssp_flag[issp]-1;
          }
        }
      }
      // reset...
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
          const int issp = sspobf_v[pob];
          assert((issp >= 0)&&(issp < subSurface->nsp));
          if (ssp_flag[issp] < 0)
            ssp_flag[issp] = -ssp_flag[issp]-1;
          assert((ssp_flag[issp] >= 0)&&(ssp_flag[issp] <= mpi_size));
        }
      }
      assert(my_count0 == my_sp_count[ifz]);
      subSurface->updateSpData(ssp_flag2,MAX_DATA);

      int ind = 0;
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        for (int tob = sstobf_i[ibf]; tob != sstobf_i[ibf+1]; ++tob) {
          const int isst = sstobf_v[tob];
          assert((isst >= 0)&&(isst < subSurface->nst));
          if (sst_flag[isst] == mpi_rank) {
            sst_flag[isst] = -1;
            FOR_I3 {
              const int issp = subSurface->spost[isst][i];
              assert((issp >= 0)&&(issp < subSurface->nsp));
              assert(ssp_flag2[issp] >= 0);
              ibuf[ind++] = ssp_flag2[issp]+1; // 1-indexed
            }
          }
        }
      }
      assert(ind == 3*my_st_count[ifz]);
      // reset...
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        for (int tob = sstobf_i[ibf]; tob != sstobf_i[ibf+1]; ++tob) {
          const int isst = sstobf_v[tob];
          assert((isst >= 0)&&(isst < subSurface->nst));
          if (sst_flag[isst] == -1)
            sst_flag[isst] = mpi_rank;
          assert((sst_flag[isst] >= 0)&&(sst_flag[isst] <= mpi_size));
        }
      }

      writeChunkedData<int>(fh,offset+my_st_disp[ifz]*3*int_size,ibuf,3*my_st_count[ifz],mpi_comm);
      offset += int_size*3*st_count[ifz];

    }
    delete[] ibuf;
    delete[] ssp_flag2;
    delete[] sst_flag_fz;
    delete[] st_count;
    delete[] my_st_count;
    delete[] my_st_disp;

    MPI_File_set_size(fh,offset);
    MPI_File_close(&fh);

  }

  // now lets write out the scalar and vector files for this step...
  if (mpi_rank == 0)
    cout << "StaticSolver::writeFlaggedBfZonesDualEnsightVarFiles()" << endl;

  if ((dw->nscalars_cv + dw->nscalars_bf + dw->nvectors_cv + dw->nvectors_bf) > 0) {

    double *inv_ssp_wgt_sum = new double[subSurface->nsp];
    for (int issp = 0; issp < subSurface->nsp; ++issp)
      inv_ssp_wgt_sum[issp] = 0.0;
    for (int ifz = 0; ifz < nfz; ++ifz) {
      const int iz = znofz[ifz];
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
          const int issp = sspobf_v[pob];
          inv_ssp_wgt_sum[issp] += sspobf_wgt[pob];
        }
      }
    }
    subSurface->updateSpData(inv_ssp_wgt_sum,ADD_DATA);
    for (int issp = 0; issp < subSurface->nsp; ++issp) {
      if (inv_ssp_wgt_sum[issp] > 0.0)
        inv_ssp_wgt_sum[issp] = 1.0/inv_ssp_wgt_sum[issp];
    }
    dumpRange(inv_ssp_wgt_sum,subSurface->nsp,"inv_ssp_wgt_sum");

    if (dw->nscalars_cv > 0||dw->nscalars_bf > 0) {
      double *ssp_data = new double[subSurface->nsp];
      for (int isc = 0; isc < dw->nscalars_cv; ++isc) {

        // open file...
        sprintf(cbuf,"%s.%08d.sca%02d",dw->name.c_str(),step,isc); // see case file
        MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
        MPI_File_delete(cbuf,MPI_INFO_NULL);
        MPI_File fh;
        MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

        // write header...
        MPI_Offset offset = 0;
        if (mpi_rank == 0) {
          sprintf(cbuf,"%-80s",dw->scalar_descs_cv[isc].c_str()); // description, hashes??
          MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        }
        offset += 80*sizeof(char);

        // get data on subsurface points...
        for (int issp = 0; issp < subSurface->nsp; ++issp)
          ssp_data[issp] = 0.0;
        CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->scalar_names_cv[isc]);
        assert((var->getType() == DN_DATA)&&(var->getUnindexedTopology() == CV_DATA));
        double* var_d = var->getDNptr();
        for (int ifz = 0; ifz < nfz; ++ifz) {
          const int iz = znofz[ifz];
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            const int icv = cvobf[ibf];
            for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
              const int issp = sspobf_v[pob];
              ssp_data[issp] += sspobf_wgt[pob]*var_d[icv];
            }
          }
        }
        subSurface->updateSpData(ssp_data,ADD_DATA);
        for (int issp = 0; issp < subSurface->nsp; ++issp)
          ssp_data[issp] *= inv_ssp_wgt_sum[issp];
        for (int ifz = 0; ifz < nfz; ++ifz) {
          const int iz = znofz[ifz];
          int* const ssp_flag = &ssp_flag_fz[ifz*subSurface->nsp];

          if (mpi_rank == 0) {
            sprintf(cbuf,"%-80s","part");
            MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
            int part = 1+ifz;
            MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
            sprintf(cbuf,"%-80s","coordinates");
            MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          }
          offset += int_size+2*80*sizeof(char);

          // write data...
          int ind = 0;
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
              const int issp = sspobf_v[pob];
              assert((issp >= 0)&&(issp < subSurface->nsp));
              if (ssp_flag[issp] == mpi_rank) {
                ssp_flag[issp] = -1;
                fbuf[ind++] = (float)ssp_data[issp];
              }
            }
          }
          // reset...
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
              const int issp = sspobf_v[pob];
              assert((issp >= 0)&&(issp < subSurface->nsp));
              if (ssp_flag[issp] == -1)
                ssp_flag[issp] = mpi_rank;
              assert((ssp_flag[issp] >= 0)&&(ssp_flag[issp] <= mpi_size));
            }
          }
          assert(ind == my_sp_count[ifz]);

          writeChunkedData<float>(fh,offset+my_sp_disp[ifz]*sizeof(float),fbuf,my_sp_count[ifz],mpi_comm);
          offset += sizeof(float)*sp_count[ifz];
        }

        // close file...
        MPI_File_set_size(fh,offset);
        MPI_File_close(&fh);
      }
      if (dw->nscalars_bf > 0) {
        for (int isc = 0; isc < dw->nscalars_bf/dw->nzones_bf; ++isc) {

          // open file...
          sprintf(cbuf,"%s.%08d.sca%02d",dw->name.c_str(),step,isc+dw->nscalars_cv); // see case file
          MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
          MPI_File_delete(cbuf,MPI_INFO_NULL);
          MPI_File fh;
          MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

          // write header...
          MPI_Offset offset = 0;
          if (mpi_rank == 0) {
            sprintf(cbuf,"%-80s",dw->scalar_descs_bf[isc].c_str()); // description, hashes??
            MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          }
          offset += 80*sizeof(char);

          for (int issp = 0; issp < subSurface->nsp; ++issp)
            ssp_data[issp] = 0.0;
          for (int izone = 0; izone < dw->nzones_bf; ++izone) {
            CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->scalar_names_bf[izone*(dw->nscalars_bf/dw->nzones_bf)+isc]);
            const int iz = var->getIndex();
            assert((var->getType() == DN_DATA)&&(var->getUnindexedTopology() == BF_DATA));
            double* var_d = var->getDNptr();
            for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
              assert(zone_bf[ibf] == iz);
              for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
                const int issp = sspobf_v[pob];
                ssp_data[issp] += sspobf_wgt[pob]*var_d[ibf-bfZoneVec[iz].ibf_f];
              }
            }
          }
          subSurface->updateSpData(ssp_data,ADD_DATA);
          for (int issp = 0; issp < subSurface->nsp; ++issp)
            ssp_data[issp] *= inv_ssp_wgt_sum[issp];

          for (int ifz = 0; ifz < nfz; ++ifz) {
            const int iz = znofz[ifz];
            int* const ssp_flag = &ssp_flag_fz[ifz*subSurface->nsp];

            if (mpi_rank == 0) {
              sprintf(cbuf,"%-80s","part");
              MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
              int part = 1+ifz;
              MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
              sprintf(cbuf,"%-80s","coordinates");
              MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
            }
            offset += int_size+2*80*sizeof(char);

            // write data...
            int ind = 0;
            for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
              assert(zone_bf[ibf] == iz);
              for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
                const int issp = sspobf_v[pob];
                assert((issp >= 0)&&(issp < subSurface->nsp));
                if (ssp_flag[issp] == mpi_rank) {
                  ssp_flag[issp] = -1;
                  fbuf[ind++] = (float)ssp_data[issp];
                }
              }
            }
            // reset...
            for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
              assert(zone_bf[ibf] == iz);
              for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
                const int issp = sspobf_v[pob];
                assert((issp >= 0)&&(issp < subSurface->nsp));
                if (ssp_flag[issp] == -1)
                  ssp_flag[issp] = mpi_rank;
                assert((ssp_flag[issp] >= 0)&&(ssp_flag[issp] <= mpi_size));
              }
            }
            assert(ind == my_sp_count[ifz]);

            writeChunkedData<float>(fh,offset+my_sp_disp[ifz]*sizeof(float),fbuf,my_sp_count[ifz],mpi_comm);
            offset += sizeof(float)*sp_count[ifz];
          }

          // close file...
          MPI_File_set_size(fh,offset);
          MPI_File_close(&fh);
        }
      }
      delete[] ssp_data;
    }
    if (dw->nvectors_cv > 0 || dw->nvectors_bf > 0) {
      double (*ssp_data)[3] = new double[subSurface->nsp][3];
      for (int ivec = 0; ivec < dw->nvectors_cv; ++ivec) {

        // open file...
        sprintf(cbuf,"%s.%08d.vec%02d",dw->name.c_str(),step,ivec); // see case file
        MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
        MPI_File_delete(cbuf,MPI_INFO_NULL);
        MPI_File fh;
        MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

        // write header...
        MPI_Offset offset = 0;
        if (mpi_rank == 0) {
          sprintf(cbuf,"%-80s",dw->vector_descs_cv[ivec].c_str()); // description, hashes??
          MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        }
        offset += 80*sizeof(char);

        // get data on subsurface points...
        for (int issp = 0; issp < subSurface->nsp; ++issp)
          FOR_I3 ssp_data[issp][i] = 0.0;
        CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->vector_names_cv[ivec]);
        assert((var->getType() == DN3_DATA)&&(var->getUnindexedTopology() == CV_DATA));
        double (*var_d3)[3] = var->getDN3ptr();
        for (int ifz = 0; ifz < nfz; ++ifz) {
          const int iz = znofz[ifz];
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            const int icv = cvobf[ibf];
            for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
              const int issp = sspobf_v[pob];
              FOR_I3 ssp_data[issp][i] += sspobf_wgt[pob]*var_d3[icv][i];
            }
          }
        }
        subSurface->updateSpData(ssp_data,ADD_DATA);
        for (int issp = 0; issp < subSurface->nsp; ++issp)
          FOR_I3 ssp_data[issp][i] *= inv_ssp_wgt_sum[issp];

        for (int ifz = 0; ifz < nfz; ++ifz) {
          const int iz = znofz[ifz];
          int* const ssp_flag = &ssp_flag_fz[ifz*subSurface->nsp];

          if (mpi_rank == 0) {
            sprintf(cbuf,"%-80s","part");
            MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
            int part = 1+ifz;
            MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
            sprintf(cbuf,"%-80s","coordinates");
            MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          }
          offset += int_size+2*80*sizeof(char);

          // write data...
          FOR_I3 {
            int ind = 0;
            for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
              assert(zone_bf[ibf] == iz);
              for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
                const int issp = sspobf_v[pob];
                assert((issp >= 0)&&(issp < subSurface->nsp));
                if (ssp_flag[issp] == mpi_rank) {
                  ssp_flag[issp] = -1;
                  fbuf[ind++] = (float)ssp_data[issp][i];
                }
              }
            }
            // reset...
            for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
              assert(zone_bf[ibf] == iz);
              for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
                const int issp = sspobf_v[pob];
                assert((issp >= 0)&&(issp < subSurface->nsp));
                if (ssp_flag[issp] == -1)
                  ssp_flag[issp] = mpi_rank;
                assert((ssp_flag[issp] >= 0)&&(ssp_flag[issp] <= mpi_size));
              }
            }
            assert(ind == my_sp_count[ifz]);

            writeChunkedData<float>(fh,offset+(i*sp_count[ifz]+my_sp_disp[ifz])*sizeof(float),fbuf,my_sp_count[ifz],mpi_comm);
          }
          offset += sizeof(float)*3*sp_count[ifz];
        }

        // close file...
        MPI_File_set_size(fh,offset);
        MPI_File_close(&fh);
      }
      if (dw->nvectors_bf > 0) {
        for (int ivec = 0; ivec < dw->nvectors_bf/dw->nzones_bf; ++ivec) {

          // open file...
          sprintf(cbuf,"%s.%08d.vec%02d",dw->name.c_str(),step,ivec+dw->nvectors_cv); // see case file
          MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
          MPI_File_delete(cbuf,MPI_INFO_NULL);
          MPI_File fh;
          MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

          // write header...
          MPI_Offset offset = 0;
          if (mpi_rank == 0) {
            sprintf(cbuf,"%-80s",dw->vector_descs_bf[ivec].c_str()); // description, hashes??
            MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          }
          offset += 80*sizeof(char);

          for (int issp = 0; issp < subSurface->nsp; ++issp)
            FOR_I3 ssp_data[issp][i] = 0.0;
          for (int izone = 0; izone < dw->nzones_bf; ++izone) {
            CtiRegister::CtiData * var = CtiRegister::getCtiData(dw->vector_names_bf[izone*(dw->nvectors_bf/dw->nzones_bf)+ivec]);
            const int iz = var->getIndex();
            assert((var->getType() == DN3_DATA)&&(var->getUnindexedTopology() == BF_DATA));
            double (*var_d3)[3] = var->getDN3ptr();
            for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
              assert(zone_bf[ibf] == iz);
              for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
                const int issp = sspobf_v[pob];
                FOR_I3 ssp_data[issp][i] += sspobf_wgt[pob]*var_d3[ibf-bfZoneVec[iz].ibf_f][i];
              }
            }
          }
          subSurface->updateSpData(ssp_data,ADD_DATA);
          for (int issp = 0; issp < subSurface->nsp; ++issp)
            FOR_I3 ssp_data[issp][i] *= inv_ssp_wgt_sum[issp];

          for (int ifz = 0; ifz < nfz; ++ifz) {
            const int iz = znofz[ifz];
            int* const ssp_flag = &ssp_flag_fz[ifz*subSurface->nsp];

            if (mpi_rank == 0) {
              sprintf(cbuf,"%-80s","part");
              MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
              int part = 1+ifz;
              MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
              sprintf(cbuf,"%-80s","coordinates");
              MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
            }
            offset += int_size+2*80*sizeof(char);

            // write data...
            FOR_I3 {
              int ind = 0;
              for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
                assert(zone_bf[ibf] == iz);
                for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
                  const int issp = sspobf_v[pob];
                  assert((issp >= 0)&&(issp < subSurface->nsp));
                  if (ssp_flag[issp] == mpi_rank) {
                    ssp_flag[issp] = -1;
                    fbuf[ind++] = (float)ssp_data[issp][i];
                  }
                }
              }
              // reset...
              for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
                assert(zone_bf[ibf] == iz);
                for (int pob = sspobf_i[ibf]; pob != sspobf_i[ibf+1]; ++pob) {
                  const int issp = sspobf_v[pob];
                  assert((issp >= 0)&&(issp < subSurface->nsp));
                  if (ssp_flag[issp] == -1)
                    ssp_flag[issp] = mpi_rank;
                  assert((ssp_flag[issp] >= 0)&&(ssp_flag[issp] <= mpi_size));
                }
              }
              assert(ind == my_sp_count[ifz]);

              writeChunkedData<float>(fh,offset+(i*sp_count[ifz]+my_sp_disp[ifz])*sizeof(float),fbuf,my_sp_count[ifz],mpi_comm);
            }
            offset += sizeof(float)*3*sp_count[ifz];
          }

          // close file...
          MPI_File_set_size(fh,offset);
          MPI_File_close(&fh);
        }
      }
      delete[] ssp_data;
    }
    delete[] inv_ssp_wgt_sum;
  }

  // cleanup
  delete[] my_sp_count;
  delete[] my_sp_disp;
  delete[] sp_count;
  delete[] ssp_flag_fz;
  delete[] cbuf;
  delete[] fbuf;
}

void StaticSolver::writeFlaggedBfZonesEnsight(DataWriter* dw,const int step,const double time) {

  // update case file with info about this step...
  updateEnsightCaseFile(dw,step,time);

  int* no_flag = new int[nno];

  FOR_INO no_flag[ino] = mpi_size;

  int nfz = 0;
  for (map<const string,int>::const_iterator it = bfZoneNameMap.begin(); it != bfZoneNameMap.end(); ++it) {
    const int iz = it->second; // zone index
    if (dw->bfzone_flag[iz]) {
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int ino = noobf_v[nob];
          assert((ino >= 0)&&(ino < nno));
          no_flag[ino] = mpi_rank;
        }
        int nnof = noobf_i[ibf+1]-noobf_i[ibf];
        assert( nnof >= 3 ); //check later ...
      }
      nfz += 1;
    }
  }

  int *znofz = new int[nfz]; // zn of flagged zn
  nfz = 0;
  for (map<const string,int>::const_iterator it = bfZoneNameMap.begin(); it != bfZoneNameMap.end(); ++it) {
    const int iz = it->second; // zone index
    if (dw->bfzone_flag[iz]) {
      znofz[nfz++] = iz;
    }
  }

  updateNoData(no_flag,MIN_NO_PERIODIC_DATA); // smallest rank wins

  // nodes b/w zones are multiply written if the proc own it

  bool *visited_no = new bool[nno]; FOR_INO visited_no[ino] = false;
  int (*my_count)[3] = new int[nfz][3]; // local ncount,ecount and necount
  for (int ifz = 0; ifz < nfz; ++ifz) {
    const int iz = znofz[ifz];
    // get zone node,element,node-element counts...
    my_count[ifz][0] = 0;
    my_count[ifz][1] = bfZoneVec[iz].nbf;
    my_count[ifz][2] = 0;
    for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
      assert(zone_bf[ibf] == iz);
      my_count[ifz][2] += noobf_i[ibf+1]-noobf_i[ibf];
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino = noobf_v[nob];
        assert((ino >= 0)&&(ino < nno));
        // first time this node was visited by this zone
        if (!visited_no[ino]) {
          visited_no[ino] = true;
          if (no_flag[ino] == mpi_rank) {
            my_count[ifz][0] += 1; // this node is dumped by us
          }
        }
      }
    }
    // reset visited_no...
    for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
      assert(zone_bf[ibf] == iz);
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino = noobf_v[nob];
        assert((ino >= 0)&&(ino < nno));
        visited_no[ino] = false;
      }
    }
  }

  // figure out our nodal and element displacement...

  int (*my_disp)[3] = new int[nfz][3];
  MPI_Scan((int*)my_count,(int*)my_disp,3*nfz,MPI_INT,MPI_SUM,mpi_comm);
  for (int ifz = 0; ifz < nfz; ++ifz) {
    my_disp[ifz][0] -= my_count[ifz][0];
    my_disp[ifz][1] -= my_count[ifz][1];
    my_disp[ifz][2] -= my_count[ifz][2];
  }

  // everyone gets the counts...

  int (*count)[3] = new int[nfz][3];
  MPI_Allreduce((int*)my_count,(int*)count,3*nfz,MPI_INT,MPI_SUM,mpi_comm);

  int max_my_count[3] = {0,0,0};
  for (int ifz = 0; ifz < nfz; ++ifz) {
    FOR_I3 max_my_count[i] = max(max_my_count[i],my_count[ifz][i]);
  }

  // geom and scalar files need these array, so just allocate here...
  char * cbuf = new char[256];
  float * fbuf = new float[max_my_count[0]];

  // for static solver the mesh is static, so we only need to build it once.
  // This assumes updateEnsightCasFile is called first!
  if (dw->ensight_casefile_nsteps == 1) {

    // lets build the geo file used for all transient data associated to this ensight case...
    COUT1("StaticSolver::writeFlaggedBfZonesEnsightGeoFile()");

    // open the file...

    string filename = dw->name + ".geo";
    MiscUtils::mkdir_for_file_collective(filename.c_str()); // allow the user to have specified a subdirectory...
    MPI_File_delete((char*)filename.c_str(),MPI_INFO_NULL);
    MPI_File fh;
    MPI_File_open(mpi_comm,(char*)filename.c_str(),MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    // ==============================
    // header...
    // ==============================

    MPI_Offset offset = 0;
    if (mpi_rank == 0) {
      // write out header

      sprintf(cbuf,"%-80s","C Binary");
      MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Geometry file"); // should probably include hashes here
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Produced by Cascade Technologies, Inc.");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","node id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","element id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
    }
    offset = 5*80*sizeof(char);

    // go through all flagged bf zones...

    int *ibuf = new int[max_my_count[1]];
    int *ibuf_2 = new int[max_my_count[2]];
    int* no_flag2 = new int[nno];
    for (int ifz = 0; ifz < nfz; ++ifz) {
      const int iz = znofz[ifz];

      // ==============================
      // part header...
      // ==============================

      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s","part");
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        int part = 1+ifz;
        MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
        sprintf(cbuf,"%-80s",bfZoneVec[iz].getName().c_str());
        MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      }
      offset += int_size+2*80*sizeof(char);

      // ==============================
      // part coordinates...
      // ==============================

      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s","coordinates"); // convex polyhedron format...
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        MPI_File_write_at(fh,offset+80*sizeof(char),&count[ifz][0],1,MPI_INT,MPI_STATUS_IGNORE);
      }
      offset += 80*sizeof(char)+int_size*1;

      // write out x[3]...
      FOR_I3 {
        int ind = 0;
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);
          for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
            const int ino = noobf_v[nob];
            assert((ino >= 0)&&(ino < nno));
            // first time this node was visited by this zone
            if (!visited_no[ino]) {
              visited_no[ino] = true;
              if (no_flag[ino] == mpi_rank) {
                fbuf[ind++] = (float)x_no[ino][i]; // x[N],y[N],z[N]
              }
            }
          }
        }
        // reset visited_no...
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);
          for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
            const int ino = noobf_v[nob];
            assert((ino >= 0)&&(ino < nno));
            visited_no[ino] = false;
          }
        }
        assert(ind == my_count[ifz][0]);

        writeChunkedData<float>(fh,offset+(i*count[ifz][0]+my_disp[ifz][0])*sizeof(float),fbuf,my_count[ifz][0],mpi_comm);
      }
      offset += sizeof(float)*3*count[ifz][0];

      // ==============================
      // connectivity...
      // ==============================

      sprintf(cbuf,"%-80s","nsided"); // polygon...
      if (mpi_rank == 0) {
        // write out element block header...
        MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        MPI_File_write_at(fh,offset+80*sizeof(char),&count[ifz][1],1,MPI_INT,MPI_STATUS_IGNORE);
      }
      offset += 80*sizeof(char)+int_size*1;

      // put global index in no_flag2...
      int my_count0 = 0;
      FOR_INO no_flag2[ino] = -2;
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int ino = noobf_v[nob];
          assert((ino >= 0)&&(ino < nno));
          // first time this node was visited by this zone
          if (!visited_no[ino]) {
            visited_no[ino] = true;
            if (no_flag[ino] == mpi_rank) {
              no_flag2[ino] = my_disp[ifz][0]+my_count0++; // this node is dumped by us
            }
            else if (no_flag[ino] < mpi_size) {
              no_flag2[ino] = -1; // this node gets dumped, but not by us
            }
            else {
              assert(no_flag[ino] == mpi_size);
              no_flag2[ino] = -2; // nobody owns this node
            }
          }
        }
      }
      // reset visited_no...
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int ino = noobf_v[nob];
          assert((ino >= 0)&&(ino < nno));
          visited_no[ino] = false;
        }
      }
      assert(my_count0 == my_count[ifz][0]);
      updateNoData(no_flag2,MAX_NO_PERIODIC_DATA);

      int ind = 0;
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);

        // counts...
        ibuf[ibf-bfZoneVec[iz].ibf_f] = noobf_i[ibf+1]-noobf_i[ibf];

        // connectivity...
        for (int nof = noobf_i[ibf]; nof != noobf_i[ibf+1]; ++nof) {
          ibuf_2[ind++] = 1+no_flag2[noobf_v[nof]]; // 1-indexed
        }
      }
      assert(ind == my_count[ifz][2]);

      writeChunkedData<int>(fh,offset+my_disp[ifz][1]*int_size,ibuf,my_count[ifz][1],mpi_comm);
      offset += int_size*count[ifz][1];
      writeChunkedData<int>(fh,offset+my_disp[ifz][2]*int_size,ibuf_2,my_count[ifz][2],mpi_comm);
      offset += int_size*count[ifz][2];
    }
    delete[] ibuf;
    delete[] ibuf_2;
    delete[] no_flag2;

    MPI_File_set_size(fh,offset);
    MPI_File_close(&fh);

  }

  // now lets write out the scalar and vector files for this step...
  COUT1("StaticSolver::writeFlaggedBfZonesEnsightVarFiles()");

  if (dw->nscalars_cv > 0||dw->nscalars_bf > 0) {
    double *no_data = new double[nno];
    for (int isc = 0; isc < dw->nscalars_cv; ++isc) {

      // open file...
      sprintf(cbuf,"%s.%08d.sca%02d",dw->name.c_str(),step,isc); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);
      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->scalar_descs_cv[isc].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      }
      offset += 80*sizeof(char);

      // get data on nodes...
      setNoDN(no_data,dw->scalar_names_cv[isc]);
      for (int ifz = 0; ifz < nfz; ++ifz) {
        const int iz = znofz[ifz];

        if (mpi_rank == 0) {
          sprintf(cbuf,"%-80s","part");
          MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          int part = 1+ifz;
          MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
          sprintf(cbuf,"%-80s","coordinates");
          MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        }
        offset += int_size+2*80*sizeof(char);

        // write data...

        int ind = 0;
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);
          for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
            const int ino = noobf_v[nob];
            assert((ino >= 0)&&(ino < nno));
            // first time this node was visited by this zone
            if (!visited_no[ino]) {
              visited_no[ino] = true;
              if (no_flag[ino] == mpi_rank) {
                fbuf[ind++] = (float)no_data[ino]; // x[N],y[N],z[N]
              }
            }
          }
        }
        // reset visited_no...
        for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
          assert(zone_bf[ibf] == iz);
          for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
            const int ino = noobf_v[nob];
            assert((ino >= 0)&&(ino < nno));
            visited_no[ino] = false;
          }
        }
        assert(ind == my_count[ifz][0]);

        writeChunkedData<float>(fh,offset+my_disp[ifz][0]*sizeof(float),fbuf,my_count[ifz][0],mpi_comm);

        offset += sizeof(float)*count[ifz][0];
      }

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);
    }
    if (dw->nscalars_bf > 0) {
      vector<string> var_subvec_bf(dw->nzones_bf);
      for (int isc = 0; isc < dw->nscalars_bf/dw->nzones_bf; ++isc) {

        // open file...
        sprintf(cbuf,"%s.%08d.sca%02d",dw->name.c_str(),step,isc+dw->nscalars_cv); // see case file
        MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
        MPI_File_delete(cbuf,MPI_INFO_NULL);
        MPI_File fh;
        MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

        // write header...
        MPI_Offset offset = 0;
        if (mpi_rank == 0) {
          sprintf(cbuf,"%-80s",dw->scalar_descs_bf[isc].c_str()); // description, hashes??
          MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        }
        offset += 80*sizeof(char);

        for (int izone = 0; izone < dw->nzones_bf; ++izone) {
          var_subvec_bf[izone] = dw->scalar_names_bf[izone*(dw->nscalars_bf/dw->nzones_bf)+isc];
        }
        setBfDN(no_data,var_subvec_bf);

        for (int ifz = 0; ifz < nfz; ++ifz) {
          const int iz = znofz[ifz];

          if (mpi_rank == 0) {
            sprintf(cbuf,"%-80s","part");
            MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
            int part = 1+ifz;
            MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
            sprintf(cbuf,"%-80s","coordinates");
            MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          }
          offset += int_size+2*80*sizeof(char);

          // write data...

          int ind = 0;
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
              const int ino = noobf_v[nob];
              assert((ino >= 0)&&(ino < nno));
              // first time this node was visited by this zone
              if (!visited_no[ino]) {
                visited_no[ino] = true;
                if (no_flag[ino] == mpi_rank) {
                  fbuf[ind++] = (float)no_data[ino]; // x[N],y[N],z[N]
                }
              }
            }
          }
          // reset visited_no...
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
              const int ino = noobf_v[nob];
              assert((ino >= 0)&&(ino < nno));
              visited_no[ino] = false;
            }
          }
          assert(ind == my_count[ifz][0]);

          writeChunkedData<float>(fh,offset+my_disp[ifz][0]*sizeof(float),fbuf,my_count[ifz][0],mpi_comm);
          offset += sizeof(float)*count[ifz][0];
        }

        // close file...
        MPI_File_set_size(fh,offset);
        MPI_File_close(&fh);
      }
    }
    delete[] no_data;
  }
  if (dw->nvectors_cv > 0 || dw->nvectors_bf > 0) {
    double (*no_data)[3] = new double[nno][3];
    for (int ivec = 0; ivec < dw->nvectors_cv; ++ivec) {

      // open file...
      sprintf(cbuf,"%s.%08d.vec%02d",dw->name.c_str(),step,ivec); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);
      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->vector_descs_cv[ivec].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      }
      offset += 80*sizeof(char);

      setNoDN3(no_data,dw->vector_names_cv[ivec]);
      for (int ifz = 0; ifz < nfz; ++ifz) {
        const int iz = znofz[ifz];

        if (mpi_rank == 0) {
          sprintf(cbuf,"%-80s","part");
          MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          int part = 1+ifz;
          MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
          sprintf(cbuf,"%-80s","coordinates");
          MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        }
        offset += int_size+2*80*sizeof(char);

        // write data...
        FOR_I3 {
          int ind = 0;
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
              const int ino = noobf_v[nob];
              assert((ino >= 0)&&(ino < nno));
              // first time this node was visited by this zone
              if (!visited_no[ino]) {
                visited_no[ino] = true;
                if (no_flag[ino] == mpi_rank) {
                  fbuf[ind++] = (float)no_data[ino][i];
                }
              }
            }
          }
          // reset visited_no...
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
              const int ino = noobf_v[nob];
              assert((ino >= 0)&&(ino < nno));
              visited_no[ino] = false;
            }
          }
          assert(ind == my_count[ifz][0]);

          writeChunkedData<float>(fh,offset+(i*count[ifz][0]+my_disp[ifz][0])*sizeof(float),fbuf,my_count[ifz][0],mpi_comm);
        }
        offset += sizeof(float)*3*count[ifz][0];
      }

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);
    }
    if (dw->nvectors_bf > 0) {
      vector<string> var_subvec_bf(dw->nzones_bf);
      for (int ivec = 0; ivec < dw->nvectors_bf/dw->nzones_bf; ++ivec) {

        // open file...
        sprintf(cbuf,"%s.%08d.vec%02d",dw->name.c_str(),step,ivec+dw->nvectors_cv); // see case file
        MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
        MPI_File_delete(cbuf,MPI_INFO_NULL);
        MPI_File fh;
        MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

        // write header...
        MPI_Offset offset = 0;
        if (mpi_rank == 0) {
          sprintf(cbuf,"%-80s",dw->vector_descs_bf[ivec].c_str()); // description, hashes??
          MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
        }
        offset += 80*sizeof(char);

        for (int izone = 0; izone < dw->nzones_bf; ++izone) {
          var_subvec_bf[izone] = dw->vector_names_bf[izone*(dw->nvectors_bf/dw->nzones_bf)+ivec];
        }
        setBfDN3(no_data,var_subvec_bf);
        for (int ifz = 0; ifz < nfz; ++ifz) {
          const int iz = znofz[ifz];

          if (mpi_rank == 0) {
            sprintf(cbuf,"%-80s","part");
            MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
            int part = 1+ifz;
            MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
            sprintf(cbuf,"%-80s","coordinates");
            MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
          }
          offset += int_size+2*80*sizeof(char);

          // write data...
          FOR_I3 {
            int ind = 0;
            for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
              assert(zone_bf[ibf] == iz);
              for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
                const int ino = noobf_v[nob];
                assert((ino >= 0)&&(ino < nno));
                // first time this node was visited by this zone
                if (!visited_no[ino]) {
                  visited_no[ino] = true;
                  if (no_flag[ino] == mpi_rank) {
                    fbuf[ind++] = (float)no_data[ino][i];
                  }
                }
              }
            }
            // reset visited_no...
            for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
              assert(zone_bf[ibf] == iz);
              for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
                const int ino = noobf_v[nob];
                assert((ino >= 0)&&(ino < nno));
                visited_no[ino] = false;
              }
            }
            assert(ind == my_count[ifz][0]);

            writeChunkedData<float>(fh,offset+(i*count[ifz][0]+my_disp[ifz][0])*sizeof(float),fbuf,my_count[ifz][0],mpi_comm);
          }
          offset += sizeof(float)*3*count[ifz][0];
        }

        // close file...
        MPI_File_set_size(fh,offset);
        MPI_File_close(&fh);
      }
    }
    delete[] no_data;
  }

  // cleanup
  delete[] my_count;
  delete[] my_disp;
  delete[] count;
  delete[] visited_no;
  delete[] znofz;
  delete[] no_flag;
  delete[] cbuf;
  delete[] fbuf;
}

void StaticSolver::writeFlaggedLpsEnsight(DataWriter* dw,const int step,const double time) {

  // update case file with info about this step...
  updateEnsightCaseFile(dw,step,time);

  const int nlp = *lpHelperVec[dw->lp_index].size_ptr;
  int8 my_count = 0;
  for (int ip = 0; ip < nlp; ++ip) {
    if (dw->lp_flag[ip] != 0)
      ++my_count;
  }
  int8 count;
  MPI_Allreduce(&my_count,&count,1,MPI_INT8,MPI_SUM,mpi_comm);

  int8 my_disp;
  MPI_Scan(&my_count,&my_disp,1,MPI_INT8,MPI_SUM,mpi_comm);
  my_disp -= my_count;

  if (mpi_rank == 0)
    cout << " > global flagged particles: " << count << endl;

  FOR_RANK {
    if (mpi_rank == rank)
      cout << my_disp << " " << my_count << endl;
    MPI_Barrier(mpi_comm);
  }

  int np = (int)my_count;

  // geom and scalar files need these array, so just allocate here...
  char * cbuf = new char[256];
  float * fbuf = new float[np];

  {

    // lets build the geo file used for all transient data associated to this ensight case...
    COUT1("StaticSolver::writeFlaggedLpsEnsightGeoFile()");

    // open the file...

    sprintf(cbuf,"%s.%08d.geo",dw->name.c_str(),step);
    MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
    MPI_File_delete(cbuf,MPI_INFO_NULL);
    MPI_File fh;
    MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    // ==============================
    // header...
    // ==============================

    MPI_Offset offset = 0;
    if (mpi_rank == 0) {
      // write out header

      sprintf(cbuf,"%-80s","C Binary");
      MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Geometry file"); // should probably include hashes here
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","Produced by Cascade Technologies, Inc.");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","node id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      sprintf(cbuf,"%-80s","element id assign");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
    }
    offset = 5*80*sizeof(char);

    // all ranks write there flagged data...

    offset += mpi_rank*(int_size+2*80*sizeof(char)); // part headers
    offset += mpi_rank*(80*sizeof(char)+int_size*1); // coordinates, nlp
    offset += sizeof(float)*3*my_disp; // xp,yp,zp
    offset += mpi_rank*(80*sizeof(char)+int_size*1); // element type, nlp
    offset += int_size*my_disp;

    // ==============================
    // part header...
    // ==============================

    sprintf(cbuf,"%-80s","part");
    MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
    int part = 1+mpi_rank;
    MPI_File_write_at(fh,offset+80*sizeof(char),&part,1,MPI_INT,MPI_STATUS_IGNORE);
    sprintf(cbuf,"%-80d",mpi_rank);
    MPI_File_write_at(fh,offset+int_size+80*sizeof(char),cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
    offset += int_size+2*80*sizeof(char);

    // ==============================
    // coordinates...
    // ==============================

    sprintf(cbuf,"%-80s","coordinates");
    MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
    MPI_File_write_at(fh,offset+80*sizeof(char),&np,1,MPI_INT,MPI_STATUS_IGNORE);
    offset += 80*sizeof(char)+int_size*1;

    {
      CtiRegister::CtiData * data = CtiRegister::getRegisteredCtiData(lpHelperVec[dw->lp_index].name+":xp",false);
      FOR_I3 {
        int ind = 0;
        for (int ip = 0; ip < nlp; ++ip) {
          if (dw->lp_flag[ip] != 0)
            fbuf[ind++] = (float)data->dn3(ip,i); // x[N],y[N],z[N]
        }
        assert(ind == np);
        writeChunkedData<float>(fh,offset,fbuf,np,mpi_comm);
        offset += sizeof(float)*np;
      }
    }

    // ==============================
    // connectivity...
    // ==============================

    sprintf(cbuf,"%-80s","point"); // point format...
    // write out element block header...
    MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
    MPI_File_write_at(fh,offset+80*sizeof(char),&np,1,MPI_INT,MPI_STATUS_IGNORE);
    offset += 80*sizeof(char)+int_size*1;

    int* ibuf = new int[np];
    {
      int ind = 0;
      for (int ip = 0; ip < nlp; ++ip) {
        if (dw->lp_flag[ip] != 0) {
          ibuf[ind] = ind+1;
          ++ind;
        }
      }
      assert(ind == np);
    }
    writeChunkedData<int>(fh,offset,ibuf,np,mpi_comm);
    delete[] ibuf;
    offset += int_size*my_count;

    // give everyone the same offset in end to set size

    offset = 5*80*sizeof(char);
    offset += mpi_size*(int_size+2*80*sizeof(char)); // part headers
    offset += mpi_size*(80*sizeof(char)+int_size); // coordinates, np
    offset += sizeof(float)*3*count; // xp,yp,zp
    offset += mpi_size*(80*sizeof(char)+int_size); // element type, np
    offset += int_size*count;

    MPI_File_set_size(fh,offset);
    MPI_File_close(&fh);

  }

  // now lets write out the scalar and vector files for this step...
  COUT1("StaticSolver::writeFlaggedLpsEnsightVarFiles()");

  if (dw->nscalars_lp > 0) {
    for (int isc = 0; isc < dw->nscalars_lp; ++isc) {

      // open file...
      sprintf(cbuf,"%s.%08d.sca%02d",dw->name.c_str(),step,isc); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);
      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->scalar_descs_lp[isc].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      }
      offset = 80*sizeof(char);
      offset += mpi_rank*(int_size+2*80*sizeof(char)); // part header
      offset += sizeof(float)*my_disp; // part scalar nodal record

      sprintf(cbuf,"%-80s","part");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      int part = 1+mpi_rank;
      MPI_File_write_at(fh,offset,&part,1,MPI_INT,MPI_STATUS_IGNORE);
      offset += int_size;
      sprintf(cbuf,"%-80s","coordinates");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);

      // write data...
      CtiRegister::CtiData * data = CtiRegister::getCtiData(dw->scalar_names_lp[isc]);
      {
        int ind = 0;
        for (int ip = 0; ip < nlp; ++ip) {
          if (dw->lp_flag[ip] != 0)
            fbuf[ind++] = (float)data->dn(ip);
        }
        assert(ind == np);
      }
      writeChunkedData<float>(fh,offset,fbuf,np,mpi_comm);
      offset += sizeof(float)*np;

      // give everyone the same offset to size size

      offset = 80*sizeof(char);
      offset += mpi_size*(int_size+2*80*sizeof(char)); // part header
      offset += sizeof(float)*count; // part scalar nodal record

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);

    }
  }
  if (dw->nvectors_lp > 0) {
    for (int ivec = 0; ivec < dw->nvectors_lp; ++ivec) {

      // open file...
      sprintf(cbuf,"%s.%08d.vec%02d",dw->name.c_str(),step,ivec); // see case file
      MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
      MPI_File_delete(cbuf,MPI_INFO_NULL);
      MPI_File fh;
      MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // write header...
      MPI_Offset offset = 0;
      if (mpi_rank == 0) {
        sprintf(cbuf,"%-80s",dw->vector_descs_lp[ivec].c_str()); // description, hashes??
        MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      }
      offset += 80*sizeof(char);

      offset += mpi_rank*(int_size+2*80*sizeof(char)); // part header
      offset += 3*sizeof(float)*my_disp; // part vector nodal record

      sprintf(cbuf,"%-80s","part");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);
      int part = 1+mpi_rank;
      MPI_File_write_at(fh,offset,&part,1,MPI_INT,MPI_STATUS_IGNORE);
      offset += int_size;
      sprintf(cbuf,"%-80s","coordinates");
      MPI_File_write_at(fh,offset,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
      offset += 80*sizeof(char);

      // write data...
      CtiRegister::CtiData * data = CtiRegister::getCtiData(dw->vector_names_lp[ivec]);
      FOR_I3 {
        int ind = 0;
        for (int ip = 0; ip < nlp; ++ip) {
          if (dw->lp_flag[ip] != 0)
            fbuf[ind++] = (float)data->dn3(ip,i);
        }
        assert(ind == np);
        writeChunkedData<float>(fh,offset,fbuf,np,mpi_comm);
        offset += sizeof(float)*np;
      }

      // give everyone the same offset to size size

      offset = 80*sizeof(char);
      offset += mpi_size*(int_size+2*80*sizeof(char)); // part header
      offset += 3*sizeof(float)*count; // part scalar nodal record

      // close file...
      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);
    }
  }

  // cleanup...
  delete[] cbuf;
  delete[] fbuf;
}

void StaticSolver::writeSimpleTriVecStl(DataWriter* dw,const vector<SimpleTriWithWgts>& triVec,const int step,const double time) {

  if (mpi_rank == 0)
    cout << "StaticSolver::writeSimpleTriVecStl()" << endl;

  // open file...
  char cbuf[256];
  sprintf(cbuf,"%s.%08d.stl",dw->name.c_str(),step);
  MiscUtils::mkdir_for_file_collective(cbuf); // allow the user to have specified a subdirectory...
  MPI_File_delete(cbuf,MPI_INFO_NULL);
  MPI_File fh;
  MPI_File_open(mpi_comm,cbuf,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

  int my_count = triVec.size();
  int count;
  MPI_Allreduce(&my_count,&count,1,MPI_INT,MPI_SUM,mpi_comm);
  int my_disp;
  MPI_Scan(&my_count,&my_disp,1,MPI_INT,MPI_SUM,mpi_comm);
  my_disp -= my_count;

  // ==============================
  // header...
  // ==============================

  MPI_Offset offset = 0;
  if (mpi_rank == 0) {
    sprintf(cbuf,"%-80s","Binary STL file");
    MPI_File_write(fh,cbuf,80,MPI_CHAR,MPI_STATUS_IGNORE);
    offset += 80*sizeof(char);
    MPI_File_write_at(fh,offset,&count,1,MPI_INT,MPI_STATUS_IGNORE);
  }
  offset = int_size+80*sizeof(char);

  // ==============================
  // tri data (normal[3 floats],x0[3 floats],x1[3 floats],x2[3 floats],attribute[2 bytes]...
  // ==============================

  assert(sizeof(unsigned char) == 1);
  assert(sizeof(float) == 4); // could generalize
  const int nchar = 50*my_count; // 12 floats and 1 short = 50 bytes
  char* data = new char[nchar];
  int ichar = 0;
  for (int ist = 0; ist < my_count; ++ist) {
    FOR_I3 FOR_J4 data[ichar++] = '\0';
    FOR_I3 {
      const float my_float = (float)triVec[ist].x0[i];
      FOR_J4 data[ichar++] = ((unsigned char*)&my_float)[j];
    }
    FOR_I3 {
      const float my_float = (float)triVec[ist].x1[i];
      FOR_J4 data[ichar++] = ((unsigned char*)&my_float)[j];
    }
    FOR_I3 {
      const float my_float = (float)triVec[ist].x2[i];
      FOR_J4 data[ichar++] = ((unsigned char*)&my_float)[j];
    }
    FOR_I2 data[ichar++] = '\0';
  }
  assert(ichar == 50*my_count);
  writeChunkedData<char>(fh,offset+50*my_disp*sizeof(char),data,50*my_count,mpi_comm);
  delete[] data;
  offset += sizeof(char)*50*count;

  MPI_File_set_size(fh,offset);
  MPI_File_close(&fh);

}

#undef IHN_VAL
