#ifndef _CHT_HPP_
#define _CHT_HPP_
#include "Prcomm.hpp"
#include "SimpleFunc.hpp"

class ChtBoundaryData {
public:
  int np; // the send-side points associated with the nodes touching a particular CHT bf zone
  int zone; // fluid side zone
  double * Tp; // pointer into a buffer for temperature stored by the bc...
  double * qp; // pointer into a buffer for heat flux x area stored by the bc...
  int * send_rank; // send side is the fluid side
  int np_recv;
  int * recv_rank; // recv side is the fem side // TODO: call this recv_index?
  int * recv_ino;  // "
  ChtBoundaryData(double * Tp,double * qp,const int np,const int zone) {
    this->Tp = Tp;
    this->qp = qp;
    this->np = np;
    this->zone = zone;
    send_rank = NULL;
    recv_rank = NULL;
    recv_ino = NULL;
  }
  void clear() {
    DELETE(send_rank);
    DELETE(recv_rank);
    DELETE(recv_ino);
  }
};

enum ChtMode {
  INIT,
  TSTART,
  TFINISH,
  QSTART,
  QFINISH
};

class Cht {
private:

  ChtMode mode;
  int *send_count,*send_disp;
  int *recv_count,*recv_disp;
  Adt<double> * adt;
  Adt<double> * bbox_adt;
  double tol;
  double (*grad_te)[4][3];
  double * vol_te;
  //double * vol_no;
  double * laplace_k_noono;
  uint8 * rbi_g; // nodal ghost rbi_g
  vector<CvPrcomm> noPrcommVec; // TODO: call CvPrcomm->Prcomm, and Prcomm->PrcommFull
  // mpi requests map the void data pointer to a set of
  // non-blocking requests that have yet to be completed...
  // this map relates update*DataStart() and update*DataFinish() non-blocking
  // ghost update routines...
  map<const void*,MpiRequestStuff*> mpiRequestMap;
  vector<ChtBoundaryData> boundaryDataVec;
  // communication between the no-based data (FEM/CHT side) and
  // the sp-based data (Flow solver bc side) has its own non-blocking
  // communication...
  class NoSpComm {
  public:
    int rank;
    int count,disp;
  };
  vector<NoSpComm> noCommVec;
  vector<NoSpComm> spCommVec;
  int no_buf_size;
  int sp_buf_size;
  double * no_buf_double;
  double * sp_buf_double;
  MPI_Request * spRequestArray;
  MPI_Request * noRequestArray;

public:

  // for striping io
  int8 nno_global;
  int8* ino_global;
  int8* noora_striped;
  DistributedDataExchanger* dde_striped;

  int * noono_i;
  int * noono_v;

  int nno; // active nodes
  int nno_g; // active+ghost nodes (nno_g >= nno)
  int nno_i; // active internal nodes: edge nbrs are also active (nno_i <= nno). Used for latency hiding
  int nte; // all tets
  int nte_i; // tets touching local nodes only (nte_i <= nte_ib)
  int nte_ib; // tets touching local and ghost nodes that we own (nte_ib <= nte)
  //int nst; // all surface tri's
  //int nst_i; // surface tri's touching local nodes only (nst_i <= nst_ib)
  //int nst_ib; // surface tri's touching local and ghost nodes that we own (nst_ib <= nst)
  int (*noote)[4]; // node-of-tet
  //int (*noost)[3]; // node-of-surface-tri
  //int *znost;
  int *znote; // volume zone of tet
  double (*x_no)[3];
  // property simple functions: index matches znote...
  vector<SimpleFunc*> rhoFuncVec;
  vector<SimpleFunc*> cpFuncVec;
  vector<SimpleFunc*> kFuncVec;
  //SimpleFunc *rho_func,*cp_func,*k_func;
  //double * rho_no;
  //double * cp_no;
  //double * k_no;
  double *rho_cp_vol_no; // product of rho*cp*vol at the nodes
  double alpha;
  double * T_no;
  double * q_no;
  double * rhs_no[3]; // rhs's for RK3
  double *A; // operator

  int write_cht_interval;
  string write_cht_prefix;

  Cht() {

    mode = INIT;
    send_count = NULL;
    send_disp = NULL;
    recv_count = NULL;
    recv_disp = NULL;
    adt = NULL;
    bbox_adt = NULL;
    tol = 0.0;
    noono_i = NULL;
    noono_v = NULL;
    grad_te = NULL;
    vol_te = NULL;
    //vol_no = NULL;
    laplace_k_noono = NULL;
    rbi_g = NULL;

    // no<->sp communication...
    no_buf_size = sp_buf_size = 0;
    no_buf_double = NULL;
    sp_buf_double = NULL;
    spRequestArray = NULL;
    noRequestArray = NULL;

    nno_global = 0;
    ino_global = NULL;
    noora_striped = NULL;
    dde_striped = NULL;

    nno = nno_g = nno_i = nte = nte_i = nte_ib = 0;
    //nst = nst_i = nst_ib = 0;
    noote = NULL;
    //noost = NULL;
    //znost = NULL;
    znote = NULL;
    x_no = NULL;
    T_no = NULL;
    q_no = NULL;
    rhs_no[0] = NULL;
    A = NULL;

    // set the properties...
    // cannot do this until mesh is read in because there
    // may be multiple volume zones with different props...
    /*
    Param * rho_param = getParam("CHT.RHO");
    Param * cp_param = getParam("CHT.CP");
    Param * k_param = getParam("CHT.K");
    if ((rho_param == NULL)||(cp_param == NULL)||(k_param == NULL)) {
      CERR("CHT missing one or more temperature-dependent properties: CHT.RHO, CHT.CP, CHT.K");
    }
    rho_func = processSimpleFunc(rho_param);
    cp_func = processSimpleFunc(cp_param);
    k_func = processSimpleFunc(k_param);
    */

    // nodal-based versions of rho,cp,k...
    //rho_no = NULL;
    //cp_no = NULL;
    //k_no = NULL;
    // now just the product of all 3...
    rho_cp_vol_no = NULL;

    // alpha is used to adjust the time-scale ratio. This number
    // multiples the solid cp to (significantly) reduce the thermal inertia,
    // while keeping the Biot number constant...
    alpha = getDoubleParam("CHT.ALPHA");
    write_cht_interval = -1; // -1 means not parsed, 0 means parsed and not present, 1,2,3... means parsed and prefix set
    write_cht_prefix = "result"; // default prefix

  }

  ~Cht() {

    DELETE(send_count);
    DELETE(send_disp);
    DELETE(recv_count);
    DELETE(recv_disp);
    if (adt != NULL) {
      delete adt;
      adt = NULL;
    }
    if (bbox_adt != NULL) {
      delete bbox_adt;
      bbox_adt = NULL;
    }

    DELETE(noono_i);
    DELETE(noono_v);
    DELETE(grad_te);
    DELETE(vol_te);
    //DELETE(vol_no);
    DELETE(laplace_k_noono);
    DELETE(rbi_g);

    DELETE(no_buf_double);
    DELETE(sp_buf_double);
    DELETE(spRequestArray);
    DELETE(noRequestArray);

    DELETE(ino_global);
    DELETE(noora_striped);
    if (dde_striped) {
      delete dde_striped;
      dde_striped = NULL;
    }

    DELETE(noote);
    //DELETE(noost);
    //DELETE(znost);
    DELETE(znote);
    DELETE(x_no);
    DELETE(T_no);
    DELETE(q_no);
    DELETE(rhs_no[0]); // this is right: see trick during new below
    DELETE(A);

    for (int ii = 0, size = rhoFuncVec.size(); ii < size; ++ii) {
      delete rhoFuncVec[ii];
    }
    for (int ii = 0, size = cpFuncVec.size(); ii < size; ++ii) {
      delete cpFuncVec[ii];
    }
    for (int ii = 0, size = kFuncVec.size(); ii < size; ++ii) {
      delete kFuncVec[ii];
    }

    delete[] rho_cp_vol_no;
    //DELETE(rho_no);
    //DELETE(cp_no);
    //DELETE(k_no);

    // boundaryDataVec does not clean up after itself...
    for (int ii = 0; ii < boundaryDataVec.size(); ++ii)
      boundaryDataVec[ii].clear();
  }

  // ===========================================================
  // ghost node data update routines...
  // ===========================================================

#define UPDATE_START updateNoDataStart
#define UPDATE_FINISH updateNoDataFinish
#define UPDATE updateNoData
#define PRCOMM_VEC noPrcommVec

  // updateNoData(int * s...

#define T int
#define PACK_BUF pack_buf_int
#define MPI_T MPI_INT
#define UPDATE_TAG 12121
#include "updateDN.hpp"
#undef T
#undef PACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateCvData(double * s and double (*s)[3]...

#define T double
#define PACK_BUF pack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 12122
#include "updateDN.hpp"
#include "updateDN3.hpp"
#undef T
#undef PACK_BUF
#undef MPI_T
#undef UPDATE_TAG

#undef UPDATE_START
#undef UPDATE_FINISH
#undef UPDATE
#undef PRCOMM_VEC

  // ===========================================================
  // end of ghost node data update routines
  // ===========================================================

  void report(const int step,const double time) {

    // NOTE: this report is only approximate because the
    // temperature-dependence of energy is not considered

    if (rho_cp_vol_no != NULL) {

      double my_buf[2] = { 0.0, 0.0 };
      double my_minmax[2] = { HUGE_VAL, HUGE_VAL };

      FOR_INO {
	my_buf[0] += rho_cp_vol_no[ino]*alpha;
	my_buf[1] += rho_cp_vol_no[ino]*alpha*T_no[ino];
	my_minmax[0] = min(my_minmax[0],T_no[ino]);
	my_minmax[1] = min(my_minmax[1],-T_no[ino]);
      }

      double buf[2];
      MPI_Reduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      double minmax[2];
      MPI_Reduce(my_minmax,minmax,2,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
      if (mpi_rank == 0) {
	cout << " > CHT summary: step, time, energy(approx), T(min,avg,max): " <<
	  step << " " << time << " " <<
	  buf[1] << " " << minmax[0] << " " << buf[1]/buf[0] << " " << -minmax[1] << endl;
      }

    }
    else {

      if (mpi_rank == 0) {
	cout << " > CHT: report not available" << endl;
      }

    }

  }

  void registerBoundaryData(const double (* const xp)[3],double* Tp,double* qp,int* ssp_flag,int& np,const int nssp,const int zone) {

    assert(mode == INIT);

    // all Cht bcs call this routine collectively (i.e. even if they have no data), so,
    // on the first time, use this as a chance to prepare the adt...

    if (adt == NULL) {

      if (mpi_rank == 0) cout << " > building x_no adt for CHT and computing node matching tolerance..." << endl;

      // put the local nodes in an adt...

      assert(bbox_adt == NULL);
      adt = new Adt<double>(nno,x_no,x_no);

      // the top leaf of this local adt stores the local bbox...
      double my_bbmin[3],my_bbmax[3];
      adt->getBbox(my_bbmin,my_bbmax);
      double (*bbmin)[3] = new double[mpi_size][3];
      double (*bbmax)[3] = new double[mpi_size][3];
      MPI_Allgather(my_bbmin,3,MPI_DOUBLE,bbmin,3,MPI_DOUBLE,mpi_comm);
      MPI_Allgather(my_bbmax,3,MPI_DOUBLE,bbmax,3,MPI_DOUBLE,mpi_comm);
      bbox_adt = new Adt<double>(mpi_size,bbmin,bbmax);
      delete[] bbmin;
      delete[] bbmax;

      // also use this as an opportunity to calculate eps. Recall
      // nte_ib are the unique tets that we own...
      double my_d2_min = HUGE_VAL;
      for (int ite = 0; ite < nte_ib; ++ite) {
        for (int i = 0; i < 3; ++i) {
          for (int j = i+1; j < 4; ++j) {
            assert(noote[ite][i] != noote[ite][j]);
            const double d2 = DIST2(x_no[noote[ite][i]],x_no[noote[ite][j]]);
            assert(d2 > 0.0);
            my_d2_min = min(my_d2_min,d2);
          }
        }
      }
      MPI_Allreduce(&my_d2_min,&tol,1,MPI_DOUBLE,MPI_MIN,mpi_comm);
      // take 1% of this minimum spacing...
      const double eps = getDoubleParam("CHT_EPS",0.01);
      tol = eps*sqrt(tol);
      if (mpi_rank == 0) cout << " > node matching tol for CHT: " << tol << endl;

    }

    double my_d2_match = 0.0;

    // add boundary data to the registered data...

    boundaryDataVec.push_back(ChtBoundaryData(Tp,qp,np,zone));
    ChtBoundaryData * cbd = &boundaryDataVec.back();

    // build and store the communication pattern for these points...
    // note that ultimately we want one message only between each
    // pair that may involve multiple zones, so this communication
    // pattern gets translated into paired communicators only after
    // all CHT boundary conditions have been registered, and the
    // interleaving of messages can be properly handled...

    if (send_count == NULL) send_count = new int[mpi_size];
    if (send_disp == NULL) send_disp = new int[mpi_size];

    FOR_RANK send_count[rank] = 0;

    vector<pair<int,int> > rankDataPairVec;
    vector<int> rankVec;
    double bbmin[3],bbmax[3];
    for (int ip = 0; ip < np; ++ip) {
      assert(rankVec.empty());
      FOR_I3 bbmin[i] = xp[ip][i]-tol;
      FOR_I3 bbmax[i] = xp[ip][i]+tol;
      bbox_adt->buildListForBBox(rankVec,bbmin,bbmax);
      // allow for fluid points to be untouched...
      //assert(!rankVec.empty()); // should always find one, unless there is a problem matching CHT
      for (int ii = 0; ii < rankVec.size(); ++ii) {
        rankDataPairVec.push_back(pair<int,int>(rankVec[ii],ip));
        send_count[rankVec[ii]] += 3;
      }
      rankVec.clear();
    }

    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
    assert(send_count_sum == rankDataPairVec.size()*3);

    // pack the xp's to get the correspondence...
    double * send_buf_double = new double[send_count_sum];
    for (int isend = 0; isend < rankDataPairVec.size(); ++isend) {
      const int rank = rankDataPairVec[isend].first;
      const int ip = rankDataPairVec[isend].second;
      FOR_I3 send_buf_double[send_disp[rank]+i] = xp[ip][i];
      send_disp[rank] += 3;
    }

    // rewind...
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

    if (recv_count == NULL) recv_count = new int[mpi_size];
    if (recv_disp == NULL) recv_disp = new int[mpi_size];

    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    assert(recv_count_sum%3 == 0);

    double * recv_buf_double = new double[recv_count_sum];
    MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
                  recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
    delete[] send_buf_double;

    int * recv_buf_int = new int[recv_count_sum/3];
    cbd->np_recv = 0;
    for (int irecv = 0; irecv < recv_count_sum/3; ++irecv) {
      // look for a match in x_no for this point. This returns
      // -1 if no match found...
      recv_buf_int[irecv] = getMatchingPoint(recv_buf_double+irecv*3);
      if (recv_buf_int[irecv] >= 0) {
        const double d2 = DIST2(recv_buf_double+irecv*3,x_no[recv_buf_int[irecv]]);
        my_d2_match = max(my_d2_match,d2);
        //cout << "Got match on rank " << mpi_rank << " " << sqrt(d2) << endl;
        ++cbd->np_recv;
      }
    }
    delete[] recv_buf_double;

    FOR_RANK {
      send_count[rank] /= 3;
      send_disp[rank] /= 3;
      recv_count[rank] /= 3;
      recv_disp[rank] /= 3;
    }

    // now loop through and record the rank and matching point on the recv side
    // based on this matching. In future exchanges, the matched points will be
    // send by the owning processor...
    assert(cbd->recv_rank == NULL);
    cbd->recv_rank = new int[cbd->np_recv];
    assert(cbd->recv_ino == NULL);
    cbd->recv_ino = new int[cbd->np_recv];
    int ip_recv = 0;
    FOR_RANK {
      for (int irecv = recv_disp[rank]; irecv < recv_disp[rank]+recv_count[rank]; ++irecv) {
        if (recv_buf_int[irecv] >= 0) {
          cbd->recv_rank[ip_recv] = rank;
          cbd->recv_ino[ip_recv] = recv_buf_int[irecv];
          ++ip_recv;
        }
      }
    }
    assert(ip_recv == cbd->np_recv);

    // send these matches back to the send side (the fluid side) to
    // confirm every point is matched once and only once...
    int * send_buf_int = new int[send_count_sum/3];
    MPI_Alltoallv(recv_buf_int,recv_count,recv_disp,MPI_INT,
                  send_buf_int,send_count,send_disp,MPI_INT,
                  mpi_comm);
    delete[] recv_buf_int;

    // now we should have one and only one match (positive int) in
    // send_buf_int for each ip in rank. For now, we just record the
    // rank for each point in send_index. The
    assert(cbd->send_rank == NULL);
    cbd->send_rank = new int[np];
    for (int ip = 0; ip < np; ++ip)
      cbd->send_rank[ip] = -1;
    for (int isend = 0; isend < rankDataPairVec.size(); ++isend) {
      const int rank = rankDataPairVec[isend].first;
      const int ip = rankDataPairVec[isend].second;
      // the matching ino on rank is then...
      const int ino = send_buf_int[send_disp[rank]++];
      if (ino >= 0) {
        // this was a match on the recv side, so store the rank...
        assert(cbd->send_rank[ip] == -1);
        cbd->send_rank[ip] = rank;
      }
    }
    delete[] send_buf_int;

    // make sure all were set: confirms each Tp will go to one and only one rank...
    int my_nmiss = 0;
    for (int ip = 0; ip < np; ++ip) {
      // allow for fluid points to be untouched...
      //assert(cbd->send_rank[ip] >= 0);
      if (!(cbd->send_rank[ip] >= 0)) {
        assert(cbd->send_rank[ip] == -1);
        ++my_nmiss;
      }
    }
    int my_buf[2] = {my_nmiss,np};
    int buf[2];
    MPI_Reduce(my_buf,buf,2,MPI_INT,MPI_SUM,0,mpi_comm);
    if (mpi_rank == 0)
      cout << " > could not match " << buf[0] << " of " << buf[1] << " points." << endl;

    // now we can reduce the match and report...
    double d2_match;
    MPI_Reduce(&my_d2_match,&d2_match,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    if (mpi_rank == 0)
      cout << " > matched points within dist: " << sqrt(d2_match) << endl;

    // reorder points, skipping unmatched...
    if (my_nmiss > 0) {
      const int np0 = np;
      np = 0;
      int* ipoip0 = new int[np0];
      for (int ip0 = 0; ip0 < np0; ++ip0) {
        if (cbd->send_rank[ip0] >= 0) {
          ipoip0[ip0] = np;
          cbd->send_rank[np] = cbd->send_rank[ip0];
          ++np;
        }
        else {
          assert(cbd->send_rank[ip0] == -1);
          ipoip0[ip0] = -1;
        }
      }
      assert(np == np0-my_nmiss);
      cbd->np = np;
      for (int issp = 0; issp < nssp; ++issp) {
        if (ssp_flag[issp] >= 0) {
          assert(ssp_flag[issp] < np0);
          ssp_flag[issp] = ipoip0[ssp_flag[issp]];
        }
        else {
          assert(ssp_flag[issp] == -1);
        }
      }
      delete[] ipoip0;
    }

  }

  /*
  void setZnost(const int nzones) {

    // get all zones that surface point touches to all nodes
    set<int>* spZoneSet = new set<int>[nno_g];
    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np_recv = boundaryDataVec[ii].np_recv;
      for (int ip_recv = 0; ip_recv < np_recv; ++ip_recv) {
        const int ino = boundaryDataVec[ii].recv_ino[ip_recv];
        spZoneSet[ino].insert(boundaryDataVec[ii].zone);
      }
    }

    int* send_count = new int[mpi_size];
    int* send_disp = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;
    int* send_buf = new int[2*(nno_g-nno)];
    for (int iter = 0; iter < 2; ++iter) {
      for (int ino = nno; ino < nno_g; ++ino) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[ino-nno]);
        if (iter == 0) {
          send_count[rank] += 2;
        }
        else {
          send_buf[send_disp[rank]++] = index;
          send_buf[send_disp[rank]++] = ino;
        }
      }
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    }

    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    int * recv_buf = new int[recv_count_sum];
    MPI_Alltoallv(send_buf,send_count,send_disp,MPI_INT,
                  recv_buf,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf;

    FOR_RANK send_count[rank] = 0;
    for (int iter = 0; iter < 2; ++iter) {
      FOR_RANK {
        int irecv = recv_disp[rank];
        while (irecv < recv_disp[rank]+recv_count[rank]) {
          const int ino = recv_buf[irecv++];
          const int ino_nbr = recv_buf[irecv++];
          if (!spZoneSet[ino].empty()) {
            if (iter == 0) {
              send_count[rank] += 2+spZoneSet[ino].size();
            }
            else {
              send_buf[send_disp[rank]++] = ino_nbr;
              send_buf[send_disp[rank]++] = spZoneSet[ino].size();
              for (set<int>::iterator iter = spZoneSet[ino].begin(); iter != spZoneSet[ino].end(); ++iter)
                send_buf[send_disp[rank]++] = *iter;
            }
          }
        }
        assert(irecv == (recv_disp[rank]+recv_count[rank]));
      }
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
      if (iter == 0)
        send_buf = new int[send_disp[mpi_size-1]+send_count[mpi_size-1]];
    }
    delete[] recv_buf;

    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    recv_buf = new int[recv_count_sum];
    MPI_Alltoallv(send_buf,send_count,send_disp,MPI_INT,
                  recv_buf,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_count;
    delete[] send_disp;
    delete[] send_buf;

    {
      int irecv = 0;
      while (irecv < recv_count_sum) {
        const int ino = recv_buf[irecv++]; assert((ino >= nno)&&(ino < nno_g));
        const int count = recv_buf[irecv++];
        for (int ii = 0; ii < count; ++ii)
          spZoneSet[ino].insert(recv_buf[irecv++]);
      }
      assert(irecv == recv_count_sum);
    }
    delete[] recv_count;
    delete[] recv_disp;
    delete[] recv_buf;

    // build stono_i/v...
    int * stono_i = new int[nno_g+1];
    for (int ino = 0; ino < nno_g; ++ino)
      stono_i[ino+1] = 0;
    for (int ist = 0; ist < nst_ib; ++ist) {
      FOR_I3 {
        const int ino = noost[ist][i];
        ++stono_i[ino+1];
      }
    }
    stono_i[0] = 0;
    for (int ino = 0; ino < nno_g; ++ino)
      stono_i[ino+1] += stono_i[ino];
    const int stono_s = stono_i[nno_g];
    int * stono_v = new int[stono_s];
    for (int ist = 0; ist < nst_ib; ++ist) {
      FOR_I3 {
        const int ino = noost[ist][i];
        stono_v[stono_i[ino]++] = ist;
      }
    }
    for (int ino = nno_g-1; ino > 0; --ino)
      stono_i[ino] = stono_i[ino-1];
    stono_i[0] = 0;

    vector<int>* stZoneVec = new vector<int>[nst_ib];
    for (int ino = 0; ino < nno_g; ++ino) {
      for (int son = stono_i[ino]; son != stono_i[ino+1]; ++son) {
        const int ist = stono_v[son];
        for (set<int>::iterator iter = spZoneSet[ino].begin(); iter != spZoneSet[ino].end(); ++iter)
          stZoneVec[ist].push_back(*iter);
      }
    }
    delete[] spZoneSet;
    delete[] stono_i;
    delete[] stono_v;

    // set znost...
    for (int ist = 0; ist < nst_ib; ++ist) {
      if (stZoneVec[ist].size() < 3) {
        // just set beyond the number of bfZones...
        znost[ist] = nzones;
      }
      else {
        // find most often and set to znost...
        sort(stZoneVec[ist].begin(),stZoneVec[ist].end());
        int max_count = 1, curr_count = 1;
        znost[ist] = stZoneVec[ist][0];
        for (int ii = 1, lim = stZoneVec[ist].size(); ii < lim; ++ii) {
          if (stZoneVec[ist][ii] == stZoneVec[ist][ii-1]) {
            ++curr_count;
          }
          else {
            if (curr_count > max_count) {
              max_count = curr_count;
              znost[ist] = stZoneVec[ist][ii-1];
            }
            curr_count = 1;
          }
        }
        if (curr_count > max_count) {
          max_count = curr_count;
          znost[ist] = stZoneVec[ist][stZoneVec[ist].size()-1];
        }
        assert(max_count >= 3); // if we get here, then we should reset zone to nzones
      }
    }
    delete[] stZoneVec;

  }
  */

private:

  int getMatchingPoint(const double * const x) {
    // returns the point ino from the fem nodes, or -1 if no match
    // is found. tol should be set to return 1 and only one point...
    assert(adt); // also means tol has been computed
    const double bbmin[3] = { x[0]-tol, x[1]-tol, x[2]-tol };
    const double bbmax[3] = { x[0]+tol, x[1]+tol, x[2]+tol };
    vector<int> intVec;
    adt->buildListForBBox(intVec,bbmin,bbmax);
    if (!intVec.empty()) {
      assert(intVec.size() == 1);
      const int ino = intVec[0];
      //const double d2 = DIST2(x,x_no[ino]);
      //cout << "Got match on rank " << mpi_rank << " " << sqrt(d2) << endl;
      return ino;
    }
    return -1;
  }

  void buildCommStuff() {

    assert(recv_count);
    FOR_RANK recv_count[rank] = 0;

    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np_recv = boundaryDataVec[ii].np_recv;
      for (int ip_recv = 0; ip_recv < np_recv; ++ip_recv) {
        const int rank = boundaryDataVec[ii].recv_rank[ip_recv];
        ++recv_count[rank];
      }
    }

    assert(recv_disp);
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    no_buf_size = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    assert(no_buf_double == NULL); no_buf_double = new double[no_buf_size];

    // now populate noCommVec with the various counts...
    assert(noCommVec.empty());
    int ii = 0;
    FOR_RANK {
      if (recv_count[rank] > 0) {
        ++ii;
      }
    }
    noCommVec.resize(ii);
    ii = 0;
    FOR_RANK {
      if (recv_count[rank] > 0) {
        noCommVec[ii].rank = rank;
        noCommVec[ii].count = recv_count[rank];
        noCommVec[ii].disp = recv_disp[rank];
        ++ii;
      }
    }
    assert(ii == noCommVec.size());

    assert(noRequestArray == NULL);
    noRequestArray = new MPI_Request[noCommVec.size()];

    // each boundaryData member needs its recv_rank converted to
    // an index for direct access to the buffer...

    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np_recv = boundaryDataVec[ii].np_recv;
      for (int ip_recv = 0; ip_recv < np_recv; ++ip_recv) {
        const int rank = boundaryDataVec[ii].recv_rank[ip_recv];
        //const int ino = boundaryDataVec[ii].recv_ino[ip_recv];
        boundaryDataVec[ii].recv_rank[ip_recv] = recv_disp[rank]++;
      }
    }

    // and the sp (i.e. fluid) side...

    assert(send_count);
    FOR_RANK send_count[rank] = 0;

    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np = boundaryDataVec[ii].np;
      for (int ip = 0; ip < np; ++ip) {
        const int rank = boundaryDataVec[ii].send_rank[ip];
        ++send_count[rank];
      }
    }

    assert(send_disp);
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    sp_buf_size = send_disp[mpi_size-1] + send_count[mpi_size-1];

    assert(sp_buf_double == NULL); sp_buf_double = new double[sp_buf_size];

    // now populate noCommVec with the various counts...
    assert(spCommVec.empty());
    ii = 0;
    FOR_RANK {
      if (send_count[rank] > 0) {
        ++ii;
      }
    }
    spCommVec.resize(ii);
    ii = 0;
    FOR_RANK {
      if (send_count[rank] > 0) {
        spCommVec[ii].rank = rank;
        spCommVec[ii].count = send_count[rank];
        spCommVec[ii].disp = send_disp[rank];
        ++ii;
      }
    }
    assert(ii == spCommVec.size());

    assert(spRequestArray == NULL);
    spRequestArray = new MPI_Request[spCommVec.size()];

    // check on send_count...
    // TODO: remove this after everything is working
    int * send_count_check = new int[mpi_size];
    MPI_Alltoall(recv_count,1,MPI_INT,send_count_check,1,MPI_INT,mpi_comm);
    FOR_RANK assert(send_count_check[rank] == send_count[rank]);
    delete[] send_count_check;

    // unpack the send buf into the appropriate boundary data locations
    // on the fluid side...

    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np = boundaryDataVec[ii].np;
      for (int ip = 0; ip < np; ++ip) {
        const int rank = boundaryDataVec[ii].send_rank[ip];
        boundaryDataVec[ii].send_rank[ip] = send_disp[rank]++;
      }
    }

    // now delete the send_count, recv_count, etc...
    delete[] send_count; send_count = NULL;
    delete[] send_disp;  send_disp = NULL;
    delete[] recv_count; recv_count = NULL;
    delete[] recv_disp;  recv_disp = NULL;

    // check complexity...
    int my_buf[2] = { (int)noCommVec.size(), (int)spCommVec.size() };
    int buf[2];
    MPI_Reduce(my_buf,buf,2,MPI_INT,MPI_MAX,0,mpi_comm);
    if (mpi_rank == 0) cout << "Cht::buildCommStuff(): max noCommVec.size(): " << buf[0] << " spCommVec.size(): " << buf[1] << endl;

  }

public:

  void updateTStart() {

    // =======================================================
    // send T_no from the FEM solver to the sp's of the
    // flow solver bf's...
    // =======================================================

    // fem data is on the recv side...

    if (mode == INIT) {
      buildCommStuff();
      // also delete adt, bbox_adt...
      assert(adt);
      delete adt; adt = NULL;
      assert(bbox_adt);
      delete bbox_adt; bbox_adt = NULL;
    }
    else {
      assert(mode == QFINISH);
    }

    // the sp side is the recv side this time...

    assert(sp_buf_double);
    assert(spRequestArray);
    for (int ii = 0; ii < spCommVec.size(); ++ii) {
      // post irecv -- directly into sp_buf...
      MPI_Irecv(sp_buf_double+spCommVec[ii].disp,spCommVec[ii].count,
                MPI_DOUBLE,spCommVec[ii].rank,19671,mpi_comm,spRequestArray+ii);
    }

    // pack and send the node data: i.e. temperature...

    assert(no_buf_double);
    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np_recv = boundaryDataVec[ii].np_recv;
      for (int ip_recv = 0; ip_recv < np_recv; ++ip_recv) {
        const int ibuf = boundaryDataVec[ii].recv_rank[ip_recv];
        const int ino = boundaryDataVec[ii].recv_ino[ip_recv];
        no_buf_double[ibuf] = T_no[ino];
      }
    }

    assert(noRequestArray);
    for (int ii = 0; ii < noCommVec.size(); ++ii) {
      // post irecv -- directly into sp_buf...
      MPI_Issend(no_buf_double+noCommVec[ii].disp,noCommVec[ii].count,
                 MPI_DOUBLE,noCommVec[ii].rank,19671,mpi_comm,noRequestArray+ii);
    }

    mode = TSTART;

  }

  void updateTFinish() {

    assert(mode == TSTART);

    // now we wait for all messages to be sent and received...

    MPI_Waitall(noCommVec.size(),noRequestArray,MPI_STATUSES_IGNORE);

    // and finally, wait for the recv's...

    MPI_Waitall(spCommVec.size(),spRequestArray,MPI_STATUSES_IGNORE);

    // and unpack...

    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np = boundaryDataVec[ii].np;
      double * Tp = boundaryDataVec[ii].Tp;
      for (int ip = 0; ip < np; ++ip) {
        const int ibuf = boundaryDataVec[ii].send_rank[ip];
        Tp[ip] = sp_buf_double[ibuf];
      }
    }

    mode = TFINISH;

  }

  void updateT() {

    updateTStart();
    updateTFinish();

  }

  void updateChtTemperatureToFlowBcsOld() {

    // this is the original blocking routine. You cannot have
    // called buildCommStuff because it modifies the recv_rank
    // arrays...

    // fem data is on the recv side...
    // for now, use all-to-all communication...

    assert(recv_count);
    FOR_RANK recv_count[rank] = 0;

    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np_recv = boundaryDataVec[ii].np_recv;
      for (int ip_recv = 0; ip_recv < np_recv; ++ip_recv) {
        const int rank = boundaryDataVec[ii].recv_rank[ip_recv];
        ++recv_count[rank];
      }
    }

    assert(recv_disp);
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    // and pack the data...

    double * recv_buf_double = new double[recv_count_sum];
    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np_recv = boundaryDataVec[ii].np_recv;
      for (int ip_recv = 0; ip_recv < np_recv; ++ip_recv) {
        const int rank = boundaryDataVec[ii].recv_rank[ip_recv];
        const int ino = boundaryDataVec[ii].recv_ino[ip_recv];
        recv_buf_double[recv_disp[rank]++] = T_no[ino];
      }
    }

    // and rewind recv_disp...
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

    // prepare the send side (which will be RECIEVING the temperature -- sorry
    // for any confusion ;)

    assert(send_count);
    FOR_RANK send_count[rank] = 0;

    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np = boundaryDataVec[ii].np;
      for (int ip = 0; ip < np; ++ip) {
        const int rank = boundaryDataVec[ii].send_rank[ip];
        ++send_count[rank];
      }
    }

    assert(send_disp);
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

    // check on send_count...
    // TODO: remove this after everything working
    int * send_count_check = new int[mpi_size];
    MPI_Alltoall(recv_count,1,MPI_INT,send_count_check,1,MPI_INT,mpi_comm);
    FOR_RANK assert(send_count_check[rank] == send_count[rank]);
    delete[] send_count_check;

    double * send_buf_double = new double[send_count_sum];
    MPI_Alltoallv(recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,
                  send_buf_double,send_count,send_disp,MPI_DOUBLE,
                  mpi_comm);
    delete[] recv_buf_double;

    // unpack the send buf into the appropriate boundary data locations
    // on the fluid side...

    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np = boundaryDataVec[ii].np;
      double * Tp = boundaryDataVec[ii].Tp;
      for (int ip = 0; ip < np; ++ip) {
        const int rank = boundaryDataVec[ii].send_rank[ip];
        Tp[ip] = send_buf_double[send_disp[rank]++];
      }
    }
    delete[] send_buf_double;

  }

  void updateQStart() {

    // =======================================================
    // send q_sp from the bf's of the flow solver to the
    // q_no source of the FEM solver...
    // =======================================================

    assert(mode == TFINISH);

    // the no side is the recv side this time...

    assert(no_buf_double);
    assert(noRequestArray);
    for (int ii = 0; ii < noCommVec.size(); ++ii) {
      // post irecv...
      MPI_Irecv(no_buf_double+noCommVec[ii].disp,noCommVec[ii].count,
                MPI_DOUBLE,noCommVec[ii].rank,19672,mpi_comm,noRequestArray+ii);
    }

    // pack and send the heat flux data...

    assert(sp_buf_double);
    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np = boundaryDataVec[ii].np;
      double * qp = boundaryDataVec[ii].qp;
      for (int ip = 0; ip < np; ++ip) {
        const int ibuf = boundaryDataVec[ii].send_rank[ip];
        sp_buf_double[ibuf] = qp[ip];
      }
    }

    assert(spRequestArray);
    for (int ii = 0; ii < spCommVec.size(); ++ii) {
      // post irecv...
      MPI_Issend(sp_buf_double+spCommVec[ii].disp,spCommVec[ii].count,
                 MPI_DOUBLE,spCommVec[ii].rank,19672,mpi_comm,spRequestArray+ii);
    }

    mode = QSTART;

  }

  void updateQFinish() {

    assert(mode == QSTART);

    // now we wait for all messages to be sent and received...

    MPI_Waitall(spCommVec.size(),spRequestArray,MPI_STATUSES_IGNORE);

    // and finally, wait for the recv's...

    MPI_Waitall(noCommVec.size(),noRequestArray,MPI_STATUSES_IGNORE);

    // and unpack...

    FOR_INO q_no[ino] = 0.0;
    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np_recv = boundaryDataVec[ii].np_recv;
      for (int ip_recv = 0; ip_recv < np_recv; ++ip_recv) {
        const int ibuf = boundaryDataVec[ii].recv_rank[ip_recv];
        const int ino = boundaryDataVec[ii].recv_ino[ip_recv];
        q_no[ino] += no_buf_double[ibuf];
      }
    }

    mode = QFINISH;

  }

  void updateQ() {

    updateQStart();
    updateQFinish();

  }

  void updateFlowBcQToChtOld() {

    // prepare the send side (which will be RECIEVING the temperature -- sorry
    // for any confusion ;)

    assert(send_count);
    FOR_RANK send_count[rank] = 0;

    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np = boundaryDataVec[ii].np;
      for (int ip = 0; ip < np; ++ip) {
        const int rank = boundaryDataVec[ii].send_rank[ip];
        ++send_count[rank];
      }
    }

    assert(send_disp);
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

    // unpack the send buf into the appropriate boundary data locations
    // on the fluid side...

    double * send_buf_double = new double[send_count_sum];
    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np = boundaryDataVec[ii].np;
      double * qp = boundaryDataVec[ii].qp;
      for (int ip = 0; ip < np; ++ip) {
        const int rank = boundaryDataVec[ii].send_rank[ip];
        send_buf_double[send_disp[rank]++] = qp[ip];
      }
    }

    // rewind disp...
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

    // prepare the recv side stuff...

    assert(recv_count);
    FOR_RANK recv_count[rank] = 0;

    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np_recv = boundaryDataVec[ii].np_recv;
      for (int ip_recv = 0; ip_recv < np_recv; ++ip_recv) {
        const int rank = boundaryDataVec[ii].recv_rank[ip_recv];
        ++recv_count[rank];
      }
    }

    assert(recv_disp);
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    // check on recv_count...
    // TODO: remove this after everything working
    int * recv_count_check = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count_check,1,MPI_INT,mpi_comm);
    FOR_RANK assert(recv_count_check[rank] == recv_count[rank]);
    delete[] recv_count_check;

    // exchange...
    double * recv_buf_double = new double[recv_count_sum];
    MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
                  recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,
                  mpi_comm);
    delete[] send_buf_double;

    // and unpack the data into the fem source term...

    FOR_INO q_no[ino] = 0.0;
    for (int ii = 0; ii < boundaryDataVec.size(); ++ii) {
      const int np_recv = boundaryDataVec[ii].np_recv;
      for (int ip_recv = 0; ip_recv < np_recv; ++ip_recv) {
        const int rank = boundaryDataVec[ii].recv_rank[ip_recv];
        const int ino = boundaryDataVec[ii].recv_ino[ip_recv];
        q_no[ino] += recv_buf_double[recv_disp[rank]++];
      }
    }
    delete[] recv_buf_double;

  }

  void writeCht(const int step) {
    char filename[128];
    sprintf(filename,"%s.%08d.cht",write_cht_prefix.c_str(),step);
    if (mpi_rank == 0) cout << " Writing cht file \"" << filename << "\"..." << endl;
    writeData(filename);
    if (mpi_rank == 0) {
      unlink("cht");
      symlink(filename,"cht");
    }
    MPI_Barrier(mpi_comm);
  }

  void processWriteCht(Param* param,const int step,const bool killfile_request) {

    if (param->size() == 0) {
      // ----------------------------------------------------
      // WRITE_CHT
      // ----------------------------------------------------
      if (killfile_request) {
        writeCht(step);
        return;
      }
      else {
        CWARN(" > WRITE_CHT syntax should be:\n" <<
        "WRITE_CHT INTERVAL <interval> [NAME <prefix>] or \n" <<
        "WRITE_CHT <interval>");
        if (write_cht_interval == -1) write_cht_interval = 0; // means no WRITE_CHT
        return;
      }
    }
    else if (param->size() == 1) {
      // ----------------------------------------------------
      // WRITE_CHT 1000
      // ----------------------------------------------------
      // if the interval has never been set, then reset it...
      const int interval = param->getInt(0);
      if (interval <= 0) {
        CWARN(" > WRITE_CHT <interval> expects positive integer");
        if (write_cht_interval == -1) write_cht_interval = 0; // means no WRITE_CHT
        return;
      }
      write_cht_interval = interval;
      // if this happens to be from a killfile, then write the current
      // cht filename if the step is divisible by write_cht_interval. Otherwise
      // this just sets the step...
      if (killfile_request && (step%write_cht_interval == 0))
      writeCht(step);
      return;
    }
    else {
      // ----------------------------------------------------
      // WRITE_CHT NAME blah INTERVAL 1000
      // ----------------------------------------------------
      string name;
      bool b_name = false;
      int interval;
      bool b_interval = false;
      int iarg = 0;
      while (iarg < param->size()) {
        const string token = param->getString(iarg++);
        if (token == "NAME") {
          name = param->getString(iarg++);
          b_name = true;
        }
        else if (token == "INTERVAL") {
          interval = param->getInt(iarg++);
          b_interval = true;
        }
        else {
          CWARN(" > unrecognized token in WRITE_CHT: " << token << "; skipping...");
        }
      }
      if ((!b_interval)||(interval <= 0)) {
        CWARN(" > check WRITE_CHT syntax:\n" <<
        "WRITE_CHT INTERVAL <interval> [NAME <prefix>] or \n" <<
        "WRITE_CHT <interval>");
        if (write_cht_interval == -1) write_cht_interval = 0; // means no WRITE_CHT
        return;
      }
      write_cht_interval = interval;
      if (b_name) write_cht_prefix = name;
      if (killfile_request && (step%write_cht_interval == 0))
      writeCht(step);
      return;
    }

  }

  void init() {

    if (mpi_rank == 0) cout << "Cht::init()..." << endl;

    int8 * noora = NULL;
    MiscUtils::buildXora(noora,nno);

    assert(noora[0] == 0);
    assert(noora[mpi_rank+1]-noora[mpi_rank] == nno);

    int8 *ino_global2 = new int8[nno];
    for (int ino = 0; ino < nno; ++ino) ino_global2[ino] = noora[mpi_rank] + ino; // global no index -- need this

    int nno_orig = nno;
    repartXcvPadt(x_no,ino_global2,nno,mpi_comm);

    // x_no should be balanced now. This reorder puts nearby data closer in memory...
    reorderXcv(x_no,ino_global2,nno); // note: new nno

    DistributedDataExchanger* dde = new DistributedDataExchanger(ino_global2,nno,noora);
    delete[] ino_global2; ino_global2 = NULL;

    // rbi == "rank-bits-index"...
    uint8* rbi = new uint8[nno];
    for (int ino = 0; ino < nno; ++ino)
      rbi[ino] = BitUtils::packRankBitsIndex(mpi_rank,0,ino);

    uint8 * rbi_orig = new uint8[nno_orig];
    dde->push(rbi_orig,rbi);
    delete[] rbi;
    delete dde; dde = NULL;

    // we now have the rank and index to send all nno_orig stuff in rbi_orig...
    // now pull the rbi_orig out to the tets...

    assert(ino_global2 == NULL); ino_global2 = new int8[nte*4];
    for (int ite = 0; ite < nte; ++ite)
      FOR_I4 ino_global2[ite*4+i] = noote[ite][i];
    delete[] noote; noote = NULL;

    assert(dde == NULL); dde = new DistributedDataExchanger(ino_global2,nte*4,noora);
    delete[] ino_global2; ino_global2 = NULL;
    delete[] noora; // reused by nst below if present

    uint8 (*rbiote)[4] = new uint8[nte][4];
    dde->pull((uint8*)rbiote,rbi_orig);
    delete dde; dde = NULL;

    // and send the tets as 4 rbi's in the final distribution. This
    // also set the ghost data. We also use this exchange to send the
    // tet volume zone, it present. Otherwise, we send zero...

    if (send_count == NULL) send_count = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;

    // here we temporarily use recv_count as a flag...
    if (recv_count == NULL) recv_count = new int[mpi_size];
    FOR_RANK recv_count[rank] = -1;
    for (int ite = 0; ite < nte; ++ite) {
      FOR_I4 {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbiote[ite][i]);
        if (recv_count[rank] != ite) {
          send_count[rank] += 5;
          recv_count[rank] = ite;
        }
      }
    }

    if (send_disp == NULL) send_disp = new int[mpi_size];
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
    int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

    uint8 * send_buf_uint8 = new uint8[send_count_sum];
    for (int ite = 0; ite < nte; ++ite) {
      FOR_I4 {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbiote[ite][i]);
        if (recv_count[rank] != nte+ite) { // use nte offset here to save resetting flag to -1 ;)
          FOR_J4 send_buf_uint8[send_disp[rank]+j] = rbiote[ite][j];
	  send_buf_uint8[send_disp[rank]+4] = ((znote == NULL) ? 0 : znote[ite]);
          send_disp[rank] += 5;
          recv_count[rank] = nte+ite;
        }
      }
    }
    delete[] rbiote;
    delete[] znote; znote = NULL;

    // rewind...
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    if (recv_disp == NULL) recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    assert(recv_count_sum%5 == 0);

    uint8 * recv_buf_uint8 = new uint8[recv_count_sum];
    MPI_Alltoallv(send_buf_uint8,send_count,send_disp,MPI_UINT8,
                  recv_buf_uint8,recv_count,recv_disp,MPI_UINT8,mpi_comm);
    delete[] send_buf_uint8; send_buf_uint8 = NULL;

    map<const uint8,int> rbiMap;
    nno_g = nno; // recall nno stores the new node distribution
    nte = recv_count_sum/5;
    int * flag_te = new int[nte];
    FOR_ITE flag_te[ite] = 0;
    assert(noote == NULL); noote = new int[nte][4];
    assert(znote == NULL); znote = new int[nte];
    for (int irecv = 0; irecv < recv_count_sum; irecv += 5) {
      const int ite = irecv/5;
      FOR_I4 {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,recv_buf_uint8[irecv+i]);
        if ((rank == mpi_rank)&&(bits == 0)) {
          assert((index >= 0)&&(index < nno));
          noote[ite][i] = index;
        }
        else {
          // mark this tet as a boundary tet:
          // 2: a lower rank also has this tet
          // 1: we are the lowest rank with this tet...
          if (rank < mpi_rank) flag_te[ite] = 2;
          else flag_te[ite] = max(flag_te[ite],1);
          map<const uint8,int>::iterator iter = rbiMap.find(recv_buf_uint8[irecv+i]);
          if (iter == rbiMap.end()) {
            // this is the first time we have seen this node...
            rbiMap[recv_buf_uint8[irecv+i]] = noote[ite][i] = nno_g++;
          }
          else {
            // this node is already indexed...
            noote[ite][i] = iter->second;
          }
        }
      }
      znote[ite] = recv_buf_uint8[irecv+4];
    }
    delete[] recv_buf_uint8; recv_buf_uint8 = NULL;

    // reorder the tets in 3 ranges:
    // 0: internal,
    // 1: internal-boundary-owned,
    // 2: internal-boundary-owned-by-another-rank

    int counts[3] = { 0, 0, 0 };
    FOR_ITE ++counts[flag_te[ite]];

    int nte_old = nte;
    int (*noote_old)[4] = noote;
    noote = new int[nte][4];
    int *znote_old = znote;
    znote = new int[nte];
    nte_i = 0;
    nte_ib = counts[0];
    nte = counts[0]+counts[1];
    // also we will be flagging the nodes for local reorder...
    int * order = new int[nno];
    FOR_INO order[ino] = 0;
    for (int ite = 0; ite < nte_old; ++ite) {
      switch (flag_te[ite]) {
      case 0:
        FOR_I4 noote[nte_i][i] = noote_old[ite][i];
	znote[nte_i] = znote_old[ite];
        ++nte_i;
        break;
      case 1:
        FOR_I4 {
          noote[nte_ib][i] = noote_old[ite][i];
          if (noote_old[ite][i] < nno) {
            order[noote_old[ite][i]] = 1;
          }
        }
	znote[nte_ib] = znote_old[ite];
        ++nte_ib;
        break;
      case 2:
        FOR_I4 {
          noote[nte][i] = noote_old[ite][i];
          if (noote_old[ite][i] < nno) {
            order[noote_old[ite][i]] = 1;
          }
        }
	znote[nte] = znote_old[ite];
        ++nte;
        break;
      default:
        assert(0);
      }
    }
    assert(nte == nte_old);
    delete[] noote_old;
    delete[] znote_old;
    delete[] flag_te;

    // ===================================
    // load in the properties...
    assert(rhoFuncVec.empty());
    assert(cpFuncVec.empty());
    assert(kFuncVec.empty());
    int my_nzn = 0;
    for (int ite = 0; ite < nte; ++ite) {
      assert(znote[ite] >= 0);
      my_nzn = max(my_nzn,znote[ite]+1);
    }
    int nzn;
    MPI_Allreduce(&my_nzn,&nzn,1,MPI_INT,MPI_MAX,mpi_comm);
    if (nzn == 1) {
      // for one zone only, look for
      // CHT.RHO, CHT.CP, CHT.K
      if (mpi_rank == 0) cout << " > CHT.RHO..." << endl;
      rhoFuncVec.push_back(processSimpleFunc(getParam("CHT.RHO")));
      if (mpi_rank == 0) cout << " > CHT.CP..." << endl;
      cpFuncVec.push_back(processSimpleFunc(getParam("CHT.CP")));
      if (mpi_rank == 0) cout << " > CHT.K..." << endl;
      kFuncVec.push_back(processSimpleFunc(getParam("CHT.K")));
    }
    else {
      // for multiple zones, use indexing...
      // CHT.0.RHO, CHT.0.CP, CHT.0.K
      // CHT.1.RHO, CHT.1.CP, CHT.1.K
      char name[128];
      for (int izn = 0; izn < nzn; ++izn) {
	sprintf(name,"CHT.%d.RHO",izn);
	if (mpi_rank == 0) cout << " > " << name << "..." << endl;
	rhoFuncVec.push_back(processSimpleFunc(getParam(name)));
	sprintf(name,"CHT.%d.CP",izn);
	if (mpi_rank == 0) cout << " > " << name << "..." << endl;
	cpFuncVec.push_back(processSimpleFunc(getParam(name)));
	sprintf(name,"CHT.%d.K",izn);
	if (mpi_rank == 0) cout << " > " << name << "..." << endl;
	kFuncVec.push_back(processSimpleFunc(getParam(name)));
      }
    }

    // now pull the rbi_orig out to the surface tri's...
    /*
    assert(ino_global2 == NULL); ino_global2 = new int8[nst*3];
    for (int ist = 0; ist < nst; ++ist)
      FOR_I3 ino_global2[ist*3+i] = noost[ist][i];
    delete[] noost; noost = NULL;

    assert(dde == NULL); dde = new DistributedDataExchanger(ino_global2,nst*3,noora);
    delete[] ino_global2;
    delete[] noora;

    uint8 (*rbiost)[3] = new uint8[nst][3];
    dde->pull((uint8*)rbiost,rbi_orig);
    delete dde;

    // and send the surface tri's as 3 rbi's in the final distribution...

    FOR_RANK send_count[rank] = 0;

    // here we temporarily use recv_count as a flag...
    FOR_RANK recv_count[rank] = -1;
    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbiost[ist][i]);
        if (recv_count[rank] != ist) {
          send_count[rank] += 3;
          recv_count[rank] = ist;
        }
      }
    }

    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
    send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

    assert(send_buf_uint8 == NULL); send_buf_uint8 = new uint8[send_count_sum];
    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbiost[ist][i]);
        if (recv_count[rank] != nst+ist) { // use nst offset here to save resetting flag to -1 ;)
          FOR_J3 send_buf_uint8[send_disp[rank]+j] = rbiost[ist][j];
          send_disp[rank] += 3;
          recv_count[rank] = nst+ist;
        }
      }
    }
    delete[] rbiost;

    // rewind...
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    assert(recv_count_sum%3 == 0);

    assert(recv_buf_uint8 == NULL); recv_buf_uint8 = new uint8[recv_count_sum];
    MPI_Alltoallv(send_buf_uint8,send_count,send_disp,MPI_UINT8,
                  recv_buf_uint8,recv_count,recv_disp,MPI_UINT8,mpi_comm);
    delete[] send_buf_uint8;

    nst = recv_count_sum/3;
    int * flag_st = new int[nst];
    FOR_IST flag_st[ist] = 0;
    assert(noost == NULL); noost = new int[nst][3];
    for (int irecv = 0; irecv < recv_count_sum; irecv += 3) {
      const int ist = irecv/3;
      FOR_I3 {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,recv_buf_uint8[irecv+i]);
        if ((rank == mpi_rank)&&(bits == 0)) {
          assert((index >= 0)&&(index < nno));
          noost[ist][i] = index;
        }
        else {
          // mark this tri as a boundary tri:
          // 2: a lower rank also has this tri
          // 1: we are the lowest rank with this tri...
          if (rank < mpi_rank) flag_st[ist] = 2;
          else flag_st[ist] = max(flag_st[ist],1);
          map<const uint8,int>::iterator iter = rbiMap.find(recv_buf_uint8[irecv+i]);
          // this node is already indexed...
          assert(iter != rbiMap.end());
          noost[ist][i] = iter->second;
        }
      }
    }
    delete[] recv_buf_uint8;

    // reorder the surface tri's in 3 ranges:
    // 0: internal,
    // 1: internal-boundary-owned,
    // 2: internal-boundary-owned-by-another-rank

    FOR_I3 counts[i] = 0;
    FOR_IST ++counts[flag_st[ist]];

    int nst_old = nst;
    int (*noost_old)[3] = noost;
    noost = new int[nst][3];
    nst_i = 0;
    nst_ib = counts[0];
    nst = counts[0]+counts[1];
    for (int ist = 0; ist < nst_old; ++ist) {
      switch (flag_st[ist]) {
      case 0:
        FOR_I3 noost[nst_i][i] = noost_old[ist][i];
        ++nst_i;
        break;
      case 1:
        FOR_I3 noost[nst_ib][i] = noost_old[ist][i];
        ++nst_ib;
        break;
      case 2:
        FOR_I3 noost[nst][i] = noost_old[ist][i];
        ++nst;
        break;
      default:
        assert(0);
      }
    }
    assert(nst == nst_old);
    delete[] noost_old;
    delete[] flag_st;
    */

    // at this point, the active nodes (ino < nno) on the boundary are flagged with
    // order == 1, and order == 0 elsewhere...
    counts[0] = 0;
    counts[1] = 0;
    FOR_INO ++counts[order[ino]];

    int nno_old = nno;
    double (*x_no_old)[3] = x_no;
    x_no = new double[nno_g][3];
    nno_i = 0;
    nno = counts[0];
    for (int ino = 0; ino < nno_old; ++ino) {
      switch(order[ino]) {
      case 0:
        FOR_I3 x_no[nno_i][i] = x_no_old[ino][i];
        order[ino] = nno_i++;
        break;
      case 1:
        FOR_I3 x_no[nno][i] = x_no_old[ino][i];
        order[ino] = nno++;
        break;
      default:
        assert(0);
      }
    }
    assert(nno == nno_old);
    delete[] x_no_old;

    // reorder ghost data so it is contiguous in rank/b/index. This is exactly the current
    // order of the map...

    int * order_g = new int[nno_g-nno];
    assert(rbi_g == NULL); rbi_g = new uint8[nno_g-nno];
    int ino_g = nno;
    for (map<const uint8,int>::iterator iter = rbiMap.begin(); iter != rbiMap.end(); ++iter) {
      order_g[iter->second-nno] = ino_g;
      rbi_g[ino_g-nno] = iter->first;
      ++ino_g;
    }
    assert(ino_g == nno_g);
    rbiMap.clear();

    // change the ghost order in noote and noost...
    for (int ite = 0; ite < nte; ++ite) {
      FOR_I4 {
        const int ino = noote[ite][i];
        if (ino < nno) {
          noote[ite][i] = order[ino];
        }
        else {
          noote[ite][i] = order_g[ino-nno];
        }
      }
    }
    /*
    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
        const int ino = noost[ist][i];
        if (ino < nno) {
          noost[ist][i] = order[ino];
        }
        else {
          noost[ist][i] = order_g[ino-nno];
        }
      }
    }
    */
    delete[] order_g;

    // and finally use rbi_g to build the node communicator...
    // TODO: put this in a common, protected location...

    buildPrcomm(noPrcommVec,nno,nno_g,rbi_g);
    reorderPrcommPackVec(noPrcommVec,order,nno);

    // reorder the index in rbi_g...
    {

      int* send_count = new int[mpi_size];
      int* send_disp = new int[mpi_size];
      FOR_RANK send_count[rank] = 0;
      int* send_buf = new int[nno_g-nno];
      for (int iter = 0; iter < 2; ++iter) {
        for (int ino = nno; ino < nno_g; ++ino) {
          int rank,bits,index;
          BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[ino-nno]);
          if (iter == 0)
            ++send_count[rank];
          else
            send_buf[send_disp[rank]++] = index;
        }
        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
      }

      int * recv_count = new int[mpi_size];
      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      int * recv_disp = new int[mpi_size];
      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      int * recv_buf = new int[recv_count_sum];
      MPI_Alltoallv(send_buf,send_count,send_disp,MPI_INT,
                    recv_buf,recv_count,recv_disp,MPI_INT,mpi_comm);

      for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
        const int ino0 = recv_buf[irecv];
        recv_buf[irecv] = order[ino0];
      }

      MPI_Alltoallv(recv_buf,recv_count,recv_disp,MPI_INT,
                    send_buf,send_count,send_disp,MPI_INT,mpi_comm);
      delete[] recv_buf;
      delete[] recv_count;
      delete[] recv_disp;

      for (int ino = nno; ino < nno_g; ++ino) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[ino-nno]);
        rbi_g[ino-nno] = BitUtils::packRankBitsIndex(rank,bits,send_buf[send_disp[rank]++]);
      }
      delete[] send_buf;
      delete[] send_count;
      delete[] send_disp;

    }

    // update x_no into ghost in prep for operator build...

    updateNoData(x_no);

    // finally, bring over T_no and ino_global from the original global index striping....

    FOR_RANK send_count[rank] = 0;
    for (int ino = 0; ino < nno_orig; ++ino) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_orig[ino]);
      ++send_count[rank];
    }
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
    const int send_count_sum2 = send_disp[mpi_size-1] + send_count[mpi_size-1];
    int * send_buf_int = new int[send_count_sum2*2];
    double * send_buf_double = new double[send_count_sum2];
    for (int ino = 0; ino < nno_orig; ++ino) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_orig[ino]);
      assert(bits == 0);
      send_buf_int[send_disp[rank]*2  ] = index;
      send_buf_int[send_disp[rank]*2+1] = ino;
      send_buf_double[send_disp[rank]] = T_no[ino];
      ++send_disp[rank];
    }
    delete[] rbi_orig;
    if (T_no != NULL) delete[] T_no; T_no = NULL;
    // rewind...
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
    // prepare recv side...
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    const int recv_count_sum2 = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    int * recv_buf_int = new int[2*recv_count_sum2];
    FOR_RANK {
      send_count[rank] *= 2;
      send_disp[rank] *= 2;
      recv_count[rank] *= 2;
      recv_disp[rank] *= 2;
    }
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
                  recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int;
    FOR_RANK {
      send_count[rank] /= 2;
      send_disp[rank] /= 2;
      recv_count[rank] /= 2;
      recv_disp[rank] /= 2;
    }

    double * recv_buf_double = new double[recv_count_sum2];
    MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
                  recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
    delete[] send_buf_double;
    // now unpack into T_no...
    T_no = new double[nno_g];
    assert(ino_global == NULL); ino_global = new int8[nno]; // need ghosts?
    assert(noora_striped != NULL);
    for (int ino = 0; ino < nno_g; ++ino) T_no[ino] = HUGE_VAL;
    FOR_RANK {
      int irecv = recv_disp[rank];
      while (irecv < recv_disp[rank]+recv_count[rank]) {
        const int ino = recv_buf_int[irecv*2  ]; assert((ino >= 0)&&(ino < nno));
        const int ino_striped = recv_buf_int[irecv*2+1]; assert((ino >= 0)&&(ino < nno));
        const int ino_new = order[ino]; assert((ino_new >= 0)&&(ino_new < nno));
        assert(T_no[ino_new] == HUGE_VAL);
        T_no[ino_new] = recv_buf_double[irecv];
        ino_global[ino_new] = ino_striped + noora_striped[rank];
        // HACK check: should be exact...
        assert(T_no[ino_new] == (x_no[ino_new][0]+x_no[ino_new][1]+x_no[ino_new][2]));
        ++irecv;
      }
    }
    delete[] recv_buf_int;
    delete[] recv_buf_double;
    delete[] order;
    assert(dde_striped == NULL); dde_striped = new DistributedDataExchanger(ino_global,nno,noora_striped);

    MiscUtils::dumpRange(&nno,1,"cht nodes");

    // ====================================================================
    // intial conditions...
    // ====================================================================

    FOR_INO T_no[ino] = HUGE_VAL;
    updateNoData(T_no);

    assert(q_no == NULL);
    q_no = new double[nno];
    FOR_INO q_no[ino] = HUGE_VAL;

    /*
    // some examples of serial file io leveraging ordering...
    double * var = new double[nno_g];
    for (int ino = 0; ino < nno_g; ++ino) var[ino] = x_no[ino][0] + x_no[ino][1] + x_no[ino][2]; // dummy scalar

    // write all tets and nodes...
    char filename[128];
    sprintf(filename,"tets.%08d.dat",mpi_rank);
    writeTecplotSerial(filename,var,x_no,nno_g,noote,nte);

    // write active tets and all nodes. Note some nodes may not have a tet, because
    // the division between active and inactive tets is not respected in the
    // ghost nodes...
    sprintf(filename,"tets-active.%08d.dat",mpi_rank);
    writeTecplotSerial(filename,var,x_no,nno_g,noote,nte_ib);

    // write internal tets and active nodes...
    sprintf(filename,"tets-internal.%08d.dat",mpi_rank);
    writeTecplotSerial(filename,var,x_no,nno,noote,nte_i);
    delete[] var;
    */

    // check on tet ordering...

    for (int ite = 0; ite < nte_i; ++ite) {
      FOR_I4 {
        const int ino = noote[ite][i];
        assert((ino >= 0)&&(ino < nno));
      }
    }

    buildNoono();
    buildOperators();

    assert(rhs_no[0] == NULL);
    rhs_no[0] = new double[nno*3];
    rhs_no[1] = rhs_no[0] + nno;
    rhs_no[2] = rhs_no[0] + nno*2;

    assert(A == NULL);
    A = new double[noono_i[nno]];

    /*
    assert(znost == NULL);
    znost = new int[nst];
    FOR_IST znost[ist] = -1; // unset
    */

  }

  void rk3Step(const double rk_wgt[3], const double dt, const int rk_stage) {

    // rk_stage is passed as 1,2,3...
    assert(rk_stage >= 1);
    assert(!rhoFuncVec.empty());
    assert(!cpFuncVec.empty());
    assert(!kFuncVec.empty());
    assert(vol_te);
    assert(grad_te);

    //if (rho_no == NULL) rho_no = new double[nno];
    //if (cp_no == NULL) cp_no = new double[nno];
    //if (k_no == NULL) k_no = new double[nno_g]; // only k needs to be valid in the ghost region

    if (rk_stage == 1) {
      // on the first stage, re-compute rho,cp,k as a funtion of the current T...
      // TODO: look at bcs that do something similar. Why is there a separate function?
      if (rho_cp_vol_no == NULL) rho_cp_vol_no = new double[nno];
      if (laplace_k_noono == NULL) laplace_k_noono = new double[noono_i[nno]];
      FOR_INO rho_cp_vol_no[ino] = 0.0;
      for (int non = 0; non < noono_i[nno]; ++non)
	laplace_k_noono[non] = 0.0;
      FOR_ITE {
	const int ino0 = noote[ite][0];
	const int ino1 = noote[ite][1];
	const int ino2 = noote[ite][2];
	const int ino3 = noote[ite][3];
	// compute the tet temperature: TODO: could eval these props a little
	// more efficiently if this is important...
	const double T_avg = 0.25*(T_no[ino0]+T_no[ino1]+T_no[ino2]+T_no[ino3]);
	double rho_te,cp_te,k_te;
	const int izn = ( (znote == NULL) ? 0 : znote[ite] );
	rhoFuncVec[izn]->eval(&rho_te,&T_avg,1);
	cpFuncVec[izn]->eval(&cp_te,&T_avg,1);
	kFuncVec[izn]->eval(&k_te,&T_avg,1);
	// distribute 1/4 to each node...
	if (ino0 < nno) rho_cp_vol_no[ino0] += 0.25*rho_te*cp_te*vol_te[ite];
	if (ino1 < nno) rho_cp_vol_no[ino1] += 0.25*rho_te*cp_te*vol_te[ite];
	if (ino2 < nno) rho_cp_vol_no[ino2] += 0.25*rho_te*cp_te*vol_te[ite];
	if (ino3 < nno) rho_cp_vol_no[ino3] += 0.25*rho_te*cp_te*vol_te[ite];
	// and build the laplace operator using this k_te...
	FOR_I4 {
	  const int ino = noote[ite][i];
	  if (ino < nno) {
	    const int non = noono_i[ino]; assert(noono_v[non] == ino);
	    FOR_J4 {
	      const int ino_nbr = noote[ite][j];
	      int non_nbr = non;
	      while (noono_v[non_nbr] != ino_nbr)
		++non_nbr;
	      assert(non_nbr < noono_i[ino+1]); // make sure we found ino_nbr...
	      // integration-by-parts flips the sign of this term...
	      laplace_k_noono[non_nbr] -= DOT_PRODUCT(grad_te[ite][i],grad_te[ite][j])*vol_te[ite]*k_te;
	    }
	  }
	}
      }
    }
    assert(rho_cp_vol_no);
    assert(laplace_k_noono);

    //if ((mpi_rank == 0)) cout << "Cht::rk3Step: " << COUT_VEC(rk_wgt) << " dt: " << dt << " rk_stage: " << rk_stage << endl;

    // use rhs_no[rk_stage-1] to store the rhs for this substep...
    // note: this is NOT the rhs that we need to store for subsequent rk substeps. This is
    // the rhs of the system:
    //
    // A*T_no = rhs
    //
    // the rhs req'd for subsequent rk substeps is reversed out below (see "really weird")...

    FOR_INO rhs_no[rk_stage-1][ino] = rk_wgt[rk_stage-1]*q_no[ino] + rho_cp_vol_no[ino]*alpha*T_no[ino]/dt;
    for (int i = 0; i < rk_stage-1; ++i) {
      assert(rk_wgt[i] != 0.0);
      FOR_INO rhs_no[rk_stage-1][ino] += rk_wgt[i]*rhs_no[i][ino];
    }

    for (int ino = 0; ino < nno; ++ino) {
      const int non_f = noono_i[ino];
      assert(noono_v[non_f] == ino);
      A[non_f] = rho_cp_vol_no[ino]*alpha/dt;
      for (int non = non_f+1; non != noono_i[ino+1]; ++non) {
        const int ino_nbr = noono_v[non];
        A[non] = -rk_wgt[rk_stage-1]*laplace_k_noono[non];
        A[non_f] -= A[non];
      }
    }

    // store T_no in q_no (q_no not needed anymore)...
    FOR_INO q_no[ino] = T_no[ino];

    // solve A*T = rhs
    // TODO: put this in a solver area eventually...
    const double T_zero = getDoubleParam("TEMPERATURE_ZERO",1.0E-8);
    const int T_maxiter = getIntParam("TEMPERATURE_MAXITER",2000);
    const bool T_verbose = getBoolParam("TEMPERATURE_VERBOSE",false);
    solveNoCg(T_no,A,rhs_no[rk_stage-1],T_zero,T_maxiter,T_verbose);

    // now store the final rhs for this level...
    // this is only necessary for stages 1 and 2...
    // this is "really weird", i.e. not great, but ensures conservation...
    if (rk_stage < 3) {
      FOR_INO rhs_no[rk_stage-1][ino] = rho_cp_vol_no[ino]*alpha*(T_no[ino]-q_no[ino])/dt; // note q contains T
      for (int i = 0; i < rk_stage-1; ++i) {
        assert(rk_wgt[i] != 0.0);
        FOR_INO rhs_no[rk_stage-1][ino] -= rk_wgt[i]*rhs_no[i][ino];
      }
      FOR_INO rhs_no[rk_stage-1][ino] /= rk_wgt[rk_stage-1];
    }

    //MiscUtils::dumpRange(T_no,nno,"T_no");

  }

  void writeTecplot(const int step) {

    // handshake version of ascii tecplot write...

    char filename[128];
    sprintf(filename,"T_no.%08d.dat",step);

    int my_buf[2] = { nno, nte_ib };
    int buf[2];
    MPI_Reduce(my_buf,buf,2,MPI_INT,MPI_SUM,0,mpi_comm);

    int offset;
    MPI_Scan(&nno,&offset,1,MPI_INT,MPI_SUM,mpi_comm);
    offset -= nno;

    int * flag_no = new int[nno_g];
    FOR_INO flag_no[ino] = ino + offset;
    updateNoData(flag_no);

    // ====================
    // header and nodes first...
    // ====================

    FILE * fp;
    if ( mpi_rank == 0 ) {
      fp = fopen(filename,"w");
      assert(fp != NULL);
      fprintf(fp,"TITLE = \"%s\"\n",filename);
      fprintf(fp,"VARIABLES = \"X\"\n");
      fprintf(fp,"\"Y\"\n");
      fprintf(fp,"\"Z\"\n");
      fprintf(fp,"\"VAR\"\n");
      fprintf(fp,"ZONE T=\"%s\"\n",filename);
      fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TETRAHEDRON\n",buf[0],buf[1]);
    }
    else {
      int dummy;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,MPI_STATUS_IGNORE);
      fp = fopen(filename,"a");
      assert(fp != NULL);
    }

    FOR_INO fprintf(fp,"%18.15le %18.15le %18.15le %18.15le\n",x_no[ino][0],x_no[ino][1],x_no[ino][2],T_no[ino]);

    fclose(fp);

    if ( mpi_rank < mpi_size-1 ) {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }
    MPI_Barrier(mpi_comm);

    // ====================
    // then noote...
    // ====================

    if ( mpi_rank == 0 ) {
      fp = fopen(filename,"a");
      assert(fp != NULL);
    }
    else {
      int dummy;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1235,mpi_comm,MPI_STATUS_IGNORE);
      fp = fopen(filename,"a");
      assert(fp != NULL);
    }

    for (int ite = 0; ite < nte_ib; ++ite) {
      fprintf(fp,"%d %d %d %d\n",
              flag_no[noote[ite][0]]+1,
              flag_no[noote[ite][1]]+1,
              flag_no[noote[ite][2]]+1,
              flag_no[noote[ite][3]]+1);
    }
    delete[] flag_no;

    fclose(fp);

    if ( mpi_rank < mpi_size-1 ) {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1235,mpi_comm);
    }
    MPI_Barrier(mpi_comm);

    delete[] flag_no;

  }

  void writeTecplotByZone(const int step) {

    // handshake version of ascii tecplot write...
    
    int my_nzn = 0;
    FOR_ITE {
      assert(znote[ite] >= 0);
      my_nzn = max(my_nzn,znote[ite]+1);
    }
    int nzn;
    MPI_Allreduce(&my_nzn,&nzn,1,MPI_INT,MPI_MAX,mpi_comm);

    int * flag_no = new int[nno_g];

    for (int izn = 0; izn < nzn; ++izn) {

      char filename[128];
      sprintf(filename,"cht%d.%08d.dat",izn,step);
      COUT1(" > Cht:writeTecplotByZone --> zone[" << izn << "]...");

      int my_buf[2] = { 0 , 0 };
      for (int ino = 0; ino < nno_g ; ++ino) flag_no[ino] = 0; 
      for (int ite = 0; ite < nte_ib; ++ite) {
        if ( izn == ((znote == NULL) ? 0 : znote[ite]) ) {
          my_buf[1] += 1; 
          FOR_I4 { 
            const int ino = noote[ite][i]; 
            if (flag_no[ino] == 0) {
              my_buf[0] += 1;
              flag_no[ino] = 1;
            }
          }
        }
      }
      int buf[2];
      MPI_Reduce(my_buf,buf,2,MPI_INT,MPI_SUM,0,mpi_comm);

      int offset;
      MPI_Scan(&my_buf[0],&offset,1,MPI_INT,MPI_SUM,mpi_comm);
      offset -= my_buf[0];

      // ====================
      // header and nodes first...
      // ====================

      FILE * fp;
      if ( mpi_rank == 0 ) {
        fp = fopen(filename,"w");
        assert(fp != NULL);
        fprintf(fp,"TITLE = \"%s\"\n",filename);
        fprintf(fp,"VARIABLES = \"X\", \"Y\", \"Z\", \"T\"\n");
        fprintf(fp,"ZONE T=\"%s\"\n",filename);
        fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TETRAHEDRON\n",buf[0],buf[1]);
      }
      else {
        int dummy;
        MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,MPI_STATUS_IGNORE);
        fp = fopen(filename,"a");
        assert(fp != NULL);
      }

      int loc = 0; 
      for (int ino = 0; ino < nno_g ; ++ino) flag_no[ino] = 0; 
      for (int ite = 0; ite < nte_ib; ++ite) {
        if ( izn == ((znote == NULL) ? 0 : znote[ite]) ) { 
          FOR_I4 { 
            const int ino = noote[ite][i]; 
            if (flag_no[ino] == 0) {
              fprintf(fp,"%18.15le %18.15le %18.15le %18.15le\n",x_no[ino][0],x_no[ino][1],x_no[ino][2],T_no[ino]);
              flag_no[ino] = ++loc; // ...location in node list (1-indexed)
            }
          }
        }
      }

      fclose(fp);

      if ( mpi_rank < mpi_size-1 ) {
        int dummy = 1;
        MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
      }
      MPI_Barrier(mpi_comm);

      // ====================
      // then noote...
      // ====================
      
      if ( mpi_rank == 0 ) {
        fp = fopen(filename,"a");
        assert(fp != NULL);
      }
      else {
        int dummy;
        MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1235,mpi_comm,MPI_STATUS_IGNORE);
        fp = fopen(filename,"a");
        assert(fp != NULL);
      }

      for (int ite = 0; ite < nte_ib; ++ite) {
        if ( izn == ((znote == NULL) ? 0 : znote[ite]) ) {
          FOR_I4 assert(flag_no[noote[ite][i]] > 0); 
          fprintf(fp,"%d %d %d %d\n",
                  flag_no[noote[ite][0]]+offset,
                  flag_no[noote[ite][1]]+offset,
                  flag_no[noote[ite][2]]+offset,
                  flag_no[noote[ite][3]]+offset);
        }
      }

      fclose(fp);

      if ( mpi_rank < mpi_size-1 ) {
        int dummy = 1;
        MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1235,mpi_comm);
      }
      MPI_Barrier(mpi_comm);

    }
    delete[] flag_no;
  }

  void writeTecplotSerial(const string& filename,const double * const var_no,const double (*const x_no)[3],const int nno,const int (*const noote)[4],const int nte) {

    FILE * fp = fopen(filename.c_str(),"w");
    fprintf(fp,"TITLE = \"%s\"\n",filename.c_str());
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");
    fprintf(fp,"\"VAR\"\n");

    fprintf(fp,"ZONE T=\"%s\"\n",filename.c_str());
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TETRAHEDRON\n",nno,nte);

    FOR_INO {
      fprintf(fp,"%18.15le %18.15le %18.15le %18.15le\n",x_no[ino][0],x_no[ino][1],x_no[ino][2],var_no[ino]);
    }

    for (int ite = 0; ite < nte; ++ite) {
      fprintf(fp,"%d %d %d %d\n",
              noote[ite][0]+1,
              noote[ite][1]+1,
              noote[ite][2]+1,
              noote[ite][3]+1);
    }

    fclose(fp);

  }

  void buildNoono() {

    if (mpi_rank == 0) cout << "Cht::buildNoono()..." << endl;

    assert(noote != NULL);
    assert(noono_i == NULL);
    assert(noono_v == NULL);

    // build teono_i/v...

    int * teono_i = new int[nno+1];
    FOR_INO teono_i[ino+1] = 0;
    FOR_ITE {
      FOR_I4 {
        const int ino = noote[ite][i];
        if (ino < nno)
          ++teono_i[ino+1];
      }
    }
    teono_i[0] = 0;
    FOR_INO teono_i[ino+1] += teono_i[ino];
    const int teono_s = teono_i[nno];
    int * teono_v = new int[teono_s];
    FOR_ITE {
      FOR_I4 {
        const int ino = noote[ite][i];
        if (ino < nno)
          teono_v[teono_i[ino]++] = ite;
      }
    }
    for (int ino = nno-1; ino > 0; --ino)
      teono_i[ino] = teono_i[ino-1];
    teono_i[0] = 0;

    // now use that to build noono_i/v...

    noono_i = new int[nno+1];
    FOR_INO noono_i[ino+1] = 0;

    int * no_flag = new int[nno_g];
    for (int ino = 0; ino < nno_g; ++ino) no_flag[ino] = -1;

    FOR_INO {
      ++noono_i[ino+1]; // diagonal...
      no_flag[ino] = ino;
      for (int ton = teono_i[ino]; ton != teono_i[ino+1]; ++ton) {
        const int ite = teono_v[ton];
        FOR_I4 {
          const int ino_nbr = noote[ite][i];
          if (no_flag[ino_nbr] != ino) {
            ++noono_i[ino+1];
            no_flag[ino_nbr] = ino;
          }
        }
      }
    }
    noono_i[0] = 0;
    FOR_INO noono_i[ino+1] += noono_i[ino];
    const int noono_s = noono_i[nno];
    noono_v = new int[noono_s];

    FOR_INO {
      noono_v[noono_i[ino]++] = ino; // diagonal...
      no_flag[ino] = ino+nno; // nno offset trick
      for (int ton = teono_i[ino]; ton != teono_i[ino+1]; ++ton) {
        const int ite = teono_v[ton];
        FOR_I4 {
          const int ino_nbr = noote[ite][i];
          if (no_flag[ino_nbr] != ino+nno) {
            noono_v[noono_i[ino]++] = ino_nbr;
            no_flag[ino_nbr] = ino+nno;
          }
        }
      }
    }

    for (int ino = nno-1; ino > 0; --ino)
      noono_i[ino] = noono_i[ino-1];
    noono_i[0] = 0;

    delete[] no_flag;
    delete[] teono_i;
    delete[] teono_v;

    // checks...

    // [0:nno_i) nodes only have [0:nno) nbrs...

    for (int ino = 0; ino < nno_i; ++ino) {
      for (int non = noono_i[ino]; non != noono_i[ino+1]; ++non) {
        const int ino_nbr = noono_v[ino];
        assert((ino_nbr >= 0)&&(ino_nbr < nno));
      }
    }

  }

  void buildOperators() {

    if (mpi_rank == 0) cout << "Cht::buildOperators()..." << endl;

    // grad_te and vol_te,vol_no...

    assert(grad_te == NULL);
    assert(vol_te == NULL);
    //assert(vol_no == NULL);
    grad_te = new double[nte][4][3];
    vol_te = new double[nte];
    //vol_no = new double[nno];

    //FOR_INO vol_no[ino] = 0.0;

    FOR_ITE {
      const int ino0 = noote[ite][0];
      const int ino1 = noote[ite][1];
      const int ino2 = noote[ite][2];
      const int ino3 = noote[ite][3];
      // compute the inward normal...
      const double n0[3] = TRI_NORMAL_2(x_no[ino1],x_no[ino3],x_no[ino2]);
      const double n1[3] = TRI_NORMAL_2(x_no[ino0],x_no[ino2],x_no[ino3]);
      const double n2[3] = TRI_NORMAL_2(x_no[ino0],x_no[ino3],x_no[ino1]);
      const double n3[3] = TRI_NORMAL_2(x_no[ino0],x_no[ino1],x_no[ino2]);
      const double vol6 = SIGNED_TET_VOLUME_6(x_no[ino0],x_no[ino1],x_no[ino2],x_no[ino3]);
      assert(vol6 > 0.0);
      FOR_I3 grad_te[ite][0][i] = n0[i]/vol6;
      FOR_I3 grad_te[ite][1][i] = n1[i]/vol6;
      FOR_I3 grad_te[ite][2][i] = n2[i]/vol6;
      FOR_I3 grad_te[ite][3][i] = n3[i]/vol6;
      vol_te[ite] = vol6/6.0;
      // distribute 1/4 to each node...
      //if (ino0 < nno) vol_no[ino0] += vol6/24.0;
      //if (ino1 < nno) vol_no[ino1] += vol6/24.0;
      //if (ino2 < nno) vol_no[ino2] += vol6/24.0;
      //if (ino3 < nno) vol_no[ino3] += vol6/24.0;
    }

    // laplace_noono...

    /*
    assert(laplace_noono == NULL);
    laplace_noono = new double[noono_i[nno]];
    for (int non = 0; non < noono_i[nno]; ++non)
      laplace_noono[non] = 0.0;

    FOR_ITE {
      FOR_I4 {
        const int ino = noote[ite][i];
        if (ino < nno) {
          const int non = noono_i[ino]; assert(noono_v[non] == ino);
          FOR_J4 {
            const int ino_nbr = noote[ite][j];
            int non_nbr = non;
            while (noono_v[non_nbr] != ino_nbr)
              ++non_nbr;
            assert(non_nbr < noono_i[ino+1]); // make sure we found ino_nbr...
            // integration-by-parts flips the sign of this term...
            laplace_noono[non_nbr] -= DOT_PRODUCT(grad_te[ite][i],grad_te[ite][j])*vol_te[ite];
          }
        }
      }
    }
    */

  }

  void calcTeGrad(double (*grad_phi)[3],const double * const phi) {

    FOR_ITE {
      FOR_J3 grad_phi[ite][j] = 0.0;
      FOR_I4 {
        const int ino = noote[ite][i];
        FOR_J3 {
          grad_phi[ite][j] += grad_te[ite][i][j]*phi[ino];
        }
      }
    }

  }

  void checkTeGrad() {

    double *phi = new double[nno_g];
    const double grad_phi_exact[3] = { 1.21321, 3.2124, 1.437829 };
    for (int ino = 0; ino < nno_g; ++ino) phi[ino] = 1.1378 + DOT_PRODUCT(grad_phi_exact,x_no[ino]);

    double (*grad_phi)[3] = new double[nte][3];
    calcTeGrad(grad_phi,phi);

    double ltwo_err = 0.0; // Squared integral norm: sqrt(\int err^2 dV)
    double lmax_err = -1.0; int lmax_err_ite = -1;
    FOR_ITE {
      //cout << "grad: " << COUT_NOC(grad_phi[ite]) << " exact: " << COUT_VEC(grad_phi_exact) << endl;
      //getchar();
      const double err_tmp3[3] = DIFF(grad_phi[ite],grad_phi_exact);
      const double ite_err = MAG(err_tmp3);
      ltwo_err += ite_err*ite_err*vol_te[ite];
      if(ite_err > lmax_err) { lmax_err = ite_err; lmax_err_ite = ite; }
    }
    ltwo_err = sqrt(ltwo_err); lmax_err *= vol_te[lmax_err_ite];
    const int ino0 = noote[lmax_err_ite][0]; const int ino1 = noote[lmax_err_ite][1];
    const int ino2 = noote[lmax_err_ite][2]; const int ino3 = noote[lmax_err_ite][3];
    const double lmax_err_pos[3] = { (x_no[ino0][0]+x_no[ino1][0]+x_no[ino2][0]+x_no[ino3][0])/4.0,
                                     (x_no[ino0][1]+x_no[ino1][1]+x_no[ino2][1]+x_no[ino3][1])/4.0,
                                     (x_no[ino0][2]+x_no[ino1][2]+x_no[ino2][2]+x_no[ino3][2])/4.0 };

    cout << "checkTeGrad: max_err = " << lmax_err << " at ite = " << lmax_err_ite << " near " << COUT_VEC(lmax_err_pos) << endl;
    cout << "checkTeGrad: two_err = " << ltwo_err << endl;

    delete[] phi;
    delete[] grad_phi;

  }

  void checkLaplace() {

    // this routine has not thought through parallel yet...
    cout << "checkLaplace()" << endl;
    assert(laplace_k_noono);
    assert(mpi_size == 0);

    double max_assym = 0.0;
    FOR_INO {
      for (int non = noono_i[ino]+1; non != noono_i[ino+1]; ++non) {
        const int ino_nbr = noono_v[non];
        assert(ino_nbr != ino);
        // look for ino in ino_nbr's row...
        int non_nbr = noono_i[ino_nbr]+1; // no need to consider the diagonal
        while (noono_v[non_nbr] != ino)
          ++non_nbr;
        assert(non_nbr < noono_i[ino_nbr+1]);
        assert(non_nbr != non);
        max_assym = max(max_assym,fabs(laplace_k_noono[non_nbr]-laplace_k_noono[non]));
      }
    }

    cout << " > max_assym: " << max_assym << endl;

    double *phi = new double[nno];
    FOR_INO phi[ino] = (x_no[ino][0]-0.5)*(x_no[ino][0]-0.5);

    double *lap_phi = new double[nno];
    FOR_INO {
      lap_phi[ino] = laplace_k_noono[noono_i[ino]]*phi[ino];
      for (int non = noono_i[ino]+1; non != noono_i[ino+1]; ++non) {
        lap_phi[ino] += laplace_k_noono[non]*phi[noono_v[non]];
      }
      //lap_phi[ino] /= vol_no[ino];
    }

    writeTecplotSerial("lap.dat",lap_phi,x_no,nno,noote,nte);
    // CWH: add a max/l2 check here some day...
    // TODO: need to do flux boundary closure

    delete[] phi;
    delete[] lap_phi;

  }

  // TODO: put this somewhere. It has no class dependecies!
  // TODO: put this somewhere. It has no class dependecies!
  // TODO: put this somewhere. It has no class dependecies!
  // TODO: put this somewhere. It has no class dependecies!
  // TODO: put this somewhere. It has no class dependecies!
  // TODO: put this somewhere. It has no class dependecies!
  // AND USE IT IN PMOVE: it is identical (except rbi_g is
  // vector, so pass ptr to first element)!

  void buildPrcomm(vector<CvPrcomm>& cvPrcommVec,const int ncv,const int ncv_g,uint8 * rbi_g) {

    // this routine uses the names "cv, ncv, icv, cvPrcommVec", but works for nodes, etc. or any other data
    // topology where paired communicators need to be built to update ghost data...

    if (mpi_rank == 0) cout << "buildPrcomm()" << endl;

    // ------------------------------------------------
    // build the paired communicator -- this should be
    // symmetric in terms of rank, but possibly not count...
    // ------------------------------------------------

    assert(cvPrcommVec.empty());

    int rank_current = -1;
    int bits_current = -1;
    CvPrcomm * prcomm = NULL;
    for (int icv = ncv;  icv < ncv_g; ++icv) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv-ncv]);
      if (rank > rank_current) {
        // we just switched ranks. If there was a prcomm, then
        // complete the size...
        if (prcomm)
          prcomm->unpack_size = icv-prcomm->unpack_offset;
        rank_current = rank;
        bits_current = -1;
        cvPrcommVec.push_back(CvPrcomm());
        prcomm = &cvPrcommVec.back();
        prcomm->rank = rank;
        prcomm->unpack_offset = icv;
      }
      else {
        assert(rank_current == rank);
        assert(prcomm);
      }
      // just tests monotoncity of the bits...
      if (bits > bits_current) {
        bits_current = bits;
      }
      else {
        assert(bits_current == bits);
      }
    }
    // we just finished. If there was a prcomm, then
    // complete the size...
    if (prcomm) prcomm->unpack_size = ncv_g-prcomm->unpack_offset;

    /*
      MPI_Pause("1");
      FOR_RANK {
      if (mpi_rank == rank) {
      cout << "rank " << rank << " communicates with ";
      for (int ii = 0, ii_end=cvPrcommVec.size(); ii < ii_end; ++ii)
      cout << " " << cvPrcommVec[ii].rank;
      cout << endl;
      }
      MPI_Pause("ok");
      }
    */

    // finally, we need to send/recv the indices and bits to the pack side
    // and build the packVecs. Note that these are not necessarily symmetric by construction
    // because internal faces may have been removed based on tolerance on one vd,
    // and not on the other...

    /*
      int * send_count = new int[mpi_size];
      FOR_RANK send_count[rank] = 0;
      for (int ii = 0, ii_end=cvPrcommVec.size(); ii < ii_end; ++ii) {
      assert(cvPrcommVec[ii].unpack_size > 0);
      send_count[cvPrcommVec[ii].rank] = cvPrcommVec[ii].unpack_size;
      }
      int * recv_count = new int[mpi_size];
      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
      for (int ii = 0, ii_end=cvPrcommVec.size(); ii < ii_end; ++ii) {
      cvPrcommVec[ii].pack_size = recv_count[cvPrcommVec[ii].rank];
      recv_count[cvPrcommVec[ii].rank] = 0; // zero to add any non-zeros below...
      }
      // add any that are missing...
      FOR_RANK if (recv_count[rank] != 0) {
      cvPrcommVec.push_back(CvPrcomm());
      prcomm = &cvPrcommVec.back();
      prcomm->rank = rank;
      prcomm->unpack_size = 0;
      prcomm->unpack_offset = ncv;
      prcomm->pack_size = recv_count[rank];
      }
      delete[] send_count;
      delete[] recv_count;
    */

    MPI_Request * sendRequestArray = new MPI_Request[cvPrcommVec.size()];
    MPI_Request * recvRequestArray = new MPI_Request[cvPrcommVec.size()];

    for (int ii = 0, ii_end=cvPrcommVec.size(); ii < ii_end; ++ii) {

      // post irecv...
      MPI_Irecv(&(cvPrcommVec[ii].pack_size),1,MPI_INT,cvPrcommVec[ii].rank,12345,mpi_comm,&(recvRequestArray[ii]));

      // and the send...
      MPI_Issend(&(cvPrcommVec[ii].unpack_size),1,MPI_INT,cvPrcommVec[ii].rank,12345,mpi_comm,&(sendRequestArray[ii]));

    }

    MPI_Waitall(cvPrcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
    MPI_Waitall(cvPrcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

    // now send from the unpack side to the pack side...

    int pack_size = 0;
    for (int ii = 0, ii_end=cvPrcommVec.size(); ii < ii_end; ++ii)
      pack_size += cvPrcommVec[ii].pack_size;

    uint8 * recv_rbi_g = new uint8[pack_size];

    pack_size = 0;
    for (int ii = 0, ii_end=cvPrcommVec.size(); ii < ii_end; ++ii) {

      // post irecv...
      MPI_Irecv(recv_rbi_g+pack_size,cvPrcommVec[ii].pack_size,MPI_UINT8,
                cvPrcommVec[ii].rank,12346,mpi_comm,&(recvRequestArray[ii]));
      pack_size += cvPrcommVec[ii].pack_size;

      // and the send...
      MPI_Issend(rbi_g+(cvPrcommVec[ii].unpack_offset-ncv),cvPrcommVec[ii].unpack_size,MPI_UINT8,
                 cvPrcommVec[ii].rank,12346,mpi_comm,&(sendRequestArray[ii]));

    }

    MPI_Waitall(cvPrcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
    MPI_Waitall(cvPrcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

    delete[] sendRequestArray;
    delete[] recvRequestArray;

    // now build the packVec and periodicity...

    double R[9], t[3];
    pack_size = 0;
    for (int ii = 0, ii_end=cvPrcommVec.size(); ii < ii_end; ++ii) {
      assert(cvPrcommVec[ii].packVec.empty());
      cvPrcommVec[ii].packVec.resize(cvPrcommVec[ii].pack_size);
      CvPrcomm::Transform * transform = NULL;
      int bits_current = -1;
      for (int i = 0; i < cvPrcommVec[ii].pack_size; ++i) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,recv_rbi_g[pack_size+i]);
        assert(rank == mpi_rank);
        if (bits > bits_current) {
          // bits are about to change, so complete any translate and rotate...
          if (transform) transform->end = i;
          // look for new rotations/translations...
          const bool has_R = false; // getPeriodicR(R,bits); <- TODO
          const bool has_t = false; // getPeriodicT(t,bits);
          cvPrcommVec[ii].transformVec.push_back(CvPrcomm::Transform(has_R,R,has_t,t,bits,i));
          transform = &(cvPrcommVec[ii].transformVec.back());
          bits_current = bits;
        }
        else {
          assert(bits_current == bits);
        }
        assert((index >= 0)&&(index < ncv));
        cvPrcommVec[ii].packVec[i] = index;
      }
      if (transform) transform->end = cvPrcommVec[ii].pack_size;
      pack_size += cvPrcommVec[ii].pack_size;
    }

    delete[] recv_rbi_g;

  }

  void reorderPrcommPackVec(vector<CvPrcomm>& cvPrcommVec,const int * const order,const int ncv) {

    for (int ii = 0; ii < cvPrcommVec.size(); ++ii) {
      for (int i = 0, limit = cvPrcommVec[ii].packVec.size(); i < limit; ++i) {
        const int icv = cvPrcommVec[ii].packVec[i];
        assert((icv >= 0)&&(icv < ncv));
        cvPrcommVec[ii].packVec[i] = order[icv];
      }
    }

  }

  // TODO: put this in a solver area eventually...

  int solveNoCg(double * phi,const double * const A,const double * const rhs,
                const double zero,const int maxiter,const bool verbose) {

    // we need the following work arrays...
    double * res      = new double[nno];
    double * v        = new double[nno];
    double * p        = new double[nno_g];
    double * inv_diag = new double[nno];

    // initialize...
    for (int ino = 0; ino < nno; ++ino)  {
      assert(noono_v[noono_i[ino]] == ino); // confirm diagonal first in noono_i/v CSR
      assert(A[noono_i[ino]] != 0.0); // diag cannot be zero for diag-preconditioned cg.
      inv_diag[ino] = 1.0/A[noono_i[ino]];
    }

    for (int ino = 0; ino < nno_g; ++ino)
      p[ino] = 0.0;
    double rho = 1.0;

    // calculate the residual in rhs format...
    for (int ino = 0; ino < nno; ++ino) {
      res[ino] = rhs[ino] - A[noono_i[ino]]*phi[ino];
      for (int non = noono_i[ino]+1; non < noono_i[ino+1]; ++non) {
        const int ino_nbr = noono_v[non];
        res[ino] -= A[non]*phi[ino_nbr];
      }
    }

    // diagonal precon/compute normalized residual...
    double my_res_max = 0.0;
    for (int ino = 0; ino < nno; ++ino) {
      v[ino] = res[ino]*inv_diag[ino];
      my_res_max = max( my_res_max, fabs(v[ino]) );
    }
    double res_max = 0.0;
    if (verbose) {
      MPI_Reduce(&my_res_max,&res_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      if ((mpi_rank == 0)&&verbose) cout << " > initial res_max: " << res_max << endl;
    }

    int iter = 0;
    int done = 0;
    while (done == 0) {

      ++iter;

      assert(rho != 0.0);
      const double rho_prev = rho;
      //if (fabs(rho_prev) < 1.0E-20)
      //rho_prev = -1.0E-20; // -1.0E-20? seems to help

      double my_rho = 0.0;
      for (int ino = 0; ino < nno; ++ino)
        my_rho += res[ino]*v[ino];
      MPI_Allreduce(&my_rho,&rho,1,MPI_DOUBLE,MPI_SUM,mpi_comm);

      assert(rho_prev != 0.0);
      const double beta = rho/rho_prev;
      for (int ino = 0; ino < nno; ++ino)
        p[ino] = v[ino] + beta*p[ino];

      // v = [Ap]{p}...
      // with some attempt to hide latency...

      updateNoDataStart(p);

      for (int ino = 0; ino < nno_i; ++ino) {
        v[ino] = A[noono_i[ino]]*p[ino];
        for (int non = noono_i[ino]+1; non < noono_i[ino+1]; ++non) {
          const int ino_nbr = noono_v[non];
          assert(ino_nbr < nno);
          v[ino] += A[non]*p[ino_nbr];
        }
      }

      updateNoDataFinish(p);

      for (int ino = nno_i; ino < nno; ++ino) {
        v[ino] = A[noono_i[ino]]*p[ino];
        for (int non = noono_i[ino]+1; non < noono_i[ino+1]; ++non) {
          const int ino_nbr = noono_v[non];
          v[ino] += A[non]*p[ino_nbr];
        }
      }

      double my_gamma = 0.0;
      for (int ino = 0; ino < nno; ++ino)
        my_gamma += p[ino]*v[ino];
      double gamma;
      MPI_Allreduce(&my_gamma,&gamma,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
      assert(gamma != 0.0);
      assert(gamma == gamma); // nan check

      const double alpha = rho/gamma;

      // check if we are done...
      if (iter%3 == 0) {

        for (int ino = 0; ino < nno_g; ++ino) {
          assert(p[ino] == p[ino]); // nan check
          phi[ino] += alpha*p[ino];
        }

        // recompute the residual...
        for (int ino = 0; ino < nno; ++ino) {
          res[ino] = rhs[ino] - A[noono_i[ino]]*phi[ino];
          for (int non = noono_i[ino]+1; non < noono_i[ino+1]; ++non) {
            const int ino_nbr = noono_v[non];
            res[ino] -= A[non]*phi[ino_nbr];
          }
        }

        for (int ino = 0; ino < nno; ++ino)
          v[ino] = res[ino]*inv_diag[ino];

        // compute the max (L-infinity) normalized residual...
	my_res_max = 0.0;
        for (int ino = 0; ino < nno; ++ino)
          my_res_max = max( my_res_max, fabs(v[ino]) );
        MPI_Reduce(&my_res_max,&res_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
        if (mpi_rank == 0) {
          // only share the last half of the convergence behaviour...
          if (verbose || (iter > maxiter/2))
            cout << " > solveCgLocal iter " << iter << " res_max " << res_max << endl;
          if (res_max <= zero) {
            if (verbose) cout << "-> Successfully converged error to " << res_max << endl;
            done = 1;
          }
          else if (iter > maxiter) {
            cout << "Warning: solveCgLocal did not converge after " << maxiter <<
              " iters, res_max: " << res_max << endl;
            done = 2;
          }
        }
        MPI_Bcast(&done,1,MPI_INT,0,mpi_comm);

      }
      else {

        // update full phi...
        for (int ino = 0; ino < nno_g; ++ino)
          phi[ino] += alpha*p[ino];

        for (int ino = 0; ino < nno; ++ino) {
          // on the other iterations, use this approximation to update
          // the unreduced residual...
          res[ino] -= alpha*v[ino];
          // still need to compute v, diag precon for next iteration...
          v[ino] = res[ino]*inv_diag[ino];
        }

      }

      //MiscUtils::dumpRange(phi,nno,"PHI RANGE");
      //getchar();

    }

    delete[] res;
    delete[] v;
    delete[] p;
    delete[] inv_diag;

    // let the calling routine know if we were successful...
    return done;

  }

  void read_tles(const string& filename) {

    assert(0);

    /*

    // TODO some of this should be persistent, but for now lets organize the data just like FluentMsh
    int8 nte_global = 0; // global number of tets
    //int8 nno_global = 0; // global number of nodes
    int8 nst_global = 0; // global number of surface-tris
    int8 (*noote_global)[4] = NULL; // node-of-tet global index
    //int8 *teost_global = NULL;      // tet-of-surface-tri global index (ist ordered)
    int8 (*noost_global)[3] = NULL; // node-of-surface-tri global index (spost ordered)
    //int8 *noora = NULL; // node-of-rank
    int8 *teora = NULL; // tet-of-rank
    int8 *stora = NULL; // surface-tri-of-rank

    COUT1("read_tles(): " << filename);

    char dummy[128]; assert(filename.length() < 128);
    sprintf(dummy,"%s",filename.c_str());
    MPI_File fh;
    int ierr = MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
    if (ierr != 0) {
      CERR("read_tles: file not found: \"" << filename << "\"");
    }

    // recall file starts: { TETS_IO_MAGIC_NUMBER, TETS_IO_VERSION, nte_global, nno_global, nst_global };
    int8 ibuf[5];
    MPI_Offset offset = 0;
    if (mpi_rank == 0)
      MPI_File_read_at(fh,offset,ibuf,5,MPI_INT8,MPI_STATUS_IGNORE);
    offset += 5*int8_size;
    MPI_Bcast(ibuf,5,MPI_INT8,0,mpi_comm);
    assert(ibuf[0] == TETS_IO_MAGIC_NUMBER); // implement byteswapping if you hit this
    assert(ibuf[1] == TETS_IO_VERSION);
    assert(nte_global == 0);
    nte_global = ibuf[2];
    assert(nno_global == 0);
    nno_global = ibuf[3];
    assert(nst_global == 0);
    nst_global = ibuf[4];

    if (mpi_rank == 0)
      cout << " > nte_global: " << nte_global << ", nno_global: " << nno_global << ", nst_global: " << nst_global << endl;

    // xora's...

    assert(noora_striped == NULL);
    MiscUtils::calcThresholdDist(noora_striped,nno_global,mpi_size,DIST_THRESHOLD);
    assert(noora_striped[mpi_rank+1]-noora_striped[mpi_rank] < TWO_BILLION);
    assert(nno == 0);
    nno = int(noora_striped[mpi_rank+1]-noora_striped[mpi_rank]);

    assert(teora == NULL);
    MiscUtils::calcThresholdDist(teora,nte_global,mpi_size,DIST_THRESHOLD);
    assert(teora[mpi_rank+1]-teora[mpi_rank] < TWO_BILLION);
    assert(nte == 0);
    nte = int(teora[mpi_rank+1]-teora[mpi_rank]);

    assert(stora == NULL);
    MiscUtils::calcThresholdDist(stora,nst_global,mpi_size,DIST_THRESHOLD);
    assert(stora[mpi_rank+1]-stora[mpi_rank] < TWO_BILLION);
    assert(nst == 0);
    nst = int(stora[mpi_rank+1]-stora[mpi_rank]);

    // x_no...

    assert(x_no == NULL); x_no = new double[nno][3];
    MPI_File_read_at_all(fh,offset+noora_striped[mpi_rank]*3*double_size,(double*)x_no,nno*3,MPI_DOUBLE,MPI_STATUS_IGNORE);
    offset += nno_global*3*double_size;

    // noote_global...

    assert(noote_global == NULL); noote_global = new int8[nte][4];
    MPI_File_read_at_all(fh,offset+teora[mpi_rank]*4*int8_size,(int8*)noote_global,nte*4,MPI_INT8,MPI_STATUS_IGNORE);
    offset += nte_global*4*int8_size;

    // teost_global...

    //assert(teost_global == NULL); teost_global = new int8[nst];
    //MPI_File_read_at_all(fh,offset+stora[mpi_rank]*int8_size,teost_global,nst,MPI_INT8,MPI_STATUS_IGNORE);
    offset += nst_global*int8_size;

    // noost_global...

    assert(noost_global == NULL); noost_global = new int8[nst][3];
    MPI_File_read_at_all(fh,offset+stora[mpi_rank]*3*int8_size,(int8*)noost_global,nst*3,MPI_INT8,MPI_STATUS_IGNORE);
    offset += nst_global*3*int8_size;

    // TODO read in znost?

    // and close...

    MPI_File_close(&fh);

    // copying into integer representations (should that be the default?)

    assert(nno_global < TWO_BILLION);
    assert(noote == NULL); noote = new int[nte][4];
    for (int ite = 0; ite < nte; ++ite)
      FOR_I4 noote[ite][i] = (int)noote_global[ite][i];
    assert(noost == NULL); noost = new int[nst][3];
    for (int ist = 0; ist < nst; ++ist)
      FOR_I3 noost[ist][i] = (int)noost_global[ist][i];

    // cleanup...

    DELETE(noote_global);
    //DELETE(teost_global);
    DELETE(noost_global);
    //DELETE(noora);
    DELETE(teora);
    DELETE(stora);

    */

  }

  void readData(const string& filename) {

    assert(ino_global);
    assert(dde_striped);
    assert(noora_striped);

    COUT1("readData(): " << filename);

    char dummy[128]; assert(filename.length() < 128);
    sprintf(dummy,"%s",filename.c_str());
    MPI_File fh;
    int ierr = MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
    if (ierr != 0) {
      CERR("readData: file not found: \"" << filename << "\"");
    }

    int8 nno_global_check;
    MPI_Offset offset = 0;
    if (mpi_rank == 0)
      MPI_File_read_at(fh,offset,&nno_global_check,1,MPI_INT8,MPI_STATUS_IGNORE);
    offset += int8_size;
    MPI_Bcast(&nno_global_check,1,MPI_INT8,0,mpi_comm);
    assert(nno_global_check == nno_global);

    // T_no...

    int nno_striped = int(noora_striped[mpi_rank+1]-noora_striped[mpi_rank]);
    double *dbuf_striped = new double[nno_striped];
    MPI_File_read_at_all(fh,offset+noora_striped[mpi_rank]*double_size,dbuf_striped,nno_striped,MPI_DOUBLE,MPI_STATUS_IGNORE);
    offset += nno_global*double_size;
    dde_striped->pull(T_no,dbuf_striped);
    delete[] dbuf_striped;
    updateNoData(T_no);

    // and close...

    MPI_File_close(&fh);

  }

  void writeData(const string& filename) {

    assert(ino_global);
    assert(dde_striped);
    assert(noora_striped);

    COUT1("writeData(): " << filename);

    char fname[128];
    sprintf(fname,"%s",filename.c_str());
    MPI_File_delete(fname,MPI_INFO_NULL);

    MPI_File fh;
    MPI_File_open(mpi_comm,fname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    if (mpi_rank == 0)
      MPI_File_write(fh,&nno_global,1,MPI_INT8,MPI_STATUS_IGNORE);
    MPI_Offset offset = int8_size;

    // T_no...

    int nno_striped = int(noora_striped[mpi_rank+1]-noora_striped[mpi_rank]);
    double *dbuf_striped = new double[nno_striped];
    dde_striped->push(dbuf_striped,T_no);
    MPI_File_write_at_all(fh,offset+noora_striped[mpi_rank]*double_size,dbuf_striped,nno_striped,MPI_DOUBLE,MPI_STATUS_IGNORE);
    delete[] dbuf_striped;
    offset += nno_global*double_size;

    // and close...

    MPI_File_set_size(fh,offset);
    MPI_File_close(&fh);

  }

};

#endif
