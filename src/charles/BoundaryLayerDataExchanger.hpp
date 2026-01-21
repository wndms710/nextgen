#ifndef BOUNDARY_LAYER_DATA_EXCHANGER_HPP
#define BOUNDARY_LAYER_DATA_EXCHANGER_HPP

template<typename SolverT>
class BoundaryLayerDataExchanger {


public:
  //data associated with local ibf's
  double *r1_bl;  //value of some r1 data at these points...to be computed 

private:
  
  //data associated with local ibf's
  
  BfZone * zone_ptr;
  SolverT* solver;
  int nbl;            //number of points in boundardy layer (per bf)
  double l_bl;            //distance bl mesh extends from bf
  double (*x_bl)[3];   //list of points of size nbl*zone_ptr->nbf  DAP TODO: maybe don't need to store this, location can always be computed

  //data assciated with local icv's
  int npt;             //number of points matched to cv's on my rank
  double (*x_pt)[3];   //the points matched to cv's on my rank
  int *icv_pt;         //indices of the cv's matched to x_pt on my rank
  // structures for sending pt buffers back to bl buffers
  int *send_count_r1;
  int *send_disp_r1;
  int *recv_count_r1;
  int *recv_disp_r1;

  int *blopt; 
  int blopt_s;

  //DAP TODO: consider interpolating to point location rather than just returning icv value...
  //int *wtopt_i; // for interpolating from noocv to pt
  //pair<int,double> *wtopt_v; 
public:

  BoundaryLayerDataExchanger() {
    zone_ptr = NULL;
    solver   = NULL;
    nbl = 0;
    l_bl = 0.0;
    x_bl = NULL;
    r1_bl = NULL;
     
    npt = 0;
    x_pt = NULL;
    icv_pt = NULL;

    send_count_r1 = NULL;
    send_disp_r1  = NULL;
    recv_count_r1 = NULL;
    recv_disp_r1  = NULL;
  
    blopt = NULL;
    blopt_s = 0;

    //wtopt_i = NULL;
    //wtopt_v = NULL;
  }

  BoundaryLayerDataExchanger(BfZone* _zone_ptr, SolverT* _solver, const int _nbl, const double _l_bl) {
    zone_ptr = _zone_ptr;
    solver = _solver;
    nbl = _nbl;
    l_bl = _l_bl;
    x_bl = NULL;
    r1_bl = NULL;
     
    npt = 0;
    x_pt = NULL;
    icv_pt = NULL;

    send_count_r1 = NULL;
    send_disp_r1  = NULL;
    recv_count_r1 = NULL;
    recv_disp_r1  = NULL;

    blopt = NULL;
    blopt_s = 0;
   
    //wtopt_i = NULL;
    //wtopt_v = NULL;
    //
    
    init();
  }

  ~BoundaryLayerDataExchanger() {
    DELETE(x_bl);
    DELETE(r1_bl);
    DELETE(x_pt);
    DELETE(icv_pt);
    DELETE(send_count_r1);
    DELETE(send_disp_r1);
    DELETE(recv_count_r1);
    DELETE(recv_disp_r1);
    DELETE(x_pt);
    DELETE(icv_pt);
    DELETE(blopt);
    //DELETE(wtopt_i);
    //DELETE(wtopt_v);
  }


  void init(){
    assert(zone_ptr!=NULL&&nbl>0&&l_bl>0.0);
    assert(x_bl==NULL);

    //build bl point grid
    x_bl = new double[zone_ptr->nbf * nbl][3];
    double d_bl = l_bl/nbl;
    for (int ibf=0; ibf<zone_ptr->nbf; ++ibf){
      double mag_n = MAG(zone_ptr->n_bf[ibf]);
      double unit_n[3];
      FOR_I3 unit_n[i] = zone_ptr->n_bf[ibf][i]/mag_n;
      for (int ibl=0; ibl<nbl; ++ibl){
        FOR_I3 x_bl[ibf*nbl + ibl][i] = zone_ptr->x_bf[ibf][i] - unit_n[i]*(ibl+0.5)*d_bl;
      }
    }

    // ------------------------------------------------------------------------
    // now send the bl points to the ranks owning the cv's that contain them
    // ------------------------------------------------------------------------
    
    // build the bbox...

    // make sure we have the stuff we need...

    if (solver->cvAdt == NULL) solver->buildCvAdt();

    assert(send_count_r1==NULL);
    send_count_r1 = new int[mpi_size];
    FOR_RANK send_count_r1[rank] = 0;

    int * send_count_r2 = new int[mpi_size];
    FOR_RANK send_count_r2[rank] = 0;

    vector<int> intVec;
    for (int ibf=0; ibf<zone_ptr->nbf; ++ibf){
      for (int ibl = 0; ibl < nbl; ++ibl) {
        assert(intVec.empty());
        solver->cvBboxAdt->buildListForPoint(intVec,&(x_bl[ibf*nbl+ibl][0]));
        for (int ii = 0; ii < intVec.size(); ++ii) {
          const int rank = intVec[ii]; assert((rank >= 0)&&(rank < mpi_size));
          ++send_count_r1[rank];
          send_count_r2[rank] += 3;
        }
        intVec.clear();
      }
    }
    assert(send_disp_r1==NULL);
    send_disp_r1 = new int[mpi_size];
    int * send_disp_r2 = new int[mpi_size];
    send_disp_r1[0] = 0;
    send_disp_r2[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank){
      send_disp_r1[rank] = send_count_r1[rank-1] + send_disp_r1[rank-1];
      send_disp_r2[rank] = send_count_r2[rank-1] + send_disp_r2[rank-1];
    }
    const int send_count_r1_sum = send_disp_r1[mpi_size-1] + send_count_r1[mpi_size-1];
    const int send_count_r2_sum = send_disp_r2[mpi_size-1] + send_count_r2[mpi_size-1];

    double * send_buf_r2 = new double[send_count_r2_sum];
    int * send_buf_r1 = new int[send_count_r1_sum];
    for (int ibf=0; ibf<zone_ptr->nbf; ++ibf){
      for (int ibl = 0; ibl < nbl; ++ibl) {
        assert(intVec.empty());
        solver->cvBboxAdt->buildListForPoint(intVec,&(x_bl[ibf*nbl+ibl][0]));
        for (int ii = 0; ii < intVec.size(); ++ii) {
          const int rank = intVec[ii]; assert((rank >= 0)&&(rank < mpi_size));
          send_buf_r2[send_disp_r2[rank]  ] = x_bl[ibf*nbl+ibl][0];
          send_buf_r2[send_disp_r2[rank]+1] = x_bl[ibf*nbl+ibl][1];
          send_buf_r2[send_disp_r2[rank]+2] = x_bl[ibf*nbl+ibl][2];
          send_buf_r1[send_disp_r1[rank]] = ibf*nbl+ibl;

          ++send_disp_r1[rank];
          send_disp_r2[rank] += 3;
        }
        intVec.clear();
      }
    }
    delete[] x_bl; x_bl = NULL;  //DAP TODO consider making this a local pointer...

    // rewind...

    send_disp_r1[0] = 0;
    send_disp_r2[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank){
      send_disp_r1[rank] = send_count_r1[rank-1] + send_disp_r1[rank-1];
      send_disp_r2[rank] = send_count_r2[rank-1] + send_disp_r2[rank-1];
    }

    // now send...

    assert(recv_count_r1==NULL);
    recv_count_r1 = new int[mpi_size];
    MPI_Alltoall(send_count_r1,1,MPI_INT,recv_count_r1,1,MPI_INT,mpi_comm);

    int * recv_count_r2 = new int[mpi_size];
    MPI_Alltoall(send_count_r2,1,MPI_INT,recv_count_r2,1,MPI_INT,mpi_comm);

    assert(recv_disp_r1==NULL);
    recv_disp_r1 = new int[mpi_size];
    recv_disp_r1[0] = 0;
    int * recv_disp_r2 = new int[mpi_size];
    recv_disp_r2[0] = 0;

    for (int rank = 1; rank < mpi_size; ++rank){
      recv_disp_r1[rank] = recv_count_r1[rank-1] + recv_disp_r1[rank-1];
      recv_disp_r2[rank] = recv_count_r2[rank-1] + recv_disp_r2[rank-1];
    }
    const int recv_count_r1_sum = recv_disp_r1[mpi_size-1] + recv_count_r1[mpi_size-1];
    const int recv_count_r2_sum = recv_disp_r2[mpi_size-1] + recv_count_r2[mpi_size-1];

    double * recv_buf_r2 = new double[recv_count_r2_sum];
    int * recv_buf_r1 = new int[recv_count_r1_sum];

    MPI_Alltoallv(send_buf_r1,send_count_r1,send_disp_r1,MPI_INT,
                  recv_buf_r1,recv_count_r1,recv_disp_r1,MPI_INT,mpi_comm);
    delete[] send_buf_r1;
    //delete[] send_disp_r1; send_disp_r1 = NULL;
    //delete[] send_count_r1; send_count_r1 = NULL;
   
    MPI_Alltoallv(send_buf_r2,send_count_r2,send_disp_r2,MPI_DOUBLE,
                  recv_buf_r2,recv_count_r2,recv_disp_r2,MPI_DOUBLE,mpi_comm);
    delete[] send_buf_r2;
    delete[] send_disp_r2; 
    delete[] send_count_r2;
   

    // now unpack and set...
    
    assert(recv_count_r2_sum%3 == 0);
    int * flag = new int[recv_count_r2_sum/3]; // -1 
    for (int ipt = 0; ipt < recv_count_r2_sum/3; ++ipt) flag[ipt] = -1; // when flag[ipt] >= 0, x_pt[ipt] lies on this rank

    assert(npt==0);
    for (int irecv = 0; irecv < recv_count_r2_sum; irecv += 3) {
      const int ipt = irecv/3;
      double xp[3]; FOR_I3 xp[i] = recv_buf_r2[irecv+i];
      assert(intVec.empty());
      solver->cvAdt->buildListForPoint(intVec,xp);              
      double dist_sq = HUGE_VAL;
      for (int ii = 0; ii < intVec.size(); ++ii) {
        const int this_icv = intVec[ii]; assert((this_icv >= 0)&&(this_icv < solver->ncv));
        const double this_dist_sq = DIST2(xp,solver->x_vv[this_icv]);
        if (this_dist_sq <= dist_sq && this_dist_sq <= (1.0+1.0E-12)*solver->r_vv[this_icv]*solver->r_vv[this_icv]) {
          flag[ipt] = this_icv;
          dist_sq = this_dist_sq;
        }
      }
      if (flag[ipt] >= 0) {
        const int icv = flag[ipt];
        // need to check against ghosts as well...
        if (icv > solver->ncv_i) {
          for (int coc = solver->cvocv_i[icv]; coc != solver->cvocv_i[icv+1]; ++coc) {
            const int icv_nbr = solver->cvocv_v[coc];
            if (icv_nbr >= solver->ncv) {
              const double this_dist_sq = DIST2(xp,solver->x_vv[icv_nbr]); 
              if (this_dist_sq < dist_sq) {
                flag[ipt] = -1;
              }
              else if (this_dist_sq == dist_sq) {
                // lower rbi owns it...
                int rank,bits,index;
                BitUtils::unpackRankBitsIndex(rank,bits,index,solver->rbi_g[icv_nbr-solver->ncv]);
                if (mpi_rank > rank) {
                  flag[ipt] = -1;
                }
              }
            }
          }
        }
        // if we still own the point...
        if (flag[ipt] >= 0) 
          ++npt;
      }
      intVec.clear();
    }
    assert(x_pt==NULL&&icv_pt==NULL);
    x_pt = new double[npt][3];
    icv_pt = new int[npt];
    int * ibl_pt = new int[npt];
    int ipt = 0;
    FOR_RANK send_count_r1[rank] = 0;
    int theRank = 0;
   
    for (int irecv = 0; irecv < recv_count_r2_sum; irecv += 3) {
      if (irecv - recv_disp_r2[theRank] >= recv_count_r2[theRank]) {
        ++theRank;
        while(recv_count_r2[theRank]==0&&theRank<mpi_size-1)
          ++theRank;
      }
      const int ind = irecv/3;
      if (flag[ind] >= 0) {
        FOR_I3 x_pt[ipt][i] = recv_buf_r2[ind*3+i];
        icv_pt[ipt] = flag[ind];
        ibl_pt[ipt] = recv_buf_r1[ind];
        ++send_count_r1[theRank];
        ++ipt;
      }
    }
    assert(ipt <= recv_count_r2_sum/3);
    delete[] flag;
    delete[] recv_count_r2;
    delete[] recv_disp_r2;
    delete[] recv_buf_r2;
    delete[] recv_buf_r1; 

    //Prepare for data to be returned to original ranks, 
    //build send_count_r1/send_disp_r1/recv_count_r1/recv_disp_r1 for a buffer of size icv_pt
    send_disp_r1[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp_r1[rank] = send_count_r1[rank-1] + send_disp_r1[rank-1];

    MPI_Alltoall(send_count_r1,1,MPI_INT,recv_count_r1,1,MPI_INT,mpi_comm);

    recv_disp_r1[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp_r1[rank] = recv_count_r1[rank-1] + recv_disp_r1[rank-1];
    assert(blopt_s==0);
    blopt_s = recv_disp_r1[mpi_size-1] + recv_count_r1[mpi_size-1];
    //send ibl_pt indices back...
    assert(blopt==NULL);
    blopt = new int[blopt_s];

    MPI_Alltoallv(ibl_pt,send_count_r1,send_disp_r1,MPI_INT,
                  blopt ,recv_count_r1,recv_disp_r1,MPI_INT,mpi_comm);
    delete[] ibl_pt;

    //WOOHOOO! We're finally ready, we have...
    //n_pt, x_pt, icv_pt to look up local data at bl points
    //send_count,send_disp,recv_count,recv_disp,blopt_s to exhange this additional r1 data
    //blopt[blopt_s] to assign exchanged data to r1_bl; NOTE! not all entries in r1_bl will necessarily be found...
    
    //check how many points were not found...
    //TODO often duplicates are found and more points are returned than requested
    //     consider an additional distance exchange and nearest distance check on the original sending ranks?
    int my_buf[2] = {zone_ptr->nbf*nbl,blopt_s};
    int buf[2];
    MPI_Reduce(my_buf,buf,2,MPI_INT,MPI_SUM,0,mpi_comm);
    if (mpi_rank==0){
      cout << " > BoundaryLayerPointExchanger::init() BfZone " << zone_ptr->getName() << " found " << buf[1] << " of " << buf[0] << " requested points" << endl;
    }
  }

  void computeBlFromPt(double * bl_delta, const double &pTotal_threshold){

    if (r1_bl==NULL){
      r1_bl = new double[nbl*zone_ptr->nbf];
    }
    for (int ibf=0; ibf<zone_ptr->nbf; ++ibf){
      for (int ibl=0; ibl<nbl; ++ibl){
        r1_bl[ibf*nbl+ibl] = -HUGE_VAL;
      }
    }

    //TODO check if there's a better way to get these...
    CtiRegister::CtiData * u_avg   = CtiRegister::getCtiData("u_avg");
    CtiRegister::CtiData * p_avg   = CtiRegister::getCtiData("p_avg");
    CtiRegister::CtiData * rho_avg = CtiRegister::getCtiData("rho_avg");
    
    double * send_buf_r1 = new double[npt];
    double * recv_buf_r1 = new double[blopt_s];
    if (u_avg&&p_avg&&rho_avg){
      if (mpi_rank==0) cout << " > BoundaryLayerPointExchanger::computeBlFromPt() BfZone " << zone_ptr->getName() << " using stat data" << endl;
      for (int ipt=0; ipt<npt; ++ipt){
        double u[3] = {u_avg->dn3(icv_pt[ipt],0),u_avg->dn3(icv_pt[ipt],1),u_avg->dn3(icv_pt[ipt],2)};
        double p = p_avg->dn(icv_pt[ipt]);
        double rho = rho_avg->dn(icv_pt[ipt]);
        double u_mag2 = DOT_PRODUCT(u,u);
        send_buf_r1[ipt] = p + 0.5*rho*u_mag2;
      }
    }
    else{
      if (mpi_rank==0) cout << " > BoundaryLayerPointExchanger::computeBlFromPt() BfZone " << zone_ptr->getName() << " using inst data" << endl;
      for (int ipt=0; ipt<npt; ++ipt){
        double u_mag2 = DOT_PRODUCT(solver->u[icv_pt[ipt]],solver->u[icv_pt[ipt]]);
        send_buf_r1[ipt] = solver->p[icv_pt[ipt]] + 0.5*solver->rho[icv_pt[ipt]]*u_mag2;
      }
    }
    assert(send_count_r1&&send_disp_r1&&recv_count_r1&&recv_disp_r1);
    MPI_Alltoallv(send_buf_r1,send_count_r1,send_disp_r1,MPI_DOUBLE,
                  recv_buf_r1,recv_count_r1,recv_disp_r1,MPI_DOUBLE,mpi_comm);
    delete[] send_buf_r1;
   
    //DAPDAP
    //FOR_RANK{
    //  if (mpi_rank == rank){
    //     int max_b = -1;
    //     for (int blopt_i = 0; blopt_i<blopt_s; ++blopt_i){
    //       if (max_b < blopt[blopt_i])
    //         max_b = blopt[blopt_i];
    //     }
    //     cout << "Rank " << rank << ": nbl*zone_ptr->nbf " << nbl*zone_ptr->nbf << ",  max(blopt[blopt_i]) " << max_b << endl;
    //  }
    //  MPI_Barrier(mpi_comm);
    //}
    for (int blopt_i = 0; blopt_i<blopt_s; ++blopt_i){
      assert(blopt[blopt_i]>=0);
      assert(blopt[blopt_i]<(nbl*zone_ptr->nbf));
    }

    for (int blopt_i = 0; blopt_i<blopt_s; ++blopt_i)
      r1_bl[blopt[blopt_i]] = recv_buf_r1[blopt_i];

    delete[] recv_buf_r1;


    double dl = l_bl/nbl;
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
      bl_delta[ibf] = 0.0;
      for (int ibl = 0; ibl<nbl; ++ibl){
        if (r1_bl[ibf*nbl+ibl]>(pTotal_threshold) || (r1_bl[ibf*nbl+ibl]==-HUGE_VAL))
          break;
        bl_delta[ibf] = (ibl+0.5)*dl;//value of bl_delta
      }
    }
    delete[] r1_bl; r1_bl = NULL;
    
    //MiscUtils::dumpRange(r1_bl,nbl*zone_ptr->nbf,zone_ptr->getName()+":BL_Pt");
  }


};
#endif
