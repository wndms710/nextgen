#include "Multigrid.hpp"
#include "SvdNew.hpp"

namespace Multigrid {

  void colorCvsPadt(int8 *color_cv,int8 &ncolors,const double (*const x_cv)[3],const int8* const icv_global,
      const int ncv,const int8 ncv_global) {

    const int debug_rank = getIntParam("DEBUG_RANK",-1);

    if (mpi_rank == 0)
      cout << "Multigrid::colorCvsPadt(), requested colors: " << ncolors << endl;

    // get base 2 size 

    int mpi_size_tmp = 1;
    while (mpi_size_tmp*2 <= mpi_size)
      mpi_size_tmp *= 2;

    if (ncolors < mpi_size_tmp) {
      CWARN("Min number of colors is 2^(floor(log2(mpi_size))). Setting ncolors to " << mpi_size_tmp << ".");
      ncolors = mpi_size_tmp;
    }

    int *send_count = new int[mpi_size];
    int *send_disp  = new int[mpi_size];
    int *recv_count = new int[mpi_size];
    int *recv_disp  = new int[mpi_size];

    // first get mpi_size_tmp colors using the padt partitioner (w/o weighting)

    if (mpi_size_tmp == mpi_size) {

      int* cv_part = new int[ncv];
      calcCvPartPadt(cv_part,x_cv,ncv,mpi_comm);
      assert(color_cv);
      FOR_ICV color_cv[icv] = (int8)cv_part[icv];
      delete[] cv_part;

    }
    else {

      // split communicator to separate those ranks

      int mpi_key;
      if (mpi_rank < mpi_size_tmp)
        mpi_key = 0;
      else
        mpi_key = 1;

      MPI_Comm mpi_comm_tmp;
      MPI_Comm_split(mpi_comm,mpi_key,mpi_rank,&mpi_comm_tmp);

      if (mpi_key == 0) {
        int mpi_rank_check;
        int mpi_size_check;
        MPI_Comm_rank(mpi_comm_tmp, &mpi_rank_check);
        MPI_Comm_size(mpi_comm_tmp, &mpi_size_check);
        assert(mpi_size_check == mpi_size_tmp);
        assert(mpi_rank_check == mpi_rank);
      }
      else {
        MPI_Comm_free(&mpi_comm_tmp); 
      }

      // collect x_cv on those ranks w/in that size

      int8* cvora_tmp = new int8[mpi_size+1];
      {
        int8* cvora_tmp_ = NULL;
        MiscUtils::calcUniformDist(cvora_tmp_,ncv_global,mpi_size_tmp);
        for (int rank = 0; rank <= mpi_size_tmp; ++rank)
          cvora_tmp[rank] = cvora_tmp_[rank];
        delete[] cvora_tmp_;
        for (int rank = mpi_size_tmp+1; rank <= mpi_size; ++rank)
          cvora_tmp[rank] = cvora_tmp[rank-1];
      }
      DistributedDataExchanger *dde_tmp = new DistributedDataExchanger(icv_global,ncv,cvora_tmp,mpi_comm);
      const int ncv_tmp = cvora_tmp[mpi_rank+1]-cvora_tmp[mpi_rank];
      delete[] cvora_tmp;
      double (*x_cv_tmp)[3] = new double[ncv_tmp][3];
      dde_tmp->push(x_cv_tmp,x_cv);

      // calcCvPartPadt on those ranks

      int* cv_part_tmp = new int[ncv_tmp];
      if (mpi_rank < mpi_size_tmp) {
        calcCvPartPadt(cv_part_tmp,x_cv_tmp,ncv_tmp,mpi_comm_tmp);
        MPI_Comm_free(&mpi_comm_tmp); 
      }
      else {
        assert(ncv_tmp == 0);
      }
      delete[] x_cv_tmp;

      // distribute color back

      int* cv_part = new int[ncv];
      dde_tmp->pull(cv_part,cv_part_tmp);
      delete[] cv_part_tmp;
      delete dde_tmp;

      assert(color_cv);
      FOR_ICV color_cv[icv] = (int8)cv_part[icv];
      delete[] cv_part;

    }

    int8 current_ncolors = mpi_size_tmp;

    // now get the remaining colors by round robining 

    if (mpi_rank == 0)
      cout << " > current/total: " << current_ncolors << "/" << ncolors << endl;

    while (current_ncolors < ncolors) {

      FOR_RANK send_count[rank] = 0;
      double* send_buf_double;
      int8* send_buf_int8;
      int send_count_sum;
      for (int iter = 0; iter < 2; ++iter) {
        FOR_ICV {
          const int rank = color_cv[icv]%mpi_size; 
          if (iter == 0) {
            ++send_count[rank];
          }
          else {
            FOR_I3 send_buf_double[send_disp[rank]*3+i] = x_cv[icv][i];
            send_buf_int8[send_disp[rank]*2+0] = color_cv[icv];
            send_buf_int8[send_disp[rank]*2+1] = BitUtils::packRankBitsIndex(mpi_rank,0,icv); // rank included for checking
            ++send_disp[rank];
          }
        }

        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

        if (iter == 0) {
          send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
          send_buf_double = new double[send_count_sum*3];
          send_buf_int8   = new int8[send_count_sum*2];
        }
      }

      // exchange

      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      double * recv_buf_double = new double[recv_count_sum*3];
      FOR_RANK {
        send_count[rank] *= 3;
        send_disp[rank]  *= 3;
        recv_count[rank] *= 3;
        recv_disp[rank]  *= 3;
      }
      MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
          recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
      delete[] send_buf_double; 

      int8 * recv_buf_int8 = new int8[recv_count_sum*2];
      FOR_RANK {
        send_count[rank] = (send_count[rank]/3)*2;
        send_disp[rank]  = (send_disp[rank]/3)*2;
        recv_count[rank] = (recv_count[rank]/3)*2;
        recv_disp[rank]  = (recv_disp[rank]/3)*2;
      }
      MPI_Alltoallv(send_buf_int8,send_count,send_disp,MPI_INT8,
          recv_buf_int8,recv_count,recv_disp,MPI_INT8,mpi_comm);

      FOR_RANK {
        send_count[rank] /= 2;
        send_disp[rank]  /= 2;
        recv_count[rank] /= 2;
        recv_disp[rank]  /= 2;
      }

      // get the range in each direction 

      const int max_ncolors_rr = current_ncolors/mpi_size+1; // max colors on round-robined rank
      double (*bbox)[6] = new double[max_ncolors_rr][6];
      for (int icolor = 0; icolor < max_ncolors_rr; ++icolor) 
        FOR_I6 bbox[icolor][i] = HUGE_VAL;
      int ncolors_rr = 0;
      for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
        const int icolor = recv_buf_int8[irecv*2+0]/mpi_size;
        assert((icolor >= 0)&&(icolor < max_ncolors_rr));
        const double * x = &recv_buf_double[irecv*3];
        if (bbox[icolor][0] == HUGE_VAL) {
          ++ncolors_rr;
          if (mpi_rank == debug_rank)
            cout << " > recvd: " << recv_buf_int8[irecv*2+0] << " " << icolor << endl;
        }
        FOR_I3 bbox[icolor][i] = min(bbox[icolor][i],x[i]);
        FOR_I3 bbox[icolor][i+3] = min(bbox[icolor][i+3],-x[i]);
      }
      if (mpi_rank == debug_rank)
        cout << " > ncolors_rr, max_ncolors_rr: " << ncolors_rr << " " << max_ncolors_rr << endl; cout.flush();
      //assert((ncolors_rr == max_ncolors_rr)||(ncolors_rr == max_ncolors_rr-1));

      const double eps1 = 1.1234E-6; // include a little y and z to break sort ties (and hopefully not cause any)...
      const double eps2 = 2.7531E-6; // include a little y and z to break sort ties (and hopefully not cause any)...

      // which direction is the largest?

      int * dir = new int[ncolors_rr];
      for (int icolor = 0; icolor < ncolors_rr; ++icolor) {
        if ( (-bbox[icolor][3]-bbox[icolor][0]) >= max(-bbox[icolor][4]-bbox[icolor][1],-bbox[icolor][5]-bbox[icolor][2]) ) {
          dir[icolor] = 0;
        }
        else if ( (-bbox[icolor][4]-bbox[icolor][1]) >= max(-bbox[icolor][3]-bbox[icolor][0],-bbox[icolor][5]-bbox[icolor][2]) ) {
          dir[icolor] = 1;
        }
        else {
          dir[icolor] = 2;
        }
        if (mpi_rank == debug_rank)
          cout << " > bbmin,-bbmax: " << COUT_VEC(bbox[icolor]) << " " << COUT_VEC(&bbox[icolor][3]) << " " << dir[icolor] << endl;
      }

      vector<pair<double,int> > * diVec = new vector<pair<double,int> >[ncolors_rr];
      for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
        const int icolor = recv_buf_int8[irecv*2+0]/mpi_size;
        const double * x = &recv_buf_double[irecv*3];
        const int i = dir[icolor];
        diVec[icolor].push_back(pair<double,int8>(x[i]+eps1*x[(i+1)%3]+eps2*x[(i+2)%3],irecv));
      }
      delete[] recv_buf_double;

      for (int icolor = 0; icolor < ncolors_rr; ++icolor) {

        // get the new color and if it exists, bisect the data in the chosen direction

        int8 new_color = icolor*mpi_size+mpi_rank+current_ncolors;
        if (mpi_rank == debug_rank)
          cout << " > ic,size,new_color,old_color,current,total: " << icolor << " " << diVec[icolor].size() << " " << new_color << " " << icolor*mpi_size+mpi_rank << " " << current_ncolors << " " << ncolors << endl;
        if (new_color < ncolors) {

          sort(diVec[icolor].begin(),diVec[icolor].end());
          const int i = dir[icolor];
          double x0 = bbox[icolor][i] + eps1*bbox[icolor][(i+1)%3] + eps2*bbox[icolor][(i+2)%3];
          double x1 = -bbox[icolor][i+3] - eps1*bbox[icolor][(i+1)%3+3] - eps2*bbox[icolor][(i+2)%3+3]; // the bbox[icolor] max is negative...
          double xmid;
          int8 count[2];

          int lim = diVec[icolor].size();
          int i0 = 0;
          int i1 = lim-1;
          int iter = 0;
          while (1) {

            xmid = 0.5*(x0+x1);

            int i0_ = i0;
            int i1_ = i1;
            if (lim == 0) {
              count[0] = 0;
              count[1] = 0;
            }
            else if (diVec[icolor][lim-1].first <= xmid) {
              count[0] = lim;
              count[1] = 0;
            }
            else if (!(diVec[icolor][0].first <= xmid)) {
              count[0] = 0;
              count[1] = lim;
            }
            else {
              while (i0_ < i1_-1) {
                const int imid = (i0_+i1_)/2;
                if (diVec[icolor][imid].first <= xmid)
                  i0_ = imid;
                else
                  i1_ = imid;
              }
              count[0] = i0_+1;
              count[1] = lim-i0_-1;
            }

            if ( (count[1] > count[0]) && ((count[1]-1) >= (count[0]+1)) ) {
              x0 = xmid;
              i0 = i0_;
            }
            else if ( (count[1] < count[0]) && ((count[1]+1) <= (count[0]-1)) ) {
              x1 = xmid;
              i1 = i1_;
            }
            else {
              break;
            }

            ++iter;
            if (iter > 50) {
              cout << " Warning: 50 iters exceed in bisection routine: count[0]: " << count[0] << " count[1]: " << count[1] << " x0: " << x0 << " xmid: " << xmid << " x1: " << x1 << endl;
              if (iter > 75)
                break;  // can't seem to balance. Just move ahead anyways.
            }

          }
          if (mpi_rank == debug_rank)
            cout << " > iter: " << iter << " count[0]: " << count[0] << " count[1]: " << count[1] << " x0: " << x0 << " xmid: " << xmid << " x1: " << x1 << endl;

          // update color in recv_buf_int8

          for (int ii = 0; ii < lim; ++ii) {
            if (diVec[icolor][ii].first <= xmid) {
              const int irecv = diVec[icolor][ii].second;
              recv_buf_int8[irecv*2+0] = new_color;
            }
          }
          diVec[icolor].clear();

        }

      }
      delete[] bbox;
      delete[] dir;
      delete[] diVec;

      // update total current number of colors and send back the new colors

      current_ncolors = min(current_ncolors*2,ncolors);

      FOR_RANK {
        send_count[rank] *= 2;
        send_disp[rank]  *= 2;
        recv_count[rank] *= 2;
        recv_disp[rank]  *= 2;
      }
      MPI_Alltoallv(recv_buf_int8,recv_count,recv_disp,MPI_INT8,
          send_buf_int8,send_count,send_disp,MPI_INT8,mpi_comm);
      delete[] recv_buf_int8;

      for (int isend = 0; isend < send_count_sum; ++isend) {
        const int8 color = send_buf_int8[isend*2+0];
        assert((color >= 0)&&(color < current_ncolors));
        //const int8 rbi   = send_buf_int8[isend*2+1];
        int rank,bits,icv;
        BitUtils::unpackRankBitsIndex(rank,bits,icv,send_buf_int8[isend*2+1]);
        assert(rank == mpi_rank);
        assert(bits == 0);
        assert((icv >= 0)&&(icv < ncv));
        color_cv[icv] = color;
      }
      delete[] send_buf_int8;

      if (mpi_rank == 0)
        cout << " > current/total: " << current_ncolors << "/" << ncolors << endl;

    }

    delete[] send_count;
    delete[] send_disp;
    delete[] recv_count;
    delete[] recv_disp;

  }

  void updateCvData(int8* color_cv,const uint8* const rbicv_g,const uint8* const rbicv_g2,const int ncv,const int ncv_g,
      const int ncv_g2) {
    assert(rbicv_g);
    assert(rbicv_g2||(ncv_g==ncv_g2));
    
    int* send_count = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;
    int* send_disp = new int[mpi_size];
    int *send_buf_int = new int[ncv_g2-ncv];
    for (int iter = 0; iter < 2; ++iter) {

      for (int icv = ncv; icv < ncv_g; ++icv) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbicv_g[icv-ncv]);
        assert((rank != mpi_rank)||bits);
        if (iter == 0)
          ++send_count[rank];
        else
          send_buf_int[send_disp[rank]++] = index;
      }
      for (int icv = ncv_g; icv < ncv_g2; ++icv) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbicv_g2[icv-ncv_g]);
        assert((rank != mpi_rank)||bits);
        if (iter == 0)
          ++send_count[rank];
        else
          send_buf_int[send_disp[rank]++] = index;
      }

      // rewind...

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      assert((send_disp[mpi_size-1] + send_count[mpi_size-1]) == (ncv_g2-ncv));
    }

    // setup recv side stuff...

    int* recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    int* recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    // exchange...

    int* recv_buf_int = new int[recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
        recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; 

    // pack color...

    int8* recv_buf_int8 = new int8[recv_count_sum];
    for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
      const int icv = recv_buf_int[irecv];
      assert((icv >= 0)&&(icv < ncv));
      recv_buf_int8[irecv] = color_cv[icv];
    }
    delete[] recv_buf_int;

    // reverse exchange...

    int8* send_buf_int8 = new int8[ncv_g2-ncv];
    MPI_Alltoallv(recv_buf_int8,recv_count,recv_disp,MPI_INT8,
        send_buf_int8,send_count,send_disp,MPI_INT8,mpi_comm);
    delete[] recv_buf_int8; 

    // unpack color...

    for (int icv = ncv; icv < ncv_g; ++icv) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbicv_g[icv-ncv]);
      assert((rank != mpi_rank)||bits);
      color_cv[icv] = send_buf_int8[send_disp[rank]++];
    }
    for (int icv = ncv_g; icv < ncv_g2; ++icv) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbicv_g2[icv-ncv_g]);
      assert((rank != mpi_rank)||bits);
      color_cv[icv] = send_buf_int8[send_disp[rank]++];
    }

    // cleanup...
    delete[] send_buf_int8;
    delete[] send_count;
    delete[] send_disp;
    delete[] recv_count;
    delete[] recv_disp;

  }

  void splitOrphanedColors(int8 * color_cv,int8 &ncolor,const int* const cvocv_i,const int* const cvocv_v,
      const int8* const icv_global,const uint8* const rbicv_g,const uint8* const rbicv_g2,const int ncv,const int ncv_g,
      const int ncv_g2,const int8 ncv_global) {
    assert(rbicv_g);
    assert(rbicv_g2||(ncv_g==ncv_g2));

    // assumes color_cv is filled in the ghosts

    if (mpi_rank == 0) 
      cout << "Multigrid::splitOrphanedColors()" << endl;

    const int8 ncolor0 = ncolor;
    int8* color_cv2 = new int8[ncv_g2];
    FOR_ICV color_cv2[icv] = icv_global[icv];
    for (int icv = ncv; icv < ncv_g2; ++icv)
      color_cv2[icv] = -1;
    updateCvData(color_cv2,rbicv_g,rbicv_g2,ncv,ncv_g,ncv_g2);
    for (int icv = ncv; icv < ncv_g2; ++icv)
      assert((color_cv2[icv] >= 0)&&(color_cv2[icv] < ncv_global));

    int done = 0;
    int iter = 0;
    while (done == 0) {
      ++iter;

      int my_count = 0;
      FOR_ICV {
        int8 min_color2 = color_cv2[icv];
        for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          if (color_cv[icv] == color_cv[icv_nbr]) 
            min_color2 = min(min_color2,color_cv2[icv_nbr]);
        }
        if (min_color2 != color_cv2[icv]) {
          color_cv2[icv] = min_color2;
          ++my_count;
        }
      }
      updateCvData(color_cv2,rbicv_g,rbicv_g2,ncv,ncv_g,ncv_g2);

      int count;
      MPI_Allreduce(&my_count,&count,1,MPI_INT,MPI_SUM,mpi_comm);

      if ((iter > int(ncv_global/ncolor0))&&(mpi_rank == 0))
        cout << " > iter, ungrouped: " << iter << " " << count << endl;

      if (count == 0)
        done = 1;
    }

    // round robining on color to reindex

    int* send_count = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;
    int* send_disp = new int[mpi_size];

    int8 *send_buf_int8 = new int8[ncv_g2];
    for (int iter = 0; iter < 2; ++iter) {

      FOR_ICV_G2 {
        const int rank = color_cv2[icv]%mpi_size;
        assert((rank >= 0)&&(rank < mpi_size));
        if (iter == 0) {
          ++send_count[rank];
        }
        else {
          send_buf_int8[send_disp[rank]++] = color_cv2[icv];
        }
      }

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      assert(ncv_g2 == (send_count[mpi_size-1] + send_disp[mpi_size-1]));
    }

    // setup recv side stuff...

    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    // exchange...

    int8* recv_buf_int8 = new int8[recv_count_sum];
    MPI_Alltoallv(send_buf_int8,send_count,send_disp,MPI_INT8,
        recv_buf_int8,recv_count,recv_disp,MPI_INT8,mpi_comm);

    vector<pair<int8,int> > colorIrecvVec(recv_count_sum);
    for (int irecv = 0; irecv < recv_count_sum; ++irecv) 
      colorIrecvVec[irecv] = pair<int8,int>(recv_buf_int8[irecv],irecv);

    sort(colorIrecvVec.begin(),colorIrecvVec.end());
    int8 my_ncolor = 0;
    int current_color = -1;
    for (int ii = 0; ii < recv_count_sum; ++ii) {
      assert(current_color <= colorIrecvVec[ii].first);
      if (current_color < colorIrecvVec[ii].first) {
        ++my_ncolor;
        current_color = colorIrecvVec[ii].first;
      }
    }
    MPI_Allreduce(&my_ncolor,&ncolor,1,MPI_INT8,MPI_SUM,mpi_comm);
    assert(ncolor >= ncolor0);

    if (mpi_rank == 0) 
      cout << " > ncolor0: " << ncolor0 << ", ncolor: " << ncolor << endl;

    // reindex colors to be b/w 0 and ncolor-1...

    int8 color_offset;
    MPI_Scan(&my_ncolor,&color_offset,1,MPI_INT8,MPI_SUM,mpi_comm);
    color_offset -= my_ncolor;

    my_ncolor = -1;
    current_color = -1;
    for (int ii = 0; ii < recv_count_sum; ++ii) {
      assert(current_color <= colorIrecvVec[ii].first);
      if (current_color < colorIrecvVec[ii].first) {
        current_color = colorIrecvVec[ii].first;
        ++my_ncolor;
      }
      recv_buf_int8[colorIrecvVec[ii].second] = my_ncolor+color_offset;
    }

    // send back updated color

    MPI_Alltoallv(recv_buf_int8,recv_count,recv_disp,MPI_INT8,
        send_buf_int8,send_count,send_disp,MPI_INT8,mpi_comm);
    delete[] recv_buf_int8; 
    delete[] recv_count;
    delete[] recv_disp;

    FOR_ICV_G2 {
      const int rank = color_cv2[icv]%mpi_size;
      assert((rank >= 0)&&(rank < mpi_size));
      color_cv[icv] = send_buf_int8[send_disp[rank]++];
      assert((color_cv[icv] >= 0)&&(color_cv[icv] < ncolor));
    }
    delete[] send_buf_int8; 
    delete[] send_disp;
    delete[] send_count;
    delete[] color_cv2;

  }

  void AlgebraicCoarseGrid::buildCcIPrcomm() {

    if (mpi_rank == 0)
      cout << "AlgebraicCoarseGrid::buildCcIPrcomm()" << endl;

    // ------------------------------------------------
    // build the inactive to active communicator
    // ------------------------------------------------

    assert(ccIPrcommVec.empty());

    // build pack side using rbicc_i...
    {
      int rank_current = -1;
      int bits_current = -1;
      Prcomm * prcomm = NULL;
      for (int icc = ncc_a; icc < ncc; ++icc) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbicc_i[icc-ncc_a]);
        //assert(bits == 0);
        if (rank > rank_current) {
          // we just switched ranks. If there was a prcomm, then
          // complete the size...
          if (prcomm) {
            prcomm->pack_size = icc-prcomm->pack_offset;
            prcomm->packVec.resize(prcomm->pack_size);
            for (int icc2 = prcomm->pack_offset; icc2 < icc; ++icc2) {
              prcomm->packVec[icc2-prcomm->pack_offset] = icc2;
            }
          }
          rank_current = rank;
          bits_current = -1;
          ccIPrcommVec.push_back(Prcomm());
          prcomm = &ccIPrcommVec.back();
          prcomm->rank = rank;
          prcomm->pack_offset = icc;
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
      if (prcomm) {
        prcomm->pack_size = ncc-prcomm->pack_offset;
        prcomm->packVec.resize(prcomm->pack_size);
        for (int icc2 = prcomm->pack_offset; icc2 < ncc; ++icc2) {
          prcomm->packVec[icc2-prcomm->pack_offset] = icc2;
        }
      }
    }

    // we need to get the rbi's for the inactive cells that go with our split active cells...

    // get rbi's for the inactive cells that contribute to our active cells 
    int* send_count = new int[mpi_size];
    int* send_disp = new int[mpi_size];
    int* send_buf_int = new int[3*(ncc-ncc_a)]; // nbr index,index
    FOR_RANK send_count[rank] = 0;
    for (int iter = 0; iter < 2; ++iter) {
      for (int icc = ncc_a; icc < ncc; ++icc) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbicc_i[icc-ncc_a]); 
        //assert(rank != mpi_rank);
        //assert(bits == 0);
        assert((rank != mpi_rank)||(bits != 0));
        if (iter == 0) {
          send_count[rank] += 3;
        }
        else {
          send_buf_int[send_disp[rank]++] = index;
          send_buf_int[send_disp[rank]++] = icc;
          send_buf_int[send_disp[rank]++] = bits;
        }
      }
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
      assert((send_count[mpi_size-1] + send_disp[mpi_size-1]) == 3*(ncc-ncc_a));
    }

    // setup recv side stuff

    int* recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    //assert(recv_count[mpi_rank] == 0); 

    int* recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    int* recv_buf_int = new int[recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
        recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; 
    delete[] send_count;
    delete[] send_disp;

    assert(recv_count_sum%3 == 0);
    vector<pair<int8,int> > rbi_index_pair_vec(recv_count_sum/3);
    FOR_RANK {
      int irecv = recv_disp[rank];
      assert(irecv%3 == 0);
      while (irecv < recv_disp[rank]+recv_count[rank]) {
        assert(recv_buf_int[irecv] < ncc_a); // touches my active cell
        rbi_index_pair_vec[irecv/3] = pair<int8,int>(BitUtils::packRankBitsIndex(rank,recv_buf_int[irecv+2],recv_buf_int[irecv+1]),recv_buf_int[irecv]);
        irecv += 3;
      }
    }
    delete[] recv_buf_int;
    delete[] recv_count;
    delete[] recv_disp;
    sort(rbi_index_pair_vec.begin(),rbi_index_pair_vec.end());

    {
      int rank_current = -1;
      int bits_current = -1;
      Prcomm * prcomm = NULL;
      for (int ii = 0, ii_end=rbi_index_pair_vec.size(); ii < ii_end; ++ii) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_index_pair_vec[ii].first);
        //assert(rank != mpi_rank);
        //assert(bits == 0);
        assert((rank != mpi_rank)||(bits != 0));
        if (rank > rank_current) {
          // we just switched ranks. If there was a prcomm, then
          // complete the size...
          if (prcomm) {
            prcomm->unpack_size = ii-prcomm->unpack_offset;
            prcomm->unpackVec.resize(prcomm->unpack_size);
            for (int ii2 = prcomm->unpack_offset; ii2 < ii; ++ii2) {
              prcomm->unpackVec[ii2-prcomm->unpack_offset] = rbi_index_pair_vec[ii2].second;
            }
          }
          rank_current = rank;
          bits_current = -1;
          vector<Prcomm>::iterator it;
          for (it = ccIPrcommVec.begin(); it != ccIPrcommVec.end(); ++it) {
            if (it->rank == rank) {
              prcomm = &(*it);
              break;
            }
          }
          if (it == ccIPrcommVec.end()) {
            ccIPrcommVec.push_back(Prcomm());
            prcomm = &ccIPrcommVec.back();
            prcomm->rank = rank;
            prcomm->pack_offset = 0;
            prcomm->pack_size = 0;
          }
          prcomm->unpack_offset = ii;
        }
        else {
          assert(rank_current == rank);
          assert(prcomm);
        }
        if (bits > bits_current) {
          bits_current = bits;
        }
        else {
          assert(bits_current == bits);
        }
      }
      // we just finished. If there was a prcomm, then
      // complete the size...
      if (prcomm) {
        prcomm->unpack_size = rbi_index_pair_vec.size()-prcomm->unpack_offset;
        prcomm->unpackVec.resize(prcomm->unpack_size);
        for (int ii2 = prcomm->unpack_offset, ii2_end=rbi_index_pair_vec.size(); ii2 < ii2_end; ++ii2) {
          prcomm->unpackVec[ii2-prcomm->unpack_offset] = rbi_index_pair_vec[ii2].second;
        }
      }
    }

  }

  void AlgebraicCoarseGrid::buildCcPrcomm() {

    if (mpi_rank == 0)
      cout << "AlgebraicCoarseGrid::buildCcPrcomm()" << endl;

    // ------------------------------------------------
    // build the active to inactive/ghost communicator
    // ------------------------------------------------

    assert(ccPrcommVec.empty());

    // we need to get the rbi's for the ghost cells that go with our active cells...

    int* send_count = new int[mpi_size];
    int* send_disp = new int[mpi_size];
    int* send_buf_int = new int[2*(ncc_g-ncc)]; // nbr index,index
    FOR_RANK send_count[rank] = 0;
    for (int iter = 0; iter < 2; ++iter) {
      for (int icc = ncc; icc < ncc_g; ++icc) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbicc_g[icc-ncc]); 
        assert(rank != mpi_rank); 
        assert(bits == 0);
        if (iter == 0) {
          send_count[rank] += 2;
        }
        else {
          send_buf_int[send_disp[rank]++] = index;
          send_buf_int[send_disp[rank]++] = icc;
        }
      }
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
      assert((send_count[mpi_size-1] + send_disp[mpi_size-1]) == 2*(ncc_g-ncc));
    }

    // setup recv side stuff

    int* recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    assert(recv_count[mpi_rank] == 0); 

    int* recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    int* recv_buf_int = new int[recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
        recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; 
    delete[] send_count;
    delete[] send_disp;

    assert(recv_count_sum%2 == 0);
    vector<pair<int8,int> > rbi_index_pair_vec(recv_count_sum/2);
    set<int> auxRankSet;
    FOR_RANK {
      int irecv = recv_disp[rank];
      assert(irecv%2 == 0);
      while (irecv < recv_disp[rank]+recv_count[rank]) {
        assert(recv_buf_int[irecv] < ncc_a); // touches my active cell
        rbi_index_pair_vec[irecv/2] = pair<int8,int>(BitUtils::packRankBitsIndex(rank,0,recv_buf_int[irecv+1]),recv_buf_int[irecv]);
        auxRankSet.insert(rank);
        irecv += 2;
      }
    }
    delete[] recv_buf_int;
    delete[] recv_count;
    delete[] recv_disp;
    sort(rbi_index_pair_vec.begin(),rbi_index_pair_vec.end());

    // throw ccIPrcommVec ranks into set 
    for (int ii = 0, lim = ccIPrcommVec.size(); ii < lim; ++ii)
      auxRankSet.insert(ccIPrcommVec[ii].getRank());

    // build unpack side using rbicc_g...
    {
      set<int>::iterator iter = auxRankSet.begin();
      int rank_current = -1;
      int bits_current = -1;
      Prcomm * prcomm = NULL;
      for (int icc = ncc; icc < ncc_g; ++icc) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbicc_g[icc-ncc]);
        assert(bits == 0);
        if (rank > rank_current) {
          // we just switched ranks. If there was a prcomm, then
          // complete the size...
          if (prcomm) 
            prcomm->unpack_size = icc-prcomm->unpack_offset;
          while ((iter != auxRankSet.end())&&((*iter) <= rank)) {
            if ((*iter) < rank) {
              // add an empty Prcomm...
              ccPrcommVec.push_back(Prcomm());
              ccPrcommVec.back().rank = *iter;
              ccPrcommVec.back().unpack_offset = icc;
              ccPrcommVec.back().unpack_size = 0;
            }
            ++iter;
          }
          rank_current = rank;
          bits_current = -1;
          ccPrcommVec.push_back(Prcomm());
          prcomm = &ccPrcommVec.back();
          prcomm->rank = rank;
          prcomm->unpack_offset = icc;
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
      if (prcomm) 
        prcomm->unpack_size = ncc_g-prcomm->unpack_offset;
      // add any not represented in Prcomm at the end...
      while (iter != auxRankSet.end()) {
        // add an empty Prcomm...
        ccPrcommVec.push_back(Prcomm());
        ccPrcommVec.back().rank = *iter;
        ccPrcommVec.back().unpack_offset = ncc_g;
        ccPrcommVec.back().unpack_size = 0;
        ++iter;
      }
      assert(iter == auxRankSet.end());
    }
    auxRankSet.clear();

    {
      int rank_current = -1;
      int bits_current = -1;
      Prcomm * prcomm = NULL;
      for (int ii = 0, ii_end=rbi_index_pair_vec.size(); ii < ii_end; ++ii) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_index_pair_vec[ii].first);
        assert(rank != mpi_rank); 
        assert(bits == 0);
        if (rank > rank_current) {
          // we just switched ranks. If there was a prcomm, then
          // complete the size...
          if (prcomm) {
            prcomm->pack_size = ii-prcomm->pack_offset;
            prcomm->packVec.resize(prcomm->pack_size);
            for (int ii2 = prcomm->pack_offset; ii2 < ii; ++ii2) {
              prcomm->packVec[ii2-prcomm->pack_offset] = rbi_index_pair_vec[ii2].second;
            }
          }
          rank_current = rank;
          bits_current = -1;
          vector<Prcomm>::iterator it;
          for (it = ccPrcommVec.begin(); it != ccPrcommVec.end(); ++it) {
            if (it->rank == rank) {
              prcomm = &(*it);
              break;
            }
          }
          assert(it != ccPrcommVec.end()); // added above so should NOT add any new ones
          prcomm->pack_offset = ii;
        }
        else {
          assert(rank_current == rank);
          assert(prcomm);
        }
        if (bits > bits_current) {
          bits_current = bits;
        }
        else {
          assert(bits_current == bits);
        }
      }
      // we just finished. If there was a prcomm, then
      // complete the size...
      if (prcomm) {
        prcomm->pack_size = rbi_index_pair_vec.size()-prcomm->pack_offset;
        prcomm->packVec.resize(prcomm->pack_size);
        for (int ii2 = prcomm->pack_offset, ii2_end=rbi_index_pair_vec.size(); ii2 < ii2_end; ++ii2) {
          prcomm->packVec[ii2-prcomm->pack_offset] = rbi_index_pair_vec[ii2].second;
        }
      }
    }

  }

  void AlgebraicCoarseGrid::calcCcGrad(double (*dpdx)[3],const double *p) {

    for (int icc = 0; icc < ncc_a; ++icc) {
      const int coc_f = ccocc_i[icc];
      FOR_I3 dpdx[icc][i] = 0.0; 
      //FOR_I3 dpdx[icc][i] = ccocc_grad_coeff[coc_f][i] * p[icc];
      for (int coc = coc_f+1; coc != ccocc_i[icc+1]; ++coc) {
        const int icc_nbr = ccocc_v[coc];
        FOR_I3 dpdx[icc][i] += ccocc_grad_coeff[coc][i] * (p[icc_nbr] - p[icc]);
        //FOR_I3 dpdx[icc][i] += ccocc_grad_coeff[coc][i] * p[icc_nbr];
      }
    }

  }

  void AlgebraicCoarseGrid::prolongCcData(double* phi_cv,const double* phi_cc) {

    // make sure that you call updateCcData to update phi_cc prior to calling this routine to get it in inactive cells. 

    for (int icc = 0; icc < ncc; ++icc) {
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        phi_cv[icv] = phi_cc[icc];
      }
    }

  }

  void AlgebraicCoarseGrid::prolongCcDataAndUpdateGuess(double* phi_cv,const double* phi_cc) {

    // make sure that you call updateCcData to update phi_cc prior to calling this routine to get it in inactive cells. 

    for (int icc = 0; icc < ncc; ++icc) {
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        phi_cv[icv] += phi_cc[icc];
      }
    }

  }

  void AlgebraicCoarseGrid::prolongCcData(double* phi_cv,const double (*x_cv)[3],const double* phi_cc,const double (*dphi_cc)[3]) {

    // make sure that you call updateCcData to update phi_cc prior to calling this routine to get it in inactive cells. 

    for (int icc = 0; icc < ncc; ++icc) {
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        const double dx[3] = DIFF(x_cv[icv],x_cc[icc]);
        phi_cv[icv] = phi_cc[icc]+DOT_PRODUCT(dphi_cc[icc],dx);
      }
    }

  }

  void AlgebraicCoarseGrid::prolongCcDataAndUpdateGuess(double* phi_cv,const double (*x_cv)[3],const double* phi_cc,const double (*dphi_cc)[3]) {

    // make sure that you call updateCcData to update phi_cc prior to calling this routine to get it in inactive cells. 

    for (int icc = 0; icc < ncc; ++icc) {
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        const double dx[3] = DIFF(x_cv[icv],x_cc[icc]);
        phi_cv[icv] += phi_cc[icc]+DOT_PRODUCT(dphi_cc[icc],dx);
      }
    }

  }

  void AlgebraicCoarseGrid::prolongCcData(double (*u_cv)[3],const double (*u_cc)[3]) {

    // make sure that you call updateCcData to update u_cc prior to calling this routine to get it in inactive cells. 

    for (int icc = 0; icc < ncc; ++icc) {
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        FOR_I3 u_cv[icv][i] = u_cc[icc][i];
      }
    }

  }

  void AlgebraicCoarseGrid::prolongCcDataAndUpdateGuess(double (*u_cv)[3],const double (*u_cc)[3],const double relax) {

    // make sure that you call updateCcData to update u_cc prior to calling this routine to get it in inactive cells. 

    for (int icc = 0; icc < ncc; ++icc) {
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        FOR_I3 u_cv[icv][i] += relax*u_cc[icc][i];
      }
    }

  }

  void AlgebraicCoarseGrid::restrictCcData(double* phi_cc,const double* phi_cv) {

    // restrict inactive cells first and send to active

    for (int icc = ncc_a; icc < ncc; ++icc) {
      phi_cc[icc] = 0.0;
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        phi_cc[icc] += vol_cv[icv]*phi_cv[icv];
      }
    }
    updateCcIDataStart(phi_cc); 

    // restrict the active cells

    for (int icc = 0; icc < ncc_a; ++icc) {
      phi_cc[icc] = 0.0;
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        phi_cc[icc] += vol_cv[icv]*phi_cv[icv];
      }
    }
    updateCcIDataFinish(phi_cc);

    // normalize

    for (int icc = 0; icc < ncc_a; ++icc)
      phi_cc[icc] *= inv_vol[icc];

    // you probably will need to call updateCcData to have phi_cc in inactive/ghosts after calling this routine 

  }

  void AlgebraicCoarseGrid::restrictExtrinsicCcData(double* phiV_cc,const double* phiV_cv) {

    // restrict inactive cells first and send to active

    for (int icc = ncc_a; icc < ncc; ++icc) {
      phiV_cc[icc] = 0.0;
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        phiV_cc[icc] += phiV_cv[icv];
      }
    }
    updateCcIDataStart(phiV_cc); 

    // restrict the active cells

    for (int icc = 0; icc < ncc_a; ++icc) {
      phiV_cc[icc] = 0.0;
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        phiV_cc[icc] += phiV_cv[icv];
      }
    }
    updateCcIDataFinish(phiV_cc);

    // you probably will need to call updateCcData to have phiV_cc in inactive/ghosts after calling this routine 

  }

  void AlgebraicCoarseGrid::restrictCcData(double (*u_cc)[3],const double (*u_cv)[3]) {

    for (int icc = ncc_a; icc < ncc; ++icc) {
      FOR_I3 u_cc[icc][i] = 0.0;
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        FOR_I3 u_cc[icc][i] += vol_cv[icv]*u_cv[icv][i];
      }
    }
    updateCcIDataStart(u_cc); 

    for (int icc = 0; icc < ncc_a; ++icc) {
      FOR_I3 u_cc[icc][i] = 0.0;
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        FOR_I3 u_cc[icc][i] += vol_cv[icv]*u_cv[icv][i];
      }
    }
    updateCcIDataFinish(u_cc); 

    for (int icc = 0; icc < ncc_a; ++icc)
      FOR_I3 u_cc[icc][i] *= inv_vol[icc];

    // you probably will need to call updateCcData to have phi_cc in inactive/ghosts after calling this routine 

  }

  void AlgebraicCoarseGrid::restrictExtrinsicCcData(double (*uV_cc)[3],const double (*uV_cv)[3]) {

    for (int icc = ncc_a; icc < ncc; ++icc) {
      FOR_I3 uV_cc[icc][i] = 0.0;
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        FOR_I3 uV_cc[icc][i] += uV_cv[icv][i];
      }
    }
    updateCcIDataStart(uV_cc); 

    for (int icc = 0; icc < ncc_a; ++icc) {
      FOR_I3 uV_cc[icc][i] = 0.0;
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        FOR_I3 uV_cc[icc][i] += uV_cv[icv][i];
      }
    }
    updateCcIDataFinish(uV_cc); 

    // you probably will need to call updateCcData to have phi_cc in inactive/ghosts after calling this routine 

  }

  void AlgebraicCoarseGrid::init(const double* vol_cv,const double (*const x_cv)[3],const double* const A_cv,
      const int* const cvocv_i,const int* const cvocv_v,const int8* const icv_global,const uint8* const rbicv_g,
      const uint8* const rbicv_g2,const int ncv,const int ncv_g,const int ncv_g2,const int8 ncv_global,const int icg,
      const double agglomeration_factor,const bool split_orphaned_colors,const double amg_coeff) {

    // icg 0 behaves slightly differently, so store the icg...

    this->icg = icg;
    this->amg_coeff = amg_coeff;

    // need to keep copy of volume for restriction operator...

    this->vol_cv = vol_cv;
    
    // color base grid based to get new global indices

    ncc_global = int8(double(ncv_global)/agglomeration_factor); 

    if (mpi_rank == 0)
      cout << "AlgebraicCoarseGrid::init(), level: " << icg << ", ncc_global: " << ncc_global << ", amg_coeff: " << amg_coeff << endl;

    int8* icc_global_cv = new int8[ncv_g2];

    colorCvsPadt(icc_global_cv,ncc_global,x_cv,icv_global,ncv,ncv_global);
    //FILE * fp = fopen("color_cv1.dat","w");
    //FOR_ICV fprintf(fp,"%lld %lld %lf %lf %lf\n",icc_global_cv[icv],icv_global[icv],x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
    //fclose(fp);

    updateCvData(icc_global_cv,rbicv_g,rbicv_g2,ncv,ncv_g,ncv_g2);
    if (split_orphaned_colors) {
      splitOrphanedColors(icc_global_cv,ncc_global,cvocv_i,cvocv_v,icv_global,
          rbicv_g,rbicv_g2,ncv,ncv_g,ncv_g2,ncv_global); // NULL because we are using compact grid only
    }

    if (mpi_rank == 0)
      cout << " final ncc_global: " << ncc_global << endl;

    map<const int8,int> globalMap;
    assert(ncc == 0);
    for (int icv = 0; icv < ncv; ++icv) {
      map<const int8,int>::iterator iter = globalMap.find(icc_global_cv[icv]);
      if (iter == globalMap.end()) 
        globalMap[icc_global_cv[icv]] = ncc++;
    }
    assert(globalMap.size() == ncc);

    // figure out which coarse cells (colors) are split amongst ranks

    assert(icc_global == NULL); icc_global = new int8[ncc]; 
    for (map<const int8,int>::iterator iter = globalMap.begin(); iter != globalMap.end(); ++iter) 
      icc_global[iter->second] = iter->first;

    // send to striping so we can rectify shared/ghost icc's

    int8 *ccora_striped = NULL;
    MiscUtils::calcUniformDist(ccora_striped,ncc_global,mpi_size);

    int * send_count = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;
    for (int icc = 0; icc < ncc; ++icc) {
      const int rank = MiscUtils::getRankInXora(icc_global[icc],ccora_striped);
      send_count[rank] += 3; // mpi_rank,icc,icc_striped
    }

    int * send_disp = new int[mpi_size];
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

    int* send_buf_int = new int[send_count_sum]; 
    for (int icc = 0; icc < ncc; ++icc) {
      const int rank = MiscUtils::getRankInXora(icc_global[icc],ccora_striped);
      send_buf_int[send_disp[rank]++] = mpi_rank;
      send_buf_int[send_disp[rank]++] = icc;
      send_buf_int[send_disp[rank]++] = icc_global[icc]-ccora_striped[rank];
    }

    // rewind...

    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

    // setup recv side stuff

    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    int * recv_buf_int = new int[recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
        recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; send_buf_int = NULL;

    const int nccs = ccora_striped[mpi_rank+1]-ccora_striped[mpi_rank];
    int* ccoccs_i = new int[nccs+1];
    for (int iccs = 0; iccs < nccs; ++iccs) 
      ccoccs_i[iccs+1] = 0;
    int8* ccoccs_v_global = NULL;
    for (int iter = 0; iter < 2; ++iter) {
      for (int irecv = 0; irecv < recv_count_sum; irecv += 3) {
        const int iccs = recv_buf_int[irecv+2];
        assert((iccs >= 0)&&(iccs < nccs));
        if (iter == 0) 
          ++ccoccs_i[iccs+1];
        else 
          ccoccs_v_global[ccoccs_i[iccs]++] = BitUtils::packRankBitsIndex(recv_buf_int[irecv+0],0,recv_buf_int[irecv+1]);
      }
      if (iter == 0) {
        ccoccs_i[0] = 0;
        for (int iccs = 0; iccs < nccs; ++iccs) ccoccs_i[iccs+1] += ccoccs_i[iccs];
        ccoccs_v_global = new int8[ccoccs_i[nccs]];
      }
      else {
        // rewind
        for (int iccs = nccs; iccs > 0; --iccs)
          ccoccs_i[iccs] = ccoccs_i[iccs-1];
        ccoccs_i[0] = 0;
      }
    }
    delete[] recv_buf_int; recv_buf_int = NULL;

    FOR_RANK send_count[rank] = 0;
    for (int iter = 0; iter < 2; ++iter) {
      for (int iccs = 0; iccs < nccs; ++iccs) {
        const int nbro = ccoccs_i[iccs+1]-ccoccs_i[iccs]-1; // bros have same icc_global (exclude yourself)
        for (int cocs = ccoccs_i[iccs]; cocs < ccoccs_i[iccs+1]; ++cocs) {
          int rank,bits,index; 
          BitUtils::unpackRankBitsIndex(rank,bits,index,ccoccs_v_global[cocs]);
          assert(bits == 0);
          if (iter == 0) {
            send_count[rank] += 2+2*nbro; // icc,nbro,(rank_bro,icc_bro)
          }
          else {
            send_buf_int[send_disp[rank]++] = index; 
            send_buf_int[send_disp[rank]++] = nbro;
            for (int cocs_ = ccoccs_i[iccs]; cocs_ < ccoccs_i[iccs+1]; ++cocs_) {
              if (cocs_ != cocs) {
                int rank_,bits_,index_; 
                BitUtils::unpackRankBitsIndex(rank_,bits_,index_,ccoccs_v_global[cocs_]);
                assert(bits_ == 0);
                send_buf_int[send_disp[rank]++] = rank_; 
                send_buf_int[send_disp[rank]++] = index_;
              }
            }
          }
        }
      }

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      if (iter == 0) {
        send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
        assert(send_buf_int == NULL); 
        send_buf_int = new int[send_count_sum];
      }
    }
    delete[] ccoccs_i;
    delete[] ccoccs_v_global;

    // setup recv side stuff

    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

    recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
        recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; send_buf_int = NULL;

    int* brocc_i = new int[ncc+1];
    for (int icc = 0; icc < ncc; ++icc) 
      brocc_i[icc+1] = 0;
    int8* brocc_v_global = NULL;
    for (int iter = 0; iter < 2; ++iter) {
      FOR_RANK {
        int irecv = recv_disp[rank];
        while (irecv < recv_disp[rank]+recv_count[rank]) {
          const int icc = recv_buf_int[irecv++];
          const int nbro = recv_buf_int[irecv++];
          if (iter == 0) {
            brocc_i[icc+1] += nbro;
            irecv += 2*nbro;
          }
          else {
            for (int ibro = 0; ibro < nbro; ++ibro) {
              const int rank = recv_buf_int[irecv++];
              const int index = recv_buf_int[irecv++];
              brocc_v_global[brocc_i[icc]++] = BitUtils::packRankBitsIndex(rank,0,index);
            }
          }
        }
      }
      if (iter == 0) {
        brocc_i[0] = 0;
        for (int icc = 0; icc < ncc; ++icc) brocc_i[icc+1] += brocc_i[icc];
        brocc_v_global = new int8[brocc_i[ncc]];
      }
      else {
        // rewind
        for (int icc = ncc; icc > 0; --icc)
          brocc_i[icc] = brocc_i[icc-1];
        brocc_i[0] = 0;
      }
    }
    delete[] recv_buf_int; recv_buf_int = NULL;

    // now figure out which cells are split/active (min rank wins)

    vector<pair<int8,int> > rbiVec(ncc);
    assert(ncc_a == 0);
    assert(ncc_in == 0);
    for (int icc = 0; icc < ncc; ++icc) {
      rbiVec[icc].second = icc;
      int min_rank = mpi_rank;
      if (brocc_i[icc+1]-brocc_i[icc] == 0) {
        rbiVec[icc].first = -2; // just set to -2 if internal
        ++ncc_in; // not shared so this is an internal coarse cell
      }
      else {
        rbiVec[icc].first = -1; // just set to -1 if active on this rank
        for (int boc = brocc_i[icc]; boc < brocc_i[icc+1]; ++boc) {
          int rank,bits,index;
          BitUtils::unpackRankBitsIndex(rank,bits,index,brocc_v_global[boc]);
          assert(rank != mpi_rank);
          if (rank < min_rank) {
            rbiVec[icc].first = brocc_v_global[boc];
            min_rank = rank;
          }
        }
        if (min_rank == mpi_rank) {
          assert(rbiVec[icc].first == -1);
          ++ncc_a; // min rank wins, so this is the active coarse cell
        }
        else {
          assert(rbiVec[icc].first >= 0);
        }
      }
    }
    ncc_a += ncc_in;
    delete[] brocc_i;
    delete[] brocc_v_global; // do we need to keep/reorder these for active communicator

    // now we need to reorder coarse cells like internal,active,inactive (rbi ordered)
    sort(rbiVec.begin(),rbiVec.end());

    int *reorder_cc = new int[ncc];
    for (int ii = 0; ii < ncc; ++ii)
      reorder_cc[rbiVec[ii].second] = ii;
    const int ncc_check = ncc;
    const int ncc_in_check = ncc_in;
    const int ncc_a_check = ncc_a;
    ncc = ncc_a;
    ncc_a = ncc_in;
    ncc_in = 0;
    for (int icc = 0; icc < ncc_check; ++icc) {
      if (rbiVec[icc].first == -2)
        ++ncc_in;
      else if (rbiVec[icc].first == -1)
        ++ncc_a;
      else 
        ++ncc;
    }
    assert(ncc == ncc_check);
    assert(ncc_in == ncc_in_check);
    assert(ncc_a == ncc_a_check);

    int8 my_ncc = (int8)ncc_a;
    int8 ncc_global_check;
    MPI_Reduce(&my_ncc,&ncc_global_check,1,MPI_INT8,MPI_SUM,0,mpi_comm);
    if (mpi_rank == 0) 
      assert(ncc_global_check == ncc_global);

    assert(rbicc_i == NULL); rbicc_i = new uint8[ncc-ncc_a];
    for (map<const int8,int>::iterator iter = globalMap.begin(); iter != globalMap.end(); ++iter) {
      const int icc_old = iter->second; 
      const int icc = reorder_cc[icc_old];
      icc_global[icc] = iter->first;
      iter->second = icc;
      if (icc >= ncc_a) { 
        assert(rbiVec[icc].first >= 0);
        rbicc_i[icc-ncc_a] = rbiVec[icc].first;
      }
      else if (icc >= ncc_in) {
        assert(rbiVec[icc].first == -1);
      }
      else {
        assert(rbiVec[icc].first == -2);
      }
    }
    rbiVec.clear();

    // we need to update rbicc_i to respect the reorder on nbr rank

    FOR_RANK send_count[rank] = 0;
    int * inactive_index = new int[ncc-ncc_a];
    for (int iter = 0; iter < 2; ++iter) {
      for (int icc = ncc_a; icc < ncc; ++icc) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbicc_i[icc-ncc_a]);
        if (iter == 0) {
          ++send_count[rank];
        }
        else {
          inactive_index[send_disp[rank]] = icc-ncc_a;
          send_buf_int[send_disp[rank]++] = index;
        }
      }

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      if (iter == 0) {
        send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
        assert(send_buf_int == NULL); 
        send_buf_int = new int[send_count_sum];
      }

    }

    // setup recv side stuff

    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

    recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
        recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);

    for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
      const int icc_old = recv_buf_int[irecv]; assert((icc_old >= 0)&&(icc_old < ncc)); 
      const int icc = reorder_cc[icc_old]; assert((icc >= ncc_in)&&(icc < ncc_a)); // should be active 
      recv_buf_int[irecv] = icc;
    }
    delete[] reorder_cc; reorder_cc = NULL;

    // send back

    MPI_Alltoallv(recv_buf_int,recv_count,recv_disp,MPI_INT,
        send_buf_int,send_count,send_disp,MPI_INT,mpi_comm);
    delete[] recv_buf_int; recv_buf_int = NULL;

    // overwrite rbicc_i

    FOR_RANK {
      int isend = send_disp[rank];
      while (isend < send_disp[rank]+send_count[rank]) {
        const int8 rbi = BitUtils::packRankBitsIndex(rank,0,send_buf_int[isend]);
        rbicc_i[inactive_index[isend]] = rbi;
        ++isend;
      }
    }
    delete[] send_buf_int; send_buf_int = NULL;
    delete[] inactive_index;

    assert(cvocc_i == NULL); cvocc_i = new int[ncc+1];
    for (int icc = 0; icc < ncc; ++icc) 
      cvocc_i[icc+1] = 0;
    for (int iter = 0; iter < 2; ++iter) {
      for (int icv = 0; icv < ncv; ++icv) {
        map<const int8,int>::iterator iter2 = globalMap.find(icc_global_cv[icv]);
        assert(iter2 != globalMap.end());
        const int icc = iter2->second;
        if (iter == 0) 
          ++cvocc_i[icc+1];
        else 
          cvocc_v[cvocc_i[icc]++] = icv;
      }
      if (iter == 0) {
        cvocc_i[0] = 0;
        for (int icc = 0; icc < ncc; ++icc) 
          cvocc_i[icc+1] += cvocc_i[icc];
        assert(ncv == cvocc_i[ncc]);
        assert(cvocc_v == NULL); cvocc_v = new int[cvocc_i[ncc]];
      }
      else {
        // rewind
        for (int icc = ncc; icc > 0; --icc)
          cvocc_i[icc] = cvocc_i[icc-1];
        cvocc_i[0] = 0;
      }
    }

    MiscUtils::dumpRange(&ncc_in_check,1,"ncc internal");
    MiscUtils::dumpRange(&ncc_a,1,"ncc active");
    MiscUtils::dumpRange(&ncc,1,"ncc");

    // now lets get the ghosts

    assert(ccocv == NULL); ccocv = new int[ncv_g2]; // fill ghosts once we have them
    for (int icv = 0; icv < ncv_g2; ++icv) 
      ccocv[icv] = -1;
    for (int icc = 0; icc < ncc; ++icc) {
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        ccocv[icv] = icc;
      }
    }
    for (int icv = 0; icv < ncv; ++icv) 
      assert((ccocv[icv] >= 0)&&(ccocv[icv] < ncc));
    ncc_g = ncc;
    assert(globalMap.size() == ncc); 
    for (int icv = 0; icv < ncv; ++icv) {
      assert(icc_global[ccocv[icv]] == icc_global_cv[icv]);
      assert(icv == cvocv_v[cvocv_i[icv]]);
      for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc]; assert((icv_nbr >= 0)&&(icv_nbr < ncv_g2));
        if (icv_nbr >= ncv) {
          if (icc_global_cv[icv] != icc_global_cv[icv_nbr]) {
            map<const int8,int>::iterator iter = globalMap.find(icc_global_cv[icv_nbr]);
            if (iter == globalMap.end()) {
              ccocv[icv_nbr] = ncc_g;
              globalMap[icc_global_cv[icv_nbr]] = ncc_g++; // add in the traditional ghosts
            }
            else {
              ccocv[icv_nbr] = iter->second;
            }
          }
          else {
            ccocv[icv_nbr] = ccocv[icv];
          }
        }
      }
    }
    delete[] icc_global_cv;
    assert(globalMap.size() == ncc_g);
    // ccocv may NOT touch all ghosts on later levels...
    //for (int icv = 0; icv < ncv_g2; ++icv) 
    //  assert((ccocv[icv] >= 0)&&(ccocv[icv] < ncc_g));
    for (int coc = 0; coc < cvocv_i[ncv]; ++coc) {
      const int icv_nbr = cvocv_v[coc]; assert((icv_nbr >= 0)&&(icv_nbr < ncv_g2));
      const int icc_nbr = ccocv[icv_nbr]; assert((icc_nbr >= 0)&&(icc_nbr < ncc_g));
    }

    // add inactive rbi's 
    assert(rbicc_g == NULL); rbicc_g = new uint8[ncc_g-ncc];
    for (int icc = ncc; icc < ncc_g; ++icc) 
      rbicc_g[icc-ncc] = numeric_limits<uint8>::max();
    int8* rbicc_cv = new int8[ncv_g2];
    for (int icv = 0; icv < ncv_g2; ++icv)
      rbicc_cv[icv] = -1;
    for (int icc = 0; icc < ncc_a; ++icc) {
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        rbicc_cv[icv] = BitUtils::packRankBitsIndex(mpi_rank,0,icc);
      }
    }
    for (int icc = ncc_a; icc < ncc; ++icc) {
      for (int coc = cvocc_i[icc]; coc < cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        rbicc_cv[icv] = rbicc_i[icc-ncc_a];
      }
    }
    for (int icv = 0; icv < ncv; ++icv)
      assert(rbicc_cv[icv] >= 0);
    updateCvData(rbicc_cv,rbicv_g,rbicv_g2,ncv,ncv_g,ncv_g2);
    for (int icv = ncv; icv < ncv_g2; ++icv)
      assert(rbicc_cv[icv] >= 0);
    for (int icv = ncv; icv < ncv_g2; ++icv) {
      const int icc = ccocv[icv];
      if (icc >= ncc) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbicc_cv[icv]);
        if (rank != mpi_rank) {
          if (rbicc_g[icc-ncc] == numeric_limits<uint8>::max()) {
            rbicc_g[icc-ncc] = rbicc_cv[icv];
          }
          // update with min rank rbi (active)
          else if (rbicc_g[icc-ncc] != rbicc_cv[icv]) {
            int rank0,bits0,index0;
            BitUtils::unpackRankBitsIndex(rank0,bits0,index0,rbicc_g[icc-ncc]); 
            //assert(rank != rank0); 
            if (rank < rank0)
              rbicc_g[icc-ncc] = rbicc_cv[icv];
            else if ((rank == rank0)&&(index < index0))
              rbicc_g[icc-ncc] = rbicc_cv[icv];
          }
        }
      }
    }
    delete[] rbicc_cv;

    delete[] icc_global; icc_global = new int8[ncc_g];
    for (map<const int8,int>::iterator iter = globalMap.begin(); iter != globalMap.end(); ++iter)
      icc_global[iter->second] = iter->first;

    set<int>* nbrSet = new set<int>[ncc]; 
    for (int icv = 0; icv < ncv; ++icv) {
      const int icc = ccocv[icv]; assert((icc >= 0)&&(icc < ncc));
      for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc]; assert((icv_nbr >= 0)&&(icv_nbr < ncv_g2));
        const int icc_nbr = ccocv[icv_nbr]; assert((icc_nbr >= 0)&&(icc_nbr < ncc_g));
        if (icc != icc_nbr) {
          if (icc < ncc)
            nbrSet[icc].insert(icc_nbr);
          if (icc_nbr < ncc)
            nbrSet[icc_nbr].insert(icc);
        }
      }
    }

    // pack inactive icc active nbrs and send them to the associated active icc 

    FOR_RANK send_count[rank] = 0;
    for (int iter = 0; iter < 2; ++iter) {
      for (int icc = ncc_a; icc < ncc; ++icc) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbicc_i[icc-ncc_a]);
        if (iter == 0) {
          send_count[rank] += 2+4*nbrSet[icc].size(); // index,ncoc,{nbr,rank_nbr,nbr_striped,nbr_rank_striped}
        }
        else {
          send_buf_int[send_disp[rank]++] = index;
          send_buf_int[send_disp[rank]++] = nbrSet[icc].size(); 
          for (set<int>::iterator it = nbrSet[icc].begin(); it != nbrSet[icc].end(); ++it) {
            const int icc_nbr = *it;
            const int rank_striped = MiscUtils::getRankInXora(icc_global[icc_nbr],ccora_striped);
            if (icc_nbr >= ncc_a) {
              int rank_nbr,bits_nbr,index_nbr;
              if (icc_nbr < ncc) 
                BitUtils::unpackRankBitsIndex(rank_nbr,bits_nbr,index_nbr,rbicc_i[icc_nbr-ncc_a]);
              else 
                BitUtils::unpackRankBitsIndex(rank_nbr,bits_nbr,index_nbr,rbicc_g[icc_nbr-ncc]);
              send_buf_int[send_disp[rank]++] = index_nbr; 
              send_buf_int[send_disp[rank]++] = rank_nbr;
            }
            else {
              send_buf_int[send_disp[rank]++] = icc_nbr;
              send_buf_int[send_disp[rank]++] = mpi_rank;
            }
            send_buf_int[send_disp[rank]++] = icc_global[icc_nbr]-ccora_striped[rank_striped]; 
            send_buf_int[send_disp[rank]++] = rank_striped;
          }
        }
      }

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      if (iter == 0) {
        send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
        assert(send_buf_int == NULL); 
        send_buf_int = new int[send_count_sum];
      }
    }
    //delete[] ccocc0_i;
    //delete[] ccocc0_v;

    // setup recv side stuff

    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

    recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
        recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; send_buf_int = NULL;

    // first count unique ghosts,store ghost rbi's and add to set (there could be duplicates)

    const int ncc_g0 = ncc_g;
    for (int iter = 0; iter < 2; ++iter) {
      FOR_RANK {
        int irecv = recv_disp[rank];
        while (irecv < recv_disp[rank]+recv_count[rank]) {
          const int icc = recv_buf_int[irecv++]; assert((icc >= 0)&&(icc < ncc_a));
          const int ncoc = recv_buf_int[irecv++];
          for (int coc = 0; coc < ncoc; ++coc) {
            const int icc_nbr = recv_buf_int[irecv++]; 
            const int rank_nbr = recv_buf_int[irecv++]; 
            const int8 rbi = BitUtils::packRankBitsIndex(rank_nbr,0,icc_nbr);
            const int icc_nbr_striped = recv_buf_int[irecv++];
            const int rank_striped = recv_buf_int[irecv++];
            const int8 icc_nbr_global = icc_nbr_striped+ccora_striped[rank_striped];
            assert(icc_nbr_global != icc_global[icc]);
            map<const int8,int>::iterator iter2 = globalMap.find(icc_nbr_global);
            if (iter == 0) {
              if (iter2 == globalMap.end()) {
                nbrSet[icc].insert(ncc_g);
                globalMap[icc_nbr_global] = ncc_g++;
              }
              else {
                nbrSet[icc].insert(iter2->second);
              }
            }
            else {
              if (iter2->second >= ncc) {
                if (rbicc_g[iter2->second-ncc] == numeric_limits<uint8>::max()) {
                  rbicc_g[iter2->second-ncc] = rbi;
                }
                // update with min rank rbi (active)
                else if (rbicc_g[iter2->second-ncc] != rbi) {
                  int rank0,bits0,index0;
                  BitUtils::unpackRankBitsIndex(rank0,bits0,index0,rbicc_g[iter2->second-ncc]); 
                  if (rank_nbr < rank0)
                    rbicc_g[iter2->second-ncc] = rbi;
                }
              }
            }
          }
        }
      }

      if (iter == 0) {
        // add inactive rbi's 
        uint8* rbicc_g0 = rbicc_g;
        rbicc_g = new uint8[ncc_g-ncc];
        for (int icc = ncc; icc < ncc_g0; ++icc)
          rbicc_g[icc-ncc] = rbicc_g0[icc-ncc];
        delete[] rbicc_g0;
        for (int icc = ncc_g0; icc < ncc_g; ++icc) 
          rbicc_g[icc-ncc] = numeric_limits<uint8>::max();
      }
    }
    delete[] recv_buf_int; recv_buf_int = NULL;
    for (int icc = ncc; icc < ncc_g; ++icc) 
      assert(rbicc_g[icc-ncc] != numeric_limits<uint8>::max());

    // add in new ghosts to icc_global

    delete[] icc_global; icc_global = new int8[ncc_g];
    for (map<const int8,int>::iterator iter = globalMap.begin(); iter != globalMap.end(); ++iter)
      icc_global[iter->second] = iter->first;

    rbiVec.resize(ncc_g-ncc);
    for (int icc = ncc; icc < ncc_g; ++icc)
      rbiVec[icc-ncc] = pair<int8,int>(rbicc_g[icc-ncc],icc);
    sort(rbiVec.begin(),rbiVec.end());

    // reorder ncc:ncc_g data
    assert(reorder_cc == NULL); reorder_cc = new int[ncc_g-ncc];
    int8* icc_global0 = icc_global;
    icc_global = new int8[ncc_g];
    for (int icc = 0; icc < ncc; ++icc) 
      icc_global[icc] = icc_global0[icc];
    for (int ii = 0, lim = ncc_g-ncc; ii < lim; ++ii) {
      const int icc = ii+ncc; 
      const int icc_old = rbiVec[ii].second; assert((icc_old >= ncc)&&(icc_old < ncc_g));
      map<const int8,int>::iterator iter = globalMap.find(icc_global0[icc_old]);
      assert(iter != globalMap.end());
      iter->second = icc;
      rbicc_g[ii] = rbiVec[ii].first;
      reorder_cc[icc_old-ncc] = icc;
      icc_global[icc] = iter->first;
    }
    delete[] icc_global0;
    rbiVec.clear();

    // update ccocv and cvocc_i/v
    for (int icv = 0; icv < ncv_g2; ++icv) {
      if (ccocv[icv] >= ncc)
        ccocv[icv] = reorder_cc[ccocv[icv]-ncc];
    }
    for (int icc = 0; icc < ncc; ++icc) cvocc_i[icc+1] = 0;
    for (int iter = 0; iter < 2; ++iter) {
      for (int icv = 0; icv < ncv; ++icv) {
        const int icc = ccocv[icv];
        if (iter == 0) 
          ++cvocc_i[icc+1];
        else 
          cvocc_v[cvocc_i[icc]++] = icv;
      }
      if (iter == 0) {
        cvocc_i[0] = 0;
        for (int icc = 0; icc < ncc; ++icc) cvocc_i[icc+1] += cvocc_i[icc];
      }
      else {
        // rewind
        for (int icc = ncc; icc > 0; --icc)
          cvocc_i[icc] = cvocc_i[icc-1];
        cvocc_i[0] = 0;
      }
    }

    // rbi check...
    {
      FOR_RANK send_count[rank] = 0;
      for (int iter = 0; iter < 2; ++iter) {
        for (int icc = ncc_a; icc < ncc_g; ++icc) {
          int rank,bits,index;
          if (icc < ncc)
            BitUtils::unpackRankBitsIndex(rank,bits,index,rbicc_i[icc-ncc_a]); 
          else
            BitUtils::unpackRankBitsIndex(rank,bits,index,rbicc_g[icc-ncc]); 
          if (iter == 0) {
            ++send_count[rank];
          }
          else {
            send_buf_int[send_disp[rank]++] = index;
          }
        }
        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

        if (iter == 0) {
          send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
          assert(send_buf_int == NULL); 
          send_buf_int = new int[send_count_sum];
        }
      }

      // setup recv side stuff

      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

      recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
      assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum];
      MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
          recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
      delete[] send_buf_int; send_buf_int = NULL;

      for (int irecv = 0; irecv < recv_count_sum; ++irecv)
        assert(recv_buf_int[irecv] < ncc_a);
      delete[] recv_buf_int; recv_buf_int = NULL;
    }

    // build ccocc using nbrSet

    assert(ccocc_i == NULL); ccocc_i = new int[ncc_a+1];
    for (int icc = 0; icc < ncc_a; ++icc) 
      ccocc_i[icc+1] = 1+nbrSet[icc].size();
    ccocc_i[0] = 0;
    for (int icc = 0; icc < ncc_a; ++icc) ccocc_i[icc+1] += ccocc_i[icc];
    assert(ccocc_v == NULL); ccocc_v = new int[ccocc_i[ncc_a]];
    for (int icc = 0; icc < ncc_a; ++icc) {
      ccocc_v[ccocc_i[icc]++] = icc;
      for (set<int>::iterator iter = nbrSet[icc].begin(); iter != nbrSet[icc].end(); ++iter) {
        if (*iter < ncc)
          ccocc_v[ccocc_i[icc]++] = *iter;
        else 
          ccocc_v[ccocc_i[icc]++] = reorder_cc[*iter-ncc];
      }
      nbrSet[icc].clear();
    }
    delete[] reorder_cc;
    delete[] nbrSet;
    for (int icc = ncc_a; icc > 0; --icc)
      ccocc_i[icc] = ccocc_i[icc-1];
    ccocc_i[0] = 0;

    // build ccIPrcomm so we can check it...

    buildCcIPrcomm();

    // cc geom data

    assert(x_cc == NULL); x_cc = new double[ncc_g][3];
    assert(vol_cc == NULL); vol_cc = new double[ncc_g];
    for (int icc = 0; icc < ncc; ++icc) {
      FOR_I3 x_cc[icc][i] = 0.0;
      vol_cc[icc] = 0.0;
      for (int coc = cvocc_i[icc]; coc != cvocc_i[icc+1]; ++coc) {
        const int icv = cvocc_v[coc];
        FOR_I3 x_cc[icc][i] += vol_cv[icv]*x_cv[icv][i];
        vol_cc[icc] += vol_cv[icv];
      }
    }

    double * vol_cc_check = new double[ncc];
    double (*x_cc_check)[3] = new double[ncc][3];
    for (int icc = 0; icc < ncc; ++icc) {
      vol_cc_check[icc] = vol_cc[icc];
      FOR_I3 x_cc_check[icc][i] = x_cc[icc][i];
    }
    DistributedDataExchanger * dde_cc_striped = new DistributedDataExchanger(icc_global,ncc,ccora_striped);
    delete[] ccora_striped;
    double (*x_ccs)[3] = new double[nccs][3];
    for (int iccs = 0; iccs < nccs; ++iccs)
      FOR_I3 x_ccs[iccs][i] = 0.0;
    dde_cc_striped->push(x_ccs,x_cc_check,ADD_DATA);
    dde_cc_striped->pull(x_cc_check,x_ccs);
    delete[] x_ccs;
    double* vol_ccs = new double[nccs];
    for (int iccs = 0; iccs < nccs; ++iccs)
      vol_ccs[iccs] = 0.0;
    dde_cc_striped->push(vol_ccs,vol_cc_check,ADD_DATA);
    dde_cc_striped->pull(vol_cc_check,vol_ccs);
    delete dde_cc_striped;
    delete[] vol_ccs;
    for (int icc = 0; icc < ncc; ++icc) {
      assert(vol_cc_check[icc] > 0);
      FOR_I3 x_cc_check[icc][i] /= vol_cc_check[icc];
    }

    updateCcIData(vol_cc);
    updateCcIData(x_cc);
    for (int icc = 0; icc < ncc_a; ++icc) {
      assert(vol_cc[icc] > 0.0);
      FOR_I3 x_cc[icc][i] /= vol_cc[icc];
    }

    for (int icc = 0; icc < ncc_a; ++icc) {
      FOR_I3 x_cc_check[icc][i] -= x_cc[icc][i];
      vol_cc_check[icc] -= vol_cc[icc];
    }
    MiscUtils::dumpRange(vol_cc_check,ncc_a,"vol_cc check");
    delete[] vol_cc_check;
    MiscUtils::dumpRange(x_cc_check,ncc_a,"x_cc check");

    double my_buf[8]; FOR_I8 my_buf[i] = 0;
    for (int icv = 0; icv < ncv; ++icv) {
      my_buf[0] += vol_cv[icv];
      FOR_I3 my_buf[1+i] += vol_cv[icv]*x_cv[icv][i];
    }
    for (int icc = 0; icc < ncc_a; ++icc) {
      my_buf[4] += vol_cc[icc];
      FOR_I3 my_buf[5+i] += vol_cc[icc]*x_cc[icc][i];
    }
    double buf[8];
    MPI_Reduce(my_buf,buf,8,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

    if (mpi_rank == 0) {
      FOR_I3 buf[1+i] /= buf[0];
      FOR_I3 buf[5+i] /= buf[4];
      cout << " > vol diff: " << buf[0]-buf[4] << ", centroid diff: " << buf[5]-buf[1] << " " << buf[6]-buf[2] << " " << buf[7]-buf[3] << endl;
    }

    buildCcPrcomm();
    for (int icc = ncc_a; icc < ncc_g; ++icc) {
      vol_cc[icc] = 1.0E+20;
      FOR_I3 x_cc[icc][i] = 1.0E+20;
    }
    updateCcData(vol_cc);
    updateCcData(x_cc); // no bits for now
    MiscUtils::dumpRange(vol_cc,ncc_g,"vol_cc");
    MiscUtils::dumpRange(x_cc,ncc_a,"x_cc");

    delete[] send_count;
    delete[] send_disp;
    delete[] recv_count;
    delete[] recv_disp;

    assert(inv_vol == NULL);
    inv_vol = new double[ncc_g];
    for (int icc = 0; icc < ncc_g; ++icc)
      inv_vol[icc] = 1.0/vol_cc[icc];

    restrictCcData(x_cc_check,x_cv);
    for (int icc = 0; icc < ncc_a; ++icc) 
      FOR_I3 x_cc_check[icc][i] -= x_cc[icc][i];
    MiscUtils::dumpRange(x_cc_check,ncc_a,"x_cc check2");
    delete[] x_cc_check;

    // now build coarse operator...
    
    assert(A_cc == NULL);
    A_cc = new double[ccocc_i[ncc_a]];
    assert(inv_diag == NULL); 
    inv_diag = new double[ncc_g];
    update(A_cv,cvocv_i,cvocv_v,rbicv_g,rbicv_g2,ncv,ncv_g,ncv_g2);

    assert(err == NULL); 
    err = new double[ncc_g];
    assert(res == NULL); 
    res = new double[ncc_g];

    {
      // build gradient operator...

      assert(ccocc_grad_coeff == NULL);
      ccocc_grad_coeff = new double[ccocc_i[ncc_a]][3];
      const double sigma_epsilon = getDoubleParam("SIGMA_EPSILON",0.45);
      int (*flag_cc)[3] = new int[ncc_a][3];
      for (int icc = 0; icc < ncc_a; ++icc) {

        double A_grad[9];
        FOR_I3 FOR_J3 A_grad[j*3+i] = 0.0;
        const int coc_f = ccocc_i[icc];
        for (int coc = coc_f+1; coc < ccocc_i[icc+1]; ++coc) {
          const int inb = ccocc_v[coc]; 
          const double dx[3] = DIFF(x_cc[inb],x_cc[icc]); // periodicity?
          const double dr = MAG(dx);
          //double wi = dr*dr*dr; // should be 2nd order on certain configurations
          double wi = dr*dr; // O(1) sigma
          //double wi = 1.0;
          assert(wi > 0.0);
          wi = 1.0/wi;
          FOR_I3 FOR_J3 {
            const double coeff = wi*dx[i]*dx[j];
            A_grad[j*3+i] += coeff;
            if (i != j)
              A_grad[i*3+j] += coeff;
          }
        }

        // use svd to invert A_grad to handle rank deficiency...

        double U[9], V[9], sigma[3];
        double svd_eps_tol;
        calcSvd33(U,V,sigma,A_grad,svd_eps_tol);

        // now supply the regularized inverse .. let V <--- V\Sigma^+

        double sigma_inv[3];
        FOR_I3 {
          assert( sigma[i] >= 0.0);
          if (sigma[i] < sigma_epsilon) {
            sigma_inv[i] = 0.0; // cant trust the gradient in this dir..
            flag_cc[icc][i] = 1;
          }
          else { 
            sigma_inv[i] = 1.0/sigma[i];
            flag_cc[icc][i] = 0;
          }
        }
        //cout << "SIGMA: " << COUT_VEC(sigma) << endl;
        FOR_I3 FOR_J3 V[j*3+i] *= sigma_inv[j];

        // finally supply the inverse for the gradient coefficients ..

        FOR_I3 ccocc_grad_coeff[coc_f][i] = 0.0;
        for (int coc = coc_f+1; coc < ccocc_i[icc+1]; ++coc) {
          const int inb = ccocc_v[coc]; 
          const double dx[3] = DIFF(x_cc[inb],x_cc[icc]); // periodicity?
          const double dr = MAG(dx);
          //double wi = dr*dr*dr; // should be 2nd order on certain configurations
          double wi = dr*dr; // O(1) sigma
          //double wi = 1.0;
          assert(wi > 0.0);
          wi = 1.0/wi;
          double tmp[3];
          FOR_I3 tmp[i] = wi*(U[i*3+0]*dx[0] + U[i*3+1]*dx[1] + U[i*3+2]*dx[2]);
          FOR_I3 ccocc_grad_coeff[coc][i] = V[0*3+i]*tmp[0] + V[1*3+i]*tmp[1] + V[2*3+i]*tmp[2];
          FOR_I3 ccocc_grad_coeff[coc_f][i] -= ccocc_grad_coeff[coc][i];
        }

      }

      {
        const double grad_check[3] = { 1.1234, -1.3243, 1.5321 }; // some order-1 gradient
        double * phi               = new double[ncc_g];
        double (*grad_phi)[3]      = new double[ncc_a][3];
        for (int icc = 0; icc < ncc_a; ++icc) {
          phi[icc] = DOT_PRODUCT(x_cc[icc],grad_check);
          //phi[icc] = 1.12;
        }
        updateCcData(phi);

        calcCcGrad(grad_phi,phi);
        delete[] phi;

        for (int icc = 0; icc < ncc_a; ++icc) {
          //cout << COUT_VEC(x_cc[icc]) << " " << COUT_VEC(grad_phi[icc]) << " ";
          FOR_I3 {
            grad_phi[icc][i] /= grad_check[i];
            grad_phi[icc][i] -= 1.0;
          }
          //cout << COUT_VEC(grad_phi[icc]) << endl;
        }
        MiscUtils::dumpRange(grad_phi,ncc_a,"grad err");
        int8 my_count = 0;
        for (int icc = 0; icc < ncc_a; ++icc) {
          FOR_I3 {
            if (flag_cc[icc][i] == 1) {
              grad_phi[icc][i] = 0.0;
              ++my_count;
            }
          }
        }
        delete[] flag_cc;
        MiscUtils::dumpRange(grad_phi,ncc_a,"grad err unflagged");
        delete[] grad_phi;
        int8 count;
        MPI_Reduce(&my_count,&count,1,MPI_INT8,MPI_SUM,0,mpi_comm);
        if (mpi_rank == 0) 
          cout << " > fraction of cc dir's flagged: " << double(count)/double(3*ncc_global) << endl;
      }
    }
    
  }
 
  void AlgebraicCoarseGrid::update(const double* const A_cv,const int* const cvocv_i,const int* const cvocv_v,
      const uint8* const rbicv_g,const uint8* const rbicv_g2,const int ncv,const int ncv_g,const int ncv_g2) {
    assert(A_cv);
    assert(cvocv_i);
    assert(cvocv_v);
    assert(A_cc);

    for (int coc = 0; coc < ccocc_i[ncc_a]; ++coc)
      A_cc[coc] = 0.0;

    map<pair<uint8,uint8>,double> Amap; 
    map<pair<uint8,uint8>,double>::iterator it; 
    FOR_ICV {
      const int icc = ccocv[icv]; assert((icc >= 0)&&(icc < ncc));
      // we own A_cc row...
      if (icc < ncc_a) {
        for (int voc = cvocv_i[icv]; voc < cvocv_i[icv+1]; ++voc) {
          const int icv_nbr = cvocv_v[voc]; assert((icv_nbr >= 0)&&(icv_nbr < ncv_g2));
          const int icc_nbr = ccocv[icv_nbr]; assert((icc_nbr >= 0)&&(icc_nbr < ncc_g));
          // find coc...
          int coc = ccocc_i[icc];
          while (ccocc_v[coc] != icc_nbr)
            ++coc;
          assert(coc < ccocc_i[icc+1]);
          A_cc[coc] += A_cv[voc]; 
        }
      }
      else {
        const uint8 rbi = rbicc_i[icc-ncc_a];
        for (int voc = cvocv_i[icv]; voc < cvocv_i[icv+1]; ++voc) {
          const int icv_nbr = cvocv_v[voc]; assert((icv_nbr >= 0)&&(icv_nbr < ncv_g2));
          const int icc_nbr = ccocv[icv_nbr]; assert((icc_nbr >= 0)&&(icc_nbr < ncc_g));
          uint8 rbi_nbr;
          if (icc_nbr < ncc_a) 
            rbi_nbr = BitUtils::packRankBitsIndex(mpi_rank,0,icc_nbr);
          else if (icc_nbr < ncc)
            rbi_nbr = rbicc_i[icc_nbr-ncc_a];
          else 
            rbi_nbr = rbicc_g[icc_nbr-ncc];
          // find map entry...
          it = Amap.find(pair<uint8,uint8>(rbi,rbi_nbr));
          if (it == Amap.end()) {
            Amap[pair<uint8,uint8>(rbi,rbi_nbr)] = A_cv[voc];
          }
          else {
            it->second += A_cv[voc];
          }
        }
      }
    }

    int* send_count = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;
    int* send_disp = new int[mpi_size];
    int* send_buf_int = NULL;
    double* send_buf_double = NULL;
    for (int iter = 0; iter < 2; ++iter) {
      for (it = Amap.begin(); it != Amap.end(); ++it) {
        int rank,bits,index; BitUtils::unpackRankBitsIndex(rank,bits,index,it->first.first);
        if (iter == 0) {
          ++send_count[rank];
        }
        else {
          int rank_nbr,bits_nbr,index_nbr; BitUtils::unpackRankBitsIndex(rank_nbr,bits_nbr,index_nbr,it->first.second);
          send_buf_int[send_disp[rank]*3  ] = index; 
          send_buf_int[send_disp[rank]*3+1] = BitUtils::packRankBits(rank_nbr,bits_nbr);
          send_buf_int[send_disp[rank]*3+2] = index_nbr;
          send_buf_double[send_disp[rank]] = it->second; 
          ++send_disp[rank];
        }
      }

      // rewind...

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      if (iter == 0) {
        const int send_count_sum = send_disp[mpi_size-1]+send_count[mpi_size-1];
        send_buf_int = new int[3*send_count_sum];
        send_buf_double = new double[send_count_sum];
      }
    }
    Amap.clear();

    // setup recv side stuff...

    int* recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    int* recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    // exchange...
    
    double* recv_buf_double = new double[recv_count_sum];
    MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
        recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
    delete[] send_buf_double; 

    FOR_RANK {
      send_count[rank] *= 3;
      send_disp[rank] *= 3;
      recv_count[rank] *= 3;
      recv_disp[rank] *= 3;
    }

    int* recv_buf_int = new int[3*recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
        recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; 
    delete[] send_count;
    delete[] send_disp;
    delete[] recv_count;
    delete[] recv_disp;

    for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
      const int icc = recv_buf_int[irecv*3]; assert((icc >= 0)&&(icc < ncc_a));
      int rank_nbr,bits_nbr; BitUtils::unpackRankBits(rank_nbr,bits_nbr,recv_buf_int[irecv*3+1]);  
      const int icc_nbr = recv_buf_int[irecv*3+2];
      int coc = -1;
      if ((bits_nbr == 0)&&(rank_nbr == mpi_rank)) {
        for (coc = ccocc_i[icc]; coc < ccocc_i[icc+1]; ++coc) {
          const int icc_nbr_check = ccocc_v[coc]; assert((icc_nbr_check >= 0)&&(icc_nbr_check < ncc_g));
          if (icc_nbr_check == icc_nbr)
            break;
        }
      }
      else {
        const uint8 rbi_nbr = BitUtils::packRankBitsIndex(rank_nbr,bits_nbr,icc_nbr);
        for (coc = ccocc_i[icc]; coc < ccocc_i[icc+1]; ++coc) {
          const int icc_nbr_check = ccocc_v[coc]; assert((icc_nbr_check >= 0)&&(icc_nbr_check < ncc_g));
          if (icc_nbr_check >= ncc) {
            if (rbi_nbr == rbicc_g[icc_nbr_check-ncc])
              break;
          }
          else if (icc_nbr_check >= ncc_a) {
            if (rbi_nbr == rbicc_i[icc_nbr_check-ncc_a])
              break;
          }
        }
      }
      assert((coc >= ccocc_i[icc])&&(coc < ccocc_i[icc+1]));
      A_cc[coc] += recv_buf_double[irecv];
    }
    delete[] recv_buf_int;
    delete[] recv_buf_double;

    //for (int icc = 0; icc < ncc_a; ++icc) {
    //  cout << COUT_VEC(x_cc[icc]) << " : ";
    //  for (int coc = ccocc_i[icc]; coc != ccocc_i[icc+1]; ++coc) 
    //    cout << A_cc[coc] << " ";
    //  cout << endl;
    //}
    //cout << endl;
    //getchar();

    for (int icc = 0; icc < ncc_a; ++icc) {
      const int coc_f = ccocc_i[icc]; 
      double offd_sum = 0.0;
      for (int coc = coc_f+1; coc < ccocc_i[icc+1]; ++coc) {
        const int icc_nbr = ccocc_v[coc]; assert((icc_nbr >= 0)&&(icc_nbr < ncc_g));
        offd_sum += A_cc[coc];
        A_cc[coc] *= amg_coeff;
      }
      A_cc[coc_f] += (1.0-amg_coeff)*offd_sum;
    }

    for (int icc = 0; icc < ncc_a; ++icc)
      inv_diag[icc] = 1.0/A_cc[ccocc_i[icc]];
    updateCcData(inv_diag);

    //for (int icc = 0; icc < ncc_a; ++icc) {
    //  cout << COUT_VEC(x_cc[icc]) << " : ";
    //  for (int coc = ccocc_i[icc]; coc != ccocc_i[icc+1]; ++coc) 
    //    cout << A_cc[coc] << " ";
    //  cout << endl;
    //}
    //cout << endl;
    //getchar();
  }

  void AlgebraicCoarseGrid::calcCcResidual(double *res,const double* phi,const double *rhs) {

    for (int icc = 0; icc < ncc_a; ++icc) {
      const int coc_f = ccocc_i[icc];
      res[icc] = rhs[icc] - A_cc[coc_f]*phi[icc];
      for (int coc = coc_f+1; coc != ccocc_i[icc+1]; ++coc) {
        const int icc_nbr = ccocc_v[coc];
        res[icc] -= A_cc[coc]*phi[icc_nbr];
      }
    }

  }

  void AlgebraicCoarseGrid::smoothCcJacobi(double* phi,double *res,const double *rhs,const int nsmooth,const double relax,const double zero) {

    int iter = 0;
    int done = 0;
    while (done == 0) {
      ++iter;

      calcCcResidual(res,phi,rhs);

      for (int icc = 0; icc < ncc_a; ++icc) 
        phi[icc] += relax*res[icc]*inv_diag[icc];
      updateCcData(phi);

      if (zero > 0.0) {
        if (iter%3 == 0) {
          calcCcResidual(res,phi,rhs);
          double my_res_max = 0.0;
          for (int icc = 0; icc < ncc_a; ++icc) 
            my_res_max = max(my_res_max,fabs(res[icc]*inv_diag[icc]));
          double res_max;
          MPI_Reduce(&my_res_max,&res_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
          if (mpi_rank == 0) {
            if (res_max <= zero) {
              //cout << " > coarse grid solve converged after " << iter << " iters, res_max: " << res_max << endl;
              done = 1;
            }
            else if (iter >= nsmooth) {
              cout << " > Warning coarse grid solve did not converge after " << nsmooth << " iters, res_max: " << res_max << endl;
              done = 2;
            }
          }
          MPI_Bcast(&done,1,MPI_INT,0,mpi_comm);
        }
      }
      else if (iter == nsmooth) {
        done = 2;
      }
    }

  }

  void AlgebraicCoarseGrid::smoothCcSgs(double* phi,double* res,const double *rhs,const int nsmooth,const double relax,const double zero) {

    // note that this is only symmetric Gauss-Seidel locally

    int iter = 0;
    int done = 0;
    while (done == 0) {
      ++iter;

      for (int icc = 0; icc < ncc_a; ++icc) {
        const double phi0 = phi[icc];
        phi[icc] = rhs[icc];
        for (int coc = ccocc_i[icc]+1; coc != ccocc_i[icc+1]; ++coc) {
          const int icc_nbr = ccocc_v[coc];
          phi[icc] -= A_cc[coc]*phi[icc_nbr];
        }
        phi[icc] = phi0 + relax*(phi[icc]*inv_diag[icc]-phi0);
      }
      updateCcData(phi);
      for (int icc = ncc_a-1; icc >= 0; --icc) {
        const double phi0 = phi[icc];
        phi[icc] = rhs[icc];
        for (int coc = ccocc_i[icc]+1; coc != ccocc_i[icc+1]; ++coc) {
          const int icc_nbr = ccocc_v[coc];
          phi[icc] -= A_cc[coc]*phi[icc_nbr];
        }
        phi[icc] = phi0 + relax*(phi[icc]*inv_diag[icc]-phi0);
      }
      updateCcData(phi);

      if (zero > 0.0) {
        if (iter%3 == 0) {
          calcCcResidual(res,phi,rhs);
          double my_res_max = 0.0;
          for (int icc = 0; icc < ncc_a; ++icc) 
            my_res_max = max(my_res_max,fabs(res[icc]*inv_diag[icc]));
          double res_max;
          MPI_Reduce(&my_res_max,&res_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
          if (mpi_rank == 0) {
            if (res_max <= zero) {
              //cout << " > coarse grid solve converged after " << iter << " iters, res_max: " << res_max << endl;
              done = 1;
            }
            else if (iter >= nsmooth) {
              cout << " > Warning coarse grid solve did not converge after " << nsmooth << " iters, res_max: " << res_max << endl;
              done = 2;
            }
          }
          MPI_Bcast(&done,1,MPI_INT,0,mpi_comm);
        }
      }
      else if (iter == nsmooth) {
        done = 2;
      }

    }

  }

  void AlgebraicCoarseGrid::smoothCcGs(double* phi,double* res,const double *rhs,const int nsmooth,const double relax,const double zero) {

    // note that this is only Gauss-Seidel locally

    int iter = 0;
    int done = 0;
    while (done == 0) {
      ++iter;

      for (int icc = 0; icc < ncc_a; ++icc) {
        const double phi0 = phi[icc];
        phi[icc] = rhs[icc];
        for (int coc = ccocc_i[icc]+1; coc != ccocc_i[icc+1]; ++coc) {
          const int icc_nbr = ccocc_v[coc];
          phi[icc] -= A_cc[coc]*phi[icc_nbr];
        }
        phi[icc] = phi0 + relax*(phi[icc]*inv_diag[icc]-phi0);
      }
      updateCcData(phi);

      if (zero > 0.0) {
        if (iter%3 == 0) {
          calcCcResidual(res,phi,rhs);
          double my_res_max = 0.0;
          for (int icc = 0; icc < ncc_a; ++icc) 
            my_res_max = max(my_res_max,fabs(res[icc]*inv_diag[icc]));
          double res_max;
          MPI_Reduce(&my_res_max,&res_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
          if (mpi_rank == 0) {
            if (res_max <= zero) {
              //cout << " > coarse grid solve converged after " << iter << " iters, res_max: " << res_max << endl;
              done = 1;
            }
            else if (iter >= nsmooth) {
              cout << " > Warning coarse grid solve did not converge after " << nsmooth << " iters, res_max: " << res_max << endl;
              done = 2;
            }
          }
          MPI_Bcast(&done,1,MPI_INT,0,mpi_comm);
        }
      }
      else if (iter == nsmooth) {
        done = 2;
      }

    }

  }

  void AlgebraicCoarseGrid::smoothCcCg(double *phi,double *res,double* v,double *p,const double *rhs,const int nsmooth,const double zero) {

    // assume we come in with a consistent initial condition...

    for (int icc = 0; icc < ncc_a; ++icc)
      p[icc] = 0.0;
    double rho = 1.0;

    // calculate the residual in rhs format...
    calcCcResidual(res,phi,rhs);

    // diagonal precon/compute normalized residual...
    for (int icc = 0; icc < ncc_a; ++icc)
      v[icc] = res[icc]*inv_diag[icc];

    int iter = 0;
    int done = 0;
    while (done == 0) {
      ++iter;

      double rho_prev = rho;
      double my_rho = 0.0;
      for (int icc = 0; icc < ncc_a; ++icc)
        my_rho += res[icc]*v[icc];
      MPI_Allreduce(&my_rho,&rho,1,MPI_DOUBLE,MPI_SUM,mpi_comm);

      double beta = rho/rho_prev;
      for (int icc = 0; icc < ncc_a; ++icc)
        p[icc] = v[icc] + beta*p[icc];
      updateCcData(p);

      // v = [Ap]{p}...
      for (int icc = 0; icc < ncc_a; ++icc) {
        v[icc] = A_cc[ccocc_i[icc]]*p[icc];
        for (int coc = ccocc_i[icc]+1; coc < ccocc_i[icc+1]; ++coc) {
          const int icc_nbr = ccocc_v[coc];
          v[icc] += A_cc[coc]*p[icc_nbr];
        }
      }

      double my_gamma = 0.0;
      for (int icc = 0; icc < ncc_a; ++icc)
        my_gamma += p[icc]*v[icc];
      double gamma;
      MPI_Allreduce(&my_gamma,&gamma,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
      const double alpha = rho/gamma;

      // update full phi including ghosts...
      for (int icc = 0; icc < ncc_g; ++icc)
        phi[icc] += alpha*p[icc];

      for (int icc = 0; icc < ncc_a; ++icc) {
        // the unreduced residual...
        res[icc] -= alpha*v[icc];
        // still need to compute v, diag precon for next iteration...
        v[icc] = res[icc]*inv_diag[icc];
      }

      if (zero > 0.0) {
        if (iter%3 == 0) {
          double my_res_max = 0.0;
          for (int icc = 0; icc < ncc_a; ++icc) 
            my_res_max = max(my_res_max,fabs(v[icc]));
          double res_max;
          MPI_Reduce(&my_res_max,&res_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
          if (mpi_rank == 0) {
            if (res_max <= zero) {
              //cout << " > coarse grid solve converged after " << iter << " iters, res_max: " << res_max << endl;
              done = 1;
            }
            else if (iter >= nsmooth) {
              cout << " > Warning coarse grid solve did not converge after " << nsmooth << " iters, res_max: " << res_max << endl;
              done = 2;
            }
          }
          MPI_Bcast(&done,1,MPI_INT,0,mpi_comm);
        }
      }
      else if (iter == nsmooth) {
        done = 2;
      }

    }

  }

  void AlgebraicCoarseGrid::smoothCcGmres(double* phi,double *res,double * w,double *z,const double *rhs,const int nsmooth,const double zero) {

    const int m = nsmooth; 
    double * H = new double[m*m];
    double * v = new double[ncc_g*m];
    double * cs = new double[m+1];
    double * sn = new double[m+1];
    double * s = new double[m+1];
    double * y = new double[m+1];
    double * inv_sqrt_diag = new double[ncc_g];

    for (int icc = 0; icc < ncc_g; ++icc) {
      //assert(inv_diag[icc] < 0.0);
      inv_sqrt_diag[icc] = sqrt(fabs(inv_diag[icc]));
      //inv_sqrt_diag[icc] = 1.0;
    }

    // calculate the residual in rhs format...
    calcCcResidual(res,phi,rhs);

    int iter = 0;
    int done = 0;
    while (done == 0) {
      ++iter;

      for (int icc = 0; icc < ncc_a; ++icc)
        res[icc] *= inv_sqrt_diag[icc];
      //for (int icc = 0; icc < ncc_a; ++icc)
      //  res[icc] *= inv_diag[icc];
      double beta = sqrt(dot(res,res));

      for (int i = 0; i < m+1; ++i) { 
        cs[i] = 0.0;
        sn[i] = 0.0;
        s[i] = 0.0;
      }
      for (int i = 0; i < m*m ; ++i) 
        H[i] = 0.0;
      for (int i = 0; i < m*ncc_g; ++i) 
        v[i] = 0.0;
      for (int icc = 0; icc < ncc_a; ++icc) 
        v[(0)*ncc_g+icc] = res[icc]/beta;
      updateCcData(v);
      s[0] = beta;

      for (int i = 0; i < m-1; ++i) {

        for (int icc = 0; icc < ncc_g; ++icc)
          z[icc] = v[(i)*ncc_g+icc]*inv_sqrt_diag[icc];
        matvec(w,A_cc,z);
        for (int icc = 0; icc < ncc_a; ++icc)
          w[icc] *= inv_sqrt_diag[icc];
        //matvec(w,A_cc,&v[(i)*ncc_g]);
        //for (int icc = 0; icc < ncc_a; ++icc)
        //  w[icc] *= inv_diag[icc];

        // gram-schmidt...
        for (int k = 0; k <= i; ++k) { 
          H[(i)*m+k] = dot(w,&v[(k)*ncc_g]);
          for (int icc = 0; icc < ncc_a; ++icc) 
            w[icc] -= H[(i)*m+k] * v[(k)*ncc_g+icc];
        }

        H[(i)*m+i+1] = sqrt(dot(w,w));
        for (int icc = 0; icc < ncc_a; ++icc)  
          v[(i+1)*ncc_g+icc] = w[icc]/(H[(i)*m+i+1]);
        updateCcData(&v[(i+1)*ncc_g]);

        for (int k = 0; k < i; ++k) 
          applyPlaneRotation(H[(i)*m+k],H[(i)*m+k+1],cs[k],sn[k]);

        generatePlaneRotation(H[(i)*m+i],H[(i)*m+i+1],cs[i],sn[i]);
        applyPlaneRotation(H[(i)*m+i],H[(i)*m+i+1],cs[i],sn[i]);
        applyPlaneRotation(s[i],s[i+1],cs[i],sn[i]);
      }

      // update gmres...
      for (int i = 0; i < m+1; ++i) 
        y[i] = s[i];

      // back substitution... 
      for (int i = m-2; i >= 0; --i) { 
        y[i] /= H[(i)*m+i];
        for (int j = i-1; j >= 0; --j) 
          y[j] -= H[(i)*m+j] * y[i];
      }

      for (int icc = 0; icc < ncc_g; ++icc)
        z[icc] = 0.0;
      for (int j = 0; j <= m-2; ++j) 
        for (int icc = 0; icc < ncc_g; ++icc) 
          z[icc] += y[j] * v[(j)*ncc_g+icc];
      for (int icc = 0; icc < ncc_g; ++icc)
        phi[icc] += z[icc]*inv_sqrt_diag[icc];
      //for (int j = 0; j <= m-2; ++j) 
      //  for (int icc = 0; icc < ncc_g; ++icc) 
      //    phi[icc] += y[j] * v[(j)*ncc_g+icc];

      calcCcResidual(res,phi,rhs);
      
      if (zero > 0.0) {
        double my_res_max = 0.0;
        for (int icc = 0; icc < ncc_a; ++icc) 
          my_res_max = max(my_res_max,fabs(res[icc]*inv_diag[icc]));
        double res_max;
        MPI_Reduce(&my_res_max,&res_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
        if (mpi_rank == 0) {
          if (res_max <= zero) {
            //cout << " > coarse grid solve converged after " << iter << " iters, res_max: " << res_max << endl;
            done = 1;
          }
          else if (iter >= nsmooth) {
            cout << " > Warning coarse grid solve did not converge after " << nsmooth << " iters, res_max: " << res_max << endl;
            done = 2;
          }
        }
        MPI_Bcast(&done,1,MPI_INT,0,mpi_comm);
      }
      else {
        done = 2;
      }
      
    }

    delete[] H;
    delete[] v;
    delete[] cs;
    delete[] sn;
    delete[] s;
    delete[] y;
    delete[] inv_sqrt_diag;
  }

  void AlgebraicMultigrid::init(double* vol_cv,double (* x_cv)[3],double* A_cv,int* cvocv_i,int* cvocv_v,
      int8* icv_global,uint8* rbicv_g,uint8* rbicv_g2,const int ncv,const int ncv_g,const int ncv_g2,
      const int8 ncv_global,const int ncg,const double agglomeration_factor,const bool split_orphaned_colors,
      const double amg_coeff) {

    if (!(ncg >= 0)) {
      CERR("The number of coarse grids must be non-negative to use multigrid.");
    }

    if (mpi_rank == 0)
      cout << "AlgebraicMultigrid::init(), nlevel (coarse + base grids): " << ncg+1 << endl;

    assert(acg == NULL);
    this->ncg = ncg;
    acg = new AlgebraicCoarseGrid[ncg+1];

    // copy base grid into level 0
    {
      acg[0].icg = 0;
      acg[0].ncc_global = ncv_global;
      acg[0].ncc_a = ncv;
      acg[0].ncc = ncv_g;
      acg[0].ncc_g = ncv_g2;
      acg[0].ccocc_i = cvocv_i;
      acg[0].ccocc_v = cvocv_v;
      acg[0].A_cc = A_cv;
      acg[0].vol_cc = vol_cv;
      acg[0].x_cc = x_cv;
      acg[0].icc_global = icv_global;
      acg[0].rbicc_i = rbicv_g;
      acg[0].rbicc_g = rbicv_g2;
      acg[0].buildCcIPrcomm();
      acg[0].buildCcPrcomm();
      assert(acg[0].inv_vol == NULL);
      acg[0].inv_vol = new double[acg[0].ncc_g];
      for (int icc = 0; icc < acg[0].ncc_g; ++icc)
        acg[0].inv_vol[icc] = 1.0/acg[0].vol_cc[icc];
      assert(acg[0].inv_diag == NULL); 
      acg[0].inv_diag = new double[acg[0].ncc_g];
      for (int icc = 0; icc < acg[0].ncc_a; ++icc)
        acg[0].inv_diag[icc] = 1.0/acg[0].A_cc[acg[0].ccocc_i[icc]];
      acg[0].updateCcData(acg[0].inv_diag);
    }

    {
      // build gradient operator...

      assert(acg[0].ccocc_grad_coeff == NULL);
      acg[0].ccocc_grad_coeff = new double[acg[0].ccocc_i[acg[0].ncc_a]][3];
      const double sigma_epsilon = getDoubleParam("SIGMA_EPSILON",0.45);
      int (*flag_cc)[3] = new int[acg[0].ncc_a][3];
      for (int icc = 0; icc < acg[0].ncc_a; ++icc) {

        double A_grad[9];
        FOR_I3 FOR_J3 A_grad[j*3+i] = 0.0;
        const int coc_f = acg[0].ccocc_i[icc];
        for (int coc = coc_f+1; coc < acg[0].ccocc_i[icc+1]; ++coc) {
          const int inb = acg[0].ccocc_v[coc]; 
          const double dx[3] = DIFF(acg[0].x_cc[inb],acg[0].x_cc[icc]); // periodicity?
          const double dr = MAG(dx);
          //double wi = dr*dr*dr; // should be 2nd order on certain configurations
          double wi = dr*dr; // O(1) sigma
          //double wi = 1.0;
          assert(wi > 0.0);
          wi = 1.0/wi;
          FOR_I3 FOR_J3 {
            const double coeff = wi*dx[i]*dx[j];
            A_grad[j*3+i] += coeff;
            if (i != j)
              A_grad[i*3+j] += coeff;
          }
        }

        // use svd to invert A_grad to handle rank deficiency...

        double U[9], V[9], sigma[3];
        double svd_eps_tol;
        calcSvd33(U,V,sigma,A_grad,svd_eps_tol);

        // now supply the regularized inverse .. let V <--- V\Sigma^+

        double sigma_inv[3];
        FOR_I3 {
          assert( sigma[i] >= 0.0);
          if (sigma[i] < sigma_epsilon) {
            sigma_inv[i] = 0.0; // cant trust the gradient in this dir..
            flag_cc[icc][i] = 1;
          }
          else { 
            sigma_inv[i] = 1.0/sigma[i];
            flag_cc[icc][i] = 0;
          }
        }
        //cout << "SIGMA: " << COUT_VEC(sigma) << endl;
        FOR_I3 FOR_J3 V[j*3+i] *= sigma_inv[j];

        // finally supply the inverse for the gradient coefficients ..

        FOR_I3 acg[0].ccocc_grad_coeff[coc_f][i] = 0.0;
        for (int coc = coc_f+1; coc < acg[0].ccocc_i[icc+1]; ++coc) {
          const int inb = acg[0].ccocc_v[coc]; 
          const double dx[3] = DIFF(acg[0].x_cc[inb],acg[0].x_cc[icc]); // periodicity?
          const double dr = MAG(dx);
          //double wi = dr*dr*dr; // should be 2nd order on certain configurations
          double wi = dr*dr; // O(1) sigma
          //double wi = 1.0;
          assert(wi > 0.0);
          wi = 1.0/wi;
          double tmp[3];
          FOR_I3 tmp[i] = wi*(U[i*3+0]*dx[0] + U[i*3+1]*dx[1] + U[i*3+2]*dx[2]);
          FOR_I3 acg[0].ccocc_grad_coeff[coc][i] = V[0*3+i]*tmp[0] + V[1*3+i]*tmp[1] + V[2*3+i]*tmp[2];
          FOR_I3 acg[0].ccocc_grad_coeff[coc_f][i] -= acg[0].ccocc_grad_coeff[coc][i];
        }

      }

      {
        const double grad_check[3] = { 1.1234, -1.3243, 1.5321 }; // some order-1 gradient
        double * phi               = new double[acg[0].ncc_g];
        double (*grad_phi)[3]      = new double[acg[0].ncc_a][3];
        for (int icc = 0; icc < acg[0].ncc_a; ++icc) {
          phi[icc] = DOT_PRODUCT(acg[0].x_cc[icc],grad_check);
          //phi[icc] = 1.12;
        }
        acg[0].updateCcData(phi);

        acg[0].calcCcGrad(grad_phi,phi);
        delete[] phi;

        //double my_err = 0.0;
        for (int icc = 0; icc < acg[0].ncc_a; ++icc) {
          //my_err += DIST(grad_phi[icc],grad_check);
          //cout << COUT_VEC(acg[0].x_cc[icc]) << " " << COUT_VEC(grad_phi[icc]) << " ";
          FOR_I3 {
            grad_phi[icc][i] /= grad_check[i];
            grad_phi[icc][i] -= 1.0;
          }
          //cout << COUT_VEC(grad_check) << " " << endl;
        }
        //cout << " AVG_GRAD_ERROR: " << my_err/acg[0].ncc_a << endl;
        MiscUtils::dumpRange(grad_phi,acg[0].ncc_a,"grad err");
        int8 my_count = 0;
        for (int icc = 0; icc < acg[0].ncc_a; ++icc) {
          FOR_I3 {
            if (flag_cc[icc][i] == 1) {
              grad_phi[icc][i] = 0.0;
              ++my_count;
            }
          }
        }
        delete[] flag_cc;
        MiscUtils::dumpRange(grad_phi,acg[0].ncc_a,"grad err unflagged");
        delete[] grad_phi;
        int8 count;
        MPI_Reduce(&my_count,&count,1,MPI_INT8,MPI_SUM,0,mpi_comm);
        if (mpi_rank == 0) 
          cout << " > fraction of cc dir's flagged: " << double(count)/double(3*acg[0].ncc_global) << endl;
      }
    }

    for (int icg = 0; icg < ncg; ++icg) {
      acg[icg+1].init(acg[icg].vol_cc,acg[icg].x_cc,acg[icg].A_cc,acg[icg].ccocc_i,acg[icg].ccocc_v,acg[icg].icc_global,
          acg[icg].rbicc_i,acg[icg].rbicc_g,acg[icg].ncc_a,acg[icg].ncc,acg[icg].ncc_g,acg[icg].ncc_global,icg+1,
          agglomeration_factor,split_orphaned_colors,amg_coeff);
    }

  }

  void AlgebraicMultigrid::update(double* A_cv) {
    
    acg[0].A_cc = A_cv;
    for (int icc = 0; icc < acg[0].ncc_a; ++icc)
      acg[0].inv_diag[icc] = 1.0/acg[0].A_cc[acg[0].ccocc_i[icc]];
    acg[0].updateCcData(acg[0].inv_diag);

    for (int icg = 0; icg < ncg; ++icg) {
      acg[icg+1].update(acg[icg].A_cc,acg[icg].ccocc_i,acg[icg].ccocc_v,acg[icg].rbicc_i,acg[icg].rbicc_g,
          acg[icg].ncc_a,acg[icg].ncc,acg[icg].ncc_g);
    }

  }

  void AlgebraicMultigrid::solve(double* phi,const double* rhs,const double zero,const string cycle,const int maxcycle,
      const bool b_linear_prolongation,const int nsmooth,const int dnsmooth,const string smoother,const double relax,const int maxiter_coarsest,
      const string solver_coarsest,const bool b_verbose) {

    double* tmp = new double[acg[0].ncc_g];
    double* p = new double[acg[0].ncc_g];
    double* v = new double[acg[0].ncc_g];
    double (*grad_err)[3] = new double[acg[0].ncc][3];

    int iter = 0;
    int done = 0;
    while (done == 0) {
      iter++;

      //cout << "icg: 0" << endl;

      if (ncg > 0) {
        if (smoother == "SGS")
          acg[0].smoothCcSgs(phi,tmp,rhs,nsmooth,relax);
        else if (smoother == "GS")
          acg[0].smoothCcGs(phi,tmp,rhs,nsmooth,relax);
        else if (smoother == "JACOBI")
          acg[0].smoothCcJacobi(phi,tmp,rhs,nsmooth,relax);
        else if (smoother == "CG")
          acg[0].smoothCcCg(phi,tmp,p,v,rhs,nsmooth);
        else if (smoother == "GMRES")
          acg[0].smoothCcGmres(phi,tmp,p,v,rhs,nsmooth);
        else
          acg[0].smoothCcJacobi(phi,tmp,rhs,nsmooth,relax);
      }
      else {
        if (smoother == "SGS")
          acg[0].smoothCcSgs(phi,tmp,rhs,nsmooth,relax,zero);
        else if (smoother == "GS")
          acg[0].smoothCcGs(phi,tmp,rhs,nsmooth,relax,zero);
        else if (smoother == "JACOBI")
          acg[0].smoothCcJacobi(phi,tmp,rhs,nsmooth,relax,zero);
        else if (smoother == "CG")
          acg[0].smoothCcCg(phi,tmp,p,v,rhs,nsmooth,zero);
        else if (smoother == "GMRES")
          acg[0].smoothCcGmres(phi,tmp,p,v,rhs,nsmooth,zero);
        else
          acg[0].smoothCcJacobi(phi,tmp,rhs,nsmooth,relax,zero);
      }

      // compute residual 

      acg[0].calcCcResidual(tmp,phi,rhs);

      // check if done

      double my_res_max = 0.0;
      for (int icc = 0; icc < acg[0].ncc_a; ++icc)
        my_res_max = max(my_res_max,fabs(tmp[icc]*acg[0].inv_diag[icc]));
      double res_max;
      MPI_Reduce(&my_res_max,&res_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      if (mpi_rank == 0) {
        if (b_verbose) 
          cout << " > solve iter, res_max: " << iter << " " << res_max << endl;
        if (res_max < zero) {
          done = 1;
        }
        else if (iter > maxcycle) {
          cout << " > Warning: solve did not converge after " << maxcycle << " cycles, res_max: " << res_max << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done,1,MPI_INT,0,mpi_comm);

      if (done != 0)
        break;

      if (ncg == 0)
        continue;

      // restrict residual 

      //acg[1].restrictCcData(acg[1].res,tmp); 
      acg[1].restrictExtrinsicCcData(acg[1].res,tmp); 
      acg[1].updateCcData(acg[1].res);

      for (int icg = 1; icg < ncg; ++icg) {

        //cout << "icg: " << icg << endl;

        // smooth

        for (int icc = 0; icc < acg[icg].ncc_g; ++icc) 
          acg[icg].err[icc] = 0.0;

        const int this_nsmooth = nsmooth+icg*dnsmooth;
        if (smoother == "SGS")
          acg[icg].smoothCcSgs(acg[icg].err,tmp,acg[icg].res,this_nsmooth,relax); 
        else if (smoother == "GS")
          acg[icg].smoothCcGs(acg[icg].err,tmp,acg[icg].res,this_nsmooth,relax); 
        else if (smoother == "JACOBI")
          acg[icg].smoothCcJacobi(acg[icg].err,tmp,acg[icg].res,this_nsmooth,relax); 
        else if (smoother == "CG")
          acg[icg].smoothCcCg(acg[icg].err,tmp,p,v,acg[icg].res,this_nsmooth);
        else if (smoother == "GMRES")
          acg[icg].smoothCcGmres(acg[icg].err,tmp,p,v,acg[icg].res,this_nsmooth);
        else
          acg[icg].smoothCcJacobi(acg[icg].err,tmp,acg[icg].res,this_nsmooth,relax); 

        // compute residual 

        acg[icg].calcCcResidual(tmp,acg[icg].err,acg[icg].res);

        // restrict residual 

        //acg[icg+1].restrictCcData(acg[icg+1].res,tmp); 
        acg[icg+1].restrictExtrinsicCcData(acg[icg+1].res,tmp); 
        acg[icg+1].updateCcData(acg[icg+1].res);

      }

      for (int icc = 0; icc < acg[ncg].ncc_g; ++icc) 
        acg[ncg].err[icc] = 0.0;

      // use maxiter_coarsest to toggle b/w smooth and solve...

      //cout << "icg: " << ncg << endl;

      if (maxiter_coarsest > 0) {

        // solve on coarsest grid
        
        if (solver_coarsest == "SGS")
          acg[ncg].smoothCcSgs(acg[ncg].err,tmp,acg[ncg].res,maxiter_coarsest,relax,zero);
        else if (solver_coarsest == "GS")
          acg[ncg].smoothCcGs(acg[ncg].err,tmp,acg[ncg].res,maxiter_coarsest,relax,zero);
        else if (solver_coarsest == "JACOBI")
          acg[ncg].smoothCcJacobi(acg[ncg].err,tmp,acg[ncg].res,maxiter_coarsest,relax,zero);
        else if (solver_coarsest == "CG")
          acg[ncg].smoothCcCg(acg[ncg].err,tmp,p,v,acg[ncg].res,maxiter_coarsest,zero);
        else if (solver_coarsest == "GMRES")
          acg[ncg].smoothCcGmres(acg[ncg].err,tmp,p,v,acg[ncg].res,maxiter_coarsest,zero);
        else
          acg[ncg].smoothCcJacobi(acg[ncg].err,tmp,acg[ncg].res,maxiter_coarsest,relax,zero);

      }
      else {

        // smooth on coarsest grid
        
        const int this_nsmooth = nsmooth+ncg*dnsmooth;
        if (smoother == "SGS")
          acg[ncg].smoothCcSgs(acg[ncg].err,tmp,acg[ncg].res,this_nsmooth,relax);
        else if (smoother == "GS")
          acg[ncg].smoothCcGs(acg[ncg].err,tmp,acg[ncg].res,this_nsmooth,relax);
        else if (smoother == "JACOBI")
          acg[ncg].smoothCcJacobi(acg[ncg].err,tmp,acg[ncg].res,this_nsmooth,relax);
        else if (smoother == "CG")
          acg[ncg].smoothCcCg(acg[ncg].err,tmp,p,v,acg[ncg].res,this_nsmooth);
        else if (smoother == "GMRES")
          acg[ncg].smoothCcGmres(acg[ncg].err,tmp,p,v,acg[ncg].res,this_nsmooth);
        else
          acg[ncg].smoothCcJacobi(acg[ncg].err,tmp,acg[ncg].res,this_nsmooth,relax);

      }

      if (cycle == "F_CYCLE") {

        for (int jcg = ncg; jcg > 1; --jcg) {

          // prolong/smooth from ncg to jcg...

          for (int icg = ncg; icg >= jcg; --icg) {

            // prolong and update guess

            if (b_linear_prolongation) {
              acg[icg].calcCcGrad(grad_err,acg[icg].err);
              acg[icg].updateCcIDataReverse(grad_err);  
              acg[icg].prolongCcDataAndUpdateGuess(acg[icg-1].err,acg[icg-1].x_cc,acg[icg].err,grad_err);
              acg[icg-1].updateCcData(acg[icg-1].err);
            }
            else {
              //acg[icg].updateCcIDataReverse(acg[icg].err); 
              acg[icg].prolongCcDataAndUpdateGuess(acg[icg-1].err,acg[icg].err);
              acg[icg-1].updateCcData(acg[icg-1].err);
            }

            // smooth

            if (icg > jcg) {

              //cout << "icg: " << icg-1 << ", jcg: " << jcg << endl;

              const int this_nsmooth = nsmooth+(icg-1)*dnsmooth;
              if (smoother == "SGS")
                acg[icg-1].smoothCcSgs(acg[icg-1].err,tmp,acg[icg-1].res,this_nsmooth,relax);
              else if (smoother == "GS")
                acg[icg-1].smoothCcGs(acg[icg-1].err,tmp,acg[icg-1].res,this_nsmooth,relax);
              else if (smoother == "JACOBI")
                acg[icg-1].smoothCcJacobi(acg[icg-1].err,tmp,acg[icg-1].res,this_nsmooth,relax);
              else if (smoother == "CG")
                acg[icg-1].smoothCcCg(acg[icg-1].err,tmp,p,v,acg[icg-1].res,this_nsmooth);
              else if (smoother == "GMRES")
                acg[icg-1].smoothCcGmres(acg[icg-1].err,tmp,p,v,acg[icg-1].res,this_nsmooth);
              else
                acg[icg-1].smoothCcJacobi(acg[icg-1].err,tmp,acg[icg-1].res,this_nsmooth,relax);
            }

          }

          // restrict/smooth from jcg-1 to ncg

          for (int icg = jcg-1; icg < ncg; ++icg) {

            //cout << "icg: " << icg << ", jcg: " << jcg << endl;

            if (icg > jcg-1) {
              for (int icc = 0; icc < acg[icg].ncc_g; ++icc) 
                acg[icg].err[icc] = 0.0;
            }

            // smooth

            const int this_nsmooth = nsmooth+(icg-1)*dnsmooth;
            if (smoother == "SGS")
              acg[icg].smoothCcSgs(acg[icg].err,tmp,acg[icg].res,this_nsmooth,relax); 
            else if (smoother == "GS")
              acg[icg].smoothCcGs(acg[icg].err,tmp,acg[icg].res,this_nsmooth,relax); 
            else if (smoother == "JACOBI")
              acg[icg].smoothCcJacobi(acg[icg].err,tmp,acg[icg].res,this_nsmooth,relax); 
            else if (smoother == "CG")
              acg[icg].smoothCcCg(acg[icg].err,tmp,p,v,acg[icg].res,this_nsmooth);
            else if (smoother == "GMRES")
              acg[icg].smoothCcGmres(acg[icg].err,tmp,p,v,acg[icg].res,this_nsmooth);
            else
              acg[icg].smoothCcJacobi(acg[icg].err,tmp,acg[icg].res,this_nsmooth,relax); 

            // compute residual 

            acg[icg].calcCcResidual(tmp,acg[icg].err,acg[icg].res);

            // restrict residual 

            //acg[icg+1].restrictCcData(acg[icg+1].res,tmp); 
            acg[icg+1].restrictExtrinsicCcData(acg[icg+1].res,tmp); 
            acg[icg+1].updateCcData(acg[icg+1].res);

          }

          //cout << "icg: " << ncg << ", jcg: " << jcg << endl;

          for (int icc = 0; icc < acg[ncg].ncc_g; ++icc) 
            acg[ncg].err[icc] = 0.0;

          // use maxiter_coarsest to toggle b/w smooth and solve...

          if (maxiter_coarsest > 0) {

            // solve on coarsest grid

            if (solver_coarsest == "SGS")
              acg[ncg].smoothCcSgs(acg[ncg].err,tmp,acg[ncg].res,maxiter_coarsest,relax,zero);
            else if (solver_coarsest == "GS")
              acg[ncg].smoothCcGs(acg[ncg].err,tmp,acg[ncg].res,maxiter_coarsest,relax,zero);
            else if (solver_coarsest == "JACOBI")
              acg[ncg].smoothCcJacobi(acg[ncg].err,tmp,acg[ncg].res,maxiter_coarsest,relax,zero);
            else if (solver_coarsest == "CG")
              acg[ncg].smoothCcCg(acg[ncg].err,tmp,p,v,acg[ncg].res,maxiter_coarsest,zero);
            else if (solver_coarsest == "GMRES")
              acg[ncg].smoothCcGmres(acg[ncg].err,tmp,p,v,acg[ncg].res,maxiter_coarsest,zero);
            else
              acg[ncg].smoothCcJacobi(acg[ncg].err,tmp,acg[ncg].res,maxiter_coarsest,relax,zero);

          }
          else {

            // smooth on coarsest grid

            const int this_nsmooth = nsmooth+ncg*dnsmooth;
            if (smoother == "SGS")
              acg[ncg].smoothCcSgs(acg[ncg].err,tmp,acg[ncg].res,this_nsmooth,relax);
            else if (smoother == "GS")
              acg[ncg].smoothCcGs(acg[ncg].err,tmp,acg[ncg].res,this_nsmooth,relax);
            else if (smoother == "JACOBI")
              acg[ncg].smoothCcJacobi(acg[ncg].err,tmp,acg[ncg].res,this_nsmooth,relax);
            else if (smoother == "CG")
              acg[ncg].smoothCcCg(acg[ncg].err,tmp,p,v,acg[ncg].res,this_nsmooth);
            else if (smoother == "GMRES")
              acg[ncg].smoothCcGmres(acg[ncg].err,tmp,p,v,acg[ncg].res,this_nsmooth);
            else
              acg[ncg].smoothCcJacobi(acg[ncg].err,tmp,acg[ncg].res,this_nsmooth,relax);

          }

        }

      }

      for (int icg = ncg; icg > 1; --icg) {

        //cout << "icg: " << icg-1 << endl;

        // prolong and update guess

        if (b_linear_prolongation) {
          acg[icg].calcCcGrad(grad_err,acg[icg].err);
          acg[icg].updateCcIDataReverse(grad_err);  
          acg[icg].prolongCcDataAndUpdateGuess(acg[icg-1].err,acg[icg-1].x_cc,acg[icg].err,grad_err);
          acg[icg-1].updateCcData(acg[icg-1].err);
        }
        else {
          //acg[icg].updateCcIDataReverse(acg[icg].err);  // shouldn't this already happen in smooth/updateCcData???
          acg[icg].prolongCcDataAndUpdateGuess(acg[icg-1].err,acg[icg].err); 
          acg[icg-1].updateCcData(acg[icg-1].err);
        }

        // smooth

        const int this_nsmooth = nsmooth+(icg-1)*dnsmooth;
        if (smoother == "SGS")
          acg[icg-1].smoothCcSgs(acg[icg-1].err,tmp,acg[icg-1].res,this_nsmooth,relax);
        else if (smoother == "GS")
          acg[icg-1].smoothCcGs(acg[icg-1].err,tmp,acg[icg-1].res,this_nsmooth,relax);
        else if (smoother == "JACOBI")
          acg[icg-1].smoothCcJacobi(acg[icg-1].err,tmp,acg[icg-1].res,this_nsmooth,relax);
        else if (smoother == "CG")
          acg[icg-1].smoothCcCg(acg[icg-1].err,tmp,p,v,acg[icg-1].res,this_nsmooth);
        else if (smoother == "GMRES")
          acg[icg-1].smoothCcGmres(acg[icg-1].err,tmp,p,v,acg[icg-1].res,this_nsmooth);
        else
          acg[icg-1].smoothCcJacobi(acg[icg-1].err,tmp,acg[icg-1].res,this_nsmooth,relax);

      }

      // prolong and update guess

      if (b_linear_prolongation) {
        acg[1].calcCcGrad(grad_err,acg[1].err);
        acg[1].updateCcIDataReverse(grad_err);  
        acg[1].prolongCcDataAndUpdateGuess(phi,acg[0].x_cc,acg[1].err,grad_err);
        acg[0].updateCcData(phi);
      }
      else {
        //acg[1].updateCcIDataReverse(acg[1].err); // shouldn't this already happen in smooth/updateCcData???
        acg[1].prolongCcDataAndUpdateGuess(phi,acg[1].err);
        acg[0].updateCcData(phi);
      }

    }

    delete[] tmp;
    delete[] p;
    delete[] v;
    delete[] grad_err;

  }
  
}
