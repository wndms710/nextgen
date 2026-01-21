// compressAndSend.hpp
{                        
                        
  // =======================================
  FOR_RANK send_count[rank] = 0;
  // we can significantly condense the rays before sending. Condense and count first...
  sort(windowVec.begin(),windowVec.end());
  for (int ilevel = 0; ilevel <= level_max; ++ilevel) 
    for (int iwindow = 0; iwindow <= nwi; ++iwindow) 
      level_count[(nwi+1)*ilevel+iwindow] = 0;
  int8 index_current = -1;
  for (int ii = 0, limit = windowVec.size(); ii < limit; ++ii) {
    const int8 index = windowVec[ii].index;
    if (index != index_current) {
      // when the index changes, the level should have been returned to zero...
      for (int ilevel = 0; ilevel <= level_max; ++ilevel) 
        for (int iwindow = 0; iwindow <= nwi; ++iwindow) 
          assert(level_count[(nwi+1)*ilevel+iwindow] == 0);
      index_current = index;
    }
    if (windowVec[ii].level < 0) {
      //const int il2 = -windowVec[ii].level-1; 
      //const int il2 = ((-windowVec[ii].level-1)&LEVEL_MASK_BITS);
      const int il2 = (nwi+1)*((-windowVec[ii].level-1)&LEVEL_MASK_BITS)+((-windowVec[ii].level-1)>>HCP_WINDOW_SHIFT_BITS);
      assert((il2 >= max(level,1))&&(il2 < (nwi+1)*(level_max+1)));
      --level_count[il2];
      if (level_count[il2] == 0) {
        // this is a keeper. It is the last "out" point
        const int rank = index%int8(mpi_size);
        ++send_count[rank];
      }
      else {
        // discard this one...
        assert(level_count[il2] > 0);
        windowVec[ii].level = 0; // this zero will be used to skip this entry during pack
      }
    }
    else {
      assert(windowVec[ii].level > 0);
      //const int il2 = windowVec[ii].level-1; 
      //assert((il2 >= max(level,1))&&(il2 <= level_max));
      const int il2 = (nwi+1)*((windowVec[ii].level-1)&LEVEL_MASK_BITS)+((windowVec[ii].level-1)>>HCP_WINDOW_SHIFT_BITS);
      assert((il2 >= max(level,1))&&(il2 < (nwi+1)*(level_max+1)));
      ++level_count[il2];
      if (level_count[il2] == 1) {
        // this is the first "in" point: keep...
        const int rank = index%int8(mpi_size);
        ++send_count[rank];
      }
      else {
        assert(level_count[il2] > 1);
        windowVec[ii].level = 0; // this zero will be used to skip this entry during pack
      }
    }
  }
  // level_count returned to zero...
  for (int ilevel = 0; ilevel <= level_max; ++ilevel) 
    for (int iwindow = 0; iwindow <= nwi; ++iwindow) 
      assert(level_count[(nwi+1)*ilevel+iwindow] == 0);
  // get ready to send...
  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
  const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
  int8 * send_buf_int8 = new int8[send_count_sum*2]; // for ints we need the index AND the level
  double * send_buf_double = new double[send_count_sum];
  for (int ii = 0, limit = windowVec.size(); ii < limit; ++ii) {
    if (windowVec[ii].level != 0) {
      const int8 index = windowVec[ii].index;
      const int rank = index%int8(mpi_size);
      send_buf_int8[send_disp[rank]*2]   = index;
      send_buf_int8[send_disp[rank]*2+1] = windowVec[ii].level; // the +/- level+1
      send_buf_double[send_disp[rank]]  = windowVec[ii].xp;
      ++send_disp[rank];
    }
  }
  // no longer reqd...
  windowVec.clear();
  // rewind...
  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
  // set up recv-side stuff...
  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
  const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
  assert(recv_buf_double == NULL);
  recv_buf_double = new double[recv_count_sum];
  MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
                recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
  delete[] send_buf_double;
  FOR_RANK {
    send_count[rank] *= 2;
    send_disp[rank] *= 2;
    recv_count[rank] *= 2;
    recv_disp[rank] *= 2;
  }
  int8 * recv_buf_int8 = new int8[recv_count_sum*2];
  MPI_Alltoallv(send_buf_int8,send_count,send_disp,MPI_INT8,
                recv_buf_int8,recv_count,recv_disp,MPI_INT8,mpi_comm);
  delete[] send_buf_int8;
  // push into the ray storage and then compress it (this will be
  // the maximum compression because we are completely responsible for
  // this ray's info).
  // here we can also confirm that every ray we recv'd has 
  // some level 0 events (ray_buf working)...
  for (int ii = 0; ii < my_ray_flag_size; ++ii)
    my_ray_flag[ii] = 0;
  int jj_prev = -1;
  for (int ii = 0, limit = inoutVec[level].size(); ii < limit; ++ii) {
    const int8 index = inoutVec[level][ii].index;
    assert((index-int8(mpi_rank))%int8(mpi_size) == 0ll);
    const int jj = (index-int8(mpi_rank))/int8(mpi_size);
    assert((jj >= 0)&&(jj < my_ray_flag_size));
    assert(jj >= jj_prev); // check sort
    jj_prev = jj;
    ++my_ray_flag[jj];
  }
  for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
    const int8 index = recv_buf_int8[irecv*2];
    assert(index >= 0);
    assert((index-int8(mpi_rank))%int8(mpi_size) == 0ll);
    const int jj = (index-int8(mpi_rank))/int8(mpi_size);
    assert((jj >= 0)&&(jj < my_ray_flag_size));
    if (my_ray_flag[jj] > 0) {
      // count only those rays that have existing surface intersections
      // and/or other entries...
      ++my_ray_flag[jj];
    }
    else {
      // this is refinement window information that is not
      // required because the ray does not pass through the 
      // geometry. Set the index to -1...
      recv_buf_int8[irecv*2] = -1;
    }
  }
  // my_ray_flag now has the counts associated with each level
  // in it. Turn it into the jj offset in each level...
  for (int ii = 1; ii < my_ray_flag_size; ++ii)
    my_ray_flag[ii] += my_ray_flag[ii-1];
  // now go backwards through the inoutVec[level], making room
  // for the new info...
  const int8 old_size = inoutVec[level].size();
  inoutVec[level].resize(my_ray_flag[my_ray_flag_size-1]);
  for (int ii = old_size-1; ii >= 0; --ii) {
    const int8 index = inoutVec[level][ii].index;
    const int jj = (index-int8(mpi_rank))/int8(mpi_size);
    --my_ray_flag[jj];
    const int ii_new = my_ray_flag[jj];
    assert(ii_new >= ii);
    inoutVec[level][ii_new] = inoutVec[level][ii];
  }
  for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
    const int8 index = recv_buf_int8[irecv*2];
    if (index >= 0) {
      assert((index-int8(mpi_rank))%int8(mpi_size) == 0ll);
      const int jj = (index-int8(mpi_rank))/int8(mpi_size);
      --my_ray_flag[jj];
      const int ii_new = my_ray_flag[jj];
      inoutVec[level][ii_new].index = index;
      inoutVec[level][ii_new].level = recv_buf_int8[irecv*2+1];
      inoutVec[level][ii_new].xp = recv_buf_double[irecv];
    }
  }
  // finally sort locally amongst each index and eliminate inout points
  // that provide no new information (i.e. in for a particular level when 
  // we are already in, or out before the final out)...
  vector<InOutData>::iterator iter_end = inoutVec[level].begin();
  while (iter_end != inoutVec[level].end()) {
    vector<InOutData>::iterator iter_begin = iter_end++;
    while ((iter_end != inoutVec[level].end())&&(iter_end->index == iter_begin->index))
      ++iter_end;
    sort(iter_begin,iter_end);
    for (int ilevel = 0; ilevel <= level_max; ++ilevel) 
      for (int iwindow = 0; iwindow <= nwi; ++iwindow) 
        level_count[(nwi+1)*ilevel+iwindow] = 0;

    //for (vector<InOutData>::iterator iter2 = iter_begin; iter2 != iter_end; ++iter2)
    //  cout << iter2->level << " ";
    //cout << endl;
    
    for (vector<InOutData>::iterator iter = iter_begin; iter != iter_end; ++iter) {
      if (iter->level > 0) {
        //const int il2 = iter->level-1;
        //assert((il2 >= 0)&&(il2 <= level_max));
        const int il2 = (nwi+1)*((iter->level-1)&LEVEL_MASK_BITS)+((iter->level-1)>>HCP_WINDOW_SHIFT_BITS);
        assert((il2 >= 0)&&(il2 < (nwi+1)*(level_max+1)));
        ++level_count[il2];
        /*
        // going from 0 to 1 is a keeper. anything larger, eliminate...
        if (level_count[il2] > 1) {
          assert(il2 > 0); // level 0 has already been handled
          iter->level = 0;
        }
        else {
          if (!(level_count[il2] == 1)) 
            cout << "il2: " << il2 << " level_count[il2]: " << level_count[il2] << " " << iter->level << endl;
          assert(level_count[il2] == 1);
        }
        */
        // going from 0 to 1 is a keeper. anything else, eliminate...
        if (level_count[il2] != 1) {
          iter->level = 0;
          if (level_count[il2] > 1)
            assert(il2 > 0); // level 0 has already been handled
        }
      }
      else {
        assert(iter->level < 0);
        //const int il2 = -iter->level-1;
        //assert((il2 >= 0)&&(il2 <= level_max));
        const int il2 = (nwi+1)*((-iter->level-1)&LEVEL_MASK_BITS)+((-iter->level-1)>>HCP_WINDOW_SHIFT_BITS);
        assert((il2 >= 0)&&(il2 < (nwi+1)*(level_max+1)));
        --level_count[il2];
        /*
        // going from 1 to 0 is a keeper. anything larger, eliminate...
        if (level_count[il2] > 0) {
          assert(il2 > 0); // level 0 has already been handled
          iter->level = 0;
        }
        else {
          if (!(level_count[il2] == 0)) 
            cout << "il2: " << il2 << " level_count[il2]: " << level_count[il2] << endl;
          assert(level_count[il2] == 0); 
        }
        */
        // going from 1 to 0 is a keeper. anything else, eliminate...
        if (level_count[il2] != 0) {
          iter->level = 0;
          if (level_count[il2] > 0)
            assert(il2 > 0); // level 0 has already been handled
        }
      }
    }
    for (int ilevel = 0; ilevel <= level_max; ++ilevel) 
      for (int iwindow = 0; iwindow <= nwi; ++iwindow) 
        assert(level_count[(nwi+1)*ilevel+iwindow] == 0);
  }
  // now all the zeros should be removed...
  int new_size = 0;
  for (int ii = 0, limit = inoutVec[level].size(); ii < limit; ++ii) {
    if (inoutVec[level][ii].level != 0) {
      inoutVec[level][new_size] = inoutVec[level][ii];
      ++new_size;
    }
  }
  inoutVec[level].resize(new_size);
  // and cleanup...
  delete[] recv_buf_int8; recv_buf_int8 = NULL; 
  delete[] recv_buf_double; recv_buf_double = NULL;
  // =======================================
                        

}
