{
    
  if (mpi_rank == 0)
    cout << " > compressing nodes";
    
  for (int isend = 0; isend < send_int8_count_sum; ++isend) {
    const int ino_local = send_int8_ino_local[isend];
    if (ino_local >= 0) {

      /*
	if (debug_cnn) {
	if (ino_global[ino_local] == ino_debug) {
	cout << "packing + ino_debug into send buf on rank: " << mpi_rank << endl;
	}
	else if (ino_global[ino_local] == -ino_debug-2) {
	cout << "packing - ino_debug into send buf on rank: " << mpi_rank << endl;
	}
	}
      */

      send_int8_buf[isend] = ino_global[ino_local];
      assert(send_int8_buf[isend] != -1); // may be -ve now
    }
    else {
      assert(ino_local == -1);
      send_int8_buf[isend] = -1;
    }
  }
    
  int done = 0;
  while (done != 1) {
      
    // assume done...
    int my_done = 1; 

    if (mpi_rank == 0) {
      cout << ".";
      cout.flush();
    }
      
    MPI_Alltoallv(send_int8_buf,send_int8_count,send_int8_disp,MPI_INT8,
		  recv_int8_buf,recv_int8_count,recv_int8_disp,MPI_INT8,
		  mpi_comm);
      
    for (int irecv = 0; irecv < recv_int8_count_sum; ++irecv) {
      const int ino_local = recv_int8_ino_local[irecv];
      if (ino_local >= 0) { // valid node
	
	// modify check/set/done - this worked but resulted in potentially too many iterations...
	/*
	  if (recv_int8_buf[irecv] != ino_global[ino_local]) {
	  my_done = 0;
	  if (recv_int8_buf[irecv] < ino_global[ino_local]) {
	  ino_global[ino_local] = recv_int8_buf[irecv];
	  }
	  else {
	  recv_int8_buf[irecv] = ino_global[ino_local];
	  }
	  }
	*/

	/*
	  if (debug_cnn) {
	  if (ino_global[ino_local] == ino_debug) {
	  cout << "comparing recv_int8_buf[irecv]: " << recv_int8_buf[irecv] << " to + ino_debug on rank: " << mpi_rank << endl;
	  }
	  else if (ino_global[ino_local] == -ino_debug-2) {
	  cout << "comparing recv_int8_buf[irecv]: " << recv_int8_buf[irecv] << " to - ino_debug on rank: " << mpi_rank << endl;
	  }
	  }
	*/
      
	if (recv_int8_buf[irecv] < ino_global[ino_local]) {
	  ino_global[ino_local] = recv_int8_buf[irecv];
	  //my_done = 0; // skip my_done on recv side -- should be OK?
	}
	
      }
      else {
	assert(ino_local == -1);
	assert(recv_int8_buf[irecv] == -1);
      }
    }

    // pack recv buf after changing EVERONE possible in above loop...

    for (int irecv = 0; irecv < recv_int8_count_sum; ++irecv) {
      const int ino_local = recv_int8_ino_local[irecv];
      if (ino_local >= 0) { // valid node
	recv_int8_buf[irecv] = ino_global[ino_local];
	assert(recv_int8_buf[irecv] != -1);
      }
    }
    
    MPI_Alltoallv(recv_int8_buf,recv_int8_count,recv_int8_disp,MPI_INT8,
		  send_int8_buf,send_int8_count,send_int8_disp,MPI_INT8,
		  mpi_comm);
      
    for (int isend = 0; isend < send_int8_count_sum; ++isend) {
      const int ino_local = send_int8_ino_local[isend];
      if (ino_local >= 0) {
	
	// modify check/set/done - see note above...
	/*
	  if (send_int8_buf[isend] != ino_global[ino_local]) {
	  my_done = 0;
	  if (send_int8_buf[isend] < ino_global[ino_local]) {
	  ino_global[ino_local] = send_int8_buf[isend];
	  }
	  else {
	  send_int8_buf[isend] = ino_global[ino_local];
	  }
	  }
	*/

	/*
	  if (debug_cnn) {
	  if (ino_global[ino_local] == ino_debug) {
	  cout << "comparing send_int8_buf[isend]: " << send_int8_buf[isend] << " to + ino_debug on rank: " << mpi_rank << endl;
	  }
	  else if (ino_global[ino_local] == -ino_debug-2) {
	  cout << "comparing send_int8_buf[isend]: " << send_int8_buf[isend] << " to - ino_debug on rank: " << mpi_rank << endl;
	  }
	  }
	*/
	
	if (send_int8_buf[isend] < ino_global[ino_local]) {
	  ino_global[ino_local] = send_int8_buf[isend];
	  my_done = 0;
	}
	
      }
      else {
	assert(ino_local == -1);
	assert(send_int8_buf[isend] == -1);
      }
    }

    MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);

    if (done != 1) {
      
      // only if we are not done, pack the send buf for next iteration...
      
      for (int isend = 0; isend < send_int8_count_sum; ++isend) {
	const int ino_local = send_int8_ino_local[isend];
	if (ino_local >= 0) {
	  send_int8_buf[isend] = ino_global[ino_local];
	  assert(send_int8_buf[isend] != -1);
	}
      }

    }
    
    /*
      if (debug_cnn)       
      MPI_Pause("\ndone cnn iter\n");
    */
      
  }
  
  if (mpi_rank == 0)
    cout << "OK" << endl;
    
}
