{
    
  if (mpi_rank == 0)
    cout << " > compressing no_bits";

  // outside the loop we set first, then pack second...
    
  for (int isend = 0; isend < send_int8_count_sum; ++isend) {
    const int ino_local = send_int8_ino_local[isend];
    if (ino_local >= 0) {
      assert(ino_local < nno_local_reduced);
      const int fa_bits = send_int8_faono_bits[isend];
      //assert((fa_bits == 0)||(fa_bits == 1)||(fa_bits == 4)||(fa_bits == 16)); // i.e. loser face nbrs have the first bit in the pair
      no_bits[ino_local] |= fa_bits;
    }
    else {
      assert(ino_local == -1);
      assert(send_int8_faono_bits[isend] == 0);
    }
  }

  for (int isend = 0; isend < send_int8_count_sum; ++isend) {
    const int ino_local = send_int8_ino_local[isend];
    if (ino_local >= 0) {
      assert(ino_local < nno_local_reduced);
      const int fa_bits = send_int8_faono_bits[isend];
      // send back with the bit (should be just one) in fa_bits cleared, and 
      // the inverse bit set...
      send_int_buf[isend] = ( (no_bits[ino_local] & (~fa_bits)) | BitUtils::flipPeriodicBits(fa_bits) );
      assert(send_int_buf[isend] != -1);
    }
    else {
      assert(ino_local == -1);
      assert(send_int8_faono_bits[isend] == 0);
      send_int_buf[isend] = -1;
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
      
    MPI_Alltoallv(send_int_buf,send_int8_count,send_int8_disp,MPI_INT,
		  recv_int_buf,recv_int8_count,recv_int8_disp,MPI_INT,
		  mpi_comm);
      
    for (int irecv = 0; irecv < recv_int8_count_sum; ++irecv) {
      const int ino_local = recv_int8_ino_local[irecv];
      if (ino_local >= 0) { // valid node
	assert(ino_local < nno_local_reduced);
	assert(recv_int_buf[irecv] >= 0);
	no_bits[ino_local] |= recv_int_buf[irecv];
	//my_done = 0; // skip my_done on recv side -- should be OK?
      }
      else {
	assert(ino_local == -1);
	assert(recv_int_buf[irecv] == -1);
      }
    }

    // pack recv buf after changing EVERONE possible in above loop...
    // note that recv side has no fa_bits. This is unnecessary. send
    // side info is here after just one iter... 
    for (int irecv = 0; irecv < recv_int8_count_sum; ++irecv) {
      const int ino_local = recv_int8_ino_local[irecv];
      if (ino_local >= 0) { // valid node
	assert(ino_local < nno_local_reduced);
	recv_int_buf[irecv] = no_bits[ino_local];
	assert(recv_int_buf[irecv] != -1);
      }
    }
    
    MPI_Alltoallv(recv_int_buf,recv_int8_count,recv_int8_disp,MPI_INT,
		  send_int_buf,send_int8_count,send_int8_disp,MPI_INT,
		  mpi_comm);
      
    for (int isend = 0; isend < send_int8_count_sum; ++isend) {
      const int ino_local = send_int8_ino_local[isend];
      if (ino_local >= 0) {
	assert(ino_local < nno_local_reduced);
	assert(send_int_buf[isend] >= 0);
	const int fa_bits = send_int8_faono_bits[isend];
	const int this_no_bits = ( (send_int_buf[isend] & ~BitUtils::flipPeriodicBits(fa_bits)) | fa_bits );
	if (this_no_bits != no_bits[ino_local]) {
	  no_bits[ino_local] |= this_no_bits;
	  my_done = 0;
	}
      }
      else {
	assert(ino_local == -1);
	assert(send_int_buf[isend] == -1);
	assert(send_int8_faono_bits[isend] == 0);
      }
    }
    
    MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);
    
    if (done != 1) {
      
      // only if we are not done, pack the send buf for next iteration...
      
      for (int isend = 0; isend < send_int8_count_sum; ++isend) {
	const int ino_local = send_int8_ino_local[isend];
	if (ino_local >= 0) {
	  assert(ino_local < nno_local_reduced);
	  const int fa_bits = send_int8_faono_bits[isend];
	  send_int_buf[isend] = ( (no_bits[ino_local] & (~fa_bits)) | BitUtils::flipPeriodicBits(fa_bits) );
	  assert(send_int_buf[isend] != -1);
	}
      }

    }
    
  }
  
  if (mpi_rank == 0)
    cout << "OK" << endl;
  
}
