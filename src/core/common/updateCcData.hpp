// =====================================================================
// these routines use the mpi_comm to exchange primitive data.
// We assume BOTH ccIPrcomm and ccPrcomm are built. These
// routines exchange ALL inactive/ghost data
// =====================================================================

void updateCcDataStart(T * s,const int action = REPLACE_DATA) {

  // here T is a class that should provide pack and unpack routines...
  
  // all ccPrcomm's are non-local at this point -- no periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateCv2DataStart(T*): s is already mapped. updateCv2DataFinish required on s.");
  
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;

  // count requests and buf sizes. This is a way to 
  // step forward simultaneously through the ccIPrcommVec and 
  // ccPrcommVec, which should be monotonic in rank by 
  // construction...
  
  int rank_check = -1; // monotonicity check...
  int unpack_size = 0;
  int pack_size = 0;
  int ii = 0;
  int ii2 = 0;
  int rank_count = 0;
  Prcomm * ccIPrcomm;
  Prcomm * ccPrcomm;
  int rank;
  while ((ii != ccIPrcommVec.size())||(ii2 != ccPrcommVec.size())) {
    
    ccIPrcomm = NULL;
    ccPrcomm = NULL;
    rank = -1;
    
    if (ii == ccIPrcommVec.size()) {
      assert(ii2 != ccPrcommVec.size());
      ccPrcomm = &ccPrcommVec[ii2++];
      rank = ccPrcomm->getRank();
    }
    else if (ii2 == ccPrcommVec.size()) {
      assert(ii != ccIPrcommVec.size());
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
    }
    else if (ccIPrcommVec[ii].getRank() < ccPrcommVec[ii2].getRank()) {
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
    }
    else if (ccIPrcommVec[ii].getRank() > ccPrcommVec[ii2].getRank()) {
      ccPrcomm = &ccPrcommVec[ii2++];
      rank = ccPrcomm->getRank();
    }
    else {
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
      ccPrcomm = &ccPrcommVec[ii2++];
      assert(rank == ccPrcomm->getRank());
    }
    assert(rank != -1);
    assert(rank > rank_check);
    rank_check = rank;
    assert(ccIPrcomm || ccPrcomm); // could be one or both
    assert(ccPrcomm); // I think we insist this on build -- could simplify the above logic
    if (ccIPrcomm) {
      // note that we use the ccIPrcomm in reverse in this exchange (pack <-> unpack)
      unpack_size += ccIPrcomm->pack_size;
      pack_size += ccIPrcomm->unpackVec.size();
    }
    if (ccPrcomm) {
      unpack_size += ccPrcomm->unpack_size;
      pack_size += ccPrcomm->packVec.size();
    }
    ++rank_count;
  }
  
  assert(rank_count == ccPrcommVec.size()); // ccPrcommVec has at least the ranks from ccIPrcommVec
  assert(unpack_size == (ncc_g-ncc_a));
  
  // allocate request arrays...
  
  mrs->sendRequestArray = new MPI_Request[rank_count];
  mrs->recvRequestArray = new MPI_Request[rank_count];
  
  // and pack and send...
  
  assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
  
  unpack_size = 0;
  pack_size = 0;
  ii = 0;
  ii2 = 0;
  rank_count = 0;
  while ((ii != ccIPrcommVec.size())||(ii2 != ccPrcommVec.size())) {

    ccIPrcomm = NULL;
    ccPrcomm = NULL;
    rank = -1;
    
    if (ii == ccIPrcommVec.size()) {
      assert(ii2 != ccPrcommVec.size());
      ccPrcomm = &ccPrcommVec[ii2++];
      rank = ccPrcomm->getRank();
    }
    else if (ii2 == ccPrcommVec.size()) {
      assert(ii != ccIPrcommVec.size());
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
    }
    else if (ccIPrcommVec[ii].getRank() < ccPrcommVec[ii2].getRank()) {
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
    }
    else if (ccIPrcommVec[ii].getRank() > ccPrcommVec[ii2].getRank()) {
      ccPrcomm = &ccPrcommVec[ii2++];
      rank = ccPrcomm->getRank();
    }
    else {
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
      ccPrcomm = &ccPrcommVec[ii2++];
      assert(rank == ccPrcomm->getRank());
    }
    assert(rank != -1);
    assert(ccIPrcomm || ccPrcomm); // could be one or both
    assert(ccPrcomm); // see above

    // post irecv...

    const int unpack_size0 = unpack_size;
    const int pack_size0 = pack_size;
    
    // pack...
    
    if (ccIPrcomm) {
      for (int i = 0, n = ccIPrcomm->unpackVec.size(); i < n; ++i) {
	const int icc = ccIPrcomm->unpackVec[i];
	assert((icc >= 0)&&(icc < ncc_a));
	mrs->PACK_BUF[pack_size+i] = s[icc];
      }
      unpack_size += ccIPrcomm->pack_size;
      pack_size += ccIPrcomm->unpackVec.size();
    }
    
    if (ccPrcomm) {
      for (int i = 0, n = ccPrcomm->packVec.size(); i < n; ++i) {
	const int icc = ccPrcomm->packVec[i];
	assert((icc >= 0)&&(icc < ncc_a));
	mrs->PACK_BUF[pack_size+i] = s[icc];
      }
      unpack_size += ccPrcomm->unpack_size;
      pack_size += ccPrcomm->packVec.size();
    }

    // send/recv...
    
    MPI_Irecv(mrs->UNPACK_BUF+unpack_size0,unpack_size-unpack_size0,
	      MPI_T,rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[rank_count]));
    
    MPI_Issend(mrs->PACK_BUF+pack_size0,pack_size-pack_size0,
	       MPI_T,rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[rank_count]));

    ++rank_count;
    
  }

  assert(rank_count == ccPrcommVec.size()); 

  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
  
}  

void updateCcDataFinish(T * s,const int action = REPLACE_DATA) {
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateCcDataFinish(T*): s is not mapped.");
  
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
   
  // now we wait for all messages to be sent and received...
    
  MPI_Waitall(ccPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);

  // as soon as we are all sent, we can clear the buffer...

  delete[] mrs->PACK_BUF;  mrs->PACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;
    
  // and finally, wait for the recv's...
    
  MPI_Waitall(ccPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);

  // once all recv's available, we need to unpack...
  
  int unpack_size = 0;
  int icc = ncc_a;
  int icc_g = ncc;
  int ii = 0;
  int ii2 = 0;
  Prcomm * ccIPrcomm;
  Prcomm * ccPrcomm;
  int rank;
  while ((ii != ccIPrcommVec.size())||(ii2 != ccPrcommVec.size())) {
    
    ccPrcomm = NULL;
    ccIPrcomm = NULL;
    rank = -1;
    
    if (ii == ccIPrcommVec.size()) {
      assert(ii2 != ccPrcommVec.size());
      ccPrcomm = &ccPrcommVec[ii2++];
      rank = ccPrcomm->getRank();
    }
    else if (ii2 == ccPrcommVec.size()) {
      assert(ii != ccIPrcommVec.size());
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
    }
    else if (ccIPrcommVec[ii].getRank() < ccPrcommVec[ii2].getRank()) {
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
    }
    else if (ccIPrcommVec[ii].getRank() > ccPrcommVec[ii2].getRank()) {
      ccPrcomm = &ccPrcommVec[ii2++];
      rank = ccPrcomm->getRank();
    }
    else {
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
      ccPrcomm = &ccPrcommVec[ii2++];
      assert(rank == ccPrcomm->getRank());
    }
    assert(rank != -1);
    assert(ccPrcomm || ccIPrcomm); // could be one or both

    // unpack...
    
    switch (action) {
    case REPLACE_DATA:
      {
        if (ccIPrcomm) {
          // remember pack is unpack
          for (int i = 0; i < ccIPrcomm->pack_size; ++i)
            s[icc++] = mrs->UNPACK_BUF[unpack_size+i];
          unpack_size += ccIPrcomm->pack_size;
        }
        if (ccPrcomm) {
          for (int i = 0; i < ccPrcomm->unpack_size; ++i)
            s[icc_g++] = mrs->UNPACK_BUF[unpack_size+i];
          unpack_size += ccPrcomm->unpack_size;
        }
      }
      break;
    default:
      assert(0);
    }
  }
  
  // make sure we ended up where we expected...
  
  assert(icc == ncc);
  assert(icc_g == ncc_g);
  
  // cleanup...
  
  delete[] mrs->UNPACK_BUF; mrs->UNPACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->recvRequestArray;  mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
    
}

// a blocking version of the complete update...

void updateCcData(T * s,const int action = REPLACE_DATA) {
  updateCcDataStart(s,action);
  updateCcDataFinish(s,action);
}

void updateCcDataStart(T (*s)[3],const int action = REPLACE_ROTATE_DATA) {

  // here T is a class that should provide pack and unpack routines...
  
  // all ccPrcomm's are non-local at this point -- no periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateCv2DataStart(T*): s is already mapped. updateCv2DataFinish required on s.");
  
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;

  // count requests and buf sizes. This is a way to 
  // step forward simultaneously through the ccIPrcommVec and 
  // ccPrcommVec, which should be monotonic in rank by 
  // construction...
  
  int rank_check = -1; // monotonicity check...
  int unpack_size = 0;
  int pack_size = 0;
  int ii = 0;
  int ii2 = 0;
  int rank_count = 0;
  Prcomm * ccIPrcomm;
  Prcomm * ccPrcomm;
  int rank;
  while ((ii != ccIPrcommVec.size())||(ii2 != ccPrcommVec.size())) {
    
    ccIPrcomm = NULL;
    ccPrcomm = NULL;
    rank = -1;
    
    if (ii == ccIPrcommVec.size()) {
      assert(ii2 != ccPrcommVec.size());
      ccPrcomm = &ccPrcommVec[ii2++];
      rank = ccPrcomm->getRank();
    }
    else if (ii2 == ccPrcommVec.size()) {
      assert(ii != ccIPrcommVec.size());
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
    }
    else if (ccIPrcommVec[ii].getRank() < ccPrcommVec[ii2].getRank()) {
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
    }
    else if (ccIPrcommVec[ii].getRank() > ccPrcommVec[ii2].getRank()) {
      ccPrcomm = &ccPrcommVec[ii2++];
      rank = ccPrcomm->getRank();
    }
    else {
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
      ccPrcomm = &ccPrcommVec[ii2++];
      assert(rank == ccPrcomm->getRank());
    }
    assert(rank != -1);
    assert(rank > rank_check);
    rank_check = rank;
    assert(ccIPrcomm || ccPrcomm); // could be one or both
    assert(ccPrcomm); // I think we insist this on build -- could simplify the above logic
    if (ccIPrcomm) {
      // note that we use the ccIPrcomm in reverse in this exchange (pack <-> unpack)
      unpack_size += ccIPrcomm->pack_size*3;
      pack_size += ccIPrcomm->unpackVec.size()*3;
    }
    if (ccPrcomm) {
      unpack_size += ccPrcomm->unpack_size*3;
      pack_size += ccPrcomm->packVec.size()*3;
    }
    ++rank_count;
  }
  
  assert(rank_count == ccPrcommVec.size()); // ccPrcommVec has at least the ranks from ccIPrcommVec
  assert(unpack_size == 3*(ncc_g-ncc_a));
  
  // allocate request arrays...
  
  mrs->sendRequestArray = new MPI_Request[rank_count];
  mrs->recvRequestArray = new MPI_Request[rank_count];
  
  // and pack and send...
  
  assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
  
  unpack_size = 0;
  pack_size = 0;
  ii = 0;
  ii2 = 0;
  rank_count = 0;
  while ((ii != ccIPrcommVec.size())||(ii2 != ccPrcommVec.size())) {

    ccIPrcomm = NULL;
    ccPrcomm = NULL;
    rank = -1;
    
    if (ii == ccIPrcommVec.size()) {
      assert(ii2 != ccPrcommVec.size());
      ccPrcomm = &ccPrcommVec[ii2++];
      rank = ccPrcomm->getRank();
    }
    else if (ii2 == ccPrcommVec.size()) {
      assert(ii != ccIPrcommVec.size());
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
    }
    else if (ccIPrcommVec[ii].getRank() < ccPrcommVec[ii2].getRank()) {
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
    }
    else if (ccIPrcommVec[ii].getRank() > ccPrcommVec[ii2].getRank()) {
      ccPrcomm = &ccPrcommVec[ii2++];
      rank = ccPrcomm->getRank();
    }
    else {
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
      ccPrcomm = &ccPrcommVec[ii2++];
      assert(rank == ccPrcomm->getRank());
    }
    assert(rank != -1);
    assert(ccIPrcomm || ccPrcomm); // could be one or both
    assert(ccPrcomm); // see above

    // post irecv...

    const int unpack_size0 = unpack_size;
    const int pack_size0 = pack_size;
    
    // pack...
    // TODO bits/transforms
    
    if (ccIPrcomm) {
      for (int i = 0, n = ccIPrcomm->unpackVec.size(); i < n; ++i) {
	const int icc = ccIPrcomm->unpackVec[i];
	assert((icc >= 0)&&(icc < ncc_a));
	mrs->PACK_BUF[pack_size+3*i  ] = s[icc][0];
	mrs->PACK_BUF[pack_size+3*i+1] = s[icc][1];
	mrs->PACK_BUF[pack_size+3*i+2] = s[icc][2];
      }
      unpack_size += ccIPrcomm->pack_size*3;
      pack_size += ccIPrcomm->unpackVec.size()*3;
    }
    
    if (ccPrcomm) {
      for (int i = 0, n = ccPrcomm->packVec.size(); i < n; ++i) {
	const int icc = ccPrcomm->packVec[i];
	assert((icc >= 0)&&(icc < ncc_a));
	mrs->PACK_BUF[pack_size+3*i  ] = s[icc][0];
	mrs->PACK_BUF[pack_size+3*i+1] = s[icc][1];
	mrs->PACK_BUF[pack_size+3*i+2] = s[icc][2];
      }
      unpack_size += ccPrcomm->unpack_size*3;
      pack_size += ccPrcomm->packVec.size()*3;
    }

    // send/recv...
    
    MPI_Irecv(mrs->UNPACK_BUF+unpack_size0,unpack_size-unpack_size0,
	      MPI_T,rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[rank_count]));
    
    MPI_Issend(mrs->PACK_BUF+pack_size0,pack_size-pack_size0,
	       MPI_T,rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[rank_count]));

    ++rank_count;
    
  }

  assert(rank_count == ccPrcommVec.size()); 

  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
  
}  

void updateCcDataFinish(T (*s)[3],const int action = REPLACE_ROTATE_DATA) {
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateCcDataFinish(T*): s is not mapped.");
  
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
   
  // now we wait for all messages to be sent and received...
    
  MPI_Waitall(ccPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);

  // as soon as we are all sent, we can clear the buffer...

  delete[] mrs->PACK_BUF;  mrs->PACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;
    
  // and finally, wait for the recv's...
    
  MPI_Waitall(ccPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);

  // once all recv's available, we need to unpack...
  
  int unpack_size = 0;
  int icc = ncc_a;
  int icc_g = ncc;
  int ii = 0;
  int ii2 = 0;
  Prcomm * ccIPrcomm;
  Prcomm * ccPrcomm;
  int rank;
  while ((ii != ccIPrcommVec.size())||(ii2 != ccPrcommVec.size())) {
    
    ccPrcomm = NULL;
    ccIPrcomm = NULL;
    rank = -1;
    
    if (ii == ccIPrcommVec.size()) {
      assert(ii2 != ccPrcommVec.size());
      ccPrcomm = &ccPrcommVec[ii2++];
      rank = ccPrcomm->getRank();
    }
    else if (ii2 == ccPrcommVec.size()) {
      assert(ii != ccIPrcommVec.size());
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
    }
    else if (ccIPrcommVec[ii].getRank() < ccPrcommVec[ii2].getRank()) {
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
    }
    else if (ccIPrcommVec[ii].getRank() > ccPrcommVec[ii2].getRank()) {
      ccPrcomm = &ccPrcommVec[ii2++];
      rank = ccPrcomm->getRank();
    }
    else {
      ccIPrcomm = &ccIPrcommVec[ii++];
      rank = ccIPrcomm->getRank();
      ccPrcomm = &ccPrcommVec[ii2++];
      assert(rank == ccPrcomm->getRank());
    }
    assert(rank != -1);
    assert(ccPrcomm || ccIPrcomm); // could be one or both

    // unpack...
    
    if (ccIPrcomm) {
      // remember pack is unpack
      for (int i = 0; i < ccIPrcomm->pack_size; ++i) {
        s[icc][0] = mrs->UNPACK_BUF[unpack_size+3*i  ];
        s[icc][1] = mrs->UNPACK_BUF[unpack_size+3*i+1];
        s[icc][2] = mrs->UNPACK_BUF[unpack_size+3*i+2];
        ++icc;
      }
      unpack_size += ccIPrcomm->pack_size*3;
    }
    if (ccPrcomm) {
      for (int i = 0; i < ccPrcomm->unpack_size; ++i) {
        s[icc_g][0] = mrs->UNPACK_BUF[unpack_size+3*i  ];
        s[icc_g][1] = mrs->UNPACK_BUF[unpack_size+3*i+1];
        s[icc_g][2] = mrs->UNPACK_BUF[unpack_size+3*i+2];
        ++icc_g;
      }
      unpack_size += ccPrcomm->unpack_size*3;
    }

  }
  
  // make sure we ended up where we expected...
  
  assert(icc == ncc);
  assert(icc_g == ncc_g);
  
  // cleanup...
  
  delete[] mrs->UNPACK_BUF; mrs->UNPACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->recvRequestArray;  mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
    
}

// a blocking version of the complete update...

void updateCcData(T (*s)[3],const int action = REPLACE_ROTATE_DATA) {
  updateCcDataStart(s,action);
  updateCcDataFinish(s,action);
}
