// =====================================================================
// these routines use the mpi_comm to exchange primitive data.
// We assume BOTH cvPrcomm and cv2Prcomm are built. These
// routines exchange ALL ghost data: i.e. no need to call updateCvData 
// AND updateCv2Data. Just call updateCv2Data.
// =====================================================================

void updateCv2DataReverseStart(T * s,const int action = MIN_DATA) {

  // here T is a class that should provide pack and unpack routines...
  
  // all cv2Prcomm's are non-local at this point -- no periodicity...
  // and perhaps it does not matter. We should use the buffer infrastructure
  // anyways, even with periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateCv2DataStart(T*): s is already mapped. updateCv2DataFinish required on s.");
  
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;

  // count requests and buf sizes. This is a way to 
  // step forward simultaneously through the cvPrcommVec and 
  // cv2PrcommVec, which should be monotonic in rank by 
  // construction...
  
  int rank_check = -1; // monotonicity check...
  //const int data_size = T::data_size();
  int unpack_size = 0;
  int pack_size = 0;
  int ii = 0;
  int ii2 = 0;
  int rank_count = 0;
  CvPrcomm * cvPrcomm;
  CvPrcomm * cv2Prcomm;
  int rank;
  while ((ii != cvPrcommVec.size())||(ii2 != cv2PrcommVec.size())) {
    
    cvPrcomm = NULL;
    cv2Prcomm = NULL;
    rank = -1;
    
    if (ii == cvPrcommVec.size()) {
      assert(ii2 != cv2PrcommVec.size());
      cv2Prcomm = &cv2PrcommVec[ii2++];
      rank = cv2Prcomm->getRank();
    }
    else if (ii2 == cv2PrcommVec.size()) {
      assert(ii != cvPrcommVec.size());
      cvPrcomm = &cvPrcommVec[ii++];
      rank = cvPrcomm->getRank();
    }
    else if (cvPrcommVec[ii].getRank() < cv2PrcommVec[ii2].getRank()) {
      cvPrcomm = &cvPrcommVec[ii++];
      rank = cvPrcomm->getRank();
    }
    else if (cvPrcommVec[ii].getRank() > cv2PrcommVec[ii2].getRank()) {
      cv2Prcomm = &cv2PrcommVec[ii2++];
      rank = cv2Prcomm->getRank();
    }
    else {
      cvPrcomm = &cvPrcommVec[ii++];
      rank = cvPrcomm->getRank();
      cv2Prcomm = &cv2PrcommVec[ii2++];
      assert(rank == cv2Prcomm->getRank());
    }
    assert(rank != -1);
    assert(rank > rank_check);
    rank_check = rank;
    assert(cvPrcomm || cv2Prcomm); // could be one or both
    assert(cv2Prcomm); // I think we insist this on build -- could simplify the above logic
    if (cvPrcomm) {
      unpack_size += cvPrcomm->unpack_size;
      pack_size += cvPrcomm->packVec.size();
    }
    if (cv2Prcomm) {
      unpack_size += cv2Prcomm->unpack_size;
      pack_size += cv2Prcomm->packVec.size();
    }
    ++rank_count;
  }
  
  assert(rank_count == cv2PrcommVec.size()); // cv2PrcommVec has at least the ranks from cvPrcommVec
  assert(unpack_size == (ncv_g2-ncv));
  
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
  int icv_g = ncv;
  int icv_g2 = ncv_g;
  rank_count = 0;
  while ((ii != cvPrcommVec.size())||(ii2 != cv2PrcommVec.size())) {

    cvPrcomm = NULL;
    cv2Prcomm = NULL;
    rank = -1;
    
    if (ii == cvPrcommVec.size()) {
      assert(ii2 != cv2PrcommVec.size());
      cv2Prcomm = &cv2PrcommVec[ii2++];
      rank = cv2Prcomm->getRank();
    }
    else if (ii2 == cv2PrcommVec.size()) {
      assert(ii != cvPrcommVec.size());
      cvPrcomm = &cvPrcommVec[ii++];
      rank = cvPrcomm->getRank();
    }
    else if (cvPrcommVec[ii].getRank() < cv2PrcommVec[ii2].getRank()) {
      cvPrcomm = &cvPrcommVec[ii++];
      rank = cvPrcomm->getRank();
    }
    else if (cvPrcommVec[ii].getRank() > cv2PrcommVec[ii2].getRank()) {
      cv2Prcomm = &cv2PrcommVec[ii2++];
      rank = cv2Prcomm->getRank();
    }
    else {
      cvPrcomm = &cvPrcommVec[ii++];
      rank = cvPrcomm->getRank();
      cv2Prcomm = &cv2PrcommVec[ii2++];
      assert(rank == cv2Prcomm->getRank());
    }
    assert(rank != -1);
    assert(cvPrcomm || cv2Prcomm); // could be one or both
    assert(cv2Prcomm); // see above

    // post irecv...

    const int unpack_size0 = unpack_size;
    const int pack_size0 = pack_size;
    
    // pack...
    
    if (cvPrcomm) {
      for (int i = 0, n = cvPrcomm->unpack_size; i < n; ++i) {
	mrs->UNPACK_BUF[unpack_size+i] = s[icv_g++];
      }
      unpack_size += cvPrcomm->unpack_size;
      pack_size += cvPrcomm->packVec.size();
    }
    
    if (cv2Prcomm) {
      for (int i = 0, n = cv2Prcomm->unpack_size; i < n; ++i) {
	mrs->UNPACK_BUF[unpack_size+i] = s[icv_g2++];
      }
      unpack_size += cv2Prcomm->unpack_size;
      pack_size += cv2Prcomm->packVec.size();
    }

    // send/recv...
    
    MPI_Irecv(mrs->PACK_BUF+pack_size0,pack_size-pack_size0,
	      MPI_T,rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[rank_count]));
    
    MPI_Issend(mrs->UNPACK_BUF+unpack_size0,unpack_size-unpack_size0,
	       MPI_T,rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[rank_count]));

    ++rank_count;
    
  }

  // make sure we ended up where we expected...
  
  assert(icv_g == ncv_g);
  assert(icv_g2 == ncv_g2);
  assert(rank_count == cv2PrcommVec.size()); 

  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
  
}  

void updateCv2DataReverseFinish(T * s,const int action = MIN_DATA) {
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateCv2DataFinish(T*): s is not mapped.");
  
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
   
  // now we wait for all messages to be sent and received...
    
  MPI_Waitall(cv2PrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);

  // as soon as we are all sent, we can clear the buffer...

  delete[] mrs->UNPACK_BUF; mrs->UNPACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;
    
  // and finally, wait for the recv's...
    
  MPI_Waitall(cv2PrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);

  // once all recv's available, we need to unpack...
  
  //const int data_size = T::data_size();
  int pack_size = 0;
  int ii = 0;
  int ii2 = 0;
  CvPrcomm * cvPrcomm;
  CvPrcomm * cv2Prcomm;
  int rank;
  while ((ii != cvPrcommVec.size())||(ii2 != cv2PrcommVec.size())) {
    
    cvPrcomm = NULL;
    cv2Prcomm = NULL;
    rank = -1;
    
    if (ii == cvPrcommVec.size()) {
      assert(ii2 != cv2PrcommVec.size());
      cv2Prcomm = &cv2PrcommVec[ii2++];
      rank = cv2Prcomm->getRank();
    }
    else if (ii2 == cv2PrcommVec.size()) {
      assert(ii != cvPrcommVec.size());
      cvPrcomm = &cvPrcommVec[ii++];
      rank = cvPrcomm->getRank();
    }
    else if (cvPrcommVec[ii].getRank() < cv2PrcommVec[ii2].getRank()) {
      cvPrcomm = &cvPrcommVec[ii++];
      rank = cvPrcomm->getRank();
    }
    else if (cvPrcommVec[ii].getRank() > cv2PrcommVec[ii2].getRank()) {
      cv2Prcomm = &cv2PrcommVec[ii2++];
      rank = cv2Prcomm->getRank();
    }
    else {
      cvPrcomm = &cvPrcommVec[ii++];
      rank = cvPrcomm->getRank();
      cv2Prcomm = &cv2PrcommVec[ii2++];
      assert(rank == cv2Prcomm->getRank());
    }
    assert(rank != -1);
    assert(cvPrcomm || cv2Prcomm); // could be one or both

    // unpack...
    
    switch (action) {
    case MIN_DATA:
      {
        if (cvPrcomm) {
          for (int i = 0; i < cvPrcomm->packVec.size(); ++i) {
            const int icv = cvPrcomm->packVec[i];
            assert((icv >= 0)&&(icv < ncv));
            s[icv] = min(s[icv],mrs->PACK_BUF[pack_size+i]);
          }
          pack_size += cvPrcomm->packVec.size();
        }
        if (cv2Prcomm) {
          for (int i = 0; i < cv2Prcomm->packVec.size(); ++i) {
            const int icv = cv2Prcomm->packVec[i];
            assert((icv >= 0)&&(icv < ncv));
            s[icv] = min(s[icv],mrs->PACK_BUF[pack_size+i]);
          }
          pack_size += cv2Prcomm->packVec.size();
        }
      }
      break;
    case ADD_DATA:
      {
        if (cvPrcomm) {
          for (int i = 0; i < cvPrcomm->packVec.size(); ++i) {
            const int icv = cvPrcomm->packVec[i];
            assert((icv >= 0)&&(icv < ncv));
            s[icv] += mrs->PACK_BUF[pack_size+i];
          }
          pack_size += cvPrcomm->packVec.size();
        }
        if (cv2Prcomm) {
          for (int i = 0; i < cv2Prcomm->packVec.size(); ++i) {
            const int icv = cv2Prcomm->packVec[i];
            assert((icv >= 0)&&(icv < ncv));
            s[icv] += mrs->PACK_BUF[pack_size+i];
          }
          pack_size += cv2Prcomm->packVec.size();
        }
      }
      break;
    default:
      assert(0);
    }
  }
  
  // cleanup...
  
  delete[] mrs->PACK_BUF;  mrs->PACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->recvRequestArray;  mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
    
}

// a blocking version of the complete update...

void updateCv2DataReverse(T * s,const int action = MIN_DATA) {
  updateCv2DataReverseStart(s,action);
  updateCv2DataReverseFinish(s,action);
}

