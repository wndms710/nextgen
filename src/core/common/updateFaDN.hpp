// =====================================================================
// these routines use the mpi_comm to exchange data... 
// =====================================================================

void updateFaDataStart(T * s,const int action) {

  // all faPrcomm's are non-local at this point -- no periodicity...
  // and perhaps it does not matter. We should use the buffer infrastructure
  // anyways, even with periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateFaDataStart(T*): s is already mapped. updateFaDataFinish required on s.");
  
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;
  
  mrs->sendRequestArray = new MPI_Request[faPrcommVec.size()];
  mrs->recvRequestArray = new MPI_Request[faPrcommVec.size()];
  //mrs->statusArray = new MPI_Status[faPrcommVec.size()];
    
  int unpack_size = 0;
  int pack_size = 0;
  for (int ii = 0; ii < faPrcommVec.size(); ++ii) {
    unpack_size += faPrcommVec[ii].unpackVec.size();
    pack_size += faPrcommVec[ii].packVec.size();
  }
  
  // and pack and send...
  
  assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
  
  unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < faPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv(mrs->UNPACK_BUF+unpack_size,faPrcommVec[ii].unpackVec.size(),
	      MPI_T,faPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
    
    unpack_size += faPrcommVec[ii].unpackVec.size();
    
    // pack and send...
    
    for (int i = 0; i < faPrcommVec[ii].packVec.size(); ++i) {
      const int ifa = faPrcommVec[ii].packVec[i];
      assert((ifa >= nfa_i)&&(ifa < nfa));
      mrs->PACK_BUF[pack_size+i] = s[ifa];
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,faPrcommVec[ii].packVec.size(),
	       MPI_T,faPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += faPrcommVec[ii].packVec.size();

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateFaDataFinish(T * s,const  int action) {

  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateFaDataFinish(T*): s is not mapped.");
  
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
  
  // now we wait for all messages to be sent and received...
  
  MPI_Waitall(faPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);

  // as soon as we are all sent, we can clear the buffer...
  
  delete[] mrs->PACK_BUF;         mrs->PACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;

    // and finally, wait for the recv's...
    
  MPI_Waitall(faPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);
  
  // once all recv's available, we need to unpack...
  
  switch (action) {
  case MIN_NO_PERIODIC_DATA:
    {
      int unpack_size = 0;
      for (int ii = 0; ii < faPrcommVec.size(); ++ii) {
	for (int i = 0; i < faPrcommVec[ii].unpackVec.size(); ++i) {
	  const int ifa = faPrcommVec[ii].unpackVec[i];
	  assert((ifa >= nfa_i)&&(ifa < nfa));
          const int icv1 = cvofa[ifa][1];
          assert((icv1 >= ncv)&&(icv1 < ncv_g));
          int rank,bits,index;
          BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv1-ncv]);
          if (bits == 0)
            s[ifa] = min(s[ifa],mrs->UNPACK_BUF[unpack_size+i]);
	}
	unpack_size += faPrcommVec[ii].unpackVec.size();
      }
    }
    break;
  case MAX_NO_PERIODIC_DATA:
    {
      int unpack_size = 0;
      for (int ii = 0; ii < faPrcommVec.size(); ++ii) {
	for (int i = 0; i < faPrcommVec[ii].unpackVec.size(); ++i) {
	  const int ifa = faPrcommVec[ii].unpackVec[i];
	  assert((ifa >= nfa_i)&&(ifa < nfa));
          const int icv1 = cvofa[ifa][1];
          assert((icv1 >= ncv)&&(icv1 < ncv_g));
          int rank,bits,index;
          BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv1-ncv]);
          if (bits == 0)
            s[ifa] = max(s[ifa],mrs->UNPACK_BUF[unpack_size+i]);
	}
	unpack_size += faPrcommVec[ii].unpackVec.size();
      }
    }
    break;
  case ADD_DATA:
    {
      int unpack_size = 0;
      for (int ii = 0; ii < faPrcommVec.size(); ++ii) {
	for (int i = 0; i < faPrcommVec[ii].unpackVec.size(); ++i) {
	  const int ifa = faPrcommVec[ii].unpackVec[i];
	  assert((ifa >= nfa_i)&&(ifa < nfa));
	  s[ifa] += mrs->UNPACK_BUF[unpack_size+i];
	}
	unpack_size += faPrcommVec[ii].unpackVec.size();
      }
    }
    break;
  case SUBTRACT_DATA:
    {
      int unpack_size = 0;
      for (int ii = 0; ii < faPrcommVec.size(); ++ii) {
	for (int i = 0; i < faPrcommVec[ii].unpackVec.size(); ++i) {
	  const int ifa = faPrcommVec[ii].unpackVec[i];
	  assert((ifa >= nfa_i)&&(ifa < nfa));
	  s[ifa] -= mrs->UNPACK_BUF[unpack_size+i];
	}
	unpack_size += faPrcommVec[ii].unpackVec.size();
      }
    }
    break;
  case REPLACE_DATA:
    {
      int unpack_size = 0;
      for (int ii = 0; ii < faPrcommVec.size(); ++ii) {
	for (int i = 0; i < faPrcommVec[ii].unpackVec.size(); ++i) {
	  const int ifa = faPrcommVec[ii].unpackVec[i];
	  assert((ifa >= nfa_i)&&(ifa < nfa));
	  s[ifa] = mrs->UNPACK_BUF[unpack_size+i];
	}
	unpack_size += faPrcommVec[ii].unpackVec.size();
      }
    }
    break;
  default:
    assert(0);
  }
  
  // cleanup...
  
  delete[] mrs->UNPACK_BUF; mrs->UNPACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
  
}

void updateFaData(T * s,const int action) {
  updateFaDataStart(s,action);
  updateFaDataFinish(s,action);
}

