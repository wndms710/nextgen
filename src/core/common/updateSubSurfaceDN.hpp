// =====================================================================
// these routines use the mpi_comm to exchange data... 
// =====================================================================

void updateSpDataStart(T * s, const UpdateAction action = ADD_DATA) {

  // all spPrcomm's are non-local at this point -- no periodicity...
  // and perhaps it does not matter. We should use the buffer infrastructure
  // anyways, even with periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateSpDataStart(T*): s is already mapped. updateSpDataFinish required on s.");
  
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;
  
  mrs->sendRequestArray = new MPI_Request[spPrcommVec.size()];
  mrs->recvRequestArray = new MPI_Request[spPrcommVec.size()];
  //mrs->statusArray = new MPI_Status[spPrcommVec.size()];
    
  int unpack_size = 0;
  int pack_size = 0;
  for (int ii = 0; ii < spPrcommVec.size(); ++ii) {
    unpack_size += spPrcommVec[ii].unpackVec.size();
    pack_size += spPrcommVec[ii].packVec.size();
  }
  
  // and pack and send...
  
  assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
  
  unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < spPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv(mrs->UNPACK_BUF+unpack_size,spPrcommVec[ii].unpackVec.size(),
	      MPI_T,spPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
    
    unpack_size += spPrcommVec[ii].unpackVec.size();
    
    // pack and send...
    
    for (int i = 0; i < spPrcommVec[ii].packVec.size(); ++i) {
      const int isp = spPrcommVec[ii].packVec[i];
      assert((isp >= 0)&&(isp < nsp));
      mrs->PACK_BUF[pack_size+i] = s[isp];
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,spPrcommVec[ii].packVec.size(),
	       MPI_T,spPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += spPrcommVec[ii].packVec.size();

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateSpDataFinish(T * s, const UpdateAction action = ADD_DATA) {

  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateSpDataFinish(T*): s is not mapped.");
  
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
  
  // now we wait for all messages to be sent and received...
  
  MPI_Waitall(spPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);

  // as soon as we are all sent, we can clear the buffer...
  
  delete[] mrs->PACK_BUF;         mrs->PACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;

    // and finally, wait for the recv's...
    
  MPI_Waitall(spPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);
  
  // once all recv's available, we need to unpack...
  
  int unpack_size = 0;
  switch (action) {
  case ADD_DATA:
    for (int ii = 0; ii < spPrcommVec.size(); ++ii) {
      
      for (int i = 0; i < spPrcommVec[ii].unpackVec.size(); ++i) {
        const int isp = spPrcommVec[ii].unpackVec[i];
        assert((isp >= 0)&&(isp < nsp));
        s[isp] += mrs->UNPACK_BUF[unpack_size+i]; // default for nodal update is addition
      }
      
      unpack_size += spPrcommVec[ii].unpackVec.size();

    }
    break;
  case MIN_DATA:
    for (int ii = 0; ii < spPrcommVec.size(); ++ii) {
      
      for (int i = 0; i < spPrcommVec[ii].unpackVec.size(); ++i) {
        const int isp = spPrcommVec[ii].unpackVec[i];
        assert((isp >= 0)&&(isp < nsp));
        s[isp] = min(s[isp],mrs->UNPACK_BUF[unpack_size+i]);
      }
      
      unpack_size += spPrcommVec[ii].unpackVec.size();

    }
    break;
  case MAX_DATA:
    for (int ii = 0; ii < spPrcommVec.size(); ++ii) {
      
      for (int i = 0; i < spPrcommVec[ii].unpackVec.size(); ++i) {
        const int isp = spPrcommVec[ii].unpackVec[i];
        assert((isp >= 0)&&(isp < nsp));
        s[isp] = max(s[isp],mrs->UNPACK_BUF[unpack_size+i]);
      }
      
      unpack_size += spPrcommVec[ii].unpackVec.size();

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

void updateSpData(T * s, const UpdateAction action = ADD_DATA) {
  updateSpDataStart(s,action);
  updateSpDataFinish(s,action);
}

void updateStDataStart(T * s, const UpdateAction action = ADD_DATA) {

  // all stPrcomm's are non-local at this point -- no periodicity...
  // and perhaps it does not matter. We should use the buffer infrastructure
  // anyways, even with periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateStDataStart(T*): s is already mapped. updateStDataFinish required on s.");
  
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;
  
  mrs->sendRequestArray = new MPI_Request[stPrcommVec.size()];
  mrs->recvRequestArray = new MPI_Request[stPrcommVec.size()];
  //mrs->statusArray = new MPI_Status[stPrcommVec.size()];
    
  int unpack_size = 0;
  int pack_size = 0;
  for (int ii = 0; ii < stPrcommVec.size(); ++ii) {
    unpack_size += stPrcommVec[ii].unpackVec.size();
    pack_size += stPrcommVec[ii].packVec.size();
  }
  
  // and pack and send...
  
  assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
  
  unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < stPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv(mrs->UNPACK_BUF+unpack_size,stPrcommVec[ii].unpackVec.size(),
	      MPI_T,stPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
    
    unpack_size += stPrcommVec[ii].unpackVec.size();
    
    // pack and send...
    
    for (int i = 0; i < stPrcommVec[ii].packVec.size(); ++i) {
      const int ist = stPrcommVec[ii].packVec[i];
      assert((ist >= 0)&&(ist < nst));
      mrs->PACK_BUF[pack_size+i] = s[ist];
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,stPrcommVec[ii].packVec.size(),
	       MPI_T,stPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += stPrcommVec[ii].packVec.size();

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateStDataFinish(T * s, const UpdateAction action = ADD_DATA) {

  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateStDataFinish(T*): s is not mapped.");
  
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
  
  // now we wait for all messages to be sent and received...
  
  MPI_Waitall(stPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);

  // as soon as we are all sent, we can clear the buffer...
  
  delete[] mrs->PACK_BUF;         mrs->PACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;

    // and finally, wait for the recv's...
    
  MPI_Waitall(stPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);
  
  // once all recv's available, we need to unpack...
  
  int unpack_size = 0;
  switch (action) {
  case ADD_DATA:
    for (int ii = 0; ii < stPrcommVec.size(); ++ii) {
      
      for (int i = 0; i < stPrcommVec[ii].unpackVec.size(); ++i) {
        const int ist = stPrcommVec[ii].unpackVec[i];
        assert((ist >= 0)&&(ist < nst));
        s[ist] += mrs->UNPACK_BUF[unpack_size+i]; // default for nodal update is addition
      }
      
      unpack_size += stPrcommVec[ii].unpackVec.size();

    }
    break;
  case MIN_DATA:
    for (int ii = 0; ii < stPrcommVec.size(); ++ii) {
      
      for (int i = 0; i < stPrcommVec[ii].unpackVec.size(); ++i) {
        const int ist = stPrcommVec[ii].unpackVec[i];
        assert((ist >= 0)&&(ist < nst));
        s[ist] = min(s[ist],mrs->UNPACK_BUF[unpack_size+i]);
      }
      
      unpack_size += stPrcommVec[ii].unpackVec.size();

    }
    break;
  case MAX_DATA:
    for (int ii = 0; ii < stPrcommVec.size(); ++ii) {
      
      for (int i = 0; i < stPrcommVec[ii].unpackVec.size(); ++i) {
        const int ist = stPrcommVec[ii].unpackVec[i];
        assert((ist >= 0)&&(ist < nst));
        s[ist] = max(s[ist],mrs->UNPACK_BUF[unpack_size+i]);
      }
      
      unpack_size += stPrcommVec[ii].unpackVec.size();

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

void updateStData(T * s, const UpdateAction action = ADD_DATA) {
  updateStDataStart(s,action);
  updateStDataFinish(s,action);
}
