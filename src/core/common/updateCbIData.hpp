// =====================================================================
// these routines use the mpi_comm to exchange data... 
// =====================================================================

void updateCbIDataStart(T * s, const UpdateAction action = ADD_NO_PERIODIC_DATA) {

  // all cbIPrcomm's are non-local at this point -- no periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateCbIDataStart(T*): s is already mapped. updateCbIDataFinish required on s.");
  
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;
  
  mrs->sendRequestArray = new MPI_Request[cbIPrcommVec.size()];
  mrs->recvRequestArray = new MPI_Request[cbIPrcommVec.size()];
    
  int unpack_size = 0;
  for (int ii = 0; ii < cbIPrcommVec.size(); ++ii) 
    unpack_size += cbIPrcommVec[ii].unpackVec.size();
  
  // and pack and send...
  
  assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  
  unpack_size = 0;
  for (int ii = 0; ii < cbIPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv(mrs->UNPACK_BUF+unpack_size,cbIPrcommVec[ii].unpackVec.size(),
	      MPI_T,cbIPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
    
    unpack_size += cbIPrcommVec[ii].unpackVec.size();
    
    // pack and send...
    
    MPI_Issend(s+cbIPrcommVec[ii].pack_offset,cbIPrcommVec[ii].pack_size,
	       MPI_T,cbIPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateCbIDataFinish(T * s, const UpdateAction action = ADD_NO_PERIODIC_DATA) {

  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateCbIDataFinish(T*): s is not mapped.");
  
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
  
  // now we wait for all messages to be sent and received...
  
  MPI_Waitall(cbIPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);

  // as soon as we are all sent, we can clear the buffer...
  
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;

    // and finally, wait for the recv's...
    
  MPI_Waitall(cbIPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);
  
  // once all recv's available, we need to unpack...
  
  int unpack_size = 0;
  if (action == ADD_NO_PERIODIC_DATA) {
    for (int ii = 0; ii < cbIPrcommVec.size(); ++ii) {
      
      for (int i = 0; i < cbIPrcommVec[ii].unpackVec.size(); ++i) {
        const int icb = cbIPrcommVec[ii].unpackVec[i];
        assert((icb >= 0)&&(icb < ncb_a));
        s[icb] += mrs->UNPACK_BUF[unpack_size+i]; // default for active update is addition
      }
      
      unpack_size += cbIPrcommVec[ii].unpackVec.size();

    }
  }
  else {
    CERR("Unsupported UpdateAction in updateCbIData");
  }

  // cleanup...
  
  delete[] mrs->UNPACK_BUF; mrs->UNPACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
  
}

void updateCbIData(T * s, const UpdateAction action = ADD_NO_PERIODIC_DATA) {
  updateCbIDataStart(s,action);
  updateCbIDataFinish(s,action);
}

void updateCbIDataStart(T (*s)[3], const UpdateAction action = ADD_NO_PERIODIC_DATA) {

  // all cbIPrcomm's are non-local at this point -- no periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateCbIDataStart(T*): s is already mapped. updateCbIDataFinish required on s.");
  
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;
  
  mrs->sendRequestArray = new MPI_Request[cbIPrcommVec.size()];
  mrs->recvRequestArray = new MPI_Request[cbIPrcommVec.size()];
    
  int unpack_size = 0;
  for (int ii = 0; ii < cbIPrcommVec.size(); ++ii) 
    unpack_size += cbIPrcommVec[ii].unpackVec.size()*3;
  
  // and pack and send...
  
  assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  
  unpack_size = 0;
  for (int ii = 0; ii < cbIPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv(mrs->UNPACK_BUF+unpack_size,cbIPrcommVec[ii].unpackVec.size()*3,
	      MPI_T,cbIPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
    
    unpack_size += cbIPrcommVec[ii].unpackVec.size()*3;
    
    // pack and send...
    
    MPI_Issend(((double*)s)+cbIPrcommVec[ii].pack_offset*3,cbIPrcommVec[ii].pack_size*3,
	       MPI_T,cbIPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateCbIDataFinish(T (*s)[3], const UpdateAction action = ADD_NO_PERIODIC_DATA) {

  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateCbIDataFinish(T*): s is not mapped.");
  
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
  
  // now we wait for all messages to be sent and received...
  
  MPI_Waitall(cbIPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);

  // as soon as we are all sent, we can clear the buffer...
  
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;

    // and finally, wait for the recv's...
    
  MPI_Waitall(cbIPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);
  
  // once all recv's available, we need to unpack...
  
  int unpack_size = 0;
  if (action == ADD_NO_PERIODIC_DATA) {
    for (int ii = 0; ii < cbIPrcommVec.size(); ++ii) {
      
      for (int i = 0; i < cbIPrcommVec[ii].unpackVec.size(); ++i) {
        const int icb = cbIPrcommVec[ii].unpackVec[i];
        assert((icb >= 0)&&(icb < ncb_a));
        s[icb][0] += mrs->UNPACK_BUF[unpack_size+3*i  ]; // default for active update is addition
        s[icb][1] += mrs->UNPACK_BUF[unpack_size+3*i+1]; // default for active update is addition
        s[icb][2] += mrs->UNPACK_BUF[unpack_size+3*i+2]; // default for active update is addition
      }
      
      unpack_size += cbIPrcommVec[ii].unpackVec.size()*3;

    }
  }
  else {
    CERR("Unsupported UpdateAction in updateCbIData");
  }

  // cleanup...
  
  delete[] mrs->UNPACK_BUF; mrs->UNPACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
  
}

void updateCbIData(T (*s)[3], const UpdateAction action = ADD_NO_PERIODIC_DATA) {
  updateCbIDataStart(s,action);
  updateCbIDataFinish(s,action);
}
