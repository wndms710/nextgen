// =====================================================================
// these routines use the mpi_comm to exchange data... 
// =====================================================================

void updateCcIDataReverseStart(T * s) {

  // all ccIPrcomm's are non-local at this point -- no periodicity...
  // and perhaps it does not matter. We should use the buffer infrastructure
  // anyways, even with periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateCcIDataReverseStart(T*): s is already mapped. updateCvDataFinish required on s.");
  
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;
  
  mrs->sendRequestArray = new MPI_Request[ccIPrcommVec.size()];
  mrs->recvRequestArray = new MPI_Request[ccIPrcommVec.size()];
    
  int pack_size = 0;
  for (int ii = 0; ii < ccIPrcommVec.size(); ++ii) {
    pack_size += ccIPrcommVec[ii].unpackVec.size(); // note its bawckwards!
  }

  // and pack and send...
    
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
  
  pack_size = 0;
  for (int ii = 0; ii < ccIPrcommVec.size(); ++ii) {
      
    // post irecv -- directly into s...
    MPI_Irecv(s+ccIPrcommVec[ii].pack_offset,ccIPrcommVec[ii].pack_size,
	      MPI_T,ccIPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
      
    
    // pack and send...
    
    for (int i = 0; i < ccIPrcommVec[ii].unpackVec.size(); ++i) {
      const int icc = ccIPrcommVec[ii].unpackVec[i];
      assert((icc >= 0)&&(icc < ncc_a));
      mrs->PACK_BUF[pack_size+i] = s[icc];
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,ccIPrcommVec[ii].unpackVec.size(),
	       MPI_T,ccIPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += ccIPrcommVec[ii].unpackVec.size();

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateCcIDataReverseFinish(T * s) {

  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateCcIDataReverseFinish(T*): s is not mapped.");
    
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
   
  // now we wait for all messages to be sent and received...
    
  MPI_Waitall(ccIPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);

  // as soon as we are all sent, we can clear the buffer...

  delete[] mrs->PACK_BUF;         mrs->PACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;
    
  // and finally, wait for the recv's...
    
  MPI_Waitall(ccIPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);
    
  // no unpacking into ghost region required...

  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
    
}

void updateCcIDataReverse(T * s) {
  updateCcIDataReverseStart(s);
  updateCcIDataReverseFinish(s);
}

void updateCcIDataReverseStart(T (*s)[3]) {

  // all ccIPrcomm's are non-local at this point -- no periodicity...
  // and perhaps it does not matter. We should use the buffer infrastructure
  // anyways, even with periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateCcIDataReverseStart(T*): s is already mapped. updateCvDataFinish required on s.");
  
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;
  
  mrs->sendRequestArray = new MPI_Request[ccIPrcommVec.size()];
  mrs->recvRequestArray = new MPI_Request[ccIPrcommVec.size()];
    
  int pack_size = 0;
  for (int ii = 0; ii < ccIPrcommVec.size(); ++ii) {
    pack_size += 3*ccIPrcommVec[ii].unpackVec.size(); // note its bawckwards!
  }

  // and pack and send...
    
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
  
  pack_size = 0;
  for (int ii = 0; ii < ccIPrcommVec.size(); ++ii) {
      
    // post irecv -- directly into s...
    MPI_Irecv(s+ccIPrcommVec[ii].pack_offset,3*ccIPrcommVec[ii].pack_size,
	      MPI_T,ccIPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
      
    
    // pack and send...
    
    for (int i = 0; i < ccIPrcommVec[ii].unpackVec.size(); ++i) {
      const int icc = ccIPrcommVec[ii].unpackVec[i];
      assert((icc >= 0)&&(icc < ncc_a));
      FOR_J3 mrs->PACK_BUF[pack_size+3*i+j] = s[icc][j];
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,3*ccIPrcommVec[ii].unpackVec.size(),
	       MPI_T,ccIPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += 3*ccIPrcommVec[ii].unpackVec.size();

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateCcIDataReverseFinish(T (*s)[3]) {

  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateCcIDataReverseFinish(T*): s is not mapped.");
    
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
   
  // now we wait for all messages to be sent and received...
    
  MPI_Waitall(ccIPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);

  // as soon as we are all sent, we can clear the buffer...

  delete[] mrs->PACK_BUF;         mrs->PACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;
    
  // and finally, wait for the recv's...
    
  MPI_Waitall(ccIPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);
    
  // no unpacking into ghost region required...

  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
    
}

void updateCcIDataReverse(T (*s)[3]) {
  updateCcIDataReverseStart(s);
  updateCcIDataReverseFinish(s);
}
