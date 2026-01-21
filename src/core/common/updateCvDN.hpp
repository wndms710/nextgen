// =====================================================================
// these routines use the mpi_comm to exchange data... 
// =====================================================================

void updateCvDataStart(T * s) {

  // all cvPrcomm's are non-local at this point -- no periodicity...
  // and perhaps it does not matter. We should use the buffer infrastructure
  // anyways, even with periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateCvDataStart(T*): s is already mapped. updateCvDataFinish required on s.");
  
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;
  
  mrs->sendRequestArray = new MPI_Request[cvPrcommVec.size()];
  mrs->recvRequestArray = new MPI_Request[cvPrcommVec.size()];
  //mrs->statusArray = new MPI_Status[cvPrcommVec.size()];
    
  //int unpack_size = 0;
  int pack_size = 0;
  for (int ii = 0; ii < cvPrcommVec.size(); ++ii) {
    //unpack_size += cvPrcommVec[ii].unpackVec.size();
    pack_size += cvPrcommVec[ii].packVec.size();
  }

  // and pack and send...
    
  //assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
  
  //unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < cvPrcommVec.size(); ++ii) {
      
    // post irecv -- directly into s...
    MPI_Irecv(s+cvPrcommVec[ii].unpack_offset,cvPrcommVec[ii].unpack_size,
	      MPI_T,cvPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
      
    //unpack_size += cvPrcommVec[ii].unpackVec.size();
    
    // pack and send...
    
    for (int i = 0; i < cvPrcommVec[ii].packVec.size(); ++i) {
      const int icv = cvPrcommVec[ii].packVec[i];
      //assert((icv >= 0)&&(icv < ncv_d));
      mrs->PACK_BUF[pack_size+i] = s[icv];
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,cvPrcommVec[ii].packVec.size(),
	       MPI_T,cvPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += cvPrcommVec[ii].packVec.size();

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateCvDataFinish(T * s) {

  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateCvDataFinish(T*): s is not mapped.");
    
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
   
  // now we wait for all messages to be sent and received...
    
  MPI_Waitall(cvPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);

  // as soon as we are all sent, we can clear the buffer...

  delete[] mrs->PACK_BUF;         mrs->PACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;
    
  // and finally, wait for the recv's...
    
  MPI_Waitall(cvPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);
    
  // no unpacking into ghost region required...

  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
    
}

void updateCvData(T * s) {
  updateCvDataStart(s);
  updateCvDataFinish(s);
}
