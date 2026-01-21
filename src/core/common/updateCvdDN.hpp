// =====================================================================
// these routines use the mpi_comm to exchange data... 
// =====================================================================

void updateCvdDataSeparateGhostsStart(T * s,T *sg) {

  // all cvdPrcomm's are non-local at this point -- no periodicity...
  // and perhaps it does not matter. We should use the buffer infrastructure
  // anyways, even with periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateCvdDataSeparateGhostsStart(T*,T*): s/sg is already mapped. updateCvdDataSeparateGhostsFinish required on s/sg.");
  
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;
  
  mrs->sendRequestArray = new MPI_Request[cvdPrcommVec.size()];
  mrs->recvRequestArray = new MPI_Request[cvdPrcommVec.size()];
  //mrs->statusArray = new MPI_Status[cvdPrcommVec.size()];
    
  //int unpack_size = 0;
  int pack_size = 0;
  for (int ii = 0; ii < cvdPrcommVec.size(); ++ii) {
    //unpack_size += cvdPrcommVec[ii].unpackVec.size();
    pack_size += cvdPrcommVec[ii].packVec.size();
  }

  // and pack and send...
    
  //assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
  
  //unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < cvdPrcommVec.size(); ++ii) {
      
    // post irecv -- directly into sg...
    MPI_Irecv((sg-ncv)+cvdPrcommVec[ii].unpack_offset,cvdPrcommVec[ii].unpack_size,
	      MPI_T,cvdPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
      
    //unpack_size += cvdPrcommVec[ii].unpackVec.size();
    
    // pack and send...
    
    for (int i = 0; i < cvdPrcommVec[ii].packVec.size(); ++i) {
      const int icv = cvdPrcommVec[ii].packVec[i];
      assert((icv >= 0)&&(icv < ncv));
      mrs->PACK_BUF[pack_size+i] = s[icv];
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,cvdPrcommVec[ii].packVec.size(),
	       MPI_T,cvdPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += cvdPrcommVec[ii].packVec.size();

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateCvdDataSeparateGhostsFinish(T * s,T * sg) {

  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateCvdDataSeparateGhostsFinish(T*,T*): s/sg is not mapped.");
    
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
   
  // now we wait for all messages to be sent and received...
    
  MPI_Waitall(cvdPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);

  // as soon as we are all sent, we can clear the buffer...

  delete[] mrs->PACK_BUF;         mrs->PACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;
    
  // and finally, wait for the recv's...
    
  MPI_Waitall(cvdPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);
    
  // no unpacking into ghost region required...

  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
    
}

void updateCvdDataSeparateGhosts(T * s,T * sg) {
  updateCvdDataSeparateGhostsStart(s,sg);
  updateCvdDataSeparateGhostsFinish(s,sg);
}

void updateCvdDataStart(T *s) {
  updateCvdDataSeparateGhostsStart(s,s+ncv);
}

void updateCvdDataFinish(T *s) {
  updateCvdDataSeparateGhostsFinish(s,s+ncv);
}

void updateCvdData(T *s) {
  updateCvdDataStart(s);
  updateCvdDataFinish(s);
}
