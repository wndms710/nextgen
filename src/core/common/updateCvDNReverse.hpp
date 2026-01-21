// =====================================================================
// these routines use the mpi_comm to exchange data... 
// =====================================================================

void updateCvDataReverseStart(T * s,const int action) {

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
    
  //int unpack_size = 0;
  int pack_size = 0;
  for (int ii = 0; ii < cvPrcommVec.size(); ++ii) {
    pack_size += cvPrcommVec[ii].packVec.size();
  }

  // and pack and send...
    
  //assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
  
  //unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < cvPrcommVec.size(); ++ii) {
      
    // post irecv in pack buf...
    MPI_Irecv(mrs->PACK_BUF+pack_size,cvPrcommVec[ii].packVec.size(),
	      MPI_T,cvPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
    pack_size += cvPrcommVec[ii].packVec.size();
    
    // send back ghost data from s...
    MPI_Issend(s+cvPrcommVec[ii].unpack_offset,cvPrcommVec[ii].unpack_size,
	      MPI_T,cvPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
      
  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateCvDataReverseFinish(T * s,const int action) {

  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateCvDataFinish(T*): s is not mapped.");
    
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
   
  // now we wait for all messages to be sent and received...
    
  MPI_Waitall(cvPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);
  
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;
    
  // and finally, wait for the recv's...
    
  MPI_Waitall(cvPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);
    
  // once all recv's available, we need to (un)pack...
  
  switch (action) {
  case BITWISE_OR_DATA:
    {
      int pack_size = 0;
      for (int ii = 0; ii < cvPrcommVec.size(); ++ii) {
        for (int i = 0; i < cvPrcommVec[ii].packVec.size(); ++i) {
          const int icv = cvPrcommVec[ii].packVec[i];
          assert((icv >= 0)&&(icv < ncv)); // should be internal
          s[icv] |= mrs->PACK_BUF[pack_size+i];
        }
        pack_size += cvPrcommVec[ii].packVec.size();
      }
    }
    break;
  default:
    assert(0);
  }
  
  // cleanup...
  
  delete[] mrs->PACK_BUF;         mrs->PACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
    
}

void updateCvDataReverse(T * s,const int action) {
  updateCvDataReverseStart(s,action);
  updateCvDataReverseFinish(s,action);
}
