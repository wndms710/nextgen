// =====================================================================
// these routines use the mpi_comm to exchange data... 
// =====================================================================

void updateNoDataStart(T * s, const UpdateAction action = ADD_NO_PERIODIC_DATA) {

  // all noPrcomm's are non-local at this point -- no periodicity...
  // and perhaps it does not matter. We should use the buffer infrastructure
  // anyways, even with periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateNoDataStart(T*): s is already mapped. updateNoDataFinish required on s.");
  
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;
  
  mrs->sendRequestArray = new MPI_Request[noPrcommVec.size()];
  mrs->recvRequestArray = new MPI_Request[noPrcommVec.size()];
  //mrs->statusArray = new MPI_Status[noPrcommVec.size()];
    
  int unpack_size = 0;
  int pack_size = 0;
  for (int ii = 0; ii < noPrcommVec.size(); ++ii) {
    unpack_size += noPrcommVec[ii].unpackVec.size();
    pack_size += noPrcommVec[ii].packVec.size();
  }
  
  // and pack and send...
  
  assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
  
  unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < noPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv(mrs->UNPACK_BUF+unpack_size,noPrcommVec[ii].unpackVec.size(),
	      MPI_T,noPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
    
    unpack_size += noPrcommVec[ii].unpackVec.size();
    
    // pack and send...
    
    for (int i = 0; i < noPrcommVec[ii].packVec.size(); ++i) {
      const int ino = noPrcommVec[ii].packVec[i];
      assert((ino >= 0)&&(ino < nno));
      mrs->PACK_BUF[pack_size+i] = s[ino];
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,noPrcommVec[ii].packVec.size(),
	       MPI_T,noPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += noPrcommVec[ii].packVec.size();

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateNoDataFinish(T * s, const UpdateAction action = ADD_NO_PERIODIC_DATA) {

  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateNoDataFinish(T*): s is not mapped.");
  
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
  
  // now we wait for all messages to be sent and received...
  
  MPI_Waitall(noPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);

  // as soon as we are all sent, we can clear the buffer...
  
  delete[] mrs->PACK_BUF;         mrs->PACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;

    // and finally, wait for the recv's...
    
  MPI_Waitall(noPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);
  
  // once all recv's available, we need to unpack...
  
  int unpack_size = 0;
  if (action == ADD_NO_PERIODIC_DATA) {
    for (int ii = 0; ii < noPrcommVec.size(); ++ii) {
      
      for (int i = 0; i < noPrcommVec[ii].unpackVec.size(); ++i) {
        const int ino = noPrcommVec[ii].unpackVec[i];
        assert((ino >= 0)&&(ino < nno));
        s[ino] += mrs->UNPACK_BUF[unpack_size+i]; // default for nodal update is addition
      }
      
      unpack_size += noPrcommVec[ii].unpackVec.size();

    }
  }
  else if (action == MIN_NO_PERIODIC_DATA) { 
    for (int ii = 0; ii < noPrcommVec.size(); ++ii) {
      
      for (int i = 0; i < noPrcommVec[ii].unpackVec.size(); ++i) {
        const int ino = noPrcommVec[ii].unpackVec[i];
        assert((ino >= 0)&&(ino < nno));
        s[ino] = min(s[ino],mrs->UNPACK_BUF[unpack_size+i]); 
      }
      
      unpack_size += noPrcommVec[ii].unpackVec.size();

    }
  }
  else if (action == MAX_NO_PERIODIC_DATA) { 
    for (int ii = 0; ii < noPrcommVec.size(); ++ii) {
      
      for (int i = 0; i < noPrcommVec[ii].unpackVec.size(); ++i) {
        const int ino = noPrcommVec[ii].unpackVec[i];
        assert((ino >= 0)&&(ino < nno));
        s[ino] = max(s[ino],mrs->UNPACK_BUF[unpack_size+i]); 
      }
      
      unpack_size += noPrcommVec[ii].unpackVec.size();

    }
  }
  else {
    CERR("Unsupported UpdateAction in updateNoData");
  }

  // cleanup...
  
  delete[] mrs->UNPACK_BUF; mrs->UNPACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
  
}

void updateNoData(T * s, const UpdateAction action = ADD_NO_PERIODIC_DATA) {
  updateNoDataStart(s,action);
  updateNoDataFinish(s,action);
}

void updateNoDataStart(T (*s)[3], const UpdateAction action = ADD_NO_PERIODIC_DATA) {

  // all noPrcomm's are non-local at this point -- no periodicity...
  // and perhaps it does not matter. We should use the buffer infrastructure
  // anyways, even with periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateNoDataStart(T*): s is already mapped. updateNoDataFinish required on s.");
    
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;
  
  mrs->sendRequestArray = new MPI_Request[noPrcommVec.size()];
  mrs->recvRequestArray = new MPI_Request[noPrcommVec.size()];
  //mrs->statusArray = new MPI_Status[noPrcommVec.size()];
    
  int unpack_size = 0;
  int pack_size = 0;
  for (int ii = 0; ii < noPrcommVec.size(); ++ii) {
    unpack_size += noPrcommVec[ii].unpackVec.size()*3;
    pack_size += noPrcommVec[ii].packVec.size()*3;
  }

  // and pack and send...
    
  assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
    
  unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < noPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv(mrs->UNPACK_BUF+unpack_size,noPrcommVec[ii].unpackVec.size()*3,
	      MPI_T,noPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));

    unpack_size += noPrcommVec[ii].unpackVec.size()*3;
    
    // pack and send...
    
    for (int i = 0; i < noPrcommVec[ii].packVec.size(); ++i) {
      const int ino = noPrcommVec[ii].packVec[i];
      assert((ino >= 0)&&(ino < nno));
      mrs->PACK_BUF[pack_size+3*i  ] = s[ino][0];
      mrs->PACK_BUF[pack_size+3*i+1] = s[ino][1];
      mrs->PACK_BUF[pack_size+3*i+2] = s[ino][2];
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,noPrcommVec[ii].packVec.size()*3,
	       MPI_T,noPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += noPrcommVec[ii].packVec.size()*3;

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateNoDataFinish(T (*s)[3], const UpdateAction action = ADD_NO_PERIODIC_DATA) {

  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateNoDataFinish(T*): s is not mapped.");
  
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
  
  // now we wait for all messages to be sent and received...
  
  MPI_Waitall(noPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);
  
  // as soon as we are all sent, we can clear the buffer...
  
  delete[] mrs->PACK_BUF;         mrs->PACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;
  
  // and finally, wait for the recv's...
  
  MPI_Waitall(noPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);
  
  // unpack...

  int unpack_size = 0;
  if (action == ADD_NO_PERIODIC_DATA) {
    for (int ii = 0; ii < noPrcommVec.size(); ++ii) {
      
      // unpack...
      
      for (int i = 0; i < noPrcommVec[ii].unpackVec.size(); ++i) {
        const int ino = noPrcommVec[ii].unpackVec[i];
        assert((ino >= 0)&&(ino < nno));
        s[ino][0] += mrs->UNPACK_BUF[unpack_size+3*i  ];
        s[ino][1] += mrs->UNPACK_BUF[unpack_size+3*i+1];
        s[ino][2] += mrs->UNPACK_BUF[unpack_size+3*i+2];
      }
      
      unpack_size += noPrcommVec[ii].unpackVec.size()*3;
      
    }
  }
  else if (action == MIN_NO_PERIODIC_DATA) { 
    for (int ii = 0; ii < noPrcommVec.size(); ++ii) {
      
      // unpack...
      
      for (int i = 0; i < noPrcommVec[ii].unpackVec.size(); ++i) {
        const int ino = noPrcommVec[ii].unpackVec[i];
        assert((ino >= 0)&&(ino < nno));
        s[ino][0] = min(s[ino][0],mrs->UNPACK_BUF[unpack_size+3*i  ]);
        s[ino][1] = min(s[ino][1],mrs->UNPACK_BUF[unpack_size+3*i+1]);
        s[ino][2] = min(s[ino][2],mrs->UNPACK_BUF[unpack_size+3*i+2]);
      }
      
      unpack_size += noPrcommVec[ii].unpackVec.size()*3;
      
    }
  }
  else if (action == MAX_NO_PERIODIC_DATA) { 
    for (int ii = 0; ii < noPrcommVec.size(); ++ii) {
      
      // unpack...
      
      for (int i = 0; i < noPrcommVec[ii].unpackVec.size(); ++i) {
        const int ino = noPrcommVec[ii].unpackVec[i];
        assert((ino >= 0)&&(ino < nno));
        s[ino][0] = max(s[ino][0],mrs->UNPACK_BUF[unpack_size+3*i  ]);
        s[ino][1] = max(s[ino][1],mrs->UNPACK_BUF[unpack_size+3*i+1]);
        s[ino][2] = max(s[ino][2],mrs->UNPACK_BUF[unpack_size+3*i+2]);
      }
      
      unpack_size += noPrcommVec[ii].unpackVec.size()*3;
      
    }
  }
  else {
    CERR("Unsupported UpdateAction in updateNoData");
  }
  
  // cleanup...

  delete[] mrs->UNPACK_BUF;       mrs->UNPACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
    
}

void updateNoData(T (*s)[3], const UpdateAction action = ADD_NO_PERIODIC_DATA) {
  updateNoDataStart(s,action);
  updateNoDataFinish(s,action);
}

void updateNoDataCheckFinish(T (*s)[3]) {

  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateNoDataFinish(T*): s is not mapped.");
  
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
  
  // now we wait for all messages to be sent and received...
  
  MPI_Waitall(noPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);
  
  // as soon as we are all sent, we can clear the buffer...
  
  delete[] mrs->PACK_BUF;         mrs->PACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;
  
  // and finally, wait for the recv's...
  
  MPI_Waitall(noPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);
  
  // unpack...

  int unpack_size = 0;
  for (int ii = 0; ii < noPrcommVec.size(); ++ii) {
    
    // unpack...
    
    for (int i = 0; i < noPrcommVec[ii].unpackVec.size(); ++i) {
      const int ino = noPrcommVec[ii].unpackVec[i];
      assert((ino >= 0)&&(ino < nno));
      assert(s[ino][0] == mrs->UNPACK_BUF[unpack_size+3*i  ]);
      assert(s[ino][1] == mrs->UNPACK_BUF[unpack_size+3*i+1]);
      assert(s[ino][2] == mrs->UNPACK_BUF[unpack_size+3*i+2]);
    }
    
    unpack_size += noPrcommVec[ii].unpackVec.size()*3;
    
  }
  
  // cleanup...

  delete[] mrs->UNPACK_BUF;       mrs->UNPACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
    
}

void updateNoDataCheck(T (*s)[3]) {
  updateNoDataStart(s);
  updateNoDataCheckFinish(s);
}
