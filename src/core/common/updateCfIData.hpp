// =====================================================================
// these routines use the mpi_comm to exchange data... 
// =====================================================================

void updateCfIDataStart(T * s, const UpdateAction action = ADD_NO_PERIODIC_DATA) {

  // all cfIPrcomm's are non-local at this point -- no periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateCfIDataStart(T*): s is already mapped. updateCfIDataFinish required on s.");
  
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;
  
  mrs->sendRequestArray = new MPI_Request[cfIPrcommVec.size()];
  mrs->recvRequestArray = new MPI_Request[cfIPrcommVec.size()];
  //mrs->statusArray = new MPI_Status[cfIPrcommVec.size()];
    
  int unpack_size = 0;
  int pack_size = 0;
  for (int ii = 0; ii < cfIPrcommVec.size(); ++ii) {
    unpack_size += cfIPrcommVec[ii].unpackVec.size();
    pack_size += cfIPrcommVec[ii].packVec.size();
  }
  
  // and pack and send...
  
  assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
  
  unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < cfIPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv(mrs->UNPACK_BUF+unpack_size,cfIPrcommVec[ii].unpackVec.size(),
	      MPI_T,cfIPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
    
    unpack_size += cfIPrcommVec[ii].unpackVec.size();
    
    // pack and send...
    
    for (int i = 0; i < cfIPrcommVec[ii].packVec.size(); ++i) {
      const int icf = max(cfIPrcommVec[ii].packVec[i],-cfIPrcommVec[ii].packVec[i]-1); // remember signed to flip 
      assert((icf >= ncf_aa)&&(icf < ncf));
      mrs->PACK_BUF[pack_size+i] = s[icf];
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,cfIPrcommVec[ii].packVec.size(),
	       MPI_T,cfIPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += cfIPrcommVec[ii].packVec.size();

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateSignedCfIDataStart(T * s, const UpdateAction action = ADD_NO_PERIODIC_DATA) {

  // all cfIPrcomm's are non-local at this point -- no periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateCfIDataStart(T*): s is already mapped. updateCfIDataFinish required on s.");
  
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;
  
  mrs->sendRequestArray = new MPI_Request[cfIPrcommVec.size()];
  mrs->recvRequestArray = new MPI_Request[cfIPrcommVec.size()];
  //mrs->statusArray = new MPI_Status[cfIPrcommVec.size()];
    
  int unpack_size = 0;
  int pack_size = 0;
  for (int ii = 0; ii < cfIPrcommVec.size(); ++ii) {
    unpack_size += cfIPrcommVec[ii].unpackVec.size();
    pack_size += cfIPrcommVec[ii].packVec.size();
  }
  
  // and pack and send...
  
  assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
  
  unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < cfIPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv(mrs->UNPACK_BUF+unpack_size,cfIPrcommVec[ii].unpackVec.size(),
	      MPI_T,cfIPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
    
    unpack_size += cfIPrcommVec[ii].unpackVec.size();
    
    // pack and send...
    
    for (int i = 0; i < cfIPrcommVec[ii].packVec.size(); ++i) {
      if (cfIPrcommVec[ii].packVec[i] >= 0) {
        const int icf = cfIPrcommVec[ii].packVec[i];
        assert((icf >= ncf_aa)&&(icf < ncf));
        mrs->PACK_BUF[pack_size+i] = s[icf];
      }
      else {
        const int icf = -cfIPrcommVec[ii].packVec[i]-1;
        assert((icf >= ncf_aa)&&(icf < ncf));
        mrs->PACK_BUF[pack_size+i] = -s[icf]; // flip data
      }
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,cfIPrcommVec[ii].packVec.size(),
	       MPI_T,cfIPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += cfIPrcommVec[ii].packVec.size();

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateCfIDataFinish(T * s, const UpdateAction action = ADD_NO_PERIODIC_DATA) {

  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateCfIDataFinish(T*): s is not mapped.");
  
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
  
  // now we wait for all messages to be sent and received...
  
  MPI_Waitall(cfIPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);

  // as soon as we are all sent, we can clear the buffer...
  
  delete[] mrs->PACK_BUF;         mrs->PACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;

    // and finally, wait for the recv's...
    
  MPI_Waitall(cfIPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);
  
  // once all recv's available, we need to unpack...
  
  int unpack_size = 0;
  if (action == ADD_NO_PERIODIC_DATA) {
    for (int ii = 0; ii < cfIPrcommVec.size(); ++ii) {
      
      for (int i = 0; i < cfIPrcommVec[ii].unpackVec.size(); ++i) {
        const int icf = cfIPrcommVec[ii].unpackVec[i];
        assert((icf >= 0)&&(icf < ncf_ai));
        s[icf] += mrs->UNPACK_BUF[unpack_size+i]; // default for update is addition
      }
      
      unpack_size += cfIPrcommVec[ii].unpackVec.size();

    }
  }
  else {
    CERR("Unsupported UpdateAction in updateCfIData");
  }

  // cleanup...
  
  delete[] mrs->UNPACK_BUF; mrs->UNPACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
  
}

void updateCfIData(T * s, const UpdateAction action = ADD_NO_PERIODIC_DATA) {
  updateCfIDataStart(s,action);
  updateCfIDataFinish(s,action);
}

void updateSignedCfIData(T * s, const UpdateAction action = ADD_NO_PERIODIC_DATA) {
  updateSignedCfIDataStart(s,action);
  updateCfIDataFinish(s,action);
}

void updateCfIDataStart(T (*s)[3], const UpdateAction action = ADD_NO_PERIODIC_DATA) {

  // all cfIPrcomm's are non-local at this point -- no periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateCfIDataStart(T*): s is already mapped. updateCfIDataFinish required on s.");
    
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;
  
  mrs->sendRequestArray = new MPI_Request[cfIPrcommVec.size()];
  mrs->recvRequestArray = new MPI_Request[cfIPrcommVec.size()];
  //mrs->statusArray = new MPI_Status[cfIPrcommVec.size()];
    
  int unpack_size = 0;
  int pack_size = 0;
  for (int ii = 0; ii < cfIPrcommVec.size(); ++ii) {
    unpack_size += cfIPrcommVec[ii].unpackVec.size()*3;
    pack_size += cfIPrcommVec[ii].packVec.size()*3;
  }

  // and pack and send...
    
  assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
    
  unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < cfIPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv(mrs->UNPACK_BUF+unpack_size,cfIPrcommVec[ii].unpackVec.size()*3,
	      MPI_T,cfIPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));

    unpack_size += cfIPrcommVec[ii].unpackVec.size()*3;
    
    // pack and send...
    
    for (int i = 0; i < cfIPrcommVec[ii].packVec.size(); ++i) {
      const int icf = max(cfIPrcommVec[ii].packVec[i],-cfIPrcommVec[ii].packVec[i]-1); // remember signed to flip 
      assert((icf >= ncf_aa)&&(icf < ncf));
      mrs->PACK_BUF[pack_size+3*i  ] = s[icf][0];
      mrs->PACK_BUF[pack_size+3*i+1] = s[icf][1];
      mrs->PACK_BUF[pack_size+3*i+2] = s[icf][2];
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,cfIPrcommVec[ii].packVec.size()*3,
	       MPI_T,cfIPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += cfIPrcommVec[ii].packVec.size()*3;

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateSignedCfIDataStart(T (*s)[3], const UpdateAction action = ADD_NO_PERIODIC_DATA) {

  // all cfIPrcomm's are non-local at this point -- no periodicity...
    
  // make sure s is not already in the map...
  
  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter != mpiRequestMap.end())
    CERR("updateCfIDataStart(T*): s is already mapped. updateCfIDataFinish required on s.");
    
  MpiRequestStuff * mrs = new MpiRequestStuff; // mrs == MpiRequestStuff
  mpiRequestMap[(void*)s] = mrs;
  
  mrs->sendRequestArray = new MPI_Request[cfIPrcommVec.size()];
  mrs->recvRequestArray = new MPI_Request[cfIPrcommVec.size()];
  //mrs->statusArray = new MPI_Status[cfIPrcommVec.size()];
    
  int unpack_size = 0;
  int pack_size = 0;
  for (int ii = 0; ii < cfIPrcommVec.size(); ++ii) {
    unpack_size += cfIPrcommVec[ii].unpackVec.size()*3;
    pack_size += cfIPrcommVec[ii].packVec.size()*3;
  }

  // and pack and send...
    
  assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
    
  unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < cfIPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv(mrs->UNPACK_BUF+unpack_size,cfIPrcommVec[ii].unpackVec.size()*3,
	      MPI_T,cfIPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));

    unpack_size += cfIPrcommVec[ii].unpackVec.size()*3;
    
    // pack and send...
    
    for (int i = 0; i < cfIPrcommVec[ii].packVec.size(); ++i) {
      if (cfIPrcommVec[ii].packVec[i] >= 0) {
        const int icf = cfIPrcommVec[ii].packVec[i];
        assert((icf >= ncf_aa)&&(icf < ncf));
        mrs->PACK_BUF[pack_size+3*i  ] = s[icf][0];
        mrs->PACK_BUF[pack_size+3*i+1] = s[icf][1];
        mrs->PACK_BUF[pack_size+3*i+2] = s[icf][2];
      }
      else {
        const int icf = -cfIPrcommVec[ii].packVec[i]-1;
        assert((icf >= ncf_aa)&&(icf < ncf));
        mrs->PACK_BUF[pack_size+3*i  ] = -s[icf][0]; // flip data
        mrs->PACK_BUF[pack_size+3*i+1] = -s[icf][1];
        mrs->PACK_BUF[pack_size+3*i+2] = -s[icf][2];
      }
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,cfIPrcommVec[ii].packVec.size()*3,
	       MPI_T,cfIPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += cfIPrcommVec[ii].packVec.size()*3;

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateCfIDataFinish(T (*s)[3], const UpdateAction action = ADD_NO_PERIODIC_DATA) {

  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateCfIDataFinish(T*): s is not mapped.");
  
  MpiRequestStuff * mrs = iter->second;
  mpiRequestMap.erase(iter);
  
  // now we wait for all messages to be sent and received...
  
  MPI_Waitall(cfIPrcommVec.size(),mrs->sendRequestArray,MPI_STATUSES_IGNORE);
  
  // as soon as we are all sent, we can clear the buffer...
  
  delete[] mrs->PACK_BUF;         mrs->PACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->sendRequestArray; mrs->sendRequestArray = NULL;
  
  // and finally, wait for the recv's...
  
  MPI_Waitall(cfIPrcommVec.size(),mrs->recvRequestArray,MPI_STATUSES_IGNORE);
  
  // unpack...

  int unpack_size = 0;
  if (action == ADD_NO_PERIODIC_DATA) {
    for (int ii = 0; ii < cfIPrcommVec.size(); ++ii) {
      
      // unpack...
      
      for (int i = 0; i < cfIPrcommVec[ii].unpackVec.size(); ++i) {
        const int icf = cfIPrcommVec[ii].unpackVec[i];
        assert((icf >= 0)&&(icf < ncf_ai));
        s[icf][0] += mrs->UNPACK_BUF[unpack_size+3*i  ];
        s[icf][1] += mrs->UNPACK_BUF[unpack_size+3*i+1];
        s[icf][2] += mrs->UNPACK_BUF[unpack_size+3*i+2];
      }
      
      unpack_size += cfIPrcommVec[ii].unpackVec.size()*3;
      
    }
  }
  else {
    CERR("Unsupported UpdateAction in updateCfIData");
  }
  
  // cleanup...

  delete[] mrs->UNPACK_BUF;       mrs->UNPACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
    
}

void updateCfIData(T (*s)[3], const UpdateAction action = ADD_NO_PERIODIC_DATA) {
  updateCfIDataStart(s,action);
  updateCfIDataFinish(s,action);
}

void updateSignedCfIData(T (*s)[3], const UpdateAction action = ADD_NO_PERIODIC_DATA) {
  updateSignedCfIDataStart(s,action);
  updateCfIDataFinish(s,action);
}
