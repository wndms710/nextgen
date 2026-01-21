// =====================================================================
// these routines use the mpi_comm to exchange data... 
// =====================================================================

void updateSpDataStart(T (*s)[3], const UpdateAction action = ADD_ROTATE_DATA) {

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
    unpack_size += spPrcommVec[ii].unpackVec.size()*3;
    pack_size += spPrcommVec[ii].packVec.size()*3;
  }

  // and pack and send...
    
  assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
    
  unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < spPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv(mrs->UNPACK_BUF+unpack_size,spPrcommVec[ii].unpackVec.size()*3,
	      MPI_T,spPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));

    unpack_size += spPrcommVec[ii].unpackVec.size()*3;
    
    // pack and send...
    
    // any periodic transformations are applied on the pack side...
    switch (action) {
    case ADD_ROTATE_DATA: // default transformation
      for (int jj = 0, limit = spPrcommVec[ii].transformVec.size(); jj < limit; ++jj) {
        if (spPrcommVec[ii].transformVec[jj].has_R) {
          for (int i = spPrcommVec[ii].transformVec[jj].start; i != spPrcommVec[ii].transformVec[jj].end; ++i) {
            const int isp = spPrcommVec[ii].packVec[i];
            MiscUtils::matVecMult(mrs->PACK_BUF+pack_size+3*i,spPrcommVec[ii].transformVec[jj].R,s[isp]);
          }
        }
        else {
          for (int i = spPrcommVec[ii].transformVec[jj].start; i != spPrcommVec[ii].transformVec[jj].end; ++i) {
            const int isp = spPrcommVec[ii].packVec[i];
            FOR_J3 mrs->PACK_BUF[pack_size+3*i+j] = s[isp][j];
          }
        }
      }
      break;
    case ADD_TRANSLATE_DATA: 
      for (int jj = 0, limit = spPrcommVec[ii].transformVec.size(); jj < limit; ++jj) {
        if ((spPrcommVec[ii].transformVec[jj].has_R)&&(spPrcommVec[ii].transformVec[jj].has_t)) {
	  for (int i = spPrcommVec[ii].transformVec[jj].start; i != spPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int isp = spPrcommVec[ii].packVec[i];
            double tmp[3];
	    MiscUtils::matVecMult(tmp,spPrcommVec[ii].transformVec[jj].R,s[isp]);
	    MiscUtils::vecVecAdd(mrs->PACK_BUF+pack_size+3*i,spPrcommVec[ii].transformVec[jj].t,tmp);
	  }
        }
        else if (spPrcommVec[ii].transformVec[jj].has_R) {
	  for (int i = spPrcommVec[ii].transformVec[jj].start; i != spPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int isp = spPrcommVec[ii].packVec[i];
	    MiscUtils::matVecMult(mrs->PACK_BUF+pack_size+3*i,spPrcommVec[ii].transformVec[jj].R,s[isp]);
	  }
	}
	else if (spPrcommVec[ii].transformVec[jj].has_t) {
	  for (int i = spPrcommVec[ii].transformVec[jj].start; i != spPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int isp = spPrcommVec[ii].packVec[i];
	    MiscUtils::vecVecAdd(mrs->PACK_BUF+pack_size+3*i,spPrcommVec[ii].transformVec[jj].t,s[isp]);
	  }
	}
	else {
	  for (int i = spPrcommVec[ii].transformVec[jj].start; i != spPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int isp = spPrcommVec[ii].packVec[i];
	    FOR_J3 mrs->PACK_BUF[pack_size+3*i+j] = s[isp][j];
	  }
	}
      }
      break;
    case ADD_DATA: 
      for (int i = 0; i < spPrcommVec[ii].packVec.size(); ++i) {
        const int isp = spPrcommVec[ii].packVec[i];
        assert((isp >= 0)&&(isp < nsp));
        mrs->PACK_BUF[pack_size+3*i  ] = s[isp][0];
        mrs->PACK_BUF[pack_size+3*i+1] = s[isp][1];
        mrs->PACK_BUF[pack_size+3*i+2] = s[isp][2];
      }
      break;
    default:
      assert(0);
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,spPrcommVec[ii].packVec.size()*3,
	       MPI_T,spPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += spPrcommVec[ii].packVec.size()*3;

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateSpDataFinish(T (*s)[3], const UpdateAction action = ADD_ROTATE_DATA) {

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
  
  // unpack...

  int unpack_size = 0;
  switch (action) {
  case ADD_DATA:
  case ADD_ROTATE_DATA:
  case ADD_TRANSLATE_DATA:
    for (int ii = 0; ii < spPrcommVec.size(); ++ii) {
      
      // unpack...
      
      for (int i = 0; i < spPrcommVec[ii].unpackVec.size(); ++i) {
        const int isp = spPrcommVec[ii].unpackVec[i];
        assert((isp >= 0)&&(isp < nsp));
        s[isp][0] += mrs->UNPACK_BUF[unpack_size+3*i  ];
        s[isp][1] += mrs->UNPACK_BUF[unpack_size+3*i+1];
        s[isp][2] += mrs->UNPACK_BUF[unpack_size+3*i+2];
      }
      
      unpack_size += spPrcommVec[ii].unpackVec.size()*3;
      
    }
    break;
  default:
    assert(0);
  }
  
  // cleanup...

  delete[] mrs->UNPACK_BUF;       mrs->UNPACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
    
}

void updateSpData(T (*s)[3], const UpdateAction action = ADD_ROTATE_DATA) {
  updateSpDataStart(s,action);
  updateSpDataFinish(s,action);
}

void updateSpDataCheckFinish(T (*s)[3]) {

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
  
  // unpack...

  double my_max_dist2 = 0.0;
  for (int iter = 0; iter < 2; ++iter) {
    int unpack_size = 0;
    for (int ii = 0; ii < spPrcommVec.size(); ++ii) {
      
      // unpack...
      
      for (int i = 0; i < spPrcommVec[ii].unpackVec.size(); ++i) {
        const int isp = spPrcommVec[ii].unpackVec[i];
        assert((isp >= 0)&&(isp < nsp));
        my_max_dist2 = max(my_max_dist2,DIST2(s[isp],mrs->UNPACK_BUF+unpack_size+3*i));
        //assert(s[isp][0] == mrs->UNPACK_BUF[unpack_size+3*i  ]);
        //assert(s[isp][1] == mrs->UNPACK_BUF[unpack_size+3*i+1]);
        //assert(s[isp][2] == mrs->UNPACK_BUF[unpack_size+3*i+2]);
      }
      
      unpack_size += spPrcommVec[ii].unpackVec.size()*3;
      
    }
  }
  double max_dist2;
  MPI_Reduce(&my_max_dist2,&max_dist2,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
  if (mpi_rank == 0) 
    cout << " > sp communicator x check (should be zero): " << sqrt(max_dist2) << endl;

  // cleanup...

  delete[] mrs->UNPACK_BUF;       mrs->UNPACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
    
}

void updateSpDataCheck(T (*s)[3]) {
  updateSpDataStart(s,ADD_TRANSLATE_DATA);
  updateSpDataCheckFinish(s);
}

void updateStDataStart(T (*s)[3], const UpdateAction action = ADD_ROTATE_DATA) {

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
    unpack_size += stPrcommVec[ii].unpackVec.size()*3;
    pack_size += stPrcommVec[ii].packVec.size()*3;
  }

  // and pack and send...
    
  assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
    
  unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < stPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv(mrs->UNPACK_BUF+unpack_size,stPrcommVec[ii].unpackVec.size()*3,
	      MPI_T,stPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));

    unpack_size += stPrcommVec[ii].unpackVec.size()*3;
    
    // pack and send...
    
    // any periodic transformations are applied on the pack side...
    switch (action) {
    case ADD_ROTATE_DATA: // default transformation
      for (int jj = 0, limit = stPrcommVec[ii].transformVec.size(); jj < limit; ++jj) {
        if (stPrcommVec[ii].transformVec[jj].has_R) {
          for (int i = stPrcommVec[ii].transformVec[jj].start; i != stPrcommVec[ii].transformVec[jj].end; ++i) {
            const int ist = stPrcommVec[ii].packVec[i];
            MiscUtils::matVecMult(mrs->PACK_BUF+pack_size+3*i,stPrcommVec[ii].transformVec[jj].R,s[ist]);
          }
        }
        else {
          for (int i = stPrcommVec[ii].transformVec[jj].start; i != stPrcommVec[ii].transformVec[jj].end; ++i) {
            const int ist = stPrcommVec[ii].packVec[i];
            FOR_J3 mrs->PACK_BUF[pack_size+3*i+j] = s[ist][j];
          }
        }
      }
      break;
    case ADD_TRANSLATE_DATA: 
      for (int jj = 0, limit = stPrcommVec[ii].transformVec.size(); jj < limit; ++jj) {
        if ((stPrcommVec[ii].transformVec[jj].has_R)&&(stPrcommVec[ii].transformVec[jj].has_t)) {
	  for (int i = stPrcommVec[ii].transformVec[jj].start; i != stPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int ist = stPrcommVec[ii].packVec[i];
            double tmp[3];
	    MiscUtils::matVecMult(tmp,stPrcommVec[ii].transformVec[jj].R,s[ist]);
	    MiscUtils::vecVecAdd(mrs->PACK_BUF+pack_size+3*i,stPrcommVec[ii].transformVec[jj].t,tmp);
	  }
        }
        else if (stPrcommVec[ii].transformVec[jj].has_R) {
	  for (int i = stPrcommVec[ii].transformVec[jj].start; i != stPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int ist = stPrcommVec[ii].packVec[i];
	    MiscUtils::matVecMult(mrs->PACK_BUF+pack_size+3*i,stPrcommVec[ii].transformVec[jj].R,s[ist]);
	  }
	}
	else if (stPrcommVec[ii].transformVec[jj].has_t) {
	  for (int i = stPrcommVec[ii].transformVec[jj].start; i != stPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int ist = stPrcommVec[ii].packVec[i];
	    MiscUtils::vecVecAdd(mrs->PACK_BUF+pack_size+3*i,stPrcommVec[ii].transformVec[jj].t,s[ist]);
	  }
	}
	else {
	  for (int i = stPrcommVec[ii].transformVec[jj].start; i != stPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int ist = stPrcommVec[ii].packVec[i];
	    FOR_J3 mrs->PACK_BUF[pack_size+3*i+j] = s[ist][j];
	  }
	}
      }
      break;
    case ADD_DATA: 
      for (int i = 0; i < stPrcommVec[ii].packVec.size(); ++i) {
        const int ist = stPrcommVec[ii].packVec[i];
        assert((ist >= 0)&&(ist < nst));
        mrs->PACK_BUF[pack_size+3*i  ] = s[ist][0];
        mrs->PACK_BUF[pack_size+3*i+1] = s[ist][1];
        mrs->PACK_BUF[pack_size+3*i+2] = s[ist][2];
      }
      break;
    default:
      assert(0);
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,stPrcommVec[ii].packVec.size()*3,
	       MPI_T,stPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += stPrcommVec[ii].packVec.size()*3;

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateStDataFinish(T (*s)[3], const UpdateAction action = ADD_ROTATE_DATA) {

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
  
  // unpack...

  int unpack_size = 0;
  switch (action) {
  case ADD_DATA:
  case ADD_ROTATE_DATA:
  case ADD_TRANSLATE_DATA:
    for (int ii = 0; ii < stPrcommVec.size(); ++ii) {
      
      // unpack...
      
      for (int i = 0; i < stPrcommVec[ii].unpackVec.size(); ++i) {
        const int ist = stPrcommVec[ii].unpackVec[i];
        assert((ist >= 0)&&(ist < nst));
        s[ist][0] += mrs->UNPACK_BUF[unpack_size+3*i  ];
        s[ist][1] += mrs->UNPACK_BUF[unpack_size+3*i+1];
        s[ist][2] += mrs->UNPACK_BUF[unpack_size+3*i+2];
      }
      
      unpack_size += stPrcommVec[ii].unpackVec.size()*3;
      
    }
    break;
  default:
    assert(0);
  }
  
  // cleanup...

  delete[] mrs->UNPACK_BUF;       mrs->UNPACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
    
}

void updateStData(T (*s)[3], const UpdateAction action = ADD_ROTATE_DATA) {
  updateStDataStart(s,action);
  updateStDataFinish(s,action);
}

void updateStDataCheckFinish(T (*s)[3]) {

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
  
  // unpack...

  double my_max_dist2 = 0.0;
  for (int iter = 0; iter < 2; ++iter) {
    int unpack_size = 0;
    for (int ii = 0; ii < stPrcommVec.size(); ++ii) {
      
      // unpack...
      
      for (int i = 0; i < stPrcommVec[ii].unpackVec.size(); ++i) {
        const int ist = stPrcommVec[ii].unpackVec[i];
        assert((ist >= 0)&&(ist < nst));
        my_max_dist2 = max(my_max_dist2,DIST2(s[ist],mrs->UNPACK_BUF+unpack_size+3*i));
        //assert(s[ist][0] == mrs->UNPACK_BUF[unpack_size+3*i  ]);
        //assert(s[ist][1] == mrs->UNPACK_BUF[unpack_size+3*i+1]);
        //assert(s[ist][2] == mrs->UNPACK_BUF[unpack_size+3*i+2]);
      }
      
      unpack_size += stPrcommVec[ii].unpackVec.size()*3;
      
    }
  }
  double max_dist2;
  MPI_Reduce(&my_max_dist2,&max_dist2,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
  if (mpi_rank == 0) 
    cout << " > st communicator x check (should be zero): " << sqrt(max_dist2) << endl;

  // cleanup...

  delete[] mrs->UNPACK_BUF;       mrs->UNPACK_BUF = NULL; // nullify for desructor check
  delete[] mrs->recvRequestArray; mrs->recvRequestArray = NULL;
  
  // the destructor checks that all stuff has been nullified...
  
  delete mrs;
    
}

void updateStDataCheck(T (*s)[3]) {
  updateStDataStart(s,ADD_TRANSLATE_DATA);
  updateStDataCheckFinish(s);
}
