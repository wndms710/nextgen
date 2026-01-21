// =====================================================================
// these routines use the mpi_comm to exchange data... 
// =====================================================================

void updateFaDataStart(T (*s)[3],const int action) {

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
    unpack_size += faPrcommVec[ii].unpackVec.size()*3;
    pack_size += faPrcommVec[ii].packVec.size()*3;
  }

  // and pack and send...
    
  assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
    
  unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < faPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv(mrs->UNPACK_BUF+unpack_size,faPrcommVec[ii].unpackVec.size()*3,
	      MPI_T,faPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));

    unpack_size += faPrcommVec[ii].unpackVec.size()*3;
    
    // pack and send...

    // any periodic transformations are applied on the pack side...
    switch (action) {
    case ADD_DATA:
    case SUBTRACT_DATA:
    case REPLACE_DATA:
      for (int i = 0; i < faPrcommVec[ii].packVec.size(); ++i) {
        const int ifa = faPrcommVec[ii].packVec[i];
        assert((ifa >= nfa_i)&&(ifa < nfa));
        mrs->PACK_BUF[pack_size+3*i  ] = s[ifa][0];
        mrs->PACK_BUF[pack_size+3*i+1] = s[ifa][1];
        mrs->PACK_BUF[pack_size+3*i+2] = s[ifa][2];
      }
      break;
    case ADD_ROTATE_DATA:
    case SUBTRACT_ROTATE_DATA:
    case REPLACE_ROTATE_DATA: // default transformation
      for (int jj = 0, limit = faPrcommVec[ii].transformVec.size(); jj < limit; ++jj) {
	if (faPrcommVec[ii].transformVec[jj].has_R) {
	  for (int i = faPrcommVec[ii].transformVec[jj].start; i != faPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int ifa = faPrcommVec[ii].packVec[i];
	    MiscUtils::matVecMult(mrs->PACK_BUF+pack_size+3*i,faPrcommVec[ii].transformVec[jj].R,s[ifa]);
	  }
	}
	else {
	  for (int i = faPrcommVec[ii].transformVec[jj].start; i != faPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int ifa = faPrcommVec[ii].packVec[i];
	    FOR_J3 mrs->PACK_BUF[pack_size+3*i+j] = s[ifa][j];
	  }
	}
      }
      break;
    case ADD_TRANSLATE_DATA:
    case SUBTRACT_TRANSLATE_DATA:
    case REPLACE_TRANSLATE_DATA:
      for (int jj = 0, limit = faPrcommVec[ii].transformVec.size(); jj < limit; ++jj) {
        if ((faPrcommVec[ii].transformVec[jj].has_R)&&(faPrcommVec[ii].transformVec[jj].has_t)) {
	  for (int i = faPrcommVec[ii].transformVec[jj].start; i != faPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int ifa = faPrcommVec[ii].packVec[i];
            double tmp[3];
	    MiscUtils::matVecMult(tmp,faPrcommVec[ii].transformVec[jj].R,s[ifa]);
	    MiscUtils::vecVecAdd(mrs->PACK_BUF+pack_size+3*i,faPrcommVec[ii].transformVec[jj].t,tmp);
	  }
        }
        else if (faPrcommVec[ii].transformVec[jj].has_R) {
	  for (int i = faPrcommVec[ii].transformVec[jj].start; i != faPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int ifa = faPrcommVec[ii].packVec[i];
	    MiscUtils::matVecMult(mrs->PACK_BUF+pack_size+3*i,faPrcommVec[ii].transformVec[jj].R,s[ifa]);
	  }
	}
	else if (faPrcommVec[ii].transformVec[jj].has_t) {
	  for (int i = faPrcommVec[ii].transformVec[jj].start; i != faPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int ifa = faPrcommVec[ii].packVec[i];
	    MiscUtils::vecVecAdd(mrs->PACK_BUF+pack_size+3*i,faPrcommVec[ii].transformVec[jj].t,s[ifa]);
	  }
	}
	else {
	  for (int i = faPrcommVec[ii].transformVec[jj].start; i != faPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int ifa = faPrcommVec[ii].packVec[i];
	    FOR_J3 mrs->PACK_BUF[pack_size+3*i+j] = s[ifa][j];
	  }
	}
      }
      break;
    default:
      assert(0); // if you get hwre, add cases above: e.g. ADD_ROTATE...
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,faPrcommVec[ii].packVec.size()*3,
	       MPI_T,faPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += faPrcommVec[ii].packVec.size()*3;

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateFaDataFinish(T (*s)[3],const int action) {

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
  
  // unpack...

  switch (action) {
  case MIN_NO_PERIODIC_DATA:
    {
      int unpack_size = 0;
      for (int ii = 0; ii < faPrcommVec.size(); ++ii) {
	// unpack...
	for (int i = 0; i < faPrcommVec[ii].unpackVec.size(); ++i) {
	  const int ifa = faPrcommVec[ii].unpackVec[i];
	  assert((ifa >= nfa_i)&&(ifa < nfa));
	  s[ifa][0] = min(s[ifa][0],mrs->UNPACK_BUF[unpack_size+3*i  ]);
	  s[ifa][1] = min(s[ifa][1],mrs->UNPACK_BUF[unpack_size+3*i+1]);
	  s[ifa][2] = min(s[ifa][2],mrs->UNPACK_BUF[unpack_size+3*i+2]);
	}
	unpack_size += faPrcommVec[ii].unpackVec.size()*3;
      }
    }
    break;
  case MAX_NO_PERIODIC_DATA:
    {
      int unpack_size = 0;
      for (int ii = 0; ii < faPrcommVec.size(); ++ii) {
	// unpack...
	for (int i = 0; i < faPrcommVec[ii].unpackVec.size(); ++i) {
	  const int ifa = faPrcommVec[ii].unpackVec[i];
	  assert((ifa >= nfa_i)&&(ifa < nfa));
	  s[ifa][0] = max(s[ifa][0],mrs->UNPACK_BUF[unpack_size+3*i  ]);
	  s[ifa][1] = max(s[ifa][1],mrs->UNPACK_BUF[unpack_size+3*i+1]);
	  s[ifa][2] = max(s[ifa][2],mrs->UNPACK_BUF[unpack_size+3*i+2]);
	}
	unpack_size += faPrcommVec[ii].unpackVec.size()*3;
      }
    }
    break;
  case ADD_ROTATE_DATA:
  case ADD_TRANSLATE_DATA:
  case ADD_DATA:
    {
      int unpack_size = 0;
      for (int ii = 0; ii < faPrcommVec.size(); ++ii) {
	// unpack...
	for (int i = 0; i < faPrcommVec[ii].unpackVec.size(); ++i) {
	  const int ifa = faPrcommVec[ii].unpackVec[i];
	  assert((ifa >= nfa_i)&&(ifa < nfa));
	  s[ifa][0] += mrs->UNPACK_BUF[unpack_size+3*i  ];
	  s[ifa][1] += mrs->UNPACK_BUF[unpack_size+3*i+1];
	  s[ifa][2] += mrs->UNPACK_BUF[unpack_size+3*i+2];
	}
	unpack_size += faPrcommVec[ii].unpackVec.size()*3;
      }
    }
    break;
  case SUBTRACT_ROTATE_DATA:
  case SUBTRACT_TRANSLATE_DATA:
  case SUBTRACT_DATA:
    {
      int unpack_size = 0;
      for (int ii = 0; ii < faPrcommVec.size(); ++ii) {
	// unpack...
	for (int i = 0; i < faPrcommVec[ii].unpackVec.size(); ++i) {
	  const int ifa = faPrcommVec[ii].unpackVec[i];
	  assert((ifa >= nfa_i)&&(ifa < nfa));
	  s[ifa][0] -= mrs->UNPACK_BUF[unpack_size+3*i  ];
	  s[ifa][1] -= mrs->UNPACK_BUF[unpack_size+3*i+1];
	  s[ifa][2] -= mrs->UNPACK_BUF[unpack_size+3*i+2];
	}
	unpack_size += faPrcommVec[ii].unpackVec.size()*3;
      }
    }
    break;
  case REPLACE_ROTATE_DATA:
  case REPLACE_TRANSLATE_DATA:
  case REPLACE_DATA:
    {
      int unpack_size = 0;
      for (int ii = 0; ii < faPrcommVec.size(); ++ii) {
	// unpack...
	for (int i = 0; i < faPrcommVec[ii].unpackVec.size(); ++i) {
	  const int ifa = faPrcommVec[ii].unpackVec[i];
	  assert((ifa >= nfa_i)&&(ifa < nfa));
	  s[ifa][0] = mrs->UNPACK_BUF[unpack_size+3*i  ];
	  s[ifa][1] = mrs->UNPACK_BUF[unpack_size+3*i+1];
	  s[ifa][2] = mrs->UNPACK_BUF[unpack_size+3*i+2];
	}
	unpack_size += faPrcommVec[ii].unpackVec.size()*3;
      }
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

void updateFaData(T (*s)[3],const int action) {
  updateFaDataStart(s,action);
  updateFaDataFinish(s,action);
}

