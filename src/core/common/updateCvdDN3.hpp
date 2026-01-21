// =====================================================================
// these routines use the mpi_comm to exchange data... 
// =====================================================================

void updateCvdDataSeparateGhostsStart(T (*s)[3],T (*sg)[3],const int action = REPLACE_ROTATE_DATA) {

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
    //unpack_size += cvdPrcommVec[ii].unpackVec.size()*3;
    pack_size += cvdPrcommVec[ii].packVec.size()*3;
  }

  // and pack and send...
    
  //assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
    
  //unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < cvdPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv((sg-ncv)+cvdPrcommVec[ii].unpack_offset,cvdPrcommVec[ii].unpack_size*3,
	      MPI_T,cvdPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
    
    // pack and send...
    
    // any periodic transformations are applied on the pack side...
    switch (action) {
    case REPLACE_ROTATE_DATA: // default transformation
      for (int jj = 0, limit = cvdPrcommVec[ii].transformVec.size(); jj < limit; ++jj) {
	// never translate...
	if (cvdPrcommVec[ii].transformVec[jj].has_R) {
	  for (int i = cvdPrcommVec[ii].transformVec[jj].start; i != cvdPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int icv = cvdPrcommVec[ii].packVec[i];
	    MiscUtils::matVecMult(mrs->PACK_BUF+pack_size+3*i,cvdPrcommVec[ii].transformVec[jj].R,s[icv]);
	  }
	}
	else {
	  for (int i = cvdPrcommVec[ii].transformVec[jj].start; i != cvdPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int icv = cvdPrcommVec[ii].packVec[i];
	    FOR_J3 mrs->PACK_BUF[pack_size+3*i+j] = s[icv][j];
	  }
	}
      }
      break;
    case REPLACE_TRANSLATE_DATA:
      for (int jj = 0, limit = cvdPrcommVec[ii].transformVec.size(); jj < limit; ++jj) {
        if ((cvdPrcommVec[ii].transformVec[jj].has_R)&&(cvdPrcommVec[ii].transformVec[jj].has_t)) {
	  for (int i = cvdPrcommVec[ii].transformVec[jj].start; i != cvdPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int icv = cvdPrcommVec[ii].packVec[i];
            double tmp[3];
	    MiscUtils::matVecMult(tmp,cvdPrcommVec[ii].transformVec[jj].R,s[icv]);
	    MiscUtils::vecVecAdd(mrs->PACK_BUF+pack_size+3*i,cvdPrcommVec[ii].transformVec[jj].t,tmp);
	  }
        }
        else if (cvdPrcommVec[ii].transformVec[jj].has_R) {
	  for (int i = cvdPrcommVec[ii].transformVec[jj].start; i != cvdPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int icv = cvdPrcommVec[ii].packVec[i];
	    MiscUtils::matVecMult(mrs->PACK_BUF+pack_size+3*i,cvdPrcommVec[ii].transformVec[jj].R,s[icv]);
	  }
	}
	else if (cvdPrcommVec[ii].transformVec[jj].has_t) {
	  for (int i = cvdPrcommVec[ii].transformVec[jj].start; i != cvdPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int icv = cvdPrcommVec[ii].packVec[i];
	    MiscUtils::vecVecAdd(mrs->PACK_BUF+pack_size+3*i,cvdPrcommVec[ii].transformVec[jj].t,s[icv]);
	  }
	}
	else {
	  for (int i = cvdPrcommVec[ii].transformVec[jj].start; i != cvdPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int icv = cvdPrcommVec[ii].packVec[i];
	    FOR_J3 mrs->PACK_BUF[pack_size+3*i+j] = s[icv][j];
	  }
	}
      }
      break;

    case REPLACE_DATA:

      for (int jj = 0, limit = cvdPrcommVec[ii].transformVec.size(); jj < limit; ++jj) {
        for (int i = cvdPrcommVec[ii].transformVec[jj].start; i != cvdPrcommVec[ii].transformVec[jj].end; ++i) {
          const int icv = cvdPrcommVec[ii].packVec[i];
          FOR_J3 mrs->PACK_BUF[pack_size+3*i+j] = s[icv][j];
        }
      }

      break;
    default:
      assert(0);
    }
    
    MPI_Issend(mrs->PACK_BUF+pack_size,cvdPrcommVec[ii].packVec.size()*3,
	       MPI_T,cvdPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += cvdPrcommVec[ii].packVec.size()*3;

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateCvdDataSeparateGhostsFinish(T (*s)[3],T (*sg)[3]) {

  map<const void*,MpiRequestStuff*>::iterator iter = mpiRequestMap.find((void*)s);
  if (iter == mpiRequestMap.end())
    CERR("updateCvdDataSeparateGhostsFinish(T*): s/sg is not mapped.");
  
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

void updateCvdDataSeparateGhosts(T (*s)[3],T (*sg)[3],const int action = REPLACE_ROTATE_DATA) {
  updateCvdDataSeparateGhostsStart(s,sg,action);
  updateCvdDataSeparateGhostsFinish(s,sg);
}

void updateCvdDataStart(T (*s)[3],const int action = REPLACE_ROTATE_DATA) {
  updateCvdDataSeparateGhostsStart(s,s+ncv,action);
}

void updateCvdDataFinish(T (*s)[3]) {
  updateCvdDataSeparateGhostsFinish(s,s+ncv);
}

void updateCvdData(T (*s)[3],const int action = REPLACE_ROTATE_DATA) {
  updateCvdDataStart(s,action);
  updateCvdDataFinish(s);
}
