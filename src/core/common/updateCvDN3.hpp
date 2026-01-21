// =====================================================================
// these routines use the mpi_comm to exchange data... 
// =====================================================================

void updateCvDataStart(T (*s)[3],const int action = REPLACE_ROTATE_DATA) {

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
    //unpack_size += cvPrcommVec[ii].unpackVec.size()*3;
    pack_size += cvPrcommVec[ii].packVec.size()*3;
  }

  // and pack and send...
    
  //assert(mrs->UNPACK_BUF == NULL); mrs->UNPACK_BUF = new T[unpack_size];
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
    
  //unpack_size = 0;
  pack_size = 0;
  for (int ii = 0; ii < cvPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv(s+cvPrcommVec[ii].unpack_offset,cvPrcommVec[ii].unpack_size*3,
	      MPI_T,cvPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
    
    // pack - done below now...
    
    // any periodic transformations are applied on the pack side...
    switch (action) {
    case REPLACE_ROTATE_DATA: // default transformation
      for (int jj = 0, limit = cvPrcommVec[ii].transformVec.size(); jj < limit; ++jj) {
	// never translate...
	if (cvPrcommVec[ii].transformVec[jj].has_R) {
	  for (int i = cvPrcommVec[ii].transformVec[jj].start; i != cvPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int icv = cvPrcommVec[ii].packVec[i];
	    MiscUtils::matVecMult(mrs->PACK_BUF+pack_size+3*i,cvPrcommVec[ii].transformVec[jj].R,s[icv]);
	  }
	}
	else {
	  for (int i = cvPrcommVec[ii].transformVec[jj].start; i != cvPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int icv = cvPrcommVec[ii].packVec[i];
	    FOR_J3 mrs->PACK_BUF[pack_size+3*i+j] = s[icv][j];
	  }
	}
      }
      break;
    case REPLACE_TRANSLATE_DATA:
      for (int jj = 0, limit = cvPrcommVec[ii].transformVec.size(); jj < limit; ++jj) {
        if ((cvPrcommVec[ii].transformVec[jj].has_R)&&(cvPrcommVec[ii].transformVec[jj].has_t)) {
	  for (int i = cvPrcommVec[ii].transformVec[jj].start; i != cvPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int icv = cvPrcommVec[ii].packVec[i];
            double tmp[3];
	    MiscUtils::matVecMult(tmp,cvPrcommVec[ii].transformVec[jj].R,s[icv]);
	    MiscUtils::vecVecAdd(mrs->PACK_BUF+pack_size+3*i,cvPrcommVec[ii].transformVec[jj].t,tmp);
	  }
        }
        else if (cvPrcommVec[ii].transformVec[jj].has_R) {
	  for (int i = cvPrcommVec[ii].transformVec[jj].start; i != cvPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int icv = cvPrcommVec[ii].packVec[i];
	    MiscUtils::matVecMult(mrs->PACK_BUF+pack_size+3*i,cvPrcommVec[ii].transformVec[jj].R,s[icv]);
	  }
	}
	else if (cvPrcommVec[ii].transformVec[jj].has_t) {
	  for (int i = cvPrcommVec[ii].transformVec[jj].start; i != cvPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int icv = cvPrcommVec[ii].packVec[i];
	    MiscUtils::vecVecAdd(mrs->PACK_BUF+pack_size+3*i,cvPrcommVec[ii].transformVec[jj].t,s[icv]);
	  }
	}
	else {
	  for (int i = cvPrcommVec[ii].transformVec[jj].start; i != cvPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int icv = cvPrcommVec[ii].packVec[i];
	    FOR_J3 mrs->PACK_BUF[pack_size+3*i+j] = s[icv][j];
	  }
	}
      }
      break;

    case REPLACE_DATA:

      for (int jj = 0, limit = cvPrcommVec[ii].transformVec.size(); jj < limit; ++jj) {
        for (int i = cvPrcommVec[ii].transformVec[jj].start; i != cvPrcommVec[ii].transformVec[jj].end; ++i) {
          const int icv = cvPrcommVec[ii].packVec[i];
          FOR_J3 mrs->PACK_BUF[pack_size+3*i+j] = s[icv][j];
        }
      }

      break;
    default:
      assert(0);
    }
    
    // send...
    MPI_Issend(mrs->PACK_BUF+pack_size,cvPrcommVec[ii].packVec.size()*3,
	       MPI_T,cvPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += cvPrcommVec[ii].packVec.size()*3;

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateCvDataFinish(T (*s)[3]) {

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

void updateCvData(T (*s)[3],const int action = REPLACE_ROTATE_DATA) {
  updateCvDataStart(s,action);
  updateCvDataFinish(s);
}

