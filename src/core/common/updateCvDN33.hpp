// =====================================================================
// these routines use the mpi_comm to exchange data... 
// =====================================================================

void updateCvDataStart(T (*s)[3][3],const int action = REPLACE_ROTATE_DATA) {

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
    
  int pack_size = 0;
  for (int ii = 0; ii < cvPrcommVec.size(); ++ii) {
    pack_size += cvPrcommVec[ii].packVec.size()*9;
  }

  // and pack and send...
    
  assert(mrs->PACK_BUF == NULL);   mrs->PACK_BUF   = new T[pack_size];    
    
  pack_size = 0;
  for (int ii = 0; ii < cvPrcommVec.size(); ++ii) {
      
    // post irecv...
    MPI_Irecv(s+cvPrcommVec[ii].unpack_offset,cvPrcommVec[ii].unpack_size*9,
	      MPI_T,cvPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->recvRequestArray[ii]));
    
    // any periodic transformations are applied on the pack side...
    switch (action) {
    case REPLACE_ROTATE_DATA: // default transformation
    case REPLACE_TRANSLATE_DATA: // for the tensor, we're not applying the vector translation..
      for (int jj = 0, limit = cvPrcommVec[ii].transformVec.size(); jj < limit; ++jj) {
	// never translate...
	if (cvPrcommVec[ii].transformVec[jj].has_R) {
          
          // recall that the R is stored in row-major format ... 
          
          double tmp[9];
          
          for (int i = cvPrcommVec[ii].transformVec[jj].start; i != cvPrcommVec[ii].transformVec[jj].end; ++i) {
           
            const int icv = cvPrcommVec[ii].packVec[i];

            for (int kk = 0; kk < 9; ++kk) 
              mrs->PACK_BUF[pack_size+9*i+kk] = 0.0;

            // compute T R^T 

            for (int j =0; j < 3; ++j) { 
              for (int k = 0; k < 3; ++k) { 
                tmp[3*j+k] = 0.0;
                for (int m = 0; m < 3; ++m) 
                  tmp[3*j+k] += s[icv][j][m] * cvPrcommVec[ii].transformVec[jj].R[3*k+m];
              }
            }

            // now pack R (TR^T)

            for (int j = 0; j < 3; ++j) 
              for (int k =0; k < 3; ++k) 
                for (int m =0; m < 3; ++m) 
                  mrs->PACK_BUF[pack_size+9*i+3*j+k] += cvPrcommVec[ii].transformVec[jj].R[3*j+m]*tmp[3*m+k];
          }
        }
        else {
	  for (int i = cvPrcommVec[ii].transformVec[jj].start; i != cvPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int icv = cvPrcommVec[ii].packVec[i];
            for (int j = 0; j < 3; ++j) { 
              for (int k =0; k < 3; ++k) { 
                mrs->PACK_BUF[pack_size+9*i+3*j+k] = s[icv][j][k];
              }
            }
	  }
	}
      }
      break;
    case REPLACE_DATA:

      for (int jj = 0, limit = cvPrcommVec[ii].transformVec.size(); jj < limit; ++jj) {
        for (int i = cvPrcommVec[ii].transformVec[jj].start; i != cvPrcommVec[ii].transformVec[jj].end; ++i) {
	    const int icv = cvPrcommVec[ii].packVec[i];
            for (int j = 0; j < 3; ++j) { 
              for (int k =0; k < 3; ++k) { 
                mrs->PACK_BUF[pack_size+9*i+3*j+k] = s[icv][j][k];
              }
            }
	  }
	}
      break;
    default:
      assert(0);
    }
    
    // send...
    MPI_Issend(mrs->PACK_BUF+pack_size,cvPrcommVec[ii].packVec.size()*9,
	       MPI_T,cvPrcommVec[ii].rank,UPDATE_TAG,mpi_comm,&(mrs->sendRequestArray[ii]));
    
    pack_size += cvPrcommVec[ii].packVec.size()*9;

  }
    
  // everything is packed now, and the request stuff is mapped. We can go back to doing work
  // if we want and call finish on the same pointer in a bit...
    
}  

void updateCvDataFinish(T (*s)[3][3]) {

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

void updateCvData(T (*s)[3][3],const int action = REPLACE_ROTATE_DATA) {
  updateCvDataStart(s,action);
  updateCvDataFinish(s);
}

