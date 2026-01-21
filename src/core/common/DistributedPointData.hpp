#ifndef DISTRIBUTEDPOINTDATA_HPP
#define DISTRIBUTEDPOINTDATA_HPP

#include "MiscUtils.hpp" // for dumpRange

// ======================================================================
// DistributedPointData.hpp
// 
// a parallel class for point data distributed in a striped fashion across
// processors.
//
// History:
//  author: Frank Ham - April 2013
//  May 2013 - scalar/vector decision made upfront to avoid large mem seg faults
// ======================================================================

// ======================================================================
// Copyright (c) 2013 Cascade Technologies Inc.
// All Rights Reserved. www.cascadetechnologies.com
//
// This source is proprietary software and cannot be used, redistributed,
// or modified except under the terms of Cascade's Software License 
// Agreement. 
//
// THIS SOURCE CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY 
// OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO 
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
// ======================================================================

typedef double *double_ptr;
typedef double (*double3_ptr)[3];

class PointData {
public:
  int n;
  double (*x)[3];
  vector< std::pair<string,double_ptr> > doubleScalarData;
  vector< std::pair<string,double3_ptr> > doubleVectorData;
  PointData() {
    x = NULL;
  }
  ~PointData() {
    DELETE(x);
    for (int iscalar = 0; iscalar < doubleScalarData.size(); ++iscalar)
      delete[] doubleScalarData[iscalar].second;
    for (int ivector = 0; ivector < doubleVectorData.size(); ++ivector)
      delete[] doubleVectorData[ivector].second;
  }

  void dumpRanges() {
    dumpRange(x,n,"X");
    for (int iscalar = 0; iscalar < doubleScalarData.size(); ++iscalar)
      dumpRange(doubleScalarData[iscalar].second,n,doubleScalarData[iscalar].first);
    for (int ivector = 0; ivector < doubleVectorData.size(); ++ivector)
      dumpRange(doubleVectorData[ivector].second,n,doubleVectorData[ivector].first);
  }

  void rename(const string& old_name,const string& new_name) {
    for (int iscalar = 0; iscalar < doubleScalarData.size(); ++iscalar)
      if (doubleScalarData[iscalar].first == old_name)
	doubleScalarData[iscalar].first = new_name;
    for (int ivector = 0; ivector < doubleVectorData.size(); ++ivector)
      if (doubleVectorData[ivector].first == old_name)
	doubleVectorData[ivector].first = new_name;
  }
 
};

class DistributedPointData : public PointData {
public:
  
  int * xora;
  
  DistributedPointData() {
    xora = NULL;
  }
  
  ~DistributedPointData() {
    DELETE(xora);
  }
  
  void readIpFile(const string& filename);
 
};


class AsciiPointData : public PointData { 

public: 
  int * xora ; 

  AsciiPointData() : xora(NULL) {} 
  ~AsciiPointData() { 
    DELETE(xora) ; 
  } 


  void readAsciiFile(const string& filename) ; 

}; 

#include "FileParser.hpp"

// define this to skip the data reading and just populate all scalars
// with 2, and all vectors with {3,4,5}...
//#define DEBUG_DPD

int safe_unmap(map<int,int>& m, int key) { 

  assert ( m.find(key) != m.end()) ; 
  return m[key] ; 
} 

void AsciiPointData::readAsciiFile(const string& filename) { 

  // 
  // first two lines of the data file are 
  // 
  // <n-global> 
  // <var names> 
  // the variables must be written as VECTOR-<name> and 
  // SCALAR-<name> so that we can distinguish between 
  // scalar/vector data. 
  // 

  int nglobal ; 
  FileParser * fp = NULL; 
  vector<string> varName ; 


  if ( mpi_rank == 0 ) { 
    cout << " > reading : " << filename << " ... " << endl ; 

    fp = new FileParser(filename, " \t\n", "[]"); 
    string token ; 
    FileParserTokenType fptt ; 

    fptt    = fp->getNextToken(token,FILEPARSER_STANDARD_TOKEN) ;
    nglobal = fp->stringToInt(token); 
    cout << " > nglobal = " << nglobal << endl ; 

    fp->getNextToken(token,FILEPARSER_DELIMITER_TOKEN) ; 

    // get the variables ..
    while ( (fp->getNextToken(token)) != FILEPARSER_DELIMITER_TOKEN) { 
      cout << " variable : " << token << endl ; 
      varName.push_back(token); 
    }
  }


  // share the header info with everyone .. 
  MPI_Bcast(&nglobal,1,MPI_INT,0,mpi_comm) ; 
  
  MPI_Bcast_stringVec(varName,0,mpi_comm); 
  const int nvar = varName.size() ; 

  map<int,int> varType ; // 0 for scalar, 1 for vector ... 
  map<int,int> varIndex ;

  int nscalar = 0 ; 
  int nvector = 0 ; 
  int ivar    = 0 ; 

  int var_pos = 0 ; 
  for (int ivar = 0 ; ivar < nvar ; ++ivar ) { 
    
    if ( varName[ivar].compare(0,7,"VECTOR-") == 0)  { 
      const string vec_name = varName[ivar].substr(7); 
      varName[ivar] = vec_name ;
      varType[var_pos]  = 1;
      varIndex[var_pos] = nvector++;
      var_pos += 3 ; 
    }
    else { 
      assert ( varName[ivar].compare(0,7,"SCALAR-") == 0) ; 
      const string vec_name = varName[ivar].substr(7); 
      varName[ivar]    = vec_name ; 
      varType[var_pos] = 0; 
      varIndex[var_pos] = nscalar++;
      var_pos += 1; 
    }
  }

  assert (xora == NULL); buildUniformXora(xora,nglobal); 
  assert ( x == NULL ) ; 

  //
  // allocate
  //
  n = xora[mpi_rank+1] - xora[mpi_rank]; 
  x = new double[n][3]; 
  doubleScalarData.resize(nscalar); 
  doubleVectorData.resize(nvector); 

  var_pos = 0 ; 
  int iscalar = 0 ; 
  int ivector = 0 ; 
  for (int ivar =0; ivar < nvar ; ++ivar) { 
    int vt = safe_unmap(varType,var_pos); 
    if ( vt == 0 ) {
      assert ( safe_unmap(varIndex,var_pos) == iscalar ) ; 
      doubleScalarData[iscalar].first  = varName[ivar] ; 
      doubleScalarData[iscalar].second = new double[n]; 
      ++iscalar ; 
      var_pos += 1 ;
    }
    else if ( vt == 1 ) { 
      assert ( safe_unmap(varIndex,var_pos) == ivector ) ; 
      doubleVectorData[ivector].first  = varName[ivar] ; 
      doubleVectorData[ivector].second = new double[n][3]; 
      ++ivector ; 
      var_pos += 3 ; 
    }
    else { 
      assert(0); 
    }
  }//ivar 


  assert ( iscalar == nscalar ) ; 
  assert ( ivector == nvector ) ; 
  varName.clear() ; 

  const int len_line = nscalar + 3*nvector ; 

  double * tmp_buf = new double[n*len_line] ; 

  int ierr = 0 ; 
  if ( mpi_rank == 0 ) { 

    int n_max = 0 ; 
    FOR_RANK n_max = max(n_max,xora[rank+1]-xora[rank]) ; 
    double * dbuf = new double[n_max*len_line]; 
    string token ; 
    FileParserTokenType fptt ; 

    FOR_RANK { 
      //cout << rank << "   " << xora[rank] << "    " << xora[rank+1] << endl ;
      for (int ii = xora[rank]; ii != xora[rank+1] ; ++ii) { 
        const int ioffset = len_line * ( ii - xora[rank] ) ; 
        for (int j = 0 ; j < len_line ; ++j ) { 
          fptt = fp->getNextToken(token,FILEPARSER_STANDARD_TOKEN) ; 
          dbuf[ioffset+j] = fp->stringToDouble(token); 
        }
        //assert ( fp->getNextToken(token) == FILEPARSER_DELIMITER_TOKEN) ; 
      }//ii

      if ( rank == 0 ) { 
        for (int j = 0 ; j < n*len_line; ++j) 
          tmp_buf[j] = dbuf[j]; 
      } 
      else {
        const int nsend = (xora[rank+1] - xora[rank])*len_line; 
        MPI_Send(dbuf,nsend,MPI_DOUBLE,rank,12345, mpi_comm); 
      } 
    }

    delete[] dbuf ; 
  }
  else { 
    // post the receive ..
    MPI_Status status ; 
    MPI_Recv(tmp_buf, (xora[mpi_rank+1]-xora[mpi_rank])*len_line,MPI_DOUBLE,0,12345,mpi_comm, &status); 
  } 


  // you now have data that you need to put into the appropriate buffers .. 
  for (int ii = xora[mpi_rank] ; ii != xora[mpi_rank+1]; ++ii) { 
    const int ioffset = len_line * ( ii - xora[mpi_rank] ) ; 
    int var_pos = 0 ; 
    int vt, vi ; // varType, varIndex .. 
    while ( var_pos < len_line ) { 
      vt = safe_unmap(varType,var_pos); 
      if ( vt == 0 ) { 
        vi = safe_unmap(varIndex,var_pos) ; 
        doubleScalarData[vi].second[ii-xora[mpi_rank]] = tmp_buf[ioffset+var_pos] ; 
        var_pos += 1; 
      } 
      else { 
        assert ( vt == 1) ; 
        vi = safe_unmap(varIndex,var_pos); 
        FOR_K3 
          doubleVectorData[vi].second[ii-xora[mpi_rank]][k] = tmp_buf[ioffset+var_pos+k]; 
        var_pos += 3 ; 
      }
    }
  }//ii 


  // x should be in the list of vectors.. we'll set that now .. 
  ivar = 0 ; 
  for (ivar = 0 ; ivar < doubleVectorData.size() ; ++ivar) { 
    if ( doubleVectorData[ivar].first == "X") { 
      for (int ii = xora[mpi_rank]; ii != xora[mpi_rank+1]; ++ii ) { 
        const int i = ii-xora[mpi_rank]; 
        FOR_J3 x[i][j] = doubleVectorData[ivar].second[i][j] ; 
      }//ii 
      break; 
    }//X
  }//ivar

  assert ( ivar != doubleVectorData.size() ); 

  delete[] tmp_buf ; 
  if ( fp != NULL ) 
    delete fp ; 

  MPI_Barrier(mpi_comm) ; 


} 


void DistributedPointData::readIpFile(const string& filename) {
  
  int buf[2];
  FileParser * fp;
  vector<string> varName;

  if (mpi_rank == 0) {
    
    cout << " > reading \"" << filename << "\"..." << endl;
    
    try {

      fp = new FileParser(filename," \t","()"); // filename, whitespace, braces
      string token; FileParserTokenType fptt;
      
      // some junk at the start: 3 3 ?
      fptt = fp->getNextToken(token,FILEPARSER_STANDARD_TOKEN);
      fptt = fp->getNextToken(token,FILEPARSER_STANDARD_TOKEN);
      
      // the global size...
      fptt = fp->getNextToken(token,FILEPARSER_STANDARD_TOKEN); 
      buf[0] = fp->stringToInt(token);
      cout << " > nglobal=" << buf[0] << endl;
      
      // the number of vars...
      fptt = fp->getNextToken(token,FILEPARSER_STANDARD_TOKEN); 
      buf[1] = fp->stringToInt(token);
      cout << " > nvar=" << buf[1] << endl;
      
      // the var names...
      varName.resize(buf[1]);
      for (int ivar = 0; ivar < buf[1]; ++ivar) {
	fptt = fp->getNextToken(varName[ivar],FILEPARSER_STANDARD_TOKEN);
	cout << " > variable=\"" << varName[ivar] << "\"" << endl;
      }

    }
    catch(int e) {

      // use nglobal to indicate an error...
      
      buf[0] = -1;
      
    }
    
  }
  
  // share the header with everyone...
  
  MPI_Bcast(buf, 2, MPI_INT, 0, mpi_comm);
  const int nglobal = buf[0];
  const int nvar = buf[1];

  if (nglobal == -1)
    CERR("Parsing problem in IP file.");
  
  // this routine bcast's a vector of strings to everyone...
 
  MPI_Bcast_stringVec(varName,0,mpi_comm);
  assert(varName.size() == nvar);
  
  // based on the names, decide if any of the data should be colapsed into vectors...
  
  vector<int> varType(nvar); // 0 = scalar, 1 = vector
  vector<int> varIndex(nvar); // scalar or vector index, depending on type
  vector<int> varIndex2(nvar); // second index for vectors: i.e. 0,1,2
  
  int nscalar = 0;
  int nvector = 0;
  int ivar = 0;
  while (ivar < varName.size()) {
    bool got_vector = false;
    if ((varName[ivar].compare(0,2,"x-") == 0)&&(ivar < varName.size()-2)) {
      const string vec_name = varName[ivar].substr(2);
      if ( (varName[ivar+1].compare(0,2,"y-") == 0) && (vec_name == varName[ivar+1].substr(2)) ) {
	if ( (varName[ivar+2].compare(0,2,"z-") == 0) && (vec_name == varName[ivar+2].substr(2)) ) {
	  // convert these ones to a vector...
	  if (mpi_rank == 0)
	    cout << " > interpreting scalar data \"" << varName[ivar] << "\", \"" << 
	      varName[ivar+1] << "\" and \"" <<  varName[ivar+2] << "\" as vector data \"" << 
	      vec_name << "\"" << endl;
	  got_vector = true;
	  varType[ivar]   = 1; varIndex[ivar]   = nvector; varIndex2[ivar]   = 0; 
	  varType[ivar+1] = 1; varIndex[ivar+1] = nvector; varIndex2[ivar+1] = 1;
	  varType[ivar+2] = 1; varIndex[ivar+2] = nvector; varIndex2[ivar+2] = 2;
	  ++nvector;
	  ivar += 3;
	}
      }
    }
    if (!got_vector) {
      varType[ivar] = 0; // scalar
      varIndex[ivar] = nscalar;
      varIndex2[ivar] = -1;
      ++nscalar;
      ++ivar;
    }
  }
  assert(ivar == nvar);
  assert(nscalar + 3*nvector == nvar);

  // everyone can now build xora...
  
  assert(xora == NULL); buildUniformXora(xora,nglobal);
  
  // everyone can now compute and allocate their local data size...
  
  n = xora[mpi_rank+1] - xora[mpi_rank]; // local n
  assert(x == NULL);
  x = new double[n][3];
  doubleScalarData.resize(nscalar);
  doubleVectorData.resize(nvector);
  int iscalar = 0;
  int ivector = 0;
  for (int ivar = 0; ivar < nvar; ++ivar) {
    if (varType[ivar] == 0) {
      assert(varIndex[ivar] == iscalar);
      doubleScalarData[iscalar].first = varName[ivar];
      doubleScalarData[iscalar].second = new double[n];
      ++iscalar;
    }
    else if (varType[ivar] == 1) {
      if (varIndex2[ivar] == 0) {
	assert(varIndex[ivar] == ivector);
	doubleVectorData[ivector].first = varName[ivar].substr(2); // strip "x-"
	doubleVectorData[ivector].second = new double[n][3];
	++ivector;
      }
    }
  }
  assert(iscalar == nscalar);
  assert(ivector == nvector);
  varName.clear();

  // parse the data section with error checking...
  
  int ierr = 0;
  if (mpi_rank == 0) {
    
    // rank 0 reads and sends...
    
    int n_max = 0; FOR_RANK n_max = max(n_max,xora[rank+1]-xora[rank]);
    double * dbuf = new double[n_max];
    string token; FileParserTokenType fptt;
    
    // coordinate data...

    FOR_J3 {

      try {
	if (ierr == 0) {
	  cout << " > x[" << j << "]";
#ifndef DEBUG_DPD
	  fptt = fp->getNextToken(token,FILEPARSER_DELIMITER_TOKEN); 
	  assert(token == "(");
#endif
	}
      }
      catch (int e) {
	ierr = -1;
      }

      FOR_RANK {
	try {
	  if (ierr == 0) {
	    for (int i = xora[rank]; i != xora[rank+1]; ++i) {
	      if (i%(nglobal/10)==0) {
		cout << ".";
		cout.flush();
	      }
#ifndef DEBUG_DPD
	      fptt = fp->getNextToken(token,FILEPARSER_STANDARD_TOKEN); 
	      dbuf[i-xora[rank]] = fp->stringToDouble(token);
#else
	      dbuf[i-xora[rank]] = 1.0;
#endif
	    }
	  }
	}
	catch (int e) {
	  ierr = -1;
	}
	if (rank == 0) {
	  for (int i = 0; i < n; ++i) {
	    x[i][j] = dbuf[i];
	  }
	}
	else {
	  // only send if this n > 0...
	  if (xora[rank+1]-xora[rank] > 0)
	    MPI_Ssend(dbuf,xora[rank+1]-xora[rank],MPI_DOUBLE,rank,12345,mpi_comm);
	}
      }
      
      try {
	if (ierr == 0) {
#ifndef DEBUG_DPD
	  fptt = fp->getNextToken(token,FILEPARSER_DELIMITER_TOKEN);
	  assert(token == ")");
#endif
	  cout << "OK" << endl;
      	}
      }
      catch (int e) {
	ierr = -1;
      }
    
    }

    // =========================================
    // named data...
    // =========================================
    
    for (int ivar = 0; ivar < nvar; ++ivar) {

      if (varType[ivar] == 0) {
	
	const int iscalar = varIndex[ivar];

	try {
	  if (ierr == 0) {
	    cout << " > " << doubleScalarData[iscalar].first;
#ifndef DEBUG_DPD
	    fptt = fp->getNextToken(token,FILEPARSER_DELIMITER_TOKEN);
	    assert(token == "(");
#endif
	  }
	}
	catch (int e) {
	  ierr = -1;
	}
	
	FOR_RANK {
	  try {
	    if (ierr == 0) {
	      for (int i = xora[rank]; i != xora[rank+1]; ++i) {
		if (i%(nglobal/10)==0) {
		  cout << ".";
		  cout.flush();
		}
#ifndef DEBUG_DPD
		fptt = fp->getNextToken(token,FILEPARSER_STANDARD_TOKEN); 
		dbuf[i-xora[rank]] = fp->stringToDouble(token);
#else
		dbuf[i-xora[rank]] = 2.0;
#endif	     
	      }
	    }
	  }
	  catch (int e) {
	    ierr = -1;
	  }
	  if (rank == 0) {
	    for (int i = 0; i < n; ++i) {
	      doubleScalarData[iscalar].second[i] = dbuf[i];
	    }
	  }
	  else {
	    if (xora[rank+1]-xora[rank] > 0)
	      MPI_Ssend(dbuf,xora[rank+1]-xora[rank],MPI_DOUBLE,rank,12346,mpi_comm);
	  }
	  
	}
      
	try {
	  if (ierr == 0) {
#ifndef DEBUG_DPD
	    fptt = fp->getNextToken(token,FILEPARSER_DELIMITER_TOKEN);
	    assert(token == ")");
#endif
	    cout << "OK" << endl;
	  }
	}
	catch (int e) {
	  ierr = -1;
	}
      
      }
      else {
	
	assert(varType[ivar] == 1);
	const int ivector = varIndex[ivar];
	const int jvector = varIndex2[ivar];
	
	try {
	  if (ierr == 0) {
	    cout << " > " << doubleVectorData[ivector].first << "[" << jvector << "]";
#ifndef DEBUG_DPD
	    fptt = fp->getNextToken(token,FILEPARSER_DELIMITER_TOKEN);
	    assert(token == "(");
#endif
	  }
	}
	catch (int e) {
	  ierr = -1;
	}
	
	FOR_RANK {
	  try {
	    if (ierr == 0) {
	      for (int i = xora[rank]; i != xora[rank+1]; ++i) {
		if (i%(nglobal/10)==0) {
		  cout << ".";
		  cout.flush();
		}
#ifndef DEBUG_DPD
		fptt = fp->getNextToken(token,FILEPARSER_STANDARD_TOKEN); 
		dbuf[i-xora[rank]] = fp->stringToDouble(token);
#else
		dbuf[i-xora[rank]] = 3.0 + double(jvector);
#endif	     
	      }
	    }
	  }
	  catch (int e) {
	    ierr = -1;
	  }
	  if (rank == 0) {
	    for (int i = 0; i < n; ++i) {
	      doubleVectorData[ivector].second[i][jvector] = dbuf[i];
	    }
	  }
	  else {
	    if (xora[rank+1]-xora[rank] > 0)
	      MPI_Ssend(dbuf,xora[rank+1]-xora[rank],MPI_DOUBLE,rank,12347,mpi_comm);
	  }
	  
	}
      
	try {
	  if (ierr == 0) {
#ifndef DEBUG_DPD
	    fptt = fp->getNextToken(token,FILEPARSER_DELIMITER_TOKEN);
	    assert(token == ")");
#endif
	    cout << "OK" << endl;
	  }
	}
	catch (int e) {
	  ierr = -1;
	}
      
      }
      
    }
    
    delete[] dbuf;
    delete fp;
    
  }
  else {
    
    // other ranks recv...
    
    if (n > 0) {
      
      MPI_Status status;
      
      // coordinate data....
      
      double * dbuf = new double[n];
      FOR_J3 {
	MPI_Recv(dbuf,n,MPI_DOUBLE,0,12345,mpi_comm,&status);
	for (int i = 0; i < n; ++i) {
	  x[i][j] = dbuf[i];
	}
      }

      // named data...
      
      for (int ivar = 0; ivar < nvar; ++ivar) {
	
	if (varType[ivar] == 0) {
	  const int iscalar = varIndex[ivar];
	  MPI_Recv(doubleScalarData[iscalar].second,n,MPI_DOUBLE,0,12346,mpi_comm,&status);
	}
	else {
	  assert(varType[ivar] == 1);
	  const int ivector = varIndex[ivar];
	  const int jvector = varIndex2[ivar];
	  MPI_Recv(dbuf,n,MPI_DOUBLE,0,12347,mpi_comm,&status);
	  for (int i = 0; i < n; ++i) {
	    doubleVectorData[ivector].second[i][jvector] = dbuf[i];
	  }
	}
	
      }
      
      delete[] dbuf;
      
    }
    
  }
  
  // exchange ierr from mpi_rank==0 process...
  
  MPI_Bcast(&ierr, 1, MPI_INT, 0, mpi_comm);
  if (ierr == -1)
    CERR("ip file data problem");
  
  if (mpi_rank == 0)
    cout << " > done." << endl;
  
}

#endif
