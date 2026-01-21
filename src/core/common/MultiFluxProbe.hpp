#ifndef _MULTI_FLUX_PROBE_HPP_
#define _MULTI_FLUX_PROBE_HPP_

enum MultiFluxProbeFormat {
  SINGLE_FILE,
  MULTI_FILE,
};

class MultiFluxProbe {
public:

  string name,header;
  int interval;
  MultiFluxProbeFormat format;
  vector<string> varNameVec; 
  int probe_rank0;
  int probe_ngr;
  int * probe_supergroup_last;
  double * probe_proj_area;
  double (*probe_x)[3];
  int probe_pack_buf_size;
  int probe_pack_index_size;
  int (*probe_pack_index)[2];  
  int probe_unpack_buf_size;
  int probe_active_rank_count;
  int * probe_raoar;
  int * probe_ioar_i;
  int * probe_ioar_v;
  std::stringstream * ss;
  vector<pair<int,int> > ssVec;
  
  MultiFluxProbe() {
    
    interval = -1;
    format = SINGLE_FILE; // default is single file
    probe_rank0 = -1;
    probe_ngr = 0;
    probe_supergroup_last = NULL;
    probe_proj_area = NULL;
    probe_x = NULL;
    probe_pack_buf_size = 0;
    probe_pack_index_size = 0;
    probe_pack_index = NULL;
    probe_unpack_buf_size = 0;
    probe_active_rank_count = 0;
    probe_raoar = NULL;
    probe_ioar_i = NULL;
    probe_ioar_v = NULL;
    ss = NULL;
     
  }

  ~MultiFluxProbe() {

    DELETE(probe_supergroup_last);
    DELETE(probe_proj_area);
    DELETE(probe_x);
    DELETE(probe_pack_index);
    DELETE(probe_raoar);
    DELETE(probe_ioar_i);
    DELETE(probe_ioar_v);
    if (ss != NULL) delete ss;

  }

  void flush() {
    
    if (probe_rank0 == mpi_rank) {
      assert(ss);
      ofstream outFile;
      char filename[128];
      if (format == SINGLE_FILE) {
	// in single-file mode, just open and append...
	sprintf(filename,"%s.mfp",name.c_str());
        if (interval == -1) {
          const string tmp_filename = MiscUtils::makeTmpPrefix(filename);
          outFile.open(tmp_filename.c_str(),ofstream::trunc);
          assert(outFile.is_open());
          outFile << ss->rdbuf();
          outFile.close();
          remove(filename);
          rename(tmp_filename.c_str(),filename);
        }
        else {
          outFile.open(filename,ofstream::app);
          assert(outFile.is_open());
          outFile << ss->rdbuf();
          outFile.close();
        }
	ss->str(string()); // this clears ss, believe it or not...
	// the ssVec should be empty...
	assert(ssVec.empty());
      }
      else {
	assert(format == MULTI_FILE);
	// in multi-file mode, write a new file for each time step...
	int size_max = 0;
	int pos1 = 0;
	for (int ii = 0; ii < ssVec.size(); ++ii) {
	  const int pos0 = pos1;
	  pos1 = ssVec[ii].second;
	  size_max = max(size_max,pos1-pos0);
	}
	char * cbuf = new char[size_max+1];
	pos1 = 0;
	for (int ii = 0; ii < ssVec.size(); ++ii) {
	  const int pos0 = pos1;
	  pos1 = ssVec[ii].second;
	  ss->rdbuf()->pubseekpos(pos0);
	  ss->rdbuf()->sgetn(cbuf,pos1-pos0);
	  cbuf[pos1-pos0] = '\0';
	  const int step = ssVec[ii].first;
	  sprintf(filename,"%s.%08d.mfp",name.c_str(),step);
          const string tmp_filename = MiscUtils::makeTmpPrefix(filename);
          outFile.open(tmp_filename.c_str(),ofstream::trunc);
          assert(outFile.is_open());
          outFile << header;
          outFile << cbuf;
          outFile.close();
          remove(filename);
          rename(tmp_filename.c_str(),filename);
	}
	delete[] cbuf;
	ss->str(string());
	ssVec.clear();
      }
    }

  }

  void doProbe(const int step = -1,const double time = 0.0) {
    
    // should have some valid vars...

    const int nvar = varNameVec.size();
    assert(nvar > 0);

    // get/evaluate the variables required for the probe. They 
    // should all be valid because they are checked when the 
    // probe is initialized in initProbe...
    
    vector<CtiRegister::CtiData*> dataVec(nvar);
    for (int ii = 0; ii < nvar; ++ii) {
      dataVec[ii] = CtiRegister::getCtiData(varNameVec[ii]);
      assert((dataVec[ii] != NULL)&&(dataVec[ii]->getTopology() == SIGNED_FA_DATA)&&(dataVec[ii]->getType() == DN_DATA));
    }

    double * pack_buf = new double[probe_pack_buf_size*nvar];
    for (int ipack = 0; ipack < probe_pack_buf_size*nvar; ++ipack)
      pack_buf[ipack] = 0.0;
      
    for (int ii = 0; ii < probe_pack_index_size; ++ii) {
      const int ipack = probe_pack_index[ii][0];
      assert((ipack >= 0)&&(ipack < probe_pack_buf_size)); 
      int ifa = probe_pack_index[ii][1];
      double fa_sign = 1.0;
      if (ifa < 0) {
	ifa = -ifa-1;
	fa_sign = -1.0;
      }
      for (int jj = 0; jj < nvar; ++jj) {
        pack_buf[ipack*nvar+jj] += fa_sign*dataVec[jj]->dn(ifa);
      }
    }

    double * unpack_buf = NULL;
    if (mpi_rank == probe_rank0)
      unpack_buf = new double[probe_ngr*nvar];
    
    probe_reduce(pack_buf,unpack_buf,nvar);
    
    delete[] pack_buf;
    
    // done reduce -- now report...
    
    if (mpi_rank == probe_rank0) {
      
      int supergroup_index = 0;
      int supergroup_ngr = 0;
      double *total_buf = new double[4+nvar];
      double *supergroup_buf = new double[4+nvar];
      for (int ii = 0; ii < 4+nvar; ++ii) {
        supergroup_buf[ii] = 0.0;
        total_buf[ii] = 0.0;
      }
      for (int igr = 0; igr < probe_ngr; ++igr) {
	// add supergroup and total bufs...
	supergroup_buf[0] += probe_x[igr][0]*probe_proj_area[igr];
	supergroup_buf[1] += probe_x[igr][1]*probe_proj_area[igr];
	supergroup_buf[2] += probe_x[igr][2]*probe_proj_area[igr];
	supergroup_buf[3] += probe_proj_area[igr];
	total_buf[0] += probe_x[igr][0]*probe_proj_area[igr];
	total_buf[1] += probe_x[igr][1]*probe_proj_area[igr];
	total_buf[2] += probe_x[igr][2]*probe_proj_area[igr];
	total_buf[3] += probe_proj_area[igr];
        for (int ivar = 0; ivar < nvar; ++ivar) {
          supergroup_buf[4+ivar] += unpack_buf[igr*nvar+ivar];
          total_buf[4+ivar] += unpack_buf[igr*nvar+ivar];
        }
	// and report...
	if (probe_supergroup_last[igr] == 1) {
	  assert(supergroup_ngr == 0);
	  // this is a probe that is the ONLY member of the current supergroup. report it separately with single indexing...
	  // unless it is the only group. In this case, we will only report the total (i.e. no indexing)...
	  if (probe_ngr > 1) {
            // if this is the last one, then report the supergroup as well...
	    *ss << "probe-" << (supergroup_index+1) << " " << step << " " << time << " " << 
	      probe_x[igr][0] << " " << probe_x[igr][1] << " " << probe_x[igr][2] << " " << probe_proj_area[igr];
            for (int ivar = 0; ivar < nvar; ++ivar) 
              *ss << " " << unpack_buf[igr*nvar+ivar]; 
	    *ss << endl;
	    // and as a sum...
	    *ss << "probe-sum-" << (supergroup_index+1) << " " << step << " " << time << " " << 
	      probe_x[igr][0] << " " << probe_x[igr][1] << " " << probe_x[igr][2] << " " << probe_proj_area[igr];
            for (int ivar = 0; ivar < nvar; ++ivar) 
              *ss << " " << unpack_buf[igr*nvar+ivar]; 
	    *ss << endl;

	  }
	  ++supergroup_index;
	  // leave supergroup_ngr = 0
          for (int ii = 0; ii < 4+nvar; ++ii) {
            supergroup_buf[ii] = 0.0;
          }
	}
	else {
	  // this is a probe that is part of the current supergroup. report it separately with double indexing...
          *ss << "probe-" << (supergroup_index+1) << "-" << igr << " " << step << " " << time << " " <<
	    probe_x[igr][0] << " " << probe_x[igr][1] << " " << probe_x[igr][2] << " " << probe_proj_area[igr];
          for (int ivar = 0; ivar < nvar; ++ivar) 
            *ss << " " << unpack_buf[igr*nvar+ivar];
          *ss << endl;
	  // if this is the last one, then report the supergroup as well...
	  if (probe_supergroup_last[igr] > 1) {
	    *ss << "probe-sum-" << (supergroup_index+1) << " " << step << " " << time << " " << 
	      supergroup_buf[0]/supergroup_buf[3] << " " << supergroup_buf[1]/supergroup_buf[3] << " " << supergroup_buf[2]/supergroup_buf[3] << " " << supergroup_buf[3];
            for (int ivar = 0; ivar < nvar; ++ivar) 
	      *ss << " " << supergroup_buf[4+ivar];
            *ss << endl;
	    ++supergroup_index;
	    supergroup_ngr = 0;
            for (int ii = 0; ii < 4+nvar; ++ii) {
              supergroup_buf[ii] = 0.0;
            }
	  }
	  else {
	    ++supergroup_ngr;
	  }
	}
      }
      // and finally report the total...
      if (probe_ngr == 1) {
	*ss << "probe " << step << " " << time << " " << 
	  total_buf[0]/total_buf[3] << " " << 
	  total_buf[1]/total_buf[3] << " " << 
	  total_buf[2]/total_buf[3] << " " << 
	  total_buf[3];
      }
      else {
	*ss << "probe-total " << step << " " << time << " " << 
	  total_buf[0]/total_buf[3] << " " << 
	  total_buf[1]/total_buf[3] << " " << 
	  total_buf[2]/total_buf[3] << " " << 
	  total_buf[3];
      }
      for (int ivar = 0; ivar < nvar; ++ivar) {
        *ss << " " << total_buf[4+ivar];
      }
      *ss << endl;
      delete[] unpack_buf;
      delete[] total_buf;
      delete[] supergroup_buf;
      
      //cout << ss.rdbuf();

      if (format == MULTI_FILE) {
	ss->seekp(0, ios::end);
	ssVec.push_back(pair<int,int>(step,ss->tellp()));
      }

    }
    
    //MPI_Barrier(mpi_comm);
    
  }

  // formatting example...
  /*
   *ss << setfill(' ') << setw(17) << "total-probe " << "     " << " ";
   *ss << setfill(' ') << setw(8) << step << " " << setfill(' ') << setw(16) << time << " ";
   FOR_I3 *ss << setfill(' ') << setw(16) << total_buf[i]/total_buf[3] << " ";
   *ss << setfill(' ') << setw(16) << total_buf[3];
   for (int ivar = 0; ivar < nvar; ++ivar) 
   *ss << setfill(' ') << setw(16) << total_buf[4+ivar] << " ";
   *ss << endl;
   */
  
  void probe_reduce(double * pack_buf,double * unpack_buf,const int nd) {
      
    if (mpi_rank == probe_rank0) {
	
      for (int i = 0; i < probe_ngr*nd; ++i)
	unpack_buf[i] = 0.0;
	
      // pull data from active ranks...
	
      double * dbuf = new double[probe_unpack_buf_size*nd];
      for (int iar = 0; iar < probe_active_rank_count; ++iar) {
	const int rank = probe_raoar[iar]; assert((rank >= 0)&&(rank < mpi_size));
	if (rank == mpi_rank) {
	  for (int ioa = probe_ioar_i[iar]; ioa != probe_ioar_i[iar+1]; ++ioa) {
	    const int igr = probe_ioar_v[ioa]; assert((igr >= 0)&&(igr < probe_ngr));
	    for (int id = 0; id < nd; ++id) 
	      unpack_buf[igr*nd+id] += pack_buf[(ioa-probe_ioar_i[iar])*nd+id];
	  }
	}
	else {
	  MPI_Recv(dbuf,(probe_ioar_i[iar+1]-probe_ioar_i[iar])*nd,MPI_DOUBLE,rank,42426,mpi_comm,MPI_STATUS_IGNORE);
	  for (int ioa = probe_ioar_i[iar]; ioa != probe_ioar_i[iar+1]; ++ioa) {
	    const int igr = probe_ioar_v[ioa]; assert((igr >= 0)&&(igr < probe_ngr));
	    for (int id = 0; id < nd; ++id) 
	      unpack_buf[igr*nd+id] += dbuf[(ioa-probe_ioar_i[iar])*nd+id];
	  }
	}
      }
      delete[] dbuf;
	
    }
    else if (probe_pack_buf_size > 0) {
	
      MPI_Ssend(pack_buf,probe_pack_buf_size*nd,MPI_DOUBLE,probe_rank0,42426,mpi_comm);
	
    }
      
  }
    
};

#endif
