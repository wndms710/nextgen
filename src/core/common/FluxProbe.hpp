#ifndef FLUXPROBE_HPP
#define FLUXPROBE_HPP

class FluxProbe {
public: 

  string name;
  int interval;
  vector<string> var_vec;
  int probe_rank0;
  stringstream * ss;
  vector<int> face_vec;

  double probe_x[3];
  double projected_area;

  FluxProbe() { 

    name        = "";
    interval    = -1;
    probe_rank0 = -1;
    ss          = NULL;

    projected_area = 0.0;
    for (int i =0; i < 3; ++i)
      probe_x[i] = 0.0;
  }

  ~FluxProbe() {
    if ( ss != NULL) delete ss;
  }

  void doProbe(const int step = -1, const double time = 0.0) { 

    const int nvar = var_vec.size();
    assert( nvar > 0);

    vector<CtiRegister::CtiData*> data_vec(nvar);
    for (int ii = 0; ii < nvar; ++ii) { 
      data_vec[ii] = CtiRegister::getCtiData(var_vec[ii]);
      assert( (data_vec[ii] != NULL) && (data_vec[ii]->getTopology() == SIGNED_FA_DATA) &&
              (data_vec[ii]->getType() == DN_DATA));
    }

    double * pack_buf = new double[nvar];
    for (int ivar = 0; ivar < nvar; ++ivar) 
      pack_buf[ivar] = 0.0;

    const int nff = face_vec.size();
    for (int ii = 0; ii < nff; ++ii) { 

      int ifa = face_vec[ii];
      double fa_sign = 1.0;
      if ( ifa < 0) { 
        ifa = -ifa-1;
        fa_sign = -1.0;
      }

      for (int ivar = 0; ivar < nvar; ++ivar) 
        pack_buf[ivar] += fa_sign*data_vec[ivar]->dn(ifa);

    }

    double * unpack_buf = NULL;
    if ( mpi_rank == probe_rank0) 
      unpack_buf = new double[nvar];

    probe_reduce(pack_buf,unpack_buf,nvar);

    delete[] pack_buf;

    if ( mpi_rank == probe_rank0) { 
      
      *ss << "probe " << step << " " << time << " " << 
        probe_x[0] << " " << 
        probe_x[1] << " " << 
        probe_x[2] << " " << 
        projected_area;

      for (int ivar = 0; ivar < nvar; ++ivar) 
        *ss << "  " << unpack_buf[ivar];
      *ss << endl;

      delete[] unpack_buf;
    }
  }


  void flush() {

    if ( mpi_rank == probe_rank0) { 

      assert( ss);
      ofstream out_file;
      char filename[128];
      sprintf(filename,"%s.fp",name.c_str());
      if (interval == -1) {
        const string tmp_filename = MiscUtils::makeTmpPrefix(filename);
        out_file.open(tmp_filename.c_str(),ofstream::trunc);
        assert(out_file.is_open());
        out_file << ss->rdbuf();
        out_file.close();
        remove(filename);
        rename(tmp_filename.c_str(),filename);
      }
      else {
        out_file.open(filename,ofstream::app);
        assert( out_file.is_open());
        out_file << ss->rdbuf();
        out_file.close();
      }
      ss->str(string()); // clears ss
    }
  }

  void probe_reduce(double * pack_buf, double * unpack_buf, const int nd) { 

    // since there was only one group, it may be more efficient to reduce
    // the data via a straight reduction (log(p) scaling) as opposed to 
    // the approach taken by the multifluxprobe; probe_reduce function is 
    // left in case we want to change this later.. 

    MPI_Reduce(pack_buf,unpack_buf,nd,MPI_DOUBLE,MPI_SUM,probe_rank0,mpi_comm);

  }
};

#endif
