#ifndef VOLUMETRICPROBE_HPP
#define VOLUMETRICPROBE_HPP

class VolumetricProbe {
public: 

  string name;
  int interval;
  vector<string> var_vec;
  int probe_rank0;
  stringstream * ss;
  vector<int> cv_vec;

  double probe_x[3];
  double probe_vol;

  VolumetricProbe() { 

    name        = "";
    interval    = -1;
    probe_rank0 = -1;
    ss          = NULL;

    probe_vol = 0.0;
    FOR_I3 probe_x[i] = 0.0;
  }

  ~VolumetricProbe() {
    if ( ss != NULL) delete ss;
  }

  void doProbe(const int step = -1, const double time = 0.0) { 

    const int nvar = var_vec.size();
    assert( nvar > 0);

    vector<CtiRegister::CtiData*> data_vec(nvar);
    for (int ii = 0; ii < nvar; ++ii) { 
      data_vec[ii] = CtiRegister::getCtiData(var_vec[ii]);
      assert( (data_vec[ii] != NULL) && (data_vec[ii]->getTopology() == CV_DATA) &&
              (data_vec[ii]->getType() == DN_DATA));
    }

    double * pack_buf = new double[nvar];
    for (int ivar = 0; ivar < nvar; ++ivar) 
      pack_buf[ivar] = 0.0;

    CtiRegister::CtiData* vol_data = CtiRegister::getRegisteredCtiData("vol_cv");
    const int nc = cv_vec.size();
    for (int ii = 0; ii < nc; ++ii) { 

      int icv = cv_vec[ii];

      for (int ivar = 0; ivar < nvar; ++ivar) 
        pack_buf[ivar] += vol_data->dn(icv)*data_vec[ivar]->dn(icv);

    }

    double * unpack_buf = NULL;
    if ( mpi_rank == probe_rank0) 
      unpack_buf = new double[nvar];

    probe_reduce(pack_buf,unpack_buf,nvar);

    delete[] pack_buf;

    if ( mpi_rank == probe_rank0) { 

      *ss << step << " " << time;
      for (int ivar = 0; ivar < nvar; ++ivar) 
        *ss << "  " << unpack_buf[ivar]/probe_vol;
      *ss << endl;

      delete[] unpack_buf;
    }
  }


  void doSMDProbe(const int step = -1, const double time = 0.0) { 

    const int nvar = var_vec.size();
    int vol_index = -1;
    int sur_index = -1;
    assert( nvar > 0);

    vector<CtiRegister::CtiData*> data_vec(nvar);
    for (int ii = 0; ii < nvar; ++ii) { 
      data_vec[ii] = CtiRegister::getCtiData(var_vec[ii]);
      assert( (data_vec[ii] != NULL) && (data_vec[ii]->getTopology() == CV_DATA) &&
              (data_vec[ii]->getType() == DN_DATA));
      if (var_vec[ii] == "vof" || var_vec[ii] == "vof_avg" ) vol_index = ii;
      if (var_vec[ii] == "lsd" || var_vec[ii] == "lsd_avg" ) sur_index = ii;
    }

    assert(vol_index >= 0);
    assert(sur_index >= 0);

    double * pack_buf = new double[nvar];
    for (int ivar = 0; ivar < nvar; ++ivar) 
      pack_buf[ivar] = 0.0;

    CtiRegister::CtiData* vol_data = CtiRegister::getRegisteredCtiData("vol_cv");
    const int nc = cv_vec.size();
    for (int ii = 0; ii < nc; ++ii) { 

      int icv = cv_vec[ii];

      for (int ivar = 0; ivar < nvar; ++ivar)  
        pack_buf[ivar] += vol_data->dn(icv)*data_vec[ivar]->dn(icv);

    }

    double * unpack_buf = NULL;
    if ( mpi_rank == probe_rank0) 
      unpack_buf = new double[nvar];

    probe_reduce(pack_buf,unpack_buf,nvar);

    delete[] pack_buf;

    assert(nvar == 2);

    if ( mpi_rank == probe_rank0) { 

      *ss << step << " " << time;
      for (int ivar = 0; ivar < nvar; ++ivar) 
        *ss << "  " << unpack_buf[ivar]/probe_vol;
      *ss << " " << 6.0*unpack_buf[vol_index]/unpack_buf[sur_index] << endl;
      // SMD = 6*vof/lsd ....

      delete[] unpack_buf;
    }
  }


  void flush() {

    if ( mpi_rank == probe_rank0) { 

      assert( ss);
      ofstream out_file;
      char filename[128];
      sprintf(filename,"%s.svp",name.c_str());
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
