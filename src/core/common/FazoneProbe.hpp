#ifndef FAZONEPROBE_HPP
#define FAZONEPROBE_HPP

class FazoneProbe {
public:

  string name;
  int interval;
  vector<string> var_vec;
  int probe_rank0;
  stringstream * ss;
  string zone;
  vector<int> bf_vec;

  double probe_x[3];
  double probe_area;

  FazoneProbe() {

    name        = "";
    zone        = "";
    interval    = -1;
    probe_rank0 = -1;
    ss          = NULL;

    probe_area = 0.0;
    FOR_I3 probe_x[i] = 0.0;
  }

  FazoneProbe(const string _zone) : zone(_zone) {

    name        = "";
    interval    = -1;
    probe_rank0 = -1;
    ss          = NULL;

    probe_area = 0.0;
    FOR_I3 probe_x[i] = 0.0;
  }

  ~FazoneProbe() {
    if ( ss != NULL) delete ss;
  }

  void setZone(const string _zone) {
    zone = _zone;
  }

  void doProbe(const int step = -1, const double time = 0.0) {

    const int nvar = var_vec.size();
    assert( nvar > 0);

    vector<CtiRegister::CtiData*> data_vec(nvar);
    for (int ii = 0; ii < nvar; ++ii) {
      data_vec[ii] = CtiRegister::getCtiData(var_vec[ii]);

      // this boundary may not exist on this rank, so allow NULL
      if (data_vec[ii] != NULL) {
        assert( (data_vec[ii]->getUnindexedTopology() == BF_DATA) &&
              (data_vec[ii]->getType() == DN_DATA));
      }
    }


    double * pack_buf = new double[nvar];
    for (int ivar = 0; ivar < nvar; ++ivar)
      pack_buf[ivar] = 0.0;

    char area_name[128];
    sprintf(area_name,"%s:area_bf",zone.c_str());

    CtiRegister::CtiData* area_data = CtiRegister::getRegisteredCtiData(area_name);
    for (int ii = 0, nf=bf_vec.size(); ii < nf; ++ii) {

      int ibf = bf_vec[ii];

      for (int ivar = 0; ivar < nvar; ++ivar)
        if (data_vec[ivar] != NULL) pack_buf[ivar] += area_data->dn(ibf)*data_vec[ivar]->dn(ibf);
    }

    double * unpack_buf = NULL;
    if ( mpi_rank == probe_rank0)
      unpack_buf = new double[nvar];

    probe_reduce(pack_buf,unpack_buf,nvar);

    delete[] pack_buf;

    if ( mpi_rank == probe_rank0) {

      *ss << step << " " << time;
      for (int ivar = 0; ivar < nvar; ++ivar)
        *ss << "  " << unpack_buf[ivar] << "  " << unpack_buf[ivar]/probe_area;
      *ss << endl;

      delete[] unpack_buf;
    }
  }


  void flush() {

    if ( mpi_rank == probe_rank0) {

      assert( ss);
      ofstream out_file;
      char filename[128];
      sprintf(filename,"%s.fzp",name.c_str());
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
