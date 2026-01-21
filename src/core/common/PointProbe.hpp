#ifndef _POINT_PROBE_HPP_
#define _POINT_PROBE_HPP_

class PointProbe {
public:

  CtiRegister::CtiDataProducer * solver;  // back pointer to who created you

  string name; //,header;
  int interval;

  // the varVec cna store a persistent pointer to registered CtiData* for efficiency.
  // a null second indicates an expression that needs to be evaluated...
  // vector<pair<string,CtiRegister::CtiData*> > varVec; // on all ranks
  vector<pair<string,vector<CtiRegister::CtiData*> > > varVec; // on all ranks
  int rank0; // rank responsible for writing
  int np;  // number of probe points
  int ns;  // number of data/scalar values per point

  vector<pair<string,std::stringstream*> > ssVec; // rank0 only

  bool b_surface;
  bool b_cht;

  int pack_buf_size;
  int * pack_buf_index;  // stores icv or ibf (for surface snapped)
  int * pack_buf_zone;  // stores izone for surface snapped points
  double * pack_buf_wgt;

  int active_rank_count;
  int * raoar; // rank-of-active-rank
  int * ioar_i; // index-of-active-rank i
  int * ioar_v; // index-of-active-rank v (i.e. what point probe to add contribution to)

  int unpack_buf_size;

  PointProbe() {

    assert(0); // march 2019

    solver = NULL;

    interval = -1;
    rank0 = -1;
    np = -1;
    ns = 0;

    pack_buf_size = 0;
    pack_buf_index = NULL;
    pack_buf_zone = NULL;
    pack_buf_wgt = NULL;

    raoar = NULL;
    ioar_i = NULL;
    ioar_v = NULL;

    b_surface = false;
    b_cht = false;

    unpack_buf_size = 0;

  }

  PointProbe(CtiRegister::CtiDataProducer * _solver) {

    solver = _solver;

    interval = -1;
    rank0 = -1;
    np = -1;
    ns = 0;

    pack_buf_size = 0;
    pack_buf_index = NULL;
    pack_buf_zone = NULL;
    pack_buf_wgt = NULL;

    raoar = NULL;
    ioar_i = NULL;
    ioar_v = NULL;

    b_surface = false;
    b_cht = false;

    unpack_buf_size = 0;

  }

  ~PointProbe() {

    for (int ii = 0; ii < ssVec.size(); ++ii) delete ssVec[ii].second; // the ss ptr

    DELETE(pack_buf_index);
    DELETE(pack_buf_zone);
    DELETE(pack_buf_wgt);

    DELETE(raoar);
    DELETE(ioar_i);
    DELETE(ioar_v);

  }

  void doProbe(const int step,const double time) {
    // should have some valid vars...

    const int nvar = varVec.size();
    assert(nvar > 0);
    assert(ns > 0);

    // get/evaluate the variables required for the probe. They
    // should all be valid because they are checked when the
    // probe is initialized in initProbe...

    double * pack_buf_double = new double[pack_buf_size*ns];

    int is = 0; // number of scalars, same as number of files...
    for (int ii = 0; ii < nvar; ++ii) {

      CtiRegister::CtiData * data = NULL;

      // check if variable requires zonename expansion
      if (varVec[ii].first.find("*:") != string::npos) {
        int current_izone = -1;
        int active_zone_index = 0;
        int d_type = 0;
        for (int ipack = 0; ipack < pack_buf_size; ++ipack) {
          // only update the CtiData* when a new zone is used; pack_bufs should have been sorted by zone in init
          if (pack_buf_zone[ipack] != current_izone) {
            current_izone = pack_buf_zone[ipack];
            if (varVec[ii].second[active_zone_index] == NULL) {
              // unregistered data so we'll need to create string accessor
              string zone_name;
              const bool found = solver->getSpecifierByIndex(zone_name,current_izone);  // init checks that is valid for all packed zones
              assert(found);
              string data_name = varVec[ii].first;
              data_name = MiscUtils::replaceAll(data_name,"*:",(zone_name+":"));
              data = CtiRegister::getUnregisteredCtiData(data_name);
            }
            else {
              data = varVec[ii].second[active_zone_index];
            }
            ++active_zone_index;
          }

          assert((data != NULL) && (data->getUnindexedTopology() == BF_DATA));

          // pack the buffer
          if (data->getType() == DN_DATA) {
            const int index = pack_buf_index[ipack];  // icv or ibf
            assert((index >= 0)&&(index < data->size()));
            pack_buf_double[pack_buf_size*is + ipack] = pack_buf_wgt[ipack] * data->dn(index);
            d_type |= 1;
          }
          else if (data->getType() == DN3_DATA) {
            FOR_I3 {
              const int index = pack_buf_index[ipack];
              assert((index >= 0)&&(index < data->size()));
              pack_buf_double[pack_buf_size*(is+i) + ipack] = pack_buf_wgt[ipack]*data->dn3(index,i);
              d_type |= 2;
            }
          }
          else {
            assert(0);
          }
        }
        int my_d_type = d_type;
        MPI_Allreduce(&my_d_type,&d_type, 1, MPI_INT, MPI_MAX, mpi_comm);
        assert((d_type==1)||(d_type==2));
        is += (d_type == 1) ? 1:3;  // update scalar count
      }
      else {
        // traditional variable without zonename expansion, i.e., single data pointer
        assert(int(varVec[ii].second.size()) == 1);

        if (varVec[ii].second[0] != NULL) data = varVec[ii].second[0]; // data is registered
        else data = CtiRegister::getUnregisteredCtiData(varVec[ii].first); // data needs to be evaluated

        //assert((data != NULL) && ((data->getUnindexedTopology()==BF_DATA)||(data->getUnindexedTopology()==CV_DATA)) && ((data->getType()==DN_DATA)||(data->getType()==DN3_DATA)));
        assert((data != NULL) && ((data->getType()==DN_DATA)||(data->getType()==DN3_DATA)));
  

        if (data->getType() == DN_DATA) {
          for (int ipack = 0; ipack < pack_buf_size; ++ipack) {
            const int index = pack_buf_index[ipack];  // icv or ibf or ino
            assert((index >= 0)&&(index < data->size()));
            pack_buf_double[pack_buf_size*is + ipack] = pack_buf_wgt[ipack]*data->dn(index);
          }
          is += 1;
        }
        else if (data->getType() == DN3_DATA) {
          FOR_I3 {
            for (int ipack = 0; ipack < pack_buf_size; ++ipack) {
              const int index = pack_buf_index[ipack];
              assert((index >= 0)&&(index < data->size()));
              pack_buf_double[pack_buf_size*(is+i) + ipack] = pack_buf_wgt[ipack]*data->dn3(index,i);
            }
          }
          is += 3;
        }
        else {
          assert(0);
        }
      }
    }
    assert(is==ns);

    //cout << " - the write rank is: " << rank0 << ", I am: " << mpi_rank << endl;

    if (mpi_rank == rank0) {

      double * probe_buf = new double[np*ns]; // number of points * number of scalars
      for (int i = 0; i < np*ns; ++i) probe_buf[i] = 0.0;

      // pull data from active ranks...

      double * unpack_buf_double = NULL;
      for (int iar = 0; iar < active_rank_count; ++iar) {
        const int rank = raoar[iar]; assert((rank >= 0)&&(rank < mpi_size));
        const int unpack_size = ioar_i[iar+1]-ioar_i[iar];
        if (rank == mpi_rank) {
          for (int ioa = ioar_i[iar]; ioa != ioar_i[iar+1]; ++ioa) {
            const int ip = ioar_v[ioa]; assert((ip >= 0)&&(ip < np));
            for (int is = 0; is < ns; ++is) probe_buf[is*np+ip] += pack_buf_double[is*unpack_size+ioa-ioar_i[iar]];
          }
        }
        else {
          if (unpack_buf_double == NULL) {
            // allocate unpack_buf_double so it is large enough for the biggest exchange...
            // TODO: should pre-compute unpack_size_max
            int unpack_size_max = 0;
            for (int iar2 = 0; iar2 < active_rank_count; ++iar2) {
              if (raoar[iar2] != mpi_rank) {
                unpack_size_max = max(unpack_size_max,ioar_i[iar2+1]-ioar_i[iar2]);
              }
            }
            assert(unpack_size_max > 0);
            unpack_buf_double = new double[unpack_size_max*ns];
          }
          MPI_Recv(unpack_buf_double,unpack_size*ns,MPI_DOUBLE,rank,52528,mpi_comm,MPI_STATUS_IGNORE);
          for (int ioa = ioar_i[iar]; ioa != ioar_i[iar+1]; ++ioa) {
            const int ip = ioar_v[ioa]; assert((ip >= 0)&&(ip < np));
            for (int is = 0; is < ns; ++is) probe_buf[is*np+ip] += unpack_buf_double[is*unpack_size+ioa-ioar_i[iar]];
          }
        }
      }
      DELETE(unpack_buf_double);

      // take a look at the data unrolled...
      /*
        for (int i = 0; i < np*ns; ++i) {
          cout << "i : " << i << " probe_buf[i]: " << probe_buf[i] << endl;
        }
      */

      // write data to stringstream...

      assert(ssVec.size() == ns);
      for (int is = 0; is < ns; ++is) {
        *ssVec[is].second << step << " " << time << " " << np;
        for (int ip = 0; ip < np; ++ip)
          *ssVec[is].second << " " << probe_buf[is*np+ip];
        *ssVec[is].second << endl;
      }

      delete[] probe_buf;
    }
    else {
      if (pack_buf_size > 0)
        MPI_Ssend(pack_buf_double,pack_buf_size*ns,MPI_DOUBLE,rank0,52528,mpi_comm);
    }

    DELETE(pack_buf_double);
  }

  void flush(const bool verbose = true) {

    if (mpi_rank == rank0) {

      ofstream outFile;
      for (int ii = 0; ii < ssVec.size(); ++ii) {
        ssVec[ii].second->rdbuf()->pubseekpos(0);  // set read position to start of stringbuf - only seemed to be required for certain compilers...
        if (ssVec[ii].second->rdbuf()->in_avail() == 0) {
          // if the stringstream is empty, report...
          if (verbose) cout << "Warning: PointProbe \"" << name << "\" requesting flush of an empty buffer. Reduce flush rate." << endl;
        }
        else {
          if (interval == -1) {
            const string filename = MiscUtils::cleanFilename(ssVec[ii].first,name.size());
            const string tmp_filename = MiscUtils::makeTmpPrefix(filename);
            outFile.open(tmp_filename.c_str(),ofstream::trunc);
            assert(outFile.is_open());
            outFile << ssVec[ii].second->rdbuf();
            outFile.close();
            remove(filename.c_str());
            rename(tmp_filename.c_str(),filename.c_str());
          }
          else {
            outFile.open((MiscUtils::cleanFilename(ssVec[ii].first,name.size())).c_str(),ofstream::app);
            assert(outFile.is_open());
            outFile << ssVec[ii].second->rdbuf();
            outFile.close();
          }
          ssVec[ii].second->str(string()); // this clears ss, believe it or not...
        }
      }
    }

  }

};

#endif
