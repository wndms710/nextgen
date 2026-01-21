#include "StaticSolver.hpp"
#include "GeomUtils.hpp"

/*
bool StaticSolver::getPeriodicR(double R[9],const int bits) const {

  // set I...
  R[0] = 1.0; R[1] = 0.0; R[2] = 0.0;
  R[3] = 0.0; R[4] = 1.0; R[5] = 0.0;
  R[6] = 0.0; R[7] = 0.0; R[8] = 1.0;

  // loop on bit pairs...
  bool has_R = false;
  FOR_I3 {
    if (bits & (1<<(2*i))) {
      // even bit pair set...
      assert(!(bits & (1<<(2*i+1)))); // only one pair set
      assert(periodicTransformVec.size() > i); // make sure transform vec supports this bit pair
      double this_R[9];
      if (periodicTransformVec[i].getR(this_R)) {
        has_R = true;
        double Rtmp[9] = { 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0,
          0.0, 0.0, 0.0 };
        FOR_M3 {
          FOR_J3 {
            FOR_K3 {
              Rtmp[3*m+j] += R[3*m+k]*this_R[3*k+j];
            }
          }
        }
        for (int ii =0; ii < 9; ++ii)
          R[ii] = Rtmp[ii];
      }
    }
    else if (bits & (1<<(2*i+1))) {
      // odd bit pair set...
      assert(periodicTransformVec.size() > i); // make sure transform vec supports this bit pair
      double this_R[9];
      if (periodicTransformVec[i].getInvR(this_R)) {
        has_R = true;
        double Rtmp[9] = { 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0,
          0.0, 0.0, 0.0 };
        FOR_M3 {
          FOR_J3 {
            FOR_K3 {
              Rtmp[3*m+j] += R[3*m+k]*this_R[3*k+j];
            }
          }
        }
        for (int ii =0; ii < 9; ++ii)
          R[ii] = Rtmp[ii];
      }
    }
  }
  return has_R;
}

bool StaticSolver::getPeriodicT(double t[3],const int bits) const {

  // zero t...
  t[0] = 0.0; t[1] = 0.0; t[2] = 0.0;

  // loop on bit pairs...
  bool has_t = false;
  FOR_I3 {
    if (bits & (1<<(2*i))) {
      // even bit pair set...
      assert(!(bits & (1<<(2*i+1)))); // only one pair set
      assert(periodicTransformVec.size() > i); // make sure transform vec supports this bit pair
      double this_t[3];
      if (periodicTransformVec[i].getT(this_t)) {
        has_t = true;
        FOR_J3 t[j] += this_t[j];
      }
    }
    else if (bits & (1<<(2*i+1))) {
      // odd bit pair set...
      assert(int(periodicTransformVec.size()) > i); // make sure transform vec supports this bit pair
      double this_t[3];
      if (periodicTransformVec[i].getInvT(this_t)) {
        has_t = true;
        FOR_J3 t[j] += this_t[j];
      }
    }
  }
  return has_t;
}
*/

void StaticSolver::readData(const string& filename) {

  CtiRegister::clearAllDataFlags();
  CtiRegister::readData(filename);
  if (lpHelperVec.size() > 0) {

    // need to (re)build the dde stuff using the current snapshots lpocv_i_global...
    clearLpData();
    readLpocvAndInitDdeStuff(filename);
    CtiRegister::readLpData(filename);

    // set particle icv...
    for (int ii = 0, lim = lpHelperVec.size(); ii < lim; ++ii) {

      CtiRegister::CtiData * data = CtiRegister::getRegisteredCtiData(lpHelperVec[ii].name+":icv",false);
      assert((data)&&(data->getType() == IN_DATA)&&(data->getUnindexedTopology() == LP_DATA));
      const int nlp = data->size();
      for (int ilp = 0; ilp < nlp; ++ilp)
        data->in(ilp) = lpHelperVec[ii].cvolp[ilp];

      DELETE(lpHelperVec[ii].cvolp);
    }
  }

}

void StaticSolver::readLpocvAndInitDdeStuff(const string filename) {

  COUT1("readLpocvAndInitDdeStuff()");

  MPI_File fh;
  char dummy[128];
  sprintf(dummy,"%s",filename.c_str());
  int ierr = MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

  if (ierr != 0) {
    CERR("cannot open restart file: " << filename);
  }

  int itmp[2] = { 0, 0 };
  if (mpi_rank == 0) {
    // first 2 ints are: 0. magic number, 1. io version
    MPI_File_read(fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
  }
  MPI_Bcast(itmp,2,MPI_INT,0,mpi_comm);

  bool byte_swap = false;
  if (itmp[0] != UGP_IO_MAGIC_NUMBER+1) {
    ByteSwap::byteSwap(itmp,2);
    if (itmp[0] != UGP_IO_MAGIC_NUMBER+1) {
      CERR("result file does not start as expected.");
    }
    COUT1(" > file requires byte swapping.");
    byte_swap = true;
  }

  int io_version = itmp[1];
  if (io_version != 5) {
    CERR("result file version not 5: " << io_version);
  }

  MPI_Offset offset = 8; // 2 ints

  Header header;
  int done = 0;
  //int record_count = 0;
  while (done != 1) {

    if (mpi_rank == 0) {
      MPI_File_read_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      if (byte_swap) ByteSwap::byteSwapHeader(&header,1);
    }
    MPI_Bcast(&header,1,MPI_Header,0,mpi_comm);

    switch (header.id) {

      case UGP_IO_LPOCV_I:

        {
          if (mpi_rank == 0) {
            cout << " > LPOCV_I \"" << header.name << "\"" << endl;
            cout.flush();
          }

          const string name = header.name;
          size_t colon_pos = name.find(":");
          assert(colon_pos != string::npos); // names must be prepended with particle class name
          const string lp_name  = name.substr(0,colon_pos);

          // find/build lp helper vec entry...

          map<const string,int>::iterator iter2 = lpHelperNameMap.find(lp_name);
          assert(iter2 != lpHelperNameMap.end());
          const int lp_index = iter2->second;
          assert((lp_index >= 0)&&(lp_index < lpHelperVec.size()));
          assert(lpHelperVec[lp_index].name  == lp_name);
          assert(lpHelperVec[lp_index].index == lp_index);

          // clear arrays, (re)building them here...

          lpHelperVec[lp_index].clear();

          assert(ncv_global+1 == ByteSwap::getInt8FromLswMswPair(header.idata+0));
          lpHelperVec[lp_index].nlp_global = ByteSwap::getInt8FromLswMswPair(header.idata+2);

          // read in lpocv_i_global...

          const int ncv_striped = cvora_striped[mpi_rank+1]-cvora_striped[mpi_rank];
          int8* lpocv_i_global = new int8[ncv_striped+1];
          readChunkedData<int8>(fh,offset+header_size+cvora_striped[mpi_rank]*int8_size,lpocv_i_global,ncv_striped+1,byte_swap,mpi_comm);
          int nlp_striped = lpocv_i_global[ncv_striped]-lpocv_i_global[0];

          {
            // setup send side stuff...

            int * send_count = new int[mpi_size];
            FOR_RANK send_count[rank] = 0;
            for (int icv_striped = 0; icv_striped < ncv_striped; ++icv_striped) {
              int rank,bits,icv;
              BitUtils::unpackRankBitsIndex(rank,bits,icv,rbi_sm[icv_striped]);
              assert(bits == 0);
              // icv,count and ilps...
              send_count[rank] += 2 + lpocv_i_global[icv_striped+1]-lpocv_i_global[icv_striped];
            }

            int * send_disp = new int[mpi_size];
            send_disp[0] = 0;
            for (int rank = 1; rank < mpi_size; ++rank)
              send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
            int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

            // build lpora...

            MiscUtils::buildXora(lpHelperVec[lp_index].lpora_cv_striped,nlp_striped);

            int * send_buf_int = new int[send_count_sum];
            for (int icv_striped = 0; icv_striped < ncv_striped; ++icv_striped) {
              int rank,bits,icv;
              BitUtils::unpackRankBitsIndex(rank,bits,icv,rbi_sm[icv_striped]);
              assert(bits == 0);
              // pack icv (promoted to int8) and ip_globals associated to icv...
              send_buf_int[send_disp[rank]++] = icv;
              send_buf_int[send_disp[rank]++] = lpocv_i_global[icv_striped+1]-lpocv_i_global[icv_striped];
              for (int loc = lpocv_i_global[icv_striped]; loc < lpocv_i_global[icv_striped+1]; ++loc)
                send_buf_int[send_disp[rank]++] = loc-lpHelperVec[lp_index].lpora_cv_striped[mpi_rank];
            }
            delete[] lpocv_i_global;

            // rewind...

            send_disp[0] = 0;
            for (int rank = 1; rank < mpi_size; ++rank)
              send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

            // setup recv side stuff...

            int * recv_count = new int[mpi_size];
            MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

            int * recv_disp = new int[mpi_size];
            recv_disp[0] = 0;
            for (int rank = 1; rank < mpi_size; ++rank)
              recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
            int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

            // exchange...

            int* recv_buf_int = new int[recv_count_sum];
            MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
                recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
            delete[] send_buf_int; send_buf_int = NULL;
            delete[] send_count;
            delete[] send_disp;

            int nlp = 0;
            FOR_RANK {
              int irecv = recv_disp[rank];
              while (irecv < recv_disp[rank]+recv_count[rank]) {
                //const int icv = recv_buf_int[irecv++];
                ++irecv;
                const int nloc = recv_buf_int[irecv++];
                nlp += nloc;
                irecv += nloc;
              }
            }

            *lpHelperVec[lp_index].size_ptr = nlp; // update registered data size

            int8* ilp_global = new int8[nlp];
            lpHelperVec[lp_index].cvolp = new int[nlp];
            nlp = 0;
            FOR_RANK {
              int irecv = recv_disp[rank];
              while (irecv < recv_disp[rank]+recv_count[rank]) {
                const int icv = recv_buf_int[irecv++];
                const int nloc = recv_buf_int[irecv++];
                for (int loc = 0; loc < nloc; ++loc) {
                  ilp_global[nlp] = recv_buf_int[irecv++]+lpHelperVec[lp_index].lpora_cv_striped[rank];
                  lpHelperVec[lp_index].cvolp[nlp] = icv;
                  ++nlp;
                }
              }
            }
            assert(*lpHelperVec[lp_index].size_ptr == nlp);
            delete[] recv_buf_int; recv_buf_int = NULL;
            delete[] recv_count;
            delete[] recv_disp;

            // build dde...
            lpHelperVec[lp_index].dde_cv_striped = new DistributedDataExchanger(ilp_global,nlp,lpHelperVec[lp_index].lpora_cv_striped);
            delete[] ilp_global;

            resizeData();
          }
        }

        break;

      case UGP_IO_EOF:

        done = 1;
        break;

        /*
           default:
           if (mpi_rank == 0) cout << " > header: " << header.id << " \"" << header.name << "\", skipping." << endl;
           */

    }

    offset += header.skip;

  }

  MPI_File_close(&fh);

}

void StaticSolver::redistReorderData(StripedMesh* sm) {

  COUT1("StaticSolver::redistReorderData()");

  // I_DATA...

  for (list<pair<CtiRegister::CtiData*,int> >::iterator iter = sm->iList.begin(); iter != sm->iList.end(); ++iter) {
    assert(iter->first != NULL);
    assert(iter->first->getType() == I_DATA);
    iter->first->i() = iter->second; // set value here
    iter->first->setFlag(1);
  }

  // D_DATA...

  for (list<pair<CtiRegister::CtiData*,double> >::iterator iter = sm->dList.begin(); iter != sm->dList.end(); ++iter) {
    assert(iter->first != NULL);
    assert(iter->first->getType() == D_DATA);
    iter->first->d() = iter->second;
    iter->first->setFlag(1);
  }

  // IN_DATA...

  for (list<InData>::iterator iter = sm->inList.begin(); iter != sm->inList.end(); ++iter) {
    assert(iter->ctiData != NULL);
    iter->ctiData->redistReorderData(iter->data);
  }

  // DN_DATA...

  for (list<DnData>::iterator iter = sm->dnList.begin(); iter != sm->dnList.end(); ++iter) {
    assert(iter->ctiData != NULL);
    iter->ctiData->redistReorderData(iter->data);
  }

  // DN3_DATA...

  for (list<Dn3Data>::iterator iter = sm->dn3List.begin(); iter != sm->dn3List.end(); ++iter) {
    assert(iter->ctiData != NULL);
    iter->ctiData->redistReorderData(iter->data);
  }

  // set particle icv...
  for (int ii = 0, lim = lpHelperVec.size(); ii < lim; ++ii) {

    CtiRegister::CtiData * data = CtiRegister::getRegisteredCtiData(lpHelperVec[ii].name+":icv",false);
    assert((data)&&(data->getType() == IN_DATA)&&(data->getUnindexedTopology() == LP_DATA));
    const int nlp = data->size();
    for (int ilp = 0; ilp < nlp; ++ilp)
      data->in(ilp) = lpHelperVec[ii].cvolp[ilp];

    DELETE(lpHelperVec[ii].cvolp);
  }

  // orient face data to respect local face direction...

  flipRWSignedFaData(); // CI why this?

}

void StaticSolver::flipRWSignedFaData() {
  // use ifa global sign to determine if we should flip registered FA/DN3_DATA with rw bits

  for (map<const string,CtiRegister::CtiData>::iterator iter = CtiRegister::registeredDataMap.begin(); iter != CtiRegister::registeredDataMap.end(); ++iter) {
    if ( iter->second.getTopology() == SIGNED_FA_DATA &&
        (iter->second.checkBit(CAN_WRITE_DATA) || iter->second.checkBit(WRITE_DATA) || iter->second.checkBit(READ_DATA)) ) {
      assert(iter->second.size() == nfa);
      if (iter->second.getType() == DN_DATA) {
        FOR_IFA {
          if (ifa_global[ifa] < 0)
            iter->second.dn(ifa) = -iter->second.dn(ifa);
        }
      }
      else if (iter->second.getType() == DN3_DATA) {
        FOR_IFA {
          if (ifa_global[ifa] < 0)
            FOR_I3 iter->second.dn3(ifa,i) = -iter->second.dn3(ifa,i);
        }
      }
      else {
        assert(0);
      }
    }
  }

}

void StaticSolver::registerBoundaryConditions() {

  COUT1("StaticSolver::registerBoundaryConditions()");

  FOR_IZONE(bfZoneVec) {

    assert(bfZoneVec[izone].area_global > 0.0); // must be called after the zone is initialized
    CtiRegister::_registerData(bfZoneVec[izone].index,bfZoneVec[izone].getName()+":index",NO_READWRITE_DATA);
    CtiRegister::_registerData(bfZoneVec[izone].area_global,bfZoneVec[izone].getName()+":area_global",NO_READWRITE_DATA);
    CtiRegister::_registerData(bfZoneVec[izone].area_over_delta_global,bfZoneVec[izone].getName()+":area_over_delta_global",NO_READWRITE_DATA);
    CtiRegister::_registerData(bfZoneVec[izone].x_global,bfZoneVec[izone].getName()+":x_global",NO_READWRITE_DATA); // can fxn know diff b/w double[3] and *double?
    CtiRegister::_registerData(bfZoneVec[izone].n_global,bfZoneVec[izone].getName()+":n_global",NO_READWRITE_DATA);

    // these CAN be written to lets say a snapshot...
    bfZoneVec[izone].registerBfData(bfZoneVec[izone].area_bf,"area_bf",CAN_WRITE_DATA);
    bfZoneVec[izone].registerBfData(bfZoneVec[izone].area_over_delta_bf,"area_over_delta_bf",CAN_WRITE_DATA);
    bfZoneVec[izone].registerBfData(bfZoneVec[izone].x_bf,"x_bf",CAN_WRITE_DATA|X_DATA); // X_DATA sets this as the "x" for this "n" type
    bfZoneVec[izone].registerBfData(bfZoneVec[izone].n_bf,"n_bf",CAN_WRITE_DATA);

    // also register a global index and name for this "n" type...
    CtiRegister::_registerGlobalIndex(bfZoneVec[izone].ibf_global,bfZoneVec[izone].nbf);
    CtiRegister::_registerName(bfZoneVec[izone].getName(),bfZoneVec[izone].nbf);

  }

}

void StaticSolver::buildCvAdt() {

  COUT1("StaticSolver::buildCvAdt()");

  assert(x_vv);
  assert(r_vv);
  assert(cvBboxAdt == NULL);
  assert(cvAdt == NULL);

  double (*bbmin)[3] = new double[ncv][3];
  double (*bbmax)[3] = new double[ncv][3];

  FOR_ICV {
    FOR_I3 bbmin[icv][i] = 1.0E+20;
    FOR_I3 bbmax[icv][i] = -1.0E+20;
  }

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
    const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
    for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      const int ino = noofa_v[nof]; assert((ino >= 0)&&(ino < nno));
      FOR_I3 bbmin[icv0][i] = min(bbmin[icv0][i],x_no[ino][i]);
      FOR_I3 bbmax[icv0][i] = max(bbmax[icv0][i],x_no[ino][i]);
      FOR_I3 bbmin[icv1][i] = min(bbmin[icv1][i],x_no[ino][i]);
      FOR_I3 bbmax[icv1][i] = max(bbmax[icv1][i],x_no[ino][i]);
    }
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
    for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      const int ino = noofa_v[nof]; assert((ino >= 0)&&(ino < nno));
      FOR_I3 bbmin[icv0][i] = min(bbmin[icv0][i],x_no[ino][i]);
      FOR_I3 bbmax[icv0][i] = max(bbmax[icv0][i],x_no[ino][i]);
    }
  }

  FOR_IBF {
    const int icv = cvobf[ibf]; assert((icv >= 0)&&(icv < ncv));
    for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
      const int ino = noobf_v[nob]; assert((ino >= 0)&&(ino < nno));
      FOR_I3 bbmin[icv][i] = min(bbmin[icv][i],x_no[ino][i]);
      FOR_I3 bbmax[icv][i] = max(bbmax[icv][i],x_no[ino][i]);
    }
  }

  // expand all bboxes by some small fraction of the local r_vv...

  double bbmin_local[3] = { 1.0E+20, 1.0E+20, 1.0E+20 };
  double bbmax_local[3] = { -1.0E+20, -1.0E+20, -1.0E+20 };
  FOR_ICV {
    FOR_I3 {
      bbmin[icv][i] -= 1.0E-6*r_vv[icv];
      bbmin_local[i] = min(bbmin_local[i],bbmin[icv][i]);
    }
    FOR_I3 {
      bbmax[icv][i] += 1.0E-6*r_vv[icv];
      bbmax_local[i] = max(bbmax_local[i],bbmax[icv][i]);
    }
  }

  cvAdt = new Adt<double>(ncv,bbmin,bbmax);

  delete[] bbmin;
  delete[] bbmax;

  double (*bbmin_global)[3] = new double[mpi_size][3];
  MPI_Allgather((double*)bbmin_local,3,MPI_DOUBLE,(double*)bbmin_global,3,MPI_DOUBLE,mpi_comm);
  double (*bbmax_global)[3] = new double[mpi_size][3];
  MPI_Allgather((double*)bbmax_local,3,MPI_DOUBLE,(double*)bbmax_global,3,MPI_DOUBLE,mpi_comm);

  cvBboxAdt = new Adt<double>(mpi_size,bbmin_global,bbmax_global);

  delete[] bbmin_global;
  delete[] bbmax_global;

}

void StaticSolver::averageBfToNo(double * v_no,const CtiRegister::CtiData * data) {

  assert(data->getType() == DN_DATA || data->getType() == IN_DATA);
  assert(data->getUnindexedTopology() == BF_DATA);

  double * wgt = new double[nno];

  FOR_INO {
    v_no[ino] = 0.0;
    wgt[ino] = 0.0;
  }

  const int iz = data->getIndex();
  if (data->getType() == DN_DATA) {
    for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
      assert(zone_bf[ibf] == iz);
      const double v = data->dn(ibf-bfZoneVec[iz].ibf_f);
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino = noobf_v[nob];
        assert(ino >= 0 && ino < nno);
        v_no[ino] += v;
        wgt[ino] += 1.0;
      }
    }
  }
  else {
    assert(data->getType() == IN_DATA);
    for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
      assert(zone_bf[ibf] == iz);
      const double v = data->in(ibf-bfZoneVec[iz].ibf_f);
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino = noobf_v[nob];
        assert(ino >= 0 && ino < nno);
        v_no[ino] += v;
        wgt[ino] += 1.0;
      }
    }
  }

  // need to do node update...
  updateNoData(v_no);
  updateNoData(wgt);

  FOR_INO {
    if (wgt[ino] > 0.0) {
      v_no[ino] /= wgt[ino];
    }
  }
  delete[] wgt;

}

void StaticSolver::averageBfToNo(double * v_no,const vector<CtiRegister::CtiData*> data_vec) {

  double * wgt = new double[nno];

  FOR_INO {
    v_no[ino] = 0.0;
    wgt[ino] = 0.0;
  }

  for (int id = 0, lim = data_vec.size(); id < lim; ++id) {
    CtiRegister::CtiData* data = data_vec[id];
    assert(data->getType() == DN_DATA || data->getType() == IN_DATA);
    assert(data->getUnindexedTopology() == BF_DATA);
    const int iz = data->getIndex();
    if (data->getType() == DN_DATA) {
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        const double v = data->dn(ibf-bfZoneVec[iz].ibf_f);
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int ino = noobf_v[nob];
          assert(ino >= 0 && ino < nno);
          v_no[ino] += v;
          wgt[ino] += 1.0;
        }
      }
    }
    else {
      assert(data->getType() == IN_DATA);
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        const double v = data->in(ibf-bfZoneVec[iz].ibf_f);
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int ino = noobf_v[nob];
          assert(ino >= 0 && ino < nno);
          v_no[ino] += v;
          wgt[ino] += 1.0;
        }
      }
    }
  }

  // need to do node update...
  updateNoData(v_no);
  updateNoData(wgt);

  FOR_INO {
    if (wgt[ino] > 0.0) {
      v_no[ino] /= wgt[ino];
    }
  }
  delete[] wgt;

}

void StaticSolver::averageBfToNo(double (*v_no)[3],const CtiRegister::CtiData * data) {

  assert(data->getType() == DN3_DATA);
  assert(data->getUnindexedTopology() == BF_DATA);

  double * wgt = new double[nno];

  FOR_INO {
    FOR_I3 v_no[ino][i] = 0.0;
    wgt[ino] = 0.0;
  }

  const int iz = data->getIndex();
  for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
    assert(zone_bf[ibf] == iz);
    double v[3]; FOR_I3 v[i] = data->dn3(ibf-bfZoneVec[iz].ibf_f,i);
    for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
      const int ino = noobf_v[nob];
      assert(ino >= 0 && ino < nno);
      FOR_I3 v_no[ino][i] += v[i];
      wgt[ino] += 1.0;
    }
  }

  // need to do node update...
  updateNoData(v_no);
  updateNoData(wgt);

  FOR_INO {
    if (wgt[ino] > 0.0) {
      FOR_I3 v_no[ino][i] /= wgt[ino];
    }
  }
  delete[] wgt;

}

void StaticSolver::averageBfToNo(double (*v_no)[3],const vector<CtiRegister::CtiData*> data_vec) {

  double * wgt = new double[nno];

  FOR_INO {
    FOR_I3 v_no[ino][i] = 0.0;
    wgt[ino] = 0.0;
  }

  for (int id = 0, lim = data_vec.size(); id < lim; ++id) {
    CtiRegister::CtiData* data = data_vec[id];
    assert(data->getType() == DN3_DATA);
    assert(data->getUnindexedTopology() == BF_DATA);
    const int iz = data->getIndex();
    for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
      assert(zone_bf[ibf] == iz);
      double v[3]; FOR_I3 v[i] = data->dn3(ibf-bfZoneVec[iz].ibf_f,i);
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino = noobf_v[nob];
        assert(ino >= 0 && ino < nno);
        FOR_I3 v_no[ino][i] += v[i];
        wgt[ino] += 1.0;
      }
    }
  }

  // need to do node update...
  updateNoData(v_no);
  updateNoData(wgt);

  FOR_INO {
    if (wgt[ino] > 0.0) {
      FOR_I3 v_no[ino][i] /= wgt[ino];
    }
  }
  delete[] wgt;

}

void StaticSolver::averageCvToNo(double * v_no,const CtiRegister::CtiData * iso) {

  assert(iso->getType() == DN_DATA || iso->getType() == IN_DATA);
  assert(iso->getUnindexedTopology() == CV_DATA);

  double * wgt = new double[nno];

  FOR_INO {
    v_no[ino] = 0.0;
    wgt[ino] = 0.0;
  }

  if (iso->getType() == DN_DATA) {

    FOR_IBF {
      const int icv = cvobf[ibf];
      const double v = iso->dn(icv);
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino = noobf_v[nob];
        assert(ino >= 0 && ino < nno);
        v_no[ino] += v;
        wgt[ino] += 1.0;
      }
    }

    FOR_INTERNAL_IFA {
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      const double vsum = iso->dn(icv0) + iso->dn(icv1);
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino = noofa_v[nof];
        assert(ino >= 0 && ino < nno);
        v_no[ino] += vsum;
        wgt[ino] += 2.0;
      }
    }

    FOR_INTERPROC_IFA {
      const int icv0 = cvofa[ifa][0];
      const double v = iso->dn(icv0);
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino = noofa_v[nof];
        assert(ino >= 0 && ino < nno);
        v_no[ino] += v;
        wgt[ino] += 1.0;
      }
    }

  }
  else {
    assert(iso->getType() == IN_DATA);

    FOR_IBF {
      const int icv = cvobf[ibf];
      const double v = iso->in(icv);
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino = noobf_v[nob];
        assert(ino >= 0 && ino < nno);
        v_no[ino] += v;
        wgt[ino] += 1.0;
      }
    }

    FOR_INTERNAL_IFA {
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      const double vsum = iso->in(icv0) + iso->in(icv1);
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino = noofa_v[nof];
        assert(ino >= 0 && ino < nno);
        v_no[ino] += vsum;
        wgt[ino] += 2.0;
      }
    }

    FOR_INTERPROC_IFA {
      const int icv0 = cvofa[ifa][0];
      const double v = iso->in(icv0);
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino = noofa_v[nof];
        assert(ino >= 0 && ino < nno);
        v_no[ino] += v;
        wgt[ino] += 1.0;
      }
    }
  }

  // need to do node update...

  updateNoData(v_no);
  updateNoData(wgt);

  FOR_INO {
    assert(wgt[ino] != 0.0);
    v_no[ino] /= wgt[ino];
  }

  delete[] wgt;

}

void StaticSolver::averageCvToNo(double (*v_no)[3],const CtiRegister::CtiData * iso) {

  assert(iso->getType() == DN3_DATA);
  assert(iso->getUnindexedTopology() == CV_DATA);

  double * wgt = new double[nno];

  FOR_INO {
    FOR_I3 v_no[ino][i] = 0.0;
    wgt[ino] = 0.0;
  }

  FOR_IBF {
    const int icv = cvobf[ibf];
    double v[3]; FOR_I3 v[i] = iso->dn3(icv,i);
    for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
      const int ino = noobf_v[nob];
      FOR_I3 v_no[ino][i] += v[i];
      wgt[ino] += 1.0;
    }
  }

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    double vsum[3]; FOR_I3 vsum[i] = iso->dn3(icv0,i) + iso->dn3(icv1,i);
    for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      const int ino = noofa_v[nof];
      FOR_I3 v_no[ino][i] += vsum[i];
      wgt[ino] += 2.0;
    }
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    double v[3]; FOR_I3 v[i] = iso->dn3(icv0,i);
    for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      const int ino = noofa_v[nof];
      FOR_I3 v_no[ino][i] += v[i];
      wgt[ino] += 1.0;
    }
  }

  // need to do node update...

  updateNoData(v_no);
  updateNoData(wgt);

  FOR_INO {
    assert(wgt[ino] != 0.0);
    FOR_I3 v_no[ino][i] /= wgt[ino];
  }

  delete[] wgt;

}

void StaticSolver::averageCvToNo(double * v_no,const double * v_cv) {

  double * wgt = new double[nno];

  FOR_INO {
    v_no[ino] = 0.0;
    wgt[ino] = 0.0;
  }

  FOR_IBF {
    const int icv = cvobf[ibf];
    for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
      const int ino = noobf_v[nob];
      assert(ino >= 0 && ino < nno);
      v_no[ino] += v_cv[icv];
      wgt[ino] += 1.0;
    }
  }

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    const double vsum = v_cv[icv0]+v_cv[icv1];
    for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      const int ino = noofa_v[nof];
      assert(ino >= 0 && ino < nno);
      v_no[ino] += vsum;
      wgt[ino] += 2.0;
    }
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      const int ino = noofa_v[nof];
      assert(ino >= 0 && ino < nno);
      v_no[ino] += v_cv[icv0];
      wgt[ino] += 1.0;
    }
  }

  // need to do node update...

  updateNoData(v_no);
  updateNoData(wgt);

  FOR_INO {
    assert(wgt[ino] != 0.0);
    v_no[ino] /= wgt[ino];
  }

  delete[] wgt;

}

void StaticSolver::averageCvToNo(double (*v_no)[3],const double (*v_cv)[3]) {

  double * wgt = new double[nno];

  FOR_INO {
    FOR_I3 v_no[ino][i] = 0.0;
    wgt[ino] = 0.0;
  }

  FOR_IBF {
    const int icv = cvobf[ibf];
    for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
      const int ino = noobf_v[nob];
      FOR_I3 v_no[ino][i] += v_cv[icv][i];
      wgt[ino] += 1.0;
    }
  }

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    double vsum[3]; FOR_I3 vsum[i] = v_cv[icv0][i]+v_cv[icv1][i];
    for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      const int ino = noofa_v[nof];
      FOR_I3 v_no[ino][i] += vsum[i];
      wgt[ino] += 2.0;
    }
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      const int ino = noofa_v[nof];
      FOR_I3 v_no[ino][i] += v_cv[icv0][i];
      wgt[ino] += 1.0;
    }
  }

  // need to do node update...

  updateNoData(v_no);
  updateNoData(wgt);

  FOR_INO {
    assert(wgt[ino] != 0.0);
    FOR_I3 v_no[ino][i] /= wgt[ino];
  }

  delete[] wgt;

}

void StaticSolver::buildIso(vector<pair<SimpleTri,int> >& triVec,const string& iso_var_name,const double iso_var_value,int * cv_hide_flag) {

  // convert the iso_var_name into a node-base DN...

  double * iso_var_no = new double[nno];
  setNoDN(iso_var_no,iso_var_name);

  buildIso(triVec,iso_var_no,iso_var_value,cv_hide_flag);

  delete[] iso_var_no;

}

void StaticSolver::buildIso(vector<pair<SimpleTri,int> >& triVec,const string& iso_var_name,const CtiRegister::CtiData* iso_var,const double iso_var_value,int * cv_hide_flag) {

  // convert the iso_var into a node-base DN...

  double * iso_var_no = new double[nno];
  setNoDN(iso_var_no,iso_var_name,iso_var);

  buildIso(triVec,iso_var_no,iso_var_value,cv_hide_flag);

  delete[] iso_var_no;

}

void StaticSolver::buildIso(vector<pair<SimpleTri,int> >& triVec,const double * iso_var_no,const double iso_var_value,int * cv_hide_flag) {

  MiscUtils::dumpRange(iso_var_no,nno,"iso_var_no");

  // build buffer to hold cv, bf, and fa data that might be necessary to build the iso...

  int * cv_flag = new int[ncv];
  setCvFlagFromIsoVar(cv_flag,iso_var_no,iso_var_value);  // could pass cv_hide flag here instead...?

  int * bf_flag = new int[nbf];
  int * fa_flag = new int[nfa];
  int ndata = setFlagsToIndexTmpData(cv_flag,bf_flag,fa_flag);

  double * iso_var_tmp = new double[ndata-nno];
  double (*x_tmp)[3] = new double[ndata-nno][3];
  double * wgt_tmp = new double[ndata-nno];

  for (int idata = 0; idata < ndata-nno; ++idata) {
    iso_var_tmp[idata] = 0.0;
    FOR_I3 x_tmp[idata][i] = 0.0;
    wgt_tmp[idata] = 0.0;
  }

  const bool cvs_hidden = !(cv_hide_flag == NULL);
  FOR_IBF {
    const int icv = cvobf[ibf];
    const bool cv_ignore = (cvs_hidden && cv_hide_flag[icv]);
    if ((cv_flag[icv] >= 0)&&(!cv_ignore)) {
      const int idata_cv = cv_flag[icv]-nno;
      assert(bf_flag[ibf] >= 0);
      const int idata_bf = bf_flag[ibf]-nno;
      double wgt_bf = 0.0;
      int ino1 = noobf_v[noobf_i[ibf+1]-1];
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino0 = ino1;
        ino1 = noobf_v[nob];
        const double wgt_ed = DIST(x_no[ino0],x_no[ino1]);
        iso_var_tmp[idata_bf] += wgt_ed*(iso_var_no[ino0] + iso_var_no[ino1]);
        FOR_I3 x_tmp[idata_bf][i] += wgt_ed*(x_no[ino0][i] + x_no[ino1][i]);
        wgt_bf += wgt_ed;
      }
      iso_var_tmp[idata_bf] /= 2.0*wgt_bf;
      FOR_I3 x_tmp[idata_bf][i] /= 2.0*wgt_bf;
      // also add to the cv...
      wgt_bf = MAG(n_bf[ibf]);
      iso_var_tmp[idata_cv] += wgt_bf*iso_var_tmp[idata_bf];
      FOR_I3 x_tmp[idata_cv][i] += wgt_bf*x_tmp[idata_bf][i];
      wgt_tmp[idata_cv] += wgt_bf;
    }
  }

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    const bool cv_ignore = (cvs_hidden && (cv_hide_flag[icv0] || cv_hide_flag[icv1]));
    if ( (!cv_ignore)&&((cv_flag[icv0] >= 0)||(cv_flag[icv1] >= 0)) ) {
      assert(fa_flag[ifa] >= 0);
      const int idata_fa = fa_flag[ifa]-nno;
      double wgt_fa = 0.0;
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        const double wgt_ed = DIST(x_no[ino0],x_no[ino1]);
        iso_var_tmp[idata_fa] += wgt_ed*(iso_var_no[ino0] + iso_var_no[ino1]);
        FOR_I3 x_tmp[idata_fa][i] += wgt_ed*(x_no[ino0][i] + x_no[ino1][i]);
        wgt_fa += wgt_ed;
      }
      iso_var_tmp[idata_fa] /= 2.0*wgt_fa;
      FOR_I3 x_tmp[idata_fa][i] /= 2.0*wgt_fa;
      // also add to the cv...
      wgt_fa = MAG(n_fa[ifa]);
      if (cv_flag[icv0] >= 0) {
        const int idata_cv0 = cv_flag[icv0]-nno;
        iso_var_tmp[idata_cv0] += wgt_fa*iso_var_tmp[idata_fa];
        FOR_I3 x_tmp[idata_cv0][i] += wgt_fa*x_tmp[idata_fa][i];
        wgt_tmp[idata_cv0] += wgt_fa;
      }
      if (cv_flag[icv1] >= 0) {
        const int idata_cv1 = cv_flag[icv1]-nno;
        iso_var_tmp[idata_cv1] += wgt_fa*iso_var_tmp[idata_fa];
        FOR_I3 x_tmp[idata_cv1][i] += wgt_fa*x_tmp[idata_fa][i];
        wgt_tmp[idata_cv1] += wgt_fa;
      }
    }
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    const bool cv_ignore = (cvs_hidden && cv_hide_flag[icv0]);
    if ((cv_flag[icv0] >= 0) && (!cv_ignore)) {
      const int idata_cv0 = cv_flag[icv0]-nno;
      assert(fa_flag[ifa] >= 0);
      const int idata_fa = fa_flag[ifa]-nno;
      double wgt_fa = 0.0;
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        const double wgt_ed = DIST(x_no[ino0],x_no[ino1]);
        iso_var_tmp[idata_fa] += wgt_ed*(iso_var_no[ino0] + iso_var_no[ino1]);
        FOR_I3 x_tmp[idata_fa][i] += wgt_ed*(x_no[ino0][i] + x_no[ino1][i]);
        wgt_fa += wgt_ed;
      }
      iso_var_tmp[idata_fa] /= 2.0*wgt_fa;
      FOR_I3 x_tmp[idata_fa][i] /= 2.0*wgt_fa;
      // also add to the cv...
      wgt_fa = MAG(n_fa[ifa]);
      iso_var_tmp[idata_cv0] += wgt_fa*iso_var_tmp[idata_fa];
      FOR_I3 x_tmp[idata_cv0][i] += wgt_fa*x_tmp[idata_fa][i];
      wgt_tmp[idata_cv0] += wgt_fa;
    }
  }

  // we now have variables at all faces, cvs and nodes, so now record the tets...

  FOR_ICV {
    const bool cv_ignore = (cvs_hidden && cv_hide_flag[icv]);
    if ((cv_flag[icv] >= 0)&&(!cv_ignore)) {
      const int idata_cv = cv_flag[icv]-nno;
      assert(wgt_tmp[idata_cv] > 0.0);
      iso_var_tmp[idata_cv] /= wgt_tmp[idata_cv];
      FOR_I3 x_tmp[idata_cv][i] /= wgt_tmp[idata_cv];
    }
  }

  delete[] wgt_tmp;

  FOR_IBF {
    const int icv = cvobf[ibf];
    const bool cv_ignore = (cvs_hidden && cv_hide_flag[icv]);
    if ((cv_flag[icv] >= 0)&&(!cv_ignore)) {
      const int idata_cv = cv_flag[icv]-nno;
      assert(bf_flag[ibf] >= 0);
      const int idata_bf = bf_flag[ibf]-nno;
      int ino1 = noobf_v[noobf_i[ibf+1]-1];
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino0 = ino1;
        ino1 = noobf_v[nob];
        tetCutPlane(triVec,
            iso_var_tmp[idata_cv]-iso_var_value,iso_var_tmp[idata_bf]-iso_var_value,iso_var_no[ino0]-iso_var_value,iso_var_no[ino1]-iso_var_value,
            x_tmp[idata_cv],x_tmp[idata_bf],x_no[ino0],x_no[ino1]);
      }
    }
  }

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    const bool cv_ignore = (cvs_hidden && (cv_hide_flag[icv0] || cv_hide_flag[icv1]));
    if ( (!cv_ignore)&&(((cv_flag[icv0] >= 0)||(cv_flag[icv1] >= 0))) ) {
      assert(fa_flag[ifa] >= 0);
      const int idata_fa = fa_flag[ifa]-nno;
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        if (cv_flag[icv0] >= 0) {
          const int idata_cv0 = cv_flag[icv0]-nno;
          tetCutPlane(triVec,
              iso_var_tmp[idata_cv0]-iso_var_value,iso_var_tmp[idata_fa]-iso_var_value,iso_var_no[ino0]-iso_var_value,iso_var_no[ino1]-iso_var_value,
              x_tmp[idata_cv0],x_tmp[idata_fa],x_no[ino0],x_no[ino1]);
        }
        if (cv_flag[icv1] >= 0) {
          const int idata_cv1 = cv_flag[icv1]-nno;
          tetCutPlane(triVec,
              iso_var_tmp[idata_cv1]-iso_var_value,iso_var_tmp[idata_fa]-iso_var_value,iso_var_no[ino1]-iso_var_value,iso_var_no[ino0]-iso_var_value,
              x_tmp[idata_cv1],x_tmp[idata_fa],x_no[ino1],x_no[ino0]);
        }
      }
    }
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    const bool cv_ignore = (cvs_hidden && cv_hide_flag[icv0]);
    if ((!cv_ignore)&&(cv_flag[icv0] >= 0)) {
      const int idata_cv0 = cv_flag[icv0]-nno;
      assert(fa_flag[ifa] >= 0);
      const int idata_fa = fa_flag[ifa]-nno;
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        tetCutPlane(triVec,
            iso_var_tmp[idata_cv0]-iso_var_value,iso_var_tmp[idata_fa]-iso_var_value,iso_var_no[ino0]-iso_var_value,iso_var_no[ino1]-iso_var_value,
            x_tmp[idata_cv0],x_tmp[idata_fa],x_no[ino0],x_no[ino1]);
      }
    }
  }

  delete[] x_tmp;
  delete[] iso_var_tmp;
  delete[] fa_flag;
  delete[] bf_flag;
  delete[] cv_flag;

}

void StaticSolver::buildIsoWithIcv(vector<pair<SimpleTri,int> >& triVec,const double * iso_var_no,const double iso_var_value) {

  //MiscUtils::dumpRange(iso_var_no,nno,"iso_var_no");

  // build buffer to hold cv, bf, and fa data that might be necessary to build the iso...

  int * cv_flag = new int[ncv];
  setCvFlagFromIsoVar(cv_flag,iso_var_no,iso_var_value);

  int * bf_flag = new int[nbf];
  int * fa_flag = new int[nfa];
  int ndata = setFlagsToIndexTmpData(cv_flag,bf_flag,fa_flag);

  double * iso_var_tmp = new double[ndata-nno];
  double (*x_tmp)[3] = new double[ndata-nno][3];
  double * wgt_tmp = new double[ndata-nno];

  for (int idata = 0; idata < ndata-nno; ++idata) {
    iso_var_tmp[idata] = 0.0;
    FOR_I3 x_tmp[idata][i] = 0.0;
    wgt_tmp[idata] = 0.0;
  }

  FOR_IBF {
    const int icv = cvobf[ibf];
    if (cv_flag[icv] >= 0) {
      const int idata_cv = cv_flag[icv]-nno;
      assert(bf_flag[ibf] >= 0);
      const int idata_bf = bf_flag[ibf]-nno;
      double wgt_bf = 0.0;
      int ino1 = noobf_v[noobf_i[ibf+1]-1];
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino0 = ino1;
        ino1 = noobf_v[nob];
        const double wgt_ed = DIST(x_no[ino0],x_no[ino1]);
        iso_var_tmp[idata_bf] += wgt_ed*(iso_var_no[ino0] + iso_var_no[ino1]);
        FOR_I3 x_tmp[idata_bf][i] += wgt_ed*(x_no[ino0][i] + x_no[ino1][i]);
        wgt_bf += wgt_ed;
      }
      iso_var_tmp[idata_bf] /= 2.0*wgt_bf;
      FOR_I3 x_tmp[idata_bf][i] /= 2.0*wgt_bf;
      // also add to the cv...
      wgt_bf = MAG(n_bf[ibf]);
      iso_var_tmp[idata_cv] += wgt_bf*iso_var_tmp[idata_bf];
      FOR_I3 x_tmp[idata_cv][i] += wgt_bf*x_tmp[idata_bf][i];
      wgt_tmp[idata_cv] += wgt_bf;
    }
  }

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    if ((cv_flag[icv0] >= 0)||(cv_flag[icv1] >= 0)) {
      assert(fa_flag[ifa] >= 0);
      const int idata_fa = fa_flag[ifa]-nno;
      double wgt_fa = 0.0;
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        const double wgt_ed = DIST(x_no[ino0],x_no[ino1]);
        iso_var_tmp[idata_fa] += wgt_ed*(iso_var_no[ino0] + iso_var_no[ino1]);
        FOR_I3 x_tmp[idata_fa][i] += wgt_ed*(x_no[ino0][i] + x_no[ino1][i]);
        wgt_fa += wgt_ed;
      }
      iso_var_tmp[idata_fa] /= 2.0*wgt_fa;
      FOR_I3 x_tmp[idata_fa][i] /= 2.0*wgt_fa;
      // also add to the cv...
      wgt_fa = MAG(n_fa[ifa]);
      if (cv_flag[icv0] >= 0) {
        const int idata_cv0 = cv_flag[icv0]-nno;
        iso_var_tmp[idata_cv0] += wgt_fa*iso_var_tmp[idata_fa];
        FOR_I3 x_tmp[idata_cv0][i] += wgt_fa*x_tmp[idata_fa][i];
        wgt_tmp[idata_cv0] += wgt_fa;
      }
      if (cv_flag[icv1] >= 0) {
        const int idata_cv1 = cv_flag[icv1]-nno;
        iso_var_tmp[idata_cv1] += wgt_fa*iso_var_tmp[idata_fa];
        FOR_I3 x_tmp[idata_cv1][i] += wgt_fa*x_tmp[idata_fa][i];
        wgt_tmp[idata_cv1] += wgt_fa;
      }
    }
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    if (cv_flag[icv0] >= 0) {
      const int idata_cv0 = cv_flag[icv0]-nno;
      assert(fa_flag[ifa] >= 0);
      const int idata_fa = fa_flag[ifa]-nno;
      double wgt_fa = 0.0;
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        const double wgt_ed = DIST(x_no[ino0],x_no[ino1]);
        iso_var_tmp[idata_fa] += wgt_ed*(iso_var_no[ino0] + iso_var_no[ino1]);
        FOR_I3 x_tmp[idata_fa][i] += wgt_ed*(x_no[ino0][i] + x_no[ino1][i]);
        wgt_fa += wgt_ed;
      }
      iso_var_tmp[idata_fa] /= 2.0*wgt_fa;
      FOR_I3 x_tmp[idata_fa][i] /= 2.0*wgt_fa;
      // also add to the cv...
      wgt_fa = MAG(n_fa[ifa]);
      iso_var_tmp[idata_cv0] += wgt_fa*iso_var_tmp[idata_fa];
      FOR_I3 x_tmp[idata_cv0][i] += wgt_fa*x_tmp[idata_fa][i];
      wgt_tmp[idata_cv0] += wgt_fa;
    }
  }

  // we now have variables at all faces, cvs and nodes, so now record the tets...

  FOR_ICV {
    if (cv_flag[icv] >= 0) {
      const int idata_cv = cv_flag[icv]-nno;
      assert(wgt_tmp[idata_cv] > 0.0);
      iso_var_tmp[idata_cv] /= wgt_tmp[idata_cv];
      FOR_I3 x_tmp[idata_cv][i] /= wgt_tmp[idata_cv];
    }
  }

  delete[] wgt_tmp;

  FOR_IBF {
    const int icv = cvobf[ibf];
    if (cv_flag[icv] >= 0) {
      const int idata_cv = cv_flag[icv]-nno;
      assert(bf_flag[ibf] >= 0);
      const int idata_bf = bf_flag[ibf]-nno;
      int ino1 = noobf_v[noobf_i[ibf+1]-1];
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino0 = ino1;
        ino1 = noobf_v[nob];
        const int triVec_size0 = triVec.size();
        tetCutPlane(triVec,
            iso_var_tmp[idata_cv]-iso_var_value,iso_var_tmp[idata_bf]-iso_var_value,iso_var_no[ino0]-iso_var_value,iso_var_no[ino1]-iso_var_value,
            x_tmp[idata_cv],x_tmp[idata_bf],x_no[ino0],x_no[ino1]);
        for (int ii = triVec_size0; ii < triVec.size(); ++ii) {
          triVec[ii].second = icv;
        }
      }
    }
  }

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    if ((cv_flag[icv0] >= 0)||(cv_flag[icv1] >= 0)) {
      assert(fa_flag[ifa] >= 0);
      const int idata_fa = fa_flag[ifa]-nno;
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        if (cv_flag[icv0] >= 0) {
          const int idata_cv0 = cv_flag[icv0]-nno;
          const int triVec_size0 = triVec.size();
          tetCutPlane(triVec,
              iso_var_tmp[idata_cv0]-iso_var_value,iso_var_tmp[idata_fa]-iso_var_value,iso_var_no[ino0]-iso_var_value,iso_var_no[ino1]-iso_var_value,
              x_tmp[idata_cv0],x_tmp[idata_fa],x_no[ino0],x_no[ino1]);
          for (int ii = triVec_size0; ii < triVec.size(); ++ii) {
            triVec[ii].second = icv0;
          }
        }
        if (cv_flag[icv1] >= 0) {
          const int idata_cv1 = cv_flag[icv1]-nno;
          const int triVec_size0 = triVec.size();
          tetCutPlane(triVec,
                      iso_var_tmp[idata_cv1]-iso_var_value,iso_var_tmp[idata_fa]-iso_var_value,iso_var_no[ino1]-iso_var_value,iso_var_no[ino0]-iso_var_value,
                      x_tmp[idata_cv1],x_tmp[idata_fa],x_no[ino1],x_no[ino0]);
          for (int ii = triVec_size0; ii < triVec.size(); ++ii) {
            triVec[ii].second = icv1;
          }
        }
      }
    }
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    if (cv_flag[icv0] >= 0) {
      const int idata_cv0 = cv_flag[icv0]-nno;
      assert(fa_flag[ifa] >= 0);
      const int idata_fa = fa_flag[ifa]-nno;
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        const int triVec_size0 = triVec.size();
        tetCutPlane(triVec,
            iso_var_tmp[idata_cv0]-iso_var_value,iso_var_tmp[idata_fa]-iso_var_value,iso_var_no[ino0]-iso_var_value,iso_var_no[ino1]-iso_var_value,
                    x_tmp[idata_cv0],x_tmp[idata_fa],x_no[ino0],x_no[ino1]);
        for (int ii = triVec_size0; ii < triVec.size(); ++ii) {
          triVec[ii].second = icv0;
        }
      }
    }
  }

  delete[] x_tmp;
  delete[] iso_var_tmp;
  delete[] fa_flag;
  delete[] bf_flag;
  delete[] cv_flag;

}

void StaticSolver::buildIsoWithData(vector<pair<SimpleTriWithData,int> >& triVec,const string& iso_var_name,const double iso_var_value,const string& data_name,int * cv_hide_flag) {

  // convert the iso_var_name and data_name into a node-base DN's...

  double * iso_var_no = new double[nno];
  setNoDN(iso_var_no,iso_var_name);

  double * data_no = new double[nno];
  setNoDN(data_no,data_name);

  buildIsoWithData(triVec,iso_var_no,iso_var_value,data_no,cv_hide_flag);

  delete[] data_no;
  delete[] iso_var_no;

}

// pass ctidata's if we already have them
void StaticSolver::buildIsoWithData(vector<pair<SimpleTriWithData,int> >& triVec,const string& iso_var_name,const CtiRegister::CtiData* iso_var,const double iso_var_value,const string& data_name, const CtiRegister::CtiData* data,int * cv_hide_flag) {

  // convert the iso_var and data into a node-base DN's...

  double * iso_var_no = new double[nno];
  setNoDN(iso_var_no,iso_var_name,iso_var);

  double * data_no = new double[nno];
  setNoDN(data_no,data_name,data);

  buildIsoWithData(triVec,iso_var_no,iso_var_value,data_no,cv_hide_flag);

  delete[] data_no;
  delete[] iso_var_no;

}

void StaticSolver::buildIsoWithData(vector<pair<SimpleTriWithData,int> >& triVec,const double * const iso_var_no,const double iso_var_value,const double * const data_no,int * cv_hide_flag) {

  //MiscUtils::dumpRange(iso_var_no,nno,"iso_var_no");
  //MiscUtils::dumpRange(data_no,nno,"data_no");

  // build buffer to hold cv, bf, and fa data that might be necessary to build the iso...

  int * cv_flag = new int[ncv];
  setCvFlagFromIsoVar(cv_flag,iso_var_no,iso_var_value);

  int * bf_flag = new int[nbf];
  int * fa_flag = new int[nfa];
  int ndata = setFlagsToIndexTmpData(cv_flag,bf_flag,fa_flag);

  double * iso_var_tmp = new double[ndata-nno];
  double (*x_tmp)[3] = new double[ndata-nno][3];
  double * wgt_tmp = new double[ndata-nno];
  double * data_tmp = new double[ndata-nno];

  // now loop again and work the tets that touch cvs that contain the isosurface...

  for (int idata = 0; idata < ndata-nno; ++idata) {
    iso_var_tmp[idata] = 0.0;
    FOR_I3 x_tmp[idata][i] = 0.0;
    data_tmp[idata] = 0.0;
    wgt_tmp[idata] = 0.0;
  }

  const bool cvs_hidden = !(cv_hide_flag == NULL);
  FOR_IBF {
    const int icv = cvobf[ibf];
    const bool cv_ignore = (cvs_hidden && cv_hide_flag[icv]);
    if ((!cv_ignore)&&(cv_flag[icv] >= 0)) {
      const int idata_cv = cv_flag[icv]-nno;
      assert(bf_flag[ibf] >= 0);
      const int idata_bf = bf_flag[ibf]-nno;
      double wgt_bf = 0.0;
      int ino1 = noobf_v[noobf_i[ibf+1]-1];
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino0 = ino1;
        ino1 = noobf_v[nob];
        const double wgt_ed = DIST(x_no[ino0],x_no[ino1]);
        iso_var_tmp[idata_bf] += wgt_ed*(iso_var_no[ino0] + iso_var_no[ino1]);
        FOR_I3 x_tmp[idata_bf][i] += wgt_ed*(x_no[ino0][i] + x_no[ino1][i]);
        data_tmp[idata_bf] += wgt_ed*(data_no[ino0] + data_no[ino1]);
        wgt_bf += wgt_ed;
      }
      iso_var_tmp[idata_bf] /= 2.0*wgt_bf;
      FOR_I3 x_tmp[idata_bf][i] /= 2.0*wgt_bf;
      data_tmp[idata_bf] /= 2.0*wgt_bf;
      // also add to the cv...
      wgt_bf = MAG(n_bf[ibf]);
      iso_var_tmp[idata_cv] += wgt_bf*iso_var_tmp[idata_bf];
      FOR_I3 x_tmp[idata_cv][i] += wgt_bf*x_tmp[idata_bf][i];
      data_tmp[idata_cv] += wgt_bf*data_tmp[idata_bf];
      wgt_tmp[idata_cv] += wgt_bf;
    }
  }

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    const bool cv_ignore = (cvs_hidden && ( cv_hide_flag[icv0] || cv_hide_flag[icv1]));
    if ((!cv_ignore) && ( ((cv_flag[icv0] >= 0)||(cv_flag[icv1] >= 0))) ) {
      assert(fa_flag[ifa] >= 0);
      const int idata_fa = fa_flag[ifa]-nno;
      double wgt_fa = 0.0;
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        const double wgt_ed = DIST(x_no[ino0],x_no[ino1]);
        iso_var_tmp[idata_fa] += wgt_ed*(iso_var_no[ino0] + iso_var_no[ino1]);
        FOR_I3 x_tmp[idata_fa][i] += wgt_ed*(x_no[ino0][i] + x_no[ino1][i]);
        data_tmp[idata_fa] += wgt_ed*(data_no[ino0] + data_no[ino1]);
        wgt_fa += wgt_ed;
      }
      iso_var_tmp[idata_fa] /= 2.0*wgt_fa;
      FOR_I3 x_tmp[idata_fa][i] /= 2.0*wgt_fa;
      data_tmp[idata_fa] /= 2.0*wgt_fa;
      // also add to the cv...
      wgt_fa = MAG(n_fa[ifa]);
      if (cv_flag[icv0] >= 0) {
        const int idata_cv0 = cv_flag[icv0]-nno;
        iso_var_tmp[idata_cv0] += wgt_fa*iso_var_tmp[idata_fa];
        FOR_I3 x_tmp[idata_cv0][i] += wgt_fa*x_tmp[idata_fa][i];
        data_tmp[idata_cv0] += wgt_fa*data_tmp[idata_fa];
        wgt_tmp[idata_cv0] += wgt_fa;
      }
      if (cv_flag[icv1] >= 0) {
        const int idata_cv1 = cv_flag[icv1]-nno;
        iso_var_tmp[idata_cv1] += wgt_fa*iso_var_tmp[idata_fa];
        FOR_I3 x_tmp[idata_cv1][i] += wgt_fa*x_tmp[idata_fa][i];
        data_tmp[idata_cv1] += wgt_fa*data_tmp[idata_fa];
        wgt_tmp[idata_cv1] += wgt_fa;
      }
    }
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    const bool cv_ignore = (cvs_hidden && cv_hide_flag[icv0]);
    if ((!cv_ignore)&&(cv_flag[icv0] >= 0)) {
      const int idata_cv0 = cv_flag[icv0]-nno;
      assert(fa_flag[ifa] >= 0);
      const int idata_fa = fa_flag[ifa]-nno;
      double wgt_fa = 0.0;
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        const double wgt_ed = DIST(x_no[ino0],x_no[ino1]);
        iso_var_tmp[idata_fa] += wgt_ed*(iso_var_no[ino0] + iso_var_no[ino1]);
        FOR_I3 x_tmp[idata_fa][i] += wgt_ed*(x_no[ino0][i] + x_no[ino1][i]);
        data_tmp[idata_fa] += wgt_ed*(data_no[ino0] + data_no[ino1]);
        wgt_fa += wgt_ed;
      }
      iso_var_tmp[idata_fa] /= 2.0*wgt_fa;
      FOR_I3 x_tmp[idata_fa][i] /= 2.0*wgt_fa;
      data_tmp[idata_fa] /= 2.0*wgt_fa;
      // also add to the cv...
      wgt_fa = MAG(n_fa[ifa]);
      iso_var_tmp[idata_cv0] += wgt_fa*iso_var_tmp[idata_fa];
      FOR_I3 x_tmp[idata_cv0][i] += wgt_fa*x_tmp[idata_fa][i];
      data_tmp[idata_cv0] += wgt_fa*data_tmp[idata_fa];
      wgt_tmp[idata_cv0] += wgt_fa;
    }
  }

  // we now have variables at all faces, cvs and nodes, so now record the tets...

  FOR_ICV {
    const bool cv_ignore = (cvs_hidden && cv_hide_flag[icv]);
    if ((!cv_ignore)&&(cv_flag[icv] >= 0)) {
      const int idata_cv = cv_flag[icv]-nno;
      assert(wgt_tmp[idata_cv] > 0.0);
      iso_var_tmp[idata_cv] /= wgt_tmp[idata_cv];
      FOR_I3 x_tmp[idata_cv][i] /= wgt_tmp[idata_cv];
      data_tmp[idata_cv] /= wgt_tmp[idata_cv];
    }
  }

  delete[] wgt_tmp;

  FOR_IBF {
    const int icv = cvobf[ibf];
    const bool cv_ignore = (cvs_hidden && cv_hide_flag[icv]);
    if ((!cv_ignore)&&(cv_flag[icv] >= 0)) {
      const int idata_cv = cv_flag[icv]-nno;
      assert(bf_flag[ibf] >= 0);
      const int idata_bf = bf_flag[ibf]-nno;
      int ino1 = noobf_v[noobf_i[ibf+1]-1];
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino0 = ino1;
        ino1 = noobf_v[nob];
        tetCutPlaneData(triVec,
            iso_var_tmp[idata_cv]-iso_var_value,iso_var_tmp[idata_bf]-iso_var_value,iso_var_no[ino0]-iso_var_value,iso_var_no[ino1]-iso_var_value,
            x_tmp[idata_cv],x_tmp[idata_bf],x_no[ino0],x_no[ino1],
            data_tmp[idata_cv],data_tmp[idata_bf],data_no[ino0],data_no[ino1]);
      }
    }
  }

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    const bool cv_ignore = (cvs_hidden && ( cv_hide_flag[icv0]||cv_hide_flag[icv1]));
    if ((!cv_ignore)&&(((cv_flag[icv0] >= 0)||(cv_flag[icv1] >= 0))) ) {
      assert(fa_flag[ifa] >= 0);
      const int idata_fa = fa_flag[ifa]-nno;
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        if (cv_flag[icv0] >= 0) {
          const int idata_cv0 = cv_flag[icv0]-nno;
          tetCutPlaneData(triVec,
              iso_var_tmp[idata_cv0]-iso_var_value,iso_var_tmp[idata_fa]-iso_var_value,iso_var_no[ino0]-iso_var_value,iso_var_no[ino1]-iso_var_value,
              x_tmp[idata_cv0],x_tmp[idata_fa],x_no[ino0],x_no[ino1],
              data_tmp[idata_cv0],data_tmp[idata_fa],data_no[ino0],data_no[ino1]);

        }
        if (cv_flag[icv1] >= 0) {
          const int idata_cv1 = cv_flag[icv1]-nno;
          tetCutPlaneData(triVec,
              iso_var_tmp[idata_cv1]-iso_var_value,iso_var_tmp[idata_fa]-iso_var_value,iso_var_no[ino1]-iso_var_value,iso_var_no[ino0]-iso_var_value,
              x_tmp[idata_cv1],x_tmp[idata_fa],x_no[ino1],x_no[ino0],
              data_tmp[idata_cv1],data_tmp[idata_fa],data_no[ino1],data_no[ino0]);
        }
      }
    }
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    const bool cv_ignore = (cvs_hidden && cv_hide_flag[icv0]);
    if ((!cv_ignore)&&(cv_flag[icv0] >= 0)) {
      const int idata_cv0 = cv_flag[icv0]-nno;
      assert(fa_flag[ifa] >= 0);
      const int idata_fa = fa_flag[ifa]-nno;
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        tetCutPlaneData(triVec,
            iso_var_tmp[idata_cv0]-iso_var_value,iso_var_tmp[idata_fa]-iso_var_value,iso_var_no[ino0]-iso_var_value,iso_var_no[ino1]-iso_var_value,
            x_tmp[idata_cv0],x_tmp[idata_fa],x_no[ino0],x_no[ino1],
            data_tmp[idata_cv0],data_tmp[idata_fa],data_no[ino0],data_no[ino1]);
      }
    }
  }

  delete[] data_tmp;
  delete[] x_tmp;
  delete[] iso_var_tmp;
  delete[] fa_flag;
  delete[] bf_flag;
  delete[] cv_flag;

}

void StaticSolver::buildIsoWithWgtsAndData(vector<SimpleTriWithWgts>& triVec,const double * const iso_var_no,const double iso_var_value,double *&dn_tmp, double *&dn_no, const int nscalars, double (*&dn3_tmp)[3], double (*&dn3_no)[3], const int nvectors,int * cv_hide_flag) {

  MiscUtils::dumpRange(iso_var_no,nno,"iso_var_no");
  for (int ivar = 0; ivar < nscalars; ++ivar) {
    stringstream ss;
    ss << ivar;
    string name = "dn_no[" + ss.str() + "]";
    MiscUtils::dumpRange(dn_no+ivar*nno,nno,name);
  }
  for (int ivar = 0; ivar < nvectors; ++ivar) {
    stringstream ss;
    ss << ivar;
    string name = "dn3_no[" + ss.str() + "]";
    MiscUtils::dumpRange(dn3_no+ivar*nno,nno,name);
  }

  // shuffle dn_no to group collocated data
  double *_dn_no = dn_no;
  if (nscalars > 0) {
    assert(dn_no != NULL); // already allocated/filled
    dn_no = new double[nno*nscalars];
    FOR_INO {
      for (int ivar = 0; ivar < nscalars; ++ivar) {
        dn_no[ino*nscalars+ivar] = _dn_no[nno*ivar+ino];
      }
    }
  }
  if (_dn_no != NULL) delete[] _dn_no;
  // shuffle dn_no to group collocated data
  double (*_dn3_no)[3] = dn3_no;
  if (nvectors > 0) {
    assert(dn3_no != NULL); // already allocated/filled
    dn3_no = new double[nno*nvectors][3];
    FOR_INO {
      for (int ivar = 0; ivar < nvectors; ++ivar) {
        FOR_I3 dn3_no[ino*nvectors+ivar][i] = _dn3_no[nno*ivar+ino][i];
      }
    }
  }
  if (_dn3_no != NULL) delete[] _dn3_no;


  // build buffer to hold cv, bf, and fa data that might be necessary to build the iso...

  int * cv_flag = new int[ncv];
  setCvFlagFromIsoVar(cv_flag,iso_var_no,iso_var_value);

  int * bf_flag = new int[nbf];
  int * fa_flag = new int[nfa];
  int ndata = setFlagsToIndexTmpData(cv_flag,bf_flag,fa_flag);

  // now loop again and work the tets that touch cvs that contain the isosurface...

  double * iso_var_tmp = new double[ndata-nno];
  double (*x_tmp)[3] = new double[ndata-nno][3];
  double * wgt_tmp = new double[ndata-nno];
  assert(dn_tmp == NULL); dn_tmp = new double[(ndata-nno)*nscalars];
  assert(dn3_tmp == NULL); dn3_tmp = new double[(ndata-nno)*nvectors][3];

  for (int idata = 0; idata < ndata-nno; ++idata) {
    iso_var_tmp[idata] = 0.0;
    FOR_I3 x_tmp[idata][i] = 0.0;
    wgt_tmp[idata] = 0.0;
    for (int ivar = 0; ivar < nscalars; ++ivar)
      dn_tmp[idata*nscalars+ivar] = 0.0;
    for (int ivar = 0; ivar < nvectors; ++ivar)
      FOR_I3 dn3_tmp[idata*nvectors+ivar][i] = 0.0;
  }

  const bool cvs_hidden = !(cv_hide_flag == NULL);
  FOR_IBF {
    const int icv = cvobf[ibf];
    const bool cv_ignore = (cvs_hidden && cv_hide_flag[icv]);
    if ((!cv_ignore)&&(cv_flag[icv] >= 0)) {
      const int idata_cv = cv_flag[icv]-nno;
      assert(bf_flag[ibf] >= 0);
      const int idata_bf = bf_flag[ibf]-nno;
      double wgt_bf = 0.0;
      int ino1 = noobf_v[noobf_i[ibf+1]-1];
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino0 = ino1;
        ino1 = noobf_v[nob];
        const double wgt_ed = DIST(x_no[ino0],x_no[ino1]);
        iso_var_tmp[idata_bf] += wgt_ed*(iso_var_no[ino0] + iso_var_no[ino1]);
        for (int ivar = 0; ivar < nscalars; ++ivar)
          dn_tmp[idata_bf*nscalars+ivar] += wgt_ed*(dn_no[ino0*nscalars+ivar] + dn_no[ino1*nscalars+ivar]);
        for (int ivar = 0; ivar < nvectors; ++ivar)
          FOR_I3 dn3_tmp[idata_bf*nvectors+ivar][i] += wgt_ed*(dn3_no[ino0*nvectors+ivar][i] + dn3_no[ino1*nvectors+ivar][i]);
        FOR_I3 x_tmp[idata_bf][i] += wgt_ed*(x_no[ino0][i] + x_no[ino1][i]);
        wgt_bf += wgt_ed;
      }
      iso_var_tmp[idata_bf] /= 2.0*wgt_bf;
      for (int ivar = 0; ivar < nscalars; ++ivar)
        dn_tmp[idata_bf*nscalars+ivar] /= 2.0*wgt_bf;
      for (int ivar = 0; ivar < nvectors; ++ivar)
        FOR_I3 dn3_tmp[idata_bf*nvectors+ivar][i] /= 2.0*wgt_bf;
      FOR_I3 x_tmp[idata_bf][i] /= 2.0*wgt_bf;
      // also add to the cv...
      wgt_bf = MAG(n_bf[ibf]);
      iso_var_tmp[idata_cv] += wgt_bf*iso_var_tmp[idata_bf];
      for (int ivar = 0; ivar < nscalars; ++ivar)
        dn_tmp[idata_cv*nscalars+ivar] += wgt_bf*dn_tmp[idata_bf*nscalars+ivar];
      for (int ivar = 0; ivar < nvectors; ++ivar)
        FOR_I3 dn3_tmp[idata_cv*nvectors+ivar][i] += wgt_bf*dn3_tmp[idata_bf*nvectors+ivar][i];
      FOR_I3 x_tmp[idata_cv][i] += wgt_bf*x_tmp[idata_bf][i];
      wgt_tmp[idata_cv] += wgt_bf;
    }
  }

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    const bool cv_ignore = (cvs_hidden && (cv_hide_flag[icv0] || cv_hide_flag[icv1]));
    if ((!cv_ignore)&&(((cv_flag[icv0] >= 0)||(cv_flag[icv1] >= 0))) ) {
      assert(fa_flag[ifa] >= 0);
      const int idata_fa = fa_flag[ifa]-nno;
      double wgt_fa = 0.0;
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        const double wgt_ed = DIST(x_no[ino0],x_no[ino1]);
        iso_var_tmp[idata_fa] += wgt_ed*(iso_var_no[ino0] + iso_var_no[ino1]);
        for (int ivar = 0; ivar < nscalars; ++ivar)
          dn_tmp[idata_fa*nscalars+ivar] += wgt_ed*(dn_no[ino0*nscalars+ivar] + dn_no[ino1*nscalars+ivar]);
        for (int ivar = 0; ivar < nvectors; ++ivar)
          FOR_I3 dn3_tmp[idata_fa*nvectors+ivar][i] += wgt_ed*(dn3_no[ino0*nvectors+ivar][i] + dn3_no[ino1*nvectors+ivar][i]);
        FOR_I3 x_tmp[idata_fa][i] += wgt_ed*(x_no[ino0][i] + x_no[ino1][i]);
        wgt_fa += wgt_ed;
      }
      iso_var_tmp[idata_fa] /= 2.0*wgt_fa;
      for (int ivar = 0; ivar < nscalars; ++ivar)
        dn_tmp[idata_fa*nscalars+ivar] /= 2.0*wgt_fa;
      for (int ivar = 0; ivar < nvectors; ++ivar)
        FOR_I3 dn3_tmp[idata_fa*nvectors+ivar][i] /= 2.0*wgt_fa;
      FOR_I3 x_tmp[idata_fa][i] /= 2.0*wgt_fa;
      // also add to the cv...
      wgt_fa = MAG(n_fa[ifa]);
      if (cv_flag[icv0] >= 0) {
        const int idata_cv0 = cv_flag[icv0]-nno;
        iso_var_tmp[idata_cv0] += wgt_fa*iso_var_tmp[idata_fa];
        for (int ivar = 0; ivar < nscalars; ++ivar)
          dn_tmp[idata_cv0*nscalars+ivar] += wgt_fa*dn_tmp[idata_fa*nscalars+ivar];
        for (int ivar = 0; ivar < nvectors; ++ivar)
          FOR_I3 dn3_tmp[idata_cv0*nvectors+ivar][i] += wgt_fa*dn3_tmp[idata_fa*nvectors+ivar][i];
        FOR_I3 x_tmp[idata_cv0][i] += wgt_fa*x_tmp[idata_fa][i];
        wgt_tmp[idata_cv0] += wgt_fa;
      }
      if (cv_flag[icv1] >= 0) {
        const int idata_cv1 = cv_flag[icv1]-nno;
        iso_var_tmp[idata_cv1] += wgt_fa*iso_var_tmp[idata_fa];
        for (int ivar = 0; ivar < nscalars; ++ivar)
          dn_tmp[idata_cv1*nscalars+ivar] += wgt_fa*dn_tmp[idata_fa*nscalars+ivar];
        for (int ivar = 0; ivar < nvectors; ++ivar)
          FOR_I3 dn3_tmp[idata_cv1*nvectors+ivar][i] += wgt_fa*dn3_tmp[idata_fa*nvectors+ivar][i];
        FOR_I3 x_tmp[idata_cv1][i] += wgt_fa*x_tmp[idata_fa][i];
        wgt_tmp[idata_cv1] += wgt_fa;
      }
    }
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    const bool cv_ignore = (cvs_hidden && cv_hide_flag[icv0]);
    if (!(cv_ignore)&&(cv_flag[icv0] >= 0)) {
      const int idata_cv0 = cv_flag[icv0]-nno;
      assert(fa_flag[ifa] >= 0);
      const int idata_fa = fa_flag[ifa]-nno;
      double wgt_fa = 0.0;
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        const double wgt_ed = DIST(x_no[ino0],x_no[ino1]);
        iso_var_tmp[idata_fa] += wgt_ed*(iso_var_no[ino0] + iso_var_no[ino1]);
        for (int ivar = 0; ivar < nscalars; ++ivar)
          dn_tmp[idata_fa*nscalars+ivar] += wgt_ed*(dn_no[ino0*nscalars+ivar] + dn_no[ino1*nscalars+ivar]);
        for (int ivar = 0; ivar < nvectors; ++ivar)
          FOR_I3 dn3_tmp[idata_fa*nvectors+ivar][i] += wgt_ed*(dn3_no[ino0*nvectors+ivar][i] + dn3_no[ino1*nvectors+ivar][i]);
        FOR_I3 x_tmp[idata_fa][i] += wgt_ed*(x_no[ino0][i] + x_no[ino1][i]);
        wgt_fa += wgt_ed;
      }
      iso_var_tmp[idata_fa] /= 2.0*wgt_fa;
      for (int ivar = 0; ivar < nscalars; ++ivar)
        dn_tmp[idata_fa*nscalars+ivar] /= 2.0*wgt_fa;
      for (int ivar = 0; ivar < nvectors; ++ivar)
        FOR_I3 dn3_tmp[idata_fa*nvectors+ivar][i] /= 2.0*wgt_fa;
      FOR_I3 x_tmp[idata_fa][i] /= 2.0*wgt_fa;
      // also add to the cv...
      wgt_fa = MAG(n_fa[ifa]);
      iso_var_tmp[idata_cv0] += wgt_fa*iso_var_tmp[idata_fa];
      for (int ivar = 0; ivar < nscalars; ++ivar)
        dn_tmp[idata_cv0*nscalars+ivar] += wgt_fa*dn_tmp[idata_fa*nscalars+ivar];
      for (int ivar = 0; ivar < nvectors; ++ivar)
        FOR_I3 dn3_tmp[idata_cv0*nvectors+ivar][i] += wgt_fa*dn3_tmp[idata_fa*nvectors+ivar][i];
      FOR_I3 x_tmp[idata_cv0][i] += wgt_fa*x_tmp[idata_fa][i];
      wgt_tmp[idata_cv0] += wgt_fa;
    }
  }

  // we now have variables at all faces, cvs and nodes, so now record the tets...

  FOR_ICV {
    const bool cv_ignore = (cvs_hidden && cv_hide_flag[icv]);
    if ((!cv_ignore)&&(cv_flag[icv] >= 0)) {
      const int idata_cv = cv_flag[icv]-nno;
      assert(wgt_tmp[idata_cv] > 0.0);
      iso_var_tmp[idata_cv] /= wgt_tmp[idata_cv];
      for (int ivar = 0; ivar < nscalars; ++ivar)
        dn_tmp[idata_cv*nscalars+ivar] /= wgt_tmp[idata_cv];
      for (int ivar = 0; ivar < nvectors; ++ivar)
        FOR_I3 dn3_tmp[idata_cv*nvectors+ivar][i] /= wgt_tmp[idata_cv];
      FOR_I3 x_tmp[idata_cv][i] /= wgt_tmp[idata_cv];
    }
  }
  delete[] wgt_tmp;

  FOR_IBF {
    const int icv = cvobf[ibf];
    const bool cv_ignore = (cvs_hidden && cv_hide_flag[icv]);
    if ((!cv_ignore)&&(cv_flag[icv] >= 0)) {
      const int idata_cv = cv_flag[icv]-nno;
      assert(bf_flag[ibf] >= 0);
      const int idata_bf = bf_flag[ibf]-nno;
      int ino1 = noobf_v[noobf_i[ibf+1]-1];
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino0 = ino1;
        ino1 = noobf_v[nob];
        tetCutPlaneWgts(triVec,
            iso_var_tmp[idata_cv]-iso_var_value,iso_var_tmp[idata_bf]-iso_var_value,iso_var_no[ino0]-iso_var_value,iso_var_no[ino1]-iso_var_value,
            x_tmp[idata_cv],x_tmp[idata_bf],x_no[ino0],x_no[ino1],
            cv_flag[icv],bf_flag[ibf],ino0,ino1);
      }
    }
  }

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    const bool cv_ignore = (cvs_hidden && (cv_hide_flag[icv0] || cv_hide_flag[icv1]));
    if ((!cv_ignore)&&(((cv_flag[icv0] >= 0)||(cv_flag[icv1] >= 0))) ) {
      assert(fa_flag[ifa] >= 0);
      const int idata_fa = fa_flag[ifa]-nno;
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        if (cv_flag[icv0] >= 0) {
          const int idata_cv0 = cv_flag[icv0]-nno;
          tetCutPlaneWgts(triVec,
              iso_var_tmp[idata_cv0]-iso_var_value,iso_var_tmp[idata_fa]-iso_var_value,iso_var_no[ino0]-iso_var_value,iso_var_no[ino1]-iso_var_value,
              x_tmp[idata_cv0],x_tmp[idata_fa],x_no[ino0],x_no[ino1],
              cv_flag[icv0],fa_flag[ifa],ino0,ino1);
        }
        if (cv_flag[icv1] >= 0) {
          const int idata_cv1 = cv_flag[icv1]-nno;
          tetCutPlaneWgts(triVec,
              iso_var_tmp[idata_cv1]-iso_var_value,iso_var_tmp[idata_fa]-iso_var_value,iso_var_no[ino0]-iso_var_value,iso_var_no[ino1]-iso_var_value,
              x_tmp[idata_cv1],x_tmp[idata_fa],x_no[ino0],x_no[ino1],
              cv_flag[icv1],fa_flag[ifa],ino0,ino1);
        }
      }
    }
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    const bool cv_ignore = (cvs_hidden && cv_hide_flag[icv0]);
    if ((!cv_ignore)&&(cv_flag[icv0] >= 0)) {
      const int idata_cv0 = cv_flag[icv0]-nno;
      assert(fa_flag[ifa] >= 0);
      const int idata_fa = fa_flag[ifa]-nno;
      int ino1 = noofa_v[noofa_i[ifa+1]-1];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noofa_v[nof];
        tetCutPlaneWgts(triVec,
            iso_var_tmp[idata_cv0]-iso_var_value,iso_var_tmp[idata_fa]-iso_var_value,iso_var_no[ino0]-iso_var_value,iso_var_no[ino1]-iso_var_value,
            x_tmp[idata_cv0],x_tmp[idata_fa],x_no[ino0],x_no[ino1],
            cv_flag[icv0],fa_flag[ifa],ino0,ino1);
      }
    }
  }

  delete[] x_tmp;
  delete[] iso_var_tmp;
  delete[] fa_flag;
  delete[] bf_flag;
  delete[] cv_flag;

}

void StaticSolver::buildChtIso(vector<pair<SimpleTri,int> >& triVec,const CtiData * iso_var_no,const double iso_var_value) {
  int (*noote)[4] = NULL;
  int nte;
  const bool b_got_tets = iso_var_no->getTets(noote,nte);
  assert(b_got_tets);
  assert(noote != NULL);
  double (*x_no_te)[3] = NULL;
  const bool b_got_x = iso_var_no->getX(x_no_te);
  assert(x_no_te != NULL);
  assert(iso_var_no->getType() == DN_DATA);
  double * data_no = iso_var_no->getDNptr();

  const int tstart=triVec.size();
  for (int ite=0; ite<nte; ++ite) {
    chtTetCutPlane(triVec,
      data_no[noote[ite][0]]-iso_var_value,
      data_no[noote[ite][1]]-iso_var_value,
      data_no[noote[ite][2]]-iso_var_value,
      data_no[noote[ite][3]]-iso_var_value,
      x_no_te[noote[ite][0]],
      x_no_te[noote[ite][1]],
      x_no_te[noote[ite][2]],
      x_no_te[noote[ite][3]]
    );
  }

  // adjust pair int so that all mesh edges draw
  // if (triVec.size() > tstart) {
  //   for (int ite=tstart, tend=triVec.size(); ite<tend; ++ite) {
  //     triVec[ite].second=7;  // naive draw all edges; eventually need to go into the tetCutPlane and return appropriate edges to draw to account for two-tri pushes with internal edge that shouldn't be highlighted
  //   }
  // }
}

void StaticSolver::buildChtIsoWithData(vector<pair<SimpleTriWithData,int> >& triVec,const CtiData * iso_var_no,const double iso_var_value,const CtiData * data_var_no) {
  int (*noote)[4] = NULL;
  int nte;
  const bool b_got_tets = iso_var_no->getTets(noote,nte);
  assert(b_got_tets);
  assert(noote != NULL);
  double (*x_no_te)[3] = NULL;
  const bool b_got_x = iso_var_no->getX(x_no_te);
  assert(x_no_te != NULL);
  assert(iso_var_no->getType() == DN_DATA);
  double * data_no = iso_var_no->getDNptr();
  double * data_var = data_var_no->getDNptr();

  const int tstart=triVec.size();
  for (int ite=0; ite<nte; ++ite) {
    chtTetCutPlaneData(triVec,
      data_no[noote[ite][0]]-iso_var_value,
      data_no[noote[ite][1]]-iso_var_value,
      data_no[noote[ite][2]]-iso_var_value,
      data_no[noote[ite][3]]-iso_var_value,
      x_no_te[noote[ite][0]],
      x_no_te[noote[ite][1]],
      x_no_te[noote[ite][2]],
      x_no_te[noote[ite][3]],
      data_var[noote[ite][0]],
      data_var[noote[ite][1]],
      data_var[noote[ite][2]],
      data_var[noote[ite][3]]
    );
  }

  // adjust pair int so that all mesh edges draw
  // if (triVec.size() > tstart) {
  //   for (int ite=tstart, tend=triVec.size(); ite<tend; ++ite) {
  //     triVec[ite].second=7;  // naive draw all edges; eventually need to go into the tetCutPlane and return appropriate edges to draw to account for two-tri pushes with internal edge that shouldn't be highlighted
  //   }
  // }
}

void StaticSolver::setBfDN(double * no_data,const string& name,double *bf_data) {
  CtiRegister::CtiData * var = CtiRegister::getCtiData(name);
  setBfDN(no_data,name,var,bf_data);
}

void StaticSolver::setBfDN(double * no_data,const string& name,const CtiRegister::CtiData * var,double *bf_data) {
  if (var==NULL){
    CERR("name \"" << name << "\" cannot be evaluated");
  }
  else {
    if (no_data != NULL)
      averageBfToNo(no_data,var);
    if (bf_data != NULL) {
      if (var->getType() == DN_DATA) {
        const int topo = var->getUnindexedTopology();
        if (topo == BF_DATA) {
          const int iz = var->getIndex();
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            bf_data[ibf] = var->dn(ibf-bfZoneVec[iz].ibf_f);;
          }
        }
        else if (topo == CV_DATA) {
          const int iz = var->getIndex();
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            const int icv = cvobf[ibf];
            bf_data[ibf] = var->dn(icv);
          }
        }
        else {
          CERR("name \"" << name << "\" is not BF/CV DATA");
        }
      }
      else if (var->getType() == IN_DATA) {
        const int topo = var->getUnindexedTopology();
        if (topo == BF_DATA) {
          const int iz = var->getIndex();
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            bf_data[ibf] = double(var->in(ibf-bfZoneVec[iz].ibf_f));;
          }
        }
        else if (topo == CV_DATA) {
          const int iz = var->getIndex();
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            const int icv = cvobf[ibf];
            bf_data[ibf] = double(var->in(icv));
          }
        }
        else {
          CERR("name \"" << name << "\" is not BF/CV DATA");
        }
      }
      else {
        CERR("name \"" << name << "\" cannot be evaluated to BF DN");
      }
    }
  }
}

void StaticSolver::setBfDN(double * no_data,const vector<string>& name_vec,double * bf_data) {
  const int nvar = name_vec.size();
  vector<CtiRegister::CtiData*> var_vec(nvar);
  for (int ivar = 0; ivar < nvar; ++ivar)
    var_vec[ivar] = CtiRegister::getCtiData(name_vec[ivar]);
  setBfDN(no_data,name_vec,var_vec,bf_data);
}

void StaticSolver::setBfDN(double * no_data,const vector<string>& name_vec,const vector<CtiRegister::CtiData *>& var_vec,double *bf_data) {

  const int nvar = name_vec.size();
  for (int ivar = 0; ivar < nvar; ++ivar) {
    if (var_vec[ivar] == NULL) {
      CERR("name \"" << name_vec[ivar] << "\" cannot be evaluated");
    }
  }
  if (no_data != NULL)
    averageBfToNo(no_data,var_vec);
  if (bf_data != NULL) {
    for (int ivar = 0; ivar < nvar; ++ivar) {
      CtiRegister::CtiData* var = var_vec[ivar];
      if (var->getType() == DN_DATA) {
        const int topo = var->getUnindexedTopology();
        if (topo == BF_DATA) {
          const int iz = var->getIndex();
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            bf_data[ibf] = var->dn(ibf-bfZoneVec[iz].ibf_f);
          }
        }
        else {
          CERR("name \"" << name_vec[ivar] << "\" is not BF_DATA");
        }
      }
      else if (var->getType() == IN_DATA) {
        const int topo = var->getUnindexedTopology();
        if (topo == BF_DATA) {
          const int iz = var->getIndex();
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            bf_data[ibf] = double(var->in(ibf-bfZoneVec[iz].ibf_f));;
          }
        }
        else {
          CERR("name \"" << name_vec[ivar] << "\" is not BF_DATA");
        }
      }
      else {
        CERR("name \"" << name_vec[ivar] << "\" cannot be evaluated to BF DN");
      }
    }
  }
}

void StaticSolver::setBfDN3(double (*no_data)[3],const string& name,double (*bf_data)[3]) {
  CtiRegister::CtiData * var = CtiRegister::getCtiData(name);
  setBfDN3(no_data,name,var,bf_data);
}

void StaticSolver::setBfDN3(double (*no_data)[3],const string& name,const CtiRegister::CtiData * var,double (*bf_data)[3]) {

  if (var==NULL){
    CERR("name \"" << name << "\" cannot be evaluated");
  }
  else {
    if (no_data != NULL)
      averageBfToNo(no_data,var);
    if (bf_data != NULL) {
      if (var->getType() == DN3_DATA) {
        const int topo = var->getUnindexedTopology();
        if (topo == BF_DATA) {
          const int iz = var->getIndex();
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            FOR_I3 bf_data[ibf][i] = var->dn3(ibf-bfZoneVec[iz].ibf_f,i);
          }
        }
        else if (topo == CV_DATA) {
          FOR_IBF {
            const int icv = cvobf[ibf];
            FOR_I3 bf_data[ibf][i] = var->dn3(icv,i);
          }
        }
        else {
          CERR("name \"" << name << "\" is not BF/CV DATA");
        }
      }
      else {
        CERR("name \"" << name << "\" cannot be evaluated to BF DN3");
      }
    }
  }
}

void StaticSolver::setBfDN3(double (*no_data)[3],const vector<string>& name_vec,double (*bf_data)[3]) {
  const int nvar = name_vec.size();
  vector<CtiRegister::CtiData*> var_vec(nvar);
  for (int ivar = 0; ivar < nvar; ++ivar)
    var_vec[ivar] = CtiRegister::getCtiData(name_vec[ivar]);
  setBfDN3(no_data,name_vec,var_vec,bf_data);
}

void StaticSolver::setBfDN3(double (*no_data)[3],const vector<string>& name_vec,const vector<CtiRegister::CtiData *>& var_vec,double (*bf_data)[3]) {

  const int nvar = name_vec.size();
  for (int ivar = 0; ivar < nvar; ++ivar) {
    if (var_vec[ivar] == NULL) {
      CERR("name \"" << name_vec[ivar] << "\" cannot be evaluated");
    }
  }
  if (no_data != NULL)
    averageBfToNo(no_data,var_vec);
  if (bf_data != NULL) {
    for (int ivar = 0; ivar < nvar; ++ivar) {
      CtiRegister::CtiData* var = var_vec[ivar];
      if (var->getType() == DN3_DATA) {
        const int topo = var->getUnindexedTopology();
        if (topo == BF_DATA) {
          const int iz = var->getIndex();
          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
            assert(zone_bf[ibf] == iz);
            FOR_I3 bf_data[ibf][i] = var->dn3(ibf-bfZoneVec[iz].ibf_f,i);
          }
        }
        else {
          CERR("name \"" << name_vec[ivar] << "\" is not BF_DATA");
        }
      }
      else {
        CERR("name \"" << name_vec[ivar] << "\" cannot be evaluated to BF DN3");
      }
    }
  }
}

void StaticSolver::setNoDN(double * no_data,const string& name) {
  CtiRegister::CtiData * var = CtiRegister::getCtiData(name);
  setNoDN(no_data,name,var);
}

void StaticSolver::setNoDN(double * no_data,const string& name,const CtiRegister::CtiData * var) {

  if (var==NULL){
    CERR("name \"" << name << "\" cannot be evaluated");
  }
  else if (var->getType() == I_DATA) {
    double val = double(var->i());
    FOR_INO no_data[ino] = val;
  }
  else if (var->getType() == D_DATA) {
    double val = var->d();
    FOR_INO no_data[ino] = val;
  }
  else if (var->getType() == DN_DATA) {
    const int topo = var->getUnindexedTopology();
    if (topo == CV_DATA) {
      averageCvToNo(no_data,var);
    }
    else if (topo == NO_DATA) {
      FOR_INO no_data[ino] = var->dn(ino);
    }
    else {
      CERR("name \"" << name << "\" is not CV_DATA or NO_DATA");
    }
  }
  else if (var->getType() == IN_DATA) {
    const int topo = var->getUnindexedTopology();
    if (topo == CV_DATA) {
      averageCvToNo(no_data,var); // can handle IN_DATA too
    }
    else if (topo == NO_DATA) {
      FOR_INO no_data[ino] = double(var->dn(ino));
    }
    else {
      CERR("name \"" << name << "\" is not CV_DATA or NO_DATA");
    }
  }
  else {
    CERR("name \"" << name << "\" cannot be evaluated to NO DN");
  }

}

void StaticSolver::setNoDN3(double (*no_data)[3],const string& name) {

  CtiRegister::CtiData * var = CtiRegister::getCtiData(name);
  if (var==NULL){
    CERR("name \"" << name << "\" cannot be evaluated");
  }
  else if (var->getType() == DN3_DATA) {
    const int topo = var->getUnindexedTopology();
    if (topo == CV_DATA) {
      averageCvToNo(no_data,var);
    }
    else if (topo == NO_DATA) {
      FOR_INO FOR_I3 no_data[ino][i] = var->dn3(ino,i);
    }
    else {
      CERR("name \"" << name << "\" is not CV_DATA or NO_DATA");
    }
  }
  else {
    CERR("name \"" << name << "\" cannot be evaluated to NO DN3");
  }

}

void StaticSolver::setCvFlagFromIsoVar(int * cv_flag,const double * iso_var_no,const double iso_var_value) {

  // at this point, we want to slice the tets associated with the
  // voronoi diagram. To figure out which tets to slice, we can
  // identify cvs that have node values on both sides of the isosurface...

  FOR_ICV cv_flag[icv] = 0;

  FOR_IBF {
    const int icv = cvobf[ibf];
    for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
      const int ino = noobf_v[nob];
      if (iso_var_no[ino] < iso_var_value) {
        cv_flag[icv] |= 1;
      }
      else if (iso_var_no[ino] > iso_var_value) {
        cv_flag[icv] |= 2;
      }
      else {
        // equals zero!
        cv_flag[icv] = 3;
        break;
      }
      if (cv_flag[icv] == 3)
        break;
    }
  }

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      const int ino = noofa_v[nof];
      if (iso_var_no[ino] < iso_var_value) {
        cv_flag[icv0] |= 1;
        cv_flag[icv1] |= 1;
      }
      else if (iso_var_no[ino] > iso_var_value) {
        cv_flag[icv0] |= 2;
        cv_flag[icv1] |= 2;
      }
      else {
        // equals zero!
        cv_flag[icv0] = 3;
        cv_flag[icv1] = 3;
        break;
      }
      if ((cv_flag[icv0] == 3)&&(cv_flag[icv1] == 3))
        break;
    }
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      const int ino = noofa_v[nof];
      if (iso_var_no[ino] < iso_var_value) {
        cv_flag[icv0] |= 1;
      }
      else if (iso_var_no[ino] > iso_var_value) {
        cv_flag[icv0] |= 2;
      }
      else {
        // equals zero!
        cv_flag[icv0] = 3;
        break;
      }
      if (cv_flag[icv0] == 3)
        break;
    }
  }
}

int StaticSolver::setFlagsToIndexTmpData(int * cv_flag,int * bf_flag,int * fa_flag) {

  // count the size of nno + nbf + nfa + ncv for the data that potentially
  // contribute to the isosurface...

  int ndata = nno;
  FOR_ICV {
    if (cv_flag[icv] == 3) {
      cv_flag[icv] = ndata++;
    }
    else {
      cv_flag[icv] = -1;
    }
  }

  FOR_IBF {
    const int icv = cvobf[ibf];
    if (cv_flag[icv] >= 0) {
      bf_flag[ibf] = ndata++;
    }
    else {
      bf_flag[ibf] = -1;
    }
  }

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    if ((cv_flag[icv0] >= 0)||(cv_flag[icv1] >= 0)) {
      fa_flag[ifa] = ndata++;
    }
    else {
      fa_flag[ifa] = -1;
    }
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    if (cv_flag[icv0] >= 0) {
      fa_flag[ifa] = ndata++;
    }
    else {
      fa_flag[ifa] = -1;
    }
  }

  return ndata;
}

void StaticSolver::flagCvsAdjacentToZones(int8* cv_flag,bool* bfzone_flag,const int index) {
  // for flagged zones, flag the adjacent cvs
  for (map<const string,int>::const_iterator it = bfZoneNameMap.begin(); it != bfZoneNameMap.end(); ++it) {
    const int iz = it->second; // zone index
    cout << "rank " << mpi_rank << ": zone,flag: " << it->first << " " << bfzone_flag[iz] << endl;
    if (bfzone_flag[iz]) {
      for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) {
        assert(zone_bf[ibf] == iz);
        cv_flag[cvobf[ibf]] = index;
      }
    }
  }
}

// flag cvs using signed-distance computed at nodes
void StaticSolver::flagCvsUnderIso(int8* cv_flag, const double *sd_no,const int index) {

  for (int icv = 0; icv < ncv; ++icv) {
    bool outside = false;
    for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
      const int ifa = faocv_v[foc];
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1];++nof) {
        const int ino = noofa_v[nof];
        if (sd_no[ino] > 1.0E-6*r_vv[icv]) { // add a lil eps based on voronoi diagram
          outside = true;
          break;
        }
      }
      if (outside) break;
    }
    // if all points are inside surface (+eps), flag cell
    if (!outside)
      cv_flag[icv] = index;
  }

}

void StaticSolver::initBfZoneVec(StripedMesh * sm) {

  assert(bfZoneVec.empty());
  assert(bfZoneNameMap.empty());

  bfZoneVec.resize(sm->bfZoneVec.size());
  FOR_IZONE(bfZoneVec) {
    // the bfZoneNameMap allows reverse lookup of the zone info...
    //bfZoneNameMap[sm->bfZoneVec[izone].name] = izone;
    bfZoneVec[izone].setName(sm->bfZoneVec[izone].name);
    bfZoneNameMap[bfZoneVec[izone].getName()] = izone;
    // set the index -- sometimes it is helpful to know what zone you are
    // (for example, data registration)...
    assert(bfZoneVec[izone].index == -1); bfZoneVec[izone].index = izone;
    // copy over global indexing information...
    bfZoneVec[izone].nbf_global   = sm->bfZoneVec[izone].nbf_global;
    bfZoneVec[izone].ibf_f_global = sm->bfZoneVec[izone].ibf_f_global;
    // and global geometry...
    bfZoneVec[izone].area_global              = sm->bfZoneVec[izone].area_global;
    bfZoneVec[izone].area_over_delta_global   = sm->bfZoneVec[izone].area_over_delta_global;
    // outwarn normal with area magnitude (but not neccessarily area_global unless planar)...
    bfZoneVec[izone].n_global[0]              = sm->bfZoneVec[izone].n_global[0];
    bfZoneVec[izone].n_global[1]              = sm->bfZoneVec[izone].n_global[1];
    bfZoneVec[izone].n_global[2]              = sm->bfZoneVec[izone].n_global[2];
    // area-weighted center of mass...
    bfZoneVec[izone].x_global[0]              = sm->bfZoneVec[izone].x_global[0];
    bfZoneVec[izone].x_global[1]              = sm->bfZoneVec[izone].x_global[1];
    bfZoneVec[izone].x_global[2]              = sm->bfZoneVec[izone].x_global[2];
  }

}

void StaticSolver::loadBalance(StripedMesh * sm) {

  COUT1("StaticSolver::loadBalance()");

  // at this point the *lb_cost int8's should be set, and we are going to
  // distribute the cvs to balance this agregate cost. To prevent the load
  // balancer from failing, there needs to be some cost...

  assert((cv_lb_cost > 0)||(fa_lb_cost > 0)||(ef_lb_cost > 0));

  // report lb cost model settings...

  if (mpi_rank == 0) {
    cout << " > cv_lb_cost: " << cv_lb_cost << " fa_lb_cost: " << fa_lb_cost << " ef_lb_cost: " << ef_lb_cost << endl;
    FOR_IZONE(bfZoneVec)
      cout << " > lb_cost for zone \"" << bfZoneVec[izone].getName() << "\" " << bfZoneVec[izone].lb_cost << endl;
    for (int ii = 0, lim = lpHelperVec.size(); ii < lim; ++ii)
      cout << " > lb_cost for lp \"" << lpHelperVec[ii].name << "\" " << lpHelperVec[ii].lb_cost << endl;
  }

  // the cv-based load balancer requires 2 data values in each cv: the global index and the wgt

  int8 (*cv_global_wgt)[2] = new int8[sm->ncv][2];
  for (int icv = 0; icv < sm->ncv; ++icv) {
    cv_global_wgt[icv][0] = sm->cvora[mpi_rank] + icv; // global cv index -- need this
    cv_global_wgt[icv][1] = cv_lb_cost;
  }

  // if we have particles and the user supplied a positive particle cost, augment the cv cost...
  for (int ii_sm = 0, lim = sm->lpHelperVec.size(); ii_sm < lim; ++ii_sm) {
    // find lp helper...
    map<const string,int>::iterator iter = lpHelperNameMap.find(sm->lpHelperVec[ii_sm].name);
    assert(iter != lpHelperNameMap.end());
    const int ii = iter->second;
    for (int icv = 0; icv < sm->ncv; ++icv)
      cv_global_wgt[icv][1] += lpHelperVec[ii].lb_cost*(sm->lpHelperVec[ii_sm].lpocv_i_global[icv+1]-sm->lpHelperVec[ii_sm].lpocv_i_global[icv]);
  }

  // if we have solution based load balanced it will be in sm->cv_lb_wgt (defined by solvers)...
  if (sm->cv_lb_wgt) {
    for (int icv = 0; icv < sm->ncv; ++icv) {
      cv_global_wgt[icv][1] += sm->cv_lb_wgt[icv];
    }
  }

  // loop through faces and sum face-based part of cv cost...

  int * send_count = new int[mpi_size];
  FOR_RANK send_count[rank] = 0;

  // primary faces...

  for (int ifa = 0; ifa < sm->nfa; ++ifa) {
    FOR_I2 {
      const int8 icv = sm->cvofa_global[ifa][i]&MASK_55BITS;
      assert((icv >= 0)&&(icv < sm->ncv_global));
      const int rank = MiscUtils::getRankInXora(icv,sm->cvora);
      if (rank != mpi_rank)
        send_count[rank] += 2;
    }
  }

  // extended faces...

  for (int ief = 0; ief < sm->nef; ++ief) {
    FOR_I2 {
      const int8 icv = sm->cvoef_global[ief][i]&MASK_52BITS;
      assert((icv >= 0)&&(icv < sm->ncv_global));
      const int rank = MiscUtils::getRankInXora(icv,sm->cvora);
      if (rank != mpi_rank)
        send_count[rank] += 2;
    }
  }

  // boundary faces...

  for (int ibf = 0; ibf < sm->nbf; ++ibf) {
    const int8 icv = sm->cvobf_global[ibf];
    assert((icv >= 0)&&(icv < sm->ncv_global));
    const int rank = MiscUtils::getRankInXora(icv,sm->cvora);
    if (rank != mpi_rank)
      send_count[rank] += 2;
  }

  int * send_disp = new int[mpi_size];
  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
  const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

  int8 * send_buf_int8 = new int8[send_count_sum];

  for (int ifa = 0; ifa < sm->nfa; ++ifa) {
    FOR_I2 {
      const int8 icv = sm->cvofa_global[ifa][i]&MASK_55BITS;
      assert((icv >= 0)&&(icv < sm->ncv_global));
      const int rank = MiscUtils::getRankInXora(icv,sm->cvora);
      if (rank == mpi_rank) {
        const int icv_local = icv - sm->cvora[mpi_rank];
        cv_global_wgt[icv_local][1] += fa_lb_cost;
      }
      else {
        send_buf_int8[send_disp[rank]  ] = icv;
        send_buf_int8[send_disp[rank]+1] = fa_lb_cost;
        send_disp[rank] += 2;
      }
    }
  }

  // extended faces...

  for (int ief = 0; ief < sm->nef; ++ief) {
    FOR_I2 {
      const int8 icv = sm->cvoef_global[ief][i]&MASK_52BITS;
      assert((icv >= 0)&&(icv < sm->ncv_global));
      const int rank = MiscUtils::getRankInXora(icv,sm->cvora);
      if (rank == mpi_rank) {
        const int icv_local = icv - sm->cvora[mpi_rank];
        cv_global_wgt[icv_local][1] += ef_lb_cost;
      }
      else {
        send_buf_int8[send_disp[rank]  ] = icv;
        send_buf_int8[send_disp[rank]+1] = ef_lb_cost;
        send_disp[rank] += 2;
      }
    }
  }

  // boundary faces...

  for (int ibf = 0; ibf < sm->nbf; ++ibf) {
    const int izone = sm->zone_bf[ibf];
    assert((izone >= 0)&&(izone < int(bfZoneVec.size())));
    const int8 bf_lb_cost = bfZoneVec[izone].lb_cost;
    const int8 icv = sm->cvobf_global[ibf];
    assert((icv >= 0)&&(icv < sm->ncv_global));
    const int rank = MiscUtils::getRankInXora(icv,sm->cvora);
    if (rank == mpi_rank) {
      const int icv_local = icv - sm->cvora[mpi_rank];
      cv_global_wgt[icv_local][1] += bf_lb_cost;
    }
    else {
      send_buf_int8[send_disp[rank]  ] = icv;
      send_buf_int8[send_disp[rank]+1] = bf_lb_cost;
      send_disp[rank] += 2;
    }
  }

  // rewind...

  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

  // exchange...

  int * recv_count = new int[mpi_size];
  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

  int * recv_disp = new int[mpi_size];
  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
  const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

  int8 * recv_buf_int8 = new int8[recv_count_sum];
  MPI_Alltoallv(send_buf_int8,send_count,send_disp,MPI_INT8,
      recv_buf_int8,recv_count,recv_disp,MPI_INT8,mpi_comm);
  delete[] send_buf_int8;
  delete[] send_count;
  delete[] send_disp;
  delete[] recv_count;
  delete[] recv_disp;

  // now unpack into cost...

  for (int irecv = 0; irecv < recv_count_sum; irecv += 2) {
    const int8 icv = recv_buf_int8[irecv];
    assert((icv >= sm->cvora[mpi_rank])&&(icv < sm->cvora[mpi_rank+1]));
    const int icv_local = icv - sm->cvora[mpi_rank];
    cv_global_wgt[icv_local][1] += recv_buf_int8[irecv+1];
  }
  delete[] recv_buf_int8;

  /*
     {
     int8 my_check_buf[2] = { sm->ncv, 0 };
     for (int icv = 0; icv < sm->ncv; ++icv)
     my_check_buf[1] += cv_global_wgt[icv][1];
     int8 buf_sum[2];
     MPI_Reduce(my_check_buf,buf_sum,2,MPI_INT8,MPI_SUM,0,mpi_comm);
     if (mpi_rank == 0)
     cout << "CHECK ncv: " << buf_sum[0] << " total cost: " << buf_sum[1] << endl;
     }
     */

  ncv = sm->ncv;
  repartXcvWeightedPadt(sm->x_vv,cv_global_wgt,ncv,mpi_comm); // NOTE: changes sm->x_vv,cv_global_wgt,ncv !

  // check final load balance...
  int8 my_buf[2] = { ncv, 0 };
  for (int icv = 0; icv < ncv; ++icv)
    my_buf[1] += cv_global_wgt[icv][1];

  int8 buf_min[2];
  MPI_Reduce(my_buf,buf_min,2,MPI_INT8,MPI_MIN,0,mpi_comm);
  int8 buf_max[2];
  MPI_Reduce(my_buf,buf_max,2,MPI_INT8,MPI_MAX,0,mpi_comm);
  int8 buf_sum[2];
  MPI_Reduce(my_buf,buf_sum,2,MPI_INT8,MPI_SUM,0,mpi_comm);

  if (mpi_rank == 0) {
    const double avg_cost = double(buf_sum[1])/double(mpi_size) ;
    cout << " > load balance: ncv min:avg:max: " <<
      buf_min[0] << ":" << double(buf_sum[0])/double(mpi_size) << ":" << buf_max[0] <<
      " cost ratio min:avg:max: " << double(buf_min[1]/avg_cost) << ":1.0:" << double(buf_max[1])/avg_cost << endl;
  }


  // we need one more reordering step... to delineate internal vs proc boundary cells
  // ie a cell is internal if it has no cell neighbors who are ghosts ...

  {

    // set a initial ordering based on the locations of the cells...

    assert( icv_global == NULL);
    icv_global = new int8[ncv];
    for (int icv = 0; icv < ncv; ++icv)
      icv_global[icv] = cv_global_wgt[icv][0];
    delete[] cv_global_wgt;
    reorderXcv(sm->x_vv,icv_global,ncv);

    DistributedDataExchanger*  dde_cv_tmp = new DistributedDataExchanger(icv_global,ncv,sm->cvora);

    uint8* rbi = new uint8[ncv];
    for (int icv = 0; icv < ncv; ++icv)
      rbi[icv] = BitUtils::packRankBitsIndex(mpi_rank,0,icv);

    assert( rbi_sm == NULL); rbi_sm = new uint8[sm->ncv];
    dde_cv_tmp->push(rbi_sm,rbi);
    delete[] rbi;
    delete dde_cv_tmp;

    // we need to pull the cvofa_global associated with this distribution...

    int * bits = new int[sm->nfa];

    for (int ifa =0; ifa < sm->nfa; ++ifa) {
      // the bits should only be in sm->cvofa_global[ifa][1]...
      bits[ifa] = int(sm->cvofa_global[ifa][1]>>55);
      if (bits[ifa]) {
        sm->cvofa_global[ifa][1] &= MASK_55BITS;
      }
    }

    uint8 (*rbiofa_sm)[2] = new uint8[sm->nfa][2];

    DistributedDataExchanger*  dde_fa = new DistributedDataExchanger((int8*)sm->cvofa_global,sm->nfa*2,sm->cvora);
    dde_fa->pull((uint8*)rbiofa_sm,rbi_sm);
    delete dde_fa;
    DELETE(rbi_sm);

    for (int ifa = 0; ifa < sm->nfa; ++ifa) {
      if (bits[ifa]) {
        sm->cvofa_global[ifa][1] |= (int8(bits[ifa])<<55);
        int rank1,bits1,index1;
        BitUtils::unpackRankBitsIndex(rank1,bits1,index1,rbiofa_sm[ifa][1]);
        assert(bits1 == 0); // no periodicity for now
        rbiofa_sm[ifa][1] = BitUtils::packRankBitsIndex(rank1,bits[ifa],index1);
      }
    }

    delete[] bits;

    int * send_count   = new int[mpi_size];
    int * send_disp    = new int[mpi_size];
    int * send_buf_int = NULL;

    FOR_RANK send_count[rank] = 0;

    for (int iter = 0; iter < 2; ++iter) {

      for (int ifa =0; ifa < sm->nfa; ++ifa) {

        int rank0,bits0,index0;
        BitUtils::unpackRankBitsIndex(rank0,bits0,index0,rbiofa_sm[ifa][0]);
        assert(bits0 == 0);
        int rank1,bits1,index1;
        BitUtils::unpackRankBitsIndex(rank1,bits1,index1,rbiofa_sm[ifa][1]);

        if ( bits1 != 0) {

          if ( iter == 0) {
            ++send_count[rank0];
          } else {
            send_buf_int[send_disp[rank0]++] = index0;
          }

        } else if ( rank0 != rank1) {

          if ( iter == 0) {
            ++send_count[rank0];
          } else {
            send_buf_int[send_disp[rank0]++] = index0;
          }

          if ( iter == 0) {
            ++send_count[rank1];
          } else {
            send_buf_int[send_disp[rank1]++] = index1;
          }
        }
      }

      // rebuidl the disp on both iterations..
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      if ( iter == 0 )
        send_buf_int = new int[send_disp[mpi_size-1]+send_count[mpi_size-1]];

    }//iter

    delete[] rbiofa_sm;

    int * recv_count = new int[mpi_size];
    int * recv_disp  = new int[mpi_size];

    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

    const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    int * recv_buf_int = new int[recv_disp[mpi_size-1] + recv_count[mpi_size-1]];

    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
        recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);

    delete[] send_buf_int;
    delete[] send_disp;
    delete[] send_count;
    delete[] recv_disp;
    delete[] recv_count;

    int * cv_flag = new int[ncv];
    for (int icv = 0; icv < ncv; ++icv)
      cv_flag[icv] = 0;

    for (int ii = 0; ii < recv_count_sum; ++ii) {
      cv_flag[recv_buf_int[ii]]--;
    }

    delete[] recv_buf_int;

    const int ncv_old = ncv;
    ncv               = 0;
    for (int iter = 0; iter < 2; ++iter) {

      for (int icv = 0; icv < ncv_old; ++icv) {

        if ( iter == 0 ) { // number the internal cells...

          if ( cv_flag[icv] == 0) {
            cv_flag[icv] = ncv++;
          }

        } else {

          if ( cv_flag[icv] < 0) { // number the cells with ghost cv nbrs..
            cv_flag[icv] = ncv++;
          }
        }
      }

      if ( iter == 0 )
        ncv_i = ncv;
    }

    assert( ncv == ncv_old);

    // reorder icv_global and x_vv...
    int8* icv_global_new = new int8[ncv];
    double (*x_vv_new)[3] = new double[ncv][3];

    for (int icv = 0; icv < ncv; ++icv) {
      icv_global_new[cv_flag[icv]] = icv_global[icv];
      for (int i =0; i < 3; ++i)
        x_vv_new[cv_flag[icv]][i] = sm->x_vv[icv][i];
    }

    delete[] cv_flag;
    delete[] icv_global; icv_global = icv_global_new;
    delete[] sm->x_vv;    sm->x_vv  = x_vv_new;

  }


  // and reorder on the rank...

  /*
     assert(icv_global == NULL);
     icv_global = new int8[ncv];
     for (int icv = 0; icv < ncv; ++icv)
     icv_global[icv] = cv_global_wgt[icv][0];
     delete[] cv_global_wgt;
     reorderXcv(sm->x_vv,icv_global,ncv);
     */


  // take a look...

  if (checkParam("WRITE_LOAD_BALANCE")) {
    char filename[128];
    sprintf(filename,"load_balance.%06d.dat",mpi_rank);
    FILE * fp = fopen(filename,"w");
    FOR_ICV {
      fprintf(fp,"%18.16e %18.16e %18.16e\n",sm->x_vv[icv][0],sm->x_vv[icv][1],sm->x_vv[icv][2]);
    }
    fclose(fp);
  }

}

void StaticSolver::redistReorderMesh(StripedMesh * sm) {

  COUT1("StaticSolver::redistReorderMesh()");

  // must be called AFTER loadBalance...

  assert(icv_global);
  assert(sm->cvora);
  assert(sm->cvofa_global);

  // set global counts from sm...
  ncv_global = sm->ncv_global;
  nbf_global = sm->nbf_global;
  nfa_global = sm->nfa_global;
  nno_global = sm->nno_global;
  nef_global = sm->nef_global;
  nno_pb_global = sm->nno_pb_global;

  // copy of periodic transform vec from sm...
  //periodicTransformVec = sm->periodicTransformVec;

  // build the striped dde -- used elsewhere for cv reading/writing...

  assert(dde_cv_striped == NULL);
  dde_cv_striped = new DistributedDataExchanger(icv_global,ncv,sm->cvora);
  // also take a copy of the sm->cvora...
  assert(cvora_striped == NULL);
  cvora_striped = new int8[mpi_size+1];
  for (int rank = 0; rank <= mpi_size; ++rank) cvora_striped[rank] = sm->cvora[rank];

  // push the local rbi back to sm's striped distribution of cvs.

  uint8 * rbi = new uint8[ncv];
  FOR_ICV rbi[icv] = BitUtils::packRankBitsIndex(mpi_rank,0,icv);

  assert(rbi_sm == NULL); rbi_sm = new uint8[sm->ncv];
  dde_cv_striped->push(rbi_sm,rbi);
  delete[] rbi;

  // ===================================================================
  // boundary faces...
  // ===================================================================

  // the bf's need to know where they are going...

  uint8 *rbiobf_sm = new uint8[sm->nbf];
  {
    DistributedDataExchanger dde(sm->cvobf_global,sm->nbf,sm->cvora);
    dde.pull(rbiobf_sm,rbi_sm);
  }

  // now count...

  int * send_count = new int[mpi_size];
  FOR_RANK send_count[rank] = 0;

  int * send_count2 = new int[mpi_size]; // nodes
  FOR_RANK send_count2[rank] = 0;

  // at this point, we will be processing the surface stuff as well, if
  // present in the striped mesh...

  assert(subSurface == NULL);
  if (sm->b_surface) {
    subSurface = new SubSurface();
    subSurface->nst_global = sm->surface_nst_global;
    subSurface->nsp_global = sm->surface_nsp_global;
  }

  for (int ibf = 0; ibf < sm->nbf; ++ibf) {
    int rank,bits,index;
    BitUtils::unpackRankBitsIndex(rank,bits,index,rbiobf_sm[ibf]);
    assert(bits == 0);
    const int nnob = sm->noobf_i[ibf+1]-sm->noobf_i[ibf];
    ++send_count[rank];
    send_count2[rank] += 2 + nnob; // ibf_global + icv0_global + nodes
    // for the surface, use the int8 buf. this is ist | (bits<<52)...
    if (sm->b_surface) {
      const int nstob = sm->surface_sbobf_i[ibf+1]-sm->surface_sbobf_i[ibf];
      send_count2[rank] += nstob;
      if (sm->surface_spobf_i != NULL) {
        const int nspob = sm->surface_spobf_i[ibf+1]-sm->surface_spobf_i[ibf];
        send_count2[rank] += nspob*2; // include uint8 v and double wgts here
      }
    }
  }

  // allocate...

  int * send_disp = new int[mpi_size];
  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
  int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

  int * send_buf_int = new int[send_count_sum*5]; // add space for nstob and nspob, even when surface is not available
  double * send_buf_double = new double[send_count_sum*17]; // Gij_bf[3][3]

  int * send_disp2 = new int[mpi_size];
  send_disp2[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp2[rank] = send_count2[rank-1] + send_disp2[rank-1];
  int send_count2_sum = send_disp2[mpi_size-1] + send_count2[mpi_size-1];

  uint8 * send_buf_uint8 = new uint8[send_count2_sum]; // for global node indices

  // pack...

  for (int ibf = 0; ibf < sm->nbf; ++ibf) {
    int rank,bits,index;
    BitUtils::unpackRankBitsIndex(rank,bits,index,rbiobf_sm[ibf]);
    assert(bits == 0);
    const int nnob = sm->noobf_i[ibf+1]-sm->noobf_i[ibf];
    // send to rank...
    send_buf_int[send_disp[rank]*5  ] = sm->zone_bf[ibf];
    send_buf_int[send_disp[rank]*5+1] = index; // cvobf[ibf] on the final partition
    send_buf_int[send_disp[rank]*5+2] = nnob;
    send_buf_int[send_disp[rank]*5+3] = 0; // nstob: set below when b_surface == true
    send_buf_int[send_disp[rank]*5+4] = 0; // nspob: set below when b_surface == true
    // geometry...
    send_buf_double[send_disp[rank]*17   ]  = sm->n_bf[ibf][0];
    send_buf_double[send_disp[rank]*17+1 ]  = sm->n_bf[ibf][1];
    send_buf_double[send_disp[rank]*17+2 ]  = sm->n_bf[ibf][2];
    send_buf_double[send_disp[rank]*17+3 ]  = sm->x_bf[ibf][0];
    send_buf_double[send_disp[rank]*17+4 ]  = sm->x_bf[ibf][1];
    send_buf_double[send_disp[rank]*17+5 ]  = sm->x_bf[ibf][2];
    send_buf_double[send_disp[rank]*17+6 ]  = sm->area_bf[ibf];
    send_buf_double[send_disp[rank]*17+7 ]  = sm->area_over_delta_bf[ibf];
    send_buf_double[send_disp[rank]*17+8 ]  = sm->Gij_bf[ibf][0][0];
    send_buf_double[send_disp[rank]*17+9 ]  = sm->Gij_bf[ibf][0][1];
    send_buf_double[send_disp[rank]*17+10]  = sm->Gij_bf[ibf][0][2];
    send_buf_double[send_disp[rank]*17+11]  = sm->Gij_bf[ibf][1][0];
    send_buf_double[send_disp[rank]*17+12]  = sm->Gij_bf[ibf][1][1];
    send_buf_double[send_disp[rank]*17+13]  = sm->Gij_bf[ibf][1][2];
    send_buf_double[send_disp[rank]*17+14]  = sm->Gij_bf[ibf][2][0];
    send_buf_double[send_disp[rank]*17+15]  = sm->Gij_bf[ibf][2][1];
    send_buf_double[send_disp[rank]*17+16]  = sm->Gij_bf[ibf][2][2];
    // int8 stuff...
    send_buf_uint8[send_disp2[rank]++] = sm->bfora[mpi_rank] + ibf; // ibf_global
    send_buf_uint8[send_disp2[rank]++] = sm->cvobf_global[ibf]; // icv_global for checking
    // nodes: just send the global node number...
    for (int nob = sm->noobf_i[ibf]; nob != sm->noobf_i[ibf+1]; ++nob)
      send_buf_uint8[send_disp2[rank]++] = sm->noobf_v_global[nob];
    // surface stuff...
    if (sm->b_surface) {
      // stobf...
      const int nstob = sm->surface_sbobf_i[ibf+1]-sm->surface_sbobf_i[ibf];
      send_buf_int[send_disp[rank]*5+3] = nstob;
      for (int stob = sm->surface_sbobf_i[ibf]; stob != sm->surface_sbobf_i[ibf+1]; ++stob) {
        send_buf_uint8[send_disp2[rank]++] = sm->surface_sbobf_v_global[stob];
      }
      // spobf...
      if (sm->surface_spobf_i != NULL) {
        const int nspob = sm->surface_spobf_i[ibf+1]-sm->surface_spobf_i[ibf];
        send_buf_int[send_disp[rank]*5+4] = nspob;
        assert(sizeof(uint8) == sizeof(double));
        memcpy(send_buf_uint8+send_disp2[rank],sm->surface_spobf_v_global+sm->surface_spobf_i[ibf],nspob*sizeof(uint8));
        send_disp2[rank] += nspob;
        memcpy(send_buf_uint8+send_disp2[rank],sm->surface_spobf_wgt_global+sm->surface_spobf_i[ibf],nspob*sizeof(double));
        send_disp2[rank] += nspob;
      }
    }
    ++send_disp[rank];
  }

  delete[] rbiobf_sm;

  // rewind disp...

  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

  // exchange...

  int * recv_count = new int[mpi_size];
  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

  int * recv_disp = new int[mpi_size];
  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
  int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

  FOR_RANK {
    send_count[rank] *= 5;
    send_disp[rank] *= 5;
    recv_count[rank] *= 5;
    recv_disp[rank] *= 5;
  }

  int * recv_buf_int = new int[recv_count_sum*5];
  MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
      recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
  delete[] send_buf_int; send_buf_int = NULL;

  FOR_RANK {
    send_count[rank] = (send_count[rank]/5)*17;
    send_disp[rank] = (send_disp[rank]/5)*17;
    recv_count[rank] = (recv_count[rank]/5)*17;
    recv_disp[rank] = (recv_disp[rank]/5)*17;
  }

  double * recv_buf_double = new double[recv_count_sum*17];
  MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
      recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
  delete[] send_buf_double; send_buf_double = NULL;

  // and the int8 stuff...

  // rewind...

  send_disp2[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp2[rank] = send_count2[rank-1] + send_disp2[rank-1];

  int * recv_count2 = new int[mpi_size];
  MPI_Alltoall(send_count2,1,MPI_INT,recv_count2,1,MPI_INT,mpi_comm);

  int * recv_disp2 = new int[mpi_size];
  recv_disp2[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp2[rank] = recv_count2[rank-1] + recv_disp2[rank-1];
  int recv_count2_sum = recv_disp2[mpi_size-1] + recv_count2[mpi_size-1];

  uint8 * recv_buf_uint8 = new uint8[recv_count2_sum];
  MPI_Alltoallv(send_buf_uint8,send_count2,send_disp2,MPI_UINT8,
      recv_buf_uint8,recv_count2,recv_disp2,MPI_UINT8,mpi_comm);
  delete[] send_buf_uint8; send_buf_uint8 = NULL;

  // unpack and sort on the recv side...
  // bf's are grouped as follows:
  //

  // these are designed to produce a set of sorted bfs:
  // 1. zone
  // 2. icv (local)
  // 3. ibf_global (through irecv -- tricky!)

  vector< pair<pair<int,int>,pair<int,int> > > bfVec;  // first.first = zone, first.second = icv (local, for sort), second.first = irecv, second.second = irecv2

  int noobf_s = 0;
  int sbobf_s = 0;
  int spobf_s = 0;
  int irecv2 = 0;
  for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
    // grab just enough data to figure out what kind of bf this is...
    const int izone = recv_buf_int[irecv*5  ];
    const int icv   = recv_buf_int[irecv*5+1]; assert((icv >= 0)&&(icv < ncv));
    const int nnob  = recv_buf_int[irecv*5+2];
    const int nstob = recv_buf_int[irecv*5+3];
    const int nspob = recv_buf_int[irecv*5+4];
    bfVec.push_back( pair<pair<int,int>,pair<int,int> >(pair<int,int>(izone,icv),pair<int,int>(irecv,irecv2)) );
    noobf_s += nnob;
    sbobf_s += nstob;
    spobf_s += nspob;
    irecv2  += 2 + nnob + nstob + nspob*2;
  }
  assert(irecv2 == recv_count2_sum);

  // sort...

  sort(bfVec.begin(),bfVec.end());

  // and set bf data...

  nbf = bfVec.size();

  assert(ibf_global == NULL); ibf_global = new int8[nbf];
  assert(zone_bf == NULL); zone_bf = new int[nbf];
  assert(cvobf == NULL); cvobf = new int[nbf];
  assert(n_bf == NULL); n_bf = new double[nbf][3];
  assert(x_bf == NULL); x_bf = new double[nbf][3];
  assert(area_bf == NULL); area_bf = new double[nbf];
  assert(area_over_delta_bf == NULL); area_over_delta_bf = new double[nbf];
  assert(Gij_bf == NULL); Gij_bf = new double[nbf][3][3];
  assert(noobf_i == NULL); noobf_i = new int[nbf+1]; noobf_i[0] = 0;
  int8 * noobf_v_global = new int8[noobf_s];

  assert(sstobf_i == NULL);
  uint8 * sbobf_v_global = NULL;
  if (sm->b_surface) {
    sstobf_i = new int[nbf+1]; sstobf_i[0] = 0;
    sbobf_v_global = new uint8[sbobf_s];
  }

  assert(sspobf_i == NULL);
  uint8 * spobf_v_global = NULL;
  double * spobf_wgt_global = NULL;
  // we need to support restarts that only read in the surface tris. This
  // was available in dev in late 2018 into 2019 to support lsp boundary
  // conditions...
  if ((sm->b_surface)&&(sm->surface_spobf_i != NULL)) {
    sspobf_i = new int[nbf+1]; sspobf_i[0] = 0;
    spobf_v_global = new uint8[spobf_s];
    spobf_wgt_global = new double[spobf_s];
  }
  else {
    assert(spobf_s == 0);
  }

  // unpack sorted data...

  int izone_current = -1;
  for (int ii = 0, nii = bfVec.size(); ii < nii; ++ii) {
    const int irecv = bfVec[ii].second.first;
    const int irecv2 = bfVec[ii].second.second;
    // the igr/icv...
    const int izone = bfVec[ii].first.first;
    if (izone != izone_current) {
      if (izone_current != -1) {
        bfZoneVec[izone_current].ibf_l = ii-1;
        bfZoneVec[izone_current].nbf = bfZoneVec[izone_current].ibf_l - bfZoneVec[izone_current].ibf_f + 1;
      }
      izone_current = izone;
      bfZoneVec[izone_current].ibf_f = ii;
    }
    const int icv = bfVec[ii].first.second;
    // unpack...
    assert(izone  == recv_buf_int[irecv*5  ]);
    assert(icv    == recv_buf_int[irecv*5+1]); assert((icv >= 0)&&(icv < ncv));
    const int nnob = recv_buf_int[irecv*5+2];
    const int nstob = recv_buf_int[irecv*5+3];
    const int nspob = recv_buf_int[irecv*5+4];
    ibf_global[ii] = recv_buf_uint8[irecv2];
    zone_bf[ii] = izone;
    cvobf[ii] = icv;
    assert(icv_global[icv] == recv_buf_uint8[irecv2+1]); // good check
    n_bf[ii][0] = recv_buf_double[irecv*17  ];
    n_bf[ii][1] = recv_buf_double[irecv*17+1];
    n_bf[ii][2] = recv_buf_double[irecv*17+2];
    x_bf[ii][0] = recv_buf_double[irecv*17+3];
    x_bf[ii][1] = recv_buf_double[irecv*17+4];
    x_bf[ii][2] = recv_buf_double[irecv*17+5];
    area_bf[ii] = recv_buf_double[irecv*17+6];
    area_over_delta_bf[ii] = recv_buf_double[irecv*17+7];
    Gij_bf[ii][0][0] = recv_buf_double[irecv*17+8 ];
    Gij_bf[ii][0][1] = recv_buf_double[irecv*17+9 ];
    Gij_bf[ii][0][2] = recv_buf_double[irecv*17+10];
    Gij_bf[ii][1][0] = recv_buf_double[irecv*17+11];
    Gij_bf[ii][1][1] = recv_buf_double[irecv*17+12];
    Gij_bf[ii][1][2] = recv_buf_double[irecv*17+13];
    Gij_bf[ii][2][0] = recv_buf_double[irecv*17+14];
    Gij_bf[ii][2][1] = recv_buf_double[irecv*17+15];
    Gij_bf[ii][2][2] = recv_buf_double[irecv*17+16];
    // nodes...
    noobf_i[ii+1] = noobf_i[ii] + nnob;
    for (int nob = 0; nob < nnob; ++nob)
      noobf_v_global[noobf_i[ii]+nob] = recv_buf_uint8[irecv2+2+nob];
    // surface tris and pts if b_surface...
    if (sm->b_surface) {
      // st stuff...
      sstobf_i[ii+1] = sstobf_i[ii] + nstob;
      for (int stob = 0; stob < nstob; ++stob)
        sbobf_v_global[sstobf_i[ii]+stob] = recv_buf_uint8[irecv2+2+nnob+stob];
      // and sp stuff: use memcpy here because we are pushing doubles through the uint8 buf...
      if (sm->surface_spobf_i != NULL) {
        sspobf_i[ii+1] = sspobf_i[ii] + nspob;
        memcpy(spobf_v_global+sspobf_i[ii],recv_buf_uint8+irecv2+2+nnob+nstob,nspob*sizeof(uint8));
        memcpy(spobf_wgt_global+sspobf_i[ii],recv_buf_uint8+irecv2+2+nnob+nstob+nspob,nspob*sizeof(double));
      }
      else {
        assert(nspob == 0);
      }
    }
  }
  if (izone_current != -1) {
    bfZoneVec[izone_current].ibf_l = bfVec.size()-1;
    bfZoneVec[izone_current].nbf = bfZoneVec[izone_current].ibf_l - bfZoneVec[izone_current].ibf_f + 1;
  }

  // check ibf indexing within zones...
  int8 nbf_check = 0; // use int8 here for reduction -- not necessary locally
  FOR_IZONE(bfZoneVec) {
    assert( bfZoneVec[izone].ibf_l - bfZoneVec[izone].ibf_f + 1 ==  bfZoneVec[izone].nbf );
    nbf_check += bfZoneVec[izone].nbf;
    for (int ibf = bfZoneVec[izone].ibf_f; ibf <= bfZoneVec[izone].ibf_l; ++ibf) {
      assert(zone_bf[ibf] == izone);
      // also, the global indices for this bf should be in the
      // global range set in the bfZoneVec...
      assert((ibf_global[ibf] >= bfZoneVec[izone].ibf_f_global)&&(ibf_global[ibf] < bfZoneVec[izone].ibf_f_global+bfZoneVec[izone].nbf_global));
    }
    // set the bfZone's ptrs into the solver nbf data with zero indexing...
    assert(bfZoneVec[izone].ibf_global == NULL);         bfZoneVec[izone].ibf_global         = ibf_global         + bfZoneVec[izone].ibf_f;
    assert(bfZoneVec[izone].area_bf == NULL);            bfZoneVec[izone].area_bf            = area_bf            + bfZoneVec[izone].ibf_f;
    assert(bfZoneVec[izone].area_over_delta_bf == NULL); bfZoneVec[izone].area_over_delta_bf = area_over_delta_bf + bfZoneVec[izone].ibf_f;
    assert(bfZoneVec[izone].n_bf == NULL);               bfZoneVec[izone].n_bf               = n_bf               + bfZoneVec[izone].ibf_f;
    assert(bfZoneVec[izone].x_bf == NULL);               bfZoneVec[izone].x_bf               = x_bf               + bfZoneVec[izone].ibf_f;
    assert(bfZoneVec[izone].cvobf == NULL);              bfZoneVec[izone].cvobf              = cvobf              + bfZoneVec[izone].ibf_f;
  }
  assert(nbf_check == nbf);

  int8 nbf_global_check;
  MPI_Reduce(&nbf_check,&nbf_global_check,1,MPI_INT8,MPI_SUM,0,mpi_comm);
  if (mpi_rank == 0) {
    assert(nbf_global == nbf_global_check);
  }

  bfVec.clear();

  delete[] recv_buf_int;    recv_buf_int = NULL;
  delete[] recv_buf_double; recv_buf_double = NULL;
  delete[] recv_buf_uint8;  recv_buf_uint8 = NULL;

  // ===================================================================
  // compact faces...
  // ===================================================================

  // the faces need to know where they are going...

  uint8 (*rbiofa_sm)[2] = new uint8[sm->nfa][2];
  {

    // for periodic, need to copy this guy and remove bits, then reapply them...
    int * bits = new int[sm->nfa];

    for (int ifa = 0; ifa < sm->nfa; ++ifa) {
      // the bits should only be in sm->cvofa_global[ifa][1]...
      bits[ifa] = int(sm->cvofa_global[ifa][1]>>55);
      if (bits[ifa]) {
        sm->cvofa_global[ifa][1] &= MASK_55BITS;
      }
    }

    DistributedDataExchanger dde((int8*)sm->cvofa_global,sm->nfa*2,sm->cvora); // will fail when periodic bits are in cvofa_global
    dde.pull((uint8*)rbiofa_sm,rbi_sm);

    for (int ifa = 0; ifa < sm->nfa; ++ifa) {
      if (bits[ifa]) {
        sm->cvofa_global[ifa][1] |= (int8(bits[ifa])<<55);
        int rank1,bits1,index1;
        BitUtils::unpackRankBitsIndex(rank1,bits1,index1,rbiofa_sm[ifa][1]);
        assert(bits1 == 0); // no periodicity for now
        rbiofa_sm[ifa][1] = BitUtils::packRankBitsIndex(rank1,bits[ifa],index1);
      }
    }

    delete[] bits;

  }

  // now count...

  FOR_RANK send_count[rank] = 0;

  FOR_RANK send_count2[rank] = 0;

  for (int ifa = 0; ifa < sm->nfa; ++ifa) {
    int rank0,bits0,index0;
    BitUtils::unpackRankBitsIndex(rank0,bits0,index0,rbiofa_sm[ifa][0]);
    assert(bits0 == 0);
    int rank1,bits1,index1;
    BitUtils::unpackRankBitsIndex(rank1,bits1,index1,rbiofa_sm[ifa][1]);
    const int nnof = sm->noofa_i[ifa+1]-sm->noofa_i[ifa];
    if ((rank0 == rank1)||(bits1 != 0)) { // if the bits are set, then only send once...
      // send once to rank0...
      ++send_count[rank0];
      send_count2[rank0] += 2 + nnof;
    }
    else {
      ++send_count[rank0];
      ++send_count[rank1];
      send_count2[rank0] += 2 + nnof;
      send_count2[rank1] += 2 + nnof;
    }

  }

  // allocate...

  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
  send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

  assert(send_buf_int == NULL); send_buf_int = new int[send_count_sum*5];
  assert(send_buf_double == NULL); send_buf_double = new double[send_count_sum*6];

  send_disp2[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp2[rank] = send_count2[rank-1] + send_disp2[rank-1];
  send_count2_sum = send_disp2[mpi_size-1] + send_count2[mpi_size-1];

  assert(send_buf_uint8 == NULL);
  int8 * send_buf_int8 = new int8[send_count2_sum]; // for global node indices

  // pack...

  for (int ifa = 0; ifa < sm->nfa; ++ifa) {
    int rank0,bits0,index0;
    BitUtils::unpackRankBitsIndex(rank0,bits0,index0,rbiofa_sm[ifa][0]);
    assert(bits0 == 0);
    int rank1,bits1,index1;
    BitUtils::unpackRankBitsIndex(rank1,bits1,index1,rbiofa_sm[ifa][1]);
    const int nnof = sm->noofa_i[ifa+1]-sm->noofa_i[ifa];
    if ((rank0 == rank1)||(bits1 != 0)) { // in this case, and bits present means the face goes to just one rank...
      // send once to rank0==rank1...
      send_buf_int[send_disp[rank0]*5  ] = index0; // cvofa[ifa][0] on the final partition
      send_buf_int[send_disp[rank0]*5+1] = rank1;  // nbr info
      send_buf_int[send_disp[rank0]*5+2] = bits1;
      send_buf_int[send_disp[rank0]*5+3] = index1;
      send_buf_int[send_disp[rank0]*5+4] = nnof;
      // geometry...
      send_buf_double[send_disp[rank0]*6  ]  = sm->n_fa[ifa][0];
      send_buf_double[send_disp[rank0]*6+1]  = sm->n_fa[ifa][1];
      send_buf_double[send_disp[rank0]*6+2]  = sm->n_fa[ifa][2];
      send_buf_double[send_disp[rank0]*6+3]  = sm->x_fa[ifa][0];
      send_buf_double[send_disp[rank0]*6+4]  = sm->x_fa[ifa][1];
      send_buf_double[send_disp[rank0]*6+5]  = sm->x_fa[ifa][2];
      ++send_disp[rank0];
      // int8 stuff...
      send_buf_int8[send_disp2[rank0]++] = sm->faora[mpi_rank] + ifa; // ifa_global
      send_buf_int8[send_disp2[rank0]++] = sm->cvofa_global[ifa][0]; // icv0_global for checking
      // nodes: just send the global node number...
      for (int nof = sm->noofa_i[ifa]; nof != sm->noofa_i[ifa+1]; ++nof)
        send_buf_int8[send_disp2[rank0]++] = sm->noofa_v_global[nof];
    }
    else {
      // send to rank0...
      send_buf_int[send_disp[rank0]*5  ] = index0; // cvofa[ifa][0] on the final partition
      send_buf_int[send_disp[rank0]*5+1] = rank1;  // nbr info
      send_buf_int[send_disp[rank0]*5+2] = bits1; assert(bits1 == 0);
      send_buf_int[send_disp[rank0]*5+3] = index1;
      send_buf_int[send_disp[rank0]*5+4] = nnof;
      // geometry...
      send_buf_double[send_disp[rank0]*6  ]  = sm->n_fa[ifa][0];
      send_buf_double[send_disp[rank0]*6+1]  = sm->n_fa[ifa][1];
      send_buf_double[send_disp[rank0]*6+2]  = sm->n_fa[ifa][2];
      send_buf_double[send_disp[rank0]*6+3]  = sm->x_fa[ifa][0];
      send_buf_double[send_disp[rank0]*6+4]  = sm->x_fa[ifa][1];
      send_buf_double[send_disp[rank0]*6+5]  = sm->x_fa[ifa][2];
      ++send_disp[rank0];
      // int8 stuff...
      send_buf_int8[send_disp2[rank0]++] = sm->faora[mpi_rank] + ifa; // ifa_global
      send_buf_int8[send_disp2[rank0]++] = sm->cvofa_global[ifa][0]; // icv0_global for checking
      // nodes: just send the global node number...
      for (int nof = sm->noofa_i[ifa]; nof != sm->noofa_i[ifa+1]; ++nof)
        send_buf_int8[send_disp2[rank0]++] = sm->noofa_v_global[nof];
      // and rank1 - flip...
      send_buf_int[send_disp[rank1]*5  ] = index1; // cvofa[ifa][0] on the final partition
      send_buf_int[send_disp[rank1]*5+1] = rank0;  // nbr info
      send_buf_int[send_disp[rank1]*5+2] = bits0;
      send_buf_int[send_disp[rank1]*5+3] = index0;
      send_buf_int[send_disp[rank1]*5+4] = nnof;
      // geometry...
      send_buf_double[send_disp[rank1]*6  ]  = -sm->n_fa[ifa][0];
      send_buf_double[send_disp[rank1]*6+1]  = -sm->n_fa[ifa][1];
      send_buf_double[send_disp[rank1]*6+2]  = -sm->n_fa[ifa][2];
      send_buf_double[send_disp[rank1]*6+3]  = sm->x_fa[ifa][0];
      send_buf_double[send_disp[rank1]*6+4]  = sm->x_fa[ifa][1];
      send_buf_double[send_disp[rank1]*6+5]  = sm->x_fa[ifa][2];
      ++send_disp[rank1];
      // int8 stuff...
      send_buf_int8[send_disp2[rank1]++] = -(sm->faora[mpi_rank] + ifa) - 1; // ifa_global -1 indexed
      send_buf_int8[send_disp2[rank1]++] = sm->cvofa_global[ifa][1]; // icv0_global for checking
      // nodes: just send the FLIPPED global node number...
      for (int nof = sm->noofa_i[ifa+1]-1; nof >= sm->noofa_i[ifa]; --nof)
        send_buf_int8[send_disp2[rank1]++] = sm->noofa_v_global[nof];
    }
  }

  delete[] rbiofa_sm;

  // rewind disp...

  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

  // exchange...

  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
  recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

  FOR_RANK {
    send_count[rank] *= 5;
    send_disp[rank] *= 5;
    recv_count[rank] *= 5;
    recv_disp[rank] *= 5;
  }

  assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum*5];
  MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
      recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
  delete[] send_buf_int; send_buf_int = NULL;

  FOR_RANK {
    send_count[rank] = (send_count[rank]/5)*6;
    send_disp[rank] = (send_disp[rank]/5)*6;
    recv_count[rank] = (recv_count[rank]/5)*6;
    recv_disp[rank] = (recv_disp[rank]/5)*6;
  }

  assert(recv_buf_double == NULL); recv_buf_double = new double[recv_count_sum*6];
  MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
      recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
  delete[] send_buf_double; send_buf_double = NULL;

  // and the int8 stuff...

  // rewind...

  send_disp2[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp2[rank] = send_count2[rank-1] + send_disp2[rank-1];

  MPI_Alltoall(send_count2,1,MPI_INT,recv_count2,1,MPI_INT,mpi_comm);

  recv_disp2[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp2[rank] = recv_count2[rank-1] + recv_disp2[rank-1];
  recv_count2_sum = recv_disp2[mpi_size-1] + recv_count2[mpi_size-1];

  assert(recv_buf_uint8 == NULL);
  int8 * recv_buf_int8 = new int8[recv_count2_sum];
  MPI_Alltoallv(send_buf_int8,send_count2,send_disp2,MPI_INT8,
      recv_buf_int8,recv_count2,recv_disp2,MPI_INT8,mpi_comm);
  delete[] send_buf_int8; send_buf_int8 = NULL;

  // unpack and sort on the recv side...
  // faces are grouped as follows:
  //

  // these are designed to produce a set of sorted faces that are soted by:
  // 1. group
  // 2. icv (local)
  // 3. ifa_global (through irecv -- tricky!)

  vector< pair<int,pair<int,int> > > internalFaceVec; // first.first = icv (local, for sort), second.first = irecv, second.second = irecv2
  vector< pair<int,pair<int,int> > > interprocessorFaceVec;

  int noofa_s = 0;
  irecv2 = 0;
  for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
    // grab just enough data to figure out what kind of face this is...
    const int icv = recv_buf_int[irecv*5  ]; assert((icv >= 0)&&(icv < ncv));
    const int rank = recv_buf_int[irecv*5+1];
    const int bits = recv_buf_int[irecv*5+2]; assert((bits >= 0)&&(bits < (1<<6)));
    //const int icv_nbr = recv_buf_int[irecv*5+3];
    const int nnof = recv_buf_int[irecv*5+4];
    if ((rank == mpi_rank)&&(bits == 0)) {
      // internal face...
      internalFaceVec.push_back( pair<int,pair<int,int> >(icv,pair<int,int>(irecv,irecv2)) );
      noofa_s += nnof;
    }
    else {
      // interprocessor face...
      interprocessorFaceVec.push_back( pair<int,pair<int,int> >(icv,pair<int,int>(irecv,irecv2)) );
      noofa_s += nnof;
    }
    irecv2 += 2 + nnof;
  }
  assert(irecv2 == recv_count2_sum);

  // sort...

  sort(internalFaceVec.begin(),internalFaceVec.end());
  sort(interprocessorFaceVec.begin(),interprocessorFaceVec.end());

  // and set face data...

  nfa_i = internalFaceVec.size();
  nfa = nfa_i + interprocessorFaceVec.size();

  // allocate face data...

  assert(ifa_global == NULL); ifa_global = new int8[nfa];
  assert(cvofa == NULL); cvofa = new int[nfa][2];
  assert(n_fa == NULL); n_fa = new double[nfa][3];
  assert(x_fa == NULL); x_fa = new double[nfa][3];
  assert(noofa_i == NULL); noofa_i = new int[nfa+1]; noofa_i[0] = 0;

  // these ones built later -- no longer coming from striped mesh...
  assert(area_over_delta_fa == NULL);
  assert(group_fa == NULL);

  int8 * noofa_v_global = new int8[noofa_s];

  // unpack sorted data...

  // internal first...

  for (int ii = 0, nii = internalFaceVec.size(); ii < nii; ++ii) {
    const int irecv = internalFaceVec[ii].second.first;
    const int irecv2 = internalFaceVec[ii].second.second;
    // the icv...
    const int icv = internalFaceVec[ii].first;
    // unpack...
    assert(icv       == recv_buf_int[irecv*5  ]); assert((icv >= 0)&&(icv < ncv));
    const int rank    = recv_buf_int[irecv*5+1]; assert(rank == mpi_rank);
    const int bits    = recv_buf_int[irecv*5+2]; assert(bits == 0);
    const int icv_nbr = recv_buf_int[irecv*5+3]; assert((icv_nbr >= 0)&&(icv_nbr < ncv));
    const int nnof    = recv_buf_int[irecv*5+4];
    ifa_global[ii] = recv_buf_int8[irecv2];
    cvofa[ii][0] = icv;
    cvofa[ii][1] = icv_nbr;
    assert(icv_global[icv] == recv_buf_int8[irecv2+1]); // good check
    n_fa[ii][0] = recv_buf_double[irecv*6  ];
    n_fa[ii][1] = recv_buf_double[irecv*6+1];
    n_fa[ii][2] = recv_buf_double[irecv*6+2];
    x_fa[ii][0] = recv_buf_double[irecv*6+3];
    x_fa[ii][1] = recv_buf_double[irecv*6+4];
    x_fa[ii][2] = recv_buf_double[irecv*6+5];
    // finally nodes...
    noofa_i[ii+1] = noofa_i[ii] + nnof;
    for (int nof = 0; nof < nnof; ++nof)
      noofa_v_global[noofa_i[ii]+nof] = recv_buf_int8[irecv2+2+nof];
  }

  internalFaceVec.clear();

  map<const uint8,int> rbiMap;
  ncv_g = ncv;

  for (int ii = 0, nii = interprocessorFaceVec.size(); ii < nii; ++ii) {
    const int irecv = interprocessorFaceVec[ii].second.first;
    const int irecv2 = interprocessorFaceVec[ii].second.second;
    // the icv...
    const int icv = interprocessorFaceVec[ii].first;
    // unpack...
    assert(icv       == recv_buf_int[irecv*5  ]); assert((icv >= 0)&&(icv < ncv));
    const int rank    = recv_buf_int[irecv*5+1];
    const int bits    = recv_buf_int[irecv*5+2];
    const int icv_nbr = recv_buf_int[irecv*5+3];
    const int nnof    = recv_buf_int[irecv*5+4];
    ifa_global[ii+nfa_i] = recv_buf_int8[irecv2];
    cvofa[ii+nfa_i][0] = icv;
    const uint8 rbi_nbr = BitUtils::packRankBitsIndex(rank,bits,icv_nbr);
    map<const uint8,int>::iterator iter = rbiMap.find(rbi_nbr);
    if (iter == rbiMap.end()) {
      rbiMap[rbi_nbr] = cvofa[ii+nfa_i][1] = ncv_g++;
    }
    else {
      cvofa[ii+nfa_i][1] = iter->second;
    }
    assert(icv_global[icv] == recv_buf_int8[irecv2+1]); // good check
    n_fa[ii+nfa_i][0] = recv_buf_double[irecv*6  ];
    n_fa[ii+nfa_i][1] = recv_buf_double[irecv*6+1];
    n_fa[ii+nfa_i][2] = recv_buf_double[irecv*6+2];
    x_fa[ii+nfa_i][0] = recv_buf_double[irecv*6+3];
    x_fa[ii+nfa_i][1] = recv_buf_double[irecv*6+4];
    x_fa[ii+nfa_i][2] = recv_buf_double[irecv*6+5];
    // finally nodes...
    noofa_i[ii+nfa_i+1] = noofa_i[ii+nfa_i] + nnof;
    for (int nof = 0; nof < nnof; ++nof)
      noofa_v_global[noofa_i[ii+nfa_i]+nof] = recv_buf_int8[irecv2+2+nof];
  }

  assert(area_fa == NULL); area_fa = new double[nfa];
  FOR_IFA area_fa[ifa] = MAG(n_fa[ifa]);

  interprocessorFaceVec.clear();

  delete[] recv_buf_int;    recv_buf_int = NULL;
  delete[] recv_buf_double; recv_buf_double = NULL;
  delete[] recv_buf_int8;  recv_buf_int8 = NULL;

  // faces should be done!

  // ---------------------------------------------------------------
  // do nodes -- start with boundaries to list boundary nodes first
  // ---------------------------------------------------------------

  map<const int8,int> globalMap;
  nno = 0;

  assert(noobf_v == NULL); noobf_v = new int[noobf_s];
  for (int nob = 0; nob < noobf_s; ++nob) {
    const int8 ino_global = noobf_v_global[nob];
    assert((ino_global >= 0)&&(ino_global < nno_global));
    map<const int8,int>::iterator iter = globalMap.find(ino_global);
    if (iter == globalMap.end()) {
      globalMap[ino_global] = noobf_v[nob] = nno++;
    }
    else {
      noobf_v[nob] = iter->second;
    }
  }

  MiscUtils::dumpRange(noobf_v_global,noobf_s,"noobf_v_global");

  delete[] noobf_v_global;

  // remember the boundary node count...

  nno_b = nno;

  // then the faces...

  assert(noofa_v == NULL); noofa_v = new int[noofa_s];
  for (int nof = 0; nof < noofa_s; ++nof) {
    const int8 ino_global = noofa_v_global[nof];
    assert((ino_global >= 0)&&(ino_global < nno_global));
    map<const int8,int>::iterator iter = globalMap.find(ino_global);
    if (iter == globalMap.end()) {
      globalMap[ino_global] = noofa_v[nof] = nno++;
    }
    else {
      noofa_v[nof] = iter->second;
    }
  }

  MiscUtils::dumpRange(noofa_v_global,noofa_s,"noofa_v_global");

  delete[] noofa_v_global;

  // record the global node number -- not so important on Voronoi meshes,
  // but whatever...

  assert(nno_global == sm->noora[mpi_size]);
  assert(ino_global == NULL); ino_global = new int8[nno];
  for (map<const int8,int>::iterator iter = globalMap.begin(); iter != globalMap.end(); ++iter) {
    ino_global[iter->second] = iter->first;
    if ((iter->first < 0)||(iter->first >= nno_global)) {
      cout << "rank: " << mpi_rank << " BAD ino_global: " << iter->first << " nno_global: " << nno_global << endl;
    }
  }

  globalMap.clear();

  assert(x_no == NULL); x_no = new double[nno][3];
  {
    DistributedDataExchanger dde(ino_global,nno,sm->noora);
    dde.pull(x_no,sm->x_no);
  }

  // compare the normal and the node-based normal...

  double my_buf[4] = { 0.0, 0.0, 0.0, 0.0 };

  for (int ibf = 0; ibf < nbf; ++ibf) {
    double n_bf_check[3] = { 0.0, 0.0, 0.0 };
    const int ino = noobf_v[noobf_i[ibf]];
    int ino1 = noobf_v[noobf_i[ibf+1]-1];
    for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
      const int ino0 = ino1;
      ino1 = noobf_v[nob];
      if ((ino != ino0)&&(ino != ino1)) {
        const double this_n[3] = TRI_NORMAL_2(x_no[ino],x_no[ino0],x_no[ino1]);
        FOR_I3 n_bf_check[i] += 0.5*this_n[i];
      }
    }
    const double n_mag = MAG(n_bf[ibf]);
    const double n_mag_check = MAG(n_bf_check);
    //if ( (n_mag > 0.0) && (n_mag_check > 0.0))
    // XXX please revisit--- there appear to be faces that are failing the finite check
    //
    my_buf[0] = max( my_buf[0], fabs(	n_mag_check - n_mag ) );
    const double dp_check = DOT_PRODUCT(n_bf[ibf],n_bf_check)/(n_mag*n_mag_check);
    my_buf[1] = max( my_buf[1], fabs(	1.0 - dp_check )*n_mag );

  }

  for (int ifa = 0; ifa < nfa; ++ifa) {
    double n_fa_check[3] = { 0.0, 0.0, 0.0 };
    const int ino = noofa_v[noofa_i[ifa]];
    int ino1 = noofa_v[noofa_i[ifa+1]-1];
    for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      const int ino0 = ino1;
      ino1 = noofa_v[nof];
      if ((ino != ino0)&&(ino != ino1)) {
        const double this_n[3] = TRI_NORMAL_2(x_no[ino],x_no[ino0],x_no[ino1]);
        FOR_I3 n_fa_check[i] += 0.5*this_n[i];
      }
    }
    const double n_mag = MAG(n_fa[ifa]);
    const double n_mag_check = MAG(n_fa_check);
    my_buf[2] = max( my_buf[2], fabs(	n_mag_check - n_mag ) );
    const double dp_check = DOT_PRODUCT(n_fa[ifa],n_fa_check)/(n_mag*n_mag_check);
    my_buf[3] = max( my_buf[3], fabs(	1.0 - dp_check )*n_mag );
  }

  double buf[4];
  MPI_Reduce(my_buf,buf,4,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
  if (mpi_rank == 0) {
    cout << " > node bf normal check (should be small): " << buf[0] << " " << buf[1] << endl;
    cout << " > node fa normal check (should be small): " << buf[2] << " " << buf[3] << endl;
  }

  delete[] send_count2;
  delete[] send_disp2;
  delete[] recv_count2;
  delete[] recv_disp2;

  // =====================================
  // surface tris nst (and points nsp)...
  // =====================================

  assert(sstobf_v == NULL);
  assert(sspobf_v == NULL);
  assert(sspobf_wgt == NULL);
  if (sm->b_surface) {

    // TODO: would like to use uint8 here, but some of the
    // infrastructure (e.g. dde) does not support it...
    assert(globalMap.empty());

    //sm->dumpSurface();
    assert(subSurface->nst == 0);
    assert(subSurface->nsp == 0);

    assert(sbobf_v_global);
    sstobf_v = new int[sbobf_s]; // a subsurface tri reference

    for (int stob = 0; stob < sbobf_s; ++stob) {
      const int8 ist_bits_global = sbobf_v_global[stob];
      map<const int8,int>::iterator iter = globalMap.find(ist_bits_global);
      if (iter == globalMap.end()) {
        globalMap[ist_bits_global] = sstobf_v[stob] = subSurface->nst++;
      }
      else {
        sstobf_v[stob] = iter->second;
      }
    }

    delete[] sbobf_v_global; sbobf_v_global = NULL;

    // record the global surface tri index...

    assert(subSurface->ist_global_and_bits == NULL); // TODO: this should be private with accessors
    subSurface->ist_global_and_bits = new int8[subSurface->nst];
    for (map<const int8,int>::iterator iter = globalMap.begin(); iter != globalMap.end(); ++iter) {
      subSurface->ist_global_and_bits[iter->second] = iter->first;
    }

    globalMap.clear();

    // sp's...
    // protect this with the sm->surface_spobf_i indicator, because some mles files
    // written for particles did not write the spobf_i/v/wgt...

    if (sm->surface_spobf_i != NULL) {

      // this is a newer stitch output that includes spobf_i/v/wgt...

      assert(spobf_v_global);
      assert(spobf_wgt_global);
      sspobf_v = new int[spobf_s];
      sspobf_wgt = spobf_wgt_global; // don't delete spobf_wgt_global.

      for (int spob = 0; spob < spobf_s; ++spob) {
        const int8 isp_bits_global = spobf_v_global[spob];
        map<const int8,int>::iterator iter = globalMap.find(isp_bits_global);
        if (iter == globalMap.end()) {
          globalMap[isp_bits_global] = sspobf_v[spob] = subSurface->nsp++;
        }
        else {
          sspobf_v[spob] = iter->second;
        }
      }
      delete[] spobf_v_global; spobf_v_global = NULL;

      // leave globalMap and delay build of subSurface->isp_global_and_bits to below
      // because it is the same with or without spobf_i/v/wgt...

    }

    // now use the dde to bring over surface tri stuff from striped surface...
    int * bits = new int[max(subSurface->nst,subSurface->nsp)];

    // strip off any bits and use subSurface->ist_global_and_bits as an
    // index into sm->surface_stora...
    for (int ist = 0; ist < subSurface->nst; ++ist) {
      bits[ist] = int(subSurface->ist_global_and_bits[ist]>>52);
      if (bits[ist])
        subSurface->ist_global_and_bits[ist] &= MASK_52BITS;
    }

    DistributedDataExchanger * dde = new DistributedDataExchanger(subSurface->ist_global_and_bits,subSurface->nst,sm->surface_stora);

    for (int ist = 0; ist < subSurface->nst; ++ist) {
      if (bits[ist])
        subSurface->ist_global_and_bits[ist] |= (int8(bits[ist])<<52);
    }

    assert(subSurface->znost == NULL); subSurface->znost = new int[subSurface->nst];
    dde->pull(subSurface->znost,sm->surface_znost);

    // and the global node numbers...

    int (*spost_global)[3] = new int[subSurface->nst][3]; // should be int8 eventually?
    dde->pull(spost_global,sm->surface_spost);

    delete dde; dde = NULL;

    // use the global map to build a local node index...

    assert(subSurface->spost == NULL); subSurface->spost = new int[subSurface->nst][3];
    for (int ist = 0; ist < subSurface->nst; ++ist) {
      const int st_bits = int(subSurface->ist_global_and_bits[ist]>>52);
      FOR_I3 {
        int8 isp_global_and_bits = spost_global[ist][i];
        if (st_bits) isp_global_and_bits |= int8(st_bits)<<52;
        map<const int8,int>::iterator iter = globalMap.find(isp_global_and_bits);
        if (iter == globalMap.end()) {
          // this must be an old mles file...
          //assert(sm->surface_spobf_i == NULL);
          globalMap[isp_global_and_bits] = subSurface->spost[ist][i] = subSurface->nsp++;
        }
        else {
          subSurface->spost[ist][i] = iter->second;
        }
      }
    }
    delete[] spost_global;

    // record the global surface isp index...

    assert(subSurface->isp_global_and_bits == NULL); // TODO: this should be private with accessors
    subSurface->isp_global_and_bits = new int8[subSurface->nsp];
    for (map<const int8,int>::iterator iter = globalMap.begin(); iter != globalMap.end(); ++iter) {
      subSurface->isp_global_and_bits[iter->second] = iter->first;
    }

    globalMap.clear();

    // don't know nsp size for old mesh until now
    //if (sm->surface_spobf_i == NULL) {
      delete[] bits;
      bits = new int[subSurface->nsp];
    //}

    // and finally, bring over the coordinates...

    // strip off any bits and use subSurface->isp_global_and_bits as an
    // index into sm->surface_spora...
    for (int isp = 0; isp < subSurface->nsp; ++isp) {
      bits[isp] = int(subSurface->isp_global_and_bits[isp]>>52);
      if (bits[isp])
        subSurface->isp_global_and_bits[isp] &= MASK_52BITS;
    }

    dde = new DistributedDataExchanger(subSurface->isp_global_and_bits,subSurface->nsp,sm->surface_spora);

    assert(subSurface->xp == NULL); subSurface->xp = new double[subSurface->nsp][3];
    dde->pull(subSurface->xp,sm->surface_xp);

    delete dde;

    for (int isp = 0; isp < subSurface->nsp; ++isp) {
      if (bits[isp]) {
        subSurface->isp_global_and_bits[isp] |= (int8(bits[isp])<<52);
        PeriodicData::periodicTranslate(subSurface->xp[isp],1,bits[isp]);
      }
    }

    delete[] bits;

    // take a look...

    if (checkParam("WRITE_SUBSURFACE"))
      subSurface->writeTecplot();

  }

  // ===================================================================
  // extended faces...
  // ===================================================================

  // the ef's need to know where they are going...

  uint8 (*rbioef_sm)[2] = new uint8[sm->nef][2];
  {

    int * bits = new int[sm->nef];

    for (int ief = 0; ief < sm->nef; ++ief) {
      bits[ief] = int(sm->cvoef_global[ief][1]>>52);
      if (bits[ief]) {
        sm->cvoef_global[ief][1] &= MASK_52BITS;
      }
    }

    DistributedDataExchanger dde((int8*)sm->cvoef_global,sm->nef*2,sm->cvora);
    dde.pull((uint8*)rbioef_sm,rbi_sm);

    for (int ief = 0; ief < sm->nef; ++ief) {
      if (bits[ief]) {
        sm->cvoef_global[ief][1] |= (int8(bits[ief])<<52);
        int rank1,bits1,index1;
        BitUtils::unpackRankBitsIndex(rank1,bits1,index1,rbioef_sm[ief][1]);
        assert(bits1 == 0); // no periodicity for now
        rbioef_sm[ief][1] = BitUtils::packRankBitsIndex(rank1,bits[ief],index1);
      }
    }

    delete[] bits;

  }

  // cleanup...

  //delete[] rbi_sm; // keep for io?

  // now count...

  FOR_RANK send_count[rank] = 0;

  for (int ief = 0; ief < sm->nef; ++ief) {
    int rank0,bits0,index0;
    BitUtils::unpackRankBitsIndex(rank0,bits0,index0,rbioef_sm[ief][0]);
    assert(bits0 == 0);
    int rank1,bits1,index1;
    BitUtils::unpackRankBitsIndex(rank1,bits1,index1,rbioef_sm[ief][1]);
    if ((rank0 == rank1)&&(bits1 == 0)) { // note ef's are different -- we don't have periodic shadows
      // send once to rank0...
      ++send_count[rank0];
    }
    else {
      ++send_count[rank0];
      ++send_count[rank1];
    }
  }

  // allocate...

  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
  send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

  assert(send_buf_int == NULL); send_buf_int = new int[send_count_sum*4];
  assert(send_buf_double == NULL); send_buf_double = new double[send_count_sum*6];

  // pack...

  for (int ief = 0; ief < sm->nef; ++ief) {
    int rank0,bits0,index0;
    BitUtils::unpackRankBitsIndex(rank0,bits0,index0,rbioef_sm[ief][0]);
    assert(bits0 == 0);
    int rank1,bits1,index1;
    BitUtils::unpackRankBitsIndex(rank1,bits1,index1,rbioef_sm[ief][1]);
    if ((rank0 == rank1)&&(bits1 == 0)) {
      // send once to rank0==rank1...
      send_buf_int[send_disp[rank0]*4  ] = index0; // cvoef[ief][0] on the final partition
      send_buf_int[send_disp[rank0]*4+1] = rank1;  // nbr info
      send_buf_int[send_disp[rank0]*4+2] = bits1;
      send_buf_int[send_disp[rank0]*4+3] = index1;
      // geometry...
      send_buf_double[send_disp[rank0]*6  ]  = sm->n_ef[ief][0];
      send_buf_double[send_disp[rank0]*6+1]  = sm->n_ef[ief][1];
      send_buf_double[send_disp[rank0]*6+2]  = sm->n_ef[ief][2];
      send_buf_double[send_disp[rank0]*6+3]  = sm->c_ef[ief][0];
      send_buf_double[send_disp[rank0]*6+4]  = sm->c_ef[ief][1];
      send_buf_double[send_disp[rank0]*6+5]  = sm->c_ef[ief][2];
      ++send_disp[rank0];
    }
    else {
      // send once to rank0==rank1...
      send_buf_int[send_disp[rank0]*4  ] = index0; // cvoef[ief][0] on the final partition
      send_buf_int[send_disp[rank0]*4+1] = rank1;  // nbr info
      send_buf_int[send_disp[rank0]*4+2] = bits1;
      send_buf_int[send_disp[rank0]*4+3] = index1;
      // geometry...
      send_buf_double[send_disp[rank0]*6  ]  = sm->n_ef[ief][0];
      send_buf_double[send_disp[rank0]*6+1]  = sm->n_ef[ief][1];
      send_buf_double[send_disp[rank0]*6+2]  = sm->n_ef[ief][2];
      send_buf_double[send_disp[rank0]*6+3]  = sm->c_ef[ief][0];
      send_buf_double[send_disp[rank0]*6+4]  = sm->c_ef[ief][1];
      send_buf_double[send_disp[rank0]*6+5]  = sm->c_ef[ief][2];
      ++send_disp[rank0];
      // and rank1 - flip...
      send_buf_int[send_disp[rank1]*4  ] = index1; // cvoef[ief][0] on the final partition
      send_buf_int[send_disp[rank1]*4+1] = rank0;  // nbr info
      const int inv_bits1 = BitUtils::flipPeriodicBits(bits1);
      send_buf_int[send_disp[rank1]*4+2] = inv_bits1;
      send_buf_int[send_disp[rank1]*4+3] = index0;
      // geometry...
      send_buf_double[send_disp[rank1]*6  ]  = -sm->n_ef[ief][0];
      send_buf_double[send_disp[rank1]*6+1]  = -sm->n_ef[ief][1];
      send_buf_double[send_disp[rank1]*6+2]  = -sm->n_ef[ief][2];
      send_buf_double[send_disp[rank1]*6+3]  = sm->c_ef[ief][0];
      send_buf_double[send_disp[rank1]*6+4]  = sm->c_ef[ief][1];
      send_buf_double[send_disp[rank1]*6+5]  = sm->c_ef[ief][2];
      if (inv_bits1) PeriodicData::periodicRotate(send_buf_double+send_disp[rank1]*6,2,inv_bits1);
      ++send_disp[rank1];
    }
  }

  delete[] rbioef_sm;

  // rewind disp...

  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

  // exchange...

  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
  recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

  FOR_RANK {
    send_count[rank] *= 4;
    send_disp[rank] *= 4;
    recv_count[rank] *= 4;
    recv_disp[rank] *= 4;
  }

  assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum*4];
  MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
      recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
  delete[] send_buf_int; send_buf_int = NULL;

  FOR_RANK {
    send_count[rank] = (send_count[rank]/4)*6;
    send_disp[rank] = (send_disp[rank]/4)*6;
    recv_count[rank] = (recv_count[rank]/4)*6;
    recv_disp[rank] = (recv_disp[rank]/4)*6;
  }

  assert(recv_buf_double == NULL); recv_buf_double = new double[recv_count_sum*6];
  MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
      recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
  delete[] send_buf_double; send_buf_double = NULL;

  // unpack and sort on the recv side...
  // ef's are grouped as follows:
  //
  // these are designed to produce a set of sorted ef's:
  // 1. icv (local)
  // 2. ief_global (through irecv -- tricky!)

  vector< pair<int,int> > internalEfVec;    // first = icv (local, for sort), second.first = irecv, second.second = irecv2
  vector< pair<int,int> > interprocessorEfVec;

  for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
    // grab just enough data to figure out what kind of ef this is...
    const int icv  = recv_buf_int[irecv*4  ]; assert((icv >= 0)&&(icv < ncv));
    const int rank = recv_buf_int[irecv*4+1];
    const int bits = recv_buf_int[irecv*4+2];
    if ((rank == mpi_rank)&&(bits == 0)) {
      // internal face...
      internalEfVec.push_back( pair<int,int>(icv,irecv) );
    }
    else {
      // interprocessor face...
      interprocessorEfVec.push_back( pair<int,int>(icv,irecv) );
    }
  }

  // sort...

  sort(internalEfVec.begin(),internalEfVec.end());
  sort(interprocessorEfVec.begin(),interprocessorEfVec.end());

  // and set ef data...

  nef_i = internalEfVec.size();
  nef = nef_i + interprocessorEfVec.size();

  // allocate ef data...

  assert(cvoef == NULL); cvoef = new int[nef][2];
  assert(n_ef == NULL);  n_ef = new double[nef][3];
  assert(c_ef == NULL);  c_ef = new double[nef][3];

  // done later now...

  assert(group_ef == NULL);

  // unpack sorted data...

  // internal first...

  for (int ii = 0, nii = internalEfVec.size(); ii < nii; ++ii) {
    const int irecv = internalEfVec[ii].second;
    // the icv...
    const int icv = internalEfVec[ii].first;
    // unpack...
    assert(icv       == recv_buf_int[irecv*4  ]); assert((icv >= 0)&&(icv < ncv));
    const int rank    = recv_buf_int[irecv*4+1]; assert(rank == mpi_rank);
    const int bits    = recv_buf_int[irecv*4+2]; assert(bits == 0);
    const int icv_nbr = recv_buf_int[irecv*4+3]; assert((icv_nbr >= 0)&&(icv_nbr < ncv));
    cvoef[ii][0] = icv;
    cvoef[ii][1] = icv_nbr;
    n_ef[ii][0] = recv_buf_double[irecv*6  ];
    n_ef[ii][1] = recv_buf_double[irecv*6+1];
    n_ef[ii][2] = recv_buf_double[irecv*6+2];
    c_ef[ii][0] = recv_buf_double[irecv*6+3];
    c_ef[ii][1] = recv_buf_double[irecv*6+4];
    c_ef[ii][2] = recv_buf_double[irecv*6+5];
  }

  internalEfVec.clear();

  // then interproc...
  // note that the extended faces will touch cells in the current ghost
  // range, as well as new cells. These are added as ncv_g2 using the
  // same Map...

  ncv_g2 = ncv_g;

  for (int ii = 0, nii = interprocessorEfVec.size(); ii < nii; ++ii) {
    const int irecv = interprocessorEfVec[ii].second;
    // the icv...
    const int icv = interprocessorEfVec[ii].first;
    // unpack...
    assert(icv       == recv_buf_int[irecv*4  ]); assert((icv >= 0)&&(icv < ncv));
    const int rank    = recv_buf_int[irecv*4+1];
    const int bits    = recv_buf_int[irecv*4+2];
    const int icv_nbr = recv_buf_int[irecv*4+3];
    cvoef[ii+nef_i][0] = icv;
    const uint8 rbi_nbr = BitUtils::packRankBitsIndex(rank,bits,icv_nbr);
    map<const uint8,int>::iterator iter = rbiMap.find(rbi_nbr);
    if (iter == rbiMap.end()) {
      rbiMap[rbi_nbr] = cvoef[ii+nef_i][1] = ncv_g2++;
    }
    else {
      cvoef[ii+nef_i][1] = iter->second;
    }
    n_ef[ii+nef_i][0] = recv_buf_double[irecv*6  ];
    n_ef[ii+nef_i][1] = recv_buf_double[irecv*6+1];
    n_ef[ii+nef_i][2] = recv_buf_double[irecv*6+2];
    c_ef[ii+nef_i][0] = recv_buf_double[irecv*6+3];
    c_ef[ii+nef_i][1] = recv_buf_double[irecv*6+4];
    c_ef[ii+nef_i][2] = recv_buf_double[irecv*6+5];
  }

  interprocessorEfVec.clear();

  delete[] recv_buf_int; recv_buf_int = NULL;
  delete[] recv_buf_double; recv_buf_double = NULL;
  delete[] recv_buf_int8; recv_buf_int8 = NULL;

  // efs should be done!

  // ===================================================
  // gcl checks....
  // ===================================================

  {

    double (*gcl)[3] = new double[ncv][3];
    int * count = new int[ncv];

    // compact gcl...

    FOR_ICV {
      FOR_I3 gcl[icv][i] = 0.0;
      count[icv] = 0;
    }

    FOR_IBF {
      const int icv = cvobf[ibf]; assert((icv >= 0)&&(icv < ncv));
      FOR_I3 gcl[icv][i] += n_bf[ibf][i];
    }

    FOR_INTERNAL_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
      FOR_I3 gcl[icv0][i] += n_fa[ifa][i];
      FOR_I3 gcl[icv1][i] -= n_fa[ifa][i];
      ++count[icv0];
      ++count[icv1];
    }

    FOR_INTERPROC_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g)); // must be in compact ghosts
      FOR_I3 gcl[icv0][i] += n_fa[ifa][i];
      ++count[icv0];
    }

    MiscUtils::dumpRange(gcl,ncv,"GCL compact faces");
    MiscUtils::dumpRange(count,ncv,"compact faces per cv");

    // extended face gcl...

    if (nef_global > 0) {

      FOR_ICV {
        FOR_I3 gcl[icv][i] = 0.0;
        count[icv] = 0;
      }

      FOR_IBF {
        const int icv = cvobf[ibf]; assert((icv >= 0)&&(icv < ncv));
        FOR_I3 gcl[icv][i] += n_bf[ibf][i];
      }

      FOR_INTERNAL_IEF {
        const int icv0 = cvoef[ief][0]; assert((icv0 >= 0)&&(icv0 < ncv));
        const int icv1 = cvoef[ief][1]; assert((icv1 >= 0)&&(icv1 < ncv));
        FOR_I3 gcl[icv0][i] += n_ef[ief][i];
        FOR_I3 gcl[icv1][i] -= n_ef[ief][i];
        ++count[icv0];
        ++count[icv1];
      }

      FOR_INTERPROC_IEF {
        const int icv0 = cvoef[ief][0]; assert((icv0 >= 0)&&(icv0 < ncv));
        const int icv1 = cvoef[ief][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g2)); // must be in either compact or extended ghosts
        FOR_I3 gcl[icv0][i] += n_ef[ief][i];
        ++count[icv0];
      }

      MiscUtils::dumpRange(gcl,ncv,"GCL extended faces");
      MiscUtils::dumpRange(count,ncv,"extended faces per cv");

    }

    delete[] gcl;
    delete[] count;

  }

  // ================================================
  // cv-based geometry...
  // copy x_vv from striped mesh. For this quantity,
  // we will also set the ghost values...
  // ================================================

  assert(x_vv == NULL);
  x_vv = new double[ncv_g2][3];
  FOR_ICV FOR_I3 x_vv[icv][i] = sm->x_vv[icv][i]; // recall x_vv already reordered in sm
  for (int icv = ncv; icv < ncv_g; ++icv) FOR_I3 x_vv[icv][i] = 1.0E+20; // we need the communicators to populate this -- see below

  assert(r_vv == NULL);
  r_vv = new double[ncv_g];
  assert(sm->r_vv);
  dde_cv_striped->pull(r_vv,sm->r_vv);
  for (int icv = ncv; icv < ncv_g; ++icv) r_vv[icv] = HUGE_VAL;

  assert(x_cv == NULL);
  //x_cv = new double[ncv_g][3];
  x_cv = new double[ncv_g2][3]; // 2-level ghosts used by FlowSolver and VofSolver
  assert(sm->x_cv);
  dde_cv_striped->pull(x_cv,sm->x_cv);

  assert(vol_cv == NULL);
  //vol_cv = new double[ncv_g];
  vol_cv = new double[ncv_g2]; // 2-level ghosts required for a_sgs
  assert(sm->vol_cv);
  dde_cv_striped->pull(vol_cv,sm->vol_cv);

  assert( inv_vol == NULL);
  inv_vol = new double[ncv];
  for (int icv = 0; icv < ncv; ++icv)
    inv_vol[icv] = 1.0/vol_cv[icv];

  // ================================================
  // now build the communicators: Cv first...
  // ================================================

  // the Map will automatically sort the ghost data into the rank-bits-index order.
  // Because a single map was used, the first and second layers will be
  // inter-mixed. We want these separated in the communicators: to allow reductions on
  // either first-layer cvs or on both first-and-second-layer cvs...

  assert(int(rbiMap.size()) == ncv_g2-ncv);
  int * icv_new_g2 = new int[ncv_g2-ncv];
  assert(rbi_g == NULL); rbi_g = new uint8[ncv_g-ncv];
  assert(rbi_g2 == NULL); rbi_g2 = new uint8[ncv_g2-ncv_g];
  int icv_g = ncv;
  int icv_g2 = ncv_g;
  for (map<const uint8,int>::const_iterator iter = rbiMap.begin(); iter != rbiMap.end(); ++iter) {
    // the second should contain the current ghost icv...
    assert(iter->second >= ncv);
    int icv_new;
    if (iter->second < ncv_g) {
      icv_new = icv_g++; // level 1 ghost
      rbi_g[icv_new-ncv] = iter->first;
    }
    else {
      assert(iter->second < ncv_g2);
      icv_new = icv_g2++; // level 2 ghost
      rbi_g2[icv_new-ncv_g] = iter->first;
    }
    icv_new_g2[iter->second-ncv] = icv_new;
  }

  rbiMap.clear();

  // change cvofa and cvoef based on this new ghost ordering...

  FOR_INTERPROC_IFA {
    const int icv1 = cvofa[ifa][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g)); // must be in compact ghosts
    const int icv1_new = icv_new_g2[icv1-ncv]; assert((icv1_new >= ncv)&&(icv1_new < ncv_g)); // same
    cvofa[ifa][1] = icv1_new;
  }

  FOR_INTERPROC_IEF {
    const int icv1 = cvoef[ief][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g2)); // must be in either compact or extended ghosts
    const int icv1_new = icv_new_g2[icv1-ncv]; assert((icv1_new >= ncv)&&(icv1_new < ncv_g2)); // same
    cvoef[ief][1] = icv1_new;
  }

  delete[] icv_new_g2;

  // build level-1 cvPrcomm

  buildCvPrcomm();

  // update x_vv...

  updateCvData(x_vv,REPLACE_TRANSLATE_DATA);
  updateCvData(x_cv,REPLACE_TRANSLATE_DATA);
  updateCvData(vol_cv);
  updateCvData(r_vv);
  MiscUtils::dumpRange(vol_cv,ncv_g,"vol_cv");

  // a minimal check is confirming that the ghost data has been over-written...

  for (int icv = ncv; icv < ncv_g; ++icv) FOR_I3 assert(x_vv[icv][i] != 1.0E+20);

  // build per_R and per_t (needed by lsp)...

  buildPeriodicRt();

  // we can now build area_over_delta_fa...

  assert(area_over_delta_fa == NULL);
  area_over_delta_fa = new double[nfa];

  // this version is based on voronoi forming point locations.
  // tests indicate it is better to use the below based on cvs centroids...
  /*
     FOR_IFA {
     const int icv0 = cvofa[ifa][0];
     const int icv1 = cvofa[ifa][1];
     const double area = MAG(n_fa[ifa]);
     const double delta = DIST(x_vv[icv0],x_vv[icv1]);
     area_over_delta_fa[ifa] = area/delta;
     }
     */

  // this version based on cv centroids...
  FOR_IFA {
    area_over_delta_fa[ifa] = 0.0;
    const double area2 = DOT_PRODUCT(n_fa[ifa],n_fa[ifa]);
    if (area2 > 0.0) {
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      const double dx[3] = DIFF(x_cv[icv1],x_cv[icv0]);
      const double dp = DOT_PRODUCT(dx,n_fa[ifa]); // this has magnitude delta*area
      if (dp > 0.0) {
        area_over_delta_fa[ifa] = area2/dp;
      }
    }
  }

  // build level-2 cvPrcomm2

  buildCv2Prcomm();
  updateCv2Data(x_vv,REPLACE_TRANSLATE_DATA);
  updateCv2Data(x_cv,REPLACE_TRANSLATE_DATA);
  updateCv2Data(vol_cv,REPLACE_DATA);

  // extended operator check...

  if (nef_global > 0)
    checkGradExtended();

  // finally fa and ef grouping...

  assert(group_fa == NULL);
  group_fa = new int[nfa];
  FOR_IFA group_fa[ifa] = -1;

  assert(group_ef == NULL);
  group_ef = new int[nef];
  FOR_IEF group_ef[ief] = -1;

  // build fa dde and striping for face writing...

  assert(dde_fa_striped == NULL);
  int8 *abs_ifa_global = new int8[nfa];
  FOR_IFA abs_ifa_global[ifa] = max(ifa_global[ifa],-ifa_global[ifa]-1);
  dde_fa_striped = new DistributedDataExchanger(abs_ifa_global,nfa,sm->faora);
  delete[] abs_ifa_global; // ifa_global is signed w/ -1 indexing for loser interproc faces!
  // also take a copy of the sm->faora...
  assert(faora_striped == NULL);
  faora_striped = new int8[mpi_size+1];
  for (int rank = 0; rank <= mpi_size; ++rank) faora_striped[rank] = sm->faora[rank];

  // particle stuff...

  assert(rbi_sm);
  for (int ii_sm = 0, lim = sm->lpHelperVec.size(); ii_sm < lim; ++ii_sm) {
    // find lp helper...
    map<const string,int>::iterator iter = lpHelperNameMap.find(sm->lpHelperVec[ii_sm].name);
    assert(iter != lpHelperNameMap.end());
    const int ii = iter->second;

    // setup send side stuff...

    FOR_RANK send_count[rank] = 0;
    int ncv_striped = cvora_striped[mpi_rank+1]-cvora_striped[mpi_rank];
    for (int icv_striped = 0; icv_striped < ncv_striped; ++icv_striped) {
      int rank,bits,icv;
      BitUtils::unpackRankBitsIndex(rank,bits,icv,rbi_sm[icv_striped]);
      assert(bits == 0);
      // icv,count and ilps...
      send_count[rank] += 2 + sm->lpHelperVec[ii_sm].lpocv_i_global[icv_striped+1]-sm->lpHelperVec[ii_sm].lpocv_i_global[icv_striped];
    }

    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

    // build lpora...

    MiscUtils::buildXora(lpHelperVec[ii].lpora_cv_striped,sm->lpHelperVec[ii_sm].nlp);

    assert(send_buf_int == NULL); send_buf_int = new int[send_count_sum];
    for (int icv_striped = 0; icv_striped < ncv_striped; ++icv_striped) {
      int rank,bits,icv;
      BitUtils::unpackRankBitsIndex(rank,bits,icv,rbi_sm[icv_striped]);
      assert(bits == 0);
      // pack icv (promoted to int8) and ip_globals associated to icv...
      send_buf_int[send_disp[rank]++] = icv;
      send_buf_int[send_disp[rank]++] = sm->lpHelperVec[ii_sm].lpocv_i_global[icv_striped+1]-sm->lpHelperVec[ii_sm].lpocv_i_global[icv_striped];
      for (int loc = sm->lpHelperVec[ii_sm].lpocv_i_global[icv_striped]; loc < sm->lpHelperVec[ii_sm].lpocv_i_global[icv_striped+1]; ++loc)
        send_buf_int[send_disp[rank]++] = loc-lpHelperVec[ii].lpora_cv_striped[mpi_rank];
    }

    // rewind...

    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

    // setup recv side stuff...

    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    // exchange...

    assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
        recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; send_buf_int = NULL;

    int nlp = 0;
    FOR_RANK {
      int irecv = recv_disp[rank];
      while (irecv < recv_disp[rank]+recv_count[rank]) {
        //const int icv = recv_buf_int[irecv++];
        ++irecv;
        const int nloc = recv_buf_int[irecv++];
        nlp += nloc;
        irecv += nloc;
      }
    }

    *lpHelperVec[ii].size_ptr = nlp; // update registered data size

    int8* ilp_global = new int8[nlp];
    lpHelperVec[ii].cvolp = new int[nlp];
    nlp = 0;
    FOR_RANK {
      int irecv = recv_disp[rank];
      while (irecv < recv_disp[rank]+recv_count[rank]) {
        const int icv = recv_buf_int[irecv++];
        const int nloc = recv_buf_int[irecv++];
        for (int loc = 0; loc < nloc; ++loc) {
          ilp_global[nlp] = recv_buf_int[irecv++]+lpHelperVec[ii].lpora_cv_striped[rank];
          lpHelperVec[ii].cvolp[nlp] = icv;
          ++nlp;
        }
      }
    }
    assert(*lpHelperVec[ii].size_ptr == nlp);
    delete[] recv_buf_int; recv_buf_int = NULL;

    // build dde...
    lpHelperVec[ii].dde_cv_striped = new DistributedDataExchanger(ilp_global,nlp,lpHelperVec[ii].lpora_cv_striped);
    delete[] ilp_global;
  }

  delete[] send_count;
  delete[] send_disp;
  delete[] recv_count;
  delete[] recv_disp;

}

void StaticSolver::buildPeriodicRt() {

  // TODO: replace this with a more fundamental handling of
  // periodicity for particles, e.g. as part of PeriodicData
  // itself...

  assert(rbi_g != NULL);
  assert(iRt_g == NULL); iRt_g = new int[ncv_g-ncv][2];

  double R[9],t[3];
  int count_pR = 0, count_pt = 0;
  for (int iter = 0; iter < 2; ++iter) {
    int rank_current = -1;
    int bits_current = -1;
    for (int icv = ncv; icv < ncv_g; ++icv) {
      if (iter == 0) FOR_I2 iRt_g[icv-ncv][i] = -1; // init to no R and no t
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv-ncv]);
      if (rank > rank_current) {
        bits_current = -1;
      }
      else {
        assert(rank == rank_current);
      }
      if (bits > bits_current) {
        // look for new rotations/translations...
        if (iter == 0) {
          const bool has_R = PeriodicData::getPeriodicR(R,bits);
          const bool has_t = PeriodicData::getPeriodicT(t,bits);
          if (has_R) {
            iRt_g[icv-ncv][0] = count_pR++;
          }
          if (has_t) {
            iRt_g[icv-ncv][1] = count_pt++;
          }
        }
        else {
          assert(iter == 1);
          const int ipR = iRt_g[icv-ncv][0];
          const int ipt = iRt_g[icv-ncv][1];
          if (ipR >= 0) {
            const bool has_R = PeriodicData::getPeriodicR(R,bits);
            FOR_I3 FOR_J3 per_R[ipR][3*i+j] = R[3*j+i]; // transpose
            assert(has_R);
          }
          if (ipt >= 0) {
            const bool has_t = PeriodicData::getPeriodicT(t,bits);
            FOR_I3 per_t[ipt][i] = -t[i]; // negate
            assert(has_t);
          }
        }
        bits_current = bits;
      }
      else {
        assert(bits_current == bits);
      }
    }

    if (iter == 0) {
      assert(per_R == NULL); per_R = new double[count_pR][9];
      assert(per_t == NULL); per_t = new double[count_pt][3];
    }
  }

}

void StaticSolver::groupExtendedFaces() const {

  // create vector of pairs of area and ief for ef's with c approx 0

  const double c2_over_n2_ratio_zero = 1.0E-24;
  const double matched_face_tol = 1.0E-11;
  const double small_face_tol = 1.0E-6;

  int8 my_count[6] = { nef, 0, 0, 0, 0, 0 };

  vector<pair<double,int> > magAndIndexVec;
  for (int ief = 0; ief < nef; ++ief) {
    // check magnitudes...
    const double n_mag2 = DOT_PRODUCT(n_ef[ief],n_ef[ief]);
    const double c_mag2 = DOT_PRODUCT(c_ef[ief],c_ef[ief]);
    // only consider faces where c is zero wrt n_mag. Use the same factor
    // as used when determining zero faces in the vd tightened by an additional
    // few orders of magnitude...
    if (c_mag2 <= c2_over_n2_ratio_zero*n_mag2) {
      magAndIndexVec.push_back(pair<double,int>(sqrt(n_mag2),ief));
      ++my_count[1];
    }
    else if (c_mag2 < 10.0*c2_over_n2_ratio_zero*n_mag2) {
      ++my_count[2];
    }
    else if (c_mag2 < 100.0*c2_over_n2_ratio_zero*n_mag2) {
      ++my_count[3];
    }
    else if (c_mag2 < 1000.0*c2_over_n2_ratio_zero*n_mag2) {
      ++my_count[4];
    }
    else if (c_mag2 < 10000.0*c2_over_n2_ratio_zero*n_mag2) {
      ++my_count[5];
    }
  }

  int8 count[6];
  MPI_Reduce(my_count,count,6,MPI_INT8,MPI_SUM,0,mpi_comm);
  if (mpi_rank == 0) {
    cout << " extended faces flagged for potential grouping: " << count[1] << " out of " << count[0] << " (" << double(count[1])/double(count[0])*100.0 << "%)" << endl;
    cout << " sensitivity check: additional extended faces flagged if c2_over_n2_ratio_zero reduced: " << count[2] << " " << count[3] << " " << count[4] << " " << count[5] << endl;
  }

  // sort c approx 0 faces by area...

  sort(magAndIndexVec.begin(),magAndIndexVec.end());

  // using the sorted vec, bin the areas based off first element and count number in bin.
  // group (A-A_0) < eps_A * A_0.

  vector<pair<double,int> > myNormalKindVec; // first: A_0, second: count for (A-A0) < eps_A * A_0
  for (int ii = 0,nii = magAndIndexVec.size(); ii < nii; ++ii) {
    if ( myNormalKindVec.empty() || ( magAndIndexVec[ii].first - myNormalKindVec.back().first > matched_face_tol*myNormalKindVec.back().first ) ) {
      myNormalKindVec.push_back(pair<double,int>(magAndIndexVec[ii].first,1));
    }
    else {
      myNormalKindVec.back().second += 1; // the count
    }
  }

  /*
     for (int ii = 0,nii = myNormalKindVec.size(); ii < nii; ++ii) {
     if (ii == 0) {
     cout << " > group: " << ii << " area: " << myNormalKindVec[ii].first << " count: " << myNormalKindVec[ii].second << endl;
     }
     else {
     cout << " > group: " << ii << " area: " << myNormalKindVec[ii].first << " count: " << myNormalKindVec[ii].second <<
     " diff with prev: " << (myNormalKindVec[ii].first-myNormalKindVec[ii-1].first)/myNormalKindVec[ii-1].first << endl;
     }
     }
     */

  // now repeat the grouping using the median value.
  // group (A-A_med) < eps_A * A_med. TODO

  int my_ngr = myNormalKindVec.size();
  int ngr;
  //MPI_Reduce(&my_ngr,&ngr,1,MPI_INT,MPI_SUM,0,mpi_comm);
  //if (mpi_rank == 0)
  //  cout << "Initial extended face groups: " << ngr << endl;

  // only keep groups greater than a threshold...

  vector<double> magVec;
  int8 my_final_count = 0;
  const int ngr_threshold = getIntParam("GROUP_MEMBER_THRESHOLD", 100); // TODO move
  for (int ii = 0,nii = myNormalKindVec.size(); ii < nii; ++ii) {
    // save groups with counts above 100. Note that at this point we
    // only have the areas grouped. In each area, there could be
    // 12 or 14 directions. So assuming we need atleast a few faces
    // in each area, we threshold at 100...
    //
    // also force the areas to larger than some epsilon ....
    if ( (myNormalKindVec[ii].second > ngr_threshold) && (myNormalKindVec[ii].first > small_face_tol*myNormalKindVec[nii-1].first) ) {
      //cout << " > FINAL extended face group: " << magVec.size() << " " << myNormalKindVec[ii].first << " count: " << myNormalKindVec[ii].second << endl;
      magVec.push_back(myNormalKindVec[ii].first);
      my_final_count += myNormalKindVec[ii].second;
    }
  }

  int8 final_count;
  MPI_Reduce(&my_final_count,&final_count,1,MPI_INT8,MPI_SUM,0,mpi_comm);
  if (mpi_rank == 0)
    cout << "Number of grouped extended faces: " << final_count << endl;

  my_ngr = magVec.size();
  MPI_Reduce(&my_ngr,&ngr,1,MPI_INT,MPI_MAX,0,mpi_comm);
  if (mpi_rank == 0)
    cout << "Max number of serial area-based extended face groups: " << ngr << endl;

  double my_max_dp = 0.0;
  int my_max_gr = 0;
  if (my_ngr > 0) {

    // accuracy that a unit-normal dot product must be aligned to be considered the same face...

    const double dp_eps = 1.0E-12;

    vector<double> doubleVec;
    int idx_offset = 0;
    int igr = 0;
    for (int ii = 0,nii = magAndIndexVec.size(); ii < nii; ++ii) {
      while ((igr < int(magVec.size())-1)&&(magVec[igr+1] <= magAndIndexVec[ii].first)) {
        ++igr;
        doubleVec.push_back(0.0); doubleVec.push_back(0.0); doubleVec.push_back(0.0);
        idx_offset = doubleVec.size();
      }
      assert(igr < (1<<21)); // igr gets shifted by 10 below, and cannot exceed the int range...
      // it is possible that a group included in the magVec just fit, so its min value was (1+1*matched_face_tol)*magVec[igr]. Members
      // in that group may have been up to a full matched_face_tol away, so here we use 2x to include these candidates...
      if ( (magAndIndexVec[ii].first - magVec[igr] >= 0.0) && (magAndIndexVec[ii].first - magVec[igr] <= 2.0*matched_face_tol*magVec[igr]) ) {
        // get the unit normal...
        const int ief = magAndIndexVec[ii].second;
        assert(group_ef[ief] == -1);
        const double unit_n[3] = {
          n_ef[ief][0]/magAndIndexVec[ii].first,
          n_ef[ief][1]/magAndIndexVec[ii].first,
          n_ef[ief][2]/magAndIndexVec[ii].first
        };
        // check the dot product against existing vectors in this group...
        int idx, idx_end;
        for (idx = idx_offset, idx_end = doubleVec.size(); idx < idx_end; idx += 3) {
          const double dp = unit_n[0]*doubleVec[idx] + unit_n[1]*doubleVec[idx+1] + unit_n[2]*doubleVec[idx+2];
          if (dp > 1.0-dp_eps) {
            my_max_dp = max(my_max_dp,1.0-dp);
            // this is aligned with this vector...
            assert(idx%3 == 0);
            // the bit rules for the group are:
            // msb: group, then member in group, finally bit0 = sign...
            group_ef[ief] = (igr<<10)|(((idx-idx_offset)/3)<<1);
            break;
          }
          else if (dp < -1.0+dp_eps) {
            my_max_dp = max(my_max_dp,1.0+dp);
            // this is exactly misaligned with this vector...
            assert(idx%3 == 0);
            group_ef[ief] = (igr<<10)|(((idx-idx_offset)/3)<<1)|1; // 0-bit is the sign bit
            break;
          }
        }
        if (idx == int(doubleVec.size())) {
          assert(idx%3 == 0);
          //cout << (idx-idx_offset)/3 << endl;
          if (!((idx-idx_offset)/3 < (1<<9)))
            cout << "Error: (idx-idx_offset)/3: " << (idx-idx_offset)/3 << endl;
          assert((idx-idx_offset)/3 < (1<<9)); // might have to increase one day from from 16 to 32
          my_max_gr = max(my_max_gr,(idx-idx_offset)/3);
          group_ef[ief] = (igr<<10)|(((idx-idx_offset)/3)<<1);
          doubleVec.push_back(unit_n[0]);
          doubleVec.push_back(unit_n[1]);
          doubleVec.push_back(unit_n[2]);
        }
      }
    }

  }

  /*
     for (int ief = 0; ief < nef; ++ief) {
     cout << "GN_EF" << group_ef[ief] << " ";
     FOR_I3 cout << n_ef[ief][i] << " ";
     cout << endl;
     }
     */

  int max_gr;
  MPI_Reduce(&my_max_gr,&max_gr,1,MPI_INT,MPI_MAX,0,mpi_comm);
  if (mpi_rank == 0)
    cout << "Max number of serial orientation-based extended face sub-groups (expect 4:12 in 2d and 8:16 in 3d): " << max_gr+1 << endl;

  double max_dp;
  MPI_Reduce(&my_max_dp,&max_dp,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
  if (mpi_rank == 0)
    cout << "Actual tol for dp matching for extended faces: " << max_dp << endl;

}

void StaticSolver::groupCompactFaces() const {

  // create vector of pairs of area and ifa for fa's with c approx 0

  const double c4_over_n2_ratio_zero = 1.0E-24; // using a length scale rather than area (like extended)
  const double matched_face_tol = 1.0E-11;
  const double small_face_tol = 1.0E-6;

  int8 my_count[6] = { nfa, 0, 0, 0, 0, 0 };

  vector<pair<double,int> > magAndIndexVec;
  for (int ifa = 0; ifa < nfa; ++ifa) {
    // check magnitudes...
    const double n_mag2 = DOT_PRODUCT(n_fa[ifa],n_fa[ifa]);
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    double c_fa[3]; FOR_I3 c_fa[i] = 0.5*(x_cv[icv0][i]+x_cv[icv1][i])-x_fa[ifa][i];
    const double c_mag2 = DOT_PRODUCT(c_fa,c_fa);;
    const double c_mag4 = c_mag2*c_mag2;
    // only consider faces where c is zero wrt delta.
    if (c_mag4 <= c4_over_n2_ratio_zero*n_mag2) {
      magAndIndexVec.push_back(pair<double,int>(sqrt(n_mag2),ifa));
      ++my_count[1];
    }
    else if (c_mag4 < 10.0*c4_over_n2_ratio_zero*n_mag2) {
      ++my_count[2];
    }
    else if (c_mag4 < 100.0*c4_over_n2_ratio_zero*n_mag2) {
      ++my_count[3];
    }
    else if (c_mag4 < 1000.0*c4_over_n2_ratio_zero*n_mag2) {
      ++my_count[4];
    }
    else if (c_mag4 < 10000.0*c4_over_n2_ratio_zero*n_mag2) {
      ++my_count[5];
    }
  }

  int8 count[6];
  MPI_Reduce(my_count,count,6,MPI_INT8,MPI_SUM,0,mpi_comm);
  if (mpi_rank == 0) {
    cout << " compact faces flagged for potential grouping: " << count[1] << " out of " << count[0] << " (" << double(count[1])/double(count[0])*100.0 << "%)" << endl;
    cout << " sensitivity check: additional compact faces flagged if c4_over_n2_ratio_zero reduced: " << count[2] << " " << count[3] << " " << count[4] << " " << count[5] << endl;
  }

  // sort c approx 0 faces by area...

  sort(magAndIndexVec.begin(),magAndIndexVec.end());

  // using the sorted vec, bin the areas based off first element and count number in bin.
  // group (A-A_0) < eps_A * A_0.

  vector<pair<double,int> > myNormalKindVec; // first: A_0, second: count for (A-A0) < eps_A * A_0
  for (int ii = 0,nii = magAndIndexVec.size(); ii < nii; ++ii) {
    if ( myNormalKindVec.empty() || ( magAndIndexVec[ii].first - myNormalKindVec.back().first > matched_face_tol*myNormalKindVec.back().first ) ) {
      myNormalKindVec.push_back(pair<double,int>(magAndIndexVec[ii].first,1));
    }
    else {
      myNormalKindVec.back().second += 1; // the count
    }
  }

  /*
     for (int ii = 0,nii = myNormalKindVec.size(); ii < nii; ++ii) {
     if (ii == 0) {
     cout << " > group: " << ii << " area: " << myNormalKindVec[ii].first << " count: " << myNormalKindVec[ii].second << endl;
     }
     else {
     cout << " > group: " << ii << " area: " << myNormalKindVec[ii].first << " count: " << myNormalKindVec[ii].second <<
     " diff with prev: " << (myNormalKindVec[ii].first-myNormalKindVec[ii-1].first)/myNormalKindVec[ii-1].first << endl;
     }
     }
     */

  // now repeat the grouping using the median value.
  // group (A-A_med) < eps_A * A_med. TODO

  int my_ngr = myNormalKindVec.size();
  int ngr;
  //MPI_Reduce(&my_ngr,&ngr,1,MPI_INT,MPI_SUM,0,mpi_comm);
  //if (mpi_rank == 0)
  // cout << "Initial compact face groups: " << ngr << endl;

  // only keep groups greater than a threshold...

  vector<double> magVec;
  int8 my_final_count = 0;
  const int ngr_threshold = getIntParam("GROUP_MEMBER_THRESHOLD", 100)/2; // does it make sense to halve it?
  for (int ii = 0,nii = myNormalKindVec.size(); ii < nii; ++ii) {
    // save groups with counts above 100. Note that at this point we
    // only have the areas grouped. In each area, there could be
    // 12 or 14 directions. So assuming we need atleast a few faces
    // in each area, we threshold at 100...
    //
    // also force the areas to larger than some epsilon ....
    if ( (myNormalKindVec[ii].second > ngr_threshold) && (myNormalKindVec[ii].first > small_face_tol*myNormalKindVec[nii-1].first) ) {
      magVec.push_back(myNormalKindVec[ii].first);
      my_final_count += myNormalKindVec[ii].second;
    }
  }

  int8 final_count;
  MPI_Reduce(&my_final_count,&final_count,1,MPI_INT8,MPI_SUM,0,mpi_comm);
  if (mpi_rank == 0)
    cout << "Final count of grouped compact faces: " << final_count << endl;

  my_ngr = magVec.size();
  MPI_Reduce(&my_ngr,&ngr,1,MPI_INT,MPI_MAX,0,mpi_comm);
  if (mpi_rank == 0)
    cout << "Max number of serial area-based compact face groups: " << ngr << endl;

  double my_max_dp = 0.0;
  int my_max_gr = 0;
  if (my_ngr > 0) {

    // accuracy that a unit-normal dot product must be aligned to be considered the same face...

    const double dp_eps = 1.0E-12;

    vector<double> doubleVec;
    int idx_offset = 0;
    int igr = 0;
    for (int ii = 0,nii = magAndIndexVec.size(); ii < nii; ++ii) {
      while ((igr < int(magVec.size())-1)&&(magVec[igr+1] <= magAndIndexVec[ii].first)) {
        ++igr;
        doubleVec.push_back(0.0); doubleVec.push_back(0.0); doubleVec.push_back(0.0);
        idx_offset = doubleVec.size();
      }
      assert(igr < (1<<21));
      // it is possible that a group included in the magVec just fit, so its min value was (1+1*matched_face_tol)*magVec[igr]. Members
      // in that group may have been up to a full matched_face_tol away, so here we use 2x to include these candidates...
      if ( (magAndIndexVec[ii].first - magVec[igr] >= 0.0) && (magAndIndexVec[ii].first - magVec[igr] <= 2.0*matched_face_tol*magVec[igr]) ) {
        // get the unit normal...
        const int ifa = magAndIndexVec[ii].second;
        assert(group_fa[ifa] == -1);
        const double unit_n[3] = {
          n_fa[ifa][0]/magAndIndexVec[ii].first,
          n_fa[ifa][1]/magAndIndexVec[ii].first,
          n_fa[ifa][2]/magAndIndexVec[ii].first
        };
        // check the dot product against existing vectors in this group...
        int idx, idx_end;
        for (idx = idx_offset, idx_end=doubleVec.size(); idx < idx_end; idx += 3) {
          const double dp = unit_n[0]*doubleVec[idx] + unit_n[1]*doubleVec[idx+1] + unit_n[2]*doubleVec[idx+2];
          if (dp > 1.0-dp_eps) {
            my_max_dp = max(my_max_dp,1.0-dp);
            // this is aligned with this vector...
            assert(idx%3 == 0);
            // the bit rules for the group are:
            // msb: group, then member in group, finally bit0 = sign...
            group_fa[ifa] = (igr<<10)|(((idx-idx_offset)/3)<<1);
            break;
          }
          else if (dp < -1.0+dp_eps) {
            my_max_dp = max(my_max_dp,1.0+dp);
            // this is exactly misaligned with this vector...
            assert(idx%3 == 0);
            group_fa[ifa] = (igr<<10)|(((idx-idx_offset)/3)<<1)|1; // 0-bit is the sign bit
            break;
          }
        }
        if (idx == int(doubleVec.size())) {
          assert(idx%3 == 0);
          assert((idx-idx_offset)/3 < (1<<9));
          my_max_gr = max(my_max_gr,(idx-idx_offset)/3);
          group_fa[ifa] = (igr<<10)|(((idx-idx_offset)/3)<<1);
          doubleVec.push_back(unit_n[0]);
          doubleVec.push_back(unit_n[1]);
          doubleVec.push_back(unit_n[2]);
        }
      }
    }

  }

  /*
     for (int ifa = 0; ifa < nfa; ++ifa) {
     cout << "GX_FA" << group_fa[ifa] << " ";
     FOR_I3 cout << x_fa[ifa][i] << " ";
     FOR_I3 cout << n_fa[ifa][i] << " ";
     cout << endl;
     }
     */

  int max_gr;
  MPI_Reduce(&my_max_gr,&max_gr,1,MPI_INT,MPI_MAX,0,mpi_comm);
  if (mpi_rank == 0)
    cout << "Max number of serial orientation-based compact face sub-groups (expect 2:6 in 2d and 4:12 in 3d): " << max_gr+1 << endl;

  double max_dp;
  MPI_Reduce(&my_max_dp,&max_dp,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
  if (mpi_rank == 0)
    cout << "Actual tol for dp matching for compact faces: " << max_dp << endl;

}

void StaticSolver::groupExtendedFacesOld() const {

  const double c2_over_n2_ratio_zero = 1.0E-24;
  const double matched_face_tol = 1.0E-11;

  int8 my_count[6] = { nef, 0, 0, 0, 0, 0 };

  vector<pair<double,int> > magAndIndexVec;
  for (int ief = 0; ief < nef; ++ief) {
    // check magnitudes...
    const double n_mag2 = DOT_PRODUCT(n_ef[ief],n_ef[ief]);
    const double c_mag2 = DOT_PRODUCT(c_ef[ief],c_ef[ief]);
    // only consider faces where c is zero wrt n_mag. Use the same factor
    // as used when determining zero faces in the vd tightened by an additional
    // few orders of magnitude...
    if (c_mag2 <= c2_over_n2_ratio_zero*n_mag2) {
      magAndIndexVec.push_back(pair<double,int>(sqrt(n_mag2),ief));
      ++my_count[1];
    }
    else if (c_mag2 < 10.0*c2_over_n2_ratio_zero*n_mag2) {
      ++my_count[2];
    }
    else if (c_mag2 < 100.0*c2_over_n2_ratio_zero*n_mag2) {
      ++my_count[3];
    }
    else if (c_mag2 < 1000.0*c2_over_n2_ratio_zero*n_mag2) {
      ++my_count[4];
    }
    else if (c_mag2 < 10000.0*c2_over_n2_ratio_zero*n_mag2) {
      ++my_count[5];
    }
  }

  int8 count[6];
  MPI_Reduce(my_count,count,6,MPI_INT8,MPI_SUM,0,mpi_comm);
  if (mpi_rank == 0) {
    cout << " extended faces flagged for potential grouping: " << count[1] << " out of " << count[0] << " (" << double(count[1])/double(count[0])*100.0 << "%)" << endl;
    cout << " sensitivity check: additional extended faces flagged if c2_over_n2_ratio_zero reduced: " << count[2] << " " << count[3] << " " << count[4] << " " << count[5] << endl;
  }

  sort(magAndIndexVec.begin(),magAndIndexVec.end());

  vector<pair<double,int> > myNormalKindVec;
  for (int ii = 0,nii = magAndIndexVec.size(); ii < nii; ++ii) {
    if ( myNormalKindVec.empty() || ( magAndIndexVec[ii].first - myNormalKindVec.back().first > matched_face_tol*myNormalKindVec.back().first ) ) {
      myNormalKindVec.push_back(pair<double,int>(magAndIndexVec[ii].first,1));
    }
    else {
      myNormalKindVec.back().second += 1; // the count
    }
  }

  // let rank 0 gather everyone's groups and build a common version...

  int my_ngr = myNormalKindVec.size();
  int * ngora = NULL;
  if (mpi_rank == 0) ngora = new int[mpi_size];
  MPI_Gather(&my_ngr,1,MPI_INT,ngora,1,MPI_INT,0,mpi_comm);

  double * my_double_buf = new double[my_ngr];
  int * my_int_buf = new int[my_ngr];
  for (int igr = 0; igr < my_ngr; ++igr) {
    my_double_buf[igr] = myNormalKindVec[igr].first;
    my_int_buf[igr] = myNormalKindVec[igr].second;
  }

  int * disp_full = NULL;
  double * double_buf_full = NULL;
  int * int_buf_full = NULL;
  int full_size;
  if (mpi_rank == 0) {
    // check for overflow...
    int8 full_size_check = 0;
    FOR_RANK full_size_check += ngora[rank];
    assert(full_size_check < TWO_BILLION);
    disp_full = new int[mpi_size];
    disp_full[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      disp_full[rank] = disp_full[rank-1] + ngora[rank-1];
    full_size = disp_full[mpi_size-1] + ngora[mpi_size-1];
    double_buf_full = new double[full_size];
    int_buf_full = new int[full_size];
  }

  MPI_Gatherv(my_double_buf,my_ngr,MPI_DOUBLE,double_buf_full,ngora,disp_full,MPI_DOUBLE,0,mpi_comm);
  delete[] my_double_buf;

  MPI_Gatherv(my_int_buf,my_ngr,MPI_INT,int_buf_full,ngora,disp_full,MPI_INT,0,mpi_comm);
  delete[] my_int_buf; my_int_buf = NULL;

  // now build the final groups and counts on the rank0...

  vector<double> magVec;

  if (mpi_rank == 0) {

    vector<pair<double,int> > sortVec(full_size);
    int8 count_check = 0;
    for (int jj = 0; jj < full_size; ++jj) {
      sortVec[jj].first = double_buf_full[jj];
      sortVec[jj].second = jj;
      count_check += int_buf_full[jj];
    }
    assert(count_check == count[1]);

    sort(sortVec.begin(),sortVec.end());

    vector<pair<double,int8> > normalKindVec;
    for (int ii = 0; ii < full_size; ++ii) {
      const int jj = sortVec[ii].second;
      assert(sortVec[ii].first == double_buf_full[jj]);
      // increase the matched_face_tol here by an order of magnitude to combine groups
      if ( normalKindVec.empty() || ( sortVec[ii].first - normalKindVec.back().first > 10.0*matched_face_tol*normalKindVec.back().first ) ) {
        normalKindVec.push_back(pair<double,int>(double_buf_full[jj],int8(int_buf_full[jj])));
      }
      else {
        normalKindVec.back().second += int8(int_buf_full[jj]);
      }
    }
    delete[] double_buf_full; double_buf_full = NULL;
    delete[] int_buf_full; int_buf_full = NULL;

    cout << "Reduced groups: " << normalKindVec.size() << endl;

    for (int ii = 0,nii = normalKindVec.size(); ii < nii; ++ii) {
      if (ii == 0) {
        cout << " > group: " << ii << " area: " << normalKindVec[ii].first << " count: " << normalKindVec[ii].second << endl;
      }
      else {
        cout << " > group: " << ii << " area: " << normalKindVec[ii].first << " count: " << normalKindVec[ii].second <<
          " diff with prev: " << (normalKindVec[ii].first-normalKindVec[ii-1].first)/normalKindVec[ii-1].first << endl;
      }
    }

    // and reduce...

    int8 final_count = 0;
    const int ngr_threshold = getIntParam("GROUP_MEMBER_THRESHOLD", 100);
    for (int ii = 0,nii = normalKindVec.size(); ii < nii; ++ii) {
      // save groups with counts above 100. Note that at this point we
      // only have the areas grouped. In each area, there could be
      // 12 or 14 directions. So assuming we need atleast a few faces
      // in each area, we threshold at 100...
      if (normalKindVec[ii].second > ngr_threshold) {
        cout << " > FINAL extended face group: " << magVec.size() << " " << normalKindVec[ii].first << " count: " << normalKindVec[ii].second << endl;
        magVec.push_back(normalKindVec[ii].first);
        final_count += normalKindVec[ii].second;
      }
    }

    cout << "final count of grouped extended faces: " << final_count << endl;

    my_ngr = magVec.size();

  }

  // everybody needs the full list of potential groups...

  MPI_Bcast(&my_ngr,1,MPI_INT,0,mpi_comm);

  if (my_ngr > 0) {

    if (mpi_rank != 0) magVec.resize(my_ngr);
    MPI_Bcast(&magVec[0],my_ngr,MPI_DOUBLE,0,mpi_comm);

    // accuracy that a unit-normal dot product must be aligned to be considered the same face...

    const double dp_eps = 1.0E-12;

    vector<double> doubleVec;
    int idx_offset = 0;
    int igr = 0;
    double my_max_dp = 0.0;
    for (int ii = 0,nii = magAndIndexVec.size(); ii < nii; ++ii) {
      while ((igr < int(magVec.size())-1)&&(magVec[igr+1] <= magAndIndexVec[ii].first)) {
        ++igr;
        doubleVec.push_back(0.0); doubleVec.push_back(0.0); doubleVec.push_back(0.0);
        idx_offset = doubleVec.size();
      }
      // it is possible that a group included in the magVec just fit, so its min value was (1+10*matched_face_tol)*magVec[igr]. Members
      // in that group may have been up to a full matched_face_tol away, so here we use 11x to include these candidates...
      if ( (magAndIndexVec[ii].first - magVec[igr] >= 0.0) && (magAndIndexVec[ii].first - magVec[igr] <= 11.0*matched_face_tol*magVec[igr]) ) {
        // get the unit normal...
        const int ief = magAndIndexVec[ii].second;
        assert(group_ef[ief] == -1);
        const double unit_n[3] = {
          n_ef[ief][0]/magAndIndexVec[ii].first,
          n_ef[ief][1]/magAndIndexVec[ii].first,
          n_ef[ief][2]/magAndIndexVec[ii].first
        };
        // check the dot product against existing vectors in this group...
        int idx,idx_end;
        for (idx = idx_offset, idx_end=doubleVec.size(); idx < idx_end; idx += 3) {
          const double dp = unit_n[0]*doubleVec[idx] + unit_n[1]*doubleVec[idx+1] + unit_n[2]*doubleVec[idx+2];
          if (dp > 1.0-dp_eps) {
            my_max_dp = max(my_max_dp,1.0-dp);
            // this is aligned with this vector...
            assert(idx%3 == 0);
            group_ef[ief] = idx/3;
            break;
          }
          else if (dp < -1.0+dp_eps) {
            my_max_dp = max(my_max_dp,1.0+dp);
            // this is exactly misaligned with this vector...
            assert(idx%3 == 0);
            group_ef[ief] = -idx/3-2; // temporarily use -2 indexing
            break;
          }
        }
        if (idx == int(doubleVec.size())) {
          assert(idx%3 == 0);
          assert((idx-idx_offset)/3 < (1<<4));
          group_ef[ief] = idx/3;
          doubleVec.push_back(unit_n[0]);
          doubleVec.push_back(unit_n[1]);
          doubleVec.push_back(unit_n[2]);
        }
      }
    }
    // complete any unfinished groups, and terminate with 0,0,0
    while (igr < int(magVec.size())) {
      ++igr;
      doubleVec.push_back(0.0); doubleVec.push_back(0.0); doubleVec.push_back(0.0);
    }

    double max_dp;
    MPI_Reduce(&my_max_dp,&max_dp,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);

    // now switch the meaning of my_ngr...

    my_ngr = doubleVec.size();
    MPI_Gather(&my_ngr,1,MPI_INT,ngora,1,MPI_INT,0,mpi_comm);

    if (mpi_rank == 0) {
      // check for overflow...
      int8 full_size_check = 0;
      FOR_RANK full_size_check += ngora[rank];
      assert(full_size_check < TWO_BILLION);
      assert(disp_full != NULL);
      disp_full[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        disp_full[rank] = disp_full[rank-1] + ngora[rank-1];
      full_size = disp_full[mpi_size-1] + ngora[mpi_size-1];
      assert(double_buf_full == NULL); double_buf_full = new double[full_size];
      assert(full_size%3 == 0);
      assert(int_buf_full == NULL); int_buf_full = new int[full_size/3];
      for (int ii = 0; ii < full_size/3; ++ii)
        int_buf_full[ii] = -1;
    }

    MPI_Gatherv(&doubleVec[0],my_ngr,MPI_DOUBLE,double_buf_full,ngora,disp_full,MPI_DOUBLE,0,mpi_comm);

    if (mpi_rank == 0) {
      double max_dp_rank0 = 0.0;
      for (int igr = 0, igr_end=magVec.size(); igr < igr_end; ++igr) {
        doubleVec.clear();
        FOR_RANK {
          assert(disp_full[rank]%3 == 0);
          int idx_i = disp_full[rank]/3;
          double unit_n[3];
          FOR_I3 unit_n[i] = double_buf_full[disp_full[rank]++];
          while (!((unit_n[0] == 0.0)&&(unit_n[1] == 0.0)&&(unit_n[2] == 0.0))) {
            int idx,idx_end;
            for (idx = 0, idx_end=doubleVec.size(); idx < idx_end; idx += 3) {
              const double dp = unit_n[0]*doubleVec[idx] + unit_n[1]*doubleVec[idx+1] + unit_n[2]*doubleVec[idx+2];
              if (dp > 1.0-dp_eps) {
                // this is aligned with this vector...
                // the bit rules for the group are:
                // msb: group, then member in group, finally bit0 = sign...
                int_buf_full[idx_i] = (igr<<5)|((idx/3)<<1);
                max_dp_rank0 = max(max_dp_rank0,1.0-dp);
                break;
              }
              else if (dp < -1.0+dp_eps) {
                // this is exactly misaligned with this vector...
                int_buf_full[idx_i] = (igr<<5)|((idx/3)<<1)|1; // 0-bit is the sign bit
                max_dp_rank0 = max(max_dp_rank0,1.0+dp);
                break;
              }
            }
            if (idx == int(doubleVec.size())) {
              assert(idx/3 < (1<<4));
              int_buf_full[idx_i] = (igr<<5)|((idx/3)<<1);
              doubleVec.push_back(unit_n[0]);
              doubleVec.push_back(unit_n[1]);
              doubleVec.push_back(unit_n[2]);
            }
            // grab the next unit normal...
            assert(disp_full[rank]%3 == 0);
            idx_i = disp_full[rank]/3;
            FOR_I3 unit_n[i] = double_buf_full[disp_full[rank]++];
          }
        }
        cout << "group: " << igr << " doubleVec.size()/3: " << doubleVec.size()/3 << endl;
      }

      cout << "actual tol for dp matching: " << max_dp << " " << max_dp_rank0 << endl;

      // check and reset disp for return int...
      int disp_check = 0;
      FOR_RANK {
        disp_check += ngora[rank];
        assert(disp_full[rank] == disp_check);
        assert(ngora[rank]%3 == 0);
        ngora[rank] /= 3;
        if (rank == 0) {
          disp_full[rank] = 0;
        }
        else {
          disp_full[rank] = disp_full[rank-1] + ngora[rank-1];
        }
      }

      // now disp_full and ngora are set to scatterv the int...

    }

    assert(my_ngr%3 == 0);
    my_ngr /= 3;

    assert(my_int_buf == NULL);
    my_int_buf = new int[my_ngr];

    MPI_Scatterv(int_buf_full,ngora,disp_full,MPI_INT,my_int_buf,my_ngr,MPI_INT,0,mpi_comm);

    // now on each processor, modify the index...

    for (int ii = 0,nii = magAndIndexVec.size(); ii < nii; ++ii) {
      const int ief = magAndIndexVec[ii].second;
      if (group_ef[ief] >= 0) {
        const int idx = group_ef[ief];
        assert((idx >= 0)&&(idx < my_ngr));
        assert(my_int_buf[idx] >= 0);
        group_ef[ief] = my_int_buf[idx];
      }
      else if (group_ef[ief] <= -2) {
        const int idx = -group_ef[ief]-2;
        assert((idx >= 0)&&(idx < my_ngr));
        assert(my_int_buf[idx] >= 0);
        group_ef[ief] = my_int_buf[idx];
        if (my_int_buf[idx] & 1)
          group_ef[ief] = (my_int_buf[idx] & (~1));
        else
          group_ef[ief] = (my_int_buf[idx] | 1);
      }
    }

    delete[] my_int_buf;

    /*
       int8 * gr_count = NULL;
       if (mpi_rank == 0) gr_count = new int8[my_ngr];
       MPI_Reduce(my_gr_count,gr_count,my_ngr,MPI_INT8,MPI_SUM,0,mpi_comm);
       if (mpi_rank == 0) {
       int8 count_check = 0;
       for (int igr = 0; igr < my_ngr; ++igr) {
       count_check += gr_count[igr];
       cout << "igr: " << igr << " count: " << gr_count[igr] << endl;
       }
       cout << " should match final count of grouped faces: " << count_check << endl;
       delete[] gr_count;
       }
       delete[] my_gr_count;
       */

  }

  if (mpi_rank == 0) {
    delete[] ngora;
    delete[] disp_full;
    delete[] double_buf_full;
    delete[] int_buf_full;
  }

  /*
     FOR_INTERPROC_IEF {
     const int icv0 = cvoef[ief][0]; assert((icv0 >= 0)&&(icv0 < ncv));
     const int icv1 = cvoef[ief][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g2));
     int rank,bits,index;
     if (icv1 < ncv_g) {
     BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv1-ncv]);
     }
     else {
     assert(icv1 < ncv_g2);
     BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g2[icv1-ncv_g]);
     }
// flip???  orientation problems in FlowSolver.hpp
if (mpi_rank < rank) {
if (group_ef[ief] & 1) group_ef[ief] = (group_ef[ief] & (~1));
else                   group_ef[ief] = (group_ef[ief] | 1);
}
}
*/

/*
// resort the extended faces...

vector< pair<int,int> > efVec(nef);
for (int ief = 0; ief < nef; ++ief) {
const int igr = group_ef[ief];
efVec[ief] = pair<int,int>(igr,ief);
}

// sort...

sort(efVec.begin(),efVec.end());

int * group_ef_old = group_ef;
int (*cvoef_old)[2] = cvoef;
double (*n_ef_old)[3] = n_ef;
double (*c_ef_old)[3] = c_ef;
delete[] group_ef; group_ef = new int[nef];
delete[] cvoef; cvoef = new int[nef][2];
delete[] n_ef; n_ef = new double[nef][3];
delete[] c_ef; c_ef = new double[nef][3];
for (int ief = 0; ief > nef; ++ief) {
const int ief_old = efVec[ief].second;
group_ef[ief] = group_ef_old[ief_old];
FOR_I2 cvoef[ief][i] = cvoef_old[ief_old][i];
FOR_I3 n_ef[ief][i] = n_ef_old[ief_old][i];
FOR_I3 c_ef[ief][i] = c_ef_old[ief_old][i];
}

efVec.clear();
*/

}

void StaticSolver::calcCvGrad(double (*__restrict__ dpdx)[3],const double *__restrict__ p) {

  for (int icv = 0; icv < ncv; ++icv) {

    const int coc_f = cvocv_i[icv];
    for (int i = 0; i < 3; ++i) {
      dpdx[icv][i] = 0.0; //cvocv_grad_coeff[coc_f][i] * p[icv];
    }

    for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      for (int i =0; i < 3; ++i) {
        dpdx[icv][i] += cvocv_grad_coeff[coc][i] * (p[icv_nbr] - p[icv]);
      }
    }

  }

}


void StaticSolver::calcCvGrad(double (*__restrict__ dudx)[3][3],const double (*__restrict__ u)[3]) {

  for (int icv = 0; icv < ncv; ++icv) {

    const int coc_f = cvocv_i[icv];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        dudx[icv][i][j] = 0.0; //cvocv_grad_coeff[coc_f][j]* u[icv][i];
      }
    }

    for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      for (int i =0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          dudx[icv][i][j] += cvocv_grad_coeff[coc][j] * (u[icv_nbr][i] - u[icv][i]);
        }
      }
    }

  }

}

void StaticSolver::_calcCvGradRange(double (*__restrict__ dudx)[3][3],const double (*__restrict__ u)[3],
    const int icv_start, const int icv_end) {

  for (int icv = icv_start; icv < icv_end; ++icv) {

    const int coc_f = cvocv_i[icv];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        dudx[icv][i][j] = 0.0; //cvocv_grad_coeff[coc_f][j]* u[icv][i];
      }
    }

    for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      for (int i =0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          dudx[icv][i][j] += cvocv_grad_coeff[coc][j] * (u[icv_nbr][i] - u[icv][i]);
        }
      }
    }

  }
}

void StaticSolver::checkCvGradCoeff(int * cv_flag) {

  COUT1("StaticSolver::checkCvGradCoeff()");

  const double grad_check[3] = { 1.1234, -1.3243, 1.5321 }; // some order-1 gradient
  double * phi               = new double[ncv_g];
  double (*grad_phi)[3]      = new double[ncv][3];

  for (int icv = 0; icv < ncv_g; ++icv)
    phi[icv] = DOT_PRODUCT(x_cv[icv],grad_check);

  calcCvGrad(grad_phi,phi);

  for (int icv = 0; icv < ncv; ++icv) {
    FOR_I3 {
      grad_phi[icv][i] /= grad_check[i];
      grad_phi[icv][i] -= 1.0;
    }
  }

  MiscUtils::dumpRange(grad_phi,ncv,"cv grad coeff check (should be order one)");

  if ( cv_flag) {

    for (int icv = 0; icv < ncv; ++icv) {
      if ( cv_flag[icv] != 0 ) {
        for (int i =0; i < 3; ++i)
          grad_phi[icv][i] = 0.0;
      }
    }

    MiscUtils::dumpRange(grad_phi,ncv,"cv grad coeff check -- internal only (should be zero)");

  }

  delete[] grad_phi;
  delete[] phi;
}

void StaticSolver::calcCvGradCoeff() {

  COUT1("StaticSolver::calcCvGradCoeff()");

  // this routine solves for the gradient using
  // vol * grad(phi) = sum( phi_face * A_face )
  // where phi_face = 0.5*(phi + phi_nbr) + grad(phi)*(x_fa - 0.5*(x_cv + x_cv_nbr))
  // this results in a system as follows:
  //
  // [A]{grad_phi} = sum(0.5*(phi + phi_nbr)*A_face)
  //
  // which we can solve with A^{-1} as
  //
  // grad_phi = A^{-1}*sum(0.5*(phi + phi_nbr)*A_face)
  //
  // the main challenge with this approach is the case of limited data where A
  // can become singular,
  // which we must treat stably by reducing the gradient, ultimately to
  // zero. To handle this, we treat boundary faces in a second pass
  // that considers the absolute projection of other areas on the face...

  assert(cvocv_i != NULL);
  assert(cvocv_v != NULL);
  assert(bfocv_i != NULL);
  assert(bfocv_v != NULL);
  assert(cvocv_grad_coeff == NULL);
  cvocv_grad_coeff = new double[cvocv_i[ncv]][3];

  assert(Gij_bf);

  const double sigma_epsilon = getDoubleParam("SIGMA_THRESH", 4.5e-1);

  int * cv_flag  = new int[ncv];

  FOR_ICV {

    const int coc_f = cvocv_i[icv];
    const int coc_l = cvocv_i[icv+1]-1;
    for (int coc = coc_f; coc <= coc_l; ++coc)
      FOR_I3 cvocv_grad_coeff[coc][i] = 0.0;

    double A[9];
    for (int ii = 0; ii < 9; ++ii)
      A[ii] = 0.0;

    // start with volume on the diagonal
    for (int i =0; i < 3; ++i)
      A[(i)*3+i] = vol_cv[icv];

    double gcl[3] = { 0.0, 0.0, 0.0 }; // for checking

    for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
      const int ibf = bfocv_v[boc];
      for (int j =0; j < 3; ++j) {
        for (int i =0; i < 3; ++i) {
          A[(j)*3+i] += Gij_bf[ibf][i][j];
        }
      }
      FOR_I3 gcl[i] += n_bf[ibf][i];
    }

    for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
      const int ifa = faocv_v[foc];
      double fa_sign;
      int icv_nbr;
      if ( cvofa[ifa][0] == icv ) {
        fa_sign = 1.0;
        icv_nbr = cvofa[ifa][1];
      }
      else {
        assert( cvofa[ifa][1] == icv );
        fa_sign = -1.0;
        icv_nbr = cvofa[ifa][0];
      }
      // also get the nbr associated with this icv_nbr...
      const int coc = coc_f+foc-faocv_i[icv]+1;
      assert(cvocv_v[coc] == icv_nbr);
      double r[3]; FOR_I3 r[i] = x_fa[ifa][i] - 0.5*(x_cv[icv][i] + x_cv[icv_nbr][i]);

      for (int j =0; j < 3; ++j) {
        for (int i =0; i < 3; ++i) {
          A[(j)*3+i] -= r[j]*n_fa[ifa][i]*fa_sign;
        }
      }

      FOR_I3 cvocv_grad_coeff[coc][i]   += 0.5*n_fa[ifa][i]*fa_sign;
      FOR_I3 cvocv_grad_coeff[coc_f][i] -= 0.5*n_fa[ifa][i]*fa_sign;
      FOR_I3 gcl[i] += n_fa[ifa][i]*fa_sign;
    }

    // gcl contains a delta(area), so its square should be compared to a vol^(4/3)
    if (DOT_PRODUCT(gcl,gcl) > 1.0E-12*pow(vol_cv[icv],4.0/3.0))
      cout << "Warning: got large gcl: " << COUT_VEC(gcl) << " at location: " << COUT_VEC(x_cv[icv]) << endl;

    // rescale the matrices with the inv_vol to help the conditioning of the system...
    const double inv_vol_cv = 1.0/vol_cv[icv];
    for (int j = 0; j < 3; ++j) {
      for (int i =0; i < 3; ++i) {
        A[(j)*3+i] *= inv_vol_cv;
      }
    }

    for (int coc = coc_f; coc <= coc_l; ++coc)
      for (int i =0; i < 3; ++i)
        cvocv_grad_coeff[coc][i] *= inv_vol_cv;

    double U[9], V[9], sigma[3];
    double svd_eps_tol;
    calcSvd33(U,V,sigma,A,svd_eps_tol);

    // the singular values are not necessarily sorted from the prior routine..
    // so we need to loop and build the thresholding here...

    cv_flag[icv] = 0; // assume that the grad is ok unless found otherwise..

    double sigma_inv[3];
    for (int i =0; i < 3; ++i) {
      assert( sigma[i] >= 0.0);
      if ( sigma[i] < sigma_epsilon) {
        sigma_inv[i] = 0.0; // cant trust the gradient in this dir..
        cv_flag[icv] = 1;
      } else {
        sigma_inv[i] = 1.0/sigma[i];
      }
    }

    // now supply the regularized inverse .. let V <--- V\Sigma^+

    for (int j = 0; j < 3; ++j) {
      for (int i =0; i < 3; ++i) {
        V[(j)*3+i] *= sigma_inv[j];
      }
    }

    // finally supply the inverse for the gradient coefficients ..

    for (int coc = coc_f; coc <= coc_l; ++coc) {

      double b_tmp[3];
      for (int i =0; i < 3; ++i) {
        b_tmp[i] = 0.0;
        for (int j =0; j < 3; ++j) {
          b_tmp[i] += A[(j)*3+i]*cvocv_grad_coeff[coc][j];
        }
      }


      double tmp[3];
      // V*U^T[r]
      for (int i =0; i <3 ; ++i)
        tmp[i] = U[(i)*3+0]*(cvocv_grad_coeff[coc][0] - b_tmp[0]) +
          U[(i)*3+1]*(cvocv_grad_coeff[coc][1] - b_tmp[1]) +
          U[(i)*3+2]*(cvocv_grad_coeff[coc][2] - b_tmp[2]) ;

      for (int i =0; i < 3; ++i)
        cvocv_grad_coeff[coc][i] += V[0*3+i]*tmp[0] +
          V[1*3+i]*tmp[1] +
          V[2*3+i]*tmp[2];
    }
  }

  delete[] Gij_bf; Gij_bf = NULL; // not needed for anything else but this gradient

  checkCvGradCoeff(cv_flag);

  delete[] cv_flag;

}

void StaticSolver::checkGradExtended() {

  COUT1("StaticSolver::checkGradExtended()");

  const double grad_check[3] = { 1.123, 2.134, -1.3423 };
  //const double grad_check[3] = { 1.1234, -1.3243, 1.5321 }; // some order-1 gradient

  // x_cv is available in ncv_g2...
  double * phi = new double[ncv_g2];
  for (int icv = 0; icv < ncv_g2; ++icv) phi[icv] = DOT_PRODUCT(x_cv[icv],grad_check);

  double (*grad_phi)[3] = new double[ncv][3];
  FOR_ICV FOR_I3 grad_phi[icv][i] = 0.0;

  /*
     double (*x_fa)[3] = new double[nfa][3];
     for (int ifa = 0; ifa < nfa; ++ifa) {
     const int ino = noofa_v[noofa_i[ifa]];
     int ino1 = noofa_v[noofa_i[ifa+1]-1];
     double area_fa = 0.0;
     for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
     const int ino0 = ino1;
     ino1 = noofa_v[nof];
     if ((ino != ino0)&&(ino != ino1)) {
     const double this_n[3] = TRI_NORMAL_2(x_no[ino],x_no[ino0],x_no[ino1]);
     const double n_mag  = MAG(this_n);
     area_fa += n_mag;
     FOR_I3 x_fa[ifa][i] += n_mag*(x_no[ino][i] + x_no[ino0][i] + x_no[ino1][i]);
     }
     }
     FOR_I3 x_fa[ifa][i] /= 3.0*area_fa;
     }

     FOR_INTERNAL_IFA {
     const int icv0 = cvofa[ifa][0];
     const int icv1 = cvofa[ifa][1];
     const double phi_fa = DOT_PRODUCT(x_fa[ifa],grad_check);
     double flux[3]; FOR_I3 flux[i] = phi_fa*n_fa[ifa][i];
     FOR_I3 grad_phi[icv0][i] += flux[i];
     FOR_I3 grad_phi[icv1][i] -= flux[i];
     }
     FOR_INTERPROC_IFA {
     const int icv0 = cvofa[ifa][0];
     const double phi_fa = DOT_PRODUCT(x_fa[ifa],grad_check);
     double flux[3]; FOR_I3 flux[i] = phi_fa*n_fa[ifa][i];
     FOR_I3 grad_phi[icv0][i] += flux[i];
     }
     delete[] x_fa;
     */

  FOR_INTERNAL_IEF {
    const int icv0 = cvoef[ief][0]; assert((icv0 >= 0)&&(icv0 < ncv));
    const int icv1 = cvoef[ief][1]; assert((icv1 >= 0)&&(icv1 < ncv));
    double flux[3];	FOR_I3 flux[i] = 0.5*n_ef[ief][i]*(phi[icv0]+phi[icv1]) + 0.5*c_ef[ief][i]*(phi[icv1]-phi[icv0]);
    FOR_I3 grad_phi[icv0][i] += flux[i];
    FOR_I3 grad_phi[icv1][i] -= flux[i];
  }

  FOR_INTERPROC_IEF {
    const int icv0 = cvoef[ief][0]; assert((icv0 >= 0)&&(icv0 < ncv));
    const int icv1 = cvoef[ief][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g2));
    double flux[3];	FOR_I3 flux[i] = 0.5*n_ef[ief][i]*(phi[icv0]+phi[icv1]) + 0.5*c_ef[ief][i]*(phi[icv1]-phi[icv0]);
    FOR_I3 grad_phi[icv0][i] += flux[i];
  }

  // and close the gradient at the boundary...

  FOR_IBF {
    const int icv = cvobf[ibf];
    const double phi_bf = DOT_PRODUCT(x_bf[ibf],grad_check);
    FOR_I3 grad_phi[icv][i] += phi_bf*n_bf[ibf][i];
  }

  FOR_ICV {
    FOR_I3 grad_phi[icv][i] /= vol_cv[icv]*grad_check[i]; // to produce 1
    FOR_I3 grad_phi[icv][i] -= 1.0; // to produce an error
  }

  MiscUtils::dumpRange(grad_phi,ncv,"extended grad (should be order one)");

  // check the internal only contribution.  note that internal here means
  // that your nbrs cannot touch boundaries either...

  int * cv_flag = new int[ncv_g2];

  for (int icv = 0; icv < ncv_g2; ++icv)
    cv_flag[icv] = 0;

  for (int ibf = 0; ibf < nbf; ++ibf)
    cv_flag[cvobf[ibf]] = 1;

  updateCv2Data(cv_flag);

  int8 ncv_internal = 0;

  // cvocv_{i,v} have not been built yet... we'll have to go through cvofa.
  assert( cvofa != NULL);

  int * cv_flag2 = new int[ncv];
  for (int icv = 0; icv < ncv; ++icv)
    cv_flag2[icv] = cv_flag[icv];

  for (int ifa = 0; ifa < nfa_i; ++ifa) {
    const int icv0 = cvofa[ifa][0]; assert( (icv0 >= 0)&&(icv0 < ncv));
    const int icv1 = cvofa[ifa][1]; assert( (icv1 >= 0)&&(icv1 < ncv));

    cv_flag2[icv0] = max(cv_flag2[icv0], cv_flag[icv1]);
    cv_flag2[icv1] = max(cv_flag2[icv1], cv_flag[icv0]);
  }

  for (int ifa = nfa_i; ifa < nfa; ++ifa) {
    const int icv0 = cvofa[ifa][0]; assert( (icv0 >= 0)&&(icv0 < ncv));
    const int icv1 = cvofa[ifa][1]; assert( (icv1 >= ncv) && (icv1 < ncv_g));

    cv_flag2[icv0] = max(cv_flag2[icv0], cv_flag[icv1]);
  }

  for (int icv = 0; icv < ncv; ++icv) {
    if ( cv_flag2[icv] == 0 ) {
      ++ncv_internal;
    }
    else {
      FOR_I3 grad_phi[icv][i] = 0.0;
    }
  }

  MiscUtils::dumpRange(grad_phi,ncv,"extended grad -- internal only (should be zero)");

  int8 ncv_internal_global;
  MPI_Reduce(&ncv_internal,&ncv_internal_global,1,MPI_INT8,MPI_SUM,0,mpi_comm);
  if ( mpi_rank == 0 )
    cout << " > internal cv fraction : " << double(ncv_internal_global)/double(ncv_global) << endl;

  delete[] cv_flag;
  delete[] cv_flag2;
  delete[] grad_phi;
  delete[] phi;

}

void StaticSolver::buildCvPrcomm() {

  COUT1("StaticSolver::buildCvPrcomm()");

  // ------------------------------------------------
  // build the paired communicator -- this should be
  // symmetric in terms of rank, but possibly not count...
  // ------------------------------------------------

  assert(cvPrcommVec.empty());

  int rank_current = -1;
  int bits_current = -1;
  CvPrcomm * prcomm = NULL;
  for (int icv = ncv;  icv < ncv_g; ++icv) {
    int rank,bits,index;
    BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv-ncv]);
    if (rank > rank_current) {
      // we just switched ranks. If there was a prcomm, then
      // complete the size...
      if (prcomm)
        prcomm->unpack_size = icv-prcomm->unpack_offset;
      rank_current = rank;
      bits_current = -1;
      cvPrcommVec.push_back(CvPrcomm());
      prcomm = &cvPrcommVec.back();
      prcomm->rank = rank;
      prcomm->unpack_offset = icv;
    }
    else {
      assert(rank_current == rank);
      assert(prcomm);
    }
    // just tests monotoncity of the bits...
    if (bits > bits_current) {
      bits_current = bits;
    }
    else {
      assert(bits_current == bits);
    }
  }
  // we just finished. If there was a prcomm, then
  // complete the size...
  if (prcomm) prcomm->unpack_size = ncv_g-prcomm->unpack_offset;

  // finally, we need to send/recv the indices and bits to the pack side
  // and build the packVecs. Note that these should be symmetric by construction
  // because internal faces always have 2 cvs, but counts may not be the
  // same because more than one face can be on a single cv...

  MPI_Request * sendRequestArray = new MPI_Request[cvPrcommVec.size()];
  MPI_Request * recvRequestArray = new MPI_Request[cvPrcommVec.size()];
  for (int ii = 0, ii_end=cvPrcommVec.size(); ii < ii_end; ++ii) {

    // post irecv...
    MPI_Irecv(&(cvPrcommVec[ii].pack_size),1,MPI_INT,cvPrcommVec[ii].rank,12345,mpi_comm,&(recvRequestArray[ii]));

    // and the send...
    MPI_Issend(&(cvPrcommVec[ii].unpack_size),1,MPI_INT,cvPrcommVec[ii].rank,12345,mpi_comm,&(sendRequestArray[ii]));

  }

  MPI_Waitall(cvPrcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
  MPI_Waitall(cvPrcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

  // now send from the unpack side to the pack side...

  int pack_size = 0;
  for (int ii = 0, ii_end=cvPrcommVec.size(); ii < ii_end; ++ii)
    pack_size += cvPrcommVec[ii].pack_size;

  uint8 * recv_rbi_g = new uint8[pack_size];

  pack_size = 0;
  for (int ii = 0, ii_end=cvPrcommVec.size(); ii < ii_end; ++ii) {

    // post irecv...
    MPI_Irecv(recv_rbi_g+pack_size,cvPrcommVec[ii].pack_size,MPI_UINT8,
        cvPrcommVec[ii].rank,12346,mpi_comm,&(recvRequestArray[ii]));
    pack_size += cvPrcommVec[ii].pack_size;

    // and the send...
    MPI_Issend(rbi_g+(cvPrcommVec[ii].unpack_offset-ncv),cvPrcommVec[ii].unpack_size,MPI_UINT8,
        cvPrcommVec[ii].rank,12346,mpi_comm,&(sendRequestArray[ii]));

  }

  MPI_Waitall(cvPrcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
  MPI_Waitall(cvPrcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

  delete[] sendRequestArray;
  delete[] recvRequestArray;

  // now build the packVec and periodicity...

  double R[9], t[3];
  pack_size = 0;
  for (int ii = 0, ii_end=cvPrcommVec.size(); ii < ii_end; ++ii) {
    assert(cvPrcommVec[ii].packVec.empty());
    cvPrcommVec[ii].packVec.resize(cvPrcommVec[ii].pack_size);
    CvPrcomm::Transform * transform = NULL;
    int bits_current = -1;
    for (int i = 0; i < cvPrcommVec[ii].pack_size; ++i) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,recv_rbi_g[pack_size+i]);
      assert(rank == mpi_rank);
      if (bits > bits_current) {
        // bits are about to change, so complete any translate and rotate...
        if (transform) transform->end = i;
        // look for new rotations/translations...
        const bool has_R = PeriodicData::getPeriodicR(R,bits);
        const bool has_t = PeriodicData::getPeriodicT(t,bits);
        cvPrcommVec[ii].transformVec.push_back(CvPrcomm::Transform(has_R,R,has_t,t,bits,i));
        transform = &(cvPrcommVec[ii].transformVec.back());
        bits_current = bits;
      }
      else {
        assert(bits_current == bits);
      }
      assert((index >= 0)&&(index < ncv));
      cvPrcommVec[ii].packVec[i] = index;
    }
    if (transform) transform->end = cvPrcommVec[ii].pack_size;
    pack_size += cvPrcommVec[ii].pack_size;
  }

  delete[] recv_rbi_g;

}

void StaticSolver::buildCv2Prcomm() {

  COUT1("StaticSolver::buildCv2Prcomm()");

  assert(rbi_g2);

  // ------------------------------------------------
  // build the paired communicator for this level 2 ghost
  // data. See discussion below. This can be non-symmetric
  // unless built carefully.
  // ------------------------------------------------

  // building the cv2PrcommVec is a bit more difficult than the cvPrcommVec
  // because it is possible for the cv2PrcommVec on its own to be non-symmetric.
  // This is because ghost parts of the pair-wise data can already be in
  // the ncv_g range for other reasons. So the cv2PrcommVec needs to have AT LEAST
  // the rank coverage of the cvPrcomm, and potentially more...

  int ii_cvPrcomm = 0;
  assert(cv2PrcommVec.empty());

  int rank_current = -1;
  int bits_current = -1;
  CvPrcomm * prcomm = NULL;
  for (int icv = ncv_g;  icv < ncv_g2; ++icv) {
    int rank,bits,index;
    BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g2[icv-ncv_g]);
    if (rank > rank_current) {
      // we just switched ranks. If there was a prcomm, then
      // complete the size...
      if (prcomm)
        prcomm->unpack_size = icv-prcomm->unpack_offset;
      while ((ii_cvPrcomm < int(cvPrcommVec.size()))&&(cvPrcommVec[ii_cvPrcomm].getRank() <= rank)) {
        if (cvPrcommVec[ii_cvPrcomm].getRank() < rank) {
          // add an empty CvPrcomm...
          cv2PrcommVec.push_back(CvPrcomm());
          cv2PrcommVec.back().rank = cvPrcommVec[ii_cvPrcomm].getRank();
          cv2PrcommVec.back().unpack_offset = icv;
          cv2PrcommVec.back().unpack_size = 0;
        }
        ++ii_cvPrcomm;
      }
      rank_current = rank;
      bits_current = -1;
      cv2PrcommVec.push_back(CvPrcomm());
      prcomm = &cv2PrcommVec.back();
      prcomm->rank = rank;
      prcomm->unpack_offset = icv;
    }
    else {
      assert(rank_current == rank);
      assert(prcomm);
    }
    if (bits > bits_current) {
      bits_current = bits;
    }
    else {
      assert(bits_current == bits);
    }
  }
  // we just finished. If there was a prcomm, then
  // complete the size...
  if (prcomm)
    prcomm->unpack_size = ncv_g2-prcomm->unpack_offset;
  // add any not represented in cvPrcomm at the end...
  while (ii_cvPrcomm < int(cvPrcommVec.size())) {
    // add an empty CvPrcomm...
    cv2PrcommVec.push_back(CvPrcomm());
    cv2PrcommVec.back().rank = cvPrcommVec[ii_cvPrcomm].getRank();
    cv2PrcommVec.back().unpack_offset = ncv_g2;
    cv2PrcommVec.back().unpack_size = 0;
    ++ii_cvPrcomm;
  }
  assert(ii_cvPrcomm == int(cvPrcommVec.size()));

  // take a look at pairing...
  /*
     FOR_RANK {
     if (mpi_rank == rank) {
     cout << "rank: " << mpi_rank << endl;
     for (int ii = 0; ii < cv2PrcommVec.size(); ++ii) {
     cout << " > cv2PrcommVec[ii].rank: " << cv2PrcommVec[ii].rank << endl;
     }
     }
     MPI_Pause("press to return see next rank");
     }
     MPI_Pause("OK");
     */

  // finally, we need to send/recv the indices and bits to the pack side
  // and build the packVecs. Note that these should be symmetric by construction
  // for reasons discussed above. We essentially forced it to be true...

  MPI_Request * sendRequestArray = new MPI_Request[cv2PrcommVec.size()];
  MPI_Request * recvRequestArray = new MPI_Request[cv2PrcommVec.size()];
  for (int ii = 0, ii_end=cv2PrcommVec.size(); ii < ii_end; ++ii) {

    // post irecv...
    MPI_Irecv(&(cv2PrcommVec[ii].pack_size),1,MPI_INT,cv2PrcommVec[ii].rank,21345,mpi_comm,&(recvRequestArray[ii]));

    // and the send...
    MPI_Issend(&(cv2PrcommVec[ii].unpack_size),1,MPI_INT,cv2PrcommVec[ii].rank,21345,mpi_comm,&(sendRequestArray[ii]));

  }

  MPI_Waitall(cv2PrcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
  MPI_Waitall(cv2PrcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

  // now send from the unpack side to the pack side...

  int pack_size = 0;
  for (int ii = 0, ii_end=cv2PrcommVec.size(); ii < ii_end; ++ii)
    pack_size += cv2PrcommVec[ii].pack_size;

  uint8 * recv_rbi_g2 = new uint8[pack_size];
  pack_size = 0;
  for (int ii = 0, ii_end=cv2PrcommVec.size(); ii < ii_end; ++ii) {

    // post irecv...
    MPI_Irecv(recv_rbi_g2+pack_size,cv2PrcommVec[ii].pack_size,MPI_UINT8,
        cv2PrcommVec[ii].rank,21346,mpi_comm,&(recvRequestArray[ii]));
    pack_size += cv2PrcommVec[ii].pack_size;

    // and the send...
    MPI_Issend(rbi_g2+(cv2PrcommVec[ii].unpack_offset-ncv_g),cv2PrcommVec[ii].unpack_size,MPI_UINT8,
        cv2PrcommVec[ii].rank,21346,mpi_comm,&(sendRequestArray[ii]));

  }

  MPI_Waitall(cv2PrcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
  MPI_Waitall(cv2PrcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

  delete[] sendRequestArray;
  delete[] recvRequestArray;

  // now build the packVec and periodicity...

  double R[9],t[3];
  pack_size = 0;
  for (int ii = 0, ii_end=cv2PrcommVec.size(); ii < ii_end; ++ii) {
    assert(cv2PrcommVec[ii].packVec.empty());
    cv2PrcommVec[ii].packVec.resize(cv2PrcommVec[ii].pack_size);
    CvPrcomm::Transform * transform = NULL;
    int bits_current = -1;
    for (int i = 0; i < cv2PrcommVec[ii].pack_size; ++i) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,recv_rbi_g2[pack_size+i]);
      assert(rank == mpi_rank);
      if (bits > bits_current) {
        // bits are about to change, so complete any translate and rotate...
        if (transform) transform->end = i;
        // look for new rotations/translations...
        const bool has_R = PeriodicData::getPeriodicR(R,bits);
        const bool has_t = PeriodicData::getPeriodicT(t,bits);
        cv2PrcommVec[ii].transformVec.push_back(CvPrcomm::Transform(has_R,R,has_t,t,bits,i));
        transform = &(cv2PrcommVec[ii].transformVec.back());
        bits_current = bits;
      }
      else {
        assert(bits_current == bits);
      }
      assert((index >= 0)&&(index < ncv));
      cv2PrcommVec[ii].packVec[i] = index;
    }
    if (transform) transform->end = cv2PrcommVec[ii].pack_size;
    pack_size += cv2PrcommVec[ii].pack_size;
  }

  delete[] recv_rbi_g2;

}

void StaticSolver::buildPrcommSymmetric(vector<CvPrcomm>& prcommVec,uint8 * rbi_g,const int ncv,const int ncv_g) {

  COUT1("buildPrcommSymmetric()");

  // ------------------------------------------------
  // build the paired communicator -- this should be
  // symmetric in terms of rank, but possibly not count...
  // ------------------------------------------------

  assert(prcommVec.empty());

  int rank_current = -1;
  int bits_current = -1;
  CvPrcomm * prcomm = NULL;
  for (int icv = ncv; icv < ncv_g; ++icv) {
    int rank,bits,index;
    BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv-ncv]);
    if (rank > rank_current) {
      // we just switched ranks. If there was a prcomm, then
      // complete the size...
      if (prcomm)
        prcomm->unpack_size = icv-prcomm->unpack_offset;
      rank_current = rank;
      bits_current = -1;
      prcommVec.push_back(CvPrcomm());
      prcomm = &prcommVec.back();
      prcomm->rank = rank;
      prcomm->unpack_offset = icv;
    }
    else {
      assert(rank_current == rank);
      assert(prcomm);
    }
    // just tests monotoncity of the bits...
    if (bits > bits_current) {
      bits_current = bits;
    }
    else {
      assert(bits_current == bits);
    }
  }
  // we just finished. If there was a prcomm, then
  // complete the size...
  if (prcomm) prcomm->unpack_size = ncv_g-prcomm->unpack_offset;

  // finally, we need to send/recv the indices and bits to the pack side
  // and build the packVecs. Note that these should be symmetric by construction
  // because internal faces always have 2 cvs, but counts may not be the
  // same because more than one face can be on a single cv...

  MPI_Request * sendRequestArray = new MPI_Request[prcommVec.size()];
  MPI_Request * recvRequestArray = new MPI_Request[prcommVec.size()];
  for (int ii = 0, ii_end=prcommVec.size(); ii < ii_end; ++ii) {

    // post irecv...
    MPI_Irecv(&(prcommVec[ii].pack_size),1,MPI_INT,prcommVec[ii].rank,12345,mpi_comm,&(recvRequestArray[ii]));

    // and the send...
    MPI_Issend(&(prcommVec[ii].unpack_size),1,MPI_INT,prcommVec[ii].rank,12345,mpi_comm,&(sendRequestArray[ii]));

  }

  MPI_Waitall(prcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
  MPI_Waitall(prcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

  // now send from the unpack side to the pack side...

  int pack_size = 0;
  for (int ii = 0, ii_end=prcommVec.size(); ii < ii_end; ++ii)
    pack_size += prcommVec[ii].pack_size;

  uint8 * recv_rbi_g = new uint8[pack_size];

  pack_size = 0;
  for (int ii = 0, ii_end=prcommVec.size(); ii < ii_end; ++ii) {

    // post irecv...
    MPI_Irecv(recv_rbi_g+pack_size,prcommVec[ii].pack_size,MPI_UINT8,
        prcommVec[ii].rank,12346,mpi_comm,&(recvRequestArray[ii]));
    pack_size += prcommVec[ii].pack_size;

    // and the send...
    MPI_Issend(rbi_g+(prcommVec[ii].unpack_offset-ncv),prcommVec[ii].unpack_size,MPI_UINT8,
        prcommVec[ii].rank,12346,mpi_comm,&(sendRequestArray[ii]));

  }

  MPI_Waitall(prcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
  MPI_Waitall(prcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

  delete[] sendRequestArray;
  delete[] recvRequestArray;

  // now build the packVec and periodicity...

  double R[9], t[3];
  pack_size = 0;
  for (int ii = 0, ii_end=prcommVec.size(); ii < ii_end; ++ii) {
    assert(prcommVec[ii].packVec.empty());
    prcommVec[ii].packVec.resize(prcommVec[ii].pack_size);
    CvPrcomm::Transform * transform = NULL;
    int bits_current = -1;
    for (int i = 0; i < prcommVec[ii].pack_size; ++i) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,recv_rbi_g[pack_size+i]);
      assert(rank == mpi_rank);
      if (bits > bits_current) {
        // bits are about to change, so complete any translate and rotate...
        if (transform) transform->end = i;
        // look for new rotations/translations...
        const bool has_R = PeriodicData::getPeriodicR(R,bits);
        const bool has_t = PeriodicData::getPeriodicT(t,bits);
        prcommVec[ii].transformVec.push_back(CvPrcomm::Transform(has_R,R,has_t,t,bits,i));
        transform = &(prcommVec[ii].transformVec.back());
        bits_current = bits;
      }
      else {
        assert(bits_current == bits);
      }
      assert((index >= 0)&&(index < ncv));
      prcommVec[ii].packVec[i] = index;
    }
    if (transform) transform->end = prcommVec[ii].pack_size;
    pack_size += prcommVec[ii].pack_size;
  }

  delete[] recv_rbi_g;

}

void StaticSolver::buildPrcomm(vector<CvPrcomm>& prcommVec,uint8 * rbi_g,const set<int>& rankSet,const int ncv,const int ncv_g) {

  COUT1("buildPrcomm()");

  // ------------------------------------------------
  // build communicator without assuming rank symmetry
  // ------------------------------------------------

  assert(prcommVec.empty());

  // build pack side using rbi_g...

  set<int>::iterator iter = rankSet.begin();
  int rank_current = -1;
  int bits_current = -1;
  CvPrcomm * prcomm = NULL;
  for (int icv = ncv; icv < ncv_g; ++icv) {
    int rank,bits,index;
    BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv-ncv]);
    if (rank > rank_current) {
      // we just switched ranks. If there was a prcomm, then
      // complete the size...
      if (prcomm)
        prcomm->unpack_size = icv-prcomm->unpack_offset;
      while ((iter != rankSet.end())&&((*iter) <= rank)) {
        if ((*iter) < rank) {
          // add an empty Prcomm...
          prcommVec.push_back(Prcomm());
          prcommVec.back().rank = *iter;
          prcommVec.back().unpack_offset = icv;
          prcommVec.back().unpack_size = 0;
        }
        ++iter;
      }
      rank_current = rank;
      bits_current = -1;
      prcommVec.push_back(CvPrcomm());
      prcomm = &prcommVec.back();
      prcomm->rank = rank;
      prcomm->unpack_offset = icv;
    }
    else {
      assert(rank_current == rank);
      assert(prcomm);
    }
    // just tests monotoncity of the bits...
    if (bits > bits_current) {
      bits_current = bits;
    }
    else {
      assert(bits_current == bits);
    }
  }
  // we just finished. If there was a prcomm, then
  // complete the size...
  if (prcomm) prcomm->unpack_size = ncv_g-prcomm->unpack_offset;
  // add any not represented in Prcomm at the end...
  while (iter != rankSet.end()) {
    // add an empty Prcomm...
    prcommVec.push_back(Prcomm());
    prcommVec.back().rank = *iter;
    prcommVec.back().unpack_offset = ncv_g;
    prcommVec.back().unpack_size = 0;
    ++iter;
  }
  assert(iter == rankSet.end());

  // finally, we need to send/recv the indices and bits to the pack side
  // and build the packVecs. Note that these should be symmetric by construction
  // because internal faces always have 2 cvs, but counts may not be the
  // same because more than one face can be on a single cv...

  MPI_Request * sendRequestArray = new MPI_Request[prcommVec.size()];
  MPI_Request * recvRequestArray = new MPI_Request[prcommVec.size()];
  for (int ii = 0, ii_end=prcommVec.size(); ii < ii_end; ++ii) {

    // post irecv...
    MPI_Irecv(&(prcommVec[ii].pack_size),1,MPI_INT,prcommVec[ii].rank,12345,mpi_comm,&(recvRequestArray[ii]));

    // and the send...
    MPI_Issend(&(prcommVec[ii].unpack_size),1,MPI_INT,prcommVec[ii].rank,12345,mpi_comm,&(sendRequestArray[ii]));

  }

  MPI_Waitall(prcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
  MPI_Waitall(prcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

  // now send from the unpack side to the pack side...

  int pack_size = 0;
  for (int ii = 0, ii_end=prcommVec.size(); ii < ii_end; ++ii)
    pack_size += prcommVec[ii].pack_size;

  uint8 * recv_rbi_g = new uint8[pack_size];

  pack_size = 0;
  for (int ii = 0, ii_end=prcommVec.size(); ii < ii_end; ++ii) {

    // post irecv...
    MPI_Irecv(recv_rbi_g+pack_size,prcommVec[ii].pack_size,MPI_UINT8,
        prcommVec[ii].rank,12346,mpi_comm,&(recvRequestArray[ii]));
    pack_size += prcommVec[ii].pack_size;

    // and the send...
    MPI_Issend(rbi_g+(prcommVec[ii].unpack_offset-ncv),prcommVec[ii].unpack_size,MPI_UINT8,
        prcommVec[ii].rank,12346,mpi_comm,&(sendRequestArray[ii]));

  }

  MPI_Waitall(prcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
  MPI_Waitall(prcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

  delete[] sendRequestArray;
  delete[] recvRequestArray;

  // now build the packVec and periodicity...

  double R[9], t[3];
  pack_size = 0;
  for (int ii = 0, ii_end=prcommVec.size(); ii < ii_end; ++ii) {
    assert(prcommVec[ii].packVec.empty());
    prcommVec[ii].packVec.resize(prcommVec[ii].pack_size);
    CvPrcomm::Transform * transform = NULL;
    int bits_current = -1;
    for (int i = 0; i < prcommVec[ii].pack_size; ++i) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,recv_rbi_g[pack_size+i]);
      assert(rank == mpi_rank);
      if (bits > bits_current) {
        // bits are about to change, so complete any translate and rotate...
        if (transform) transform->end = i;
        // look for new rotations/translations...
        const bool has_R = PeriodicData::getPeriodicR(R,bits);
        const bool has_t = PeriodicData::getPeriodicT(t,bits);
        prcommVec[ii].transformVec.push_back(CvPrcomm::Transform(has_R,R,has_t,t,bits,i));
        transform = &(prcommVec[ii].transformVec.back());
        bits_current = bits;
      }
      else {
        assert(bits_current == bits);
      }
      assert((index >= 0)&&(index < ncv));
      prcommVec[ii].packVec[i] = index;
    }
    if (transform) transform->end = prcommVec[ii].pack_size;
    pack_size += prcommVec[ii].pack_size;
  }

  delete[] recv_rbi_g;

}

void StaticSolver::buildXfa() {

  COUT1("StaticSolver::buildXfa()");

  // eventually x_fa should be read in from the StripedMesh, but for now
  // build it out of the nodes...

  assert(noofa_i);
  assert(noofa_v);
  assert(x_no);
  assert(n_fa);

  assert(x_fa == NULL);
  x_fa = new double[nfa][3];

  FOR_IFA {

    const double dx[3] = DIFF(x_vv[cvofa[ifa][1]],x_vv[cvofa[ifa][0]]);
    assert(DOT_PRODUCT(dx,n_fa[ifa]) > 0.0);

    FOR_I3 x_fa[ifa][i] = 0.0;
    double wgt = 0.0;

    double n_test[3] = { 0.0, 0.0, 0.0 };

    const int ino0 = noofa_v[noofa_i[ifa]];
    int ino1 = ino0;
    for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      const int ino2 = ino1;
      ino1 = noofa_v[nof];
      if ((ino1 != ino0)&&(ino2 != ino0)) {
        const double this_n[3] = TRI_NORMAL_2(x_no[ino0],x_no[ino1],x_no[ino2]);

        FOR_I3 n_test[i] += this_n[i];

        const double this_wgt = DOT_PRODUCT(this_n,n_fa[ifa]);
        FOR_I3 x_fa[ifa][i] += this_wgt*(x_no[ino0][i]+x_no[ino1][i]+x_no[ino2][i]);
        wgt += this_wgt;
      }
    }

    FOR_I3 n_test[i] = n_test[i]/2.0 + n_fa[ifa][i];
    if (MAG(n_test) > 1.0E-08) {
      cout << "wgt: " << wgt << " compare n_test: " << COUT_VEC(n_test) << " to n_fa: " << COUT_VEC(n_fa[ifa]) << endl;
    }

    //assert(wgt >= 0);
    FOR_I3 x_fa[ifa][i] /= wgt*3.0;

  }

}

void StaticSolver::buildFaocv() {

  COUT1("StaticSolver::buildFaocv()");

  assert(faocv_i == NULL);
  assert(faocv_v == NULL);
  assert(cvofa != NULL);

  faocv_i = new int[ncv+1];

  // count...

  FOR_ICV faocv_i[icv+1] = 0;

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    ++faocv_i[icv0+1];
    const int icv1 = cvofa[ifa][1];
    ++faocv_i[icv1+1];
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    ++faocv_i[icv0+1];
  }

  // allocate...

  faocv_i[0] = 0;
  FOR_ICV faocv_i[icv+1] += faocv_i[icv];
  const int faocv_s = faocv_i[ncv];

  faocv_v = new int[faocv_s];

  // and set...

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    faocv_v[faocv_i[icv0]++] = ifa;
    faocv_v[faocv_i[icv1]++] = ifa;
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    assert( icv0 >= ncv_i);
    faocv_v[faocv_i[icv0]++] = ifa;
  }

  // rewind...

  for (int icv = ncv-1; icv > 0; --icv)
    faocv_i[icv] = faocv_i[icv-1];
  faocv_i[0] = 0;


  // check the internal cell count ..

  for (int icv = 0; icv < ncv_i; ++icv) {
    for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
      const int ifa = faocv_v[foc];
      assert( (ifa < nfa_i) && (ifa >= 0));
      if ( cvofa[ifa][0] == icv) {
        assert( (cvofa[ifa][1] < ncv) && (cvofa[ifa][1] >= 0));
      } else {
        assert( cvofa[ifa][1] == icv);
        assert( (cvofa[ifa][0] < ncv) && (cvofa[ifa][0] >= 0));
      }
    }
  }

}

void StaticSolver::buildBfocv() {

  COUT1("StaticSolver::buildBfocv()");

  assert(bfocv_i == NULL);
  assert(bfocv_v == NULL);
  assert(cvobf != NULL);

  bfocv_i = new int[ncv+1];

  // count...

  FOR_ICV bfocv_i[icv+1] = 0;

  FOR_IBF {
    const int icv = cvobf[ibf];
    ++bfocv_i[icv+1];
  }

  // allocate...

  bfocv_i[0] = 0;
  FOR_ICV bfocv_i[icv+1] += bfocv_i[icv];
  const int bfocv_s = bfocv_i[ncv];

  bfocv_v = new int[bfocv_s];

  // and set...

  FOR_IBF {
    const int icv = cvobf[ibf];
    bfocv_v[bfocv_i[icv]++] = ibf;
  }

  // rewind...

  for (int icv = ncv-1; icv > 0; --icv)
    bfocv_i[icv] = bfocv_i[icv-1];
  bfocv_i[0] = 0;

}

void StaticSolver::buildCvocv() {

  COUT1("StaticSolver::buildCvocv()");

  assert(cvocv_i == NULL);
  assert(cvocv_v == NULL);
  assert(cvofa != NULL);

  cvocv_i = new int[ncv+1];

  // count...

  FOR_ICV cvocv_i[icv+1] = 1; // diagonal

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    ++cvocv_i[icv0+1];
    const int icv1 = cvofa[ifa][1];
    ++cvocv_i[icv1+1];
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    ++cvocv_i[icv0+1];
  }

  // allocate...

  cvocv_i[0] = 0;
  FOR_ICV cvocv_i[icv+1] += cvocv_i[icv];
  const int cvocv_s = cvocv_i[ncv];

  cvocv_v = new int[cvocv_s];

  // and set...

  FOR_ICV {
    cvocv_v[cvocv_i[icv]++] = icv; // diagonal
  }

  FOR_INTERNAL_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    cvocv_v[cvocv_i[icv0]++] = icv1;
    cvocv_v[cvocv_i[icv1]++] = icv0;
  }

  FOR_INTERPROC_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    cvocv_v[cvocv_i[icv0]++] = icv1;
  }

  // rewind...

  for (int icv = ncv-1; icv > 0; --icv)
    cvocv_i[icv] = cvocv_i[icv-1];
  cvocv_i[0] = 0;


  for (int icv = 0; icv < ncv_i; ++icv) {
    for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      assert( (icv_nbr >= 0) && (icv_nbr < ncv));
    }
  }

}

void StaticSolver::buildCvocv2() { 

  if (mpi_rank == 0) 
    cout << "buildCvocv2()" << endl;

  assert( cvocv2_i == NULL);
  assert( cvocv2_v == NULL);
  cvocv2_i = new int[ncv+1];
  for (int icv = 0; icv <= ncv; ++icv) 
    cvocv2_i[icv] = 0;

  for (int iter =0; iter < 2; ++iter) { 

    for (int ief = 0; ief < nef; ++ief) { 

      const int icv0 = cvoef[ief][0];
      const int icv1 = cvoef[ief][1];

      if ( iter == 0) { 

        ++cvocv2_i[icv0+1];
        if ( icv1 < ncv) 
          ++cvocv2_i[icv1+1];

      } else { 

        // for this face, icv1 is a neighbor of icv0 and vice-versa .. 

        cvocv2_v[cvocv2_i[icv0]] = icv1;
        ++cvocv2_i[icv0];

        if ( icv1 < ncv) { 
          cvocv2_v[cvocv2_i[icv1]] = icv0;
          ++cvocv2_i[icv1];
        }

      }

    }

    if ( iter == 0 ) { 

      cvocv2_i[0] = 0;
      for (int icv = 1; icv <= ncv; ++icv) 
        cvocv2_i[icv] += cvocv2_i[icv-1]+1; // add one for the diagonal

      cvocv2_v = new int[cvocv2_i[ncv]];

      // set the diagonal in the cvocv2_i operator .. 

      for (int icv = 0; icv < ncv; ++icv) {

        cvocv2_v[cvocv2_i[icv]] = icv;
        ++cvocv2_i[icv];

      }
    }

  }

  // rewind .. 

  for (int icv = ncv; icv > 0; --icv) 
    cvocv2_i[icv] = cvocv2_i[icv-1];
  cvocv2_i[0] = 0;

}

void StaticSolver::buildNoPrcomm() {

  COUT1("StaticSolver::buildNoPrcomm()");

  assert(ino_global);
  FOR_INO assert((ino_global[ino] >= 0)&&(ino_global[ino] < nno_global));

  int8 * noora_s = NULL; // s == striped
  MiscUtils::buildUniformXora(noora_s,nno_global);
  const int nno_s = noora_s[mpi_rank+1] - noora_s[mpi_rank];

  DistributedDataExchanger dde(ino_global,nno,noora_s);

  int * rbi_i_s = new int[nno_s+1];
  for (int ino_s = 0; ino_s < nno_s; ++ino_s)
    rbi_i_s[ino_s+1] = 0;
  dde.push_one(rbi_i_s+1,ADD_DATA);

  // pull the count back to the load-balanced distribution...

  int * rbi_i = new int[nno+1];
  dde.pull(rbi_i+1,rbi_i_s+1);

  // now build the csr on the striped distribution...

  rbi_i_s[0] = 0;
  for (int ino_s = 0; ino_s < nno_s; ++ino_s)
    rbi_i_s[ino_s+1] += rbi_i_s[ino_s];
  uint8 * rbi_v_s = new uint8[rbi_i_s[nno_s]];

  uint8 * rbi = new uint8[nno];
  FOR_INO rbi[ino] = BitUtils::packRankBitsIndex(mpi_rank,0,ino);
  dde.push_to_csr(rbi_i_s,rbi_v_s,rbi);
  delete[] rbi;

  /*
  // ====================================
  // costly, complete check...
  // ====================================
  {

  int * no_flag = new int[nno];
  FOR_INO no_flag[ino] = 0;

  FOR_RANK {

  int nron = rbi_i_s[nno_s];
  MPI_Bcast(&nron,1,MPI_INT,rank,mpi_comm);

  // this is bad form, but split the behavior around MPI_Bcast...

  uint8 buf[2];
  if (mpi_rank == rank) {
  for (int ino_s = 0; ino_s < nno_s; ++ino_s) {
  for (int ron = rbi_i_s[ino_s]; ron != rbi_i_s[ino_s+1]; ++ron) {
  buf[0] = noora_s[mpi_rank] + ino_s;
  buf[1] = rbi_v_s[ron];
  MPI_Bcast(buf,2,MPI_UINT8,rank,mpi_comm);
  int check_rank,check_bits,check_index;
  BitUtils::unpackRankBitsIndex(check_rank,check_bits,check_index,buf[1]);
  assert(check_bits == 0);
  if (check_rank == mpi_rank) {
  assert(ino_global[check_index] == buf[0]);
  assert(no_flag[check_index] == 0);
  no_flag[check_index] = 1;
  }
  }
  }
  }
  else {
  for (int ron = 0; ron < nron; ++ron) {
  MPI_Bcast(buf,2,MPI_UINT8,rank,mpi_comm);
  int check_rank,check_bits,check_index;
  BitUtils::unpackRankBitsIndex(check_rank,check_bits,check_index,buf[1]);
  assert(check_bits == 0);
  if (check_rank == mpi_rank) {
  assert(ino_global[check_index] == buf[0]);
  assert(no_flag[check_index] == 0);
  no_flag[check_index] = 1;
  }
  }
  }

  if (mpi_rank == 0)
  cout << " > rank " << rank << " looks good" << endl;

  }

  FOR_INO assert(no_flag[ino] == 1);

  if (mpi_rank == 0)
  cout << " > ALL nodes represented ONCE in striped csr." << endl;

  delete[] no_flag;
  }
  */

  // and pull the csr back to the load balanced distribution...
  // NOTE: this is not working, so revisit later.

  /*
     rbi_i[0] = 0;
     FOR_INO {
     assert(rbi_i[ino+1] >= 1);
     rbi_i[ino+1] += rbi_i[ino];
     }
     uint8 * rbi_v = new uint8[rbi_i[nno]];
     dde.pull_csr(rbi_i,rbi_v,rbi_i_s,rbi_v_s,ino_global,noora_s); // does no work -- do it the hard way
     */

  delete[] noora_s;

  // FOR NOW, do this the hard way, but also send less data...
  // go from recv (the striped side) to send -- here we send everybody
  // except out own entry...

  rbi_i[0] = 0;
  FOR_INO {
    assert(rbi_i[ino+1] >= 1);
    rbi_i[ino+1] -= 1; // remove one for ourselves. We still get our bit in the exchange and can handle periodicity?
    rbi_i[ino+1] += rbi_i[ino];
  }
  uint8 * rbi_v = new uint8[rbi_i[nno]];

  int * recv_count = new int[mpi_size];
  FOR_RANK recv_count[rank] = 0;

  for (int ino_s = 0; ino_s < nno_s; ++ino_s) {
    const int nron = rbi_i_s[ino_s+1] - rbi_i_s[ino_s];
    if (nron >= 2) {
      for (int ron = rbi_i_s[ino_s]; ron != rbi_i_s[ino_s+1]; ++ron) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_v_s[ron]);
        recv_count[rank] += nron; // everyone gets one less than the number, but also gets the index
      }
    }
  }

  int * recv_disp = new int[mpi_size];
  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
  const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

  uint8 * recv_buf = new uint8[recv_count_sum];

  for (int ino_s = 0; ino_s < nno_s; ++ino_s) {
    const int nron = rbi_i_s[ino_s+1] - rbi_i_s[ino_s];
    if (nron >= 2) {
      for (int ron = rbi_i_s[ino_s]; ron != rbi_i_s[ino_s+1]; ++ron) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_v_s[ron]);
        recv_buf[recv_disp[rank]++] = BitUtils::packRankBitsIndex(nron-1,0,index); // use a trick here to store the count where rank goes...
        for (int ron2 = rbi_i_s[ino_s]; ron2 != rbi_i_s[ino_s+1]; ++ron2) if (ron2 != ron) {
          recv_buf[recv_disp[rank]++] = rbi_v_s[ron2];
        }
      }
    }
  }

  // clear the striped csr...

  delete[] rbi_v_s;
  delete[] rbi_i_s;

  // rewind...

  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

  int * send_count = new int[mpi_size];
  MPI_Alltoall(recv_count,1,MPI_INT,send_count,1,MPI_INT,mpi_comm);

  int * send_disp = new int[mpi_size];
  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
  const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
  //assert(send_count_sum == rbi_i[nno]+nno);

  uint8 * send_buf = new uint8[send_count_sum];
  MPI_Alltoallv(recv_buf,recv_count,recv_disp,MPI_UINT8,
      send_buf,send_count,send_disp,MPI_UINT8,
      mpi_comm);
  delete[] recv_buf;
  delete[] recv_count;
  delete[] recv_disp;
  delete[] send_count;
  delete[] send_disp;

  int isend = 0;
  while (isend < send_count_sum) {
    int nron,bits,ino;
    BitUtils::unpackRankBitsIndex(nron,bits,ino,send_buf[isend++]);
    assert(nron > 0);
    assert(bits == 0);
    assert((ino >= 0)&&(ino < nno));
    assert(rbi_i[ino+1]-rbi_i[ino] == nron);
    for (int ron = 0; ron < nron; ++ron) {
      rbi_v[rbi_i[ino]+ron] = send_buf[isend++];
    }
  }
  assert(isend == send_count_sum);
  // delete[] send_buf; // we use this below

  // now we have a csr structure associated with every node on the
  // load balanced distribution. We built this the hard way and
  // do not have ourselves in the loop any more, but will verify later
  // but exchanging x_no...

  vector<pair<uint8,int> > rbi_index_pair_vec(rbi_i[nno]);
  FOR_INO {
    for (int ron = rbi_i[ino]; ron != rbi_i[ino+1]; ++ron) {
      rbi_index_pair_vec[ron].first = rbi_v[ron];
      rbi_index_pair_vec[ron].second = ino;
    }
  }

  delete[] rbi_i;
  delete[] rbi_v;

  sort(rbi_index_pair_vec.begin(),rbi_index_pair_vec.end());

  assert(noPrcommVec.empty());
  assert(send_count_sum >= int(rbi_index_pair_vec.size()));

  int rank_current = -1;
  int bits_current = -1;
  Prcomm * prcomm = NULL;
  for (int ii = 0, ii_end=rbi_index_pair_vec.size(); ii < ii_end; ++ii) {
    send_buf[ii] = rbi_index_pair_vec[ii].first;
    int rank,bits,index;
    BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_index_pair_vec[ii].first);
    assert(rank != mpi_rank);
    assert(bits == 0);
    if (rank > rank_current) {
      // we just switched ranks. If there was a prcomm, then
      // complete the size...
      if (prcomm) {
        prcomm->unpack_size = ii-prcomm->unpack_offset;
        prcomm->unpackVec.resize(prcomm->unpack_size);
        for (int ii2 = prcomm->unpack_offset; ii2 < ii; ++ii2) {
          prcomm->unpackVec[ii2-prcomm->unpack_offset] = rbi_index_pair_vec[ii2].second;
        }
      }
      rank_current = rank;
      bits_current = -1;
      noPrcommVec.push_back(Prcomm());
      prcomm = &noPrcommVec.back();
      prcomm->rank = rank;
      prcomm->unpack_offset = ii;
    }
    else {
      assert(rank_current == rank);
      assert(prcomm);
    }
    if (bits > bits_current) {
      bits_current = bits;
    }
    else {
      assert(bits_current == bits);
    }
  }
  // we just finished. If there was a prcomm, then
  // complete the size...
  if (prcomm) {
    prcomm->unpack_size = rbi_index_pair_vec.size()-prcomm->unpack_offset;
    prcomm->unpackVec.resize(prcomm->unpack_size);
    for (int ii2 = prcomm->unpack_offset, ii2_end=rbi_index_pair_vec.size(); ii2 < ii2_end; ++ii2) {
      prcomm->unpackVec[ii2-prcomm->unpack_offset] = rbi_index_pair_vec[ii2].second;
    }
  }

  // finally, we need to send/recv the indices and bits to the pack side
  // and build the packVecs. Note that these should be symmetric by construction
  // because internal faces always have 2 cvs, but counts may not be the
  // same because more than one face can be on a single cv...

  MPI_Request * sendRequestArray = new MPI_Request[noPrcommVec.size()];
  MPI_Request * recvRequestArray = new MPI_Request[noPrcommVec.size()];
  for (int ii = 0,ii_end=noPrcommVec.size(); ii < ii_end; ++ii) {

    // post irecv...
    MPI_Irecv(&(noPrcommVec[ii].pack_size),1,MPI_INT,noPrcommVec[ii].rank,12345,mpi_comm,&(recvRequestArray[ii]));

    // and the send...
    MPI_Issend(&(noPrcommVec[ii].unpack_size),1,MPI_INT,noPrcommVec[ii].rank,12345,mpi_comm,&(sendRequestArray[ii]));

  }

  MPI_Waitall(noPrcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
  MPI_Waitall(noPrcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

  // now send from the unpack side to the pack side...

  int pack_size = 0;
  for (int ii = 0,ii_end=noPrcommVec.size(); ii < ii_end; ++ii)
    pack_size += noPrcommVec[ii].pack_size;

  uint8 * recv_rbi = new uint8[pack_size];
  pack_size = 0;
  for (int ii = 0,ii_end=noPrcommVec.size(); ii < ii_end; ++ii) {

    // post irecv...
    MPI_Irecv(recv_rbi+pack_size,noPrcommVec[ii].pack_size,MPI_UINT8,
        noPrcommVec[ii].rank,12346,mpi_comm,&(recvRequestArray[ii]));
    pack_size += noPrcommVec[ii].pack_size;

    // and the send...
    MPI_Issend(send_buf+noPrcommVec[ii].unpack_offset,noPrcommVec[ii].unpack_size,MPI_UINT8,
        noPrcommVec[ii].rank,12346,mpi_comm,&(sendRequestArray[ii]));

  }

  MPI_Waitall(noPrcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
  MPI_Waitall(noPrcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

  delete[] send_buf;
  delete[] sendRequestArray;
  delete[] recvRequestArray;

  // now build the packVec (and periodicity in the future)...

  pack_size = 0;
  for (int ii = 0,ii_end=noPrcommVec.size(); ii < ii_end; ++ii) {
    assert(noPrcommVec[ii].packVec.empty());
    assert(noPrcommVec[ii].pack_size == noPrcommVec[ii].unpack_size);
    noPrcommVec[ii].packVec.resize(noPrcommVec[ii].pack_size);
    for (int i = 0; i < noPrcommVec[ii].pack_size; ++i) {
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,recv_rbi[pack_size+i]);
      assert(rank == mpi_rank);
      assert(bits == 0); // for now
      assert((index >= 0)&&(index < nno));
      noPrcommVec[ii].packVec[i] = index;
      if (i > 0) assert( noPrcommVec[ii].packVec[i] > noPrcommVec[ii].packVec[i-1] ); // should be valid for same bits
    }
    pack_size += noPrcommVec[ii].pack_size;
  }

  delete[] recv_rbi;

  // =====================================================
  // test that the node coordinates are exactly matched...
  // =====================================================

  updateNoDataCheck(x_no);

}

void StaticSolver::buildFaPrcomm() {

  COUT1("StaticSolver::buildFaPrcomm()");

  // ------------------------------------------------
  // build the paired communicator -- this should be
  // symmetric in terms of rank AND count...
  // ------------------------------------------------

  assert(faPrcommVec.empty());
  assert(rbi_g);

  vector<pair<uint8,int> > rbiVec(nfa-nfa_i);
  FOR_INTERPROC_IFA {
    const int icv1 = cvofa[ifa][1];
    assert((icv1 >= ncv)&&(icv1 < ncv_g));
    rbiVec[ifa-nfa_i].first = rbi_g[icv1-ncv];
    rbiVec[ifa-nfa_i].second = ifa;
  }

  sort(rbiVec.begin(),rbiVec.end());

  int rank_current = -1;
  int bits_current = -1;
  Prcomm * prcomm = NULL;
  for (int ii = 0, ii_end=rbiVec.size(); ii < ii_end; ++ii) {
    int rank,bits,index;
    BitUtils::unpackRankBitsIndex(rank,bits,index,rbiVec[ii].first);
    if (rank > rank_current) {
      // we just switched ranks. If there was a prcomm, then
      // complete the size...
      if (prcomm) {
        prcomm->unpack_size = ii-prcomm->unpack_offset;
        assert(prcomm->unpackVec.empty());
        prcomm->unpackVec.resize(prcomm->unpack_size);
        for (int ii2 = prcomm->unpack_offset; ii2 < ii; ++ii2) {
          prcomm->unpackVec[ii2-prcomm->unpack_offset] = rbiVec[ii2].second; // the ifa
        }
      }
      rank_current = rank;
      bits_current = -1;
      faPrcommVec.push_back(Prcomm());
      prcomm = &faPrcommVec.back();
      prcomm->rank = rank;
      prcomm->unpack_offset = ii;
    }
    else {
      // monotonicity check...
      assert(rank_current == rank);
      assert(prcomm);
    }
    if (bits > bits_current) {
      bits_current = bits;
    }
    else {
      assert(bits_current == bits);
    }
  }
  // we just finished. If there was a prcomm, then
  // complete the size...
  if (prcomm) {
    prcomm->unpack_size = rbiVec.size()-prcomm->unpack_offset;
    assert(prcomm->unpackVec.empty());
    prcomm->unpackVec.resize(prcomm->unpack_size);
    for (int ii2 = prcomm->unpack_offset,ii2_end=rbiVec.size(); ii2 < ii2_end; ++ii2) {
      prcomm->unpackVec[ii2-prcomm->unpack_offset] = rbiVec[ii2].second; // the ifa
    }
  }

  rbiVec.clear();

  // finally, we need to send/recv the indices and bits to the pack side
  // and build the packVecs. Note that these should be symmetric by construction
  // because faces are paired...

  MPI_Request * sendRequestArray = new MPI_Request[faPrcommVec.size()];
  MPI_Request * recvRequestArray = new MPI_Request[faPrcommVec.size()];
  for (int ii = 0, ii_end=faPrcommVec.size(); ii < ii_end; ++ii) {

    // post irecv...
    MPI_Irecv(&(faPrcommVec[ii].pack_size),1,MPI_INT,faPrcommVec[ii].rank,22345,mpi_comm,&(recvRequestArray[ii]));

    // and the send...
    MPI_Issend(&(faPrcommVec[ii].unpack_size),1,MPI_INT,faPrcommVec[ii].rank,22345,mpi_comm,&(sendRequestArray[ii]));

  }

  MPI_Waitall(faPrcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
  MPI_Waitall(faPrcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

  // now send from the unpack side to the pack side...

  uint8 * send_buf = new uint8[nfa-nfa_i];
  int pack_size = 0;
  for (int ii = 0, ii_end=faPrcommVec.size(); ii < ii_end; ++ii) {
    assert(faPrcommVec[ii].pack_size == faPrcommVec[ii].unpack_size);
    assert(int(faPrcommVec[ii].unpackVec.size()) == faPrcommVec[ii].unpack_size);
    assert(faPrcommVec[ii].unpack_offset == pack_size);
    for (int jj = 0; jj < faPrcommVec[ii].unpack_size; ++jj) {
      const int ifa = faPrcommVec[ii].unpackVec[jj];
      const int icv1 = cvofa[ifa][1];
      assert((icv1 >= ncv)&&(icv1 < ncv_g));
      int rank,bits,index;
      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv1-ncv]);
      assert(faPrcommVec[ii].rank == rank);
      assert((bits >= 0)&&(bits < (1<<6)));
      // trick here we use the rank part of the uint8 pack to hold the icv on the other side...
      send_buf[pack_size+jj] = BitUtils::packRankBitsIndex(index,bits,cvofa[ifa][0]);
    }
    pack_size += faPrcommVec[ii].pack_size;
  }
  assert(pack_size == nfa-nfa_i);

  uint8 * recv_rbi = new uint8[pack_size];
  pack_size = 0;
  for (int ii = 0, ii_end=faPrcommVec.size(); ii < ii_end; ++ii) {

    // post irecv...
    MPI_Irecv(recv_rbi+pack_size,faPrcommVec[ii].pack_size,MPI_UINT8,
        faPrcommVec[ii].rank,22346,mpi_comm,&(recvRequestArray[ii]));
    pack_size += faPrcommVec[ii].pack_size;

    // and the send...
    MPI_Issend(send_buf+faPrcommVec[ii].unpack_offset,faPrcommVec[ii].unpack_size,MPI_UINT8,
        faPrcommVec[ii].rank,22346,mpi_comm,&(sendRequestArray[ii]));

  }

  MPI_Waitall(faPrcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
  MPI_Waitall(faPrcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

  delete[] send_buf;
  delete[] sendRequestArray;
  delete[] recvRequestArray;

  // now build the packVec and periodicity...

  double R[9],t[3];
  pack_size = 0;
  for (int ii = 0, ii_end=faPrcommVec.size(); ii < ii_end; ++ii) {
    CvPrcomm::Transform * transform = NULL;
    int bits_current = -1;
    assert(faPrcommVec[ii].packVec.empty());
    assert(faPrcommVec[ii].pack_size == faPrcommVec[ii].unpack_size);
    faPrcommVec[ii].packVec.resize(faPrcommVec[ii].pack_size);
    for (int i = 0; i < faPrcommVec[ii].pack_size; ++i) {
      // note trick with rbi above -- rank holds icv0...
      int icv0,bits,index;
      BitUtils::unpackRankBitsIndex(icv0,bits,index,recv_rbi[pack_size+i]);
      assert((icv0 >= 0)&&(icv0 < ncv));
      assert((bits >= 0)&&(bits < (1<<6)));
      // periodic bits...
      if (bits > bits_current) {
        // bits are about to change, so complete any translate and rotate...
        if (transform) transform->end = i;
        // look for new rotations/translations...
        const bool has_R = PeriodicData::getPeriodicR(R,bits);
        const bool has_t = PeriodicData::getPeriodicT(t,bits);
        faPrcommVec[ii].transformVec.push_back(CvPrcomm::Transform(has_R,R,has_t,t,bits,i));
        transform = &(faPrcommVec[ii].transformVec.back());
        bits_current = bits;
      }
      else {
        assert(bits_current == bits);
      }
      // now loop through faces of icv0 to find rbi_g that matches. This is the face...
      const int inv_bits = BitUtils::flipPeriodicBits(bits);
      const uint8 rbi_nbr = BitUtils::packRankBitsIndex(faPrcommVec[ii].rank,inv_bits,index);
      int foc;
      for (foc = faocv_i[icv0]; foc != faocv_i[icv0+1]; ++foc) {
        const int ifa = faocv_v[foc];
        if (ifa >= nfa_i) {
          const int icv1 = cvofa[ifa][1];
          assert((icv1 >= ncv)&&(icv1 < ncv_g));
          if (rbi_g[icv1-ncv] == rbi_nbr) {
            faPrcommVec[ii].packVec[i] = ifa;
            break;
          }
        }
      }
      // ensure we found it...
      assert(foc != faocv_i[icv0+1]);
    }
    if (transform) transform->end = faPrcommVec[ii].pack_size;
    pack_size += faPrcommVec[ii].pack_size;
  }

  delete[] recv_rbi;

  {

    // check face correspondence...

    // 1:1 first...

    int * fa_flag = new int[nfa];
    FOR_INTERNAL_IFA fa_flag[ifa] = 2;
    FOR_INTERPROC_IFA fa_flag[ifa] = 1;
    updateFaData(fa_flag,ADD_DATA);
    FOR_IFA assert(fa_flag[ifa] == 2);
    delete[] fa_flag;

    // then x...

    assert(x_fa);
    double (*x_fa_check)[3] = new double[nfa][3];
    FOR_IFA FOR_I3 x_fa_check[ifa][i] = x_fa[ifa][i];
    updateFaData(x_fa_check,REPLACE_TRANSLATE_DATA);
    double my_d2_max = 0.0;
    FOR_IFA {
      const double d2 = DIST2(x_fa_check[ifa],x_fa[ifa]);
      my_d2_max = max(my_d2_max,d2);
    }
    delete[] x_fa_check;
    double d2_max;
    MPI_Reduce(&my_d2_max,&d2_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    if (mpi_rank == 0) cout << " > face communicator x check (should be zero): " << sqrt(d2_max) << endl;

    // then n...

    assert(n_fa);
    double (*n_fa_check)[3] = new double[nfa][3];
    FOR_IFA FOR_I3 n_fa_check[ifa][i] = n_fa[ifa][i];
    updateFaData(n_fa_check,REPLACE_ROTATE_DATA);
    my_d2_max = 0.0;
    FOR_INTERNAL_IFA {
      FOR_I3 assert(n_fa_check[ifa][i] == n_fa[ifa][i]); // these should be untouched
    }
    FOR_INTERPROC_IFA {
      double d2 = 0.0;
      FOR_I3 d2 += (n_fa_check[ifa][i]+n_fa[ifa][i])*(n_fa_check[ifa][i]+n_fa[ifa][i]); // inter-proc should be flipped
      my_d2_max = max(my_d2_max,d2);
    }
    delete[] n_fa_check;
    MPI_Reduce(&my_d2_max,&d2_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    if (mpi_rank == 0) cout << " > face communicator n check (should be zero): " << sqrt(d2_max) << endl;

    /*
       FOR_RANK {
       if (rank == mpi_rank) {
    // HACK XXX
    FOR_INTERPROC_IFA {
    cout << COUT_VEC(x_fa[ifa]) << " " << ifa_global[ifa] << endl; cout.flush();
    }

    }
    MPI_Pause("next");
    }
    */

  }

}

void StaticSolver::writeLpocvAndInitDdeStuff(MPI_File &fh,MPI_Offset &offset,const int lp_index) {

  lpHelperVec[lp_index].clear();

  assert(cvora_striped);
  CtiRegister::CtiData * data = CtiRegister::getRegisteredCtiData(lpHelperVec[lp_index].name + ":icv",false);
  assert(data);

  int * send_count = new int[mpi_size];
  FOR_RANK send_count[rank] = 0;

  const int nlp = data->size();
  int * pack_index = new int[nlp];
  for (int ilp = 0; ilp < nlp; ++ilp) {
    const int icv = data->in(ilp);
    assert((icv >= 0)&&(icv < ncv));
    const int rank = MiscUtils::getRankInXora(getIcvGlobal(icv),cvora_striped);
    send_count[rank] += 1;
    // temporarily store the rank in the pack...
    pack_index[ilp] = rank;
  }

  int * send_disp = new int[mpi_size];
  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
  const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
  assert(send_count_sum == nlp);

  // build the send_buf_int with the unpack cv index...

  int * send_buf_int = new int[send_count_sum];
  for (int ilp = 0; ilp < nlp; ++ilp) {
    const int icv = data->in(ilp);
    assert((icv >= 0)&&(icv < ncv));
    // recall the rank is in the pack index...
    const int rank = pack_index[ilp];
    const int icv_striped = getIcvGlobal(icv) - cvora_striped[rank];
    // switch the flag to the pack index in send_buf...
    pack_index[ilp] = send_disp[rank];
    send_buf_int[send_disp[rank]++] = icv_striped;
  }

  // rewind...

  send_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

  // setup recv side stuff...

  int * recv_count = new int[mpi_size];
  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

  int * recv_disp = new int[mpi_size];
  recv_disp[0] = 0;
  for (int rank = 1; rank < mpi_size; ++rank)
    recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
  const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

  int * recv_buf_int = new int[recv_count_sum];
  MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
      recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);

  // the recv_buf_int contains the cv index where the data needs to go...

  const int ncv_striped = cvora_striped[mpi_rank+1]-cvora_striped[mpi_rank];
  int8* lpocv_i_global = new int8[ncv_striped+1];
  for (int icv_striped = 0; icv_striped < ncv_striped; ++icv_striped)
    lpocv_i_global[icv_striped+1] = 0;

  for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
    const int icv_striped = recv_buf_int[irecv];
    assert((icv_striped >= 0)&&(icv_striped < ncv_striped));
    ++lpocv_i_global[icv_striped+1];
  }

  // turn this into csr...

  lpocv_i_global[0] = 0;
  for (int icv_striped = 0; icv_striped < ncv_striped; ++icv_striped)
    lpocv_i_global[icv_striped+1] += lpocv_i_global[icv_striped];
  assert(recv_count_sum == lpocv_i_global[ncv_striped]);

  // will send back new global index...

  for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
    const int icv_striped = recv_buf_int[irecv];
    assert((icv_striped >= 0)&&(icv_striped < ncv_striped));
    recv_buf_int[irecv] = lpocv_i_global[icv_striped]++;
  }

  // rewind...

  for (int icv_striped = ncv_striped-1; icv_striped > 0; --icv_striped)
    lpocv_i_global[icv_striped] = lpocv_i_global[icv_striped-1];
  lpocv_i_global[0] = 0;

  // build lpora_cv_striped...

  assert(lpocv_i_global[ncv_striped] < TWO_BILLION); // could relax this fairly easily
  assert(lpHelperVec[lp_index].lpora_cv_striped == NULL);
  const int nlp_striped = lpocv_i_global[ncv_striped];
  MiscUtils::buildXora(lpHelperVec[lp_index].lpora_cv_striped,nlp_striped);
  lpHelperVec[lp_index].nlp_global = lpHelperVec[lp_index].lpora_cv_striped[mpi_size];

  // now increase all our lpocv_i_global's by the displacement...

  for (int icv_striped = 0; icv_striped <= ncv_striped; ++icv_striped)
    lpocv_i_global[icv_striped] += lpHelperVec[lp_index].lpora_cv_striped[mpi_rank];

  MPI_Alltoallv(recv_buf_int,recv_count,recv_disp,MPI_INT,
      send_buf_int,send_count,send_disp,MPI_INT,mpi_comm);
  delete[] recv_buf_int;
  delete[] send_count;
  delete[] send_disp;
  delete[] recv_count;
  delete[] recv_disp;

  int8* ilp_global = new int8[nlp];
  for (int ilp = 0; ilp < nlp; ++ilp) {
    const int icv = data->in(ilp);
    assert((icv >= 0)&&(icv < ncv));
    const int rank = MiscUtils::getRankInXora(getIcvGlobal(icv),cvora_striped);
    ilp_global[ilp] = send_buf_int[pack_index[ilp]]+lpHelperVec[lp_index].lpora_cv_striped[rank];
  }
  delete[] send_buf_int;
  delete[] pack_index;

  lpHelperVec[lp_index].dde_cv_striped = new DistributedDataExchanger(ilp_global,nlp,lpHelperVec[lp_index].lpora_cv_striped);
  assert(lpHelperVec[lp_index].dde_cv_striped != NULL);
  delete[] ilp_global;

  // write lpocv...

  assert(cvora_striped);
  assert(getNcvGlobal() == cvora_striped[mpi_size]);

  if ( mpi_rank == 0 ) {
    Header header;
    string name = lpHelperVec[lp_index].name + ":lpocv_i";
    cout << " > LPOCV_I " << name << endl;
    sprintf(header.name,"%s",name.c_str());
    header.id     = UGP_IO_LPOCV_I;
    header.skip = header_size + (getNcvGlobal()+1)*int8_size;
    ByteSwap::setLswMswPairForInt8(header.idata+0,getNcvGlobal()+1);
    ByteSwap::setLswMswPairForInt8(header.idata+2,lpHelperVec[lp_index].nlp_global);
    MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
  }

  int n = cvora_striped[mpi_rank+1]-cvora_striped[mpi_rank];
  if (mpi_rank == mpi_size-1) ++n; // last rank writes the 1 extra entry in lpocv_i...

  writeChunkedData<int8>(fh,offset+header_size+cvora_striped[mpi_rank]*int8_size,lpocv_i_global,n,mpi_comm);
  delete[] lpocv_i_global;
  offset += header_size + (getNcvGlobal()+1)*int8_size;

}

void StaticSolver::writeData(const string& filename) {

  const double wtime0 = MPI_Wtime();

  // orient face data to respect global face direction...
  flipRWSignedFaData();

  // create tmp file to write too...
  MiscUtils::mkdir_for_file_collective(filename,0);
  string tmp_filename = MiscUtils::makeTmpPrefix(filename);
  MPI_File_delete(tmp_filename.c_str(),MPI_INFO_NULL);
  MPI_File fh;
  MPI_File_open(mpi_comm,tmp_filename.c_str(),MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

  if ( mpi_rank == 0 )  {
    int itmp[2] = {UGP_IO_MAGIC_NUMBER+1, 5}; // use a different magic number? -- to indicate file type as data, not restart
    cout << " > UGP_IO_VERSION: " << itmp[1] << endl;
    MPI_File_write(fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
  }

  MPI_Offset offset = int_size*2;

  if (mpi_rank == 0) {
    {
      //snapshot hash constructed from the mles hash and sles hash from which the calc began
      //includes all variable names as well as single valued variable values.
      stringstream ss;  ss << RestartHashUtilities::mlesHash << RestartHashUtilities::slesHash;
      for (map<const string,CtiRegister::CtiData>::iterator iter = CtiRegister::registeredDataMap.begin(); iter != CtiRegister::registeredDataMap.end(); ++iter) {
        if (iter->second.checkBit(WRITE_DATA)) {
          ss << iter->first; // add name to hash
          if (iter->second.getType() == D_DATA)
            ss << iter->second.d();
          if (iter->second.getType() == I_DATA)
            ss << iter->second.i();
        }
      }
      RestartHashUtilities::myHash.init(ss,RestartHashUtilities::sha1hashlength);
      //cout << "debug: sles hash string: " << ss.str() << endl;
      cout << " > hash " << RestartHashUtilities::myHash << endl;

      offset += RestartHashUtilities::writeHashToRestartHeader(fh);
    }

    //-----------
    // params...
    //-----------
    {
      cout << " > params" << endl;

      std::stringstream ss;
      CTI_Dump_solver_info(ss);
      dumpParams(ss);

      Header header;
      sprintf(header.name,"UGP_IO_PARAMS");
      header.id = UGP_IO_PARAMS;
      header.idata[0] = ss.str().length();
      header.skip = header_size + header.idata[0];

      MPI_File_write(fh,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      MPI_File_write(fh,(void*)ss.str().c_str(),header.idata[0],MPI_CHAR,MPI_STATUS_IGNORE);

      offset += header.skip;
    }
  }
  //sync everyone else at this offset...
  MPI_Bcast(&offset,1,MPI_OFFSET_DATATYPE,0,mpi_comm);

  // write all the lpocv's first...

  for (int ii = 0, lim = lpHelperVec.size(); ii < lim; ++ii)
    writeLpocvAndInitDdeStuff(fh,offset,ii);

  // write other registered data...

  CtiRegister::writeData(fh,offset);

  if ( mpi_rank == 0 ) {
    cout << " > EOF ... " << endl;
    Header header;
    header.id = UGP_IO_EOF;
    sprintf(header.name,"EOF");
    header.skip = header_size;
    MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
  }

  offset += header_size;

  MPI_File_set_size(fh,offset);
  MPI_File_close(&fh);

  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0) {
    const double seconds = MPI_Wtime()-wtime0;
    cout << " > write size: " << double(offset)/1.0E+9 << " [GB], write rate: " << double(offset)/(1.0E+9*seconds) << " [GB/s]" << endl;
  }

  // orient face data to respect local face direction...
  flipRWSignedFaData();

  // if we successfully wrote the tmp file, rename/replace the file with it...
  if (mpi_rank == 0) {
    remove(filename.c_str());
    rename(tmp_filename.c_str(),filename.c_str());
  }

}

void StaticSolver::registerStats(Param * param,const bool b_init) {

  // the list of stats are managed by CtiRegister, but we build the list here
  // because we may need to use the solver bcs to learn how to read/write the
  // stats to the disk for certain types...

  // b_init is false for stats that are known at init time because they
  // are registered data that may be read from teh restart. pass true
  // for stats that are requested during the run...

  COUT1("StaticSolver::registerStats");

  for (int iarg = 0; iarg < param->size(); ++iarg) {

    const string vname = param->getString(iarg);

    // catch wild-card expansion for boundary data here
    vector<string> var_names;
    if (vname.find("*:") != string::npos) {
      // insert expansion with all boundary zones, the data pointer checking below will filter
      // zones which don't 'have' this variable
      for (int izn=0, nzn=bfZoneVec.size(); izn<nzn; ++izn) {
        string data_name = vname;
        data_name = MiscUtils::replaceAll(data_name,"*:",(bfZoneVec[izn].getName()+":"));
        var_names.push_back(data_name);
      }
    }
    else {
      // if not a wild-card string, simply add this variable
      var_names.push_back(vname);
    }

    bool found = false;
    for (int v=0, vmax=var_names.size(); v<vmax; ++v) {
      CtiRegister::CtiData * data = CtiRegister::getCtiData(var_names[v],false);
      if (data != NULL) {
        if (mpi_rank == 0) cout << " > add stats for " << data->getTypeAsString() << "|" << data->getTopologyAsString() << " \"" << var_names[v] << "\"" << endl;

        // stats require the ability to do IO. Look for the dde stuff...
        if ((data->getType() >= DN_DATA)&&(!data->hasDdeStuff())) {
          CERR("this scalar/vector data has no dde stuff -- need to figure out how to build");
        }

        CtiRegister::registerStats(var_names[v],data,b_init);
        found = true;
      }
    }

    // only report issues if a single variable was attempted to register and failed; otherwise could be filtering of *:<var> and we don't want to report all of these
    if ((mpi_rank == 0) && (!found)) {
      // in the case of cht: stats, they are delayed because cht: registration is delayed...
      if (MiscUtils::startsWith(vname,"cht:")&&(!b_init)) {
	cout << " > Warning: CHT stats registration is delayed until later in the initialization. They will be available, but reset on each run" << endl;
      }
      else {
	cout << " > Warning: cannot evaluate STATS \"" << vname << "\" to any registered variables. Skipping." << endl;
      }
    }

  }

}

void StaticSolver::interpFromData() {

  if (checkParam("INTERP_FROM_RESTART")) {

    // preallocate
    int * send_count = new int[mpi_size];
    int * send_disp = new int[mpi_size];
    int * recv_count = new int[mpi_size];
    int * recv_disp = new int[mpi_size];
    double * fa_dist2 = new double[nfa];
    double * fa_dn_data = new double[nfa];
    double (*fa_dn3_data)[3] = new double[nfa][3];
    bool * cv_interp = new bool[ncv];
    FOR_ICV cv_interp[icv] = false; // flag if interpolated
    // for nearest nbr
    double *min_dist2 = new double[ncv];
    FOR_ICV min_dist2[icv] = HUGE_VAL;

    stringstream ss_hash; //build a string to identify the interpolated data in a hash

    int ndc = 0;
    FOR_PARAM_MATCHING("INTERP_FROM_RESTART") ndc++;

    int idc = 0;
    bool b_reset = false; // reset time and step to 0 (true if any INTERP_FROM_RESTART requests it)
    FOR_PARAM_MATCHING("INTERP_FROM_RESTART") {

      // put the data into the data container...

      int iarg = 0;
      const string filename = param->getString(iarg++);
      string secondary_filename = "";

      // certain files expect to be passed in pairs; parse this here so can do additional parameter parsing
      // later without interruption
      if (iarg < param->size()) {
        if ((filename == "mles")||(filename.find(".mles") != string::npos)) {
          // this version may not be necessary going forward. x_vv is now included in sles
          secondary_filename = param->getString(iarg++);
          if ((secondary_filename != "sles")&&(secondary_filename.find(".sles") == string::npos)) {
            CERR("expecting \"*.sles\" file after \"mles\", but \"" << secondary_filename << "\" does not have the correct file extension");
          }
        }
        else if ((filename.find(".cas") != string::npos) || (filename.find(".msh") != string::npos)) {
          secondary_filename = param->getString(iarg++);
          if (secondary_filename.find(".dat") == string::npos) {
            CERR("expecting \"*.dat\" file after \"cas/msh\", but \"" << secondary_filename << "\" does not have the correct file extension");
          }
        }
        else if (filename.find(".pbin") != string::npos) {
          secondary_filename = param->getString(iarg++);
          CERR("expecting \"*.dat\" file after \"pbin\", but \"" << secondary_filename << "\" does not have the correct file extension");
        }
      }

      // ===============================================================
      // we can pass hints to the data container. Use the currently
      // registered data...
      // ===============================================================
      vector<string> nameVec;

      CtiRegister::setRegisteredCvDnAndDn3Names(nameVec);
      set<string> varSet;
      for (int ii = 0, ii_end=nameVec.size(); ii < ii_end; ++ii)  {
        // skip registered coordinate data: it should never be interpolated!...
        if ((nameVec[ii] != "x_cv")&&(nameVec[ii] != "x_vv"))
          varSet.insert(nameVec[ii]);
      }

      nameVec.clear();
      CtiRegister::setRegisteredDAndINames(nameVec);
      for (int ii = 0, ii_end=nameVec.size(); ii < ii_end; ++ii)  {
        //if (mpi_rank == 0) cout << "XXXXXXXXXX d and i: " << nameVec[ii] << endl;
        varSet.insert(nameVec[ii]);
      }

      // additional filtering of VARS to interp
      // cull varSet based on user input
      // for now must be located immediately after file names
      if (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg));
        if ((token == "VAR")||(token == "VARS")) {
          ++iarg;
          // only insert listed VARS
          vector<string> varsVec;
          set<string> tmpVarSet;
          tmpVarSet.swap(varSet);
          MiscUtils::splitCsv(varsVec,param->getString(iarg++));
          set<string>::iterator it;
          for (int i=0,end=varsVec.size(); i < end; ++i) {
            it = tmpVarSet.find(varsVec[i]);
            if (it != tmpVarSet.end()) varSet.insert(varsVec[i]);
          }
        }
        else if ((token == "!VAR")||(token == "!VARS")||(token == "EXCLUDE_VAR")||(token == "EXCLUDE_VARS")) {
          ++iarg;
          // remove listed VARS
          vector<string> varsVec;
          MiscUtils::splitCsv(varsVec,param->getString(iarg++));
          for (int i=0,end=varsVec.size(); i < end; ++i) {
            varSet.erase(varsVec[i]);
          }
        }
      }

      if (varSet.empty()) {
        CWARN("no valid variables specified for interpolation; skipping...");
        continue;  // go to next INTERP_FROM_RESTART found
      }

      // now build the data container...
      DataContainer dc;

      if (filename.find(".les") != string::npos) {
        //dc.initFromOldRestart(filename,"u,T,step,time,dt"); // another example of passing hints
        int ierr = dc.initFromOldRestart(filename,varSet);
        if (ierr != 0) {
          CERR("initFromOldRestart failed");
        }
      }
      else if ((filename == "sles")||(filename.find(".sles") != string::npos)) {
        int ierr = dc.initFromSles(filename,varSet);
        if (ierr != 0) {
          CERR("DataContainer::initFromSles() failed");
        }
        // ss_hash stuff
        ss_hash << RestartHashUtilities::slesInterpHash;
        RestartHashUtilities::slesInterpHash.clear();
      }
      else if ((filename == "mles")||(filename.find(".mles") != string::npos)) {
        // this version may not be necessary going forward. x_vv is now included in sles
        assert((secondary_filename == "sles")||(secondary_filename.find(".sles") != string::npos));
        int ierr = dc.initFromRestart(filename,secondary_filename,varSet);
        if (ierr != 0) {
          CERR("DataContainer::initFromRestart failed");
        }
        //add interpolated file hashes to ss_hash;
        ss_hash << RestartHashUtilities::mlesInterpHash << RestartHashUtilities::slesInterpHash;
        RestartHashUtilities::mlesInterpHash.clear();
        RestartHashUtilities::slesInterpHash.clear();
      }
      else if ((filename.find(".cas") != string::npos) || (filename.find(".msh") != string::npos)) {
        int ierr = dc.initFromFluent(filename,secondary_filename,varSet);
        if (ierr != 0) {
          CERR("DataContainer::initFromFluent failed");
        }
      }
      else if (filename.find(".pbin") != string::npos) {
        assert(secondary_filename.find(".dat") != string::npos);
        int ierr = dc.initFromPbinFluent(filename,secondary_filename,varSet);
        if (ierr != 0) {
          CERR("DataContainer::initFromPbinFluent failed");
        }
      }
      else if (filename.find(".ip") != string::npos) {
        int ierr = dc.initFromFluentIp(filename,varSet);
        if (ierr != 0) {
          CERR("DataContainer::initFromFluentIp() failed");
        }
      }
      else {
        CERR("Unrecognized file extension in INTERP_FROM_RESTART");
      }

      // TODO pbin, fluent, ASCII

      // store hash in RestartHashUtilities::slesHash
      RestartHashUtilities::slesHash.clear();
      RestartHashUtilities::slesHash.init(ss_hash,RestartHashUtilities::sha1hashlength);

      int double_count = prepareInterp(dc,idc == 0);
      double_count += 3; // for the x_vv's

      // if the dc requires transform, process these in order...

      while (iarg < param->size()) {
        const string token = param->getString(iarg++);
        if (token == "ROTATE_X") {
          const double degrees = param->getDouble(iarg++);
          if (mpi_rank == 0) cout << " > ROTATE_X " << degrees << " degrees" << endl;
          dc.rotate(0,degrees);
          MiscUtils::dumpRange(dc.x_cv,dc.ncv,"dc.x_cv");
        }
        else if (token == "ROTATE_Y") {
          const double degrees = param->getDouble(iarg++);
          if (mpi_rank == 0) cout << " > ROTATE_Y " << degrees << " degrees" << endl;
          dc.rotate(1,degrees);
          MiscUtils::dumpRange(dc.x_cv,dc.ncv,"dc.x_cv");
        }
        else if (token == "ROTATE_Z") {
          const double degrees = param->getDouble(iarg++);
          if (mpi_rank == 0) cout << " > ROTATE_Z " << degrees << " degrees" << endl;
          dc.rotate(2,degrees);
          MiscUtils::dumpRange(dc.x_cv,dc.ncv,"dc.x_cv");
        }
        else if (token == "ROTATE") {
          double point[3]; FOR_I3 point[i] = param->getDouble(iarg++);
          double axis[3]; FOR_I3 axis[i] = param->getDouble(iarg++);
          const double degrees = param->getDouble(iarg++);
          dc.rotate(point,axis,degrees);
          MiscUtils::dumpRange(dc.x_cv,dc.ncv,"dc.x_cv");
        }
        else if (token == "ROTATE_VECTOR_DATA_X") {
          const double degrees = param->getDouble(iarg++);
          if (mpi_rank == 0) cout << " > ROTATE_VECTOR_DATA_X " << degrees << " degrees" << endl;
          dc.rotate(0,degrees,false);
        }
        else if (token == "ROTATE_VECTOR_DATA_Y") {
          const double degrees = param->getDouble(iarg++);
          if (mpi_rank == 0) cout << " > ROTATE_VECTOR_DATA_Y " << degrees << " degrees" << endl;
          dc.rotate(1,degrees,false);
        }
        else if (token == "ROTATE_VECTOR_DATA_Z") {
          const double degrees = param->getDouble(iarg++);
          if (mpi_rank == 0) cout << " > ROTATE_VECTOR_DATA_Z " << degrees << " degrees" << endl;
          dc.rotate(2,degrees,false);
        }
        else if (token == "ROTATE_VECTOR_DATA") {
          double point[3]; FOR_I3 point[i] = param->getDouble(iarg++);
          double axis[3]; FOR_I3 axis[i] = param->getDouble(iarg++);
          const double degrees = param->getDouble(iarg++);
          dc.rotate(point,axis,degrees,false);
        }
        else if (token == "TRANSLATE") {
          const double dx = param->getDouble(iarg++);
          const double dy = param->getDouble(iarg++);
          const double dz = param->getDouble(iarg++);
          if (mpi_rank == 0) cout << " > TRANSLATE " << dx << " " << dy << " " << dz << endl;
          dc.translate(dx,dy,dz);
          MiscUtils::dumpRange(dc.x_cv,dc.ncv,"dc.x_cv");
        }
        else if (token == "MIRROR") {
          double point[3]; FOR_I3 point[i] = param->getDouble(iarg++);
          double axis[3]; FOR_I3 axis[i] = param->getDouble(iarg++);
          dc.mirror(point,axis);
          MiscUtils::dumpRange(dc.x_cv,dc.ncv,"dc.x_cv");
        }
        else if (token == "RESET") {
          b_reset = true;
        }
        else {
          if (mpi_rank == 0) cout << "Warning: skipping unrecognized INTERP_FROM_RESTART token: " << token << endl;
        }
      }

      // build the bbox...

      // make sure we have the stuff we need...

      if (cvAdt == NULL) buildCvAdt();

      // ========================================
      // we now have everything we need...
      // ========================================

      if (mpi_rank == 0)
        cout << " > packing coarse grid data..." << endl;

      FOR_RANK send_count[rank] = 0;

      vector<int> intVec;
      for (int icv_ = 0; icv_ < dc.ncv; ++icv_) {
        assert(intVec.empty());
        cvBboxAdt->buildListForPoint(intVec,dc.x_cv[icv_]);
        for (int ii = 0, ii_end=intVec.size(); ii < ii_end; ++ii) {
          const int rank = intVec[ii]; assert((rank >= 0)&&(rank < mpi_size));
          ++send_count[rank];
        }
        intVec.clear();
      }

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
      const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

      MiscUtils::dumpRange(&send_count_sum,1,"packed data size");

      double * send_buf = new double[3*send_count_sum];
      int *cvois = new int[send_count_sum];
      for (int icv_ = 0; icv_ < dc.ncv; ++icv_) {
        assert(intVec.empty());
        cvBboxAdt->buildListForPoint(intVec,dc.x_cv[icv_]);
        for (int ii = 0, ii_end=intVec.size(); ii < ii_end; ++ii) {
          const int rank = intVec[ii]; assert((rank >= 0)&&(rank < mpi_size));
          cvois[send_disp[rank]] = icv_;
          send_buf[send_disp[rank]*3  ] = dc.x_cv[icv_][0];
          send_buf[send_disp[rank]*3+1] = dc.x_cv[icv_][1];
          send_buf[send_disp[rank]*3+2] = dc.x_cv[icv_][2];
          ++send_disp[rank];
        }
        intVec.clear();
      }

      // rewind...

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      // now send...

      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      MiscUtils::dumpRange(&recv_count_sum,1,"unpacked data size");

      double * recv_buf = new double[3*recv_count_sum];
      FOR_RANK {
        send_count[rank] *= 3;
        recv_count[rank] *= 3;
        send_disp[rank] *= 3;
        recv_disp[rank] *= 3;
      }
      MPI_Alltoallv(send_buf,send_count,send_disp,MPI_DOUBLE,
          recv_buf,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);

      // now unpack and set closest point index...

      int *irocv = new int[ncv];
      FOR_ICV irocv[icv] = -1;
      for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
        double xp[3]; FOR_I3 xp[i] = recv_buf[irecv*3+i];
        assert(intVec.empty());
        cvAdt->buildListForPoint(intVec,xp);
        for (int ii = 0, ii_end=intVec.size(); ii < ii_end; ++ii) {
          const int icv = intVec[ii]; assert((icv >= 0)&&(icv < ncv));
          const double dist2 = DIST2(xp,x_vv[icv]);
          if ((dist2 < min_dist2[icv])&&(dist2 <= r_vv[icv]*r_vv[icv])) { // HACK -- adjust tol here for exact point match
            min_dist2[icv] = dist2;
            irocv[icv] = irecv;
          }
        }
        intVec.clear();
      }

      // now set...
      
      for (list<Dn3Data>::iterator iter = dc.dn3List.begin(); iter != dc.dn3List.end(); ++iter) {
        if (iter->ctiData != NULL) {
          for (int isend = 0; isend < send_count_sum; ++isend) {
            const int icv_ = cvois[isend]; assert((icv_ >= 0)&&(icv_ < dc.ncv));
            FOR_I3 send_buf[isend*3+i] = iter->data[icv_][i];
          }
          MPI_Alltoallv(send_buf,send_count,send_disp,MPI_DOUBLE,
              recv_buf,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
          FOR_ICV {
            const int irecv = irocv[icv];
            if (irecv >= 0) {
              assert(irecv < recv_count_sum);
              FOR_I3 iter->ctiData->dn3(icv,i) = recv_buf[irecv*3+i];
            }
          }
        }
      }
      
      FOR_RANK {
        send_count[rank] /= 3;
        recv_count[rank] /= 3;
        send_disp[rank] /= 3;
        recv_disp[rank] /= 3;
      }

      for (list<DnData>::iterator iter = dc.dnList.begin(); iter != dc.dnList.end(); ++iter) {
        if (iter->ctiData != NULL) {
          for (int isend = 0; isend < send_count_sum; ++isend) {
            const int icv_ = cvois[isend]; assert((icv_ >= 0)&&(icv_ < dc.ncv));
            send_buf[isend] = iter->data[icv_];
          }
          MPI_Alltoallv(send_buf,send_count,send_disp,MPI_DOUBLE,
              recv_buf,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
          FOR_ICV {
            const int irecv = irocv[icv];
            if (irecv >= 0) {
              assert(irecv < recv_count_sum);
              iter->ctiData->dn(icv) = recv_buf[irecv];
            }
          }
        }
      }
      delete[] cvois;
      delete[] recv_buf;
      delete[] send_buf;
      
      int8 my_count[2] = { ncv, 0 };
      double my_dist_max = 0.0;
      FOR_ICV {
        const int irecv = irocv[icv];
        if (irecv >= 0) {
          assert(irecv < recv_count_sum);
          ++my_count[1];
          my_dist_max = max(my_dist_max,sqrt(min_dist2[icv])/r_vv[icv]); // remember the normalized distance
          // keep track of interpolated cells from all data sources
          cv_interp[icv] = true;
        }
      }
      delete[] irocv;

      int8 count[2];
      MPI_Reduce(my_count,count,2,MPI_INT8,MPI_SUM,0,mpi_comm);
      if (mpi_rank == 0)
        cout << " > set " << count[1] << " out of " << count[0] << " cvs." << endl;

      // check normalization and prox...

      double dist_max;
      MPI_Reduce(&my_dist_max,&dist_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      if (mpi_rank == 0)
        cout << " > dist_max: " << dist_max << endl;

      // only march after finishing all interpolations
      if (idc == (ndc-1) ) {

        // now use a marching algorithm to set any vars still unset...
        // here we use a face-based apporach to handle inter-processor data to allow
        // the registered CtiData to have arbitrary stride...

        if (mpi_rank == 0) cout << " > marching data:" << endl;

        // data set by the above process gets -1...
        FOR_ICV if (cv_interp[icv]) min_dist2[icv] = -1.0;
        delete[] cv_interp;

        {

          for (list<DnData>::iterator iter = dc.dnList.begin(); iter != dc.dnList.end(); ++iter) {
            if (iter->ctiData != NULL) {

              if (mpi_rank == 0) cout << " > " << iter->name << "...";

              // reset min_dist2...

              FOR_ICV if (min_dist2[icv] != -1.0) min_dist2[icv] = HUGE_VAL;

              int done = 0;
              while (done == 0) {

                if (mpi_rank == 0) {
                  cout << ".";
                  cout.flush();
                }

                int my_done = 1;

                FOR_INTERPROC_IFA {
                  const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
                  const int icv1 = cvofa[ifa][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g));
                  if (min_dist2[icv0] <= -1.0) {
                    fa_dist2[ifa] = DIST2(x_vv[icv1],x_vv[icv0]);
                    fa_dn_data[ifa] = iter->ctiData->dn(icv0);
                  }
                  else {
                    assert(min_dist2[icv0] == HUGE_VAL);
                    fa_dist2[ifa] = -1.0;
                    fa_dn_data[ifa] = 0.0;
                  }
                }

                updateFaData(fa_dist2,REPLACE_DATA);
                updateFaData(fa_dn_data,REPLACE_DATA);

                FOR_INTERPROC_IFA {
                  const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
                  if (min_dist2[icv0] >= 0.0) {
                    if ((fa_dist2[ifa] >= 0.0)&&(fa_dist2[ifa] < min_dist2[icv0])) {
                      // this nbr is wins...
                      min_dist2[icv0] = fa_dist2[ifa];
                      iter->ctiData->dn(icv0) = fa_dn_data[ifa];
                    }
                  }
                }

                FOR_INTERNAL_IFA {
                  const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
                  const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
                  const double dist2 = DIST2(x_vv[icv1],x_vv[icv0]);
                  if ((min_dist2[icv0] >= 0.0)&&(min_dist2[icv1] <= -1.0)) {
                    if (dist2 < min_dist2[icv0]) {
                      min_dist2[icv0] = dist2;
                      iter->ctiData->dn(icv0) = iter->ctiData->dn(icv1);
                    }
                  }
                  else if ((min_dist2[icv0] <= -1.0)&&(min_dist2[icv1] >= 0.0)) {
                    if (dist2 < min_dist2[icv1]) {
                      min_dist2[icv1] = dist2;
                      iter->ctiData->dn(icv1) = iter->ctiData->dn(icv0);
                    }
                  }
                }

                FOR_ICV {
                  if ((min_dist2[icv] >= 0.0)&&(min_dist2[icv] < HUGE_VAL)) {
                    min_dist2[icv] = -2.0; // use -2 to indicate a value that has been set by the marching
                  }
                  else if (min_dist2[icv] == HUGE_VAL) {
                    // any HUGE_VALs left indicate we are not done yet...
                    my_done = 0;
                  }
                  else {
                    assert(min_dist2[icv] <= -1.0);
                  }
                }

                MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);

              }

              if (mpi_rank == 0) cout << "OK" << endl;

            }
          }

        }

        {

          for (list<Dn3Data>::iterator iter = dc.dn3List.begin(); iter != dc.dn3List.end(); ++iter) {
            if (iter->ctiData != NULL) {

              if (mpi_rank == 0) cout << " > " << iter->name << "...";

              // reset min_dist2...

              FOR_ICV if (min_dist2[icv] != -1.0) min_dist2[icv] = HUGE_VAL;

              int done = 0;
              while (done == 0) {

                if (mpi_rank == 0) {
                  cout << ".";
                  cout.flush();
                }

                int my_done = 1;

                FOR_INTERPROC_IFA {
                  const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
                  const int icv1 = cvofa[ifa][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g));
                  if (min_dist2[icv0] <= -1.0) {
                    fa_dist2[ifa] = DIST2(x_vv[icv1],x_vv[icv0]);
                    FOR_I3 fa_dn3_data[ifa][i] = iter->ctiData->dn3(icv0,i);
                  }
                  else {
                    assert(min_dist2[icv0] == HUGE_VAL);
                    fa_dist2[ifa] = -1.0;
                    FOR_I3 fa_dn3_data[ifa][i] = 0.0;
                  }
                }

                updateFaData(fa_dist2,REPLACE_DATA);
                updateFaData(fa_dn3_data,REPLACE_DATA);

                FOR_INTERPROC_IFA {
                  const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
                  if (min_dist2[icv0] >= 0.0) {
                    if ((fa_dist2[ifa] >= 0.0)&&(fa_dist2[ifa] < min_dist2[icv0])) {
                      // this nbr is wins...
                      min_dist2[icv0] = fa_dist2[ifa];
                      FOR_I3 iter->ctiData->dn3(icv0,i) = fa_dn3_data[ifa][i];
                    }
                  }
                }

                FOR_INTERNAL_IFA {
                  const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
                  const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
                  const double dist2 = DIST2(x_vv[icv1],x_vv[icv0]);
                  if ((min_dist2[icv0] >= 0.0)&&(min_dist2[icv1] <= -1.0)) {
                    if (dist2 < min_dist2[icv0]) {
                      min_dist2[icv0] = dist2;
                      FOR_I3 iter->ctiData->dn3(icv0,i) = iter->ctiData->dn3(icv1,i);
                    }
                  }
                  else if ((min_dist2[icv0] <= -1.0)&&(min_dist2[icv1] >= 0.0)) {
                    if (dist2 < min_dist2[icv1]) {
                      min_dist2[icv1] = dist2;
                      FOR_I3 iter->ctiData->dn3(icv1,i) = iter->ctiData->dn3(icv0,i);
                    }
                  }
                }

                FOR_ICV {
                  if ((min_dist2[icv] >= 0.0)&&(min_dist2[icv] < HUGE_VAL)) {
                    min_dist2[icv] = -2.0; // use -2 to indicate a value that has been set by the marching
                  }
                  else if (min_dist2[icv] == HUGE_VAL) {
                    // any HUGE_VALs left indicate we are not done yet...
                    my_done = 0;
                  }
                  else {
                    assert(min_dist2[icv] <= -1.0);
                  }
                }

                MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);

              }

              if (mpi_rank == 0) cout << "OK" << endl;

            }
          }
        }
      }
      ++idc;

      for (list<DnData>::iterator iter = dc.dnList.begin(); iter != dc.dnList.end(); ++iter) {
        if (iter->ctiData != NULL) {
          dumpRange(iter->ctiData->getDNptr(),iter->ctiData->size(),iter->name);
        }
      }
      for (list<Dn3Data>::iterator iter = dc.dn3List.begin(); iter != dc.dn3List.end(); ++iter) {
        if (iter->ctiData != NULL) {
          dumpRange(iter->ctiData->getDN3ptr(),iter->ctiData->size(),iter->name);
        }
      }
    }

    // reset time and step if requested. we put this here because if the user only
    // reset the time on one of the INTERP_FROM_RESTART calls, we still want to
    // reset the time and step.
    if (b_reset) {
      COUT1(" > reseting time and step to 0");
      if ( CtiRegister::CtiData* data = CtiRegister::getRegisteredCtiData("step") ) {
        data->i() = 0;
      }
      if ( CtiRegister::CtiData* data = CtiRegister::getRegisteredCtiData("time") ) {
        data->d() = 0.0;
      }
    }

    delete[] fa_dn3_data;
    delete[] fa_dn_data;
    delete[] fa_dist2;
    delete[] min_dist2;
    delete[] recv_disp;
    delete[] recv_count;
    delete[] send_disp;
    delete[] send_count;
  }
}

int StaticSolver::prepareInterp(DataContainer& dc,bool first) {

  using namespace CtiRegister;

  COUT1("StaticSolver::prepareInterp()");

  // look for the simulations registered data that matches the data in the data container,
  // zero its value and set its flag...

  if (first)
    clearAllDataFlags();

  // actually set the I data and D data here.

  // I_DATA...

  for (list<pair<string,int> >::iterator iter = dc.iList.begin(); iter != dc.iList.end(); ++iter) {
    CtiData * ctiData = getRegisteredCtiData(iter->first);
    if (ctiData == NULL) {

      if (mpi_rank == 0) cout << " > did not match i \"" << iter->first << "\"" << endl;

    }
    else {

      if (mpi_rank == 0) cout << " > matched i \"" << iter->first << "\"" << endl;

      assert(ctiData->getType() == I_DATA);
      ctiData->i() = iter->second; // set value here
      ctiData->setFlag(1);

    }
  }

  // D_DATA...

  for (list<pair<string,double> >::iterator iter = dc.dList.begin(); iter != dc.dList.end(); ++iter) {
    CtiData * ctiData = getRegisteredCtiData(iter->first);
    if (ctiData == NULL) {

      if (mpi_rank == 0) cout << " > did not match d \"" << iter->first << "\"" << endl;

    }
    else {

      if (mpi_rank == 0) cout << " > matched d \"" << iter->first << "\"" << endl;

      assert(ctiData->getType() == D_DATA);
      ctiData->d() = iter->second;
      ctiData->setFlag(1);

    }
  }

  // =======================================================================
  // for field data, we just get the ctiData stuff connected...
  // =======================================================================

  int double_count = 0;

  for (list<DnData>::iterator iter = dc.dnList.begin(); iter != dc.dnList.end(); ++iter) {
    assert(iter->ctiData == NULL);
    iter->ctiData = getRegisteredCtiData(iter->name);
    if (iter->ctiData == NULL) {

      if (mpi_rank == 0) cout << " > did not match dn \"" << iter->name << "\"" << endl;

    }
    else {

      if (mpi_rank == 0) cout << " > matched dn \"" << iter->name << "\"" << endl;

      assert(iter->ctiData->getTopology() == CV_DATA);
      assert(iter->ctiData->getType() == DN_DATA);
      assert(iter->ctiData->size() == ncv);
      iter->ctiData->setFlag(1);
      if (first) FOR_ICV iter->ctiData->dn(icv) = 0.0;

      double_count += 1;

    }
  }

  for (list<Dn3Data>::iterator iter = dc.dn3List.begin(); iter != dc.dn3List.end(); ++iter) {
    assert(iter->ctiData == NULL);
    iter->ctiData = getRegisteredCtiData(iter->name);
    if (iter->ctiData == NULL) {

      if (mpi_rank == 0) cout << " > did not match dn3 \"" << iter->name << "\"" << endl;

    }
    else {

      if (mpi_rank == 0) cout << " > matched dn3 \"" << iter->name << "\"" << endl;

      assert(iter->ctiData->getTopology() == CV_DATA);
      assert(iter->ctiData->getType() == DN3_DATA);
      assert(iter->ctiData->size() == ncv);
      iter->ctiData->setFlag(1);
      if (first) FOR_ICV FOR_I3 iter->ctiData->dn3(icv,i) = 0.0;

      double_count += 3;

    }
  }

  return double_count;

}

void StaticSolver::buildCvLaplacian(double * A) {

  assert(A);
  assert(cvocv_i);
  assert(cvocv_v);

  for (int coc = 0; coc < cvocv_i[ncv]; ++coc)
    A[coc] = 0.0;

  for (int ifa = 0; ifa < nfa_i; ++ifa) {
    const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
    const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
    const double coeff = area_over_delta_fa[ifa];
    {
      int coc = cvocv_i[icv0];
      assert(cvocv_v[coc] == icv0);
      A[coc] -= coeff;
      ++coc;
      while (cvocv_v[coc] != icv1)
        ++coc;
      assert(coc < cvocv_i[icv0+1]);
      A[coc] += coeff;
    }
    {
      int coc = cvocv_i[icv1];
      assert(cvocv_v[coc] == icv1);
      A[coc] -= coeff;
      ++coc;
      while (cvocv_v[coc] != icv0)
        ++coc;
      assert(coc < cvocv_i[icv1+1]);
      A[coc] += coeff;
    }
  }

  for (int ifa = nfa_i; ifa < nfa; ++ifa) {
    const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
    const int icv1 = cvofa[ifa][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g));
    const double coeff = area_over_delta_fa[ifa];
    {
      int coc = cvocv_i[icv0];
      assert(cvocv_v[coc] == icv0);
      A[coc] -= coeff;
      ++coc;
      while (cvocv_v[coc] != icv1)
        ++coc;
      assert(coc < cvocv_i[icv0+1]);
      A[coc] += coeff;
    }
  }

}

void StaticSolver::buildCv2Laplacian(double * A) {

  assert(A);
  assert(cvocv2_i);
  assert(cvocv2_v);

  for (int coc = 0; coc < cvocv2_i[ncv]; ++coc)
    A[coc] = 0.0;

  for (int ief = 0; ief < nef_i; ++ief) {
    const int icv0 = cvoef[ief][0]; assert((icv0 >= 0)&&(icv0 < ncv));
    const int icv1 = cvoef[ief][1]; assert((icv1 >= 0)&&(icv1 < ncv));
    // approximate area over delta for extended face...
    const double coeff = sqrt(DOT_PRODUCT(n_ef[ief],n_ef[ief])/DIST2(x_cv[icv1],x_cv[icv0])); 
    {
      int coc = cvocv2_i[icv0];
      assert(cvocv2_v[coc] == icv0);
      A[coc] -= coeff;
      ++coc;
      while (cvocv2_v[coc] != icv1)
        ++coc;
      assert(coc < cvocv2_i[icv0+1]);
      A[coc] += coeff;
    }
    {
      int coc = cvocv2_i[icv1];
      assert(cvocv2_v[coc] == icv1);
      A[coc] -= coeff;
      ++coc;
      while (cvocv2_v[coc] != icv0)
        ++coc;
      assert(coc < cvocv2_i[icv1+1]);
      A[coc] += coeff;
    }
  }

  for (int ief = nef_i; ief < nef; ++ief) {
    const int icv0 = cvoef[ief][0]; assert((icv0 >= 0)&&(icv0 < ncv));
    const int icv1 = cvoef[ief][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g2));
    // approximate area over delta for extended face...
    const double coeff = sqrt(DOT_PRODUCT(n_ef[ief],n_ef[ief])/DIST2(x_cv[icv1],x_cv[icv0]));
    {
      int coc = cvocv2_i[icv0];
      assert(cvocv2_v[coc] == icv0);
      A[coc] -= coeff;
      ++coc;
      while (cvocv2_v[coc] != icv1)
        ++coc;
      assert(coc < cvocv2_i[icv0+1]);
      A[coc] += coeff;
    }
  }

}

void StaticSolver::calcCvResidual(double *res,const double* phi,const double *A,const double *rhs) {

  for (int icv = 0; icv < ncv; ++icv) {
    const int coc_f = cvocv_i[icv];
    res[icv] = rhs[icv] - A[coc_f]*phi[icv];
    for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      res[icv] -= A[coc]*phi[icv_nbr];
    }
  }

}

int StaticSolver::solveCvCg(double * phi,const double * const A,const double * const rhs,const double zero,const int maxiter,const bool verbose) {

  // assume we come in with a consistent initial condition...

  // we need the following work arrays...

  double * res      = new double[ncv];
  double * v        = new double[ncv];
  double * p        = new double[ncv_g];
  double * inv_diag = new double[ncv];

  // initialize...
  for (int icv = 0; icv < ncv; ++icv)
    inv_diag[icv] = 1.0/A[cvocv_i[icv]];

  for (int icv = 0; icv < ncv; ++icv)
    p[icv] = 0.0;
  double rho = 1.0;

  // calculate the residual in rhs format...
  for (int icv = 0; icv < ncv; ++icv) {
    res[icv] = rhs[icv] - A[cvocv_i[icv]]*phi[icv];
    for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      res[icv] -= A[coc]*phi[icv_nbr];
    }
  }

  // diagonal precon/compute normalized residual...
  for (int icv = 0; icv < ncv; ++icv)
    v[icv] = res[icv]*inv_diag[icv];

  int iter = 0;
  int done = 0;
  while (done == 0) {

    ++iter;

    double rho_prev = rho;
    if (fabs(rho_prev) < 1.0E-20)
      rho_prev = -1.0E-20; // -1.0E-20? seems to help

    double my_rho = 0.0;
    for (int icv = 0; icv < ncv; ++icv)
      my_rho += res[icv]*v[icv];
    MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    double beta = rho/rho_prev;
    for (int icv = 0; icv < ncv; ++icv)
      p[icv] = v[icv] + beta*p[icv];
    updateCvData(p);

    // v = [Ap]{p}...
    for (int icv = 0; icv < ncv; ++icv) {
      v[icv] = A[cvocv_i[icv]]*p[icv];
      for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        v[icv] += A[coc]*p[icv_nbr];
      }
    }

    double my_gamma = 0.0;
    for (int icv = 0; icv < ncv; ++icv)
      my_gamma += p[icv]*v[icv];
    double gamma;
    MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    if (fabs(gamma) < 1.0E-20)
      gamma = 1.0E-20;

    const double alpha = rho/gamma;

    // check if we are done...
    if (iter%3 == 0) {

      for (int icv = 0; icv < ncv; ++icv)
        phi[icv] += alpha*p[icv];
      updateCvData(phi);

      // recompute the residual...
      for (int icv = 0; icv < ncv; ++icv) {
        res[icv] = rhs[icv] - A[cvocv_i[icv]]*phi[icv];
        for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          res[icv] -= A[coc]*phi[icv_nbr];
        }
      }

      for (int icv = 0; icv < ncv; ++icv)
        v[icv] = res[icv]*inv_diag[icv];

      // compute the max (L-infinity) normalized residual...
      double  my_res_max = 0.0;
      for (int icv = 0; icv < ncv; ++icv)
        my_res_max = max( my_res_max, fabs(v[icv]) );
      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank == 0) {
        // only share the last half of the convergence behaviour...
        if (verbose || (iter > maxiter/2))
          cout << " > solveCvCg iter " << iter << " res_max " << res_max << endl;
        if (res_max <= zero) {
          //cout << "-> Successfully converged error to " << res_max << endl;
          done = 1;
        }
        else if (iter > maxiter) {
          cout << "Warning: solveCvCg did not converge after " << maxiter <<
            " iters, res_max: " << res_max << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

    }
    else {

      // update full phi including ghosts...
      for (int icv = 0; icv < ncv_g; ++icv)
        phi[icv] += alpha*p[icv];

      for (int icv = 0; icv < ncv; ++icv) {
        // on the other iterations, use this approximation to update
        // the unreduced residual...
        res[icv] -= alpha*v[icv];
        // still need to compute v, diag precon for next iteration...
        v[icv] = res[icv]*inv_diag[icv];
      }

    }

  }

  delete[] res;
  delete[] v;
  delete[] p;
  delete[] inv_diag;

  // let the calling routine know if we were successful...
  return( done == 1 );

}

int StaticSolver::solveCvCg(double (*u)[3],const double * const A,const double (* const rhs)[3],const double zero,const int maxiter,const bool verbose) {

  // assume we come in with a consistent initial condition...

  // we need the following work arrays...

  double (*res)[3]      = new double[ncv][3];
  double (*v)[3]        = new double[ncv][3];
  double (*p)[3]        = new double[ncv_g][3];
  double *inv_diag = new double[ncv];

  // initialize...
  for (int icv = 0; icv < ncv; ++icv)
    inv_diag[icv] = 1.0/A[cvocv_i[icv]];

  for (int icv = 0; icv < ncv; ++icv)
    FOR_I3 p[icv][i] = 0.0;
  double rho[3] = {1.0,1.0,1.0};

  // calculate the residual in rhs format...
  for (int icv = 0; icv < ncv; ++icv) {
    FOR_I3 res[icv][i] = rhs[icv][i] - A[cvocv_i[icv]]*u[icv][i];
    for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      FOR_I3 res[icv][i] -= A[coc]*u[icv_nbr][i];
    }
  }

  // diagonal precon/compute normalized residual...
  for (int icv = 0; icv < ncv; ++icv)
    FOR_I3 v[icv][i] = res[icv][i]*inv_diag[icv];

  int iter = 0;
  int done = 0;
  while (done == 0) {

    ++iter;

    double rho_prev[3] = {rho[0],rho[1],rho[2]};
    FOR_I3 {
      if (fabs(rho_prev[i]) < 1.0E-20)
        rho_prev[i] = -1.0E-20; // -1.0E-20? seems to help
    }

    double my_rho[3] = {0.0,0.0,0.0};
    for (int icv = 0; icv < ncv; ++icv)
      FOR_I3 my_rho[i] += res[icv][i]*v[icv][i];
    MPI_Allreduce(my_rho,rho,3,MPI_DOUBLE,MPI_SUM,mpi_comm);

    double beta[3]; FOR_I3 beta[i] = rho[i]/rho_prev[i];
    for (int icv = 0; icv < ncv; ++icv)
      FOR_I3 p[icv][i] = v[icv][i] + beta[i]*p[icv][i];
    updateCvData(p);

    // v = [Ap]{p}...
    for (int icv = 0; icv < ncv; ++icv) {
      FOR_I3 v[icv][i] = A[cvocv_i[icv]]*p[icv][i];
      for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        FOR_I3 v[icv][i] += A[coc]*p[icv_nbr][i];
      }
    }

    double my_gamma[3] = {0.0,0.0,0.0};
    for (int icv = 0; icv < ncv; ++icv)
      FOR_I3 my_gamma[i] += p[icv][i]*v[icv][i];
    double gamma[3];
    MPI_Allreduce(my_gamma,gamma,3,MPI_DOUBLE,MPI_SUM,mpi_comm);
    FOR_I3 {
      if (fabs(gamma[i]) < 1.0E-20)
        gamma[i] = 1.0E-20;
    }

    double alpha[3]; FOR_I3 alpha[i] = rho[i]/gamma[i];

    // check if we are done...
    if (iter%3 == 0) {

      for (int icv = 0; icv < ncv; ++icv)
        FOR_I3 u[icv][i] += alpha[i]*p[icv][i];
      updateCvData(u);

      // recompute the residual...
      for (int icv = 0; icv < ncv; ++icv) {
        FOR_I3 res[icv][i] = rhs[icv][i] - A[cvocv_i[icv]]*u[icv][i];
        for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          FOR_I3 res[icv][i] -= A[coc]*u[icv_nbr][i];
        }
      }

      for (int icv = 0; icv < ncv; ++icv)
        FOR_I3 v[icv][i] = res[icv][i]*inv_diag[icv];

      // compute the max (L-infinity) normalized residual...
      double  my_res_max = 0.0;
      for (int icv = 0; icv < ncv; ++icv)
        FOR_I3 my_res_max = max( my_res_max, fabs(v[icv][i]) );
      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank == 0) {
        // only share the last half of the convergence behaviour...
        if (verbose || (iter > maxiter/2))
          cout << " > solveCvCg iter " << iter << " res_max " << res_max << endl;
        if (res_max <= zero) {
          //cout << "-> Successfully converged error to " << res_max << endl;
          done = 1;
        }
        else if (iter > maxiter) {
          cout << "Warning: solveCvCg did not converge after " << maxiter <<
            " iters, res_max: " << res_max << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

    }
    else {

      // update full phi including ghosts...
      for (int icv = 0; icv < ncv_g; ++icv)
        FOR_I3 u[icv][i] += alpha[i]*p[icv][i];

      for (int icv = 0; icv < ncv; ++icv) {
        // on the other iterations, use this approximation to update
        // the unreduced residual...
        FOR_I3 res[icv][i] -= alpha[i]*v[icv][i];
        // still need to compute v, diag precon for next iteration...
        FOR_I3 v[icv][i] = res[icv][i]*inv_diag[icv];
      }

    }

  }

  delete[] res;
  delete[] v;
  delete[] p;
  delete[] inv_diag;

  // let the calling routine know if we were successful...
  return( done == 1 );

}

int StaticSolver::solveCvJacobi(double * phi,const double *A,const double *rhs,
    const double zero,const double relax,const int maxiter,const bool verbose) {

  double (*res)      = new double[ncv];

  int iter = 0;
  int done = 0;
  while (done == 0) {

    iter++;

    // calculate the residual...
    FOR_ICV {
      const int coc_f = cvocv_i[icv];
      res[icv] = rhs[icv] - A[coc_f]*phi[icv];
      for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        res[icv] -= A[coc]*phi[icv_nbr];
      }
      res[icv] /= A[coc_f];
    }

    // update the active u's...
    FOR_ICV phi[icv]   += relax*res[icv];

    // and ghosts...
    updateCvData(phi);

    // check...
    double my_res_max = 0.0;
    FOR_ICV my_res_max = max(my_res_max,fabs(res[icv]));
    double res_max;
    MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
    if (mpi_rank == 0) {
      if ((verbose)||(iter > maxiter/2))
        cout << " > solveCvJacobi iter, res_max: " << iter << " " << res_max << endl;
      if (res_max < zero) {
        done = 1;
      }
      else if (iter > maxiter) {
        cout << " > Warning: solveCvJacobi did not converge after " << maxiter << " iters, res_max: " << res_max << endl;
        done = 2;
      }
    }
    MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

  }

  delete[] res;

  return( done == 1 );

}

void StaticSolver::calcCvResidual(double (*res)[3],const double (*u)[3],const double *A,const double (*rhs)[3]) {

  for (int icv = 0; icv < ncv; ++icv) {
    const int coc_f = cvocv_i[icv];
    FOR_I3 res[icv][i] = rhs[icv][i] - A[coc_f]*u[icv][i];
    for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      FOR_I3 res[icv][i] -= A[coc]*u[icv_nbr][i];
    }
  }

}

int StaticSolver::solveCvJacobi(double (*__restrict__ u)[3],const double *__restrict__ A,const double (*__restrict__ rhs)[3],
    const double zero,const double relax,const int maxiter,const bool verbose) {

  double (*res)[3] = new double[ncv][3];

  double * inv_diag = new double[ncv];
  for (int icv = 0; icv < ncv; ++icv) {
    const int coc_f = cvocv_i[icv];
    inv_diag[icv] = 1.0/A[coc_f];
  }

  int iter = 0;
  int done = 0;
  while (done == 0) {

    iter++;

    // calculate the residual...
    FOR_ICV {
      const int coc_f = cvocv_i[icv];
      const int coc_e = cvocv_i[icv+1];
      FOR_I3 res[icv][i] = rhs[icv][i] - A[coc_f]*u[icv][i];
      for (int coc = coc_f+1; coc < coc_e; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        FOR_I3 res[icv][i] -= A[coc]*u[icv_nbr][i];
      }
      //FOR_I3 res[icv][i] /= A[coc_f];
      FOR_I3 res[icv][i] *= inv_diag[icv];
    }

    // update the active u's...
    FOR_ICV FOR_I3 u[icv][i]   += relax*res[icv][i];

    // and ghosts...
    updateCvData(u,REPLACE_ROTATE_DATA);

    // check...
    if ( iter%3 == 0) {
      double my_res_max = 0.0;
      FOR_ICV FOR_I3 my_res_max = max(my_res_max,fabs(res[icv][i]));
      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank == 0) {
        if ((verbose)||(iter > maxiter/2))
          cout << " > solveCvJacobi iter, res_max: " << iter << " " << res_max << endl;
        if (res_max < zero) {
          done = 1;
        }
        else if (iter > maxiter) {
          cout << " > Warning: solveCvJacobi did not converge after " << maxiter << " iters, res_max: " << res_max << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
    }

  }

  delete[] res;
  delete[] inv_diag;

  return( done == 1 );

}

void StaticSolver::calcCvResidual(double (*res)[3],const double (*u)[3],const double *A,const double (*A_diag)[3],const double (*rhs)[3]) {

  for (int icv = 0; icv < ncv; ++icv) {
    const int coc_f = cvocv_i[icv];
    FOR_I3 res[icv][i] = rhs[icv][i] - (A[coc_f]+A_diag[icv][i])*u[icv][i];
    for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      FOR_I3 res[icv][i] -= A[coc]*u[icv_nbr][i];
    }
  }

}

int StaticSolver::solveCvJacobi(double (*__restrict__ u)[3],const double *__restrict__ A,const double (*__restrict__ A_diag)[3],const double (*__restrict__ rhs)[3],
    const double zero,const double relax,const int maxiter,const bool verbose) {

  double (*res)[3] = new double[ncv][3];

  double (*inv_diag)[3] = new double[ncv][3];
  for (int icv = 0; icv < ncv; ++icv) {
    const int coc_f = cvocv_i[icv];
    FOR_I3 inv_diag[icv][i] = 1.0/(A[coc_f]+A_diag[icv][i]);
  }

  int iter = 0;
  int done = 0;
  while (done == 0) {

    iter++;

    // calculate the residual...
    FOR_ICV {
      const int coc_f = cvocv_i[icv];
      const int coc_e = cvocv_i[icv+1];
      FOR_I3 res[icv][i] = rhs[icv][i] - (A[coc_f]+A_diag[icv][i])*u[icv][i];
      for (int coc = coc_f+1; coc < coc_e; ++coc) {
        const int icv_nbr = cvocv_v[coc];
        FOR_I3 res[icv][i] -= A[coc]*u[icv_nbr][i];
      }
      //FOR_I3 res[icv][i] /= A[coc_f];
      FOR_I3 res[icv][i] *= inv_diag[icv][i];
    }

    // update the active u's...
    FOR_ICV FOR_I3 u[icv][i]   += relax*res[icv][i];

    // and ghosts...
    updateCvData(u,REPLACE_ROTATE_DATA);

    // check...
    if ( iter%3 == 0) {
      double my_res_max = 0.0;
      FOR_ICV FOR_I3 my_res_max = max(my_res_max,fabs(res[icv][i]));
      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank == 0) {
        if ((verbose)||(iter > maxiter/2))
          cout << " > solveCvJacobi iter, res_max: " << iter << " " << res_max << endl;
        if (res_max < zero) {
          done = 1;
        }
        else if (iter > maxiter) {
          cout << " > Warning: solveCvJacobi did not converge after " << maxiter << " iters, res_max: " << res_max << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
    }

  }

  delete[] res;
  delete[] inv_diag;

  return( done == 1 );

}


int StaticSolver::solveCvPatr(double *phi,const double * A,const double * At, const double *rhs,
    const double zero,const double relax1,const double relax2,const int maxiter,const bool verbose) {

  // assume we come in with a consistent initial condition...

  // we need the following work arrays...
  double *r   = new double[ncv_g];
  double *p   = new double[ncv_g];
  double *inv_diag = new double[ncv_g];
  double *Ar  = new double[ncv];
  double *Ap  = new double[ncv];

  // normalization...
  FOR_ICV inv_diag[icv] = 1.0/A[cvocv_i[icv]];
  updateCvData(inv_diag);

  // residial: r = b-A*x...
  FOR_ICV {
    r[icv] = A[cvocv_i[icv]]*phi[icv] - rhs[icv];
    for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      r[icv] += A[coc]*phi[icv_nbr];
    }
    r[icv] *= inv_diag[icv];
  }

  int iter = 0;
  int done = 0;
  while (done == 0) {

    // check if we are done...
    if (iter%10 == 0) {

      // compute the max (L-infinity) normalized residual...
      double  my_res_max = 0.0;
      FOR_ICV  my_res_max = max( my_res_max, fabs(r[icv]) );
      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank == 0) {
        // only share the last half of the convergence behaviour...
        if (verbose || (iter > maxiter/2))
          cout << " > solveCvPatr iter " << iter << " res_max " << res_max << endl;
        if (res_max <= zero) {
          //cout << "-> Successfully converged error to " << res_max << endl;
          done = 1;
        }
        else if (iter > maxiter) {
          cout << "Warning: solveCvPatr did not converge after " << maxiter <<
            " iters, res_max: " << res_max << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

    }

    // if we are not done, get the search direction and mag, and update guess...
    if (done == 0) {
      ++iter;

      // need residual in ghosts...
      updateCvData(r);

      // search direction: p = [A^T]{r}...
      FOR_ICV {
        p[icv] = At[cvocv_i[icv]]*inv_diag[icv]*r[icv];
        for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          p[icv] += At[coc]*inv_diag[icv_nbr]*r[icv_nbr];
        }
      }
      updateCvData(p);

      // second search direction is r...

      // get work arrays: Ap and Ar...
      FOR_ICV {
        Ar[icv] = A[cvocv_i[icv]]*r[icv];
        Ap[icv] = A[cvocv_i[icv]]*p[icv];
        for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          Ar[icv] += A[coc]*r[icv_nbr];
          Ap[icv] += A[coc]*p[icv_nbr];
        }
        Ar[icv] *= inv_diag[icv];
        Ap[icv] *= inv_diag[icv];
      }

      // calc alpha and beta...
      double my_buf[5]= {0.0,0.0,0.0,0.0,0.0};
      FOR_ICV {
        my_buf[0] += Ap[icv]*Ap[icv];
        my_buf[1] += Ap[icv]*Ar[icv];
        my_buf[2] += Ar[icv]*Ar[icv];
        my_buf[3] -=  r[icv]*Ap[icv];
        my_buf[4] -=  r[icv]*Ar[icv];
      }
      double buf[5];
      MPI_Allreduce((double*)my_buf,(double*)buf,5,MPI_DOUBLE,MPI_SUM,mpi_comm);
      double alpha,beta;
      const double inv_det = 1.0/(buf[0]*buf[2] - buf[1]*buf[1] + 1.0E-40);
      assert(inv_det > 0.0);
      alpha = relax1*inv_det*( buf[2]*buf[3] - buf[1]*buf[4]);
      beta  = relax2*inv_det*(-buf[1]*buf[3] + buf[0]*buf[4]);

      // phi += alpha*p + beta*r...
      // r += alpha*Ap + beta*Ar...
      FOR_ICV_G phi[icv] += alpha*p[icv]  + beta*r[icv];
      FOR_ICV r[icv] += alpha*Ap[icv] + beta*Ar[icv];

    }

  }

  // cleanup...
  delete[] r;
  delete[] p;
  delete[] inv_diag;
  delete[] Ar;
  delete[] Ap;

  // let the calling routine know if we were successful...
  return( done == 1 );

}

int StaticSolver::solveCvPatr(double (*u)[3],const double * A,const double * At, const double (*rhs)[3],
    const double zero,const double relax1,const double relax2,const int maxiter,const bool verbose) {
  // uses search directions p=A^T*res,res...

  // assume we come in with a consistent initial condition...

  // we need the following work arrays...
  double (*r)[3]   = new double[ncv_g][3];
  double (*p)[3]   = new double[ncv_g][3];
  double *inv_diag = new double[ncv_g];
  double (*Ar)[3]  = new double[ncv][3];
  double (*Ap)[3]  = new double[ncv][3];

  // normalization...
  FOR_ICV inv_diag[icv] = 1.0/A[cvocv_i[icv]];
  updateCvData(inv_diag);

  // residial: r = b-A*x...
  FOR_ICV {
    FOR_I3 r[icv][i] = A[cvocv_i[icv]]*u[icv][i] - rhs[icv][i];
    for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      FOR_I3 r[icv][i] += A[coc]*u[icv_nbr][i];
    }
    FOR_I3 r[icv][i] *= inv_diag[icv];
  }

  int iter = 0;
  int done = 0;
  while (done == 0) {

    // check if we are done...
    if (iter%10 == 0) {

      // compute the max (L-infinity) normalized residual...
      double  my_res_max = 0.0;
      FOR_ICV FOR_I3 my_res_max = max( my_res_max, fabs(r[icv][i]) );
      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank == 0) {
        // only share the last half of the convergence behaviour...
        if (verbose || (iter > maxiter/2))
          cout << " > solveCvPatr iter " << iter << " res_max " << res_max << endl;
        if (res_max <= zero) {
          //cout << "-> Successfully converged error to " << res_max << endl;
          done = 1;
        }
        else if (iter > maxiter) {
          cout << "Warning: solveCvPatr did not converge after " << maxiter <<
            " iters, res_max: " << res_max << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

    }

    // if we are not done, get the search direction and mag, and update guess...
    if (done == 0) {
      ++iter;

      // need residual in ghosts...
      updateCvData(r,REPLACE_ROTATE_DATA);

      // search direction: p = [A^T]{r}...
      FOR_ICV {
        FOR_I3 p[icv][i] = At[cvocv_i[icv]]*inv_diag[icv]*r[icv][i];
        for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          FOR_I3 p[icv][i] += At[coc]*inv_diag[icv_nbr]*r[icv_nbr][i];
        }
      }
      updateCvData(p,REPLACE_ROTATE_DATA);

      // second search direction is r...

      // get work arrays: Ap and Ar...
      FOR_ICV {
        FOR_I3 {
          Ar[icv][i] = A[cvocv_i[icv]]*r[icv][i];
          Ap[icv][i] = A[cvocv_i[icv]]*p[icv][i];
        }
        for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          FOR_I3 {
            Ar[icv][i] += A[coc]*r[icv_nbr][i];
            Ap[icv][i] += A[coc]*p[icv_nbr][i];
          }
        }
        FOR_I3 {
          Ar[icv][i] *= inv_diag[icv];
          Ap[icv][i] *= inv_diag[icv];
        }
      }

      // calc alpha and beta...
      //                      0   1   2
      double my_buf[5][3] = {{0.0,0.0,0.0},  // Ap^T*Ap
        {0.0,0.0,0.0},  // Ap^T*Ar
        {0.0,0.0,0.0},  // Ar^T*Ar
        {0.0,0.0,0.0},  //  r^T*Ap
        {0.0,0.0,0.0}}; //  r^T*Ar
      FOR_ICV {
        FOR_I3 {
          my_buf[0][i] += Ap[icv][i]*Ap[icv][i];
          my_buf[1][i] += Ap[icv][i]*Ar[icv][i];
          my_buf[2][i] += Ar[icv][i]*Ar[icv][i];
          my_buf[3][i] -=  r[icv][i]*Ap[icv][i];
          my_buf[4][i] -=  r[icv][i]*Ar[icv][i];
        }
      }
      double buf[5][3];
      MPI_Allreduce((double*)my_buf,(double*)buf,15,MPI_DOUBLE,MPI_SUM,mpi_comm);
      double alpha[3],beta[3];
      FOR_I3 {
        const double inv_det = 1.0/(buf[0][i]*buf[2][i] - buf[1][i]*buf[1][i] + 1.0E-40);
        assert(inv_det > 0.0);
        alpha[i] = relax1*inv_det*( buf[2][i]*buf[3][i] - buf[1][i]*buf[4][i]);
        beta[i]  = relax2*inv_det*(-buf[1][i]*buf[3][i] + buf[0][i]*buf[4][i]);
      }

      // u += alpha*p + beta*r...
      // r += alpha*Ap + beta*Ar...
      FOR_ICV_G FOR_I3 u[icv][i] += alpha[i]*p[icv][i]  + beta[i]*r[icv][i];
      FOR_ICV FOR_I3 r[icv][i] += alpha[i]*Ap[icv][i] + beta[i]*Ar[icv][i];

    }

  }

  // cleanup...
  delete[] r;
  delete[] p;
  delete[] inv_diag;
  delete[] Ar;
  delete[] Ap;

  // let the calling routine know if we were successful...
  return( done == 1 );

}

int StaticSolver::solveCvPatr(double (*u)[3],const double * A,const double * At, const double (*A_diag)[3], const double (*rhs)[3],
    const double zero,const double relax1,const double relax2,const int maxiter,const bool verbose) {
  // uses search directions p=A^T*res,res...

  // assume we come in with a consistent initial condition...

  // we need the following work arrays...
  double (*r)[3]   = new double[ncv_g][3];
  double (*p)[3]   = new double[ncv_g][3];
  double (*inv_diag)[3] = new double[ncv_g][3];
  double (*Ar)[3]  = new double[ncv][3];
  double (*Ap)[3]  = new double[ncv][3];

  // normalization...
  FOR_ICV FOR_I3 inv_diag[icv][i] = 1.0/(A[cvocv_i[icv]]+A_diag[icv][i]);
  updateCvData(inv_diag);

  // residial: r = b-A*x...
  FOR_ICV {
    FOR_I3 r[icv][i] = (A[cvocv_i[icv]]+A_diag[icv][i])*u[icv][i] - rhs[icv][i];
    for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
      const int icv_nbr = cvocv_v[coc];
      FOR_I3 r[icv][i] += A[coc]*u[icv_nbr][i];
    }
    FOR_I3 r[icv][i] *= inv_diag[icv][i];
  }

  int iter = 0;
  int done = 0;
  while (done == 0) {

    // check if we are done...
    if (iter%10 == 0) {

      // compute the max (L-infinity) normalized residual...
      double  my_res_max = 0.0;
      FOR_ICV FOR_I3 my_res_max = max( my_res_max, fabs(r[icv][i]) );
      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank == 0) {
        // only share the last half of the convergence behaviour...
        if (verbose || (iter > maxiter/2))
          cout << " > solveCvPatr iter " << iter << " res_max " << res_max << endl;
        if (res_max <= zero) {
          //cout << "-> Successfully converged error to " << res_max << endl;
          done = 1;
        }
        else if (iter > maxiter) {
          cout << "Warning: solveCvPatr did not converge after " << maxiter <<
            " iters, res_max: " << res_max << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

    }

    // if we are not done, get the search direction and mag, and update guess...
    if (done == 0) {
      ++iter;

      // need residual in ghosts...
      updateCvData(r,REPLACE_ROTATE_DATA);

      // search direction: p = [A^T]{r}...
      FOR_ICV {
        FOR_I3 p[icv][i] = (At[cvocv_i[icv]]+A_diag[icv][i])*inv_diag[icv][i]*r[icv][i];
        for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          FOR_I3 p[icv][i] += At[coc]*inv_diag[icv_nbr][i]*r[icv_nbr][i];
        }
      }
      updateCvData(p,REPLACE_ROTATE_DATA);

      // second search direction is r...

      // get work arrays: Ap and Ar...
      FOR_ICV {
        FOR_I3 {
          Ar[icv][i] = (A[cvocv_i[icv]]+A_diag[icv][i])*r[icv][i];
          Ap[icv][i] = (A[cvocv_i[icv]]+A_diag[icv][i])*p[icv][i];
        }
        for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          FOR_I3 {
            Ar[icv][i] += A[coc]*r[icv_nbr][i];
            Ap[icv][i] += A[coc]*p[icv_nbr][i];
          }
        }
        FOR_I3 {
          Ar[icv][i] *= inv_diag[icv][i];
          Ap[icv][i] *= inv_diag[icv][i];
        }
      }

      // calc alpha and beta...
      //                      0   1   2
      double my_buf[5][3] = {{0.0,0.0,0.0},  // Ap^T*Ap
        {0.0,0.0,0.0},  // Ap^T*Ar
        {0.0,0.0,0.0},  // Ar^T*Ar
        {0.0,0.0,0.0},  //  r^T*Ap
        {0.0,0.0,0.0}}; //  r^T*Ar
      FOR_ICV {
        FOR_I3 {
          my_buf[0][i] += Ap[icv][i]*Ap[icv][i];
          my_buf[1][i] += Ap[icv][i]*Ar[icv][i];
          my_buf[2][i] += Ar[icv][i]*Ar[icv][i];
          my_buf[3][i] -=  r[icv][i]*Ap[icv][i];
          my_buf[4][i] -=  r[icv][i]*Ar[icv][i];
        }
      }
      double buf[5][3];
      MPI_Allreduce((double*)my_buf,(double*)buf,15,MPI_DOUBLE,MPI_SUM,mpi_comm);
      double alpha[3],beta[3];
      FOR_I3 {
        const double inv_det = 1.0/(buf[0][i]*buf[2][i] - buf[1][i]*buf[1][i] + 1.0E-40);
        assert(inv_det > 0.0);
        alpha[i] = relax1*inv_det*( buf[2][i]*buf[3][i] - buf[1][i]*buf[4][i]);
        beta[i]  = relax2*inv_det*(-buf[1][i]*buf[3][i] + buf[0][i]*buf[4][i]);
      }

      // u += alpha*p + beta*r...
      // r += alpha*Ap + beta*Ar...
      FOR_ICV_G FOR_I3 u[icv][i] += alpha[i]*p[icv][i]  + beta[i]*r[icv][i];
      FOR_ICV FOR_I3 r[icv][i] += alpha[i]*Ap[icv][i] + beta[i]*Ar[icv][i];

    }

  }

  // cleanup...
  delete[] r;
  delete[] p;
  delete[] inv_diag;
  delete[] Ar;
  delete[] Ap;

  // let the calling routine know if we were successful...
  return( done == 1 );

}

int StaticSolver::solveCvTim(double (*u)[3],const double *A,const double *As,const double (*rhs)[3],
    const double zero,const double tau_bound,const double relax,const int maxiter,const bool verbose) {

  // NOTE: assumes diag(A) = 1 (i.e. A is actually DinvA). As is the skew part of DinvA

  double (*res)[3] = new double[ncv][3];

  const double tau = relax*tau_bound;
  double (*w)[3] = new double[ncv][3];

  int iter = 0;
  int done = 0;
  while (done == 0) {

    iter++;

    // calculate the residual...
    calcCvResidual(res,u,A,rhs);

    // check...
    if ( iter%3 == 0) {
      double my_res_max = 0.0;
      FOR_ICV FOR_I3 my_res_max = max(my_res_max,fabs(res[icv][i]));
      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank == 0) {
        if ((verbose)||(iter > maxiter/2))
          cout << " > solveCvTim iter, res_max: " << iter << " " << res_max << endl;
        if (res_max < zero) {
          done = 1;
        }
        else if (iter > maxiter) {
          cout << " > Warning: solveCvTim did not converge after " << maxiter << " iters, res_max: " << res_max << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
    }

    if (done == 0) {
      FOR_ICV {
        FOR_I3 w[icv][i] = u[icv][i] + tau*(rhs[icv][i] - A[cvocv_i[icv]]*u[icv][i]);
        for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          if (icv_nbr > icv) {
            FOR_I3 w[icv][i] -= tau*A[coc]*u[icv_nbr][i];
          }
          else {
            assert(icv_nbr < icv);
            FOR_I3 w[icv][i] -= tau*((A[coc]-2*As[coc])*u[icv_nbr][i] + 2.0*As[coc]*w[icv_nbr][i]);
          }
        }
      }
      FOR_ICV FOR_I3 u[icv][i] = w[icv][i];
      updateCvData(u,REPLACE_ROTATE_DATA);
    }

  }

  delete[] res;
  delete[] w;

  return( done == 1 );

}

void StaticSolver::calcCvGradLeastSquares(double (*dphidx)[3],double * phi) {

  // computes the cv-base gradient assuming the face-based gradients
  // in the NORMAL direction are given by
  //
  // dpdn = (phi[icv1] - phi[icv0])*area_over_delta_fa[ifa]/area_fa;
  //

  // phi should be valid in the ghost cvs...

  double (*sAninj_diag)[3] = new double[ncv][3];
  double (*sAninj_offd)[3] = new double[ncv][3];
  double (*sAnidpdn)[3] = new double[ncv][3];

  FOR_ICV {
    FOR_I3 sAninj_diag[icv][i] = 0.0;
    FOR_I3 sAninj_offd[icv][i] = 0.0;
    FOR_I3 sAnidpdn[icv][i] = 0.0;
  }

  FOR_IFA {
    const double area_fa = MAG(n_fa[ifa]);
    if (area_fa > 0.0) {
      const double unit_n[3] = {
        n_fa[ifa][0]/area_fa,
        n_fa[ifa][1]/area_fa,
        n_fa[ifa][2]/area_fa
      };
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      FOR_I3 sAninj_diag[icv0][i] += area_fa*unit_n[i]*unit_n[i];
      FOR_I3 sAninj_offd[icv0][i] += area_fa*unit_n[(i+1)%3]*unit_n[(i+2)%3];
      const double dpdn = (phi[icv1] - phi[icv0])*area_over_delta_fa[ifa]/area_fa;
      FOR_I3 sAnidpdn[icv0][i] += area_fa*unit_n[i]*dpdn;
      if (icv1 < ncv) {
        // sign flips, but everything contains 2 terms, so disregard...
        FOR_I3 sAninj_diag[icv1][i] += area_fa*unit_n[i]*unit_n[i];
        FOR_I3 sAninj_offd[icv1][i] += area_fa*unit_n[(i+1)%3]*unit_n[(i+2)%3];
        FOR_I3 sAnidpdn[icv1][i] += area_fa*unit_n[i]*dpdn;
      }
    }
  }

  // for boundary faces, assume zero gradient...

  FOR_IBF {
    const double area_bf = MAG(n_bf[ibf]);
    if (area_bf > 0.0) {
      const double unit_n[3] = {
        n_bf[ibf][0]/area_bf,
        n_bf[ibf][1]/area_bf,
        n_bf[ibf][2]/area_bf
      };
      const int icv = cvobf[ibf];
      FOR_I3 sAninj_diag[icv][i] += area_bf*unit_n[i]*unit_n[i];
      FOR_I3 sAninj_offd[icv][i] += area_bf*unit_n[(i+1)%3]*unit_n[(i+2)%3];
      // assume dpdn == 0, so nothing to add...
    }
  }

  // now invert to build the gradient...

  FOR_ICV {

    const double denom =
      sAninj_diag[icv][0]*sAninj_diag[icv][1]*sAninj_diag[icv][2] +
      2.0*sAninj_offd[icv][0]*sAninj_offd[icv][1]*sAninj_offd[icv][2] -
      sAninj_diag[icv][0]*sAninj_offd[icv][0]*sAninj_offd[icv][0] -
      sAninj_diag[icv][1]*sAninj_offd[icv][1]*sAninj_offd[icv][1] -
      sAninj_diag[icv][2]*sAninj_offd[icv][2]*sAninj_offd[icv][2];
    assert(denom != 0.0);

    dphidx[icv][0] =
      ( (sAninj_diag[icv][1]*sAninj_diag[icv][2]-sAninj_offd[icv][0]*sAninj_offd[icv][0])*sAnidpdn[icv][0] +
        (sAninj_offd[icv][0]*sAninj_offd[icv][1]-sAninj_offd[icv][2]*sAninj_diag[icv][2])*sAnidpdn[icv][1] +
        (sAninj_offd[icv][0]*sAninj_offd[icv][2]-sAninj_offd[icv][1]*sAninj_diag[icv][1])*sAnidpdn[icv][2] )/denom;

    dphidx[icv][1] =
      ( (sAninj_diag[icv][2]*sAninj_diag[icv][0]-sAninj_offd[icv][1]*sAninj_offd[icv][1])*sAnidpdn[icv][1] +
        (sAninj_offd[icv][1]*sAninj_offd[icv][2]-sAninj_offd[icv][0]*sAninj_diag[icv][0])*sAnidpdn[icv][2] +
        (sAninj_offd[icv][1]*sAninj_offd[icv][0]-sAninj_offd[icv][2]*sAninj_diag[icv][2])*sAnidpdn[icv][0] )/denom;

    dphidx[icv][2] =
      ( (sAninj_diag[icv][0]*sAninj_diag[icv][1]-sAninj_offd[icv][2]*sAninj_offd[icv][2])*sAnidpdn[icv][2] +
        (sAninj_offd[icv][2]*sAninj_offd[icv][0]-sAninj_offd[icv][1]*sAninj_diag[icv][1])*sAnidpdn[icv][0] +
        (sAninj_offd[icv][2]*sAninj_offd[icv][1]-sAninj_offd[icv][0]*sAninj_diag[icv][0])*sAnidpdn[icv][1] )/denom;
  }

  delete[] sAninj_diag;
  delete[] sAninj_offd;
  delete[] sAnidpdn;

}

void StaticSolver::calcCvGradLeastSquares(double (*dphidx)[3],double * phi,double * minus_dphidn_dA_bf) {

  // here we are using the relationship rhou = -dphidn

  // computes the cv-base gradient assuming the face-based gradients
  // in the NORMAL direction are given by
  //
  // dpdn = (phi[icv1] - phi[icv0])*area_over_delta_fa[ifa]/area_fa;
  //

  // phi should be valid in the ghost cvs...

  double (*sAninj_diag)[3] = new double[ncv][3];
  double (*sAninj_offd)[3] = new double[ncv][3];
  double (*sAnidpdn)[3] = new double[ncv][3];

  FOR_ICV {
    FOR_I3 sAninj_diag[icv][i] = 0.0;
    FOR_I3 sAninj_offd[icv][i] = 0.0;
    FOR_I3 sAnidpdn[icv][i] = 0.0;
  }

  FOR_IFA {
    const double area_fa = MAG(n_fa[ifa]);
    if (area_fa > 0.0) {
      const double unit_n[3] = {
        n_fa[ifa][0]/area_fa,
        n_fa[ifa][1]/area_fa,
        n_fa[ifa][2]/area_fa
      };
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      FOR_I3 sAninj_diag[icv0][i] += area_fa*unit_n[i]*unit_n[i];
      FOR_I3 sAninj_offd[icv0][i] += area_fa*unit_n[(i+1)%3]*unit_n[(i+2)%3];
      const double dpdn_dA = (phi[icv1] - phi[icv0])*area_over_delta_fa[ifa];
      FOR_I3 sAnidpdn[icv0][i] += unit_n[i]*dpdn_dA;
      if (icv1 < ncv) {
        // sign flips, but everything contains 2 terms, so disregard...
        FOR_I3 sAninj_diag[icv1][i] += area_fa*unit_n[i]*unit_n[i];
        FOR_I3 sAninj_offd[icv1][i] += area_fa*unit_n[(i+1)%3]*unit_n[(i+2)%3];
        FOR_I3 sAnidpdn[icv1][i] += unit_n[i]*dpdn_dA;
      }
    }
  }

  // for boundary faces, use passed normal gradient to close...

  FOR_IBF {
    const double area_bf = MAG(n_bf[ibf]);
    if (area_bf > 0.0) {
      const double unit_n[3] = {
        n_bf[ibf][0]/area_bf,
        n_bf[ibf][1]/area_bf,
        n_bf[ibf][2]/area_bf
      };
      const int icv = cvobf[ibf];
      FOR_I3 sAninj_diag[icv][i] += area_bf*unit_n[i]*unit_n[i];
      FOR_I3 sAninj_offd[icv][i] += area_bf*unit_n[(i+1)%3]*unit_n[(i+2)%3];
      FOR_I3 sAnidpdn[icv][i]    -= unit_n[i]*minus_dphidn_dA_bf[ibf];
    }
  }

  // now invert to build the gradient...

  FOR_ICV {

    const double denom =
      sAninj_diag[icv][0]*sAninj_diag[icv][1]*sAninj_diag[icv][2] +
      2.0*sAninj_offd[icv][0]*sAninj_offd[icv][1]*sAninj_offd[icv][2] -
      sAninj_diag[icv][0]*sAninj_offd[icv][0]*sAninj_offd[icv][0] -
      sAninj_diag[icv][1]*sAninj_offd[icv][1]*sAninj_offd[icv][1] -
      sAninj_diag[icv][2]*sAninj_offd[icv][2]*sAninj_offd[icv][2];
    if (denom == 0.0)
      cout << mpi_rank << " about to fail calcCvGradLeastSquares: " << COUT_VEC(x_cv[icv]) << endl;
    assert(denom != 0.0);

    dphidx[icv][0] =
      ( (sAninj_diag[icv][1]*sAninj_diag[icv][2]-sAninj_offd[icv][0]*sAninj_offd[icv][0])*sAnidpdn[icv][0] +
        (sAninj_offd[icv][0]*sAninj_offd[icv][1]-sAninj_offd[icv][2]*sAninj_diag[icv][2])*sAnidpdn[icv][1] +
        (sAninj_offd[icv][0]*sAninj_offd[icv][2]-sAninj_offd[icv][1]*sAninj_diag[icv][1])*sAnidpdn[icv][2] )/denom;

    dphidx[icv][1] =
      ( (sAninj_diag[icv][2]*sAninj_diag[icv][0]-sAninj_offd[icv][1]*sAninj_offd[icv][1])*sAnidpdn[icv][1] +
        (sAninj_offd[icv][1]*sAninj_offd[icv][2]-sAninj_offd[icv][0]*sAninj_diag[icv][0])*sAnidpdn[icv][2] +
        (sAninj_offd[icv][1]*sAninj_offd[icv][0]-sAninj_offd[icv][2]*sAninj_diag[icv][2])*sAnidpdn[icv][0] )/denom;

    dphidx[icv][2] =
      ( (sAninj_diag[icv][0]*sAninj_diag[icv][1]-sAninj_offd[icv][2]*sAninj_offd[icv][2])*sAnidpdn[icv][2] +
        (sAninj_offd[icv][2]*sAninj_offd[icv][0]-sAninj_offd[icv][1]*sAninj_diag[icv][1])*sAnidpdn[icv][0] +
        (sAninj_offd[icv][2]*sAninj_offd[icv][1]-sAninj_offd[icv][0]*sAninj_diag[icv][0])*sAnidpdn[icv][1] )/denom;
  }

  delete[] sAninj_diag;
  delete[] sAninj_offd;
  delete[] sAnidpdn;

}

void StaticSolver::setupFilterMoments() { 

  // note that x_vv is available to you in ncv_g but not ncv_g2

  assert( faocv_i);
  assert( faocv_v);
  assert( bfocv_i);
  assert( bfocv_v);
  assert( filter_mom2 == NULL);

  filter_mom2 = new double[ncv];

  for (int icv = 0; icv < ncv; ++icv) { 

    // approximate the tensor variables with the scalar representation 
    // from a spherical ball filter (mom2 in theory could be expanded 
    // to include all 9 tensorial components ... ) 

    filter_mom2[icv]  = 0.0;
    double wgt        = 0.0;

    for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) { 

      const int ifa     = faocv_v[foc];
      const int ii      = foc - faocv_i[icv];

      // cvocv has a diagonal entry, so this needs shifted.

      const int icv_nbr = cvocv_v[cvocv_i[icv]+ii+1]; 

      if ( cvofa[ifa][0] == icv ) { 

        assert( cvofa[ifa][1] == icv_nbr);

      } else { 

        assert( cvofa[ifa][1] == icv);
        assert( cvofa[ifa][0] == icv_nbr);

      }

      const double area    = MAG(n_fa[ifa]);
      const double delta   = DIST(x_vv[icv_nbr],x_vv[icv]);
      const double vol_    = area*delta/6.0;

      filter_mom2[icv]    += vol_ * delta * delta; 
      wgt                 += vol_ ; 

    }


    for (int bof = bfocv_i[icv]; bof != bfocv_i[icv+1]; ++bof) { 

      const int ibf = bfocv_v[bof];
      const double area  = area_bf[ibf]; 
      const double delta = area_bf[ibf]/area_over_delta_bf[ibf];

      const double vol_  = delta*area/3.0;
      filter_mom2[icv]  += vol_ * delta * delta; 
      wgt               += vol_; 

    } 


    // normalized mom2... 

    const double del_cv = pow(vol_cv[icv], 1.0/3.0);
    filter_mom2[icv]   /= (wgt * del_cv * del_cv ) ;


  }

  MiscUtils::dumpRange(filter_mom2, ncv, "filter mom2"); 

} 

void StaticSolver::filterCvR1(double* phif, const double* phi) {

  // the following routine assumes that the ghost values have
  // already been synced in phi.  requires that the first level
  // ghosts are populated.

  // the weight here should correspond to the volume of the cells on
  // the interior cells and for boundaries where the boundary face is
  // indeed planar.

  double * wgt = new double[ncv];

  for (int icv = 0; icv < ncv; ++icv) {
    wgt[icv]  = 0.0;
    phif[icv] = 0.0;
  }

  for (int ifa = 0; ifa < nfa; ++ifa) {

    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];

    // consider precomputing the following...
    const double area    = MAG(n_fa[ifa]);
    const double delta   = DIST(x_vv[icv1],x_vv[icv0]);

    const double phi_avg = 0.5*(phi[icv0] + phi[icv1]);
    const double vol_    = area*delta/6.0;
    const double c_avg   = phi_avg*vol_;

    phif[icv0]          += c_avg;
    wgt[icv0]           += vol_;
    if ( icv1 < ncv) {
      phif[icv1]        += c_avg;
      wgt[icv1]         += vol_;
    }

  }

  for (int ibf = 0; ibf < nbf; ++ibf) {

    const int icv      = cvobf[ibf];
    //const double dx[3] = DIFF(x_vv[icv],x_bf[ibf]);
    //const double delta = abs(DOT_PRODUCT(dx,n_bf[ibf]))/MAG(n_bf[ibf]);
    const double area  = area_bf[ibf];
    const double delta = area_bf[ibf]/area_over_delta_bf[ibf];

    const double vol_  = delta*area/3.0;
    phif[icv]         += phi[icv]*vol_;
    wgt[icv]          += vol_;
  }


  for (int icv = 0; icv < ncv; ++icv) {
    //phif[icv] /= vol_cv[icv];
    phif[icv]   /= wgt[icv];
  }


  // check the consistency of the weight with the volume on interior cells...
  /*
     for (int icv = 0 ; icv < ncv; ++icv) {
     wgt[icv] -= vol_cv[icv];
     wgt[icv] /= vol_cv[icv];
     }

     MiscUtils::dumpRange(wgt,ncv, "wgt vol error");

     for (int ibf = 0; ibf < nbf; ++ibf)
     wgt[cvobf[ibf]] = 0.0;

     MiscUtils::dumpRange(wgt,ncv, "wgt vol error (internal only)");
     */

  delete[] wgt;

}

void StaticSolver::filterCvR2(double (*phif)[3], const double (*phi)[3]) {

  double * wgt = new double[ncv];

  for (int icv = 0; icv < ncv; ++icv) {
    wgt[icv] = 0.0;
    for (int i = 0; i < 3; ++i)
      phif[icv][i] = 0.0;
  }

  for (int ifa = 0; ifa < nfa; ++ifa) {

    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];

    // consider precomputing the following...
    const double area    = MAG(n_fa[ifa]);
    const double delta   = DIST(x_vv[icv1],x_vv[icv0]);

    double c_avg[3];
    for (int i = 0; i < 3; ++i)
      c_avg[i] = 0.5*(phi[icv0][i] + phi[icv1][i]);

    const double fax = area*delta/6.0;

    wgt[icv0]       += fax;
    for (int i = 0; i < 3; ++i)
      phif[icv0][i] += c_avg[i]*fax;

    if ( icv1 < ncv) {
      wgt[icv1]     += fax;
      for (int i = 0; i < 3; ++i)
        phif[icv1][i] += c_avg[i]*fax;
    }
  }

  for (int ibf = 0; ibf < nbf; ++ibf) {

    const int icv      = cvobf[ibf];
    //const double dx[3] = DIFF(x_vv[icv],x_bf[ibf]);
    //const double area  = MAG(n_bf[ibf]);
    //const double delta = abs(DOT_PRODUCT(dx,n_bf[ibf]))/area;

    const double area    = area_bf[ibf];
    const double delta   = area_bf[ibf]/area_over_delta_bf[ibf];

    const double fax = area*delta/3.0;
    wgt[icv]        += fax;
    for (int i = 0; i < 3; ++i)
      phif[icv][i]  += phi[icv][i]*fax;
  }


  for (int icv = 0; icv < ncv; ++icv) {
    for (int i = 0; i < 3; ++i) {
      //phif[icv][i] /= vol_cv[icv];
      phif[icv][i] /= wgt[icv];
    }
  }

  delete[] wgt;

}

int StaticSolver::pointIsInside(const double xp[3], int &icv_ret) {

  icv_ret = -1;

  // returns:
  // 1: is inside   icv_ret: icv_closest
  // 0: is outside  icv_ret: -1
  ensureCvAdt();

  vector<int> cvList;
  cvAdt->buildListForPoint(cvList, xp);

  int icv_closest = -1;
  double d2_closest;
  for (int ivec = 0; ivec < cvList.size(); ivec++) {
    const int icv = cvList[ivec];
    const double d2 = DIST2(x_vv[icv], xp);
    if ((icv_closest==-1)||(d2<d2_closest)) {
      d2_closest = d2;
      icv_closest = icv;
    }
  }
  if (icv_closest==-1) {
    // could not find a cell to own the point
    return 0;
  }

  assert((icv_closest>=0)&&(icv_closest<ncv));
  // find the closest neighbor cell to the point
  int icv_nbr_closest = -1;
  double d2_nbr_closest;
  for (int coc = cvocv_i[icv_closest]+1; coc != cvocv_i[icv_closest+1]; coc++) {
    const int icv_nbr = cvocv_v[coc];
    assert(icv_nbr!=icv_closest);
    const double d2_nbr = DIST2(x_vv[icv_nbr], xp);
    if ((icv_nbr_closest==-1)||(d2_nbr < d2_nbr_closest)) {
      d2_nbr_closest = d2_nbr;
      icv_nbr_closest = icv_nbr;
    }
  }

  if (d2_nbr_closest < d2_closest) {
    // we got a closer neighbor
    return 0;
  }

  // now we should have icv_closest with the point inside of it
  assert((icv_closest>=0)&&(icv_closest<ncv));

  // check if the closest cell has a boundary surface and decide if it is inside
  const int my_nbf = bfocv_i[icv_closest+1] - bfocv_i[icv_closest];
  if (my_nbf==0) {
    // no boundary faces so the point is inside
    icv_ret = icv_closest;
    return 1;
  }
  else {

    // and we need the closest ibf,ist_ss...
    // We need a length scale...
    const double tol2 = 1.0E-10*pow(vol_cv[icv_closest],2.0/3.0);
    int ibf_closest = -1;
    int ist_ss_closest;
    double d2_closest,dist_closest;
    //double dx_closest[3],normal_closest[3];
    for (int boc = bfocv_i[icv_closest]; boc != bfocv_i[icv_closest+1]; ++boc) {
      const int ibf = bfocv_v[boc];
      for (int sob = sstobf_i[ibf]; sob != sstobf_i[ibf+1]; ++sob) {
        const int ist_ss = sstobf_v[sob];

        // the corner coords of this sub-surface tri are (recall that the
        // subsurface tri already has any periodic transfrom taken care of)...
        const double * const x0 = subSurface->xp[subSurface->spost[ist_ss][0]];
        const double * const x1 = subSurface->xp[subSurface->spost[ist_ss][1]];
        const double * const x2 = subSurface->xp[subSurface->spost[ist_ss][2]];

        double xp_tri[3]; MiscUtils::getClosestPointOnTriRobust(xp_tri,xp,x0,x1,x2);

        const double dx[3] = DIFF(xp_tri, xp);
        const double d2 = DOT_PRODUCT(dx,dx);

        if ((ibf_closest == -1)||(d2 < d2_closest+tol2)) {
          const double normal[3] = TRI_NORMAL_2(x0,x1,x2);
          const double normal_mag = MAG(normal);
          assert(normal_mag != 0); // the subsurface should not have zero-area tris!
          const double dist = DOT_PRODUCT(dx,normal)/normal_mag;
          if ((ibf_closest == -1)||(d2 < d2_closest-tol2)||(fabs(dist) > fabs(dist_closest))) {
            ibf_closest = ibf;
            ist_ss_closest = ist_ss;
            d2_closest = d2;
            //FOR_I3 dx_closest[i] = dx[i];
            //FOR_I3 normal_closest[i] = normal[i];
            dist_closest = dist;
          }
        }
      }
    }
    assert(ist_ss_closest >= 0);
    if (dist_closest >= 0.0) {
      // dist with positive sign means we are inside the fluid volume...
      icv_ret = icv_closest;
      return 1;
    } else {
      // dist with negative sign means we are outside the fluid volume...
      return 0;
    }
  }

  // we should never reach here
  assert(0);
  return 0;
}

int StaticSolver::initPlaneInjector(double * &cdf_area_rank, double * &cdf_area_tri, vector<pair<SimpleTri,int> >& triVec, const double * const xp_plane, const double * const np_plane) {

  triVec.clear();

  double e1[3], e2[3];
  MiscUtils::getBestE1E2FromE0(e1, e2, np_plane);

  IF_RANK0 cout << "   > injector GEOM = PLANE xp: " << COUT_VEC(xp_plane) << " np: " << COUT_VEC(np_plane) << endl;

  double * iso_var_no = new double[nno];

  for (int ino = 0; ino < nno; ++ino) {
    iso_var_no[ino] =
      (x_no[ino][0]-xp_plane[0])*np_plane[0] +
      (x_no[ino][1]-xp_plane[1])*np_plane[1] +
      (x_no[ino][2]-xp_plane[2])*np_plane[2];
  }

  buildIsoWithIcv(triVec,iso_var_no,0.0);

  delete[] iso_var_no;

  // check if any triangles are selected
  int my_ntri = triVec.size();
  int ntri_global;
  MPI_Allreduce(&my_ntri, &ntri_global, 1, MPI_INT, MPI_SUM, mpi_comm);
  if (ntri_global == 0) return 1;

  //GeomUtils::writeTecplot("triVec.dat",triVec);

  // put the selected cvs in groups based on their participation in trivec...

  // start by putting every cv in its own group. Recall the cv associated with
  // every tri in the triVec is stored in the second...
  int * cv_flag = new int[ncv_g];
  for (int icv = 0; icv < ncv; ++icv) cv_flag[icv] = -1;
  int ngr = 0;
  for (int ii = 0; ii < triVec.size(); ++ii) {
    const int icv = triVec[ii].second;
    assert((icv >= 0)&&(icv < ncv));
    if (cv_flag[icv] == -1) cv_flag[icv] = ngr++;
  }

  // offset the cv_flag so we are all globally unique...
  int offset;
  MPI_Scan(&ngr,&offset,1,MPI_INT,MPI_SUM,mpi_comm);
  offset -= ngr; assert(offset >= 0);

  int icv_closest = -1;
  double d2_closest = HUGE_VAL;
  for (int icv = 0; icv < ncv; ++icv) {
    if (cv_flag[icv] >= 0) {
      cv_flag[icv] += offset;
      const double d2 = DIST2(x_vv[icv],xp_plane);
      if (d2 < d2_closest) {
        icv_closest = icv;
        d2_closest = d2;
      }
    }
  }

  int done = 0;
  while (done == 0) {

    int my_done = 1;

    updateCvData(cv_flag);

    // TODO: this makes no use of ordering and internal-only, so 100's of iterations
    // may be necessary for large problems...
    // e.g. loop through interproc faces and see if any local cvs need to be updated,etc...

    for (int icv = 0; icv < ncv; ++icv) {
      if (cv_flag[icv] >= 0) {
        for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          if ((cv_flag[icv_nbr] >= 0)&&(cv_flag[icv_nbr] < cv_flag[icv])) {
            cv_flag[icv] = cv_flag[icv_nbr];
            my_done = 0;
          }
        }
      }
    }

    MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);

  }

  // at this point the group indexing will be unique, but not continuous.
  // this does not matter. We just want to save the group of the
  // closest icv...

  DoubleInt my_di, di;
  my_di.this_double = d2_closest;
  my_di.this_int    = mpi_rank;
  MPI_Allreduce(&my_di,&di,1,MPI_DOUBLE_INT,MPI_MINLOC, mpi_comm);

  int igr_closest;
  if (di.this_int == mpi_rank) {
    assert(icv_closest >= 0);
    igr_closest = cv_flag[icv_closest];
    assert(igr_closest >= 0);
  }
  MPI_Bcast(&igr_closest,1,MPI_INT,di.this_int,mpi_comm);

  // clean the triVec...
  int ii_new = 0;
  for (int ii = 0; ii < triVec.size(); ++ii) {
    const int icv = triVec[ii].second;
    if (cv_flag[icv] == igr_closest) {
      if (ii_new != ii) {
        FOR_I3 triVec[ii_new].first.x0[i] = triVec[ii].first.x0[i];
        FOR_I3 triVec[ii_new].first.x1[i] = triVec[ii].first.x1[i];
        FOR_I3 triVec[ii_new].first.x2[i] = triVec[ii].first.x2[i];
        triVec[ii_new].second = icv;
      }
      ++ii_new;
    }
  }
  triVec.resize(ii_new);

  delete[] cv_flag;

  //GeomUtils::writeTecplot("triVec2.dat",triVec);

  // =====================================
  // also build a cdf of tris and ranks...
  // =====================================

  const double np_mag = MAG(np_plane); assert(np_mag > 0.0); assert(fabs(np_mag-1.0)<1e-6);
  assert(cdf_area_tri==NULL);
  cdf_area_tri = new double [triVec.size()];
  double my_area = 0.0;
  for (int ii = 0; ii < triVec.size(); ++ii) {
    const double this_n[3] = TRI_NORMAL_2(triVec[ii].first.x0,triVec[ii].first.x1,triVec[ii].first.x2);
    const double area = 0.5*DOT_PRODUCT(this_n,np_plane);
    my_area += area;
    cdf_area_tri[ii] = fabs(area);
    if (ii > 0) cdf_area_tri[ii] += cdf_area_tri[ii-1];
  }

  for (int ii = 0; ii < triVec.size(); ++ii) {
    cdf_area_tri[ii] /= cdf_area_tri[triVec.size()-1];
  }

  assert(cdf_area_rank == NULL);
  if (mpi_rank == 0) {
    cdf_area_rank = new double[mpi_size];
  }
  MPI_Gather(&my_area,1,MPI_DOUBLE,cdf_area_rank,1,MPI_DOUBLE,0,mpi_comm);
  if (mpi_rank == 0) {
    double area_sum = 0.0;
    FOR_RANK {
      area_sum += cdf_area_rank[rank];
      cdf_area_rank[rank] = area_sum;
    }
    assert(area_sum>0);
    FOR_RANK cdf_area_rank[rank] /= area_sum;
    cout << "   > plane area: " << area_sum << endl;
  }

  //IF_RANK0 {
  //  FOR_RANK {
  //    cout << rank << " cdf_area_rank: " << cdf_area_rank[rank] << endl;
  //  }
  //}
  //MPI_Pause("HERE");

  return 0;
}

int StaticSolver::initZoneInjector(double * &cdf_area_tri, int (* &spost_zn)[3] , double (* &xp_zn)[3], int& nst_zn, int& nsp_zn, const string zone_name) {

  map<const string, int>::iterator iter = bfZoneNameMap.find(zone_name);
  if (iter == bfZoneNameMap.end()) {
    return 1;
  }

  int zone_ind = bfZoneNameMap[zone_name];
  if (mpi_rank == 0) {
    cout << "   > injector GEOM = ZONE name: " << zone_name << " id (bfZone index): " << zone_ind << " " << endl;
  }

  int my_nst = 0; // number of triangles in this zone
  for (int ist = 0; ist < subSurface->nst; ist++) {
    if (subSurface->znost[ist] == zone_ind)
      my_nst++;
  }

  int (* my_spost)[3] = new int [my_nst][3];
  map<int,int> spMap; // a map from subsurface isp to the zone isp
  int my_nsp = 0; // number of distinct points
  int ist = 0;
  if (my_nst > 0) {
    for (int ist_ss = 0; ist_ss < subSurface->nst; ist_ss++) {
      if (subSurface->znost[ist_ss] == zone_ind) {
        FOR_I3 {
          int isp_ss = subSurface->spost[ist][i];
          map<int,int>::iterator it = spMap.find(isp_ss);
          if (it == spMap.end()) { // ip_old is a new point, add it to the map
            spMap[isp_ss] = my_nsp++;
          }
          my_spost[ist][i] = spMap[isp_ss];
        }
        ist++;
      }
    }
  }
  assert(ist==my_nst);

  double (* my_xp)[3] = new double [my_nsp][3];
  for (map<int,int>::iterator it = spMap.begin(); it!=spMap.end(); it++) {
    int isp_ss = it->first;
    int isp = it->second;
    FOR_I3 my_xp[isp][i] = subSurface->xp[isp_ss][i];
  }

  // communicate nst, nsp
  int * nst_ranks = new int [mpi_size]; // nst for each processor
  int * nsp_ranks = new int [mpi_size]; // nsp for each processor
  MPI_Allgather(&my_nst,1,MPI_INT,nst_ranks,1,MPI_INT,mpi_comm);
  MPI_Allgather(&my_nsp,1,MPI_INT,nsp_ranks,1,MPI_INT,mpi_comm);

  nst_zn = 0;
  nsp_zn = 0;
  for (int iproc = 0; iproc < mpi_size; iproc++) {
    nst_zn += nst_ranks[iproc];
    nsp_zn += nsp_ranks[iproc];
  }

  // communicate spost
  // correct my_spost, shift the isp to accomodate for other processors
  int shift = 0;
  for (int iproc = 0; iproc < mpi_rank; iproc++)
    shift += nsp_ranks[iproc];
  for (int ist = 0; ist < my_nst; ist++)
    FOR_I3 my_spost[ist][i] += shift;

  int * disp = new int [mpi_size];
  disp[0] = 0;
  for (int iproc = 1; iproc < mpi_size; iproc++)
    disp[iproc] = disp[iproc-1] + nst_ranks[iproc-1]*3;
  assert((disp[mpi_size-1]+nst_ranks[mpi_size-1]*3)==(nst_zn*3));

  int * rcv_count = new int [mpi_size];
  for (int iproc = 0; iproc < mpi_size; iproc++)
    rcv_count[iproc] = nst_ranks[iproc]*3;

  assert(spost_zn == NULL);
  spost_zn = new int [nst_zn][3];
  MPI_Allgatherv((int*)my_spost,my_nst*3,MPI_INT,(int*)spost_zn,rcv_count,disp,MPI_INT,mpi_comm);

  // communicate xp
  disp[0] = 0;
  for (int iproc = 1; iproc < mpi_size; iproc++)
    disp[iproc] = disp[iproc-1] + nsp_ranks[iproc-1]*3;
  assert((disp[mpi_size-1]+nsp_ranks[mpi_size-1]*3)==(nsp_zn*3));

  for (int iproc = 0; iproc < mpi_size; iproc++)
    rcv_count[iproc] = nsp_ranks[iproc]*3;

  assert(xp_zn == NULL);
  xp_zn = new double[nsp_zn][3];
  MPI_Allgatherv((double*)my_xp,my_nsp*3,MPI_DOUBLE,(double*)xp_zn,rcv_count,disp,MPI_DOUBLE,mpi_comm);

  assert(cdf_area_tri == NULL);
  cdf_area_tri = new double [nst_zn+1];
  cdf_area_tri[0] = 0.;
  double avg_normal[3] = {0,0,0};
  for (int ist = 0; ist < nst_zn; ist++) {
    const double * const x0 = xp_zn[spost_zn[ist][0]];
    const double * const x1 = xp_zn[spost_zn[ist][1]];
    const double * const x2 = xp_zn[spost_zn[ist][2]];

    // twice the tri normal
    double normal[3] = TRI_NORMAL_2(x0,x2,x1);
    double area = MAG(normal) / 2.;
    cdf_area_tri[ist+1] = cdf_area_tri[ist] + area;

    FOR_I3 avg_normal[i] += normal[i];
  }

  FOR_I3 avg_normal[i] /= (2.*cdf_area_tri[nst_zn]);

  for (int ist = 0; ist < nst_zn+1; ist++)
    cdf_area_tri[ist] /= cdf_area_tri[nst_zn];

  IF_RANK0 {
    cout << "   > zone info, nst_zn: " << nst_zn << " nsp_zn: " << nsp_zn << " average normal: " << COUT_VEC(avg_normal) << endl;
    //  // dump as tecplot...
    //  FILE * fp = fopen("tris_zone.dat","w");
    //  assert(fp != NULL);
    //  fprintf(fp,"TITLE = \"tris\"\n");
    //  fprintf(fp,"VARIABLES = \"X\"\n");
    //  fprintf(fp,"\"Y\"\n");
    //  fprintf(fp,"\"Z\"\n");
    //  // zone header
    //  fprintf(fp,"ZONE T=\"blah\"\n");
    //  fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nsp_zn,nst_zn);

    //  // order should be good...
    //  for (int isp = 0; isp < nsp_zn; ++isp) {
    //    fprintf(fp,"%lf %lf %lf\n",xp_zn[isp][0],xp_zn[isp][1],xp_zn[isp][2]);
    //  }

    //  for (int ist = 0; ist < nst_zn; ++ist) {
    //    fprintf(fp,"%d %d %d\n",spost_zn[ist][0]+1,spost_zn[ist][1]+1,spost_zn[ist][2]+1);
    //  }

    //  fclose(fp);

    //  for (int ist = 0; ist < nst_zn; ist++)
    //    cout << "ist: " << ist << " isp: " << spost_zn[ist][0] << " " << spost_zn[ist][1] << " " <<spost_zn[ist][2] << " " << COUT_VEC(xp_zn[spost_zn[ist][0]]) << " " << COUT_VEC(xp_zn[spost_zn[ist][1]]) << " " << COUT_VEC(xp_zn[spost_zn[ist][2]]) << endl;
    //  for (int ist = 0; ist < nst_zn+1; ist++)
    //    cout << "ist: " << ist << " cdf_area_tri: " << cdf_area_tri[ist] << endl;
  }

  delete[] my_spost;
  delete[] my_xp;
  delete[] nst_ranks;
  delete[] nsp_ranks;
  delete[] rcv_count;
  delete[] disp;
  return 0;

}

int StaticSolver::addStructuredParticlesInBox(vector<LpDataBase>& lpDataVec, const double xBox[6], const int npx, const int npy, const int npz, const double dp, const string injectorName) {

  int my_np_count = 0;

  double x_min = xBox[0]; // <xmin>
  double x_max = xBox[1]; // <xmax>
  double y_min = xBox[2]; // <ymin>
  double y_max = xBox[3]; // <ymax>
  double z_min = xBox[4]; // <zmin>
  double z_max = xBox[5]; // <zmax>

  int np_total = npx * npy * npz;
  double dx = (x_max - x_min)/double(npx);
  double dy = (y_max - y_min)/double(npy);
  double dz = (z_max - z_min)/double(npz);

  for (int ip = 0; ip < np_total; ip++) {
    // set xp_seed
    double xp_seed[3];
    int iz = ip / (npx*npy);
    int iy = (ip-npx*npy*iz) / npx;
    int ix = ip % npx;
    xp_seed[0] = x_min + (double(ix)+0.5)*dx;
    xp_seed[1] = y_min + (double(iy)+0.5)*dy;
    xp_seed[2] = z_min + (double(iz)+0.5)*dz;
    assert(xp_seed[0]>=x_min);
    assert(xp_seed[0]<=x_max);
    assert(xp_seed[1]>=y_min);
    assert(xp_seed[1]<=y_max);
    assert(xp_seed[2]>=z_min);
    assert(xp_seed[2]<=z_max);

    int icv_closest;
    const int isInside = pointIsInside(xp_seed, icv_closest);
    if ( isInside ) {
      assert( (icv_closest >= 0) && (icv_closest < ncv));
      lpDataVec.resize(lpDataVec.size()+1);
      int iback = lpDataVec.size()-1;
      lpDataVec[iback].icv = icv_closest;
      FOR_I3 lpDataVec[iback].xp[i] = xp_seed[i];
      lpDataVec[iback].dp = dp;
      my_np_count++;
    }

  }

  // we are finished
  int count;
  MPI_Reduce(&my_np_count,&count,1,MPI_INT,MPI_SUM,0,mpi_comm);
  IF_RANK0 cout << " > injecting " << count << " tracers for STRUCTURED injector: " << injectorName << endl;

  return count;
}

int StaticSolver::addEachCellParticlesInBox(vector<LpDataBase>& lpDataVec, const double xBox[6], const int np_in_cell, const double dp, string injectorName) {

  assert(np_in_cell>0);

  int my_np_count = 0;

  double x_min = xBox[0]; // <xmin>
  double x_max = xBox[1]; // <xmax>
  double y_min = xBox[2]; // <ymin>
  double y_max = xBox[3]; // <ymax>
  double z_min = xBox[4]; // <zmin>
  double z_max = xBox[5]; // <zmax>

  int8 * cvora = NULL;
  buildXora(cvora,ncv);
  assert( cvora[mpi_rank+1] - cvora[mpi_rank] == ncv );

  for (int icv = 0; icv < ncv; icv++) {
    if ((x_vv[icv][0] >= x_min) and (x_vv[icv][0] <= x_max)) {
      if ((x_vv[icv][1] >= y_min) and (x_vv[icv][1] <= y_max)) {
        if ((x_vv[icv][2] >= z_min) and (x_vv[icv][2] <= z_max)) {
          int icv_closest;
          const int isInside = pointIsInside(x_vv[icv], icv_closest);
          if ( isInside ) {
            assert(icv == icv_closest);
            for (int ip = 0; ip < np_in_cell; ip++) {
              lpDataVec.resize(lpDataVec.size()+1);
              int iback = lpDataVec.size()-1;
              lpDataVec[iback].icv = icv_closest;
              FOR_I3 lpDataVec[iback].xp[i] = x_vv[icv][i];
              lpDataVec[iback].dp = dp;
              my_np_count++;
            }
          }
        }
      }
    }
  }
  delete[] cvora;

  // we are finished
  int count;
  MPI_Reduce(&my_np_count,&count,1,MPI_INT,MPI_SUM,0,mpi_comm);
  IF_RANK0 cout << " > injecting " << count << " tracers for EACHCELL injector: " << injectorName << endl;

  return count;
}

void StaticSolver::computeBfDistanceFromFlaggedCvs(double * wall_dist, const BfZone * zone_ptr, const double * cv_flag) const {
  // wall distance should already be allocated?
  const int nbf = zone_ptr->nbf;
  FOR_IBF {
    wall_dist[ibf] = -1.0;
  }

  // push flagged cells (flag contains wall-dist) distances to their bfs
  int myCount = 0;  // number of ibfs with wall distance set
  FOR_IBF {
    int icv = zone_ptr->cvobf[ibf];
    if (cv_flag[icv] > 0.0) {
      assert(cv_flag[icv]!=HUGE_VAL);
      wall_dist[ibf] = cv_flag[icv];
      myCount++;
    }
  }

  // All bfs should have either the distance to the wall, or a negative value.
  // We have also counted the positive values. We need to communicate their distance to the wall and coordinates to all processors.
  double * myValues = new double[4*myCount];  // share x,y,z and wall-distance

  //    Now populate
  int check = 0;
  FOR_IBF {
    if (wall_dist[ibf] > 0.0) {
      const int index0 = 4*check;
      FOR_I3 myValues[index0 + i] = zone_ptr->x_bf[ibf][i];
      myValues[index0 + 3] = wall_dist[ibf];
      ++check;
    }
  }
  assert(check == myCount);

  // share counts and offsets
  int totalCount = 0;
  int * numbers = new int[mpi_size];
  int * displ = new int[mpi_size];
  MPI_Allgather(&myCount,1,MPI_INT,numbers,1,MPI_INT,mpi_comm);
  FOR_RANK {
    totalCount += numbers[rank];
    numbers[rank] *= 4;
    if(rank==0) displ[rank]=0;
    else displ[rank] = displ[rank-1] + numbers[rank-1];
  }

  double * allValues = new double[4*totalCount];
  MPI_Allgatherv(myValues,4*myCount,MPI_DOUBLE,allValues,numbers,displ,MPI_DOUBLE,mpi_comm);
  DELETE(numbers);
  DELETE(displ);
  DELETE(myValues);

  // Now we have all the boundary corner faces and their distance to the wall. The next step is to loop over all boundary faces and set the length.
  FOR_IBF {
    if (wall_dist[ibf] < 0.0) {
      wall_dist[ibf] = HUGE_VAL;

      double x_d[3];
      for(int iv=0; iv < totalCount; ++iv) {
        const int index0 = 4*iv;
        FOR_I3 x_d[i] = allValues[index0 + i];  // location of that ibf

        double distance = DIST(x_d,zone_ptr->x_bf[ibf]);
        distance += allValues[index0 + 3];  // wall_distance + distance to that ibf
        if (distance < wall_dist[ibf]) wall_dist[ibf] = distance;
      }
    }
  }
  DELETE(allValues);
}
