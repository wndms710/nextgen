#ifndef CONDITIONALPROBE_HPP
#define CONDITIONALPROBE_HPP

#include "SimpleGeom.hpp"

class ConditionalProbe {
public:

  string name;
  int write_interval;
  int sample_interval;
  bool reset_after_write;

  SimpleGeom * geom;
  vector<pair<string,CtiRegister::CtiData*> > var_vec;

  pair<string,CtiRegister::CtiData*> wvar;

  pair<string,CtiRegister::CtiData*> cvar;
  int nbin;
  int ibin_min,ibin_max;
  double range[2];
  bool b_range;

  pair<string,CtiRegister::CtiData*> cvar_2;
  int nbin_2;
  int ibin_min_2,ibin_max_2;
  double range_2[2];
  bool b_range_2;

  double * wgt_buf;
  double * avg_buf;
  double * rms_buf;
  double * min_buf;
  double * max_buf;

  list<pair<string,stringstream*> > cpFileList;

  ConditionalProbe() {

    name        = "";
    write_interval  = -1; // -1 used to indicate was never set
    sample_interval = -1;
    reset_after_write = false; // default is to accumulate stats indefinately

    geom = NULL;

    nbin      = 100;    // default bin size
    ibin_min  = 0;      // default to full range
    ibin_max  = nbin+1; // +2 to encompass out of range data
    b_range   = false;  // by default range will be set by first doProbe...

    // for 2d binning...
    nbin_2       = 100;      // default bin size
    ibin_min_2   = 0;        // default to full range
    ibin_max_2   = nbin_2+1; // +2 to encompass out of range data
    b_range_2    = false;    // by default range will be set by first doProbe...
    cvar_2.first = "";

    wvar.first  = "";   // unit weighting when we don't get wgt var...
    wvar.second = NULL;

    wgt_buf = NULL;
    avg_buf = NULL;
    rms_buf = NULL;
    min_buf = NULL;
    max_buf = NULL;

  }

  ~ConditionalProbe() {

    if (geom) delete geom;

    DELETE(wgt_buf);
    DELETE(avg_buf);
    DELETE(rms_buf);
    DELETE(min_buf);
    DELETE(max_buf);

    for (list<pair<string,stringstream*> >::iterator iter = cpFileList.begin(); iter != cpFileList.end(); ++iter)
      delete iter->second;

  }

  void doProbe(const int step,const double time) {

    static int rank0 = 0; // used to round-robin the writing...

    // conditional variable...
    CtiRegister::CtiData * cvar_data = cvar.second;
    if (cvar_data == NULL) cvar_data = CtiRegister::getUnregisteredCtiData(cvar.first);
    assert(cvar_data);

    // wgt...
    double *wvar_ptr = NULL;
    CtiRegister::CtiData * wvar_data = wvar.second;
    if (wvar_data != NULL) {
      // registered data
      wvar_ptr = wvar_data->getDNptr();
    }
    else {
      if (wvar.first != "") {
        // unregistered data
        wvar_data = CtiRegister::getUnregisteredCtiData(wvar.first);
        assert(wvar_data);
        wvar_ptr = wvar_data->getDNptr();
      }
      else {
        // constant wgt (default)
        wvar_ptr = new double[cvar_data->size()];
        for (int ii = 0; ii < cvar_data->size(); ++ii) wvar_ptr[ii] = 1.0;
      }
    }
    assert(wvar_ptr);

    // if geom is present, then the user is restricting the analysis using
    // the pointIsInside capability of the geom...
    bool * is_inside = NULL;
    if (geom) {
      // for now, must be cv data, because we don't know how to access coordinate
      // data otherwise. This could be made part of the CONDITIONAL_PROBE...
      CtiRegister::CtiData * x_cv_data;
      if (cvar_data->getUnindexedTopology() == BF_DATA) {
        const string zone_name = CtiRegister::getFirstSpecifier(cvar.first);
        x_cv_data = CtiRegister::getRegisteredCtiData(zone_name+":x_bf");
      }
      else if (cvar_data->getName() == "cht") {
        //cvar_data->getX(x_cv_data);  // doesn't support ctiData * 
        x_cv_data = CtiRegister::getRegisteredCtiData("cht:x_no");
      }
      else x_cv_data = CtiRegister::getRegisteredCtiData("x_cv");

      assert(x_cv_data);
      assert(x_cv_data->getType() == DN3_DATA);
//      if (x_cv_data->size() != cvar_data->size()) {
//        if (mpi_rank == 0) {
//          cout << "## x_cv_data->size(): " << x_cv_data->size() << endl;
//           cout << "## cvar_data->size(): " << cvar_data->size() << endl; 
//        }
//      }
      assert(x_cv_data->size() == cvar_data->size());
      is_inside = new bool[cvar_data->size()];
      double x_cv[3];
      for (int i = 0; i < cvar_data->size(); ++i) {
      	FOR_J3 x_cv[j] = x_cv_data->dn3(i,j);
      	is_inside[i] = geom->pointIsInside(x_cv);
      }
    }

    double my_range[2] = { 1.0E+20, 1.0E+20 };
    int my_count = 0;
    for (int i = 0; i < cvar_data->size(); ++i) {
      if ((is_inside == NULL)||(is_inside[i])) {
        ++my_count;
        const double c = cvar_data->dn(i);
        if (c < my_range[0]) my_range[0] = c;
        if (-c < my_range[1]) my_range[1] = -c;
      }
    }
    if (!b_range) {
      MPI_Allreduce(my_range,range,2,MPI_DOUBLE,MPI_MIN,mpi_comm);
      range[1] = -range[1]; // flip the max
      b_range = true;
    }
    my_range[1] = -my_range[1];

    if (range[0] == range[1]) {
      CWARN("Conditional variable " << cvar.first << " is constant. Skipping...");
      return;
    }

    if ( my_count == 0) {

      // due to the presence of the geom restriction, it is possible that this
      // rank has no participating data at all..

      ibin_min  = 0;
      ibin_max  = -1;

    } else {

      if (my_range[0] < range[0]) {
        // below lower bound...
        ibin_min = 0;
      }
      else if (my_range[0] > range[1]) {
        // above upper bound...
        ibin_min = nbin+1;
      }
      else {
        // should be [1:nbin] inclusive. need min operation to include dn == range[1] in last valid bin...
        ibin_min = min(nbin,1+(int)floor(double(nbin)*(my_range[0]-range[0])/(range[1]-range[0])));
      }
      if (my_range[1] < range[0]) {
        // below lower bound...
        ibin_max = 0;
      }
      else if (my_range[1] > range[1]) {
        // above upper bound...
        ibin_max = nbin+1;
      }
      else {
        // should be [1:nbin] inclusive. need min operation to include dn == range[1] in last valid bin...
        ibin_max = min(nbin,1+(int)floor(double(nbin)*(my_range[1]-range[0])/(range[1]-range[0])));
      }
      //cout << mpi_rank << " " << ibin_min << ":" << ibin_max << endl; cout.flush();
    }

    // for each variable we are conditionally probing, we add the avg,rms,min,max

    const int nvar = var_vec.size();
    assert(nvar > 0);

    if (cvar_2.first == "") {

      // 1d conditional probe...

      if (wgt_buf == NULL) {

        assert(avg_buf == NULL);
        assert(rms_buf == NULL);
        assert(min_buf == NULL);
        assert(max_buf == NULL);

        // the size of these buffers is nbin+2 per scalar nvar...

        //wgt_buf = new double[(nbin+2)*3]; // stores wgt, wgt*c and wgt*c*c
        //avg_buf = new double[(nbin+2)*nvar];
        //rms_buf = new double[(nbin+2)*nvar];
        //min_buf = new double[(nbin+2)*nvar];
        //max_buf = new double[(nbin+2)*nvar];
        wgt_buf = new double[(ibin_max-ibin_min+1)*3]; // stores wgt, wgt*c and wgt*c*c
        avg_buf = new double[(ibin_max-ibin_min+1)*nvar];
        rms_buf = new double[(ibin_max-ibin_min+1)*nvar];
        min_buf = new double[(ibin_max-ibin_min+1)*nvar];
        max_buf = new double[(ibin_max-ibin_min+1)*nvar];

        resetBuffers();

      }

      // set the bin index for each piece of data...

      int * ibin_buf = new int[cvar_data->size()];
      for (int i = 0; i < cvar_data->size(); ++i) {
        if ((is_inside == NULL)||(is_inside[i])) {
          const double c = cvar_data->dn(i);
          if (c < range[0]) {
            // below lower bound...
            ibin_buf[i] = 0;
          }
          else if (c > range[1]) {
            // above upper bound...
            ibin_buf[i] = nbin+1;
          }
          else {
            // should be [1:nbin] inclusive. need min operation to include dn == range[1] in last valid bin...
            ibin_buf[i] = min(nbin,1+(int)floor(double(nbin)*(c-range[0])/(range[1]-range[0])));
          }
          // and accumulate in wgt_buf...
          //wgt_buf[ibin_buf[i]]            += wvar_ptr[i];
          //wgt_buf[(nbin+2)+ibin_buf[i]]   += wvar_ptr[i]*c;
          //wgt_buf[(nbin+2)*2+ibin_buf[i]] += wvar_ptr[i]*c*c;
          assert((ibin_buf[i] >= ibin_min)&&(ibin_buf[i] <= ibin_max));
          wgt_buf[ibin_buf[i]-ibin_min]                         += wvar_ptr[i];
          wgt_buf[(ibin_max-ibin_min+1)+ibin_buf[i]-ibin_min]   += wvar_ptr[i]*c;
          wgt_buf[(ibin_max-ibin_min+1)*2+ibin_buf[i]-ibin_min] += wvar_ptr[i]*c*c;
        }
      }

      // now cycle through and add our stats into the bufs...

      for (int ivar = 0; ivar < nvar; ++ivar) {

        CtiRegister::CtiData * var_data = var_vec[ivar].second;
        if (var_data == NULL) var_data = CtiRegister::getUnregisteredCtiData(var_vec[ivar].first);
        assert(var_data);
        assert(var_data->size() == cvar_data->size());
        for (int i = 0; i < var_data->size(); ++i) {
          if ((is_inside == NULL)||(is_inside[i])) {
            const double v = var_data->dn(i);
            const int ibin = ibin_buf[i];
            //avg_buf[ivar*(nbin+2)+ibin] += wvar_ptr[i]*v;
            //rms_buf[ivar*(nbin+2)+ibin] += wvar_ptr[i]*v*v;
            //min_buf[ivar*(nbin+2)+ibin] = min(min_buf[ivar*(nbin+2)+ibin],v);
            //max_buf[ivar*(nbin+2)+ibin] = max(max_buf[ivar*(nbin+2)+ibin],v);
            assert((ibin_buf[i] >= ibin_min)&&(ibin_buf[i] <= ibin_max));
            avg_buf[ivar*(ibin_max-ibin_min+1)+ibin-ibin_min] += wvar_ptr[i]*v;
            rms_buf[ivar*(ibin_max-ibin_min+1)+ibin-ibin_min] += wvar_ptr[i]*v*v;
            min_buf[ivar*(ibin_max-ibin_min+1)+ibin-ibin_min] = min(min_buf[ivar*(ibin_max-ibin_min+1)+ibin-ibin_min],v);
            max_buf[ivar*(ibin_max-ibin_min+1)+ibin-ibin_min] = max(max_buf[ivar*(ibin_max-ibin_min+1)+ibin-ibin_min],v);
          }
        }

      }
      delete[] ibin_buf;
      if (is_inside) delete[] is_inside;

      if ((step == -1)||(step%write_interval == 0)) {

        // reduce and write into probe_rank0...

        // for now, we reduce EVERYTHING. probably a better way to do this some day...

        double * wgt_buf_global = NULL;
        double * avg_buf_global = NULL;
        double * rms_buf_global = NULL;
        double * min_buf_global = NULL;
        double * max_buf_global = NULL;

        if (mpi_rank == rank0) {

          wgt_buf_global = new double[(nbin+2)*3];
          avg_buf_global = new double[(nbin+2)*nvar];
          rms_buf_global = new double[(nbin+2)*nvar];
          min_buf_global = new double[(nbin+2)*nvar];
          max_buf_global = new double[(nbin+2)*nvar];

        }

        // and do reductions...

        //MPI_Reduce(wgt_buf,wgt_buf_global,(nbin+2)*3,MPI_DOUBLE,MPI_SUM,rank0,mpi_comm);
        //MPI_Reduce(avg_buf,avg_buf_global,(nbin+2)*nvar,MPI_DOUBLE,MPI_SUM,rank0,mpi_comm);
        //MPI_Reduce(rms_buf,rms_buf_global,(nbin+2)*nvar,MPI_DOUBLE,MPI_SUM,rank0,mpi_comm);
        //MPI_Reduce(min_buf,min_buf_global,(nbin+2)*nvar,MPI_DOUBLE,MPI_MIN,rank0,mpi_comm);
        //MPI_Reduce(max_buf,max_buf_global,(nbin+2)*nvar,MPI_DOUBLE,MPI_MAX,rank0,mpi_comm);

        // init global bufs...
        if (mpi_rank == rank0) {
          for (int i = 0, limit = (nbin+2)*3; i < limit; ++i) {
            wgt_buf_global[i] = 0.0;
          }
          for (int i = 0, limit = (nbin+2)*var_vec.size(); i < limit; ++i) {
            avg_buf_global[i] = 0.0;
            rms_buf_global[i] = 0.0;
            min_buf_global[i] = 1.0E+20;
            max_buf_global[i] = -1.0E+20;
          }
        }

        // gather local ranges...
        int my_ibin_range[2] = {ibin_min,ibin_max};
        int (*ibin_range)[2] = NULL;
        if (mpi_rank == rank0)
          ibin_range = new int[mpi_size][2];
        MPI_Gather(my_ibin_range,2,MPI_INT,(int*)ibin_range,2,MPI_INT,rank0,mpi_comm);
        int * count = NULL;
        int * disp = NULL;
        double *recv_buf_double = NULL;
        if (mpi_rank == rank0) {
          count = new int[mpi_size];
          FOR_RANK count[rank] = ibin_range[rank][1]-ibin_range[rank][0]+1;
          disp  = new int[mpi_size];
          disp[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank)
            disp[rank] = disp[rank-1] + count[rank-1];
          int recv_count_sum = disp[mpi_size-1] + count[mpi_size-1]; // int8?
          recv_buf_double = new double[max(3,nvar)*recv_count_sum];
        }

        // now gather the double data and populate
        if (mpi_rank == rank0) {
          FOR_RANK {
            count[rank] *= 3;
            disp[rank]  *= 3;
          //cout << rank << " " << ibin_range[rank][0] << ":" << ibin_range[rank][1] << " " <<  count[rank] << " " << disp[rank] << endl;
          }
        }
        MPI_Gatherv(wgt_buf,(ibin_max-ibin_min+1)*3,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            FOR_I3 {
              for (int ii = 0; ii < count[rank]/3; ++ii) {
                wgt_buf_global[(nbin+2)*i+ii+ibin_range[rank][0]] += recv_buf_double[disp[rank]+ii+(count[rank]/3)*i];
              }
            }
          }
        }
        if (mpi_rank == rank0) {
          FOR_RANK {
            count[rank] = (count[rank]/3)*nvar;
            disp[rank]  = (disp[rank]/3)*nvar;
          }
        }
        MPI_Gatherv(avg_buf,(ibin_max-ibin_min+1)*nvar,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            for (int ivar = 0; ivar < nvar; ++ivar) {
              for (int ii = 0; ii < count[rank]/nvar; ++ii) {
                avg_buf_global[(nbin+2)*ivar+ii+ibin_range[rank][0]] += recv_buf_double[disp[rank]+ii+(count[rank]/nvar)*ivar];
              }
            }
          }
        }
        MPI_Gatherv(rms_buf,(ibin_max-ibin_min+1)*nvar,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            for (int ivar = 0; ivar < nvar; ++ivar) {
              for (int ii = 0; ii < count[rank]/nvar; ++ii) {
                rms_buf_global[(nbin+2)*ivar+ii+ibin_range[rank][0]] += recv_buf_double[disp[rank]+ii+(count[rank]/nvar)*ivar];
              }
            }
          }
        }
        MPI_Gatherv(min_buf,(ibin_max-ibin_min+1)*nvar,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            for (int ivar = 0; ivar < nvar; ++ivar) {
              for (int ii = 0; ii < count[rank]/nvar; ++ii) {
                min_buf_global[(nbin+2)*ivar+ii+ibin_range[rank][0]] = min(min_buf_global[(nbin+2)*ivar+ii+ibin_range[rank][0]],recv_buf_double[disp[rank]+ii+(count[rank]/nvar)*ivar]);
              }
            }
          }
        }
        MPI_Gatherv(max_buf,(ibin_max-ibin_min+1)*nvar,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            for (int ivar = 0; ivar < nvar; ++ivar) {
              for (int ii = 0; ii < count[rank]/nvar; ++ii) {
                max_buf_global[(nbin+2)*ivar+ii+ibin_range[rank][0]] = max(max_buf_global[(nbin+2)*ivar+ii+ibin_range[rank][0]],recv_buf_double[disp[rank]+ii+(count[rank]/nvar)*ivar]);
              }
            }
          }
        }
        if (mpi_rank == rank0) {
          delete[] ibin_range;
          delete[] count;
          delete[] disp;
          delete[] recv_buf_double;
        }

        if (mpi_rank == rank0) {

          char filename[128];
          if (step == -1) {
            sprintf(filename,"%s.cp",name.c_str());
          }
          else {
            sprintf(filename,"%s.%08d.cp",name.c_str(),step);
          }

          cpFileList.push_back(pair<string,stringstream*>(string(filename),new stringstream()));
          stringstream* ss = cpFileList.back().second;

          // add the header...

          *ss << setprecision(8);
          *ss << "# step " << step << " time " << time << endl;
          if (wvar.first == "")
            *ss << "# 1:bin 2:count 3:" << cvar.first << "_avg 4:" << cvar.first << "_rms";
          else
            *ss << "# 1:bin 2:total 3:" << cvar.first << "_avg 4:" << cvar.first << "_rms";
          for (int i = 0; i < int(var_vec.size()); ++i)
            *ss << " " << 5+i*4 << ":" << var_vec[i].first << "_avg " <<
              6+i*4 << ":" << var_vec[i].first << "_rms " <<
              7+i*4 << ":" << var_vec[i].first << "_min " <<
              8+i*4 << ":" << var_vec[i].first << "_max";
          *ss << endl;

          for (int ibin = 0; ibin < nbin+2; ++ibin) {
            if (wgt_buf_global[ibin] > 0.0) {
              *ss << ibin << " " << wgt_buf_global[ibin];
              // this bin contains valid data. Report the bin center using the
              // c average and rms...
              const double cavg = wgt_buf_global[(nbin+2)+ibin]/wgt_buf_global[ibin];
              const double crms = sqrt(max(0.0,wgt_buf_global[(nbin+2)*2+ibin]/wgt_buf_global[ibin] - cavg*cavg));
              *ss << " " << cavg << " " << crms;
              for (int ivar = 0; ivar < nvar; ++ivar) {
                const double avg = avg_buf_global[(nbin+2)*ivar+ibin]/wgt_buf_global[ibin];
                const double rms = sqrt(max(0.0,rms_buf_global[(nbin+2)*ivar+ibin]/wgt_buf_global[ibin] - avg*avg));
                *ss << " " << avg << " " << rms << " " << min_buf_global[(nbin+2)*ivar+ibin] << " " << max_buf_global[(nbin+2)*ivar+ibin];
              }
              *ss << endl;
            }
          }

          // cleanup...

          delete[] wgt_buf_global;
          delete[] avg_buf_global;
          delete[] rms_buf_global;
          delete[] min_buf_global;
          delete[] max_buf_global;

        }

        // advance rank0 for the next one...

        ++rank0;
        if (rank0 == mpi_size)
          rank0 = 0;

        // and rest the averages if reset_after_write...

        if (reset_after_write)
          resetBuffers();

      }

    }
    else {

      // 2d conditional probe...

      CtiRegister::CtiData * cvar_2_data = cvar_2.second;
      if (cvar_2_data == NULL) cvar_2_data = CtiRegister::getUnregisteredCtiData(cvar_2.first);
      assert(cvar_2_data);

      double my_range_2[2] = { 1.0E+20, 1.0E+20 };
      int my_count_2 = 0;
      for (int i = 0; i < cvar_2_data->size(); ++i) {
        if ((is_inside == NULL)||(is_inside[i])) {
          ++my_count_2;
          const double c = cvar_2_data->dn(i);
          if (c < my_range_2[0]) my_range_2[0] = c;
          if (-c < my_range_2[1]) my_range_2[1] = -c;
        }
      }
      if (!b_range_2) {
        MPI_Allreduce(my_range_2,range_2,2,MPI_DOUBLE,MPI_MIN,mpi_comm);
        range_2[1] = -range_2[1]; // flip the max
        b_range_2 = true;
      }
      my_range_2[1] = -my_range_2[1];

      if (range_2[0] == range_2[1]) {
        CWARN("Conditional variable " << cvar_2.first << " is constant. Skipping...");
        return;
      }

      if ( my_count_2 == 0) {

        // due to the presence of the geom restriction, it is possible that this
        // rank has no participating data at all..

        ibin_min_2  = 0;
        ibin_max_2  = -1;

      } else {
        if (my_range_2[0] < range_2[0]) {
          // below lower bound...
          ibin_min_2 = 0;
        }
        else if (my_range_2[0] > range_2[1]) {
          // above upper bound...
          ibin_min_2 = nbin_2+1;
        }
        else {
          // should be [1:nbin_2] inclusive. need min operation to include dn == range_2[1] in last valid bin...
          ibin_min_2 = min(nbin_2,1+(int)floor(double(nbin_2)*(my_range_2[0]-range_2[0])/(range_2[1]-range_2[0])));
        }
        if (my_range_2[1] < range_2[0]) {
          // below lower bound...
          ibin_max_2 = 0;
        }
        else if (my_range_2[1] > range_2[1]) {
          // above upper bound...
          ibin_max_2 = nbin_2+1;
        }
        else {
          // should be [1:nbin_2] inclusive. need min operation to include dn == range_2[1] in last valid bin...
          ibin_max_2 = min(nbin_2,1+(int)floor(double(nbin_2)*(my_range_2[1]-range_2[0])/(range_2[1]-range_2[0])));
        }
      }

      if (wgt_buf == NULL) {

        assert(avg_buf == NULL);
        assert(rms_buf == NULL);
        assert(min_buf == NULL);
        assert(max_buf == NULL);

        // the size of these buffers is (nbin+2)*(nbin_2+2) per scalar nvar...

        //wgt_buf = new double[(nbin+2)*(nbin_2+2)*5]; // stores wgt, wgt*c, wgt*c*c, wgt*c_2 and wgt*c_2*c_2
        //avg_buf = new double[(nbin+2)*(nbin_2+2)*nvar];
        //rms_buf = new double[(nbin+2)*(nbin_2+2)*nvar];
        //min_buf = new double[(nbin+2)*(nbin_2+2)*nvar];
        //max_buf = new double[(nbin+2)*(nbin_2+2)*nvar];
        wgt_buf = new double[(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*5]; // stores wgt,wgt*c,wgt*c*c,wgt*c_2,wgt*c_2*c_2
        avg_buf = new double[(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*nvar];
        rms_buf = new double[(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*nvar];
        min_buf = new double[(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*nvar];
        max_buf = new double[(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*nvar];

        resetBuffers();

      }

      // set the bin index for each piece of data...

      assert(cvar_2_data->size() == cvar_data->size());
      int * ibin_buf = new int[cvar_data->size()];
      int * jbin_buf = new int[cvar_2_data->size()];
      for (int ii = 0; ii < cvar_data->size(); ++ii) {
        if ((is_inside == NULL)||(is_inside[ii])) {
          const double c = cvar_data->dn(ii);
          if (c < range[0]) {
            // below lower bound...
            ibin_buf[ii] = 0;
          }
          else if (c > range[1]) {
            // above upper bound...
            ibin_buf[ii] = nbin+1;
          }
          else {
            // should be [1:nbin] inclusive. need min operation to include dn == range[1] in last valid bin...
            ibin_buf[ii] = min(nbin,1+(int)floor(double(nbin)*(c-range[0])/(range[1]-range[0])));
          }
          const double c_2 = cvar_2_data->dn(ii);
          if (c_2 < range_2[0]) {
            // below lower bound...
            jbin_buf[ii] = 0;
          }
          else if (c_2 > range_2[1]) {
            // above upper bound...
            jbin_buf[ii] = nbin_2+1;
          }
          else {
            // should be [1:nbin_2] inclusive. need min operation to include dn == range_2[1] in last valid bin...
            jbin_buf[ii] = min(nbin_2,1+(int)floor(double(nbin_2)*(c_2-range_2[0])/(range_2[1]-range_2[0])));
          }
          // and accumulate in wgt_buf...
          //const int kbin = ibin_buf[ii]*(nbin_2+2)+jbin_buf[ii];
          //wgt_buf[                      kbin] += wvar_ptr[ii];
          //wgt_buf[(nbin+2)*(nbin_2+2)  +kbin] += wvar_ptr[ii]*c;
          //wgt_buf[(nbin+2)*(nbin_2+2)*2+kbin] += wvar_ptr[ii]*c*c;
          //wgt_buf[(nbin+2)*(nbin_2+2)*3+kbin] += wvar_ptr[ii]*c_2;
          //wgt_buf[(nbin+2)*(nbin_2+2)*4+kbin] += wvar_ptr[ii]*c_2*c_2;
          assert((ibin_buf[ii] >= ibin_min)&&(ibin_buf[ii] <= ibin_max));
          assert((jbin_buf[ii] >= ibin_min_2)&&(jbin_buf[ii] <= ibin_max_2));
          const int index = (ibin_buf[ii]-ibin_min)*(ibin_max_2-ibin_min_2+1)+jbin_buf[ii]-ibin_min_2;
          wgt_buf[                                                  index] += wvar_ptr[ii];
          wgt_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)  +index] += wvar_ptr[ii]*c;
          wgt_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*2+index] += wvar_ptr[ii]*c*c;
          wgt_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*3+index] += wvar_ptr[ii]*c_2;
          wgt_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*4+index] += wvar_ptr[ii]*c_2*c_2;
        }
      }

      // now cycle through and add our stats into the bufs...

      for (int ivar = 0; ivar < nvar; ++ivar) {

        CtiRegister::CtiData * var_data = var_vec[ivar].second;
        if (var_data == NULL) var_data = CtiRegister::getUnregisteredCtiData(var_vec[ivar].first);
        assert(var_data);
        assert(var_data->size() == cvar_data->size());
        for (int ii = 0; ii < var_data->size(); ++ii) {
          if ((is_inside == NULL)||(is_inside[ii])) {
            const double v = var_data->dn(ii);
            //const int index = ivar*(nbin+2)*(nbin_2+2)+ibin_buf[ii]*(nbin_2+2)+jbin_buf[ii];
            //avg_buf[index] += wvar_ptr[ii]*v;
            //rms_buf[index] += wvar_ptr[ii]*v*v;
            //min_buf[index] = min(min_buf[index],v);
            //max_buf[index] = max(max_buf[index],v);
            assert((ibin_buf[ii] >= ibin_min)&&(ibin_buf[ii] <= ibin_max));
            assert((jbin_buf[ii] >= ibin_min_2)&&(jbin_buf[ii] <= ibin_max_2));
            const int index = (ibin_buf[ii]-ibin_min)*(ibin_max_2-ibin_min_2+1)+jbin_buf[ii]-ibin_min_2;
            avg_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*ivar+index] += wvar_ptr[ii]*v;
            rms_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*ivar+index] += wvar_ptr[ii]*v*v;
            min_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*ivar+index] = min(min_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*ivar+index],v);
            max_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*ivar+index] = max(max_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*ivar+index],v);
          }
        }

      }
      delete[] ibin_buf;
      delete[] jbin_buf;
      if (is_inside) delete[] is_inside;

      if ((step == -1)||(step%write_interval == 0)) {

        // reduce and write into probe_rank0...

        // for now, we reduce EVERYTHING. probably a better way to do this some day...

        double * wgt_buf_global = NULL;
        double * avg_buf_global = NULL;
        double * rms_buf_global = NULL;
        double * min_buf_global = NULL;
        double * max_buf_global = NULL;

        if (mpi_rank == rank0) {

          wgt_buf_global = new double[(nbin+2)*(nbin_2+2)*5];
          avg_buf_global = new double[(nbin+2)*(nbin_2+2)*nvar];
          rms_buf_global = new double[(nbin+2)*(nbin_2+2)*nvar];
          min_buf_global = new double[(nbin+2)*(nbin_2+2)*nvar];
          max_buf_global = new double[(nbin+2)*(nbin_2+2)*nvar];

        }

        // and do reductions...

        //MPI_Reduce(wgt_buf,wgt_buf_global,(nbin+2)*(nbin_2+2)*5,MPI_DOUBLE,MPI_SUM,rank0,mpi_comm);
        //MPI_Reduce(avg_buf,avg_buf_global,(nbin+2)*(nbin_2+2)*nvar,MPI_DOUBLE,MPI_SUM,rank0,mpi_comm);
        //MPI_Reduce(rms_buf,rms_buf_global,(nbin+2)*(nbin_2+2)*nvar,MPI_DOUBLE,MPI_SUM,rank0,mpi_comm);
        //MPI_Reduce(min_buf,min_buf_global,(nbin+2)*(nbin_2+2)*nvar,MPI_DOUBLE,MPI_MIN,rank0,mpi_comm);
        //MPI_Reduce(max_buf,max_buf_global,(nbin+2)*(nbin_2+2)*nvar,MPI_DOUBLE,MPI_MAX,rank0,mpi_comm);

        // init global bufs...
        if (mpi_rank == rank0) {
          for (int i = 0, limit = (nbin+2)*(nbin_2+2)*5; i < limit; ++i) {
            wgt_buf_global[i] = 0.0;
          }
          for (int i = 0, limit = (nbin+2)*(nbin_2+2)*var_vec.size(); i < limit; ++i) {
            avg_buf_global[i] = 0.0;
            rms_buf_global[i] = 0.0;
            min_buf_global[i] = 1.0E+20;
            max_buf_global[i] = -1.0E+20;
          }
        }

        // gather local ranges...
        int my_ibin_range[4] = {ibin_min,ibin_max,ibin_min_2,ibin_max_2};
        int (*ibin_range)[4] = NULL;
        if (mpi_rank == rank0)
          ibin_range = new int[mpi_size][4];
        MPI_Gather(my_ibin_range,4,MPI_INT,(int*)ibin_range,4,MPI_INT,rank0,mpi_comm);
        int * count = NULL;
        int * disp = NULL;
        double *recv_buf_double = NULL;
        if (mpi_rank == rank0) {
          count = new int[mpi_size];
          FOR_RANK count[rank] = (ibin_range[rank][1]-ibin_range[rank][0]+1)*(ibin_range[rank][3]-ibin_range[rank][2]+1);
          disp  = new int[mpi_size];
          disp[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank)
            disp[rank] = disp[rank-1] + count[rank-1];
          int recv_count_sum = disp[mpi_size-1] + count[mpi_size-1]; // int8?
          recv_buf_double = new double[max(5,nvar)*recv_count_sum];
        }

        // now gather the double data and populate
        if (mpi_rank == rank0) {
          FOR_RANK {
            count[rank] *= 5;
            disp[rank]  *= 5;
            //cout << rank << " " << count[rank] << " " << disp[rank] << " " << ibin_range[rank][0] << " " << ibin_range[rank][1] << " " << ibin_range[rank][2] << " " << ibin_range[rank][3] << endl;
          }
        }
        MPI_Gatherv(wgt_buf,(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*5,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            FOR_I5 {
              for (int ii = 0; ii < count[rank]/5; ++ii) {
                const int ix = ibin_range[rank][0]+ii/(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int iy = ibin_range[rank][2]+ii%(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int index = (nbin_2+2)*ix+iy;
                wgt_buf_global[(nbin_2+2)*(nbin+2)*i+index] += recv_buf_double[disp[rank]+ii+(count[rank]/5)*i];
              }
            }
          }
        }
        if (mpi_rank == rank0) {
          FOR_RANK {
            count[rank] = (count[rank]/5)*nvar;
            disp[rank]  = (disp[rank]/5)*nvar;
          }
        }
        MPI_Gatherv(avg_buf,(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*nvar,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            for (int ivar = 0; ivar < nvar; ++ivar) {
              for (int ii = 0; ii < count[rank]/nvar; ++ii) {
                const int ix = ibin_range[rank][0]+ii/(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int iy = ibin_range[rank][2]+ii%(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int index = (nbin_2+2)*ix+iy;
                avg_buf_global[(nbin_2+2)*(nbin+2)*ivar+index] += recv_buf_double[disp[rank]+ii+(count[rank]/nvar)*ivar];
              }
            }
          }
        }
        MPI_Gatherv(rms_buf,(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*nvar,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            for (int ivar = 0; ivar < nvar; ++ivar) {
              for (int ii = 0; ii < count[rank]/nvar; ++ii) {
                const int ix = ibin_range[rank][0]+ii/(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int iy = ibin_range[rank][2]+ii%(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int index = (nbin_2+2)*ix+iy;
                rms_buf_global[(nbin_2+2)*(nbin+2)*ivar+index] += recv_buf_double[disp[rank]+ii+(count[rank]/nvar)*ivar];
              }
            }
          }
        }
        MPI_Gatherv(min_buf,(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*nvar,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            for (int ivar = 0; ivar < nvar; ++ivar) {
              for (int ii = 0; ii < count[rank]/nvar; ++ii) {
                const int ix = ibin_range[rank][0]+ii/(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int iy = ibin_range[rank][2]+ii%(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int index = (nbin_2+2)*ix+iy;
                min_buf_global[(nbin_2+2)*(nbin+2)*ivar+index] = min(min_buf_global[(nbin_2+2)*(nbin+2)*ivar+index],recv_buf_double[disp[rank]+ii+(count[rank]/nvar)*ivar]);
              }
            }
          }
        }
        MPI_Gatherv(max_buf,(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*nvar,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            for (int ivar = 0; ivar < nvar; ++ivar) {
              for (int ii = 0; ii < count[rank]/nvar; ++ii) {
                const int ix = ibin_range[rank][0]+ii/(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int iy = ibin_range[rank][2]+ii%(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int index = (nbin_2+2)*ix+iy;
                max_buf_global[(nbin_2+2)*(nbin+2)*ivar+index] = max(max_buf_global[(nbin_2+2)*(nbin+2)*ivar+index],recv_buf_double[disp[rank]+ii+(count[rank]/nvar)*ivar]);
              }
            }
          }
        }
        if (mpi_rank == rank0) {
          delete[] ibin_range;
          delete[] count;
          delete[] disp;
          delete[] recv_buf_double;
        }

        if (mpi_rank == rank0) {

          char filename[128];
          if (step == -1) {
            sprintf(filename,"%s.cp",name.c_str());
          }
          else {
            sprintf(filename,"%s.%08d.cp",name.c_str(),step);
          }

          cpFileList.push_back(pair<string,stringstream*>(string(filename),new stringstream()));
          stringstream* ss = cpFileList.back().second;

          // add the header...

          *ss << setprecision(8);
          *ss << "# step " << step << " time " << time << endl;
          if (wvar.first == "")
            *ss << "# 1:ibin 2:jbin 3:count 4:";
          else
            *ss << "# 1:ibin 2:jbin 3:total 4:";
          *ss << cvar.first << "_avg 5:" << cvar.first << "_rms 6:"
              << cvar_2.first << "_avg 7:" << cvar_2.first << "_rms";
          for (int i = 0; i < int(var_vec.size()); ++i)
            *ss << " " << 8+i*4 << ":" << var_vec[i].first << "_avg " <<
              9+i*4 << ":" << var_vec[i].first << "_rms " <<
              10+i*4 << ":" << var_vec[i].first << "_min " <<
              11+i*4 << ":" << var_vec[i].first << "_max";
          *ss << endl;

          for (int ibin = 0; ibin < nbin+2; ++ibin) {
            for (int jbin = 0; jbin < nbin_2+2; ++jbin) {
              const int kbin = ibin*(nbin_2+2)+jbin;
              if (wgt_buf_global[kbin] > 0.0) {
                *ss << ibin << " " << jbin << " " << wgt_buf_global[kbin];
                // this bin contains valid data. Report the bin center using the
                // c average and rms...
                const double c_avg = wgt_buf_global[(nbin+2)*(nbin_2+2)+kbin]/wgt_buf_global[kbin];
                const double c_rms = sqrt(max(0.0,wgt_buf_global[(nbin+2)*(nbin_2+2)*2+kbin]/wgt_buf_global[kbin] - c_avg*c_avg));
                *ss << " " << c_avg << " " << c_rms;
                const double c_2_avg = wgt_buf_global[(nbin+2)*(nbin_2+2)*3+kbin]/wgt_buf_global[kbin];
                const double c_2_rms = sqrt(max(0.0,wgt_buf_global[(nbin+2)*(nbin_2+2)*4+kbin]/wgt_buf_global[kbin] - c_2_avg*c_2_avg));
                *ss << " " << c_2_avg << " " << c_2_rms;
                for (int ivar = 0; ivar < nvar; ++ivar) {
                  const int index = ivar*(nbin+2)*(nbin_2+2)+kbin;
                  const double avg = avg_buf_global[index]/wgt_buf_global[kbin];
                  const double rms = sqrt(max(0.0,rms_buf_global[index]/wgt_buf_global[kbin] - avg*avg));
                  *ss << " " << avg << " " << rms << " " << min_buf_global[index] << " " << max_buf_global[index];
                }
                *ss << endl;
              }
            }
          }

          // cleanup...

          delete[] wgt_buf_global;
          delete[] avg_buf_global;
          delete[] rms_buf_global;
          delete[] min_buf_global;
          delete[] max_buf_global;

        }

        // advance rank0 for the next one...

        ++rank0;
        if (rank0 == mpi_size)
          rank0 = 0;

        // and rest the averages if reset_after_write...

        if (reset_after_write)
          resetBuffers();

      }
    }

    if (wvar.first == "")
      delete[] wvar_ptr;
  }

  void doSMDProbe(const int step,const double time) {

    static int rank0 = 0; // used to round-robin the writing...

    // conditional variable...
    CtiRegister::CtiData * cvar_data = cvar.second;
    if (cvar_data == NULL) cvar_data = CtiRegister::getUnregisteredCtiData(cvar.first);
    assert(cvar_data);

    // wgt...
    double *wvar_ptr = NULL;
    CtiRegister::CtiData * wvar_data = wvar.second;
    if (wvar_data != NULL) {
      // registered data
      wvar_ptr = wvar_data->getDNptr();
    }
    else {
      if (wvar.first != "") {
        // unregistered data
        wvar_data = CtiRegister::getUnregisteredCtiData(wvar.first);
        assert(wvar_data);
        wvar_ptr = wvar_data->getDNptr();
      }
      else {
        // constant wgt (default)
        wvar_ptr = new double[cvar_data->size()];
        for (int ii = 0; ii < cvar_data->size(); ++ii) wvar_ptr[ii] = 1.0;
      }
    }
    assert(wvar_ptr);

    // if geom is present, then the user is restricting the analysis using
    // the pointIsInside capability of the geom...
    bool * is_inside = NULL;
    if (geom) {
      // for now, must be cv data, because we don't know how to access coordinate
      // data otherwise. This could be made part of the CONDITIONAL_PROBE...
      CtiRegister::CtiData * x_cv_data;
      if (cvar_data->getUnindexedTopology() == BF_DATA) {
        const string zone_name = CtiRegister::getFirstSpecifier(cvar.first);
        x_cv_data = CtiRegister::getRegisteredCtiData(zone_name+":x_bf");
      }
      else if (cvar_data->getName() == "cht") {
        //cvar_data->getX(x_cv_data);  // doesn't support ctiData * 
        x_cv_data = CtiRegister::getRegisteredCtiData("cht:x_no");
      }
      else x_cv_data = CtiRegister::getRegisteredCtiData("x_cv");

      assert(x_cv_data);
      assert(x_cv_data->getType() == DN3_DATA);
//      if (x_cv_data->size() != cvar_data->size()) {
//        if (mpi_rank == 0) {
//          cout << "## x_cv_data->size(): " << x_cv_data->size() << endl;
//           cout << "## cvar_data->size(): " << cvar_data->size() << endl; 
//        }
//      }
      assert(x_cv_data->size() == cvar_data->size());
      is_inside = new bool[cvar_data->size()];
      double x_cv[3];
      for (int i = 0; i < cvar_data->size(); ++i) {
      	FOR_J3 x_cv[j] = x_cv_data->dn3(i,j);
      	is_inside[i] = geom->pointIsInside(x_cv);
      }
    }

    double my_range[2] = { 1.0E+20, 1.0E+20 };
    int my_count = 0;
    for (int i = 0; i < cvar_data->size(); ++i) {
      if ((is_inside == NULL)||(is_inside[i])) {
        ++my_count;
        const double c = cvar_data->dn(i);
        if (c < my_range[0]) my_range[0] = c;
        if (-c < my_range[1]) my_range[1] = -c;
      }
    }
    if (!b_range) {
      MPI_Allreduce(my_range,range,2,MPI_DOUBLE,MPI_MIN,mpi_comm);
      range[1] = -range[1]; // flip the max
      b_range = true;
    }
    my_range[1] = -my_range[1];

    if (range[0] == range[1]) {
      CWARN("Conditional variable " << cvar.first << " is constant. Skipping...");
      return;
    }

    if ( my_count == 0) {

      // due to the presence of the geom restriction, it is possible that this
      // rank has no participating data at all..

      ibin_min  = 0;
      ibin_max  = -1;

    } else {

      if (my_range[0] < range[0]) {
        // below lower bound...
        ibin_min = 0;
      }
      else if (my_range[0] > range[1]) {
        // above upper bound...
        ibin_min = nbin+1;
      }
      else {
        // should be [1:nbin] inclusive. need min operation to include dn == range[1] in last valid bin...
        ibin_min = min(nbin,1+(int)floor(double(nbin)*(my_range[0]-range[0])/(range[1]-range[0])));
      }
      if (my_range[1] < range[0]) {
        // below lower bound...
        ibin_max = 0;
      }
      else if (my_range[1] > range[1]) {
        // above upper bound...
        ibin_max = nbin+1;
      }
      else {
        // should be [1:nbin] inclusive. need min operation to include dn == range[1] in last valid bin...
        ibin_max = min(nbin,1+(int)floor(double(nbin)*(my_range[1]-range[0])/(range[1]-range[0])));
      }
      //cout << mpi_rank << " " << ibin_min << ":" << ibin_max << endl; cout.flush();
    }

    // for each variable we are conditionally probing, we add the avg,rms,min,max

    const int nvar = var_vec.size();
    assert(nvar > 0);
    int vol_var= -1;
    int sur_var= -1;
    for (int i = 0; i < nvar; ++i) {
      if (var_vec[i].first == "vof" || var_vec[i].first == "vof_avg") 
        vol_var = i;
      if (var_vec[i].first == "lsd" || var_vec[i].first == "lsd_avg") 
        sur_var = i;
    }
    
    assert(vol_var >= 0);
    assert(sur_var >= 0);
    

    if (cvar_2.first == "") {

      // 1d conditional probe...

      if (wgt_buf == NULL) {

        assert(avg_buf == NULL);
        assert(rms_buf == NULL);
        assert(min_buf == NULL);
        assert(max_buf == NULL);

        // the size of these buffers is nbin+2 per scalar nvar...

        //wgt_buf = new double[(nbin+2)*3]; // stores wgt, wgt*c and wgt*c*c
        //avg_buf = new double[(nbin+2)*nvar];
        //rms_buf = new double[(nbin+2)*nvar];
        //min_buf = new double[(nbin+2)*nvar];
        //max_buf = new double[(nbin+2)*nvar];
        wgt_buf = new double[(ibin_max-ibin_min+1)*3]; // stores wgt, wgt*c and wgt*c*c
        avg_buf = new double[(ibin_max-ibin_min+1)*nvar];
        rms_buf = new double[(ibin_max-ibin_min+1)*nvar];
        min_buf = new double[(ibin_max-ibin_min+1)*nvar];
        max_buf = new double[(ibin_max-ibin_min+1)*nvar];

        resetBuffers();

      }

      // set the bin index for each piece of data...

      int * ibin_buf = new int[cvar_data->size()];
      for (int i = 0; i < cvar_data->size(); ++i) {
        if ((is_inside == NULL)||(is_inside[i])) {
          const double c = cvar_data->dn(i);
          if (c < range[0]) {
            // below lower bound...
            ibin_buf[i] = 0;
          }
          else if (c > range[1]) {
            // above upper bound...
            ibin_buf[i] = nbin+1;
          }
          else {
            // should be [1:nbin] inclusive. need min operation to include dn == range[1] in last valid bin...
            ibin_buf[i] = min(nbin,1+(int)floor(double(nbin)*(c-range[0])/(range[1]-range[0])));
          }
          // and accumulate in wgt_buf...
          //wgt_buf[ibin_buf[i]]            += wvar_ptr[i];
          //wgt_buf[(nbin+2)+ibin_buf[i]]   += wvar_ptr[i]*c;
          //wgt_buf[(nbin+2)*2+ibin_buf[i]] += wvar_ptr[i]*c*c;
          assert((ibin_buf[i] >= ibin_min)&&(ibin_buf[i] <= ibin_max));
          wgt_buf[ibin_buf[i]-ibin_min]                         += wvar_ptr[i];
          wgt_buf[(ibin_max-ibin_min+1)+ibin_buf[i]-ibin_min]   += wvar_ptr[i]*c;
          wgt_buf[(ibin_max-ibin_min+1)*2+ibin_buf[i]-ibin_min] += wvar_ptr[i]*c*c;
        }
      }

      // now cycle through and add our stats into the bufs...

      for (int ivar = 0; ivar < nvar; ++ivar) {

        CtiRegister::CtiData * var_data = var_vec[ivar].second;
        if (var_data == NULL) var_data = CtiRegister::getUnregisteredCtiData(var_vec[ivar].first);
        assert(var_data);
        assert(var_data->size() == cvar_data->size());
        for (int i = 0; i < var_data->size(); ++i) {
          if ((is_inside == NULL)||(is_inside[i])) {
            const double v = var_data->dn(i);
            const int ibin = ibin_buf[i];
            //avg_buf[ivar*(nbin+2)+ibin] += wvar_ptr[i]*v;
            //rms_buf[ivar*(nbin+2)+ibin] += wvar_ptr[i]*v*v;
            //min_buf[ivar*(nbin+2)+ibin] = min(min_buf[ivar*(nbin+2)+ibin],v);
            //max_buf[ivar*(nbin+2)+ibin] = max(max_buf[ivar*(nbin+2)+ibin],v);
            assert((ibin_buf[i] >= ibin_min)&&(ibin_buf[i] <= ibin_max));
            avg_buf[ivar*(ibin_max-ibin_min+1)+ibin-ibin_min] += wvar_ptr[i]*v;
            rms_buf[ivar*(ibin_max-ibin_min+1)+ibin-ibin_min] += wvar_ptr[i]*v*v;
            min_buf[ivar*(ibin_max-ibin_min+1)+ibin-ibin_min] = min(min_buf[ivar*(ibin_max-ibin_min+1)+ibin-ibin_min],v);
            max_buf[ivar*(ibin_max-ibin_min+1)+ibin-ibin_min] = max(max_buf[ivar*(ibin_max-ibin_min+1)+ibin-ibin_min],v);
          }
        }

      }
      delete[] ibin_buf;
      if (is_inside) delete[] is_inside;

      if ((step == -1)||(step%write_interval == 0)) {

        // reduce and write into probe_rank0...

        // for now, we reduce EVERYTHING. probably a better way to do this some day...

        double * wgt_buf_global = NULL;
        double * avg_buf_global = NULL;
        double * rms_buf_global = NULL;
        double * min_buf_global = NULL;
        double * max_buf_global = NULL;

        if (mpi_rank == rank0) {

          wgt_buf_global = new double[(nbin+2)*3];
          avg_buf_global = new double[(nbin+2)*nvar];
          rms_buf_global = new double[(nbin+2)*nvar];
          min_buf_global = new double[(nbin+2)*nvar];
          max_buf_global = new double[(nbin+2)*nvar];

        }

        // and do reductions...

        //MPI_Reduce(wgt_buf,wgt_buf_global,(nbin+2)*3,MPI_DOUBLE,MPI_SUM,rank0,mpi_comm);
        //MPI_Reduce(avg_buf,avg_buf_global,(nbin+2)*nvar,MPI_DOUBLE,MPI_SUM,rank0,mpi_comm);
        //MPI_Reduce(rms_buf,rms_buf_global,(nbin+2)*nvar,MPI_DOUBLE,MPI_SUM,rank0,mpi_comm);
        //MPI_Reduce(min_buf,min_buf_global,(nbin+2)*nvar,MPI_DOUBLE,MPI_MIN,rank0,mpi_comm);
        //MPI_Reduce(max_buf,max_buf_global,(nbin+2)*nvar,MPI_DOUBLE,MPI_MAX,rank0,mpi_comm);

        // init global bufs...
        if (mpi_rank == rank0) {
          for (int i = 0, limit = (nbin+2)*3; i < limit; ++i) {
            wgt_buf_global[i] = 0.0;
          }
          for (int i = 0, limit = (nbin+2)*var_vec.size(); i < limit; ++i) {
            avg_buf_global[i] = 0.0;
            rms_buf_global[i] = 0.0;
            min_buf_global[i] = 1.0E+20;
            max_buf_global[i] = -1.0E+20;
          }
        }

        // gather local ranges...
        int my_ibin_range[2] = {ibin_min,ibin_max};
        int (*ibin_range)[2] = NULL;
        if (mpi_rank == rank0)
          ibin_range = new int[mpi_size][2];
        MPI_Gather(my_ibin_range,2,MPI_INT,(int*)ibin_range,2,MPI_INT,rank0,mpi_comm);
        int * count = NULL;
        int * disp = NULL;
        double *recv_buf_double = NULL;
        if (mpi_rank == rank0) {
          count = new int[mpi_size];
          FOR_RANK count[rank] = ibin_range[rank][1]-ibin_range[rank][0]+1;
          disp  = new int[mpi_size];
          disp[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank)
            disp[rank] = disp[rank-1] + count[rank-1];
          int recv_count_sum = disp[mpi_size-1] + count[mpi_size-1]; // int8?
          recv_buf_double = new double[max(3,nvar)*recv_count_sum];
        }

        // now gather the double data and populate
        if (mpi_rank == rank0) {
          FOR_RANK {
            count[rank] *= 3;
            disp[rank]  *= 3;
          //cout << rank << " " << ibin_range[rank][0] << ":" << ibin_range[rank][1] << " " <<  count[rank] << " " << disp[rank] << endl;
          }
        }
        MPI_Gatherv(wgt_buf,(ibin_max-ibin_min+1)*3,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            FOR_I3 {
              for (int ii = 0; ii < count[rank]/3; ++ii) {
                wgt_buf_global[(nbin+2)*i+ii+ibin_range[rank][0]] += recv_buf_double[disp[rank]+ii+(count[rank]/3)*i];
              }
            }
          }
        }
        if (mpi_rank == rank0) {
          FOR_RANK {
            count[rank] = (count[rank]/3)*nvar;
            disp[rank]  = (disp[rank]/3)*nvar;
          }
        }
        MPI_Gatherv(avg_buf,(ibin_max-ibin_min+1)*nvar,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            for (int ivar = 0; ivar < nvar; ++ivar) {
              for (int ii = 0; ii < count[rank]/nvar; ++ii) {
                avg_buf_global[(nbin+2)*ivar+ii+ibin_range[rank][0]] += recv_buf_double[disp[rank]+ii+(count[rank]/nvar)*ivar];
              }
            }
          }
        }
        MPI_Gatherv(rms_buf,(ibin_max-ibin_min+1)*nvar,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            for (int ivar = 0; ivar < nvar; ++ivar) {
              for (int ii = 0; ii < count[rank]/nvar; ++ii) {
                rms_buf_global[(nbin+2)*ivar+ii+ibin_range[rank][0]] += recv_buf_double[disp[rank]+ii+(count[rank]/nvar)*ivar];
              }
            }
          }
        }
        MPI_Gatherv(min_buf,(ibin_max-ibin_min+1)*nvar,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            for (int ivar = 0; ivar < nvar; ++ivar) {
              for (int ii = 0; ii < count[rank]/nvar; ++ii) {
                min_buf_global[(nbin+2)*ivar+ii+ibin_range[rank][0]] = min(min_buf_global[(nbin+2)*ivar+ii+ibin_range[rank][0]],recv_buf_double[disp[rank]+ii+(count[rank]/nvar)*ivar]);
              }
            }
          }
        }
        MPI_Gatherv(max_buf,(ibin_max-ibin_min+1)*nvar,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            for (int ivar = 0; ivar < nvar; ++ivar) {
              for (int ii = 0; ii < count[rank]/nvar; ++ii) {
                max_buf_global[(nbin+2)*ivar+ii+ibin_range[rank][0]] = max(max_buf_global[(nbin+2)*ivar+ii+ibin_range[rank][0]],recv_buf_double[disp[rank]+ii+(count[rank]/nvar)*ivar]);
              }
            }
          }
        }
        if (mpi_rank == rank0) {
          delete[] ibin_range;
          delete[] count;
          delete[] disp;
          delete[] recv_buf_double;
        }

        if (mpi_rank == rank0) {

          char filename[128];
          if (step == -1) {
            sprintf(filename,"%s.scp",name.c_str());
          }
          else {
            sprintf(filename,"%s.%08d.scp",name.c_str(),step);
          }

          cpFileList.push_back(pair<string,stringstream*>(string(filename),new stringstream()));
          stringstream* ss = cpFileList.back().second;

          // add the header...

          int index_vol = -1;
          int index_sur = -1;
          for (int i = 0; i < nvar; ++i) {
            if (var_vec[i].first == "vof" || var_vec[i].first == "vof_avg") 
              index_vol = i;
            if (var_vec[i].first == "lsd" || var_vec[i].first == "lsd_avg") 
              index_sur = i;
          }

          assert(index_vol >= 0);
          assert(index_sur >= 0);

          *ss << setprecision(8);
          *ss << "# step " << step << " time " << time << endl;
          if (wvar.first == "")
            *ss << "# 1:bin 2:count 3:" << cvar.first << "_avg 4:" << cvar.first << "_rms";
          else
            *ss << "# 1:bin 2:total 3:" << cvar.first << "_avg 4:" << cvar.first << "_rms";
          for (int i = 0; i < int(var_vec.size()); ++i)
            *ss << " " << 5+i*4 << ":" << var_vec[i].first << "_avg " <<
              6+i*4 << ":" << var_vec[i].first << "_rms " <<
              7+i*4 << ":" << var_vec[i].first << "_min " <<
              8+i*4 << ":" << var_vec[i].first << "_max ";
          *ss << 9+(nvar-1)*4 << ":" <<  "smd_avg ";
          *ss << endl;

          for (int ibin = 0; ibin < nbin+2; ++ibin) {
            if (wgt_buf_global[ibin] > 0.0) {
              *ss << ibin << " " << wgt_buf_global[ibin];
              // this bin contains valid data. Report the bin center using the
              // c average and rms...
              const double cavg = wgt_buf_global[(nbin+2)+ibin]/wgt_buf_global[ibin];
              const double crms = sqrt(max(0.0,wgt_buf_global[(nbin+2)*2+ibin]/wgt_buf_global[ibin] - cavg*cavg));
              *ss << " " << cavg << " " << crms;
              for (int ivar = 0; ivar < nvar; ++ivar) {
                const double avg = avg_buf_global[(nbin+2)*ivar+ibin]/wgt_buf_global[ibin];
                const double rms = sqrt(max(0.0,rms_buf_global[(nbin+2)*ivar+ibin]/wgt_buf_global[ibin] - avg*avg));
                *ss << " " << avg << " " << rms << " " << min_buf_global[(nbin+2)*ivar+ibin] << " " << max_buf_global[(nbin+2)*ivar+ibin];
              }
              if (fabs(avg_buf_global[(nbin+2)*sur_var+ibin]) > 1.0E-16) 
                *ss << " "<<  6*avg_buf_global[(nbin+2)*vol_var+ibin]/avg_buf_global[(nbin+2)*sur_var+ibin];
              else 
                *ss << " "<<  0.0;

              *ss << endl;
            }
          }

          // cleanup...

          delete[] wgt_buf_global;
          delete[] avg_buf_global;
          delete[] rms_buf_global;
          delete[] min_buf_global;
          delete[] max_buf_global;

        }

        // advance rank0 for the next one...

        ++rank0;
        if (rank0 == mpi_size)
          rank0 = 0;

        // and rest the averages if reset_after_write...

        if (reset_after_write)
          resetBuffers();

      }

    }
    else {

      // 2d conditional probe...

      CtiRegister::CtiData * cvar_2_data = cvar_2.second;
      if (cvar_2_data == NULL) cvar_2_data = CtiRegister::getUnregisteredCtiData(cvar_2.first);
      assert(cvar_2_data);

      double my_range_2[2] = { 1.0E+20, 1.0E+20 };
      int my_count_2 = 0;
      for (int i = 0; i < cvar_2_data->size(); ++i) {
        if ((is_inside == NULL)||(is_inside[i])) {
          ++my_count_2;
          const double c = cvar_2_data->dn(i);
          if (c < my_range_2[0]) my_range_2[0] = c;
          if (-c < my_range_2[1]) my_range_2[1] = -c;
        }
      }
      if (!b_range_2) {
        MPI_Allreduce(my_range_2,range_2,2,MPI_DOUBLE,MPI_MIN,mpi_comm);
        range_2[1] = -range_2[1]; // flip the max
        b_range_2 = true;
      }
      my_range_2[1] = -my_range_2[1];

      if (range_2[0] == range_2[1]) {
        CWARN("Conditional variable " << cvar_2.first << " is constant. Skipping...");
        return;
      }

      if ( my_count_2 == 0) {

        // due to the presence of the geom restriction, it is possible that this
        // rank has no participating data at all..

        ibin_min_2  = 0;
        ibin_max_2  = -1;

      } else {
        if (my_range_2[0] < range_2[0]) {
          // below lower bound...
          ibin_min_2 = 0;
        }
        else if (my_range_2[0] > range_2[1]) {
          // above upper bound...
          ibin_min_2 = nbin_2+1;
        }
        else {
          // should be [1:nbin_2] inclusive. need min operation to include dn == range_2[1] in last valid bin...
          ibin_min_2 = min(nbin_2,1+(int)floor(double(nbin_2)*(my_range_2[0]-range_2[0])/(range_2[1]-range_2[0])));
        }
        if (my_range_2[1] < range_2[0]) {
          // below lower bound...
          ibin_max_2 = 0;
        }
        else if (my_range_2[1] > range_2[1]) {
          // above upper bound...
          ibin_max_2 = nbin_2+1;
        }
        else {
          // should be [1:nbin_2] inclusive. need min operation to include dn == range_2[1] in last valid bin...
          ibin_max_2 = min(nbin_2,1+(int)floor(double(nbin_2)*(my_range_2[1]-range_2[0])/(range_2[1]-range_2[0])));
        }
      }

      if (wgt_buf == NULL) {

        assert(avg_buf == NULL);
        assert(rms_buf == NULL);
        assert(min_buf == NULL);
        assert(max_buf == NULL);

        // the size of these buffers is (nbin+2)*(nbin_2+2) per scalar nvar...

        //wgt_buf = new double[(nbin+2)*(nbin_2+2)*5]; // stores wgt, wgt*c, wgt*c*c, wgt*c_2 and wgt*c_2*c_2
        //avg_buf = new double[(nbin+2)*(nbin_2+2)*nvar];
        //rms_buf = new double[(nbin+2)*(nbin_2+2)*nvar];
        //min_buf = new double[(nbin+2)*(nbin_2+2)*nvar];
        //max_buf = new double[(nbin+2)*(nbin_2+2)*nvar];
        wgt_buf = new double[(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*5]; // stores wgt,wgt*c,wgt*c*c,wgt*c_2,wgt*c_2*c_2
        avg_buf = new double[(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*nvar];
        rms_buf = new double[(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*nvar];
        min_buf = new double[(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*nvar];
        max_buf = new double[(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*nvar];

        resetBuffers();

      }

      // set the bin index for each piece of data...

      assert(cvar_2_data->size() == cvar_data->size());
      int * ibin_buf = new int[cvar_data->size()];
      int * jbin_buf = new int[cvar_2_data->size()];
      for (int ii = 0; ii < cvar_data->size(); ++ii) {
        if ((is_inside == NULL)||(is_inside[ii])) {
          const double c = cvar_data->dn(ii);
          if (c < range[0]) {
            // below lower bound...
            ibin_buf[ii] = 0;
          }
          else if (c > range[1]) {
            // above upper bound...
            ibin_buf[ii] = nbin+1;
          }
          else {
            // should be [1:nbin] inclusive. need min operation to include dn == range[1] in last valid bin...
            ibin_buf[ii] = min(nbin,1+(int)floor(double(nbin)*(c-range[0])/(range[1]-range[0])));
          }
          const double c_2 = cvar_2_data->dn(ii);
          if (c_2 < range_2[0]) {
            // below lower bound...
            jbin_buf[ii] = 0;
          }
          else if (c_2 > range_2[1]) {
            // above upper bound...
            jbin_buf[ii] = nbin_2+1;
          }
          else {
            // should be [1:nbin_2] inclusive. need min operation to include dn == range_2[1] in last valid bin...
            jbin_buf[ii] = min(nbin_2,1+(int)floor(double(nbin_2)*(c_2-range_2[0])/(range_2[1]-range_2[0])));
          }
          // and accumulate in wgt_buf...
          //const int kbin = ibin_buf[ii]*(nbin_2+2)+jbin_buf[ii];
          //wgt_buf[                      kbin] += wvar_ptr[ii];
          //wgt_buf[(nbin+2)*(nbin_2+2)  +kbin] += wvar_ptr[ii]*c;
          //wgt_buf[(nbin+2)*(nbin_2+2)*2+kbin] += wvar_ptr[ii]*c*c;
          //wgt_buf[(nbin+2)*(nbin_2+2)*3+kbin] += wvar_ptr[ii]*c_2;
          //wgt_buf[(nbin+2)*(nbin_2+2)*4+kbin] += wvar_ptr[ii]*c_2*c_2;
          assert((ibin_buf[ii] >= ibin_min)&&(ibin_buf[ii] <= ibin_max));
          assert((jbin_buf[ii] >= ibin_min_2)&&(jbin_buf[ii] <= ibin_max_2));
          const int index = (ibin_buf[ii]-ibin_min)*(ibin_max_2-ibin_min_2+1)+jbin_buf[ii]-ibin_min_2;
          wgt_buf[                                                  index] += wvar_ptr[ii];
          wgt_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)  +index] += wvar_ptr[ii]*c;
          wgt_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*2+index] += wvar_ptr[ii]*c*c;
          wgt_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*3+index] += wvar_ptr[ii]*c_2;
          wgt_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*4+index] += wvar_ptr[ii]*c_2*c_2;
        }
      }

      // now cycle through and add our stats into the bufs...

      for (int ivar = 0; ivar < nvar; ++ivar) {

        CtiRegister::CtiData * var_data = var_vec[ivar].second;
        if (var_data == NULL) var_data = CtiRegister::getUnregisteredCtiData(var_vec[ivar].first);
        assert(var_data);
        assert(var_data->size() == cvar_data->size());
        for (int ii = 0; ii < var_data->size(); ++ii) {
          if ((is_inside == NULL)||(is_inside[ii])) {
            const double v = var_data->dn(ii);
            //const int index = ivar*(nbin+2)*(nbin_2+2)+ibin_buf[ii]*(nbin_2+2)+jbin_buf[ii];
            //avg_buf[index] += wvar_ptr[ii]*v;
            //rms_buf[index] += wvar_ptr[ii]*v*v;
            //min_buf[index] = min(min_buf[index],v);
            //max_buf[index] = max(max_buf[index],v);
            assert((ibin_buf[ii] >= ibin_min)&&(ibin_buf[ii] <= ibin_max));
            assert((jbin_buf[ii] >= ibin_min_2)&&(jbin_buf[ii] <= ibin_max_2));
            const int index = (ibin_buf[ii]-ibin_min)*(ibin_max_2-ibin_min_2+1)+jbin_buf[ii]-ibin_min_2;
            avg_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*ivar+index] += wvar_ptr[ii]*v;
            rms_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*ivar+index] += wvar_ptr[ii]*v*v;
            min_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*ivar+index] = min(min_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*ivar+index],v);
            max_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*ivar+index] = max(max_buf[(ibin_max_2-ibin_min_2+1)*(ibin_max-ibin_min+1)*ivar+index],v);
          }
        }

      }
      delete[] ibin_buf;
      delete[] jbin_buf;
      if (is_inside) delete[] is_inside;

      if ((step == -1)||(step%write_interval == 0)) {

        // reduce and write into probe_rank0...

        // for now, we reduce EVERYTHING. probably a better way to do this some day...

        double * wgt_buf_global = NULL;
        double * avg_buf_global = NULL;
        double * rms_buf_global = NULL;
        double * min_buf_global = NULL;
        double * max_buf_global = NULL;

        if (mpi_rank == rank0) {

          wgt_buf_global = new double[(nbin+2)*(nbin_2+2)*5];
          avg_buf_global = new double[(nbin+2)*(nbin_2+2)*nvar];
          rms_buf_global = new double[(nbin+2)*(nbin_2+2)*nvar];
          min_buf_global = new double[(nbin+2)*(nbin_2+2)*nvar];
          max_buf_global = new double[(nbin+2)*(nbin_2+2)*nvar];

        }

        // and do reductions...

        //MPI_Reduce(wgt_buf,wgt_buf_global,(nbin+2)*(nbin_2+2)*5,MPI_DOUBLE,MPI_SUM,rank0,mpi_comm);
        //MPI_Reduce(avg_buf,avg_buf_global,(nbin+2)*(nbin_2+2)*nvar,MPI_DOUBLE,MPI_SUM,rank0,mpi_comm);
        //MPI_Reduce(rms_buf,rms_buf_global,(nbin+2)*(nbin_2+2)*nvar,MPI_DOUBLE,MPI_SUM,rank0,mpi_comm);
        //MPI_Reduce(min_buf,min_buf_global,(nbin+2)*(nbin_2+2)*nvar,MPI_DOUBLE,MPI_MIN,rank0,mpi_comm);
        //MPI_Reduce(max_buf,max_buf_global,(nbin+2)*(nbin_2+2)*nvar,MPI_DOUBLE,MPI_MAX,rank0,mpi_comm);

        // init global bufs...
        if (mpi_rank == rank0) {
          for (int i = 0, limit = (nbin+2)*(nbin_2+2)*5; i < limit; ++i) {
            wgt_buf_global[i] = 0.0;
          }
          for (int i = 0, limit = (nbin+2)*(nbin_2+2)*var_vec.size(); i < limit; ++i) {
            avg_buf_global[i] = 0.0;
            rms_buf_global[i] = 0.0;
            min_buf_global[i] = 1.0E+20;
            max_buf_global[i] = -1.0E+20;
          }
        }

        // gather local ranges...
        int my_ibin_range[4] = {ibin_min,ibin_max,ibin_min_2,ibin_max_2};
        int (*ibin_range)[4] = NULL;
        if (mpi_rank == rank0)
          ibin_range = new int[mpi_size][4];
        MPI_Gather(my_ibin_range,4,MPI_INT,(int*)ibin_range,4,MPI_INT,rank0,mpi_comm);
        int * count = NULL;
        int * disp = NULL;
        double *recv_buf_double = NULL;
        if (mpi_rank == rank0) {
          count = new int[mpi_size];
          FOR_RANK count[rank] = (ibin_range[rank][1]-ibin_range[rank][0]+1)*(ibin_range[rank][3]-ibin_range[rank][2]+1);
          disp  = new int[mpi_size];
          disp[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank)
            disp[rank] = disp[rank-1] + count[rank-1];
          int recv_count_sum = disp[mpi_size-1] + count[mpi_size-1]; // int8?
          recv_buf_double = new double[max(5,nvar)*recv_count_sum];
        }

        // now gather the double data and populate
        if (mpi_rank == rank0) {
          FOR_RANK {
            count[rank] *= 5;
            disp[rank]  *= 5;
            //cout << rank << " " << count[rank] << " " << disp[rank] << " " << ibin_range[rank][0] << " " << ibin_range[rank][1] << " " << ibin_range[rank][2] << " " << ibin_range[rank][3] << endl;
          }
        }
        MPI_Gatherv(wgt_buf,(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*5,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            FOR_I5 {
              for (int ii = 0; ii < count[rank]/5; ++ii) {
                const int ix = ibin_range[rank][0]+ii/(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int iy = ibin_range[rank][2]+ii%(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int index = (nbin_2+2)*ix+iy;
                wgt_buf_global[(nbin_2+2)*(nbin+2)*i+index] += recv_buf_double[disp[rank]+ii+(count[rank]/5)*i];
              }
            }
          }
        }
        if (mpi_rank == rank0) {
          FOR_RANK {
            count[rank] = (count[rank]/5)*nvar;
            disp[rank]  = (disp[rank]/5)*nvar;
          }
        }
        MPI_Gatherv(avg_buf,(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*nvar,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            for (int ivar = 0; ivar < nvar; ++ivar) {
              for (int ii = 0; ii < count[rank]/nvar; ++ii) {
                const int ix = ibin_range[rank][0]+ii/(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int iy = ibin_range[rank][2]+ii%(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int index = (nbin_2+2)*ix+iy;
                avg_buf_global[(nbin_2+2)*(nbin+2)*ivar+index] += recv_buf_double[disp[rank]+ii+(count[rank]/nvar)*ivar];
              }
            }
          }
        }
        MPI_Gatherv(rms_buf,(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*nvar,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            for (int ivar = 0; ivar < nvar; ++ivar) {
              for (int ii = 0; ii < count[rank]/nvar; ++ii) {
                const int ix = ibin_range[rank][0]+ii/(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int iy = ibin_range[rank][2]+ii%(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int index = (nbin_2+2)*ix+iy;
                rms_buf_global[(nbin_2+2)*(nbin+2)*ivar+index] += recv_buf_double[disp[rank]+ii+(count[rank]/nvar)*ivar];
              }
            }
          }
        }
        MPI_Gatherv(min_buf,(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*nvar,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            for (int ivar = 0; ivar < nvar; ++ivar) {
              for (int ii = 0; ii < count[rank]/nvar; ++ii) {
                const int ix = ibin_range[rank][0]+ii/(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int iy = ibin_range[rank][2]+ii%(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int index = (nbin_2+2)*ix+iy;
                min_buf_global[(nbin_2+2)*(nbin+2)*ivar+index] = min(min_buf_global[(nbin_2+2)*(nbin+2)*ivar+index],recv_buf_double[disp[rank]+ii+(count[rank]/nvar)*ivar]);
              }
            }
          }
        }
        MPI_Gatherv(max_buf,(ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*nvar,MPI_DOUBLE,recv_buf_double,count,disp,MPI_DOUBLE,rank0,mpi_comm);
        if (mpi_rank == rank0) {
          FOR_RANK {
            for (int ivar = 0; ivar < nvar; ++ivar) {
              for (int ii = 0; ii < count[rank]/nvar; ++ii) {
                const int ix = ibin_range[rank][0]+ii/(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int iy = ibin_range[rank][2]+ii%(ibin_range[rank][3]-ibin_range[rank][2]+1);
                const int index = (nbin_2+2)*ix+iy;
                max_buf_global[(nbin_2+2)*(nbin+2)*ivar+index] = max(max_buf_global[(nbin_2+2)*(nbin+2)*ivar+index],recv_buf_double[disp[rank]+ii+(count[rank]/nvar)*ivar]);
              }
            }
          }
        }
        if (mpi_rank == rank0) {
          delete[] ibin_range;
          delete[] count;
          delete[] disp;
          delete[] recv_buf_double;
        }

        if (mpi_rank == rank0) {

          char filename[128];
          if (step == -1) {
            sprintf(filename,"%s.scp",name.c_str());
          }
          else {
            sprintf(filename,"%s.%08d.scp",name.c_str(),step);
          }

          cpFileList.push_back(pair<string,stringstream*>(string(filename),new stringstream()));
          stringstream* ss = cpFileList.back().second;

          // add the header...


          *ss << setprecision(8);
          *ss << "# step " << step << " time " << time << endl;
          if (wvar.first == "")
            *ss << "# 1:ibin 2:jbin 3:count 4:";
          else
            *ss << "# 1:ibin 2:jbin 3:total 4:";
          *ss << cvar.first << "_avg 5:" << cvar.first << "_rms 6:"
              << cvar_2.first << "_avg 7:" << cvar_2.first << "_rms";
          for (int i = 0; i < int(var_vec.size()); ++i)
            *ss << " " << 8+i*4 << ":" << var_vec[i].first << "_avg " <<
              9+i*4 << ":" << var_vec[i].first << "_rms " <<
              10+i*4 << ":" << var_vec[i].first << "_min " <<
              11+i*4 << ":" << var_vec[i].first << "_max ";
          *ss << 12 + (nvar-1)*4 << ":" << "smd_avg";
          *ss << endl;

          for (int ibin = 0; ibin < nbin+2; ++ibin) {
            for (int jbin = 0; jbin < nbin_2+2; ++jbin) {
              const int kbin = ibin*(nbin_2+2)+jbin;
              if (wgt_buf_global[kbin] > 0.0) {
                *ss << ibin << " " << jbin << " " << wgt_buf_global[kbin];
                // this bin contains valid data. Report the bin center using the
                // c average and rms...
                const double c_avg = wgt_buf_global[(nbin+2)*(nbin_2+2)+kbin]/wgt_buf_global[kbin];
                const double c_rms = sqrt(max(0.0,wgt_buf_global[(nbin+2)*(nbin_2+2)*2+kbin]/wgt_buf_global[kbin] - c_avg*c_avg));
                *ss << " " << c_avg << " " << c_rms;
                const double c_2_avg = wgt_buf_global[(nbin+2)*(nbin_2+2)*3+kbin]/wgt_buf_global[kbin];
                const double c_2_rms = sqrt(max(0.0,wgt_buf_global[(nbin+2)*(nbin_2+2)*4+kbin]/wgt_buf_global[kbin] - c_2_avg*c_2_avg));
                *ss << " " << c_2_avg << " " << c_2_rms;
                for (int ivar = 0; ivar < nvar; ++ivar) {
                  const int index = ivar*(nbin+2)*(nbin_2+2)+kbin;
                  const double avg = avg_buf_global[index]/wgt_buf_global[kbin];
                  const double rms = sqrt(max(0.0,rms_buf_global[index]/wgt_buf_global[kbin] - avg*avg));
                  *ss << " " << avg << " " << rms << " " << min_buf_global[index] << " " << max_buf_global[index];
                }
                const int vol_index = vol_var*(nbin+2)*(nbin_2+2)+kbin;
                const int sur_index = sur_var*(nbin+2)*(nbin_2+2)+kbin;
                if ( fabs(avg_buf_global[sur_index]) > 1.0E-16) 
                  *ss << " " << 6.0*avg_buf_global[vol_index]/avg_buf_global[sur_index] << endl;
                else 
                  *ss << " " << 0.0 << endl; ;
              }
            }
          }

          // cleanup...

          delete[] wgt_buf_global;
          delete[] avg_buf_global;
          delete[] rms_buf_global;
          delete[] min_buf_global;
          delete[] max_buf_global;

        }

        // advance rank0 for the next one...

        ++rank0;
        if (rank0 == mpi_size)
          rank0 = 0;

        // and rest the averages if reset_after_write...

        if (reset_after_write)
          resetBuffers();

      }
    }

    if (wvar.first == "")
      delete[] wvar_ptr;
  }


  void resetBuffers() {

    assert(wgt_buf);
    assert(avg_buf);
    assert(rms_buf);
    assert(min_buf);
    assert(max_buf);

    if (cvar_2.first == "") {

      // 1d conditional probe...

      /*
      for (int i = 0, limit = (nbin+2)*3; i < limit; ++i) {
        wgt_buf[i] = 0.0;
      }

      for (int i = 0, limit = (nbin+2)*var_vec.size(); i < limit; ++i) {
        avg_buf[i] = 0.0;
        rms_buf[i] = 0.0;
        min_buf[i] = 1.0E+20;
        max_buf[i] = -1.0E+20;
      }
      */

      for (int i = 0, limit = (ibin_max-ibin_min+1)*3; i < limit; ++i) {
        wgt_buf[i] = 0.0;
      }

      for (int i = 0, limit = (ibin_max-ibin_min+1)*var_vec.size(); i < limit; ++i) {
        avg_buf[i] = 0.0;
        rms_buf[i] = 0.0;
        min_buf[i] = 1.0E+20;
        max_buf[i] = -1.0E+20;
      }

    }
    else {

      // 2d conditional probe...

      /*
      for (int i = 0, limit = (nbin+2)*(nbin_2+2)*5; i < limit; ++i) {
        wgt_buf[i] = 0.0;
      }

      for (int i = 0, limit = (nbin+2)*(nbin_2+2)*var_vec.size(); i < limit; ++i) {
        avg_buf[i] = 0.0;
        rms_buf[i] = 0.0;
        min_buf[i] = 1.0E+20;
        max_buf[i] = -1.0E+20;
      }
      */

      for (int i = 0, limit = (ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*5; i < limit; ++i) {
        wgt_buf[i] = 0.0;
      }

      for (int i = 0, limit = (ibin_max-ibin_min+1)*(ibin_max_2-ibin_min_2+1)*var_vec.size(); i < limit; ++i) {
        avg_buf[i] = 0.0;
        rms_buf[i] = 0.0;
        min_buf[i] = 1.0E+20;
        max_buf[i] = -1.0E+20;
      }

    }

  }

  void flush() {

    // any files to be written are stored in cpFileList. Could be on any ranks...

    for (list<pair<string,stringstream*> >::iterator iter = cpFileList.begin(); iter != cpFileList.end(); ++iter) {
      ofstream out_file;
      if ((sample_interval == -1) && (write_interval == -1)) {
        const string tmp_filename = MiscUtils::makeTmpPrefix(iter->first);
        out_file.open(tmp_filename.c_str(),ofstream::trunc);
        assert(out_file.is_open());
        out_file << iter->second->rdbuf();
        out_file.close();
        remove(iter->first.c_str());
        rename(tmp_filename.c_str(),iter->first.c_str());
      }
      else {
        out_file.open(iter->first.c_str(),ofstream::out);
        assert(out_file.is_open());
        out_file << iter->second->rdbuf();
        out_file.close();
      }
      delete iter->second;
    }
    cpFileList.clear();

  }

};

#endif
