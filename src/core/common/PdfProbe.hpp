#ifndef PDFPROBE_HPP
#define PDFPROBE_HPP

#include "SimpleGeom.hpp"

class PdfProbe {
public:

  string name,header;
  int write_interval;
  int sample_interval;
  bool reset_after_write;

  SimpleGeom * geom;
  pair<string,CtiRegister::CtiData*> var;

  pair<string,CtiRegister::CtiData*> wvar;

  int nbin,nbin_target;
  double range[2];
  bool b_range;
  double * count; // weighted count! (e.g. 1.0,volume,area_bf,...)

  list<pair<string,stringstream*> > pdfpFileList;

  PdfProbe() {

    name = "";
    write_interval  = -1; // used to determine whether param specified value or not
    sample_interval = -1;
    reset_after_write = false; // default is to accumulate stats indefinately

    geom = NULL;

    nbin = nbin_target = 10;   // default bin size
    b_range = false; // by default range will be set by first doProbe...
    range[0] = range[1] = 0.0;
    count = NULL;

    wvar.first  = "";   // unit weighting when we don't get wgt var...
    wvar.second = NULL;

  }

  ~PdfProbe() {


    if (geom) delete geom;

    for (list<pair<string,stringstream*> >::iterator iter = pdfpFileList.begin(); iter != pdfpFileList.end(); ++iter)
      delete iter->second;

    DELETE(count);

  }

  void doProbe(const int step,const double time) {

    static int rank0 = 0; // used to round-robin the writing...

    // conditional variable...
    CtiRegister::CtiData * var_data = var.second;
    if (var_data == NULL) var_data = CtiRegister::getUnregisteredCtiData(var.first);
    assert(var_data);

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
        wvar_ptr = new double[var_data->size()];
        for (int ii = 0; ii < var_data->size(); ++ii) wvar_ptr[ii] = 1.0;
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
      if (var_data->getUnindexedTopology() == BF_DATA) {
        // const size_t colon_pos = var.first.find_first_of(":");
        // const string delimiters = " =+-/*(^%$,";
        // const size_t start_pos = var.first.substr(0,colon_pos).find_last_of(delimiters);
        // const string zone_name = var.first.substr(start_pos+1,colon_pos-start_pos-1);
        const string zone_name = CtiRegister::getFirstSpecifier(var.first);
        x_cv_data = CtiRegister::getRegisteredCtiData(zone_name+":x_bf");
      }
      else x_cv_data = CtiRegister::getRegisteredCtiData("x_cv");

      assert(x_cv_data);
      assert(x_cv_data->getType() == DN3_DATA);
      assert(x_cv_data->size() == var_data->size());
      is_inside = new bool[var_data->size()];
      double x_cv[3];
      for (int i = 0; i < var_data->size(); ++i) {
	FOR_J3 x_cv[j] = x_cv_data->dn3(i,j);
	is_inside[i] = geom->pointIsInside(x_cv);
      }
    }

    if (count == NULL) {
      count = new double[nbin];
      resetBuffers();
    }

    if (!b_range) {

      // this is the first time, and b_range is not set, so set the range...
      double my_data_minmax[2] = { HUGE_VAL, HUGE_VAL };
      for (int i = 0; i < var_data->size(); ++i) {
        if ((is_inside == NULL)||(is_inside[i])) {
            my_data_minmax[0] = min(my_data_minmax[0],var_data->dn(i));
            my_data_minmax[1] = min(my_data_minmax[1],-var_data->dn(i));
          }
      }
      double data_minmax[2];
      MPI_Allreduce(my_data_minmax,data_minmax,2,MPI_DOUBLE,MPI_MIN,mpi_comm);

      double drange = 0.01*(-data_minmax[1]-data_minmax[0]);
      if (drange == 0.0) drange = 0.01; // avoid /0 errors
      range[0] = data_minmax[0] - drange;
      range[1] = -data_minmax[1] + drange;
      b_range = true;

      if (mpi_rank == rank0) cout << " > pdf probe range set to: " << range[0] << " " << range[1] << endl;
    }

    // add the data...

    for (int i = 0; i < var_data->size(); ++i) {
      if ((is_inside == NULL)||(is_inside[i])) {
        const int ibin = int((var_data->dn(i) - range[0])/(range[1]-range[0])*double(nbin));
        if ((ibin >= 0)&&(ibin < nbin))
          count[ibin] += wvar_ptr[i];
      }
    }

    if (is_inside) delete[] is_inside;

    if ((step == -1)||(step%write_interval == 0)) {

      // reduce and write into probe_rank0...

      // for now, we reduce EVERYTHING. probably a better way to do this some day...

      double * count_global = NULL;

      if (mpi_rank == rank0) {

        count_global = new double[nbin];

      }

      // and do reductions...

      MPI_Reduce(count,count_global,nbin,MPI_DOUBLE,MPI_SUM,rank0,mpi_comm);

      if (mpi_rank == rank0) {

        char filename[128];
        if (step == -1) {
          sprintf(filename,"%s.pdfp",name.c_str());
        }
        else {
          sprintf(filename,"%s.%08d.pdfp",name.c_str(),step);
        }

        pdfpFileList.push_back(pair<string,stringstream*>(string(filename),new stringstream()));
        stringstream* ss = pdfpFileList.back().second;

        // add the header...

        *ss << setprecision(8);
        *ss << "# step " << step << " time " << time << endl;

        if (wvar.first == "")
          *ss << "# 1:bin 2:count 3:bin center" << endl;
        else
          *ss << "# 1:bin 2:total 3:bin center" << endl;
        for (int ibin = 0; ibin < nbin; ++ibin)
          *ss << ibin << " " << count_global[ibin] << " " << range[0]+(range[1]-range[0])*(double(ibin)+0.5)/double(nbin) << endl;
        //*ss << "# 1:bin 2:count 3:bin center 4:bin width" << endl;
        //for (int ibin = 0; ibin < nbin+2; ++ibin) {
        //  *ss << ibin << " " << (int)count_global[ibin] << " " << range[0]+(range[1]-range[0])*(double(ibin)+0.5)/double(nbin) << " "
        //      << " " << (range[1]-range[0])/double(nbin) << endl;
        //}

        // cleanup...

        delete[] count_global;

      }

      // advance rank0 for the next one...

      ++rank0;
      if (rank0 == mpi_size)
        rank0 = 0;

      // and rest the averages if reset_after_write...

      if (reset_after_write)
        resetBuffers();

    }

    if (wvar.first == "")
      delete[] wvar_ptr;
  }

  void resetBuffers() {

    assert(count);
    for (int i = 0; i < nbin; ++i)
      count[i] = 0.0;

  }

  void flush() {

    // any files to be written are stored in pdfpFileList. Could be on any ranks...

    for (list<pair<string,stringstream*> >::iterator iter = pdfpFileList.begin(); iter != pdfpFileList.end(); ++iter) {
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
      }
      out_file.close();
      delete iter->second;
    }
    pdfpFileList.clear();

  }

};

#endif
