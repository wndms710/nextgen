
#include "CTI.hpp"
#include "MiscUtils.hpp"
#include "Dmd.hpp"
#include "ByteSwap.hpp"
#include "MpiTimer.hpp"
#include <set>

using namespace CTI;
using namespace MiscUtils;

class FileHash {
public:
  FILE* fp;
  int8 offset;
  double * data;

  FileHash() : fp(NULL), offset(0), data(NULL) {}
};

void my_fseek(FILE * fp, int8 offset) {

  // not very efficient, but works on BG
  int8 cur_offset = offset;
  int8 max_offset = 1000000000ll;

  fseek(fp,0,SEEK_SET);
  while (cur_offset > max_offset) {
    fseek(fp,max_offset,SEEK_CUR);
    cur_offset -= max_offset;
  }
  fseek(fp,cur_offset,SEEK_CUR);
}

class ModeDescriptor {
public:
  double freq;
  double damp;
  double amp_n;  // amplitude at the final instant..
  double amp_0;  // amplitude at the first instant
  int idx;

  ModeDescriptor(const double _li, const double _lr, const double _an, const double _a0, const int _i) {
    freq  = _li;
    damp  = _lr;
    amp_n = _an;
    amp_0 = _a0;
    idx   = _i;
  }

  bool operator<(const ModeDescriptor& _d) const {
    return (amp_n < _d.amp_n);
  }
};

class Modal {
public:

  // snapshot correlation matrix setup

  double * A_corr;    // temporal correlation matrix..
  int ns;             // number of snapshots
  string name_prefix; // snapshot name prefix
  string var_name;    // variable name..
  int istart,iend;    // start and ending steps for the snapshot..
  int istep;          // step index..
  int* snora;         // snapshot of rank
  vector<FileHash> fh; // snapshot file pointers..
  int8 ncv_global;
  double var_ref;     // ability to remove a single global reference value..
  double var_norm;    // normalize the variable with a single global value...
  double dt_samp;     // time step between snapshots (assumed constant in the analysis)

  // dmd decomposition inputs and output

  int nnz;
  double * lambda_i;
  double * lambda_r;
  complex<double>* Vr;
  double * Wk;
  double * amp;

  Modal() {
    A_corr      = NULL;
    ns          = 0;
    name_prefix = "";
    istart      = -1;
    iend        = -1;
    istep       = -1;
    snora       = NULL;
    ncv_global  = -1;
    var_ref     = 0.0;
    var_norm    = 0.0;

    dt_samp     = 0.0;
    nnz         = 0;
    lambda_i    = NULL;
    lambda_r    = NULL;
    Vr          = NULL;
    Wk          = NULL;
    amp         = NULL;

  }

  void init() {

    Param * param = getParam("SNAPSHOT_MODES");
    int iarg = 0;

    // default reference values assume no normalization..

    var_ref  = 0.0;
    var_norm = 1.0;

    while ( iarg < param->size()) {
      const string token = param->getString(iarg++);

      if ( token == "NAME") {
        name_prefix = param->getString(iarg++);
      } else if ( token == "START") {
        istart      = param->getInt(iarg++);
      } else if ( token == "END") {
        iend        = param->getInt(iarg++);
      } else if ( token == "STEP") {
        istep       = param->getInt(iarg++);
      } else if ( token == "VAR") {
        var_name    = param->getString(iarg++);
      } else if ( token == "REMOVE_REF") {
        var_ref     = param->getDouble(iarg++);
      } else if ( token == "NORMALIZE_VAR") {
        var_norm    = param->getDouble(iarg++);
      } else {
        CERR(" > unrecognized token in SNAPSHOT_MODES " << token);
      }
    }

    ns = (iend - istart)/istep;
    assert( A_corr == NULL);
    A_corr = new double[ns*ns];
    for (int ii = 0; ii < ns*ns; ++ii)
      A_corr[ii] = 0.0;

    //MiscUtils::buildUniformXora(snora,ns);
    MiscUtils::calcThresholdDist(snora,ns,mpi_size,1);
    const int nfl = snora[mpi_rank+1] - snora[mpi_rank];

    // tell everyone to go ahead and open the files that they are responsible for..

    MiscUtils::dumpRange(&nfl,1,"nfl");
    if ( mpi_rank == 0) 
      cout << " > number of snapshots : " << ns << endl;
    

    fh.resize(nfl);

    // compute the dt_samp from the first two snapshots.. the first snapshot
    // is guaranteed to live on rank 0, but the next snapshot may live elsewhere..

    double time01[2] = {0.0,0.0};
    for (int ifl = 0; ifl < nfl; ++ifl) {

      char filename[128];
      const int step = (snora[mpi_rank]+ifl)*istep + istart;
      sprintf(filename,"%s.%08d.sles",name_prefix.c_str(),step);
      fh[ifl].fp = fopen(filename,"rb");

      // we need to advance to the start of the record..
      fh[ifl].offset = sizeof(int)*2;
      Header header;

      int done   = 0;
      ncv_global = -1;

      bool found_time  = false;
      bool found_var   = false;

      while ( done != 1) {

        my_fseek(fh[ifl].fp,fh[ifl].offset);
        fread(&header,sizeof(Header),1,fh[ifl].fp);
        
        if ( mpi_rank == 0)
          cout << " found record: " << header.name << endl;

        switch (header.id) {

        case UGP_IO_D0:
          {
            if ( strcmp(header.name,"time") == 0) {

              found_time = true;

              if ( snora[mpi_rank]+ifl == 0) {
                time01[0] = header.rdata[0];
                cout << " time t0 found: " << header.rdata[0] << endl;
                cout.flush();
              } else if ( snora[mpi_rank]+ifl == 1) {
                time01[1] = header.rdata[0];
                cout << " time t1 found: " << header.rdata[1] << endl;
                cout.flush();
              }
            }
          }
          break;

        case UGP_IO_CV_D1:
          {
            if ( strcmp(header.name,var_name.c_str()) == 0) {

              found_var = true;

              if ( ncv_global == -1) {
                ncv_global = header.idata[0];
              } else {
                assert( ncv_global == header.idata[0]);
              }

              // advance the offset to the start of the data record..
              fh[ifl].offset += header_size;
              my_fseek(fh[ifl].fp,fh[ifl].offset);
              done = 1;
            }
          }
          break;

        case UGP_IO_EOF:
          {
            cout << "mpi_rank : " << mpi_rank << " ; unable to find " << var_name << endl;
            done = 1;
            throw(1);
          }
        }

        if ( done == 0 ) {
          // if we found the record, we actually want to store the offset for fast rewind
          fh[ifl].offset += header.skip;
        }
      }
    }//for(ifl)

    // everyone needs ncv_global, it is possible that some ranks actually do not
    // have any snapshots to read, but they will still help in the correlation calc.

    {
      int8 ncv_global_max;
      MPI_Allreduce(&ncv_global,&ncv_global_max,1,MPI_INT8,MPI_MAX,mpi_comm);

      if ( ncv_global == -1) {

        assert (fh.size() == 0);
        ncv_global = ncv_global_max;

      } else {
        assert( ncv_global_max == ncv_global);
      }
    }

    // compute the dt_samp from the first two entries..

    {
      double tmp[2];
      MPI_Allreduce(time01,tmp,2,MPI_DOUBLE,MPI_SUM,mpi_comm);
      dt_samp = tmp[1] - tmp[0];

      if ( mpi_rank == 0 )
        cout << " dt_samp (estimated) : " << dt_samp << endl;
    }

  }

  void buildCorrelationMatrix() {

    if ( mpi_rank == 0)
      cout << " buildCorrelationMatrix() " << endl;

    MPI_Barrier(mpi_comm);
    double wtime0 = MPI_Wtime();

    double * A_ = new double[ns*ns];
    for (int ii = 0; ii < ns*ns; ++ii)
      A_[ii] = 0.0;

    const int local_mem   = getIntParam("LOCAL_MEM", 1024*1024*8); // 8MB record...
    const int nfl         = fh.size();
    const int8 nq          = int8(local_mem)*int8(mpi_size)/int8(8*ns); // contiguous cvs read at once..

    if ( mpi_rank == 0 ) 
      cout << " > nq calc: " << local_mem << "   " << mpi_size << "   " << ns << "   " << local_mem*mpi_size << "    " << nq << endl;

    assert(nq > 0);

    for (int ifl = 0; ifl < nfl; ++ifl) {
      assert( fh[ifl].data == NULL);
      fh[ifl].data = new double[nq];

      //setvbuf( fh[ifl].fp, (char *)NULL, _IONBF, 0 ); // turn off the buffering
    }

    double io_time     = 0.0;
    double build_time  = 0.0;
    double comm_time   = 0.0;
    double comp_time   = 0.0;

    MpiTimer timer;

    int8 icv_global = 0;
    while ( icv_global < ncv_global) {

     const int nn  = (int) min(ncv_global-icv_global,int8(nq));
     
     if ( mpi_rank == 0 )
        cout << " > starting icv " << icv_global << " of  " << ncv_global << "    " << nq << "    " << nn << endl;

     timer.start("init");

      for (int ifl = 0; ifl < nfl; ++ifl) {

        double io_elapsed = MPI_Wtime();
        fread(fh[ifl].data,sizeof(double),nn,fh[ifl].fp);
        io_elapsed = MPI_Wtime() - io_elapsed;
        io_time   += io_elapsed;

        // remove the specified reference and renormalize the variable.

        const double tmp = 1.0/var_norm;
        for (int i =0 ; i <nn; ++i) {
          fh[ifl].data[i] = (fh[ifl].data[i] - var_ref)*tmp;
        }

        // HACK .. fill the data buffer with known values
        /*
        for (int ii = 0; ii < nn; ++ii)
          fh[ifl].data[ii] = double(snora[mpi_rank] + ifl); // snapshot index..
        */
      }

      MPI_Barrier(mpi_comm);
      timer.split("io");

      
      double build_elapsed = MPI_Wtime();
      
      // at the end of the transposition each rank will
      // hold a chunk of (nchunk) cvs for all time..
      // for the global icv range you need to offset with the icv_global


      int * cvora = NULL; MiscUtils::buildUniformXora(cvora,nn);
      const int nchunk = cvora[mpi_rank+1] - cvora[mpi_rank];

      // the amount of data that has been read by this rank is nn* nfl
      // after the transposition you should have ns*nchunk

      double * send_buf = new double[nn * nfl];
      int * send_count  = new int[mpi_size];
      int * send_disp   = new int[mpi_size];

      // first pack the send_buf; the data it sends to a rank
      // needs to be contiguous wrt to the time coordinate for a
      // given cv.

      for (int rank = 0; rank < mpi_size; ++rank)
        send_count[rank] = nfl*(cvora[rank+1] - cvora[rank]);

      int jj = 0;
      for (int ii = 0; ii < nn; ++ii) {
        for (int ifl = 0; ifl < nfl; ++ifl) {
          send_buf[jj++] = fh[ifl].data[ii];
        }
      }

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

      timer.split("pack"); 

      // set up the receive side of the data ..

      double * recv_buf = new double[ns * nchunk];
      int * recv_count  = new int[mpi_size];
      int * recv_disp   = new int[mpi_size];

      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_disp[rank-1] + recv_count[rank-1];

      // check that the counts are exactly what we believe they should
      // be and then perform the exchange.

      const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
      assert( recv_count_sum == ns*nchunk);

      MPI_Alltoallv(send_buf,send_count,send_disp,MPI_DOUBLE,
                    recv_buf,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);


      timer.split(" communication"); 
      // now need one more reorder of the recv_buf to get going..
      // recall that the time is contiguous across the ranks..

      double * recv_buf2 = new double[ns*nchunk];
      jj = 0;
      for (int ii = 0; ii < nchunk; ++ii) {
        for (int rank = 0; rank < mpi_size; ++rank) {
          const int ns_rank = snora[rank+1] - snora[rank];
          const int i_offset = recv_disp[rank]+ii*ns_rank;
          for (int is = 0; is < ns_rank; ++is) {
            recv_buf2[jj++] = recv_buf[recv_disp[rank]+ii*ns_rank+is];
            //recv_buf2[jj++] = recv_buf[i_offset+is];
          }
        }
      }

      timer.split("unpack");

      const double tmp_time = MPI_Wtime();
      comm_time += tmp_time - build_elapsed;

      // XXX HACK ... when the buffer was set with known values, lets
      // ensure that we got exactly what we were expecting..
      /*
      for (int k = 0; k < nchunk; ++k) {
        for (int i =0; i < ns; ++i) {
        assert( recv_buf2[ns*k+i] == double(i));
        }
      }
      */

      // reduce the local data into the temporary correlation matrix
      // note we are only going to be build the lower triangle since
      // the matrix is symmetric

      /*
      for (int i =0; i < ns; ++i) {
        for (int j = 0; j <= i; ++j) {
          double sum = 0.0;
          for (int k = 0; k < nchunk; ++k)
            sum += recv_buf[ns*k+i]*recv_buf[ns*k+j];
        }
      }
      */

      // reverse the order of the loops to exploit data contiguity..

      const double tmp  = 1.0/double(ncv_global);
      for (int k =0; k < nchunk; ++k) {
        for (int i =0; i < ns; ++i) {
          for (int j = 0; j <= i; ++j) {
            A_[(j)*ns+i] += recv_buf2[ns*k+i]*recv_buf2[ns*k+j]*tmp;
          }
        }
      }

      timer.split("compute");

      const double tmp_time2 = MPI_Wtime();
      comp_time += tmp_time2 - tmp_time;

      build_elapsed = tmp_time2 - build_elapsed;
      build_time   += build_elapsed;

      delete[] send_buf;
      delete[] send_count;
      delete[] send_disp;
      delete[] recv_buf;
      delete[] recv_buf2;
      delete[] recv_count;
      delete[] recv_disp;
      delete[] cvora;

      icv_global += int8(nn);
      timer.accumulate();
    }

    // reduce the correlation matrix across all of the nodes.  this
    // could be an MPI_Reduce, but we'll let everyone keep it anyway

    MPI_Allreduce(A_,A_corr,ns*ns,MPI_DOUBLE,MPI_SUM,mpi_comm);
    delete[] A_;

    MPI_Barrier(mpi_comm);
    double wtime = MPI_Wtime();

    if ( mpi_rank == 0) { 
      cout << " > correlation matrix built in [secs] : " << wtime - wtime0 << endl;
    }

    timer.report();

    MiscUtils::dumpRange(&io_time,1,"io_time");
    MiscUtils::dumpRange(&build_time,1,"build_time");
    MiscUtils::dumpRange(&comm_time,1,"comm_time");
    MiscUtils::dumpRange(&comp_time,1,"comp_time");

    // re-populate the upper triangle of the correlation matrix

    for (int i =0; i < ns; ++i) {
      for (int j = i+1; j < ns; ++j) {
        A_corr[(j)*ns+i] = A_corr[(i)*ns+j];
      }
    }

    // since the correlation matrix was so cumbersome to build, write
    // it out to the disk for safe-keeping

    if ( mpi_rank == 0) {

      FILE* fp = fopen("A.dat", "wb");
      fwrite(&ns,sizeof(int),1,fp);
      fwrite(A_corr,sizeof(double),ns*ns,fp);
      fclose(fp);

    }

    // cleanup the data buffers

    for (int ifl = 0; ifl < nfl; ++ifl) {
      delete[] fh[ifl].data;
      fh[ifl].data = NULL;
    }

  }

  void readCorrelationMatrix(Param* param) {

    FILE * fp = NULL;
    if ( mpi_rank == 0 ) {

      const string fname = param->getString(0);
      fp                 = fopen(fname.c_str(),"rb");

      if (!fp) {
        cerr << " > unable to open file: " << fname << endl;
        throw(1);
      }

      fread(&ns,sizeof(int),1,fp);

    }

    MPI_Bcast(&ns,1,MPI_INT,0,mpi_comm);

    assert( A_corr != NULL);

    if ( mpi_rank == 0) {

      fread(A_corr,sizeof(double),ns*ns,fp);
      fclose(fp);

    }

    MPI_Bcast(A_corr,ns*ns,MPI_DOUBLE,0,mpi_comm);

  }

  void dumpDmd() const {

    FILE * fp = fopen("freq_amp.dat", "w");
    set<ModeDescriptor> descriptors;

    for (int i =0; i < nnz; ++i) {

      if ( lambda_i[i] >= 0.0) {

        complex<double> lambda_c(lambda_r[i]*dt_samp,lambda_i[i]*2.0*M_PI*dt_samp);
        complex<double> mu_c     = exp(lambda_c);
        complex<double> mu_n     = pow(mu_c,ns);
        const double a_n         = amp[i]*std::abs(mu_n);
        descriptors.insert( ModeDescriptor(lambda_i[i],lambda_r[i], a_n, amp[i], i));
      }
    }

    for (set<ModeDescriptor>::reverse_iterator it = descriptors.rbegin(); it != descriptors.rend(); ++it) {
      fprintf(fp, "%12.8g    %12.8g    %12.8g    %12.8g    %d\n",
              it->freq, it->amp_0, it->damp, it->amp_n, it->idx);
    }

    fclose(fp);
  }

  void computeDmdFreqAndAmp() {

    // can overwrite the estimate sampling time if desired..

    if ( checkParam("DT_SAMP")) {
      dt_samp = getDoubleParam("DT_SAMP");
    }

    if ( mpi_rank == 0 ) {

      assert( amp      == NULL);
      assert( Wk       == NULL);
      assert( lambda_r == NULL);
      assert( lambda_i == NULL);
      assert( Vr       == NULL);
      assert( A_corr   != NULL);
      assert( dt_samp > 0.0);

      nnz = dmdAnalysis(ns, A_corr, lambda_r, lambda_i, Wk, Vr, amp, dt_samp);

      dumpDmd();
    }

    MPI_Bcast(&nnz, 1, MPI_INT, 0, mpi_comm);

    if ( mpi_rank != 0 ) {

      Wk       = new double[(ns-1)*(ns-1)];
      Vr       = new complex<double>[nnz*nnz];
      lambda_i = new double[nnz];
      lambda_r = new double[nnz];

    }

    MPI_Bcast(Wk, (ns-1)*(ns-1), MPI_DOUBLE, 0, mpi_comm);
    MPI_Bcast(Vr, nnz*nnz      , MPI_DOUBLE_COMPLEX, 0, mpi_comm);
    MPI_Bcast(lambda_i, nnz    , MPI_DOUBLE, 0, mpi_comm);
    MPI_Bcast(lambda_r, nnz    , MPI_DOUBLE, 0, mpi_comm);

  }

  void extractDmdModes(Param* param) {

    assert( param->getString(0) == "FREQUENCY");
    const double freq_nom = param->getDouble(1);

    // find the nearest mode to this frequency.

    int imode   = -1;

    if ( mpi_rank == 0 ) {

      double dist = 1e+16;
      for (int imo = 0; imo < nnz; ++imo) {

        const double d = abs(lambda_i[imo] - freq_nom);
        if ( d < dist) {
          dist  = d;
          imode = imo;
        }
      }

      cout << "  > Working on mode with frequency: " << lambda_i[imode] << endl;
    }

    MPI_Bcast(&imode,1,MPI_INT,0,mpi_comm);

    // we are going to write the mode in chunks so we need to
    // set up the files for the real and imaginary components

    char fname[128];
    sprintf(fname, "mode.%.3f.sles", lambda_i[imode]);

    MPI_File file;
    MPI_File_delete(fname,MPI_INFO_NULL);
    MPI_File_open(mpi_comm,fname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&file);

    if ( mpi_rank == 0 )  {
      int itmp[2] = {UGP_IO_MAGIC_NUMBER+1, 5};
      MPI_File_write(file,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
    }

    MPI_Offset offset = int_size * 2;

    // determine the chunk size of the cvs..

    const int local_mem   = getIntParam("LOCAL_MEM", 1024*1024*8); // 8MB record...
    const int nfl         = fh.size();
    const int nchunk      = local_mem/8;
    double * my_buf       = new double[nchunk];
    double * buf          = new double[nchunk];

    // resize the data buffer accordingly
    for (int ifl = 0; ifl < nfl; ++ifl) {

      if ( fh[ifl].data != NULL) {
        delete[] fh[ifl].data;
      }

      fh[ifl].data = new double[nchunk];
    }

    for (int iter = 0; iter < 2; ++iter) {

      if ( mpi_rank == 0 ) {
        if ( iter == 0 )
          cout << " >>  working on real part of the mode " << endl;
        else
          cout << " >>  working on imag part of the mode " << endl;
      }

      // rewind the file pointers to the beginning of the data record...

      for (int ifl = 0; ifl < nfl; ++ifl) {
        my_fseek(fh[ifl].fp,fh[ifl].offset);
      }

      if ( mpi_rank == 0 ) {
        Header header;

        if ( iter == 0 ) {
          sprintf(header.name,"mode-real");
        } else {
          sprintf(header.name,"mode-imag");
        }

        header.id = UGP_IO_CV_D1; // for now, this is all that's supported...
	header.skip = header_size + ncv_global * double_size;
	ByteSwap::setLswMswPairForInt8(header.idata+0,ncv_global);
	MPI_File_write_at(file,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      offset += header_size;

      int8 icv_global = 0;
      while ( icv_global < ncv_global) {

        if ( mpi_rank == 0 )
          cout << " > working on icv " << icv_global << " of " << ncv_global << endl;

        const int nn = min(int(ncv_global - icv_global), nchunk);

        // read the chunks of data for the files we need..

        for (int ii = 0; ii < nn; ++ii)
          my_buf[ii] = 0.0;

        for (int ifl = 0; ifl < nfl; ++ifl) {

          // the range of the reconstruct only extends to ns-1, so skip if you are the last snapshot
          if ( snora[mpi_rank] + ifl == ns-1)
            continue;

          fread(fh[ifl].data,sizeof(double),nn,fh[ifl].fp);
          complex<double> c_coeff;

          reconstructCoeff(c_coeff, Wk, Vr, nnz, ns, imode, snora[mpi_rank]+ifl);
          const double coeff = (iter == 0) ? c_coeff.real() : c_coeff.imag();
          for (int ii = 0; ii < nn; ++ii)
            my_buf[ii] += coeff * fh[ifl].data[ii];

        }

        MPI_Reduce(my_buf,buf,nn,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

        // mpi rank 0 writes this chunk bit

        if ( mpi_rank == 0 ) {
          MPI_File_write_at(file,offset,buf,nn,MPI_DOUBLE,MPI_STATUS_IGNORE);
        }

        offset += int8(nn)* double_size;
        icv_global += int8(nn);

      }

    }

    if ( mpi_rank == 0 ) {
      cout << " > EOF" << endl;
      Header header;
      header.id = UGP_IO_EOF;
      sprintf(header.name,"EOF");
      header.skip = header_size;
      MPI_File_write_at(file,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
    }

    offset += header_size;

    MPI_File_set_size(file,offset);
    MPI_File_close(&file);

    delete[] my_buf;
    delete[] buf;

    // cleanup the data buffers again

    for (int ifl = 0; ifl < nfl; ++ifl) {
      delete[] fh[ifl].data;
      fh[ifl].data = NULL;
    }

  }

  void run() {

    if ( Param * param = getParam("CORRELATION_MATRIX" )) {

      // allow all the ranks to get the correlation matrix
      readCorrelationMatrix(param);

    } else {

      // otherwise the correlation matrix needs to be built.
      buildCorrelationMatrix();

    }

    // to start, we are going to compute the dmd; we can
    // always enable the calculation of other modal decompositions
    // later.

    computeDmdFreqAndAmp();

    FOR_PARAM_MATCHING("WRITE_DMD_MODE") {
      extractDmdModes(&(*param));
    }

  }

  ~Modal() {
    DELETE(A_corr);
    DELETE(snora);
    DELETE(Wk);
    DELETE(Vr);
    DELETE(lambda_r);
    DELETE(lambda_i);
    DELETE(amp);
  }
};

int main(const int argc, char* argv[]) {

  try {

    CTI_Init(argc,argv,"modal.in");

    {
      Modal m;
      m.init();
      m.run();

    }

    CTI_Finalize();

  } catch(...) {

    CTI_Abort();

  }

  return 0;
}
