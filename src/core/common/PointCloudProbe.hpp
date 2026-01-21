#ifndef POINTCLOUD_PROBE_HPP
#define POINTCLOUD_PROBE_HPP 

enum PrecisionType {
  FLOAT_PRECISION,
  DOUBLE_PRECISION
};

class PointCloudProbe {
public:
  
  // this guy could have a hash or some indication of where it came from...
  
  int8 npt_global;
  int npt; // local 
  double (*x_pt)[3];
  int *icv_pt;
  int8 *ipt0_pt;
  string name; 
  int interval; 
  vector<pair<string,CtiRegister::CtiData*> > var_vec;
  int8 * ptora_striped;
  DistributedDataExchanger * dde_striped; // pull and push from the striped dist (the one for reading/writing)
  int *wtopt_i; // for interpolating from noocv to pt
  pair<int,double> *wtopt_v; 
  int precision; 
  bool b_ascii;

  PointCloudProbe() {
    npt_global = 0;
    ptora_striped = NULL;
    dde_striped = NULL;
    npt = 0;
    x_pt = NULL;
    icv_pt = NULL;
    ipt0_pt = NULL;
    name = "";
    interval = -1; 
    wtopt_i = NULL;
    wtopt_v = NULL;
    precision = DOUBLE_PRECISION;
    b_ascii = true;
  }

  ~PointCloudProbe() {
    DELETE(ptora_striped);
    if (dde_striped != NULL) delete dde_striped;
    DELETE(x_pt);
    DELETE(icv_pt);
    DELETE(ipt0_pt);
    var_vec.clear();
    DELETE(wtopt_i);
    DELETE(wtopt_v);
  }

  double atod(char * token) {
    double d;
    sscanf(token, "%lf", &d);
    return (d);
  }

  int getNextToken(char * token,int& pos,const char * buf, const int max_pos) {
    int token_pos = 0;
    while (pos < max_pos) {
      char c = buf[pos++];
      //cout << "getNextToken pos: " << pos << " c: \"" << c << "\"" << endl;
      if ( (c == ' ') || (c == '\t') || (c == '\n') || (c == ',') || (c == '\0') || c == 13 || c == 31 || c == EOF) {
        if (token_pos > 0) {
          // we have a valid token, so quit...
          break;
        } 
        else {
          // nothing yet - keep going...
          continue;
        }
      } 
      else { 
        token[token_pos++] = c;
      }
    }
    // null terminate token...
    if (token_pos != 0) token[token_pos] = '\0';
    return token_pos; // returns size
  }

  // restricted to no header and 3 columns of double data.
  // entries can be of any length though
  void readPointsAscii(vector<double>& x_vec,const string& filename) {
    
    COUT1("PointCloudProbe::readPointsAscii(\"" << filename << "\")...");
      
    // open...

    char dummy[128]; assert(filename.length() < 128);
    sprintf(dummy,"%s",filename.c_str());
    MPI_File fh;
    int ierr = MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
    if (ierr != 0) {
      CERR("cannot open restart file: " << filename);
    }

    // read in relevant chunk of file into "chunk",
    // which starts at location in the file globalstart
    // and has size mysize 

    // figure out who reads what
    MPI_Offset filesize;
    MPI_File_get_size(fh, &filesize);
    //filesize--;  /* get rid of text file eof */
    int8 mysize = filesize/mpi_size;
    MPI_Offset globalstart = mpi_rank * mysize;
    MPI_Offset globalend   = globalstart + mysize - 1;
    if (mpi_rank == mpi_size-1) globalend = filesize-1;
    mysize = globalend - globalstart + 1;

    // allocate memory 
    char *chunk = new char[mysize];

    // everyone reads in their part 
    readChunkedData<char>(fh, globalstart, chunk, mysize, false,mpi_comm);

    // count and index entries that you own

    int8 *ofonl = NULL;
    int8 nnl;
    for (int iter = 0; iter < 2; ++iter) {
      nnl = 0;
      for (int i = 0; i < mysize; ++i) {
        if (chunk[i] == '\n' || chunk[i] == EOF) {
          if (iter == 0) 
            nnl++;
          else
            ofonl[nnl++] = i+globalstart;
        }
      }
      if (iter == 0)
        ofonl = new int8[nnl];
    }
    delete[] chunk; chunk = NULL;

    // scan to determine which global entry/newlines you own...

    int8 nnl_global;
    MPI_Scan(&nnl,&nnl_global,1,MPI_INT8,MPI_SUM,mpi_comm);
    int8 my_inl_disp = nnl_global - nnl;
    MPI_Bcast(&nnl_global,1,MPI_INT8,mpi_size-1,mpi_comm); // everybody needs nl_cnt_global to set offset...

    // build striping
    int8 *nlora_striped = NULL;
    MiscUtils::calcThresholdDist(nlora_striped,nnl_global,mpi_size,DIST_THRESHOLD);



    // build dde 
    int8 *my_inl_global = new int8[nnl];
    const int nnl_striped = nlora_striped[mpi_rank+1]-nlora_striped[mpi_rank];
    int8 *ofonl_striped = new int8[nnl_striped+1]; ofonl_striped[0] = 0;

    int used_rank = 0;
    if (nnl_striped > 0) used_rank = 1;
    int nrank;
    MPI_Scan(&used_rank,&nrank,1,MPI_INT,MPI_SUM,mpi_comm);
    MPI_Bcast(&nrank,1,MPI_INT,mpi_size-1,mpi_comm); // everybody needs to know used striped ranks ...

    // send offsete of new-line to striped...
    {
      for (int i = 0; i < nnl; ++i) my_inl_global[i] = my_inl_disp+i; 
      DistributedDataExchanger dde(my_inl_global,nnl,nlora_striped);
      dde.push(ofonl_striped+1,ofonl);
    }
    delete[] my_inl_global;
    DELETE(ofonl);

    // send up your last new line...
    MPI_Request send_req,recv_req;
    if (mpi_rank < nrank-1) 
      MPI_Isend(&ofonl_striped[nnl_striped],1,MPI_INT8,mpi_rank+1,1234,MPI_COMM_WORLD,&send_req);
    if ((mpi_rank > 0) && (mpi_rank < nrank))
      MPI_Irecv(&ofonl_striped[0],1,MPI_INT8,mpi_rank-1,1234,MPI_COMM_WORLD,&recv_req);
    if (mpi_rank < nrank-1) 
      MPI_Wait(&send_req,MPI_STATUSES_IGNORE);
    if ((mpi_rank > 0) && (mpi_rank < nrank))
      MPI_Wait(&recv_req,MPI_STATUSES_IGNORE);
    if (mpi_rank == 0) ofonl_striped[0] = ofonl_striped[1];

    // now try reading again

    globalstart = ofonl_striped[0]+1;
    globalend = ofonl_striped[nlora_striped[mpi_rank+1]-nlora_striped[mpi_rank]];
    delete[] ofonl_striped;
    delete[] nlora_striped;
    //if (mpi_rank == mpi_size-1) globalend = filesize-1;
    mysize = globalend - globalstart + 1;
    if (mysize <= 0) {
      mysize = globalend = globalstart = 0;
    }
    assert(chunk == NULL); chunk = new char[mysize];

    // everyone reads in their part 
    readChunkedData<char>(fh, globalstart, chunk, mysize, false,mpi_comm);
    
    /*
    // for ranks > 0 send back the first line, so that the prev rank can complete their
    // full line. rank 0 has the one line header ("<npts> <nvars>\n"). we don't need it here, 
    // so just skip it...

    int locstart=0;
    while (chunk[locstart] != '\n') 
      locstart++;
    locstart++;
    
    // send/recv size...
    MPI_Request send_req,recv_req;
    if (mpi_rank > 0) 
      MPI_Isend(&locstart,1,MPI_INT,mpi_rank-1,1234,MPI_COMM_WORLD,&send_req);
    int nrecv = 0;
    if (mpi_rank < mpi_size-1)
      MPI_Irecv(&nrecv,1,MPI_INT,mpi_rank+1,1234,MPI_COMM_WORLD,&recv_req);
    if (mpi_rank > 0) 
      MPI_Wait(&send_req,MPI_STATUSES_IGNORE);
    if (mpi_rank < mpi_size-1) 
      MPI_Wait(&recv_req,MPI_STATUSES_IGNORE);
    
    // send/recv chars....
    MPI_Request send_req2,recv_req2;
    if (mpi_rank > 0) 
      MPI_Isend(chunk,locstart,MPI_CHAR,mpi_rank-1,2345,MPI_COMM_WORLD,&send_req2);
    char *overlap = new char[nrecv];
    if (mpi_rank < mpi_size-1) 
      MPI_Irecv(overlap,nrecv,MPI_CHAR,mpi_rank+1,2345,MPI_COMM_WORLD,&recv_req2);
    if (mpi_rank > 0) 
      MPI_Wait(&send_req2,MPI_STATUSES_IGNORE);
    if (mpi_rank < mpi_size-1) 
      MPI_Wait(&recv_req2,MPI_STATUSES_IGNORE);

    // build full char (for checking)...
    char* full_chunk = new char[mysize-locstart+nrecv]; 
    for (int i = locstart; i < mysize; ++i) full_chunk[i-locstart] = chunk[i];
    for (int i = 0; i < nrecv; ++i) full_chunk[i+mysize-locstart] = overlap[i];
    delete[] overlap;
    delete[] chunk;

    // get doubles from char tokens...
    char token[256];
    int pos = 0;
    assert(x_vec.size() == 0);
    while (pos < (nrecv+mysize-locstart)) { 
      if (getNextToken(token,pos,full_chunk,nrecv+mysize-locstart)) {
        x_vec.push_back(atod(token));
        //cout << mpi_rank << " " << x_vec.back() << endl;
      }
    }
    assert(x_vec.size()%3 == 0); 
    delete[] full_chunk;
    */

    // get doubles from char tokens...
    char token[256];
    int pos = 0;
    assert(x_vec.size() == 0);
    while (pos < mysize) { 
      if (getNextToken(token,pos,chunk,mysize)) {
        x_vec.push_back(atod(token));
        //if (mpi_rank == 0) cout << mpi_rank << " " << x_vec.back() << endl;
      }
    }
    assert(x_vec.size()%3 == 0); 
    delete[] chunk;

    // and close...
    MPI_File_close(&fh);

  }

  // to support boundary data io, we need to build a dde for any boundary data that
  // gets registered.
  void initDdeStuff() {
    assert(dde_striped == NULL);
    assert(ptora_striped == NULL);

    // get global info...
    int8 * ptora = NULL;
    MiscUtils::buildXora(ptora,npt);
    npt_global = ptora[mpi_size];
    int8 * my_ipt_global = new int8[npt];
    for (int ipt = 0; ipt < npt; ++ipt)
      my_ipt_global[ipt] = ipt + ptora[mpi_rank];
    delete[] ptora;

    // everyone build ptora based on our common global index...
    MiscUtils::calcThresholdDist(ptora_striped,npt_global,mpi_size,DIST_THRESHOLD);
    dde_striped = new DistributedDataExchanger(my_ipt_global,npt,ptora_striped);
    delete[] my_ipt_global;
  }

  // write out pbin without delta (double for consistency with other pbins)...
  void writePointsBinary() const {

    assert(x_pt != NULL);
    assert(ipt0_pt != NULL);
    assert(ptora_striped != NULL);
    assert(dde_striped != NULL);

    string filename = name+".pbin";

    COUT1("PointCloudProbe::writePointsBinary(): " << filename);
    
    // 1. pack the stided data into a contiguous buffer and
    // push to stiped data using the xora/dde...
      
    const int nrecv = ptora_striped[mpi_rank+1]-ptora_striped[mpi_rank];

    char dummy[128]; assert(filename.length() < 128);
    sprintf(dummy,"%s",filename.c_str());
    MPI_File_delete(dummy,MPI_INFO_NULL);

    MPI_File fh;
    MiscUtils::mkdir_for_file_collective(dummy); // allow the user to have specified a subdirectory...
    MPI_File_open(mpi_comm,dummy,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    // header...

    if (mpi_rank == 0) {
      int8 ibuf[4] = { POINTS_IO_MAGIC_NUMBER, POINTS_IO_VERSION, npt_global, 2 }; // last entry means include index
      MPI_File_write_at(fh,0,ibuf,4,MPI_INT8,MPI_STATUS_IGNORE);
    }
    int8 header_size = 4*int8_size;

    // x_pt...

    double (*recv_buf_dble)[3] = new double[nrecv][3];
    dde_striped->push(recv_buf_dble,x_pt); 
    writeChunkedData<double>(fh,header_size+ptora_striped[mpi_rank]*3*double_size,(double*)recv_buf_dble,3*nrecv,mpi_comm);
    delete[] recv_buf_dble;

    // ipt0_ipt...
    
    int8 (*recv_buf_int8) = new int8[nrecv];
    dde_striped->push(recv_buf_int8,ipt0_pt); 
    writeChunkedData<int8>(fh,header_size+npt_global*3*double_size+ptora_striped[mpi_rank]*int8_size,recv_buf_int8,nrecv,mpi_comm);
    delete[] recv_buf_int8;

    // and close...

    MPI_File_close(&fh);

  }

  // write out pbin without delta (double for consistency with other pbins)...
  void writePointsAscii() const {
    assert(b_ascii);

    assert(x_pt != NULL);
    assert(ipt0_pt != NULL);
    assert(ptora_striped != NULL);
    assert(dde_striped != NULL);

    string filename = name+".pxyz";

    COUT1("PointCloudProbe::writePointsAscii(): " << filename);
    
    // 1. pack the stided data into a contiguous buffer and
    // push to stiped data using the xora/dde...
      
    const int nrecv = ptora_striped[mpi_rank+1]-ptora_striped[mpi_rank];

    char dummy[128]; assert(filename.length() < 128);
    sprintf(dummy,"%s",filename.c_str());
    MPI_File_delete(dummy,MPI_INFO_NULL);

    MPI_File fh;
    MiscUtils::mkdir_for_file_collective(dummy); // allow the user to have specified a subdirectory...
    MPI_File_open(mpi_comm,dummy,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    // header...

    char tmp[4*23+1];
    if (mpi_rank == 0) {
      sprintf(tmp,"%22s %22s %22s %22s\n","index","x","y","z");
      MPI_File_write_at(fh,0,tmp,4*23,MPI_CHAR,MPI_STATUS_IGNORE);
    }
    int8 header_size = 4*23;

    // x_pt...

    double (*recv_buf_dble)[3] = new double[nrecv][3];
    dde_striped->push(recv_buf_dble,x_pt); 
    
    char* carray = new char[23*4*nrecv+1];
    for (int irecv = 0; irecv < nrecv; ++irecv) {
      sprintf(tmp,"% 18.15e ",recv_buf_dble[irecv][0]);  
      for (int i = 0; i < 23; ++i) carray[23*(4*irecv+1)+i] = tmp[i];
      sprintf(tmp,"% 18.15e ",recv_buf_dble[irecv][1]);  
      for (int i = 0; i < 23; ++i) carray[23*(4*irecv+2)+i] = tmp[i];
      sprintf(tmp,"% 18.15e\n",recv_buf_dble[irecv][2]);  
      for (int i = 0; i < 23; ++i) carray[23*(4*irecv+3)+i] = tmp[i];
    }
    delete[] recv_buf_dble;

    // ipt0_ipt...
    
    int8 (*recv_buf_int8) = new int8[nrecv];
    dde_striped->push(recv_buf_int8,ipt0_pt); 
    for (int irecv = 0; irecv < nrecv; ++irecv) {
      sprintf(tmp,"% 22lld ",recv_buf_int8[irecv]); // first column
      for (int i = 0; i < 23; ++i) carray[23*(4*irecv+0)+i] = tmp[i];
    }
    delete[] recv_buf_int8;

    writeChunkedData<char>(fh,header_size+ptora_striped[mpi_rank]*4*23,carray,4*23*nrecv,mpi_comm);
    delete[] carray;

    // and close...

    MPI_File_close(&fh);

  }

  void writeReadMe(Param * param, const string& infile) {

    const string filename = name+".README";
    char dummy[128]; assert(filename.length() < 128);
    sprintf(dummy,"%s",filename.c_str());
    MiscUtils::mkdir_for_file_collective(dummy); // allow the user to have specified a subdirectory...
    if (mpi_rank == 0) {
      ofstream outFile;
      outFile.open(filename.c_str(),ofstream::trunc);
      assert(outFile.is_open());
      outFile << setprecision(8);

      outFile << "# " << param->str() << endl;
      outFile << "# mles hash id " << RestartHashUtilities::mlesHash << endl;
      outFile << "# sles hash id " << RestartHashUtilities::slesHash << endl;
      outFile << endl;

      outFile << "input file: " << infile << endl; 
      outFile << "number of points (npt): " << npt_global << endl;
      outFile << "points io magic number for endianness check (num) : " << POINTS_IO_MAGIC_NUMBER << endl; 
      outFile << "points io version (ver): " << POINTS_IO_VERSION << endl; 
      outFile << "points io delta flag (flag): " << 0 << endl; 
      outFile << "point cloud data precision (prec): " << precision; 
      if (precision == FLOAT_PRECISION) 
        outFile << " => float" << endl;
      else 
        outFile << " => double" << endl;
      outFile << "number of scalars (nsc): " << var_vec.size() << endl;
      outFile << endl;

      outFile << "the scalars are:" << endl;
      for (int i = 0; i < var_vec.size(); ++i) {
        outFile << "   sc_" << i+1 << " " << var_vec[i].first << endl;
      }
      outFile << endl;

      outFile << "point cloud coordinates written in binary format (.pbin): " << name+".pbin" << endl;
      outFile << "    num[long] ver[long] npt[long] flg[long] " << endl;
      outFile << "    x_1[double] y_1[double] z_1[double] " << endl;
      outFile << "    x_2[double] y_2[double] z_2[double] " << endl;
      outFile << "    ..." << endl;
      outFile << "    x_npt[double] y_npt[double] z_npt[double] " << endl;
      outFile << "    i_1[long] " << endl;
      outFile << "    i_2[long] " << endl;
      outFile << "    ..." << endl;
      outFile << "    i_npt[long] " << endl;
      outFile << endl; 

      outFile << "point cloud data written in binary fomat: " << name + "XXXXXXXX.pcd" << ", where XXXXXXXX is the step number." << endl;
      outFile << "    num[long] ver[long] npt[long] nsc[long] prec[long]" << endl;
      if (precision == FLOAT_PRECISION) {
        outFile << "    sc_1(x_1,y_1,z_1)[float] sc_1(x_2,y_2,z_2)[float] ... sc_1(x_npt,y_npt,z_npt)[float] " << endl;
        outFile << "    sc_2(x_1,y_1,z_1)[float] sc_2(x_2,y_2,z_2)[float] ... sc_2(x_npt,y_npt,z_npt)[float] " << endl;
        outFile << "    ..." << endl;
        outFile << "    sc_nsc(x_1,y_1,z_1)[float] sc_nsc(x_2,y_2,z_2)[float] ... sc_nsc(x_npt,y_npt,z_npt)[float] " << endl;
      }
      else {
        outFile << "    sc_1(x_1,y_1,z_1)[double] sc_1(x_2,y_2,z_2)[double] ... sc_1(x_npt,y_npt,z_npt)[double] " << endl;
        outFile << "    sc_2(x_1,y_1,z_1)[double] sc_2(x_2,y_2,z_2)[double] ... sc_2(x_npt,y_npt,z_npt)[double] " << endl;
        outFile << "    ..." << endl;
        outFile << "    sc_nsc(x_1,y_1,z_1)[double] sc_nsc(x_2,y_2,z_2)[double] ... sc_nsc(x_npt,y_npt,z_npt)[double] " << endl;
      }
      outFile << endl; 

      outFile.close();
    }
  }

};
#endif
