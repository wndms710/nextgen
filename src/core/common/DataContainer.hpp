#ifndef DATACONTAINER_HPP
#define DATACONTAINER_HPP

#include "RestartHashUtilities.hpp"
#include "CtiRegister.hpp"

class InData {
public:
  string name;
  int * data;
  CtiRegister::CtiData * ctiData; // link to corresponding registered ctiData
  InData(const string& name_) {
    name = name_;
    data = NULL;
    ctiData = NULL;
  }
  ~InData() {
    DELETE(data);
  }
};

class LnData {
public:
  string name;
  int8 * data;
  CtiRegister::CtiData * ctiData; // link to corresponding registered ctiData
  LnData(const string& name_) {
    name = name_;
    data = NULL;
    ctiData = NULL;
  }
  ~LnData() {
    DELETE(data);
  }
};

class DnData {
public:
  string name;
  double * data;
  CtiRegister::CtiData * ctiData; // link to corresponding registered ctiData
  DnData(const string& name_) {
    name = name_;
    data = NULL;
    ctiData = NULL;
  }
  ~DnData() {
    DELETE(data);
  }
};

class Dn3Data {
public:
  string name;
  double (*data)[3];
  CtiRegister::CtiData * ctiData; // link to corresponding registered ctiData
  Dn3Data(const string& name_) {
    name = name_;
    data = NULL;
    ctiData = NULL;
  }
  ~Dn3Data() {
    DELETE(data);
  }
};

// down here because it needs to know about DnData/Dn3Data
#include "CasDatReader.hpp"

class DataContainer {
public:

  // this guy could have a hash or some indication of where it came from...

  int8 ncv_global;
  int8 * cvora;
  int ncv; // local
  double (*x_cv)[3];
  list<pair<string,int> > iList; // int values
  list<pair<string,double> > dList; // double values
  list<DnData> dnList; // scalar fields
  list<Dn3Data> dn3List; // vector fields
  double * delta; // required for writing pbins

  DataContainer() {
    ncv_global = 0;
    cvora = NULL;
    ncv = 0;
    x_cv = NULL;
    delta = NULL;
  }

  ~DataContainer() {
    DELETE(cvora);
    DELETE(x_cv);
    DELETE(delta);
  }

  void rotate(const int iaxis,const double degrees,const bool b_rot_coord = true) {

    // rotate data around iaxis

    const double rad = degrees*M_PI/180.0;
    const double cos_t = cos(rad);
    const double sin_t = sin(rad);
    const int i1 = (iaxis+1)%3;
    const int i2 = (iaxis+2)%3;

    // rotate coordinate data...

    if (b_rot_coord) {
      assert(x_cv);
      for (int icv = 0; icv < ncv; ++icv) {
        const double x1 = x_cv[icv][i1]*cos_t - x_cv[icv][i2]*sin_t;
        const double x2 = x_cv[icv][i2]*cos_t + x_cv[icv][i1]*sin_t;
        x_cv[icv][i1] = x1;
        x_cv[icv][i2] = x2;
      }
    }

    // also rotate all vector quantities...

    for (list<Dn3Data>::iterator iter = dn3List.begin(); iter != dn3List.end(); ++iter) {
      for (int icv = 0; icv < ncv; ++icv) {
	const double v1 = iter->data[icv][i1]*cos_t - iter->data[icv][i2]*sin_t;
	const double v2 = iter->data[icv][i2]*cos_t + iter->data[icv][i1]*sin_t;
	iter->data[icv][i1] = v1;
	iter->data[icv][i2] = v2;
      }
    }

  }

  void rotate(const double point[3],const double _axis[3],const double degrees,const bool b_rot_coord = true) {
    double axis[3] = {_axis[0],_axis[1],_axis[2]};
    NORMALIZE(axis);
    const double angle = degrees/180.0*M_PI;  // radians

    const double cp_mat[3][3] = {
      {0.0,-axis[2],axis[1]},
      {axis[2],0.0,-axis[0]},
      {-axis[1],axis[0],0.0}
    };

    // create rotation matrix
    // R = cos0*I + sin0*[axis]_x + (1-cos0)*(u tensor-product u)
    double R[3][3];
    FOR_I3 {
      R[i][i] = cos(angle);
      FOR_J3 {
        if (i != j) R[i][j] = sin(angle)*cp_mat[i][j];
        R[i][j] += (1-cos(angle))*axis[i]*axis[j];
      }
    }

    // rotate coordinate data...

    if (b_rot_coord) {
      assert(x_cv);
      FOR_ICV {
        const double x[3] = DIFF(x_cv[icv],point);
        const double temp[3] = MATRIX_PRODUCT(R,x);
        FOR_I3 x_cv[icv][i] = temp[i] + point[i];
      }
    }

    // also rotate all vector quantities...

    for (list<Dn3Data>::iterator iter = dn3List.begin(); iter != dn3List.end(); ++iter) {
      FOR_ICV {
        const double temp[3] = MATRIX_PRODUCT(R,iter->data[icv]);
        FOR_I3 iter->data[icv][i] = temp[i];
      }
    }
  }

  void translate(const double dx, const double dy, const double dz) {

    // translate coordinate data...

    assert(x_cv);
    for (int icv = 0; icv < ncv; ++icv) {
      x_cv[icv][0] += dx;
      x_cv[icv][1] += dy;
      x_cv[icv][2] += dz;
    }

  }

  void mirror(const double point[3],const double _axis[3]) {
    double axis[3] = {_axis[0],_axis[1],_axis[2]};
    NORMALIZE(axis);

    // mirror coordinate data...

    if (true) {  // in case want to only do this to vectors and not positions...
      assert(x_cv);
      FOR_ICV {
        double x_plane[3];
        int ierr = MiscUtils::computeLinePlaneIntersection(x_plane,x_cv[icv],axis,point,axis);
        if (ierr == 0) {
          FOR_J3 x_cv[icv][j] = x_plane[j] + (x_plane[j]-x_cv[icv][j]);
        }
      }
    }

    // also mirror vector quantities...
    for (list<Dn3Data>::iterator iter = dn3List.begin(); iter != dn3List.end(); ++iter) {
      FOR_ICV {
        double in_plane[3];  // component of vector orthogonal to plane
        double orth_plane[3];  // component of vector normal to plane (gets inverted)
        FOR_I3 {
          orth_plane[i] = iter->data[icv][i]*axis[i];
          in_plane[i] = iter->data[icv][i] - orth_plane[i];

          iter->data[icv][i] = in_plane[i] - orth_plane[i];
        }
      }
    }
  }

  void buildOrCheckCvora(const int8 ncv_global_check) {

    if (cvora == NULL) {
      ncv_global = ncv_global_check;
      if (mpi_rank == 0) cout << " > ncv_global: " << ncv_global << endl;
      MiscUtils::calcUniformDist(cvora,ncv_global,mpi_size);
      assert(cvora[mpi_rank+1]-cvora[mpi_rank] < TWO_BILLION);
      ncv = cvora[mpi_rank+1]-cvora[mpi_rank];
    }
    else {
      assert(ncv_global == ncv_global_check);
      assert(ncv == cvora[mpi_rank+1]-cvora[mpi_rank]);
    }

  }

  void initPointsFromFluent(const string& cas_filename, int byte_swap = 0) {

    COUT1("DataContainer::initPointsFromFluent(\"" << cas_filename << "\")...");

    // 1.0 build x_cv/delta from cas...
    CasDatReader cdr;
    cdr.setForceByteSwap(byte_swap);
    cdr.initFromCas(x_cv,delta,cvora,ncv,ncv_global,cas_filename);
  }

  void writePointsFromFluent(const string& cas_filename) {

    COUT1("DataContainer::writePointsFromFluent(\"" << cas_filename << "\")...");

    // 1.0 build x_cv/delta from cas...
    CasDatReader cdr;
    cdr.initFromCas(x_cv,delta,cvora,ncv,ncv_global,cas_filename);

    size_t last_dot = cas_filename.find_last_of('.');
    string filename = cas_filename.substr(0,last_dot) + ".pbin";
    writePointsBinary(filename);

  }

  int initFromFluent(const string& cas_filename,const string& dat_filename,const string& csVarNames) {

    // csVarNames == comma-separated var names...

    // turn csVarNames into a set of strings...

    set<string> varNameSet;
    MiscUtils::splitCsv(varNameSet,csVarNames);

    return initFromFluent(cas_filename,dat_filename,varNameSet);

  }

  int initFromFluent(const string& cas_filename,const string& dat_filename,const set<string>& varNameSet) {

    COUT1("DataContainer::initFromFluent(\"" << cas_filename << "\",\"" << dat_filename << "\")...");

    // 1.0 build x_cv/delta from cas...
    CasDatReader cdr;
    cdr.initFromCas(x_cv,delta,cvora,ncv,ncv_global,cas_filename);

    // 2.0 get registered cv data from dat...
    cdr.initFromDat(dnList,dn3List,dat_filename,varNameSet,cvora,ncv,ncv_global);
    return 0;

  }

  int initFromPbinFluent(const string& pbin_filename,const string& dat_filename,const string& csVarNames) {

    // csVarNames == comma-separated var names...

    // turn csVarNames into a set of strings...

    set<string> varNameSet;
    MiscUtils::splitCsv(varNameSet,csVarNames);

    return initFromPbinFluent(pbin_filename,dat_filename,varNameSet);

  }

  int initFromPbinFluent(const string& pbin_filename,const string& dat_filename,const set<string>& varNameSet) {

    COUT1("DataContainer::initFromPbinFluent(\"" << pbin_filename << "\",\"" << dat_filename << "\")...");
    COUT1(" > mixed pbin/dat importing is no longer supported; binary cas/dat has been included to speed up importing directly from Fluent files");
    // 1.0 extract x_cv (and delta from pbin)
    // readPointsBinary(pbin_filename);
    //
    // // 2.0 get registered cv data from dat...
    // CasDatReader cdr;
    // cdr.initFromDat(dnList,dn3List,dat_filename,varNameSet,cvora,ncv,ncv_global);
    return 0;

  }

  int initFromFluentIp(const string& filename,const string& csVarNames) {

    // csVarNames == comma-separated var names...

    // turn csVarNames into a set of strings...

    set<string> varNameSet;
    MiscUtils::splitCsv(varNameSet,csVarNames);

    return initFromFluentIp(filename,varNameSet);

  }

  int initFromFluentIp(const string& filename,const set<string>& varNameSet) {

    COUT1("DataContainer::initFromFluentIp(\"" << filename << "\")...");

    // read header info on rank 0
    int version = -1;
    int dim;
    int np;
    int nvar;
    vector<string> varNameVec;
    int pos;
    if (mpi_rank == 0) {
      ifstream ifp(filename.c_str());
      assert(!ifp.fail());
      ifp >> version;

      if ((version != 5) && (version != 4)) {
        COUT1("only binary formats (version 4 or 5) of Fluent IP files are supported currently");
        return -1;
      }

      ifp >> dim;
      if (dim != 3) {
        COUT1("only 3D Fluent IP files are supported currently");
        return -2;
      }

      ifp >> np;
      ifp >> nvar;
      varNameVec.push_back("x");
      varNameVec.push_back("y");
      varNameVec.push_back("z");

      for (int ivar = 0; ivar < nvar; ++ivar) {
        string var;
        ifp >> var;
        varNameVec.push_back(var);
      }
      nvar += 3; // account for x,y,z, but don't want in the above loop

      pos = ifp.tellg();
      cout << " > version: " << version << ", dimension: " << dim << ", points: " << np << ", variables: " << nvar << ", header_size (bytes): " << pos << endl;
      cout << " > variables: ";
      for (int ii = 0, lim = varNameVec.size(); ii < lim; ++ii) {
        cout << varNameVec[ii];;
        if (ii < lim-1)
          cout << ", ";
      }
      cout << endl;
    }

    // share header info
    int int_buf_size = 4+varNameVec.size();
    MPI_Bcast(&int_buf_size,1,MPI_INT,0,mpi_comm);
    int* int_buf = new int[int_buf_size];
    int char_buf_size = 0;
    if (mpi_rank == 0) {
      int_buf[0] = np;
      int_buf[1] = nvar;
      int_buf[2] = pos;
      int_buf[3] = version;
      for (int ivar = 0; ivar < nvar; ++ivar) {
        char_buf_size += varNameVec[ivar].length()+1;
        int_buf[4+ivar] = varNameVec[ivar].length()+1; // +1 to include null character
      }
    }
    MPI_Bcast(int_buf,int_buf_size,MPI_INT,0,mpi_comm);
    MPI_Bcast(&char_buf_size,1,MPI_INT,0,mpi_comm);
    char* char_buf = new char[char_buf_size];
    if (mpi_rank == 0) {
      int char_buf_disp = 0;
      for (int ivar = 0; ivar < nvar; ++ivar) {
        strcpy(char_buf+char_buf_disp,varNameVec[ivar].c_str());
        char_buf_disp += varNameVec[ivar].length()+1;
      }
      assert(char_buf_disp == char_buf_size);
    }
    MPI_Bcast(char_buf,char_buf_size,MPI_CHAR,0,mpi_comm);
    if (mpi_rank != 0) {
      int int_buf_disp = 0;
      np = int_buf[int_buf_disp++];
      nvar = int_buf[int_buf_disp++];
      pos = int_buf[int_buf_disp++];
      version = int_buf[int_buf_disp++];
      int char_buf_disp = 0;
      for (int ivar = 0; ivar < nvar; ++ivar) {
        const int length = int_buf[int_buf_disp++];
        varNameVec.push_back(string(char_buf+char_buf_disp,length-1));
        char_buf_disp += length;
      }
      assert(int_buf_disp == int_buf_size);
      assert(char_buf_disp == char_buf_size);
    }
    delete[] int_buf;
    delete[] char_buf;

    assert(ncv_global == 0); ncv_global = np;
    assert(cvora == NULL);
    MiscUtils::calcUniformDist(cvora,ncv_global,mpi_size);
    assert(ncv == 0); ncv = cvora[mpi_rank+1] - cvora[mpi_rank];
    assert(x_cv == NULL); x_cv = new double[ncv][3];
    vector<pair<double *,int> > dataVec(nvar);
    FOR_I3 {
      dataVec[i].first = (double*)x_cv;
      dataVec[i].second = i;
    }
    for (int ivar = 3; ivar < nvar; ++ivar) {
      if ((varNameVec[ivar] == "pressure")&&(varNameSet.find("p") != varNameSet.end())) {
        dnList.push_back(DnData("p"));
        assert(dnList.back().data == NULL);
        dnList.back().data = new double[ncv];
        dataVec[ivar].first = dnList.back().data;
        dataVec[ivar].second = -1; // scalar
      }
      // assumes velocity chunks are ordered and contiguous
      else if ((varNameVec[ivar] == "x-velocity")&&(varNameSet.find("u") != varNameSet.end())) {
        dn3List.push_back(Dn3Data("u"));
        assert(dn3List.back().data == NULL);
        dn3List.back().data = new double[ncv][3];
        dataVec[ivar].first = (double*)dn3List.back().data;
        dataVec[ivar].second = 0;
      }
      else if ((varNameVec[ivar] == "y-velocity")&&(varNameSet.find("u") != varNameSet.end())) {
        // vector data must come in order for now
        dataVec[ivar].first = (double*)dn3List.back().data;
        dataVec[ivar].second = 1;
      }
      else if ((varNameVec[ivar] == "z-velocity")&&(varNameSet.find("u") != varNameSet.end())) {
        dataVec[ivar].first = (double*)dn3List.back().data;
        dataVec[ivar].second = 2;
      }
      else if ((varNameVec[ivar] == "temperature")&&(varNameSet.find("T") != varNameSet.end())) {
        dnList.push_back(DnData("T"));
        assert(dnList.back().data == NULL);
        dnList.back().data = new double[ncv];
        dataVec[ivar].first = dnList.back().data;
        dataVec[ivar].second = -1;
      }
      //else if ((varNameVec[ivar] == "enthalpy")&&(varNameSet.find("h") != varNameSet.end())) {
      //  dnList.push_back(DnData("h"));
      //  assert(dnList.back().data == NULL);
      //  dnList.back().data = new double[ncv];
      //  dataVec[ivar].first = dnList.back().data;
      //  dataVec[ivar].second = -1;
      //}
      else if ((varNameVec[ivar] == "fmean")&&(varNameSet.find("Z") != varNameSet.end())) {
        dnList.push_back(DnData("Z"));
        assert(dnList.back().data == NULL);
        dnList.back().data = new double[ncv];
        dataVec[ivar].first = dnList.back().data;
        dataVec[ivar].second = -1;
      }
      //else if ((varNameVec[ivar] == "premixc")&&(varNameSet.find("C") != varNameSet.end())) {
      //  dnList.push_back(DnData("C")); // needs to be rescaled by solver
      //  assert(dnList.back().data == NULL);
      //  dnList.back().data = new double[ncv];
      //  dataVec[ivar].first = dnList.back().data;
      //  dataVec[ivar].second = -1;
      //}
      else {
        dataVec[ivar].first = NULL; // skip
        dataVec[ivar].second = -1;
      }
    }

    char dummy[128]; assert(filename.length() < 128);
    sprintf(dummy,"%s",filename.c_str());
    MPI_File fh;
    MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
    MPI_Offset offset = pos+2; // skip to first double record


    bool byte_swap = false;  // ip file doesn't specify anything in header to indicate swap vs no-swap

    if (version == 5) {
      double *double_buf = new double[ncv];
      // x y z
      for (int ivar = 0; ivar < 3; ++ivar) {
        readChunkedData<double>(fh,offset+cvora[mpi_rank]*double_size,double_buf,ncv,byte_swap,mpi_comm);
        FOR_ICV x_cv[icv][ivar] = double_buf[icv];
        MiscUtils::dumpRange(double_buf,ncv,varNameVec[ivar]);
        offset += double_size*ncv_global+33;
      }
      // remaining variables
      for (int ivar = 3; ivar < nvar; ++ivar) {
        if (dataVec[ivar].first != NULL) {
          if (dataVec[ivar].second == -1) {
            readChunkedData<double>(fh,offset+cvora[mpi_rank]*double_size,dataVec[ivar].first,ncv,byte_swap,mpi_comm);
            MiscUtils::dumpRange(dataVec[ivar].first,ncv,varNameVec[ivar]);
          }
          else {
            assert((dataVec[ivar].second >= 0)&&(dataVec[ivar].second < 3));
            readChunkedData<double>(fh,offset+cvora[mpi_rank]*double_size,double_buf,ncv,byte_swap,mpi_comm);
            FOR_ICV ((double(*)[3])dataVec[ivar].first)[icv][dataVec[ivar].second] = double_buf[icv];
            MiscUtils::dumpRange(double_buf,ncv,varNameVec[ivar]);
          }
        }
        offset += double_size*ncv_global+33;
      }
      delete[] double_buf;
    }
    else if (version == 4) {
      float *float_buf = new float[ncv];
      // x y z
      for (int ivar = 0; ivar < 3; ++ivar) {
        readChunkedData<float>(fh,offset+cvora[mpi_rank]*float_size,float_buf,ncv,byte_swap,mpi_comm);
        FOR_ICV x_cv[icv][ivar] = double(float_buf[icv]);
        MiscUtils::dumpRange(float_buf,ncv,varNameVec[ivar]);
        offset += float_size*ncv_global+33;
      }
      // remaining variables
      for (int ivar = 3; ivar < nvar; ++ivar) {
        if (dataVec[ivar].first != NULL) {
          if (dataVec[ivar].second == -1) {
            readChunkedData<float>(fh,offset+cvora[mpi_rank]*float_size,float_buf,ncv,byte_swap,mpi_comm);
            FOR_ICV dataVec[ivar].first[icv] = double(float_buf[icv]);
            MiscUtils::dumpRange(dataVec[ivar].first,ncv,varNameVec[ivar]);
          }
          else {
            assert((dataVec[ivar].second >= 0)&&(dataVec[ivar].second < 3));
            readChunkedData<float>(fh,offset+cvora[mpi_rank]*float_size,float_buf,ncv,byte_swap,mpi_comm);
            FOR_ICV ((double(*)[3])dataVec[ivar].first)[icv][dataVec[ivar].second] = double(float_buf[icv]);
            MiscUtils::dumpRange(float_buf,ncv,varNameVec[ivar]);
          }
        }
        offset += float_size*ncv_global+33;
      }
      delete[] float_buf;
    }
    else {
      assert(0);  // should have been caught above
    }

    MPI_File_close(&fh);

    return 0;
  }

  int initFromRestart(const string& mles_filename,const string& sles_filename,const string& csVarNames) {

    // csVarNames == comma-separated var names...

    // turn csVarNames into a set of strings...

    set<string> varNameSet;
    MiscUtils::splitCsv(varNameSet,csVarNames);

    return initFromRestart(mles_filename,sles_filename,varNameSet);

  }

  int initFromRestart(const string& mles_filename,const string& sles_filename,const set<string>& varNameSet) {

    COUT1("DataContainer::initFromRestart(\"" << mles_filename << "\",\"" << sles_filename << "\")...");

    bool b_foundHash = false;

    // 1.0 get x_cv from mles...

    MPI_File mles_fh;
    char dummy[128];
    sprintf(dummy,"%s",mles_filename.c_str());
    int ierr = MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&mles_fh);

    if (ierr != 0) {
      if (mpi_rank == 0) cout << "Error: cannot open mles file: " << mles_filename << endl;
      return -1;
    }

    int itmp[2] = { 0, 0 };
    if (mpi_rank == 0) {
      // first 2 ints are: 0. magic number, 1. io version
      MPI_File_read(mles_fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
    }
    MPI_Bcast(itmp,2,MPI_INT,0,mpi_comm);

    bool byte_swap = false;
    if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
      ByteSwap::byteSwap(itmp,2);
      if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
	if (mpi_rank == 0) cout << "Error: mles file does not start as expected." << endl;
	return -1;
      }
      if (mpi_rank == 0) cout << " > file requires byte swapping." << endl;
      byte_swap = true;
    }

    int io_version = itmp[1];
    if (!(io_version >= 5)) {
      if (mpi_rank == 0) cout << "Error: restart file version not >= 5: " << io_version << endl;
      return -1;
    }

    Header header;
    MPI_Offset offset = 8; // 2 ints
    int done = 0;
    while (done != 1) {

      if (mpi_rank == 0) {
	MPI_File_read_at(mles_fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
	if (byte_swap) ByteSwap::byteSwapHeader(&header,1);
      }
      MPI_Bcast(&header,1,MPI_Header,0,mpi_comm);

      switch (header.id) {

      case UGP_IO_HASHIDS:

        //read hashes from mles to perform consistency check with sles
        b_foundHash = true;
        RestartHashUtilities::mlesInterpReadHashes(mles_fh, offset, header); //will bCast to all ranks
        if (mpi_rank==0){
          cout << " > interp mles hash: " << RestartHashUtilities::mlesInterpHash << endl;
        }
        break;

      case UGP_IO_I0:

        if (mpi_rank == 0) cout << " > skipping I0    \"" << header.name << "\" value: " << header.idata[0] << endl;
	break;

      case UGP_IO_D0:

	if (mpi_rank == 0) cout << " > skipping D0    \"" << header.name << "\" value: " << header.rdata[0] << endl;
	break;

      case UGP_IO_CV_D1:

        if (mpi_rank == 0) cout << " > skipping CV_D1 \"" << header.name << "\"" << endl;
	break;

      case UGP_IO_CV_D2:

	if (strcmp(header.name,"x_vv") == 0) {

	  buildOrCheckCvora(ByteSwap::getInt8FromLswMswPair(header.idata+0));

	  if (mpi_rank == 0) cout << " > reading  CV_D2 \"" << header.name << "\"";

	  assert(x_cv == NULL);
	  x_cv = new double[ncv][3];

	  readChunkedData<double>(mles_fh,offset+header_size+cvora[mpi_rank]*double_size*3,(double*)x_cv,ncv*3,byte_swap,mpi_comm);

	  MiscUtils::dumpRange(x_cv,ncv,"");

	}
	else {

	  if (mpi_rank == 0) cout << " > skipping CV_D2 \"" << header.name << "\"" << endl;

	}
	break;

      case UGP_IO_EOF:

	done = 1;
        break;

	/*
	  default:
	  if (mpi_rank == 0) cout << " > skipping header: " << header.id << " \"" << header.name << "\"." << endl;
	*/

      }

      offset += header.skip;

    }

    MPI_File_close(&mles_fh);

    if (x_cv == NULL) {
      if (mpi_rank == 0) cout << " > Warning: DataContainer has not initialized x_cv" << endl;
      return -1;
    }

    if (!b_foundHash){ //set a default mles hash if one was not found
      std::stringstream ss;
      ss << "Missing mles Hash";
      RestartHashUtilities::mlesInterpHash.init(ss,RestartHashUtilities::sha1hashlength);
      if (mpi_rank==0) cout << " > setting default interp mles hash id " <<  RestartHashUtilities::mlesInterpHash << endl;
    }
    b_foundHash = false;

    // 2.0 get registered sles data...

    MPI_File sles_fh;
    sprintf(dummy,"%s",sles_filename.c_str());
    ierr = MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&sles_fh);

    if (ierr != 0) {
      if (mpi_rank == 0) cout << "Error: cannot open sles file: " << sles_filename << endl;
      return -1;
    }

    itmp[0] = 0; itmp[1] = 0;
    if (mpi_rank == 0) {
      // first 2 ints are: 0. magic number, 1. io version
      MPI_File_read(sles_fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
    }
    MPI_Bcast(itmp,2,MPI_INT,0,mpi_comm);

    byte_swap = false;
    //cout << UGP_IO_MAGIC_NUMBER << " " << itmp[0] << endl;
    if (itmp[0] != UGP_IO_MAGIC_NUMBER+1) {
      ByteSwap::byteSwap(itmp,2);
      if (itmp[0] != UGP_IO_MAGIC_NUMBER+1) {
	if (mpi_rank == 0) cout << "Error: sles file does not start as expected." << endl;
	return -1;
      }
      if (mpi_rank == 0) cout << " > file requires byte swapping." << endl;
      byte_swap = true;
    }

    io_version = itmp[1];
    if (!(io_version >= 5)) {
      if (mpi_rank == 0) cout << "Error: restart file version not >= 5: " << io_version << endl;
      return -1;
    }

    offset = 8; // 2 ints
    done = 0;
    while (done != 1) {

      if (mpi_rank == 0) {
	MPI_File_read_at(sles_fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
	if (byte_swap) ByteSwap::byteSwapHeader(&header,1);
      }
      MPI_Bcast(&header,1,MPI_Header,0,mpi_comm);

      switch (header.id) {
      case UGP_IO_HASHIDS:
        b_foundHash = true;
        int slesNonMatchFlag;
        //read two hash id's from sles file. First id identifies
        //the current file, store in slesHash.  Second  id
        //identifies the "parent" mles file.  Store in mlesHash
        if (mpi_rank==0)  cout << "RestartHashUtilities::slesInterpReadHashes()" << endl;
        RestartHashUtilities::slesInterpReadHashes(sles_fh, offset, header); //will bCast to all ranks
        if (mpi_rank==0){
          slesNonMatchFlag = RestartHashUtilities::slesInterpConsistencyCheck(); //-1 not enough info, 0 match, 1 no match
          cout << " > Found interp sles hash: " << RestartHashUtilities::slesInterpHash << endl;
          if (slesNonMatchFlag==0){
            cout << " >  with interp mles hash: " << RestartHashUtilities::mlesInterpHash << endl;
          }
          else{
            cout << " > sles expects interp mles hash: " << RestartHashUtilities::myHash << endl;
            cout << " >      current interp mles hash: " << RestartHashUtilities::mlesInterpHash << endl;
          }
        }
        MPI_Bcast(&slesNonMatchFlag,1,MPI_INT,0,mpi_comm);
        if (slesNonMatchFlag>0){
           CERR("interp sles data file is not a match for the existing interp mles file");
        }
        break;

      case UGP_IO_I0:

	if (varNameSet.find(header.name) != varNameSet.end()) {

	  if (mpi_rank == 0) cout << " > reading  I0    \"" << header.name << "\" value: " << header.idata[0] << endl;

	  iList.push_back(pair<string,int>(header.name,header.idata[0]));

	}
	else {

	  if (mpi_rank == 0) cout << " > skipping I0    \"" << header.name << "\" value: " << header.idata[0] << endl;

	}
	break;

      case UGP_IO_D0:

	if (varNameSet.find(header.name) != varNameSet.end()) {

	  if (mpi_rank == 0) cout << " > reading  D0    \"" << header.name << "\" value: " << header.rdata[0] << endl;

	  dList.push_back(pair<string,double>(header.name,header.rdata[0]));

	}
	else {

	  if (mpi_rank == 0) cout << " > skipping D0    \"" << header.name << "\" value: " << header.rdata[0] << endl;

	}
	break;

      case UGP_IO_CV_D1:

	if (varNameSet.find(header.name) != varNameSet.end()) {

	  buildOrCheckCvora(ByteSwap::getInt8FromLswMswPair(header.idata+0));

	  if (mpi_rank == 0) cout << " > reading  CV_D1 \"" << header.name << "\"";

	  dnList.push_back(DnData(header.name));
	  assert(dnList.back().data == NULL);
	  dnList.back().data = new double[ncv];

	  readChunkedData<double>(sles_fh,offset+header_size+cvora[mpi_rank]*double_size,dnList.back().data,ncv,byte_swap,mpi_comm);

	  MiscUtils::dumpRange(dnList.back().data,ncv,"");

	}
	else {

	  if (mpi_rank == 0) cout << " > skipping CV_D1 \"" << header.name << "\"" << endl;

	}
	break;

      case UGP_IO_CV_D2:

	if (varNameSet.find(header.name) != varNameSet.end()) {

	  buildOrCheckCvora(ByteSwap::getInt8FromLswMswPair(header.idata+0));

	  if (mpi_rank == 0) cout << " > reading  CV_D2 \"" << header.name << "\"";

	  dn3List.push_back(Dn3Data(header.name));
	  assert(dn3List.back().data == NULL);
	  dn3List.back().data = new double[ncv][3];

	  readChunkedData<double>(sles_fh,offset+header_size+cvora[mpi_rank]*double_size*3,(double*)dn3List.back().data,ncv*3,byte_swap,mpi_comm);

	  MiscUtils::dumpRange(dn3List.back().data,ncv,"");

	}
	else {

	  if (mpi_rank == 0) cout << " > skipping CV_D2 \"" << header.name << "\"" << endl;

	}
	break;

      case UGP_IO_EOF:

	done = 1;
        break;

	/*
	  default:
	  if (mpi_rank == 0) cout << " > skipping header: " << header.id << " \"" << header.name << "\"." << endl;
	*/

      }

      offset += header.skip;

    }

    MPI_File_close(&sles_fh);

    if (!b_foundHash){ //set a default sles hash if one was not found
      std::stringstream ss;
      ss << "Missing sles Hash";
      RestartHashUtilities::slesInterpHash.init(ss,RestartHashUtilities::sha1hashlength);
      if (mpi_rank==0) cout << " > setting default interp sles hash id " <<  RestartHashUtilities::slesInterpHash << endl;
    }

    return 0;

  }

  int initFromSles(const string& sles_filename,const set<string>& varNameSet) {

    COUT1("DataContainer::initFromSles(\"" << sles_filename << "\")...");

    bool b_foundHash = false;

    // registered sles data...

    MPI_File sles_fh;
    char dummy[128];
    sprintf(dummy,"%s",sles_filename.c_str());
    int ierr = MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&sles_fh);

    if (ierr != 0) {
      if (mpi_rank == 0) cout << "Error: cannot open sles file: " << sles_filename << endl;
      return -1;
    }

    int itmp[2] = { 0, 0 };
    if (mpi_rank == 0) {
      // first 2 ints are: 0. magic number, 1. io version
      MPI_File_read(sles_fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
    }
    MPI_Bcast(itmp,2,MPI_INT,0,mpi_comm);

    bool byte_swap = false;
    //cout << UGP_IO_MAGIC_NUMBER << " " << itmp[0] << endl;
    if (itmp[0] != UGP_IO_MAGIC_NUMBER+1) {
      ByteSwap::byteSwap(itmp,2);
      if (itmp[0] != UGP_IO_MAGIC_NUMBER+1) {
	if (mpi_rank == 0) cout << "Error: sles file does not start as expected." << endl;
	return -1;
      }
      if (mpi_rank == 0) cout << " > file requires byte swapping." << endl;
      byte_swap = true;
    }

    const int io_version = itmp[1];
    if (!(io_version >= 5)) {
      if (mpi_rank == 0) cout << "Warning: sles file version not >= 5: " << io_version << endl;
    }

    Header header;
    MPI_Offset offset = 8; // 2 ints
    int done = 0;
    while (done != 1) {

      if (mpi_rank == 0) {
	MPI_File_read_at(sles_fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
	if (byte_swap) ByteSwap::byteSwapHeader(&header,1);
      }
      MPI_Bcast(&header,1,MPI_Header,0,mpi_comm);

      switch (header.id) {

      case UGP_IO_HASHIDS:

	b_foundHash = true;
        //read two hash id's from sles file. First id identifies
        //the current file, store in slesInterpHash.  Second  id
        //identifies the "parent" mles file and is not used, and thus cleared
	// from myHash below...
        if (mpi_rank==0)  cout << "RestartHashUtilities::slesInterpReadHashes()" << endl;
        RestartHashUtilities::slesInterpReadHashes(sles_fh, offset, header); //will bCast to all ranks
	RestartHashUtilities::myHash.clear();
        break;

      case UGP_IO_I0:

	if (varNameSet.find(header.name) != varNameSet.end()) {

	  if (mpi_rank == 0) cout << " > reading  I0    \"" << header.name << "\" value: " << header.idata[0] << endl;

	  iList.push_back(pair<string,int>(header.name,header.idata[0]));

	}
	else {

	  if (mpi_rank == 0) cout << " > skipping I0    \"" << header.name << "\" value: " << header.idata[0] << endl;

	}
	break;

      case UGP_IO_D0:

	if (varNameSet.find(header.name) != varNameSet.end()) {

	  if (mpi_rank == 0) cout << " > reading  D0    \"" << header.name << "\" value: " << header.rdata[0] << endl;

	  dList.push_back(pair<string,double>(header.name,header.rdata[0]));

	}
	else {

	  if (mpi_rank == 0) cout << " > skipping D0    \"" << header.name << "\" value: " << header.rdata[0] << endl;

	}
	break;

      case UGP_IO_CV_D1:

	if (varNameSet.find(header.name) != varNameSet.end()) {

	  buildOrCheckCvora(ByteSwap::getInt8FromLswMswPair(header.idata+0));

	  if (mpi_rank == 0) cout << " > reading  CV_D1 \"" << header.name << "\"";

	  dnList.push_back(DnData(header.name));
	  assert(dnList.back().data == NULL);
	  dnList.back().data = new double[ncv];

	  readChunkedData<double>(sles_fh,offset+header_size+cvora[mpi_rank]*double_size,dnList.back().data,ncv,byte_swap,mpi_comm);

	  MiscUtils::dumpRange(dnList.back().data,ncv,"");

	}
	else {

	  if (mpi_rank == 0) cout << " > skipping CV_D1 \"" << header.name << "\"" << endl;

	}
	break;

      case UGP_IO_CV_F1:

	if (varNameSet.find(header.name) != varNameSet.end()) {

	  buildOrCheckCvora(ByteSwap::getInt8FromLswMswPair(header.idata+0));

	  if (mpi_rank == 0) cout << " > reading  CV_F1 \"" << header.name << "\"";

	  dnList.push_back(DnData(header.name));
	  assert(dnList.back().data == NULL);
	  dnList.back().data = new double[ncv];

	  float * fdata = new float[ncv];
	  readChunkedData<float>(sles_fh,offset+header_size+cvora[mpi_rank]*float_size,fdata,ncv,byte_swap,mpi_comm);
	  for (int icv = 0; icv < ncv; ++icv)
	    dnList.back().data[icv] = double(fdata[icv]);
	  delete[] fdata;

	  MiscUtils::dumpRange(dnList.back().data,ncv,"");

	}
	else {

	  if (mpi_rank == 0) cout << " > skipping CV_F1 \"" << header.name << "\"" << endl;

	}
	break;

      case UGP_IO_CV_D2:

	if (strcmp(header.name,"x_vv") == 0) {

	  // this is th record containing x_vv data...

	  buildOrCheckCvora(ByteSwap::getInt8FromLswMswPair(header.idata+0));

	  if (mpi_rank == 0) cout << " > reading  CV_D2 \"" << header.name << "\"";

	  // we push this into the data container's x_cv

	  assert(x_cv == NULL);
	  x_cv = new double[ncv][3];

	  readChunkedData<double>(sles_fh,offset+header_size+cvora[mpi_rank]*double_size*3,(double*)x_cv,ncv*3,byte_swap,mpi_comm);

	  MiscUtils::dumpRange(x_cv,ncv,"");

	}
	else if (varNameSet.find(header.name) != varNameSet.end()) {

	  buildOrCheckCvora(ByteSwap::getInt8FromLswMswPair(header.idata+0));

	  if (mpi_rank == 0) cout << " > reading  CV_D2 \"" << header.name << "\"";

	  dn3List.push_back(Dn3Data(header.name));
	  assert(dn3List.back().data == NULL);
	  dn3List.back().data = new double[ncv][3];

	  readChunkedData<double>(sles_fh,offset+header_size+cvora[mpi_rank]*double_size*3,(double*)dn3List.back().data,ncv*3,byte_swap,mpi_comm);

	  MiscUtils::dumpRange(dn3List.back().data,ncv,"");

	}
	else {

	  if (mpi_rank == 0) cout << " > skipping CV_D2 \"" << header.name << "\"" << endl;

	}
	break;

      case UGP_IO_CV_F2:

	if (strcmp(header.name,"x_vv") == 0) {

	  // this is th record containing x_vv data...

	  buildOrCheckCvora(ByteSwap::getInt8FromLswMswPair(header.idata+0));

	  if (mpi_rank == 0) cout << " > reading  CV_F2 \"" << header.name << "\"";

	  // we push this into the data container's x_cv

	  assert(x_cv == NULL);
	  x_cv = new double[ncv][3];

	  float (*fdata)[3] = new float[ncv][3];
	  readChunkedData<float>(sles_fh,offset+header_size+cvora[mpi_rank]*float_size*3,(float*)fdata,ncv*3,byte_swap,mpi_comm);
	  for (int icv = 0; icv < ncv; ++icv)
	    FOR_I3 x_cv[icv][i] = double(fdata[icv][i]);
	  delete[] fdata;

	  MiscUtils::dumpRange(x_cv,ncv,"");

	}
	else if (varNameSet.find(header.name) != varNameSet.end()) {

	  buildOrCheckCvora(ByteSwap::getInt8FromLswMswPair(header.idata+0));

	  if (mpi_rank == 0) cout << " > reading  CV_F2 \"" << header.name << "\"";

	  dn3List.push_back(Dn3Data(header.name));
	  assert(dn3List.back().data == NULL);
	  dn3List.back().data = new double[ncv][3];

	  float (*fdata)[3] = new float[ncv][3];
	  readChunkedData<float>(sles_fh,offset+header_size+cvora[mpi_rank]*float_size*3,(float*)fdata,ncv*3,byte_swap,mpi_comm);
	  for (int icv = 0; icv < ncv; ++icv)
	    FOR_I3 dn3List.back().data[icv][i] = double(fdata[icv][i]);
	  delete[] fdata;

	  MiscUtils::dumpRange(dn3List.back().data,ncv,"");

	}
	else {

	  if (mpi_rank == 0) cout << " > skipping CV_F2 \"" << header.name << "\"" << endl;

	}
	break;

      case UGP_IO_EOF:

	done = 1;
        break;

	/*
	  default:
	  if (mpi_rank == 0) cout << " > skipping header: " << header.id << " \"" << header.name << "\"." << endl;
	*/

      }

      offset += header.skip;

    }

    MPI_File_close(&sles_fh);

    if (!b_foundHash){ //set a default sles hash if one was not found
      std::stringstream ss;
      ss << "Missing sles Hash";
      RestartHashUtilities::slesInterpHash.init(ss,RestartHashUtilities::sha1hashlength);
      if (mpi_rank==0) cout << " > setting default interp sles hash id " <<  RestartHashUtilities::slesInterpHash << endl;
    }

    if (x_cv == NULL) {
      if (mpi_rank == 0) cout << " > Warning: DataContainer has not initialized x_cv. You will need the matching mles file." << endl;
      return -1;
    }

    return 0;

  }

  int initFromOldRestart(const string& filename,const string& csVarNames) {

    // csVarNames == comma-separated var names...

    // turn csVarNames into a set of strings...

    set<string> varNameSet;
    MiscUtils::splitCsv(varNameSet,csVarNames);

    return initFromOldRestart(filename,varNameSet);

  }

  int initFromOldRestart(const string& filename,const set<string>& varNameSet) {

    COUT1("DataContainer::initFromOldRestart(\"" << filename << "\")...");

    MPI_File fh;
    char dummy[128];
    sprintf(dummy,"%s",filename.c_str());
    int ierr = MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

    if (ierr != 0) {
      if (mpi_rank == 0) cout << "Error: cannot open old restart file: " << filename << endl;
      return -1;
    }

    int itmp[2] = { 0, 0 };
    if (mpi_rank == 0) {
      // first 2 ints are: 0. magic number, 1. io version
      MPI_File_read(fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
    }
    MPI_Bcast(itmp,2,MPI_INT,0,mpi_comm);

    bool byte_swap = false;
    if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
      ByteSwap::byteSwap(itmp,2);
      if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
	if (mpi_rank == 0) cout << "Error: restart file does not start as expected." << endl;
	return -1;
      }
      if (mpi_rank == 0) cout << " > file requires byte swapping." << endl;
      byte_swap = true;
    }

    const int io_version = itmp[1];
    if (!((io_version == 3)||(io_version == 4))) {
      if (mpi_rank == 0) cout << "Error: restart file version not 3 or 4: " << io_version << endl;
      return -1;
    }

    map<string,string> varNameMap; //if io3 fill map with varNameMap[VARNAME] = varName
    for (set<string>::iterator it = varNameSet.begin(); it!=varNameSet.end(); ++it){
      string varName = (*it);
      if (io_version==3){
        transform(varName.begin(), varName.end(), varName.begin(), ::toupper);
      }
      varNameMap[varName] = *it;
    }


    Header header;
    MPI_Offset offset = 8; // 2 ints
    int done = 0;
    while (done != 1) {

      if (mpi_rank == 0) {
	MPI_File_read_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
	if (byte_swap) ByteSwap::byteSwapHeader(&header,1);
      }
      MPI_Bcast(&header,1,MPI_Header,0,mpi_comm);

      switch (header.id) {
      case UGP_IO_I0:

	if (varNameMap.find(header.name) != varNameMap.end()) {

	  if (mpi_rank == 0) cout << " > reading  I0    \"" << header.name << "\" value: " << header.idata[0] << endl;

	  iList.push_back(pair<string,int>(varNameMap[header.name],header.idata[0]));

	}
	else {

	  if (mpi_rank == 0) cout << " > skipping I0    \"" << header.name << "\" value: " << header.idata[0] << endl;

	}
	break;

      case UGP_IO_D0:

	if (varNameMap.find(header.name) != varNameMap.end()) {

	  if (mpi_rank == 0) cout << " > reading  D0    \"" << header.name << "\" value: " << header.rdata[0] << endl;

	  dList.push_back(pair<string,double>(varNameMap[header.name],header.rdata[0]));

	}
	else {

	  if (mpi_rank == 0) cout << " > skipping D0    \"" << header.name << "\" value: " << header.rdata[0] << endl;

	}
	break;

      case UGP_IO_CV_D1:

        if (varNameMap.find(header.name) != varNameMap.end()) {

	  buildOrCheckCvora(ByteSwap::getInt8FromLswMswPair(header.idata+0));

	  if (mpi_rank == 0) cout << " > reading  CV_D1 \"" << header.name << "\"";

	  dnList.push_back(DnData(varNameMap[header.name]));
	  assert(dnList.back().data == NULL);
	  dnList.back().data = new double[ncv];

	  readChunkedData<double>(fh,offset+header_size+cvora[mpi_rank]*double_size,dnList.back().data,ncv,byte_swap,mpi_comm);

	  MiscUtils::dumpRange(dnList.back().data,ncv,"");

	}
	else {

	  if (mpi_rank == 0) cout << " > skipping CV_D1 \"" << header.name << "\"" << endl;

	}
	break;

      case UGP_IO_CV_D2:

	if (strcmp(header.name,"x_vv") == 0 || strcmp(header.name,"X_VV") == 0) {

          int8 ncv_global_check;
          if (io_version==3){
            ncv_global_check = header.idata[0]; //for io v3 idata[0] = ncv_global, idata[1]=3
          }
          else{
            ncv_global_check = ByteSwap::getInt8FromLswMswPair(header.idata+0);
          }
          buildOrCheckCvora(ncv_global_check);


	  if (mpi_rank == 0) cout << " > reading  CV_D2 \"" << header.name << "\"";

	  assert(x_cv == NULL);
	  x_cv = new double[ncv][3];

	  readChunkedData<double>(fh,offset+header_size+cvora[mpi_rank]*double_size*3,(double*)x_cv,ncv*3,byte_swap,mpi_comm);

	  MiscUtils::dumpRange(x_cv,ncv,"");

	}
	else if (varNameMap.find(header.name) != varNameMap.end()) {

          int8 ncv_global_check;
          if (io_version==3){
            ncv_global_check = header.idata[0]; //for io v3 idata[0] = ncv_global, idata[1]=3
          }
          else{
            ncv_global_check = ByteSwap::getInt8FromLswMswPair(header.idata+0);
          }
	  buildOrCheckCvora(ncv_global_check);

	  if (mpi_rank == 0) cout << " > reading  CV_D2 \"" << header.name << "\"";

	  dn3List.push_back(Dn3Data(varNameMap[header.name]));
	  assert(dn3List.back().data == NULL);
	  dn3List.back().data = new double[ncv][3];

	  readChunkedData<double>(fh,offset+header_size+cvora[mpi_rank]*double_size*3,(double*)dn3List.back().data,ncv*3,byte_swap,mpi_comm);

	  MiscUtils::dumpRange(dn3List.back().data,ncv,"");

	}
	else {

	  if (mpi_rank == 0) cout << " > skipping CV_D2 \"" << header.name << "\"" << endl;

	}
 	break;

      case UGP_IO_EOF:

	done = 1;
        break;

	/*
	  default:
	  if (mpi_rank == 0) cout << " > skipping header: " << header.id << " \"" << header.name << "\"." << endl;
	*/

      }

      offset += header.skip;

    }

    MPI_File_close(&fh);

    /*
      MPI_Barrier(mpi_comm);
      if (mpi_rank == 0) {
      const double seconds = MPI_Wtime() - wtime0;
      cout << " > read_mles summary: ncv_global: " << ncv_global <<
      " file size: " << double(offset)/1.0E+9 << " [GB] read rate: " << double(offset)/1.0E+9/seconds << " [GB/s]" << endl;
      }
    */

    if (x_cv == NULL) {
      if (mpi_rank == 0) cout << " > Warning: DataContainer has not initialized x_cv" << endl;
      return -1;
    }

    return 0;

  }


  void writePointsBinary(const string& filename) const {

    COUT1("writePointsBinary(): " << filename);
    assert(delta != NULL); // make sure we populated delta (maybe we can make it optional)

    char dummy[128]; assert(filename.length() < 128);
    sprintf(dummy,"%s",filename.c_str());
    MPI_File_delete(dummy,MPI_INFO_NULL);

    MPI_File fh;
    MPI_File_open(mpi_comm,dummy,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    if (mpi_rank == 0) {
      int8 ibuf[4] = { POINTS_IO_MAGIC_NUMBER, POINTS_IO_VERSION, ncv_global, 1 }; // last is include delta or not maybe?
      MPI_File_write_at(fh,0,ibuf,4,MPI_INT8,MPI_STATUS_IGNORE);
    }

    // x_cv...

    int8 my_ncv = ncv;
    int8 my_disp;
    MPI_Scan(&my_ncv,&my_disp,1,MPI_INT8,MPI_SUM,mpi_comm);
    assert( (mpi_rank != mpi_size-1) || (my_disp == ncv_global) );
    MPI_Offset offset = (my_disp-my_ncv)*3*8 + 4*8;
    writeChunkedData<double>(fh,offset,(double*)x_cv,ncv*3,mpi_comm);

    // delta...

    offset = ncv_global*3*8 + 4*8 + (my_disp-my_ncv)*8;
    writeChunkedData<double>(fh,offset,delta,ncv,mpi_comm);

    // and close...

    MPI_File_close(&fh);

  }

  void readPointsBinary(const string& filename) {

    COUT1("readPointsBinary(): " << filename);

    char dummy[128]; assert(filename.length() < 128);
    sprintf(dummy,"%s",filename.c_str());
    MPI_File fh;
    MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

    // recall file starts: { POINTS_IO_MAGIC_NUMBER, POINTS_IO_VERSION, ncv_global, 1 }; // last is include delta or not maybe?
    int8 ibuf[4];
    if (mpi_rank == 0)
      MPI_File_read_at(fh,0,ibuf,4,MPI_INT8,MPI_STATUS_IGNORE);
    MPI_Bcast(ibuf,4,MPI_INT8,0,mpi_comm);
    const bool byte_swap = false;
    assert(ibuf[0] == POINTS_IO_MAGIC_NUMBER); // implement byteswapping if you hit this
    assert(ibuf[1] == POINTS_IO_VERSION);
    ncv_global = ibuf[2];
    assert(ibuf[3] == 1); // where we can specify if delta exists

    MiscUtils::calcUniformDist(cvora,ncv_global,mpi_size);
    assert(cvora[mpi_rank+1]-cvora[mpi_rank] < TWO_BILLION);
    ncv = int(cvora[mpi_rank+1]-cvora[mpi_rank]);

    // x_cv...

    assert(x_cv == NULL); x_cv = new double[ncv][3];
    MPI_Offset offset = cvora[mpi_rank]*3*8 + 4*8;
    readChunkedData<double>(fh,offset,(double*)x_cv,ncv*3,byte_swap,mpi_comm);

    // delta...

    assert(delta == NULL); delta = new double[ncv];
    offset = ncv_global*3*8 + 4*8 + cvora[mpi_rank]*8;
    readChunkedData<double>(fh,offset,delta,ncv,byte_swap,mpi_comm);

    // and close...

    MPI_File_close(&fh);

    if (mpi_rank == 0) cout << " > ncv_global: " << ncv_global << endl;
    MiscUtils::dumpRange(x_cv,ncv,"x_cv from pbin");
    MiscUtils::dumpRange(delta,ncv,"delta from pbin");
  }

};

#endif
