#ifndef _POINTS_HPP_
#define _POINTS_HPP_

#include "CTI.hpp"
using namespace CTI;

#include "MiscUtils.hpp"
using namespace MiscUtils;

#include "Adt.hpp"
#include "ByteSwap.hpp"
#include "CommonIo.hpp"

class Points {
public:

  // Notes: Should Points maintain a ipora[] structure
  // based on its current striping? Most of the points don't require
  // this, and it can be heavy when there are MANY ranks, so just store
  // np_global, and build the ipora if/when required...

  int8 np_global;
  int np;
  double (*xp)[3];
  double *delta;

  // Jan 2020: introduced a 
  bool b_hcp;
  int hcp_packing;
  double hcp_delta_factor,hcp_delta;
  int *level;
  
  int *flag;

  // for fast point location...

  Adt<double> * xpBbAdt;
  Adt<double> * xpAdt;

private:

  void init() {
    np_global = -1;
    np = -1;
    xp = NULL;
    delta = NULL;
    b_hcp = false;
    level = NULL;
    flag = NULL;
    xpBbAdt = NULL;
    xpAdt = NULL;
  }

public:

  Points() {
    init();
  }

  Points(const string& filename) {
    init();
    readBinary(filename);
  }

  void init(const int np) {
    this->np = np;
    int8 np_int8 = np;
    MPI_Allreduce(&np_int8,&np_global,1,MPI_INT8,MPI_SUM,mpi_comm);
    xp = new double[np][3];
    delta = new double[np];
  }

  // copy constructor
  Points(const Points& other) {
    init();

    np = other.np;
    np_global = other.np_global;
    assert(xp == NULL); xp = new double[np][3];
    for (int ip=0; ip < np; ++ip) {
      FOR_I3 xp[ip][i] = other.xp[ip][i];
    }

    assert(delta == NULL); delta = new double[np];
    for (int ip=0; ip < np; ++ip) delta[ip] = other.delta[ip];
  }

  ~Points() {
    DELETE(xp);
    DELETE(delta);
    DELETE(level);
    DELETE(flag);
    if (xpBbAdt != NULL) delete xpBbAdt;
    if (xpAdt != NULL) delete xpAdt;
  }

  int8 getNpGlobal() const {
    assert(np_global >= 0);
    return np_global;
  }

  void setHcpData(const int hcp_packing_,const double hcp_delta_factor_,const double hcp_delta_) {
    b_hcp = true;
    this->hcp_packing      = hcp_packing_;
    this->hcp_delta_factor = hcp_delta_factor_;
    this->hcp_delta        = hcp_delta_;
  }

  void getHcpData(int &hcp_packing_,double &hcp_delta_factor_,double &hcp_delta_) const {
    assert(b_hcp);
    hcp_packing_      = this->hcp_packing;
    hcp_delta_factor_ = this->hcp_delta_factor;
    hcp_delta_        = this->hcp_delta;
  }
  
  void add(const vector<Points*>& ptsPtrVec) {

    np_global = 0;
    np = 0;

    for (int ii = 0; ii < ptsPtrVec.size(); ++ii) {
      if (ptsPtrVec[ii]) {
        np_global += ptsPtrVec[ii]->np_global;
        np += ptsPtrVec[ii]->np;
      }
    }

    assert(xp == NULL);    xp = new double[np][3];
    assert(delta == NULL); delta = new double[np];

    int ip = 0;
    for (int ii = 0; ii < ptsPtrVec.size(); ++ii) {
      if (ptsPtrVec[ii]) {
	for (int this_ip = 0; this_ip < ptsPtrVec[ii]->np; ++this_ip) {
	  FOR_I3 xp[ip][i] = ptsPtrVec[ii]->xp[this_ip][i];
	  delta[ip]        = ptsPtrVec[ii]->delta[this_ip];
	  ++ip;
	}
      }
    }
    assert(ip == np);

  }

  void reorder() {
    // reorders points' data based on the current order in flag...
    assert(xp);
    assert(delta);
    assert(flag);
    MiscUtils::reorder(xp,flag,np);
    MiscUtils::reorder(delta,flag,np);
  }

  bool checkAdt() {

    return( xpAdt != NULL );

  }

  void initAdt() {

    // no ghost data treatment here -- just the owned points...

    assert(xpBbAdt == NULL);
    assert(xpAdt == NULL);

    double bbmin[3] = { 1.0E+20, 1.0E+20, 1.0E+20 };
    double bbmax[3] = { -1.0E+20, -1.0E+20, -1.0E+20 };
    for (int ip = 0; ip < np; ++ip) {
      FOR_I3 bbmin[i] = min(bbmin[i],xp[ip][i]);
      FOR_I3 bbmax[i] = max(bbmax[i],xp[ip][i]);
    }

    double (*bbminGlobal)[3] = new double[mpi_size][3];
    MPI_Allgather(bbmin,3,MPI_DOUBLE,bbminGlobal,3,MPI_DOUBLE,mpi_comm);
    double (*bbmaxGlobal)[3] = new double[mpi_size][3];
    MPI_Allgather(bbmax,3,MPI_DOUBLE,bbmaxGlobal,3,MPI_DOUBLE,mpi_comm);

    xpBbAdt = new Adt<double>(mpi_size,bbminGlobal,bbmaxGlobal);

    delete[] bbminGlobal;
    delete[] bbmaxGlobal;

    // we also need the global data limits...
    /*
      double xpBbminGlobal[3],xpBbmaxGlobal[3]; // the global bounding box for active xp points
      MPI_Allreduce(bbmin,xpBbminGlobal,3,MPI_DOUBLE,MPI_MIN,mpi_comm);
      MPI_Allreduce(bbmax,xpBbmaxGlobal,3,MPI_DOUBLE,MPI_MAX,mpi_comm);
      if (mpi_rank == 0)
      cout << " > global point location bounding box: " << COUT_VEC(xpBbminGlobal) << ":" << COUT_VEC(xpBbmaxGlobal) << endl;
    */

    // finally we can build the xAdt from all points from 0 to nv...
    // recall that periodicity is handled by querying the volume for points
    // at the transformed location -- this guarantees that all potential
    // periodic neighbors can be found,rather than just copying a few "layers"
    // and sending them over to include in a larged box...

    xpAdt = new Adt<double>(np,xp,xp);

  }

  void resetOrClearAfterMove() {

    // we just moved, so...
    clearAdt();

  }

  void clearAdt() {

    if (xpBbAdt != NULL) {
      delete xpBbAdt;
      xpBbAdt = NULL;
    }

    if (xpAdt != NULL) {
      delete xpAdt;
      xpAdt = NULL;
    }

  }

  void translate(const double dx[3]) {
    COUT1("Points::translate: " << COUT_VEC(dx));
    for (int ip = 0; ip < np; ++ip)
      FOR_I3 xp[ip][i] += dx[i];
  }

  void rotate(const double (&point)[3], const double (&_axis)[3], const double angle_deg) {
    double axis[3] = {_axis[0],_axis[1],_axis[2]};
    NORMALIZE(axis);
    const double angle = angle_deg/180.0*M_PI;  // radians

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

    for (int ip = 0; ip < np; ++ip) {
      double _xp[3];
      FOR_J3 _xp[j] = xp[ip][j] - point[j];
      const double temp[3] = MATRIX_PRODUCT(R,_xp);
      FOR_J3 xp[ip][j] = temp[j] + point[j];
    }
  }

  void scale(const double dx[3]) {
    COUT1("Points::scale: " << COUT_VEC(dx));
    for (int ip = 0; ip < np; ++ip) {
      FOR_I3 xp[ip][i] *= dx[i];
      delta[ip] *= max(dx[0],max(dx[1],dx[2]));
    }
  }

  void writeBinary(const string& filename) {
    
    COUT1("Points::writeBinary(): " << filename);

    // check that np_global has been set...
    
    if (np_global == -1) {
      int8 np_int8 = np;
      MPI_Allreduce(&np_int8,&np_global,1,MPI_INT8,MPI_SUM,mpi_comm);
    }
    
    char dummy[128]; assert(filename.length() < 128);
    sprintf(dummy,"%s",filename.c_str());
    MPI_File_delete(dummy,MPI_INFO_NULL);

    MPI_File fh;
    MPI_File_open(mpi_comm,dummy,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    MPI_Offset offset;
    if (b_hcp) {
      
      // the b_hcp format should have level information, and NOT delta...
      assert(xp);
      assert(delta == NULL);
      assert(level);
      
      if (mpi_rank == 0) {
        
        int8 i8buf[4] = { POINTS_IO_MAGIC_NUMBER, POINTS_IO_VERSION, np_global, 2 }; // use a 2 to trigger this kind
        MPI_File_write(fh,i8buf,4,MPI_INT8,MPI_STATUS_IGNORE);
        
        // and write the hcp info...
        int ibuf[2] = { hcp_packing, 0 }; // extra int here to keep alignment
        MPI_File_write(fh,ibuf,2,MPI_INT,MPI_STATUS_IGNORE);
        double dbuf[2] = { hcp_delta_factor, hcp_delta };
        MPI_File_write(fh,dbuf,2,MPI_DOUBLE,MPI_STATUS_IGNORE);

      }
      
      offset = 4*8 + 2*4 + 2*8;

    }
    else {

      assert(xp);
      assert(delta);
      assert(level == NULL);
      
      if (mpi_rank == 0) {

        int8 i8buf[4] = { POINTS_IO_MAGIC_NUMBER, POINTS_IO_VERSION, np_global, 1 }; 
        MPI_File_write(fh,i8buf,4,MPI_INT8,MPI_STATUS_IGNORE);

      }
      
      offset = 4*8;
      
    }

    // xp...

    int8 my_np = np;
    int8 my_disp;
    MPI_Scan(&my_np,&my_disp,1,MPI_INT8,MPI_SUM,mpi_comm);
    assert( (mpi_rank != mpi_size-1) || (my_disp == np_global) );
    writeChunkedData<double>(fh,offset+(my_disp-my_np)*3*8,(double*)xp,3*np,mpi_comm);

    offset += np_global*3*8;

    // level or delta...

    if (b_hcp) {

      writeChunkedData<int>(fh,offset+(my_disp-my_np)*4,level,np,mpi_comm);

    }
    else {

      writeChunkedData<double>(fh,offset+(my_disp-my_np)*8,delta,np,mpi_comm);
      
    }

    // and close...

    MPI_File_close(&fh);

  }

  void writeResult(FILE * fp) const {

    assert(fp);

    int ibuf[2] = {
      POINTS_IO_STAMP,
      np
    };
    fwrite(ibuf,sizeof(int),2,fp);
    fwrite(&np_global,sizeof(int8),1,fp);
    fwrite(xp,sizeof(double),np*3,fp);
    fwrite(delta,sizeof(double),np,fp);

  }

  void readRestart(FILE * fp) {

    assert(fp);

    int ibuf[2];
    fread(ibuf,sizeof(int),2,fp);
    assert(ibuf[0] == POINTS_IO_STAMP);
    np = ibuf[1];

    fread(&np_global,sizeof(int8),1,fp);

    assert(xp == NULL);
    xp = new double[np][3];
    fread(xp,sizeof(double),np*3,fp);

    assert(delta == NULL);
    delta = new double[np];
    fread(delta,sizeof(double),np,fp);

  }

  void readBinary(const string& filename) {

    COUT1("Points::readBinary(): " << filename);

    char dummy[128]; assert(filename.length() < 128);
    sprintf(dummy,"%s",filename.c_str());
    MPI_File fh;
    int ierr = MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
    if (ierr != 0) {
      CERR("Points::readBinary: file not found: \"" << filename << "\"");
    }

    // recall file starts: { POINTS_IO_MAGIC_NUMBER, POINTS_IO_VERSION, np_global, 1 }; // last is include delta or not maybe?
    int8 i8buf[4];
    if (mpi_rank == 0)
      MPI_File_read(fh,i8buf,4,MPI_INT8,MPI_STATUS_IGNORE);
    MPI_Bcast(i8buf,4,MPI_INT8,0,mpi_comm);
    bool byte_swap = false;
    if (i8buf[0] != POINTS_IO_MAGIC_NUMBER) {
      ByteSwap::byteSwap(i8buf,4);
      if (i8buf[0] != POINTS_IO_MAGIC_NUMBER) {
	CERR("Error: pbin file does not start as expected.")
      }
      if (mpi_rank == 0) cout << " > file requires byte swapping." << endl;
      byte_swap = true;
    }
    assert(i8buf[1] == POINTS_IO_VERSION);
    np_global = i8buf[2];
    
    MPI_Offset offset = 4*8;

    if (mpi_rank == 0)
      cout << " > np_global: " << np_global << endl;
    
    // the last int8 triggers any special points file type...
    
    if (i8buf[3] == 2) {
      
      // this is a new hcp points version with some additional info about hcp points in the header...
      int ibuf[2];
      double dbuf[2];
      if (mpi_rank == 0) {
        MPI_File_read(fh,ibuf,2,MPI_INT,MPI_STATUS_IGNORE);
        MPI_File_read(fh,dbuf,2,MPI_DOUBLE,MPI_STATUS_IGNORE);
      }
      MPI_Bcast(ibuf,2,MPI_INT,0,mpi_comm);
      MPI_Bcast(dbuf,3,MPI_DOUBLE,0,mpi_comm);
      
      b_hcp = true;
      hcp_packing = ibuf[0];
      hcp_delta_factor = dbuf[0];
      hcp_delta = dbuf[1];
      if (mpi_rank == 0) cout << " > pbin is hcp format: hcp_delta: " << hcp_delta << endl;
      
      offset += 2*4 + 2*8;

    }
    else {
      
      assert(i8buf[3] == 1);
      
    }
    
    int8 * ipora = NULL;
    MiscUtils::calcUniformDist(ipora,np_global,mpi_size);
    assert(ipora[mpi_rank+1]-ipora[mpi_rank] < TWO_BILLION);
    assert(np == -1);
    np = int(ipora[mpi_rank+1]-ipora[mpi_rank]);
    
    // xp...

    assert(xp == NULL); xp = new double[np][3];
    readChunkedData<double>(fh,offset+ipora[mpi_rank]*3*8,(double*)xp,3*np,byte_swap,mpi_comm);
    offset += np_global*3*8;

    if (b_hcp) {

      // hcp format gets a level next...
      assert(level == NULL);
      level = new int[np];
      readChunkedData<int>(fh,offset+ipora[mpi_rank]*4,level,np,byte_swap,mpi_comm);
      offset += np_global*4;
      
    }
    else {
      
      // regular pbin format reads delta...
      assert(delta == NULL); 
      delta = new double[np];
      readChunkedData<double>(fh,offset+ipora[mpi_rank]*8,delta,np,byte_swap,mpi_comm);
      offset += np_global*8;

    }

    // and close...

    MPI_File_close(&fh);
    
    if (mpi_rank == 0) cout << " > filesize: " << offset << endl; 
    
    // cleanup...
    
    delete[] ipora;
    
    // and report...

    double my_buf[8];
    for (int j = 0; j < 8; ++j) my_buf[j] = HUGE_VAL; // something large
    
    for (int ip = 0; ip < np; ++ip) {
      FOR_J3 {
        my_buf[2*j] = min(my_buf[2*j],xp[ip][j]);
        my_buf[2*j+1] = min(my_buf[2*j+1],-xp[ip][j]);
      }
      if (delta) {
        my_buf[6] = min(my_buf[6],delta[ip]);
        my_buf[7] = min(my_buf[7],-delta[ip]);
      }
      else {
        my_buf[6] = min(my_buf[6],double(level[ip]));
        my_buf[7] = min(my_buf[7],double(-level[ip]));
      }
    }

    double buf[8];
    MPI_Allreduce(my_buf, buf, 8, MPI_DOUBLE, MPI_MIN, mpi_comm);
    if (mpi_rank == 0) {
      cout << " > xp range: ";
      FOR_J3 {
        cout << j << ": " << buf[2*j] << ":" << -buf[2*j+1] << ", ";
      }
      cout << endl;
      if (delta) {
        cout << " > delta range: " << buf[6] << ":" << -buf[7] << endl;
      }
      else {
        assert(level);
        cout << " > level range: " << buf[6] << ":" << -buf[7] << endl;
      }
    }
    
  }

  void addFromBinary(const string& filename) {
    // TODO maybe combine this with above?

    COUT1("addFromBinary(): " << filename);

    char dummy[128]; assert(filename.length() < 128);
    sprintf(dummy,"%s",filename.c_str());
    MPI_File fh;
    int ierr = MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
    if (ierr != 0) {
      CERR("Points::addFromBinary: file not found: \"" << filename << "\"");
    }

    // recall file starts: { POINTS_IO_MAGIC_NUMBER, POINTS_IO_VERSION, np_global, 1 }; // last is include delta or not maybe?
    int8 ibuf[4];
    if (mpi_rank == 0)
      MPI_File_read_at(fh,0,ibuf,4,MPI_INT8,MPI_STATUS_IGNORE);
    MPI_Bcast(ibuf,4,MPI_INT8,0,mpi_comm);
    bool byte_swap = false;
    assert(ibuf[0] == POINTS_IO_MAGIC_NUMBER); // implement byteswapping if you hit this
    assert(ibuf[1] == POINTS_IO_VERSION);
    np_global += ibuf[2];
    assert(ibuf[3] == 1);

    if (mpi_rank == 0)
      cout << " > np_global: " << np_global << endl;

    int8 * ipora = NULL;
    MiscUtils::calcUniformDist(ipora,ibuf[2],mpi_size);
    assert(ipora[mpi_rank+1]-ipora[mpi_rank] < TWO_BILLION);
    const int np0 = np;
    const int this_np = int(ipora[mpi_rank+1]-ipora[mpi_rank]);
    np += this_np;

    // xp...

    double (*xp0)[3] = xp;
    xp = new double[np][3];
    for (int ip = 0; ip < np0; ++ip)
      FOR_I3 xp[ip][i] = xp0[ip][i];
    delete[] xp0;
    MPI_Offset offset = ipora[mpi_rank]*3*8 + 4*8;
    readChunkedData<double>(fh,offset,(double*)(xp)+3*np0,3*this_np,byte_swap,mpi_comm);

    // delta...

    double *delta0 = delta;
    delta = new double[np];
    for (int ip = 0; ip < np0; ++ip)
      delta[ip] = delta0[ip];
    delete[] delta0;
    offset = ibuf[2]*3*8 + 4*8 + ipora[mpi_rank]*8;
    readChunkedData<double>(fh,offset,delta+np0,this_np,byte_swap,mpi_comm);

    // cleanup...

    delete[] ipora;

    // and report...

    double my_buf[8];
    for (int j = 0; j < 8; ++j) my_buf[j] = HUGE_VAL; // something large

    for (int ip = 0; ip < np; ++ip) {
      FOR_J3 {
        my_buf[2*j] = min(my_buf[2*j],xp[ip][j]);
        my_buf[2*j+1] = min(my_buf[2*j+1],-xp[ip][j]);
      }
      my_buf[6] = min(my_buf[6],delta[ip]);
      my_buf[7] = min(my_buf[7],-delta[ip]);
    }

    double buf[8];
    MPI_Allreduce(my_buf, buf, 8, MPI_DOUBLE, MPI_MIN, mpi_comm);
    if (mpi_rank == 0) {
      cout << " > xp range: ";
      FOR_J3 {
        cout << j << ": " << buf[2*j] << ":" << -buf[2*j+1] << ", ";
      }
      cout << endl;
      cout << " > delta range: " << buf[6] << ":" << -buf[7] << endl;
    }

    // and close...

    MPI_File_close(&fh);

  }
  
  void readLes(const string& filename) {

    char dummy[128];
    sprintf(dummy,"%s",filename.c_str());
    MPI_File fh;
    MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
    
    // assume we do not have to byte_swap...
    bool byte_swap = false;

    int itmp[2];
    if (mpi_rank == 0) MPI_File_read(fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
    MPI_Bcast(itmp,2,MPI_INT,0,mpi_comm);

    if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
      ByteSwap::byteSwap(itmp,2);
      if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
	CERR("readLes: file does not start as expected. Exiting.");
      }
      COUT2(" > file requires byte swapping.");
      byte_swap = true;
    }

    // initialize the offset...

    MPI_Offset offset = int_size*2;
    int8 * ipora = NULL;
    assert(xp == NULL);
    assert(delta == NULL);
    
    Header header;
    int done = 0;
    while (done < 2) {
      
      if (mpi_rank == 0) MPI_File_read_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      MPI_Bcast(&header,1,MPI_Header,0,mpi_comm);
      if (byte_swap) ByteSwap::byteSwapHeader(&header,1);

      switch (header.id) {
      case UGP_IO_CV_D2:
	if (strcmp(header.name,"x_vv") == 0) {
	  int8 np_global_check = ByteSwap::getInt8FromLswMswPair(header.idata+0);
	  if (ipora == NULL) {
	    np_global = np_global_check;
	    MiscUtils::calcUniformDist(ipora,np_global,mpi_size);
	    assert(ipora[mpi_rank+1]-ipora[mpi_rank] < TWO_BILLION);
	    assert(np == -1);
	    np = int(ipora[mpi_rank+1]-ipora[mpi_rank]);
	  }
	  assert(np_global_check == ipora[mpi_size]);
	  assert(xp == NULL); xp = new double[np][3];
	  readChunkedData<double>(fh,offset+header_size+ipora[mpi_rank]*3*8,(double*)xp,3*np,byte_swap,mpi_comm);
	}
	break;
	
      case UGP_IO_CV_D1:
	if (strcmp(header.name,"r_vv") == 0) {
	  int8 np_global_check = ByteSwap::getInt8FromLswMswPair(header.idata+0);
	  if (ipora == NULL) {
	    np_global = np_global_check;
	    MiscUtils::calcUniformDist(ipora,np_global,mpi_size);
	    assert(ipora[mpi_rank+1]-ipora[mpi_rank] < TWO_BILLION);
	    assert(np == -1);
	    np = int(ipora[mpi_rank+1]-ipora[mpi_rank]);
	  }
	  assert(np_global_check == ipora[mpi_size]);
	  assert(delta == NULL); delta = new double[np];
	  readChunkedData<double>(fh,offset+header_size+ipora[mpi_rank]*8,delta,np,byte_swap,mpi_comm);
	  // increase delta because it must be typically double r_vv...
	  for (int ip = 0; ip < np; ++ip)
	    delta[ip] *= 2.01;
	}
	break;
	
      case UGP_IO_EOF:
	done = 2;
	break;
	
      }
      
      offset += header.skip;
      
    }
    
    MPI_File_close(&fh);
    
    assert(xp);
    assert(delta);

    // cleanup...

    delete[] ipora;

    // and report...

    double my_buf[8];
    for (int j = 0; j < 8; ++j) my_buf[j] = HUGE_VAL; // something large

    for (int ip = 0; ip < np; ++ip) {
      FOR_J3 {
        my_buf[2*j] = min(my_buf[2*j],xp[ip][j]);
        my_buf[2*j+1] = min(my_buf[2*j+1],-xp[ip][j]);
      }
      my_buf[6] = min(my_buf[6],delta[ip]);
      my_buf[7] = min(my_buf[7],-delta[ip]);
    }

    double buf[8];
    MPI_Allreduce(my_buf, buf, 8, MPI_DOUBLE, MPI_MIN, mpi_comm);
    if (mpi_rank == 0) {
      cout << " > xp range: ";
      FOR_J3 {
        cout << j << ": " << buf[2*j] << ":" << -buf[2*j+1] << ", ";
      }
      cout << endl;
      cout << " > delta range: " << buf[6] << ":" << -buf[7] << endl;
    }
    
  }
  
  void initFromRestart(const string& filename) {

    // read the mesh points in from the passed restart file, distributed
    // uniformly across ranks...

    using namespace ByteSwap;

    COUT1("Points::initFromRestart("<<filename<<")");

    MPI_File fh;
    {
      char dummy[128];
      sprintf(dummy,"%s",filename.c_str());
      MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
    }

    // first 2 ints are:
    // 0. magic nummber
    // 1. i/o version

    // assume we do not have to byte_swap...
    bool byte_swap = false;

    int itmp[2];
    if (mpi_rank == 0) MPI_File_read(fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
    MPI_Bcast(itmp,2,MPI_INT,0,mpi_comm);

    if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
      ByteSwap::byteSwap(itmp,2);
      if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
	CERR("readRestart: restart file does not start as expected. Exiting.");
      }
      COUT2(" > file requires byte swapping.");
      byte_swap = true;
    }

    assert(UGP_IO_VERSION == 4); // incre
    const int io_version = itmp[1];
    assert((io_version <= 3));
    COUT2(" > io_version " << io_version);

    // initialize the offset...

    MPI_Offset offset = int_size*2;
    int8 * ipora = NULL;

    Header header;
    int done = 0;
    while (done < 2) {

      if (mpi_rank == 0) MPI_File_read_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      MPI_Bcast(&header,1,MPI_Header,0,mpi_comm);
      if (byte_swap) ByteSwap::byteSwapHeader(&header,1);

      switch (header.id) {
      case UGP_IO_NO_FA_CV_COUNTS:
	{
	  assert(np_global == -1);
	  np_global = ByteSwap::getInt8FromLswMswPair(header.idata+4);
	  if (mpi_rank == 0)
	    cout << " > np_global: " << np_global << endl;

	  // use globals to build the striped readers...

	  assert(ipora == NULL);
	  MiscUtils::calcUniformDist(ipora,np_global,mpi_size);
	  assert(ipora[mpi_rank+1]-ipora[mpi_rank] < TWO_BILLION);
	  assert(np == -1);
	  np = int(ipora[mpi_rank+1]-ipora[mpi_rank]);
	}

	break;

      case UGP_IO_CV_CHECK:
	assert(ipora != NULL);
	{
	  int * cv_check = new int[np];
	  MPI_File_read_at_all(fh,offset+header_size+ipora[mpi_rank]*4,
			       cv_check,np,MPI_INT,MPI_STATUS_IGNORE);
	  if (byte_swap) ByteSwap::byteSwap(cv_check,np);
	  for (int icv = 0; icv < np; ++icv)
	    assert(cv_check[icv] == ipora[mpi_rank]+icv);
	  delete[] cv_check;
	}
	if (mpi_rank == 0)
	  cout << " > cv check OK." << endl;
	break;

      case UGP_IO_CV_D2:
	if ( (strcmp(header.name,"X_VV") == 0) || (strcmp(header.name,"X_CV") == 0) ) {
	  assert(ipora != NULL);
	  assert(xp == NULL);
	  xp = new double[np][3];
	  MPI_File_read_at_all(fh,offset+header_size+ipora[mpi_rank]*24,
			       xp,np*3,MPI_DOUBLE,MPI_STATUS_IGNORE);
	  if (byte_swap) ByteSwap::byteSwap(xp,np*3);
	  ++done;
	}
	break;

      case UGP_IO_CV_D1:
	if ( (strcmp(header.name,"VOL_VV") == 0) || (strcmp(header.name,"vol_cv") == 0) ) {
          COUT1(" > setting delta from cell volume");
	  assert(ipora != NULL);
	  assert(delta == NULL);
	  delta = new double[np];
	  MPI_File_read_at_all(fh,offset+header_size+ipora[mpi_rank]*8,
			       delta,np,MPI_DOUBLE,MPI_STATUS_IGNORE);
	  if (byte_swap) ByteSwap::byteSwap(delta,np);

          //estimate delta from vol^1/3
          double power = 1.0/3.0;
          for (int ip = 0; ip < np; ++ip) {
            delta[ip] = pow(delta[ip],power);
          }
	  ++done;
	}
	break;

      case UGP_IO_EOF:
	done = 2;
	break;

      }

      offset += header.skip;

    }

    MPI_File_close(&fh);

    // =======================================================
    // check that the data is as we expect...
    // =======================================================

    assert(np_global >= 0);
    assert(ipora); delete[] ipora;
    if (xp == NULL)
      CERR("restart file " << filename << " does not contain X_CV. add\n\n--REGISTER_CV_R2 X_CV --INIT X_CV=x_cv\n\nto your prepro.in");

    // perterb points randomly in x,y...

    if (Param * param = getParam("RANDOM_XY")) {
      const double delta = param->getDouble();
      COUT1("RANDOM_XY: " << delta);
      for (int ip = 0; ip < np; ++ip) {
	FOR_I2 {
	  xp[ip][i] += 2.0*delta*(double(rand())/double(RAND_MAX)-0.5);
	}
      }
    }

    // and finally set delta...

    double my_buf[6];
    FOR_J6 my_buf[j] = HUGE_VAL; // something large

    for (int ip = 0; ip < np; ++ip) {
      FOR_J3 {
        my_buf[2*j] = min(my_buf[2*j],xp[ip][j]);
        my_buf[2*j+1] = min(my_buf[2*j+1],-xp[ip][j]);
      }
    }

    double buf[6];
    MPI_Allreduce(my_buf, buf, 6, MPI_DOUBLE, MPI_MIN, mpi_comm);
    if (mpi_rank == 0) {
      cout << " > xp range: ";
      FOR_J3 {
        cout << j << ": " << buf[2*j] << ":" << -buf[2*j+1] << ", ";
      }
      cout << endl;
    }

    if (delta == NULL) {
      COUT1(" > using constant delta estimated from max bbox dimension");


      // we choose a small fraction of the maximum length to set delta...

      const double default_delta = 1.0E-3*max( -buf[1]-buf[0], max( -buf[3]-buf[2], -buf[5]-buf[4] ) );
      COUT1(" >   default_delta = " << default_delta);

      delta = new double[np];
      for (int ip = 0; ip < np; ++ip) delta[ip] = default_delta;
    }

  }

  void writeTecplot(const string& filename) const {

    COUT1("Points::writeTecplot(\"" << filename << "\")...");

    FILE * fp;
    if ( mpi_rank == 0 ) {
      fp = fopen(filename.c_str(),"w");
      assert(fp != NULL);
      fprintf(fp,"TITLE = \"flagged faces\"\n");
      fprintf(fp,"VARIABLES = \"X\"\n");
      fprintf(fp,"\"Y\"\n");
      fprintf(fp,"\"Z\"\n");
      fprintf(fp,"\"DELTA\"\n");
      fprintf(fp,"\"RANK\"\n");
    }
    else {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
    }

    //cout << "HACK ******************** writing mpi_rank in delta ****************************** HACK" << endl;

    for (int ip = 0; ip < np; ++ip) {
      fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %d\n",
	      xp[ip][0],
	      xp[ip][1],
	      xp[ip][2],
	      delta[ip]*2.0/3.0,mpi_rank);
    }

    fclose(fp);

    if ( mpi_rank < mpi_size-1 ) {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }

    MPI_Barrier(mpi_comm);

  }

  void writeTecplot(const string& filename,const int * const flag) const {

    COUT1("Points::writeTecplot(\"" << filename << "\")...");

    FILE * fp;
    if ( mpi_rank == 0 ) {
      fp = fopen(filename.c_str(),"w");
      assert(fp != NULL);
      fprintf(fp,"TITLE = \"flagged faces\"\n");
      fprintf(fp,"VARIABLES = \"X\"\n");
      fprintf(fp,"\"Y\"\n");
      fprintf(fp,"\"Z\"\n");
      fprintf(fp,"\"DELTA\"\n");
    }
    else {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
    }

    //cout << "HACK ******************** writing mpi_rank in delta ****************************** HACK" << endl;

    for (int ip = 0; ip < np; ++ip) {
      if (flag[ip]) {
	fprintf(fp,"%18.15le %18.15le %18.15le %18.15le\n",
		xp[ip][0],
		xp[ip][1],
		xp[ip][2],
		delta[ip]*2.0/3.0); //double(mpi_rank));
      }
    }

    fclose(fp);

    if ( mpi_rank < mpi_size-1 ) {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }

    MPI_Barrier(mpi_comm);

  }

  int getClosestPoint(const double xp_check[3]) const {

    double d2_closest = 1.0E+200;
    int ip_closest = -1;
    for (int ip = 0; ip < np; ++ip) {
      const double this_d2 = DIST2(xp[ip],xp_check);
      if (this_d2 < d2_closest) {
	d2_closest = this_d2;
	ip_closest = ip;
      }
    }

    // we need to get the minimum value AND the winning rank: use
    // MPI's MPI_MINLOC capability...

    DoubleInt myDi,di;
    myDi.this_double = d2_closest;
    myDi.this_int = mpi_rank;
    MPI_Allreduce(&myDi,&di,1,MPI_DOUBLE_INT,MPI_MINLOC,mpi_comm);

    if (di.this_int == mpi_rank) {
      assert(ip_closest >= 0);
      return ip_closest;
    }

    // we do not own the closest...
    return -1;

  }

  void makeCart(const double x0[3],const double x1[3],const int n[3]) {

    FOR_I3 assert(n[i] > 0);
    FOR_I3 assert(x1[i] > x0[i]);
    assert(xp == NULL);
    assert(delta == NULL);

    np_global = n[0]*n[1]*n[2];
    assert(np_global > 0);

    // round-robin the points...

    // count first...

    np = 0;
    int round_robin_rank = 0;
    for (int i = 0; i < n[0]; ++i) {
      for (int j = 0; j < n[1]; ++j) {
        for (int k = 0; k < n[2]; ++k) {
          if (mpi_rank == round_robin_rank)
            ++np;
        }
      }
      ++round_robin_rank;
      if (round_robin_rank == mpi_size)
        round_robin_rank = 0;
    }

    // then set...

    //cout << "x0: " << COUT_VEC(x0) << " x1: " << COUT_VEC(x1) << endl;

    xp = new double[np][3];
    delta = new double[np];
    int ip = 0;
    round_robin_rank = 0;
    for (int i = 0; i < n[0]; ++i) {
      for (int j = 0; j < n[1]; ++j) {
        for (int k = 0; k < n[2]; ++k) {
          if (mpi_rank == round_robin_rank) {
            const double wi = (double(i)+0.5)/double(n[0]);
            xp[ip][0] = wi*x1[0] + (1.0-wi)*x0[0];
            const double wj = (double(j)+0.5)/double(n[1]);
            xp[ip][1] = wj*x1[1] + (1.0-wj)*x0[1];
            const double wk = (double(k)+0.5)/double(n[2]);
            xp[ip][2] = wk*x1[2] + (1.0-wk)*x0[2];
            // and for delta, use the diagonal...
            delta[ip] = 1.1*sqrt( (x1[0]-x0[0])*(x1[0]-x0[0])/double(n[0]*n[0]) +
                                  (x1[1]-x0[1])*(x1[1]-x0[1])/double(n[1]*n[1]) +
                                  (x1[2]-x0[2])*(x1[2]-x0[2])/double(n[2]*n[2]) );

            //if ((i==0)&&(j==0)&&(k==0)) {
            //  cout << "min xp[ip]: " << COUT_VEC(xp[ip]) << endl;
            //}
            //else if ((i==n[0]-1)&&(j==n[1]-1)&&(k==n[2]-1)) {
            //  cout << "max xp[ip]: " << COUT_VEC(xp[ip]) << endl;
            //}

            ++ip;
          }
        }
      }
      ++round_robin_rank;
      if (round_robin_rank == mpi_size)
        round_robin_rank = 0;
    }
    assert(ip == np);

  }

  void makeAnnularCyl(const double x0[3],const double x1[3],const double r0,const double r1,const int n[3]) {

    FOR_I3 assert(n[i] > 0);
    assert(r1 > r0);
    assert(xp == NULL);
    assert(delta == NULL);

    np_global = n[0]*n[1]*n[2];
    assert(np_global > 0);

    // round-robin the points...

    // count first...

    np = 0;
    int round_robin_rank = 0;
    for (int i = 0; i < n[0]; ++i) {
      for (int j = 0; j < n[1]; ++j) {
        for (int k = 0; k < n[2]; ++k) {
          if (mpi_rank == round_robin_rank)
            ++np;
        }
      }
      ++round_robin_rank;
      if (round_robin_rank == mpi_size)
        round_robin_rank = 0;
    }

    // then set...

    const double dx[3] = DIFF(x1,x0);
    const double L = MAG(dx);
    double e1[3],e2[3];
    MiscUtils::getBestE1E2FromE0(e1,e2,dx);

    xp = new double[np][3];
    delta = new double[np];
    int ip = 0;
    round_robin_rank = 0;
    for (int i = 0; i < n[0]; ++i) {
      for (int j = 0; j < n[1]; ++j) {
        for (int k = 0; k < n[2]; ++k) {
          if (mpi_rank == round_robin_rank) {
            const double wi = (double(i)+0.5)/double(n[0]);
            const double dl = L/double(n[0]);
            const double wj = (double(j)+0.5)/double(n[1]);
            const double r = r1*wj + r0*(1.0-wj);
            const double dr = (r1-r0)/double(n[1]);
            const double wk = (double(k)+0.5)/double(n[2]);
            FOR_I3 xp[ip][i] = x0[i] + wi*dx[i] + r*cos(2.0*M_PI*wk)*e1[i] + r*sin(2.0*M_PI*wk)*e2[i];
            // and for delta, use the diagonal...
            delta[ip] = 1.1*sqrt(dl*dl + dr*dr + 4.0*M_PI*M_PI*r*r/double(n[2]*n[2]));
            ++ip;
          }
        }
      }
      ++round_robin_rank;
      if (round_robin_rank == mpi_size)
        round_robin_rank = 0;
    }
    assert(ip == np);

  }

  void makeHexcore(const double x0[3],const double x1[3],const double r,const int n[3]) {

    FOR_I3 assert(n[i] > 0);
    assert(r > 0.0);
    assert(xp == NULL);
    assert(delta == NULL);

    int ni = n[0];
    int nr = n[1];
    int nt = n[2];

    // ensure integer counts per direction are compatible
    if (nr%2 != 0) {
      nr++;
      COUT1(" > radial spacing must be a multiple of 2; incrementing 'nr' to " << nr);
    }
    if (nt%4 != 0) {
      nt += 4 - (nt%4);
      COUT1(" > theta spacing must be a multiple of 4; incrementing 'nt' to " << nt);
    }
    np_global = int8(ni)*( int8(nr*nt/2) + int8(nt*nt/16) );
    COUT1(" > generating " << np_global << " points");
    assert(np_global > 0);


    // then set...
    const double e0[3] = DIFF(x1,x0);
    const double L = MAG(e0);
    double e1[3],e2[3];
    MiscUtils::getBestE1E2FromE0(e1,e2,e0);
    COUT2(" > e1,e2 vectors: " << COUT_VEC(e1) << " " << COUT_VEC(e2));

    // assign points
    int8 * ipora = NULL;
    buildUniformXora(ipora,np_global);
    assert( ipora[mpi_rank+1]-ipora[mpi_rank] < TWO_BILLION );
    np = int( ipora[mpi_rank+1]-ipora[mpi_rank] );

    xp = new double[np][3];
    delta = new double[np];

    int8 ip_global = 0;
    for (int ii = 0; ii < ni; ++ii) {
      double x[3];
      FOR_I3 x[i] = x0[i] + ((double(ii)+0.5)/double(ni))*e0[i];
      // 4 near-wall sectors...
      for (int t0 = 0; t0 < 4; ++t0) {
        const double theta0 = (double(t0)+0.5)/2.0*M_PI;
        const double theta1 = (double(t0)+1.5)/2.0*M_PI;
        double xws[3];
        FOR_I3 xws[i] = 0.5*r*(e1[i]*cos(theta1) + e2[i]*sin(theta1));
        double xwn[3];
        FOR_I3 xwn[i] = 1.0*r*(e1[i]*cos(theta1) + e2[i]*sin(theta1));
        double xes[3];
        FOR_I3 xes[i] = 0.5*r*(e1[i]*cos(theta0) + e2[i]*sin(theta0));
        double xen[3];
        FOR_I3 xen[i] = 1.0*r*(e1[i]*cos(theta0) + e2[i]*sin(theta0));
        for (int j = 0; j < nt/4; ++j) {
          for (int k = 0; k < nr/2; ++k) {
            if ( (ipora[mpi_rank] <= ip_global)&&(ip_global < ipora[mpi_rank+1]) ) {
              const int ip = ip_global - ipora[mpi_rank];
              assert((ip >= 0)&&(ip < np));

              const double a = (double(j)+0.5)/double(nt/4);
              const double b = (double(k)+0.5)/double(nr/2);
              const double theta = a*theta1 + (1.0-a)*theta0;
              double xw[3]; FOR_I3 xw[i] = spacingFunction(b,1.0,xws[i],xwn[i]);
              double xe[3]; FOR_I3 xe[i] = spacingFunction(b,1.0,xes[i],xen[i]);
              double xs[3];
              FOR_I3 xs[i] = 0.25*r*cos(theta)*e1[i] + 0.25*r*sin(theta)*e2[i] + (0.5*a*xws[i] + 0.5*(1.0-a)*xes[i])*1.0;
              double xn[3];
              FOR_I3 xn[i] = r*cos(theta)*e1[i] + r*sin(theta)*e2[i];
              FOR_I3 xp[ip][i] = x[i] + (1.0-b)*xs[i] + b*xn[i] + (1.0-a)*xw[i] + a*xe[i] - ( a*b*xen[i] + a*(1.0-b)*xes[i] + (1.0-a)*b*xwn[i] + (1.0-a)*(1.0-b)*xws[i] );

              // now perturb in each direction to get deltav...
              double d_max = L/double(ni);
              {
                const double a = (double(j)    )/double(nt/4);
                const double b = (double(k)+0.5)/double(nr/2);
                const double theta = a*theta1 + (1.0-a)*theta0;
                double xw[3]; FOR_I3 xw[i] = spacingFunction(b,1.0,xws[i],xwn[i]);
                double xe[3]; FOR_I3 xe[i] = spacingFunction(b,1.0,xes[i],xen[i]);
                double xs[3];
                FOR_I3 xs[i] = 0.25*r*cos(theta)*e1[i] + 0.25*r*sin(theta)*e2[i] + (0.5*a*xws[i] + 0.5*(1.0-a)*xes[i])*1.0;
                double xn[3];
                FOR_I3 xn[i] = r*cos(theta)*e1[i] + r*sin(theta)*e2[i];
                double xdeltav[3]; FOR_I3 xdeltav[i] = x[i] + (1.0-b)*xs[i] + b*xn[i] + (1.0-a)*xw[i] + a*xe[i] - ( a*b*xen[i] + a*(1.0-b)*xes[i] + (1.0-a)*b*xwn[i] + (1.0-a)*(1.0-b)*xws[i] );
                d_max = max(d_max,2.0*DIST(xdeltav,xp[ip]));
              }
              {
                const double a = (double(j)+1.0)/double(nt/4);
                const double b = (double(k)+0.5)/double(nr/2);
                const double theta = a*theta1 + (1.0-a)*theta0;
                double xw[3]; FOR_I3 xw[i] = spacingFunction(b,1.0,xws[i],xwn[i]);
                double xe[3]; FOR_I3 xe[i] = spacingFunction(b,1.0,xes[i],xen[i]);
                double xs[3];
                FOR_I3 xs[i] = 0.25*r*cos(theta)*e1[i] + 0.25*r*sin(theta)*e2[i] + (0.5*a*xws[i] + 0.5*(1.0-a)*xes[i])*1.0;
                double xn[3];
                FOR_I3 xn[i] = r*cos(theta)*e1[i] + r*sin(theta)*e2[i];
                double xdeltav[3]; FOR_I3 xdeltav[i] = x[i] + (1.0-b)*xs[i] + b*xn[i] + (1.0-a)*xw[i] + a*xe[i] - ( a*b*xen[i] + a*(1.0-b)*xes[i] + (1.0-a)*b*xwn[i] + (1.0-a)*(1.0-b)*xws[i] );
                d_max = max(d_max,2.0*DIST(xdeltav,xp[ip]));
              }
              {
                const double a = (double(j)+0.5)/double(nt/4);
                const double b = (double(k)    )/double(nr/2);
                const double theta = a*theta1 + (1.0-a)*theta0;
                double xw[3]; FOR_I3 xw[i] = spacingFunction(b,1.0,xws[i],xwn[i]);
                double xe[3]; FOR_I3 xe[i] = spacingFunction(b,1.0,xes[i],xen[i]);
                double xs[3];
                FOR_I3 xs[i] = 0.25*r*cos(theta)*e1[i] + 0.25*r*sin(theta)*e2[i] + (0.5*a*xws[i] + 0.5*(1.0-a)*xes[i])*1.0;
                double xn[3];
                FOR_I3 xn[i] = r*cos(theta)*e1[i] + r*sin(theta)*e2[i];
                double xdeltav[3]; FOR_I3 xdeltav[i] = x[i] + (1.0-b)*xs[i] + b*xn[i] + (1.0-a)*xw[i] + a*xe[i] - ( a*b*xen[i] + a*(1.0-b)*xes[i] + (1.0-a)*b*xwn[i] + (1.0-a)*(1.0-b)*xws[i] );
                d_max = max(d_max,2.0*DIST(xdeltav,xp[ip]));
              }
              {
                double xw[3]; FOR_I3 xw[i] = spacingFunction(b,1.0,xws[i],xwn[i]);
                double xe[3]; FOR_I3 xe[i] = spacingFunction(b,1.0,xes[i],xen[i]);
                double xs[3];
                FOR_I3 xs[i] = 0.25*r*cos(theta)*e1[i] + 0.25*r*sin(theta)*e2[i] + (0.5*a*xws[i] + 0.5*(1.0-a)*xes[i])*1.0;
                double xn[3];
                FOR_I3 xn[i] = r*cos(theta)*e1[i] + r*sin(theta)*e2[i];
                double xdeltav[3]; FOR_I3 xdeltav[i] = x[i] + (1.0-b)*xs[i] + b*xn[i] + (1.0-a)*xw[i] + a*xe[i] - ( a*b*xen[i] + a*(1.0-b)*xes[i] + (1.0-a)*b*xwn[i] + (1.0-a)*(1.0-b)*xws[i] );
                d_max = max(d_max,2.0*DIST(xdeltav,xp[ip]));
              }
              delta[ip] = 1.1*d_max;

              //cout << "XXXXXX: " << pts->xp[ip][0] << " " << pts->xp[ip][1] << " " << pts->xp[ip][2] << endl;

            }
            ++ip_global;
          }
        }
      }
      // and the center sector...
      double xws[3];
      FOR_I3 xws[i] = -0.5*r/sqrt(2.0)*e1[i] -0.5*r/sqrt(2.0)*e2[i];
      double xwn[3];
      FOR_I3 xwn[i] = -0.5*r/sqrt(2.0)*e1[i] +0.5*r/sqrt(2.0)*e2[i];
      double xes[3];
      FOR_I3 xes[i] = 0.5*r/sqrt(2.0)*e1[i] -0.5*r/sqrt(2.0)*e2[i];
      double xen[3];
      FOR_I3 xen[i] = 0.5*r/sqrt(2.0)*e1[i] +0.5*r/sqrt(2.0)*e2[i];
      for (int j = 0; j < nt/4; ++j) {
        for (int k = 0; k < nt/4; ++k) {
          if ( (ipora[mpi_rank] <= ip_global)&&(ip_global < ipora[mpi_rank+1]) ) {

            const int ip = ip_global - ipora[mpi_rank];
            assert((ip >= 0)&&(ip < np));

            const double a = (double(j)+0.5)/double(nt/4);
            const double b = (double(k)+0.5)/double(nt/4);
            const double theta_a = a*0.5/2.0*M_PI + (1.0-a)*1.5/2.0*M_PI;
            double xs[3];
            FOR_I3 xs[i] = 0.25*r*cos(theta_a)*e1[i] -0.25*r*sin(theta_a)*e2[i] + (0.5*a*xes[i] + 0.5*(1.0-a)*xws[i])*1.0;
            double xn[3];
            FOR_I3 xn[i] = 0.25*r*cos(theta_a)*e1[i] + 0.25*r*sin(theta_a)*e2[i] + (0.5*a*xen[i] + 0.5*(1.0-a)*xwn[i])*1.0;
            const double theta_b = b*1.5/2.0*M_PI + (1.0-b)*2.5/2.0*M_PI;
            double xw[3];
            FOR_I3 xw[i] = 0.25*r*cos(theta_b)*e1[i] + 0.25*r*sin(theta_b)*e2[i] + (0.5*b*xwn[i] + 0.5*(1.0-b)*xws[i])*1.0;
            double xe[3];
            FOR_I3 xe[i] = -0.25*r*cos(theta_b)*e1[i] + 0.25*r*sin(theta_b)*e2[i] + (0.5*b*xen[i] + 0.5*(1.0-b)*xes[i])*1.0;
            FOR_I3 xp[ip][i] = x[i] + (1.0-b)*xs[i] + b*xn[i] + (1.0-a)*xw[i] + a*xe[i] - ( a*b*xen[i] + a*(1.0-b)*xes[i] + (1.0-a)*b*xwn[i] + (1.0-a)*(1.0-b)*xws[i] );

            // now perturb in each direction to get deltav...
            double d_max = L/double(ni);
            {
              const double a = (double(j)    )/double(nt/4);
              const double b = (double(k)+0.5)/double(nt/4);
              const double theta_a = a*0.5/2.0*M_PI + (1.0-a)*1.5/2.0*M_PI;
              double xs[3];
              FOR_I3 xs[i] = 0.25*r*cos(theta_a)*e1[i] -0.25*r*sin(theta_a)*e2[i] + (0.5*a*xes[i] + 0.5*(1.0-a)*xws[i])*1.0;
              double xn[3];
              FOR_I3 xn[i] = 0.25*r*cos(theta_a)*e1[i] + 0.25*r*sin(theta_a)*e2[i] + (0.5*a*xen[i] + 0.5*(1.0-a)*xwn[i])*1.0;
              const double theta_b = b*1.5/2.0*M_PI + (1.0-b)*2.5/2.0*M_PI;
              double xw[3];
              FOR_I3 xw[i] = 0.25*r*cos(theta_b)*e1[i] + 0.25*r*sin(theta_b)*e2[i] + (0.5*b*xwn[i] + 0.5*(1.0-b)*xws[i])*1.0;
              double xe[3];
              FOR_I3 xe[i] = -0.25*r*cos(theta_b)*e1[i] + 0.25*r*sin(theta_b)*e2[i] + (0.5*b*xen[i] + 0.5*(1.0-b)*xes[i])*1.0;
              double xdeltav[3]; FOR_I3 xdeltav[i] = x[i] + (1.0-b)*xs[i] + b*xn[i] + (1.0-a)*xw[i] + a*xe[i] - ( a*b*xen[i] + a*(1.0-b)*xes[i] + (1.0-a)*b*xwn[i] + (1.0-a)*(1.0-b)*xws[i] );
              d_max = max(d_max,2.0*DIST(xdeltav,xp[ip]));
            }
            {
              const double a = (double(j)+1.0)/double(nt/4);
              const double b = (double(k)+0.5)/double(nt/4);
              const double theta_a = a*0.5/2.0*M_PI + (1.0-a)*1.5/2.0*M_PI;
              double xs[3];
              FOR_I3 xs[i] = 0.25*r*cos(theta_a)*e1[i] -0.25*r*sin(theta_a)*e2[i] + (0.5*a*xes[i] + 0.5*(1.0-a)*xws[i])*1.0;
              double xn[3];
              FOR_I3 xn[i] = 0.25*r*cos(theta_a)*e1[i] + 0.25*r*sin(theta_a)*e2[i] + (0.5*a*xen[i] + 0.5*(1.0-a)*xwn[i])*1.0;
              const double theta_b = b*1.5/2.0*M_PI + (1.0-b)*2.5/2.0*M_PI;
              double xw[3];
              FOR_I3 xw[i] = 0.25*r*cos(theta_b)*e1[i] + 0.25*r*sin(theta_b)*e2[i] + (0.5*b*xwn[i] + 0.5*(1.0-b)*xws[i])*1.0;
              double xe[3];
              FOR_I3 xe[i] = -0.25*r*cos(theta_b)*e1[i] + 0.25*r*sin(theta_b)*e2[i] + (0.5*b*xen[i] + 0.5*(1.0-b)*xes[i])*1.0;
              double xdeltav[3]; FOR_I3 xdeltav[i] = x[i] + (1.0-b)*xs[i] + b*xn[i] + (1.0-a)*xw[i] + a*xe[i] - ( a*b*xen[i] + a*(1.0-b)*xes[i] + (1.0-a)*b*xwn[i] + (1.0-a)*(1.0-b)*xws[i] );
              d_max = max(d_max,2.0*DIST(xdeltav,xp[ip]));
            }
            {
              const double a = (double(j)+0.5)/double(nt/4);
              const double b = (double(k)    )/double(nt/4);
              const double theta_a = a*0.5/2.0*M_PI + (1.0-a)*1.5/2.0*M_PI;
              double xs[3];
              FOR_I3 xs[i] = 0.25*r*cos(theta_a)*e1[i] -0.25*r*sin(theta_a)*e2[i] + (0.5*a*xes[i] + 0.5*(1.0-a)*xws[i])*1.0;
              double xn[3];
              FOR_I3 xn[i] = 0.25*r*cos(theta_a)*e1[i] + 0.25*r*sin(theta_a)*e2[i] + (0.5*a*xen[i] + 0.5*(1.0-a)*xwn[i])*1.0;
              const double theta_b = b*1.5/2.0*M_PI + (1.0-b)*2.5/2.0*M_PI;
              double xw[3];
              FOR_I3 xw[i] = 0.25*r*cos(theta_b)*e1[i] + 0.25*r*sin(theta_b)*e2[i] + (0.5*b*xwn[i] + 0.5*(1.0-b)*xws[i])*1.0;
              double xe[3];
              FOR_I3 xe[i] = -0.25*r*cos(theta_b)*e1[i] + 0.25*r*sin(theta_b)*e2[i] + (0.5*b*xen[i] + 0.5*(1.0-b)*xes[i])*1.0;
              double xdeltav[3]; FOR_I3 xdeltav[i] = x[i] + (1.0-b)*xs[i] + b*xn[i] + (1.0-a)*xw[i] + a*xe[i] - ( a*b*xen[i] + a*(1.0-b)*xes[i] + (1.0-a)*b*xwn[i] + (1.0-a)*(1.0-b)*xws[i] );
              d_max = max(d_max,2.0*DIST(xdeltav,xp[ip]));
            }
            {
              const double a = (double(j)+0.5)/double(nt/4);
              const double b = (double(k)+1.0)/double(nt/4);
              const double theta_a = a*0.5/2.0*M_PI + (1.0-a)*1.5/2.0*M_PI;
              double xs[3];
              FOR_I3 xs[i] = 0.25*r*cos(theta_a)*e1[i] -0.25*r*sin(theta_a)*e2[i] + (0.5*a*xes[i] + 0.5*(1.0-a)*xws[i])*1.0;
              double xn[3];
              FOR_I3 xn[i] = 0.25*r*cos(theta_a)*e1[i] + 0.25*r*sin(theta_a)*e2[i] + (0.5*a*xen[i] + 0.5*(1.0-a)*xwn[i])*1.0;
              const double theta_b = b*1.5/2.0*M_PI + (1.0-b)*2.5/2.0*M_PI;
              double xw[3];
              FOR_I3 xw[i] = 0.25*r*cos(theta_b)*e1[i] + 0.25*r*sin(theta_b)*e2[i] + (0.5*b*xwn[i] + 0.5*(1.0-b)*xws[i])*1.0;
              double xe[3];
              FOR_I3 xe[i] = -0.25*r*cos(theta_b)*e1[i] + 0.25*r*sin(theta_b)*e2[i] + (0.5*b*xen[i] + 0.5*(1.0-b)*xes[i])*1.0;
              double xdeltav[3]; FOR_I3 xdeltav[i] = x[i] + (1.0-b)*xs[i] + b*xn[i] + (1.0-a)*xw[i] + a*xe[i] - ( a*b*xen[i] + a*(1.0-b)*xes[i] + (1.0-a)*b*xwn[i] + (1.0-a)*(1.0-b)*xws[i] );
              d_max = max(d_max,2.0*DIST(xdeltav,xp[ip]));
            }
            delta[ip] = 1.1*d_max;
            //cout << "XXXXXX: " << pts->xp[ip][0] << " " << pts->xp[ip][1] << " " << pts->xp[ip][2] << endl;
          }
          ++ip_global;
        }
      }
    }

    delete[] ipora;
    assert( ip_global == np_global );
  }

  void makeDlrt(const double x0[3],const double x1[3], double r0[2], double r1[2],const double d[3]) {

    // ============================================================================
    // Annulus/cyl mesh points...
    //
    // dx = axis aligned resolution
    // dr = radius resolution
    // dt = theta resolution (arc-length of circle spacing)
    //
    // ============================================================================

    FOR_I3 assert(d[i] > 0);
    assert(r0[1] > 0.0);
    assert(r1[1] > 0.0);
    assert(xp == NULL);
    assert(delta == NULL);

    // then set...
    const double e0[3] = DIFF(x1,x0);
    const double L = MAG(e0);
    double e1[3],e2[3];
    MiscUtils::getBestE1E2FromE0(e1,e2,e0);
    COUT2(" > e1,e2 vectors: " << COUT_VEC(e1) << " " << COUT_VEC(e2));

    const double SMALL = 1E-12;
    if (r0[0] < SMALL) r0[0] = 0.0;
    if (r1[0] < SMALL) r1[0] = 0.0;

    // determine total number of points to generate
    const int nx = max(1, int(ceil(L/d[0])) );
    assert(nx > 0);
    COUT2(" > points along axis: " << nx);
    int (*npoxi)[2] = new int[nx][2];  // at each x-index 0: nr, 1:nt
    int nt_max = 0;  // store largest nt for pointer max allocation later

    // first pass to count
    np_global = 0;
    for (int xi = 0; xi < nx; ++xi) {
      // min/max radius at this dx
      const double factor = (double(xi) + 0.5)/double(nx);
      const double ri_min = r0[0] + factor*(r1[0]-r0[0]);
      const double ri_max = r0[1] + factor*(r1[1]-r0[1]);

      npoxi[xi][0] = max(1,int(ceil(fabs(ri_max-ri_min) / d[1])));
      assert(npoxi[xi][0] > 0);

      // count points per radius level
      npoxi[xi][1] = 0;
      for (int ir = 0, nr = npoxi[xi][0]; ir < nr; ++ir) {
        const double this_r = spacingFunction((double(ir)+0.5)/double(nr),1.0,ri_min,ri_max);
        const int ntheta = max(3,int(ceil(2.0*M_PI*this_r / d[2])));
        assert(ntheta > 0);
        nt_max = max(nt_max,ntheta);
        npoxi[xi][1] += ntheta;
      }
      np_global += int8(npoxi[xi][1]);
    }

    // allocate for parallel points build
    int8 * ipora = NULL;
    buildUniformXora(ipora,np_global);
    assert( ipora[mpi_rank+1]-ipora[mpi_rank] < TWO_BILLION );
    np = int( ipora[mpi_rank+1]-ipora[mpi_rank] );

    assert(xp == NULL);    xp = new double[np][3];
    assert(delta == NULL); delta = new double[np];

    int8 ip_global = 0;
    double (*xp_temp)[3] = new double[nt_max][3];

    for (int xi = 0; xi < nx; ++xi) {
      // min/max radius at this dx
      const double factor = (double(xi) + 0.5)/double(nx);
      const double ri_min = r0[0] + factor*(r1[0]-r0[0]);
      const double ri_max = r0[1] + factor*(r1[1]-r0[1]);

      double this_xc[3];
      FOR_J3 {
        this_xc[j] = x0[j] + factor*(e0[j]);
      }

      for (int ir = 0, nr = npoxi[xi][0]; ir < nr; ++ir) {
        const double this_r = spacingFunction((double(ir)+0.5)/double(nr),1.0,ri_min,ri_max);

        // populate the entire ring, but only pull in the points you own
        const int nt = max(3,int(ceil(2.0*M_PI*this_r / d[2])));

        const bool stagger = (ir%2 == 1);

        MiscUtils::createCirclePts(xp_temp,0,this_xc,e0,this_r,nt,stagger);  // stagger starting point if odd nt

        for (int ti = 0; ti < nt; ++ti) {
          if ( (ipora[mpi_rank] <= ip_global)&&(ip_global < ipora[mpi_rank+1]) ) {
            const int ip = ip_global - ipora[mpi_rank];

            FOR_I3 xp[ip][i] = xp_temp[ti][i];

            // set delta based on max of x,r,t spacing
            double dmax = max(d[0],max(d[2],d[1]));
            delta[ip] = 1.1 * dmax;
          }
          ++ip_global;
        }//ti
      }//ri
    }//xi

    // cleanup
    DELETE(xp_temp);
    DELETE(npoxi);

    delete[] ipora;
    assert( ip_global == np_global );

  } //initDxDrDtPoints

};

#endif
