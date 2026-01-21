#include "CtiRegister.hpp"
#include "ByteSwap.hpp"
#include "MiscUtils.hpp"
#include <map>
#include "RestartHashUtilities.hpp"

namespace CtiRegister {

  bool cti_eval_verbosity = true;

  void setEvalVerbosity(const bool b) {
    cti_eval_verbosity = b;
  }

  CtiData::CtiData() {
    //cout << "CtiData(): " << this << endl;
    reset();
  }

  CtiData::CtiData(const CtiData& data) {
    //cout << "CtiData(): " << this << " copy of: " << &data << " " << data << endl;
    copy(data);
  }

  CtiData::~CtiData() {
    //cout << "~CtiData(): " << this << " mem_flag: " << mem_flag << endl;
    clear();
  }

  void CtiData::clear() {
    if (mem_flag) {
      const int dt = getType();
      switch (dt) {
      case S_DATA:
        assert(size_ptr == NULL);
        delete (string*)data_ptr; data_ptr = NULL;
        break;
      case I_DATA:
        assert(size_ptr == NULL);
        delete (int*)data_ptr; data_ptr = NULL;
	break;
      case D_DATA:
        assert(size_ptr == NULL);
        delete (double*)data_ptr; data_ptr = NULL;
	break;
      case D3_DATA:
        assert(size_ptr == NULL);
        delete[] (double*)data_ptr; data_ptr = NULL;
	break;
      case DN_DATA:
        assert( size_ptr != NULL);
        //delete size_ptr; // now just points directly to &n
        size_ptr = NULL;
        // this check because of delayed allocation. Get rid of it when STATS is changed
        if (data_ptr != NULL) {
          delete[] (double*)data_ptr;
          data_ptr = NULL;
        }
	break;
      case DN3_DATA:
        assert( size_ptr != NULL);
        //delete size_ptr;
        size_ptr = NULL;
        // this check because of delayed allocation. Get rid of it when STATS is changed
        if (data_ptr != NULL) {
          delete[] (double(*)[3])data_ptr;
          data_ptr = NULL;
        }
	break;
      case DN33_DATA:
        assert( size_ptr != NULL);
        //delete size_ptr;
        size_ptr = NULL;
        // this check because of delayed allocation. Get rid of it when STATS is changed
        if (data_ptr != NULL) {
          delete[] (double(*)[3][3])data_ptr;
          data_ptr = NULL;
        }
	break;
      case IN_DATA:
        assert( size_ptr != NULL);
        assert( data_ptr != NULL);
        //delete size_ptr;
        size_ptr = NULL;
	delete[] (int*)data_ptr; data_ptr = NULL;
	break;
      default:
	assert(size_ptr == NULL);
	assert(data_ptr == NULL);
      }
    }
  }

  bool CtiData::empty() const {
    if ((data_ptr == NULL) &&
	(size_ptr == NULL) &&
        (mem_flag == false) &&
        (access_flag == false) &&
	(bits == 0) &&
	(stride == 0) &&
	(offset == 0) &&
	(flag == 0) &&
	(xora_ptr == NULL) &&
	(dde_ptr == NULL))
      return true;
    else
      return false;
  }

  void CtiData::reset() {
    data_ptr = NULL;
    size_ptr = NULL;
    mem_flag = false;
    access_flag = false;
    bits     = 0;
    stride   = 0;
    flag     = 0;
    offset   = 0;
    xora_ptr = NULL;
    dde_ptr  = NULL;
  }

  void CtiData::copy(const CtiData& data) {
    data_ptr = data.data_ptr;
    size_ptr = data.size_ptr;
    mem_flag = data.mem_flag;
    access_flag = data.access_flag;
    bits     = data.bits;
    stride   = data.stride;
    offset   = data.offset;
    flag     = data.flag;
    xora_ptr = data.xora_ptr;
    dde_ptr  = data.dde_ptr;
  }

  void CtiData::setType(const int type) {
    bits |= type;
  }

  int CtiData::getType() const {
    // for now, the type is stored as a bit. This is a bit expensive,
    // because types are mutually exclusive. Fix later...
    if (bits & D_DATA) {
      assert(!(bits & (D3_DATA|DN_DATA|DN3_DATA|DN33_DATA|I_DATA|IN_DATA)));
      return D_DATA;
    }
    else if (bits & D3_DATA) {
      assert(!(bits & (D_DATA|DN_DATA|DN3_DATA|DN33_DATA|I_DATA|IN_DATA)));
      return D3_DATA;
    }
    else if (bits & DN_DATA) {
      assert(!(bits & (D_DATA|D3_DATA|DN3_DATA|DN33_DATA|I_DATA|IN_DATA)));
      return DN_DATA;
    }
    else if (bits & DN3_DATA) {
      assert(!(bits & (D_DATA|D3_DATA|DN_DATA|DN33_DATA|I_DATA|IN_DATA)));
      return DN3_DATA;
    }
    else if (bits & DN33_DATA) {
      assert(!(bits & (D_DATA|D3_DATA|DN_DATA|DN3_DATA|I_DATA|IN_DATA)));
      return DN33_DATA;
    }
    else if (bits & I_DATA) {
      assert(!(bits & (D_DATA|D3_DATA|DN_DATA|DN3_DATA|DN33_DATA|IN_DATA)));
      return I_DATA;
    }
    else if (bits & IN_DATA) {
      assert(!(bits & (D_DATA|D3_DATA|DN_DATA|DN3_DATA|DN33_DATA|I_DATA)));
      return IN_DATA;
    }
    else if (bits & S_DATA) {
      assert(!(bits & (D_DATA|D3_DATA|DN_DATA|DN3_DATA|DN33_DATA|IN_DATA|I_DATA)));
      return S_DATA;
    }
    else {
      return 0;
    }
  }

  string CtiData::getTypeAsString() const {
    if (bits & D_DATA) {
      assert(!(bits & (D3_DATA|DN_DATA|DN3_DATA|DN33_DATA|I_DATA|IN_DATA)));
      return "D_DATA";
    }
    else if (bits & D3_DATA) {
      assert(!(bits & (D_DATA|DN_DATA|DN3_DATA|DN33_DATA|I_DATA|IN_DATA)));
      return "D3_DATA";
    }
    else if (bits & DN_DATA) {
      assert(!(bits & (D_DATA|D3_DATA|DN3_DATA|DN33_DATA|I_DATA|IN_DATA)));
      return "DN_DATA";
    }
    else if (bits & DN3_DATA) {
      assert(!(bits & (D_DATA|D3_DATA|DN_DATA|DN33_DATA|I_DATA|IN_DATA)));
      return "DN3_DATA";
    }
    else if (bits & DN33_DATA) {
      assert(!(bits & (D_DATA|D3_DATA|DN_DATA|DN3_DATA|I_DATA|IN_DATA)));
      return "DN33_DATA";
    }
    else if (bits & I_DATA) {
      assert(!(bits & (D_DATA|D3_DATA|DN_DATA|DN3_DATA|DN33_DATA|IN_DATA)));
      return "I_DATA";
    }
    else if (bits & IN_DATA) {
      assert(!(bits & (D_DATA|D3_DATA|DN_DATA|DN3_DATA|DN33_DATA|I_DATA)));
      return "IN_DATA";
    }
    else if (bits & S_DATA) {
      assert(!(bits & (D_DATA|D3_DATA|DN_DATA|DN3_DATA|DN33_DATA|IN_DATA|I_DATA)));
      return "S_DATA";
    }
    else {
      return "UNKNOWN";
    }
  }

  // the topology is defined by one of the data types
  // like cv, bf, (lp eventually), AND also an index.
  // this prevents doing operations on data that are
  // not part of the same topology.

  void CtiData::setTopology(const int topo) {
    // insist on valid topo...
    assert(topo & (CV_DATA|BF_DATA|FA_DATA|EF_DATA|NO_DATA|LP_DATA));
    bits |= topo; // could include index
  }


  void CtiData::setTopology(const int unindexed_topo,const int index) {
    // topo here needs to be unindexed
    // index can be particle class, bf zone index, or fa sign boolean (for now)
    assert((unindexed_topo == BF_DATA)||(unindexed_topo == LP_DATA)||(unindexed_topo == FA_DATA)||(unindexed_topo == NO_DATA));
    if (unindexed_topo == FA_DATA)
      assert((index == 0)||(index == 1));
    bits |= unindexed_topo;
    bits |= (index<<INDEX_SHIFT);
  }

  int CtiData::getTopology() const {
    // for now, the topology is stored as a bit. This is a bit expensive,
    // because, like types, topologies are mutually exclusive. Fix later...
    if (bits & CV_DATA) {
      assert(!(bits & (BF_DATA|FA_DATA|EF_DATA|NO_DATA|LP_DATA) ));
      assert((bits>>INDEX_SHIFT) == 0);
      return CV_DATA;
    }
    else if (bits & BF_DATA) {
      assert(!(bits & (CV_DATA|FA_DATA|EF_DATA|NO_DATA|LP_DATA) ));
      return bits & TOPOLOGY_MASK;
    }
    else if (bits & FA_DATA) {
      assert(!(bits & (CV_DATA|BF_DATA|EF_DATA|NO_DATA|LP_DATA) ));
      assert(((bits>>INDEX_SHIFT) == 0)||((bits>>INDEX_SHIFT) == 1));
      return bits & TOPOLOGY_MASK;
    }
    else if (bits & EF_DATA) {
      assert(!(bits & (CV_DATA|BF_DATA|FA_DATA|NO_DATA|LP_DATA) ));
      assert((bits>>INDEX_SHIFT) == 0);
      return EF_DATA;
    }
    else if (bits & NO_DATA) {
      assert(!(bits & (CV_DATA|BF_DATA|FA_DATA|EF_DATA|LP_DATA) ));
      return bits & TOPOLOGY_MASK;
    }
    else if (bits & LP_DATA) {
      assert(!(bits & (CV_DATA|BF_DATA|FA_DATA|EF_DATA|NO_DATA) ));
      return bits & TOPOLOGY_MASK;
    }
    else {
      return 0; // may be value data?
    }
  }

  string CtiData::getTopologyAsString() const {
    // for now, the topology is stored as a bit. This is a bit expensive,
    // because, like types, topologies are mutually exclusive. Fix later...
    if (bits & CV_DATA) {
      assert(!(bits & (BF_DATA|FA_DATA|EF_DATA|NO_DATA|LP_DATA) ));
      assert((bits>>INDEX_SHIFT) == 0);
      return "CV_DATA";
    }
    else if (bits & BF_DATA) {
      assert(!(bits & (CV_DATA|FA_DATA|EF_DATA|NO_DATA|LP_DATA) ));
      stringstream ss;
      ss << "BF_DATA " << (bits>>INDEX_SHIFT);
      return ss.str();
    }
    else if (bits & FA_DATA) {
      assert(!(bits & (CV_DATA|BF_DATA|EF_DATA|NO_DATA|LP_DATA) ));
      assert(((bits>>INDEX_SHIFT) == 0)||((bits>>INDEX_SHIFT) == 1));
      if ((bits>>INDEX_SHIFT) == 0)
        return "UNSIGNED_FA_DATA";
      else
        return "SIGNED_FA_DATA";
    }
    else if (bits & EF_DATA) {
      assert(!(bits & (CV_DATA|BF_DATA|FA_DATA|NO_DATA|LP_DATA) ));
      assert((bits>>INDEX_SHIFT) == 0);
      return "EF_DATA";
    }
    else if (bits & NO_DATA) {
      assert(!(bits & (CV_DATA|BF_DATA|FA_DATA|EF_DATA|LP_DATA) ));
      stringstream ss;
      ss << "NO_DATA " << (bits>>INDEX_SHIFT);
      return ss.str();
    }
    else if (bits & LP_DATA) {
      assert(!(bits & (CV_DATA|BF_DATA|FA_DATA|EF_DATA|NO_DATA) ));
      stringstream ss;
      ss << "LP_DATA " << (bits>>INDEX_SHIFT);
      return ss.str();
    }
    else {
      assert((bits>>INDEX_SHIFT) == 0);
      return "VALUE_DATA"; // may be value data?
    }
  }

  int CtiData::getUnindexedTopology() const {
    // for now, the topology is stored as a bit. This is a bit expensive,
    // because, like types, topologies are mutually exclusive. Fix later...
    if (bits & CV_DATA) {
      assert(!(bits & (BF_DATA|FA_DATA|EF_DATA|NO_DATA|LP_DATA) ));
      assert((bits>>INDEX_SHIFT) == 0);
      return CV_DATA;
    }
    else if (bits & BF_DATA) {
      assert(!(bits & (CV_DATA|FA_DATA|EF_DATA|NO_DATA|LP_DATA) ));
      return BF_DATA;
    }
    else if (bits & FA_DATA) {
      assert(!(bits & (CV_DATA|BF_DATA|EF_DATA|NO_DATA|LP_DATA) ));
      assert(((bits>>INDEX_SHIFT) == 0)||((bits>>INDEX_SHIFT) == 1));
      return FA_DATA;
    }
    else if (bits & EF_DATA) {
      assert(!(bits & (CV_DATA|BF_DATA|FA_DATA|NO_DATA|LP_DATA) ));
      assert((bits>>INDEX_SHIFT) == 0);
      return EF_DATA;
    }
    else if (bits & NO_DATA) {
      assert(!(bits & (CV_DATA|BF_DATA|FA_DATA|EF_DATA|LP_DATA) ));
      return NO_DATA;
    }
    else if (bits & LP_DATA) {
      assert(!(bits & (CV_DATA|BF_DATA|FA_DATA|EF_DATA|NO_DATA) ));
      return LP_DATA;
    }
    else {
      return 0; // may be value data?
    }
  }

  int CtiData::getIndex() const {
    if (bits & CV_DATA) {
      assert(!(bits & (BF_DATA|FA_DATA|EF_DATA|NO_DATA|LP_DATA) ));
      assert((bits>>INDEX_SHIFT) == 0);
      return 0;
    }
    else if (bits & BF_DATA) {
      assert(!(bits & (CV_DATA|FA_DATA|EF_DATA|NO_DATA|LP_DATA) ));
      return bits>>INDEX_SHIFT;
    }
    else if (bits & FA_DATA) {
      assert(!(bits & (CV_DATA|BF_DATA|EF_DATA|NO_DATA|LP_DATA) ));
      assert(((bits>>INDEX_SHIFT) == 0)||((bits>>INDEX_SHIFT)==1));
      return bits>>INDEX_SHIFT;
    }
    else if (bits & EF_DATA) {
      assert(!(bits & (CV_DATA|BF_DATA|FA_DATA|NO_DATA|LP_DATA) ));
      assert((bits>>INDEX_SHIFT) == 0);
      return 0;
    }
    else if (bits & NO_DATA) {
      assert(!(bits & (CV_DATA|BF_DATA|FA_DATA|EF_DATA|LP_DATA) ));
      return bits>>INDEX_SHIFT;
    }
    else if (bits & LP_DATA) {
      assert(!(bits & (CV_DATA|BF_DATA|FA_DATA|EF_DATA|NO_DATA) ));
      return bits>>INDEX_SHIFT;
    }
    else {
      return 0; // may be value data?
    }
  }


  void CtiData::setDdeStuff(int8 * &_xora,DistributedDataExchanger * &_dde) {
    assert(xora_ptr == NULL);
    assert(dde_ptr == NULL);
    xora_ptr = &_xora;
    dde_ptr = &_dde;
  }
  void CtiData::initData() {
    assert(size_ptr != NULL);
    assert(mem_flag);
    if (getType() == IN_DATA) {
      this->alloc_in();
    }
    else if (getType() == DN_DATA) {
      this->alloc_dn();
    }
    else if (getType() == DN3_DATA) {
      this->alloc_dn3();
    }
    else if (getType() == DN33_DATA) {
      this->alloc_dn33();
    }
  }

  void CtiData::writeData(const string& name,MPI_File& fh,MPI_Offset& offset) {

    //if (mpi_rank == 0) cout << "CtiRegister::writeData \"" << name << "\"" << endl;

    assert((bits & WRITE_DATA)||(bits & CAN_WRITE_DATA));

    // check the dde if present....
    // TODO: get rid of this later...

    if (dde_ptr) {
      assert(xora_ptr);
      assert(size_ptr);
      assert((*size_ptr) >= 0);
      assert((*dde_ptr) != NULL);
      assert((*dde_ptr)->my_nda >= 0);
      assert((*size_ptr) == (*dde_ptr)->my_nda);
      // check 1-to-1 correspondence...
      int8 n_i8 = *size_ptr;
      int8 n_sum;
      MPI_Reduce(&n_i8,&n_sum,1,MPI_INT8,MPI_SUM,0,mpi_comm);
      int * send_buf = new int[*size_ptr];
      for (int i = 0; i < (*size_ptr); i++)
	send_buf[i] = 101010;
      const int nrecv = (*(*xora_ptr+mpi_rank+1))-(*(*xora_ptr+mpi_rank));
      int *recv_buf = new int[nrecv];
      for (int i = 0; i < nrecv; ++i)
	recv_buf[i] = 1;
      (*dde_ptr)->pull(send_buf,recv_buf);
      for (int i = 0; i < (*size_ptr); i++)
	assert(send_buf[i] == 1);
      for (int i = 0; i < nrecv; ++i)
	recv_buf[i] = 101;
      (*dde_ptr)->push(recv_buf,send_buf,REPLACE_DATA);
      for (int i = 0; i < nrecv; ++i)
	assert(recv_buf[i] == 1);
      delete[] recv_buf;
      delete[] send_buf;
    }

    // the data type can be determined from the bits...

    const int data_type = getType();
    if (data_type == D_DATA) {

      // a single real...

      if ( mpi_rank == 0 ) {
	cout << " > D0 " << name << "=" << d() << endl;
	Header header;
	header.id = UGP_IO_D0;
	sprintf(header.name,"%s",name.c_str());
	header.skip = header_size;
	header.rdata[0] = d();
	MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }
      offset += header_size;

    }
    else if (data_type == DN_DATA) {

      assert(xora_ptr);
      assert(dde_ptr);
      assert(data_ptr);
      assert(size_ptr);

      // 1. pack the strided data into a contiguous buffer and
      // push to stiped data using the xora/dde...

      double *send_buf = new double[*size_ptr];
      for (int ii = 0; ii < (*size_ptr); ++ii)
        send_buf[ii] = dn(ii);

      const int nrecv = (*(*xora_ptr+mpi_rank+1))-(*(*xora_ptr+mpi_rank));
      double *recv_buf = new double[nrecv];

      // hack -- fill recv buf to check dde push...
      for (int irecv = 0; irecv < nrecv; ++irecv)
        recv_buf[irecv] = 1.0E+20;

      (*dde_ptr)->push(recv_buf,send_buf);

      delete[] send_buf;

      // now write the recv_buf as the correct datatype...

      if ( mpi_rank == 0 ) {
	Header header;
	sprintf(header.name,"%s",name.c_str());
	const int topo = getUnindexedTopology();
	if (topo == CV_DATA) {
	  header.id = UGP_IO_CV_D1;
	  cout << " > CV_D1 " << name << endl;
	}
	else if (topo == BF_DATA) {
	  header.id = UGP_IO_BF_D1;
	  cout << " > BF_D1 " << name << endl;
	}
	else if (topo == FA_DATA) {
	  header.id = UGP_IO_FA_D1;
	  cout << " > FA_D1 " << name << endl;
	}
	else if (topo == LP_DATA) {
	  header.id = UGP_IO_LP_D1;
	  cout << " > LP_D1 " << name << endl;
	}
	else {
	  assert(0);
	}
	header.skip = header_size + (*(*xora_ptr+mpi_size))*double_size;
	ByteSwap::setLswMswPairForInt8(header.idata+0,(*(*xora_ptr+mpi_size)));
	MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      writeChunkedData<double>(fh,offset+header_size+(*(*xora_ptr+mpi_rank))*double_size,recv_buf,nrecv,mpi_comm);

      delete[] recv_buf;

      offset += header_size + (*(*xora_ptr+mpi_size))*double_size;

    }
    else if (data_type == DN3_DATA) {

      assert(xora_ptr);
      assert(dde_ptr);
      assert(data_ptr);
      assert(size_ptr);

      // 1. pack the stided data into a contiguous buffer and
      // push to stiped data using the xora/dde...

      double (*send_buf)[3] = new double[*size_ptr][3];
      for (int ii = 0; ii < (*size_ptr); ++ii)
	FOR_I3 send_buf[ii][i] = dn3(ii,i);

      const int nrecv = (*(*xora_ptr+mpi_rank+1))-(*(*xora_ptr+mpi_rank));
      double (*recv_buf)[3] = new double[nrecv][3];

      // hack -- fill recv buf to check dde push...
      for (int irecv = 0; irecv < nrecv; ++irecv)
	FOR_I3 recv_buf[irecv][i] = 1.0E+20;

      (*dde_ptr)->push(recv_buf,send_buf); // TODO I think we could avoid repacking into send_buf with a push_strided_d3 routine in dde

      delete[] send_buf;

      // now write the recv_buf as the correct datatype...

      if ( mpi_rank == 0 ) {
	Header header;
	sprintf(header.name,"%s",name.c_str());
	const int topo = getUnindexedTopology();
	if (topo == CV_DATA) {
	  header.id = UGP_IO_CV_D2;
	  cout << " > CV_D2 " << name << endl;
	}
	else if (topo == BF_DATA) {
	  header.id = UGP_IO_BF_D2;
	  cout << " > BF_D2 " << name << endl;
	}
	else if (topo == FA_DATA) {
	  header.id = UGP_IO_FA_D2;
	  cout << " > FA_D2 " << name << endl;
	}
	else if (topo == LP_DATA) {
	  header.id = UGP_IO_LP_D2;
	  cout << " > LP_D2 " << name << endl;
	}
	else {
	  assert(0);
	}
	header.skip   = header_size + (*(*xora_ptr+mpi_size))*double_size*3;
	ByteSwap::setLswMswPairForInt8(header.idata+0,(*(*xora_ptr+mpi_size)));
	header.idata[2] = 3;
	MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      writeChunkedData<double>(fh,offset+header_size+(*(*xora_ptr+mpi_rank))*double_size*3,(double*)recv_buf,nrecv*3,mpi_comm);

      delete[] recv_buf;

      offset += header_size + (*(*xora_ptr+mpi_size))*double_size*3;

    }
    else if (data_type == I_DATA) {

      // a single int...

      if ( mpi_rank == 0 ) {
	cout << " > I0 " << name << "=" << i() << endl;
	Header header;
	header.id = UGP_IO_I0;
	sprintf(header.name,"%s",name.c_str());
	header.skip = header_size;
	header.idata[0] = i();
	MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }
      offset += header_size;

    }
    else if (data_type == IN_DATA) {

      assert(xora_ptr);
      assert(dde_ptr);
      assert(data_ptr);
      assert(size_ptr);

      // 1. pack the strided data into a contiguous buffer and
      // push to stiped data using the xora/dde...

      int *send_buf = new int[*size_ptr];
      for (int ii = 0; ii < (*size_ptr); ++ii) {
        send_buf[ii] = in(ii);
      }

      const int nrecv = (*(*xora_ptr+mpi_rank+1))-(*(*xora_ptr+mpi_rank));
      int *recv_buf = new int[nrecv];

      // hack -- fill recv buf to check dde push...
      for (int irecv = 0; irecv < nrecv; ++irecv)
        recv_buf[irecv] = INT_MAX;

      (*dde_ptr)->push(recv_buf,send_buf);

      delete[] send_buf;

      // now write the recv_buf as the correct datatype...

      if ( mpi_rank == 0 ) {
	Header header;
	sprintf(header.name,"%s",name.c_str());
	const int topo = getUnindexedTopology();
	if (topo == CV_DATA) {
	  header.id     = UGP_IO_CV_I1;
	  cout << " > CV_I1 " << name << endl;
	}
	else if (topo == BF_DATA) {
	  header.id     = UGP_IO_BF_I1;
	  cout << " > BF_I1 " << name << endl;
	}
	else if (topo == FA_DATA) {
	  header.id     = UGP_IO_FA_I1;
	  cout << " > FA_I1 " << name << endl;
	}
	else if (topo == LP_DATA) {
	  header.id     = UGP_IO_LP_I1;
	  cout << " > LP_I1 " << name << endl;
	}
	else {
	  assert(0);
	}
	header.skip = header_size + (*(*xora_ptr+mpi_size))*int_size;
	ByteSwap::setLswMswPairForInt8(header.idata+0,(*(*xora_ptr+mpi_size)));
	MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      writeChunkedData<int>(fh,offset+header_size+(*(*xora_ptr+mpi_rank))*int_size,recv_buf,nrecv,mpi_comm);

      delete[] recv_buf;

      offset += header_size + (*(*xora_ptr+mpi_size))*int_size;

    }
    else {

      cout << "don't know how to write this yet" << endl;
      assert(0);

    }

  }

  void CtiData::redistReorderData(int* data) {
    assert(getType() == IN_DATA);
    assert(xora_ptr);
    assert(dde_ptr);
    assert(data_ptr);
    assert(size_ptr);

    int *send_buf = new int[*size_ptr];
    (*dde_ptr)->pull(send_buf,data);

    for (int ii = 0; ii < *size_ptr; ++ii)
      in(ii) = send_buf[ii];
    delete[] send_buf;

    setFlag(1);
  }
  void CtiData::redistReorderData(double* data) {
    assert(getType() == DN_DATA);
    assert(xora_ptr);
    assert(dde_ptr);
    assert(data_ptr);
    assert(size_ptr);

    double *send_buf = new double[*size_ptr];
    (*dde_ptr)->pull(send_buf,data);

    for (int ii = 0; ii < *size_ptr; ++ii)
      dn(ii) = send_buf[ii];
    delete[] send_buf;

    setFlag(1);
  }
  void CtiData::redistReorderData(double (*data)[3]) {
    assert(getType() == DN3_DATA);
    assert(xora_ptr);
    assert(dde_ptr);
    assert(data_ptr);
    assert(size_ptr);

    double (*send_buf)[3] = new double[*size_ptr][3];
    (*dde_ptr)->pull(send_buf,data);

    for (int ii = 0; ii < *size_ptr; ++ii)
      FOR_I3 dn3(ii,i) = send_buf[ii][i];
    delete[] send_buf;

    setFlag(1);
  }

  void CtiData::readData(MPI_File& fh,const Header& header,const MPI_Offset& offset,const bool byte_swap) {

    if (mpi_rank == 0) cout << "CtiRegister::readData \"" << header.name << "\"" << endl;

    assert(bits & READ_DATA);

    // check the dde if present....
    // TODO: get rid of this later...

    if (dde_ptr) {
      assert(xora_ptr);
      assert(size_ptr);
      assert((*size_ptr) == (*dde_ptr)->my_nda);
      // check 1-to-1 correspondence...
      int8 n_i8 = (*size_ptr);
      int8 n_sum;
      MPI_Reduce(&n_i8,&n_sum,1,MPI_INT8,MPI_SUM,0,mpi_comm);
      int * send_buf = new int[*size_ptr];
      for (int i = 0; i < (*size_ptr); i++)
	send_buf[i] = 101010;
      const int nrecv = (*(*xora_ptr+mpi_rank+1))-(*(*xora_ptr+mpi_rank));
      int *recv_buf = new int[nrecv];
      for (int i = 0; i < nrecv; ++i)
	recv_buf[i] = 1;
      (*dde_ptr)->pull(send_buf,recv_buf);
      for (int i = 0; i < (*size_ptr); i++)
	assert(send_buf[i] == 1);
      for (int i = 0; i < nrecv; ++i)
	recv_buf[i] = 101;
      (*dde_ptr)->push(recv_buf,send_buf,REPLACE_DATA);
      for (int i = 0; i < nrecv; ++i)
	assert(recv_buf[i] == 1); // see note above for change
      delete[] recv_buf;
      delete[] send_buf;
    }

    // the data type can be determined from the bits...

    const int data_type = getType();
    if (data_type == DN_DATA) {

      assert(xora_ptr);
      assert(dde_ptr);
      assert(data_ptr);
      assert(size_ptr);

      // the d1 size should be in Header...
      const int8 n_global = ByteSwap::getInt8FromLswMswPair(header.idata+0);
      assert(n_global == (*(*xora_ptr+mpi_size)));

      if (header.id == UGP_IO_CV_D1) {
	assert(getUnindexedTopology() == CV_DATA);
      }
      else if (header.id == UGP_IO_BF_D1) {
	assert(getUnindexedTopology() == BF_DATA);
      }
      else if (header.id == UGP_IO_FA_D1) {
	assert(getUnindexedTopology() == FA_DATA);
      }
      else if (header.id == UGP_IO_LP_D1) {
	assert(getUnindexedTopology() == LP_DATA);
      }
      else {
	assert(0);
      }

      // read from the file into strided data...

      const int nrecv = (*(*xora_ptr+mpi_rank+1))-(*(*xora_ptr+mpi_rank));
      double *recv_buf = new double[nrecv];
      readChunkedData<double>(fh,offset+header_size+(*(*xora_ptr+mpi_rank))*double_size,recv_buf,nrecv,byte_swap,mpi_comm);

      double *send_buf = new double[*size_ptr];
      (*dde_ptr)->pull(send_buf,recv_buf); // TODO could have similar pull routines for strided data...
      delete[] recv_buf;

      for (int ii = 0; ii < (*size_ptr); ++ii)
        dn(ii) = send_buf[ii]; // ptr[ii*stride+i];

      delete[] send_buf;

    }
    else if (data_type == DN3_DATA) {

      assert(xora_ptr);
      assert(dde_ptr);
      assert(data_ptr);
      assert(size_ptr);

      // the d2 size should be in Header...
      const int8 n_global = ByteSwap::getInt8FromLswMswPair(header.idata+0);
      assert(n_global == (*(*xora_ptr+mpi_size)));
      assert(header.idata[2] == 3);

      if (header.id == UGP_IO_CV_D2) {
	assert(getUnindexedTopology() == CV_DATA);
      }
      else if (header.id == UGP_IO_BF_D2) {
	assert(getUnindexedTopology() == BF_DATA);
      }
      else if (header.id == UGP_IO_FA_D2) {
	assert(getUnindexedTopology() == FA_DATA);
      }
      else if (header.id == UGP_IO_LP_D2) {
	assert(getUnindexedTopology() == LP_DATA);
      }
      else {
	assert(0);
      }

      // read from the file into strided data...

      const int nrecv = (*(*xora_ptr+mpi_rank+1))-(*(*xora_ptr+mpi_rank));
      double (*recv_buf)[3] = new double[nrecv][3];
      readChunkedData<double>(fh,offset+header_size+(*(*xora_ptr+mpi_rank))*double_size*3,(double*)recv_buf,nrecv*3,byte_swap,mpi_comm);

      double (*send_buf)[3] = new double[*size_ptr][3];
      (*dde_ptr)->pull(send_buf,recv_buf); // TODO could have similar pull routines for strided data...
      delete[] recv_buf;

      for (int ii = 0; ii < (*size_ptr); ++ii)
	FOR_I3 dn3(ii,i) = send_buf[ii][i]; // ptr[ii*stride+i];
      delete[] send_buf;

    }
    else if (data_type == IN_DATA) {

      assert(xora_ptr);
      assert(dde_ptr);
      assert(data_ptr);
      assert(size_ptr);

      // the d1 size should be in Header...
      const int8 n_global = ByteSwap::getInt8FromLswMswPair(header.idata+0);
      assert(n_global == (*(*xora_ptr+mpi_size)));

      if (header.id == UGP_IO_CV_I1) {
	assert(getUnindexedTopology() == CV_DATA);
      }
      else if (header.id == UGP_IO_BF_I1) {
	assert(getUnindexedTopology() == BF_DATA);
      }
      else if (header.id == UGP_IO_FA_I1) {
	assert(getUnindexedTopology() == FA_DATA);
      }
      else if (header.id == UGP_IO_LP_I1) {
	assert(getUnindexedTopology() == LP_DATA);
      }
      else {
	assert(0);
      }

      // read from the file into strided data...

      const int nrecv = (*(*xora_ptr+mpi_rank+1))-(*(*xora_ptr+mpi_rank));
      int *recv_buf = new int[nrecv];
      readChunkedData<int>(fh,offset+header_size+(*(*xora_ptr+mpi_rank))*int_size,recv_buf,nrecv,byte_swap,mpi_comm);

      int *send_buf = new int[*size_ptr];
      (*dde_ptr)->pull(send_buf,recv_buf); // TODO could have similar pull routines for strided data...
      delete[] recv_buf;

      for (int ii = 0; ii < (*size_ptr); ++ii)
        in(ii) = send_buf[ii];

      delete[] send_buf;

    }
    else {

      assert(0);

    }

  }

  map<const string,CtiData> registeredDataMap;
  map<const string,CtiData> currentDataMap;

  class NData {
  public:
    double (**x_ptr)[3];
    double (**l_ptr);
    int8 (**index_ptr);
    string name;
    int (**noote_ptr)[4];
    int *nte_ptr;
    NData() {
      x_ptr = NULL;
      l_ptr = NULL;
      index_ptr = NULL;
      name = "UNKNOWN";
      noote_ptr = NULL;
      nte_ptr = NULL;
    }
  };
  map<const int*,NData> nDataMap;

  bool CtiData::getX(double (*&x)[3]) const {
    map<const int*,NData>::iterator it = nDataMap.find(size_ptr);
    if (it == nDataMap.end()) {
      return false;
    }
    else if (it->second.x_ptr == NULL) {
      return false;
    }
    else {
      x = *it->second.x_ptr;
      return true;
    }
  }

  bool CtiData::getL(double *&delta) const {
    map<const int*,NData>::iterator it = nDataMap.find(size_ptr);
    if (it == nDataMap.end()) {
      return false;
    }
    else if (it->second.l_ptr == NULL) {
      return false;
    }
    else {
      delta = *it->second.l_ptr;
      return true;
    }
  }

  bool CtiData::getGlobalIndex(int8 *&index) const {
    map<const int*,NData>::iterator it = nDataMap.find(size_ptr);
    if (it == nDataMap.end()) {
      return false;
    }
    else if (it->second.index_ptr == NULL) {
      return false;
    }
    else {
      index = *it->second.index_ptr;
      return true;
    }
  }

  string CtiData::getName() const {
    map<const int*,NData>::iterator it = nDataMap.find(size_ptr);
    if (it == nDataMap.end()) {
      return "NOT_IN_NDATA";
    }
    else {
      return it->second.name;
    }
  }

  bool CtiData::hasTets() const {
    map<const int*,NData>::iterator it = nDataMap.find(size_ptr);
    if (it == nDataMap.end()) {
      return false;
    }
    else if (it->second.noote_ptr == NULL) {
      assert(it->second.nte_ptr == NULL);
      return false;
    }
    else {
      assert(it->second.nte_ptr != NULL);
      return true;
    }
  }

  bool CtiData::getTets(int (*&noote)[4],int &nte) const {
    map<const int*,NData>::iterator it = nDataMap.find(size_ptr);
    if (it == nDataMap.end()) {
      return false;
    }
    else if (it->second.noote_ptr == NULL) {
      assert(it->second.nte_ptr == NULL);
      return false;
    }
    else {
      noote = *it->second.noote_ptr;
      assert(it->second.nte_ptr != NULL);
      nte = *it->second.nte_ptr;
      return true;
    }
  }

  //map<const string,string> defineMap;
  list<CtiDataProducer*> ctiDataProducerList;

  void addCtiDataProducer(CtiDataProducer * cdp) {
    ctiDataProducerList.push_back(cdp);
  }
  void removeLastCtiDataProducer() {
    ctiDataProducerList.pop_back();
  }

  void dumpCurrentData() {
    if (mpi_rank == 0) {
      cout << "\n====================== currentDataMap ===========================" << endl;
      for (map<const string,CtiData>::iterator iter = currentDataMap.begin(); iter != currentDataMap.end(); ++iter)
	cout << "\"" << iter->first << "\" " << &(iter->second) << " " << iter->second << endl;
      cout << "==================== end of currentDataMap ==========================" << endl;
    }
  }

  void dumpRegisteredData() {
    if (mpi_rank == 0) {
      cout << "\n====================== registeredDataMap ===========================" << endl;
      for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter) {
	cout << " > \"" << iter->first << "\" type: " << iter->second.getTypeAsString()
          << " topo: " << iter->second.getTopologyAsString() << " flag: " << iter->second.getFlag() << endl;
      }
      cout << "==================== end of registeredDataMap ==========================" << endl;
    }
  }

  enum CtiDataTokenType {
    CTI_DATA_VARIABLE_NAME_TOKEN,
    CTI_DATA_FUNCTION_NAME_TOKEN,
    CTI_DATA_NUMBER_TOKEN,
    CTI_DATA_OPERATOR_TOKEN,
    CTI_DATA_EXPRESSION_TOKEN,
    CTI_DATA_COMMA_TOKEN,
    CTI_DATA_NOTHING_TOKEN,
    CTI_DATA_STRING_TOKEN,
    CTI_DATA_ERROR_TOKEN
  };

  const string cti_data_first_vars = ":_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  const string cti_data_all_vars = cti_data_first_vars+"1234567890.";
  const string cti_data_operators = "~=*/+-^<>!|&";

  const string cti_data_numbers = "1234567890";
  const string cti_data_numbers_plus_decimal = "1234567890.";
  const string cti_data_numbers_plus_e = cti_data_numbers_plus_decimal+"Ee";

  // ===============================================================================
  // CtiData evaluation routines...
  // ===============================================================================

  void _throwError(const string& message, const string& args, const string::size_type pos) {
    /*
    if (mpi_rank == 0) {
      cout << "\n\n=======================================================================\n" <<
	"Error: " << message << " \"" << args << "\"" << endl;
      for (unsigned int i = 0; i < message.length()+9+pos; ++i)
	cout << " ";
      cout << "^\n=======================================================================\n" << endl;
    }
    */
    if (cti_eval_verbosity) {
      CWARN("CtiData evaluation - " << message << " \"" << args << "\"");
    }
    throw(0);
  }

  void _advanceToPairedBracket(string::size_type& pos, const string& line) {
    int open_count = 1;
    string::size_type pos0 = pos;
    while (open_count != 0) {
      const string brackets = "()";
      pos = line.find_first_of(brackets,pos);
      if (pos == string::npos)
	_throwError("brackets not balanced",line,pos0-1);
      if (line[pos] == '(')
	++open_count;
      else
	--open_count;
      ++pos;
    }
  }

  CtiDataTokenType _parseNextToken(string& token, string::size_type& pos, const string& line) {

    pos = line.find_first_not_of(" \t",pos);
    if (pos == string::npos)
      return(CTI_DATA_NOTHING_TOKEN);

    //cout << "at \"" << line[pos] << "\"" << endl;

    if (cti_data_first_vars.find(line[pos]) != string::npos) {

      string::size_type pos0 = pos;
      pos = line.find_first_not_of(cti_data_all_vars,pos0+1);
      if (pos == string::npos) {
	token = line.substr(pos0);
	return(CTI_DATA_VARIABLE_NAME_TOKEN);
      }
      else {
	token = line.substr(pos0,pos-pos0);
	// this could be a function - check the next char. If it is a "(" then we have a function...
	if ((pos < line.length())&&(line[pos] == '(')) {
	  ++pos;
          // check that it actually isn't stats of a function, which is a var.
          // e.g. mach()_avg or log(p/rho)_avg.
          // we will assume for now that ')' will be followed by '_'
          string::size_type pos_tmp = pos;
          _advanceToPairedBracket(pos_tmp,line);
	  if ((pos_tmp < line.length())&&(line[pos_tmp] == '_')) {
            pos = line.find_first_not_of(cti_data_all_vars,pos_tmp+1);
	    token = line.substr(pos0,pos-pos0);
	    return(CTI_DATA_VARIABLE_NAME_TOKEN);
          }
          else {
	    return(CTI_DATA_FUNCTION_NAME_TOKEN);
          }
	}
	else {
	  return(CTI_DATA_VARIABLE_NAME_TOKEN);
	}
      }

    }
    else if (cti_data_operators.find(line[pos]) != string::npos) {

      string::size_type pos0 = pos;
      pos = line.find_first_not_of(cti_data_operators,pos0+1);
      if (pos == string::npos) {
	token = line.substr(pos0);
      }
      else {
	token = line.substr(pos0,pos-pos0);
      }
      return(CTI_DATA_OPERATOR_TOKEN);

    }
    else if (cti_data_numbers_plus_decimal.find(line[pos]) != string::npos) {

      string::size_type pos0 = pos;
      pos = line.find_first_not_of(cti_data_numbers_plus_e,pos0+1);
      if (pos == string::npos) {
	// make sure the last number is NOT an E/e. This means the
	// string did not properly complete the exponential...
	if ((line[line.length()-1] == 'E')||(line[line.length()-1] == 'e')) {
	  if (mpi_rank == 0)
	    cout << "Error: exponential not completed properly: \"" << line.substr(pos0) << "\"" << endl;
	  throw(0);
	}
	token = line.substr(pos0);
      }
      else {
	// look for a trailing E...
	if ((line[pos-1] == 'E')||(line[pos-1] == 'e')) {
	  if ((line[pos] == '+')||(line[pos] == '-')) {
	    // use 1.0E+20 notation...
	    if (line.length() == pos+1) {
	      if (mpi_rank == 0)
		cout << "Error parsing exponential: \"" << line.substr(pos0,pos-pos0) << "\"" << endl;
	      throw(0);
	    }
	    string::size_type pos1 = pos+1; // position of first part of number...
	    pos = line.find_first_not_of(cti_data_numbers,pos1);
	    if (pos == pos1) {
	      if (mpi_rank == 0)
		cout << "Error parsing exponential: \"" << line.substr(pos0,pos-pos0) << "\"" << endl;
	      throw(0);
	    }
	  }
	  else {
	    if (mpi_rank == 0)
	      cout << "Error parsing exponential: \"" << line.substr(pos0,pos-pos0) << "\"" << endl;
	    throw(0);
	  }
	}
	token = line.substr(pos0,pos-pos0);
      }
      return(CTI_DATA_NUMBER_TOKEN);

    }
    else if (line[pos] == '"') {

      string::size_type pos0 = pos;
      ++pos;
      while (pos < line.length()) {
	if (line[pos] == '"') {
	  token = line.substr(pos0+1,pos-pos0-1);
	  ++pos; // skip closing quote -- we have already handled it.
	  return(CTI_DATA_STRING_TOKEN);
	}
	++pos;
      }

      // if we got here, we did not find a closing quote, so throw an error...
      if (mpi_rank == 0)
	cout << "Error: could not find closing quote" << endl;
      throw(0);

    }
    else if (line[pos] == '(') {

      string::size_type pos0 = pos;
      ++pos;
      _advanceToPairedBracket(pos,line);

      // make sure there is something inside the brackets...
      assert(pos > pos0+1);

      // check that it actually isn't stats of a function, which is a var.
      // e.g. (r*u)_avg
      // we will assume for now that ')' will be followed by '_'
      if ((pos < line.length())&&(line[pos] == '_')) {
        pos = line.find_first_not_of(cti_data_all_vars,pos+1);
        token = line.substr(pos0,pos-pos0);
        return(CTI_DATA_VARIABLE_NAME_TOKEN);
      }
      else {
        token = line.substr(pos0+1,pos-pos0-2);
        return(CTI_DATA_EXPRESSION_TOKEN);
      }

    }
    else if (line[pos] == ',') {
      ++pos;
      return(CTI_DATA_COMMA_TOKEN);

    }
    else if (line[pos] == ')') {
      _throwError("unmatched closing bracket",line,pos);

    }
    else {
      _throwError("unrecognized char",line,pos);
    }

    // should NEVER get here...
    return(CTI_DATA_ERROR_TOKEN);

  }

  template <typename SimpleFunc2>
  CtiDataError CtiDataSimpleFunc2(CtiData& v, const CtiData * arg1, const CtiData * arg2, SimpleFunc2 simpleFunc2, const bool b_eval_func){

    const int type1 = arg1->getType();
    const int type2 = arg2->getType();

    if (type1 == I_DATA) {
      if (type2 == I_DATA) {
	v.new_i();
	if ( b_eval_func) v.i() = simpleFunc2(arg1->i(),arg2->i());
	return(CTI_DATA_OK);
      }
      else if (type2 == IN_DATA) {
	v.new_in(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.in(i) = simpleFunc2(arg1->i(),arg2->in(i));
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == D_DATA) {
	v.new_d();
	if ( b_eval_func) v.d() = simpleFunc2(double(arg1->i()),arg2->d());
	return(CTI_DATA_OK);
      }
      else if (type2 == D3_DATA) {
	v.new_d3();
        if ( b_eval_func) FOR_I3 v.d3(i) = simpleFunc2(double(arg1->i()),arg2->d3(i));
	return(CTI_DATA_OK);
      }
      else if (type2 == DN_DATA) {
	v.new_dn(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = simpleFunc2(double(arg1->i()),arg2->dn(i));
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == DN3_DATA) {
	v.new_dn3(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = simpleFunc2(double(arg1->i()),arg2->dn3(i,j));
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == DN33_DATA) {
	v.new_dn33(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = simpleFunc2(double(arg1->i()),arg2->dn33(i,j,k));
          }
	}
	return(CTI_DATA_OK);
      }
      else {
	return(CTI_DATA_NOT_VALID);
      }
    }
    else if (type1 == IN_DATA) {
      if (type2 == I_DATA) {
	v.new_in(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.in(i) = simpleFunc2(arg1->in(i),arg2->i());
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == IN_DATA) {
        v.new_in(*arg1,*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.in(i) = simpleFunc2(arg1->in(i),arg2->in(i));
          }
        }
	return(CTI_DATA_OK);
      }
      else if (type2 == D_DATA) {
	v.new_dn(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = simpleFunc2(double(arg1->in(i)),arg2->d());
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == D3_DATA) {
	v.new_dn3(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = simpleFunc2(double(arg1->in(i)),arg2->d3(j));
          }
        }
	return(CTI_DATA_OK);
      }
      else if (type2 == DN_DATA) {
        v.new_dn(*arg1,*arg2); // checks handled internally...
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = simpleFunc2(double(arg1->in(i)),arg2->dn(i));
          }
        }
	return(CTI_DATA_OK);
      }
      else if (type2 == DN3_DATA) {
	v.new_dn3(*arg1,*arg2); // checks handled internally...
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = simpleFunc2(double(arg1->in(i)),arg2->dn3(i,j));
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == DN33_DATA) {
	v.new_dn33(*arg1,*arg2); // checks handled internally...
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = simpleFunc2(double(arg1->in(i)),arg2->dn33(i,j,k));
          }
	}
	return(CTI_DATA_OK);
      }
      else {
	return(CTI_DATA_NOT_VALID);
      }
    }
    else if (type1 == D_DATA) {
      if (type2 == I_DATA) {
	v.new_d();
	if ( b_eval_func) v.d() = simpleFunc2(arg1->d(),double(arg2->i()));
	return(CTI_DATA_OK);
      }
      else if (type2 == IN_DATA) {
	v.new_dn(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = simpleFunc2(arg1->d(),double(arg2->in(i)));
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == D_DATA) {
	v.new_d();
	if ( b_eval_func) v.d() = simpleFunc2(arg1->d(),arg2->d());
	return(CTI_DATA_OK);
      }
      else if (type2 == D3_DATA) {
	v.new_d3();
        if ( b_eval_func) FOR_I3 v.d3(i) = simpleFunc2(arg1->d(),arg2->d3(i));
	return(CTI_DATA_OK);
      }
      else if (type2 == DN_DATA) {
	v.new_dn(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = simpleFunc2(arg1->d(),arg2->dn(i));
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == DN3_DATA) {
	v.new_dn3(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = simpleFunc2(arg1->d(),arg2->dn3(i,j));
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == DN33_DATA) {
	v.new_dn33(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = simpleFunc2(arg1->d(),arg2->dn33(i,j,k));
          }
	}
	return(CTI_DATA_OK);
      }
      else {
	return(CTI_DATA_NOT_VALID);
      }
    }
    else if (type1 == D3_DATA) {
      if (type2 == I_DATA) {
	v.new_d3();
	if ( b_eval_func) FOR_I3 v.d3(i) = simpleFunc2(arg1->d3(i),double(arg2->i()));
	return(CTI_DATA_OK);
      }
      else if (type2 == IN_DATA) {
	v.new_dn3(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = simpleFunc2(arg1->d3(j),double(arg2->in(i)));
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == D_DATA) {
	v.new_d3();
	if ( b_eval_func) FOR_I3 v.d3(i) = simpleFunc2(arg1->d3(i),arg2->d());
	return(CTI_DATA_OK);
      }
      else if (type2 == D3_DATA) {
	v.new_d3();
        if ( b_eval_func) FOR_I3 v.d3(i) = simpleFunc2(arg1->d3(i),arg2->d3(i));
	return(CTI_DATA_OK);
      }
      else if (type2 == DN_DATA) {
	v.new_dn3(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = simpleFunc2(arg1->d3(j),arg2->dn(i));
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == DN3_DATA) {
	v.new_dn3(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = simpleFunc2(arg1->d3(j),arg2->dn3(i,j));
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == DN33_DATA) {
	v.new_dn33(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = simpleFunc2(arg1->d3(j),arg2->dn33(i,j,k));
          }
	}
	return(CTI_DATA_OK);
      }
      else {
	return(CTI_DATA_NOT_VALID);
      }
    }
    else if (type1 == DN_DATA) {
      if (type2 == I_DATA) {
	v.new_dn(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = simpleFunc2(arg1->dn(i),double(arg2->i()));
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == IN_DATA) {
        v.new_dn(*arg1,*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = simpleFunc2(arg1->dn(i),double(arg2->in(i)));
          }
        }
	return(CTI_DATA_OK);
      }
      else if (type2 == D_DATA) {
	v.new_dn(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = simpleFunc2(arg1->dn(i),arg2->d());
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == D3_DATA) {
	v.new_dn3(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = simpleFunc2(arg1->dn(i),arg2->d3(j));
          }
        }
	return(CTI_DATA_OK);
      }
      else if (type2 == DN_DATA) {
        v.new_dn(*arg1,*arg2); // checks handled internally...
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = simpleFunc2(arg1->dn(i),arg2->dn(i));
          }
        }
	return(CTI_DATA_OK);
      }
      else if (type2 == DN3_DATA) {
	v.new_dn3(*arg1,*arg2); // checks handled internally...
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = simpleFunc2(arg1->dn(i),arg2->dn3(i,j));
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == DN33_DATA) {
	v.new_dn33(*arg1,*arg2); // checks handled internally...
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = simpleFunc2(arg1->dn(i),arg2->dn33(i,j,k));
          }
	}
	return(CTI_DATA_OK);
      }
      else {
	return(CTI_DATA_NOT_VALID);
      }
    }
    else if (type1 == DN3_DATA) {
      if (type2 == I_DATA) {
	v.new_dn3(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = simpleFunc2(arg1->dn3(i,j),double(arg2->i()));
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == IN_DATA) {
	v.new_dn3(*arg1,*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = simpleFunc2(arg1->dn3(i,j),double(arg2->in(i)));
          }
        }
	return(CTI_DATA_OK);
      }
      else if (type2 == D_DATA) {
	v.new_dn3(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = simpleFunc2(arg1->dn3(i,j),arg2->d());
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == D3_DATA) {
	v.new_dn3(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = simpleFunc2(arg1->dn3(i,j),arg2->d3(j));
          }
        }
	return(CTI_DATA_OK);
      }
      else if (type2 == DN_DATA) {
	v.new_dn3(*arg1,*arg2); // checks handled internally...
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = simpleFunc2(arg1->dn3(i,j),arg2->dn(i));
          }
        }
	return(CTI_DATA_OK);
      }
      else if (type2 == DN3_DATA) {
	v.new_dn3(*arg1,*arg2); // checks handled internally...
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = simpleFunc2(arg1->dn3(i,j),arg2->dn3(i,j));
          }
        }
	return(CTI_DATA_OK);
      }
      else if (type2 == DN33_DATA) {
	v.new_dn33(*arg1,*arg2); // checks handled internally...
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = simpleFunc2(arg1->dn3(i,j),arg2->dn33(i,j,k));
          }
        }
	return(CTI_DATA_OK);
      }
      else {
	return(CTI_DATA_NOT_VALID);
      }
    }
    else if (type1 == DN33_DATA) {
      if (type2 == I_DATA) {
	v.new_dn33(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = simpleFunc2(arg1->dn33(i,j,k),double(arg2->i()));
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == IN_DATA) {
	v.new_dn33(*arg1,*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = simpleFunc2(arg1->dn33(i,j,k),double(arg2->in(i)));
          }
        }
	return(CTI_DATA_OK);
      }
      else if (type2 == D_DATA) {
	v.new_dn33(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = simpleFunc2(arg1->dn33(i,j,k),arg2->d());
          }
	}
	return(CTI_DATA_OK);
      }
      else if (type2 == D3_DATA) {
	v.new_dn33(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = simpleFunc2(arg1->dn33(i,j,k),arg2->d3(j));
          }
        }
	return(CTI_DATA_OK);
      }
      else if (type2 == DN_DATA) {
	v.new_dn33(*arg1,*arg2); // checks handled internally...
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = simpleFunc2(arg1->dn33(i,j,k),arg2->dn(i));
          }
        }
	return(CTI_DATA_OK);
      }
      else if (type2 == DN3_DATA) {
	v.new_dn33(*arg1,*arg2); // checks handled internally...
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = simpleFunc2(arg1->dn33(i,j,k),arg2->dn3(i,j));
          }
        }
	return(CTI_DATA_OK);
      }
      else if (type2 == DN33_DATA) {
	v.new_dn33(*arg1,*arg2); // checks handled internally...
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = simpleFunc2(arg1->dn33(i,j,k),arg2->dn33(i,j,k));
          }
        }
	return(CTI_DATA_OK);
      }
      else {
	return(CTI_DATA_NOT_VALID);
      }
    }
    else {
      assert(0);
      return(CTI_DATA_NOT_VALID);
    }

  }

  double opPlusFunctionDouble(const double &a, const double &b){return a+b;}
  double opMinusFunctionDouble(const double &a, const double &b){return a-b;}
  double opMultiplyFunctionDouble(const double &a, const double &b){return a*b;}
  double opDivideFunctionDouble(const double &a, const double &b){return a/b;}
  double opPowerFunctionDouble(const double &a, const double &b){return pow(a,b);}
  double opLessFunctionDouble(const double &a, const double &b){return (double)(a<b);}
  double opLessEqualFunctionDouble(const double &a, const double &b){return (double)(a<=b);}
  double opGreaterFunctionDouble(const double &a, const double &b){return (double)(a>b);}
  double opGreaterEqualFunctionDouble(const double &a, const double &b){return (double)(a>=b);}
  double opEqualFunctionDouble(const double &a, const double &b){return (double)(a==b);}
  double opNotEqualFunctionDouble(const double &a, const double &b){return (double)(a!=b);}
  double opOrFunctionDouble(const double &a, const double &b){return (double)((int)a||(int)b);}
  double opAndFunctionDouble(const double &a, const double &b){return (double)((int)a&&(int)b);}
  int opPlusFunctionInt(const int &a, const int &b){return a+b;}
  int opMinusFunctionInt(const int &a, const int &b){return a-b;}
  int opMultiplyFunctionInt(const int &a, const int &b){return a*b;}
  int opDivideFunctionInt(const int &a, const int &b){return a/b;}
  int opPowerFunctionInt(const int &a, const int &b){
    // pow(int,int) not always available...
    //return pow(a,b);
    if (b < 0) return 0; // cannot return an int, except for a = 1, -1
    else if (b == 0) return 1;
    else {
      int b_copy = b;
      int result = 1;
      while (b_copy > 0) {
        result *= a;
        --b_copy;
      }
      return result;
    }
  }
  int opLessFunctionInt(const int &a, const int &b){return (int)(a<b);}
  int opLessEqualFunctionInt(const int &a, const int &b){return (int)(a<=b);}
  int opGreaterFunctionInt(const int &a, const int &b){return (int)(a>b);}
  int opGreaterEqualFunctionInt(const int &a, const int &b){return (int)(a>=b);}
  int opEqualFunctionInt(const int &a, const int &b){return (int)(a==b);}
  int opNotEqualFunctionInt(const int &a, const int &b){return (int)(a!=b);}
  int opOrFunctionInt(const int &a, const int &b){return (int)(a||b);}
  int opAndFunctionInt(const int &a, const int &b){return (int)(a&&b);}

  CtiDataError _evalOperator(CtiData& v, const CtiData * arg1, const CtiData * arg2, const string& op,
                             const bool b_eval_func) {

    bool b_double = false;
    if ( (arg1->getType() == D_DATA)||(arg1->getType() == D3_DATA)||(arg1->getType() == DN_DATA)||(arg1->getType() == DN3_DATA)||(arg1->getType() == DN33_DATA)||
         (arg2->getType() == D_DATA)||(arg2->getType() == D3_DATA)||(arg2->getType() == DN_DATA)||(arg2->getType() == DN3_DATA)||(arg2->getType() == DN33_DATA) ) {
      b_double = true;
    }

    if (op == "+"){
      if (b_double)
        CtiDataSimpleFunc2(v,arg1,arg2,&opPlusFunctionDouble,b_eval_func);
      else
        CtiDataSimpleFunc2(v,arg1,arg2,&opPlusFunctionInt,b_eval_func);
    }
    else if (op == "-") {
      if (b_double)
        CtiDataSimpleFunc2(v,arg1,arg2,&opMinusFunctionDouble,b_eval_func);
      else
        CtiDataSimpleFunc2(v,arg1,arg2,&opMinusFunctionInt,b_eval_func);
    }
    else if (op == "*") {
      if (b_double)
        CtiDataSimpleFunc2(v,arg1,arg2,&opMultiplyFunctionDouble,b_eval_func);
      else
        CtiDataSimpleFunc2(v,arg1,arg2,&opMultiplyFunctionInt,b_eval_func);
    }
    else if (op == "/") {
      if (b_double)
        CtiDataSimpleFunc2(v,arg1,arg2,&opDivideFunctionDouble,b_eval_func);
      else
        CtiDataSimpleFunc2(v,arg1,arg2,&opDivideFunctionInt,b_eval_func);
    }
    else if ((op == "**")||(op == "^")) {
      if (b_double)
        CtiDataSimpleFunc2(v,arg1,arg2,&opPowerFunctionDouble,b_eval_func);
      else
        CtiDataSimpleFunc2(v,arg1,arg2,&opPowerFunctionInt,b_eval_func);
    }
    else if (op == "<") {
      if (b_double)
        CtiDataSimpleFunc2(v,arg1,arg2,&opLessFunctionDouble,b_eval_func);
      else
        CtiDataSimpleFunc2(v,arg1,arg2,&opLessFunctionInt,b_eval_func);
    }
    else if (op == "<~") {
      if (b_double)
        CtiDataSimpleFunc2(v,arg1,arg2,&opLessEqualFunctionDouble,b_eval_func);
      else
        CtiDataSimpleFunc2(v,arg1,arg2,&opLessEqualFunctionInt,b_eval_func);
    }
    else if (op == ">") {
      if (b_double)
        CtiDataSimpleFunc2(v,arg1,arg2,&opGreaterFunctionDouble,b_eval_func);
      else
        CtiDataSimpleFunc2(v,arg1,arg2,&opGreaterFunctionInt,b_eval_func);
    }
    else if (op == ">~") {
      if (b_double)
        CtiDataSimpleFunc2(v,arg1,arg2,&opGreaterEqualFunctionDouble,b_eval_func);
      else
        CtiDataSimpleFunc2(v,arg1,arg2,&opGreaterEqualFunctionInt,b_eval_func);
    }
    else if (op == "~~") {
      if (b_double)
        CtiDataSimpleFunc2(v,arg1,arg2,&opEqualFunctionDouble,b_eval_func);
      else
        CtiDataSimpleFunc2(v,arg1,arg2,&opEqualFunctionInt,b_eval_func);
    }
    else if (op == "!~") {
      if (b_double)
        CtiDataSimpleFunc2(v,arg1,arg2,&opNotEqualFunctionDouble,b_eval_func);
      else
        CtiDataSimpleFunc2(v,arg1,arg2,&opNotEqualFunctionInt,b_eval_func);
    }
    else if (op == "||") {
      if (b_double)
        CtiDataSimpleFunc2(v,arg1,arg2,&opOrFunctionDouble,b_eval_func);
      else
        CtiDataSimpleFunc2(v,arg1,arg2,&opOrFunctionInt,b_eval_func);
    }
    else if (op == "&&") {
      if (b_double)
        CtiDataSimpleFunc2(v,arg1,arg2,&opAndFunctionDouble,b_eval_func);
      else
        CtiDataSimpleFunc2(v,arg1,arg2,&opAndFunctionInt,b_eval_func);
    }
    else {
      CERR("unrecognized operator: " << (*arg1) << " " << op << " " << (*arg2));
    }

    return(CTI_DATA_OK);

  }

  void _processStack(CtiData& v,vector<pair<bool,CtiData*> >& valueStack,vector<string>& opStack,
                     const bool b_eval_func) {

    //cout << "CtiRegister::_processStack..." << endl;

    /*
    cout << "\n\nvalueStack: " << endl;
    for (int i = 0; i < valueStack.size(); ++i)
      cout << " > bool: " << valueStack[i].first << " data: " << *(valueStack[i].second) << endl;
    cout << "opStack: " << endl;
    for (int i = 0; i < opStack.size(); ++i)
      cout << " > " << opStack[i] << endl;
    cout << "\n" << endl;
    */

    if (valueStack.size() == 1) {
      pair<bool,CtiData*> value = valueStack.back(); valueStack.pop_back();
      // copy value into v...
      v.copy(*(value.second));
      if (value.first) {
	// true means the valueStack allocated the memory for this CtiData*, so
	// since we just copied into v, we reset before delete to leave any pointers
	// valid...
	value.second->reset();
	delete value.second;
      }
      else {
	// in this case, the valueStack contained a CtiData ptr that it did not
	// allocate. In this case, we must set the mem_flag in the copy to false,
	// so it does not reset the memory when deleted...
	v.setMemFlag(false);
      }
    }
    else {
      assert(valueStack.size() > 1);
      pair<bool,CtiData*> value1 = valueStack.back(); valueStack.pop_back();
      while (1) {
	string op = opStack.back(); opStack.pop_back();
	pair<bool,CtiData*> value0 = valueStack.back(); valueStack.pop_back();
	if (_evalOperator(v,value0.second,value1.second,op,b_eval_func) != CTI_DATA_OK) {
	  // the valueStack will be cleaned when this error is caught,
	  // but we have already popped 2 values off, and need to clean them here...
	  if (value0.first) delete value0.second;
	  if (value1.first) delete value1.second;
	  throw(0);
	}
	if (value0.first) delete value0.second;
	if (opStack.size() == 0) {
	  if (value1.first) delete value1.second;
	  break;
	}
	// if we made it here, then we need to put v into value1, and reset v...
	if (value1.first) delete value1.second;
	value1.first = true;
	value1.second = new CtiData();
	value1.second->copy(v);
	v.reset();
      }
    }
  }

  CtiDataError _maxCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() == 1) {

      // 2. grab the arg...
      list<CtiData>::iterator arg = args.begin();

      // 3. and apply the function...
      if (arg->getType() == IN_DATA) {
	v.new_i();
	if ( b_eval_func) {
	  int my_max_v = -INT_MAX;
	  for (int i = 0; i < arg->size(); ++i) {
	    my_max_v = max(my_max_v,arg->in(i));
	  }
	  int max_v;
	  MPI_Allreduce(&my_max_v,&max_v,1,MPI_INT,MPI_MAX,mpi_comm);
	  v.i() = max_v;
	}
	return(CTI_DATA_OK);
      }
      else if (arg->getType() == DN_DATA) {
	v.new_d();
	if ( b_eval_func) {
	  double my_max_v = -HUGE_VAL;
	  for (int i = 0; i < arg->size(); ++i) {
	    my_max_v = max(my_max_v,arg->dn(i));
	  }
	  double max_v;
	  MPI_Allreduce(&my_max_v,&max_v,1,MPI_DOUBLE,MPI_MAX,mpi_comm);
	  v.d() = max_v;
	}
	return(CTI_DATA_OK);
      }
      return(CTI_DATA_NOT_VALID);

    }
    else if (args.size() == 2) {

      // the max of 2 args take the larger...

      list<CtiData>::iterator arg1 = args.begin();
      list<CtiData>::iterator arg2 = arg1; arg2++;

      if (arg1->getType() == I_DATA) {
        if (arg2->getType() == I_DATA) {
          v.new_i();
          if ( b_eval_func) {
            v.i() = max(arg1->i(),arg2->i());
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == IN_DATA) {
          v.new_in(*arg2);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.in(i) = max(arg1->i(),arg2->in(i));
            }
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == D_DATA) {
          v.new_d();
          if ( b_eval_func) {
            v.d() = max(double(arg1->i()),arg2->d());
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == DN_DATA) {
          v.new_dn(*arg2);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = max(double(arg1->i()),arg2->dn(i));
            }
          }
          return(CTI_DATA_OK);
        }
        return(CTI_DATA_NOT_VALID);
      }
      else if (arg1->getType() == IN_DATA) {
        if (arg2->getType() == I_DATA) {
          v.new_in(*arg1);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.in(i) = max(arg1->in(i),arg2->i());
            }
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == IN_DATA) {
          v.new_in(*arg1,*arg2);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.in(i) = max(arg1->in(i),arg2->in(i));
            }
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == D_DATA) {
          v.new_dn(*arg1);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = max(double(arg1->in(i)),arg2->d());
            }
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == DN_DATA) {
          v.new_dn(*arg1,*arg2);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = max(double(arg1->in(i)),arg2->dn(i));
            }
          }
          return(CTI_DATA_OK);
        }
        return(CTI_DATA_NOT_VALID);
      }
      else if (arg1->getType() == D_DATA) {
        if (arg2->getType() == I_DATA) {
          v.new_d();
          if ( b_eval_func) {
            v.d() = max(arg1->d(),double(arg2->i()));
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == IN_DATA) {
          v.new_dn(*arg2);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = max(arg1->d(),double(arg2->in(i)));
            }
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == D_DATA) {
          v.new_d();
          if ( b_eval_func) {
            v.d() = max(arg1->d(),arg2->d());
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == DN_DATA) {
          v.new_dn(*arg2);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = max(arg1->d(),arg2->dn(i));
            }
          }
          return(CTI_DATA_OK);
        }
        return(CTI_DATA_NOT_VALID);
      }
      else if (arg1->getType() == DN_DATA) {
        if (arg2->getType() == I_DATA) {
          v.new_dn(*arg1);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = max(arg1->dn(i),double(arg2->i()));
            }
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == IN_DATA) {
          v.new_dn(*arg1,*arg2);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = max(arg1->dn(i),double(arg2->in(i)));
            }
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == D_DATA) {
          v.new_dn(*arg1);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = max(arg1->dn(i),arg2->d());
            }
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == DN_DATA) {
          v.new_dn(*arg1,*arg2);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = max(arg1->dn(i),arg2->dn(i));
            }
          }
          return(CTI_DATA_OK);
        }
        return(CTI_DATA_NOT_VALID);
      }
      return(CTI_DATA_NOT_VALID);

    }
    else {

      return(CTI_DATA_ARG_COUNT);

    }

  }

  CtiDataError _minCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() == 1) {

      // 2. grab the arg...
      list<CtiData>::iterator arg = args.begin();

      // 3. and apply the function...
      if (arg->getType() == IN_DATA) {
	v.new_i();
	if ( b_eval_func) {
	  int my_min_v = INT_MAX;
	  for (int i = 0; i < arg->size(); ++i) {
	    my_min_v = min(my_min_v,arg->in(i));
	  }
	  int min_v;
	  MPI_Allreduce(&my_min_v,&min_v,1,MPI_INT,MPI_MIN,mpi_comm);
	  v.i() = min_v;
	}
	return(CTI_DATA_OK);
      }
      else if (arg->getType() == DN_DATA) {
	v.new_d();
	if ( b_eval_func) {
	  double my_min_v = HUGE_VAL;
	  for (int i = 0; i < arg->size(); ++i) {
	    my_min_v = min(my_min_v,arg->dn(i));
	  }
	  double min_v;
	  MPI_Allreduce(&my_min_v,&min_v,1,MPI_DOUBLE,MPI_MIN,mpi_comm);
	  v.d() = min_v;
	}
	return(CTI_DATA_OK);
      }
      return(CTI_DATA_NOT_VALID);

    }
    else if (args.size() == 2) {

      // the min of 2 args take the smaller...

      list<CtiData>::iterator arg1 = args.begin();
      list<CtiData>::iterator arg2 = arg1; arg2++;

      if (arg1->getType() == I_DATA) {
        if (arg2->getType() == I_DATA) {
          v.new_i();
          if ( b_eval_func) {
            v.i() = min(arg1->i(),arg2->i());
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == IN_DATA) {
          v.new_in(*arg2);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.in(i) = min(arg1->i(),arg2->in(i));
            }
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == D_DATA) {
          v.new_d();
          if ( b_eval_func) {
            v.d() = min(double(arg1->i()),arg2->d());
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == DN_DATA) {
          v.new_dn(*arg2);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = min(double(arg1->i()),arg2->dn(i));
            }
          }
          return(CTI_DATA_OK);
        }
        return(CTI_DATA_NOT_VALID);
      }
      else if (arg1->getType() == IN_DATA) {
        if (arg2->getType() == I_DATA) {
          v.new_in(*arg1);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.in(i) = min(arg1->in(i),arg2->i());
            }
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == IN_DATA) {
          v.new_in(*arg1,*arg2);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.in(i) = min(arg1->in(i),arg2->in(i));
            }
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == D_DATA) {
          v.new_dn(*arg1);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = min(double(arg1->in(i)),arg2->d());
            }
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == DN_DATA) {
          v.new_dn(*arg1,*arg2);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = min(double(arg1->in(i)),arg2->dn(i));
            }
          }
          return(CTI_DATA_OK);
        }
        return(CTI_DATA_NOT_VALID);
      }
      else if (arg1->getType() == D_DATA) {
        if (arg2->getType() == I_DATA) {
          v.new_d();
          if ( b_eval_func) {
            v.d() = min(arg1->d(),double(arg2->i()));
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == IN_DATA) {
          v.new_dn(*arg2);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = min(arg1->d(),double(arg2->in(i)));
            }
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == D_DATA) {
          v.new_d();
          if ( b_eval_func) {
            v.d() = min(arg1->d(),arg2->d());
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == DN_DATA) {
          v.new_dn(*arg2);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = min(arg1->d(),arg2->dn(i));
            }
          }
          return(CTI_DATA_OK);
        }
        return(CTI_DATA_NOT_VALID);
      }
      else if (arg1->getType() == DN_DATA) {
        if (arg2->getType() == I_DATA) {
          v.new_dn(*arg1);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = min(arg1->dn(i),double(arg2->i()));
            }
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == IN_DATA) {
          v.new_dn(*arg1,*arg2);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = min(arg1->dn(i),double(arg2->in(i)));
            }
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == D_DATA) {
          v.new_dn(*arg1);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = min(arg1->dn(i),arg2->d());
            }
          }
          return(CTI_DATA_OK);
        }
        else if (arg2->getType() == DN_DATA) {
          v.new_dn(*arg1,*arg2);
          if ( b_eval_func) {
            for (int i = 0; i < v.size(); ++i) {
              v.dn(i) = min(arg1->dn(i),arg2->dn(i));
            }
          }
          return(CTI_DATA_OK);
        }
        return(CTI_DATA_NOT_VALID);
      }
      return(CTI_DATA_NOT_VALID);

    }
    else {

      return(CTI_DATA_ARG_COUNT);

    }

  }

  CtiDataError _sumCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();

    // 3. and apply the function...
    if (arg->getType() == IN_DATA) {
      v.new_i();
      if ( b_eval_func) {
        int my_sum_v = 0;
        for (int i = 0; i < arg->size(); ++i) {
          my_sum_v += arg->in(i);
        }
        int sum_v;
        MPI_Allreduce(&my_sum_v,&sum_v,1,MPI_INT,MPI_SUM,mpi_comm);
        v.i() = sum_v;
      }
      return(CTI_DATA_OK);
    }
    else if (arg->getType() == DN_DATA) {
      v.new_d();
      if ( b_eval_func) {
        double my_sum_v = 0.0;
        for (int i = 0; i < arg->size(); ++i) {
          my_sum_v += arg->dn(i);
        }
        double sum_v;
        MPI_Allreduce(&my_sum_v,&sum_v,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
        v.d() = sum_v;
      }
      return(CTI_DATA_OK);
    }
    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _avgCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();

    // 3. and apply the function...
    if (arg->getType() == IN_DATA) {
      v.new_d();
      if ( b_eval_func) {
        double my_buf_v[2] = {0.0,0.0}; // sum, count
        for (int i = 0; i < arg->size(); ++i) {
          my_buf_v[0] += double(arg->in(i));
        }
        my_buf_v[1] = arg->size();
        double buf_v[2];
        MPI_Allreduce(my_buf_v,buf_v,2,MPI_DOUBLE,MPI_SUM,mpi_comm);
        v.d() = buf_v[0]/buf_v[1];
      }
      return(CTI_DATA_OK);
    }
    else if (arg->getType() == DN_DATA) {
      v.new_d();
      if ( b_eval_func) {
        double my_buf_v[2] = {0.0,0.0}; // sum, count
        for (int i = 0; i < arg->size(); ++i) {
          my_buf_v[0] += arg->dn(i);
        }
        my_buf_v[1] = arg->size();
        double buf_v[2];
        MPI_Allreduce(my_buf_v,buf_v,2,MPI_DOUBLE,MPI_SUM,mpi_comm);
        v.d() = buf_v[0]/buf_v[1];
      }
      return(CTI_DATA_OK);
    }
    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _distCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 2)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the args...
    list<CtiData>::iterator arg1 = args.begin();
    list<CtiData>::iterator arg2 = arg1; arg2++;

    // 3. and apply the function...
    if (arg1->getType() == D3_DATA) {
      if (arg2->getType() == D3_DATA) {
        v.new_d();
        if ( b_eval_func) {
          v.d() = sqrt( (arg1->d3(0)-arg2->d3(0))*(arg1->d3(0)-arg2->d3(0))
                      + (arg1->d3(1)-arg2->d3(1))*(arg1->d3(1)-arg2->d3(1))
                      + (arg1->d3(2)-arg2->d3(2))*(arg1->d3(2)-arg2->d3(2)) );
        }
        return(CTI_DATA_OK);
      }
      else if (arg2->getType() == DN3_DATA) {
        v.new_dn(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = sqrt( (arg1->d3(0)-arg2->dn3(i,0))*(arg1->d3(0)-arg2->dn3(i,0))
                          + (arg1->d3(1)-arg2->dn3(i,1))*(arg1->d3(1)-arg2->dn3(i,1))
                          + (arg1->d3(2)-arg2->dn3(i,2))*(arg1->d3(2)-arg2->dn3(i,2)) );
          }
        }
        return(CTI_DATA_OK);
      }
    }
    else if (arg1->getType() == DN3_DATA) {
      if (arg2->getType() == D3_DATA) {
        v.new_dn(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = sqrt( (arg2->d3(0)-arg1->dn3(i,0))*(arg2->d3(0)-arg1->dn3(i,0))
                          + (arg2->d3(1)-arg1->dn3(i,1))*(arg2->d3(1)-arg1->dn3(i,1))
                          + (arg2->d3(2)-arg1->dn3(i,2))*(arg2->d3(2)-arg1->dn3(i,2)) );
          }
        }
        return(CTI_DATA_OK);
      }
      else if (arg2->getType() == DN3_DATA) {
        v.new_dn(*arg2);
        if ( b_eval_func ) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = sqrt( (arg1->dn3(i,0)-arg2->dn3(i,0))*(arg1->dn3(i,0)-arg2->dn3(i,0))
                          + (arg1->dn3(i,1)-arg2->dn3(i,1))*(arg1->dn3(i,1)-arg2->dn3(i,1))
                          + (arg1->dn3(i,2)-arg2->dn3(i,2))*(arg1->dn3(i,2)-arg2->dn3(i,2)) );
          }
        }
        return(CTI_DATA_OK);
      }
    }
    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _d3CtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 3)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg1 = args.begin();
    list<CtiData>::iterator arg2 = arg1; arg2++;
    list<CtiData>::iterator arg3 = arg2; arg3++;

    // 3. and apply the function...
    if (arg1->getType() == D_DATA && arg2->getType() == D_DATA && arg3->getType() == D_DATA) {
      v.new_d3();
      if ( b_eval_func) {
        v.d3(0) = arg1->d();
        v.d3(1) = arg2->d();
        v.d3(2) = arg3->d();
      }
      return(CTI_DATA_OK);
    }
    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _dotCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 2)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg1 = args.begin();
    list<CtiData>::iterator arg2 = arg1; arg2++;

    // 3. and apply the function...
    if (arg1->getType() == D3_DATA) {
      if (arg2->getType() == D3_DATA) {
        v.new_d();
        if ( b_eval_func) v.d() = arg1->d3(0)*arg2->d3(0) + arg1->d3(1)*arg2->d3(1) + arg1->d3(2)*arg2->d3(2);
        return(CTI_DATA_OK);
      }
      else if (arg2->getType() == DN3_DATA) {
        v.new_dn(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = arg1->d3(0)*arg2->dn3(i,0) + arg1->d3(1)*arg2->dn3(i,1) + arg1->d3(2)*arg2->dn3(i,2);
          }
        }
        return(CTI_DATA_OK);
      }
      else if (arg2->getType() == DN33_DATA) {
        v.new_dn3(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = arg1->d3(0)*arg2->dn33(i,0,j) + arg1->d3(1)*arg2->dn33(i,1,j) + arg1->d3(2)*arg2->dn33(i,2,j);
          }
        }
        return(CTI_DATA_OK);
      }
    }
    else if (arg1->getType() == DN3_DATA) {
      if (arg2->getType() == D3_DATA) {
        v.new_dn(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = arg2->d3(0)*arg1->dn3(i,0) + arg2->d3(1)*arg1->dn3(i,1) + arg2->d3(2)*arg1->dn3(i,2);
          }
        }
        return(CTI_DATA_OK);
      }
      else if (arg2->getType() == DN3_DATA) {
        v.new_dn(*arg1,*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = arg1->dn3(i,0)*arg2->dn3(i,0) + arg1->dn3(i,1)*arg2->dn3(i,1) + arg1->dn3(i,2)*arg2->dn3(i,2);
          }
        }
        return(CTI_DATA_OK);
      }
      else if (arg2->getType() == DN33_DATA) {
        v.new_dn3(*arg1,*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = arg1->dn3(i,0)*arg2->dn33(i,0,j) + arg1->dn3(i,1)*arg2->dn33(i,1,j) + arg1->dn3(i,2)*arg2->dn33(i,2,j);
          }
        }
        return(CTI_DATA_OK);
      }
    }
    else if (arg1->getType() == DN33_DATA) {
      if (arg2->getType() == D3_DATA) {
        v.new_dn3(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = arg1->dn33(i,j,0)*arg2->d3(0) + arg1->dn33(i,j,1)*arg2->d3(1) + arg1->dn33(i,j,2)*arg2->d3(2);
          }
        }
        return(CTI_DATA_OK);
      }
      else if (arg2->getType() == DN3_DATA) {
        v.new_dn3(*arg1,*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = arg1->dn33(i,j,0)*arg2->dn3(i,0) + arg1->dn33(i,j,1)*arg2->dn3(i,1) + arg1->dn33(i,j,2)*arg2->dn3(i,2);
          }
        }
        return(CTI_DATA_OK);
      }
    }
    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _ddotCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 2)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg1 = args.begin();
    list<CtiData>::iterator arg2 = arg1; arg2++;

    // 3. and apply the function...
    if (arg1->getType() == DN33_DATA) {
      if (arg2->getType() == DN33_DATA) {
        v.new_dn(*arg1,*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = arg1->dn33(i,0,0)*arg2->dn33(i,0,0) + arg1->dn33(i,0,1)*arg2->dn33(i,1,0) + arg1->dn33(i,0,2)*arg2->dn33(i,2,0) +
                      arg1->dn33(i,1,0)*arg2->dn33(i,0,1) + arg1->dn33(i,1,1)*arg2->dn33(i,1,1) + arg1->dn33(i,1,2)*arg2->dn33(i,2,1) +
                      arg1->dn33(i,2,0)*arg2->dn33(i,0,2) + arg1->dn33(i,2,1)*arg2->dn33(i,1,2) + arg1->dn33(i,2,2)*arg2->dn33(i,2,2);
          }
        }
        return(CTI_DATA_OK);
      }
    }
    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _crossCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 2)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg1 = args.begin();
    list<CtiData>::iterator arg2 = arg1; arg2++;

    // 3. and apply the function...
    if (arg1->getType() == D3_DATA) {
      if (arg2->getType() == D3_DATA) {
        v.new_d3();
        if ( b_eval_func) {
          v.d3(0) = arg1->d3(1)*arg2->d3(2) - arg1->d3(2)*arg2->d3(1);
          v.d3(1) = arg1->d3(2)*arg2->d3(0) - arg1->d3(0)*arg2->d3(2);
          v.d3(2) = arg1->d3(0)*arg2->d3(1) - arg1->d3(1)*arg2->d3(0);
        }
        return(CTI_DATA_OK);
      }
      else if (arg2->getType() == DN3_DATA) {
        v.new_dn3(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn3(i,0) = arg1->d3(1)*arg2->dn3(i,2) - arg1->d3(2)*arg2->dn3(i,1);
            v.dn3(i,1) = arg1->d3(2)*arg2->dn3(i,0) - arg1->d3(0)*arg2->dn3(i,2);
            v.dn3(i,2) = arg1->d3(0)*arg2->dn3(i,1) - arg1->d3(1)*arg2->dn3(i,0);
          }
        }
        return(CTI_DATA_OK);
      }
    }
    else if (arg1->getType() == DN3_DATA) {
      if (arg2->getType() == D3_DATA) {
        v.new_dn3(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn3(i,0) = arg1->dn3(i,1)*arg2->d3(2) - arg1->dn3(i,2)*arg2->d3(1);
            v.dn3(i,1) = arg1->dn3(i,2)*arg2->d3(0) - arg1->dn3(i,0)*arg2->d3(2);
            v.dn3(i,2) = arg1->dn3(i,0)*arg2->d3(1) - arg1->dn3(i,1)*arg2->d3(0);
          }
        }
        return(CTI_DATA_OK);
      }
      else if (arg2->getType() == DN3_DATA) {
        v.new_dn3(*arg1,*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn3(i,0) = arg1->dn3(i,1)*arg2->dn3(i,2) - arg1->dn3(i,2)*arg2->dn3(i,1);
            v.dn3(i,1) = arg1->dn3(i,2)*arg2->dn3(i,0) - arg1->dn3(i,0)*arg2->dn3(i,2);
            v.dn3(i,2) = arg1->dn3(i,0)*arg2->dn3(i,1) - arg1->dn3(i,1)*arg2->dn3(i,0);
          }
        }
        return(CTI_DATA_OK);
      }
    }
    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _magCtiData(CtiData& v,list<CtiData>& args, const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();

    // 3. and apply the function...
    if (arg->getType() == D3_DATA) {
      v.new_d();
      if ( b_eval_func) v.d() = sqrt(arg->d3(0)*arg->d3(0) + arg->d3(1)*arg->d3(1) + arg->d3(2)*arg->d3(2));
      return(CTI_DATA_OK);
    }
    else if (arg->getType() == DN3_DATA) {
      v.new_dn(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          v.dn(i) = sqrt(arg->dn3(i,0)*arg->dn3(i,0) + arg->dn3(i,1)*arg->dn3(i,1) + arg->dn3(i,2)*arg->dn3(i,2));
        }
      }
      return(CTI_DATA_OK);
    }
    else if (arg->getType() == DN33_DATA) {
      v.new_dn(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          v.dn(i) = sqrt(0.5*(
                arg->dn33(i,0,0)*arg->dn33(i,0,0) + arg->dn33(i,0,1)*arg->dn33(i,0,1) + arg->dn33(i,0,2)*arg->dn33(i,0,2) +
                arg->dn33(i,1,0)*arg->dn33(i,1,0) + arg->dn33(i,1,1)*arg->dn33(i,1,1) + arg->dn33(i,1,2)*arg->dn33(i,1,2) +
                arg->dn33(i,2,0)*arg->dn33(i,2,0) + arg->dn33(i,2,1)*arg->dn33(i,2,1) + arg->dn33(i,2,2)*arg->dn33(i,2,2)));
        }
      }
      return(CTI_DATA_OK);
    }
    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _dyadicCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 2)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg1 = args.begin();
    list<CtiData>::iterator arg2 = arg1; arg2++;

    // 3. and apply the function...
    if (arg1->getType() == D3_DATA) {
      if (arg2->getType() == DN3_DATA) {
        v.new_dn33(*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = arg1->d3(j)*arg2->dn3(i,k);
          }
        }
        return(CTI_DATA_OK);
      }
    }
    else if (arg1->getType() == DN3_DATA) {
      if (arg2->getType() == D3_DATA) {
        v.new_dn33(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = arg1->dn3(i,j)*arg2->d3(k);
          }
        }
        return(CTI_DATA_OK);
      }
      else if (arg2->getType() == DN3_DATA) {
        v.new_dn33(*arg1,*arg2);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = arg1->dn3(i,j)*arg2->dn3(i,k);
          }
        }
        return(CTI_DATA_OK);
      }
    }
    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _transposeCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();

    // 3. and apply the function...
    if (arg->getType() == DN33_DATA) {
      v.new_dn33(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          FOR_J3 FOR_K3 v.dn33(i,k,j) = arg->dn33(i,j,k);
        }
      }
    }
    return(CTI_DATA_NOT_VALID);

  }

  template<double (*T)(double)>
  CtiDataError _funcOfCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();

    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == D_DATA) {
      v.new_d();
      if ( b_eval_func) v.d() = T(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      if ( b_eval_func) FOR_I3 v.d3(i) = T(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          v.dn(i) = T(arg->dn(i));
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          FOR_J3 v.dn3(i,j) = T(arg->dn3(i,j));
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN33_DATA) {
      v.new_dn33(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          FOR_J3 FOR_K3 v.dn33(i,j,k) = T(arg->dn33(i,j,k));
        }
      }
      return(CTI_DATA_OK);
    }

    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _absCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();

    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == I_DATA) {
      v.new_i();
      if ( b_eval_func) v.i() = abs(arg->i());
      return(CTI_DATA_OK);
    }
    else if (datatype == IN_DATA) {
      v.new_in(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          v.in(i) = abs(arg->in(i));
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == D_DATA) {
      v.new_d();
      if ( b_eval_func) v.d() = fabs(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      if ( b_eval_func) FOR_I3 v.d3(i) = fabs(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          v.dn(i) = fabs(arg->dn(i));
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          FOR_J3 v.dn3(i,j) = fabs(arg->dn3(i,j));
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN33_DATA) {
      v.new_dn33(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          FOR_J3 FOR_K3 v.dn33(i,j,k) = fabs(arg->dn33(i,j,k));
        }
      }
      return(CTI_DATA_OK);
    }

    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _sgnCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();

    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == I_DATA) {
      v.new_i();
      if ( b_eval_func) {
        if (arg->i() > 0)
          v.i() = 1;
        else if (arg->i() < 0)
          v.i() = -1;
        else
          v.i() = 0;
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == IN_DATA) {
      v.new_in(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          if (arg->in(i) > 0)
            v.in(i) = 1;
          else if (arg->in(i) < 0)
            v.in(i) = -1;
          else
            v.in(i) = 0;
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == D_DATA) {
      v.new_d();
      if ( b_eval_func) v.d() = SGN(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      if ( b_eval_func) FOR_I3 v.d3(i) = SGN(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          v.dn(i) = SGN(arg->dn(i));
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          FOR_J3 v.dn3(i,j) = SGN(arg->dn3(i,j));
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN33_DATA) {
      v.new_dn33(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          FOR_J3 FOR_K3 v.dn33(i,j,k) = SGN(arg->dn33(i,j,k));
        }
      }
      return(CTI_DATA_OK);
    }

    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _intCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();

    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == I_DATA) {
      v.new_i();
      if ( b_eval_func) v.i() = arg->i();
      return(CTI_DATA_OK);
    }
    else if (datatype == IN_DATA) {
      v.new_in(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          v.in(i) = arg->in(i);
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == D_DATA) {
      v.new_i();
      if ( b_eval_func) v.i() = int(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_in(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          v.in(i) = int(arg->dn(i));
        }
      }
      return(CTI_DATA_OK);
    }
    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _doubleCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();

    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == I_DATA) {
      v.new_d();
      if ( b_eval_func) v.d() = double(arg->i());
      return(CTI_DATA_OK);
    }
    else if (datatype == IN_DATA) {
      v.new_dn(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          v.dn(i) = double(arg->in(i));
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == D_DATA) {
      v.new_d();
      if ( b_eval_func) v.d() = arg->d();
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          v.dn(i) = arg->dn(i);
        }
      }
      return(CTI_DATA_OK);
    }
    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _compCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if ((args.size() != 2)&&(args.size() != 3))
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg1 = args.begin();
    list<CtiData>::iterator arg2 = arg1; ++arg2;

    // 3. and apply the function...
    if (arg2->getType() == D_DATA) {
      if (args.size() == 2) {

        // convert to ints...
        int index = -1;
        if (arg2->d() == 0.0) index = 0;
        else if (arg2->d() == 1.0) index = 1;
        else if (arg2->d() == 2.0) index = 2;

        if ( (index >= 0) && (index < 3) ) {
          if (arg1->getType() == D3_DATA) {
            v.new_d();
            if ( b_eval_func) v.d() = arg1->d3(index);
            return(CTI_DATA_OK);
          }
          else if (arg1->getType() == DN3_DATA) {
            v.new_dn(*arg1);
            if ( b_eval_func) {
              for (int i = 0; i < v.size(); ++i) {
                v.dn(i) = arg1->dn3(i,index);
              }
            }
            return(CTI_DATA_OK);
          }
        }
      }
      else if (args.size() == 3) {
        list<CtiData>::iterator arg3 = arg2; ++arg3;
        if (arg3->getType() == D_DATA) {

          // convert to ints...
          int index[2] = {-1,-1};
          if (arg2->d() == 0.0) index[0] = 0;
          else if (arg2->d() == 1.0) index[0] = 1;
          else if (arg2->d() == 2.0) index[0] = 2;
          if (arg3->d() == 0.0) index[1] = 0;
          else if (arg3->d() == 1.0) index[1] = 1;
          else if (arg3->d() == 2.0) index[1] = 2;

          if ( (index[0] >= 0) && (index[0] < 3) && (index[1] >= 0) && (index[1] < 3) ) {
            if (arg1->getType() == DN33_DATA) {
              v.new_dn(*arg1);
              if ( b_eval_func) {
                for (int i = 0; i < v.size(); ++i) {
                  v.dn(i) = arg1->dn33(i,index[0],index[1]);
                }
              }
              return(CTI_DATA_OK);
            }
          }
        }
      }
    }
    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _powCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 2)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg1 = args.begin();
    list<CtiData>::iterator arg2 = arg1; ++arg2;

    // 3. and apply the function...
    if ( arg2->getType() == D_DATA ) {
      if (arg1->getType() == D_DATA) {
        v.new_d();
        if ( b_eval_func) v.d() = pow(arg1->d(),arg2->d());
        return(CTI_DATA_OK);
      }
      else if (arg1->getType() == D3_DATA) {
        v.new_d3();
        if ( b_eval_func) FOR_I3 v.d3(i) = pow(arg1->d3(i),arg2->d());
        return(CTI_DATA_OK);
      }
      else if (arg1->getType() == DN_DATA) {
        v.new_dn(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            v.dn(i) = pow(arg1->dn(i),arg2->d());
          }
        }
        return(CTI_DATA_OK);
      }
      else if (arg1->getType() == DN3_DATA) {
        v.new_dn3(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 v.dn3(i,j) = pow(arg1->dn3(i,j),arg2->d());
          }
        }
        return(CTI_DATA_OK);
      }
      else if (arg1->getType() == DN33_DATA) {
        v.new_dn33(*arg1);
        if ( b_eval_func) {
          for (int i = 0; i < v.size(); ++i) {
            FOR_J3 FOR_K3 v.dn33(i,j,k) = pow(arg1->dn33(i,j,k),arg2->d());
          }
        }
        return(CTI_DATA_OK);
      }
    }
    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _atan2CtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {

    // 1. check arg count...
    if (args.size() != 2)
      return(CTI_DATA_ARG_COUNT);

    // 2. grab the arg...
    list<CtiData>::iterator arg1 = args.begin();
    list<CtiData>::iterator arg2 = arg1; ++arg2;

    // 3. and apply the function...
    if ((arg1->getType()==D_DATA )&&(arg2->getType()==D_DATA)) {
      v.new_d();
      if ( b_eval_func) v.d() = atan2(arg1->d(),arg2->d());
      return(CTI_DATA_OK);
    }
    else if ((arg1->getType()==DN_DATA)&&(arg2->getType()==DN_DATA)) {
      v.new_dn(*arg1,*arg2);
      if ( b_eval_func) {
        for (int i = 0; i < v.size(); ++i) {
          v.dn(i) = atan2(arg1->dn(i),arg2->dn(i));
        }
      }
      return(CTI_DATA_OK);
    }
    return(CTI_DATA_NOT_VALID);

  }

  CtiDataError _funcEval(CtiData& v, const string& name, list<CtiData>& args, const bool b_eval_func) {

    //  CERR("_funcEval not implemented for name: " << name);
    //enum CtiDataError {
    //  CTI_DATA_OK,
    //  CTI_DATA_ARG_COUNT,
    //  CTI_DATA_NOT_FOUND,
    //  CTI_DATA_NOT_VALID
    //};

    if (name == "mag") {

      // ----------------------------------------------
      // v = mag(u)
      // returns the magnitude of u in each position
      // ----------------------------------------------

      return _magCtiData(v,args,b_eval_func);

    }
    else if (name == "transpose") {

      // ----------------------------------------------
      // v = transpose(u)
      // returns the transpose of tensor u in each position
      // ----------------------------------------------

      return _transposeCtiData(v,args,b_eval_func);

    }
    else if (name == "dyadic") {

      // ----------------------------------------------
      // v = dyadic(u1,u2)
      // returns the dyadic product of vectors u1 and u2 in each position
      // ----------------------------------------------

      return _dyadicCtiData(v,args,b_eval_func);
    }
    else if (name == "dot") {

      // ----------------------------------------------
      // v = dot(u1,u2)
      // returns the dot product of u1 and u2 in each position
      // ----------------------------------------------

      return _dotCtiData(v,args,b_eval_func);

    }
    else if (name == "ddot") {

      // ----------------------------------------------
      // v = ddot(u1,u2)
      // returns the double dot product of tensors u1 and u2 in each position
      // ----------------------------------------------

      return _ddotCtiData(v,args,b_eval_func);

    }
    else if (name == "cross") {

      // ----------------------------------------------
      // v = cross(u1,u2)
      // returns the cross product of u1 and u2 in each position
      // ----------------------------------------------

      return _crossCtiData(v,args,b_eval_func);

    }
    else if (name == "d3") {

      // ----------------------------------------------
      // v = d3(a,b,c)
      // returns the D3_DATA representation of doubles a b c
      // ----------------------------------------------

      return _d3CtiData(v,args,b_eval_func);

    }
    else if (name == "dist") {

      // ----------------------------------------------
      // v = dist(u1,u2)
      // returns the dist b/w u1 and u2 in each position
      // ----------------------------------------------

      return _distCtiData(v,args,b_eval_func);

    }
    else if (name == "sqrt") {

      // ----------------------------------------------
      // v = sqrt(u)
      // returns the sqrt of u in each position
      // ----------------------------------------------

      return _funcOfCtiData<sqrt>(v,args,b_eval_func);
    }
    else if (name == "exp") {

      // ----------------------
      // ----------------------------------------------
      // v = exp(u)
      // returns the exp of u in each position
      // ----------------------
      // ----------------------------------------------

      return _funcOfCtiData<exp>(v,args,b_eval_func);

    }
    else if (name == "log") {

      // ----------------------------------------------
      // v = log(u)
      // returns the (natural) log of u in each position
      // ----------------------------------------------

      return _funcOfCtiData<log>(v,args,b_eval_func);

    }
    else if (name == "log10") {

      // ----------------------------------------------
      // v = log10(u)
      // returns the (common) log10 of u in each position
      // ----------------------------------------------

      return _funcOfCtiData<log10>(v,args,b_eval_func);

    }
    else if (name == "cos") {

      // ----------------------------------------------
      // v = cos(u)
      // returns the cos of u in each position
      // ----------------------------------------------

      return _funcOfCtiData<cos>(v,args,b_eval_func);

    }
    else if (name == "sin") {

      // ----------------------------------------------
      // v = sin(u)
      // returns the sin of u in each position
      // ----------------------------------------------

      return _funcOfCtiData<sin>(v,args,b_eval_func);

    }
    else if (name == "tan") {

      // ----------------------------------------------
      // v = tan(u)
      // returns the tan of u in each position
      // ----------------------------------------------

      return _funcOfCtiData<tan>(v,args,b_eval_func);

    }
    else if (name == "cosh") {

      // ----------------------------------------------
      // v = cosh(u)
      // returns the cosh of u in each position
      // ----------------------------------------------

      return _funcOfCtiData<cosh>(v,args,b_eval_func);

    }
    else if (name == "sinh") {

      // ----------------------------------------------
      // v = sinh(u)
      // returns the sinh of u in each position
      // ----------------------------------------------

      return _funcOfCtiData<sinh>(v,args,b_eval_func);

    }
    else if (name == "tanh") {

      // ----------------------------------------------
      // v = tanh(u)
      // returns the tanh of u in each position
      // ----------------------------------------------

      return _funcOfCtiData<tanh>(v,args,b_eval_func);

    }
    else if (name == "acos") {

      // ----------------------------------------------
      // v = acos(u)
      // returns the acos of u in each position
      // ----------------------------------------------

      return _funcOfCtiData<acos>(v,args,b_eval_func);

    }
    else if (name == "asin") {

      // ----------------------------------------------
      // v = asin(u)
      // returns the asin of u in each position
      // ----------------------------------------------

      return _funcOfCtiData<asin>(v,args,b_eval_func);

    }
    else if (name == "atan") {

      // ----------------------------------------------
      // v = atan(u)
      // returns the atan of u in each position
      // ----------------------------------------------

      return _funcOfCtiData<atan>(v,args,b_eval_func);

    }
    else if (name == "ceil") {

      // ----------------------------------------------
      // v = ceil(u)
      // returns the ceil of u in each position
      // ----------------------------------------------

      return _funcOfCtiData<ceil>(v,args,b_eval_func);

    }
    else if (name == "floor") {

      // ----------------------------------------------
      // v = floor(u)
      // returns the floor of u in each position
      // ----------------------------------------------

      return _funcOfCtiData<floor>(v,args,b_eval_func);

    }
    else if (name == "abs") {

      // ----------------------------------------------
      // v = abs(u)
      // returns the abs of u in each position
      // ----------------------------------------------

      return _absCtiData(v,args,b_eval_func);

    }
    else if (name == "sgn") {

      // ----------------------------------------------
      // v = sgn(u)
      // returns the sign of u in each position
      // ----------------------------------------------

      return _sgnCtiData(v,args,b_eval_func);

    }
    else if (name == "comp") {

      // ----------------------------------------------
      // v = comp(u,i) or comp(u,i,j)
      // returns component i(,j) of u
      // ----------------------------------------------

      return _compCtiData(v,args,b_eval_func);

    }
    else if (name == "pow") {

      // ----------------------------------------------
      // v = pow(u,i)
      // returns u raised to the power i at each position
      // ----------------------------------------------

      return _powCtiData(v,args,b_eval_func);

    }
    else if (name == "atan2") {
      // ----------------------------------------------
      // v = atan2(y,x)
      // returns the 4-quadrant arctangent of (y,x)
      // ----------------------------------------------

      return _atan2CtiData(v,args,b_eval_func);

    }
    else if (name == "max") {
      // ----------------------------------------------
      // v = max(u)
      // returns max element of u
      // ----------------------------------------------

      return _maxCtiData(v,args,b_eval_func);

    }
    else if (name == "min") {
      // ----------------------------------------------
      // v = min(u)
      // returns min element of u
      // ----------------------------------------------

      return _minCtiData(v,args,b_eval_func);
    }
    else if (name == "sum") {
      // ----------------------------------------------
      // v = sum(u)
      // returns sum of u's elements
      // ----------------------------------------------

      return _sumCtiData(v,args,b_eval_func);
    }
    else if (name == "avg") {
      // ----------------------------------------------
      // v = avg(u)
      // returns avg of u's elements
      // ----------------------------------------------

      return _avgCtiData(v,args,b_eval_func);
    }
    else if (name == "double") {
      // ----------------------------------------------
      // v = double(u)
      // return double cast u's elements
      // ----------------------------------------------

      return _doubleCtiData(v,args,b_eval_func);

    }
    else if (name == "int") {
      // ----------------------------------------------
      // v = int(u)
      // return int cast u's elements
      // ----------------------------------------------

      return _intCtiData(v,args,b_eval_func);

    }
    else {

      // we did not find the function in our standard set of functions. Try to
      // find it in and CtiDataProducers that are registered with us....

      for (list<CtiDataProducer*>::iterator iter = ctiDataProducerList.begin(); iter != ctiDataProducerList.end(); ++iter) {
	const CtiDataError ierr = (*iter)->funcEvalCtiData(v,name,args,b_eval_func);
	if ((ierr == CTI_DATA_OK)||(ierr == CTI_DATA_ARG_COUNT)||(ierr == CTI_DATA_NOT_VALID))
	  return ierr;
	assert(ierr == CTI_DATA_NOT_FOUND);
      }

      //CERR("_funcEval not implemented for name: " << name);
      return(CTI_DATA_NOT_FOUND);

    }

  }

  int _getOpRank(const string& op) {
    if ((op == "+=")||(op == "-=")||(op == "*=")||(op == "/=")||(op == "=")||(op == "~"))
      return(0);
    else if ((op == "||")||(op == "&&"))
      return(1);
    else if ((op == "<")||(op == "<~")||(op == ">")||(op == ">~")||(op == "~~")||(op == "!~"))
      return(2);
    else if ((op == "+")||(op == "-"))
      return(3);
    else if ((op == "*")||(op == "/"))
      return(4);
    else if ((op == "**")||(op == "^"))
      return(5);
    else {
      if (mpi_rank == 0)
	cout << "Error: unrecognized op: " << op << endl;
      throw(0);
    }
  }

  string::size_type _eval(CtiData& v,const bool b_eval_func, const string& line,string::size_type pos = 0) {

    vector<pair<bool,CtiData*> > valueStack; // bool is true CtiData* needs deletion
    vector<string> opStack;

    try {

      // start to parse...
      while (1) {

	string token;
	const string::size_type pos0 = pos;
	CtiDataTokenType tt = _parseNextToken(token,pos,line);

	switch (tt) {
	case CTI_DATA_NOTHING_TOKEN:

	  // ===========================================================
	  // we got nothing more, so process anything left on the stack
	  // into "v" and return...
	  // ===========================================================
	  if (valueStack.size() < opStack.size()+1)
	    _throwError("expecting value",line,pos0);
	  else if (valueStack.size() > opStack.size()+1)
	    _throwError("expecting operator",line,pos0);
	  assert(pos == string::npos);
	  _processStack(v,valueStack,opStack,b_eval_func);
	  assert(valueStack.size() == 0);
	  assert(opStack.size() == 0);
	  return pos;

	case CTI_DATA_COMMA_TOKEN:

	  // ===========================================================
	  // we got a comma, so process anything left on the stack
	  // into "v" and return...
	  // ===========================================================
	  // process anything left in the stack and return...
	  assert(valueStack.size() == opStack.size()+1);
	  assert(pos < string::npos);
	  _processStack(v,valueStack,opStack,b_eval_func);
	  assert(valueStack.size() == 0);
	  assert(opStack.size() == 0);
	  return pos;

	case CTI_DATA_VARIABLE_NAME_TOKEN:

	  // ===========================================================
	  // we got a variable name (alpha-numeric). Try and figure
	  // out what it refers to...
	  // ===========================================================

	  {

	    if (valueStack.size() > opStack.size())
	      _throwError("expecting operator",line,pos0);

	    map<const string,CtiData>::iterator iter = registeredDataMap.find(token);
	    if (iter != registeredDataMap.end()) {
	      //cout << "setting from token: " << token << " data type: " << iter->second.getTypeAsString() << endl;
	      valueStack.push_back(pair<bool,CtiData*>(false,&(iter->second)));
	    }
            else {
              // could actually live in current data map now that we added varEvalCtiData functionality
              iter = currentDataMap.find(token);
              if (iter != currentDataMap.end()) {
                //cout << "setting from token: " << token << " data type: " << iter->second.getTypeAsString() << endl;
                valueStack.push_back(pair<bool,CtiData*>(false,&(iter->second)));
              }
              else {
                // varEvalCtiData is defined by some CtiDataProducers (e.g. LesImageMapper)
                list<CtiDataProducer*>::iterator iter;
                for (iter = ctiDataProducerList.begin(); iter != ctiDataProducerList.end(); ++iter) {
                  if (CtiData * cti_data = (*iter)->varEvalCtiData(token)) {
                    valueStack.push_back(pair<bool,CtiData*>(false,cti_data));
                    break;
                  }
                }
                if (iter == ctiDataProducerList.end()) {
                  _throwError("unrecognized token",line,pos0);
                  /*
                    map<const string,string>::iterator iter = defineMap.find(token);
                    if (iter != defineMap.end()) {
                    valueStack.push_back(pair<bool,CtiData*>(true,new CtiData()));
                    _eval(*(valueStack.back().second),iter->second);
                    }
                    else {
                    _throwError("unrecognized token",line,pos0);
                    }
                  */
                }
              }
	    }
	  }
	  break;

	case CTI_DATA_EXPRESSION_TOKEN:

	  // ===========================================================
	  // we got a separate expression, so we need to evaluate...
	  // ===========================================================

	  if (valueStack.size() > opStack.size())
	    _throwError("expecting operator",line,pos0);
	  // an expression requires some temporary memory...
	  valueStack.push_back(pair<bool,CtiData*>(true,new CtiData()));
	  _eval(*(valueStack.back().second),b_eval_func,token);
	  break;

	case CTI_DATA_FUNCTION_NAME_TOKEN:

	  // ===========================================================
	  // we know this is a function because it ended in an open bracket. Advance to
	  // the paired close bracket, then look for the function in the map of
	  // functions. If it exists, then put that on the valueStack, otherwise
	  // and eval the function argument...
	  // ===========================================================
	  {

	    if (valueStack.size() > opStack.size())
	      _throwError("expecting operator",line,pos0);

	    string::size_type pos1 = pos;
	    _advanceToPairedBracket(pos,line);

	    map<const string,CtiData>::iterator iter = currentDataMap.find(line.substr(pos0,pos-pos0));
	    if (iter != currentDataMap.end()) {

	      //cout << "_eval CTI_DATA_FUNCTION_NAME_TOKEN: found \"" << iter->first << "\" in currentDataMap." << endl;
	      valueStack.push_back(pair<bool,CtiData*>(false,&(iter->second)));

	    }
	    else {

	      pair<const string,CtiData> key_value_pair(line.substr(pos0,pos-pos0),CtiData());

	      // evaluate the argument or arguments into a list (could be multiple args,
	      // also supports no args)...

	      string subline = line.substr(pos1,pos-pos1-1);
	      list<CtiData> args;
	      string::size_type pos = subline.find_first_not_of(" \t",0);
	      while (pos != string::npos) {
		args.push_back(CtiData());
		pos = _eval(args.back(),b_eval_func,subline,pos);
	      }

	      // insert into the current data map before eval -- this should be ok

	      //cout << "XXXX: inserting \"" << key_value_pair.first << "\" as a function in the currentDataMap" << endl;

	      pair<map<const string,CtiData>::iterator,bool> return_pair = currentDataMap.insert(key_value_pair);
	      assert(return_pair.second); // should insert successfully
	      assert(key_value_pair.second.empty()); // should be empty - don't use anymore, and his destructor cannot affect the inserted CtiData

	      // check token against solver specific functions (i.e chemtable) as well as
	      // CtiDataMachine functions like mag(), comp()
	      CtiDataError ierr = _funcEval(return_pair.first->second,token,args,b_eval_func);

	      // recall possible returned values...
	      // enum CtiDataError {
	      // CTI_DATA_OK,
	      // CTI_DATA_ARG_COUNT,
	      // CTI_DATA_NOT_FOUND,
	      // CTI_DATA_NOT_VALID
              // }
	      switch (ierr) {
	      case CTI_DATA_OK:
		// fine -- add the valueStack...
		valueStack.push_back(pair<bool,CtiData*>(false,&(return_pair.first->second)));
		break;
	      case CTI_DATA_ARG_COUNT:
		currentDataMap.erase(return_pair.first);
		_throwError("wrong arg count in funcEval",line,pos0);
	      case CTI_DATA_NOT_FOUND:
		currentDataMap.erase(return_pair.first);
		_throwError("could not find function in funcEval",line,pos0);
	      case CTI_DATA_NOT_VALID:
		currentDataMap.erase(return_pair.first);
		_throwError("invalid use of function in funcEval",line,pos0);
	      default:
		currentDataMap.erase(return_pair.first);
		_throwError("unexpected ierr in funcEval",line,pos0);
	      }

	    }

	  }
	  break;

	case CTI_DATA_NUMBER_TOKEN:

	  // ===========================================================
	  // make sure op and value stack are balanced...
	  // ===========================================================
	  {
	    if (valueStack.size() > opStack.size())
	      _throwError("expecting operator",line,pos0);
	    valueStack.push_back(pair<bool,CtiData*>(true,new CtiData()));
	    valueStack.back().second->new_d();
	    valueStack.back().second->d() = atof(token.c_str());
	  }
	  break;

	case CTI_DATA_STRING_TOKEN:
    {
      valueStack.push_back(pair<bool,CtiData*>(true,new CtiData()));
      valueStack.back().second->new_s();
      valueStack.back().second->s() = token.c_str();
    }
    break;
	  // ===========================================================
	  // this is a string, which should not have any associated ops
	  // or other values (i.e. you cannot operate on a string)...
	  // ===========================================================

	    /*
      {

	    // a string should not be part of a complex expression...

	    assert(valueStack.size() == 0);
	    assert(opStack.size() == 0);

	    // a string should only be followed by a comma, or a closing bracket, which
	    // will have been stripped off by the bracket pairing routine, so we check
	    // for a comma, or nothing, and return the string in v...

	    string token2;
	    CtiDataTokenType tt2 = _parseNextToken(token,pos,line);
	    if ((tt2 == CTI_DATA_COMMA_TOKEN)||(tt2 == CTI_DATA_NOTHING_TOKEN)) {
	    v.hash = CTI_DATA_STRING_BIT;
	    v.str = token;
	    return pos;
	    }
	    else {
	    _throwError("unsupported use of string",line,pos0);
	    }



	    assert(0);

	  }
    */
	  break;

	case CTI_DATA_OPERATOR_TOKEN:

	  // ===========================================================
	  // ===========================================================
	  {

	    // when we get a new operator, we are able to evaluate any of the previous ones
	    // that have the same or higher rank. This means that the value/op stack should
	    // only ever be > 2/1 if there are a series of operators of increaasing rank, e.g.
	    // expressions like
	    //
	    // 2+3/4**5
	    //
	    // have to be fully stacked before they can be evalulated

	    int this_op_rank = _getOpRank(token);
	    //cout << "parsed op: " << token << " with rank: " << this_op_rank << endl;
	    if (opStack.size() > 0) {
	      string top_op = opStack.back();
	      int top_op_rank = _getOpRank(top_op);
	      //cout << "compare to top operator: " << top_op << " with rank: " << top_op_rank << endl;
	      if (top_op_rank == this_op_rank) {
		// the top op is the same rank as this new op. To maintain L-to-R evaluation,
		// we process the two values and the top operator in the stack right now...
		opStack.pop_back();
		// use the top operator to act on the top 2 values in the valueStack, and
		// replace the top value with the answer...
		pair<bool,CtiData*> value1 = valueStack.back(); valueStack.pop_back();
		pair<bool,CtiData*> value0 = valueStack.back(); valueStack.pop_back();
		valueStack.push_back(pair<bool,CtiData*>(true,new CtiData()));
		if (_evalOperator(*(valueStack.back().second),value0.second,value1.second,top_op,b_eval_func) != CTI_DATA_OK) {
		  // delete the two popped values...
		  if (value0.first) delete value0.second;
		  if (value1.first) delete value1.second;
		  _throwError("_evalOperator error",line,pos0);
		}
		if (value0.first) delete value0.second;
		if (value1.first) delete value1.second;
	      }
	      else if (top_op_rank > this_op_rank) {
		// we need to process the whole stack and then push that evaluation
		// back on the stack...
		CtiData * value = new CtiData();
		try {
		  _processStack(*value,valueStack,opStack,b_eval_func);
		}
		catch (int e) {
		  delete value;
		  throw(0);
		}
		assert(valueStack.size() == 0);
		assert(opStack.size() == 0);
		valueStack.push_back(pair<bool,CtiData*>(true,value));
	      }
	    }
	    else if ((valueStack.size() == 0)&&((token=="-")||(token=="+"))) {
	      // both value and op stacks are empty, so this is a prefix negative
	      // or positive, so add a zero and treat the operator as binary...
	      valueStack.push_back(pair<bool,CtiData*>(true,new CtiData()));
	      valueStack.back().second->new_d();
	      valueStack.back().second->d() = 0.0;
	    }
	    if (valueStack.size() <= opStack.size())
	      _throwError("expecting value",line,pos0);
	    opStack.push_back(token);
	  }
	  break;

	default:
	  assert(0);
	  break;
	}

	// stacks are now...
	/*
	  cout << "\n\nvalueStack: " << endl;
	  for (int i = 0; i < valueStack.size(); ++i)
	  cout << " > " << *valueStack[i] << endl;
	  cout << "opStack: " << endl;
	  for (int i = 0; i < opStack.size(); ++i)
	  cout << " > " << opStack[i] << endl;
	  cout << "\n" << endl;
	*/

      }

    }
    catch (int e) {
      // any error in this routine needs to delete any undeleted members of the valueStack...
      for (unsigned int i = 0; i < valueStack.size(); ++i)
	if (valueStack[i].first)
	  delete valueStack[i].second;
      throw(0);
    }

    // logic problem
    assert(0);

  }

  CtiData * getUnregisteredCtiData(const string& expression, const bool b_eval_func) {

    // same as getCtiData, but do NOT look at registered data map...

    // NOTE: if you change this routine, you should check/change getUnregisteredCtiData
    // as well, which should be identical but not traverse the registered data...

    //cout << "getCtiData(\"" << expression << "\").." << endl;

    // start with the currentDataMap. It gets flushed every timestep, but holds
    // values for expressions that have just been evaluated...

    map<const string,CtiData>::iterator iter = currentDataMap.find(expression);
    if (iter != currentDataMap.end()) {

      //cout << "getCtiData: found \"" << expression << "\" in currentDataMap." << endl;
      //cout << "getCtiData: returning: " << &iter->second << endl;
      return &iter->second;

    }
    else {

      // behavior of the get*CtiData is to return NULL when an expression cannot be found or evaluated...

      try {

	// must need to evaluate this expression. Here we need to be careful with memory...
	//cout << "getCtiData: \"" << expression << "\" not found in registeredDataMap. must need to evaluate..." << endl;

	// construct the key-value pair and use this second's memory for the eval...
	pair<const string,CtiData> key_value_pair(expression,CtiData());
	_eval(key_value_pair.second,b_eval_func,expression);

	//cout << "getCtiData: back from _eval, about to push into map: " << key_value_pair.second << endl;

	// once an expression has been properly evaluated without throwing an exception, put it in the currentDataMap...
        const bool mem_flag = key_value_pair.second.getMemFlag();
	key_value_pair.second.setMemFlag(false);
	pair<map<const string,CtiData>::iterator,bool> return_pair = currentDataMap.insert(key_value_pair); //pair<const string,CtiData>(expression,data));
	key_value_pair.second.setMemFlag(mem_flag);

	//cout << "getCtiData: back from map, return_pair.second: " << return_pair.second << " data: " << key_value_pair.second << endl;

	if (return_pair.second) {
	  // the data was properly inserted, so this prevents data's destructor from cleaning up
	  // any memory when it goes out of scope...
	  //cout << "getCtiData: insert was good. Clear mem flag for: " << &key_value_pair.second << endl;
	  key_value_pair.second.setMemFlag(false);
          return_pair.first->second.setMemFlag(mem_flag);
	}
	else {

	  //cout << "getCtiData: insert FAILED. Leave mem flag " << key_value_pair.second.getMemFlag() << " for: " << &key_value_pair.second << endl;

	}

	//cout << "getCtiData: returning: " << &(return_pair.first->second) << endl;
	return &(return_pair.first->second); // return the ptr

      }
      catch(int ierr) {

	return NULL;

      }

    }

  }

  CtiData * getRegisteredCtiData(const string& expression,const bool verbose) { // verbose default == true

    // request registered data specifically

    //cout << "getRegisteredCtiData(\"" << expression << "\").." << endl;

    map<const string,CtiData>::iterator iter = registeredDataMap.find(expression);
    if (iter != registeredDataMap.end()) {

      //cout << "getRegistedCtiData: found \"" << expression << "\" in registeredDataMap." << endl;
      //cout << "getRegistedCtiData: returning: " << &iter->second << endl;
      return &iter->second;

    }
    else {

      if (verbose) {
	CWARN("registered data not found for \"" << expression << "\". Returning NULL.");
      }

      return NULL;

    }

  }

  CtiData * getCtiData(const string& expression, const bool b_eval_func) {

    // NOTE: if you change this routine, you should check/change getUnregisteredCtiData
    // as well, which should be identical but not traverse the registered data...

    //cout << "getCtiData(\"" << expression << "\").." << endl;

    // start with the currentDataMap. It gets flushed every timestep, but holds
    // values for expressions that have just been evaluated...

    map<const string,CtiData>::iterator iter = currentDataMap.find(expression);
    if (iter != currentDataMap.end()) {

      //cout << "getCtiData: found \"" << expression << "\" in currentDataMap." << endl;
      //cout << "getCtiData: returning: " << &iter->second << endl;
      return &iter->second;

    }
    else {

      //cout << "getCtiData: \"" << expression << "\" not found in currentDataMap. Checking registeredDataMap..." << endl;

      // next try the registered data map...

      iter = registeredDataMap.find(expression);
      if (iter != registeredDataMap.end()) {

	//cout << "getCtiData: found \"" << expression << "\" in registeredDataMap." << endl;
	//cout << "getCtiData: returning: " << &iter->second << endl;
	return &iter->second;

      }
      else {

	// behavior of the getCtiData is to return NULL when an expression cannot be found or evaluated...

	try {

          // must need to evaluate this expression. Here we need to be careful with memory...
         // cout << "getCtiData: \"" << expression << "\" not found in registeredDataMap. must need to evaluate..." << endl;

          // construct the key-value pair and use this second's memory for the eval...
          pair<const string,CtiData> key_value_pair(expression,CtiData());
          _eval(key_value_pair.second,b_eval_func,expression);

         // cout << "getCtiData: back from _eval, about to push into map: " << key_value_pair.second << endl;

          // once an expression has been properly evaluated without throwing an exception, put it in the currentDataMap...
          // we have to do this to prevent the extra copies' destructors from deleting the data. clang compilers did this.
          // TODO: maybe go to pointers to CtiData
          const bool mem_flag = key_value_pair.second.getMemFlag();
          key_value_pair.second.setMemFlag(false);
	  pair<map<const string,CtiData>::iterator,bool> return_pair = currentDataMap.insert(key_value_pair); //pair<const string,CtiData>(expression,data));
	  key_value_pair.second.setMemFlag(mem_flag);

         // cout << "getCtiData: back from map, return_pair.second: " << return_pair.second << " data: " << key_value_pair.second << endl;

          if (return_pair.second) {
            // the data was properly inserted, so this prevents data's destructor from cleaning up
            // any memory when it goes out of scope...
           // cout << "getCtiData: insert was good. Clear mem flag for: " << &key_value_pair.second << endl;
            key_value_pair.second.setMemFlag(false);
            return_pair.first->second.setMemFlag(mem_flag);

           // cout << "getCtiData: back from map, return_pair.second: " << return_pair.second << " data: " << key_value_pair.second << endl;
          }
          else {

           // cout << "getCtiData: insert FAILED. Leave mem flag " << key_value_pair.second.getMemFlag() << " for: " << &key_value_pair.second << endl;

          }

         // cout << "getCtiData: returning: " << &(return_pair.first->second) << endl;
          return &(return_pair.first->second); // return the ptr

	}
	catch(int ierr) {

	  return NULL;

	}

      }

    }

  }

  bool isRegisteredData(CtiData * data) {

    for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter)
      if (&iter->second == data)
	return true;

    return false;

  }

  // ==================================
  // defines...
  // ==================================

  /*
  void addDefine(const string& key,const string& value) {

    map<const string,string>::iterator iter = defineMap.find(key);
    if (iter != defineMap.end()) {
      if (mpi_rank == 0) cout << "CtiRegister::addDefine: modifying \"" << key << "\" from \"" << iter->second << "\" to \"" << value << "\"" << endl;
      iter->second = value;
    }
    else {
      if (mpi_rank == 0) cout << "CtiRegister::addDefine: \"" << key << "\" = \"" << value << "\"" << endl;
      defineMap[key] = value;
    }

  }
  */

  // =================================
  // interacting with registered data...
  // =================================

  void setRegisteredCvDnAndDn3Names(vector<string>& nameVec) {

    assert(nameVec.empty());

    for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter) {
      if ( ((iter->second.getType() == DN_DATA)||(iter->second.getType() == DN3_DATA)) && (iter->second.getUnindexedTopology() == CV_DATA) ) {
	nameVec.push_back(iter->first);
      }
    }

  }

  void setRegisteredDAndINames(vector<string>& nameVec) {

    assert(nameVec.empty());

    for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter) {
      if ((iter->second.getType() == D_DATA)||(iter->second.getType() == I_DATA)) {
	nameVec.push_back(iter->first);
      }
    }

  }

  void setRegisteredCvAndNoDNNames(vector<string>& nameVec) {

    assert(nameVec.empty());

    for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter) {
      if ( (iter->second.getType() == DN_DATA) && ((iter->second.getUnindexedTopology() == CV_DATA)||(iter->second.getUnindexedTopology() == NO_DATA)) ) {
	nameVec.push_back(iter->first);
      }
    }

  }
  void setRegisteredCvAndNoDN3Names(vector<string>& nameVec) {

    assert(nameVec.empty());

    for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter) {
      if ( (iter->second.getType() == DN3_DATA) && ((iter->second.getUnindexedTopology() == CV_DATA)||(iter->second.getUnindexedTopology() == NO_DATA)) ) {
	nameVec.push_back(iter->first);
      }
    }

  }

  void setRegisteredBfDN3Names(vector<string>& nameVec) {

    assert(nameVec.empty());

    for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter) {
      if ( (iter->second.getType() == DN3_DATA) && (iter->second.getUnindexedTopology() == BF_DATA) ) {
	nameVec.push_back(iter->first);
      }
    }

  }
  void setRegisteredBfDNNames(vector<string>& nameVec) {

    assert(nameVec.empty());

    for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter) {
      if ( (iter->second.getType() == DN_DATA) && (iter->second.getUnindexedTopology() == BF_DATA) )  {
	nameVec.push_back(iter->first);
      }
    }

  }

  void setRegisteredCvBfAndNoDNNames(vector<string>& nameVec) {

    assert(nameVec.empty());

    for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter) {
      if ( (iter->second.getType() == DN_DATA) &&
           ((iter->second.getUnindexedTopology() == CV_DATA)||(iter->second.getUnindexedTopology() == NO_DATA)||(iter->second.getUnindexedTopology() == BF_DATA)) ) {
	nameVec.push_back(iter->first);
      }
    }

  }
  void setRegisteredCvBfAndNoDN3Names(vector<string>& nameVec) {

    assert(nameVec.empty());

    for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter) {
      if ( (iter->second.getType() == DN3_DATA) &&
           ((iter->second.getUnindexedTopology() == CV_DATA)||(iter->second.getUnindexedTopology() == NO_DATA)||(iter->second.getUnindexedTopology() == BF_DATA)) ) {
	nameVec.push_back(iter->first);
      }
    }

  }

  void setRegisteredFuncEvalNames(vector<string>& nameVec) {
    for (list<CtiDataProducer*>::iterator iter = ctiDataProducerList.begin(); iter != ctiDataProducerList.end(); ++iter) {
      nameVec.insert(nameVec.end(),(*iter)->funcEvalList.begin(),(*iter)->funcEvalList.end());
    }
  }

  // ===========================================
  // int value I
  // ===========================================

  void _registerData(int& val,const string& name,const uint rw_bits) {

    if (mpi_rank == 0) cout << " > CtiRegister: int \"" << name << "\"" << endl;

    // map it...

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);

    // uniquesness check...

    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    // build the crd...

    CtiData * data = &(return_pair.first->second);
    data->setDataPtr(val);
    data->setBits(rw_bits);
    data->setType(I_DATA);

  }

  // ===================================================================================
  // int scalar IN registration
  // ===================================================================================

  void _registerData(int *&val,const string& name,const int topology,const uint rw_bits,int &n) {

    if (mpi_rank == 0) cout << " > CtiRegister: int array \"" << name << "\" without io" << endl;

    // map it...

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);

    // uniquesness check...

    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    // build the crd...

    CtiData * data = &(return_pair.first->second);
    data->setDataPtr(val);
    data->setStride(1);
    data->setBits(rw_bits);
    data->setType(IN_DATA);
    data->setTopology(topology);
    data->setSizePtr(n);

    // ensure that dde is provided whenever i/o is needed...

    assert(!( (rw_bits & READ_DATA) || (rw_bits & WRITE_DATA) || (rw_bits & CAN_WRITE_DATA) ));

  }
  void _registerData(int *&val,const string& name,const int topology,const uint rw_bits,int &n,int8 *&xora,DistributedDataExchanger *&dde) {

    if (mpi_rank == 0) cout << " > CtiRegister: int array \"" << name << "\" with io" << endl;

    // map it...

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);

    // uniquesness check...

    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    // build the crd...

    CtiData * data = &(return_pair.first->second);
    data->setDataPtr(val);
    data->setStride(1);
    data->setBits(rw_bits);
    data->setType(IN_DATA);
    data->setTopology(topology);
    data->setSizePtr(n);

    // dde consistency check...

    assert((rw_bits & READ_DATA) || (rw_bits & WRITE_DATA) || (rw_bits & CAN_WRITE_DATA));
    //assert(xora);
    //assert(dde);
    //assert(n == dde->my_nda);
    data->setDdeStuff(xora,dde);

  }

  // ==================================
  // data registration...
  // double value D
  // ==================================

  void _registerData(double& val,const string& name,const uint rw_bits) {

    if (mpi_rank == 0) cout << " > CtiRegister: double \"" << name << "\"" << endl;

    // map it...

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);

    // uniquesness check...

    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    // build the crd...

    CtiData * data = &(return_pair.first->second);
    data->setDataPtr(val);
    data->setBits(rw_bits);
    data->setType(D_DATA);

  }

  // ==================================
  // data registration...
  // double[3] value D3
  // ==================================

  void _registerData(double val[3],const string& name,const uint rw_bits) {

    if (mpi_rank == 0) cout << " > CtiRegister: double[3] \"" << name << "\"" << endl;

    // map it...

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);

    // uniquesness check...

    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    // build the crd...

    CtiData * data = &(return_pair.first->second);
    data->setDataPtr(val);
    data->setBits(rw_bits);
    data->setType(D3_DATA);

  }

  // ===================================================================================
  // double scalar DN registration
  // ===================================================================================

  void _registerData(double *&val,const string& name,const int topology,const uint rw_bits,int &n) {

    if (mpi_rank == 0) cout << " > CtiRegister: double array \"" << name << "\" without io X" << endl;

    // map it...

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);

    // uniquesness check...

    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    // build the crd...

    CtiData * data = &(return_pair.first->second);
    data->setDataPtr(val);
    data->setStride(2);
    data->setBits(rw_bits);
    data->setType(DN_DATA);
    data->setTopology(topology);
    data->setSizePtr(n);

    // ensure that dde is provided whenever i/o is needed...

    assert(!( (rw_bits & READ_DATA) || (rw_bits & WRITE_DATA) || (rw_bits & CAN_WRITE_DATA) ));

    // n-map stuff...

    map<const int*,NData>::iterator it = nDataMap.find(&n);
    if (it == nDataMap.end()) {
      pair<const int*,NData> key_value_pair = pair<const int*,NData>(&n,NData());
      pair<map<const int*,NData>::iterator,bool> return_pair = nDataMap.insert(key_value_pair);
      assert(return_pair.second);
      it = return_pair.first;
    }
    if (rw_bits & L_DATA) {
      assert(it->second.l_ptr == NULL);
      it->second.l_ptr = &val;
      if (mpi_rank == 0) cout << " > setting this data as lengthscale L" << endl;
    }

  }

  void _registerData(double *&val,const string& name,const int topology,const uint rw_bits,int& n,int8 * &xora,DistributedDataExchanger * &dde) {

    if (mpi_rank == 0) cout << " > CtiRegister: double array \"" << name << "\" with io" << endl;

    // map it...

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);

    // uniquesness check...

    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    // build the crd...

    CtiData * data = &(return_pair.first->second);
    data->setDataPtr(val);
    data->setStride(2);
    data->setBits(rw_bits);
    data->setType(DN_DATA);
    data->setTopology(topology);
    data->setSizePtr(n);

    // dde consistency check...

    assert((rw_bits & READ_DATA) || (rw_bits & WRITE_DATA) || (rw_bits & CAN_WRITE_DATA));
    //assert(xora);
    //assert(dde);
    //assert(n == dde->my_nda);
    data->setDdeStuff(xora,dde);

  }

  // trying to eliminate the need for topology for these registrations...
  /*
  void _registerData(double *&val,const string& name,const uint rw_bits,int &n) {
    // for now, just use NO_DATA...
    _registerData(val,name,NO_DATA,rw_bits,n);
  }
  */

  // ===================================================================================
  // double vector DN3 registration
  // ===================================================================================

  void _registerData(double (*&val)[3],const string& name,const int topology,const uint rw_bits,int &n) {

    if (mpi_rank == 0) cout << " > CtiRegister: double[3] array \"" << name << "\" without io X" << endl;

    // map it...

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);

    // uniqueness check...

    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    CtiData * data = &(return_pair.first->second);
    data->setDataPtr(val);
    data->setStride(6);
    data->setBits(rw_bits);
    data->setType(DN3_DATA);
    data->setTopology(topology);
    data->setSizePtr(n);

    // ensure that dde is provided whenever i/o is needed...

    assert(!( (rw_bits & READ_DATA) || (rw_bits & WRITE_DATA) || (rw_bits & CAN_WRITE_DATA) ));

    // n-map stuff...

    map<const int*,NData>::iterator it = nDataMap.find(&n);
    if (it == nDataMap.end()) {
      pair<const int*,NData> key_value_pair = pair<const int*,NData>(&n,NData());
      pair<map<const int*,NData>::iterator,bool> return_pair = nDataMap.insert(key_value_pair);
      assert(return_pair.second);
      it = return_pair.first;
    }
    if (rw_bits & X_DATA) {
      assert(it->second.x_ptr == NULL);
      it->second.x_ptr = &val;
      if (mpi_rank == 0) cout << " > setting this data as position X" << endl;
    }

  }
  void _registerData(double (*&val)[3],const string& name,const int topology,const uint rw_bits,int& n,int8 * &xora,DistributedDataExchanger * &dde) {

    if (mpi_rank == 0) cout << " > CtiRegister: double[3] array \"" << name << "\" with io" << endl;

    // map it...

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);

    // uniqueness check...

    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    CtiData * data = &(return_pair.first->second);
    data->setDataPtr(val);
    data->setStride(6);
    data->setBits(rw_bits);
    data->setType(DN3_DATA);
    data->setTopology(topology);
    data->setSizePtr(n);

    // dde consistency check...

    assert((rw_bits & READ_DATA) || (rw_bits & WRITE_DATA) || (rw_bits & CAN_WRITE_DATA));
    //assert(xora);
    //assert(dde);
    //assert(n == dde->my_nda);
    data->setDdeStuff(xora,dde);

    // n-map stuff...

    map<const int*,NData>::iterator it = nDataMap.find(&n);
    if (it == nDataMap.end()) {
      pair<const int*,NData> key_value_pair = pair<const int*,NData>(&n,NData());
      pair<map<const int*,NData>::iterator,bool> return_pair = nDataMap.insert(key_value_pair);
      assert(return_pair.second);
      it = return_pair.first;
    }
    if (rw_bits & X_DATA) {
      assert(it->second.x_ptr == NULL);
      it->second.x_ptr = &val;
      if (mpi_rank == 0) cout << " > setting this data as position X" << endl;
    }

  }

  // trying to eliminate the need for topology for these registrations...
  /*
  void _registerData(double (*&val)[3],const string& name,const uint rw_bits,int &n) {
    // for now, just use NO_DATA...
    _registerData(val,name,NO_DATA,rw_bits,n);
  }
  */

  // ============================================================================
  // name-based registration where CtiData manages memory...
  // ============================================================================

  // I_DATA, D_DATA and D3_DATA
  void _registerData(const string& name,const int datatype,const uint rw_bits) {

    // this routine is used to register data when the calling process does not manage memory,
    // so CtiRegister manages the memory in the CtiData that lives in the registeredDataMap....

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);
    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    CtiData * data = &(return_pair.first->second);
    assert(data->empty());

    if (datatype == I_DATA) {

      if (mpi_rank == 0) cout << " > CtiRegister: int \"" << name << "\"";
      data->new_i();

    }
    else if (datatype == D_DATA) {

      if (mpi_rank == 0) cout << " > CtiRegister: double \"" << name << "\"";
      data->new_d();

    }
    else if (datatype == D3_DATA) {

      if (mpi_rank == 0) cout << " > CtiRegister: double[3] \"" << name << "\"";
      data->new_d3();

    }
    else {

      CERR("unsupported datatype: " << datatype << " for name \"" << name << "\"");

    }

    // check rw_bits stuff...

    if ((rw_bits & READ_DATA) || (rw_bits & WRITE_DATA) || (rw_bits & CAN_WRITE_DATA)) {
      if (mpi_rank == 0) cout << " with io" << endl;
      data->setBits(rw_bits);
    }
    else {
      if (mpi_rank == 0) cout << " without io" << endl;
      assert(rw_bits == 0);
    }

  }

  // IN_DATA, DN_DATA and DN3_DATA
  void _registerData(const string& name,const int topology,const int datatype,const uint rw_bits,int &n) {
    // this routine is used to register data when the calling process does not manage memory,
    // so CtiRegister manages the memory in the CtiData that lives in the registeredDataMap....

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);
    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    CtiData * data = &(return_pair.first->second);
    assert(data->empty());

    if (datatype == IN_DATA) {

      if (mpi_rank == 0) cout << " > CtiRegister: int array \"" << name << "\"";
      data->new_in(topology,n);

    }
    else if (datatype == DN_DATA) {

      if (mpi_rank == 0) cout << " > CtiRegister: double array \"" << name << "\"";
      data->new_dn(topology,n);

    }
    else if (datatype == DN3_DATA) {

      if (mpi_rank == 0) cout << " > CtiRegister: double[3] array \"" << name << "\"";
      data->new_dn3(topology,n);

    }
    else {

      CERR("unsupported datatype: " << datatype << " for name \"" << name << "\"");

    }

    // make sure we don't need io...

    assert(!((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA)));

  }

  void _registerData(const string& name,const int topology,const int datatype,const uint rw_bits,int &n,int8 * &xora,DistributedDataExchanger * &dde) {

    // this routine is used to register data when the calling process does not manage memory,
    // so CtiRegister manages the memory in the CtiData that lives in the registeredDataMap....

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);
    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    CtiData * data = &(return_pair.first->second);
    assert(data->empty());

    if (datatype == IN_DATA) {

      if (mpi_rank == 0) cout << " > CtiRegister: int array \"" << name << "\"";
      //data->new_in(topology,n);
      data->setStride(1);

    }
    else if (datatype == DN_DATA) {

      if (mpi_rank == 0) cout << " > CtiRegister: double array \"" << name << "\"";
      //data->new_dn(topology,n);
      data->setStride(2);

    }
    else if (datatype == DN3_DATA) {

      if (mpi_rank == 0) cout << " > CtiRegister: double[3] array \"" << name << "\"";
      //data->new_dn3(topology,n);
      data->setStride(6);

    }
    else {

      CERR("unsupported datatype: " << datatype << " for name \"" << name << "\"");

    }

    // check rw_bits stuff...
    assert((rw_bits & READ_DATA)||(rw_bits & WRITE_DATA)||(rw_bits & CAN_WRITE_DATA));
    if (mpi_rank == 0) cout << " with io" << endl;

    // here we just set everything, but we DO NOT allocate yet, because the size has not been set...
    data->setType(datatype);
    data->setTopology(topology);
    data->setBits(rw_bits);
    data->setDdeStuff(xora,dde);
    data->setMemFlag(true);
    data->setAccessFlag(true);
    data->setSizePtr(n);
  }
  void _initData() {
    COUT1("CtiRegister::_initData()");
    CtiRegister::_initStats(); // must be done separately so size_ptr can get its own value (to no screw with the ptr it referenced)
    for (map<const string,CtiData>::iterator it = registeredDataMap.begin(); it != registeredDataMap.end(); ++it) {
      // if registration owns the memory and it hasn't allocated it yet (because the size wasn't known)...
      if ((it->second.getMemFlag())&&(it->second.getDataPtr() == NULL)) {
        it->second.initData();
      }
    }
  }

  void _registerGlobalIndex(int8 * &index,int &n) {
    // for now, we insist that this call is not the first call to
    // registered data of a certain "n"...
    map<const int*,NData>::iterator it = nDataMap.find(&n);
    assert(it != nDataMap.end()); // expect data of this "n" to exist
    assert(it->second.index_ptr == NULL); // should set the index only once per "n" data type
    it->second.index_ptr = &index;
  }

  void _registerName(const string& name,int &n) {
    // for now, we insist that this call is not the first call to
    // registered data of a certain "n"...
    map<const int*,NData>::iterator it = nDataMap.find(&n);
    assert(it != nDataMap.end()); // expect data of this "n" to exist
    assert(it->second.name == "UNKNOWN"); // should set the name only once per "n" data type
    it->second.name = name;
  }

  void _registerNoote(int (*&noote)[4],int &nte,int &n) {
    // NODE_OF_TET...
    // for now, we insist that this call is not the first call to
    // registered data of a certain "n"...
    map<const int*,NData>::iterator it = nDataMap.find(&n);
    assert(it != nDataMap.end()); // expect data of this "n" to exist
    assert(it->second.noote_ptr == NULL);
    it->second.noote_ptr = &noote;
    assert(it->second.nte_ptr == NULL);
    it->second.nte_ptr = &nte;
  }

  string getFirstSpecifier(const string var,const bool filter_protected) {
    const size_t colon_pos = var.find_first_of(":");
    if (filter_protected) {
      const string delimiters = " =+-/*(^%$,";  // protected chars in expression parsing
      const size_t start_pos = var.substr(0,colon_pos).find_last_of(delimiters);
      return var.substr(start_pos+1,colon_pos-start_pos-1);
    }
    else {
      return var.substr(0,colon_pos);
    }
  }

  string getSecondSpecifier(const string var) {
    const size_t colon_pos = var.find_first_of(":");
    return var.substr(colon_pos+1,var.size()-colon_pos);
  }

  int getDataFlag(const string& name) {

    map<const string,CtiData>::iterator iter = registeredDataMap.find(name);
    if (iter == registeredDataMap.end()) {
      CERR("data name not found in registeredDataMap: \"" << name << "\"");
    }
    return iter->second.getFlag();

  }

  void setDataFlag(const string& name,const int val) {

    map<const string,CtiData>::iterator iter = registeredDataMap.find(name);
    if (iter == registeredDataMap.end()) {
      CERR("data name not found in registeredDataMap: \"" << name << "\"");
    }
    iter->second.setFlag(val);

  }

  bool checkDataFlag(const string& name) {

    return getDataFlag(name) != 0;

  }

  void clearAllDataFlags() {

    for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter) {
      iter->second.setFlag(0);
    }

  }

  void readData(const string& filename) {

    COUT1("CtiRegister::readData: " << filename);

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

    bool b_foundHash = false;
    MPI_Offset offset = 8; // 2 ints
    Header header;
    map<const string,CtiData>::iterator iter;
    int done = 0;
    int record_count = 0;
    while (done != 1) {

      if (mpi_rank == 0) {
	MPI_File_read_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
	if (byte_swap) ByteSwap::byteSwapHeader(&header,1);
      }
      MPI_Bcast(&header,1,MPI_Header,0,mpi_comm);

      switch (header.id) {

      case UGP_IO_D0:

	if (mpi_rank == 0) {
	  cout << " > D0 \"" << header.name << "\" (=" << header.rdata[0] << ")";
	  cout.flush();
	}

	iter = registeredDataMap.find(header.name);
	if ( (iter == registeredDataMap.end()) || (iter->second.getType() != D_DATA) || (!iter->second.checkBit(READ_DATA)) ) {
	  if (mpi_rank == 0) cout << " skipping." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " setting." << endl;
	  iter->second.d() = header.rdata[0];
	  iter->second.setFlag(1);
	}

	break;

      case UGP_IO_I0:

	if (mpi_rank == 0) {
	  cout << " > I0 \"" << header.name << "\" (=" << header.idata[0] << ")";
	  cout.flush();
	}

	iter = registeredDataMap.find(header.name);
	if ( (iter == registeredDataMap.end()) || (iter->second.getType() != I_DATA) || (!iter->second.checkBit(READ_DATA)) ) {
	  if (mpi_rank == 0) cout << " skipping." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " setting." << endl;
	  iter->second.i() = header.idata[0];
	  iter->second.setFlag(1);
	}

	break;

      case UGP_IO_CV_I1:

	if (mpi_rank == 0) {
	  cout << " > CV_I1 \"" << header.name << "\"";
	  cout.flush();
	}

	iter = registeredDataMap.find(header.name);
	if ( (iter == registeredDataMap.end()) || (iter->second.getType() != IN_DATA) ||
             (iter->second.getUnindexedTopology() != CV_DATA) || (!iter->second.checkBit(READ_DATA)) ) {
	  if (mpi_rank == 0) cout << " skipping." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " setting." << endl;
	  iter->second.readData(fh,header,offset,byte_swap);
	  iter->second.setFlag(1);
	}

	break;

      case UGP_IO_BF_I1:

	if (mpi_rank == 0) {
	  cout << " > BF_I1 \"" << header.name << "\"";
	  cout.flush();
	}

	iter = registeredDataMap.find(header.name);
	if ( (iter == registeredDataMap.end()) || (iter->second.getType() != IN_DATA) ||
             (iter->second.getUnindexedTopology() != BF_DATA) || (!iter->second.checkBit(READ_DATA)) ) {
	  if (mpi_rank == 0) cout << " skipping." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " setting." << endl;
	  iter->second.readData(fh,header,offset,byte_swap);
	  iter->second.setFlag(1);
	}

	break;

      case UGP_IO_FA_I1:

	if (mpi_rank == 0) {
	  cout << " > FA_I1 \"" << header.name << "\"";
	  cout.flush();
	}

	iter = registeredDataMap.find(header.name);
	if ( (iter == registeredDataMap.end()) || (iter->second.getType() != IN_DATA) ||
             (iter->second.getUnindexedTopology() != FA_DATA) || (!iter->second.checkBit(READ_DATA)) ) {
	  if (mpi_rank == 0) cout << " skipping." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " setting." << endl;
	  iter->second.readData(fh,header,offset,byte_swap);
	  iter->second.setFlag(1);
	}

	break;

      case UGP_IO_CV_D1:

	if (mpi_rank == 0) {
	  cout << " > CV_D1 \"" << header.name << "\"";
	  cout.flush();
	}

	iter = registeredDataMap.find(header.name);
	if ( (iter == registeredDataMap.end()) || (iter->second.getType() != DN_DATA) ||
             (iter->second.getUnindexedTopology() != CV_DATA) || (!iter->second.checkBit(READ_DATA)) ) {
	  if (mpi_rank == 0) cout << " skipping." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " setting." << endl;
	  iter->second.readData(fh,header,offset,byte_swap);
	  iter->second.setFlag(1);
	}

	break;

      case UGP_IO_BF_D1:

	if (mpi_rank == 0) {
	  cout << " > BF_D1 \"" << header.name << "\"";
	  cout.flush();
	}

	iter = registeredDataMap.find(header.name);
	if ( (iter == registeredDataMap.end()) || (iter->second.getType() != DN_DATA) ||
             (iter->second.getUnindexedTopology() != BF_DATA) || (!iter->second.checkBit(READ_DATA)) ) {
	  if (mpi_rank == 0) cout << " skipping." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " setting." << endl;
	  iter->second.readData(fh,header,offset,byte_swap);
	  iter->second.setFlag(1);
	}

	break;

      case UGP_IO_FA_D1:

	if (mpi_rank == 0) {
	  cout << " > FA_D1 \"" << header.name << "\"";
	  cout.flush();
	}

	iter = registeredDataMap.find(header.name);
	if ( (iter == registeredDataMap.end()) || (iter->second.getType() != DN_DATA) ||
             (iter->second.getUnindexedTopology() != FA_DATA) || (!iter->second.checkBit(READ_DATA)) ) {
	  if (mpi_rank == 0) cout << " skipping." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " setting." << endl;
	  iter->second.readData(fh,header,offset,byte_swap);
	  iter->second.setFlag(1);
	}

	break;

      case UGP_IO_CV_D2:

	if (mpi_rank == 0) {
	  cout << " > CV_D2 \"" << header.name << "\"";
	  cout.flush();
	}

	iter = registeredDataMap.find(header.name);
	if ( (iter == registeredDataMap.end()) || (iter->second.getType() != DN3_DATA) ||
             (iter->second.getUnindexedTopology() != CV_DATA) || (!iter->second.checkBit(READ_DATA)) ) {
	  if (mpi_rank == 0) cout << " skipping." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " setting." << endl;
	  iter->second.readData(fh,header,offset,byte_swap);
	  iter->second.setFlag(1);
	}

	break;

      case UGP_IO_BF_D2:

	if (mpi_rank == 0) {
	  cout << " > BF_D2 \"" << header.name << "\"";
	  cout.flush();
	}

	iter = registeredDataMap.find(header.name);
	if ( (iter == registeredDataMap.end()) || (iter->second.getType() != DN3_DATA) ||
             (iter->second.getUnindexedTopology() != BF_DATA) || (!iter->second.checkBit(READ_DATA)) ) {
	  if (mpi_rank == 0) cout << " skipping." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " setting." << endl;
	  iter->second.readData(fh,header,offset,byte_swap);
	  iter->second.setFlag(1);
	}

	break;

      case UGP_IO_FA_D2:

	if (mpi_rank == 0) {
	  cout << " > FA_D2 \"" << header.name << "\"";
	  cout.flush();
	}

	iter = registeredDataMap.find(header.name);
	if ( (iter == registeredDataMap.end()) || (iter->second.getType() != DN3_DATA) ||
             (iter->second.getUnindexedTopology() != FA_DATA) || (!iter->second.checkBit(READ_DATA)) ) {
	  if (mpi_rank == 0) cout << " skipping." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " setting." << endl;
	  iter->second.readData(fh,header,offset,byte_swap);
	  iter->second.setFlag(1);
	}

        break;

      case UGP_IO_HASHIDS:

        b_foundHash = true;
        int slesNonMatchFlag;
        //read two hash id's from sles file. First id identifies
        //the current file, store in slesHash.  Second  id
        //identifies the "parent" mles file.  Store in mlesHash
        if (mpi_rank==0)  cout << "RestartHashUtilities::slesReadHashes()" << endl;
        RestartHashUtilities::slesReadHashes(fh, offset, header); //will bCast to all ranks
        if (mpi_rank==0){
          slesNonMatchFlag = RestartHashUtilities::slesConsistencyCheck(); //-1 not enough info, 0 match, 1 no match
          cout << " > Found sles hash: " << RestartHashUtilities::slesHash << endl;
          if (slesNonMatchFlag==0){
            cout << " >  with mles hash: " << RestartHashUtilities::mlesHash << endl;
          }
          else{
            cout << " > sles expects mles hash: " << RestartHashUtilities::myHash << endl;
            cout << " >      current mles hash: " << RestartHashUtilities::mlesHash << endl;
          }
        }
        MPI_Bcast(&slesNonMatchFlag,1,MPI_INT,0,mpi_comm);
        if (slesNonMatchFlag>0){
	  CERR("sles data file is not a match for the existing mesh.\n" <<
	       "Consider using INTERP_FROM to initialize with this data.");
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

      ++record_count;
      if (record_count > 1000) {
	CERR("data file record count exceeds 1000. file is probably corrupt: \"" << filename << "\"");
      }

    }

    MPI_File_close(&fh);

    if (!b_foundHash) { //set a default sles hash if one was not found
      std::stringstream ss;
      ss << "Missing sles Hash";
      RestartHashUtilities::slesHash.init(ss,RestartHashUtilities::sha1hashlength);
      if (mpi_rank==0) cout << " > setting default sles hash id " <<  RestartHashUtilities::slesHash << endl;
    }

  }

  void readLpData(const string& filename) {

    COUT1("CtiRegister::readLpData: " << filename);

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
    map<const string,CtiData>::iterator iter;
    int done = 0;
    int record_count = 0;
    while (done != 1) {

      if (mpi_rank == 0) {
	MPI_File_read_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
	if (byte_swap) ByteSwap::byteSwapHeader(&header,1);
      }
      MPI_Bcast(&header,1,MPI_Header,0,mpi_comm);

      switch (header.id) {

      case UGP_IO_LP_I1:

	if (mpi_rank == 0) {
	  cout << " > LP_I1 \"" << header.name << "\"";
	  cout.flush();
	}

	iter = registeredDataMap.find(header.name);
	if ( (iter == registeredDataMap.end()) || (iter->second.getType() != IN_DATA) ||
             (iter->second.getUnindexedTopology() != LP_DATA) || (!iter->second.checkBit(READ_DATA)) ) {
	  if (mpi_rank == 0) cout << " skipping." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " setting." << endl;
	  iter->second.readData(fh,header,offset,byte_swap);
	  iter->second.setFlag(1);
	}

	break;

      case UGP_IO_LP_D1:

	if (mpi_rank == 0) {
	  cout << " > LP_D1 \"" << header.name << "\"";
	  cout.flush();
	}

	iter = registeredDataMap.find(header.name);
        if (iter == registeredDataMap.end()) {
          cout << " cannot find " << endl;

        }
	if ( (iter == registeredDataMap.end()) || (iter->second.getType() != DN_DATA) ||
             (iter->second.getUnindexedTopology() != LP_DATA) || (!iter->second.checkBit(READ_DATA)) ) {
	  if (mpi_rank == 0) cout << " skipping." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " setting." << endl;
	  iter->second.readData(fh,header,offset,byte_swap);
	  iter->second.setFlag(1);
	}

	break;

      case UGP_IO_LP_D2:

	if (mpi_rank == 0) {
	  cout << " > LP_D2 \"" << header.name << "\"";
	  cout.flush();
	}

	iter = registeredDataMap.find(header.name);
	if ( (iter == registeredDataMap.end()) || (iter->second.getType() != DN3_DATA) ||
             (iter->second.getUnindexedTopology() != LP_DATA) || (!iter->second.checkBit(READ_DATA)) ) {
	  if (mpi_rank == 0) cout << " skipping." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " setting." << endl;
	  iter->second.readData(fh,header,offset,byte_swap);
	  iter->second.setFlag(1);
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

      ++record_count;
      if (record_count > 1000) {
	CERR("data file record count exceeds 1000. file is probably corrupt: \"" << filename << "\"");
      }

    }

    MPI_File_close(&fh);

  }

  void writeData(const string& filename) {

    //if (mpi_rank == 0) cout << "writeData: " << filename << endl;

    const double wtime0 = MPI_Wtime();

    char fname[128];
    sprintf(fname,"%s",filename.c_str());
    MPI_File_delete(fname,MPI_INFO_NULL);

    MPI_File fh;
    MPI_File_open(mpi_comm,fname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    if ( mpi_rank == 0 )  {
      int itmp[2] = {UGP_IO_MAGIC_NUMBER+1, 5}; // use a different magic number? -- to indicate file type as data, not restart
      cout << " > UGP_IO_VERSION: " << itmp[1] << endl;
      MPI_File_write(fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
    }

    MPI_Offset offset = int_size*2;

    if (mpi_rank == 0){
      //snapshot hash constructed from the mles hash and sles hash from which the calc began
      //includes all variable names as well as single valued variable values.
      stringstream ss;  ss << RestartHashUtilities::mlesHash << RestartHashUtilities::slesHash;
      for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter) {
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
    }
    offset += RestartHashUtilities::writeHashToRestartHeader(fh); //update offset on all ranks

    for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter) {
      if (iter->second.checkBit(WRITE_DATA)) {
	iter->second.writeData(iter->first,fh,offset);
      }
    }

    if ( mpi_rank == 0 ) {
      cout << " > EOF" << endl;
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

  }

  void writeData(MPI_File &fh, MPI_Offset &offset) {

    for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter) {
      if (iter->second.checkBit(WRITE_DATA)) {
        iter->second.writeData(iter->first,fh,offset);
      }
    }

  }

  // =======================================================================
  // stats infrastructure
  // =======================================================================

  vector<Stats> statsVec;

  void registerStats(const string& vname,CtiData * data,const bool b_init) {

    COUT1("CtiRegister::registerStats(): " << vname << " b_init: " << b_init);

    const int dtype = data->getType();
    if (dtype == D_DATA) {

      // make sure there are no conflicts with the registered data names....
      if ((registeredDataMap.find(vname+"_wgt") != registeredDataMap.end())||
      (registeredDataMap.find(vname+"_avg") != registeredDataMap.end())||
      (registeredDataMap.find(vname+"_rms") != registeredDataMap.end())) {
        CWARN("STATS for " << vname << " already exist.");
        return;
      }

      // stats on D_DATA requires wgt, avg and rms...
      statsVec.push_back(Stats(vname,D_DATA));
      Stats * stats = &statsVec.back();

      pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(pair<const string,CtiData>(vname+"_wgt",CtiData()));
      assert(return_pair.second); // make sure insert was successful
      stats->wgt = &return_pair.first->second;
      stats->wgt->new_d();
      stats->wgt->setBits(READWRITE_DATA);

      return_pair = registeredDataMap.insert(pair<const string,CtiData>(vname+"_avg",CtiData()));
      assert(return_pair.second); // make sure insert was successful
      stats->avg = &return_pair.first->second;
      stats->avg->new_d();
      stats->avg->setBits(READWRITE_DATA);

      return_pair = registeredDataMap.insert(pair<const string,CtiData>(vname+"_rms",CtiData()));
      assert(return_pair.second); // make sure insert was successful
      stats->rms = &return_pair.first->second;
      stats->rms->new_d();
      stats->rms->setBits(READWRITE_DATA);

      // if the source is registered data, then we take a reference directly to it. This
      // makes the processing of certain stats faster because we do not have to invoke
      // the data evaluation process...

      for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter) {
        if (&iter->second == data) {
          stats->source = data;
          break;
        }
      }

      // the option to init this particular stats is included to allow the
      // user to add stats to a running job.

      stats->wgt->d() = HUGE_VAL;  // default value to eval whether new stat or not
      if (b_init) {
        stats->wgt->d() = 0.0;
        stats->avg->d() = 0.0;
        stats->rms->d() = 0.0;
      }

    }
    else if (dtype == DN_DATA) {

      // make sure there are no conflicts with the registered data names....
      if ((registeredDataMap.find(vname+"_wgt") != registeredDataMap.end())||
	  (registeredDataMap.find(vname+"_avg") != registeredDataMap.end())||
	  (registeredDataMap.find(vname+"_rms") != registeredDataMap.end())) {
	CWARN("STATS for " << vname << " already exist.");
	return;
      }

      // stats on DN_DATA requires wgt, avg and rms...
      bool b_io = true; // as a HACK to support CHT, we allow for stats without io for now...
      if (!data->hasDdeStuff()) {
	if (mpi_rank == 0) cout << "WARNING: stats requested for var: " << vname << " will not support io" << endl;
	b_io = false;
      }
      statsVec.push_back(Stats(vname,DN_DATA));
      Stats * stats = &statsVec.back();

      pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(pair<const string,CtiData>(vname+"_wgt",CtiData()));
      assert(return_pair.second); // make sure insert was successful
      stats->wgt = &return_pair.first->second;
      stats->wgt->new_d();
      if (b_io) stats->wgt->setBits(READWRITE_DATA);

      return_pair = registeredDataMap.insert(pair<const string,CtiData>(vname+"_avg",CtiData()));
      assert(return_pair.second); // make sure insert was successful
      stats->avg = &return_pair.first->second;
      stats->avg->new_dn(*data,b_init);
      if (b_io) stats->avg->setBits(READWRITE_DATA);

      return_pair = registeredDataMap.insert(pair<const string,CtiData>(vname+"_rms",CtiData()));
      assert(return_pair.second); // make sure insert was successful
      stats->rms = &return_pair.first->second;
      stats->rms->new_dn(*data,b_init);
      if (b_io) stats->rms->setBits(READWRITE_DATA);

      for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter) {
        if (&iter->second == data) {
          stats->source = data;
          break;
        }
      }

      stats->wgt->d() = HUGE_VAL;  // default value to eval whether new stat or not
      if (b_init) {
        stats->wgt->d() = 0.0;
        for (int i = 0; i < stats->avg->size(); ++i) stats->avg->dn(i) = 0.0;
        for (int i = 0; i < stats->rms->size(); ++i) stats->rms->dn(i) = 0.0;
      }

    }
    else if (dtype == DN3_DATA) {

      // make sure there are no conflicts with the registered data names....
      if ((registeredDataMap.find(vname+"_wgt") != registeredDataMap.end())||
      (registeredDataMap.find(vname+"_avg") != registeredDataMap.end())||
      (registeredDataMap.find(vname+"_rms") != registeredDataMap.end())||
      (registeredDataMap.find(vname+"_rey") != registeredDataMap.end())) {
        CWARN("STATS for " << vname << " already exist.");
        return;
      }

      // stats on D2 requires wgt, avg and rms and rey...
      bool b_io = true; // as a HACK to support CHT, we allow for stats without io for now...
      if (!data->hasDdeStuff()) {
	if (mpi_rank == 0) cout << "WARNING: stats requested for var: " << vname << " will not support io" << endl;
	b_io = false;
      }
      statsVec.push_back(Stats(vname,DN3_DATA));
      Stats * stats = &statsVec.back();

      pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(pair<const string,CtiData>(vname+"_wgt",CtiData()));
      assert(return_pair.second); // make sure insert was successful
      stats->wgt = &return_pair.first->second;
      stats->wgt->new_d();
      if (b_io) stats->wgt->setBits(READWRITE_DATA);

      return_pair = registeredDataMap.insert(pair<const string,CtiData>(vname+"_avg",CtiData()));
      assert(return_pair.second); // make sure insert was successful
      stats->avg = &return_pair.first->second;
      stats->avg->new_dn3(*data,b_init);
      if (b_io) stats->avg->setBits(READWRITE_DATA);

      return_pair = registeredDataMap.insert(pair<const string,CtiData>(vname+"_rms",CtiData()));
      assert(return_pair.second); // make sure insert was successful
      stats->rms = &return_pair.first->second;
      stats->rms->new_dn3(*data,b_init);
      if (b_io) stats->rms->setBits(READWRITE_DATA);

      return_pair = registeredDataMap.insert(pair<const string,CtiData>(vname+"_rey",CtiData()));
      assert(return_pair.second); // make sure insert was successful
      stats->rey = &return_pair.first->second;
      stats->rey->new_dn3(*data,b_init);
      if (b_io) stats->rey->setBits(READWRITE_DATA);

      for (map<const string,CtiData>::iterator iter = registeredDataMap.begin(); iter != registeredDataMap.end(); ++iter) {
        if (&iter->second == data) {
          stats->source = data;
          break;
        }
      }

      stats->wgt->d() = HUGE_VAL;  // default value to eval whether new stat or not
      if (b_init) {
        stats->wgt->d() = 0.0;
        for (int i = 0; i < stats->avg->size(); ++i) FOR_J3 stats->avg->dn3(i,j) = 0.0;
        for (int i = 0; i < stats->rms->size(); ++i) FOR_J3 stats->rms->dn3(i,j) = 0.0;
        for (int i = 0; i < stats->rey->size(); ++i) FOR_J3 stats->rey->dn3(i,j) = 0.0;
      }

    }
    else {
      CERR("data type for \"" << vname << "\" not supported for stats: " << data->getTypeAsString());
    }

  }

  void _initStats() {

    for (int ii = 0,ii_end=statsVec.size(); ii < ii_end; ++ii) {
      const int dtype = statsVec[ii].datatype;
      if (dtype == D_DATA) {
        statsVec[ii].wgt->d() = 0.0;
        statsVec[ii].avg->d() = 0.0;
        statsVec[ii].rms->d() = 0.0;
      }
      else if (dtype == DN_DATA) {
        statsVec[ii].wgt->d() = 0.0;
        statsVec[ii].avg->alloc_dn();
        for (int i = 0; i < statsVec[ii].avg->size(); ++i) statsVec[ii].avg->dn(i) = 0.0;
        statsVec[ii].rms->alloc_dn();
        for (int i = 0; i < statsVec[ii].rms->size(); ++i) statsVec[ii].rms->dn(i) = 0.0;
      }
      else if (dtype == DN3_DATA) {
        statsVec[ii].wgt->d() = 0.0;
        statsVec[ii].avg->alloc_dn3();
        for (int i = 0; i < statsVec[ii].avg->size(); ++i) FOR_J3 statsVec[ii].avg->dn3(i,j) = 0.0;
        statsVec[ii].rms->alloc_dn3();
        for (int i = 0; i < statsVec[ii].rms->size(); ++i) FOR_J3 statsVec[ii].rms->dn3(i,j) = 0.0;
        statsVec[ii].rey->alloc_dn3();
        for (int i = 0; i < statsVec[ii].rey->size(); ++i) FOR_J3 statsVec[ii].rey->dn3(i,j) = 0.0;
      }
      else {
        assert(0);
      }
    }
  }

  void resetStats() {

    for (int ii = 0,ii_end=statsVec.size(); ii < ii_end; ++ii) {
      const int dtype = statsVec[ii].datatype;
      if (dtype == D_DATA) {
	statsVec[ii].wgt->d() = 0.0;
	statsVec[ii].avg->d() = 0.0;
	statsVec[ii].rms->d() = 0.0;
      }
      else if (dtype == DN_DATA) {
	statsVec[ii].wgt->d() = 0.0;
	for (int i = 0; i < statsVec[ii].avg->size(); ++i) statsVec[ii].avg->dn(i) = 0.0;
	for (int i = 0; i < statsVec[ii].rms->size(); ++i) statsVec[ii].rms->dn(i) = 0.0;
      }
      else if (dtype == DN3_DATA) {
	statsVec[ii].wgt->d() = 0.0;
	for (int i = 0; i < statsVec[ii].avg->size(); ++i) FOR_J3 statsVec[ii].avg->dn3(i,j) = 0.0;
	for (int i = 0; i < statsVec[ii].rms->size(); ++i) FOR_J3 statsVec[ii].rms->dn3(i,j) = 0.0;
	for (int i = 0; i < statsVec[ii].rey->size(); ++i) FOR_J3 statsVec[ii].rey->dn3(i,j) = 0.0;
      }
      else {
	assert(0);
      }
    }

  }

  void clearStatsBits(const uint _bits) {

    for (int ii =0,ii_end = statsVec.size(); ii < ii_end; ++ii) {

      assert( statsVec[ii].wgt != NULL);
      statsVec[ii].wgt->clearBits(_bits);

      if ( statsVec[ii].avg )
        statsVec[ii].avg->clearBits(_bits);

      if ( statsVec[ii].rms )
        statsVec[ii].rms->clearBits(_bits);

      if ( statsVec[ii].rey )
        statsVec[ii].rey->clearBits(_bits);

    }
  }

  void updateStats(const double wgt,const bool verbose) {

    // this can happen if you try to update the stats w/o setting the timestep...
    if (wgt == 0.0) return;

    //COUT1("CtiRegister::updateStats: wgt=" << wgt);

    for (int ii = 0,ii_end=statsVec.size(); ii < ii_end; ++ii) {

      // recall that when the source points to registered data, it is persistent, and
      // thus does not need to invoke the getCtiData call. Since we know that any
      // NULL source must be unregiastered, we use the "getUnregisteredCtiData" routine
      // simply to save the time associated with traversing the map of registered data,
      // because we know it will not be found in there.

      CtiData * source = statsVec[ii].source;
      if (source == NULL) {
	source = getUnregisteredCtiData(statsVec[ii].name);
	assert(source);
      }
      const int dtype = statsVec[ii].datatype;
      assert(source->getType() == dtype);
      const double tmp = 1.0 / ( statsVec[ii].wgt->d() + wgt );
      if (dtype == D_DATA) {
        statsVec[ii].rms->d() = ( statsVec[ii].wgt->d() * ( statsVec[ii].rms->d() * statsVec[ii].rms->d() + statsVec[ii].avg->d() * statsVec[ii].avg->d() ) +
				  wgt * ( source->d() * source->d() ) ) * tmp;
        statsVec[ii].avg->d() = ( statsVec[ii].wgt->d() * statsVec[ii].avg->d() + wgt * source->d() ) * tmp;
        statsVec[ii].rms->d() = sqrt( fabs( statsVec[ii].rms->d() - statsVec[ii].avg->d()*statsVec[ii].avg->d() ) );
        statsVec[ii].wgt->d() += wgt;
      }
      else if (dtype == DN_DATA) {
        assert(source->size() == statsVec[ii].avg->size());
        assert(source->size() == statsVec[ii].rms->size());
        double * rms = statsVec[ii].rms->getDNptr();
        double * avg = statsVec[ii].avg->getDNptr();
        double * src = source->getDNptr();
        const int size = source->size();
        const double wgt0 = statsVec[ii].wgt->d();
        for (int i = 0; i < size; ++i) {
          rms[i] = ( wgt0 * ( rms[i]*rms[i] + avg[i]*avg[i] ) + wgt*src[i]*src[i] ) * tmp;
          avg[i] = ( wgt0*avg[i] + wgt*src[i] ) * tmp;
        }
        for (int i = 0; i < size; ++i) {
          rms[i] = sqrt( fabs(rms[i]-avg[i]*avg[i]) );
        }
        statsVec[ii].wgt->d() += wgt;
      }
      else if (dtype == DN3_DATA) {
        assert(source->size() == statsVec[ii].avg->size());
        assert(source->size() == statsVec[ii].rms->size());
        assert(source->size() == statsVec[ii].rey->size());
        double (*rey)[3] = statsVec[ii].rey->getDN3ptr();
        double (*rms)[3] = statsVec[ii].rms->getDN3ptr();
        double (*avg)[3] = statsVec[ii].avg->getDN3ptr();
        double (*src)[3] = source->getDN3ptr();
        const int size = source->size();
        const double wgt0 = statsVec[ii].wgt->d();
        for (int i = 0; i < size; ++i) {
          FOR_J3 rms[i][j] = ( wgt0 * ( rms[i][j]*rms[i][j] + avg[i][j]*avg[i][j] ) + wgt*src[i][j]*src[i][j] ) * tmp;
          rey[i][0] = ( wgt0 * ( rey[i][0] + avg[i][1]*avg[i][2] ) + wgt * ( src[i][1]*src[i][2] ) ) * tmp;
          rey[i][1] = ( wgt0 * ( rey[i][1] + avg[i][2]*avg[i][0] ) + wgt * ( src[i][2]*src[i][0] ) ) * tmp;
          rey[i][2] = ( wgt0 * ( rey[i][2] + avg[i][0]*avg[i][1] ) + wgt * ( src[i][0]*src[i][1] ) ) * tmp;
          FOR_J3 avg[i][j] = ( wgt0*avg[i][j] + wgt*src[i][j] ) * tmp;
        }
        for (int i = 0; i < size; ++i) {
          FOR_J3 rms[i][j] = sqrt( fabs( rms[i][j] - avg[i][j]*avg[i][j] ) );
          rey[i][0] -= avg[i][1] * avg[i][2];
          rey[i][1] -= avg[i][2] * avg[i][0];
          rey[i][2] -= avg[i][0] * avg[i][1];
        }
        statsVec[ii].wgt->d() += wgt;
      }
      else {
	assert(0);
      }

      if ((verbose)&&(mpi_rank == 0)) {
	cout << " > updateStats() \"" << statsVec[ii].name << "\", total averaging time: " << statsVec[ii].wgt->d() << endl;
      }

    }

  }

  void clearCurrentData() {
    currentDataMap.clear();
  }

} // namespace CtiRegister
