#ifndef COMMONIO_HPP
#define COMMONIO_HPP

#include "Common.hpp"
#include "DistributedDataExchanger.hpp"
#include "ByteSwap.hpp"
#include "Utils.hpp"

namespace IoParams {
  static const int8 n_chnk = (1ll<<26);
};

inline void buildDistAndDde(int8*& xora_node, int8*& xora_full, int& nx_s,
                            DistributedDataExchanger* &dde,
                            const int8* ix_global, const int n,
                            const int8 nx_global, const int nthresh = -1) {
  assert(xora_node == NULL);
  assert(xora_full == NULL);
  assert(dde       == NULL);

  if ( mpi_rank_shared == 0 ) {
    if ( nthresh > 0 )
      MiscUtils::calcThresholdDist(xora_node,nx_global,mpi_size_internode,nthresh);
    else
      MiscUtils::calcUniformDist(xora_node,nx_global,mpi_size_internode);

    nx_s = xora_node[mpi_rank_internode+1] - xora_node[mpi_rank_internode];
    assert( nx_s < TWO_BILLION);
  } else {
    nx_s = 0;
  }

  MiscUtils::buildXora(xora_full,nx_s);
  dde = new DistributedDataExchanger(ix_global,n,xora_full);
}

inline void writeHeaderRank0(MPI_File &fh, MPI_Offset& offset, const string& name,
                             const int8 skip, const int id, const int8 n_global,const int K= -1) {
  if ( mpi_rank == 0 ) {
    Header header;
    sprintf(header.name,"%s",name.c_str());
    header.id    = id;
    header.skip  = skip;
    ByteSwap::setLswMswPairForInt8(header.idata+0,n_global);
    if ( K >= 0) header.idata[2] = K;
    MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
  }
}

template<typename T>
void readChunkedData(MPI_File& fh, const MPI_Offset offset, T* my_buf, int8 my_count,
                     const bool byte_swap, MPI_Comm& mpi_comm) {

  MPI_Datatype MPI_T     = MpiStuff::getMpiDatatype<T>();
  const int8 T_size      = int8(sizeof(T));
  MPI_Offset this_offset = offset;
  int8 n_remain          = my_count;

  // during the chunked data reading, if everyone has something to read, then
  // we will issue collective read calls.  otherwise they will issue non-collective
  // calls.

  int done = 0;
  while ( done == 0 ) {

    int my_done = 1;

    // check to see if everyone is going to be reading ...

    const int8 this_count = min(n_remain,IoParams::n_chnk); assert(this_count >= 0);
    const int8 n_read     = my_count-n_remain;

    if ( this_count > 0)
      my_done = 0;

    MPI_File_read_at_all(fh,this_offset,my_buf+n_read,this_count,MPI_T,MPI_STATUS_IGNORE);

    n_remain    -= this_count; assert(n_remain >= 0);
    this_offset += this_count*T_size;

    MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);

  }

  assert(this_offset == (offset+T_size*my_count));

  if ( byte_swap) ByteSwap::byteSwap(my_buf,my_count);

}

template<typename T>
void writeChunkedData(MPI_File& fh, MPI_Offset offset, T* my_buf, int8 my_count,
                      MPI_Comm& mpi_comm) {

  MPI_Datatype MPI_T      = MpiStuff::getMpiDatatype<T>();
  const int8 T_size       = int8(sizeof(T));

  MPI_Offset this_offset  = offset;
  int8 n_remain           = my_count;

  while ( true ) {

    int8 this_count       = min(n_remain,IoParams::n_chnk); assert(this_count >= 0);
    const int8 n_read     = my_count-n_remain;

    int8 my_counts[2] = {this_count,-this_count};
    int8 counts[2];
    MPI_Allreduce(my_counts,counts,2,MPI_INT8,MPI_MAX,mpi_comm);

    counts[1] = -counts[1];

    // everyone is done because the max of the counts is 0...
    if ( counts[0] == 0)
      break;

#ifdef SERIAL_IO
    MPI_File_write_at_all(fh,this_offset,my_buf+n_read,this_count,MPI_T,MPI_STATUS_IGNORE);
#else
    if (  counts[1] > 0) {
      // collective operation because everyone has something to write...
      MPI_File_write_at_all(fh,this_offset,my_buf+n_read,this_count,MPI_T,MPI_STATUS_IGNORE);
    } else if ( this_count > 0) {
      MPI_File_write_at(fh,this_offset,my_buf+n_read,this_count,MPI_T,MPI_STATUS_IGNORE);
    }
#endif

    n_remain             -= this_count; assert( n_remain >= 0);
    this_offset          += this_count*T_size;
  }

  assert(this_offset == (offset+T_size*my_count));

}

template<typename T>
void writeChunkedDataAscii(MPI_File& fh,MPI_Offset& offset, T* my_buf, int8 my_count,MPI_Comm& mpi_comm,const string fmt) {
  // due to unknown file offsets based on write length of formatted chars, we
  // use a 2-pass system to write the chunked ASCII data. The first pass counts
  // the number of chars to write, while the second does the parallel writing

  // first pass to simply count the chars to be written by each rank
  MPI_Offset my_offset  = 0;
  int8 n_remain         = my_count;

  while ( true ) {
    int cbuf_n = 0;

    const int8 this_count = min(n_remain,IoParams::n_chnk); assert(this_count >= 0);
    const int8 n_read     = my_count-n_remain;

    int8 my_counts[2] = {this_count,-this_count};
    int8 counts[2];
    MPI_Allreduce(my_counts,counts,2,MPI_INT8,MPI_MAX,mpi_comm);

    counts[1] = -counts[1];

    // everyone is done because the max of the counts is 0...
    if ( counts[0] == 0) break;

    MPI_Offset cbuf_size = 40*this_count;  // accounts for space-delimiters and end-of-line char; related to max digits possible with fmt
    char * cbuf = new char[cbuf_size];

    for (int ii=0,end=int(this_count); ii<end; ++ii) {
      cbuf_n += sprintf(cbuf+cbuf_n,fmt.c_str(),my_buf[n_read+ii]);
    }

    n_remain             -= this_count; assert( n_remain >= 0);
    my_offset          += cbuf_n;

    DELETE(cbuf);
  }

  // my_offset now contains the number of chars I will write, so socialize
  int8 my_offset8 = my_offset;  // rank-local number of chars to write
  int8 final_offset;
  MPI_Allreduce(&my_offset8,&final_offset,1,MPI_INT8,MPI_SUM,mpi_comm);

  // on systems where every rank can write use offsets to do parallel I/O

  int8 my_disp;
  MPI_Scan(&my_offset8,&my_disp,1,MPI_INT8,MPI_SUM,mpi_comm);

  my_offset = offset + MPI_Offset(my_disp-int8(my_offset));

  // second pass to actually write
  n_remain           = my_count;
  while ( true ) {
    int cbuf_n = 0;

    int8 this_count       = min(n_remain,IoParams::n_chnk); assert(this_count >= 0);
    const int8 n_read     = my_count-n_remain;

    int8 my_counts[2] = {this_count,-this_count};
    int8 counts[2];
    MPI_Allreduce(my_counts,counts,2,MPI_INT8,MPI_MAX,mpi_comm);

    counts[1] = -counts[1];

    // everyone is done because the max of the counts is 0...
    if ( counts[0] == 0) break;

    MPI_Offset cbuf_size = 40*this_count;  // accounts for space-delimiters and end-of-line char; related to max digits possible with fmt
    char * cbuf = new char[cbuf_size];

    for (int ii=0,end=int(this_count); ii<end; ++ii) {
      cbuf_n += sprintf(cbuf+cbuf_n,fmt.c_str(),my_buf[n_read+ii]);
    }

#ifdef SERIAL_IO
    MPI_File_write_at_all(fh,my_offset,cbuf,cbuf_n,MPI_CHAR,MPI_STATUS_IGNORE);
#else
    if (  counts[1] > 0) {
      // collective operation because everyone has something to write...
      MPI_File_write_at_all(fh,my_offset,cbuf,cbuf_n,MPI_CHAR,MPI_STATUS_IGNORE);
    } else if ( this_count > 0) {
      MPI_File_write_at(fh,my_offset,cbuf,cbuf_n,MPI_CHAR,MPI_STATUS_IGNORE);
    }
#endif

    n_remain             -= this_count; assert( n_remain >= 0);
    my_offset          += cbuf_n;

    DELETE(cbuf);
  }

  // make sure everybody has final offset after all ranks write
  offset += MPI_Offset(final_offset);
}

template<typename T>
void readXR1(MPI_File& fh, MPI_Offset& offset, T* buf, const int8* xora,
             const int n,const bool byte_swap, MPI_Comm& mpi_comm) {

  assert( xora != NULL);
  const int8 T_size  = int8(sizeof(T));
  MPI_Offset my_offset = offset+header_size+xora[mpi_rank_internode]*T_size;
  readChunkedData(fh,my_offset,buf,n,IoParams::n_chnk,byte_swap,mpi_comm);
}

template<typename T, int K>
void readXR2Stride(MPI_File& fh, MPI_Offset& offset, T (*buf)[K], const int8* xora,
             const int n,const bool byte_swap, MPI_Comm& mpi_comm) {

  assert( xora != NULL);
  const int8 T_size  = int8(sizeof(T));
  MPI_Offset my_offset = offset+header_size+xora[mpi_rank_internode]*int8(K)*T_size;
  readChunkedData(fh,my_offset,(T*)buf,K*n,IoParams::n_chnk,byte_swap,mpi_comm);
}

template<typename T>
void readXR2(MPI_File& fh, MPI_Offset& offset, T (*buf)[3], const int8* xora,
             const int n,const bool byte_swap, MPI_Comm& mpi_comm) {
  readXR2Stride<T,3>(fh,offset,buf,xora,n,byte_swap,mpi_comm);
}

template<typename T>
void io_xcheck(MPI_File& fh, MPI_Offset& offset, const int8* xora,const int n,
               const char* msg, const bool byte_swap,MPI_Comm& mpi_comm) {

  T* x_check         = new T[n];
  readXR1<T>(fh,offset,x_check,xora,n,byte_swap,mpi_comm);

  for (int i =0; i < n ; ++i)
    assert( x_check[i] == xora[mpi_rank_internode]+i);
  delete[] x_check;

  if ( mpi_rank == 0 )
    cout << msg << endl;
}

template<typename T>
void writeXR1(MPI_File& fh, MPI_Offset& offset,const T* buf, const int id,
              const string& name, DistributedDataExchanger* dde,
              const int8* daora_node, const int8 n_global,MPI_Comm& mpi_comm) {

  T* buf_s          = new T[dde->nda]; // striped data.
  const int8 T_size = int8(sizeof(T));
  const int8 skip   = header_size + T_size*n_global;
  dde->push(buf_s,buf);

  if ( mpi_rank_shared == 0 ) {

    writeHeaderRank0(fh,offset,name,skip,id,n_global);

    const int n_s = int(daora_node[mpi_rank_internode+1]-daora_node[mpi_rank_internode]);
    assert( n_s == dde->nda);

    MPI_Offset my_offset = offset+header_size+daora_node[mpi_rank_internode]*T_size;
    writeChunkedData(fh,my_offset,buf_s,n_s,IoParams::n_chnk,mpi_comm);
  }

  offset += skip;
  delete[] buf_s;
}

template<typename T,int K>
void writeXR2Stride(MPI_File& fh, MPI_Offset& offset, const T (*buf)[K], const int id,
                    const string& name, DistributedDataExchanger* dde,
                    const int8* daora_node, const int8 n_global,MPI_Comm& mpi_comm) {

  T (*buf_s)[K]     = new T[dde->nda][K];
  const int8 T_size = int8(sizeof(T));
  const int8 skip   = header_size + int8(K)*T_size*n_global;
  dde->push(buf_s,buf);

  if ( mpi_rank_shared == 0) {

    writeHeaderRank0(fh,offset,name,skip,id,n_global,K);

    const int n_s = int(daora_node[mpi_rank_internode+1]-daora_node[mpi_rank_internode]);
    assert( n_s == dde->nda);

    MPI_Offset my_offset = offset+header_size+daora_node[mpi_rank_internode]*int8(K)*T_size;
    writeChunkedData(fh,my_offset,(T*)buf_s,K*n_s,IoParams::n_chnk,mpi_comm);
  }

  offset += skip;
  delete[] buf_s;
}

template<typename T>
void writeXR2(MPI_File& fh, MPI_Offset& offset, const T (*buf)[3], const int id,
              const string& name, DistributedDataExchanger* dde,
              const int8* daora_node, const int8 n_global,MPI_Comm& mpi_comm) {
  writeXR2Stride<T,3>(fh,offset,buf,id,name,dde,daora_node,n_global,mpi_comm);
}

#endif
