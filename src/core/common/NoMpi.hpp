#ifndef NOMPI_HPP
#define NOMPI_HPP

#include <iostream>
#include <assert.h>
#include <cstring>
#include <map>
#include <stdio.h>

using namespace std;

/// dummy MPI defines
/// defines are size of types !?!?!?!?!
#define MPI_CHAR               1
#define MPI_UNSIGNED_CHAR      1
#define MPI_SHORT              2
#define MPI_INT                4
#define MPI_FLOAT              4
#define MPI_BYTE               1
#define MPI_INTEGER            4
#define MPI_DOUBLE             8
#define MPI_DOUBLE_PRECISION   8
#define MPI_LONG               8
#define MPI_LONG_LONG          8
#define MPI_UNSIGNED           4
#define MPI_UNSIGNED_LONG      8
#define MPI_UNSIGNED_LONG_LONG 8
#define MPI_UNSIGNED_SHORT     2
#define MPI_DOUBLE_INT         12
#define MPI_DOUBLE_COMPLEX     16
#define MPI_UNDEFINED          -1
#define MPI_OFFSET             8

#define MPI_MAX 1
#define MPI_SUM 2
#define MPI_MIN 3
#define MPI_MINLOC 4
#define MPI_LOR 5

#define MPI_COMM_WORLD 1
#define MPI_COMM_SELF  1
#define MPI_INFO_NULL 1
#define MPI_COMM_NULL -1

#define MPI_MODE_RDONLY              2  // ADIO_RDONLY
#define MPI_MODE_RDWR                8  // ADIO_RDWR
#define MPI_MODE_WRONLY              4  // ADIO_WRONLY
#define MPI_MODE_CREATE              1  // ADIO_CREATE
#define MPI_MODE_EXCL               64  // ADIO_EXCL
#define MPI_MODE_DELETE_ON_CLOSE    16  // ADIO_DELETE_ON_CLOSE
#define MPI_MODE_UNIQUE_OPEN        32  // ADIO_UNIQUE_OPEN
#define MPI_MODE_APPEND            128  // ADIO_APPEND
#define MPI_MODE_SEQUENTIAL        256  // ADIO_SEQUENTIAL
#define MPI_SEEK_SET  1
#define MPI_SEEK_CUR  1

#define MPI_MAX_ERROR_STRING 64
#define MPI_MAX_PROCESSOR_NAME 256

// start of addr space, unsafe
#define MPI_BOTTOM (char*)0

#define MPI_STATUS_IGNORE 0
#define MPI_REQUEST_NULL 0
#define MPI_STATUSES_IGNORE 0

#define MPI_IN_PLACE 0

/// typedefs
typedef int MPI_Comm;
typedef int MPI_Status;
typedef int MPI_Request;
typedef int MPI_Datatype;
typedef int MPI_Info;
typedef int MPI_Op;
typedef int MPI_Fint;
typedef int MPI_Group;

// mpi_win class is used to spoof
// the one-sided communication RMA
// protocol
struct MPI_Win {
  void * base ;
  int disp_unit ;

} ;

#define MPI_DATATYPE_NULL ((MPI_Datatype)0)

#ifdef INT8_IS_LONG_INT
typedef long int MPI_Offset;
typedef long int MPI_Aint;
#else
typedef long long int MPI_Offset;
typedef long long int MPI_Aint;
#endif

#define MPI_AINT 8

typedef FILE * MPI_File;

namespace NoMpiBuffer { 
  extern map<int,char*> mpi_buffers;
  extern map<int,pair<char*,int> > recv_req;
};

int MPI_Init(int * argc, char *** argv);
int MPI_Comm_dup(MPI_Comm comm, MPI_Comm * mycomm);
int MPI_Comm_free(MPI_Comm * mycomm) ;
int MPI_Comm_split(MPI_Comm comm,int color,int key,MPI_Comm * comm_out);
int MPI_Comm_rank(MPI_Comm comm, int * mpi_rank);
int MPI_Comm_size(MPI_Comm mycomm, int * mpi_size);
int MPI_Comm_group( MPI_Comm comm, MPI_Group* group) ;
MPI_Fint MPI_Comm_c2f(MPI_Comm comm);
int MPI_Barrier(MPI_Comm comm);
int MPI_Finalize();
int MPI_Abort(MPI_Comm comm,int errorcode);
int MPI_Waitall(int count, MPI_Request *requests, MPI_Status *status);
int MPI_Waitany(int count, MPI_Request *requests, int *index, MPI_Status *status);
int MPI_Wait(MPI_Request *request, MPI_Status * status); 

int MPI_File_open(MPI_Comm comm,const char *fname, int mode, int info, MPI_File *fp);
int MPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype,
		      const char * datarep, int info);
int MPI_File_close(MPI_File *fh);
int MPI_File_delete(const char *fname,MPI_Info info);

int MPI_Type_free(MPI_Datatype *datatype);
int MPI_Type_indexed(int count, int *blocklens, int *indices, MPI_Datatype old_type,
		     MPI_Datatype *newtype);
int MPI_Type_hindexed(int count, int *blocklens, MPI_Aint *indices, MPI_Datatype old_type,
		      MPI_Datatype *newtype);
int MPI_Type_create_hindexed(int count, int *blocklens, MPI_Aint *indices, MPI_Datatype old_type,
			     MPI_Datatype *newtype);
int MPI_Type_commit(MPI_Datatype *datatype);
int MPI_Type_struct(int count, int *blocklens, MPI_Aint *indices, MPI_Datatype *old_types,
		    MPI_Datatype *newtype);
int MPI_Type_create_struct(int count, int *blocklens, MPI_Aint *indices, MPI_Datatype *old_types,
			   MPI_Datatype *newtype);
int MPI_Type_get_extent(MPI_Datatype datatype, MPI_Aint *lb,MPI_Aint *extent);
int MPI_Type_size(MPI_Datatype datatype, int *size);
int MPI_File_seek(MPI_File fh, MPI_Offset offset,int whence);
int MPI_Address(void *location, MPI_Aint *address);
int MPI_Get_address(void *location, MPI_Aint *address);
int MPI_Get_processor_name(char *name,int *resultlen);

int MPI_Error_string(int ierr,char * message,int *length);

double MPI_Wtime();

int MPI_File_set_size(MPI_File fh,MPI_Offset size);
int MPI_File_get_size(MPI_File fh, MPI_Offset * size); 
int MPI_Testall(int count, MPI_Request array_of_requests[], int* flag, MPI_Status array_of_statuses[]);
int MPI_Get_count( const MPI_Status *status, MPI_Datatype datatype, int *count );
int MPI_Request_free(MPI_Request* req) ;

// one-sided communication
template<class T>
int MPI_Win_create( T* base, int size, int disp_unit, MPI_Info info, MPI_Comm comm, MPI_Win * win) {

  // no checking of the end of the window copy
  // using size argument

  win->base       = (void*) base ;
  win->disp_unit  = disp_unit ;
  return 0 ;
}

template<class T>
int MPI_Put( T* origin_addr, int count, MPI_Datatype o_type, int target_rank,
             MPI_Aint target_disp, int target_count, MPI_Datatype t_type, MPI_Win& win) {

  char * target_start = (char*)win.base +  target_disp * win.disp_unit ;
  int nbytes          = count * o_type ;
  // no checking to ensure origin and target are of same size
  memcpy ( (void*) target_start, (void*) origin_addr, nbytes ) ;

  return 0;
}


template<class T>
int MPI_Get( T* origin_addr, int count, MPI_Datatype o_type, int target_rank,
             MPI_Aint target_disp, int target_count, MPI_Datatype t_type, MPI_Win& win) {

  // note that the api is declared to be consistent with the MPI docs, but
  // data is copied from target TO origin in MPI_Get
  char * target_start = (char*)win.base + target_disp * win.disp_unit ;
  int nbytes          = count * o_type ;
  // no checking to ensure origin and target are of same size
  memcpy ( (void*) origin_addr, (void *) target_start, nbytes) ;
  
  return(0) ;
}

int MPI_Win_fence(int assert, MPI_Win& win);
int MPI_Win_free( MPI_Win * win) ;


//
// implementation of DUMMY-MPI functions when compiling without mpi,
// function templates have to be implemented in the header file
//

///
/// MPI imitator fuction is using C style for io stream
/// instead of using sizeof(T) for size specification in fread -
/// had to use DEFINE of datatype, see functions 'read_no_vector' and 'read_cvofa' ...
template<class T>
int MPI_File_read_all(MPI_File fp, T *buf, int count, MPI_Datatype datatype,MPI_Status *status) {
  fread(buf, datatype, count, fp);
  return 0;
}

template<class T>
int MPI_File_read(MPI_File fp, T *buf, int count, MPI_Datatype datatype,MPI_Status *status) {
  fread(buf, datatype, count, fp);
  return 0;
}

template<class T>
int MPI_File_write_all(MPI_File fp, T *buf, int count, MPI_Datatype datatype,MPI_Status *status) {
  fwrite(buf, datatype, count, fp);
  return 0;
}

template<class T>
int MPI_File_write(MPI_File fp, T *buf, int count, MPI_Datatype datatype,MPI_Status *status) {
  fwrite(buf, datatype, count, fp);
  return 0;
}

template<class T>
int MPI_File_write_at(MPI_File fp, MPI_Offset offset, T *buf, int count, MPI_Datatype datatype,MPI_Status *status) {
  MPI_File_seek(fp,offset,0);
  fwrite(buf, datatype, count, fp);
  return 0;
}

template<class T>
int MPI_File_write_at_all(MPI_File fp, MPI_Offset offset, T *buf, int count, MPI_Datatype datatype,MPI_Status *status) {
  MPI_File_seek(fp,offset,0);
  fwrite(buf, datatype, count, fp);
  return 0;
}

template<class T>
int MPI_File_read_at(MPI_File fp, MPI_Offset offset, T *buf, int count, MPI_Datatype datatype,MPI_Status *status) {
  fseek(fp,offset,SEEK_SET);
  fread(buf,datatype,count,fp) ;
  return 0;
}

template<class T>
int MPI_File_read_at_all(MPI_File fp, MPI_Offset offset, T *buf, int count, MPI_Datatype datatype,MPI_Status *status) {
  fseek(fp,offset,SEEK_SET);
  fread(buf,datatype,count,fp) ;
  return 0;
}

template<class T>
int MPI_Reduce(T * local, T * global, int n, int type, int action, int rank, MPI_Comm comm) {
  for (int i = 0; i < n; i++)
    global[i] = local[i];
  return 0;
}

template<class T>
int MPI_Reduce(const int mpi_in_place, T * global, int n, int type, int action, int rank, MPI_Comm comm) {
  assert(mpi_in_place == MPI_IN_PLACE);
  return 0;
}

template<class T>
int MPI_Reduce(T (*local)[2], T (*global)[2], int n, int type, int action, int rank, MPI_Comm comm) {
  
  // &(global[0]) cannot be used here.  n could be zero.
  // this applies for all of the following... 
  
  assert(n%2 == 0);
  for (int i = 0; i < n/2; i++)
    for (int j = 0; j < 2; j++)
      global[i][j] = local[i][j];
  return 0;
}
template<class T>
int MPI_Reduce(T (*local)[3], T (*global)[3], int n, int type, int action, int rank, MPI_Comm comm) {
  assert(n%3 == 0);
  for (int i = 0; i < n/3; i++)
    for (int j = 0; j < 3; j++)
      global[i][j] = local[i][j];
  return 0;
}
template<class T>
int MPI_Reduce(T (*local)[8], T (*global)[8], int n, int type, int action, int rank, MPI_Comm comm) {
  assert(n%8 == 0);
  for (int i = 0; i < n/8; i++)
    for (int j = 0; j < 8; j++)
      global[i][j] = local[i][j];
  return 0;
}

template<class T>
int MPI_Allreduce(T *local, T *global, int n, int type, int action, MPI_Comm comm) {
  for (int i = 0; i < n; i++)
    global[i] = local[i];
  return 0;
}

template<class T>
int MPI_Allreduce(int in_place, T *global, int n, int type, int action, MPI_Comm comm) {
  return 0;
}

template<class T>
int MPI_Send(T *data, int count, int type, int rank, int tag, MPI_Comm comm) {
  
  // in NoMpi mode, the only way this will work is if an Irecv has 
  // already been posted... 
 
  using namespace NoMpiBuffer;

  int ierr = -1;

  for (map<int,pair<char*,int> >::iterator it = recv_req.begin(); 
       it != recv_req.end(); ++it) { 

    // looking to match with the tag we supplied .. 

    if ( it->first == tag) {
      assert( ierr == -1); // this will check no duplicated tags.. 
      memcpy(it->second.first,(void*)data,sizeof(char)*it->second.second);
      recv_req.erase(it++); // this request has been serviced... 
      ierr = 0;
    }
  }

  assert( ierr == 0);
  return 0;
}

template<class T>
int MPI_Ssend(T *data, int count, int type, int rank, int tag, MPI_Comm comm) {
  return MPI_Send(data,count,type,rank,tag,comm);
}

template<class T>
int MPI_Recv(T *data, int count, MPI_Datatype datatype, int rank, int tag,
	     MPI_Comm comm, MPI_Status * status) {
 
  using namespace NoMpiBuffer;

  // in NoMpi mode, the only way this will work is if an Isend/Issend has 
  // already been posted... 
  
  int ierr = -1;

  for (map<int,char*>::iterator it = mpi_buffers.begin(); it != mpi_buffers.end(); ++it) { 

    if ( it->first == tag) { 
      
      assert( ierr == -1);

      // matching tag in the posted isends ... 

      memcpy(data,it->second,sizeof(T)*count);
      delete[] it->second;
      mpi_buffers.erase(it++);
      ierr = 0;
    }
  }

  assert( ierr = 0);
  return 0;
}

template<class T>
int MPI_Sendrecv(T * sbuf, int ns, int type_s, int rank_s, int tag_s, T * rbuf,
		 int nr, int type_r, int rank_r, int tag_r, MPI_Comm
		 comm, MPI_Status * status) {
  
  
  assert( type_s == type_r);
  assert( ns == nr);
  assert( rank_s == rank_r);

  for (int i = 0; i < ns; ++i) 
    rbuf[i] = sbuf[i];
  
  return 0;
}

template<class T>
int MPI_Alltoall(T *sendbuf, int sendcount, MPI_Datatype sendtype, T *recvbuf,int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
  for (int i = 0; i < sendcount; ++i) 
    recvbuf[i] = sendbuf[i];
  return 0;
}

template<class T>
int MPI_Alltoallv(T *sendbuf, int *sendcnts, int *sdispls, MPI_Datatype sendtype,
    T *recvbuf, int *recvcnts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm) {

  for (int i = 0; i < sendcnts[0]; i++)
    recvbuf[i] = sendbuf[i];
  return 0;
}

template<class T>
int MPI_Alltoallv(T (*sendbuf)[3], int *sendcnts, int *sdispls, MPI_Datatype sendtype,
		  T (*recvbuf)[3], int *recvcnts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm) {
  assert( sendcnts[0]%3 == 0 );
  for (int i = 0; i < sendcnts[0]/3; ++i) {
    recvbuf[i][0] = sendbuf[i][0];
    recvbuf[i][1] = sendbuf[i][1];
    recvbuf[i][2] = sendbuf[i][2];
  }
  return 0;
}

template<class T>
int MPI_Alltoallv(T *sendbuf, int *sendcnts, int *sdispls, MPI_Datatype sendtype,
		  T (*recvbuf)[3], int *recvcnts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm) {
  assert( sendcnts[0]%3 == 0 );
  for (int i = 0; i < sendcnts[0]/3; ++i) {
    recvbuf[i][0] = sendbuf[3*i];
    recvbuf[i][1] = sendbuf[3*i+1];
    recvbuf[i][2] = sendbuf[3*i+2];
  }
  return 0;
}

template<class T>
int MPI_Gather(T *sendbuf, int sendcnt, MPI_Datatype sendtype,
               T *recvbuf, int recvcount, MPI_Datatype recvtype,
               int root, MPI_Comm comm) {
  assert(sendcnt == recvcount);
  assert(root == 0);
  for (int i = 0; i < sendcnt; ++i) 
    recvbuf[i] = sendbuf[i];
  return 0;
}

template<class T>
int MPI_Gather(T (*buf1)[3], int n1, int type1, T (* buf2)[3], int n2, int type2,
	       int root, MPI_Comm comm) {
  assert( n1%3 == 0 );
  assert( n1 == n2);
  for (int i = 0; i < n1/3; ++i) {
    buf2[i][0] = buf1[0];
    buf2[i][1] = buf1[1];
    buf2[i][2] = buf1[2];
  }
  return 0;
}

template<class T>
int MPI_Gatherv(T *sendbuf, int sendcnt, MPI_Datatype sendtype,
                T *recvbuf, int *recvcnts, int *displs,
                MPI_Datatype recvtype,
                int root, MPI_Comm comm) {
  assert(sendcnt == recvcnts[0]);
  assert(root == 0);
  for (int i = 0; i < sendcnt; ++i) 
    recvbuf[i] = sendbuf[i];
  return 0;
}

template<class T>
int MPI_Allgatherv(T *sendbuf, int sendcnt, MPI_Datatype sendtype,
		   T *recvbuf, int *recvcnts, int *displs,
		   MPI_Datatype recvtype,
		   MPI_Comm comm) {
  assert(sendcnt == recvcnts[0]);
  for (int i = 0; i < sendcnt; ++i) 
    recvbuf[i] = sendbuf[i];
  return 0;
}

template<class T>
int MPI_Allgather(T * buf1, int n1, int type1, T * buf2, int n2, int type2,MPI_Comm comm) {
  for (int i = 0; i < n1; ++i) 
    buf2[i] = buf1[i];
  return 0;
}

template<class T>
int MPI_Allgather(T (buf1)[3], int n1, int type1, T (* buf2)[3], int n2, int type2,MPI_Comm comm) {
  assert( n1%3 == 0 );
  for (int i = 0; i < n1/3; ++i) {
    buf2[i][0] = buf1[0];
    buf2[i][1] = buf1[1];
    buf2[i][2] = buf1[2];
  }
  return 0;
}

template<class T>
int MPI_Scatter(T *sendbuf, int sendcnt, MPI_Datatype sendtype,
                T* recvbuf, int recvcnt, MPI_Datatype recvtype,
                int root, MPI_Comm comm) {
  for (int i = 0; i < recvcnt; ++i) 
    recvbuf[i] = sendbuf[i];
  return 0;
}


template<class T>
int MPI_Scatterv(T *sendbuf, int *sendcnts, int *displs, MPI_Datatype sendtype,
                 T *recvbuf, int recvcnt, MPI_Datatype recvtype,
                 int root, MPI_Comm comm) {
  assert(sendcnts[0] == recvcnt);
  assert(root == 0);
  for (int i = 0; i < recvcnt; ++i) 
    recvbuf[i] = sendbuf[i];
  return 0;
}

template<class T>
int MPI_Scan(T * buf1, T * buf2, int n, int type, int action, MPI_Comm comm) {
  for (int i = 0; i < n; ++i) 
    buf2[i] = buf1[i];
  return 0;
}

template<class T>
int MPI_Exscan(T * buf1, T * buf2, int n, int type, int action, MPI_Comm comm) {
  for (int i = 0; i < n; ++i) 
    buf2[i] = 0;
  return 0;
}

template<class T>
int MPI_Bcast(T * data, int n, int type, int rank, MPI_Comm comm) {
  return 0;
}

template<class T>
int MPI_Issend(T *data, int count, MPI_Datatype datatype, int rank, int tag,MPI_Comm comm, MPI_Request * request) {
  
  // copy the send buffer into a memory address .. 
  // recall MPI_Datatype is keeping a size of the data members

  char * buf = new char[count*int(datatype)];
  memcpy(buf,(char*)data,sizeof(char)*count*int(datatype));
  assert( NoMpiBuffer::mpi_buffers.find(tag) == NoMpiBuffer::mpi_buffers.end());
  NoMpiBuffer::mpi_buffers[tag] = buf;
  
  return 0;
}

template<class T>
int MPI_Irecv(T *data, int count, MPI_Datatype datatype, int rank, int tag,MPI_Comm comm, MPI_Request * request) {
  
  assert( NoMpiBuffer::recv_req.find(tag) == NoMpiBuffer::recv_req.end());
  NoMpiBuffer::recv_req[tag] = pair<char*,int>((char*)data,count*int(datatype));

  return 0;
}

template<class T>
int MPI_Isend(T *data, int count, MPI_Datatype datatype, int rank, int tag,MPI_Comm comm, MPI_Request * request ) {
 
  char * buf = new char[count*int(datatype)];
  memcpy(buf,(char*)data,sizeof(char)*count*int(datatype));
  assert( NoMpiBuffer::mpi_buffers.find(tag) == NoMpiBuffer::mpi_buffers.end());
  NoMpiBuffer::mpi_buffers[tag] = buf;
  
  return 0;
}

#endif
