#ifdef NO_MPI

#include <iostream>
#include <assert.h>
#include <string.h>
#include <cstdio>
#include "NoMpi.hpp"
#include <sys/time.h>
#include <unistd.h>

using namespace std;

namespace NoMpiBuffer { 
  map<int,char*> mpi_buffers;
  map<int,pair<char*,int> > recv_req;
};

int MPI_Comm_dup(MPI_Comm comm, MPI_Comm * mycomm) {
  *mycomm = comm;
  return (0);
}

int MPI_Comm_rank(MPI_Comm comm, int * mpi_rank) {
  *mpi_rank = 0;
  return (0);
}

int MPI_Comm_size(MPI_Comm mycomm, int * mpi_size) {
  *mpi_size = 1;
  return (0);
}

int MPI_Comm_group(MPI_Comm comm, MPI_Group * group) {
  *group = 0 ;
  return (0) ;
}

int MPI_File_open(MPI_Comm comm,const char *fname, int mode, int info, MPI_File *fp) {

  char cmode[10];
  switch (mode) {
    case MPI_MODE_RDONLY:
    sprintf(cmode, "rb");
    break;
    case MPI_MODE_WRONLY:
    sprintf(cmode, "wb");   // only binary
    break;
    case MPI_MODE_WRONLY | MPI_MODE_CREATE:
    sprintf(cmode, "wb");   // only binary
    break;
    default:
    std::cerr << "Error: cannot figure out mode" << std::endl;
    break;
  }

  if ((*fp=fopen(fname, cmode)) == NULL) {
    printf("could not open %s\n", fname);
    return -1;
  }

  return 0;
}

int MPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype,const char * datarep, int info) {
  fseek (fh ,disp ,SEEK_SET);
  return 0;
}

int MPI_File_close(MPI_File *fp) {
  fclose(*fp);
  return 0;
}

int MPI_File_delete(const char *fname,MPI_Info info) {
  // remove returns 0 if successfull...
  return remove(fname);
}

int MPI_Type_free(MPI_Datatype *datatype) {
  *datatype = MPI_DATATYPE_NULL;
  return 0;
}

int MPI_Type_indexed(int count, int *blocklens, int *indices, MPI_Datatype old_type, MPI_Datatype *newtype) {
  *newtype = 1; // something other than MPI_DATATYPE_NULL
  return 0;
}

int MPI_Type_hindexed(int count, int *blocklens, MPI_Aint *indices, MPI_Datatype old_type, MPI_Datatype *newtype) {
  *newtype = 1; // something other than MPI_DATATYPE_NULL
  return 0;
}

int MPI_Type_create_hindexed(int count, int *blocklens, MPI_Aint *indices, MPI_Datatype old_type, MPI_Datatype *newtype) {
  *newtype = 1; // something other than MPI_DATATYPE_NULL
  return 0;
}

int MPI_Type_commit(MPI_Datatype *datatype) {
  assert( *datatype != MPI_DATATYPE_NULL );
  return 0;
}

int MPI_Type_struct(int count, int *blocklens, MPI_Aint *indices, MPI_Datatype *old_types, MPI_Datatype *newtype) {
  // here we use the datatype to store the size...
  *newtype = 0;
  for (int i=0; i<count; i++)
    *newtype += blocklens[i]*old_types[i];
  return 0;
}

int MPI_Type_create_struct(int count, int *blocklens, MPI_Aint *indices, MPI_Datatype *old_types, MPI_Datatype *newtype) {
  return MPI_Type_struct(count,blocklens,indices,old_types,newtype);
}

int MPI_Type_get_extent(MPI_Datatype datatype, MPI_Aint *lb,MPI_Aint *extent) {
  *lb = 0; // lower bound
  *extent = (long int)datatype;
  return 0;
}

int MPI_Type_size(MPI_Datatype datatype, int *size) {
  *size = (long int)datatype;
  return 0;
}

double MPI_Wtime() {
  static bool first = true ;
  static bool disabled = false;
  static timeval start;
  if ( first ) {
    if (gettimeofday(&start,NULL) < 0){
      cout << "Warning: timing disabled, problem with gettimeofday()" << endl;
      disabled = true;
    }
    first = false ;
  }
  if (disabled){
    return 0.0;
  }

  //microsecond resolution
  double start_seconds = start.tv_sec + start.tv_usec*1e-6;
  timeval now;
  gettimeofday(&now,NULL);

  double tmp = now.tv_sec + now.tv_usec*1e-6 - (start_seconds) ;
  return tmp;
}

int MPI_File_get_size(MPI_File fh, MPI_Offset * size) { 
  fseek(fh,0L,SEEK_END);
  const int sz = ftell(fh);
  rewind(fh);
  return sz;
}

int MPI_File_set_size(MPI_File fh,MPI_Offset size) {
  return 0;
}

int MPI_File_seek(MPI_File fh, MPI_Offset offset,int whence) {
  fseek(fh,offset,SEEK_SET);
  return 0;
}

int MPI_Address(void *location, MPI_Aint *address) {
  *address = (MPI_Aint)location;
  return 0;
}

int MPI_Get_address(void *location, MPI_Aint *address) { return MPI_Address(location,address); }

int MPI_Get_processor_name(char *name,int *resultlen) {
  //unlink real mpi routine, serial gethostname will not update resultlen, 
  //so require the user to pass in the max length of name.
  //if gethostname returns successfull, decrease resultlen by one
  //as another indication to calling routine of success.
  assert(*resultlen==MPI_MAX_PROCESSOR_NAME);
  int ierr =  gethostname(name,*resultlen);
  if (ierr==0) //success
    *resultlen = MPI_MAX_PROCESSOR_NAME-1;
  return ierr;
  //sprintf(name,"one");
  //*resultlen = 3;
  //return 0;
}

int MPI_Error_string(int ierr,char * message,int *length) {
  sprintf(message,"NoMpi Error %d",ierr);
  *length = strlen(message);
  return 0;
}

int MPI_Testall(int count, MPI_Request array_of_requests[], int* flag, MPI_Status array_of_statuses[]) {
  *flag = 1;
  return 0;
}

int MPI_Get_count( const MPI_Status *status, MPI_Datatype datatype, int *count ) {
  *count = 0;
  return 0;
}

int MPI_Init(int * argc, char *** argv) {
  return (0);
}

int MPI_Comm_free(MPI_Comm * mycomm) {
  return 0;
}

int MPI_Comm_split(MPI_Comm comm,int color,int key,MPI_Comm * comm_out) {
  return MPI_Comm_dup(comm,comm_out);
}

int MPI_Barrier(MPI_Comm comm) {
  return (0);
}

int MPI_Finalize() { return (0); }

int MPI_Abort(MPI_Comm comm,int errorcode) {
  return(errorcode);
}

int MPI_Wait(MPI_Request *request, MPI_Status * status) { 
  return(0);
}


int MPI_Waitall(int count, MPI_Request * requests, MPI_Status * statuses) {

  using namespace NoMpiBuffer;

  // if there are any recv_request pending, please complete them.. 

  for (map<int,pair<char*,int> >::iterator it = recv_req.begin(); it != recv_req.end(); ++it) { 

    map<int,char*>::iterator send_it = mpi_buffers.find(it->first);
    if ( send_it != mpi_buffers.end()) { 
      memcpy(it->second.first,send_it->second,sizeof(char)*it->second.second);
      delete[] send_it->second;
      mpi_buffers.erase(send_it);
      it->second.first = NULL; // to note completion.. 

    } else { 
      assert(0); // unable to find an associated send.  not good.
    }
  }

  // ensure that everyone is completed 

  assert( mpi_buffers.size() == 0);
  recv_req.clear();

  return (0);
}

int MPI_Waitany(int count, MPI_Request * requests, int *index, MPI_Status * statuses) {

  // this functionality in nompi mode is presently broken, because 
  // we do not have a facility where we are keeping track of requests...
  assert(0); 
  return (0);
}

MPI_Fint MPI_Comm_c2f(MPI_Comm comm) { return((MPI_Fint)comm); }

int MPI_Request_free(MPI_Request* req) {
  return 0;
}

int MPI_Win_fence(int assert,MPI_Win& win) {
  return 0;
}

int MPI_Win_free(MPI_Win* win) {
  return 0;
}
#endif
