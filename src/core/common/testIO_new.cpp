
#include <mpi.h>
#include <cstdio>
#include <iostream>
#include <assert.h>
#include <map>

using namespace std;

int mpi_rank;
int mpi_size;
MPI_Comm mpi_comm;

#define checkMpiError(error) \
  if ( error != MPI_SUCCESS) { \
    cout << " mpi error returned : " << error << "   " \
    << "   " << __FILE__ << ":" << __LINE__ << endl; \
    throw(0);\
  }


void calcUniformDist(int * &xod, const int nx, const int ndist) {
  if (xod == NULL) xod = new int[ndist+1];
  xod[0] = 0;
  for (int id = 1; id <= ndist; id++) {
    xod[id] = (int)((double)id/(double)ndist*(double)nx + 0.5);
    assert(xod[id] >= xod[id-1]);
  }
  assert(xod[ndist] == nx);
}


#ifdef SERIAL_IO

// in the mpich implementation, MPI_File is an int...  

int ifile = 0;
const int max_file_desc = 32;
int file_descp[max_file_desc];
map<MPI_File,FILE*> file_map;

int MPI_File_open(MPI_Comm comm, const char *filename, int amode,
    MPI_Info info, MPI_File * fh) { 

  if ( mpi_rank == 0) 
    cout << " intercepted an mpi file open call " << endl;

  // create the file handle for everyone ; this is a 
  // collective operation... 

  
  (*fh) = (MPI_File)(&file_descp[ifile]);
  ++ifile;

  FILE * fp = NULL;

  if ( mpi_rank == 0) { 

    if ( (amode & MPI_MODE_WRONLY)) { 

      fp = fopen(filename, "wb");

    } else if ( (amode & MPI_MODE_RDONLY)) { 

      fp = fopen(filename, "rb");

    } else { 

      assert(0);

    }

  } 

  file_map[*fh] = fp;

  return MPI_SUCCESS;
} 


int MPI_File_read(MPI_File fh, void *buf, int count,
    MPI_Datatype datatype, MPI_Status *status) { 
  

  if ( mpi_rank == 0) 
    cout << " intercepted a read call... " << endl;

  FILE * fp; 
  map<MPI_File,FILE*>::iterator it = file_map.find(fh);
  assert( it != file_map.end());

  fp = it->second;
  assert( fp != NULL); 

  int nbytes = count;
  if ( datatype == MPI_INT) 
    nbytes *= sizeof(int);
  else if ( datatype == MPI_DOUBLE) 
    nbytes *= sizeof(double);
  else if ( datatype == MPI_FLOAT)
    nbytes *= sizeof(float);
  else
    assert(0);

  fread(buf,sizeof(char),nbytes,fp);

  return MPI_SUCCESS;
} 

int MPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void *buf,
    int count, MPI_Datatype datatype, MPI_Status * status) { 

  if ( mpi_rank == 0) 
    cout << " intercepted a read at all call... " << endl;

  int * send_count = new int[mpi_size];
  int * send_disp  = new int[mpi_size];
  
  MPI_Gather(&count,1,MPI_INT,send_count,1,MPI_INT,0,mpi_comm);

  FILE * fp; 
  map<MPI_File,FILE*>::iterator it = file_map.find(fh);
  assert( it != file_map.end());

  fp = it->second;

  int count_global = 0;
  int nbytes;
  void* tmp_buf = NULL;
  
  if ( mpi_rank == 0) { 

    for (int rank = 0; rank < mpi_size; ++rank) 
      count_global += send_count[rank];
    
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank) 
      send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

    nbytes = count_global;
    if ( datatype == MPI_INT) 
      nbytes *= sizeof(int);
    else if ( datatype == MPI_DOUBLE) 
      nbytes *= sizeof(double);
    else if ( datatype == MPI_FLOAT)
      nbytes *= sizeof(float);
    else
      assert(0);

    tmp_buf = (void*)(new char[nbytes]);
  } 

  if ( mpi_rank == 0) { 

    assert( fp != NULL); 
    fread(tmp_buf,sizeof(char),nbytes,fp);

  }

  MPI_Scatterv(tmp_buf,send_count,send_disp,datatype,buf,count,datatype,0,mpi_comm);

  delete[] send_count;
  delete[] send_disp;
  if ( tmp_buf != NULL) { 
    
    char* tmp = (char*)(tmp_buf);
    delete[] tmp;

  }

  return MPI_SUCCESS;

} 

int MPI_File_read_at(MPI_File fh, MPI_Offset offset, void *buf,
    int count, MPI_Datatype datatype, MPI_Status * status) { 

  if ( mpi_rank == 0) 
    cout << " intercepted a read at call... " << endl;

  FILE * fp; 
  map<MPI_File,FILE*>::iterator it = file_map.find(fh);
  assert( it != file_map.end());

  fp = it->second;
  assert( fp != NULL); 

  int nbytes; 
  if ( datatype == MPI_INT) 
    nbytes = sizeof(int);
  else if ( datatype == MPI_DOUBLE) 
    nbytes = sizeof(double);
  else if ( datatype == MPI_FLOAT)
    nbytes = sizeof(float);
  else
    assert(0);
  nbytes *= count;

  fseek(fp,offset,SEEK_SET);
  fread(buf,sizeof(char),nbytes,fp);

  return MPI_SUCCESS;
} 

int MPI_File_write(MPI_File fh, const void *buf, int count,
    MPI_Datatype datatype, MPI_Status * status) { 

  if ( mpi_rank == 0) 
    cout << " intercepted a write call... " << endl;

  FILE * fp; 
  map<MPI_File,FILE*>::iterator it = file_map.find(fh);
  assert( it != file_map.end());

  fp = it->second;
  assert( fp != NULL); 

  int nbytes = count;
  if ( datatype == MPI_INT) 
    nbytes *= sizeof(int);
  else if ( datatype == MPI_DOUBLE) 
    nbytes *= sizeof(double);
  else if ( datatype == MPI_FLOAT)
    nbytes *= sizeof(float);
  else
    assert(0);


  fwrite((void*)buf,sizeof(char),nbytes,fp);

  return MPI_SUCCESS;

} 

int MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, const void *buf,
    int count, MPI_Datatype datatype, MPI_Status * status) { 


  if ( mpi_rank == 0) 
    cout << " intercepted a write at all call... " << endl;


  int * recv_count = new int[mpi_size];
  int * recv_disp  = new int[mpi_size];
  
  MPI_Gather(&count,1,MPI_INT,recv_count,1,MPI_INT,0,mpi_comm);

  FILE * fp; 
  map<MPI_File,FILE*>::iterator it = file_map.find(fh);
  assert( it != file_map.end());

  fp = it->second;

  int count_global = 0;
  int nbytes;
  void* tmp_buf = NULL;
  
  if ( mpi_rank == 0) { 

    for (int rank = 0; rank < mpi_size; ++rank) 
      count_global += recv_count[rank];
    
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank) 
      recv_disp[rank] = recv_disp[rank-1] + recv_count[rank-1];

    nbytes = count_global;
    if ( datatype == MPI_INT) 
      nbytes *= sizeof(int);
    else if ( datatype == MPI_DOUBLE) 
      nbytes *= sizeof(double);
    else if ( datatype == MPI_FLOAT)
      nbytes *= sizeof(float);
    else
      assert(0);

    tmp_buf = (void*)(new char[nbytes]);
  } 


  MPI_Gatherv(buf,count,datatype,tmp_buf,recv_count,recv_disp,datatype,0,mpi_comm);


  if ( mpi_rank == 0) { 

    assert( fp != NULL);
    fwrite(tmp_buf,sizeof(char),nbytes,fp);

  }

  delete[] recv_count;
  delete[] recv_disp;
  if ( tmp_buf != NULL) { 
    
    char* tmp = (char*)(tmp_buf);
    delete[] tmp;

  }

  return MPI_SUCCESS;

} 

int MPI_File_write_at(MPI_File fh, MPI_Offset offset, const void *buf,
    int count, MPI_Datatype datatype, MPI_Status * status) { 

  if ( mpi_rank == 0) 
    cout << " intercepted a write at call... " << endl;

  FILE * fp; 
  map<MPI_File,FILE*>::iterator it = file_map.find(fh);
  assert( it != file_map.end());

  fp = it->second;
  assert(fp != NULL);

  int nbytes; 
  if ( datatype == MPI_INT) 
    nbytes = sizeof(int);
  else if ( datatype == MPI_DOUBLE) 
    nbytes = sizeof(double);
  else if ( datatype == MPI_FLOAT)
    nbytes = sizeof(float);
  else
    assert(0);
  nbytes *= count;

  fseek(fp,offset,SEEK_SET);
  fwrite((void*)buf,sizeof(char),nbytes,fp);

  return MPI_SUCCESS;

} 

int MPI_File_close(MPI_File * fh) { 

  if ( mpi_rank == 0) 
    cout << " intercepted a close call... " << endl;

  FILE * fp; 
  map<MPI_File,FILE*>::iterator it = file_map.find(*fh);
  assert( it != file_map.end());

  fp = it->second;

  if ( fp != NULL) { 
    fclose(fp);
  }

  return MPI_SUCCESS;
} 


#endif

int main(int argc, char* argv[]) { 

  //=============================================
  // configuration options 
  // imode_write = {0,1} 
  //             = 0    means a single rank writes
  //             = 1    means everyone writes a chunk
  //             = 2    means rank writes at an offset
  // imode_read  = {0,1}
  //             = 0    means a single rank reads
  //             = 1    means everyone reads a chunk
  //             = 2    means rank reads at an offset
  //=============================================

  int imode_write = 1; 
  int imode_read  = 1;

  if ( argc == 3) {
    imode_write = atoi(argv[1]);
    imode_read = atoi(argv[2]);
  }

  checkMpiError(MPI_Init(&argc,&argv));

  //MPI_Comm mpi_comm;
  //int mpi_rank, mpi_size;

  checkMpiError(MPI_Comm_dup(MPI_COMM_WORLD,&mpi_comm));
  checkMpiError(MPI_Comm_rank(mpi_comm,&mpi_rank));
  checkMpiError(MPI_Comm_size(mpi_comm,&mpi_size));

  //cout << " i am : " << mpi_rank << " of " << mpi_size << endl;
  //cout.flush();


  // make a buffer that everyone has access to...

  int n = 10000; 
  int * buf = new int[n];
  for (int i = 0; i < n; ++i) 
    buf[i] = i;

  
  // write the buf to the disk (either everyone or just rank 0..)

  MPI_File fh;

  checkMpiError(
      MPI_File_open(mpi_comm, "test_file.dat", MPI_MODE_CREATE | MPI_MODE_WRONLY,
        MPI_INFO_NULL, &fh));
  
  if ( imode_write == 0) { 

    if ( mpi_rank == 0) 
      cout << " > writing from rank 0 . " << endl;

    if ( mpi_rank == 0) 
      checkMpiError(MPI_File_write(fh,buf,n,MPI_INT,MPI_STATUS_IGNORE));

  } else if (imode_write == 1) { 

    if ( mpi_rank == 0) 
      cout << " > collectively writing . " << endl;

    int * iora = NULL; calcUniformDist(iora,n,mpi_size);
    int my_n   = iora[mpi_rank+1] - iora[mpi_rank];

    MPI_Offset offset = sizeof(int) * iora[mpi_rank];

    checkMpiError( 
        MPI_File_write_at_all(fh,offset,buf+iora[mpi_rank],my_n,
          MPI_INT,MPI_STATUS_IGNORE));
    
    delete[] iora;

  } else {

    assert( imode_write == 2);

    if (mpi_rank == 0)
      cout << " > writing from rank 0 at an offset. " << endl;

    if (mpi_rank == 0)
      checkMpiError(MPI_File_write_at(fh,0,buf,n,MPI_INT,MPI_STATUS_IGNORE));

  }

  checkMpiError(MPI_File_close(&fh));
  
  MPI_Barrier(mpi_comm);

  if ( mpi_rank == 0) 
    cout << " Successfully wrote test_file.dat" << endl;

  MPI_Barrier(mpi_comm);

  
  // clobber the buffer..

  for (int i =0; i < n; ++i) 
    buf[i] = -1;

  checkMpiError(
      MPI_File_open(mpi_comm, "test_file.dat", MPI_MODE_RDONLY,
        MPI_INFO_NULL, &fh));
  
  if ( imode_read == 0) { 

    if ( mpi_rank == 0) 
      cout << " > reading from rank 0. " << endl;

    if ( mpi_rank == 0) 
      checkMpiError( MPI_File_read(fh,buf,n,MPI_INT,MPI_STATUS_IGNORE));

    // broadcast the info to everyone .. 
    
    checkMpiError(MPI_Bcast(buf,n,MPI_INT,0,mpi_comm));

  } else if (imode_read == 1) { 

    if ( mpi_rank == 0) 
      cout << " > collectively reading. " << endl;

    int * iora = NULL; calcUniformDist(iora,n,mpi_size);
    int my_n   = iora[mpi_rank+1] - iora[mpi_rank];

    MPI_Offset offset = sizeof(int) * iora[mpi_rank];

    int * tmp_buf = new int[mpi_size*my_n];

    checkMpiError( 
        MPI_File_read_at_all(fh,offset,tmp_buf,my_n,
          MPI_INT,MPI_STATUS_IGNORE));

    // a little bit of a shuffle.

    int * send_count = new int[mpi_size];
    for (int rank = 0; rank < mpi_size; ++rank)
      send_count[rank] = my_n;

    int * send_disp = new int[mpi_size];
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank) 
      send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

    for (int rank = 1; rank < mpi_size; ++rank) { 
      for (int i =0; i < my_n; ++i) 
        tmp_buf[rank*my_n+i] = tmp_buf[i];
    }

    int * recv_count = new int[mpi_size];
    for (int rank = 0; rank < mpi_size; ++rank) 
      recv_count[rank] = iora[rank+1] - iora[rank];

    int * recv_disp  = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank) 
      recv_disp[rank] = recv_disp[rank-1] + recv_count[rank-1];

    checkMpiError( 
        MPI_Alltoallv(tmp_buf, send_count, send_disp, MPI_INT,
                      buf, recv_count, recv_disp, MPI_INT, mpi_comm));


    delete[] send_count;
    delete[] recv_count;
    delete[] send_disp;
    delete[] recv_disp;
    delete[] tmp_buf;
    delete[] iora;
  
  } else {

    assert( imode_read == 2);

    if ( mpi_rank == 0) 
      cout << " > reading from rank 0 at an offset. " << endl;

    if ( mpi_rank == 0) 
      checkMpiError( MPI_File_read_at(fh,0,buf,n,MPI_INT,MPI_STATUS_IGNORE));

    // broadcast the info to everyone .. 
    
    checkMpiError(MPI_Bcast(buf,n,MPI_INT,0,mpi_comm));

  }
  
  checkMpiError( MPI_File_close(&fh));
  
  MPI_Barrier(mpi_comm);

  if ( mpi_rank == 0) 
    cout << " Successfully read test_file.dat" << endl;


  MPI_Barrier(mpi_comm);

  
  for (int i = 0; i < n; ++i) 
    assert( buf[i] == i);


  MPI_Barrier(mpi_comm);

  if ( mpi_rank == 0) 
    cout << " Buffer passes check. " << endl;

  MPI_Barrier(mpi_comm);


  delete[] buf;

  checkMpiError(MPI_Comm_free(&mpi_comm));
  checkMpiError(MPI_Finalize());

  return 0;

} 
