#include "CTI.hpp"
using namespace CTI;
#include "MiscUtils.hpp"
using namespace MiscUtils;

int main(int argc,char * argv[]) {
    
  // initialize MPI...
  CTI_Init(argc,argv,"testIO.in");

  const string mode = getStringParam("MODE","ALL");
  const uint8 count = getInt8Param("COUNT",1000000);

  if (mode == "ALL" || mode == "SERIAL_WRITE") {
    if (mpi_rank == 0) {

      cout << " > serial writing... " << endl;
      cout << "   > count: " << count << endl;
      remove("dummy.dat");
      FILE * fp = fopen("dummy.dat","wb");
      assert(fp != NULL);
      double *buf = new double[count];
      for (int i = 0; i < count; ++i) {
        buf[i] = 1.0-2.0*rand()/double(RAND_MAX);
        if (i < 10 && i < count) 
          cout << buf[i] << " ";
      }
      cout << endl;
      fwrite(&count,sizeof(uint8),1,fp);
      fwrite(buf,sizeof(double),count,fp);
      delete[] buf; 
      fclose(fp); 

    }
  }

  MPI_Barrier(mpi_comm);

  if (mode == "ALL" || mode == "PARALLEL_READ_WRITE" || mode == "PARALLEL_READ") {
  
    if (mpi_rank == 0)
      cout << " > parallel reading... " << endl;
    MPI_File fh;
    int ierr = MPI_File_open(mpi_comm,"dummy.dat",MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
    if (ierr)
      cout << " > unable to open dummy.dat on rank " << mpi_rank << endl;
    uint8 count;
    if (mpi_rank == 0) {
      ierr = MPI_File_read_at(fh,0,&count,1,MPI_UINT8,MPI_STATUS_IGNORE);
      if (ierr)
        cout << " > unable to read dummy.dat on rank " << mpi_rank << endl;
    }
    MPI_Bcast(&count,1,MPI_UINT8,0,mpi_comm);
    int8* xora = NULL;
    MiscUtils::calcUniformDist(xora,count,mpi_size);
    assert(xora[mpi_rank+1]-xora[mpi_rank] < TWO_BILLION);
    int my_count = int(xora[mpi_rank+1]-xora[mpi_rank]);
    double* buf = new double[my_count];
    MPI_Offset offset = xora[mpi_rank]*1*8 + 1*8;
    ierr = MPI_File_read_at_all(fh,offset,buf,my_count,MPI_DOUBLE,MPI_STATUS_IGNORE);
    if (ierr)
      cout << " > unable to read dummy.dat on rank " << mpi_rank << endl;
    if (mpi_rank == 0) {
      for (int i = 0; i < min(my_count,10); ++i)
        cout << buf[i] << " ";
      cout << endl;
    }
    ierr = MPI_File_close(&fh);
    fh = NULL;
    if (ierr)
      cout << " > unable to close dummy.dat on rank " << mpi_rank << endl;

    MPI_Barrier(mpi_comm);

    if (mode != "PARALLEL_READ") {

      if (mpi_rank == 0)
        cout << " > parallel writing..." << endl;
      int my_ierr = MPI_File_delete("dummy.dat",MPI_INFO_NULL);
      MPI_Reduce(&my_ierr,&ierr,1,MPI_INT,MPI_MIN,0,mpi_comm);
      if (mpi_rank == 0) {
        if (ierr)
          cout << " > unable to delete dummy.dat" << endl;
      }
      assert(fh == NULL);
      ierr = MPI_File_open(mpi_comm,"dummy.dat",MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);
      if (ierr)
        cout << " > unable to open dummy.dat on rank " << mpi_rank << endl;
      assert(ierr == MPI_SUCCESS);
      if (mpi_rank == 0) {
        ierr = MPI_File_write_at(fh,0,&count,1,MPI_UINT8,MPI_STATUS_IGNORE);
        if (ierr)
          cout << " > unable to write dummy.dat on rank " << mpi_rank << endl;
      }
      ierr = MPI_File_write_at_all(fh,offset,buf,my_count,MPI_DOUBLE,MPI_STATUS_IGNORE);
      if (ierr)
        cout << " > unable to write dummy.dat on rank " << mpi_rank << endl;
      ierr = MPI_File_close(&fh);
      if (ierr)
        cout << " > unable to close dummy.dat on rank " << mpi_rank << endl;
    }
    delete[] xora;
    delete[] buf;

  }

  MPI_Barrier(mpi_comm);

  if (mode == "ALL" || mode == "SERIAL_READ") {
    if (mpi_rank == 0) {

      cout << " > serial reading... " << endl;
      FILE* fp = fopen("dummy.dat","rb");
      assert(fp != NULL);
      uint8 count;
      fread(&count,sizeof(uint8),1,fp);
      cout << "   > count: " << count << endl;
      double *buf = new double[count];
      fread(buf,sizeof(double),count,fp); 
      for (int i = 0; i < min(count,uint8(10)); ++i) {
        cout << buf[i] << " ";
      }
      cout << endl;
      delete[] buf;
      fclose(fp);

    }
  }

  if (mode == "ALL" || mode == "SHM") {
    if (mpi_rank == 0)
      cout << " > checking shm stuff..." << endl;
    double *buf = NULL;
    CTI_Mmap(buf,count);
    if (mpi_rank_shared == 0) {
      for (int i = 0; i < count; ++i) 
        buf[i] = 1.0-2.0*rand()/double(RAND_MAX);
    }
    MPI_Barrier(mpi_comm_shared);
    CTI_Munmap(buf,count);
    buf = NULL;
    if (mpi_rank == 0)
      cout << " > checking shm rw stuff..." << endl;
    CTI_Mmap_rw(buf,count); // _rw
    int* xora = NULL;
    MiscUtils::calcUniformDist(xora,count,mpi_size_shared);
    int my_count = int(xora[mpi_rank_shared+1]-xora[mpi_rank_shared]);
    for (int i = 0; i < my_count; ++i)
      buf[i+xora[mpi_rank_shared]] = 1.0-2.0*rand()/double(RAND_MAX);
    delete[] xora;
    MPI_Barrier(mpi_comm_shared);
    CTI_Munmap(buf,count);
  }

  MPI_Barrier(mpi_comm);

  CTI_Finalize();

  return 0;
}
