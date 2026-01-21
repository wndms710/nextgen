// =========================================
//
// hold
//
// A parallel placeholder for facilitating interactive
// behavior in a queueing system.
//
// =========================================
 
#include <iostream>
#include <unistd.h> // it is used for usleep 
#include <mpi.h>
#include "Logger.hpp"

using namespace std;

int main(int argc,char * argv[]) {

  Logger * logger = NULL;

  MPI_Init(&argc,&argv);
 
  try {

    MPI_Comm mpi_comm;
    int mpi_rank,mpi_size;

    MPI_Comm_dup(MPI_COMM_WORLD, &mpi_comm);
    MPI_Comm_rank(mpi_comm, &mpi_rank);
    MPI_Comm_size(mpi_comm, &mpi_size);
 
    if (mpi_rank==0){
      int iarg0 = 1;
      while (iarg0 < argc) {
        if (strcmp(argv[iarg0],"--with-logging") == 0) {
          cout << "Found --with-logging, enabling central logging" << endl;
#ifdef WITH_SQLITE
          //solver_name,ncores,
          logger = new LoggerSqlite(argv[0],mpi_size);
#else
          logger = new LoggerAscii(argv[0],mpi_size); 
#endif
          logger->setKillFilename("killhold");
          logger->insertSolverAction(SIMULATION_INIT);
        }
        ++iarg0;
      }
    }
    
    if (mpi_rank == 0)
      cout << "hold started on " << mpi_size << " processes..." << endl;
 
    if (logger) logger->insertSolverAction(SIMULATION_RUN);

    int done = 0;
    while (done != 1) {
      
      // look for killhold to stop...
      if (mpi_rank == 0) {
	cout << " > checking for killhold..." << endl;
	const int ierr = MPI_File_delete("killhold",MPI_INFO_NULL);
	if (ierr == 0)
	  done = 1;
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

      // do something costly...
      if (done != 1) 
	usleep(5000000);
    }
  } catch(...) {
    cerr << "unhandled Exception.\n" << endl;
    MPI_Finalize();
    return(-1);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  if (logger) {
    logger->insertSolverAction(SIMULATION_FINAL);
    delete logger;
    logger = NULL;
  }

  return 0;
}

