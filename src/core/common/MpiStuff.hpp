#ifndef MPISTUFF_HPP
#define MPISTUFF_HPP

#include "Common.hpp"

#ifdef NO_MPI
# include "NoMpi.hpp"
#else
# include <mpi.h>
#endif

#ifdef INT8_IS_LONG_INT
# define MPI_OFFSET_DATATYPE MPI_LONG
# define MPI_INT8 MPI_LONG
# define MPI_UINT8 MPI_UNSIGNED_LONG
#else
# define MPI_OFFSET_DATATYPE MPI_LONG_LONG
# define MPI_INT8 MPI_LONG_LONG
# define MPI_UINT8 MPI_UNSIGNED_LONG_LONG
#endif

#define MPI_INT2 MPI_SHORT
#define MPI_UINT2 MPI_UNSIGNED_SHORT

#ifdef WITH_SHM
#include <errno.h>     // POSIX error codes
#include <sys/mman.h>  // shm_open, shm_unlink, mmap, memory management defines
#include <fcntl.h>     // file access defines
#include <unistd.h>    // for ftruncate()
#define SHM_NAME "/cti_shm"
#endif

namespace MpiStuff {
  
  // these namespace members are declared "extern" so we can include
  // the namespace definition in a header file included by
  // multiple routines
  extern int mpi_rank;
  extern int mpi_size;
  extern MPI_Comm mpi_comm;
  
  extern double mpi_wtime_0;
  
  // the shared communicator connects cores associated with a specific node. Its
  // size should be the number of cores per node, and does not have to be
  // homgeneous across the resource...
  extern int mpi_rank_shared;
  extern int mpi_size_shared;
  extern MPI_Comm mpi_comm_shared;
  
  // the internode communicator connects all mpi_rank_shared == 0 procs together
  // Its size should be the number of nodes. It is not valid for mpi_rank_shared != 0...
  extern int mpi_rank_internode;
  extern int mpi_size_internode;
  extern MPI_Comm mpi_comm_internode;

  extern int mpi_rank_internode_actual; // known by all ranks on node 
  
  // this is a size mpi_size_internode array that 
  // contains the rank (int mpi_comm) of the first 
  // rank on each node. This is useful when some algorithm 
  // needs to cycle amongst the nodes...
  extern int * rank_of_rank_internode;
  
  // type sizes...
  extern MPI_Offset int_size,
    float_size,
    int8_size,
    uint8_size,
    double_size,
    header_size;
  
  extern MPI_Datatype MPI_Header;

  extern int ishm; // index to help increment unique filenames 
  extern map<void*,SharedMemDescriptor> shm_map;

  //extern int phys_cores_per_node;
  
  void initMpiStuff();
  void initMpiStuff(MPI_Comm& comm);
  void initMpiStuffCommon();
  void destroyMpiStuff();
  
  void MPI_Pause(const char * message);
  void MPI_File_pause(const char * message);
  void MPI_Sync(const char * message,const bool reset = false);
  void MPI_Run_check();
  int MPI_Type_indexed_clean(int n_zones, int *my_zone_count, int *my_zone_disp,
                             MPI_Datatype old_type, MPI_Datatype *new_type);
  int MPI_Type_indexed_clean(int n_zones, int *my_zone_count, int8 *my_zone_disp,
                             MPI_Datatype old_type, MPI_Datatype *new_type);
  void MPI_Check(const int ierr,const char * message);
  int MPI_Bcast_string(std::string& str,int rank,MPI_Comm& comm);
  int MPI_Bcast_stringVec(std::vector<std::string>& stringVec,int rank,MPI_Comm& comm);
  
  inline int getStatusCount(MPI_Status& status, MPI_Datatype type) {
    int count ;
    MPI_Get_count(&status, type, &count);
    return count;
  }
  
  void checkMemoryUsage() ;
  
  template<class T>
  MPI_Datatype getMpiDatatype();
  size_t sizeOfMpiDatatype(MPI_Datatype datatype);
  
  // unfortunately, the following could not be templated, so just use defines instead...

#define T unsigned char
# include "Mmap_Munmap_defs.hpp"
#undef T
  
#define T int
# include "Mmap_Munmap_defs.hpp"
#undef T
  
#define T uint
# include "Mmap_Munmap_defs.hpp"
#undef T
  
#define T float
# include "Mmap_Munmap_defs.hpp"
#undef T

#define T double
# include "Mmap_Munmap_defs.hpp"
#undef T

  /*
    #define T uint32_t
    # include "Mmap_Munmap_defs.hpp"
    #undef T
  */

#define T int8
# include "Mmap_Munmap_defs.hpp"
#undef T

#define T uint8
# include "Mmap_Munmap_defs.hpp"
#undef T

#define T uint2
# include "Mmap_Munmap_defs.hpp"
#undef T

  // parallel-adt partitioning...

  void setPartKind(const string& name); // supports "PADT", "BAD"
  void repartXcvPadt(double (*&x_cv)[3],int8 *&icv_global,int &ncv,const MPI_Comm& mpi_comm,const int bits = 7); // 111 -> all directions
  void repartXcvWeightedPadt(double (*&x_cv)[3],int8 (*&cv_global_wgt)[2],int &ncv,const MPI_Comm& mpi_comm,const int bits = 7);
  void calcCvPartPadt(int * cv_part,const double (*x_cv_orig)[3],const int ncv,const MPI_Comm& mpi_comm,const int bits = 7);

  // and reordering...
  
  void reorderXcv(double (*x_cv)[3],int8 *icv_global,const int ncv);

}

#endif
