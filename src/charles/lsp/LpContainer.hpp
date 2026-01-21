#ifndef __LPCONTAINER_HPP__
#define __LPCONTAINER_HPP__

#include "LpData.hpp"

// =============================
// flag of Lp for recycling
// =============================
//enum {
//  TRASH = 10,
//  KEEP
//};

template<class T>
class LpContainer {
  // ======================================================
  // every class that manages lagrangian particles (lp)
  // should inherit this container class that makes
  // many common routines available.
  // ======================================================
protected:
  int np_max;
public:
  T * lp;
  int np;
  int * lpocv_i;
  LpContainer() {
    lp = NULL;
    np = np_max = 0;
    lpocv_i = NULL;
  }
  virtual ~LpContainer() {
    DELETE(lp);
    DELETE(lpocv_i);
  }
  void resize(int np_new) {
    if (np_new > np_max) {
      T * lp_old = lp; // store the old one -- even if null...
      np_max = np_new + max(256,np_new/5); // pad by 20% of np_new
      lp = new T[np_max];
      if (lp_old != NULL) {
        for (int ip = 0; ip < np; ++ip) {
          lp[ip] = lp_old[ip];
        }
        delete[] lp_old;
      }
    }
    np = np_new;
  }
  int size() const {
    return np;
  }
  int size_max() const {
    return np_max;
  }
  virtual void report() const {
    int8 np_int8 = np;
    int8 np_global;
    MPI_Reduce(&np_int8,&np_global,1,MPI_INT8,MPI_SUM,0,mpi_comm);
    int np_minmax[2] = { np, -np };
    int np_minmax_global[2];
    MPI_Reduce(np_minmax,np_minmax_global,2,MPI_INT,MPI_MIN,0,mpi_comm);
    if (mpi_rank == 0) 
      cout << " > np_global: " << np_global << 
        " min,avg,max per rank: " << np_minmax_global[0] << 
        " " << double(np_global)/double(mpi_size) << 
        " " << -np_minmax_global[1] << endl;
  }

  // this allocates lp for the first time, 
  // assuming np is read from restart or even 0
  void init() {
    assert(lp == NULL);
    np_max = np;
    lp = new T [np];
  }

  int getLpSize() const {

    int my_size = this->size();
    int total_size;
    MPI_Allreduce(&my_size, &total_size, 1, MPI_INT, MPI_SUM, mpi_comm);

    return total_size;
  }

  // recycling routine
  int recycleLp() {
    int ip_new = 0;
    for (int ip = 0; ip < np; ++ip) {
      if (lp[ip].flag == KEEP) {
        if (ip != ip_new) {
          //lp[ip_new].copy(lp[ip]);
          lp[ip_new] = lp[ip];
        }
        ip_new++;
      } else if (lp[ip].flag == TRASH) {
      } else {
        CERR("ERROR:particle flag not recognized" << " , particle: " << ip << " , flag: " << lp[ip].flag << "\n");
      }
    }
    int myNTrash = np - ip_new;
    int NTrash;
    MPI_Reduce(&myNTrash,&NTrash,1,MPI_INT,MPI_SUM,0,mpi_comm);
    resize(ip_new);
    return NTrash;
  }

};

#endif
