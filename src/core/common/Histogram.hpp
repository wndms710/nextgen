#ifndef _HISTOGRAM_HPP_
#define _HISTOGRAM_HPP_

#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <assert.h>

using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::min;
using std::vector;

#include "MpiStuff.hpp"

class Histogram {

private:

  int nbin;
  double * my_count;
  double * count;
  double * cdf;
  double * bins;
  double x0,x1;
  bool rangeSet;
  bool logFlag;
  bool intBins;
  MPI_Comm comm;

public:

  Histogram() {
    _initEmpty();
    comm = MpiStuff::mpi_comm;
  }

  Histogram(MPI_Comm _comm) {
    _initEmpty();
    comm = _comm;
  }

  void _initEmpty() {
    nbin     = 100;
    my_count = NULL;
    count    = NULL;
    cdf      = NULL;
    bins     = NULL;
    intBins  = false;
    logFlag  = false;
    rangeSet = false;
  }

  Histogram(MPI_Comm _comm,const int _nbin) {
    _initEmptyWithBin(_nbin);
    comm = _comm;
  }

  void _initEmptyWithBin(const int _nbin) {
    nbin     = _nbin;
    my_count = NULL;
    count    = NULL;
    cdf      = NULL;
    bins     = NULL;
    intBins  = false;
    logFlag  = false;
    rangeSet = false;
  }

  Histogram(const double x0,const double x1) {
    _initRange(x0,x1);
    comm = MpiStuff::mpi_comm;
  }

  Histogram(const double x0, const double x1, MPI_Comm _comm) {
    _initRange(x0,x1);
    comm = _comm;
  }

  void _initRange(const double _x0, const double _x1) {
    nbin     = 100;
    rangeSet = true;
    this->x0 = _x0;
    this->x1 = _x1;
    my_count = NULL;
    count    = NULL;
    cdf      = NULL;
    bins     = NULL;
    intBins  = false;
    logFlag  = false;
  }

  // Histogram(const double x0,const double x1,const int nbin) {
  //   _initRangeBin(x0,x1,nbin);
  //   comm = MpiStuff::mpi_comm;
  // }

  Histogram(const double x0,const double x1, const int nbin, MPI_Comm _comm) {
    _initRangeBin(x0,x1,nbin);
    comm = _comm;
  }

  void _initRangeBin(const double _x0, const double _x1, const int _nbin) {
    this->nbin = _nbin;
    rangeSet   = true;
    this->x0   = _x0;
    this->x1   = _x1;
    my_count   = NULL;
    count      = NULL;
    cdf        = NULL;
    bins       = NULL;
    intBins    = false;
    logFlag    = false;
  }

  ~Histogram() {
    DELETE(my_count);
    DELETE(count);
    DELETE(cdf);
    DELETE(bins);
  }

  void setLog(const bool logFlag) {
    this->logFlag = logFlag;
    if (logFlag && rangeSet)
    x0 = max(x0,1.0E-20); // cannot be zero or negative
  }

  void setIntBins(const bool _bins) {
    this->intBins = _bins;
  }

  void add(const double * var,const int nv) {

    if (!rangeSet) {
      rangeSet = true;
      double my_buf[2] = { 1.0E+20, 1.0E+20 };
      for (int iv = 0; iv < nv; ++iv) {
        my_buf[0] = min(my_buf[0],var[iv]);
        my_buf[1] = min(my_buf[1],-var[iv]);
      }
      double buf[2];
      if ( comm != MPI_COMM_NULL)
      MPI_Allreduce(my_buf,buf,2,MPI_DOUBLE,MPI_MIN,comm);
      else {
        buf[0] = my_buf[0];
        buf[1] = my_buf[1];
      }
      x0 = buf[0];
      x1 = -buf[1] + 1.0E-12*(-buf[1]-buf[0]);

      if ( x1-x0 < 1.0E-12) {
        assert( abs(x1-x0) < 1.0E-12);
        x1 = x0 + 1.0E-12;
      }

      if (logFlag) {
        // make sure range is positive...
        x0 = max(x0,1.0E-20);
      }
      assert( x1 > x0 );
    }

    if (my_count == NULL) {
      my_count = new double[nbin+2];
      for (int i = 0; i < nbin+2; ++i)
      my_count[i] = 0.0;
    }

    for (int iv = 0; iv < nv; ++iv) {
      if (var[iv] < x0)
      my_count[0] += 1.0;
      else if (var[iv] >= x1)
      my_count[nbin+1] += 1.0;
      else if (logFlag) {
        int index = int((log(var[iv])-log(x0))/(log(x1)-log(x0))*double(nbin));
        assert((index >= 0)&&(index < nbin));
        my_count[index+1] += 1.0;
      }
      else {
        int index = int((var[iv]-x0)/(x1-x0)*double(nbin));
        if (!((index >= 0)&&(index < nbin))) {
          cout << "got index: " << index << " for var_v: " << var[iv] << " in range: " << x0 << ":" << x1 << endl;
        }
        else {
          assert((index >= 0)&&(index < nbin));
          my_count[index+1] += 1.0;
        }
      }
    }
  }

  double binOfCdf(const double _val) const {
    if ( _val <= cdf[0]) {
      return bins[0];
    } else if ( _val >= cdf[nbin-1]) {
      return bins[nbin-1];
    } else {
      int r; MiscUtils::findInterval(r,_val,cdf,nbin);
      assert( (cdf[r+1] >= _val)&&(cdf[r] <= _val));
      return bins[r];
    }
  }

  void reduce() {

    int mpi_rank   = 0;
    if ( comm != MPI_COMM_NULL)
    MPI_Comm_rank(comm,&mpi_rank);

    // everyone is going to keep a copy of the final counts...
    if ( count == NULL) count = new double[nbin+2];
    if ( cdf   == NULL) cdf   = new double[nbin];
    if ( bins  == NULL) bins  = new double[nbin];
    if ( comm != MPI_COMM_NULL) {
      MPI_Allreduce(my_count,count,nbin+2,MPI_DOUBLE,MPI_SUM,mpi_comm);
    } else {
      for (int i =0; i < nbin+2; ++i)
      count[i] = my_count[i];
    }

    double count_sum = 0.0;
    for (int i =0; i < nbin+2; ++i)
    count_sum += count[i];

    double this_cdf = count[0];
    for (int i = 0; i < nbin; ++i) {
      if (logFlag)
      bins[i] = exp(log(x0) + (double(i)+0.5)/double(nbin)*(log(x1)-log(x0)));
      else
      bins[i] =  x0 + (double(i)+0.5)/double(nbin)*(x1-x0);

      cdf[i] = (this_cdf + 0.5*count[i+1])/count_sum;
      this_cdf += count[i+1];
    }
  }

  void write(const char* filename) {

    int mpi_rank   = 0;
    if ( comm != MPI_COMM_NULL)
    MPI_Comm_rank(comm,&mpi_rank);

    this->reduce();

    if (mpi_rank == 0) {
      ofstream ofile;
      ofile.open(filename);
      ofile << "# bin count cummulative-count cdf" << endl;
      double count_sum = 0.0;
      for (int i =0; i < nbin+2; ++i)
      count_sum += count[i];

      for (int i = 0; i < nbin; ++i) {
        if (this->intBins) ofile << floor(bins[i]);
        else ofile << bins[i];

        ofile << " " << count[i+1] << " " << cdf[i]*count_sum << " " << cdf[i] << endl;
      }
      ofile.close();
    }
  }
};
#endif
