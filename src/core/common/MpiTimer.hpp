#ifndef MPITIMER_HPP
#define MPITIMER_HPP

#include "MpiStuff.hpp"
using namespace MpiStuff;

class MpiTimer {
public:

  vector<double> buf;
  vector<double> buf_sum;
  vector<string> tags;
  
  void addTag(const string& tag) {
    if ( tag == "") { 
      char arr[16]; sprintf(arr,"%d",int(buf.size())-1);
      tags.push_back(string(arr));
    } else { 
      tags.push_back(tag);
    }
  }
  
  void start(const string tag ="") {
    buf.clear();
    tags.clear();

    // do not disturb the accumulation in buf_sum..
    buf.push_back(MPI_Wtime());
    addTag(tag);
  }
  
  void split(const string tag="") {
    buf.push_back(MPI_Wtime());
    addTag(tag);
  }
  
  void stop() { split(); }
  
  // report the elapsed time from the start to split j
  double elapsed(const int j1, const int j0) const{
    return buf[j1] - buf[j0];
  }
    
  double elapsed(const int j) const { return elapsed(j,0); }
  double elapsed() const { return elapsed(1,0) ; }

  void accumulate() {

    if ( buf_sum.size() != buf.size()) { 
      assert(buf_sum.size() == 0);
      buf_sum.resize(buf.size());
      for (int i =0; i <int(buf.size()); ++i) 
        buf_sum[i] = 0.0;
    }

    const int vec_size = int(buf.size());
    for (int i =1; i < vec_size ; ++i) {
      buf_sum[i-1] += elapsed(i,i-1); 
    }
  }
  
  double norm_interval(const int i) const { 
    double sum = 0.0;
    for (int j =0; j < int(buf_sum.size())-1; ++j) 
      sum += buf_sum[j];
    return buf_sum[i-1]/sum; 
  }
  
  void report() {
    
    // user is responsible for accumulating the data while the calculation runs
    // aggregate the timing data from all of the ranks and report and average.. 
    const int vec_size = int(buf.size());
    if ( vec_size > 1 ) { 

      double * tmp = new double[vec_size-1];
      MPI_Reduce(&(buf_sum[0]),tmp,vec_size-1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      if ( mpi_rank == 0) { 
        
        double sum = 0.0;
        for (int i =0; i < vec_size-1; ++i) 
          sum += tmp[i];
        
        assert( int(tags.size()) == vec_size);
        cout << " ============== timing report ===================== " << endl;
        for (int i =0; i < vec_size-1; ++i) { 
          cout << " > frac time spent in " << tags[i+1] << "  : " << tmp[i]/sum << endl;
        }
        cout << " ================================================== " << endl;
      }
      
      delete[] tmp;

    } else { 
      if ( mpi_rank == 0 )
        cout << " no timing information to report. " << endl;
    }
  }
};
#endif
