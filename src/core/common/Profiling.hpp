#ifndef PROFILING_HPP
#define PROFILING_HPP

#include "papi.h"
#include "MpiStuff.hpp"

using namespace MpiStuff;  

namespace PapiProfiling {
  
  static int events[6] = {PAPI_L1_DCM, //L1 data cache miss
                          PAPI_L1_ICM, //L1 inst cache miss
                          PAPI_L2_DCM, //L2 data cache miss
                          PAPI_L2_ICM, //L2 inst cache miss
                          PAPI_L3_DCM, //L3 data cache miss
                          PAPI_L3_ICM};
  
  static const int ncounters = 6;
  long long profile_counts[6];

  inline int papi_start() { 
    if ( PAPI_num_counters() < 6) 
      return -1;
    if ( PAPI_start_counters(events, 6) != PAPI_OK) 
      return -2;
    return 0;
  }

  inline int papi_finish() { 
    if ( PAPI_read_counters(profile_counts, 6) != PAPI_OK)
      return -1;
    return 0;
  }

  void eventToStr(char * filename, const int i) { 
    switch(i) { 
    case 0: 
      sprintf(filename, "L1_data_cache_miss.dat"); 
      break; 
    case 1: 
      sprintf(filename, "L1_inst_cache_miss.dat"); 
      break;
    case 2: 
      sprintf(filename, "L2_data_cache_miss.dat"); 
      break; 
    case 3: 
      sprintf(filename, "L2_inst_cache_miss.dat"); 
      break; 
    case 4: 
      sprintf(filename, "L3_data_cache_miss.dat"); 
      break; 
    case 5: 
      sprintf(filename, "L3_inst_cache_miss.dat"); 
      break; 
    }
  }


};


#endif
