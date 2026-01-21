#ifndef COMMON_HPP
#define COMMON_HPP

//void srand(const int seed);

#include <iostream>
#include <cstdio>
#include <cstdlib> 
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <string>
#include <cstring>
#include <vector>
#include <assert.h>
#include <list>
#include <math.h> 
#include <algorithm> 
#include <functional> 
#include <map>
#include <set>
#include <deque>
#include <sys/stat.h> 
#include <errno.h>
#include <unistd.h> 
#include <iomanip>
#include <complex>
#include <utility>
#include <climits>
#include <limits>
#include <stdint.h>

#include "Macros.hpp" 
#include "Defs.hpp"
#include "Structs.hpp"

using namespace std;

// assorted type defs
typedef complex<double> dcomplex;

template <class T1,class T2,class T3>
class triple {
public:
  T1 first;
  T2 second;
  T3 third;
  triple() {} 
  triple(const T1 first,const T2 second,const T3 third) {
    this->first = first;
    this->second = second;
    this->third = third;
  }
};

template<typename T>
inline T getMaxLimit();

template<>
inline int getMaxLimit<int>() {
  return INT_MAX;
}

template<>
inline uint getMaxLimit<uint>() {
  return UINT_MAX;
}

template<>
inline int8 getMaxLimit<int8>() {
  return LLONG_MAX;
}

template<>
inline uint8 getMaxLimit<uint8>() {
  return ULLONG_MAX;
}

template<>
inline float getMaxLimit<float>() {
  return HUGE_VALF;
}

template<>
inline double getMaxLimit<double>() {
  return HUGE_VAL;
}

#endif
