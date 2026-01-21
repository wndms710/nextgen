#ifndef BYTESWAP_HPP
#define BYTESWAP_HPP

#include "Common.hpp"
#include <assert.h>

namespace ByteSwap {

  inline int byteSwap(const int old_value) {
    int new_value = 0;
    unsigned char * cold = (unsigned char *)(&old_value);
    unsigned char * cnew = (unsigned char *)(&new_value);
    cnew[0] = cold[3];
    cnew[1] = cold[2];
    cnew[2] = cold[1];
    cnew[3] = cold[0];
    return (new_value);
  }

  inline float byteSwap(const float old_value) {
    float new_value = 0;
    unsigned char * cold = (unsigned char *)(&old_value);
    unsigned char * cnew = (unsigned char *)(&new_value);
    cnew[0] = cold[3];
    cnew[1] = cold[2];
    cnew[2] = cold[1];
    cnew[3] = cold[0];
    return (new_value);
  }

  inline int8 byteSwap(const int8 old_value) {
    int8 new_value = 0;
    unsigned char * cold = (unsigned char *)(&old_value);
    unsigned char * cnew = (unsigned char *)(&new_value);
    cnew[0] = cold[7];
    cnew[1] = cold[6];
    cnew[2] = cold[5];
    cnew[3] = cold[4];
    cnew[4] = cold[3];
    cnew[5] = cold[2];
    cnew[6] = cold[1];
    cnew[7] = cold[0];
    return (new_value);
  }

  inline double byteSwap(const double old_value) {
    double new_value = 0;
    unsigned char * cold = (unsigned char *)(&old_value);
    unsigned char * cnew = (unsigned char *)(&new_value);
    cnew[0] = cold[7];
    cnew[1] = cold[6];
    cnew[2] = cold[5];
    cnew[3] = cold[4];
    cnew[4] = cold[3];
    cnew[5] = cold[2];
    cnew[6] = cold[1];
    cnew[7] = cold[0];
    return (new_value);
  }

  template<typename T>
  inline T byteSwap(const T old_value) {
    T new_value = 0;
    unsigned char * cold = (unsigned char *)(&old_value);
    unsigned char * cnew = (unsigned char *)(&new_value);
    for (int index=0; index < int(sizeof(T)); ++index) {
      cnew[index] = cold[int(sizeof(T))-1-index];
    }
    return (new_value);
  }

  template<typename T>
  inline void byteSwap(T* ptr, const int n) {
    for (int i =0; i < n ; ++i)
      ptr[i] = byteSwap(ptr[i]);
  }

  template<typename T, int K>
  inline void byteSwap( T (*ptr)[K], const int n) {
    for (int i =0; i < n; ++i)
      for (int k=0; k < K; ++k)
        ptr[i][k] = byteSwap(ptr[i][k]);
  }

  template<typename T>
  inline void byteSwap( T (*ptr)[2], const int n) {
    byteSwap<T,2>(ptr,n);
  }

  template<typename T>
  inline void byteSwap( T (*ptr)[3], const int n) {
    byteSwap<T,3>(ptr,n);
  }

  template<typename T>
  inline void byteSwap(T (*ptr)[3][3], const int n) {
    for (int i =0; i < n ; ++i)
      for (int j =0; j < 3 ; ++j)
        for (int k=0; k <3 ; ++k)
          ptr[i][j][k] = byteSwap(ptr[i][j][k]);
  }

  inline void byteSwapHeader(Header * header, const int n) {
    for (int i = 0; i < n; i++) {
      header[i].id = byteSwap(header[i].id);
      header[i].skip = byteSwap(header[i].skip);
      byteSwap(header[i].idata, 16);
      byteSwap(header[i].rdata, 12);
      byteSwap(header[i].ui8data, 4);
    }
  }

  // note that these 2 routines split an int8 into 2 int's using a most-significant
  // and least-significant word. Strictly speaking, these are not most and least
  // significant words because the lsw is only considered to be the bottom 31 bits
  // rather than 32 (i.e. we don't use the sign).
  inline int8 getInt8FromLswMswPair(const int lswmsw[]) {

    int8 lsw = (int8)lswmsw[0];
    int8 msw = (int8)lswmsw[1];

    // do that flipping thing for the case where restart file writing
    // wasn't fixed yet -- i.e. JWN cases in Oct/Nov 2011...
    if (lsw < 0) {
      // hopefully msw is zero...
      assert( msw == 0 );
      lsw += 4294967296ll;
      assert( lsw > 0 );
      return(lsw);
    }

    return( (msw<<31)+lsw );
  }

  inline void setLswMswPairForInt8(int lswmsw[],const int8 value) {

    int8 msw = value >> 31; // bit shifting should get rid of lsw stuff
    int8 lsw = value - (msw<<31);

    lswmsw[0] = (int)lsw;
    lswmsw[1] = (int)msw;
  }
}
#endif
