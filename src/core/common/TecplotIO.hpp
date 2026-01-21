#ifndef TECPLOTIO_HPP
#define TECPLOTIO_HPP

#include <iostream>
#include <assert.h>
#include <string.h>

namespace TecplotIO {
  
  inline int writeIntsTecplot(char * cbuf,const int * ibuf,const int n) {
    assert(sizeof(int) == 4);
    memcpy(cbuf,ibuf,4*n);
    return( 4*n );
  }
  
  inline int writeFloatsTecplot(char * cbuf,const float * fbuf,const int n) {
    assert(sizeof(float) == 4);
    memcpy(cbuf,fbuf,4*n);
    return( 4*n );
  }
  
  inline int writeCharsTecplot(char * cbuf, const char * message) {
    int n = strlen(message);
    for (int i = 0; i < n; ++i) {
      int ichar = (int)message[i];
      memcpy(cbuf+4*i,&ichar,4);
    }
    int ichar = 0;
    memcpy(cbuf+4*n,&ichar,4);
    return(4*(n+1));
  }
  
  inline int writeCharsTecplot(char * cbuf, const std::string& message) { return writeCharsTecplot(cbuf, message.c_str()); }

  inline int writeCharsTecplot(char * cbuf, const char * message1,const char * message2) {
    int n1 = strlen(message1);
    for (int i = 0; i < n1; ++i) {
      int ichar = (int)message1[i];
      memcpy(cbuf+4*i,&ichar,4);
    }
    int n2 = strlen(message2);
    for (int i = 0; i < n2; ++i) {
      int ichar = (int)message2[i];
      memcpy(cbuf+4*(n1+i),&ichar,4);
    }
    int ichar = 0;
    memcpy(cbuf+4*(n1+n2),&ichar,4);
    return(4*(n1+n2+1));
  }
  
  inline int writeCharsTecplot(char * cbuf, const std::string& message1,const char * message2) { return writeCharsTecplot(cbuf, message1.c_str(), message2); }
}
#endif
