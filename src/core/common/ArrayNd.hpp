// =============================
// N dimentional array operators
// =============================

#ifndef ARRAYND_HPP
#define ARRAYND_HPP

#include <iostream>
using namespace std;
#include <assert.h>


// ======================
// array data structure
// ======================
template<typename T>
class CTIarray {

public:
  
  T * data;
  int n1, n2, n3, n4, n5, n6;
  size_t size;
  int flag;

  bool extern_alloc;

public: 

  // =================
  // 1D
  // =================

  // constructor
  CTIarray(const int n1) {
    assert( n1 > 0 );
    this->n1 = n1;
    this->n2 = -1;
    this->n3 = -1;
    this->n4 = -1;
    this->n5 = -1;
    this->n6 = -1;
    size     = size_t(n1);
    data     = new T[size];
    flag     = 1;
    int tmp  = 0;
    for (unsigned int i = 0; i < size; ++i)
      data[i] = (T)tmp;
    extern_alloc = false;
  }  

  // array acces
  inline T& operator()(const int i1) {
    int index = i1;
    assert( flag == 1 );
#ifdef CTI_DEBUG
    assert( (i1 >= 0)    && (i1 < n1)      );
    assert( (index >= 0) && (index < size) );
#endif
    return( data[index] );
  }
  
  inline T& operator()(const int i1) const {
    int index = i1;
    assert( flag == 1 );
#ifdef CTI_DEBUG
    assert( (i1 >= 0)    && (i1 < n1)      );
    assert( (index >= 0) && (index < size) );
#endif   
    return( data[index] );
  }


  

  // =================
  // 2D
  // =================

  // constructor
  CTIarray(const int n1, const int n2) {
    assert( n1 > 0 );
    assert( n2 > 0 );
    this->n1 = n1;
    this->n2 = n2;
    this->n3 = -1;
    this->n4 = -1;
    this->n5 = -1;
    this->n6 = -1;
    size     = size_t(n1)*n2;
    data     = new T[size];
    flag     = 2;
    int tmp  = 0;
    for (unsigned int i = 0; i < size; ++i)
      data[i] = (T)tmp;
    extern_alloc = false;
  }  

  // array acces
  inline T& operator()(const int i1, const int i2) {
    int index = i1*n2 + i2;
    assert( flag == 2 );
#ifdef CTI_DEBUG
    assert( (i1 >= 0)    && (i1 < n1)      );
    assert( (i2 >= 0)    && (i2 < n2)      );
    assert( (index >= 0) && (index < size) );
#endif
    return( data[index] );
  }
  
  inline T& operator()(const int i1, const int i2) const {
    int index = i1*n2 + i2;
    assert( flag == 2 );
#ifdef CTI_DEBUG
    assert( (i1 >= 0)    && (i1 < n1)      );
    assert( (i2 >= 0)    && (i2 < n2)      );
    assert( (index >= 0) && (index < size) );
#endif
    return( data[index] );
  }



  // =================
  // 3D
  // =================

  // constructor
  CTIarray(const int n1, const int n2, const int n3) {
    assert( n1 > 0 );
    assert( n2 > 0 );
    assert( n3 > 0 );
    this->n1 = n1;
    this->n2 = n2;
    this->n3 = n3;
    this->n4 = -1;
    this->n5 = -1;
    this->n6 = -1;
    size     = size_t(n1)*n2*n3;
    data     = new T[size];
    flag     = 3;
    int tmp  = 0;
    for (unsigned int i = 0; i < size; ++i)
      data[i] = (T)tmp;
    extern_alloc = false;
  }  

  // array acces
  T& operator()(const int i1, const int i2, const int i3) {
    int index = (i1*n2 + i2)*n3 + i3;
    assert( flag == 3 );
#ifdef CTI_DEBUG
    assert( (i1 >= 0)    && (i1 < n1)      );
    assert( (i2 >= 0)    && (i2 < n2)      );
    assert( (i3 >= 0)    && (i3 < n3)      );
    assert( (index >= 0) && (index < size) );
#endif
    return( data[index] );
  }

// array acces
  T& operator()(const int i1, const int i2, const int i3) const {
    int index = (i1*n2 + i2)*n3 + i3;
    assert( flag == 3 );
#ifdef CTI_DEBUG
    assert( (i1 >= 0)    && (i1 < n1)      );
    assert( (i2 >= 0)    && (i2 < n2)      );
    assert( (i3 >= 0)    && (i3 < n3)      );
    assert( (index >= 0) && (index < size) );
#endif
    return( data[index] );
  }




  // =================
  // 4D
  // =================

  // constructor
  CTIarray(const int n1, const int n2, const int n3, const int n4) {
    assert( n1 > 0 );
    assert( n2 > 0 );
    assert( n3 > 0 );
    assert( n4 > 0 );
    this->n1 = n1;
    this->n2 = n2;
    this->n3 = n3;
    this->n4 = n4;
    this->n5 = -1;
    this->n6 = -1;
    size     = size_t(n1)*n2*n3*n4;
    data     = new T[size];
    flag     = 4;
    int tmp  = 0;
    for (unsigned int i = 0; i < size; ++i)
      data[i] = (T)tmp;
    extern_alloc = false;
  }  
 
  CTIarray(T * data, const int n1, const int n2, const int n3, const int n4) {
    assert( n1 > 0 );
    assert( n2 > 0 );
    assert( n3 > 0 );
    assert( n4 > 0 );
    this->n1 = n1;
    this->n2 = n2;
    this->n3 = n3;
    this->n4 = n4;
    size     = size_t(n1)*n2*n3*n4;
    this->data = data;
    flag     = 4;
    extern_alloc = true;
  }  

  // array acces
  T& operator()(const int i1, const int i2, const int i3, const int i4) {
    int index = ((i1*n2 + i2)*n3 + i3)*n4 + i4;
    assert( flag == 4 );
#ifdef CTI_DEBUG
    assert( (i1 >= 0)    && (i1 < n1)      );
    assert( (i2 >= 0)    && (i2 < n2)      );
    assert( (i3 >= 0)    && (i3 < n3)      );
    assert( (i4 >= 0)    && (i4 < n4)      );
    assert( (index >= 0) && (index < size) );
#endif
    return( data[index] );
  }

// array acces
  T& operator()(const int i1, const int i2, const int i3, const int i4) const {
    int index = ((i1*n2 + i2)*n3 + i3)*n4 + i4;
    assert( flag == 4 );
#ifdef CTI_DEBUG
    assert( (i1 >= 0)    && (i1 < n1)      );
    assert( (i2 >= 0)    && (i2 < n2)      );
    assert( (i3 >= 0)    && (i3 < n3)      );
    assert( (i4 >= 0)    && (i4 < n4)      );
    assert( (index >= 0) && (index < size) );
#endif
    return( data[index] );
  }




  // =================
  // 5D
  // =================

  // constructor
  CTIarray(const int n1, const int n2, const int n3, const int n4, const int n5) {
    assert( n1 > 0 );
    assert( n2 > 0 );
    assert( n3 > 0 );
    assert( n4 > 0 );
    assert( n5 > 0 );
    this->n1 = n1;
    this->n2 = n2;
    this->n3 = n3;
    this->n4 = n4;
    this->n5 = n5;
    this->n6 = -1;
    size     = size_t(n1)*n2*n3*n4*n5;
    data     = new T[size];
    flag     = 5;
    int tmp  = 0;
    for (int i = 0; i < size; ++i)
      data[i] = (T)tmp;
    extern_alloc = false;
  }  

  // array acces
  T& operator()(const int i1, const int i2, const int i3, const int i4, const int i5) {
    int index = (((i1*n2 + i2)*n3 + i3)*n4 + i4)*n5 + i5;
    assert( flag == 5 );
#ifdef CTI_DEBUG
    assert( (i1 >= 0)    && (i1 < n1)      );
    assert( (i2 >= 0)    && (i2 < n2)      );
    assert( (i3 >= 0)    && (i3 < n3)      );
    assert( (i4 >= 0)    && (i4 < n4)      );
    assert( (i5 >= 0)    && (i5 < n5)      );
    assert( (index >= 0) && (index < size) );
#endif
    return( data[index] );
  }

 T& operator()(const int i1, const int i2, const int i3, const int i4, const int i5) const {
    int index = (((i1*n2 + i2)*n3 + i3)*n4 + i4)*n5 + i5;
    assert( flag == 5 );
#ifdef CTI_DEBUG
    assert( (i1 >= 0)    && (i1 < n1)      );
    assert( (i2 >= 0)    && (i2 < n2)      );
    assert( (i3 >= 0)    && (i3 < n3)      );
    assert( (i4 >= 0)    && (i4 < n4)      );
    assert( (i5 >= 0)    && (i5 < n5)      );
    assert( (index >= 0) && (index < size) );
#endif
    return( data[index] );
  }

  
  
  // =================
  // 6D
  // =================

  // constructor
  CTIarray(const int n1, const int n2, const int n3, const int n4, const int n5, const int n6) {
    assert( n1 > 0 );
    assert( n2 > 0 );
    assert( n3 > 0 );
    assert( n4 > 0 );
    assert( n5 > 0 );
    assert( n6 > 0 );
    this->n1 = n1;
    this->n2 = n2;
    this->n3 = n3;
    this->n4 = n4;
    this->n5 = n5;
    this->n6 = n6;
    size     = size_t(n1)*n2*n3*n4*n5*n6;
    data     = new T[size];
    flag     = 6;
    int tmp  = 0;
    for (int i = 0; i < size; ++i)
      data[i] = (T)tmp;
    extern_alloc = false;
  }  

  // array acces
  T& operator()(const int i1, const int i2, const int i3, const int i4, const int i5, const int i6) {
    int index = ((((i1*n2 + i2)*n3 + i3)*n4 + i4)*n5 + i5)*n6 + i6;
    assert( flag == 6 );
#ifdef CTI_DEBUG
    assert( (i1 >= 0)    && (i1 < n1)      );
    assert( (i2 >= 0)    && (i2 < n2)      );
    assert( (i3 >= 0)    && (i3 < n3)      );
    assert( (i4 >= 0)    && (i4 < n4)      );
    assert( (i5 >= 0)    && (i5 < n5)      );
    assert( (i6 >= 0)    && (i6 < n6)      );
    assert( (index >= 0) && (index < size) );
#endif
    return( data[index] );
  }

  T& operator()(const int i1, const int i2, const int i3, const int i4, const int i5, const int i6) const {
    int index = ((((i1*n2 + i2)*n3 + i3)*n4 + i4)*n5 + i5)*n6 + i6;
    assert( flag == 6 );
#ifdef CTI_DEBUG
    assert( (i1 >= 0)    && (i1 < n1)      );
    assert( (i2 >= 0)    && (i2 < n2)      );
    assert( (i3 >= 0)    && (i3 < n3)      );
    assert( (i4 >= 0)    && (i4 < n4)      );
    assert( (i5 >= 0)    && (i5 < n5)      );
    assert( (i6 >= 0)    && (i6 < n6)      );
    assert( (index >= 0) && (index < size) );
#endif
    return( data[index] );
  }

  // copy constructor
  CTIarray(const CTIarray& other) {
    n1   = other.n1;
    n2   = other.n2;
    n3   = other.n3;
    n4   = other.n4;
    n5   = other.n5;
    n6   = other.n6;
    size = other.size;
    flag = other.flag;
    extern_alloc = other.extern_alloc;

    if (extern_alloc) {
      data = other.data;
    }
    else {
      data = new T[size];
      for (unsigned int i = 0; i < size; ++i)
	data[i] = other.data[i];
    }
  }

  // array size
  size_t getSize() {
    return(size);
  }

  // array dimention
  int getDim() {
    return(flag);
  }
  
  // destructor...  
  ~CTIarray(){
    if (extern_alloc) {
      // can't do anything here!
    }
    else {
      DELETE(data) ; 
    }
  }
  

private:
  void placeholder();

};





#endif
