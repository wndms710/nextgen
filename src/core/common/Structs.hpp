#ifndef STRUCTS_HPP
#define STRUCTS_HPP

//
// common lightweight data structures
//
typedef struct {
  int i1,i2;
} IntPair;

class IntPairCompare1 : public std::binary_function<IntPair&, IntPair&, bool> {
  public:
    bool operator()(const IntPair& a,const IntPair& b) {
      return( a.i1 < b.i1 );
    }
};

// same as above, but produces a stable sort behavior...
class IntPairCompare12 : public std::binary_function<IntPair&, IntPair&, bool> {
  public:
    bool operator()(const IntPair& a,const IntPair& b) {
      return( (a.i1 < b.i1) || ((a.i1 == b.i1)&&(a.i2 < b.i2)) );
    }
};

typedef struct {
  int8 i1;
  int i2;
} Int8IntPair;

class Int8IntPairCompare1 : public std::binary_function<Int8IntPair&, Int8IntPair&, bool> {
public:
  bool operator()(const Int8IntPair& a,const Int8IntPair& b) {
    return( a.i1 < b.i1 );
  }
};

// same as above, but produces a stable sort behavior...
class Int8IntPairCompare12 : public std::binary_function<Int8IntPair&, Int8IntPair&, bool> {
 public:
  bool operator()(const Int8IntPair& a,const Int8IntPair& b) {
    return( (a.i1 < b.i1) || ((a.i1 == b.i1)&&(a.i2 < b.i2)) );
  }
};

class Int8Pair {
public:
  int8 i1,i2;

  bool operator<(const Int8Pair& _p2) const {
    return ( (i1 < _p2.i1) && (i2 < _p2.i2));
  }
};

template <class T0, class T1, class T2>
class Triple {
public:
  T0 first;
  T1 second;
  T2 third;

  Triple() {}
  Triple(const Triple& other) {
    this->first = other.first;
    this->second = other.second;
    this->third = other.third;
  }
  Triple(const T0 _first,const T1 _second,const T2 _third) {
    this->first = _first;
    this->second = _second;
    this->third = _third;
  }
  ~Triple() {}

  // overload operators so comparisons/sorts work as desired
  inline bool operator < (const Triple& other) const {
    return( (this->first < other.first) ||
            ((this->first == other.first)&&(this->second < other.second)) ||
            ((this->first == other.first)&&(this->second == other.second)&&(this->third < other.third)) );
  }
  inline bool operator == (const Triple& other) const {
    return( (this->first == other.first) && (this->second == other.second) && (this->third == other.third) );
  }
  inline bool operator > (const Triple& other) const {
    return( (this->first > other.first) ||
            ((this->first == other.first)&&(this->second > other.second)) ||
            ((this->first == other.first)&&(this->second == other.second)&&(this->third > other.third)) );
  }
  inline bool operator >= (const Triple& other) const {
    return ((this > other) || (this == other));
  }
  inline bool operator <= (const Triple& other) const {
    return ((this < other) || (this == other));
  }
};

class IntTriple {
public:
  int i1,i2,i3;
  IntTriple() {
  }
  IntTriple(const int i1,const int i2,const int i3) {
    this->i1 = i1;
    this->i2 = i2;
    this->i3 = i3;
  }
  int operator[](const int index) {
    if (index == 2) return this->i3;
    else if (index == 1) return this->i2;
    else return this->i1;
  }
};

class IntTripleCompare1 : public std::binary_function<IntTriple&, IntTriple&, bool> {
 public:
  bool operator()(const IntTriple& a,const IntTriple& b) {
    return( a.i1 < b.i1 );
  }
};

class IntTripleCompare3 : public std::binary_function<IntTriple&, IntTriple&, bool> {
 public:
  bool operator()(const IntTriple& a,const IntTriple& b) {
    return( a.i3 < b.i3 );
  }
};

class IntTripleCompare12 : public std::binary_function<IntTriple&, IntTriple&, bool> {
 public:
  bool operator()(const IntTriple& a,const IntTriple& b) {
    return( (a.i1 < b.i1) || ((a.i1 == b.i1)&&(a.i2 < b.i2)) );
  }
};

class IntTripleCompare23 : public std::binary_function<IntTriple&, IntTriple&, bool> {
 public:
  bool operator()(const IntTriple& a,const IntTriple& b) {
    return( (a.i2 < b.i2) || ((a.i2 == b.i2)&&(a.i3 < b.i3)) );
  }
};

// class IntTripleCompare123 : public std::binary_function<IntTriple&, IntTriple&, bool> {
//  public:
//   bool operator()(const IntTriple& a,const IntTriple& b) {
//     return( (a.i1 < b.i1) ||
//             ((a.i1 == b.i1)&&(a.i2 < b.i2)) ||
//             ((a.i1 == b.i1)&&(a.i2 == b.i2)&&(a.i3 < b.i3)) );
//   }
// };

class IntTripleCompare123 {
 public:
  inline bool operator()(const IntTriple& a,const IntTriple& b) {
    return( (a.i1 < b.i1) ||
            ((a.i1 == b.i1)&&(a.i2 < b.i2)) ||
            ((a.i1 == b.i1)&&(a.i2 == b.i2)&&(a.i3 < b.i3)) );
  }
};

typedef struct {
  int i1,i2,i3,i4;
} IntQuad;

class IntQuadCompare1 : public std::binary_function<IntQuad&, IntQuad&, bool> {
 public:
  bool operator()(const IntQuad& a,const IntQuad& b) {
    return( a.i1 < b.i1 );
  }
};

class IntQuadCompare4 : public std::binary_function<IntQuad&, IntQuad&, bool> {
 public:
  bool operator()(const IntQuad& a,const IntQuad& b) {
    return( a.i4 < b.i4 );
  }
};

class IntQuadCompare12 : public std::binary_function<IntQuad&, IntQuad&, bool> {
 public:
  bool operator()(const IntQuad& a,const IntQuad& b) {
    return( (a.i1 < b.i1) || ((a.i1 == b.i1)&&(a.i2 < b.i2)) );
  }
};

class IntQuadCompare123 : public std::binary_function<IntQuad&, IntQuad&, bool> {
 public:
  bool operator()(const IntQuad& a,const IntQuad& b) {
    return( (a.i1 < b.i1) ||
            ((a.i1 == b.i1)&&(a.i2 < b.i2)) ||
            ((a.i1 == b.i1)&&(a.i2 == b.i2)&&(a.i3 < b.i3)) );
  }
};

class IntQuadCompare1234 : public std::binary_function<IntQuad&, IntQuad&, bool> {
 public: 
  bool operator()(const IntQuad& a,const IntQuad& b) {
    return(  (a.i1  < b.i1) || 
            ((a.i1 == b.i1)&&(a.i2  < b.i2)) ||
            ((a.i1 == b.i1)&&(a.i2 == b.i2)&&(a.i3  < b.i3)) ||
            ((a.i1 == b.i1)&&(a.i2 == b.i2)&&(a.i3 == b.i3)&&(a.i4 < b.i4)) );
  }
};

inline bool compareFirst2(const std::pair<int8,int>& a,const std::pair<int8,int>& b) {
  // no duplicates: interestingly, this can fail for the subsequent
  // vector sort in the FlaggedFaceWriter, even though there should not be any duplicates...
  //assert( a.first != b.first );
  //assert( a.second != b.second );
  return( a.first < b.first );
}

struct PlanarCheck {
  bool bPlanar;

  PlanarCheck(int * noofa_v, double x_no[][3], const int nof_f, const double tol) {
    const int ino0 = noofa_v[nof_f];
    const int ino1 = noofa_v[nof_f+1];
    const int ino2 = noofa_v[nof_f+2];
    const int ino3 = noofa_v[nof_f+3];
    const double v1[3] = {
      x_no[ino1][0] - x_no[ino0][0],
      x_no[ino1][1] - x_no[ino0][1],
      x_no[ino1][2] - x_no[ino0][2] };
    const double v2[3] = {
      x_no[ino2][0] - x_no[ino0][0],
      x_no[ino2][1] - x_no[ino0][1],
      x_no[ino2][2] - x_no[ino0][2] };

    const double normal[3] = {
      v1[1]*v2[2] - v1[2]*v2[1],
      v1[2]*v2[0] - v1[0]*v2[2],
      v1[0]*v2[1] - v1[1]*v2[0]};

    // check the distance of the 4th node off the plane of the other 3...
    const double delta =
      (x_no[ino3][0] - x_no[ino0][0])*normal[0] +
      (x_no[ino3][1] - x_no[ino0][1])*normal[1] +
      (x_no[ino3][2] - x_no[ino0][2])*normal[2];

    bPlanar = fabs(delta) < tol*pow(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2], 0.75);
  }
};

class ColoredTriMesh {
public :
  int ntr ;
  int ntr_local ;
  int nno ;
  int (*nootr)[3] ;
  double (*x_no)[3];
  unsigned short* fa_color ;

  int * ugp_fa_flag ;
  int * recv_counts ;
  int * recv_disp ;

  ColoredTriMesh() : ntr(0), ntr_local(0), nno(0), nootr(NULL), x_no(NULL), fa_color(NULL),
     ugp_fa_flag(NULL), recv_counts(NULL), recv_disp(NULL) {}
  ~ColoredTriMesh() {
    DELETE(nootr);
    DELETE(x_no );
    DELETE(fa_color);
    DELETE(ugp_fa_flag);
    DELETE(recv_counts);
    DELETE(recv_disp) ;
  }
};

typedef struct {
    int v0;
    int v1;
    int v2;
} TriConnectivity;

typedef struct {
    double x;
    double y;
    double z;
    double v;
} TriData;

typedef std::pair<double,double> CTI_double2 ;

// used by mpi's DOUBLE_INT...

typedef struct { double this_double; int this_int; } DoubleInt;


class SharedMemDescriptor { 
public: 
  size_t nbytes;   // n*sizeof(T) stored directly.. 
  std::string filename;
  int fd_shm;

  SharedMemDescriptor() : nbytes(0), filename(""), fd_shm(-1) {}
  SharedMemDescriptor(const size_t _n, const std::string& _filename, const int _fd_shm) : 
    nbytes(_n), filename(_filename), fd_shm(_fd_shm) {}

};

inline std::string stringify_shm_name(const int i) { 
  std::stringstream ss; 
  ss << "/cti_shm";
  ss << i;
  return ss.str();
}

inline std::string stringify_shm_name(const int i,const int j) { 
  std::stringstream ss; 
  ss << "/cti_shm";
  ss << i;
  ss << "_";
  ss << j;
  return ss.str();
}

#endif
