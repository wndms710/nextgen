#ifndef _PRCOMM_HPP_
#define _PRCOMM_HPP_

class CvPrcomm {
public:
  int rank;
  vector<int> packVec;
  int pack_offset,pack_size; // should be same as packVec.size()

  // we unpack directly into the ghost data, so no vector<int> unpackVec;
  int unpack_offset,unpack_size;

  // the packVec is divided into pack "ranges" associated with
  // different periodic bcs...
  //int packRangeSize;
  //int * packRangeIndex;
  //int * packRangeBits;

  // for variable-sized pack/unpack...
  int npack_v,nunpack_v;
  int pack_offset_v,unpack_offset_v;
  
  // note: translations and rotations are stored in the transformVec that spans
  // all data in the pack vec, even if both are NULL...

  class Transform {
  public:
    int bits;
    bool has_R;
    bool has_t;
    double R[9];
    double t[3];
    int start,end; // range is [start:end), i.e. does NOT include end
    Transform(const bool has_R_,const double * R_,const bool has_t_,const double * t_,const int bits_,const int start_) {
      has_R = has_R_;
      if (has_R_) {
	for (int i = 0; i < 9; ++i) R[i] = R_[i];
      }
      has_t = has_t_;
      if (has_t_) {
	FOR_I3 t[i] = t_[i];
      }
      bits = bits_;
      start = start_;
    }
  };
  vector<Transform> transformVec;

  CvPrcomm() {
    //packRangeSize  = 0;
    //packRangeIndex = NULL;
    //packRangeBits  = NULL;
  }

  virtual ~CvPrcomm() {
    //DELETE(packRangeIndex);
    //DELETE(packRangeBits);
  }

  int getRank() const { return rank; }
  
  void copy(const CvPrcomm& other) {
    // this copy routine is used when the CvPrcomm is duplicated in the StaticSolverNS namespace.
    // it should be eliminated when we go to an entirely namespace-based static solver...
    rank = other.rank;
    packVec = other.packVec; 
    assert(packVec.size() == other.packVec.size());
    pack_offset = other.pack_offset;
    pack_size = other.pack_size; // should be same as packVec.size()
    unpack_offset = other.unpack_offset;
    unpack_size = other.unpack_size;
    //packRangeSize = other.packRangeSize;
    //assert(other.packRangeIndex == NULL);
    //assert(other.packRangeBits == NULL);
    npack_v = other.npack_v;
    nunpack_v = other.nunpack_v;
    pack_offset_v = other.pack_offset_v;
    unpack_offset_v = other.unpack_offset_v;
    transformVec = other.transformVec;
    assert(transformVec.size() == other.transformVec.size());
  }

};

// the general Prcomm has both packVec and unpackVec...

class Prcomm : public CvPrcomm {
public:
  vector<int> unpackVec;
};

class MpiRequestStuff {
public:
  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  int nsend,nrecv; // used to allocate above requests when prcomm is asymmetric (i.e. some zeros)
  int * pack_buf_int;
  int * unpack_buf_int;
  int8 * pack_buf_int8;
  int8 * unpack_buf_int8;
  double * pack_buf_double;
  double * unpack_buf_double;
  void * pack_buf_v;
  void * unpack_buf_v;
  MpiRequestStuff() {
    sendRequestArray  = NULL;
    recvRequestArray  = NULL;
    pack_buf_int      = NULL;
    unpack_buf_int    = NULL;
    pack_buf_int8     = NULL;
    unpack_buf_int8   = NULL;
    pack_buf_double   = NULL;
    unpack_buf_double = NULL;
    pack_buf_v        = NULL;
    unpack_buf_v      = NULL;
  }
  ~MpiRequestStuff() {
    // just check that we are all cleaned up...
    assert(sendRequestArray  == NULL);
    assert(recvRequestArray  == NULL);
    assert(pack_buf_int      == NULL);
    assert(pack_buf_int8     == NULL);
    assert(pack_buf_double   == NULL);
    assert(pack_buf_v        == NULL);
    assert(unpack_buf_int    == NULL);
    assert(unpack_buf_int8   == NULL);
    assert(unpack_buf_double == NULL);
    assert(unpack_buf_v      == NULL);
  }
};

#endif
