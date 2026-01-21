#include "PeriodicData.hpp"

namespace PeriodicData {

  vector<PeriodicTransform> periodicTransformVec;

  bool isPeriodic() {
    return !periodicTransformVec.empty();
  }

  bool checkPeriodicBitPair(const int i) {
    return i < int(periodicTransformVec.size());
  }

  bool getPeriodicR(double R[9],const int bits) {

    // set I...
    R[0] = 1.0; R[1] = 0.0; R[2] = 0.0;
    R[3] = 0.0; R[4] = 1.0; R[5] = 0.0;
    R[6] = 0.0; R[7] = 0.0; R[8] = 1.0;

    // loop on bit pairs...
    bool has_R = false;
    FOR_I3 {
      if (bits & (1<<(2*i))) {
	// even bit pair set...
	assert(!(bits & (1<<(2*i+1)))); // only one pair set
	assert(int(periodicTransformVec.size()) > i); // make sure transform vec supports this bit pair
	double this_R[9];
	if (periodicTransformVec[i].getR(this_R)) {
	  has_R = true;
	  double Rtmp[9] = { 0.0, 0.0, 0.0,
			     0.0, 0.0, 0.0,
			     0.0, 0.0, 0.0 };
	  FOR_M3 {
	    FOR_J3 {
	      FOR_K3 {
		Rtmp[3*m+j] += R[3*m+k]*this_R[3*k+j];
	      }
	    }
	  }
	  for (int ii =0; ii < 9; ++ii)
	    R[ii] = Rtmp[ii];
	}
      }
      else if (bits & (1<<(2*i+1))) {
	// odd bit pair set...
	assert(int(periodicTransformVec.size()) > i); // make sure transform vec supports this bit pair
	double this_R[9];
	if (periodicTransformVec[i].getInvR(this_R)) {
	  has_R = true;
	  double Rtmp[9] = { 0.0, 0.0, 0.0,
			     0.0, 0.0, 0.0,
			     0.0, 0.0, 0.0 };
	  FOR_M3 {
	    FOR_J3 {
	      FOR_K3 {
		Rtmp[3*m+j] += R[3*m+k]*this_R[3*k+j];
	      }
	    }
	  }
	  for (int ii =0; ii < 9; ++ii)
	    R[ii] = Rtmp[ii];
	}
      }
    }
    return has_R;
  }

  bool getPeriodicT(double t[3],const int bits) {

    // zero t...
    t[0] = 0.0; t[1] = 0.0; t[2] = 0.0;

    // loop on bit pairs...
    bool has_t = false;
    FOR_I3 {
      if (bits & (1<<(2*i))) {
	// even bit pair set...
	assert(!(bits & (1<<(2*i+1)))); // only one pair set
	assert(int(periodicTransformVec.size()) > i); // make sure transform vec supports this bit pair
	double this_t[3];
	if (periodicTransformVec[i].getT(this_t)) {
	  has_t = true;
	  FOR_J3 t[j] += this_t[j];
	}
      }
      else if (bits & (1<<(2*i+1))) {
	// odd bit pair set...
	assert(int(periodicTransformVec.size()) > i); // make sure transform vec supports this bit pair
	double this_t[3];
	if (periodicTransformVec[i].getInvT(this_t)) {
	  has_t = true;
	  FOR_J3 t[j] += this_t[j];
	}
      }
    }
    return has_t;
  }

  void periodicTranslate(double (*xp_t)[3],const int n,const int bits) {
    for (int bit_pair = 0; bit_pair < 3; ++bit_pair) {
      if (bits & (1<<(2*bit_pair))) {
	assert((bits & (1<<(2*bit_pair+1)))==0);
	assert(int(periodicTransformVec.size()) > bit_pair);
	periodicTransformVec[bit_pair].translate(xp_t,n);
      }
      else if (bits & (1<<(2*bit_pair+1))) {
	assert(int(periodicTransformVec.size()) > bit_pair);
	periodicTransformVec[bit_pair].inv_translate(xp_t,n);
      }
    }
    // these ones are the doubles...
    for (int bit_pair = 3; bit_pair < 6; ++bit_pair) {
      if (bits & (1<<(2*bit_pair))) {
	assert((bits & (1<<(2*bit_pair+1)))==0);
	assert(int(periodicTransformVec.size()) > bit_pair-3);
	periodicTransformVec[bit_pair-3].translate(xp_t,n);
	periodicTransformVec[bit_pair-3].translate(xp_t,n);
      }
      else if (bits & (1<<(2*bit_pair+1))) {
	assert(int(periodicTransformVec.size()) > bit_pair-3);
	periodicTransformVec[bit_pair-3].inv_translate(xp_t,n);
	periodicTransformVec[bit_pair-3].inv_translate(xp_t,n);
      }
    }

  }

  void periodicTranslate(double *xp_t,const int n,const int bits) {
    periodicTranslate((double(*)[3])xp_t,n,bits);
  }

  void periodicRotate(double (*xp_t)[3],const int n,const int bits) {
    for (int bit_pair = 0; bit_pair < 3; ++bit_pair) {
      if (bits & (1<<(2*bit_pair))) {
	assert((bits & (1<<(2*bit_pair+1)))==0);
	assert(int(periodicTransformVec.size()) > bit_pair);
	periodicTransformVec[bit_pair].rotate(xp_t,n);
      }
      else if (bits & (1<<(2*bit_pair+1))) {
	assert(int(periodicTransformVec.size()) > bit_pair);
	periodicTransformVec[bit_pair].inv_rotate(xp_t,n);
      }
    }
    for (int bit_pair = 3; bit_pair < 6; ++bit_pair) {
      if (bits & (1<<(2*bit_pair))) {
	assert((bits & (1<<(2*bit_pair+1)))==0);
	assert(int(periodicTransformVec.size()) > bit_pair-3);
	periodicTransformVec[bit_pair-3].rotate(xp_t,n);
	periodicTransformVec[bit_pair-3].rotate(xp_t,n);
      }
      else if (bits & (1<<(2*bit_pair+1))) {
	assert(int(periodicTransformVec.size()) > bit_pair-3);
	periodicTransformVec[bit_pair-3].inv_rotate(xp_t,n);
	periodicTransformVec[bit_pair-3].inv_rotate(xp_t,n);
      }
    }
  }

  void periodicRotate(double *xp_t,const int n,const int bits) {
    periodicRotate((double(*)[3])xp_t,n,bits);
  }

  void scalePeriodicTransformVec(double s[3]) {

    for(vector<PeriodicTransform>::iterator it = periodicTransformVec.begin();
        it != periodicTransformVec.end(); ++it) {
      it->scale(s);
    }

  }

  void buildPeriodicRtVec(vector<PeriodicRt>& periodicRtVec) {
    switch (periodicTransformVec.size()) {
    case 0:
      periodicRtVec.push_back(PeriodicRt());
      break;
    case 1:
      periodicRtVec.push_back(PeriodicRt());
      // HACK: remove copies for single periodicity to check nasa37 (not now)...
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<0)); // ie 10
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<0)); // ie 10
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<1)); // ie 01
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<1)); // ie 01
      break;
    case 2:
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<0)); // ie 1000
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<0)); // ie 1000
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<1)); // ie 0100
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<1)); // ie 0100
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<2)); // ie 0010
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<2)); // ie 0010
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<3)); // ie 0001
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<3)); // ie 0001
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<0)|(1<<2)); // ie 1010
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<0)|(1<<2)); // ie 1010
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<1)|(1<<2)); // ie 0110
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<1)|(1<<2)); // ie 0110
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<0)|(1<<3)); // ie 1001
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<0)|(1<<3)); // ie 1001
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<1)|(1<<3)); // ie 0101
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<1)|(1<<3)); // ie 0101
      break;
    case 3:
      // 0...
      periodicRtVec.push_back(PeriodicRt());
      // 1...
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<0)); // ie 100000
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<0)); // ie 100000
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<1)); // ie 010000
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<1)); // ie 010000
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<2)); // ie 001000
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<2)); // ie 001000
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<3)); // ie 000100
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<3)); // ie 000100
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<4)); // ie 000010
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<4)); // ie 000010
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<5)); // ie 000001
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<5)); // ie 000001
      // 2...
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<0)|(1<<2)); // ie 101000
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<0)|(1<<2)); // ie 101000
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<1)|(1<<2)); // ie 011000
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<1)|(1<<2)); // ie 011000
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<0)|(1<<3)); // ie 100100
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<0)|(1<<3)); // ie 100100
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<1)|(1<<3)); // ie 010100
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<1)|(1<<3)); // ie 010100
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<0)|(1<<4)); // ie 100010
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<0)|(1<<4)); // ie 100010
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<1)|(1<<4)); // ie 010010
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<1)|(1<<4)); // ie 010010
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<0)|(1<<5)); // ie 100001
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<0)|(1<<5)); // ie 100001
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<1)|(1<<5)); // ie 010001
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<1)|(1<<5)); // ie 010001
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<2)|(1<<4)); // ie 001010
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<2)|(1<<4)); // ie 001010
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<3)|(1<<4)); // ie 000110
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<3)|(1<<4)); // ie 000110
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<2)|(1<<5)); // ie 001001
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<2)|(1<<5)); // ie 001001
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<3)|(1<<5)); // ie 000101
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<3)|(1<<5)); // ie 000101
      // 3...
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<0)|(1<<2)|(1<<4)); 
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<0)|(1<<2)|(1<<4));
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<1)|(1<<2)|(1<<4)); 
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<1)|(1<<2)|(1<<4));
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<0)|(1<<3)|(1<<4)); 
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<0)|(1<<3)|(1<<4));
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<1)|(1<<3)|(1<<4)); 
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<1)|(1<<3)|(1<<4));
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<0)|(1<<2)|(1<<5)); 
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<0)|(1<<2)|(1<<5));
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<1)|(1<<2)|(1<<5)); 
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<1)|(1<<2)|(1<<5));
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<0)|(1<<3)|(1<<5)); 
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<0)|(1<<3)|(1<<5));
      periodicRtVec.push_back(PeriodicRt());
      periodicRtVec.back().b_R = getPeriodicR(periodicRtVec.back().R,(1<<1)|(1<<3)|(1<<5)); 
      periodicRtVec.back().b_t = getPeriodicT(periodicRtVec.back().t,(1<<1)|(1<<3)|(1<<5));
      break;
    default:
      assert(0);
    }
  }

  void buildBitVec(vector<int>& bitsVec) {

    assert(bitsVec.empty());
    
    switch (periodicTransformVec.size()) {
    case 0:
      bitsVec.push_back(0);
      break;
    case 1:
      bitsVec.push_back(0);
      bitsVec.push_back((1<<0));
      bitsVec.push_back((1<<1));
      break;
    case 2:
      bitsVec.push_back(0);
      bitsVec.push_back((1<<0));
      bitsVec.push_back((1<<1));
      bitsVec.push_back((1<<2));
      bitsVec.push_back((1<<3));
      bitsVec.push_back((1<<2)|(1<<0));
      bitsVec.push_back((1<<3)|(1<<0));
      bitsVec.push_back((1<<2)|(1<<1));
      bitsVec.push_back((1<<3)|(1<<1));
      break;
    case 3:
      bitsVec.push_back(0);
      bitsVec.push_back((1<<0));
      bitsVec.push_back((1<<1));
      bitsVec.push_back((1<<2));
      bitsVec.push_back((1<<3));
      bitsVec.push_back((1<<4));
      bitsVec.push_back((1<<5));
      bitsVec.push_back((1<<0)|(1<<2));
      bitsVec.push_back((1<<0)|(1<<3));
      bitsVec.push_back((1<<0)|(1<<4));
      bitsVec.push_back((1<<0)|(1<<5));
      bitsVec.push_back((1<<1)|(1<<2));
      bitsVec.push_back((1<<1)|(1<<3));
      bitsVec.push_back((1<<1)|(1<<4));
      bitsVec.push_back((1<<1)|(1<<5));
      bitsVec.push_back((1<<2)|(1<<4));
      bitsVec.push_back((1<<2)|(1<<5));
      bitsVec.push_back((1<<3)|(1<<4));
      bitsVec.push_back((1<<3)|(1<<5));
      bitsVec.push_back((1<<0)|(1<<2)|(1<<4));
      bitsVec.push_back((1<<0)|(1<<2)|(1<<5));
      bitsVec.push_back((1<<0)|(1<<3)|(1<<4));
      bitsVec.push_back((1<<0)|(1<<3)|(1<<5));
      bitsVec.push_back((1<<1)|(1<<2)|(1<<4));
      bitsVec.push_back((1<<1)|(1<<2)|(1<<5));
      bitsVec.push_back((1<<1)|(1<<3)|(1<<4));
      bitsVec.push_back((1<<1)|(1<<3)|(1<<5));
      break;
    default:
      assert(0);
    }
    
  }

} // namespace PeriodicData
