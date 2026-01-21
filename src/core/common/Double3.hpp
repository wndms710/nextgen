#ifndef _DOUBLE3_HPP_
#define _DOUBLE3_HPP_

class Double3 {
public:
  double data[3];
  Double3() {
    FOR_I3 data[i] = 0.0;
  }    
  Double3(const double data[3]) {
    FOR_I3 this->data[i] = data[i];
  }
  Double3(const double data[3],const double factor) {
    FOR_I3 this->data[i] = data[i]*factor;
  }
  inline int operator[] (int i) const {
    assert((i >= 0)&&(i < 3));
    return data[i];
  }
};

#endif
