#ifndef _PLANE_DATA_HPP_
#define _PLANE_DATA_HPP_

template <class T>
class PlaneData {
public:
  T center[3];
  T normal[3];
  PlaneData() {}
  PlaneData(const T center[3],const T normal[3]) {
    FOR_I3 {
      this->center[i] = center[i];
      this->normal[i] = normal[i];
    }
  }
};

#endif
