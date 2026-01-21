#ifndef _DOUBLE_VERTEX_HPP_
#define _DOUBLE_VERTEX_HPP_

class DoubleVertex {
public:
  double x,y,z;
  DoubleVertex() {}
  DoubleVertex(const double xyz[3]) {
    x = xyz[0];
    y = xyz[1];
    z = xyz[2];
  }
  bool operator<(const DoubleVertex& other) const {
    return (x < other.x) || ((x == other.x) && ((y < other.y) || ((y == other.y)&&(z < other.z))));
  }
  bool operator!=(const DoubleVertex& other) const {
    return (x != other.x) || (y != other.y) || (z != other.z);
  }
  bool operator==(const DoubleVertex& other) const {
    return (x == other.x) && (y == other.y) && (z == other.z);
  }
};

#endif
