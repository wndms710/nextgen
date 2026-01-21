#ifndef _FLOAT_VERTEX_HPP_
#define _FLOAT_VERTEX_HPP_

class FloatVertex {
public:
  float x,y,z;
  FloatVertex() {}
  FloatVertex(const float xyz[3]) {
    x = xyz[0];
    y = xyz[1];
    z = xyz[2];
  }
  bool operator<(const FloatVertex& other) const { 
    return (x < other.x) || ((x == other.x) && ((y < other.y) || ((y == other.y)&&(z < other.z))));
  }
  bool operator!=(const FloatVertex& other) const { 
    return (x != other.x) || (y != other.y) || (z != other.z);
  }
  bool operator==(const FloatVertex& other) const { 
    return (x == other.x) && (y == other.y) && (z == other.z);
  }
};

#endif
