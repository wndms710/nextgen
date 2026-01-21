#ifndef _SIMPLE_POINT_HPP_
#define _SIMPLE_POINT_HPP_

class SimplePoint {
public:
  double x[3];
  SimplePoint() {
    this->x[0] = 0.0;
    this->x[1] = 0.0;
    this->x[2] = 0.0;
  }
  SimplePoint(const SimplePoint& other) {
    this->x[0] = other.x[0];
    this->x[1] = other.x[1];
    this->x[2] = other.x[2];
  }
  SimplePoint(const double x[3]) {
    this->x[0] = x[0];
    this->x[1] = x[1];
    this->x[2] = x[2];
  }
  SimplePoint(const double x,const double y,const double z) {
    this->x[0] = x;
    this->x[1] = y;
    this->x[2] = z;
  }
  bool operator<(const SimplePoint& other) const {
    return (x[0] < other.x[0]) || ((x[0] == other.x[0]) && ((x[1] < other.x[1]) || ((x[1] == other.x[1])&&(x[2] < other.x[2]))));
  }
  bool operator!=(const SimplePoint& other) const {
    return (x[0] != other.x[0]) || (x[1] != other.x[1]) || (x[2] != other.x[2]);
  }
  bool operator==(const SimplePoint& other) const {
    return (x[0] == other.x[0]) && (x[1] == other.x[1]) && (x[2] == other.x[2]);
  }
};

#endif
