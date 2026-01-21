#ifndef _ST_ZONE_HPP_
#define _ST_ZONE_HPP_

class StZone {
private:
  string name;
  uint2 periodic_bits;
public:
  int flag; // multi-purpose flag
  StZone() {
    this->name = "";
    this->periodic_bits = 0;
    this->flag = 0;
  }
  StZone(const string& name) {
    this->name = name;
    this->periodic_bits = 0;
    this->flag = 0;
  }
  StZone(const StZone& other) {
    name = other.name;
    periodic_bits = other.periodic_bits;
    flag = other.flag;
  }
  string getName() const {
    return name;
  }
  void setName(const string& name) {
    this->name = name;
  }
  void setName(const char * name) {
    this->name = name;
  }
  uint2 getPeriodicBits() const {
    return periodic_bits;
  }
  void setPeriodicBits(const uint2 periodic_bits_) {
    periodic_bits = periodic_bits_;
  }
  bool isPeriodic() const {
    return periodic_bits != 0;
  }
  bool isBoundary() const {
    return periodic_bits == 0;
  }
};

#endif
