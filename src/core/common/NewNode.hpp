#ifndef _NEW_NODE_HPP_
#define _NEW_NODE_HPP_

// temporary container to hold new node information
class NewNode {
public:
  double xsp[3];
  int isp;  // keep track if this is on top of an existing node
  double data_d;  // double data
  NewNode() {
    FOR_I3 xsp[i] = 0.0;
    isp = -1;
    data_d = 0.0;
  }
  NewNode(const NewNode& other) {
    FOR_I3 xsp[i] = other.xsp[i];
    isp = other.isp;
    data_d = other.data_d;
  }
  NewNode(const double x,const double y,const double z) {
    xsp[0] = x;
    xsp[1] = y;
    xsp[2] = z;
    isp = -1;
    data_d = 0.0;
  }
  NewNode(const double _xsp[3]) {
    FOR_I3 xsp[i] = _xsp[i];
    isp = -1;
    data_d = 0.0;
  }
  void setIsp(const int isp_) {
    isp = isp_;
  }
  int getIsp() const {
    return isp;
  }
};

#endif
