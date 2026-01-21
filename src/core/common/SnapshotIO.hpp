#ifndef SNAPSHOT_HPP
#define SNAPSHOT_HPP

class Snapshot {
public: 

  string name;  // filename prefix
  int interval; 
  vector<pair<string,CtiRegister::CtiData*> > var_vec;

  Snapshot() : name(""), interval(-1) {
    var_vec.clear();
  }
};
#endif
