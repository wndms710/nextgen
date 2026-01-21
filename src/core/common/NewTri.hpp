#ifndef _NEW_TRI_HPP_
#define _NEW_TRI_HPP_

// temporary container to hold new tri information
class NewTri {
public:
  int spost[3];
  int znost,szost;
  NewTri() {
    spost[0] = -1;
    spost[1] = -1;
    spost[2] = -1;
    znost = -1;
    szost = -1;
  }
  NewTri(const NewTri& other) {
    FOR_I3 spost[i] = other.spost[i];
    znost = other.znost;
    szost = other.szost;
  }
  NewTri(const int sp0,const int sp1,const int sp2) {
    spost[0] = sp0;
    spost[1] = sp1;
    spost[2] = sp2;
    znost = -1;
    szost = -1;
  }
  NewTri(const int sp0,const int sp1,const int sp2,const int zone,const int subzone) {
    spost[0] = sp0;
    spost[1] = sp1;
    spost[2] = sp2;
    znost = zone;
    szost = subzone;
  }
  NewTri(const int _spost[3],const int zone,const int subzone) {
    FOR_I3 spost[i] = _spost[i];
    znost = zone;
    szost = subzone;
  }
};

#endif
