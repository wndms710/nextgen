#ifndef _PLANE_IMAGE_MAP_
#define _PLANE_IMAGE_MAP_

class CvPlaneMap {
public:
  float xp[3];
  float np[3];
private:
  vector<int8> cvopl_vec; // unique, ordered, cvs
  vector<double> x_vv_vec;
  vector<double> r_vv_vec;
  bool b_mapSet;
public:

  CvPlaneMap() {
    FOR_I3 xp[i] = 0.0;
    FOR_I3 np[i] = 0.0;
    b_mapSet = false;
  }

  CvPlaneMap(const float _xp[3], const float _np[3]) {
    init(_xp,_np);
  }

  void init(const float _xp[3], const float _np[3]){
    float np_mag = MAG(_np);
    if (b_mapSet){ //check if new _xp, _np matches the existing
      FOR_I3 {
        if (xp[i] != _xp[i]){
          b_mapSet=false;
          break;
        }

        if (np[i] != _np[i]/np_mag){
          b_mapSet=false;
          break;
        }
      }
      if (!b_mapSet){
        cout << "Resetting CvPlaneMap data" << endl;
        cvopl_vec.clear();
        x_vv_vec.clear();
        r_vv_vec.clear();
        FOR_I3 xp[i] = _xp[i];
        FOR_I3 np[i] = _np[i]/np_mag;
      }
    }
    else{
      FOR_I3 xp[i] = _xp[i];
      FOR_I3 np[i] = _np[i]/np_mag;
    }
  }

  bool checkMapFlag(){
    return b_mapSet;
  }

  void setMapFlag(){
    b_mapSet = true;
    //report size of vector data
    COUT1("CvPlaneMap Data Report: ncv " << cvopl_vec.size() << " [" << cvopl_vec.size()*40/(1024.0*1024.0) << " MBs]");
  }

  int size(){
    return cvopl_vec.size();
  }

  void push_back(const int8 &icv, const double x_vv[3], const double &r_vv){
    cvopl_vec.push_back(icv);
    FOR_I3 x_vv_vec.push_back(x_vv[i]);
    r_vv_vec.push_back(r_vv);
  }

  int8 getIcvGlobal(const int8 &icv){
    assert(icv<int8(cvopl_vec.size()));
    return cvopl_vec[icv];
  }

  void getXvv(double x_vv[3],const int &icv){
    assert(icv<int8(x_vv_vec.size())/3);
    FOR_I3 x_vv[i] = x_vv_vec[icv*3+i];
  }

  void getXcs(double xc_s[3],const int &icv){
    double x_vv[3];
    getXvv(x_vv,icv);
    const double dist =
       (x_vv[0]-xp[0])*np[0] +
       (x_vv[1]-xp[1])*np[1] +
       (x_vv[2]-xp[2])*np[2];
    FOR_I3 xc_s[i] = x_vv[i] - dist*np[i];
  }

  double getRs2(const int &icv){
    assert(icv<int(r_vv_vec.size()));
    double r_vv = r_vv_vec[icv];
    double x_vv[3];
    getXvv(x_vv,icv);
    const double dist =
       (x_vv[0]-xp[0])*np[0] +
       (x_vv[1]-xp[1])*np[1] +
       (x_vv[2]-xp[2])*np[2];
    return r_vv*r_vv - dist*dist;
  }

  ~CvPlaneMap() {

  }


};


#endif
