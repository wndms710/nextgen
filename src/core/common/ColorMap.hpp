//
//  ColorMap.hpp
//

#ifndef _ColorMap_hpp
#define _ColorMap_hpp

#include "Common.hpp"
#include <fstream>

class ColorMap_Base {

protected:

  bool b_dataFromRgb;

public:
  ColorMap_Base() : b_dataFromRgb(false) { }

  virtual ~ColorMap_Base() { }

  virtual string getName(const bool _b_flip)=0;

  //override if colormap goes by multiple names...
  virtual bool isAlias(const string &nameWithoutInverted) const {
    return false;
  };
  virtual int getImageType()=0;

  virtual inline void calcColorMapRGB(float _rgb[3], const float &_phi01, const bool _b_flip) const {
    //You should not be calling this funtion if your colormap is properly defined...
    assert(0);
  }

  virtual inline float calcPhiFromRgb(const unsigned char _rgb[3], const bool _b_flip) const {
    //please check b_dataFromRgb before calling this function...
    assert(0);
  }


  inline void calcColorMapRGB(unsigned char _rgb[3], const float &phi01, const bool _b_flip) const {
    float _rgb_f[3];
    calcColorMapRGB(_rgb_f, phi01, _b_flip);

    for(int i=0; i<3; ++i) _rgb[i] = (unsigned char)_rgb_f[i];
  }

  inline bool get_b_dataFromRgb() {return b_dataFromRgb;}

};

class ColorMap_HOT_METAL: public ColorMap_Base {
  static const int mapSize = 5;
  static const float local_r[mapSize], local_g[mapSize], local_b[mapSize];

  float colormap_r[8];
  float colormap_g[8];
  float colormap_b[8];

public:
  ColorMap_HOT_METAL() {
    for (int i=0;i<mapSize;i++){
      colormap_r[i] = local_r[i];
      colormap_g[i] = local_g[i];
      colormap_b[i] = local_b[i];
    }
    b_dataFromRgb = true;
  }

  string getName(const bool _b_flip) {
    if (_b_flip)
      return "INVERTED_HOT_METAL";
    else
      return "HOT_METAL";
  }

  int getImageType() { return 0; }

  inline void calcColorMapRGB(float _rgb[3], const float &_phi01, const bool _b_flip) const {
    float phi01;
    if (_b_flip) phi01 = (1.0f-_phi01);
    else        phi01 = _phi01;

    int bin = (int) floor(phi01 / 0.250001);
    assert(bin < mapSize-1 && bin>=0);
    float w1 = (phi01 - 0.250001f*float(bin)) / 0.250001f;
    float w2 = 1.0f - w1;
    _rgb[0] = (float)( (w2*colormap_r[bin] + w1*colormap_r[bin+1])  );
    _rgb[1] = (float)( (w2*colormap_g[bin] + w1*colormap_g[bin+1])  );
    _rgb[2] = (float)( (w2*colormap_b[bin] + w1*colormap_b[bin+1])  );
  }

  inline float calcPhiFromRgb(const unsigned char _rgb[3], const bool _b_flip) const {
    int bin;
    float w1;
    uint8_t r = _rgb[0];
    uint8_t g = _rgb[1];
    uint8_t b = _rgb[2];

    if (b>0) {
      bin=3;
      w1 = b/255.0f;
    }
    else if (g>0) {
      bin=2;
      w1 = g/255.0f;
    }
    else if (r>127) {
      bin=1;
      w1 = (r-127.0f)/(255.0f-127.0f);
    }
    else {
      bin=0;
      w1 = (r-2.0f)/(127.0f-2.0f);
    }

    float phi = min(max((w1+bin)*0.250001f,0.0f),1.0f);

    if (_b_flip) return (1.0f-phi);
    else        return phi;
  }

};

class ColorMap_BLUE_RED_RAINBOW: public ColorMap_Base {
  static const int mapSize = 5;
  static const float local_r[mapSize], local_g[mapSize], local_b[mapSize];

  float colormap_r[8];
  float colormap_g[8];
  float colormap_b[8];

public:
  ColorMap_BLUE_RED_RAINBOW() {
    for (int i=0;i<mapSize;i++){
      colormap_r[i] = local_r[i]*((float)255.);
      colormap_g[i] = local_g[i]*((float)255.);
      colormap_b[i] = local_b[i]*((float)255.);
    }
    b_dataFromRgb = true;
  }

  string getName(const bool _b_flip) {
    if (_b_flip)
      return "INVERTED_BLUE_RED_RAINBOW";
    else
      return "BLUE_RED_RAINBOW";
  }

  int getImageType() { return 0; }

  inline void calcColorMapRGB(float _rgb[3], const float &_phi01, const bool _b_flip) const {
    float phi01;
    if (_b_flip) phi01 = (1.0f-_phi01);
    else        phi01 = _phi01;

    int bin = (int) floor(phi01 / 0.250001);
    assert(bin < mapSize-1 && bin>=0);
    float w1 = (phi01 - 0.250001f*float(bin)) / 0.250001f;
    float w2 = 1.0f - w1;
    _rgb[0] = (float) ( (w2*colormap_r[bin] + w1*colormap_r[bin+1])  );
    _rgb[1] = (float) ( (w2*colormap_g[bin] + w1*colormap_g[bin+1])  );
    _rgb[2] = (float) ( (w2*colormap_b[bin] + w1*colormap_b[bin+1])  );
  }

  inline float calcPhiFromRgb(const unsigned char _rgb[3], const bool _b_flip) const {
    int bin;
    float w1;
    uint8_t r = _rgb[0];
    uint8_t g = _rgb[1];
    uint8_t b = _rgb[2];
    if (r >= 255) {
      bin=3;
      w1 = 1.0f-(g/255.0f);
    }

    else if (r>0) {
      bin=2;
      w1 = r/255.0f;
    }
    else if (g>=255) {
      bin=1;
      w1 = 1.0f-(b/255.0f);
    }
    else {
      bin=0 ;
      w1 = (g/255.0f);
    }

    float phi = min(max((w1+bin)*0.250001f,0.0f),1.0f);

    if (_b_flip) return (1.0f-phi);
    else        return phi;
  }


};

class ColorMap_COOL_WARM: public ColorMap_Base {
  static const int mapSize = 5;
  static const float local_r[mapSize], local_g[mapSize], local_b[mapSize];

  float colormap_r[8];
  float colormap_g[8];
  float colormap_b[8];

public:
  ColorMap_COOL_WARM() {
    for (int i=0;i<mapSize;i++){
      colormap_r[i] = local_r[i]*((float)255.);
      colormap_g[i] = local_g[i]*((float)255.);
      colormap_b[i] = local_b[i]*((float)255.);
    }
    b_dataFromRgb = true;
  }

  string getName(const bool _b_flip) {
    if (_b_flip)
      return "INVERTED_COOL_WARM";
    else
      return "COOL_WARM";
  }


  int getImageType() { return 0; }

  inline void calcColorMapRGB(float _rgb[3], const float &_phi01, const bool _b_flip) const {
    float phi01;
    if (_b_flip) phi01 = (1.0f-_phi01);
    else        phi01 = _phi01;

    int bin = (int) floor(phi01 / 0.50001);
    assert(bin < mapSize-1 && bin>=0);
    float w1 = (phi01 - 0.50001f*float(bin)) / 0.50001f;
    float w2 = 1.0f - w1;
    _rgb[0] = (float) ( (w2*colormap_r[bin] + w1*colormap_r[bin+1])  );
    _rgb[1] = (float) ( (w2*colormap_g[bin] + w1*colormap_g[bin+1])  );
    _rgb[2] = (float) ( (w2*colormap_b[bin] + w1*colormap_b[bin+1])  );
  }

  inline float calcPhiFromRgb(const unsigned char _rgb[3], const bool _b_flip) const {
    int bin;
    float w1;
    uint8_t r = _rgb[0];
    //uint8_t g = _rgb[1];
    uint8_t b = _rgb[2];

    if (r<191) {
      bin=0;
      w1 = r/191.25f;
    }
    else {
      bin=1;
      w1 = 1.0-b/191.25f;
    }

    float phi = min(max((w1+bin)*0.50001f,0.0f),1.0f);

    if (_b_flip) return (1.0f-phi);
    else        return phi;
  }

};

class ColorMap_RAINBOW_DESATURATED: public ColorMap_Base {
  static const int mapSize = 8;
  static const float local_r[mapSize], local_g[mapSize], local_b[mapSize];

  float colormap_r[8];
  float colormap_g[8];
  float colormap_b[8];

public:
  ColorMap_RAINBOW_DESATURATED() {
    for (int i=0;i<mapSize;i++){
      colormap_r[i] = local_r[i]*((float)255.);
      colormap_g[i] = local_g[i]*((float)255.);
      colormap_b[i] = local_b[i]*((float)255.);
    }
  }

  string getName(const bool _b_flip) {
    if (_b_flip)
      return "INVERTED_RAINBOW_DESATURATED";
    else
      return "RAINBOW_DESATURATED";
  }

  int getImageType() { return 0; }

  inline void calcColorMapRGB(float _rgb[3], const float &_phi01, const bool _b_flip) const {
    float phi01;
    if (_b_flip) phi01 = (1.0f-_phi01);
    else        phi01 = _phi01;

    int bin = (int) floor(phi01 / 0.142858);
    assert(bin < mapSize-1 && bin>=0);
    float w1 = (phi01 - 0.142858f*float(bin)) / 0.142858f;
    float w2 = 1.0f - w1;
    _rgb[0] = (float) ( (w2*colormap_r[bin] + w1*colormap_r[bin+1])  );
    _rgb[1] = (float) ( (w2*colormap_g[bin] + w1*colormap_g[bin+1])  );
    _rgb[2] = (float) ( (w2*colormap_b[bin] + w1*colormap_b[bin+1])  );
  }
};

class ColorMap_GRAYSCALE: public ColorMap_Base {

public:

  ColorMap_GRAYSCALE() {
    b_dataFromRgb = true;
  }

  string getName(const bool _b_flip) {
    if (_b_flip)
      return "INVERTED_GRAYSCALE_RGB";
    else
      return "GRAYSCALE_RGB";
  }

  bool isAlias(const string &nameWithoutInverted) const {
    if (nameWithoutInverted=="GRAYSCALE"||nameWithoutInverted=="GREYSCALE")
      return true;
    else
      return false;
  }

  int getImageType() { return 0; }
  inline void calcColorMapRGB(float _rgb[3], const float &_phi01, const bool _b_flip) const {
    float phi01;
    if (_b_flip) phi01 = (1.0f-_phi01);
    else        phi01 = _phi01;

    _rgb[0] = (float)( phi01 * float(255) + 0.5f );
    _rgb[1] = (float)( phi01 * float(255) + 0.5f );
    _rgb[2] = (float)( phi01 * float(255) + 0.5f );
  }

  inline float calcPhiFromRgb(const unsigned char _rgb[3], const bool _b_flip) const {
    int _rgb_i = (int) _rgb[0];
    float phi = _rgb_i/255.0f;
    if (_b_flip) return (1.0f-phi);
    else        return phi;
  }

};

class ColorMap_USERDEFINED: public ColorMap_Base {

private:

  string name;
  int mapSize;
  float * phi_ref;
  uint8_t (*colormap_rgb)[3];

public:

  ColorMap_USERDEFINED() : name("CUSTOM"), mapSize(0), phi_ref(NULL), colormap_rgb(NULL) {
    b_dataFromRgb = false;
  }

  bool readColorMapFile(const string &filename){

    if (filename.substr(filename.length()-4) == ".lut") {
      cout << " > processing LUT file" << endl;
      ifstream infile(filename.c_str());
      if (infile.is_open()) {
        size_t start = filename.find_last_of("/\\");
        if (start == string::npos) start = 0;
        else ++start;
        const size_t end = filename.length() - 4 - start;
        name = filename.substr(start,end);

        mapSize = 256;
        phi_ref = new float[mapSize];
        colormap_rgb = new uint8_t[mapSize][3];

        int irow=0;
        string line;
        while (std::getline(infile, line) && irow<mapSize )
        {
          istringstream iss(line);
          float r_f, g_f, b_f;
          if (!(iss >> r_f >> g_f >> b_f)) {
             cerr << "Unable to process colormap line " << irow << " in " << filename << endl;
             return false;
          }
          phi_ref[irow] = float(irow)/float(mapSize);
          colormap_rgb[irow][0] = r_f;
          colormap_rgb[irow][1] = g_f;
          colormap_rgb[irow][2] = b_f;
          cout << "[" << irow << "] " << phi_ref[irow] << " " << r_f << " " << g_f << " " << b_f << endl;
          ++irow;
        }
      }
      return true;
    }


    ifstream infile(filename.c_str());
    if (infile.is_open()) {
      getline(infile,name);  //first line contains colormap name
      string line;
      getline(infile,line);  //second line contains the color scaling 1 or 255
      int colorScale = atoi(line.c_str());
      if (colorScale!=1&&colorScale!=255){
        cerr << "Invalid colormap color scaling " << colorScale << " for " << filename << ". Accepted values are 1 and 255" << endl;
        return false;
      }
      getline(infile,line);  //third line contains number of control points
      mapSize = atoi(line.c_str());
      if (mapSize <= 0){
         cerr << "Invalid colormap control point size in " << filename << " " << mapSize << endl;
         return false;
      }
      phi_ref = new float[mapSize];
      colormap_rgb = new uint8_t[mapSize][3];

      int irow = 0;
      //cout << "USERDEFINED: " << name << " " << mapSize << endl;
      //cout << "USERDEFINED: " << phi_ref[i] << " " << colormap_rgb[i][0] << " "  << colormap_rgb[i][1] << " " << colormap_rgb[i][2] << endl;
      while (std::getline(infile, line) && irow<mapSize )
      {
        istringstream iss(line);
        float r_f, g_f, b_f;
        if (!(iss >> phi_ref[irow] >> r_f >> g_f >> b_f)) {
           cerr << "Unable to process colormap line " << irow+3 << " in " << filename << endl;
           return false;
        }
        //cout << phi_ref[irow] << " " << r_int << " "  << g_int << " "  << b_int << endl;
        if (colorScale==1){
          colormap_rgb[irow][0] = r_f*255.0f + 0.5f;
          colormap_rgb[irow][1] = g_f*255.0f + 0.5f;
          colormap_rgb[irow][2] = b_f*255.0f + 0.5f;
        }
        else{
          colormap_rgb[irow][0] = r_f;
          colormap_rgb[irow][1] = g_f;
          colormap_rgb[irow][2] = b_f;
        }
        ++irow;
      }
      if (irow<mapSize){
         cerr << "Unable able to read full colormap record; only " << irow << " row(s) of " << mapSize << " available in " << filename << endl;
         return false;
      }
      cout << "Read user-defined colormap " << name << " with " << mapSize << " control points " << endl;
      cout << "  first row " << setw(16) << phi_ref[0]         << setw(6) << (int) colormap_rgb[0][0]         << setw(6) << (int) colormap_rgb[0][1]         << setw(6) << (int) colormap_rgb[0][2] << endl;
      cout << "   last row " << setw(16) << phi_ref[mapSize-1] << setw(6) << (int) colormap_rgb[mapSize-1][0] << setw(6) << (int) colormap_rgb[mapSize-1][1] << setw(6) << (int) colormap_rgb[mapSize-1][2] << endl;

      return true;
    }
    cerr << "Unable to open colormap file: " << filename << endl;

    return false;
  }

  ~ColorMap_USERDEFINED() {
    cout << "~ColorMap_USERDEFINED() " << name << endl; //TODO remove

    if (phi_ref){
      delete[] phi_ref;
      phi_ref=NULL;
    }
    if (colormap_rgb){
      delete[] colormap_rgb;
      colormap_rgb=NULL;
    }
  }

  string getName(const bool _b_flip) {
    if (_b_flip)
      return "INVERTED_"+name;
    else
      return name;
  }

  int getImageType() { return 0; }

  //output _rgb should be pixel scaled (not 0-1)
  inline void calcColorMapRGB(float _rgb[3], const float &_phi01, const bool _b_flip) const {
    assert(phi_ref&&colormap_rgb);

    float phi01;
    if (_b_flip) phi01 = (1.0f-_phi01);
    else        phi01 = _phi01;

    //check for _phi01 outside control points
    //(i.e. user did not specify
    // color values for phi=0 and phi=1)
    if (phi01<phi_ref[0]){
      FOR_I3 _rgb[i] = colormap_rgb[0][i];
    }
    else if (phi01>=phi_ref[mapSize-1]){
      FOR_I3 _rgb[i] = colormap_rgb[mapSize-1][i];
    }
    else{
      //linear interpolate color between control points
      for (int bin = 0; bin<mapSize-1; ++bin){
        if (phi01<phi_ref[bin+1]){
          float w0 = min(max((phi01-phi_ref[bin])/(phi_ref[bin+1]-phi_ref[bin]),0.0f),1.0f);
          FOR_I3{
            int bin_size = colormap_rgb[bin+1][i] - colormap_rgb[bin][i];
            //apply 0.5 shift to map w0<0.5/bin_size to colormap_rgb[bin][i]
            //            and to map w0>(bin_size-0.5)/bin_size to colormap_rgb[bin+1][i]
            _rgb[i] = colormap_rgb[bin][i] + w0*bin_size + 0.5f; //float is trucated to a uchar in calling process
          }
          //cout << "w0 " << w0 << " PHI " << phi01 << " rgb " << _rgb[0] << " " << _rgb[1] << " " << _rgb[2] << endl;
          break;
        }
      }
    }
  }

};



class ColorMap {

private:

  string my_colormap_str;
  bool b_flip;

public:
  static map<const string,ColorMap_Base*> colormapsMap;

  static void initStandardColorMaps(){
    assert(colormapsMap.find("GRAYSCALE_RGB")==colormapsMap.end());
    assert(colormapsMap.find("HOT_METAL")==colormapsMap.end());
    assert(colormapsMap.find("BLUE_RED_RAINBOW")==colormapsMap.end());
    assert(colormapsMap.find("COOL_WARM")==colormapsMap.end());
    assert(colormapsMap.find("RAINBOW_DESATURATED")==colormapsMap.end());

    colormapsMap["GRAYSCALE_RGB"]       = new ColorMap_GRAYSCALE();
    colormapsMap["HOT_METAL"]           = new ColorMap_HOT_METAL();
    colormapsMap["BLUE_RED_RAINBOW"]    = new ColorMap_BLUE_RED_RAINBOW();
    colormapsMap["COOL_WARM"]           = new ColorMap_COOL_WARM();
    colormapsMap["RAINBOW_DESATURATED"] = new ColorMap_RAINBOW_DESATURATED();
  }

  static void addUserDefinedColorMap(const string &filename){
    if (colormapsMap.size()==0){
      initStandardColorMaps();
    }

    ColorMap_USERDEFINED * theUDColorMap = new ColorMap_USERDEFINED();

    if (theUDColorMap->readColorMapFile(filename)){
      if(colormapsMap.find(theUDColorMap->getName(false))!=colormapsMap.end()){
        cout << "Warning: User defined colormap " << theUDColorMap->getName(false) << " already exists, ignoring" << endl;
      }
      else {
        colormapsMap[theUDColorMap->getName(false)] = theUDColorMap;
      }
    }
    else{
      cout << "Warning: problem adding user defined colormap, ignoring file " << filename << endl;
    }

  }

  static bool checkColorMapName(const string &colorMapName){
    bool b_isAvailableColorMap=false;

    string nameWithoutInverted;
    if (colorMapName.length()>9 && colorMapName.find("INVERTED_") == 0) {
      nameWithoutInverted = colorMapName.substr(9);
    }
    else{
      nameWithoutInverted = colorMapName;
    }

    for (map<const string,ColorMap_Base*>::iterator it = colormapsMap.begin(); it!=colormapsMap.end(); ++it){
      if (it->first == nameWithoutInverted){
        b_isAvailableColorMap = true;
        break;
      }
      else if (it->second->isAlias(nameWithoutInverted)){
        b_isAvailableColorMap = true;
        break;
      }
    }

    return b_isAvailableColorMap;
  }

  static void getValidColorMapNames(stringstream &ss){
    if (colormapsMap.size()>0){
      map<const string,ColorMap_Base*>::iterator it = colormapsMap.begin();
      ss << it->first << ", INVERTED_" << it->first;
      if (colormapsMap.size()>1){
        int count = 0;
        for (it = colormapsMap.begin(); it!=colormapsMap.end(); ++it){
          if (count++>0)
            ss << ", " << it->first << ", INVERTED_" << it->first;
        }
      }
    }
  }

  static void clearColorMaps(){
    for (map<const string,ColorMap_Base*>::iterator it = colormapsMap.begin(); it!=colormapsMap.end(); ++it){
      if (it->second){
        delete it->second;
        it->second = NULL;
      }
    }
    colormapsMap.clear();
  }

  ColorMap() : my_colormap_str("GRAYSCALE_RGB"), b_flip(false) {
    if (colormapsMap.size()==0){
      initStandardColorMaps();
    }
  }

  ColorMap(const ColorMap& _copy) {
    set_type(_copy.getName());
  }

  ColorMap(const string& token) : my_colormap_str("GRAYSCALE_RGB"), b_flip(false){
    if (colormapsMap.size()==0){
      initStandardColorMaps();
    }
    set_type(token);
  }

  ~ColorMap() {
  }

  inline bool set_type(const string& colorMapName) {

    //if token begins with "INVERTED_", the set b_flip to true

    //check if valid colormap name, return string
    //stripped of "INVERTED_" prefix if present and
    //set flip bool
    string nameWithoutInverted;
    if (colorMapName.length()>9 && colorMapName.find("INVERTED_") == 0) {
      nameWithoutInverted = colorMapName.substr(9);
      b_flip = true;
    }
    else{
      nameWithoutInverted = colorMapName;
      b_flip = false;
    }

    bool b_isAvailableColorMap=false;
    for (map<const string,ColorMap_Base*>::iterator it = colormapsMap.begin(); it!=colormapsMap.end(); ++it){
      if (it->first == nameWithoutInverted){
        b_isAvailableColorMap = true;
        break;
      }
      else if (it->second->isAlias(nameWithoutInverted)){
        b_isAvailableColorMap = true;
        nameWithoutInverted = it->first; //re-assign alias to default name
        break;
      }
    }

    if (b_isAvailableColorMap){
      my_colormap_str = nameWithoutInverted; //if an alias was used, it will be set to default
    }
    else{
      //cout << "Warning: invalid colormap named " << colorMapName << ", defaulting to GRAYSCALE" << endl;;
      my_colormap_str = "GRAYSCALE_RGB";
      b_flip = false;
    }
    return b_isAvailableColorMap;

  }

  inline string getName() const {
    return colormapsMap[my_colormap_str]->getName(b_flip);
  }

  inline int getImageType() const {
    return colormapsMap[my_colormap_str]->getImageType();
  }

  inline void calcColorMapRGB(unsigned char _rgb[3], const float &phi01) const {
    colormapsMap[my_colormap_str]->calcColorMapRGB(_rgb, phi01, b_flip);
  }

  inline bool get_b_dataFromRgb() const {return colormapsMap[my_colormap_str]->get_b_dataFromRgb();}

  inline float calcPhiFromRgb(const unsigned char _rgb[3]) const { return colormapsMap[my_colormap_str]->calcPhiFromRgb(_rgb,b_flip);}

  bool operator==(const ColorMap &other) const {
    return (getName()==other.getName());
  }

  bool operator!=(const ColorMap &other) const {
    return (getName()!=other.getName());
  }


};


#endif
