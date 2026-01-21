#ifndef IMAGE_METADATA_HPP
#define IMAGE_METADATA_HPP

// ======================================================================
// ImageMetadata.hpp
//
// Class defining additional data stored in the header of PNG images
// written by the Cascade solvers
//
//   -Metadata that is suitable for direct display is saved as uncompressed
//    text metadata
//   -Metadata utilized primarily by the WebUI is saved in a JsonObject
//    as compressed text metadata
//
// History:
//  author: David Philips
// ======================================================================

// ======================================================================
// Copyright (c) 2015 Cascade Technologies Inc.
// All Rights Reserved. www.cascadetechnologies.com
//
// This source is proprietary software and cannot be used, redistributed,
// or modified except under the terms of Cascade's Software License
// Agreement.
//
// THIS SOURCE CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY
// OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
// ======================================================================
using namespace std;

// 65531 - first "special" zone index - need to clarify exactly what this means

class ImageMetadata {

public:
  enum {
    ZONEID_MIN    = 0,
    ZONEID_MAX    = 32766,
    MESH_DATA     = 32767, // includes planes,isos,surfaces and particles W/O data
    OEDGEID_MIN   = 32768, // also used as offset for surface boundary data
    OEDGEID_MAX   = 49151,
    MERDGEID_MIN  = 49152,
    MERDGEID_MAX  = 65531,
    ISO_DATA      = 65532,
    PARTICLE_DATA = 65533,
    VOLUME_DATA   = 65534,
    BACKGROUND    = 65535
  };

private:
  map<string,string> textMap;  //User viewable data

public:

  double transformMat[16];  //column based, 4x4 transformation matrix
  vector<string> zoneNames;
  vector<uint2> zoneIds;
  vector<uint2> subZoneMin;
  vector<uint2> subZoneMax;

  vector<float> planeXp;
  vector<float> planeNp;

  ImageMetadata() {
    init();
  }

  void init() {
    //Init textMap with empty strings for all valid tEXt fields,
    //will be useful for reverse process of initializing an
    //ImageMetadata object from a read image.
    textMap["WRITE_IMAGE"]            = ""; //command used to generate the image
    textMap["RESTART_HASH_ID"]        = ""; //identifies restart file
    textMap["TIME"]                   = ""; //simulation time
    textMap["VAR_ON_SURFACE"]         = ""; //simulation surface variable u, p, etc.
    textMap["VAR_ON_PARTICLE"]        = ""; //simulation particle variable u, T, etc.
    textMap["VAR_ON_ISO"]             = ""; //simulation particle variable u, T, etc.
    textMap["VAR"]                    = ""; //simulation variable u, p, etc.
    textMap["COLORMAP"]               = ""; //name of colormap: Grayscale, Hot Metal, etc.
    textMap["COLORMAP_SURFACE"]       = ""; //name of colormap: Grayscale, Hot Metal, etc.
    textMap["COLORMAP_PARTICLE"]      = ""; //name of colormap: Grayscale, Hot Metal, etc.
    textMap["COLORMAP_ISO"]           = ""; //name of colormap: Grayscale, Hot Metal, etc.
    textMap["UNIT_SYSTEM_NAME"]       = ""; //name of the unit system used by the computation
    textMap["RANGE_ON_SURFACE_MIN"]   = ""; //surface var contour map minimum given in refUnitsSystem
    textMap["RANGE_ON_SURFACE_MAX"]   = ""; //surface var contour map maximum given in refUnitsSystem
    textMap["RANGE_ON_PARTICLE_MIN"]  = ""; //particle var contour map minimum given in refUnitsSystem
    textMap["RANGE_ON_PARTICLE_MAX"]  = ""; //particle var contour map maximum given in refUnitsSystem
    textMap["RANGE_ON_ISO_MIN"]       = ""; //iso var contour map minimum given in refUnitsSystem
    textMap["RANGE_ON_ISO_MAX"]       = ""; //iso var contour map maximum given in refUnitsSystem
    textMap["RANGE_MIN"]              = ""; //contour map minimum given in refUnitsSystem
    textMap["RANGE_MAX"]              = ""; //contour map maximum given in refUnitsSystem
    textMap["LENGTH_SCALE"]           = ""; //pixel dimension in simulation coordinates
    textMap["RGB_MODE"]               = ""; //DATA, NORMALS
    textMap["CAM_DEPTH"]              = ""; //depth offset in pixels of the camera from closest surface
    textMap["LIGHT_DATA"]             = "FALSE"; //tell app to use normal when coloring plane data
    textMap["HASH_MLES"]              = ""; //Hash from mesh or surface
    textMap["HASH_SLES"]              = ""; //Hash from solution data used to init

    //initialize transformMat to identity matrix
    transformMat[0] = 1.0;
    transformMat[1] = 0.0;
    transformMat[2] = 0.0;
    transformMat[3] = 0.0;
    transformMat[4] = 0.0;
    transformMat[5] = 1.0;
    transformMat[6] = 0.0;
    transformMat[7] = 0.0;
    transformMat[8] = 0.0;
    transformMat[9] = 0.0;
    transformMat[10] = 1.0;
    transformMat[11] = 0.0;
    transformMat[12] = 0.0;
    transformMat[13] = 0.0;
    transformMat[14] = 0.0;
    transformMat[15] = 1.0;
  }

  void resetAll(){
    textMap.clear();
    zoneNames.clear();
    zoneIds.clear();
    subZoneMin.clear();
    subZoneMax.clear();
    planeXp.clear();
    planeNp.clear();
    init();
  }

  //Set methods for User viewable data (tEXt chunks)
  void setWriteImageParam(const string &_paramStr){
    textMap["WRITE_IMAGE"] = _paramStr;
  }

  void setLightData(const bool light = false) {
    if (light)
      textMap["LIGHT_DATA"] = "TRUE";
    else
      textMap["LIGHT_DATA"] = "FALSE";
  }

  void setParentRestartHash(const string &_parentRestartHash){
    textMap["RESTART_HASH_ID"]  = _parentRestartHash;
  }
  void setTime(const double &_time){
    char strbuf[64]={'\0'};
    sprintf(strbuf,"%.15g",_time);
    textMap["TIME"]             = strbuf;
  }
  void setVarOnSurfaceId(const string &_varOnSurfaceId){
    textMap["VAR_ON_SURFACE"]   = _varOnSurfaceId;
  }
  void setVarOnParticleId(const string &_varOnParticleId){
    textMap["VAR_ON_PARTICLE"]   = _varOnParticleId;
  }
  void setVarOnIsoId(const string &_varOnIsoId){
    textMap["VAR_ON_ISO"]   = _varOnIsoId;
  }
  void setVarId(const string &_varId){
    textMap["VAR"]              = _varId;
  }
  void setColorMapName(const string &_colorMapName){
    textMap["COLORMAP"]         = _colorMapName;
  }
  void setColorMapSurfaceName(const string &_colorMapName){
    textMap["COLORMAP_SURFACE"]         = _colorMapName;
  }
  void setColorMapParticlesName(const string &_colorMapName){
    textMap["COLORMAP_PARTICLE"]         = _colorMapName;
  }
  void setColorMapIsoName(const string &_colorMapName){
    textMap["COLORMAP_ISO"]         = _colorMapName;
  }
  void setRefUnitSystem(const string &_refUnitSystem){
    textMap["UNIT_SYSTEM_NAME"] = _refUnitSystem;
  }
  void setRangeMin(const double &_rangeMin){
    char strbuf[64]={'\0'};
    sprintf(strbuf,"%.15g",_rangeMin);
    textMap["RANGE_MIN"] = strbuf;
  }
  void setRangeMax(const double &_rangeMax){
    char strbuf[64]={'\0'};
    sprintf(strbuf,"%.15g",_rangeMax);
    textMap["RANGE_MAX"] = strbuf;
  }
  void setRangeOnSurfaceMin(const double &_rangeOnSurfaceMin){
    char strbuf[64]={'\0'};
    sprintf(strbuf,"%.15g",_rangeOnSurfaceMin);
    textMap["RANGE_ON_SURFACE_MIN"] = strbuf;
  }
  void setRangeOnSurfaceMax(const double &_rangeOnSurfaceMax){
    char strbuf[64]={'\0'};
    sprintf(strbuf,"%.15g",_rangeOnSurfaceMax);
    textMap["RANGE_ON_SURFACE_MAX"] = strbuf;
  }
  void setRangeOnParticleMin(const double &_rangeOnParticleMin){
    char strbuf[64]={'\0'};
    sprintf(strbuf,"%.15g",_rangeOnParticleMin);
    textMap["RANGE_ON_PARTICLE_MIN"] = strbuf;
  }
  void setRangeOnParticleMax(const double &_rangeOnParticleMax){
    char strbuf[64]={'\0'};
    sprintf(strbuf,"%.15g",_rangeOnParticleMax);
    textMap["RANGE_ON_PARTICLE_MAX"] = strbuf;
  }
  void setRangeOnIsoMin(const double &_rangeOnIsoMin){
    char strbuf[64]={'\0'};
    sprintf(strbuf,"%.15g",_rangeOnIsoMin);
    textMap["RANGE_ON_ISO_MIN"] = strbuf;
  }
  void setRangeOnIsoMax(const double &_rangeOnIsoMax){
    char strbuf[64]={'\0'};
    sprintf(strbuf,"%.15g",_rangeOnIsoMax);
    textMap["RANGE_ON_ISO_MAX"] = strbuf;
  }

  void setLengthScale(const double &_lengthScale){
    char strbuf[64]={'\0'};
    sprintf(strbuf,"%.15g",_lengthScale);
    textMap["LENGTH_SCALE"]     = strbuf;
  }
  void setRgbMode(const int &_rgbMode){
    string rgbModeStr;
    switch (_rgbMode){
      case 0:
        rgbModeStr = "DATA";
        break;
      case 1:
        rgbModeStr = "NORMALS";
        break;
      default:
        rgbModeStr = "";
        break;
    }
    textMap["RGB_MODE"]         = rgbModeStr;
  }
  void setCamDepth(const float &_camDepth){
    char strbuf[64]={'\0'};
    sprintf(strbuf,"%.15g",_camDepth);
    textMap["CAM_DEPTH"]        = strbuf;
  }
  void setHashMles(const string &_mlesHash){
    textMap["HASH_MLES"]         = _mlesHash;
  }
  void setHashSles(const string &_slesHash){
    textMap["HASH_SLES"]         = _slesHash;
  }

  //set methods for webui (json) data (zTXt chunk)
  void addZoneNameId(const string &name, const uint2 &zoneid) {
    zoneNames.push_back(name);
    zoneIds.push_back(zoneid);
  }

  void clearZoneNameId() {
    zoneNames.clear();
    zoneIds.clear();
  }

  void addPlanePointAndNormal(const float xp[3], const float np[3]){
    for (int i=0; i<3; ++i) {
      planeXp.push_back(xp[i]);
      planeNp.push_back(np[i]);
    }
  }

  void clearPlanePointAndNormal() {
    planeXp.clear();
    planeNp.clear();
  }

  //PngImage Call 1 - pass uncompressed text metadata,
  //                  only set values will be written;
  void buildTextMetadataVector(vector<string> &textKeyValue) {
    textKeyValue.clear();
    for (map<string, string>::iterator it = textMap.begin(); it!=textMap.end(); it++){
      if (it->second != ""){
        textKeyValue.push_back(it->first);
        textKeyValue.push_back(it->second);
      }
    }
  }

  //PngImage Call 2 - put transformMat and zone info into a
  //                  json formatted string.  JSONMetadata will
  //                  always be written with at least the
  //                  transformMat (default is identity matrix)
  void buildZTextMetadataVector(vector<string> &ztextKeyValue) {
    ztextKeyValue.clear();
    ztextKeyValue.push_back("JSONMetadata");
    ztextKeyValue.push_back(getJsonMetadataStr());
    //cout << getJsonMetadataStr() << endl;
  }


  //INIT ImageMetadata object after reading png image metadata
  void setTextFromImage(const char * key, const char * value, const bool verbose) {
    string keyStr = key;
    string valStr = value;
    if ( textMap.find(keyStr) != textMap.end() ){
      textMap[keyStr] = valStr;
    }
    else if (keyStr == "JSONMetadata") {

      int tmk = valStr.find("transformMat");
      int tm0 = valStr.find("[",tmk);  //start of transformMat data
      int tm1 = valStr.find("]",tm0);  //end of transformMat data

      string tm_str = valStr.substr(tm0+1,tm1-tm0-1);

      int prev, pos;
      vector<string> tm_vec;
      prev = tm_str.find_first_not_of(" ,\n\t",0);
      while ((pos= tm_str.find_first_of(",",prev)) != int(std::string::npos))
      {
        //pos set to comma location, prev set to start of number
        if (pos > prev){
          tm_vec.push_back(tm_str.substr(prev, pos-prev));
        }
        prev=tm_str.find_first_not_of(" ,\n\t",pos);
      }
      if (prev < int(tm_str.length())) {  //no comma after last entry in matrix
        pos=tm_str.find_first_of(" ,\n\t",prev);
        tm_vec.push_back(tm_str.substr(prev, pos-prev));
      }
      assert(tm_vec.size()==16); //16 elements in transformMat
      for (int i =0; i<16; ++i){
        //cout << "transformMat: " << tm_vec[i] << endl;
        transformMat[i] = atof(tm_vec[i].c_str());
      }

      //parse zoneIds and zoneNames
      int znk = valStr.find("zoneNames");
      if (znk != int(std::string::npos)){
        zoneNames.clear();
        zoneIds.clear();
        int zn0 = valStr.find("[",znk);  //start of zoneNames data
        int zn1 = valStr.find("]",zn0);  //end of zoneNames data
        string zn_str = valStr.substr(zn0+1,zn1-zn0-1);

        prev = zn_str.find_first_not_of(" ,\n\t",0);
        while ((pos= zn_str.find_first_of(",",prev)) != int(std::string::npos))
        {
          //pos set to comma location, prev set to start of number
          if (pos > prev){
            zoneNames.push_back(zn_str.substr(prev+1, pos-prev-2));
          }
          prev=zn_str.find_first_not_of(" ,\n\t",pos);
        }
        if (prev < int(zn_str.length())) {  //no comma after last entry in matrix
          pos=zn_str.find_first_of(" ,\n\t",prev);
          zoneNames.push_back(zn_str.substr(prev+1, pos-prev-2));
        }

        znk = valStr.find("zoneIds");
        zn0 = valStr.find("[",znk);  //start of zoneIds data
        zn1 = valStr.find("]",zn0);  //end of zoneIds data
        zn_str = valStr.substr(zn0+1,zn1-zn0-1);

        vector<string> zn_vec;
        prev = zn_str.find_first_not_of(" ,\n\t",0);
        while ((pos= zn_str.find_first_of(",",prev)) != int(std::string::npos))
        {
          //pos set to comma location, prev set to start of number
          if (pos > prev){
            zn_vec.push_back(zn_str.substr(prev, pos-prev));
          }
          prev=zn_str.find_first_not_of(" ,\n\t",pos);
        }
        if (prev < int(zn_str.length())) {  //no comma after last entry in matrix
          pos=zn_str.find_first_of(" ,\n\t",prev);
          zn_vec.push_back(zn_str.substr(prev, pos-prev));
        }
        for (int i =0,i_max=zn_vec.size(); i<i_max; ++i){
          zoneIds.push_back((uint2)atof(zn_vec[i].c_str()));
        }
      }

      //parse planeXp and planeNp
      int plk = valStr.find("planeXp");
      if (plk != int(std::string::npos)){
        planeXp.clear();
        planeNp.clear();

        vector<string> pl_vec;

        int pl0 = valStr.find("[",plk);  //start of planeXp data
        int pl1 = valStr.find("\n]",pl0);  //end of planeXp data
        string pl_str = valStr.substr(pl0+1,pl1-pl0-1);

        prev = pl_str.find_first_not_of(" ,\n\t[]",0);
        int count=0;
        string search_char = ",";
        while ((pos= pl_str.find_first_of(search_char,prev)) != int(std::string::npos))
        {
          //pos set to comma location, prev set to start of number
          if (pos > prev){
            pl_vec.push_back(pl_str.substr(prev, pos-prev));
          }
          prev=pl_str.find_first_not_of(" ,\n\t[]",pos);

          ++count;
          if (count%3==2)
            search_char = "]";
          else
            search_char = ",";
        }
        for (int i =0,i_max=pl_vec.size(); i<i_max; ++i){
          planeXp.push_back(atof(pl_vec[i].c_str()));
        }

        pl_vec.clear();

        plk = valStr.find("planeNp");
        pl0 = valStr.find("[",plk);  //start of planeNp data
        pl1 = valStr.find("\n]",pl0);  //end of planeNp data
        pl_str = valStr.substr(pl0+1,pl1-pl0-1);

        prev = pl_str.find_first_not_of(" ,\n\t[]",0);
        count=0;
        search_char = ",";
        while ((pos= pl_str.find_first_of(search_char,prev)) != int(std::string::npos))
        {
          //pos set to comma location, prev set to start of number
          if (pos > prev){
            pl_vec.push_back(pl_str.substr(prev, pos-prev));
          }
          prev=pl_str.find_first_not_of(" ,\n\t[]",pos);

          ++count;
          if (count%3==2)
            search_char = "]";
          else
            search_char = ",";
        }
        for (int i =0,i_max=pl_vec.size(); i<i_max; ++i){
          planeNp.push_back(atof(pl_vec[i].c_str()));
        }

      }


    }
    else if (verbose)
      cout << "Ignoring tEXt metadata: " << key << " " << value << endl;
  }

  bool hasText(const string &keyStr){
     if (textMap.find(keyStr) == textMap.end())
       return false;

     return textMap[keyStr] != "";

  }

  void dumpInfo(ostream& ofile){
    ofile << "*****teXt metadata*****" << endl;
    for (map<string,string>::iterator it=textMap.begin(); it!=textMap.end(); ++it)
      ofile << it->first << " = " << it->second << endl;
    ofile << "*****ztXt metadata*****" << endl;
    ofile << getJsonMetadataStr() << endl;
  }

  //getters for access after reading
  string getWriteImageParam(){
    return(textMap["WRITE_IMAGE"]);
  }
  double getRangeMin(){
    return(atof(textMap["RANGE_MIN"].c_str()));
  }
  double getRangeMax(){
    return(atof(textMap["RANGE_MAX"].c_str()));
  }
  double getSurfRangeMin(){
    return(atof(textMap["RANGE_ON_SURFACE_MIN"].c_str()));
  }
  double getSurfRangeMax(){
    return(atof(textMap["RANGE_ON_SURFACE_MAX"].c_str()));
  }

  double getPartRangeMin(){
    return(atof(textMap["RANGE_ON_PARTICLE_MIN"].c_str()));
  }
  double getPartRangeMax(){
    return(atof(textMap["RANGE_ON_PARTICLE_MAX"].c_str()));
  }

  double getIsoRangeMin(){
    return(atof(textMap["RANGE_ON_ISO_MIN"].c_str()));
  }
  double getIsoRangeMax(){
    return(atof(textMap["RANGE_ON_ISO_MAX"].c_str()));
  }
  string getVarId(){
    return textMap["VAR"];
  }
  string getSurfVarId(){
    return textMap["VAR_ON_SURFACE"];
  }
  string getPartVarId(){
    return textMap["VAR_ON_PARTICLE"];
  }
  string getIsoVarId(){
    return textMap["VAR_ON_ISO"];
  }
  double getTime(){
    return(atof(textMap["TIME"].c_str()));
  }
  double getLengthScale(){
    return(atof(textMap["LENGTH_SCALE"].c_str()));
  }
  double getCamDepth(){
    return(atof(textMap["CAM_DEPTH"].c_str()));
  }
  string getColorMapName(){
    return(textMap["COLORMAP"]);
  }
  string getSurfColorMapName(){
    return(textMap["COLORMAP_SURFACE"]);
  }
  string getPartColorMapName(){
    return(textMap["COLORMAP_PARTICLE"]);
  }
  string getIsoColorMapName(){
    return(textMap["COLORMAP_ISO"]);
  }
  int getRgbMode(){
    if (textMap["RGB_MODE"]=="NORMALS")
      return 1;
    else
      return 0;
  }

  private:

  string getJsonMetadataStr(){

    // write the JSON file...
    stringstream ss;
    ss << "{\n\"transformMat\":[\n";
    ss << setprecision(12) <<transformMat[0] <<",\n";
    ss << setprecision(12) <<transformMat[1] <<",\n";
    ss << setprecision(12) <<transformMat[2] <<",\n";
    ss << setprecision(12) <<transformMat[3] <<",\n";
    ss << setprecision(12) <<transformMat[4] <<",\n";
    ss << setprecision(12) <<transformMat[5] <<",\n";
    ss << setprecision(12) <<transformMat[6] <<",\n";
    ss << setprecision(12) <<transformMat[7] <<",\n";
    ss << setprecision(12) <<transformMat[8] <<",\n";
    ss << setprecision(12) <<transformMat[9] <<",\n";
    ss << setprecision(12) <<transformMat[10]<<",\n";
    ss << setprecision(12) <<transformMat[11]<<",\n";
    ss << setprecision(12) <<transformMat[12]<<",\n";
    ss << setprecision(12) <<transformMat[13]<<",\n";
    ss << setprecision(12) <<transformMat[14]<<",\n";
    ss << setprecision(12) <<transformMat[15]<<"\n";

    if (zoneIds.size()>0 && zoneNames.size()==zoneIds.size()){
      ss <<"],\n";

      ss << "\"zoneIds\":[\n";
      for (int izi = 0,limit=zoneIds.size()-1; izi<limit; ++izi){
        int zoneIdInt = zoneIds[izi];
        ss << zoneIdInt << ",\n";
      }
      ss << zoneIds[zoneIds.size()-1] << "\n";
      ss <<"],\n";

      ss << "\"zoneNames\":[\n";
      for (int izi = 0,limit=zoneNames.size()-1; izi<limit; ++izi){
        ss << "\"" << zoneNames[izi] << "\",\n";
      }
      ss << "\"" << zoneNames[zoneNames.size()-1] << "\"\n";
    }
    if (planeXp.size()>0 && planeXp.size()==planeNp.size()){
      ss <<"],\n";

      ss << "\"planeXp\":[\n";
      for (int ixp = 0,limit=planeXp.size()-1; ixp<limit; ++ixp){
        if (ixp%3==0)
          ss << "[";

        ss << planeXp[ixp];

        if (ixp%3==2)
          ss << "],\n";
        else
          ss << ",";
      }
      ss << planeXp[planeXp.size()-1] << "]\n";
      ss <<"],\n";

      ss << "\"planeNp\":[\n";
      for (int inp = 0,limit=planeNp.size()-1; inp<limit; ++inp){
        if (inp%3==0)
          ss << "[";

        ss << planeNp[inp];

        if (inp%3==2)
          ss << "],\n";
        else
          ss << ",";
      }
      ss << planeNp[planeNp.size()-1] << "]\n";
    }
    ss <<"]\n";

    ss <<"}\n";

    //cout << "***********************" << endl;
    //cout << ss.str().length();
    //cout << ss.str();
    //cout << "***********************" << endl;

    return ss.str();

  }




};


#endif
