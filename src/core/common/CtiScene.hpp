#ifndef _CTI_SCENE_HPP_
#define _CTI_SCENE_HPP_

#include "IntFlag.hpp"
#include "Params.hpp"
#include "WriteImageData.hpp"
#include "CtiCanvas.hpp"
#include "ImageMetadata.hpp"
#include "LesImageMapper.hpp"
#include "CvPlaneMap.hpp"

#include "SimpleTri.hpp" // this is a double-based tri. Any other type uses type name, e.g. IntTri, FlaotTri...
#include "RestartHashUtilities.hpp"

#include "GeomUtils.hpp"

#include "CtiRegister.hpp"
#include "WebUI.hpp"

// needs to be data producer for serial surface les reader...
// TODO reorganize lim and cim to encapsulate this...
class CtiScene : public CtiRegister::CtiDataProducer {

private:

  // had to expose this for drawing in solver's processImage...

  WriteImageData wid; //kind of like the camera settings
  CtiCanvas * canvas;

  ImageMetadata * imd;
  float rmax;
  bool b_rmax;
  float bbox[6];
  bool b_bbox;
  float center[3];
  bool b_center;

  bool b_cached;

  vector<string> zoneNames;
  vector<bool> bf_hide;  // better memory compression than array

  CvImageMap cim;
  CvImageMap cimS;
  LesImageMapper * lim;
  string filename_volumeData; // mles or les
  string filename_volumeData_2; // sles (if available)
  string lengthscaleName;

  // added for varEvalCtiData (try to think of a better way when we restructure this)
  CvImageMap* current_cim_ptr;
  LesImageMapper * current_lim_ptr;

public:

  // you can initialize the scene from the WRITE_IMAGE param
  // directly. This builds the wid from the param.

  CtiScene(Param *param) : wid(param) {
    init();
  }

  // alternatively, if you have a wid already available,
  // use this initialization. It simply copies the pod
  // from the passed wid_ into the scene's wid...

  CtiScene(const WriteImageData& wid_) : wid(wid_) {
    init();
  }

  virtual ~CtiScene(){
    if (canvas!=NULL){
      delete canvas;
      canvas = NULL;
    }
    if (imd!=NULL){
      delete imd;
      imd = NULL;
    }
    if (lim!=NULL){
      delete lim;
      lim = NULL;
    }
  }

private:

  void init() {

    canvas = NULL;
    imd    = NULL;
    b_rmax = false;
    b_center = false;
    b_bbox = false;
    b_cached = false;

    lim = NULL;
    filename_volumeData="";
    filename_volumeData_2="";
    lengthscaleName="r_vv";

  }

public:

  bool skip(const int &step){
    return !( (step%wid.interval == 0) || (wid.interval == -1) ); // -1 is one-off
  }

  void setRmax(const float &_rmax){
    rmax = _rmax;
    b_rmax = true;
  }

  void setCenter(const float &_xc,const float &_yc,const float &_zc){
    center[0] = _xc;
    center[1] = _yc;
    center[2] = _zc;
    b_center = true;
  }

  void setBoundingBbox(const float &_xmin,const float &_xmax,const float &_ymin,const float &_ymax,const float &_zmin,const float &_zmax){
    bbox[0] = _xmin;
    bbox[1] = _xmax;
    bbox[2] = _ymin;
    bbox[3] = _ymax;
    bbox[4] = _zmin;
    bbox[5] = _zmax;
    b_bbox = true;
  }

  void addZoneName(const string &zoneName){
    zoneNames.push_back(zoneName);

    // user may have hidden zones by name instead of index
    // we want to convert to index to reuse existing implementation
    // wid has already been initialized, so simply update
    // hiddenZonesSet if this zoneName found in hidden set
    // slower approach if a long zone list, but solver won't have to explicitly tell scene when it is done adding zones
    // if (wid.b_hidezonenames) {
    //   set<string>::const_iterator it=wid.hiddenZoneNamesSet.find(zoneName);
    //   if (it != wid.hiddenZoneNamesSet.end()) {
    //     // zone found in hidden set, so add its index to hiddenZonesSet
    //     wid.b_hidezones = true;
    //     wid.hiddenZonesSet.insert(zoneNames.size()-1);
    //
    //     //TODO: potentially erase it from hiddenZoneNamesSet here to make future searches smaller... worth it?
    //   }
    // }
  }

  void enablePlanarData(const string &filename, const string &tempMeshLengthscaleName){
    filename_volumeData=filename;
    lengthscaleName = tempMeshLengthscaleName; //TODO once we're not using DELTA_VV naming anymore, remove this
  }

  void enablePlanarData(const string &filename, const string &filename_2,const string &tempMeshLengthscaleName){
    filename_volumeData=filename;
    filename_volumeData_2=filename_2;
    lengthscaleName = tempMeshLengthscaleName; //TODO once we're not using DELTA_VV naming anymore, remove this
  }

  void convertHiddenZoneNamesToIndices() {
    // should run this once all zone names have been added
    // so can populate hiddenZones vec with zone indices of named zones
    // solver must explicitly call before any functions requiring knowledge of wid.hiddenZonesSet are used

    if (wid.b_hidezonenames) {
      for (set<string>::const_iterator it=wid.hiddenZoneNamesSet.begin(); it != wid.hiddenZoneNamesSet.end(); ++it) {
        vector<string>::iterator vit = find(zoneNames.begin(), zoneNames.end(), *it);
        if (vit != zoneNames.end()) {
          // found, so insert index
          wid.b_hidezones = true;
          wid.hiddenZonesSet.insert(int(std::distance(zoneNames.begin(),vit)));
        }
      }
    }
  }

  void initCanvas(){

    // fill in any missing stuff required to set the view details...
    // the view requires
    //  > target[3]
    //  > camera[3]
    //  > up[3]
    //  > width
    //  > size

    imd = new ImageMetadata();

    if (wid.b_geom_plane) {
      assert(b_center);
      assert(b_bbox);
      // this is a plane request...
      // if it is also a plane_frac, we now have the bounding box, and can set the xp...
      if (wid.b_geom_plane_frac) {
        FOR_I3 wid.xp[i] = center[i] + wid.np[i]*(wid.geom_plane_frac-0.5)*(bbox[2*i+1]-bbox[2*i]);
        COUT1("computed PLANE_FRAC xp = " << COUT_VEC(wid.xp));
      }
      //add plane info to image metadata
      imd->addPlanePointAndNormal(wid.xp,wid.np);
    }

    if (!wid.b_target) {
      assert(b_center);
      if (wid.b_geom_plane) {
        const double dx[3] = DIFF(wid.xp,center);
        const double dp = DOT_PRODUCT(dx,wid.np);
        FOR_I3 wid.target[i] = center[i] + dx[i] - dp*wid.np[i];
      }
      else {
        FOR_I3 wid.target[i] = center[i];
      }
    }

    if (!wid.b_camera) {
      // if geom (i.e. PLANE has been set, then use this to set the camera)...
      if (wid.b_geom_plane) {
        FOR_I3 wid.camera[i] = wid.target[i] + rmax*wid.np[i];
      }
      else {
        assert(b_rmax);
        // user has not specified the camera. Put the camera
        // in the default location:
        wid.camera[0] = wid.target[0] + rmax*0.75;
        wid.camera[1] = wid.target[1] + rmax*sqrt(3.0)*0.25;
        wid.camera[2] = wid.target[2] + rmax*0.5;
      }
    }

    if (!wid.b_up) {
      float forward[3] = DIFF(wid.target,wid.camera);
      const float forward_mag = MAG(forward);
      FOR_I3 forward[i] /= forward_mag;
      // try to make e0 either +x or +y...
      float e0_approx[3];
      if (fabs(forward[0]) <= fabs(forward[1])) {
        e0_approx[0] = 1.0;
        e0_approx[1] = 0.0;
        e0_approx[2] = 0.0;
      }
      else {
        e0_approx[0] = 0.0;
        e0_approx[1] = 1.0;
        e0_approx[2] = 0.0;
      }
      wid.up[0] = CROSS_PRODUCT_0(e0_approx,forward);
      wid.up[1] = CROSS_PRODUCT_1(e0_approx,forward);
      wid.up[2] = CROSS_PRODUCT_2(e0_approx,forward);
    }

    if (!wid.b_width) {
      assert(b_rmax);
      // user has not specified the width.
      // the width is based upon a best fit of the
      // domain bounding box to the current image canvas.
      // If the image target is not the bounding box center
      // additional padding is added to ensure the full
      // geometry is included in the image (with the target
      // at the center of the image).
      assert(b_center);

      // First compute maximum bounding box extent in the
      // e0 and e1 image directions

      float e2[3] = DIFF(wid.target,wid.camera);
      const float e2_mag = MAG(e2);
      assert(e2_mag > 0.0); // camera and target cannot be coincident
      FOR_I3 e2[i] = e2[i]/(e2_mag);

      float e0[3] = CROSS_PRODUCT(e2,wid.up);
      const float e0_mag = MAG(e0); assert(e0_mag > 0.0); // problem with up
      FOR_I3 e0[i] = e0[i]/(e0_mag);

      float e1[3] = CROSS_PRODUCT(e2,e0);
      const float e1_mag = MAG(e1); assert(e1_mag > 0.0);
      FOR_I3 e1[i] = e1[i]/(e1_mag); //should already be a unit vector...

      float e0_min[3], e0_max[3], e1_min[3], e1_max[3];
      FOR_I3{
        if (e0[i] > 0){
          e0_min[i] = bbox[2*i];
          e0_max[i] = bbox[2*i+1];
        }
        else{
          e0_min[i] = bbox[2*i+1];
          e0_max[i] = bbox[2*i];
        }
        if (e1[i] > 0){
          e1_min[i] = bbox[2*i];
          e1_max[i] = bbox[2*i+1];
        }
        else{
          e1_min[i] = bbox[2*i+1];
          e1_max[i] = bbox[2*i];
        }

      }
      float e0_diff[3] = DIFF(e0_max,e0_min);
      float e1_diff[3] = DIFF(e1_max,e1_min);


      //maximum bounding box height and width projected into the image plane
      float bbox_proj[2];
      bbox_proj[0] = fabs(DOT_PRODUCT(e0,e0_diff));
      bbox_proj[1] = fabs(DOT_PRODUCT(e1,e1_diff));


      //Compute the target to center distance and project along
      //the image e0, e1 axes
      float offset[3] = DIFF(wid.target,center);
      float offset_e0 = fabs(DOT_PRODUCT(offset,e0));
      float offset_e1 = fabs(DOT_PRODUCT(offset,e1));

      //Minimum width/height to capture the bounding box plus target offset
      //from the bounding box center
      float width  = 2.f*offset_e0 + bbox_proj[0];
      float height = 2.f*offset_e1 + bbox_proj[1];

      //Use aspect ratio to determine which dimension is more restrictive
      //Add 5 percent padding
      if (wid.size[0]*height > wid.size[1]*width){
        //cout << "Height limited image" << endl;
        assert(wid.size[1]>0);
        wid.width = (height*wid.size[0]/wid.size[1])*1.05;
      }
      else{
        //cout << "Width limited image" << endl;
        wid.width = width*1.05;
      }

      COUT1("computed WIDTH = " << wid.width);

    }

    imd->setWriteImageParam(wid.paramStr);
    imd->setHashMles(RestartHashUtilities::mlesHash.getAsciiHash());
    imd->setHashSles(RestartHashUtilities::slesHash.getAsciiHash());

    if (wid.b_var_on_particle){
      imd->setVarOnParticleId(wid.var_on_particle);
    }
    else{
      imd->setVarOnParticleId("OFF");
    }

    if (wid.b_var_on_surface){
      imd->setVarOnSurfaceId(wid.var_on_surface);
    }
    else{
      imd->setVarOnSurfaceId("OFF");
    }

    if (wid.b_var_on_iso){
      imd->setVarOnIsoId(wid.var_on_iso);
    }
    else{
      imd->setVarOnIsoId("OFF");
    }

    // when there is no ISO_VAR and there is an GEOM_ISO use VAR
    if ((wid.b_geom_plane&&wid.b_var)||(wid.b_geom_iso&&wid.b_var&&!wid.b_var_on_iso)) {
      imd->setVarId(wid.var);
    }
    else{
      imd->setVarId("OFF");
    }

    for (int index=0,end=zoneNames.size(); index<end; ++index){
      imd->zoneIds.push_back(index);
      imd->zoneNames.push_back(zoneNames[index]);
    }

    // this should be enough to initialize the canvas...
    canvas = new CtiCanvas(wid.target,wid.camera,wid.up,wid.width,wid.size,imd);

    // add any blanking planes from the wid to the canvas...
    if (wid.b_blank) {
      for (int i=0, limit=wid.vBlankPlaneData.size(); i<limit; ++i) {
        canvas->addBlankPlane(wid.vBlankPlaneData[i].center,wid.vBlankPlaneData[i].normal);
        //store blank plane info in image metadata
        imd->addPlanePointAndNormal(wid.vBlankPlaneData[i].center,wid.vBlankPlaneData[i].normal);
      }
    }

  }

  // ==========================================================
  // sometimes a drawing routine needs to add the requested
  // GEOM PLANE ... as a blanking plane...
  // ==========================================================

  bool isClientRequest() const {
    // for now, use the rgb_mode to detect if this is a client request...
    return wid.rgb_mode == NORMALS;
  }

  double getWidth() const {
    return wid.width;
  }

  double getVectorDensity() const {
    return wid.vector_density;
  }

  double getVectorScale() const {
    return wid.vector_scale;
  }

  int getInterval() const {
    return wid.interval;
  }

  bool hasName() const {
    return wid.b_name;
  }

  void setWriteRank(const int write_rank) {
    assert(canvas);
    canvas->write_rank = write_rank;
  }

  int getWriteRank() const {
    return canvas->write_rank;
  }

  void cacheImage() {
    if (wid.b_range){
      canvas->setDataRange(wid.range[0],wid.range[1]);
    }
    if (wid.b_range_on_surface){
      canvas->setSurfaceDataRange(wid.range_on_surface[0],wid.range_on_surface[1]);
    }
    if (wid.b_range_on_particle){
      canvas->setParticleDataRange(wid.range_on_particle[0],wid.range_on_particle[1]);
    }
    if (wid.b_range_on_iso){
      canvas->setIsoDataRange(wid.range_on_iso[0],wid.range_on_iso[1]);
    }
    if (wid.range[2] != 0.f){
      canvas->setDataRangeBins((int)wid.range[2]);
    }
    if (wid.range_on_surface[2] != 0.f){
      canvas->setSurfaceDataRangeBins((int)wid.range_on_surface[2]);
    }
    if (wid.range_on_iso[2] != 0.f){
      canvas->setIsoDataRangeBins((int)wid.range_on_iso[2]);
    }
    canvas->cacheImage();
    b_cached = true;
  }

  bool isCached() const {
    return b_cached;
  }

  void setStep(const int step) {
    wid.step = step;
  }

  void flushImage() {
    assert(b_cached);
    char filename[128];
    sprintf(filename,"%s.%08d.png",wid.name.c_str(),wid.step);
    canvas->flushImage(filename,wid.rgb_mode);
    // cleanup...
    canvas->clearPBuf();
  }

  bool hasGeomPrimitive() {
    return wid.b_geom_primitive;
  }

  string getPrimitiveName() {
    assert(wid.b_geom_primitive);
    return wid.primitive_name;
  }

  void getPrimitiveData(double primitive_data[8]) {
    assert(wid.b_geom_primitive);
    FOR_I8 primitive_data[i] = wid.primitive_data[i];
  }

  bool hasGeomIso() const {
    return wid.b_geom_iso;
  }

  double getIsoValue(const int ii) const {
    assert(wid.b_geom_iso);
    return double(wid.iso_values[ii]);
  }

  string getIsoVar(const int ii) const {
    assert(wid.b_geom_iso);
    return wid.iso_vars[ii];
  }

  string getVolvisVar() const {
    assert(wid.b_volvis);
    return wid.volvis_var;
  }

  int getIsoCount() const {
    const int n = wid.iso_vars.size();
    assert(n == int(wid.iso_values.size()));
    return n;
  }

  bool hasVolvis() const {
    return wid.b_volvis;
  }
  bool hasVolvisSurface() const {
    return wid.b_volvis_surface;
  }
  string getVolvisType() const {
    return wid.volvis_type;
  }

  bool hasGeomPlane() const {
    return wid.b_geom_plane;
  }

  void getGeomPlaneXpAndNp(double _xp[3],double _np[3]) const {
    assert(wid.b_geom_plane);
    FOR_I3 _xp[i] = wid.xp[i];
    FOR_I3 _np[i] = wid.np[i];
  }

  //Geom blanking planes get added to beginning
  //of blanking vector, serial planar vis assumes
  //first entry is its own blanking plane
  // DP: why first? this makes the methods more complicated. i.e. insert vs push_back, first vs .back(), etc
  void addGeomPlaneAsBlankPlane() {
    assert(wid.b_geom_plane);
    assert(canvas);
    canvas->insertBlankPlane(wid.xp,wid.np);
  }

  void removeGeomPlaneAsBlankPlane() {
    assert(wid.b_geom_plane);
    assert(canvas);
    canvas->deleteFirstBlankPlane();
  }

  void setGeomPlaneAsDataPlane() {
    assert(wid.b_geom_plane);
    assert(canvas);
    canvas->setDataPlane(wid.xp,wid.np);
    canvas->setBlankDataPlane(wid.b_blank_data_plane);
  }
  void removeDataPlane() {
    assert(canvas);
    canvas->removeDataPlane();
  }

  void getBlankPlanes(vector<PlaneData<float> >& blank_planes) const {
    canvas->getBlankPlanes(blank_planes);
  }

  void deleteBlankPlanes() {
    canvas->deleteBlankPlanes();
  }

  void insertBlankPlanes(const vector<PlaneData<float> >& blank_planes) {
    canvas->insertBlankPlanes(blank_planes);
  }

  bool blankDataPlane() {
    return wid.b_blank_data_plane;
  }

  bool cellFlooded() {
    return wid.b_cell_flooded;
  }

  bool hasVar() const {
    return wid.b_var;
  }
  bool hasIsoDataVar() const {
    return wid.b_var_on_iso;
  }

  string getVar() const {
    assert(wid.b_var);
    return wid.var;
  }
  string getIsoDataVar() const {
    assert(wid.b_var_on_iso);
    return wid.var_on_iso;
  }

  bool hasSolverSurface() const {
    return wid.b_solver_surface;
  }
  string getSolverSurface() const {
    assert(wid.b_solver_surface);
    return wid.solver_surface;
  }

  bool hasHiddenSubzones() const {
    return wid.b_hidesubzones;
  }

  void getHiddenSubzones(int * sz_flag,const int nsz) {
    for (int isz = 0; isz < nsz; ++isz) sz_flag[isz] = 0;
    flagHiddenSubzones(sz_flag,nsz);
  }

  void flagHiddenSubzones(int * sz_flag,const int nsz) {
    for (set<int>::iterator iter = wid.hiddenSubzonesSet.begin(); iter != wid.hiddenSubzonesSet.end(); ++iter) {
      // -1 reserved to indicate all
      if (*iter == -1) {
        // assert(wid.hiddenZonesSet.size() == 1);
        for (int isz = 0; isz < nsz; ++isz) {
          sz_flag[isz] = 1;
        }
        return;  // no need to process others
      }
      else {
        if ((*iter >= 0)&&(*iter < nsz)) sz_flag[*iter] = 1;
        else {
          if (mpi_rank == 0) cout << "hidden subzone index: " << *iter << " is out of range; ignoring" << endl;
        }
      }
    }
  }

  bool isSubzoneHidden(const int index) const {
    if ((index != -1) && isSubzoneHidden(-1)) return true;  // hide all subzones flag

    set<int>::iterator it = wid.hiddenSubzonesSet.find(index);
    return (it != wid.hiddenSubzonesSet.end());
  }

  bool hasHiddenZones() const {
    return (wid.b_hidezones);
  }

  void getHiddenZones(int * zn_flag,const int nzn) {
    for (int izn = 0; izn < nzn; ++izn) zn_flag[izn] = 0;
    flagHiddenZones(zn_flag,nzn);
  }

  void flagHiddenZones(int * zn_flag,const int nzn) {
    for (set<int>::iterator iter = wid.hiddenZonesSet.begin(); iter != wid.hiddenZonesSet.end(); ++iter) {
      // -1 reserved to indicate all
      if (*iter == -1) {
        // assert(wid.hiddenZonesSet.size() == 1);
        for (int izn = 0; izn < nzn; ++izn) {
          zn_flag[izn] = 1;
        }
        return;  // no need to process others
      }
      else {
        if ((*iter >= 0)&&(*iter < nzn)) zn_flag[*iter] = 1;
        else {
          if (mpi_rank == 0) cout << "hidden zone index: " << *iter << " is out of range; ignoring" << endl;
        }
      }
    }
  }

  bool isZoneHidden(const int index) const {
    if ((index != -1) && isZoneHidden(-1)) return true;  // hide all zones flag

    set<int>::iterator it = wid.hiddenZonesSet.find(index);
    return (it != wid.hiddenZonesSet.end());
  }

  void ensureHiddenBfFlag(const int nbf) {
    if (nbf != int(bf_hide.size())) {
      bf_hide.resize(nbf);
      for (int ibf=0; ibf<nbf; ++ibf) bf_hide[ibf] = false;
    }
  }

  void setHiddenBfsFromHiddenZones(const int * zone_bf,const int nbf) {
    ensureHiddenBfFlag(nbf);
    if (wid.b_hidesubzones || wid.b_hidezones) {
      const int nzn = zoneNames.size();
      int * hide_zone = new int[nzn];
      for (int izn = 0; izn < nzn; ++izn) hide_zone[izn] = 0;  // default is show

      // in case special catch for all-zones was used
      // from solver zones/subzones are interchangeable, so treat identically
      flagHiddenZones(hide_zone,nzn);
      flagHiddenSubzones(hide_zone,nzn);

      // old approach; may be faster b/c no set lookups...
      // bool * b_hide = NULL;
      // const int nzn = zoneNames.size();
      // if (wid.b_hidesubzones || wid.b_hidezones) {
      //   b_hide = new bool[nzn];
      //   for (int izn = 0; izn < nzn; ++izn)
      //   b_hide[izn] = false; // default to show
      //
      //   set<int>::iterator iter;
      //   for (iter = wid.hiddenSubzonesSet.begin(); iter != wid.hiddenSubzonesSet.end(); ++iter) {
      //     if ((*iter >= 0)&&(*iter < nzn)) {
      //       b_hide[*iter] = true;
      //     }
      //     else {
      //       if (mpi_rank == 0) cout << "WARNING: izn out of bounds: " << *iter << " (nzn=" << nzn << ") File: " << __FILE__ << ", line: " << __LINE__ << endl;
      //     }
      //   }
      //   for (iter = wid.hiddenZonesSet.begin(); iter != wid.hiddenZonesSet.end(); ++iter) {
      //     if (*iter == -1) {
      //       assert(wid.hiddenZonesSet.size() == 1);
      //       for (int izn = 0; izn < nzn; ++izn)
      //       b_hide[izn] = true;
      //     }
      //     else {
      //       if ((*iter >= 0)&&(*iter < nzn)) {
      //         b_hide[*iter] = true;
      //       }
      //       else {
      //         if (mpi_rank == 0) cout << "WARNING: izn out of bounds: " << *iter << " (nzn=" << nzn << ") File: " << __FILE__ << ", line: " << __LINE__ << endl;
      //       }
      //     }
      //   }
      // }

      // now hidden zones are known; create bf_flag
      for (int ibf=0; ibf<nbf; ++ibf) {
        bf_hide[ibf] = (bf_hide[ibf] || hide_zone[zone_bf[ibf]]);  // if set to hide from somewhere else we obey
      }
      DELETE(hide_zone);
    }
  }

  void setHiddenBfsFromCvFlag(const int * cv_flag,const int * cvobf,const int nbf) {
    ensureHiddenBfFlag(nbf);

    // hide bfs attached to flagged cvs
    for (int ibf=0; ibf<nbf; ++ibf) {
      bf_hide[ibf] = (bf_hide[ibf] || (cv_flag[cvobf[ibf]] == 1));  // if set to hide from somewhere else we obey
    }
  }

  float getCanvasPixelDepth(const float xp[3]) const {
    assert(canvas);
    return canvas->calcPixelDepth(xp);
  }

  float getCanvasPixelDepth(const double xp[3]) const {
    assert(canvas);
    return canvas->calcPixelDepth(xp);
  }

  // ==========================================================

  bool hasParticleVar() const {
    if (!wid.b_var_on_particle)
      return false;
    return true;
  }
  string getParticleVar() const {
    assert(wid.b_var_on_particle);
    return wid.var_on_particle;
  }

  double getParticleMagnification() {
    if (wid.b_particle_magnification)
      return wid.particle_magnification;
    else
      return 1.0;
  }

  bool getParticleSize(double &dp) {
    if (wid.b_particle_size) {
      dp = wid.particle_size;
      return true;
    }
    else {
      return false;
    }
  }

  bool getSurfaceVar(string &surfaceVarStr){
    if (!wid.b_var_on_surface)
      return false;

    surfaceVarStr = wid.var_on_surface;
    return true;
  }

  bool hasEdges() {
    return wid.showOpenEdges;
  }

  int hasDynamicEdges() {
    if (wid.showDynamicEdges) {
      return int(wid.dEdgeType);
    }
    else return 0;  // don't do any dynamic edge viz
  }

  bool hasSurfaceMesh() const {
    return wid.showSurfaceMesh;
  }

  bool showBfSurface() const {
    return wid.showBfSurface;
  }

  void addBoundaryFaces(const double (*x_no)[3],const double (*x_bf)[3],const double (*n_bf)[3],const double *data_no,const int *noobf_i,const int *noobf_v,const int * zone_bf,const int nbf,const set<int>& zones_with_data) {
    // to be used when no data or if data_no is specified

    assert(canvas!=NULL);

    const int nzn = zoneNames.size();
    bool *b_data = new bool[nzn];
    if (zones_with_data.size() && (data_no != NULL) ) {
      for (int izn = 0; izn < nzn; ++izn) {
        b_data[izn] = (zones_with_data.find(izn) != zones_with_data.end());
      }
    }
    else if (data_no != NULL) {
      // no data zones specified but data exists, so assume on all zones
      for (int izn = 0; izn < nzn; ++izn) b_data[izn] = true;
    }
    else {
      //zones_with_data not specified; assume all don't have
      for (int izn = 0; izn < nzn; ++izn) b_data[izn] = false;
    }

    const uchar mesh = (wid.showSurfaceMesh) ? 1:0;  // only draw the first edge

    if ((!wid.b_transform_cyl_x)&&(!wid.b_transform_cyl_z)) {
      // Cartesian image

      double data[3];
      FOR_IBF {
        double data_bf = 0.0;
        const int izn = zone_bf[ibf];
        assert((izn >= 0)&&(izn < 32768));

        if (b_data[izn]) {  // data on the surface

          // build data at the x_bf...
          double wgt_bf = 0.0;
          int ino1 = noobf_v[noobf_i[ibf+1]-1];
          for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
            const int ino0 = ino1;
            ino1 = noobf_v[nob];
            const double this_wgt = DIST(x_no[ino0],x_no[ino1]);
            data_bf += this_wgt*(data_no[ino0]+data_no[ino1]);
            wgt_bf += this_wgt;
          }
          assert(wgt_bf > 0.0);
          data_bf /= 2.0*wgt_bf;

          // triangulate voronoi face
          ino1 = noobf_v[noobf_i[ibf+1]-1];
          for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
            const int ino0 = ino1;
            ino1 = noobf_v[nob];
            data[0] = data_no[ino0];
            data[1] = data_no[ino1];
            data[2] = data_bf;
            canvas->addSurfaceTri(x_no[ino0],x_no[ino1],x_bf[ibf],izn,bf_hide[ibf],mesh,false,&n_bf[ibf][0],data);
          }

        }
        else {

          // no data on this surface tri...
          // triangulate voronoi face
          int ino1 = noobf_v[noobf_i[ibf+1]-1];
          for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
            const int ino0 = ino1;
            ino1 = noobf_v[nob];
            canvas->addSurfaceTri(x_no[ino0],x_no[ino1],x_bf[ibf],izn,bf_hide[ibf]||wid.b_hide_surfaces_without_data,mesh,false,n_bf[ibf]);
          }

        }
      }
    }
    else if (wid.b_transform_cyl_x || wid.b_transform_cyl_z) {

      // to allow allocation once, figure out the largest boundary nnob size
      // outside of the FOR_IBF loop...

      double x_bf_t[3];
      int nnob_max = 0;
      FOR_IBF nnob_max = max(nnob_max,noobf_i[ibf+1]-noobf_i[ibf]);
      double (*x_no_t)[3] = new double[nnob_max][3];

      double data[3];
      FOR_IBF {
	const int izn = zone_bf[ibf];
	assert((izn >= 0)&&(izn < 32768));
	//TODO have wid(?) do the transform for you...
	coordinateTransform(x_bf_t,x_bf[ibf]);

	const int nob_f = noobf_i[ibf];
	for (int nob = nob_f; nob != noobf_i[ibf+1]; ++nob) {
	  const int ino = noobf_v[nob];
	  coordinateTransformAndAdjustTheta(x_no_t[nob-nob_f], x_no[ino], x_bf_t);
	}

        // triangulate voronoi face

	if (b_data[izn]) {

          // build data at the x_bf...
          double data_bf = 0.0;
	  double wgt_bf = 0.0;
	  int ino1 = noobf_v[noobf_i[ibf+1]-1];
	  for (int nob = nob_f; nob != noobf_i[ibf+1]; ++nob) {
	    const int ino0 = ino1;
	    ino1 = noobf_v[nob];
	    const double this_wgt = DIST(x_no[ino0],x_no[ino1]);
	    data_bf += this_wgt*(data_no[ino0]+data_no[ino1]);
	    wgt_bf += this_wgt;
	  }
	  assert(wgt_bf > 0.0);
	  data_bf /= 2.0*wgt_bf;

          int nob1 = noobf_i[ibf+1]-1-nob_f;
          ino1 = noobf_v[noobf_i[ibf+1]-1];
          for (int nob = nob_f; nob != noobf_i[ibf+1]; ++nob) {
            const int nob0 = nob1;
            const int ino0 = ino1;
            nob1 = nob-nob_f;
            ino1 = noobf_v[nob];
            data[0] = data_no[ino0];
            data[1] = data_no[ino1];
            data[2] = data_bf;
            canvas->addSurfaceTri(x_no_t[nob0],x_no_t[nob1],x_bf_t,izn,bf_hide[ibf],mesh,false,NULL,data);
          }

        }
        else {

          // no data on this surface...
          int nob1 = noobf_i[ibf+1]-1-nob_f;
          int ino1 = noobf_v[noobf_i[ibf+1]-1];
          for (int nob = nob_f; nob != noobf_i[ibf+1]; ++nob) {
            const int nob0 = nob1;
            const int ino0 = ino1;
            nob1 = nob-nob_f;
            ino1 = noobf_v[nob];
            canvas->addSurfaceTri(x_no_t[nob0],x_no_t[nob1],x_bf_t,izn,bf_hide[ibf]|wid.b_hide_surfaces_without_data,mesh,false);
          }

        }

      } // FOR_IBF

      delete[] x_no_t;

    }

    delete[] b_data;

  }

  void addBoundaryTris(const double (*x_no)[3],const double * data_no,const int nno,const int (*noote)[4],const int nte,const int index) {

    assert(canvas!=NULL);

    assert(data_no != NULL);
    MiscUtils::dumpRange(data_no,nno,"addBoundaryTris range");

    const uchar mesh = (wid.showSurfaceMesh) ? 7:0;  // draw all3 edges or none

    if ((!wid.b_transform_cyl_x)&&(!wid.b_transform_cyl_z)) {

      // Cartesian image

      double data[3];

      // for now, render ALL tet tris...

      for (int ite = 0; ite < nte; ++ite) {
        const int ino0 = noote[ite][0]; assert((ino0 >= 0)&&(ino0 < nno));
        const int ino1 = noote[ite][1]; assert((ino1 >= 0)&&(ino1 < nno));
        const int ino2 = noote[ite][2]; assert((ino2 >= 0)&&(ino2 < nno));
        const int ino3 = noote[ite][3]; assert((ino3 >= 0)&&(ino3 < nno));
        // render all 4 faces...
        // 0...
        data[0] = data_no[ino0];
        data[1] = data_no[ino1];
        data[2] = data_no[ino2];
        canvas->addSurfaceTri(x_no[ino0],x_no[ino1],x_no[ino2],index,false,mesh,false,NULL,data);
        // 1...
        data[0] = data_no[ino1];
        data[1] = data_no[ino3];
        data[2] = data_no[ino2];
        canvas->addSurfaceTri(x_no[ino1],x_no[ino3],x_no[ino2],index,false,mesh,false,NULL,data);
        // 2...
        data[0] = data_no[ino2];
        data[1] = data_no[ino3];
        data[2] = data_no[ino0];
        canvas->addSurfaceTri(x_no[ino2],x_no[ino3],x_no[ino0],index,false,mesh,false,NULL,data);
        // 3...
        data[0] = data_no[ino0];
        data[1] = data_no[ino3];
        data[2] = data_no[ino1];
        canvas->addSurfaceTri(x_no[ino0],x_no[ino3],x_no[ino1],index,false,mesh,false,NULL,data);
      }

    }
    else if (wid.b_transform_cyl_x || wid.b_transform_cyl_z) {

      // Cylindrical unwrap

      double data[3];
      double x_no_t[4][3];

      // for now, render ALL tet tris...

      for (int ite = 0; ite < nte; ++ite) {
        const int ino0 = noote[ite][0]; assert((ino0 >= 0)&&(ino0 < nno));
        const int ino1 = noote[ite][1]; assert((ino1 >= 0)&&(ino1 < nno));
        const int ino2 = noote[ite][2]; assert((ino2 >= 0)&&(ino2 < nno));
        const int ino3 = noote[ite][3]; assert((ino3 >= 0)&&(ino3 < nno));
	coordinateTransform(x_no_t[0],x_no[ino0]);
	coordinateTransform(x_no_t[1],x_no[ino1]);
	coordinateTransform(x_no_t[2],x_no[ino2]);
	coordinateTransform(x_no_t[3],x_no[ino3]);
        // render all 4 faces...
        // 0...
        data[0] = data_no[ino0];
        data[1] = data_no[ino1];
        data[2] = data_no[ino2];
	adjustThetaTri(x_no_t[0],x_no_t[1],x_no_t[2]);
        canvas->addSurfaceTri(x_no_t[0],x_no_t[1],x_no_t[2],index,false,mesh,false,NULL,data);
        // 1...
        data[0] = data_no[ino1];
        data[1] = data_no[ino3];
        data[2] = data_no[ino2];
	adjustThetaTri(x_no_t[1],x_no_t[3],x_no_t[2]);
        canvas->addSurfaceTri(x_no_t[1],x_no_t[3],x_no_t[2],index,false,mesh,false,NULL,data);
        // 2...
        data[0] = data_no[ino2];
        data[1] = data_no[ino3];
        data[2] = data_no[ino0];
	adjustThetaTri(x_no_t[2],x_no_t[3],x_no_t[0]);
        canvas->addSurfaceTri(x_no_t[2],x_no_t[3],x_no_t[0],index,false,mesh,false,NULL,data);
        // 3...
        data[0] = data_no[ino0];
        data[1] = data_no[ino3];
        data[2] = data_no[ino1];
	adjustThetaTri(x_no_t[0],x_no_t[3],x_no_t[1]);
        canvas->addSurfaceTri(x_no_t[0],x_no_t[3],x_no_t[1],index,false,mesh,false,NULL,data);
      }

    }

  }

  void addBoundaryFacesWithVectors(const double (*x_no)[3],const double (*x_bf)[3],const double (*n_bf)[3],const double (*bf_data)[3],const int *noobf_i,const int *noobf_v,const int *znobf,const double *area_bf,const int nbf) {

    // vectors will be drawn on top of no-data surfaces
    // first add surfaces and compute areas

    assert(canvas!=NULL);
    assert(b_rmax);

    const uchar mesh = (wid.showSurfaceMesh) ? 1:0;  // draw only on first edge

    double my_area = 0.0;
    double my_mag_max = 0.0;
    FOR_IBF {
      const int izn = znobf[ibf];
      assert((izn >= 0)&&(izn < 32768));

      // triangulate voronoi face and draw as a surface to the canvas
      int ino1 = noobf_v[noobf_i[ibf+1]-1];
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino0 = ino1;
        ino1 = noobf_v[nob];
        canvas->addSurfaceTri(x_no[ino0],x_no[ino1],x_bf[ibf],izn,bf_hide[ibf],mesh,false,&n_bf[ibf][0]);
      }
      if (!bf_hide[ibf]) {
        my_area += area_bf[ibf];
        my_mag_max = max(my_mag_max,MAG(bf_data[ibf]));
      }
    }

    double area;
    MPI_Allreduce(&my_area,&area,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
    double mag_max;
    MPI_Allreduce(&my_mag_max,&mag_max,1,MPI_DOUBLE,MPI_MAX,mpi_comm);

    if ((area > 0.0)&&(mag_max > 0.0)) {

      // user/app controls (defaults are 1.0)...
      const double factor = 0.2*(double)rmax*getVectorScale()/mag_max;
      const double target_area = area/(250.0*getVectorDensity());

      if (factor > 0.0) {
        vector<double> x0Vec;
        vector<double> dxVec;
        vector<int> zoneVec;
        double current_area = 0.0;
        FOR_IBF {

          const int izn = znobf[ibf];
          assert((izn >= 0)&&(izn < 32768));
          if (!bf_hide[ibf]) {
            current_area += area_bf[ibf];
            if (current_area > target_area) {
              current_area = 0.0; // reset current representative area

              FOR_I3 x0Vec.push_back(x_bf[ibf][i]);
              FOR_I3 dxVec.push_back(bf_data[ibf][i]);
              zoneVec.push_back(izn);
            }
          }
        }
        if (x0Vec.size() > 0) addSurfaceVectors((double(*)[3])(&x0Vec[0]),(double(*)[3])(&dxVec[0]),&zoneVec[0],(int)x0Vec.size()/3,factor);
      }
    }
  }

  void addBoundaryFacesWithPartialVectors(const double (*x_no)[3],const double (*x_bf)[3],const double (*n_bf)[3],const double (*bf_data)[3],
					  const int *noobf_i,const int *noobf_v,const int *znobf,const double *area_bf,const int nbf,const set<int>& zones_with_data) {

    // TODO: no transform implemented here...

    assert(b_rmax);
    assert(canvas!=NULL);

    const int nzn = zoneNames.size();
    bool *b_data = new bool[nzn];
    for (int izn = 0; izn < nzn; ++izn)
      b_data[izn] = (zones_with_data.find(izn) != zones_with_data.end());

    double my_area = 0.0;
    double my_mag_max = 0.0;

    const uchar mesh = (wid.showSurfaceMesh) ? 1:0;  // render mesh on first edge

    FOR_IBF {
      const int izn = znobf[ibf];
      assert((izn >= 0)&&(izn < 32768));

      // triangulate voronoi face and draw as non-data surface
      int ino1 = noobf_v[noobf_i[ibf+1]-1];
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
        const int ino0 = ino1;
        ino1 = noobf_v[nob];
        canvas->addSurfaceTri(x_no[ino0],x_no[ino1],x_bf[ibf],izn,bf_hide[ibf],mesh,false,&n_bf[ibf][0]);
      }

      if (b_data[izn]&&(!bf_hide[ibf])) {
        my_area += area_bf[ibf];
        my_mag_max = max(my_mag_max,MAG(bf_data[ibf]));
      }
    }

    double area;
    MPI_Allreduce(&my_area,&area,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
    double mag_max;
    MPI_Allreduce(&my_mag_max,&mag_max,1,MPI_DOUBLE,MPI_MAX,mpi_comm);

    if ((area > 0.0)&&(mag_max > 0.0)) {

      // user/app controls (defaults are 1.0)...
      const double factor = 0.2*(double)rmax*getVectorScale()/mag_max;
      const double target_area = area/(250.0*getVectorDensity());

      if (factor > 0.0) {

        vector<double> x0Vec;
        vector<double> dxVec;
        vector<int> zoneVec;
        double current_area = 0.0;
        FOR_IBF {

          const int izn = znobf[ibf];
          assert((izn >= 0)&&(izn < 32768));
          if (b_data[izn]&&(!bf_hide[ibf])) {
            current_area += area_bf[ibf];
            if (current_area > target_area) {
              current_area = 0.0; // reset current representative area

              FOR_I3 x0Vec.push_back(x_bf[ibf][i]);
              FOR_I3 dxVec.push_back(bf_data[ibf][i]);
              zoneVec.push_back(izn);

            }
          }
        }
        if (x0Vec.size() > 0)
          addSurfaceVectors((double(*)[3])(&x0Vec[0]),(double(*)[3])(&dxVec[0]),&zoneVec[0],(int)x0Vec.size()/3,factor);
      }
    }

    delete[] b_data;
  }

  void addInternalFace(const double (*x_no)[3],const double x_fa[3],const int nnof,const int * noofa_v) {

    assert(canvas != NULL);

    // triangulate voronoi face
    const uchar mesh = (wid.showSurfaceMesh) ? 1:0;  // render mesh on first edge

    if (wid.b_transform_cyl_x || wid.b_transform_cyl_z) {
      double x_fa_t[3];
      coordinateTransform(x_fa_t, x_fa);

      double (*x_no_t)[3] = new double[nnof][3];
      for (int nof = 0; nof < nnof; ++nof) {
	const int ino = noofa_v[nof];
	coordinateTransformAndAdjustTheta(x_no_t[nof], x_no[ino], x_fa_t);
      }

      int ino1 = nnof-1;
      for (int nof = 0; nof < nnof; ++nof) {
	const int ino0 = ino1;
	ino1 = nof;
	canvas->addInternalTri(x_no_t[ino0],x_no_t[ino1],x_fa_t,32767,mesh,true);
      }

      delete[] x_no_t;
    }
    else {
      int ino1 = noofa_v[nnof-1];
      for (int nof = 0; nof < nnof; ++nof) {
	const int ino0 = ino1;
	ino1 = noofa_v[nof];
	canvas->addInternalTri(x_no[ino0],x_no[ino1],x_fa,32767,mesh,true);
      }
    }
  }

  void addInternalFaceFlip(const double (*x_no)[3],const double x_fa[3],const int nnof,const int * noofa_v) {

    assert(canvas != NULL);

    // triangulate voronoi face
    const uchar mesh = (wid.showSurfaceMesh) ? 1:0;  // render mesh on first edge

    if (wid.b_transform_cyl_x || wid.b_transform_cyl_z) {
      double x_fa_t[3];
      coordinateTransform(x_fa_t, x_fa);

      double (*x_no_t)[3] = new double[nnof][3];
      for (int nof = 0; nof < nnof; ++nof) {
	const int ino = noofa_v[nof];
	coordinateTransformAndAdjustTheta(x_no_t[nof], x_no[ino], x_fa_t);
      }

      int ino1 = 0;
      for (int nof = nnof-1; nof >= 0; --nof) {
	const int ino0 = ino1;
	ino1 = nof;
	canvas->addInternalTri(x_no_t[ino0],x_no_t[ino1],x_fa_t,32767,mesh,true);
      }

      delete[] x_no_t;
    }
    else {
      int ino1 = noofa_v[0];
      for (int nof = nnof-1; nof >= 0; --nof) {
	const int ino0 = ino1;
	ino1 = noofa_v[nof];
	canvas->addInternalTri(x_no[ino0],x_no[ino1],x_fa,32767,mesh,true);
      }
    }
  }

  void addSurfaceTris(const double (*xsp)[3],const int (*spost)[3],const int * znost,const int8 * ist_global_and_bits,const int nst) {

    // simplest routine to add tris to a mesh...

    assert(canvas != NULL);

    // the set is sorted, so we want the last entry. Set has no back(), so
    // do this with the reverse iterator...

    int hidden_isz_max = -1;
    if (!wid.hiddenSubzonesSet.empty())
      hidden_isz_max = *wid.hiddenSubzonesSet.rbegin();

    assert(hidden_isz_max < 65535);
    int * hidden_isz_flag = new int[hidden_isz_max+1];
    for (int isz = 0; isz <= hidden_isz_max; ++isz)
      hidden_isz_flag[isz] = 0;

    for (set<int>::iterator it = wid.hiddenSubzonesSet.begin(); it != wid.hiddenSubzonesSet.end(); ++it) {
      const int isz = *it;
      if ((isz >= 0)&&(isz <= hidden_isz_max)) hidden_isz_flag[isz] = 1;
      else {
        if (mpi_rank == 0) cout << "hidden subzone index: " << isz << " is out of range; ignoring" << endl;
      }
    }
    for (set<int>::iterator it = wid.hiddenZonesSet.begin(); it != wid.hiddenZonesSet.end(); ++it) {
      const int isz = *it;
      if ((isz >= 0)&&(isz <= hidden_isz_max)) hidden_isz_flag[isz] = 1;
      else {
        if (mpi_rank == 0) cout << "hidden zone index: " << isz << " is out of range; ignoring" << endl;
      }
    }

    // add the surface tris...

    const uchar mesh = (wid.showSurfaceMesh) ? 7:0;  // show all the edges of the tri

    for (int ist = 0; ist < nst; ++ist) {
      const int index = znost[ist];
      assert((index >= 0)&&(index < 32768));
      const int bits = int(ist_global_and_bits[ist]>>52);
      const bool hide = bits||!((index > hidden_isz_max)||(hidden_isz_flag[index] == 0));

      // default is to use xsp node locations
      double x[3][3];
      FOR_I3 {
	FOR_J3 {
	  x[i][j] = xsp[spost[ist][i]][j];
	}
      }
      // adjust if a transform is specified
      if (wid.b_transform_cyl_x || wid.b_transform_cyl_z) {
	FOR_I3 {
	  const int isp = spost[ist][i];
	  coordinateTransform(x[i], xsp[isp]);
	}
	adjustThetaTri(x[0],x[1],x[2]);
      }

      canvas->addSurfaceTri(x[0],x[1],x[2],index,hide,mesh);
    }

    delete[] hidden_isz_flag;
  }

  void addSurfaceTris(const double (*xsp)[3],const int (*spost)[3],const int * znost,const int8 * ist_global_and_bits,const int nst,set<int>& skipZoneSet) {

    assert(canvas != NULL);

    int hidden_isz_max = -1;
    if (!wid.hiddenSubzonesSet.empty())
      hidden_isz_max = *wid.hiddenSubzonesSet.rbegin();

    assert(hidden_isz_max < 65535);
    int * hidden_isz_flag = new int[hidden_isz_max+1];
    for (int isz = 0; isz <= hidden_isz_max; ++isz)
      hidden_isz_flag[isz] = 0;

    for (set<int>::iterator it = wid.hiddenSubzonesSet.begin(); it != wid.hiddenSubzonesSet.end(); ++it) {
      const int isz = *it;
      if ((isz >= 0)&&(isz <= hidden_isz_max)) hidden_isz_flag[isz] = 1;
      else {
        if (mpi_rank == 0) cout << "hidden subzone index: " << isz << " is out of range; ignoring" << endl;
      }
    }
    for (set<int>::iterator it = wid.hiddenZonesSet.begin(); it != wid.hiddenZonesSet.end(); ++it) {
      const int isz = *it;
      if ((isz >= 0)&&(isz <= hidden_isz_max)) hidden_isz_flag[isz] = 1;
      else {
        if (mpi_rank == 0) cout << "hidden zone index: " << isz << " is out of range; ignoring" << endl;
      }
    }

    // add the surface tris...

    const uchar mesh = (wid.showSurfaceMesh) ? 7:0;  // show all edges

    for (int ist = 0; ist < nst; ++ist) {
      const int index = znost[ist];
      assert((index >= 0)&&(index < 32768));
      if (skipZoneSet.find(index) == skipZoneSet.end()) {
        const int bits = int(ist_global_and_bits[ist]>>52);
        const bool hide = bits||!((index > hidden_isz_max)||(hidden_isz_flag[index] == 0));

	// default is to use xsp node locations
	double x[3][3];
	FOR_I3 {
	  FOR_J3 {
	    x[i][j] = xsp[spost[ist][i]][j];
	  }
	}
	// adjust if a transform is specified
	if (wid.b_transform_cyl_x || wid.b_transform_cyl_z) {
	  FOR_I3 {
	    const int isp = spost[ist][i];
	    coordinateTransform(x[i], xsp[isp]);
	  }
	  adjustThetaTri(x[0],x[1],x[2]);
	}

	canvas->addSurfaceTri(x[0],x[1],x[2],index,hide,mesh);
      }
    }

    delete[] hidden_isz_flag;
  }

  void addSurfaceTris(const double (*xsp)[3],const int (*spost)[3],const int * znost,const int nst) {

    // simplest routine to add tris to a mesh...
    // currently this is only used byStitch; Charles uses version with bits&rank, Surfer uses version with szozn_i

    assert(canvas != NULL);

    // the set is sorted, so we want the last entry. Set has no back(), so
    // do this with the reverse iterator...

    int hidden_isz_max = -1;
    if (!wid.hiddenSubzonesSet.empty()) hidden_isz_max = *wid.hiddenSubzonesSet.rbegin();
    if (!wid.hiddenZonesSet.empty()) hidden_isz_max = max(hidden_isz_max,*wid.hiddenZonesSet.rbegin());

    assert(hidden_isz_max < 65535);
    int * hidden_isz_flag = new int[hidden_isz_max+1];
    for (int isz = 0; isz <= hidden_isz_max; ++isz)
      hidden_isz_flag[isz] = 0;

    for (set<int>::iterator it = wid.hiddenSubzonesSet.begin(); it != wid.hiddenSubzonesSet.end(); ++it) {
      const int isz = *it;
      if ((isz >= 0)&&(isz <= hidden_isz_max)) hidden_isz_flag[isz] = 1;
      else {
        if (mpi_rank == 0) cout << "hidden subzone index: " << isz << " is out of range; ignoring" << endl;
      }
    }
    for (set<int>::iterator it = wid.hiddenZonesSet.begin(); it != wid.hiddenZonesSet.end(); ++it) {
      const int isz = *it;
      if ((isz >= 0)&&(isz <= hidden_isz_max)) hidden_isz_flag[isz] = 1;
      else {
        if (mpi_rank == 0) cout << "hidden zone index: " << isz << " is out of range; ignoring" << endl;
      }
    }

    const uchar mesh = (wid.showSurfaceMesh) ? 7:0; // either show all edges or none

    // add the surface tris...
    for (int ist = 0; ist < nst; ++ist) {
      const int index = znost[ist];
      assert((index >= 0)&&(index < 32768));
      const bool hide = !((index > hidden_isz_max)||(hidden_isz_flag[index] == 0));

      // default is to use xsp node locations
      double x[3][3];
      FOR_I3 {
        FOR_J3 {
          x[i][j] = xsp[spost[ist][i]][j];
        }
      }
      // adjust if a transform is specified
      if (wid.b_transform_cyl_x || wid.b_transform_cyl_z) {
        FOR_I3 {
          const int isp = spost[ist][i];
	  coordinateTransform(x[i], xsp[isp]);
        }
	adjustThetaTri(x[0],x[1],x[2]);
      }

      canvas->addSurfaceTri(x[0],x[1],x[2],index,hide,mesh);
    }

    delete[] hidden_isz_flag;
  }

  void addSurfaceTrisWithData(const double (*xsp)[3],const double *dsp,const int (*spost)[3],const int * znost,const int nst) {

    // simplest routine to add tris to a mesh...
    assert(canvas != NULL);

    // the set is sorted, so we want the last entry. Set has no back(), so
    // do this with the reverse iterator...

    int hidden_isz_max = -1;
    if (!wid.hiddenSubzonesSet.empty())
      hidden_isz_max = *wid.hiddenSubzonesSet.rbegin();

    assert(hidden_isz_max < 65535);
    int * hidden_isz_flag = new int[hidden_isz_max+1];
    for (int isz = 0; isz <= hidden_isz_max; ++isz)
      hidden_isz_flag[isz] = 0;

    for (set<int>::iterator it = wid.hiddenSubzonesSet.begin(); it != wid.hiddenSubzonesSet.end(); ++it) {
      const int isz = *it;
      if ((isz >= 0)&&(isz <= hidden_isz_max)) hidden_isz_flag[isz] = 1;
      else {
        if (mpi_rank == 0) cout << "hidden subzone index: " << isz << " is out of range; ignoring" << endl;
      }
    }
    for (set<int>::iterator it = wid.hiddenZonesSet.begin(); it != wid.hiddenZonesSet.end(); ++it) {
      const int isz = *it;
      if ((isz >= 0)&&(isz <= hidden_isz_max)) hidden_isz_flag[isz] = 1;
      else {
        if (mpi_rank == 0) cout << "hidden zone index: " << isz << " is out of range; ignoring" << endl;
      }
    }

    // add the surface tris...
    const uchar mesh = (wid.showSurfaceMesh) ? 7:0;  // all three edges
    double * null_ptr = NULL;  // use tri nodes to compute normal

    for (int ist = 0; ist < nst; ++ist) {
      const int index = znost[ist];
      assert((index >= 0)&&(index < 32768));
      const bool hide = !((index > hidden_isz_max)||(hidden_isz_flag[index] == 0));
      const double data[3] = {dsp[spost[ist][0]],dsp[spost[ist][1]],dsp[spost[ist][2]]};
      canvas->addSurfaceTri(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],index,hide,mesh,false,null_ptr,&data[0]);
    }

    delete[] hidden_isz_flag;
  }

  void addSurfaceTrisWithDataFromCim(const double (*xsp)[3],const int (*spost)[3],const int * znost,const IntFlag &szost,const IntFlag &szozn_i,const int nst) {
    UNUSED(znost);

    assert(canvas!=NULL);
    // add the surface tris...

    // flag stores which zones should be hidden (i.e., not added)
    const int nzn = szozn_i.getLength() - 1;
    const int nsz = szozn_i[nzn];
    IntFlag show_sz_flag(nsz);
    show_sz_flag.setAll(1);
    if (wid.b_hidesubzones || wid.b_hidezones) {
      set<int>::iterator it;
      for (it = wid.hiddenZonesSet.begin(); it != wid.hiddenZonesSet.end(); ++it) {
        if (*it == -1) {
          assert(wid.hiddenZonesSet.size() == 1);
          for (int isz = 0; isz < nsz; ++isz) {
            show_sz_flag[isz] = 0;
          }
        }
        else {
          const int izn = *it;
          if ((izn >= 0)&&(izn < nzn)) {
            for (int iszn = szozn_i[izn]; iszn < szozn_i[izn+1]; ++iszn) show_sz_flag[iszn] = 0;
          }
          else {
            if (mpi_rank == 0) cout << "WARNING: izn out of bounds: " << izn << " (nzn=" << nzn << ") File: " << __FILE__ << ", line: " << __LINE__ << endl;
          }
        }
      }
      for (it = wid.hiddenSubzonesSet.begin(); it != wid.hiddenSubzonesSet.end(); ++it) {
        const int iszn = *it;
        if ((iszn >= 0)&&(iszn < nsz)) {
          show_sz_flag[iszn] = 0;
        }
        else {
          if (mpi_rank == 0) cout << "WARNING: isnz out of bounds: " << iszn << " (nsz=" << nsz << ") File: " << __FILE__ << ", line: " << __LINE__ << endl;
        }
      }
    }

    const uchar mesh = (wid.showSurfaceMesh) ? 7:0;  // highlight all edges

    // everybody renders the tris
    for (int ist = 0; ist < nst; ++ist) {
      const int index = szost[ist];  // subzone index
      canvas->addSurfaceTri(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],index,!(index<nsz&&index>=0&&show_sz_flag[index]),mesh);
    }

    if (wid.b_var_on_surface){
      if (wid.b_range_on_surface){
        canvas->setSurfaceDataRange(wid.range_on_surface[0],wid.range_on_surface[1]);
      }

      if (filename_volumeData != "") {

        //Use LesImageMapper to read surface data...
        assert(lim==NULL);  //this will happen before plane cims are built...
        if (filename_volumeData_2=="")
          lim = new LesImageMapper(filename_volumeData);
        else
          lim = new LesImageMapper(filename_volumeData,filename_volumeData_2);

        //read or build a CvImageMap...

        //pass blank planes for unique cim name generations
        //must include geom plane blanking if active
        vector<PlaneData<float> > blank_planes;
        canvas->getBlankPlanes(blank_planes);

        string cimName = lim->buildCvImageMapFilenameSurface(canvas,bbox,nst,blank_planes,show_sz_flag);
        int readCim = cimS.read(cimName.c_str());
        if (readCim<0){
          //cimS was not found, must build it

          //Get a list of zones that "won" the depth test to speed up (limit) the file io
          IntFlag bfzone_flag(nzn); bfzone_flag.setAll(0);
          canvas->getActiveZones(bfzone_flag,szozn_i);

          lim->buildCvImageMapSurface(cimS,canvas,bfzone_flag,lengthscaleName);
          cimS.write(cimName.c_str());
        }
      }
    }
  }

  void addSurfaceTrisWithData(const double (*xsp)[3],const int (*spost)[3],const IntFlag &szost,const IntFlag &szozn_i,const int nst,const double* surface_data = NULL,const int surface_data_type = 0) {

    assert(canvas!=NULL);

    // flag stores which zones should be hidden (i.e., not added)...

    const int nzn = szozn_i.getLength() - 1;
    const int nsz = szozn_i[nzn];
    IntFlag show_sz_flag(nsz);
    show_sz_flag.setAll(1);
    if (wid.b_hidesubzones || wid.b_hidezones) {
      set<int>::iterator it;
      for (it = wid.hiddenZonesSet.begin(); it != wid.hiddenZonesSet.end(); ++it) {
        if (*it == -1) {
          assert(wid.hiddenZonesSet.size() == 1);
          for (int isz = 0; isz < nsz; ++isz) {
            show_sz_flag[isz] = 0;
          }
        }
        else {
          const int izn = *it;
          if ((izn >= 0)&&(izn < nzn)) {
            for (int iszn = szozn_i[izn]; iszn < szozn_i[izn+1]; ++iszn) show_sz_flag[iszn] = 0;
          }
          else {
            if (mpi_rank == 0) cout << "WARNING: izn out of bounds: " << izn << " (nzn=" << nzn << ") File: " << __FILE__ << ", line: " << __LINE__ << endl;
          }
        }
      }
      for (it = wid.hiddenSubzonesSet.begin(); it != wid.hiddenSubzonesSet.end(); ++it) {
        const int iszn = *it;
        if ((iszn >= 0)&&(iszn < nsz)) {
          show_sz_flag[iszn] = 0;
        }
        else {
          if (mpi_rank == 0) cout << "WARNING: isnz out of bounds: " << iszn << " (nsz=" << nsz << ") File: " << __FILE__ << ", line: " << __LINE__ << endl;
        }
      }
    }

    double * null_ptr = NULL;  // compute normal based on tri nodes pushed into the canvas
    const uchar mesh = (wid.showSurfaceMesh) ? 7:0;  // if shown show all edges of tri
    const bool b_has_surface_data = (surface_data != NULL)&&(surface_data_type != 0)&&wid.b_var_on_surface;


    if(wid.b_proj_cyl&&xsp){
      //collapse surface about cyl-theta to R-Z
      int nsp = -1; //compute nsp
      for (int ist = 0; ist < nst; ++ist) {
        FOR_I3 nsp = max(nsp,spost[ist][i]);   
      }
      ++nsp;//one more than max index required
      double (*xsp_local)[3] = new double[nsp][3];
      double axis_mag = MAG(wid.proj_cyl_np);
      const double axis_unit[3] = {wid.proj_cyl_np[0]/axis_mag, wid.proj_cyl_np[1]/axis_mag, wid.proj_cyl_np[2]/axis_mag};
      for (int isp=0; isp<nsp; ++isp){
        xsp_local[isp][0] = 0.0;
        double x_rel[3] = DIFF(xsp[isp],wid.proj_cyl_xp);
        xsp_local[isp][2] = DOT_PRODUCT(x_rel,axis_unit);
        xsp_local[isp][1] = sqrt(max(0.0,DOT_PRODUCT(x_rel,x_rel)-xsp_local[isp][2]*xsp_local[isp][2]));
      }
      canvas->setProjPlane();
      if (b_has_surface_data) {
        if (wid.b_range_on_surface) canvas->setSurfaceDataRange(wid.range_on_surface[0],wid.range_on_surface[1]);

        if (surface_data_type == 1) {
          for (int ist = 0; ist < nst; ++ist) {
            const int index = szost[ist];
            assert((index >= 0)&&(index < 32768));
            const double data[3] = {surface_data[spost[ist][0]],surface_data[spost[ist][1]],surface_data[spost[ist][2]]};
            canvas->addSurfaceTriThickened(xsp_local[spost[ist][0]],xsp_local[spost[ist][1]],xsp_local[spost[ist][2]],index,!(index<nsz&&index>=0&&show_sz_flag[index]),mesh,false,null_ptr,&data[0]);
          }
        }
        else if (surface_data_type == 2) {
          for (int ist = 0; ist < nst; ++ist) {
            const int index = szost[ist];
            assert((index >= 0)&&(index < 32768));
            const double * data = &surface_data[ist];
            canvas->addSurfaceTriThickened(xsp_local[spost[ist][0]],xsp_local[spost[ist][1]],xsp_local[spost[ist][2]],index,!(index<nsz&&index>=0&&show_sz_flag[index]),mesh,false,null_ptr,data,true);
          }
        }
        else {
          CWARN("unrecognized surface data type");
          assert(0); // TODO surface_data_type 3 and 4
        }
      }
      else {
        // no surface data
        for (int ist = 0; ist < nst; ++ist) {
          const int index = szost[ist];  // subzone index
          assert((index >= 0)&&(index < 32768));
          canvas->addSurfaceTriThickened(xsp_local[spost[ist][0]],xsp_local[spost[ist][1]],xsp_local[spost[ist][2]],index,!(index<nsz&&index>=0&&show_sz_flag[index]),mesh);
        }
      }

      delete[] xsp_local;
    }
    else {
      if (b_has_surface_data) {
        if (wid.b_range_on_surface) canvas->setSurfaceDataRange(wid.range_on_surface[0],wid.range_on_surface[1]);

        if (surface_data_type == 1) {
          for (int ist = 0; ist < nst; ++ist) {
            const int index = szost[ist];
            assert((index >= 0)&&(index < 32768));
            const double data[3] = {surface_data[spost[ist][0]],surface_data[spost[ist][1]],surface_data[spost[ist][2]]};
            canvas->addSurfaceTri(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],index,!(index<nsz&&index>=0&&show_sz_flag[index]),mesh,false,null_ptr,&data[0]);
          }
        }
        else if (surface_data_type == 2) {
          for (int ist = 0; ist < nst; ++ist) {
            const int index = szost[ist];
            assert((index >= 0)&&(index < 32768));
            const double * data = &surface_data[ist];
            canvas->addSurfaceTri(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],index,!(index<nsz&&index>=0&&show_sz_flag[index]),mesh,false,null_ptr,data,true);
          }
        }
        else {
          CWARN("unrecognized surface data type");
          assert(0); // TODO surface_data_type 3 and 4
        }
      }
      else {
        // no surface data
        for (int ist = 0; ist < nst; ++ist) {
          const int index = szost[ist];  // subzone index
          assert((index >= 0)&&(index < 32768));
          canvas->addSurfaceTri(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],index,!(index<nsz&&index>=0&&show_sz_flag[index]),mesh);
        }
      }
    }
  }

  void addPlaneDataRvvPositive(const double * phi,const double (*x_vv)[3],const double *r_vv,const int n,const double factor) {
    // factor is a number that multiples the passed r_vv to determine the max sphere of potential influence of the
    // the vd's...
    if (wid.b_transform_cyl_x) {
      double (*x_vv_t)[3] = new double[n][3];
      double *r_vv_t = new double[n];
      for (int i = 0; i < n; ++i) {
        if (r_vv[i] <= 0.0) {
          r_vv_t[i] = r_vv[i];
        }
        else {
          // (x,y,z) -> (x,r,r0*theta)
	  coordinateTransform(x_vv_t[i], x_vv[i]);
          // for r less than wid.r0, the space will be stretched in the
          // azimuthal direction. Increase r_vv proportionally to avoid missing pixels,
          // up to a maximum factor of 2 (otherwise gets expensive)...
          if (x_vv_t[i][1] > 0.5*wid.r0) {
            r_vv_t[i] = r_vv[i]*max(wid.r0/x_vv_t[i][1],1.0);
          }
          else {
            r_vv_t[i] = r_vv[i]*2.0;
          }
        }
      }
      canvas->addPlaneDataRvvPositive(phi,x_vv_t,r_vv_t,n,factor);
      delete[] x_vv_t;
      delete[] r_vv_t;
    }
    else if (wid.b_transform_cyl_z) {
      double (*x_vv_t)[3] = new double[n][3];
      double *r_vv_t = new double[n];
      for (int i = 0; i < n; ++i) {
        if (r_vv[i] <= 0.0) {
          r_vv_t[i] = r_vv[i];
        }
        else {
          // (x,y,z) -> (r,r0*theta,z)
	  coordinateTransform(x_vv_t[i], x_vv[i]);
          // for r less than wid.r0, the space will be stretched in the
          // azimuthal direction. Increase r_vv proportionally to avoid missing pixels,
          // up to a maximum factor of 2 (otherwise gets expensive)...
          if (x_vv_t[i][0] > 0.5*wid.r0) {
            r_vv_t[i] = r_vv[i]*max(wid.r0/x_vv_t[i][0],1.0);
          }
          else {
            r_vv_t[i] = r_vv[i]*2.0;
          }
        }
      }
      canvas->addPlaneDataRvvPositive(phi,x_vv_t,r_vv_t,n,factor);
      delete[] x_vv_t;
      delete[] r_vv_t;
    }
    else {
      canvas->addPlaneDataRvvPositive(phi,x_vv,r_vv,n,factor);
    }
  }

  void addPartialSurfaceData(const double * phi,const double (*x_vv)[3],const double *r_vv,const int * zone_vv,const int n,const double factor) {
    // factor is a number that multiples the passed r_vv to determine the max sphere of potential influence of the
    // the vd's...
    if (wid.b_transform_cyl_x) {
      double (*x_vv_t)[3] = new double[n][3];
      double *r_vv_t = new double[n];
      for (int i = 0; i < n; ++i) {
	// (x,y,z) -> (x,r,r0*theta)
	coordinateTransform(x_vv_t[i], x_vv[i]);
        // for r less than wid.r0, the space will be stretched in the
        // azimuthal direction. Increase r_vv proportionally to avoid missing pixels,
        // up to a maximum factor of 2 (otherwise gets expensive)...
        if (x_vv_t[i][1] > 0.5*wid.r0) {
          r_vv_t[i] = r_vv[i]*max(wid.r0/x_vv_t[i][1],1.0);
        }
        else {
          r_vv_t[i] = r_vv[i]*2.0;
        }
      }
      canvas->addPartialSurfaceData(phi,x_vv_t,r_vv_t,zone_vv,n,factor);
      delete[] x_vv_t;
      delete[] r_vv_t;
    }
    else if (wid.b_transform_cyl_z) {
      double (*x_vv_t)[3] = new double[n][3];
      double *r_vv_t = new double[n];
      for (int i = 0; i < n; ++i) {
	// (x,y,z) -> (r,r0*theta,z)
	coordinateTransform(x_vv_t[i], x_vv[i]);
        // for r less than wid.r0, the space will be stretched in the
        // azimuthal direction. Increase r_vv proportionally to avoid missing pixels,
        // up to a maximum factor of 2 (otherwise gets expensive)...
        if (x_vv_t[i][0] > 0.5*wid.r0) {
          r_vv_t[i] = r_vv[i]*max(wid.r0/x_vv_t[i][0],1.0);
        }
        else {
          r_vv_t[i] = r_vv[i]*2.0;
        }
      }
      canvas->addPartialSurfaceData(phi,x_vv_t,r_vv_t,zone_vv,n,factor);
      delete[] x_vv_t;
      delete[] r_vv_t;
    }
    else {
      canvas->addPartialSurfaceData(phi,x_vv,r_vv,zone_vv,n,factor);
    }

  }

  void addSurfaceData(const double * phi,const double (*x_vv)[3],const double *r_vv,const int n,const double factor) {
    // factor is a number that multiples the passed r_vv to determine the max sphere of potential influence of the
    // the vd's...
    if (wid.b_transform_cyl_x) {
      double (*x_vv_t)[3] = new double[n][3];
      double *r_vv_t = new double[n];
      for (int i = 0; i < n; ++i) {
	// (x,y,z) -> (x,r,r0*theta)
	coordinateTransform(x_vv_t[i], x_vv[i]);
        // for r less than wid.r0, the space will be stretched in the
        // azimuthal direction. Increase r_vv proportionally to avoid missing pixels,
        // up to a maximum factor of 2 (otherwise gets expensive)...
        if (x_vv_t[i][1] > 0.5*wid.r0) {
          r_vv_t[i] = r_vv[i]*max(wid.r0/x_vv_t[i][1],1.0);
        }
        else {
          r_vv_t[i] = r_vv[i]*2.0;
        }
      }
      canvas->addSurfaceData(phi,x_vv_t,r_vv_t,n,factor);
      delete[] x_vv_t;
      delete[] r_vv_t;
    }
    else if (wid.b_transform_cyl_z) {
      double (*x_vv_t)[3] = new double[n][3];
      double *r_vv_t = new double[n];
      for (int i = 0; i < n; ++i) {
	// (x,y,z) -> (r,r0*theta,z)
	coordinateTransform(x_vv_t[i], x_vv[i]);
        // for r less than wid.r0, the space will be stretched in the
        // azimuthal direction. Increase r_vv proportionally to avoid missing pixels,
        // up to a maximum factor of 2 (otherwise gets expensive)...
        if (x_vv_t[i][0] > 0.5*wid.r0) {
          r_vv_t[i] = r_vv[i]*max(wid.r0/x_vv_t[i][0],1.0);
        }
        else {
          r_vv_t[i] = r_vv[i]*2.0;
        }
      }
      canvas->addSurfaceData(phi,x_vv_t,r_vv_t,n,factor);
      delete[] x_vv_t;
      delete[] r_vv_t;
    }
    else {
      canvas->addSurfaceData(phi,x_vv,r_vv,n,factor);
    }

  }

  void rmPbdsWithoutSurfaceData() {
    canvas->rmPbdsWithoutSurfaceData();
  }

  void addSurfaceMeshRvvPositive(const double (*x_vv)[3],const double *r_vv,const int *i_vv,const int *ivv_global,const int n,const double factor = 1.0) {
    // i_vv is used to store stuff like level and window index?
    canvas->addSurfaceMeshRvvPositive(x_vv,r_vv,i_vv,ivv_global,n,factor);
  }

  void addPlaneMeshRvvPositive(const double (*x_vv)[3],const double *r_vv,const int *i_vv,const int *ivv_global,const int n,const double factor = 1.0) {
    // i_vv is used to store stuff like level and window index
    if (wid.b_transform_cyl_x) {
      double (*x_vv_t)[3] = new double[n][3];
      double *r_vv_t = new double[n];
      for (int i = 0; i < n; ++i) {
        if (r_vv[i] <= 0.0) {
          r_vv_t[i] = r_vv[i];
        }
        else {
          // (x,y,z) -> (x,r,r0*theta)
	  coordinateTransform(x_vv_t[i], x_vv[i]);
          // for r less than wid.r0, the space will be stretched in the
          // azimuthal direction. Increase r_vv proportionally to avoid missing pixels,
          // up to a maximum factor of 2 (otherwise gets expensive)...
          if (x_vv_t[i][1] > 0.5*wid.r0) {
            r_vv_t[i] = r_vv[i]*max(wid.r0/x_vv_t[i][1],1.0);
          }
          else {
            r_vv_t[i] = r_vv[i]*2.0;
          }
        }
      }
      canvas->addPlaneMeshRvvPositive(x_vv_t,r_vv_t,i_vv,ivv_global,n,factor);
      delete[] x_vv_t;
      delete[] r_vv_t;
    }
    else if (wid.b_transform_cyl_z) {
      double (*x_vv_t)[3] = new double[n][3];
      double *r_vv_t = new double[n];
      for (int i = 0; i < n; ++i) {
        if (r_vv[i] <= 0.0) {
          r_vv_t[i] = r_vv[i];
        }
        else {
          // (x,y,z) -> (r,r0*theta,z)
	  coordinateTransform(x_vv_t[i], x_vv[i]);
          // for r less than wid.r0, the space will be stretched in the
          // azimuthal direction. Increase r_vv proportionally to avoid missing pixels,
          // up to a maximum factor of 2 (otherwise gets expensive)...
          if (x_vv_t[i][0] > 0.5*wid.r0) {
            r_vv_t[i] = r_vv[i]*max(wid.r0/x_vv_t[i][0],1.0);
          }
          else {
            r_vv_t[i] = r_vv[i]*2.0;
          }
        }
      }
      canvas->addPlaneMeshRvvPositive(x_vv_t,r_vv_t,i_vv,ivv_global,n,factor);
      delete[] x_vv_t;
      delete[] r_vv_t;
    }
    else {
      canvas->addPlaneMeshRvvPositive(x_vv,r_vv,i_vv,ivv_global,n,factor);
    }
  }

  void addPointsMeshRvvPositive(const double (*x_vv)[3],const double *r_vv,const int *i_vv,const int n,const double factor = 1.0) {
    // i_vv is used to store stuff like level and window index
    canvas->addPointsMeshRvvPositive(x_vv,r_vv,i_vv,n,factor);
  }

  void addSimpleInternalTrisWithMeshOnEdge(const vector<pair<SimpleTri,int> >& triVec) {
    assert(canvas!=NULL);

    // add the internal tris...
    for (int it = 0, limit = triVec.size(); it < limit; ++it) {
      uchar mesh = 0;
      if (wid.showSurfaceMesh && (triVec[it].second != -1)) mesh=triVec[it].second;

      //    if (triVec[it].second < 3) mesh |= (1 << triVec[it].second);
      //    else mesh=7;  // draw all edges (tet mesh portions)
      // }

      if (wid.b_transform_cyl_x || wid.b_transform_cyl_z) {
	double x[3][3];
	FOR_I3 x[0][i] = triVec[it].first.x0[i];
	FOR_I3 x[1][i] = triVec[it].first.x1[i];
	FOR_I3 x[2][i] = triVec[it].first.x2[i];
	coordinateTransform(x[0], triVec[it].first.x0);
	coordinateTransform(x[1], triVec[it].first.x1);
	coordinateTransform(x[2], triVec[it].first.x2);
	adjustThetaTri(x[0],x[1],x[2]);
	canvas->addInternalTri(x[0],x[1],x[2],32767,mesh);
      }
      else {
	canvas->addInternalTri(triVec[it].first.x0,triVec[it].first.x1,triVec[it].first.x2,32767,mesh);
      }
    }
  }

  void addSimpleInternalTrisWithMeshOnFirstEdge(const vector<SimpleTri>& triVec,const bool on_data_plane = false) {
    addSimpleInternalTrisWithMesh(triVec,1,on_data_plane);
  }

  void addSimpleInternalTrisWithMesh(const vector<SimpleTri>& triVec,const int _mesh = 7,const bool on_data_plane = false,const int zone=32767) {
    // currently defaults to protected internal face zone
    assert(canvas!=NULL);

    const uchar mesh = (wid.showSurfaceMesh) ? _mesh:0;

    // add the internal tris...
    for (int it = 0, limit = triVec.size(); it < limit; ++it) {
      if (wid.b_transform_cyl_x || wid.b_transform_cyl_z) {
	double x[3][3];
	FOR_I3 x[0][i] = triVec[it].x0[i];
	FOR_I3 x[1][i] = triVec[it].x1[i];
	FOR_I3 x[2][i] = triVec[it].x2[i];
	coordinateTransform(x[0], triVec[it].x0);
	coordinateTransform(x[1], triVec[it].x1);
	coordinateTransform(x[2], triVec[it].x2);
	adjustThetaTri(x[0],x[1],x[2]);
	canvas->addInternalTri(x[0],x[1],x[2],zone,mesh,on_data_plane);
      }
      else {
	canvas->addInternalTri(triVec[it].x0,triVec[it].x1,triVec[it].x2,zone,mesh,on_data_plane);
      }
    }
  }

  void addSimpleSurfaceTrisWithMesh(const vector<pair<SimpleTriWithData,int> >& triVec,const int _mesh = 7,const bool on_data_plane = false) {
    // tri data includes the normal to use
    assert(canvas != NULL);

    // the set is sorted, so we want the last entry. Set has no back(), so
    // do this with the reverse iterator...

    // in Stitch both subzones and zones are redundant
    // to comprehend zone-name based lsitings, we nee dto process both zones and subzones
    int hidden_isz_max = -1;
    if (!wid.hiddenSubzonesSet.empty())
      hidden_isz_max = *wid.hiddenSubzonesSet.rbegin();

    if (!wid.hiddenZonesSet.empty())
      hidden_isz_max = max(hidden_isz_max,*wid.hiddenZonesSet.rbegin());

    assert(hidden_isz_max < 65535);
    int * hidden_isz_flag = new int[hidden_isz_max+1];
    for (int isz = 0; isz <= hidden_isz_max; ++isz) hidden_isz_flag[isz] = 0;

    for (set<int>::iterator it = wid.hiddenSubzonesSet.begin(); it != wid.hiddenSubzonesSet.end(); ++it) {
      const int isz = *it;
      if ((isz >= 0)&&(isz <= hidden_isz_max)) hidden_isz_flag[isz] = 1;
      else {
        if (mpi_rank == 0) cout << "hidden subzone index: " << isz << " is out of range; ignoring" << endl;
      }
    }

    for (set<int>::iterator it = wid.hiddenZonesSet.begin(); it != wid.hiddenZonesSet.end(); ++it) {
      const int isz = *it;
      if ((isz >= 0)&&(isz <= hidden_isz_max)) hidden_isz_flag[isz] = 1;
      else {
        if (mpi_rank == 0) cout << "hidden zone index: " << isz << " is out of range; ignoring" << endl;
      }
    }

    // add the surface tris...
    const uchar mesh = (wid.showSurfaceMesh) ? _mesh:0;

    for (int it = 0, limit = triVec.size(); it < limit; ++it) {
      const int index = triVec[it].second;
      assert((index >= 0)&&(index < 32768));
      const bool hide = !((index > hidden_isz_max)||(hidden_isz_flag[index] == 0));

      const double normal[3] = {triVec[it].first.d0,triVec[it].first.d1,triVec[it].first.d2};
      canvas->addSurfaceTri(triVec[it].first.x0,triVec[it].first.x1,triVec[it].first.x2,index,hide,mesh,on_data_plane,&normal[0]);
    }

    delete[] hidden_isz_flag;
  }

  void addSimpleSurfaceTrisWithMesh(const vector<pair<SimpleTri,int> >& triVec,const int _mesh = 7,const bool on_data_plane = false) {

    assert(canvas != NULL);

    // the set is sorted, so we want the last entry. Set has no back(), so
    // do this with the reverse iterator...

    int hidden_isz_max = -1;
    if (!wid.hiddenSubzonesSet.empty())
      hidden_isz_max = *wid.hiddenSubzonesSet.rbegin();

    if (!wid.hiddenZonesSet.empty())
      hidden_isz_max = max(hidden_isz_max,*wid.hiddenZonesSet.rbegin());

    assert(hidden_isz_max < 65535);
    int * hidden_isz_flag = new int[hidden_isz_max+1];
    for (int isz = 0; isz <= hidden_isz_max; ++isz) hidden_isz_flag[isz] = 0;

    for (set<int>::iterator it = wid.hiddenSubzonesSet.begin(); it != wid.hiddenSubzonesSet.end(); ++it) {
      const int isz = *it;
      if ((isz >= 0)&&(isz <= hidden_isz_max)) hidden_isz_flag[isz] = 1;
      else {
        if (mpi_rank == 0) cout << "hidden subzone index: " << isz << " is out of range; ignoring" << endl;
      }
    }

    for (set<int>::iterator it = wid.hiddenZonesSet.begin(); it != wid.hiddenZonesSet.end(); ++it) {
      const int isz = *it;
      if ((isz >= 0)&&(isz <= hidden_isz_max)) hidden_isz_flag[isz] = 1;
      else {
        if (mpi_rank == 0) cout << "hidden zone index: " << isz << " is out of range; ignoring" << endl;
      }
    }

    // add the surface tris...
    const uchar mesh = (wid.showSurfaceMesh) ? _mesh:0;

    for (int it = 0, limit = triVec.size(); it < limit; ++it) {
      const int index = triVec[it].second;
      assert((index >= 0)&&(index < 32768));
      const bool hide = !((index > hidden_isz_max)||(hidden_isz_flag[index] == 0));
      canvas->addSurfaceTri(triVec[it].first.x0,triVec[it].first.x1,triVec[it].first.x2,index,hide,mesh,on_data_plane);
    }

    delete[] hidden_isz_flag;
  }

  // vector vizualation
  void addInternalVectors(const double (*x0)[3],const double (*dx)[3],const int index,const int n,const double scale) {
    assert(canvas != NULL);
    for (int i = 0; i < n; ++i) canvas->addVector(ISOSURFACE_VECTOR,x0[i],dx[i],scale,index);
  }
  void addSurfaceVectors(const double (*x0)[3],const double (*dx)[3],const int *zones,const int n,const double scale) {
    assert(canvas != NULL);
    for (int i = 0; i < n; ++i) canvas->addVector(SURFACE_VECTOR,x0[i],dx[i],scale,zones[i]);
  }

  void addSimpleInternalTrisWithData(const vector<pair<SimpleTriWithData,int> >& triVec, const int index) {
    assert(canvas != NULL);

    if (index == 65532) {
      // isosurface data...
      float * range = NULL;
      if (wid.b_clip_iso) range = &wid.range_on_iso[0];

      for (int it = 0, limit = triVec.size(); it < limit; ++it) {
        const double data[3] = {triVec[it].first.d0,triVec[it].first.d1,triVec[it].first.d2};

        uchar mesh = 0;
        if (wid.showSurfaceMesh && (triVec[it].second != -1)) mesh=triVec[it].second;
        //    if (triVec[it].second < 3) mesh |= (1 << triVec[it].second);
        //    else mesh=7;  // draw all edges (tet mesh portions)
        // }

        //uchar mesh = 0;
        //if (wid.showSurfaceMesh && (triVec[it].second != -1)) mesh |= (1 << triVec[it].second);

	if (wid.b_transform_cyl_x || wid.b_transform_cyl_z) {

	  double x_no0[3][3], x_no_t[3][3];
	  FOR_I3 x_no0[0][i] = triVec[it].first.x0[i];
	  FOR_I3 x_no0[1][i] = triVec[it].first.x1[i];
	  FOR_I3 x_no0[2][i] = triVec[it].first.x2[i];
	  for (int ino = 0; ino < 3; ++ino) {
	    if (ino > 0)
	      coordinateTransformAndAdjustTheta(x_no_t[ino], x_no0[ino], x_no_t[0]);
	    else
	      coordinateTransform(x_no_t[ino], x_no0[ino]);
	  }

	  canvas->addInternalTri(x_no_t[0],x_no_t[1],x_no_t[2],index,mesh,true,&data[0]);
	}
	else {
	  canvas->addInternalTri(triVec[it].first.x0,triVec[it].first.x1,triVec[it].first.x2,index,mesh,false,&data[0],false,range);
	}
      }
    }
    else if (index == 65534) {
      // plane data...
      for (int it = 0, limit = triVec.size(); it < limit; ++it) {
        const double data[3] = {triVec[it].first.d0,triVec[it].first.d1,triVec[it].first.d2};

        uchar mesh = 0;
        if (wid.showSurfaceMesh && (triVec[it].second != -1)) mesh=triVec[it].second;
        // mesh |= (1 << triVec[it].second);

	if (wid.b_transform_cyl_x || wid.b_transform_cyl_z) {

	  double x_no0[3][3], x_no_t[3][3];
	  FOR_I3 x_no0[0][i] = triVec[it].first.x0[i];
	  FOR_I3 x_no0[1][i] = triVec[it].first.x1[i];
	  FOR_I3 x_no0[2][i] = triVec[it].first.x2[i];
	  for (int ino = 0; ino < 3; ++ino) {
	    if (ino > 0)
	      coordinateTransformAndAdjustTheta(x_no_t[ino], x_no0[ino], x_no_t[0]);
	    else
	      coordinateTransform(x_no_t[ino], x_no0[ino]);
	  }

	  canvas->addInternalTri(x_no_t[0],x_no_t[1],x_no_t[2],index,mesh,true,&data[0]);
	}
	else {
	  canvas->addInternalTri(triVec[it].first.x0,triVec[it].first.x1,triVec[it].first.x2,index,mesh,true,&data[0]);
	}
      }
    }
    else {
      // unsupported data...
      assert(0);
    }
  }

  void addParticle(const double xp[3],const double dp) {
    canvas->addParticle(xp,dp,wid.b_no_particle_blanking);
  }

  void addParticle(const double xp[3],const double dp,const double vp) {
    canvas->addParticle(xp,dp,vp,wid.b_no_particle_blanking);
  }

  void addParticles(const int np,const double (*xp)[3],const double * dp,const double factor = 1.0) {
    for (int ip = 0; ip < np; ++ip)
      canvas->addParticle(xp[ip],dp[ip]*factor,wid.b_no_particle_blanking);
  }

  void addParticles(const int np,const double (*xp)[3],const double * dp,const double * vp,const double factor = 1.0) {
    for (int ip = 0; ip < np; ++ip)
      canvas->addParticle(xp[ip],dp[ip]*factor,vp[ip],wid.b_no_particle_blanking);
  }

  void prepVolvis() {

    canvas->setVolVis();
    canvas->setVolvisType(wid.volvis_type);
    if (wid.b_volvis_aux_data) canvas->setVolvisAuxData(wid.volvis_aux_data[0],wid.volvis_aux_data[1]);
    if (wid.b_volvis_surface) canvas->setVolvisSurface();
  }

  void addSurfaceTrisVolVis(const double (*xsp)[3],const int (*spost)[3],const int * znost,const int nst) {

    assert(canvas);

    const int nzn = zoneNames.size();
    IntFlag show_zn_flag(nzn);
    if (isSubzoneHidden(-1)) {
      // hide all flag
      show_zn_flag.setAll(0);
    }
    else {
      show_zn_flag.setAll(1);
      if (wid.b_hidesubzones || wid.b_hidezones) {
        set<int>::iterator it;
        for (set<int>::iterator it = wid.hiddenSubzonesSet.begin(); it != wid.hiddenSubzonesSet.end(); ++it) {
          const int iszn = *it;
          if ((iszn >= 0)&&(iszn < nzn)) show_zn_flag[iszn] = 0;
          else {
            if (mpi_rank == 0) cout << "hidden subzone index: " << iszn << " is out of range; ignoring" << endl;
          }
        }
        for (set<int>::iterator it = wid.hiddenZonesSet.begin(); it != wid.hiddenZonesSet.end(); ++it) {
          const int izn = *it;
          if ((izn >= 0)&&(izn < nzn)) show_zn_flag[izn] = 0;
          else {
            if (mpi_rank == 0) cout << "hidden zone index: " << izn << " is out of range; ignoring" << endl;
          }
        }
      }
    }

    for (int ist = 0; ist < nst; ++ist) {
      const int izn = znost[ist];
      assert((izn >= 0)&&(izn < nzn));
      // for VolVis, we add EVERY tri (even periodic, izn >= nzn)..
      canvas->addSurfaceTriVolVis(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],show_zn_flag[izn]==1); // arg is show
    }

    // if requested, include the surface...
    if (wid.b_volvis_surface) {
      const uchar mesh = (wid.showSurfaceMesh) ? 7:0;  // render all three edges

      for (int ist = 0; ist < nst; ++ist) {
        const int izn = znost[ist];
        canvas->addSurfaceTri(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],znost[ist],show_zn_flag[izn]==0,mesh);
      }
    }

  }

  void addPeriodicFaceVolVis(const double (*x_no)[3],const double x_fa[3],const int nnof,const int *noofa_v) {
    assert(canvas);

    // triangulate voronoi face
    int ino1 = noofa_v[nnof-1];
    for (int nof = 0; nof < nnof; ++nof) {
      const int ino0 = ino1;
      ino1 = noofa_v[nof];
      canvas->addSurfaceTriVolVis(x_no[ino0],x_no[ino1],x_fa,false); // periodic faces are always no show (user has no control of them)
    }

  }

  void addBoundaryFacesVolVis(const double (*x_no)[3],const double (*x_bf)[3],const double (*n_bf)[3],const int *zone_bf,const int *noobf_i,const int *noobf_v,const int nbf) {

    assert(canvas);

    // triangulate voronoi faces
    // only add shown tris for path integration start/stops
    FOR_IBF {
      int ino1 = noobf_v[noobf_i[ibf+1]-1];
      for (int nof = noobf_i[ibf]; nof != noobf_i[ibf+1]; ++nof) {
        const int ino0 = ino1;
        ino1 = noobf_v[nof];
        // add if not hidden
        canvas->addSurfaceTriVolVis(x_no[ino0],x_no[ino1],x_bf[ibf],!bf_hide[ibf],n_bf[ibf][0],n_bf[ibf][1],n_bf[ibf][2]); // periodic faces are always no show (user has no control of them)
      }
    }

    // if requested, include the surface rendering too...
    if (wid.b_volvis_surface) {
      const uchar mesh = (wid.showSurfaceMesh) ? 1:0;  // draw edges on first edge

      FOR_IBF {
        int ino1 = noobf_v[noobf_i[ibf+1]-1];
        for (int nof = noobf_i[ibf]; nof != noobf_i[ibf+1]; ++nof) {
          const int ino0 = ino1;
          ino1 = noobf_v[nof];
          canvas->addSurfaceTri(x_no[ino0],x_no[ino1],x_bf[ibf],zone_bf[ibf],bf_hide[ibf],mesh,false,&n_bf[ibf][0]);
        }
      }
    }
  }

  void addVoronoiPointVolVis(const double x_vv[3],const double r_vv,const float v) {

    assert(canvas);
    canvas->addVoronoiPointVolVis(x_vv,r_vv,v);

  }

  //this routine will build or read a CvImageMap for the requested data plane
  //if a mesh image is requested it will build additional data if the CvImageMap was
  //read from disk.
  void initDataPlaneCim(CvPlaneMap& cvPlaneMap) {
    assert(canvas!=NULL);
    if (wid.b_geom_plane&&wid.b_var&&(filename_volumeData!="")){

      if (!lim){  //if also visualizing surface data, the lim will already be allocated
        if (filename_volumeData_2 == "")
          lim = new LesImageMapper(filename_volumeData);
        else
          lim = new LesImageMapper(filename_volumeData,filename_volumeData_2);
      }

      //Read cim from disk if available, filename construction does not require canvas mask be built
      string cimName = lim->buildCvImageMapFilename(canvas,wid.xp,wid.np);
      //cout << " > Unique CvImageMap file name: " << cimName << endl;
      int readCim = cim.read(cimName.c_str());
      if (readCim<0){
        //cim is not available from disk, we must construct a mask and build a cim
        cvPlaneMap.init(wid.xp,wid.np); //will release old cv data if new xp,np
        lim->buildCvImageMap(cim,canvas,cvPlaneMap,lengthscaleName);
        cim.write(cimName.c_str());  //write the cim in case this view is ever requested again
      }
      /*
      // required call for DAP mesh viz. Currently using stitch mesh viz.
      if (wid.var=="MESH"||wid.var=="mesh"){
      //if requesting a mesh view, set the mesh_buf pointer in lim.
      //if the cim was read, this will be computed prior to being set.
      //if the cim was just computed, d2_buf is reassigned to mesh_buf.
      //lim can continue to be used to compute cims for surface data.
      lim->computeCvPixelD2(cim, canvas, wid.xp, wid.np);
      }
      */
    }
  }

  void addEdges(const double (*xsp)[3],const int (*spost)[3],const int (*teost)[3],const int *znost,const IntFlag &szost,const IntFlag &szozn_i,const int nst) {
    assert(canvas!=NULL);

    // flag stores which zones should be hidden (i.e., not added)
    const int nzn = szozn_i.getLength() - 1;
    const int nsz = szozn_i[nzn];
    IntFlag show_sz_flag(nsz);
    show_sz_flag.setAll(1);
    if (wid.b_hidesubzones || wid.b_hidezones) {
      set<int>::iterator it;
      for (it = wid.hiddenZonesSet.begin(); it != wid.hiddenZonesSet.end(); ++it) {
        if (*it == -1) {
          assert(wid.hiddenZonesSet.size() == 1);
          for (int isz = 0; isz < nsz; ++isz) {
            show_sz_flag[isz] = 0;
          }
        }
        else {
          const int izn = *it;

          if ((izn >= 0)&&(izn < nzn)) {
            for (int iszn = szozn_i[izn]; iszn < szozn_i[izn+1]; ++iszn) show_sz_flag[iszn] = 0;
          }
          else {
            if (mpi_rank == 0) cout << "WARNING: izn out of bounds: " << izn << " (nzn=" << nzn << ") File: " << __FILE__ << ", line: " << __LINE__ << endl;
          }
        }
      }
      for (it = wid.hiddenSubzonesSet.begin(); it != wid.hiddenSubzonesSet.end(); ++it) {
        const int iszn = *it;
        if ((iszn >= 0)&&(iszn < nsz)) {
          show_sz_flag[iszn] = 0;
        }
        else {
          if (mpi_rank == 0) cout << "WARNING: iszn out of bounds: " << iszn << " (nsz=" << nsz << ") File: " << __FILE__ << ", line: " << __LINE__ << endl;
        }
      }
    }

    // add edges if requested...
    if (wid.showOpenEdges) {
      canvas->setShowEdges(true);
      if (wid.showSurfaceMesh) {
        for (int ist = 0; ist < nst; ++ist) {
          const int index = szost[ist];  // subzone index
          if (index<nsz&&index>=0&&show_sz_flag[index]) {
            for (int i = 0; i < 3; ++i) {
              if (teost[ist][i] != 0) {
                const int isp0 = spost[ist][i];
                const int isp1 = spost[ist][(i+1)%3];
                const double dx01[3] = DIFF(xsp[isp1],xsp[isp0]);
                const double dx02[3] = DIFF(xsp[spost[ist][(i+2)%3]],xsp[isp0]);
                const double l01_sq = DOT_PRODUCT(dx01,dx01);
                if (l01_sq > 0.0) {
                  const double proj = DOT_PRODUCT(dx01,dx02)/l01_sq;
                  // dx here is the highlighted edge displacement relative to the actual edge location
                  double dx[3]; // edge unit normal in dir of i+2
                  FOR_J3 dx[j] = dx02[j]-proj*dx01[j];
                  const double l = MAG(dx);

                  if (l <= 0.0) FOR_J3 dx[j] = 0.0;  // linear tri; don't displace b/c proper direction is undefined
                  else FOR_J3 dx[j] /= l;

                  if (teost[ist][i] == -1) {
                    FOR_J3 dx[j] = 0.0;
                    canvas->addEdge(xsp[isp0],xsp[isp1],dx,32768,true);
                  }
                  else if (teost[ist][i] < -1) {
                    canvas->addEdge(xsp[isp0],xsp[isp1],dx,-teost[ist][i]-1+32768,true);
                  }
                  else {
                    canvas->addEdge(xsp[isp0],xsp[isp1],dx,teost[ist][i]-1+49152,true);
                  }
                }
              }
            }
          }
        }
      }
      else {
        for (int ist = 0; ist < nst; ++ist) {
          const int index = szost[ist];  // subzone index
          if (index<nsz&&index>=0&&show_sz_flag[index]) {
            for (int i = 0; i < 3; ++i) {
              if (teost[ist][i] != 0) {
                const int isp0 = spost[ist][i];
                const int isp1 = spost[ist][(i+1)%3];
                const double dx01[3] = DIFF(xsp[isp1],xsp[isp0]);
                const double dx02[3] = DIFF(xsp[spost[ist][(i+2)%3]],xsp[isp0]);
                const double l01_sq = DOT_PRODUCT(dx01,dx01);
                if (l01_sq > 0.0) {
                  const double proj = DOT_PRODUCT(dx01,dx02)/l01_sq;
                  // dx here is the highlighted edge displacement relative to the actual edge location
                  double dx[3]; // edge unit normal in dir of i+2
                  FOR_J3 dx[j] = dx02[j]-proj*dx01[j];
                  const double l = MAG(dx);

                  if (l <= 0.0) FOR_J3 dx[j] = 0.0;  // linear tri; don't dislpace b/c proper direction is undefined
                  else FOR_J3 dx[j] /= l;

                  if (teost[ist][i] == -1) {
                    FOR_J3 dx[j] = 0.0;
                    canvas->addEdge(xsp[isp0],xsp[isp1],dx,32768,false);  // protected for misaligned nbr edges
                  }
                  else if (teost[ist][i] < -1) {
                    canvas->addEdge(xsp[isp0],xsp[isp1],dx,-teost[ist][i]-1+32768,false);
                  }
                  else {
                    canvas->addEdge(xsp[isp0],xsp[isp1],dx,teost[ist][i]-1+49152,false);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // void writeImage() {
  //   if (wid.b_range){
  //     canvas->setDataRange(wid.range[0],wid.range[1]);
  //   }
  //   if (wid.b_range_on_surface){
  //     canvas->setSurfaceDataRange(wid.range_on_surface[0],wid.range_on_surface[1]);
  //   }
  //   if (wid.b_range_on_particle){
  //     canvas->setParticleDataRange(wid.range_on_particle[0],wid.range_on_particle[1]);
  //   }
  //   if (wid.b_range_on_iso){
  //     canvas->setIsoDataRange(wid.range_on_iso[0],wid.range_on_iso[1]);
  //   }
  //   if (wid.range[2] != 0.f){
  //     canvas->setDataRangeBins((int)wid.range[2]));
  //   }
  //   if (wid.range_on_surface[2] != 0.f){
  //     canvas->setSurfaceDataRangeBins((int)wid.range_on_surface[2]));
  //   }
  //   if (wid.range_on_iso[2] != 0.f){
  //     canvas->setIsoDataRangeBins((int)wid.range_on_iso[2]));
  //   }
  //   canvas->writeImage(wid.name+".png",wid.rgb_mode);
  // }

  void writeImage(const int index=-1) {
    char filename[128];
    // add "app" param parsing for compatibility with new app interactions...
    if (wid.b_app) {
      sprintf(filename,"%s",wid.app.c_str());
    }
    else if (index != -1) {
      sprintf(filename,"%s.%08d.png",wid.name.c_str(),index);
    }
    else {
      sprintf(filename,"%s.png",wid.name.c_str());
    }

    if (wid.b_range){
      canvas->setDataRange(wid.range[0],wid.range[1]);
    }
    if (wid.b_range_on_surface){
      canvas->setSurfaceDataRange(wid.range_on_surface[0],wid.range_on_surface[1]);
    }
    if (wid.b_range_on_particle){
      canvas->setParticleDataRange(wid.range_on_particle[0],wid.range_on_particle[1]);
    }
    if (wid.b_range_on_iso){
      canvas->setIsoDataRange(wid.range_on_iso[0],wid.range_on_iso[1]);
    }
    if (wid.range[2] != 0.f){
      canvas->setDataRangeBins((int)wid.range[2]);
    }
    if (wid.range_on_surface[2] != 0.f){
      canvas->setSurfaceDataRangeBins((int)wid.range_on_surface[2]);
    }
    if (wid.range_on_iso[2] != 0.f){
      canvas->setIsoDataRangeBins((int)wid.range_on_iso[2]);
    }
    canvas->writeImage(filename,wid.rgb_mode);
  }

  void setImdTime(const double time = 0.0) {
    imd->setTime(time);
  }

  void addCimDataAndWriteImage(){
    assert(canvas!=NULL);

    if (wid.b_snapshot){
      writeImagesFromSnapshotData();
    }
    else{ //no snapshots, just a single image
      //Add serial surface data if present
      sparseReadAndAddSurfaceData();
      //Add serial volume data if present
      sparseReadAndAddVolumeData();
      this->writeImage();
    }
  }

  void writeImagesFromSnapshotData(){
    bool first = true;

    if (wid.b_range){
      canvas->setDataRange(wid.range[0],wid.range[1]);
    }
    if (wid.b_range_on_surface){
      canvas->setSurfaceDataRange(wid.range_on_surface[0],wid.range_on_surface[1]);
    }
    if (wid.range[2] != 0.f){
      canvas->setDataRangeBins((int)wid.range[2]);
    }
    if (wid.range_on_surface[2] != 0.f){
      canvas->setSurfaceDataRangeBins((int)wid.range_on_surface[2]);
    }

    for (int iS=wid.snapshot_f; iS<=wid.snapshot_l; iS+=wid.snapshot_inc) {

      // always clear current data when user requests for a new snapshot (because we have NEW current data)...
      CtiRegister::clearCurrentData();

      std::stringstream ss;
      ss << wid.snapshot_name << "." << std::setw(8) << std::setfill('0') << iS << ".sles";
      LesImageMapper limSnapshot(filename_volumeData,ss.str());
      current_lim_ptr = &limSnapshot;

      if (first) {
        if (cimS.ncv_unique > 0) {
          current_cim_ptr = &cimS;
          CtiRegister::CtiData * data = CtiRegister::getCtiData(wid.var_on_surface);
          if ((data == NULL)||(data->getUnindexedTopology() != CV_DATA)||(data->getType() != DN_DATA)) {
            COUT1(" > Warning: requested expression " << wid.var << " does not eval to CVDN data. Skipping.");
          }
          else {
            canvas->addCimSurfaceData(cimS,data->getDNptr());
          }
        }
        if (cim.ncv_unique > 0) {
          current_cim_ptr = &cim;
          CtiRegister::CtiData * data = CtiRegister::getCtiData(wid.var);
          if ((data == NULL)||(data->getUnindexedTopology() != CV_DATA)||(data->getType() != DN_DATA)) {
            COUT1(" > Warning: requested expression " << wid.var << " does not eval to CVDN data. Skipping.");
          }
          else {
            canvas->addCimPlaneData(cim,data->getDNptr());
          }
        }
        canvas->cacheImage();
        first = false;
      }
      else {
        if (cimS.ncv_unique > 0) {
          current_cim_ptr = &cimS;
          CtiRegister::CtiData * data = CtiRegister::getCtiData(wid.var_on_surface);
          if ((data == NULL)||(data->getUnindexedTopology() != CV_DATA)||(data->getType() != DN_DATA)) {
            COUT1(" > Warning: requested expression " << wid.var << " does not eval to CVDN data. Skipping.");
          }
          else {
            canvas->updateCimSurfaceData(cimS,data->getDNptr());
          }
        }
        if (cim.ncv_unique > 0) {
          current_cim_ptr = &cim;
          CtiRegister::CtiData * data = CtiRegister::getCtiData(wid.var);
          if ((data == NULL)||(data->getUnindexedTopology() != CV_DATA)||(data->getType() != DN_DATA)) {
            COUT1(" > Warning: requested expression " << wid.var << " does not eval to CVDN data. Skipping.");
          }
          else {
            canvas->updateCimPlaneData(cim,data->getDNptr());
          }
        }
      }

      ss.str(std::string());
      ss.clear();
      ss << wid.name << "." << std::setw(8) << std::setfill('0') << iS << ".png";
      imd->setTime(limSnapshot.getTime());  //update simulation time in metadata
      canvas->flushImage(ss.str(),wid.rgb_mode);
    } //for all snapshots...
    canvas->clearPBuf();
  }

  void sparseReadAndAddSurfaceData() {
    //if (wid.b_var_on_surface && (surfaceData==NULL || !(wid.var_on_surface=="G"||wid.var_on_surface=="G_ST"||wid.var_on_surface=="G_SP") ) )
    if (cimS.ncv_unique>0){

      // need to make sure cimS is current...
      current_cim_ptr = &cimS;
      current_lim_ptr = lim;

      if (wid.var_on_surface == "LIC") {
        //Surface LIC should not be functional - DP
        double (*vol_xcc)[3] = new double[cimS.ncv_unique][3];

        int8 data_offset_xcc = lim->getRfpCvD2Offset("x_cv");
        if (data_offset_xcc!=-1){
          lim->sparseReadCvR2(vol_xcc,cimS,"x_cv",data_offset_xcc);
        }
        else if (filename_volumeData_2 != ""){
          data_offset_xcc = lim->getSfpCvD2Offset("x_cv");
          assert(data_offset_xcc!=-1);
          lim->sparseReadCvR2_2(vol_xcc,cimS,"x_cv",data_offset_xcc);
        }
        else{
          assert(0);
        }
        double (*vol_vec)[3] = new double[cimS.ncv_unique][3];
        int8 data_offset = lim->getRfpCvD2Offset("u");
        if (data_offset!=-1){
          lim->sparseReadCvR2(vol_vec,cimS,"u",data_offset);
        }
        else if (filename_volumeData_2 != ""){
          data_offset = lim->getSfpCvD2Offset("u");
          assert(data_offset!=-1);
          lim->sparseReadCvR2_2(vol_vec,cimS,"u",data_offset);
        }
        else{
          assert(0);
        }
        const int ni = canvas->getNi();
        const int nj = canvas->getNj();
        float * lic_image = new float[ni*nj];

	//      canvas->projectVectorOnSurface(cimS, vol_vec);
        canvas->projectVectorOnSurfaceAveraged(cimS, vol_vec);
	//      canvas->projectVectorOnSphere(cimS,vol_vec,vol_xcc);
        computeLIC(cimS, lic_image,vol_vec,vol_xcc);
        canvas->addCimSurfaceLIC(cimS, lic_image);

        delete[] lic_image;
        delete[] vol_vec;
        delete[] vol_xcc;

      }
      else {
        //sparse read data and add to canvas
        CtiRegister::CtiData * data = CtiRegister::getCtiData(wid.var_on_surface);
        if ((data == NULL)||(data->getUnindexedTopology() != CV_DATA)||(data->getType() != DN_DATA)) {
          COUT1(" > Warning: requested expression " << wid.var << " does not eval to CVDN data. Skipping.");
        }
        else {
          canvas->addCimSurfaceData(cimS,data->getDNptr());
        }
      }

    }
  }

  void sparseReadAndAddVolumeData(){

    cout << "sparseReadAndAddVolumeData()..." << endl;

    if (cim.ncv_unique > 0) {
      cout << "  cim.ncv_unique > 0..." << endl;

      // need to make sure cim is current...
      current_cim_ptr = &cim;
      current_lim_ptr = lim;

      assert(lim!=NULL);
      imd->setTime(lim->getTime()); //set simulation time in image metadata

      if (wid.b_range){
        canvas->setDataRange(wid.range[0],wid.range[1]);
      }
      if (wid.range[2] != 0.f){
        canvas->setDataRangeBins((int)wid.range[2]);
      }

      if (wid.var=="MESH"||wid.var=="mesh"){
        canvas->addCimPlaneMesh(cim,lim->mesh_buf);
      }
      else if (wid.var=="LIC") {
        double (*vol_xcc)[3] = new double[cim.ncv_unique][3];
        int8 data_offset_xcc = lim->getRfpCvD2Offset("x_cv");
        if (data_offset_xcc!=-1){
          lim->sparseReadCvR2(vol_xcc,cim,"x_cv",data_offset_xcc);
        }
        else if (filename_volumeData_2 != ""){
          data_offset_xcc = lim->getSfpCvD2Offset("x_cv");
          assert(data_offset_xcc!=-1);
          lim->sparseReadCvR2_2(vol_xcc,cim,"x_cv",data_offset_xcc);
        }
        else{
          data_offset_xcc = lim->getRfpCvD2Offset("x_vv"); //legacy data
          if (data_offset_xcc!=-1){
            lim->sparseReadCvR2(vol_xcc,cim,"x_vv",data_offset_xcc);
          }
          else{
            assert(0);
          }
        }
        double (*vol_vec)[3] = new double[cim.ncv_unique][3];
        int8 data_offset = lim->getRfpCvD2Offset("u");
        if (data_offset!=-1){
          lim->sparseReadCvR2(vol_vec,cim,"u",data_offset);
        }
        else if (filename_volumeData_2 != ""){
          data_offset = lim->getSfpCvD2Offset("u");
          assert(data_offset!=-1);
          lim->sparseReadCvR2_2(vol_vec,cim,"u",data_offset);
        }
        else{
          assert(0);
        }

        const int ni = canvas->getNi();
        const int nj = canvas->getNj();
        float * lic_image = new float[ni*nj];

        canvas->projectVectorOnPlane(cim, vol_vec, wid.np);
        computeLIC(cim, lic_image,vol_vec,vol_xcc);
        canvas->addCimPlaneLIC(cim, lic_image);

        delete[] lic_image;
        delete[] vol_vec;
        delete[] vol_xcc;
      }
      else {
        CtiRegister::CtiData * data = CtiRegister::getCtiData(wid.var);
        if ((data == NULL)||(data->getUnindexedTopology() != CV_DATA)||(data->getType() != DN_DATA)) {
          COUT1(" > Warning: requested expression " << wid.var << " does not eval to CVDN data. Skipping.");
        }
        else {
          cout << "  canvas->addCimPlaneData(cim,data->getDNptr())" << endl;
          canvas->addCimPlaneData(cim,data->getDNptr());
        }
      }
    }
  }

  // perform line integral convolution
  void DoConvolution(int width,int height,float* vector,float* mag,int* mask,unsigned char* image) {
    //TODO
    // for now LIC parameters are hard-coded, but in the future we would like these
    // to be user defined, even perhaps the vector-variables used in the convolution
    float LENGTH_RANGE[2] = {1.E10f,-1.E10f};
    const float LENGTH_MAX = 75.f;
    const float LENGTH_MIN = 10.f;
    int DISCRETE_FILTER_SIZE  = 2048;
    float LINE_SQUARE_CLIP_MAX = 100000.0;
    float VMIN = 0.001000;
    COUT2("    > filter size: " << DISCRETE_FILTER_SIZE);
    COUT2("    > line square clip: " << LINE_SQUARE_CLIP_MAX);
    // float VMIN = 0.050000;

    int  vec_id;      ///ID in the VECtor buffer (for the input flow field)
    int  advDir;      ///ADVection DIRection (0: positive;  1: negative)
    int  advcts;      ///number of ADVeCTion stepS per direction (a step counter)
    int  ADVCTS = (int)(LENGTH_MAX * 3); ///MAXIMUM number of advection steps per direction to break dead loops

    float vctr_x;      ///x-component  of the VeCToR at the forefront point
    float vctr_y;      ///y-component  of the VeCToR at the forefront point
    float clp0_x;      ///x-coordinate of CLiP point 0 (current)
    float clp0_y;      ///y-coordinate of CLiP point 0 (current)
    float clp1_x;      ///x-coordinate of CLiP point 1 (next   )
    float clp1_y;      ///y-coordinate of CLiP point 1 (next   )
    float samp_x;      ///x-coordinate of the SAMPle in the current pixel
    float samp_y;      ///y-coordinate of the SAMPle in the current pixel
    float tmpLen;      ///TeMPorary LENgth of a trial clipped-segment
    float segLen;      ///SEGment   LENgth
    float curLen;      ///CURrent   LENgth of the streamline
    float prvLen;      ///PReVious  LENgth of the streamline
    float W_ACUM;      ///ACcuMulated Weight from the seed to the current streamline forefront
    float texVal;      ///TEXture VALue
    float smpWgt;      ///WeiGhT of the current SaMPle
    float t_acum[2];     ///two ACcUMulated composite Textures for the two directions, perspectively
    float w_acum[2];     ///two ACcUMulated Weighting values   for the two directions, perspectively
    float* wgtLUT = NULL;    ///WeiGhT Look Up Table pointing to the target filter LUT

    unsigned char* noise  = (unsigned char* ) malloc( sizeof(unsigned char) * width * height     );
    float*         p_LUT0 = (float*   ) malloc( sizeof(float        ) * DISCRETE_FILTER_SIZE);
    float*         p_LUT1 = (float*   ) malloc( sizeof(float        ) * DISCRETE_FILTER_SIZE);

    // store filter
    for(int i = 0;  i < DISCRETE_FILTER_SIZE;  i ++)  p_LUT0[i] = p_LUT1[i] = i;

    // make noise...
    for( int j = 0;   j < height;  j ++) {
      for( int i = 0;   i < width;  i ++) {
        int  r = rand();
        r = (  (r & 0xff) + ( (r & 0xff00) >> 8 )  ) & 0xff;
        noise[j * width + i] = (unsigned char) r;
      }
    }

    ///for each pixel in the 2D output LIC image///
    for( int j = 0;  j < height;  j ++) {
      for( int i = 0;  i < width;  i ++) {
        if (mask[j * width + i] == 0) {

          // float krnlen = LENGTH_MAX;
          float krnlen = ceil(mag[j * width + i] * LENGTH_MAX);
          krnlen = max(LENGTH_MIN,min(LENGTH_MAX,krnlen));
          LENGTH_RANGE[0] = min(LENGTH_RANGE[0],krnlen);
          LENGTH_RANGE[1] = max(LENGTH_RANGE[1],krnlen);

          float len2ID = (DISCRETE_FILTER_SIZE - 1) / krnlen; ///map a curve LENgth TO an ID in the LUT

          ///init the composite texture accumulators and the weight accumulators///
          t_acum[0] = t_acum[1] = w_acum[0] = w_acum[1] = 0.0f;

          ///for either advection direction///
          for(advDir = 0;  advDir < 2;  ++advDir) {
            ///init the step counter, curve-length measurer, and streamline seed///
            advcts = 0;
            curLen = 0.0f;
            clp0_x = i + 0.5f;
            clp0_y = j + 0.5f;

            ///access the target filter LUT///
            wgtLUT = (advDir == 0) ? p_LUT0 : p_LUT1;

            ///until the streamline is advected long enough or a tightly  spiralling center / focus is encountered///
            while( (curLen < krnlen) && (advcts < ADVCTS) ) {
              ///access the vector at the sample///
              vec_id = ( (int)(clp0_y) * width + (int)(clp0_x) ) << 1;

              // oe possibility is to skip the update of the velocity if the current point is masked!
              if (mask[(int)(clp0_y) * width + (int)(clp0_x)] == 0) {
                vctr_x = vector[vec_id    ];
                vctr_y = vector[vec_id + 1];
              }

              ///in case of a critical point///
              if( (vctr_x == 0.0f) && (vctr_y == 0.0f) ) {
                t_acum[advDir] = (advcts == 0) ? 0.0f : t_acum[advDir];
                w_acum[advDir] = (advcts == 0) ? 1.0f : w_acum[advDir];
                break;
              }

              ///negate the vector for the backward-advection case///
              vctr_x = (advDir == 0) ? vctr_x : -vctr_x;
              vctr_y = (advDir == 0) ? vctr_y : -vctr_y;

              ///clip the segment against the pixel boundaries --- find the shorter from the two clipped segments///
              segLen = LINE_SQUARE_CLIP_MAX;
              segLen = (vctr_x < -VMIN) ? ((int)(      clp0_x      )-clp0_x)/vctr_x : segLen;
              segLen = (vctr_x >  VMIN) ? ((int)((int)(clp0_x)+1.5f)-clp0_x)/vctr_x : segLen;
              segLen = (vctr_y < -VMIN) ? (((tmpLen = ((int)(      clp0_y)      -clp0_y)/vctr_y) < segLen) ? tmpLen : segLen) : segLen;
              segLen = (vctr_y >  VMIN) ? (((tmpLen = ((int)((int)(clp0_y)+1.5f)-clp0_y)/vctr_y) < segLen) ? tmpLen : segLen) : segLen;

              ///update the curve-length measurers///
              prvLen = curLen;
              curLen+= segLen;
              segLen+= 0.0004f;

              ///check if the filter has reached either end///
              segLen = (curLen > krnlen) ? ( (curLen = krnlen) - prvLen ) : segLen;

              ///obtain the next clip point///
              clp1_x = clp0_x + vctr_x * segLen;
              clp1_y = clp0_y + vctr_y * segLen;

              ///obtain the middle point of the segment as the texture-contributing sample///
              samp_x = (clp0_x + clp1_x) * 0.5f;
              samp_y = (clp0_y + clp1_y) * 0.5f;

              ///obtain the texture value of the sample///
              const int sample_loc = min( width*height - 1, int( round(samp_y)*width + round(samp_x)) );
              texVal = noise[sample_loc];

              ///update the accumulated weight and the accumulated composite texture (texture x weight)///
              W_ACUM = wgtLUT[ (int)(curLen * len2ID) ];
              smpWgt = W_ACUM - w_acum[advDir];
              w_acum[advDir]  = W_ACUM;
              t_acum[advDir] += texVal * smpWgt;

              ///update the step counter and the "current" clip point///
              advcts ++;
              clp0_x = clp1_x;
              clp0_y = clp1_y;

              ///check if the streamline has gone beyond the flow field///
              if( clp0_x < 0.0f || clp0_x >= width || clp0_y < 0.0f || clp0_y >= height)  break;
            }
          }

          ///normalize the accumulated composite texture///
          texVal = (t_acum[0] + t_acum[1]) / (w_acum[0] + w_acum[1]);

          ///clamp the texture value against the displayable intensity range [0, 255]
          texVal = (texVal <   0.0f) ?   0.0f : texVal;
          texVal = (texVal > 255.0f) ? 255.0f : texVal;
          image[j * width + i] = (unsigned char) texVal;
        }
      }
    }
    free(noise); noise = NULL;
    free(p_LUT0); p_LUT0 = NULL;
    free(p_LUT1); p_LUT1 = NULL;

    COUT2("    > integration length range: [" << LENGTH_RANGE[0] << "," << LENGTH_RANGE[1] << "]");
  }

  void computeLIC(const CvImageMap &cim, float *lic_image,const double (*vol_vec)[3], const double (*vol_xcc)[3]){

    float x0[3];  //simulation coordinates, center of the view
    float e0[3];  //NOTE, not unit vectors, have length ni/width
    float e1[3];
    float e2[3];
    canvas->getX0(x0);
    canvas->getE0(e0);
    canvas->getE1(e1);
    canvas->getE2(e2);
    const int ni = canvas->getNi();
    const int nj = canvas->getNj();

    // GI: construct the mask....and for now store the fields in arrays...
    float *gi_vector  = new float[ni*nj*2];
    float *gi_mag  = new float[ni*nj];
    float vmag_max = -1.0f;
    int *gi_mask = new int[ni*nj];
    unsigned char* gi_image = new unsigned char[ni*nj];

    // initialize the map...mask is ON everywhere, unless there are cells...
    for (int ij=0; ij<ni*nj; ij++) gi_mask[ij] = 1;

    for (int icv = 0; icv < cim.ncv_unique; ++icv) {
      for (int poc = cim.poucv_i[icv]; poc != cim.poucv_i[icv+1]; ++poc) {
        unsigned int index = cim.getIndex(poc);
        unsigned int repeat_count = cim.getRepeatCount(poc);
        assert( (index | (repeat_count<<(32-CIM_REPEAT_BITS))) == cim.poucv_v[poc]);
        for (unsigned int offset = 0; offset <= repeat_count; ++offset) {
          const int image_ij = index+offset;
          const int image_j = (int) image_ij/ni;
          const int image_i = image_ij%ni;

          // GI: store the active pixels..(deactivae the mask...)
          gi_mask[image_j*ni+image_i] = 0;
          // GI: store the vector and normalize....
          int  index = (image_j * ni + image_i) << 1;

          // convert xcc and ut into image coordinates...
          double p[3],pv[3];

          FOR_I3 p[i]=vol_xcc[icv][i];
          FOR_I3 pv[i]=vol_xcc[icv][i]+vol_vec[icv][i];

          const double i_p = ( (p[0]-x0[0])*e0[0] + (p[1]-x0[1])*e0[1] + (p[2]-x0[2])*e0[2] );
          const double j_p = ( (p[0]-x0[1])*e1[0] + (p[1]-x0[1])*e1[1] + (p[2]-x0[2])*e1[2] );
          const double i_pv = ( (pv[0]-x0[0])*e0[0] + (pv[1]-x0[1])*e0[1] + (pv[2]-x0[2])*e0[2] );
          const double j_pv = ( (pv[0]-x0[1])*e1[0] + (pv[1]-x0[1])*e1[1] + (pv[2]-x0[2])*e1[2] );

          double u = i_pv-i_p;
          double v = j_pv-j_p;
          double vmag = sqrt(u*u+v*v);
          gi_mag[image_j*ni+image_i] = float(vmag);
          if (vmag < 1.0E-10) {
            u=1.0E-10;
            v=1.0E-10;
            gi_mag[image_j*ni+image_i] = 1.0E-10;
          }
          else { u = u/vmag; v=v/vmag; }
          vmag_max = max(vmag_max,gi_mag[image_j*ni+image_i]);

          gi_vector[index    ] =  u;
          gi_vector[index + 1] =  v;

        }
      }
    }

    // scale velocity magnitude by max in image
    for (int ij=0; ij<ni*nj; ij++) gi_mag[ij] /= vmag_max;

    COUT1(" > computing line integral convolution (LIC)");
    DoConvolution(ni,nj,gi_vector,gi_mag,gi_mask,gi_image);

    for (int image_i=0;image_i<ni;++image_i){
      for (int image_j=0;image_j<nj;++image_j){
        lic_image[image_j * ni + image_i]=  (double)(gi_image[image_j * ni + image_i])/255.0;
      }
    }

    delete[] gi_mask;
    delete[] gi_vector;
    delete[] gi_mag;
    delete[] gi_image;
    // GI: done...

  }

  // function added for CtiRegister
  CtiRegister::CtiData * varEvalCtiData(const string& name) {

    if (mpi_rank == 0) cout << " > CtiScene::varEvalCtiData(\"" << name << "\")..." << endl;

    pair<const string,CtiRegister::CtiData> key_value_pair(name,CtiRegister::CtiData());
    pair<map<const string,CtiRegister::CtiData>::iterator,bool> return_pair = CtiRegister::currentDataMap.insert(key_value_pair);
    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(return_pair.second); // should insert successfully
    assert(key_value_pair.second.empty()); // should be empty - don't use anymore, and his destructor cannot affect the inserted CtiData

    CtiRegister::CtiData * data = &(return_pair.first->second);
    assert(data->empty());

    int8 data_offset = current_lim_ptr->getRfpCvD1Offset(name);
    if (data_offset!=-1){
      data->new_dn(CV_DATA,current_cim_ptr->ncv_unique); // not really CV_DATA (actually a subset)
      current_lim_ptr->sparseReadCvR1(data->getDNptr(),*current_cim_ptr,name,data_offset);
    }
    else{
      data_offset = current_lim_ptr->getRfpCvD2Offset(name);
      if (data_offset!=-1){
        data->new_dn3(CV_DATA,current_cim_ptr->ncv_unique);
        current_lim_ptr->sparseReadCvR2(data->getDN3ptr(),*current_cim_ptr,name,data_offset);
      }
      else if (filename_volumeData_2 != "") {
        data_offset = current_lim_ptr->getSfpCvD1Offset(name);
        if (data_offset!=-1){
          data->new_dn(CV_DATA,current_cim_ptr->ncv_unique);
          current_lim_ptr->sparseReadCvR1_2(data->getDNptr(),*current_cim_ptr,name,data_offset);
        }
        else{
          data_offset = current_lim_ptr->getSfpCvD2Offset(name);
          if (data_offset!=-1) {
            data->new_dn3(CV_DATA,current_cim_ptr->ncv_unique);
            current_lim_ptr->sparseReadCvR2_2(data->getDN3ptr(),*current_cim_ptr,name,data_offset);
          }
          else {
            COUT1(" > Warning, could not find var " << name << " in *.sles file");
            CtiRegister::currentDataMap.erase(return_pair.first);
          }
        }
      }
      else {
        COUT1(" > Warning, could not find var " << name << " in *.mles,*.les file");
        CtiRegister::currentDataMap.erase(return_pair.first);
      }
    }

    return data;
  }

  virtual CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,const string& name,list<CtiRegister::CtiData>& args,
                                                    const bool b_eval_func) {

    return CtiRegister::CTI_DATA_NOT_FOUND;

  }

  void coordinateTransform(double x_t[3], const double x[3]) {

    if (wid.b_transform_cyl_x) {
      GeomUtils::coordinateTransform(x_t,x,wid.r0,wid.cos_theta0,wid.sin_theta0,wid.cyl_center,GeomUtils::TRANSFORM_CYL_X);
    }
    else if (wid.b_transform_cyl_z) {
      GeomUtils::coordinateTransform(x_t,x,wid.r0,wid.cos_theta0,wid.sin_theta0,wid.cyl_center,GeomUtils::TRANSFORM_CYL_Z);
    }
    else {
      FOR_I3 x_t[i] = x[i];
    }

  }

  void coordinateTransformAndAdjustTheta(double x_t[3], const double x[3], const double x_ref[3]) {

    if (wid.b_transform_cyl_x) {
      GeomUtils::coordinateTransform(x_t,x,wid.r0,wid.cos_theta0,wid.sin_theta0,wid.cyl_center,GeomUtils::TRANSFORM_CYL_X);
      if (x_t[2] < x_ref[2]) {
	if ((x_t[2]+2.0*M_PI*wid.r0)-x_ref[2] < x_ref[2]-x_t[2])
	  x_t[2] += 2.0*M_PI*wid.r0;
      }
      else if (x_t[2] > x_ref[2]) {
	if (x_ref[2]-(x_t[2]-2.0*M_PI*wid.r0) < x_t[2]-x_ref[2])
	  x_t[2] -= 2.0*M_PI*wid.r0;
      }
    }
    else if (wid.b_transform_cyl_z) {
      GeomUtils::coordinateTransform(x_t,x,wid.r0,wid.cos_theta0,wid.sin_theta0,wid.cyl_center,GeomUtils::TRANSFORM_CYL_Z);
      if (x_t[1] < x_ref[1]) {
	if ((x_t[1]+2.0*M_PI*wid.r0)-x_ref[1] < x_ref[1]-x_t[1])
	  x_t[1] += 2.0*M_PI*wid.r0;
      }
      else if (x_t[1] > x_ref[1]) {
	if (x_ref[1]-(x_t[1]-2.0*M_PI*wid.r0) < x_t[1]-x_ref[1])
	  x_t[1] -= 2.0*M_PI*wid.r0;
      }
    }
    else {
      FOR_I3 x_t[i] = x[i];
    }

  }

  void adjustThetaTri(double x0[3], double x1[3], double x2[3]) {

    if (wid.b_transform_cyl_x) {
      // ensure that the tri does not cross over the -180/+180 boundary
      // find the 2 closest nodes in theta...
      const double dt_01 = x1[2]-x0[2];
      const double dt_12 = x2[2]-x1[2];
      const double dt_20 = x0[2]-x2[2];
      if (fabs(dt_01) < min(fabs(dt_12),fabs(dt_20))) {
	// 0,1 are close, so try moving 2 and see what happens...
	if (dt_12 < 0.0) {
	  if (fabs(dt_12+2.0*M_PI*wid.r0) < -dt_12)
	    x2[2] += 2.0*M_PI*wid.r0;
	}
	else {
	  if (fabs(dt_12-2.0*M_PI*wid.r0) < dt_12)
	    x2[2] -= 2.0*M_PI*wid.r0;
	}
      }
      else if (fabs(dt_12) < fabs(dt_20)) {
	// 1,2 are close, so try moving 0 and see what happens...
	if (dt_20 < 0.0) {
	  if (fabs(dt_20+2.0*M_PI*wid.r0) < -dt_20)
	    x0[2] += 2.0*M_PI*wid.r0;
	}
	else {
	  if (fabs(dt_20-2.0*M_PI*wid.r0) < dt_20)
	    x0[2] -= 2.0*M_PI*wid.r0;
	}
      }
      else {
	// 2,0 are closest. Test 1...
	if (dt_01 < 0.0) {
	  if (fabs(dt_01+2.0*M_PI*wid.r0) < -dt_01)
	    x1[2] += 2.0*M_PI*wid.r0;
	}
	else {
	  if (fabs(dt_01-2.0*M_PI*wid.r0) < dt_01)
	    x1[2] -= 2.0*M_PI*wid.r0;
	}
      }
    }
    else if (wid.b_transform_cyl_z) {
      // ensure that the tri does not cross over the -180/+180 boundary
      // find the 2 closest nodes in theta...
      const double dt_01 = x1[1]-x0[1];
      const double dt_12 = x2[1]-x1[1];
      const double dt_20 = x0[1]-x2[1];
      if (fabs(dt_01) < min(fabs(dt_12),fabs(dt_20))) {
	// 0,1 are close, so try moving 2 and see what happens...
	if (dt_12 < 0.0) {
	  if (fabs(dt_12+2.0*M_PI*wid.r0) < -dt_12)
	    x2[1] += 2.0*M_PI*wid.r0;
	}
	else {
	  if (fabs(dt_12-2.0*M_PI*wid.r0) < dt_12)
	    x2[1] -= 2.0*M_PI*wid.r0;
	}
      }
      else if (fabs(dt_12) < fabs(dt_20)) {
	// 1,2 are close, so try moving 0 and see what happens...
	if (dt_20 < 0.0) {
	  if (fabs(dt_20+2.0*M_PI*wid.r0) < -dt_20)
	    x0[1] += 2.0*M_PI*wid.r0;
	}
	else {
	  if (fabs(dt_20-2.0*M_PI*wid.r0) < dt_20)
	    x0[1] -= 2.0*M_PI*wid.r0;
	}
      }
      else {
	// 2,0 are closest. Test 1...
	if (dt_01 < 0.0) {
	  if (fabs(dt_01+2.0*M_PI*wid.r0) < -dt_01)
	    x1[1] += 2.0*M_PI*wid.r0;
	}
	else {
	  if (fabs(dt_01-2.0*M_PI*wid.r0) < dt_01)
	    x1[1] -= 2.0*M_PI*wid.r0;
	}
      }
    }

  }

  void getTransformInfo(int& transform_type, double& r0_, double& cos_theta0_, double& sin_theta0_, double cyl_center_[3]) {

    if (wid.b_transform_cyl_x)
      transform_type = 1;
    else if (wid.b_transform_cyl_z)
      transform_type = 3;
    else
      transform_type = 0;

    if (transform_type > 0) {
      r0_ = wid.r0;
      cos_theta0_ = wid.cos_theta0;
      sin_theta0_ = wid.sin_theta0;
      FOR_I3 cyl_center_[i] = wid.cyl_center[i];
    }
    else {
      // no transform. set to HUGE_VAL
      r0_ = HUGE_VAL;
      cos_theta0_ = HUGE_VAL;
      sin_theta0_ = HUGE_VAL;
      FOR_I3 cyl_center_[i] = HUGE_VAL;
    }
  }

};


#endif
