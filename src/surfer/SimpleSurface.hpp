#ifndef _SIMPLE_SURFACE_HPP_
#define _SIMPLE_SURFACE_HPP_

#include <stack>
#include <sstream>
#include "IntFlag.hpp"
#include "Teost.hpp"
#include "NonManifoldTriFlag.hpp"
#include "DoubleVertex.hpp"
#include "CTI.hpp"
#include "MiscUtils.hpp"
#include "MinHeap.hpp"
#include "RestartHashUtilities.hpp"
#include "Adt.hpp"
#include "Adt2d.hpp"
#include "DirectedGraph.hpp"
#include "SimpleGeom.hpp"
#include "WebUI.hpp"
#include "GeomUtils.hpp"
#include "../stitch/IntersectionStuff.hpp"

using namespace CTI;

#include "PeriodicData.hpp"

#include "NewTri.hpp"
#include "NewNode.hpp"

class SimpleSurface : public SimpleGeom {

private:

  // the area-weighted centroid of the surface
  bool b_centroid;
  double centroid[3];

  // the maximum distance from the centroid to any node on the surface
  bool b_rmax;
  double rmax;

  // the x,y,z bounding box (min/max extent) of the surface
  bool b_boundingBox;
  double boundingBox[6]; // xmin,xmax,ymin,yamx,zmin,zmax

  // indicates whether non-manifold issues were detected during teost construction
  NonManifoldTriFlag nmTriData;

  // tri-edge data structure
  bool b_teost;
  Teost teost;

  bool b_open_edge_groups;
  int n_open_edge_groups;

  bool b_pbi_bits; // periodic bits set for ALL transforms in periodicTransformVec
  uint8 * pbi; // periodic-bits and index: top 12 bits are periodic bits, remaining 52 are the index of the original point

public:

  // in PeriodicData namespace now...
  //vector<PeriodicTransform> periodicTransformVec; // need access for JSON

  /* ======================================================
  // tri nodes (spost) and edges (in teost_data for example)
  // are indexed as folows...
  //
  //      2                    1
  //      |\                   |\
  //      | \        flip      | \
  //    e2|  \e1     --->    e0|  \e1
  //      |   \                |    \
  //      |    \               |    \
  //      0-----1              0-----2
  //        e0                   e2
  //
  //     normal out           normal into
  //     of screen              screen
  //     (rh-rule)
  //
  // ====================================================== */

  // surface data
  int nsp;
  double (*xsp)[3];
  int nst;
  int (*spost)[3];
  int *znost;

  // material data
  enum material_type {
    FLUID,
    MAT_MIN_VAL=FLUID,
    SOLID,
    MAT_MAX_VAL=SOLID
  };

  class VolumeMaterial {
  private:
    string name;
    material_type type;  // 0: fluid, 1: solid (others open for different cv-zone types, i.e., sponge, porous-media)

  public:
    set<int> boundary_zones;  // store which zones have been assigned to me

    VolumeMaterial (const string& name_, const material_type type_=FLUID) {
      name = name_;
      type = type_;
      boundary_zones.clear();
    }

    VolumeMaterial() {
      name = "";
      type = FLUID;  // default to fluid
      boundary_zones.clear();
    }

    ~VolumeMaterial() {}

    string getName() const { return name; }
    void setName(const string& name_) { name = name_; }

    material_type getType() const {return type; }
    void setType(const material_type type_) {
      type = type_;
    }

    void addZone(const int zone_index) {
      boundary_zones.insert(zone_index);
    }
    bool rmZone(const int zone_index) {
      bool found = false;
      set<int>::iterator it;
      for (it=boundary_zones.begin(); it!=boundary_zones.end(); ++it) {
        if (*it == zone_index) {
          found = true;
          break;
        }
      }

      if (found) {
        boundary_zones.erase(it);
      }
      return found;
    }
    void clearZones() {
      boundary_zones.clear();
    }
  };

  vector<VolumeMaterial> materialVec;

  // zone data
  class SurfaceZone {
  private:

    string name;
    int level;
    uint2 periodic_bits; // bits for any associated periodicity
    string boundary_condition;

    vector<bool> b_metadata;  // centroid, normal, area
    double xc[3];
    double normal[3];
    double area;

    int material[2];  // index into material vec

  public:
    int flag; // multi-purpose flag
    SurfaceZone(const string& _name) {
      name = _name;
      MiscUtils::eraseAllSubStr(name, "\"");
      level = 0;
      periodic_bits = 0;
      boundary_condition = "WALL_ADIABATIC";
      material[0] = 0;
      material[1] = -1;  // nothing
      clearMetadata();
    }
    SurfaceZone() {
      name = "";
      level = 0;
      periodic_bits = 0;
      boundary_condition = "WALL_ADIABATIC";
      material[0] = 0;
      material[1] = -1;  // nothing
      clearMetadata();
    }

    string getName() const { return name; }
    void setName(const string& _name) {
      name = _name;
      MiscUtils::eraseAllSubStr(name, "\"");
    }

    int getLevel() const { return level; }
    void setLevel(const int& _level) { level = _level; }

    int silver() {
      return getSilverMaterial();
    }
    int getSilverMaterial() const {
      return material[0];
    }
    int gold() {
      return getGoldMaterial();
    }
    int getGoldMaterial() const {
      return material[1];
    }
    bool setMaterial(const int mat_index, const int side) {
      bool set = false;
      if ((side==0) || (side==1)) {
        material[side] = mat_index;
        set = true;
      }
      return set;
    }
    bool setSilver(const int mat_index) {
      return setMaterial(mat_index,0);
    }
    bool setGold(const int mat_index) {
      return setMaterial(mat_index,1);
    }
    void resetMaterial() {
      setSilver(0);
      setGold(-1);
    }

    bool getCentroid(double (&_xc)[3]) const {
      if (!b_metadata[0]) return false;

      FOR_I3 _xc[i] = xc[i];
      return true;
    }
    void setCentroid(double (&_xc)[3]) {
      FOR_I3 xc[i] = _xc[i];
      b_metadata[0] = true;
    }

    bool getNormal(double (&_normal)[3]) const {
      if (!b_metadata[1]) return false;

      FOR_I3 _normal[i] = normal[i];
      return true;
    }
    void setNormal(double (&_normal)[3]) {
      FOR_I3 normal[i] = _normal[i];
      b_metadata[1] = true;
    }

    bool getArea(double &_area) const {
      if (!b_metadata[2]) return false;

      _area = area;
      return true;
    }
    void setArea(double& _area) {
      area = _area;
      b_metadata[2] = true;
    }

    bool hasMetadata() const {
      for (int i=0, limit=b_metadata.size(); i<limit; ++i) {
        if (!b_metadata[i]) return false;
      }
      // all values have been set
      return true;
    }

    void clearMetadata() {
      b_metadata.clear();
      b_metadata.resize(3,false);
    }

    string getBC() const {return boundary_condition; }
    void setBC(const string& _bc) {
      boundary_condition = _bc;
    }
    bool isBCPeriodic() const {
      const string prefix = boundary_condition.substr(0,5);
      if (prefix == "PER1_" || prefix == "PER2_" || prefix == "PER3_") return true;
      else return false;
    }

    uint2 getPeriodicBits() const { return periodic_bits; }
    void setPeriodicBits(const uint2& _periodic_bits) { periodic_bits = _periodic_bits; }
    void dump() const {
      cout << "SurfaceZone \"" << name << "\" level: " << level;
      if (periodic_bits) {
        cout << " periodic_bits: ";
        for (int i = 5; i >= 0; --i) {
          if (periodic_bits & (1<<i)) cout << "1";
          else cout << "0";
        }
      }
      cout << endl;
    }
  };
  bool b_zone_data;
  vector<SurfaceZone> zoneVec;
  double zone_level_hcp_delta;

  // subzone data
  class SubzoneData {
  public:
    double area;
    double normal[3];
    double xc[3];
    //double r_rms; // another moment?

    int zone;  // parent zone

    // mesh related values
    bool has_gap_h;
    double gap_h;
    bool has_curvature;
    double curvature;
    int level;

    SubzoneData() {
      zero();
    }
    void zero() {
      area = 0.0;
      normal[0] = 0.0;
      normal[1] = 0.0;
      normal[2] = 0.0;
      xc[0] = 0.0;
      xc[1] = 0.0;
      xc[2] = 0.0;
      has_gap_h = false;
      gap_h = 0.0;
      has_curvature = false;
      curvature = 0.0;
      //r_rms = 0.0;
      level = 0;  // default is to background resolution
      zone = -1;
    }
    void dump() const {
      cout << " > surface area: " << area << endl;
      cout << " > normal vec: " << COUT_VEC(normal) << endl;
      cout << " > mag(normal): " << MAG(normal) << endl;
      cout << " > x_centroid: " << COUT_VEC(xc) << endl;
      cout << " > parent zone (index): " << zone << endl;
      if (has_gap_h) cout << " > gap height: " << gap_h << endl;
      if (has_curvature) cout << " > curvature (10% CDF): " << curvature << endl;
      if (level) cout << " > suggested Level: " << level << endl;
    }
  };

  // edge groups of interest data
  class HalfEdgeGroupData {
  public:
    int ned;
    // store the first and last edge (ist/i) using compact half-edge math
    pair<bool,uint> end_ed[2]; // bool: does it exist, if yes uint: ist (30 bits) | i (2-bits)
    double length;
    double xc[3];
    double mean_d;
    set<int> zones;

    HalfEdgeGroupData() {
      zero();
    }

    void zero() {
      ned = 0;
      end_ed[0] = pair<bool,uint> (false,0);
      end_ed[1] = pair<bool,uint> (false,0);
      length = 0.0;
      xc[0] = 0.0;
      xc[1] = 0.0;
      xc[2] = 0.0;
      mean_d = 0.0;
      zones.clear();
    }

    void dump() const {
      cout << "halfEdgeGroupData: edges: " << ned << ", length: " << length;
      cout << ", xc: " << COUT_VEC(xc) << ", mean_dist: " << mean_d;
      const int ist0_i = end_ed[0].second & 3;
      const int ist0 = (end_ed[0].second >> 2) & MASK_30BITS;
      const int ist1_i = end_ed[1].second & 3;
      const int ist1 = (end_ed[1].second >> 2) & MASK_30BITS;
      cout << ", end_ed[0] (ist,i): " << ist0 << "," << ist0_i << " end_ed[1] (ist,i): " << ist1 << "," << ist1_i;
      if (((end_ed[0].first)&&(!end_ed[1].first))||
          ((!end_ed[0].first)&&(end_ed[1].first))) {
        cout << ", This group had problems during indexing process.";
      }
      //assert(end_ed[1] != -1); // XXX BAD
      cout << endl;
    }
  };

  // these are always available...
  int nsz;
  IntFlag szost; // subzone of surface tri
  IntFlag szozn_i; // subzone of zone index range. Note: nsz == szozn_i[zoneVec.size()]

  bool b_subzone_data;
  vector<SubzoneData> subzoneDataVec;

  // container for some topo and geom of open edge group (build w/ open edge group)
  vector<HalfEdgeGroupData> openEdgeGroupDataVec;
  map<uint,int> oe_to_group;

  double feature_cos;
  // container to hold group indexing for dynamically defined edges-of-interest
  enum DYNAMIC_EDGE_TYPE {
    NONE,
    MULTI,
    METAL_BOUNDARY,
    FLUID_BOUNDARY,
    ZONE_BOUNDARY,
    SUBZONE_BOUNDARY,
    FEATURE,
    OPEN,
    NME,
    SELECTED_BOUNDARY
  };
  bool b_eoi;
  DYNAMIC_EDGE_TYPE eoi_type;
  map<uint,int> eoi_to_group;
  vector<HalfEdgeGroupData> eoiGroupDataVec;

  vector<int> selectedSubzoneVec;  // container for what zones to highlight in UI
  bool b_update_hidden_subzones;
  vector<int> hiddenSubzoneVec;  // when subzone indices change, populate this

  bool b_stosz_szosz;
  int * stosz_i;
  int * stosz_v;
  int * szosz_i;
  int * szosz_v;

  int n_closed_edge_groups; // number of closed open edge groups (needs to be persistent?)

  // various flags...
  IntFlag st_flag;
  IntFlag sp_flag;
  IntFlag zone_flag;
  IntFlag sz_flag;

  bool b_sposp;
  int * sposp_i;  // point-point connectivity
  int * sposp_v;  // not built by default
  double *sp_buf;
  double *st_buf;
  double (*sp_buf3)[3];
  double (*st_buf3)[3];

  SimpleSurface() {

    b_centroid = false;
    b_rmax = false;
    b_boundingBox = false;

    nsp = nst = 0;
    xsp = NULL;
    spost = NULL;

    b_teost = false;
    znost = NULL;

    b_open_edge_groups = false;
    n_open_edge_groups = 0;
    n_closed_edge_groups = 0;

    clearMaterials();

    zone_level_hcp_delta = 0.0;

    nsz = 0;
    b_subzone_data = false;
    b_zone_data = false;
    b_update_hidden_subzones = false;

    b_stosz_szosz = false;
    stosz_i = NULL;
    stosz_v = NULL;
    szosz_i = NULL;
    szosz_v = NULL;

    b_sposp = false;
    sposp_i = NULL;
    sposp_v = NULL;

    sp_buf = NULL;
    st_buf = NULL;
    sp_buf3 = NULL;
    st_buf3 = NULL;

    b_pbi_bits = false;
    pbi = NULL;

    stAdt2d = NULL;

  }

  ~SimpleSurface() {

    DELETE(xsp);
    DELETE(spost);
    DELETE(znost);

    DELETE(stosz_i);
    DELETE(stosz_v);
    DELETE(szosz_i);
    DELETE(szosz_v);
    DELETE(sposp_i);
    DELETE(sposp_v);
    DELETE(sp_buf);
    DELETE(st_buf);
    DELETE(sp_buf3);
    DELETE(st_buf3);

    DELETE(pbi);

    if (stAdt2d != NULL) delete stAdt2d;

  }

  void check() {
    cout << "SimpleSurface::check(): nst: " << nst << " nsp: " << nsp << " nzn: " << zoneVec.size() << " nsz: " << nsz; cout.flush();

    /*
    assert(nsz == szozn_i[zoneVec.size()]);
    for (int ist = 0; ist < nst; ++ist) {
      const int izn = znost[ist];
      assert((izn >= 0)&&(izn < int(zoneVec.size())));
      const int isz = szost[ist];
      assert((isz >= szozn_i[izn])&&(isz < szozn_i[izn+1]));
    }
    */

    cout << " OK" << endl;

    // hidden subzone stuff
    // cout << "hidden subzones: ";
    // for (int isz=0,end=hiddenSubzoneVec.size(); isz<end; ++isz) cout << hiddenSubzoneVec[isz] << ",";
    // cout << endl;
  }

  void clear() {
    COUT1("SimpleSurface::clear()");

    zoneVec.clear();
    clearMaterials();

    b_centroid = false;
    b_rmax = false;
    b_boundingBox = false;

    nsp = nst = 0;
    zone_level_hcp_delta = 0.0;

    DELETE(xsp);
    DELETE(spost);
    DELETE(znost);

    b_teost = false;

    b_open_edge_groups = false;
    n_open_edge_groups = 0;
    n_closed_edge_groups = 0;

    feature_cos = cos((180.0-150.0)*M_PI/180.0);
    b_eoi = false;
    eoi_type = NONE;

    // clear EdgeGroup data
    eoi_to_group.clear();
    eoiGroupDataVec.clear();
    openEdgeGroupDataVec.clear();

    st_flag.clear();
    sp_flag.clear();
    zone_flag.clear();
    szost.clear();
    szozn_i.clear();

    // clear Surface data
    clearZoneData();
    clearSubzoneData();
    clearDynamicEdgeGroups();

    b_sposp = false;
    DELETE(sposp_i);
    DELETE(sposp_v);

    b_stosz_szosz = false;
    DELETE(stosz_i);
    DELETE(stosz_v);
    DELETE(szosz_i);
    DELETE(szosz_v);

    DELETE(sp_buf);
    DELETE(st_buf);
    DELETE(sp_buf3);
    DELETE(st_buf3);

    PeriodicData::periodicTransformVec.clear();
    DELETE(pbi);
  }

  int init(Param * param) {

    COUT1("SimpleSurface::init()");

    // we are about to initialize or add a surface, so clear
    // anything that is currently set based on these...

    clearSubzoneData();
    clearZoneData();
    clearDynamicEdgeGroups();
    b_centroid = false;
    b_rmax = false;
    b_boundingBox = false;
    clearTeost();

    int ierr = 0;

    int iarg = 0;
    string token = MiscUtils::toUpperCase(param->getString(iarg++));
    if (token == "STL_GROUP") {
      while (iarg != param->size()) {
        const string name = param->getString(iarg++);
        COUT1(" > about to read stl filename: " << name);
        ierr = addStl(name);
        if (ierr != 0) COUT1(" ! skipping file " << name);
      }
      if (nsp && nst) ierr = 0;  // at least some files were read; so proceed...
      //TODO decide if we want to exit here or not when only partial file-list is read

      // need to merge stls at their collocated nodes (and cleanup tris)
      mergeCollocatedNodes();
      deleteTrisWithCollocatedNodes();
    }
    else if (token == "STL") {
      const string name = param->getString(iarg++);
      COUT1(" > about to read stl filename: " << name);
      ierr = addStl(name);
      // need to merge collocated nodes (and cleanup tris)
      mergeCollocatedNodes();
      deleteTrisWithCollocatedNodes();
    }
    else if (token == "MSH" || token == "CAS" || token == "FLUENT") {
      const string name = param->getString(iarg++);
      int force_byte_swap = -1;
      while (iarg != param->size()) {
        const string tok2 = MiscUtils::toUpperCase(param->getString(iarg++));
        if (tok2 == "BYTE_SWAP") {
          const string val = MiscUtils::toUpperCase(param->getString(iarg++));
          if (val == "TRUE") force_byte_swap = 1;
          else if (val == "FALSE") force_byte_swap = 0;
          else {
            CWARN(" ! unrecognized BYTE_SWAP value " << val << "; skipping");
          }
        }
      }
      COUT1(" > about to read msh filename: " << name);
      ierr = addMsh(name,force_byte_swap);
      deleteIsolatedNodes(); // Pointwise seemed to have a lot of these
    }
    else if (token == "PLY") {
      const string name = param->getString(iarg++);
      COUT1(" > about to read ply filename: " << name);
      ierr = addPly(name);
    }
    else if (token == "UGRID") {
      const string name = param->getString(iarg++);
      COUT1(" > about to read ugrid filename: " << name);
      if (iarg < param->size()) {
        const string bc_file = param->getString(iarg++);
        ierr = addUgrid(name,bc_file);
      }
      else {
        ierr = addUgrid(name);
      }
    }
    else if (token == "OBJ") {
      const string name = param->getString(iarg++);
      COUT1(" > about to read obj filename: " << name);
      ierr = addObj(name);
    }
    else if ((token == "BIN")||(token == "SBIN")) {
      const string name = param->getString(iarg++);
      COUT1(" > about to read surface binary filename \"" << name << "\"...");
      ierr = addBinary(name);
    }
    else if (token == "RESTART") {
      const string name = param->getString(iarg++);
      COUT1(" > about to read restart file: " << name);
      ierr = addRestartSurface(name);
    }
    else if ((token == "TEC")||(token == "TECPLOT")) {
      const string name = param->getString(iarg++);
      COUT1(" > about to read tecplot file: " << name);
      ierr = addTecplot(name);
    }
    else if (token == "PART") {
      const string name = param->getString(iarg++);
      COUT1(" > about to read surface from part file: " << name);
      ierr = addPartSurface(name);
    }
    else if (token == "GEOM") {
      const string geom_type = param->getString(iarg++);
      if (geom_type == "PLANE") ierr = constructPlane(iarg,param);
      else {
        CWARN("unrecognized GEOM surface type \"" << geom_type << "\"; skipping");
        ierr = -1;
      }
    }
    else if ((token == "SIMPLE_BOX")||(token == "BOX")) {
      // SIMPLE_BOX x0 x1 y0 y1 z0 z1...
      double x0[3];
      double x1[3];
      x0[0] = param->getDouble(iarg++);
      x1[0] = param->getDouble(iarg++);
      x0[1] = param->getDouble(iarg++);
      x1[1] = param->getDouble(iarg++);
      x0[2] = param->getDouble(iarg++);
      x1[2] = param->getDouble(iarg++);
      bool flip = false;
      while (iarg < param->size()) {
	token = MiscUtils::toUpperCase(param->getString(iarg++));
	if (token == "FLIP") {
	  flip = true;
	}
	else {
	  cout << "unrecognized BOX param: " << token << endl;
	}
      }
      ierr = addBox(x0,x1,flip);
    }
    else if ((token == "CIRCLE")) {
      // CIRCLE x y z nx ny nz r N
      double xc[3];
      double np[3];
      double r;
      int n=128;
      xc[0] = param->getDouble(iarg++);
      xc[1] = param->getDouble(iarg++);
      xc[2] = param->getDouble(iarg++);
      np[0] = param->getDouble(iarg++);
      np[1] = param->getDouble(iarg++);
      np[2] = param->getDouble(iarg++);
      r = param->getDouble(iarg++);

      if (iarg < param->size()) {
        if (param->getString(iarg) == "N") {
          ++iarg;
          n = param->getInt(iarg++);
        }
      }

      ierr = addCircle(xc,np,r,n);
    }
    else if ((token == "ANNULUS")) {
      // ANNULUS x y z nx ny nz r0 r1 N
      double xc[3];
      double np[3];
      double r0,r1;
      int n=128;
      xc[0] = param->getDouble(iarg++);
      xc[1] = param->getDouble(iarg++);
      xc[2] = param->getDouble(iarg++);
      np[0] = param->getDouble(iarg++);
      np[1] = param->getDouble(iarg++);
      np[2] = param->getDouble(iarg++);
      r0 = param->getDouble(iarg++);
      r1 = param->getDouble(iarg++);

      if (iarg < param->size()) {
        if (param->getString(iarg) == "N") {
          ++iarg;
          n = param->getInt(iarg++);
        }
      }

      ierr = addAnnulus(xc,np,r0,r1,n);
    }
    else if ((token == "SIMPLE_DISK")||(token == "DISK")||(token == "PIPE")||(token == "SIMPLE_PIPE")||(token == "CYLINDER")||(token == "CYL")) {
      // SIMPLE_DISK x0 y0 z0 x1 y1 z1 r...
      double x0[3];
      double x1[3];
      double r;
      int n = 128;
      token = MiscUtils::toUpperCase(param->getString(iarg));
      if (token == "X0") ++iarg;
      x0[0] = param->getDouble(iarg++);
      x0[1] = param->getDouble(iarg++);
      x0[2] = param->getDouble(iarg++);
      token = MiscUtils::toUpperCase(param->getString(iarg));
      if (token == "X1") ++iarg;
      x1[0] = param->getDouble(iarg++);
      x1[1] = param->getDouble(iarg++);
      x1[2] = param->getDouble(iarg++);
      token = MiscUtils::toUpperCase(param->getString(iarg));
      if ((token == "R")||(token == "RADIUS")) ++iarg;
      r = param->getDouble(iarg++);
      // additional optional tokens...
      bool b_flip = false;
      while (iarg < param->size()) {
	token = MiscUtils::toUpperCase(param->getString(iarg++));
	if (token == "N") {
          n = param->getInt(iarg++);
        }
	else if (token == "FLIP") {
	  b_flip = true;
	}
	else {
	  cout << "Warning: skipping unrecognized token: " << token << endl;
	}
      }
      ierr = addDisk(x0,x1,r,n,b_flip);
    }
    else if ((token == "SIMPLE_TCONE")||(token == "TCONE")) {
      // SIMPLE_TCONE x0 y0 z0 r0 x1 y1 z1 r1...
      double x0[3];
      double x1[3];
      double r0;
      double r1;
      int n=128;
      x0[0] = param->getDouble(iarg++);
      x0[1] = param->getDouble(iarg++);
      x0[2] = param->getDouble(iarg++);
      r0 = param->getDouble(iarg++);
      x1[0] = param->getDouble(iarg++);
      x1[1] = param->getDouble(iarg++);
      x1[2] = param->getDouble(iarg++);
      r1 = param->getDouble(iarg++);
      if (iarg < param->size()) {
        if (param->getString(iarg) == "N") {
          ++iarg;
          n = param->getInt(iarg++);
        }
      }

      ierr = addTcone(x0,x1,r0,r1,n);
    }
    else if ((token == "SIMPLE_ANNULUS")||(token == "ANNULUS")) {
      // ANNULUS x0 y0 z0 r00 r01 x1 y1 z1 r10 r11...
      double x0[3];
      double x1[3];
      double r0[2];
      double r1[2];
      int n = 128;
      x0[0] = param->getDouble(iarg++);
      x0[1] = param->getDouble(iarg++);
      x0[2] = param->getDouble(iarg++);
      r0[0] = param->getDouble(iarg++);
      r0[1] = param->getDouble(iarg++);
      x1[0] = param->getDouble(iarg++);
      x1[1] = param->getDouble(iarg++);
      x1[2] = param->getDouble(iarg++);
      r1[0] = param->getDouble(iarg++);
      r1[1] = param->getDouble(iarg++);
      if (iarg < param->size()) {
        if (param->getString(iarg) == "N") {
          ++iarg;
          n = param->getInt(iarg++);
        }
      }

      ierr = addAnnularTcone(x0,x1,r0[0],r0[1],r1[0],r1[1],n);
    }
    else if ((token == "SIMPLE_SPHERE")||(token == "SPHERE")) {
      // SIMPLE_SPHERE x y z r...
      double x[3];
      double r;
      x[0] = param->getDouble(iarg++);
      x[1] = param->getDouble(iarg++);
      x[2] = param->getDouble(iarg++);
      r = param->getDouble(iarg++);
      int n = 4; // default sphere edge segment count 4
      if (iarg < param->size()) {
        token = param->getString(iarg);
        if ((token == "N")||(token == "n")) {
          ++iarg;
          n = param->getInt(iarg++);
        }
      }
      ierr = addSphere(x,r,n,false); // last arg is flip: default is fluid inside
    }
    else if (token == "HEMISPHERE") {
      // HEMISPHERE x y z nx ny nz r...
      double xp[3];
      double np[3];
      double rp;
      xp[0] = param->getDouble(iarg++);
      xp[1] = param->getDouble(iarg++);
      xp[2] = param->getDouble(iarg++);
      np[0] = param->getDouble(iarg++);
      np[1] = param->getDouble(iarg++);
      np[2] = param->getDouble(iarg++);
      rp = param->getDouble(iarg++);
      int ntheta = 128; // default hemisphere edge segmentation: 128
      if (iarg < param->size()) {
        token = param->getString(iarg);
        if ((token == "N")||(token == "NTHETA")||(token == "n")||(token == "ntheta")) {
          ++iarg;
          ntheta = param->getInt(iarg++);
        }
      }
      ierr = addHemisphere(xp,np,rp,ntheta,false); // last arg is flip: default is fluid inside
    }
    else if (token == "PISTON_FF") {
      // PISTON_FF dr0 dz0 dr1 dz1
      double dr0 = param->getDouble(iarg++);
      double dz0 = param->getDouble(iarg++);
      double dr1 = param->getDouble(iarg++);
      double dz1 = param->getDouble(iarg++);
      ierr = addPistonFF(dr0,dz0,dr1,dz1);
    }
    else if (token == "LIFTED") {
      // LIFTED DN <double>
      bool b_dn = false;
      double dn;
      vector<int> subzone_indices;
      while (iarg < param->size()) {
        token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (parseSpecifiedSubzones(subzone_indices,token,param,iarg)) {
          // no-op
        }
        else if (token == "DN") {
          b_dn = true;
          dn = param->getDouble(iarg++);
        }
        else {
          cout << "Error: unrecognized LIFTED_SURFACE token: " << token << endl;
          ierr = -1;
        }
      }
      if (!b_dn)
        ierr = -1;
      if (ierr == 0) {
        // the lifted surface acts on the flagged subzones...
        sz_flag.resize(nsz);
        if (subzone_indices.empty()) {
          sz_flag.setAll(1);
        }
        else {
          sz_flag.setAll(0);
          for (vector<int>::iterator it=subzone_indices.begin(); it!=subzone_indices.end(); ++it)
            sz_flag[*it] = 1;
        }
        ierr = addLiftedSurface(dn);
      }
    }
    else if (token == "LOFTED") {
      vector<string> profileVec;
      bool b_nl = false;
      int nl; // number of points in the lofted direction. The other direction will be set to be approximately isotropic
      int iarg = 1;
      while (iarg < param->size()) {
        token = MiscUtils::toUpperCase(param->getString(iarg++));
        if ((token == "PROFILE")||(token == "PROFILES")) {
          // assume the rest of the tokens are filenames...
          while (iarg < param->size()) {
            token = param->getString(iarg++);
            cout << " > adding profile: " << token << endl;
            profileVec.push_back(token);
          }
        }
        else if (token == "NL") {
          b_nl = true;
          nl = param->getInt(iarg++);
          cout << " > NL: " << nl << endl;
        }
        else {
          cout << "WARNING: skipping unrecognized token: " << token << endl;
        }
      }
      if (profileVec.size() < 2) {
        cout << "Error: LOFTED requires atleast 2 PROFILE files" << endl;
        ierr = -1;
      }
      if (!b_nl) {
        cout << "Error: point count in lofted direction NL not specified" << endl;
        ierr = -1;
      }
      if (ierr == 0) {
        ierr = addLofted(profileVec,nl);
      }
    }
    else if (token == "REVOLVE_SPLINE") {
      vector<double> dVec;
      double axis[3] = { 1.0, 0.0, 0.0 };
      double origin[3] = { 0.0, 0.0, 0.0 };
      int ntheta = 128;
      int iarg = 1;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));
	if (token == "PTS") {
	  // PTS can be followed by a series of points, or a filename.
	  const string filename = param->getString(iarg++);
	  ierr = GeomUtils::readXY(dVec,filename);
	}
	else if (token == "AXIS") {
	  FOR_I3 axis[i] = param->getDouble(iarg++);
	}
	else if (token == "ORIGIN") {
	  FOR_I3 origin[i] = param->getDouble(iarg++);
	}
	else if (token == "NTHETA") {
          ntheta = param->getInt(iarg++);
        }
        else {
          cout << "WARNING: skipping unrecognized token: " << token << endl;
        }
      }
      if (ierr == 0) {
	ierr = addRevolveSpline(dVec,axis,origin,ntheta);
      }
    }
    else if (token == "EXTRUDE_SPLINE") {
      vector<double> dVec;
      double x0[3] = { 0.0, 0.0, 0.0 };
      double x1[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
      int nspan = -1;
      int iarg = 1;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));
	if (token == "PTS") {
	  // PTS can be followed by a series of points, or a filename.
	  const string filename = param->getString(iarg++);
	  ierr = GeomUtils::readXY(dVec,filename);
	}
	else if ((token == "X0")||(token == "START")) {
	  FOR_I3 x0[i] = param->getDouble(iarg++);
	}
	else if ((token == "X1")||(token == "END")) {
	  FOR_I3 x1[i] = param->getDouble(iarg++);
        }
	else if ((token == "NSPAN")||(token == "N_SPAN")) {
          nspan = param->getInt(iarg++);
        }
        else {
          cout << "WARNING: skipping unrecognized EXTRUDE_SPLINE token: " << token << endl;
        }
      }
      if (x1[0] == HUGE_VAL) {
        cout << "WARNING: EXTRUDE_SPLINE missing X1. Skipping..." << endl;
        ierr = 1;
      }
      if (nspan <= 0) {
        cout << "WARNING: EXTRUDE_SPLINE needs a positive NSPAN. Skipping..." << endl;
        ierr = 1;
      }
      if (ierr == 0) {
	ierr = addExtrudeSpline(dVec,x0,x1,nspan);
      }
    }
    else if (token == "OFFSET") {
      // is this used anymore? I think it got replaced with lifted surface above.
      // SURF OFFSET DN <double>
      token = param->getString(iarg++);
      assert((token == "DN")||(token == "dn"));
      double dn = param->getDouble(iarg++);
      ierr = addOffset(dn);
    }
    else if ((token == "SIMPLE_REVOLVE")||(token == "REVOLVE")) {

      // example:
      // SIMPLE_REVOLVE AXIS <nx> <ny> <nz> START <x> <y> <z> PROFILE <filename> N_THETA <n> CLOSE MAX_AR <aspect-ratio>
      // e.g. SURF SIMPLE_REVOLVE AXIS 1 0 0 START 0 0 0 PROFILE 2d.dat N_THETA 128 MAX_AR 2.0
      // where 2d.dat looks like:
      // # 2d points
      // 0 0.5
      // 0.5 0.5
      // 1.0 1.0
      // 3.0 1.0

      // revolve reads profile from file now
      double origin[3] = { 0.0, 0.0, 0.0 };  // starting point
      double axis[3];
      bool b_axis = false;
      double (*ir_profile)[2] = NULL;
      int n_vals = 0;
      int n_theta = 0;
      bool b_cap = false;
      double maxAR = -1.0;

      while (iarg != param->size()) {
        token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (token == "AXIS") {
          b_axis = true;
          FOR_I3 axis[i] = param->getDouble(iarg++);
        }
        else if ((token == "START")||(token == "ORIGIN")) {
          FOR_I3 origin[i] = param->getDouble(iarg++);
        }
        else if (token == "PROFILE") {
          n_vals = MiscUtils::read2DAsciiTable(ir_profile,param->getString(iarg++));
        }
        else if ((token == "N_THETA")||(token == "NTHETA")) {
          n_theta = param->getInt(iarg++);
        }
        else if (token == "CLOSE") {
          b_cap = true;
        }
        else if (token == "MAX_AR") {
          maxAR = param->getDouble(iarg++);
        }
        else {
          CWARN("unrecognized REVOLVE parameter: " << token << "; skipping");
        }
      }

      if (b_axis && n_vals && (n_theta > 2)) {
        ierr = addRevolve(axis,origin,ir_profile,n_vals,n_theta,b_cap,maxAR);
      }
      else {
        CWARN("not all required REVOLVE parameters were specified: sample usage" <<
              "\nSURF REVOLVE AXIS <nx> <ny> <nz> START <x> <y> <z> PROFILE <filename> N_THETA <n> CLOSE MAX_AR <aspect-ratio>");
        ierr = -1;
      }
      DELETE(ir_profile);
    }
    else if ((token == "EXTRUDE")) {
      double x0[3];  // starting point
      double x1[3];  // ending point
      bool b_start = false;
      bool b_end = false;
      bool b_n = false;
      bool b_d = false;
      double (*xyz_profile)[3] = NULL;
      int n_vals = 0;
      int n_span = 1;
      double normal[3];
      double dist;
      bool closed = false;
      int skip = 1;

      while (iarg != param->size()) {
        const string tok = MiscUtils::toUpperCase(param->getString(iarg++));
        if (tok == "START") {
          b_start = true;
          FOR_I3 x0[i] = param->getDouble(iarg++);
        }
        else if (tok == "END") {
          b_end = true;
          FOR_I3 x1[i] = param->getDouble(iarg++);
        }
        else if (tok == "SKIP") {
          skip = param->getInt(iarg++);
        }
        else if (tok == "N") {
          b_n = true;
          FOR_I3 normal[i] = param->getDouble(iarg++);
          NORMALIZE(normal);
        }
        else if (tok == "D") {
          b_d = true;
        dist = param->getDouble(iarg++);
        }
        else if (tok == "PROFILE_XY" || tok == "PROFILE_XZ" || tok == "PROFILE_YZ") {
          if (xyz_profile != NULL) {
            WUI(WARN,"cannot specify multiple profiles; skipping");
            continue;
          }
          else {
            double (*profile2d)[2] = NULL;
            n_vals = MiscUtils::read2DAsciiTable(profile2d,param->getString(iarg++));
            int a,b,c;
            if (tok == "PROFILE_XY") {
              a = 0;
              b = 1;
              c = 2;
            }
            else if (tok == "PROFILE_XZ") {
              a = 0;
              b = 2;
              c = 1;
            }
            else {
              a = 1;
              b = 2;
              c = 0;
            }
            xyz_profile = new double[n_vals][3];
            for (int ival=0; ival<n_vals; ++ival) {
              xyz_profile[ival][a] = profile2d[ival][0];
              xyz_profile[ival][b] = profile2d[ival][1];
              xyz_profile[ival][c] = 0.0;
            }
          }
        }
        else if (tok == "PROFILE_3D") {
          if (xyz_profile != NULL) {
            WUI(WARN,"cannot specify multiple profiles; skipping");
          }
          else {
            n_vals = MiscUtils::read3DAsciiTable(xyz_profile,param->getString(iarg++));
          }
        }
        else if (tok == "CLOSED") {
          closed = true;
        }
        else if (tok == "N_SPAN" || tok == "NSPAN") {
          n_span = param->getInt(iarg++);
        }
        else {
          CWARN("unrecognized EXTRUDE parameter: " << tok << "; skipping");
        }
      }

      if (xyz_profile != NULL) {
        bool axis_set = false;
        if (b_start && b_end) {
          axis_set = true;
          COUT2(" > using axis determined from START/END parameters");
        }
        else if (b_n && b_d) {
          axis_set = true;
          COUT2(" > using axis determined from N/D parameters");
          FOR_I3 x0[i] = xyz_profile[0][i];
          FOR_I3 x1[i] = x0[i] + dist*normal[i];
        }
        else {
          WUI(WARN,"extrusion axis not properly set; skipping");
          ierr = -1;
        }

        if (skip > 1) {
          double (*xyz_profile0)[3] = xyz_profile;
          xyz_profile = new double[n_vals/skip][3];
          for (int i_new=0,i_old=0; i_old < n_vals; ++i_new,i_old+=skip) {
            FOR_I3 xyz_profile[i_new][i] = xyz_profile0[i_old][i];
          }
          n_vals /= skip;
          DELETE(xyz_profile0);
        }

        if (axis_set) {
          if (n_vals && n_span) {
            ierr = addExtrudedProfile(x0,x1,xyz_profile,n_vals,n_span,closed);
          }
          else {
            CWARN("not all EXTRUDE parameters were correctly specified; skipping");
            ierr = -1;
          }
        }
        DELETE(xyz_profile);
      }
      else {
        WUI(WARN,"profile to extrude was not properly parsed; skipping");
        ierr = -1;
      }
    }
    else if (token == "NACA") {
      bool parsed = false;
      bool symmetric = false;
      int xx = -1;
      int m = -1;
      int p = -1;
      int n_chord = 100;
      int n_span = -1;
      double chord = 1.0;
      double span = 0.10;
      bool b_sharp = false;
      const string naca_type = param->getString(iarg++);
      if (naca_type.length() != 4) {
        CWARN("only 4-digit airfoils are currently supported; skipping");
        ierr = -1;
      }
      else {
        const string first_two_digits = naca_type.substr(0,2);
        if (first_two_digits == "00") {
          symmetric = true;
          const string last_two_digits = naca_type.substr(2,2);
          xx = atoi(last_two_digits.c_str());
          if (xx) parsed = true;  // xx is non-zero so string-to-int conversion was successful
        }
        else {
          const string first_digit = naca_type.substr(0,1);
          const string second_digit = naca_type.substr(1,1);
          const string last_two_digits = naca_type.substr(2,2);
          m = atoi(first_digit.c_str());
          p = atoi(second_digit.c_str());
          xx = atoi(last_two_digits.c_str());
          if (m && p && xx) parsed = true;  // if non-zero string-to-int conversion was successful
        }
      }

      while (iarg != param->size()) {
        const string tok = param->getString(iarg++);
        if (tok == "N_CHORD") {
          n_chord = param->getInt(iarg++);
        }
        else if (tok == "N_SPAN") {
          n_span = param->getInt(iarg++);
        }
        else if (tok == "CHORD") {
          chord = param->getDouble(iarg++);
        }
        else if (tok == "SPAN") {
          span = param->getDouble(iarg++);
        }
        else if (tok == "SHARP") {
          b_sharp = true;
        }
        else {
          CWARN("unrecognized NACA parameter \"" << tok << "\"; skipping");
          ierr = -1;
        }
      }

      if (parsed) {
        if (symmetric) {
          ierr = addNaca00Airfoil(xx,n_chord,chord,n_span,span,b_sharp);
        }
        else {
          ierr = addNaca4DigitAirfoil(m,p,xx,n_chord,chord,n_span,span);
        }
      }
      else {
        CWARN("NACA parameters did not sufficiently describe an airfoil surface");
        ierr = -1;
      }
    }
    else if (token == "GRID") {
      // -----------------------------------------------------------------------
      // adds a grid of turbulence sphere primitives that can be used to trip
      // turbulence at a specified boundary...
      //
      // SURF GRID NAME = annular_grid D 0.0004 L 0.001 X -0.09 0 0 N 1 0 0 SUBZONES 8,7
      // SURF GRID NAME = annular_grid D 0.0004 L 0.001 X -0.09 0 0 N 1 0 0
      // -----------------------------------------------------------------------
      st_flag.setLength(nst);
      st_flag.setAll(1); // consider all tris, unless user specifies participating subzones
      bool b_diam = false;
      double diam;
      bool b_spacing = false;
      double spacing;
      bool b_xc = false;
      double xc[3];
      bool b_normal = false;
      double normal[3];
      int iarg = 1;
      bool b_name = false;
      string name;
      while (iarg < param->size()) {
        string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (token == "NAME") {
          b_name = true;
          name = param->getString(iarg++);
        }
        else if ((token == "D")||(token == "DIAM")||(token == "DIAMETER")) {
          b_diam = true;
          diam = param->getDouble(iarg++);
        }
        else if ((token == "L")||(token == "SPACING")) {
          b_spacing = true;
          spacing = param->getDouble(iarg++);
        }
        else if ((token == "X")||(token == "XC")) {
          b_xc = true;
          FOR_I3 xc[i] = param->getDouble(iarg++);
        }
        else if ((token == "N")||(token == "NORMAL")) {
          b_normal = true;
          FOR_I3 normal[i] = param->getDouble(iarg++);
        }
        else if (token == "SUBZONES") {
          st_flag.setAll(0);
          vector<int> subzonesVec;
          MiscUtils::splitCsv(subzonesVec,param->getString(iarg++));
          sz_flag.setLength(nsz);
          sz_flag.setAll(0);
          for (int ii = 0, limit=subzonesVec.size(); ii < limit; ++ii) {
            const int isz = subzonesVec[ii]; assert((isz >= 0)&&(isz < nsz));
            sz_flag[isz] = 1;
          }
          for (int ist = 0; ist < nst; ++ist)
          if (sz_flag[szost[ist]] == 1)
          st_flag[ist] = 1;
        }
        else {
          WUI(WARN,"unrecognized SURF GRID token: " << token << "; skipping. Example:" <<
              "\nSURF GRID NAME=gr D=0.005 L=0.01 X=0.2 0 0 N=1 0 0 SUBZONES 8,9");
          ierr = -1;
        }
      }
      if (ierr == 0) {
        if (!(b_diam && b_spacing && b_xc && b_normal && b_name)) {
          WUI(WARN,"SURF GRID missing information; skipping. Example:" <<
              "\nSURF GRID NAME=gr D=0.005 L=0.01 X=0.2 0 0 N=1 0 0 SUBZONES 8,9");
          ierr = -1;
        }
        else {
          ierr = addGridForFlaggedTris(diam,spacing,xc,normal,name);
          WUI(INFO,"SURF GRID added successfully");
        }
      }
    }
    else {
      CERR("unrecognized SURF token: " << token);
    }

    if (ierr != 0) return ierr;

    // the subzone indexing (but not subzone data) is always available.
    // it is initialized to the zone indexing: i.e. 1 subzone per zone. Note that
    // if we are adding some new surface, this will reset ALL subzoning, even
    // if it existed before the add...
    setSubzonesToZones();

    return 0;  // success
  }

  // routines helpful for parsing "ZONES" paramaters...

  bool parseSpecifiedZones(vector<int>& zone_indices,const string& token,Param * param,int& iarg) const {

    bool parsed = false;  // indicates whether token was used here
    if ((token == "ZONE_NAME") || (token == "ZONE_NAMES") || (token == "ZONE") || (token == "ZONES") || (token == "ZNS")) {
      parsed = true;
      const string zonesCsv = param->getString(iarg++);
      vector<string> zonesVec;
      MiscUtils::splitCsv(zonesVec,zonesCsv);
      for (vector<string>::iterator it=zonesVec.begin(); it!=zonesVec.end(); ++it) {
        const int nzn=zoneVec.size();
        int zone_index;
        for (zone_index=0; zone_index < nzn; ++zone_index) {
          if (*it == zoneVec[zone_index].getName()) break;
        }
        if (zone_index != int(zoneVec.size())) {
          zone_indices.push_back(zone_index);
        }
        else COUT2(" > skipping invalid zone name: " << *it);
      }
    }
    else if ((token == "ZONE_INDICES") || (token == "ZONE_IDS") || (token == "ZNIDS")) {
      parsed = true;
      const string zonesCsv = param->getString(iarg++);
      vector<string> zonesVec;
      MiscUtils::splitCsv(zonesVec,zonesCsv);
      const int nzn=zoneVec.size();
      for (vector<string>::iterator it=zonesVec.begin(); it!=zonesVec.end(); ++it) {
        //const int zone_index = stoi(zonesVec[i]);
        const int zone_index = atoi((*it).c_str());
        if (zone_index >= 0 && zone_index < nzn) {
          zone_indices.push_back(zone_index);
        }
        else COUT2(" > skipping invalid zone index: " << zone_index);
      }
    }
    else if (token == "ALL") {
      parsed = true;
      zone_indices.clear();  // to avoid redundant if any of above specified
      for (int izn=0,nzn=zoneVec.size(); izn < nzn; ++izn) zone_indices.push_back(izn);

    }
    return parsed;
  }

  bool parseSpecifiedSubzones(vector<int>& subzone_indices,const string& token,Param * param,int& iarg) const {

    bool parsed = false;  // indicates whether token was used here
    if ((token == "SUBZONES") || (token == "SUBZONE") || (token == "SZNS") || (token == "SZN")) {
      parsed = true;
      const string subzonesCsv = param->getString(iarg++);
      vector<int> subzonesVec;
      MiscUtils::splitCsv(subzonesVec,subzonesCsv);
      for (vector<int>::iterator it=subzonesVec.begin(); it!=subzonesVec.end(); ++it) {
        const int isz = *it;
        if ((isz >= 0)&&(isz < nsz)) subzone_indices.push_back(isz);
        else cout << "Warning: skipping invalid subzone index " << isz << " (nsz=" << nsz << ")" << endl;
      }
    }
    else if ((token == "ZONE_NAME") || (token == "ZONE_NAMES") || (token == "ZONE") || (token == "ZONES") || (token == "ZNS")) {
      parsed = true;
      const string zonesCsv = param->getString(iarg++);
      vector<string> zonesVec;
      MiscUtils::splitCsv(zonesVec,zonesCsv);
      for (vector<string>::iterator it=zonesVec.begin(); it!=zonesVec.end(); ++it) {
        const int nzn=zoneVec.size();
        int zone_index;
        for (zone_index=0; zone_index < nzn; ++zone_index) {
          if (*it == zoneVec[zone_index].getName()) break;
        }
        if (zone_index != int(zoneVec.size())) {
          for (int isz=szozn_i[zone_index]; isz<szozn_i[zone_index+1]; isz++) subzone_indices.push_back(isz);
        }
        else COUT2(" > skipping invalid zone name: " << *it);
      }
    }
    else if ((token == "ZONE_INDICES") || (token == "ZONE_IDS") || (token == "ZNIDS")) {
      parsed = true;
      const string zonesCsv = param->getString(iarg++);
      vector<string> zonesVec;
      MiscUtils::splitCsv(zonesVec,zonesCsv);
      const int nzn=zoneVec.size();
      for (vector<string>::iterator it=zonesVec.begin(); it!=zonesVec.end(); ++it) {
        //const int zone_index = stoi(zonesVec[i]);
        const int zone_index = atoi((*it).c_str());
        if (zone_index >= 0 && zone_index < nzn) {
          for (int isz=szozn_i[zone_index]; isz<szozn_i[zone_index+1]; isz++) subzone_indices.push_back(isz);
        }
        else COUT2(" > skipping invalid zone index: " << zone_index);
      }
    }
    else if (token == "ALL") {
      parsed = true;
      subzone_indices.clear();  // to avoid redundant if any of above specified
      for (int isz=0; isz < szozn_i[zoneVec.size()]; ++isz) subzone_indices.push_back(isz);

    }
    return parsed;
  }

  bool parseSpecifiedEdges(vector<int>& edge_indices,const string token,Param * param,int& iarg) {
      bool parsed = false;  // indicates whether token was used here
      if ((token == "ZONES") || (token == "SUBZONES") || (token == "LOOPS")) {
        parsed = true;
        const string zonesCsv = param->getString(iarg++);
        vector<string> zonesVec;
        MiscUtils::splitCsv(zonesVec,zonesCsv);
        for (vector<string>::iterator it=zonesVec.begin(); it!=zonesVec.end(); ++it) {
          edge_indices.push_back(atoi((*it).c_str()));
        }
      }
      else if (token == "ALL") {
        parsed = true;
        edge_indices.clear();  // parsed but empty means all

      }
      return parsed;
    }

  bool flagZones(const string& csv) {
    for (int izn = 0, size = zoneVec.size(); izn < size; ++izn)
      zoneVec[izn].flag =  0;
    vector<string> tokens;
    MiscUtils::splitCsv(tokens,csv);
    bool worked = true;
    for (int ii = 0; ii < tokens.size(); ++ii) {
      // allow for a wildcard...
      bool found = false;
      for (int izn = 0, size = zoneVec.size(); izn < size; ++izn) {
	if (MiscUtils::strcmp_wildcard(zoneVec[izn].getName(),tokens[ii])) {
	  zoneVec[izn].flag = 1;
	  found = true;
	}
      }
      if (!found) {
	// maybe it was a zone index...
	int izn;
	if (from_string<int>(izn,tokens[ii],std::dec)&&(izn >= 0)&&(izn < zoneVec.size())) {
	  zoneVec[izn].flag = 1;
	}
	else {
	  WUI(WARN,"no matching zone found for " << tokens[ii]);
	  worked = false;
	}
      }
    }
    return worked;
  }

  // ============================================
  // these routines in SimpleSurface_teost.cpp...
  // ============================================

public:
  void buildTeost();
  void clearTeost();
  void ensureTeost();
  bool gotTeost() const;
  int getTriNbrDataFull(int& ist_nbr,int& i_nbr,int& orient_nbr,const int ist,const int i) const;
  bool getTriNbrData(int& ist_nbr,int& i_nbr,int& orient_nbr,const int ist,const int i) const;
  void setTriNbrData(const int ist_nbr,const int i_nbr,const int orient_nbr,const int ist,const int i);
  void setTriNbrOpen(const int ist,const int i);
  void setTriNbrMulti(const int ime,const int orient_ime,const int ist,const int i);
  void setTriNbrToSameAs(const int ist_ref,const int i_ref,const int ist,const int i,const bool b_flip_orient = false);
  void getTriNbrs(int (&ist_nbr)[3],const int ist) const;
  bool getAlignedTriNbr(int& ist_nbr,const int ist,const int i) const;
  bool isNbrAligned(const int ist,const int i) const;
  bool isEdgeOpen(const int ist,const int i) const;
  bool isEdgeMulti(const int ist,const int i) const;
  bool isEdgeMulti(int& ime,int& orient_ime,const int ist,const int i) const;
private:
  void setOpenEdgeGroup(const int igr,const int ist,const int i);
public:
  bool getOpenEdgeGroup(int& igr,const int ist,const int i) const;
  bool isEdgeZoneBoundary(const int ist, const int i) const;
  bool isEdgeNmeBoundary(const int ist, const int i) const;
  bool isEdgeSubzoneBoundary(const int ist, const int i) const;
  bool isEdgeSelectionBoundary(const int ist, const int i) const;
  bool isEdgeMisaligned(const int ist, const int i) const;
  // bool isEdgeMetalBoundary(const int ist, const int i) const;
  // bool isEdgeFluidBoundary(const int ist, const int i) const;
  double edgeCreaseCos(const int ist,const int i) const;
  double edgeCreaseCos(const int ist0,const int ist1, const bool aligned) const;
  bool isEdgeCrease(const int ist,const int i) const;
  bool isEdgeCrease(const int ist0,const int ist1, const bool aligned) const;
  bool isEdgeCrease(const double unit_n[3],const int ist1, const bool aligned) const;
  bool isEdgeFeature(const int ist, const int i) const;
  void multiEdgesToDynamicEdges();
  void flagTrisTouchingMultiEdges();
  void flagTrisTouchingMultiEdges(set<int>& specificMultiEdges);
  void flagTrisTouchingOpenEdges();
  void flagNmeTris(const vector<int>& zone_indices,const int nme_bits,const bool b_add);
  void flagTrisTouchingFlaggedTris();
  void flagTrisTouchingFlaggedTrisInSameZone();
  void flagTrisWithAreaLessThan(const double area);
  void flagTrisWithSubzoneAreaLessThan(const double area);

  // these new selection functions allow you to add to existing flagging...
  void flagTrisInZones(const vector<int>& zone_indices,const bool b_add);
  void flagTrisTouchingZones(const vector<int>& zone_indices,const bool b_add);
  void countFlaggedTris();
  void invertFlaggedTris();

  // ============================================
  // these routines in SimpleSurface_stl.cpp...
  // ============================================

  int addStl(const string& filename);  // returns 0 on success...
  void writeSelectedTrisToStl(const string& filename, const bool single_zone);
  void writeSelectedTrisToStlSolid(FILE * fp, const string& solidname, IntFlag& my_st_flag) const;
  void writeSelectedTrisToStlSolid(FILE * fp, const string& solidname) const;

  // ============================================
  // these routines in SimpleSurface_msh.cpp...
  // ============================================

  int addMsh(const string& filename,const int force_byte_swap);  // returns 0 on success...

  // ============================================
  // these routines in SimpleSurface_ply.cpp...
  // ============================================

  int addPly(const string& filename);  // returns 0 on success...

  // ============================================
  // these routines in SimpleSurface_obj.cpp...
  // ============================================

  int addObj(const string& filename);  // returns 0 on success...

  // ============================================
  // these routines in SimpleSurface_ugrid.cpp...
  // ============================================

  int addUgrid(const string& filename,const string& bc_filename = "");  // returns 0 on success...


  // ============================================
  // SimpleSurface_primitives.cpp
  // ============================================

  // this version produces all spheres, including those that intersect the boundary...
  int addGridForFlaggedTris(const double diam,const double spacing,const double xc[3],const double normal[3],const string& zone_name);

  int addPlane(const double xp[3],const double np[3],const double width,const double height,const int _nx, const int _ny);
  int addCircle(const double xp[3],const double np[3],const double r,const int n);
  int addAnnulus(const double xp[3],const double np[3],const double r0,const double r1,const int n);
  int addBox(const double x0[3],const double x1[3],const bool flip = false);
  int addDisk(const double x0[3],const double x1[3],const double r,const int n,const bool b_flip = false);
  int addTcone(const double x0[3],const double x1[3],const double rad0,const double rad1,const int n);
  int addAnnularTcone(const double x0[3],const double x1[3],const double rad00,const double rad01,const double rad10,const double rad11,const int n);
  int addRevolve(const double rAxis[3],const double x0[3],double (*ir_vals)[2],const int n_vals,const int _n_theta, bool capSurf, const double maxAR);
  int addExtrudedProfile(const double x0[3], const double x1[3], double (*xyz_vals)[3], const int n_vals, const int _n_span, const bool closed);

  // these should be put into a geom namespace...
  void facetCircleToPoint(int (* spost)[3], int * znost, const int indexPt, const int indexCircle, const int st0, const int zoneId, const int n, const bool flip);
  void facetCircleToCircle(int (* spost)[3], int * znost, const int indexC0, const int indexC1, const int st0, const int zoneId, const int n, const bool closed,const bool flip);

  int addSphere(const double xc0[3],const double r,const int n,const bool flip);
  int addSphereWithBl(const double xc0[3],const double r,const double dn,const double dt_target,const int nlayers);
  int addHemisphere(const double xp[3],const double np[3],const double rp,const int ntheta,const bool flip);
  int addPistonFF(const double dr0,const double dz0,const double dr1,const double dz1);
  int addLiftedSurface(const double dn);
  int addOffset(const double dn);
  int alignNormalsUsingBloatedSurface(const double delta,const bool seed_with_bloat);

private:

  int constructPlane(int& iarg,Param * param);

  int addLofted(const vector<string>& profileVec,const int nr);
  int addRevolveSpline(vector<double>& dVec,const double axis[3],const double origin[3],const int ntheta);
  int addExtrudeSpline(vector<double>& dVec,const double x0[3],const double x1[3],const int nspan);

public:

  void interfaceFlaggedZones(const double xc[3],const double axis[3]);
  void makeFF(const double xc[3],const double axis[3],const double dn,const double ds); // dn == standoff distance from moving surface, ds == target segment length

private:

  int addRevolve(vector<double>& segVec,const double xc[3],const double axis[3],const int ntheta);

private:

  // ============================================
  // SimpleSurface_naca.cpp
  // ============================================

  int addNaca00Airfoil(const int xx,const int _n_panels,const double chord,const int _n_span,const double span,const bool b_sharp=false);
  int addNaca4DigitAirfoil(const int m, const int p,const int xx,const int _n_panels,const double chord,const int _n_span,const double span);

public:

  // ============================================
  // these routines in SimpleSurface_intersect.cpp
  // ============================================

  void intersectFlaggedSubzones(const string& mode);

private:

  void processTriTriIntersection(map<const pair<int,int>,int>& edgeMap,vector<pair<int,int> >& stoedVec,vector<IntersectionData>& intersectionVec,
				 const int * const idata,const double * const ddata,const int ist0,const int ist1);

  // ============================================
  // these routines in SimpleSurface_zipper.cpp
  // ============================================

public:

  void zipOpenEdges(const vector<int>& group_indices,const double edge_factor,const double delta_max);
  void zipOpenEdges(const double edge_factor,const double delta_max);

private:

  class NodeNodeMerge {
  public:
    bool active;
    int isp0,isp1;
    double d2;
    NodeNodeMerge() {
      active = false;
    }
    NodeNodeMerge(const int isp0_,const int isp1_,const double d2_) {
      active = true;
      isp0 = isp0_;
      isp1 = isp1_;
      d2   = d2_;
    }
  };

  class EdgeNodeMerge {
  public:
    bool active;
    int ioe,isp;
    double t,d2;
    EdgeNodeMerge() {
      active = false;
    }
    EdgeNodeMerge(const int ioe_,const int isp_,const double t_,const double d2_) {
      active = true;
      ioe = ioe_;
      isp = isp_;
      t   = t_;
      d2  = d2_;
    }
    bool operator < (const EdgeNodeMerge& other) const {
      return (ioe < other.ioe) || ((ioe == other.ioe)&&(t < other.t));
    }
  };

  void zipOpenEdges(vector<uint8>& zipTeostVec,const double edge_factor,const double delta_max); // zip up the nodes and edges of these flagged tri-edges
  int setZipMergeVecs(vector<NodeNodeMerge>& nodeNodeMergeVec,vector<EdgeNodeMerge>& edgeNodeMergeVec,vector<uint8>& zipTeostVec,const double edge_factor,const double delta_max);
  void splitEdgesAndAddNodeNodeMerges(vector<NodeNodeMerge>& nodeNodeMergeVec,vector<EdgeNodeMerge>& edgeNodeMergeVec,vector<uint8>& zipTeostVec);
  void mergeNodes(vector<NodeNodeMerge>& nodeNodeMergeVec);

  void buildZipMergeVecs(vector<pair<int,int> >& nodeNodeMergeVec,vector<pair<pair<int,double>,int> >& edgeNodeMergeVec,const vector<uint8>& zipTeostVec);
  void splitZipEdgesAndCreateNodeNodeMerges(vector<pair<int,int> >& nodeNodeMergeVec,vector<pair<pair<int,double>,int> >& edgeNodeMergeVec,const vector<uint8>& zipTeostVec);
  void compressZipEdgeNodes(const vector<pair<int,int> >& nodeNodeMergeVec);

public:

  // ============================================
  // these routines in SimpleSurface_sbin.cpp...
  // ============================================

  int addBinary(const string& filename);  // returns 0 on success...
  void writeBinary(const string& filename) const;
  void writeSelectedTrisToBinary(const string& filename,const IntFlag& st_flag) const;
  void writeSubzonesToBinary(const string& filename,const vector<int>& subzonesVec);
  int addPartSurface(const string& filename);  // returns 0 on success...

  // ============================================
  // these routines in SimpleSurface_restart.cpp...
  // ============================================

  int addRestartSurface(const string& filename);  // returns 0 on success...

  // ============================================
  // these routines in SimpleSurface_tecplot.cpp...
  // ============================================

  void writeTecplot(const string& filename) const;
  void writeFlaggedZonesToTecplot(const string& filename);
  void writeSelectedFacesByZoneToTecplot(const string& filename, const IntFlag& st_flag,const bool b_flag=false) const;
  void writeSelectedFacesByZoneToTecplot(const string& filename, const IntFlag& st_flag,const double x0[3]) const;
  void writeSpDataToTecplot(const string& filename,const double *sp_buf);
  void writeSpDataToTecplot(const string& filename,const double (*sp_buf3)[3]);

  int addTecplot(const string& filename);

  // ============================================
  // these routines in SimpleSurface_input_files.cpp...
  // ============================================

  void writeStitchInputFile(const string& filename, const bool withHelp) const;
  vector<string> buildStitchInputFile(const bool withHelp) const;
  void writeCharlesInputFile(const string& filename, const int eos, const bool withHelp) const;
  vector<string> buildCharlesInputFile(const int eos,const bool withHelp) const;
  void buildIdealGasInputOptions(vector<string>& input_file,const bool withHelp) const;
  void buildPremixedInputOptions(vector<string>& input_file,const bool withHelp) const;
  void buildNonPremixedInputOptions(vector<string>& input_file,const bool withHelp) const;

  void writeFluxProbesCylinder(const string& name,const int isz,const int isz_prev,const int isz_next,const bool b_write_image,const double up[3]);
  void suggestMultiFluxProbeFlaggedZones();

  // ============================================
  // these routines in SimpleSurface_data.cpp...
  // ============================================

  double zoneArea(const int izn) {
    // need to rethink zone geometry management...
    if ((izn < 0)||(izn >= zoneVec.size()))
      return 0.0;
    double area = 0.0;
    for (int ist = 0; ist < nst; ++ist) {
      if (znost[ist] == izn) {
	const double n2[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
	area += MAG(n2);
      }
    }
    return 0.5*area;
  }

  void flagNodesOfFlaggedTris();
  void buildSubzoneIndexingFromZone();
  void dumpZones() const;
  int isManifold();
  void calcRmax();
  double getRmax();

  void ensureCentroid();
  void calcCentroid();
  void calcCenterOfMass(double xcom[3]);
  void calcFlaggedTrisCentroid(double (&_centroid)[3]);
  void getCentroid(double x[3]);

  void ensureBoundingBox();
  void calcBoundingBox();
  void getBoundingBox(double _bBox[6]);

  double getBoundingBoxRmax();
  void getBoundingBoxCenter(double _bBoxCenter[3]);

  void getAreaAndVolume(double& volume,double& area);

  void flipTrisRandom();
  void deleteIsolatedNodes();

  void deleteFlaggedTris();
  int countAndIndexDisjointSurfaceGroups(double& volume,double& area, double (&gcl)[3], vector<double>&groupVolVec,const bool detailed);
  void buildSubzoneData();
  void ensureSubzoneData();
  void clearSubzoneData();

  void buildZoneData();
  void ensureZoneData();
  void clearZoneData();

  void updateSubzoneData(const int isz);

  void buildStoszSzosz();
  void ensureStoszSzosz();
  void clearStoszSzosz();

  void reportZoneData(const bool subzone_info);
  void showSurfaceDiagnostics(const bool detailed,const bool subzone_info);
  void computeSurfaceCurvatureGaussBonnet(const int n_filter);
  void computeSurfaceCurvatureOld();
  void computeSubzoneCurvature(const double radius_cdf,const bool dump_histogram);
  void computeSubzoneCurvatureForSubzone(const int isz,double (*sp_normal)[3],double * curv_st);
  void buildSzoznFromLocalSzost();
  void buildSzoznFromGlobalSzost();
  void globalSzostToLocal();
  void localSzostToGlobal();
  void deleteSelectedSubzones(const vector<int>& subzoneVec);

  // routines for calculating geodesic distance from a point using FMM : O(NlogN)
  void fmmLoop(const int * stosp_i, const int * stosp_v, MinHeap& trialHeap);
  void calcGeoDistFromPoint(const double x[3]);
  void calcGeoDistFromFlaggedNodes();
  void getSubzoneBoundingBoxCenterAndDiagonal(double _bBox[6], double _bBoxCenter[3],double &_diag,IntFlag show_sz_flag);

  // these routines are used to provide stitch inputs (estimating the local grid scale)...
  void checkTriNbrUnitN(int *seed_st,const int ist,const int iseed,const int ist_seed,const double n_st_seed[3],
			const double n_mag_seed,const double dp_tol,int count);
  void getGapThicknesses(const double delta, const bool gap_viz);
  void getGapThicknessesCheap(const double planar_angle, const double delta, const bool gap_viz);

  void growNspData(const int nsp_new,const int nsp_old);
  void growNstData(const int nsp_new,const int nsp_old);
  void resizeNspData(const int nsp_new,const int nsp_old);
  void resizeNstData(const int nsp_new,const int nsp_old);

  // ============================================
  // these routines in SimpleSurface_periodic.cpp
  // ============================================

  void clearPeriodicity();
  void setTeostVec(vector<uint8>& pTeostVec,const vector<int>& open_edge_groups) const;
  bool setTeostVecs(vector<uint8>& p0TeostVec,vector<uint8>& p1TeostVec) const;
  bool setTeostVecs(vector<uint8>& p0TeostVec,vector<uint8>& p1TeostVec,const vector<int>& edges0,const vector<int>& edges1) const;
  int setPeriodicFlaggedZones(PeriodicTransform& pt,const bool b_force,const double crease_angle);
  int setPeriodicOpenEdgeGroups(const vector<int>& edges0,const vector<int>& edges1,PeriodicTransform& pt);
  int setPeriodic(vector<uint8>& p0TeostVec,vector<uint8>& p1TeostVec,const int ipt,const bool b_force,const double crease_angle);
  int resetPreviousPeriodic();
  bool checkLoopOrthogonalToX(const vector<uint8>& edgeVec) const;
  ///
  void computeOpenEdgeLoopNormal(double (&n0)[3],const int igroup) const;
  void computeEdgeLoopWeightCentroidNormal(double& wgt,double (&x0)[3],double (&n0)[3],const vector<uint>& openEdgeVec) const;
  ///
  void computeEdgeLoopWeightCentroidNormal(double& wgt,double (&x0)[3],double (&n0)[3],const vector<uint8>& openEdgeVec) const;
  int setPeriodicMergeVecs(vector<NodeNodeMerge>& node0Node1MergeVec,vector<EdgeNodeMerge>& edge0Node1MergeVec,vector<EdgeNodeMerge>& edge1Node0MergeVec,const vector<uint8>& p0TeostVec,const vector<uint8>& p1TeostVec,const int ipt);
  int setPeriodicMergeVecsForce(vector<NodeNodeMerge>& node0Node1MergeVec,vector<EdgeNodeMerge>& edge0Node1MergeVec,vector<EdgeNodeMerge>& edge1Node0MergeVec,const vector<uint8>& p0TeostVec,const vector<uint8>& p1TeostVec,const int ipt,const double crease_angle);
  void splitEdgesAndAddNodeNodeMerges(vector<NodeNodeMerge>& node0Node1MergeVec,vector<EdgeNodeMerge>& edge0Node1MergeVec,vector<EdgeNodeMerge>& edge1Node0MergeVec,vector<uint8>& p0TeostVec,vector<uint8>& p1TeostVec,const int ipt);
  int addPeriodicTransform(PeriodicTransform& pt);
  void mergeNodesAndSetPeriodicPbi(vector<NodeNodeMerge>& node0Node1MergeVec,const int ipt);

  int zipOpenEdgesChtFlaggedZones(const double crease_angle);

  // ============================================
  // these routines in SimpleSurface_transform.cpp...
  // ============================================

  void translateFlaggedTris(const double dx[3]);
  void translateFlaggedNodes(const double dx[3]);
  void mirrorFlaggedTris(const double xc[3],const double np[3]);
  void mirrorFlaggedNodes(const double xc[3],const double np[3]);
  void translateFlaggedTrisNormal(const double dn, const bool b_flagged_only);
  void rotateFlaggedTris(const double axis[3],const double point[3],const double angle_deg);
  void rotateFlaggedNodes(const double axis[3],const double point[3],const double angle_deg);
  void scaleFlaggedTris(const double sx[3],const double x0[3],const bool b_norm);
  void scaleFlaggedTrisRadially(const double factor,const double x0[3],const double axis[3],bool b_norm);
  void scaleFlaggedNodes(const double sx[3],const double x0[3],const bool b_norm);
  void scaleFlaggedNodesRadially(const double factor,const double x0[3],const double axis[3],const bool b_norm);

  void forceOrthogonalityOpenEdges(const vector<int>& open_edge_groups,const double dir[3]);
  void mirrorSelectedSubzones(const vector<int>& subzonesVec,const double xp[3],const double np[3]);
  void translateSelectedSubzones(const vector<int>& subzonesVec,const double dx[3]);
  void translateSelectedSubzonesNormal(const vector<int>& subzonesVec,const double dn);
  void rotateSelectedSubzones(const vector<int>& subzonesVec,const double axis[3],const double point[3],const double angle_deg);
  void scaleSelectedSubzones(const vector<int>& subzonesVec,const double sx[3],const double x0[3],const bool b_scale_centroid,const bool b_norm);
  void scaleSelectedSubzonesRadially(const vector<int>& subzonesVec,const double factor,const double x0[3],const double axis[3],const bool b_scale_centroid,const bool b_norm);

  void translateSurface(const double dx[3]);
  void rotateSurface(const double axis[3],const double point[3],const double angle_deg);
  void scaleSurface(const double sx[3],const double x0[3],const bool b_scale_centroid,const bool b_norm);

  void translateOpenEdges(const vector<int>& edge_indices, const double dx[3],const int n_duplicates);
  void rotateOpenEdges(const vector<int>& edge_indices, const double axis[3],const double point[3],const double angle_deg,const int n_duplicates);
  void scaleOpenEdges(const vector<int>& edge_indices, const double sx[3], const double x0[3], const bool b_scale_norm,const bool b_centroids,const int n_duplicates);
  void copySelectedSubzones(const vector<int>& subzonesVec,const int n_copies);
  void copyTranslateSelectedSubzones(const vector<int>& subzonesVec,const int n_copies,const double dx[3]);
  void copyRotateSelectedSubzones(const vector<int>& subzonesVec,const int n_copies,const double _axis[3],const double point[3],const double angle_deg);
  void roughenFlaggedSubzones(const double range[2],const int n_filter);

  // ============================================
  // these routines in SimpleSurface_geom.cpp...
  // ============================================

private:

  Adt2d<int> * stAdt2d; // the adt for integer (y,z) bounding boxes for the tris (ODD)
  double d2i_factor; // double-to-int factor used to convert y,z to integer (used with centroid)

public:
  void setFeatureAngle(const double degrees);
  double getFeatureAngle() const;

  void calcGcl(double gcl[3], const bool include_open_edge_groups = true);

  void ensureSposp();
  void buildSposp();
  void clearSposp();

  bool pointInTriangle(const double a[3], const double b[3], const double c[3], const double p[3], const double eps = 0.0);
  int howIsPointOnTriangle(const double a[3], const double b[3], const double c[3], const double p[3],const double dist_tol);

  // these routines are used to triangulate a 2d polygon using ear clipping: cost ~ O(n^2)
  inline double triangleSignedArea(const double xa, const double ya, const double xb, const double yb, const double xc, const double yc);
  bool between(const double xa, const double ya, const double xb, const double yb, const double xc, const double yc);
  bool collinear(const double xa, const double ya, const double xb, const double yb, const double xc, const double yc);
  bool diagonal(const int im1, const int ip1, const int n, const int* prev_node, const int* next_node, const double* x, const double* y);
  bool diagonalie(const int im1, const int ip1, const int n, const int* next_node, const double* x, const double* y);
  bool inCone(const int im1, const int ip1, const int n, const int* prev_node, const int* next_node, const double* x, const double* y);
  bool intersect(const double xa, const double ya, const double xb, const double yb, const double xc, const double yc, const double xd, const double yd);
  bool intersectProp(const double xa, const double ya, const double xb, const double yb, const double xc, const double yc, const double xd, const double yd);
  double polygonSignedArea(const int n, const double* x, const double* y);
  bool polygonTriangulate(int* triangles, const int n, const double* x, const double* y);

  bool pointIsInside(const double xp[3]);

  // ============================================
  // these routines in SimpleSurface_edges.cpp...
  // ============================================

  uint packEdge(const uint ist,const uchar i) const;
  uint packEdge(const int ist,const int i) const;
  void unpackEdge(uint& ist,uchar& i,const uint edge) const;
  void buildOpenEdgeGroups(const double crease_angle = 0.0);
  int getNOpenEdgeGroups() const;
  void ensureOpenEdgeGroups(const double crease_angle = 0.0); // default is to only break at multiloop pts
  void clearOpenEdgeGroups();
  int reportOpenEdgeGroupData(const bool detailed);

  void findNestedHalfEdgeLoops(vector<pair<int,int> >& gr_pairs,const double tol,const vector<int>& group_indices,const vector<HalfEdgeGroupData>& halfEdgeGroupDataVec) const;

  void ensureDynamicEdgeGroups(const DYNAMIC_EDGE_TYPE requested_type);
  void ensureDynamicEdgeGroups();
  void clearDynamicEdgeGroups();
  void buildDynamicEdgeGroups();
  // void buildMetalZoneEdgeGroups();
  // void buildFluidZoneEdgeGroups();
  void buildZoneEdgeGroups();
  void buildSubzoneEdgeGroups();
  void buildSelectedEdgeGroups();
  void buildFeatureEdgeGroups();
  void buildOpenEdgeGroups2();
  void buildNmeEdgeGroups();
  void selectEdgesInWindow(const double window[4][3],const bool b_strictly_inside,const bool b_open);
  void orderHalfEdgesInVec(vector<uint>& edgeVec);
  void orderEdgesInGroup(vector<uint>& edgeVec,const int igroup,const map<uint,int>& edge_to_group,const vector<HalfEdgeGroupData>& heGroupDataVec);
  void orderEdgesInGroup(vector<uint>& edgeVec,vector<int>& group_indices,const map<uint,int>& edge_to_group,const vector<HalfEdgeGroupData>& heGroupDataVec);

  bool getHalfEdgeGroup(int& igr,const uint ist,const uchar i);
  int buildHalfEdgeGroups(vector<HalfEdgeGroupData>& halfEdgeGroupDataVec, map<uint,int>& he_to_gr, bool (SimpleSurface::*halfEdgeCriterion)(const int,const int) const);

  void closeHalfEdgeLoops(vector<int>& group_indices, const bool ear_clipping,const bool open);
  void closeHalfEdgeLoop(vector<NewNode>& xsp_new,vector<NewTri>& spoed_new,const int igroup, const bool ear_clipping,const map<uint,int>& edge_to_group,const vector<HalfEdgeGroupData>& heGroupDataVec);

  void getNodesFromHalfEdges(vector<int>& nodes,const vector<uint>& edges,const bool loop) const;
  double getNodesAndLengthFromHalfEdges(vector<pair<int,double> >& nodes,const vector<uint>& edges,const bool loop) const;
  void facetGap(vector<NewTri>& newTris,const vector<int>& isp0Vec,const vector<int>& isp1Vec,const bool loop) const;
  void facetGap(vector<NewTri>& newTris,const vector<int>& isp0Vec,const vector<int>& isp1Vec,const double (*xsp)[3],const bool loop) const;
  void facetGapFromHalfEdges(vector<NewTri>& newTris,const vector<uint>& ise0Vec,const vector<uint>& ise1Vec,const bool loop) const;
  void sampleTriUniformRandom(vector<NewNode>& v_samples,const int n_samples,const int ist,const double delta) const;
  void sampleTriR2(vector<NewNode>& v_samples,const int n_samples,const int ist,const double delta,const int useEdge=0) const;
  void sampleTriR2AR(vector<NewNode>& v_samples,const int n_samples,const int ist,const double delta,const int useEdge=0) const;
  double computeTriAspectRatio(const int ist) const;
  int init2DVorPointsOnGroup(vector<NewNode>& internal_nodes,const double delta,const int np_fixed,const double factor,const int igr, const bool b_keep_edge_nodes,const int init_type) const;
  void reDiscretizeEdges(const vector<int>& edge_indices,const double delta,const int type,const bool b_keep_edge_nodes);

  // ============================================
  // these routines in SimpleSurface_imprint.cpp...
  // ============================================

  void imprintPlane(const double x_plane[3],const double n_plane[3],const double node_tol,const bool splitSubZones);
  void imprintCyl(const double x_cyl[3],const double n_cyl[3],const double r_cyl,const bool splitSubzones);
  // the common routine where tris with sp_dist==0 will be cut...
  void imprint(const double * const sp_dist,const bool splitSubzones);


  // ============================================
  // these routines in SimpleSurface_repair.cpp...
  // ============================================

  void deleteTrisWithCollocatedNodes();
  void deleteTrisWithIdenticalNodes();
  void mergeCollocatedNodes();
  void mergeCollocatedNodes(const double eps);
  string replaceSpacesWithUnderscores(const string& newname);

  void buildHalfEdgeVec(vector<uint>& halfEdgeVec,const map<const int,int> groupToLoop,const int iloop);

  void createPlanarOpenEdgeLoopTris(vector<NewTri>& newTrisVec,const vector<int>& group_indices);
  void createPlanarHalfEdgeLoopTris(vector<NewTri>& newTrisVec,const vector<int>& group_indices);
  void createPlanarEdgeLoopTris(vector<NewTri>& newTrisVec, vector<pair<int8,char> >& closeTeostVec);
  void closePlanarHalfEdgeLoops(vector<int>& group_indices,vector<HalfEdgeGroupData>& halfEdgeGroupDataVec, const bool nested, const double nested_tol=0.001,const string zonename="",const bool b_open = false);
  bool triangulateMarchingFront(vector<NewTri>& newTrisVec,set< std::pair<int,int> >& frontSet,const double (*x)[2],const int n);

  void closeHalfEdgeDonutLoops(vector<pair<int,int> >& donutLoopsVec,const bool open);
  void closeHalfEdgeDonutLoop(vector<NewTri>& newTris,const int group0, const int group1, const int izone, const int isz,const bool open);
  void splitMultiEdges();

private:
  void closeHalfEdgeLoopWith3Nodes(vector<NewTri>& spoed_new,vector<uint>& edgeVec);
  void closeHalfEdgeLoopWith4Nodes(vector<NewTri>& spoed_new,vector<uint>& edgeVec);
  void closeHalfEdgeLoopEarClipping(vector<NewTri>& spoed_new,vector<uint>& edgeVec);
  void closeHalfEdgeLoopMeanVisible(vector<NewNode>& xsp_new,vector<NewTri>& spoed_new,vector<uint>& edgeVec);

public:
  void flipFlaggedTris();
  void separateOverlappingTris(const double dn);

  void buildNonManifoldData();
  void ensureNonManifoldData();
  void clearNonManifoldData();
  bool hasNonManifoldData();

  void setNonManifoldCheckProperties(const double dist_tol,const double angle_tol_degrees,const bool check_self_intersections,const string output_type,const int n_samples,const int nbr_layers);
  void diagnoseManifold();
  void flagNonManifoldTris();
  void flagNonManifoldTris(const int nme_bits);
  double findMinNonzeroEdgeLength();
  void flagNonManifoldTrisAndNeighborhood(IntFlag& st_show,bool (SimpleSurface::*nmtCriterion)(const int) const,const int nbr_layers) const;
  bool isMultiNeighbor(const int ist) const;
  bool isNmeAdjacent(const int ist) const;
  bool isLinear(const int ist) const;
  bool isIntersecting(const int ist) const;
  bool isOverlapping(const int ist) const;
  bool isImprinting(const int ist) const;
  void searchForAndCategorizeIntersectingTris(const double dist_tol,const double angle_tol_degrees);
  void categorizeEdgeNeighbor(const int ist,const int c_ist,const bool aligned,const double angle_tol);
  void categorizeNodeNeighbor(const int ist,const int c_ist,const int (&shareNode)[3][2],const IntFlag& signToIst,const IntFlag& signToCIst,const bool isCoplanar,const double dist_tol);
  void categorizeMultiNodeButNotEdgeNeighbor(const int ist,const int c_ist,const double angle_tol);
  void categorizeCoplanarVicinityTri(const int ist,const int c_ist,const double dist_tol);
  void categorizeVicinityTri(const int ist,const int c_ist,IntFlag& signToIst,IntFlag& signToCIst,const double dist_tol);
  int checkForValidTriPiercing(const int ist,const int c_ist) const;
  int categorizeOneNodeInPlaneVicinityTri(const int ist,const int c_ist,IntFlag& signToIst,const double dist_tol);
  bool categorizeTwoNodesInPlaneVicinityTri(const int ist,const int c_ist,IntFlag& signToIst,const double dist_tol);
  void openLinearTris(IntFlag& st_flag);

  // structures used for open edge loop closing
  struct oegNode {
    int isp;
    int match;
    double minD2;  // can be reused as any lengthscale
    double tol;
    oegNode() {
      isp = -1;
      match = -1;
      minD2 = -1.0;
      tol = -1.0;
    };
    oegNode(const int _isp) {
      isp = _isp;
      match = -1;
      minD2 = -1.0;
      tol = -1.0;
    };
    oegNode(const int _isp, const double d2) {
      isp = _isp;
      match = -1;
      minD2 = d2;
      tol = -1.0;
    };
    int setMinD2(const double d2) {
      if (minD2 < 0) {
        minD2 = d2;
        return 1;  // minD2 was set
      }

      if (d2 > 0.0 && d2 < minD2) {
        minD2 = d2;
        return 1;  // minD2 was set
      }
      else {
        return 0;
      }
    };
    void clearMinD2() {
      minD2 = -1.0;
    };
    void setMatch(const int isp_match, const double d2) {
      int wasSet = setMinD2(d2);
      if (wasSet) match = isp_match;
    };
    void setTol(const double _tol) {
      tol = _tol;
    }
  };

  struct oegEdge {
    int ied;
    double d2;
    int isp0;
    int isp1;
    oegEdge(const int _ied, const int _isp0, const int _isp1, const double _d2) {
      ied = _ied;
      isp0 = _isp0;
      isp1 = _isp1;
      d2 = _d2;
    };
  };

  void mergeNodesOfOpenEdges(const double merge_tol);
  void mergeNodesOfOpenEdgeGroups(const int loop0, const int loop1, const double merge_tol);
  void mergeNodesOfOpenEdgeGroupsOnSubzones(const int zone0, const int zone1, const double merge_tol);
  void mergeNodesOfOpenEdgeVectors(vector<oegEdge>& edges0, vector<oegEdge>& edges1, const double merge_tol);

  void flipSelectedSubzones(const vector<int>& subzonesVec);
  void flipSelectedGroups(const vector<int>& groupsVec);
  void flipSurface();
  void alignNormalsAuto();
  void alignNormals(const int st_seed0);
  void deleteLinearTris();
  void splitLinearTriNeighbor();
  void deleteMultiNeighborTris();

  // surface mesh refinement...

  void refineFlaggedTris(const double delta);

private:

  // supports refineFlaggedTris...
  void recursivelyRefineTri(const int ist,const double delta,map<const pair<int,int>,int>& newNodeMap,int& nsp_max,int& nst_max);

public:

  void checkSelfIntersectionsNew();

  // ============================================
  // these routines in SimpleSurface_zoning.cpp...
  // ============================================

  void migrateSubzonesToZones(const vector<int> zone_indices);
  void setSubzonesToZones();
  void combineZonesWithSameName();
  int addNewZone(const string& zonename); // returns new zone index
  void splitFlaggedSubzonesIntoTris();
  void splitFlaggedSubzonesInSphere(const double x[3],const double r);
  void splitFlaggedSubzonesIntoDisjointSubzones();
  void splitFlaggedSubzonesAtCreaseAngle(const double crease_angle);
  void splitFlaggedSubzonesByMarch(const double* split_data, const string& dest_zone);
  void splitFlaggedSubzonesByCoordinate(const int isz_mode,const double* split_data,const string& dest_zone);
  void splitFlaggedSubzonesByPlane(const double x[3],const double n[3]);
  int getZoneIndex(const string& zonename) const;
  bool pruneEmptyZonesAndSubzones();
  void renameZone(const string& name,const string& newname);
  void renameZone(const int index, const string& newname);
  void moveFlaggedSubzonesToZone(const int izn_dest);
  void moveFlaggedTrisToZone(const int izn_dest);
  void flagNodesFromSubzoneVec(const vector<int>& subzonesVec,const bool b_exclude_boundaries = false);
  void flagTrisFromSubzoneVec(const vector<int>& subzonesVec);
  void flagNodesAndTrisFromSubzoneVec(const vector<int>& subzonesVec);
  void indexNodesAndTrisFromSubzoneVec(int& sp_count, int& st_count,const vector<int>& subzonesVec);
  void flagZonesFromZoneNameVec(const vector<string>& zoneNamesVec);
  void flagZonesFromZoneIndexVec(const vector<int>& zone_indices);
  void selectSimilarSubzone(const int sz_index,const bool bSameZone,const bool byArea, const bool byNormal, const double diffTol, const double normalTol);
  void selectSimilarOpenEdgeGroup(const int loop_index,const bool byLength,const bool byAngleSum, const bool byCentroidMean, const double diffTol, const double angleSumTol, const double centroidMeanTol);

  // also in SimpleSurface_zoning.cpp because it is like splitFlaggedSubzonesIntoDisjointSubzones...
  void makeInjectorForFlaggedSubzones();
  void characterizePipe(double &l_on_d,double theta[2],const int isz);
  void calcClosestCylinderForFlaggedTris(double xc[3],double nc[3],double& r);
  void calcClosestCylinderForFlaggedTris(double nc[3],double xc[3]);

  // put these into a "SimpleSurface_part.hpp eventually...
  void makePartBlForFlaggedSubzones(const string& name,const double dn,const double dt,const int n);
  void writeTecplotFF(const string& name,const double (* const xsp_ff)[3],const int nsp_ff,const int (* const spost_ff)[3],const int nst_ff);
  void writeSbinFF(const string& name,const double (* const xsp_ff)[3],const int nsp_ff,const int (* const spost_ff)[3],const int nst_ff);

  // try and get some info on the flagged zones...
  void queryFlaggedSubzones();

  void selectZonesInWindow(const double window[4][3],const bool b_strictly_inside);

  double identifyFeatureXp(vector<NewNode>& fixedXp,vector<pair<int,int> >& featureEdgeVec,vector<pair<int,int> >& adj_tris,vector<pair<int,int> >& isp_to_nearest_fixed,const vector<int>& edge_indices,const double delta,const bool no_adj_feature,const bool b_keep_edge_nodes);
  double identifyFeatureXp(vector<NewNode>& fixedXp,vector<pair<int,int> >& featureEdgeVec,vector<pair<int,int> >& adj_tris,vector<pair<int,int> >& isp_to_nearest_fixed,const int (*stost)[3],const double delta,const bool no_adj_feature,const bool b_keep_edge_nodes);
  double rediscretizeEdges(vector<NewNode>& fixedXp,vector<pair<int,int> >& adj_tris,vector<pair<int,int> >& isp_to_nearest_fixed,const multimap<int,int>& ispToEdgeMap,const vector<pair<int,int> >& featureEdgeVec,const double delta,const bool no_adj_feature,const bool edge_only,const bool b_keep_edge_nodes);
  void buildStostForFlaggedZones(int (**stost)[3]);
  int groupFlaggedTrisByFeatures(const int (*stost)[3]);
  void makeSelectedTrisDelaunay(const int iter_limit);
  int triangulateVorPoints(vector<NewTri>& triVec,const double (*xp)[3],const int * stoxp,const int * nboxp_i,const vector<int>& nboxp_v,const int np) const;
  void updateSurfaceWithRetessellation(const vector<NewNode>& new_nodes,const int np_fixed,const vector<NewTri>& new_adj_tris,const vector<NewTri>& new_tris);
  void zipAdjacentTris(vector<NewTri>& new_tris,vector<pair<int,int> >& adj_tris,const multimap <pair<int,int>,int>& edge_to_fixedNodes,const map<int,pair<int,bool> >& isp_to_fixedNodes,vector<NewNode>& fixedNodes);
  void reDiscretizeZones(const vector<int>& zone_indices,const double delta,const double edge_factor,const double seed_factor,const int n_lloyd,const int type,const bool b_delaunay,const bool b_keep_edge_nodes,const int init_type,const bool b_power_diagram,const double growth_factor,const double growth_power);

  void splitZonesAfterCopyAll(const int n);

  // ============================================
  // these routines in SimpleSurface_material.cpp...
  // ============================================

  void clearMaterials();
  void ensureMaterialZones();
  void reportMaterialInfo();

  // ============================================
  // these routines in SimpleSurface_morph.cpp...
  // ============================================

  void morphTranslateFlaggedZones(const double dx[3]);
  bool checkMorphedZone(const int zone);
  void morphTranslate(const set<int>& moved, const set<int>& freed, const double dx[3], const int nrelax, const bool use_bump);
  void morphRotateCylinder(const int cyl_id, const int top_id, const int base_id, const double phi,const int nrelax);
  void morphRotateHalfCylinder(const int cyl_id, const int top_id, const int base_id, const int side_id, const double phi,const int nrelax);
  void morphCorner(const int top_id, const int side_id, const double dx, const double dl,
		   const bool use_bump, const double dw);
  void morphStretch(const int top, const double dn, const bool use_bump, const double dl);
  void morphCylinder(const int cyl_id, const int top_id, const int bot_id, const double L,
		     const double theta0, const double a, const double b, const int n);
  void morphCylinderZ(const int cyl_id, const int top_id, const int bot_id, const double x0, const double y0,
		      const double r0, const double L, const double theta0, const double a, const double b, const int n);
  void morphCylinderR(const int cyl_id,const double factor);
  void morphGaussianBump();
  void smoothMorph(double (*dx)[3]);
  void nsmoothFlaggedSubzones(const int niter,const double relax,const bool b_poisson);
  void nsmoothNearPoint(const int niter,const double x[3],const double r);
  void morphRotateSelectedSubzones(const vector<int>& subzonesVec,const double axis[3], const double point[3],const double angle_deg);

  // ============================================
  // these routines in SimpleSurface_seed.cpp...
  // ============================================
  void testSeed(Param * param,const bool help);
  void testNonmanifoldSeeds(Adt<double> * stAdt,const int n_samples);
  bool isTriValidBasedOnSeeds(vector<int>& neighborhood_tris,const int ist, Adt<double> * stAdt, const int n_samples);

  // ============================================
  // these routines in SimpleSurface_vor2d.cpp...
  // ============================================
  void retri();
  void lloydIterateInternalVorPoints(double (*xp)[3],int * stoxp, double * deltap,const int np_fixed,const int np_local,int * nboxp_i,vector<int>& nboxp_v,const int (*stost)[3],const int maxiter,const double delta,const double zero,const int igr,const bool b_power_diagram,const double growth_factor,const double growth_power) const;

};

#endif
