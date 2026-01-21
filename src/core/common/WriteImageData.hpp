#ifndef _WRITE_IMAGE_DATA_HPP_
#define _WRITE_IMAGE_DATA_HPP_

#include "MiscUtils.hpp"
#include "PlaneData.hpp"

enum RGB_MODE {
  DATA,
  NORMALS
};

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

class WriteImageData {
public:

  string paramStr;

  int interval;
  int step; // used to create solver written image filenames
  
  string app;
  bool b_app;

  RGB_MODE rgb_mode;

  string name; // file prefix
  bool b_name;

  float geom_plane_frac; // fractional distance through geometry in direction "np"
  bool b_geom_plane_frac;

  float xp[3],np[3]; // plane normal and point
  bool b_geom_plane;

  // iso stuff...

  vector<float> iso_values;
  vector<string> iso_vars;
  bool b_geom_iso;
  bool b_var_on_iso;
  string var_on_iso;
  bool b_range_on_iso;
  float range_on_iso[3];
  bool b_clip_iso;

  // iso primitive stuff...
  bool b_geom_primitive;
  string primitive_name;
  double primitive_data[8];

  // vector viz stuff..
  float vector_scale;
  float vector_density;

  bool b_up;
  float up[3]; // up vector: has to have non-zero cross with np

  bool b_width;
  float width; // width of image in simulation units

  bool b_range;
  float range[3];

  bool b_size;
  int size[2]; // pixel size: width x height

  bool b_target;
  float target[3];

  bool b_camera;
  float camera[3];

  bool b_hidezones;
  bool b_hidezonenames;
  bool b_hidesubzones;
  std::set<int> hiddenZonesSet;
  std::set<string> hiddenZoneNamesSet;
  std::set<int> hiddenSubzonesSet;

  bool b_blank;
  std::vector< PlaneData<float> > vBlankPlaneData;
  bool b_blank_data_plane;
  bool b_no_particle_blanking;

  bool b_var;
  string var;

  bool b_edgeIndex;
  int edgeIndex;

  bool showBfSurface;
  bool b_var_on_surface;
  string var_on_surface;

  bool b_range_on_surface;
  float range_on_surface[3];

  bool b_particle_magnification;
  double particle_magnification;

  bool b_particle_size;
  double particle_size;

  bool b_var_on_particle;
  string var_on_particle;

  bool b_range_on_particle;
  float range_on_particle[2];

  bool b_snapshot;
  string snapshot_name;
  int snapshot_f;
  int snapshot_inc;
  int snapshot_l;

  bool showSurfaceMesh;
  bool showInteriorMesh;
  bool showOpenEdges;
  bool showDynamicEdges;
  DYNAMIC_EDGE_TYPE dEdgeType;

  bool  b_volvis;
  string volvis_var;
  string volvis_type; // linear: dI(phi,ds) ~ max(min(phi,phi_max),phi_min)*ds,
                      // gaussian: dI(phi,ds) ~ exp(-0.5*((phi-phi_avg)/phi_std)^2)*dx
  float volvis_aux_data[2]; // min/max (linear), avg/std (gaussian)
  bool b_volvis_aux_data;
  bool b_volvis_surface;

  bool b_solver_surface; // like PLIC
  string solver_surface;

  bool b_cell_flooded;

  bool b_transform_cyl_x;
  bool b_transform_cyl_z;
  double r0,cos_theta0,sin_theta0;
  double cyl_center[3];

  bool b_proj_cyl;
  double proj_cyl_xp[3];
  double proj_cyl_np[3];

  bool b_hide_surfaces_without_data; // when true, surfaces without the requested data are hidden

  void init() {

    interval  = -1; // CI changed this to indicate one-offs vs persistent images
    step = 0;
    b_app = false;
    rgb_mode = DATA;
    name = "";
    b_name  = false;
    b_geom_plane_frac  = false;
    b_geom_plane = false;
    b_geom_iso = false;
    b_geom_primitive = false;
    b_up    = false;
    b_width = false;
    b_range = false;
    b_range_on_surface = false;
    b_var   = false;
    b_size = false;
    b_hidezones = false;
    b_hidezonenames = false;
    b_hidesubzones = false;
    size[0] = 640;
    size[1] = 480;
    b_target = false;
    b_camera = false;
    b_blank = false;
    b_blank_data_plane = true; // default to blanking in front
    b_no_particle_blanking = false; // default to blank particles
    b_edgeIndex = false;
    showBfSurface = true; // default to using bf surface instead of sbin
    b_var_on_surface = false;
    b_particle_magnification = false;
    b_particle_size = false;
    b_range_on_particle = false;
    b_var_on_particle = false;
    b_range_on_iso = false;
    b_clip_iso = false;
    b_var_on_iso = false;
    b_snapshot = false;

    FOR_I3 {
      range_on_iso[i] = 0.f;
      range[i] = 0.f;
      range_on_surface[i] = 0.f;
    }

    showSurfaceMesh = false;
    showInteriorMesh = false;
    showOpenEdges = false;
    showDynamicEdges = false;
    dEdgeType = MULTI;

    b_volvis = false;
    b_volvis_aux_data = false;
    b_volvis_surface = false;
    volvis_type = "X-RAY";

    b_solver_surface = false;
    solver_surface = "";

    vector_scale = 1.0;
    vector_density = 1.0;

    b_cell_flooded = false;
    b_transform_cyl_x = false;
    b_transform_cyl_z = false;

    // when using the cyl transform for images, the user can specify a
    // non-zero axis by adding the tokens:
    // CYL_CENTER <x0> <y0> <z0>
    cyl_center[0] = 0.0;
    cyl_center[1] = 0.0;
    cyl_center[2] = 0.0;

    b_proj_cyl = false;
    proj_cyl_xp[0] = 0.0;
    proj_cyl_xp[1] = 0.0;
    proj_cyl_xp[2] = 0.0;
    proj_cyl_np[0] = 0.0;
    proj_cyl_np[1] = 0.0;
    proj_cyl_np[2] = 0.0;


    b_hide_surfaces_without_data = false;

  }

  WriteImageData() {

    init();

  }

  WriteImageData(Param * param) {

    init(param);

  }

  int init(Param * param,const int step_ = 0) { // underscore because step is a class variable - why? DP

    // the wid init takes the step_ with a default value of
    // zero. It is simply used for fast return when step_ is
    // known and the step_%interval can be checked quickly
    // and the whole initialization skipped.

    init();

    try {

      int iarg = 0;
      while (iarg < param->size()) {
        string token = param->getString(iarg++);
        if (token == "NAME") {
          name = param->getString(iarg++);
          b_name = true;
        }
        else if (token == "INTERVAL") {
          interval = param->getInt(iarg++);
          if (interval > 0) {
            if (step_%interval != 0) return -1; // we are not going to write, so just skip
          }
          else {
            CWARN(" > WRITE_IMAGE INTERVAL expects a positive integer");
            // invalid value entered, so treat as unset (default value)
            interval = -1;
          }
        }
	else if (token == "APP") {
	  app = param->getString(iarg++);
	  b_app = true;
	}
        else if (token == "RGB_MODE"){
          token = param->getString(iarg++);
          if (token == "DATA") {
            rgb_mode = DATA;
          }
          else if (token == "NORMALS") {
            rgb_mode = NORMALS;
          }
          else {
            if (mpi_rank == 0) cout << "Warning: unrecognized WRITE_IMAGE RGB_MODE token \"" << token << "\", using default" << endl;
          }
        }
        else if (token == "GEOM") {
          token = param->getString(iarg++);
          if (token == "Z_PLANE_FRAC") {
            geom_plane_frac = param->getDouble(iarg++);
            np[0] = 0.0;
            np[1] = 0.0;
            np[2] = 1.0;
            b_geom_plane = true;
            b_geom_plane_frac = true;
          }
          else if (token == "Y_PLANE_FRAC") {
            geom_plane_frac = param->getDouble(iarg++);
            np[0] = 0.0;
            np[1] = -1.0;
            np[2] = 0.0;
            b_geom_plane = true;
            b_geom_plane_frac = true;
          }
          else if (token == "X_PLANE_FRAC") {
            geom_plane_frac = param->getDouble(iarg++);
            np[0] = 1.0;
            np[1] = 0.0;
            np[2] = 0.0;
            b_geom_plane = true;
            b_geom_plane_frac = true;
          }
          else if (token == "PLANE_FRAC") {
            np[0] = param->getDouble(iarg++);
            np[1] = param->getDouble(iarg++);
            np[2] = param->getDouble(iarg++);
            geom_plane_frac = param->getDouble(iarg++);
            b_geom_plane = true;
            b_geom_plane_frac = true;
          }
          else if (token == "PLANE") {
            xp[0] = param->getDouble(iarg++);
            xp[1] = param->getDouble(iarg++);
            xp[2] = param->getDouble(iarg++);
            np[0] = param->getDouble(iarg++);
            np[1] = param->getDouble(iarg++);
            np[2] = param->getDouble(iarg++);
            // make np a unit normal...
            const double mag = MAG(np);
            if (mag <= 0.0) {
              if (mpi_rank == 0) cout << "Warning: invalid PLANE normal direction specified; skipping PLANE visualization" << endl;
            }
            else {
              FOR_I3 np[i] /= mag;
              b_geom_plane = true;
            }
          }
          else if (token == "ISO") {
            b_geom_iso = true;
            while (token == "ISO") {
              iso_vars.push_back(param->getString(iarg++));
              iso_values.push_back(param->getDouble(iarg++));
              if (iarg >= param->size()){
                ++iarg;
                break;
              }
              token = param->getString(iarg++);
            }
            --iarg;
          }
          // the following simple primitives are handled like iso's...
          else if ( (token == "CYL_X") || (token == "CYL_Y") || (token == "CYL_Z") ) {
            b_geom_primitive = true;
            primitive_name = token;
            // x0[3] and x1[3] determined from bbox and primitive name
            primitive_data[0] = param->getDouble(iarg++); // r
          }
          else if ( token == "CYLINDER" ) {
            b_geom_primitive = true;
            primitive_name = token;
            primitive_data[0] = param->getDouble(iarg++); // x0
            primitive_data[1] = param->getDouble(iarg++); // y0
            primitive_data[2] = param->getDouble(iarg++); // z0
            primitive_data[3] = param->getDouble(iarg++); // x1
            primitive_data[4] = param->getDouble(iarg++); // y1
            primitive_data[5] = param->getDouble(iarg++); // z1
            primitive_data[6] = param->getDouble(iarg++); // r
          }
          else if ( (token == "TCONE_X") || (token == "TCONE_Y") || (token == "TCONE_Z") ) {
            b_geom_primitive = true;
            primitive_name = token;
            // x0[3] and x1[3] determined from bbox and primitive name
            primitive_data[0] = param->getDouble(iarg++); // r0
            primitive_data[1] = param->getDouble(iarg++); // r1
          }
          else if ( token == "TCONE" ) {
            b_geom_primitive = true;
            primitive_name = token;
            primitive_data[0] = param->getDouble(iarg++); // x0
            primitive_data[1] = param->getDouble(iarg++); // y0
            primitive_data[2] = param->getDouble(iarg++); // z0
            primitive_data[3] = param->getDouble(iarg++); // r0
            primitive_data[4] = param->getDouble(iarg++); // x1
            primitive_data[5] = param->getDouble(iarg++); // y1
            primitive_data[6] = param->getDouble(iarg++); // z1
            primitive_data[7] = param->getDouble(iarg++); // r1
          }
          else if ( token == "SPHERE" ) {
            b_geom_primitive = true;
            primitive_name = token;
            primitive_data[0] = param->getDouble(iarg++); // xc
            primitive_data[1] = param->getDouble(iarg++); // yc
            primitive_data[2] = param->getDouble(iarg++); // zc
            primitive_data[3] = param->getDouble(iarg++); // r
          }
          else {
            if (mpi_rank == 0) cout << "Warning: skipping unrecognized WRITE_IMAGE GEOM token \"" << token << "\"" << endl;
          }
        }
        else if (token == "CAMERA") {
          if (b_camera) {
            if (mpi_rank == 0) cout << "Warning: CAMERA already specified; overwriting values" << endl;
          }
          camera[0] = param->getFloat(iarg++);
          camera[1] = param->getFloat(iarg++);
          camera[2] = param->getFloat(iarg++);
          b_camera = true;
        }
        else if (token == "TARGET") {
          if (b_target) {
            if (mpi_rank == 0) cout << "Warning: TARGET already specified; overwriting values" << endl;
          }
          target[0] = param->getFloat(iarg++);
          target[1] = param->getFloat(iarg++);
          target[2] = param->getFloat(iarg++);
          b_target = true;
        }
        else if (token == "VAR") {
          var = param->getString(iarg++);
          b_var = true;
          if (var == "MESH" || var == "mesh") showInteriorMesh = true;
        }
        else if (token == "MESH_ON_SURFACE") {
          showSurfaceMesh = true;
        }
        else if (token == "RANGE_ON_PARTICLE") {
          range_on_particle[0] = param->getFloat(iarg++);
          range_on_particle[1] = param->getFloat(iarg++);
          b_range_on_particle = true;
        }
        else if (token == "VAR_ON_PARTICLE") {
          b_var_on_particle = true;
          var_on_particle = param->getString(iarg++);
        }
        else if (token == "PARTICLE_MAGNIFICATION") {
          particle_magnification = param->getDouble(iarg++);
          b_particle_magnification = true;
        }
        else if (token == "PARTICLE_SIZE") {
          particle_size = param->getDouble(iarg++);
          b_particle_size = true;
        }
        else if (token == "RANGE_ON_ISO") {
          range_on_iso[0] = param->getFloat(iarg++);
          range_on_iso[1] = param->getFloat(iarg++);
          if (iarg<param->size() && param->getString(iarg)=="CLIP"){
            b_clip_iso = true;
            ++iarg;
          }
          b_range_on_iso = true;
        }
        else if (token == "NBIN_ISO") {
          range_on_iso[2] = (float)param->getInt(iarg++);
          if ((range_on_iso[2]<2) || (range_on_iso[2]>256)) {
            if (mpi_rank == 0) cout << "Warning: only NBIN_SIO values 2 <= NBIN_SIO <= 256 are processed; skipping contour levels" << endl;
            range_on_iso[2] = 0.f;  // only process if valid
          }
        }
        else if (token == "CLIP_ISO") {
          b_clip_iso = true;
        }
        else if (token == "VECTOR_SCALE") {
          vector_scale = param->getFloat(iarg++);
        }
        else if (token == "VECTOR_DENSITY") {
          vector_density = param->getFloat(iarg++);
        }
        else if (token == "VAR_ON_ISO") {
          b_var_on_iso = true;
          // VAR_ON_ISO is used instead of VAR when available
          // this is due to ISO being a type of GEOM which uses VAR
          // will have to be careful with logic and intended behavior
          var_on_iso = param->getString(iarg++);
        }
        else if (token == "RANGE_ON_SURFACE") {
          range_on_surface[0] = param->getFloat(iarg++);
          range_on_surface[1] = param->getFloat(iarg++);
          b_range_on_surface = true;
        }
        else if (token == "NBIN_SURFACE") {
          range_on_surface[2] = (float)param->getInt(iarg++);
          if ((range_on_surface[2]<2) || (range_on_surface[2]>256)) {
            if (mpi_rank == 0) cout << "Warning: only NBIN_SURFACE values 2 <= NBIN_SURFACE <= 256 are processed; skipping contour levels" << endl;
            range_on_surface[2] = 0.f;  // only process if valid
          }
        }
        else if (token == "VAR_ON_SURFACE") {
          b_var_on_surface = true;
          var_on_surface = param->getString(iarg++);
          /*
            if (val=="X"||val=="x")
            var_on_surface = 0;
            else if (val=="Y"||val=="y")
            var_on_surface = 1;
            else if (val=="Z"||val=="z")
            var_on_surface = 2;
            else{
            cout << "Warning: unrecognized VAR_ON_SURFACE " << val << ", valid options are X, Y, or Z" << endl;
            b_var_on_surface = false;
            --iarg;
            }
          */
        }
        else if (token == "NO_GEOM_BLANKING") {
          // prevent data plane from being a blank plane...
          b_blank_data_plane = false;
        }
        else if (token == "NO_PARTICLE_BLANKING") {
          // prevent data plane from being a blank plane...
          b_no_particle_blanking = true;
        }
        else if (token == "CELL_FLOODED") {
          b_cell_flooded = true;
        }
        else if (token == "HIDE_ZONES_WITHOUT_DATA") {
          b_hide_surfaces_without_data = true;
        }
        else if (token == "NO_BF_SURFACE") {
          // use the triangulated surface instead of the boundary face in charles viz
          showBfSurface = false;
        }
        else if (token == "SOLVER_SURFACE") {
          b_solver_surface = true;
          solver_surface = param->getString(iarg++);
        }
        else if (token == "SHOW_OPEN_EDGES") {
          showOpenEdges = true;
        }
        else if (token == "HIGHLIGHT_DYN_EDGES") {
        showDynamicEdges = true;
        token = param->getString(iarg++);
        if (token == "MULTI") dEdgeType=MULTI;
        else if (token == "METAL_BOUNDARY") dEdgeType=METAL_BOUNDARY;
        else if (token == "FLUID_BOUNDARY") dEdgeType=FLUID_BOUNDARY;
        else if (token == "ZONE_BOUNDARY") dEdgeType=ZONE_BOUNDARY;
        else if (token == "SUBZONE_BOUNDARY") dEdgeType=SUBZONE_BOUNDARY;
        else if (token == "FEATURE") dEdgeType=FEATURE;
        else if (token == "OPEN") dEdgeType=OPEN;
        else if (token == "NME") dEdgeType=NME;
        else if (token == "SELECTED_BOUNDARY") dEdgeType=SELECTED_BOUNDARY;
        else {
          cout << "\nWARNING: unrecognized edge type: " << token << "; skipping\n"<< endl;
          dEdgeType=NONE;
          showDynamicEdges = false;
        }
      }
        else if (token == "SIZE") {
          int img_width = param->getInt(iarg++);
          int img_height = param->getInt(iarg++);

          if ( (img_width < 1) || (img_width >= 16000) ) {
            if (mpi_rank == 0) cout << "Warning: unsupported image width: " << img_width << "px; using default of 640px" << endl;
            img_width = 640;
          }

          if ( (img_height < 1) || (img_height >= 16000) ) {
            if (mpi_rank == 0) cout << "Warning: unsupported image height: " << img_height << "px; using default of 480px" << endl;
            img_height = 480;
          }

          size[0] = img_width;
          size[1] = img_height;
          b_size = true;
        }
        else if (token == "UP") {
          up[0] = param->getFloat(iarg++);
          up[1] = param->getFloat(iarg++);
          up[2] = param->getFloat(iarg++);
          b_up = true;
        }
        else if (token == "BLANK") {
          if (!b_blank) b_blank = true;
          token = param->getString(iarg++);
          if (token == "PLANE") {
            const float center[3] = {float(param->getDouble(iarg)),
                                     float(param->getDouble(iarg+1)),
                                     float(param->getDouble(iarg+2))};
            iarg += 3;
            const float normal[3] = {float(param->getDouble(iarg)),
                                     float(param->getDouble(iarg+1)),
                                     float(param->getDouble(iarg+2))};
            iarg += 3;
            vBlankPlaneData.push_back( PlaneData<float>(center,normal) );
          }
          else {
            if (mpi_rank == 0) cout << "Warning: skipping unrecognized BLANK token \"" << token << "\"" << endl;
          }
        }
        else if (token == "WIDTH") {
          width = param->getFloat(iarg++);
          b_width = true;
        }
        else if (token == "RANGE") {
          range[0] = param->getFloat(iarg++);
          range[1] = param->getFloat(iarg++);
          b_range = true;
        }
        else if (token == "NBIN") {
          range[2] = (float)param->getInt(iarg++);
          if ((range[2]<2) || (range[2]>256)) {
            if (mpi_rank == 0) cout << "Warning: only NBIN values 2 <= NBIN <= 256 are processed; skipping contour levels" << endl;
            range[2] = 0.f;  // only process if valid
          }
        }
        else if (token == "EDGE") {
          edgeIndex = param->getInt(iarg++);
          if ((edgeIndex < 32768) || (edgeIndex >= 65535)) {
            if (mpi_rank == 0) cout << "Warning: invalid edge index " << edgeIndex << " specified; skipping" << endl;
          }
          else b_edgeIndex = true;
        }
        else if (token == "HIDE_SUBZONES") {
          token = param->getString(iarg++);
          if (token == "ALL") {
            b_hidesubzones = true;
            hiddenSubzonesSet.clear();
            hiddenSubzonesSet.insert(-1);
          }
          else {
            vector<string> hiddenSubzonesVec;
            MiscUtils::splitCsv(hiddenSubzonesVec,token);
            if (!hiddenSubzonesVec.empty()) {
              b_hidesubzones = true;
              hiddenSubzonesSet.clear();
              for (int i=0, limit=hiddenSubzonesVec.size(); i<limit; ++i) {
                hiddenSubzonesSet.insert(atoi(hiddenSubzonesVec[i].c_str()));
              }
            }
          }
        }
        else if (token == "HIDE_ZONES") {
          // DP/ME: can this be cleaned out?
          token = param->getString(iarg++);
          if (token == "ALL") {
            b_hidezones = true;
            hiddenZonesSet.clear();
            hiddenZonesSet.insert(-1); // TODO: is this supported?
          }
          else {
            vector<string> hiddenZonesVec;
            MiscUtils::splitCsv(hiddenZonesVec,token);
            if (hiddenZonesVec.size()) {
              b_hidezones = true;
              hiddenZonesSet.clear();
              for (int i=0, limit=hiddenZonesVec.size(); i<limit; ++i) {
                hiddenZonesSet.insert(atoi(hiddenZonesVec[i].c_str()));
              }
            }
          }
        }
        else if (token == "HIDE_ZONES_NAMED") {
          token = param->getString(iarg++);
          vector<string> hiddenZonesVec;
          MiscUtils::splitCsv(hiddenZonesVec,token);
          if (hiddenZonesVec.size()) {
            b_hidezonenames = true;
            hiddenZoneNamesSet.clear();
            hiddenZoneNamesSet.insert(hiddenZonesVec.begin(),hiddenZonesVec.end());
            // CLEAN
            // for (int i=0, limit=hiddenZonesVec.size(); i<limit; ++i) {
            //   hiddenZoneNamesSet.insert(hiddenZonesVec[i]);
            // }
          }
        }
        else if (token == "VOLVIS_VAR") {
          b_volvis = true;
          volvis_var = param->getString(iarg++);
        }
        else if (token == "VOLVIS_SURFACE") {
          b_volvis_surface = true;
        }
        else if (token == "VOLVIS_AUX_DATA") {
          b_volvis_aux_data = true;
          volvis_aux_data[0] = param->getFloat(iarg++);
          volvis_aux_data[1] = param->getFloat(iarg++);
        }
        else if (token == "VOLVIS_TYPE") {
          token = param->getString(iarg++);
          if (token == "X-RAY" || token == "LINEAR" || token == "GAUSSIAN" || token == "SURFACE") {
            volvis_type = token;
          }
          else {
            if (mpi_rank == 0) cout << "Warning: unrecognized VOLVIS_TYPE token \"" << token << "\", using default" << endl;
          }
          // x-ray doesn't require var, so ensure b_volvis set here
          if (volvis_type == "X-RAY" || volvis_type == "SURFACE") {
            b_volvis = true;
            volvis_var = "none";
          }
        }
        else if (token == "SNAPSHOT") {
          b_snapshot = true;
          snapshot_name = param->getString(iarg++);
          snapshot_f    = param->getInt(iarg++);
          snapshot_inc  = param->getInt(iarg++);
          snapshot_l    = param->getInt(iarg++);
        }
        else if (token == "TRANSFORM_CYL_X") {
          // usage:
          // TRANSFORM_CYL_X <r0> <theta0-radians>
          // projects image and solution data into CYL_X space (x,y,z)->(x,r,r0*theta)
          b_transform_cyl_x = true;
          r0 = param->getDouble(iarg++);
          const double theta0 = param->getDouble(iarg++);
          cos_theta0 = cos(theta0);
          sin_theta0 = sin(theta0);
          //if (mpi_rank == 0) cout << " > XXXXXXXXXXXXXXXXXXXXXXXX TRANSFORM_CYL_X: r0=" << r0 << ", theta0=" << theta0 << endl;
        }
        else if (token == "TRANSFORM_CYL_Z") {
          // usage:
          // TRANSFORM_CYL_Z <r0> <theta0-radians>
          // projects image and solution data into CYL_Z space (x,y,z)->(r,r0*theta,z)
          b_transform_cyl_z = true;
          r0 = param->getDouble(iarg++);
          const double theta0 = param->getDouble(iarg++);
          cos_theta0 = cos(theta0);
          sin_theta0 = sin(theta0);
          //if (mpi_rank == 0) cout << " > XXXXXXXXXXXXXXXXXXXXXXXXX TRANSFORM_CYL_Z: r0=" << r0 << ", theta0=" << theta0 << endl;
        }
        else if (token == "PROJECT_CYL"){
          FOR_I3 proj_cyl_xp[i] = param->getDouble(iarg++);
          FOR_I3 proj_cyl_np[i] = param->getDouble(iarg++);
          b_proj_cyl = true;

          b_geom_iso = true;
          iso_vars.push_back("comp(x_no,0)");
          iso_values.push_back(0.0);
          b_var_on_iso = true;
          var_on_iso   = "D";

          range_on_iso[0] = 0.0;
          range_on_iso[1] = 1.0;
          b_range_on_iso = true;
        }
	else if (token == "CYL_CENTER") {
	  // for use in combination with one of the "TRANSFORM_CYL_*" tokens to offset the axis...
	  cyl_center[0] = param->getDouble(iarg++);
	  cyl_center[1] = param->getDouble(iarg++);
	  cyl_center[2] = param->getDouble(iarg++);
	}
        else {
          if (mpi_rank == 0) cout << "Warning: skipping unrecognized WRITE_IMAGE token \"" << token << "\"" << endl;
        }
      }

      // and finally take a copy of the parameter...
      paramStr = param->str();

    }
    catch(int ierr) {

      return -2; // can't use this wid

    }

    return 0;

  }

};

#endif
