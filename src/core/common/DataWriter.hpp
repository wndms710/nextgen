#ifndef _DATA_WRITER_HPP_
#define _DATA_WRITER_HPP_

#include "TecplotIO.hpp"

enum TopoType {
  ISO_TOPO,
  BF_TOPO,
  BF_TOPO_EXPLICIT,
  CV_TOPO,
  LP_TOPO,
};

enum FormatType {
  TECPLOT_FORMAT,
  VTK_FORMAT,
  PVTK_FORMAT, // 3D only (reverts to VTK_FORMAT if 2D)
  ENSIGHT_FORMAT,
  FLUENT_BP_FORMAT,
  ASCII_FORMAT,
  STL_FORMAT, // iso surface geometry only (could add bf stuff)
};

enum GeomType {
  USERDEF_GEOM,
  SIMPLE_GEOM,
  ALL_GEOM,
  ALL_WITH_ZONES_GEOM, // ENSIGHT only for now
  FILE_GEOM,
};

// simple container for writing data to tecplot, vtk, ensight, etc...
class DataWriter {
public:

  string param_str; // param
  string name; // name prefix
  int interval;
  int format; // data format (tecplot, vtk, etc.) to be written
  int topo; // data topology to be written
  int geom; // data geometry to be written
  SimpleGeom* simple_geom; 
  double iso_var_value; // used for ISO_TOPO
  string iso_var_name; // used for USERDEF_GEOM ISO_TOPO
  double *iso_var; // used for ISO_TOPO (on cv's when in DUAL mode)
  bool* bfzone_flag; // used for BF_TOPO
  string sbin_filename; // used for FILE_GEOM ISO_TOPO
  int8 *cv_flag; // used for CV_TOPO
  bool b_dual; // write dual tet format instead of poly

  // used for LP_TOPO
  int lp_index;
  int *lp_flag;

  bool write_from_rank0; // shuffle data to rank 0 for write (stupid NFS)

  int ensight_casefile_loc; // used for updating ensight casefile nsteps
  int ensight_casefile_nsteps; // used by ensight data writer
  vector<string> scalar_descs_cv; // needed to remove reserved chars in ensight
  vector<string> vector_descs_cv;
  vector<string> scalar_descs_bf;
  vector<string> vector_descs_bf;
  vector<string> scalar_descs_lp;
  vector<string> vector_descs_lp;

  // these are the variables names written in a form understandable to each of data writer functions
  vector<string> scalar_names_cv;
  vector<string> vector_names_cv;
  vector<string> scalar_names_bf;
  vector<string> vector_names_bf;
  vector<string> scalar_zones_bf;
  vector<string> vector_zones_bf;
  vector<string> scalar_names_lp;
  vector<string> vector_names_lp;
  int nscalars_cv;
  int nvectors_cv;
  int nscalars_bf;
  int nvectors_bf;
  int nzones_bf;
  int nvectors_lp;
  int nscalars_lp;

  DataWriter() {
    name = "";
    geom = -1;
    simple_geom = NULL;
    topo = -1;
    interval = -1;
    format = TECPLOT_FORMAT; // default to tecplot
    ensight_casefile_loc = -1;
    ensight_casefile_nsteps = 0;
    bfzone_flag = NULL;
    nscalars_cv = 0;
    nvectors_cv = 0;
    nscalars_bf = 0;
    nvectors_bf = 0;
    nzones_bf = 0;
    sbin_filename = "";
    cv_flag = NULL;
    b_dual = false;
    iso_var = NULL;
    write_from_rank0 = false;
    lp_flag = NULL;
    lp_index = -1;
  }

  ~DataWriter() {
    DELETE(bfzone_flag);
    if (simple_geom) {
      delete simple_geom;
      simple_geom = NULL;
    }
    DELETE(cv_flag);
    DELETE(iso_var);
    DELETE(lp_flag);
  }
};

#endif
