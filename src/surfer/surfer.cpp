#include "CTI.hpp"
using namespace CTI;

// make this part of CTI eventually...
#include "IntFlag.hpp"
#include "SimpleSurface.hpp"
#include "CtiScene.hpp"
#include "KillfileReader.hpp"
#include "Histogram.hpp"
#include "LesImageMapper.hpp"
#include "SimplePoint.hpp"
#include "WebUI.hpp"

// Sept 2020: see MORPH_TRANSLATE and morphTranslateFlaggedZones(const double dx[3]) for
// recent implementation of surface morphing. Need unified cleanup.

class Surfer {

private:

  KillfileReader kfr;
  list<string> journal;

  bool b_skip_image;
  string skip_image_name;

public:

  SimpleSurface ss;

  string filetype;

  string lastTemplate;
  vector<pair<int,SimpleSurface::SubzoneData> > lastVec;

  int surface_data_type;
  string surface_data_name;

  CvPlaneMap cvPlaneMap; //persist cvs associated with a plane to speed up serial planar vis

  Surfer() {
    clear();
  }

  ~Surfer() {
    lastVec.clear();
  }

  void init() {
    COUT1("surfer::init()");
    if (mpi_rank==0) logger->setKillFilename("killsurfer");
  }

  void clear() {
    COUT1("surfer::clear()");
    surface_data_type = 0; // 0:none set,1:sp,2:st,3:sp3,4:st3
    b_skip_image = false;
    skip_image_name = "";
    surface_data_name = "";
    filetype = "unkown";
    ss.clear();
    lastVec.clear();
    srand(2); // for consistent random seeding
  }

  void run() {
    COUT1("surfer::run()");
    FOR_ALL_PARAM {
      param->touch();
      processParam(&(*param));
    }
  }

  void runInteractive(Param * param,const bool help) {

    COUT1("starting interactive session...");

    if (mpi_rank==0) logger->insertSolverAction(SIMULATION_RUN);

    // clear the webUIOutput before starting the interactive session:
    // this will prevent the deluge of INFO/WARNING messages that can occur
    // when the session starts...
    WebUI::webUIOutput.clear();

    if (param->size()>0){
      double timeout_mins = param->getDouble(0);
      kfr.setHoldWithTimeout(true,timeout_mins);
      //WebUI::webUIOutput.writeOK(true);
    }
    else{
      kfr.setHold(true);
      //WebUI::webUIOutput.writeOK(true);
    }

    while (Param * param = kfr.getParam("killsurfer")) {
      string token = MiscUtils::toUpperCase(param->getName());
      if ((token != "INTERACTIVE")&&(token != "I")) processParam(param);
    }

  }

  inline bool hasSurface() {
    return !((ss.xsp == NULL)||(ss.spost==NULL)||(ss.nst==0)||(ss.nsp==0));
  }

  void processParam(Param * param) {

    COUT1("\nprocessParam \"" << param->str() << "\"");
    const double wtime = MPI_Wtime();
    journal.push_back(param->str());

    bool b_help = false;

    string token = MiscUtils::toUpperCase(param->getName());
    if (token == "HELP") {
      b_help = true;
      if (param->size() == 0) {
        helpHelp();
        return;
      }
      token = MiscUtils::toUpperCase(param->getString(0));  // first token is keyword for help
    }

    if (token == "SURF") {
      processSurf(param,b_help);
    }
    else if (token == "CLEAR") {
      processClear(param,b_help);
    }
    else if ((token == "INTERACTIVE")||(token == "I")) {
      runInteractive(param,b_help);
    }
    else if ((token == "STOP")||(token == "STOP!")) {
      // TODO: is there a better way to do this other than throw(0)?...
      if (b_help) {
        WUI(INFO,"stop or stop! immediately exits Surfer");
      }
      else {
        WUI(INFO,"got stop or stop!: stopping surfer and exiting");
        throw(0);
      }
    }
    else if (token == "CTI_LICENSE") {
      COUT1("license location set by input file");
    }
    else if (token == "WRITE_JOURNAL") {
      processWriteJournal(param,b_help);
    }
    else if (token == "WRITE_IMAGE") {
      processWriteImage(param,b_help);
    }
    else if (token == "WRITE_JSON") {
      processWriteJson(param,b_help);
    }
    else if ((token == "SPLIT_FLAG") || (token == "SPLIT_ZONE") || (token == "SPLIT_ZONES") || (token == "SPLIT_SUBZONES")) {
      processSplitFlag(param,b_help);
    }
    else if ((token == "SUBZONES_TO_ZONES")) {
      processSubzonesToIndexedZones(param,b_help);
    }
    else if (token == "MAKE_INJECTOR") {
      processMakeInjector(param,b_help);
    }
    else if (token == "MAKE_PART") {
      processMakePart(param,b_help);
    }
    else if (token == "WINDOW_SELECT") {
      processWindowSelect(param,b_help);
    }
    else if (token == "FLIP") {
      processFlip(param,b_help);
    }
    else if (token == "ALIGN_NORMALS") {
      processAlignNormals(param,b_help);
    }
    else if (token == "MOVE_TO_ZONE") {
      processMoveToZone(param,b_help);
    }
    else if (token == "COMBINE_ZONES_WITH_SAME_NAME") {
      ss.combineZonesWithSameName();
    }
    else if (token == "RENAME_ZONE") {
      processRenameZone(param,b_help);
    }
    else if (token == "SET_ZONE_BC") {
      processSetZoneBC(param,b_help);
    }
    else if (token == "DELETE") {
      processDeleteSelectedSubzones(param,b_help);
    }
    else if (token == "TRANSLATE") {
      processTranslateSelectedSubzones(param,b_help);
    }
    else if (token == "SNAP_TO_COM") {
      processSnapToCOM(param,b_help);
    }
    else if (token == "ROTATE") {
      processRotateSelectedSubzones(param,b_help);
    }
    else if (token == "COPY") {
      processCopySelectedSubzones(param,b_help);
    }
    else if (token == "SPLIT_AFTER_COPY") {
      processSplitZonesAfterCopyAll(param,b_help);
    }
    else if (token == "SCALE") {
      processScaleSelectedSubzones(param,b_help);
    }
    else if (token == "MIRROR") {
      processMirrorSelectedSubzones(param,b_help);
    }
    else if (token == "IMPRINT") {
      processImprint(param,b_help);
    }
    else if (token == "GET_PIPE_DATA") {
      processPipes(param,b_help);
    }
    else if (token == "GET_SURFACE_CURVATURE") {
      processSurfaceCurvature(param,b_help);
    }
    else if (token == "GET_GAP_THICKNESSES") {
      processGetGapThicknesses(param,b_help);
    }
    else if (token == "MORPH_TRANSLATE") {
      processMorphTranslate(param,b_help);
    }
    else if (token == "MORPH_ROTATE_CYLINDER") {
      processMorphRotateCylinder(param,b_help);
    }
    else if (token == "MORPH_ROTATE_HALF_CYLINDER") {
      processMorphRotateHalfCylinder(param,b_help);
    }
    else if (token == "MORPH_CORNER") {
      processMorphCorner(param,b_help);
    }
    else if (token == "MORPH_STRETCH") {
      processMorphStretch(param,b_help);
    }
    else if (token == "MORPH_CYLINDER") {
      processMorphCylinder(param,b_help);
    }
    else if (token == "MORPH_CYLINDER_R") {
      processMorphCylinderR(param,b_help);
    }
    else if (token == "MORPH_GAUSSIAN_BUMP") {
      processMorphGaussianBump(param,b_help);
    }
    else if (token == "REFINE_TRIS") {
      processRefineTris(param,b_help);
    }
    else if (token == "APPLY_SIMILAR") {
      processApplySimilar(param,b_help);
    }
    else if (token == "SELECT_SIMILAR") {
      processSelectSimilar(param,b_help);
    }
    else if (token == "SUGGEST_MESH") {
      processSuggestMesh(param,b_help);
    }
    else if (token == "WRITE_STITCH_INPUT") {
      processWriteStitch(param,b_help);
    }
    else if (token == "SHOW_STITCH_INPUT") {
      processShowStitch(param,b_help);
    }
    else if (token == "WRITE_CHARLES_INPUT") {
      processWriteCharles(param,b_help);
    }
    else if (token == "SHOW_CHARLES_INPUT") {
      processShowCharles(param,b_help);
    }
    else if (token == "EDGE_LENGTH_PDF") {
      processEdgeLengthPdf(param,b_help);
    }
    else if (token == "RUN_DIAGNOSTICS") {
      processRunDiagnostics(param,b_help);
    }
    else if ((token == "ENSURE_TEOST")||(token == "BUILD_TEOST")) {
      processEnsureTeost(param,b_help);
    }
    else if (token == "ENSURE_OPEN_EDGE_GROUPS") {
      processEnsureOpenEdgeGroups(param,b_help);
    }
    else if (token == "REGROUP_OPEN_EDGES") {
      processRegroupOpenEdges(param,b_help);
    }
    else if (token == "CALC_GCL") {
      processCalcGcl(param,b_help);
    }
    else if (token == "DELETE_TRIS_WITH_COLLOCATED_NODES") {
      processDeleteTrisWithCollocatedNodes(param,b_help);
    }
    else if (token == "MERGE_COLLOCATED_NODES") {
      processMergeCollocatedNodes(param,b_help);
    }
    else if (token == "CLOSE_OPEN_EDGE_LOOPS") {
      processCloseOpenEdgeLoops(param,b_help);
    }
    else if (token == "ZIP_OPEN_EDGES") {
      processZipOpenEdges(param,b_help);
    }
    else if (token == "ZIP_OPEN_EDGES_CHT") {
      processZipOpenEdgesCht(param,b_help);
    }
    else if (token == "MERGE_OPEN_EDGE_NODES") {
      processMergeOpenEdgeNodes(param,b_help);
    }
    else if (token == "SET_FEATURE_ANGLE") {
      processSetFeatureAngle(param,b_help);
    }
    else if (token == "GROUP_DYNAMIC_EDGES") {
      processGroupDynamicEdges(param,b_help);
    }
    else if (token == "DISCRETIZE_EDGES") {
      processDiscretizeEdges(param,b_help);
    }
    else if (token == "ENSURE_DELAUNAY") {
      processEnsureDelaunayOnZones(param,b_help);
    }
    else if (token == "DISCRETIZE_ZONES") {
      processDiscretizeZones(param,b_help);
    }
    else if (token == "ROUGHEN") {
      processRoughenSelectedSubzones(param,b_help);
    }
    else if (token == "SET_MATERIAL") {
      processSetMaterial(param);
    }
    else if (token == "CLOSE_EDGE_LOOPS") {
      processCloseEdgeLoops(param,b_help);
    }
    else if (token == "DELETE_ISOLATED_NODES") {
      processDeleteIsolatedNodes(param,b_help);
    }
    else if (token == "REPAIR_LINEAR_TRIS") {
      processRepairLinearTris(param,b_help);
    }
    else if (token == "REPAIR_MULTI_NBR_TRIS") {
      processRepairMultiNeighborTris(param,b_help);
    }
    else if (token == "SPLIT_MULTIEDGES") {
      ss.splitMultiEdges();
    }
    else if (token == "SHOW_GEODESIC_DISTANCE") {
      processShowGeoDistFromPoint(param,b_help);
    }
    else if (token == "SLICE_GEN") {
      processSliceGen(param,b_help);
    }
    else if ((token == "WRITE_SBIN")||(token == "WRITE_BIN")) {
      processWriteSbin(param,b_help);
    }
    else if (token == "WRITE_STL") {
      processWriteStl(param,b_help);
    }
    else if (token == "WRITE_TECPLOT") {
      processWriteTecplot(param,b_help);
    }
    else if (token == "SET_PERIODIC") {
      processSetPeriodic(param,b_help);
    }
    else if (token == "CLEAR_PERIODIC") {
      processClearPeriodic(param,b_help);
    }
    else if (token == "DUMP_ZONES") {
      processDumpZones(param,b_help);
    }
    else if (token == "CLEAR_ST_FLAG") {
      ss.st_flag.setLength(ss.nst);
      ss.st_flag.setAll(0);
    }
    else if (token == "REPORT_ZONE_DATA") {
      processReportZoneData(param,b_help);
    }
    else if (token == "SET_ST_FLAG") {
      processSetStFlag(param,b_help);
    }
    else if (token == "WRITE_FLUX_PROBES_CYLINDER") {
      processWriteFluxProbesCylinder(param,b_help);
    }
    else if (token == "MULTI_FLUX_PROBE") {
      processMultiFluxProbe(param,b_help);
    }
    else if (token == "CALC_CLOSEST_CYL") {
      processCalcClosestCylinder(param,b_help);
    }
    else if (token == "SCALE_OPEN_EDGES") {
      processScaleOpenEdges(param,b_help);
    }
    else if (token == "TRANSLATE_OPEN_EDGES") {
      processTranslateOpenEdges(param,b_help);
    }
    else if (token == "ROTATE_OPEN_EDGES") {
      processRotateOpenEdges(param,b_help);
    }
    else if (token == "FORCE_ORTHOGONALITY_OPEN_EDGES") {
      processForceOrthogonalityOpenEdges(param,b_help);
    }
    else if (token == "NSMOOTH") {
      processNsmooth(param,b_help);
    }
    else if (token == "INTERSECT") {
      processIntersect(param,b_help);
    }
    else if (token == "HOLE_PAIRING") {
      processHolePairing(param, b_help);
    }
    else if (token == "TEST_SURFACE_RETRI") {
      // testing surface re-triangulation using the 2d voronoi diagram built
      // from cutting the surface...
      ss.retri();
    }
    else if ((token == "JOURNAL")||(token == "HISTORY")) {
      processHistory(param,b_help);
    }
    else if (token == "REWIND") {
      processRewind(param,b_help);
    }
    else if (token == "UNDO") {
      processUndo(param,b_help);
    }
    else if (token == "FLAG_TRIS_TOUCHING_MULTIEDGES") {
      if (param->size() == 0) {
        ss.flagTrisTouchingMultiEdges();
      }
      else if (param->size() == 1) {
        // expecting a comma-delimited list of multiedge indices...
        set<int> meSet;
        MiscUtils::splitCsv(meSet,param->getString());
        ss.flagTrisTouchingMultiEdges(meSet);
      }
      else {
        cout << "unsupported size: " << param->size() << endl;
      }
    }
    else if (token == "FLAG_NON_MANIFOLD_TRIS") {
      if (!ss.hasNonManifoldData()) {
        WUI(WARN,"Please call RUN_DIAGNOSTICS with selected params before calling FLAG_NON_MANIFOLD_TRIS. Returning.");
        helpRunDiagnostics();
      }
      else {
        if (param->size() == 0) ss.flagNonManifoldTris();
        else {
          int nme_bits = 0;
          int iarg = 0;
          while (iarg < param->size()) {
            const string token = MiscUtils::toUpperCase(param->getString(iarg++));
            if (token == "TYPE" || token == "TYPES") {
              vector<string> nmeType;
              MiscUtils::splitCsv(nmeType,param->getString(iarg++));
              for (int i=0,end=nmeType.size(); i<end; ++i) {
                // should correspond to enum in NonManifoldTriFlag.hpp
                // NM_MULTI_EDGE_ADJACENT = (1<<0),
                // NM_COLLAPSED = (1<<1),
                // NM_LINEAR = (1<<2),
                // NM_SHARED_EDGES = (1<<3),
                // NM_INTERSECT = (1<<4),
                // NM_OVERLAP = (1<<5),
                // NM_IMPRINT = (1<<6),
                if (nmeType[i] == "MULTI_EDGE_ADJACENT") nme_bits |= 1<<0;
                else if (nmeType[i] == "COLLAPSED") nme_bits |= 1<<1;
                else if (nmeType[i] == "LINEAR") nme_bits |= 1<<2;
                else if (nmeType[i] == "SHARED_EDGES") nme_bits |= 1<<3;
                else if (nmeType[i] == "INTERSECT") nme_bits |= 1<<4;
                else if (nmeType[i] == "OVERLAP") nme_bits |= 1<<5;
                else if (nmeType[i] == "IMPRINT") nme_bits |= 1<<6;
                else {
                  WUI(WARN,"unrecognized non-manifold tri type: " << nmeType[i] << "; skipping");
                }
              }
              if (nme_bits) ss.flagNonManifoldTris(nme_bits);
              else {
                WUI(WARN,"no tris of the specified type(s) were found; skipping");
              }
            }
            else {
              WUI(WARN,"unrecognized nme type token: " << token << "; skipping");
            }
          }
        }
      }
    }
    else if (token == "FLAG_TRIS_TOUCHING_OPEN_EDGES") {
      ss.flagTrisTouchingOpenEdges();
    }
    else if (token == "FLAG_TRIS_TOUCHING_FLAGGED_TRIS") {
      if (param->size() == 0) ss.flagTrisTouchingFlaggedTris();
      else {
        int iarg = 0;
        while (iarg < param->size()) {
          const string token = MiscUtils::toUpperCase(param->getString(iarg++));
          if (token == "SAME_ZONE" || token == "SAME_ZONES") ss.flagTrisTouchingFlaggedTrisInSameZone();
          else {
            WUI(WARN,"unrecognized flag tris token: " << token << "; skipping");
          }
        }
      }
    }
    else if (token == "FLAG_TRIS_WITH_AREA_LESS_THAN") {
      if (param->size() == 0) {
        cout << "expecting FLAG_TRIS_WITH_AREA_LESS_THAN <double>" << endl;
      }
      else {
        const double area = param->getDouble();
        ss.flagTrisWithAreaLessThan(area);
      }
    }
    else if (token == "FLAG_TRIS_WITH_SUBZONE_AREA_LESS_THAN") {
      if (param->size() == 0) {
        cout << "expecting FLAG_TRIS_WITH_SUBZONE_AREA_LESS_THAN <double>" << endl;
      }
      else {
        const double area = param->getDouble();
        ss.flagTrisWithSubzoneAreaLessThan(area);
      }
    }
    else if (token == "DELETE_FLAGGED_TRIS") {
      ss.deleteFlaggedTris();
    }
    else if (token == "WRITE_FLAGGED_TRIS_TECPLOT") {
      bool b_x0 = false;
      double x0[3];
      string filename = "flagged_tris.dat";
      int iarg = 0;
      while (iarg < param->size()) {
        token = param->getString(iarg++);
        if (token == "X0") {
          b_x0 = true;
          FOR_I3 x0[i] = param->getDouble(iarg++);
        }
        else {
          // assume this is the name:
          filename = token;
        }
      }
      if (b_x0) {
        ss.writeSelectedFacesByZoneToTecplot(filename,ss.st_flag,x0);
      }
      else {
        ss.writeSelectedFacesByZoneToTecplot(filename,ss.st_flag);
      }
    }
    else if (token == "FLAG_TRIS") {
      processFlagTris(param,b_help);
    }
    else if (token == "FLAG_TRIS_TOUCHING") {
      processFlagTrisTouching(param,b_help);
    }
    else if (token == "COUNT_FLAGGED_TRIS") {
      ss.countFlaggedTris();
    }
    else if (token == "INVERT_FLAGGED_TRIS") {
      ss.invertFlaggedTris();
    }
    else if (token == "DELETE_TRIS_WITH_IDENTICAL_NODES") {
      ss.deleteTrisWithIdenticalNodes();
    }
    else if (token == "SEPARATE_OVERLAPPING_TRIS") {
      processSeparateOverlappingTris(param,b_help);
    }
    else if (token == "DEBUG_TRI_TRI") {
      processDebugTriTri(param,b_help);
    }
    else if (token == "QUERY") {
      processQuery(param,b_help);
    }
    else if (token == "FLATTEN_ROTOR") {
      processFlattenRotor(param,b_help);
    }
    else if (token == "TEST_CHT") {
      processTestCht(param,b_help);
    }
    else if (token == "THREE_POINTS") {
      processThreePoints(param,b_help);
    }
    else if (token == "INTERFACE") {
      processInterface(param,b_help);
    }
    else if (token == "MAKE_FF") {
      processMakeFF(param,b_help);
    }
    //-------------------
    // below methods are candidates for deprecation/removal
    //-------------------
    else if (token == "FLIP_TRIS_RANDOM") {
      processFlipTrisRandom(param,b_help);
    }
    else if (token == "ADD_OPEN_EDGE_LOOP") {
      processAddOpenEdgeLoop(param,b_help);
    }
    else if (token == "TEST_SEED") {
      ss.testSeed(param,b_help);
    }
    else if (token == "BBOX") {
      if (b_help) {
        WUI(INFO,"BBOX reports the boundary box of the current surface(s)");
      }
      else {
        double bbox[6];
        ss.getBoundingBox(bbox);
        WUI(INFO,"bounding box x: " << bbox[0] << " " << bbox[1] << " y: " << bbox[2] << " " << bbox[3] << " z: " << bbox[4] << " " << bbox[5] <<
            "\nbounding box dimensions dx: " << bbox[1]-bbox[0] << " dy: " << bbox[3]-bbox[2] << " dz: " << bbox[5]-bbox[4]);
      }
    }
    else {
      if (b_help) {
        WUI(WARN,"HELP not available for unrecognized param: " << token);
      }
      else {
	if ((token == "VERBOSE") || (token == "DEFINE")) {
	  COUT1(" > processed");
	}
	else {
	  WUI(WARN,"skipping unrecognized param: " << token);
	}
      }
    }

    // do a little checking...
    if (ss.nsp || ss.nst) ss.check();

    // report timing...
    COUT1("processParam elapsed time " << MPI_Wtime()-wtime << " [s]\n");
  }

  void helpHelp() {

    WUI(WARN,"HELP must be followed by a parameter. HELP is available for\n" <<
        "  SURF           CLEAR           SPLIT_FLAG      FLIP\n" <<
        "  ALIGN_NORMALS  MOVE_TO_ZONE    WRITE_IMAGE     RENAME_ZONE\n" <<
        "  DELETE         TRANSLATE       SNAP_TO_COM     ROTATE\n" <<
        "  COPY           SCALE           MIRROR          IMPRINT\n" <<
        "  REFINE_TRIS    RUN_DIAGNOSTICS ZIP_OPEN_EDGES  WRITE_SBIN\n" <<
        "  WRITE_STL      WRITE_TECPLOT   SET_PERIODIC    CLEAR_PERIODIC\n" <<
        "  CLOSE_OPEN_EDGE_LOOPS WRITE_CHARLES_INPUT\n" <<
        "  WRITE_JOURNAL  WRITE_JSON      WRITE_STITCH_INPUT");

  }

  void processClear(Param * param,const bool help) {
    UNUSED(param);
    if (help) {
      WUI(INFO,"removes all surfaces from the active workspace. there are no associated parameters\n");
      return;
    }

    try {
      clear();
    }
    catch(int e) {
      if (e == 0) {
        WUI(WARN,"I don't know how, but you failed in CLEAR");
      }
      else {
        WUI(WARN,"I don't know how, but you failed in CLEAR");
      }
    }
  }

  void processSurf(Param * param,const bool help) {

    if (help) {
      // allow request for specific help: e.g. HELP SURF GRID
      if (param->size() > 1)
        helpSurf(MiscUtils::toUpperCase(param->getString(1)));
      else
        helpSurf();
      return;
    }

    try {
      bool isEmpty = false;
      int ierr = 0;

      if (param->size() == 0) {
        // no params specified...
        isEmpty = true;
        CWARN("single entry SURF param command: " << param->str());
      }

      ierr = ss.init(param);
      if (!isEmpty) filetype = "surface";

      if (ierr == 0) {
        assert(ss.znost);
        WebUI::webUIOutput.ensureImage();
      }
    }
    catch(int e) {
      if (e == 0) {
        WUI(WARN,"Expecting SURF <file-type> <file-name>");
      }
      else {
        WUI(WARN,"Expecting SURF <file-type> <file-name>");
      }
      helpSurf();
    }
  }

  void helpSurf() {
    WUI(INFO,"SURF imports a faceted surface from a file. Subsequent commands will append the new surface to the existing surface.\n" <<
        "Supported file formats include: SBIN, STL, MSH/CAS, PLY, OBJ. Examples:\n"
        "  SURF SBIN awesome_car.sbin\n" <<
        "  SURF STL pyramids.stl\n" <<
        "SURF can also be used to create or add geometric primitives such as spheres and cylinders. For example:\n"
        "  SURF SPHERE <x> <y> <z> <r>\n"
        "  SURF PIPE <x0> <y0> <z0> <x1> <y1> <z1> <r>\n"
        "For more info see [$CWIKB:surfer_import]"
        );
  }

  void helpSurf(const string& token) {
    if (token == "GRID") {
      WUI(INFO,"SURF GRID creates a Cartesian array of spheres suitable for generating turbulence near an inlet. Example:\n" <<
          "  SURF GRID NAME=grid D=0.002 L=0.005 X=1.1 0 0 N=1 0 0 SUBZONES 7,8\n" <<
          "If the SUBZONES token is omitted, a grid will be added inside all inlets cut by the specified plane.\n" <<
          "Normally you should either delete spheres touching the boundaries or use the INTERSECT command to produce a\n"
          "single watertight surface including the spheres.");
    }
    else if (token == "LOFTED") {
      WUI(INFO,"SURF LOFTED creates a lofted surface from curves. Example:\n" <<
          "  SURF LOFTED NL 10 PROFILES file1 file2 [file3...]\n" <<
          "where the files are ASCII files of 3D points.");
    }
    else {
      WUI(WARN,"Specific SURF help not available for " << token);
      helpSurf();
    }
  }

  void processWriteJournal(Param * param,const bool b_help) {

    if (b_help) {
      helpWriteJournal();
    }
    else {

      if (mpi_rank == 0) {
        assert(!journal.empty());
        string name = "journal.in"; // TODO: if this name exists, don't overwrite. Turn this to journal(1).in or something
        bool b_verbose = false;
        int iarg = 0;
        while (iarg < param->size()) {
          string token = param->getString(iarg++);
          if ((token == "VERBOSE")||(token == "verbose")) {
            b_verbose = true;
          }
          else {
            // assume this is the name...
            name = token;
          }
        }
        FILE * fp = fopen(name.c_str(),"w");
        int count = 0;
        for (list<string>::iterator iter = journal.begin(); iter != journal.end(); ++iter) {
          // default is to write active commands only...
          if (b_verbose||isActiveCommand(*iter)) {
            fprintf(fp,"%s\n",iter->c_str());
            ++count;
          }
        }
        fclose(fp);
        WUI(INFO,count << " journal params written to file " << name);
      }
    }

  }

  void helpWriteJournal() const {
    WUI(INFO,
        "WRITE_JOURNAL exports command history to file. Examples:\n" <<
        "WRITE_JOURNAL # exports to journal.in\n" <<
        "WRITE_JOURNAL my_journal.in\n" <<
        "WRITE_JOURNAL VERBOSE # includes ALL commands");
  }

  void processWriteSbin(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpWriteSbin();
      return;
    }

    int iarg = 0;
    string filename = "surface.sbin";
    vector<int> subzone_indices;
    bool b_st_flag = false;
    string first_arg;
    while (iarg < param->size()) {
      if (iarg == 0) first_arg = param->getString(iarg);
      string token = MiscUtils::toUpperCase(param->getString(iarg++));

      if (parseSpecifiedSubzones(subzone_indices,token,param,iarg)) {
        //no op; parsed info
      }
      else if (token == "NAME") {
        filename = param->getString(iarg++);
      }
      else if (token == "ST_FLAG") {
        b_st_flag = true;
      }
      else {
        // if empty and was first string, assume is filename
        if (iarg == 1) filename = first_arg;
      }
    }

    if (int(subzone_indices.size()) == ss.nsz) subzone_indices.clear();  // so that ALL also writes periodic info


    if (!subzone_indices.empty()) {
      ss.writeSubzonesToBinary(filename,subzone_indices);
    }
    else if (b_st_flag) {
      ss.writeSelectedTrisToBinary(filename,ss.st_flag);
    }
    else {
      // no subzones were specified so assume entire surface
      ss.writeBinary(filename);
    }

  }
  void helpWriteSbin() {
    WUI(INFO,"writes specified portions of the surface to a Cascade binary file format\n" <<
        "  examples:\n" <<
        "    WRITE_SBIN my_surface.sbin \n" <<
        "    WRITE_SBIN ZONE_NAMES pipe,nozzle,fins NAME model_a.sbin\n" <<
        "  more info: [$CWIKB:surfer_export]"
        );
  }

  void processWriteStl(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpWriteStl();
      return;
    }

    int iarg = 0;
    string name = "surface.stl";
    vector<int> subzone_indices;
    ss.st_flag.resize(ss.nst);
    ss.st_flag.setAll(0);
    bool b_single_zone = false;

    string first_arg;
    while (iarg < param->size()) {
      if (iarg == 0) first_arg = param->getString(iarg);
      const string token = MiscUtils::toUpperCase(param->getString(iarg++));

      if (parseSpecifiedSubzones(subzone_indices,token,param,iarg)) {
        //no op; parsed info
      }
      else if (token == "NAME") {
        name = param->getString(iarg++);
      }
      else if (token == "SINGLE_ZONE") {
        b_single_zone = true;
      }
      else {
        // if empty and was first string, assume is filename
        if (iarg == 1) name = first_arg;
        else {
          WUI(WARN,"unrecognized STL token \"" << token << "\"; skipping");
        }
      }
    }

    if (!subzone_indices.empty()) {
      ss.flagTrisFromSubzoneVec(subzone_indices);  // flagtris appropriately
    }
    else {
      // default is to flag all tris
      for (int ist=0; ist < ss.nst; ++ist) ss.st_flag[ist] = 1;
    }

    ss.writeSelectedTrisToStl(name,b_single_zone);
  }
  void helpWriteStl() {
    WUI(INFO,"writes specified portions of the surface to an ASCII STL file\n" <<
        "  examples:\n" <<
        "    WRITE_STL my_surface.stl \n" <<
        "    WRITE_STL ZONE_NAMES pipe,nozzle,fins NAME model_a.stl\n" <<
        "  more info: [$CWIKB:surfer_export]"
        );
  }

  void processWriteTecplot(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpWriteTecplot();
      return;
    }

    int iarg = 0;
    string name = "surface.dat";
    vector<int> zone_indices;

    string first_arg;  // so not harmed by uppercase
    while (iarg < param->size()) {
      if (iarg == 0) first_arg = param->getString(iarg);
      const string token = MiscUtils::toUpperCase(param->getString(iarg++));

      if (parseSpecifiedZones(zone_indices,token,param,iarg)) {
        //no op; parsed info
      }
      else if (token == "NAME") {
        name = param->getString(iarg++);
      }
      else {
        // if empty and was first string, assume is filename
        if (iarg == 1) name = first_arg;
      }
    }

    if (!zone_indices.empty()) {
      ss.flagZonesFromZoneIndexVec(zone_indices);
      ss.writeFlaggedZonesToTecplot(name);
    }
    else ss.writeTecplot(name);
  }
  void helpWriteTecplot() {
    WUI(INFO,"writes specified portions of the surface to a Tecplot ascii .dat file\n" <<
        "  examples:\n" <<
        "    WRITE_TECPLOT my_surface.dat \n" <<
        "    WRITE_TECPLOT ZONE_NAMES pipe,nozzle,fins NAME model_a.dat\n" <<
        "  more info: [$CWIKB:surfer_export]"
        );
  }

  void processSetPeriodic(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpSetPeriodic();
      return;
    }

    bool got_periodic_transform = false;
    PeriodicTransform pt;

    bool got_zones = false;
    vector<string> zone0NameVec,zone1NameVec;

    bool got_zone_ids = false;
    bool got_edges = false;
    vector<int> zone0idVec,zone1idVec;

    bool b_force = false; // new "force" mode walks graph of edges between corners and loops
    double crease_angle = 120.0;

    // NOTE: probably should replace all periodicitiy setting with this force mode eventually
    // but for now make it optional to handle particularly tough cases where the current local
    // checks fail...

    int iarg = 0;
    while (iarg < param->size()) {
      const string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if ((token == "ZONE")||(token == "ZONES")) {
        // ZONES is followed by 2 csv params, e.g. ZONES x0,t0 x1,t1
        if (iarg > param->size()-2) {
          WUI(WARN,"SET_PERIODIC ZONES: expecting two sets of comma-separated zone names; skipping");
          return;
        }
        else if (got_edges) {
          WUI(WARN,"SET_PERIODIC connot currently use both EDGES and ZONES simultaneously");
          return;
        }
        MiscUtils::splitCsv(zone0NameVec,param->getString(iarg++));
        MiscUtils::splitCsv(zone1NameVec,param->getString(iarg++));
        got_zones = true;
      }
      else if (token == "ZONE_IDS") {
        if (iarg > param->size()-2) {
          CWARN("SET_PERIODIC ZONE_IDS: expecting two sets of comma-separated zone ids; skipping");
          return;
        }
        else if (got_edges) {
          WUI(WARN,"SET_PERIODIC connot currently use both EDGES and ZONES simultaneously");
          return;
        }
        MiscUtils::splitCsv(zone0idVec,param->getString(iarg++));
        MiscUtils::splitCsv(zone1idVec,param->getString(iarg++));
        got_zone_ids = true;
      }
      else if (token == "EDGES") {
        if (iarg > param->size()-2) {
          CWARN("SET_PERIODIC EDGES: expecting two sets of comma-separated open edge ids; skipping");
          return;
        }
        else if (got_zones || got_zone_ids) {
          WUI(WARN,"SET_PERIODIC connot currently use both EDGES and ZONES or ZONE_IDS simultaneously");
          return;
        }
        MiscUtils::splitCsv(zone0idVec,param->getString(iarg++));
        MiscUtils::splitCsv(zone1idVec,param->getString(iarg++));
        got_edges = true;
      }
      else if (token == "CART") {
        // CART <dx> <dy> <dz>
        if (iarg > param->size()-3) {
          cout << "ERROR: SET_PERIODIC CART: expecting 3 components of CART transformation. Check syntax." << endl;
          return;
        }
        double dx[3]; FOR_I3 dx[i] = param->getDouble(iarg++);
        pt.setCart(dx);
        got_periodic_transform = true;
      }
      else if (token == "CYL_X") {
        // CYL_X <degrees>
        if (iarg > param->size()-1) {
          cout << "ERROR: SET_PERIODIC CYL_X: expecting degrees of CYL_X transformation. Check syntax." << endl;
          return;
        }
        const double degrees = param->getDouble(iarg++);
        pt.setCylx(degrees);
        got_periodic_transform = true;
      }
      else if (token == "CYL_Y") {
        // CYL_Y <degrees>
        if (iarg > param->size()-1) {
          cout << "ERROR: SET_PERIODIC CYL_Y: expecting degrees of CYL_Y transformation. Check syntax." << endl;
          return;
        }
        const double degrees = param->getDouble(iarg++);
        pt.setCyly(degrees);
        got_periodic_transform = true;
      }
      else if (token == "CYL_Z") {
        // CYL_Z <degrees>
        if (iarg > param->size()-1) {
          cout << "ERROR: SET_PERIODIC CYL_Z: expecting degrees of CYL_Z transformation. Check syntax." << endl;
          return;
        }
        const double degrees = param->getDouble(iarg++);
        pt.setCylz(degrees);
        got_periodic_transform = true;
      }
      else if (token == "SUGGEST") {
        pt.setKind(PERIODIC_TRANSFORM_NULL);
        got_periodic_transform = true;
      }
      else if (token == "FORCE") {
        b_force = true;
      }
      else if (token == "CREASE_ANGLE") {
        crease_angle = param->getDouble(iarg++);
      }
      else {
        cout << "ERROR: SET_PERIODIC: unrecognized token \"" << token << "\". Check syntax." << endl;
        return;
      }
    }

    // use the subzones to build a teost-like list of periodic edges...

    if (!got_zones && !got_zone_ids && !got_edges) {
      cout << "ERROR: SET_PERIODIC: missing ZONES, ZONE_IDS, or EDGES. Check syntax." << endl;
      return;
    }

    if (!got_periodic_transform) {
      cout << "ERROR: SET_PERIODIC: missing transform definition: CART <dx> <dy> <dz>, CYL_X <degrees>, CYL_Y <degrees>, or CYL_Z <degrees>. Check syntax." << endl;
      return;
    }

    int ierr = 0;
    if (got_zones || got_zone_ids) {
      ss.zone_flag.setLength(ss.zoneVec.size());
      ss.zone_flag.setAll(0);

      if (got_zones) {
        for (int ii = 0, ii_end=zone0NameVec.size(); ii < ii_end; ++ii) {
          const int izone = ss.getZoneIndex(zone0NameVec[ii]);
          if (izone == -1) {
            cout << "ERROR: SET_PERIODIC: cannot find zone named \"" << zone0NameVec[ii] << "\"" << endl;
            ierr = 1;
            continue;
          }
          assert((izone >= 0)&&(izone < int(ss.zoneVec.size())));
          if (ss.zoneVec[izone].getPeriodicBits() != 0) ierr = 2;
          ss.zone_flag[izone] = 1;
        }
        for (int ii = 0, ii_end=zone1NameVec.size(); ii < ii_end; ++ii) {
          const int izone = ss.getZoneIndex(zone1NameVec[ii]);
          if (izone == -1) {
            cout << "ERROR: SET_PERIODIC: cannot find zone named \"" << zone1NameVec[ii] << "\"" << endl;
            ierr = 1;
            continue;
          }
          assert((izone >= 0)&&(izone < int(ss.zoneVec.size())));
          if (ss.zoneVec[izone].getPeriodicBits() != 0) ierr = 2;
          ss.zone_flag[izone] = 2;
        }
      }
      else if (got_zone_ids) {
        for (int ii = 0, ii_end=zone0idVec.size(); ii < ii_end; ++ii) {
          const int izone = zone0idVec[ii];
          assert((izone >= 0)&&(izone < int(ss.zoneVec.size())));
          if (ss.zoneVec[izone].getPeriodicBits() != 0) ierr = 2;
          ss.zone_flag[izone] = 1;
        }
        for (int ii = 0, ii_end=zone1idVec.size(); ii < ii_end; ++ii) {
          const int izone = zone1idVec[ii];
          assert((izone >= 0)&&(izone < int(ss.zoneVec.size())));
          if (ss.zoneVec[izone].getPeriodicBits() != 0) ierr = 2;
          ss.zone_flag[izone] = 2;
        }
      }

      if (ierr) {
        COUT1(" > Periodicity already set for (some of) these zones. To reapply periodicity run CLEAR_PERIODIC first.");
        if (ierr == 2) {
          WUI(WARN,"Probably a simple oversight, but periodicity is already set for some of these zones. Clear the current periodicity before attempting to reapply");
        }
        return;
      }

      ierr = ss.setPeriodicFlaggedZones(pt,b_force,crease_angle);

    }
    else {
      assert(got_edges);
      ss.ensureOpenEdgeGroups();
      // convert from UI edge index to 0-indexed group
      for (vector<int>::iterator it=zone0idVec.begin(),end=zone0idVec.end(); it!=end; ) {
        if (*it < OPEN_E_SZ_MIN || *it > OPEN_E_SZ_MAX) {
          WUI(WARN,"removing invalid open edge group " << *it << " from list");
          it = zone0idVec.erase(it);
        }
        else {
          *it -= OPEN_E_SZ_MIN;
          ++it;
        }
      }
      for (vector<int>::iterator it=zone1idVec.begin(),end=zone1idVec.end(); it!=end; ) {
        if (*it < OPEN_E_SZ_MIN || *it > OPEN_E_SZ_MAX) {
          WUI(WARN,"removing invalid open edge group " << *it << " from list");
          it = zone1idVec.erase(it);
        }
        else {
          *it -= OPEN_E_SZ_MIN;
          ++it;
        }
      }

      if (zone0idVec.empty() || zone1idVec.empty()) {
        WUI(ERR,"no valid edges were specified on one of the periodic sides; skipping");
        return;
      }

      ierr = ss.setPeriodicOpenEdgeGroups(zone0idVec,zone1idVec,pt);
    }

    if (ierr == 0) {
      WUI(INFO,"Things that leave will now come back! Periodicity set successfully");
    }
    else if (ierr == 2) {}  // user suggestion
    else {
      WUI(ERR,"Bad news friend, setting periodicity failed. We took evasive action, a side effect of which required resetting all periodic connectivity. Please re-apply all previously applied periodicities.");
      ss.clearPeriodicity();
    }
  }
  void helpSetPeriodic() {
    WUI(INFO,"sets periodic transformation(s) for a single surface\n" <<
        "  examples:\n" <<
        "    SET_PERIODIC ZONES per_a0,per_a1 per_b0,per_b1 CART 0 0 5 \n" <<
        "    SET_PERIODIC ZONES perA per_B CYL_X 30.0 \n" <<
        "    SET_PERIODIC ZONES perA per_B SUGGEST (get a suggested periodic transform for the provided zones) \n" <<
        "  more info: [$CWIKB:surfer_periodicity]"
        );
  }

  void processClearPeriodic(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    UNUSED(param);
    if (help) {
      WUI(INFO,"clears all applied periodic transformations\n" <<
          "  more info: [$CWIKB:surfer_periodicity]"
          );
      return;
    }

    ss.clearPeriodicity();
  }

  void processDumpZones(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    UNUSED(param);
    if (help) {
      WUI(INFO,"writes zone information & metadata to stdout\n"
          );
      return;
    }

    ss.dumpZones();
  }

  void processReportZoneData(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"writes zone or subzone information & metadata to stdout\n" <<
          "  examples:\n" <<
          "    REPORT_ZONE_DATA SUBZONES (data is written per-subzone)"
          );
      return;
    }

    int iarg = 0;
    bool b_subzone = false;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "SUBZONES" || token == "SUBZONE") b_subzone = true;
      else {
        WUI(WARN,"unrecognized REPORT_ZONE_DATA token: " << token << "; skipping");
      }
    }
    ss.reportZoneData(b_subzone);
  }

  void processSetStFlag(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"set the internal tri flag to specific portions of the surface\n"
          );
      return;
    }

    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "ZONES") {
        // select zones based on zone name provided in comma-delimited list
        vector<string> zoneNameVec;
        MiscUtils::splitCsv(zoneNameVec,param->getString(iarg++));
        ss.zone_flag.setLength(ss.zoneVec.size());
        ss.zone_flag.setAll(0);
        assert(0); // need to finish this one
      }
      else if (token == "ZONE_IDS") {
        // select zones based on zone index provided in comma-delimited list
        vector<int> zoneIdVec;
        MiscUtils::splitCsv(zoneIdVec,param->getString(iarg++));
        ss.zone_flag.setLength(ss.zoneVec.size());
        ss.zone_flag.setAll(0);
        const int nzn = ss.zoneVec.size();
        for (int ii = 0, ii_end=zoneIdVec.size(); ii < ii_end; ++ii) {
          const int izone = zoneIdVec[ii];
          if ((izone < 0)||(izone >= nzn)) {
            cout << "Warning: SET_ST_FLAG ZONE_IDS out of range: " << izone << endl;
          }
          else {
            cout << " > ZONE_IDS: " << izone << ", selecting st's in zone: " << ss.zoneVec[izone].getName() << endl;
            ss.zone_flag[izone] = 1;
          }
        }
        for (int ist = 0; ist < ss.nst; ++ist) {
          const int izone = ss.znost[ist];
          if (ss.zone_flag[izone] == 1) ss.st_flag[ist] = 1;
        }
      }
      else if (token == "TRI_IDS") {
        // select tris based on tri index provided in comma-delimitted list
        vector<int> triIdVec;
        MiscUtils::splitCsv(triIdVec,param->getString(iarg++));
        for (int ii = 0, ii_end=triIdVec.size(); ii < ii_end; ++ii) {
          const int ist = triIdVec[ii];
          if ((ist < 0)||(ist >= ss.nst)) {
            cout << "Warning: SET_ST_FLAG TRI_IDS out of range: " << ist << endl;
          }
          else {
            ss.st_flag[ist] = 1;
          }
        }
      }
      else if (token == "MARCH_ZONES") {
        // select zones adjacent to existing flagged tris...
        ss.ensureTeost();
        ss.zone_flag.setLength(ss.zoneVec.size());
        ss.zone_flag.setAll(0);
        for (int ist = 0; ist < ss.nst; ++ist) {
          if (ss.st_flag[ist] == 1) {
            FOR_I3 {
              int ist_nbr,i_nbr,orient_nbr;
              if (ss.getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
                if (ss.st_flag[ist_nbr] == 0) {
                  const int izone_nbr = ss.znost[ist_nbr];
                  ss.zone_flag[izone_nbr] = 1;
                }
              }
            }
          }
        }
        for (int izone = 0, nzn=ss.zoneVec.size(); izone < nzn; ++izone) {
          if (ss.zone_flag[izone] == 1) cout << " > MARCH_ZONES: " << izone << ", selecting st's in zone: " << ss.zoneVec[izone].getName() << endl;
        }
        for (int ist = 0; ist < ss.nst; ++ist) {
          const int izone = ss.znost[ist];
          if (ss.zone_flag[izone] == 1) ss.st_flag[ist] = 1;
        }
      }
      else if (token == "GEOM") {
        token = param->getString(iarg++);
        if (token == "SPHERE") {
          double xs[3]; FOR_I3 xs[i] = param->getDouble(iarg++);
          double r = param->getDouble(iarg++);
          cout << " > SET_ST_FLAG: selecting sts inside sphere at x: " << COUT_VEC(xs) << " r: " << r << endl;
          int count = 0;
          for (int ist = 0; ist < ss.nst; ++ist) {
            if (ist%100000 == 0) cout << " > done " << ist << " out of " << ss.nst << endl;
            if (MiscUtils::getPointToTriDist2(xs,ss.xsp[ss.spost[ist][0]],ss.xsp[ss.spost[ist][1]],ss.xsp[ss.spost[ist][2]]) <= r*r) {
              ss.st_flag[ist] = 1;
              ++count;
            }
          }
          cout << " > count: " << count << endl;
        }
        else {
          cout << "ERROR: SET_ST_FLAG: unrecognized GEOM: \"" << token << "\". Check syntax." << endl;
          return;
        }
      }
      else {
        cout << "ERROR: SET_ST_FLAG: unrecognized token \"" << token << "\". Check syntax." << endl;
        return;
      }
    }

  }

  void processWriteStitch(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"writes a Stitch.in template, creating plce-holders for zone-based refinement windows\n" <<
          "  examples:\n" <<
          "    WRITE_STITCH_INPUT NAME=stitch.case_a.in \n"
          );
      return;
    }

    int iarg = 0;
    string filename = "stitch.template.in";
    bool bWithHelp = false;

    while (iarg < param->size()) {
      const string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "NAME") {
        filename = param->getString(iarg++);
      }
      else if (token == "LEAN") {
        bWithHelp = false;
      }
      // else if (token == "GUIDANCE") {
      //     const string guidance = param->getString(iarg++);
      //     if (guidance == "CURVATURE") {
      //       // suggest levels by curvature...
      //       // perhaps put levels in zone_flag?
      //     }
      //   }
    }
    ss.writeStitchInputFile(filename,bWithHelp);
  }

  void processShowStitch(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"send a templated stitch.in to the Cascade App; only meant for use from the App\n"
          );
      return;
    }

    int iarg = 0;
    bool bWithHelp = false;

    while (iarg < param->size()) {
      const string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "LEAN") {
        bWithHelp = false;
      }
      // else if (token == "GUIDANCE") {
      //     const string guidance = param->getString(iarg++);
      //     if (guidance == "CURVATURE") {
      //       // suggest levels by curvature...
      //       // perhaps put levels in zone_flag?
      //     }
      //   }
    }
    // TODO: update this format eventually?...
    WebUI::webUIOutput.add_(UIMessage_(STITCH_IN,ss.buildStitchInputFile(bWithHelp)));
  }

  void processWriteCharles(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"writes a Charles.in template, creating plce-holders for boundary conditions\n" <<
          "  examples:\n" <<
          "    WRITE_CHARLES_INPUT NAME=charles.case_a.in \n"
          );
      return;
    }

    // default values
    string filename = "charles.template.in";
    int eos = -1;
    bool bWithHelp = true;  // verbosity and text guidance in template input file

    int iarg = 0;
    while (iarg < param->size()) {
      const string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "NAME") {
        filename = param->getString(iarg++);
      }
      else if (token == "EOS") {
        const string _eos = MiscUtils::toUpperCase(param->getString(iarg++));
        if (_eos == "IDEAL_GAS") {
          eos = 0;
        }
        else if (_eos == "PREMIXED") {
          eos = 1;
        }
        else if (_eos == "NONPREMIXED") {
          eos = 2;
        }
        else {
          CWARN("unrecognized equation of state: " << _eos);
        }
      }
      else if (token == "LEAN") {
        bWithHelp = false;
      }
    }
    ss.writeCharlesInputFile(filename,eos,bWithHelp);
  }

  void processShowCharles(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"send a templated charles.in to the Cascade App; only meant for use from the App\n"
          );
      return;
    }

    // default values
    int eos = -1;
    bool bWithHelp = true;  // verbosity and text guidance in template input file

    int iarg = 0;
    while (iarg < param->size()) {
      const string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "EOS") {
        const string _eos = MiscUtils::toUpperCase(param->getString(iarg++));
        if (_eos == "IDEAL_GAS") {
          eos = 0;
        }
        else if (_eos == "PREMIXED_FPV") {
          eos = 1;
        }
        else if (_eos == "NONPREMIXED_FPV") {
          eos = 2;
        }
        else {
          CWARN("unrecognized equation of state: " << _eos);
        }
      }
      else if (token == "LEAN") {
        bWithHelp = false;
      }
    }
    // TODO: update this format...
    WebUI::webUIOutput.add_(UIMessage_(CHARLES_IN,ss.buildCharlesInputFile(eos,bWithHelp)));
  }

  void processEdgeLengthPdf(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    UNUSED(param);

    if (help) {
      WUI(INFO,"writes a histogram of the surface tri edge lengths to: edge_length.dat\n"
          );
      return;
    }

    double * edge_length = new double[ss.nst*3];
    for (int ist = 0; ist < ss.nst; ++ist) {
      FOR_I3 {
        edge_length[ist*3+i] = DIST(ss.xsp[ss.spost[ist][i]],ss.xsp[ss.spost[ist][(i+1)%3]]);
      }
    }

    Histogram h(MPI_COMM_NULL); // serial?
    h.setLog(true);
    h.add(edge_length,ss.nst*3);
    delete[] edge_length;

    h.write("edge_length.dat");
  }

  void processRunDiagnostics(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpRunDiagnostics();
      return;
    }

    int iarg=0;
    double dist_tol = 1.0e-10;
    double angle_tol_degrees = 0.5;
    bool check_self_intersections = false;
    bool check_self_intersections_new = false;
    bool detailed = false;
    bool subzone_info = false;
    string format = "none";
    int n_samples = 0;
    int nbr_layers = 0;
    while (iarg < param->size()) {
      const string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "DIST_TOL") {
        dist_tol = param->getDouble(iarg++);
      }
      else if (token == "ANGLE_TOL") {
        angle_tol_degrees = param->getDouble(iarg++);
      }
      else if (token == "CHECK_SELF_INTERSECTIONS") {
        check_self_intersections = true;
      }
      else if (token == "CHECK_SELF_INTERSECTIONS_NEW") {
        check_self_intersections_new = true;
      }
      else if (token == "WITH_DETAILS") {
        detailed = true;
      }
      else if (token == "WITH_SUBZONE_INFO") {
        subzone_info = true;
      }
      else if (token == "OUTPUT_FORMAT") {
        format = param->getString(iarg++);
      }
      else if (token == "SEED_SAMPLES") {
        n_samples = param->getInt(iarg++);
      }
      else if (token == "NBR_LAYERS") {
        nbr_layers = param->getInt(iarg++);
      }
      else {
        WUI(WARN,"unrecognized token: " << token);
      }
    }

    if (check_self_intersections_new) {
      ss.checkSelfIntersectionsNew();
    }

    ss.showSurfaceDiagnostics(detailed,subzone_info);
    ss.setNonManifoldCheckProperties(dist_tol,angle_tol_degrees,check_self_intersections,format,n_samples,nbr_layers);
    ss.diagnoseManifold();
  }

  void helpRunDiagnostics() {
    WUI(INFO,"run surface diagnostics that evaluate the manifold-ness of the surface\n" <<
        "  examples:\n" <<
        "    RUN_DIAGNOSTICS\n" <<
        "    RUN_DIAGNOSTICS WITH_DETAILS\n" <<
        "    RUN_DIAGNOSTICS CHECK_SELF_INTERSECTIONS\n" <<
        "    RUN_DIAGNOSTICS CHECK_SELF_INTERSECTIONS_NEW\n" <<
        "    RUN_DIAGNOSTICS CHECK_SELF_INTERSECTIONS SEED_SAMPLES 200\n" <<
        "  more info: [$CWIKB:surfer_diagnostics]"
        );
  }

  void processFlipTrisRandom(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    UNUSED(param);
    if (help) {
      WUI(INFO,"a debugging command that flips random tri normals; you shouldn't need or use this"
          );
      return;
    }

    try {
      ss.flipTrisRandom();
      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem flipping random tris");
    }
  }

  void processDeleteTrisWithCollocatedNodes(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    UNUSED(param);
    if (help) {
      WUI(INFO,"surface repair utility that deletes redundant tris; no arguments required"
          );
      return;
    }

    try {
      ss.deleteTrisWithCollocatedNodes();
      WebUI::webUIOutput.ensureImage();
      WebUI::webUIOutput.ensureMenu();
    }
    catch (int e) {
      WUI(WARN,"failure while deleting collocated tris");
    }
  }

  void processMergeCollocatedNodes(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"merges surface points; when no arguments provided, merges identical points. You can also\n" <<
          "specify a distance tolerance, and all points closer than that distance will be merged."
          );
      return;
    }

    try {
      // allow the user to specify a length scale...
      if (param->size() > 0) {
        const double eps = param->getDouble();
        ss.mergeCollocatedNodes(eps);
      }
      else {
        ss.mergeCollocatedNodes();
      }
      WebUI::webUIOutput.ensureImage();
      WebUI::webUIOutput.ensureMenu();
    }
    catch (int e) {
      WUI(WARN,"failure while merging collocated nodes");
    }
  }

  void processEnsureTeost(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    UNUSED(param);
    if (help) {
      WUI(INFO,"debugging tool to clear and build the edge data structure Teost; you shouldn't need or use this"
          );
      return;
    }

    try {
      ss.clearTeost();
      ss.ensureTeost();
    }
    catch (int e) {
      WUI(WARN,"failure while building teost (half-edge data structure)");
    }
  }

  void processAddOpenEdgeLoop(Param * param,const bool help) {

    if (help) {
      WUI(INFO,"debugging utility; you shouldn't need or use this"
          );
      return;
    }

    int ned_loop = 2;
    int iarg = 0;
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "NEDGES") {
        ned_loop = param->getInt(iarg++);
      }
    }
    if (!((ned_loop == 2)||(ned_loop == 3))) {
      CWARN("Unsupported number of edges in loop. Skipping");
      return;
    }

    // need teost
    ss.ensureTeost();
    if (ned_loop == 2) {
      // find one edge to open
      for (int ist = ss.nst; ist < ss.nst; ++ist) {
        bool found = false;
        FOR_I3 {
          if (!ss.isEdgeOpen(ist,i)) {
            int ist_nbr, i_nbr, orient_nbr;
            if (ss.getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
              found = true;
              ss.setTriNbrOpen(ist,i);
              ss.setTriNbrOpen(ist_nbr,i_nbr);
            }
          }

          if (found) break;
        }

        if (found) break;
      }
    }
    else if (ned_loop == 3) {
      // find a single tri to open edges on
      for (int ist = 0; ist < ss.nst; ++ist) {
        bool found = true;
        FOR_I3 if (ss.isEdgeOpen(ist,i)) found = false;
        if (!found) continue;
        FOR_I3 {
          if (!ss.isEdgeOpen(ist,i)) {
            int ist_nbr, i_nbr, orient_nbr;
            if (ss.getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
              ss.setTriNbrOpen(ist,i);
              ss.setTriNbrOpen(ist_nbr,i_nbr);
            }
          }
        }
        if (found) break;
      }
    }
    ss.clearOpenEdgeGroups();
  }

  void processEnsureOpenEdgeGroups(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    UNUSED(param);
    if (help) {
      WUI(INFO,"debugging utility; you shouldn't need or use this"
          );
      return;
    }

    // add a LOOP to a param for timing...
    /*
      int nloop = 1;
      int iarg=0;
      while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "LOOP") {
      nloop = param->getInt(iarg++);
      }
      else {
      cout << "Warning: processIndexOpenEdges: unrecognized token: " << token << endl;
      }
      }
    */

    ss.ensureOpenEdgeGroups();
  }

  // forces regrouping open edges using user-defined crease-angle
  void processRegroupOpenEdges(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpRegroupOpenEdges();
      return;
    }

    try {
      int iarg = 0;
      double crease_angle = 150.0;
      while (iarg < param->size()) {
        const string token = param->getString(iarg++);
        if (token == "CREASE_ANGLE") {
          crease_angle = param->getDouble(iarg++);
        }
        else {
          cout << "Warning: processRegroupOpenEdges: unrecognized token: " << token << endl;
        }
      }
      ss.clearOpenEdgeGroups();
      ss.ensureOpenEdgeGroups(crease_angle);

      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem while re-grouping open edges");
      helpRegroupOpenEdges();
    }
  }
  void helpRegroupOpenEdges() {
    WUI(INFO,"allows decimating open edge loops into sub-groups based on crease angle\n" <<
        "  examples:\n" <<
        "    REGROUP_OPEN_EDGES CREASE_ANGLE 165"
        );
  }

  void processCalcGcl(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"diagnostic that reports the surface's gcl"
          );
      return;
    }

    bool include_open_edge_groups = true;
    int iarg=0;
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "INCLUDE_OPEN_EDGE_GROUPS") {
        include_open_edge_groups = param->getBool(iarg++);
      }
    }
    double gcl[3];
    ss.calcGcl(gcl,include_open_edge_groups);
    WUI(INFO,"GCL: " << COUT_VEC(gcl));
  }

  void processCloseOpenEdgeLoops(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpCloseOpenEdgeLoops();
      return;
    }

    try {
      ss.ensureOpenEdgeGroups();
      const int ngr = ss.getNOpenEdgeGroups();

      if (ngr == 0) {
        cout << "There are no open edge groups. Returning." << endl;
        return;
      }
      //cout << ngr << endl;

      vector<int> edge_indices;

      bool ear_clipping = false; // default to mean_visible
      bool nested = false;
      double nested_tol = 0.001;
      bool donut = false;
      bool b_plane_fmm = false;
      bool b_parsed = false;

      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (!b_parsed && parseSpecifiedOpenEdges(edge_indices,token,param,iarg)) {
          // prune removes invalid entries AND gets loops in 0-indexed values
          b_parsed = true;
        }
        else if (token == "EAR_CLIPPING") {
          // if requested, use ear-clipping algorithm directly
          ear_clipping = true;
        }
        else if (token == "DONUT") {
          donut = true;
        }
        else if (token == "PLANAR_FMM") {
          b_plane_fmm = true;
        }
        else if (token == "NESTED") {
          nested = true;
        }
        else if (token == "NESTED_TOL") {
          nested_tol = param->getDouble(iarg++);
        }
        else {
          cout << "Warning: unrecognized CLOSE_OPEN_EDGE_LOOPS param: " << token << ". Skipping" << endl;
        }
      }
      if (ear_clipping) {
        // triangulate using ear clipping algorithm
        cout << "Closing open edge loops using ear clipping algorithm" << endl;
      }
      else if (donut) {
        cout << "Closing open edge loops using donut delauney algorithm." << endl;
      }
      else {
        // triangulate by triangulating agains the loop centroid
        cout << "Closing open edge loops using algorithm which assumes the loop center is visible to all vertices." << endl;
      }

      if (edge_indices.empty()) {
        // default is to work on all
        edge_indices.reserve(ss.getNOpenEdgeGroups());
        for (int iedge=0,lim=ss.getNOpenEdgeGroups(); iedge<lim; ++iedge) edge_indices.push_back(iedge);
      }

      if (b_plane_fmm) {
        // bool true at end indicates closing of open edges vs. eoi
        if (nested) ss.closePlanarHalfEdgeLoops(edge_indices,ss.openEdgeGroupDataVec,nested,nested_tol,"Closed_Nested_Open_Loops",true);
        else ss.closePlanarHalfEdgeLoops(edge_indices,ss.openEdgeGroupDataVec,nested,nested_tol,"",true);
      }
      else if (donut) {
        vector<pair<int,int> > donut_loops;
        if (nested) {
          ss.findNestedHalfEdgeLoops(donut_loops,nested_tol,edge_indices,ss.openEdgeGroupDataVec);
          if (donut_loops.empty()) {
            WUI(WARN,"DONUT loops were not detected in this selection; skipping");
            return;
          }
        }
        else {
          // for now only process if 2 passed in
          if (edge_indices.size() != 2) {
            WUI(WARN,"DONUT closing currently only supports passing 2 loops; skipping");
            return;
          }
          donut_loops.push_back(pair<int,int> (edge_indices[0],edge_indices[1]));
        }

        ss.closeHalfEdgeDonutLoops(donut_loops,true);
      }
      else {
        ss.closeHalfEdgeLoops(edge_indices,ear_clipping,true);
      }
      WebUI::webUIOutput.ensureImage();
      WebUI::webUIOutput.ensureMenu();
    }
    catch (int e) {
      WUI(WARN,"problem closing open edges");
      helpCloseOpenEdgeLoops();
    }
  }
  void helpCloseOpenEdgeLoops() {
    WUI(INFO,"allows user to try and fill/cap a specified set of open edge loops, i.e., holes\n" <<
        "  examples:\n" <<
        "    CLOSE_OPEN_EDGE_LOOPS (tries to fill all assuming simple holes)\n" <<
        "    CLOSE_OPEN_EDGE_LOOPS DONUT (close annular & concentric loops)\n" <<
        "  more info: [$CWIKB:surfer_repair]"
        );
  }

  void processCloseEdgeLoops(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpCloseEdgeLoops();
      return;
    }

    if (ss.eoiGroupDataVec.empty()) {
      COUT1("there are no dynamic edge groups; skipping");
      return;
    }

    try {
      vector<int> edge_indices;

      bool ear_clipping = false; // default to mean_visible
      bool b_donut = false;
      bool b_nested = false;
      double nested_tol = 0.001;
      bool b_planar_fmm = false;
      bool b_parsed = false;
      int iarg = 0;
      while (iarg < param->size()) {
        const string token = param->getString(iarg++);
        if (!b_parsed && parseSpecifiedDynamicEdges(edge_indices,token,param,iarg)) {
          // prune removes invalid entries AND gets loops in 0-indexed values
          b_parsed = true;
        }
        else if (token == "EAR_CLIPPING") {
          // if requested, use ear-clipping algorithm directly
          ear_clipping = true;
        }
        else if (token == "PLANAR_FMM") {
          b_planar_fmm = true;
        }
        else if (token == "NESTED") {
          b_nested = true;
        }
        else if (token == "NESTED_TOL") {
          nested_tol = param->getDouble(iarg++);
        }
        else if ((token == "DONUT")||(token == "ANNULUS")||(token == "RING")||(token == "CONCENTRIC")) {
          b_donut = true;
        }
        else {
          cout << "Warning: unrecognized CLOSE_EDGE_LOOPS param: " << token << ". Skipping" << endl;
        }
      }

      if (b_planar_fmm) {
        if (b_nested) ss.closePlanarHalfEdgeLoops(edge_indices,ss.eoiGroupDataVec,b_nested,nested_tol,"Closed_Nested_Loops");
        else ss.closePlanarHalfEdgeLoops(edge_indices,ss.eoiGroupDataVec,b_nested);
      }
      else if (b_donut) {
        vector<pair<int,int> > donut_loops;
        if (b_nested) {
          ss.findNestedHalfEdgeLoops(donut_loops,nested_tol,edge_indices,ss.eoiGroupDataVec);
          if (donut_loops.empty()) {
            WUI(WARN,"DONUT loops were not detected in this selection; skipping");
            return;
          }
        }
        else {
          // for now only process if 2 passed in
          if (edge_indices.size() != 2) {
            WUI(WARN,"DONUT closing currently only supports passing 2 loops; skipping");
            return;
          }
          donut_loops.push_back(pair<int,int> (edge_indices[0],edge_indices[1]));
        }

        ss.closeHalfEdgeDonutLoops(donut_loops,false);  // false = dynamic edges
      }
      else {
        ss.closeHalfEdgeLoops(edge_indices,ear_clipping,false);  // false = dynamic edges
      }
      WebUI::webUIOutput.ensureImage();
      WebUI::webUIOutput.ensureMenu();
    }
    catch (int e) {
      WUI(WARN,"problem closing edges");
      helpCloseEdgeLoops();
    }
  }
  void helpCloseEdgeLoops() {
    WUI(INFO,"allows user to try and fill/cap a specified set of dynamic edge loops,\n" <<
        "  examples:\n" <<
        "    CLOSE_EDGE_LOOPS\n" <<
        "    CLOSE_EDGE_LOOPS ZONES 41952,41953 DONUT\n" <<
        "    CLOSE_EDGE_LOOPS ZONES 41952,41953,41954,41955 NESTED DONUT\n" <<
        "  more info: [$CWIKB:surfer_repair]"
        );
  }

  void processZipOpenEdges(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpZipOpenEdges();
      return;
    }


    try {
      ss.ensureOpenEdgeGroups();
      const int ngr = ss.getNOpenEdgeGroups();

      if (ngr == 0) {
        cout << "There are no open edge groups. Returning." << endl;
        return;
      }

      vector<int> edge_indices;
      bool b_parsed = false;
      double edge_factor = 0.75;
      double delta_max = 1.0E+20;
      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));

        if (!b_parsed && parseSpecifiedOpenEdges(edge_indices,token,param,iarg)) {
          // prune removes invalid entries AND gets loops in 0-indexed values
          b_parsed = true;
        }
        else if (token == "EDGE_FACTOR") {
          edge_factor = param->getDouble(iarg++);
        }
        else if (token == "DELTA_MAX") {
          delta_max = param->getDouble(iarg++);
        }
        else {
          cout << "Warning: unrecognized ZIP_OPEN_EDGES param: " << token << ". Skipping" << endl;
        }
      }

      if (edge_indices.empty())
        ss.zipOpenEdges(edge_factor,delta_max);
      else
        ss.zipOpenEdges(edge_indices,edge_factor,delta_max);

      WebUI::webUIOutput.ensureImage();
      WebUI::webUIOutput.ensureMenu();
    }
    catch (int e) {
      WUI(WARN,"problem during zip open edges");
      helpZipOpenEdges();
    }
  }
  void helpZipOpenEdges() {
    WUI(INFO,"allows user to glue nearly adjacent open-edge loops together. By default this routine processes all open edges\n" <<
        "  examples:\n" <<
        "    ZIP_OPEN_EDGES\n" <<
        "    ZIP_OPEN_EDGES LOOPS 32679,32670,32671\n" <<
        "  more info: [$CWIKB:surfer_repair]"
        );
  }

  void processZipOpenEdgesCht(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }
    if (help) {
      WUI(INFO,"allows user to glue nearly adjacent open-edge loops together, and is guaranteed NOT to \n" <<
	  "move the nodes of the first zone listed. Examples:\n" <<
	  "    ZIP_OPEN_EDGES_CHT ZONES x0 y0\n" <<
	  "  more info: [$CWIKB:surfer_repair]"
	  );
      return;
    }

    try {

      bool got_zones = false;
      vector<string> zone0NameVec,zone1NameVec;

      bool got_zone_ids = false;
      vector<int> zone0idVec,zone1idVec;

      double crease_angle = 120.0;

      int iarg = 0;
      while (iarg < param->size()) {
	const string token = MiscUtils::toUpperCase(param->getString(iarg++));
	if ((token == "ZONE")||(token == "ZONES")) {
	  // ZONES is followed by 2 csv params, e.g. ZONES x0,t0 x1,t1
	  if (iarg > param->size()-2) {
	    WUI(WARN,"ZIP_OPEN_EDGES_CHT ZONES: expecting two sets of comma-separated zone names; skipping");
	    return;
	  }
	  MiscUtils::splitCsv(zone0NameVec,param->getString(iarg++));
	  MiscUtils::splitCsv(zone1NameVec,param->getString(iarg++));
	  got_zones = true;
	}
	else if (token == "ZONE_IDS") {
	  if (iarg > param->size()-2) {
	    WUI(WARN,"ZIP_OPEN_EDGES_CHT ZONE_IDS: expecting two sets of comma-separated zone ids; skipping");
	    return;
	  }
	  MiscUtils::splitCsv(zone0idVec,param->getString(iarg++));
	  MiscUtils::splitCsv(zone1idVec,param->getString(iarg++));
	  got_zone_ids = true;
	}
	else if (token == "CREASE_ANGLE") {
	  crease_angle = param->getDouble(iarg++);
	}
	else {
	  WUI(WARN,"ZIP_OPEN_EDGES_CHT: unrecognized token \"" << token << "\". Check syntax.");
	  return;
	}
      }

      if (!got_zones && !got_zone_ids) {
	WUI(WARN,"ZIP_OPEN_EDGES_CHT: missing ZONES or ZONE_IDS. Check syntax.");
	return;
      }

      ss.zone_flag.setLength(ss.zoneVec.size());
      ss.zone_flag.setAll(0);
      int ierr = 0;
      if (got_zones) {
        for (int ii = 0, ii_end=zone0NameVec.size(); ii < ii_end; ++ii) {
          const int izone = ss.getZoneIndex(zone0NameVec[ii]);
          if (izone == -1) {
            WUI(WARN,"ZIP_OPEN_EDGES_CHT: cannot find zone named \"" << zone0NameVec[ii] << "\"");
	    ierr = -1;
	    continue;
          }
	  assert((izone >= 0)&&(izone < int(ss.zoneVec.size())));
          ss.zone_flag[izone] = 1;
        }
        for (int ii = 0, ii_end=zone1NameVec.size(); ii < ii_end; ++ii) {
          const int izone = ss.getZoneIndex(zone1NameVec[ii]);
          if (izone == -1) {
            WUI(WARN,"ZIP_OPEN_EDGES_CHT: cannot find zone named \"" << zone1NameVec[ii] << "\"");
            ierr = -1;
            continue;
          }
          assert((izone >= 0)&&(izone < int(ss.zoneVec.size())));
          ss.zone_flag[izone] = 2;
        }
      }
      else if (got_zone_ids) {
        for (int ii = 0, ii_end=zone0idVec.size(); ii < ii_end; ++ii) {
          const int izone = zone0idVec[ii];
          if ((izone < 0)||(izone >= int(ss.zoneVec.size()))) {
            WUI(WARN,"ZIP_OPEN_EDGES_CHT: zone index " << izone << " out of range.");
	    ierr = -1;
	    continue;
	  }
          ss.zone_flag[izone] = 1;
        }
        for (int ii = 0, ii_end=zone1idVec.size(); ii < ii_end; ++ii) {
          const int izone = zone1idVec[ii];
          if ((izone < 0)||(izone >= int(ss.zoneVec.size()))) {
            WUI(WARN,"ZIP_OPEN_EDGES_CHT: zone index " << izone << " out of range.");
	    ierr = -1;
	    continue;
	  }
          ss.zone_flag[izone] = 2;
        }
      }

      if (ierr != 0)
        return;

      ierr = ss.zipOpenEdgesChtFlaggedZones(crease_angle);

      WebUI::webUIOutput.ensureImage();
      WebUI::webUIOutput.ensureMenu();
    }
    catch (int e) {
      WUI(WARN,"problem during zip open edges");
    }
  }

  void processMergeOpenEdgeNodes(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpMergeOpenEdgeNodes();
      return;
    }

    try {
      ss.ensureOpenEdgeGroups();
      const int ngr = ss.getNOpenEdgeGroups();

      if (ngr == 0) {
        cout << "There are no open edge groups. Returning." << endl;
        return;
      }

      vector<int> loop_indices;
      vector<int> zone_indices;
      double merge_tol = -1.0;
      bool all_oe_nodes = true;

      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (token == "EDGE_LOOPS") {
          all_oe_nodes = false;
          const string splitZonesCsv = param->getString(iarg++);
          vector<string> splitZonesVec;
          MiscUtils::splitCsv(splitZonesVec,splitZonesCsv);
          for (int i = 0, limit = splitZonesVec.size(); i < limit; ++i) {
            const int loop_index = (atoi(splitZonesVec[i].c_str())-OPEN_E_SZ_MIN);
            if (loop_index < 0 || loop_index >= ngr) {
              CWARN("loop index out-of-bounds: " << loop_index);
            }
            else {
              loop_indices.push_back(loop_index);
            }
          }
        }
        else if (token == "SUBZONES") {
          all_oe_nodes = false;
          const string splitZonesCsv = param->getString(iarg++);
          vector<string> splitZonesVec;
          MiscUtils::splitCsv(splitZonesVec,splitZonesCsv);
          for (int i = 0, limit = splitZonesVec.size(); i < limit; ++i) {
            const int zone_index = atoi(splitZonesVec[i].c_str());
            if (zone_index < 0 || zone_index > ss.nsz) {
              CWARN("subzone index out-of-bounds: " << zone_index);
            }
            else {
              zone_indices.push_back(zone_index);
            }
          }
        }
        else if (token == "TOL") {
          merge_tol = param->getDouble(iarg++);
        }
        else {
          cout << "Warning: unrecognized MERGE_OPEN_EDGE_NODES param: " << token << ". Skipping" << endl;
        }
      }

      if (!zone_indices.empty()) {
        if (!loop_indices.empty()) {
          CWARN("cannot merge both edge loops and zone loops at the same time; skipping");
          return;
        }
        if (zone_indices.size() != 2) {
          CWARN("must select exactly two zones to merge; skipping");
          return;
        }
        //HACK hard code to first two zones selected
        ss.mergeNodesOfOpenEdgeGroupsOnSubzones(zone_indices[0],zone_indices[1],merge_tol);
      }
      else if (!loop_indices.empty()) {
        if (!zone_indices.empty()) {
          CWARN("cannot merge both edge loops and zone loops at the same time; skipping");
          return;
        }
        if (loop_indices.size() != 2) {
          CWARN("must select exactly two loops to merge; skipping");
          return;
        }
        ss.mergeNodesOfOpenEdgeGroups(loop_indices[0],loop_indices[1],merge_tol);
      }
      else if (all_oe_nodes) {
        // merge all open edge nodes based on tolerance
        ss.mergeNodesOfOpenEdges(merge_tol);
      }
      else {
        CWARN("Unrecognized description of which nodes to merge; skipping");
      }
      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"failure during merge of open edge nodes");
      helpMergeOpenEdgeNodes();
    }
  }
  void helpMergeOpenEdgeNodes() {
    WUI(INFO,"merges nearly collocated (within a user specified tolerance) nodes on open-edge loops\n" <<
        "  examples:\n" <<
        "    MERGE_OPEN_EDGE_NODES TOL 0.02\n" <<
        "  more info: [$CWIKB:surfer_repair]"
        );
  }

  void processSetFeatureAngle(Param * param,const bool help) {
    if (help) {
      WUI(INFO,"sets feature crease angle that is used when grouping dynamic edges by FEATURE; takes no parameters\n"
          );
      return;
    }

    try {
      ss.setFeatureAngle(param->getDouble(0));
    }
    catch (int e) {
      WUI(WARN,"somehow setting the feature angle failed");
    }
  }

  void processGroupDynamicEdges(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpGroupDynamicEdges();
      return;
    }

    try {
      // ensure dynamic edges are grouped as desired
      SimpleSurface::DYNAMIC_EDGE_TYPE eType;
      double feat_angle = -1.0;
      vector<int> subzone_indices;
      int iarg = 0;
      while (iarg < param->size()) {
        const string token = param->getString(iarg++);
        if (token == "NONE") eType = SimpleSurface::NONE;
        else if (token == "MULTI") eType = SimpleSurface::MULTI;
        // else if (token == "METAL_BOUNDARY") eType = SimpleSurface::METAL_BOUNDARY;
        // else if (token == "FLUID_BOUNDARY") eType = SimpleSurface::FLUID_BOUNDARY;
        else if (token == "ZONE_BOUNDARY") eType = SimpleSurface::ZONE_BOUNDARY;
        else if (token == "SUBZONE_BOUNDARY") eType = SimpleSurface::SUBZONE_BOUNDARY;
        else if (token == "FEATURE") eType = SimpleSurface::FEATURE;
        else if (token == "OPEN") eType = SimpleSurface::OPEN;
        else if (token == "SELECTED_BOUNDARY") eType = SimpleSurface::SELECTED_BOUNDARY;
        else if (token == "NME") eType = SimpleSurface::NME;
        else if (token == "FEATURE_ANGLE") feat_angle = param->getDouble(iarg++);
        else if (parseSpecifiedSubzones(subzone_indices,token,param,iarg)) {
          // tris parsed, so add to selected tris
          ss.selectedSubzoneVec.clear();
          ss.selectedSubzoneVec = subzone_indices;
        }
        else {
          CWARN("unrecognized dynamic edge type: " << token << "; skipping");
          return;
        }
      }
      if (feat_angle > 0.0) ss.setFeatureAngle(feat_angle);

      ss.ensureDynamicEdgeGroups(eType);

      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"failure while re-grouping dynamic edges");
      helpGroupDynamicEdges();
    }
  }
  void helpGroupDynamicEdges() {
    WUI(INFO,"modifies how dynamic edges are registered, i.e., by what criteria\n" <<
        "  examples:\n" <<
        "    GROUP_DYNAMIC_EDGES MULTI\n" <<
        "    GROUP_DYNAMIC_EDGES ZONE_BOUNDARY\n" <<
        "    GROUP_DYNAMIC_EDGES SUBZONE_BOUNDARY\n" <<
        "    GROUP_DYNAMIC_EDGES FEATURE\n" <<
        "    GROUP_DYNAMIC_EDGES OPEN\n" <<
        "  more info: [$CWIKB:surfer_edges]"
        );
  }

  void processDiscretizeEdges(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpDiscretizeEdges();
      return;
    }

    try {
      vector<int> edge_indices;
      bool b_parsed = false;
      double delta = 0.0;
      int type = 0;  // later can be used for non-uniform spacing if desired...?
      bool b_keep_edge_nodes = false;

      int iarg = 0;
      while (iarg < param->size()) {
        string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (!b_parsed && parseSpecifiedDynamicEdges(edge_indices,token,param,iarg)) {
          // prunes invalid indices AND gets loops in 0-indexed values
          b_parsed = true;
        }
        else if (token == "DELTA") {
          delta = param->getDouble(iarg++);
        }
        else if (token == "TYPE") {
          type = param->getInt(iarg++);
        }
        else if (token == "KEEP_EDGE_NODES") {
          b_keep_edge_nodes = true;
        }
        else {
          WUI(WARN,"unrecognized token " << token << "; skipping");
        }
      }
      ss.reDiscretizeEdges(edge_indices,delta,type,b_keep_edge_nodes);
      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem rediscretizing edges");
      helpDiscretizeEdges();
    }
  }
  void helpDiscretizeEdges() {
    WUI(INFO,"re discretize edges of a dynamic edge loop\n" <<
        "where the edges are discretized uniformly into an integer number of new edges;\n" <<
        "thus actual length scale will be les than or equal to specified one.\n" <<
        "  examples:\n" <<
        "    DISCRETIZE_EDGES ZONES 41952 DELTA 0.05\n" <<
        "  more info: [$CWIKB:surfer_cht]"
        );
  }

  void processEnsureDelaunayOnZones(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpEnsureDelaunayOnZones();
      return;
    }

    try {
      vector<int> subzone_indices;
      int iters = 10;

      int iarg = 0;
      while (iarg < param->size()) {
        string token = param->getString(iarg++);
        if (parseSpecifiedSubzones(subzone_indices,token,param,iarg)) {
          // prune non-dynamic edge indices
          vector<int>::iterator it=subzone_indices.begin();
          while (it!=subzone_indices.end()) {
            if ((*it > SURFACE_SZ_MAX)) {
              COUT1(" > pruned specified zone "<< *it<< " because it is not a valid zone");
              it = subzone_indices.erase(it);
            }
            else {
              ++it;
            }
          }
        }
        else if (token == "ITERS" || token == "NITERS") {
          iters = param->getInt(iarg++);
        }
        else {
          WUI(WARN,"unrecognized token " << token << "; skipping");
        }
      }

      if (subzone_indices.empty()) {
        WUI(WARN,"must select subset of zones to test edge flipping on; skipping");
        throw(1);
      }
      else {
        ss.flagTrisFromSubzoneVec(subzone_indices);  // flag tris appropriately
      }
      ss.makeSelectedTrisDelaunay(iters);
      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem encountered during deaulany edge flipping");
      helpEnsureDelaunayOnZones();
    }
  }
  void helpEnsureDelaunayOnZones() {
    WUI(INFO,"process triangles on selected surface regions and perform edge flipping to make triangulation Delaunay\n" <<
        "  examples:\n" <<
        "    ENSURE_DELAUNAY ZONE_NAME z0 ITERS 20\n" <<
        "  more info: [$CWIKB:surfer_cht]"
        );
  }

  void processDiscretizeZones(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpDiscretizeZones();
      return;
    }

    try {
      vector<int> zone_indices;
      double delta = 0.0;
      double edge_factor = 0.8;
      double seed_factor = 0.9;
      int type = 0;  // later can be used for non-uniform spacing if desired...?
      int n_lloyd = 25;
      bool b_delaunay = false;
      bool b_keep_edge_nodes = false;
      bool b_power_diagram = false;
      double growth_factor = 1.0;
      double growth_power = 2.0;
      int init_type = 0;

      int iarg = 0;
      while (iarg < param->size()) {
        string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (parseSpecifiedZones(zone_indices,token,param,iarg)) {
          // prune non-dynamic edge indices
          vector<int>::iterator it=zone_indices.begin();
          while (it!=zone_indices.end()) {
            if ((*it > SURFACE_SZ_MAX)) {
              COUT1(" > pruned specified zone "<< *it<< " because it is not a valid zone");
              it = zone_indices.erase(it);
            }
            else {
              ++it;
            }
          }
        }
        else if (token == "DELTA") {
          delta = param->getDouble(iarg++);
        }
        else if (token == "EDGE_FACTOR") {
          edge_factor = param->getDouble(iarg++);
        }
        else if (token == "SEED_FACTOR") {
          seed_factor = param->getDouble(iarg++);
        }
        else if (token == "ENSURE_DELAUNAY") {
          b_delaunay = true;
        }
        else if (token == "ITERS" || token == "NITERS") {
          n_lloyd = param->getInt(iarg++);
        }
        else if (token == "TYPE") {
          type = param->getInt(iarg++);
          // value of 1: ignore adj tri feature termini
        }
        else if (token == "INIT_TYPE") {
          init_type = param->getInt(iarg++);
          // 0: uniform random, 1: R2 quasi-random
        }
        else if (token == "KEEP_EDGE_NODES") {
          b_keep_edge_nodes = true;
        }
        else if (token == "POWER_DIAGRAM") {
          b_power_diagram = true;
        }
        else if (token == "GROWTH_FACTOR") {
          growth_factor = param->getDouble(iarg++);
        }
        else if (token == "GROWTH_POWER") {
          growth_power = param->getDouble(iarg++);
        }
        else {
          WUI(WARN,"unrecognized token " << token << "; skipping");
        }
      }

      if (delta <= 0.0) {
        WUI(WARN,"invalid DELTA specified; skipping");
        throw(1);
      }
      if (zone_indices.empty()) {
        WUI(WARN,"must select subset of zones to tessellate; skipping");
        throw(1);
      }

      if (b_power_diagram) {
        seed_factor *= sqrt(2.0*(growth_factor/(growth_power+1)));
      }

      ss.reDiscretizeZones(zone_indices,delta,edge_factor,seed_factor,n_lloyd,type,b_delaunay,b_keep_edge_nodes,init_type,b_power_diagram,growth_factor,growth_power);
      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem encountered during rediscretization of zones");
      helpDiscretizeZones();
    }

  }
  void helpDiscretizeZones() {
    WUI(INFO,"discretize selected zones by a specific length-scale (retessellation)\n" <<
        "  examples:\n" <<
        "    DISCRETIZE_ZONES ZONE_NAME z0 DELTA 0.05 ITERS 20\n" <<
        "  more info: [$CWIKB:surfer_cht]"
        );
  }

  void processRoughenSelectedSubzones(Param* param,const bool b_help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (b_help) {
      helpRoughenSelectedSubzones();
      return;
    }

    try {
      vector<int> subzone_indices;
      double range[2];
      bool b_range = false;
      int n_filter = 10;

      int iarg = 0;
      while (iarg < param->size()) {
        string token = param->getString(iarg++);
        if (parseSpecifiedSubzones(subzone_indices,token,param,iarg)) {
          // prune subzones
          vector<int>::iterator it=subzone_indices.begin();
          while (it!=subzone_indices.end()) {
            if ((*it > SURFACE_SZ_MAX)) {
              COUT1(" > pruned specified zone "<< *it<< " because it is not a valid zone");
              it = subzone_indices.erase(it);
            }
            else {
              ++it;
            }
          }
        }
        else if (token == "RANGE") {
          range[0] = param->getDouble(iarg++); // min
          range[1] = param->getDouble(iarg++); // max
          if (range[1] < range[0]) {
            WUI(WARN,"MIN value smaller than MAX; swapping");
            std::swap(range[0],range[1]);
          }
          else if (range[1] == range[0]) {
            WUI(WARN,"MIN and MAX values cannot be identical");
            throw(1);
          }
          b_range = true;
        }
        else if (token == "NFILTER") {
          n_filter = param->getInt(iarg++);
        }
        else {
          WUI(WARN,"unrecognized token " << token << "; skipping");
        }
      }

      if (!b_range) {
        WUI(WARN,"must specify RANGE <min> <max> to roughen; skipping");
        throw(1);
      }
      else {
        ss.sz_flag.resize(ss.nsz);
        if (subzone_indices.empty()) {
          ss.sz_flag.setAll(1);
        }
        else {
          ss.sz_flag.setAll(0);
          for (vector<int>::iterator it=subzone_indices.begin(); it!=subzone_indices.end(); ++it)
            ss.sz_flag[*it] = 1;
        }
        ss.roughenFlaggedSubzones(range,n_filter);
        WebUI::webUIOutput.ensureImage();
      }
    }
    catch (int e) {
      WUI(WARN,"problem encountered during roughen selected subzones");
      helpRoughenSelectedSubzones();
    }

  }
  void helpRoughenSelectedSubzones() {
    WUI(INFO,"roughen selected zones by a uniform distribution\n" <<
        "  examples:\n" <<
        "    ROUGHEN ZONE_NAME z0 RANGE -0.05 0.05 # NFILTER 10 (default)\n" <<
        "    ROUGHEN ZONE_NAME z0 RANGE -0.05 0.05 NFILTER 0\n");
  }

  void processSetMaterial(Param * param) {}

  void processScaleOpenEdges(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpScaleOpenEdges();
      return;
    }

    try {
      bool b_norm = true;    // Whether the scaling factors are normalized or are distances
      double sx[3] = {0.0,0.0,0.0}; // Permits non-uniform scaling.  If zero, restrict scaling to non-zero directions only
      double x0[3] = {0.0,0.0,0.0};
      bool b_centroids = false;
      int n_duplicates = 0;
      vector<int> edge_indices;
      bool b_parsed = false;
      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (!b_parsed && parseSpecifiedOpenEdges(edge_indices,token,param,iarg)) {
          // prunes invalid indices AND gets loops in 0-indexed values
          b_parsed = true;
        }
        else if (token == "ABOUT") {
          FOR_I3 x0[i] = param->getDouble(iarg++);
        }
        else if (token == "ABOUT_CENTROIDS" || token == "ABOUT_CENTROID") {
          b_centroids = true;
        }
        else if (token == "DUPLICATES" || token == "DUPLICATE" || token == "EXTRUDE") {
          n_duplicates = 1;
        }
        else if (token == "N_DUPLICATES" || token == "N_DUPLICATE" || token == "N_EXTRUDE") {
          n_duplicates = param->getInt(iarg++);
        }
        else if (token == "DISTANCE") {
          b_norm = false;
          sx[0] = param->getDouble(iarg++);
          sx[1] = sx[0];
          sx[2] = sx[0];
        }
        else if (token == "FACTOR") {
          sx[0] = param->getDouble(iarg++);
          sx[1] = sx[0];
          sx[2] = sx[0];
        }
        else if (token == "SX") {
          FOR_I3 sx[i] = param->getDouble(iarg++);
        }
        else if (token == "X") {
          sx[0] = param->getDouble(iarg++);
        }
        else if (token == "Y") {
          sx[1] = param->getDouble(iarg++);
        }
        else if (token == "Z") {
          sx[2] = param->getDouble(iarg++);
        }
        else {
          WUI(WARN,"unrecognized SCALE_OPEN_EDGES param: " << token << "; skipping");
        }
      }

      ss.scaleOpenEdges(edge_indices, sx, x0,b_norm,b_centroids,n_duplicates);
      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem scaling open edges");
      helpScaleOpenEdges();
    }
  }

  void helpScaleOpenEdges() {
    WUI(INFO,"allows user to scale/resize (all) open-edge vertices directed about the vector defined by the origin\n" <<
        "  examples:\n" <<
        "      SCALE_OPEN_EDGES SX 10E-6 0.0 10E-6 ABOUT 0.0 0.0 0.0\n" <<
        "      SCALE_OPEN_EDGES SX 1.1   0.0 1.1   ABOUT 0.0 0.0 0.0 NORM\n" <<
        "      SCALE_OPEN_EDGES FACTOR 1.1 ABOUT 0.0 0.0 0.0\n" <<
        "  more info: [$CWIKB:surfer_repair]"
        );
  }

  void processTranslateOpenEdges(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpTranslateOpenEdges();
      return;
    }

    try {
      double dx[3] = {0.0,0.0,0.0};
      int n_duplicates = 0;
      vector<int> edge_indices;
      bool b_parsed = false;
      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (!b_parsed && parseSpecifiedOpenEdges(edge_indices,token,param,iarg)) {
          // prunes invalid indices AND gets loops in 0-indexed values
          b_parsed = true;
        }
        else if (token == "DUPLICATES" || token == "DUPLICATE" || token == "EXTRUDE") {
          n_duplicates = 1;
        }
        else if (token == "N_DUPLICATES" || token == "N_DUPLICATE" || token == "N_EXTRUDE") {
          n_duplicates = param->getInt(iarg++);
        }
        else if (token == "DX") {
          FOR_I3 dx[i] = param->getDouble(iarg++);
        }
        else if (token == "X") {
          dx[0] = param->getDouble(iarg++);
        }
        else if (token == "Y") {
          dx[1] = param->getDouble(iarg++);
        }
        else if (token == "Z") {
          dx[2] = param->getDouble(iarg++);
        }
        else {
          WUI(WARN,"unrecognized TRANSLATE_OPEN_EDGES param: " << token << "; skipping");
        }
      }

      ss.translateOpenEdges(edge_indices,dx,n_duplicates);
      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem translating open edges");
      helpTranslateOpenEdges();
    }
  }
  void helpTranslateOpenEdges() {
    WUI(INFO,"allows user to translate open-edge vertices\n" <<
        "  examples:\n" <<
        "      TRANSLATE_OPEN_EDGES DX 1.1 0.0 1.1 EXTRUDE\n" <<
        "  more info: [$CWIKB:surfer_repair]"
        );
  }

  void processRotateOpenEdges(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpRotateOpenEdges();
      return;
    }

    try {
      double axis[3] = {0.0,0.0,0.0};
      double point[3] = {0.0,0.0,0.0};
      double angle = 0.0;
      int n_duplicates = 0;
      vector<int> edge_indices;
      bool b_parsed = false;
      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (!b_parsed && parseSpecifiedOpenEdges(edge_indices,token,param,iarg)) {
          // prunes invalid indices AND gets loops in 0-indexed values
          b_parsed = true;
        }
        else if (token == "DUPLICATES" || token == "DUPLICATE" || token == "EXTRUDE") {
          n_duplicates = 1;
        }
        else if (token == "N_DUPLICATES" || token == "N_DUPLICATE" || token == "N_EXTRUDE") {
          n_duplicates = param->getInt(iarg++);
        }
        else if (token == "AXIS") {
          FOR_I3 axis[i] = param->getDouble(iarg++);
        }
        else if (token == "X") {
          axis[0] = 1.0;
          axis[1] = 0.0;
          axis[2] = 0.0;
        }
        else if (token == "Y") {
          axis[0] = 0.0;
          axis[1] = 1.0;
          axis[2] = 0.0;
        }
        else if (token == "Z") {
          axis[0] = 0.0;
          axis[1] = 0.0;
          axis[2] = 1.0;
        }
        else if (token == "POINT") {
          FOR_I3 point[i] = param->getDouble(iarg++);
        }
        else if (token == "ANGLE") {
          angle = param->getDouble(iarg++);
        }
        else {
          WUI(WARN,"unrecognized ROTATE_OPEN_EDGES param: " << token << "; skipping");
        }
      }

      ss.rotateOpenEdges(edge_indices,axis,point,angle,n_duplicates);
      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem rotating open edges");
      helpRotateOpenEdges();
    }
  }
  void helpRotateOpenEdges() {
    WUI(INFO,"allows user to rotate open-edge vertices\n" <<
        "  examples:\n" <<
        "      ROTATE_OPEN_EDGES AXIS 0 0 1 POINT 0 0 0 EXTRUDE\n" <<
        "  more info: [$CWIKB:surfer_repair]"
        );
  }

  void processForceOrthogonalityOpenEdges(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpForceOrthogonalityOpenEdges();
      return;
    }

    try {
      double dir[3] = {0.0,0.0,0.0};
      vector<int> edge_indices;
      bool b_parsed = false;
      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (!b_parsed && parseSpecifiedOpenEdges(edge_indices,token,param,iarg)) {
          // prunes invalid indices AND gets loops in 0-indexed values
          b_parsed = true;
        }
        else if (token == "DIR") {
          FOR_I3 dir[i] = param->getDouble(iarg++);
        }
        else if (token == "X") {
          dir[0] = 1.0;
        }
        else if (token == "Y") {
          dir[1] = 1.0;
        }
        else if (token == "Z") {
          dir[2] = 1.0;
        }
        else {
          WUI(WARN,"unrecognized FORCE_ORTHOGONALITY_OPEN_EDGES param: " << token << "; skipping");
        }
      }

      ss.forceOrthogonalityOpenEdges(edge_indices,dir);
      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem enforcing orthogonality of open edges");
      helpForceOrthogonalityOpenEdges();
    }
  }
  void helpForceOrthogonalityOpenEdges() {
    WUI(INFO,"allows user to adjust nodes on open edge loop(s) such that they are orthognal to a specified direction\n" <<
        "  examples:\n" <<
        "      FORCE_ORTHOGONALITY_OPEN_EDGES DIR 1 0 0\n" <<
        "  more info: [$CWIKB:surfer_repair]"
        );
  }

  void processNsmooth(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"NSMOOTH applies a simple smoothing algorithm to a specified region\n" <<
          "  examples:\n" <<
          "      NSMOOTH 10 SUBZONES 1,4\n" <<
          "      NSMOOTH 15 X 0.1 0.2 2 R 0.5"
	  );
      return;
    }

    try {
      int iarg = 0;
      const int niter = param->getInt(iarg++);
      vector<int> subzone_indices;
      double relax = 1.0;
      bool b_poisson = false;
      double x[3];
      bool b_x;
      double r;
      bool b_r;
      while (iarg < param->size()) {
        string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (parseSpecifiedSubzones(subzone_indices,token,param,iarg)) {
          //no op; parsed info
        }
        else if (token == "RELAX") {
          relax = param->getDouble(iarg++);
        }
        else if (token == "POISSON") {
          b_poisson = true;
        }
	else if (token == "X") {
	  FOR_I3 x[i] = param->getDouble(iarg++);
	  b_x = true;
	}
	else if (token == "R") {
	  r = param->getDouble(iarg++);
	  b_r = true;
	}
        else {
          WUI(WARN,"unrecognized NSMOOTH token: " << token << "; skipping");
        }
      }

      if (b_x && b_r) {
	ss.nsmoothNearPoint(niter,x,r);
      }
      else {
	ss.sz_flag.resize(ss.nsz);
	if (subzone_indices.empty()) {
	  ss.sz_flag.setAll(1);
	}
	else {
	  ss.sz_flag.setAll(0);
	  for (vector<int>::iterator it=subzone_indices.begin(); it!=subzone_indices.end(); ++it)
	    ss.sz_flag[*it] = 1;
	}
	ss.nsmoothFlaggedSubzones(niter,relax,b_poisson);
      }
      WebUI::webUIOutput.ensureImage();
      WUI(INFO,"NSMOOTH applied successfully");
    }
    catch (int e) {
      WUI(WARN,"problem with NSMOOTH");
    }

  }

  void processDeleteIsolatedNodes(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"deletes nodes that have been imported but are not participating in any valid tri\n" <<
          "  more info: [$CWIKB:surfer_repair]"
          );
      return;
    }

    try {
      ss.deleteIsolatedNodes();
    }
    catch (int e) {
      WUI(WARN,"problem deleting isolated nodes");
    }
  }

  void processRepairLinearTris(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpRepairLinearTris();
      return;
    }

    try {
      bool b_open = false;
      bool b_del = false;
      bool b_tri_split = false;

      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (token == "BY") {
          const string method = param->getString(iarg++);
          if (method == "DELETION") b_del = true;
          else if (method == "OPENING") b_open = true;
          else if (method == "TRI_SPLIT") b_tri_split = true;
          else CWARN("invalid option \"" << method << "\"; skipping");
        }
        else {
          CWARN("unrecognized REPAIR_LINEAR_TRIS string: " << token << ", skipping");
        }
      }

      if (b_open) ss.openLinearTris(ss.st_flag);
      else if (b_del) ss.deleteLinearTris();
      else if (b_tri_split) ss.splitLinearTriNeighbor();

      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem repairing linear tris");
      helpRepairLinearTris();
    }
  }
  void helpRepairLinearTris() {
    WUI(INFO,"allows user mutiple approaches to repairing \"linear\" tris\n" <<
        "  examples:\n" <<
        "    REPAIR_LINEAR_TRIS BY DELETION (removes offending tris)\n" <<
        "    REPAIR_LINEAR_TRIS BY OPENING (spreading nodes apart)\n" <<
        "    REPAIR_LINEAR_TRIS BY TRI_SPLIT (split opposing tri)\n" <<
        "  more info: [$CWIKB:surfer_repair]"
        );
  }

  void processRepairMultiNeighborTris(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"allows user to delete tris adjacent to multi-edges\n" <<
          "  more info: [$CWIKB:surfer_repair]"
          );
      return;
    }

    try {
      bool b_del = false;

      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (token == "BY") {
          const string method = param->getString(iarg++);
          if (method == "DELETION") b_del = true;
          else CWARN("invalid option \"" << method << "\"; skipping");
        }
        else {
          CWARN("unrecognized REPAIR_MULTI_NBR_TRIS string: " << token << ", skipping");
        }
      }

      if (b_del) ss.deleteMultiNeighborTris();

      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem repairing multi-edge neighbors");
    }
  }

  void processShowGeoDistFromPoint(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"computes geodesic distance at all tris from a seed location on the surface\n" <<
          "  examples:\n" <<
          "    SHOW_GEODESIC_DISTANCE FROM_POINT 1.2 0 0"
          );
      return;
    }

    // expects x y and z
    int iarg = 0;
    assert(MiscUtils::toUpperCase(param->getString(iarg++)) == "FROM_POINT");
    double x[3]; FOR_I3 x[i] = param->getDouble(iarg++);
    cout << "Performing Eikonal solve over surface from point " << COUT_VEC(x) << " using fast marching method." << endl;
    ss.calcGeoDistFromPoint(x);

    surface_data_type = 1;
    surface_data_name = "geodesic_D";

    WebUI::webUIOutput.ensureImage();
  }

  void processDeleteSelectedSubzones(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpDeleteSelectedSubzones();
      return;
    }

    try {
      // expects a comma delimited list of subzone indices
      vector<int> subzone_indices;

      int iarg = 0;

      ss.st_flag.resize(ss.nst);
      ss.st_flag.setAll(0);
      bool b_tri_ids = false;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (parseSpecifiedSubzones(subzone_indices,token,param,iarg)) {
          //no op; parsed info
        }
        else if (token == "TRI_IDS") {
          const string byCsv = param->getString(iarg++);
          vector<int> byVec;
          MiscUtils::splitCsv(byVec,byCsv);
          for (int ii = 0,lim = byVec.size(); ii < lim; ++ii) {
            if ((byVec[ii] >= 0)&&(byVec[ii] < ss.nst)) {
              ss.st_flag[byVec[ii]] = 1;
              b_tri_ids = true;
            }
          }
        }
        else {
          WUI(WARN,"unrecognized DELETE param: " << token << ", skipping");
        }
      }

      if (b_tri_ids) {
        ss.deleteFlaggedTris();
      }
      else if ( !subzone_indices.empty() ) {
        ss.deleteSelectedSubzones(subzone_indices);
        WebUI::webUIOutput.ensureMenu();  // in case zones changed
        WebUI::webUIOutput.ensureImage();
      }
      else {
        WUI(WARN,"no surface(s) was selected for deleting; skipping");
      }
    }
    catch (int e) {
      WUI(WARN,"problem deleting subzones");
    }
  }
  void helpDeleteSelectedSubzones() {
    WUI(INFO,"allows users to delete user-specified portions of the surface\n" <<
        "  examples:\n" <<
        "    DELETE ZONE_NAMES top,bottom\n" <<
        "    DELETE SUBZONES 0,3,2,4"
        );
  }

  void processTranslateSelectedSubzones(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpTranslateSelectedSubzones();
      return;
    }

    try {
      // expects a comma delimited list of subzone indices

      // container which persists selected/highlgihting in UI
      // this transform doesn't change subzones, just their metadata so can leverage to persist selection
      ss.selectedSubzoneVec.clear();
      double dx[3]       = {0.0,0.0,0.0};
      double dn          = 0.0;
      int translate_type = -1;
      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));

        if (parseSpecifiedSubzones(ss.selectedSubzoneVec,token,param,iarg)) {
          //no op; parsed info
        }
        else if (token == "DX") {
          FOR_I3 dx[i]   = param->getDouble(iarg++);
          translate_type = 0;
        }
        else if (token == "X") {
          dx[0] = param->getDouble(iarg++);
          translate_type = 0;
        }
        else if (token == "Y") {
          dx[1] = param->getDouble(iarg++);
          translate_type = 0;
        }
        else if (token == "Z") {
          dx[2] = param->getDouble(iarg++);
          translate_type = 0;
        }
        else if ( token == "DN") {
          dn = param->getDouble(iarg++);
          translate_type = 1;
        }
      }

      if ( translate_type == 0 ) {
        // this is a cartesian translation with the dx vector...
        if ( !ss.selectedSubzoneVec.empty() ) ss.translateSelectedSubzones(ss.selectedSubzoneVec,dx);
        else ss.translateSurface(dx);
      }
      else if ( translate_type == 1) {
        // normal translation of the surface ...
        if ( ss.selectedSubzoneVec.empty() ) {
          // default to entire surface
          for (int isz=0; isz<ss.nsz; ++isz) ss.selectedSubzoneVec.push_back(isz);
        }
        ss.translateSelectedSubzonesNormal(ss.selectedSubzoneVec,dn);
      }
      else {
        WUI(WARN," > unrecognized translation type: " << translate_type);
      }

      WebUI::webUIOutput.ensureImage();
      WUI(WARN,"Transform could have made your surface camera shy. Please Reset Scene.");
    }
    catch (int e) {
      WUI(WARN,"failure to translate selection");
      helpTranslateSelectedSubzones();
    }
  }
  void helpTranslateSelectedSubzones() {
    WUI(INFO,"translates specified portions of the surface by a cartesian offset\n" <<
        "  examples:\n" <<
        "    TRANSLATE ALL DX 1 0 0\n" <<
        "    TRANSLATE SUBZONES 4,5,6 DX 1 4 0\n" <<
        "  more info: [$CWIKB:surfer_transform]"
        );
  }

  void processSnapToCOM(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"translates surface so that (0,0,0) is at the center of mass based on surface nodes"
          );
      return;
    }

    try {
      // snap-to-center-of-mass
      double x_com[3];
      ss.calcCenterOfMass(x_com);
      FOR_I3 x_com[i] = -x_com[i];
      ss.translateSurface(x_com);
      ss.calcCenterOfMass(x_com);  // run again to confirm 0...

      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem snapping to center-of-mass");
    }
  }

  void processRotateSelectedSubzones(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpRotateSelectedSubzones();
      return;
    }

    try {
      // expects a comma delimited list of subzone indices
      // container which persists selected/highlgihting in UI
      // this transform doesn't change subzones, just their metadata so can leverage to persist selection
      ss.selectedSubzoneVec.clear();

      double axis[3] = {0.0,0.0,0.0};
      double point[3] = {0.0,0.0,0.0};
      double angle = 0.0;
      bool b_smooth = false;

      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));

        if (parseSpecifiedSubzones(ss.selectedSubzoneVec,token,param,iarg)) {
          //no op; parsed info
        }
        else if (token == "AXIS") {
          FOR_I3 axis[i] = param->getDouble(iarg++);
        }
        else if (token == "X") {
          axis[0] = 1.0;
          axis[1] = 0.0;
          axis[2] = 0.0;
        }
        else if (token == "Y") {
          axis[0] = 0.0;
          axis[1] = 1.0;
          axis[2] = 0.0;
        }
        else if (token == "Z") {
          axis[0] = 0.0;
          axis[1] = 0.0;
          axis[2] = 1.0;
        }
        else if (token == "POINT") {
          FOR_I3 point[i] = param->getDouble(iarg++);
        }
        else if (token == "ANGLE") {
          angle = param->getDouble(iarg++);
        }
        else if (token == "SMOOTH") {
          b_smooth = true;
        }
        else {
          CWARN("ROTATE: unrecognized token \"" << token << "\"");
        }
      }

      if (angle == 0.0) {
        // CWARN("a rotation angle was not specified; skipping");
        WUI(WARN,"a rotation angle was not specified; skipping this operation");
        return;
      }

      if (MAG(axis) == 0.0) {
        // CWARN("an invalid rotation axis was specified; skipping");
        WUI(WARN,"Whoops, you've specified an invalid rotation axis. Turn around (get it?!) and try again.");
        return;
      }

      if ( !ss.selectedSubzoneVec.empty() ) {
        if (b_smooth)
          ss.morphRotateSelectedSubzones(ss.selectedSubzoneVec,axis,point,angle);
        else
          ss.rotateSelectedSubzones(ss.selectedSubzoneVec,axis,point,angle);
      }
      else {
        ss.rotateSurface(axis,point,angle);
      }

      WebUI::webUIOutput.ensureImage();
      WUI(WARN,"Transform could have made your surface camera shy. Please Reset Scene.");
    }
    catch (int e) {
      WUI(WARN,"problem rotating selected surface");
      helpRotateSelectedSubzones();
    }
  }
  void helpRotateSelectedSubzones() {
    WUI(INFO,"rotates specified portions of the surface about an arbitrary axis\n" <<
        "  examples:\n" <<
        "    ROTATE ALL AXIS 1 0 0 POINT 13 4 5 ANGLE 30.0\n" <<
        "    ROTATE X ZONE_NAMES top,left,right ANGLE 30.0\n" <<
        "  more info: [$CWIKB:surfer_transform]"
        );
  }

  void processCopySelectedSubzones(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpCopySelectedSubzones();
      return;
    }

    try {
      // expects a comma delimited list of subzone indices
      // container which persists selected/highlgihting in UI
      // this transform doesn't change subzones, just their metadata so can leverage to persist selection
      ss.selectedSubzoneVec.clear();
      bool b_all = true;
      bool b_dx = false;
      double dx[3] = {0.0,0.0,0.0};  // use as axis or translation
      double point[3] = {0.0,0.0,0.0};
      double angle = 0.0;
      int transform_type = -1;
      int copies = 1;

      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (parseSpecifiedSubzones(ss.selectedSubzoneVec,token,param,iarg)) {
          b_all = false;
        }
        else if (token == "TRANSLATE") {
          transform_type = 0;
        }
        else if (token == "ROTATE") {
          transform_type = 1;
        }
        else if ((token == "AXIS")||(token == "DX")) {
          b_dx = true;
          FOR_I3 dx[i] = param->getDouble(iarg++);
        }
        else if (token == "X") {
          dx[0] = 1.0;
          dx[1] = 0.0;
          dx[2] = 0.0;
        }
        else if (token == "Y") {
          dx[0] = 0.0;
          dx[1] = 1.0;
          dx[2] = 0.0;
        }
        else if (token == "Z") {
          dx[0] = 0.0;
          dx[1] = 0.0;
          dx[2] = 1.0;
        }
        else if (token == "POINT") {
          FOR_I3 point[i] = param->getDouble(iarg++);
        }
        else if (token == "ANGLE") {
          angle = param->getDouble(iarg++);
        }
        else if (token == "COPIES") {
          copies = param->getInt(iarg++);
        }
        else if (token == "ALL") {
          b_all = true;
        }
        else {
          CWARN("COPY: unrecognized token \"" << token << "\"");
        }
      }

      // error checking...

      if (b_all) {
        ss.selectedSubzoneVec.clear();
        for (int isz=0; isz<ss.nsz; ++isz) ss.selectedSubzoneVec.push_back(isz);
      }
      if (copies < 1) {
        WUI(WARN,"Whoops, you've specified too few copies; skipping.");
        return;
      }
      if (transform_type == -1) {
        if (b_dx) transform_type = 0;
        else {
          WUI(WARN,"you must specify TRANSLATE or ROTATE; skipping.");
          return;
        }
      }

      if (transform_type == 0) {
        if ((!b_dx)||(MAG(dx) == 0.0)) {
          WUI(WARN,"DX was not specified or was zero; skipping.");
          return;
        }
        if ( !ss.selectedSubzoneVec.empty() ) ss.copyTranslateSelectedSubzones(ss.selectedSubzoneVec,copies,dx);
        else CWARN("no surface(s) was selected for translating; skipping");
      }
      else if (transform_type == 1) {
        if (angle == 0.0) {
          WUI(WARN,"a rotation angle was not specified; skipping this operation");
          return;
        }
        if (MAG(dx) == 0.0) {
          WUI(WARN,"Whoops, you've specified an invalid rotation axis. Turn around (get it?!) and try again.");
          return;
        }
        if ( !ss.selectedSubzoneVec.empty() ) ss.copyRotateSelectedSubzones(ss.selectedSubzoneVec,copies,dx,point,angle);
        else CWARN("no surface(s) was selected for rotating; skipping");
      }
      WebUI::webUIOutput.ensureImage();
      WebUI::webUIOutput.ensureMenu();
      WUI(WARN,"Transform could have made your surface camera shy. Please Reset Scene.");
    }
    catch (int e) {
      WUI(WARN,"problem copying selected surface");
      helpCopySelectedSubzones();
    }
  }
  void helpCopySelectedSubzones() {
    WUI(INFO,"copy specified portions of the surface and either translates or rotates\n" <<
        "  examples:\n" <<
        "    COPY ROTATE AXIS 1 0 0 POINT 13 4 5 ANGLE 30.0 ALL\n" <<
        "    COPY TRANSLATE ZONE_NAMES top,left,right DX 5 0 0 COPIES 3\n" <<
        "  more info: [$CWIKB:surfer_transform]"
        );
  }

  void processSplitZonesAfterCopyAll(Param * param,const bool help) {
    // like copy below, but acts on entire mesh always and adds new zone names with
    // index suffixes...
    try {
      int n = param->getInt();
      cout << " > got N " << n << endl;
      ss.splitZonesAfterCopyAll(n);
    }
    catch (int e) {
      WUI(WARN,"problem with splitZonesAfterCopyAll");
    }
  }

  void processScaleSelectedSubzones(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpScaleSelectedSubzones();
      return;
    }

    try {
      // expects a comma delimited list of subzone indices
      // container which persists selected/highlgihting in UI
      // this transform doesn't change subzones, just their metadata so can leverage to persist selection
      ss.selectedSubzoneVec.clear();
      //bool b_all = true;
      bool b_about = false;
      bool b_scale_centroid = false;
      bool b_radial = false;
      bool b_axis = false;
      bool b_norm = true;
      //int ierr = 0;
      double sx[3] = {1.0,1.0,1.0};
      double scale_centroid[3] = {0.0,0.0,0.0};
      double axis[3] = {0.0,0.0,0.0};
      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));

        if (parseSpecifiedSubzones(ss.selectedSubzoneVec,token,param,iarg)) {
          //no op; parsed info
        }
        else if (token == "ABOUT") {
          b_about = true;
          FOR_I3 scale_centroid[i] = param->getDouble(iarg++);
        }
        else if (token == "CENTROID") {
          b_scale_centroid = true;
        }
        else if (token == "RADIAL") {
          b_radial = true;
        }
        else if (token == "AXIS") {
          b_axis = true;
          FOR_I3 axis[i] = param->getDouble(iarg++);
        }
        else if (token == "FACTOR") {
          sx[0] = param->getDouble(iarg++);
          sx[1] = sx[0];
          sx[2] = sx[0];
        }
        else if (token == "DISTANCE") {
          b_norm = false;
          sx[0] = param->getDouble(iarg++);
          sx[1] = sx[0];
          sx[2] = sx[0];
        }
        else if (token == "SX") {
          FOR_I3 sx[i] = param->getDouble(iarg++);
        }
        else if (token == "X") {
          sx[0] = param->getDouble(iarg++);
        }
        else if (token == "Y") {
          sx[1] = param->getDouble(iarg++);
        }
        else if (token == "Z") {
          sx[2] = param->getDouble(iarg++);
        }
        else {
          WUI(WARN,"unrecognized SCALE token: " << token << "; skipping");
        }
      }

      if (b_scale_centroid && b_about) {
        WUI(WARN,"Cannot specify both \"ABOUT <x> <y> <z>\" and CENTROID; using the ABOUT reference point for scaling");
        b_scale_centroid = false;
      }

      if (b_radial && !b_axis) {
        WUI(WARN,"An AXIS direction must be specified in order to scale radially; syntax uses a point (via ABOUT or CENTRIOD) and AXIS; skipping");
        return;
      }

      if ((sx[0] == 1.0) && (sx[1] == 1.0) && (sx[2] == 1.0)) {
        // CWARN("scaling by 1 in all dimensions results in no changes to the surface; skipping");
        WUI(WARN,"So you've decided to scale the surface by 1.0 in all dimensions...which does nothing...so I'm just going to skip this command...");
        return;
      }

      if (!ss.selectedSubzoneVec.empty()) {
        if (b_radial) ss.scaleSelectedSubzonesRadially(ss.selectedSubzoneVec,sx[0],scale_centroid,axis,b_scale_centroid,b_norm);
        else ss.scaleSelectedSubzones(ss.selectedSubzoneVec,sx,scale_centroid,b_scale_centroid,b_norm);
      }
      else {
        if (b_radial) {
          ss.selectedSubzoneVec.reserve(ss.nsz);
          for (int isz=0; isz<ss.nsz; ++isz) ss.selectedSubzoneVec.push_back(isz);
          ss.scaleSelectedSubzonesRadially(ss.selectedSubzoneVec,sx[0],scale_centroid,axis,b_scale_centroid,b_norm);
        }
        else ss.scaleSurface(sx,scale_centroid,b_scale_centroid,b_norm);
      }

      WebUI::webUIOutput.ensureImage();
      WUI(WARN,"Transform could have made your surface camera shy. Please Reset Scene.");
    }
    catch (int e) {
      WUI(WARN,"problem scaling selected surface");
      helpScaleSelectedSubzones();
    }
  }
  void helpScaleSelectedSubzones() {
    WUI(INFO,"scales specified portions of the surface about an point\n" <<
      "  examples:\n" <<
      "    SCALE ALL SX 2 1.2 1.2 (control over x,y,z scaling)\n" <<
      "    SCALE ALL DISTANCE 0.04 (scale a specific distance from reference point)\n" <<
      "    SCALE ZONE_NAMES tube RADIAL ABOUT 0 0 0 AXIS 1 0 0 DISTANCE 0.04 (scales specified surface a distance of 0.04m about x-axis)\n" <<
      "    SCALE ZONE_NAMES top,left,right FACTOR 0.5 (single parameter defines scaling factor)\n" <<
      "  more info: [$CWIKB:surfer_transform]"
    );
  }

  void processMirrorSelectedSubzones(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpMirrorSelectedSubzones();
      return;
    }

    try {
      // expects a comma delimited list of subzone indices
      // container which persists selected/highlgihting in UI
      // this transform doesn't change subzones, just their metadata so can leverage to persist selection
      ss.selectedSubzoneVec.clear();

      double xp[3] = {0.0,0.0,0.0};
      double np[3] = {1.0,0.0,0.0};
      bool b_plane = false;

      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));

        if (parseSpecifiedSubzones(ss.selectedSubzoneVec,token,param,iarg)) {
          //no op; parsed info
        }
        else if (token == "PLANE") {
          b_plane = true;
          xp[0] = param->getDouble(iarg++);
          xp[1] = param->getDouble(iarg++);
          xp[2] = param->getDouble(iarg++);
          np[0] = param->getDouble(iarg++);
          np[1] = param->getDouble(iarg++);
          np[2] = param->getDouble(iarg++);
        }
        else if (token == "X") {
          b_plane = true;
        }
        else if (token == "Y") {
          b_plane = true;
          np[0] = 0.0;
          np[1] = 1.0;
          np[2] = 0.0;
        }
        else if (token == "Z") {
          b_plane = true;
          np[0] = 0.0;
          np[1] = 0.0;
          np[2] = 1.0;
        }
      }

      if (!b_plane) {
        CWARN("mirroring plane not specified; skipping");
        WUI(WARN,"mirroring plane wasn't properly specified; so I'm just... not going to do... anything... Yeah.");
        return;
      }

      if (ss.selectedSubzoneVec.empty() ) {
        for (int isz=0; isz<ss.nsz; ++isz) ss.selectedSubzoneVec.push_back(isz);
      }

      ss.mirrorSelectedSubzones(ss.selectedSubzoneVec,xp,np);
      ss.flipSelectedSubzones(ss.selectedSubzoneVec);  // mirroring causes surface normals to invert; so correct for this

      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem mirroring selected surface");
      helpMirrorSelectedSubzones();
    }
  }
  void helpMirrorSelectedSubzones() {
    WUI(INFO,"mirrors specified portions of the surface about an arbitrary plane\n" <<
        "  examples:\n" <<
        "    MIRROR ALL PLANE 0 0 0 1 0 0 (about yz-axis)\n" <<
        "    MIRROR ZONE_NAMES the_box PLANE 0 0.4 0 1 1 0 (about arbitrary plane)\n" <<
        "  more info: [$CWIKB:surfer_transform]"
        );
  }

  void processImprint(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpImprint();
      return;
    }

    try {
      vector<int> subzone_indices;
      bool split = false;

      bool b_plane = false;
      double xp[3] = {0.0,0.0,0.0};
      double np[3];

      double tol = 0.0;

      bool b_cyl = false;
      double rp;

      bool b_window = false;
      double window[4][3];

      int iarg = 0;
      while (iarg < param->size()) {
        string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (parseSpecifiedSubzones(subzone_indices,token,param,iarg)) {
          //no op; parsed info
        }
        else if (token == "PLANE") {
          b_plane = true;
          FOR_I3 xp[i] = param->getDouble(iarg++);
          FOR_I3 np[i] = param->getDouble(iarg++);
        }
        else if (token == "CYL") {
          b_cyl = true;
          FOR_I3 xp[i] = param->getDouble(iarg++);
          FOR_I3 np[i] = param->getDouble(iarg++);
        }
        else if (token == "CYL_X") {
          b_cyl = true;
          np[0] = 1.0; np[1] = 0.0; np[2] = 0.0;
        }
        else if (token == "CYL_Y") {
          b_cyl = true;
          np[0] = 0.0; np[1] = 1.0; np[2] = 0.0;
        }
        else if (token == "CYL_Z") {
          b_cyl = true;
          np[0] = 0.0; np[1] = 0.0; np[2] = 1.0;
        }
        else if (token == "POINT") {
          FOR_I3 xp[i] = param->getDouble(iarg++);
        }
        else if (token == "R") {
          rp = param->getDouble(iarg++);
        }
        else if (token == "WINDOW") {
          FOR_I3 window[0][i] = param->getDouble(iarg++);
          FOR_I3 window[1][i] = param->getDouble(iarg++);
          FOR_I3 window[2][i] = param->getDouble(iarg++);
          FOR_I3 window[3][i] = param->getDouble(iarg++);
          b_window = true;
        }
        else if (token == "SPLIT_ZONES") {
          split = true;
        }
        else if (token == "TOL") {
          tol = param->getDouble(iarg++);
        }
        else {
          CWARN("unrecognized IMPRINT param \"" << token << "\"; skipping");
        }
      }

      // flag subzones
      ss.sz_flag.resize(ss.nsz);
      if (subzone_indices.empty()) {
        ss.sz_flag.setAll(1);
      }
      else {
        ss.sz_flag.setAll(0);
        for (vector<int>::iterator it=subzone_indices.begin(); it!=subzone_indices.end(); ++it) ss.sz_flag[*it] = 1;
      }

      if (b_plane) {
        ss.imprintPlane(xp,np,tol,split);
      }
      else if (b_cyl) {
        ss.imprintCyl(xp,np,rp,split);
      }
      else if (b_window) {

        const double dx01[3] = DIFF(window[1],window[0]);
        const double dx12[3] = DIFF(window[2],window[1]);
        const double n[3] = CROSS_PRODUCT(dx01,dx12);

        const int nzn0 = ss.zoneVec.size();
        int *nsz_zn0 = new int[nzn0];
        for (int izn = 0; izn < nzn0; ++izn)
          nsz_zn0[izn] = ss.szozn_i[izn+1]-ss.szozn_i[izn];

        vector<int>* flagged_local_subzones0 = new vector<int>[nzn0];
        for (vector<int>::iterator it = subzone_indices.begin(); it != subzone_indices.end(); ++it) {
          const int isz = *it;
          int izn = -1;
          for (izn = 0; izn < nzn0; ++izn) {
            if (ss.szozn_i[izn+1] > isz)
              break;
          }
          assert((izn >= 0)&&(izn < nzn0));
          flagged_local_subzones0[izn].push_back(isz-ss.szozn_i[izn]);
        }

        int ii0 = 3;
        for (int ii = 0; ii < 4; ++ii) {

          // get plane orthogonal to this edge and n
          const double dx_ii[3] = DIFF(window[ii],window[ii0]);
          const double n_ii[3] = CROSS_PRODUCT(dx_ii,n);
          ss.imprintPlane(window[ii],n_ii,tol,split);

          if (ii < 3) {

            // reflag subzones. could do a better job if need be by maintaining sz_flag
            ss.sz_flag.resize(ss.nsz);
            if (subzone_indices.empty()) {
              ss.sz_flag.setAll(1);
            }
            else {
              ss.sz_flag.setAll(0);
              for (int izn = 0; izn < nzn0; ++izn) {
                // flag old subzones...
                for (int ii = 0, lim = flagged_local_subzones0[izn].size(); ii < lim; ++ii)
                  ss.sz_flag[ss.szozn_i[izn]+flagged_local_subzones0[izn][ii]] = 1;
                // flag new subzones...
                const int nsz = ss.szozn_i[izn+1]-ss.szozn_i[izn];
                for (int isz = nsz_zn0[izn]; isz < nsz; ++isz)
                  ss.sz_flag[ss.szozn_i[izn]+isz] = 1;
              }
              // flag new zones...
              const int nzn = ss.zoneVec.size();
              for (int izn = nzn0; izn < nzn; ++izn) {
                for (int isz = ss.szozn_i[izn]; isz != ss.szozn_i[izn+1]; ++isz)
                  ss.sz_flag[isz] = 1;
              }
            }

          }

          ii0 = ii;
        }
        delete[] nsz_zn0;
        delete[] flagged_local_subzones0;

      }
      else {
        helpImprint();
      }

      WebUI::webUIOutput.ensureImage();
      WebUI::webUIOutput.ensureMenu();
    }
    catch (int e) {
      WUI(WARN,"problem during imprint");
      helpImprint();
    }
  }
  void helpImprint() {
    WUI(INFO,"introduce edges by intersecting a specified geometry with portions of the existing surface\n" <<
        "  examples:\n" <<
        "    IMPRINT PLANE 0 0 0 1 0 0 ZONE_NAMES top,bottom\n" <<
        "    IMPRINT PLANE 0 0 0 0 1 0 TOL 1.0e-12 SPLIT_ZONES\n" <<
        "    IMPRINT WINDOW 0 1 0 0 0 0 1 0 0 1 1 1\n" <<
        "    IMPRINT CYL 0 0 0 0 0 1\n" <<
        "  more info: [$CWIKB:surfer_zoning]"
        );
  }

  void processIntersect(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available to intersect; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"INTERSECT intersects the triangles of 2 surfaces, joining them along their common edge(s)\n" <<
	  "and retaining only that portion of the resulting surface that is correctly walkable.\n" <<
	  "  examples:\n" <<
	  "    INTERSECT ZONES x0,y0 cylinder\n" <<
	  "    INTERSECT ZONES x0,y0 cylinder MODE SUBTRACT # default\n" <<
	  "    INTERSECT ZONES x0,y0 cylinder MODE ADD\n" <<
	  "  more info: knowledge base link here"
	  );
      return;
    }

    bool got_sz_ids = false;
    vector<int> sz0idVec,sz1idVec;

    string mode = "SUBTRACT"; // default mode is subtract

    int iarg = 0;
    while (iarg < param->size()) {
      const string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (parseSpecifiedSubzones(sz0idVec,token,param,iarg)){
        if (iarg > param->size()-1) {
          WUI(WARN,"INTERSECT: expecting two sets of comma-separated zone names or indices; skipping");
          return;
        }
        if (parseSpecifiedSubzones(sz1idVec,token,param,iarg)){
          got_sz_ids = true;
        }
      }
      else if (token == "MODE") {
        mode = param->getString(iarg++);
        cout << " > got mode: " << mode << endl;
      }
      else {
        WUI(WARN,"INTERSECT: skipping unrecognized token: " << token);
      }
    }

    if (!got_sz_ids) {
      cout << "ERROR: INTERSECT: missing SUBZONES. Check syntax." << endl;
      return;
    }

    ss.sz_flag.resize(ss.nsz);
    ss.sz_flag.setAll(-1);
    for (int ii = 0,nii=sz0idVec.size(); ii < nii; ++ii) {
      const int isz = sz0idVec[ii];
      if ((isz < 0)&&(isz >= ss.nsz)) {
        WUI(WARN,"INTERSECT: subzone index out of range: " << isz << "; skipping");
        return;
      }
      assert(ss.sz_flag[isz] == -1);
      ss.sz_flag[isz] = 0;
    }
    for (int ii = 0,nii=sz1idVec.size(); ii < nii; ++ii) {
      const int isz = sz1idVec[ii];
      if ((isz < 0)&&(isz >= ss.nsz)) {
        WUI(WARN,"INTERSECT: subzone index out of range: " << isz << "; skipping");
        return;
      }
      assert(ss.sz_flag[isz] == -1);
      ss.sz_flag[isz] = 1;
    }

    // now some tris are flagged 0, some are flagged 1, and most are -1. Intersect 0's and 1's...

    ss.intersectFlaggedSubzones(mode);

    WebUI::webUIOutput.ensureImage();
    WebUI::webUIOutput.ensureMenu();
  }

  void processHolePairing(Param * param, const bool help) {
    int iarg = 0;
    const double P_min  = param->getDouble(iarg++);
    const double P_max  = param->getDouble(iarg++);
    double ax[3]; FOR_I3 ax[i] = param->getDouble(iarg++);
    const double tol = param->getDouble(iarg++);

    cout << "Attempting tube reconstruction, candidate edge loops will have perimeters between " << P_min << " and " << P_max << endl;
    cout << "Pipe axis is " << COUT_VEC(ax) << endl;
    cout << "Matching tolerance is " << tol << endl;

    try {
      bool keep_going = true;
      int its = -1;
      ss.ensureOpenEdgeGroups();
      const int ngr_init = ss.getNOpenEdgeGroups();
      while (keep_going) {
        its++;
        cout << "Keep_going iteration # " << its << endl;
        ss.ensureOpenEdgeGroups();
        const int ngr = ss.getNOpenEdgeGroups();

        if (ngr == 0) {
          WUI(WARN,"There are no open edge groups. Returning.");;
          return;
        }//(ngr == 0)

        // Loop through all of the groups, searching for a valid candidate loop
        cout << "it = " << its << ", ngr = " << ngr << endl;
        for (int igr0 = 0;  igr0 < ngr; igr0++) {
          if ((ss.openEdgeGroupDataVec[igr0].length > P_min) && (ss.openEdgeGroupDataVec[igr0].length < P_max)) {
            const double pc0[3] = {ss.openEdgeGroupDataVec[igr0].xc[0] - ax[0],
                                   ss.openEdgeGroupDataVec[igr0].xc[1] - ax[1],
                                   ss.openEdgeGroupDataVec[igr0].xc[2] - ax[2]};
            const double pc1[3] = {ss.openEdgeGroupDataVec[igr0].xc[0] + ax[0],
                                   ss.openEdgeGroupDataVec[igr0].xc[1] + ax[1],
                                   ss.openEdgeGroupDataVec[igr0].xc[2] + ax[2]};

            cout << "Hole pairing candidate found!"                               << endl \
                 << "Group "      << igr0                                         << ", " \
                 << "Perim = "    << ss.openEdgeGroupDataVec[igr0].length         << ", " \
                 << "Centroid = " << COUT_VEC(ss.openEdgeGroupDataVec[igr0].xc)   << endl \
                 << "Expecting partner loop center to be at either:" << COUT_VEC(pc0) << " or " << COUT_VEC(pc1) << endl;

            int partner = -1;
            double err = HUGE_VAL;

            // Look for its match among the remaining open edge loops
            for (int igr1 = igr0+1; igr1 < ngr; igr1++) {
              if ((ss.openEdgeGroupDataVec[igr1].length > P_min) && (ss.openEdgeGroupDataVec[igr1].length < P_max)) {
                const double my_dist = min(DIST(ss.openEdgeGroupDataVec[igr1].xc, pc0), DIST(ss.openEdgeGroupDataVec[igr1].xc, pc1));
                if (my_dist < err) {
                  err = my_dist;
                  partner  = igr1;
                }//Improvement found
              }//(Perimeter guard
            }//(igr1 < ngr)

            if ((partner > 0) && ( err < tol)) {
              cout << "Valid matching partner found. Donut closing loops " << igr0 << " and " << partner << endl;
              vector <pair<int, int> > donut_loops;
              donut_loops.push_back(pair<int,int> (igr0, partner));

              ss.closeHalfEdgeDonutLoops(donut_loops, true);
              break;
            }//valid candidate, close loop and restart
          }// Valid perimeter
        }//(igr0 < ngr)

        if (its > ngr_init) {
          WUI(WARN,"Hole repair timed out");
          keep_going = false;
        }
      }//(keep_going)
    }//try
    catch (int e) {
      WUI(WARN, "Hole pairing failed!\n");
      helpProcessHolePairing();
    }
  }//processHolePairing()

  void helpProcessHolePairing() {
    WUI(INFO,"Finds two similar open edge loops and forms a tube to close them\n" <<
        "Arguments are minimum perimeter, maximum perimter, and expected tube centerline vector\n" <<
        "example:\n" <<
        "HOLE_PAIRING 1E-3 5E-3 0.0 1.0 0.0 ");
  }//helpProcessHolePairing()

  bool isActiveCommand(const string& cmd) {
    // grow this list of inactive commmands over time...
    // it is used by WRITE_JOURNAL, UNDO, REWIND, HISTORY/JOURNAL...
    if (MiscUtils::startsWith(cmd,"WRITE_IMAGE")||MiscUtils::startsWith(cmd,"write_image")||
        MiscUtils::startsWith(cmd,"INTERACTIVE")||MiscUtils::startsWith(cmd,"interactive")||
        (cmd == "I")||(cmd == "i")||
        MiscUtils::startsWith(cmd,"HELP")||MiscUtils::startsWith(cmd,"help")||
        MiscUtils::startsWith(cmd,"WRITE_JSON")||MiscUtils::startsWith(cmd,"write_json")||
        // ???...
        //MiscUtils::startsWith(cmd,"WRITE_SBIN")||MiscUtils::startsWith(cmd,"write_sbin")||
        //MiscUtils::startsWith(cmd,"RUN_DIAGNOSTICS")||MiscUtils::startsWith(cmd,"run_diagnostics")||
        MiscUtils::startsWith(cmd,"WRITE_JOURNAL")||MiscUtils::startsWith(cmd,"write_journal")||
        MiscUtils::startsWith(cmd,"WRITE_TECPLOT")||MiscUtils::startsWith(cmd,"write_tecplot")||
        MiscUtils::startsWith(cmd,"UNDO")||MiscUtils::startsWith(cmd,"undo")||
        MiscUtils::startsWith(cmd,"HISTORY")||MiscUtils::startsWith(cmd,"history")||
        MiscUtils::startsWith(cmd,"JOURNAL")||MiscUtils::startsWith(cmd,"journal")||
        MiscUtils::startsWith(cmd,"REWIND")||MiscUtils::startsWith(cmd,"rewind")||
        MiscUtils::startsWith(cmd,"BBOX")||MiscUtils::startsWith(cmd,"bbox"))
      return false;
    return true;
  }

  void processHistory(Param * param,const bool help) {

    if (help) {
      WUI(INFO,"HISTORY reports an indexed version of the parameter history for use with rewind");
      return;
    }

    WUI(INFO,"Current history (use REWIND <index> to return to a previous state):");
    int count = 0;
    for (list<string>::iterator iter = journal.begin(); iter != journal.end(); ++iter) {
      // skip all inactive commands...
      if (isActiveCommand(*iter)) {
        ++count;
        WUI(INFO,count << " " << *iter);
      }
    }

  }

  void processRewind(Param * param,const bool help) {

    if (help) {
      WUI(INFO,"REWIND <index> rewinds the state to the index reported in the active journal HISTORY");
      return;
    }

    // try to parse the index...
    int index;
    try {
      index = param->getInt();
    }
    catch(int e) {
      WUI(WARN,"REWIND expects a single index");
      return;
    }

    if (index < 0) {
      WUI(WARN,"REWIND <index> out of range. Use HISTORY to show valid command indices");
      return;
    }

    if (index == 0) {
      // index == 0 means just clear...
      Param journal_param("CLEAR");
      processParam(&journal_param);
    }
    else {
      int count = 0;
      list<string>::iterator iter_start = journal.begin();
      int count_start = 0;
      for (list<string>::iterator iter = journal.begin(); iter != journal.end(); ++iter) {
        // skip all inactive commands...
        if (isActiveCommand(*iter)) {
          if ((*iter == "CLEAR")||(*iter == "clear")) {
            iter_start = iter;
            ++iter_start; // advance past the clear. We add a new clear below
            count_start = count + 1;
          }
          ++count;
          if (count == index)
            break;
        }
      }
      if (index > count) {
        WUI(WARN,"REWIND <index> out of range. Use HISTORY to show valid command indices");
        return;
      }
      if (count_start == index) {
        // must have just been clear...
        Param journal_param("CLEAR");
        processParam(&journal_param);
      }
      else {
        assert(count_start < index);
        cout << " > REWIND will process CLEAR then all of:" << endl;
        list<string> undoList;
        count = count_start;
        for (list<string>::iterator iter = iter_start; iter != journal.end(); ++iter) {
          // skip all inactive commands...
          if (isActiveCommand(*iter)) {
            undoList.push_back(*iter);
            cout << "   > " << *iter << endl;
            ++count;
            if (count == index)
              break;
          }
        }
        Param journal_param("CLEAR");
        processParam(&journal_param);
        for (list<string>::iterator iter = undoList.begin(); iter != undoList.end(); ++iter) {
          if (createParamFromString(journal_param,*iter)) {
            processParam(&journal_param);
          }
          else {
            assert(0);
          }
        }
      }
    }
    WUI(INFO,"REWIND successful");

  }

  void processUndo(Param * param,const bool help) {

    if (help) {
      WUI(INFO,"UNDO uses the current journal to reset and redo everything except the last command");
      return;
    }

    // to undo, we are actually going to clear everything, then redo the current journal,
    // skipping passive commands...

    // look in the current journal for the most recent clear. This will be the case for
    // multiple undo commands... If the last command the user wants to undo is clear, then
    // we need a little more logic here...

    list<string>::iterator iter_start = journal.begin();
    for (list<string>::iterator iter = journal.begin(); iter != journal.end(); ++iter) {
      if ((*iter == "CLEAR")||(*iter == "clear")) {
        iter_start = iter;
        ++iter_start; // advance past the clear. We add a new clear below
      }
    }

    // now build the list of active commands...

    cout << " > UNDO will process CLEAR then all but the last of:" << endl;
    list<string> undoList;
    int count = 0;
    for (list<string>::iterator iter = iter_start; iter != journal.end(); ++iter) {
      if (isActiveCommand(*iter)) {
        undoList.push_back(*iter);
        cout << "   > " << *iter << endl;
        ++count;
      }
    }

    if (count <= 1) {
      WUI(WARN,"Nothing to UNDO");
      return;
    }

    Param journal_param("CLEAR");
    processParam(&journal_param);

    // now do everything but the last one...

    for (list<string>::iterator iter = undoList.begin(); iter != undoList.end(); ++iter) {
      if (createParamFromString(journal_param,*iter)) {
        processParam(&journal_param);
      }
      else {
        assert(0);
      }
      --count;
      if (count == 1)
        break;
    }

    WUI(INFO,"UNDO successfully applied");

  }

  void processApplySimilar(Param * param,const bool help) {
    if (help) {
      WUI(INFO,"obsolete functionality; you shouldn't need or use this"
          );
      return;
    }

    double area_tol = 0.02;
    int iarg = 0;
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "AREA_TOL") {
        area_tol = param->getDouble(iarg++);
        cout << "Using AREA_TOL: " << area_tol << endl;
      }
      else {
        cout << "Warning: unrecognized APPLY_SIMILAR param: " << token << ". Skipping" << endl;
      }
    }

    cout << "[ACTION] processApplySimilar: lastTemplate: \"" << lastTemplate << "\"" << endl;

    if (lastVec.size() == 2) {

      // ======================================
      // this is a 2-subsurface command...
      // ======================================

      assert(ss.szosz_i);
      assert(ss.szosz_v);

      // look for matching normal (within X% of area)...

      vector< pair<int,int> > similarVec;
      for (int isz=0, nsz=ss.subzoneDataVec.size(); isz < nsz; ++isz) {
        if (lastVec[0].first != isz) {
          if (fabs(lastVec[0].second.area - ss.subzoneDataVec[isz].area) < area_tol*lastVec[0].second.area) {
            //cout << "got an area match!" << endl;
            // now check normal alignment...
            const float dp = DOT_PRODUCT(lastVec[0].second.normal,ss.subzoneDataVec[isz].normal);
            const float mag1 = MAG(lastVec[0].second.normal);
            const float mag2 = MAG(ss.subzoneDataVec[isz].normal);
            if (dp > 0.95*mag1*mag2) {
              //cout << "got a normal match!" << endl;
              // now look for a connecting surface that matches lastVec[1]...
              int isz_match = -1;
              for (int sos = ss.szosz_i[isz]; sos != ss.szosz_i[isz+1]; ++sos) {
                const int isz_nbr = ss.szosz_v[sos];
                if (fabs(lastVec[1].second.area - ss.subzoneDataVec[isz_nbr].area) < area_tol*lastVec[1].second.area) {
                  if (isz_match == -1) {
                    isz_match = isz_nbr;
                  }
                  else {
                    //cout << "got multiple matches -- ambiguous -- skip" << endl;
                    isz_match = -2;
                  }
                }
              }
              if (isz_match >= 0) {
                similarVec.push_back(pair<int,int>(isz,isz_match));
              }
            }
          }
        }
      }

      cout << "similarVec.size(): " << similarVec.size() << endl;

      for (int ii = 0, ii_end = similarVec.size(); ii < ii_end; ++ii) {
        size_t pos0 = lastTemplate.find("%0");
        assert(pos0 != string::npos);
        size_t pos1 = lastTemplate.find("%1");
        assert(pos1 != string::npos);
        std::stringstream ss;
        ss << lastTemplate.substr(0,pos0);
        ss << similarVec[ii].first;
        ss << lastTemplate.substr(pos0+2,pos1-pos0-2);
        ss << similarVec[ii].second;
        ss << lastTemplate.substr(pos1+2);
        kfr.push_back(ss.str()); // here we really want push front
      }

    }
    else if (lastVec.size() == 1) {

      // ======================================
      // this is a 1-subsurface command...
      // ======================================

      vector<int> similarVec;
      for (int isz = 0,nsz=ss.subzoneDataVec.size(); isz < nsz; ++isz) {
        if (lastVec[0].first != isz) {
          if (fabs(lastVec[0].second.area - ss.subzoneDataVec[isz].area) < area_tol*lastVec[0].second.area) {
            //cout << "got an area match!" << endl;
            // now check normal alignment...
            const float dp = DOT_PRODUCT(lastVec[0].second.normal,ss.subzoneDataVec[isz].normal);
            const float mag1 = MAG(lastVec[0].second.normal);
            const float mag2 = MAG(ss.subzoneDataVec[isz].normal);
            if (dp > 0.95*mag1*mag2) {
              similarVec.push_back(isz);
            }
          }
        }
      }

      cout << "similarVec.size(): " << similarVec.size() << endl;

      for (int ii = 0, ii_end=similarVec.size(); ii < ii_end; ++ii) {
        size_t pos0 = lastTemplate.find("%0");
        assert(pos0 != string::npos);
        std::stringstream ss;
        ss << lastTemplate.substr(0,pos0);
        ss << similarVec[ii];
        ss << lastTemplate.substr(pos0+2);
        kfr.push_back(ss.str()); // here we really want push front
      }

    }
    else {

      cout << "Don't know how to handle lastVec.size(): " << lastVec.size() << endl;

    }

  }

  void processSelectSimilar(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"allows users to identify similar zones or edge loops to the specified zone or loop\n" <<
          "  examples:\n" <<
          "    SELECT_SIMILAR 45"
          );
      return;
    }

    // optional parameters to tweak resulting "similar-ness"
    bool isEdge = false;

    double diffTol = 0.02;  // allowable difference in surface areas/loop-length
    bool bySameZone = false;  // only select subzones with the same parent zone

    bool byArea = false;
    bool byNormal = false;
    double normalTol = 0.95;  // alignment of normals for similarity

    bool byLength = false;
    bool byCentroidMean = false;
    double centroidMeanTol = 0.02;
    bool byAngleSum = false;
    double angleSumTol = 0.95;

    int iarg = 0;
    int sz_index = param->getInt(iarg++);
    if (sz_index < 0) {
      CWARN("subzone index: " << sz_index << " invalid. Skipping");
      return;
    }
    else if (sz_index >= 0 && sz_index < ss.nsz) {
      // valid subzone
      isEdge = false;
    }
    else if (sz_index >= ss.nsz && sz_index < OPEN_E_SZ_MIN) {
      // neither a valid zone or edge loop
      CWARN("subzone index: " << sz_index << " invalid. Skipping");
      return;
    }
    else if (sz_index >= OPEN_E_SZ_MIN && sz_index < OPEN_E_SZ_MIN+ss.getNOpenEdgeGroups()) {
      // valid edge loop
      isEdge = true;
      sz_index -= OPEN_E_SZ_MIN;  // 0-indexed
    }
    else {
      CWARN("subzone index: " << sz_index << " invalid. Skipping");
      return;
    }

    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "DIFF_TOL") {
        diffTol = param->getDouble(iarg++);
      }
      else if (token == "NORMAL_TOL") {
        normalTol = param->getDouble(iarg++);
      }
      else if (token == "ANGLE_SUM_TOL") {
        angleSumTol = param->getDouble(iarg++);
      }
      else if (token == "CENTROID_MEAN_TOL") {
        centroidMeanTol = param->getDouble(iarg++);
      }
      else if (token == "BY") {
        const string byCsv = param->getString(iarg++);
        vector<string> byVec;
        MiscUtils::splitCsv(byVec,byCsv);

        for (int i = 0, limit = byVec.size(); i < limit; ++i) {
          const string criterion = byVec[i].c_str();
          if (criterion == "ZONE") {
            bySameZone = true;
          }
          else if (criterion == "AREA") {
            byArea = true;
          }
          else if (criterion == "NORMAL") {
            byNormal = true;
          }
          else if (criterion == "LENGTH") {
            byLength = true;
          }
          else if (criterion == "ANGLE_SUM") {
            byAngleSum = true;
          }
          else if (criterion == "CENTROID_MEAN") {  // maybe useful in both contexts?
            byCentroidMean = true;
          }
          else {
            CWARN("unrecognized BY criterion: " << criterion);
          }
        }
      }
      else {
        cout << "Warning: unrecognized APPLY_SIMILAR param: " << token << ". Skipping" << endl;
      }
    }

    COUT1("[ACTION] processSelectSimilar");
    COUT2(" > looking for similar " << ((isEdge) ? "edge loops":"subzones"));
    COUT2(" > searching only within same zone: " << ((bySameZone) ? "true":"false"));
    if (isEdge) {
      if (byLength) COUT2(" > length difference tolerance: " << diffTol);
      if (byAngleSum) COUT2(" > angle sum difference tolerance: " << angleSumTol);
      if (byCentroidMean) COUT2(" > mean dist. to centroid difference tolerance: " << centroidMeanTol);
    }
    else {
      if (byArea)   COUT2(" > area difference tolerance: " << diffTol);
      if (byNormal) COUT2(" > normal alignment tolerance: " << normalTol);
      // if (byCentroidMean) COUT2(" > mean dist. to centroid difference tolerance: " << centroidMeanTol);
    }

    if (isEdge && (byLength || byAngleSum || byCentroidMean)) {
      ss.selectSimilarOpenEdgeGroup(sz_index,byLength,byAngleSum,byCentroidMean,diffTol,angleSumTol,centroidMeanTol);
    }
    else if (!isEdge && (byArea || byNormal)) {
      // look for similar surface tri subzones
      ss.selectSimilarSubzone(sz_index,bySameZone,byArea,byNormal,diffTol,normalTol);
    }
    else {
      CWARN("criteria for selection not properly specified. Skipping");
    }
  }

  void processPipes(Param * param,const bool help) {
    surface_data_type = 2;
    surface_data_name = "";
    vector<int> zone_indices;
    int iarg = 0;
    while (iarg < param->size()) {
      const string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (parseSpecifiedZones(zone_indices,token,param,iarg)) {
        //no op; parsed info
      }
      else if ((token == "L_ON_D")||(token == "THETA0")||(token == "THETA1")) {
        surface_data_name = token;
      }
      else {
        CWARN("unrecognized token: " << token << "; skipping");
      }
    }
    if (surface_data_name == "") {
      CWARN("unrecognized pipe viz, should be L_ON_D, THETA0 or THETA1. returning...")
      return;
    }
    if (zone_indices.size() != 1) {
      CWARN("unrecognized single pipe-filled zone. returning...")
      return;
    }

    const int izone = zone_indices[0];
    ss.setSubzonesToZones();

    // split subzones into disjoint pipes...
    ss.sz_flag.resize(ss.nsz);
    ss.sz_flag.setAll(0);
    ss.sz_flag[izone] = 1;
    ss.splitFlaggedSubzonesIntoDisjointSubzones();

    if (ss.st_buf == NULL) ss.st_buf = new double[ss.nst];
    for (int ist = 0; ist < ss.nst; ++ist) 
      ss.st_buf[ist] = 0.0;

    ss.clearStoszSzosz();
    ss.ensureStoszSzosz();

    cout << " > disjoint-groups/pipes for zone: " << ss.szozn_i[izone+1]-ss.szozn_i[izone] << endl;

    for (int isz=ss.szozn_i[izone]; isz<ss.szozn_i[izone+1]; ++isz) {
      cout << " > tris for disjoint-group/pipe: " << ss.stosz_i[isz+1]-ss.stosz_i[isz] << endl;
      double l_on_d,theta[2]; // length over diameter, angles b/w pipe axis and face normals
      ss.characterizePipe(l_on_d,theta,isz);
      cout << " > l_on_d: " << l_on_d << " theta[2]: " << theta[0] << " " << theta[1] << endl;
      for (int toz = ss.stosz_i[isz]; toz != ss.stosz_i[isz+1]; ++toz) {
        const int ist = ss.stosz_v[toz];
        if (surface_data_name == "L_ON_D") {
          ss.st_buf[ist] = l_on_d;
        }
        else if (surface_data_name == "THETA0") {
          ss.st_buf[ist] = theta[0];
        }
        else if (surface_data_name == "THETA1") {
          ss.st_buf[ist] = theta[1];
        }
      }
    }

  }

  // this fella is one part of estimating a grid scale
  void processGetGapThicknesses(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"compute gap thicknesses around the surface"
          );
      return;
    }

    int iarg = 0;
    double planar_angle = 175.0; // use 175 deg by default
    double delta = 0.05; // must specify, this is the level0 delta (V/N)^(1/3)
    bool gap_viz = true; // no viz by default
    // bool cheap = false; // use cheap estimate ?
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      // angle for grouping almost planar tris...
      if (token == "PLANAR_ANGLE") {
        const double _planar_angle = param->getDouble(iarg++);
        if (_planar_angle >= 0.0 && _planar_angle <= 180.0) planar_angle = _planar_angle;
        else CWARN("PLANAR_ANGLE must be between [0,180]; using default value of " << planar_angle);
      }
      // delta for free stream...
      else if (token == "HCP_DELTA") {
        const double _delta = param->getDouble(iarg++);
        if (_delta > 0.0) delta = _delta;
        else CWARN("must specify positive HCP_DELTA; using default value of " << delta);
      }
      // vizualize gap
      else if (token == "VIZ") {
        gap_viz = param->getBool(iarg++);
      }
      // use cheap version?
      // else if (token == "CHEAP") {
      //   cheap = param->getBool(iarg++);
      // }
      else {
        CWARN("unrecognized token " << token << "; skipping");
      }
    }

    COUT1(" > computing surface gap thicknesses:")
      COUT1("    - hcp_delta estimate: " << delta);
    //if (cheap) {
    //  COUT1("    - planar_angle for coloring: " << planar_angle);
    //  ss.getGapThicknessesCheap(planar_angle,delta,gap_viz);
    // }
    //else
    ss.getGapThicknesses(delta,gap_viz);

    surface_data_type = 2;
    surface_data_name = "gap_D";
  }

  void processSurfaceCurvature(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"evaluate the surface curvature"
          );
      return;
    }

    int iarg = 0;
    int n_filter = 0;
    int k_method = 0;
    bool dump_histogram = false;
    double cdf_limit = 0.10;

    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "N_FILTER") {
        n_filter = param->getInt(iarg++);
      }
      else if (token == "METHOD") {
        const string method = param->getString(iarg++);
        if (method == "GAUSSBONNET") {
          k_method = 2;
        }
        else if (method == "OLD") {
          k_method = 1;
        }
        else if (method == "OLD_BY_SZ") {
          k_method = 0;
        }
      }
      else if (token == "DUMP_HIST") {
        dump_histogram = true;
      }
      else if (token == "CDF_LIMIT") {
        const double _cdf_limit = param->getDouble(iarg++);
        if (_cdf_limit >= 0.0 && _cdf_limit <= 1.0) cdf_limit = _cdf_limit;
        else CWARN("CFD_LIMIT must be between [0.0,1.0]; using default value of " << cdf_limit);
      }
      else {
        CWARN("unrecognized token: " << token << "; skipping");
      }
    }

    COUT1(" > computing surface curvature:")
      COUT1("    - method: " << (k_method==2 ? "Gauss-Bonnet":"normal-variation"));
    COUT1("    - grouping: " << (k_method==0 ? "by subzone":"entire surface"));
    if (dump_histogram) COUT1("    - writing histograms");
    if (k_method == 0) COUT1("    - cdf_limit (for r_curvature): " << cdf_limit);
    if (n_filter) COUT1("    - node-tri filtering passes: " << n_filter);


    if (k_method == 2) {
      ss.computeSurfaceCurvatureGaussBonnet(n_filter);
      surface_data_type = 1;
      surface_data_name = "Kappa_GB";
    }
    else if (k_method == 1) {
      ss.computeSurfaceCurvatureOld();
      surface_data_type = 2;
      surface_data_name = "Kappa_O";
    }
    else if (k_method == 0) {
      ss.computeSubzoneCurvature(cdf_limit,dump_histogram);
      surface_data_type = 2;
      surface_data_name = "Kappa_SZ";
    }
    else {
      COUT1(" > unrecognized curvature method specified; skipping");
    }
  }

  void processSuggestMesh(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"suggest mesh resolution length-scales per zone based on properties of the surface\n" <<
          "  examples:\n" <<
          "    SUGGEST_MESH HCP_DELTA 0.5"
          );
      return;
    }

    int iarg = 0;
    bool with_gap = false;
    double delta = 0.050;
    double crease_angle = 150.0;
    double cdf_limit = 0.10;
    double planar_angle = 175.0;

    UIMessage_ warn_msg(WARN); // TODO: update this to WUI format?
    UIMessage_ info_msg(INFO);
    info_msg.setClosable(true);
    // bool cheap = false; // default to full check
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "HCP_DELTA") {
        const double _delta = param->getDouble(iarg++);
        if (_delta > 0.0) delta = _delta;
        else {
          CWARN("HCP_DELTA must be greater than zero; using default value of " << delta);
          warn_msg.addLine(stringstream().flush() << "HCP_DELTA must be greater than zero; using default value of " << delta);
        }
      }
      else if (token == "CREASE_ANGLE") {
        const double _crease_angle = param->getDouble(iarg++);
        if (_crease_angle >= 0.0 || _crease_angle <= 180.0) crease_angle = _crease_angle;
        else {
          CWARN("CREASE_ANGLE must be between [0,180]; using default value of " << crease_angle);
          warn_msg.addLine("");
          warn_msg.addLine(stringstream().flush() << "CREASE_ANGLE must be between [0,180]; using default value of " << crease_angle);
        }
      }
      else if (token == "PLANAR_ANGLE") {
        const double _planar_angle = param->getDouble(iarg++);
        if (_planar_angle >= 0.0 || _planar_angle <= 180.0) planar_angle = _planar_angle;
        else {
          CWARN("PLANAR_ANGLE must be between [0,180]; using default value of " << planar_angle);
          warn_msg.addLine("");
          warn_msg.addLine(stringstream().flush() << "PLANAR_ANGLE must be between [0,180]; using default value of " << planar_angle);
        }
      }
      else if (token == "CDF_LIMIT") {
        const double _cdf_limit = param->getDouble(iarg++);
        if (_cdf_limit >= 0.0 && _cdf_limit <= 1.0) cdf_limit = _cdf_limit;
        else {
          CWARN("CFD_LIMIT must be between [0.0,1.0]; using default value of " << cdf_limit);
          warn_msg.addLine("");
          warn_msg.addLine(stringstream().flush() << "CFD_LIMIT must be between [0.0,1.0]; using default value of " << cdf_limit);
        }
      }
      else if (token == "WITH_GAP_THICKNESS") {
        with_gap = true;
      }
      else {
        CWARN("unrecognized token: " << token << "; skipping");
      }
    }

    if (warn_msg.size() > 0) WebUI::webUIOutput.add_(warn_msg);

    COUT1(" > computing mesh resolution estimates given:")
      COUT1("    - hcp_delta (Level 0): " << delta);
    COUT1("    - crease_angle (for sub-zoning): " << crease_angle);
    COUT1("    - cdf_limit (for r_curvature): " << cdf_limit);
    if (with_gap) COUT1("    - planar_angle (for gap_thickness): " << planar_angle);

    info_msg.addLine("Computing mesh resolution estimates given:");
    info_msg.addLine(stringstream().flush() << " - hcp_delta (Level 0): " << delta);
    info_msg.addLine(stringstream().flush() << " - crease_angle (for sub-zoning): " << crease_angle);
    info_msg.addLine(stringstream().flush() << " - cdf_limit (for r_curvature): " << cdf_limit);
    info_msg.addLine(stringstream().flush() << " - planar_angle (for gap_thickness): " << planar_angle);

    // split surface into subzones by crease angle
    assert(!ss.szost.isNull()); // subzone ids
    ss.zone_flag.setLength(ss.zoneVec.size());
    assert(ss.nsz == ss.szozn_i[ss.zoneVec.size()]);
    ss.sz_flag.setLength(ss.nsz);
    ss.sz_flag.setAll(1); // default is to split all
    ss.splitFlaggedSubzonesAtCreaseAngle(crease_angle);
    ss.clearSubzoneData();
    ss.clearZoneData();
    ss.ensureSubzoneData();
    ss.clearStoszSzosz();
    ss.ensureStoszSzosz();

    //   assert(0);
    // this routine is now the subzone based, but you should be able to use it
    // by resetting any subzone'ing to zones...
    //ss.splitSubzonesAtCreaseAngle(ss.szost,ss.zone_flag,crease_angle);  // produces zone-local sz-indexing
    // ss.buildSubzoneStuffFromZoneLocalSzFlag();
    //ss.buildSzoznFromLocalSzost();
    //ss.localSzostToGlobal();

    if (with_gap) ss.getGapThicknesses(delta,false); // no vizualization

    // compute r_curvature per-subzone
    ss.computeSubzoneCurvature(cdf_limit,false);  // no histogram output

    // now dump subzone information
    for (int isz=0; isz<ss.nsz; ++isz) {
      // cout << ss.subzoneDataVec[isz].curvature << " " << ss.subzoneDataVec[isz].gap_h << endl;
      // compute suggested level
      double dist_criteria;
      if (with_gap) dist_criteria = min(ss.subzoneDataVec[isz].curvature,ss.subzoneDataVec[isz].gap_h);
      else dist_criteria = ss.subzoneDataVec[isz].curvature;
      ss.subzoneDataVec[isz].level = (int)ceil(log2(delta/dist_criteria));
      if (ss.subzoneDataVec[isz].level < 0) ss.subzoneDataVec[isz].level = 0;
      COUT2("subzone[" << isz << "] metadata:");
      ss.subzoneDataVec[isz].dump();
    }

    // aggregate subzones based on similar information; perform zone suggestions
    // ss.makeSzostZoneLocal();
    //ss.globalSzostToLocal();

    set< pair<int,int> > level_sz_set;
    vector<int> lvosz;

    int current_sz=-1;

    for (int izone=0, nzn=ss.zoneVec.size(); izone<nzn; ++izone) {
      // COUT2("zone[" << izone << "]");

      for (int isz=ss.szozn_i[izone]; isz<ss.szozn_i[izone+1]; ++isz) {
        // COUT2("level_sz_set new pair: ("<< ss.subzoneDataVec[isz].level << ","<< isz <<")")
        level_sz_set.insert(pair<int,int> (ss.subzoneDataVec[isz].level,isz));
      }

      // agglomerate szs with the same level
      int current_level = -1;
      for (set< pair<int,int> >::iterator it=level_sz_set.begin(); it!=level_sz_set.end(); ++it) {
        if (current_level != it->first) {
          current_level = it->first;
          ++current_sz;
          lvosz.push_back(current_level);
        }

        for (int toz = ss.stosz_i[it->second]; toz != ss.stosz_i[it->second+1]; ++toz) {
          const int ist = ss.stosz_v[toz];
          ss.szost[ist] = current_sz;  // zone global indexing
        }
      }

      // correct
      ss.szozn_i[izone] = current_sz+1;

      level_sz_set.clear();
    }
    // rewind...
    for (int izone=ss.zoneVec.size(); izone > 0; --izone) ss.szozn_i[izone] = ss.szozn_i[izone-1];
    ss.szozn_i[0] = 0;

    ss.nsz = current_sz+1;
    assert(ss.szozn_i[ss.zoneVec.size()] == ss.nsz);

    // TODO maybe remove subzone data vec and just carry levels myself
    ss.clearStoszSzosz();
    ss.clearSubzoneData();
    ss.clearZoneData();
    ss.buildSubzoneData();

    // update level in sz metadata and parent zone metadata
    for (int isz=0; isz<ss.nsz; ++isz) {
      ss.subzoneDataVec[isz].level = lvosz[isz];
      const int zone_level = ss.zoneVec[ss.subzoneDataVec[isz].zone].getLevel();
      ss.zoneVec[ss.subzoneDataVec[isz].zone].setLevel(max(zone_level,lvosz[isz]));
    }
    ss.zone_level_hcp_delta = delta;

    // color st_buf with szs
    if (ss.st_buf == NULL) ss.st_buf = new double[ss.nst];

    for (int ist=0; ist<ss.nst; ++ist) {
      ss.st_buf[ist] = double(ss.subzoneDataVec[ss.szost[ist]].level);
    }
    surface_data_type = 2;
    surface_data_name = "Level";

    info_msg.addLine("");
    if (ss.nsz == int(ss.zoneVec.size())) {
      COUT2(" > current surface zones appear adequate for mesh resolution specification");
      info_msg.addLine("Looking good, the currently defined surface zones appear adequate for use with mesh resolution specification");
    }
    else {
      assert(ss.nsz > int(ss.zoneVec.size()));
      COUT2(" > analysis suggests zones should be further split to optimize mesh resolution specification;");
      COUT2("   inspect suggested zone decimation by looking at current subzones");
      info_msg.addLine("Our analysis suggests that the current zones should be further split to optimize mesh resolution specification. We have re-subzoned the surface with our suggestions on how to further decimate the current zones.");
      info_msg.addLine("");
      info_msg.addLine("Don't be shy, go inspect our suggestions for yourself. Trust me, they're probably good...at least OK. Worst case they're awful.");
    }
    WebUI::webUIOutput.add_(info_msg);
    ss.selectedSubzoneVec.clear();  // new subzones are generated so will lose old ones anyway
  }

  void processMorphTranslate(Param * param,const bool help) {

    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
      return;
    }

    if (help) {
      WUI(INFO,"MORPH_TRANSLATE ZONE <zone,another-zone,...> DX <dx> <dy> <dz>");
      return;
    }

    bool b_zones = false;
    bool b_dx = false;
    double dx[3];

    int iarg = 0;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if ((token == "ZONE")||(token == "ZONES")) {
	if (!ss.flagZones(param->getString(iarg++))) {
	  return;
	}
	b_zones = true;
      }
      else if (token == "DX") {
	FOR_I3 dx[i] = param->getDouble(iarg++);
	b_dx = true;
      }
      else {
	WUI(WARN,"unrecognized MORPH_TRANSLATE token: " << token);
      }
    }

    int ierr = 0;
    if (!b_zones) {
      WUI(WARN,"missing ZONES <comma-delimited-zone-names-or-numbers>");
      ierr = -1;
    }
    if (!b_dx) {
      WUI(WARN,"missing DX <dx> <dy> <dz>");
      ierr = -1;
    }
    if (ierr == -1)
      return;

    ss.morphTranslateFlaggedZones(dx);
    WUI(INFO,"morph applied successfully");
    WebUI::webUIOutput.ensureImage();

    // morphTranslateFlaggedZones leaves the morph factor in ss.sp_buf at the nodes...
    surface_data_type = 1; // sp
    surface_data_name = "sp_buf";

    return;





    // this stuff left for now...
    // this stuff left for now...
    // this stuff left for now...

    MPI_Pause("OKO:");

    iarg = 0;
    set<int> moved;
    set<int> freed;
    //double dx[3];
    int nrelax = 1000;
    bool use_bump = false;
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "MOVED") {
        const string splitZonesCsv = param->getString(iarg++);
        vector<string> splitZonesVec;
        MiscUtils::splitCsv(splitZonesVec,splitZonesCsv);
        for (int i = 0, limit = splitZonesVec.size(); i < limit; ++i) {
          const int zone = atoi(splitZonesVec[i].c_str());
          if (zone < 0 || zone >= ss.nsz) {
            CWARN("zone index out-of-bounds: " << zone);
          }
          else {
            moved.insert(zone);
          }
        }
      }
      else if (token == "FREED") {
        const string splitZonesCsv = param->getString(iarg++);
        vector<string> splitZonesVec;
        MiscUtils::splitCsv(splitZonesVec,splitZonesCsv);
        for (int i = 0, limit = splitZonesVec.size(); i < limit; ++i) {
          const int zone = atoi(splitZonesVec[i].c_str());
          if (zone < 0 || zone >= ss.nsz) {
            CWARN("zone index out-of-bounds: " << zone);
          }
          else {
            freed.insert(zone);
          }
        }
      }
      else if (token == "DX") {
        FOR_I3 dx[i] = param->getDouble(iarg++);
      }
      else if (token == "NRELAX") {
        nrelax = param->getInt(iarg++);
      }
      else if (token == "USE_BUMP") {
        use_bump = param->getBool(iarg++);
      }
      else {
        cout << "Warning: unrecognized MORPH_TRANSLATE_X param: " << token << ". Skipping." << endl;
      }
    }

    ss.morphTranslate(moved, freed, dx, nrelax, use_bump);

  }

  void processMorphRotateCylinder(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"morph capabilities or under development"
          );
      return;
    }

    int iarg = 0;
    int cyl_id = -1, top_id = -1, base_id = -1;
    double phi = 0.0; // default to no op
    int nrelax = 1000;
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "CYL" || token == "CYLINDER") {
        cyl_id = param->getInt(iarg++);
      }
      else if (token == "TOP") {
        top_id = param->getInt(iarg++);
      }
      else if (token == "BASE" || token == "BOTTOM") {
        base_id = param->getInt(iarg++);
      }
      else if (token == "PHI") {
        phi = -M_PI*param->getDouble(iarg++)/180.0;
      }
      else {
        cout << "Warning: unrecognized MORPH_ROTATE_CYLINDER param: " << token << ". Skipping." << endl;
        return;
      }
    }
    ss.morphRotateCylinder(cyl_id, top_id, base_id, phi, nrelax);

  }

  void processMorphRotateHalfCylinder(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"morph capabilities or under development"
          );
      return;
    }

    int iarg = 0;
    int cyl_id = -1, top_id = -1, base_id = -1, side_id = -1;
    double phi = 0.0; // default to no op
    int nrelax = 1000;
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "CYL" || token == "CYLINDER") {
        cyl_id = param->getInt(iarg++);
      }
      else if (token == "TOP") {
        top_id = param->getInt(iarg++);
      }
      else if (token == "BASE" || token == "BOT" ||  token == "BOTTOM") {
        base_id = param->getInt(iarg++);
      }
      else if (token == "SIDE") {
        side_id = param->getInt(iarg++);
      }
      else if (token == "PHI") {
        phi = -M_PI*param->getDouble(iarg++)/180.0;
      }
      else if (token == "NRELAX") {
        nrelax = param->getInt(iarg++);
      }
      else {
        cout << "Warning: unrecognized MORPH_ROTATE_CYLINDER param: " << token << ". Skipping." << endl;
        return;
      }
    }
    ss.morphRotateHalfCylinder(cyl_id, top_id, base_id, side_id, phi, nrelax);
  }

  void processMorphCorner(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    // TODO move to two step paradigm: select similar and then apply
    ss.ensureSubzoneData();

    int iarg = 0;
    int side = -1, top = -1;
    double dx, dw, dl;
    bool use_bump = false; // use linear unless bump specified
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "TOP") {
        top = param->getInt(iarg++);
      }
      else if (token == "SIDE") {
        side = param->getInt(iarg++);
      }
      else if (token == "DX") {
        dx = param->getDouble(iarg++);
      }
      else if (token == "DL") {
        dl = param->getDouble(iarg++);
        // init dw using dx
        dw = 2.0*dx;
      }
      else if (token == "USE_BUMP") {
        use_bump = param->getBool(iarg++);
      }
      else if (token == "DW") {
        dw = param->getDouble(iarg++);
      }
      else {
        cout << "Warning: unrecognized MORPH_CORNER param: " << token << ". Skipping." << endl;
      }
    }

    // store the "last" stuff...

    lastTemplate = "MORPH_CORNER TOP %0 SIDE %1";
    {
      stringstream oss; oss << dx;
      lastTemplate += " DX " + oss.str();
    }
    {
      stringstream oss; oss << dl;
      lastTemplate += " DL " + oss.str();
    }
    {
      stringstream oss; oss << use_bump;
      lastTemplate += " USE_BUMP " + oss.str();
    }
    {
      stringstream oss; oss << dw;
      lastTemplate += " DW " + oss.str();
    }
    cout << lastTemplate << endl;
    lastVec.resize(2);
    lastVec[0].first = top;
    lastVec[0].second = ss.subzoneDataVec[top];
    lastVec[1].first = side;
    lastVec[1].second = ss.subzoneDataVec[side];

    ss.morphCorner(top, side, dx, dl, use_bump, dw);

  }

  void processMorphStretch(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"morph capabilities or under development"
          );
      return;
    }

    // needed for lastVec
    ss.ensureSubzoneData();

    int iarg = 0;
    int top = -1;
    double dn, dl;
    bool use_bump = false; // use linear unless bump specified
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "TOP") {
        top = param->getInt(iarg++);
      }
      else if (token == "DN") {
        dn = param->getDouble(iarg++);
        // initialize dl usng dn
        dl = 2.0*dn;
      }
      else if (token == "USE_BUMP") {
        use_bump = param->getBool(iarg++);
      }
      else if (token == "DL") {
        dl = param->getDouble(iarg++);
      }
      else {
        cout << "Warning: unrecognized MORPH_STRETCH param: " << token << ". Skipping." << endl;
      }
    }

    // store the "last" stuff...

    lastTemplate = "MORPH_STRETCH TOP %0";
    {
      stringstream oss; oss << dn;
      lastTemplate += " DN " + oss.str();
    }
    {
      stringstream oss; oss << use_bump;
      lastTemplate += " USE_BUMP " + oss.str();
    }
    {
      stringstream oss; oss << dl;
      lastTemplate += " DL " + oss.str();
    }
    cout << lastTemplate << endl;
    lastVec.resize(1);
    lastVec[0].first = top;
    lastVec[0].second = ss.subzoneDataVec[top];

    ss.morphStretch(top, dn, use_bump, dl);

  }

  // Morph cylinder into a rotated superellipse using a  bump function
  void processMorphCylinder(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"morph capabilities or under development"
          );
      return;
    }

    int iarg = 0;
    int cyl_id = -1, top_id = -1, bot_id = -1;
    double x0, y0, r0, L, theta0, a, b;
    int n;
    bool z_cylinder = false; // default to arbitrary cylinder
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "CYL" || token == "CYLINDER") {
        cyl_id = param->getInt(iarg++);
      }
      else if (token == "TOP") {
        top_id = param->getInt(iarg++);
      }
      else if (token == "BOT" || token == "BOTTOM") {
        bot_id = param->getInt(iarg++);
      }
      else if (token == "CYL_X0") {
        x0 = param->getDouble(iarg++);
        z_cylinder = true;
      }
      else if (token == "CYL_Y0") {
        y0 = param->getDouble(iarg++);
      }
      else if (token == "CYL_R0") {
        r0 = param->getDouble(iarg++);
      }
      else if (token == "BUMP_L") {
        L = param->getDouble(iarg++);
      }
      else if (token == "SE_THETA0") {
        theta0 = M_PI/180.0*param->getDouble(iarg++);
      }
      else if (token == "SE_A") {
        a = param->getDouble(iarg++);
      }
      else if (token == "SE_B") {
        b = param->getDouble(iarg++);
      }
      else if (token == "SE_N") {
        n = param->getInt(iarg++);
        assert(n > 0);
      }
      else {
        cout << "Warning: unrecognized MORPH_CYLINDER param: " << token << ". Skipping." << endl;
        return;
      }
    }

    if (z_cylinder) {
      // morph cylinder into superellipse : (x/a)^n+(y/b)^n = 1
      ss.morphCylinderZ(cyl_id, top_id, bot_id, x0, y0, r0, L, theta0, a, b, n);
    }
    else {
      // place holder for now (just does some checks to see if it actually is a cylinder)
      ss.morphCylinder(cyl_id, top_id, bot_id, L, theta0, a, b, n);
    }

  }

  // Morph cylinder radius by some factor...
  void processMorphCylinderR(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"morph capabilities or under development"
          );
      return;
    }

    int iarg = 0;
    int cyl_id = -1;
    bool b_factor = false;
    double factor;

    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "SUBZONE") {
        cyl_id = param->getInt(iarg++);
      }
      else if (token == "FACTOR") {
        b_factor = true;
        factor = param->getDouble(iarg++);
      }
      else {
        cout << "Warning: unrecognized MORPH_CYLINDER_R param: " << token << ". Skipping." << endl;
        return;
      }
    }

    if ((cyl_id == -1)||(!b_factor)) {
      cout << "Warning: incorrect use of MORPH_CYLINDER_R: example\nMORPH_CYLINDER_R SUBZONE 9 FACTOR 1.02" << endl;
      return;
    }

    ss.morphCylinderR(cyl_id,factor);

  }

  void processMorphGaussianBump(Param * param,const bool help) {

    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
      return;
    }

    if (help) {
      WUI(INFO,"MORPH_GAUSSIAN_BUMP ZONE <zone,another-zone,...>");
      return;
    }

    bool b_zones = false;
    int iarg = 0;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if ((token == "ZONE")||(token == "ZONES")) {
	if (!ss.flagZones(param->getString(iarg++))) {
	  return;
	}
	b_zones = true;
      }
      else {
	WUI(WARN,"unrecognized MORPH_GAUSSIAN_BUMP token: " << token);
      }
    }

    int ierr = 0;
    if (!b_zones) {
      WUI(WARN,"missing ZONES <comma-delimited-zone-names-or-numbers>");
      ierr = -1;
    }
    if (ierr == -1)
      return;

    ss.morphGaussianBump();
    WUI(INFO,"morph applied successfully");
    WebUI::webUIOutput.ensureImage();

    return;

  }

  // Refine tris specified subzones to a specified length scale...
  void processRefineTris(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpRefineTris();
      return;
    }

    try {
      int iarg = 0;
      //bool b_st_flag = false;
      double delta = -1.0;

      // and zero the st_flag...
      ss.st_flag.setLength(ss.nst);
      ss.st_flag.setAll(0);

      vector<int> subzone_indices;

      double crease_angle = 0.0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));

        if (parseSpecifiedSubzones(subzone_indices,token,param,iarg)) {
          //no op; parsed info
        }
        else if (token == "DELTA") {
          delta = param->getDouble(iarg++);
        }
        else if (token == "CREASE_ANGLE") {
          crease_angle = param->getDouble(iarg++);
        }
        else {
          COUT2("Warning: unrecognized REFINE_TRIS param: " << token << "; skipping.");
          return;
        }
      }

      if (crease_angle == 0.0) {
        if (delta <= 0.0) {
          WUI(WARN,"DELTA was not specified or is not positive in REFINE_TRIS; skipping.");
          return;
        }

        if (subzone_indices.empty()) {
          // flag all tris
          for (int isz=0; isz < ss.nsz; ++isz) subzone_indices.push_back(isz);
        }

        // build the subzone to st structures...
        ss.ensureStoszSzosz();
        ss.flagTrisFromSubzoneVec(subzone_indices);
      }
      else {
        // flag tris at creases
        ss.ensureTeost();
        ss.sp_flag.setLength(ss.nsp);
        ss.sp_flag.setAll(0);
        const double dp_tol = cos((180.0-crease_angle)*M_PI/180.0);
        delta = HUGE_VAL;
        for (int ist = 0; ist < ss.nst; ++ist) {
          const double n[3] = TRI_NORMAL_2(ss.xsp[ss.spost[ist][0]],ss.xsp[ss.spost[ist][1]],ss.xsp[ss.spost[ist][2]]);
          const double area2 = MAG(n);
          FOR_I3 {
            int ist_nbr;
            if (ss.getAlignedTriNbr(ist_nbr,ist,i)) {
              const double n_nbr[3] = TRI_NORMAL_2(ss.xsp[ss.spost[ist_nbr][0]],ss.xsp[ss.spost[ist_nbr][1]],ss.xsp[ss.spost[ist_nbr][2]]);
              const double area_nbr2 = MAG(n_nbr);
              if (DOT_PRODUCT(n,n_nbr) < dp_tol*area2*area_nbr2) {
                ss.sp_flag[ss.spost[ist][i]] = 1;
                ss.sp_flag[ss.spost[ist][(i+1)%3]] = 1;
                const double vol_6 = fabs(SIGNED_TET_VOLUME_6(ss.xsp[ss.spost[ist_nbr][0]],ss.xsp[ss.spost[ist_nbr][1]],ss.xsp[ss.spost[ist_nbr][2]],ss.xsp[ss.spost[ist][(i+2)%3]]));
                delta = min(delta,vol_6/area_nbr2);
              }
            }
          }
        }
        for (int ist = 0; ist < ss.nst; ++ist) {
          FOR_I3 {
            if (ss.sp_flag[ss.spost[ist][i]] == 1) {
              ss.st_flag[ist] = 1;
              break;
            }
          }
        }
      }
      ss.refineFlaggedTris(delta);

      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem refining tris");
    }
  }
  void helpRefineTris() {
    WUI(INFO,"decimate selected tris until their edges are of a particular length scale"
        );
  }

  void processWriteJson(Param * param,const bool help) {
    if (help) {
      WUI(INFO,"a method to write surface data to a Json file. This is meant to be used by the Cascade App to communicate with surfer"
          );
      return;
    }

    try {
      bool b_hidesubzones = false;
      std::set<int> hiddenSubZonesSet;
      std::set<int> hiddenZonesSet;
      uint json_bits = 0;

      cout << "JUST GOT WRITE_JSON command: " << param->getString() << endl;
      assert(param->getString(0) == "NAME");
      string prefix = "surfer";

      //Check for hide zone information for bounding box
      //computation
      int iarg = 0;
      while (iarg < param->size()) {
        string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (token == "NAME") {
          prefix = param->getString(iarg++);
        }
        else if (token == "HIDE_SUBZONES") {
          token = param->getString(iarg++);

          vector<string> hiddenSubZonesVec;
          MiscUtils::splitCsv(hiddenSubZonesVec,token);
          if (hiddenSubZonesVec.size()) {
            b_hidesubzones = true;
            hiddenSubZonesSet.clear();
            for (int i=0, limit=hiddenSubZonesVec.size(); i<limit; ++i) {
              hiddenSubZonesSet.insert(atoi(hiddenSubZonesVec[i].c_str()));
            }
          }
        }
        else if (token == "HIDE_ZONES") {
          token = MiscUtils::toUpperCase(param->getString(iarg++));
          if (token == "ALL") {
            hiddenZonesSet.clear();
            for (int izn = 0, nzn = ss.zoneVec.size(); izn < nzn; ++izn) {
              hiddenZonesSet.insert(izn);
            }
          }
          else {
            vector<string> hiddenZonesVec;
            MiscUtils::splitCsv(hiddenZonesVec,token);
            if (hiddenZonesVec.size()) {
              b_hidesubzones = true;
              hiddenZonesSet.clear();
              for (int i=0, limit=hiddenZonesVec.size(); i<limit; ++i) {
                hiddenZonesSet.insert(atoi(hiddenZonesVec[i].c_str()));
              }
            }
          }
        }
        else if (token == "BITS") {
          json_bits = uint(param->getInt(iarg++));
        }
        else {
          if (mpi_rank == 0) cout << "Warning: skipping unrecognized WRITE_JSON token \"" << token << "\"" << endl;
        }
      }


      // write the JSON file...
      const string filename = prefix + ".json";
      MiscUtils::mkdir_for_file(filename);
      const string tmp_filename = MiscUtils::makeTmpPrefix(filename);
      FILE * fp = fopen(tmp_filename.c_str(),"w");
      assert(fp);

      double bbox[6];
      double xbb[3] = {0.0, 0.0, 0.0};
      double diag   = 0.0;
      bool b_surface_loaded = true;
      if (ss.nsp == 0 && ss.nst == 0){ //no surface has been loaded
        b_surface_loaded = false;
      }
      else{
        if (!b_hidesubzones) {
          ss.getBoundingBox(bbox);
          ss.getBoundingBoxCenter(xbb);
          diag = 2.0*ss.getBoundingBoxRmax();
        }
        else {
          // flag stores which zones should be hidden (i.e., not added)
          const int nzn = ss.szozn_i.getLength() - 1;
          const int nsz = ss.szozn_i[nzn];
          IntFlag show_sz_flag(nsz);
          show_sz_flag.setAll(1);
          if (b_hidesubzones) {
            std::set<int>::iterator it;
            // flag subzones based on hidden zones
            for (int izn = 0; izn < nzn; ++izn) {
              it = hiddenZonesSet.find(izn);
              if ( it!=hiddenZonesSet.end() ) {
                for (int iszn = ss.szozn_i[izn]; iszn < ss.szozn_i[izn+1]; ++iszn) show_sz_flag[iszn] = 0;  // hides tri from image
              }
            }

            for (int iszn = 0; iszn < nsz; ++iszn) {
              if (show_sz_flag[iszn] == 0) continue;  // skip if already hidden

              it = hiddenSubZonesSet.find(iszn);
              if ( it!=hiddenSubZonesSet.end() ) show_sz_flag[iszn] = 0;  // hides tri from image
            }
          }

          ss.getSubzoneBoundingBoxCenterAndDiagonal(bbox,xbb,diag,show_sz_flag);
        }
      }
      fprintf(fp,"{\n\"solver\":\"surfer\"");
      fprintf(fp,",\n\"b_surface_loaded\":%s",(b_surface_loaded?"true":"false"));
      const int manifold_state = ss.isManifold();
      switch (manifold_state) {
      case -1:
        fprintf(fp,",\n\"manifold_state\":\"unknown\"");
        break;
      case 1:
        fprintf(fp,",\n\"manifold_state\":\"true\"");
        break;
      default:
        fprintf(fp,",\n\"manifold_state\":\"false\"");
        break;
      }
      fprintf(fp,",\n\"BoundingBox\":[%f,%f,%f,%f,%f,%f]",bbox[0],bbox[1],bbox[2],bbox[3],bbox[4],bbox[5]);
      fprintf(fp,",\n\"BoundingBoxCentroid\":[%f,%f,%f]",xbb[0],xbb[1],xbb[2]);
      fprintf(fp,",\n\"BoundingBoxDiagonal\":%f",diag);
      fprintf(fp,",\n\"LightLocation\":[0.039503,0,1.01189,-0.530497,-0.57,0.581893]");

      if (!WebUI::webUIOutput.empty() || WebUI::webUIOutput.skipImage() || WebUI::webUIOutput.buildMenu()) {
        if (!WebUI::webUIOutput.empty()) assert(WebUI::webUIOutput.message_flag);  // consistency check

        if (WebUI::webUIOutput.skipImage()) {
          // the processParam routine that may have pushed messages into this JSON
          // may have set the newImage to false. If this is the case, we can skip
          // the image writing, which should share the same prefix as this JSON...
          b_skip_image = true;
          skip_image_name = filename;
        }
        WebUI::webUIOutput.writeJson(fp);  // only writes to the json the fields that are available
        WebUI::webUIOutput.clear();
      }

      // bit-based criterion
      if (json_bits & (1<<0)) {
        // siData for ui autocompletion
        writeAutoCompleteJson(fp);
      }

      fprintf(fp,",\n\"filetype\": \"%s\"", filetype.c_str());

      const int nzn=ss.zoneVec.size();
      // check if an empty surfer or not
      if (nzn>0) {
        fprintf(fp,",\n\"zoneNames\":[\n");
        for (int izn = 0; izn<nzn; ++izn) {
          if (izn == 0) fprintf(fp," \"%s\"",ss.zoneVec[izn].getName().c_str());
          else fprintf(fp,",\n \"%s\"",ss.zoneVec[izn].getName().c_str());
        }
        fprintf(fp,"\n]");
        fprintf(fp,",\n\"zoneSubzones\":[\n");
        assert((nzn+1) == ss.szozn_i.getLength());
        for (int i = 0, limit = nzn+1; i < limit; ++i) {
          if (i == 0) fprintf(fp," \"%d\"",ss.szozn_i[i]);
          else fprintf(fp,",\n \"%d\"",ss.szozn_i[i]);
        }
        fprintf(fp,"\n]");

        if (ss.selectedSubzoneVec.size()) {
          fprintf(fp,",\n\"selectedSubzones\":[\n");
          for (int i = 0, limit = ss.selectedSubzoneVec.size(); i < limit; ++i) {
            if (i == 0) fprintf(fp," \"%d\"",ss.selectedSubzoneVec[i]);
            else fprintf(fp,",\n \"%d\"",ss.selectedSubzoneVec[i]);
          }
          fprintf(fp,"\n]");
        }
        if (ss.b_update_hidden_subzones) {
          IntFlag sz_flag(ss.nsz);
          sz_flag.setAll(0);

          for (vector<int>::iterator it=ss.hiddenSubzoneVec.begin(); it!=ss.hiddenSubzoneVec.end(); ++it) {
            sz_flag[*it] = 1;
          }

          // flag fully hidden zones
          IntFlag zn_flag(ss.zoneVec.size());
          zn_flag.setAll(1);  // default to hidden
          for (int izn=0,nzn=ss.zoneVec.size(); izn<nzn; ++izn) {
            for (int isz=ss.szozn_i[izn],limit=ss.szozn_i[izn+1]; isz<limit; ++isz) {
              if (sz_flag[isz] == 0) {
                // sz not hidden, so zone cannot be hidden
                zn_flag[izn] = 0;
              }
              sz_flag[isz] = izn;  // hold zone of sz
            }
          }

          if (zn_flag.count()) {
            // need to hide full zones
            fprintf(fp,",\n\"hiddenZones\":[\n");
            // cout << "wrote hiddenZones to JSON: ";
            bool first = true;
            for (int i = 0, limit = zn_flag.size(); i < limit; ++i) {
              if (zn_flag[i]) {
                if (first) {
                  fprintf(fp," \"%d\"",i);
                  first = false;
                }
                else fprintf(fp,",\n \"%d\"",i);
                // cout << i << ",";
              }
            }
            // cout << endl;
            fprintf(fp,"\n]");

            // prune from hiddenSubzonesVec so not doubly hidden
            vector<int>::iterator it=ss.hiddenSubzoneVec.begin();
            while (it!=ss.hiddenSubzoneVec.end()) {
              if (zn_flag[sz_flag[*it]]) it = ss.hiddenSubzoneVec.erase(it);  // deletes and points to next el, so don't increment
              else ++it;
            }
          }

          if (ss.hiddenSubzoneVec.size()) {
            fprintf(fp,",\n\"hiddenSubzones\":[\n");
            // cout << "wrote hiddenSubzones to JSON: ";
            bool first = true;
            for (int i = 0, limit = ss.hiddenSubzoneVec.size(); i < limit; ++i) {
              if (first) {
                fprintf(fp," \"%d\"",ss.hiddenSubzoneVec[i]);
                first = false;
              }
              else fprintf(fp,",\n \"%d\"",ss.hiddenSubzoneVec[i]);
              // cout << ss.hiddenSubzoneVec[i] << ",";
            }
            // cout << endl;
            fprintf(fp,"\n]");
            ss.hiddenSubzoneVec.clear();
          }

          ss.b_update_hidden_subzones = false;
        }
        fprintf(fp,",\n\"zoneBoundaryConditions\":[\n");
        for (int izn = 0; izn<nzn; ++izn) {
          if (izn == 0) fprintf(fp," \"%s\"",ss.zoneVec[izn].getBC().c_str());
          else fprintf(fp,",\n \"%s\"",ss.zoneVec[izn].getBC().c_str());
        }
        fprintf(fp,"\n]");
        // periodic info
        const int npt = PeriodicData::periodicTransformVec.size();
        if (npt > 0) {
          fprintf(fp,",\n\"zonePeriodicBits\":[\n");
          for (int izn = 0; izn < nzn; ++izn) {
            if (izn == 0) fprintf(fp," \"%d\"",ss.zoneVec[izn].getPeriodicBits());
            else fprintf(fp,",\n \"%d\"",ss.zoneVec[izn].getPeriodicBits());
          }
          fprintf(fp,"\n]");
          fprintf(fp,",\n\"periodicTransformKind\":[\n");
          for (int ipt = 0; ipt < npt; ++ipt) {
            if (ipt == 0) fprintf(fp," \"%d\"",PeriodicData::periodicTransformVec[ipt].getKind());
            else fprintf(fp,",\n \"%d\"",PeriodicData::periodicTransformVec[ipt].getKind());
          }
          fprintf(fp,"\n]");
          fprintf(fp,",\n\"periodicTransformData\":[\n");
          for (int ipt = 0; ipt < npt; ++ipt) {
            double pt_data[3]; PeriodicData::periodicTransformVec[ipt].getData(pt_data);
            if (ipt == 0) fprintf(fp," [%f,%f,%f]",pt_data[0],pt_data[1],pt_data[2]);
            else fprintf(fp,",\n [%f,%f,%f]",pt_data[0],pt_data[1],pt_data[2]);
          }
          fprintf(fp,"\n]");
        }
        if (surface_data_name != "") {
          fprintf(fp,",\n\"boundaryVarsAndFxns\":{\n");
          fprintf(fp," \"D1\":[\n");
          if ((surface_data_type == 1)||(surface_data_type == 2))
            fprintf(fp,"  \"%s\"",surface_data_name.c_str());
          fprintf(fp,"\n ],\n");
          fprintf(fp," \"D2\":[\n");
          if ((surface_data_type == 3)||(surface_data_type == 4))
            fprintf(fp,"  \"%s\"",surface_data_name.c_str());
          fprintf(fp,"\n  ]");
          fprintf(fp,"\n }");
        }
      }  // if nzn>0
      fprintf(fp,"\n}\n");
      fclose(fp);
      remove(filename.c_str());
      rename(tmp_filename.c_str(),filename.c_str());

      // clear this vector so subsequent calls don't populate it
      ss.selectedSubzoneVec.clear();
    }
    catch (int e) {
      WUI(WARN,"problem writing json");
    }
  }

  void processWriteImage(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"allows a user to export a png image of the current surface\n" <<
          "  examples:\n" <<
          "    WRITE_IMAGE NAME images/surface GEOM Y_PLANE_FRAC 0.5 SIZE 1200 600 VAR SURFACE\n" <<
          "  more info: [$CWIKB:write_image]"
          );
      return;
    }

    try {
      // note - any INTERVAL param in the WRITE_IMAGE command will be ignored here...
      WriteImageData wid;
      if (wid.init(param) == -2) return;

      // when interacting with the client, sometimes an image is not required,
      // even though one was requested.
      if (b_skip_image && wid.b_name && (skip_image_name == wid.name)) {
        // skip this immage...
        if (mpi_rank == 0) cout << " > skipping this image" << endl;
        b_skip_image = false;
        return;
      }

      CtiScene scene(wid);

      scene.setRmax(ss.getBoundingBoxRmax());
      double bbox[6];
      ss.getBoundingBox(bbox);
      scene.setBoundingBbox(bbox[0],bbox[1],bbox[2],bbox[3],bbox[4],bbox[5]);
      double center[3];
      ss.getBoundingBoxCenter(center);
      scene.setCenter(center[0],center[1],center[2]);

      for (int izn=0,nzn=ss.zoneVec.size(); izn<nzn; ++izn) scene.addZoneName(ss.zoneVec[izn].getName());
      scene.convertHiddenZoneNamesToIndices();

      // set hidden subzone info if present
      ss.hiddenSubzoneVec.clear();
      if (scene.hasHiddenSubzones()) {
        int * sz_flag_tmp = new int[ss.nsz];
        scene.getHiddenSubzones(sz_flag_tmp,ss.nsz);
        for (int isz=0; isz<ss.nsz; ++isz) {
          if (sz_flag_tmp[isz]) ss.hiddenSubzoneVec.push_back(isz);
        }
        DELETE(sz_flag_tmp);
      }

      // treat zones as subzones so hiding gets processed by one routine
      // report back new fully hidden zones in json
      if (scene.hasHiddenZones()) {
        const int nzn = ss.zoneVec.size();
        int * zn_flag_tmp = new int[nzn];
        scene.getHiddenZones(zn_flag_tmp,nzn);
        for (int izn=0; izn<nzn; ++izn) {
          if (zn_flag_tmp[izn]) {
            for (int isz=ss.szozn_i[izn], limit=ss.szozn_i[izn+1]; isz<limit; ++isz) ss.hiddenSubzoneVec.push_back(isz);
          }
        }
        DELETE(zn_flag_tmp);
      }

      scene.initCanvas();

      if (scene.hasGeomPlane()) {
        scene.setGeomPlaneAsDataPlane();
        if (scene.blankDataPlane()) scene.addGeomPlaneAsBlankPlane(); // should be changed to back in scene/canvas
      }

      // and render...
      // if (ss.eoi_type != SimpleSurface::NONE) {
      //   // special rendering if dynamic edges have been requested edge rendered...
      //   ss.ensureDynamicEdgeGroups(ss.eoi_type);
      //
      //   //HACK currently znost is not used within scene, only subzone
      //   // we can create an int to carry the zone edges (or dynamic edges) that should be highlighted/focused
      //   int * focus_edges = new int[ss.nst];
      //   for (int ist = 0; ist < ss.nst; ++ist) focus_edges[ist] = 0;
      //
      //   // update focus edges with dynamic edge info
      //   for (map<pair<int,int>,int>::iterator it=ss.eoi_to_group.begin(); it!=ss.eoi_to_group.end(); ++it) {
      //     assert(it->second >= 0);  //negative group encountered
      //     const int ist = it->first.first;
      //     const int edge = it->first.second;
      //     assert(ist < ss.nst); // invalid  ist
      //     assert((edge < 3) &&(edge >= 0));  // bad edge index
      //     focus_edges[ist] |= (1 << edge);  // flag zone boundary edge
      //   }
      //
      //   scene.addSurfaceTrisFocusEdge(ss.xsp, ss.spost, focus_edges, ss.szost, ss.szozn_i, ss.nst);  // replace znost with focus edge indicator....
      //   DELETE(focus_edges);
      // }
      // else {
      //   ss.clearDynamicEdgeGroups();
      //   scene.addSurfaceTris(ss.xsp, ss.spost, ss.znost, ss.szost, ss.szozn_i, ss.nst);
      // }

      string surface_var;
      if (scene.getSurfaceVar(surface_var)&&(surface_var == surface_data_name)&&(surface_data_type != 0)) {

        if (surface_data_type == 1) {
          scene.addSurfaceTrisWithData(ss.xsp,ss.spost,ss.szost,ss.szozn_i,ss.nst,ss.sp_buf,surface_data_type);
        }
        else if (surface_data_type == 2) {
          scene.addSurfaceTrisWithData(ss.xsp,ss.spost,ss.szost,ss.szozn_i,ss.nst,ss.st_buf,surface_data_type);
        }
        else {
          // TODO sp_buf3 and st_buf3
          scene.addSurfaceTrisWithData(ss.xsp,ss.spost,ss.szost,ss.szozn_i,ss.nst,NULL,0);
        }

      }
      else {
        scene.addSurfaceTrisWithData(ss.xsp,ss.spost,ss.szost,ss.szozn_i,ss.nst,NULL,0);
      }

      if (scene.hasEdges()) {
        // TODO: add in room for periodic edges
        ss.ensureOpenEdgeGroups();

        // edge indexing for selection is as follows:
        // 0: regular edge (aligned)
        // -1: misaligned nbr edge
        // -2,-3,-4,...: open edge groups are -2 indexed
        // 1,2,3,...: multi-edges are +1 indexed
        // we need to decide how to do other edge showing:
        // * imprinted edges
        // * zone boundary edges
        int (*teost)[3] = new int[ss.nst][3];
        for (int ist = 0; ist < ss.nst; ++ist) {
          FOR_I3 {
            int igr;
            if (ss.getOpenEdgeGroup(igr,ist,i)) {
              teost[ist][i] = -igr-2;
            }
            else if (ss.isEdgeMisaligned(ist,i)) {
              teost[ist][i] = -1;
            }
            else {
              teost[ist][i] = 0;
            }
          }
        }

        // special behavior for 1+ indexed edges
        if (ss.eoi_type != SimpleSurface::NONE) {
          ss.ensureDynamicEdgeGroups(ss.eoi_type);
          // render dynamic edges selectable
          for (map<uint,int>::iterator it=ss.eoi_to_group.begin(); it!=ss.eoi_to_group.end(); ++it) {
            assert(it->second >= 0);  //negative group encountered
            uint ist; uchar edge;
            ss.unpackEdge(ist,edge,it->first);
            assert(int(ist) < ss.nst); // invalid  ist
            assert(edge < 3);  // bad edge index
            teost[ist][edge] = it->second+1;  // +1 indexed
          }
        }
        // else {
        //   for (int ist = 0; ist < ss.nst; ++ist) {
        //     FOR_I3 {
        //       int ime, orient_ime;
        //       if (ss.isEdgeMulti(ime,orient_ime,ist,i)) {
        //         teost[ist][i] = ime+1;
        //       }
        //     }
        //   }
        // }

        scene.addEdges(ss.xsp, ss.spost, teost, ss.znost, ss.szost, ss.szozn_i, ss.nst);
        delete[] teost;
      }

      scene.writeImage();

    }
    catch (int e) {
      WUI(WARN,"problem writing image");
    }
  }

  void processSplitFlag(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpSplitFlag();
      return;
    }

    try {
      assert(!ss.szost.isNull()); // subzone ids
      assert(ss.nsz == ss.szozn_i[ss.zoneVec.size()]);
      ss.sz_flag.setLength(ss.nsz);
      ss.sz_flag.setAll(0); // default is to split none

      string destination_zone = "";
      string source_zone      = "";
      int split_mode            = 0;
      double split_data[6]    = {0.0,0.0,0.0,0.0,0.0};
      //bool b_all = false;
      vector<int> subzone_indices;

      int iarg = 0;
      while (iarg < param->size()) {
        string token = MiscUtils::toUpperCase(param->getString(iarg++));

        if (parseSpecifiedSubzones(subzone_indices,token,param,iarg)) {
          //no op; parsed info
        }
        else if (token == "DESTINATION_ZONE") {
          destination_zone = param->getString(iarg++);
        }
        else if (token == "BY_ANGLE") {
          split_data[0] = param->getDouble(iarg++);
          split_mode |= (1<<7);
        }
        else if (token == "BY_DISJOINT") {
          split_mode |= (1<<8);
        }
        else if (token == "BY_COORDINATE")  {
          split_mode |= (1<<9);
          token = MiscUtils::toUpperCase(param->getString(iarg++));
          if ( token == "X") {
            split_mode |= (1<<1);   // keeping track of the splitting via bits...
          }
          else if ( token == "Y") {
            split_mode |= (1<<2);
          }
          else if ( token == "Z") {
            split_mode |= (1<<3);
          }
          else {
            CWARN( " > unknown zone coordinate direction : " << token );
            return;
          }

          token = MiscUtils::toUpperCase(param->getString(iarg++));
          if ( token == "ABOVE") {
            split_mode |= (1<<4);
          }
          else if ( token == "BELOW") {
            split_mode |= (1<<5);
          }
          else {
            CWARN(" > Unrecognized token " << token << ". You must specify ABOVE or BELOW to split BY_COORDINATE; skipping");
            return;
          }
          split_data[0] = param->getDouble(iarg++);

        }
        else if ( token == "MARCH_FROM_POINT") {

          split_mode |= (1<<10);

          FOR_I3 split_data[i] = param->getDouble(iarg++);

          split_data[3] = 0.865; // default angle tolerance set... (150 degree crease)
          while (iarg < param->size()) {
            const string tok = MiscUtils::toUpperCase(param->getString(iarg++));
            if (tok == "ANGLE_TOL") {
              const double angle = param->getDouble(iarg++);
              split_data[3]      = cos(angle*M_PI/180.0);
            }
            else {
              CWARN("unrecognized  MARCH_FROM_POINT parameter " << tok << "; skipping");
            }
          }
        }
        else if (token == "BY_TRI") {
          split_mode |= (1<<11);
        }
        else if (token == "IN_SPHERE") {
          split_data[0] = param->getDouble(iarg++); // x
          split_data[1] = param->getDouble(iarg++); // y
          split_data[2] = param->getDouble(iarg++); // z
          split_data[3] = param->getDouble(iarg++); // r
          split_mode |= (1<<12);
        }
        else if (token == "BY_PLANE") {
          split_data[0] = param->getDouble(iarg++); // x
          split_data[1] = param->getDouble(iarg++); // y
          split_data[2] = param->getDouble(iarg++); // z
          split_data[3] = param->getDouble(iarg++); // nx
          split_data[4] = param->getDouble(iarg++); // ny
          split_data[5] = param->getDouble(iarg++); // nz
          split_mode |= (1<<13);
        }
        else {
          CWARN("Unknown SPLIT parameter " << token << "; skipping");
        }
      }

      // flag subzones
      ss.sz_flag.resize(ss.nsz);
      ss.sz_flag.setAll(0);
      for (vector<int>::iterator it=subzone_indices.begin(); it!=subzone_indices.end(); ++it) ss.sz_flag[*it] = 1;
      if (subzone_indices.empty()) ss.sz_flag.setAll(1);

      if ((ss.sz_flag.count() == 0) && (split_mode != 0)) {
        COUT2(" > no valid subzones were selected for splitting; skipping");
        return;
      }

      if (split_mode & (1<<7)) {
        ss.splitFlaggedSubzonesAtCreaseAngle(split_data[0]);
      }
      else if (split_mode & (1<<8)) {
        ss.splitFlaggedSubzonesIntoDisjointSubzones();
      }
      else if (split_mode & (1<<9)) {
        ss.splitFlaggedSubzonesByCoordinate(split_mode,split_data,destination_zone);
      }
      else if (split_mode & (1<<10)) {
        ss.splitFlaggedSubzonesByMarch(split_data,destination_zone);
      }
      else if (split_mode & (1<<11)) {
        ss.splitFlaggedSubzonesIntoTris();
      }
      else if (split_mode & (1<<12)) {
        ss.splitFlaggedSubzonesInSphere(split_data,split_data[3]); // x[3],radius
      }
      else if (split_mode & (1<<13)) {
        ss.splitFlaggedSubzonesByPlane(split_data,split_data+3); // x[3],n[3]
      }
      else {
        ss.setSubzonesToZones();
      }
      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem splitting subzones");
      helpSplitFlag();
    }
  }
  void helpSplitFlag() {
    WUI(INFO,"splits surface zones into smaller subzones based on differing criteria\n" <<
        "  examples:\n" <<
        "    SPLIT_ZONE ALL BY_ANGLE 155.0 (based on crease angle)\n" <<
        "    SPLIT_ZONE ZONE_NAMES left,right BY_DISJOINT (by disconnected pieces)"
        );
  }

  void processMakeInjector(Param * param,const bool help) {
    if (help) {
      WUI(INFO,"generate a specific surface for injector geometries"
          );
      return;
    }

    assert(!ss.szost.isNull()); // subzone ids
    assert(ss.nsz == ss.szozn_i[ss.zoneVec.size()]);

    vector<int> subzone_indices;
    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "SUBZONES" || token == "FOR_ZONES") {
        const string subzonesCsv = param->getString(iarg++);
        MiscUtils::splitCsv(subzone_indices,subzonesCsv);
      }
      else {
        CWARN("Unknown SPLIT parameter " << token << "; skipping");
      }
    }

    // flag subzones
    ss.sz_flag.resize(ss.nsz);
    ss.sz_flag.setAll(0);
    for (vector<int>::iterator it=subzone_indices.begin(); it!=subzone_indices.end(); ++it) {
      assert((*it >= 0)&&(*it < ss.nsz));
      ss.sz_flag[*it] = 1;
    }
    ss.makeInjectorForFlaggedSubzones();

  }

  void processWindowSelect(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"allows Cascade App to query subzones within a fixed region of the App screen"
          );
      return;
    }

    try {
      double window[4][3];
      bool b_window = false;
      bool b_strictly_inside = false;
      bool b_edges = false;
      bool b_open = false;

      int iarg = 0;
      while (iarg < param->size()) {
        string token = param->getString(iarg++);
        if (token == "WINDOW") {
          FOR_I3 window[0][i] = param->getDouble(iarg++);
          FOR_I3 window[1][i] = param->getDouble(iarg++);
          FOR_I3 window[2][i] = param->getDouble(iarg++);
          FOR_I3 window[3][i] = param->getDouble(iarg++);
          b_window = true;
        }
        else if (token == "SUBZONES") {
          // add already selected subzones to vec
          const string subzonesCsv = param->getString(iarg++);
          MiscUtils::splitCsv(ss.selectedSubzoneVec,subzonesCsv);
        }
        else if (token == "STRICTLY_INSIDE") {
          b_strictly_inside = true;
        }
        else if (token == "DYNAMIC_EDGES") {
          b_edges = true;
        }
        else if (token == "OPEN_EDGES") {
          b_open = true;
          b_edges = true;
        }
        else {
          CWARN("Unknown WINDOW_SELECT parameter " << token << "; skipping");
        }
      }

      if (b_window) {
        if (b_edges)
          ss.selectEdgesInWindow(window,b_strictly_inside,b_open);
        else
          ss.selectZonesInWindow(window,b_strictly_inside);
      }
      else {
        CWARN("WINDOW_SELECT requires WINDOW x0 y0 z0 x1 y1 z1 x2 y2 z2 x3 y3 z3. returning...");
      }
    }
    catch (int e) {
      WUI(WARN,"problem during WINDOW_SELECT");
    }
  }

  void processMakePart(Param * param,const bool help) {
    if (help) {
      WUI(INFO,"developmental part construction"
          );
      return;
    }
    // the first param needs to be the type of part...

    if (param->size() < 1) {
      CWARN("MAKE_PART expecting part keyword. Skipping");
      return;
    }

    const string part_keyword = param->getString();
    if (part_keyword == "BL") {
      processMakePartBl(param,false);
    }
    else {
      CWARN("MAKE_PART: unrecognized part keyword: " << part_keyword << ". Skipping");
      return;
    }

  }

  void processMakePartBl(Param * param,const bool help) {
    if (help) {
      WUI(INFO,"developmental part generation"
          );
      return;
    }
    cout << "processMakePartBl" << endl;

    bool b_subzones = false;
    //bool b_name = false;
    string name;
    bool b_dn = false;
    double dn;
    bool b_dt = false;
    double dt;
    bool b_n = false;
    int n;

    assert(!ss.szost.isNull()); // subzone ids
    assert(ss.nsz == ss.szozn_i[ss.zoneVec.size()]);
    vector<int> subzone_indices;
    int iarg = 1;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "SUBZONES") {
        // expect a comma-delimited set of integers...
        b_subzones = true;
        const string subzonesCsv = param->getString(iarg++);
        MiscUtils::splitCsv(subzone_indices,subzonesCsv);
      }
      else if (token == "DN") {
        b_dn = true;
        dn = param->getDouble(iarg++);
      }
      else if (token == "DT") {
        b_dt = true;
        dt = param->getDouble(iarg++);
      }
      else if ((token == "N")||(token == "NLAYERS")) {
        b_n = true;
        n = param->getInt(iarg++);
      }
      else if (token == "NAME") {
        //b_name = true;
        name = param->getString(iarg++);
      }
      else {
        CWARN("Unknown MAKE_PART BL parameter " << token << "; skipping");
      }
    }

    int ierr = 0;
    if (!b_dn) {
      cout << "Error: MAKE_PART BL missing normal spacing DN <double>" << endl;
      ierr = -1;
    }
    if (!b_dt) {
      cout << "Error: MAKE_PART BL missing tangential spacing DT <double>" << endl;
      ierr = -1;
    }
    if (!b_n) {
      cout << "Error: MAKE_PART BL missing N <int>" << endl;
      ierr = -1;
    }
    if (!b_subzones) {
      cout << "Error: MAKE_PART BL missing SUBZONES <comma-deliminted-indices>" << endl;
      ierr = -1;
    }
    if (ierr != 0) {
      cout << "Skipping MAKE_PART BL. Fix errors and rerun" << endl;
      return;
    }

    // flag subzones
    ss.sz_flag.resize(ss.nsz);
    ss.sz_flag.setAll(0);
    for (vector<int>::iterator it=subzone_indices.begin(); it!=subzone_indices.end(); ++it) {
      assert((*it >= 0)&&(*it < ss.nsz));
      ss.sz_flag[*it] = 1;
    }
    ss.makePartBlForFlaggedSubzones(name,dn,dt,n);

  }

  void processFlagTris(Param * param,const bool b_help) {

    if (b_help) {
      WUI(INFO,"flags surface tris in requested zones\n" <<
          "  examples:\n" <<
          "    FLAG_TRIS ZONES x0,x1\n" <<
          "  more info: [$CWIKB:surfer_tri_flagging]"
          );
      return;
    }

    vector<int> zone_indices;
    bool b_add = false;

    int iarg = 0;
    while (iarg < param->size()) {
      const string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (parseSpecifiedZones(zone_indices,token,param,iarg)) {
        // no op
      }
      else if (token == "ADD") {
        b_add = true;
      }
      else {
        WUI(WARN,"unrecognized FLAG_TRIS token " << token << ". Check syntax and try again.");
        return;
      }
    }

    if (zone_indices.empty()) {
      WUI(WARN,"no zones or subzones selected");
      return;
    }

    ss.flagTrisInZones(zone_indices,b_add);

  }

  void processFlagTrisTouching(Param * param,const bool b_help) {

    if (b_help) {
      WUI(INFO,"flags surface tris touching requested zones at any edge or node\n" <<
          "  examples:\n" <<
          "    FLAG_TRIS_TOUCHING ZONES x0,x1\n" <<
          "  more info: [$CWIKB:surfer_tri_flagging]"
          );
      return;
    }

    vector<int> zone_indices;
    bool b_add = false;

    int iarg = 0;
    while (iarg < param->size()) {
      const string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (parseSpecifiedZones(zone_indices,token,param,iarg)) {
        // no op
      }
      else if (token == "ADD") {
        b_add = true;
      }
      else {
        WUI(WARN,"unrecognized FLAG_TRIS_TOUCHING token " << token << ". Check syntax and try again.");
        return;
      }
    }

    if (zone_indices.empty()) {
      WUI(WARN,"no zones or subzones selected");
      return;
    }

    ss.flagTrisTouchingZones(zone_indices,b_add); // last bool is add

  }

  void processSeparateOverlappingTris(Param * param,const bool b_help) {

    if (b_help) {
      helpSeparateOverlappingTris();
      return;
    }

    if (!ss.hasNonManifoldData()) {
      WUI(WARN,"Please call RUN_DIAGNOSTICS with selected params before calling SEPARATE_OVERLAPPING_TRIS. Returning.");
      helpRunDiagnostics();
      return;
    }

    try {

      const double length_scale = ss.getRmax();
      double dn = 1.0E-8*length_scale;

      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (token == "DN") {
          dn = param->getDouble(iarg++);
        }
        else {
          WUI(WARN,"unrecognized SEPARATE_OVERLAPPING_TRIS token " << token << ". skipping");
        }
      }

      ss.separateOverlappingTris(dn);

      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem during SEPARATE_OVERLAPPING_TRIS");
      helpSeparateOverlappingTris();
    }

  }

  void helpSeparateOverlappingTris() const {
    WUI(INFO,"moves overlapping tris in normal direction to separate them\n" <<
        "  examples:\n" <<
        "    SEPARATE_OVERLAPPING_TRIS\n" <<
        "    SEPARATE_OVERLAPPING_TRIS DN 1.0E-12 # defaults\n" <<
        "  more info: [$CWIKB:surfer_repair]"
        );
  }

  void processQueryHelp() {
    WUI(INFO,"QUERY provides geometric info on zones or subzones\n" <<
	"  examples:\n" <<
	"    QUERY SUBZONE 3\n"
	"    QUERY ZONE a,b,c"
	);
  }

  void processQuery(Param * param,const bool help) {
    if (help) {
      processQueryHelp();
      return;
    }
    bool b_zones = false;
    bool b_subzones = false;
    vector<int> subzone_indices;
    int iarg = 0;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if ((token == "ZONE")||(token == "ZONES")) {
	// should be another token that is a comma-delimited list of zone names and/or numbers...
	if (iarg >= param->size()) {
	  WUI(WARN,"QUERY ZONE expects a comma-delimited list of zone names and/or numbers");
	  processQueryHelp();
	  return;
	}
	else {
	  for (int izn = 0; izn < ss.zoneVec.size(); ++izn)
	    ss.zoneVec[izn].flag = 0;
	  vector<string> zoneNameVec;
	  MiscUtils::splitCsv(zoneNameVec,param->getString(iarg++));
	  int ierr = 0;
	  for (int ii = 0; ii < zoneNameVec.size(); ++ii) {
	    // allow for a wildcard...
	    bool found = false;
	    for (int izn = 0; izn < ss.zoneVec.size(); ++izn) {
	      if (MiscUtils::strcmp_wildcard(ss.zoneVec[izn].getName(),zoneNameVec[ii])) {
		ss.zoneVec[izn].flag = 1;
		found = true;
	      }
	    }
	    if (!found) {
	      // maybe it was a zone index...
	      int izn;
	      if (from_string<int>(izn,zoneNameVec[ii],std::dec)&&(izn >= 0)&&(izn < ss.zoneVec.size())) {
		ss.zoneVec[izn].flag = 1;
	      }
	      else {
		WUI(WARN,"no matching ZONE found for " << zoneNameVec[ii]);
		ierr = -1;
	      }
	    }
	  }
	  if (ierr != 0) {
	    processQueryHelp();
	    return;
	  }
	}
	b_zones = true;
      }
      else if ((token == "SUBZONE")||(token == "SUBZONES")) {
        // expect a comma-delimited set of integers...
        b_subzones = true;
        const string subzonesCsv = param->getString(iarg++);
        MiscUtils::splitCsv(subzone_indices,subzonesCsv);
      }
      else {
	WUI(WARN,"unrecognized QUERY token " << token);
	processQueryHelp();
	return;
      }
    }
    if (b_zones) {
      double area_sum = 0.0;
      int count = 0;
      for (int izn = 0; izn < ss.zoneVec.size(); ++izn) {
	if (ss.zoneVec[izn].flag == 1) {
	  double area = ss.zoneArea(izn);
	  WUI(INFO,"zone " << izn << " " << ss.zoneVec[izn].getName() << " has area " << area);
	  area_sum += area;
	  ++count;
	}
      }
      if (count > 1) {
	WUI(INFO,"total area all " << count << " zones: " << area_sum);
      }
    }
    else if (b_subzones) {
      ss.sz_flag.resize(ss.nsz);
      ss.sz_flag.setAll(0);
      for (vector<int>::iterator it=subzone_indices.begin(); it!=subzone_indices.end(); ++it) {
	assert((*it >= 0)&&(*it < ss.nsz));
	ss.sz_flag[*it] = 1;
      }
      ss.queryFlaggedSubzones();
    }
  }

  void processFlattenRotor(Param * param,const bool help) {
    if (help) {
      WUI(INFO,"FLATTEN_ROTOR writes all points transformed in x-r for disk building in moving mesh simulations.\n" <<
	  "  examples:\n" <<
	  "    FLATTEN_ROTOR X0 6.356035 0 0.15643 X1 6.437815 0 0.15643"
          );
      return;
    }

    bool b_x0 = false;
    double x0[3];
    bool b_x1 = false;
    double x1[3];

    int iarg = 0;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "X0") {
	FOR_I3 x0[i] = param->getDouble(iarg++);
	b_x0 = true;
      }
      else if (token == "X1") {
	FOR_I3 x1[i] = param->getDouble(iarg++);
	b_x1 = true;
      }
      else {
	WUI(WARN,"unrecognized FLATTEN_ROTOR token " << token << ". skipping");
      }
    }

    if ((!b_x0)||(!b_x1)) {
      WUI(WARN,"FLATTEN_ROTOR missing X0 or X1");
      return;
    }

    const double dx[3] = DIFF(x1,x0);
    if ((dx[1] == 0)&&(dx[2] == 0)) {
      // x-aligned...
      FILE * fp = fopen("flatten_rotor.dat","w");
      for (int isp = 0; isp < ss.nsp; ++isp) {
	if ((ss.xsp[isp][0] > x0[0])&&(ss.xsp[isp][0] < x1[0])) {
	  const double dy = ss.xsp[isp][1] - x0[1];
	  const double dz = ss.xsp[isp][2] - x0[2];
	  const double r = sqrt(dy*dy+dz*dz);
	  fprintf(fp,"%18.15e %18.15e\n",ss.xsp[isp][0],r);
	}
      }
      fclose(fp);
    }
    else {
      assert(0);
    }

    cout << "take a look at flatten_rotor.dat" << endl;

  }

  void processFlip(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpFlip();
      return;
    }

    try {
      // expects a comma delimited list of subzone indices
      vector<int> subzone_indices;
      vector<int> group_indices;

      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));

        if (parseSpecifiedSubzones(subzone_indices,token,param,iarg)) {
          //no op; parsed info
        }
        else if (token == "GROUPS") {
          const string groupsCsv = param->getString(iarg++);
          vector<string> groupsVec;
          MiscUtils::splitCsv(groupsVec,groupsCsv);
          for (int i = 0, limit = groupsVec.size(); i < limit; ++i) {
            const int group_index = atoi(groupsVec[i].c_str());
            group_indices.push_back(group_index);
          }
        }
      }

      if ( !group_indices.empty() && !subzone_indices.empty() ) {
        WUI(WARN,"cannot specify both subzones and groups for flipping; skipping");
        return;
      }

      if ( !group_indices.empty() )  ss.flipSelectedGroups(group_indices);
      else if ( !subzone_indices.empty() ) ss.flipSelectedSubzones(subzone_indices);
      else ss.flipSurface();

      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem during flip");
      helpFlip();
    }
  }
  void helpFlip() {
    WUI(INFO,"inverts the normals on a spefied portion of the surface\n" <<
        "  examples:\n" <<
        "    FLIP (entire surface)\n" <<
        "    FLIP ZONE_NAMES top,bottom (portion of surface)\n" <<
        "    FLIP GROUPS 3,4,7 (disjoint parts of surface)\n" <<
        "  more info: [$CWIKB:surfer_repair]"
        );
  }

  void processAlignNormals(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpAlignNormals();
      return;
    }

    enum AlignAlgo {
      SEED,
      AUTO,
      BLOAT,
      BLOAT_SEED,
      NONE
    };

    try {
      ss.ensureStoszSzosz();
      AlignAlgo align_algo = NONE;
      int st_seed = 0;
      double delta = 0.0;

      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (token == "SEED") {
          align_algo = (align_algo == AUTO) ? AUTO:SEED;  // persists if AUTO mode, keep it that way
          const int sz_index = param->getInt(iarg++);
          if (sz_index < 0 || sz_index >= ss.nsz) {
            CWARN("invalid subzone index: " << sz_index << "; skipping");
          }
          else {
            COUT2(" > aligning tris to subzone: " << sz_index);
            st_seed = ss.stosz_v[ss.stosz_i[sz_index]];
            assert(st_seed >= 0 && st_seed < ss.nst);
          }
        }
        else if (token == "AUTO") {
          align_algo = AUTO;
        }
        else if (token == "BLOAT") {
          align_algo = BLOAT;
        }
        else if (token == "DELTA") {
          delta = param->getDouble(iarg++);
        }
        else if (token == "BLOAT_SEED") {
          align_algo = BLOAT_SEED;
        }
        else {
          CWARN("unrecognized ALIGN_NORMALS token: " << token << "; skipping");
        }
      }

      if ((delta == 0.0)&&((align_algo == BLOAT)||(align_algo == BLOAT_SEED))) {
        double volume,area;
        ss.getAreaAndVolume(volume,area);
        delta = sqrt(2.0*area/double(ss.nst)); // use avg area if unspecified
      }

      switch (align_algo) {
      case SEED:
        ss.alignNormals(st_seed);
        break;
      case AUTO:
        ss.alignNormals(st_seed);  // orients normals by group first, then volume based flipping
        ss.alignNormalsAuto();
        break;
      case BLOAT:
        ss.alignNormalsUsingBloatedSurface(delta,false);
        break;
      case BLOAT_SEED:
        ss.alignNormalsUsingBloatedSurface(delta,true);
        break;
      default:
        ss.alignNormals(st_seed);
      }

      WebUI::webUIOutput.ensureImage();
    }
    catch (int e) {
      WUI(WARN,"problem during align normals");
      helpAlignNormals();
    }
  }
  void helpAlignNormals() {
    WUI(INFO,"automatically marches out and aligns surface normals to the specified subzone tri\n" <<
        "  examples:\n" <<
        "    ALIGN_NORMALS SEED 23 \n" <<
        "  more info: [$CWIKB:surfer_repair]"
        );
  }

  void processRenameZone(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpRenameZone();
      return;
    }

    try {
      int iarg = 0;
      int izone = -1;
      string oldname = param->getString(iarg++);
      if (MiscUtils::toUpperCase(oldname) == "INDEX") {
        izone = param->getInt(iarg++);

        if (izone < 0 || izone >= int(ss.zoneVec.size())) {
          CWARN("out-of-bounds INDEX " << izone << "; skipping rename");
          return;
        }
        oldname = ss.zoneVec[izone].getName();
      }


      string newname = param->getString(iarg++);
      while (iarg < param->size()) newname += " " + param->getString(iarg++);

      if (izone == -1) ss.renameZone(oldname,newname);
      else ss.renameZone(izone,newname);

      WebUI::webUIOutput.ensureMenu();
    }
    catch (int e) {
      WUI(WARN,"problem renaming zone");
      helpRenameZone();
    }
  }
  void helpRenameZone() {
    WUI(INFO,"renames an existing zone\n" <<
        "  examples:\n" <<
        "    RENAME_ZONE top very_top (by name) \n" <<
        "    RENAME_ZONE INDEX 13 very_top (by zone index)"
        );
  }

  void processSetZoneBC(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"sets the boundary condition for a particilar zone; only used by the App to fill in boundary condition template information"
          );
      return;
    }

    int iarg = 0;
    string bc_name;
    vector<int> zone_indices;

    //const int nzn = ss.zoneVec.size();

    while (iarg < param->size()) {
      string token = param->getString(iarg++);

      if (parseSpecifiedZones(zone_indices,token,param,iarg)) {
        //no op; parsed info
      }
      else if (token == "BC") {
        bc_name = param->getString(iarg++);
      }
      else {
        CWARN("unrecognized token \"" << token << "\"; skipping");
      }
    }

    if (bc_name.length() == 0) {
      WUI(WARN,"no destination boundary condition was presecribed; skipping");
      return;
    }

    if (!zone_indices.empty()) {
      for (vector<int>::iterator it=zone_indices.begin(); it!=zone_indices.end(); ++it) {
        ss.zoneVec[*it].setBC(bc_name);
      }
    }
  }

  void processSubzonesToIndexedZones(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    try {

      vector<int> zone_indices;

      int iarg = 0;
      while (iarg < param->size()) {
        const string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (parseSpecifiedZones(zone_indices,token,param,iarg)) {
          // no op
        }
        else {
          WUI(WARN,"unrecognized SUBZONES_TO_ZONES token " << token << "; skipping");
        }
      }
      ss.migrateSubzonesToZones(zone_indices);

      WebUI::webUIOutput.ensureImage();
      WebUI::webUIOutput.ensureMenu();
    }
    catch (int e) {
      WUI(WARN,"problem making each subzone a zone");
    }
  }

  void processMoveToZone(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      helpMoveToZone();
      return;
    }

    try {
      int iarg = 0;

      //bool b_all = false;
      int izone = -1;
      string zonename = "new_zone";
      vector<int> subzone_indices;

      while (iarg < param->size()) {
        string token = MiscUtils::toUpperCase(param->getString(iarg++));

        if (parseSpecifiedSubzones(subzone_indices,token,param,iarg)) {
          //no op; parsed info
        }
        else if (token == "NAME") {
          zonename = param->getString(iarg++);
        }
        else if (token == "INDEX") {
          izone = param->getInt(iarg++);
          if (izone < 0 || izone >= int(ss.zoneVec.size())) {
            CWARN("invalid zone index " << izone << " specified");
            izone = -1;
          }
        }
        else {
          CWARN("unrecognized token \"" << token << "\"; skipping");
        }
      }

      if (izone == -1) {
        izone = ss.getZoneIndex(zonename);
        if (izone == -1) {
          izone = ss.addNewZone(zonename);
        }
      }
      assert(izone != -1);

      assert(ss.nsz == ss.szozn_i[ss.zoneVec.size()]);
      ss.sz_flag.setLength(ss.nsz);
      if (!subzone_indices.empty()) {
        ss.sz_flag.setAll(0);
        for (vector<int>::iterator it=subzone_indices.begin(); it!=subzone_indices.end(); ++it) {
          const int isz = *it;
          if (isz < 0 || isz >= ss.nsz) {
            CWARN("subzone index out-of-bounds: " << isz);
          }
          else {
            ss.sz_flag[isz] = 1;
          }
        }
      }
      ss.moveFlaggedSubzonesToZone(izone);

      WebUI::webUIOutput.ensureImage();
      WebUI::webUIOutput.ensureMenu();
    }
    catch (int e) {
      WUI(WARN,"problem moving subzones to zone");
      helpMoveToZone();
    }
  }
  void helpMoveToZone() {
    WUI(INFO,"moves the specified portions of a surface to a new or existing zone\n" <<
        "  examples:\n" <<
        "    MOVE_TO_ZONE NAME flaps ZONE_NAMES flap_a,flap_b,flap_c"
        );
  }

  void processSliceGen(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"generates image slices along an axis between specified subzones"
          );
      return;
    }

    string token;

    int iarg = 0;

    token = param->getString(iarg++);
    assert(token == "START");
    const int start = param->getInt(iarg++);

    token = param->getString(iarg++);
    assert(token == "END");
    const int end = param->getInt(iarg++);

    token = param->getString(iarg++);
    assert(token == "N");
    const int n = param->getInt(iarg++);

    // start and end refer to zones in the subzoneDataVec...

    assert( (start >= 0) && (start < int(ss.subzoneDataVec.size())) );
    assert( (end >= 0) && (end < int(ss.subzoneDataVec.size())) );

    // check if both start and end are nearly normal...

    double x_start[3];
    FOR_I3 x_start[i] = ss.subzoneDataVec[start].xc[i];
    const double n_start_mag = MAG(ss.subzoneDataVec[start].normal);
    if (n_start_mag < 0.95*ss.subzoneDataVec[start].area)
      cout << "Warning: START may not be planar" << endl;
    double n_start[3];
    FOR_I3 n_start[i] = ss.subzoneDataVec[start].normal[i]/n_start_mag;

    double x_end[3];
    FOR_I3 x_end[i] = ss.subzoneDataVec[end].xc[i];
    const double n_end_mag = MAG(ss.subzoneDataVec[end].normal);
    if (n_end_mag < 0.95*ss.subzoneDataVec[end].area)
      cout << "Warning: END may not be planar" << endl;
    double n_end[3];
    FOR_I3 n_end[i] = ss.subzoneDataVec[end].normal[i]/n_end_mag;

    // now build the bezier curve between the points with the specified end-point normals...

    const double dx[3] = DIFF(x_end,x_start);

    // check normal alignment -- should be pointing out...
    assert(DOT_PRODUCT(n_start,dx) < 0.0);
    assert(DOT_PRODUCT(n_end,dx) > 0.0);

    const double dn[3] = DIFF(n_end,n_start);
    const double alpha = 0.5*DOT_PRODUCT(dx,dn)/(2.0-DOT_PRODUCT(n_start,n_end));

    // {
    //   FILE * fp = fopen("pts.dat","w");
    //   fprintf(fp,"%lf %lf %lf\n",x_start[0],x_start[1],x_start[2]);
    //   fprintf(fp,"%lf %lf %lf\n",x_start[0]-alpha*n_start[0],x_start[1]-alpha*n_start[1],x_start[2]-alpha*n_start[2]);
    //   fprintf(fp,"%lf %lf %lf\n",x_end[0]-alpha*n_end[0],x_end[1]-alpha*n_end[1],x_end[2]-alpha*n_end[2]);
    //   fprintf(fp,"%lf %lf %lf\n",x_end[0],x_end[1],x_end[2]);
    //   fclose(fp);
    // }

    // write out points...

    FILE * fp = fopen("killfile.slice","w");
    for (int ip = 0; ip <= n; ++ip) {
      const double t = max(1.0E-5,min(1.0-1.0E-5,double(ip)/double(n))); // modify first and last planes slightly
      const double omt = 1.0-t;
      double xp[3]; FOR_J3 xp[j] = omt*omt*omt*x_start[j] + 3.0*omt*omt*t*(x_start[j]-alpha*n_start[j]) + 3.0*omt*t*t*(x_end[j]-alpha*n_end[j]) + t*t*t*x_end[j];
      // n is dxdt...
      double np[3]; FOR_J3 np[j] = -3.0*omt*omt*alpha*n_start[j] + 6.0*omt*t*(dx[j]-alpha*dn[j]) + 3.0*t*t*alpha*n_end[j];
      const double np_mag = MAG(np); assert(np_mag > 0.0);
      FOR_I3 np[i] /= np_mag;

      fprintf(fp,"WRITE_IMAGE NAME images/un/ax%03d SIZE 2000 2000 WIDTH 0.6 GEOM PLANE %lf %lf %lf %lf %lf %lf VAR (%lf*comp(u_avg,0))+(%lf*comp(u_avg,1))+(%lf*comp(u_avg,2)) RANGE -30.48 304.8\n",
	      ip,xp[0],xp[1],xp[2],np[0],np[1],np[2],np[0],np[1],np[2]);

      fprintf(fp,"WRITE_IMAGE NAME images/Tavg/ax%03d SIZE 2000 2000 WIDTH 0.6 GEOM PLANE %lf %lf %lf %lf %lf %lf VAR T_avg RANGE 644.26 2144.26\n",
	      ip,xp[0],xp[1],xp[2],np[0],np[1],np[2]);

      fprintf(fp,"WRITE_IMAGE NAME images/Trms/ax%03d SIZE 2000 2000 WIDTH 0.6 GEOM PLANE %lf %lf %lf %lf %lf %lf VAR T_rms RANGE 0 444.44\n",
	      ip,xp[0],xp[1],xp[2],np[0],np[1],np[2]);

      fprintf(fp,"WRITE_IMAGE NAME images/co/ax%03d SIZE 2000 2000 WIDTH 0.6 GEOM PLANE %lf %lf %lf %lf %lf %lf VAR log10(X_CO_DRY_avg) RANGE -10 0\n",
	      ip,xp[0],xp[1],xp[2],np[0],np[1],np[2]);

      fprintf(fp,"WRITE_IMAGE NAME images/nox/ax%03d SIZE 2000 2000 WIDTH 0.6 GEOM PLANE %lf %lf %lf %lf %lf %lf VAR log10(X_NOX_DRY15_avg) RANGE -10 -4\n",
	      ip,xp[0],xp[1],xp[2],np[0],np[1],np[2]);

    }
    fclose(fp);
    cout << "just wrote killfile.slice" << endl;

    getchar();

  }

  void processWriteFluxProbesCylinder(Param * param,const bool help) {
    if (help) {
      WUI(INFO,"debugging utility"
          );
      return;
    }

    bool b_name = false;
    string name;

    int isz = -1;
    int isz_prev = -1;
    int isz_next = -1;

    bool b_write_image = false;
    double up[3];

    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "NAME") {
        b_name = true;
        name = param->getString(iarg++);
      }
      else if (token == "SUBZONE") {
        isz = param->getInt(iarg++);
      }
      else if (token == "SUBZONE_PREV") {
        isz_prev = param->getInt(iarg++);
      }
      else if (token == "SUBZONE_NEXT") {
        isz_next = param->getInt(iarg++);
      }
      else if (token == "WRITE_IMAGE") {
        b_write_image = true;
        token = param->getString(iarg++);
        assert(token == "UP");
        FOR_I3 up[i] = param->getDouble(iarg++);
      }
      else {
        cout << "Error: unrecognized processWriteFluxProbesCylinder param: " << token << endl;
        return;
      }
    }

    assert(b_name);
    assert(isz != -1);
    assert(isz_prev != -1);
    assert(isz_next != -1);

    ss.writeFluxProbesCylinder(name,isz,isz_prev,isz_next,b_write_image,up);

  }

  void processMultiFluxProbeHelp() {
    WUI(INFO,"MULTI_FLUX_PROBE provides help describing multi-flux-probe commands for the solvers,\n" <<
	"which use zones to indicate probe locations and direction\n" <<
	"  examples:\n" <<
	"    MULTI_FLUX_PROBE ZONE liner_holes"
	);
  }

  void processMultiFluxProbe(Param * param,const bool help) {

    if (help) {
      processMultiFluxProbeHelp();
      return;
    }

    bool b_zones = false;
    bool b_error = false;

    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if ((token == "ZONE")||(token == "ZONES")) {
	// should be another token that is a comma-delimited list of zone names and/or numbers...
	if (iarg >= param->size()) {
	  WUI(WARN,"MULTI_FLUX_PROBE ZONE expects a comma-delimited list of zone names and/or numbers");
	  processMultiFluxProbeHelp();
	  return;
	}
	else {
	  for (int izn = 0; izn < ss.zoneVec.size(); ++izn)
	    ss.zoneVec[izn].flag = 0;
	  vector<string> zoneNameVec;
	  MiscUtils::splitCsv(zoneNameVec,param->getString(iarg++));
	  int ierr = 0;
	  for (int ii = 0; ii < zoneNameVec.size(); ++ii) {
	    // allow for a wildcard...
	    bool found = false;
	    for (int izn = 0; izn < ss.zoneVec.size(); ++izn) {
	      if (MiscUtils::strcmp_wildcard(ss.zoneVec[izn].getName(),zoneNameVec[ii])) {
		ss.zoneVec[izn].flag = 1;
		found = true;
	      }
	    }
	    if (!found) {
	      // maybe it was a zone index...
	      int izn;
	      if (from_string<int>(izn,zoneNameVec[ii],std::dec)&&(izn >= 0)&&(izn < ss.zoneVec.size())) {
		ss.zoneVec[izn].flag = 1;
	      }
	      else {
		WUI(WARN,"no matching ZONE found for " << zoneNameVec[ii]);
		ierr = -1;
	      }
	    }
	  }
	  if (ierr != 0) {
	    b_error = true;
	  }
	}
	b_zones = true;
      }
      else {
	WUI(WARN,"unrecognized MULTI_FLUX_PROBE token: " << token);
	b_error = true;
      }
    }

    if (b_error) {
      WUI(INFO,"Fix syntax problems and rerun");
      return;
    }
    if (!b_zones) {
      WUI(WARN,"missing ZONE or ZONES param");
      return;
    }

    ss.suggestMultiFluxProbeFlaggedZones();

  }

  void processCalcClosestCylinder(Param * param,const bool help) {
    if (!hasSurface()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"calculates closest cylinder"
          );
      return;
    }

    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if ((token == "ZONES")||(token == "ZONE")) {
        // select zones based on zone name provided in comma-delimited list
        vector<string> zoneNameVec;
        MiscUtils::splitCsv(zoneNameVec,param->getString(iarg++));
        ss.zone_flag.setLength(ss.zoneVec.size());
        ss.zone_flag.setAll(0);
        // would be better to use a zone name map structure
        for (int ii = 0, ii_end=zoneNameVec.size(); ii < ii_end; ++ii) {
          string zoneName = zoneNameVec[ii];
          bool b_zoneFound = false;
          for (int izone = 0,nzn=ss.zoneVec.size(); izone < nzn; izone++) {
            if (zoneName == ss.zoneVec[izone].getName()) {
              cout << " > ZONE: " << izone << ", selecting st's in zone: " << ss.zoneVec[izone].getName() << endl;
              ss.zone_flag[izone] = 1;
              b_zoneFound = true;
              break;
            }
          }
          if (not b_zoneFound) {
            cout << "Warning: CALC_CLOSEST_CYL ZONE \"" << zoneName << "\" not found! " << endl;
          }
        }
        ss.st_flag.setLength(ss.nst);
        ss.st_flag.setAll(0);
        for (int ist = 0; ist < ss.nst; ++ist) {
          const int izone = ss.znost[ist];
          if (ss.zone_flag[izone] == 1)
            ss.st_flag[ist] = 1;
        }
      }
      else if (token == "ZONE_IDS") {
        // select zones based on zone index provided in comma-delimited list
        vector<int> zoneIdVec;
        MiscUtils::splitCsv(zoneIdVec,param->getString(iarg++));
        ss.zone_flag.setLength(ss.zoneVec.size());
        ss.zone_flag.setAll(0);
        const int nzn = ss.zoneVec.size();
        for (int ii = 0, ii_end=zoneIdVec.size(); ii < ii_end; ++ii) {
          const int izone = zoneIdVec[ii];
          if ((izone < 0)||(izone >= nzn)) {
            cout << "Warning: CALC_CLOSEST_CYL ZONE_IDS out of range: " << izone << endl;
          }
          else {
            cout << " > ZONE_IDS: " << izone << ", selecting st's in zone: " << ss.zoneVec[izone].getName() << endl;
            ss.zone_flag[izone] = 1;
          }
        }
        ss.st_flag.setLength(ss.nst);
        ss.st_flag.setAll(0);
        for (int ist = 0; ist < ss.nst; ++ist) {
          const int izone = ss.znost[ist];
          if (ss.zone_flag[izone] == 1)
            ss.st_flag[ist] = 1;
        }
      }
      else if (token == "SUBZONE_IDS") {
        vector<int> subZoneIdVec;
        MiscUtils::splitCsv(subZoneIdVec,param->getString(iarg++));
        if (!subZoneIdVec.empty()) ss.flagTrisFromSubzoneVec(subZoneIdVec);  // flag tris appropriately
      }
    }
    double normal[3];
    double xc[3];
    ss.calcClosestCylinderForFlaggedTris(normal, xc);
    cout << "\n=====================================================\n"
         << " > Closest Cylinder Normal: " << COUT_VEC(normal) << endl
         << " > Area Averaged Position: " << COUT_VEC(xc) << endl
         << "=====================================================\n\n";
  }

  void writeAutoCompleteJson(FILE * fp) {
    fprintf(fp,",\n\"siData\": {");
    fprintf(fp,"\n  \"kbVersion\": \"%s\"",CTI::cti_docs_version);
    const string siData =
#include "surfer.siData"
      ;
    fprintf(fp,",\n  %s\n",siData.c_str());
    fprintf(fp,"}");
  }

  // these parsing helper routines moved to ss, because they are used there too...

  bool parseSpecifiedZones(vector<int>& zone_indices,const string& token,Param * param,int& iarg) {
    return ss.parseSpecifiedZones(zone_indices,token,param,iarg);
  }

  bool parseSpecifiedSubzones(vector<int>& subzone_indices,const string& token,Param * param,int& iarg) {
    return ss.parseSpecifiedSubzones(subzone_indices,token,param,iarg);
  }

  bool parseSpecifiedDynamicEdges(vector<int>& zone_indices,const string& token,Param * param,int& iarg) {
    ss.ensureDynamicEdgeGroups();
    bool parsed=ss.parseSpecifiedEdges(zone_indices,token,param,iarg);

    if (!parsed) {
      return false;
    }
    else {
      // if none specified, interpret as "ALL"
      if (zone_indices.empty()) {
        for (int izone=0,lim=ss.eoiGroupDataVec.size(); izone<lim; ++izone) zone_indices.push_back(izone+DYN_E_SZ_MIN);
      }

      vector<int>::iterator it=zone_indices.begin();
      while (it!=zone_indices.end()) {
        if ((*it < DYN_E_SZ_MIN) || (*it > DYN_E_SZ_MAX)) {
          COUT1(" > pruned dynamic-edge zone "<< *it << " because it is not a dynamic edge");
          it = zone_indices.erase(it);
        }
        else {
          *it -= DYN_E_SZ_MIN;
          ++it;
        }
      }
      return true;
    }
  }

  bool parseSpecifiedOpenEdges(vector<int>& zone_indices,const string& token,Param * param,int& iarg) {
    ss.ensureOpenEdgeGroups();
    bool parsed=ss.parseSpecifiedEdges(zone_indices,token,param,iarg);

    if (!parsed) {
      return false;
    }
    else {
      // if none specified, interpret as "ALL"
      if (zone_indices.empty()) {
        for (int izone=0,lim=ss.getNOpenEdgeGroups(); izone<lim; ++izone) zone_indices.push_back(izone+OPEN_E_SZ_MIN);
      }

      vector<int>::iterator it=zone_indices.begin();
      while (it!=zone_indices.end()) {
        if ((*it < OPEN_E_SZ_MIN) || (*it > OPEN_E_SZ_MAX)) {
          COUT1(" > pruned open-edge zone "<< *it<< " because it is not an open edge");
          it = zone_indices.erase(it);
        }
        else {
          *it -= OPEN_E_SZ_MIN;
          ++it;
        }
      }
      return true;
    }
  }

  void dumpVol(double x0[3],double x1[3],double x2[3],double x3[3]) {

    cout << "dumpVol: " << endl;
    cout << SIGNED_TET_VOLUME_6(x0,x1,x2,x3) << endl;
    cout << SIGNED_TET_VOLUME_6(x0,x2,x3,x1) << endl;
    cout << SIGNED_TET_VOLUME_6(x0,x3,x1,x2) << endl;
    cout << SIGNED_TET_VOLUME_6(x1,x0,x3,x2) << endl;
    cout << SIGNED_TET_VOLUME_6(x1,x2,x0,x3) << endl;
    cout << SIGNED_TET_VOLUME_6(x1,x3,x2,x0) << endl;
    cout << SIGNED_TET_VOLUME_6(x2,x0,x1,x3) << endl;
    cout << SIGNED_TET_VOLUME_6(x2,x1,x3,x0) << endl;
    cout << SIGNED_TET_VOLUME_6(x2,x3,x0,x1) << endl;
    cout << SIGNED_TET_VOLUME_6(x3,x0,x2,x1) << endl;
    cout << SIGNED_TET_VOLUME_6(x3,x1,x0,x2) << endl;
    cout << SIGNED_TET_VOLUME_6(x3,x2,x1,x0) << endl;

  }

  void checkTetVolExact() {

    cout << "checkTetVolExact()" << endl;

    double x[4][3];

    double vol[24];
    bool pmatch[24][24];
    bool mmatch[24][24];

    for (int ii = 0; ii < 24; ++ii) {
      for (int jj = 0; jj < 24; ++jj) {
        pmatch[ii][jj] = true;
        mmatch[ii][jj] = true;
      }
    }


    int iter = 0;
    while (1) {

      ++iter;
      if (iter%1000 == 0)
        cout << " > done " << iter << endl;

      FOR_I3 x[0][i] = double(rand())/double(RAND_MAX)-0.5;
      FOR_I3 x[1][i] = double(rand())/double(RAND_MAX)-0.5;
      FOR_I3 x[2][i] = double(rand())/double(RAND_MAX)-0.5;
      FOR_I3 x[3][i] = double(rand())/double(RAND_MAX)-0.5;

      int ii = 0;
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          if (j != i) {
            for (int k = 0; k < 4; ++k) {
              if ((k != i)&&(k != j)) {
                for (int l = 0; l < 4; ++l) {
                  if ((l != i)&&(l != j)&&(l != k)) {
                    vol[ii] = SIGNED_TET_VOLUME_6(x[i],x[j],x[k],x[l]);
                    cout << ii << " ijkl: " << i << j << k << l << " vol: " << vol[ii] << endl;
                    ++ii;
                  }
                }
              }
            }
          }
        }
      }
      assert(ii == 24);

      dumpVol(x[0],x[1],x[2],x[3]);

      for (int ii = 0; ii < 24; ++ii) {
        for (int jj = 0; jj < 24; ++jj) {
          if (pmatch[ii][jj]) {
            if (vol[ii] != vol[jj]) {
              pmatch[ii][jj] = false;
              cout << "reseting pmatch" << endl;
            }
          }
          if (mmatch[ii][jj]) {
            if (vol[ii] != -vol[jj]) {
              mmatch[ii][jj] = false;
              cout << "reseting mmatch" << endl;
            }
          }
        }
      }

      cout << "current match: " << endl;

      for (int ii = 0; ii < 24; ++ii) {
        cout << "ii " << ii << " matches: ";
        for (int jj = 0; jj < 24; ++jj) {
          if (ii == jj) {
            assert(pmatch[ii][jj]);
            assert(!mmatch[ii][jj]);
          }
          else if (pmatch[ii][jj]) {
            cout << jj << " ";
          }
          else if (mmatch[ii][jj]) {
            cout << "-" << jj << " ";
          }
        }
        cout << endl;
      }

      getchar();

    }

  }

  bool processDebugTriTri(Param * param,const bool help) {

    cout << "\n================================================================================\n" <<
      "WORKING ON TRITRI_DEBUG..." << endl;

    bool b_hold = false;
    int iarg = 0;
    if (param->getString(iarg) == "HOLD") {
      b_hold = true;
      ++iarg;
    }

    while (iarg < param->size()) {

      cout << " > " << param->getString(iarg) << endl;
      double x0[3][3];
      double x1[3][3];
      double * x0p[3]; FOR_I3 x0p[i] = x0[i];
      double * x1p[3]; FOR_I3 x1p[i] = x1[i];
      MiscUtils::readTriTriBin(x0p,x1p,param->getString(iarg));
      ++iarg;

      int idata[6];
      double ddata[4];

      //int ierr = MiscUtils::calcTriTriIntersection(idata,ddata,x0p,x1p);
      //if (ierr == -2) {

      // this code checks the uniqueness of hte tet volume calc -- end result:
      // there are 12 unique volumes, as illustrated in dumpVol.
      //checkTetVolExact();

      // Notes:
        //
        // the exact math does NOT appear to comprehend ALL cases and this can lead to
        // some unrealizable geometries. I think we may have to do all 12 unique
        // tet volume calculations per tet (as in dumpVol). For example, uncomment
        // thcode below and run this routine on the ./debug_tritri/tritri.nearly_coplanar_edge.bin
        // binary triangle pair that is part of this source...
        //
        // ./surfer.exe -i junk --DEBUG_TRI_TRI ./debug_tritri/tritri.nearly_coplanar_edge.bin
        //
        // (the "-i junk" simply ensures you do not read the surfer.in file in this directory)

        /*
          cout << "looking at vol: " << endl;
          FOR_I3 {
          cout << "\ntri 0, point " << i << endl;
          // the volume of the tet from a particular x0 to tri x1...
          dumpVol(x0p[i],x1p[0],x1p[1],x1p[2]);
          }
          FOR_I3 {
          cout << "\ntri 1, point " << i << endl;
          // the volume of the tet from a particular x0 to tri x1...
          dumpVol(x1p[i],x0p[0],x0p[1],x0p[2]);
          }
          getchar();
        */

        // run the suite of permutations...
        FOR_I3 FOR_J3 {
          x0p[0] = x0[i];
          x0p[1] = x0[(i+1)%3];
          x0p[2] = x0[(i+2)%3];
          x1p[0] = x1[j];
          x1p[1] = x1[(j+1)%3];
          x1p[2] = x1[(j+2)%3];
          MiscUtils::calcTriTriIntersection(idata,ddata,x0p,x1p);
          if (b_hold) getchar();
          x1p[0] = x0[i];
          x1p[1] = x0[(i+1)%3];
          x1p[2] = x0[(i+2)%3];
          x0p[0] = x1[j];
          x0p[1] = x1[(j+1)%3];
          x0p[2] = x1[(j+2)%3];
          MiscUtils::calcTriTriIntersection(idata,ddata,x0p,x1p);
          if (b_hold) getchar();
          x0p[0] = x0[(i+2)%3];
          x0p[1] = x0[(i+1)%3];
          x0p[2] = x0[i];
          x1p[0] = x1[j];
          x1p[1] = x1[(j+1)%3];
          x1p[2] = x1[(j+2)%3];
          MiscUtils::calcTriTriIntersection(idata,ddata,x0p,x1p);
          if (b_hold) getchar();
          x0p[0] = x0[(i+2)%3];
          x0p[1] = x0[(i+1)%3];
          x0p[2] = x0[i];
          x1p[0] = x1[(j+2)%3];
          x1p[1] = x1[(j+1)%3];
          x1p[2] = x1[j];
          MiscUtils::calcTriTriIntersection(idata,ddata,x0p,x1p);
          if (b_hold) getchar();
          x1p[0] = x0[(i+2)%3];
          x1p[1] = x0[(i+1)%3];
          x1p[2] = x0[i];
          x0p[0] = x1[j];
          x0p[1] = x1[(j+1)%3];
          x0p[2] = x1[(j+2)%3];
          MiscUtils::calcTriTriIntersection(idata,ddata,x0p,x1p);
          if (b_hold) getchar();
          x1p[0] = x0[(i+2)%3];
          x1p[1] = x0[(i+1)%3];
          x1p[2] = x0[i];
          x0p[0] = x1[(j+2)%3];
          x0p[1] = x1[(j+1)%3];
          x0p[2] = x1[j];
          MiscUtils::calcTriTriIntersection(idata,ddata,x0p,x1p);
          if (b_hold) getchar();
        }

        /*
      }
      else {
        cout << " > success: returned n: " << ierr << endl;
      }
        */

    }

    throw(0);

  }

  void processTestCht(Param * param,const bool help) {

    if (help) {
      WUI(INFO,"TEST_CHT tests which surface zones connected to the solid in the tet_bin file. Example:\n" <<
	  "  TEST_CHT cht.tet_bin");
      return;
    }

    FILE * fp;
    if ((param->size() == 0)||((fp = fopen(param->getString().c_str(),"rb")) == NULL)) {
      WUI(WARN,"TEST_CHT requires a valid tet_bin file. Example:\n" <<
	  "  TEST_CHT cht.tet_bin");
      return;
    }

    int tet_bin_version = 0;
    int nno,nte;
    fread(&nno,sizeof(int),1,fp); // can be used to indicate the version that includes the tet zone
    if (nno == -1) {
      tet_bin_version = 1;
      fread(&nno,sizeof(int),1,fp);
    }
    assert(nno > 0);
    fread(&nte,sizeof(int),1,fp);

    cout << " > tet_bin version: " << tet_bin_version << " nno: " << nno << " nte: " << nte << endl;

    double (*x_no)[3] = new double[nno][3];
    fread(x_no,sizeof(double),nno*3,fp);

    int (*noote)[4] = new int[nte][4];
    fread(noote,sizeof(int),nte*4,fp);

    int *znote = new int[nte];
    int nzn = 1;
    if (tet_bin_version == 0) {
      for (int ite = 0; ite < nte; ++ite) {
	znote[ite] = 0;
      }
    }
    else {
      assert(tet_bin_version == 1);
      fread(znote,sizeof(int),nte,fp);
      for (int ite = 0; ite < nte; ++ite) {
	nzn = max(nzn,znote[ite]+1);
      }
    }

    fclose(fp);

    double *zone_vol = new double[nzn];
    int *zone_nte = new int[nzn];
    for (int izn = 0; izn < nzn; ++izn) {
      zone_vol[izn] = 0.0;
      zone_nte[izn] = 0;
    }

    double eps = HUGE_VAL;
    for (int ite = 0; ite < nte; ++ite) {
      const int ino0 = noote[ite][0];
      const int ino1 = noote[ite][1];
      const int ino2 = noote[ite][2];
      const int ino3 = noote[ite][3];
      const double vol = SIGNED_TET_VOLUME_6(x_no[ino0],x_no[ino1],x_no[ino2],x_no[ino3])/6.0;
      const int izn = znote[ite];
      zone_vol[izn] += vol;
      zone_nte[izn] += 1;
      // also calc eps from edge lengths...
      for (int i0 = 0; i0 < 3; ++i0) {
	const int ino0 = noote[ite][i0];
	for (int i1 = i0+1; i1 < 4; ++i1) {
	  const int ino1 = noote[ite][i1];
	  const double d2 = DIST2(x_no[ino0],x_no[ino1]);
	  eps = min(eps,d2);
	}
      }
    }
    eps = sqrt(eps);
    cout << " > minimum tet edge length: " << eps << endl;

    for (int izn = 0; izn < nzn; ++izn) {
      cout << " > zone: " << izn << " nte: " << zone_nte[izn] << " volume: " << zone_vol[izn] << endl;
    }

    // put the x_no in an adt...

    Adt<double> * adt = new Adt<double>(nno,x_no,x_no);
    vector<int> intVec;

    std::stringstream exact_ss;
    std::stringstream partial_ss;
    std::stringstream none_ss;

    ss.sp_flag.setLength(ss.nsp);
    for (int isp = 0; isp < ss.nsp; ++isp)
      ss.sp_flag[isp] = -1;
    for (int isz = 0, nsz = ss.zoneVec.size(); isz < nsz; ++isz) {
      for (int ist = 0; ist < ss.nst; ++ist) {
	if (ss.znost[ist] == isz) {
	  FOR_I3 {
	    const int isp = ss.spost[ist][i];
	    ss.sp_flag[isp] = isz;
	  }
	}
      }
      int count = 0;
      int matched = 0;
      double d2_max = 0.0;
      for (int isp = 0; isp < ss.nsp; ++isp) {
	if (ss.sp_flag[isp] == isz) {
	  ++count;
	  assert(intVec.empty());
	  adt->buildListForSphere(intVec,ss.xsp[isp],0.1*eps);
	  if (!intVec.empty()) {
	    assert(intVec.size() == 1);
	    const int ino = intVec[0];
	    const double d2 = DIST2(ss.xsp[isp],x_no[ino]);
	    d2_max = max(d2_max,d2);
	    ++matched;
	    intVec.clear();
	  }
	}
      }
      if (matched == count) {
	exact_ss << " > surface zone " << ss.zoneVec[isz].getName() << " matched " <<  matched << " out of " << count << " nodes, dist: " << sqrt(d2_max) << endl;
      }
      else if (matched > 0) {
	partial_ss << " > surface zone " << ss.zoneVec[isz].getName() << " matched " <<  matched << " out of " << count << " nodes, dist: " << sqrt(d2_max) << endl;
      }
      else {
	none_ss << " > surface zone " << ss.zoneVec[isz].getName() << " matched " <<  matched << " out of " << count << " nodes, dist: " << sqrt(d2_max) << endl;
      }
    }

    cout << "\n===========================================================================\n" <<
      "The following zones exactly matched all nodes" << endl;
    cout << exact_ss.str() << endl;

    cout << "\n===========================================================================\n" <<
      "The following zones partially matched" << endl;
    cout << partial_ss.str() << endl;

    cout << "\n===========================================================================\n" <<
      "The following zones did not match at any nodes" << endl;
    cout << none_ss.str() << endl;

    delete adt;
    delete[] zone_nte;
    delete[] zone_vol;
    delete[] x_no;
    delete[] noote;
    delete[] znote;

  }

  void processThreePoints(Param * param,const bool help) {

    if (help) {
      WUI(INFO,"THREE_POINTS returns the center and normal formed by three points. Example:\n" <<
	  "  THREE_POINTS [X0] <x0> <y0> <z0> [X1] <x1> <y1> <z1> [X2] <x2> <y2> <z2>");
      return;
    }

    int iarg = 0;
    if (param->getString(iarg)[0] == 'X') ++iarg;
    double x0[3];
    FOR_I3 x0[i] = param->getDouble(iarg++);
    if (param->getString(iarg)[0] == 'X') ++iarg;
    double x1[3];
    FOR_I3 x1[i] = param->getDouble(iarg++);
    if (param->getString(iarg)[0] == 'X') ++iarg;
    double x2[3];
    FOR_I3 x2[i] = param->getDouble(iarg++);

    cout << " > x0: " << COUT_VEC(x0) << endl;
    cout << " > x1: " << COUT_VEC(x1) << endl;
    cout << " > x2: " << COUT_VEC(x2) << endl;

    cout << " > center: " <<
      (x0[0]+x1[0]+x2[0])/3.0 << " " <<
      (x0[1]+x1[1]+x2[1])/3.0 << " " <<
      (x0[2]+x1[2]+x2[2])/3.0 << endl;

    const double n2[3] = TRI_NORMAL_2(x0,x1,x2);
    cout << " > normal: " << 0.5*n2[0] << " " << 0.5*n2[1] << " " << 0.5*n2[2] << endl;

  }

  void processInterface(Param * param,const bool help) {

    if (help) {
      WUI(INFO,"INTERFACE builds a interface of revolution around a specified zone. Example:\n" <<
	  "  INTERFACE ZONES <name> MOTION ROTATION XC <x> <y> <z> AXIS <ax> <ay> <az>");
      return;
    }

    bool b_zones = false;
    bool b_motion = false;
    bool b_rotation = false;
    bool b_axis = false;
    double axis[3];
    bool b_xc = false;
    double xc[3];

    int iarg = 0;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if ((token == "ZONE")||(token == "ZONES")) {
	// should be another token that is a comma-delimited list of zone names and/or numbers...
	if (iarg >= param->size()) {
	  WUI(WARN,"INTERFACE ZONE expects a comma-delimited list of zone names and/or numbers");
	  return;
	}
	else {
	  for (int izn = 0; izn < ss.zoneVec.size(); ++izn)
	    ss.zoneVec[izn].flag = 0;
	  vector<string> zoneNameVec;
	  MiscUtils::splitCsv(zoneNameVec,param->getString(iarg++));
	  int ierr = 0;
	  for (int ii = 0; ii < zoneNameVec.size(); ++ii) {
	    // allow for a wildcard...
	    bool found = false;
	    for (int izn = 0; izn < ss.zoneVec.size(); ++izn) {
	      if (MiscUtils::strcmp_wildcard(ss.zoneVec[izn].getName(),zoneNameVec[ii])) {
		ss.zoneVec[izn].flag = 1;
		found = true;
	      }
	    }
	    if (!found) {
	      // maybe it was a zone index...
	      int izn;
	      if (from_string<int>(izn,zoneNameVec[ii],std::dec)&&(izn >= 0)&&(izn < ss.zoneVec.size())) {
		ss.zoneVec[izn].flag = 1;
	      }
	      else {
		WUI(WARN,"no matching ZONE found for " << zoneNameVec[ii]);
		ierr = -1;
	      }
	    }
	  }
	  if (ierr != 0) {
	    return;
	  }
	}
	b_zones = true;
      }
      else if (token == "MOTION") {
	b_motion = true;
      }
      else if (token == "ROTATION") {
	b_rotation = true;
      }
      else if (token == "XC") {
	FOR_I3 xc[i] = param->getDouble(iarg++);
	b_xc = true;
      }
      else if (token == "AXIS") {
	FOR_I3 axis[i] = param->getDouble(iarg++);
	b_axis = true;
      }
      else {
	WUI(WARN,"unrecognized INTERFACE token \"" << token << "\"");
	return;
      }
    }

    assert(b_zones);
    assert(b_motion);
    assert(b_rotation);
    assert(b_xc);
    assert(b_axis);

    ss.interfaceFlaggedZones(xc,axis);
    WebUI::webUIOutput.ensureImage();

  }

  void processMakeFF(Param * param,const bool help) {

    if (help) {
      WUI(INFO,"MAKE_FF builds a farfield surface around a moving part. Example:\n" <<
	  "  MAKE_FF ZONES <name> MOTION ROTATION XC <x> <y> <z> AXIS <ax> <ay> <az>");
      return;
    }
    
    bool b_zones = false;
    bool b_motion = false;
    bool b_rotation = false;
    bool b_axis = false;
    double axis[3];
    bool b_xc = false;
    double xc[3];
    bool b_dn = false;
    double dn;
    bool b_ds = false;
    double ds;
    
    int iarg = 0;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if ((token == "ZONE")||(token == "ZONES")) {
	// should be another token that is a comma-delimited list of zone names and/or numbers...
	if (iarg >= param->size()) {
	  WUI(WARN,"MAKE_FF ZONE expects a comma-delimited list of zone names and/or numbers");
	  return;
	}
	else {
	  for (int izn = 0; izn < ss.zoneVec.size(); ++izn)
	    ss.zoneVec[izn].flag = 0;
	  vector<string> zoneNameVec;
	  MiscUtils::splitCsv(zoneNameVec,param->getString(iarg++));
	  int ierr = 0;
	  for (int ii = 0; ii < zoneNameVec.size(); ++ii) {
	    // allow for a wildcard...
	    bool found = false;
	    for (int izn = 0; izn < ss.zoneVec.size(); ++izn) {
	      if (MiscUtils::strcmp_wildcard(ss.zoneVec[izn].getName(),zoneNameVec[ii])) {
		ss.zoneVec[izn].flag = 1;
		found = true;
	      }
	    }
	    if (!found) {
	      // maybe it was a zone index...
	      int izn;
	      if (from_string<int>(izn,zoneNameVec[ii],std::dec)&&(izn >= 0)&&(izn < ss.zoneVec.size())) {
		ss.zoneVec[izn].flag = 1;
	      }
	      else {
		WUI(WARN,"no matching ZONE found for " << zoneNameVec[ii]);
		ierr = -1;
	      }
	    }
	  }
	  if (ierr != 0) {
	    return;
	  }
	}
	b_zones = true;
      }
      else if (token == "MOTION") {
	b_motion = true;
      }
      else if (token == "ROTATION") {
	b_rotation = true;
      }
      else if (token == "XC") {
	FOR_I3 xc[i] = param->getDouble(iarg++);
	b_xc = true;
      }
      else if (token == "AXIS") {
	FOR_I3 axis[i] = param->getDouble(iarg++);
	b_axis = true;
      }
      else if (token == "DN") {
	dn = param->getDouble(iarg++);
	b_dn = true;
      }
      else if (token == "DS") {
	ds = param->getDouble(iarg++);
	b_ds = true;
      }
      else {
	WUI(WARN,"unrecognized INTERFACE token \"" << token << "\"");
	return;
      }
    }
    
    assert(b_zones);
    //assert(b_motion);
    //assert(b_rotation);
    assert(b_xc);
    assert(b_axis);
    assert(b_dn);

    // if the segment distance has not been provided, some fraction of the standoff
    // distance is a good guess...
    if (!b_ds) {
      ds = 0.25*dn;
    }
    
    ss.makeFF(xc,axis,dn,ds);

    // what is required here...
    ss.setSubzonesToZones();
    WebUI::webUIOutput.ensureImage();

  }


};


int main(int argc, char * argv[]) {

  try {

    CTI_Init(argc,argv,"surfer.in");

    {
      Surfer surfer;
      surfer.init();
      surfer.run();
    }

    CTI_Finalize();

  }
  catch (int e) {
    if (e >= 0) CTI_Finalize();
    else CTI_Abort();
  }
  catch(...) {
    CTI_Abort();
  }

  return(0);

}
