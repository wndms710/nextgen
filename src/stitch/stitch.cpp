#include "CTI.hpp"
using namespace CTI;

#include "KillfileReader.hpp"
#include "WebUI.hpp"
#include "StitchNS.hpp"

class Histogram2 {
private:

  int nbin;
  bool b_x0x1;
  double x0,x1;
  double * count;

public:

  Histogram2(const double x0,const double x1) {
    this->x0 = x0;
    this->x1 = x1;
    b_x0x1 = true;
    // clear everything else...
    nbin = 100;
    count = NULL;
  }

  ~Histogram2() {
    if (count != NULL) delete[] count;
  }

  void add(const double val) {
    assert(b_x0x1);
    if (val < x0) {
      cout << "val out of range: " << val << endl;
    }
    else if (val >= x1) {
      cout << "val out of range: " << val << endl;
    }
    else {
      const int ibin = int((val-x0)/(x1-x0)*double(nbin));
      if ((ibin >= 0)&&(ibin < nbin))
	add(ibin,1.0);
    }
  }

  void add(const int ibin,const double wgt) {
    assert((ibin >= 0)&&(ibin < nbin));
    if (count == NULL) {
      count = new double[nbin];
      for (int ibin_ = 0; ibin_ < nbin; ++ibin_)
        count[ibin_] = 0.0;
    }
    count[ibin] += wgt;
  }

  void write(const string& filename,const int rank = 0) {
    // reduce and write histogram to a given rank...
    assert(count);
    double * count_reduced = NULL;
    if (mpi_rank == rank) count_reduced = new double[nbin];
    MPI_Reduce(count,count_reduced,nbin,MPI_DOUBLE,MPI_SUM,rank,mpi_comm);
    if (mpi_rank == rank) {
      double sum_wgt = 0.0;
      for (int ibin = 0; ibin < nbin; ++ibin)
        sum_wgt += count_reduced[ibin];
      const double bin_width = (x1-x0)/double(nbin);
      FILE * fp = fopen(filename.c_str(),"w");
      fprintf(fp,"# histogram x0=%lf x1=%lf nbin=%d bin-width=%lf sum(count)=%lf\n",x0,x1,nbin,bin_width,sum_wgt);
      fprintf(fp,"# 1:bin 2:xbin 3:count 4:normalized-count 5:cfd 6:normalized-cdf\n");
      double sum = 0.0;
      for (int ibin = 0; ibin < nbin; ++ibin) {
        sum += 0.5*count_reduced[ibin];
        fprintf(fp,"%d %18.15e %18.15e %18.15e %18.15e %18.15e\n",
                ibin,
                x0+(x1-x0)*(double(ibin)+0.5)/double(nbin),
                count_reduced[ibin],
                count_reduced[ibin]/(sum_wgt*bin_width),
                sum,
                sum/sum_wgt);
        sum += 0.5*count_reduced[ibin];
      }
      fclose(fp);
      delete[] count_reduced;
    }
  }

  void dump() {
    // reduce and dump histogram to rank 0 stdout...
    const int nmax = 120;
    assert(count);
    double * count_reduced = NULL;
    if (mpi_rank == 0) count_reduced = new double[nbin];
    MPI_Reduce(count,count_reduced,nbin,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    if (mpi_rank == 0) {
      double max_count = 0.0;
      for (int ibin = 0; ibin < nbin; ++ibin)
        max_count = max(max_count,count_reduced[ibin]);
      for (int ibin = 0; ibin < nbin; ++ibin) {
        const double xbin = x0+(x1-x0)*(double(ibin)+0.5)/double(nbin);
        cout << setw(7) << xbin << "|";
        const int n = int(count_reduced[ibin]/max_count*double(nmax));
        for (int i = 0; i < n; ++i)
          cout << "#";
        cout << endl;
      }
      delete[] count_reduced;
    }
  }

};

class Stitch {

private:

  list<string> journal;
  KillfileReader kfr;

  bool interactive;
  //bool b_skip_image;
  //string skip_image_name;

public:

  Stitch() {
    interactive = false; // true when stitch is in interactive mode
  }

  void init() {
    COUT1("stitch::init()");
    if (mpi_rank==0) logger->setKillFilename("killstitch");
    MPI_Sync("init()");
  }

  void run() {
    COUT1("stitch::run()");
    FOR_ALL_PARAM {
      param->touch();
      processParam(&(*param));
    }
  }

  void processParam(Param * param) {

    double wtime0;
    if (mpi_rank == 0) {
      journal.push_back(param->str());
      cout << "\nprocessParam \"" << param->str() << "\"" << endl;
      wtime0 = MPI_Wtime();
    }

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

    if ((token == "INTERACTIVE")||(token == "I")) {
      if (!interactive) runInteractive(param);
    }
    else if ((token == "STOP")||(token == "STOP!")) {
      // TODO: is there a better way to do this other than throw(0)?...
      if (b_help) {
        WUI(INFO,"stop or stop! immediately stops stitch and exits");
      }
      else {
        WUI(INFO,"got stop or stop!: stopping stitch and exiting");
        throw(0);
      }
    }
    else if (token == "CTI_LICENSE") {
      COUT1("license location set by input file");
    }
    else if (token == "HCP_PACKING") {
      StitchNS::hcp.processHcpPacking(param,b_help);
    }
    else if (token == "HCP_DELTA") {
      StitchNS::hcp.processHcpDelta(param,b_help);
    }
    else if (token == "HCP_X0") {
      StitchNS::hcp.processHcpX0(param,b_help);
    }
    else if (token == "HCP_DX0") {
      StitchNS::hcp.processHcpDX0(param,b_help);
    }
    else if (token == "DXOST_FACTOR") {
      StitchNS::hcp.processDxostFactor(param,b_help);
    }
    else if (token == "NLAYERS_DEFAULT") {
      StitchNS::hcp.processNlayersDefault(param,b_help);
    }
    else if (token == "HCP_BIT_MAX") {
      StitchNS::hcp.processHcpBitMax(param,b_help);
    }
    else if (token == "HCP_WINDOW") {
      StitchNS::hcp.processHcpWindow(param,b_help);
    }
    else if (token == "REPORT_LEVELS") {
      StitchNS::hcp.processReportLevels(param,b_help);
    }
    else if (token == "COUNT_POINTS") {
      StitchNS::hcp.processCountPoints(param,b_help);
    }
    else if ((token == "Q")||(token == "QUERY")) {
      processQuery(param,b_help);
    }
    else if (token == "PART") {
      PartData::processPart(param,b_help);
      MPI_Sync("done PART");
      StitchNS::hcp.clearPointBooleans();
    }
    else if (token == "INTERFACE") {
      PartData::processInterface(param,b_help);
    }
    else if (token == "FORCE_PLANAR") {
      StitchNS::processForcePlanar(param,b_help);
    }
    else if (token == "PISTON_FF") {
      PartData::processPistonFF(param); // TODO
      StitchNS::hcp.clearPointBooleans();
    }
    else if (token == "WRITE_SBIN") {
      PartData::processWriteSbin(param,b_help);
    }
    else if (token == "WRITE_PBIN") {
      PartData::processWritePbin(param,b_help);
    }
    else if (token == "WRITE_HCP_PTS") {
      StitchNS::hcp.processWritePts(param);
    }
    else if (token == "READ_HCP_PTS") {
      StitchNS::hcp.processReadPts(param);
    }
    else if (token == "QUERY_LENGTH_SCALES") {
      PartData::queryLengthScales(); // TODO
    }
    else if ( (token == "WRITE_MLES")      ||
              (token == "WRITE_RESTART")   ||
              (token == "WRITE_RESTART_V5")||
              (token == "WRITE_RESTART_V4")||
              (token == "WRITE_RESTART_V3") ) {
      StitchNS::processWriteMles(param,b_help);
    }
    else if (token == "RUN_DIAGNOSTICS") {
      PartData::runDiagnostics(param,b_help);
    }
    else if (token == "CLEAR_SMOOTHING") {
      StitchNS::processClearSmoothing(param,b_help);
    }
    else if (token == "SMOOTH_LIMIT") {
      StitchNS::processSmoothLimit(param,b_help);
    }
    else if (token == "SMOOTH_MODE") {
      StitchNS::processSmoothMode(param,b_help);
    }
    else if (token == "CREASE_ANGLE_DEGREES") {
      StitchNS::processCreaseAngleDegrees(param,b_help);
    }
    else if (token == "RC_LOCAL_THRESHOLD") {
      StitchNS::processRcLocalThreshold(param,b_help);
    }
    else if (token == "RC_LIMIT") {
      StitchNS::processRcLimit(param,b_help);
    }
    else if (token == "LOAD_BALANCE_MODE") {
      PartData::processLoadBalanceMode(param,b_help);
    }
    else if (token == "NSMOOTH") {
      processNsmooth(param,b_help);
    }
    else if ((token == "STATUS")||(token == "VERSION")) {
      processStatus(param,b_help);
    }
    else if (token == "WRITE_JSON") {
      StitchNS::processWriteJson(param,b_help);
    }
    else if (token == "WRITE_IMAGE") {
      StitchNS::processWriteImage(param,b_help);
    }
    else if (token == "WRITE_JOURNAL") {
      processWriteJournal(param,b_help);
    }
    else if (token == "WRITE_PARAMS") {
      processWriteParams(param,b_help);
    }
    else if (token == "WRITE_PART") {
      processWritePart(param,b_help);
    }
    else if (token == "WRITE_SINGLE_PART") {
      processWriteSinglePart(param,b_help);
    }
    else if (token == "FLAG_ZONES_WITH_PART") {
      processFlagZonesWithPart(param,b_help);
    }
    else if (token == "WRITE_PTS_TECPLOT") {
      processWritePtsTecplot(param,b_help);
    }
    else if (token == "SHOW") {
      StitchNS::processShow(param,b_help);
    }
    else {
      if (b_help) {
        WUI(WARN,"HELP not available for unrecognized param: " << token);
      }
      else {
        WUI(WARN,"skipping unrecognized param: " << token);
      }
    }

    // report timing...
    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0) cout << "processParam elapsed time " << MPI_Wtime()-wtime0 << " [s]\n" << endl;

  }

  void runInteractive(Param * _param) {

    COUT1("starting interactive session...");

    if (mpi_rank==0) logger->insertSolverAction(SIMULATION_RUN);

    assert(!interactive);
    interactive = true;

    // clear the webUIOutput before starting the interactive session:
    // this will prevent the deluge of INFO/WARNING messages that can occur
    // when the session starts...
    WebUI::webUIOutput.clear();

    if (_param->size() > 0) {
      double timeout_mins = _param->getDouble(0);
      kfr.setHoldWithTimeout(true,timeout_mins);
    }
    else{
      kfr.setHold(true);
    }

    while (Param * param = kfr.getParam("killstitch")) {
      processParam(param);
    }

    interactive = false;

  }

  void helpHelp() {

    WUI(WARN,"HELP must be followed by a parameter. HELP is available for\n" <<
        "  PART          HCP_DELTA     STATUS          VERSION\n" <<
        "  NSMOOTH       HCP_WINDOW    WRITE_IMAGE     HOLD" <<
        "  WRITE_JOURNAL WRITE_JSON");

  }

  void processStatus(Param *param,const bool b_help) {
    if (b_help) {
      WUI(INFO,"STATUS prints some information about your current stitch run.");
    }
    else {
      WUI(INFO,"stitch running on " << mpi_size << " cores, " << mpi_size_internode << " nodes.\n" <<
          "version: " << cti_core_version << " built: " << cti_core_date);
    }
  }

  void processNsmooth(Param * param,const bool b_help) {
    if (b_help) {
      helpNsmooth();
    }
    else {
      try {

        int nsmooth = param->getInt(0);
        WUI(INFO,"NSMOOTH set to " << nsmooth);

        // make sure hcpPointsBuilder points are up-to-date...
        StitchNS::hcp.ensurePoints();

        // if hcpPointBuilder changed anything, then we start over...
        if (PartData::hcpPts == NULL) {
          StitchNS::hcp.writeHcpPointsToPartData();
          StitchNS::clearSmoothing(); // must be called before clear points
          PartData::clearPoints();
          StitchNS::finalize();
        }
        // otherwise, we check to see if the smoothing has been cleared...
        else if (!PartData::b_vp) {
          StitchNS::clearSmoothing();
          StitchNS::finalize();
        }

        // make sure we have points copies from PartData...
        PartData::ensurePoints();

        // smooth points
        StitchNS::rebuildVd(nsmooth);

        // set bool that vd's are available 
        PartData::b_vp = true;

        // update image after smoothing...
        WebUI::webUIOutput.ensureImage();

      }
      catch(int e) {
        if (e == 0) {
          WUI(WARN,"expecting NSMOOTH <int>");
        }
        else if (e == -1) {
          // failure in nsmooth
          throw(-1);
        }
        else {
          WUI(WARN,"invalid param for NSMOOTH: " << param->getString() <<
              "\nexpecting NSMOOTH <int>");
        }
        helpNsmooth();
      }
    }
  }

  void helpNsmooth() const {
    WUI(INFO,
        "NSMOOTH performs Lloyd iteration on current HCP points. Examples:\n" <<
        "NSMOOTH 0\n" <<
        "NSMOOTH 20\n" <<
        "for more detail see [$CWIKB:stitch_smoothing]");
  }
  
  void processWritePtsTecplot(Param * param,const bool b_help) {
    if (b_help) {
      WUI(INFO,"WRITE_PTS_TECPLOT <tecplot-filename>");
    }
    else {
      try {
	const string filename = param->getString(0);

        // make sure hcpPointsBuilder points are up-to-date...
        StitchNS::hcp.ensurePoints();

        // if hcpPointBuilder changed anything, then we start over...
        if (PartData::hcpPts == NULL) {
          StitchNS::hcp.writeHcpPointsToPartData();
          StitchNS::clearSmoothing(); // must be called before clear points
          PartData::clearPoints();
          StitchNS::finalize();
        }
        // otherwise, we check to see if the smoothing has been cleared...
        else if (!PartData::b_vp) {
          StitchNS::clearSmoothing();
          StitchNS::finalize();
        }

        // make sure we have points copies from PartData...
        PartData::ensurePoints();
	
	GeomUtils::writePtsTecplot(filename,PartData::x_vd,PartData::ncv);

      }
      catch(int e) {
	WUI(WARN,"failed WRITE_PTS_TECPLOT <filename>");
      }
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
        bool b_quiet = false;
        int iarg = 0;
        while (iarg < param->size()) {
          string token = param->getString(iarg++);
          if ((token == "VERBOSE")||(token == "verbose")) {
            b_verbose = true;
          }
          else if ((token == "QUIET")||(token == "quiet")) {
            b_quiet = true; //suppress WUI output when called by GUI
          }
          else {
            // assume this is the name...
            name = token;
          }
        }
        FILE * fp = fopen(name.c_str(),"w");
        int count = 0;
        for (list<string>::iterator iter = journal.begin(); iter != journal.end(); ++iter) {
          // default is to skip all WRITE_IMAGE and WRITE_JSON lines...
          if ( (!b_verbose) && ((iter->find("WRITE_IMAGE") == 0)||(iter->find("write_image") == 0)||
                                (iter->find("WRITE_JSON") == 0)||(iter->find("write_json") == 0)||
                                (iter->find("INTERACTIVE") == 0)||(iter->find("interactive") == 0)||
                                (iter->find("WRITE_JOURNAL") == 0)||(iter->find("write_journal") == 0)||
                                (iter->find("WRITE_PARAMS") == 0)||(iter->find("write_params") == 0)) )
            continue;
          fprintf(fp,"%s\n",iter->c_str());
          ++count;
        }
        fclose(fp);
        if (!b_quiet)
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

  void processWriteParams(Param * param,const bool b_help) {

    if (b_help) {
      WUI(INFO,
          "WRITE_PARAMS writes out the current param list. Examples:\n" <<
          "WRITE_PARAMS");
    }
    else {

      // only rank 0 write the json...
      string filename = "params.in"; // TODO: if this name exists, don't overwrite. Turn this to journal(1).in or something

      try {

        int iarg = 0;
        while (iarg < param->size()) {
          string token = param->getString(iarg++);
          if (token == "NAME") {
            const string prefix = param->getString(iarg++);
            filename = prefix+".in";
          }
        }

      }
      catch(int e) {
        WUI(WARN,"expecting WRITE_PARAMS NAME <prefix>");
        return;
      }

      if (mpi_rank == 0) {

        string tmp_filename = "." + filename;
        FILE * fp = fopen(tmp_filename.c_str(),"w");
        assert(fp);

        stringstream ss;
        FOR_ALL_PARAM {
          param->dump(ss);
        }

        fprintf(fp,"%s\n",ss.str().c_str());
        fclose(fp);
        rename(tmp_filename.c_str(),filename.c_str());

      }

    }

  }

  void processQuery(Param * param,const bool b_help) {

    if (b_help || (param->size() == 0)) {
      WUI(INFO,
          "QUERY (or Q for short) returns information about the current model\n" <<
          "  examples:\n" <<
          "    Q PARTS           # prints the list of parts\n" <<
          "    Q ZONES           # prints the list of zones\n" <<
          "    Q BBOX            # prints the bboxes for all zones\n" <<
          "    Q ZONEID 0,1      # prints analysis of a list of comma-separated zone ids\n" <<
          "    Q BBOX ZONEID 0,1 # prints the bboxes for a list of comma-separated zone ids");
      return;
    }

    bool b_bbox = false;
    vector<int> zoneidVec;
    int ierr = 0;

    int iarg = 0;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "PARTS") {
        PartData::queryParts();
      }
      else if ((token == "ZONE")||(token == "ZONES")) {
        // if there is another argument, then parse it and send it...
        if (iarg < param->size()) {
          const string names = param->getString(iarg++);
          PartData::queryZones(names);
        }
        else {
          PartData::queryZones();
        }
      }
      else if ((token == "HCP_WINDOW")||(token == "HCP_WINDOWS")) {
        StitchNS::hcp.queryHcpWindows();
      }
      else if (token == "BBOX") {
        b_bbox = true;
      }
      else if (token == "ZONEID") {
        // expect a comma-delimited set of integers...
        const string zonesCsv = param->getString(iarg++);
        vector<int> tmp;
        MiscUtils::splitCsv(tmp,zonesCsv);
        for (int ii = 0; ii < tmp.size(); ++ii) {
          const int izone = tmp[ii];
          if ((izone >= 0)&&(izone < PartData::zoneVec.size())) {
            zoneidVec.push_back(izone);
          }
          else {
            WUI(WARN,"ZONE index out of range: " << izone << ". 0 <= izone < " << PartData::zoneVec.size());
            ierr = 1;
          }
        }
        if(tmp.size() <= 0)
        {
          WUI(WARN,"Zone ID CSV string [ZONEID " << zonesCsv << "] must have valid zone id values separated by commas, e.g. [ZONEID 0], [ZONEID 0,1], etc.");
          ierr = 1;
        }
      }
      else {
        WUI(WARN,"Skipping unsupported QUERY: " << token);
      }
    }

    if (ierr == 0) {
      if (b_bbox) {
        double bbmin[3],bbmax[3];
        if (zoneidVec.empty()) {
          // dump the entire part bounding box if no zones provided...
          PartData::getBbox(bbmin,bbmax);

          WUI(INFO,
              "BBOX X0 " << bbmin[0] << " " << bbmin[1] << " " << bbmin[2] << "\n"
              "     X1 " << bbmax[0] << " " << bbmax[1] << " " << bbmax[2]);
        }
        else {
          // user has provided specific zones...
          PartData::getBboxForZoneids(bbmin,bbmax,zoneidVec);

          std::ostringstream zoneIdStr("[",std::ostringstream::ate);
          for(int izn=0; izn < zoneidVec.size()-1; ++izn) zoneIdStr << zoneidVec[izn] << ",";
          zoneIdStr << zoneidVec.back() << "]";
          WUI(INFO,
              "BBOX for zoneids: " << zoneIdStr.str() << "\n"
              "                  X0 " << bbmin[0] << " " << bbmin[1] << " " << bbmin[2] << "\n"
              "                  X1 " << bbmax[0] << " " << bbmax[1] << " " << bbmax[2]);
        }
      }
      else {
        //
      }
    }
  }

  void processWritePart(Param * param,const bool b_help) {

    // look for help...
    if (b_help) {
      WUI(INFO,
          "WRITE_PART writes the current surface(s) and points to a flattened part file\n" <<
          "  example:\n" <<
          "    WRITE_PART mypart.part\n"
          "for more detail see [$CWIKB:stitch_export]");
      return;
    }

    // the param requires the part filename...
    if (param->size() == 0) {
      WUI(WARN,"Error: WRITE_PART requires a filename. e.g. WRITE_PART fanblades.part");
      return;
    }

    const string filename = param->getString(0);

    // make sure hcpPointsBuilder points are up-to-date...
    StitchNS::hcp.ensurePoints();

    // if hcpPointBuilder changed anything, then we start over...
    if (PartData::hcpPts == NULL) {
      StitchNS::hcp.writeHcpPointsToPartData();
      StitchNS::clearSmoothing(); // must be called before clear points
      PartData::clearPoints();
      StitchNS::finalize();
    }
    // otherwise, we check to see if the smoothing has been cleared...
    else if (!PartData::b_vp) {
      StitchNS::clearSmoothing();
      StitchNS::finalize();
    }

    // make sure we have points copies from PartData...
    PartData::ensurePoints();

    {

      char dummy[128];
      sprintf(dummy,"%s",filename.c_str());
      MPI_File fh;
      MPI_Offset offset = 0;
      MPI_File_open(mpi_comm,dummy,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      // ===============================================================
      // header
      // ===============================================================

      int itmp[5] = { PART_IO_MAGIC_NUMBER, PART_IO_VERSION, 0, 0, 0 };
      for (int ipart = 0; ipart < PartData::partVec.size(); ++ipart) {
        if (PartData::partVec[ipart]->surface && PartData::partVec[ipart]->b_write_surface) {
          // look for the keyword FF...
          // TODO: for now, just assume all surface...
          itmp[2] = 1;
        }
        if (PartData::partVec[ipart]->ff_surface && PartData::partVec[ipart]->b_write_ff) {
          // look for the keyword FF...
          // TODO: for now, just assume all surface...
          itmp[3] = 1;
        }
        if (PartData::partVec[ipart]->pts) {
          // look for the keyword FF...
          // TODO: for now, just assume all surface...
          itmp[4] = 1;
        }
        if (PartData::hcpPts) itmp[4] = 1;
      }
      assert(itmp[2] == 1); // surface
      // assert(itmp[3] == 0); // ff_surface
      if (itmp[3] == 1) {
        WUI(WARN,"currently WRITE_PART does not support Parts with farfield surface (may use NO_WRITE_FF to ignore ff); skipping");
        MPI_File_close(&fh);
        throw(0);
      }
      assert(itmp[4] == 1); // pts (any pts, not just hcpPts)

      if (mpi_rank == 0) {

        cout << " > surface: " << itmp[2] << ", ff_surface: " << itmp[3] << ", pts: " << itmp[4] << endl;

        MPI_File_write(fh,itmp,5,MPI_INT,MPI_STATUS_IGNORE);
        offset += int_size*5;

        if (itmp[2] == 1) {

          // ===============================================================
          // surface
          // ===============================================================

          int surface_itmp[3] = { 0, 0, 0 };
	  int count = 0;
          for (int ipart = 0; ipart < PartData::partVec.size(); ++ipart) {
            if (PartData::partVec[ipart]->surface) {
              // HACK: to get going, assume only ipart 0 is written...
              //if (ipart != 0) {
              //  cout << "Warning: assuming part: " << PartData::partVec[ipart]->name << " is FF" << endl;
              //  continue;
              //}
	      if (PartData::partVec[ipart]->b_write_surface) {
		++count;
		if (count > 1) {
		  cout << "> surface_itmp..." << endl;
		  cout << "  Warning: currently only support surfaces from one part; skipping part: " << PartData::partVec[ipart]->name << endl;
		  continue;
		}
		SurfaceShm * surface = PartData::partVec[ipart]->surface;
		// TODO: here we should separate out anything FF?...
		// for now, just take it all...
		surface_itmp[0] += surface->zoneVec.size();
		surface_itmp[1] += surface->nsp;
		surface_itmp[2] += surface->nst;
	      }
	      else {
		cout << "Surface of part: " << PartData::partVec[ipart]->name << " is not written" << endl;
	      }
            }
          }

          cout << " > surface nzones: " << surface_itmp[0] << ", nsp: " << surface_itmp[1] << ", nst: " << surface_itmp[2] << endl;

          MPI_File_write(fh,surface_itmp,3,MPI_INT,MPI_STATUS_IGNORE);
          offset += int_size*3;

          // compressed surface zones...

          int nzone = 0;
	  count = 0;
          for (int ipart = 0; ipart < PartData::partVec.size(); ++ipart) {
            if (PartData::partVec[ipart]->surface) {
              // HACK: to get going, assume only ipart 0 is written...
              //if (ipart != 0) {
              //  if (mpi_rank == 0) cout << "Warning: assuming part: " << PartData::partVec[ipart]->name << " is FF" << endl;
              //  continue;
              //}
	      if (PartData::partVec[ipart]->b_write_surface) {
		++count;
		if (count > 1) {
		  cout << "> zone names..." << endl;
		  cout << "  Warning: currently only support surfaces from one part; skipping part: " << PartData::partVec[ipart]->name << endl;
		  continue;
		}
		SurfaceShm * surface = PartData::partVec[ipart]->surface;
		for (int izone = 0; izone < surface->zoneVec.size(); ++izone) {
		  int length = surface->zoneVec[izone].getName().length();
		  MPI_File_write(fh,&length,1,MPI_INT,MPI_STATUS_IGNORE);
		  offset += int_size;
		  MPI_File_write(fh,(char*)surface->zoneVec[izone].getName().c_str(),length,MPI_CHAR,MPI_STATUS_IGNORE);
		  offset += length;
		  cout << " > zone " << nzone << " name \"" << surface->zoneVec[izone].getName() << "\"" << endl;
		  ++nzone;
		}
	      }
            }
          }
          assert(nzone == surface_itmp[0]);

          // compressed surface xsp...
          // we could do this in parallel, but serial for now...

          int nsp = 0;
	  count = 0;
          for (int ipart = 0; ipart < PartData::partVec.size(); ++ipart) {
            if (PartData::partVec[ipart]->surface) {
              // HACK: to get going, assume only ipart 0 is written...
              //if (ipart != 0) {
              //  if (mpi_rank == 0) cout << "Warning: assuming part: " << PartData::partVec[ipart]->name << " is FF" << endl;
              //  continue;
              //}
	      if (PartData::partVec[ipart]->b_write_surface) {
		++count;
		if (count > 1) {
		  cout << "> xsp..." << endl;
		  cout << "  Warning: currently only support surfaces from one part; skipping part: " << PartData::partVec[ipart]->name << endl;
		  continue;
		}
		SurfaceShm * surface = PartData::partVec[ipart]->surface;
		MPI_File_write(fh,(double*)surface->xsp,surface->nsp*3,MPI_DOUBLE,MPI_STATUS_IGNORE);
		nsp += surface->nsp;
	      }
            }
          }
          assert(nsp == surface_itmp[1]);
          offset += double_size*nsp*3;

          // compressed surface spost...
          int nst = 0;
	  count = 0;
          for (int ipart = 0; ipart < PartData::partVec.size(); ++ipart) {
            if (PartData::partVec[ipart]->surface) {
              // HACK: to get going, assume only ipart 0 is written...
              //if (ipart != 0) {
              //  if (mpi_rank == 0) cout << "Warning: assuming part: " << PartData::partVec[ipart]->name << " is FF" << endl;
              //  continue;
              //}
	      if (PartData::partVec[ipart]->b_write_surface) {
		++count;
		if (count > 1) {
		  cout << "> spost..." << endl;
		  cout << "  Warning: currently only support surfaces from one part; skipping part: " << PartData::partVec[ipart]->name << endl;
		  continue;
		}
		SurfaceShm * surface = PartData::partVec[ipart]->surface;
		MPI_File_write(fh,(int*)surface->spost,surface->nst*3,MPI_INT,MPI_STATUS_IGNORE);
		nst += surface->nst;
	      }
            }
          }
          assert(nst == surface_itmp[2]);
          offset += int_size*nst*3;

          nst = 0;
	  count = 0;
          for (int ipart = 0; ipart < PartData::partVec.size(); ++ipart) {
            if (PartData::partVec[ipart]->surface) {
              // HACK: to get going, assume only ipart 0 is written...
              //if (ipart != 0) {
              //  if (mpi_rank == 0) cout << "Warning: assuming part: " << PartData::partVec[ipart]->name << " is FF" << endl;
              //  continue;
              //}
	      if (PartData::partVec[ipart]->b_write_surface) {
		++count;
		if (count > 1) {
		  cout << "> znost..." << endl;
		  cout << "  Warning: currently only support surfaces from one part; skipping part: " << PartData::partVec[ipart]->name << endl;
		  continue;
		}
		SurfaceShm * surface = PartData::partVec[ipart]->surface;
		MPI_File_write(fh,(int*)surface->znost,surface->nst,MPI_INT,MPI_STATUS_IGNORE);
		nst += surface->nst;
	      }
            }
          }
          assert(nst == surface_itmp[2]);
          offset += int_size*nst;

        }

        if (itmp[3] == 1) {

          // ===============================================================
          // ff_surface
          // ===============================================================

          assert(0);

        }

      }

      MPI_Bcast(&offset,1,MPI_INT8,0,mpi_comm);

      // ===============================
      // finally points...
      // these come from the PartData...
      // ===============================

      if (itmp[4] == 1) {

        int my_np = PartData::ncv;
        int my_disp;
        MPI_Scan(&my_np,&my_disp,1,MPI_INT,MPI_SUM,mpi_comm);

        int np_global = my_disp;
        MPI_Bcast(&np_global,1,MPI_INT,mpi_size-1,mpi_comm);
        assert(np_global > 0);

        if (mpi_rank == 0) {
          cout << " > total point count: " << np_global << endl;
          MPI_File_write(fh,&np_global,1,MPI_INT,MPI_STATUS_IGNORE);
        }
        offset += int_size;

        // now xp...

        writeChunkedData<double>(fh,offset+double_size*(my_disp-my_np)*3,(double*)PartData::x_vd,my_np*3,mpi_comm);
        //MPI_File_write_at_all(fh,offset,xp,np*3,MPI_DOUBLE,MPI_STATUS_IGNORE);
        offset += double_size*np_global*3;

        // delta...

        writeChunkedData<double>(fh,offset+double_size*(my_disp-my_np),PartData::delta_vd,my_np,mpi_comm);
        //MPI_File_write_at_all(fh,offset,delta,np,MPI_DOUBLE,MPI_STATUS_IGNORE);
        offset += double_size*np_global;

      }

      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);

    }

    WUI(INFO," > file " << filename << " written successfully");

  }

  void processWriteSinglePart(Param * param,const bool b_help) {

    // look for help...
    if (b_help) {
      WUI(INFO,
          "WRITE_SINGLE_PART writes out the named part to a part file\n" <<
          "  example:\n" <<
          "    WRITE_SINGLE_PART bl\n"
          "for more detail see [$CWIKB:stitch_export]");
      return;
    }

    // the param requires the part filename...
    if (param->size() == 0) {
      WUI(WARN,"Error: WRITE_SINGLE_PART requires a filename. e.g. WRITE_SINGLE_PART bl");
      return;
    }

    const string partname = param->getString(0);
    int ipart = -1;
    for (ipart = 0; ipart < PartData::partVec.size(); ++ipart) {
      if (PartData::partVec[ipart]->getName() == partname)
        break;
    }

    if (ipart == PartData::partVec.size()) {
      WUI(WARN,"Error: WRITE_SINGLE_PART part was not found");
      return;
    }

    const string filename = partname+".part";
    PartData::partVec[ipart]->writeBinary(filename);
    WUI(INFO," > file " << filename << " written successfully");

  }

  void processFlagZonesWithPart(Param * param,const bool b_help) {

    if (b_help || (param->size() == 0)) {
      WUI(INFO,"FLAG_ZONES_WITH_PART ZONES <zones list> PART <name> sets the surface flag with part index.");
      return;
    }

    bool b_zone_indices = false;
    set<int> zone_indices;
    int ipart = -1;

    int iarg = 0;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if ((token == "FAZONE")||(token == "ZONE")||(token == "ZONES")) {
        b_zone_indices = true;
        const string zoneNamesCsv = param->getString(iarg++);
        vector<string> tmp;
        MiscUtils::splitCsv(tmp,zoneNamesCsv);
        for (int ii = 0; ii < tmp.size(); ++ii) {
          const int izone = PartData::getZoneIndex(tmp[ii]); // get the unique zoneIndex in the partVec
          if (izone >= 0) {
            zone_indices.insert(izone);
          }
          else {
            WUI(WARN,"FLAG_ZONES ZONE not recognized: " << tmp[ii]);
            return;
          }
        }

      }
      else if (token == "PART") {
        const string part_name = param->getString(iarg++);
        for (ipart = 0; ipart < PartData::partVec.size(); ++ipart) {
          if (PartData::partVec[ipart]->getName() == part_name)
            break;
        }
        if (ipart == PartData::partVec.size()) {
          WUI(WARN,"Error: FLAG_ZONES_WITH_PART PART not recognized: " << part_name);
          return;
        }
      }
    }

    if (ipart == -1) {
      WUI(WARN,"Error: FLAG_ZONES_WITH_PART requires PART <part-name>: ");
      return;
    }
    if (!b_zone_indices) {
      WUI(WARN,"Error: FLAG_ZONES_WITH_PART requires ZONES <zone-names>: ");
      return;
    }

    if (mpi_rank_shared == 0) {
      for (set<int>::iterator it = zone_indices.begin(); it != zone_indices.end(); ++it) {
        pair<int,int> partZonePair = PartData::getPartSurfaceZonePair(*it);
        if (PartData::partVec[partZonePair.first]->surface) {
          const int izone = partZonePair.second;
          for (int ist = 0; ist < PartData::partVec[partZonePair.first]->surface->nst; ++ist) {
            if (PartData::partVec[partZonePair.first]->surface->znost[ist] == izone) {
              PartData::partVec[partZonePair.first]->setPartForSt(ipart,ist);
            }
          }
        }
      }
    }
    MPI_Barrier(mpi_comm_shared);

    assert(!PartData::partVec[ipart]->b_link_surface);
    assert(PartData::partVec[ipart]->ipart_link.empty());
    assert(PartData::partVec[ipart]->ist_offset_link.empty());

    for (set<int>::iterator it = zone_indices.begin(); it != zone_indices.end(); ++it) {
      pair<int,int> partZonePair = PartData::getPartSurfaceZonePair(*it);
      if (PartData::partVec[partZonePair.first]->surface) {
	// record the associated parts
	PartData::partVec[ipart]->b_link_surface = true;
	PartData::partVec[ipart]->ipart_link.push_back(partZonePair.first);
	if (PartData::partVec[ipart]->ist_offset_link.empty()) {
	  PartData::partVec[ipart]->ist_offset_link.push_back(0);
	}
	else {
	  const int back = PartData::partVec[ipart]->ist_offset_link.back();
	  PartData::partVec[ipart]->ist_offset_link.push_back(back+PartData::partVec[partZonePair.first]->surface->nst);
	}
      }
    }

    WUI(INFO," > flagged zones with part");

  }

};

int main(int argc, char * argv[]) {

  try {

    CTI_Init(argc,argv,"stitch.in");

    {

      Stitch stitch;
      stitch.init();
      stitch.run();
      
      StitchNS::finalize();
      PartData::clear();
      PartData::clear_new(); // for new data: should be moved to Cascade NS some day
      
    }

    CTI_Finalize();

  }
  catch (int e) {
    // integer exceptions are thrown from some parts of the code.
    // negative exceptions are thrown when only one rank is having a problem,
    // so these require a CTI_Abort to free the resources. Positive exceptions
    // e.g. throw(0), throw(1) are associated with points where everyone is
    // synchronized, and we can call CTI_Finalize and shut down cleanly.
    if (e >= 0) { 
      CTI_Finalize();
    }
    else {
      cout << "catch " << e << " on rank: " << mpi_rank << ": calling abort..." << endl;
      CTI_Abort();
      CTI_Finalize();
    }
  }
  catch(...) {
    // don't know what this could be...
    cout << "catch unknown error on rank: " << mpi_rank << ": calling abort..." << endl;
    CTI_Abort();
    CTI_Finalize();
  }

  return(0);

}
