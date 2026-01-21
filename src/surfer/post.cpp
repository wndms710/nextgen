#include "CTI.hpp"
using namespace CTI;

#include "SimpleRestart.hpp"

#include "IntFlag.hpp"
#include "CtiScene.hpp"
#include "KillfileReader.hpp"
#include "LesImageMapper.hpp"
#include "SimplePoint.hpp"
#include "WebUI.hpp"

class Post {

private:

  KillfileReader kfr;
  list<string> journal;

  bool b_skip_image;
  string skip_image_name;

public:

  string filetype;

  SimpleRestart rs;

  CvPlaneMap cvPlaneMap; //persist cvs associated with a plane to speed up serial planar vis

  Post() {
    clear();
  }

  ~Post() {
  }

  void init() {
    COUT1("post::init()");
    if (mpi_rank==0) logger->setKillFilename("killpost");
  }

  void clear() {
    COUT1("post::clear()");
    b_skip_image = false;
    skip_image_name = "";
    rs.clear();
  }

  void run() {
    COUT1("post::run()");
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

    while (Param * param = kfr.getParam("killpost")) {
      string token = MiscUtils::toUpperCase(param->getName());
      if ((token != "INTERACTIVE")&&(token != "I")) processParam(param);
    }

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

    if (token == "RESTART") {
      processRestart(param,b_help);
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
        WUI(INFO,"stop or stop! immediately exits Post");
      }
      else {
        WUI(INFO,"got stop or stop!: stopping post and exiting");
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
    else if (token == "WRITE_PARAMS") {
      processWriteParams(param,b_help);
    }
    else if ((token == "JOURNAL")||(token == "HISTORY")) {
      processHistory(param,b_help);
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

    //// do a little checking...
    //if (ss.nsp || ss.nst) ss.check();

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

  void processRestart(Param * param,const bool help) {
    if (help) {
      helpRestart();
      return;
    }

    try {
      int ierr = rs.init(param);
      if (ierr!=0) {
        helpRestart();
        return;
      }
    }
    catch(int e) {
      if (e == 0) {
        WUI(WARN,"Expecting RESTART <restart.mles> [<result.sles>]");
      }
      else {
        WUI(WARN,"Expecting RESTART <restart.mles> [<result.sles>]");
      }
      helpRestart();
    }
  }

  void helpRestart() {
    WUI(INFO,"RESTART initializes post processing from an mles file and optionally an sles file.\n" <<
        "The expected syntax is RESTART <restart.mles> [<result.sles>]. The domain surface\n" <<
        "is read by post immediately while volume data is simply indexed for fast visualization.\n"
        );
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

  void processExtractPts(Param * param,const bool help) {
    if (!rs.isInit()) {
      WUI(WARN,"no surface available currently; please use SURF to import or construct one");
    }

    if (help) {
      WUI(INFO,"extracts points from a Cascade mesh or restart file and writes them to Cascade's binary pbin format\n"
      );
      return;
    }

    bool b_restart = false;
    string name;

    int iarg = 0;
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "RESTART") {
        name = param->getString(iarg++);
        b_restart = true;
      }
      else {
        cout << "ERROR: EXTRACT_PTS: unrecognized token \"" << token << "\". Check syntax." << endl;
        return;
      }
    }

    vector<SimplePoint> pointVec;
    if (b_restart) {

      // for extracting points from a restart file, use the LesImageMapper...
      LesImageMapper lim(name);
      lim.getPointsInsideGeom(pointVec,&(rs.ss));


    }
    else {
      cout << "ERROR: Skipping EXTRACT_PTS: Check syntax." << endl;
      return;
    }

    // ==========================================
    // finally write the pbin file....
    // ==========================================

    cout << " > pointVec.size(): " << pointVec.size() << endl;

    /*
      FILE * fp = fopen("points.dat","w");
      for (int ip = 0, limit = pointVec.size(); ip < limit; ++ip) {
      fprintf(fp,"%18.15e %18.15e %18.15e\n",pointVec[ip].x[0],pointVec[ip].x[1],pointVec[ip].x[2]);
      }
      fclose(fp);
      cout << "take a look" << endl;
      getchar();
    */

    cout << " > writing file \"extract_pts.pbin\"..." << endl;
    FILE * fp = fopen("extract_pts.pbin","wb");

    int8 ibuf[4] = { POINTS_IO_MAGIC_NUMBER, POINTS_IO_VERSION, (int8)pointVec.size(), 1 }; // last is include delta or not maybe?
    fwrite(ibuf,sizeof(int8),4,fp);

    // xp...
    double *dbuf = new double[pointVec.size()*3];
    for (int ip = 0, limit = pointVec.size(); ip < limit; ++ip)
      FOR_I3 dbuf[ip*3+i] = pointVec[ip].x[i];
    fwrite(dbuf,sizeof(double),pointVec.size()*3,fp);

    // delta...
    for (int ip = 0, limit = pointVec.size(); ip < limit; ++ip)
      dbuf[ip] = 1.0E-6;
    fwrite(dbuf,sizeof(double),pointVec.size(),fp);

    fclose(fp);

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

  bool isActiveCommand(const string& cmd) {
    // grow this list of inactive commmands over time...
    // it is used by WRITE_JOURNAL, UNDO, REWIND, HISTORY/JOURNAL...
    if (MiscUtils::startsWith(cmd,"WRITE_IMAGE")||MiscUtils::startsWith(cmd,"write_image")||
        MiscUtils::startsWith(cmd,"INTERACTIVE")||MiscUtils::startsWith(cmd,"interactive")||
        (cmd == "I")||(cmd == "i")||
        MiscUtils::startsWith(cmd,"HELP")||MiscUtils::startsWith(cmd,"help")||
        MiscUtils::startsWith(cmd,"WRITE_JSON")||MiscUtils::startsWith(cmd,"write_json")||
        MiscUtils::startsWith(cmd,"WRITE_SBIN")||MiscUtils::startsWith(cmd,"write_sbin")||
        MiscUtils::startsWith(cmd,"WRITE_TECPLOT")||MiscUtils::startsWith(cmd,"write_tecplot")||
        MiscUtils::startsWith(cmd,"UNDO")||MiscUtils::startsWith(cmd,"undo")||
        MiscUtils::startsWith(cmd,"HISTORY")||MiscUtils::startsWith(cmd,"history")||
        MiscUtils::startsWith(cmd,"JOURNAL")||MiscUtils::startsWith(cmd,"journal")||
        MiscUtils::startsWith(cmd,"REWIND")||MiscUtils::startsWith(cmd,"rewind")||
        MiscUtils::startsWith(cmd,"RUN_DIAGNOSTICS")||MiscUtils::startsWith(cmd,"run_diagnostics")||
        MiscUtils::startsWith(cmd,"BBOX")||MiscUtils::startsWith(cmd,"bbox")||
        MiscUtils::startsWith(cmd,"WRITE_JOURNAL")||MiscUtils::startsWith(cmd,"write_journal"))
      return false;
    return true;
  }

  void processWriteJson(Param * param,const bool help) {
    if (help) {
      WUI(INFO,"a method to write surface data to a Json file. This is meant to be used by the Cascade App to communicate with post"
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
      string filename = "post";

      //Check for hide zone information for bounding box
      //computation
      int iarg = 0;
      while (iarg < param->size()) {
        string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (token == "NAME") {
          filename = param->getString(iarg++);
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
            for (int izn = 0, nzn = rs.ss.zoneVec.size(); izn < nzn; ++izn) {
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
      FILE * fp = fopen((filename+".json").c_str(),"w");
      assert(fp);

      double bbox[6];
      double xbb[3] = {0.0, 0.0, 0.0};
      double diag   = 0.0;
      bool b_surface_loaded = true;
      if (!rs.isInit()){ //no restart has been loaded
        b_surface_loaded = false;
      }
      else{
        if (!b_hidesubzones) {
          rs.ss.getBoundingBox(bbox);
          rs.ss.getBoundingBoxCenter(xbb);
          diag = 2.0*rs.ss.getBoundingBoxRmax();
        }
        else {
          // flag stores which zones should be hidden (i.e., not added)
          const int nzn = rs.ss.szozn_i.getLength() - 1;
          const int nsz = rs.ss.szozn_i[nzn];
          IntFlag show_sz_flag(nsz);
          show_sz_flag.setAll(1);
          if (b_hidesubzones) {
            std::set<int>::iterator it;
            // flag subzones based on hidden zones
            for (int izn = 0; izn < nzn; ++izn) {
              it = hiddenZonesSet.find(izn);
              if ( it!=hiddenZonesSet.end() ) {
                for (int iszn = rs.ss.szozn_i[izn]; iszn < rs.ss.szozn_i[izn+1]; ++iszn) show_sz_flag[iszn] = 0;  // hides tri from image
              }
            }

            for (int iszn = 0; iszn < nsz; ++iszn) {
              if (show_sz_flag[iszn] == 0) continue;  // skip if already hidden

              it = hiddenSubZonesSet.find(iszn);
              if ( it!=hiddenSubZonesSet.end() ) show_sz_flag[iszn] = 0;  // hides tri from image
            }
          }

          rs.ss.getSubzoneBoundingBoxCenterAndDiagonal(bbox,xbb,diag,show_sz_flag);
        }
      }
      fprintf(fp,"{\n\"solver\":\"post\"");
      fprintf(fp,",\n\"b_surface_loaded\":%s",(b_surface_loaded?"true":"false"));
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

      const int nzn=rs.ss.zoneVec.size();
      // check if an empty post or not
      if (nzn>0) {
        fprintf(fp,",\n\"zoneNames\":[\n");
        for (int izn = 0; izn<nzn; ++izn) {
          if (izn == 0) fprintf(fp," \"%s\"",rs.ss.zoneVec[izn].getName().c_str());
          else fprintf(fp,",\n \"%s\"",rs.ss.zoneVec[izn].getName().c_str());
        }
        fprintf(fp,"\n]");
        fprintf(fp,",\n\"zoneSubzones\":[\n");
        assert((nzn+1) == rs.ss.szozn_i.getLength());
        for (int i = 0, limit = nzn+1; i < limit; ++i) {
          if (i == 0) fprintf(fp," \"%d\"",rs.ss.szozn_i[i]);
          else fprintf(fp,",\n \"%d\"",rs.ss.szozn_i[i]);
        }
        fprintf(fp,"\n]");

        if (rs.ss.selectedSubzoneVec.size()) {
          fprintf(fp,",\n\"selectedSubzones\":[\n");
          for (int i = 0, limit = rs.ss.selectedSubzoneVec.size(); i < limit; ++i) {
            if (i == 0) fprintf(fp," \"%d\"",rs.ss.selectedSubzoneVec[i]);
            else fprintf(fp,",\n \"%d\"",rs.ss.selectedSubzoneVec[i]);
          }
          fprintf(fp,"\n]");
        }
        if (rs.ss.b_update_hidden_subzones) {
          IntFlag sz_flag(rs.ss.nsz);
          sz_flag.setAll(0);

          for (vector<int>::iterator it=rs.ss.hiddenSubzoneVec.begin(); it!=rs.ss.hiddenSubzoneVec.end(); ++it) {
            sz_flag[*it] = 1;
          }

          // flag fully hidden zones
          IntFlag zn_flag(rs.ss.zoneVec.size());
          zn_flag.setAll(1);  // default to hidden
          for (int izn=0,nzn=rs.ss.zoneVec.size(); izn<nzn; ++izn) {
            for (int isz=rs.ss.szozn_i[izn],limit=rs.ss.szozn_i[izn+1]; isz<limit; ++isz) {
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
            vector<int>::iterator it=rs.ss.hiddenSubzoneVec.begin();
            while (it!=rs.ss.hiddenSubzoneVec.end()) {
              if (zn_flag[sz_flag[*it]]) it = rs.ss.hiddenSubzoneVec.erase(it);  // deletes and points to next el, so don't increment
              else ++it;
            }
          }

          if (rs.ss.hiddenSubzoneVec.size()) {
            fprintf(fp,",\n\"hiddenSubzones\":[\n");
            // cout << "wrote hiddenSubzones to JSON: ";
            bool first = true;
            for (int i = 0, limit = rs.ss.hiddenSubzoneVec.size(); i < limit; ++i) {
              if (first) {
                fprintf(fp," \"%d\"",rs.ss.hiddenSubzoneVec[i]);
                first = false;
              }
              else fprintf(fp,",\n \"%d\"",rs.ss.hiddenSubzoneVec[i]);
              // cout << ss.hiddenSubzoneVec[i] << ",";
            }
            // cout << endl;
            fprintf(fp,"\n]");
            rs.ss.hiddenSubzoneVec.clear();
          }

          rs.ss.b_update_hidden_subzones = false;
        }
        fprintf(fp,",\n\"zoneBoundaryConditions\":[\n");
        for (int izn = 0; izn<nzn; ++izn) {
          if (izn == 0) fprintf(fp," \"%s\"",rs.ss.zoneVec[izn].getBC().c_str());
          else fprintf(fp,",\n \"%s\"",rs.ss.zoneVec[izn].getBC().c_str());
        }
        fprintf(fp,"\n]");
        // periodic info
        const int npt = PeriodicData::periodicTransformVec.size();
        if (npt > 0) {
          fprintf(fp,",\n\"zonePeriodicBits\":[\n");
          for (int izn = 0; izn < nzn; ++izn) {
            if (izn == 0) fprintf(fp," \"%d\"",rs.ss.zoneVec[izn].getPeriodicBits());
            else fprintf(fp,",\n \"%d\"",rs.ss.zoneVec[izn].getPeriodicBits());
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
        fprintf(fp,",\n\"volumeVars\":{\n");
        fprintf(fp," \"D1\":[\n");
        map<const string, int8>::iterator it;
        bool first = true;
        // if velocity is available; register LIC as a D1 var
        it = rs.cvD2OffsetMap.find("u");
        bool bLIC = false;
        if (it != rs.cvD2OffsetMap.end())
          bLIC = true;
        else{
          it = rs.cvD2OffsetMap.find("U");
          if (it != rs.cvD2OffsetMap.end())
            bLIC = true;
        }
        for (it = rs.cvD1OffsetMap.begin(); it != rs.cvD1OffsetMap.end(); ++it) {
          if (first) {
            fprintf(fp,"  \"mesh\"");
            if (bLIC) fprintf(fp,",\n  \"LIC\"");
            first = false;
          }
          if (first) {
            fprintf(fp,"  \"%s\"",it->first.c_str());
            first = false;
          }
          else fprintf(fp,",\n  \"%s\"",it->first.c_str());
        }
        fprintf(fp,"\n ],\n");
        fprintf(fp," \"D2\":[\n");
        first = true;
        for (it = rs.cvD2OffsetMap.begin(); it != rs.cvD2OffsetMap.end(); ++it) {
          if (first) {
            fprintf(fp,"  \"mag(%s)\"",it->first.c_str());
            first = false;
          }
          else fprintf(fp,",\n  \"mag(%s)\"",it->first.c_str());
        }
        fprintf(fp,"\n  ]");
        fprintf(fp,"\n }");
      }  // if nzn>0
      fprintf(fp,"\n}\n");
      fclose(fp);

      // clear this vector so subsequent calls don't populate it
      rs.ss.selectedSubzoneVec.clear();
    }
    catch (int e) {
      WUI(WARN,"problem writing json");
    }
  }

  void writeAutoCompleteJson(FILE * fp) {
    fprintf(fp,",\n\"siData\": {");
    fprintf(fp,"\n  \"kbVersion\": \"%s\"",CTI::cti_docs_version);
    const string siData =
#include "post.siData"
    ;
    fprintf(fp,",\n  %s\n",siData.c_str());
    fprintf(fp,"}");
  }

  void processWriteParams(Param * param,const bool help) {
    if (help) {
      WUI(INFO,"a method to write parameter information to a Json file. This is meant to be used by the Cascade App to communicate with post"
      );
      return;
    }

    assert(param->getString(0) == "NAME");
    const string prefix = param->getString(1);

    // write the input file
    string filename = prefix+".in";
    string tmp_filename = "." + filename;
    FILE * fp = fopen(tmp_filename.c_str(),"w");
    assert(fp);

    if (!rs.isInit()) {
      WUI(WARN,"no simulation result available currently; please use \"RESTART <restart.mles> [result.sles]\" to initialize.");
    }
    else{
      string inFile;
      if (rs.filename_data=="")
        inFile = rs.filename_mesh;
      else
        inFile = rs.filename_data;

      char * param_buf = NULL;
      if(rs.readResultParams(param_buf,inFile)){
        assert(param_buf);
        fwrite(param_buf,sizeof(char), strlen(param_buf), fp);
        delete[] param_buf;
      }
    }

    fclose(fp);
    rename(tmp_filename.c_str(),filename.c_str());
  }

  void processWriteImage(Param * param,const bool help) {
    if (!rs.isInit()) {
      WUI(WARN,"no simulation result available currently; please use \"RESTART <restart.mles> [result.sles]\" to initialize.");
    }

    if (help) {
      WUI(INFO,"allows a user to export a png image of the current restart\n" <<
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

      scene.setRmax(rs.ss.getBoundingBoxRmax());
      double bbox[6];
      rs.ss.getBoundingBox(bbox);
      scene.setBoundingBbox(bbox[0],bbox[1],bbox[2],bbox[3],bbox[4],bbox[5]);
      double center[3];
      rs.ss.getBoundingBoxCenter(center);
      scene.setCenter(center[0],center[1],center[2]);

      // stuff for expression parsing...
      CtiRegister::clearCurrentData();         // reading in new data, so clear old
      CtiRegister::addCtiDataProducer(&scene); // this is necessary for registration to access varEvalCtiData (io) and funcEvalCtiData

      if (rs.filename_data == "") scene.enablePlanarData(rs.filename_mesh,"r_vv");
      else scene.enablePlanarData(rs.filename_mesh,rs.filename_data,"r_vv");

      for (int izn=0,nzn=rs.ss.zoneVec.size(); izn<nzn; ++izn) scene.addZoneName(rs.ss.zoneVec[izn].getName());
      scene.convertHiddenZoneNamesToIndices();
      
      // set hidden subzone info if present
      rs.ss.hiddenSubzoneVec.clear();
      if (scene.hasHiddenSubzones()) {
        int * sz_flag_tmp = new int[rs.ss.nsz];
        scene.getHiddenSubzones(sz_flag_tmp,rs.ss.nsz);
        for (int isz=0; isz<rs.ss.nsz; ++isz) {
          if (sz_flag_tmp[isz]) rs.ss.hiddenSubzoneVec.push_back(isz);
        }
        DELETE(sz_flag_tmp);
      }

      // treat zones as subzones so hiding gets processed by one routine
      // report back new fully hidden zones in json
      if (scene.hasHiddenZones()) {
        const int nzn = rs.ss.zoneVec.size();
        int * zn_flag_tmp = new int[nzn];
        scene.getHiddenZones(zn_flag_tmp,nzn);
        for (int izn=0; izn<nzn; ++izn) {
          if (zn_flag_tmp[izn]) {
            for (int isz=rs.ss.szozn_i[izn], limit=rs.ss.szozn_i[izn+1]; isz<limit; ++isz) rs.ss.hiddenSubzoneVec.push_back(isz);
          }
        }
        DELETE(zn_flag_tmp);
      }

      scene.initCanvas();

      if (scene.hasGeomPlane()) {
        scene.setGeomPlaneAsDataPlane();
        if (scene.blankDataPlane()) scene.addGeomPlaneAsBlankPlane(); // should be changed to back in scene/canvas
      }

      scene.addSurfaceTrisWithDataFromCim(rs.ss.xsp, rs.ss.spost, rs.ss.znost, rs.ss.szost, rs.ss.szozn_i, rs.ss.nst);
      // scene.addSurfaceTrisWithData(rs.ss.xsp, rs.ss.spost, rs.ss.szost, rs.ss.szozn_i, rs.ss.nst);  // if no data* passed then no data rendered

      if (scene.hasGeomPlane()) {
        scene.initDataPlaneCim(cvPlaneMap);
      }

      scene.addCimDataAndWriteImage();

      CtiRegister::removeLastCtiDataProducer();
    }
    catch (int e) {
      WUI(WARN,"problem writing image");
    }
  }

};


int main(int argc, char * argv[]) {

  try {

    CTI_Init(argc,argv,"post.in");

    {
      Post post;
      post.init();
      post.run();
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
