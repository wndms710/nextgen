#ifndef BASICPOSTPRO_HPP
#define BASICPOSTPRO_HPP

class BasicPostpro : public StaticSolver { 
public:

  //double * phi;

  BasicPostpro() {
  
    // if you customize this instance, you can register your data here.. 
  
    //phi = NULL; registerCvData(phi, "PHI", READWRITE_DATA);

  }

  ~BasicPostpro() { 

    // you are responsible for deallocating any customized data you allocated

    //DELETE(phi);
  } 

  virtual void initData() {
    
    // allocate any registered data (or unregistered data) here

    // for the test example, we are going to allocate a registered variable 
    // with (first) level ghosts.

    //phi = new double[ncv_g];

  } 

  virtual void initMin() { 

    // no need to build the extended face instrastructure here..

    StaticSolver::init(INIT_COMPACT_FACES | INIT_CV_GRAD);
    StaticSolver::initFromParams();

    if (mpi_rank == 0)
      logger->setKillFilename("killcharles");

  }

  void syncPostState() { 

    // default behavior is operating on registered data that does not 
    // any ghost information, so this routine is empty.  however, if you
    // have registered data with ghost information, you would need to 
    // updateCvData ghosts

    //updateCvData(phi);
  }

  void runPost() { 

    // for the solvers running inside of post mode- solutions are not 
    // advanced and instead we either operate on the result file pair 
    // that was loaded or process a series of snapshots

    // tell the central log a solver is running...
    
    if (mpi_rank==0) logger->insertSolverAction(SIMULATION_RUN);

    // we need to register particles here to prevent StaticSolver from trying to init unsized particle data

    bool b_particles = false;
    FOR_ALL_PARAM {
      if ( ( (param->name == "REGISTER_LP_R1") || (param->name == "REGISTER_LP_DN") ) && param->touch() ) {
        b_particles = true;

        if (checkParam("SNAPSHOT")) {
          for (int iarg = 0; iarg < param->size(); ++iarg)
            registerLpDN( param->getString(iarg) );
        }

      } else if ( ( (param->name == "REGISTER_LP_R2") || (param->name == "REGISTER_LP_DN3") ) && param->touch() ) {
        b_particles = true;

        if (checkParam("SNAPSHOT")) {
          for (int iarg = 0; iarg < param->size(); ++iarg)
            registerLpDN3( param->getString(iarg) );
        }

      } 
    }

    // needed by app to control contex menu
    if (b_particles)
      added_data_fields.insert("particles");
        
    initialHook();

    if ( Param* param = getParam("SNAPSHOT")) { 

      if ( (param->getString(0) != "NAME") || (param->getString(2) != "RANGE")) { 
        CERR( " Invalid snapshot syntax; SNAPSHOT NAME <prefix> RANGE <start> <inc> <last>");
      }

      const string prefix         = param->getString(1);
      const int snapshot_first    = param->getInt(3);
      const int snapshot_inc      = param->getInt(4);
      const int snapshot_last     = param->getInt(5);

      // if we are in snapshot mode, we also need to check for the intended behavior 
      // for statistical averaging behavior if stats are already present in the snapshot
      // files.  we can either take ensemble avg data over the snapshots (default) or 
      // we can allow the user to step through their averages.

      const bool b_overwrite_stats = getBoolParam("OVERWRITE_STATS",true);

      // reset the stats so that we can take ensemble statistics from the available data
      
      double snapshot_stats_wgt = 1.0;   // ensemble weighting

      if ( b_overwrite_stats) { 

        // we need to clear the read bits for any stats data and reset the stats ..
        // so our own ensemble averaging data is not clobbered by the snapshot files
      
        CtiRegister::resetStats();
        CtiRegister::clearStatsBits(READ_DATA);

      } else { 

        snapshot_stats_wgt = 0.0;  // to prevent updating the stats info.. 

      }

      // loop over the snapshots and give access to the step processing .. 

      int done = 0;
      int step = 0;
      double time = 0.0;

      for (int snapshot = snapshot_first; snapshot <= snapshot_last; snapshot += snapshot_inc) { 

        CtiRegister::clearCurrentData(); 

        char filename[128];
        MiscUtils::buildIndexedFilename(filename,prefix.c_str(),snapshot,"sles");
        if ( snapshot != snapshot_first) RestartHashUtilities::clearSlesHash();
        readData(filename);

        // we should have a valid step and time at this stage from the snapshot file
        // we will report to stdout just to confirm

        if ( mpi_rank == 0) 
          cout << " > processing snapshot: " << filename << endl; 

        // solver specific sync of the data and population of primitive data fields..

        syncPostState();
        
        temporalHook();
        
        // interactive viewing of snapshots

        if ((snapshot > snapshot_first)&&(checkParam("INTERACTIVE")))
          kfr.setHold(true);

        // process step is called to activate the probing, imaging, etc.
        // if the step and time are registered, we will use their values in the output

        if ( CtiRegister::CtiData* s = CtiRegister::getRegisteredCtiData("step",false)) 
          step = s->i();
       
        if ( CtiRegister::CtiData* t = CtiRegister::getRegisteredCtiData("time",false))
          time = t->d();
        
        done = processStep(step,time,snapshot_stats_wgt,false);

        if ( done == -2) 
          return;
        else if ( done == -1) 
          break;

      }

      flushProbes();
      flushImages();
      finalHook();

      if ( (done >= -1) && (checkParam("RESULT")) ) { 

        string result_name = getStringParam("RESULT");
        MiscUtils::mkdir_for_file(result_name.c_str());

        if ( mpi_rank == 0) 
          cout << " Writing solution file \"" << result_name << "\"..." << endl;
        CtiRegister::writeData(result_name);

      }

    } else { 

      // run the post operations on a single piece of data that has been read in... 

      CtiRegister::clearCurrentData();

      syncPostState();

      temporalHook();

      // set to 1 in order to be able to write result (step = 0 is skipped) if step is not registered
      int step = 1;
      if ( CtiRegister::CtiData* s = CtiRegister::getRegisteredCtiData("step",false)) 
        step = s->i();
      
      double time = 0.0;
      if ( CtiRegister::CtiData* t = CtiRegister::getRegisteredCtiData("time",false))
        time = t->d();

      processStep(step,time,0.0,false);

      flushProbes();
      flushImages();
      finalHook();

      if ( checkParam("RESULT")) { 

        string result_name = getStringParam("RESULT");
        MiscUtils::mkdir_for_file(result_name.c_str());

        if ( mpi_rank == 0) 
          cout << " Writing solution file \"" << result_name << "\"..." << endl;
        CtiRegister::writeData(result_name);
      }

    }

  }

  // additional user defined hooks 

  virtual void initialHook() { 

    // called before the snapshot loop is performed, for instance

    //FOR_ICV_G phi[icv] = 0.0;

  }

  virtual void temporalHook() {
  
    // called either once or after each snapshot is read in.. 

    //CtiRegister::CtiData* p_data = CtiRegister::getCtiData("p"); // assuming it was registered from the input file
    //assert( p_data != NULL);
    //assert( p_data->getTopology() == CV_DATA);
    //FOR_ICV phi[icv] += p_data->dn(icv);
    
  }
  
  virtual void finalHook() {
  
    // this is called prior to exit
  
  }

};

#endif
