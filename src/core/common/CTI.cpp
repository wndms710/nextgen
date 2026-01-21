#include "CTI.hpp"
#include "CtiRegister.hpp"

/*
void srand(const int seed) {
  std::cout << "CALLING srand with seed: " << seed << std::endl;
  std::srand(2);
  MpiStuff::MPI_Pause("OK");
}
*/

namespace CTI {

  // already in CTI.hpp...
  //using namespace MpiStuff;
  //using namespace Params;
  //using namespace MiscUtils;

  int cti_debug = -1;

#if defined(CTI_VERBOSE3) || defined(CTI_VERBOSE2)
  bool cti_verbose = true;
#else
  // default logging level
  bool cti_verbose = false;
#endif

  std::stringstream cti_journal;

  int argc_copy;
  char **argv_copy;

  MpiGuardian * guard;

  Logger * logger = NULL;

  void CTI_Init(int argc,char * argv[],const char * input_file) {

    // take a copy of the args...
    argc_copy = argc;
    argv_copy = argv;

    // initialize MPI...
    MPI_Init(&argc, &argv);

    // and the MpiStuff namespace...
    initMpiStuff();
    
    // output number formatting...
    cout << setprecision(8);
    
#ifdef WITH_SHM
    // the internode communicator links all (mpi_rank_shared == 0) cores, i.e. one per node...
    if (mpi_rank_shared == 0) {
      int my_ibuf[2] = { -mpi_size_shared, mpi_size_shared };
      int ibuf[2];
      MPI_Reduce(my_ibuf,ibuf,2,MPI_INT,MPI_MAX,0,mpi_comm_internode);
      if (mpi_rank == 0) {
	assert(mpi_rank_internode == 0);
	cout << argv[0] << " running on " << mpi_size << " cores, " << mpi_size_internode << " nodes";
	if (-ibuf[0]==ibuf[1]) cout << " (" << ibuf[1] << " cores per node)" << endl;
	else cout << " (" << -ibuf[0] << ":" << ibuf[1] << " cores per node)" << endl;
	cout << "version: " << cti_core_version << " built: " << cti_core_date << endl;
      }
    }
#else
    COUT1(argv[0] << " running on " << mpi_size << " cores");
    COUT1("version: " << cti_core_version << " built: " << cti_core_date);
#endif

// NFS file systems might need this
#ifdef SERIAL_IO
  if (mpi_rank == 0)
    cout << "Warning: shuffling data to rank0 for serial io!" << endl;
#endif

    /*
    // look for the special -i parameter.. and use that as the
    // input file if available ...
    // string input_file ;
    bool bDashIFound = false;
    int iargc = 1;
    {
      while (iargc < argc) {

        // look for a -i param, eventually you can look for other stuff too ...
        if ( (strlen(argv[iargc]) >= 2) && (argv[iargc][0] == '-') && (argv[iargc][1] == 'i') ) {

          // the filename should be the next argument ...
          if ( (iargc == argc-1) && (mpi_rank == 0))
            cout << " >> Found an input file argument, but no input file is specified ... " << endl ;

          assert ( iargc < argc-1 ) ;

          ++iargc ;

          input_file.assign(argv[iargc]) ;
	  bDashIFound = true;
          break;
        }

        ++iargc ;
      }//while
    }

    // if there is no input file specified, then use the default name ..
    if ( !bDashIFound ) iargc = 0 ; // rewind ...
    
    // and the Params namespace...
    // pass all the parameters that are AFTER the special single - keywords ..
    */
    
    initParams(argc,argv,input_file);
    
    // Verbosity: the VERBOSE param controls the verbosity. Valid uses:
    // VERBOSE 
    // VERBOSE true
    // VERBOSE false
    // VERBOSE 1
    
    cti_verbose = false; // default is false...
    if (Param * param = getParam("VERBOSE")) {
      if (param->tokens.size() > 0) {
	// expect a boolean...
	cti_verbose = param->getBool();
      }
      else {
	cti_verbose = true;
      }
    }

    if (checkParam("RANK_CHECK")) {
      MPI_Run_check();
    }
    
#ifndef NO_GUARDIAN
    // Check license after params initialization so that a user can specify an alternate
    // license file location in their input file
    if (mpi_rank== 0)
      cout << "checking license..." << endl;

    //Will throw an error if license is invalid
    guard = new MpiGuardian();
    guard->initialize();
    assert(guard->isValid());

    if (mpi_rank== 0){
      cout << "========================== license summary ============================\n" << 
	((guard->isFingerprintCheckEnabled())?" > License Verification Successful!":" > License Accepted with Warnings") <<
	"\n >   User:              " << guard->getLicUsr() <<
	"\n >   Organization:      " << guard->getLicOrg() << 
	"\n >   Serial Number:     " << guard->getLicSer() <<
	"\n >   Expiration Date:   " << guard->getLicExp() <<
	" >   Fingerprint Check: " << ((guard->isFingerprintCheckEnabled())?"Found Matching ID!":"Warning - Disabled!" ) <<
	"\n====================== end of license summary =========================" << endl;
    }
#endif

    // and journal...
    if (mpi_rank == 0) {
      CTI_Timestamp_journal();
      cti_journal << argv[0] << " " << mpi_size << " " << cti_core_version << " " << cti_core_date << endl;
    }

#ifdef FPE_TRAP
    feenableexcept(FE_DIVBYZERO|FE_INVALID) ;
#endif

    if(mpi_rank==0){
#ifdef WITH_SQLITE
      //solver_name,ncores,
      logger = new LoggerSqlite(argv_copy[0],mpi_size);
#else
      logger = new LoggerAscii(argv_copy[0],mpi_size); 
#endif
      logger->insertSolverAction(SIMULATION_INIT);
    }

    // and seed the rand with something different on each core...
    // NOTE: this is put here for the purposes of repeatability, however some mpi versions
    // can call their own "srand" with another unknown seed, disrupting this effort to seed the
    // rand function the same way each time. It is neccessary to delay this srand call until just 
    // before the rand routine is used (e.g. srand has been added to FlowSolver's init() routine to make lsp 
    // injection repeatable)...
    srand(mpi_rank+2); // 1 is reserved to reinitialized to its initial value, seems to be srand(0), so start at 2
    
  }

  void CTI_Finalize() {
    if (guard != NULL)
      delete guard;

    MPI_Barrier(MPI_COMM_WORLD);

    COUT1("CTI_Finalize()");

    if (mpi_rank==0 && logger!=NULL) {
      logger->insertSolverAction(SIMULATION_FINAL);
      delete logger;
      logger = NULL;
    }

    // do this no matter what verbosity level
    // report the final params and parameter usage for archival and debugging purposes...
    //if (cti_verbose) {
    if (mpi_rank == 0)
      dumpParams();
    dumpParamUsage();
    //}

    destroyMpiStuff();
    MPI_Finalize();
  }

  void CTI_Abort() {
    MPI_Abort(mpi_comm,-1);
  }

  void CTI_Dump_solver_info(ostream& ofile) {
    // argument info...
    ofile << "# solver: " << argv_copy[0] << " running on " << mpi_size << " cores" << endl;
    ofile << "# version: " << cti_core_version << " built: " << cti_core_date << endl;
  }

  void CTI_Make_timestamp(char timestamp[],const tm * time_tm) {
    strftime(timestamp,CTI_TIMESTAMP_SIZE,"%y%m%d:%H%M%S",time_tm);
  }

  void CTI_Timestamp_journal() {
    // only rank 0 should journal...
    assert(mpi_rank == 0);

    time_t raw_time = time(0);
    tm * now = gmtime(&raw_time);
    char timestamp[CTI_TIMESTAMP_SIZE];
    CTI_Make_timestamp(timestamp,now);
    cti_journal << timestamp << ":";
  }
}

#ifndef NO_GUARDIAN
//MpiGuardian Implementation...

bool MpiGuardian::getFingerprintUNAndSN(string &username,string &serialNo){

  using namespace MpiStuff;
  getUName(username);
  getSerialNumber(serialNo);

  //If unable to get username on one or more nodes,
  //skip system fingerprint check
  foundUsername=true;
  if (username == "N/A"){
    myLicError = 1;
  }
  MPI_Allreduce(&myLicError,&licError,1,MPI_INT,MPI_SUM,mpi_comm);
  if (licError > 0){
    foundUsername = false;
    licError = 0;
    myLicError = 0;
  }
  return foundUsername;
}

void MpiGuardian::printLicenseError(const string &message){
  using namespace MpiStuff;
  MPI_Allreduce(&myLicError,&licError,1,MPI_INT,MPI_SUM,mpi_comm);
  if (licError>0){
    if (mpi_rank==0){
      if (foundUsername){
        std::cout <<
        "\n\n************************ LICENSE ERROR *****************************\n > " <<
        message << "\n" <<
        "\n       The following resource username can be used to generate" <<
        "\n           a license at support.cascadetechnologies.com\n\n" <<
          "             > " << sysUN <<
        "\n\n********************* END OF LICENSE ERROR *************************\n" << std::endl;
      }
      else {
      std::cout <<
      "\n\n************************ LICENSE ERROR *****************************\n > " <<
      message <<
      "\n********************* END OF LICENSE ERROR *************************\n" << std::endl;
      }
    }
    throw(0);
  }
}

void MpiGuardian::initialize(){
  using namespace MpiStuff;
  //Set path for license file, order preference:
  //param, env var, compilation default
  bool foundPath = false;
  string filePath;
  using namespace Params;
  FOR_PARAM_MATCHING ("CTI_LICENSE"){
    if (param->tokens.size()>0){
      filePath = param->getString(0);

      if (mpi_rank==0)
        cout << " > Found param CTI_LICENSE=" << filePath << endl;

      foundPath = true;
    }
  }

  if (!foundPath){
    char* tPath;
    tPath = getenv("CTI_LICENSE");
    if (tPath!=NULL){
      filePath = tPath;

      if (mpi_rank==0)
        cout << " > Found env CTI_LICENSE=" << filePath << endl;;

      foundPath = true;
    }
    else{
      myLicError = 1;
    }
    MPI_Allreduce(&myLicError,&licError,1,MPI_INT,MPI_SUM,mpi_comm);
    if (licError > 0 && licError < mpi_size){
      myLicError = 1;
      printLicenseError("CTI_LICENSE environment variable not found on all ranks.\nSpecify CTI_LICENSE=<path/license.dat> parameter in your solver *.in file");
    }
    else
      myLicError = 0;
  }

  if (!foundPath){
    filePath = ctiLicenseFile;
    if (mpi_rank==0)
      cout << " > Using default path, CTI_LICENSE=" << filePath << endl;
  }

  Guardian::initialize(filePath);
}
#else
// dummy mpiguardian implementation...
bool MpiGuardian::getFingerprintUNAndSN(string &username,string &serialNo){ return true; }
void MpiGuardian::printLicenseError(const string &message){}
void MpiGuardian::initialize(){}
#endif
