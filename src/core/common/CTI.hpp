#ifndef CTI_HPP
#define CTI_HPP

/* *****************************************************************************
 *
 * CTI is a wrapper namespace for a collection of common functionality
 * including: MPI-related stuff, Parameter-related stuff, and misc utilities
 * that are nice to have.
 *
 * Include it in your solver/program as follows:
 *
 * include "CTI.hpp"
 * using namespace CTI;
 *
 ***************************************************************************** */

#ifdef FPE_TRAP

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <fenv.h>
#endif

#include "Common.hpp"
#include "MpiStuff.hpp"
#include "Params.hpp"
#include "Guardian.hpp"
#include "Logger.hpp"

using namespace MpiStuff;
using namespace Params;

namespace CTI {

  extern const char* cti_core_version;
  extern const char* cti_core_date ;
  extern const char* cti_docs_version;

  extern int cti_debug;    // various uses throughout code
  extern bool cti_verbose; // either true or false, controlled by VERBOSE on the input file line
  extern std::stringstream cti_journal;

  extern Logger * logger;

  void CTI_Init(int argc,char * argv[],const char *name);
  void CTI_Finalize();
  void CTI_Abort();
  void CTI_Dump_solver_info(ostream& ofile);

  // timestamp routines...
  void CTI_Make_timestamp(char timestamp[],const tm * time_tm);
  void CTI_Timestamp_journal();
}

class MpiGuardian : public Guardian {

protected:
  bool getFingerprintUNAndSN(string &username,string &serialNo);

public:
  void printLicenseError(const string &message);
  void initialize();
};

#endif
