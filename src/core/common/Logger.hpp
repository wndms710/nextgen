#ifndef LOGGER_HPP
#define LOGGER_HPP

// ======================================================================
// Logger.hpp
//
// This class writes information about launched simulations to
// a centralized log file for use by the WebUI.  When connecting to
// a resource, the WebUI will check this log and be able to connect
// users to their running simulations without additional manual
// entry of run directories, solver names, etc.
//
// The log is written to a user's home directory by default.  For
// systems where compute nodes cannot access the home directory,
// (ex: Titan) an alternate path can be set with the environment 
// variable CTI_RDB.  When using an alternate path it must be configured 
// once in the resources tab of the WebUI.  If any problems are 
// encountered initializing the log file writing it is disabled for the
// duration of a simulation.
//
// Logger writes an ASCII formated csv file by default with a row 
// appended for each solver action.  There is basic checking the log
// file begins as expected but there is no error checking of the
// log file contents prior to writing or protection against concurrent 
// writes.
//
// Alternately, logger can write to an sqlite3 database file.  This 
// guarantees inserted rows are recorded in their entirety or not at
// all. It requires building the code against libsqlite3 (found on many
// systems) or including sqlite3.c and sqlite3.h in this repo (ommitted
// because sqlite3.c is a ~6MB source file, also a C compiler must be
// specified in Makefile.in).  
//
// History:
//  author: David Philips - Feb 2017
// ======================================================================


#include "Common.hpp"

#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include <time.h>

enum LoggerAction{
 SIMULATION_INIT,
 SIMULATION_RUN,
 SIMULATION_FINAL
};

class Logger {

protected:

  //path to db file
  string db_filename;
  bool bEnableLogging;

  //data
  string solver_name;
  int    ncores;
  string cwd_path;
  string username;
  string hash_id;
  string killfile;
  string rank0Host;

public:
  
  Logger(const string _solver_name, const int _ncores);
  virtual ~Logger() {};

  virtual int enableLogging() = 0; //return 0 for ok, 1 for error
  virtual int insertSolverAction(LoggerAction _action) = 0;
  
  string getTimestamp();
  void setKillFilename(const string& _killfilename);
  void setHashId(const string& _hash_id);
};

class LoggerAscii: public Logger {

private:
  //file pointer
  FILE *fp;

public:
   LoggerAscii (const string _solver_name, const int _ncores);
 
   int enableLogging();
   int insertSolverAction(LoggerAction _action);
  
};

#ifdef WITH_SQLITE
#undef Z
#include <sqlite3.h>

class LoggerSqlite: public Logger {

private:
  //database pointer
  sqlite3 *db;

public:
  LoggerSqlite(const string _solver_name, const int _ncores);

  int enableLogging();
  int insertSolverAction(LoggerAction _action);

};
#endif

#endif
