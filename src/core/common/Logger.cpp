#include "Logger.hpp"
#include "MpiStuff.hpp"

namespace LoggerData {
  static const int nColumns=8;
  static const string columnNames[nColumns] = {"timestamp",
                                        "hash_id",
                                        "username",
                                        "solver_name",
                                        "ncores",
                                        "rank0_host:cwd_path",
                                        "killfile",
                                        "action"};
  static const string columnTypes[nColumns] = {"TEXT NOT NULL",
                                               "TEXT NOT NULL",
                                               "TEXT NOT NULL",
                                               "TEXT NOT NULL",
                                               "INT NOT NULL",
                                               "TEXT NOT NULL",
                                               "TEXT NOT NULL",
                                               "INT NOT NULL"};
}

Logger::Logger(const string _solver_name, const int _ncores) : db_filename(""), bEnableLogging(false), solver_name(_solver_name), ncores(_ncores), cwd_path(""), username(""), hash_id("0000000000000000000000000000"), killfile("killsolver"), rank0Host("") {
 
  //initial database field entries
  struct passwd *pw = getpwuid(getuid());
  if (pw){
    username = pw->pw_name;
  }
  else{
    username = "N/A";
  }

  //  working directory
  char cwd[1024];
  if (getcwd(cwd,sizeof(cwd))){
    cwd_path = cwd;
  }
  else{
    cout << " > Warning: WebUI run logging cannot determine current working directory" << endl;
  }

  //Logger is currently run from rank0, so this will be the rank0 host though not enforced here...
  char hostnameChar[MPI_MAX_PROCESSOR_NAME];
  for (int ic = 0; ic<MPI_MAX_PROCESSOR_NAME; ++ic)
    hostnameChar[ic] = '\0';

  int len = MPI_MAX_PROCESSOR_NAME;
  MPI_Get_processor_name(hostnameChar,&len);
  if (len < MPI_MAX_PROCESSOR_NAME && hostnameChar[len] == '\0'){
    rank0Host = hostnameChar;
  }
  else{
    cout << " > Warning: WebUI run logging cannot determine rank0 host name" << endl;
  }

  //Set path to logger db and name
  //  1. check for CTI_DB env var
  //  2. use default: (a) $HOME env+cti_log (b) get home dir
 
  char* tPath = NULL;
  tPath = getenv("CTI_RDB");
  if (!tPath){
    tPath = getenv("HOME"); //home dir
    if (!tPath){
      if (pw){ 
        tPath = pw->pw_dir; //home dir
      }
      else{
        cout << " > Warning: failed to find path to home dir for WebUI run logging. " << endl;
        cout << " >          Specify alternate path with env var CTI_RDB=/<path>/cti_rdb.db" << endl;
        cout << " >          Run logging DISABLED" << endl;
        return;
      }
    }
    db_filename = tPath;
    db_filename += "/cti_rdb.db";
    cout << " > Using home directory for WebUI run log ... " << db_filename << endl;
  } 
  else{
    db_filename = tPath;
    cout << " > Using env CTI_RDB for WebUI run log ... " << db_filename << endl;
  }

}

string Logger::getTimestamp(){
  if (!bEnableLogging)
    return "";

  time_t now;
  time(&now);
  struct tm * timeUTC;
  timeUTC = gmtime(&now);
  char buf[80];
  strftime(buf,80,"%F_%T",timeUTC);
  return buf;
}

void Logger::setKillFilename(const string& _killfilename){
  killfile = _killfilename;
}
void Logger::setHashId(const string& _hash_id){
  hash_id = _hash_id;
}


LoggerAscii::LoggerAscii(const string _solver_name, const int _ncores) : Logger(_solver_name,_ncores), fp(NULL) {
  //Check for existance of log, if none present perform any setup needed,
  //output specific error messages within checkRunDb
  int error_code = enableLogging();
  if (error_code){
    cout << " >          WebUI run logging DISABLED" << endl;
    bEnableLogging = false;
  }
  else{
    bEnableLogging = true;
  }

}

int LoggerAscii::enableLogging(){

  string keyStr="CascadeRunLog";
  //check if log file exists
  fp = fopen(db_filename.c_str(),"r");
  if (!fp){
    //file doesn't exist, try to create one
    cout << " > Logger: using new db file, creating header" << endl;
    fp = fopen(db_filename.c_str(),"w");
    if (!fp){
      cout << " > Warning: unable to create new WebUI log file" << endl;
      return 1;
    }
    fwrite(keyStr.c_str(),sizeof(char),keyStr.length(),fp);
    if (ferror(fp)){
      cout << " > Warning: unable to write to WebUI log file" << endl;
      fclose(fp);
      return 1;
    }
    //write column headers
    string header="\n";
    for (int i=0;i<LoggerData::nColumns-1;++i)
      header+=LoggerData::columnNames[i]+",";
    header+=LoggerData::columnNames[LoggerData::nColumns-1]+"\n";
    fwrite(header.c_str(),sizeof(char),header.length(),fp);
    if (ferror(fp)){
      cout << " > Warning: unable to write to WebUI log file" << endl;
      fclose(fp);
      return 1;
    }
    fclose(fp);
  }
  else{
    //check file starts as expected
    char buf[14];
    memset(buf,'\0',sizeof(char)*14);
    fread(buf,sizeof(char),13,fp);
    fclose(fp);
    if (strcmp(keyStr.c_str(),buf)!=0){
      cout << " > Warning: WebUI log file does not start as expected" << endl;
      //cout << " >          Expecting: " << keyStr << endl;
      //cout << " >          Found:     " << buf << endl;
      return 1;
    }
  }
  return 0;
}

int LoggerAscii::insertSolverAction(LoggerAction _action){
  if (!bEnableLogging)
    return 0;

  //Add row to ascii log fp with specified action
  /* Open file for (append) writing */
  fp = fopen(db_filename.c_str(), "a");
  if (!fp){
    cout << " > Warning: unable to append to WebUI log file" << endl;
    cout << " >          WebUI run logging DISABLED" << endl;
    bEnableLogging = false;
    return 1;
  }
  else{
    /* Create CSV entry*/
    //timestamp,hash_id,username,solver_name,ncores,rank0Host:cwd_path,killfile,action
    stringstream ss;
    ss << getTimestamp() << "," << hash_id << "," << username << ","; 
    ss << solver_name << "," << ncores << ",";
    if (rank0Host!=""){
      //timestamp,hash_id,username,solver_name,ncores,cwd_path,killfile,action
      ss << rank0Host << ":";
    }
    ss << cwd_path << "," << killfile << "," << _action << "\n";

    string csvStr = ss.str(); 

    /* Write CSV entry */
    fwrite(csvStr.c_str(),sizeof(char),csvStr.length(),fp);
    if (ferror(fp)){
      cout << " > Warning: unable to write to WebUI log file" << endl;
      cout << " >          WebUI run logging DISABLED" << endl;
      if (rank0Host!="")
        cout << " >          rank_0_Host ... " << rank0Host << endl;

      bEnableLogging = false;
      fclose(fp);
      return 1;
    }
    fclose(fp);    

    if (rank0Host!="")
      cout << " > Appended to WebUI log successfully, rank_0_Host ... " << rank0Host << endl;
    else 
      cout << " > Appended to WebUI log successfully" << endl;

    return 0;
  }
}


#ifdef WITH_SQLITE
static int logger_callback(void *data, int argc, char **argv, char **azColName){
   int i;
   for(i=0; i<argc; i++){
      printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
   }
   printf("\n");
   return 0;
}


LoggerSqlite::LoggerSqlite(const string _solver_name, const int _ncores) : db(NULL), Logger(_solver_name,_ncores) {

  //Check for existance of db, if none present create needed schema
  //enable logging if table is ready to accept inserts
  int error_code = enableLogging();
  if (error_code){
    cout << " >          WebUI run logging DISABLED" << endl;
    bEnableLogging = false;
  }
  else{
    bEnableLogging = true;
  }

};


int LoggerSqlite::enableLogging(){

  // Check if db file can be opened
  // Create the db schema if it doesn't exist
  int error_code = sqlite3_open(db_filename.c_str(), &db);
  if( error_code ){
    cout << " > Warning: can't open database: " << sqlite3_errmsg(db) << endl;
    return error_code;
  }
  else{
    //Database was opened successfully
    //Check for table existance
    string sql = "SELECT name FROM sqlite_master WHERE type='table' AND name='solver_history'";
    sqlite3_stmt *stmt; //binary sql query
    error_code = sqlite3_prepare_v2(db,sql.c_str(),-1,&stmt, 0);
    if (error_code == SQLITE_OK && stmt!=NULL){

      error_code = sqlite3_step(stmt);
      bool bCreateTable = true;
      if (error_code == SQLITE_ROW){
        //table exists
        sqlite3_clear_bindings(stmt);
        sqlite3_reset(stmt);
        bCreateTable = false;
      }
      error_code = sqlite3_finalize(stmt);
  
      if (bCreateTable){
        //create table
        cout << " > Logger: using new db file, creating table 'solver_history'" << endl;

        sql = "CREATE TABLE solver_history(";
        // "timestamp       TEXT NOT NULL," 
        // .....
        // "action            INT NOT NULL);";
        for (int i=0;i<LoggerData::nColumns-1;++i){
          sql+= LoggerData::columnNames[i]+"     ";
          sql+= LoggerData::columnTypes[i]+", ";
        }
        sql+=LoggerData::columnNames[LoggerData::nColumns-1]+"     ";
        sql+=LoggerData::columnTypes[LoggerData::nColumns-1]+");";


        /* Execute SQL statement */
        char *zErrMsg;
        error_code = sqlite3_exec(db, sql.c_str(), logger_callback, 0, &zErrMsg);
        if( error_code != SQLITE_OK ){
          cout << "Warning, " <<  zErrMsg << endl;
          sqlite3_free(zErrMsg);
          return error_code;
        }else{
          cout << " > Table created successfully" << endl;
        }
      }
    }
    else{
      cout << " > Warning, problem checking cti_rdb tables" << endl;
      return 1;
    }
    //close for now
    sqlite3_close(db);
    return 0;
  }

}


int LoggerSqlite::insertSolverAction(LoggerAction _action){
  if (!bEnableLogging)
    return 0;

  //Add row to sqlite db with specified action
  /* Open database */
  int error_code = sqlite3_open(db_filename.c_str(), &db);
  if( error_code ){
    cout << " > Warning: can't open WebUI log database: " << sqlite3_errmsg(db) << endl;
    return(0);
  }else{
    cout << " > Opened WebUI log database" << endl;

    /* Create SQL statement */
    stringstream ss;
    ss << "INSERT INTO solver_history("; 
    for (int i=0;i<LoggerData::nColumns-1;++i){
      ss << LoggerData::columnNames[i] << ",";
    }
    ss << LoggerData::columnNames[LoggerData::nColumns-1] << ") VALUES (";

    //timestamp,hash_id,username,solver_name,ncores,rank0Host:cwd_path,killfile,action
    ss << "'" << getTimestamp() << "', '" << hash_id << "', '" << username;
    ss << "', '" << solver_name << "'," << ncores << ", '";
    if (rank0Host!=""){
      ss << rank0Host << ":";
    }
    ss << cwd_path << "', '" << killfile << "', " << _action  << "); ";
    string sql = ss.str(); 

    /* Execute SQL statement */
    char *zErrMsg;
    error_code = sqlite3_exec(db, sql.c_str(), logger_callback, 0, &zErrMsg);
    if( error_code != SQLITE_OK ){
      cout << " > Warning: can't insert row into WebUI log database: " << zErrMsg << endl;
      sqlite3_free(zErrMsg);
    }
    else {
      cout << " > Row inserted successfully into WebUI log database" << endl;
    }
    sqlite3_close(db);    
  }
}
#endif

