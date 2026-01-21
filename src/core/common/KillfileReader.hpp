#ifndef _KILLFILE_READER_HPP_
#define _KILLFILE_READER_HPP_

#include "Params.hpp"
#include <queue>
#include "MpiStuff.hpp"

using namespace MpiStuff;

class KillfileReader {

  // ========================================================
  // This class implements the method getParam("killfile")
  // that returns a param for parsing to a calling process.
  //
  // In the calling process, it should be used something like
  // this:
  //
  //    KillfileReader kfr;
  //    ...
  //    while (Param * param = kfr.getParam("killsolver")) {
  //      // code to process param...
  //    }
  //
  // When hold is false (default) getParam will return NULL
  // if there is no killfile and all params have been consumed.
  // When hold is true, getParam will block until the killfile
  // is written.
  //
  // Note that this class no longer implicity processes "stop" or
  // "stop!". These are returned as param's for the calling
  // process to handle.
  //
  // Also note: if the user provides an empty killsolver file
  // (e.g. from a "touch killfile" command), the kfr will return
  // the param "stop" to allow the calling process to implement
  // the expected behavior in response to "stop".
  // ========================================================

private:

  enum KillFileInfo {
    KILL_FILE_NOT_FOUND,
    EMPTY_KILL_FILE_FOUND,
    KILL_FILE_FOUND,
  };

  bool hold;
  double hold_timeout_mins;
  queue<Param*> paramQueue;
  Param* current;

  // NOTE: journal stuff removed from KillfileReader. Journaling in general
  // will be MORE than just the params killfilereader touches, so handle in 
  // the process that owns the killfilereader...
  // for journal
  //bool journal;
  //string journal_filename;
  //list<string> strList; // opted to write to file immediately

public:

  KillfileReader() {
    hold = false;
    //journal = false;
    //journal_filename = "";
    hold_timeout_mins = -1;
    current = NULL;
  }

  ~KillfileReader() {
    if (current) delete current;
    while (!paramQueue.empty()) {
      delete paramQueue.front();
      paramQueue.pop();
    }
    //strList.clear();
  }

  void setHold(const bool hold) {
    // in some cases
    this->hold = hold;
  }

  bool getHold() {
    return this->hold;
  }

  /*
  void setJournal(const bool journal,const string& journal_filename = "killsurfer.journal") {
    // in some cases
    this->journal = journal;
    this->journal_filename = journal_filename;
    if (journal)
      remove(journal_filename.c_str());
  }
  */

  void setHoldWithTimeout(const bool hold, const double timeout_mins) {
    this->hold = hold;
    this->hold_timeout_mins = timeout_mins;
  }

  void push_back(const string& str) {

    cout << "[KillfileReader] ADDING: \"" << str << "\"" << endl;
    paramQueue.push(new Param());
    Params::createParamFromString(*paramQueue.back(),str);
  }

  Param * getParam(const string& filename) {

    // if the "current" param is not NULL, then this means it was returned in a
    // previous call to this routine, and we can delete it now...

    if (current) {
      //cout << "GOT VALID CURRENT! deleting" << endl;
      delete current;
      current = NULL;
      // if we are not in hold mode, do NOT try to read the killfile
      // again. Touching the filesystem has some performance impact,
      // so we avoid it when not in hold mode...
      if (!hold) {
        if (paramQueue.empty()) return NULL;
        else {
          current = paramQueue.front();
          paramQueue.pop();
          return current;
        }
      }
    }

    // now parse the passed filename and add params to the paramQueue. We do
    // this part first even if there are current entries in the paramQueue because
    // we want the reading to be greedy, especially in hold mode...

    int status;
    vector<string> lines;
    if (mpi_rank == 0) {

      KillFileInfo kfi = readAndRemoveKillFile(lines,filename);

      // if:
      // 1. we did not find a file, AND
      // 2. there is no existing queue of params to return, AND
      // 3. we are in hold mode,
      // then loop until we get a killfile...

      double mpi_wtime_0 = MPI_Wtime();
      while ((kfi == KILL_FILE_NOT_FOUND)&&(paramQueue.empty())&&(hold)) {
        //cout << "got kfi: " << kfi << endl;

        //Exit the hold loop if a timeout is set and exceed.
        //Note, this only changes the value of hold on rank 0 but
        //this should lead to a collective status==-1 and a return NULL
        if ( hold_timeout_mins>=0 && (MPI_Wtime()-mpi_wtime_0  > hold_timeout_mins*60) ){
           cout << "Reached timeout limit of " << hold_timeout_mins << " min(s), hold off" << endl;
           //kfi == KILL_FILE_NOT_FOUND, paramQueue.empty() == true
           hold = false;
           break;
        }
        usleep(500000);
        kfi = readAndRemoveKillFile(lines,filename);
      }

      if (kfi == KILL_FILE_NOT_FOUND) {
        if (paramQueue.empty()) {
          assert(!hold);
          // this is most common: no file was found, no params to process, not in hold mode...
          status = -1;
        }
        else {
          // we did not find a file, but there are params to return...
          status = -2;
        }
      }
      else if (kfi == EMPTY_KILL_FILE_FOUND) {
        // this is the same as a stop.
        lines.push_back("stop");
        status = 5; // 4 + 1
      }
      else {
        assert(kfi == KILL_FILE_FOUND);
        assert(!lines.empty());
        status = 0;
        for (int i = 0, i_end=lines.size(); i < i_end; ++i) status += lines[i].length()+1;
      }

    }

    MPI_Bcast(&status, 1, MPI_INT, 0, mpi_comm);

    // ----------------------------------------------------------
    // for a negative status, we can take action immediately...
    // ----------------------------------------------------------

    if (status == -1) {
      // this is most common: no file was found, no params to process, not in hold mode...
      return NULL;
    }
    else if (status == -2) {
      // we did not find a file, but there are still params in the queue...
      current = paramQueue.front();
      paramQueue.pop();
      return current;
    }

    // ----------------------------------------------------------
    // for a positive status, we need to exchange the new strings
    // that were read and turn them into parameters on ALL ranks.
    // ----------------------------------------------------------

    assert(status > 0);
    char * cbuf = new char[status];
    if (mpi_rank == 0) {
      int offset = 0;
      for (int i = 0,i_end=lines.size(); i < i_end; ++i) {
        for (int j = 0, j_end=lines[i].length(); j < j_end; ++j) cbuf[offset++] = lines[i].at(j);
        cbuf[offset++] = '\0';
      }
      assert(offset == status);
    }
    MPI_Bcast(cbuf,status,MPI_CHAR,0,mpi_comm);

    // now everyone has the cbuf. Use it to queue params...
    int last = -1;
    while (last < status-1) {
      const int first = last+1;
      last = first + strlen(cbuf+first);
      // take a look (use the last rank to confirm parallelism)...
      //if (mpi_rank == mpi_size-1) {
      //  cout << "got status: " << status << " first: " << first << " and last: " << last << endl;
      //  cout << "current string is \"" << (cbuf+first) << "\"" << endl;
      //}
      // ----------------------------------------------------------------------------
      // at this point, we have the param line in the char buffer. Look for special
      // behaviors first...
      //
      // special keywords are:
      //
      // "hold on"  : turns hold mode on
      // "hold off" : turns hold mode off
      // "flush"    : flushes unprocessed params - maybe just "const" params
      //
      // Other lines get interpreted as standard params
      // ----------------------------------------------------------------------------
      if (strcmp(cbuf+first,"hold on") == 0) {
        if (mpi_rank == 0) cout << "hold on" << endl;
        hold = true;
      }
      else if (strcmp(cbuf+first,"hold off") == 0) {
        if (mpi_rank == 0) cout << "hold off" << endl;
        hold = false;
      }
      else if (strcmp(cbuf+first,"flush") == 0) {
        if (mpi_rank == 0) cout << "flush" << endl;
        // clear the queue...
        while (!paramQueue.empty()) {
          delete paramQueue.front();
          paramQueue.pop();
        }
      }
      else {
        // assum this is a new param...
        paramQueue.push(new Param());
        Params::createParamFromString(*paramQueue.back(),cbuf+first);
        // look at the param...
        //if (mpi_rank == mpi_size-1) {
        //  Param * param = paramQueue.back();
        //  cout << "param->getName(): \"" << param->getName() << "\"" << endl;
        //  for (int i = 0; i < param->size(); ++i)
        //    cout << " > param->getString(i): \"" <<  param->getString(i) << "\"" << endl;
        //}
        //MPI_Pause("OK");
      }
    }
    delete[] cbuf;

    if (paramQueue.empty()) {
      if (hold) {
        // call this routine recursively. Note that this does create a solver
        // stack memory vulnerability where a user repeatedly submits "hold on".
        // This user is a loser and deserves the consequences.
        return getParam(filename);
      }
      else {
        return NULL;
      }
    }
    else {
      current = paramQueue.front();
      paramQueue.pop();
      return current;
    }

  }

private:

  KillFileInfo readAndRemoveKillFile(vector<string>& lines,const string& filename) {

    // rank 0 only...
    assert(mpi_rank == 0);

    //cout << " ==================== TOUCHING FILE =======================" << endl;

    ifstream ifile;
    ifile.open(filename.c_str());
    if (ifile.is_open()) {
      parseParamFile(lines,ifile);
      ifile.close();
      remove(filename.c_str());
      if (lines.empty()) return(EMPTY_KILL_FILE_FOUND);
      else return(KILL_FILE_FOUND);
    }
    else {
      return(KILL_FILE_NOT_FOUND);
    }

  }

  void parseParamFile(vector<string>& lines,const string& filename) {

    ifstream ifile;
    ifile.open(filename.c_str());
    if (ifile.is_open()) {
      parseParamFile(lines,ifile);
      ifile.close();
    }

  }

  void parseParamFile(vector<string>& lines,ifstream& ifile) {

    // ====================================================
    // robust param file parser. Returns the passed, opened
    // param file parsed as a vector of strings.
    // - Continuations are all processed to produce single lines,
    // - starting and trailing whitespace is eliminated
    // - comments and empty lines are removed (# comment format)
    // ====================================================

    bool multiline = false;
    while (1) {

      string line;
      getline(ifile,line);

      if (ifile.fail()) break;

      //cout << "\n\nWorking on line \"" << line << "\"..." << endl;

      // look for an empty string...
      const std::size_t start = line.find_first_not_of(" \t\f\v\n\r");
      if (start == string::npos) {
        //cout << "empty" << endl;
        if (multiline) {
          //cout << "Warning: param file format has continuation to empty line. Ending line continuation." << endl;
          const std::size_t end = lines.back().find_last_not_of(" \t\f\v\n\r");
          if (end < lines.back().length()-1) lines.back().resize(end+1);
          multiline = false;
        }
        continue;
      }
      // look for full comment line...
      const std::size_t hash = line.find_first_of("#");
      if (hash == start) {
        //cout << "full comment" << endl;
        if (multiline) {
          //cout << "Warning: param file format has continuation to comment line. Ending line continuation." << endl;
          const std::size_t end = lines.back().find_last_not_of(" \t\f\v\n\r");
          if (end < lines.back().length()-1) lines.back().resize(end+1);
          multiline = false;
        }
        continue;
      }
      else if (hash != string::npos) {
        //cout << "trailing comment: truncating line at hash: " << hash << endl;
        line.resize(hash);
      }

      // set end to eliminate trailing whitespaces...
      std::size_t end = line.find_last_not_of(" \t\f\v\n\r");

      /*
      if (end < line.length()-1) {
      cout << "looks like there is trailing whitespace" << endl;
    }
    */

    bool activate_multiline = false;
    if (line.at(end) == '\\') {
      //cout << "LOOKS LIKE MULTILINE " << start << " " << end << endl;
      if (start == end) continue; // must be just a "\" in a single line
      activate_multiline = true;
      --end; // do not add continuation character "\" to lines
    }

    assert(end != string::npos);
    assert(start <= end);

    if (multiline) {
      assert(!lines.empty());
      lines.back().append(line,0,end+1);
    }
    else {
      lines.push_back(line.substr(start,end-start+1));
    }

    multiline = activate_multiline;

    // report lines so far...
    /*
    cout << "XXXXXXXXXX lines XXXXXXXXXXX" << endl;
    for (int i = 0; i < lines.size(); ++i)
    cout << "\"" << lines[i] << "\"" << endl;
    cout << "XXXXXX end of lines XXXXXXXX" << endl;
    getchar();
    */
  }

  // if we came out with multiline, then one last correction...
  if (multiline) {
    //cout << "Warning: param file format has continuation to EOF. Ending line continuation." << endl;
    assert(!lines.empty());
    const std::size_t end = lines.back().find_last_not_of(" \t\f\v\n\r");
    if (end+1 < lines.back().length()) lines.back().resize(end+1);
    multiline = false;
  }

  // report final lines...
  /*
  cout << "XXXXXXXX final lines XXXXXXX" << endl;
  for (int i = 0; i < lines.size(); ++i)
  cout << "\"" << lines[i] << "\"" << endl;
  cout << "XXXXXX end of lines XXXXXXXX" << endl;
  getchar();
  */

}

};

#endif
