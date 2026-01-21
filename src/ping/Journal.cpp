#include "Journal.hpp"

namespace Journal {

  int step = -1;
  int iJournalLine = 0;
  int oJournalLine = 0;
  vector<string> iJournal;
  vector<string> oJournal;

  //return true to stop parsing additional lines
  //return false to indicate additional lines are needed to complete the param
  bool createMultiLineParamFromString(Param& param,const string& line) {
    bool inParamFlag = param.checkFlag();

    // look for a hash - this could be an entire line, or a comment at the end of a line...

    string::size_type hashPos = line.find_first_of("#",0);
    string trimmed_line = line.substr(0,hashPos);

    string::size_type contPos = trimmed_line.find_first_of("\\",0);
    if (contPos != string::npos)
      trimmed_line = line.substr(0,contPos);

    // start to parse...

    string delimiters = " \t=";
    string::size_type lastPos = trimmed_line.find_first_not_of(delimiters,0); // skip delimiters at beginning

    // check for an empty line...
    if ( !(lastPos == string::npos) ){

      // if a parameter is being continued, then we add this string's tokens to the
      // current parameter. Otherwise, we start a new parameter...

      string::size_type pos = trimmed_line.find_first_of(delimiters,lastPos); // find first "delimiter"

      if (param.checkFlag()) {
        param.tokens.push_back(trimmed_line.substr(lastPos,pos-lastPos));
      }
      else {
        param.setName(trimmed_line.substr(lastPos,pos-lastPos));
      }

      // and continue parsing...
      lastPos = trimmed_line.find_first_not_of(delimiters,pos);
      pos = trimmed_line.find_first_of(delimiters,lastPos);
      while (pos != string::npos || lastPos != string::npos) {
        param.tokens.push_back(trimmed_line.substr(lastPos,pos-lastPos));
        lastPos = trimmed_line.find_first_not_of(delimiters,pos); // skip delimiters
        pos = trimmed_line.find_first_of(delimiters,lastPos);
      }

      // set param flag if we found a "\" before at the end...
      if (contPos != string::npos)
        param.setFlag();
      else{
        param.clearFlag();
        return true;
      }
    }
    else{
      if (inParamFlag){
        cout << "Warning: Empty line found after multiline param" << endl;
        return true;
      }
    }
    return false;
  }

  void readInputJournal(const string &filename){
    std::ifstream infile(filename.c_str());
    readInputJournal(infile);
  }

  void readInputJournal(std::istream &infile){
    std::string line;
    while(std::getline(infile, line)){
      iJournal.push_back(line);
    }
    iJournalLine=0;

    if (iJournal.size() == 0) {
	CWARN("could not find input file or file is empty; skipping");
    }
  }



  bool isEndOfInputJournal(){
    if (iJournalLine >= int(iJournal.size()))
      return true;
    else
      return false;
  }



  bool parseParamFromJournal(Param &param){
    int i;
    for (i = iJournalLine; i < int(iJournal.size()); i++){
      if (createMultiLineParamFromString(param,iJournal[i])){
        iJournalLine = i+1;
        break;
      }
    }
    param.clearFlag();
    if (i>= int(iJournal.size())){
      iJournalLine = i;
      return false; //no more lines in journal
    }
    return true; //more lines in journal
  }

  void logCommandLine(const int &_step){
    while (oJournalLine < iJournalLine){
      logCommandLine(iJournal[oJournalLine++],_step);
    }
  }

  void logCommandLine(const string &line, const int &_step){
    if (_step > step){
      step = _step;
      char step_comment[20];
      sprintf(step_comment,"#STEP %i",step);
      oJournal.push_back("#-------------------");
      oJournal.push_back(step_comment);
      oJournal.push_back("#-------------------");
    }
    oJournal.push_back(line);
  }

  void dumpOutputJournal(){
    if (mpi_rank==0){
      cout << "====================== Journal of Commands =========================" << endl;
      for (int i=0; i< int(oJournal.size()); i++){
        cout << oJournal[i] << endl;
      }
      cout << "==================== End Journal of Commands =======================" << endl;
    }
  }

}
