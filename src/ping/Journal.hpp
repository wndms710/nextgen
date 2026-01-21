#ifndef _JOURNAL_HPP_
#define _JOURNAL_HPP_

#include "CTI.hpp"
using namespace CTI;

namespace Journal {
 
  bool createMultiLineParamFromString(Param& param,const string& line);
  void readInputJournal(std::istream &infile);
  void readInputJournal(const string &filename);
  bool isEndOfInputJournal();
  bool parseParamFromJournal(Param &param);

  void logParams(const int &step, const string message);
  void logParams(const vector<string> inputs, const int &step, const string message);

  void logCommandLine(const int &_step);
  void logCommandLine(const string &line, const int &_step);
  void dumpOutputJournal();
}

#endif 
