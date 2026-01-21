#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <sstream>

using namespace std;
#include "Defs.hpp"

class Param {

private:
  int count;
  bool flag;

public:
  string name;
  vector<string> tokens;
  Param() {
    count = 0;
    flag = false;
  }
  Param(const string& name) {
    setName(name);
    count = 0;
    flag = false;
  }
  void init(const Param& param) {
    setName(param.getName());
    count = 0;
    flag = false;
    tokens = param.tokens;
  }
  void setName(const string& name) {
    this->name = name;
  }
  string getName() const {
    return(name);
  }
  void setFlag() { flag = true; }
  void clearFlag() { flag = false; }
  bool checkFlag() const { return(flag); }
  void dump(ostream& ofile = cout);
  string str() const;
  void writeToFile(FILE * fp);
  int size() const { return( tokens.size() ); }
  string getString(const int i = 0) const;
  float getFloat(const int i = 0) const;
  bool isConstDouble(const int i = 0) const;
  double getDouble(const int i = 0) const;
  int getInt(const int i = 0) const;
  int8 getInt8(const int i = 0) const;
  bool getBool(const int i = 0) const;
  // recording parameter accesses...
  bool touch() { ++count; return(true); } // return true here for use in FOR_PARAM_MATCHING(S)
  int getTouchCount() const { return(count); }
  // calling process can pop the last string off...
  void pop() {
    if (tokens.size() > 0)
      tokens.resize(tokens.size()-1);
  }
  void popFront() {
    const string new_name = getString(0);
    if (tokens.size() > 0) tokens.erase(tokens.begin());
    setName(new_name);
  }
  int applyDefines(const map<const string,string>& defineMap);
};

// a useful MACRO for looping on all the parameters with the

// same name "S", where S is some string...
//#define FOR_PARAM(S) for (setCurrentParam(S);checkCurrentParam();setCurrentParamNext(S))

// March 2019: discarding the idea of a global currentParam, because it can get
// messed up when used recursively, in favor of something like this...
//typedef list<Param>::iterator ParamIter;

#define FOR_PARAM_MATCHING(S) for (list<Param>::iterator param = paramList_begin(); param != paramList_end(); ++param) \
    if ( (param->name == (S)) && param->touch() )

/*
 * #define FOR_PARAM_MATCHING_IN_LIST(PARAM_LIST,S) for (ParamIter param = PARAM_LIST.begin(); param != PARAM_LIST.end(); ++param) \
 *    if ( (param->name == (S)) && param->touch() )
 */

#define FOR_ALL_PARAM for (list<Param>::iterator param = paramList_begin(); param != paramList_end(); ++param)

// March 2019: parameter list no longer exposed to solvers. Use Params methods below
/*
  namespace ParamsData {
  extern list<Param> paramList_;
  }
*/

namespace Params {

  // param list accessor methods...
  list<Param>::iterator paramList_begin();
  list<Param>::iterator paramList_end();

  // these namespace members are declared "extern" so we can include
  // the namespace definition in a header file included by
  // multiple routines

  // =========================================
  // methods for setting up parameters
  // =========================================

  /**
   * initialize parameters from args (first) then from a file. In general,
   * this is all you need to call.
   */
  void initParams(int argc,char * argv[],const char * input_file);

  /**
   * add the parameters from the args using --NAME to
   * indicate the start of a new parameter line.
   */
  //void addParamsFromArgs(int argc,char * argv[],const bool verbose);

  /**
   * add the parameters from an input file
   */
  //void addParamsFromFile(list<Param>& paramList,const char * input_file,const bool verbose);

  /**
   * add a parameter from a string. This takes the boolean argument "continued" to
   * either add the tokens to the current Param (if true), or start a new parameter (if false).
   */
  //bool addParamFromString(list<Param>& paramList,const string& line,const bool continued = false);

  /**
   * adds a new param to the param list from the param ADD_PARAM...
   */
  Param * addParamFromAddParam(Param * add_param);

  /**
   * add a param object to the paramList vector
   */
  void addParam(Param * add_param);

  /**
   * remove any param(s) matching tokens after RM_PARAM...
   */
  bool rmParam(Param * rm_param);

  /**
   * turn a line into a parameter...
   */
  
  bool createParamFromString(Param& param,const string& line);

  // report the parameters. Note: in parallel, you should only call this
  // from rank 0...
  void dumpParams(ostream& ofile = cout);

  // report parameter usage based on access counts...
  void dumpParamUsage();

  // parameter substitution...
  //void replaceParamsSubString(const string& oldSubString,const string& newSubString);

  // =========================================
  // methods for accessing parameters
  // =========================================
  Param * getParam(const string &name);
  Param * getParamNoTouch(const string &name);
  bool checkParam(const string& name); // returns true if the requested param exists, and sets it to the current param
  void clearAllParamFlags();
  void deleteFlaggedParams();
  string getStringParam(const string& name);
  string getStringParam(const string& name,const string& default_value);
  float getFloatParam(const string& name);
  float getFloatParam(const string& name,const float default_value);
  double getDoubleParam(const string& name);
  double getDoubleParam(const string& name,const double default_value);
  int getIntParam(const string& name);
  int getIntParam(const string& name,const int default_value);
  int8 getInt8Param(const string& name);
  int8 getInt8Param(const string& name,const int8 default_value);
  bool getBoolParam(const string& name);
  bool getBoolParam(const string& name,const bool default_value);

}

#endif
