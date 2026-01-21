#include "Params.hpp"
#include "MpiStuff.hpp"
#include "Macros.hpp"
#include "CtiRegister.hpp"

#include <fstream>
#include <assert.h>
#include <iomanip>
#include <string.h>
#include <stack>

void Param::dump(ostream& ofile) {
  ofile << name;
  for (unsigned int i = 0; i < tokens.size(); ++i)
    ofile << " " << tokens[i];
  ofile << endl;
}

string Param::str() const {
  std::stringstream ss;
  ss << name;
  for (unsigned int i = 0; i < tokens.size(); ++i)
    ss << " " << tokens[i];
  return(ss.str());
}

void Param::writeToFile(FILE *fp) {
  fprintf(fp, "%s =",name.c_str());
  for (unsigned int i = 0; i < tokens.size(); ++i)
    fprintf(fp, " %s",tokens[i].c_str());
  fprintf(fp, "\n");
}

string Param::getString(const int i) const {
  using namespace MpiStuff;
  if (int(tokens.size()) <= i) {
    if (mpi_rank == 0)
      cerr << "\n\n\n********************************************************************\n" <<
	"Error: no string token at position " << i << " of param " << name <<
	"\n********************************************************************\n\n\n" << endl;
    throw(0);
  }
  return(tokens[i]);
}

// =======================================================
// a note on error handling: these routines throw 0
// when a parameter is missing, and 1 when it is present
// but does not evaluate properly. These are both positive
// int errors consistent with the throwing of positive from
// collective operations. It is assumed that parameter parsing
// is (exclusively?) collective.
// =======================================================

float Param::getFloat(const int i) const {
  using namespace MpiStuff;
  if (int(tokens.size()) <= i) {
    if (mpi_rank == 0)
      cerr << "\n\n\n********************************************************************\n" <<
	"Error: no float token at position " << i << " of param " << name << "." <<
	"\n********************************************************************\n\n\n" << endl;
    throw(0);
  }
  float value;
  if (!from_string<float>(value,tokens[i],std::dec)) {
    CtiRegister::CtiData * var = CtiRegister::getCtiData(tokens[i]);
    if ((var == NULL)||(var->getType() != D_DATA)) {
      if (mpi_rank == 0)
	cerr << "\n\n\n********************************************************************\n" <<
	  "Error: value at pos: "<< i << " of param " << name << " is not float: " << tokens[i] <<
	  "\n********************************************************************\n\n\n" << endl;
      throw(1);
    }
    value = float(var->d());
  }
  return (value);
}

bool Param::isConstDouble(const int i) const {
  using namespace MpiStuff;
  if (int(tokens.size()) <= i) {
    if (mpi_rank == 0)
      cerr << "\n\n\n********************************************************************\n" <<
	"Error: no double token at position " << i << " of param " << name << "." <<
	"\n********************************************************************\n\n\n" << endl;
    throw(0);
  }
  double value;
  return from_string<double>(value,tokens[i],std::dec);
}

double Param::getDouble(const int i) const {
  using namespace MpiStuff;
  if (int(tokens.size()) <= i) {
    if (mpi_rank == 0)
      cerr << "\n\n\n********************************************************************\n" <<
	"Error: no double token at position " << i << " of param " << name << "." <<
	"\n********************************************************************\n\n\n" << endl;
    throw(0);
  }
  double value;
  if (!from_string<double>(value,tokens[i],std::dec)) {
    CtiRegister::CtiData * var = CtiRegister::getCtiData(tokens[i]);
    if ((var == NULL)||(var->getType() != D_DATA)) {
      if (mpi_rank == 0)
	cerr << "\n\n\n********************************************************************\n" <<
	  "Error: value at pos: "<< i << " of param " << name << " is not double: " << tokens[i] <<
	  "\n********************************************************************\n\n\n" << endl;
      throw(1);
    }
    value = var->d();
  }
  return (value);
}

// TODO: CI should we revisit parsing of int's and int8's now that we have int-registration
// in CtiRegister?...

int Param::getInt(const int i) const {
  using namespace MpiStuff;
  if (int(tokens.size()) <= i) {
    if (mpi_rank == 0)
      cerr << "\n\n\n********************************************************************\n" <<
	"Error: no int token at position " << i << " of param " << name << "." <<
	"\n********************************************************************\n\n\n" << endl;
    throw(0);
  }
  int value;
  if (!from_string<int>(value,tokens[i],std::dec)) {
    CtiRegister::CtiData * var = CtiRegister::getCtiData(tokens[i]);
    // note that the number_token evaluation if it succeeds, will return a value that is a double type
    if ( (var == NULL) || ( var->getType() != D_DATA)) {
      if (mpi_rank == 0)
        cerr << "\n\n\n********************************************************************\n" <<
          "Error: value at pos: "<< i << " of param " << name << " is not int: " << tokens[i] <<
          "\n********************************************************************\n\n\n" << endl;
      throw(1);
    }
    value = int(var->d());
  }
  return (value);
}

int8 Param::getInt8(const int i) const {
  using namespace MpiStuff;
  if (int(tokens.size()) <= i) {
    if (mpi_rank == 0)
      cerr << "\n\n\n********************************************************************\n" <<
	"Error: no int8 token at position " << i << " of param " << name << "." <<
	"\n********************************************************************\n\n\n" << endl;
    throw(0);
  }
  int8 value;
  if (!from_string<int8>(value,tokens[i],std::dec)) {
    CtiRegister::CtiData * var = CtiRegister::getCtiData(tokens[i]);
    // note that the number_token evaluation if it succeeds, will return a value that is a double type
    if ( (var == NULL) || ( var->getType() != D_DATA)) {
      if (mpi_rank == 0)
        cerr << "\n\n\n********************************************************************\n" <<
          "Error: value at pos: "<< i << " of param " << name << " is not int8: " << tokens[i] <<
          "\n********************************************************************\n\n\n" << endl;
      throw(1);
    }
    value = int8(var->d());
  }
  return (value);
}

bool Param::getBool(const int i) const {
  using namespace MpiStuff;
  if (int(tokens.size()) <= i) {
    // Make the default behavior of the presence of a parameter with no "true" or 1 equal to true... 
    if ((i == 0)&&(tokens.empty()))
      return true;
    if (mpi_rank == 0)
      cerr << "\n\n\n********************************************************************\n" <<
	"Error: no bool token at position " << i << " of param " << name << "." <<
	"\n********************************************************************\n\n\n" << endl;
    throw(0);
  }
  else if ((tokens[i] == "1")||(tokens[i] == "true")||(tokens[i] == "TRUE"))
    return(true);
  else if ((tokens[i] == "0")||(tokens[i] == "false")||(tokens[i] == "FALSE"))
    return(false);
  else {
    if (mpi_rank == 0)
      cerr << "\n\n\n********************************************************************\n" <<
	"Error: value at pos: "<< i << " of param " << name << " is not bool: " << tokens[i] <<
	"\n********************************************************************\n\n\n" << endl;
    throw(1);
  }
}

int Param::applyDefines(const map<const string,string>& defineMap) {

  int ierr = 0;

  // apply defines to the param name as well...
  std::size_t pos0 = name.find("$(");
  while (pos0 != string::npos) {
    std::size_t pos1 = name.find_first_of(")",pos0+2);
    if (pos1 == string::npos) {
      if (mpi_rank == 0) cerr << "Error: Param::applyDefines: closing bracket not found for DEFINE in name \"" << name << "\"" << endl;
      ierr = -1; pos0 = string::npos;
    }
    else {
      // grab the substring...
      const string var = name.substr(pos0+2,pos1-pos0-2);
      map<const string,string>::const_iterator iter = defineMap.find(var);
      if (iter == defineMap.end()) {
	if (mpi_rank == 0) cerr << "Error: Param::applyDefines: DEFINE not found for " << var << endl;
	ierr = -1; pos0 = string::npos;
      }
      else {
	// replace...
	name.replace(pos0,pos1-pos0+1,iter->second);
	// and look again, from the start...
	pos0 = name.find("$(");
      }
    }
  }

  for (int i = 0, tkn_size=tokens.size(); i < tkn_size; ++i) {
    int count = 0;
    std::size_t pos0 = tokens[i].find("$(");
    while (pos0 != string::npos) {
      std::size_t pos1 = tokens[i].find_first_of(")",pos0+2);
      if (pos1 == string::npos) {
	if (mpi_rank == 0) cerr << "Error: Param::applyDefines: closing bracket not found for DEFINE in token \"" << tokens[i] << "\"" << endl;
	ierr = -1; pos0 = string::npos;
      }
      else {
	// grab the substring...
	const string var = tokens[i].substr(pos0+2,pos1-pos0-2);
	map<const string,string>::const_iterator iter = defineMap.find(var);
	if (iter == defineMap.end()) {
	  if (mpi_rank == 0) cerr << "Error: Param::applyDefines: DEFINE not found for " << var << endl;
	  ierr = -1; pos0 = string::npos;
	}
	else {
	  // at this point, check how many times we have done substitution on this token. It is
	  // possible to get into recursive loops that can bring the solver down. e.g.
	  // DEFINE a $(b)
	  // DEFINE b $(a)*sqrt(2)
	  //
	  // F. Ham:
	  // March 2019: I think this recursive example is not valid anymore because the DEFINES are now
          // processed in the order they are received, so the first define makes no sense, because $(b)
          // does not exist yet.
	  ++count;
	  if (count > 1000) {
	    if (mpi_rank == 0) cerr << "Error: Param::applyDefines: DEFINE " << var << " " << iter->second << " is involved in infinite recursion." << endl;
	    ierr = -1; pos0 = string::npos;
	  }
	  else {
	    // replace...
	    tokens[i].replace(pos0,pos1-pos0+1,iter->second);
	    // and look again, from the start...
	    pos0 = tokens[i].find("$(");
	  }
	}
      }
    }
  }
  return ierr;

}

namespace Params {

  list<Param> activeParamList;
  map<const string,string> defineMap;
  set<string> unknownParamSet;

  list<Param>::iterator paramList_begin() {
    return activeParamList.begin();
  }

  list<Param>::iterator paramList_end() {
    return activeParamList.end();
  }

  bool addParamFromString(list<Param>& paramList,const string& line,const bool continued = false) {

    // look for a hash - this could be an entire line, or a comment at the end of a line...

    string::size_type hashPos = line.find_first_of("#",0);
    string trimmed_line = line.substr(0,hashPos);

    string::size_type contPos = trimmed_line.find_first_of("\\",0);
    if (contPos != string::npos)
      trimmed_line = line.substr(0,contPos);

    // start to parse...

    string delimiters = " \t=\r";
    string::size_type lastPos = trimmed_line.find_first_not_of(delimiters,0); // skip delimiters at beginning

    // check for an empty line...
    if (lastPos == string::npos)
      return(false); // breaks any continuation...

    // if a parameter is being continued, then we add this string's tokens to the
    // current parameter. Otherwise, we start a new parameter...

    string::size_type pos = trimmed_line.find_first_of(delimiters,lastPos); // find first "delimiter"

    Param * param = NULL;
    if (continued) {
      assert(!paramList.empty());
      // we are adding to the last param...
      param = &paramList.back();
      param->tokens.push_back(trimmed_line.substr(lastPos,pos-lastPos));
    }
    else {
      paramList.push_back(Param(trimmed_line.substr(lastPos,pos-lastPos)));
      param = &paramList.back();
    }

    // and continue parsing...
    lastPos = trimmed_line.find_first_not_of(delimiters,pos);
    pos = trimmed_line.find_first_of(delimiters,lastPos);
    while ((pos != string::npos)||(lastPos != string::npos)) {
      param->tokens.push_back(trimmed_line.substr(lastPos,pos-lastPos));
      lastPos = trimmed_line.find_first_not_of(delimiters,pos); // skip delimiters
      pos = trimmed_line.find_first_of(delimiters,lastPos);
    }

    // return true if we found a "\" before at the end...
    return(contPos != string::npos);

  }

  void addParamsFromCbuf(list<Param>& paramList,const char * const cbuf,const int count) {

    // the cbuf is a convenient way to transfer params with am MPI_Bcast. It is used
    // to allow only rank 0 to read an input or killfile.
    //
    // Each param token is separated by a single '\0' character, and each
    // separate parameter is separated by a double '\0\0'...

    int ic = 0;
    int ic_prev = 0;
    bool new_param = true;
    while (ic < count) {
      if (cbuf[ic] == '\0') {
        //cout << "got space: ic: " << ic << " ic_prev: " << ic_prev << endl;
        // a space means the current token is done...
        assert(ic > ic_prev);
        // there is a token...
        if (new_param) {
          paramList.push_back(Param(cbuf+ic_prev));
          //cout << "just added new param: " << paramList.back().name << endl;
          new_param = false;
        }
        else {
          // there are some tokens...
          paramList.back().tokens.push_back(cbuf+ic_prev);
          //cout << "tokens: \"" << paramList.back().tokens.back() << "\"" << endl;
        }
        ++ic;
        if (cbuf[ic] == '\0') {
          ++ic;
          new_param = true;
        }
        ic_prev = ic;
      }
      else {
        //cout << "got char: " << cbuf[ic] << endl;
        ++ic;
      }
    }
    // we should have ended with new_param set and ic_prev == ic == count
    // i.e. we ended with double '\0\0'...
    assert(ic == count);
    assert(ic_prev == count);
    assert(new_param);

  }

  int addParamsFromFile(list<Param>& paramList,const char * input_file,const bool verbose) {

    // return value:
    // -1: file not found
    // 0: empty file found
    // >0: new params added

    int count = 0;
    char * cbuf = NULL;
    if (mpi_rank == 0) {

      // only rank0 reads...

      const string delimiters = " \t=";
      const bool paramListWasEmpty = paramList.empty();
      list<Param>::iterator begin = paramList.end();
      if (!paramListWasEmpty) {
        // if the current paramList is not empty, store a pointer to the
        // last entry to facilitate exchange of new params with other ranks...
        --begin;
      }

      ifstream ifile;
      ifile.open(input_file);
      if (!ifile.is_open()) {
        if (verbose) cout <<
                       "# *********************************************************************\n" <<
                       "# Warning: input file \"" << input_file << "\" is not present or unreadable. Skipping\n" <<
                       "# *********************************************************************" << endl;
        count = -1;
      }
      else {
        // here we make use of the fact that the user can use "\" to
        // continue a parameter line in the solver...
        bool continued = false;
        string line;
        while (ifile.good()) {
          getline(ifile,line);
          if (verbose) cout << line << endl;
          continued = addParamFromString(paramList,line,continued);
        }
        ifile.close();

        // loop through the new params and count the char buf size req'd
        // to send them...

        if (paramListWasEmpty) {
          begin = paramList.begin();
        }
        else {
          ++begin;
        }

        count = 0;
        for (list<Param>::iterator iter = begin; iter != paramList.end(); ++iter) {
          count += iter->name.size() + 2; // one for the space after name and the final space...
          for (int it = 0, nit=iter->tokens.size(); it < nit; ++it) {
            count += iter->tokens[it].size() + 1;
          }
        }

        //cout << "got count: " << count << endl;

        if ((mpi_size > 1)&&(count > 0)) {
          assert(cbuf == NULL);
          cbuf = new char[count];
          int ic = 0;
          for (list<Param>::iterator iter = begin; iter != paramList.end(); ++iter) {
            memcpy(cbuf+ic,iter->name.c_str(),iter->name.size());
            ic += iter->name.size();
            cbuf[ic++] = '\0';
            for (int it = 0,nit=iter->tokens.size(); it < nit; ++it) {
              memcpy(cbuf+ic,iter->tokens[it].c_str(),iter->tokens[it].size());
              ic += iter->tokens[it].size();
              cbuf[ic++] = '\0';
            }
            cbuf[ic++] = '\0';
          }
          assert(ic == count);

          /*
          // check that new params can be correctly constructed from the cbuf...
          list<Param> newParamList;
          addParamsFromCbuf(newParamList,cbuf,count);
          cout << "about to check" << endl;
          list<Param>::iterator iter = begin;
          list<Param>::iterator new_iter = newParamList.begin();
          while (iter != paramList.end()) {
          assert(new_iter != newParamList.end());
          assert(iter->name == new_iter->name);
          assert(iter->tokens.size() == new_iter->tokens.size());
          for (int it = 0; it < iter->tokens.size(); ++it) {
          assert(iter->tokens[it] == new_iter->tokens[it]);
          }
          ++iter;
          ++new_iter;
          }
          assert(iter == paramList.end());
          assert(new_iter == newParamList.end());
          cout << "LOOKS GOOD" << endl;
          */

        }

      }

    }

    // now update paramList on other ranks...

    if (mpi_size > 1) {

      MPI_Bcast(&count,1,MPI_INT,0,mpi_comm);

      if (count > 0) {

        if (mpi_rank != 0) {
          assert(cbuf == NULL);
          cbuf = new char[count];
        }

        // bcast form rank 0...
        MPI_Bcast(cbuf,count,MPI_CHAR,0,mpi_comm);

        // every EXCEPT rank 0 add to their paramList...
        if (mpi_rank != 0)
          addParamsFromCbuf(paramList,cbuf,count);

        // everyone cleanup...
        delete[] cbuf;
        cbuf = NULL;

      }

    }

    // check...
    assert(cbuf == NULL);
    return count;

  }

  /*
    void dumpDefineMap() {
    cout << "defineMap: ";
    for (map<const string,string>::iterator iter = ParamsData::defineMap.begin(); iter != ParamsData::defineMap.end(); ++iter)
    cout << iter->first << "=" << iter->second << ", ";
    cout << endl;
    }
  */

  class ForLoopData {
  public:
    string var;
    // for now, assume index advances towards "last" by "inc",
    // but other looping strategies would be possible.
    int index,last,inc;
  };

  void expandLoopsAndDefines(list<Param>& newParamList,const list<Param>& paramList) {

    // we use CtiRegister's evaluation infrastructure to provide some guidance on what the
    // expressions can evaluate to, but many evaluations will fail, so turn off
    // verbosity. We turn it on again below...
    CtiRegister::setEvalVerbosity(false);

    stack<pair<list<Param>::const_iterator,ForLoopData> > forStack;

    Param current_param;
    list<Param>::const_iterator iter__ = paramList.begin();
    while (iter__ != paramList.end()) {
      // take a copy and apply current defines to this param, no matter what it is. This
      // allows the use of defined variables in any parameter, including FOR, etc...
      current_param.init(*iter__);
      if (current_param.applyDefines(defineMap) != 0) {
	CERR("defines failing for param: " << iter__->str());
      }
      // check if this param has any tokens that we can try to eval as
      // doubles, and if we can, report the full param again with these
      // evaluations...
      bool got_eval_double = false;
      std::stringstream ss;
      // only needed on rank 0...
      if (mpi_rank == 0) {
        ss << "# eval check: " << current_param.name;
        for (int ii = 0; ii < current_param.size(); ++ii) {
          //cout << ">>>>>>> working on param: " << current_param.tokens[ii];
          double value;
          if (!from_string<double>(value,current_param.tokens[ii],std::dec)) {
            CtiRegister::CtiData * var = CtiRegister::getCtiData(current_param.tokens[ii]);
            if ((var == NULL)||(var->getType() != D_DATA)) {
              ss << " " << current_param.tokens[ii];
              //cout << " string" << endl;
            }
            else {
              got_eval_double = true;
              ss << " " << var->d();
              //cout << " eval double" << var->d() << endl;
            }
          }
          else {
            ss << " " << value;
            //cout << " regular double " << value << endl;
          }
        }
        //CtiRegister::dumpCurrentData();
        //cout << "HERE IS ss: " << ss.str() << endl;
      }
      // ======================================================================
      // some params have special behavior (FOR/ENDFOR, DEFINE/UNDEF) and will
      // not necessarily appear in the final param list...
      // ======================================================================
      string token = MiscUtils::toUpperCase(current_param.getName());
      if (token == "DEFINE") {
	// ==========================================================
        // DEFINE A 1.234
        // DEFINE IMAGE_STUFF INTERVAL=10
	// ==========================================================
        if (current_param.size() != 2) {
	  CERR("DEFINE syntax problem: " << iter__->str());
	}
	else if (defineMap.find(current_param.tokens[0]) != defineMap.end()) {
	  if (mpi_rank == 0) cout <<
                               "# *********************************************************************\n" <<
                               "# Warning: " << current_param.tokens[0] << " already in use. Skipping " << iter__->str() <<
                               "\n# *********************************************************************" << endl;
	}
	else {
          assert(current_param.size() == 2);
          // removed parenthesis to prevent string substitution issues
          defineMap[current_param.tokens[0]] = current_param.tokens[1];
          if (mpi_rank == 0) {
            cout << "# DEFINE " << current_param.tokens[0] << " " << current_param.tokens[1] << endl;
            if (got_eval_double) cout << ss.str() << endl;
          }
        }
      }
      else if (token == "UNDEF") {
        if (current_param.size() < 1) {
	  CERR("UNDEF syntax problem: " << iter__->str());
	}
        else {
          map<const string,string>::iterator id = defineMap.find(current_param.tokens[0]);
          if (id == defineMap.end()) {
            if (mpi_rank == 0) cout <<
                                 "# *********************************************************************\n" <<
                                 "# Warning: " << current_param.tokens[0] << " not currently defined. Skipping " << iter__->str() <<
                                 "\n# *********************************************************************" << endl;
          }
          else {
            if (mpi_rank == 0) cout << "# UNDEF " << current_param.tokens[0] << endl;
            defineMap.erase(id);
          }
        }
      }
      else if ((token == "$FOR")||(token == "$ENDFOR")) {
        CERR("old $for/$endfor syntax no longer supported. Use FOR/ENDFOR. Example:\n" <<
             "FOR I = <start> TO <end> [BY <inc>]\n" <<
             "# some params you want repeated, accessiogn the loop variable as $(I) (just like DEFINE)...\n" <<
             "ENDFOR\n");
      }
      else if (token == "$RETURN") {
        CERR("old $return no longer supported.");
      }
      else if (token == "FOR") {
	// ==========================================================
	// for loop syntax:
	// FOR I 1 TO 10
	// FOR J 7 TO 11 BY 2
	// ==========================================================
        if (forStack.empty() || (forStack.top().first != iter__)) {
          if (mpi_rank == 0) {
            cout << "# " << current_param.str() << endl;
            if (got_eval_double) cout << ss.str() << endl;
          }
          // this is the first time in this loop, so add its loop variable as a define...
          ForLoopData fld;
          fld.var = current_param.getString(0);
	  if (defineMap.find(fld.var) != defineMap.end()) {
	    CERR("for loop variable already in use: " << fld.var);
	  }
	  fld.index = current_param.getInt(1);
	  token = MiscUtils::toUpperCase(current_param.getString(2));
	  if (token != "TO") {
	    CERR("for loop syntax problem: missing TO: " << current_param.str());
	  }
	  fld.last = current_param.getInt(3);
	  fld.inc = 1;
	  if (current_param.size() > 4) {
	    token = MiscUtils::toUpperCase(current_param.getString(4));
	    if (token != "BY") {
	      CERR("for loop syntax problem: expecting BY: " << current_param.str());
	    }
	    fld.inc = current_param.getInt(5);
	  }
	  // check that first, last and inc are consistent...
	  if (fld.inc >= 1) {
	    if (fld.last < fld.index) {
	      CERR("for loop indexing problem: last < first with inc positive: " << current_param.str());
	    }
	    else if ((fld.last > fld.index)&&((fld.last-fld.index)%fld.inc != 0)) {
	      CERR("for loop indexing problem: (last-first)%inc != 0: " << current_param.str());
	    }
	  }
	  else if (fld.inc <= -1) {
	    if (fld.last > fld.index) {
	      CERR("for loop indexing problem: last > first with inc negative: " << current_param.str());
	    }
	    else if ((fld.last < fld.index)&&((fld.index-fld.last)%fld.inc != 0)) {
	      CERR("for loop indexing problem: (first-last)%inc != 0: " << current_param.str());
	    }
	  }
	  else {
	    CERR("for loop indexing problem: inc != 0: " << current_param.str());
	  }
	  forStack.push(pair<list<Param>::const_iterator,ForLoopData>(iter__,fld));
          // define map takes a string, so just use the original string instead of
          // using a ss to convert the int index (currently at the loop start)
          // back to a string...
	  defineMap[fld.var] = current_param.getString(1);
          if (mpi_rank == 0) cout << "# setting loop variable " << fld.var << "=" << fld.index << endl;
	}
      }
      else if (token == "ENDFOR") {
	// the "FOR" associated with this "ENDFOR" should be top of stack...
	if (forStack.empty()) {
	  CERR("Params::expandLoopsAndDefines: FOR/ENDFOR pairs not matched");
	}
	// we have already parsed the param, so we skip error-checking here...
	const string var = forStack.top().second.var;
	// this var should be in the defineMap...
	map<const string,string>::iterator id = defineMap.find(var);
	assert(id != defineMap.end());
        if (forStack.top().second.index == forStack.top().second.last) {
	  // we are done with this for loop...
	  forStack.pop();
	  // also remove our loop variable from the current defineMap...
	  defineMap.erase(id);
          if (mpi_rank == 0) cout << "# ENDFOR loop variable " << var << endl;
	}
	else {
	  // inc loop index...
	  forStack.top().second.index += forStack.top().second.inc;
	  // and modify define...
	  stringstream ss2;
	  ss2 << forStack.top().second.index;
	  id->second = ss2.str();
	  // reset loop iter to the FOR position...
	  iter__ = forStack.top().first;
          if (mpi_rank == 0) cout << "# setting loop variable " << var << "=" << id->second << endl;
	  continue;
	}
      }
      else {
        newParamList.push_back(current_param);
        if (mpi_rank == 0) {
          // push the current param to cout...
          cout << current_param.str() << endl;
          if (got_eval_double) cout << ss.str() << endl;
        }
      }
      // proceed to the next param...
      ++iter__;
    }

    // check that FOR/ENDFOR loops paired...
    if (!forStack.empty()) {
      CERR("Params::expandLoopsAndDefines: FOR/ENDFOR pairs not matched");
    }

    // and turn it back...
    CtiRegister::setEvalVerbosity(true);
    CtiRegister::clearCurrentData();

  }

  void initParams(int argc,char * argv[],const char * input_file) {

    // should only call this once...

    assert(activeParamList.empty());

    if (mpi_rank == 0)
      cout << "================= init params: output of raw sources ==================" << endl;

    // loop through the args and look for params specified with --KEYWORD VALUE
    // and one or more files specified with -i <input_file>...

    list<Param> paramList;
    bool custom_input_file = false;
    bool first_param_from_args = true;

    string line;
    assert(line.empty());

    int iargc = 1;
    while (iargc < argc) {

      if ((strlen(argv[iargc]) == 2)&&(argv[iargc][0] == '-')&&(argv[iargc][1] == 'i')) {

	// this is a command line input file. Flush any params in line...
	if (!line.empty()) {
	  if (first_param_from_args) {
	    if (mpi_rank == 0) cout << "# ------------------ param(s) from command line args ------------------" << endl;
	    first_param_from_args = false;
	  }
	  if (mpi_rank == 0) cout << line << endl;
	  addParamFromString(paramList,line);
	  line.clear();
	  assert(line.empty());
	}

	// next argument should be the input_file...
	++iargc;
	if (iargc == argc) {
	  CERR("command line formatting: expecting -i <input_file>");
	}

	custom_input_file = true; // this will skip default input_file reading at the end...
	first_param_from_args = true; // if there are more args after this -i <input_file>, we get the args comment again
        if (mpi_rank == 0) cout << "# ------------------ param(s) from input file: " << argv[iargc] << "..." << endl;
	addParamsFromFile(paramList,argv[iargc],mpi_rank==0); // TODO: params should be mpi-aware

      }
      else if ((strlen(argv[iargc]) >= 3)&&(argv[iargc][0] == '-')&&(argv[iargc][1] == '-')) {

	// if there is already something in line, add it...
	if (!line.empty()) {
	  if (first_param_from_args) {
	    if (mpi_rank == 0) cout << "# ------------------ param(s) from command line args ------------------" << endl;
	    first_param_from_args = false;
	  }
	  if (mpi_rank == 0) cout << line << endl;
	  addParamFromString(paramList,line);
	}

	// start the next...
	line = &(argv[iargc][2]);

      }
      else {

	if (line.empty()) {
	  CERR("command line formatting: use --KEYWORD to start a new param");
	}

	line.append(" ");
	line.append(argv[iargc]);

      }

      iargc++;

    }

    // complete the last command line param, if present...
    if (!line.empty()) {
      if (first_param_from_args) {
	if (mpi_rank == 0) cout << "# ------------------ param(s) from command line args ------------------" << endl;
	first_param_from_args = false;
      }
      if (mpi_rank == 0) cout << line << endl;
      addParamFromString(paramList,line);
    }

    // and read from params from the default input file, if no custom one
    // was specified...
    if (!custom_input_file) {
      if (mpi_rank == 0) cout << "# ------------------ param(s) from input file: " << input_file << "..." << endl;
      addParamsFromFile(paramList,input_file,mpi_rank==0); // TODO: params should be mpi-aware
    }

    if (mpi_rank == 0) {
      cout <<
        "=================== init params: end of raw sources ===================\n" <<
        "================ init params: expand loops and defines ================" << endl;
    }

    assert(activeParamList.empty());
    expandLoopsAndDefines(activeParamList,paramList);

    if (mpi_rank == 0)
      cout << "======================== end of init params ===========================" << endl;

  }

  Param * addParamFromAddParam(Param * add_param) {
    // here the user has prefixed the param tokens with the
    // ADD_PARAM keyword. For example:
    // ADD_PARAM WRITE_IMAGE INTERVAL 100 NAME x0 [...]
    //
    // shift the first token to become the keyword and add this param...
    // also: return pointer to param for use by calling function
    // if needed
    assert(add_param->name == "ADD_PARAM");
    if (add_param->size() > 0) {
      activeParamList.push_back(Param(add_param->tokens[0]));
      Param * param = &(activeParamList.back());
      for (int i = 1; i < add_param->size(); ++i)
        param->tokens.push_back(add_param->tokens[i]);
      return param;
    }
    return NULL;
  }

  void addParam(Param * add_param) {
    activeParamList.push_back(Param(add_param->name));
    Param * param = &(activeParamList.back());
    param->tokens.resize(add_param->size());
    for (int i = 0; i < add_param->size(); ++i)
      param->tokens[i] = add_param->tokens[i];
  }

  bool rmParam(Param * rm_param) {
    // remove any params that EXACTLY match the passed tokens. At any point,
    // the rm_param can contain a "*", at which point the params are considered
    // to have matched if they have matched to that point...
    assert(rm_param->name == "RM_PARAM");
    if (rm_param->size() > 0) {
      list<Param>::iterator param = activeParamList.begin();
      while (param != activeParamList.end()) {
        //cout << "WORKING on param: " << param->str() << endl;
        bool matched = false;
        if (param->name == rm_param->tokens[0]) {
          // matched the name...
          //cout << " > name matched!... " << endl;
          for (int i = 1; i < rm_param->size(); ++i) {
            if (rm_param->tokens[i] == "*") {
              //cout << " > got * - matched" << endl;
              matched = true;
              break;
            }
            else if (i <= param->size()) {
              if (rm_param->tokens[i] != param->tokens[i-1]) {
                //cout << " > mismatch: " << rm_param->tokens[i] << " != " << param->tokens[i-1] << endl;
                matched = false;
                break;
              }
              else if ((i == param->size())&&(i == rm_param->size()-1)) {
                //cout << " > all matched" << endl;
                matched = true;
                break;
              }
            }
            else {
              //cout << " > param out of tokens" << endl;
              matched = false;
              break;
            }
          }
        }
        if (matched) {
          //cout << " > matched: removing..." << endl;
          list<Param>::iterator param_copy = param;
          ++param;
          activeParamList.erase(param_copy);
          return true;
        }
        else {
          //cout << " > did not match: leaving..." << endl;
          ++param;
        }
      }
    }
    return false;
  }

  /*
    void addParam(Param& param) {
    using namespace ParamsData;
    activeParamList.push_back(param);
    }
  */

  bool createParamFromString(Param& param,const string& line) {

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
    if (lastPos == string::npos)
      return(false); // no parameter created

    // if a parameter is being continued, then we add this string's tokens to the
    // current parameter. Otherwise, we start a new parameter...

    string::size_type pos = trimmed_line.find_first_of(delimiters,lastPos); // find first "delimiter"
    param.setName(trimmed_line.substr(lastPos,pos-lastPos));

    // and parse tokens...
    param.tokens.clear();
    lastPos = trimmed_line.find_first_not_of(delimiters,pos);
    pos = trimmed_line.find_first_of(delimiters,lastPos);
    while (pos != string::npos || lastPos != string::npos) {
      param.tokens.push_back(trimmed_line.substr(lastPos,pos-lastPos));
      lastPos = trimmed_line.find_first_not_of(delimiters,pos); // skip delimiters
      pos = trimmed_line.find_first_of(delimiters,lastPos);
    }

    // parameter was created...
    return(true);

  }

  void dumpParams(ostream& ofile) {

    // this is a non-collective routine. Call just from rank 0

    assert( mpi_rank == 0 ); // generally should be called from rank 0 only...

    if (&ofile == &cout)
      ofile << "=========================== param summary =============================" << endl;

    for (list<Param>::iterator p = activeParamList.begin(); p != activeParamList.end(); ++p)
      p->dump(ofile);

    if (&ofile == &cout)
      ofile << "======================== end of param summary =========================" << endl;

  }

  void dumpParamUsage() {

    // this is a collective routine...

    int np = activeParamList.size();
    int * my_touch_count = new int[2*np];

    int i = 0;
    for (list<Param>::iterator p = activeParamList.begin(); p != activeParamList.end(); ++p) {
      my_touch_count[i++] = p->getTouchCount();
      my_touch_count[i++] = -p->getTouchCount();
    }

    int * touch_count = NULL;
    if (mpi_rank == 0)
      touch_count = new int[2*np];

    MPI_Reduce(my_touch_count,touch_count,2*np,MPI_INT,MPI_MIN,0,mpi_comm);

    delete[] my_touch_count;

    if (mpi_rank == 0) {

      cout << "============================ param usage ==============================\n" <<
	"The following params were never accessed:" << endl;

      i = 0;
      for (list<Param>::iterator p = activeParamList.begin(); p != activeParamList.end(); ++p) {
        int touch_min = touch_count[i++];
        int touch_max = -touch_count[i++];
        if ((touch_min == 0)&&(touch_max == 0))
          cout << p->name << endl;
      }

      cout << "-----------------------------------------------------------------------\n" <<
	"The following params were accessed differently on different ranks:" << endl;

      i = 0;
      for (list<Param>::iterator p = activeParamList.begin(); p != activeParamList.end(); ++p) {
        int touch_min = touch_count[i++];
        int touch_max = -touch_count[i++];
        if (touch_min != touch_max)
          cout << p->name << ", access count range: " << touch_min << " to " << touch_max << endl;
      }

      cout << "------------------------------------------------------------------------\n" <<
	"The following params were accessed uniformly across ranks:" << endl;

      i = 0;
      for (list<Param>::iterator p = activeParamList.begin(); p != activeParamList.end(); ++p) {
        int touch_min = touch_count[i++];
        int touch_max = -touch_count[i++];
        if ((touch_min > 0)&&(touch_min == touch_max))
          cout << p->name << ", access count: " << touch_min << endl;
      }

      delete[] touch_count;

      cout << "-----------------------------------------------------------------------\n" <<
	"The solver checked or requested the following params that were not specified: " << endl;

      for (set<string>::iterator up = unknownParamSet.begin(); up != unknownParamSet.end(); ++up)
        cout << *up << endl;

      cout << "======================== end of param usage ===========================" << endl;
    }

  }

  /*
    void replaceParamsSubString(const string& oldSubString,const string& newSubString) {

    using namespace ParamsData;

    for (list<Param>::iterator p = paramList.begin(); p != paramList.end(); ++p) {
    for (unsigned int i = 0; i < p->tokens.size(); ++i) {
    size_t np = 0;
    while ((np = p->tokens[i].find(oldSubString)) != string::npos) {
    string before = p->tokens[i].substr(0,np);
    string after  = p->tokens[i].substr(np+oldSubString.size(),string::npos);
    if (after.size() == 0) {
    if (before.size() == 0) {
    p->tokens[i] = newSubString;
    }
    else {
    p->tokens[i] = before+newSubString;
    }
    break;
    }
    else if (before.size() == 0) {
    p->tokens[i] = newSubString+after;
    }
    else {
    p->tokens[i] = before+newSubString+after;
    }
    np += newSubString.size() - oldSubString.size();
    }
    }
    }

    }
  */

  // =========================================
  // methods for accessing parameters
  // =========================================

  Param * getParam(const string &name) {

    for (list<Param>::iterator param = activeParamList.begin(); param != activeParamList.end(); ++param) {
      if (param->name == name) {
        param->touch();
	return &*param;
      }
    }

    unknownParamSet.insert(name);

    return NULL;

  }

  Param * getParamNoTouch(const string &name) {

    for (list<Param>::iterator param = activeParamList.begin(); param != activeParamList.end(); ++param) {
      if (param->name == name) {
	return &*param;
      }
    }

    return NULL;

  }

  bool checkParam(const string& name) {

    Param * param = getParam(name);

    return param != NULL;

  }

  void clearAllParamFlags() {

    for (list<Param>::iterator param = activeParamList.begin(); param != activeParamList.end(); ++param)
      param->clearFlag();

  }

  void deleteFlaggedParams() {

    // loop through params delete any that have their flag set...
    list<Param>::iterator param = activeParamList.begin();
    while (param != activeParamList.end()) {
      list<Param>::iterator param_copy = param;
      ++param;
      if (param_copy->checkFlag()) {
        activeParamList.erase(param_copy);
      }
    }

  }

  string getStringParam(const string& name) {

    if (Param * param = getParam(name)) {
      return param->getString();
    }
    else {
      CERR("cannot find string param: " << name);
    }

  }

  string getStringParam(const string& name,const string& default_value) {

    // here use the "no touch" version of getParam so we can
    // insert the default value into the unknownParamSet...
    if (Param * param = getParamNoTouch(name)) {
      param->touch();
      return param->getString();
    }
    else {
      std::stringstream ss;
      ss << name << ", default used: " << default_value;
      unknownParamSet.insert(ss.str());
      return default_value;
    }

  }

  float getFloatParam(const string& name) {

    if (Param * param = getParam(name)) {
      return param->getFloat();
    }
    else {
      CERR("cannot find float parameter " << name);
    }

  }

  float getFloatParam(const string& name,const float default_value) {

    if (Param * param = getParamNoTouch(name)) {
      param->touch();
      return param->getFloat();
    }
    else {
      std::stringstream ss;
      ss << name << ", default used: " << default_value;
      unknownParamSet.insert(ss.str());
      return default_value;
    }

  }

  double getDoubleParam(const string& name) {

    if (Param * param = getParam(name)) {
      return param->getDouble();
    }
    else {
      CERR("cannot find double parameter " << name);
    }

  }

  double getDoubleParam(const string& name,const double default_value) {

    if (Param * param = getParamNoTouch(name)) {
      param->touch();
      return param->getDouble();
    }
    else {
      std::stringstream ss;
      ss << name << ", default used: " << default_value;
      unknownParamSet.insert(ss.str());
      return default_value;
    }

  }

  int getIntParam(const string& name) {

    if (Param * param = getParam(name)) {
      return param->getInt();
    }
    else {
      CERR("cannot find int parameter " << name);
    }

  }

  int getIntParam(const string& name,const int default_value) {

    if (Param * param = getParamNoTouch(name)) {
      param->touch();
      return param->getInt();
    }
    else {
      std::stringstream ss;
      ss << name << ", default used: " << default_value;
      unknownParamSet.insert(ss.str());
      return default_value;
    }

  }

  int8 getInt8Param(const string& name) {

    if (Param * param = getParam(name)) {
      return param->getInt8();
    }
    else {
      CERR("cannot find int8 parameter " << name);
    }

  }

  int8 getInt8Param(const string& name,const int8 default_value) {

    if (Param * param = getParamNoTouch(name)) {
      param->touch();
      return param->getInt8();
    }
    else {
      std::stringstream ss;
      ss << name << ", default used: " << default_value;
      unknownParamSet.insert(ss.str());
      return default_value;
    }

  }

  bool getBoolParam(const string& name) {

    if (Param * param = getParam(name)) {
      return param->getBool();
    }
    else {
      CERR("cannot find int8 parameter " << name);
    }

  }

  bool getBoolParam(const string& name,const bool default_value) {

    if (Param * param = getParamNoTouch(name)) {
      param->touch();
      return param->getBool();
    }
    else {
      std::stringstream ss;
      ss << name << ", default used: " << default_value;
      unknownParamSet.insert(ss.str());
      return default_value;
    }

  }

}
