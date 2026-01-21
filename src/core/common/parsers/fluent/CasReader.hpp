#ifndef FLUENT_CAS_READER_HPP
#define FLUENT_CAS_READER_HPP

#include "Common.hpp"
#include "CTI.hpp"
using namespace CTI;
#include "parsers/fluent/FluentParser.hpp"

class CasReader : public FluentParser {
public:

  CasReader() : FluentParser() {};

  CasReader(const string& filename) : FluentParser(filename) {};

  ~CasReader() {};


  void readBinaryFaces(int * d, const int nNodes) {

    // node list + two cells
    int nvars = nNodes + 2;

    for (int i=0; i<nvars; ++i) {
      d[i] = getNextBinary<int,int>();
    }
  }

  void readBinaryInts(int * d, const int nvars) {

    for (int i=0; i<nvars; ++i) {
      d[i] = getNextBinary<int,int>();
    }
  }

  template<typename T>
  void readBinaryNodes(double * d, const int dim) {
    for (int i=0; i<dim; ++i) {
      d[i] = getNextBinary<double,T>();
    }
  }

};
  #endif
