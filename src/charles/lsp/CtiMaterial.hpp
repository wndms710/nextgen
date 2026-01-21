#ifndef CTI_MATERIAL_HPP
#define CTI_MATERIAL_HPP

#include "Common.hpp" 
#include "MpiStuff.hpp"
#include "Params.hpp" 
using namespace std; 
using namespace MpiStuff; 
using namespace Params; 

class CtiMaterial {

protected:

  string name;
  string info_prefix;

public:

  CtiMaterial() {}

  virtual ~CtiMaterial() {}

  virtual void setInfoPrefix(const string& s) = 0;
  virtual void info() = 0;
  virtual void info(const double T0,const double T1,const int n) = 0;
  string getName() {return name;}

  // common methods
  virtual double calcRho(const double T) = 0;
  virtual double calcCp(const double T) = 0;

};

#endif
