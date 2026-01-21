#ifndef CTI_SOLID_HPP
#define CTI_SOLID_HPP

#include "CtiMaterial.hpp"

class CtiSolid : public CtiMaterial {

  // =============================================================================
  // CtiSolid is an abstract class that describes the solid properties necessary
  // for solid particle tracking modeling...
  // =============================================================================

public:

 // parameters
  double rho, cp, k, Tmelt, Hmelt;
  double yieldStress;
  double elasticity;

  CtiSolid(const string& name);
  void setInfoPrefix(const string& s);
  void info();
  void info(const double T0,const double T1,const int n);
    
  // solid properties
  virtual double calcCp(const double T) = 0;
  virtual double calcK(const double T) = 0;
  virtual double calcRho(const double T) = 0;
  virtual double calcHmelt(const double T) = 0;

  double calcYieldStress() {return yieldStress;}
  double calcElasticity() {return elasticity;}

  virtual ~CtiSolid() {} 
};

class Dust : public CtiSolid {
  // ====================================
  // a simple solid dust based on Quartz
  // ====================================

  public:

  Dust(const string& name);
  double calcCp(const double T);
  double calcK(const double T);
  double calcRho(const double T);
  double calcHmelt(const double T);

  virtual ~Dust() {}
};

class UserSolid : public CtiSolid {
  // ====================================
  // a simple solid dust based on Quartz
  // ====================================

  public:

  UserSolid(const string& name);
  double calcCp(const double T);
  double calcK(const double T);
  double calcRho(const double T);
  double calcHmelt(const double T);

  virtual ~UserSolid() {}
};

CtiSolid * newCtiSolid(Param * param,int &iarg);

#endif

