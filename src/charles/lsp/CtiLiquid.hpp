#ifndef CTI_LIQUID_HPP
#define CTI_LIQUID_HPP

#include "CtiMaterial.hpp"

class CtiLiquid : public CtiMaterial {

  // =============================================================================
  // CtiLiquid is an abstract class that describes the liquid properties necessary
  // for liquid spray modeling...
  // =============================================================================

public:

  // parameters
  double R_UNIVERSAL;
  double Pref, Tboil, MW;

  CtiLiquid(const string& name, double Pref);
  void setInfoPrefix(const string& s);
  void info();
  void info(const double T0,const double T1,const int n);
    
  // liquid properties
  virtual double calcCp(const double T) = 0;
  virtual double calcRho(const double T) = 0;
  virtual double calcHvap(const double T) = 0;
  virtual double calcSigma(const double T) = 0;

  // vapor properties
  virtual double calcPv(const double T) = 0;
  virtual double calcKv(const double T) = 0;
  virtual double calcDv(const double T) = 0;
  virtual double calcMuv(const double T) = 0; 
  virtual double calcCpv(const double T) = 0;
  virtual double calcRhov(const double T) = 0;

  virtual ~CtiLiquid() {} 
};


class YawsLiquid : public CtiLiquid {

  // =============================================================================
  // see: Carl L. Yaws, Chemical Properties Handbook, McGraw-Hill, NY 1999.
  // =============================================================================

private:

  double Ar,Br,nr;
  double AL,nL;
  double AC,BC,CC,DC;
  double AV,BV,CV,DV,EV;
  double Ak,Bk,Ck;
  double AP,BP,CP,DP,EP;
  double Asigma,Tsigma;
  double Tcrit;
  
public:
  
  YawsLiquid(const string& name, double Pref);
  double calcRho(const double T);
  double calcHvap(const double T);
  double calcCp(const double T);
  double calcCpv(const double T);
  double calcKv(const double T);
  double calcPv(const double T);
  double calcSigma(const double T);
  double calcDv(const double T);
  double calcMuv(const double T);
  double calcRhov(const double T);
};

class MHB98Liquid : public CtiLiquid {

  // =============================================================================
  // see: Miller et al. 1998, Inter. J. of MultiPhase Flow 24 1024-1055
  // =============================================================================
  
private:

  double AL,BL,CL,DL,nL;
  double AV,BV,CV,DV;
  double Ak,Bk,Ck,Dk;
  double Am, Bm;
  double rho,Hvap,cp,cpv,kv;//,pv,sigma,dv,muv;

public:

  MHB98Liquid(const string& name, double Pref);
  
  double calcRho(const double T);
  double calcHvap(const double T);
  double calcCp(const double T);
  double calcCpv(const double T);
  double calcKv(const double T);
  double calcDv(const double T);
  double calcPv(const double T);
  double calcSigma(const double T);
  double calcMuv(const double T);
  double calcRhov(const double T);
};


class MHB98decane : public CtiLiquid {

  // =============================================================================
  // see: Miller et al. 1998, Inter. J. of MultiPhase Flow 24 1024-1055
  // =============================================================================

public:

  MHB98decane(const string& name, double Pref); 
  double calcRho(const double T);
  double calcHvap(const double T);
  double calcCp(const double T);
  double calcCpv(const double T);
  double calcKv(const double T);
  double calcDv(const double T);
  double calcPv(const double T);
  double calcSigma(const double T);
  double calcMuv(const double T);
  double calcRhov(const double T);

};

class VLSATLiquid : public CtiLiquid {

private:

  double Pcrit,Tcrit,rho_c; 
  double PV0,PV1;
  double DV0,DV1;
  double DL0,DL1,DL2; 
  double HV0,HV1,HV2; 
  double ST0,ST1,ST2;
  double KV0,KV1,KV2;
  double VV0,VV1,VV2,VV3; 
  double CV0,CV1,CV2,CV3,CV4; 
  double CL0,CL1,CL2,CL3,CL4; 

public:

  VLSATLiquid(const string& name, double Pref);

  double calcPv(const double T);
  double calcHvap(const double T);
  double calcSigma(const double T);
  double calcRho(const double T);
  double calcRhov(const double T);
  double calcCp(const double T);
  double calcCpv(const double T);
  double calcKv(const double T);
  double calcMuv(const double T);
  double calcDv(const double T);

};

class MHB98heptane : public CtiLiquid {

  // =============================================================================
  // see: Miller et al. 1998, Inter. J. of MultiPhase Flow 24 1024-1055
  // =============================================================================

private:

  double AL,BL,CL,DL,nL;
  double AV,BV,CV,DV;
  double Ak,Bk,Ck,Dk;
  double Am, Bm, Cm, Dm;
  double Ad, Bd, Cd, Dd;
  double rho,Hvap,cp,cpv,kv,Dv;//,pv,sigma,dv,muv;

public:

  MHB98heptane(const string& name, double Pref); 
  double calcRho(const double T);
  double calcHvap(const double T);
  double calcCp(const double T);
  double calcCpv(const double T);
  double calcKv(const double T);
  double calcDv(const double T);
  double calcPv(const double T);
  double calcSigma(const double T);
  double calcMuv(const double T);
  double calcRhov(const double T);

};

// example of how to add additional liquids using the 
// abstract class CtiLiquid.
// 1. add a new class that implements all the virtual 
//    functions in CtiLiquid: e.g. calcRho(T), etc...
// 2. modify the routine "newCtiLiquid()" (defined in CtiLiquid.cpp)
//    to return your liquid when requested.

class UserLiquid : public CtiLiquid {

private:

  double rho,Hvap,cp,cpv,kv,pv,muv,sigma;
  
public:
  
  UserLiquid(const string& name, double Pref);
  double calcRho(const double T);
  double calcHvap(const double T);
  double calcCp(const double T);
  double calcCpv(const double T);
  double calcKv(const double T);
  double calcPv(const double T);
  double calcSigma(const double T);
  double calcDv(const double T);
  double calcMuv(const double T);
  double calcRhov(const double T);
  
};

//CtiLiquid * newCtiLiquid(const string& name, double Pref);

// Note: if you modify/add a new liquid class, then this routine
// needs to return it when requested by name.

CtiLiquid * newCtiLiquid(Param * param,int &iarg);

#endif

