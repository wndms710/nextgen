#include "CtiSolid.hpp"
#include "MiscUtils.hpp"

CtiSolid::CtiSolid(const string& name) : CtiMaterial() {
  this->name = name;
  info_prefix = "";
  
  yieldStress = getDoubleParam(name+".YIELD_STRESS",-1);
  elasticity = getDoubleParam(name+".ELASTICITY",-1);;
}

void CtiSolid::setInfoPrefix(const string& s) {
  info_prefix = s+" ";
}

void CtiSolid::info() {
  cout << "======================= solid info ========================" << endl;
  cout << "CtiSolid: " << name << endl;
  cout << "Tmelt [K]: " << Tmelt << endl;
  cout << "density at 300K [kg/m3]: " << calcRho(300) << endl;
  cout << "heat capacity (cp) at 300K [J/kg/K]: " << calcCp(300) << endl;
  cout << "heat conductivity (k) at 300K [J/m/K/s]: " << calcK(300) << endl;
  cout << "=================== end of solid info ======================" << endl;
}

void CtiSolid::info(const double T0,const double T1,const int n = 25) {
  const int nset = 12;
  cout << "# CtiSolid " << name << ", Tmelt=" << Tmelt << "[K]" << endl;
  cout << setw(nset) << left << 
          "# name" << setw(nset) << left <<
          "T" << setw(nset) << left <<
          "rho" << setw(nset) << left <<     
          "Cp" << setw(nset) << left <<
          "k" << setw(nset) << left <<       
          "Hmelt" << endl;
  cout << setw(nset) << left <<
          "# [none]" << setw(nset) << left <<
          "[K]" << setw(nset) << left <<
          "[kg/m3]" << setw(nset) << left <<
          "[J/(kg*K)]" << setw(nset) << left <<
          "[W/(m*K)]" << setw(nset) << left <<
          "[J/kg]" << endl;
  for (int i = 0; i < n; ++i) {
    const double T = T0 + double(i)/double(n-1)*(T1-T0);
    cout << setw(nset) << left <<
            info_prefix << setw(nset) << left <<
            T << setw(nset) << left << 
            calcRho(T) << setw(nset) << left <<
            calcCp(T) << setw(nset) << left <<
            calcK(T) << setw(nset) << left <<
            calcHmelt(T) << endl;
  }
}

Dust::Dust(const string& name) : CtiSolid(name) {
  IF_RANK0 cout << "Dust()" << endl;

  // http://www.mt-berlin.com/frames_cryst/descriptions/quartz%20.htm
  rho = 2200.;          // kg/m3
  cp = 1400.;            // J/kg/K
  k = 1.2;             // J/s/m/K
  Tmelt = 1715. + 273.; // K
  Hmelt = 0.;
}

double Dust::calcRho(const double T) {
  return(rho);
}

double Dust::calcCp(const double T) {
  return(cp);
}

double Dust::calcK(const double T) {
  return(k);
}

double Dust::calcHmelt(const double T) {
  return(Hmelt);
}

UserSolid::UserSolid(const string& name) : CtiSolid(name) {
  IF_RANK0 cout << "UserSolid()" << endl;

  rho     = getDoubleParam(name+".RHO");
  cp      = getDoubleParam(name+".CP");
  k       = getDoubleParam(name+".K");
  Tmelt   = getDoubleParam(name+".TMELT");
  Hmelt = 0.;
}

double UserSolid::calcRho(const double T) {
  return(rho);
}

double UserSolid::calcCp(const double T) {
  return(cp);
}

double UserSolid::calcK(const double T) {
  return(k);
}

double UserSolid::calcHmelt(const double T) {
  return(Hmelt);
}

CtiSolid * newCtiSolid(Param * param,int &iarg) {

  // first token should be a recognized name...
  
  if (iarg >= param->size()) {
    CERR("LP.MATERIAL SOLID missing id. Possible choices include:\n" << 
         "USER, DUST");
  }
  const string id = param->getString(iarg++);
  string name = id;
  
  // parse anything else...
  
  while (iarg < param->size()) {
    string token = MiscUtils::toUpperCase(param->getString(iarg++));
    if (token == "NAME") {
      name = param->getString(iarg++);
    }
    else {
      if (mpi_rank == 0) cout << " > skipping unrecognized LP.MATERIAL token: " << token << endl;
    }
  }
  
  // this function gets material props based on the passed id. It can  
  // get smarter as the number and type of supported liquids grows.
  if (id == "DUST") {
    return(new Dust(name));
  }
  else if (id == "USER") {
    return(new UserSolid(name));
  } 
  else {
    CERR("unrecognized LP.MATERIAL id: " << id << ". Possible choices include:\n" << 
         "USER, DUST");
  }
  
  // should never get here...
  return NULL;
  
}
