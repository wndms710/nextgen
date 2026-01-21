// =========
// Flamelet 
// =========

#ifndef FLAMELET_HPP
#define FLAMELET_HPP

#include "Params.hpp"
using namespace Params;
#include "ArrayNd.hpp"
#include "properties.hpp"
#include "PdfConvolution.hpp" 

#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"

#include <fstream>
#include <sstream>
#include <list>
#include <map>
//#include <tr1/unordered_map>
//using std::tr1::unordered_map;

// pdf types
enum PRESUMED_PDF_FUNCTIONS {
  BETA_PDF
};

// define a type for sorting purposes 
typedef std::pair<int,double> pairIntDbl; 

inline bool isSecondSmallerPair(const pairIntDbl& lhs, const pairIntDbl& rhs ){
  return(lhs.second < rhs.second);
}

// =================================================
//  Test if progress variable is larger or smaller
// =================================================
template<class T>
bool isProgSmallerFlamelet(const T* flamelet1, const T* flamelet2) {
  
  double Cmax1 = flamelet1->Cmax;
  double Cmax2 = flamelet2->Cmax;

  return(Cmax2>Cmax1);
}

template<class T>
bool isPhiSmallerFlamelet(const T* flamelet1, const T* flamelet2) {
  
  double phi1 = flamelet1->phi;
  double phi2 = flamelet2->phi;

  return(phi2>phi1);
    
}

// ============================
// class that handles flamelets
// ============================
class Flamelet {

public: 
  string filename;
  double pressure;
  double Cmin;
  double Cmax;
  double Tmax;

  // vector of data
  vector<string>   varNameVec;
  vector<double *> varDataVec;
 
  int    npoints;  

  // HACK for dual-stream 
  double z2;

  // make it protected so that a derived class 
  // can have access to these var
protected:
  
  int    nspecies;
  
  // vector of species
  vector<string>  speciesNameVec;
  
  // progress variable
  vector<string> progSpeciesNameVec;
  vector<double> progSpeciesWeightVec;

  // header
  vector<string> headerVec;
 
  // pdf function
  PRESUMED_PDF_FUNCTIONS pdf_function;
 
  Cantera::IdealGasMix *mix;
  Cantera::Transport *trans;

private:

  // pdf
  double *pdf;
  double pdfMean;
  double pdfVar;
  double *x_pdf;
#ifndef NEW_CHEM
  vector<pairIntDbl> x_pdfPair;
#endif
  int *idx_pdf; // for sorted pdf ..
  double xmin_pdf;
  double xmax_pdf;
  double invLx_pdf;

  map<string, double *> varHash;
  //unordered_map<string, double *> varHash;
    
public:
 
  Flamelet() {
    pdf = NULL;
    // initialize the mean and var to someting huge and negative
    pdfMean  = -1.0e+20;
    pdfVar   = -1.0e+20;
    
    xmin_pdf = 1.0e+20;
    xmax_pdf = -1.0e+20;
    x_pdf    = NULL;
  }

  virtual ~Flamelet() {

    for (int ivar=0; ivar<varDataVec.size(); ++ivar)
      delete[] varDataVec[ivar];
    
#ifdef NEW_CHEM
    DELETE(pdf); 
    DELETE(x_pdf); 
#else
    if ( pdf != NULL )
      delete[] pdf;
    x_pdf = NULL;
#endif

  }

  // =======================
  // initialize the flamelet
  // =======================
  void init(const string& filename) {
    // make sure that the flamelet has not alreay been initialized
    assert(varNameVec.size() == 0);
    assert(varDataVec.size() == 0);
    
    // ======================
    // read the flamelet file
    // ====================== 
    readFlameletFile(filename);
   
    // ==============
    // clean up names
    // ==============
    cleanupNamesAndInitSpeciesNames();   

    // ==============================
    // clean up negative source terms
    // ==============================
    cleanupNegativeMassFractions(); 
   
    // =====================================
    // progress variable species and weights
    // =====================================
    addProgAndSrcProg();

    // =====================================
    // divide source terms by rho (and cp if necessary)
    // =====================================
    divSrcByRho();
    divSrcByCp();
 
    // ================
    // add lambda / cp
    // ================
    addLocp();
    addGamma();

    // ==============
    // check for NAN
    // ==============
    if (checkParam("FLAMELET.CHECK_NAN"))
      checkForNaN();
    
    // ======================================
    // dump variables range
    // ======================================
    if (checkParam("FLAMELET.DUMP_VARS_RANGE")) 
      dumpVarsRange();
      
    // =========
    // pdf stuff
    // =========
    if (Param * param = getParam("FLAMELET.PRESUMED_PDF")) { 
      string pdf_name = param->getString(0);
      //cout << "    presumed pdf:                            " << pdf_name << endl;
      if (pdf_name == "BETA_FUNCTION") 
	      pdf_function = BETA_PDF;
      else 
	      CERR_S("unrecognized presumed pdf " << pdf_name);
    }
    else {
      //cout << "    presumed pdf:                            " << "BETA_FUNCTION"<< endl;
      pdf_function = BETA_PDF;
    }

    // if the PDF is not allocated, allocate the memory
    if ( pdf == NULL )
      pdf = new double[npoints];

    string mechanism = getStringParam("CHEMTABLE.MECHANISM");
    string transport = getStringParam("CHEMTABLE.MIXING");
    if (mechanism=="GRI3.0") mechanism = "gri30.xml";
    if (transport=="MIXTURE_AVERAGED")    transport = "Mix";
    else if (transport=="MULTICOMPONENT") transport = "Multi";
    mix   = new Cantera::IdealGasMix(mechanism);
    trans = newTransportMgr(transport, mix, 0);  // XXX check this
  }

  // ====================================
  // get the maximum value for a variable
  // ====================================
  double getVarMax(const string& varname){
    // variable
    double * var = getVarPtr(varname);
    assert(npoints>0);
    // get the max
    double varmax = var[0];
    for (int i=1; i<npoints; ++i)
      varmax = max(varmax,var[i]);
    
    return(varmax);
  }
  
  // ====================================
  // get the minimum value for a variable
  // ====================================
  double getVarMin(const string& varname){
    // variable
    double * var = getVarPtr(varname);
    assert(npoints>0);
    // get the min
    double varmin = var[0];
    for (int i=1; i<npoints; ++i)
      varmin = min(varmin,var[i]);
    
    return(varmin);
  }

  // =========================================
  // check if a variable is loaded
  // =========================================
  bool checkVar(const string& varname) {

    for (int i=0; i<varNameVec.size(); ++i)
      if ( varNameVec[i] == varname ) 
	return(true);

    return(false);
  }

  // =========================================
  // find a variable in the flamelet data list
  // and return a pointer to the data
  // =========================================
  double * getVarPtr(const string& varname) {
    int count = 0;

    for (int i=0; i<varNameVec.size(); ++i)
      if ( varNameVec[i] == varname ) 
	break;
      else
	++count;
    
    // if the variable is not found
    if ( count == varNameVec.size() ) 
      CERR_S("could not find variable " << varname);

    // return the data pointer
    return(varDataVec[count]);
  }
 
  // ====================================================
  // add variables to the flamelet and return the pointer
  // ====================================================
  double * addVarAndReturnDataPtr(const string& varname) {

    // first make sure it is not in the list
    for (int i=0; i<varNameVec.size(); ++i)
      if ( varNameVec[i] == varname )
	CERR_S("variable " << varname << " is already in the varNameVec");
    varNameVec.push_back(varname);
    double * tmpvec = new double[npoints];
    varDataVec.push_back(tmpvec);
    for (int i=0; i<npoints; ++i)
      tmpvec[i] = 0.0;
    return(tmpvec);
  }
  
  // =============================
  // add variables to the flamelet
  // =============================
  void addVar(const string& varname) { 
    double * tmpvec = addVarAndReturnDataPtr(varname);
  }

  // =======================================
  // dump the flamelet into a tecplot format
  // =======================================
  void writeTecplotfile(const string& filename) {
    string tecplotfilename = filename;
    tecplotfilename.append(".dat");
    cout << "  > writing flamelet data to the Tecplot file: "<< tecplotfilename << endl;
    FILE * fp = fopen(tecplotfilename.c_str(),"w");
    assert(fp);
    fprintf(fp,"title = \"flamelet\"\n");
    fprintf(fp,"variables = ");
    for (int var=0; var<varNameVec.size(); ++var)
      fprintf(fp,"\"%s\"\n",varNameVec[var].c_str());
    fprintf(fp,"zone i = %d, j = %d, k = %d, DATAPACKING=POINT\n",npoints,1,1);
    for (int i=0; i<npoints; ++i){
      for (int var=0; var<varNameVec.size(); ++var)
	fprintf(fp,"%lf ",varDataVec[var][i]);
      fprintf(fp,"\n");
    }

    fclose(fp);
  }

  // ==================================
  // add the pdf intergation coordinate
  // ==================================
  void setIndependentCoor(const string& coorName){
    
#ifdef NEW_CHEM
    assert(npoints>0);

    const char * name = coorName.c_str(); 
    PdfConvolution::buildIndepCoord( npoints, getVarPtr(name), 
				     name, getVarMin(name), getVarMax(name) ); 

    // 
    // now store the state from PdfConvolution 
    // for use later (needs to be a deep copy) 
    //
 
    /*
    DELETE(x_pdf);
    DELETE(idx_pdf);

    memcpy( x_pdf, PdfConvolution::x_pdf, npoints * sizeof(double) ); 
    memcpy( idx_pdf, PdfConvolution::idx, npoints * sizeof(int)    ); 

    xmin_pdf  = PdfConvolution::xmin;
    invLx_pdf = PdfConvolution::invLx; 
    */
#else
    assert(npoints>0);
    assert(x_pdfPair.size() == 0 );

    // set the variable
    assert(x_pdf == NULL);
    x_pdf     = getVarPtr(coorName);
    xmin_pdf  = getVarMin(coorName)-1.0e-16;
    xmax_pdf  = getVarMax(coorName)+1.0e-16;
    assert( xmax_pdf > xmin_pdf);
    invLx_pdf = 1.0 / (xmax_pdf - xmin_pdf);

    for (int i=0; i<npoints; ++i){
      double x = ( x_pdf[i] - xmin_pdf ) * invLx_pdf;
      x_pdfPair.push_back(pairIntDbl(i,x));
    }
    assert( x_pdfPair.size() == npoints);
    // now sort the vector
    std::sort(x_pdfPair.begin(),x_pdfPair.end(),isSecondSmallerPair);
#endif

  }

  // =================================
  // compute convolution with the pdf
  // =================================
  double getValGivenMeanAndVar(const double mean, const double var, const string& varname) {
    
#ifdef NEW_CHEM
    // first update the pdf given the zmean and zvar
    PdfConvolution::updatePdf("Z",mean,var);
    
    // compute the convolution
    return(PdfConvolution::convolveWithPdf(varname, getVarPtr(varname), npoints));
#else
    // first update the pdf given the zmean and zvar
    updatePDF(mean,var);

    // compute the convolution
    return(convolveWithPDF(varname));
#endif

  }

  // ======================================================
  // perturb energy of a flamelet and measure the quantites
  // ======================================================
  void addMixtureERandPerturbedVariables(Mixture* &mixture_in, const double rhoDeltaE) {
    assert(1==2);
  }

  void setMixFromFlamelet(const int i) {
    Cantera::vector_fp y(nspecies);
    size_t nsp = mix->nSpecies();

    for (int isp=0; isp<nsp; isp++) {
      double *Y = getVarPtr("Y_"+mix->speciesName(isp));
      y[isp] = Y[i];
    }

    double *T = getVarPtr("T");
    mix->setState_TPY(T[i], pressure, &y[0]);
  }

  void addMixtureERandPerturbedVariables(const double rhoDeltaE) {
    size_t nsp = mix->nSpecies();
    
    // make sure that at this point "mixture" does not point to null
    assert(mix != NULL);
    // check for consistency
    assert(speciesNameVec.size() == nsp);
    assert(nspecies == nsp);

    // ======================================
    // tolerance for the relative error check
    double tolerance = 0.02; // i.e. 2%
    // ======================================

    // XXX all the checks can be combined to minimize mixture state setting
    // +++++++++++++++++
    // check temperature
    // +++++++++++++++++
    {
      double *H   = getVarPtr("TotalEnthalpy");
      double *T   = getVarPtr("T");
      double Terr = 0.0;
      for (int i = 0; i < npoints; ++i) {
        setMixFromFlamelet(i);
        mix->setState_HP(H[i], pressure);
	double Tf       = T[i];
	double Tm       = mix->temperature();
	double thisTerr = fabs(Tf-Tm) / min(Tf,Tm);
	Terr            = max(Terr,thisTerr);
	T[i]            = Tm;
      }
      double tol = 0.1 * tolerance;
      if (Terr>tol)
	CERR_S("for flamelet " << filename <<
	  "\n large difference between flamelet temperature\n and mixture temperature is observed: \n  Max(error) = "
	       << Terr*100.0 << "%" );
    }
    
    // =============================
    // add constant gass (R_u/M)
    //   R_u: universal gas constant
    // =============================
    {
      double * tmpvec = addVarAndReturnDataPtr("R");
      double *T = getVarPtr("T");

      for (int i = 0; i < npoints; ++i) {
        setMixFromFlamelet(i);
	tmpvec[i] = Cantera::GasConstant / mix->meanMolecularWeight();
      }
    }

    // =======
    // add R*T
    // =======
    {
      double *T       = getVarPtr("T");
      double *R       = getVarPtr("R");
      double *tmpvec = addVarAndReturnDataPtr("RT");
      for (int i = 0; i < npoints; ++i)
	tmpvec[i] = R[i] * T[i];
    }

    // +++++++++++++
    // check density
    // +++++++++++++
    {
      double *RHO   = getVarPtr("rho");
      double *T     = getVarPtr("T");
      double *R     = getVarPtr("R");
      double RHOerr = 0.0;
      for (int i = 0; i < npoints; ++i) {
	double RHOf       = RHO[i];
	double RHOm       = pressure / (R[i]*T[i]);
	double thisRHOerr = fabs(RHOf-RHOm) / min(RHOf,RHOm);
	RHOerr            = max(RHOerr,thisRHOerr);
	RHO[i]            = RHOm;
      }
      double tol = 0.1 * tolerance;
      if (RHOerr>tol) 
	CERR_S("for flamelet " << filename <<
	  "\n large difference between flamelet density\n and mixture density is observed: \n  Max(error) = "
	       << RHOerr*100.0 << "%");
    }
    
    // ================================
    // add internal energy (E = H - RT)
    // ================================
    {
      double *H      = getVarPtr("TotalEnthalpy");
      double *RT     = getVarPtr("RT");
      double *tmpvec = addVarAndReturnDataPtr("e");
      for (int i = 0; i < npoints; ++i)
	tmpvec[i] = H[i] - RT[i];
    }

    // +++++++++
    // check mu
    // +++++++++
    {
      double *MU   = getVarPtr("mu");
      double MUerr = 0.0;
      for (int i = 0; i < npoints; ++i) {
        setMixFromFlamelet(i);
	double MUf       = MU[i];
	double MUm       = trans->viscosity();  //mixtureVec[i].GetMixMul();
	double thisMUerr = fabs(MUf-MUm) / min(MUf,MUm);
	MUerr            = max(MUerr,thisMUerr);
	MU[i]            = MUm;
      }
      double tol = 0.1 * tolerance;
      if (MUerr > tol)
	CERR_S("for flamelet " << filename <<
	  "\n large difference between flamelet viscosity\n and mixture viscosity is observed: \n  Max(error) = "
	       << MUerr*100.0 << "%");
    }

    // +++++++++
    // check cp
    // +++++++++    
    {
      double *T     = getVarPtr("T");
      double *CP    = getVarPtr("cp");
      double CPerr  = 0.0;
      for (int i = 0; i < npoints; ++i) {
        setMixFromFlamelet(i);
	double CPf       = CP[i];
	double CPm       = mix->cp_mass();  //mixtureVec[i].ComputeMixCp(T[i]);
	double thisCPerr = fabs(CPf-CPm) / min(CPf,CPm);
	CPerr            = max(CPerr,thisCPerr);
	CP[i]            = CPm;
      }
      double tol = tolerance;
      if (CPerr > tol)
	CERR_S("for flamelet " << filename <<
	  "\n large difference between flamelet cp\n and mixture cp is observed: \n  Max(error) = "
	       << CPerr*100.0 << "%");
    }

    // ++++++++++++
    // check lambda
    // ++++++++++++
    {
      double *L    = getVarPtr("lambda");
      double Lerr  = 0.0;
      for (int i = 0; i < npoints; ++i) {
        setMixFromFlamelet(i);
//        Cantera::Transport *trans = newTransportMgr(transport, mix);
	double Lf       = L[i];
	double Lm       = trans->thermalConductivity();  //mixtureVec[i].GetMixLambda();
	double thisLerr = fabs(Lf-Lm) / min(Lf,Lm);
	Lerr            = max(Lerr,thisLerr);
	L[i]            = Lm;
      }
      double tol = tolerance;
      if (Lerr>tol) 
	CERR_S("for flamelet " << filename <<
	  "\n large difference between flamelet  Lambda\n and mixture Lambda is observed: \n  Max(error) = "
	       << Lerr*100.0 << "%");
    }

    // ++++++++++++
    // replace locp
    // ++++++++++++
    {
      double *L    = getVarPtr("lambda");
      double *Cp   = getVarPtr("cp");
      double *Locp = getVarPtr("locp");
      for (int i = 0; i < npoints; ++i) {
	Locp[i] = L[i] / Cp[i];
      }
    }

    // ===================
    // perturbed quantites
    // ===================

    // positive perturbations
    double *Ep    = addVarAndReturnDataPtr("Ep");
    double *Tp    = addVarAndReturnDataPtr("Tp");
    double *RTp   = addVarAndReturnDataPtr("RTp");
    double *MUp   = addVarAndReturnDataPtr("MUp");
    double *LOCPp = addVarAndReturnDataPtr("LOCPp");

    // negative perturbations
    double *Em    = addVarAndReturnDataPtr("Em");
    double *Tm    = addVarAndReturnDataPtr("Tm");
    double *RTm   = addVarAndReturnDataPtr("RTm");
    double *MUm   = addVarAndReturnDataPtr("MUm");
    double *LOCPm = addVarAndReturnDataPtr("LOCPm");

    // get the original vars
    double *T    = getVarPtr("T");
    double *RHO  = getVarPtr("rho");
    double *R    = getVarPtr("R");
    double *E    = getVarPtr("e");
   
    // for all the points
    for (int i=0; i<npoints; ++i ) {
      setMixFromFlamelet(i);
//      Cantera::Transport *trans = newTransportMgr(transport, mix);

      // perturbed thermo peroperties
      Ep[i] = E[i] + rhoDeltaE/RHO[i];
//      mix->setState_HP(Ep[i]+pressure/RHO[i], pressure);  // assuming dh = de
      setMixState_UP(Ep[i], pressure);
      Tp[i]    = mix->temperature();
      MUp[i]   = trans->viscosity();
      LOCPp[i] = trans->thermalConductivity() / mix->cp_mass();
      RTp[i] = Tp[i] * R[i];  // CFPV assumes Yk do not change

      Em[i] = E[i] - rhoDeltaE/RHO[i];
//      mix->setState_HP(Em[i]+pressure/RHO[i], pressure);  // assuming dh = de
      setMixState_UP(Em[i], pressure);
      Tm[i]    = mix->temperature();
      MUm[i]   = trans->viscosity();
      LOCPm[i] = trans->thermalConductivity() / mix->cp_mass();
      RTm[i] = Tm[i] * R[i];  // CFPV assumes Yk do not change
    }

  }

  void setMixState_UP(const double e, const double p, const double tol = 1.0e-4) {
    double Terror = 1.0e20;

    while (Terror >= tol) {
      double Told = mix->temperature();
      mix->setState_HP(e+p/mix->density(), p);
      Terror = fabs(mix->temperature() - Told);
    }
  }

public:

  // ==============
  // check for NAN
  // ==============
  void checkForNaN() {
    for (int j=0; j<varNameVec.size(); ++j) {
      double* tmpVec = getVarPtr(varNameVec[j]);
      for (int i=0; i<npoints; ++i) 
	if ( tmpVec[i] != tmpVec[i] )
	  CERR_S("in flamelet file "  << filename << " variable " << varNameVec[j] << " is NaN at location i = " << i);
    }
  }

  // ======================================
  // dump variables range
  // ======================================
  void dumpVarsRange(){
    cout << "\n===========================================" << endl;
    cout << " flamelet file " << filename << endl;
    for (int i=0; i<varNameVec.size(); ++i)
      cout << "  > " << varNameVec[i] << " range: (" << getVarMin(varNameVec[i]) << "," << getVarMax(varNameVec[i]) << ")" << endl;
    cout << "===========================================\n" << endl; 
  }

  // ======================
  // read the flamelet file
  // ======================  
  void readFlameletFile(const string& filename) {

    ifstream ifp;

    // open the file
    ifp.open(filename.c_str());
    if ( ifp.fail() ) 
      CERR_S("could not open flamelet file: " << filename )
        
    // ================
    // header
    // ================
    {
      // make sure the size is zero
      assert(headerVec.size() == 0);
      string buffer_str; 
      while(getline(ifp,buffer_str)) {
	string header, dummy;      
	istringstream buf(buffer_str);
	buf >> header;
	if ( !header.compare("body") )
	  break;
	else {
	  headerVec.push_back(buffer_str);
	  // initialize grids size and number of species
	  if ( !header.compare("numOfSpecies") )
	    buf >> dummy >> nspecies;
	  else if ( !header.compare("gridPoints") )
	    buf >> dummy >> npoints;
	}
      }
    }//header

    //=================
    // XXX HACK TO GET Z2 from filename ..
    //================
    { 
      string::size_type found = filename.find("Z2_"); 
      if (found != filename.npos) { 
	string str2 = filename.substr(found+3,found+10); 
	z2          = atof(str2.c_str()); 
	//cout << " FOUND Z2 = " << z2 << endl;
      } 
      else { 
	z2 = 0.0; 
      }
    }//Z2 hack

    // ===============
    // variables
    // ===============  
    string buffer_str;
    // number of lines  
    int tmpint = npoints % 5;      
    int nLines = npoints/5;
    if ( tmpint > 0 ) ++nLines;
    
    while(getline(ifp,buffer_str)) {
      int count = 0;
      // variable name
      istringstream buf(buffer_str);
      string dummy_str;
      buf >> dummy_str;
      // if we are at the trailer location
      if (dummy_str == "trailer")
	      break;

      if (checkParam("DUMP_NAMES")) 
        cout << dummy_str << endl;
      // 
      // decide if we are going to read 
      // the data 
      // 
      bool isStoredVar;
      //if ( checkParam("CHEMTABLE.VARS")) { 
      if (Param * param = getParam("CHEMTABLE.VARS")) { 
	      isStoredVar     = false;
        
        for (int i = 0; i < param->size(); i++) { 
          if (param->getString(i) == dummy_str) { 
	          isStoredVar = true;
	          break;
	        }
	      }//i 
      }
      else 
	      isStoredVar = true; // default is yes

      if ( isStoredVar ) { 
	      varNameVec.push_back(dummy_str);
	
	      // and then the data
	      double * tmpvec = new double[npoints];
	      varDataVec.push_back(tmpvec);
	
	      for ( int iLines=0; iLines<nLines; ++iLines ) {
	        getline(ifp,buffer_str);
	        istringstream thisbuf(buffer_str);
	        for ( int iCol=0; iCol<5; ++iCol) {
	          double tmp;
	          if ( !(thisbuf >> tmp) )  break;
	          tmpvec[count] = tmp;
	          ++count;	
	        }
	      }//ilines 

	      assert(count==npoints);
      } 
      else { 
	      // skip nlines .. 
	      for (int iLines=0; iLines<nLines; ++iLines)
	      getline(ifp,buffer_str);
      }
    } // while(getline(ifp,buffer_str))
    
    // done with the file
    ifp.close();

    // check the both variable name and variable data vectors have the same size
    assert(varNameVec.size() == varDataVec.size());
    
  }//readFlameletFile()

  // =============================
  // clean up names from the table
  // =============================
  void cleanupNamesAndInitSpeciesNames() {

    // ==================
    // mass fractions
    // ==================
    int tmpint = 0;
    assert(speciesNameVec.size()==0);
    for (int ivar=0; ivar<varNameVec.size(); ++ivar) {
      size_t pos;
      pos = varNameVec[ivar].find("massfraction-");
      if ((int)pos != -1) {
	varNameVec[ivar].erase(pos,13);
	speciesNameVec.push_back(varNameVec[ivar]);	
	string name = "Y_" + varNameVec[ivar];
	varNameVec[ivar] = name;
	++tmpint;
      }
    }
    nspecies = tmpint;
    //assert(nspecies==tmpint);
    assert(nspecies==speciesNameVec.size());

    // ==================
    // mole fractions
    // ==================
    tmpint = 0;
    for (int ivar=0; ivar<varNameVec.size(); ++ivar) {
      size_t pos;
      pos = varNameVec[ivar].find("molefraction-");
      if ((int)pos != -1) {
	varNameVec[ivar].erase(pos,13);
	string name = "X_" + varNameVec[ivar];
	varNameVec[ivar] = name;
	++tmpint;
      }
    }
    
    // ==============
    // source terms
    // ==============
    for (int ivar=0; ivar<varNameVec.size(); ++ivar) {
      string tmp = varNameVec[ivar];
      size_t pos;
      pos = tmp.find("ProdRate-");
      if ((int)pos != -1) {
        tmp.erase(pos,9);
	varNameVec[ivar] = "src_" + tmp;
      }
      // 01/2013 -- backward compatibility
      // older versions of flamemaster omitted the "-" from "ProdRate-"
      pos = tmp.find("ProdRate");
      if ((int)pos != -1) {
        tmp.erase(pos,8);
	varNameVec[ivar] = "src_" + tmp;
      }
    }

    // ========================================
    // density, temperature, and mass flow rate
    // ========================================

    for (int ivar = 0; ivar < varNameVec.size(); ++ivar) {
      string tmp = varNameVec[ivar];
      if ( tmp == "density" )
	varNameVec[ivar] = "rho";
      else if ( tmp == "temperature" )
	varNameVec[ivar] = "T";
      else if ( tmp == "massflowrate" )
	varNameVec[ivar] = "mdot";
    }

    // ================================
    // NOx, enthalpy, cp, radiation, heat release
    // ================================

    for (int ivar = 0; ivar < varNameVec.size(); ++ivar) {
      string tmp = varNameVec[ivar];
      if ( tmp == "TotalEnthalpy2" )
        varNameVec[ivar] = "h";
      else if ( tmp == "src_Pos-NO" )
        varNameVec[ivar] = "src_posNO";
      else if ( tmp == "src_Neg-NO" )
        varNameVec[ivar] = "src_negNO";
      else if ( tmp == "GasRadiation" )
        // 01/2013 -- backward compatibility
        // "GasRadiation" superceded by "RadiationSource" in new flamemaster
        varNameVec[ivar] = "src_rad";
      else if ( tmp == "RadiationSource" )
        varNameVec[ivar] = "src_rad";
      else if ( tmp == "Cp" || tmp == "cp" )
        varNameVec[ivar] = "cp";
      else if ( tmp == "HeatRelease" )
        varNameVec[ivar] = "heatrelease";
      else if ( tmp == "lambdaOverCp" )
        varNameVec[ivar] = "locp";
      else if ( tmp == "entropy" )
        varNameVec[ivar] = "s";

      // info for droplet evaporation model
      else if ( tmp == "Tref" )
        varNameVec[ivar] = "lsp_Tref";
      else if ( tmp == "rho_fuel_Tref" )
        varNameVec[ivar] = "lsp_rho_fuel";
      else if ( tmp == "cp_fuel_Tref" )
        varNameVec[ivar] = "lsp_cp_fuel";
      else if ( tmp == "mu_fuel_Tref" )
        varNameVec[ivar] = "lsp_mu_fuel";
      else if ( tmp == "cond_fuel_Tref" )
        varNameVec[ivar] = "lsp_cond_fuel";
      else if ( tmp == "rho_nonfuel_Tref" )
        varNameVec[ivar] = "lsp_rho_nonfuel";
      else if ( tmp == "cp_nonfuel_Tref" )
        varNameVec[ivar] = "lsp_cp_nonfuel";
      else if ( tmp == "mu_nonfuel_Tref" )
        varNameVec[ivar] = "lsp_mu_nonfuel";
      else if ( tmp == "cond_nonfuel_Tref" )
        varNameVec[ivar] = "lsp_cond_nonfuel";
      else if ( tmp == "mw_nonfuel" )
        varNameVec[ivar] = "lsp_mw_nonfuel";

      // molecular weight
      else if ( tmp == "W" )
        varNameVec[ivar] = "mw";

    }

  }

  // ================================
  // clean up negative mass fractions
  // ================================
  void cleanupNegativeMassFractions() {

    // ==================
    // mass fractions
    // ==================
    for (int ivar=0; ivar<varNameVec.size(); ++ivar) {
      size_t pos;
      pos = varNameVec[ivar].find("Y_");
      if ((int)pos != -1) {
	for (int i=0; i<npoints; ++i)
	  if ( varDataVec[ivar][i] < 0.0 ) {
	    //cout << "WARNING: in flamelet " << filename << " mass fraction species" << 
	    //  varNameVec[ivar] << " = " <<  varDataVec[ivar][i] << " set to zero" << endl;
	    varDataVec[ivar][i] = 0.0;
	  }
      }
    }

  }
  
  // ================================================
  // add progress variable and its source term to the 
  // table. note that the source term in progress 
  // variable equation is not multiplied by density
  // ================================================
  void addProgAndSrcProg() {

    assert(progSpeciesNameVec.size() == 0 );
    assert(progSpeciesWeightVec.size() == 0 );
   
    //if (checkParam("PROGRESS_VARIABLE.SPECIES+WEIGHTS")) 
    if (Param * param = getParam("PROGRESS_VARIABLE.SPECIES+WEIGHTS"))
      for (int i = 0; i < param->size()/2; ++i) { 
        progSpeciesNameVec.push_back(param->getString(2*i));
        progSpeciesWeightVec.push_back(param->getDouble(2*i+1));
      }
    else 
      CERR_S("could not find a consistent definition for progress variable\n parameter PROGRESS_VARIABLE.SPECIES+WEIGHTS should be provided");
    assert(progSpeciesNameVec.size() == progSpeciesWeightVec.size() );

    // add the progress variable
    {
      double * tmpvec = addVarAndReturnDataPtr("prog");
      // loop over all the progress variable species
      for (int iprog=0; iprog<progSpeciesNameVec.size(); ++iprog){
	double * thisY = getVarPtr("Y_"+progSpeciesNameVec[iprog]);
	for (int i=0; i<npoints; ++i)
	  tmpvec[i] += progSpeciesWeightVec[iprog] * thisY[i];
      }

      // normalization of tiny numbers .. machine prec issues
      Cmax = -1.0e+16; 
      for (int i=0; i < npoints; ++i) {
	if ( fabs(tmpvec[i]) < 1.0e-16) 
	  tmpvec[i] = 0.0; 
	Cmax = fmax(tmpvec[i], Cmax); 
      }

    }

    // and the source term
    {  
      double * tmpvec = addVarAndReturnDataPtr("src_prog"); 
      // loop over all the progress variable species
      for (int iprog=0; iprog<progSpeciesNameVec.size(); ++iprog){
	double * this_src_Y = getVarPtr("src_"+progSpeciesNameVec[iprog]);
	for (int i=0; i<npoints; ++i)
	  tmpvec[i] += progSpeciesWeightVec[iprog] * this_src_Y[i];
      }
    }

  }

  // ===================================================
  // divide source terms by density
  // ===================================================
  void divSrcByRho() {
    
    double * tmpvec = NULL;
    double * rho = getVarPtr("rho");
    for (int ivar=0; ivar<varNameVec.size(); ++ivar) {
      string tmp = varNameVec[ivar];
      if ( tmp == "src_NO" ) 
         tmpvec = varDataVec[ivar];
      else if ( tmp == "src_NO2" ) 
         tmpvec = varDataVec[ivar];
      else if ( tmp == "src_negNO" ) 
         tmpvec = varDataVec[ivar];
      else if ( tmp == "src_posNO" ) 
         tmpvec = varDataVec[ivar];
      else if ( tmp == "src_rad" ) 
         tmpvec = varDataVec[ivar];
      else if ( tmp == "src_prog" ) 
         tmpvec = varDataVec[ivar];
      else
         continue;

      assert( tmpvec !=NULL);
  
      for (int i=0; i<npoints; ++i) {
         tmpvec[i] /= rho[i];
      }
    }
  }

  // ===================================================
  // divide source terms by heat capacity (for T-eqn)
  // ===================================================
  void divSrcByCp() {

    double * tmpvec = NULL;
    double * cp = getVarPtr("cp");
    for (int ivar=0; ivar<varNameVec.size(); ++ivar) {
      string tmp = varNameVec[ivar];
      if ( tmp == "src_rad" )
         tmpvec = varDataVec[ivar];
      else
         continue;

      assert( tmpvec !=NULL);

      for (int i=0; i<npoints; ++i) {
         tmpvec[i] /= cp[i];
      }
    }
  }

  void addGamma() {
    
    double R = kR; // universal gas constant
    double * gamma = addVarAndReturnDataPtr("gamma"); 
    double * cp = getVarPtr("cp"); 
    double * mw = getVarPtr("mw"); 
    for (int i=0; i<npoints; ++i) {
      gamma[i] = 1.0/(1.0-R/mw[i]/cp[i]); // gamma=Cp/Cv, Cv=Cp-R
    }

  }

  // ================
  // add lambda / cp
  // ================
  void addLocp() {
    
    if ( checkVar("locp") )
      return;

    double * tmpvec = addVarAndReturnDataPtr("locp"); 
    // get lambda and cp
    double * lambda =  getVarPtr("lambda"); 
    double * cp     =  getVarPtr("cp"); 
    for (int i=0; i<npoints; ++i) {
      assert(cp[i]>0.0);
      tmpvec[i] = lambda[i] / cp[i];
    }

  }
   
  // ====================================
  // compute the convolution with the pdf
  // ===================================
  double convolveWithPDF(const string& varname) {

    // get the variable
    double * var; // = getVarPtr(varname);
    
    if (varHash.size() != varDataVec.size()) {
      varHash.clear();
      assert(varHash.size() == 0);
      for (int ivar = 0; ivar < varNameVec.size(); ++ivar) 
	varHash[varNameVec[ivar]] = varDataVec[ivar];
      assert(varHash.size() == varDataVec.size());
    }

    map<string,double*>::const_iterator found = varHash.find(varname);
    //unordered_map<string,double*>::const_iterator found = varHash.find(varname);
    if (found != varHash.end()) var = found->second;
    else CERR_S("could not find variable " << varname << " in flamelet");
    
    double sum = 0.0;
    if ( varname == "rho" ) {
      for (int i=0; i<npoints; ++i) {
          sum += pdf[i] / var[i];
//          if (var[i]<1.0e-2)
//            cout << "var is almost 0: " << var[i] << "; " << pdf[i] << endl;
      }
//      cout << "varName = " << varname << "; sum = " << sum << endl;  // XXX ars
      assert(sum>0.0);
      return(1.0/sum);
    }
    else {
      for (int i=0; i<npoints; ++i)
	sum += pdf[i] * var[i];
      return(sum);
    }
  
  }
    
  // ===================================
  // update pdf for the current flamelet
  // ===================================
  void updatePDF(const double mean, const double var) {

    assert(pdf!=NULL);

    switch (pdf_function) {
    case BETA_PDF:
      updateBetaPDF(mean,var);
      break;
    default:
      CERR_S("unrecognized pdf function");
    }
    
  }

  // ========================================
  // update beta pdf for the current flamelet
  // ========================================
  void  updateBetaPDF(const double mean_in, const double var_in) {

    // ==============================================================
    // Beta distribution is parametrized by a and b so that
    //   pdf = Gamma(a+b)/(Gamma(a)*Gamma(b)) * x^(a-1) * (1-x)^(b-1)
    // where
    //   mean = a / (a+b)
    //   variance = (a*b) / ((a+b)^2 * (a+b+1))
    // which gives
    //   a = mean * ( (mean*(1-mean))/variance - 1)
    //   b = a / mean - a
    // The pdf is computed as e^[log(pdf)] to increase accuracy
    // This procedure requires the variables to be normalized 
    // between 0 and 1. 
    // ==============================================================   

#ifdef NEW_CHEM
    assert(0); // this has been moved to PdfConvolution.hpp
#else
    // assert that mean and variance are consistent
    //assert(mean>=0.0);
    //assert(mean<=1.0);
    //assert(var>=0.0);
    //assert(var<=(mean*(1.0-mean)));


    assert(x_pdf != NULL);
    assert(x_pdfPair.size() == npoints);

    // rescale mean and variance
    double mean = (mean_in-xmin_pdf) * invLx_pdf;
    double var  = var_in * invLx_pdf*invLx_pdf;

    // recompute the Beta PDF only if mean and var have changed
    if ((mean != pdfMean) || (var != pdfVar)) {

      // reset pdfMean and pdfVar
      pdfMean = mean;
      pdfVar  = var;
      
      // set the initial value to zero
      for (int i=0; i<npoints; ++i)
        pdf[i] = 0.0;
      

      // limiting behaviour
      // if var is zero or mean is zero or one
      // the pdf becomes a delta function
      // define with a tolerance

      double tol = 1.0e-10;
      if ( mean < tol ) {
	
	pdf[0] = 1.0;
	return;
      
      }
      else if ( (1.0-mean) < tol ) {
	
	pdf[npoints-1] = 1.0;
	return;

      }      
      else if ( var <= tol ) {
	
	int count = 0;
	for (int i=0; i<npoints-1; ++i)
	  if ( mean<(x_pdfPair[i+1].second) ) 
	    break;
	  else
	    ++count;
	if ( count == (npoints-1) ) 
	  CERR_S("mean is not in the x_pdf range");
	assert( (x_pdfPair[count+1].second) > (x_pdfPair[count].second) );
	assert( mean >= (x_pdfPair[count].second) );
	double wgt = ( x_pdfPair[count+1].second - mean ) / 
	  ( x_pdfPair[count+1].second - x_pdfPair[count].second );
	pdf[x_pdfPair[count].first]   = wgt;
	assert(wgt<=1.0);
	pdf[x_pdfPair[count+1].first] = 1.0 - wgt;       
	return;
      
      }
      // Impossible cases (i.e., var > mean*(1-mean)): two delta pdf at x_pdf=0 and x_pdf=1 
      // note that we already checked for  var < mean*(1-mean)   
      else if ( var >= mean * (1.0 - mean)) {
      
	pdf[x_pdfPair[0].first]         = 1.0 - mean;
        pdf[x_pdfPair[npoints-1].first] = mean;
        return;
      
      }

      else {

	double a      = mean * (mean * (1.0 - mean) / var - 1.0);
	assert( a != 0.0 );
	double b      = a / mean - a;
	assert( b != 0.0 );
	double factor = lgamma(a+b) - lgamma(a) - lgamma(b);
	// Left BC: explicit integration
	double dx = 0.5 * (x_pdfPair[1].second - x_pdfPair[0].second);
	if ( dx < tol ) 
	  pdf[x_pdfPair[0].first] = 0.0;
	else {
	  double tmp = a * log(dx) + factor;
	  pdf[x_pdfPair[0].first] = exp(tmp) / a;
	}
	// Right BC: explicit integration
	dx = 0.5 * (x_pdfPair[npoints-1].second - x_pdfPair[npoints-2].second);
	if ( dx < tol ) 
	  pdf[x_pdfPair[npoints-1].first] = 0.0;
	else {
	  double tmp = b * log(dx) + factor;
	  pdf[x_pdfPair[npoints-1].first] = exp(tmp) / b;
	}
	// other points
	for (int i=1; i<npoints-1; ++i) {
	  dx     = 0.5 * (x_pdfPair[i+1].second - x_pdfPair[i-1].second);
	  if ( dx < tol )
	    pdf[x_pdfPair[i].first] = 0.0;
	  else {
	    double tmp = (a - 1.0) * log(x_pdfPair[i].second) + 
	      (b - 1.0) * log(1.0 - x_pdfPair[i].second) + factor;
	    pdf[x_pdfPair[i].first] = exp(tmp) * dx;
	  }
	}

      }

      // normalize pdf
      double sum = 0.0;
      for  (int i=0; i<npoints; ++i)
        sum += pdf[i];
      assert(sum>0.0);
      double one_over_sum = 1.0 / sum;
      for  (int i=0; i<npoints; ++i)
        pdf[i] *= one_over_sum;

    } //  if ((mean != pdfMean) || (var != pdfVar))
#endif  
  }

};

// =====================
// premixed flamelet
// =====================
class PremixedFlamelet: public Flamelet {

public: 
 
  double phi;
  double Zst;
  double Tunburnt;
  double laminarFlameSpeed;
  double laminarFlameThickness;
  vector<string> unburntSpeciesVec;
  vector<double> unburntMassfractionVec;
  
  PremixedFlamelet() : Flamelet() {
    //cout << "PremixedFlamelet():" << endl;

    // initialize the values
    phi                   = -1.0e+20;
    Zst                   = -1.0e+20;
    Tmax                  = -1.0e+20;
    Tunburnt              = -1.0e+20;
    pressure              = -1.0e+20;
    Cmin                  = -1.0e+20;
    Cmax                  = -1.0e+20;
    laminarFlameSpeed     = -1.0e+20;
    laminarFlameThickness = -1.0e+20;

    // XXX not necessary when using Cantera
//    // later we will use the NULL for checking
//    mixture = NULL;
  }
  
  virtual ~PremixedFlamelet(){
    //cout << "~PremixedFlamelet()" << endl;
  }

  // ==================================
  // write flamelet into a tecplot file
  // ==================================
  void writeTecplot(const string& filename){
    writeTecplotfile(filename);
  }

  // ==========
  // initialize
  // ==========
  void init(const string& filename){

    this->filename = filename;

    Flamelet::init(filename);

    // add rho * src_prog
    double *src = addVarAndReturnDataPtr("rho_src_prog");
    double *rho = getVarPtr("rho");
    double *src_prog = getVarPtr("src_prog");
    for (int i=0; i<npoints; ++i)
      src[i] = rho[i] * src_prog[i];

    // cumulative integral of (rho * csrc) dx
    double *x = getVarPtr("y");
    double* int_src = addVarAndReturnDataPtr("int_rho_src");
    int_src[0] = 0.0;
    for (int i=1; i<npoints; ++i)
      int_src[i] = int_src[i-1] + 0.5*(x[i]-x[i-1])*(src[i]+src[i-1]);

    // add heat release integral
    double *q = getVarPtr("heatrelease");
    double *int_q = addVarAndReturnDataPtr("int_heatrelease");
    int_q[0] = 0.0;
    for (int i = 1; i < npoints; ++i)
      int_q[i] = int_q[i-1] + 0.5*(x[i]-x[i-1])*(q[i]+q[i-1]);

    assert( unburntSpeciesVec.size()      == 0 );
    assert( unburntMassfractionVec.size() == 0 );
    parseHeaderAndInitParams();

    setIndependentCoor("prog");

    // compute min and max progress variable
    Cmin = getVarMin("prog");
    //assert(Cmin<=1.0);
    //assert(Cmin>=0.0);
    Cmax = getVarMax("prog");
    //assert(Cmax<=1.0);
    //assert(Cmax>=0.0);

    // do some checks
    assert( Tmax                     >  0.0 ); 
    assert( Tunburnt                 >  0.0 ); 
    assert( pressure                 >  0.0 );
    //assert( Cmin                     >= 0.0 );
    //assert( Cmax                     <= 1.0 );
    assert( Cmin                     <= Cmax);      
    assert( laminarFlameSpeed        >= 0.0 ); 
    assert( laminarFlameThickness    >  0.0 ); 
    assert( unburntSpeciesVec.size() >  0   );
    assert( unburntSpeciesVec.size() == unburntMassfractionVec.size() );

    // dump a tecplot file for the flamelet
    if ( checkParam("FLAMELET.TECPLOT") )
      writeTecplotfile(filename);

  }

  // ========================
  // dump info on the screen
  // ========================
  void info() {
    
    cout << " > Premixed Flamelets info:" << endl;
    cout << "       progress variable definition:      C = ";
    for (int i=0; i<progSpeciesNameVec.size()-1; ++i) 
      cout << progSpeciesWeightVec[i] << "*" << progSpeciesNameVec[i] << " + ";
    int tmpint = progSpeciesNameVec.size()-1;
    cout << progSpeciesWeightVec[tmpint] << "*" << progSpeciesNameVec[tmpint] << endl;
    cout << "       background pressure (Pa):          " << pressure << endl;
    cout << "       unburnt temperature (K):           " << Tunburnt << endl;
    cout << "       number of species:                 " << nspecies << endl;
    cout << "       unburnt species:                   ";
    for (int i=0; i<unburntSpeciesVec.size()-1; ++i)
      cout << unburntSpeciesVec[i] << ", ";
    cout << unburntSpeciesVec[unburntSpeciesVec.size()-1] << endl;
    
  }

  // ===================
  // add NOx source term
  // ===================
  void addNOxSrc() {
    
    if ( checkVar("src_NOX") )
      return;

    double *tmpvec = addVarAndReturnDataPtr("src_NOX");
    double *srcNO  = getVarPtr("src_NO");
    double *srcNO2 = getVarPtr("src_NO2");
    for (int i=0; i<npoints; ++i) {
      tmpvec[i] = srcNO[i] + srcNO2[i];
    }

  }

  void addNOxSrcSplit() {
    
    // NOTE: this implementation is not robust to kinetic mechanisms that
    //       don't include NO, NO2 or where they are named differently
    
    double *s_t = addVarAndReturnDataPtr("src_nox_therm"); 
    double *I_p = addVarAndReturnDataPtr("int_src_nox_prompt"); 

    double *x      = getVarPtr("y");
    double *rho    = getVarPtr("rho");
    double *prog   = getVarPtr("prog"); 
    double *srcNO  = getVarPtr("src_NO"); 
    double *srcNO2 = getVarPtr("src_NO2"); 
    double *s_p    = new double[npoints];

    // split "prompt" and "thermal" NOx
    const double c_pt  = 0.99; // splitting location
    const double delta = 0.01; // transition width

    for (int i=0; i<npoints; ++i) {
      const double cc = (prog[i]-prog[0])/(prog[npoints-1]-prog[0]+1e-15);
      double wt;
      if      (cc < c_pt-delta/2.0) 
        wt = 0.0; 
      else if (cc > c_pt+delta/2.0) 
        wt = 1.0; 
      else
        wt = 0.5*(1.0+sin(M_PI*(cc-c_pt)/delta));

      const double rho_src = rho[i]*(srcNO[i]+srcNO2[i]); // total   [=] (kg NOx)/m3/sec
      const double st = rho_src * wt;                     // thermal [=] (kg NOx)/m3/sec
      s_p[i] = rho_src - st;                              // prompt  [=] (kg NOx)/m3/sec
      s_t[i] = st/rho[i];                                 // thermal [=] 1/sec
    }

    // cumulative integral of (rho * src_prompt) dx
    I_p[0] = 0.0; 
    for (int i=1; i<npoints; ++i)
      I_p[i] = I_p[i-1] + 0.5*(x[i]-x[i-1])*(s_p[i]+s_p[i-1]); 

    delete[] s_p; 

  }

  void addNOxSrcSplitHCN() {
    
    // NOTE: this implementation is not robust to kinetic mechanisms that
    //       don't include NO, NO2, HCN or where they are named differently
    
    double *s_t = addVarAndReturnDataPtr("src_nox_therm"); 
    double *I_p = addVarAndReturnDataPtr("int_src_nox_prompt"); 

    double *x      = getVarPtr("y");
    double *rho    = getVarPtr("rho");
    double *prog   = getVarPtr("prog"); 
    double *srcNO  = getVarPtr("src_NO"); 
    double *srcNO2 = getVarPtr("src_NO2");
    double *yHCN   = getVarPtr("Y_HCN");

    // use HCN concentration to differentiate prompt and thermal NOx
    const double Ymax = getVarMax("Y_HCN");
    const double Ylim = Ymax/10.0;
    int    imax;
    for (imax=0; imax<npoints-1; ++imax)
      if (yHCN[imax]==Ymax) break;

    int ilim; 
    for (ilim=imax; ilim<npoints-1; ++ilim)
      if (yHCN[ilim+1]<Ylim) break; 
    ilim = min(ilim,npoints-2); 

    // split prompt and thermal
    double *s_p = new double[npoints];
    for (int i=0; i<npoints; ++i) {

      const double rho_src = rho[i]*(srcNO[i]+srcNO2[i]);
      if (i>ilim) { 
        s_p[i] = 0.0;
        s_t[i] = rho_src/rho[i]; // thermal [=] 1/sec
      }
      else { 
        s_p[i] = rho_src; // prompt [=] (kg NOx)/m3/sec
        s_t[i] = 0.0; 
      }
    }
    
    // cumulative integral of (rho * src_prompt) dx
    I_p[0] = 0.0; 
    for (int i=1; i<npoints; ++i) 
      I_p[i] = I_p[i-1] + 0.5*(x[i]-x[i-1])*(s_p[i]+s_p[i-1]); 

    delete[] s_p; 

  }

  // =================
  // add mole fraction
  // =================
  void addMoleFraction(string species) {

    double *X   = addVarAndReturnDataPtr("X_"+species);
    double *Y   = getVarPtr("Y_"+species);
    double *mw  = getVarPtr("mw");

    // get species molecular weight
    Species sp;
    char *name_char = new char[species.length() + 1];
    strcpy(name_char, species.c_str());
    sp.SetSpeciesName(name_char);
    sp.ComputeSpeciesConstants(1.0);
    double W_sp = sp.GetSpeciesMolarMass();

    for (int i=0; i<npoints; ++i) {
      X[i] = Y[i]*mw[i]/W_sp;
    }

    delete[] name_char;
  }

  // =====================
  // add dry mole fraction
  // =====================
  void addDryMoleFraction(string species) {

    double *X_dry = addVarAndReturnDataPtr("X_"+species+"_DRY");
    double *Y     = getVarPtr("Y_"+species);
    double *mw    = getVarPtr("mw");
    double *Y_H2O = getVarPtr("Y_H2O");

    // get species molecular weight
    Species sp;
    char *name_char = new char[species.length() + 1];
    strcpy(name_char, species.c_str());
    sp.SetSpeciesName(name_char);
    sp.ComputeSpeciesConstants(1.0);
    double W_sp  = sp.GetSpeciesMolarMass();
    double W_H2O = 2.0*AWT_H + AWT_O;

    for (int i=0; i<npoints; ++i) {
      const double X     = Y[i]    *mw[i]/W_sp;
      const double X_H2O = Y_H2O[i]*mw[i]/W_H2O;
      X_dry[i] = X / (1.0 - X_H2O);
    }

    delete[] name_char;

  }

  // =======================
  // add laminar flame speed
  // =======================
  void addLaminarFlameSpeed() {

    if ( checkVar("sL") )
      return;

    double * tmpvec = addVarAndReturnDataPtr("sL"); 
    for (int i=0; i<npoints; ++i) {
      tmpvec[i] = laminarFlameSpeed;

    }
  }

  // ===========================
  // add laminar flame thickness
  // ===========================
  void addLaminarFlameThickness() {

    if ( checkVar("lF") )
      return;

    double * tmpvec = addVarAndReturnDataPtr("lF"); 
    for (int i=0; i<npoints; ++i) {
      tmpvec[i] = laminarFlameThickness;
    }
  }

private:

  // =====================================
  // parse the header of the flamelet file
  // =====================================
  void parseHeaderAndInitParams() {
    
    for (int i=0; i<headerVec.size(); ++i){
      
      istringstream buf(headerVec[i]);    
      string header;
      buf >> header;
      string dummy;
 
      if ( !header.compare("body") )
	CERR_S("should not be here")
      else if ( !header.compare("pressure") ) {
	buf >> dummy >> pressure;
	pressure *= 100000.0; // convert from bar to Pa
      } 
      
      else if ( !header.compare("fuel-air-equivalence-ratio") )
	buf >> dummy >> phi;
      
      else if ( !header.compare("mixture-fraction-stoichiometric") )
        buf >> dummy >> Zst;

      else if ( !header.compare("Tmax") )
	buf >> dummy >> Tmax;
      
      else if ( !header.compare("Temperature") )
	buf >> dummy >> Tunburnt;
      
      else if ( !header.compare("burningVelocity") ) { 
	buf >> dummy >> laminarFlameSpeed;
        laminarFlameSpeed *= 0.01; // convert from cm/s to m/s
      }
      
      else if ( !header.compare("FlameThickness") )
	buf >> dummy >> laminarFlameThickness;
      
      else if ( !header.compare("numOfSpecies") ) {
	int inttmp = 0;
	buf >> dummy >> inttmp;
	assert(nspecies == inttmp);
      }
      
      else if ( !header.compare("gridPoints") ) {
	int inttmp = 0;
	buf >> dummy >> inttmp;
	assert(npoints == inttmp);
      }      
      
      else {
	size_t pos = header.find("Massfraction-");
	if ( (int)pos != -1 ) {
	  header.erase(pos,13);
	  unburntSpeciesVec.push_back(header);
	  double tmp;
	  buf >> dummy >> tmp;
	  unburntMassfractionVec.push_back(tmp);
	}
      }

    } //  for (int i=0; i<headerVec.size(); ++i)

  }

};

// =====================
// non premixed flamelet
// =====================
class NonPremixedFlamelet: public Flamelet {

public:
  
  double chi_st;
  double Tfuel;
  double Toxi;
 
public:
  
  NonPremixedFlamelet() : Flamelet() {

    //cout << "NonPremixedFlamelet():" << endl;

    // initialize the values
    chi_st    = -1.0e+20;
    Tmax      = -1.0e+20;
    pressure  = -1.0e+20;
    Tfuel     = -1.0e+20;
    Toxi      = -1.0e+20;
    Cmin      = -1.0e+20;
    Cmax      = -1.0e+20;

    // XXX not necessary when using Cantera
//    // later we will use the NULL for checking
//    mixture = NULL;
  }
  
  virtual ~NonPremixedFlamelet(){
    //cout << "~NonPremixedFlamelet()" << endl;
  }

  // ==========
  // initialize
  // ==========
  void init(const string& filename){

    this->filename = filename;

    // initialize flamelet
    Flamelet::init(filename);

    // parse the deader and set parameters
    parseHeaderAndInitParams();

    // set the independent coordinate for convolution
    setIndependentCoor("Z");

    // compute min and max progress variable
    Cmin = getVarMin("prog");
    assert(Cmin<=1.0);
    if( Cmin<0 ) {
      cout << "WARNING: negative progress variable " << Cmin << " in flamelet file: " << filename << endl;
      writeTecplotfile(filename);
    }
    Cmax = getVarMax("prog");
    //assert(Cmax<=1.0);
    //assert(Cmax>=0.0);
   
    // do some checks
//    assert( chi_st   > 0.0 );
    assert( Tmax     > 0.0 ); 
    assert( pressure > 0.0 );
    assert( Tfuel    > 0.0 );
    assert( Toxi     > 0.0 );
    //assert( Cmin     >= 0.0 );
    //assert( Cmax     <= 1.0 );
    assert( Cmin     <= Cmax);      

    // dump a tecplot file for the flamelet
    if ( checkParam("FLAMELET.TECPLOT") )
      writeTecplotfile(filename);

  }

  // ========================
  // dump info on the screen
  // ========================
  void info() {
    
    cout << " > Non-premixed Flamelets info:" << endl;
    string pdf_name;
    if ( pdf_function == BETA_PDF ) 
      pdf_name = "BETA_FUNCTION";
    else 
      CERR_S("Error: unknown pdf function");
    cout << "       presumed pdf:                      " << pdf_name << endl;
    cout << "       progress variable definition:      C = ";
    for (int i=0; i<progSpeciesNameVec.size()-1; ++i) 
      cout << progSpeciesWeightVec[i] << "*" << progSpeciesNameVec[i] << " + ";
    int tmpint = progSpeciesNameVec.size()-1;
    cout << progSpeciesWeightVec[tmpint] << "*" << progSpeciesNameVec[tmpint] << endl;
    cout << "       fuel temperature (K):              " << Tfuel << endl;
    cout << "       oxidizer temperature (K):          " << Toxi << endl;
    cout << "       background pressure (Pa):          " << pressure << endl;
    cout << "       number of species:                 " << nspecies << endl;
    
  }

private:

  // =====================================
  // parse the header of the flamelet file
  // =====================================
  void parseHeaderAndInitParams() {
    
    int Tfuelread = 0;
    for (int i=0; i<headerVec.size(); ++i){
      
      istringstream buf(headerVec[i]);    
      string header;
      buf >> header;
      string dummy;
      
      if ( !header.compare("chi_st") )
	buf >> dummy >> chi_st;

      else if ( !header.compare("chi_ref") )
	buf >> dummy >> chi_st;

      
      else if ( !header.compare("Tmax") )
	buf >> dummy >> Tmax;
      
      else if ( !header.compare("pressure") ) {
	buf >> dummy >> pressure;
	pressure *= 100000.0; // convert from bar to Pa
      }
      
      // first comes fuel temperature
      else if ( !header.compare("Temperature") && Tfuelread!=1 ) { 
	buf >> dummy >> Tfuel;
	Tfuelread = 1;
      }
      
      // then comes oxidizer temperature
      else if ( !header.compare("Temperature") && Tfuelread!=0) 
	buf >> dummy >> Toxi;
      
      else if ( !header.compare("numOfSpecies") ) {
	int inttmp = 0;
	buf >> dummy >> inttmp;
	assert(nspecies == inttmp);
      }
      
      else if ( !header.compare("gridPoints") ) {
	int inttmp = 0;
	buf >> dummy >> inttmp;
	assert(npoints == inttmp);
      }      
      
      else if ( !header.compare("body") )
	CERR_S("shoould not be here")
      
    } //  for (int i=0; i<headerVec.size(); ++i){
    
  }
  
};

#endif
