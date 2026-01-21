// Test inflow turbulence program
// ==============================
// ==============================

#ifndef CREATECHEMTABLE_HPP
#define CREATECHEMTABLE_HPP

#include "Common.hpp"
#include "Flamelet.hpp"
#include "ArrayNd.hpp"
#include "candle.hpp"  // ars

extern void stretchgrid(const double* xp, const double* qp, const int np, double* &x, const int nx) {

  // build a 1D non-uniform grid from control points
  assert( qp[0]    == 0.0 );
  assert( qp[np-1] == 1.0 );

  const double ks = 0.2; // fraction of interval to apply smoothing (ks < 0.5)
  const double dq = 1.0/double(nx-1); // uniform spacing for q-grid

  double x0, x1;
  double q0, q1;

  double q  = 0; 
  int    ix = 0; 

  // build x from q=0 to first transition
  q0 = qp[1]-ks*(qp[1]-qp[0]);
  for (int n=0; n<nx; ++n) {
    x[ix] = xp[0]+(xp[1]-xp[0])*(q-qp[0])/(qp[1]-qp[0]);
    q+=dq; ++ix; if (q > q0) break;
  }

  // loop over transitions ...
  // build x in transition region and the neighboring interval
  for (int i=1; i<np-1; ++i) {

    x0 = xp[i]-ks*(xp[i]-xp[i-1]);
    q0 = qp[i]-ks*(qp[i]-qp[i-1]);
    x1 = xp[i]+ks*(xp[i+1]-xp[i]);
    q1 = qp[i]+ks*(qp[i+1]-qp[i]);

    // transition between intervals
    double w, xs, xl; 

    for (int n=0; n<nx; ++n) {
      if (q < qp[i]) { 
        xs = xp[i-1]+(xp[i]-xp[i-1])*(q-qp[i-1])/(qp[i]-qp[i-1]);
        w  = 0.5*(q-q0)/(qp[i]-q0);
      }
      else { 
        xs = xp[i]+(xp[i+1]-xp[i])*(q-qp[i])/(qp[i+1]-qp[i]);
        w  = 0.5*(q-q1)/(qp[i]-q1);
      }
      xl = x0+(x1-x0)*(q-q0)/(q1-q0);
      x[ix] = (1.0-w)*xs + w*xl;
      q+=dq; ++ix; if (q > q1) break;
    }

    // interior of next interval
    if (i == np-2)
      q0 = qp[np-1]+1.0e-15;
    else
      q0 = qp[i+1]-ks*(qp[i+1]-qp[i]);

    for (int n=0; n<nx; ++n) {
      x[ix] = xp[i]+(xp[i+1]-xp[i])*(q-qp[i])/(qp[i+1]-qp[i]);
      q+=dq; ++ix; if (q > q0) break;
    }
  }
  x[0]=xp[0]; x[nx-1]=xp[np-1]; // enforce endpoints
}


// ====================
// base chemtable class
// ====================
class AbstractChemtable {

public:

  string tablename;
  string tabletype;
  string filename;

public:

  AbstractChemtable() {
    cout << "AbstractChemtable()" << endl;
    // table name
    tablename = getStringParam("CHEMTABLE.NAME", "table");
  }
   
  virtual ~AbstractChemtable(){
    cout << "~AbstractChemtable()" << endl;
  }
  
  virtual void init() = 0;
  
  virtual void build() = 0;
  
  virtual void finalize() = 0;

  // ars
  void createPremixedFlamelets() {
    string fuel;//     = getStringParam("CHEMTABLE.FUEL");
    string oxidizer;// = getStringParam("CHEMTABLE.OXIDIZER");
    double pTab     = getDoubleParam("CHEMTABLE.P");
    string mech     = getStringParam("CHEMTABLE.MECHANISM", "gri30.xml");
    string mix      = getStringParam("CHEMTABLE.MIXING", "MIXTURE_AVERAGED");
    int nfl         =    getIntParam("CHEMTABLE.N_FLAMELETS", 100);
    double phi1     = getDoubleParam("CHEMTABLE.PHI_1", 0.5);
    double phi2     = getDoubleParam("CHEMTABLE.PHI_2", 1.5);
    int loglvl      =    getIntParam("CHEMTABLE.LOGLEVEL", 0);
    int soret       =    getIntParam("CHEMTABLE.SORET", 0);
    double Lewis    = getDoubleParam("CHEMTABLE.LEWIS_NUMBER", 1.0);
    int createInitFlamelet = getIntParam("CHEMTABLE.FLAMELET_GENERATE_INIT", 1);
    string init_file       = getStringParam("CHEMTABLE.FLAMELET_INIT_FILE", "init.xml");
    int onlyPhi1ToPh2      = getIntParam("CHEMTABLE.FLAMELET_GENERATE_PHI1_PHI2", 0);
    double Tf;
    double Tox;

    if (Param *param = getParam("CHEMTABLE.FUEL")) {
      fuel = param->getString(0);

      for (int ip=1; ip<param->size(); ip++) {
        if ( param->getString(ip).substr(0, 2) != "T:")
          fuel += " " + param->getString(ip);
        else {
          Tf = atof( param->getString(ip).substr(2).c_str() );
        }
      }
    }

    if (Param *param = getParam("CHEMTABLE.OXIDIZER")) {
      oxidizer = param->getString(0);

      for (int ip=1; ip<param->size(); ip++) {
        if ( param->getString(ip).substr(0, 2) != "T:")
          oxidizer += " " + param->getString(ip);
        else {
          Tox = atof( param->getString(ip).substr(2).c_str() );
        }
      }
    }

    cout << "\n ===== Premixed flamelet solver parameters ===== " << endl;
    cout << "\tFUEL \t\t= " << fuel << endl;
    cout << "\tOXIDIZER \t= " << oxidizer << endl;
    cout << "\tP \t\t= " << pTab << endl;
    cout << "\tMECHANISM \t= " << mech << endl;
    cout << "\tMIXING \t\t= " << mix << endl;
    cout << "\tN \t\t= " << nfl << endl;

//    Flame f("CH4", "Air", 1.05e5, 298.0, "gri30.xml", "Mix", 0);
//    Premixed f(fuel, oxid,   p,    tf, to, mech, tran, Lewis, loglev, createInitFlamelet, "init_p23_0", writeTecplot);
    Premixed f(fuel, oxidizer, pTab, Tf, Tox, mech, mix, Lewis, loglvl, createInitFlamelet, init_file, onlyPhi1ToPh2);//, soret);

    f.createFlamelets(phi1, phi2, nfl);
  }

  void createDiffusionFlamelets() {
    string fuel;//     = getStringParam("CHEMTABLE.FUEL");
    string oxidizer;// = getStringParam("CHEMTABLE.OXIDIZER");
    double pTab     = getDoubleParam("CHEMTABLE.P");
    string mech     = getStringParam("CHEMTABLE.MECHANISM", "gri30.xml");
    string mix      = getStringParam("CHEMTABLE.MIXING", "MIXTURE_AVERAGED");
    int nfl         =    getIntParam("CHEMTABLE.N_FLAMELETS", 100);
    int loglvl      =    getIntParam("CHEMTABLE.LOGLEVEL", 0);
//    int soret       =    getIntParam("CHEMTABLE.SORET", 0);
    double Lewis    = getDoubleParam("CHEMTABLE.LEWIS_NUMBER", 1.0);
    int ng          =    getIntParam("CHEMTABLE.NP_INIT", 51);
    double aInit    = getDoubleParam("CHEMTABLE.A_INIT", 0.5);
    double Tf;
    double Tox;

    if (Param *param = getParam("CHEMTABLE.FUEL")) {
      fuel = param->getString(0);

      for (int ip=1; ip<param->size(); ip++) {
        if ( param->getString(ip).substr(0, 2) != "T:")
          fuel += " " + param->getString(ip);
        else {
          Tf = atof( param->getString(ip).substr(2).c_str() );
        }
      }
    }

    if (Param *param = getParam("CHEMTABLE.OXIDIZER")) {
      oxidizer = param->getString(0);

      for (int ip=1; ip<param->size(); ip++) {
        if ( param->getString(ip).substr(0, 2) != "T:")
          oxidizer += " " + param->getString(ip);
        else {
          Tox = atof( param->getString(ip).substr(2).c_str() );
        }
      }
    }

    cout << "\n ===== Diffusion flamelet solver parameters ===== " << endl;
    cout << "\tFUEL \t\t= " << fuel << endl;
    cout << "\tOXIDIZER \t= " << oxidizer << endl;
    cout << "\tP \t\t= " << pTab << endl;
    cout << "\tMECHANISM \t= " << mech << endl;
    cout << "\tMIXING \t\t= " << mix << endl;
//    cout << "\tN \t\t= " << nfl << endl;

    Diffusion f(fuel, oxidizer, pTab, Tf, Tox, mech, mix, Lewis, loglvl, 1, ng, "init.xml", aInit);
//    Diffusion f(fuel, oxidizer, pTab, Tf, Tox, mech, mix, Lewis, loglvl, 0, "init2.xml", 1841.71866);
    f.createFlamelets(nfl);
  }


};

// ==================================
// Premixed flamelet table for the 
// variable density flow solver: VIDA
// ==================================
class VidaChemtablePFPVCart1D: public AbstractChemtable {

private:

  PremixedFlamelet* flamelet;
  vector<string> tableVarsNameVec;
  int nvars;
  double* C;
  int nC;
  double dC, Cmin, Cmax;
  double* weight;
  int (*index)[2];


public:

  VidaChemtablePFPVCart1D(const string& tabletype) {
    cout << "VidaChemtablePFPVCart1D()" << endl;
    this->tabletype = tabletype;
    filename = tablename + ".VidaPFPVCart1D.chem"; // tablename is initialized at the base class level
    flamelet = NULL;
    nvars    = 0;
    C        = NULL;
    nC       = 0;
    weight   = NULL;
    index    = NULL;
  }

  virtual ~VidaChemtablePFPVCart1D() {
    cout << "~VidaChemtablePFPVCart1D()" << endl;
    if ( flamelet != NULL) 
      delete flamelet;
    if ( C != NULL ) 
      delete[] C;
    if ( weight != NULL ) 
      delete[] weight;
    if ( index != NULL ) 
      delete[] index;

  }


  void init() {

    // read the list and initialize flamelet
    readFlameletFileList();
              
    assert(tableVarsNameVec.size() == 0);
    // a vida premixed table requires at least 8 variables
    //  (1) density,      
    //  (2) temperature,
    //  (3) viscosity, 
    //  (4) diffusivity,  
    //  (5) progress variable, 
    //  (6) progress variable source term,
    //  (7) axial coordinate
    //  (8) mass flow rate
    tableVarsNameVec.push_back("rho");
    tableVarsNameVec.push_back("T");
    tableVarsNameVec.push_back("mu");
    tableVarsNameVec.push_back("locp");
    tableVarsNameVec.push_back("prog");
    tableVarsNameVec.push_back("src_prog");
    tableVarsNameVec.push_back("y");
    tableVarsNameVec.push_back("mdot");
    if (checkParam("NOX_MODEL")) {
      tableVarsNameVec.push_back("src_NOX");
      flamelet->addNOxSrc();
    }
    // and add the user defined species
    //if (checkParam("CHEMTABLE.SPECIES")) 
    if (Param * param = getParam("CHEMTABLE.SPECIES"))
      for (int i = 0; i < param->size(); ++i) 
        tableVarsNameVec.push_back("Y_" + param->getString(i));
    // print on screen
    cout << " > species to be included in the table: " << endl;
    for (int i = 0; i<tableVarsNameVec.size(); ++i) 
      cout << "       " << tableVarsNameVec[i] << endl;

    nvars = tableVarsNameVec.size();
    assert(nvars>=8);
  }
  
  // =========================
  // build the data for table
  // =========================
  void build() {
    buildCgrid();
    findIndexAndWeight();
  }
  
  
  // ===================================
  // write the table into a binary file
  // and dump some info on screen
  // ===================================
  void finalize(){
    cout << " finalize()" << endl;    

    // write the table to a binary file
    writeBinary(filename);
    
    if ( checkParam("CHEMTABLE.TECPLOT") )
      writeTecplotFile();
  }

private:


  // ==============================
  // read flamelet file list
  // ==============================
  void readFlameletFileList() {

    string filename = getStringParam("CHEMTABLE.FLAMELET_LIST");        
    cout << " > reading the flamelet list file: " << filename << endl;

    // open the file to get the size
    ifstream ifp;
    // open the file
    ifp.open(filename.c_str());
    if ( ifp.fail() ) 
      CERR_S("could not open flamelet file list: " << filename );
    string buffer_str;
    int nfl = 0;
    while(getline(ifp, buffer_str)) {
      if (buffer_str.empty())
	break;
      ++nfl;
    }
    ifp.close();
    cout << "  found " << nfl << " flamelet files" << endl;

    if ( nfl != 1 ) 
      CERR_S("for VidaPFPVCart1D, you should only provide one flamelet file");

    cout << " > reading the flamelet files ..." << endl;
 
    // open the file again
    ifp.open(filename.c_str());
    if ( ifp.fail() ) 
      CERR_S("could not open flamelet file list: " << filename );
    int ifl = 0;
    int tmpint = max((int)(0.1*nfl),1);
    while(getline(ifp, buffer_str)) {
      if (buffer_str.empty())
	break;
      istringstream buf(buffer_str);
      string thisfile;
      buf >> thisfile;
      // read the flamelet file
      cout << " > reading the flamelet file: " << thisfile << endl;
      flamelet = new PremixedFlamelet();
      flamelet->init(thisfile);
      ++ifl;
    }
    ifp.close();
    
    // check
    assert(ifl == nfl);
  }
  

  
  // =================================
  // write the table to a tecplot file
  // =================================
  void writeTecplotFile() {

    string filename = tablename + ".dat";

    cout << " > writing table to the tecplot file: " << filename << endl;
   
    double* prog = flamelet->getVarPtr("prog");

    vector<double*> tmpVec;
    for (int ivar=0; ivar<tableVarsNameVec.size(); ++ivar) {
      double* tmp = flamelet->getVarPtr(tableVarsNameVec[ivar]);
      tmpVec.push_back(tmp);
    }
           
    FILE * fp = fopen(filename.c_str(),"w");
    fprintf(fp,"title = \"chemtable\"\n");
    fprintf(fp,"variables = \"C\"\n");
    for (int ivar=0; ivar<tableVarsNameVec.size(); ++ivar)
      fprintf(fp,"\"%s\"\n",tableVarsNameVec[ivar].c_str());
    fprintf(fp,"zone i = %d, DATAPACKING=POINT\n",nC);
    for (int i=0; i<nC; ++i) {
      fprintf(fp,"%g",C[i]);
      for (int ivar=0; ivar<tableVarsNameVec.size(); ++ivar) {
	double* tmp = tmpVec[ivar];
	double thistmp = weight[i]*tmp[index[i][0]] + (1.0-weight[i])*tmp[index[i][1]];
	fprintf(fp," %g",thistmp);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
    
  }

  // ============================
  // build progress variable grid
  // ============================

  void buildCgrid() {
    nC    = getIntParam("CART_CHEMTABLE.nC");
    assert(nC>=2);
    assert( C == NULL);
    C     = new double[nC];

    Cmin = flamelet->getVarMin("prog");
    Cmax = flamelet->getVarMax("prog");
    assert(Cmax>Cmin);

    dC = (Cmax - Cmin) / double(nC-1);
    for (int i=0; i<nC; ++i )
      C[i] = Cmin + double(i) * dC;
     
    assert(fabs(C[0]-Cmin)<1.0e-14);
    assert(fabs(C[nC-1]-Cmax)<1.0e-14);
  }


  // =============================
  // build the interpolation stuff
  // =============================

  void findIndexAndWeight(){

    assert(index==NULL);
    index = new int[nC][2];
    assert(weight==NULL);
    weight = new double[nC];

    double* prog = flamelet->getVarPtr("prog");
    

    // not very good but live with it for now
    double Ctarget;   
    for (int j=0; j<nC; ++j) {
      Ctarget = C[j];
      Ctarget = min(Ctarget,Cmax-1.0e-14);
      Ctarget = max(Ctarget,Cmin);
      int ub = 0;
      while ( (Ctarget>=prog[ub]) && (ub<(flamelet->npoints-1)) )
	++ub;
      int lb = ub-1;      
      assert(lb>=0);
      assert(ub<=(flamelet->npoints-1));
      assert(prog[lb]<prog[ub]);
      double L = prog[ub] - prog[lb];
      assert(L>0.0);
      index[j][0] = lb;
      index[j][1] = ub;
      weight[j]   = (prog[ub] - Ctarget) / L;
      assert(weight[j]<=1.0);
      assert(weight[j]>=0.0);	     
    }
  }
  
  
  // ===============================
  // write a core IO compatible file
  // ===============================
  void writeBinary(const string& filename){
    
    cout << " > writing vida compatible table to the binary file: " << filename << endl;

    // get the progress variable
    double* prog = flamelet->getVarPtr("prog");
    
    // open the file
    FILE * fp = NULL;
    fp=fopen(filename.c_str(),"wb");
    if (fp == NULL) 
      CERR_S(" cannot open chemtable file: " << filename);
      
    assert(fp);
    

    // use the header in the core IO ...
   

    // (1) write ugp magic number and version
    {
      int tmpint[2]={CHEM_IO_MAGIC_NUMBER,CHEM_IO_VERSION};
      size_t dummy= fwrite(tmpint,sizeof(int),2,fp);
      assert(dummy==2);
    }
    
    // (2) write the table type and dimensions
    {
      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      size_t dummy = tabletype.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and bye skip
      header.id = UGP_IO_CT_CART_1D;
      header.skip = sizeof(Header);
      // table dimensions
      header.idata[0] = nC; 
      header.idata[1] = nvars;
      // table pressure and temperature
      header.rdata[0] = flamelet->pressure;
      header.rdata[1] = flamelet->laminarFlameSpeed;
      header.rdata[2] = flamelet->laminarFlameThickness; 
      // write into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
    }

    // (3) write "C" coordinate
    {
      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      string tmp    = "C";
      size_t dummy = tmp.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and bye skip
      header.id = UGP_IO_CT_COOR;
      header.skip = nC * sizeof(double) + sizeof(Header);
      // table dimensions
      header.idata[0] = nC; 
      // C coordinate
      header.rdata[0] = dC;
      header.rdata[1] = Cmin;
      header.rdata[2] = Cmax;
      // write the into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
      // and the coordinate
      dummy= fwrite(C,sizeof(double),nC,fp);
      assert(dummy==nC);
    }

    // (4) write the data
    for(int ivar=0; ivar<nvars; ++ivar){
      
      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      size_t dummy = tableVarsNameVec[ivar].copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and bye skip
      header.id = UGP_IO_CT_DATA;

      double* tmpVec = flamelet->getVarPtr(tableVarsNameVec[ivar]);
      header.skip = nC * sizeof(double) + sizeof(Header);
      // table dimensions
      header.idata[0] = nC; 
      // write the header into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
      // and the data
      for (int i=0; i<nC; ++i ) {
	double tmp = weight[i]*tmpVec[index[i][0]] + (1.0-weight[i])*tmpVec[index[i][1]];
	dummy= fwrite(&tmp,sizeof(double),1,fp);
	assert(dummy==1);
      }
    }
    
    // end of file header
    {
      Header header;
      header.id = UGP_IO_EOF;
      sprintf(header.name,"EOF");
      header.skip = sizeof(Header);
      // write the header into the file      
      size_t dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);  
    }
    
    // close the file
    fclose(fp);
  
  }

};


// =====================================
// 2D (Z,C) premixed flamelet table
// for Vida w/ thickened flame model
// =====================================
class VidaChemtablePFPVCart2D: public AbstractChemtable {

protected:

  list<PremixedFlamelet*> flameletList;
  vector<string> tableVarsNameVec;

  int nz; 
  int nc; 
  int nvars;
  double pressure;
  double * mixfrac;
  double * Z;
  double * C;
  CTIarray<double> * table;

public:

  VidaChemtablePFPVCart2D(const string& tabletype) {

    cout << "VidaChemtablePFPVCart2D()" << endl;
    this->tabletype = tabletype;
    filename = tablename + ".VidaPFPVCart2D.chem"; // tablename is initialized at the base class level

    nz       = 0;
    nc       = 0;
    nvars    = 0;
    pressure = 0.0; 
    mixfrac  = NULL;
    Z        = NULL;
    C        = NULL;
    table    = NULL;

  }

  virtual ~VidaChemtablePFPVCart2D() {

    cout << "~VidaChemtablePFPVCart2D()" << endl;

    if ( mixfrac != NULL ) delete[] mixfrac;
    if ( Z       != NULL ) delete[] Z;
    if ( C       != NULL ) delete[] C;
    if ( table   != NULL ) delete   table; 

    for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f)
      if ((*f) != NULL ) delete (*f);

  }

  void init() {

    cout << " init()" << endl; 

    // =======================================
    // generate flamelets
    // =======================================
    if (!checkParam("CHEMTABLE.FLAMELET_LIST"))
      createPremixedFlamelets();

    //========================================
    // get grid info
    //========================================

    nz = getIntParam("CART_CHEMTABLE.nZ"); 
    Z  = new double[nz];

    nc = getIntParam("CART_CHEMTABLE.nC");  
    C  = new double[nc];

    // read files and sort by phi
    readAndSortFlamelets();

    //========================================
    // initialize variables
    //========================================

    assert(tableVarsNameVec.size() == 0);
    // a 2D vida premixed table requires at least 8 variables
    //  (1) density,      
    //  (2) temperature,
    //  (3) viscosity, 
    //  (4) diffusivity,  
    //  (5) progress variable, 
    //  (6) progress variable source term,
    //  (7) laminar flame speed
    //  (8) laminar flame thickness
    tableVarsNameVec.push_back("rho");
    tableVarsNameVec.push_back("T");
    tableVarsNameVec.push_back("mu");
    tableVarsNameVec.push_back("locp");
    tableVarsNameVec.push_back("prog");
    tableVarsNameVec.push_back("src_prog");
    tableVarsNameVec.push_back("sL");
    tableVarsNameVec.push_back("lF");
    tableVarsNameVec.push_back("mw");
    for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f) {
      (*f)->addLaminarFlameThickness();
      (*f)->addLaminarFlameSpeed();
    }
    if (checkParam("NOX_MODEL")) {
      tableVarsNameVec.push_back("src_NOX");
      for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f)
        (*f)->addNOxSrc();
    }

    // and add the user defined species
    // options:
    // Y_SPECIES    : mass fraction of SPECIES
    // X_SPECIES    : mole fraction of SPECIES
    // X_SPECIES_DRY: dry mole fraction of SPECIES
    // SPECIES      : mass, mole, and dry mole fractions of SPECIES
    //if (checkParam("CHEMTABLE.SPECIES")) 
    if (Param * param = getParam("CHEMTABLE.SPECIES")) {
      for (int i = 0; i < param->size(); ++i) {
        string var_name = param->getString(i);

        if (var_name[0] == 'X' || var_name[0] == 'Y') {
          tableVarsNameVec.push_back(var_name);

          if (var_name[0] == 'X') {
//            cout << " >>>>>>> species name: " << var_name.substr(2,var_name.length()-1) << endl;
            for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f) {
              if ( var_name.substr(var_name.length()-4) == "_DRY" || var_name.substr(var_name.length()-4) == "_dry" )
                (*f)->addDryMoleFraction(var_name.substr(2,var_name.length()-6));
              else
                (*f)->addMoleFraction(var_name.substr(2,var_name.length()-2));
            }
          }
        }//(var_name[0] == "X" || "Y")
        else {
          tableVarsNameVec.push_back("Y_"+var_name);         // add mass fraction
          tableVarsNameVec.push_back("X_"+var_name);         // add mole fraction
          if (var_name != "H2O") tableVarsNameVec.push_back("X_"+var_name+"_DRY");  // add dry mole fraction

          for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f) {
            (*f)->addMoleFraction(var_name);
            if (var_name != "H2O") (*f)->addDryMoleFraction(var_name);
          }
        }//(var_name[0] != "X" || "Y")
      }//(i < param->size())
    }//CHEMTABLE.SPECIES

    // print on screen
    cout << " > variables to be included in the table: " << endl;
    for (int i = 0; i<tableVarsNameVec.size(); ++i) 
      cout << "       " << tableVarsNameVec[i] << endl;

    nvars = tableVarsNameVec.size();
    assert(nvars>=8);

    // report min/max for each variable
    for (int i=0; i<nvars; i++){ 
      double vmin =  1.0e+20; 
      double vmax = -1.0e+20; 
      for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f){
        vmin=min(vmin,(*f)->getVarMin(tableVarsNameVec[i])); 
        vmax=max(vmax,(*f)->getVarMax(tableVarsNameVec[i])); 
      }
      cout << " > " << tableVarsNameVec[i] << ": min = " << vmin << ", max = " << vmax << endl;
    }

  }
  
  // =========================
  // build the data for table
  // =========================
  void build() {

    cout << " build()" << endl; 
    buildZgrid();
    buildCgrid();
    interpData(); 
  }
  
  // ===================================
  // write the table into a binary file
  // and dump some info on screen
  // ===================================
  void finalize(){
    
    cout << " finalize()" << endl;

    // write the table to a binary file
    writeBinary(filename);
    
    if (checkParam("CHEMTABLE.TECPLOT"))
      writeTecplotFile();
    
  }

private:

  void readAndSortFlamelets() { 

    assert(flameletList.size()==0);

    string filename = getStringParam("CHEMTABLE.FLAMELET_LIST", "FlameletList.txt");
    cout << " > reading flamelet list file: " << filename << endl;

    // open the file to get the size
    ifstream ifp;
    // open the file
    ifp.open(filename.c_str());
    if ( ifp.fail() )
      CERR_S("could not open flamelet file list: " << filename );
    string buffer_str;
    int nfl = 0;
    while(getline(ifp, buffer_str)) {
      if (buffer_str.empty())
        break;
      ++nfl;
    }
    ifp.close();
    cout << "   found " << nfl << " flamelet files" << endl;

    if ( nfl < 2 )
      CERR_S("VidaPFPVCart2D requires more than one flamelet file");

    cout << " > reading flamelet files ..." << endl;

    // open the file again
    ifp.open(filename.c_str());
    if ( ifp.fail() )
      CERR_S("could not open flamelet file list: " << filename );

    cout << "   progress:"; cout.flush(); 
    int ifl = 0;
    int tmpint = max((int)(0.1*nfl),1);
    while(getline(ifp, buffer_str)) {
      if ( (ifl+1)%tmpint == 0 ) {
        cout << " .." << (ifl+1)*10/tmpint << "%";
        cout.flush(); 
      }
      if (buffer_str.empty())
        break;
      istringstream buf(buffer_str);
      string thisfile;
      buf >> thisfile;

      // read the flamelet file
      PremixedFlamelet * flamelet = new PremixedFlamelet();
      flamelet->init(thisfile);
      flameletList.push_back(flamelet); 
      ++ifl;
    }
    ifp.close();
    cout << endl;
    assert(ifl == nfl);

    // sort based on phi
    flameletList.sort(isPhiSmallerFlamelet<PremixedFlamelet>); 

    // dump flamelet list if requested
    if ( checkParam("CHEMTABLE.DUMP_FLAMELETLIST") ) {
      cout << " > flamelet list ..." << endl;
      for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f) 
        cout << (*f)->filename << " ";
      cout << endl;
    }

    // dump some info on screen
    if (cti_verbose) {
      list<PremixedFlamelet*>::iterator f=flameletList.begin();
      (*f)->info();
    }
    
    {
      // check that the pressure is the same for all flamelets and store it
      list<PremixedFlamelet*>::iterator f=flameletList.begin();
      pressure = (*f)->pressure;
      for (f=++(flameletList.begin()); f!=flameletList.end(); ++f) {
        assert(pressure == (*f)->pressure );
      }
    }

  }

  void buildZgrid() { 

    cout << " > building Z grid ..." << endl;

    //--------------------------------------------
    // construct mixture fraction from flamelets
    //--------------------------------------------

    int nf = flameletList.size(); 
    mixfrac = new double[nf];

    int iz = 0;
    for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f) {

      // get unburned mixture fuel mass fraction
      // XXXX WARNING: ignoring O2,N2 is not robust for all mechanisms (i.e. could have different names)
      //               there could also be trace O2 or N2 in the fuel (i.e. real natural gas) 

      double Yfuel=0.0; 
      for (int i=0; i<(*f)->unburntSpeciesVec.size(); ++i)
        if ( (*f)->unburntSpeciesVec[i] != "O2"
          && (*f)->unburntSpeciesVec[i] != "N2" )
          Yfuel += (*f)->unburntMassfractionVec[i];
      mixfrac[iz] = Yfuel; 
      ++iz; 
      //cout << "   Z[" << iz-1 << "] = " << mixfrac[iz-1] << endl; 
    } 
    assert(iz==nf);

    //--------------------------------------------
    // build mixture fraction grid for table
    //--------------------------------------------

    double zmin=mixfrac[0];
    double zmax=mixfrac[flameletList.size()-1];

    //if (checkParam("CART_CHEMTABLE.STRETCH_Z")) 
    if (Param * param = getParam("CART_CHEMTABLE.STRETCH_Z")) {
      int np=4;
      double *xp = new double[np];
      xp[0] = zmin;
      xp[1] = max(param->getDouble(0), zmin + 1.0E-8);
      xp[2] = min(param->getDouble(1), zmax-1.0E-8);
      xp[3] = zmax;
      double frac  = param->getDouble(2);

      assert( (xp[0] < xp[1]) &&
              (xp[1] < xp[2]) &&
              (xp[2] < xp[3]) &&
              (frac  >  0.0 ) &&
              (frac  <  1.0 ) );

      double *qp = new double[np];
      qp[0] = 0.0;
      qp[1] = (1.0-frac)*(xp[1]-xp[0])/(xp[1]-xp[0] + xp[3]-xp[2]);
      qp[2] = qp[1]+frac;
      qp[3] = 1.0;

      stretchgrid(xp,qp,np,Z,nz);
      delete [] xp;
      delete [] qp;
    }
    else // uniform grid
      for (int i=0; i<nz; i++)
        Z[i] = zmin + (zmax-zmin)*double(i)/double(nz-1);

    //for (int i=0; i<nz; i++) cout << Z[i] << endl; 
  }

  void buildCgrid() {
  
    cout << " > building C grid ..." << endl;

    //--------------------------------------------
    // build progress variable grid for table
    //--------------------------------------------

    double cmin =  1.0e+20; 
    double cmax = -1.0e+20; 
    for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f) {
      cmin = min(cmin,(*f)->Cmin); 
      cmax = max(cmax,(*f)->Cmax); 
    }
  
    //if (checkParam("CART_CHEMTABLE.STRETCH_C")) 
    if (Param * param = getParam("CART_CHEMTABLE.STRETCH_C")) {
      int np=4;
      double *xp = new double[np];
      xp[0] = cmin;
      xp[1] = max(param->getDouble(0), cmin + 1.0E-8);
      xp[2] = min(param->getDouble(1), cmax - 1.0E-8);
      xp[3] = cmax;
      double frac  = param->getDouble(2);

      assert( (xp[0] < xp[1]) &&
              (xp[1] < xp[2]) &&
              (xp[2] < xp[3]) &&
              (frac  >  0.0 ) &&
              (frac  <  1.0 ) );

      double *qp = new double[np];
      qp[0] = 0.0;
      qp[1] = (1.0-frac)*(xp[1]-xp[0])/(xp[1]-xp[0] + xp[3]-xp[2]);
      qp[2] = qp[1]+frac;
      qp[3] = 1.0;

      stretchgrid(xp,qp,np,C,nc);
      delete [] xp;
      delete [] qp;
    }
    else // uniform grid
      for (int i=0; i<nc; i++)
        C[i] = cmin + (cmax-cmin)*double(i)/double(nc-1);

    //for (int i=0; i<nc; i++) cout << C[i] << endl; 
  }

  void interpData() { 

    cout << " > interpolating data ..." << endl;

    assert(table == NULL);
    table = new CTIarray<double>(nz,nc,nvars);
    int nf = flameletList.size(); 

    int iprog=-1; 
    int irho =-1; 
    int isrc =-1; 
    int isL  =-1;
    int ilF  =-1;
    for (int ivar=0; ivar<nvars; ++ivar)
      if (tableVarsNameVec[ivar].compare("prog") == 0) iprog=ivar;
      else if (tableVarsNameVec[ivar].compare("rho") == 0) irho=ivar;
      else if (tableVarsNameVec[ivar].compare("src_prog") == 0) isrc=ivar;
      else if (tableVarsNameVec[ivar].compare("sL") == 0) isL=ivar;
      else if (tableVarsNameVec[ivar].compare("lF") == 0) ilF=ivar;
    assert(iprog != -1);
    assert(irho  != -1);
    assert(isrc  != -1);
    assert(isL   != -1);
    assert(ilF   != -1);

    for (int i=0; i<nz; ++i) {

      // get flamelets that bracket current z value
      int iz; 
      for (iz=0; iz<nf-1; ++iz) 
        if (mixfrac[iz+1] > Z[i]) break;
      iz = min(iz,nf-2);
      double wz = max(0.0,min(1.0,(mixfrac[iz+1]-Z[i])/(mixfrac[iz+1]-mixfrac[iz]))); 

      list<PremixedFlamelet*>::iterator f0 = flameletList.begin(); advance(f0,iz  );
      list<PremixedFlamelet*>::iterator f1 = flameletList.begin(); advance(f1,iz+1);
      int n0 = (*f0)->npoints;
      int n1 = (*f1)->npoints;

      // oversample flamelets (if not the same size) and interpolate to z grid
      int np=max(n0,n1);
      double data[np][nvars];

      if (n0 < np) { 
        for (int j=0; j<np; ++j) { 
          double x = double(j)/double(np-1); 
          int   ii = min(n0-2,int(x*double(n0-1)));
          double w = max(0.0,min(1.0,double(ii+1)-x*double(n0-1))); 
           
          for (int ivar=0; ivar<nvars; ++ivar) { 
            if (ivar == irho ) // use inverse
              data[j][ivar] = (     w /(*f0)->getVarPtr(tableVarsNameVec[ivar])[ii  ]
                              +(1.0-w)/(*f0)->getVarPtr(tableVarsNameVec[ivar])[ii+1])*wz; 
            else if (ivar == isrc || ivar == isL || ivar == ilF) // use logarithm
              data[j][ivar] = (log(max(1.0e-16,(*f0)->getVarPtr(tableVarsNameVec[ivar])[ii  ]))*w
                              +log(max(1.0e-16,(*f0)->getVarPtr(tableVarsNameVec[ivar])[ii+1]))*(1.0-w))*wz; 
            else
              data[j][ivar] = ((*f0)->getVarPtr(tableVarsNameVec[ivar])[ii  ]*w
                              +(*f0)->getVarPtr(tableVarsNameVec[ivar])[ii+1]*(1.0-w))*wz; 
          }
        }
      }
      else {
        for (int j=0; j<np; ++j) { 
          for (int ivar=0; ivar<nvars; ++ivar) {
            if (ivar == irho ) // use inverse
              data[j][ivar] = wz/(*f0)->getVarPtr(tableVarsNameVec[ivar])[j];
            else if (ivar == isrc || ivar == isL || ivar == ilF) // use logarithm
              data[j][ivar] = log(max(1.0e-16,(*f0)->getVarPtr(tableVarsNameVec[ivar])[j]))*wz;
            else
              data[j][ivar] = (*f0)->getVarPtr(tableVarsNameVec[ivar])[j]*wz;
          }
        }
      }
      // NOTE: data now contains wz*f0

      if (n1 < np) {
        for (int j=0; j<np; ++j) { 
          double x = double(j)/double(np-1); 
          int   ii = min(n1-2,int(x*double(n1-1)));
          double w = max(0.0,min(1.0,double(ii+1)-x*double(n1-1))); 
           
          for (int ivar=0; ivar<nvars; ++ivar) {
            if (ivar == irho ) // use inverse
              data[j][ivar] += (     w /(*f1)->getVarPtr(tableVarsNameVec[ivar])[ii  ]
                               +(1.0-w)/(*f1)->getVarPtr(tableVarsNameVec[ivar])[ii+1])*(1.0-wz); 
            else if (ivar == isrc || ivar == isL || ivar == ilF) // use logarithm
              data[j][ivar] += (log(max(1.0e-16,(*f1)->getVarPtr(tableVarsNameVec[ivar])[ii  ]))*w
                               +log(max(1.0e-16,(*f1)->getVarPtr(tableVarsNameVec[ivar])[ii+1]))*(1.0-w))*(1.0-wz); 
            else
              data[j][ivar] += ((*f1)->getVarPtr(tableVarsNameVec[ivar])[ii  ]*w
                               +(*f1)->getVarPtr(tableVarsNameVec[ivar])[ii+1]*(1.0-w))*(1.0-wz); 
          }
        }
      }
      else {
        for (int j=0; j<np; ++j) { 
          for (int ivar=0; ivar<nvars; ++ivar)
            if (ivar == irho ) // use inverse
              data[j][ivar] += (1.0-wz)/(*f1)->getVarPtr(tableVarsNameVec[ivar])[j]; 
            else if (ivar == isrc || ivar == isL || ivar == ilF) // use logarithm
              data[j][ivar] += log(max(1.0e-16,(*f1)->getVarPtr(tableVarsNameVec[ivar])[j]))*(1.0-wz); 
            else
              data[j][ivar] += (*f1)->getVarPtr(tableVarsNameVec[ivar])[j]*(1.0-wz); 
        }
      }
      // NOTE: data now contains wz*f0+(1-wz)*f1, i.e. a flamelet at the correct z value

      // interpolate data to c grid
      for (int j=0; j<nc; ++j) { 
        int ic;
        for (ic=0; ic<np-1; ++ic)
          if (data[ic+1][iprog] > C[j]) break;
        ic = min(ic,np-2); 

        double wc = 0.0;
        double denom = data[ic+1][iprog]-data[ic][iprog]; 
        if (denom > 1.0e-10)
          wc = max(0.0,min(1.0,(data[ic+1][iprog]-C[j])/denom));

        for (int ivar=0; ivar<nvars; ++ivar) { 
          (*table)(i,j,ivar) = wc*data[ic][ivar] + (1.0-wc)*data[ic+1][ivar];

          if (ivar == irho)
            (*table)(i,j,ivar) = 1.0/(*table)(i,j,ivar); 
          else if (ivar == isrc || ivar == isL || ivar == ilF)
            (*table)(i,j,ivar) = exp((*table)(i,j,ivar)); 
        }

      }
      
    }
  }

  // =================================
  // write the table to a tecplot file
  // =================================
  void writeTecplotFile() {

    string filename = tablename + ".dat";

    cout << " > writing table to tecplot file: " << filename << endl;
   
    FILE * fp = fopen(filename.c_str(),"w");
    fprintf(fp,"title = \"chemtable\"\n");
    fprintf(fp,"variables = \"Z\"\n\"C\"\n");
    for (int ivar=0; ivar<tableVarsNameVec.size(); ++ivar)
      fprintf(fp,"\"%s\"\n",tableVarsNameVec[ivar].c_str());
    fprintf(fp,"zone i = %d, j = %d, DATAPACKING=POINT\n",nz,nc);
    for (int j=0; j<nc; ++j) {
      for (int i=0; i<nz; ++i) {
        fprintf(fp,"%g %g",Z[i],C[j]);
        for (int ivar=0; ivar<tableVarsNameVec.size(); ++ivar)
	  fprintf(fp," %g",(*table)(i,j,ivar));
        fprintf(fp,"\n");
      }
    }
    fclose(fp);
  }

  // ===============================
  // write a core IO compatible file
  // ===============================
  void writeBinary(const string& filename){
    
    cout << " > writing vida compatible table to the binary file: " << filename << endl;

    // open the file
    FILE * fp = NULL;
    fp=fopen(filename.c_str(),"wb");
    if (fp == NULL) 
      CERR_S(" cannot open chemtable file: " << filename);
    assert(fp);

    // use the header in the core IO ...

    // (1) write ugp magic number and version
    {
      int tmpint[2]={CHEM_IO_MAGIC_NUMBER,CHEM_IO_VERSION};
      size_t dummy= fwrite(tmpint,sizeof(int),2,fp);
      assert(dummy==2);
    }
    
    // (2) write the table type and dimensions
    {
      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      size_t dummy = tabletype.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and bye skip
      header.id = UGP_IO_CT_CART_2D;
      header.skip = sizeof(Header);
      // table dimensions
      header.idata[0] = nz; 
      header.idata[1] = nc; 
      header.idata[2] = nvars;
      // table pressure
      header.rdata[0] = pressure;
      // write into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
    }

    // (3) write "Z" coordinate
    {
      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
        header.name[i] = '\0';  
      string tmp    = "Z";
      size_t dummy = tmp.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and bye skip
      header.id = UGP_IO_CT_COOR;
      header.skip = nz * sizeof(double) + sizeof(Header);
      // table dimensions
      header.idata[0] = nz; 
      // write the into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
      // and the coordinate
      dummy= fwrite(Z,sizeof(double),nz,fp);
      assert(dummy==nz);
    }

    // (4) write "C" coordinate
    {
      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
        header.name[i] = '\0';  
      string tmp   = "C";
      size_t dummy = tmp.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and bye skip
      header.id = UGP_IO_CT_COOR;
      header.skip = nc * sizeof(double) + sizeof(Header);
      // table dimensions
      header.idata[0] = nc; 
      // write the into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
      // and the coordinate
      dummy= fwrite(C,sizeof(double),nc,fp);
      assert(dummy==nc);
    }

    // (5) write the data
    for(int ivar=0; ivar<nvars; ++ivar){

      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
        header.name[i] = '\0';  
      size_t dummy = tableVarsNameVec[ivar].copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and bye skip
      header.id = UGP_IO_CT_DATA;
      header.skip = nz*nc * sizeof(double) + sizeof(Header);
      // table dimensions
      header.idata[0] = nz*nc; 
      header.idata[1] = nz; 
      header.idata[2] = nc; 
      // write the header into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
      // and the data
      for (int iz=0; iz<nz; ++iz)
        for (int ic=0; ic<nc; ++ic) {
          double tmp = (*table)(iz,ic,ivar);
          dummy= fwrite(&tmp,sizeof(double),1,fp);
          assert(dummy==1);
        }
    }
        
    // end of file header
    {
      Header header;
      header.id = UGP_IO_EOF;
      sprintf(header.name,"EOF");
      header.skip = sizeof(Header);
      // write the header into the file      
      size_t dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);  
    }
    
    // close the file
    fclose(fp);
  
  }

}; //VidaChemtablePFPVCart2D


// ===============================
// Flamelet table for the variable
// density flow solver: VIDA
// ===============================
class ChemtableNFPVCart3D: public AbstractChemtable {

private:
 
  double pressure;
  double Tfuel;
  double Toxi;

  CTIarray<double> * table;

protected:

  list<NonPremixedFlamelet*> flameletList;
 
  vector<string> tableVarsNameVec;
  int nz;
  int nzvar;
  int nc;
  int nvars;

  double * Z;
  double * Zvar;
  double * C;

  CTIarray<double>* data4d;


public:

  ChemtableNFPVCart3D() {

    cout << "ChemtableNFPVCart3D()" << endl;
    cout << " building a 3D cartesian table using non-premixed fpv model ... " << endl;

    Z    = NULL;
    Zvar = NULL;
    C    = NULL;

    data4d  = NULL;
    table   = NULL;

  }

  virtual ~ChemtableNFPVCart3D() {
    
    cout << "~ChemtableNFPVCart3D()" << endl;

    if ( Z    != NULL ) delete[] Z;
    if ( Zvar != NULL ) delete[] Zvar;
    if ( C    != NULL ) delete[] C;
    
    for (list<NonPremixedFlamelet*>::iterator flamelet=flameletList.begin(); flamelet!=flameletList.end(); ++flamelet)
      if ( (*flamelet) != NULL ) delete (*flamelet);
    
    if ( data4d  != NULL ) delete data4d;
    if ( table   != NULL ) delete table;
    
  }
 

  // =====================================
  // variables to be included in the table
  // =====================================   
  virtual void initTableVarsNameVec() = 0;
  
  // =======================
  // compute the convolution
  // =======================
  virtual void calcConvolution() = 0;  

  // =============================
  // initialize table and its data
  // =============================
  void init() {

    // =======================================
    // generate flamelets
    // =======================================
    if (!checkParam("CHEMTABLE.FLAMELET_LIST"))
      createDiffusionFlamelets();

    // =====================================
    // variables to be included in the table
    // =====================================   
    initTableVarsNameVec();
 
     
    // ==============================
    // initialize grid for Z and Zvar
    // ==============================
    nz    = getIntParam("CART_CHEMTABLE.nZ");
    assert( Z == NULL);
    Z     = new double[nz]; 
    
    nzvar = getIntParam("CART_CHEMTABLE.nZvar");
    assert( Zvar == NULL);
    Zvar  = new double[nzvar]; 
    

    // =======================
    // and get dimension for C
    // will be build after we 
    // have its bounds
    // =======================
    nc    = getIntParam("CART_CHEMTABLE.nC");
    assert( C == NULL);
    C     = new double[nc];

    cout << " > table dimensions: (nvars,nz,nzvar,nc) = (" <<nvars<<","<<nz<<","<<nzvar<<"," <<nc<<")"<< endl;

  }


  // ===================================
  // read flamelet files and build table
  // ===================================
  void build(){

    cout << " build()" << endl;

    // read the flamelet files and sort based in maximum progress variable
    readFlameletFilesAndSort();

    // build grid for Z and Zvar
    // need Z and Zvar for cleanup
    buildZgrid();
    buildZvargrid();

    // cleanup the list based on unique mapping in C direction
    cleanupFlameletList();
    
    // build the c grid based on min and max progress variable
    buildCgrid(); 
    
    // convolve with the presumed pdf
    calcConvolution();

    //dumpPointCloud(); //lee
    
    // and build the c mapping
    buildCmapping();

  }

  
  // ===================================
  // write the table into a binary file
  // and dump some info on screen
  // ===================================
  void finalize(){
    
    cout << " finalize()" << endl;    

    // write the table to a binary file

    writeBinary(filename);

    if ( checkParam("CHEMTABLE.TECPLOT") )
      dumpTecplotFile();
    
    info();

  }
  
  // ===============================
  // dump table infomation on screen
  // ===============================
  void info() {

    cout << "\n====================== chemtable info ======================"<<endl;
    cout << "   type:                     " << tabletype << endl;
    cout << "   dimension (nz,nzvar,nc):  " << "(" << nz << "," << nzvar << "," << nc << ")" << endl;
    cout << "   reference pressure:       " << pressure << endl;
    cout << "   oxidizer temperature:     " << Toxi << endl;
    cout << "   fuel temperature:         " << Tfuel << endl;
    
    cout << "   coordinates:" << endl;
    double varmin = 1.0e20;
    double varmax = 1.0e-20;
    getVarMinAndMax(varmin,varmax,Z,nz);
    cout << "      ";
    cout.width(8);      
    cout << left << "Z";
    cout.width(8);
    cout << "   (min,max):  " << "(" << varmin << "," << varmax << ")" << endl;

    getVarMinAndMax(varmin,varmax,Zvar,nzvar);
    cout << "      ";
    cout.width(8);      
    cout << left << "Zvar";
    cout.width(8);
    cout << "   (min,max):  " << "(" << varmin << "," << varmax << ")" << endl;

    getVarMinAndMax(varmin,varmax,C,nc);
    cout << "      ";
    cout.width(8);      
    cout << left << "C";
    cout.width(8);
    cout << "   (min,max):  " << "(" << varmin << "," << varmax << ")" << endl;
    
    cout << "   variables (nvars = " << nvars << "):" << endl;
    for (int ivar=0; ivar<tableVarsNameVec.size(); ++ivar) {
      double varmin = 1.0e20;
      double varmax = 1.0e-20;
      getVarMinAndMax(varmin,varmax,&(*table)(ivar,0,0,0),nz*nzvar*nc);
      cout << "      ";
      cout.width(8);      
      cout << left << tableVarsNameVec[ivar];
      cout.width(8);
      cout << "   (min,max):  " << "(" << varmin << "," << varmax << ")" << endl;
    }

    cout << "============================================================\n"<<endl;

  }


private:

  
  // ========================================
  // report min and max of the table variable
  // ========================================
  void getVarMinAndMax(double& varmin,double& varmax,const double* data, const int n){
    
    varmin = data[0];
    varmax = data[0];
    for (int i=1; i<n; ++i){
      double tmp = data[i];	  
      varmin = min(varmin,tmp);
      varmax = max(varmax,tmp);
    }

  }


protected:

  // ===============================
  // write a core IO compatible file
  // ===============================
  void writeBinary(const string& filename){
    
    cout << " > writing table to the binary file: " << filename << endl;

    // open the file
    FILE * fp = NULL;
    fp=fopen(filename.c_str(),"wb");
    if (fp == NULL) 
      CERR_S(" cannot open chemtable file: " << filename);
    
    assert(fp);
    
    // use the header in the core IO ...
   

    // (1) write ugp magic number and version
    {
      int tmpint[2]={CHEM_IO_MAGIC_NUMBER,CHEM_IO_VERSION};
      size_t dummy= fwrite(tmpint,sizeof(int),2,fp);
      assert(dummy==2);
    }
    

    // (2) write the table type and dimensions
    {
      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      size_t dummy = tabletype.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and bye skip
      header.id = UGP_IO_CT_CART_3D;
      header.skip = sizeof(Header);
      // table dimensions
      header.idata[0] = nz; 
      header.idata[1] = nzvar;
      header.idata[2] = nc;
      header.idata[3] = nvars;
      // table pressure and temperature
      header.rdata[0] = pressure;
      header.rdata[1] = Tfuel;
      header.rdata[2] = Toxi;
      // write into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
    }

    // (3) write "Z" coordinate
    {
      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      string tmp    = "Z";
      size_t dummy = tmp.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and bye skip
      header.id = UGP_IO_CT_COOR;
      header.skip = nz * sizeof(double) + sizeof(Header);
      // table dimensions
      header.idata[0] = nz; 
      // write the into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
      // and the coordinate
      dummy= fwrite(Z,sizeof(double),nz,fp);
      assert(dummy==nz);
    }

    // (4) write "Zvar" coordinate
    {
      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      string tmp   = "Zvar";
      size_t dummy = tmp.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and bye skip
      header.id = UGP_IO_CT_COOR;
      header.skip = nzvar * sizeof(double) + sizeof(Header);
      // table dimensions
      header.idata[0] = nzvar; 
      // write the into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
      // and the coordinate
      dummy= fwrite(Zvar,sizeof(double),nzvar,fp);
      assert(dummy==nzvar);
    }

    // (5) write "C" coordinate
    {
      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      string tmp   = "C";
      size_t dummy = tmp.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and bye skip
      header.id = UGP_IO_CT_COOR;
      header.skip = nc * sizeof(double) + sizeof(Header);
      // table dimensions
      header.idata[0] = nc; 
      // write the into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
      // and the coordinate
      dummy= fwrite(C,sizeof(double),nc,fp);
      assert(dummy==nc);
    }
    


    // (6) write the data
    for(int ivar=0; ivar<nvars; ++ivar){

      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      size_t dummy = tableVarsNameVec[ivar].copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and bye skip
      header.id = UGP_IO_CT_DATA;
      header.skip = nz*nzvar*nc * sizeof(double) + sizeof(Header);
      // table dimensions
      header.idata[0] = nz*nzvar*nc; 
      header.idata[1] = nz; 
      header.idata[2] = nzvar; 
      header.idata[3] = nc; 
      // write the header into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
      // and the data
      for (int iz=0; iz<nz; ++iz)
	for (int izvar=0; izvar<nzvar; ++izvar)
	  for (int ic=0; ic<nc; ++ic) {
	    double tmp = (*table)(ivar,iz,izvar,ic);
	    dummy= fwrite(&tmp,sizeof(double),1,fp);
	    assert(dummy==1);
	  }
      
    }

    // end of file header
    {
      Header header;
      header.id = UGP_IO_EOF;
      sprintf(header.name,"EOF");
      header.skip = sizeof(Header);
      // write the header into the file      
      size_t dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);  
    }

    // close the file
    fclose(fp);
  
  }

private:

  // ================================
  // write the table to a binary file
  // ================================
  void writeBinaryVidaOld(){

    // length of the string
    // this should be unique for all the files
    // either on the write side or the read side
    int stringLength = 64;
    char * tmpchar = new char[stringLength];
    //char tmpchar[64];
    size_t dummy; 


    string filename = tablename + ".VidaNFPVCart3D.chem";
    cout << " > writing table to the vida compatible binary file: " << filename << endl;
    
    // open the file
    FILE * fp = NULL;
    fp=fopen(filename.c_str(),"wb");
    if (fp == NULL)
      CERR_S(" cannot open chemtable file: " << filename);
      
    assert(fp);

    // (1) write the string length 
    {        
      size_t dummy= fwrite(&stringLength,sizeof(int),1,fp);
      assert(dummy==1);
    }


    // (2) write the table type
    {
      char * tmpchar = new char[stringLength];
      for(int i=0; i<stringLength; ++i)
	tmpchar[i] = '\0';	
      dummy = tabletype.copy(tmpchar,stringLength);
      assert(dummy<stringLength);
      assert(tmpchar[dummy] == '\0');
      dummy = fwrite(tmpchar,sizeof(char),stringLength,fp);
      assert(dummy == stringLength);
      delete[] tmpchar;
    }

    

    // (3) write the name of the independent variables
    {
      for(int i=0; i<stringLength; ++i)
	tmpchar[i] = '\0';
      // Z mean
      string varname = "Z";
      dummy = varname.copy(tmpchar,stringLength);
      assert(dummy<stringLength);
      assert(tmpchar[dummy] == '\0');
      dummy = fwrite(tmpchar,sizeof(char),stringLength,fp);
      assert(dummy == stringLength);
      // Z variance
      for(int i=0; i<stringLength; ++i)
	tmpchar[i] = '\0';
      varname = "Zvar";
      dummy = varname.copy(tmpchar,stringLength);
      assert(dummy<stringLength);
      assert(tmpchar[dummy] == '\0');
      dummy = fwrite(tmpchar,sizeof(char),stringLength,fp);
      assert(dummy == stringLength);
      // Progress variable
      for(int i=0; i<stringLength; ++i)
	tmpchar[i] = '\0';
      varname = "C";
      dummy = varname.copy(tmpchar,stringLength);
      assert(dummy<stringLength);
      assert(tmpchar[dummy] == '\0');
      dummy = fwrite(tmpchar,sizeof(char),stringLength,fp);
      assert(dummy == stringLength);
    }

    
    // (4) write table dimensions
    {    
      // Z mean    
      dummy= fwrite(&nz,sizeof(int),1,fp);
      assert(dummy==1);
      // Z variance
      dummy= fwrite(&nzvar,sizeof(int),1,fp);
      assert(dummy==1);
      // progress variable
      dummy= fwrite(&nc,sizeof(int),1,fp);
      assert(dummy==1);
    }

    // (5) write table coordinates
    {    
      // Z mean    
      dummy= fwrite(Z,sizeof(double),nz,fp);
      assert(dummy==nz);
      // Z variance
      dummy= fwrite(Zvar,sizeof(double),nzvar,fp);
      assert(dummy==nzvar);
      // progress variable
      dummy= fwrite(C,sizeof(double),nc,fp);
      assert(dummy==nc);
    }

    // (6) write the variables
    for(int ivar=0; ivar<nvars; ++ivar){
      // variable name  
      for(int i=0; i<stringLength; ++i)
	tmpchar[i] = '\0';   
      dummy = tableVarsNameVec[ivar].copy(tmpchar,stringLength);
      assert(dummy<stringLength);
      assert(tmpchar[dummy] == '\0');
      dummy = fwrite(tmpchar,sizeof(char),stringLength,fp);
      assert(dummy == stringLength);
      // and the data
      for (int ic=0; ic<nc; ++ic)
	for (int izvar=0; izvar<nzvar; ++izvar)
	  for (int iz=0; iz<nz; ++iz) {
	    double tmp = (*table)(ivar,iz,izvar,ic);
	    dummy= fwrite(&tmp,sizeof(double),1,fp);
	    assert(dummy==1);
	  }
    }


    delete[] tmpchar;
    fclose(fp);
    
  }

 
  // =================================
  // write the table to a tecplot file
  // =================================
  void dumpTecplotFile() {

    string filename = tablename + ".dat";

    cout << " > writing table to the tecplot file: " << filename << endl;
          
    FILE * fp = fopen(filename.c_str(),"w");
    fprintf(fp,"title = \"chemtable\"\n");
    fprintf(fp,"variables = \"Z\"\n\"Zvar\"\n\"C\"\n");
    for (int ivar=0; ivar<tableVarsNameVec.size(); ++ivar)
      fprintf(fp,"\"%s\"\n",tableVarsNameVec[ivar].c_str());
    fprintf(fp,"zone i = %d, j = %d, k = %d, DATAPACKING=POINT\n",nz,nzvar,nc);
    for (int k=0; k <nc; ++k) 
      for (int j=0; j <nzvar; ++j) 
	for (int i=0; i <nz; ++i){ 
	  fprintf(fp,"%g %g %g",Z[i],Zvar[j],C[k]);
	  for (int ivar=0; ivar<tableVarsNameVec.size(); ++ivar)
	    fprintf(fp," %g",(*table)(ivar,i,j,k));
	  fprintf(fp,"\n");
	}
    fclose(fp);
    

  }

 
  // ==============================================
  // read the flamelets files and sort base on 
  // maximum progress variable
  // ==============================================
  void readFlameletFilesAndSort() {

    assert(flameletList.size() == 0 );

    string filename = getStringParam("CHEMTABLE.FLAMELET_LIST", "FlameletList.txt");
    cout << " > reading the flamelet list file: " << filename << endl;

    // open the file to get the size
    ifstream ifp;
    // open the file
    ifp.open(filename.c_str());
    if ( ifp.fail() ) 
      CERR_S("could not open flamelet file list: " << filename );
    string buffer_str;
    int nfl = 0;
    while(getline(ifp, buffer_str)) {
      if (buffer_str.empty())
	break;
      ++nfl;
    }
    ifp.close();
    cout << "  found " << nfl << " flamelet files" << endl;
    if ( nfl <1 )
      CERR_S("at least one flamelet file is required" );
	cout << " > reading the flamelet files ..." << endl;
 
    // open the file again
    ifp.open(filename.c_str());
    if ( ifp.fail() )
      CERR_S("could not open flamelet file list: " << filename );
	
    cout << "  progress:  ";
    cout.flush();
    int ifl = 0;
    int tmpint = max((int)(0.1*nfl),1);
    while(getline(ifp, buffer_str)) {
      if ( (ifl+1)%tmpint == 0 ) {
	cout << (ifl+1)*10/tmpint << "% ... ";
	cout.flush();
      }
      if (buffer_str.empty())
	break;
      istringstream buf(buffer_str);
      string thisfile;
      buf >> thisfile;
      // read the flamelet file
      NonPremixedFlamelet * flamelet = new NonPremixedFlamelet();
      flamelet->init(thisfile);
      flameletList.push_back(flamelet);
      ++ifl;
    }
    ifp.close();
    cout << endl;
    // check
    assert(ifl == nfl);
  
    // and sort based on Cmax
    flameletList.sort(isProgSmallerFlamelet<NonPremixedFlamelet>);

    // dump the list
    if ( checkParam("CHEMTABLE.DUMP_FLAMELETLIST") ) {
      cout << " > flamelet list before clean up ..." << endl;
      for (list<NonPremixedFlamelet*>::iterator flamelet=flameletList.begin(); flamelet!=flameletList.end(); ++flamelet) 
	cout << (*flamelet)->filename << " ";
      cout << endl;
    }

    // dump some info on screen
    if (cti_verbose) {
      list<NonPremixedFlamelet*>::iterator flamelet=flameletList.begin();
      (*flamelet)->info();
    }
    
    {
      // check that the pressure,and temperature at the fuel and oxidizer are the same for all the flamelets and store that pressure
      list<NonPremixedFlamelet*>::iterator flamelet=flameletList.begin();
      pressure = (*flamelet)->pressure;
      Tfuel    = (*flamelet)->Tfuel;
      Toxi     = (*flamelet)->Toxi;
      for (flamelet=++(flameletList.begin()); flamelet!=flameletList.end(); ++flamelet) {
	assert(pressure == (*flamelet)->pressure );
	assert(Tfuel == (*flamelet)->Tfuel );
	assert(Toxi == (*flamelet)->Toxi );
      }
    }
    
  }


  // ===========
  // Z mean grid
  // ===========
  void buildZgrid() {

    cout << " > building the Z grid ..." << endl;

    if (Param * param = getParam("CART_CHEMTABLE.STRETCH_Z")) {
      int np=4;
      double *xp = new double[np];
      xp[0] = 0.0;
      xp[1] = max(param->getDouble(0), 0.0 + 1.0E-8);
      xp[2] = min(param->getDouble(1), 1.0 - 1.0E-8);
      xp[3] = 1.0;
      double frac  = param->getDouble(2);

      assert( (xp[0] < xp[1]) &&
              (xp[1] < xp[2]) &&
              (xp[2] < xp[3]) &&
              (frac  >  0.0 ) &&
              (frac  <  1.0 ) );

      double *qp = new double[np];
      qp[0] = 0.0;
      qp[1] = (1.0-frac)*(xp[1]-xp[0])/(xp[1]-xp[0] + xp[3]-xp[2]);
      qp[2] = qp[1]+frac;
      qp[3] = 1.0;

      stretchgrid(xp,qp,np,Z,nz);
      delete [] xp;
      delete [] qp;
    }
    else { // use "average" Z mesh from all flamelets

      double *z = new double[nz]; 
      for (int iz=0; iz<nz; ++iz) Z[iz] = 0.0;

      for (list<NonPremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f) { 

        double *zz = (*f)->getVarPtr("Z");
        int     nn = (*f)->npoints;
        double  dq = 1.0/double(nn-1);

        if (nn != nz) {
          // resize grid to nz points, preserve relative spacing
          for (int iz=0; iz<nz; ++iz) { 
            double q  = double(iz)/double(nz-1);  
            int    i  = fmin(nn-2,floor(q/dq)); 
            double w  = double(i+1) - q/dq;
            z[iz]     = w*zz[i] + (1.0-w)*zz[i+1];
          }
        }
        else
          for (int iz=0; iz<nz; ++iz) z[iz] = zz[iz]; 

        // sum for averaging
        for (int iz=0; iz<nz; ++iz) Z[iz] += z[iz]; 
      }
      for (int iz=0; iz<nz; ++iz) Z[iz] /= double(flameletList.size());

      delete [] z; 
    }
    //else // uniform grid
    //  for (int i=0; i<nz; i++)
    //    Z[i] = double(i)/double(nz-1);
  }
  

  // ===============
  // Z variance grid
  // ===============
  void buildZvargrid() {
    
    cout << " > building the Zvar grid ..." << endl;   

    // Quadratic mesh between 0 and 0.25 (variance cannot be larger than Z*(1-Z))
    for (int i=0; i<nzvar; i++)
      Zvar[i] = 0.25 * pow((double)i / (double)(nzvar-1), 2.0);
        
  }  

  
  // ==============================================
  // clean up data list by making the 
  // mapping unique in (Z,Zvar,C) coordinate system
  // for zvar == 0   <=== ???
  // ==============================================
   void cleanupFlameletList() {

    cout << " > cleaning up Flameletlist ... " << endl;

    vector<string> deletedFlameletVec;

    int nflamelet0 = flameletList.size();
    int count =0;

    list<NonPremixedFlamelet*>::iterator flamelet   = flameletList.begin();
    list<NonPremixedFlamelet*>::iterator flameletp1 = flamelet;
    ++flameletp1;
    cout << "  progress:  ";
    cout.flush();
    int whilecount = 1; 
    int nfl = flameletList.size();
    int tmpint = max((int)(0.1*nfl),1);
    while ( flameletp1!=flameletList.end() ) {
      if ( (whilecount+1)%tmpint == 0 ) {
	cout << (whilecount+1)*10/tmpint << "% ... ";
	cout.flush();
      }
      double error = 0.0;
      double tol2  = 1.0e-4;
      double tol1  = 1.0e-6;
      for (int iz=0; iz<nz; ++iz){
	double Ci    = (*flamelet)->getValGivenMeanAndVar(Z[iz],Z[0],"prog");	
	double Cip1  = (*flameletp1)->getValGivenMeanAndVar(Z[iz],Z[0],"prog");
	double Cdiff = Cip1 - Ci;
	double Cmean = max(0.5*(Cip1+Ci),tol1);
	error = max(-Cdiff/Cmean,error);
      }
      
      //cout << "DEBUG: " << whilecount << " flamelet   : " << (*flamelet)->filename << endl;
      //cout << "       " << whilecount << " flameletp1 : " << (*flameletp1)->filename << endl;
      //cout << "       " << whilecount << " error      : " << error << endl;
      if ( error>tol2 ) {
        /*
	cout << "   >> inconsistency between  the flamelets: " << endl;
	cout << "        " << (*flamelet)->filename << endl;
	cout << "          Tmax:    " << (*flamelet)->Tmax << endl;
	cout << "          Cmax:    " << (*flamelet)->Cmax << endl;
	cout << "          chi_st:  " << (*flamelet)->chi_st << endl;
	cout << "        " << (*flameletp1)->filename << endl;
	cout << "          Tmax:    " << (*flameletp1)->Tmax << endl;
	cout << "          Cmax:    " << (*flameletp1)->Cmax << endl;
	cout << "          chi_st:  " << (*flameletp1)->chi_st << endl;
	cout << "   >> with maximum relative error(%) = " << error*100.0 << endl;
	cout << " *** removing flamelet file " << (*flameletp1)->filename << endl;
        */
	deletedFlameletVec.push_back((*flameletp1)->filename);
	delete((*flameletp1));
	flameletp1 = flameletList.erase(flameletp1);
	++count;
      }
      else{
	flamelet = flameletp1;
	++flameletp1;
      }
      ++whilecount;
    }
    cout << endl;

    // check
    assert((flameletList.size()+count)==nflamelet0);
    cout << "  >> " << count << " flamelets were deleted from the list: " << endl; 
    for (int i=0; i<deletedFlameletVec.size(); ++i)
      cout << "     " << deletedFlameletVec[i] << endl;

    //    if (count!=0)
    //  cleanupFlameletList();
   
    // dump the list
    if ( checkParam("CHEMTABLE.DUMP_FLAMELETLIST") ) { 
      cout << " > flamelet list after clean up ..." << endl;
      for (list<NonPremixedFlamelet*>::iterator flamelet=flameletList.begin(); flamelet!=flameletList.end(); ++flamelet) 
	cout << (*flamelet)->filename << " ";
      cout << endl;
    }
    
   }


  // ======================
  // progress variable grid
  // ======================
  void buildCgrid() {

    cout << " > building the C grid ..." << endl;

    assert(flameletList.size() > 0 );

    // first get the min and max values for C
    double Cmin = 1.0e20;
    double Cmax = -1.0e20;
    for (list<NonPremixedFlamelet*>::iterator flamelet=flameletList.begin(); flamelet!=flameletList.end(); ++flamelet) {
      Cmin = min(Cmin,(*flamelet)->Cmin);
      Cmax = max(Cmax,(*flamelet)->Cmax);
    }
    if (fabs(Cmin)<1.0e-14) Cmin = 0.0;
    assert(Cmax>Cmin);

    if (Param * param = getParam("CART_CHEMTABLE.STRETCH_C")) {
      int np=4;
      double *xp = new double[np];
      xp[0] = Cmin;
      xp[1] = max(param->getDouble(0), Cmin + 1.0E-8);
      xp[2] = min(param->getDouble(1), Cmax - 1.0E-8);
      xp[3] = Cmax;
      double frac  = param->getDouble(2);

      assert( (xp[0] < xp[1]) &&
              (xp[1] < xp[2]) &&
              (xp[2] < xp[3]) &&
              (frac  >  0.0 ) &&
              (frac  <  1.0 ) );

      double *qp = new double[np];
      qp[0] = 0.0;
      qp[1] = (1.0-frac)*(xp[1]-xp[0])/(xp[1]-xp[0] + xp[3]-xp[2]);
      qp[2] = qp[1]+frac;
      qp[3] = 1.0;

      stretchgrid(xp,qp,np,C,nc);
      delete [] xp;
      delete [] qp;
    }
    else // uniform grid
      for (int i=0; i<nc; i++)
        C[i] = Cmin + double(i)/double(nc-1)*(Cmax-Cmin);

    cout << "   using (Cmin,Cmax) = (" << C[0] << "," << C[nc-1] << ") to build the C coordinate" << endl;

  }



  // ======================
  // build the C mapping
  // ======================
  void buildCmapping(){

    cout << " > building C mapping ..." << endl;

    assert(table == NULL);
    table = new CTIarray<double>(nvars,nz,nzvar,nc);
    
    // get the progress variable index
    int prog_index = -1;
    for (int i = 0; i<tableVarsNameVec.size(); ++i) 
      if ( tableVarsNameVec[i] == "prog" ) 
	prog_index = i;
    if ( prog_index < 0 ) 
      CERR_S(" could not find prog");
        
    // number of flamelets
    int nfl = flameletList.size();
           
    // take advantage of flamelets being sorted ...
    cout << "  progress:  ";
    cout.flush();
    int tmpint = max((int)(0.1*(double)nzvar),1);
    for (int izvar=0; izvar<nzvar; ++izvar) {
      if ( (izvar+1)%tmpint == 0 ) {
	cout << (izvar+1)*10/ tmpint<< "% ... ";
	cout.flush();
      }
      for (int iz=0; iz<nz; ++iz) {

	if ( nfl == 1 ) {
	  for (int ivar=0; ivar<nvars; ++ivar) 
	    for (int ic=0; ic<nc; ++ic)
	      (*table)(ivar,iz,izvar,ic) = (*data4d)(ivar,iz,izvar,0);
	}
	else{
	  // sort based on progress variable ...
	  // first store the ifl and prog
	  vector<pairIntDbl> myPairVec;
	  for (int ifl=0; ifl<nfl; ++ifl){
	    double prog = (*data4d)(prog_index,iz,izvar,ifl);
	    myPairVec.push_back(pairIntDbl(ifl,prog));
	  }
	  assert(myPairVec.size() == nfl);
	  // now sort the list
	  std::sort(myPairVec.begin(),myPairVec.end(),isSecondSmallerPair);
	  
	  // first find the index
	  int ic  = 0;
	  int ifl = 1;
	  // check
	  assert(C[0]<=myPairVec[0].second);
	  // now for all Cs
	  while ( ic != nc ) {
	    // find the flamelet that C is just smaller than the flamelet's prog
	    while ( C[ic] > myPairVec[ifl].second ) {
	      if ( ifl == (nfl-1) ) break;
	      ++ifl;
	    }
	    // check
	    assert(ifl>0);
	    assert(ifl<=(nfl-1));
	    // clipp C to flamelets min and max C
	    double thisC = C[ic];
	    thisC = max(myPairVec[0].second,thisC);
	    thisC = min(myPairVec[nfl-1].second,thisC);
	    // check
	    assert(thisC>=myPairVec[ifl-1].second);
	    assert(thisC<=myPairVec[ifl].second);
	    
	    // now C is between myPairVec[ifl-1].second <= C     <= myPairVec[ifl].second
	    // and              myPairVec[0].second     <= thisC <= myPairVec[nfl-1].second
	    // the actual indecise are
	    int iflm1 = myPairVec[ifl-1].first;
	    int iflp1 = myPairVec[ifl].first;
	    // just a check
	    double Cp1 = (*data4d)(prog_index,iz,izvar,iflp1);
	    double Cm1 = (*data4d)(prog_index,iz,izvar,iflm1);
	    assert(thisC<=Cp1);
	    assert(thisC>=Cm1);
	    double Cdiff = Cp1 - Cm1;
	    double wm1;
	    double wp1;
	    if ( Cdiff > 0.0 ) {
	      wm1 = (Cp1-thisC) / Cdiff;
	      wp1 = 1.0 - wm1;
	    }
	    else if ( Cdiff == 0.0 ) {
	      wm1 = 0.5;
	      wp1 = 0.5;
	    }
	    else {
	      CERR_S(" Cdiff should be non-negative");
	    }
	    
	    // now interpolate
	    for (int ivar=0; ivar<nvars; ++ivar) 
	      (*table)(ivar,iz,izvar,ic) = wm1 * (*data4d)(ivar,iz,izvar,iflm1) + wp1 * (*data4d)(ivar,iz,izvar,iflp1);
	  
	    ++ic;
	  } // while ( ic != nc )
	} //if ( nfl == 1 ) {	

      } // for (int iz=0; iz<nz; ++iz)
    } //  for (int izvar=0; izvar<nzvar; ++izvar)
    cout << endl; 

  }

  void dumpPointCloud(){ //lee

    cout << " > dumping flamelet point cloud ..." << endl;

    // get variable indices
    int prog_index = -1;
    int rho_index  = -1;
    int src_index  = -1;
    int tmp_index  = -1;
    int yco_index  = -1;
    for (int i = 0; i<tableVarsNameVec.size(); ++i) {
      if ( tableVarsNameVec[i] == "prog" ) {
	prog_index = i;
        cout << "   > found prog, i = " << i << endl;
      }
      else if ( tableVarsNameVec[i] == "rho" ) {
	rho_index = i;
        cout << "   > found rho,  i = " << i << endl;
      }
      else if ( tableVarsNameVec[i] == "src_prog" ) {
	src_index = i;
        cout << "   > found src,  i = " << i << endl;
      }
      else if ( tableVarsNameVec[i] == "Y_CO" ) {
	yco_index = i;
        cout << "   > found yco,  i = " << i << endl;
      }
      else if ( tableVarsNameVec[i] == "T" ) {
	tmp_index = i;
        cout << "   > found T,    i = " << i << endl;
      }
    }
        
    // number of flamelets
    int nfl = flameletList.size();

    // output file
    ofstream fout;
    fout.open("flameletPointCloud.dat");
    fout.precision(10);
           
    // take advantage of flamelets being sorted ...
    cout << "  progress:  ";
    cout.flush();
    int tmpint = max((int)(0.1*(double)nzvar),1);
    for (int izvar=0; izvar<nzvar; ++izvar) {
      if ( (izvar+1)%tmpint == 0 ) {
	cout << (izvar+1)*10/ tmpint<< "% ... ";
	cout.flush();
      }
      for (int iz=0; iz<nz; ++iz) {

        if ( Zvar[izvar] > Z[iz]*(1.0-Z[iz]) ) continue;

	for (int ifl=0; ifl<nfl; ++ifl){
	  double prog = (*data4d)(prog_index,iz,izvar,ifl);
          double rho  = (*data4d)( rho_index,iz,izvar,ifl);
          double src  = (*data4d)( src_index,iz,izvar,ifl);
          double tmp  = (*data4d)( tmp_index,iz,izvar,ifl);
          double yco  = (*data4d)( yco_index,iz,izvar,ifl);

          fout << Z[iz] << "  " << Zvar[izvar] << "  " << prog << "  "
               << rho << "  " << src << "  " << tmp << "  " << yco << "  " 
               << ifl << endl; 
	} // ifl

      } // iz
    } //  izvar

    fout.close();
    cout << endl; 

    cout << "shutting down" << endl;
    assert(0);

  }
  
};




// =====================================
// Chemtable compatible with vida solver
// =====================================

class VidaChemtableNFPVCart3D: public ChemtableNFPVCart3D {

private:
  bool nox_model;

public: 

  VidaChemtableNFPVCart3D(const string& tabletype) : ChemtableNFPVCart3D() {
    cout << "VidaChemtableNFPVCart3D()" << endl;
    this->tabletype = tabletype;
    filename = tablename + ".VidaNFPVCart3D.chem"; // tablename is initialized at the base class level
    nox_model = checkParam("NOX_MODEL");
  }
  
  virtual ~VidaChemtableNFPVCart3D() {
    
    cout << "~VidaChemtableNFPVCart3D()" << endl;
    
  }
  
  
  // =====================================
  // variables to be included in the table
  // =====================================   
  
  virtual void initTableVarsNameVec() {
    
    assert(tableVarsNameVec.size() == 0);
    
    // a vida table requires at least 8 variables
    //  (1)  density
    //  (2)  temperature
    //  (3)  viscosity
    //  (4)  diffusivity  
    //  (5)  progress variable
    //  (6)  progress variable source term
    //  (7)  heat capacity ratio (gamma)
    //  (8)  molecular weight
    tableVarsNameVec.push_back("rho");
    tableVarsNameVec.push_back("T");
    tableVarsNameVec.push_back("mu");
    tableVarsNameVec.push_back("locp");
    tableVarsNameVec.push_back("prog");
    tableVarsNameVec.push_back("src_prog");
    tableVarsNameVec.push_back("gamma");
    tableVarsNameVec.push_back("mw");
    // add vars for Ihme & Pitsch NOx model
    if (nox_model) {
      tableVarsNameVec.push_back("src_posNO");      // production term
      tableVarsNameVec.push_back("src_negNO");      // consumption term
      tableVarsNameVec.push_back("src_negNOdivNO"); // consumption term divided by mass fraction
      tableVarsNameVec.push_back("Y_NO");           // mass fraction
    }
    // add vars for liquid fuel
    if (checkParam("LIQUID_FUEL")) {
      //tableVarsNameVec.push_back("lsp_rho_fuel");
      //tableVarsNameVec.push_back("lsp_cp_fuel");
      //tableVarsNameVec.push_back("lsp_mu_fuel");
      //tableVarsNameVec.push_back("lsp_cond_fuel");
      tableVarsNameVec.push_back("lsp_rho_nonfuel");
      tableVarsNameVec.push_back("lsp_cp_nonfuel");
      tableVarsNameVec.push_back("lsp_mu_nonfuel");
      tableVarsNameVec.push_back("lsp_cond_nonfuel");
      tableVarsNameVec.push_back("lsp_mw_nonfuel");
      tableVarsNameVec.push_back("lsp_Tref");
    }
    // and add the user defined species
    if (checkParam("CHEMTABLE.SPECIES")) 
    if (Param * param = getParam("CHEMTABLE.SPECIES")) 
      for (int i = 0; i<param->size(); ++i) 
	      tableVarsNameVec.push_back("Y_"+param->getString(i));

    // print on screen
    cout << " > species to be included in the table: " << endl;
    for (int i = 0; i<tableVarsNameVec.size(); ++i) 
      cout << "       " << tableVarsNameVec[i] << endl;

    nvars = tableVarsNameVec.size();
    assert(nvars>=6);

  }
  
  
  // =======================
  // compute the convolution
  // =======================
  virtual void calcConvolution() {

    const double SMALL = 1.0E-10;
    cout << " > computing convolution with presumed pdf ... " << endl;
  
    assert(data4d == NULL);
    data4d = new CTIarray<double>(nvars,nz,nzvar,flameletList.size());

    // for all the flamelets
    int ifl = 0;
    cout << "  progress:  ";
    cout.flush();
    int tmpint = max((int)(0.1*(double)flameletList.size()),1);
    for (list<NonPremixedFlamelet*>::iterator flamelet=flameletList.begin(); flamelet!=flameletList.end(); ++flamelet){
      if ( (ifl+1)%tmpint == 0 ) {
	cout << (ifl+1)*10/tmpint << "% ... ";
	cout.flush();
      }
      for (int izvar=0; izvar<nzvar; ++izvar) {
	double thisZvar = Zvar[izvar];
      	for (int iz=0; iz<nz; ++iz) {
	  double thisZ = Z[iz];

          double srcNegNO = 0.0;
          if (nox_model) {
            // source term for Ihme & Pitsch 2008 NOx model
            double NO=(*flamelet)->getValGivenMeanAndVar(thisZ,thisZvar,"Y_NO");
            srcNegNO=(*flamelet)->getValGivenMeanAndVar(thisZ,thisZvar,"src_negNO");
            if ( NO <= SMALL ) srcNegNO = 0.0;
            else  srcNegNO /= NO;
          }

          for (int ivar=0; ivar<nvars; ++ivar) {
            string varname = tableVarsNameVec[ivar];
//            cout << "$$$$$$$$$$$$$$$$ I am before 2499; var = " << varname << "; flamelet = " << (*flamelet)->filename << endl;  // XXX ars
            double value = 1.0e+24;
            if      ( varname == "src_negNOdivNO") value = srcNegNO;
            else value =(*flamelet)->getValGivenMeanAndVar(thisZ,thisZvar,varname);
            // and set the data
            (*data4d)(ivar,iz,izvar,ifl) = value;

          /* old stuff
	  for (int ivar=0; ivar<nvars; ++ivar) {
	    string varname = tableVarsNameVec[ivar];
	    (*data4d)(ivar,iz,izvar,ifl) = 
	      (*flamelet)->getValGivenMeanAndVar(thisZ,thisZvar,varname);
          */

	  } // for (int ivar=0; ivar<nvars; ++ivar)
	} // for (int iz=0; iz<nz; ++iz) 
      } // for (int izvar=0; izvar<nzvar; ++izvar)

      ++ifl;
    }
    cout << endl;
    assert(ifl==flameletList.size());

    // we do not need the flamelets anymore
    for (list<NonPremixedFlamelet*>::iterator flamelet=flameletList.begin(); flamelet!=flameletList.end(); ++flamelet) {
      delete (*flamelet);
      (*flamelet) = NULL;
    }
    
  }

};



// ======================================
// Chemtable compatible with chris solver
// ======================================

class ChrisChemtableNFPVCart3D: public ChemtableNFPVCart3D {

private:

  Mixture * mixture;

public: 
  
  ChrisChemtableNFPVCart3D(const string& tabletype) : ChemtableNFPVCart3D() {
    cout << "ChrisChemtableNFPVCart3D()" << endl;
    this->tabletype = tabletype;
    filename = tablename + ".ChrisNFPVCart3D.chem";
    mixture = NULL;

  }
  
  virtual ~ChrisChemtableNFPVCart3D() {
    
    cout << "~ChrisChemtableNFPVCart3D()" << endl;
    
    if ( mixture != NULL )
      delete mixture;
    
  }


  // =====================================
  // variables to be included in the table
  // =====================================   
  
  virtual void initTableVarsNameVec() {
    
    assert(tableVarsNameVec.size() == 0);
    
    // a chris table requires at least 12 variables
    //  (1)  density,      
    //  (2)  temperature,
    //  (3)  gas constant,
    //  (4)  internal energy,
    //  (5)  progress variable,
    //  (6)  progress variable source term,
    //  (7)  specific heat ratio,
    //  (9)  coefficient for specific heat ratio,
    //  (9)  viscosity,
    //  (10) coefficent for the viscosity,
    //  (11) diffusivity, and
    //  (12) coefficient for diffusivity,

    tableVarsNameVec.push_back("rho");
    tableVarsNameVec.push_back("T");
    tableVarsNameVec.push_back("R");
    tableVarsNameVec.push_back("e");
    tableVarsNameVec.push_back("prog");
    tableVarsNameVec.push_back("src_prog");
    
    tableVarsNameVec.push_back("gamma");
    tableVarsNameVec.push_back("a_gamma");
    
    tableVarsNameVec.push_back("mu");
    tableVarsNameVec.push_back("a_mu");
    
    tableVarsNameVec.push_back("locp");
    tableVarsNameVec.push_back("a_locp");

    // add vars for liquid fuel
    // NOTE: Thermo quantities are evaluated at Tref, so no adjustments are
    //       necessary. Eventually, we should correct for pressure.
    if (checkParam("LIQUID_FUEL")) {
      //tableVarsNameVec.push_back("lsp_rho_fuel");
      //tableVarsNameVec.push_back("lsp_cp_fuel");
      //tableVarsNameVec.push_back("lsp_mu_fuel");
      //tableVarsNameVec.push_back("lsp_cond_fuel");
      tableVarsNameVec.push_back("lsp_rho_nonfuel");
      tableVarsNameVec.push_back("lsp_cp_nonfuel");
      tableVarsNameVec.push_back("lsp_mu_nonfuel");
      tableVarsNameVec.push_back("lsp_cond_nonfuel");
      tableVarsNameVec.push_back("lsp_mw_nonfuel");
      tableVarsNameVec.push_back("lsp_Tref");
    }
    
    // and add the user defined species
    if (Param * param = getParam("CHEMTABLE.SPECIES"))
      for (int i = 0; i<param->size(); i++)
	      tableVarsNameVec.push_back("Y_" + param->getString(i));

    // print on screen
    cout << " > species to be included in the table: " << endl;
    for (int i = 0; i<tableVarsNameVec.size(); ++i) 
      cout << "       " << tableVarsNameVec[i] << endl;

    nvars = tableVarsNameVec.size();
    assert(nvars>=6);

  }


  // =======================
  // compute the convolution
  // =======================
  virtual void calcConvolution() {

    cout << " > computing convolution with presumed pdf ... " << endl;
  
    // perturbation in the energy
    double rhoDeltaE = 5000.0;
    cout << "  >> perturbing flamelets internal energy (rhoe) by  Delta(rhoe) = " << rhoDeltaE << endl;

    assert(data4d == NULL);
    data4d = new CTIarray<double>(nvars,nz,nzvar,flameletList.size());
    
    // for all the flamelets
    int ifl = 0;
    int tmpint = max((int)(0.1*(double)flameletList.size()),1);
    for (list<NonPremixedFlamelet*>::iterator flamelet=flameletList.begin(); flamelet!=flameletList.end(); ++flamelet){
      if ( (ifl+1)%tmpint == 0 ) {
	cout << (ifl+1)*10/tmpint << "% ... ";
	cout.flush();
      }
      
      // load the mixture and build the necessary elements
      (*flamelet)->addMixtureERandPerturbedVariables(rhoDeltaE);

      if ( ifl == 0 )  
	cout << "  convolution progress:  ";

      for (int izvar=0; izvar<nzvar; ++izvar) {
	double thisZvar = Zvar[izvar];
      	for (int iz=0; iz<nz; ++iz) {
	  double thisZ = Z[iz];
	  // compute the variables
	  double rho, T, R, e, prog, src_prog, gamma, a_gamma, mu, a_mu, locp, a_locp;
	  getVarsGivenZandZvar(rho,T,R,e,prog,src_prog,gamma,a_gamma,mu,a_mu,
			       locp,a_locp,thisZ,thisZvar,rhoDeltaE,(*flamelet));
	  // now for the variables
	  for (int ivar=0; ivar<nvars; ++ivar) {
	    string varname = tableVarsNameVec[ivar];
	    double value = 1.0e+24;
	    if      ( varname == "rho"      ) value = rho;
	    else if ( varname == "T"        ) value = T;
	    else if ( varname == "R"        ) value = R;
	    else if ( varname == "e"        ) value = e;
	    else if ( varname == "prog"     ) value = prog;
	    else if ( varname == "src_prog" ) value = src_prog;
	    else if ( varname == "gamma"    ) value = gamma;
	    else if ( varname == "a_gamma"  ) value = a_gamma;
	    else if ( varname == "mu"       ) value = mu;
	    else if ( varname == "a_mu"     ) value = a_mu;
	    else if ( varname == "locp"     ) value = locp;
	    else if ( varname == "a_locp"   ) value = a_locp;
	    else value =(*flamelet)->getValGivenMeanAndVar(thisZ,thisZvar,varname);
	    // and set the data
	    (*data4d)(ivar,iz,izvar,ifl) = value;
	    
	  } // for (int ivar=0; ivar<nvars; ++ivar)
	} // for (int iz=0; iz<nz; ++iz) 
      } // for (int izvar=0; izvar<nzvar; ++izvar)
      
      ++ifl;
    }
    cout << endl;
    assert(ifl==flameletList.size());
    
    // we do not need the flamelets anymore
    for (list<NonPremixedFlamelet*>::iterator flamelet=flameletList.begin(); flamelet!=flameletList.end(); ++flamelet) {
      delete (*flamelet);
      (*flamelet) = NULL;
    }
    
  }



  // ================================================
  // compute variables from flamelet given Z and Zvar
  // ================================================
  inline void getVarsGivenZandZvar(double& rho, double& T, double& R, double& e, 
				   double& prog, double& src_prog, double& gamma, 
				   double& a_gamma, double& mu, double& a_mu, 
				   double& locp, double& a_locp, 
				   const double& Z, const double& Zvar, 
				   const double& rhoDeltaE,  NonPremixedFlamelet * flamelet) {
    
    // get the variables that are already in the table 
    rho           = flamelet->getValGivenMeanAndVar(Z,Zvar,"rho");
    const double deltaE = rhoDeltaE/rho;
    R             = flamelet->getValGivenMeanAndVar(Z,Zvar,"R");
    double RT     = flamelet->getValGivenMeanAndVar(Z,Zvar,"RT");
    T             = RT / R;
    e             = flamelet->getValGivenMeanAndVar(Z,Zvar,"e");
    prog          = flamelet->getValGivenMeanAndVar(Z,Zvar,"prog");
    src_prog      = flamelet->getValGivenMeanAndVar(Z,Zvar,"src_prog");
    mu            = flamelet->getValGivenMeanAndVar(Z,Zvar,"mu");
    locp          = flamelet->getValGivenMeanAndVar(Z,Zvar,"locp");
    
    // get the plus and minus perturbatiosn
    double RTp    = flamelet->getValGivenMeanAndVar(Z,Zvar,"RTp");
    double Tp     = RTp / R;
    double MUp    = flamelet->getValGivenMeanAndVar(Z,Zvar,"MUp");
    double LOCPp  = flamelet->getValGivenMeanAndVar(Z,Zvar,"LOCPp");
    
    double RTm    = flamelet->getValGivenMeanAndVar(Z,Zvar,"RTm");
    double Tm     = RTm / R;
    double MUm    = flamelet->getValGivenMeanAndVar(Z,Zvar,"MUm");
    double LOCPm  = flamelet->getValGivenMeanAndVar(Z,Zvar,"LOCPm");

    
    // derivatives: de/dT and d^2e/dT^2 at T=T
    double dTm    = T  - Tm;
    double dTp    = Tp - T;
    double dT     = Tp - Tm;
    double dedT   = ( dTm*dTm*(e+deltaE) + (dTp-dTm)*dT*e - dTp*dTp*(e-deltaE) ) / (dTp*dTm*dT);
    double d2edT2 = 2.0 * ( dTm*(e+deltaE) - dT*e + dTp*(e-deltaE) ) / (dTp*dTm*dT);
        
    // compute coefficients
    gamma         = R / dedT + 1.0;
    a_gamma       = - d2edT2 * (gamma-1.0) * (gamma-1.0) / R;

    double dTm_log = log(T/Tm);
    double dTp_log = log(Tp/T);
    double dT_log = log(Tp/Tm);

    /* DAP: corrected tabulation of a_mu, a_locp */
    a_mu          = ( dTm_log*dTm_log*log(MUp)   + (dTp_log-dTm_log)*dT_log*log(mu)   - dTp_log*dTp_log*log(MUm)   ) / (dTp_log*dTm_log*dT_log);
    a_locp        = ( dTm_log*dTm_log*log(LOCPp)   + (dTp_log-dTm_log)*dT_log*log(locp)   - dTp_log*dTp_log*log(LOCPm)   ) / (dTp_log*dTm_log*dT_log);
/*
    // this is adapted from Vincent's code. I think it is wrong
    // it should take the log of the MUp and ...
    a_mu          = ( dTm*dTm*MUp   + (dTp-dTm)*dT*mu   - dTp*dTp*MUm   ) / (dTp*dTm*dT);
    a_locp        = ( dTm*dTm*LOCPp + (dTp-dTm)*dT*locp - dTp*dTp*LOCPm ) / (dTp*dTm*dT);
*/
      
  }

};


// =============================================
// charles 2d (Z,C) premixed chemtable
// =============================================

class CharlesChemtablePFPVCart2D: public AbstractChemtable {

protected:

  list<PremixedFlamelet*> flameletList;
  vector<string> tableVarsNameVec;

  int nz; 
  int nc; 
  int nvars;
  double pressure;
  double * mixfrac;
  double * Z;
  double * C;
  CTIarray<double> * table;

public:

  CharlesChemtablePFPVCart2D(const string& tabletype) {

    cout << "CharlesChemtablePFPVCart2D()" << endl;
    this->tabletype = tabletype;
    filename = tablename + ".CharlesPFPVCart2D.chem"; // tablename is initialized at the base class level

    nz       = 0;
    nc       = 0;
    nvars    = 0;
    pressure = 0.0; 
    mixfrac  = NULL;
    Z        = NULL;
    C        = NULL;
    table    = NULL;

  }

  virtual ~CharlesChemtablePFPVCart2D() {

    cout << "~CharlesChemtablePFPVCart2D()" << endl;

    if ( mixfrac != NULL ) delete[] mixfrac;
    if ( Z       != NULL ) delete[] Z;
    if ( C       != NULL ) delete[] C;
    if ( table   != NULL ) delete   table; 

    for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f)
      if ((*f) != NULL ) delete (*f);

  }

  void init() {

    cout << " init()" << endl;

    // =======================================
    // generate flamelets
    // =======================================
    if (!checkParam("CHEMTABLE.FLAMELET_LIST"))
      createPremixedFlamelets();

    //========================================
    // get grid info
    //========================================

    nz = getIntParam("CART_CHEMTABLE.nZ"); 
    Z  = new double[nz];

    nc = getIntParam("CART_CHEMTABLE.nC");  
    C  = new double[nc];

    // read files and sort by phi
    readAndSortFlamelets();

    //========================================
    // initialize variables
    //========================================

    assert(tableVarsNameVec.size() == 0);

    // 2D charles premixed table requires the following variables
    tableVarsNameVec.push_back("rho");             // density
    tableVarsNameVec.push_back("T");               // temperature
    tableVarsNameVec.push_back("R");               // gas constant
    tableVarsNameVec.push_back("e");               // internal energy
    tableVarsNameVec.push_back("prog");            // progress variable
    tableVarsNameVec.push_back("src_prog");        // progress variable source term
    tableVarsNameVec.push_back("gamma");           // heat capacity ratio
    tableVarsNameVec.push_back("a_gamma");         // heat capacity coefficient
    tableVarsNameVec.push_back("mu");              // viscosity
    tableVarsNameVec.push_back("a_mu");            // viscosity coefficient
    tableVarsNameVec.push_back("locp");            // diffusivity
    tableVarsNameVec.push_back("a_locp");          // diffusivity coefficient
    tableVarsNameVec.push_back("sL");              // laminar burning velocity
    tableVarsNameVec.push_back("lF");              // laminar flame thickness
    tableVarsNameVec.push_back("mw");              // molecular weight
    tableVarsNameVec.push_back("int_rho_src");     // cumulative integral of (rho * csrc) dx
    tableVarsNameVec.push_back("int_heatrelease"); // cumulative integral of heat release (qdot * dx)
    tableVarsNameVec.push_back("s");               // entropy
    
    for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f) {
      (*f)->addLaminarFlameThickness();
      (*f)->addLaminarFlameSpeed();
    }

    if (checkParam("NOX_MODEL")) {
      tableVarsNameVec.push_back("src_nox_therm");
      tableVarsNameVec.push_back("int_src_nox_prompt");
      for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f)
        (*f)->addNOxSrcSplit();
    }

    // add thermo vars and coefficients
    perturbEnergy();

    // and add the user defined species
    // options:
    // Y_SPECIES    : mass fraction of SPECIES
    // X_SPECIES    : mole fraction of SPECIES
    // X_SPECIES_DRY: dry mole fraction of SPECIES
    // SPECIES      : mass, mole, and dry mole fractions of SPECIES
    if (Param * param = getParam("CHEMTABLE.SPECIES")) { 
      for (int i = 0; i<param->size(); ++i) {
        string var_name = param->getString(i);

        if (var_name[0] == 'X' || var_name[0] == 'Y') {
          tableVarsNameVec.push_back(var_name);

          if (var_name[0] == 'X') {
            for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f) {
              if ( var_name.substr(var_name.length()-4) == "_DRY" || var_name.substr(var_name.length()-4) == "_dry" )
                (*f)->addDryMoleFraction(var_name.substr(2,var_name.length()-6));
              else
                (*f)->addMoleFraction(var_name.substr(2,var_name.length()-2));
            }
          }//"X"
        }//(var_name[0] == "X" || "Y")
        else {
          tableVarsNameVec.push_back("Y_"+var_name);         // add mass fraction
          tableVarsNameVec.push_back("X_"+var_name);         // add mole fraction
          if (var_name != "H2O") tableVarsNameVec.push_back("X_"+var_name+"_DRY");  // add dry mole fraction

          for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f) {
            (*f)->addMoleFraction(var_name);
            if (var_name != "H2O") (*f)->addDryMoleFraction(var_name);
          }
        }//(var_name[0] != "X" || "Y")
      }//(i < param->size())
    }//"CHEMTABLE.SPECIES"

    nvars = tableVarsNameVec.size();
    assert(nvars>=15);

    // report min/max for each variable
    cout << " > variables to be included in the table: " << endl;
    for (int i=0; i<nvars; i++){ 
      double vmin =  1.0e+20; 
      double vmax = -1.0e+20; 
      for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f){
        vmin=min(vmin,(*f)->getVarMin(tableVarsNameVec[i])); 
        vmax=max(vmax,(*f)->getVarMax(tableVarsNameVec[i])); 
      }
      cout << "   > " << tableVarsNameVec[i] << ": min = " << vmin << ", max = " << vmax << endl;
    }

  }
  
  // =========================
  // build the data for table
  // =========================
  void build() {

    cout << " build()" << endl; 
    buildZgrid();
    buildCgrid();
    interpData(); 
  }
  
  // ===================================
  // write the table into a binary file
  // and dump some info on screen
  // ===================================
  void finalize(){
    
    cout << " finalize()" << endl;

    // write the table to a binary file
    writeBinary(filename);
    
    if (checkParam("CHEMTABLE.TECPLOT"))
      writeTecplotFile();
    
  }

private:

  void readAndSortFlamelets() { 

    assert(flameletList.size()==0);

    string filename = getStringParam("CHEMTABLE.FLAMELET_LIST", "FlameletList.txt");
    cout << " > reading flamelet list file: " << filename << endl;

    // open the file to get the size
    ifstream ifp;
    // open the file
    ifp.open(filename.c_str());
    if ( ifp.fail() )
      CERR_S("could not open flamelet file list: " << filename );
    string buffer_str;
    int nfl = 0;
    while(getline(ifp, buffer_str)) {
      if (buffer_str.empty())
        break;
      ++nfl;
    }
    ifp.close();
    cout << "   found " << nfl << " flamelet files" << endl;

    if ( nfl < 2 )
      CERR_S("CharlesPFPVCart2D requires more than one flamelet file");

    cout << " > reading flamelet files ..." << endl;

    // open the file again
    ifp.open(filename.c_str());
    if ( ifp.fail() )
      CERR_S("could not open flamelet file list: " << filename );

    int ifl = 0;
    int interval = max(int(0.05*nfl),1);
    while(getline(ifp, buffer_str)) {
      if ( (ifl+1) % interval == 0 )
        cout << "\r   progress ... " << (ifl+1)*5/interval << "%" << flush; 
      if (buffer_str.empty())
        break;
      istringstream buf(buffer_str);
      string thisfile;
      buf >> thisfile;

      // read the flamelet file
      PremixedFlamelet * flamelet = new PremixedFlamelet();
      flamelet->init(thisfile);
      flameletList.push_back(flamelet); 
      ++ifl;
    }
    ifp.close();
    cout << endl;
    assert(ifl == nfl);

    // sort based on phi
    flameletList.sort(isPhiSmallerFlamelet<PremixedFlamelet>); 

    // dump flamelet list if requested
    if ( checkParam("CHEMTABLE.DUMP_FLAMELETLIST") ) {
      cout << " > flamelet list ..." << endl;
      for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f) 
        cout << (*f)->filename << " ";
      cout << endl;
    }

    // dump some info on screen
    if (cti_verbose) {
      list<PremixedFlamelet*>::iterator f=flameletList.begin();
      (*f)->info();
    }
    
    {
      // check that the pressure is the same for all flamelets and store it
      list<PremixedFlamelet*>::iterator f=flameletList.begin();
      pressure = (*f)->pressure;
      for (f=++(flameletList.begin()); f!=flameletList.end(); ++f) {
        assert(pressure == (*f)->pressure );
      }
    }

  }

  void perturbEnergy() {

    const double deltaRhoE = 5000.0;
    Mixture * mixture = NULL; 

    cout << " > perturbing internal energy by deltaRhoE = " << deltaRhoE << endl;

    int count = 0;
    int step = max(int(0.05*flameletList.size()),1); 

    for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f) {

      (*f)->addMixtureERandPerturbedVariables(deltaRhoE);

      double *gamma   = (*f)->getVarPtr("gamma");
      double *a_gamma = (*f)->addVarAndReturnDataPtr("a_gamma");
      double *a_locp  = (*f)->addVarAndReturnDataPtr("a_locp");
      double *a_mu    = (*f)->addVarAndReturnDataPtr("a_mu");

      for (int i = 0; i<(*f)->npoints; ++i) {

        // get baseline vars
        const double R     = (*f)->getVarPtr("R")[i]; 
        const double T     = (*f)->getVarPtr("RT")[i]/R;
        const double e     = (*f)->getVarPtr("e")[i]; 
        const double mu    = (*f)->getVarPtr("mu")[i]; 
        const double locp  = (*f)->getVarPtr("locp")[i]; 

        // get +/- perturbations
        const double Tp    = (*f)->getVarPtr("RTp")[i]/R; 
        const double MUp   = (*f)->getVarPtr("MUp")[i];
        const double LOCPp = (*f)->getVarPtr("LOCPp")[i];

        const double Tm    = (*f)->getVarPtr("RTm")[i]/R; 
        const double MUm   = (*f)->getVarPtr("MUm")[i];
        const double LOCPm = (*f)->getVarPtr("LOCPm")[i];
	
        // derivatives: de/dT and d^2e/dT^2 at T=T
        const double deltaE  = deltaRhoE/(*f)->getVarPtr("rho")[i];
        const double dTm     = T  - Tm;
        const double dTp     = Tp - T;
        const double dT      = Tp - Tm;
        const double dedT    = (dTm*dTm*(e+deltaE) + (dTp-dTm)*dT*e - dTp*dTp*(e-deltaE))/(dTp*dTm*dT);
        const double d2edT2  = 2.0*(dTm*(e+deltaE) - dT*e + dTp*(e-deltaE))/(dTp*dTm*dT);
        const double dTm_log = log(T/Tm);
        const double dTp_log = log(Tp/T);
        const double dT_log  = log(Tp/Tm);
        
        // compute coefficients
        gamma[i]   = R / dedT + 1.0;
        a_gamma[i] = - d2edT2 * (gamma[i]-1.0) * (gamma[i]-1.0) / R;
        a_locp[i]  = (dTm_log*dTm_log*log(LOCPp) + (dTp_log-dTm_log)*dT_log*log(locp) - dTp_log*dTp_log*log(LOCPm))/(dTp_log*dTm_log*dT_log);
        a_mu[i]    = (dTm_log*dTm_log*log(MUp)   + (dTp_log-dTm_log)*dT_log*log(mu)   - dTp_log*dTp_log*log(MUm)  )/(dTp_log*dTm_log*dT_log);

      }

      ++count;
      if (count % step == 0)
        cout << "\r   progress ... " << count*5/step << "% " << flush; 
    }
    cout << endl; 
  }
  
  void buildZgrid() { 

    cout << " > building Z grid ..." << endl;

    //--------------------------------------------
    // construct mixture fraction from flamelets
    //--------------------------------------------

    int nf = flameletList.size(); 
    mixfrac = new double[nf];

    // Mixture fraction is computed from equivalence ratio and Zst
    // Both of these values are stored in premixed flamelets
    int iz = 0;
    for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f) {
      double FlameletZst = (*f)->Zst;
      mixfrac[iz] = 1.0 / ( 1.0 + (1.0-FlameletZst)/FlameletZst / ((*f)->phi+1.0e-20) );
      ++iz; 
    } 
    assert(iz == nf);

    //--------------------------------------------
    // build mixture fraction grid for table
    //--------------------------------------------

    double zmin = mixfrac[0];
    double zmax = mixfrac[flameletList.size()-1];

    if (Param * param = getParam("CART_CHEMTABLE.STRETCH_Z")) {
      int np=4;
      double *xp = new double[np];
      xp[0] = zmin;
      xp[1] = max(param->getDouble(0), zmin + 1.0E-8);
      xp[2] = min(param->getDouble(1), zmax - 1.0E-8);
      xp[3] = zmax;
      double frac  = param->getDouble(2);

      assert( (xp[0] < xp[1]) &&
              (xp[1] < xp[2]) &&
              (xp[2] < xp[3]) &&
              (frac  >  0.0 ) &&
              (frac  <  1.0 ) );

      double *qp = new double[np];
      qp[0] = 0.0;
      qp[1] = (1.0-frac)*(xp[1]-xp[0])/(xp[1]-xp[0] + xp[3]-xp[2]);
      qp[2] = qp[1]+frac;
      qp[3] = 1.0;

      stretchgrid(xp, qp, np, Z, nz);
      delete [] xp;
      delete [] qp;
    }
    else // uniform grid
      for (int i=0; i<nz; i++)
        Z[i] = zmin + (zmax-zmin)*double(i)/double(nz-1);

    //for (int i=0; i<nz; i++) cout << Z[i] << endl; 
  }

  void buildCgrid() {
  
    cout << " > building C grid ..." << endl;

    //--------------------------------------------
    // build progress variable grid for table
    //--------------------------------------------

    // NOTE: C is normalized between [0,1] in this table!!!
  
    if (Param * param = getParam("CART_CHEMTABLE.STRETCH_C")) { 
      int np=4;
      double *xp = new double[np];
      xp[0] = 0.0;
      xp[1] = max(param->getDouble(0), 0.0 + 1.0E-8);
      xp[2] = min(param->getDouble(1), 1.0 - 1.0E-8); 
      xp[3] = 1.0;
      double frac  = param->getDouble(2);

      assert( (xp[0] < xp[1]) &&
              (xp[1] < xp[2]) &&
              (xp[2] < xp[3]) &&
              (frac  >  0.0 ) &&
              (frac  <  1.0 ) );

      double *qp = new double[np];
      qp[0] = 0.0;
      qp[1] = (1.0-frac)*(xp[1]-xp[0])/(xp[1]-xp[0] + xp[3]-xp[2]);
      qp[2] = qp[1]+frac;
      qp[3] = 1.0;

      stretchgrid(xp, qp, np, C, nc);
      delete [] xp;
      delete [] qp;
    }//"CART_CHEMTABLE.STRETCH_C"
    else // uniform grid
      for (int i=0; i<nc; i++)
        C[i] = double(i)/double(nc-1);

    //for (int i=0; i<nc; i++) cout << C[i] << endl; 
  }

  void interpData() { 

    cout << " > interpolating data ..." << endl;

    assert(table == NULL);
    table = new CTIarray<double>(nz,nc,nvars);

    int isrc = -1; 
    int isL  = -1;
    int ilF  = -1;
    for (int ivar=0; ivar<nvars; ++ivar)
      if      (tableVarsNameVec[ivar].compare("src_prog") == 0) isrc=ivar;
      else if (tableVarsNameVec[ivar].compare("sL") == 0) isL=ivar;
      else if (tableVarsNameVec[ivar].compare("lF") == 0) ilF=ivar;
    assert(isrc != -1);
    assert(isL  != -1);
    assert(ilF  != -1);

    const double EPS = 1.0e-15; 
    const int nf = flameletList.size(); 

    // interpolate each flamelet to normalized progvar coordinate
      
    int i = 0;
    CTIarray<double> *data = new CTIarray<double>(nf,nc,nvars);
    for (list<PremixedFlamelet*>::iterator f=flameletList.begin(); f!=flameletList.end(); ++f) {

      const int np = (*f)->npoints; 
      double *prog = (*f)->getVarPtr("prog"); 
      const double cmin = (*f)->Cmin;
      const double crng = (*f)->Cmax-cmin+EPS;
       
      for (int j=0; j<nc; ++j) { 

        int ic;
        for (ic=0; ic<np-1; ++ic)
          if ((prog[ic+1]-cmin)/crng > C[j]) break;
        ic = min(ic,np-2); 

        double wc = 0.0;
        double denom = prog[ic+1]-prog[ic];
        if (denom > EPS)
          wc = max(0.0,min(1.0,(prog[ic+1]-cmin-C[j]*crng)/denom));

        for (int ivar=0; ivar<nvars; ++ivar) { 
          double *var = (*f)->getVarPtr(tableVarsNameVec[ivar]); 
          if (ivar == isrc || ivar == isL || ivar == ilF) // use log scale
            (*data)(i,j,ivar) = log(max(EPS,var[ic  ]))*wc
                                + log(max(EPS,var[ic+1]))*(1.0-wc); 
          else
            (*data)(i,j,ivar) = var[ic]*wc + var[ic+1]*(1.0-wc);
        }

      }
      ++i;
    }

    // interpolate to mixture fraction coordinate

    for (i=0; i<nz; ++i) {

      int iz; 
      for (iz=0; iz<nf-1; ++iz) 
        if (mixfrac[iz+1] > Z[i]) break;
      iz = min(iz,nf-2);
      double wz = max(0.0,min(1.0,(mixfrac[iz+1]-Z[i])/(mixfrac[iz+1]-mixfrac[iz]))); 

      for (int j=0; j<nc; ++j) { 
        for (int ivar=0; ivar<nvars; ++ivar) { 
          (*table)(i,j,ivar) = (*data)(iz,j,ivar)*wz + (*data)(iz+1,j,ivar)*(1.0-wz);

          if (ivar == isrc || ivar == isL || ivar == ilF)
            (*table)(i,j,ivar) = exp((*table)(i,j,ivar)); 
        }
      }
      (*table)(i,nc-1,isrc) = 0.0; // zero source term at max progvar
    }
    delete data; 
  }

  // =================================
  // write the table to a tecplot file
  // =================================
  void writeTecplotFile() {

    string filename = tablename + ".dat";

    cout << " > writing table to tecplot file: " << filename << endl;
   
    FILE * fp = fopen(filename.c_str(),"w");
    fprintf(fp,"title = \"chemtable\"\n");
    fprintf(fp,"variables = \"Z\"\n\"C\"\n");
    for (int ivar=0; ivar<tableVarsNameVec.size(); ++ivar)
      fprintf(fp,"\"%s\"\n",tableVarsNameVec[ivar].c_str());
    fprintf(fp,"zone i = %d, j = %d, DATAPACKING=POINT\n",nz,nc);
    for (int j=0; j<nc; ++j) {
      for (int i=0; i<nz; ++i) {
        fprintf(fp,"%g %g",Z[i],C[j]);
        for (int ivar=0; ivar<tableVarsNameVec.size(); ++ivar)
	  fprintf(fp," %g",(*table)(i,j,ivar));
        fprintf(fp,"\n");
      }
    }
    fclose(fp);
  }

protected:
  // ===============================
  // write a core IO compatible file
  // ===============================
  void writeBinary(const string& filename){
    
    cout << " > writing charles compatible table to the binary file: " << filename << endl;

    // open the file
    FILE * fp = NULL;
    fp=fopen(filename.c_str(),"wb");
    if (fp == NULL) 
      CERR_S(" cannot open chemtable file: " << filename);
    assert(fp);

    // use the header in the core IO ...

    // (1) write ugp magic number and version
    {
      int tmpint[2]={CHEM_IO_MAGIC_NUMBER,CHEM_IO_VERSION};
      size_t dummy= fwrite(tmpint,sizeof(int),2,fp);
      assert(dummy==2);
    }
    
    // (2) write the table type and dimensions
    {
      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      size_t dummy = tabletype.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and bye skip
      header.id = UGP_IO_CT_CART_2D;
      header.skip = sizeof(Header);
      // table dimensions
      header.idata[0] = nz; 
      header.idata[1] = nc; 
      header.idata[2] = nvars;
      // table pressure
      header.rdata[0] = pressure;
      // write into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
    }

    // (3) write "Z" coordinate
    {
      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
        header.name[i] = '\0';  
      string tmp    = "Z";
      size_t dummy = tmp.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and bye skip
      header.id = UGP_IO_CT_COOR;
      header.skip = nz * sizeof(double) + sizeof(Header);
      // table dimensions
      header.idata[0] = nz; 
      // write the into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
      // and the coordinate
      dummy= fwrite(Z,sizeof(double),nz,fp);
      assert(dummy==nz);
    }

    // (4) write "C" coordinate
    {
      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
        header.name[i] = '\0';  
      string tmp   = "C";
      size_t dummy = tmp.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and bye skip
      header.id = UGP_IO_CT_COOR;
      header.skip = nc * sizeof(double) + sizeof(Header);
      // table dimensions
      header.idata[0] = nc; 
      // write the into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
      // and the coordinate
      dummy= fwrite(C,sizeof(double),nc,fp);
      assert(dummy==nc);
    }

    // (5) write the data
    for(int ivar=0; ivar<nvars; ++ivar){

      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
        header.name[i] = '\0';  
      size_t dummy = tableVarsNameVec[ivar].copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and bye skip
      header.id = UGP_IO_CT_DATA;
      header.skip = nz*nc * sizeof(double) + sizeof(Header);
      // table dimensions
      header.idata[0] = nz*nc; 
      header.idata[1] = nz; 
      header.idata[2] = nc; 
      // write the header into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
      // and the data
      for (int iz=0; iz<nz; ++iz)
        for (int ic=0; ic<nc; ++ic) {
          double tmp = (*table)(iz,ic,ivar);
          dummy= fwrite(&tmp,sizeof(double),1,fp);
          assert(dummy==1);
        }
    }
        
    // end of file header
    {
      Header header;
      header.id = UGP_IO_EOF;
      sprintf(header.name,"EOF");
      header.skip = sizeof(Header);
      // write the header into the file      
      size_t dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);  
    }
    
    // close the file
    fclose(fp);
  
  }

}; //CharlesChemtablePFPVCart2D



// ================================================================================ //
//                           BEGIN KDT-SPECIFIC CODE                                //
// ================================================================================ //


// inheret to add some needed POD...
class MyFlamelet : public NonPremixedFlamelet {
public:
  double my_c;
  double w0;
  double w1;
  int this_point;  
};

struct dataStruct {
  string name;
  bool error_flag;
  
  dataStruct(const string& name) {
    this->name = name;
    error_flag = TRUE; // default to monitoring error 
  }
  
  dataStruct(const string& name, const bool error_flag) {
    this->name = name;
    this->error_flag = error_flag;
  }
};

// this is just for sorting flamelets based on progress variable interpolated to common Z
inline bool flameletCompareMyC(const MyFlamelet* flamelet0, const MyFlamelet* flamelet1) {
  return(flamelet1->my_c > flamelet0->my_c);
}

// container class for holding all of the flamelets
// also, filter flamelets in this class
class NonPremixedFlameletList {  

public:  
  list<MyFlamelet*> flameletList;
  list<MyFlamelet*> flameletList_copy; // we need a copy to hold old values
  
  int nflamelets;
  int np_total;

  double cmin, cmax;
  double pressure, Tfuel, Toxi;

  double rtol;
  
  NonPremixedFlameletList() {
    cout << "NonPremixedFlameletList()" << endl;
    
    readFlameletFilesAndSort(flameletList);
    readFlameletFilesAndSort(flameletList_copy);  // in lieu of writing copy constructor, we read in twice...
    
    this->nflamelets = flameletList.size();
    assert(nflamelets > 0);
    
    np_total = 0;
    for (list<MyFlamelet*>::iterator ifl = flameletList.begin(); ifl != flameletList.end(); ++ifl) {
      this->np_total += (*ifl)->npoints;
    }
  }

  ~NonPremixedFlameletList() {
    cout << "~NonPremixedFlameletList()" << endl;
    for (list<MyFlamelet*>::iterator flamelet = flameletList.begin(); flamelet!=flameletList.end(); ++flamelet)
      if ( (*flamelet) != NULL ) delete (*flamelet);

    for (list<MyFlamelet*>::iterator flamelet = flameletList_copy.begin(); flamelet!=flameletList_copy.end(); ++flamelet)
      if ( (*flamelet) != NULL ) delete (*flamelet);
  }
  
  // add these just for readability
  list<MyFlamelet*>::iterator begin() {
    return(flameletList.begin());
  }
  list<MyFlamelet*>::iterator end() {
    return(flameletList.end());
  }

  void addMixture(Mixture * mixture, const double rhoDeltaE) {
    cout << "  >> perturbing flamelets internal energy (rhoe) by  Delta(rhoe) = " << rhoDeltaE << endl;

    // have to add mixture to both original AND copy
    for (list<MyFlamelet*>::iterator ifl = flameletList.begin(); ifl != flameletList.end(); ++ifl) 
      (*ifl)->addMixtureERandPerturbedVariables(rhoDeltaE);
    for (list<MyFlamelet*>::iterator ifl_copy = flameletList_copy.begin(); ifl_copy != flameletList_copy.end(); ++ifl_copy) 
      (*ifl_copy)->addMixtureERandPerturbedVariables(rhoDeltaE);
    
  }

  void filter(const double sigma, const vector<dataStruct>& dataStructVec) {
    assert(nflamelets > 0);
    
    assert((sigma > 0) && (sigma < 1.0));    
    COUT2_S(" > filtering flamelets in progress variable with Gaussian filter with sigma = " << sigma << "... ");
    cout << "  progress:  ";
    cout.flush();
        
    int nvars = (*flameletList.begin())->varNameVec.size();

    // get min and max to compute relative error...
    double * var_min  = new double[nvars];
    double * var_max  = new double[nvars];
    int * var_flag    = new int[nvars];
    double * var_linf = new double[nvars];

    for (int ivar = 0; ivar < nvars; ++ivar) {
      var_flag[ivar] = -1;
      var_linf[ivar] = 0.0;
      var_min[ivar] = 1e20;
      var_max[ivar] = -1e20;
      for (list<MyFlamelet*>::iterator ifl = flameletList.begin(); ifl != flameletList.end(); ++ifl) {
	var_min[ivar] = min(var_min[ivar], (*ifl)->getVarMin((*ifl)->varNameVec[ivar]));
	var_max[ivar] = max(var_max[ivar], (*ifl)->getVarMax((*ifl)->varNameVec[ivar]));
      }
      assert(var_min[ivar] <= var_max[ivar]);
    }
    
    // flag vars according to info in dataStructVec
    // this is needed so as to report filtering-induced errors only in the variables we're interested in...
    for (int idata = 0; idata < dataStructVec.size(); ++idata) {
      if (dataStructVec[idata].name != "prog") {
	bool found = false;
	for (int ivar = 0; ivar < nvars; ++ivar) {
	  if ((*flameletList.begin())->varNameVec[ivar] == dataStructVec[idata].name) {
	    assert(var_flag[ivar] == -1 && found == false);
	    var_flag[ivar] = idata;
	    found = true;
	  }
	}
	
	// not all variables will be in flamelets. chris variables, for example, are computed after the fact
	//if (found == false)
	//  CERR_S(" variable " << dataNameVec[idata] << " not found in flamelets " << endl);
      }
    }

    double * prog_flamelet = new double[nflamelets];
    double * filter_weights = new double[nflamelets];

    int tmpint = max((int)(0.1*nflamelets),1);
    int ifl_count = 0;
    for (list<MyFlamelet*>::iterator ifl = flameletList.begin(); ifl != flameletList.end(); ++ifl) {
      
      assert((*ifl)->varNameVec.size() == nvars);

      if ( (ifl_count+1)%tmpint == 0 ) {
	cout << (ifl_count+1)*10/tmpint << "% ... ";
	cout.flush();
      }
      ifl_count++;
      
      double * Z = (*ifl)->getVarPtr("Z");
      double * C = (*ifl)->getVarPtr("prog");

      for (int ipoint = 0; ipoint < (*ifl)->npoints; ++ipoint) {
	
	// get interpolation weights for all flamelets
	// these weights are used to interpolate onto common Z
	for (list<MyFlamelet*>::iterator ifl_copy = flameletList_copy.begin(); 
	     ifl_copy != flameletList_copy.end(); ++ifl_copy) {
	  double * Z2 = (*ifl_copy)->getVarPtr("Z");
	  int this_point = 0;
	  while (Z2[this_point] <= Z[ipoint] && this_point < ((*ifl_copy)->npoints-1)) ++this_point;
	  const double diff = max(Z2[this_point] - Z2[this_point-1], 1e-10);
	  const double w0 = (Z2[this_point] - Z[ipoint])/diff;
	  const double w1 = 1.0 - w0;
	  assert(w0 >= 0 && w0 <= 1.0 && w1 >= 0 && w1 <= 1.0);

	  double * C2 = (*ifl_copy)->getVarPtr("prog");
	  (*ifl_copy)->my_c = w0*C2[this_point-1] + w1*C2[this_point];	  
	  
	  (*ifl_copy)->w0         = w0;
	  (*ifl_copy)->w1         = w1;
	  (*ifl_copy)->this_point = this_point;
	}
	// sort based on progress variable at common Z 
	flameletList_copy.sort(flameletCompareMyC);
	
	{
	  int index = 0;
	  for (list<MyFlamelet*>::iterator ifl_copy = flameletList_copy.begin(); 
	       ifl_copy != flameletList_copy.end(); ++ifl_copy) {
	    
	    prog_flamelet[index] = (*ifl_copy)->my_c; // store these in array for easy access
	    filter_weights[index] = 0.0;
	    ++index;
	  }
	}
	
	const double mu = C[ipoint];
	for (int index = 0; index < (nflamelets-1); ++index) {
	  const double c0 = prog_flamelet[index];
	  const double c1 = prog_flamelet[index+1];
	  const double diff = c1 - c0;
	  
	  assert(diff >= 0.0);
	  if ( diff > 1e-9 ) {
	    // perform exact Gaussian filtering on linearly interpolated data
	    const double val1 = 0.5*(erf((c1-mu)/(sigma*sqrt(2.0))) - erf((c0-mu)/(sigma*sqrt(2.0))));
	    const double val2 = mu*val1 - sigma/sqrt(2.0*M_PI)*(exp(-(c1-mu)*(c1-mu)/(2.0*sigma*sigma)) 
								- exp(-(c0-mu)*(c0-mu)/(2.0*sigma*sigma)));
	    filter_weights[index]   +=  (c1*val1 - val2)/diff;
	    filter_weights[index+1] += -(c0*val1 - val2)/diff;
	  }
	}
	// ends are clipped and filtering is done on unbounded domain
	filter_weights[0] += 0.5*(1.0 + erf((prog_flamelet[0]-mu)/(sigma*sqrt(2.0))));
	filter_weights[nflamelets-1] += 1.0 - 
	  0.5*(1.0 + erf((prog_flamelet[nflamelets-1]-mu)/(sigma*sqrt(2.0))));
	
	{
	  // force poisitive weights and check convex combination
	  double w_sum = 0.0;
	  for (int index = 0; index < nflamelets; ++index) {	    
	    if (filter_weights[index] < 0.0) 
	      filter_weights[index] = 0.0;	
	    w_sum += filter_weights[index];
	  }
	  assert(fabs(1.0 - w_sum) < 1e-5);	
	}

	for (int ivar = 0; ivar < nvars; ++ivar) {
	  // don't filter Z and progress variable since they are independent...
	  if (((*ifl)->varNameVec[ivar].compare("Z") != 0) && ((*ifl)->varNameVec[ivar].compare("prog") != 0)) {
	    double * var = (*ifl)->varDataVec[ivar];
	    double tmp = var[ipoint];
	    var[ipoint] = 0.0;
	    int index = 0;
	    
	    // interpolate data at common Z using w0 and w1 computed before
	    // and then filter...
	    for (list<MyFlamelet*>::iterator ifl_copy = flameletList_copy.begin(); 
		 ifl_copy != flameletList_copy.end(); ++ifl_copy) {
	      const int this_point = (*ifl_copy)->this_point;
	      const double w0 = (*ifl_copy)->w0;
	      const double w1 = (*ifl_copy)->w1;
	      assert(w0 >= 0 && w0 <= 1.0 && w1 >= 0 && w1 <= 1.0);
	      
	      double * var2 = (*ifl_copy)->varDataVec[ivar];
	      var[ipoint] += filter_weights[index]*(w0*var2[this_point-1] + w1*var2[this_point]);
	      ++index;
	    }
	    
	    var_linf[ivar] = max(var_linf[ivar], fabs((tmp - var[ipoint])/(var_max[ivar] - var_min[ivar])));             
	    
	  }
	}
	
      } // for (int ipoint = 0; ipoint < (*ifl)->npoints; ++ipoint)
	     
    } // for (list<NonPremixedFlamelet*>::iterator ifl = flameletList.begin(); ifl != flameletList.end(); ++ifl) 
    
    delete[] prog_flamelet;
    delete[] filter_weights;
    
    cout << " done." << endl;
    
    if (cti_verbose) {
      cout << "==================== filtering error summary ===================="<< endl;
      cout << "   relative linf errors due to filtering in each variable "<< endl;
      for (int ivar = 0; ivar < nvars; ++ivar) {
	if (var_flag[ivar] != -1 && dataStructVec[var_flag[ivar]].error_flag) {
	  cout << "      ";
	  cout.width(12);      
	  cout << left << dataStructVec[var_flag[ivar]].name;
	  cout.width(12);
	  cout << " " << var_linf[ivar] << endl;
	}
      }
      cout << "================================================================="<< endl;
    }        

    delete[] var_min;
    delete[] var_max;
    delete[] var_flag;
  }
  
  void readFlameletFilesAndSort(list<MyFlamelet*>& my_list) {
    assert(my_list.size() == 0 );

    string filename = getStringParam("CHEMTABLE.FLAMELET_LIST");        
    cout << " > reading the flamelet list file: " << filename << endl;

    // open the file to get the size
    ifstream ifp;
    // open the file
    ifp.open(filename.c_str());
    if ( ifp.fail() )
      CERR_S("could not open flamelet file list: " << filename );
    string buffer_str;
    int nfl = 0;
    while(getline(ifp, buffer_str)) {
      if (buffer_str.empty())
	break;
      ++nfl;
    }
    ifp.close();
    cout << "  found " << nfl << " flamelet files" << endl;

    cout << " > reading the flamelet files ..." << endl;
    
    // open the file again
    ifp.open(filename.c_str());
    if ( ifp.fail() )
      CERR_S("could not open flamelet file list: " << filename );
    cout << "  progress:  ";
    cout.flush();
    int ifl = 0;
    int tmpint = max((int)(0.1*nfl),1);
    while(getline(ifp, buffer_str)) {
      if ( (ifl+1)%tmpint == 0 ) {
	cout << (ifl+1)*10/tmpint << "% ... ";
	cout.flush();
      }
      if (buffer_str.empty())
	break;
      istringstream buf(buffer_str);
      string thisfile;
      buf >> thisfile;
      // read the flamelet file
      MyFlamelet * flamelet = new MyFlamelet();
      flamelet->init(thisfile);
      my_list.push_back(flamelet);
      ++ifl;
    }
    ifp.close();
    cout << endl;
    // check
    assert(ifl == nfl);
  
    // and sort based on Cmax
    my_list.sort(isProgSmallerFlamelet<NonPremixedFlamelet>);

    // dump the list
    if ( checkParam("CHEMTABLE.DUMP_FLAMELETLIST") ) {
      cout << " > flamelet list before clean up ..." << endl;
      for (list<MyFlamelet*>::iterator flamelet=my_list.begin(); flamelet!=my_list.end(); ++flamelet) 
	cout << (*flamelet)->filename << " ";
      cout << endl;
    }

    // dump some info on screen
    if (cti_verbose) {
      list<MyFlamelet*>::iterator flamelet=my_list.begin();
      (*flamelet)->info();
    }
    
    {
      // check that the pressure,and temperature at the fuel and oxidizer are the same for all the flamelets and store that pressure
      list<MyFlamelet*>::iterator flamelet=my_list.begin();
      pressure = (*flamelet)->pressure;
      Tfuel    = (*flamelet)->Tfuel;
      Toxi     = (*flamelet)->Toxi;
      for (flamelet=++(my_list.begin()); flamelet!=my_list.end(); ++flamelet) {
	assert(pressure == (*flamelet)->pressure );
	assert(Tfuel == (*flamelet)->Tfuel );
	assert(Toxi == (*flamelet)->Toxi );
      }
    }
    
    // get the min and max values for C here
    cmin = 1.0e20;
    cmax = -1.0e20;
    for (list<MyFlamelet*>::iterator flamelet=my_list.begin(); flamelet!=my_list.end(); ++flamelet) {
      cmin = min(cmin,(*flamelet)->Cmin);
      cmax = max(cmax,(*flamelet)->Cmax);
    }
    
    if (fabs(cmin)<1.0e-14) cmin = 0.0;
    assert(cmin>=0.0);
    assert(cmax<=1.0);

  }

};

// PointCloud is unstructured data - the whole chemtable generation problem is essentially how to interpolate back and forth
// from unstructured data (flamelets) to structured data (the chemtable)
template <int DIM_K>
class PointCloud {

public:
  
  int data_m;
  int np;

  double (*x_array)[DIM_K];
  double *data_array;
  double *error_array;
  int * flag;
  
  double x_min[DIM_K];
  double x_max[DIM_K];
  double * data_min;
  double * data_max;
  
  vector<dataStruct> dataStructVec; // data names - depends on table type, user-specified species, etc..
  vector<string> xNameVec;    // coordinate names - fixed by table type...
  
  PointCloud() {
    cout << "PointCloud()" << endl;
    
    np     = -1;
    data_m = -1;
    data_min     = NULL;
    data_max     = NULL;
    x_array     = NULL;
    data_array  = NULL;
    error_array = NULL;
    flag        = NULL;   
  }
  
  virtual ~PointCloud() {
    if (data_min    != NULL)  delete[] data_min;
    if (data_max    != NULL)  delete[] data_max;
    if (x_array     != NULL)  delete[] x_array;
    if (data_array  != NULL)  delete[] data_array;
    if (error_array != NULL)  delete[] error_array;  
    if (flag        != NULL)  delete[] flag;
  }
  
  void init() {
    assert(np > 0);
    
    // derived class is required to fill the NameVec...
    assert(dataStructVec.size() > 0);
    data_m = dataStructVec.size();
    assert(xNameVec.size() == DIM_K);
    
    // ABC handles allocation... 
    assert(x_array     == NULL); x_array     = new double[np][DIM_K];
    assert(data_array  == NULL); data_array  = new double[np*data_m];
    assert(error_array == NULL); error_array = new double[np*data_m];   
    assert(flag        == NULL); flag        = new int[np];    

    assert(data_min == NULL); data_min = new double[data_m];
    assert(data_max == NULL); data_max = new double[data_m];   
    
    populatePC();
    calcXMinMax(np);
    calcDataMinMax(np);

    cout << "> number of points in cloud: " << np << endl;
  }
   
  void calcXMinMax(const int n) {
    assert(x_array != NULL);
    for (int idim = 0; idim < DIM_K; ++idim) {
      x_max[idim] = -1e20;
      x_min[idim] = 1e20;
      for (int ip = 0; ip < n; ++ip) {
	x_max[idim] = max(x_max[idim],x(ip,idim));
	x_min[idim] = min(x_min[idim],x(ip,idim));
      }
      assert(x_max[idim] > x_min[idim]);
    }
  }

  void calcDataMinMax(const int n) {
    assert(data_array != NULL && data_min != NULL && data_max != NULL);
    for (int idata = 0; idata < data_m; ++idata) {
      data_min[idata] = 1e20;
      data_max[idata] = -1e20;
      for (int ip = 0; ip < n; ++ip) {
	data_min[idata] = min(data_min[idata],data(ip,idata));
	data_max[idata] = max(data_max[idata],data(ip,idata));	
      }
      assert(data_max[idata] > data_min[idata]);
    }  
  }

  void info() {
    assert(xNameVec.size() == DIM_K);
    assert(dataStructVec.size() == data_m);
    
    extraInfo();
    cout << "\n";
    cout << "   coordinates (kdim = " << DIM_K << "):" << endl;
    for (int idim = 0; idim < DIM_K; ++idim) {
      cout << "      ";
      cout.width(12);      
      cout << left << xNameVec[idim];
      cout.width(12);
      cout << "   (min,max):  " << "(" << x_min[idim] << "," << x_max[idim] << ")" << endl;
    }
        
    cout << "   variables (nvars = " << data_m << "):" << endl;
    for (int idata = 0; idata < data_m; ++idata) {
      cout << "      ";
      cout.width(12);      
      cout << left << dataStructVec[idata].name;
      cout.width(12);
      cout << "   (min,max):  " << "(" << data_min[idata] << "," << data_max[idata] << ")" << endl;
    }
    
  }

  // allows us to dump derived class-specific info...    
  virtual void extraInfo() {}

  // allow derived classes to fill info in binary chemtable header...
  virtual void extraMainHeader(Header& header, const int ioffset, const int roffset) { }

  // various accesors...  
  inline double& x(const int ip, const int idim) {
#ifdef CTI_DEBUG
    assert(ip >= 0 && ip < np);
    assert(idim >= 0 && idim < DIM_K);
#endif
    return(x_array[ip][idim]);
  }

  inline double* x(const int ip) {
#ifdef CTI_DEBUG
    assert(ip >= 0 && ip < np);
#endif
    return(x_array[ip]);
  }

  inline double& data(const int ip, const int idata) {
#ifdef CTI_DEBUG
    assert(ip >= 0 && ip < np);
    assert(idata >= 0 && idata < data_m);
#endif
    return(data_array[ip*data_m + idata]);
  }

  inline double* data(const int ip) {
#ifdef CTI_DEBUG
    assert(ip >= 0 && ip < np);
#endif
    return(&data_array[ip*data_m]);
  }
  
  inline double& error(const int ip, const int idata) {
#ifdef CTI_DEBUG
    assert(ip >= 0 && ip < np);
    assert(idata >= 0 && idata < data_m);
#endif
    return(error_array[ip*data_m + idata]);
  }
  
  virtual int queryPC(double * dataq, const double xq[DIM_K]) = 0;

  virtual void populatePC() = 0;
    
};

class NonPremixed3DPC : public PointCloud<3> {
protected:
  
  int nzvar;
  int nclipped;
  
  NonPremixedFlameletList * flamelets;
  
public:
  
  NonPremixed3DPC() {
    cout << "NonPremixed3DPC()" << endl;    
    
    flamelets = new NonPremixedFlameletList;
    
    this->nzvar = getIntParam("KD_CHEMTABLE.nZvar", 100);
    this->nclipped = getIntParam("KD_CHEMTABLE.NCLIPPED", 20000);
    this->np = (flamelets->np_total)*nzvar + nclipped;
    
    assert(this->xNameVec.size() == 0);
    this->xNameVec.push_back("Z");
    this->xNameVec.push_back("Zvar_norm");
    this->xNameVec.push_back("C");
  
    // all nonpremixed3d chemtables can transport extra species...
    assert(dataStructVec.size() == 0);
    bool flag_species = checkParam("KD_CHEMTABLE.SPECIES_ERROR");
    //if (checkParam("CHEMTABLE.SPECIES")) 
    if (Param * param = getParam("CHEMTABLE.SPECIES")) { 
      for (int i = 0; i<param->size(); ++i) {
	      dataStruct specie("Y_" + param->getString(i), flag_species);
	      this->dataStructVec.push_back(specie);
      }//(i < param->size())
    }//"CHEMTABLE.SPECIES"

    // all non-premixed (both chris and vida) require these 6 variables...
    dataStruct rho("rho");            dataStructVec.push_back(rho);       //  (1) density,
    dataStruct T("T");                dataStructVec.push_back(T);         //  (2) temperature,
    dataStruct mu("mu");              dataStructVec.push_back(mu);        //  (3) viscosity,   
    dataStruct locp("locp");          dataStructVec.push_back(locp);      //  (4) diffusivity,  
    dataStruct prog("prog");          dataStructVec.push_back(prog);      //  (5) progress variable,   
    dataStruct src_prog("src_prog");  dataStructVec.push_back(src_prog);  //  (6) progress variable source term

    // more will be set in derived constructor...
  }
  
  virtual ~NonPremixed3DPC() {
    cout << "~NonPremixed3DPC()" << endl;
    
    delete flamelets;
  }
  
  void init() {
    
    // do filtering of flamelets here...
    double sigma = getDoubleParam("FLAMELET.FILTER_WIDTH", 5e-4);
    
    // turning off filtering is accomplished by setting sigma = 0
    if (sigma > 0.0) 
      flamelets->filter(sigma, dataStructVec);
    else 
      COUT2_S(" > flamelet filtering disabled");
    
    PointCloud<3>::init();
  }
  
  void populatePC() {
    COUT2_S("> populating point cloud... ");
    cout << " (flamelets) progress:  ";
    cout.flush();
        
    assert(x_array != NULL);
    assert(data_array != NULL);
    
    int count = 0;
    int tmpint = max((int)(0.1*(flamelets->np_total)*nzvar),1);
    for (list<MyFlamelet*>::iterator ifl = flamelets->begin(); ifl != flamelets->end(); ++ifl) {
      
      double * Z = (*ifl)->getVarPtr("Z");
      for (int ipoint = 0; ipoint < (*ifl)->npoints; ++ipoint) {
	for (int izvar = 0; izvar < nzvar; ++izvar) {
	  if ( (count+1)%tmpint == 0 ) {
	    cout << (count+1)*10/tmpint << "% ... ";
	    cout.flush();
	  }  
	  
	  x(count,0) = Z[ipoint];
	  x(count,1) = izvar*izvar/((double) (nzvar-1)*(nzvar-1));  // quadratic in normalized zvar
	  
	  const double this_zvar = x(count,1)*x(count,0)*(1.0 - x(count,0));
	  x(count,2) = (*ifl)->getValGivenMeanAndVar(x(count,0), this_zvar, "prog");
	  
	  getAllData(data(count), *ifl, x(count,0), this_zvar);
	  ++count;
	}
      }
    }
    cout << "done." << endl;
    
    assert(count == (flamelets->np_total)*nzvar);
    
    cout << " (clipped)  progress:  ";
    cout.flush();
    
    calcXMinMax(count);
    // add points in clipped region
    // let's keep it deterministic...
    srand(1);
    tmpint = max((int)(0.1*nclipped),1);
    while (count < np) {
      for (int idim = 0; idim < 3; ++idim) 
	x(count,idim) = (double) rand()/(double) RAND_MAX*(x_max[idim] - x_min[idim]);
      
      if (queryPC(data(count), x(count)) == 0) {
	// if we're clipped, increment count and move on...
	++count;
	
	if ( (count - (flamelets->np_total)*nzvar + 1)%tmpint == 0 ) {
	  cout << (count - (flamelets->np_total)*nzvar + 1)*10/tmpint << "% ... ";
	  cout.flush();
	}     
      }
    }
    COUT2_S("done.");
  
  }

  int queryPC(double * dataq, const double xq[3]) {
    
    double this_z = xq[0];
    double this_zvar = xq[1]*this_z*(1.0-this_z);

    const double this_c = xq[2];
    // try to bracket this_c with c0 and c1...
    MyFlamelet * ifl0 = NULL;
    MyFlamelet * ifl1 = NULL;
    double c0 = -1e20;
    double c1 = 1e20;
    for (list<MyFlamelet*>::iterator ifl = flamelets->begin(); ifl != flamelets->end(); ++ifl) {
      double cifl = (*ifl)->getValGivenMeanAndVar(this_z, this_zvar, "prog");
      if (cifl < this_c) {
	if ((this_c - cifl) <= (this_c - c0)) {
	  c0 = cifl;
	  ifl0 = (*ifl);
	}
      }
      else {
	if ((cifl - this_c) <= (c1 - this_c)) {
	  c1 = cifl;
	  ifl1 = (*ifl);
	}
      }
    }      
    
    // make sure we found at least one...
    assert((ifl0 != NULL) || (ifl1 != NULL));   
    
    if (ifl0 == NULL) {
      // all flamelets have c greater than given so clip...
      getAllData(dataq, ifl1, this_z, this_zvar);      
      return(0); // let calling code know that we clipped
    }
    else if (ifl1 == NULL) {
      // all flamelets have c less than given so clip...
      getAllData(dataq, ifl0, this_z, this_zvar);
      return(0); // let calling code know that we clipped
    }
    else {
      assert((c1-c0) >= 0.0);
      // linearly interpolate
      
      const double diff = max(c1 - c0, 1e-10);
      const double w0 = (c1 - this_c)/diff;
      const double w1 = 1.0 - w0;
      
      double * data0 = new double[data_m];
      double * data1 = new double[data_m];

      getAllData(data0, ifl0, this_z, this_zvar);
      getAllData(data1, ifl1, this_z, this_zvar);

      for (int idata = 0; idata < data_m; ++idata) 
	dataq[idata] = w0*data0[idata] + w1*data1[idata];

      delete[] data0;
      delete[] data1;
      return(1); // let calling code know that we interpolated
    }
  }

  void extraInfo() {   
    cout << "   reference pressure:       " << flamelets->pressure << endl;  
    cout << "   oxidizer temperature:     " << flamelets->Toxi << endl;
    cout << "   fuel temperature:         " << flamelets->Tfuel << endl;
  }

  void extraMainHeader(Header& header, const int ioffset, const int roffset) {
    header.rdata[roffset]     = flamelets->pressure;
    header.rdata[roffset + 1] = flamelets->Toxi;
    header.rdata[roffset + 2] = flamelets->Tfuel;
  }


  // this is virtual so that derived Chris can have different behavoir
  virtual void getAllData(double * this_data, MyFlamelet * ifl, const double this_z, const double this_zvar) {
    
    for (int idata = 0; idata < data_m; ++idata)
      this_data[idata] = ifl->getValGivenMeanAndVar(this_z, this_zvar, dataStructVec[idata].name); 
  }

};

class VidaNonPremixed3DPC : public NonPremixed3DPC {
  
public:
  VidaNonPremixed3DPC() {
    cout << "VidaNonPremixed3DPC()" << endl;    
    
    // a vida table requires at least 6 variables (already set in NonPremixed3DPC constructor)   
    
    NonPremixed3DPC::init(); 
  }

  ~VidaNonPremixed3DPC() {
    cout << "~VidaNonPremixed3DPC()" << endl;    
  }
  
};

class ChrisNonPremixed3DPC : public NonPremixed3DPC {

public:  
  double rhoDeltaE;
  Mixture * mixture;

  ChrisNonPremixed3DPC() {
    cout << "ChrisNonPremixed3DPC()" << endl;    
    
    // a chris table requires at least 12 variables
    // first 6 are same as vida (already set in NonPremixed3DPC constructor)...
    
    dataStruct R("R");              dataStructVec.push_back(R);       //  (7)  gas constant
    dataStruct e("e");              dataStructVec.push_back(e);       //  (8)  internal energy
    dataStruct gamma("gamma");      dataStructVec.push_back(gamma);   //  (9)  specific heat ratio  
    dataStruct a_gamma("a_gamma");  dataStructVec.push_back(a_gamma); //  (10) coefficient for specific heat ratio
    dataStruct a_mu("a_mu");        dataStructVec.push_back(a_mu);    //  (11) coefficent for the viscosity
    dataStruct a_locp("a_locp");    dataStructVec.push_back(a_locp);  //  (12) coefficient for diffusivity
           
    mixture = NULL;
    rhoDeltaE = 5000.0;  
    flamelets->addMixture(mixture, rhoDeltaE);    
    
    NonPremixed3DPC::init(); 
  }

  ~ChrisNonPremixed3DPC() {
    cout << "~ChrisNonPremixed3DPC()" << endl;    
    
    if (mixture != NULL) delete mixture;
  }
  
  void getAllData(double * this_data, MyFlamelet * ifl, const double this_z, const double this_zvar) {
    
    // get the variables that are already in the table 
    const double rho      = ifl->getValGivenMeanAndVar(this_z, this_zvar, "rho");
    const double deltaE   = rhoDeltaE/rho;
    const double R        = ifl->getValGivenMeanAndVar(this_z, this_zvar, "R");
    const double RT       = ifl->getValGivenMeanAndVar(this_z, this_zvar, "RT");
    const double T        = RT / R;
    const double e        = ifl->getValGivenMeanAndVar(this_z, this_zvar, "e");
    const double prog     = ifl->getValGivenMeanAndVar(this_z, this_zvar, "prog");
    const double src_prog = ifl->getValGivenMeanAndVar(this_z, this_zvar, "src_prog");
    const double mu       = ifl->getValGivenMeanAndVar(this_z, this_zvar, "mu");
    const double locp     = ifl->getValGivenMeanAndVar(this_z, this_zvar, "locp");
    
    // get the plus and minus perturbatiosn
    const double RTp    = ifl->getValGivenMeanAndVar(this_z, this_zvar, "RTp");
    const double Tp     = RTp / R;
    const double MUp    = ifl->getValGivenMeanAndVar(this_z, this_zvar, "MUp");
    const double LOCPp  = ifl->getValGivenMeanAndVar(this_z, this_zvar, "LOCPp");
    
    const double RTm    = ifl->getValGivenMeanAndVar(this_z, this_zvar, "RTm");
    const double Tm     = RTm / R;
    const double MUm    = ifl->getValGivenMeanAndVar(this_z, this_zvar, "MUm");
    const double LOCPm  = ifl->getValGivenMeanAndVar(this_z, this_zvar, "LOCPm");

    // derivatives: de/dT and d^2e/dT^2 at T=T
    const double dTm    = T  - Tm;
    const double dTp    = Tp - T;
    const double dT     = Tp - Tm;
    const double dedT   = ( dTm*dTm*(e + deltaE) + (dTp-dTm)*dT*e - dTp*dTp*(e - deltaE) ) / (dTp*dTm*dT);
    const double d2edT2 = 2.0 * ( dTm*(e + deltaE) - dT*e + dTp*(e - deltaE) ) / (dTp*dTm*dT);
    
    // compute coefficients
    const double gamma   = R / dedT + 1.0;
    const double a_gamma = - d2edT2 * (gamma-1.0) * (gamma-1.0) / R;

    const double dTm_log = log(T/Tm);
    const double dTp_log = log(Tp/T);
    const double dT_log = log(Tp/Tm);

    /* DAP: corrected tabulation of a_mu, a_locp */
    const double a_mu   = ( dTm_log*dTm_log*log(MUp)   + (dTp_log-dTm_log)*dT_log*log(mu)   - dTp_log*dTp_log*log(MUm)   ) / (dTp_log*dTm_log*dT_log);
    const double a_locp = ( dTm_log*dTm_log*log(LOCPp)   + (dTp_log-dTm_log)*dT_log*log(locp)   - dTp_log*dTp_log*log(LOCPm)   ) / (dTp_log*dTm_log*dT_log);
/*
    // this is adapted from Vincent's code. I think it is wrong
    // it should take the log of the MUp and ...
    const double a_mu    = ( dTm*dTm*MUp   + (dTp-dTm)*dT*mu   - dTp*dTp*MUm   ) / (dTp*dTm*dT);
    const double a_locp  = ( dTm*dTm*LOCPp + (dTp-dTm)*dT*locp - dTp*dTp*LOCPm ) / (dTp*dTm*dT);
*/
            
    for (int idata = 0; idata < data_m; ++idata) {
      if      (dataStructVec[idata].name == "rho"      ) this_data[idata] = rho;
      else if (dataStructVec[idata].name == "T"        ) this_data[idata] = T;
      else if (dataStructVec[idata].name == "R"        ) this_data[idata] = R;
      else if (dataStructVec[idata].name == "e"        ) this_data[idata] = e;
      else if (dataStructVec[idata].name == "prog"     ) this_data[idata] = prog;
      else if (dataStructVec[idata].name == "src_prog" ) this_data[idata] = src_prog;
      else if (dataStructVec[idata].name == "gamma"    ) this_data[idata] = gamma;
      else if (dataStructVec[idata].name == "a_gamma"  ) this_data[idata] = a_gamma;
      else if (dataStructVec[idata].name == "mu"       ) this_data[idata] = mu;
      else if (dataStructVec[idata].name == "a_mu"     ) this_data[idata] = a_mu;
      else if (dataStructVec[idata].name == "locp"     ) this_data[idata] = locp;
      else if (dataStructVec[idata].name == "a_locp"   ) this_data[idata] = a_locp;
      else 
	this_data[idata] = ifl->getValGivenMeanAndVar(this_z, this_zvar, dataStructVec[idata].name); 
    }
  }

};

  
// space is discretized with RANGE indicating the max
#define RANGE 2147483648 // 2^31
template <int DIM_K>
class KdTable : public AbstractChemtable {

private:

  struct Vertex {
    int8 coords[DIM_K];   // index in discretized space
    double x[DIM_K];
    double * data;
    int ive;
    Vertex * left;         // verticies are stored in binary tree for fast comparison (may not be worth the complexity)
    Vertex * right;
    
    Vertex(const int data_m) {
      assert(data_m > 0);
      data = new double[data_m];
      for (int idim = 0; idim < DIM_K; ++idim) coords[idim] = -1;
      left  = NULL;
      right = NULL;
    }
    
    ~Vertex() {
      delete[] data;
    }
    
  };
  
  // using different vertex numbering to take advantage of bit shifts -
  // in 3D...  
  //      (6)      (7)
  //      x-------x
  //     /|      /|
  // (2)/ |  (3)/ |
  //   x-------x  |
  //   |  x----|--x 
  //   | / (4) | / (5)
  //   |/      |/
  //   x-------x
  // (0)        (1)
  
  struct HyperCube {
    int ve[1<<DIM_K];  // 2^k corners, ordered according to bits, i.e. in 3D:
    // (0,0,0) = 0 is minimum corner and (1,1,1) = 7 is maximum  
    int ve_split[DIM_K][1<<(DIM_K-1)];  // k*(2^(k-1)) edges split for candidates
    int depth;
    int index;
    list<int> ipList;
   
  };
  
  // nodes corresponding to branches in kd-tree
  struct Node {
    int dim_split;
    double x_split;
    Node *left;
    Node *right;
    HyperCube * hc;
    int index;
    
    Node() {
      dim_split = -1;
      x_split = -123.4567;  // to silence valgrind and catch errors 
      left  = NULL;
      right = NULL;
      hc    = NULL;
    }
  }; 
  
  double dx[DIM_K]; // space is discretized, so dx is minimum unit of distance in each dimension 
  
  vector<Vertex *> vecPtrVe;   // redunant information, but helpful for I/O
  Vertex * root_vertex;

  vector<HyperCube *> vecPtrHC;  // redundant information, but helps with I/O 
  
  Node * root_node;
  vector<Node *> vecPtrNode;  // redundant information, but helps with I/O 

  double rtol;
  int max_nhc;
  int max_depth;
  
public:

  PointCloud<DIM_K> * pc;

  KdTable() {
    cout << "KdTable()" << endl;
        
    // RANGE must be power of 2 because we rely on bisecting the set...
    assert(!((RANGE)&(RANGE-1)));
    
    pc = NULL;
  
    // grab parameters
    rtol = getDoubleParam("KD_CHEMTABLE.RTOL", 1e-2);
    max_nhc = getIntParam("KD_CHEMTABLE.MAX_LEAVES", 160000);
    
    {
      int max_depth_def = 0;
      uint8 tmp = RANGE;
      while (tmp > 1) {
	tmp >>= 1;
	++max_depth_def;
      }
      max_depth = getIntParam("KD_CHEMTABLE.MAX_DEPTH", 8*DIM_K);
     
      if (max_depth > max_depth_def) 
	CERR_S(" max depth of " << max_depth << " exceeds integer representation");
    }

  }
  
  ~KdTable() {
    if (pc != NULL) delete pc;

    // manual cleanup is necessary....
    for (int ihc = 0; ihc < vecPtrHC.size(); ++ihc) {
      assert(vecPtrHC[ihc] != NULL);
      delete vecPtrHC[ihc];
    }
    for (int ive = 0; ive < vecPtrVe.size(); ++ive) {
      assert(vecPtrVe[ive] != NULL);
      delete vecPtrVe[ive];
    }
    for (int inode = 0; inode < vecPtrNode.size(); ++inode) {
      assert(vecPtrNode[inode] != NULL);
      delete vecPtrNode[inode];
    }
    
  }

  void init() {  
    cout << "KdTable::init()" << endl;

    // make sure we have allocated a point cloud
    assert(pc != NULL);

    for (int idim = 0; idim < DIM_K; ++idim) {
      assert(pc->x_min[idim] < pc->x_max[idim]);
      dx[idim] = (pc->x_max[idim] - pc->x_min[idim])/((double) RANGE);
    }
    
    // initialize first HyperCube
    HyperCube * first_hc = new HyperCube;
    
    // manually generate root vertex...
    root_vertex = new Vertex(pc->data_m);
    root_vertex->ive = 0;
    for (int idim = 0; idim < DIM_K; ++idim) {
      root_vertex->coords[idim] = 0;
      root_vertex->x[idim]      = pc->x_min[idim];
    }
    pc->queryPC(root_vertex->data, root_vertex->x);      
    vecPtrVe.push_back(root_vertex);
    first_hc->ve[0] = 0; 
   
    for (int ive = 1; ive < (1<<DIM_K); ++ive) {
      Vertex * new_vertex = new Vertex(pc->data_m); 
      for (int idim = 0; idim < DIM_K; ++idim) {
	if (ive&(1<<idim))	  
	  new_vertex->coords[idim] = RANGE;	
	else
	  new_vertex->coords[idim] = 0;
      }
      addVertex(new_vertex);
      first_hc->ve[ive] = new_vertex->ive;
    }
    initCandidates(first_hc);
    vecPtrHC.push_back(first_hc);
    first_hc->index = 0;
    
    // initialize root node
    root_node = new Node;
    
    root_node->hc = first_hc;
    first_hc->depth = 0;
    root_node->index = 0;
    vecPtrNode.push_back(root_node);

    // first hypercube has all of the point cloud
    for (int ip = 0; ip < pc->np; ++ip) 
      first_hc->ipList.push_back(ip);
    
    // compute error for the whole PC...
    double * dataq = new double[pc->data_m];
    for (int ip = 0; ip < pc->np; ++ip) {
      linearInterpHC(dataq, pc->x(ip), first_hc);
      
      for (int idata = 0; idata < pc->data_m; ++idata) 
	pc->error(ip,idata) = fabs((pc->data(ip,idata) - dataq[idata])/(pc->data_max[idata] - pc->data_min[idata]));
    }
    delete[] dataq;

  }

  void build() {
    cout << "KdTable::build()" << endl;
    
    list<Node*> active_leaves;
    active_leaves.push_back(root_node);
    
    double * dataq = new double[pc->data_m];
    
    // for reporting clipped leaves...
    int count_zeroed = 0;
    double max_error_zeroed = 0.0;

    int iter = 0;
    while (active_leaves.size() > 0) {
      // pop the first one
      Node * this_leaf = *(active_leaves.begin());
      HyperCube * this_hc = this_leaf->hc;
      assert(this_hc != NULL);

      // for all points inside the hypercube, find maximum error...
      double max_error = 0.0;
      int max_error_idata, max_error_ip;
      for (list<int>::iterator ilist = this_hc->ipList.begin(); ilist != this_hc->ipList.end(); ++ilist) {
	const int ip = *ilist;
#ifdef CTI_DEBUG
	//assert(pointIsInsideLeaf(this_leaf, pc->x(ip)));   
#endif
	for (int idata = 0; idata < pc->data_m; ++idata) {
	  if (pc->dataStructVec[idata].error_flag) { // check if this variable has been flagged for error monitoring
	    const double ip_error = pc->error(ip,idata);
	    if (ip_error > max_error) {
	      max_error_ip    = ip;
	      max_error_idata = idata;
	      max_error       = ip_error;
	    }
	  }
	}
      }
      
      // decide which dimension to split based on which causes the most change in the error
      // NB: if we try to force a monotonic decrease, pathological cases can ensue...
      double diff = 0.0;
      int dim_split = 0;
      for (int idim = 0; idim < DIM_K; ++idim) {
	queryHC_split(dataq, pc->x(max_error_ip), this_leaf, idim);
	double this_error = fabs((pc->data(max_error_ip,max_error_idata) - dataq[max_error_idata])/(pc->data_max[max_error_idata] - pc->data_min[max_error_idata]));
	double my_diff = fabs(this_error - max_error);
	if (my_diff > diff) {
	  diff = my_diff;
	  dim_split = idim;
	}
      }
      
      // and split...
      Node * left, * right;
      splitLeaf(left, right, this_leaf, dim_split);
      // left is old, right is new...
      assert(left->hc == this_hc);
      HyperCube * new_hc = right->hc;
      
      // erase this_leaf from list because no longer a leaf...
      active_leaves.erase(active_leaves.begin());
      
      // redistribute point list...
      assert(new_hc->ipList.size() == 0);
      {
	list<int>::iterator ilist = this_hc->ipList.begin();
	while (ilist != this_hc->ipList.end()) {
	  const int ip = *ilist;
	  if (pointIsInsideLeaf(right, pc->x(ip))) {
	    new_hc->ipList.push_back(ip);
	    this_hc->ipList.erase(ilist++);
	  }
	  else {
	    //assert(pointIsInsideLeaf(left, pc->x(ip)));
	    ++ilist;
	  }
	}
      }
      
      // recompute error...
      max_error = 0.0;
      for (list<int>::iterator ilist = this_hc->ipList.begin(); ilist != this_hc->ipList.end(); ++ilist) {
	const int ip = *ilist;
#ifdef CTI_DEBUG 
	assert(pointIsInsideLeaf(left, pc->x(ip))); 
#endif
	linearInterpHC(dataq, pc->x(ip), this_hc);
	for (int idata = 0; idata < pc->data_m; ++idata) {
	  if (pc->dataStructVec[idata].error_flag) { 
	    pc->error(ip,idata) = fabs((pc->data(ip,idata) - dataq[idata])/(pc->data_max[idata] - pc->data_min[idata]));
	    max_error = max(max_error, pc->error(ip,idata));
	  }
	}
      }
            
      // if we haven't met tolerance, add to list...      
      if (max_error >= rtol) {
	if (this_hc->depth < max_depth) {
	  active_leaves.push_back(left);
	}
	else {
	  for (list<int>::iterator ilist = this_hc->ipList.begin(); ilist != this_hc->ipList.end(); ++ilist) 
	    ++count_zeroed;
	  max_error_zeroed = max(max_error_zeroed, max_error);
	} 
      }     

      max_error = 0.0;
      for (list<int>::iterator ilist = new_hc->ipList.begin(); ilist != new_hc->ipList.end(); ++ilist) {
	const int ip = *ilist;
#ifdef CTI_DEBUG 
	assert(pointIsInsideLeaf(right, pc->x(ip)));   
#endif
	linearInterpHC(dataq, pc->x(ip), new_hc);
	for (int idata = 0; idata < pc->data_m; ++idata) {
	  if (pc->dataStructVec[idata].error_flag) { 
	    pc->error(ip,idata) = fabs((pc->data(ip,idata) - dataq[idata])/(pc->data_max[idata] - pc->data_min[idata]));
	    max_error = max(max_error,pc->error(ip,idata));
	  }
	}
      }
      
      // if we haven't met tolerance, add to list...
      if (max_error >= rtol) {
	if (new_hc->depth < max_depth) {
	  active_leaves.push_back(right);
	}
	else {
	  for (list<int>::iterator ilist = new_hc->ipList.begin(); ilist != new_hc->ipList.end(); ++ilist) 
	    ++count_zeroed;
	  max_error_zeroed = max(max_error_zeroed, max_error);
	} 
      }     

      ++iter;
      
      if (vecPtrHC.size()%1000 == 0)
	cout << "       > # of leaves, # of active leaves: " <<  vecPtrHC.size() << ", " << active_leaves.size() << endl;
      
      if (vecPtrHC.size() >= max_nhc) {
	max_error = 0.0;
	for (typename list<Node*>::iterator ileaf = active_leaves.begin(); ileaf != active_leaves.end(); ++ileaf) {
	  HyperCube * hc = (*ileaf)->hc;
	  for (list<int>::iterator ilist = hc->ipList.begin(); ilist != hc->ipList.end(); ++ilist) {
	    const int ip = *ilist;
	    for (int idata = 0; idata < pc->data_m; ++idata)
	      if (pc->dataStructVec[idata].error_flag) 
		max_error = max(max_error, pc->error(ip,idata));
	  }
	}
	cout << "       > WARNING: maximum number of leaves (" << max_nhc << ") reached, exiting iteration with relative error " << max_error << " (tolerance is " << rtol << ")" << endl; 
	break;
      }
    }

    if (count_zeroed > 0) {
      COUT2_S("      > WARNING: maximum depth of " << max_depth << " reached, ignoring " << count_zeroed << 
	      " (" << (double)100*count_zeroed/pc->np << " %) points with maximum relative error " << max_error_zeroed);
    }

    delete[] dataq;
  }
  
  void queryHC(double * data, const double xq[DIM_K]) {
    
    // first, find containing hypercube for xq 
    Node * this_leaf = findLeaf(xq);
    
    // then perform k-linear interpolation with base HC
    linearInterpHC(data, xq, this_leaf->hc);
  }
  
  void queryHC_split(double * data, const double xq[DIM_K], const Node * this_leaf, const int idim) {
    
    // construct child HC according to the split dimension argument passed...
    HyperCube tmphc = *(this_leaf->hc);
    // for each split, the query point can lie to the left or the right of the split plane...
    if (xq[idim] < vecPtrVe[tmphc.ve_split[idim][0]]->x[idim]) {
      int ve_split_count = 0;
      for (int ive = 0; ive < (1<<DIM_K); ++ive) {
	if (ive&(1<<idim)) {
	  tmphc.ve[ive] = tmphc.ve_split[idim][ve_split_count];
	  ++ve_split_count;
	}
      }
      assert(ve_split_count == (1<<(DIM_K-1)));
    }
    else {
      int ve_split_count = 0;
      for (int ive = 0; ive < (1<<DIM_K); ++ive) {
	if (!(ive&(1<<idim))) {
	  tmphc.ve[ive] = tmphc.ve_split[idim][ve_split_count];
	  ++ve_split_count;
	}
      }
      assert(ve_split_count == (1<<(DIM_K-1)));
    }
    
    // then perform k-linear interpolation with split HC containing query point
    linearInterpHC(data, xq, &tmphc);
  }
  
  void linearInterpHC(double * data, const double xq[DIM_K], const HyperCube * this_hc) {
    // k-linear interpolation of data at point xq in HyperCube hc
    double weights[1<<DIM_K];
    for (int ive = 0; ive < (1<<DIM_K); ++ive) weights[ive] = 1.0;
    
    for (int idim = 0; idim < DIM_K; ++idim) {
      const double inv_denom = 1.0/(vecPtrVe[this_hc->ve[(1<<DIM_K)-1]]->x[idim] - vecPtrVe[this_hc->ve[0]]->x[idim]);
      for (int ive = 0; ive < (1<<DIM_K); ++ive) { 
	const double tmp = (xq[idim] - vecPtrVe[this_hc->ve[0]]->x[idim])*inv_denom; 
	if (ive&(1<<idim))
	  weights[ive] *= tmp;	    
	else
	  weights[ive] *= (1.0 - tmp);
	assert(weights[ive] >= 0.0); // always interpolate
      }      
    }
    
    for (int idata = 0; idata < pc->data_m; ++idata) {
      data[idata] = 0.0;
      for (int ive = 0; ive < (1<<DIM_K); ++ive) data[idata] += weights[ive]*vecPtrVe[this_hc->ve[ive]]->data[idata];
    }
  }
  
  Node * findLeaf(const double x[DIM_K]) {
    for (int idim = 0; idim < DIM_K; ++idim) assert((x[idim] >= pc->x_min[idim]) && (x[idim] <= pc->x_max[idim]));
    
    Node * this_node = root_node;
    while (this_node->hc == NULL) {
      if (x[this_node->dim_split] < this_node->x_split)
	this_node = this_node->left;
      else
	this_node = this_node->right;
    }
    
#ifdef CTI_DEBUG 
    // make sure we've actually found it...
    assert(pointIsInsideLeaf(this_node, x));   
#endif
    
    return(this_node);
  }

  bool pointIsInsideLeaf(const Node * this_leaf, const double x[DIM_K]) {
    // checks if hypercube corresponding to given leaf contains given point    
    HyperCube * this_hc = this_leaf->hc;    
    for (int ive = 0; ive < (1<<DIM_K); ++ive) {
      for (int idim = 0; idim < DIM_K; ++idim) {
	if (ive&(1<<idim)) {
	  if (x[idim] > vecPtrVe[this_hc->ve[ive]]->x[idim]) {
	    return(false);
	  }
	}
	else {
	  if (x[idim] < vecPtrVe[this_hc->ve[ive]]->x[idim]) {
	    return(false);
	  }
	}
      }
    }
    
    return(true);
  }
 
  void splitLeaf(Node*& left, Node*& right, Node * this_leaf, const int dim_split) {
    // splits given leaf/hypercube along given dimension...
    
    assert((dim_split >= 0) && (dim_split < DIM_K));
    
    // make sure not branch..
    assert(this_leaf->hc != NULL);
    
    HyperCube * this_hc = this_leaf->hc;
    this_leaf->hc        = NULL;
    this_leaf->dim_split = dim_split;
    this_leaf->x_split   = vecPtrVe[this_hc->ve_split[dim_split][0]]->x[dim_split];
          
    HyperCube * new_hc = new HyperCube;  // new is always "right" hc (1 bit)
    int ve_split_count = 0;
     for (int ive = 0; ive < (1<<DIM_K); ++ive) {
       if (!(ive&(1<<dim_split))) {
	 new_hc->ve[ive] = this_hc->ve_split[dim_split][ve_split_count];
	 ++ve_split_count;
       }
       else 
	 new_hc->ve[ive] = this_hc->ve[ive];
     }
     assert(ve_split_count == (1<<(DIM_K-1)));
     ve_split_count = 0;
     for (int ive = 0; ive < (1<<DIM_K); ++ive) {
       if (ive&(1<<dim_split)) {
	 this_hc->ve[ive] = this_hc->ve_split[dim_split][ve_split_count];
	 ++ve_split_count;
       }
     }
     assert(ve_split_count == (1<<(DIM_K-1)));
    
     // init candidate splits (i.e. new verticies on unsplit edges) for both hypercubes after split
     initCandidates(new_hc);
     initCandidates(this_hc);
     
     new_hc->index = vecPtrHC.size();
     vecPtrHC.push_back(new_hc);
     
     // append binary tree info
     left = new Node;
     left->hc = this_hc;
     this_hc->depth += 1;
     left->index = vecPtrNode.size();
     vecPtrNode.push_back(left);

     right = new Node;
     right->hc = new_hc;
     new_hc->depth = this_hc->depth;
     right->index = vecPtrNode.size();
     vecPtrNode.push_back(right);

     this_leaf->left      = left;
     this_leaf->right     = right;
  }

  void initCandidates(HyperCube * thishc) {
    
    // initialize candidate splits
    for (int idim = 0; idim < DIM_K; ++idim) {
      int ve_count = 0;
      // loop over lower-dim hypercubes in each dimension by toggling bits and add verticies...
      for (int ive = 0; ive < (1<<DIM_K); ++ive) {
	if (ive&(1<<idim)) {
	  Vertex * new_vertex = new Vertex(pc->data_m);
	  for (int j = 0; j < DIM_K; ++j) new_vertex->coords[j] = vecPtrVe[thishc->ve[ive]]->coords[j];
	  
	  assert(vecPtrVe[thishc->ve[ive]]->coords[idim] > vecPtrVe[thishc->ve[ive^(1<<idim)]]->coords[idim]);
	  assert((vecPtrVe[thishc->ve[ive]]->coords[idim] - vecPtrVe[thishc->ve[ive^(1<<idim)]]->coords[idim])%2 == 0);
	  
	  new_vertex->coords[idim] = vecPtrVe[thishc->ve[ive^(1<<idim)]]->coords[idim] + (vecPtrVe[thishc->ve[ive]]->coords[idim] - vecPtrVe[thishc->ve[ive^(1<<idim)]]->coords[idim])/2;
	  addVertex(new_vertex);
	  thishc->ve_split[idim][ve_count] = new_vertex->ive;
	  ++ve_count;
	}
      }
      assert(ve_count == (1<<(DIM_K-1)));  
    } 
  }
    
  void addVertex(Vertex * &new_vertex) {
    for (int idim = 0; idim < DIM_K; ++idim) 
      assert((new_vertex->coords[idim] >= 0) && (new_vertex->coords[idim] <= RANGE));
    
    Vertex * this_vertex = root_vertex;
    bool found = false;
    int leaf   = 0;
    // check to make sure vertex does not already exist...
    while (!found && (leaf == 0)) {
      found = true;
      for (int idim = 0; idim < DIM_K; ++idim) {
	// ordering of verticies is based on dimenion-wise coordinate ordering
	if (this_vertex->coords[idim] != new_vertex->coords[idim]) {
	  if (this_vertex->coords[idim] < new_vertex->coords[idim]) {
	    if ( this_vertex->left != NULL ) 
	      this_vertex = this_vertex->left;
	    else 
	      leaf = -1; // -1 is left
	  }
	  else {
	    if ( this_vertex->right != NULL ) 
	      this_vertex = this_vertex->right;
	    else 
	      leaf = 1; // +1 is right
	  }
	  found = false;
	  break;
	}
      }
    }

    if (found) {
      // if found, return the found vertex 
      for (int idim = 0; idim < DIM_K; ++idim) assert(this_vertex->coords[idim] == new_vertex->coords[idim]);
      delete new_vertex;
      new_vertex = this_vertex;
    }
    else {
      // otherwise, the vertex does not exist yet, so compute x location and pushback...
      if (leaf == -1) {
	this_vertex->left = new_vertex;
      }
      else {
	assert(leaf == 1);
	this_vertex->right = new_vertex; 
      }

      for (int idim = 0; idim < DIM_K; ++idim) 
	new_vertex->x[idim] = pc->x_min[idim] + dx[idim]*new_vertex->coords[idim];   
      pc->queryPC(new_vertex->data, new_vertex->x);      
      
      new_vertex->ive = vecPtrVe.size();
      vecPtrVe.push_back(new_vertex);
    }

  }  
  
  void finalize() {
    cout << "finalize()" << endl;
    writeBinary(filename);     
    
    if (checkParam("CHEMTABLE.TECPLOT"))     
      writeTecplot();
    
    if (cti_verbose) 
      info();
  }

  void info() {
    
    cout << "\n====================== chemtable info ======================"<<endl;
    cout << "   type:                     " << tabletype << endl;
       
    cout << "\n";
    cout << "   number of verticies:      " << vecPtrVe.size() << endl;
    cout << "   number of tree nodes:     " << vecPtrNode.size() << endl;
    cout << "   number of leaves:         " << vecPtrHC.size() << endl;
    cout << "\n";
    
    // point cloud is in charge of dumping its info - kdtree shouldn't know names, pressure, etc.
    pc->info();
    cout << "============================================================\n"<<endl;

  }
  
  void writeBinary(const string& filename) {
    cout << " > writing kd table to the binary file " << filename << "... " << endl;

    // open the file
    FILE * fp = NULL;
    fp = fopen(filename.c_str(),"wb");
    if (fp == NULL) 
      CERR_S("cannot open chemtable file: " << filename);
        
    // (1) write ugp magic number and version
    {
      int tmpint[2]={CHEM_IO_MAGIC_NUMBER,CHEM_IO_VERSION};
      size_t dummy= fwrite(tmpint, sizeof(int), 2, fp);
      assert(dummy == 2);
    }
    
    // x and data name strings come from point cloud - also, are indexed by vertex ordering in stl vector        
    int nve   = vecPtrVe.size();   assert(nve > 0);
    int nhc   = vecPtrHC.size();   assert(nhc > 0); 
    int nnode = vecPtrNode.size(); assert(nnode > 0);   

    // (2) write the table type and dimensions
    {
      Header header;
      // copy the name
      for(int i=0; i<UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';
      size_t dummy = tabletype.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy<UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      // id and byte skip
      header.id = UGP_IO_CT_KD;
      header.skip = sizeof(Header);
      // table dimensions
      header.idata[0] = DIM_K; 
      header.idata[1] = pc->data_m;
      header.idata[2] = nve;
      header.idata[3] = nhc;
      header.idata[4] = nnode;
            
      header.rdata[0] = rtol;
      // allow point cloud to fill some info...
      pc->extraMainHeader(header, 5, 1);

      // write into the file      
      dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);
    }

    // use tmp to move vertex data to be contiguous in memory 
    double * tmp = new double[nve];
    
    assert(pc->xNameVec.size()    == DIM_K);
    assert(pc->dataStructVec.size() == pc->data_m);

    // (3) write vertex x coordinates
    cout << "  vertex coordinates... " << endl;
    for (int idim = 0; idim < DIM_K; ++idim) {
      Header header;
      // copy the name
      for(int i = 0; i < UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      string my_name = pc->xNameVec[idim];
      cout << "      " << my_name <<  endl;
      size_t dummy = my_name.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy < UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      
      // id and byte skip
      header.id = UGP_IO_CT_COOR;
      header.skip = nve*sizeof(double) + sizeof(Header);
      
      // number of verticies
      header.idata[0] = nve; 
      header.idata[1] = idim;
      
      // write header
      dummy = fwrite(&header, sizeof(Header), 1, fp);
      assert(dummy == 1);
      
      for (int ive = 0; ive < nve; ++ive)
	tmp[ive] = vecPtrVe[ive]->x[idim];
      
      // write coordinate
      dummy = fwrite(tmp, sizeof(double), header.idata[0], fp);
      assert(dummy == header.idata[0]);
    }
    
    // (4) write vertex data
    cout << "  vertex data... " << endl;
    for (int idata = 0; idata < pc->data_m; ++idata) {
      Header header;
      // copy the name
      for(int i = 0; i < UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      string my_name = pc->dataStructVec[idata].name;
      cout << "      " << my_name << endl;
      size_t dummy = my_name.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy < UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');

      // id and byte skip
      header.id = UGP_IO_CT_DATA;
      header.skip = nve*sizeof(double) + sizeof(Header);

      // number of verticies
      header.idata[0] = nve; 
      header.idata[1] = idata;
      
      // write header
      dummy = fwrite(&header, sizeof(Header), 1, fp);
      assert(dummy == 1);
      
      for (int ive = 0; ive < nve; ++ive)
	tmp[ive] = vecPtrVe[ive]->data[idata];
      
      // write data
      dummy = fwrite(tmp, sizeof(double), header.idata[0], fp);
      assert(dummy == header.idata[0]);
    }
    
    delete[] tmp;
    
    // temporary structure for writing out hypercube connectivity
    int (*tmp_veohc)[1<<DIM_K] = new int[nhc][1<<DIM_K];
    
    for (int ihc = 0; ihc < nhc; ++ihc) 
      for (int ive = 0; ive < (1<<DIM_K); ++ive) 
	tmp_veohc[ihc][ive] = vecPtrHC[ihc]->ve[ive];
    
    // (5) write hypercube connectivity
    cout << "  hypercube connectivity... " << endl;
    { 
      Header header;
      // copy the name
      for(int i = 0; i < UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      string my_name = "VEOHC";
      size_t dummy = my_name.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy < UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');

      // id and byte skip
      header.id = UGP_IO_CT_KD_HC;
      header.skip = (1<<DIM_K)*nhc*sizeof(int) + sizeof(Header);

      // length of structure
      header.idata[0] = nhc*(1<<DIM_K); 

      // write header
      dummy = fwrite(&header, sizeof(Header), 1, fp);
      assert(dummy == 1);
      
      // write connectivity
      dummy = fwrite(tmp_veohc, sizeof(int), header.idata[0], fp);
      assert(dummy == header.idata[0]);
    }
    delete[] tmp_veohc;
        
    // due to non-portability of writing structs (because of padding), 
    // nodes are written after serialization...
    int * tmp_dim_split = new int[nnode];
    double * tmp_x_split = new double[nnode];
    int * tmp_left  = new int[nnode];
    int * tmp_right = new int[nnode];
    int * tmp_hc    = new int[nnode];
    
    for (int inode = 0; inode < nnode; ++inode) {
      tmp_dim_split[inode] = vecPtrNode[inode]->dim_split;
      tmp_x_split[inode]   = vecPtrNode[inode]->x_split;
      if (vecPtrNode[inode]->left != NULL) 
	tmp_left[inode] = (vecPtrNode[inode]->left)->index;
      else
	tmp_left[inode] = -1;
      if (vecPtrNode[inode]->right != NULL)
	tmp_right[inode] = (vecPtrNode[inode]->right)->index;      
      else
	tmp_right[inode] = -1;
      if (vecPtrNode[inode]->hc != NULL)
	tmp_hc[inode] = (vecPtrNode[inode]->hc)->index;      
      else 
	tmp_hc[inode] = -1;
    }

    // (6) write tree 
    cout << "  tree structure... " << endl;
    { 
      Header header;
      // copy the name
      for(int i = 0; i < UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      string my_name = "DIM_SPLIT";
      size_t dummy = my_name.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy < UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      
      // id and byte skip
      header.id = UGP_IO_CT_KD_NODE;
      header.skip = nnode*sizeof(int) + sizeof(Header);
      
      // length of int data
      header.idata[0] = nnode;
      
      // write header
      dummy = fwrite(&header, sizeof(Header), 1, fp);
      assert(dummy == 1);
      
      // write node struct ints serially
      dummy  = fwrite(tmp_dim_split, sizeof(int), nnode, fp);
      assert(dummy == header.idata[0]);
    }
    delete[] tmp_dim_split;
    
    { 
      Header header;
      // copy the name
      for(int i = 0; i < UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      string my_name = "X_SPLIT";
      size_t dummy = my_name.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy < UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      
      // id and byte skip
      header.id = UGP_IO_CT_KD_NODE;
      header.skip = nnode*sizeof(double) + sizeof(Header);
      
      // length of int data
      header.idata[0] = nnode;
      
      // write header
      dummy = fwrite(&header, sizeof(Header), 1, fp);
      assert(dummy == 1);
      
      // write node struct ints serially
      dummy  = fwrite(tmp_x_split, sizeof(double), nnode, fp);
      assert(dummy == header.idata[0]);
    }
    delete[] tmp_x_split;
    
    { 
      Header header;
      // copy the name
      for(int i = 0; i < UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      string my_name = "LEFT";
      size_t dummy = my_name.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy < UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      
      // id and byte skip
      header.id = UGP_IO_CT_KD_NODE;
      header.skip = nnode*sizeof(int) + sizeof(Header);
      
      // length of int data
      header.idata[0] = nnode;
      
      // write header
      dummy = fwrite(&header, sizeof(Header), 1, fp);
      assert(dummy == 1);
      
      // write node struct ints serially
      dummy = fwrite(tmp_left, sizeof(int), nnode, fp);
      assert(dummy == header.idata[0]);
    }
    delete[] tmp_left;
    
    { 
      Header header;
      // copy the name
      for(int i = 0; i < UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      string my_name = "RIGHT";
      size_t dummy = my_name.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy < UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      
      // id and byte skip
      header.id = UGP_IO_CT_KD_NODE;
      header.skip = nnode*sizeof(int) + sizeof(Header);
      
      // length of int data
      header.idata[0] = nnode;
      
      // write header
      dummy = fwrite(&header, sizeof(Header), 1, fp);
      assert(dummy == 1);
      
      // write node struct ints serially
      dummy = fwrite(tmp_right, sizeof(int), nnode, fp);
      assert(dummy == header.idata[0]);
    }
    delete[] tmp_right;

    { 
      Header header;
      // copy the name
      for(int i = 0; i < UGP_IO_HEADER_NAME_LEN; ++i)
	header.name[i] = '\0';	
      string my_name = "HC";
      size_t dummy = my_name.copy(header.name,UGP_IO_HEADER_NAME_LEN);
      assert(dummy < UGP_IO_HEADER_NAME_LEN);
      assert(header.name[dummy] == '\0');
      
      // id and byte skip
      header.id = UGP_IO_CT_KD_NODE;
      header.skip = nnode*sizeof(int) + sizeof(Header);
      
      // length of int data
      header.idata[0] = nnode;
      
      // write header
      dummy = fwrite(&header, sizeof(Header), 1, fp);
      assert(dummy == 1);
      
      // write node struct ints serially
      dummy = fwrite(tmp_hc, sizeof(int), nnode, fp);
      assert(dummy == header.idata[0]);
    }
    delete[] tmp_hc;
   
    // end of file header
    {
      Header header;
      header.id = UGP_IO_EOF;
      sprintf(header.name,"EOF");
      header.skip = sizeof(Header);

      // write the header into the file      
      size_t dummy = fwrite(&header,sizeof(Header),1,fp);
      assert(dummy == 1);  
    }

    // close the file
    fclose(fp);
  
  }

  void writeTecplot() {
    string filename = tablename + ".dat";
    cout << " > writing table to the tecplot file: " << filename << endl;
    
    FILE* fp = fopen(filename.c_str(), "w");
    fprintf(fp,"TITLE = \"%s\"\n",filename.c_str());
    fprintf(fp,"VARIABLES = ");
    
    for(int idim = 0; idim < DIM_K; ++idim) 
      fprintf(fp,"\"%s\"\n", pc->xNameVec[idim].c_str()); 
    
    for(int idata = 0; idata < pc->data_m; ++idata)
      fprintf(fp,"\"%s\"\n", pc->dataStructVec[idata].name.c_str());
    
    fprintf(fp,"\"Zvar\"\n");
    
    fprintf(fp,"ZONE T=\"%s\"\n",filename.c_str());
    if (DIM_K == 2)
      fprintf(fp,"N=%u, E=%u, F=FEPOINT, ET=QUADRILATERAL\n",int(vecPtrVe.size()),int(vecPtrHC.size()));
    else if (DIM_K == 3)
      fprintf(fp,"N=%u, E=%u, F=FEPOINT, ET=BRICK\n",int(vecPtrVe.size()),int(vecPtrHC.size()));
    else 
      CERR_S("TECPLOT output for dimension >= 4 is not supported yet");
    
    for (int ive = 0; ive < vecPtrVe.size(); ++ive) {
      for (int idim = 0; idim < DIM_K; ++idim) 
	fprintf(fp, "%g ", vecPtrVe[ive]->x[idim]);
      
      for (int idata = 0; idata < pc->data_m; ++idata) 
	fprintf(fp, "%g ", vecPtrVe[ive]->data[idata]);
      
      double this_zvar = vecPtrVe[ive]->x[1]*(vecPtrVe[ive]->x[0]*(1.0 - vecPtrVe[ive]->x[0]));
      fprintf(fp, "%g ", this_zvar);
      
      fprintf(fp,"\n");
    }
    
    for (int ihc = 0; ihc < vecPtrHC.size(); ++ihc) {
      // 1-indexing for tecplot...
      if (DIM_K == 2)
	fprintf(fp, "%d %d %d %d\n", vecPtrHC[ihc]->ve[0]+1, vecPtrHC[ihc]->ve[1]+1, vecPtrHC[ihc]->ve[3]+1, vecPtrHC[ihc]->ve[2]+1);
      else if (DIM_K == 3)
	fprintf(fp, "%d %d %d %d %d %d %d %d\n", vecPtrHC[ihc]->ve[0]+1, vecPtrHC[ihc]->ve[1]+1, vecPtrHC[ihc]->ve[3]+1, vecPtrHC[ihc]->ve[2]+1, vecPtrHC[ihc]->ve[4]+1, vecPtrHC[ihc]->ve[5]+1, vecPtrHC[ihc]->ve[7]+1, vecPtrHC[ihc]->ve[6]+1);
      else 
	CERR_S("TECPLOT output for dimension >= 4 is not supported yet");
    }

   
    fclose(fp);
  }
    
};

class VidaChemtableNFPVKdt3D : public KdTable<3> {

public:

  VidaChemtableNFPVKdt3D(const string& tabletype) {
    cout << "VidaChemtableNFPVKdt3D()" << endl;
    this->tabletype = tabletype;
    filename = tablename + ".VidaNFPVKdt3D.chem"; 
  }

  ~VidaChemtableNFPVKdt3D() {
    cout << "~VidaChemtableNFPVKdt3D()" << endl;
  }

  void init() {
    // initialize appropriate point cloud...
    pc = new VidaNonPremixed3DPC();
    
    KdTable<3>::init();
  }
  
};

class ChrisChemtableNFPVKdt3D : public KdTable<3> {

public:

  ChrisChemtableNFPVKdt3D(const string& tabletype) {
    cout << "ChrisChemtableNFPVKdt3D()" << endl;
    this->tabletype = tabletype;
    filename = tablename + ".ChrisNFPVKdt3D.chem"; 
  }

  ~ChrisChemtableNFPVKdt3D() {
    cout << "~ChrisChemtableNFPVKdt3D()" << endl;
  }

  void init() {
    // initialize appropriate point cloud...
    pc = new ChrisNonPremixed3DPC();
    
    KdTable<3>::init();
  }

 
};

// ================================================================================ //
//                             END KDT-SPECIFIC CODE                                //
// ================================================================================ //


#endif

