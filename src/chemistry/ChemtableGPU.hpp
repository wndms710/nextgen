#ifndef CHEMTABLEGPU_HPP
#define CHEMTABLEGPU_HPP

#include "AbstractChemtable.hpp"
#include <cuda_runtime.h>

using namespace MiscUtils;

__global__
void lookupSpecialKernel(double * rout, const double * y1, const double * y2a,
    const double * y2b, const double * data, const double* x1d, const double* x2d,
    const double* x2lowerd, const double* x2upperd, const double (*invDenom1d)[4], const double (*invDenom2d)[4],
    indexMap* idxMap1d, indexMap* idxMap2d, const int n1,
    const int n2, const int n);

__global__
void lookupCubicPtrVectorKernel(double** routd, const double* y1, const double *y2,
    double** datav, const double* x1d, const double* x2d, const double* x2lowerd,
    const double* x2upperd, const double (*invDenom1d)[4], const double (*invDenom2d)[4],
    indexMap* idxMap1d, indexMap* idxMap2d, const int ndata, const int n1, const int n2, const int n);

__global__
void lookupReducedVectorKernel(double** routd, const double* y1,
    double** datav, const double* x1d, const double (*invDenom1d)[4],
    indexMap* idxMap1d, const int ndata, const int n1, const int n2, const int n);

// this kernel could be moved elsewhere, but for the purposes of keeping
// the cuda code together, its defined here currently..
__global__
void computeFaceReductionKernel(double * Z_fa_d,double * C0_fa_d,double * C1_fa_d,
                                const double * Z_cv_d, const double * C_cv_d, const int (*cvofa_d)[2],const int n);


void computeFaceReduction(double* Z_fa_d, double* C0_fa_d, double * C1_fa_d,
                          const double* Z_cv_d, const double* C_cv_d, const int (*cvofa_d)[2], const int n);

class CartesianChemtable2dGpu: public AbstractChemtable2D {

private:

  string tableName;
  int n1, n2, nvars;
  double *x1, *x2;
  double *x2lower, *x2upper;
  double (*invDenom1)[4], (*invDenom2)[4];
  indexMap *idxMap1, *idxMap2;

  //=====
  // cuda copies
  //=====
  double *x1d, *x2d;
  double *x2lowerd, *x2upperd, *x2upper2dd;
  double (*invDenom1d)[4], (*invDenom2d)[4];
  indexMap *idxMap1d, *idxMap2d;

  string *varName;
  list<ChemtableData2D> datalist;
  typedef std::pair<int,double> myPair;
  typedef std::pair<int,int> myPairInt;
  double SMALL_CLIP;

  map<string,double*> deviceVars;
  map<string,double*> cudaBuffers;

public:

  // constructor
  CartesianChemtable2dGpu(const string& _name): n1(0), n2(0), nvars(0),
  x1(NULL), x2(NULL), x2lower(NULL), x2upper(NULL), invDenom1(NULL),
  invDenom2(NULL), idxMap1(NULL), idxMap2(NULL), varName(NULL), tableName(_name) {

    if ( mpi_rank == 0 )
      cout << "CartesianChemtable2dGpu():" << endl;

    x1d = NULL;
    x2d = NULL;
    x2lowerd = NULL;
    x2upperd = NULL;
    x2upper2dd = NULL;
    invDenom1d = NULL;
    invDenom2d = NULL;
    idxMap1d = NULL;
    idxMap2d = NULL;

    deviceVars.clear();
    SMALL_CLIP = 1.01e-6;
    init();
  }

  ~CartesianChemtable2dGpu(){
    if ( mpi_rank == 0 )
      cout << "~CartesianChemtable2D():" << endl;
    DELETE(x1);
    DELETE(x2);
    DELETE(x2lower);
    DELETE(x2upper);
    DELETE(varName);
    DELETE(invDenom1);
    DELETE(invDenom2);
    delete idxMap1;
    delete idxMap2;
  }

  void info(){

    if ( mpi_rank == 0 ) {
      cout << "============== chemistry table information ===============" << endl;
      cout << " chemistry table information: " << endl;
      cout << "  > name:                     " << tableName << endl;
      //cout << "  > type:                     " << tableType << endl;
      cout << "  > reference pressure:       " << pressure << endl;
      cout << "  > dimension:                ("<< n1 <<"," << n2 << ")" << endl;
      cout << "  > first coordinate range:   ("<< x1[0] <<"," << x1[n1-1] << ")" << endl;
      cout << "  > second coordinate range:  ("<< x2[0] <<"," << x2[n2-1] << ")" << endl;
      cout << "  > variables("<< nvars << "):"<< endl;
      cout << "     ";
      for (int var=0; var<nvars; ++var) {
	cout <<  varName[var];
	if ( var == (nvars-1) )
	  cout << "."<< endl;
	else
	  cout << ", ";
      }
    }

    // for all the variables
    if ( datalist.begin() == datalist.end() ){
      if ( mpi_rank == 0 )
	cout << "  > no variables are currently loaded " <<endl;
    }
    else {
      if ( mpi_rank == 0 )
	cout << "  > the following variables are currently loaded: " << endl;

      for (list<ChemtableData2D>::iterator data = datalist.begin(); data != datalist.end(); ++data) {
	double varmin;
	double varmax;
	assert(getTableVarRange(varmin,varmax,data->name));
	if ( mpi_rank == 0 )
	  cout << "   > " << data->name << ": "<< varmin << "<"<<data->name <<"<" << varmax<< endl ;
      }
    }

    if ( mpi_rank == 0 )
      cout << "=========== end of chemistry table information ============" << endl;
  }

  int getTableVarRange(double& varmin, double& varmax, const string& name){
    // initialze
    varmin = 0.0;
    varmax = 0.0;
    // loop through the list
    for (list<ChemtableData2D>::iterator data = datalist.begin(); data != datalist.end(); ++data) {
      if ( data->name == name ) {
	varmin = (*data)(0,0);
	varmax = (*data)(0,0);
	for ( int i=0; i<n1; ++i)
	  for ( int j=0; j<n2; ++j){
	    varmax = max(varmax,(*data)(i,j));
	    varmin = min(varmin,(*data)(i,j));
	  }
	return(1);
      }
    }
    if ( mpi_rank == 0)
      cout << "***WARNING in getTableVarRange: the variable "
	   << name << " is not loaded yet/n" << " ...setting the range to 0:0 " << endl;

    return(0);
  }

  void getRhoCsrcRange(double& Smin, double& Smax) {
    ChemtableData2D * rho = getVar("rho");
    ChemtableData2D * src = getVar("src_prog");
    Smin =  1.0e+20;
    Smax = -1.0e+20;
    for (int i=0; i<n1; ++i) {
      for (int j=0; j<n2; ++j) {
        const double rhoSrc = (*rho)(i,j) * (*src)(i,j);
        Smin = min(Smin,rhoSrc);
        Smax = max(Smax,rhoSrc);
      }
    }
    return;
  }

  // ==========================================
  // write the table with the current variables
  // the one that are loaded to a tecplot file
  // ==========================================
  void writeTecplot(const string& filename) {

    // load all the variables
    for (int i=0; i<nvars; ++i)
      ChemtableData2D* tmp = getVar(varName[i]);

    // only rank 0 needs to write
    if (mpi_rank == 0) {

      cout << " > writing 2d cartesian chemtable " << tableName << " to the Tecplot file "<< filename << endl;
      cout << "  > the following variables are currently loaded: " << endl;
      for (list<ChemtableData2D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
	cout << "   > " << data->name << endl ;

      FILE * fp = fopen(filename.c_str(),"w");
      fprintf(fp,"title = \"chemtable\"\n");
      fprintf(fp,"variables = \"Z\"\n\"C\"\n");
      for (list<ChemtableData2D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
	fprintf(fp,"\"%s\"\n",data->name.c_str());
      fprintf(fp,"zone i = %d, j = %d, DATAPACKING=POINT\n",n1,n2);
      for (int j=0; j <n2; ++j)
	for (int i=0; i <n1; ++i){
	  fprintf(fp,"%lf %lf",x1[i],x2[j]);
	  for (list<ChemtableData2D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
	    fprintf(fp," %lf",(*data)(i,j));
	  fprintf(fp,"\n");
	}
      fclose(fp);
    }
  }

public:

  void lookupDataVec(vector<double*>& routPtrVec, const vector<string>& nameVec,
		     const double* x1, const double* x2, const int n) {
    //lookupLinearPtrVector(routPtrVec,nameVec,x1,x2,n);
    lookupCubicPtrVector(routPtrVec,nameVec,x1,x2,n);
  }

  void lookupDataVecReduced(vector<double*>& routPtrVec, const vector<string>& nameVec,
                            const double * y1, const int n);

  void init() {

    readHeaderAndCoordinates();
    info();
    initProgvarBounds();
    initInterpWts();

    // fast indexing finding...
    idxMap1 = new indexMap(x1,n1);
    idxMap2 = new indexMap(x2,n2);

  }

public:
  void initGpu() {

    COUT1("chemtable init gpu");

    // the cpu side has been inititalized at this point,
    // we need to now copy some of the data over the device.
    // the table coordinates, the index finding map, the
    // table variables, and precomputed interpolation weights..
    int n_bytes;

    x1d        = NULL;
    n_bytes    = n1*sizeof(double);
    cudaMalloc((void**)&x1d,n_bytes);
    cudaMemcpy(x1d,x1,n_bytes,cudaMemcpyHostToDevice);

    x2d        = NULL;
    n_bytes    = n2*sizeof(double);
    cudaMalloc((void**)&x2d,n_bytes);
    cudaMemcpy(x2d,x2,n_bytes,cudaMemcpyHostToDevice);

    x2lowerd   = NULL;
    n_bytes    = n1*sizeof(double);
    cudaMalloc((void**)&x2lowerd,n_bytes);
    cudaMemcpy(x2lowerd,x2lower,n_bytes,cudaMemcpyHostToDevice);

    x2upperd   = NULL;
    n_bytes    = n1*sizeof(double);
    cudaMalloc((void**)&x2upperd,n_bytes);
    cudaMemcpy(x2upperd,x2upper,n_bytes,cudaMemcpyHostToDevice);

    invDenom1d = NULL;
    n_bytes    = (n1-3)*4*sizeof(double);
    cudaMalloc((void**)&invDenom1d,n_bytes);
    cudaMemcpy(invDenom1d,invDenom1,n_bytes,cudaMemcpyHostToDevice);

    invDenom2d = NULL;
    n_bytes    = (n2-3)*4*sizeof(double);
    cudaMalloc((void**)&invDenom2d,n_bytes);
    cudaMemcpy(invDenom2d,invDenom2,n_bytes,cudaMemcpyHostToDevice);

    //=====================================
    // the index map is a bit of a pain, although we could simplify
    // the following using the unified memory address...
    //=====================================
    cudaMalloc((void**)&idxMap1d,sizeof(indexMap));
    cudaMalloc((void**)&idxMap2d,sizeof(indexMap));

    cudaMemcpy(idxMap1d,idxMap1,sizeof(indexMap), cudaMemcpyHostToDevice);
    cudaMemcpy(idxMap2d,idxMap2,sizeof(indexMap), cudaMemcpyHostToDevice);

    // deep copy of the xloc, iloc arrays ...
    double * xloc1d = NULL, *iloc1d = NULL;
    n_bytes = (idxMap1->nloc)*sizeof(double);
    cudaMalloc((void**)&xloc1d, n_bytes);
    cudaMalloc((void**)&iloc1d, n_bytes);
    cudaMemcpy(xloc1d,idxMap1->xloc,n_bytes,cudaMemcpyHostToDevice);
    cudaMemcpy(iloc1d,idxMap1->iloc,n_bytes,cudaMemcpyHostToDevice);

    // and link the ptrs..
    cudaMemcpy(&(idxMap1d->xloc),&xloc1d,sizeof(double*),cudaMemcpyHostToDevice);
    cudaMemcpy(&(idxMap1d->iloc),&iloc1d,sizeof(double*),cudaMemcpyHostToDevice);

    double * xloc2d = NULL, *iloc2d = NULL;
    n_bytes = (idxMap2->nloc)*sizeof(double);
    cudaMalloc((void**)&xloc2d, n_bytes);
    cudaMalloc((void**)&iloc2d, n_bytes);
    cudaMemcpy(xloc2d,idxMap2->xloc,n_bytes,cudaMemcpyHostToDevice);
    cudaMemcpy(iloc2d,idxMap2->iloc,n_bytes,cudaMemcpyHostToDevice);
    cudaMemcpy(&(idxMap2d->xloc),&xloc2d,sizeof(double*),cudaMemcpyHostToDevice);
    cudaMemcpy(&(idxMap2d->iloc),&iloc2d,sizeof(double*),cudaMemcpyHostToDevice);

    // and the chemtable variables ..
    for(list<ChemtableData2D>::iterator data = datalist.begin(); data != datalist.end(); ++data) {
      double * tmp = NULL;
      n_bytes = n1*n2*sizeof(double);
      cudaMalloc((void**)&tmp,n_bytes);
      cudaMemcpy(tmp,data->data,n_bytes,cudaMemcpyHostToDevice);
      deviceVars[data->name] = tmp;
    }

    // making a 2d version of Cmax available on the gpu to avoid
    // having to branch inside of the cuda kernels...
    {
      double* tmp_cpu = new double[n1*n2];
      for (int i = 0; i < n1 ; ++i) {
        for (int j =0; j < n2; ++j) {
          tmp_cpu[(i)*n2+j] = x2upper[i];
        }
      }

      n_bytes = n1*n2*sizeof(double);
      cudaMalloc((void**)&x2upper2dd,n_bytes);
      cudaMemcpy(x2upper2dd,tmp_cpu,n_bytes,cudaMemcpyHostToDevice);
      deviceVars["upper"] = x2upper2dd;
      delete[] tmp_cpu;
    }
  }//init

  void initProgvarBounds() {
    ChemtableData2D * c = getVar("prog");
    assert(x2lower==NULL); x2lower = new double[n1];
    assert(x2upper==NULL); x2upper = new double[n1];
    for (int i=0; i<n1; ++i) {
      x2lower[i] = (*c)(i,0);
      x2upper[i] = (*c)(i,n2-1);
    }
    return;
  }

  void initInterpWts() {
    assert(invDenom1==NULL); invDenom1 = new double[n1-3][4];
    assert(invDenom2==NULL); invDenom2 = new double[n2-3][4];
    for (int i=0; i<n1-3; ++i)
      cubicInterpDenom(invDenom1[i],&x1[i]);
    for (int i=0; i<n2-3; ++i)
      cubicInterpDenom(invDenom2[i],&x2[i]);
    return;
  }

  // ================================================
  // 2D piecewise cubic interpolation
  // ================================================

  void lookupSpecial(double* rout, const string& name,
                     const double * y1d, const double * y2ad,
                     const double * y2bd, const int n);

  void lookupSpecialFinish(double * rout, const string& name, const int n) {
    double * tmpd = cudaBuffers[name];
    cudaMemcpy(rout,tmpd, n*sizeof(double),cudaMemcpyDeviceToHost);
    cudaBuffers.erase(name);
    cudaFree(tmpd);
  }

  void lookupCubicPtrVector(vector<double*>& routPtrVec, const vector<string>& nameVec,
                            const double *y1, const double *y2,const int n);

  void loadVariables(const vector<string>& strVec) {
    for (int istr = 0; istr < strVec.size(); ++istr)
      ChemtableData2D * dummy = getVar(strVec[istr]);
  }

  // =======================
  // add a variable to table
  // =======================
  ChemtableData2D * getVar(const string& varname){

    // double check that the variable is not already in the list
    for (list<ChemtableData2D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
      if ( data->name == varname ) return(&(*data));

    // add the variable to the list
    // this needs to be done on all the processors

    datalist.push_back(ChemtableData2D(varname,n1,n2));
    ChemtableData2D * thisdata = &(datalist.back());
    assert(thisdata->name == varname);

    // read the data
    readData(thisdata);

    // return the data pointer
    return(thisdata);

  }

  void deleteVar(const string& varname){

    // double check that the variable is not already in the list
    for (list<ChemtableData2D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
      if ( data->name == varname ) {
	if ( mpi_rank == 0)
	  cout << " > deleting variable "<< data->name << " from chemtable" << endl;
	datalist.erase(data);
	return;
      }
  }

  void readHeaderAndCoordinates() {

    // read on all the ranks for now
    // open the file and go to data
    FILE * fp = NULL;
    // open the file
    const int file_err = MiscUtils::openFile(&fp,tableName,"rb");
    if (file_err != 0) throw(file_err);
    // fp = fopen(tableName.c_str(), "rb");
    // if (!fp) {
    //   CERR("cannot open chemtable file: " << tableName);
    // }

    // read the magic number and IO version and figure out the byte swap
    int byte_swap = 0;
    {
      int itmp[2];
      fread(itmp,sizeof(int),2,fp);
      if (itmp[0] != CHEM_IO_MAGIC_NUMBER) {
	ByteSwap::byteSwap(itmp,2);
	if (itmp[0] != CHEM_IO_MAGIC_NUMBER) {
	  cerr << ERRSTART << "Error: file " << tableName << " does not start as expected. aborting." << ERREND << endl;
	  throw(0);
	}
	byte_swap = 1;
      }
      if (itmp[1] != CHEM_IO_VERSION) {
	cerr << ERRSTART << "Error: in file " << tableName << "io version differs from current version: " << itmp[1] << ERREND << endl;
	throw(0);
      }
    }

    // initialize the offset...
    int8 offset = sizeof(int)*2;
    Header header;

    // find the data and read it
    int done  = 0;
    int ivar  = 0;
    int icoor = 0;
    while ( done!=1 ) {

      // go to the offset position
      my_fseek(fp,offset);
      // read the header
      fread(&header,sizeof(Header),1,fp);
      if (byte_swap)
	ByteSwap::byteSwapHeader(&header,1);

      // compare the header id
      if ( header.id == UGP_IO_CT_CART_2D ) {
	// type of the table
	tableType = header.name;
	// table dimensions
	n1    = header.idata[0]; assert(n1>1);
	n2    = header.idata[1]; assert(n2>1);
	nvars = header.idata[2];
	// allocate space for table coordinates
	assert(x1==NULL);      x1 = new double[n1];
	assert(x2==NULL);      x2 = new double[n2];
	assert(varName==NULL); varName = new string[nvars];
	// reference quantities
	pressure = header.rdata[0]; assert(pressure>0.0);
      }


      else if ( header.id == UGP_IO_CT_COOR ){
	double* tmpPtr;
	// copy the name to a string
	string name = header.name;
	if ( name == "Z" ) {
	  assert( header.idata[0] == n1 );
	  assert(x1!=NULL);
	  tmpPtr = x1;
	  ++icoor;
        }
	else if ( name == "C" ) {
	  assert( header.idata[0] == n2 );
	  assert(x2!=NULL);
	  tmpPtr = x2;
	  ++icoor;
	}
	else{
	  cerr << ERRSTART << "Error: unknown coordinate " << header.name << ERREND << endl;
	  throw(0);
	}
	// read the data
	size_t dummy= fread(&(tmpPtr[0]),sizeof(double),header.idata[0],fp);
	assert(dummy==header.idata[0]);
	// byte swap
	if (byte_swap)
	  ByteSwap::byteSwap(&(tmpPtr[0]),header.idata[0]);
      }

      else if ( header.id == UGP_IO_CT_DATA ) {
	varName[ivar] = header.name;
	++ivar;
      }

      else if ( header.id == UGP_IO_EOF ) {
	done = 1;
      }

      else {
	cerr << ERRSTART << "Error: could not find one of the two coordinates or the header " << ERREND << endl;
	throw(0);
      }

      // increase the offset
      offset += header.skip;

    } // while (done != 1)

    // close the file
    fclose(fp);

    assert(ivar==nvars);
    assert(icoor==2);

  }


  // ======================================
  // read the data associated with a header
  // ======================================
  void readData(ChemtableData2D * data) {

    // use this to broadcast an error to all processors
    int error_rank = 0;

    // read only on rank zero
    if ( mpi_rank == 0 ) {

      cout << " > reading " << data->name << " from chemtable " << tableName << ": ";

      // open the file and go to data
      FILE * fp = NULL;
      // open the file
      const int file_err = MiscUtils::openFile(&fp,tableName,"rb");
      // fp = fopen(tableName.c_str(), "rb");
      if (file_err != 0) {
        cerr << "cannot open chemtable file: " << tableName << endl;
        error_rank = 1;
      }


      // read the magic number and IO version and figure out the byte swap
      int byte_swap = 0;
      {
	int itmp[2];
	fread(itmp,sizeof(int),2,fp);
	if (itmp[0] != CHEM_IO_MAGIC_NUMBER) {
	  ByteSwap::byteSwap(itmp,2);
	  if (itmp[0] != CHEM_IO_MAGIC_NUMBER) {
	    cerr << ERRSTART << "Error: file " << tableName << " does not start as expected. aborting." << ERREND << endl;
	    error_rank = 1;
	  }
	  byte_swap = 1;
	}
	if (itmp[1] != CHEM_IO_VERSION) {
	  cerr << ERRSTART << "Error: in file " << tableName << "io version differs from current version: " << itmp[1] << ERREND << endl;
	  error_rank = 1;
	}
      }

      // initialize the offset...
      int8 offset = sizeof(int)*2;
      Header header;

      // find the data and read it
      int done = 0;
      while ( (done!=1) && (error_rank==0) ) {

	// go to the offset position
	my_fseek(fp,offset);
	// read the header
	fread(&header,sizeof(Header),1,fp);
	if (byte_swap)
	  ByteSwap::byteSwapHeader(&header,1);
	// check we are not at the end of file
	if ( header.id == UGP_IO_EOF ) {
	  cerr << ERRSTART << "Error: could not find variable " << data->name << " in the chemtable" << ERREND << endl;
	  error_rank = 1;
	}
	// compare the names
	string name = header.name;
	if ( name == data->name ) {
	  // check compatibility
	  assert( header.id       == UGP_IO_CT_DATA );
	  assert( header.idata[0] == (n1*n2) );
	  assert( header.idata[1] == n1 );
	  assert( header.idata[2] == n2 );
	  // read the data
	  size_t dummy= fread(&((*data)(0,0)),sizeof(double),header.idata[0],fp);
          assert(dummy==header.idata[0]);
	  // byte swap
	  if (byte_swap)
	    ByteSwap::byteSwap(&((*data)(0,0)),header.idata[0]);
	  // dump variable range
	  double varmin = (*data)(0,0);
	  double varmax = (*data)(0,0);
	  for ( int i=0; i<n1; ++i)
	    for ( int j=0; j<n2; ++j){
  	      varmax = max(varmax,(*data)(i,j));
	      varmin = min(varmin,(*data)(i,j));
	    }
	  cout << varmin << " < "<< data->name <<" < " << varmax << endl ;
	  // set the done
	  done = 1;
	}

	// increase the offset
	offset += header.skip;
      } // while (done != 1)

      // close the file
      if ( error_rank == 0 )
	fclose(fp);
    } // if ( mpi_rank == 0 )

    // send the error check from rank zero to all the processors
    MPI_Bcast(&error_rank,1,MPI_INT,0,mpi_comm);

    // if we have an error
    if ( error_rank == 1 )
      CERR("Problem reading file " << tableName);

    // send the data to all other ranks ...
    // here I assume that the * data is
    // already set properly on all the
    // processors

    MPI_Bcast(&((*data)(0,0)),data->getSize(),MPI_DOUBLE,0,mpi_comm);
  }
};

void initChemtableGpu(AbstractChemtable2D* &chemtable, const string& tablename);
void deleteChemtableGpu(AbstractChemtable2D* &chemtable);
#endif
