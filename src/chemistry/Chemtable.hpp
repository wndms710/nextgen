#ifndef CHEMTABLE_HPP
#define CHEMTABLE_HPP

#include "AbstractChemtable.hpp"

template <int DIM_K>
class KdChemtable {
private:
  class Verticies {
  private:
    double * x_array;
    double * data_array;

  public:
    int nve;
    int data_m;

    double x_min[DIM_K];
    double x_max[DIM_K];
    double * data_min;
    double * data_max;

    Verticies(const int nve) {
      COUT3("Verticies()");
      assert(nve > 0);

      this->nve = nve;
      x_array = new double[nve*DIM_K];

      data_m = 0;
      data_array = NULL;

      data_min = NULL;
      data_max = NULL;
    }

    ~Verticies() {
      COUT3("~Verticies()");
      if (x_array != NULL)    delete[] x_array;
      if (data_array != NULL) delete[] data_array;

      if (data_min != NULL) delete[] data_min;
      if (data_max != NULL) delete[] data_max;
    }

    void resize_data(const int data_m_new) {
      assert(data_m_new > data_m);

      if (data_array == NULL) {
	data_array = new double[nve*data_m_new];
      }
      else {
	double * data_copy = new double[nve*data_m];
	for (int i = 0; i < nve*data_m; ++i)
	  data_copy[i] = data_array[i];

	delete[] data_array;
	data_array = new double[nve*data_m_new];

	for (int ive = 0; ive < nve; ++ive) {
	  for (int idata = 0; idata < data_m; ++idata) {
	    data_array[ive*data_m_new + idata] = data_copy[ive*data_m + idata];
	  }
	}

	delete[] data_copy;
      }

      data_m = data_m_new;
    }

    void calcXMinMax() {
      for (int idim = 0; idim < DIM_K; ++idim) {
	x_max[idim] = -1e20;
	x_min[idim] = 1e20;
	for (int ive = 0; ive < nve; ++ive) {
	  x_max[idim] = max(x_max[idim],x(ive,idim));
	  x_min[idim] = min(x_min[idim],x(ive,idim));
	}
	assert(x_max[idim] > x_min[idim]);
      }
    }

    void calcDataMinMax() {
      if (data_min != NULL || data_max != NULL) {
	assert(data_min != NULL && data_max != NULL);
	delete[] data_min;
	delete[] data_max;
      }

      data_min = new double[data_m];
      data_max = new double[data_m];

      for (int idata = 0; idata < data_m; ++idata) {
	data_min[idata] = 1e20;
	data_max[idata] = -1e20;
	for (int ive = 0; ive < nve; ++ive) {
	  data_min[idata] = min(data_min[idata],data(ive,idata));
	  data_max[idata] = max(data_max[idata],data(ive,idata));
	}
	assert(data_max[idata] > data_min[idata]);
      }
    }

    inline double& x(const int ive, const int idim) {
#ifdef CTI_DEBUG
      assert((ive >= 0) && (ive < nve));
      assert((idim >= 0) && (idim < DIM_K));
#endif
      return(x_array[ive*DIM_K + idim]);
    }

    inline double& data(const int ive, const int idata) {
#ifdef CTI_DEBUG
      assert((ive >= 0) && (ive < nve));
      assert((idata >= 0) && (idata < data_m));
#endif
      return(data_array[ive*data_m + idata]);
    }

  };

  class HyperCubes {
  private:
    int * veohc_array;

  public:
    int nhc;

    HyperCubes(const int nhc) {
      COUT3("HyperCubes()");
      assert(nhc > 0); this->nhc = nhc;

      veohc_array = new int[nhc*(1<<DIM_K)];
    }

    ~HyperCubes() {
      COUT3("~HyperCubes()");

      if (veohc_array != NULL) delete[] veohc_array;
    }

    inline int& veohc(const int ihc, const int ive) {
#ifdef CTI_DEBUG
      assert((ihc >= 0) && (ihc < nhc));
      assert((ive >= 0) && (ive < (1<<DIM_K)));
#endif
      return(veohc_array[ihc*(1<<DIM_K) + ive]);
    }

  };

  struct Node {
    int dim_split;
    double x_split;
    Node * left;
    Node * right;
    int hc;
  };

protected:

  KdChemtableParallelIO * file;

  Verticies * ve;
  HyperCubes * hc;

  int nnode;
  Node * nodes;
  Node * root_node;

  vector<string> xNameVec;
  vector<string> dataNameVec;

public:

  KdChemtable() {
    COUT2("KdChemtable()");

    ve    = NULL;
    hc    = NULL;
    nodes = NULL;

  }

  virtual ~KdChemtable() {
    COUT2("~KdChemtable()");

    if (ve != NULL)    delete ve;
    if (hc != NULL)    delete hc;
    if (nodes != NULL) delete[] nodes;

    if (file != NULL) delete file;
  }


  // init gets called at derived class...
  void init(const string& filename, const string& tabletype) {
    // make sure coordinate names have been set...
    assert(xNameVec.size() == DIM_K);

    COUT1(" > loading kd table: " << filename << "... ");
    file = new KdChemtableParallelIO(filename);

    Header header;
    file->getHeaderFromFileAndBCast(header, tabletype, UGP_IO_CT_KD);
    int file_dim_k  = header.idata[0];
    if (file_dim_k != DIM_K)
      CERR("mismatch between declared dimension of table (" << DIM_K <<
	   ") and dimension in chemtable file (" << file_dim_k << ") ");

    if (mpi_rank == 0)
      assert(header.idata[1] == int( (file->serial_ptr)->dataNameVec.size()));

    ve = new Verticies(header.idata[2]);
    hc = new HyperCubes(header.idata[3]);

    nnode = header.idata[4];

    extraMainHeader(header);

    loadXCoords();

    loadHyperCubes();
    loadTree();
  }


  virtual void extraMainHeader(Header& header) { }

  void loadXCoords() {
    assert(ve != NULL);

    double * tmp_x = new double[ve->nve];
    for (int idim = 0; idim < DIM_K; ++idim) {

      file->getDataFromFileAndBCast(tmp_x, ve->nve, xNameVec[idim], UGP_IO_CT_COOR, MPI_DOUBLE);

      for (int ive = 0; ive < ve->nve; ++ive)
	ve->x(ive,idim) = tmp_x[ive];
    }
    delete[] tmp_x;

    ve->calcXMinMax();
  }

  void loadHyperCubes() {
    COUT2(" > loading hypercubes... ");
    int * tmp_veohc = new int[hc->nhc*(1<<DIM_K)];

    file->getDataFromFileAndBCast(tmp_veohc, hc->nhc*(1<<DIM_K), "VEOHC", UGP_IO_CT_KD_HC, MPI_INT);

    for (int ihc = 0; ihc < hc->nhc; ++ihc)
      for (int ive = 0; ive < (1<<DIM_K); ++ive)
	hc->veohc(ihc,ive) = tmp_veohc[ihc*(1<<DIM_K) + ive];

    delete[] tmp_veohc;
  }

  void loadTree() {
    COUT2(" > loading binary tree... ");

    int * tmp_dim_split = new int[nnode];
    file->getDataFromFileAndBCast(tmp_dim_split, nnode, "DIM_SPLIT", UGP_IO_CT_KD_NODE, MPI_INT);
    nodes = new Node[nnode];
    for (int inode = 0; inode < nnode; ++inode)
      nodes[inode].dim_split = tmp_dim_split[inode];
    delete[] tmp_dim_split;

    double * tmp_x_split = new double[nnode];
    file->getDataFromFileAndBCast(tmp_x_split, nnode, "X_SPLIT", UGP_IO_CT_KD_NODE, MPI_DOUBLE);
    for (int inode = 0; inode < nnode; ++inode)
      nodes[inode].x_split = tmp_x_split[inode];
    delete[] tmp_x_split;

    int * tmp_left  = new int[nnode];
    file->getDataFromFileAndBCast(tmp_left, nnode, "LEFT", UGP_IO_CT_KD_NODE, MPI_INT);
    for (int inode = 0; inode < nnode; ++inode) {
      if (tmp_left[inode] == -1) {
	nodes[inode].left = NULL;
      }
      else {
	assert(tmp_left[inode] >= 0 && tmp_left[inode] < nnode);
	nodes[inode].left = &nodes[tmp_left[inode]];
      }
    }
    delete[] tmp_left;

    int * tmp_right = new int[nnode];
    file->getDataFromFileAndBCast(tmp_right, nnode, "RIGHT", UGP_IO_CT_KD_NODE, MPI_INT);
    for (int inode = 0; inode < nnode; ++inode) {
      if (tmp_right[inode] == -1) {
	nodes[inode].right = NULL;
      }
      else {
	assert(tmp_right[inode] >= 0 && tmp_right[inode] < nnode);
	nodes[inode].right = &nodes[tmp_right[inode]];
      }
    }
    delete[] tmp_right;

    int * tmp_hc = new int[nnode];
    file->getDataFromFileAndBCast(tmp_hc, nnode, "HC", UGP_IO_CT_KD_NODE, MPI_INT);
    for (int inode = 0; inode < nnode; ++inode)
      nodes[inode].hc = tmp_hc[inode];
    delete[] tmp_hc;

    root_node = &nodes[0];
  }


  void loadData(const vector<string>& varNameVec) {
    // checks to see if data associated with passed strings have been loaded
    // if not, loads...

    assert(ve != NULL);

    vector<string> newVarVec;
    for (unsigned int ivar = 0; ivar < varNameVec.size(); ++ivar) {
      bool found = false;
      for (int idata = 0; idata < ve->data_m; ++idata) {
	if (varNameVec[ivar] == dataNameVec[idata]) {
	  found = true;
	  break;
	}
      }

      // push to list if not found...
      if (!found) newVarVec.push_back(varNameVec[ivar]);
    }

    if (newVarVec.size() != 0) {

      const int data_m_old = ve->data_m;
      {
	// resize arrays...
	const int data_m_new = ve->data_m + newVarVec.size();
	ve->resize_data(data_m_new);
	assert(ve->data_m == data_m_new);
      }

      double * tmp_data = new double[ve->nve];
      for (unsigned int ivar = 0;  ivar < newVarVec.size(); ++ivar) {
	if (mpi_rank == 0) cout << " > loading " << newVarVec[ivar] << "... ";

	file->getDataFromFileAndBCast(tmp_data, ve->nve, newVarVec[ivar], UGP_IO_CT_DATA, MPI_DOUBLE);

	double dmin = 1e20;
	double dmax = -1e20;
	for (int ive = 0; ive < ve->nve; ++ive) {
	  ve->data(ive,data_m_old + ivar) = tmp_data[ive];
	  dmin = min(dmin, tmp_data[ive]);
	  dmax = max(dmax, tmp_data[ive]);
	}

	if (mpi_rank == 0) cout << "[ " << dmin << ", " << dmax << " ] " << endl;

	dataNameVec.push_back(newVarVec[ivar]);
      }
      delete[] tmp_data;

      ve->calcDataMinMax();

    }


  }


  void kdlookup(vector<double*>& varPtrVec, const vector<string>& varNameVec,
		vector<double*>& xPtrVec, const int npoints) {
    // NB: kdLookup will clip xPtrVec, so it should be a copy!

    loadData(varNameVec);
    int * data_flag = new int[ve->data_m];

    for (int idata = 0; idata < ve->data_m; ++idata) data_flag[idata] = -1;

    for (unsigned int ivar = 0; ivar < varNameVec.size(); ++ivar) {
      bool found = false;
      for (int idata = 0; idata < ve->data_m; ++idata) {
	if (varNameVec[ivar] == dataNameVec[idata]) {
	  data_flag[idata] = ivar;
	  assert(found == false); // make sure no duplicates
	  found = true;
	}
      }
      // make sure everything was found
      assert(found == true);
    }

    for (int ipoint = 0; ipoint < npoints; ++ipoint) {
      // clip coordinates
      for (int idim = 0; idim < DIM_K; ++idim)
	xPtrVec[idim][ipoint] = min(max(xPtrVec[idim][ipoint],ve->x_min[idim]),ve->x_max[idim]);

      // and initialize data
      for (unsigned int ivar = 0; ivar < varNameVec.size(); ++ivar)
	varPtrVec[ivar][ipoint] = 0.0;
    }

    for (int ipoint = 0; ipoint < npoints; ++ipoint) {

      Node * this_node = root_node;
      while (this_node->hc == -1) {
	if (xPtrVec[this_node->dim_split][ipoint] < this_node->x_split)
	  this_node = this_node->left;
	else
	  this_node = this_node->right;
      }
      const int ihc = this_node->hc;

#ifdef CTI_DEBUG
      // make sure we've actually found it...
      bool found = true;
      for (int ive = 0; ive < (1<<DIM_K); ++ive) {
	for (int idim = 0; idim < DIM_K; ++idim) {
	  if (ive&(1<<idim)) {
	    if (xPtrVec[idim][ipoint] > ve->x(hc->veohc(ihc,ive),idim)) {
	    found = false;
	    break;
	    }
	  }
	  else {
	    if (xPtrVec[idim][ipoint] < ve->x(hc->veohc(ihc,ive),idim)) {
	    found = false;
	    break;
	    }
	  }
	}
      }
      assert(found);
#endif

      double weights[1<<DIM_K];
      for (int ive = 0; ive < (1<<DIM_K); ++ive) weights[ive] = 1.0;

      for (int idim = 0; idim < DIM_K; ++idim) {
	const double inv_denom = 1.0/(ve->x(hc->veohc(ihc,(1<<DIM_K) - 1),idim) - ve->x(hc->veohc(ihc,0),idim));
	for (int ive = 0; ive < (1<<DIM_K); ++ive) {
	  const double tmp = (xPtrVec[idim][ipoint] - ve->x(hc->veohc(ihc,0),idim))*inv_denom;
	  if (ive&(1<<idim))
	    weights[ive] *= tmp;
	  else
	    weights[ive] *= (1.0 - tmp);
	}
      }

#ifdef CTI_DEBUG
      // make sure weights are convex combination...
      double weight_sum = 0.0;
      for (int ive = 0; ive < (1<<DIM_K); ++ive) {
	assert(weights[ive] >= 0.0);
	weight_sum += weights[ive];
      }
      assert(fabs(weight_sum - 1.0) < 1e-9);
#endif

      for (int ive = 0; ive < (1<<DIM_K); ++ive) {
	for (int idata = 0; idata < ve->data_m; ++idata) {
	  if (data_flag[idata] >= 0) {
	    varPtrVec[data_flag[idata]][ipoint] += weights[ive]*ve->data(hc->veohc(ihc,ive),idata);
	  }
	}
      }

    }

    delete[] data_flag;
  }

  void writeKdTableTecplot(const string& filename) {
    cout << " > writing table to the tecplot file: " << filename << endl;

    assert(hc != NULL);
    assert(ve != NULL);

    FILE* fp = fopen(filename.c_str(), "w");
    fprintf(fp,"TITLE = \"%s\"\n",filename.c_str());
    fprintf(fp,"VARIABLES = ");

    for(int idim = 0; idim < DIM_K; ++idim)
      fprintf(fp,"\"%s\"\n", xNameVec[idim].c_str());

    for(int idata = 0; idata < ve->data_m; ++idata)
      fprintf(fp,"\"%s\"\n", dataNameVec[idata].c_str());

    fprintf(fp,"ZONE T=\"%s\"\n",filename.c_str());
    if (DIM_K == 2)
      fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n",ve->nve,hc->nhc);
    else if (DIM_K == 3)
      fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=BRICK\n",ve->nve,hc->nhc);
    else
      CERR_S("TECPLOT output for dimension >= 4 is not supported yet");

    for (int ive = 0; ive < ve->nve; ++ive) {
      for (int idim = 0; idim < DIM_K; ++idim)
	fprintf(fp, "%g ", ve->x(ive,idim));

      for (int idata = 0; idata < ve->data_m; ++idata)
	fprintf(fp, "%g ", ve->data(ive,idata));

      fprintf(fp,"\n");
    }

    for (int ihc = 0; ihc < hc->nhc; ++ihc) {
      // 1-indexing for tecplot...
      if (DIM_K == 2)
	fprintf(fp, "%d %d %d %d\n", hc->veohc(ihc,0) + 1, hc->veohc(ihc,1) + 1,
		hc->veohc(ihc,3) + 1, hc->veohc(ihc,2) + 1);
      else if (DIM_K == 3)
	fprintf(fp, "%d %d %d %d %d %d %d %d\n",  hc->veohc(ihc,0) + 1, hc->veohc(ihc,1) + 1,
		hc->veohc(ihc,3) + 1, hc->veohc(ihc,2) + 1, hc->veohc(ihc,4) + 1,
		hc->veohc(ihc,5) + 1, hc->veohc(ihc,7) + 1, hc->veohc(ihc,6) + 1);
      else
	CERR_S("TECPLOT output for dimension >= 4 is not supported yet");
    }

    fclose(fp);

  }

};


class KdChemtable3D :  protected KdChemtable<3>, public AbstractChemtable3D {

public:

  KdChemtable3D(const string& filename, const string& tabletype) {
    COUT2("KdChemtable3D()");

    assert(xNameVec.size() == 0);
    this->xNameVec.push_back("Z");
    this->xNameVec.push_back("Zvar_norm");
    this->xNameVec.push_back("C");

    this->tableName = filename;
    this->tableType = tabletype;

    init(filename, tabletype);

    info();
  }

  ~KdChemtable3D() {
    COUT2("~KdChemtable3D()");

    info();
  }

  bool varExists(const string& var_name) const { 
    
    return (std::find(dataNameVec.begin(),dataNameVec.end(),var_name) != dataNameVec.end());

  }

  void loadVariables(const vector<string>& strVec) {
    if (strVec.size() == 0)
      CERR("pointless loadVariables()");

    loadData(strVec);

    // XXXXXXX this should go somewhere else...
    if (mpi_rank == 0) cout << " > clipping lower bound of src_prog to prevent numerical ignition..." << endl;
    int isrc_prog = - 1;
    int iprog = -1;
    for (int idata = 0; idata < ve->data_m; ++idata) {
      if (dataNameVec[idata] == "prog") {
	assert(iprog == -1);
	iprog = idata;
      }
      else if (dataNameVec[idata] == "src_prog") {
	assert(isrc_prog == -1);
	isrc_prog = idata;
      }
    }
    assert(iprog >= 0);
    assert(isrc_prog >= 0);

    for (int ive = 0; ive < ve->nve; ++ive) {
      const double this_c = ve->x(ive,2);
      if (this_c <= 1e-6 || this_c < ve->data(ive, iprog) || ve->data(ive, isrc_prog) < 0.0) {
	ve->data(ive,isrc_prog) = 0.0;
      }
    }
  }

  // get some derived class-specific data from header
  void extraMainHeader(Header& header) {
    assert(header.name == tableType);
    pressure = header.rdata[1];
  }

  void lookupDataVecReduced(vector<double*>& routPtrVec, const vector<string>& nameVec, 
                            const double * y1, const int n) {


    // this takes that actual value of Z, a fake large value of C, and no turbulence variance.

    double * C_tmp    = new double[n]; 
    double * Zvar_tmp = new double[n];

    for (int i = 0; i < n; ++i) { 

      C_tmp[i]    = 10.0; // assuming C is weighted sum of mass fractions this is unphysically large
      Zvar_tmp[i] = 0.0;

    }

    lookupDataVec(routPtrVec,nameVec,y1,C_tmp,Zvar_tmp,n); 

    delete[] C_tmp;
    delete[] Zvar_tmp;

  }


  void lookupDataVec(vector<double*>& routPtrVec, const vector<string>& nameVec,
		     const double* x0, const double* x1, const double* x2, const int n) {

    // make a copy so we can modify...
    double * x0_copy = new double[n];
    double * x1_copy = new double[n];
    double * x2_copy = new double[n];

    for (int i = 0; i < n; ++i) {
      x0_copy[i] = x0[i];
      // in the table, x1 is a normalized zvar
      // hence, we need to normalize the given zvar...
      const double denom = x0[i]*(1.0 - x0[i]);
      if (denom != 0.0)
	x1_copy[i] = x1[i]/denom;
      else
	x1_copy[i] = 0.0;
      x2_copy[i] = x2[i];
    }

    vector<double*> rinPtrVec(3);
    rinPtrVec[0] = x0_copy;
    rinPtrVec[1] = x1_copy;
    rinPtrVec[2] = x2_copy;

    kdlookup(routPtrVec, nameVec, rinPtrVec, n);

    delete[] x0_copy;
    delete[] x1_copy;
    delete[] x2_copy;
  }

  // dummies!
  void info() {

    if ( mpi_rank == 0 ) {
      cout << "============== chemistry table information ===============" << endl;
      cout << " chemistry table information: " << endl;
      cout << "  > name:                     " << tableName << endl;
      cout << "  > type:                     " << tableType << endl;
      cout << "  > reference pressure:       " << pressure << endl;
      cout << "  > # of leaves:              " << hc->nhc << endl;
      cout << "  > coordinates and ranges:          " << endl;
      for (int idim = 0; idim < 3; ++idim)
	cout << "   > " << xNameVec[idim]  << ": "<< ve->x_min[idim] << " < " << xNameVec[idim] << " < " << ve->x_max[idim] << endl ;
      cout << "  > variables("<< (file->serial_ptr)->dataNameVec.size() << "):"<< endl;
      cout << "     ";
      for (uint var = 0; var < (file->serial_ptr)->dataNameVec.size(); ++var)
	cout << (file->serial_ptr)->dataNameVec[var] << " ";
      cout << endl;

      if ( ve->data_m == 0 ) {
	cout << "  > no variables are currently loaded " <<endl;
      }
      else {
	cout << "  > the following variables are currently loaded: " << endl;
	for (uint idata = 0; idata < dataNameVec.size(); ++idata) {
	  cout << "   > " << dataNameVec[idata]  << ": "<< ve->data_min[idata] << " < " << dataNameVec[idata] << " < " << ve->data_max[idata] << endl ;
	}
      }
    }

    if ( mpi_rank == 0 )
      cout << "=========== end of chemistry table information ============" << endl;

  }


  void writeTecplot(const string& filename) {
    writeKdTableTecplot(filename);
  }

  void writeTecplot2(const string& filename) {
    CERR("writeTecplot2 not implemented for kdtrees");
  }

  int getTableVarRange(double& varmin, double& varmax, const string& name) {
    for (int idata = 0; idata < ve->data_m; ++idata) {
      if (dataNameVec[idata] == name) {
	varmin = ve->data_min[idata];
	varmax = ve->data_max[idata];
	return(1);
      }
    }

    COUT2("***WARNING in getTableVarRange: the variable "
	  << name << " is not loaded yet/n" << " ...setting the range to 0:0 ");

    varmin = 0.0;
    varmax = 0.0;
    return(0);
  }

};

// =========================
// Cartesian chemtable class
// =========================
class CartesianChemtable1D: public AbstractChemtable1D {

private:

  double laminarFlameSpeed;
  double laminarFlameThickness;

  int n1, nvars;
  double* x1;

  double Cmin, Cmax;
  double dC;

  string *varName;

  list<ChemtableData1D> datalist;


public:

  // constructor
  CartesianChemtable1D(const string& name){
    if ( mpi_rank == 0 )
      cout << "CartesianChemtable1D():" << endl;

    // table dimensions
    n1    = 0;
    nvars = 0;

    // coordinates
    x1 = NULL;

    // array of variable names
    varName = NULL;

    // initialize table
    this->tableName =  name;
    init();

  }

  // deconstructor
  ~CartesianChemtable1D(){
    if ( mpi_rank == 0 )
      cout << "~CartesianChemtable1D():" << endl;

    //testLookup();

    //writeTecplot("test.dat");

    info();

    if (x1!=NULL) delete[] x1;
    if (varName!=NULL) delete[] varName;

  }

  virtual double getFlameThickness(){
    return(laminarFlameThickness);
  }


  virtual double getLaminarFlameSpeed() {
    return(laminarFlameSpeed);
  }

  // ================
  // dump table info
  // ================
  void info(){

    if ( mpi_rank == 0 ) {
      cout << "============== chemistry table information ===============" << endl;
      cout << " chemistry table information: " << endl;
      cout << "  > name:                     " << tableName << endl;
      cout << "  > type:                     " << tableType << endl;
      cout << "  > reference pressure:       " << pressure << endl;
      cout << "  > dimension:                ("<< n1 << ")" << endl;
      cout << "   > coordinate range:         ("<< x1[0] <<"," << x1[n1-1] << ")" << endl;
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

      for (list<ChemtableData1D>::iterator data = datalist.begin(); data != datalist.end(); ++data) {
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


  // =======================
  // find the variable range
  // =======================
  int getTableVarRange(double& varmin, double& varmax, const string& name){
    // initialze
    varmin = 0.0;
    varmax = 0.0;
    // loop through the list
    for (list<ChemtableData1D>::iterator data = datalist.begin(); data != datalist.end(); ++data) {
      if ( data->name == name ) {
	varmin = (*data)(0);
	varmax = (*data)(0);
	for ( int i=0; i<n1; ++i){
	  varmax = max(varmax,(*data)(i));
	  varmin = min(varmin,(*data)(i));
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
    ChemtableData1D * rho = getVar("rho");
    ChemtableData1D * src = getVar("src_prog");
    Smin =  1.0e+20;
    Smax = -1.0e+20;
    for (int i=0; i<n1; ++i) {
      const double rhoSrc = (*rho)(i) * (*src)(i);
      Smin = min(Smin,rhoSrc);
      Smax = max(Smax,rhoSrc);
    }
    return;
  }

  // ==========================================
  // write the table with the current variables
  // the one that are loaded to a tecplot file
  // ==========================================
  void writeTecplot(const string& filename) {

    // load all the variables
    for (int i=0; i<nvars; ++i) getVar(varName[i]);

    // only rank 0 needs to write
    if (mpi_rank == 0) {

      cout << " > writing 1d cartesian chemtable " << tableName << " to the Tecplot file "<< filename << endl;
      cout << "  > the following variables are currently loaded: " << endl;
      for (list<ChemtableData1D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
	cout << "   > " << data->name << endl ;

      FILE * fp = fopen(filename.c_str(),"w");
      fprintf(fp,"title = \"chemtable\"\n");
      fprintf(fp,"variables = \"C\"\n");
      for (list<ChemtableData1D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
	fprintf(fp,"\"%s\"\n",data->name.c_str());
      fprintf(fp,"zone i = %d, DATAPACKING=POINT\n",n1);
      for (int i=0; i <n1; ++i){
	fprintf(fp,"%lf",x1[i]);
	for (list<ChemtableData1D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
	  fprintf(fp," %lf",(*data)(i));
	fprintf(fp,"\n");
      }
      fclose(fp);
    }
  }


private:

  // ===================
  // Linear table lookup
  // ===================

  void lookupDataVec(vector<double*>& routPtrVec, const vector<string>& nameVec,
		     const double* x1, const int n) {
    lookupLinearPtrVector(routPtrVec,nameVec,x1,n);
  }



  // ====================================
  // initialize the table
  // Here we should initialize the types,
  // dimensions, and common variables
  // no major read at this stage
  // ====================================
  void init() {

    //if ( mpi_rank == 0 )
    //  cout << "initialize chemtable: "<< endl;

    // note: open the table and read it on all
    // processors. This is inefficient but since
    // a very small amount of data is involved
    // at this stage, it is not going to make
    // a big difference...


    // read header and coordinates
    readHeaderAndCoordinates();

    // dump table info
    info();

  }



  // ===========================
  // look up vector of variables
  // ===========================
  void lookupLinearPtrVector(vector<double*>& routPtrVec, const vector<string>& nameVec,
		       const double* x, const int n) {

    vector<ChemtableData1D *> dataPtrVec(nameVec.size());
    for (uint iname = 0; iname < nameVec.size(); ++iname)
      dataPtrVec[iname] = getVar(nameVec[iname]);


    double dx = x1[1] - x1[0];
    assert(fabs(dx-dC)<1.0e-14);
    for ( int i=0; i<n; ++i ) {
      int lb;
      double w;
      if ( x[i] <= x1[0] ) {
	lb = 0;
	w  = 1.0;
      }
      else if ( x[i] >= x1[n1-1] ) {
	lb = n1-2;
	w  = 0.0;
      }
      else {
	lb = int((x[i]-x1[0])/dx);
	assert(lb>=0);
	assert(lb<n1);
	lb = min(lb,(n1-2));
	assert(x[i]>=x1[lb]);
	assert(x[i]<=x1[lb+1]);
	assert(x1[lb+1]>x1[lb]);
	w = (x1[lb+1]-x[i])/(x1[lb+1]-x1[lb]);
	assert(w>=0.0);
	assert(w<=1.0);
      }
      for (uint j=0; j<routPtrVec.size(); ++j) {
	routPtrVec[j][i] =
	  w * (*(dataPtrVec[j]))(lb)     +
	  (1.0-w) * (*(dataPtrVec[j]))(lb+1);
      }

    }

  }

  void loadVariables(const vector<string>& strVec) {
    if (strVec.size() == 0)
      CERR("pointless loadVariables()");

    for (uint istr = 0; istr < strVec.size(); ++istr)
      getVar(strVec[istr]);
  }

  // =======================
  // add a variable to table
  // =======================
  ChemtableData1D * getVar(const string& varname){

    // double check that the variable is not already in the list
    for (list<ChemtableData1D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
      if ( data->name == varname ) return(&(*data));

    // add the variable to the list
    // this needs to be done on all the processors

    datalist.push_back(ChemtableData1D(varname,n1));
    ChemtableData1D * thisdata = &(datalist.back());
    assert(thisdata->name == varname);

    // read the data
    readData(thisdata);

    // return the data pointer
    return(thisdata);

  }


  // ============================
  // delete a variable from table
  // ============================
  void deleteVar(const string& varname){

    // double check that the variable is not already in the list
    for (list<ChemtableData1D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
      if ( data->name == varname ) {
	if ( mpi_rank == 0)
	  cout << " > deleting variable "<< data->name << " from chemtable" << endl;
	datalist.erase(data);
	return;
      }

  }



  // ===========================
  // read header and coordinates
  // ===========================
  void readHeaderAndCoordinates() {


    // read on all the ranks for now

    // open the file and go to data
    FILE * fp = NULL;
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
      if ( header.id == UGP_IO_CT_CART_1D ) {
	// type of the table
	tableType = header.name;
	// table dimensions
	n1    = header.idata[0]; assert(n1>1);
	nvars = header.idata[1];
	// allocate space for table coordinates
	assert(x1==NULL);      x1 = new double[n1];
	assert(varName==NULL); varName = new string[nvars];
	// reference quantities
	pressure              = header.rdata[0]; assert(pressure>0.0);
	laminarFlameSpeed     = header.rdata[1];
	laminarFlameThickness = header.rdata[2];
      }


      else if ( header.id == UGP_IO_CT_COOR ){
	double* tmpPtr;
	// copy the name to a string
	string name = header.name;
	if ( name == "C" ) {
	  assert( header.idata[0] == n1 );
	  assert(x1!=NULL);
	  tmpPtr = x1;
	  ++icoor;
	  dC =  header.rdata[0];
	  Cmin = header.rdata[1];
	  Cmax = header.rdata[2];
	}
	else{
	  cerr << ERRSTART << "Error: unknown coordinate " << header.name << ERREND << endl;
	  throw(0);
	}
	// read the data
	size_t dummy= fread(&(tmpPtr[0]),sizeof(double),header.idata[0],fp);
	assert(int(dummy)==header.idata[0]);
	// byte swap
	if (byte_swap)
	  ByteSwap::byteSwap(&(tmpPtr[0]),header.idata[0]);
	assert(fabs(x1[0]-Cmin)<1.0e-14);
	assert(fabs(x1[n1-1]-Cmax)<1.0e-14);

      }

      else if ( header.id == UGP_IO_CT_DATA ) {
	varName[ivar] = header.name;
	++ivar;
      }

      else if ( header.id == UGP_IO_EOF ) {
	done = 1;
      }

      else {
	cerr << ERRSTART << "Error: could not find one of the three coordinates or the header " << ERREND << endl;
	throw(0);
      }

      // increase the offset
      offset += header.skip;

    } // while (done != 1)

    // close the file
    fclose(fp);

    assert(ivar==nvars);
    assert(icoor==1);

  }

  // ======================================
  // read the data associated with a header
  // ======================================
  void readData(ChemtableData1D * data) {

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
	  assert( header.idata[0] == n1 );
	  // read the data
	  size_t dummy= fread(&((*data)(0)),sizeof(double),header.idata[0],fp);
	  // byte swap
	  if (byte_swap)
	    ByteSwap::byteSwap(&((*data)(0)),header.idata[0]);
	  assert(int(dummy)==header.idata[0]);
	  // dump variable range
	  double varmin = (*data)(0);
	  double varmax = (*data)(0);
	  for ( int i=0; i<n1; ++i){
	    varmax = max(varmax,(*data)(i));
	    varmin = min(varmin,(*data)(i));
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

    MPI_Bcast(&((*data)(0)),data->getSize(),MPI_DOUBLE,0,mpi_comm);


  }




public:


  // =================
  // test table lookup
  // =================
  void testLookup() {

    // Do it on rank zero only



    if ( mpi_rank == 0 )
      cout << "testing table lookup ..." << endl;
    // pick two variables in the table
    ChemtableData1D * var1 = getVar("T");
    if ( mpi_rank == 0 )
      cout << " > replacing T with linear function 2*C1 "<< endl;
    ChemtableData1D * var2 = getVar("rho");
    if ( mpi_rank == 0 )
      cout << " > replacing rho with linear function 20.1*C1"<< endl;
    // set them with linear functions
    for ( int i=0; i<n1; ++i) {
      (*var1)(i) = 2.0*x1[i];
      (*var2)(i) = 20.1*x1[i];
    }

    // Build a random array
    int myN = 100000;
    double *myC1 = new double[myN];
    //srand(10);
    for ( int i=0; i<myN; ++i ) {
      myC1[i] = x1[0] + (x1[n1-1]-x1[0])*(double)(rand())/(double)(RAND_MAX);
        }

    // result vector
    double *rout1 = new double[myN];
    double *rout2 = new double[myN];

    // do lookup
    lookup(rout1, "T", myC1,myN);
    lookup(rout2, "rho", myC1,myN);

    // compute error
    double myerror1 = 0.0;
    for ( int i=0; i<myN; ++i)
      myerror1 = max(fabs(rout1[i]-(2.0*myC1[i])),myerror1);
    double myerror2 = 0.0;
    for ( int i=0; i<myN; ++i)
      myerror2 = max(fabs(rout2[i]-(20.1*myC1[i])),myerror2);
    // take the max across processors
    double error1;
    MPI_Reduce(&myerror1,&error1,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    double error2;
    MPI_Reduce(&myerror2,&error2,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    if ( mpi_rank == 0 )
      cout << "  >> error in lookups of the two linear functions (should be zero) = " << error1 << " , " << error2 << endl;


    // deletle variables from table
    deleteVar("T");
    deleteVar("rho");

    // delete memory
    delete[] myC1;
    delete[] rout1;
    delete[] rout2;

    if ( mpi_rank == 0 )
      cout << "done with table lookup test" << endl;

  }


};

// =========================
// Cartesian chemtable class
// =========================
class CartesianChemtable2D: public AbstractChemtable2D {

private:

  string tableName;
  //string tableType;

  int n1, n2, nvars;
  double *x1, *x2;
  double *x2lower, *x2upper;
  double (*invDenom1)[4], (*invDenom2)[4];
  MiscUtils::indexMap *idxMap1, *idxMap2;

  string *varName;
  set<string> sorted_vars;

  list<ChemtableData2D> datalist;
  typedef std::pair<int,double> myPair;
  typedef std::pair<int,int> myPairInt;
  double SMALL_CLIP;

  double x2tol;

public:

  // constructor
  CartesianChemtable2D(const string& name){
    if ( mpi_rank == 0 )
      cout << "CartesianChemtable2D():" << endl;

    // table dimensions
    n1    = 0;
    n2    = 0;
    nvars = 0;

    // coordinates
    x1 = NULL;
    x2 = NULL;

    x2lower = NULL;
    x2upper = NULL;

    invDenom1 = NULL;
    invDenom2 = NULL;

    idxMap1 = NULL;
    idxMap2 = NULL;

    // array of variable names
    varName = NULL;

    // small value used for clipping
    SMALL_CLIP = 1.01e-6;

    // initialize table
    this->tableName =  name;
    //this->tableType =  type;
    init();

    // test table lookup
    //testLookup();

  }

  ~CartesianChemtable2D(){
    if ( mpi_rank == 0 )
      cout << "~CartesianChemtable2D()" << endl;

    if (x1!=NULL) delete[] x1;
    if (x2!=NULL) delete[] x2;
    if (x2lower!=NULL) delete[] x2lower;
    if (x2upper!=NULL) delete[] x2upper;
    if (varName!=NULL) delete[] varName;
    if (invDenom1!=NULL) delete[] invDenom1;
    if (invDenom2!=NULL) delete[] invDenom2;
    if (idxMap1!=NULL) delete idxMap1;
    if (idxMap2!=NULL) delete idxMap2;
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


  // =======================
  // find the variable range
  // =======================
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
    for (int i=0; i<nvars; ++i) getVar(varName[i]);

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

  // ===================
  // Linear table lookup
  // ===================
  void lookupDataVec(vector<double*>& routPtrVec, const vector<string>& nameVec,
		     const double* x1, const double* x2, const int n) {
    //lookupLinearPtrVector(routPtrVec,nameVec,x1,x2,n);
    lookupCubicPtrVector(routPtrVec,nameVec,x1,x2,n);
  }

  void lookupDataVecReducedNew(vector<double*>& routPtrVec, const vector<string>& nameVec,
                            const double *__restrict__ y1, const int (*__restrict__ cvofa)[2], const int n) {

    // this only depends on the first dimension
    vector<ChemtableData2D*> dataPtrVec(nameVec.size());
    for (unsigned int iname = 0; iname < nameVec.size(); ++iname) {
      if ( (nameVec[iname] == "UPPER")||(nameVec[iname] == "upper"))
        dataPtrVec[iname] = NULL; // this is a special keyword..
      else
        dataPtrVec[iname] = getVar(nameVec[iname]);
    }//iname ...

    for (int i =0; i < n ; ++i) {

      const int icv0   = cvofa[i][0];
      const int icv1   = cvofa[i][1];
      const double y1_ = 0.5*(y1[icv0] + y1[icv1]);

      // find the interval along the first coord ...
      double Y1    = min(x1[n1-1],max(x1[0],y1_));
      const int i1 = min(n1-4,max(0,idxMap1->findIndex(Y1)-1));

      double w1[4];
      //MiscUtils::cubicInterpWts(w1,&x1[i1],Y1);
      MiscUtils::cubicInterpWts(w1,invDenom1[i1],&x1[i1],Y1);

      // we're ready to go; the reduced data vec isnt
      // a function of the second coordinate ...
      for (unsigned int j =0; j < routPtrVec.size(); ++j) {

        double val = 0.0;
        if ( dataPtrVec[j] == NULL ) {
          for (int k =0; k < 4; ++k)
            val += w1[k]*x2upper[i1+k];
        } else {
          double * this_data = dataPtrVec[j]->data;
          const int n2_ = dataPtrVec[j]->n2;
          // the data vector is 2d but it doesn't depend on the
          // second coord, so we'll use idx 0.  going forward,
          // these data values should be precomputed in 1d...
          for (int k =0 ; k < 4 ; ++k) {
            const int offset = (i1+k)*n2_;
            val += w1[k]*this_data[offset];
          }
        }
        routPtrVec[j][i] = val;
      }//j
    }//i
  }//lookupDataVecReducedNew

  void lookupDataVecReduced(vector<double*>& routPtrVec, const vector<string>& nameVec,
                            const double * y1, const int n) {

    // this only depends on the first dimension
    vector<ChemtableData2D*> dataPtrVec(nameVec.size());
    for (unsigned int iname = 0; iname < nameVec.size(); ++iname) {
      if ( (nameVec[iname] == "UPPER")||(nameVec[iname] == "upper"))
        dataPtrVec[iname] = NULL; // this is a special keyword..
      else
        dataPtrVec[iname] = getVar(nameVec[iname]);
    }//iname ...

    for (int i =0; i < n ; ++i) {

      // find the interval along the first coord ...
      double Y1    = min(x1[n1-1],max(x1[0],y1[i]));
      const int i1 = min(n1-4,max(0,idxMap1->findIndex(Y1)-1));

      double w1[4];
      //MiscUtils::cubicInterpWts(w1,&x1[i1],Y1);
      MiscUtils::cubicInterpWts(w1,invDenom1[i1],&x1[i1],Y1);

      // we're ready to go; the reduced data vec isnt
      // a function of the second coordinate ...
      for (unsigned int j =0; j < routPtrVec.size(); ++j) {

        double val = 0.0;
        if ( dataPtrVec[j] == NULL ) {
          for (int k =0; k < 4; ++k)
            val += w1[k]*x2upper[i1+k];
        } else {
          double * this_data = dataPtrVec[j]->data;
          const int n2_ = dataPtrVec[j]->n2;
          // the data vector is 2d but it doesn't depend on the
          // second coord, so we'll use idx 0.  going forward,
          // these data values should be precomputed in 1d...
          for (int k =0 ; k < 4 ; ++k) {
            const int offset = (i1+k)*n2_;
            val += w1[k]*this_data[offset];
          }
        }
        routPtrVec[j][i] = val;
      }//j
    }//i
  }//lookupDataVecReduced...

  // ====================================
  // initialize the table
  // Here we should initialize the types,
  // dimensions, and common variables
  // no major read at this stage
  // ====================================
  void init() {

    //if ( mpi_rank == 0 )
    //  cout << "initialize chemtable: "<< endl;

    // note: open the table and read it on all
    // processors. This is inefficient but since
    // a very small amount of data is involved
    // at this stage, it is not going to make
    // a big difference...


    // read header and coordinates
    readHeaderAndCoordinates();

    // dump table info
    info();

    // get progvar upper/lower bounds for normalization
    initProgvarBounds();

    // precompute cubic weight denominator
    initInterpWts();

    // initialize fast index-finding
    idxMap1 = new MiscUtils::indexMap(x1,n1);
    idxMap2 = new MiscUtils::indexMap(x2,n2);

    // compute tolerance for fast integral lookups
    initX2Tol();
  }

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

  void initX2Tol() {

    // set max allowable error for integral lookups... i.e. I_err/I_max < errmax
    const double errmax = getDoubleParam("CHEM_LOOKUP_TOL",0.01);

    assert(x2lower!=NULL);
    assert(x2upper!=NULL);

    ChemtableData2D * I = getVar("int_rho_src");
    int iz = -1;
    double I_max = -1e+20;
    for (int i=0; i<n1; ++i)
      for (int j=0; j<n2; ++j)
        if (abs((*I)(i,j)) > I_max) { I_max = abs((*I)(i,j)); iz = i; }

    assert( iz >= 0);

    double dIdc_max = -1e+20;
    for (int j=0; j<n2-1; ++j) {
      double dc   = (x2[j+1]-x2[j])*(x2upper[iz]-x2lower[iz]); //add EPS to avoid divide-by-zero???
      double dIdc = abs(((*I)(iz,j+1)-(*I)(iz,j))/dc);
      if (dIdc > dIdc_max) dIdc_max = dIdc;
    }
    x2tol = errmax*I_max/dIdc_max;
    if (mpi_rank==0)
      cout << " > CHEM_LOOKUP_TOL = " << errmax << " -> x2 tol = " << x2tol << endl;
    return;
  }

  void initInterpWts() {
    assert(invDenom1==NULL); invDenom1 = new double[n1-3][4];
    assert(invDenom2==NULL); invDenom2 = new double[n2-3][4];
    for (int i=0; i<n1-3; ++i)
      MiscUtils::cubicInterpDenom(invDenom1[i],&x1[i]);
    for (int i=0; i<n2-3; ++i)
      MiscUtils::cubicInterpDenom(invDenom2[i],&x2[i]);
    return;
  }

  // ================================================
  // 2D piecewise cubic interpolation
  // ================================================
  class ValueInterval {
  public:
    double Y1;
    double Y2;
    int i1;
    int i2;
    double w1[4];
    double w2[4];
  };

  class ValueInterval2 : public ValueInterval {
  public:
    int i2_2;
    double Y2_2 ; // second Y2 entry ...
    double w2_2[4];
  };

  void lookupSpecialNew(double* __restrict__ rout, const string& name,
                        const double* __restrict__ y1, const double *__restrict__ y2,
                        const int (*__restrict__ cvofa)[2], const int istart, const int iend) {

    // we are going to compute rout = f(y1,y2a) - f(y1,y2b)
    const ChemtableData2D* data = getVar(name);

    // i1, w1 find..
    for (int i = istart;i < iend; ++i) {

      const int icv0    = cvofa[i][0];
      const int icv1    = cvofa[i][1];
      const double y2a_ = y2[icv0];
      const double y2b_ = y2[icv1];

      if (abs(y2a_-y2b_) < x2tol) {
        // no meaningful change in y2... skip the lookup
        rout[i] = 0.0;
      }
      else {
        // change in y2 exceeds threshold... do integral lookup

        const double y1_  = 0.5*(y1[icv0] + y1[icv1]);
        double w1[4], w2[4], w2_2[4];

        const double Y1 = min(x1[n1-1],max(x1[0],y1_));
        const int i1 = min(n1-4,max(0,idxMap1->findIndex(Y1)-1));

        MiscUtils::cubicInterpWts(w1,invDenom1[i1],&x1[i1],Y1);

        double x2l = 0.0, x2u = 0.0;

        for (int k=0; k<4; ++k)
          x2l += w1[k]*x2lower[i1+k];

        for (int k=0; k<4; ++k)
          x2u += w1[k]*x2upper[i1+k];

        const double den  = x2u - x2l;
        const double Y2   = min(1.0,max(0.0,(y2a_-x2l)/(den+1.0e-16)));
        const double Y2_2 = min(1.0,max(0.0,(y2b_-x2l)/(den+1.0e-16)));

        const int i2   = min(n2-4,max(0,idxMap2->findIndex(Y2  )-1));
        const int i2_2 = min(n2-4,max(0,idxMap2->findIndex(Y2_2)-1));

        MiscUtils::cubicInterpWts(w2,invDenom2[i2],&x2[i2],Y2);
        MiscUtils::cubicInterpWts(w2_2,invDenom2[i2_2],&x2[i2_2],Y2_2);

        double r_1[4] = {0.0,0.0,0.0,0.0};
        double r_2[4] = {0.0,0.0,0.0,0.0};
        for (int k =0; k < 4 ; ++k) {
          const int ii       = (i1+k)*data->n2;
          const int offset_1 = ii+i2;   // should i precompute the offsets... ?
          const int offset_2 = ii+i2_2;

          for (int l=0; l < 4 ; ++l)
            r_1[l] += w1[k]* data->data[offset_1+l] ;

          for (int l=0; l < 4 ; ++l)
            r_2[l] += w1[k]* data->data[offset_2+l];
        }//k

        for (int k=0; k < 4 ; ++k)
          r_1[k] *= w2[k];

        for (int k =0; k < 4 ; ++k)
          r_2[k] *= w2_2[k];


        const double sum1 = (r_1[0] + r_1[1]) + (r_1[2] + r_1[3]);
        const double sum2 = (r_2[0] + r_2[1]) + (r_2[2] + r_2[3]);
        rout[i] = std::abs(sum1-sum2) ;
      }
    }//i
  }

  void lookupSpecialNew(double* __restrict__ rout1, const string& name1,
                        double* __restrict__ rout2, const string& name2,
                        const double* __restrict__ y1, const double *__restrict__ y2,
                        const int (*__restrict__ cvofa)[2], const int istart, const int iend) {

    // we are going to compute rout = f(y1,y2a) - f(y1,y2b)

    const ChemtableData2D* data1 = getVar(name1);
    const ChemtableData2D* data2 = getVar(name2);

    // ensure that the two pieces of data have the same size ..
    assert( data1->n2 == data2->n2);

    // i1, w1 find..
    for (int i = istart;i < iend; ++i) {

      const int icv0    = cvofa[i][0];
      const int icv1    = cvofa[i][1];
      const double y2a_ = y2[icv0];
      const double y2b_ = y2[icv1];

      if (abs(y2a_-y2b_) < x2tol) {
        // no meaningful change in y2... skip the lookup
        rout1[i] = 0.0;
        rout2[i] = 0.0;
      }
      else {
        // change in y2 exceeds threshold... do integral lookup

        const double y1_  = 0.5*(y1[icv0] + y1[icv1]);
        double w1[4], w2[4], w2_2[4];

        const double Y1 = min(x1[n1-1],max(x1[0],y1_));
        const int i1 = min(n1-4,max(0,idxMap1->findIndex(Y1)-1));

        MiscUtils::cubicInterpWts(w1,invDenom1[i1],&x1[i1],Y1);

        double x2l = 0.0, x2u = 0.0;

        for (int k=0; k<4; ++k)
          x2l += w1[k]*x2lower[i1+k];

        for (int k=0; k<4; ++k)
          x2u += w1[k]*x2upper[i1+k];

        const double den  = x2u - x2l;
        const double Y2   = min(1.0,max(0.0,(y2a_-x2l)/(den+1.0e-16)));
        const double Y2_2 = min(1.0,max(0.0,(y2b_-x2l)/(den+1.0e-16)));

        const int i2   = min(n2-4,max(0,idxMap2->findIndex(Y2  )-1));
        const int i2_2 = min(n2-4,max(0,idxMap2->findIndex(Y2_2)-1));

        MiscUtils::cubicInterpWts(w2,invDenom2[i2],&x2[i2],Y2);
        MiscUtils::cubicInterpWts(w2_2,invDenom2[i2_2],&x2[i2_2],Y2_2);

        double r_11[4] = {0.0,0.0,0.0,0.0};
        double r_21[4] = {0.0,0.0,0.0,0.0};

        double r_12[4] = {0.0,0.0,0.0,0.0};
        double r_22[4] = {0.0,0.0,0.0,0.0};

        for (int k =0; k < 4 ; ++k) {

          const int ii       = (i1+k)*data1->n2;
          const int offset_1 = ii+i2;   // should i precompute the offsets... ?
          const int offset_2 = ii+i2_2;

          for (int l=0; l < 4 ; ++l) {
            r_11[l] += w1[k]* data1->data[offset_1+l] ;
            r_12[l] += w1[k]* data2->data[offset_1+l] ;
          }

          for (int l=0; l < 4 ; ++l) {
            r_21[l] += w1[k]* data1->data[offset_2+l];
            r_22[l] += w1[k]* data2->data[offset_2+l];
          }
        }//k

        for (int k=0; k < 4 ; ++k) {
          r_11[k] *= w2[k];
          r_12[k] *= w2[k];
        }

        for (int k =0; k < 4 ; ++k) {
          r_21[k] *= w2_2[k];
          r_22[k] *= w2_2[k];
        }


        const double sum11 = (r_11[0] + r_11[1]) + (r_11[2] + r_11[3]);
        const double sum21 = (r_21[0] + r_21[1]) + (r_21[2] + r_21[3]);
        rout1[i] = std::abs(sum11-sum21) ;

        const double sum12 = (r_12[0] + r_12[1]) + (r_12[2] + r_12[3]);
        const double sum22 = (r_22[0] + r_22[1]) + (r_22[2] + r_22[3]);
        rout2[i] = std::abs(sum12-sum22) ;
      }
    }//i
  }

  void lookupSpecial(double* rout, const string& name,
                     const double *__restrict__ y1, const double *__restrict__ y2a,
                     const double *__restrict__ y2b, const int n) {

    // we are going to compute rout = f(y1,y2a) - f(y1,y2b)
    const ChemtableData2D* data = getVar(name);

    // i1, w1 find..
    for (int i =0;i <n; ++i) {

      double w1[4], w2[4], w2_2[4];

      const double Y1 = min(x1[n1-1],max(x1[0],y1[i]));
      const int i1 = min(n1-4,max(0,idxMap1->findIndex(Y1)-1));

      MiscUtils::cubicInterpWts(w1,invDenom1[i1],&x1[i1],Y1);

      double x2l = 0.0, x2u = 0.0;

      for (int k=0; k<4; ++k)
        x2l += w1[k]*x2lower[i1+k];

      for (int k=0; k<4; ++k)
        x2u += w1[k]*x2upper[i1+k];

      const double den  = x2u - x2l;
      const double Y2   = min(1.0,max(0.0,(y2a[i]-x2l)/(den+1.0e-16)));
      const double Y2_2 = min(1.0,max(0.0,(y2b[i]-x2l)/(den+1.0e-16)));

      const int i2   = min(n2-4,max(0,idxMap2->findIndex(Y2  )-1));
      const int i2_2 = min(n2-4,max(0,idxMap2->findIndex(Y2_2)-1));

      MiscUtils::cubicInterpWts(w2,invDenom2[i2],&x2[i2],Y2);
      MiscUtils::cubicInterpWts(w2_2,invDenom2[i2_2],&x2[i2_2],Y2_2);

      double r_1[4] = {0.0,0.0,0.0,0.0};
      double r_2[4] = {0.0,0.0,0.0,0.0};
      for (int k =0; k < 4 ; ++k) {
        const int ii       = (i1+k)*data->n2;
        const int offset_1 = ii+i2;   // should i precompute the offsets... ?
        const int offset_2 = ii+i2_2;

        for (int l=0; l < 4 ; ++l)
          r_1[l] += w1[k]* data->data[offset_1+l] ;

        for (int l=0; l < 4 ; ++l)
          r_2[l] += w1[k]* data->data[offset_2+l];
      }//k

      for (int k=0; k < 4 ; ++k)
        r_1[k] *= w2[k];

      for (int k =0; k < 4 ; ++k)
        r_2[k] *= w2_2[k];


      const double sum1 = (r_1[0] + r_1[1]) + (r_1[2] + r_1[3]);
      const double sum2 = (r_2[0] + r_2[1]) + (r_2[2] + r_2[3]);
      rout[i] = std::abs(sum1-sum2) ;
    }//i
  }//lookupSpecial...

  void lookupCubicPtrVector(vector<double*>& routPtrVec, const vector<string>& nameVec,
                   const double * y1, const double * y2,const int n) {

    vector<ChemtableData2D*> dataPtrVec(nameVec.size());
    for (unsigned int iname = 0; iname < nameVec.size(); ++iname)
      dataPtrVec[iname] = getVar(nameVec[iname]);


    // i1, w1 find..
    for (int i =0;i <n; ++i) {

      double w1[4], w2[4];

      const double Y1  = min(x1[n1-1],max(x1[0],y1[i]));
      const int i1     = min(n1-4,max(0,idxMap1->findIndex(Y1)-1));
      MiscUtils::cubicInterpWts(w1,invDenom1[i1],&x1[i1],Y1);

      double x2l = 0.0, x2u = 0.0;
      for (int k=0; k<4; ++k)
        x2l += w1[k]*x2lower[i1+k];

      for (int k=0; k<4; ++k)
        x2u += w1[k]*x2upper[i1+k];

      const double Y2  = min(1.0,max(0.0,(y2[i]-x2l)/(x2u-x2l+1.0e-16)));
      const int i2 = min(n2-4,max(0,idxMap2->findIndex(Y2)-1));
      MiscUtils::cubicInterpWts(w2,invDenom2[i2],&x2[i2],Y2);

      const int vec_size = int(routPtrVec.size());
      for (int j=0; j < vec_size; ++j) {

        const double* data = dataPtrVec[j]->data;
        const int n2_      = dataPtrVec[j]->n2;
        double r[4] = {0.0, 0.0, 0.0, 0.0};

        for (int k =0; k<4; ++k) {
          double tmp[4];
          const double w1k = w1[k];
          // this vectorized the two loops ===================
          // and was the fastest so far...
          const int offset = (i1+k)*n2_+i2;
          for (int l=0; l <4;++l)
            tmp[l] = w1k*data[offset+l];

          for (int l = 0; l <4; ++l)
            r[l] += tmp[l];
          //================================================
        }

        for (int k=0; k <4 ;++k)
          r[k] *= w2[k];

        routPtrVec[j][i] = (r[0]+r[1])+(r[2]+r[3]);
      }
    }
  }


  // ===========================
  // look up vector of variables
  // ===========================
  void lookupLinearPtrVector(vector<double*>& routPtrVec, const vector<string>& nameVec,
		       const double* x1, const double* x2,const int n) {

    vector<ChemtableData2D*> dataPtrVec(nameVec.size());
    for (unsigned int iname = 0; iname < nameVec.size(); ++iname)
      dataPtrVec[iname] = getVar(nameVec[iname]);

    // compute the indecies and integration weights
    // first variable
    int * ti1   = new int[n]; // table index 1 (lower bound)
    double* W1 = new double[n];
    // second variable
    int * ti2   = new int[n];
    double* W2 = new double[n];
    findIndexAndComputeLinearWeights(ti1,W1,ti2,W2,x1,x2,n);

    for ( int i=0; i<n; ++i) {
      // now having the local index, get the table index
      int t1 = ti1[i];
      int t2 = ti2[i];
      // get the weights
      double w1 = W1[i];
      double w2 = W2[i];
      // linear interpolation weights
      double c1 = w1*w2;
      double c3 = w1*(1.0-w2);
      double c5 = (1.0-w1)*w2;
      double c7 = (1.0-w1)*(1.0-w2);

      // get the variable and interpolate
      for ( unsigned int j=0; j<routPtrVec.size(); ++j) {
	routPtrVec[j][i] =
	  c1 * (*(dataPtrVec[j]))(t1,t2)     +
	  c3 * (*(dataPtrVec[j]))(t1,t2+1)   +
	  c5 * (*(dataPtrVec[j]))(t1+1,t2)   +
	  c7 * (*(dataPtrVec[j]))(t1+1,t2+1);
      }
    }

    // delete the allocated memory
    delete[] ti1;
    delete[] W1;

    delete[] ti2;
    delete[] W2;

  }




  // ===================================
  // for each variable "var" compute the
  // lower bound index in the table and
  // the weight associated with the
  // integration
  // ===================================

  void findIndexAndComputeLinearWeights(int * ti1, double* w1,
					int * ti2, double* w2,
					const double* var1,
					const double* var2,
					const int nv ) {
    findIndexAndComputeLinearWeights(ti1,w1,var1,nv,x1,n1);
    findIndexAndComputeLinearWeights(ti2,w2,var2,nv,x2,n2);
  }

  void findIndexAndComputeLinearWeights(int * ti1, double* w1,
					const double* var1,
					const int nv ) {
    findIndexAndComputeLinearWeights(ti1,w1,var1,nv,x1,n1);
  }

  void findIndexAndComputeLinearWeights(int * ti, double* w,
					const double* var,
					const int nv,
					const double* tvar,
					const int nt) {
    assert(nt>1);

    // pair the variables with the index
    vector<myPair> varPair;
    for ( int i=0; i<nv; ++i )
      varPair.push_back(myPair(i,var[i]));

    // sort the vector of pairs
    lessSecond<int,double> comparator;
    std::sort(varPair.begin(),varPair.end(),comparator);

    // find the index
    int tindex = 0;
    while  ( tvar[tindex+1] == tvar[tindex] )
      ++tindex;
    double one_over_length = 1.0 / ((tvar[tindex+1] - tvar[tindex])+1.0E-20);
    for (vector<myPair>::iterator thisPair = varPair.begin(); thisPair!=varPair.end(); ++thisPair) {
      int i = thisPair->first;
      double thisvar = thisPair->second;
      if (thisvar >= tvar[nt-1]) {
	ti[i] =  nt-2;
	one_over_length = 1.0 / ((tvar[nt-1] - tvar[nt-2])+1.0E-20);
      }
      else {
	while (thisvar > tvar[tindex+1]) {
	  ++tindex;
	  one_over_length = 1.0 / ((tvar[tindex+1] - tvar[tindex])+1.0E-20);
	}
	ti[i] = tindex;
      }
      // compute integration weight
      // IMPORTANT: clip the out of bound variavles to the table limites
      // this is ad hoc but it seems we need to clip the
      // variables close to the table limit to the bondaries
      if ( thisvar < (tvar[0]+SMALL_CLIP) )     thisvar = tvar[0];
      if ( thisvar > (tvar[nt-1]-SMALL_CLIP) )  thisvar = tvar[nt-1];
      double tmp = thisvar - tvar[ti[i]];
      w[i] = 1.0 - tmp*one_over_length;
    }
    assert(tindex<(nt-1));

  }

  void loadVariables(const vector<string>& strVec) {
    if (strVec.size() == 0)
      CERR("pointless loadVariables()");

    for (unsigned int istr = 0; istr < strVec.size(); ++istr)
      getVar(strVec[istr]);
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


  // ============================
  // delete a variable from table
  // ============================
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



  // ===========================
  // read header and coordinates
  // ===========================
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
      // read the header
      if (mpi_rank == 0 ) {
        my_fseek(fp,offset);
        fread(&header,sizeof(Header),1,fp);
       }

      MPI_Bcast(&header,1,MPI_Header,0,mpi_comm);

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
	if ( mpi_rank == 0 ) {
	  size_t dummy= fread(&(tmpPtr[0]),sizeof(double),header.idata[0],fp);
	  assert(int(dummy)==header.idata[0]);
        }

        MPI_Bcast(&(tmpPtr[0]),header.idata[0],MPI_DOUBLE,0,mpi_comm);
	if (byte_swap) ByteSwap::byteSwap(&(tmpPtr[0]),header.idata[0]);
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

    // build a sorted set of the var names for fast checking to see
    // if a given var exists ....

    assert( sorted_vars.size() == 0);
    for (int ii = 0; ii < nvars; ++ii)
      sorted_vars.insert(varName[ii]);

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
          assert(int(dummy)==header.idata[0]);
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



public:


  // =================
  // test table lookup
  // =================
  void testLookup() {

    // Do it on rank zero only


    // set the clipping to zero
    double SMALL_CLIP_orig = SMALL_CLIP;
    SMALL_CLIP = 0.0;

    if ( mpi_rank == 0 )
      cout << "testing table lookup ..." << endl;
    // pick two variables in the table
    ChemtableData2D * var1 = getVar("T");
    if ( mpi_rank == 0 )
      cout << " > replacing T with linear function 2*Z + 3*C"<< endl;
    ChemtableData2D * var2 = getVar("rho");
    if ( mpi_rank == 0 )
      cout << " > replacing rho with linear function 20.1*Z + 29.0001*C"<< endl;
    // set them with linear functions
    for ( int i=0; i<n1; ++i)
      for ( int j=0; j<n2; ++j){
	(*var1)(i,j) = 2.0*x1[i] + 3.0*x2[j];
	(*var2)(i,j) = 20.1*x1[i]+29.0001*x2[j];
      }

    // Build a random array
    int myN = 100000;
    double *myZ = new double[myN];
    double *myC = new double[myN];
    //srand(10);
    for ( int i=0; i<myN; ++i ) {
      myZ[i] = x1[0] + (x1[n1-1]-x1[0])*(double)(rand())/(double)(RAND_MAX);
      myC[i] = x2[0] + (x2[n2-1]-x2[0])*(double)(rand())/(double)(RAND_MAX);
    }

    // result vector
    double *rout1 = new double[myN];
    double *rout2 = new double[myN];

    // do lookup
    lookup(rout1, "T", rout2, "rho", myZ,myC,myN);

    // compute error
    double myerror1 = 0.0;
    for ( int i=0; i<myN; ++i)
      myerror1 = max(fabs(rout1[i]-(2.0*myZ[i] + 3.0*myC[i])),myerror1);
    double myerror2 = 0.0;
    for ( int i=0; i<myN; ++i)
      myerror2 = max(fabs(rout2[i]-(20.1*myZ[i]+29.0001*myC[i])),myerror2);
    // take the max across processors
    double error1;
    MPI_Reduce(&myerror1,&error1,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    double error2;
    MPI_Reduce(&myerror2,&error2,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    if ( mpi_rank == 0 )
      cout << "  >> error in lookups of the two linear functions (should be zero) = " << error1 << " , " << error2 << endl;

    // set the clipping back to original one
    SMALL_CLIP = SMALL_CLIP_orig;

    // deletle variables from table
    deleteVar("T");
    deleteVar("rho");

    // delete memory
    delete[] myZ;
    delete[] myC;
    delete[] rout1;
    delete[] rout2;

    if ( mpi_rank == 0 )
      cout << "done with table lookup test" << endl;

  }

  bool varExists(const string& var_name) const {

    return (sorted_vars.find(var_name) != sorted_vars.end());

  }

};

// =========================
// Cartesian chemtable class
// =========================
class CartesianChemtable3D: public AbstractChemtable3D {

private:
  int n1, n2, n3, nvars;
  double*x1, *x2, *x3;

  string *varName;
  set<string> sorted_vars;

  list<ChemtableData3D> datalist;
  typedef std::pair<int,double> myPair;
  typedef std::pair<int,int> myPairInt;
  double SMALL_CLIP;

public:

  // constructor
  CartesianChemtable3D(const string& name, const string& type){
    if ( mpi_rank == 0 )
      cout << "CartesianChemtable3D():" << endl;

    // table dimensions
    n1    = 0;
    n2    = 0;
    n3    = 0;
    nvars = 0;

    // coordinates
    x1 = NULL;
    x2 = NULL;
    x3 = NULL;

    // array of variable names
    varName = NULL;

    // small value used for clipping
    SMALL_CLIP = 1.01e-6;

    // initialize table
    this->tableName = name;
    this->tableType = type;
    init();
  }

  // deconstructor
  ~CartesianChemtable3D(){
    if ( mpi_rank == 0 )
      cout << "~CartesianChemtable3D():" << endl;

    DELETE(x1);
    DELETE(x2);
    DELETE(x3);
    DELETE(varName);
  }

  void info(){

    if ( mpi_rank == 0 ) {
      cout << "============== chemistry table information ===============" << endl;
      cout << " chemistry table information: " << endl;
      cout << "  > name:                     " << tableName << endl;
      cout << "  > type:                     " << tableType << endl;
      cout << "  > reference pressure:       " << pressure << endl;
      cout << "  > dimension:                ("<< n1 <<"," << n2 <<"," << n3 << ")" << endl;
      cout << "   > first coordinate range:   ("<< x1[0] <<"," << x1[n1-1] << ")" << endl;
      cout << "   > second coordinate range:  ("<< x2[0] <<"," << x2[n2-1] << ")" << endl;
      cout << "   > third coordinate range:   ("<< x3[0] <<"," << x3[n3-1] << ")" << endl;
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

      for (list<ChemtableData3D>::iterator data = datalist.begin(); data != datalist.end(); ++data) {
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
    for (list<ChemtableData3D>::iterator data = datalist.begin(); data != datalist.end(); ++data) {
      if ( data->name == name ) {
	varmin = (*data)(0,0,0);
	varmax = (*data)(0,0,0);
	for ( int i=0; i<n1; ++i)
	  for ( int j=0; j<n2; ++j)
	    for ( int k=0; k<n3; ++k){
	      varmax = max(varmax,(*data)(i,j,k));
	      varmin = min(varmin,(*data)(i,j,k));
	    }
	return(1);
      }
    }
    if ( mpi_rank == 0)
      cout << "***WARNING in getTableVarRange: the variable "
	   << name << " is not loaded yet/n" << " ...setting the range to 0:0 " << endl;

    return(0);

  }

  // ==========================================
  // write the table with the current variables
  // the one that are loaded to a tecplot file
  // ==========================================
  void writeTecplot(const string& filename) {

    // load all the variables
    for (int i=0; i<nvars; ++i) getVar(varName[i]);

    // only rank 0 needs to write
    if (mpi_rank == 0) {

      cout << " > writing 3d cartesian chemtable " << tableName << " to the Tecplot file "<< filename << endl;
      cout << "  > the following variables are currently loaded: " << endl;
      for (list<ChemtableData3D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
	cout << "   > " << data->name << endl ;

      FILE * fp = fopen(filename.c_str(),"w");
      fprintf(fp,"title = \"chemtable\"\n");
      fprintf(fp,"variables = \"Z\"\n\"Zvar\"\n\"C\"\n");
      for (list<ChemtableData3D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
	fprintf(fp,"\"%s\"\n",data->name.c_str());
      fprintf(fp,"zone i = %d, j = %d, k = %d, DATAPACKING=POINT\n",n1,n2,n3);
      for (int k=0; k <n3; ++k)
	for (int j=0; j <n2; ++j)
	  for (int i=0; i <n1; ++i){
	    fprintf(fp,"%lf %lf %lf",x1[i],x2[j],x3[k]);
	    for (list<ChemtableData3D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
	      fprintf(fp," %lf",(*data)(i,j,k));
	    fprintf(fp,"\n");
	  }
      fclose(fp);
    }
  }

  void writeTecplot2(const string& filename) {

    // write zvar planes

    // load all the variables
    for (int i=0; i<nvars; ++i) getVar(varName[i]);

    // only rank 0 needs to write
    if (mpi_rank == 0) {

      for (int j=0; j <n2; ++j) {
        cout << " > writing Zvar plane " << j << endl;
        char index[21];
        sprintf(index,".%04d.dat",j);
        string thisfile = filename + index;
        FILE * fp = fopen(thisfile.c_str(),"w");
        fprintf(fp,"title = \"chemtable\"\n");
        fprintf(fp,"variables = \"Z\"\n\"Zvar\"\n\"C\"\n");
        for (list<ChemtableData3D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
 	  fprintf(fp,"\"%s\"\n",data->name.c_str());
        fprintf(fp,"zone i = %d, j = %d, k = %d, DATAPACKING=POINT\n",n1,1,n3);

        for (int k=0; k <n3; ++k)
	  for (int i=0; i <n1; ++i){
	    fprintf(fp,"%lf %lf %lf",x1[i],x2[j],x3[k]);
	    for (list<ChemtableData3D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
	      fprintf(fp," %lf",(*data)(i,j,k));
	    fprintf(fp,"\n");
	  }
        fclose(fp);
      }
    }
  }

private:

  void lookupDataVec(vector<double*>& routPtrVec, const vector<string>& nameVec,
		     const double* x1, const double* x2, const double* x3,const int n) {
    lookupLinearPtrVector(routPtrVec,nameVec,x1,x2,x3,n);
  }

  void init() {
    readHeaderAndCoordinates();
    info();
  }

  // ===========================
  // linear interpolation of the table data..
  // ===========================
  void lookupLinearPtrVector(vector<double*>& routPtrVec, const vector<string>& nameVec,
		       const double* x1, const double* x2, const double* x3,const int n) {

    int n1 = -1, n2 =-1, n3 = -1;
    vector<ChemtableData3D *> dataPtrVec(nameVec.size());
    for (unsigned int iname = 0; iname < nameVec.size(); ++iname) {
      dataPtrVec[iname] = getVar(nameVec[iname]);
      if ( n1 == -1) {
        n1 = dataPtrVec[iname]->n1;
        n2 = dataPtrVec[iname]->n2;
        n3 = dataPtrVec[iname]->n3;
      }
    }

    // compute the indecies and integration weights
    int * ti1   = new int[n]; // table index 1 (lower bound)
    double* W1 = new double[n];
    int * ti2   = new int[n];
    double* W2 = new double[n];
    int * ti3   = new int[n];
    double* W3 = new double[n];

    findIndexAndComputeLinearWeights(ti1,W1,ti2,W2,ti3,W3,x1,x2,x3,n);

    for ( int i=0; i<n; ++i) {
      // now having the local index, get the table index
      const int t1 = ti1[i];
      const int t2 = ti2[i];
      const int t3 = ti3[i];

      // get the weights
      const double w1 = W1[i];
      const double w2 = W2[i];
      const double w3 = W3[i];

      // linear interpolation weights
      const double c1 = w1*w2*w3;
      const double c2 = w1*w2*(1.0-w3);
      const double c3 = w1*(1.0-w2)*w3;
      const double c4 = w1*(1.0-w2)*(1.0-w3);
      const double c5 = (1.0-w1)*w2*w3;
      const double c6 = (1.0-w1)*w2*(1.0-w3);
      const double c7 = (1.0-w1)*(1.0-w2)*w3;
      const double c8 = (1.0-w1)*(1.0-w2)*(1.0-w3);


      // get the variable and interpolate
      for ( unsigned int j=0; j<routPtrVec.size(); ++j) {
        /*
	routPtrVec[j][i] =
	  c1 * (*(dataPtrVec[j]))(t1,t2,t3)     + c2 * (*(dataPtrVec[j]))(t1,t2,t3+1)      +
	  c3 * (*(dataPtrVec[j]))(t1,t2+1,t3)   + c4 * (*(dataPtrVec[j]))(t1,t2+1,t3+1)    +
	  c5 * (*(dataPtrVec[j]))(t1+1,t2,t3)   + c6 * (*(dataPtrVec[j]))(t1+1,t2,t3+1)    +
	  c7 * (*(dataPtrVec[j]))(t1+1,t2+1,t3) + c8 * (*(dataPtrVec[j]))(t1+1,t2+1,t3+1);
         */

        const double * data     = dataPtrVec[j]->data;

        const int index0  = (t1*n2 + t2)*n3 + t3;
        routPtrVec[j][i]  = c1*data[index0] + c2*data[index0+1];

        const int index1  = (t1*n2 + t2+1)*n3 + t3;
        routPtrVec[j][i] += c3*data[index1] + c4*data[index1+1];

        const int index2  = ((t1+1)*n2 + t2)*n3 + t3;
        routPtrVec[j][i] += c5*data[index2] + c6*data[index2+1];

        const int index3  = ((t1+1)*n2 + t2+1)*n3 + t3;
        routPtrVec[j][i] += c7*data[index3] + c8*data[index3+1];
      }
    }

    delete[] ti1;
    delete[] W1;
    delete[] ti2;
    delete[] W2;
    delete[] ti3;
    delete[] W3;
  }


  // ===================================
  // for each variable "var" compute the
  // lower bound index in the table and
  // the weight associated with the
  // integration
  // ===================================
  void findIndexAndComputeLinearWeights(int * ti1, double* w1, int * ti2, double* w2,
					int * ti3, double* w3, const double* var1,
					const double* var2, const double* var3, const int nv) {
    findIndexAndComputeLinearWeightsUnsorted(ti1,w1,var1,nv,x1,n1);
    findIndexAndComputeLinearWeightsUnsorted(ti2,w2,var2,nv,x2,n2);
    findIndexAndComputeLinearWeightsUnsorted(ti3,w3,var3,nv,x3,n3);
  }

  void findIndexAndComputeLinearWeights(int * ti1, double* w1, int * ti2, double* w2,
					const double* var1, const double* var2, const int nv ) {
    findIndexAndComputeLinearWeightsUnsorted(ti1,w1,var1,nv,x1,n1);
    findIndexAndComputeLinearWeightsUnsorted(ti2,w2,var2,nv,x2,n2);
  }

  void findIndexAndComputeLinearWeights(int * ti1, double* w1, const double* var1, const int nv ) {
    findIndexAndComputeLinearWeightsUnsorted(ti1,w1,var1,nv,x1,n1);
  }

  void findIndexAndComputeLinearWeightsUnsorted(int* ti, double* w, const double* var,
                                                const int nv, const double* tvar, const int nt)  {
    assert(nt > 1);

    for (int i = 0; i < nv; ++i) {
      double thisvar  = max(var[i],tvar[0]);
      if ( thisvar >= tvar[nt-1]) {
        thisvar  = tvar[nt-1];
        ti[i] = nt-2;
      } else {
        MiscUtils::findInterval(ti[i],thisvar,tvar,nt);   // replace with the idxMap
      }

      assert(ti[i] < nt-1);
      double one_over_length = 1.0/(tvar[ti[i]+1] - tvar[ti[i]]);
      w[i] = 1.0 - (thisvar - tvar[ti[i]])*one_over_length;
    }
  }

  void findIndexAndComputeLinearWeights(int * ti, double* w, const double* var, const int nv,
					const double* tvar, const int nt) {

    assert(nt>1);

    // pair the variables with the index
    vector<myPair> varPair;
    for ( int i=0; i<nv; ++i )
      varPair.push_back(myPair(i,var[i]));

    // sort the vector of pairs
    lessSecond<int,double> comparator;
    std::sort(varPair.begin(),varPair.end(),comparator);

    // find the index
    int tindex = 0;
    double one_over_length = 1.0 / (tvar[tindex+1] - tvar[tindex]);
    for (vector<myPair>::iterator thisPair = varPair.begin(); thisPair!=varPair.end(); ++thisPair) {
      int i = thisPair->first;
      double thisvar = thisPair->second;
      if (thisvar >= tvar[nt-1]) {
	ti[i] =  nt-2;
	one_over_length = 1.0 / (tvar[nt-1] - tvar[nt-2]);
      }
      else {
	while (thisvar > tvar[tindex+1]) {
	  ++tindex;
	  one_over_length = 1.0 / (tvar[tindex+1] - tvar[tindex]);
	}
	ti[i] = tindex;
      }
      // compute integration weight
      // IMPORTANT: clip the out of bound variavles to the table limits
      // this is ad hoc but it seems we need to clip the
      // variables close to the table limit to the bondaries
      if ( thisvar < (tvar[0]+SMALL_CLIP) )     thisvar = tvar[0];
      if ( thisvar > (tvar[nt-1]-SMALL_CLIP) )  thisvar = tvar[nt-1];
      double tmp = thisvar - tvar[ti[i]];
      w[i] = 1.0 - tmp*one_over_length;
    }
    assert(tindex<(nt-1));

  }

  void loadVariables(const vector<string>& strVec) {
    if (strVec.size() == 0)
      CERR("pointless loadVariables()");

    for (unsigned int istr = 0; istr < strVec.size(); ++istr)
      getVar(strVec[istr]);
  }

  // =======================
  // add a variable to table
  // =======================
  ChemtableData3D * getVar(const string& varname){

    // double check that the variable is not already in the list
    for (list<ChemtableData3D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
      if ( data->name == varname ) return(&(*data));

    // add the variable to the list
    // this needs to be done on all the processors

    datalist.push_back(ChemtableData3D(varname,n1,n2,n3));
    ChemtableData3D * thisdata = &(datalist.back());
    assert(thisdata->name == varname);

    readData(thisdata);

    // if the variable is the source term, clean up
    if ( thisdata->name == "src_prog" )
      cleanupSrcProg(thisdata);

    return thisdata;
  }


  // ============================
  // delete a variable from table
  // ============================
  void deleteVar(const string& varname){

    // double check that the variable is not already in the list
    for (list<ChemtableData3D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
      if ( data->name == varname ) {
	if ( mpi_rank == 0)
	  cout << " > deleting variable "<< data->name << " from chemtable" << endl;
	datalist.erase(data);
	return;
      }

  }


  // ================================================================
  // cleanup source term such that SRC_PROG is zero at C=0 and C=Cmax
  // ================================================================
  void cleanupSrcProg(ChemtableData3D * data){

    // check the name
    assert(data->name == "src_prog");

    if (mpi_rank == 0 )
      cout << " > cleaning up source term to assure it is zero at C=0 and C=Cmax ..."<<endl;

    // read progress variable
    int table_has_PROG = 0;
    // double check that the variable is not already in the list
    for (list<ChemtableData3D>::iterator thisdata = datalist.begin(); thisdata != datalist.end(); ++thisdata)
      if ( thisdata->name == "prog" ) table_has_PROG = 1;
    ChemtableData3D * Progdata = getVar("prog");

    // set the source term to zero if c<SMALL_CLIP or c>Cmax
    double max_SrcProg_zero = 0.0;
    for ( int i=0; i<n1; ++i)
      for ( int j=0; j<n2; ++j) {
	double Cmax = (*Progdata)(i,j,n3-1);
	for ( int k=0; k<n3; ++k)
	  if ( (x3[k]<SMALL_CLIP) || (x3[k]>=Cmax) || ((*data)(i,j,k)<0.0) ) {
	    max_SrcProg_zero = max(max_SrcProg_zero,fabs((*data)(i,j,k)));
	    (*data)(i,j,k) = 0.0;
	  }
      }

    if (mpi_rank == 0)
      cout << "   > maximum value set to zero = " << max_SrcProg_zero <<endl;

    // delete data
    if ( !table_has_PROG  ) deleteVar("prog");

  }



  // ===========================
  // read header and coordinates
  // ===========================
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
      if ( header.id == UGP_IO_CT_CART_3D ) {
	// type of the table
	tableType = header.name;
	// table dimensions
	n1    = header.idata[0]; assert(n1>1);
	n2    = header.idata[1]; assert(n2>1);
	n3    = header.idata[2]; assert(n3>1);
	nvars = header.idata[3];
	// allocate space for table coordinates
	assert(x1==NULL);      x1 = new double[n1];
	assert(x2==NULL);      x2 = new double[n2];
	assert(x3==NULL);      x3 = new double[n3];
	assert(varName==NULL); varName = new string[nvars];
	// reference quantities
	pressure   = header.rdata[0]; assert(pressure>0.0);
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
	else if ( name == "Zvar" ) {
	  assert( header.idata[0] == n2 );
	  assert(x2!=NULL);
	  tmpPtr = x2;
	  ++icoor;
	}
	else if ( name == "C" ) {
	  assert( header.idata[0] == n3 );
	  assert(x3!=NULL);
	  tmpPtr = x3;
	  ++icoor;
	}
	else{
	  cerr << ERRSTART << "Error: unknown coordinate " << header.name << ERREND << endl;
	  throw(0);
	}
	// read the data
	size_t dummy= fread(&(tmpPtr[0]),sizeof(double),header.idata[0],fp);
	assert(int(dummy)==header.idata[0]);
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
	cerr << ERRSTART << "Error: could not find one of the three coordinates or the header " << ERREND << endl;
	throw(0);
      }

      // increase the offset
      offset += header.skip;

    } // while (done != 1)

    // close the file
    fclose(fp);

    assert(ivar==nvars);
    assert(icoor==3);

    // build a sorted set of the var names for fast checking to see
    // if a given var exists ....

    assert( sorted_vars.size() == 0);
    for (int ii = 0; ii < nvars; ++ii)
      sorted_vars.insert(varName[ii]);


  }

  // ======================================
  // read the data associated with a header
  // ======================================
  void readData(ChemtableData3D * data) {

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
	  assert( header.idata[0] == (n1*n2*n3) );
	  assert( header.idata[1] == n1 );
	  assert( header.idata[2] == n2 );
	  assert( header.idata[3] == n3 );
	  // read the data
	  size_t dummy= fread(&((*data)(0,0,0)),sizeof(double),header.idata[0],fp);
	  assert(int(dummy)==header.idata[0]);
	  // byte swap
	  if (byte_swap)
	    ByteSwap::byteSwap(&((*data)(0,0,0)),header.idata[0]);
	  // dump variable range
	  double varmin = (*data)(0,0,0);
	  double varmax = (*data)(0,0,0);
	  for ( int i=0; i<n1; ++i)
	    for ( int j=0; j<n2; ++j)
	      for ( int k=0; k<n3; ++k){
		varmax = max(varmax,(*data)(i,j,k));
		varmin = min(varmin,(*data)(i,j,k));
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

    MPI_Bcast(&((*data)(0,0,0)),data->getSize(),MPI_DOUBLE,0,mpi_comm);
  }

public:

  // =================
  // test table lookup
  // =================
  void testLookup() {

    // Do it on rank zero only


    // set the clipping to zero
    double SMALL_CLIP_orig = SMALL_CLIP;
    SMALL_CLIP = 0.0;

    if ( mpi_rank == 0 )
      cout << "testing table lookup ..." << endl;
    // pick two variables in the table
    ChemtableData3D * var1 = getVar("T");
    if ( mpi_rank == 0 )
      cout << " > replacing T with linear function 2*Z + 3*Zv + 4*C"<< endl;
    ChemtableData3D * var2 = getVar("rho");
    if ( mpi_rank == 0 )
      cout << " > replacing rho with linear function 20.1*Z + 31.0*Zv + 29.0001*C"<< endl;
    // set them with linear functions
    for ( int i=0; i<n1; ++i)
      for ( int j=0; j<n2; ++j)
	for ( int k=0; k<n3; ++k) {
	  (*var1)(i,j,k) = 2.0*x1[i] + 3.0*x2[j] + 4.0*x3[k];
	  (*var2)(i,j,k) = 20.1*x1[i] - 31.0*x2[j] + 29.0001*x3[k];
	}

    // Build a random array
    int myN = 100000;
    double *myZ  = new double[myN];
    double *myZv = new double[myN];
    double *myC  = new double[myN];
    //srand(10);
    for ( int i=0; i<myN; ++i ) {
      myZ[i]  = x1[0] + (x1[n1-1]-x1[0])*(double)(rand())/(double)(RAND_MAX);
      myZv[i] = x2[0] + (x2[n2-1]-x2[0])*(double)(rand())/(double)(RAND_MAX);
      myC[i]  = x3[0] + (x3[n3-1]-x3[0])*(double)(rand())/(double)(RAND_MAX);
    }

    // result vector
    double *rout1 = new double[myN];
    double *rout2 = new double[myN];

    // do lookup
    lookup(rout1, "T", rout2, "rho", myZ,myZv,myC,myN);

    // compute error
    double myerror1 = 0.0;
    for ( int i=0; i<myN; ++i)
      myerror1 = max(fabs(rout1[i]-(2.0*myZ[i] + 3.0*myZv[i] + 4.0*myC[i])),myerror1);
    double myerror2 = 0.0;
    for ( int i=0; i<myN; ++i)
      myerror2 = max(fabs(rout2[i]-(20.1*myZ[i] - 31.0*myZv[i] + 29.0001*myC[i])),myerror2);
    // take the max across processors
    double error1;
    MPI_Reduce(&myerror1,&error1,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    double error2;
    MPI_Reduce(&myerror2,&error2,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    if ( mpi_rank == 0 )
      cout << "  >> error in lookups of the two linear functions (should be zero) = " << error1 << " , " << error2 << endl;

    // set the clipping back to original one
    SMALL_CLIP = SMALL_CLIP_orig;

    // deletle variables from table
    deleteVar("T");
    deleteVar("rho");

    // delete memory
    delete[] myZ;
    delete[] myZv;
    delete[] myC;
    delete[] rout1;
    delete[] rout2;

    if ( mpi_rank == 0 )
      cout << "done with table lookup test" << endl;

  }


  void lookupDataVecReduced(vector<double*>& routPtrVec, const vector<string>& nameVec, 
                            const double * y1, const int n) {


    // this takes that actual value of Z, a fake large value of C, and no turbulence variance.

    double * C_tmp    = new double[n]; 
    double * Zvar_tmp = new double[n];

    for (int i = 0; i < n; ++i) { 

      C_tmp[i]    = 100.0; // ... a "large" C-value, assuming C is weighted sum of mass fractions
      Zvar_tmp[i] = 0.0;

    }

    lookupDataVec(routPtrVec,nameVec,y1,Zvar_tmp,C_tmp,n); 

    delete[] C_tmp;
    delete[] Zvar_tmp;

  }

  bool varExists(const string& var_name) const {
    return (sorted_vars.find(var_name) != sorted_vars.end());
  }

};

// =========================
// Cartesian chemtable class
// =========================
class CartesianChemtable4D: public AbstractChemtable4D {

private:

  string tableName;
  string tableType;
  double  pressure;

  int n1, n2, n3, n4, nvars;
  double*x1, *x2, *x3, *x4;

  string *varName;

  list<ChemtableData4D> datalist;   // XXXXXX should be pointers someday for lower-memory
  typedef std::pair<int,double> myPair;
  typedef std::pair<int,int> myPairInt;
  double SMALL_CLIP;

  // shared memory stuff -- Now in core
  //MPI_Comm internode_comm;
  //MPI_Comm intranode_comm;
  //int chemtable_io;         // flags ranks to read/write chemtable


public:

  // constructor
  CartesianChemtable4D(const string& name){
    if ( mpi_rank == 0 )
      cout << "CartesianChemtable4D():" << endl;

    // table dimensions
    n1    = 0;
    n2    = 0;
    n3    = 0;
    n4    = 0;
    nvars = 0;

    // coordinates
    x1 = NULL;
    x2 = NULL;
    x3 = NULL;
    x4 = NULL;

    // array of variable names
    varName = NULL;


    // small value used for clipping
    SMALL_CLIP = 1.01e-6;

    // initialize table
    this->tableName =  name;
    init();

    // test table lookup
    //testLookup();

    //internode_comm = MPI_COMM_NULL;
    //intranode_comm = MPI_COMM_NULL;
    //chemtable_io = 1;
  }

  // deconstructor
  ~CartesianChemtable4D(){
    if ( mpi_rank == 0 )
      cout << "~CartesianChemtable4D():" << endl;

    //writeTecplot("test4D.dat");
    //testLookup();

    info();

    if (x1!=NULL) delete[] x1;
    if (x2!=NULL) delete[] x2;
    if (x3!=NULL) delete[] x3;
    if (x4!=NULL) delete[] x4;
    if (varName!=NULL) delete[] varName;

    // cleanup data pointers...
    for (list<ChemtableData4D>::iterator data = datalist.begin(); data != datalist.end(); ++data) {
      // defaults to delete if no shared memory...
      CTI_Munmap(data->data,data->getSize());
    }

  }



  // ================
  // dump table info
  // ================
  void info(){

    if ( mpi_rank == 0 ) {
      cout << "============== chemistry table information ===============" << endl;
      cout << " chemistry table information: " << endl;
      cout << "  > name:                     " << tableName << endl;
      cout << "  > type:                     " << tableType << endl;
      cout << "  > reference pressure:       " << pressure << endl;
      cout << "  > dimension:                ("<< n1 <<"," << n2 <<"," << n3 << "," << n4 << ")" << endl;
      cout << "   > first coordinate range:   ("<< x1[0] <<"," << x1[n1-1] << ")" << endl;
      cout << "   > second coordinate range:  ("<< x2[0] <<"," << x2[n2-1] << ")" << endl;
      cout << "   > third coordinate range:   ("<< x3[0] <<"," << x3[n3-1] << ")" << endl;
      cout << "   > forth coordinate range:   ("<< x4[0] <<"," << x4[n4-1] << ")" << endl;
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

      for (list<ChemtableData4D>::iterator data = datalist.begin(); data != datalist.end(); ++data) {
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


  // =======================
  // find the variable range
  // =======================
  int getTableVarRange(double& varmin, double& varmax, const string& name){
    // initialze
    varmin = 0.0;
    varmax = 0.0;
    // loop through the list
    for (list<ChemtableData4D>::iterator data = datalist.begin(); data != datalist.end(); ++data) {
      if ( data->name == name ) {
	varmin = (*data)(0,0,0,0);
	varmax = (*data)(0,0,0,0);
	for ( int i=0; i<n1; ++i)
	  for ( int j=0; j<n2; ++j)
	    for ( int k=0; k<n3; ++k)
	      for ( int l=0; l<n4; ++l){
		varmax = max(varmax,(*data)(i,j,k,l));
		varmin = min(varmin,(*data)(i,j,k,l));
	      }
	return(1);
      }
    }
    if ( mpi_rank == 0)
      cout << "***WARNING in getTableVarRange: the variable "
	   << name << " is not loaded yet/n" << " ...setting the range to 0:0 " << endl;

    return(0);

  }


  // ==========================================
  // write the table with the current variables
  // the one that are loaded to a tecplot file
  // ==========================================
  void writeTecplot(const string& filename) {

    // only rank 0 needs to write
    if (mpi_rank == 0) {

      cout << " > writing 4D cartesian chemtable " << tableName << " to the Tecplot file "<< filename << endl;
      cout << "  > the following variables are currently loaded: " << endl;
      for (list<ChemtableData4D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
	cout << "   > " << data->name << endl ;

      FILE * fp = fopen(filename.c_str(),"w");
      fprintf(fp,"title = \"chemtable\"\n");
      fprintf(fp,"variables = \"Z\"\n\"Zvar\"\n\"C\"\n\"H\"\n");
      for (list<ChemtableData4D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
	fprintf(fp,"\"%s\"\n",data->name.c_str());
      fprintf(fp,"zone i = %d, j = %d, k = %d, l = %d, DATAPACKING=POINT\n",n1,n2,n3,n4);
      for (int l=0; l<n4; ++l)
	for (int k=0; k <n3; ++k)
	  for (int j=0; j <n2; ++j)
	    for (int i=0; i <n1; ++i){
	      fprintf(fp,"%lf %lf %lf %lf",x1[i],x2[j],x3[k],x4[l]);
	      for (list<ChemtableData4D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
		fprintf(fp," %lf",(*data)(i,j,k,l));
	      fprintf(fp,"\n");
	    }
      fclose(fp);
    }
  }


private:

  // ===================
  // Linear table lookup
  // ===================
  void lookupDataVec(vector<double*>& routPtrVec, const vector<string>& nameVec,
		     const double* x1, const double* x2, const double* x3,const double* x4, const int n) {
    lookupLinearPtrVector(routPtrVec,nameVec,x1,x2,x3,x4,n);
  }



  // ====================================
  // initialize the table
  // Here we should initialize the types,
  // dimensions, and common variables
  // no major read at this stage
  // ====================================
  void init() {

    //if ( mpi_rank == 0 )
    //  cout << "initialize chemtable: "<< endl;

    // note: open the table and read it on all
    // processors. This is inefficient but since
    // a very small amount of data is involved
    // at this stage, it is not going to make
    // a big difference...


    // read header and coordinates
    readHeaderAndCoordinates();

    // dump table info
    info();

  }

  // ===========================
  // look up vector of variables
  // ===========================
  void lookupLinearPtrVector(vector<double*>& routPtrVec, const vector<string>& nameVec,
		       const double* x1, const double* x2, const double* x3,const double* x4, const int n) {

    vector<ChemtableData4D *> dataPtrVec(nameVec.size());
    for (unsigned int iname = 0; iname < nameVec.size(); ++iname)
      dataPtrVec[iname] = getVar(nameVec[iname]);

    // compute the indecies and integration weights
    // first variable
    int * ti1   = new int[n]; // table index 1 (lower bound)
    double* W1  = new double[n];
    // second variable
    int * ti2   = new int[n];
    double* W2  = new double[n];
    // third variable
    int * ti3   = new int[n];
    double* W3  = new double[n];
    // forth variable
    int * ti4   = new int[n];
    double* W4  = new double[n];
    findIndexAndComputeLinearWeights(ti1,W1,ti2,W2,ti3,W3,ti4,W4,x1,x2,x3,x4,n);

    for ( int i=0; i<n; ++i) {
      // now having the local index, get the table index
      int t1 = ti1[i];
      int t2 = ti2[i];
      int t3 = ti3[i];
      int t4 = ti4[i];
      // get the weights
      double w1 = W1[i];
      double w2 = W2[i];
      double w3 = W3[i];
      double w4 = W4[i];
      // linear interpolation weights
      double c1 = w4*w1*w2*w3;
      double c2 = w4*w1*w2*(1.0-w3);
      double c3 = w4*w1*(1.0-w2)*w3;
      double c4 = w4*w1*(1.0-w2)*(1.0-w3);
      double c5 = w4*(1.0-w1)*w2*w3;
      double c6 = w4*(1.0-w1)*w2*(1.0-w3);
      double c7 = w4*(1.0-w1)*(1.0-w2)*w3;
      double c8 = w4*(1.0-w1)*(1.0-w2)*(1.0-w3);

      double d1 = (1.0-w4)*w1*w2*w3;
      double d2 = (1.0-w4)*w1*w2*(1.0-w3);
      double d3 = (1.0-w4)*w1*(1.0-w2)*w3;
      double d4 = (1.0-w4)*w1*(1.0-w2)*(1.0-w3);
      double d5 = (1.0-w4)*(1.0-w1)*w2*w3;
      double d6 = (1.0-w4)*(1.0-w1)*w2*(1.0-w3);
      double d7 = (1.0-w4)*(1.0-w1)*(1.0-w2)*w3;
      double d8 = (1.0-w4)*(1.0-w1)*(1.0-w2)*(1.0-w3);

      // get the variable and interpolate
      for ( unsigned int j=0; j<routPtrVec.size(); ++j) {
	routPtrVec[j][i] =
	  c1 * (*(dataPtrVec[j]))(t1,t2,t3,t4)     + c2 * (*(dataPtrVec[j]))(t1,t2,t3+1,t4)      +
	  c3 * (*(dataPtrVec[j]))(t1,t2+1,t3,t4)   + c4 * (*(dataPtrVec[j]))(t1,t2+1,t3+1,t4)    +
	  c5 * (*(dataPtrVec[j]))(t1+1,t2,t3,t4)   + c6 * (*(dataPtrVec[j]))(t1+1,t2,t3+1,t4)    +
	  c7 * (*(dataPtrVec[j]))(t1+1,t2+1,t3,t4) + c8 * (*(dataPtrVec[j]))(t1+1,t2+1,t3+1,t4)  +
	  d1 * (*(dataPtrVec[j]))(t1,t2,t3,t4+1)     + d2 * (*(dataPtrVec[j]))(t1,t2,t3+1,t4+1)      +
	  d3 * (*(dataPtrVec[j]))(t1,t2+1,t3,t4+1)   + d4 * (*(dataPtrVec[j]))(t1,t2+1,t3+1,t4+1)    +
	  d5 * (*(dataPtrVec[j]))(t1+1,t2,t3,t4+1)   + d6 * (*(dataPtrVec[j]))(t1+1,t2,t3+1,t4+1)    +
	  d7 * (*(dataPtrVec[j]))(t1+1,t2+1,t3,t4+1) + d8 * (*(dataPtrVec[j]))(t1+1,t2+1,t3+1,t4+1);


      }
    }

    // delete the allocated memory
    delete[] ti1;
    delete[] W1;

    delete[] ti2;
    delete[] W2;

    delete[] ti3;
    delete[] W3;

    delete[] ti4;
    delete[] W4;

  }


  // ===================================
  // for each variable "var" compute the
  // lower bound index in the table and
  // the weight associated with the
  // integration
  // ===================================

  void findIndexAndComputeLinearWeights(int * ti1, double* w1,
					int * ti2, double* w2,
					int * ti3, double* w3,
					int * ti4, double* w4,
					const double* var1,
					const double* var2,
					const double* var3,
					const double* var4,
					const int nv) {
    findIndexAndComputeLinearWeights(ti1,w1,var1,nv,x1,n1);
    findIndexAndComputeLinearWeights(ti2,w2,var2,nv,x2,n2);
    findIndexAndComputeLinearWeights(ti3,w3,var3,nv,x3,n3);
    findIndexAndComputeLinearWeights(ti4,w4,var4,nv,x4,n4);
  }

  void findIndexAndComputeLinearWeights(int * ti1, double* w1,
					int * ti2, double* w2,
					int * ti3, double* w3,
					const double* var1,
					const double* var2,
					const double* var3,
					const int nv) {
    findIndexAndComputeLinearWeights(ti1,w1,var1,nv,x1,n1);
    findIndexAndComputeLinearWeights(ti2,w2,var2,nv,x2,n2);
    findIndexAndComputeLinearWeights(ti3,w3,var3,nv,x3,n3);
  }

  void findIndexAndComputeLinearWeights(int * ti1, double* w1,
					int * ti2, double* w2,
					const double* var1,
					const double* var2,
					const int nv ) {
    findIndexAndComputeLinearWeights(ti1,w1,var1,nv,x1,n1);
    findIndexAndComputeLinearWeights(ti2,w2,var2,nv,x2,n2);
  }

  void findIndexAndComputeLinearWeights(int * ti1, double* w1,
					const double* var1,
					const int nv ) {
    findIndexAndComputeLinearWeights(ti1,w1,var1,nv,x1,n1);
  }

  void findIndexAndComputeLinearWeights(int * ti, double* w,
					const double* var,
					const int nv,
					const double* tvar,
					const int nt) {

    assert(nt>1);

    // pair the variables with the index
    vector<myPair> varPair;
    for ( int i=0; i<nv; ++i )
      varPair.push_back(myPair(i,var[i]));

    // sort the vector of pairs
    lessSecond<int,double> comparator;
    std::sort(varPair.begin(),varPair.end(),comparator);

    // find the index
    int tindex = 0;
    double one_over_length = 1.0 / (tvar[tindex+1] - tvar[tindex]);
    for (vector<myPair>::iterator thisPair = varPair.begin(); thisPair!=varPair.end(); ++thisPair) {
      int i = thisPair->first;
      double thisvar = thisPair->second;
      if (thisvar >= tvar[nt-1]) {
	ti[i] =  nt-2;
	one_over_length = 1.0 / (tvar[nt-1] - tvar[nt-2]);
      }
      else {
	while (thisvar > tvar[tindex+1]) {
	  ++tindex;
	  one_over_length = 1.0 / (tvar[tindex+1] - tvar[tindex]);
	}
	ti[i] = tindex;
      }
      // compute integration weight
      // IMPORTANT: clip the out of bound variavles to the table limites
      // this is ad hoc but it seems we need to clip the
      // variables close to the table limit to the bondaries
      if ( thisvar < (tvar[0]+SMALL_CLIP) )     thisvar = tvar[0];
      if ( thisvar > (tvar[nt-1]-SMALL_CLIP) )  thisvar = tvar[nt-1];
      double tmp = thisvar - tvar[ti[i]];
      w[i] = 1.0 - tmp*one_over_length;
    }
    assert(tindex<(nt-1));

  }

  void loadVariables(const vector<string>& strVec) {
    if (strVec.size() == 0)
      CERR("pointless loadVariables()");

    for (unsigned int istr = 0; istr < strVec.size(); ++istr)
      getVar(strVec[istr]);
  }

  // =======================
  // add a variable to table
  // =======================
  ChemtableData4D * getVar(const string& varname){

    // double check that the variable is not already in the list
    for (list<ChemtableData4D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
      if ( data->name == varname ) return(&(*data));

    // this was moved down into routine...

    //datalist.push_back(ChemtableData4D(varname,n1,n2,n3,n4));
    //ChemtableData4D * thisdata = &(datalist.back());
    //assert(thisdata->name == varname);

    // read the data
    ChemtableData4D * thisdata = readData(varname);

    // if the variable is the source term, clean up
    if ( thisdata->name == "src_prog" )
      cleanupSrcProg(thisdata);

    // return the data pointer
    return(thisdata);

  }


  // ============================
  // delete a variable from table
  // ============================
  void deleteVar(const string& varname){

    // double check that the variable is not already in the list
    for (list<ChemtableData4D>::iterator data = datalist.begin(); data != datalist.end(); ++data)
      if ( data->name == varname ) {
	if ( mpi_rank == 0)
	  cout << " > deleting variable "<< data->name << " from chemtable" << endl;
	CTI_Munmap(data->data,data->getSize());
	datalist.erase(data);
	return;
      }

  }


  // ================================================================
  // cleanup source term such that SRC_PROG is zero at C=0 and C=Cmax
  // ================================================================
  void cleanupSrcProg(ChemtableData4D * data){

    // check the name
    assert(data->name == "src_prog");

    if (mpi_rank == 0 )
      cout << " > cleaning up source term to assure it is zero at C=0 and C=Cmax ..."<<endl;

    // read progress variable
    int table_has_PROG = 0;
    // double check that the variable is not already in the list
    for (list<ChemtableData4D>::iterator thisdata = datalist.begin(); thisdata != datalist.end(); ++thisdata)
      if ( thisdata->name == "prog" ) table_has_PROG = 1;
    ChemtableData4D * Progdata = getVar("prog");

    // set the source term to zero if c<SMALL_CLIP or c>Cmax
    double max_SrcProg_zero = 0.0;
    if (mpi_rank_shared == 0) {
      for ( int i=0; i<n1; ++i)
	for ( int j=0; j<n2; ++j)
	  for ( int l=0; l<n4; ++l) {
	    double Cmax = (*Progdata)(i,j,n3-1,l);
	    for ( int k=0; k<n3; ++k)
	      if ( (x3[k]<SMALL_CLIP) || (x3[k]>=Cmax) || ((*data)(i,j,k,l)<0.0) ) {
		max_SrcProg_zero = max(max_SrcProg_zero,fabs((*data)(i,j,k,l)));
		(*data)(i,j,k,l) = 0.0;
	      }
	  }
      if (mpi_rank == 0)
	cout << "   > maximum value set to zero = " << max_SrcProg_zero <<endl;
    }

    // delete data
    if ( !table_has_PROG  ) deleteVar("prog");

  }

  // ===========================
  // read header and coordinates
  // ===========================
  void readHeaderAndCoordinates() {


    // read on all the ranks for now

    // open the file and go to data
    FILE * fp = NULL;
    // open the file
    const int file_err = MiscUtils::openFile(&fp,tableName,"rb");
    if (file_err != 0) throw(file_err);
    // fp = fopen(tableName.c_str(), "rb");
    // if (!fp) {
    //   cerr << "cannot open chemtable file: " << tableName << endl;
    //   throw(0);
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
      if ( header.id == UGP_IO_CT_CART_4D ) {
	// type of the table
	tableType = header.name;
	// table dimensions
	n1    = header.idata[0]; assert(n1>1);
	n2    = header.idata[1]; assert(n2>1);
	n3    = header.idata[2]; assert(n3>1);
	n4    = header.idata[3]; assert(n4>1);
	nvars = header.idata[4];
	// allocate space for table coordinates
	assert(x1==NULL);      x1 = new double[n1];
	assert(x2==NULL);      x2 = new double[n2];
	assert(x3==NULL);      x3 = new double[n3];
	assert(x4==NULL);      x4 = new double[n4];
	assert(varName==NULL); varName = new string[nvars];
	// reference quantities
	pressure   = header.rdata[0]; assert(pressure>0.0);
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
	else if ( name == "Zvar" ) {
	  assert( header.idata[0] == n2 );
	  assert(x2!=NULL);
	  tmpPtr = x2;
	  ++icoor;
	}
	else if ( name == "C" ) {
	  assert( header.idata[0] == n3 );
	  assert(x3!=NULL);
	  tmpPtr = x3;
	  ++icoor;
	}
	else if ( name == "Enth" ) {
	  assert( header.idata[0] == n4 );
	  assert(x4!=NULL);
	  tmpPtr = x4;
	  ++icoor;
	}
        else if ( name == "Temp" ) {
          assert( header.idata[0] == n4 );
          assert(x4!=NULL);
          tmpPtr = x4;
          ++icoor;
        }
	else{
	  cerr << ERRSTART << "Error: unknown coordinate " << header.name << ERREND << endl;
	  throw(0);
	}
	// read the data
	size_t dummy= fread(&(tmpPtr[0]),sizeof(double),header.idata[0],fp);
	assert(int(dummy)==header.idata[0]);
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
	cerr << ERRSTART << "Error: could not find one of the three coordinates or the header " << ERREND << endl;
	throw(0);
      }

      // increase the offset
      offset += header.skip;

    } // while (done != 1)

    // close the file
    fclose(fp);

    assert(ivar==nvars);
    assert(icoor==4);
  }

  ChemtableData4D * readData(const string& varname) {

    // use this to broadcast an error to all processors
    int error_rank = 0;

    size_t data_size = size_t(n1)*n2*n3*n4;
    double * dataptr = NULL;
    CTI_Mmap(dataptr,data_size);

    // only rank 0 reads the data...

    if (mpi_rank == 0) {

      cout << " > reading " << varname << " from chemtable " << tableName << ": ";

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
	  cerr << ERRSTART << "Error: could not find variable " << varname << " in the chemtable" << ERREND << endl;
	  error_rank = 1;
	}
	// compare the names
	string name = header.name;
	if ( name == varname ) {
	  // check compatibility
	  assert( header.id       == UGP_IO_CT_DATA );
	  assert( header.idata[0] == (n1*n2*n3*n4) );
	  assert( header.idata[1] == n1 );
	  assert( header.idata[2] == n2 );
	  assert( header.idata[3] == n3 );
	  assert( header.idata[4] == n4 );

	  assert(int(data_size) == header.idata[0]);
	  // read the data
	  size_t dummy= fread(dataptr, sizeof(double), header.idata[0], fp);
	  assert(int(dummy)==header.idata[0]);
	  // byte swap
	  if (byte_swap)
	    ByteSwap::byteSwap(dataptr,header.idata[0]);
	  // dump variable range
	  double varmin = dataptr[0];
	  double varmax = dataptr[0];
	  for (int ii = 0; ii < header.idata[0]; ++ii) {
	    varmax = max(varmax, dataptr[ii]);
	    varmin = min(varmin, dataptr[ii]);
	  }
	  cout << varmin << " < "<< varname <<" < " << varmax << endl ;
	  // set the done
	  done = 1;
	}

	// increase the offset
	offset += header.skip;

      } // while (done != 1)

      // close the file
      if ( error_rank == 0 )
	fclose(fp);

    }

    // send the error check from rank zero to all the processors
    {
      int err;
      MPI_Allreduce(&error_rank, &err, 1, MPI_INT, MPI_MAX, mpi_comm);
      if (err == 1)
	CERR("Problem reading file " << tableName);
    }

    // broadcast to other nodes
    if (mpi_rank_shared == 0) {
      MPI_Bcast(dataptr, data_size, MPI_DOUBLE, 0, mpi_comm_internode);
    }

    MPI_Barrier(mpi_comm);

    // init chemtabledata with mapped shared memory...
    datalist.push_back(ChemtableData4D(varname, dataptr, n1, n2, n3, n4));
    return(&(datalist.back()));
  }

  /*
  // no longer used
  // ======================================
  // read the data associated with a header
  // ======================================
  ChemtableData4D * readData(const string& varname) {

    // add the variable to the list
    // this needs to be done on all the processors

    datalist.push_back(ChemtableData4D(varname,n1,n2,n3,n4));
    ChemtableData4D * thisdata = &(datalist.back());
    assert(thisdata->name == varname);

    // use this to broadcast an error to all processors
    int error_rank = 0;

    // read only on rank zero
    if ( mpi_rank == 0 ) {

      cout << " > reading " << thisdata->name << " from chemtable " << tableName << ": ";

      // open the file and go to data
      FILE * fp = NULL;
      // open the file
      fp = fopen(tableName.c_str(), "rb");
      if (!fp) {
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
	  cerr << ERRSTART << "Error: could not find variable " << thisdata->name << " in the chemtable" << ERREND << endl;
	  error_rank = 1;
	}
	// compare the names
	string name = header.name;
	if ( name == thisdata->name ) {
	  // check compatibility
	  assert( header.id       == UGP_IO_CT_DATA );
	  assert( header.idata[0] == (n1*n2*n3*n4) );
	  assert( header.idata[1] == n1 );
	  assert( header.idata[2] == n2 );
	  assert( header.idata[3] == n3 );
	  assert( header.idata[4] == n4 );
	  // read the data
	  size_t dummy= fread(&((*thisdata)(0,0,0,0)),sizeof(double),header.idata[0],fp);
	  assert(dummy==header.idata[0]);
	  // byte swap
	  if (byte_swap)
	    ByteSwap::byteSwap(&((*thisdata)(0,0,0,0)),header.idata[0]);
	  // dump variable range
	  double varmin = (*thisdata)(0,0,0,0);
	  double varmax = (*thisdata)(0,0,0,0);
	  for ( int i=0; i<n1; ++i)
	    for ( int j=0; j<n2; ++j)
	      for ( int k=0; k<n3; ++k)
		for ( int l=0; l<n4; ++l){
		  varmax = max(varmax,(*thisdata)(i,j,k,l));
		  varmin = min(varmin,(*thisdata)(i,j,k,l));
	      }
	  cout << varmin << " < "<< thisdata->name <<" < " << varmax << endl ;
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

    MPI_Bcast(&((*thisdata)(0,0,0,0)),thisdata->getSize(),MPI_DOUBLE,0,mpi_comm);

    return(thisdata);
  }
  */


public:




  // =================
  // test table lookup
  // =================
  void testLookup() {

    // Do it on rank zero only


    // set the clipping to zero
    double SMALL_CLIP_orig = SMALL_CLIP;
    SMALL_CLIP = 0.0;

    if ( mpi_rank == 0 )
      cout << "testing table lookup ..." << endl;
    // pick two variables in the table
    ChemtableData4D * var1 = getVar("T");
    if ( mpi_rank == 0 )
      cout << " > replacing T with linear function 2*Z + 3*Zv + 4*C + 5*h"<< endl;
    ChemtableData4D * var2 = getVar("rho");
    if ( mpi_rank == 0 )
      cout << " > replacing rho with linear function 20.1*Z + 31.0*Zv + 29.0001*C + 0.201*h"<< endl;
    // set them with linear functions
    for ( int i=0; i<n1; ++i)
      for ( int j=0; j<n2; ++j)
	for ( int k=0; k<n3; ++k)
	  for ( int l=0; l<n4; ++l) {
	    (*var1)(i,j,k,l) = 2.0*x1[i] + 3.0*x2[j] + 4.0*x3[k] + 5.0*x4[l] ;
	    (*var2)(i,j,k,l) = 20.1*x1[i] - 31.0*x2[j] + 29.0001*x3[k] + 0.201*x4[l];
	  }

    // Build a random array
    int myN = 100000;
    double *myZ  = new double[myN];
    double *myZv = new double[myN];
    double *myC  = new double[myN];
    double *myh  = new double[myN];
    //srand(10);
    for ( int i=0; i<myN; ++i ) {
      myZ[i]  = x1[0] + (x1[n1-1]-x1[0])*(double)(rand())/(double)(RAND_MAX);
      myZv[i] = x2[0] + (x2[n2-1]-x2[0])*(double)(rand())/(double)(RAND_MAX);
      myC[i]  = x3[0] + (x3[n3-1]-x3[0])*(double)(rand())/(double)(RAND_MAX);
      myh[i]  = x4[0] + (x4[n4-1]-x4[0])*(double)(rand())/(double)(RAND_MAX);
    }

    // result vector
    double *rout1 = new double[myN];
    double *rout2 = new double[myN];

    // do lookup
    lookup(rout1, "T", rout2, "rho", myZ,myZv,myC,myh,myN);

    // compute error
    double myerror1 = 0.0;
    for ( int i=0; i<myN; ++i)
      myerror1 = max(fabs(rout1[i]-(2.0*myZ[i] + 3.0*myZv[i] + 4.0*myC[i] + 5.0*myh[i]  )),myerror1);
    double myerror2 = 0.0;
    for ( int i=0; i<myN; ++i)
      myerror2 = max(fabs(rout2[i]-(20.1*myZ[i] - 31.0*myZv[i] + 29.0001*myC[i] + 0.201*myh[i])),myerror2);
    // take the max across processors
    double error1;
    MPI_Reduce(&myerror1,&error1,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    double error2;
    MPI_Reduce(&myerror2,&error2,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    if ( mpi_rank == 0 )
      cout << "  >> error in lookups of the two linear functions (should be zero) = " << error1 << " , " << error2 << endl;

    // set the clipping back to original one
    SMALL_CLIP = SMALL_CLIP_orig;

    // deletle variables from table
    deleteVar("T");
    deleteVar("rho");

    // delete memory
    delete[] myZ;
    delete[] myZv;
    delete[] myC;
    delete[] myh;
    delete[] rout1;
    delete[] rout2;

    if ( mpi_rank == 0 )
      cout << "done with table lookup test" << endl;

  }
};

// ===================
// delete chemtable
// ===================
void deleteChemtable(AbstractChemtable1D* &chemtable);
void deleteChemtable(AbstractChemtable3D* &chemtable);
void deleteChemtable(AbstractChemtable4D* &chemtable);


// =======================
// initialize chemtable 1D
// =======================
inline void initChemtable(AbstractChemtable1D * &chemtable, const string& tablename) {


  string tabletype = getChemtableType(tablename);

  COUT1(" > initializing 1D table: " << tablename << " with the type: " << tabletype);

  if ( (tabletype == "VIDA_PREMIXED_FPV_CART1D"))

    chemtable = new CartesianChemtable1D(tablename);

  else
    CERR("incompatible 1D table type " <<  tabletype);

}




// =======================
// initialize chemtable 2D
// =======================
inline void initChemtable(AbstractChemtable2D * &chemtable, const string& tablename) {


  string tabletype = getChemtableType(tablename);

  COUT1(" > initializing 2D table: " << tablename << " with the type: " << tabletype);

  if ( (tabletype == "VIDA_PREMIXED_FPV_CART2D") ||
       (tabletype == "CHARLES_PREMIXED_FPV_CART2D")|| (tabletype == "PREMIXED"))

    chemtable = new CartesianChemtable2D(tablename);

  else
    CERR("incompatible 2D table type " <<  tabletype);

}


// =======================
// initialize chemtable 3D
// =======================
inline void initChemtable(AbstractChemtable3D * &chemtable, const string& tablename) {


  string tabletype = getChemtableType(tablename);

  COUT1(" > initializing 3D table: " << tablename << " with the type: " << tabletype);

  if ( (tabletype == "VIDA_NON-PREMIXED_FPV_CART3D") ||
       (tabletype == "CHRIS_NON-PREMIXED_FPV_CART3D") ||
       (tabletype == "NON-PREMIXED") || 
       (tabletype == "NONPREMIXED")) {

    chemtable = new CartesianChemtable3D(tablename, tabletype);
  }
  else if ( (tabletype == "VIDA_NON-PREMIXED_FPV_KDT3D") ||
	    (tabletype == "CHRIS_NON-PREMIXED_FPV_KDT3D") ) {

    chemtable = new KdChemtable3D(tablename, tabletype);
  }
  else
    CERR("incompatible 3D table type " <<  tabletype);

}



// =======================
// initialize chemtable 4D
// =======================
inline void initChemtable(AbstractChemtable4D * &chemtable, const string& tablename) {


  string tabletype = getChemtableType(tablename);

  COUT1(" > initializing 4D table: " << tablename << " with the type: " << tabletype);

  if ( (tabletype == "VIDA_NON-PREMIXED_FPVH_CART4D") ||
       (tabletype == "VIDA_NON-PREMIXED_FPVT_CART4D") ) {

    chemtable = new CartesianChemtable4D(tablename);

  }
  else if ( (tabletype == "VIDA_NON-PREMIXED_FPVH_KDT4D") ){

    CERR("initialization for 4D kdt type " <<  tabletype << " is not implemented yet");
  }
  else
    CERR("incompatible 4D table type " <<  tabletype);

}

// ===================
// delete chemtable 2D
// ===================
inline void deleteChemtable(AbstractChemtable2D * &chemtable) {
  if ( chemtable != NULL ) delete chemtable;
}

// ===================
// delete chemtable 3D
// ===================
inline void deleteChemtable(AbstractChemtable3D * &chemtable) {
  if ( chemtable != NULL ) delete chemtable;
}

#endif
