#ifndef ABSTRACTCHEMTABLE_HPP
#define ABSTRACTCHEMTABLE_HPP

#include "Common.hpp"
#include "ArrayNd.hpp"
#include "ByteSwap.hpp"
#include "MiscUtils.hpp"
#include "CTI.hpp"

using namespace CTI;

// set the fseek to zero and go to the offset position
inline void my_fseek(FILE * fp, int8 offset) {
  fseek(fp,0,SEEK_SET);
  fseek(fp,offset,SEEK_CUR);
}

inline string getChemtableType(const string& filename){

  // open the file and go to data
  FILE* fp = NULL;
  const int file_err = MiscUtils::openFile(&fp,filename,"rb");
  if (file_err != 0) {
    CERR("cannot open chemtable file: " << filename );
  }

  // read the magic number and IO version and figure out the byte swap
  int byte_swap = 0;
  {
    int itmp[2];
    fread(itmp,sizeof(int),2,fp);
    if (itmp[0] != CHEM_IO_MAGIC_NUMBER) {
      ByteSwap::byteSwap(itmp,2);
      if (itmp[0] != CHEM_IO_MAGIC_NUMBER)
	CERR("file " << filename << " does not start as expected. aborting.");
      byte_swap = 1;
    }
    if (itmp[1] != CHEM_IO_VERSION)
      CERR("in file " << filename << "io version differs from current version: ");
  }
  // initialize the offset...
  int8 offset = sizeof(int)*2;
  Header header;
  // find the data and read it
  int done  = 0;
  while ( done!=1 ) {
    // go to the offset position
    my_fseek(fp,offset);
    // read the ehader
    fread(&header,sizeof(Header),1,fp);
    if (byte_swap)
      ByteSwap::byteSwapHeader(&header,1);
    // compare the header id
    if ( (header.id == UGP_IO_CT_CART_1D) || (header.id == UGP_IO_CT_CART_2D) || (header.id == UGP_IO_CT_CART_3D) ||  (header.id == UGP_IO_CT_KD) || (header.id == UGP_IO_CT_CART_4D) )
      break;
    if ( header.id == UGP_IO_EOF )
      CERR("could not find table type");
    // increase the offset
    offset += header.skip;
  } // while (done != 1)

  // close the file
  fclose(fp);
  return(header.name);
}


template <typename tp1, typename tp2>
class lessSecond {
  typedef std::pair<tp1,tp2> myPair;
public:
  bool operator() ( const myPair& lhs , const myPair& rhs ) {
    return(lhs.second < rhs.second);
  }
};

class AbstractChemtable1D {

public:
  string tableName;
  string tableType;
  double pressure;

  AbstractChemtable1D() {
    // initialize string info...
    tableName   = "None";
    tableType   = "None";
    // reference quantities
    pressure    = -1.0;
  }

  virtual ~AbstractChemtable1D() {};
  virtual void info() = 0;
  virtual double getFlameThickness() = 0;
  virtual double getLaminarFlameSpeed() = 0;
  virtual void writeTecplot(const string& filename) = 0;
  virtual void loadVariables(const vector<string>& strVec) = 0;
  virtual void getRhoCsrcRange(double& Smin, double& Smax) = 0;
  virtual int getTableVarRange(double& varmin, double& varmax, const string& name) = 0;

private:
  virtual void lookupDataVec(vector<double*>& routPtrVec, const vector<string>& nameVec,
		             const double* x1,const int n) = 0;
public:

  void lookup(vector<double*>& routPtrVec, const vector<string>& nameVec,
	      const double* x1, const int n) {
    lookupDataVec(routPtrVec,nameVec,x1,n);
  }

  void lookup(double* rout1,const string& name1,
	      const double* x1, const int n) {
    // build a vector pointer of the data and names
    vector<double*> routPtrvec(1);
    vector<string> namePtrvec(1);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      const double* x1, const int n) {
    // build a vector pointer of the data and names
    vector<double*> routPtrvec(2);
    vector<string> namePtrvec(2);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      const double* x1, const int n) {
    // build a vector pointer of the data and names
    vector<double*> routPtrvec(3);
    vector<string> namePtrvec(3);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      const double* x1, const int n) {
    // build a vector pointer of the data and names
    vector<double*> routPtrvec(4);
    vector<string> namePtrvec(4);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      const double* x1, const int n) {
    // build a vector pointer of the data and names
    vector<double*> routPtrvec(5);
    vector<string> namePtrvec(5);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      const double* x1, const int n) {
    // build a vector pointer of the data and names
    vector<double*> routPtrvec(6);
    vector<string> namePtrvec(6);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,n);
  }
};

class AbstractChemtable2D {

public:

  string tableName;
  string tableType;
  double  pressure;

  AbstractChemtable2D() {
    // initialize string info...
    tableName   = "None";
    tableType   = "None";
    // reference quantities
    pressure    = -1.0;
  }

  virtual ~AbstractChemtable2D() {};
  virtual void info() = 0;
  virtual void writeTecplot(const string& filename) = 0;
  virtual void loadVariables(const vector<string>& strVec) = 0;
  virtual void getRhoCsrcRange(double& Smin, double& Smax) = 0;
  virtual int getTableVarRange(double& varmin, double& varmax, const string& name) = 0;

private:

  virtual void initProgvarBounds() = 0;
  virtual void lookupDataVec(vector<double*>& routPtrVec, const vector<string>& nameVec,
		             const double* x1, const double* x2, const int n) = 0;

  virtual void lookupDataVecReduced(vector<double*>& routPtrVec, const vector<string>& nameVec,
                                    const double * x1, const int n) = 0;

  virtual void lookupDataVecReducedNew(vector<double*>& routPtrVec, const vector<string>& nameVec,
                                    const double *__restrict__ x1,
                                    const int (*__restrict__ cvofa)[2], const int n) = 0;

public:

  virtual void lookupSpecialNew(double *__restrict__ rout, const string& name,
                               const double* __restrict__ y1, const double *__restrict__ y2,
                               const int (*__restrict__ cvofa)[2], const int n1, const int n2) = 0;

  virtual void lookupSpecialNew(double *__restrict__ rout1, const string& name1,
                                double *__restrict__ rout2, const string& name2,
                               const double* __restrict__ y1, const double *__restrict__ y2,
                               const int (*__restrict__ cvofa)[2], const int n1, const int n2) = 0;

  virtual void lookupSpecial(double* rout, const string& name,
                             const double *__restrict__ y1, const double *__restrict__ y2a,
                             const double *__restrict__ y2b, const int n) = 0;

  void lookup(vector<double*>& routPtrVec, const vector<string>& nameVec,
		const double* x1, const double* x2, const int n) {
    lookupDataVec(routPtrVec,nameVec,x1,x2,n);
  }

  void lookup(double* rout1,const string& name1,
	      const double* x1,const double* x2,
	      const int n) {
    // build a vector pointer of the data and names
    vector<double*> routPtrvec(1);
    vector<string> namePtrvec(1);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,n);
  }

  void lookupReduced(double* rout1,const string& name1,
	      const double* x1, const int n) {
    // build a vector pointer of the data and names
    vector<double*> routPtrvec(1);
    vector<string> namePtrvec(1);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1;
    // lookup
    lookupDataVecReduced(routPtrvec,namePtrvec,x1,n);
  }

  void lookup(double* rout1,const string& name1,
              double* rout2,const string& name2,
	      const double* x1,const double* x2,
	      const int n) {
    // build a vector pointer of the data and names
    vector<double*> routPtrvec(2);
    vector<string> namePtrvec(2);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,n);
  }

  void lookupReduced(double* rout1,const string& name1,
                     double* rout2,const string& name2,
                     const double* x1,const int n) {
    // build a vector pointer of the data and names
    vector<double*> routPtrvec(2);
    vector<string> namePtrvec(2);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2;
    // lookup
    lookupDataVecReduced(routPtrvec,namePtrvec,x1,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      const double* x1,const double* x2,
	      const int n) {
    // build a vector pointer of the data and names
    vector<double*> routPtrvec(3);
    vector<string> namePtrvec(3);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,n);
  }

  void lookupReduced(double* rout1,const string& name1,
	             double* rout2,const string& name2,
	             double* rout3,const string& name3,
                     const double* x1, const int n) {
    // build a vector pointer of the data and names
    vector<double*> routPtrvec(3);
    vector<string> namePtrvec(3);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3;
    // lookup
    lookupDataVecReduced(routPtrvec,namePtrvec,x1,n);
  }

  void lookupReducedNew(double* rout1, const string& name1,
                        double* rout2, const string& name2,
                        const double* x1, const int (*cvofa)[2], 
                        const int n) { 
    
    // build a vector pointer of the data and names
    vector<double*> routPtrvec(2);
    vector<string> namePtrvec(2);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2;
    // lookup 
    lookupDataVecReducedNew(routPtrvec,namePtrvec,x1,cvofa,n);
  }
  
  void lookupReducedNew(double* rout1, const string& name1,
                        double* rout2, const string& name2,
                        double* rout3, const string& name3,
                        const double* x1, const int (*cvofa)[2],
                        const int n) {

    // build a vector pointer of the data and names
    vector<double*> routPtrvec(3);
    vector<string> namePtrvec(3);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3;
    // lookup
    lookupDataVecReducedNew(routPtrvec,namePtrvec,x1,cvofa,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      const double* x1,const double* x2,
	      const int n) {
    // build a vector pointer of the data and names
    vector<double*> routPtrvec(4);
    vector<string> namePtrvec(4);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      const double* x1,const double* x2,
	      const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(5);
    vector<string> namePtrvec(5);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      const double* x1,const double* x2,
	      const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(6);
    vector<string> namePtrvec(6);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      const double* x1,const double* x2,
	      const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(7);
    vector<string> namePtrvec(7);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      double* rout8,const string& name8,
	      const double* x1,const double* x2,
	      const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(8);
    vector<string> namePtrvec(8);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7; ++i;
    routPtrvec[i] = rout8; namePtrvec[i] = name8;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      double* rout8,const string& name8,
	      double* rout9,const string& name9,
	      const double* x1,const double* x2,
	      const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(9);
    vector<string> namePtrvec(9);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7; ++i;
    routPtrvec[i] = rout8; namePtrvec[i] = name8; ++i;
    routPtrvec[i] = rout9; namePtrvec[i] = name9;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      double* rout8,const string& name8,
	      double* rout9,const string& name9,
	      double* rout10,const string& name10,
	      const double* x1,const double* x2,
	      const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(10);
    vector<string> namePtrvec(10);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7; ++i;
    routPtrvec[i] = rout8; namePtrvec[i] = name8; ++i;
    routPtrvec[i] = rout9; namePtrvec[i] = name9; ++i;
    routPtrvec[i] = rout10; namePtrvec[i] = name10;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      double* rout8,const string& name8,
	      double* rout9,const string& name9,
	      double* rout10,const string& name10,
	      double* rout11,const string& name11,
	      const double* x1,const double* x2,
	      const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(11);
    vector<string> namePtrvec(11);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7; ++i;
    routPtrvec[i] = rout8; namePtrvec[i] = name8; ++i;
    routPtrvec[i] = rout9; namePtrvec[i] = name9; ++i;
    routPtrvec[i] = rout10; namePtrvec[i] = name10; ++i;
    routPtrvec[i] = rout11; namePtrvec[i] = name11;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      double* rout8,const string& name8,
	      double* rout9,const string& name9,
	      double* rout10,const string& name10,
	      double* rout11,const string& name11,
	      double* rout12,const string& name12,
	      const double* x1,const double* x2,
	      const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(12);
    vector<string> namePtrvec(12);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7; ++i;
    routPtrvec[i] = rout8; namePtrvec[i] = name8; ++i;
    routPtrvec[i] = rout9; namePtrvec[i] = name9; ++i;
    routPtrvec[i] = rout10; namePtrvec[i] = name10; ++i;
    routPtrvec[i] = rout11; namePtrvec[i] = name11; ++i;
    routPtrvec[i] = rout12; namePtrvec[i] = name12;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,n);
  }

  virtual bool varExists(const string& var_name) const = 0;
};

class AbstractChemtable3D {

public:

  string tableName;
  string tableType;
  double  pressure;

  AbstractChemtable3D() {
    // initialize string info...
    tableName   = "None";
    tableType   = "None";
    // reference quantities
    pressure    = -1.0;
  }

  virtual ~AbstractChemtable3D() {};
  virtual void info() = 0;
  virtual void writeTecplot(const string& filename) = 0;
  virtual void writeTecplot2(const string& filename) = 0;
  virtual int getTableVarRange(double& varmin, double& varmax, const string& name) = 0;
  virtual void loadVariables(const vector<string>& strVec) = 0;

private:

  virtual void lookupDataVec(vector<double*>& routPtrVec, const vector<string>& nameVec,
		     const double* x1, const double* x2, const double* x3, const int n) = 0;

public:

  void lookup(vector<double*>& routPtrVec, const vector<string>& nameVec,
		     const double* x1, const double* x2, const double* x3, const int n) {
    lookupDataVec(routPtrVec,nameVec,x1,x2,x3,n);
  }

  void lookup(double* rout1,const string& name1,
	      const double* x1,const double* x2,
	      const double* x3,const int n) {
    // build a vector pointer of the data and names
    vector<double*> routPtrvec(1);
    vector<string> namePtrvec(1);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      const double* x1,const double* x2,
	      const double* x3,const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(2);
    vector<string> namePtrvec(2);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      const double* x1,const double* x2,
	      const double* x3,const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(3);
    vector<string> namePtrvec(3);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      const double* x1,const double* x2,
	      const double* x3,const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(4);
    vector<string> namePtrvec(4);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      const double* x1,const double* x2,
	      const double* x3,const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(5);
    vector<string> namePtrvec(5);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      const double* x1,const double* x2,
	      const double* x3,const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(6);
    vector<string> namePtrvec(6);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      const double* x1,const double* x2,
	      const double* x3,const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(7);
    vector<string> namePtrvec(7);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      double* rout8,const string& name8,
	      const double* x1,const double* x2,
	      const double* x3,const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(8);
    vector<string> namePtrvec(8);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7; ++i;
    routPtrvec[i] = rout8; namePtrvec[i] = name8;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      double* rout8,const string& name8,
	      double* rout9,const string& name9,
	      const double* x1,const double* x2,
	      const double* x3,const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(9);
    vector<string> namePtrvec(9);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7; ++i;
    routPtrvec[i] = rout8; namePtrvec[i] = name8; ++i;
    routPtrvec[i] = rout9; namePtrvec[i] = name9;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      double* rout8,const string& name8,
	      double* rout9,const string& name9,
	      double* rout10,const string& name10,
	      const double* x1,const double* x2,
	      const double* x3,const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(10);
    vector<string> namePtrvec(10);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7; ++i;
    routPtrvec[i] = rout8; namePtrvec[i] = name8; ++i;
    routPtrvec[i] = rout9; namePtrvec[i] = name9; ++i;
    routPtrvec[i] = rout10; namePtrvec[i] = name10;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      double* rout8,const string& name8,
	      double* rout9,const string& name9,
	      double* rout10,const string& name10,
	      double* rout11,const string& name11,
	      const double* x1,const double* x2,
	      const double* x3,const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(11);
    vector<string> namePtrvec(11);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7; ++i;
    routPtrvec[i] = rout8; namePtrvec[i] = name8; ++i;
    routPtrvec[i] = rout9; namePtrvec[i] = name9; ++i;
    routPtrvec[i] = rout10; namePtrvec[i] = name10; ++i;
    routPtrvec[i] = rout11; namePtrvec[i] = name11;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      double* rout8,const string& name8,
	      double* rout9,const string& name9,
	      double* rout10,const string& name10,
	      double* rout11,const string& name11,
	      double* rout12,const string& name12,
	      const double* x1,const double* x2,
	      const double* x3,const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(12);
    vector<string> namePtrvec(12);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7; ++i;
    routPtrvec[i] = rout8; namePtrvec[i] = name8; ++i;
    routPtrvec[i] = rout9; namePtrvec[i] = name9; ++i;
    routPtrvec[i] = rout10; namePtrvec[i] = name10; ++i;
    routPtrvec[i] = rout11; namePtrvec[i] = name11; ++i;
    routPtrvec[i] = rout12; namePtrvec[i] = name12;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,n);
  }

  void lookupReduced( double* rout1,const string& name1, 
                      const double* x1, const int n) {
    
    //build a vector pointer of the data and names
  
    vector<double*> routPtrvec(1);
    vector<string> namePtrvec(1);
    routPtrvec[0] = rout1; namePtrvec[0] = name1;
    lookupDataVecReduced(routPtrvec,namePtrvec,x1,n);
  
  }

  virtual void lookupDataVecReduced(vector<double*>& routPtrVec, const vector<string>& nameVec,
                                    const double * x1, const int n) = 0;

  virtual bool varExists(const string& var_name) const = 0;


};

class AbstractChemtable4D {

public:

  string tableName;
  string tableType;
  double  pressure;

  AbstractChemtable4D() {
    // initialize string info...
    tableName   = "None";
    tableType   = "None";
    // reference quantities
    pressure    = -1.0;
  }

  virtual ~AbstractChemtable4D() {};
  virtual void info() = 0;
  virtual void writeTecplot(const string& filename) = 0;
  virtual int getTableVarRange(double& varmin, double& varmax, const string& name) = 0;
  virtual void loadVariables(const vector<string>& strVec) = 0;

private:

  virtual void lookupDataVec(vector<double*>& routPtrVec,
			     const vector<string>& nameVec,
			     const double* x1, const double* x2,
			     const double* x3, const double* x4,
			     const int n) = 0;

public:

  void lookup(vector<double*>& routPtrVec, const vector<string>& nameVec,
	      const double* x1, const double* x2,
	      const double* x3, const double* x4, const int n) {
    lookupDataVec(routPtrVec,nameVec,x1,x2,x3,x4,n);
  }

  void lookup(double* rout1,const string& name1,
	      const double* x1,const double* x2,
	      const double* x3,const double* x4, const int n) {
    // build a vector pointer of the data and names
    vector<double*> routPtrvec(1);
    vector<string> namePtrvec(1);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,x4,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      const double* x1,const double* x2,
	      const double* x3,const double* x4, const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(2);
    vector<string> namePtrvec(2);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,x4,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      const double* x1,const double* x2,
	      const double* x3,const double* x4, const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(3);
    vector<string> namePtrvec(3);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,x4,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      const double* x1,const double* x2,
	      const double* x3,const double* x4, const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(4);
    vector<string> namePtrvec(4);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,x4,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      const double* x1,const double* x2,
	      const double* x3,const double* x4, const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(5);
    vector<string> namePtrvec(5);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,x4,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      const double* x1,const double* x2,
	      const double* x3,const double* x4, const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(6);
    vector<string> namePtrvec(6);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,x4,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      const double* x1,const double* x2,
	      const double* x3,const double* x4, const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(7);
    vector<string> namePtrvec(7);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,x4,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      double* rout8,const string& name8,
	      const double* x1,const double* x2,
	      const double* x3,const double* x4, const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(8);
    vector<string> namePtrvec(8);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7; ++i;
    routPtrvec[i] = rout8; namePtrvec[i] = name8;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,x4,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      double* rout8,const string& name8,
	      double* rout9,const string& name9,
	      const double* x1,const double* x2,
	      const double* x3,const double* x4, const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(9);
    vector<string> namePtrvec(9);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7; ++i;
    routPtrvec[i] = rout8; namePtrvec[i] = name8; ++i;
    routPtrvec[i] = rout9; namePtrvec[i] = name9;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,x4,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      double* rout8,const string& name8,
	      double* rout9,const string& name9,
	      double* rout10,const string& name10,
	      const double* x1,const double* x2,
	      const double* x3,const double* x4, const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(10);
    vector<string> namePtrvec(10);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7; ++i;
    routPtrvec[i] = rout8; namePtrvec[i] = name8; ++i;
    routPtrvec[i] = rout9; namePtrvec[i] = name9; ++i;
    routPtrvec[i] = rout10; namePtrvec[i] = name10;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,x4,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      double* rout8,const string& name8,
	      double* rout9,const string& name9,
	      double* rout10,const string& name10,
	      double* rout11,const string& name11,
	      const double* x1,const double* x2,
	      const double* x3,const double* x4, const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(11);
    vector<string> namePtrvec(11);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7; ++i;
    routPtrvec[i] = rout8; namePtrvec[i] = name8; ++i;
    routPtrvec[i] = rout9; namePtrvec[i] = name9; ++i;
    routPtrvec[i] = rout10; namePtrvec[i] = name10; ++i;
    routPtrvec[i] = rout11; namePtrvec[i] = name11;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,x4,n);
  }

  void lookup(double* rout1,const string& name1,
	      double* rout2,const string& name2,
	      double* rout3,const string& name3,
	      double* rout4,const string& name4,
	      double* rout5,const string& name5,
	      double* rout6,const string& name6,
	      double* rout7,const string& name7,
	      double* rout8,const string& name8,
	      double* rout9,const string& name9,
	      double* rout10,const string& name10,
	      double* rout11,const string& name11,
	      double* rout12,const string& name12,
	      const double* x1,const double* x2,
	      const double* x3,const double* x4, const int n) {
    // build a vector pointer of the data and output
    vector<double*> routPtrvec(12);
    vector<string> namePtrvec(12);
    // set the vectors
    int i = 0;
    routPtrvec[i] = rout1; namePtrvec[i] = name1; ++i;
    routPtrvec[i] = rout2; namePtrvec[i] = name2; ++i;
    routPtrvec[i] = rout3; namePtrvec[i] = name3; ++i;
    routPtrvec[i] = rout4; namePtrvec[i] = name4; ++i;
    routPtrvec[i] = rout5; namePtrvec[i] = name5; ++i;
    routPtrvec[i] = rout6; namePtrvec[i] = name6; ++i;
    routPtrvec[i] = rout7; namePtrvec[i] = name7; ++i;
    routPtrvec[i] = rout8; namePtrvec[i] = name8; ++i;
    routPtrvec[i] = rout9; namePtrvec[i] = name9; ++i;
    routPtrvec[i] = rout10; namePtrvec[i] = name10; ++i;
    routPtrvec[i] = rout11; namePtrvec[i] = name11; ++i;
    routPtrvec[i] = rout12; namePtrvec[i] = name12;
    // lookup
    lookupDataVec(routPtrvec,namePtrvec,x1,x2,x3,x4,n);
  }
};

//==========================
// cartesian data classes
//==========================

class ChemtableData1D : public CTIarray<double> {
public:
  string name;
  ChemtableData1D(const string& _name, const int n1) : CTIarray<double>(n1), name(_name) {}
  ~ChemtableData1D(){}
};

class ChemtableData2D : public CTIarray<double> {
public:
  string name;
  ChemtableData2D(const string& _name, const int n1, const int n2) : CTIarray<double>(n1,n2), name(_name) {}
  ~ChemtableData2D(){}
};

class ChemtableData3D : public CTIarray<double> {
public:
  string name;
  ChemtableData3D(const string& _name, const int n1, const int n2, const int n3) :
    CTIarray<double>(n1,n2,n3), name(_name) {}
  ~ChemtableData3D(){}
};

class ChemtableData4D : public CTIarray<double> {
public:
  string name;
  ChemtableData4D(const string& _name, const int n1, const int n2, const int n3, const int n4) :
    CTIarray<double>(n1,n2,n3,n4), name(_name) {}
  ChemtableData4D(const string& _name, double * ptr, const int n1, const int n2, const int n3, const int n4) :
    CTIarray<double>(ptr,n1,n2,n3,n4), name(_name) {}
  ~ChemtableData4D(){}
};



//============================
// chemtable IO
//============================

class SerialIO {
public:

  int byte_swap;
  FILE * fp;
  string filename;

  SerialIO(const string& filename) {
    COUT3("SerialIO()");
    this->filename = filename;

    // open file for reading
    fp = NULL;
    const int file_err = MiscUtils::openFile(&fp,filename,"rb");
    if (file_err != 0) throw(file_err);
    // fp = fopen(filename.c_str(),"rb");
    // if (fp == NULL) CERR("cannot open chemtable file:" << filename);

    // read magic number and version number and set byteswap
    byte_swap = 0;
    {
      int itmp[2];
      size_t size_read = fread(itmp,sizeof(int),2,fp);
      assert(size_read == 2);

      if (itmp[0] != CHEM_IO_MAGIC_NUMBER) {
	ByteSwap::byteSwap(itmp,2);
	if (itmp[0] != CHEM_IO_MAGIC_NUMBER) {
	  cerr << ERRSTART << "Error: file " << filename << " does not start as expected. aborting." << ERREND << endl;
	  throw(0);
	}
	byte_swap = 1;
      }
      if (itmp[1] != CHEM_IO_VERSION) {
	cerr << ERRSTART << "Error: in file " << filename << "io version differs from current version: " << itmp[1] << ERREND << endl;
	throw(0);
      }
    }
  }

  ~SerialIO() {
    COUT3("~SerialIO()");
    if (fp != NULL) fclose(fp);
  }

  int getHeaderFromFile(Header &header, int8 &offset, const string data_name, const int data_id) {

    assert(fp != NULL);
    offset = sizeof(int)*2;

    // scan file for desired header
    while (1) {
      my_fseek(fp,offset);

      // read the header
      size_t size_read = fread(&header,sizeof(Header),1,fp);
      if (byte_swap) ByteSwap::byteSwapHeader(&header,1);
      assert(size_read == 1);
      if (header.id == data_id && !data_name.compare(0,data_name.length(), header.name) ) {
	offset += sizeof(Header);
	return(0);
      }
      else if (header.id == UGP_IO_EOF) {
	return(1);
      }
      offset += header.skip;
    }
  }

  template <class T>
  void getDataFromFileAtOffset(T * data, const int8 offset, const int ndata, const string data_name,
			       const int data_id) {
    assert(fp != NULL);
    my_fseek(fp,offset);
    size_t size_read = fread(data, sizeof(T), ndata, fp);
    assert(int(size_read) == ndata);
    if (byte_swap) ByteSwap::byteSwap(data, ndata);
  }
};


class ParallelIO {
public:

  SerialIO * rank0_file;
  MPI_Datatype MPI_HEADER;
  string filename;

  ParallelIO(const string& filename) {
    COUT3("ParallelIO()");
    this->filename = filename;

    // build header datatype - this is done in Ugp.cpp and lives in a namespace (can we remove this?)...
    MPI_Datatype oldtypes[5];
    int blockcounts[5];
    MPI_Aint offsets[5];

    // the name...
    offsets[0] = 0;
    oldtypes[0] = MPI_CHAR;
    blockcounts[0] = UGP_IO_HEADER_NAME_LEN;

    // the id...
    offsets[1] = offsets[0] + 1*UGP_IO_HEADER_NAME_LEN;
    oldtypes[1] = MPI_INT;
    blockcounts[1] = 1;

    // the offset skip...
    offsets[2] = offsets[1] + 4;
    oldtypes[2] = MPI_OFFSET_DATATYPE;
    blockcounts[2] = 1;

    // the 16 ints...
    offsets[3] = offsets[2] + 8;
    oldtypes[3] = MPI_INT;
    blockcounts[3] = 16;

    // the 16 reals...
    offsets[4] = offsets[3] + 4*16;
    oldtypes[4] = MPI_DOUBLE;
    blockcounts[4] = 16;

    // build the header type...
    MPI_Type_create_struct(5, blockcounts, offsets, oldtypes, &MPI_HEADER);
    MPI_Type_commit(&MPI_HEADER);

    rank0_file = NULL;
  }

  ~ParallelIO() {
    COUT3("~ParallelIO()");
    if (rank0_file != NULL) delete rank0_file;
    MPI_Type_free(&MPI_HEADER);
  }

  void getHeaderFromFileAndBCast(Header &header, const string data_name, const int data_id) {

    if (mpi_rank == 0) assert(rank0_file != NULL);
    int8 offset;
    int rank0_error;
    if (mpi_rank == 0)
      rank0_error = rank0_file->getHeaderFromFile(header, offset, data_name, data_id);

    MPI_Bcast(&rank0_error, 1, MPI_INT, 0, mpi_comm);
    if (rank0_error != 0)
      CERR("data " << data_name << " with id " << data_id << " not found in file " << filename);

    MPI_Bcast(&offset, 1, MPI_INT8, 0, mpi_comm);
    MPI_Bcast(&header, 1, MPI_HEADER, 0, mpi_comm);
  }


  template <class T>
  void getDataFromFileAndBCast(T * data, const int ndata, const string data_name,
			       const int data_id, const MPI_Datatype mpi_type) {

    if (mpi_rank == 0) assert(rank0_file != NULL);
    Header header;
    int8 offset = 0;
    int rank0_error;

    if (mpi_rank == 0)
      rank0_error = rank0_file->getHeaderFromFile(header, offset, data_name, data_id);

    MPI_Bcast(&rank0_error, 1, MPI_INT, 0, mpi_comm);
    if (rank0_error != 0)
      CERR("data " << data_name << " with id " << data_id << " not found in file " << filename);

    if (mpi_rank == 0) {
      assert(header.idata[0] == ndata);
      rank0_file->getDataFromFileAtOffset(data, offset, ndata, data_name, data_id);
    }

    MPI_Bcast(data, ndata, mpi_type, 0, mpi_comm);
  }
};

class KdChemtableSerialIO : public SerialIO {
public:

  vector<string> dataNameVec;

  KdChemtableSerialIO(const string& filename) : SerialIO(filename) {
    // fill vector will all the datanames in file

    // scan file until we hit EOF
    int8 offset = sizeof(int)*2;
    while (1) {
      my_fseek(fp,offset);

      // read the header
      Header header;
      size_t size_read = fread(&header, sizeof(Header), 1, fp);
      assert(size_read == 1);
      if (byte_swap) ByteSwap::byteSwapHeader(&header,1);
      if (header.id == UGP_IO_CT_DATA) {
	dataNameVec.push_back(header.name);
      }
      else if (header.id == UGP_IO_EOF) {
	break;
      }
      offset += header.skip;
    }
  }
};

class KdChemtableParallelIO : public ParallelIO {
public:

  KdChemtableSerialIO* serial_ptr;

  KdChemtableParallelIO(const string& filename) : ParallelIO(filename) {
    assert(rank0_file == NULL);
    serial_ptr = NULL;
    // only open file on rank 0...
    if (mpi_rank == 0) {
      serial_ptr = new KdChemtableSerialIO(filename);
      rank0_file = serial_ptr;
    }
    else {
      serial_ptr = NULL;
    }
  }

  ~KdChemtableParallelIO(){
    if ( serial_ptr != NULL ) {
      delete serial_ptr;
      rank0_file = NULL;
    }
  }
};

#endif
