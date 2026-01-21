#ifndef MISCUTILS_HPP
#define MISCUTILS_HPP

#include "Common.hpp"
#include "MpiStuff.hpp"

using namespace MpiStuff;
using namespace std;

/**
 * MiscUtils
 *
 * misc stuff that comes in handy from time to time...
 *
 */
namespace MiscUtils {

  inline void applyRotation(double out[3],const double R[9],const double in[3]) {
    FOR_I3 {
      out[i] = 0.0;
      FOR_J3 out[i] += R[i*3+j]*in[j];
    }
  }

  inline void applyInvRotation(double out[3],const double R[9],const double in[3]) {
    FOR_I3 {
      out[i] = 0.0;
      FOR_J3 out[i] += R[i+j*3]*in[j];
    }
  }

  inline void matVecMult(double out[3],const double R[9],const double in[3]) {
    FOR_I3 {
      out[i] = 0.0;
      FOR_J3 out[i] += R[i*3+j]*in[j];
    }
  }

  inline void vecVecAdd(double out[3],const double t[3],const double in[3]) {
    FOR_I3 out[i] = t[i]+in[i];
  }

  inline void invertMat(double A_inv[3][3],const double A[3][3]) {
    const double det = DETERMINANT(A); assert(det != 0.0);
    const double inv_det = 1.0/det;
    A_inv[0][0] = (A[1][1] * A[2][2] - A[2][1] * A[1][2]) * inv_det;
    A_inv[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) * inv_det;
    A_inv[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * inv_det;
    A_inv[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) * inv_det;
    A_inv[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * inv_det;
    A_inv[1][2] = (A[1][0] * A[0][2] - A[0][0] * A[1][2]) * inv_det;
    A_inv[2][0] = (A[1][0] * A[2][1] - A[2][0] * A[1][1]) * inv_det;
    A_inv[2][1] = (A[2][0] * A[0][1] - A[0][0] * A[2][1]) * inv_det;
    A_inv[2][2] = (A[0][0] * A[1][1] - A[1][0] * A[0][1]) * inv_det;
  }

  void nullTerminate(char *name);
  void nullTerminate(char *name, int len);
  string toLowerCase(const string word);
  string toUpperCase(const string word);

  void calcUniformDist(int * &xod, const int nx, const int ndist);
  void calcUniformDist(int8 * &xod, const int8 nx, const int ndist);
  void calcUniformDistAlign8(int8 * &xod, const int8 nx, const int ndist);
  //void calcUniformDist(int64_t * &xod, const int64_t nx, const int ndist);
  void calcUniformDistNoNew(int8 * xod, const int8 nx, const int ndist);

  void calcThresholdDist(int * &xod, const int nx, const int ndist, const int nmin) ;
  void calcThresholdDist(int8 * &xod, const int8 nx, const int ndist, const int nmin) ;

  void buildUniformXora(int * &Xora, const int nX_global);
  void buildUniformXora(int8 * &Xora, const int8 nX_global);
  //void buildUniformXora(int64_t * &Xora, const int64_t nX_global);
  void buildXora(int * &Xora,int nX_local);
  void buildXora(int8 * &Xora,int nX_local);
  void buildXora(int8 * &Xora,int8 nX_i8);

  void copyXora(int * &xora,const int8 * xora_global);
  int getRankInXora(const int ix,const int * xora_global);
  int getRankInXora(const int8 ix,const int8 * xora_global);
  //int getRankInXora(const double x,const double * xora_global);
  int getIntervalBisection(const double value,const double * const valueArray,const int n);
  int getIdxInXoso(const int ix,const int * xoso,const int nn);

  void dumpRange(const int * s, const int n, const string& message);
  void dumpRange(const int8 * s, const int n, const string& message);
  void dumpRange(const int (*v)[3], const int n, const string& message);
  void dumpRange(const double * s, const int n, const string& message,MPI_Comm& comm = mpi_comm);
  void dumpRange(const double (*v)[3], const int n, const string& message,MPI_Comm& comm = mpi_comm);
  void dumpRange(const float * s, const int n, const string& message,MPI_Comm& comm = mpi_comm);
  void dumpRange(const int value, const string& message);

  void setRange(const double (*v)[3], const int n, double (&buf)[3][2]);
  void dumpRange(const double (*t)[3][3], const int n, const string& message);

  // historgram in text output
  void dumpHist0(const double * var,const int n,const string& message);
  void dumpHist0(const double * var,const int n,const int nbins,const string& message,vector<double>&,vector<double>&,bool);

  void dumpSum(const int * s, const int n, const string& message);
  void dumpSum(const double * s, const int n, const string& message);
  void dumpSum(const double (*v)[3], const int n, const string& message);

  void dumpBins(const int * s, const int n, const string& message);

  void dumpNodeAvgReport(double value,const string& message);
  void dumpNodeSumReport(double value,const string& message);

  void reorder_csr(int * x_i, int * x_v, int * order, const int n);
  void reorder_csr(int * x_i, int * x_v1, int * x_v2, int * order, const int n);
  void reorder(double * var,const int * order,const int n);
  void reorder(double (*var)[3],const int * order,const int n);
  void reorder(double (*var)[3][3],const int * order,const int n);
  void reorder(int * v,const int * order,const int n);
  void reorder(int8 * v,const int * order,const int n);
  void reorder(int (*v2)[2],const int * order,const int n);

  /**
   * checks if a token is in string
   * \param string
   * \return 1 ... if token on lhs
   * \return 0 ... if no token
   * \return -1 ... if no token is within non-whitespace string
   */
  int checkHash(const string &str);
  bool parseFunction(string& function,string& args,const string& str);

  /**
   * takes a string and seperates it into tokens
   * \return tokens stored in vector !!!
   * \param str ... processed string
   * \param str ... delimiters, standard delimiter is whitespace
   */
  void tokenizeString(vector<string> &tokens, const string &str,
                      const string &delimiters = " ");

  void splitCsv(vector<string>& nameVec,const string& names);
  void splitCsv(set<string>& nameSet,const string& names);

  bool splitCsv(vector<uint2>& uint2Vec,const string& csv);

  bool splitCsv(vector<int>& intVec,const string& csv);
  bool splitCsv(set<int>& intSet,const string& names);

  bool splitCsv(vector<double>& doubleVec,const string& csv);
  bool splitCsv(set<double>& doubleSet,const string& names);

  string getFileNameExtension(const string &str, const string &delimiters);

  void buildTilePrefix(char * tilePrefix, const char * prefix, const int tile[3]);
  void buildUnindexedFilename(char * filename,const char * prefix,const char * suffix);
  void buildIndexedFilename(char * filename,const char * prefix,const int index,const char * suffix);
  void buildIndexedDirAndFilename(char * filename,const char * dir_prefix,const int dir_index,const int file_index,const char * file_suffix);

  bool fileExists(const string& filename);
  bool fileExists(const char * filename);
  int getIndexFromFilename(const string& filename);

  int openFile(FILE ** fp,const string& filename,const string& mode);
  int read2DAsciiTable(double (* &xy_vals)[2], const string filename,const int skip=1);
  int read3DAsciiTable(double (* &vals)[3], const string filename,const int skip=1);

  // resizing utilities...
  void resize(bool * &scalar,const int n_new);
  void resize(bool * &scalar,const int n_old,const int n_new);

  void resize(int * &scalar,const int n_new);
  void resize(int * &scalar,const int n_old,const int n_new);

  void resize(int8 * &scalar,const int n_new);
  void resize(int8 * &scalar,const int n_old,const int n_new);

  void resize(uint8 * &scalar,const int n_new);
  void resize(uint8 * &scalar,const int n_old,const int n_new);

  void resize(int (* &scalar)[2],const int n_new);
  void resize(int (* &scalar)[2],const int n_old,const int n_new);

  void resize(int (* &scalar)[4],const int n_new);
  void resize(int (* &scalar)[8],const int n_new);

  void resize(int (* &scalar)[4],const int n_old,const int n_new);
  void resize(int (* &scalar)[4][2],const int n_old,const int n_new);

  void resize(int (* &scalar)[6],const int n_old,const int n_new);
  void resize(int (* &scalar)[6][4],const int n_old,const int n_new);

  void resize(int (* &scalar)[8],const int n_old,const int n_new);
  void resize(int (* &scalar)[8][4],const int n_old,const int n_new);

  void resize(double * &scalar,const int n_new);
  void resize(double * &scalar,const int n_old,const int n_new);

  void resize(double (* &d3)[3],const int n_new);
  void resize(double (* &d3)[3],const int n_old,const int n_new);

  void resize(double (* &d3x3)[3][3],const int n_new);
  void resize(double (* &d3x3)[3][3],const int n_old,const int n_new);

  // random functions...

  double uniformRand(const double rmin,const double rmax);
  double randn();

  // minimum of a range...
  int rangeMin(const int * values,const int n);
  int rangeMax(const int * values,const int n);

  void removeMean(double * var,const double * norm,const int n);
  double calcMean(double * var,const double * norm,const int n);

  void newArray2d(double ** &data,const int ni,const int nj);
  double ** newArray2d(const int ni,const int nj);
  void deleteArray2d(double ** &data);
  void transposeArray2d(double ** pt,const double * const * p,const int * iora,const int * jora);

  // complex versions of above...
  void newComplexArray2d(complex<double> ** &data,const int ni,const int nj);
  complex<double> ** newComplexArray2d(const int ni,const int nj);
  void deleteComplexArray2d(complex<double> ** &data);

  void mkdir_for_file(const string& filename);
  void mkdir_for_file_collective(const string& filename, const int write_rank=0);

  string makeTmpPrefix(const string& prefix);
  string getFilenameNoPath(const string& filename_plus_path);
  string addPrefixToFilename(const string& prefix,const string& filename);
  string getPrefix(const string& filename);
  void eraseAllSubStr(string& mainStr, const string& toErase);

  // void replaceOldSubStringWithNewSubString(string& str,const string& oldSubString,const string& newSubString);

  int gaussj(double a[][3], int n, double b[][3], int m);

  bool isNan(const double& value);

  // some helpful comparisons for sorts...
  bool compareFirstSecondIntIntPair(const std::pair<int,int>& a,const std::pair<int,int>& b);
  bool compareSecondIntDoublePair(const std::pair<int,double>& a,const std::pair<int,double>& b);

  // returns true if bigstring starts with smallstring...
  bool startsWith(const string& bigString,const string& smallString);
  bool strcmp_wildcard(const string& str,const string& str_wildcard);
  bool regexMatchesString(const string& regex,const string& str);

  inline string getArg(const string& name, const string& func) {
    size_t pos = 0, pos2 = 0;
    string token;
    if ( (pos = name.find(func+"(")) != string::npos) {
      string token = name.substr(pos);
      pos = token.find("(");
      pos2 = token.find(")");
      if ( (pos != string::npos) && (pos2 != string::npos) && (pos2 >= pos+1)) {
        token = token.substr(pos+1,pos2-pos-1);
        return token;
      }
    }
    return "";
  }

  /*
    inline double fast_pow(const double x,const double p) {
    const double y = x - 1.0;
    const double xx = 1.0 + p * y * (1.0 + (p - 1.0) * y * (0.5 + (1.0 / 6.0) * (p - 2.0) * y));
    if ( abs(y) < 1.0e-04) {
    return xx;
    } else {
    return pow(x,p);
    }
    }
  */

  inline double fast_pow_posx(const double x, const double p) {
    const double pp = p*log(x);
    return 1.0 + expm1(pp);
  }

  // search
  template<class T>
  inline void findInterval( int &r, const T val, const T* arr, int n ) {
    unsigned int ilo = 0 ;
    unsigned int ihi = n-1;
    while ( ihi != ilo+1 ) {
      unsigned int imi = (ihi+ilo) / 2 ;
      if ( arr[imi] > val )
        ihi = imi ;
      else if ( arr[imi] == val) {
        ilo = imi ;
        ihi = imi+1;
      }
      else
        ilo = imi ;
    }
    r = int(ilo) ;
  }

  // sorting ..
  void quickSortByValue( int * idx, double * vals, int n) ;
  void quickSortByValue( int * idx, int *    vals, int n) ;

  // eigenvalue stuff...
  void eigenDecomposition(double A[3][3], double V[3][3], double d[3]);
  void transpose(double A[3][3]);

  // bit stuff
  void setBit(int&i,int ibit);
  void clearBit(int&i,int ibit);
  int setNextBit(int &i);
  int clearLastSetBit(int& i);
  int nSetBits(const int& i);

  void getBestE1E2FromE0(double e1[3],double e2[3],const double e0[3]);

  // geometry...
  double calcAngle(const double x0[3],const double x1[3],const double x2[3]);
  double getCosAngleBetweenTris(const double v00[3],const double v01[3],const double v02[3],const double v10[3],const double v11[3],const double v12[3]);
  double getEdgeToEdgeDist2(const double v00[3],const double v01[3], const double v10[3],const double v11[3]);
  double getPointToLineDist2(const double xp[3],const double x0[3],const double n0[3]);
  double getPointToLineDist(const double xp[3],const double x0[3],const double n0[3]);
  double getPointToEdgeDist2(const double xp[3],const double v0[3],const double v1[3]);
  double getPointToPlaneDist(const double xp[3],const double x0[3],const double n0[3]);
  void getClosestPointOnEdgeRobust(double xc[3],const double xp[3],const double v0[3],const double v1[3]);
  void getPointToEdgeDir(double dir[3],const double xp[3],const double v0[3],const double v1[3]);
  double getPointToTriDist2(const double xp[3],const double v0[3],const double v1[3],const double v2[3]);
  void getTriUnitNormal(double (&un)[3],const double v0[3],const double v1[3],const double v2[3]);
  void getClosestPointOnTriRobust(double xc[3],const double xp[3],const double v0[3],const double v1[3],const double v2[3],const bool debug = false);
  void getClosestPointOnTriInterior(double xc[3],const double xp[3],const double v0[3],const double v1[3],const double v2[3]);
  bool getClosestPointOnTri(double xc[3],const double xp[3],const double v0[3],const double v1[3],const double v2[3]);
  double spacingFunction(const double xi,const double stretch,const double x0,const double x1);
  double spacingFunctionDoubleSided(const double xi,const double stretch,const double x0,const double x1);
  void getOrthogonalVectors(double (&radDir1)[3],double (&radDir2)[3],const double xc[3],const double np[3], const double hint[3]);
  void createCirclePts(double (* xp)[3],const int index0,const double xc[3],const double np[3],const double radius, const int n,const bool stagger=true);
  void createEllipsePts(double (* xp)[3],const int index0,const double xc[3],const double np[3],const double rM,const double rm, const double primary_dir[3],const int n);
  void create3DFrom2DPtAndNormal(double (*xp)[3],const int index0,const double xc[3],const double np[3],double (*xy_vals)[2],const int n);
  int lineSphereIntersections(float (&intersections)[2][3], const float line_xc[3], const float line_n[3], const float sphere_xc[3], const float sphere_r);
  bool getBarycentricCoordinates(double lambda[4],const double xp[3],const double x0[3],
                                 const double x1[3], const double x2[3],const double x3[3]);

  // 1d cubic interp
  inline void cubicInterpWts( double w[4], double *x, const double xi ) {
    const double d0=xi-x[0];
    const double d1=xi-x[1];
    const double d2=xi-x[2];
    const double d3=xi-x[3];

    const double d01=x[0]-x[1];
    const double d02=x[0]-x[2];
    const double d03=x[0]-x[3];
    const double d12=x[1]-x[2];
    const double d13=x[1]-x[3];
    const double d23=x[2]-x[3];

    w[0]= d1*d2*d3/(d01*d02*d03);
    w[1]=-d2*d3*d0/(d01*d12*d13);
    w[2]= d3*d0*d1/(d02*d12*d23);
    w[3]=-d0*d1*d2/(d03*d13*d23);
  }

  inline void cubicInterpWts(double w[4], double inv_denom[4], double *x, const double xi) {
    const double d0=xi-x[0];
    const double d1=xi-x[1];
    const double d2=xi-x[2];
    const double d3=xi-x[3];
    w[0]=d1*d2*d3*inv_denom[0];
    w[1]=d2*d3*d0*inv_denom[1];
    w[2]=d3*d0*d1*inv_denom[2];
    w[3]=d0*d1*d2*inv_denom[3];
  }

  inline void cubicInterpDenom(double inv_denom[4], double *x) {
    const double d01=x[0]-x[1];
    const double d02=x[0]-x[2];
    const double d03=x[0]-x[3];
    const double d12=x[1]-x[2];
    const double d13=x[1]-x[3];
    const double d23=x[2]-x[3];
    inv_denom[0]= 1.0/(d01*d02*d03);
    inv_denom[1]=-1.0/(d01*d12*d13);
    inv_denom[2]= 1.0/(d02*d12*d23);
    inv_denom[3]=-1.0/(d03*d13*d23);
  }

  // fast index-finding
  class indexMap {

  public:

    int nloc;
    double *xloc, dxloc;
    double *iloc, diloc;
    double xmin, xmax, xsmall;

  public:

    indexMap(): nloc(0), xloc(NULL), iloc(NULL) {}

    indexMap(double *x, int nx) {

      // construct uniform grid for fast index-finding
      xmin = fmin(x[0],x[nx-1]);
      xmax = fmax(x[0],x[nx-1]);
      xsmall = (xmax-xmin)*1e-12;
      double dxmin = xmax-xmin;
      for (int i=0; i<nx-1; ++i)
        dxmin = fmin(dxmin,x[i+1]-x[i]);

      // xloc needs at least 2 points per interval in x
      nloc = (int)ceil(2.0*(xmax-xmin)/dxmin)+1;
      xloc = new double[nloc];
      iloc = new double[nloc];
      dxloc = (xmax-xmin)/double(nloc-1);
      for (int i=0; i<nloc; ++i)
        xloc[i] = xmin+i*dxloc;

      // locate index for each x and store result in iloc
      // (xloc,iloc) should be monotonic and pass through each (x,i)
      diloc = 1.0/dxloc;
      for (int j=0; j<nloc; ++j) {
        int ix;
        double w;
        for (ix=0; ix<nx-1; ++ix)
          if (x[ix+1]>=xloc[j]) break;
        ix = min(ix,nx-2);
        w = fmin(1.0,fmax(0.0,(x[ix+1]-xloc[j])/(x[ix+1]-x[ix])));
        iloc[j] = w*ix+(1.0-w)*(ix+1);
      }
      double dimin=nx;
      for (int i=0; i<nloc-1; ++i)
        dimin = fmin(dimin,iloc[i+1]-iloc[i]);
      int *ix = new int[nx];
      ix[0] = 0;
      for (int j=1; j<nx-1; ++j) {
        int i;
        for (i=0; i<nloc-1; ++i)
          if (xloc[i+1]>=x[j]) break;
        ix[j]=i;
      }
      ix[nx-1]=nloc;

      // correct iloc to pass through each (x,i)
      iloc[0]=0.0;
      for (int j=1; j<nx-1; ++j) {
        iloc[ix[j]  ]=double(j)+dimin*diloc*(xloc[ix[j]  ]-x[j]);
        iloc[ix[j]+1]=double(j)+dimin*diloc*(xloc[ix[j]+1]-x[j]);
      }
      iloc[nloc-1]=double(nx-1);
      delete[] ix;
    }

    ~indexMap() {
      delete[] xloc;
      delete[] iloc;
    }

    inline int findIndex(const double x) const {
      const double xi = fmin(xmax-xsmall,fmax(xmin+xsmall,x));
      const int j = (int)floor((xi-xmin)*diloc);
      return (int)floor(iloc[j]+(iloc[j+1]-iloc[j])*(xi-xloc[j])*diloc);
    }

    inline void findIndex(int *ix, double *x, const int n) const {
      for (int i=0; i<n; ++i) {
        const double xi = fmin(xmax-xsmall,fmax(xmin+xsmall,x[i]));
        const int j = (int)floor((xi-xmin)*diloc);
        ix[i] = (int)floor(iloc[j]+(iloc[j+1]-iloc[j])*(xi-xloc[j])*diloc);
      }
    }
  };

  // this parallel routine returns the index of the closest point on the rank that
  // owns the closest point, and -1 everywhere else...
  int getClosestPoint(const double xp_check[3],const double (* const xp)[3],const int np);

  // get plane coefficients from three points (i.e., triangle)
  void computePlaneCoeffsFromPoints(double (&N0)[3], double& d0, const double v0[3],const double v1[3],const double v2[3]);

  // // get projected point on a plane
  // void projectPointToPlane(double(&x_proj)[3],const double N0[3], const double d0) {
  //
  // }

  int computeLinePlaneIntersection(double (&x_intersection)[3],const double xl[3], const double nl[3],const double xp[3],const double np[3]);

  // determine how a line segment (a,b) pierces a triangle (v0,v1,v2); type of penetration indicated by return int value
  int howDoesLinePierceTri(const double v0[3],const double v1[3],const double v2[3],const double a[3],const double b[3],const double dist_tol=0.0);

  // replace all occurances of "find" by "replace" in "str" and return it
  std::string replaceAll(const std::string str ,const std::string find ,const std::string replace, const size_t start=0);  // default is to process entire string

  std::string replaceLast(const std::string str ,const std::string find ,const std::string replace);

  // clean OS/filesystem protected chars from filenames so file can be created/accessed
  std::string cleanFilename(const std::string& str, const size_t start=0);

  // compile-time endianness (see Defs.hpp)
  int getEndianness();

  // simple serial solvers...
  int solveCgSerial(double * phi,const int ncv,const double * const A,const double * const rhs,
                    const int * const cvocv_i,const int * const cvocv_v,
                    const double zero,const int maxiter,const bool verbose);

  int solveGaussSidelSerial(double * phi,const int ncv,const double * const A,const double * const rhs,
                            const int * const cvocv_i,const int * const cvocv_v,
                            const double zero,const int maxiter,const double relax,const bool verbose);

  // fast extraction of a column of data from a text file, with grep and tail options...
  // column indexing is "gnuplot" style: 1,2,3. These routines are designed to be fast and
  // do not use catch/throw error management. Instead, the calling process should look for
  // return val != 0, and this should be very rare.
  int xcol(double * &data,int &n,char const *fname,const int column);
  int xcol2(double * &data1,double * &data2,int &n,char const *fname,const int column1,const int column2);
  int xcol4(double * &data1,double * &data2,double * &data3,double * &data4,
	    int &n,char const *fname,
	    const int column1,const int column2,const int column3,const int column4);
  int grepxcol(double * &data,int &n,char const *fname,char const *expr,const int column);
  int grepxcol2(double * &data1,double * &data2,int &n,char const *fname,char const *expr,const int column1,const int column2);
  int tailxcol(double * &data,int &n,char const *fname,const int column);
  int tailxcol2(double * &data1,double * &data2,int &n,char const *fname,const int column1,const int column2);
  int tailgrepxcol(double * &data,int &n,char const *fname,char const *expr,const int column);
  int tailgrepxcol2(double * &data1,double * &data2,int &n,char const *fname,char const *expr,const int column1,const int column2);

  // ======================================================================
  // ROBUST (hopefully) double precision tri-tri intersection routines...
  // ======================================================================

  enum TriTriIntKind {
    MATCHED_INT,
    NODE_NODE_INT,
    EDGE_TRI_INT,
    TRI_EDGE_INT,
    EDGE_EDGE_INT,
    EDGE_NODE_INT,
    NODE_EDGE_INT,
  };

  class TriTriInt {
  public:

    TriTriIntKind kind;
    int idata[2];
    double ddata[2];
    string str[2];
    int match_count;

    TriTriInt(const TriTriIntKind kind,const int i0,const int i1) {
      assert(kind == NODE_NODE_INT);
      this->kind = kind;
      idata[0] = i0;
      idata[1] = i1;
      ddata[0] = HUGE_VAL;
      ddata[1] = HUGE_VAL;
      match_count = 0;
    }

    TriTriInt(const TriTriIntKind kind,const int i,const double wgt,const string& wgt_str) {
      assert((kind == EDGE_TRI_INT)||(kind == TRI_EDGE_INT));
      this->kind = kind;
      idata[0] = i;
      idata[1] = -1;
      ddata[0] = wgt;
      ddata[1] = HUGE_VAL;
      str[0] = wgt_str;
      match_count = 0;
    }

    TriTriInt(const TriTriIntKind kind,const int i0,const int i1,const double wgt0,const string& wgt0_str,const double wgt1,const string& wgt1_str) {
      assert(kind == EDGE_EDGE_INT);
      this->kind = kind;
      idata[0] = i0;
      idata[1] = i1;
      ddata[0] = wgt0;
      ddata[1] = wgt1;
      str[0] = wgt0_str;
      str[1] = wgt1_str;
      match_count = 0;
    }

    TriTriInt(const TriTriIntKind kind,const int i0,const int i1,const double wgt,const string& wgt_str) {
      assert((kind == EDGE_NODE_INT)||(kind == NODE_EDGE_INT));
      this->kind = kind;
      idata[0] = i0;
      idata[1] = i1;
      ddata[0] = wgt; // NOTE: edge wgt goes in ddata[0] for both EDGE_NODE_INT and NODE_EDGE_INT
      ddata[1] = HUGE_VAL;
      str[0] = wgt_str;
      match_count = 0;
    }

    void getX(double x[3],double *x0[3],double *x1[3]) const {

      switch (kind) {
      case NODE_NODE_INT:
        {
          double xi0[3]; FOR_I3 xi0[i] = x0[idata[0]][i];
          double xi1[3]; FOR_I3 xi1[i] = x1[idata[1]][i];
          if (DIST(xi0,xi1) > 1.0E-12) {
            cout << "Warning: NODE_NODE_INT: DIST(xi0,xi1) > 1.0E-12: " << DIST(xi0,xi1) << endl;
            //assert(0);
          }
        }
        FOR_I3 x[i] = 0.5*(x0[idata[0]][i] + x1[idata[1]][i]);
        return;
      case EDGE_TRI_INT:
        // idata[0] contains the index of the edge in tri0...
        // store wgt wrt the head of the edge so sorted intersections
        // will go from tail to head...
        FOR_I3 x[i] = ddata[0]*x0[(idata[0]+1)%3][i] + (1.0-ddata[0])*x0[idata[0]][i];
        return;
      case TRI_EDGE_INT:
        // idata[0] contains the index of the edge in tri1...
        FOR_I3 x[i] = ddata[0]*x1[(idata[0]+1)%3][i] + (1.0-ddata[0])*x1[idata[0]][i];
        return;
      case EDGE_EDGE_INT:
        {
          double xi0[3]; FOR_I3 xi0[i] = ddata[0]*x0[(idata[0]+1)%3][i] + (1.0-ddata[0])*x0[idata[0]][i];
          double xi1[3]; FOR_I3 xi1[i] = ddata[1]*x1[(idata[1]+1)%3][i] + (1.0-ddata[1])*x1[idata[1]][i];
          if (DIST(xi0,xi1) > 1.0E-12) {
            cout << "Warning: EDGE_EDGE_INT: DIST(xi0,xi1) > 1.0E-12: " << DIST(xi0,xi1) << endl;

            //assert(0);
          }
        }
        FOR_I3 x[i] = 0.5*( ddata[0]*x0[(idata[0]+1)%3][i] + (1.0-ddata[0])*x0[idata[0]][i] +
                            ddata[1]*x1[(idata[1]+1)%3][i] + (1.0-ddata[1])*x1[idata[1]][i] );
        return;
      case NODE_EDGE_INT:
        {
          double xi1[3]; FOR_I3 xi1[i] = ddata[0]*x1[(idata[1]+1)%3][i] + (1.0-ddata[0])*x1[idata[1]][i];
          if (DIST(x0[idata[0]],xi1) > 1.0E-12) {
            cout << "Warning: NODE_EDGE_INT: DIST(x0[idata[0]],xi1) > 1.0E-12: " << DIST(x0[idata[0]],xi1) << endl;
            //assert(0);
          }
          cout << "NODE_EDGE_INT dist looks good!" << DIST(x0[idata[0]],xi1) << endl;
        }
        // return the node...
        FOR_I3 x[i] = x0[idata[0]][i];
        return;
      case EDGE_NODE_INT:
        {
          double xi0[3]; FOR_I3 xi0[i] = ddata[0]*x0[(idata[0]+1)%3][i] + (1.0-ddata[0])*x0[idata[0]][i];
          if (DIST(xi0,x1[idata[1]]) > 1.0E-12) {
            cout << "Warning: EDGE_NODE_INT: DIST(xi0,x1[idata[1]]) > 1.0E-12: " << DIST(xi0,x1[idata[1]]) << endl;
            //assert(0);
          }
          cout << "EDGE_NODE_INT dist looks good!" << DIST(xi0,x1[idata[1]]) << endl;
        }
        // return the node...
        FOR_I3 x[i] = x1[idata[1]][i];
        return;
      default:
        assert(0); // not handled yet
      }

    }

    string kindString() const {
      switch (kind) {
      case NODE_NODE_INT:
        return "NODE_NODE_INT";
      case EDGE_TRI_INT:
        return "EDGE_TRI_INT";
      case TRI_EDGE_INT:
        return "TRI_EDGE_INT";
      case EDGE_EDGE_INT:
        return "EDGE_EDGE_INT";
      case NODE_EDGE_INT:
        return "NODE_EDGE_INT";
      case EDGE_NODE_INT:
        return "EDGE_NODE_INT";
      default:
        assert(0);
      }
      return "UNKNOWN_INT";
    }

    void dump(const string& message) {
      switch (kind) {
      case NODE_NODE_INT:
        cout << message << " NODE_NODE_INT " << idata[0] << " " << idata[1] << endl;
        break;
      case EDGE_TRI_INT:
        cout << message << " EDGE_TRI_INT " << idata[0] << " " << ddata[0] << endl;
        break;
      case TRI_EDGE_INT:
        cout << message << " TRI_EDGE_INT " << idata[0] << " " << ddata[0] << endl;
        break;
      case EDGE_EDGE_INT:
        cout << message << " EDGE_EDGE_INT " << idata[0] << " " << ddata[0] << " " << idata[1] << " " << ddata[1] << endl;
        break;
      case NODE_EDGE_INT:
        cout << message << " NODE_EDGE_INT " << idata[0] << " " << idata[1] << " " << ddata[0] << endl;
        break;
      case EDGE_NODE_INT:
        cout << message << " EDGE_NODE_INT " << idata[0] << " " << idata[1] << " " << ddata[0] << endl;
        break;
      default:
        assert(0);
      }
    }

    inline bool operator<(const TriTriInt& rhs) const {
      return (kind < rhs.kind) || ( (kind == rhs.kind) && ((idata[0] < rhs.idata[0]) || ((idata[0] == rhs.idata[0]) && (idata[1] < rhs.idata[1]))));
    }
    inline bool operator==(const TriTriInt& rhs) const {
      return (kind == rhs.kind) && (idata[0] == rhs.idata[0]) && (idata[1] == rhs.idata[1]);
    }
    inline bool operator!=(const TriTriInt& rhs) const {
      return (kind != rhs.kind) || (idata[0] != rhs.idata[0]) || (idata[1] != rhs.idata[1]);
    }

  };

  void getTriTriX(double x[3],double *x0[3],double *x1[3],const int kind,int * idata,double * ddata);
  void writeTriTriBin(const int index,double *x0[3],double *x1[3]);
  void writeTriTriBin(const string& filename,double *x0[3],double *x1[3]);
  void writeTriTriBinDebug(const int index,double *x0[3],double *x1[3]);
  void readTriTriBin(double *x0[3],double *x1[3],const string& filename);
  void writeTriTri(const int index,double *x0[3],double *x1[3]);

  void writeEdge(FILE * fp,const double x0[3],const double x1[3]);
  void writeEdge(const string& filename,const double x0[3],const double x1[3]);
  void writeEdge(const int ied,const double x0[3],const double x1[3]);

  int calcTriTriIntersection(int idata[6],double ddata[4],double *x0[3],double *x1[3]);
  void checkTriTri();

}

#endif
