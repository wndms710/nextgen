#ifndef FLUENT_BP_READER_HPP
#define FLUENT_BP_READER_HPP

#include "Common.hpp"
#include "CTI.hpp"
using namespace CTI;
#include "Adt.hpp"
#include "parsers/fluent/FluentParser.hpp"

enum ProfileType {
  POINT_3D,
  LINE_R,
  LINE_WD
};

class FluentBPReader : FluentParser {

private:

  int n;
  string name;
  string filename;
  string data_type;
  vector<string> varNameVec;
  vector<double*> varDataVec;

  // for now just simple nearest nbr
  int nno;
  int* index_no;
  double* wgt_no;  // for radial profiles, stores wgt between indices

  ProfileType type;
  bool use_mean;

public:

  FluentBPReader() {
    n = 0;
    use_mean = false;
    nno = 0;
    index_no = NULL;
    wgt_no = NULL;
  }

  FluentBPReader(const string& _filename): FluentParser(_filename),filename(_filename) {
    n = 0;
    use_mean = false;
    nno = 0;
    index_no = NULL;
    wgt_no = NULL;

    readData();

    if (data_type == "point") {
      type = POINT_3D;
      COUT1(" > zero-th order interpolation (nearest neighbor) being used for \"point\" format");
    }//"point"
    else if (data_type == "radial") {
      type = LINE_R;
      COUT1(" > linear interpolation being used for \"radial\" format");
    }//"radial"
    else if (data_type == "wall-distance") {
      type = LINE_WD;
      COUT1(" > linear interpolation being used for \"wall-distance\" format");
    }//"wall-distance"

  }//FluentBPReader()

  ~FluentBPReader() {
    for (int i=0, end=varDataVec.size(); i<end; ++i) {
      DELETE(varDataVec[i]);
    }
    DELETE(index_no);
    DELETE(wgt_no);
  }//~FluentBPReader()

  void init(const string& _filename) {
    filename = _filename;
    FluentParser::init(filename);

    readData();

    if (data_type == "point") {
      type = POINT_3D;
      COUT1(" > zero-th order interpolation (nearest neighbor) being used for \"point\" format");
    }//"point"
    else if (data_type == "radial") {
      type = LINE_R;
      COUT1(" > linear interpolation being used for \"radial\" format");
    }//"radial"
    else if (data_type == "wall-distance") {
      type = LINE_WD;
      COUT1(" > linear interpolation being used for \"wall-distance\" format");
    }//"wall-distance"
  }

  bool isInitialized() {return (!filename.empty());}

  // =============================================
  // check if a variable exist in the profile list
  // =============================================
  void setUseMean(const bool _use_mean) {
    use_mean = _use_mean;
    if (use_mean) COUT1(" > using mean-variables to set boundary state");
  }

  // =============================================
  // check if a variable exist in the profile list
  // =============================================
  bool checkVar(const string& name) {
    string var_name;
    if ((name == "x") || (name == "y") || (name == "z") || (name == "radius") || (name == "wall-distance") || !use_mean) {
      var_name = name;
    }
    else {
      assert(use_mean);
      var_name = "mean-"+name;
    }

    for (int i=0,end=varNameVec.size(); i<end; ++i) {
      if ( varNameVec[i] == var_name ) return(true);
    }
    return(false);

  }//checkVar()

  void ensureVar(const string& name) {
    string var_name;
    if ((name == "x") || (name == "y") || (name == "z") || (name == "radius") || (name == "wall-distance") || !use_mean) {
      var_name = name;
    }
    else {
      assert(use_mean);
      var_name = "mean-"+name;
    }

    for (int i=0,end=varNameVec.size(); i<end; ++i) {
      if ( varNameVec[i] == var_name ) return;
    }

    CERR("could not find variable " << var_name << " in profile " << filename);
  }

  // ==============================================
  // set points
  // ==============================================
  void setPointsAdt(const double (*x)[3],const double x0[3],const int nno) {
    this->nno   = nno;
    assert( index_no == NULL );
    index_no = new int[nno];

    if (type == POINT_3D) {

      // adt-based search
      assert(checkVar("x") && checkVar("y") && checkVar("z"));

      double* profile_x = getVarPtr("x");
      double* profile_y = getVarPtr("y");
      double* profile_z = getVarPtr("z");

      COUT1("    > building point cloud Adt...");

      double (*prof_xyz)[3] = new double [n][3];
      for (int ino=0; ino<n; ++ino){
        prof_xyz[ino][0] = profile_x[ino];
        prof_xyz[ino][1] = profile_y[ino];
        prof_xyz[ino][2] = profile_z[ino];
      }
      Adt<double> * adt = new Adt<double>(n,prof_xyz,prof_xyz);

      DELETE(prof_xyz);
      vector<int> inoVec;

      COUT1("    > finding closest profile point...");
      // find closest profile node for each point in x
      for (int ix=0; ix<nno; ++ix){
        inoVec.clear();
        adt->buildListForClosestPoint(inoVec,x[ix]);

        if  (inoVec.empty()) {
          assert(0);
        }
        else {
          assert(inoVec.size() == 1);
          index_no[ix] = inoVec[0];
        }
      }//(ino < nno)
      delete adt;
    }//(POINT_3D)
    else {
      assert(0);  // calling the wrong constructor for points!
    }
  }//setPoints()

  void setPoints(const double (*x)[3], const double x0[3],const int nno) {
    this->nno   = nno;
    assert( index_no == NULL );
    index_no = new int[nno];

    if (type == POINT_3D) {
      // for now simple nearest neighbor
      // with n^2 search algorithm
      assert(checkVar("x") && checkVar("y") && checkVar("z"));

      double* profile_x = getVarPtr("x");
      double* profile_y = getVarPtr("y");
      double* profile_z = getVarPtr("z");

      for (int ino=0; ino<nno; ++ino){
        double dist3[3];
        dist3[0] = x[ino][0] - profile_x[0];
        dist3[1] = x[ino][1] - profile_y[0];
        dist3[2] = x[ino][2] - profile_z[0];
        double dist_last = dist3[0]*dist3[0] + dist3[1]*dist3[1] + dist3[2]*dist3[2];
        int this_index = 0;
        for (int j = 1; j<n; ++j) {
          dist3[0] = x[ino][0] - profile_x[j];
          dist3[1] = x[ino][1] - profile_y[j];
          dist3[2] = x[ino][2] - profile_z[j];
          double dist = dist3[0]*dist3[0] + dist3[1]*dist3[1] + dist3[2]*dist3[2];
          if ( dist < dist_last ) {
            dist_last = dist;
            this_index = j;
          }
        }//(j < n)
        index_no[ino] = this_index;
      }//(ino < nno)
    }//(POINT_3D)
    else if ((type == LINE_R) || (type == LINE_WD)) {
      assert( wgt_no == NULL );
      wgt_no = new double[nno];

      double * profile_r;
      if (type == LINE_R) {
        assert(checkVar("radius"));
        profile_r = getVarPtr("radius");
      }
      else {
        assert(checkVar("wall-distance"));
        profile_r = getVarPtr("wall-distance");
      }

      for (int ino=0; ino<nno; ++ino){
        const double dx[3] = DIFF(x[ino],x0);
        const double my_r = MAG(dx);

        int l_index = findMinIndex(my_r,profile_r);
        if (l_index == -1) {
          wgt_no[ino] = 0.0;
          l_index = 0;
        }
        else if (l_index == n-1) {
          wgt_no[ino] = 0.0;
        }
        else {
          wgt_no[ino] = (my_r-profile_r[l_index])/(profile_r[l_index+1] - profile_r[l_index]);
        }
        assert((l_index >= 0) && (l_index < n));
        assert((wgt_no[ino] >= 0.0) && (wgt_no[ino] <= 1.0));
        index_no[ino] = l_index;  // always store left index
      }
    }//(LINE_R) || (LINE_WD))
  }//setPoints()

  void setPoints(const double* my_r,const int nno) {
    assert((type == LINE_R) || (type == LINE_WD));  // assumes user has computed radius value already, and setting points via this

    this->nno   = nno;
    assert( index_no == NULL );
    index_no = new int[nno];

    assert( wgt_no == NULL );
    wgt_no = new double[nno];

    double* profile_r;
    if (type == LINE_R) {
      assert(checkVar("radius"));
      profile_r = getVarPtr("radius");
    }
    else {
      assert(checkVar("wall-distance"));
      profile_r = getVarPtr("wall-distance");
    }

    for (int ino=0; ino<nno; ++ino){
      int l_index = findMinIndex(my_r[ino],profile_r);
      if (l_index == -1) {
        wgt_no[ino] = 0.0;
        l_index = 0;
      }
      else if (l_index == n-1) {
        wgt_no[ino] = 0.0;
      }
      else {
        wgt_no[ino] = (my_r[ino]-profile_r[l_index])/(profile_r[l_index+1] - profile_r[l_index]);
      }
      assert((l_index >= 0) && (l_index < n));
      assert((wgt_no[ino] >= 0.0) && (wgt_no[ino] <= 1.0));
      index_no[ino] = l_index;  // always store left index
    }//(ino < nno)
  }//setPoints()

  int findMinIndex(const double r,double * profile_r) {
      // assumes profile_r is sorted
      // return position of lower bounding element
    int left = 0;
    int right = n-1;
    int middle;

    while (left <= right) {
      middle = (int)floor(0.5*(left+right));
      if (profile_r[middle] <= r) left = middle+1;
      else right = middle-1;
    }

    return --left;
  }//findMinIndex()

  // ========
  // get data
  // ========
  double getData(const int ino, const string& name) {
    if ((type == LINE_R) || (type == LINE_WD)) {
      assert((ino>=0)&&(ino < nno));
      double* data = getVarPtr(name);
      if (ino == (nno-1)) {
        return data[index_no[ino]];
      }
      else {
        return (data[index_no[ino]] + wgt_no[ino]*(data[index_no[ino]+1] - data[index_no[ino]]));
      }
    }
    else {
      assert(type == POINT_3D);
      assert((ino>=0)&&(ino < nno));
      double* data = getVarPtr(name);
      return(data[index_no[ino]]);
    }
  }//getData()

  ProfileType getType() const {
    return type;
  }//getType()

  // ===================================
  // write the profile to a tecplot file
  // ===================================
  void dumpTecplot(const string& filename) {
    if (mpi_rank == 0) {
      FILE * fp = fopen(filename.c_str(),"w");
      fprintf(fp,"title = \"Fluent Profile\"\n");
      fprintf(fp,"variables = ");
      for (int i=0; i<varNameVec.size(); ++i) fprintf(fp,"\"%s\"\n",varNameVec[i].c_str());
      fprintf(fp,"zone i = %d, DATAPACKING=POINT\n",n);
      for (int j=0; j<n; ++j) {
        for (int i=0; i<varNameVec.size(); ++i) fprintf(fp," %lf",varDataVec[i][j]);
        fprintf(fp,"\n");
      }
      fclose(fp);
    }
  }//dumpTecplot()

  void printVarsInFile() const {
    if (varNameVec.empty()) return;

    COUT1(" > variables found in profile:");
    for (int i=0; i<varNameVec.size(); ++i) {
      COUT1("    - " << varNameVec[i]);
    }
  }//printVarsInFile()

private:
  // ============================
  // get a variable from the list
  // ============================
  double* getVarPtr(const string& name){
    string var_name;
    if ((name == "x") || (name == "y") || (name == "z") || (name == "radius") || (name == "wall-distance") || (!use_mean)) {
      var_name = name;
    }
    else {
      assert(use_mean);
      var_name = "mean-"+name;
    }

    for (int i=0; i<varNameVec.size(); ++i) {
      if (varNameVec[i] == var_name) return(varDataVec[i]);
    }

    // we should not be here
    CERR("could not find variable " << var_name << " from profile " << filename);
    return NULL;
  }//getVarPtr()


  // =============================
  // read the profile on rank zero
  // =============================
  void readData() {
    // opening of file handled in FluentParser, so now just parse

    // ==========
    // format: ((title type n-elements)
    //         (var-a <a0> <a1> <a2> ....)
    //         (var-b <b0> <b1> <b2> ....)
    //         )
    // ==========
    {
      advanceToLevel(2);
      name = getNextTokenAsString(2);
      data_type = getNextTokenAsString(2);
      n = getNextTokenAsInt(2);
      COUT1(" > profile name: " << name << ", data-type: " << data_type << ", points: " << n);

      if ((data_type != "point") && (data_type != "radial") && (data_type != "wall-distance")) {
        WUI(ERR,"problematic profile format \"" << data_type << "\" found; currently only 3D \"point\" and \"radial\" are supported. Please reformat your profile accordingly");
      }
      else if (data_type == "radial") {
        CWARN("assumes data in file is sorted by radius");
      }
    }

    // ==========
    // variables
    // ==========
    while (advanceToLevel(2)) {
      const string var_name = getNextTokenAsString(2);

      COUT1("    > reading variable " << var_name << "...");
      varNameVec.push_back(var_name);

      double * tmpdbl = new double[n];

      readVar(tmpdbl);
      varDataVec.push_back(tmpdbl);
    }
  }

  void readVar(double * var) {
    for (int i = 0; i < n; ++i) {
      var[i] = getNextTokenAsDouble(2);
    }
  }
};//class FluentBPReader
#endif
