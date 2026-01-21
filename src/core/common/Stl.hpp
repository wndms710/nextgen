#ifndef STL_HPP
#define STL_HPP

#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

#include "Adt.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::min;
using std::vector;

#define STL_ASCII_FILE 1
#define PLOT3D_ASCII_FILE 2

class Stl {

 private:

  Adt<float> * fadt;
  int nf,nv;                    // facets, vertice
  float (*normal)[3];
  // XXXXXX PAUL
  //int (*v_of_f)[3];
  int (*nb_of_f)[3];
  float bbmin[3],bbmax[3];
  char filename[128];
  vector<int> facetFlag;

  // ******Mike
  int ns;                       // number of solids defined
  vector<string> solidName;     // vector of solid block names, if present
  int *sof;                     // solid-of-facet

 protected:

  // XXXXXX PAUL
  int (*v_of_f)[3];
  float (*xv)[3];

 public:

  Stl(const std::string& filename) {

    double tol = 0.0;
    init(filename.c_str(),tol);

  }

  Stl(const char * filename) {

    double tol = 0.0;
    init(filename,tol);

  }

  Stl(const std::string& filename, float tol) {

    init(filename.c_str(), tol);

  }

  Stl(const char * filename, float tol) {

    init(filename, tol);

  }

  void init(const char * filename, float tol) {

#ifdef DEBUG
    cout << "Stl(" << filename << ")" << endl;
#endif

    fadt = NULL;

    if (strlen(filename) > 128) {
      cerr << "Error: filename to long: " << filename << endl;
      throw(-1);
    }
    strcpy(this->filename,filename);

    nf = nv = ns = 0;
    normal = NULL;
    xv = NULL;
    v_of_f = NULL;
    nb_of_f = NULL;
    sof = NULL;
    solidName.clear();
    facetFlag.clear();

    int file_format = check_stl_file(filename);
    if (file_format == 1) {
         // cout << " file is binary..trying " << endl;
         initFromFileBinary(filename);
     }
     else {
         // cout << " file is ascii..OK " << endl;
         initFromFile(filename);
     }

    compress(tol);
    repair();
  }

  ~Stl() {

#ifdef DEBUG
      cout << "~Stl()" << endl;
#endif

    if ( normal != NULL ) delete[] normal;
    if ( nb_of_f != NULL ) delete[] nb_of_f;
    if ( v_of_f != NULL ) delete[] v_of_f;
    if ( xv != NULL ) delete[] xv;

    if (fadt != NULL) delete fadt;
    if (sof != NULL) delete[] sof;
  }

  // this is used for STL_GROUP
  Stl(char * *filename, int stl_count, float tol, bool asZones) : fadt (NULL) {

#ifdef DEBUG
    cout << "Stl() group instanciation" << endl;
#endif

    strcpy(this->filename,filename[0]);

    // count first
    int nf_global = 0;
    int nv_global = 0;
    int ns_global = 0;

    // loop to count facets/vertices/solids
    for (int i=0; i < stl_count; ++i) {
      nf = nv = ns = 0;
      normal = NULL;
      xv = NULL;
      v_of_f = NULL;
      nb_of_f = NULL;
      sof = NULL;
      solidName.clear();
      facetFlag.clear();

      int file_format = check_stl_file(filename[i]);
      if (file_format == 1) {
        // cout << " file is binary..trying " << endl;
        initFromFileBinary(filename[i]);
      }
      else {
        // cout << " file is ascii..OK " << endl;
        asciiCountObjects(filename[i]);
      }

      nf_global = nf_global+nf;
      nv_global = nv_global+nv;
      ns_global = ns_global+ns;
    }

    float (*normal_global)[3];
    int (*v_of_f_global)[3];
    float (*xv_global)[3];
    int (*sof_global);
    vector<string> solidName_global;

    // dimension global arrays
    normal_global = NULL;
    xv_global = NULL;
    v_of_f_global = NULL;
    sof_global = NULL;
    normal_global = new float[nf_global][3];
    xv_global = new float[nv_global][3];
    v_of_f_global = new int[nf_global][3];
    sof_global = new int[nf_global];
    solidName_global.resize(ns_global);

    if (!asZones) {
      ns_global = 1;
      solidName_global.resize(ns_global);
      solidName_global[0] = "surface";
    }

    nf_global = 0;
    nv_global = 0;
    ns_global = 0;

    // loop to populate global arrays
    for (int i=0; i < stl_count; ++i) {
      nf = nv = ns = 0;
      normal = NULL;
      xv = NULL;
      v_of_f = NULL;
      nb_of_f = NULL;
      sof = NULL;
      solidName.clear();
      facetFlag.clear();

      cout << " > adding file " << filename[i];
      if (asZones) cout << " as zone(s)";
      cout << endl;

      int file_format = check_stl_file(filename[i]);
      if (file_format == 1) {
        initFromFileBinary(filename[i]);
      }
      else {
        initFromFile(filename[i]);
      }

      // facet structures
      for (int ifa=0; ifa < nf; ifa++) {
        int ifa_index = ifa + nf_global;
        for (int i=0; i<3; ++i) {
          normal_global[ifa_index][i] = normal[ifa][i];
          v_of_f_global[ifa_index][i] = v_of_f[ifa][i] + nv_global;
        }
        if (asZones) {
          sof_global[ifa_index] = sof[ifa] + ns_global;
        }
        else {
          sof_global[ifa_index] = 0;
        }

       }
       // vertex structures
       for (int iv=0; iv < nv; iv++) {
         int iv_index = iv + nv_global;
         for (int i=0; i<3; ++i) {
           xv_global[iv_index][i] = xv[iv][i];
         }
       }
       if (asZones) {
         // solid structures
         for (int is=0; is < ns; is++) {
           int is_index = is + ns_global;
           solidName_global[is_index] = solidName[is];
         }
       }
       nf_global = nf_global + nf;
       nv_global = nv_global + nv;
       ns_global = ns_global + ns;
    }

    // move global values to the current Stl
    nf = nv = ns = 0;
    normal = NULL;
    xv = NULL;
    v_of_f = NULL;
    nb_of_f = NULL;
    sof = NULL;
    solidName.clear();
    facetFlag.clear();

    nf = nf_global;
    nv = nv_global;
    ns = ns_global;

    if (!asZones) ns = 1;

    normal = new float[nf][3];
    xv = new float[nv][3];
    v_of_f = new int[nf][3];
    sof = new int[nf];
    solidName.resize(ns);

    for  (int ifa=0; ifa < nf; ifa++) {
      for (int i=0; i<3; ++i) {
        normal[ifa][i] = normal_global[ifa][i];
        v_of_f[ifa][i] = v_of_f_global[ifa][i];
      }
      sof[ifa] = sof_global[ifa];
    }
    for (int iv=0; iv < nv; iv++) {
      for (int i=0; i<3; ++i) {
        xv[iv][i] = xv_global[iv][i];
      }
    }
    for (int is=0; is < ns; is++) {
       solidName[is] = solidName_global[is];
    }

    delete[] xv_global;
    delete[] normal_global;
    delete[] v_of_f_global;
    delete[] sof_global;

    compress(tol);
    repair();
  }

  Stl(const char * filename, int filetype, int plot3d[6]) : fadt (NULL) {

#ifdef DEBUG
    cout << "Stl(" << filename << ")" << endl;
#endif

    if (strlen(filename) > 128) {
      cerr << "Error: filename to long: " << filename << endl;
      throw(-1);
    }
    strcpy(this->filename,filename);

    nf = nv = 0;
    normal = NULL;
    xv = NULL;
    v_of_f = NULL;
    nb_of_f = NULL;
    facetFlag.clear();

    if (filetype == STL_ASCII_FILE) {
      initFromFile(filename);
    }
    else if (filetype == PLOT3D_ASCII_FILE) {
      int istart = plot3d[0];
      int iend   = plot3d[1];
      int jstart = plot3d[2];
      int jend   = plot3d[3];
      int kstart = plot3d[4];
      int kend   = plot3d[5];
      initFromPlot3DFile(filename,istart,iend,jstart,jend,kstart,kend);
    }
    else {
      cerr << " input file type " << filetype << " not supported... " << endl;
      throw(-1);
    }

    float tol = 0.0; // TODO : Introduce tolerance in STL Plot3D
    compress(tol);
    repair();
  }

  void flagFacesInsideBbox(const float x0,const float x1,const float y0,const float y1,const float z0,const float z1) {


    cout << "flagFacesInsideBbox: " << x0 << ":" << x1 << " " << y0 << ":" << y1 << " " << z0 << ":" << z1 << endl;

    int nf_flagged = 0;
    facetFlag.resize(nf);

    for (int ifa = 0; ifa < nf; ++ifa) {
      facetFlag[ifa] = 0;
      FOR_I3 {
	const int iv = v_of_f[ifa][i];
	if ((xv[iv][0] >= x0)&&(xv[iv][0] <= x1)&&(xv[iv][1] >= y0)&&(xv[iv][1] <= y1)&&(xv[iv][2] >= z0)&&(xv[iv][2] <= z1)) {
	  facetFlag[ifa] = 1;
	  ++nf_flagged;
	  break;
	}
      }
    }

    cout << " > flagged " << nf_flagged << " out of " << nf << " faces." << endl;


  }

  void flipFacetFlag() {

    assert(int(facetFlag.size()) == nf);
    for (int ifa = 0; ifa < nf; ++ifa)
      facetFlag[ifa] = 1-facetFlag[ifa];

  }

  void deleteFlaggedFaces() {

    cout << "deleteFlaggedFaces..." << endl;

    int * v_flag = new int[nv];
    for (int iv = 0; iv < nv; ++iv)
      v_flag[iv] = 0;

    for (int ifa = 0; ifa < nf; ++ifa) {
      if (facetFlag[ifa] == 0) {
	FOR_I3 {
	  const int iv = v_of_f[ifa][i];
	  v_flag[iv] = 1;
	}
      }
    }

    const int nv_old = nv;
    nv = 0;
    for (int iv = 0; iv < nv_old; ++iv) {
      if (v_flag[iv] == 1) {
	FOR_I3 xv[nv][i] = xv[iv][i]; // copy down xv...
	v_flag[iv] = nv++;
      }
      else {
	v_flag[iv] = -1;
      }
    }

    cout << " > nv: " << nv << endl;

    const int nf_old = nf;
    nf = 0;
    for (int ifa = 0; ifa < nf_old; ++ifa) {
      if (facetFlag[ifa] == 0) {
	FOR_I3 {
	  const int iv = v_of_f[ifa][i];
	  v_of_f[nf][i] = v_flag[iv];
	}
	++nf;
      }
    }

    cout << " > nf: " << nf << endl;

    delete[] v_flag;

    facetFlag.resize(nf);

  }

  int getNf() const { return(nf); }
  int getNv() const { return(nv); }
  int getNs() const { return(ns); }

  string getSolidName(const int idx) const {
    assert((idx < ns) && (idx >= 0));
    return(solidName[idx]);
  }

  int getSOfF(const int ifacet) const {
    assert((ifacet < nf) && (ifacet >= 0));
    return(sof[ifacet]);
  }

  float getXv(int iv,int id) const {
#ifdef DEBUG
    if ((iv < 0)||(iv >= nv)||(id < 0)||(id >= 3)) {
      cerr << "Error: out of range in getXv." << endl;
      throw(-1);
    }
#endif
    return(xv[iv][id]);
  }

  int getVOfF(int ifa,int vof) const {
#ifdef DEBUG
    if ((ifa < 0)||(ifa >= nf)||(vof < 0)||(vof >= 3)) {
      cerr << "Error: out of range in getVOfF." << endl;
      throw(-1);
    }
#endif
    return(v_of_f[ifa][vof]);
  }

  void translate(float dx,float dy,float dz) {
    int iv;
    for (iv=0;iv<nv;iv++) {
      xv[iv][0] += dx;
      xv[iv][1] += dy;
      xv[iv][2] += dz;
    }
  }

  void translateX(float dx) {
    int iv;
    for (iv=0;iv<nv;iv++) {
      xv[iv][0] += dx;
    }
  }

  void translateY(float dy) {
    int iv;
    for (iv=0;iv<nv;iv++) {
      xv[iv][1] += dy;
    }
  }

  void translateZ(float dz) {
    int iv;
    for (iv=0;iv<nv;iv++) {
      xv[iv][2] += dz;
    }
  }

  void scale(float factor) {
    int iv;
    for (iv=0;iv<nv;iv++) {
      int i;
      for (i=0;i<3;i++) {
        xv[iv][i] *= factor;
      }
    }
  }

  void scaleX(float factor) {
    int iv;
    for (iv=0;iv<nv;iv++) {
      xv[iv][0] *= factor;
    }
  }

  void scaleY(float factor) {
    int iv;
    for (iv=0;iv<nv;iv++) {
      xv[iv][1] *= factor;
    }
  }

  void scaleZ(float factor) {
    int iv;
    for (iv=0;iv<nv;iv++) {
      xv[iv][2] *= factor;
    }
  }

  void rotateX(float degrees) {
    int iv;
    float rad = degrees*M_PI/180.0;
    float cos_t = cos(rad);
    float sin_t = sin(rad);
    for (iv=0;iv<nv;iv++) {
      float y = xv[iv][1]*cos_t - xv[iv][2]*sin_t;
      float z = xv[iv][2]*cos_t + xv[iv][1]*sin_t;
      xv[iv][1]  = y;
      xv[iv][2]  = z;
    }
  }

  void rotateY(float degrees) {
    int iv;
    float rad = degrees*M_PI/180.0;
    float cos_t = cos(rad);
    float sin_t = sin(rad);
    for (iv=0;iv<nv;iv++) {
      float z = xv[iv][2]*cos_t - xv[iv][0]*sin_t;
      float x = xv[iv][0]*cos_t + xv[iv][2]*sin_t;
      xv[iv][2]  = z;
      xv[iv][0]  = x;
    }
  }

  void rotateZ(float degrees) {
    int iv;
    float rad = degrees*M_PI/180.0;
    float cos_t = cos(rad);
    float sin_t = sin(rad);
    for (iv=0;iv<nv;iv++) {
      float x = xv[iv][0]*cos_t - xv[iv][1]*sin_t;
      float y = xv[iv][1]*cos_t + xv[iv][0]*sin_t;
      xv[iv][0]  = x;
      xv[iv][1]  = y;
    }
  }

  void mirrorX() {
    int iv;
    for (iv=0;iv<nv;iv++) {
      xv[iv][0] = -xv[iv][0];
    }
  }

  void mirrorY() {
    int iv;
    for (iv=0;iv<nv;iv++) {
      xv[iv][1] = -xv[iv][1];
    }
  }

  void mirrorZ() {
    int iv;
    for (iv=0;iv<nv;iv++) {
      xv[iv][2] = -xv[iv][2];
    }
  }

  void calcBbox() {

    if (nv == 0) {
      cerr << "Error: zero vertices." << endl;
      throw(-1);
    }

    int i;
    for (i=0;i<3;i++) {
      bbmin[i] = xv[0][i];
      bbmax[i] = xv[0][i];
    }

    int iv;
    for (iv=1;iv<nv;iv++) {
      for (i=0;i<3;i++) {
        bbmin[i] = min(xv[iv][i],bbmin[i]);
        bbmax[i] = max(xv[iv][i],bbmax[i]);
      }
    }

#ifdef DEBUG
    cout << "stl bounding box: " << bbmin[0] << " " << bbmax[0] << " " <<
      bbmin[1] << " " << bbmax[1] << " " <<
      bbmin[2] << " " << bbmax[2] << endl;
#endif

  }

  void calcBbox(float * bbmin,float * bbmax) {

    if (nv == 0) {
      cerr << "Error: zero vertices." << endl;
      throw(-1);
    }

    int i;
    for (i=0;i<3;i++) {
      bbmin[i] = xv[0][i];
      bbmax[i] = xv[0][i];
    }

    int iv;
    for (iv=1;iv<nv;iv++) {
      for (i=0;i<3;i++) {
        bbmin[i] = min(xv[iv][i],bbmin[i]);
        bbmax[i] = max(xv[iv][i],bbmax[i]);
      }
    }

#ifdef DEBUG
    cout << "stl bounding box: " << bbmin[0] << " " << bbmax[0] << " " <<
      bbmin[1] << " " << bbmax[1] << " " <<
      bbmin[2] << " " << bbmax[2] << endl;
#endif

  }

  void calcFacetBoundingBox(double * bb_min,double * bb_max,const int ifacet) const {
    float bb_min_float[3],bb_max_float[3];
    calcFacetBoundingBox(bb_min_float,bb_max_float,ifacet);
    bb_min[0] = (double)bb_min_float[0];
    bb_min[1] = (double)bb_min_float[1];
    bb_min[2] = (double)bb_min_float[2];
    bb_max[0] = (double)bb_max_float[0];
    bb_max[1] = (double)bb_max_float[1];
    bb_max[2] = (double)bb_max_float[2];
  }

  void calcFacetBoundingBox(float * bb_min,float * bb_max,const int ifacet) const {
    // set the bbox from the first vertex...
    int iv = v_of_f[ifacet][0];
    for (int i = 0; i < 3; ++i) {
      bb_min[i] = xv[iv][i];
      bb_max[i] = xv[iv][i];
    }
    // expand bbox based on other 2 vertices...
    for (int j = 1; j <= 2; ++j) {
      iv = v_of_f[ifacet][j];
      for (int i = 0; i < 3; ++i) {
	bb_min[i] = min(bb_min[i],xv[iv][i]);
	bb_max[i] = max(bb_max[i],xv[iv][i]);
      }
    }
  }

  void calcFacetNormal(double * normal,int ifacet) const {
    float normal_float[3];
    calcFacetNormal(normal_float,ifacet);
    normal[0] = (double)normal_float[0];
    normal[1] = (double)normal_float[1];
    normal[2] = (double)normal_float[2];
  }

  void calcFacetNormal(float * normal,int ifacet) const {
    int iv0 = v_of_f[ifacet][0];
    int iv1 = v_of_f[ifacet][1];
    int iv2 = v_of_f[ifacet][2];
    double v01[3],v02[3];
    int i;
    for (i=0;i<3;i++) {
      v01[i] = double(xv[iv1][i]) - double(xv[iv0][i]);
      v02[i] = double(xv[iv2][i]) - double(xv[iv0][i]);
    }
    normal[0] = 0.5*float(v01[1]*v02[2] - v01[2]*v02[1]);
    normal[1] = 0.5*float(v01[2]*v02[0] - v01[0]*v02[2]);
    normal[2] = 0.5*float(v01[0]*v02[1] - v01[1]*v02[0]);
  }

  void calcFacetUnitNormal(double * normal,int ifacet) const {
    float normal_float[3];
    calcFacetUnitNormal(normal_float,ifacet);
    normal[0] = (double)normal_float[0];
    normal[1] = (double)normal_float[1];
    normal[2] = (double)normal_float[2];
  }

  void calcFacetUnitNormal(float * normal,int ifacet) const {
    calcFacetNormal(normal,ifacet);
    float mag = sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );
    if (mag == 0.0) {
      cout << "Warning: normal has zero magnitude." << endl;
      normal[0] = 1.0;
      normal[1] = 0.0;
      normal[2] = 0.0;
    }
    else {
      normal[0] /= mag;
      normal[1] /= mag;
      normal[2] /= mag;
    }
  }

  void split(const double degrees) {

    cout << "Split: degrees = " << degrees << endl;
    double cos_theta = cos(degrees*M_PI/180.0);
    cout << "cos_theta: " << cos_theta << endl;

    assert(nb_of_f != NULL);
    facetFlag.resize(nf);

    for (int ifa = 0; ifa < nf; ++ifa)
      facetFlag[ifa] = 0;

    int * facetQueue = new int[nf];
    int queueSize = 0;

    int ifa_f = 0;
    int done = 0;
    int iter = 0;
    while (!done) {

      // find the first zero facet from the old flagged set.
      // This is the starting point...
      while (facetFlag[ifa_f] != 0) {
	++ifa_f;
	if (ifa_f == nf)
	  break;
      }
      if (ifa_f == nf)
	break;

      // zero all facets...
      for (int ifa = 0; ifa < nf; ++ifa)
	facetFlag[ifa] = 0;

      // pop the first onto the queue...
      facetQueue[queueSize++] = ifa_f;
      facetFlag[ifa_f] = 1;

      // and process...
      while (queueSize != 0) {

	// pop the top face off the queue...
	int ifa = facetQueue[--queueSize];
	assert(facetFlag[ifa] == 1);

	double n1[3];
	calcFacetUnitNormal(n1,ifa);

	for (int i = 0; i < 3; ++i) {
	  int ifa_nbr = nb_of_f[ifa][i];
	  if ((ifa_nbr >= 0)&&(facetFlag[ifa_nbr] == 0)) {
	    double n2[3];
	    calcFacetUnitNormal(n2,ifa_nbr);
	    if (fabs(n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2]) > cos_theta) {
	      // pop nbr on the queue...
	      facetQueue[queueSize++] = ifa_nbr;
	      facetFlag[ifa_nbr] = 1;
	    }
	  }
	}

      }

      // we should have a bunch of facets flagged...
      char newfilename[128];
      sprintf(newfilename,"%s.%02d.stl",filename,iter++);
      writeFlaggedFacets(newfilename);

    }


  }

  void write(const char * filename) {

    cout << " > write: " << filename << endl;

    FILE * fp;
    if ( (fp=fopen(filename,"w"))==NULL ) {
      cerr << "Error: cannot open file " << filename << endl;
      throw(-1);
    }

    for (int is=0; is < ns; ++is) {
      fprintf(fp,"solid %s\n",solidName[is].c_str());

      for (int ifa = 0; ifa < nf; ++ifa) {
        if (sof[ifa] == is) {
          fprintf(fp," facet normal 0 0 0\n");
          fprintf(fp,"  outer loop\n");
          for (int i = 0; i < 3; ++i) {
            fprintf(fp,"   vertex");
            int iv = v_of_f[ifa][i];
            for (int j = 0; j < 3; ++j) fprintf(fp," %f",xv[iv][j]);
            fprintf(fp,"\n");
          }
          fprintf(fp,"  endloop\n");
          fprintf(fp," endfacet\n");
        }
      }

      fprintf(fp,"endsolid %s\n", solidName[is].c_str());
    }
    fclose(fp);

  }

  void writeFlaggedFacets(const char * filename) {

    cout << " > writeFlaggedFacets: " << filename << endl;

    FILE * fp;
    if ( (fp=fopen(filename,"w"))==NULL ) {
      cerr << "Error: cannot open file " << filename << endl;
      throw(-1);
    }

    fprintf(fp,"solid\n");

    for (int ifa = 0; ifa < nf; ++ifa) {
      if (facetFlag[ifa] == 1) {
	fprintf(fp,"facet normal 1 1 1\n");
	fprintf(fp,"outer loop\n");
	for (int i = 0; i < 3; ++i) {
	  fprintf(fp,"vertex");
	  int iv = v_of_f[ifa][i];
	  for (int j = 0; j < 3; ++j)
	    fprintf(fp," %f",xv[iv][j]);
	  fprintf(fp,"\n");
	}
	fprintf(fp,"endloop\n");
	fprintf(fp,"endfacet\n");
      }
    }

    fprintf(fp,"endsolid\n");

    fclose(fp);

  }

  void initFromFile(const char * filename);
  void initFromFileBinary(const char * filename);
  void initFromPlot3DFile(const char * filename, int istart, int iend, int jstart, int jend, int kstart, int kend);
  void compress(float tol);
  void repair();
  void writeTecplot(const char * filename);
  void writeTecplotZone(FILE * fp,char * zonename);

  // stuff that uses facetFlag...
  void zeroFacetFlag() {
    if (int(facetFlag.size()) < nf)
      facetFlag.resize(nf);
    int ifa;
    for (ifa = 0; ifa < nf; ifa++)
      facetFlag[ifa] = 0;
  }

  void setFacetFlag(int ifa,int flag) {
    assert((ifa >= 0)&&(ifa <= nf)&&(int(facetFlag.size())==nf));
    facetFlag[ifa] = flag;
  }

  void writeFlaggedFacetsTecplot(const char * filename,int flag);
  void writeFlaggedFacetsTecplotZone(FILE * fp,char * zonename,int flag);

  //---------- stuff for binary STLs

#if !defined(SEEK_SET)
#define SEEK_SET 0
#endif
#define HEADER_SIZE 84
#define SIZEOF_STL_FACET 50
#define STL_MIN_FILE_SIZE 284
#define LABEL_SIZE 80

  int stl_get_number_of_facet(const char * filename);
  int check_stl_file(const char * filename);
  void asciiCountObjects(FILE * fp);
  void asciiCountObjects(const char * filename);
  void asciiReadObjects(FILE * fp);
  void asciiReadObjects(const char * filename);
  int stl_get_little_int(FILE *fp);
  float stl_get_little_float(FILE *fp);

  // point-search stuff...

  void setClosestPointAndNormal(float x[],float n[],const float xp[]) {

    // make sure the facet-adt still exists...
    assert(fadt != NULL);

    std::vector<int> bbVec;
    fadt->buildListForClosestPoint(bbVec,xp);
    assert(bbVec.size() > 0);

    float d2_min = 1.0E+20;
    int ifacet_min = -1;
    for (unsigned int ii = 0; ii < bbVec.size(); ++ii) {
      const int ifacet = bbVec[ii];
      float this_x[3];
      const float this_d2 = calcSingleFacetClosestPoint(this_x,xp,xv[v_of_f[ifacet][0]],xv[v_of_f[ifacet][1]],xv[v_of_f[ifacet][2]]);
      if (this_d2 < d2_min) {
	d2_min = this_d2;
	FOR_I3 x[i] = this_x[i];
	ifacet_min = ifacet;
      }
    }

    assert(ifacet_min != -1);
    calcFacetNormal(n,ifacet_min);

    /*
    // CHECK...
    // loop through every tri and do the calc again...
    // Note that the facet may NOT match, because 2 facets that return
    // the same distance may appear in a different order...
    {
      float d2_check = 1.0E+20;
      float x_check[3];
      int ifacet_check = -1;
      for (int ifacet = 0; ifacet < nf; ++ifacet) {
	float this_x[3];
	const float this_d2 = calcSingleFacetClosestPoint(this_x,xp,xv[v_of_f[ifacet][0]],xv[v_of_f[ifacet][1]],xv[v_of_f[ifacet][2]]);
	if (this_d2 < d2_check) {
	  d2_check = this_d2;
	  FOR_I3 x_check[i] = this_x[i];
	  ifacet_check = ifacet;
	}
      }
      // note that the facet may not match...
      assert(d2_min == d2_check);
      //assert(x_check[0] == x[0]);
      //assert(x_check[1] == x[1]);
      //assert(x_check[2] == x[2]);
    }
    */

  }

  float calcSingleFacetClosestPoint(float xc[],const float xp[],const float v0[],const float v1[],const float v2[]) const {

    //Compute Facet Triangle Edges: T(s,t) = v0 + s*edge0 + t*edge1
    double edge0[3], edge1[3];
    FOR_I3{
      edge0[i] = double(v1[i]) - double(v0[i]);
      edge1[i] = double(v2[i]) - double(v0[i]);
    }

    //Compute coefficients of distance squared function, Q(s,t)
    double diff_v0xp[3];
    double qa, qb, qc, qd, qe; //, qf;

    FOR_I3 diff_v0xp[i] = double(v0[i]) - double(xp[i]);
    qa = DOT_PRODUCT(edge0,edge0);
    qb = DOT_PRODUCT(edge0,edge1);
    qc = DOT_PRODUCT(edge1,edge1);
    qd = DOT_PRODUCT(edge0,diff_v0xp);
    qe = DOT_PRODUCT(edge1,diff_v0xp);
    //qf = DOT_PRODUCT(diff_v0xp,diff_v0xp);

    //Compute coordinates Qmin = Q(sbar/det,tbar/det)
    double det, sbar, tbar;
    det = qa*qc - qb*qb;
    if (det <= 0.0) {
      // edges of facet are colinear...
      cout << "det: " << det << endl;
      cout << "edge0: " << edge0[0] << " " << edge0[1] << " " << edge0[2] << endl;
      cout << "edge1: " << edge1[0] << " " << edge1[1] << " " << edge1[2] << endl;
      throw(0);
    }
    sbar = qb*qe - qc*qd;
    tbar = qb*qd - qa*qe;

    /***********************************
    * Identify region containing xp
    *        t
    *      \2|
    *       \|
    *        \
    *        |\
    *        | \
    *      3 | 0\  1
    *     ---|---\-----s
    *        |    \
    *      4 | 5   \  6
    ************************************/
    int region = -1;
    if (sbar + tbar <= det){
      if (sbar<0){
        if (tbar<0) region = 4;
        else region = 3;
      }
      else if (tbar<0) region = 5;
      else region = 0;
    }
    else{
      if (sbar<0) region = 2;
      else if (tbar<0) region = 6;
      else region = 1;
    }
    assert(region>-1);
    //Facet coordinates that give min distance to xp
    //will be assigned to sbar, tbar
    double invDet, numer, denom, tmp0, tmp1;
    switch (region) {
      case 0:
        invDet = 1/det;
        sbar *= invDet;
        tbar *= invDet;
        break;
      case 1:
        numer = (qc+qe-qb-qd);
        if (numer<=0) sbar = 0;
        else {
          denom = qa - 2*qb + qc; //positive
          sbar = ( numer >= denom ? 1 : numer/denom );
        }
        tbar = 1 - sbar;
        break;
      case 2:
        tmp0 = qb + qd;
        tmp1 = qc + qe;
        if (tmp1 > tmp0){ //minimum on edge s+t = 1
          numer = tmp1 - tmp0;
          denom = qa - 2*qb + qc;
          sbar = ( numer >= denom ? 1 : numer/denom);
          tbar = 1-sbar;
        }
        else { //minimum on edge s = 0
          sbar = 0;
          tbar = ( tmp1 <= 0 ? 1 : ( qe >= 0 ? 0 : -qe/qc ) );
        }
        break;
       case 3:
        sbar = 0;
        tbar = ( qe >= 0 ? 0 : ( -qe >= qc ? 1 : -qe/qc ) );
        break;
      case 4:
        if (qe < 0){ //minimum on edge s = 0
          sbar = 0;
          tbar = ( qc <= -qe ? 1 : -qe/qc );
        }
        else { //minimum on edge t = 0
          tbar = 0;
          sbar = ( qa <= -qd ? 1 : ( qd >= 0 ? 0 : -qd/qa ) );
        }
        break;
      case 5:
        tbar = 0;
        sbar = ( qd >= 0 ? 0 : ( -qd >= qa ? 1 : -qd/qa ) );
        break;
      case 6:
        tmp0 = qa + qd;
        tmp1 = qb + qe;
        if (tmp0 > tmp1){//minimum on edge s+t = 1
          numer = tmp0 - tmp1;
          denom = qa - 2*qb + qc;
          tbar = ( numer >= denom ? 1 : numer/denom );
          sbar = 1-tbar;
        }
        else{ //minimum on edge t = 0
          tbar = 0;
          sbar = ( tmp0 <= 0 ? 1 : ( qd >= 0 ? 0 : -qd/qa ) );
        }
        break;
    }

    //Return point coordinates in physical space
    FOR_I3 xc[i] = sbar*edge0[i] + tbar*edge1[i] + v0[i];

    //return distance squared
    return float( pow(double(xp[0])-double(xc[0]),2) +
		  pow(double(xp[1])-double(xc[1]),2) +
		  pow(double(xp[2])-double(xc[2]),2) );

  }

};

#endif
