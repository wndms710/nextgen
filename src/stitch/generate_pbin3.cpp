
#include <iostream>
#include <cstdio>
#include <assert.h>
#include <climits>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

typedef long long int int8;

#define POINTS_IO_VERSION 2
#define POINTS_IO_MAGIC_NUMBER 1235813
#define HUGE_INT (1<<24)

void writePointsToFile(const string& filename,const vector<double>& xpVec,const vector<double>& deltaVec) {

  FILE* fp = fopen(filename.c_str(), "wb");
  assert(fp);
  
  assert(xpVec.size() == deltaVec.size()*3);
  const int np = deltaVec.size(); 

  int8 ibuf[4] = {POINTS_IO_MAGIC_NUMBER, POINTS_IO_VERSION, np, 1};

  fwrite(ibuf, sizeof(int8),4,fp);
  fwrite(&xpVec.front(),sizeof(double),np*3,fp);
  fwrite(&deltaVec.front(),sizeof(double),np,fp);
  fclose(fp);

}

int main(const int argc, const char* argv[]) { 
  
  vector<double> xpVec;
  vector<double> deltaVec;

  // set the coarse grid size. The 2D grid will fine from -1:0 and crs from 0:1... 
  //const int n = 32;
  const int n = atoi(argv[1]);

  // 2*n+1 points spacing 1 in the interval -1.0..0.0
  for (int i = 0; i <= 6*n; ++i) {
    double x = -3.0 + (double(i)/double(6*n))*3.0+0.00001; // apply x-shift for x-periodicity
    cout << "fin: " << x << endl;
    xpVec.push_back(x); xpVec.push_back(0.0); xpVec.push_back(0.0);
    deltaVec.push_back(1.1); // TODO: too large!
  }

  // or n+1 crs...
  /*
  for (int i = 0; i <= n; ++i) {
    double x = -1.0 + (double(i)/double(n))*1.0+0.00001; // apply x-shift for x-periodicity
    cout << x << endl;
    xpVec.push_back(x); xpVec.push_back(0.0); xpVec.push_back(0.0);
    deltaVec.push_back(1.1); // TODO: too large!
  }
  */

  // n-1 points spacing 1 in the interval 0.0..1.0
  for (int i = 1; i < 3*n; ++i) {
    double x = 0.0 + (double(i)/double(3*n))*3.0+0.00001;
    cout << "crs: " << x << endl;
    xpVec.push_back(x); xpVec.push_back(0.0); xpVec.push_back(0.0);
    deltaVec.push_back(2.2);
  }
 
  cout << "about to smooth" << endl;
  getchar();
  
  // apply the smoothing here...
  int np = deltaVec.size();
  double *dx = new double[np];
  for (int iter = 0; iter < 10; ++iter) {
    for (int ip = 0; ip < np; ++ip) {
      double xm1;
      if (ip > 0) {
	xm1 = xpVec[(ip-1)*3];
      }
      else {
	xm1 = xpVec[(np-1)*3]-6.0;
      }
      double xp1;
      if (ip < np-1) {
	xp1 = xpVec[(ip+1)*3];
      }
      else {
	xp1 = xpVec[0]+6.0;
      }
      dx[ip] = 0.25*xm1 - 0.5*xpVec[ip*3] + 0.25*xp1;
      cout << "got dx: " << dx[ip] << endl;
    }
    for (int ip = 0; ip < np; ++ip) {
      xpVec[ip*3] += dx[ip];
    }
    //getchar();
  }
  delete[] dx;

  char filename[128];
  //sprintf(filename,"pts.1D-%dto%d-long.pbin",2*n,n);  
  sprintf(filename,"pts.1D-%dto%d-long-smooth.pbin",2*n,n);  
  //sprintf(filename,"pts.1D-%dto%d-smooth.pbin",2*n,n);  
  //sprintf(filename,"pts.1D-%d.pbin",n);  
  writePointsToFile(filename,xpVec,deltaVec);
  
  return 0;

}
