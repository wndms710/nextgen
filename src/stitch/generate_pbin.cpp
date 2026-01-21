
#include <iostream>
#include <cstdio>
#include <assert.h>
#include <climits>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "Adt.hpp"

using namespace std;

//
// the format of the pbin file that is passed to stitch 
// is as follows: 
// 
//     int8 ibuf[4] = {POINTS_IO_MAGIC_NUMBER,   // to check endianness 
//                     POINTS_IO_VERSION,        // version number 
//                     np_global,                // number of points, 
//                     1};                       // indicates that the delta record is present
// 
//     double (*xp)[3];                          // point locations
//
//     double *delta;                            // measure of distance between neighboring points
//
// 
// the input file format is a space delimited list of x,y,z coordinates
// 
//       x_0    y_0    z_0 
//       x_1    y_1    z_1 
//             ... 
//       x_{n-1} y_{n-1} z_{n-1}

typedef long long int int8;

#define POINTS_IO_VERSION 2
#define POINTS_IO_MAGIC_NUMBER 1235813
#define HUGE_INT (1<<24)

class Coords { 
public: 
  double x[3];

  bool operator< (const Coords& _c) const { 

    if ( x[0] == _c.x[0]) { 

      if ( x[1] == _c.x[1]) { 
        return (x[2] < _c.x[2]); 
      } else { 
        return (x[1] < _c.x[1]);
      }
    } else { 
      return (x[0] < _c.x[0]);
    }
  }
};

void readPointsFromFile(double (*&xp)[3], int&np, const char* filename) { 

  ifstream fin(filename);
  
  if ( !fin.is_open() ) { 

    cerr << " > unable to open file with points: " << filename << endl;
    cerr << " > please check the file path ... " << endl;
    throw(0);

  }


  // grab a count of the number of points (lines in the file) .. 

  np = 0;
  string line_str;
  while ( getline(fin,line_str)) 
    ++np;

  cout << " > number of points in generating file: " << np << endl;

  fin.clear();
  fin.seekg(0);

  vector<Coords> coords;
  coords.resize(np);

  for (int ip = 0; ip < np; ++ip) { 

    for (int i = 0; i < 3; ++i) 
      fin >> coords[ip].x[i];

  }

  fin.close();

  cout << " finished reading file. " << endl;
  cout << " starting to sort .. " << endl;

  sort(coords.begin(), coords.end());

  cout << " finished sort. " << endl;

  assert( xp == NULL);
  xp = new double[np][3];

  for (int ip =0; ip < np; ++ip) 
    for (int i = 0; i < 3; ++i) 
      xp[ip][i] = coords[ip].x[i];

}

double getLengthScaleFromPoints(const double (*xp)[3], const int np) { 

  //
  // this can be done much more intelligently, but for now, we are 
  // going to get a global length scale associated with bbox volumes.
  //
  
  double bbmin[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL};
  double bbmax[3] = {-HUGE_VAL,-HUGE_VAL,-HUGE_VAL};

  for (int ip = 0; ip < np; ++ip) {
    for (int i = 0; i < 3; ++i) { 
      bbmin[i] = fmin(bbmin[i], xp[ip][i]);
      bbmax[i] = fmax(bbmax[i], xp[ip][i]);
    }
  }

  const double bbox_vol = (bbmax[0] - bbmin[0])*
                          (bbmax[1] - bbmin[1])*
                          (bbmax[2] - bbmin[2]);

  assert( bbox_vol > 0.0);

  const double L_estimate = pow(bbox_vol/double(np),1.0/3.0);
  cout << " > using estimate for local length scale: " << L_estimate << endl;

  return L_estimate;

}

void setDelta2(double* &delta, const double (*xp)[3], const int np, const double L_est) { 

  cout << " > setDelta2() " << endl;

  assert( delta == NULL); delta = new double[np];

  for (int ip = 0; ip < np; ++ip) 
    delta[ip] = -1.0;

  
  cout << " > building Adt ... " ; cout.flush(); 
  Adt<double> * xp_adt = new Adt<double>(np,xp,xp);
  cout << " ok. " << endl; cout.flush();

  int * nbr_count = new int[np];
  int buf[2] = {HUGE_INT,-HUGE_INT};

  const int lb_nbr = 3;
  const int ub_nbr = 50;

  for (int ip = 0; ip < np; ++ip) { 
 
    if ( ip%(np/1000) == 0) 
      cout << " ... working on " << ip << " of : " << np << endl; 

    double L_guess;
   
    if ( delta[ip] < 0.0) 
      L_guess = L_est;
    else 
      L_guess = delta[ip];

    vector<int> nbr_vec;
    nbr_vec.clear();

    int iter           = 0;
    bool done          = false;
    const int max_iter = 100;
  
    while ( !done) { 

      // copy down the previous solutions 

      int nnbr_prev      = nbr_vec.size();
      double Lguess_prev = L_guess;

      if ( nnbr_prev <= lb_nbr ) { 

        if ( iter != 0) 
          L_guess *= M_PI*0.5;

      } else if ( nnbr_prev >= ub_nbr ) { 

        L_guess *= 0.5;

      } else { 
        
        // shouldnt get here ... 

        cout << " iter, nnbr_prev : " << iter << "    " << nnbr_prev << "   " << L_guess << endl;

        assert(0);
      }
        
      // assume we arent done .. 

      done = false;
      nbr_vec.clear();
      //xp_adt->buildListForSphere(nbr_vec,xp[ip],L_guess);
     
      double bbmin[3], bbmax[3];
      for (int i = 0; i < 3; ++i) { 
        bbmin[i] = xp[ip][i] - L_guess;
        bbmax[i] = xp[ip][i] + L_guess;
      }

      xp_adt->buildListForBBox(nbr_vec,bbmin,bbmax);

      if ( (nbr_vec.size() >= lb_nbr) && ( nbr_vec.size() <  ub_nbr)) { 

        // this guess is okay .. 

        done = true;

      } else if ( (nbr_vec.size() > ub_nbr) && (nnbr_prev > lb_nbr)) { 

        // the previous guess is okay.. 

        L_guess = Lguess_prev;
        done    = true;

      } else if ( (nbr_vec.size() < lb_nbr) && (nnbr_prev > lb_nbr)) { 

        // the previous guess is okay... 

        L_guess = Lguess_prev;
        done    = true;

      }

      ++iter;

      if ( iter == max_iter) { 
        cerr << " unable to find solution: " << nnbr_prev << "   " << nbr_vec.size() << "   " << L_guess << endl;
        throw(0);
      }
    }

    nbr_count[ip] = nbr_vec.size();
    buf[0]        = min(buf[0], nbr_count[ip]);
    buf[1]        = max(buf[1], nbr_count[ip]);
    delta[ip]     = L_guess;

    // set the delta inside of the nbr_vec as an initial guess as well... 
    
    for (int ii = 0, nn = nbr_vec.size(); ii < nn; ++ii) { 

      const int i_nbr = nbr_vec[ii];

      if ( i_nbr == ip) 
        continue;

      if ( delta[i_nbr] < 0.0) 
        delta[i_nbr] = L_guess;
      else 
        delta[i_nbr] = min(delta[i_nbr],L_guess);

    }
  }

  cout << " > range of nbr count : " << buf[0] << " , " << buf[1] << endl;

  delete[] nbr_count;
  delete xp_adt;

}

void setDelta(double* &delta, const double (*xp)[3], const int np, const double L_est) {

  assert( delta == NULL); delta = new double[np];

  Adt<double> * xp_adt = new Adt<double>(np,xp,xp);


  // we need to populate a guess for how far apart the neighbors are -- this is 
  // potentially expensive, but we'll just inflate or deflate sphere until we find 
  // an "acceptable" number of neighbors ... 

  int* nbr_count = new int[np];
  int buf[2] = {HUGE_INT, -HUGE_INT};

  const int lb_nbr = 3;
  const int ub_nbr = 50;

  for (int ip = 0; ip < np; ++ip) { 

   
    double L_guess = L_est;
    vector<int> nbr_vec;
    nbr_vec.clear();

    int iter           = 0;
    bool done          = false;
    const int max_iter = 100;
  
    while ( !done) { 

      // copy down the previous solutions 

      int nnbr_prev      = nbr_vec.size();
      double Lguess_prev = L_guess;

      if ( nnbr_prev <= lb_nbr ) { 

        if ( iter != 0) 
          L_guess *= M_PI*0.5;

      } else if ( nnbr_prev >= ub_nbr ) { 

        L_guess *= 0.5;

      } else { 
        
        // shouldnt get here ... 

        cout << " iter, nnbr_prev : " << iter << "    " << nnbr_prev << "   " << L_guess << endl;

        assert(0);
      }
        
      // assume we arent done .. 

      done = false;
      nbr_vec.clear();
      xp_adt->buildListForSphere(nbr_vec,xp[ip],L_guess);


      if ( (nbr_vec.size() > lb_nbr) && ( nbr_vec.size() <  ub_nbr)) { 

        // this guess is okay .. 

        done = true;

      } else if ( (nbr_vec.size() > ub_nbr) && (nnbr_prev > lb_nbr)) { 

        // the previous guess is okay.. 

        L_guess = Lguess_prev;
        done    = true;

      } else if ( (nbr_vec.size() < lb_nbr) && (nnbr_prev > lb_nbr)) { 

        // the previous guess is okay... 

        L_guess = Lguess_prev;
        done    = true;

      }

      ++iter;

      if ( iter == max_iter) { 
        cerr << " unable to find solution: " << nnbr_prev << "   " << nbr_vec.size() << "   " << L_guess << endl;
        throw(0);
      }
    }

    nbr_count[ip] = nbr_vec.size();
    buf[0]        = min(buf[0], nbr_count[ip]);
    buf[1]        = max(buf[1], nbr_count[ip]);
    delta[ip]     = L_guess;
  }

  cout << " > range of nbr count : " << buf[0] << " , " << buf[1] << endl;

  delete[] nbr_count;
  delete xp_adt;
}


void writePointsToFile( const double (*xp)[3], const double* delta, const int np, const char* filename) { 

  FILE* fp = fopen(filename, "wb");

  if ( !fp) { 

    cerr << " > unable to open output file: " << filename << endl;
    throw(0);

  }

  int8 ibuf[4] = {POINTS_IO_MAGIC_NUMBER, POINTS_IO_VERSION, np, 1};

  fwrite(ibuf       , sizeof(int8), 4, fp);
  fwrite((double*)xp, sizeof(double), np*3, fp);
  fwrite(delta      , sizeof(double), np  , fp);
  fclose(fp);

}

int main(const int argc, const char* argv[]) { 

  double L_est = -1.0;

  if ( argc < 2 ) { 
    
    cerr << " Usage: ./generate_pbin.exe <file-containing-points> [<approx local length scale>] "<< endl;
    return 1;
    
  }

  try { 

    const int nbr_type = 1;  // 0 - exhaustive, 1 - approximate
    int np             = 0;
    double (*xp)[3]    = NULL;
    double* delta      = NULL;
    
    readPointsFromFile(xp,np,argv[1]);
 
    if ( argc == 2) { 

      L_est = getLengthScaleFromPoints(xp,np);

    } else { 

      L_est = atof(argv[2]);
      
      if ( L_est <= 0.0) { 

        cerr << " > must specify a postive length scale if providing one; got : " << L_est << endl;
        throw(0);

      }
    }

    if ( nbr_type == 0) 
      setDelta(delta,xp,np,L_est);
    else if ( nbr_type == 1)
      setDelta2(delta,xp,np,L_est);
    else 
      assert(0);

    writePointsToFile(xp,delta,np,"points.pbin");


    if ( xp != NULL) delete[] xp;
    if ( delta != NULL) delete[] delta;

  } catch(...) { 

    cerr << " Unhandled exception -- see error message above.  Exiting... " << endl;
    return 1;

  }

  return 0;

}
