#include "CompData.hpp"

vector<vector<double> > getData(string file_name) {

  // this routime gets the data written in columns
  // ignores lines starting with #, and ignores blank lines
  // will return a vector containing columns of data
  ifstream fh (file_name.c_str());

  vector<vector<double> > Data;
  int ncol = 0;

  if (fh.is_open()) {
    string line;
    while ( getline(fh,line) ) {

      if (line[0] == '#') continue;
      if (line == "") continue;

      istringstream ss(line);

      // populate columns once
      if (Data.size()==0) {
        istringstream ss_check(line);
        string token;
        while ( getline(ss_check, token, ' ')) {
          if (token=="") continue;
          ncol++;
        }
        for (int icol = 0; icol < ncol; icol++) {
          vector<double> vec;
          Data.push_back(vec);
        }
      }

      int icol = 0;
      string token;
      // read line
      while ( getline(ss,token,' ') ) {
        if (token=="") continue;
        double num = atof(token.c_str());
        assert(icol<ncol);
        Data[icol++].push_back(num);
      }

    }
  } else {
    return Data;
  }

  fh.close();
  return Data;
}

int interpolate_index(double& x, vector<double>& x_vec) {

  // this routine finds an index in x_vec (ind) where x will be interpolated linearly between x_vec[ind] and x_vec[ind+1]

  if (x < x_vec[0])
    return 0;
  if (x >= x_vec[x_vec.size()-1])
    return x_vec.size()-2;

  // return the first index where x >= x_vec[index]
  for (int i = 0; i < x_vec.size()-1; i++) {
    if ((x >= x_vec[i]) and (x < x_vec[i+1])) {
      return i;
    }
  }
  
  return -1;
  
}

int checkData(vector<double>& x, vector<double>& y, vector<double>& x_ref, vector<double>& y_ref, const double tol) {

  // this routine checks the data (x,y) with reference (x_ref,y_ref), they don't have to be the same size.
  // returns 0: if they match
  // returns 1: if they don't match
  
  assert(x.size() == y.size());
  assert(x_ref.size() == y_ref.size());

  vector<double> error;
  for (int ip = 0; ip < x.size(); ip++) {
    int ip_match = interpolate_index(x[ip], x_ref);
    //cout << "ip_match: " << ip_match << " x: " << x[ip] << "x_ref: " << x_ref[ip_match] << " y: " << y[ip];
    if (ip_match >= 0) {
      double y_match = ( (x_ref[ip_match+1] - x[ip])*y_ref[ip_match] + (x[ip] - x_ref[ip_match])*y_ref[ip_match+1] ) / (x_ref[ip_match+1] - x_ref[ip_match]);
      //cout << " y_ref: " << y_ref[ip_match] << endl;
      error.push_back((y[ip]-y_match));
    }
  }

  if (error.size() != x.size())
    return 1;

  double rms_error = 0;
  for (int ie = 0; ie < error.size(); ie++) {
    //cout << "error: " << error[ie] << endl;
    rms_error += error[ie]*error[ie];
  }
  rms_error /= double(error.size());
  rms_error = sqrt(rms_error);
  
  // report the error, make sure to output this to a file not the terminal
  cout << " error: " << rms_error << " ";
  return (rms_error < tol) ? 0 : 1;
}

