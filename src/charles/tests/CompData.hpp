#ifndef __COMPDATA_HPP__
#define __COMPDATA_HPP__
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <assert.h>
#include <math.h>
using namespace std;

vector<vector<double> > getData(string file_name);
int checkData(vector<double>& x, vector<double>& y, vector<double>& x_ref, vector<double>& y_ref, const double tol);

#endif
