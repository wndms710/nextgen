#include "CompData.hpp"

int main(int argc, char* argv[]) {

  // how to run this: exe -i inputFile.dat -r referenceFile.dat -t tolerance

  string file_input;
  string file_ref;
  double tol = 1e-8;

  // read arguments
  int iarg = 1;
  while(iarg < argc) {
    if (string(argv[iarg])=="-i") {
      file_input = argv[++iarg];
    } else if (string(argv[iarg])=="-r") {
      file_ref = argv[++iarg];
    } else if (string(argv[iarg])=="-t") {
      tol = atof(argv[++iarg]);
    }
    iarg++;
  }

  // check files are good
  {
    ifstream f_i (file_input.c_str());
    ifstream f_r (file_ref.c_str());
    if (not f_i.good())
      return 1;
    if (not f_r.good())
      return 1;
  }

  vector<vector<double> > Data_i = getData(file_input);
  vector<vector<double> > Data_r = getData(file_ref);

  if (Data_i.size() == 0) return 1;
  if (Data_r.size() == 0) return 1;

  return (checkData(Data_i[0], Data_i[1], Data_r[0], Data_r[1], tol)) ? 1 : 0;
}
