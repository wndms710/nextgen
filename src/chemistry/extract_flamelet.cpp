
#include "Flamelet.hpp"

void reformat(const string& fout, const string& fin, const vector<string>& vars) { 

  Flamelet f; 
  f.readFlameletFile(fin);
  f.cleanupNamesAndInitSpeciesNames(); 

  cout << " > flamelet np " << f.npoints << endl;

  map<string,int> var_map;
  for (unsigned int i =0; i < f.varNameVec.size(); ++i) {
    cout << "var found: " << f.varNameVec[i] << endl;
    for (unsigned int j =0; j < vars.size(); ++j) { 
      if ( vars[j] == f.varNameVec[i]) 
        var_map[vars[j]] = i;
    }
  }

  for (unsigned int j =0; j < vars.size(); ++j) { 
    if ( var_map.find(vars[j]) == var_map.end()) { 
      cerr << "couldnt find : " << vars[j] << endl; 
      throw(0);
    }
  }

  ofstream out(fout.c_str());
  out << "# ";
  for (unsigned int i =0; i < vars.size(); ++i) 
    out << vars[i] << "   "; 
  out << endl;

  for (int i = 0; i < f.npoints; ++i) { 
    for (unsigned int j =0; j < vars.size(); ++j) { 
      const int ivar = var_map[vars[j]]; 
      out << f.varDataVec[ivar][i] << "    " ;
    }
    out << endl;
  }//i

  out.close();
}

int main(const int argc, const char* argv[]) { 
  vector<string> vars; 
  vars.push_back("y"); 
  vars.push_back("Y_CO");
  vars.push_back("Y_CO2");
  vars.push_back("Y_H2O");
  vars.push_back("src_CO");
  vars.push_back("s");
  vars.push_back("T");
  vars.push_back("TotalEnthalpy");

  const string fname(argv[1]);
  reformat("tmp.dat",fname,vars);
  return 0;
}
