
#include "Chemtable.hpp"

inline double randu() { 
  return 2.0*double(rand())/double(RAND_MAX)-1.0;
}

double lookupTest(AbstractChemtable2D* cpu_chem, double* Z_cpu,double* C_cpu,
                const int n,const string& name) { 
  
  double * val_cpu = new double[n];
  
  clock_t timer = clock(); 
  cpu_chem->lookup(val_cpu,name,Z_cpu,C_cpu,n);
  timer         = clock() - timer;
  delete[] val_cpu;

  return (double(timer))/CLOCKS_PER_SEC;
}

double lookupSpecialTest(AbstractChemtable2D* cpu_chem, 
                         double* Z_cpu, double* C_cpu0, double* C_cpu1, 
                         const int n, const string& name) { 

  double * val_cpu = new double[n];
  clock_t timer    = clock();
  cpu_chem->lookupSpecial(val_cpu,name, Z_cpu, C_cpu0, C_cpu1, n);
  timer            = clock() - timer;
  delete[] val_cpu;
  
  return (double(timer))/CLOCKS_PER_SEC;
}

void doit() { 

  AbstractChemtable2D* cpu_chem = NULL; 

  vector<string> strVec;
  strVec.push_back("rho");
  strVec.push_back("T");
  strVec.push_back("R");
  strVec.push_back("e");
  strVec.push_back("prog");
  strVec.push_back("src_prog");
  strVec.push_back("gamma");
  strVec.push_back("a_gamma");
  strVec.push_back("mu");
  strVec.push_back("a_mu");
  strVec.push_back("locp");
  strVec.push_back("a_locp");
  strVec.push_back("sL"); 
  strVec.push_back("lF"); 
  strVec.push_back("mw"); 
  strVec.push_back("int_rho_src"); 
  
  initChemtable(cpu_chem,getStringParam("CHEMTABLE"));

  cout << "======================================" << endl;
  cpu_chem->loadVariables(strVec); 
  cout << "======================================" << endl;

  // sweep through the Z space ...
  const int n          = 100;
  //const int n          = 65536;
  //const int n          = 1 << 20;
  const double C_mean  = 0.01; 
  double * Z_cpu       = new double[n];
  double * C_cpu       = new double[n];
  double * C_cpu1      = new double[n];

  for (int i =0; i < n; ++i) { 
    C_cpu[i]  = C_mean + 0.1*randu();
    C_cpu1[i] = C_mean + 0.1*randu(); 
  }
 
  const double dZ    = 0.01;
  const double Z_max = 1.0;
  double Z           = 0.2;
  
  cout << "--- Uniform Z lookup ----" << endl;
  while ( Z <= Z_max) {

    double time = 0.0;
    for (int i = 0; i < n; ++i) Z_cpu[i] = Z;

    for (int j = 0; j < 2; ++j) { 
      time += lookupTest(cpu_chem,Z_cpu,C_cpu,n,"T"); 
      time += lookupTest(cpu_chem,Z_cpu,C_cpu,n,"R");
      time += lookupTest(cpu_chem,Z_cpu,C_cpu,n,"mu");
      time += lookupSpecialTest(cpu_chem,Z_cpu,C_cpu,C_cpu1,n,"int_rho_src");
    }

    cout << "Z, time spent = " << Z << "   " << time << endl;
    Z += dZ;
  }

  cout << "--- Rand Z test -----" << endl;
  const double Z_mean = 0.1;
  for (int j = 0; j < 2; ++j) { 
    
    for (int i =0; i < n ; ++i) 
      Z_cpu[i] = Z_mean*(1.0 + randu());

    double time = 0.0;
    for (int k =0; k < 2; ++k) { 
      time       += lookupTest(cpu_chem,Z_cpu,C_cpu,n,"T");
      time       += lookupTest(cpu_chem,Z_cpu,C_cpu,n,"R");
      time       += lookupTest(cpu_chem,Z_cpu,C_cpu,n,"mu");
      time       += lookupSpecialTest(cpu_chem,Z_cpu,C_cpu,C_cpu1,n,"int_rho_src"); 
    }
    cout << "rand z iter, time spent = " << j << "    " << time << endl;
  }
} 

int main(int argc, char* argv[]) { 

  CTI_Init(argc,argv,"test_query.in");
  doit();
  cout << " done!" << endl;
  CTI_Finalize();
  return 0;
} 
