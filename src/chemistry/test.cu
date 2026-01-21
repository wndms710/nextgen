
#include "Chemtable.hpp"
#include "ChemtableGPU.hpp"

inline double randu() { 
  return 2.0*double(rand())/double(RAND_MAX)-1.0;
}


void lookupTest(AbstractChemtable2D* cpu_chem, AbstractChemtable2D* gpu_chem, 
                double* Z_cpu,double* C_cpu,double* Z_gpu,double* C_gpu,
                const int n,const string& name) { 
  
  double * val_cpu = new double[n];
  double * val_gpu = new double[n];
 
  cpu_chem->lookup(val_cpu,name,Z_cpu,C_cpu,n);
  gpu_chem->lookup(val_gpu,name,Z_gpu,C_gpu,n);

  int istart = 0;
  for (int i =istart; i <istart+5 ; ++i) 
    cout << "i, val_cpu, val_gpu = " << i << "    " << val_cpu[i] << "   " << val_gpu[i] << endl; 

  double max_diff = 0.0;
  for (int i =0; i<n; ++i) { 
    //assert(fabs(val_gpu[i]-val_cpu[i]) < 1.0e-08); 
    max_diff = max(abs(val_gpu[i]-val_cpu[i]),max_diff);
  }

  cout << " Max diff : " << max_diff << endl;
  if ( max_diff > 1.0e-08) 
    assert(0);

  cout << "Passed lookup for " << name << endl; 
  delete[] val_cpu;
  delete[] val_gpu;
}

void lookupSpecialTest(AbstractChemtable2D* cpu_chem, AbstractChemtable2D* gpu_chem,
                       double* Z_cpu, double* C_cpu0, double* C_cpu1, 
                       double* Z_gpu, double* C_gpu0, double* C_gpu1, 
                       const int n, const string& name) { 

  double * val_cpu = new double[n];
  double * val_gpu = new double[n];

  cpu_chem->lookupSpecial(val_cpu,name, Z_cpu, C_cpu0, C_cpu1, n); 
  
  CartesianChemtable2dGpu* gpu_chem_ = dynamic_cast<CartesianChemtable2dGpu*>(gpu_chem);
  gpu_chem_->lookupSpecial(val_gpu,name, Z_gpu, C_gpu0, C_gpu1, n);
  gpu_chem_->lookupSpecialFinish(val_gpu,name, n);
 
  int istart = 0;
  for (int i =istart; i <istart+5 ; ++i) 
    cout << "i, val_cpu, val_gpu = " << i << "    " << val_cpu[i] << "   " << val_gpu[i] << endl; 

  for(int i =0; i<n ; ++i) 
    assert( fabs(val_gpu[i]-val_cpu[i]) < 1.0e-08);

  cout << "Passed lookup special for " << name << endl;

  delete[] val_cpu;
  delete[] val_gpu;
}


void lookupReducedTest(AbstractChemtable2D* cpu_chem, AbstractChemtable2D* gpu_chem, 
                       double * Z_cpu, double * Z_gpu, const int n, const string& name) { 

  double * val_cpu = new double[n];
  double * val_gpu = new double[n];

  cpu_chem->lookupReduced(val_cpu,name,Z_cpu,n);
  gpu_chem->lookupReduced(val_gpu,name,Z_gpu,n);

  int istart = 0;
  for (int i =istart; i <istart+5 ; ++i) 
    cout << "i, val_cpu, val_gpu = " << i << "    " << val_cpu[i] << "   " << val_gpu[i] << endl; 
  
  double max_diff = 0.0;
  for (int i =0; i<n; ++i) { 
    //assert(fabs(val_gpu[i]-val_cpu[i]) < 1.0e-08); 
    max_diff = max(abs(val_gpu[i]-val_cpu[i]),max_diff);
  }

  cout << " Max diff : " << max_diff << endl;
  if ( max_diff > 1.0e-08) 
    assert(0);

  cout << "Passed lookup for " << name << endl; 
  delete[] val_cpu;
  delete[] val_gpu;
}


// serial test, no mpi.  
void doit() { 

  AbstractChemtable2D* cpu_chem = NULL; 
  AbstractChemtable2D* gpu_chem = NULL; 

  // premixed chemistry vars ... 
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
  initChemtableGpu(gpu_chem,getStringParam("CHEMTABLE"));

  cout << "======================================" << endl;
  cpu_chem->loadVariables(strVec); 
  cout << "======================================" << endl;
  gpu_chem->loadVariables(strVec);
  cout << "======================================" << endl;

  CartesianChemtable2dGpu* gpu_chem_ = dynamic_cast<CartesianChemtable2dGpu*>(gpu_chem);
  gpu_chem_->initGpu();

  // randomnly populate a Z, C for the test..
  const int n         = 1221; 
  const double Z_mean = 0.035; 
  const double C_mean = 0.08;
  
  double * Z_cpu      = new double[n];
  double * C_cpu      = new double[n];
  double * C_cpu1     = new double[n];

  for (int i =0; i<n; ++i) { 
    Z_cpu[i] = Z_mean + 0.02*randu();
    C_cpu[i] = C_mean + 0.06*randu();
    C_cpu1[i]= C_mean + 0.06*randu();
  }

  double *Z_gpu, *C_gpu, *C_gpu1;
  cudaMalloc((void**)&Z_gpu,n*sizeof(double)); 
  cudaMalloc((void**)&C_gpu,n*sizeof(double));
  cudaMalloc((void**)&C_gpu1,n*sizeof(double));
  cudaMemcpy(Z_gpu,Z_cpu,n*sizeof(double), cudaMemcpyHostToDevice); 
  cudaMemcpy(C_gpu,C_cpu,n*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(C_gpu1,C_cpu1,n*sizeof(double), cudaMemcpyHostToDevice); 

  // now run some lookup tests... 
  lookupTest(cpu_chem,gpu_chem,Z_cpu,C_cpu,Z_gpu,C_gpu,n,"T"); 
  lookupTest(cpu_chem,gpu_chem,Z_cpu,C_cpu,Z_gpu,C_gpu,n,"R"); 
  lookupTest(cpu_chem,gpu_chem,Z_cpu,C_cpu,Z_gpu,C_gpu,n,"mu"); 

  lookupSpecialTest(cpu_chem, gpu_chem, Z_cpu, C_cpu, C_cpu1, 
                    Z_gpu, C_gpu, C_gpu1, n, "int_rho_src"); 

  lookupReducedTest(cpu_chem,gpu_chem,Z_cpu,Z_gpu,n,"lF");
  lookupReducedTest(cpu_chem,gpu_chem,Z_cpu,Z_gpu,n,"sL");
  lookupReducedTest(cpu_chem,gpu_chem,Z_cpu,Z_gpu,n,"upper");

  cudaFree(Z_gpu); 
  cudaFree(C_gpu);
  cudaFree(C_gpu1);
} 


int main(int argc, char* argv[]) { 

  CTI_Init(argc,argv,"test.in");
  cudaSetDevice(1);
  doit();
  cout << " done!" << endl;
  CTI_Finalize();
  return 0;
} 
