
#include <iostream>
#include <map>
#include <vector>
#include <assert.h>

using namespace std;

#define I_DATA 0
#define D_DATA 1
#define D3_DATA 2


class LpRegister { 
public:

  int lp_stride;
  void** lp_ptr;
  string lp_name;
  int* lp_n_ptr;

  map<string,pair<int,int> > lp_data_map;

  LpRegister(): lp_stride(0), lp_ptr(NULL), lp_name(""), lp_n_ptr(NULL) {}

  template<class T>
  void registerLp(T * &lp,const string& name,int &n) {
    lp_ptr = (void**)&lp;
    lp_name = name;
    lp_n_ptr = &n;
    assert(sizeof(T)%sizeof(int) == 0);
    lp_stride = sizeof(T)/sizeof(int); // number of integers
  }

  template<class T>
  void registerLpData(T * &lp,double &d_data,const string& name) {
    assert(lp_ptr);
    const int offset = (&d_data)-(*(double**)lp_ptr);
    lp_data_map[lp_name+":"+name] = pair<int,int>(2*offset,D_DATA); // int offsets
  }
  
  template<class T>
  void registerLpData(T * &lp,double (&d3_data)[3],const string& name) {
    assert(lp_ptr);
    const int offset = ((double*)&d3_data)-(*(double**)lp_ptr);
    lp_data_map[lp_name+":"+name] = pair<int,int>(2*offset,D3_DATA); // int offsets
  }

  void report() { 

    for (map<string,pair<int,int> >::iterator it = lp_data_map.begin(); it != lp_data_map.end(); ++it)  
      cout << it->first << " : " << it->second.first  << "  , " << it->second.second << endl;

  }

};

class Particle { 
public:
  double T;
  double x[3];
  double u[3];
};

class ParticleContainer { 
public:

  vector<Particle> part_vec;
  Particle* part;
  int np;

  ParticleContainer() {
   
    np = 0;

    //part = new Particle[np];

    part_vec.resize(np);
    part  = &(part_vec.front());
    
  }

  ~ParticleContainer() { 
    //delete[] part;
  }

};

int main() { 

  LpRegister reg_class;
  ParticleContainer* pc = new ParticleContainer();

  reg_class.registerLp<Particle>(pc->part,"test_lp",pc->np);
  reg_class.registerLpData<Particle>(pc->part,pc->part->x,"x");
  reg_class.registerLpData<Particle>(pc->part,pc->part->u,"u");
  reg_class.registerLpData<Particle>(pc->part,pc->part->T,"T");

  reg_class.report();

  delete pc;

  return 0;
}

