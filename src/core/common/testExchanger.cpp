#include "CTI.hpp"
using namespace CTI;

#include "StaticSolver.hpp"
#include "DataExchanger.hpp"

class SampleExchanger : public StaticSolver { 

public: 

  DataExchanger* exchanger;


  SampleExchanger() : exchanger(NULL) {}

  void initData() {}

  void init() { 

    StaticSolver::init();

    // now build the exchanger.. 
    
    const string zone_name = "y0";
    BfZone* zone_ptr       = NULL;

    FOR_IZONE(bfZoneVec) { 

      if ( bfZoneVec[izone].getName() == zone_name) { 
        zone_ptr = &bfZoneVec[izone];
        break;
      }
    }

    assert( zone_ptr != NULL);

    // let's loko for a set of x points that are 1.5\Delta into the domain .. 

    const int nx      = zone_ptr->nbf;
    double (*x_ex)[3] = new double[zone_ptr->nbf][3];

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

      const double delta = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
      double unit_n[3]; 
      const double nmag  = MAG(zone_ptr->n_bf[ibf]);
      for (int i = 0; i < 3; ++i) 
        unit_n[i] = zone_ptr->n_bf[ibf][i]/nmag;

      for (int i = 0; i < 3; ++i) 
        x_ex[ibf][i] = zone_ptr->x_bf[ibf][i] - 3.0*delta*unit_n[i];

    }

    exchanger = new DataExchanger(this);
    exchanger->init(x_ex,nx);

    exchanger->checkMaxDist();


    double (*x_tmp)[3] = new double[zone_ptr->nbf][3];

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) 
      for (int i = 0; i < 3; ++i) 
        x_tmp[ibf][i] = HUGE_VAL;

    double (*grad_x)[3][3] = new double[ncv_g][3][3];
    FOR_ICV_G FOR_I3 FOR_J3 { 
      if (i == j) 
        grad_x[icv][i][j] = 1.0;
      else 
        grad_x[icv][i][j] = 0.0;
    }

    exchanger->exchangeCvDataWithGhosts(x_tmp,x_cv,grad_x);
    delete[] grad_x;

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

      for (int i = 0; i <3; ++i) {
        assert(abs(x_tmp[ibf][i] - x_ex[ibf][i]) < 1.0e-12);
      }

      /*
      for (int i = 0; i < 3; ++i) { 
        
        if ( x_tmp[ibf][i] == HUGE_VAL) { 
          cout << "uh oh: " << COUT_VEC(x_ex[ibf]) << endl;
          cout.flush();
        }

        assert(x_tmp[ibf][i] != HUGE_VAL);
      }
      */

    }

    delete[] x_tmp;
    delete[] x_ex;

  }


  void run() {}


  ~SampleExchanger() { 

    if ( exchanger != NULL) delete exchanger;

  }

};

int main(int argc, char * argv[]) {

  try {
    
    CTI_Init(argc,argv,"testExchanger.in");
    
    {
      
      SampleExchanger solver;
      solver.init();
      solver.run();
      
    }
    
    CTI_Finalize();
    
  }
  catch (int e) {
    if (e == 0)
      CTI_Finalize();
    else
      CTI_Abort();
  }
  catch(...) {
    CTI_Abort();
  }

  return(0);

}

