#include "../../catch1/catch.hpp"
#include "../StaticSolver.hpp"
#include "../DataExchanger.hpp"

// local test setup
class SampleExchanger : public StaticSolver {

public:

  DataExchanger* exchanger;

  SampleExchanger() : exchanger(NULL) {}

  void initData() {}

  void init(const string zone_name) {

    StaticSolver::init();

    // now build the exchanger..

    // const string zone_name = "y0";
    BfZone* zone_ptr       = NULL;

    FOR_IZONE(bfZoneVec) {

      if ( bfZoneVec[izone].getName() == zone_name) {
        zone_ptr = &bfZoneVec[izone];
        break;
      }
    }

    if (zone_ptr == NULL) throw(-1);

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
    }
    delete[] x_tmp;
    delete[] x_ex;
  }

  void run() {}

  ~SampleExchanger() {
    if ( exchanger != NULL) delete exchanger;
  }
};

// TEST_CASE("data exchanger init with invalid zone","[core]") {
//   SampleExchanger solver;
//   CHECK_THROWS(solver1.init("y0"));
// }

TEST_CASE("data exchanger init","[core]") {
  SampleExchanger solver;
  CHECK_NOTHROW(solver.init("inlet"));
}
