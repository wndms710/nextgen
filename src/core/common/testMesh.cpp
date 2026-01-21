#include "CTI.hpp"
using namespace CTI;

#include "StaticSolver.hpp"

class SampleSolver : public StaticSolver {
public:

  SampleSolver() {
    COUT1("SampleSolver()");
  }

  void initCtiRegisterManagedData() {}

  void initData() {}

  void init() {
    StaticSolver::init(INIT_COMPACT_FACES);
  }

  void run() {
    if (checkParam("INTERACTIVE")) kfr.setHold(true);

    COUT1("Running mesh statistical analysis:");

    // data for histograms
    double * cv_nfa = new double[ncv];
    double * cv_volRatio = new double[ncv];
    double * cv_orth_angle = new double[ncv];
    double * wallDelta = new double[nbf];
    double * wall_orth_angle = new double[nbf];
    double * nbr_distRatio = new double[ncv];
    // FOR_IFA {
    //   orth_angle[ifa] = HUGE_VAL;
    // }
    FOR_IBF {
      wallDelta[ibf] = HUGE_VAL;
      wall_orth_angle[ibf] = HUGE_VAL;
    }

    // min/max/mean for each of the above + vol_cv
    double (*my_range)[2] = new double[7][2];
    double * my_sum = new double[7];
    FOR_I6 {
      my_range[i][0] = HUGE_VAL;
      my_range[i][1] = HUGE_VAL;  // store -max here so same MPI_MIN can be used in comm
      my_sum[i] = 0.0;
    }

    COUT1(" > looping over cells to compute cell and face-based stats");
    for (int icv=0; icv < ncv; ++icv) {
      cv_nfa[icv] = faocv_i[icv+1]-faocv_i[icv];  // internal face neighbors

      // internal face count
      my_range[0][0] = min(my_range[0][0],cv_nfa[icv]);
      my_range[0][1] = min(my_range[0][1],-cv_nfa[icv]);
      my_sum[0] += cv_nfa[icv];

      // cell volume
      my_range[5][0] = min(my_range[5][0],vol_cv[icv]);
      my_range[5][1] = min(my_range[5][1],-vol_cv[icv]);
      my_sum[5] += vol_cv[icv];

      double nbrVol[2] = {HUGE_VAL,-HUGE_VAL};
      double nbrDist[2] = {HUGE_VAL,-HUGE_VAL};
      double orthAngle_area_sum[2] = {0.0,0.0};
      for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
        const int ifa = faocv_v[foc];
        const int icv_nbr = (cvofa[ifa][0] == icv) ? cvofa[ifa][1]:cvofa[ifa][0];
        assert(icv != icv_nbr);

        // volume ratio
        nbrVol[0] = min(vol_cv[icv_nbr],nbrVol[0]);
        nbrVol[1] = max(vol_cv[icv_nbr],nbrVol[1]);

        // orthogonality requires building cv-cv vector
        double cv_axis[3] = DIFF(x_cv[icv_nbr],x_cv[icv]);
        NORMALIZE(cv_axis);
        double vv_axis[3] = DIFF(x_vv[icv_nbr],x_vv[icv]);
        const double d2nbr = MAG(vv_axis);

        // nbr distance
        nbrDist[0] = min(d2nbr,nbrDist[0]);
        nbrDist[1] = max(d2nbr,nbrDist[1]);

        NORMALIZE(vv_axis);
        double dp = DOT_PRODUCT(vv_axis,cv_axis);
        if (dp < -1.0) {
          dp = -(dp+2.0);
        }
        else if (dp > 1.0) {
          dp = 2.0-dp;
        }

        // face-area weighted orthogonality
        const double area = MAG(n_fa[ifa]);
        orthAngle_area_sum[0] += area*(90.0-(acos(dp) * 180.0 / M_PI));
        orthAngle_area_sum[1] += area;
      }
      // volume ratio
      assert(nbrVol[0] != 0.0);
      cv_volRatio[icv] = nbrVol[1]/nbrVol[0];
      my_range[1][0] = min(my_range[1][0],cv_volRatio[icv]);
      my_range[1][1] = min(my_range[1][1],-cv_volRatio[icv]);
      my_sum[1] += cv_volRatio[icv];

      // nbr distance ratio
      assert(nbrDist[0] != 0.0);
      nbr_distRatio[icv] = nbrDist[1]/nbrDist[0];
      my_range[6][0] = min(my_range[1][0],nbr_distRatio[icv]);
      my_range[6][1] = min(my_range[1][1],-nbr_distRatio[icv]);
      my_sum[6] += nbr_distRatio[icv];

      // orthogonality angle
      assert(orthAngle_area_sum[1] != 0.0);
      cv_orth_angle[icv] = orthAngle_area_sum[0]/orthAngle_area_sum[1];
      my_range[2][0] = min(my_range[2][0],cv_orth_angle[icv]);
      my_range[2][1] = min(my_range[2][1],-cv_orth_angle[icv]);
      my_sum[2] += cv_orth_angle[icv];

      // boundary stuff
      for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
        const int ibf = bfocv_v[boc];
        wallDelta[ibf] = min(wallDelta[ibf],area_bf[ibf]/area_over_delta_bf[ibf]);
        my_range[3][0] = min(my_range[3][0],wallDelta[ibf]);
        my_range[3][1] = min(my_range[3][1],-wallDelta[ibf]);
        my_sum[3] += wallDelta[ibf];

        // wall orthognality
        double cv_axis[3] = DIFF(x_bf[ibf],x_cv[icv]);
        const double cv_axis_mag = MAG(cv_axis);
        if ((area_bf[ibf] > 0.0) && (cv_axis_mag > 0.0)) {
          const double unit_n[3] = {
            n_bf[ibf][0]/area_bf[ibf],
            n_bf[ibf][1]/area_bf[ibf],
            n_bf[ibf][2]/area_bf[ibf]
          };
          NORMALIZE(cv_axis);

          double dp = DOT_PRODUCT(unit_n,cv_axis);
          if (dp < -1.0) {
            dp = -(dp+2.0);
          }
          else if (dp > 1.0) {
            dp = 2.0-dp;
          }

          wall_orth_angle[ibf] = min(wall_orth_angle[ibf],90.0-(acos(dp) * 180.0 / M_PI));
        }
        else {
          wall_orth_angle[ibf] = min(wall_orth_angle[ibf],90.0);
        }
        my_range[4][0] = min(my_range[4][0],wall_orth_angle[ibf]);
        my_range[4][1] = min(my_range[4][1],-wall_orth_angle[ibf]);
        my_sum[4] += wall_orth_angle[ibf];
      }
    }

    COUT1(" > communicating information and generating histograms");
    Histogram hist_nfa(3.5,33.5,30,mpi_comm);
    hist_nfa.setIntBins(false);
    hist_nfa.add(cv_nfa,ncv);
    DELETE(cv_nfa);
    hist_nfa.reduce();
    hist_nfa.write("hist_nfa.dat");

    Histogram hist_vol(mpi_comm,50);
    hist_vol.setLog(true);
    hist_vol.setIntBins(false);
    hist_vol.add(vol_cv,ncv);
    hist_vol.reduce();
    hist_vol.write("hist_volCv.dat");

    Histogram hist_volRatio(mpi_comm,50);
    hist_volRatio.setIntBins(false);
    hist_volRatio.add(cv_volRatio,ncv);
    DELETE(cv_volRatio);
    hist_volRatio.reduce();
    hist_volRatio.write("hist_volRatio.dat");

    Histogram hist_orthAngle(59.5,90.5,30,mpi_comm);
    hist_orthAngle.setIntBins(false);
    hist_orthAngle.add(cv_orth_angle,ncv);
    DELETE(cv_orth_angle);
    hist_orthAngle.reduce();
    hist_orthAngle.write("hist_orthAngle.dat");

    Histogram hist_wallOrthAngle(29.5,90.5,60,mpi_comm);
    hist_wallOrthAngle.setIntBins(false);
    hist_wallOrthAngle.add(wall_orth_angle,nbf);
    DELETE(wall_orth_angle);
    hist_wallOrthAngle.reduce();
    hist_wallOrthAngle.write("hist_wallOrthAngle.dat");

    Histogram hist_wallDelta(mpi_comm,50);
    hist_wallDelta.setLog(true);
    hist_wallDelta.setIntBins(false);
    hist_wallDelta.add(wallDelta,nbf);
    DELETE(wallDelta);
    hist_wallDelta.reduce();
    hist_wallDelta.write("hist_wallDelta.dat");

    Histogram hist_distRatio(mpi_comm,50);
    hist_distRatio.setIntBins(false);
    hist_distRatio.add(nbr_distRatio,ncv);
    DELETE(nbr_distRatio);
    hist_distRatio.reduce();
    hist_distRatio.write("hist_nbrDistRatio.dat");

    COUT1(" > performing communication for individual metric stats");
    double (*range)[2] = NULL;
    double * sum = NULL;
    int8 * counts = NULL;
    if (mpi_rank == 0) {
      range = new double[7][2];
      sum = new double[7];
      counts = new int8[7];
    }
    int8 my_counts[7] = {int8(ncv),int8(ncv),int8(ncv),int8(nbf),int8(nbf),int8(ncv),int8(ncv)};
    MPI_Reduce(&my_counts,counts,7,MPI_INT8,MPI_SUM,0,mpi_comm);
    MPI_Reduce(my_range,range,7*2,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
    MPI_Reduce(my_sum,sum,7,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

    if (mpi_rank == 0) {
      // report out min/max/stats
      COUT1("    > nfa-per-cell              (min:mean:max) " << range[0][0] << " : " << sum[0]/double(counts[0]) << " : " << -range[0][1]);
      COUT1("    > volume-ratio-per-cell     (min:mean:max) " << range[1][0] << " : " << sum[1]/double(counts[1]) << " : " << -range[1][1]);
      COUT1("    > orthog-angle-per-cell     (min:mean:max) " << range[2][0] << " : " << sum[2]/double(counts[2]) << " : " << -range[2][1]);
      COUT1("    > wall-delta-per-ibf        (min:mean:max) " << range[3][0] << " : " << sum[3]/double(counts[3]) << " : " << -range[3][1]);
      COUT1("    > wall-orthog-angle-per-ibf (min:mean:max) " << range[4][0] << " : " << sum[4]/double(counts[4]) << " : " << -range[4][1]);
      COUT1("    > cell-vol                  (min:mean:max) " << range[5][0] << " : " << sum[5]/double(counts[5]) << " : " << -range[5][1]);
      COUT1("    > nbr-dist-ratio            (min:mean:max) " << range[6][0] << " : " << sum[6]/double(counts[6]) << " : " << -range[6][1]);
    }

  }

  virtual ~SampleSolver() {}

  void initBoundaryConditions();
};

void SampleSolver::initBoundaryConditions() {
  StaticSolver::initBoundaryConditions();
}

int main(int argc, char * argv[]) {

  try {

    CTI_Init(argc,argv,"testMesh.in");

    {

      //MPI_Run_check();

      SampleSolver solver;
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
