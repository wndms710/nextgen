#ifndef _MULTIGRID_HPP_
#define _MULTIGRID_HPP_

#include "MpiStuff.hpp"
#include "MiscUtils.hpp"
#include "DistributedDataExchanger.hpp"
#include "Utils.hpp"
#include "Prcomm.hpp"

namespace Multigrid {

  void colorCvsPadt(int8 *color_cv,int8 &ncolors,const double (*const x_cv)[3],const int8* const icv_global,
      const int ncv,const int8 ncv_global);
  void splitOrphanedColors(int8 * color_cv,int8 &ncolor,const int* const cvocv_i,const int* const cvocv_v,
      const int8* const icv_global,const uint8* const rbicv_g,const uint8* const rbicv_g2,const int ncv,const int ncv_g,
      const int ncv_g2,const int8 ncv_global);

  class AlgebraicCoarseGrid {
    public:

      int icg; // icg = 0 is the original grid

      // grid stuff...
      // cv data, cv == parent cell
      // cc data, cc == coarse cell
      
      const double *vol_cv;
      int8 ncc_global;
      int ncc_in; // internal
      int ncc_a; // +active cc's
      int ncc;   // +inactive cc's
      int ncc_g; // +ghosts cc's (note that inactive are treated as ghosts too)
      int8* icc_global;
      uint8* rbicc_i;  // rank-bits-icc inactive (for pack)
      uint8* rbicc_g; // rank-bits-icc ghost (for unpack)
      double (*x_cc)[3];
      double *vol_cc;
      double *inv_vol;
      int *ccocc_i;
      int *ccocc_v;
      double (*ccocc_grad_coeff)[3];
      int *cvocc_i; // coarse to fine
      int *cvocc_v;
      int *ccocv; // fine to coarse

      vector<Prcomm> ccIPrcommVec; // reduce inactive cc's onto active cc's
      vector<Prcomm> ccPrcommVec; // replace inactive/ghost cc's with active cc's
      map<const void*,MpiRequestStuff*> mpiRequestMap;

      // solver stuff...
      
      double* A_cc;
      double* inv_diag;
      double* err;
      double* res;
      double amg_coeff; // multiplies coefficients in coarse matrix

      AlgebraicCoarseGrid() {
        vol_cv = NULL;
        inv_vol = NULL;
        icg = 0;
        ncc_global = 0;
        ncc_in = 0;
        ncc_a = 0;
        ncc = 0;
        ncc_g = 0;
        icc_global = NULL;
        rbicc_i = NULL;
        rbicc_g = NULL;
        x_cc = NULL;
        vol_cc = NULL;
        ccocc_i = NULL;
        ccocc_v = NULL;
        ccocc_grad_coeff = NULL;
        cvocc_i = NULL;
        cvocc_v = NULL;
        ccocv = NULL;
        A_cc = NULL;
        inv_diag = NULL;
        err = NULL;
        res = NULL;
        amg_coeff = 1.0;
      }

      ~AlgebraicCoarseGrid() {
        clear();
      }

      void clear() {
        // icg 0 just copies base grid ptrs for some data...
        if (icg > 0) {
          DELETE(icc_global);
          DELETE(rbicc_i);
          DELETE(rbicc_g);
          DELETE(x_cc);
          DELETE(vol_cc);
          DELETE(ccocc_i);
          DELETE(ccocc_v);
          DELETE(ccocc_grad_coeff);
          DELETE(cvocc_i);
          DELETE(cvocc_v);
          DELETE(ccocv);
          DELETE(A_cc);
        }
        DELETE(inv_vol);
        DELETE(inv_diag);
        DELETE(res);
        DELETE(err);
      }

      void buildCcIPrcomm();

      // updateCcIData(double * s and double (*s)[3]...

#define T double
#define UNPACK_BUF unpack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 52122
#include "updateCcIData.hpp"
#undef T
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

      // updateCcIDataReverse(double * s and double (*s)[3]...

#define T double
#define PACK_BUF pack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 52126
#include "updateCcIDataReverse.hpp"
#undef T
#undef PACK_BUF
#undef MPI_T
#undef UPDATE_TAG

      void buildCcPrcomm();

      // updateCcData(double * s and double (*s)[3]...

#define T double
#define PACK_BUF pack_buf_double
#define UNPACK_BUF unpack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 52123
#include "updateCcData.hpp"
#undef T
#undef PACK_BUF
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

      void calcCcGrad(double (*dpdx)[3],const double *p);

      void prolongCcData(double* phi_cv,const double* phi_cc);
      void prolongCcDataAndUpdateGuess(double* phi_cv,const double* phi_cc);
      void prolongCcData(double* phi_cv,const double (*x_cv)[3],const double* phi_cc,const double (*dphi_cc)[3]);
      void prolongCcDataAndUpdateGuess(double* phi_cv,const double (*x_cv)[3],const double* phi_cc,const double (*dphi_cc)[3]);
      void prolongCcData(double (*u_cv)[3],const double (*u_cc)[3]);
      void prolongCcDataAndUpdateGuess(double (*u_cv)[3],const double (*u_cc)[3],const double relax = 1.0);

      void restrictCcData(double* phi_cc,const double* phi_cv);
      void restrictExtrinsicCcData(double* phiV_cc,const double* phiV_cv);
      void restrictCcData(double (*u_cc)[3],const double (*u_cv)[3]);
      void restrictExtrinsicCcData(double (*uV_cc)[3],const double (*uV_cv)[3]);

      void init(const double* const vol_cv,const double (*const x_cv)[3],const double* const A_cv,const int* const cvocv_i,
          const int* const cvocv_v,const int8* const icv_global,const uint8* const rbicv_g,const uint8* const rbicv_g2,
          const int ncv,const int ncv_g,const int ncv_g2,const int8 ncv_global,const int icg,
          const double agglomeration_factor,const bool split_orphaned_colors,const double amg_coeff);

      void update(const double* const A_cv,const int* const cvocv_i,const int* const cvocv_v,const uint8* const rbicv_g,
          const uint8* const rbicv_g2,const int ncv,const int ncv_g,const int ncv_g2);

      inline void matvec(double * w, const double * A_cc, const double * phi) const { 
        for (int icc = 0; icc < ncc_a; ++icc) {
          const int coc_f = ccocc_i[icc];
          w[icc] = A_cc[coc_f]*phi[icc];
          for (int coc = coc_f+1; coc != ccocc_i[icc+1]; ++coc) {
            const int icc_nbr = ccocc_v[coc];
            w[icc] += A_cc[coc]*phi[icc_nbr];
          }
        }
      }

      inline double dot(const double* a, const double* b) const { 
        double my_sum = 0.0; 
        for (int icc = 0; icc < ncc_a; ++icc) 
          my_sum += a[icc]*b[icc];
        double sum;
        MPI_Allreduce(&my_sum,&sum,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
        return sum;
      }

      void calcCcResidual(double *res,const double* phi,const double *rhs);
      // perform nsmooth iterations unless zero > 0.0
      void smoothCcJacobi(double* phi,double *tmp,const double *rhs,const int nsmooth,const double relax,const double zero = -1.0);
      void smoothCcSgs(double* phi,double * res,const double *rhs,const int nsmooth,const double relax,const double zero = -1.0);
      void smoothCcGs(double* phi,double * res,const double *rhs,const int nsmooth,const double relax,const double zero = -1.0);
      void smoothCcCg(double *phi,double *res,double* v,double *p,const double *rhs,const int nsmooth,const double zero = -1.0);

      inline void generatePlaneRotation(double& dx, double &dy, double& cs, double& sn) const { 

        if ( dy == 0.0) { 
          cs = 1.0;
          sn = 0.0;
        } else if ( abs(dy) > abs(dx)) { 
          double tmp = dx/dy;
          sn = 1.0/ sqrt(1.0+tmp*tmp);
          cs = tmp*sn;
        } else { 
          double tmp = dy/dx;
          cs = 1.0/sqrt(1.0 + tmp*tmp);
          sn = tmp* cs;
        }

      }

      inline void applyPlaneRotation(double& dx, double& dy, double& cs, double& sn) const { 
        double tmp = cs * dx + sn * dy;
        dy         = -sn* dx + cs* dy;
        dx         = tmp;
      }

      void smoothCcGmres(double* phi,double *res,double* w,double* z,const double *rhs,const int nsmooth,const double zero = -1.0);

  };

  class AlgebraicMultigrid {
    public:
      int ncg;
      AlgebraicCoarseGrid* acg;
      
      AlgebraicMultigrid() {
        ncg = 0;
        acg = NULL;
      }

      ~AlgebraicMultigrid() {
        clear();
      }

      void clear() {
        DELETE(acg);
      }

      void init(double* vol_cv,double (* x_cv)[3],double* A_cv,int* cvocv_i,int* cvocv_v,int8* icv_global,
          uint8* rbicv_g,uint8* rbicv_g2,const int ncv,const int ncv_g,const int ncv_g2,const int8 ncv_global,
          const int ncg,const double agglomeration_factor,const bool split_orphaned_colors,const double amg_coeff);

      void update(double* A_cv);
      void solve(double* phi,const double* rhs,const double zero,const string cycle,const int maxcycle,const bool b_linear_prolongation,
          const int nsmooth,const int dnsmooth,const string smoother,const double relax,const int maxiter_coarsest,const string solver_coarsest,
          const bool b_verbose);

  };

}

#endif
