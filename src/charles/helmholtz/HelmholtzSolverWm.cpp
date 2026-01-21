
#include "HelmholtzSolver.hpp"

// convenience macro to loop over the boundary conditions...
#define FOR_BCZONE for(vector<HelmholtzBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)

void HelmholtzSolver::computeSlipLength() {

  const double dynAlfa_deltaR = getDoubleParam("DELTAR", 1.5);
  // XXX should only compute on subset of cells...
  double elapsed_time            = MPI_Wtime() ;
  double * rho_hat               = new double[ncv_g];
  double (*uhat)[3]              = new double[ncv_g][3];
  double (*rhou)[3]              = new double[ncv_g][3];
  double (*drhoudx)[3][3]        = new double[ncv_g][3][3];
  double (*drhoudx_hat)[3][3]    = new double[ncv_g][3][3];
  double *qijqij                 = new double[ncv_g];
  double *lijqij                 = new double[ncv_g];
  double *cv_fa_area             = new double[ncv_g];

  double (*lij_d)[3]             = new double[ncv_g][3];
  double (*lij_od)[3]            = new double[ncv_g][3];
  double (*tmp_v)[3]             = new double[ncv_g][3];
  double (*wall_motion)[3]       = new double[ncv_g][3];

  FOR_ICV FOR_I3 rhou[icv][i] = rho[icv] * u[icv][i];
  updateCvData(rhou, REPLACE_ROTATE_DATA);

  FOR_ICV_G {
    FOR_I3 {
      uhat[icv][i]         = 0.0;
      lij_d[icv][i]        = 0.0;
      lij_od[icv][i]       = 0.0;
      cv_fa_area[icv]      = 0.0;
      wall_motion[icv][i]  = 0.0;

      FOR_J3 {
        drhoudx[icv][i][j] = 0.0;
        drhoudx_hat[icv][i][j] = 0.0 ;
      }
    }//i

    qijqij[icv]  = 0.0 ;
    lijqij[icv]  = 0.0 ;
    rho_hat[icv] = 0.0;
  }

  // get the favre-averaged test-filtered velocity field ..
  filterCvR1_mod(rho_hat, rho) ;

  //StaticSolver::calcCvGrad(drhoudx, rhou);
  calcCvGrad(drhoudx,rhou);
  filterCvR2_mod(uhat,rhou);
  FOR_ICV FOR_I3 uhat[icv][i] /= rho_hat[icv] ;

  updateCvData(drhoudx, REPLACE_ROTATE_DATA) ;
  filterCvR3_mod(drhoudx_hat, drhoudx);

  //dumpVectorRange(uhat, nno, "UHAT") ;

  int * cv_flag = new int[ncv];
  FOR_ICV cv_flag[icv] = -1;

  FOR_ICV_G {
    FOR_I3 lij_d[icv][i] = rho[icv] * u[icv][i] * u[icv][i];
    lij_od[icv][0]       = rho[icv] * u[icv][0] * u[icv][2];
    lij_od[icv][1]       = rho[icv] * u[icv][2] * u[icv][1];
    lij_od[icv][2]       = rho[icv] * u[icv][1] * u[icv][0];
  }


  FOR_BCZONE {
    Param * param = getParam((*it)->getName());
    if ( param->getString(0) == "WM_SLIP_WALL") {
      SlipWallModelHBc * slipWall = dynamic_cast<SlipWallModelHBc*>(*it);
      for (int ibf=0;ibf<slipWall->zone_ptr->nbf;++ibf){

        //slipWall->cdel_w[ibf] = 0.0;  //DAP: this breaks relaxation if uncommented!!!!!!

        //==================================================
        // dynamic procedure after initializations ..
        //==================================================

        //===
        // compute lij ..
        //===
        const int icv = slipWall->zone_ptr->cvobf[ibf];
        cv_flag[icv] = 0;
        FOR_I3 wall_motion[icv][i] = slipWall->u_bc[ibf][i];
        cv_fa_area[icv]      += slipWall->zone_ptr->area_bf[ibf] ;
      }
    }
  }

  filterCvR2_mod( tmp_v, lij_d) ;
  FOR_ICV {
    FOR_I3 lij_d[icv][i] = tmp_v[icv][i];
  }

  filterCvR2_mod( tmp_v, lij_od);
  FOR_ICV {
    FOR_I3 lij_od[icv][i] = tmp_v[icv][i];
  }

  // \hat{\bu_i\bu_j} - \bu_i\bu_j
  FOR_ICV {
    FOR_I3 lij_d[icv][i] -= rho[icv] * u[icv][i] * u[icv][i] ;
    lij_od[icv][0]       -= rho[icv] * u[icv][0] * u[icv][2] ;
    lij_od[icv][1]       -= rho[icv] * u[icv][2] * u[icv][1] ;
    lij_od[icv][2]       -= rho[icv] * u[icv][1] * u[icv][0] ;

    // TODO DAP confirm moving slip wm terms with STB
    // and the wall motion terms .. approximations for small density
    // perturbations...
    FOR_I3 lij_d[icv][i] -= ( rho[icv] * wall_motion[icv][i] * (uhat[icv][i]-u[icv][i]) +
        rho[icv] * wall_motion[icv][i] * (uhat[icv][i]-u[icv][i])  );
    lij_od[icv][0]       -= ( rho[icv] * wall_motion[icv][0] * (uhat[icv][2]-u[icv][2]) +
        rho[icv] * wall_motion[icv][2] * (uhat[icv][0]-u[icv][0])  );
    lij_od[icv][1]       -= ( rho[icv] * wall_motion[icv][2] * (uhat[icv][1]-u[icv][1]) +
        rho[icv] * wall_motion[icv][1] * (uhat[icv][2]-u[icv][2])  );
    lij_od[icv][2]       -= ( rho[icv] * wall_motion[icv][1] * (uhat[icv][0]-u[icv][0]) +
        rho[icv] * wall_motion[icv][0] * (uhat[icv][1]-u[icv][1])  );

  }


  double * qijqij_bf = new double[nbf];
  double * lijqij_bf = new double[nbf];
  int * bf_flag      = new int[nbf];

  for (int ibf = 0; ibf < nbf; ++ibf) { 

    qijqij_bf[ibf] = 0.0;
    lijqij_bf[ibf] = 0.0;
    bf_flag[ibf]   = 0;

  }


  FOR_BCZONE {
    Param * param = getParam((*it)->getName());
    if ( param->getString(0) == "WM_SLIP_WALL") {

      SlipWallModelHBc * slipWall = dynamic_cast<SlipWallModelHBc*>(*it);

      for (int ibf = 0; ibf < slipWall->zone_ptr->nbf; ++ibf) {

        const int icv = slipWall->zone_ptr->cvobf[ibf];

        //===========================================
        // compute qij -----
        //===========================================
        double drhoudn_hat[3], drhoudn[3];
        double qij_d[3], qij_od[3];

        double unit_n[3] = {slipWall->zone_ptr->n_bf[ibf][0],slipWall->zone_ptr->n_bf[ibf][1],slipWall->zone_ptr->n_bf[ibf][2]};
        double mag_n     = MAG(unit_n);
        FOR_I3 unit_n[i] = unit_n[i]/mag_n;
        FOR_I3 {
          drhoudn_hat[i]   = -DOT_PRODUCT(drhoudx_hat[icv][i], unit_n);
          drhoudn[i]       = -DOT_PRODUCT(drhoudx[icv][i]    , unit_n);
          qij_d[i]         = dynAlfa_deltaR* dynAlfa_deltaR* drhoudn_hat[i] * drhoudn_hat[i] / rho_hat[icv] -
            drhoudn[i]     * drhoudn[i]     / rho[icv]        ;
        }

        qij_od[0]  = dynAlfa_deltaR* dynAlfa_deltaR* drhoudn_hat[0] * drhoudn_hat[2] / rho_hat[icv] -
          drhoudn[0]     * drhoudn[2]     / rho[icv]        ;
        qij_od[1]  = dynAlfa_deltaR* dynAlfa_deltaR* drhoudn_hat[2] * drhoudn_hat[1] / rho_hat[icv] -
          drhoudn[2]     * drhoudn[1]     / rho[icv]        ;
        qij_od[2]  = dynAlfa_deltaR* dynAlfa_deltaR* drhoudn_hat[1] * drhoudn_hat[0] / rho_hat[icv] -
          drhoudn[1]     * drhoudn[0]     / rho[icv]        ;


        FOR_I3 {
          qijqij[icv] += slipWall->zone_ptr->area_bf[ibf] * (qij_d[i]      * qij_d[i]      + 2.0 * qij_od[i]      * qij_od[i] );
          lijqij[icv] += slipWall->zone_ptr->area_bf[ibf] * (lij_d[icv][i] * qij_d[i]      + 2.0 * lij_od[icv][i] * qij_od[i] );
        }
      }
    }
  }


  FOR_ICV {
    if (cv_flag[icv] == 0) {
      assert ( cv_fa_area[icv] > 0.0 ) ;
      qijqij[icv] /= cv_fa_area[icv];
      lijqij[icv] /= cv_fa_area[icv];

      // clip prior to smooth..
      lijqij[icv]  = max(0.0,lijqij[icv]);
    }
  }
  delete[] cv_flag;

  //    if (step%check_interval==0){
  //      MiscUtils::dumpRange(qijqij,ncv, "qijqij-cv");
  //      MiscUtils::dumpRange(lijqij,ncv, "lijqij-cv");
  //    }

  // put the data back on the bfs .. 

  FOR_BCZONE {
    Param * param = getParam((*it)->getName());
    if ( param->getString(0) == "WM_SLIP_WALL") {

      SlipWallModelHBc * slipWall = dynamic_cast<SlipWallModelHBc*>(*it);

      for (int ibf = 0; ibf < slipWall->zone_ptr->nbf; ++ibf) {

        const int icv = slipWall->zone_ptr->cvobf[ibf];
        qijqij_bf[ibf+slipWall->zone_ptr->ibf_f] = qijqij[icv];
        lijqij_bf[ibf+slipWall->zone_ptr->ibf_f] = lijqij[icv];
        bf_flag[ibf+slipWall->zone_ptr->ibf_f]  = 1;
      }
    }
  }

  BfFilter bf_filter(this,bf_flag);

  double * tmp_bf1 = new double[nbf];
  double * tmp_bf2 = new double[nbf];

  const int wm_smooth_iter = getIntParam("WM_SMOOTH_ITER", 50);

  for (int iter = 0; iter < wm_smooth_iter; ++iter) { 

    bf_filter.doFilter(tmp_bf1,qijqij_bf);
    bf_filter.doFilter(tmp_bf2,lijqij_bf);

    for (int ibf = 0; ibf < nbf; ++ibf) { 
      qijqij_bf[ibf] = tmp_bf1[ibf];
      lijqij_bf[ibf] = tmp_bf2[ibf];
    }
  }

  delete[] tmp_bf1;
  delete[] tmp_bf2;

  //
  // average the different vectors.. note that for terms involving, qij
  // only averaging involving the boundary can be used..
  //
#ifdef NO_AVG
  {

    double (*tmp_v_no)[3]     = new double[nno][3];
    double (*tmp_v_no_hat)[3] = new double[nno][3];
    //
    // local averaging if fnacy avg isnt available but this can only be
    // done using wall nodes for any terms using qij
    //

    if (step%check_interval==0) COUT1("> Using WALL_MODEL_AVG = LOCAL") ;

    double elapsed2 = MPI_Wtime() ;
    FOR_INO {
      FOR_I3 tmp_v_no[ino][i] = 0.0;
    }

    int * no_flag = new int[nno];
    FOR_INO no_flag[ino] = -1;


    FOR_BCZONE {
      Param * param = getParam((*it)->getName());
      if ( param->getString(0) == "WM_SLIP_WALL") {
        SlipWallModelHBc * slipWall = dynamic_cast<SlipWallModelHBc*>(*it);
        for (int ibf=0;ibf<slipWall->zone_ptr->nbf;++ibf){
          for (int nob = noobf_i[slipWall->zone_ptr->ibf_f+ibf]; nob != noobf_i[slipWall->zone_ptr->ibf_f+ibf+1]; ++nob) {
            const int ino = noobf_v[nob]; assert((ino >= 0)&&(ino < nno));
            no_flag[ino] = 1;
          }
        }
      }
    }

    updateNoData(no_flag, MAX_NO_PERIODIC_DATA);

    FOR_BCZONE {
      Param * param = getParam((*it)->getName());
      if ( param->getString(0) == "WM_SLIP_WALL") {
        SlipWallModelHBc * slipWall = dynamic_cast<SlipWallModelHBc*>(*it);
        for (int ibf=0;ibf<slipWall->zone_ptr->nbf;++ibf){
          const int nof_f = noobf_i[slipWall->zone_ptr->ibf_f+ibf];
          const int nof_l = noobf_i[slipWall->zone_ptr->ibf_f+ibf+1]-1;
          //const int nnof  = nof_l - nof_f +1;

          // store the data into tmp bufs, weighted
          // by the sub-tet volume
          int ino1 = noobf_v[nof_l] ;

          // loop on nodes-of-face in order that face normal is outward pointing
          // with nodes in rh-format...

          for (int nof = nof_f; nof <= nof_l ; ++nof) {

            // copy down previous ino1 stuff to ino0...

            int ino0 = ino1;
            ino1 = noobf_v[nof];

            double dx_no0[3], dx_no1[3];
            for (int i = 0; i < 3; i++) {
              dx_no0[i]  = x_no[ino0][i] - slipWall->zone_ptr->x_bf[ibf][i];
              dx_no1[i]  = x_no[ino1][i] - slipWall->zone_ptr->x_bf[ibf][i];
            }

            // we now have a subtri ifa-ino0-ino1 ..
            double st_normal[3];
            st_normal[0] = 0.25* ( dx_no0[1]*dx_no1[2] - dx_no0[2]*dx_no1[1] ) ;
            st_normal[1] = 0.25* ( dx_no0[2]*dx_no1[0] - dx_no0[0]*dx_no1[2] );
            st_normal[2] = 0.25* ( dx_no0[0]*dx_no1[1] - dx_no0[1]*dx_no1[0]) ;

            const double st_area = sqrt(DOT_PRODUCT(st_normal, st_normal));

            tmp_v_no[ino0][0] += st_area * qijqij[slipWall->zone_ptr->cvobf[ibf]];
            tmp_v_no[ino0][1] += st_area * lijqij[slipWall->zone_ptr->cvobf[ibf]];

            tmp_v_no[ino1][0] += st_area * qijqij[slipWall->zone_ptr->cvobf[ibf]];
            tmp_v_no[ino1][1] += st_area * lijqij[slipWall->zone_ptr->cvobf[ibf]];

            tmp_v_no[ino0][2] += st_area;
            tmp_v_no[ino1][2] += st_area;
          }
        }
      }
    }
    updateNoData( tmp_v_no,   ADD_NO_PERIODIC_DATA) ;

    FOR_INO {
      if (no_flag[ino] >= 0) {
        tmp_v_no[ino][0] /= tmp_v_no[ino][2];
        tmp_v_no[ino][1] /= tmp_v_no[ino][2];
      }
    }

    for (int iter = 0 ; iter < wm_smooth_iter ; ++iter ) {

      FOR_INO {
        FOR_I3 tmp_v_no_hat[ino][i] = 0.0 ;
      }


      // clip & smooth ..
      FOR_INO {
        tmp_v_no[ino][1] = max(tmp_v_no[ino][1],0.0);
      }

      FOR_BCZONE {
        Param * param = getParam((*it)->getName());
        if ( param->getString(0) == "WM_SLIP_WALL") {
          SlipWallModelHBc * slipWall = dynamic_cast<SlipWallModelHBc*>(*it);
          for (int ibf=0;ibf<slipWall->zone_ptr->nbf;++ibf){
            const int nof_f = noobf_i[slipWall->zone_ptr->ibf_f+ibf];
            const int nof_l = noobf_i[slipWall->zone_ptr->ibf_f+ibf+1]-1;
            const int nnof  = nof_l - nof_f +1;

            double qijqij_bf = 0.0 ;
            double lijqij_bf = 0.0 ;

            for (int nof = nof_f ; nof <= nof_l ; ++nof) {
              const int ino = noobf_v[nof];
              assert ( no_flag[ino] >= 0) ;
              qijqij_bf += tmp_v_no[ino][0];
              lijqij_bf += tmp_v_no[ino][1];
            }
            qijqij_bf /= double(nnof) ;
            lijqij_bf /= double(nnof) ;

            // store the data into tmp bufs, weighted
            // by the sub-tet volume
            int ino1 = noobf_v[nof_l] ;

            // loop on nodes-of-face in order that face normal is outward pointing
            // with nodes in rh-format...

            for (int nof = nof_f; nof <= nof_l ; ++nof) {

              // copy down previous ino1 stuff to ino0...

              int ino0 = ino1;
              ino1 = noobf_v[nof];

              double dx_no0[3], dx_no1[3];
              for (int i = 0; i < 3; i++) {
                dx_no0[i]  = x_no[ino0][i] - slipWall->zone_ptr->x_bf[ibf][i];
                dx_no1[i]  = x_no[ino1][i] - slipWall->zone_ptr->x_bf[ibf][i];
              }

              // we now have a subtri ibf-ino0-ino1 ..
              double st_normal[3];
              st_normal[0] = 0.25* ( dx_no0[1]*dx_no1[2] - dx_no0[2]*dx_no1[1] ) ;
              st_normal[1] = 0.25* ( dx_no0[2]*dx_no1[0] - dx_no0[0]*dx_no1[2] );
              st_normal[2] = 0.25* ( dx_no0[0]*dx_no1[1] - dx_no0[1]*dx_no1[0]) ;

              double st_area = sqrt( DOT_PRODUCT(st_normal,st_normal)) ;

              tmp_v_no_hat[ino0][0] += st_area * ( tmp_v_no[ino0][0] + tmp_v_no[ino1][0] + qijqij_bf)/3.0;
              tmp_v_no_hat[ino1][0] += st_area * ( tmp_v_no[ino0][0] + tmp_v_no[ino1][0] + qijqij_bf)/3.0;

              tmp_v_no_hat[ino0][1] += st_area * ( tmp_v_no[ino0][1] + tmp_v_no[ino1][1] + lijqij_bf)/3.0 ;
              tmp_v_no_hat[ino1][1] += st_area * ( tmp_v_no[ino0][1] + tmp_v_no[ino1][1] + lijqij_bf)/3.0 ;

              tmp_v_no_hat[ino0][2] += st_area ;
              tmp_v_no_hat[ino1][2] += st_area ;
            }//nof
          }
        }
      }


      updateNoData( tmp_v_no_hat,   ADD_NO_PERIODIC_DATA) ;

      // and replace the data
      FOR_INO {
        if ( no_flag[ino] >= 0 ) {
          assert ( tmp_v_no_hat[ino][2] > 0.0)  ;
          tmp_v_no[ino][0] = tmp_v_no_hat[ino][0] / tmp_v_no_hat[ino][2];
          tmp_v_no[ino][1] = tmp_v_no_hat[ino][1] / tmp_v_no_hat[ino][2];
        }
      }//ino
    }//iter

    // push back to faces
    FOR_ICV {
      qijqij[icv] = 0.0;
      lijqij[icv] = 0.0;
    }

    FOR_BCZONE {
      Param * param = getParam((*it)->getName());
      if ( param->getString(0) == "WM_SLIP_WALL") {
        SlipWallModelHBc * slipWall = dynamic_cast<SlipWallModelHBc*>(*it);
        for (int ibf=0;ibf<slipWall->zone_ptr->nbf;++ibf){
          const int nof_f = noobf_i[slipWall->zone_ptr->ibf_f+ibf];
          const int nof_l = noobf_i[slipWall->zone_ptr->ibf_f+ibf+1]-1;
          const int nnof  = nof_l - nof_f +1;
          //double x_bf[3] = { 0.0,0.0,0.0} ;

          double qijqij_bf = 0.0 ;
          double lijqij_bf = 0.0 ;

          for (int nof = nof_f ; nof <= nof_l ; ++nof) {
            const int ino = noobf_v[nof];
            assert ( no_flag[ino] >= 0) ;
            qijqij_bf += tmp_v_no[ino][0];
            lijqij_bf += tmp_v_no[ino][1];
          }
          qijqij_bf /= double(nnof) ;
          lijqij_bf /= double(nnof) ;

          const int icv = slipWall->zone_ptr->cvobf[ibf];
          qijqij[icv] += qijqij_bf * slipWall->zone_ptr->area_bf[ibf] / cv_fa_area[icv];
          lijqij[icv] += lijqij_bf * slipWall->zone_ptr->area_bf[ibf] / cv_fa_area[icv];

        }
      }
    }

    delete[] no_flag;

    elapsed2 = MPI_Wtime() - elapsed2 ;
    if ( mpi_rank == 0 && step%check_interval==0)
      cout << " dynamic alpha local avg [secs]: " << elapsed2 << endl ;

    delete[] tmp_v_no;
    delete[] tmp_v_no_hat;

  }//local averaging ...
#endif

  double cdel_w_relax = getDoubleParam("SLIP_LENGTH_RELAX",1.0);
  if (cdel_w_relax < 0.0)
    cdel_w_relax = 0.0;
  else if (cdel_w_relax > 1.0)
    cdel_w_relax = 1.0;

  //
  // now compute the dynamic coefficients and
  // set the viscosity ..
  //
  int my_clip_count[2] = {0,0} ;
  int clip_count[2] ;

  double neg_slip_length = 0.0;

  FOR_BCZONE {
    Param * param = getParam((*it)->getName());
    if ( param->getString(0) == "WM_SLIP_WALL") {
      SlipWallModelHBc * slipWall = dynamic_cast<SlipWallModelHBc*>(*it);
      for (int ibf=0;ibf<slipWall->zone_ptr->nbf;++ibf){

        //const int icv = slipWall->zone_ptr->cvobf[ibf];
        //double yy ;

        ++my_clip_count[1];

        // C_W\Delta^2 ..
        //double xx = lijqij[icv] / (qijqij[icv] + 1.0e-10);
        double xx = lijqij_bf[ibf+slipWall->zone_ptr->ibf_f] / (qijqij_bf[ibf+slipWall->zone_ptr->ibf_f] + 1.0e-10);

        // now some clipping
        if ( xx >= 0.0 )
          slipWall->cdel_w[ibf] = cdel_w_relax * sqrt( xx) + (1.0 - cdel_w_relax) * slipWall->cdel_w[ibf]; //relaxation added
        else {
          ++my_clip_count[0] ;
          neg_slip_length = max(neg_slip_length,sqrt(-xx));
          slipWall->cdel_w[ibf] = (1.0 - cdel_w_relax) * slipWall->cdel_w[ibf]; //relaxtion added
        }

      }
    }
  }//FOR_FAZONE_BOUNDARY

  if (step%check_interval==0){
    MPI_Reduce(my_clip_count, clip_count, 2, MPI_INT, MPI_SUM, 0, mpi_comm);
    if (( mpi_rank == 0 ) && (step % check_interval == 0))
      cout << " ... wm clipping : " << clip_count[0] << " of " << clip_count[1] << " ( "
        << double(clip_count[0])/double(clip_count[1])*100.0 << " %)" << endl ;

    MiscUtils::dumpRange(&neg_slip_length,1,"neg slip length");
  }


  delete[] qijqij_bf;
  delete[] lijqij_bf;
  delete[] bf_flag;

  //
  // cleanup
  //
  delete[] rho_hat ;
  delete[] uhat;
  delete[] drhoudx ;
  delete[] drhoudx_hat ;
  delete[] qijqij;
  delete[] lijqij;
  delete[] cv_fa_area;
  delete[] wall_motion;

  delete[] lij_d ;
  delete[] lij_od ;

  delete[] tmp_v;

  elapsed_time = MPI_Wtime() - elapsed_time ;
  if (( mpi_rank == 0 ) && (step % check_interval == 0))
    cout << " dynamic alpha calc [secs] : " << elapsed_time << endl ;


  // for debugging, you can manually specify a constant cdel...
  if (checkParam("CDEL")) {
    double cdel_constant = getDoubleParam("CDEL") ;
    FOR_BCZONE {
      Param * param = getParam((*it)->getName());
      if ( param->getString(0) == "WM_SLIP_WALL") {
        SlipWallModelHBc * slipWall = dynamic_cast<SlipWallModelHBc*>(*it);
        for (int ibf=0;ibf<slipWall->zone_ptr->nbf;++ibf){
          const int icv = slipWall->zone_ptr->cvobf[ibf];
          slipWall->cdel_w[ibf] = cdel_constant*0.5*pow(vol_cv[icv],1.0/3.0) ;
        }
      }
    }
  }

  //Limit Slip Length
  if (checkParam("CDEL_MAX")) {
    double cdel_constant = getDoubleParam("CDEL_MAX") ;
    FOR_BCZONE {
      Param * param = getParam((*it)->getName());
      if ( param->getString(0) == "WM_SLIP_WALL") {
        SlipWallModelHBc * slipWall = dynamic_cast<SlipWallModelHBc*>(*it);
        for (int ibf=0;ibf<slipWall->zone_ptr->nbf;++ibf){
          if (slipWall->cdel_w[ibf] > cdel_constant)
            slipWall->cdel_w[ibf] = cdel_constant ;
        }
      }
    }
  }

  //    Use QUERY_BC instead...
  //    FOR_BCZONE {
  //      Param * param = getParam((*it)->getName());
  //      if ( param->getString(0) == "WM_SLIP_WALL") {
  //        SlipWallModelHBc * slipWall = dynamic_cast<SlipWallModelHBc*>(*it);
  //        if (step%check_interval==0){
  //          dumpRange(slipWall->cdel_w,slipWall->zone_ptr->nbf,slipWall->getName()+".cdel_w") ;
  //        }
  //      }
  //    }


  delete[] rhou;
}


void HelmholtzSolver::filterCvR1_mod(double *phif, const double *phi ) {

  double (*tmp) = new double[nfa];
  double (*tmp_bf) = new double[nbf];

  for (int ifa=0; ifa<nfa; ++ifa) {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    tmp[ifa] = 0.5 * ( phi[icv0] + phi[icv1] ) ;
  }

  FOR_IBF {
    const int icv0 = cvobf[ibf];
    tmp_bf[ibf] = phi[icv0];
  }

  FOR_ICV {

    phif[icv] = 0.0;
    double area_total = 0.0;

    for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
      const int ibf = bfocv_v[boc];
      phif[icv] += area_bf[ibf]*tmp_bf[ibf];
      area_total += area_bf[ibf];
    }

    for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
      const int ifa     = faocv_v[foc];
      const double area = MAG(n_fa[ifa]);
      phif[icv] += area*tmp[ifa];
      area_total += area;
    }

    phif[icv] /= area_total;

  }

  delete[] tmp;
  delete[] tmp_bf;

}

void HelmholtzSolver::filterCvR2_mod(double (*phif)[3], const double (*phi)[3] ) {

  double (*tmp)[3] = new double[nfa][3];
  double (*tmp_bf)[3] = new double[nbf][3];

  for (int ifa=0; ifa<nfa; ++ifa) {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    FOR_I3 tmp[ifa][i] = 0.5 * ( phi[icv0][i] + phi[icv1][i] );

  }

  FOR_IBF {
    const int icv0 = cvobf[ibf];
    FOR_I3 tmp_bf[ibf][i] = phi[icv0][i];
  }

  FOR_ICV {
    FOR_I3 phif[icv][i] = 0.0;
    double area_total = 0.0;

    for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
      const int ibf = bfocv_v[boc];
      area_total += area_bf[ibf];
      FOR_I3 phif[icv][i] += area_bf[ibf]*tmp_bf[ibf][i];
    }

    for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
      const int ifa = faocv_v[foc];
      double area   = MAG(n_fa[ifa]);
      area_total   += area;
      FOR_I3 phif[icv][i] += area*tmp[ifa][i];
    }

    FOR_I3 phif[icv][i] /= area_total;
  }

  delete[] tmp;
  delete[] tmp_bf;


}

void HelmholtzSolver::filterCvR3_mod(double (*phif)[3][3], const double (*phi)[3][3] ) {

  double (*tmp)[3][3] = new double[nfa][3][3];
  double (*tmp_bf)[3][3] = new double[nbf][3][3];

  for (int ifa=0; ifa<nfa; ++ifa) {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    FOR_I3 FOR_J3 tmp[ifa][i][j] = 0.5 * ( phi[icv0][i][j] + phi[icv1][i][j] ) ;
  }

  FOR_IBF {
    const int icv0 = cvobf[ibf];
    FOR_I3 FOR_J3 tmp_bf[ibf][i][j] = phi[icv0][i][j];
  }

  FOR_ICV {
    FOR_I3 FOR_J3 phif[icv][i][j] = 0.0;
    double area_total = 0.0;

    for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
      const int ibf = bfocv_v[boc];
      area_total += area_bf[ibf];
      FOR_I3 FOR_J3 phif[icv][i][j] += area_bf[ibf]*tmp_bf[ibf][i][j];
    }

    for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
      const int ifa = faocv_v[foc];
      const double area = MAG(n_fa[ifa]);
      area_total += area;
      FOR_I3 FOR_J3 phif[icv][i][j] += area*tmp[ifa][i][j];
    }

    FOR_I3 FOR_J3 phif[icv][i][j] /= area_total;

  }

  delete[] tmp;
  delete[] tmp_bf;

}

#undef FOR_BCZONE
