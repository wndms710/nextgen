
#include "HelmholtzSolver.hpp"

void HelmholtzSolver::computeSgsNone() {

  FOR_ICV_G mu_sgs[icv] = 0.0;
}

void HelmholtzSolver::computeSgsVreman() {

  //
  // constant coeff vreman sgs ..
  //

  const double vreman_coeff = getDoubleParam("VREMAN_COEFF", 0.07);
  double (*duidxj)[3][3] = new double[ncv][3][3];
  //StaticSolver::calcCvGrad(duidxj,u); 
  calcCvGrad(duidxj,u); 

  FOR_ICV {

    double dx2 = pow( vol_cv[icv], 2.0/3.0 );
    double dy2 = dx2;
    double dz2 = dx2;

    double alpha11 = duidxj[icv][0][0];
    double alpha22 = duidxj[icv][1][1];
    double alpha33 = duidxj[icv][2][2];

    double alpha12 = duidxj[icv][0][1];
    double alpha13 = duidxj[icv][0][2];
    double alpha23 = duidxj[icv][1][2];

    double alpha21 = duidxj[icv][1][0];
    double alpha31 = duidxj[icv][2][0];
    double alpha32 = duidxj[icv][2][1];

    double beta11  = dx2*alpha11*alpha11+dy2*alpha21*alpha21+dz2*alpha31*alpha31;
    double beta12  = dx2*alpha11*alpha12+dy2*alpha21*alpha22+dz2*alpha31*alpha32;
    double beta13  = dx2*alpha11*alpha13+dy2*alpha21*alpha23+dz2*alpha31*alpha33;
    double beta22  = dx2*alpha12*alpha12+dy2*alpha22*alpha22+dz2*alpha32*alpha32;
    double beta23  = dx2*alpha12*alpha13+dy2*alpha22*alpha23+dz2*alpha32*alpha33;
    double beta33  = dx2*alpha13*alpha13+dy2*alpha23*alpha23+dz2*alpha33*alpha33;

    double B       = beta11*beta22-beta12*beta12+beta11*beta33-beta13*beta13+beta22*beta33-beta23*beta23;
    B              = (B + fabs(B))*0.5;
    double alal    =
      alpha11*alpha11+alpha22*alpha22+alpha33*alpha33 +
      alpha12*alpha12+alpha13*alpha13+alpha23*alpha23 +
      alpha21*alpha21+alpha31*alpha31+alpha32*alpha32;

    double s_mag = sqrt(B/(alal+1.0E-20)); // includes lengthscale squared too

    // mu_sgs...
    mu_sgs[icv] = rho[icv]*vreman_coeff*s_mag;

    // clip at the constant smagorinsky level with c2 = 0.5^2...
    mu_sgs[icv] = min( mu_sgs[icv], rho[icv]*0.05*dx2*sqrt(2.0*alal) );
  }

  updateCvData(mu_sgs); //set ghost cvs...

  // HACK! to match vida vor bug...
  //for (int icv = ncv; icv < ncv_g;++icv) 
  // mu_sgs[icv] = 0.0;

  delete[] duidxj;

}

void HelmholtzSolver::computeSgsSigma() {

  const double sigma_coeff = getDoubleParam("SIGMA_COEFF", 1.5);
  double (*duidxj)[3][3] = new double[ncv][3][3];
  //StaticSolver::calcCvGrad(duidxj,u); 
  calcCvGrad(duidxj,u); 

  for (int icv = 0; icv < ncv; ++icv) { 

    double A_tmp[9];
    for (int i = 0; i < 3; ++i) { 
      for (int j = 0; j < 3; ++j) 
        A_tmp[(j)*3+i] = duidxj[icv][i][j]; // store in col major format..
    }

    double U[9], V[9], sigma[3], max_err;
    calcSvd33(U,V,sigma,A_tmp,max_err);

    //
    // sort the sigma to be descending order .. 
    // there's only 3 values, so we can bubble sort pretty easily .. 
    // 

    if ( sigma[0] < sigma[1])  
      swap(sigma[0],sigma[1]);

    if ( sigma[1] < sigma[2]) 
      swap(sigma[1], sigma[2]);

    if ( sigma[0] < sigma[1]) 
      swap(sigma[0],sigma[1]);

    assert( (sigma[0] >= sigma[1]) && (sigma[0] >= sigma[2]));
    assert(  sigma[1] >= sigma[2]);

    double D_sigma = 0.0;
    if ( sigma[0] > 1.0e-12) { // protect against zero gradient case .. 

      D_sigma = sigma[2]*(sigma[0] - sigma[1])*(sigma[1] - sigma[2])/(sigma[0]*sigma[0]);

    }

    const double dx2 = pow(vol_cv[icv],2.0/3.0);
    mu_sgs[icv]      = sigma_coeff*sigma_coeff*dx2*D_sigma;

  }


  updateCvData(mu_sgs);

  delete[] duidxj;

} 

void HelmholtzSolver::computeSgsDsm() {

  //TODO follow the ideal gas procedure

  const int nfilter = 5;

  //
  // dynamic smagorinsky sgs ...
  //

  double (*dudx)[3][3] = new double[ncv_g][3][3];
  //StaticSolver::calcCvGrad(dudx,u); 
  calcCvGrad(dudx,u); 

  double * s_mag  = new double[ncv_g];
  double * lijmij = new double[ncv_g];
  double * mijmij = new double[ncv_g];

  updateCvData(rho);
  updateCvData(u);
  updateCvData(dudx);
  calcDSMStuff(s_mag,lijmij,mijmij,rho,u,dudx);
  delete[] dudx;

  // filter nfilterx2 times...

  for (int iter = 0; iter < nfilter; ++iter) {
    updateCvData(lijmij);
    filterCvR1_mod(mu_sgs,lijmij);
    updateCvData(mu_sgs);
    filterCvR1_mod(lijmij,mu_sgs);
    updateCvData(mijmij);
    filterCvR1_mod(mu_sgs,mijmij);
    updateCvData(mu_sgs);
    filterCvR1_mod(mijmij,mu_sgs);
  }

  // now calc mu_sgs, clipping to zero...

  FOR_ICV mu_sgs[icv] = max(0.0,lijmij[icv]/(mijmij[icv]+1.0E-20))*s_mag[icv];

  delete[] s_mag;
  delete[] lijmij;
  delete[] mijmij;

  updateCvData(mu_sgs); //set ghost cvs...

}

void HelmholtzSolver::calcDSMStuff(double* s_mag, double* lijmij, double* mijmij,
    const double* rho, const double (*ui)[3],
    const double (*duidxj)[3][3]) {

  // rho, u and dudx must be updated to ghosts...

  double filter_ratio_squared = 4.0;

  // memory...

  double (*sij_d)[3]     = new double[ncv_g][3];
  double (*sij_od)[3]    = new double[ncv_g][3];
  double (*ui_hat)[3]    = new double[ncv_g][3];
  double (*mij_d)[3]     = new double[ncv][3];
  double (*mij_od)[3]    = new double[ncv][3];
  double* rho_hat        = new double[ncv];
  double (*rhoui_hat)[3] = new double[ncv][3];

  // compute Sij, Smag...

  for (int icv = 0; icv < ncv_g; icv++) {

    // diagonal...
    sij_d[icv][0] = duidxj[icv][0][0];
    sij_d[icv][1] = duidxj[icv][1][1];
    sij_d[icv][2] = duidxj[icv][2][2];

    // remove trace...
    double tmp = (sij_d[icv][0] + sij_d[icv][1] + sij_d[icv][2])/3.0;
    sij_d[icv][0] -= tmp;
    sij_d[icv][1] -= tmp;
    sij_d[icv][2] -= tmp;

    // off-diagonal...
    sij_od[icv][0] = (duidxj[icv][2][1] + duidxj[icv][1][2])/2.0;
    sij_od[icv][1] = (duidxj[icv][0][2] + duidxj[icv][2][0])/2.0;
    sij_od[icv][2] = (duidxj[icv][1][0] + duidxj[icv][0][1])/2.0;

    // Smag = sqrt(2 Sij Sij)...
    s_mag[icv] = sqrt( 2.0*( sij_d[icv][0]*sij_d[icv][0] + 
          sij_d[icv][1]*sij_d[icv][1] + 
          sij_d[icv][2]*sij_d[icv][2] ) +
        4.0*( sij_od[icv][0]*sij_od[icv][0] + 
          sij_od[icv][1]*sij_od[icv][1] + 
          sij_od[icv][2]*sij_od[icv][2] ) );

    // multiply Sij by its magnitude and density in preparation for filtering...
    for (int i = 0; i < 3; i++) {
      sij_d[icv][i]  *= rho[icv]*s_mag[icv];
      sij_od[icv][i] *= rho[icv]*s_mag[icv];
    }

  }

  // filter the velocity, density, momentum, and rho|S|Sij...
  double (*rhoui)[3]     = new double[ncv_g][3];
  FOR_ICV_G FOR_I3 rhoui[icv][i] = rho[icv]*u[icv][i];
  filterCvR2_mod(rhoui_hat,rhoui);
  delete[] rhoui;
  filterCvR1_mod(rho_hat,rho);  
  filterCvR2_mod(ui_hat,ui);  
  filterCvR2_mod(mij_d,sij_d);
  filterCvR2_mod(mij_od,sij_od);

  // compute the filtered strain rate...

  updateCvData(ui_hat);
  double (*duidxj_hat)[3][3] = new double[ncv][3][3];
  //StaticSolver::calcCvGrad(duidxj_hat,ui_hat);
  calcCvGrad(duidxj_hat,ui_hat);

  // complete Mij = 2 ( hat(|S| S_ij) - (alpha^2) |S_hat| S_hat_ij )...

  for (int icv = 0; icv < ncv; icv++) {

    // diagonal...
    double sij_hat_d[3];
    sij_hat_d[0] = duidxj_hat[icv][0][0]; // filtered gradient
    sij_hat_d[1] = duidxj_hat[icv][1][1];
    sij_hat_d[2] = duidxj_hat[icv][2][2];

    // remove trace...
    double tmp = (sij_hat_d[0] + sij_hat_d[1] + sij_hat_d[2])/3.0;
    sij_hat_d[0] -= tmp;
    sij_hat_d[1] -= tmp;
    sij_hat_d[2] -= tmp;

    // off-diagonal...
    double sij_hat_od[3];
    sij_hat_od[0] = (duidxj_hat[icv][2][1] + duidxj_hat[icv][1][2])/2.0;
    sij_hat_od[1] = (duidxj_hat[icv][0][2] + duidxj_hat[icv][2][0])/2.0;
    sij_hat_od[2] = (duidxj_hat[icv][1][0] + duidxj_hat[icv][0][1])/2.0;

    double s_hat_mag = sqrt( 2.0*( sij_hat_d[0]*sij_hat_d[0] + 
          sij_hat_d[1]*sij_hat_d[1] + 
          sij_hat_d[2]*sij_hat_d[2] ) +
        4.0*( sij_hat_od[0]*sij_hat_od[0] + 
          sij_hat_od[1]*sij_hat_od[1] + 
          sij_hat_od[2]*sij_hat_od[2] ) );

    for (int i = 0; i < 3; i++) {
      mij_d[icv][i]  = 2.0*( mij_d[icv][i]  - filter_ratio_squared*rho_hat[icv]*s_hat_mag*sij_hat_d[i] );
      mij_od[icv][i] = 2.0*( mij_od[icv][i] - filter_ratio_squared*rho_hat[icv]*s_hat_mag*sij_hat_od[i] );
    }

    // compute MijMij...
    mijmij[icv] = ( mij_d[icv][0]*mij_d[icv][0] + 
        mij_d[icv][1]*mij_d[icv][1] + 
        mij_d[icv][2]*mij_d[icv][2] ) 
      + 2.0*( mij_od[icv][0]*mij_od[icv][0] + 
          mij_od[icv][1]*mij_od[icv][1] + 
          mij_od[icv][2]*mij_od[icv][2] );

  }
  delete[] duidxj_hat;

  // -------------------------------------------
  // now the Lenard term...
  // L_ij = hat(rhou_i u_j) - hat(rhou_i) hat(rhou_j) / hat(rho)
  // -------------------------------------------

  // use sij to store hat(rhou_i u_j)...

  double (*rhouiuj_d)[3] = new double[ncv_g][3];
  double (*rhouiuj_od)[3] = new double[ncv_g][3];
  FOR_ICV_G {
    rhouiuj_d[icv][0] = rho[icv]*u[icv][0]*u[icv][0];
    rhouiuj_d[icv][1] = rho[icv]*u[icv][1]*u[icv][1];
    rhouiuj_d[icv][2] = rho[icv]*u[icv][2]*u[icv][2];
    rhouiuj_od[icv][0] = rho[icv]*u[icv][2]*u[icv][1];
    rhouiuj_od[icv][1] = rho[icv]*u[icv][0]*u[icv][2];
    rhouiuj_od[icv][2] = rho[icv]*u[icv][1]*u[icv][0];
  }
  filterCvR2_mod(sij_d,rhouiuj_d);
  filterCvR2_mod(sij_od,rhouiuj_od);
  delete[] rhouiuj_d;
  delete[] rhouiuj_od;

  // assemble the Leonard term...
  // L_ij =  hat(rhou_i u_j) - hat(rhou_i) hat(rhou_j) / hat(rho)

  for (int icv = 0; icv < ncv; icv++) {

    // diagonal part...
    for (int i = 0; i < 3; i++)
      sij_d[icv][i] -= rhoui_hat[icv][i]*rhoui_hat[icv][i]/rho_hat[icv];

    // remove trace...
    double tmp = (sij_d[icv][0] + sij_d[icv][1] + sij_d[icv][2])/3.0;
    sij_d[icv][0] -= tmp;
    sij_d[icv][1] -= tmp;
    sij_d[icv][2] -= tmp;

    // off-diagonal part...
    sij_od[icv][0] -= rhoui_hat[icv][1]*rhoui_hat[icv][2]/rho_hat[icv];
    sij_od[icv][1] -= rhoui_hat[icv][2]*rhoui_hat[icv][0]/rho_hat[icv];
    sij_od[icv][2] -= rhoui_hat[icv][0]*rhoui_hat[icv][1]/rho_hat[icv];

    // compute LijMij...
    lijmij[icv] = ( sij_d[icv][0]*mij_d[icv][0] + 
        sij_d[icv][1]*mij_d[icv][1] + 
        sij_d[icv][2]*mij_d[icv][2] ) 
      + 2.0*( sij_od[icv][0]*mij_od[icv][0] + 
          sij_od[icv][1]*mij_od[icv][1] + 
          sij_od[icv][2]*mij_od[icv][2] );
  }

  // cleanup...

  delete[] sij_d;
  delete[] sij_od;
  delete[] ui_hat;
  delete[] mij_d;
  delete[] mij_od;    
  delete[] rho_hat;
  delete[] rhoui_hat;

}

void HelmholtzSolver::initASgs() {

  FOR_ICV_G a_sgs[icv] = 0.0;
  FOR_ICV_G a_sgs_sponge[icv] = 0.0;
  
  // TODO move this?
  if (Param * param = getParam("ACOUSTIC_SPONGE")) {

    double L = -1.0,n[3],x[3]; 
    bool b_geom = false;
    int iarg = 0;
    int sponge_type = 2;
    while ( iarg < param->size()) { 
      string token = param->getString(iarg++);
      if ( token == "L" ) {
        L = param->getDouble(iarg++);
      } 
      else if ( token == "GEOM") { 
        token = param->getString(iarg++);
        if (token == "PLANE") {
          FOR_I3 x[i] = param->getDouble(iarg++);
          FOR_I3 n[i] = param->getDouble(iarg++);
          b_geom = true;
        }
        else {
          CERR( " > unrecognized SPONGE GEOM token : " << token);
        }
      }
      else if ( token == "RAMP"){
        token = param->getString(iarg++);
        if (token=="CONSTANT"){
          sponge_type = 0;
        }
        else if (token=="LINEAR"){
          sponge_type = 1;
        }
        else if (token=="QUADRATIC") {
          sponge_type = 2;
        }
        else if (token=="LOGISTIC") {
          sponge_type = 3;
        }
        else{
          CERR(" > unrecognized RAMP parameter " << token << ", valid options: LINEAR (default), QUADRATIC");
        }
      }
      else {
        CERR( " > unrecognized SPONGE token : " << token);
      }

    }

    if (L < 0.0) {
      CERR(" > must specify a non-negative SPONGE L <length>");
    }

    if (!b_geom) {
      CERR(" > must specify a SPONGE GEOM. Examples: \n" << 
           "SPONGE GEOM PLANE <x> <y> <z> <nx> <ny> <nz>");
    }

    double my_buf[2] = {0.0,0.0}; // sd_max,delta_max
    FOR_ICV {
      const double dx[3] = DIFF(x_cv[icv],x);
      my_buf[0] = max(my_buf[0],DOT_PRODUCT(n,dx));
      my_buf[1] = max(my_buf[1],pow(vol_cv[icv],1.0/3.0));
    }
    double buf[2];
    MPI_Allreduce(my_buf,buf,2,MPI_DOUBLE,MPI_MAX,mpi_comm);
    
    FOR_ICV_G {
      const double dx[3] = DIFF(x_cv[icv],x);
      const double sd = DOT_PRODUCT(n,dx);
      if (sd > 0.0) {
        switch (sponge_type) {
          case 0:
            a_sgs_sponge[icv] = L*helmholtz_sos;
            break;
          case 1:
            a_sgs_sponge[icv] = L*helmholtz_sos*sd/buf[0];
            break;
          case 2:
            a_sgs_sponge[icv] = L*helmholtz_sos*sd*sd/(buf[0]*buf[0]);
            break;
          case 3:
            a_sgs_sponge[icv] = L*helmholtz_sos/(1.0+exp(-(sd-2.0*buf[1])/(0.5*buf[1]))); 
            break;
          default:
            assert(0); // handled in parsing
        }
      }
    }
    //dumpRange(a_sgs_sponge,ncv,"a_sgs_sponge");
  }
}

void HelmholtzSolver::calcASgs() {

  FOR_ICV_G a_sgs[icv] = a_sgs_sponge[icv];

  // add some subgrid acoustic modeling here...
  const double a_sgs_coeff = getDoubleParam("A_SGS_COEFF",0.0); // moving uses 0.1
  FOR_ICV_G a_sgs[icv] += a_sgs_coeff*pow(vol_cv[icv],1.0/3.0)*helmholtz_sos; 

}
