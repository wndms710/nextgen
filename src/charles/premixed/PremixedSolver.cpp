
#include "PremixedSolver.hpp"

#define FOR_BCZONE for(vector<PremixedBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)
void PremixedSolver::initData() {

  COUT1("PremixedSolver::initData()");

  FlowSolver::initData();

  assert( rho == NULL);       rho  = new double[ncv_g2];
  assert( u   == NULL);       u    = new double[ncv_g2][3];
  assert( rhoE == NULL);      rhoE = new double[ncv_g2];
  assert( Z   == NULL);       Z    = new double[ncv_g2];
  assert( C   == NULL);       C    = new double[ncv_g2];
  assert( p   == NULL);       p    = new double[ncv_g2];
  assert( T   == NULL);       T    = new double[ncv_g2];
  assert( sos == NULL);       sos  = new double[ncv_g2];
  assert( ent == NULL);       ent  = new double[ncv_g2];
  assert( h   == NULL);       h    = new double[ncv_g2];

  // since the viscous closures only involve the compact
  // faces, these structures need only be populated into
  // the first level ghosts...

  assert( mu_lam == NULL);    mu_lam  = new double[ncv_g2];
  assert( mu_sgs == NULL);    mu_sgs  = new double[ncv_g2];
  assert( loc_lam == NULL);   loc_lam = new double[ncv_g2];
  assert( loc_sgs == NULL);   loc_sgs = new double[ncv_g2];
  assert( a_sgs == NULL);     a_sgs   = new double[ncv_g2];

  assert( vv2   == NULL);     vv2     = new double[ncv_g2];

  assert( R == NULL);         R       = new double[ncv_g2];
  assert( e_p == NULL);       e_p     = new double[ncv_g2];
  assert( gamma == NULL);     gamma   = new double[ncv_g2];
  assert( a_gam == NULL);     a_gam   = new double[ncv_g2];
  assert( T_p == NULL);       T_p     = new double[ncv_g2];

  // perturbation quantities

  assert( h_prime == NULL);   h_prime = new double[ncv_g2];
  assert( s_prime == NULL);   s_prime = new double[ncv_g2];


  // face-based progress variable source terms and eff

  assert( int_rho_csrc == NULL); int_rho_csrc = new double[nfa];
  assert( e_fa         == NULL); e_fa         = new double[nfa];

  // use the align_new? command here eventually
  // has shown some impact on performance
  cv_light       = new PremixedState[ncv_g2];
  cv_compact     = new AuxGasProp[ncv_g2];

  comm_container = new PremixedComm(rho,u,rhoE,Z,C);
  sgs_comm_container = new SgsComm(mu_sgs,a_sgs,loc_sgs);

  rhs0 = new PremixedRhs[ncv];
  rhs1 = new PremixedRhs[ncv];
  rhs2 = new PremixedRhs[ncv];

  // updates nsc_transport, so must come before rhs_sc allocation


  // allocate buffers required for the scalar transport

  initScalars();

  mf           = new double[nef];

  if ( nsc_transport > 0) {

    assert ( rhs0_sc == NULL); rhs0_sc = new double[ncv*nsc_transport];
    assert ( rhs1_sc == NULL); rhs1_sc = new double[ncv*nsc_transport];
    assert ( rhs2_sc == NULL); rhs2_sc = new double[ncv*nsc_transport];

    // we are transport rho*(scalar) , so in the rk update, we
    // would need a copy of the unadvanced density

    rho_old    = new double[ncv];

  }

  // if the nox model is enabled then we'll add some memory management
  // for the source terms which are split into prmpt and thermal components

  if ( b_nox ) {

    pair< map<string,double*>::iterator,bool> retp = scalar_src_map.insert(
        pair<string,double*>("nox_prompt_src",NULL));

    pair< map<string,double*>::iterator,bool> rett = scalar_src_map.insert(
        pair<string,double*>("nox_thermal_src",NULL));

    // ensure that the assertions were successful and allocate the memory

    assert ( retp.second );
    assert ( rett.second );
    assert( retp.first->second == NULL);  retp.first->second = new double[nfa];
    assert( rett.first->second == NULL);  rett.first->second = new double[ncv];

  }

  dudx = new double[ncv_g][3][3];
}


// set the f(gr) that is used to decide when the fluxes need stabilization ...
// this calculation loops the extended faces of all the cells...

inline double gr_func(const double &x, const double& nrm) {


  const double l_width = 5.0e-03; // logistic function width param
  const double max_val = 0.1; // 0.5
  return max_val/(1.0+exp(-2.0*(x-nrm)/l_width));

}

void PremixedSolver::setFgrCvs() {

  // parameter that controls how large of a nondimensional gibbs remainder
  // is considered large.  this value is O(10^-2).  if you decrease the parameter
  // the schemes become less oscillatory as it aims to disallow entropy violations
  // the analysis would suggest that the full discrete gibbs remainder should be
  // used to set the dissipation, but it appears that only the pressure contributions
  // are being used.


  const double gr_nrm = getDoubleParam("GR_NRM", 0.03);

  // double * vv2 = new double[ncv_g2];

  for (int icv = 0; icv < ncv; ++icv)
    vv2[icv]    = 0.0;

  double * logp         = new double[ncv_g2];
  const double inv_pref = 1.0/chemtable->pressure;

  for (int icv = 0; icv < ncv; ++icv)
    logp[icv] = log(p[icv]*inv_pref);

  updateCv2DataStart(logp);

  for (int ief = 0; ief < nef_i; ++ief) {

    const int icv0 = cvoef[ief][0];
    const int icv1 = cvoef[ief][1];

    // the gibbs remainder is non-dimensionalized locally with RT_stag

    const double dp       = p[icv1] - p[icv0];
    const double v_avg    = 0.5*(1.0/rho[icv0] + 1.0/rho[icv1]);
    const double beta_avg = 0.5*(1.0/(R[icv0]*T[icv0]) + 1.0/(R[icv1]*T[icv1]));
    //const double gr       = log(p[icv1]/p[icv0]) - v_avg*dp*beta_avg;
    const double gr       = logp[icv1] - logp[icv0] - v_avg*dp*beta_avg;

    // icv0, icv1 are both valid in this loop ..

    vv2[icv0] = max(vv2[icv0],fabs(gr));
    vv2[icv1] = max(vv2[icv1],fabs(gr));

  }

  updateCv2DataFinish(logp);

  // complete in the ghosts ..

  for (int ief = nef_i; ief < nef; ++ief) {

    const int icv0 = cvoef[ief][0];
    const int icv1 = cvoef[ief][1];

    // the gibbs remainder is non-dimensionalized locally with RT_stag

    const double dp       = p[icv1] - p[icv0];
    const double v_avg    = 0.5*(1.0/rho[icv0] + 1.0/rho[icv1]);
    const double beta_avg = 0.5*(1.0/(R[icv0]*T[icv0]) + 1.0/(R[icv1]*T[icv1]));
    //const double gr       = log(p[icv1]/p[icv0]) - v_avg*dp*beta_avg;
    const double gr       = logp[icv1] - logp[icv0] - v_avg*dp*beta_avg;

    // just icv0 this time around ..

    vv2[icv0] = max(vv2[icv0],fabs(gr));
  }

  updateCv2DataStart(vv2);

  for (int icv = 0; icv < ncv; ++icv)
    cv_compact[icv].fgr = gr_func(vv2[icv],gr_nrm);

  updateCv2DataFinish(vv2);

  for (int icv = ncv; icv < ncv_g2; ++icv)
    cv_compact[icv].fgr = gr_func(vv2[icv],gr_nrm);

  if ( checkParam("NO_DISS") ) {
    for (int icv = 0; icv < ncv_g2; ++icv)
      cv_compact[icv].fgr = 0.0;
  }

  //delete[] vv2;

  delete[] logp;
}

void PremixedSolver::setFgrCvsFull() {

  // parameter that controls how large of a nondimensional gibbs remainder
  // is considered large.  this value is O(10^-2).  if you decrease the parameter
  // the schemes become less oscillatory as it aims to disallow entropy violations

  const double gr_nrm = getDoubleParam("GR_NRM", 0.1);

  double * T_stag = new double[ncv_g2];


  for (int icv = 0; icv < ncv; ++icv)  {

    vv2[icv]    = 0.0;
    T_stag[icv] = 1.0 + 0.5*(gamma[icv]-1.0)*(DOT_PRODUCT(u[icv],u[icv]))/(sos[icv]*sos[icv]);

  }

  updateCv2DataStart(T_stag);

  for (int ief = 0; ief < nef_i; ++ief) {

    const int icv0 = cvoef[ief][0];
    const int icv1 = cvoef[ief][1];

    // the gibbs remainder is non-dimensionalized locally with RT_stag

    const double Ts_avg   = 0.5*(T_stag[icv0] + T_stag[icv1]);

    const double dp       = p[icv1] - p[icv0];
    const double v_avg    = 0.5*(1.0/rho[icv0] + 1.0/rho[icv1]);
    const double R_avg    = 0.5*(R[icv0] + R[icv1]);

    const double ds       = s_prime[icv1] - s_prime[icv0];
    const double ds0      = ent[icv1]     - ent[icv0];
    const double dh       = h_prime[icv1] - h_prime[icv0];
    const double theta_f  = 0.5*(1.0/T_p[icv1] + 1.0/T_p[icv0]);
    const double theta    = 0.5*(1.0/T[icv1]   + 1.0/T[icv0]);
    const double theta_p  = theta - theta_f;

    double gr              = ds - theta*( dh - v_avg*dp) - theta_p*ds0/theta_f;
    gr                    /= (R_avg*Ts_avg);

    // icv0, icv1 are both valid in this loop ..

    vv2[icv0] = max(vv2[icv0],fabs(gr));
    vv2[icv1] = max(vv2[icv1],fabs(gr));

  }

  updateCv2DataFinish(T_stag);

  // complete in the ghosts ..

  for (int ief = nef_i; ief < nef; ++ief) {

    const int icv0 = cvoef[ief][0];
    const int icv1 = cvoef[ief][1];

    // the gibbs remainder is non-dimensionalized locally with RT_stag

    const double Ts_avg   = 0.5*(T_stag[icv0] + T_stag[icv1]);

    const double dp       = p[icv1] - p[icv0];
    const double v_avg    = 0.5*(1.0/rho[icv0] + 1.0/rho[icv1]);
    const double R_avg    = 0.5*(R[icv0] + R[icv1]);

    const double ds       = s_prime[icv1] - s_prime[icv0];
    const double ds0      = ent[icv1]     - ent[icv0];
    const double dh       = h_prime[icv1] - h_prime[icv0];
    const double theta_f  = 0.5*(1.0/T_p[icv1] + 1.0/T_p[icv0]);
    const double theta    = 0.5*(1.0/T[icv1]   + 1.0/T[icv0]);
    const double theta_p  = theta - theta_f;

    double gr              = ds - theta*( dh - v_avg*dp) - theta_p*ds0/theta_f;
    gr                    /= (R_avg*Ts_avg);

    // just icv0 this time around ..

    vv2[icv0] = max(vv2[icv0],fabs(gr));
  }

  updateCv2DataStart(vv2);

  for (int icv = 0; icv < ncv; ++icv)
    cv_compact[icv].fgr = gr_func(vv2[icv],gr_nrm);

  updateCv2DataFinish(vv2);

  for (int icv = ncv; icv < ncv_g2; ++icv)
    cv_compact[icv].fgr = gr_func(vv2[icv],gr_nrm);

  delete[] T_stag;

}


void PremixedSolver::calcRhs(PremixedRhs *__restrict__ rhs, double *__restrict__ rhs_sc, const double time, const int rk_stage) {

  for(vector<PremixedBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)
    (*it)->addBoundaryFlux(rhs);

  for (int ief = 0; ief < n_e__i; ++ief) {
    const int icv0           = fa_hashed[ief].cvofa[0];
    const int icv1           = fa_hashed[ief].cvofa[1];
    const int ief_compressed = fa_hashed[ief].ief_compressed;

    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],1,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    PremixedRhs flux;
    double u_avg[3], inv_v_avg, p_avg, h_stag, Z_avg, C_avg;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
        u_avg,inv_v_avg,h_stag,p_avg, Z_avg, C_avg);

    mf[ief]        = flux.rho;

    rhs[icv0].rho -= flux.rho;
    for (int i =0; i < 3; ++i) {
      rhs[icv0].rhou[i] -= flux.rhou[i];
    }
    rhs[icv0].rhoE -= flux.rhoE;
    rhs[icv0].rhoZ -= flux.rhoZ;
    rhs[icv0].rhoC -= flux.rhoC;

    rhs[icv1].rho += flux.rho;
    for (int i =0; i <3 ; ++i) {
      rhs[icv1].rhou[i] += flux.rhou[i];
    }
    rhs[icv1].rhoE += flux.rhoE;
    rhs[icv1].rhoZ += flux.rhoZ;
    rhs[icv1].rhoC += flux.rhoC;
  }

  // these faces have compact contributions to the viscous & modeling terms...

  for (int ief = n_e__i; ief < n_c__i; ++ief) {
    const int icv0           = fa_hashed[ief].cvofa[0];
    const int icv1           = fa_hashed[ief].cvofa[1];
    const int ief_compressed = fa_hashed[ief].ief_compressed;
    const int icf_compressed = fa_hashed[ief].icf_compressed;

    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],1,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    PremixedRhs flux;
    double u_avg[3], inv_v_avg, p_avg, h_stag, Z_avg, C_avg;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
        u_avg,inv_v_avg,h_stag,p_avg,Z_avg,C_avg);

    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    // aggregate into the flux...
    calcInternalCompactFluxNewGrad(flux,cf_geom[icf_compressed],cv_light[icv0],cv_light[icv1],
        cv_compact[icv0],cv_compact[icv1],u_avg,inv_v_avg,h_stag,
        Z_avg,C_avg, dudx[icv0], dudx[icv1]);

    // the progress variable also has a source term applied on compact faces...
    // XXX could pre-multiply e_fa*int_rho_csrc to reduce the amount of data loaded into this loop
    const int ifa                = efsorted2fa[ief]; // could check this ifa >= 0
    //assert( (ifa >=0) && (ifa < nfa_i));
    const double half_rho_csrc   = 0.5*e_fa[ifa]*int_rho_csrc[ifa]*cf_geom[icf_compressed].area;

    mf[ief]        = flux.rho;

    rhs[icv0].rho -= flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv0].rhou[i] -= flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE;
    rhs[icv0].rhoZ -= flux.rhoZ;
    rhs[icv0].rhoC += half_rho_csrc - flux.rhoC;

    rhs[icv1].rho += flux.rho;
    for (int i =0; i <3 ; ++i)
      rhs[icv1].rhou[i] += flux.rhou[i];
    rhs[icv1].rhoE += flux.rhoE;
    rhs[icv1].rhoZ += flux.rhoZ;
    rhs[icv1].rhoC += half_rho_csrc + flux.rhoC;
  }

  // faces have a c component but are purely extended (associated with the convection)..
  for (int ief = n_c__i; ief < n_ec_i; ++ief) {
    const int icv0           = fa_hashed[ief].cvofa[0];
    const int icv1           = fa_hashed[ief].cvofa[1];
    const int ief_compressed = fa_hashed[ief].ief_compressed;

    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],1,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    PremixedRhs flux;
    double u_avg[3], inv_v_avg, p_avg, h_stag, Z_avg, C_avg;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
        u_avg,inv_v_avg,h_stag,p_avg, Z_avg, C_avg);

    PremixedRhs c_flux;
    calcInternalFluxSymmetricNew(c_flux,ef_geom[ief_compressed],u_avg,inv_v_avg,
        h_stag,p_avg,Z_avg,C_avg,cv_light[icv0],cv_light[icv1]);

    // reduce the flux and c_flux
    flux.rho += c_flux.rho;
    for (int i = 0; i < 3; ++i)
      flux.rhou[i] += c_flux.rhou[i];
    flux.rhoE += c_flux.rhoE;
    flux.rhoZ += c_flux.rhoZ;
    flux.rhoC += c_flux.rhoC;

    mf[ief]        = flux.rho;

    rhs[icv0].rho -= flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv0].rhou[i] -= flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE;
    rhs[icv0].rhoZ -= flux.rhoZ;
    rhs[icv0].rhoC -= flux.rhoC;

    rhs[icv1].rho += flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv1].rhou[i] += flux.rhou[i];
    rhs[icv1].rhoE += flux.rhoE;
    rhs[icv1].rhoZ += flux.rhoZ;
    rhs[icv1].rhoC += flux.rhoC;
  }

  // have c, and are also compact faces ...
  for (int ief = n_ec_i; ief < n_cc_i; ++ief) {
    const int icv0           = fa_hashed[ief].cvofa[0];
    const int icv1           = fa_hashed[ief].cvofa[1];
    const int ief_compressed = fa_hashed[ief].ief_compressed;
    const int icf_compressed = fa_hashed[ief].icf_compressed;

    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],1,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    PremixedRhs flux;
    double u_avg[3], inv_v_avg, p_avg, h_stag, Z_avg, C_avg;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
        u_avg,inv_v_avg,h_stag,p_avg,Z_avg,C_avg);

    PremixedRhs c_flux;
    calcInternalFluxSymmetricNew(c_flux,ef_geom[ief_compressed],u_avg,inv_v_avg,
        h_stag,p_avg,Z_avg,C_avg,cv_light[icv0],cv_light[icv1]);

    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    // aggregate into the flux...
    calcInternalCompactFluxNewGrad(flux,cf_geom[icf_compressed],cv_light[icv0],cv_light[icv1],
        cv_compact[icv0],cv_compact[icv1],u_avg,inv_v_avg,h_stag,
        Z_avg,C_avg,dudx[icv0], dudx[icv1]);

    // rxn prog variable source terms ..
    const int ifa                = efsorted2fa[ief]; // could check this ifa >= 0
    //assert( (ifa >= 0) && (ifa < nfa_i));
    const double half_rho_csrc   = 0.5*e_fa[ifa]*int_rho_csrc[ifa]*cf_geom[icf_compressed].area;

    // reduce the flux and c_flux
    flux.rho += c_flux.rho;
    for (int i = 0; i < 3; ++i)
      flux.rhou[i] += c_flux.rhou[i];
    flux.rhoE += c_flux.rhoE;
    flux.rhoZ += c_flux.rhoZ;
    flux.rhoC += c_flux.rhoC;

    mf[ief]        = flux.rho;

    rhs[icv0].rho -= flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv0].rhou[i] -= flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE;
    rhs[icv0].rhoZ -= flux.rhoZ;
    rhs[icv0].rhoC += half_rho_csrc - flux.rhoC;

    rhs[icv1].rho += flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv1].rhou[i] += flux.rhou[i];
    rhs[icv1].rhoE += flux.rhoE;
    rhs[icv1].rhoZ += flux.rhoZ;
    rhs[icv1].rhoC += half_rho_csrc + flux.rhoC;
  }

  // boundary extended faces with no c, no cmpact
  for (int ief = n_cc_i; ief < n_e__b; ++ief) {
    const int icv0           = fa_hashed[ief].cvofa[0];
    const int icv1           = fa_hashed[ief].cvofa[1];
    const int ief_compressed = fa_hashed[ief].ief_compressed;

    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    PremixedRhs flux;
    double u_avg[3], inv_v_avg, p_avg, h_stag, Z_avg, C_avg;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
        u_avg,inv_v_avg,h_stag,p_avg,Z_avg,C_avg);

    mf[ief]        = flux.rho;

    rhs[icv0].rho -= flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv0].rhou[i] -= flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE;
    rhs[icv0].rhoZ -= flux.rhoZ;
    rhs[icv0].rhoC -= flux.rhoC;
  }

  // boundary extended faces with no c, but has compact closure terms
  for (int ief = n_e__b; ief < n_c__b; ++ief) {
    const int icv0           = fa_hashed[ief].cvofa[0];
    const int icv1           = fa_hashed[ief].cvofa[1];
    const int ief_compressed = fa_hashed[ief].ief_compressed;
    const int icf_compressed = fa_hashed[ief].icf_compressed;

    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    PremixedRhs flux;
    double u_avg[3], inv_v_avg,p_avg,h_stag,Z_avg,C_avg;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
        u_avg,inv_v_avg,h_stag,p_avg,Z_avg,C_avg);

    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    // aggregate into the flux...
    calcInternalCompactFluxNewGrad(flux,cf_geom[icf_compressed],cv_light[icv0],cv_light[icv1],
        cv_compact[icv0],cv_compact[icv1],u_avg,inv_v_avg,h_stag,
        Z_avg,C_avg,dudx[icv0], dudx[icv1]);

    // rxn prog variable source terms ..
    const int ifa                = efsorted2fa[ief]; // could check this ifa >= 0
    //assert( (ifa >= nfa_i) && (ifa < nfa));
    const double half_rho_csrc   = 0.5*e_fa[ifa]*int_rho_csrc[ifa]*cf_geom[icf_compressed].area;

    mf[ief]        = flux.rho;

    rhs[icv0].rho -= flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv0].rhou[i] -= flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE;
    rhs[icv0].rhoZ -= flux.rhoZ;
    rhs[icv0].rhoC += half_rho_csrc - flux.rhoC;
  }

  // proc boundary with c, but no compact closures...
  for (int ief = n_c__b; ief < n_ec_b; ++ief) {
    const int icv0           = fa_hashed[ief].cvofa[0];
    const int icv1           = fa_hashed[ief].cvofa[1];
    const int ief_compressed = fa_hashed[ief].ief_compressed;

    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    PremixedRhs flux;
    double u_avg[3],inv_v_avg,p_avg,h_stag,Z_avg,C_avg;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
        u_avg,inv_v_avg,h_stag,p_avg,Z_avg,C_avg);

    PremixedRhs c_flux;
    calcInternalFluxSymmetricNew(c_flux,ef_geom[ief_compressed],u_avg,inv_v_avg,
        h_stag,p_avg,Z_avg,C_avg,cv_light[icv0],cv_light[icv1]);

    const double frho =  flux.rho + c_flux.rho;
    mf[ief]           =  frho;

    rhs[icv0].rho -= frho;
    for (int i =0; i < 3; ++i)
      rhs[icv0].rhou[i] -= flux.rhou[i] + c_flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE + c_flux.rhoE;
    rhs[icv0].rhoZ -= flux.rhoZ + c_flux.rhoZ;
    rhs[icv0].rhoC -= flux.rhoC + c_flux.rhoC;
  }

  // proc boundary with c and compact closures
  for (int ief = n_ec_b; ief < n_cc_b; ++ief) {
    const int icv0           = fa_hashed[ief].cvofa[0];
    const int icv1           = fa_hashed[ief].cvofa[1];
    const int ief_compressed = fa_hashed[ief].ief_compressed;
    const int icf_compressed = fa_hashed[ief].icf_compressed;

    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    PremixedRhs flux;
    double u_avg[3],inv_v_avg,p_avg,h_stag,Z_avg,C_avg;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
        u_avg,inv_v_avg,h_stag,p_avg,Z_avg,C_avg);

    PremixedRhs c_flux;
    calcInternalFluxSymmetricNew(c_flux,ef_geom[ief_compressed],u_avg,inv_v_avg,
        h_stag,p_avg,Z_avg,C_avg,cv_light[icv0],cv_light[icv1]);

    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    // aggregate into the flux...
    calcInternalCompactFluxNewGrad(flux,cf_geom[icf_compressed],cv_light[icv0],cv_light[icv1],
        cv_compact[icv0],cv_compact[icv1],u_avg,inv_v_avg,h_stag,
        Z_avg,C_avg,dudx[icv0],dudx[icv1]);

    // rxn prog variable source terms ..
    const int ifa                = efsorted2fa[ief]; // could check this ifa >= 0
    //assert( (ifa >= nfa_i) && (ifa < nfa));
    const double half_rho_csrc   = 0.5*e_fa[ifa]*int_rho_csrc[ifa]*cf_geom[icf_compressed].area;

    const double frho = flux.rho + c_flux.rho;
    mf[ief]           = frho;

    rhs[icv0].rho -= frho;
    for (int i =0; i < 3; ++i)
      rhs[icv0].rhou[i] -= flux.rhou[i] + c_flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE + c_flux.rhoE;
    rhs[icv0].rhoZ -= flux.rhoZ + c_flux.rhoZ;
    rhs[icv0].rhoC += half_rho_csrc - flux.rhoC - c_flux.rhoC;
  }


  calcRhsScalars(rhs_sc, mf, time, rk_stage);

  addSourceHook(rhs,time,rk_stage);

  addBodyForceSourceTerms(rhs,rhs_sc,time,rk_stage);

}

void PremixedSolver::calcRhsScalars( double* __restrict__ rhs_sc, const double* __restrict__ mf,
    const double time, const int rk_stage) {

  // immediately return if there are no scalars to consider.

  if ( nsc_transport == 0 )
    return;

  // compute the convection and diffusion terms for the scalars -- the passive scalars
  // are all closed with a unity (laminar and turbulent) Lewis number

  for(vector<PremixedBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)
    (*it)->addScalarBoundaryFlux(rhs_sc);

  // we consider all of the extended faces that do not have viscous closures (either
  // with finite c or without).  the effect of the non-zero c has already been comprehended
  // in the calculation of the mass flux which has been stored in mf

  for (int iter = 0; iter < 2; ++iter) {

    int ief_start, ief_end;
    if ( iter == 0) {
      ief_start = 0;
      ief_end   = n_e__i;
    } else {
      ief_start = n_c__i;
      ief_end   = n_ec_i;
    }

    for (int ief = ief_start; ief < ief_end; ++ief) {

      const int icv0         = fa_hashed[ief].cvofa[0];
      const int icv1         = fa_hashed[ief].cvofa[1];

      for (int isc = 0; isc < nsc_transport; ++isc) {

        const double flux    = 0.5*(transport_scalar_vec[isc][icv0] + transport_scalar_vec[isc][icv1])*mf[ief];
        rhs_sc[isc+icv0*nsc_transport]     -= flux;
        rhs_sc[isc+icv1*nsc_transport]     += flux;

      }
    }
  }

  // these faces have compact contributions to the viscous & modeling terms...

  for (int iter = 0; iter < 2; ++iter) {

    int ief_start, ief_end;
    if ( iter == 0 ) {
      ief_start = n_e__i;
      ief_end   = n_c__i;
    } else {
      ief_start = n_ec_i;
      ief_end   = n_cc_i;
    }

    for (int ief = ief_start; ief < ief_end; ++ief) {

      const int icv0           = fa_hashed[ief].cvofa[0];
      const int icv1           = fa_hashed[ief].cvofa[1];
      const int icf_compressed = fa_hashed[ief].icf_compressed;

      const double loc_total   = 0.5*( loc_lam[icv0] + loc_sgs[icv0] + loc_lam[icv1] + loc_sgs[icv1]);
      const double d_coeff     = loc_total *cf_geom[icf_compressed].area_over_delta;

      for (int isc = 0; isc < nsc_transport; ++isc) {

        double flux            = 0.5*(transport_scalar_vec[isc][icv0] + transport_scalar_vec[isc][icv1])*mf[ief];
        flux                  -= d_coeff*(transport_scalar_vec[isc][icv1] - transport_scalar_vec[isc][icv0]);

        rhs_sc[isc+icv0*nsc_transport]  -= flux;
        rhs_sc[isc+icv1*nsc_transport]  += flux;

      }
    }
  }

  // (processor) boundary extended faces with no c, no cmpact

  for (int iter = 0; iter < 2; ++iter) {

    int ief_start, ief_end;
    if ( iter == 0) {
      ief_start = n_cc_i;
      ief_end   = n_e__b;
    } else {
      ief_start = n_c__b;
      ief_end   = n_ec_b;
    }

    for (int ief = ief_start; ief < ief_end; ++ief) {

      const int icv0         = fa_hashed[ief].cvofa[0];
      const int icv1         = fa_hashed[ief].cvofa[1];

      for (int isc = 0; isc < nsc_transport; ++isc) {

        const double flux      = 0.5*(transport_scalar_vec[isc][icv0] + transport_scalar_vec[isc][icv1])*mf[ief];
        rhs_sc[isc+icv0*nsc_transport]  -= flux;

      }
    }
  }

  // these (proc boundary) faces have compact contributions to the viscous & modeling terms...

  for (int iter = 0; iter < 2; ++iter) {

    int ief_start, ief_end;
    if ( iter == 0 ) {
      ief_start = n_e__b;
      ief_end   = n_c__b;
    } else {
      ief_start = n_ec_b;
      ief_end   = n_cc_b;
    }

    for (int ief = ief_start; ief < ief_end; ++ief) {

      const int icv0           = fa_hashed[ief].cvofa[0];
      const int icv1           = fa_hashed[ief].cvofa[1];
      const int icf_compressed = fa_hashed[ief].icf_compressed;

      const double loc_total   = 0.5*( loc_lam[icv0] + loc_sgs[icv0] + loc_lam[icv1] + loc_sgs[icv1]);
      const double d_coeff     = loc_total *cf_geom[icf_compressed].area_over_delta;

      for (int isc = 0; isc < nsc_transport; ++isc) {

        double flux            = 0.5*(transport_scalar_vec[isc][icv0] + transport_scalar_vec[isc][icv1])*mf[ief];
        flux                  -= d_coeff*(transport_scalar_vec[isc][icv1] - transport_scalar_vec[isc][icv0]);
        rhs_sc[isc+icv0*nsc_transport]  -= flux;

      }
    }
  }


  for (map<string,ScalarField*>::iterator it = scalar_map.begin(); it != scalar_map.end(); ++it) {

    const int isc = it->second->idx;

    if ( it->second->name == "Y_NOX") {

      assert( (isc >= 0) && (isc < nsc_transport));

      // since nox has such special source terms, these should have
      // been computed in the updateConservativeAndPrimitiveData ..

      double * prmpt_src = scalar_src_map["nox_prompt_src"];  assert( prmpt_src != NULL);
      double * thrml_src = scalar_src_map["nox_thermal_src"]; assert( thrml_src != NULL);

      // thermal nox source...

      for (int icv = 0; icv < ncv; ++icv)
        rhs_sc[isc+icv*nsc_transport] += vol_cv[icv]*thrml_src[icv]*rho[icv];

      // prompt nox source -- recall that this source is defined on the faces...

      for (int ifa = 0; ifa < nfa_i; ++ifa) {

        const int icv0        = cvofa[ifa][0];
        const int icv1        = cvofa[ifa][1];

        const double half_src = 0.5*MAG(n_fa[ifa])*e_fa[ifa]*prmpt_src[ifa];
        rhs_sc[isc+icv0*nsc_transport]      += half_src;
        rhs_sc[isc+icv1*nsc_transport]      += half_src;

      }

      for (int ifa = nfa_i; ifa < nfa; ++ifa) {

        const int icv0        = cvofa[ifa][0];
        const double half_src = 0.5*MAG(n_fa[ifa])*e_fa[ifa]*prmpt_src[ifa];
        rhs_sc[isc+icv0*nsc_transport]      += half_src;

      }

    } else if ( it->second->name == "RES_TIME") {

      assert((isc >= 0) && (isc < nsc_transport));

      // source term is rho ..

      for (int icv = 0; icv < ncv; ++icv)
        rhs_sc[isc+icv*nsc_transport] += vol_cv[icv]*rho[icv];
    }

  }

}

void PremixedSolver::advanceScalars(const double * rho_old, const double rk_wgt[3], const double time, const int rk_stage) {

  if ( nsc_transport == 0)
    return;

  for (int icv = 0; icv < ncv; ++icv) {

    const double tmp     = dt*inv_vol[icv]*dd[icv];
    const int ii         = icv*nsc_transport;
    const double inv_rho = 1.0/rho[icv];

    for (int isc = 0; isc < nsc_transport; ++isc) {
      // assumes that rk_wgt[i] >= rk_stage = 0..
      const double rhs_agg = tmp*(rk_wgt[0]*rhs0_sc[ii+isc] + rk_wgt[1]*rhs1_sc[ii+isc] + rk_wgt[2]*rhs2_sc[ii+isc]);
      transport_scalar_vec[isc][icv] = (rho_old[icv]*transport_scalar_vec[isc][icv] + rhs_agg)*inv_rho;
    }
  }
}

void PremixedSolver::rk3Step(const double rk_wgt[3], const double dt, const int rk_stage) {

  if ( rho_old ) {

    // if there are scalars , then we need to record a copy of the unadvanced density...
    assert( nsc_transport > 0);
    for (int icv = 0; icv < ncv; ++icv)
      rho_old[icv] = rho[icv];

  }

  for (int icv = 0; icv < ncv; ++icv) {
    PremixedRhs rhs_agg;
    const double tmp = dt*inv_vol[icv]*dd[icv];
    // assumes that rk_wgt[i] >= rk_stage = 0..
    rhs_agg.rho     = tmp*(rk_wgt[0]*rhs0[icv].rho     + rk_wgt[1]*rhs1[icv].rho     + rk_wgt[2]*rhs2[icv].rho);
    rhs_agg.rhou[0] = tmp*(rk_wgt[0]*rhs0[icv].rhou[0] + rk_wgt[1]*rhs1[icv].rhou[0] + rk_wgt[2]*rhs2[icv].rhou[0]);
    rhs_agg.rhou[1] = tmp*(rk_wgt[0]*rhs0[icv].rhou[1] + rk_wgt[1]*rhs1[icv].rhou[1] + rk_wgt[2]*rhs2[icv].rhou[1]);
    rhs_agg.rhou[2] = tmp*(rk_wgt[0]*rhs0[icv].rhou[2] + rk_wgt[1]*rhs1[icv].rhou[2] + rk_wgt[2]*rhs2[icv].rhou[2]);
    rhs_agg.rhoE    = tmp*(rk_wgt[0]*rhs0[icv].rhoE    + rk_wgt[1]*rhs1[icv].rhoE    + rk_wgt[2]*rhs2[icv].rhoE);
    rhs_agg.rhoZ    = tmp*(rk_wgt[0]*rhs0[icv].rhoZ    + rk_wgt[1]*rhs1[icv].rhoZ    + rk_wgt[2]*rhs2[icv].rhoZ);
    rhs_agg.rhoC    = tmp*(rk_wgt[0]*rhs0[icv].rhoC    + rk_wgt[1]*rhs1[icv].rhoC    + rk_wgt[2]*rhs2[icv].rhoC);

    const double rho_old = rho[icv];
    rho[icv]            += rhs_agg.rho;
    const double inv_rho = 1.0/rho[icv];
    for (int i =0; i < 3; ++i)
      u[icv][i]          = (rho_old*u[icv][i] + rhs_agg.rhou[i])*inv_rho;
    Z[icv]               = (rho_old*Z[icv] + rhs_agg.rhoZ)*inv_rho;
    C[icv]               = (rho_old*C[icv] + rhs_agg.rhoC)*inv_rho;
    rhoE[icv]           += rhs_agg.rhoE;
  }

  advanceScalars(rho_old, rk_wgt, dt, rk_stage);

}

void PremixedSolver::calcSgsAndMaterialProperties() {

  if ( sgs_model == "NONE") {
    computeSgsNone();
  } else if ( sgs_model == "VREMAN") {
    computeSgsVreman();
  } else{
    assert(0);
  }

  // March 2019: these default values (COEFF_P = 0.5, COEFF_H = 0.1)
  // are based on a variety of flows (jets, supersonic bl) and
  // minimize grid sensitivity/artifacts and reduce high-freq
  // pile-up in acoustic spectra. Change carefully!

  const double coeff_p_delta = getDoubleParam("COEFF_P", 0.5); // was 1?
  const double coeff_h_delta = getDoubleParam("COEFF_H", 0.1);

  // compute the expansion coefficients for the material properties
  // XXX should we pre-allocate these arrays for computational savings?

  double * a_mu  = new double[ncv_g2];
  double * a_loc = new double[ncv_g2];
  chemtable->lookup(mu_lam,"mu",a_mu,"a_mu",loc_lam,"locp",a_loc,"a_locp",Z,C,ncv_g2);

  // we could consider breaking this loop if we need to hide addtl latency..
  for (int icv = 0; icv < ncv_g2; ++icv) {
    const double TT           = T[icv]/T_p[icv]; // wrt to the table pressure..
    mu_lam[icv]              *= fast_pow_posx(TT,a_mu[icv]);
    loc_lam[icv]             *= fast_pow_posx(TT,a_loc[icv]);

    cv_compact[icv].mu_total  = mu_lam[icv] + mu_sgs[icv];
    cv_compact[icv].loc_total = loc_lam[icv] + loc_sgs[icv];
    cv_compact[icv].a_sgs     = a_sgs[icv];

    // add the contribution associated with the symmetric divergence stabilization..
    const double sos2           = sos[icv]*sos[icv];
    const double lambda         = MAG(u[icv]) + sos[icv]; // u+c characteristic..
    //cv_compact[icv].a_sgs      += coeff_p_delta*dcmag_hat_cv[icv]*lambda/sos2;
    //cv_compact[icv].loc_total  += coeff_h_delta*dcmag_hat_cv[icv]*rho[icv]*MAG(u[icv]);

    // add a delta based on the existing values of the sgs fctors

    const double dx2_approx     = pow(vol_cv[icv], 2.0/3.0);
    cv_compact[icv].symp_fax    = max(coeff_p_delta*dcmag_hat_cv[icv]*lambda/sos2-dx2_approx*a_sgs[icv],0.0);

    // holding the symmetric div part of the enthalpy stabliziation separately..
    // this was previously bundled with the locp contributions (commented line below)

    cv_compact[icv].symh_fax     = coeff_h_delta*dcmag_hat_cv[icv]; // rho u factor not included here...
  }

  delete[] a_mu;
  delete[] a_loc;

  // compute the model for the turbulent flame efficiency factor..

  if ( eff_model == "LAMINAR") {
    computeEfficiencyLaminar();
  } else if ( eff_model == "CHARLETTE") {
    computeEfficiencyCharlette();
  } else if ( eff_model == "OLD") {
    computeEfficiencyCharletteOld(); // for debugging purposes only
  } else {
    // this error was checked earlier..
    assert(0);
  }

  // allow for a global eff factor...

  if (Param * param = getParam("EFF_FACTOR")) {
    const double eff_factor = param->getDouble();
    for (int ifa = 0; ifa < nfa; ++ifa)
      e_fa[ifa] *= eff_factor;
  }

  if ( checkParam("USE_FULL_GR"))
    setFgrCvsFull();
  else
    setFgrCvs();

}

void PremixedSolver::updateExtendedState(PremixedState *__restrict__ cv_light, AuxGasProp *__restrict__ cv_compact,
    double *__restrict__ p, double *__restrict__ T,
    double *__restrict__ h_prime, double *__restrict__ s_prime,
    const double *__restrict__ rho, const double (*__restrict__ u)[3],
    const double *__restrict__ rhoE, const double *__restrict__ Z,
    const double *__restrict__ C, const double* __restrict__ R,
    const double *__restrict__ T_p, const double* __restrict__ e_p,
    const double *__restrict__ gamma, const double * __restrict__ a_gam,
    const double *__restrict__ ent, const int icv_start, const int icv_end) {

  const double p_cmb = chemtable->pressure;

  for (int icv = icv_start; icv < icv_end; ++icv) {
    double usq = 0.0;
    for (int i = 0 ; i < 3; ++i)
      usq += u[icv][i]*u[icv][i];

    // compute p,T..
    const double sp_vol = 1.0/rho[icv];
    const double e_cv   = rhoE[icv]*sp_vol - 0.5*usq;
    const double de     = e_cv - e_p[icv];
    const double deor   = de/R[icv];
    const double adeor  = deor*a_gam[icv];

    // series expansion about a_gam=0 for (exp(a_gam*(e-e0)/R)-1)/a_gam
    // original expr is ( exp( agamma*(e-e0)/R ) - 1.0 ) / agamma * (gamma[i]-1.0);
    const double dT     = (gamma[icv] - 1.0)*deor*(1.0+adeor*0.5000000000000000*
        (1.0+adeor*0.3333333333333333*
         (1.0+adeor*0.2500000000000000*
          (1.0+adeor*0.2000000000000000*
           (1.0+adeor*0.1666666666666667)))));

    T[icv]              = T_p[icv] + dT;
    const double gamma_ = gamma[icv] + a_gam[icv]*dT;
    const double RT_    = T[icv]*R[icv];
    p[icv]              = rho[icv]*RT_;

    sos[icv]            = sqrt(gamma_*RT_);
    h[icv]              = e_cv + RT_;

    // compute the perturbations from the flamelet state ..

    const double dp     = p[icv] - p_cmb;
    const double cp_    = gamma[icv]/(gamma[icv]-1.0)*R[icv];
    h_prime[icv]        = de + R[icv]*dT;
    s_prime[icv]        = cp_* log1p(dT/T_p[icv]) - R[icv]*log1p(dp/p_cmb);

    // pack the condensed thermo state..
    for (int i =0; i < 3; ++i)
      cv_light[icv].u[i] = u[icv][i];
    cv_light[icv].sp_vol = sp_vol;
    cv_light[icv].h      = e_cv + RT_;
    cv_light[icv].p      = p[icv];
    cv_light[icv].Z      = Z[icv];
    cv_light[icv].C      = C[icv];

    // and the rest of the gas properties..
    cv_compact[icv].half_usq = 0.5*usq;
    cv_compact[icv].rho      = rho[icv];
    cv_compact[icv].sos      = sos[icv];
  }

}

double PremixedSolver::calcCfl(double* cfl_, const double dt_) const {

  // functions returns the rank-local cfl maximum..

  bool b_memflag = false;
  double * cfl   = NULL ;
  if ( cfl_ == NULL) {
    cfl       = new double[ncv];
    b_memflag = true;
  } else {
    cfl = cfl_;
  }

  // computation of the convective cfl estimate below.  this
  // does not consider the effect of the extended faces, which
  // are omitted to try to increase the efficiency.  if this
  // cfl estimate is grossly inaccurate, it can be revisited

  for (int icv = 0; icv < ncv; ++icv)
    cfl[icv] = 0.0;

  for (int ibf = 0; ibf < nbf; ++ibf) {
    const int icv = cvobf[ibf];
    const double area = MAG(n_bf[ibf]);
    const double undA = DOT_PRODUCT(u[icv],n_bf[ibf]);
    cfl[icv] += abs(undA) + area*sos[icv];
  }

  for (int ifa = 0; ifa < nfa; ++ifa) {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    const double area = MAG(n_fa[ifa]);
    double undA = 0.5*(DOT_PRODUCT(u[icv0],n_fa[ifa]) + DOT_PRODUCT(u[icv1],n_fa[ifa]));
    const double sos_avg = 0.5*(sos[icv0] + sos[icv1]);
    cfl[icv0] += abs(undA) + area*sos_avg;
    if ( icv1 < ncv)
      cfl[icv1] += abs(undA) + area*sos_avg;
  }

  double my_cfl_max = 0.0;
  for (int icv = 0; icv < ncv; ++icv) {
    cfl[icv]  *= 0.5*dt_*inv_vol[icv];
    my_cfl_max = max(my_cfl_max,cfl[icv]);
  }

  if ( b_memflag) delete[] cfl;
  return my_cfl_max;
}


void PremixedSolver::addBodyForceSourceTerms(PremixedRhs* __restrict__ rhs, double *__restrict__ rhs_sc,
                                            const double time, const int rk_stage) {
  parseBodyForces(); //required, call this first

  if (body_force_flag){
   // if ((step%check_interval == 0 && rk_stage == 1)&&(mpi_rank == 0))
   //   cout << "Adding body force terms..." << endl;
    BodyForce * body_force;
    FOR_ICV {
      if (body_force_flag[icv] >= 0) {
        body_force = &bodyForceVec[body_force_flag[icv]];
        if (body_force->type == PMM_TYPE) {
          const double umag = MAG(u[icv]);
          double f_cv[3] = {0.0,0.0,0.0};
          FOR_I3 {
            FOR_J3 {
              f_cv[i] += vol_cv[icv]*u[icv][j]*(0.5*rho[icv]*umag*body_force->C[i][j]+mu_lam[icv]*body_force->D[i][j]);
            }
            // these may need to be cleared before function is called...
            body_force->force[i] += f_cv[i];
            body_force->work     += f_cv[i]*u[icv][i];

            // add to momentum...
            rhs[icv].rhou[i] -= f_cv[i];
            // and add term to energy...
            rhs[icv].rhoE -= f_cv[i]*u[icv][i];
          }
          body_force->volume   += vol_cv[icv];
          double m_cv[3] = CROSS_PRODUCT(x_cv[icv],f_cv);
          FOR_I3 body_force->moment[i] -= m_cv[i];
        }
        else if (body_force->type == MRF_TYPE) {
          CERR("ADD_BODY_FORCE ... TYPE=MRF ... is not supported by the premixed equation of state");
        }
        else if (body_force->type == HOOK_TYPE) {
          CERR("ADD_BODY_FORCE ... TYPE=HOOK ... is not supported by the premixed equation of state");
        }
        else {
          assert(0);  //logic problem, shouldn't get here
        }
      }
    }
  }
}


PremixedSolver::~PremixedSolver() {
  DELETE(rho);
  DELETE(u);
  DELETE(rhoE);
  DELETE(Z);
  DELETE(C);
  DELETE(p);
  DELETE(T);
  DELETE(sos);
  DELETE(ent);

  DELETE(mu_lam);
  DELETE(mu_sgs);
  DELETE(loc_lam);
  DELETE(loc_sgs);
  DELETE(a_sgs);
  DELETE(vv2);
  DELETE(rhs0);
  DELETE(rhs1);
  DELETE(rhs2);
  if ( comm_container) delete comm_container;
  if ( sgs_comm_container) delete sgs_comm_container;
  DELETE(cv_light);
  DELETE(cv_compact);

  DELETE(R);
  DELETE(h);
  DELETE(e_p);
  DELETE(gamma);
  DELETE(a_gam);
  DELETE(T_p);

  DELETE(h_prime);
  DELETE(s_prime);

  DELETE(int_rho_csrc);
  DELETE(e_fa);
  FOR_BCZONE delete *it;

  deleteChemtable(chemtable);

  DELETE(mf);
  DELETE(rhs0_sc);
  DELETE(rhs1_sc);
  DELETE(rhs2_sc);
  DELETE(rho_old);

  DELETE(dudx);
}
#undef FOR_BCZONE
