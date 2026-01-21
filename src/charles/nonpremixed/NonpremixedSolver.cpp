
#include "NonpremixedSolver.hpp"

#define FOR_BCZONE for(vector<NonpremixedBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)

void NonpremixedSolver::initData() {

  COUT1("NonpremixedSolver::initData()");

  FlowSolver::initData();

  assert(rho    == NULL);       rho  = new double[ncv_g2];
  assert(u      == NULL);       u    = new double[ncv_g2][3];
  assert(rhoE   == NULL);       rhoE = new double[ncv_g2];
  assert(Z      == NULL);       Z    = new double[ncv_g2];
  assert(Zvar   == NULL);       Zvar = new double[ncv_g2];
  assert(C      == NULL);       C    = new double[ncv_g2];
  assert(p      == NULL);       p    = new double[ncv_g2];
  assert(T      == NULL);       T    = new double[ncv_g2];
  assert(sos    == NULL);       sos  = new double[ncv_g2];
  assert(csrc   == NULL);       csrc = new double[ncv_g2];

  // since the viscous closures only involve the compact
  // faces, these structures need only be populated into 
  // the first level ghosts...

  assert(mu_lam   == NULL);     mu_lam  = new double[ncv_g2];
  assert(mu_sgs   == NULL);     mu_sgs  = new double[ncv_g2];
  assert(loc_lam  == NULL);     loc_lam = new double[ncv_g2];
  assert(loc_sgs  == NULL);     loc_sgs = new double[ncv_g2];
  assert(a_sgs    == NULL);     a_sgs   = new double[ncv_g2];

  assert(R        == NULL);     R       = new double[ncv_g2];
  assert(e_p      == NULL);     e_p     = new double[ncv_g2];
  assert(gamma    == NULL);     gamma   = new double[ncv_g2];
  assert(a_gam    == NULL);     a_gam   = new double[ncv_g2];
  assert(T_p      == NULL);     T_p     = new double[ncv_g2];

  assert(ratio    == NULL);     ratio   = new double[ncv_g2];
  assert(Tratio   == NULL);     Tratio  = new double[ncv_g2];
  assert(pratio   == NULL);     pratio  = new double[ncv_g2];

  cv_light            = new NonpremixedState[ncv_g2];
  cv_compact          = new AuxGasProp[ncv_g2];

  // note that this special comm container does not update 
  // the mixture fraction Z.  if a non-algebraic closure 
  // for the mixture fraction variance is introduced, this 
  // should be revisted.

  comm_container      = new NonpremixedComm(rho,u,rhoE,C);
  sgs_comm_container  = new SgsComm(mu_sgs,a_sgs,loc_sgs);

  rhs0                = new NonpremixedRhs[ncv];
  rhs1                = new NonpremixedRhs[ncv];
  rhs2                = new NonpremixedRhs[ncv];

  // allocate buffers required for the scalar transport

  initScalars();

  mf           = new double[nef];

  if (nsc_transport > 0) {

    assert (rhs0_sc == NULL); rhs0_sc = new double[ncv*nsc_transport];
    assert (rhs1_sc == NULL); rhs1_sc = new double[ncv*nsc_transport];
    assert (rhs2_sc == NULL); rhs2_sc = new double[ncv*nsc_transport];

    // we are transport rho*(scalar) , so in the rk update, we 
    // would need a copy of the unadvanced density

    rho_old    = new double[ncv];

  }

  dudx = new double[ncv_g][3][3];

}

inline double gr_func(const double &x, const double& nrm) {

  const double l_width = 5.0e-03;// logistic function width param
  const double max_val = 0.1; //0.5
  return max_val/(1.0 + exp(-2.0*(x-nrm)/l_width));

}

void NonpremixedSolver::setFgrCvs() { 

  // parameter that controls how large of a nondimensional gibbs remainder 
  // is considered large.  this value is O(10^-2).  if you decrease the parameter 
  // the schemes become less oscillatory as it aims to disallow entropy violations 
  // the analysis would suggest that the full discrete gibbs remainder should be 
  // used to set the dissipation, but it appears that only the pressure contributions
  // are being used.

  if ( checkParam("NO_DISS") ) { 

    for (int icv = 0; icv < ncv_g2; ++icv) 
      cv_compact[icv].fgr = 0.0;
    return ;

  }

  const double gr_nrm = getDoubleParam("GR_NRM", 0.03);

  double * vv2       = new double[ncv_g2];
  double * logp         = new double[ncv_g2];
  const double inv_pref = 1.0/chemtable->pressure; 

  for (int icv = 0; icv < ncv; ++icv) { 

    vv2[icv]    = 0.0;
    logp[icv]   = log(p[icv]*inv_pref);

  }

  updateCv2DataStart(logp);

  for (int ief = 0; ief < nef_i; ++ief) { 

    const int icv0 = cvoef[ief][0];
    const int icv1 = cvoef[ief][1];

    // the gibbs remainder is non-dimensionalized locally with RT_stag

    const double dp       = p[icv1] - p[icv0];
    const double v_avg    = 0.5*(1.0/rho[icv0] + 1.0/rho[icv1]);
    const double beta_avg = 0.5*(1.0/(R[icv0]*T[icv0]) + 1.0/(R[icv1]*T[icv1]));
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

  delete[] vv2;
  delete[] logp;

}

void NonpremixedSolver::calcRhs(NonpremixedRhs *__restrict__ rhs, double *__restrict__ rhs_sc, const double time, const int rk_stage) {

  for(vector<NonpremixedBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)  
    (*it)->addBoundaryFlux(rhs);

  for (int ief = 0; ief < n_e__i; ++ief) {
    const int icv0           = fa_hashed[ief].cvofa[0];
    const int icv1           = fa_hashed[ief].cvofa[1];
    const int ief_compressed = fa_hashed[ief].ief_compressed;

    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],1,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    NonpremixedRhs flux;
    double u_avg[3], inv_v_avg, p_avg, h_stag, Z_avg, C_avg;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
        u_avg,inv_v_avg,h_stag,p_avg, Z_avg, C_avg);

    mf[ief]        = flux.rho;

    rhs[icv0].rho -= flux.rho;
    for (int i =0; i < 3; ++i) 
      rhs[icv0].rhou[i] -= flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE;
    rhs[icv0].rhoZ -= flux.rhoZ;
    rhs[icv0].rhoC -= flux.rhoC;

    rhs[icv1].rho += flux.rho;
    for (int i =0; i <3 ; ++i) 
      rhs[icv1].rhou[i] += flux.rhou[i];
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

    NonpremixedRhs flux;
    double u_avg[3], inv_v_avg, p_avg, h_stag, Z_avg, C_avg;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
        u_avg,inv_v_avg,h_stag,p_avg,Z_avg, C_avg);

    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    // aggregate into the flux...

    calcInternalCompactFluxNewGrad(flux,cf_geom[icf_compressed],cv_light[icv0],cv_light[icv1],
        cv_compact[icv0],cv_compact[icv1],u_avg,inv_v_avg,h_stag,
        Z_avg, C_avg,dudx[icv0], dudx[icv1]);

    mf[ief]        = flux.rho;

    rhs[icv0].rho -= flux.rho;
    for (int i =0; i < 3; ++i) 
      rhs[icv0].rhou[i] -= flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE;
    rhs[icv0].rhoZ -= flux.rhoZ;
    rhs[icv0].rhoC -= flux.rhoC;

    rhs[icv1].rho += flux.rho;
    for (int i =0; i <3 ; ++i)  
      rhs[icv1].rhou[i] += flux.rhou[i];
    rhs[icv1].rhoE += flux.rhoE;
    rhs[icv1].rhoZ += flux.rhoZ;
    rhs[icv1].rhoC += flux.rhoC;
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

    NonpremixedRhs flux;
    double u_avg[3], inv_v_avg, p_avg, h_stag, Z_avg, C_avg;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
        u_avg,inv_v_avg,h_stag,p_avg, Z_avg, C_avg);

    NonpremixedRhs c_flux;
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

    NonpremixedRhs flux;
    double u_avg[3], inv_v_avg, p_avg, h_stag, Z_avg, C_avg;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
        u_avg,inv_v_avg,h_stag,p_avg,Z_avg,C_avg);

    NonpremixedRhs c_flux;
    calcInternalFluxSymmetricNew(c_flux,ef_geom[ief_compressed],u_avg,inv_v_avg,
        h_stag,p_avg,Z_avg,C_avg,cv_light[icv0],cv_light[icv1]); 

    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    // aggregate into the flux...
    calcInternalCompactFluxNewGrad(flux,cf_geom[icf_compressed],cv_light[icv0],cv_light[icv1],
        cv_compact[icv0],cv_compact[icv1],u_avg,inv_v_avg,h_stag,
        Z_avg,C_avg,dudx[icv0], dudx[icv1]);

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

  // boundary extended faces with no c, no cmpact 
  for (int ief = n_cc_i; ief < n_e__b; ++ief) {
    const int icv0           = fa_hashed[ief].cvofa[0];
    const int icv1           = fa_hashed[ief].cvofa[1];
    const int ief_compressed = fa_hashed[ief].ief_compressed;

    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    NonpremixedRhs flux;
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

    NonpremixedRhs flux;
    double u_avg[3], inv_v_avg,p_avg,h_stag,Z_avg,C_avg;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
        u_avg,inv_v_avg,h_stag,p_avg,Z_avg,C_avg);

    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    // aggregate into the flux...
    calcInternalCompactFluxNewGrad(flux,cf_geom[icf_compressed],cv_light[icv0],cv_light[icv1],
        cv_compact[icv0],cv_compact[icv1],u_avg,inv_v_avg,h_stag,
        Z_avg,C_avg,dudx[icv0], dudx[icv1]);

    mf[ief]        = flux.rho;

    rhs[icv0].rho -= flux.rho;
    for (int i =0; i < 3; ++i) 
      rhs[icv0].rhou[i] -= flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE;
    rhs[icv0].rhoZ -= flux.rhoZ;
    rhs[icv0].rhoC -= flux.rhoC;
  }

  // proc boundary with c, but no compact closures...
  for (int ief = n_c__b; ief < n_ec_b; ++ief) {
    const int icv0           = fa_hashed[ief].cvofa[0];
    const int icv1           = fa_hashed[ief].cvofa[1];
    const int ief_compressed = fa_hashed[ief].ief_compressed;

    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    NonpremixedRhs flux;
    double u_avg[3],inv_v_avg,p_avg,h_stag,Z_avg,C_avg;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
        u_avg,inv_v_avg,h_stag,p_avg,Z_avg,C_avg);

    NonpremixedRhs c_flux;
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

    NonpremixedRhs flux;
    double u_avg[3],inv_v_avg,p_avg,h_stag,Z_avg,C_avg;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
        u_avg,inv_v_avg,h_stag,p_avg,Z_avg,C_avg);

    NonpremixedRhs c_flux;
    calcInternalFluxSymmetricNew(c_flux,ef_geom[ief_compressed],u_avg,inv_v_avg,
        h_stag,p_avg,Z_avg,C_avg,cv_light[icv0],cv_light[icv1]); 

    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    // aggregate into the flux...
    calcInternalCompactFluxNewGrad(flux,cf_geom[icf_compressed],cv_light[icv0],cv_light[icv1],
        cv_compact[icv0],cv_compact[icv1],u_avg,inv_v_avg,h_stag,
        Z_avg,C_avg,dudx[icv0],dudx[icv1]);

    const double frho = flux.rho + c_flux.rho;
    mf[ief]           = frho; 

    rhs[icv0].rho -= frho; 
    for (int i =0; i < 3; ++i)
      rhs[icv0].rhou[i] -= flux.rhou[i] + c_flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE + c_flux.rhoE;
    rhs[icv0].rhoZ -= flux.rhoZ + c_flux.rhoZ;
    rhs[icv0].rhoC -= flux.rhoC + c_flux.rhoC;
  }

  for (int icv = 0; icv < ncv; ++icv) 
    rhs[icv].rhoC += vol_cv[icv]*csrc[icv]*rho[icv];

  calcRhsScalars(rhs_sc, mf, time, rk_stage);

  addSourceHook(rhs,time,rk_stage);

}//calcRhs()

void NonpremixedSolver::calcRhsScalars(double* __restrict__ rhs_sc, const double* __restrict__ mf,
    const double time, const int rk_stage) {

  // immediately return if there are no scalars to consider.

  if (nsc_transport == 0) 
    return;

  // compute the convection and diffusion terms for the scalars -- the passive scalars 
  // are all closed with a unity (laminar and turbulent) Lewis number

  for(vector<NonpremixedBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)  
    (*it)->addScalarBoundaryFlux(rhs_sc);

  // we consider all of the extended faces that do not have viscous closures (either 
  // with finite c or without).  the effect of the non-zero c has already been comprehended
  // in the calculation of the mass flux which has been stored in mf 

  for (int iter = 0; iter < 2; ++iter) {

    int ief_start, ief_end;

    if (iter == 0) {
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
    if (iter == 0) {
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

      const double loc_total   = 0.5*(loc_lam[icv0] + loc_sgs[icv0] + loc_lam[icv1] + loc_sgs[icv1]);
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
    if (iter == 0) {
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
    if (iter == 0) {
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

      const double loc_total   = 0.5*(loc_lam[icv0] + loc_sgs[icv0] + loc_lam[icv1] + loc_sgs[icv1]);
      const double d_coeff     = loc_total *cf_geom[icf_compressed].area_over_delta;

      for (int isc = 0; isc < nsc_transport; ++isc) {

        double flux            = 0.5*(transport_scalar_vec[isc][icv0] + transport_scalar_vec[isc][icv1])*mf[ief];
        flux                  -= d_coeff*(transport_scalar_vec[isc][icv1] - transport_scalar_vec[isc][icv0]);
        rhs_sc[isc+icv0*nsc_transport]  -= flux;
        
      }
    }
  }

}//calcRhsScalars()

void NonpremixedSolver::advanceScalars(const double * rho_old, const double rk_wgt[3], const double time, const int rk_stage) {

  if (nsc_transport == 0) 
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

}//advanceScalars()

void NonpremixedSolver::rk3Step(const double rk_wgt[3], const double dt, const int rk_stage) {

  if (rho_old) {

    // if there are scalars, then we need to record a copy of the unadvanced density()...
    assert(nsc_transport > 0); 
    for (int icv = 0; icv < ncv; ++icv) 
      rho_old[icv] = rho[icv];

  }

  for (int icv = 0; icv < ncv; ++icv) {
    NonpremixedRhs rhs_agg;
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

}//rk3Step()

void NonpremixedSolver::calcSgsAndMaterialProperties() {

  if (sgs_model == "NONE")          
    computeSgsNone();
  else if (sgs_model == "VREMAN")  
    computeSgsVreman();
  else
    assert(0);

  // March 2019: these default values (COEFF_P = 0.5, COEFF_H = 0.1)
  // are based on a variety of flows (jets, supersonic bl) and
  // minimize grid sensitivity/artifacts and reduce high-freq
  // pile-up in acoustic spectra. Change carefully!

  const double coeff_p_delta = getDoubleParam("COEFF_P", 0.5);
  const double coeff_h_delta = getDoubleParam("COEFF_H", 0.1);

  // compute the expansion coefficients for the material properties
  // XXX should we pre-allocate these arrays for computational savings?

  double * a_mu  = new double[ncv_g2];
  double * a_loc = new double[ncv_g2];
  chemtable->lookup(mu_lam,"mu",a_mu,"a_mu",loc_lam,"locp",a_loc,"a_locp",Z,Zvar,C,ncv_g2);

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

  setFgrCvs();

}//calcSgsAndMaterialProperties()

void NonpremixedSolver::updateExtendedState(NonpremixedState* __restrict__ cv_light, AuxGasProp* __restrict__ cv_compact,
    double *__restrict__ p,            double* __restrict__ T,
    const double *__restrict__ rho,    const double (*__restrict__ u)[3],
    const double *__restrict__ rhoE,   const double *__restrict__ Z, 
    const double *__restrict__ C,      const double* __restrict__ R,      
    const double *__restrict__ T_p,    const double* __restrict__ e_p,    
    const double *__restrict__ gamma,  const double * __restrict__ a_gam, 
    const int icv_start,               const int icv_end) {

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
    // original expr is (exp(agamma*(e-e0)/R) - 1.0) / agamma * (gamma[i]-1.0);
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

}//updateExtendedState()

double NonpremixedSolver::calcCfl(double* cfl_, const double dt_) const {

  // functions returns the rank-local cfl maximum..

  bool b_memflag = false;
  double * cfl   = NULL ;
  if (cfl_ == NULL) {
    cfl = new double[ncv];
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
    if (icv1 < ncv) 
      cfl[icv1] += abs(undA) + area*sos_avg;
  }

  double my_cfl_max = 0.0;
  for (int icv = 0; icv < ncv; ++icv) {
    cfl[icv]  *= 0.5*dt_*inv_vol[icv];
    my_cfl_max = max(my_cfl_max,cfl[icv]);
  }

  if (b_memflag) delete[] cfl;
  return my_cfl_max;
}//calcCfl()

NonpremixedSolver::~NonpremixedSolver() {
  DELETE(rho);
  DELETE(u);
  DELETE(rhoE);
  DELETE(Z);
  DELETE(Zvar);
  DELETE(C);
  DELETE(p);
  DELETE(T);
  DELETE(sos);
  DELETE(csrc);

  DELETE(ratio);
  DELETE(Tratio);
  DELETE(pratio);

  DELETE(mu_lam);
  DELETE(mu_sgs);
  DELETE(loc_lam);
  DELETE(loc_sgs);
  DELETE(a_sgs);

  DELETE(rhs0);
  DELETE(rhs1);
  DELETE(rhs2);
  if (comm_container) delete comm_container;
  if (sgs_comm_container) delete sgs_comm_container;
  DELETE(cv_light);
  DELETE(cv_compact);

  DELETE(R);
  DELETE(e_p);
  DELETE(gamma);
  DELETE(a_gam);
  DELETE(T_p);

  FOR_BCZONE delete *it;

  deleteChemtable(chemtable);

  DELETE(mf);
  DELETE(rhs0_sc);
  DELETE(rhs1_sc);
  DELETE(rhs2_sc);
  DELETE(rho_old);

  DELETE(dudx);
}//~NonpremixedSolver()
#undef FOR_BCZONE
