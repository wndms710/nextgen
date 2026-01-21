
#include "IdealGasSolver.hpp"
#include "IdealGasSolverFlux.hpp"

// convenience macro to loop over the boundary conditions...
#define FOR_BCZONE for(vector<IdealGasBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)

IdealGasSolver::~IdealGasSolver() {
  DELETE(rho);
  DELETE(u);
  DELETE(rhoE);
  DELETE(p);
  DELETE(T);
  DELETE(sos);
  DELETE(vv2);
  DELETE(fgr_sponge);
  DELETE(mu_lam);
  DELETE(mu_sgs);
  DELETE(loc_lam);
  DELETE(loc_sgs);
  DELETE(a_sgs);
  DELETE(rhs0);
  DELETE(rhs1);
  DELETE(rhs2);

  if ( comm_container) delete comm_container;
  if ( sgs_comm_container) delete sgs_comm_container;
  DELETE(cv_light);
  DELETE(cv_compact);

  FOR_BCZONE delete (*it);

  DELETE(mf);
  DELETE(rhs0_sc);
  DELETE(rhs1_sc);
  DELETE(rhs2_sc);
  DELETE(rho_old);

  DELETE(dudx);
}

void IdealGasSolver::initData() {

  FlowSolver::initData();

  assert(rho == NULL);  rho  = new double[ncv_g2];
  assert(u == NULL);    u    = new double[ncv_g2][3];
  assert(rhoE == NULL); rhoE = new double[ncv_g2];
  assert(p == NULL);    p    = new double[ncv_g2];
  assert(T == NULL);    T    = new double[ncv_g2];
  assert(sos == NULL);  sos  = new double[ncv_g2];
  assert(vv2 == NULL);  vv2  = new double[ncv_g2];

  // set conservative variables to ref values...
  // this avoids bringing the solver down when the
  // user has not specified an initial condition...

  FOR_ICV {
    rho[icv] = rho_ref;
    FOR_I3 u[icv][i] = 0.0;
    rhoE[icv] = p_ref/(gamma-1.0);
  }

  // since the viscous closures only involve the compact
  // faces, these structures need only be populated into
  // the first level ghosts, but put in both layers because
  // updates currently update both layers...

  assert(mu_lam == NULL);  mu_lam  = new double[ncv_g2];
  assert(mu_sgs == NULL);  mu_sgs  = new double[ncv_g2];
  assert(loc_lam == NULL); loc_lam = new double[ncv_g2];
  assert(loc_sgs == NULL); loc_sgs = new double[ncv_g2];
  assert(a_sgs == NULL);   a_sgs   = new double[ncv_g2];

  // packed data structures for flux calcs...

  cv_light       = new IdealGasState[ncv_g2];
  cv_compact     = new AuxGasProp[ncv_g2];

  comm_container     = new IgComm(rho,u,rhoE);
  sgs_comm_container = new SgsComm(mu_sgs,a_sgs,loc_sgs);

  rhs0           = new IdealGasRhs[ncv];
  rhs1           = new IdealGasRhs[ncv];
  rhs2           = new IdealGasRhs[ncv];

  if ( (mpi_rank == 0) && cti_verbose) {
    cout << "Gas properties: " << endl;
    cout << "  gamma        : " << gamma << endl;
    cout << "  p_ref        : " << p_ref << endl;
    cout << "  rho_ref      : " << rho_ref << endl;
    cout << "  T_ref        : " << T_ref << endl;
    cout << "  R_gas        : " << R_gas << endl;
    cout << "  mu_ref       : " << mu_ref << endl;
    cout << "  mu_power_law : " << mu_power_law << endl;
    cout << "  Pr_lam       : " << Pr_lam << endl;
  }

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

  dudx = new double[ncv_g][3][3];

  //==================
  // Lp stuff...
  //==================

  // these pointers should be set after their initialization in the solver
  if (lpTracer!=NULL) {
    lpTracer->u = this->u;
    lpTracer->ncv = ncv;
  }

  if (lpSolid!=NULL) {
    lpSolid->u = this->u;
    lpSolid->ncv = ncv;
  }

  if (lpLiquid!=NULL) {
    lpLiquid->u = this->u;
    lpLiquid->ncv = ncv;
  }

}

void IdealGasSolver::computeFgrSponge() {

  if (Param* param = getParam("DISS_SPONGE")) {
    assert(fgr_sponge == NULL); fgr_sponge = new double[ncv_g2];
    FOR_ICV fgr_sponge[icv] = 0.0;
    SimpleGeom* simple_geom = NULL;
    int iarg = 0;
    while (iarg < param->size()) {
      const string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "GEOM") {
        assert(simple_geom == NULL);
        simple_geom = newSimpleGeom(param,iarg);
        if (simple_geom) {
          FOR_ICV {
            if (simple_geom->pointIsInside(x_cv[icv])) 
              fgr_sponge[icv] = 0.5;
          }
          delete simple_geom; simple_geom = NULL;
        }
      }
      else if (token == "FLIP_GEOM") {
        assert(simple_geom == NULL);
        simple_geom = newSimpleGeom(param,iarg);
        if (simple_geom) {
          FOR_ICV {
            if (!simple_geom->pointIsInside(x_cv[icv])) 
              fgr_sponge[icv] = 0.5;
          }
          delete simple_geom; simple_geom = NULL;
        }
      }
      else {
        CWARN("unrecognized DISS_SPONGE token: " << token << ", skipping.");
      }
    }
    updateCv2Data(fgr_sponge);
  }

}

void IdealGasSolver::calcSgsAndMaterialProperties() {

  if ( sgs_model == "NONE") {
    computeSgsNone();
  } else if ( sgs_model == "VREMAN") {
    computeSgsVreman();
  } else if ( sgs_model == "SIGMA") {
    computeSgsSigma();
  } else if ( sgs_model == "DSM_LOCAL") {
    computeSgsDsm();
  } else {
    // this error is checked earlier, so we shouldnt get here..
    assert(0);
  }

  // March 2019: these default values (COEFF_P = 0.5, COEFF_H = 0.1)
  // are based on a variety of flows (jets, supersonic bl) and
  // minimize grid sensitivity/artifacts and reduce high-freq
  // pile-up in acoustic spectra. Change carefully!

  const double coeff_p_delta = getDoubleParam("COEFF_P", 0.5);
  const double coeff_h_delta = getDoubleParam("COEFF_H", 0.1);

  // set the material properties now, note that the ghost
  // values should have already been set (ncv_g2 level ghosts)

  const double inv_Tref      = 1.0/T_ref;
  const double invPr_lam     = 1.0/Pr_lam;

  for (int icv = 0; icv < ncv_g2; ++icv) {

    //mu_lam[icv]                 = mu_ref*pow(T[icv]*inv_Tref,mu_power_law);

    // the argument of this function must be > 0

    mu_lam[icv]                 = mu_ref*fast_pow_posx(T[icv]*inv_Tref,mu_power_law);

    loc_lam[icv]                = mu_lam[icv]*invPr_lam;

    cv_compact[icv].mu_total    = mu_lam[icv] + mu_sgs[icv];
    cv_compact[icv].loc_total   = loc_lam[icv] + loc_sgs[icv];
    cv_compact[icv].a_sgs       = a_sgs[icv];

    // add the contribution associated with the symmetric divergence stabilization here...

    const double sos2           = gamma*p[icv]/rho[icv];
    const double lambda         = MAG(u[icv]) + sqrt(sos2); // u+c characteristic..
    //cv_compact[icv].a_sgs      += coeff_p_delta*dcmag_hat_cv[icv]*lambda/sos2;
    //cv_compact[icv].loc_total  += coeff_h_delta*dcmag_hat_cv[icv]*rho[icv]*MAG(u[icv]);

    // add a delta based on the existing values of the sgs fctors

    const double dx2_approx     = pow(vol_cv[icv], 2.0/3.0);
    cv_compact[icv].symp_fax    = max(coeff_p_delta*dcmag_hat_cv[icv]*lambda/sos2-dx2_approx*a_sgs[icv],0.0);

    // holding the symmetric div part of the enthalpy stabliziation separately..
    // this was previously bundled with the locp contributions (commented line below)

    cv_compact[icv].symh_fax     = coeff_h_delta*dcmag_hat_cv[icv]; // rho u factor not included here...

  }

  if ( checkParam("USE_FULL_GR"))
    setFgrCvsFull();
  else
    setFgrCvs();

  // user-defined region for applying upwind dissipation...
  if (fgr_sponge) {
    for (int icv = 0 ; icv < ncv_g2; ++icv)
      cv_compact[icv].fgr = max(cv_compact[icv].fgr,fgr_sponge[icv]);
  }

}

void IdealGasSolver::calcRhsBase(IdealGasRhs *__restrict__ rhs, double *__restrict__ rhs_sc,
    const double time, const int rk_stage) {


  for(vector<IdealGasBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)
    (*it)->addBoundaryFlux(rhs);

  for (int ief = 0; ief < n_e__i; ++ief) {
    const int icv0           = fa_hashed[ief].cvofa[0];
    const int icv1           = fa_hashed[ief].cvofa[1];
    const int ief_compressed = fa_hashed[ief].ief_compressed;

    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],1,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    IdealGasRhs flux;
    double u_avg[3], inv_v_avg, p_avg, h_stag;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
                        u_avg,inv_v_avg,h_stag,p_avg);

    mf[ief]        = flux.rho;

    rhs[icv0].rho -= flux.rho;
    rhs[icv1].rho += flux.rho;
    for (int i =0; i <3 ; ++i) {
      rhs[icv0].rhou[i] -= flux.rhou[i];
      rhs[icv1].rhou[i] += flux.rhou[i];
    }
    rhs[icv0].rhoE -= flux.rhoE;
    rhs[icv1].rhoE += flux.rhoE;
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

    IdealGasRhs flux;
    double u_avg[3], inv_v_avg, p_avg, h_stag;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
                        u_avg,inv_v_avg,h_stag,p_avg);

    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    // aggregate into the flux...
    calcInternalCompactFluxNewGrad(flux,cf_geom[icf_compressed],cv_light[icv0],cv_light[icv1],
                                   cv_compact[icv0],cv_compact[icv1],u_avg,inv_v_avg,h_stag,
                                   dudx[icv0], dudx[icv1]);

    mf[ief]        = flux.rho;

    rhs[icv0].rho -= flux.rho;
    rhs[icv1].rho += flux.rho;
    for (int i =0; i <3 ; ++i) {
      rhs[icv0].rhou[i] -= flux.rhou[i];
      rhs[icv1].rhou[i] += flux.rhou[i];
    }
    rhs[icv0].rhoE -= flux.rhoE;
    rhs[icv1].rhoE += flux.rhoE;
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

    IdealGasRhs flux;
    double u_avg[3], inv_v_avg, p_avg, h_stag;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
                        u_avg,inv_v_avg,h_stag,p_avg);

    IdealGasRhs c_flux;
    calcInternalFluxSymmetricNew(c_flux,ef_geom[ief_compressed],u_avg,inv_v_avg,
        h_stag,p_avg,cv_light[icv0],cv_light[icv1]);

    // reduce the flux and c_flux
    flux.rho       += c_flux.rho;
    for (int i =0; i < 3; ++i)
      flux.rhou[i] += c_flux.rhou[i];
    flux.rhoE      += c_flux.rhoE;

    mf[ief]         = flux.rho;

    rhs[icv0].rho -= flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv0].rhou[i] -= flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE;

    rhs[icv1].rho += flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv1].rhou[i] += flux.rhou[i];
    rhs[icv1].rhoE += flux.rhoE;
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

    IdealGasRhs flux;
    double u_avg[3], inv_v_avg, p_avg, h_stag;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
                        u_avg,inv_v_avg,h_stag,p_avg);

    IdealGasRhs c_flux;
    calcInternalFluxSymmetricNew(c_flux,ef_geom[ief_compressed],u_avg,inv_v_avg,
        h_stag,p_avg,cv_light[icv0],cv_light[icv1]);

    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    // aggregate into the flux...
    calcInternalCompactFluxNewGrad(flux,cf_geom[icf_compressed],cv_light[icv0],cv_light[icv1],
                                   cv_compact[icv0],cv_compact[icv1],u_avg,inv_v_avg,h_stag,
                                   dudx[icv0], dudx[icv1]);

    // reduce the flux and c_flux
    flux.rho       += c_flux.rho;
    for (int i =0; i < 3; ++i)
      flux.rhou[i] += c_flux.rhou[i];
    flux.rhoE      += c_flux.rhoE;

    mf[ief]         = flux.rho;

    rhs[icv0].rho -= flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv0].rhou[i] -= flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE;

    rhs[icv1].rho += flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv1].rhou[i] += flux.rhou[i];
    rhs[icv1].rhoE += flux.rhoE;
  }

  // boundary extended faces with no c, no cmpact
  for (int ief = n_cc_i; ief < n_e__b; ++ief) {
    const int icv0           = fa_hashed[ief].cvofa[0];
    const int icv1           = fa_hashed[ief].cvofa[1];
    const int ief_compressed = fa_hashed[ief].ief_compressed;

    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    IdealGasRhs flux;
    double u_avg[3], inv_v_avg, p_avg, h_stag;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
                        u_avg,inv_v_avg,h_stag,p_avg);

    mf[ief]        = flux.rho;

    rhs[icv0].rho -= flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv0].rhou[i] -= flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE;
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

    IdealGasRhs flux;
    double u_avg[3], inv_v_avg,p_avg,h_stag;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
                        u_avg,inv_v_avg,h_stag,p_avg);

    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    // aggregate into the flux...
    calcInternalCompactFluxNewGrad(flux,cf_geom[icf_compressed],cv_light[icv0],cv_light[icv1],
                                   cv_compact[icv0],cv_compact[icv1],u_avg,inv_v_avg,h_stag,
                                   dudx[icv0], dudx[icv1]);

    mf[ief]        = flux.rho;

    rhs[icv0].rho -= flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv0].rhou[i] -= flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE;
  }

  // proc boundary with c, but no compact closures...
  for (int ief = n_c__b; ief < n_ec_b; ++ief) {
    const int icv0           = fa_hashed[ief].cvofa[0];
    const int icv1           = fa_hashed[ief].cvofa[1];
    const int ief_compressed = fa_hashed[ief].ief_compressed;

    cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    IdealGasRhs flux;
    double u_avg[3],inv_v_avg,p_avg,h_stag;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
                        u_avg,inv_v_avg,h_stag,p_avg);

    IdealGasRhs c_flux;
    calcInternalFluxSymmetricNew(c_flux,ef_geom[ief_compressed],u_avg,inv_v_avg,
        h_stag,p_avg,cv_light[icv0],cv_light[icv1]);

    const double frho  = flux.rho + c_flux.rho;
    mf[ief]            = frho;

    rhs[icv0].rho -= frho;
    for (int i =0; i < 3; ++i)
      rhs[icv0].rhou[i] -= flux.rhou[i] + c_flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE + c_flux.rhoE;
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

    IdealGasRhs flux;
    double u_avg[3],inv_v_avg,p_avg,h_stag;
    calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
                        u_avg,inv_v_avg,h_stag,p_avg);

    IdealGasRhs c_flux;
    calcInternalFluxSymmetricNew(c_flux,ef_geom[ief_compressed],u_avg,inv_v_avg,
        h_stag,p_avg,cv_light[icv0],cv_light[icv1]);

    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
    cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);

    // aggregate into the flux...
    calcInternalCompactFluxNewGrad(flux,cf_geom[icf_compressed],cv_light[icv0],cv_light[icv1],
                                   cv_compact[icv0],cv_compact[icv1],u_avg,inv_v_avg,h_stag,
                                   dudx[icv0], dudx[icv1]);

    const double frho    = flux.rho + c_flux.rho;
    mf[ief]              = frho;

    rhs[icv0].rho -= frho;
    for (int i =0; i < 3; ++i)
      rhs[icv0].rhou[i] -= flux.rhou[i] + c_flux.rhou[i];
    rhs[icv0].rhoE -= flux.rhoE + c_flux.rhoE;
  }

}

void IdealGasSolver::calcRhsScalars( double* __restrict__ rhs_sc, const double* __restrict__ mf,
    const double time, const int rk_stage) {

  // immediately return if there are no scalars to consider.

  if ( nsc_transport == 0 )
    return;

  // compute the convection and diffusion terms for the scalars -- the passive scalars
  // are all closed with a unity (laminar and turbulent) Lewis number

  for(vector<IdealGasBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)
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

  // there are some scalars that have source terms that should be applied

  for (map<string,ScalarField*>::iterator it = scalar_map.begin();
      it != scalar_map.end(); ++it) {

    const int isc = it->second->idx;

    // approach to adding the source terms...

    if ( it->second->name == "RES_TIME") {

      // source term is just rho ..

      for (int icv = 0; icv < ncv; ++icv)
        rhs_sc[isc+icv*nsc_transport] += vol_cv[icv]*rho[icv];

    }
  }

}

void IdealGasSolver::addRotationSourceTerms(IdealGasRhs* __restrict__ rhs, double *__restrict__ rhs_sc,
                                            const double time, const int rk_stage) {


  if ( frame_rotation != NULL) {


    for (int icv = 0; icv < ncv; ++icv) {

      // term 1: -omega x omega x r
      // centrifugal force

      double coeff[3];
      coeff[0] = frame_rotation[1]*x_cv[icv][2] - frame_rotation[2]*x_cv[icv][1];
      coeff[1] = frame_rotation[2]*x_cv[icv][0] - frame_rotation[0]*x_cv[icv][2];
      coeff[2] = frame_rotation[0]*x_cv[icv][1] - frame_rotation[1]*x_cv[icv][0];

      double msource[3];
      msource[0] = -rho[icv]*vol_cv[icv]*(frame_rotation[1]*coeff[2] - frame_rotation[2]*coeff[1]);
      msource[1] = -rho[icv]*vol_cv[icv]*(frame_rotation[2]*coeff[0] - frame_rotation[0]*coeff[2]);
      msource[2] = -rho[icv]*vol_cv[icv]*(frame_rotation[0]*coeff[1] - frame_rotation[1]*coeff[0]);

      // term 2: -2 omega x v
      // Coriolis force

      msource[0] -= 2.0*rho[icv]*vol_cv[icv]*(frame_rotation[1]*u[icv][2] - frame_rotation[2]*u[icv][1]);
      msource[1] -= 2.0*rho[icv]*vol_cv[icv]*(frame_rotation[2]*u[icv][0] - frame_rotation[0]*u[icv][2]);
      msource[2] -= 2.0*rho[icv]*vol_cv[icv]*(frame_rotation[0]*u[icv][1] - frame_rotation[1]*u[icv][0]);

      // add to momentum...

      FOR_I3 rhs[icv].rhou[i] += msource[i];

      // and add term to energy...

      rhs[icv].rhoE += DOT_PRODUCT(msource,u[icv]);
    }
  }

}

inline double gr_func(const double &x, const double &nrm) {

  //const double nrm4 = nrm*nrm*nrm*nrm;
  //return erf(x*x*x*x/nrm4);

  const double l_width = 5.0e-03; // logistic function width param
  const double max_val = 0.5;
  return max_val/(1.0+exp(-2.0*(x-nrm)/l_width));

}

void IdealGasSolver::setFgrCvsFull() {

  // parameter that controls how large of a nondimensional gibbs remainder
  // is considered large.  this value is O(10^-2).  if you decrease the parameter
  // the schemes become less oscillatory as it aims to disallow entropy violations

  const double gr_nrm = getDoubleParam("GR_NRM", 0.03);

  double * T_stag = new double[ncv_g2];
  double * ent    = new double[ncv_g2];

  const double cp         = R_gas*gamma/(gamma-1.0);
  const double gm1        = gamma - 1.0;
  const double invgm1     = 1.0/gm1;
  const double inv_pref   = 1.0/p_ref;
  const double inv_rhoref = 1.0/rho_ref;

  for (int icv = 0; icv < ncv; ++icv)  {

    vv2[icv]    = 0.0;
    T_stag[icv] = 1.0 + 0.5*gm1*(DOT_PRODUCT(u[icv],u[icv]))/(sos[icv]*sos[icv]);

    const double tmp1 = log1p((p[icv]-p_ref)*inv_pref); // log(p[icv]*inv_pref)
    const double tmp2 = log1p((rho[icv]-rho_ref)*inv_rhoref); // log(rho[icv]*inv_rhoref)
    ent[icv] = (tmp1-gamma*tmp2)*invgm1;

  }

  updateCv2DataStart(T_stag);
  updateCv2DataStart(ent);

  for (int ief = 0; ief < nef_i; ++ief) {

    const int icv0 = cvoef[ief][0];
    const int icv1 = cvoef[ief][1];

    // the gibbs remainder is non-dimensionalized locally with RT_stag

    const double Ts_avg   = 0.5*(T_stag[icv0] + T_stag[icv1]);
    const double dp       = p[icv1] - p[icv0];
    const double dh       = cp*(T[icv1] - T[icv0]);
    const double v_avg    = 0.5*(1.0/rho[icv0] + 1.0/rho[icv1]);
    const double beta_avg = 0.5*(1.0/(R_gas*T[icv0]) + 1.0/(R_gas*T[icv1]));
    const double dsor     = ent[icv1] - ent[icv0];  // the entropy already has an R factor built in...
    const double gr       = (-dsor + beta_avg*(dh - v_avg*dp))/Ts_avg;

    // icv0, icv1 are both valid in this loop ..

    vv2[icv0] = max(vv2[icv0],fabs(gr));
    vv2[icv1] = max(vv2[icv1],fabs(gr));

  }

  updateCv2DataFinish(T_stag);
  updateCv2DataFinish(ent);

  // complete in the ghosts ..

  for (int ief = nef_i; ief < nef; ++ief) {

    const int icv0 = cvoef[ief][0];
    const int icv1 = cvoef[ief][1];

    // the gibbs remainder is non-dimensionalized locally with RT_stag

    const double Ts_avg   = 0.5*(T_stag[icv0] + T_stag[icv1]);
    const double dp       = p[icv1] - p[icv0];
    const double dh       = cp*(T[icv1] - T[icv0]);
    const double v_avg    = 0.5*(1.0/rho[icv0] + 1.0/rho[icv1]);
    const double beta_avg = 0.5*(1.0/(R_gas*T[icv0]) + 1.0/(R_gas*T[icv1]));
    const double dsor     = ent[icv1] - ent[icv0];  // the entropy already has an R factor built in...
    const double gr       = (-dsor + beta_avg*(dh - v_avg*dp))/Ts_avg;

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
  delete[] ent;

}

void IdealGasSolver::rk3Step(const double rk_wgt[3], const double dt, const int rk_stage) {

  if ( rho_old ) {

    // if there are scalars , then we need to record a copy of the unadvanced density...

    assert( nsc_transport > 0);
    for (int icv = 0; icv < ncv; ++icv)
      rho_old[icv] = rho[icv];

  }

  for (int icv = 0; icv < ncv; ++icv) {
    IdealGasRhs rhs_agg;
    const double tmp = dt*inv_vol[icv];
    // assumes that rk_wgt[i] >= rk_stage = 0..
    rhs_agg.rho     = tmp*(rk_wgt[0]*rhs0[icv].rho     + rk_wgt[1]*rhs1[icv].rho     + rk_wgt[2]*rhs2[icv].rho);
    rhs_agg.rhou[0] = tmp*(rk_wgt[0]*rhs0[icv].rhou[0] + rk_wgt[1]*rhs1[icv].rhou[0] + rk_wgt[2]*rhs2[icv].rhou[0]);
    rhs_agg.rhou[1] = tmp*(rk_wgt[0]*rhs0[icv].rhou[1] + rk_wgt[1]*rhs1[icv].rhou[1] + rk_wgt[2]*rhs2[icv].rhou[1]);
    rhs_agg.rhou[2] = tmp*(rk_wgt[0]*rhs0[icv].rhou[2] + rk_wgt[1]*rhs1[icv].rhou[2] + rk_wgt[2]*rhs2[icv].rhou[2]);
    rhs_agg.rhoE    = tmp*(rk_wgt[0]*rhs0[icv].rhoE    + rk_wgt[1]*rhs1[icv].rhoE    + rk_wgt[2]*rhs2[icv].rhoE);

    const double rho_old = rho[icv];
    rho[icv]            += rhs_agg.rho;
    const double inv_rho = 1.0/rho[icv];
    for (int i =0; i < 3; ++i)
      u[icv][i]          = (rho_old*u[icv][i] + rhs_agg.rhou[i])*inv_rho;
    rhoE[icv]           += rhs_agg.rhoE;
  }

  advanceScalars(rho_old, rk_wgt, dt, rk_stage);

}//rk3Step()

void IdealGasSolver::advanceScalars(const double * rho_old, const double rk_wgt[3],
    const double time, const int rk_stage) {

  if ( nsc_transport == 0)
    return;

  for (int icv = 0; icv < ncv; ++icv) {

    const double tmp     = dt*inv_vol[icv];
    const int ii         = icv*nsc_transport;
    const double inv_rho = 1.0/rho[icv];

    for (int isc = 0; isc < nsc_transport; ++isc) {
      // assumes that rk_wgt[i] >= rk_stage = 0..
      const double rhs_agg = tmp*(rk_wgt[0]*rhs0_sc[ii+isc] + rk_wgt[1]*rhs1_sc[ii+isc] + rk_wgt[2]*rhs2_sc[ii+isc]);
      transport_scalar_vec[isc][icv] = (rho_old[icv]*transport_scalar_vec[isc][icv] + rhs_agg)*inv_rho;
    }
  }
}//advanceScalars()

void IdealGasSolver::setFgrCvs() {

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
  const double inv_pref = 1.0/p_ref;

  for (int icv = 0; icv < ncv; ++icv)
    logp[icv] = log(p[icv]*inv_pref);


  updateCv2DataStart(logp);

  for (int ief = 0; ief < nef_i; ++ief) {

    const int icv0 = cvoef[ief][0];
    const int icv1 = cvoef[ief][1];

    // the gibbs remainder is non-dimensionalized locally with RT_stag

    const double dp       = p[icv1] - p[icv0];
    const double v_avg    = 0.5*(1.0/rho[icv0] + 1.0/rho[icv1]);
    const double beta_avg = 0.5*(1.0/(R_gas*T[icv0]) + 1.0/(R_gas*T[icv1]));
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
    const double beta_avg = 0.5*(1.0/(R_gas*T[icv0]) + 1.0/(R_gas*T[icv1]));
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
  } else if ( checkParam("ALL_DISS")) {
    for (int icv = 0 ; icv < ncv_g2; ++icv)
      cv_compact[icv].fgr = 0.5;
  }

  //delete[] vv2;

  delete[] logp;
}

void IdealGasSolver::updateExtendedState(IdealGasState* __restrict__ cv_light, AuxGasProp* __restrict__ cv_compact,
    double *__restrict__ p, double * __restrict__ T, const double* __restrict__ rho,
    const double (*__restrict__ u)[3], const double* __restrict__ rhoE,
    const int icv_start, const int icv_end) {

  // need to populate the extended state variables
  // inside of the cv_light,cv_compact structures for the
  // cvs owned by this rank..

  const double gm1           = gamma - 1.0;
  const double invgm1        = 1.0/gm1;
  const double Rgogm1        = R_gas*gamma*invgm1;

  for (int icv = icv_start; icv < icv_end; ++icv) {
    double usq = 0.0;
    for (int i =0; i < 3; ++i)
      usq += u[icv][i]*u[icv][i];

    p[icv] = gm1*(rhoE[icv] - 0.5*rho[icv]*usq);
    T[icv] = p[icv]/(rho[icv]*R_gas);
    sos[icv] = sqrt(gamma*R_gas*T[icv]);

    //assert( p[icv] > 0.0);
    //assert( rho[icv] > 0.0);

    // pack the light state..
    for (int i = 0; i < 3; ++i)
      cv_light[icv].u[i] = u[icv][i];
    cv_light[icv].sp_vol = 1.0/rho[icv];
    cv_light[icv].h      = Rgogm1*T[icv];
    cv_light[icv].p      = p[icv];

    // and the rest of the extended state variables..
    cv_compact[icv].half_usq  = 0.5*usq;
    //cv_compact[icv].sor       = (tmp1-gamma*tmp2)*invgm1;
    //cv_compact[icv].invRT     = rho[icv]/p[icv];
    cv_compact[icv].rho       = rho[icv];
    cv_compact[icv].sos       = sos[icv];

  }
}

void IdealGasSolver::computeIgEntropy(double* ent) const {

  // computes sor (non-dimensionalized by the gas constant R)

  const double invgm1        = 1.0/(gamma-1.0);
  const double inv_pref      = 1.0/p_ref;
  const double inv_rhoref    = 1.0/rho_ref;

  for (int icv = 0; icv < ncv; ++icv) {

    const double tmp1 = log1p((p[icv]-p_ref)*inv_pref); // log(p[icv]*inv_pref)
    const double tmp2 = log1p((rho[icv]-rho_ref)*inv_rhoref); // log(rho[icv]*inv_rhoref)
    ent[icv] = (tmp1-gamma*tmp2)*invgm1;

    //ent[icv]   = gamma*invgm1*log(T[icv]/T_ref) - log(p[icv]/p_ref);

  }

}

double IdealGasSolver::calcCfl(double* cfl_, const double dt_) const {

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

double IdealGasSolver::calcDfl(double* cfl_, const double dt_) const {

  bool b_memflag = false;
  double * cfl = NULL;
  if (cfl_ == NULL) {
    cfl = new double[ncv];
    b_memflag = true;
  }
  else
    cfl = cfl_;

  for (int icv =0; icv < ncv; ++icv)
    cfl[icv] = 0.0;

  for (int ibf = 0; ibf < nbf; ++ibf) {

    const int icv = cvobf[ibf];
    const double nu_total = (mu_lam[icv] + mu_sgs[icv])/rho[icv];
    cfl[icv] += area_over_delta_bf[ibf]*nu_total;

  }

  for (int ifa = 0; ifa < nfa; ++ifa) {

    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];

    double nu_total = 0.5*((mu_lam[icv0] + mu_sgs[icv0])/rho[icv0] +
        (mu_lam[icv1] + mu_sgs[icv1])/rho[icv1]);

    cfl[icv0] += area_over_delta_fa[ifa]*nu_total;
    if (icv1 < ncv)
      cfl[icv1] += area_over_delta_fa[ifa]*nu_total;
  }

  double my_cfl_max = 0.0;

  for (int icv = 0; icv < ncv; ++icv) {

    cfl[icv]  *= 0.5*dt_*inv_vol[icv];
    my_cfl_max = max(my_cfl_max, cfl[icv]);

  }

  if (b_memflag) delete[] cfl;
  return my_cfl_max;

}

void IdealGasSolver::preprocessLp(LpSolidRhs * lpSolidRhs_arr[3], LpLiquidRhs * lpLiquidRhs_arr[3], LpSolidRhs * &lpSolidRhs_d, LpLiquidRhs * &lpLiquidRhs_d) {

  if (lpSolid) {
    FOR_I3 {
      lpSolidRhs_arr[i] = new LpSolidRhs [lpSolid->size()];
    }
    if (lpCoupling==ONE_WAY_COUPLED)
      lpSolidRhs_d = new LpSolidRhs [lpSolid->size()];
  }

  if (lpLiquid) {
    FOR_I3 {
      lpLiquidRhs_arr[i] = new LpLiquidRhs [lpLiquid->size()];
    }
    if (lpCoupling==ONE_WAY_COUPLED)
      lpLiquidRhs_d = new LpLiquidRhs [lpLiquid->size()];
  }

  if ((lpSolid)&&(lpLiquid)) {
    FOR_ICV {
      int npSolid = lpSolid->lpocv_i[icv+1] - lpSolid->lpocv_i[icv];
      int npLiquid = lpLiquid->lpocv_i[icv+1] - lpLiquid->lpocv_i[icv];
      npaocv_max = max(npaocv_max, (npSolid + npLiquid));
    }
  }
  else if (lpSolid) {
    FOR_ICV {
      int npSolid = lpSolid->lpocv_i[icv+1] - lpSolid->lpocv_i[icv];
      npaocv_max = max(npaocv_max, npSolid);
    }
  }
  else if (lpLiquid) {
    FOR_ICV {
      int npLiquid = lpLiquid->lpocv_i[icv+1] - lpLiquid->lpocv_i[icv];
      npaocv_max = max(npaocv_max, npLiquid);
    }
  }
  else {
    npaocv_max = 0;
  }
  assert(npaocv_max >= 0);

}

void IdealGasSolver::postprocessLp(LpSolidRhs * lpSolidRhs_arr[3], LpLiquidRhs * lpLiquidRhs_arr[3], LpSolidRhs * lpSolidRhs_d, LpLiquidRhs * lpLiquidRhs_d) {

  if (lpSolid) {
    FOR_I3 {
      delete[] lpSolidRhs_arr[i];
    }
    if (lpSolidRhs_d!=NULL) delete[] lpSolidRhs_d;
  }

  if (lpLiquid) {
    FOR_I3 {
      delete[] lpLiquidRhs_arr[i];
    }
    if (lpLiquidRhs_d!=NULL) delete[] lpLiquidRhs_d;
  }

}

void IdealGasSolver::calcRhsLpSolid_1w(LpSolidRhs * rhs, LpSolidRhs * rhs_d) {

  assert(rhs);
  assert(rhs_d);

  FOR_ICV {

    // compute fixed cv-based props here...
    // gas properties
    const double rhog = rho[icv];
    const double Tg   = T[icv];
    const double mug  = mu_lam[icv];
    const double cpg  = R_gas*gamma/(gamma-1.0);
    const double kg   = loc_lam[icv]*cpg; // loc: lambda(k) over cp
    const double ug[3] = { u[icv][0], u[icv][1], u[icv][2] };

    for (int ip = lpSolid->lpocv_i[icv]; ip != lpSolid->lpocv_i[icv+1]; ++ip) {

      // check sort...
      assert(icv == lpSolid->lp[ip].icv);

      // initialize rhs with zero
      rhs[ip].zero();
      rhs_d[ip].zero();

      //#############################################
      // npar: we compute the Rhs of one droplet,
      //       the effect of npar is taken into acount
      //       in energy transfer and mass transfer to
      //       the gas only
      //#############################################
      // make sure mass is not zero otherwise:
      // set rhs and rhs_d to zero -> no change
      // and trash the particle at the end of time step
      //#############################################
      int keep = getLpKeep(lpSolid->lp[ip].flag);
      if ((keep == TRASH) || (lpSolid->lp[ip].mp <= 0.0)) {
        setLpKeep(lpSolid->lp[ip].flag, TRASH); // make sure it will be trashed
        continue;
      }

      // compute particle properties at parcel temperature
      const int material_id = getLpMaterial(lpSolid->lp[ip].flag);
      assert(material_id>=0); assert(material_id<dustVec.size());
      CtiSolid * dust = dustVec[material_id];
      const double rhop = dust->calcRho(lpSolid->lp[ip].Tp);
      const double Dp   = lpSolid->lp[ip].dp;
      assert(Dp > 0.0);

      double du[3];
      FOR_I3 du[i] = ug[i] - lpSolid->lp[ip].up[i];
      double dus = MAG(du); // slip velocity
      // compute non-dimensional and model parameters
      const double tau = rhop*Dp*Dp/(18.*mug);
      const double Rep = rhog*Dp*dus/mug;
      const double Prg = mug*cpg/kg;
      const double Nu  = 2. + 0.552*sqrt(Rep)*pow(Prg,1./3.);

      // drag law from White, eqn 3-225
      // CD = 24.0/Rep + 6.0/(1.0+sqrt(Rep)) + 0.4 and f1 = Rep/24.0*CD but for numerical reasons when Rep=0 we add a correction to f1
      const double CD_c = 6.0/(1.0+sqrt(Rep)) + 0.4;
      const double f1 =  1.0 + Rep/24.0*CD_c;
      //const double f1 = 1.0; // only for verificatiion step where exact solution needs to be compared

      const double f2 = 1.0;

      ///////////////////////////////
      // setting rhs and rhs_d
      // rhs_total = rhs + rhs_d*var
      ///////////////////////////////
      // xp rhs
      FOR_I3 rhs[ip].xp[i] = lpSolid->lp[ip].up[i];

      // up rhs
      FOR_I3 rhs[ip].up[i] = f1/tau*(ug[i]) + grav[i];
      FOR_I3 rhs_d[ip].up[i] = -f1/tau;

      // Tp rhs
      double const coef = f2*Nu*(cpg/dust->calcCp(lpSolid->lp[ip].Tp))/(3.0*Prg*tau);
      rhs[ip].Tp = coef * (Tg);
      rhs_d[ip].Tp = -coef;

      // monitor if needed
      //if (step%check_interval==0) {
      //  cout << " > ip: "    << ip   << endl <<
      //          " > mug = "  << mug  << endl <<
      //          " > rhog = " << rhog << endl <<
      //          " > ug = " << COUT_VEC(ug) << endl;

      //  cout << " > rhop = " << rhop << endl <<
      //          " > Dp = " << Dp << endl;

      //  cout << " > tau = " << tau << endl <<
      //          " > Rep = " << Rep << endl;

      //  cout << " > f1 = " << f1 << endl;

      //}
    }
  }
}

void IdealGasSolver::calcRhsLpSolid_2w(LpSolidRhs * rhs) {

  FOR_ICV {

    // compute fixed cv-based props here...
    // gas properties
    const double rhog = rho[icv];
    const double mug  = mu_lam[icv];
    const double cpg  = R_gas*gamma/(gamma-1.0);
    const double kg   = loc_lam[icv]*cpg; // loc: lambda(k) over cp
    const double ug[3] = { u[icv][0], u[icv][1], u[icv][2] };

    for (int ip = lpSolid->lpocv_i[icv]; ip != lpSolid->lpocv_i[icv+1]; ++ip) {

      // check sort...
      assert(icv == lpSolid->lp[ip].icv);

      // initialize rhs with zero
      rhs[ip].zero();

      //#############################################
      // npar: we compute the Rhs of one droplet,
      //       the effect of npar is taken into acount
      //       in energy transfer and mass transfer to
      //       the gas only
      //#############################################
      // make sure mass is not zero otherwise:
      // set rhs k and kt to zero -> no change
      // and trash the particle at the end of time step
      //#############################################
      int keep = getLpKeep(lpSolid->lp[ip].flag);
      if ((keep == TRASH) || (lpSolid->lp[ip].mp <= 0.0)) {
        setLpKeep(lpSolid->lp[ip].flag, TRASH); // make sure it will be trashed
        lpSolid->lp[ip].k = 0.0;
        lpSolid->lp[ip].kt = 0.0;
        lpSolid->lp[ip].kb = 0.0;
        lpSolid->lp[ip].ktb = 0.0;
        continue;
      }

      // compute particle properties at parcel temperature
      const int material_id = getLpMaterial(lpSolid->lp[ip].flag);
      assert(material_id>=0); assert(material_id<dustVec.size());
      CtiSolid * dust = dustVec[material_id];
      const double rhop = dust->calcRho(lpSolid->lp[ip].Tp);
      const double Dp   = lpSolid->lp[ip].dp;
      assert(Dp > 0.0);

      double du[3];
      FOR_I3 du[i] = ug[i] - lpSolid->lp[ip].up[i];
      double dus = MAG(du); // slip velocity
      // compute non-dimensional and model parameters
      const double tau = rhop*Dp*Dp/(18.*mug);
      const double Rep = rhog*Dp*dus/mug;
      const double Prg = mug*cpg/kg;
      const double Nu  = 2. + 0.552*sqrt(Rep)*pow(Prg,1./3.);

      // drag law from White, eqn 3-225
      // CD = 24.0/Rep + 6.0/(1.0+sqrt(Rep)) + 0.4 and f1 = Rep/24.0*CD but for numerical reasons when Rep=0 we add a correction to f1
      const double CD_c = 6.0/(1.0+sqrt(Rep)) + 0.4;
      const double f1 =  1.0 + Rep/24.0*CD_c;
      //const double f1 = 1.0; // only for verificatiion step where exact solution needs to be compared

      const double f2 = 1.0;

      // xp rhs
      FOR_I3 rhs[ip].xp[i] = lpSolid->lp[ip].up[i];

      // update the variable k and kt
      lpSolid->lp[ip].k =  f1/tau;
      lpSolid->lp[ip].kt = f2*Nu*(cpg/dust->calcCp(lpSolid->lp[ip].Tp))/(3.0*Prg*tau);
      lpSolid->lp[ip].kb =  0.0;
      lpSolid->lp[ip].ktb = 0.0;

      // up rhs
      FOR_I3 rhs[ip].up[i] = grav[i];

      // Tp rhs, no heat of evaporation
      rhs[ip].Tp = 0.0;

      // monitor if needed
      //if (step%check_interval==0) {
      //  cout << " > ip: "    << ip   << endl <<
      //          " > mug = "  << mug  << endl <<
      //          " > rhog = " << rhog << endl <<
      //          " > ug = " << COUT_VEC(ug) << endl;

      //  cout << " > rhop = " << rhop << endl <<
      //          " > Dp = " << Dp << endl;

      //  cout << " > tau = " << tau << endl <<
      //          " > Rep = " << Rep << endl;

      //  cout << " > f1 = " << f1 << endl;

      //}
      assert(lpSolid->lp[ip].k>=0);
      assert(lpSolid->lp[ip].kt>=0);
    }
  }
}

void IdealGasSolver::calcRhsLpLiquid_1w(LpLiquidRhs * rhs, LpLiquidRhs * rhs_d) {

  assert(rhs);
  assert(rhs_d);

  for (int ip = 0; ip < lpLiquid->size(); ++ip) {
    // initialize rhs with zero
    rhs[ip].zero();
    rhs_d[ip].zero();

    //#############################################
    // npar: we compute the Rhs of one droplet,
    //       the effect of npar is taken into acount
    //       in energy transfer and mass transfer to
    //       the gas only
    //#############################################
    // make sure mass is not zero otherwise:
    // set rhs and rhs_d to zero -> no change
    // and trash the particle at the end of time step
    //#############################################
    int keep = getLpKeep(lpLiquid->lp[ip].flag);
    if ((keep == TRASH) || (lpLiquid->lp[ip].mp <= 0.0)) {
      setLpKeep(lpLiquid->lp[ip].flag, TRASH); // make sure it will be trashed
      continue;
    }

    // gas properties
    const int icv     = lpLiquid->lp[ip].icv;
    const double Tg   = T[icv];
    const double Pg   = p[icv];
    const double rhog = rho[icv];
    const double mug  = mu_lam[icv];
    const double cpg  = R_gas*gamma/(gamma-1.0);
    const double kg   = loc_lam[icv] * cpg; // loc: lambda(k) over cp
    //const double Yfg  = max(0.0,cv[icv].Y);
    const double Yfg  = 0.0;
    const double ug[3]= {u[icv][0], u[icv][1], u[icv][2]};
    const double MWg  = 28.97; // kg/m^3 for air

    // compute liquid properties at parcel temperature
    const int material_id = getLpMaterial(lpLiquid->lp[ip].flag);
    assert(material_id>=0); assert(material_id<fuelVec.size());
    CtiLiquid * fuel = fuelVec[material_id];
    const double Tr   = fuel->Tboil;;
    const double rhof = fuel->calcRho(lpLiquid->lp[ip].Tp);  // liquid density
    const double cpf  = fuel->calcCp(lpLiquid->lp[ip].Tp);   // liquid heat capacity
    const double Lv   = fuel->calcHvap(lpLiquid->lp[ip].Tp); // latent heat of evaporation
    const double Dp   = lpLiquid->lp[ip].dp;
    assert(Dp > 0.0);

    // equilibrium fuel composition at droplet surface
    const double nearOne = 0.9999999999;
    const double Xfseq   = min(nearOne,fuel->calcPv(lpLiquid->lp[ip].Tp)/Pg);
    const double Yfseq   = Xfseq/(Xfseq+(1.0-Xfseq)*MWg/fuel->MW);
    const double BM      = (Yfseq - Yfg)/(1.0 - Yfseq);            // Spalding mass transfer number

    // compute reference prop
    //const double Yfr = (2.*Yfseq + Yfg)/3.;  // 1/3 rule
    const double Yfr = Yfg;
    const double mur = Yfr*fuel->calcMuv(Tr) + (1.-Yfr)*mug;
    const double rhor= Yfr*fuel->calcRhov(Tr) + (1.-Yfr)*rhog;
    const double kr  = Yfr*fuel->calcKv(Tr) + (1.-Yfr)*kg;
    const double cpr = Yfr*fuel->calcCpv(Tr) + (1.-Yfr)*cpg;

    double du[3];
    FOR_I3 du[i] = ug[i] - lpLiquid->lp[ip].up[i];
    double dus = MAG(du);
    // compute non-dimensional and model parameters
    const double tau = rhof*Dp*Dp/(18.*mur);
    const double Rep = rhor*Dp*dus/mur;
    const double Prg = mur*cpr/kr;
    const double Scg = mur/(rhor*fuel->calcDv(Tr));
    const double Nu  = 2. + 0.552*sqrt(Rep)*pow(Prg,1./3.);
    const double Sh  = 2. + 0.552*sqrt(Rep)*pow(Scg,1./3.);
    const double Lk  = mur*sqrt(2.0*M_PI*lpLiquid->lp[ip].Tp*fuel->R_UNIVERSAL/fuel->MW)/(Scg*Pg);

    // Newton-iteration for evaporation rate mpdot
    const double F1 = (Sh/Scg/3.0)*(lpLiquid->lp[ip].mp/tau); // mdot + F1*ln(1+BM) = 0, see eq (4)
    const double F2 = -1.5*Prg*tau/lpLiquid->lp[ip].mp;        // beta = F2*mdot, see eq (17)
    double beta, Xfsneq, Yfsneq, BMneq, Fm;
    double mdot  = -1.0E-08;
    double mdot0 =  0.0;
    double Fm0    =  mdot0 + F1*log(1.0+BM); // initialize with equilib value
    double eps    =  1.0E-15;
    double resid  =  eps+1.0;
    int    iter   =  0;
    while (resid > eps) {
      iter += 1;
      beta   = F2*mdot;
      Xfsneq = Xfseq - (Lk/Dp*2.0)*beta;
      Yfsneq = Xfsneq/(Xfsneq+(1.0-Xfsneq)*MWg/fuel->MW);
      BMneq  = max(-nearOne,(Yfsneq-Yfg)/(1.0-Yfsneq));
      Fm     = mdot + F1*log(1.0+BMneq);
      double tmp = mdot;
      if (fabs(Fm-Fm0) < eps ) break;
      mdot = mdot - Fm/(Fm-Fm0)*(mdot-mdot0);
      resid = sqrt(pow(Fm-Fm0,2));
      Fm0 = Fm;
      mdot0 = tmp;
      if (iter > 20) cout << "ERROR: LSP Newton iteration did not converge !!!" << endl;
    }
    //cout << " > iter = " << iter << endl;
    beta = F2*mdot;

    // Dont let possitive mdot
    mdot = min(mdot, 0.0);
    assert(mdot<=0);

    const double mdot_over_mp = mdot/lpLiquid->lp[ip].mp;

    // drag law from White, eqn 3-225
    // CD = 24.0/Rep + 6.0/(1.0+sqrt(Rep)) + 0.4 and f1 = Rep/24.0*CD but for numerical reasons when Rep=0 we add a correction to f1
    const double CD_c = 6.0/(1.0+sqrt(Rep)) + 0.4;
    const double f1 =  1.0 + Rep/24.0*CD_c;
    double f2 = beta/(exp(beta)-1.0);
    //double f2 = 1.0; // only for non-evaporating case
    if (beta < 1e-16)
      f2 = 1.0;

    ///////////////////////////////
    // setting rhs and rhs_d
    // rhs_total = rhs + rhs_d*var
    ///////////////////////////////
    // xp rhs
    FOR_I3 rhs[ip].xp[i] = lpLiquid->lp[ip].up[i];

    // up rhs
    FOR_I3 rhs[ip].up[i] = f1/tau*(ug[i]);
    FOR_I3 rhs_d[ip].up[i] = -f1/tau;

    // Tp rhs
    const double theta1 = cpr/cpf;
    //rhs[ip].Tp = f2*Nu*theta1/(3.*Prg*tau)*(Tg - lpLiquid->lp[ip].Tp) + Lv/cpf*mdot_over_mp;

    double const coef = f2*Nu*theta1/(3.0*Prg*tau);
    rhs[ip].Tp = coef * (Tg) + Lv/cpf*mdot_over_mp;
    rhs_d[ip].Tp = -coef;

    // mp rhs
    rhs[ip].mp = 0.0;
    rhs_d[ip].mp = mdot_over_mp;

    // Y is the mass fraction of vapor in the gas (vapor mass/ carrier mass)
    // the lspSrc.rhoY is the evaporated vapor mass per time (integrated in vol)
    //
    // YOU SHOULD ADD THIS WHEN YOU HAVE RHOY XXX
    //
    //lspSrc[icv].rhoY += -mdot*lpLiquid->lp[ip].npar;

    // monitor if needed:
    //if (step%check_interval=0 && monitor_bool) {
    //  cout << " > ip: "    << ip   << endl <<
    //          " > Tp: "    << lpLiquid->lp[ip].Tp << endl;;
    //  cout << " > Yfg = "  << Yfg  << endl <<
    //          " > mug = "  << mug  << endl <<
    //          " > rhog = " << rhog << endl <<
    //          " > kg = "   << kg   << endl <<
    //          " > Pg = "   << Pg   << endl <<
    //          " > Tg = "   << Tg   << endl <<
    //          " > cpg = "  << cpg  << endl;

    //  cout << " > rhof = " << rhof << endl <<
    //          " > cpf = " << cpf << endl <<
    //          " > Lv = " << Lv << endl <<
    //          " > Dp = " << Dp << endl;

    //  cout << " > Xfseq = " << Xfseq << endl <<
    //          " > Yfseq = " << Yfseq << endl <<
    //          " > BM = " << BM << endl <<
    //          " > Tr " << Tr << endl;

    //  cout << " > Yfr = " << Yfr << endl <<
    //          " > mur = " << mur << endl <<
    //          " > rhor = " << rhor << endl <<
    //          " > kr = " << kr << endl <<
    //          " > cpr = " << cpr << endl;

    //  cout << " > tau = " << tau << endl <<
    //          " > Rep = " << Rep << endl <<
    //          " > Prg = " << Prg << endl <<
    //          " > Scg = " << Scg << endl <<
    //          " > Nu = " << Nu << endl <<
    //          " > Sh = " << Sh << endl <<
    //          " > Lk = " << Lk << endl;

    //  cout << "mdot: " << mdot << endl;

    //  cout << " > f1 = " << f1 << endl <<
    //          " > f2 = " << f2 << endl <<
    //          " > beta = " << beta << endl;
    //}
  }
}

void IdealGasSolver::calcRhsLpLiquid_2w(LpLiquidRhs * rhs) {

  for (int ip = 0; ip < lpLiquid->size(); ++ip) {

    // initialize rhs with zero
    rhs[ip].zero();

    //#############################################
    // npar: we compute the Rhs of one droplet,
    //       the effect of npar is taken into acount
    //       in energy transfer and mass transfer to
    //       the gas only
    //#############################################
    // make sure mass is not zero otherwise:
    // set rhs k and kt to zero -> no change
    // and trash the particle at the end of time step
    //#############################################
    int keep = getLpKeep(lpLiquid->lp[ip].flag);
    if ((keep == TRASH) || (lpLiquid->lp[ip].mp <= 0.0)) {
      setLpKeep(lpLiquid->lp[ip].flag, TRASH); // make sure it will be trashed
      lpLiquid->lp[ip].k = 0.0;
      lpLiquid->lp[ip].kt = 0.0;
      lpLiquid->lp[ip].kb = 0.0;
      lpLiquid->lp[ip].ktb = 0.0;
      continue;
    }

    // gas properties
    const int icv     = lpLiquid->lp[ip].icv;
    //const double Tg   = T[icv];
    const double Pg   = p[icv];
    const double rhog = rho[icv];
    const double mug  = mu_lam[icv];
    const double cpg  = R_gas*gamma/(gamma-1.0);
    const double kg   = loc_lam[icv] * cpg; // loc: lambda(k) over cp
    //const double Yfg  = max(0.0,cv[icv].Y);
    const double Yfg  = 0.0;
    const double ug[3]= {u[icv][0], u[icv][1], u[icv][2]};
    const double MWg  = 28.97; // kg/m^3 for air

    // compute liquid properties at parcel temperature
    const int material_id = getLpMaterial(lpLiquid->lp[ip].flag);
    assert(material_id>=0); assert(material_id<fuelVec.size());
    CtiLiquid * fuel = fuelVec[material_id];
    const double Tr   = fuel->Tboil;;
    const double rhof = fuel->calcRho(lpLiquid->lp[ip].Tp);  // liquid density
    const double cpf  = fuel->calcCp(lpLiquid->lp[ip].Tp);   // liquid heat capacity
    const double Lv   = fuel->calcHvap(lpLiquid->lp[ip].Tp); // latent heat of evaporation
    const double Dp   = lpLiquid->lp[ip].dp;
    assert(Dp > 0.0);

    // equilibrium fuel composition at droplet surface
    const double nearOne = 0.9999999999;
    const double Xfseq   = min(nearOne,fuel->calcPv(lpLiquid->lp[ip].Tp)/Pg);
    const double Yfseq   = Xfseq/(Xfseq+(1.0-Xfseq)*MWg/fuel->MW);
    const double BM      = (Yfseq - Yfg)/(1.0 - Yfseq);            // Spalding mass transfer number

    // compute reference prop
    //const double Yfr = (2.*Yfseq + Yfg)/3.;  // 1/3 rule
    const double Yfr = Yfg;
    const double mur = Yfr*fuel->calcMuv(Tr) + (1.-Yfr)*mug;
    const double rhor= Yfr*fuel->calcRhov(Tr) + (1.-Yfr)*rhog;
    const double kr  = Yfr*fuel->calcKv(Tr) + (1.-Yfr)*kg;
    const double cpr = Yfr*fuel->calcCpv(Tr) + (1.-Yfr)*cpg;

    double du[3];
    FOR_I3 du[i] = ug[i] - lpLiquid->lp[ip].up[i];
    double dus = MAG(du);
    // compute non-dimensional and model parameters
    const double tau = rhof*Dp*Dp/(18.*mur);
    const double Rep = rhor*Dp*dus/mur;
    const double Prg = mur*cpr/kr;
    const double Scg = mur/(rhor*fuel->calcDv(Tr));
    const double Nu  = 2. + 0.552*sqrt(Rep)*pow(Prg,1./3.);
    const double Sh  = 2. + 0.552*sqrt(Rep)*pow(Scg,1./3.);
    const double Lk  = mur*sqrt(2.0*M_PI*lpLiquid->lp[ip].Tp*fuel->R_UNIVERSAL/fuel->MW)/(Scg*Pg);

    // Newton-iteration for evaporation rate mpdot
    const double F1 = (Sh/Scg/3.0)*(lpLiquid->lp[ip].mp/tau); // mdot + F1*ln(1+BM) = 0, see eq (4)
    const double F2 = -1.5*Prg*tau/lpLiquid->lp[ip].mp;        // beta = F2*mdot, see eq (17)
    double beta, Xfsneq, Yfsneq, BMneq, Fm;
    double mdot  = -1.0E-08;
    double mdot0 =  0.0;
    double Fm0    =  mdot0 + F1*log(1.0+BM); // initialize with equilib value
    double eps    =  1.0E-15;
    double resid  =  eps+1.0;
    int    iter   =  0;
    while (resid > eps) {
      iter += 1;
      beta   = F2*mdot;
      Xfsneq = Xfseq - (Lk/Dp*2.0)*beta;
      Yfsneq = Xfsneq/(Xfsneq+(1.0-Xfsneq)*MWg/fuel->MW);
      BMneq  = max(-nearOne,(Yfsneq-Yfg)/(1.0-Yfsneq));
      Fm     = mdot + F1*log(1.0+BMneq);
      double tmp = mdot;
      if (fabs(Fm-Fm0) < eps ) break;
      mdot = mdot - Fm/(Fm-Fm0)*(mdot-mdot0);
      resid = sqrt(pow(Fm-Fm0,2));
      Fm0 = Fm;
      mdot0 = tmp;
      if (iter > 20) cout << "ERROR: LSP Newton iteration did not converge !!!" << endl;
    }
    //cout << " > iter = " << iter << endl;
    beta = F2*mdot;

    // Dont let possitive mdot
    mdot = min(mdot, 0.0);
    assert(mdot<=0);

    const double mdot_over_mp = mdot/lpLiquid->lp[ip].mp;

    // drag law from White, eqn 3-225
    // CD = 24.0/Rep + 6.0/(1.0+sqrt(Rep)) + 0.4 and f1 = Rep/24.0*CD but for numerical reasons when Rep=0 we add a correction to f1
    const double CD_c = 6.0/(1.0+sqrt(Rep)) + 0.4;
    const double f1 =  1.0 + Rep/24.0*CD_c;
    double f2 = beta/(exp(beta)-1.0);
    //double f2 = 1.0; // only for non-evaporating case
    if (beta < 1e-16)
      f2 = 1.0;

    // xp rhs
    FOR_I3 rhs[ip].xp[i] = lpLiquid->lp[ip].up[i];

    lpLiquid->lp[ip].k = f1/tau;
    lpLiquid->lp[ip].kt = f2*Nu*(cpr/cpf)/(3.0*Prg*tau);
    lpLiquid->lp[ip].kb =  0.0;
    lpLiquid->lp[ip].ktb = 0.0;

    assert(lpLiquid->lp[ip].k>0);
    assert(lpLiquid->lp[ip].kt>0);

    // up rhs
    FOR_I3 rhs[ip].up[i] = grav[i];

    // Tp rhs
    //const double theta1 = cpr/cpf;
    //rhs[ip].Tp = f2*Nu*theta1/(3.*Prg*tau)*(Tg - P[ip].Tp) + Lv/cpf*mdot/P[ip].mp;
    //rhs[ip].Tp = f2*Nu*theta1/(3.*Prg*tau)*(Tg) + Lv/cpf*mdot_over_mp;
    //rhs_d[ip].Tp = - f2*Nu*theta1/(3.*Prg*tau);

    //#############################################
    // only evaporation latent heat is treated explicitly
    //#############################################
    rhs[ip].Tp = Lv/cpf*mdot_over_mp;

    //#############################################
    // if you like to test a zero evaporation case, make sure
    // in addition to zeroing mass and temperature source terms
    // you make f2=1. This is the limit for non-evaporative case
    // otherwise your kt coeff is not correct
    //#############################################

    // mp rhs
    // mdot/mp is stored here. This is for implicit treatment.
    rhs[ip].mp = mdot_over_mp;

    // Y is the mass fraction of vapor in the gas (vapor mass/ carrier mass)
    // the lspSrc.rhoY is the evaporated vapor mass per time (integrated in vol)
    //
    // YOU SHOULD ADD THIS WHEN YOU HAVE RHOY XXX
    //
    //lspSrc[icv].rhoY += -mdot*lpLiquid->lp[ip].npar;

    // monitor if needed:
    //if (step%check_interval=0 && monitor_bool) {
    //  cout << " > ip: "    << ip   << endl <<
    //          " > Tp: "    << lpLiquid->lp[ip].Tp << endl;;
    //  cout << " > Yfg = "  << Yfg  << endl <<
    //          " > mug = "  << mug  << endl <<
    //          " > rhog = " << rhog << endl <<
    //          " > kg = "   << kg   << endl <<
    //          " > Pg = "   << Pg   << endl <<
    //          " > Tg = "   << Tg   << endl <<
    //          " > cpg = "  << cpg  << endl;

    //  cout << " > rhof = " << rhof << endl <<
    //          " > cpf = " << cpf << endl <<
    //          " > Lv = " << Lv << endl <<
    //          " > Dp = " << Dp << endl;

    //  cout << " > Xfseq = " << Xfseq << endl <<
    //          " > Yfseq = " << Yfseq << endl <<
    //          " > BM = " << BM << endl <<
    //          " > Tr " << Tr << endl;

    //  cout << " > Yfr = " << Yfr << endl <<
    //          " > mur = " << mur << endl <<
    //          " > rhor = " << rhor << endl <<
    //          " > kr = " << kr << endl <<
    //          " > cpr = " << cpr << endl;

    //  cout << " > tau = " << tau << endl <<
    //          " > Rep = " << Rep << endl <<
    //          " > Prg = " << Prg << endl <<
    //          " > Scg = " << Scg << endl <<
    //          " > Nu = " << Nu << endl <<
    //          " > Sh = " << Sh << endl <<
    //          " > Lk = " << Lk << endl;

    //  cout << "mdot: " << mdot << endl;

    //  cout << " > f1 = " << f1 << endl <<
    //          " > f2 = " << f2 << endl <<
    //          " > beta = " << beta << endl;
    //}
  }
}

void IdealGasSolver::rk3Step_lp_gas(const double rk_wgt[3], IdealGasRhs ** rhs_arr, LpSolidRhs ** lpSolidRhs_arr, LpLiquidRhs ** lpLiquidRhs_arr, const int rkstep) {

  // we already know the maximum number of parcels per cv, npaocv_max, so
  // allocate memory once...

  // below we build the system
  //
  // [ A_d v ]{ up,ug } = {rhs}
  // [ q   D ]
  //
  // [ A_d v ]{ Tp,rhoE } = {rhs}
  // [ q   D ]

  LpSolidState * lpS  = lpSolid ? lpSolid->lp:NULL;
  LpLiquidState * lpL = lpLiquid ? lpLiquid->lp:NULL;

  assert(npaocv_max>=0);
  double * A_d = new double[npaocv_max];
  double * v   = new double[npaocv_max];
  double * q   = new double[npaocv_max];
  double * rhs[3];
  double * X[3];
  for (int i = 0; i < 3; ++i) {
    rhs[i] = new double [npaocv_max+1];
    X[i] = new double [npaocv_max+1];
  }

  const double cp = R_gas*gamma/(gamma-1.0);
  IdealGasRhs rhs_agg;
  FOR_ICV {
    rhs_agg.zero();
    const double tmp = dt/vol_cv[icv];
    for (int irk = 0; irk < rkstep; ++irk) {
      const double wgt = tmp*rk_wgt[irk];
      rhs_agg.rho += wgt*rhs_arr[irk][icv].rho;
      FOR_I3 rhs_agg.rhou[i] += wgt*rhs_arr[irk][icv].rhou[i];
      rhs_agg.rhoE += wgt*rhs_arr[irk][icv].rhoE;
      //rhs_agg.axpy(wgt,rhs_arr[irk][icv]);
    }
    const double rho_old = rho[icv];
    rho[icv] += rhs_agg.rho;
    // SHOULD ADD LATER XXX
    //cv[icv].Y = (rho_old*cv[icv].Y + rhs_agg.rhoY) / cv[icv].rho;

    // now solve the gas and particle velocity together
    const int np_Sol = (lpSolid) ? (lpSolid->lpocv_i[icv+1]-lpSolid->lpocv_i[icv]):0;
    const int np_Liq = (lpLiquid) ? (lpLiquid->lpocv_i[icv+1]-lpLiquid->lpocv_i[icv]):0;
    const int np_SolLiq = np_Sol + np_Liq;
    if (np_SolLiq == 0) {
      // solve gas velocity/energy with regular rk3
      FOR_I3 u[icv][i] = (u[icv][i]*rho_old + rhs_agg.rhou[i])/rho[icv];
      rhoE[icv] += rhs_agg.rhoE;
    }
    else {
      assert((np_SolLiq > 0)&&(np_SolLiq <= npaocv_max));
      ////////////////////////////////////////////////
      // Solve the coupled particles/gas velocities
      ////////////////////////////////////////////////
      double D = 0.0;
      if (lpSolid) {
        for (int ip = lpSolid->lpocv_i[icv]; ip != lpSolid->lpocv_i[icv+1]; ++ip) {
          const int ii = ip - lpSolid->lpocv_i[icv];
          assert((ii >= 0)&&(ii < np_SolLiq));
          assert((ip >= 0)&&(ip < lpSolid->size()));
          assert(rk_wgt[rkstep-1] > 0.0);
          assert(lpS[ip].k >= 0.0);
          const double ki = lpS[ip].k * rk_wgt[rkstep-1];
          assert(lpS[ip].kb >= 0.0);
          const double kbi = lpS[ip].kb * rk_wgt[rkstep-1];
          A_d[ii] = ki + 1.0/dt + kbi;
          v[ii] = -ki;
          q[ii] = -ki*lpS[ip].mp;
          FOR_I3 rhs[i][ii] = lpS[ip].up[i]/dt;
          for (int irk = 0; irk < rkstep; ++irk) {
            FOR_I3 rhs[i][ii] += lpSolidRhs_arr[irk][ip].up[i]*rk_wgt[irk];
          }
          D += -q[ii];
        }
      }
      if (lpLiquid) {
        for (int ip = lpLiquid->lpocv_i[icv]; ip != lpLiquid->lpocv_i[icv+1]; ++ip) {
          const int ii = np_Sol + ip - lpLiquid->lpocv_i[icv];
          assert((ii >= 0)&&(ii < np_SolLiq));
          assert((ip >= 0)&&(ip < lpLiquid->size()));
          assert(rk_wgt[rkstep-1] > 0.0);
          assert(lpL[ip].k >= 0.0);
          const double ki = lpL[ip].k * rk_wgt[rkstep-1];
          assert(lpL[ip].kb >= 0.0);
          const double kbi = lpL[ip].kb * rk_wgt[rkstep-1];
          A_d[ii] = ki + 1.0/dt + kbi;
          v[ii] = -ki;
          q[ii] = -ki*lpL[ip].mp;
          FOR_I3 rhs[i][ii] = lpL[ip].up[i]/dt;
          for (int irk = 0; irk < rkstep; ++irk) {
            FOR_I3 rhs[i][ii] += lpLiquidRhs_arr[irk][ip].up[i]*rk_wgt[irk];
          }
          D += -q[ii];
        }
      }
      D += rho[icv]/tmp;
      FOR_I3 rhs[i][np_SolLiq] = rhs_agg.rhou[i]/tmp + rho_old*u[icv][i]/tmp;
      solveLinearSysGE(X, A_d, v, q, D, rhs, np_SolLiq);
      //checkSolution(X, np);
      //MPI_Pause("check the solution");
      ///////////////////////////////////////
      // update gas and particle velocities
      ///////////////////////////////////////
      FOR_I3 u[icv][i] = X[i][np_SolLiq];
      if (lpSolid) {
        for (int ip = lpSolid->lpocv_i[icv]; ip != lpSolid->lpocv_i[icv+1]; ++ip) {
          const int ii = ip - lpSolid->lpocv_i[icv];
          assert(ip < lpSolid->size());
          FOR_I3 lpS[ip].up[i] = X[i][ii];

          // no rk_wgt multiply to these rhs corrections:
          // update the rhs of gas rhou for last iteration
          // NOTE: rk_wgt is not multiplied
          FOR_I3 rhs_arr[rkstep-1][icv].rhou[i] += lpS[ip].k*lpS[ip].mp*(lpS[ip].up[i] - u[icv][i]);
          // update the rhs of particle up for last iteration
          FOR_I3 lpSolidRhs_arr[rkstep-1][ip].up[i] += lpS[ip].k*(u[icv][i] - lpS[ip].up[i]);
          // add the boundary piece: kb * (0-up) , this is the implicit part,
          FOR_I3 lpSolidRhs_arr[rkstep-1][ip].up[i] += lpS[ip].kb * (-lpS[ip].up[i]);
        }
      }
      if (lpLiquid) {
        for (int ip = lpLiquid->lpocv_i[icv]; ip != lpLiquid->lpocv_i[icv+1]; ++ip) {
          const int ii = np_Sol + ip - lpLiquid->lpocv_i[icv];
          assert(ip < lpLiquid->size());
          FOR_I3 lpL[ip].up[i] = X[i][ii];

          // no rk_wgt multiply to these rhs corrections:
          // update the rhs of gas rhou for last iteration
          // NOTE: rk_wgt is not multiplied
          FOR_I3 rhs_arr[rkstep-1][icv].rhou[i] += lpL[ip].k*lpL[ip].mp*(lpL[ip].up[i] - u[icv][i]);
          // update the rhs of particle up for last iteration
          FOR_I3 lpLiquidRhs_arr[rkstep-1][ip].up[i] += lpL[ip].k*(u[icv][i] - lpL[ip].up[i]);
          // add the boundary piece: kb * (0-up) , this is the implicit part,
          FOR_I3 lpLiquidRhs_arr[rkstep-1][ip].up[i] += lpL[ip].kb * (-lpL[ip].up[i]);
        }
      }
      ////////////////////////////////////////////////
      // Solve the coupled particle/gas energies
      ////////////////////////////////////////////////
      D = 0.0;
      rhs[0][np_SolLiq] = 0.0;
      const double ug2mag = DOT_PRODUCT(u[icv], u[icv]);
      if (lpSolid) {
        for (int ip = lpSolid->lpocv_i[icv]; ip != lpSolid->lpocv_i[icv+1]; ++ip) {
          const int ii = ip - lpSolid->lpocv_i[icv];
          assert(ii >= 0);  assert(ii < np_SolLiq);
          assert(ip < lpSolid->size());
          assert(rk_wgt[rkstep-1] > 0.0);
          assert(lpS[ip].kt >= 0.0);
          const double ki = lpS[ip].kt * rk_wgt[rkstep-1];
          assert(lpS[ip].ktb >= 0.0);
          const double kbi = lpS[ip].ktb * rk_wgt[rkstep-1];
          A_d[ii] = ki + 1.0/dt + kbi;
          v[ii] = -ki/(rho[icv]*cp/gamma);
          //q[ii] = -ki*(fuel->calcRho(lspVec[ip].Tp) * fuel->calcCp(lspVec[ip].Tp)) * cv[icv].vol * lspVec[ip].npar;
          const int material_id = getLpMaterial(lpS[ip].flag);
          assert(material_id>=0); assert(material_id<dustVec.size());
          CtiSolid * dust = dustVec[material_id];
          q[ii] = -ki*(dust->calcRho(lpS[ip].Tp) * dust->calcCp(lpS[ip].Tp)) * (M_PI/6.*pow(lpS[ip].dp,3)) * lpS[ip].npar;
          rhs[0][ii] = lpS[ip].Tp/dt - ki * (0.5 * ug2mag) / (cp/gamma);
          rhs[0][np_SolLiq] += -q[ii] * (0.5 * ug2mag) / (cp/gamma);
          for (int irk = 0; irk < rkstep; ++irk) {
            rhs[0][ii] += lpSolidRhs_arr[irk][ip].Tp*rk_wgt[irk];
          }
          D += -q[ii] / (rho[icv] * cp/gamma);
        }
      }
      if (lpLiquid) {
        for (int ip = lpLiquid->lpocv_i[icv]; ip != lpLiquid->lpocv_i[icv+1]; ++ip) {
          const int ii = np_Sol + ip - lpLiquid->lpocv_i[icv];
          assert(ii >= 0);  assert(ii < np_SolLiq);
          assert(ip < lpLiquid->size());
          assert(rk_wgt[rkstep-1] > 0.0);
          assert(lpL[ip].kt >= 0.0);
          const double ki = lpL[ip].kt * rk_wgt[rkstep-1];
          assert(lpL[ip].ktb >= 0.0);
          const double kbi = lpL[ip].ktb * rk_wgt[rkstep-1];
          A_d[ii] = ki + 1.0/dt + kbi;
          v[ii] = -ki/(rho[icv]*cp/gamma);
          //q[ii] = -ki*(fuel->calcRho(lspVec[ip].Tp) * fuel->calcCp(lspVec[ip].Tp)) * cv[icv].vol * lspVec[ip].npar;
          const int material_id = getLpMaterial(lpL[ip].flag);
          assert(material_id>=0); assert(material_id<fuelVec.size());
          CtiLiquid * fuel = fuelVec[material_id];
          q[ii] = -ki*(fuel->calcRho(lpL[ip].Tp) * fuel->calcCp(lpL[ip].Tp)) * (M_PI/6.*pow(lpL[ip].dp,3)) * lpL[ip].npar;
          rhs[0][ii] = lpL[ip].Tp/dt - ki * (0.5 * ug2mag) / (cp/gamma);
          rhs[0][np_SolLiq] += -q[ii] * (0.5 * ug2mag) / (cp/gamma);
          for (int irk = 0; irk < rkstep; ++irk) {
            rhs[0][ii] += lpLiquidRhs_arr[irk][ip].Tp*rk_wgt[irk];
          }
          D += -q[ii] / (rho[icv] * cp/gamma);
        }
      }
      D += 1.0/tmp;
      rhs[0][np_SolLiq] += rhs_agg.rhoE/tmp + rhoE[icv]/tmp;
      solveLinearSysGE(X[0], A_d, v, q, D, rhs[0], np_SolLiq);
      //checkSolution(X[0], np);
      ///////////////////////////////////////
      // update gas and particle temp/energy
      ///////////////////////////////////////
      rhoE[icv] = X[0][np_SolLiq];
      const double Tg = (rhoE[icv] - 0.5*rho[icv]*ug2mag)/(rho[icv]*cp/gamma);
      if (lpSolid) {
        for (int ip = lpSolid->lpocv_i[icv]; ip != lpSolid->lpocv_i[icv+1]; ++ip) {
          const int ii = ip - lpSolid->lpocv_i[icv];
          assert(ip < lpSolid->size());
          lpS[ip].Tp = X[0][ii];
          // update the rhs of gas rhoE for last iteration, NOTE: no rk_wgt should be added here
          rhs_arr[rkstep-1][icv].rhoE +=  (-q[ii]/rk_wgt[rkstep-1]) * (lpS[ip].Tp - Tg); // npar is added in q above, devide by rk_wgt to elliminate it
          // update the rhs of particle Tp for last iteration (no need to mult npar cause the equation is normalized with mass)
          lpSolidRhs_arr[rkstep-1][ip].Tp += lpS[ip].kt * ( Tg - lpS[ip].Tp);
          // add the boundary piece of the rk_wgt * ktb * (-Tp) , this is the implicit part, note that the Ts part is added already in calcRhsLsp
          lpSolidRhs_arr[rkstep-1][ip].Tp += lpS[ip].ktb * (-lpS[ip].Tp);
        }
      }
      if (lpLiquid) {
        for (int ip = lpLiquid->lpocv_i[icv]; ip != lpLiquid->lpocv_i[icv+1]; ++ip) {
          const int ii = np_Sol + ip - lpLiquid->lpocv_i[icv];
          assert(ip < lpLiquid->size());
          lpL[ip].Tp = X[0][ii];
          // update the rhs of gas rhoE for last iteration, NOTE: no rk_wgt should be added here
          rhs_arr[rkstep-1][icv].rhoE +=  (-q[ii]/rk_wgt[rkstep-1]) * (lpL[ip].Tp - Tg); // npar is added in q above, devide by rk_wgt to elliminate it
          // update the rhs of particle Tp for last iteration (no need to mult npar cause the equation is normalized with mass)
          lpLiquidRhs_arr[rkstep-1][ip].Tp += lpL[ip].kt * ( Tg - lpL[ip].Tp);
          // add the boundary piece of the rk_wgt * ktb * (-Tp) , this is the implicit part, note that the Ts part is added already in calcRhsLsp
          lpLiquidRhs_arr[rkstep-1][ip].Tp += lpL[ip].ktb * (-lpL[ip].Tp);
        }
      }
    }
  }

  // cleanup matrix memory...
  delete[] A_d;
  delete[] v;
  delete[] q;
  for (int i = 0; i < 3; ++i) {
    delete[] rhs[i];
    delete[] X[i];
  }

  //cout << "position: " << i_pos++ << endl;
  // update other particle properties, xp, mp...
  if (lpSolid) {
    double xp_rhs_agg[3];
    for (int ip = 0; ip < lpSolid->size(); ++ip) {
      FOR_I3 xp_rhs_agg[i] = 0.0;
      for (int irk = 0; irk < rkstep; ++irk) {
        const double wgt = dt*rk_wgt[irk];
        FOR_I3 xp_rhs_agg[i] += wgt * lpSolidRhs_arr[irk][ip].xp[i];
      }
      FOR_I3 lpS[ip].xp[i] += xp_rhs_agg[i];
      // make sure mass is possitive
      assert(lpS[ip].mp > 0);
    }
  }
  if (lpLiquid) {
    double xp_rhs_agg[3];
    for (int ip = 0; ip < lpLiquid->size(); ++ip) {
      FOR_I3 xp_rhs_agg[i] = 0.0;
      for (int irk = 0; irk < rkstep; ++irk) {
        const double wgt = dt*rk_wgt[irk];
        FOR_I3 xp_rhs_agg[i] += wgt * lpLiquidRhs_arr[irk][ip].xp[i];
      }
      const double tmp = dt*rk_wgt[rkstep-1];
      FOR_I3 lpL[ip].xp[i] += xp_rhs_agg[i];
      // the value of rhs for mp is mdot/mp. So add it implicitly for stability...
      // this value is not needed for any other iteration and it's overwritten in each substep
      lpL[ip].mp = (lpL[ip].mp)/(1.0 - lpLiquidRhs_arr[rkstep-1][ip].mp*tmp);
      // make sure mass is possitive
      lpL[ip].mp = max(0.0,lpL[ip].mp);
      // trash particles with zero mass
      if (lpL[ip].mp == 0.0) setLpKeep(lpL[ip].flag, TRASH);
    }
  }

}

// These solvers implement Gauss Elimination that leverages the
// diagonal structure of the matrix to minimize flop's...
inline void IdealGasSolver::solveLinearSysGE(double* X[3],const double* const A_d, const double* const v, const double* const q, const double D, double* const rhs[3], const int np) const {
  double Di = D;
  for (int ip = 0; ip<np; ++ip) {
    const double mult = -q[ip]/A_d[ip];
    Di += mult * v[ip];
    FOR_I3 rhs[i][np] += mult * rhs[i][ip];
  }
  FOR_I3 X[i][np] = rhs[i][np]/Di;
  for (int ip = 0; ip<np; ++ip) {
    FOR_I3 X[i][ip] = (rhs[i][ip] - v[ip]*X[i][np]) / A_d[ip];
  }
}

inline void IdealGasSolver::solveLinearSysGE(double* X, const double* const A_d, const double * const v, const double* const q, const double D,double* const rhs, const int np) const {
  double Di = D;
  for (int ip = 0; ip<np; ++ip) {
    const double mult = -q[ip]/A_d[ip];
    Di += mult * v[ip];
    rhs[np] += mult * rhs[ip];
  }
  X[np] = rhs[np]/Di;
  for (int ip = 0; ip<np; ++ip) {
    X[ip] = (rhs[ip] - v[ip]*X[np]) / A_d[ip];
  }
}

void IdealGasSolver::rk3Step_lpSolid(const double rk_wgt[3], LpSolidRhs ** lpSolidRhs_arr, LpSolidRhs * rhs_d, const int rkstep) {
  for (int ip=0; ip < lpSolid->size(); ++ip) {
    // build the new state
    LpSolidRhs rhs_agg;
    rhs_agg.zero();
    for (int irk=0; irk < rkstep; ++irk) {
      double wgt = dt*rk_wgt[irk];
      rhs_agg.axpy(wgt,lpSolidRhs_arr[irk][ip]);
    }
    const double temp = dt*rk_wgt[rkstep-1];
    lpSolid->lp[ip].add_im(rhs_agg,rhs_d[ip],temp);
    // correct the rhs of previous state
    FOR_I3 lpSolidRhs_arr[rkstep-1][ip].xp[i] += rhs_d[ip].xp[i]*lpSolid->lp[ip].xp[i];
    FOR_I3 lpSolidRhs_arr[rkstep-1][ip].up[i] += rhs_d[ip].up[i]*lpSolid->lp[ip].up[i];
    lpSolidRhs_arr[rkstep-1][ip].Tp += rhs_d[ip].Tp*lpSolid->lp[ip].Tp;
  }
}

void IdealGasSolver::rk3Step_lpLiquid(const double rk_wgt[3], LpLiquidRhs ** lpLiquidRhs_arr, LpLiquidRhs * rhs_d, const int rkstep) {
  for (int ip=0; ip < lpLiquid->size(); ++ip) {
    // build the new state
    LpLiquidRhs rhs_agg;
    rhs_agg.zero();
    for (int irk=0; irk < rkstep; ++irk) {
      double wgt = dt*rk_wgt[irk];
      rhs_agg.axpy(wgt,lpLiquidRhs_arr[irk][ip]);
    }
    const double temp = dt*rk_wgt[rkstep-1];
    lpLiquid->lp[ip].add_im(rhs_agg,rhs_d[ip],temp);
    // correct the rhs of previous state
    FOR_I3 lpLiquidRhs_arr[rkstep-1][ip].xp[i] += rhs_d[ip].xp[i]*lpLiquid->lp[ip].xp[i];
    FOR_I3 lpLiquidRhs_arr[rkstep-1][ip].up[i] += rhs_d[ip].up[i]*lpLiquid->lp[ip].up[i];
    lpLiquidRhs_arr[rkstep-1][ip].Tp += rhs_d[ip].Tp*lpLiquid->lp[ip].Tp;
    lpLiquidRhs_arr[rkstep-1][ip].mp += rhs_d[ip].mp*lpLiquid->lp[ip].mp;
  }
}

void IdealGasSolver::addBodyForceSourceTerms(IdealGasRhs* __restrict__ rhs, double *__restrict__ rhs_sc,
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
          CERR("ADD_BODY_FORCE ... TYPE=MRF ... is not supported by the ideal gas equation of state");
        }
        else if (body_force->type == HOOK_TYPE) {
          CERR("ADD_BODY_FORCE ... TYPE=HOOK ... is not supported by the ideal gas equation of state");
        }
        else {
          assert(0);  //logic problem, shouldn't get here
        }
      }
    }
  }
}



#undef FOR_BCZONE
