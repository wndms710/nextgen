
#include "HelmholtzSolver.hpp"

// convenience macro to loop over the boundary conditions...
#define FOR_BCZONE for(vector<HelmholtzBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)

HelmholtzSolver::HelmholtzSolver(const int icg) : FlowSolver(icg) {

  b_helmholtzCorrection = checkParam("HELMHOLTZ_CORRECTION");

  // register data with i/o

  rho   = NULL; 
  rho0  = NULL; 
  rho00 = NULL; 
  u     = NULL; 
  u0    = NULL; 
  u00   = NULL; 
  p     = NULL; 
  dpdx  = NULL; 

  mf    = NULL; 
  mf0   = NULL; 
  mf00  = NULL; 

  // data registration without i/o
  mu_lam  = NULL; 
  mu_sgs  = NULL; 
  a_sgs  = NULL; 
  a_sgs_sponge  = NULL; 
  rho_uncorrected = NULL; 

  //rhs = NULL;
  rhou_source = NULL; 
  A_diag_source = NULL; 
  b_bodyForceHook = false;

  planar_vol_ratio_cv = NULL;

  // no coarse grids/solvers on lower levels
  ncg_u = 0;
  ncg_p = 0;

  helmholtz_sos = getDoubleParam("HELMHOLTZ_SOS");
  
  frame_rotation = NULL;

  serial_state_interval = -1;
}

HelmholtzSolver::~HelmholtzSolver() {

  DELETE(rho);
  DELETE(rho0);
  DELETE(rho00);
  DELETE(u);
  DELETE(u0);
  DELETE(u00);
  DELETE(p);
  DELETE(dpdx);
  DELETE(mf);
  DELETE(mf0);
  DELETE(mf00);
  DELETE(mu_lam);
  DELETE(mu_sgs);
  DELETE(a_sgs);
  DELETE(a_sgs_sponge);
  DELETE(rho_uncorrected);

  FOR_BCZONE delete *it;

  DELETE(rhou_source);
  DELETE(A_diag_source);

  // coarse grid/data

  for (int ii = 0, lim = cgs.size(); ii < lim; ++ii)
    delete cgs[ii];
  for (int ii = 0, lim = css.size(); ii < lim; ++ii)
    delete css[ii];

  DELETE(planar_vol_ratio_cv);

  DELETE(frame_rotation);

}

void HelmholtzSolver::initData() {

  FlowSolver::initData();

  assert(rho   == NULL); rho     = new double[ncv_g];
  assert(rho0  == NULL); rho0    = new double[ncv_g];
  assert(rho00 == NULL); rho00   = new double[ncv_g];
  assert(u     == NULL); u       = new double[ncv_g][3];
  assert(u0    == NULL); u0      = new double[ncv_g][3];
  assert(u00   == NULL); u00     = new double[ncv_g][3];
  assert(p     == NULL); p       = new double[ncv_g];
  assert(dpdx  == NULL); dpdx    = new double[ncv][3];
  
  assert(mf    == NULL); mf      = new double[nfa];
  if (!b_helmholtzCorrection){
    assert(mf0   == NULL); mf0     = new double[nfa];
    assert(mf00  == NULL); mf00    = new double[nfa];
  }

  // allocate buffers required for the scalar transport

  initScalars();
    
  // data registration without i/o

  assert(mu_lam  == NULL); mu_lam = new double[ncv_g]; 
  assert(mu_sgs  == NULL); mu_sgs = new double[ncv_g];
  assert(a_sgs   == NULL);  a_sgs = new double[ncv_g];
  assert(a_sgs_sponge  == NULL);  a_sgs_sponge = new double[ncv_g];

  assert(rho_uncorrected == NULL); rho_uncorrected = new double[ncv_g];

  if ((max(ncg_u,ncg_p) > 0)&&(mpi_rank == 0))
    cout << " > initializing " << max(ncg_u,ncg_p) << " multigrid levels." << endl;
  for (int icg = 0; icg < max(ncg_u,ncg_p); ++icg) {
    CoarseGrid* cg;
    if (icg == 0) {
      cg = new CoarseGrid(this);
    }
    else {
      cg = new CoarseGrid(css[icg-1]);
    }
    cgs.push_back(cg);
    cg->initCoarseGrid(1,agglomeration_factor,split_orphaned_colors);
    HelmholtzSolver* cs = getHelmholtzSolverMg(icg); 
    cs->initFromCoarseGrid(cg);
    css.push_back(cs);
  }

  assert(planar_vol_ratio_cv == NULL); planar_vol_ratio_cv = new double[ncv_g];

}

//may be overridden by a user defined BC using
//multigrid
HelmholtzSolver* HelmholtzSolver::getHelmholtzSolverMg(const int icg){
  return(new HelmholtzSolver(icg));
}

void HelmholtzSolver::initVVOps() {

  // these are calculated from vida vor's methodology...

  FOR_ICV_G FOR_I3 x_cv[icv][i] = x_vv[icv][i];

  const double area_eps = 1.0E-20;
  int my_count = 0;
  FOR_IBF {
    FOR_I3 x_bf[ibf][i] = 0.0;
    FOR_I3 n_bf[ibf][i] = 0.0;

    // loop to get the normal and approximate x_bf...
    double x_bf_tmp[3] = { 0.0, 0.0, 0.0 };
    const int ino_f = noobf_v[noobf_i[ibf]];
    double wgt_sum = 0.0;
    int ino1 = noobf_v[noobf_i[ibf+1]-1];
    for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
      const int ino0 = ino1;
      ino1 = noobf_v[nob];
      const double wgt = DIST(x_no[ino0],x_no[ino1]);
      FOR_I3 x_bf_tmp[i] += wgt*(x_no[ino0][i]+x_no[ino1][i]);
      wgt_sum += wgt;
      const double n[3] = TRI_NORMAL_2(x_no[ino_f],x_no[ino0],x_no[ino1]);
      FOR_I3 n_bf[ibf][i] += 0.5*n[i];
    }
    assert(wgt_sum > 0.0);
    FOR_I3 x_bf_tmp[i] /= wgt_sum*2.0;

    // and get the area of this face...
    area_bf[ibf] = MAG(n_bf[ibf]);

    // get exact x_bf...
    wgt_sum = 0.0;
    ino1 = noobf_v[noobf_i[ibf+1]-1];
    for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
      const int ino0 = ino1;
      ino1 = noobf_v[nob];
      const double n[3] = TRI_NORMAL_2(x_bf_tmp,x_no[ino0],x_no[ino1]);
      const double wgt = DOT_PRODUCT(n_bf[ibf],n); // area squared weighting !?
      FOR_I3 x_bf[ibf][i] += wgt*(x_bf_tmp[i]+x_no[ino0][i]+x_no[ino1][i]);
      wgt_sum += wgt;
    }
    if (wgt_sum > 0.0) {
      FOR_I3 x_bf[ibf][i] /= wgt_sum*3.0;
    }
    else {
      // this must be a zero-area face -- so just pick a node...
      FOR_I3 x_bf[ibf][i] = x_no[ino1][i];
    }
    const int icv = cvobf[ibf];
    const double dx[3] = DIFF(x_bf[ibf],x_vv[icv]);
    const double unit_n[3] = {n_bf[ibf][0]/max(area_bf[ibf],area_eps),n_bf[ibf][1]/max(area_bf[ibf],area_eps),n_bf[ibf][2]/max(area_bf[ibf],area_eps)};
    const double dp = DOT_PRODUCT(dx,unit_n);
    const double dx_vol = pow(vol_cv[icv],1.0/3.0);
    const double dx_vol_area = vol_cv[icv]/max(area_bf[ibf],area_eps);
    const double dx_clip = min( 0.25*dx_vol , 0.25*dx_vol_area );
    if (dp < dx_clip) {
      ++my_count;
      area_over_delta_bf[ibf] = area_bf[ibf]/dx_clip;
    }
    else {
      area_over_delta_bf[ibf] = area_bf[ibf]/dp;
    }
  }

  int count;
  MPI_Reduce(&my_count,&count,1,MPI_INT,MPI_SUM,0,mpi_comm);
  if (mpi_rank == 0)
    cout << " > boundary length scale clipped in " << count << " faces." << endl;

  // recompute face geometry from nodes

  FOR_IFA {
    // zero the "face-geometry" elements...
    FOR_I3 x_fa[ifa][i] = 0.0;
    FOR_I3 n_fa[ifa][i] = 0.0;
    // loop to get total face normal and approximate x_fa...
    double x_fa_tmp[3] = { 0.0, 0.0, 0.0 };
    const int ino_f = noofa_v[noofa_i[ifa]];
    double wgt_sum = 0.0;
    int ino1 = noofa_v[noofa_i[ifa+1]-1];
    for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      const int ino0 = ino1;
      ino1 = noofa_v[nof];
      const double wgt = DIST(x_no[ino0],x_no[ino1]);
      FOR_I3 x_fa_tmp[i] += wgt*(x_no[ino0][i]+x_no[ino1][i]);
      wgt_sum += wgt;
      const double n[3] = TRI_NORMAL_2(x_no[ino_f],x_no[ino0],x_no[ino1]);
      FOR_I3 n_fa[ifa][i] += 0.5*n[i];
    }
    assert(wgt_sum > 0.0);
    FOR_I3 x_fa_tmp[i] /= wgt_sum*2.0;
    // get exact x_fa...
    wgt_sum = 0.0;
    ino1 = noofa_v[noofa_i[ifa+1]-1];
    for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      const int ino0 = ino1;
      ino1 = noofa_v[nof];
      const double n[3] = TRI_NORMAL_2(x_fa_tmp,x_no[ino0],x_no[ino1]);
      const double wgt = DOT_PRODUCT(n_fa[ifa],n); // area squared weighting !?
      FOR_I3 x_fa[ifa][i] += wgt*(x_fa_tmp[i]+x_no[ino0][i]+x_no[ino1][i]);
      wgt_sum += wgt;
    }
    if (wgt_sum > 0.0) {
      FOR_I3 x_fa[ifa][i] /= wgt_sum*3.0;
    }
    else {
      // this must be a zero-area face -- so just pick a node...
      FOR_I3 x_fa[ifa][i] = x_no[ino1][i];
    }
  }

  // area over delta based on forming points

  FOR_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    const double area = MAG(n_fa[ifa]);
    const double delta = DIST(x_vv[icv0],x_vv[icv1]);
    area_over_delta_fa[ifa] = area/delta;
  }

  // build volume at boundary cells assuming planar boundary faces (true on interior)...

  FOR_ICV planar_vol_ratio_cv[icv] = 0.0;
  FOR_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    const double dx[3] = DIFF(x_vv[icv1],x_vv[icv0]);
    const double vol = DOT_PRODUCT(n_fa[ifa],dx)/6.0;
    planar_vol_ratio_cv[icv0] += vol;
    if (icv1 < ncv) 
      planar_vol_ratio_cv[icv1] += vol;
  }
  FOR_IBF {
    const int icv = cvobf[ibf];
    assert((icv >= 0)&&(icv < ncv));
    const double dx[3] = DIFF(x_bf[ibf],x_vv[icv]);
    const double vol = DOT_PRODUCT(n_bf[ibf],dx)/3.0;
    planar_vol_ratio_cv[icv] += vol;
  }
  updateCvData(planar_vol_ratio_cv);
  FOR_ICV_G planar_vol_ratio_cv[icv] /= vol_cv[icv];

}

void HelmholtzSolver::initFromCoarseGrid(CoarseGrid* cg) {
  StaticSolver::initFromCoarseGrid(cg);

  // note that we need to allocate with room for restriction/prolongation
  assert(rho == NULL); rho = new double[cg->ncc_g];
  assert(p == NULL); p = new double[cg->ncc_g];
  assert(u == NULL); u = new double[cg->ncc_g][3];
  assert(mu_lam == NULL); mu_lam = new double[cg->ncc_g];
  assert(mu_sgs == NULL); mu_sgs = new double[cg->ncc_g];
  assert(a_sgs == NULL); a_sgs = new double[cg->ncc_g];
  assert(mf == NULL); mf = new double[cg->ncf];
  assert(A_diag_source == NULL); A_diag_source = new double[cg->ncc][3];

  FOR_IZONE(bfZoneVec) {
    const string zone_name = bfZoneVec[izone].getName();
    if ( Param* p = getParam(zone_name)) {
      const string bc_type = p->getString(0);
      if ( (bc_type == "SLIP")||(bc_type == "SYMMETRY")) {
        bcs.push_back(new SlipWallHBc(&bfZoneVec[izone],this,cg_level));
      } else if ( bc_type == "WALL") {
        bcs.push_back(new WallHBc(&bfZoneVec[izone],this,cg_level));
      } else if ( bc_type == "WM_ALG_WALL" ) {
        bcs.push_back(new AlgebraicWallModelHBc(&bfZoneVec[izone],this,cg_level));
      } else if ( bc_type == "WM_SLIP_WALL") {
        bcs.push_back(new SlipWallModelHBc(&bfZoneVec[izone],this,cg_level));
      } else if ( bc_type == "INLET") {
        bcs.push_back(new InletHBc(&bfZoneVec[izone],this,cg_level));
      } else if ( bc_type == "INLET_PROFILE") {
        bcs.push_back(new InletHBcProfile(&bfZoneVec[izone],this,cg_level));
      } else if ( bc_type == "INFLOW_TURB") {
        bcs.push_back(new InflowTurbulenceHBc(&bfZoneVec[izone],this,cg_level));
      } else if ( bc_type == "OUTLET") {
        bcs.push_back(new OutletHBc(&bfZoneVec[izone],this,cg_level));
      } else if ( bc_type == "OUTLET_VV") {
        bcs.push_back(new OutletVVHBc(&bfZoneVec[izone],this,cg_level));
      } else if ( bc_type == "SPONGE") {
        bcs.push_back(new SpongeHBc(&bfZoneVec[izone],this,cg_level));
      } else if ( bc_type == "HOOK") {
        bcs.push_back(initHookBcMg(&bfZoneVec[izone],cg_level));
      } else {
        assert(0);
      }
    } else {
      assert(0);
    }
  }

  // init with room for restriction/prolongation
  FOR_BCZONE (*it)->initFromCoarseGrid(cg);

}

void HelmholtzSolver::updateTimeLevels(){

  // reset the previous time states .. the native algorithm
  // does not necessarily need the previous time states in the 
  // ghost regions, but they are kept in case someone requests 
  // a gradient of the previous time levels for instance.

  FOR_ICV_G {
    
    FOR_I3 {
      u00[icv][i] = u0[icv][i];
      u0[icv][i]  = u[icv][i];
    }

    rho00[icv] = rho0[icv];
    rho0[icv]  = rho[icv];

  }


  if (!b_helmholtzCorrection){ // one-helmholtz algorithm uses extrapolated mf by default but may be switched off
    if (mf_init>0){
      COUT1("MF_INIT: extrapolation off for " << mf_init << " more step(s) while time history is populated");
      --mf_init;
      FOR_IFA {
        mf00[ifa]   = mf0[ifa];
        mf0[ifa]    = mf[ifa];
      }
    }
    else if (checkParam("NO_EXTRAPOLATE_MF")||n_outer_iters!=1){
      if (step%check_interval==0){
        COUT1("NO_EXTRAPOLATE_MF: extrapolation of face-based mass flow is off");
      }
      FOR_IFA {
        mf00[ifa]   = mf0[ifa];
        mf0[ifa]    = mf[ifa];
      }
    }
    else{
      // three step ab extrapolator for face-based mf...
      FOR_IFA {
        const double tmp = 3.0*mf[ifa] - 3.0*mf0[ifa] + mf00[ifa];
        mf00[ifa]   = mf0[ifa];
        mf0[ifa]    = mf[ifa];
        mf[ifa]     = tmp; 
      }
    }
  }
}

void HelmholtzSolver::updateSgs() {

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

  // acoustic modeling...
  calcASgs();

}

void HelmholtzSolver::correctContinuity(){

  if (b_helmholtzCorrection)
    correctContinuityHelmholtz();
  else 
    correctContinuityAlgebraic();
}

void HelmholtzSolver::correctContinuityHelmholtz() {

  if ( (mpi_rank==0) && (step%check_interval==0) )
    cout << "correctContinuityHelmholtz()" << endl;

  // start by correcting the rho to include the current pressure...
  // this is consistent with the extrapolated mf's that include the
  // pressure correction effects from previous time steps...

  const double inv_sos2 = 1.0/(helmholtz_sos*helmholtz_sos);

  //DAPDAP TODO - is this needed?  Is this currently the value of rho?
  FOR_ICV rho[icv] = rho_uncorrected[icv] + p[icv]*inv_sos2; 

  double * rhs = new double[ncv];
  calcContinuity(rhs,true);

  // scale rhs...
  FOR_ICV rhs[icv] *= 1.5/dt;

  // since we are solving a correction eqn then, the orig
  // field of rho should satisfy the bcs.  this is true in
  // a linearized sense for the case with no outer iters

  // now solve for phi...

  double * phi = new double[ncv_g];
  FOR_ICV_G phi[icv] = 0.0; // this is a correction equation, so initialize with zero

  if (ncg_p > 0)
    solveHelmholtzMg(phi,rhs,false);
  else {
    solveHelmholtz(phi,rhs,false);
    //solveHelmholtzPatr(phi,rhs,false); // HACK for testing MG
  }

  // correct density...

  FOR_ICV rho[icv] += phi[icv]*inv_sos2;
  updateCvData(rho); 

  // correct mf at internal faces...
  for (int ifa = 0; ifa < nfa; ++ifa) {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    mf[ifa] -= 2.0/3.0*dt*area_over_delta_fa[ifa]*(phi[icv1]-phi[icv0]);
  }
  delete[] phi;

  calcContinuity(rhs,true);

  if (step%check_interval==0)
    dumpRange(rhs,ncv,"continuity");

  delete[] rhs;
}

void HelmholtzSolver::correctContinuityAlgebraic() {

  if ( (mpi_rank==0) && (step%check_interval==0) )
    cout << "correctContinuityAlgebraic()" << endl;

  // temporarily store sum_fa mf_fa in rho...

  FOR_ICV rho[icv] = 0.0;
  FOR_BCZONE (*it)->addMassFlux(rho);
  for (int ifa = 0; ifa < nfa ; ++ifa){
    const int icv0 = cvofa[ifa][0];
    rho[icv0] += mf[ifa];
    const int icv1 = cvofa[ifa][1];
    if (icv1 < ncv)
      rho[icv1] -= mf[ifa];
  }
  const double inv_sos2 = 1.0/(helmholtz_sos*helmholtz_sos);
  for (int ifa = 0; ifa < nfa ; ++ifa){
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    const double flux = 0.5*(a_sgs[icv1]+a_sgs[icv0])*(p[icv1]-p[icv0])*inv_sos2*area_over_delta_fa[ifa];
    rho[icv0] -= flux;
    if (icv1 < ncv)
      rho[icv1] += flux;
  }

  // set rho to satisfy continuity discretely...

  // must take into account any sponge zones...
  double * sponge_sigma = new double[ncv];
  FOR_ICV sponge_sigma[icv] = 0.0;
  FOR_BCZONE {
    Param * param = getParam((*it)->getName());
    if ( param->getString(0) == "SPONGE") {
      SpongeHBc * sponge = dynamic_cast<SpongeHBc*>(*it);
      sponge->computeSigma(sponge_sigma);
    }
  }
  //FOR_ICV rho[icv] = (4.0*rho0[icv]-rho00[icv]-2.0*dt*rho[icv]/vol_cv[icv])/3.0;
  FOR_ICV rho[icv] = (2.0*rho0[icv]-0.5*rho00[icv]-dt*rho[icv]/vol_cv[icv]+dt*sponge_sigma[icv]*rho_uncorrected[icv]) /
                     (1.5 + dt*sponge_sigma[icv]);


  double * rhs = new double[ncv];
  calcContinuity(rhs,true);
  FOR_BCZONE {
    Param * param = getParam((*it)->getName());
    if ( param->getString(0) == "SPONGE") {
      SpongeHBc * sponge = dynamic_cast<SpongeHBc*>(*it);
      sponge->addContSource(rhs);
    }
  }

  if (step%check_interval==0){
    //dumpRange(sponge_sigma,ncv,"sponge_sigma");
    dumpRange(rhs,ncv,"continuity");
  }

  delete[] sponge_sigma;
  delete[] rhs;
}

void HelmholtzSolver::buildHelmholtzMg(double * &A,double* rhs,const bool b_pressure) {
  assert(A == NULL);

  const int cvocv_s = cvocv_i[ncv];
  A = new double[cvocv_s];

  const double sos2 = helmholtz_sos*helmholtz_sos;
  FOR_ICV {
    const int coc_f = cvocv_i[icv];
    A[coc_f] = -9.0/4.0*vol_cv[icv]/(sos2*dt*dt);
    for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc)
      A[coc] = 0.0;
  }

  if (b_pressure) {

    for (int ifa = 0; ifa < nfa; ++ifa) {
      const int icv0 = cvofa[ifa][0];
      const int coc00 = cvocv_i[icv0];
      const int icv1 = cvofa[ifa][1];
      const double coeff = (1.0 + 1.5*0.5*(a_sgs[icv1]+a_sgs[icv0])/(sos2*dt))*area_over_delta_fa[ifa];
      int coc01 = coc00+1;
      while (cvocv_v[coc01] != icv1) {
        ++coc01;
        assert(coc01 != cvocv_i[icv0+1]);
      }
      A[coc00] -= coeff;
      A[coc01] += coeff;

      if (icv1 < ncv) {
        const int coc11 = cvocv_i[icv1];
        int coc10 = coc11+1;
        while (cvocv_v[coc10] != icv0) {
          ++coc10;
          assert(coc10 != cvocv_i[icv1+1]);
        }
        A[coc11] -= coeff;
        A[coc10] += coeff;
      }
    }

    FOR_BCZONE (*it)->addPressureFlux(A,rhs);

  }
  else {

    for (int ifa = 0; ifa < nfa; ++ifa) {
      const int icv0 = cvofa[ifa][0];
      const int coc00 = cvocv_i[icv0];
      const int icv1 = cvofa[ifa][1];
      int coc01 = coc00+1;
      while (cvocv_v[coc01] != icv1) {
        ++coc01;
        assert(coc01 != cvocv_i[icv0+1]);
      }
      A[coc00] -= area_over_delta_fa[ifa];
      A[coc01] += area_over_delta_fa[ifa];

      if (icv1 < ncv) {
        const int coc11 = cvocv_i[icv1];
        int coc10 = coc11+1;
        while (cvocv_v[coc10] != icv0) {
          ++coc10;
          assert(coc10 != cvocv_i[icv1+1]);
        }
        A[coc11] -= area_over_delta_fa[ifa];
        A[coc10] += area_over_delta_fa[ifa];
      }
    }

  }

}

void HelmholtzSolver::solveHelmholtzPatr(double * phi,double * rhs,const bool b_pressure) {

  if ( (mpi_rank==0) && (step%check_interval==0) )
    cout << "solveHelmholtzPatr()" << endl;

  double *A = NULL;
  buildHelmholtzMg(A,rhs,b_pressure);

  const double zero = getDoubleParam("PRESSURE_ZERO", 1.0E-6);
  const int maxiter  = getIntParam("PRESSURE_MAXITER",1000);
  solveCvPatr(phi,A,A,rhs,zero,1.0,1.0,maxiter,true);

  delete[] A;

}

void HelmholtzSolver::solveHelmholtzMg(double * phi,double * rhs,const bool b_pressure) {

  if ( (mpi_rank==0) && (step%check_interval==0) )
    cout << "solveHelmholtzMg()" << endl;

  // build lhs matrices (and modify rhs's

  double *A = NULL;
  buildHelmholtzMg(A,rhs,b_pressure); 
  double *inv_diag = new double[ncv_g];
  FOR_ICV inv_diag[icv] = 1.0/A[cvocv_i[icv]];
  updateCvData(inv_diag);
  //dumpRange(inv_diag,ncv,"inv_diag 0");
  double *res = new double[ncv];
  //dumpRange(rhs,ncv,"rhs");

  // we need the following work arrays...

  double *r = new double[ncv_g];
  double *p = new double[ncv_g];
  double *Ar = new double[ncv];
  double *Ap = new double[ncv];

  vector<double *> A_cg(ncg_p);
  vector<double *> inv_diag_cg(ncg_p);
  vector<double *> err_cg(ncg_p);
  vector<double *> res_cg(ncg_p);
  for (int icg = 0; icg < ncg_p; ++icg) {

    // get a_sgs on coarse grid
    // only needed if momentum multigrid not enabled
    // on the coarse grid level... 
    if (b_pressure && (icg > (ncg_u-1))) {
      if (icg == 0) {
        cgs[icg]->restrictCcData(css[icg]->a_sgs,a_sgs);
      }
      else {
        cgs[icg]->restrictCcData(css[icg]->a_sgs,css[icg-1]->a_sgs);
      }
      // ghosts
      css[icg]->updateCvData(css[icg]->a_sgs);
    }

    A_cg[icg] = NULL;
    css[icg]->buildHelmholtzMg(A_cg[icg],NULL,b_pressure); 
    inv_diag_cg[icg] = new double[css[icg]->ncv_g];
    for (int icv = 0; icv < css[icg]->ncv; ++icv)
      inv_diag_cg[icg][icv] = 1.0/A_cg[icg][css[icg]->cvocv_i[icv]];
    css[icg]->updateCvData(inv_diag_cg[icg]);
    //dumpRange(inv_diag_cg[icg],css[icg]->ncv,"inv_diag");
    err_cg[icg] = new double[cgs[icg]->ncc_g]; // MUST be ncc_g for prolongation
    res_cg[icg] = new double[cgs[icg]->ncc_g]; // MUST be ncc_g for restriction
  }

  // solve using V-cycle 

  const double zero = getDoubleParam("PRESSURE_ZERO",1.0E-6);
  const double maxiter = getDoubleParam("PRESSURE_MAXITER",100);
  const double relax = getDoubleParam("PRESSURE_RELAX",1.0); // jacobi needs something around this

  int iter = 0;
  int done = 0;
  while (done == 0) {
    iter++;

    if (pressureSmootherVec[0] == 0)
      smoothCvJacobi(phi,r,inv_diag,A,rhs,pressureNsmoothVec[0],relax);
    else if (pressureSmootherVec[0] == 1)
      smoothCvSgs(phi,inv_diag,A,rhs,pressureNsmoothVec[0],relax);
    else if (pressureSmootherVec[0] == 2)
      smoothCvPatr(phi,r,p,Ar,Ap,inv_diag,A,A,rhs,pressureNsmoothVec[0],relax,relax);
    else
      smoothCvCg(phi,Ar,Ap,p,inv_diag,A,rhs,pressureNsmoothVec[0]);

    // compute residual 

    calcCvResidual(res,phi,A,rhs);

    // check if done

    double my_res_max = 0.0;
    FOR_ICV my_res_max = max(my_res_max,fabs(res[icv]*inv_diag[icv]));
    double res_max;
    MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
    if (mpi_rank == 0) {
      if (step%check_interval==0) 
        cout << " > solveHelmholtz iter, res_max: " << iter << " " << res_max << endl;
      if (res_max < zero) {
        done = 1;
      }
      else if (iter > maxiter) {
        cout << " > Warning: solveHelmholtz did not converge after " << maxiter << " iters, res_max: " << res_max << endl;
        done = 2;
      }
    }
    MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

    if (done != 0)
      break;

    // restrict residual 

    //cgs[0]->restrictCcData(res_cg[0],res); 
    cgs[0]->restrictExtrinsicCcData(res_cg[0],res); 
    css[0]->updateCvData(res_cg[0]);

    for (int icg = 0; icg < ncg_p-1; ++icg) {

      //if (mpi_rank == 0)
      //  cout << "icg: " << icg << endl;

      // smooth

      for (int icc = 0; icc < css[icg]->ncv_g; ++icc) 
        err_cg[icg][icc] = 0.0;

      if (pressureSmootherVec[icg+1] == 0)
        css[icg]->smoothCvJacobi(err_cg[icg],r,inv_diag_cg[icg],A_cg[icg],res_cg[icg],pressureNsmoothVec[icg+1],relax);
      else if (pressureSmootherVec[icg+1] == 1)
        css[icg]->smoothCvSgs(err_cg[icg],inv_diag_cg[icg],A_cg[icg],res_cg[icg],pressureNsmoothVec[icg+1],relax);
      else if (pressureSmootherVec[icg+1] == 2)
        css[icg]->smoothCvPatr(err_cg[icg],r,p,Ar,Ap,inv_diag_cg[icg],A_cg[icg],A_cg[icg],res_cg[icg],pressureNsmoothVec[icg+1],relax,relax);
      else
        css[icg]->smoothCvCg(err_cg[icg],Ar,Ap,p,inv_diag_cg[icg],A_cg[icg],res_cg[icg],pressureNsmoothVec[icg+1]);

      // compute residual 

      css[icg]->calcCvResidual(r,err_cg[icg],A_cg[icg],res_cg[icg]);

      // restrict residual 

      //cgs[icg+1]->restrictCcData(res_cg[icg+1],r); 
      cgs[icg+1]->restrictExtrinsicCcData(res_cg[icg+1],r); 
      css[icg+1]->updateCvData(res_cg[icg+1]);

    }

    //if (mpi_rank == 0)
    //  cout << "icg: " << ncg_p-1 << endl;

    // smooth on coarsest grid

    for (int icc = 0; icc < css[ncg_p-1]->ncv_g; ++icc) 
      err_cg[ncg_p-1][icc] = 0.0;

    if (b_solve_p_on_coarsest) {
      css[ncg_p-1]->solveCvCg(err_cg[ncg_p-1],A_cg[ncg_p-1],res_cg[ncg_p-1],zero,pressureNsmoothVec[ncg_p],false);
    }
    else {
      if (pressureSmootherVec[ncg_p] == 0)
        css[ncg_p-1]->smoothCvJacobi(err_cg[ncg_p-1],r,inv_diag_cg[ncg_p-1],A_cg[ncg_p-1],res_cg[ncg_p-1],pressureNsmoothVec[ncg_p],relax);
      else if (pressureSmootherVec[ncg_p] == 1)
        css[ncg_p-1]->smoothCvSgs(err_cg[ncg_p-1],inv_diag_cg[ncg_p-1],A_cg[ncg_p-1],res_cg[ncg_p-1],pressureNsmoothVec[ncg_p],relax);
      else if (pressureSmootherVec[ncg_p] == 2)
        css[ncg_p-1]->smoothCvPatr(err_cg[ncg_p-1],r,p,Ar,Ap,inv_diag_cg[ncg_p-1],A_cg[ncg_p-1],A_cg[ncg_p-1],res_cg[ncg_p-1],pressureNsmoothVec[ncg_p],relax,relax);
      else
        css[ncg_p-1]->smoothCvCg(err_cg[ncg_p-1],Ar,Ap,p,inv_diag_cg[ncg_p-1],A_cg[ncg_p-1],res_cg[ncg_p-1],pressureNsmoothVec[ncg_p]);
    }

    for (int icg = ncg_p-1; icg >= 1; --icg) {

      // prolong and update guess

      cgs[icg]->updateCcIDataReverse(err_cg[icg]); 
      cgs[icg]->prolongCcDataAndUpdateGuess(err_cg[icg-1],err_cg[icg]);
      css[icg-1]->updateCvData(err_cg[icg-1]);

      // smooth

      if (pressureSmootherVec[icg] == 0)
        css[icg-1]->smoothCvJacobi(err_cg[icg-1],r,inv_diag_cg[icg-1],A_cg[icg-1],res_cg[icg-1],pressureNsmoothVec[icg],relax);
      else if (pressureSmootherVec[icg] == 1)
        css[icg-1]->smoothCvSgs(err_cg[icg-1],inv_diag_cg[icg-1],A_cg[icg-1],res_cg[icg-1],pressureNsmoothVec[icg],relax);
      else if (pressureSmootherVec[icg] == 2)
        css[icg-1]->smoothCvPatr(err_cg[icg-1],r,p,Ar,Ap,inv_diag_cg[icg-1],A_cg[icg-1],A_cg[icg-1],res_cg[icg-1],pressureNsmoothVec[icg],relax,relax);
      else
        css[icg-1]->smoothCvCg(err_cg[icg-1],Ar,Ap,p,inv_diag_cg[icg-1],A_cg[icg-1],res_cg[icg-1],pressureNsmoothVec[icg]);

    }

    // prolong and update guess

    cgs[0]->updateCcIDataReverse(err_cg[0]); 
    cgs[0]->prolongCcDataAndUpdateGuess(phi,err_cg[0]);
    updateCvData(phi);

  }

  // cleanup

  delete[] A;
  delete[] inv_diag;
  delete[] res;
  delete[] r;
  delete[] p;
  delete[] Ar;
  delete[] Ap;

  for (int icg = 0; icg < ncg_p; ++icg) {
    delete[] A_cg[icg];
    delete[] inv_diag_cg[icg];
    delete[] err_cg[icg];
    delete[] res_cg[icg];
  }

}

void HelmholtzSolver::solveHelmholtz(double * p,double * rhs,const bool pressure) {

  const int cvocv_s = cvocv_i[ncv];
  double * A = new double[cvocv_s];

  const double sos2 = helmholtz_sos*helmholtz_sos;
  FOR_ICV {
    const int coc_f = cvocv_i[icv];
    A[coc_f] = -9.0/4.0*vol_cv[icv]/(sos2*dt*dt);
    for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc)
      A[coc] = 0.0;
  }


  if (pressure) {

    for (int ifa = 0; ifa < nfa; ++ifa) {
      const int icv0 = cvofa[ifa][0];
      const int coc00 = cvocv_i[icv0];
      const int icv1 = cvofa[ifa][1];
      const double coeff = (1.0 + 1.5*0.5*(a_sgs[icv1]+a_sgs[icv0])/(sos2*dt))*area_over_delta_fa[ifa];
      int coc01 = coc00+1;
      while (cvocv_v[coc01] != icv1) {
        ++coc01;
        assert(coc01 != cvocv_i[icv0+1]);
      }
      A[coc00] -= coeff;
      A[coc01] += coeff;

      if (icv1 < ncv) {
        const int coc11 = cvocv_i[icv1];
        int coc10 = coc11+1;
        while (cvocv_v[coc10] != icv0) {
          ++coc10;
          assert(coc10 != cvocv_i[icv1+1]);
        }
        A[coc11] -= coeff;
        A[coc10] += coeff;
      }

    }

    FOR_BCZONE (*it)->addPressureFlux(A,rhs);

  }
  else {

    for (int ifa = 0; ifa < nfa; ++ifa) {
      const int icv0 = cvofa[ifa][0];
      const int coc00 = cvocv_i[icv0];
      const int icv1 = cvofa[ifa][1];
      int coc01 = coc00+1;
      while (cvocv_v[coc01] != icv1) {
        ++coc01;
        assert(coc01 != cvocv_i[icv0+1]);
      }
      A[coc00] -= area_over_delta_fa[ifa];
      A[coc01] += area_over_delta_fa[ifa];

      if (icv1 < ncv) {
        const int coc11 = cvocv_i[icv1];
        int coc10 = coc11+1;
        while (cvocv_v[coc10] != icv0) {
          ++coc10;
          assert(coc10 != cvocv_i[icv1+1]);
        }
        A[coc11] -= area_over_delta_fa[ifa];
        A[coc10] += area_over_delta_fa[ifa];
      }

    }

  }

  // and solve...
  const double pressure_zero = getDoubleParam("PRESSURE_ZERO",1.0E-6);
  const int pressure_maxiter = getIntParam("PRESSURE_MAXITER",1000);
  solveCvCg(p,A,rhs,pressure_zero,pressure_maxiter,false); 

  delete[] A;

}

void HelmholtzSolver::calcContinuity(double * cont,const bool b_include_a_sgs) {

  FOR_ICV cont[icv] = vol_cv[icv]*(1.5*rho[icv] - 2.0*rho0[icv] + 0.5*rho00[icv])/dt;

  FOR_BCZONE (*it)->addMassFlux(cont);

  for (int ifa = 0; ifa < nfa ; ++ifa){
    const int icv0 = cvofa[ifa][0];
    cont[icv0] += mf[ifa];
    const int icv1 = cvofa[ifa][1];
    if (icv1 < ncv)
      cont[icv1] -= mf[ifa];
  }

  if (b_include_a_sgs) {
    const double inv_sos2 = 1.0/(helmholtz_sos*helmholtz_sos);
    for (int ifa = 0; ifa < nfa ; ++ifa){
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      const double flux = 0.5*(a_sgs[icv1]+a_sgs[icv0])*(p[icv1]-p[icv0])*inv_sos2*area_over_delta_fa[ifa];
      cont[icv0] -= flux;
      if (icv1 < ncv)
        cont[icv1] += flux;
    }
  }

}

void HelmholtzSolver::solveU() {

  if ( (mpi_rank==0) && (step%check_interval==0) )
    cout << "solveU()" << endl;

  // solve...
  // [A]*{u} = {rhs}

  const int cvocv_s = cvocv_i[ncv];
  double * A = new double[cvocv_s];
  for (int coc = 0; coc < cvocv_s; ++coc)
    A[coc] = 0.0;

  double (*rhs)[3] = new double[ncv][3];
  FOR_ICV FOR_I3 rhs[icv][i] = -dpdx[icv][i]*vol_cv[icv];

  addFrameRotationSourceTerms(rhs);

  // bcs...
  FOR_BCZONE (*it)->addMomentumFlux(A,rhs);

  const double inv_sos2 = 1.0/(helmholtz_sos*helmholtz_sos);

  // internal faces...
  for (int ifa = 0; ifa < nfa; ++ifa){

    const int icv0 = cvofa[ifa][0];
    assert((icv0 >= 0)&&(icv0 < ncv));
    const int icv1 = cvofa[ifa][1];
    assert((icv1 >= 0)&&(icv0 < ncv_g));

    const int coc00 = cvocv_i[icv0];
    int coc01 = coc00+1;
    while (cvocv_v[coc01] != icv1) {
      ++coc01;
      assert(coc01 != cvocv_i[icv0+1] );
    }

    const double mu_coeff = 0.5*(mu_lam[icv0] + mu_lam[icv1] + mu_sgs[icv0] + mu_sgs[icv1])*area_over_delta_fa[ifa]; // + fabs(cj_rhouj);

    const double asgs_coeff = 0.5*(a_sgs[icv0] + a_sgs[icv1])*inv_sos2*(p[icv1]-p[icv0])*area_over_delta_fa[ifa];

    A[coc00] += 0.5*(mf[ifa] - asgs_coeff) + mu_coeff;
    A[coc01] += 0.5*(mf[ifa] - asgs_coeff) - mu_coeff;

    if (icv1 < ncv) {

      const int coc11 = cvocv_i[icv1];
      int coc10 = coc11+1;
      while (cvocv_v[coc10] != icv0) {
        ++coc10;
        assert(coc10 != cvocv_i[icv1+1] );
      }

      A[coc11] -= 0.5*(mf[ifa] - asgs_coeff) - mu_coeff;
      A[coc10] -= 0.5*(mf[ifa] - asgs_coeff) + mu_coeff;

    }

  }

  // time term...

  FOR_ICV {
    A[cvocv_i[icv]] += 1.5*vol_cv[icv]*rho[icv]/dt;
    FOR_I3 rhs[icv][i] += vol_cv[icv]*( 2.0*rho0[icv]*u0[icv][i] - 0.5*rho00[icv]*u00[icv][i] )/dt;
  }

  // source terms...

  momentumSourceHook(A,rhs);
  if (rhou_source) {
    FOR_ICV {
      FOR_I3 {
        rhs[icv][i] += rhou_source[icv][i];
      }
    }
  }

  // and solve...

  const double mom_zero = getDoubleParam("MOMENTUM_ZERO", 1.0E-6);
  const double mom_relax = getDoubleParam("MOMENTUM_RELAX",0.7);
  const int mom_maxiter  = getIntParam("MOMENTUM_MAXITER", 100);
  if (A_diag_source) {
    //COUT1("solveCvJacobi A_diag");
    solveCvJacobi(u,A,A_diag_source,rhs,mom_zero,mom_relax,mom_maxiter,step%check_interval==0); 
  }
  else{
    solveCvJacobi(u,A,rhs,mom_zero,mom_relax,mom_maxiter,step%check_interval==0); 
  }

  // cleanup...

  delete[] A;
  delete[] rhs;

}

void HelmholtzSolver::solveUPatr() {

  if ( (mpi_rank==0) && (step%check_interval==0) )
    cout << "solveUPatr()" << endl;

  double *A = NULL;
  double *At = NULL;
  double (*rhs)[3] = new double[ncv][3];
  buildUMg(A,At,rhs);

  const double zero = getDoubleParam("MOMENTUM_ZERO", 1.0E-6);
  const int maxiter  = getIntParam("MOMENTUM_MAXITER",1000);
  const double relax = getDoubleParam("MOMENTUM_RELAX",1.0);

  if (A_diag_source) {
    //COUT1("solveCvPatr A_diag");
    solveCvPatr(u,A,At,A_diag_source,rhs,zero,relax,relax,maxiter,false);
  }
  else{
    solveCvPatr(u,A,At,rhs,zero,relax,relax,maxiter,false);
  }


  delete[] A;
  delete[] At;
  delete[] rhs;

}

void HelmholtzSolver::solveUTim() {

  if ( (mpi_rank==0) && (step%check_interval==0) )
    cout << "solveUTim()" << endl;

  if (A_diag_source){
    CERR("MOMENTUM_SOLVER TIM (triangular iterative method) not supported when using implicit body forces");
  }

  double *A = NULL;
  double *At = NULL;
  double (*rhs)[3] = new double[ncv][3];
  buildUMg(A,At,rhs);
  double *inv_diag = new double[ncv_g];
  FOR_ICV inv_diag[icv] = 1.0/A[cvocv_i[icv]];
  updateCvData(inv_diag);

  FOR_ICV {
    FOR_I3 rhs[icv][i] *= inv_diag[icv];
    for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) {
      A[coc] *= inv_diag[icv];
      At[coc] *= inv_diag[cvocv_v[coc]];
    }
  }

  double *A0 = new double[cvocv_i[ncv]];
  double *A1 = new double[cvocv_i[ncv]];
  for (int coc = 0; coc < cvocv_i[ncv]; ++coc) {
    A0[coc] = 0.5*(A[coc]+At[coc]);
    A1[coc] = 0.5*(A[coc]-At[coc]);
  }

  //double g1 = HUGE_VAL;
  double g2 = 0.0;
  double g3 = 0.0;
  FOR_ICV {
    double this_A1_inf = 0.0;
    double this_A0_inf = 0.0;
    for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) {
      this_A1_inf += fabs(A1[coc]);
      this_A0_inf += fabs(A0[coc]);
    }
    g3 = max(g3,this_A1_inf);
    g2 = max(g2,this_A0_inf);
    //g1 = min(g1,this_A0_inf);
  }
  const double tau_bound = 2.0/(2.0*g3+g2);
  //tau_opt = 4.0/(g1+g2+sqrt(pow(g2-g1+4.0*g3,2)+4.0*g1*g2));

  if (mpi_rank == 0) 
    cout << " > tau bound: " <<  tau_bound << endl;

  const double zero = getDoubleParam("MOMENTUM_ZERO", 1.0E-6);
  const int maxiter  = getIntParam("MOMENTUM_MAXITER",1000);
  const double relax = getDoubleParam("MOMENTUM_RELAX",1.0);
  solveCvTim(u,A,A1,rhs,zero,tau_bound,relax,maxiter,false);

  delete[] A;
  delete[] At;
  delete[] A0;
  delete[] A1;
  delete[] rhs;
  delete[] inv_diag;

}

void HelmholtzSolver::solveUMg() {
  assert(ncg_u >= 1);

  if ( (mpi_rank==0) && (step%check_interval==0) )
    cout << "solveUMg()" << endl;

  // build lhs matrices (and modify rhs's)

  double *A = NULL;
  double *At = NULL;
  double (*rhs)[3] = new double[ncv][3];
  buildUMg(A,At,rhs);
  double *inv_diag = new double[ncv_g];
  FOR_ICV inv_diag[icv] = 1.0/A[cvocv_i[icv]];
  updateCvData(inv_diag);
  //dumpRange(A,cvocv_i[ncv],"A");
  double (*res)[3] = new double[ncv][3];

  bool b_skew_decomp = false;
  for (int icg = 0; icg <= ncg_u; ++icg) {
    if (momentumSmootherVec[icg] == 3)
      b_skew_decomp = true;
  }

  double *A0 = NULL;
  double *A1 = NULL;
  double tau_bound;
  if (b_skew_decomp) {

    A0 = new double[cvocv_i[ncv]];
    A1 = new double[cvocv_i[ncv]];
    FOR_ICV {
      for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) {
        A0[coc] = 0.5*(inv_diag[icv]*A[coc]+inv_diag[cvocv_v[coc]]*At[coc]);
        A1[coc] = 0.5*(inv_diag[icv]*A[coc]-inv_diag[cvocv_v[coc]]*At[coc]);
      }
    }

    //double g1 = HUGE_VAL;
    double g2 = 0.0;
    double g3 = 0.0;
    FOR_ICV {
      double this_A1_inf = 0.0;
      double this_A0_inf = 0.0;
      for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) {
        this_A1_inf += fabs(A1[coc]);
        this_A0_inf += fabs(A0[coc]);
      }
      g3 = max(g3,this_A1_inf);
      g2 = max(g2,this_A0_inf);
      //g1 = min(g1,this_A0_inf);
    }
    tau_bound = 2.0/(2.0*g3+g2);
    //tau_opt = 4.0/(g1+g2+sqrt(pow(g2-g1+4.0*g3,2)+4.0*g1*g2));

    if (mpi_rank == 0) 
      cout << " > tau bound: " <<  tau_bound << endl;
  }

  // we need the following work arrays...

  double (*r)[3] = new double[ncv_g][3];
  double (*ps)[3] = new double[ncv_g][3];
  double (*Ar)[3] = new double[ncv][3];
  double (*Aps)[3] = new double[ncv][3];

  vector<double *> A_cg(ncg_u);
  vector<double *> At_cg(ncg_u);
  vector<double *> inv_diag_cg(ncg_u);
  vector<double (*)[3]> err_cg(ncg_u);
  vector<double (*)[3]> res_cg(ncg_u);
  vector<double *> A0_cg(ncg_u);
  vector<double *> A1_cg(ncg_u);
  vector<double> tau_bound_cg(ncg_u);
  for (int icg = 0; icg < ncg_u; ++icg) {

    // get rho, mu_lam, mu_sgs and mf on coarse grid
    if (icg == 0) {
      cgs[icg]->restrictCcData(css[icg]->rho,rho);
      cgs[icg]->restrictCcData(css[icg]->p,p);
      cgs[icg]->restrictCcData(css[icg]->u,u);
      cgs[icg]->restrictCcData(css[icg]->mu_lam,mu_lam);
      cgs[icg]->restrictCcData(css[icg]->mu_sgs,mu_sgs);
      cgs[icg]->restrictCcData(css[icg]->a_sgs,a_sgs);
      cgs[icg]->restrictExtrinsicSignedCfData(css[icg]->mf,mf);
    }
    else {
      cgs[icg]->restrictCcData(css[icg]->rho,css[icg-1]->rho);
      cgs[icg]->restrictCcData(css[icg]->p,css[icg-1]->p);
      cgs[icg]->restrictCcData(css[icg]->u,css[icg-1]->u);
      cgs[icg]->restrictCcData(css[icg]->mu_lam,css[icg-1]->mu_lam);
      cgs[icg]->restrictCcData(css[icg]->mu_sgs,css[icg-1]->mu_sgs);
      cgs[icg]->restrictCcData(css[icg]->a_sgs,css[icg-1]->a_sgs);
      cgs[icg]->restrictExtrinsicSignedCfData(css[icg]->mf,css[icg-1]->mf);
    }
    // need mu in ghosts to get face averages
    css[icg]->updateCvData(css[icg]->p); 
    css[icg]->updateCvData(css[icg]->mu_lam); 
    css[icg]->updateCvData(css[icg]->mu_sgs);
    css[icg]->updateCvData(css[icg]->a_sgs);

    // get boundary conditions on coarse grid
    if (icg == 0) {
      for (vector<HelmholtzBc*>::iterator it = css[icg]->bcs.begin(), it2 = bcs.begin(); 
          it != css[icg]->bcs.end(); ++it, ++it2) {
        assert(it2 != bcs.end());
        (*it)->restrictBc(cgs[icg],(*it2));
      }
    }
    else {
      for (vector<HelmholtzBc*>::iterator it = css[icg]->bcs.begin(), it2 = css[icg-1]->bcs.begin(); 
          it != css[icg]->bcs.end(); ++it, ++it2) {
        assert(it2 != css[icg-1]->bcs.end());
        (*it)->restrictBc(cgs[icg],(*it2));
      }
    }

    // build A
    A_cg[icg] = NULL;
    At_cg[icg] = NULL;
    css[icg]->buildUMg(A_cg[icg],At_cg[icg],NULL);
    inv_diag_cg[icg] = new double[css[icg]->ncv_g];
    for (int icv = 0; icv < css[icg]->ncv; ++icv)
      inv_diag_cg[icg][icv] = 1.0/A_cg[icg][css[icg]->cvocv_i[icv]];
    css[icg]->updateCvData(inv_diag_cg[icg]);
    //dumpRange(A_cg[icg],css[icg]->cvocv_i[css[icg]->ncv],"A_cg");
    err_cg[icg] = new double[cgs[icg]->ncc_g][3]; // MUST be ncc_g for prolongation
    res_cg[icg] = new double[cgs[icg]->ncc_g][3]; // MUST be ncc_g for restriction

    A0_cg[icg] = NULL;
    A1_cg[icg] = NULL;
    if (b_skew_decomp) {
      A0_cg[icg] = new double[css[icg]->cvocv_i[css[icg]->ncv]];
      A1_cg[icg] = new double[css[icg]->cvocv_i[css[icg]->ncv]];
      for (int icv = 0; icv < css[icg]->ncv; ++icv) {
        for (int coc = css[icg]->cvocv_i[icv]; coc != css[icg]->cvocv_i[icv+1]; ++coc) {
          A0_cg[icg][coc] = 0.5*(inv_diag_cg[icg][icv]*A_cg[icg][coc]+inv_diag_cg[icg][css[icg]->cvocv_v[coc]]*At_cg[icg][coc]);
          A1_cg[icg][coc] = 0.5*(inv_diag_cg[icg][icv]*A_cg[icg][coc]-inv_diag_cg[icg][css[icg]->cvocv_v[coc]]*At_cg[icg][coc]);
        }
      }

      double g2 = 0.0;
      double g3 = 0.0;
      for (int icv = 0; icv < css[icg]->ncv; ++icv) {
        double this_A1_inf = 0.0;
        double this_A0_inf = 0.0;
        for (int coc = css[icg]->cvocv_i[icv]; coc != css[icg]->cvocv_i[icv+1]; ++coc) {
          this_A1_inf += fabs(A1_cg[icg][coc]);
          this_A0_inf += fabs(A0_cg[icg][coc]);
        }
        g3 = max(g3,this_A1_inf);
        g2 = max(g2,this_A0_inf);
      }
      tau_bound_cg[icg] = 2.0/(2.0*g3+g2);

      if (mpi_rank == 0)
        cout << " > icg, tau_bound: " << icg << " " << tau_bound_cg[icg] << endl;

    }
  }

  // solve using V-cycle

  const double zero = getDoubleParam("MOMENTUM_ZERO", 1.0E-6);
  const int maxiter = getIntParam("MOMENTUM_MAXITER",100);
  const double relax = getDoubleParam("MOMENTUM_RELAX",1.0); 

  const double prolong_relax = getDoubleParam("MOMENTUM_MG_PROLONG_RELAX",1.0); 

  int iter = 0;
  int done = 0;
  while (done == 0) {
    iter++;

    if (momentumSmootherVec[0] == 0)
      smoothCvJacobi(u,r,inv_diag,A,rhs,momentumNsmoothVec[0],relax);
    else if (momentumSmootherVec[0] == 1)
      smoothCvGs(u,inv_diag,A,rhs,momentumNsmoothVec[0],relax);
    else if (momentumSmootherVec[0] == 2)
      smoothCvPatr(u,r,ps,Ar,Aps,inv_diag,A,At,rhs,momentumNsmoothVec[0],relax,relax);
    else 
      smoothCvTim(u,r,inv_diag,A,A1,rhs,momentumNsmoothVec[0],tau_bound,relax);

    // compute residual 

    calcCvResidual(res,u,A,rhs);

    // check if done

    double my_res_max = 0.0;
    FOR_ICV FOR_I3 my_res_max = max(my_res_max,fabs(res[icv][i]*inv_diag[icv]));
    double res_max;
    MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
    if (mpi_rank == 0) {
      if (step%check_interval==0)
        cout << " > solveU iter, res_max: " << iter << " " << res_max << endl;
      if (res_max < zero) {
        done = 1;
      }
      else if (iter > maxiter) {
        cout << " > Warning: solveU did not converge after " << maxiter << " iters, res_max: " << res_max << endl;
        done = 2;
      }
    }
    MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

    if (done != 0)
      break;

    // restrict residual 

    cgs[0]->restrictExtrinsicCcData(res_cg[0],res); 
    css[0]->updateCvData(res_cg[0]);

    for (int icg = 0; icg < ncg_u-1; ++icg) {

      //if (mpi_rank == 0)
      //  cout << " icg: " << icg << endl;

      // smooth

      for (int icc = 0; icc < css[icg]->ncv_g; ++icc) 
        FOR_I3 err_cg[icg][icc][i] = 0.0;

      if (momentumSmootherVec[icg+1] == 0)
        css[icg]->smoothCvJacobi(err_cg[icg],r,inv_diag_cg[icg],A_cg[icg],res_cg[icg],momentumNsmoothVec[icg+1],relax);
      else if (momentumSmootherVec[icg+1] == 1)
        css[icg]->smoothCvGs(err_cg[icg],inv_diag_cg[icg],A_cg[icg],res_cg[icg],momentumNsmoothVec[icg+1],relax);
      else if (momentumSmootherVec[icg+1] == 2)
        css[icg]->smoothCvPatr(err_cg[icg],r,ps,Ar,Aps,inv_diag_cg[icg],A_cg[icg],At_cg[icg],res_cg[icg],momentumNsmoothVec[icg+1],relax,relax);
      else
        css[icg]->smoothCvTim(err_cg[icg],r,inv_diag_cg[icg],A_cg[icg],A1_cg[icg],res_cg[icg],momentumNsmoothVec[icg+1],tau_bound_cg[icg],relax);

      // compute residual 

      css[icg]->calcCvResidual(r,err_cg[icg],A_cg[icg],res_cg[icg]);

      // restrict residual 

      cgs[icg+1]->restrictExtrinsicCcData(res_cg[icg+1],r);
      css[icg+1]->updateCvData(res_cg[icg+1]);

    }

    //if (mpi_rank == 0)
    //  cout << " icg: " << ncg_u-1 << endl;

    // smooth on coarsest grid

    for (int icc = 0; icc < css[ncg_u-1]->ncv_g; ++icc) 
      FOR_I3 err_cg[ncg_u-1][icc][i] = 0.0;

    if (momentumSmootherVec[ncg_u] == 0)
      css[ncg_u-1]->smoothCvJacobi(err_cg[ncg_u-1],r,inv_diag_cg[ncg_u-1],A_cg[ncg_u-1],res_cg[ncg_u-1],momentumNsmoothVec[ncg_u],relax);
    else if (momentumSmootherVec[ncg_u] == 1)
      css[ncg_u-1]->smoothCvGs(err_cg[ncg_u-1],inv_diag_cg[ncg_u-1],A_cg[ncg_u-1],res_cg[ncg_u-1],momentumNsmoothVec[ncg_u],relax);
    else if (momentumSmootherVec[ncg_u] == 2)
      css[ncg_u-1]->smoothCvPatr(err_cg[ncg_u-1],r,ps,Ar,Aps,inv_diag_cg[ncg_u-1],A_cg[ncg_u-1],At_cg[ncg_u-1],res_cg[ncg_u-1],momentumNsmoothVec[ncg_u],relax,relax);
    else 
      css[ncg_u-1]->smoothCvTim(err_cg[ncg_u-1],r,inv_diag_cg[ncg_u-1],A_cg[ncg_u-1],A1_cg[ncg_u-1],res_cg[ncg_u-1],momentumNsmoothVec[ncg_u],tau_bound_cg[ncg_u-1],relax);

    for (int icg = ncg_u-1; icg >= 1; --icg) {

      // prolong and update guess

      cgs[icg]->updateCcIDataReverse(err_cg[icg]); 
      cgs[icg]->prolongCcDataAndUpdateGuess(err_cg[icg-1],err_cg[icg],prolong_relax);
      css[icg-1]->updateCvData(err_cg[icg-1]);

      // smooth

      if (momentumSmootherVec[icg] == 0)
        css[icg-1]->smoothCvJacobi(err_cg[icg-1],r,inv_diag_cg[icg-1],A_cg[icg-1],res_cg[icg-1],momentumNsmoothVec[icg],relax);
      else if (momentumSmootherVec[icg] == 1)
        css[icg-1]->smoothCvGs(err_cg[icg-1],inv_diag_cg[icg-1],A_cg[icg-1],res_cg[icg-1],momentumNsmoothVec[icg],relax);
      else if (momentumSmootherVec[icg] == 2)
        css[icg-1]->smoothCvPatr(err_cg[icg-1],r,ps,Ar,Aps,inv_diag_cg[icg-1],A_cg[icg-1],At_cg[icg-1],res_cg[icg-1],momentumNsmoothVec[icg],relax,relax);
      else 
        css[icg-1]->smoothCvTim(err_cg[icg-1],r,inv_diag_cg[icg-1],A_cg[icg-1],A1_cg[icg-1],res_cg[icg-1],momentumNsmoothVec[icg],tau_bound_cg[icg-1],relax);

    }

    // prolong and update guess

    cgs[0]->updateCcIDataReverse(err_cg[0]); 
    cgs[0]->prolongCcDataAndUpdateGuess(u,err_cg[0],prolong_relax);
    updateCvData(u);

  }

  // cleanup

  delete[] A;
  delete[] At;
  if (b_skew_decomp) {
    delete[] A0;
    delete[] A1;
  }
  delete[] rhs;
  delete[] inv_diag;
  delete[] res;
  delete[] r;
  delete[] ps;
  delete[] Ar;
  delete[] Aps;

  for (int icg = 0; icg < ncg_u; ++icg) {
    delete[] A_cg[icg];
    delete[] At_cg[icg];
    if (b_skew_decomp) {
      delete[] A0_cg[icg];
      delete[] A1_cg[icg];
    }
    delete[] inv_diag_cg[icg];
    delete[] err_cg[icg];
    delete[] res_cg[icg];
  }

}

void HelmholtzSolver::solveUMgAdiag() {
  assert(ncg_u >= 1);
  assert(A_diag_source);

  if ( (mpi_rank==0) && (step%check_interval==0) )
    cout << "solveUMgAdiag()" << endl;

  // build lhs matrices (and modify rhs's)

  double *A = NULL;
  double *At = NULL;
  double (*rhs)[3] = new double[ncv][3];
  buildUMg(A,At,rhs);
  double (*inv_diag)[3] = new double[ncv_g][3];
  FOR_ICV FOR_I3 inv_diag[icv][i] = 1.0/(A[cvocv_i[icv]]+A_diag_source[icv][i]);
  updateCvData(inv_diag);
  double (*res)[3] = new double[ncv][3];

  for (int icg = 0; icg <= ncg_u; ++icg) {
    //if (momentumSmootherVec[icg] == 0){
    //  CERR("Momentum Smoother Option 0 (jacobi) not supported when using implicit body forces");
    //}else 
    if (momentumSmootherVec[icg] == 1){
      CERR("Momentum Smoother Option 1 (gauss-seidel) not supported when using implicit body forces");
    }else if (momentumSmootherVec[icg] == 3){
      CERR("Momentum Smoother Option 3 (triangular iterative method) not supported when using implicit body forces");
    }
  }


  // we need the following work arrays...

  double (*r)[3] = new double[ncv_g][3];
  double (*ps)[3] = new double[ncv_g][3];
  double (*Ar)[3] = new double[ncv][3];
  double (*Aps)[3] = new double[ncv][3];

  vector<double *> A_cg(ncg_u);
  vector<double *> At_cg(ncg_u);
  vector<double (*)[3]> inv_diag_cg(ncg_u);
  vector<double (*)[3]> err_cg(ncg_u);
  vector<double (*)[3]> res_cg(ncg_u);

  for (int icg = 0; icg < ncg_u; ++icg) {

    // get rho, mu_lam, mu_sgs and mf on coarse grid
    if (icg == 0) {
      cgs[icg]->restrictCcData(css[icg]->rho,rho);
      cgs[icg]->restrictCcData(css[icg]->p,p);
      cgs[icg]->restrictCcData(css[icg]->u,u);
      cgs[icg]->restrictCcData(css[icg]->mu_lam,mu_lam);
      cgs[icg]->restrictCcData(css[icg]->mu_sgs,mu_sgs);
      cgs[icg]->restrictCcData(css[icg]->a_sgs,a_sgs);
      cgs[icg]->restrictExtrinsicCcData(css[icg]->A_diag_source,A_diag_source);
      cgs[icg]->restrictExtrinsicSignedCfData(css[icg]->mf,mf);
    }
    else {
      cgs[icg]->restrictCcData(css[icg]->rho,css[icg-1]->rho);
      cgs[icg]->restrictCcData(css[icg]->p,css[icg-1]->p);
      cgs[icg]->restrictCcData(css[icg]->u,css[icg-1]->u);
      cgs[icg]->restrictCcData(css[icg]->mu_lam,css[icg-1]->mu_lam);
      cgs[icg]->restrictCcData(css[icg]->mu_sgs,css[icg-1]->mu_sgs);
      cgs[icg]->restrictCcData(css[icg]->a_sgs,css[icg-1]->a_sgs);
      cgs[icg]->restrictExtrinsicCcData(css[icg]->A_diag_source,css[icg-1]->A_diag_source);
      cgs[icg]->restrictExtrinsicSignedCfData(css[icg]->mf,css[icg-1]->mf);
    }
    // need mu in ghosts to get face averages
    css[icg]->updateCvData(css[icg]->p); 
    css[icg]->updateCvData(css[icg]->mu_lam); 
    css[icg]->updateCvData(css[icg]->mu_sgs);
    css[icg]->updateCvData(css[icg]->a_sgs);

    // get boundary conditions on coarse grid
    if (icg == 0) {
      for (vector<HelmholtzBc*>::iterator it = css[icg]->bcs.begin(), it2 = bcs.begin(); 
          it != css[icg]->bcs.end(); ++it, ++it2) {
        assert(it2 != bcs.end());
        (*it)->restrictBc(cgs[icg],(*it2));
      }
    }
    else {
      for (vector<HelmholtzBc*>::iterator it = css[icg]->bcs.begin(), it2 = css[icg-1]->bcs.begin(); 
          it != css[icg]->bcs.end(); ++it, ++it2) {
        assert(it2 != css[icg-1]->bcs.end());
        (*it)->restrictBc(cgs[icg],(*it2));
      }
    }

    // build A
    A_cg[icg] = NULL;
    At_cg[icg] = NULL;
    css[icg]->buildUMg(A_cg[icg],At_cg[icg],NULL);
    inv_diag_cg[icg] = new double[css[icg]->ncv_g][3];
    for (int icv = 0; icv < css[icg]->ncv; ++icv)
      FOR_I3 inv_diag_cg[icg][icv][i] = 1.0/(A_cg[icg][css[icg]->cvocv_i[icv]]+css[icg]->A_diag_source[icv][i]);
    css[icg]->updateCvData(inv_diag_cg[icg]);
    //dumpRange(A_cg[icg],css[icg]->cvocv_i[css[icg]->ncv],"A_cg");
    err_cg[icg] = new double[cgs[icg]->ncc_g][3]; // MUST be ncc_g for prolongation
    res_cg[icg] = new double[cgs[icg]->ncc_g][3]; // MUST be ncc_g for restriction

  }

  // solve using V-cycle

  const double zero = getDoubleParam("MOMENTUM_ZERO", 1.0E-6);
  const int maxiter = getIntParam("MOMENTUM_MAXITER",100);
  const double relax = getDoubleParam("MOMENTUM_RELAX",1.0); 

  const double prolong_relax = getDoubleParam("MOMENTUM_MG_PROLONG_RELAX",1.0); 

  int iter = 0;
  int done = 0;
  while (done == 0) {
    iter++;

    if (momentumSmootherVec[0] == 0)
      smoothCvJacobi(u,r,inv_diag,A,A_diag_source,rhs,momentumNsmoothVec[0],relax);
    else if (momentumSmootherVec[0] == 1)
      assert(0); //smoothCvGs(u,inv_diag,A,rhs,momentumNsmoothVec[0],relax);
    else if (momentumSmootherVec[0] == 2)
      smoothCvPatr(u,r,ps,Ar,Aps,inv_diag,A,At,A_diag_source,rhs,momentumNsmoothVec[0],relax,relax);
    else
      assert(0); //smoothCvTim(u,r,inv_diag,A,A1,rhs,momentumNsmoothVec[0],tau_bound,relax);

    // compute residual 

    calcCvResidual(res,u,A,A_diag_source,rhs);

    // check if done

    double my_res_max = 0.0;
    FOR_ICV FOR_I3 my_res_max = max(my_res_max,fabs(res[icv][i]*inv_diag[icv][i]));
    double res_max;
    MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
    if (mpi_rank == 0) {
      if (step%check_interval==0)
        cout << " > solveU iter, res_max: " << iter << " " << res_max << endl;
      if (res_max < zero) {
        done = 1;
      }
      else if (iter > maxiter) {
        cout << " > Warning: solveU did not converge after " << maxiter << " iters, res_max: " << res_max << endl;
        done = 2;
      }
    }
    MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

    if (done != 0)
      break;

    // restrict residual 

    cgs[0]->restrictExtrinsicCcData(res_cg[0],res); 
    css[0]->updateCvData(res_cg[0]);

    for (int icg = 0; icg < ncg_u-1; ++icg) {

      //if (mpi_rank == 0)
      //  cout << " icg: " << icg << endl;

      // smooth

      for (int icc = 0; icc < css[icg]->ncv_g; ++icc) 
        FOR_I3 err_cg[icg][icc][i] = 0.0;

      if (momentumSmootherVec[icg+1] == 0)
         css[icg]->smoothCvJacobi(err_cg[icg],r,inv_diag_cg[icg],A_cg[icg],css[icg]->A_diag_source,res_cg[icg],momentumNsmoothVec[icg+1],relax);
      else if (momentumSmootherVec[icg+1] == 1)
         assert(0); //css[icg]->smoothCvGs(err_cg[icg],inv_diag_cg[icg],A_cg[icg],res_cg[icg],momentumNsmoothVec[icg+1],relax);
      else if (momentumSmootherVec[icg+1] == 2)
        css[icg]->smoothCvPatr(err_cg[icg],r,ps,Ar,Aps,inv_diag_cg[icg],A_cg[icg],At_cg[icg],css[icg]->A_diag_source,res_cg[icg],momentumNsmoothVec[icg+1],relax,relax);
      else
         assert(0); //css[icg]->smoothCvTim(err_cg[icg],r,inv_diag_cg[icg],A_cg[icg],A1_cg[icg],res_cg[icg],momentumNsmoothVec[icg+1],tau_bound_cg[icg],relax);

      // compute residual 

      css[icg]->calcCvResidual(r,err_cg[icg],A_cg[icg],css[icg]->A_diag_source,res_cg[icg]);

      // restrict residual 

      cgs[icg+1]->restrictExtrinsicCcData(res_cg[icg+1],r);
      css[icg+1]->updateCvData(res_cg[icg+1]);

    }

    //if (mpi_rank == 0)
    //  cout << " icg: " << ncg_u-1 << endl;


    // smooth on coarsest grid

    for (int icc = 0; icc < css[ncg_u-1]->ncv_g; ++icc) 
      FOR_I3 err_cg[ncg_u-1][icc][i] = 0.0;

    if (momentumSmootherVec[ncg_u] == 0)
       css[ncg_u-1]->smoothCvJacobi(err_cg[ncg_u-1],r,inv_diag_cg[ncg_u-1],A_cg[ncg_u-1],css[ncg_u-1]->A_diag_source,res_cg[ncg_u-1],momentumNsmoothVec[ncg_u],relax);
    else if (momentumSmootherVec[ncg_u] == 1)
       assert(0); //css[ncg_u-1]->smoothCvGs(err_cg[ncg_u-1],inv_diag_cg[ncg_u-1],A_cg[ncg_u-1],res_cg[ncg_u-1],momentumNsmoothVec[ncg_u],relax);
    else if (momentumSmootherVec[ncg_u] == 2)
      css[ncg_u-1]->smoothCvPatr(err_cg[ncg_u-1],r,ps,Ar,Aps,inv_diag_cg[ncg_u-1],A_cg[ncg_u-1],At_cg[ncg_u-1],css[ncg_u-1]->A_diag_source,res_cg[ncg_u-1],momentumNsmoothVec[ncg_u],relax,relax);
    else 
       assert(0); //css[ncg_u-1]->smoothCvTim(err_cg[ncg_u-1],r,inv_diag_cg[ncg_u-1],A_cg[ncg_u-1],A1_cg[ncg_u-1],res_cg[ncg_u-1],momentumNsmoothVec[ncg_u],tau_bound_cg[ncg_u-1],relax);

    for (int icg = ncg_u-1; icg >= 1; --icg) {

      // prolong and update guess

      cgs[icg]->updateCcIDataReverse(err_cg[icg]); 
      cgs[icg]->prolongCcDataAndUpdateGuess(err_cg[icg-1],err_cg[icg],prolong_relax);
      css[icg-1]->updateCvData(err_cg[icg-1]);

      // smooth

      if (momentumSmootherVec[icg] == 0)
         css[icg-1]->smoothCvJacobi(err_cg[icg-1],r,inv_diag_cg[icg-1],A_cg[icg-1],css[icg-1]->A_diag_source,res_cg[icg-1],momentumNsmoothVec[icg],relax);
      else if (momentumSmootherVec[icg] == 1)
         assert(0); //css[icg-1]->smoothCvGs(err_cg[icg-1],inv_diag_cg[icg-1],A_cg[icg-1],res_cg[icg-1],momentumNsmoothVec[icg],relax);
      else if (momentumSmootherVec[icg] == 2)
         css[icg-1]->smoothCvPatr(err_cg[icg-1],r,ps,Ar,Aps,inv_diag_cg[icg-1],A_cg[icg-1],At_cg[icg-1],css[icg-1]->A_diag_source,res_cg[icg-1],momentumNsmoothVec[icg],relax,relax);
      else 
         assert(0); //css[icg-1]->smoothCvTim(err_cg[icg-1],r,inv_diag_cg[icg-1],A_cg[icg-1],A1_cg[icg-1],res_cg[icg-1],momentumNsmoothVec[icg],tau_bound_cg[icg-1],relax);

    }

    // prolong and update guess

    cgs[0]->updateCcIDataReverse(err_cg[0]); 
    cgs[0]->prolongCcDataAndUpdateGuess(u,err_cg[0],prolong_relax);
    updateCvData(u);

  }


  // cleanup

  delete[] A;
  delete[] At;
  delete[] rhs;
  delete[] inv_diag;
  delete[] res;
  delete[] r;
  delete[] ps;
  delete[] Ar;
  delete[] Aps;

  for (int icg = 0; icg < ncg_u; ++icg) {
    delete[] A_cg[icg];
    delete[] At_cg[icg];
    delete[] inv_diag_cg[icg];
    delete[] err_cg[icg];
    delete[] res_cg[icg];
  }

}


void HelmholtzSolver::buildUMg(double * &A,double * &At, double (*rhs)[3]) {
  assert(A == NULL);
  assert(At == NULL);

  const int cvocv_s = cvocv_i[ncv];
  A = new double[cvocv_s];
  for (int coc = 0; coc < cvocv_s; ++coc)
    A[coc] = 0.0;
  At = new double[cvocv_s];
  for (int coc = 0; coc < cvocv_s; ++coc)
    At[coc] = 0.0;

  if (rhs != NULL){
    FOR_ICV FOR_I3 rhs[icv][i] = -dpdx[icv][i]*vol_cv[icv];
    addFrameRotationSourceTerms(rhs);
  }

  // bcs...
  FOR_BCZONE (*it)->addMomentumFlux(A,At,rhs);

  const double inv_sos2 = 1.0/(helmholtz_sos*helmholtz_sos);

  // internal faces...
  for (int ifa = 0; ifa < nfa; ++ifa){

    const int icv0 = cvofa[ifa][0];
    assert((icv0 >= 0)&&(icv0 < ncv));
    const int icv1 = cvofa[ifa][1];
    assert((icv1 >= 0)&&(icv0 < ncv_g));

    const int coc00 = cvocv_i[icv0];
    int coc01 = coc00+1;
    while (cvocv_v[coc01] != icv1) {
      ++coc01;
      assert(coc01 != cvocv_i[icv0+1] );
    }

    const double mu_coeff = 0.5*(mu_lam[icv0] + mu_lam[icv1] + mu_sgs[icv0] + mu_sgs[icv1])*area_over_delta_fa[ifa]; 

    const double asgs_coeff = 0.5*(a_sgs[icv0] + a_sgs[icv1])*inv_sos2*(p[icv1]-p[icv0])*area_over_delta_fa[ifa];
   
    A[coc00]  +=  0.5*(mf[ifa] - asgs_coeff) + mu_coeff;
    A[coc01]  +=  0.5*(mf[ifa] - asgs_coeff) - mu_coeff;
    At[coc00] +=  0.5*(mf[ifa] - asgs_coeff) + mu_coeff;
    At[coc01] += -0.5*(mf[ifa] - asgs_coeff) - mu_coeff;

    if (icv1 < ncv) {

      const int coc11 = cvocv_i[icv1];
      int coc10 = coc11+1;
      while (cvocv_v[coc10] != icv0) {
        ++coc10;
        assert(coc10 != cvocv_i[icv1+1] );
      }

      A[coc11]  -=  0.5*(mf[ifa] - asgs_coeff) - mu_coeff;
      A[coc10]  -=  0.5*(mf[ifa] - asgs_coeff) + mu_coeff;
      At[coc11] -=  0.5*(mf[ifa] - asgs_coeff) - mu_coeff;
      At[coc10] -= -0.5*(mf[ifa] - asgs_coeff) + mu_coeff;

    }

  }

  // time term...

  FOR_ICV A[cvocv_i[icv]] += 1.5*vol_cv[icv]*rho[icv]/dt;
  if (rhs != NULL)
    FOR_ICV FOR_I3 rhs[icv][i] += vol_cv[icv]*( 2.0*rho0[icv]*u0[icv][i] - 0.5*rho00[icv]*u00[icv][i] )/dt;

  // source terms...

  momentumSourceHook(A,At,rhs);
  if ((rhou_source)&&(rhs != NULL)) {
    FOR_ICV FOR_I3 {
      rhs[icv][i] += rhou_source[icv][i];
    }
  }

}

void HelmholtzSolver::solveScalars() {

  if ( nsc_transport == 0)
    return;

  const double zero = getDoubleParam("SCALAR_ZERO", 1.0E-6);
  const double relax = getDoubleParam("SCALAR_RELAX",0.7);
  const int maxiter  = getIntParam("SCALAR_MAXITER", 1000);

  if ( (mpi_rank == 0) && (step%check_interval == 0) )
    cout << "solveScalars()" << endl;

  // solve...
  // [A]*{phi} = {rhs_phi} for each phi
  
  // uses SDIRK2 scheme...

  const int cvocv_s = cvocv_i[ncv];
  double * A = new double[cvocv_s];
  double * rhs = new double[ncv];
  double * phi0 = new double[ncv_g];
  const double inv_sos2 = 1.0/(helmholtz_sos*helmholtz_sos);
  const double gamma = 1.0-0.5*sqrt(2.0);
  const double inv_gamma_dt = 1.0/(gamma*dt);
  for (map<string,ScalarField*>::iterator it = scalar_map.begin(); it != scalar_map.end(); ++it) { 

    const ScalarField * scalar = it->second;
    const int isc = scalar->idx; 
    const double Sc_lam = scalar->Sc_lam;
    const double Sc_t   = scalar->Sc_t;

    for (int coc = 0; coc < cvocv_s; ++coc)
      A[coc] = 0.0;

    for (int ifa = 0; ifa < nfa; ++ifa) {

      const int icv0 = cvofa[ifa][0];
      assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1];
      assert((icv1 >= 0)&&(icv0 < ncv_g));

      const int coc00 = cvocv_i[icv0];
      int coc01 = coc00+1;
      while (cvocv_v[coc01] != icv1) {
        ++coc01;
        assert(coc01 != cvocv_i[icv0+1] );
      }

      const double D_coeff = 0.5*(Sc_lam*(mu_lam[icv0] + mu_lam[icv1]) + Sc_t*(mu_sgs[icv0] + mu_sgs[icv1]))*area_over_delta_fa[ifa]; 
      const double asgs_coeff = 0.5*(a_sgs[icv0] + a_sgs[icv1])*inv_sos2*(p[icv1]-p[icv0])*area_over_delta_fa[ifa];

      A[coc00] += 0.5*(mf[ifa] - asgs_coeff) + D_coeff;
      A[coc01] += 0.5*(mf[ifa] - asgs_coeff) - D_coeff;

      if (icv1 < ncv) {

        const int coc11 = cvocv_i[icv1];
        int coc10 = coc11+1;
        while (cvocv_v[coc10] != icv0) {
          ++coc10;
          assert(coc10 != cvocv_i[icv1+1] );
        }

        A[coc11] -= 0.5*(mf[ifa] - asgs_coeff) - D_coeff;
        A[coc10] -= 0.5*(mf[ifa] - asgs_coeff) + D_coeff;

      }

    }

    // first implicit step 0->gamma...
   
    FOR_ICV_G phi0[icv] = transport_scalar_vec[isc][icv];

    // time term...

    FOR_ICV {
      const double rho_tmp = (1.0-gamma)*rho0[icv]+gamma*rho[icv]; // linear interpolation
      A[cvocv_i[icv]] += vol_cv[icv]*rho_tmp*inv_gamma_dt;
      rhs[icv] = vol_cv[icv]*rho0[icv]*phi0[icv]*inv_gamma_dt;
    }

    // add boundary fluxes...

    FOR_BCZONE (*it)->addScalarBoundaryFlux(A,rhs,scalar);

    // and solve...

    solveCvJacobi(transport_scalar_vec[isc],A,rhs,zero,relax,maxiter,step%check_interval==0); 

    // extrapolation gamma->1-gamma...

    FOR_ICV_G {
      const double rho_tmp0 = (1.0-gamma)*rho0[icv]+gamma*rho[icv];
      const double rho_tmp1 = gamma*rho0[icv]+(1.0-gamma)*rho[icv];
      transport_scalar_vec[isc][icv] = ((1.0-gamma)*rho_tmp0*transport_scalar_vec[isc][icv]-(1.0-2.0*gamma)*rho0[icv]*phi0[icv])/(rho_tmp1*gamma);
    }

    // second implicit step 1-gamma->1...

    // time term...
    
    FOR_ICV {
      const double rho_tmp0 = (1.0-gamma)*rho0[icv]+gamma*rho[icv];
      const double rho_tmp1 = gamma*rho0[icv]+(1.0-gamma)*rho[icv];
      A[cvocv_i[icv]] += vol_cv[icv]*(rho[icv]-rho_tmp0)*inv_gamma_dt;
      rhs[icv] = vol_cv[icv]*rho_tmp1*transport_scalar_vec[isc][icv]*inv_gamma_dt;
    }

    // add boundary fluxes...

    FOR_BCZONE (*it)->addScalarBoundaryFlux(NULL,rhs,scalar);

    // and solve...

    solveCvJacobi(transport_scalar_vec[isc],A,rhs,zero,relax,maxiter,step%check_interval==0); 

  }

  // cleanup...

  delete[] A;
  delete[] rhs;
  delete[] phi0;

}

void HelmholtzSolver::solveScalarsPatr() {

  if ( nsc_transport == 0)
    return;

  const double zero = getDoubleParam("SCALAR_ZERO",1.0E-6);
  const double relax = getDoubleParam("SCALAR_RELAX",1.0);
  const int maxiter  = getIntParam("SCALAR_MAXITER",1000);

  if ( (mpi_rank == 0) && (step%check_interval == 0) )
    cout << "solveScalarsPatr()" << endl;

  // solve...
  // [A]*{phi} = {rhs_phi} for each phi
  
  // uses SDIRK2 scheme assuming frozen rho...

  const int cvocv_s = cvocv_i[ncv];
  double * A = new double[cvocv_s];
  double * At = new double[cvocv_s];
  double * rhs = new double[ncv];
  double * phi0 = new double[ncv_g];
  const double inv_sos2 = 1.0/(helmholtz_sos*helmholtz_sos);
  const double gamma = 1.0-0.5*sqrt(2.0);
  const double inv_gamma_dt = 1.0/(gamma*dt);
  for (map<string,ScalarField*>::iterator it = scalar_map.begin(); it != scalar_map.end(); ++it) { 

    const ScalarField * scalar = it->second;
    const int isc = scalar->idx; 
    const double Sc_lam = scalar->Sc_lam;
    const double Sc_t   = scalar->Sc_t;

    for (int coc = 0; coc < cvocv_s; ++coc)
      A[coc] = At[coc] = 0.0;

    for (int ifa = 0; ifa < nfa; ++ifa) {

      const int icv0 = cvofa[ifa][0];
      assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1];
      assert((icv1 >= 0)&&(icv0 < ncv_g));

      const int coc00 = cvocv_i[icv0];
      int coc01 = coc00+1;
      while (cvocv_v[coc01] != icv1) {
        ++coc01;
        assert(coc01 != cvocv_i[icv0+1] );
      }

      const double D_coeff = 0.5*(Sc_lam*(mu_lam[icv0] + mu_lam[icv1]) + Sc_t*(mu_sgs[icv0] + mu_sgs[icv1]))*area_over_delta_fa[ifa]; 
      const double asgs_coeff = 0.5*(a_sgs[icv0] + a_sgs[icv1])*inv_sos2*(p[icv1]-p[icv0])*area_over_delta_fa[ifa];

      A[coc00]  +=  0.5*(mf[ifa] - asgs_coeff) + D_coeff;
      A[coc01]  +=  0.5*(mf[ifa] - asgs_coeff) - D_coeff;
      At[coc00] +=  0.5*(mf[ifa] - asgs_coeff) + D_coeff;
      At[coc01] += -0.5*(mf[ifa] - asgs_coeff) - D_coeff;

      if (icv1 < ncv) {

        const int coc11 = cvocv_i[icv1];
        int coc10 = coc11+1;
        while (cvocv_v[coc10] != icv0) {
          ++coc10;
          assert(coc10 != cvocv_i[icv1+1] );
        }

        A[coc11]  -=  0.5*(mf[ifa] - asgs_coeff) - D_coeff;
        A[coc10]  -=  0.5*(mf[ifa] - asgs_coeff) + D_coeff;
        At[coc11] -=  0.5*(mf[ifa] - asgs_coeff) - D_coeff;
        At[coc10] -= -0.5*(mf[ifa] - asgs_coeff) + D_coeff;

      }

    }

    // first implicit step 0->gamma...
   
    FOR_ICV_G phi0[icv] = transport_scalar_vec[isc][icv];

    // time term...

    FOR_ICV {
      const double rho_tmp = (1.0-gamma)*rho0[icv]+gamma*rho[icv]; // linear interpolation
      A[cvocv_i[icv]]  += vol_cv[icv]*rho_tmp*inv_gamma_dt;
      At[cvocv_i[icv]] += vol_cv[icv]*rho_tmp*inv_gamma_dt;
      rhs[icv] = vol_cv[icv]*rho0[icv]*phi0[icv]*inv_gamma_dt;
    }

    // add boundary fluxes...

    FOR_BCZONE (*it)->addScalarBoundaryFlux(A,At,rhs,scalar);

    // and solve...

    solveCvPatr(transport_scalar_vec[isc],A,At,rhs,zero,relax,relax,maxiter,step%check_interval==0); 

    // extrapolation gamma->1-gamma...

    FOR_ICV_G {
      const double rho_tmp0 = (1.0-gamma)*rho0[icv]+gamma*rho[icv];
      const double rho_tmp1 = gamma*rho0[icv]+(1.0-gamma)*rho[icv];
      transport_scalar_vec[isc][icv] = ((1.0-gamma)*rho_tmp0*transport_scalar_vec[isc][icv]-(1.0-2.0*gamma)*rho0[icv]*phi0[icv])/(rho_tmp1*gamma);
    }

    // second implicit step 1-gamma->1...

    // time term...
    
    FOR_ICV {
      const double rho_tmp0 = (1.0-gamma)*rho0[icv]+gamma*rho[icv];
      const double rho_tmp1 = gamma*rho0[icv]+(1.0-gamma)*rho[icv];
      A[cvocv_i[icv]]  += vol_cv[icv]*(rho[icv]-rho_tmp0)*inv_gamma_dt; // update diagonal...
      At[cvocv_i[icv]] += vol_cv[icv]*(rho[icv]-rho_tmp0)*inv_gamma_dt; // update diagonal...
      rhs[icv] = vol_cv[icv]*rho_tmp1*transport_scalar_vec[isc][icv]*inv_gamma_dt;
    }

    // add boundary fluxes...

    FOR_BCZONE (*it)->addScalarBoundaryFlux(NULL,NULL,rhs,scalar);

    // and solve...

    solveCvPatr(transport_scalar_vec[isc],A,At,rhs,zero,relax,relax,maxiter,step%check_interval==0); 

  }

  // cleanup...

  delete[] A;
  delete[] At;
  delete[] rhs;
  delete[] phi0;

}

void HelmholtzSolver::solvePAndCorrectU() {

  if (mpi_rank==0&&step%check_interval==0)
    cout << "solvePAndCorrectU()" << endl;

  // temporarily store rho*phi in all scalars...
  // here it does not matter whether we are using the M*rho*phi formulation
  // or the V*rho*phi formulation...


  // use rhou^* to compute the new mf...
 
 
  FOR_ICV FOR_I3 u[icv][i] = u[icv][i] + 2.0/3.0*dt*dpdx[icv][i]/rho[icv];
  updateCvData(u,REPLACE_ROTATE_DATA);

  // bc mf update expects u, not rhou
  FOR_BCZONE (*it)->updateMf();
  timer.split("bczone->updateMf");

  // temporarily store rhou^* in u...
  FOR_ICV FOR_I3 u[icv][i] *= rho[icv];
  updateCvData(u,REPLACE_ROTATE_DATA);

  // use rhou^* to compute the new mf...

  for (int ifa = 0; ifa<nfa; ++ifa){
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];

    mf[ifa] = 0.5*( n_fa[ifa][0]*(u[icv0][0]+u[icv1][0]) +
                    n_fa[ifa][1]*(u[icv0][1]+u[icv1][1]) +
                    n_fa[ifa][2]*(u[icv0][2]+u[icv1][2]) );

    if (body_force_flag){
      //faces along MRF interface get half a mass flux correction, but if the
      //grid is setup properly then DOT_PRODUCT(n_fa[ifa],u_fa) should be zero...
      BodyForce * body_force0 = NULL;
      if (body_force_flag[icv0] >= 0){
        body_force0 = &bodyForceVec[body_force_flag[icv0]];
        if (body_force0->type == MRF_TYPE){
          const double r[3] = DIFF(x_fa[ifa],body_force0->x_rot);
          //face rotational velocity: omega x r
          double u_fa[3];
          u_fa[0] = body_force0->omega_rot*(body_force0->axis_rot[1]*r[2] - body_force0->axis_rot[2]*r[1]);
          u_fa[1] = body_force0->omega_rot*(body_force0->axis_rot[2]*r[0] - body_force0->axis_rot[0]*r[2]);
          u_fa[2] = body_force0->omega_rot*(body_force0->axis_rot[0]*r[1] - body_force0->axis_rot[1]*r[0]);
          mf[ifa] -= 0.5*rho[icv0]*DOT_PRODUCT(n_fa[ifa],u_fa);
        }
      }
      BodyForce * body_force1 = NULL;
      if (body_force_flag[icv1] >= 0){
        body_force1 = &bodyForceVec[body_force_flag[icv1]];
        if (body_force1->type == MRF_TYPE){
          const double r[3] = DIFF(x_fa[ifa],body_force1->x_rot);
          //face rotational velocity: omega x r
          double u_fa[3];
          u_fa[0] = body_force1->omega_rot*(body_force1->axis_rot[1]*r[2] - body_force1->axis_rot[2]*r[1]);
          u_fa[1] = body_force1->omega_rot*(body_force1->axis_rot[2]*r[0] - body_force1->axis_rot[0]*r[2]);
          u_fa[2] = body_force1->omega_rot*(body_force1->axis_rot[0]*r[1] - body_force1->axis_rot[1]*r[0]);
          mf[ifa] -= 0.5*rho[icv1]*DOT_PRODUCT(n_fa[ifa],u_fa);
        }
      }
    }
  }

  double * rhs = new double[ncv];
  calcContinuity(rhs);

  // modify continuity so it is based on the uncorrected rho...

  FOR_ICV rhs[icv] += 1.5*vol_cv[icv]*(rho_uncorrected[icv]-rho[icv])/dt;

  // mass sources

  massSourceHook(rhs);

  // scale rhs...

  FOR_ICV rhs[icv] *= 1.5/dt;

  // now solve for p...

  if (ncg_p > 0) {
    solveHelmholtzMg(p,rhs);
  }
  else {
    solveHelmholtz(p,rhs);
    //solveHelmholtzPatr(p,rhs); // HACK for testing MG
  }

  // update the bc's...
  FOR_BCZONE (*it)->updateP();

  // correct density...

  const double sos2 = helmholtz_sos*helmholtz_sos;
  FOR_ICV rho[icv] = rho_uncorrected[icv] + p[icv]/sos2;
  updateCvData(rho); 

  // correct mf...

  for (int ifa = 0; ifa < nfa; ++ifa) {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    mf[ifa] -= 2.0/3.0*dt*area_over_delta_fa[ifa]*(p[icv1]-p[icv0]);
  }

  //any BC that has addPressureFlux should also correct here...
  FOR_BCZONE (*it)->correctMf();

  calcPressureGrad(dpdx,p);

  // and correct u^{*} to u^{n+1}...

  FOR_ICV FOR_I3 u[icv][i] = (u[icv][i] - 2.0/3.0*dt*dpdx[icv][i])/rho[icv];
  updateCvData(u,REPLACE_ROTATE_DATA);

  //DAPDAP
  FOR_ICV rhs[icv] = 0.0;
  calcContinuity(rhs,true);
  FOR_BCZONE {
    Param * param = getParam((*it)->getName());
    if ( param->getString(0) == "SPONGE") {
      SpongeHBc * sponge = dynamic_cast<SpongeHBc*>(*it);
      sponge->addContSource(rhs);
    }
  }
  if (step%check_interval==0)
    dumpRange(rhs,ncv,"continuity");

  delete[] rhs;

}

void HelmholtzSolver::calcCvGrad(double (*dudx)[3][3],const double (*u)[3]) {

  FOR_ICV FOR_I3 FOR_J3 dudx[icv][i][j] = 0.0;

  FOR_IFA {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    FOR_I3 FOR_J3 dudx[icv0][i][j] += n_fa[ifa][j]*0.5*(u[icv1][i]-u[icv0][i]);
    if (icv1 < ncv)
      FOR_I3 FOR_J3 dudx[icv1][i][j] -= n_fa[ifa][j]*0.5*(u[icv0][i]-u[icv1][i]);
  }

  FOR_ICV FOR_I3 FOR_J3 dudx[icv][i][j] /= vol_cv[icv];

}

void HelmholtzSolver::calcPressureGrad(double (*dpdx)[3],const double * p) {

  // the calculation of the cv-based pressure gradient
  // active must be carefully closed at boundaries...

  FOR_ICV FOR_I3 dpdx[icv][i] = 0.0;

  for (int ifa = 0; ifa < nfa; ++ifa ) {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    const double flux[3] = {
      n_fa[ifa][0]*0.5*(p[icv0]+p[icv1]),
      n_fa[ifa][1]*0.5*(p[icv0]+p[icv1]),
      n_fa[ifa][2]*0.5*(p[icv0]+p[icv1])
    };

    FOR_I3 dpdx[icv0][i] += flux[i];
    if (icv1 < ncv)
      FOR_I3 dpdx[icv1][i] -= flux[i];
  }

  FOR_BCZONE (*it)->completePressureGrad(dpdx);

  FOR_ICV FOR_I3 dpdx[icv][i] *= inv_vol[icv];

}

double HelmholtzSolver::calcCfl(double* cfl_, const double dt_) const {

  // functions returns the rank-local cfl maximum..

  bool b_memflag = false;
  double * cfl   = NULL ;
  if ( cfl_ == NULL) {
    cfl       = new double[ncv];
    b_memflag = true;
  } else {
    cfl = cfl_;
  }


  FOR_ICV cfl[icv] = 0.0;

  for (int ibf = 0; ibf < nbf; ++ibf) { 
    const int icv = cvobf[ibf];
    cfl[icv] += rho[icv]*abs(DOT_PRODUCT(u[icv],n_bf[ibf]));
  }

  for (int ifa = 0; ifa<nfa; ++ifa){
    const int icv0 = cvofa[ifa][0];
    cfl[icv0] += fabs(mf[ifa]);
    const int icv1 = cvofa[ifa][1];
    if (icv1 < ncv)
      cfl[icv1] += fabs(mf[ifa]);
  }

  FOR_ICV cfl[icv] *= 0.5*dt_/rho[icv]/vol_cv[icv];

  double my_cfl_max = 0.0;
  FOR_ICV my_cfl_max = max(my_cfl_max,cfl[icv]);

  if ( b_memflag) delete[] cfl;
  return my_cfl_max;
}


void HelmholtzSolver::computeBodyForces() {
    parseBodyForces(); //required, call this first

  if (body_force_flag){
    if ((step%check_interval == 0)&&(mpi_rank == 0))
      cout << "Computing body force terms..." << endl;

    if (rhou_source == NULL) {
      int myHookFlag = 0;
      rhou_source = new double[ncv][3];
      int my_implicit_init_flag = 0;
      FOR_ICV {
        FOR_I3 rhou_source[icv][i] = 0.0;
        if (my_implicit_init_flag<=0 && body_force_flag[icv] >= 0) {
          BodyForce * body_force = &bodyForceVec[body_force_flag[icv]];
          if (body_force->type == PMM_TYPE){
            ++my_implicit_init_flag;
          }
          else if (body_force->type == HOOK_TYPE){  //Check for user defined body forces
            myHookFlag = 1;
          }
        }  
      }

      int hookFlag;
      MPI_Allreduce(&myHookFlag,&hookFlag,1,MPI_INT,MPI_SUM,mpi_comm);
      if (hookFlag>0){
        b_bodyForceHook = true;
        initBodyForceHook(my_implicit_init_flag);  //Check for user defined body forces
      }

      int implicit_init_flag;
      MPI_Allreduce(&my_implicit_init_flag,&implicit_init_flag,1,MPI_INT,MPI_SUM,mpi_comm);
      if (implicit_init_flag>0){
        COUT1("...using partially implicit porous media model formulation.");
        assert(!A_diag_source);
        A_diag_source = new double[ncv][3];
        FOR_ICV FOR_I3 A_diag_source[icv][i] = 0.0;
      } 
      else{
        COUT1("...using fully explicit body force formulation.");
      }

      //Report some info for MRF type
      FOR_BCZONE {
        HelmholtzBc * hbc = (*it);
        int bf_count=0;
        for (int ibf=0; ibf<hbc->zone_ptr->nbf; ++ibf){
          int icv0 = hbc->zone_ptr->cvobf[ibf];
          if (body_force_flag[icv0] >= 0) {
            BodyForce * body_force = &bodyForceVec[body_force_flag[icv0]];
            if (body_force->type == MRF_TYPE){
              ++bf_count;
            }
          }
        }
        int my_buf[2] = {bf_count,hbc->zone_ptr->nbf};
        int buf[2];
        MPI_Reduce(my_buf,buf,2,MPI_INT,MPI_SUM,0,mpi_comm);
        if (mpi_rank==0&&buf[0]>0)
          cout << "  MRF INIT: Boundary zone " << hbc->getName() << " has " << buf[0] << " of " << buf[1] << " faces bordering a MRF region" << endl;
      }
      
      //sum the net volume flow rate out of all MRF regions
      double my_qdot = 0.0;
      FOR_IFA{
        const int icv0 = cvofa[ifa][0];
        const int icv1 = cvofa[ifa][1];

        if (body_force_flag[icv0]<0&&body_force_flag[icv1]<0){
          continue;
        }
        else if (body_force_flag[icv0]>=0||body_force_flag[icv1]>=0){
          if (body_force_flag[icv0]==body_force_flag[icv1]){
            continue;
          }
          else{
            //At this point icv0 and icv1 may either both be in different
            //body force regions (strange case of bordering body force regions)
            // or only one is in a body force region
            BodyForce * body_force0 = NULL;
            if (body_force_flag[icv0] >= 0){
              body_force0 = &bodyForceVec[body_force_flag[icv0]];
              if (body_force0->type == MRF_TYPE){
                const double r[3] = DIFF(x_fa[ifa],body_force0->x_rot);
                //face rotational velocity: omega x r
                double u_fa[3];
                u_fa[0] = body_force0->omega_rot*(body_force0->axis_rot[1]*r[2] - body_force0->axis_rot[2]*r[1]);
                u_fa[1] = body_force0->omega_rot*(body_force0->axis_rot[2]*r[0] - body_force0->axis_rot[0]*r[2]);
                u_fa[2] = body_force0->omega_rot*(body_force0->axis_rot[0]*r[1] - body_force0->axis_rot[1]*r[0]);
                my_qdot += 0.5*DOT_PRODUCT(n_fa[ifa],u_fa);
              }
            }
            BodyForce * body_force1 = NULL;
            if (body_force_flag[icv1] >= 0){
              body_force1 = &bodyForceVec[body_force_flag[icv1]];
              if (body_force1->type == MRF_TYPE){
                const double r[3] = DIFF(x_fa[ifa],body_force1->x_rot);
                //face rotational velocity: omega x r
                double u_fa[3];
                u_fa[0] = body_force1->omega_rot*(body_force1->axis_rot[1]*r[2] - body_force1->axis_rot[2]*r[1]);
                u_fa[1] = body_force1->omega_rot*(body_force1->axis_rot[2]*r[0] - body_force1->axis_rot[0]*r[2]);
                u_fa[2] = body_force1->omega_rot*(body_force1->axis_rot[0]*r[1] - body_force1->axis_rot[1]*r[0]);
                my_qdot -= 0.5*DOT_PRODUCT(n_fa[ifa],u_fa);
              }
            }
          }
        }
      }
      double qdot;
      MPI_Reduce(&my_qdot,&qdot,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      if (mpi_rank==0&&qdot>0)
        cout << "  MRF INIT: Net flow rate (volume/time) due to rotation through MRF interface " << qdot << endl;
    }
    //end one time if (rhou==NULL)

    if (A_diag_source){
      //clear A_diag_source if present
      FOR_ICV FOR_I3 A_diag_source[icv][i] = 0.0;
    }

    BodyForce * body_force;
    FOR_ICV {
      FOR_I3 rhou_source[icv][i] = 0.0;
      if (body_force_flag[icv] >= 0) {
        body_force = &bodyForceVec[body_force_flag[icv]];
        if (body_force->type == PMM_TYPE) {
          const double umag = MAG(u[icv]);
          double f_cv[3] = {0.0,0.0,0.0};
          FOR_I3 {
            A_diag_source[icv][i] += vol_cv[icv]*(0.5*rho[icv]*umag*body_force->C[i][i]+mu_lam[icv]*body_force->D[i][i]);
            FOR_J3 { 
              if (i!=j){
                f_cv[i] -= vol_cv[icv]*u[icv][j]*(0.5*rho[icv]*umag*body_force->C[i][j]+mu_lam[icv]*body_force->D[i][j]);
              }
            }
            // these may need to be cleared before function is called...
            rhou_source[icv][i]  += f_cv[i];
            f_cv[i] -= A_diag_source[icv][i]*u[icv][i]; 
            body_force->force[i] -= f_cv[i];
            body_force->work     -= f_cv[i]*u[icv][i];
          }
          body_force->volume   += vol_cv[icv]; 
          double m_cv[3] = CROSS_PRODUCT(x_cv[icv],f_cv);
          FOR_I3 body_force->moment[i] -= m_cv[i];
        }
        else if (body_force->type == MRF_TYPE) {
          //omega x u
          double omega_X_r[3];
          omega_X_r[0] = body_force->omega_rot*(body_force->axis_rot[1]*u[icv][2] - body_force->axis_rot[2]*u[icv][1]);
          omega_X_r[1] = body_force->omega_rot*(body_force->axis_rot[2]*u[icv][0] - body_force->axis_rot[0]*u[icv][2]);
          omega_X_r[2] = body_force->omega_rot*(body_force->axis_rot[0]*u[icv][1] - body_force->axis_rot[1]*u[icv][0]);

          double f_rhs[3];
          FOR_I3 f_rhs[i] = -rho[icv]*vol_cv[icv]*omega_X_r[i];

          FOR_I3 rhou_source[icv][i]  += f_rhs[i];
          FOR_I3 body_force->force[i] -= f_rhs[i];
          body_force->work            -= DOT_PRODUCT(f_rhs,u[icv]);
          body_force->volume          += vol_cv[icv];
          double m_cv[3] = CROSS_PRODUCT(x_cv[icv],f_rhs);
          FOR_I3 body_force->moment[i] -= m_cv[i];
        }
        else if (body_force->type == HOOK_TYPE) {
          assert(b_bodyForceHook);
        }
        else {
          assert(0);  //logic problem, shouldn't get here
        }
      }
    }
    if (b_bodyForceHook){
      computeBodyForceHook();
    }
  }
}

void HelmholtzSolver::initBodyForceHook(int &my_implicit_init_flag){
  CERR("User must override initBodyForceHook(int &my_implicit_init_flag) to use ADD_BODY_FORCE ... TYPE=HOOK");
}

void HelmholtzSolver::computeBodyForceHook(){
  CERR("User must override computeBodyForceHook() to use ADD_BODY_FORCE ... TYPE=HOOK");
}


void HelmholtzSolver::addFrameRotationSourceTerms(double (*rhs)[3]) { 

  if ( frame_rotation != NULL) {

    for (int icv = 0; icv < ncv; ++icv) {
      
      // term 1: -omega x omega x r
      // centrifugal force
      //
      double r[3];
      FOR_I3 r[i] = x_cv[icv][i] - frame_rotation[i+3];
      
      double coeff[3];
      coeff[0] = frame_rotation[1]*r[2] - frame_rotation[2]*r[1];
      coeff[1] = frame_rotation[2]*r[0] - frame_rotation[0]*r[2];
      coeff[2] = frame_rotation[0]*r[1] - frame_rotation[1]*r[0];
      
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
      
      FOR_I3 rhs[icv][i] += msource[i];
      
    }
  }

}

#undef FOR_BCZONE
