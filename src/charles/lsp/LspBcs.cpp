#include "LspBcs.hpp"
#include "IdealGasSolver.hpp"
#include "IdealGasSolverWithLsp.hpp"
lpBaseBc::~lpBaseBc() {COUT1("~lpBaseBc()");}

void lpBaseBc::queryBc() { 
  int total_np_hit;
  MPI_Reduce(&np_hit,&total_np_hit,1,MPI_INT,MPI_SUM,0,mpi_comm);
  COUT1(" > lsp query, " << zone_ptr->getName() << " np: " << total_np_hit);
}

void lpOutletBc::applyBc(LspState& lp, double* xi, double* st_normal, double* buf) {
  lp.icv = -1;        // set to negative to not consider this lp anymore
  lp.flag = TRASH;    // set the flag to be removed in the next recycle
  np_hit++;
}

bool lpOutletBc::applyBc2(LspState& lp, double* dx, double* st_normal) {
  np_hit++;
  return false; // This return indicates that the particle should be removed
}

void lpPerfectBounceBc::applyBc(LspState& lp, double* xi, double* st_normal, double* buf) {
  const double dx[3]  = DIFF(lp.xp,xi);
  const double dp     = DOT_PRODUCT(st_normal,dx); //assert(dp > 0);
  const double nmag2  = DOT_PRODUCT(st_normal,st_normal);
  FOR_I3 lp.xp[i]    -= 2.0*fabs(dp)*st_normal[i]/nmag2;
  // and update xp0 to the xi_closest, because the remaining part of the 
  // trajectory may have additional intersections...
  FOR_I3 lp.xp0[i]    =  xi[i];
  const double un     = DOT_PRODUCT(st_normal,lp.up);
  FOR_I3 lp.up[i]    -= 2.0*un*st_normal[i]/nmag2;
  
  np_hit++;
}

bool lpPerfectBounceBc::applyBc2(LspState& lp, double* dx, double* st_normal) {
  // dx is the vector from the particle to the intersection point on the relevant tri. 
  // st_normal is the outward normal of the relevant surface tri with area magnitude...
  const double dp             = DOT_PRODUCT(dx,st_normal); // should be negative
  assert(dp <= 0.0);
  const double st_normal_mag2 = DOT_PRODUCT(st_normal,st_normal);
  assert(st_normal_mag2 > 0.0);
  FOR_I3 lp.xp[i]            += 2.0*dp*st_normal[i]/st_normal_mag2;
  const double un             = DOT_PRODUCT(lp.up,st_normal);
  FOR_I3 lp.up[i]            -= 2.0*un*st_normal[i]/st_normal_mag2;

  np_hit++;
  return true;
}

void lpIrregularBounceBc::initialHook() { 
  // Read the parameter values from the input file
  FOR_PARAM_MATCHING("LSP.BC") {
    if (param->getString(0) == zone_ptr->getName()) {
      int iarg = 2;
      while (iarg < param->size()) { 
        const string token = param->getString(iarg++);
        if      (token == "RESTITUTION")          resCoef             = param->getDouble(iarg++);
        else if (token == "MU_S")                 mu_s                = param->getDouble(iarg++);
        else if (token == "MU_K")                 mu_k                = param->getDouble(iarg++);
        else if (token == "ROUGHNESS_MAX_ANGLE")  roughness_max_angle = param->getDouble(iarg++);
        else CERR("Error: unrecognized parameter in LSP.BC for zone " << zone_ptr->getName() << " bad token: " << token);
      }
      break;
    }
  }

  // Output the paramters even if they weren't properly initialized
  COUT1(" > Assigning IRREGULAR_BOUNCE lsp bc to zone " << zone_ptr->getName()  << "\n" 
        " > RESTITUTION = "                             << resCoef              << "\n" 
        " > MU_S = "                                    << mu_s                 << "\n" 
        " > MU_K = "                                    << mu_k                 << "\n" 
        " > ROUGHNESS_MAX_ANGLE = "                     << roughness_max_angle);

  // Error if the parameters were improperly assigned...
  if ((resCoef < 0.0) || (mu_s < 0.0) || (mu_k < 0.0) || (roughness_max_angle < 0.0))
    CERR("ERROR: LSP.BC IRREGULAR BOUNCE does not have valid required parameters! (Should be nonnegative. See above\n");
}

void lpIrregularBounceBc::applyBc(LspState& lp, double* xi, double* st_normal, double* buf) {
  const double dx[3] = DIFF(lp.xp,xi);
  const double dp = DOT_PRODUCT(st_normal,dx); //assert(dp > 0);
  const double nmag2 = DOT_PRODUCT(st_normal,st_normal);
  FOR_I3 lp.xp[i] -= 2.0*fabs(dp)*st_normal[i]/nmag2;
  // and update xp0 to the xi_closest, because the remaining part of the 
  // trajectory may have additional intersections...
  FOR_I3 lp.xp0[i] = xi[i];

  // form e1 and e2
  double n[3]; FOR_I3 n[i] = st_normal[i];
  NORMALIZE(n);
  double e1[3];
  double e2[3];
  if ((fabs(n[0]) <= fabs(n[1])) and (fabs(n[0]) <= fabs(n[2]))) {
    e1[0] = 1; e1[1] = 0; e1[2] = 0;
  } else if ((fabs(n[1]) <= fabs(n[0])) and (fabs(n[1]) <= fabs(n[2]))) {
    e1[0] = 0; e1[1] = 1; e1[2] = 0;
  } else if ((fabs(n[2]) <= fabs(n[0])) and (fabs(n[2]) <= fabs(n[1]))) {
    e1[0] = 0; e1[1] = 0; e1[2] = 1;
  }
  double temp1[3] = CROSS_PRODUCT(n,e1);
  FOR_J3 e2[j] = temp1[j];
  NORMALIZE(e2);
  double temp2[3] = CROSS_PRODUCT(e2,n);
  FOR_J3 e1[j] = temp2[j];
  //assert(abs(DOT_PRODUCT(e1,e2))<1e-10);
  //assert(abs(DOT_PRODUCT(e1,n))<1e-10);
  //assert(abs(DOT_PRODUCT(e2,n))<1e-10);
  //assert((DOT_PRODUCT(e1,e1)-1.)<1e-10);
  //assert((DOT_PRODUCT(e2,e2)-1.)<1e-10);
  //assert((DOT_PRODUCT(n,n)-1.)<1e-10);

  double rand_e1 = double(rand())/double(RAND_MAX); // between 0 and 1
  double rand_e2 = sqrt(1. - rand_e1*rand_e1);
  double dir[3];
  FOR_I3 dir[i] = rand_e1*e1[i] + rand_e2*e2[i];
  //assert((DOT_PRODUCT(dir,dir)-1.)<1e-10);
  //assert(abs(DOT_PRODUCT(dir,n))<1e-10);

  // perturb the normal vector first
  //double u_orth[3];     // the orthagonal to the normal part of the velocity vector
  //FOR_I3 u_orth[i] = lp.up[i] - DOT_PRODUCT(st_normal,lp.up)*st_normal[i]/nmag2;
  //double mag_u_orth = MAG(u_orth);
  double random = (double(rand())/double(RAND_MAX) - 0.5)*2.0;
  double MAX_GAMMA = roughness_max_angle/180.0*M_PI; 
  double tan_Gamma = tan(MAX_GAMMA * random);
  //double st_normal_pert[3];     // perturbed normal vector, tan_Gamma = |dn|/|n|
  //FOR_I3 st_normal_pert[i] = st_normal[i] + sqrt(nmag2)*tan_Gamma*u_orth[i]/mag_u_orth;
  double normal_pert[3];     // perturbed normal vector, tan_Gamma = |dn|/|n|
  FOR_I3 normal_pert[i] = n[i] + tan_Gamma*dir[i];
  const double nmag2_pert = DOT_PRODUCT(normal_pert,normal_pert);
  const double un = DOT_PRODUCT(normal_pert,lp.up);
  double un_old[3];
  double ut_old[3];
  FOR_I3 un_old[i] = normal_pert[i]*un/nmag2_pert;
  FOR_I3 ut_old[i] = lp.up[i] - un_old[i];
  const double mag_un_old = MAG(un_old);
  const double mag_ut_old = MAG(ut_old);
  double un_new[3]; 
  double ut_new[3];
  FOR_I3 un_new[i] = -un_old[i]*resCoef;
  double slide = mag_ut_old - 7./2.*mu_s*(1.+resCoef)*mag_un_old; // |u_p| < 7/2 mu_s (1+e) |v_p| -> no sliding
  //FOR_I3 ut_new[i] = ut_old[i]*(1. - mu_k*(1+resCoef)*mag_un_old/mag_ut_old); // sliding
  //FOR_I3 ut_new[i] = 5./7.*ut_old[i]; // no sliding
  if (slide < 0) { // no sliding
    //cout << "impact_mode: "<< "NO slide" << " , 7./2.*mu_s*(1.+resCoef): " << 7./2.*mu_s*(1.+resCoef) << " , u/v: " << mag_ut_old/mag_un_old << " , tan_Gamma: " << tan_Gamma << endl;
    FOR_I3 ut_new[i] = 5./7.*ut_old[i];
  } else { // sliding
    //cout << "impact_mode: "<< "slide" << " , 7./2.*mu_s*(1.+resCoef): " << 7./2.*mu_s*(1.+resCoef) << " , u/v: " << mag_ut_old/mag_un_old << " , tan_Gamma: " << tan_Gamma << endl;
    FOR_I3 ut_new[i] = ut_old[i]*(1. - mu_k*(1+resCoef)*mag_un_old/mag_ut_old);
  }
  FOR_I3 lp.up[i] = ut_new[i] + un_new[i];

  np_hit++;
}

bool lpIrregularBounceBc::applyBc2(LspState& lp, double* dx, double* st_normal) {
  // dx is the vector from the particle to the intersection point on the relevant tri. 
  // st_normal is the outward normal of the relevant surface tri with area magnitude...
  const double dp                 = DOT_PRODUCT(dx,st_normal);        assert(dp <= 0.0);            // (should be negative)
  const double inv_st_normal_mag2 = 1.0/DOT_PRODUCT(st_normal,st_normal); assert(inv_st_normal_mag2 > 0.0); // (should be positive)
  FOR_I3 lp.xp[i]                += 2.0*dp*st_normal[i]*inv_st_normal_mag2;
  //const double un                 = DOT_PRODUCT(lp.up,st_normal);
  //FOR_I3 lp.up[i]                -= 2.0*un*st_normal[i]*inv_st_normal_mag2;
  
  // Form e1 and e2
  double n[3], e2[3],  e1[3] = {0.0, 0.0, 0.0};
  FOR_I3 n[i] = st_normal[i]*inv_st_normal_mag2;
  if      ((fabs(n[0]) <= fabs(n[1])) and (fabs(n[0]) <= fabs(n[2]))) e1[0] = 1.0;  
  else if ((fabs(n[1]) <= fabs(n[0])) and (fabs(n[1]) <= fabs(n[2]))) e1[1] = 1.0;   
  else if ((fabs(n[2]) <= fabs(n[0])) and (fabs(n[2]) <= fabs(n[1]))) e1[2] = 1.0;
  else assert(0);
  double temp1[3]                 = CROSS_PRODUCT(n,e1);
  FOR_J3 e2[j]                    = temp1[j];
  NORMALIZE(e2);
  double temp2[3]                 = CROSS_PRODUCT(e2,n);
  FOR_J3 e1[j]                    = temp2[j];

  double rand_e1                  = double (rand())/double(RAND_MAX); // between 0 and 1
  double rand_e2                  = sqrt(1.0 - rand_e1*rand_e1);
  double dir[3];
  FOR_I3 dir[i]                   = rand_e1*e1[i] + rand_e2*e2[i];

  double random                   = (double(rand())/double(RAND_MAX) - 0.5)*2.0;
  double MAX_GAMMA                = roughness_max_angle/180.0*M_PI;
  double tan_Gamma                = tan(MAX_GAMMA*random);

  double normal_pert[3]; // perturbed normal vector, tan_gamma = |dn|/|n|
  FOR_I3 normal_pert[i]           = n[i] + tan_Gamma*dir[i];
  const double inv_nmag2_pert     = 1.0/(DOT_PRODUCT(normal_pert, normal_pert));
  const double un                 = DOT_PRODUCT(normal_pert, lp.up);
  double un_old[3], ut_old[3];
  FOR_I3 un_old[i]                = normal_pert[i]*un*inv_nmag2_pert;
  FOR_I3 ut_old[i]                = lp.up[i] - un_old[i];
  const double mag_un_old         = MAG(un_old);
  const double mag_ut_old         = MAG(ut_old);
  double un_new[3], ut_new[3];
  FOR_I3 un_new[i]                = -un_old[i]*resCoef;
  // Sliding criterion: |u_p| < 7/2 mu_s (1 + e) |v_p| ->no_sliding
  if (mag_ut_old - 3.5*mu_s*(1.0+resCoef)*mag_un_old < 0.0) // No sliding
    FOR_I3 
      ut_new[i]                   = 5.0/7.0*ut_old[i];
  else // Sliding
    FOR_I3 
      ut_new[i] = ut_old[i]*(1.0 - mu_k*(1.0 + resCoef)*mag_un_old/mag_ut_old);

  FOR_I3 lp.up[i] = ut_new[i] + un_new[i];
  
  np_hit++;
  return true;
}

void lpInelasticBounceBc::initialHook() { 
  // Initialize stats here (no restart/registration) to set CoR values...
  updateStats();

  // Read the parameter values from the input file
  FOR_PARAM_MATCHING("LSP.BC") { 
    if (param->getString(0) == zone_ptr->getName()) { 
      assert(param->getString(1) == "INELASTIC_BOUNCE");
      assert(param->getString(2) == "MU");
      mu = param->getDouble(3);
      break;
    }
  }

  COUT1(" > Assigning INELASTIC_BOUNCE lsp bc to zone " << zone_ptr->getName()  << "\n"
        " > MU = "                                      << mu                   << "\n")

  // Error if the parameters were unassigned or improperly set
  if (mu < 0.0) CERR(" > LSP.BC " << zone_ptr->getName() << " (IRREGULAR BOUNCE) does not have required parameter MU!");

  // HACK - solid set to dust in IdealGasSolverWithLsp.hpp because we do not have solver templating
  if (this->solid == NULL) {CERR("ERROR: INELASTIC BOUNCE BC (assigned to zone " << zone_ptr->getName() 
                                  << " only works if the LSP.MATERIAL is dust (for now)!!!!")}

  COUT1(" > YIELD STRESS = "        << this->solid->calcYieldStress()   << "\n" 
        " > ELASTICITY = "          << this->solid->calcElasticity());

  if ((solid->calcYieldStress() < 0) || (solid->calcElasticity() < 0))
    CERR("Error: for INELASTIC_BOUNCE BC you must specify material yield stress and elasticity."
          << "Specify <material-name>.YIELD_STRESS, <material-name>.ELASTICITY in the input file.");
}

void lpInelasticBounceBc::applyBc(LspState& lp, double* xi, double* st_normal, double* buf) {
  const double rhop = solid->calcRho(lp.Tp);
  const double yieldStress = solid->calcYieldStress();
  const double elasticity = solid->calcElasticity();
  assert(yieldStress > 0);
  assert(elasticity > 0);

  const double nmag = MAG(st_normal);
  const double un = DOT_PRODUCT(st_normal,lp.up)/nmag;
  double up_tan[3]; FOR_I3 up_tan[i] = lp.up[i] - un*st_normal[i]/nmag;
  const double ut_mag = MAG(up_tan);
  const double E_crit = yieldStress*yieldStress / elasticity / 2.;
  const double KEn = rhop/2. * un*un;

  double un2 = un;
  double CoR = 1.;
  if (KEn > E_crit) {
    un2 = yieldStress/sqrt(rhop*elasticity);
    CoR = un2/un;
  }

  const double umag = MAG(lp.up);
  const double cosine = ((umag>0) ? (ut_mag/umag):1.);
  const double ut2 = ut_mag - mu*fabs(un)*(1. + CoR)*cosine*cosine;

  FOR_I3 lp.up[i] = un2*(-st_normal[i])/nmag + ut2*up_tan[i]/ut_mag;

  const double dx[3] = DIFF(lp.xp,xi);
  const double dp = DOT_PRODUCT(st_normal,dx); //assert(dp > 0);
  //const double nmag2 = DOT_PRODUCT(st_normal,st_normal);
  FOR_I3 lp.xp[i] -= 2.0*fabs(dp)*st_normal[i]/(nmag*nmag);
  // and update xp0 to the xi_closest, because the remaining part of the 
  // trajectory may have additional intersections...
  FOR_I3 lp.xp0[i] = xi[i];
  
  np_hit++;
}

bool lpInelasticBounceBc::applyBc2(LspState& lp, double* dx, double* st_normal) {
  // dx is the vector from the particle to the intersection point on the relevant tri. st_normal is 
  // the outward normal of the relevant surface tri with area magnitude...
  const double dp             = DOT_PRODUCT(dx,st_normal);        assert(dp <= 0.0);            // (should be negative)
  const double st_normal_mag2 = DOT_PRODUCT(st_normal,st_normal); assert(st_normal_mag2 > 0.0); // (should be positive)
  FOR_I3 lp.xp[i]            += 2.0*dp*st_normal[i]/st_normal_mag2;

  const double rhop           = solid->calcRho(lp.Tp);
  const double yieldStress    = solid->calcYieldStress(); assert(yieldStress > 0);
  const double elasticity     = solid->calcElasticity();  assert(elasticity > 0);

  const double inv_nmag       = 1.0/MAG(st_normal);
  const double un             = DOT_PRODUCT(st_normal,lp.up)*inv_nmag;
  double up_tan[3]; 
  FOR_I3 up_tan[i]            = lp.up[i] - un*st_normal[i]*inv_nmag;
  const double ut_mag         = MAG(up_tan);
  const double E_crit         = 0.5*yieldStress*yieldStress / elasticity;
  const double KEn            = 0.5*rhop*un*un;

  double un2                  = un;
  double CoR                  = 1.0;
  if (KEn > E_crit) {
    un2                       = yieldStress/sqrt(rhop*elasticity);
    CoR                       = un2/un;
  }

  // XXX Why don't we just regularize ut2 instead of this linearization of cosine stuff?
  const double umag           = MAG(lp.up);
  const double cosine         = ((umag>0) ? (ut_mag/umag):1.);
  const double ut2            = ut_mag - mu*fabs(un)*(1. + CoR)*cosine*cosine;

  FOR_I3 lp.up[i]             = un2*(-st_normal[i])*inv_nmag + ut2*up_tan[i]/ut_mag;

  // Update stats for query
  CoR_min   = min(CoR_min, CoR);
  CoR_max   = max(CoR_max, CoR);
  CoR_bar  += CoR;
  np_hit++;
  return true;
}

void lpInelasticBounceBc::queryBc() { 
  int total_np_hit;
  double CoR_min_glob, CoR_max_glob, CoR_bar_glob;
  MPI_Reduce(&np_hit,   &total_np_hit, 1, MPI_INT,    MPI_SUM, 0, mpi_comm);
  MPI_Reduce(&CoR_min,  &CoR_min_glob, 1, MPI_DOUBLE, MPI_MIN, 0, mpi_comm);
  MPI_Reduce(&CoR_max,  &CoR_max_glob, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
  MPI_Reduce(&CoR_bar,  &CoR_bar_glob, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
  COUT1(" > lsp query, " << zone_ptr->getName()
     << " np: "          << total_np_hit
     << " CoR_min: "     << CoR_min_glob
     << " CoR_bar: "     << CoR_bar_glob/(double) total_np_hit
     << " CoR_max: "     << CoR_max_glob);
}

void lpInelasticBounceBc::updateStats() { 
  // Reset stats
  //if (solver->step % solver->lsp_stats_interval == 0) {
  if (1 == 1) { 
    np_hit  = 0;
    CoR_min = HUGE_VAL;
    CoR_max = -HUGE_VAL;
    CoR_bar = 0.0;
  }
}
