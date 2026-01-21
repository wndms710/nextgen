#ifndef IDEALGASSOLVERFLUX_HPP
#define IDEALGASSOLVERFLUX_HPP

#include "States.hpp"

//===========================
// flux definitions
//===========================
inline void calcInternalFluxNew(IdealGasRhs& flux, const FaNC& fa, const IdealGasState& cv0, 
                                const IdealGasState& cv1, double u_avg[3],
                                double& inv_v_avg, double& h_stag, double& p_avg) __attribute__((always_inline)); 
 
void calcInternalFluxNew(IdealGasRhs& flux, const FaNC& fa, const IdealGasState& cv0, 
                         const IdealGasState& cv1, double u_avg[3],
                         double& inv_v_avg, double& h_stag, double& p_avg) { 
  
  inv_v_avg =   2.0/(cv0.sp_vol + cv1.sp_vol); // specific volume 
  p_avg     =   0.5*(cv1.p + cv0.p); 

  double un_avg = 0.0;
  double u0u1   = 0.0;
  for (int i =0; i < 3; ++i) {
    double uu  = 0.5*(cv0.u[i] + cv1.u[i]); 
    un_avg      += uu * fa.n[i];
    u0u1        += cv0.u[i] * cv1.u[i]; 
    u_avg[i]     = uu;
    flux.rhou[i] = p_avg*fa.n[i];
  }
  
  double Frho =     un_avg*inv_v_avg;
  h_stag      =   cv1.h + cv0.h; 
  flux.rho    =   Frho; 
  
  for (int i =0; i < 3; ++i)   
    flux.rhou[i]     += Frho*u_avg[i];

  h_stag     = 0.5*(h_stag + u0u1);
  flux.rhoE  =  Frho*h_stag;
}

inline void calcInternalFluxSymmetricNew(IdealGasRhs& flux, const FaNC& fa, double u_avg[3], 
                                         double& inv_v_avg, double& h_stag, double& p_avg, 
                                         const IdealGasState& cv0, const IdealGasState& cv1)  
  __attribute__((always_inline));

void calcInternalFluxSymmetricNew(IdealGasRhs& flux, const FaNC& fa, double u_avg[3], 
                                  double& inv_v_avg, double& h_stag, double& p_avg, 
                                  const IdealGasState& cv0, const IdealGasState& cv1) { 

  // HACK: uncomment to skip symmetric flux...
  //flux.rho = 0.0;
  //flux.rhou[0] = 0.0;
  //flux.rhou[1] = 0.0;
  //flux.rhou[2] = 0.0;
  //flux.rhoE = 0.0;
  //return;

  // this returns the flux difference dotted on c... 
  // the flux that's populated is actually a flux difference.. 
    
  // computes the symmetric part of the flux; note that 
  // this assumes that the u_avg, inv_v_avg, etc were 
  // previously computed by the skew internal flux and are 
  // defined on entry...

  double du[3];
  double duc        = 0.0;
  double uavg_c     = 0.0;
  const double dp   = cv1.p - cv0.p;
  double dH         = cv1.h - cv0.h;
  const double drho = (cv0.sp_vol - cv1.sp_vol)/(cv0.sp_vol*cv1.sp_vol);

  for (int i = 0; i < 3; ++i) { 
    const double du_ = cv1.u[i] - cv0.u[i];
    du[i]   = du_;
    duc    += du_      * fa.c[i];
    uavg_c += u_avg[i] * fa.c[i];
    dH     += u_avg[i] * du_;

    // since we've loaded in the face c, lets start building the momentum flux difference
    flux.rhou[i] = dp * fa.c[i];
  }

  const double dFrho = uavg_c* drho + inv_v_avg * duc;
  const double Frho  = uavg_c* inv_v_avg;
  flux.rho           = dFrho;
  for (int i =0; i < 3; ++i) 
    flux.rhou[i]    += Frho* du[i] + u_avg[i] * dFrho;
  flux.rhoE          = h_stag*dFrho + Frho*dH;
}


inline void calcInternalCompactFluxNewGrad(IdealGasRhs& compact_flux,const FaCompact& fa, IdealGasState& cv0,
                                       IdealGasState& cv1, AuxGasProp& cv0_compact, 
                                       AuxGasProp& cv1_compact, const double u_avg[3], const double& inv_v_avg, 
                                       const double& h_stag,  
                                       const double grad_u0[3][3], const double grad_u1[3][3]) __attribute__((always_inline));

void calcInternalCompactFluxNewGrad(IdealGasRhs& compact_flux,const FaCompact& fa, IdealGasState& cv0,
                                    IdealGasState& cv1, AuxGasProp& cv0_compact, AuxGasProp& cv1_compact, 
                                    const double u_avg[3], const double& inv_v_avg, const double& h_stag, 
                                    const double grad_u0[3][3], const double grad_u1[3][3]) { 
 
  //HACK: uncomment to skip viscous flux...
  //return;

  const double aod_half               =   0.5*fa.area_over_delta;
  const double delta                  =   fa.area/fa.area_over_delta;
  const double mu_tmp                 =   0.5*(cv0_compact.mu_total + cv1_compact.mu_total)*fa.area;
  const double mu_total_c             =   (cv0_compact.mu_total+cv1_compact.mu_total)*aod_half;
  const double loc_coeff              =   (cv0_compact.loc_total + cv1_compact.loc_total)*aod_half;
  const double dhalf_usq              =   cv1_compact.half_usq - cv0_compact.half_usq;
  const double dh                     =   cv1.h - cv0.h;

  // dun, dun2 is actually on the compact set of faces ..
  double u0_cn = 0.0;
  double u1_cn = 0.0;

  for (int i = 0; i < 3; ++i) { 
    u0_cn += cv0.u[i]* fa.unit_n[i];
    u1_cn += cv1.u[i]* fa.unit_n[i];
  }

  const double dun           =  u1_cn - u0_cn;
  const double one_third_dun =  dun/3.0;
  const double one_sixth_dun2 = (u1_cn*u1_cn - u0_cn*u0_cn)/6.0;

  /// entropy analysis based stabilization contributions..
  
  const double dp       = cv1.p      - cv0.p;
  const double fgr_avg  = 0.5*(cv1_compact.fgr + cv0_compact.fgr);
 
  const double rho_avg  = 0.5*(cv1_compact.rho + cv0_compact.rho);
  const double un_avg   = 0.5*(u0_cn + u1_cn);
  const double sos_avg  = 0.5*(cv0_compact.sos + cv1_compact.sos);
  const double Ma_n     = abs(un_avg)/sos_avg;
  const double drho     = cv1_compact.rho - cv0_compact.rho;

  const double lf_coeff = min(fgr_avg*Ma_n,0.5)*(sos_avg + abs(un_avg))*fa.area;
  //const double lf_coeff = fgr_avg*abs(un_avg)*fa.area; 
  const double lfdrho   = lf_coeff*drho;
  const double lfrhoa   = lf_coeff*rho_avg;

  double Skk_avg = 0.0;
  for (int i = 0; i < 3; ++i) 
    Skk_avg += 0.5*(grad_u0[i][i] + grad_u1[i][i]);

  /// acoustic sgs modeling... note that a_sgs is a_sgs over c^2 from the previous formulation...
  // note that a_sgs has units of cv0_compact has units of velocity/delta

  double a_sgs_coeff       = 0.5*(cv0_compact.a_sgs + cv1_compact.a_sgs)*delta*fa.area;
  a_sgs_coeff             += (cv0_compact.symp_fax + cv1_compact.symp_fax)*aod_half;
  const double a_sgs_cdp   = a_sgs_coeff*dp;

  compact_flux.rho -= lfdrho + a_sgs_cdp;

  for (int i =0; i < 3; ++i) {  
    compact_flux.rhou[i] -= mu_total_c*(cv1.u[i] - cv0.u[i] + one_third_dun*fa.unit_n[i]);
    
    compact_flux.rhou[i] -= lfdrho*u_avg[i];
    compact_flux.rhou[i] -= lfrhoa*dun*fa.unit_n[i];
    //compact_flux.rhou[i] -= lf_coeff*( u_avg[i]*drho + rho_avg*(u1_cn - u0_cn)*fa.unit_n[i]);
    
    compact_flux.rhou[i] -= u_avg[i]*a_sgs_cdp;
  }

  compact_flux.rhoE -= mu_total_c*(dhalf_usq + one_sixth_dun2);

  // split the enthalpy based stabilization out .. so we can use the local un to scale 
  // the eigenvalues... 

  double h_fax = (cv0_compact.symh_fax + cv1_compact.symh_fax)*rho_avg*abs(un_avg)*aod_half;
  h_fax        = max(h_fax - loc_coeff, 0.0);
  compact_flux.rhoE -= (loc_coeff + h_fax)*dh; 
  
  // total energy contribution of the asgs .. 
  
  compact_flux.rhoE -= a_sgs_cdp*h_stag;

  const double pv0   = cv0.p*cv0.sp_vol;
  const double pv1   = cv1.p*cv1.sp_vol;
  const double pv_avg = 0.5*(pv0 + pv1);
  const double de    = dh - (pv1 - pv0); 

  compact_flux.rhoE -= lfdrho*(h_stag - pv_avg);
  compact_flux.rhoE -= lfrhoa*de;
  compact_flux.rhoE -= lfrhoa*un_avg*dun;

  //compact_flux.rhoE -= lf_coeff*e_avg*drho + lf_coeff*(DOT_PRODUCT(u_avg,u_avg) - 0.5*(cv1_compact.half_usq + cv0_compact.half_usq))*drho;
  //compact_flux.rhoE -= lf_coeff*rho_avg*de + lf_coeff*rho_avg*un_avg*(u1_cn - u0_cn);

  // viscous transpose terms... 

  double dundx[3]           = {0.0, 0.0, 0.0};
  double one_third_dundxn   = 0.0;
  double two_third_Skk      = 0.0;
  for (int k = 0; k < 3; ++k) { 
    dundx[k] = 0.5* ( (grad_u0[0][k] + grad_u1[0][k])*fa.unit_n[0] +  
                      (grad_u0[1][k] + grad_u1[1][k])*fa.unit_n[1] +
                      (grad_u0[2][k] + grad_u1[2][k])*fa.unit_n[2] ); 

    one_third_dundxn  += dundx[k]*fa.unit_n[k];
    two_third_Skk +=   grad_u0[k][k] + grad_u1[k][k];
  }

  two_third_Skk /= 3.0;
  one_third_dundxn /= 3.0;

  for (int i =0; i < 3; ++i) {  
    const double tmp      = (dundx[i] - fa.unit_n[i]*(one_third_dundxn + two_third_Skk))*mu_tmp;
    compact_flux.rhou[i] -= tmp;  // mu_tmp = mu_avg*area... tmp reused here
    compact_flux.rhoE    -= u_avg[i]*tmp;
  }

}

// this is the variable gamma version of the flux.. 
inline void calcRiemannFlux(IdealGasRhs& flux, const double n_bf[3], 
                            const IdealGasState& cv0, const IdealGasState& cv1,
                            const double gammaL, const double gammaR) { 

  // we should consider precomputing the unit normal vector for the boundary
  const double area = MAG(n_bf);
  double unit_n[3];
  for (int i =0; i < 3; ++i) 
    unit_n[i] = n_bf[i]/area;

  // for efficiency, we could pass in the properties object as well but 
  // that would complicate the calling procedure.  we are going to recompute
  // usq and rho (from sp_vol) in the following.  if the bcs that use this
  // become a large fraction of the cost, then we can revisit..
  double unL  = DOT_PRODUCT(cv0.u, unit_n);
  double uLuL = DOT_PRODUCT(cv0.u, cv0.u); 
  double cL   = sqrt(gammaL*cv0.p*cv0.sp_vol); 
  double pL   = cv0.p;
  double HL   = cv0.h + 0.5*uLuL; 
  double rhoL = 1.0/cv0.sp_vol; 
  double EL   = rhoL*HL - cv0.p;
  
  double unR  = DOT_PRODUCT(cv1.u, unit_n);
  double uRuR = DOT_PRODUCT(cv1.u, cv1.u); 
  double cR   = sqrt(gammaR*cv1.p*cv1.sp_vol); 
  double pR   = cv1.p;
  double HR   = cv1.h + 0.5*uRuR; 
  double rhoR = 1.0/cv1.sp_vol; 
  double ER   = rhoR*HR - cv1.p;
  
  // Roe's aveaging
  double Rrho = sqrt(rhoR/rhoL);
  double tmp = 1.0/(1.0+Rrho);
  double velRoe[3];
  FOR_I3 velRoe[i] = tmp*(cv0.u[i] + cv1.u[i]*Rrho);
  double uRoe  = DOT_PRODUCT(velRoe, unit_n);
  
  double gamPdivRho = tmp*((gammaL*pL/rhoL+0.5*(gammaL-1.0)*uLuL) + (gammaR*pR/rhoR+0.5*(gammaR-1.0)*uRuR)*Rrho);
  double cRoe  = sqrt(gamPdivRho - ((gammaL+gammaR)*0.5-1.0)*0.5*DOT_PRODUCT(velRoe, velRoe));
  
  // speed of sound at L and R
  double sL = min(uRoe-cRoe, unL-cL);
  double sR = max(uRoe+cRoe, unR+cR);
  
  // speed of contact surface
  double sM = (pL-pR-rhoL*unL*(sL-unL)+rhoR*unR*(sR-unR))/(rhoR*(sR-unR)-rhoL*(sL-unL));
  
  // pressure at right and left (pR=pL) side of contact surface
  double pStar = rhoR*(unR-sR)*(unR-sM)+pR;
  
  if (sM >= 0.0) {
    if (sL > 0.0) {
      flux.rho = rhoL*unL;
      FOR_I3 flux.rhou[i] = rhoL*cv0.u[i]*unL + pL*unit_n[i];
      flux.rhoE = (EL+pL)*unL;
    } else {
      double invSLmSs = 1.0/(sL-sM);
      double sLmuL = sL-unL;
      double rhoSL = rhoL*sLmuL*invSLmSs;
      double rhouSL[3];
      FOR_I3 rhouSL[i] = (rhoL*cv0.u[i]*sLmuL+(pStar-pL)*unit_n[i])*invSLmSs;
      double eSL = (sLmuL*EL-pL*unL+pStar*sM)*invSLmSs;
      
      flux.rho = rhoSL*sM;
      FOR_I3 flux.rhou[i] = rhouSL[i]*sM + pStar*unit_n[i];
      flux.rhoE = (eSL+pStar)*sM;
    }
  } else {
    if (sR >= 0.0) {
      double invSRmSs = 1.0/(sR-sM);
      double sRmuR = sR-unR;
      double rhoSR = rhoR*sRmuR*invSRmSs;
      double rhouSR[3];
      FOR_I3 rhouSR[i] = (rhoR*cv1.u[i]*sRmuR+(pStar-pR)*unit_n[i])*invSRmSs;
      double eSR = (sRmuR*ER-pR*unR+pStar*sM)*invSRmSs;
      
      flux.rho = rhoSR*sM;
      FOR_I3 flux.rhou[i] = rhouSR[i]*sM + pStar*unit_n[i];
      flux.rhoE = (eSR+pStar)*sM;
    } else {
      flux.rho = rhoR*unR;
      FOR_I3 flux.rhou[i] = rhoR*cv1.u[i]*unR + pR*unit_n[i];
      flux.rhoE = (ER+pR)*unR;
    }
  }
  flux.rho            *= area;
  FOR_I3 flux.rhou[i] *= area;
  flux.rhoE           *= area;
}
#endif
