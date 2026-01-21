#ifndef BCCOMMON_HPP
#define BCCOMMON_HPP

#include "CtiRegister.hpp"

template<typename SolverT, typename RhsT>
class BcTemplate : public CtiRegister::CtiDataProducer {
public:

  BfZone* zone_ptr;
  SolverT* solver;

  // for passive scalar transport -- may be necessary to record
  // the boundary condition mass flux as it is computed. if you
  // want to specify more than a single value for the bc (at
  // an inlet for instance), you can do so but you will need to
  // write a special boundary condition.

  double* mf;
  double * scbc_vals;
  BcTemplate(BfZone* _zone_ptr, SolverT* _solver) {

    // in order facilitate debugging and verification of the
    // particular boundary condition types, these pointers must
    // be allowed to take a null state.

    zone_ptr = _zone_ptr;
    solver   = _solver;

    mf        = NULL;
    scbc_vals = NULL;

    addCtiDataProducer(this);

  }

  virtual ~BcTemplate() {
    DELETE(mf);
    DELETE(scbc_vals);

  }

  string getName() const {
    assert( zone_ptr) ;
    return zone_ptr->getName();
  }

  // pure virtual methods that bcs must provide
  // the addBoundaryFlux routine has access to the solver rather
  // than a particular data structure bc depending on the particular bc
  // it may be more convenient to access the packed state or the
  // continguous data (for instance, a slip wall only requires p..)

  virtual void addBoundaryFlux(RhsT* rhs) const = 0;

  // advance any boundary condition degrees of freedom; the time advance provided
  // here is an explicit rk3 to mimic the flow solver

  virtual void calcRhs(const double time, const int rk_stage) = 0;
  virtual void rk3Step(const double dt, const int rk_stage) = 0;

  // registration of boundary condition variables..
  // ACHTUNG--------------
  // in addition, the constructors below are all empty, the registration
  // hook provides a common function where variables are parsed and read
  // ---------------------

  virtual void initData() = 0;
  virtual void initialHook() = 0; // perhaps doesnt need to be pure virtual?

  // a boundary condition may be asked to store or restore its data during 
  // a run...
  
  virtual int storeData(double * buf) const = 0; // if you pass NULL to this it will simply return a count of doubles
  virtual int restoreData(double * buf) = 0;

  // each boundary condition can also define a query which will report
  // relevant information about the bc (dumpRange, integrated values) depending
  // on what is relevant for its particular case.  the default behavior is
  // that the query is empty..

  virtual void query(const string& param_str) const {}

  // force routines. force_bf must be specified by each boundary condition.

  virtual void force_bf(double (*f_bf)[9]) const = 0;

  void force(double (*f_buf)[3],double (*m_buf)[3],const double p_inf = 0.0) const {

    //f_buf: pressure, viscous, convective force components

    FOR_J3 FOR_I3 f_buf[i][j] = 0.0;

    //m_buf: compute moment about origin

    FOR_J3 FOR_I3 m_buf[i][j] = 0.0;

    double (*f_bf)[9] = new double[zone_ptr->nbf][9];
    force_bf(f_bf);

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      // put in dp form...
      FOR_I3 f_bf[ibf][3+i] -= p_inf*zone_ptr->n_bf[ibf][i];

      FOR_I3 f_buf[0][i] += f_bf[ibf][i]; 
      FOR_I3 f_buf[1][i] += f_bf[ibf][3+i];
      FOR_I3 f_buf[2][i] += f_bf[ibf][6+i];

      FOR_I3 {
        const double m_bf[3] = CROSS_PRODUCT(zone_ptr->x_bf[ibf],&f_bf[ibf][3*i]);
        FOR_J3 m_buf[i][j] += m_bf[j];
      }

    }

    delete[] f_bf;
  }

  // the default behavior is that the Bc cannot produce any data

  virtual CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
      const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

    if ( zone_ptr->isLeadingStrMatched(name)) {

      const string proj_str = zone_ptr->getName() + ":" + "proj";
      //const string xbf_str = zone_ptr->getName() + ":" + "x_bf"; // should work now w/o this - CI

      if ( name == proj_str) {

        if ( args.size() != 1)
          return CtiRegister::CTI_DATA_ARG_COUNT;

        list<CtiRegister::CtiData>::iterator arg = args.begin();
        if ( arg->getTopology() != CV_DATA)
          return CtiRegister::CTI_DATA_NOT_VALID;

        if ( arg->getType() == DN_DATA) {

          double * tmp = zone_ptr->createBfD1Data(v);

          if ( b_eval_func) {

            for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

              const int icv = zone_ptr->cvobf[ibf];
              tmp[ibf]      = arg->dn(icv);

            }
          }

          return CtiRegister::CTI_DATA_OK;

        } else if ( arg->getType() == DN3_DATA) {

          double (*tmp)[3] = zone_ptr->createBfD2Data(v);

          if ( b_eval_func) {

            for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

              const int icv = zone_ptr->cvobf[ibf];
              for (int i = 0; i < 3; ++i)
                tmp[ibf][i] = arg->dn3(icv,i);

            }
          }

          return CtiRegister::CTI_DATA_OK;

        } else {

          return CtiRegister::CTI_DATA_NOT_VALID;

        }

      }
      /*
      else if ( name == xbf_str) {

        if ( args.size() != 0)
          return CtiRegister::CTI_DATA_ARG_COUNT;

        double (*tmp)[3] = zone_ptr->createBfD2Data(v);

        if ( b_eval_func) {

          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            for (int i = 0; i < 3; ++i)
              tmp[ibf][i] = zone_ptr->x_bf[ibf][i];

          }
        }

        return CtiRegister::CTI_DATA_OK;

      }
      */
    }

    return CtiRegister::CTI_DATA_NOT_FOUND;

  }

  //-------------------------------------
  // passive scalar transport
  //-------------------------------------

  virtual void addScalarBoundaryFlux(double * rhs_sc) const {}

  void parseScalarBc(Param * param) {

    assert ( solver != NULL);

    if ( solver->nsc_transport == 0 )
      return;

    int * sc_flag = new int[solver->nsc_transport];
    for (int isc = 0; isc < solver->nsc_transport; ++isc)
      sc_flag[isc] = 0;

    assert ( scbc_vals == NULL);
    scbc_vals = new double[solver->nsc_transport];

    int iarg = 0;
    while ( iarg < param->size() ) {

      const string token = param->getString(iarg++);

      const int isc  = solver->getScalarIndex(token);
      if ( isc >= 0) {

        if ( sc_flag[isc] == 1) {

          CERR(" > multiply specified scalar bcs for var: " << token);

        } else {

          sc_flag[isc]   = 1;
          scbc_vals[isc] = param->getDouble(iarg++);

        }
      }
    }

    // check to make sure that all the scalars were set.

    int err_flag = 0;
    for (int isc = 0; isc < solver->nsc_transport; ++isc) {

      if ( sc_flag[isc] == 0) {

        if ( mpi_rank == 0 )
          cout << " > failed to provide scalar bc for : " << solver->getScalarName(isc)
               << "   at zone : " << zone_ptr->getName() << endl;

        err_flag = 1;

      } else {

        assert( sc_flag[isc] == 1);

      }
    }

    delete[] sc_flag;

    if ( err_flag == 1) { 
      CERR( " > please correct boundary condition syntax for scalar vars above ");
    }

  }

  void _addScalarBoundaryFlux(double * rhs_sc) const {

    assert ( solver != NULL);
    assert ( mf != NULL);
    assert ( scbc_vals != NULL);

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv = zone_ptr->cvobf[ibf];
      for (int isc = 0; isc < solver->nsc_transport; ++isc) {
        if ( mf[ibf] < 0.0)
          rhs_sc[isc+icv*solver->nsc_transport] -= mf[ibf] * scbc_vals[isc];
        else
          rhs_sc[isc+icv*solver->nsc_transport] -= mf[ibf] * solver->transport_scalar_vec[isc][icv];
      }
    }
  }

};


enum SpongeType {
  SPONGE_LX,
  SPONGE_LY,
  SPONGE_LZ,
  SPONGE_LRX,
  SPONGE_LRZ,
  SPONGE_LTH,
  SPONGE_UNDEFINED
};

//===========================================================
// auxiliary functions for the bc calculations
// XXX these can be moved elsewhere ..
//===========================================================

// for ideal gas, nscbc boundary conditions...
inline int calcNscbcIgSubsonicInletRhs(double* rhs, const double unit_n[3], const double delta,
				       const double rho_bc, const double u_bc[3], const double p_bc,
				       const double rho_cv, const double u_cv[3], const double p_cv,
				       const double gamma) __attribute__((always_inline)) ;

inline int calcNscbcIgSubsonicInletRhs(double* rhs, const double unit_n[3], const double delta,
				       const double rho_bc, const double u_bc[3], const double p_bc,
				       const double rho_cv, const double u_cv[3], const double p_cv,
				       const double gamma) {
  const double u0     = DOT_PRODUCT(u_cv, unit_n);
  const double un_bc  = DOT_PRODUCT(u_bc, unit_n);
  const double dudn   = (un_bc - u0)/delta;
  const double dpdn   = (p_bc -  p_cv)/delta;
  const double sos_bc = sqrt(gamma*p_bc/rho_bc);
  const double lam5   = u0 + sos_bc; 
  if (lam5 < 0.0) return -1;
  const double L5     = lam5* ( dpdn + rho_bc*sos_bc*dudn);
  
  // since the mdot, temperature are fixed, the following relations hold
  const double Ma   = un_bc/sos_bc; //signed..
  double L1         = -(1.0 +Ma*gamma)/(Ma-1.0)*L5;
  L1               /=  (1.0 + Ma/(Ma-1.0)*(gamma-1.0));
  const double L2   = 0.5*(gamma-1.0)*(L5+L1);
  rhs[0] = -(L2 + 0.5*(L5+L1))/sos_bc/sos_bc;
  return 0;
}

inline void calcNscbcIgSubsonicLIRhs(double * rhs, const double unit_n[3], const double delta,
                                     const double rho_bc, const double u_bc[3], const double p_bc,
                                     const double rho_cv, const double u_cv[3], const double p_cv,
                                     const double gamma);

inline void calcNscbcIgSubsonicLIRhs(double * rhs, const double unit_n[3], const double delta,
                                      const double rho_bc, const double u_bc[3], const double p_bc,
                                      const double rho_cv, const double u_cv[3], const double p_cv,
                                      const double gamma) {

  // In this bc, we keep rho, u, v, and w fixed
  // This gives that L2 + 0.5*(L5 + L1) = 0
  //                 L1 = L5
  //                 L3 = 0
  //                 L4 = 0

  const double u0     = DOT_PRODUCT(u_cv, unit_n);
  const double un_bc  = DOT_PRODUCT(u_bc, unit_n);
  const double dudn   = (un_bc - u0)/delta;
  const double dpdn   = (p_bc - p_cv)/delta;
  const double sos_bc = sqrt(gamma*p_bc/rho_bc);
  const double lam5   = u0 + sos_bc; assert( lam5 >= 0.0);
  const double L5     = lam5* ( dpdn + rho_bc*sos_bc*dudn);

  // Above this is computed as L5 instead of L1
  //const double lam1   = u0 - sos_bc;
  //if (lam1 > 0.0)
  // cout << "Woah! lam1 = " << lam1 << ", u0 = " << u0 << ", sos_bc = " << sos_bc << "\n" << flush;
  //assert(lam1 <= 0.0);
  //const double L1     = lam1*(dpdn - rho_bc*sos_bc*dudn);

  //since the rho and u values are fixed, the following relations hold
  //const double Ma = un_bc/sos_bc; //signed...

  const double L1       = L5;
  rhs[1]          = -0.5*(L1 + L5);

}

// for ideal gas, nscbc boundary conditions...
inline void calcNscbcIgSubsonicInletSoftRhs(double* rhs, const double unit_n[3], const double delta,
                                        const double rho_bc, const double u_bc[3], const double p_bc,
                                        const double rho_cv, const double u_cv[3], const double p_cv,
                                        const double gamma, const double T_nom, const double t_relax,
                                        const double R_gas) __attribute__((always_inline)) ;

inline void calcNscbcIgSubsonicInletSoftRhs(double* rhs, const double unit_n[3], const double delta,
                                        const double rho_bc, const double u_bc[3], const double p_bc,
                                        const double rho_cv, const double u_cv[3], const double p_cv,
                                        const double gamma, const double T_nom, const double t_relax,
                                        const double R_gas) {
  const double u0     = DOT_PRODUCT(u_cv, unit_n);
  const double un_bc  = DOT_PRODUCT(u_bc, unit_n);
  const double dudn   = (un_bc - u0)/delta;
  const double dpdn   = (p_bc -  p_cv)/delta;
  const double sos_bc = sqrt(gamma*p_bc/rho_bc);
  const double lam5   = u0 + sos_bc; assert( lam5 >= 0.0);
  const double L5     = lam5* ( dpdn + rho_bc*sos_bc*dudn);

  // note that the mach number here is signed ( likely < 0 )
  const double Ma   = un_bc/sos_bc;

  // compute a soft relaxed coefficient for the acoustic wave ..
  const double T_bc = p_bc/(rho_bc*R_gas);
  const double L1   = 0.5*(rho_bc + rho_cv)*R_gas*(T_bc - T_nom)/t_relax;

  // now compute the L2 characteristic to enforce that the mdot is constant..
  const double L2   = -1.0/(2.0*Ma)*( (Ma+1.0)*L5 + (Ma-1.0)*L1);

  // advance the rhs for the density and pressure ..
  rhs[0] = -(L2 + 0.5*(L5+L1))/sos_bc/sos_bc;
  rhs[1] = -0.5*(L5+L1);
}

inline double calcNscbcIgOutletRhs(double* rhs, const double unit_n[3], const double delta,
                                   const double rho_bc, const double u_bc[3], const double p_bc,
                                   const double rho_cv, const double u_cv[3], const double p_cv,
                                   const double gamma, const double L_ref, const double sigma,
                                   const double p_ref) __attribute__((always_inline)) ;

inline double calcNscbcIgOutletRhs(double* rhs, const double unit_n[3], const double delta,
                                   const double rho_bc, const double u_bc[3], const double p_bc,
                                   const double rho_cv, const double u_cv[3], const double p_cv,
                                   const double gamma, const double L_ref, const double sigma,
                                   const double p_ref) {

  const double un_bc  = max(0.0,DOT_PRODUCT(u_bc,unit_n));
  const double u0     = max(0.0,DOT_PRODUCT(u_cv,unit_n));
  const double sos_bc = sqrt(gamma*p_bc/rho_bc);
  const double Ma     = un_bc/sos_bc; // this is a signed Mach number..
  const double dudn   = (un_bc - u0)/delta;
  const double dpdn   = (p_bc - p_cv)/delta;
  const double L5     = (un_bc + sos_bc) * (dpdn + rho_bc*sos_bc*dudn);

  double dudn_[3];
  for (int i =0; i < 3; ++i)
    dudn_[i] = (u_bc[i] - u_cv[i])/delta;

  // the incoming acoustic characteristic ..
  double L1;
  if ( sigma >= 0.0) {
    L1           = sigma* (1.0-Ma*Ma) * sos_bc/ L_ref * (p_bc - p_ref);
  } else {
    // a negative sigma indicates a perfectly reflective condition so L1 = -L5;
    L1 = -L5;
  }

  if ( Ma >= 0.0 ) {
    const double drhodn = (rho_bc - rho_cv)/delta;
    const double L2     = un_bc * (sos_bc*sos_bc* drhodn - dpdn);

    if ( Ma >= 1.0 ) {
      // all the characteristics are going out...
      L1 = sos_bc*(Ma - 1.0)*( dpdn - rho_bc*sos_bc*dudn);
    }

    rhs[0] = -(L2 + 0.5*(L5+L1) )/sos_bc/sos_bc; // rho_bc
    rhs[1] = -0.5*(L5+L1);                       // p_bc

    // special combination of the velocity bcs to avoid coordinate decomposition...
    for (int i =0; i < 3; ++i) {
      rhs[2+i] = -(rho_bc*un_bc*dudn_[i] +
                   0.5*rho_bc*sos_bc*(1.0-Ma)*dudn*unit_n[i])/rho_bc
      -0.5*(1.0+Ma)*dpdn*unit_n[i]/rho_bc
      +0.5*L1*unit_n[i]/sos_bc/rho_bc;
    }
  }
  else {
    cout << "got Ma: " << Ma << endl;
    assert(0); // un_bc was maxed to 0.0... prohibiting any backflow.
  }

  return un_bc;
}

inline double calcNscbcIgOutletMdotWeakRhs(double* rhs, const double unit_n[3], const double delta,
                                           const double rho_bc, const double u_bc[3], const double p_bc,
                                           const double rho_cv, const double u_cv[3], const double p_cv,
                                           const double gamma, 
                                           const double rhou_avg, const double u_avg, const double rhou_target, 
                                           const double T_relax) { 

  const double un_bc  = max(0.0,DOT_PRODUCT(u_bc,unit_n));
  const double u0     = max(0.0,DOT_PRODUCT(u_cv,unit_n));
  const double sos_bc = sqrt(gamma*p_bc/rho_bc);
  const double Ma     = un_bc/sos_bc; // this is a signed Mach number..
  const double dudn   = (un_bc - u0)/delta;
  const double dpdn   = (p_bc - p_cv)/delta;
  const double L5     = (un_bc + sos_bc) * (dpdn + rho_bc*sos_bc*dudn);

  double dudn_[3];
  for (int i =0; i < 3; ++i) 
    dudn_[i] = (u_bc[i] - u_cv[i])/delta;

  // the incoming acoustic characteristic .. 
  // L1 = sigma*u_avg*(rhou_avg - rhou_target)/T, T \approx n*FTT
  
  //double L1    = abs(u_avg)*(rhou_avg - rhou_target)/T_relax;
  double L1      =  -sos_bc*(rhou_avg - rhou_target)/T_relax;

  if ( Ma >= 0.0 ) {
    const double drhodn = (rho_bc - rho_cv)/delta;
    const double L2     = un_bc * (sos_bc*sos_bc* drhodn - dpdn);
    
    if ( Ma >= 1.0 ) {
      // all the characteristics are going out...
      L1 = sos_bc*(Ma - 1.0)*( dpdn - rho_bc*sos_bc*dudn);
    }
    
    rhs[0] = -(L2 + 0.5*(L5+L1) )/sos_bc/sos_bc; // rho_bc
    rhs[1] = -0.5*(L5+L1);                       // p_bc
    
    // special combination of the velocity bcs to avoid coordinate decomposition...
    for (int i =0; i < 3; ++i) {
      rhs[2+i] = -(rho_bc*un_bc*dudn_[i] +
                   0.5*rho_bc*sos_bc*(1.0-Ma)*dudn*unit_n[i])/rho_bc
      -0.5*(1.0+Ma)*dpdn*unit_n[i]/rho_bc
      +0.5*L1*unit_n[i]/sos_bc/rho_bc;
    }
  } 
  else {
    cout << "got Ma: " << Ma << endl;
    assert(0); // un_bc was maxed to 0.0... prohibiting any backflow.
  }
  
  return un_bc;
}

inline double calcNscbcIgOutletRhsMdot(double* rhs, const double unit_n[3], const double delta,
                                       const double rho_bc, const double u_bc[3], const double p_bc,
                                       const double rho_cv, const double u_cv[3], const double p_cv,
                                       const double gamma) {

  const double un_bc  = max(0.0,DOT_PRODUCT(u_bc,unit_n));
  const double u0     = max(0.0,DOT_PRODUCT(u_bc,unit_n));
  const double sos_bc = sqrt(gamma*p_bc/rho_bc);
  const double Ma     = un_bc/sos_bc; // this is a signed Mach number..
  const double dudn   = (un_bc - u0)/delta;
  const double dpdn   = (p_bc - p_cv)/delta;
  const double L5     = (un_bc + sos_bc) * (dpdn + rho_bc*sos_bc*dudn);
  const double drhodn = (rho_bc - rho_cv)/delta;
  const double L2     = un_bc * (sos_bc*sos_bc* drhodn - dpdn);

  // the incoming acoustic characteristic ..
  double L1 = L5*(1.0+Ma)/(1.0-Ma) + 2.0*Ma/(1.0-Ma)*L2;

  if ( Ma >= 0.0 ) {

    if ( Ma >= 1.0 ) {
      // all the characteristics are going out...
      L1 = sos_bc*(Ma - 1.0)*( dpdn - rho_bc*sos_bc*dudn);
    }

    rhs[0] = -(L2 + 0.5*(L5+L1) )/sos_bc/sos_bc; // rho_bc
    rhs[1] = -0.5*(L5+L1);                       // p_bc

    // since we are enforcing a constant mdot by augmenting the velocity
    // data, there is no need to populate the rhs for the velocity field..
    for (int i =0; i < 3; ++i)
      rhs[2+i] = 0.0;

  } else {
    cout << "got Ma: " << Ma << endl;
    assert(0); // un_bc was maxed to 0.0... prohibiting any backflow.
  }

  return un_bc;
}

template <class T>
inline void force_bf_slip_wall(double (*f_bf)[9],T* bc) {

  for (int ibf = 0; ibf < bc->zone_ptr->nbf; ++ibf) {

    const int icv0 = bc->zone_ptr->cvobf[ibf];

    FOR_I3 {
      f_bf[ibf][3*0+i] = 0.0;
      f_bf[ibf][3*1+i] = bc->solver->p[icv0]*bc->zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = 0.0;
    }
  }
}

template <class T>
inline void force_bf_wall(double (*f_bf)[9],T* bc) {

  for (int ibf = 0; ibf < bc->zone_ptr->nbf; ++ibf) {

    const int icv0 = bc->zone_ptr->cvobf[ibf];
    double mu_coeff = bc->solver->mu_lam[icv0]*bc->zone_ptr->area_over_delta_bf[ibf];

    FOR_I3 {
      f_bf[ibf][3*0+i] = 0.0;
      f_bf[ibf][3*1+i] = bc->solver->p[icv0]*bc->zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = mu_coeff*(bc->solver->u[icv0][i] - bc->u_bc[i]);
    }
  }
}

// once all bcs are updated with the mrf, then this function signature 
// will replace the force_bf_wall above.

template <class T>
inline void force_bf_wall_new(double (*f_bf)[9],T* bc) {

  for (int ibf = 0; ibf < bc->zone_ptr->nbf; ++ibf) {

    const int icv0 = bc->zone_ptr->cvobf[ibf];
    double mu_coeff = bc->solver->mu_lam[icv0]*bc->zone_ptr->area_over_delta_bf[ibf];

    FOR_I3 {
      f_bf[ibf][3*0+i] = 0.0;
      f_bf[ibf][3*1+i] = bc->solver->p[icv0]*bc->zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = mu_coeff*(bc->solver->u[icv0][i] - bc->u_bc[ibf][i]);
    }
  }
}

template <class T>
inline void force_bf_cbc(double (*f_bf)[9],T* bc) {
  for (int ibf = 0; ibf < bc->zone_ptr->nbf; ++ibf) {

    const int icv0 = bc->zone_ptr->cvobf[ibf];
    FOR_I3 {
      f_bf[ibf][3*0+i] = bc->mf[ibf]*bc->solver->u[icv0][i];
      f_bf[ibf][3*1+i] = bc->solver->p[icv0]*bc->zone_ptr->n_bf[ibf][i];
      f_bf[ibf][3*2+i] = 0.0;
    }
  }
}

template <class T>
inline void force_bf_wm(double (*f_bf)[9],T* bc) {

  for (int ibf = 0; ibf < bc->zone_ptr->nbf; ++ibf) {

    const int icv0 = bc->zone_ptr->cvobf[ibf];

    // tau_wall has already been computed at this point, but we'll assume
    // that the wall stress is instantaneously aligned with the velocity vector

    double u_mag = DOT_PRODUCT(bc->solver->u[icv0],bc->solver->u[icv0]);
    u_mag = sqrt( max(0.0, u_mag));

    if ( u_mag > 0.0) {

      FOR_I3 {
        f_bf[ibf][3*0+i] = 0.0;
        f_bf[ibf][3*1+i] = bc->solver->p[icv0]*bc->zone_ptr->n_bf[ibf][i];
        f_bf[ibf][3*2+i] = bc->tau_wall[ibf]*bc->solver->u[icv0][i]/u_mag*bc->zone_ptr->area_bf[ibf];
      }

    } else {

      FOR_I3 {
        f_bf[ibf][3*0+i] = 0.0;
        f_bf[ibf][3*1+i] = bc->solver->p[icv0]*bc->zone_ptr->n_bf[ibf][i];
        f_bf[ibf][3*2+i] = 0.0;
      }

    }
  }
}

enum WmType { 
  WM_ALGEBRAIC_ISOTHERMAL,
  WM_EQUILIBRIUM_ISOTHERMAL,
  WM_ALGEBRAIC_ADIABATIC,
  WM_EQUILIBRIUM_ADIABATIC
};


#endif
