#ifndef SIMPLEHOOKBC_HPP
#define SIMPLEHOOKBC_HPP

class SimpleHookBc : public SimpleIdealGasBc { 
public: 

  SimpleHookBc(BfZone* p, IdealGasSolver* s) : SimpleIdealGasBc(p,s) {}

  void addBoundaryFlux(IdealGasRhs* rhs) const { 

    // 
    // will supply a sinusoidal variation of the 
    // x-velocity at a constant pressure, temperature
    // in time.. 
    // 
    // the solver ptr in the bc class is a back reference 
    // to the full solver state (to get data living on cvs, etc)
    // 
    // the zone_ptr ptr in the bc class gives access to structures
    // meant to exist for this particular fazone.

    assert( solver != NULL);
    const double time = solver->time;

    const double u_ref = 0.1;
    const double u_amp = 0.01;
    const double p_ref   = solver->p_ref;
    const double rho_ref = solver->rho_ref;
    const double h_ref   = (p_ref/(solver->gamma-1.0) + p_ref)/rho_ref;
   

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

      const int icv   = zone_ptr->cvobf[ibf];
    
      // ====================================================
      //  this is the bc state specification here to edit
      //=====================================================

      const double yy = zone_ptr->x_bf[ibf][1];   // y coordinate of this face

      IdealGasState bf_state; // requires a definition of u,p,h,sp_vol
      bf_state.sp_vol      = 1.0/rho_ref;
      bf_state.p           = p_ref;
      bf_state.h           = h_ref;
      bf_state.u[0]        = u_ref + u_amp*sin(2.0*M_PI*time)*cos(2.0*M_PI*yy);
      bf_state.u[1]        = 0.0;
      bf_state.u[2]        = 0.0;

      //======================================================
      //   end of specification
      //======================================================

      IdealGasRhs flux;
      calcRiemannFlux(flux,zone_ptr->n_bf[ibf], solver->cv_light[icv], 
                      bf_state,solver->gamma,solver->gamma);

      mf[ibf]       = flux.rho;
      
      rhs[icv].rho -= flux.rho;
      for (int i =0; i < 3; ++i) 
        rhs[icv].rhou[i] -= flux.rhou[i];
      rhs[icv].rhoE -= flux.rhoE;
    }
  }
};

class SimpleHookRandomBc : public SimpleIdealGasBc { 
public: 

  double tke;
  double * u_buf_p;

  SimpleHookRandomBc(BfZone* p, IdealGasSolver* s) : SimpleIdealGasBc(p,s) {
  
    tke = getDoubleParam("TRIP_TKE");
    u_buf_p = NULL; zone_ptr->registerBfData(u_buf_p,"u_buf",READWRITE_DATA);

  }

  ~SimpleHookRandomBc() { 
    DELETE(u_buf_p);
  }

  void initData() { 

    // call the parent to alloc aux memory..

    SimpleIdealGasBc::initData(); 

    assert( u_buf_p == NULL); u_buf_p = new double[zone_ptr->nbf];
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) 
      u_buf_p[ibf] = 0.0;

  }

  void addBoundaryFlux(IdealGasRhs* rhs) const { 

    // 
    // will supply a sinusoidal variation of the 
    // x-velocity at a constant pressure, temperature
    // in time.. 
    // 
    // the solver ptr in the bc class is a back reference 
    // to the full solver state (to get data living on cvs, etc)
    // 
    // the zone_ptr ptr in the bc class gives access to structures
    // meant to exist for this particular fazone.

    assert( solver != NULL);
    const double time = solver->time;

    double * u_buf       = new double[zone_ptr->nbf];

    double my_buf[2] = {0.0,0.0};

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

      const int icv     = zone_ptr->cvobf[ibf];
      const double u_pr = sqrt(2.0/3.0*tke/solver->rho[icv]); // tke = 3/2*rho*u'^2
      const double area = MAG(zone_ptr->n_bf[ibf]);

      u_buf[ibf]        = u_pr*(2.0*double(rand())/double(RAND_MAX)-1.0);
      my_buf[0]        += area;
      my_buf[1]        += area*u_buf[ibf];
    }

    //MiscUtils::dumpRange(u_buf,zone_ptr->nbf, "u_trip_inst: " + getName() );

    // remove the mean across the zone for the uncorrelated white noise... 
    
    double buf[2];
    MPI_Allreduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,mpi_comm);

    const double um = buf[1]/buf[0];
    const double eps = 0.01;

    assert( u_buf_p != NULL);
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
      //u_buf_p[ibf] = (1.0-eps)*u_buf_p[ibf] + eps*(u_buf[ibf]);
      u_buf_p[ibf] = (1.0-eps)*u_buf_p[ibf] + eps*(u_buf[ibf] - um);
    }

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 

      const int icv     = zone_ptr->cvobf[ibf];

      // ====================================================
      //  this is the bc state specification here to edit
      //=====================================================

      const double yy = zone_ptr->x_bf[ibf][1];   // y coordinate of this face

      const double uu = u_buf_p[ibf];
      const double area = MAG(zone_ptr->n_bf[ibf]);
      double unit_n[3];
      for (int i = 0; i < 3; ++i) 
        unit_n[i] = zone_ptr->n_bf[ibf][i]/area;

      IdealGasState bf_state; // requires a definition of u,p,h,sp_vol
      bf_state.sp_vol      = 1.0/solver->rho[icv];
      bf_state.p           = solver->p[icv]; 
      bf_state.h           = (solver->p[icv]/(solver->gamma-1.0) + solver->p[icv])/solver->rho[icv]; 
      bf_state.u[0]        = uu*unit_n[0]; 
      bf_state.u[1]        = uu*unit_n[1];
      bf_state.u[2]        = uu*unit_n[2];

      //======================================================
      //   end of specification
      //======================================================

      IdealGasRhs flux;
      calcRiemannFlux(flux,zone_ptr->n_bf[ibf], solver->cv_light[icv], 
                      bf_state,solver->gamma,solver->gamma);

      mf[ibf]       = flux.rho;
      
      rhs[icv].rho -= flux.rho;
      for (int i =0; i < 3; ++i) 
        rhs[icv].rhou[i] -= flux.rhou[i];
      rhs[icv].rhoE -= flux.rhoE;
    }


    //MiscUtils::dumpRange(u_buf_p,zone_ptr->nbf, "u_trip: " + getName() );

    delete[] u_buf;
  }

  void force_bf(double (*f_bf)[9]) const {force_bf_slip_wall(f_bf,this);} 
};

class SimpleHookSolver : public IdealGasSolver { 
public: 

  IdealGasBc* initHookBc(BfZone* p) {
  
    //
    // assumes that the input file required as bc specification 
    // of the form.. 
    // 
    //   x0 = HOOK 
    // 

    assert( p->getName() == "x0");
    return new SimpleHookBc(p,this);

  }
};

#endif
