#include "IdealGasSolver.hpp"
#include "RkWeights.hpp"
#include "wm/AlgebraicWM.hpp"
#include "wm/EquilibriumWM.hpp"
#include "SpongeCommon.hpp"
#include "it/InflowTurbulenceIG.hpp"
#include "bcprofile/ProfileReader.hpp"
#include "DataExchanger.hpp"
#include "IdealGasSolverFlux.hpp"
#include "oneDProfile.hpp"

using namespace RkWeights;

void SlipWallBc::addBoundaryFlux(IdealGasRhs* rhs) const {

  assert( solver != NULL);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv = zone_ptr->cvobf[ibf];
    for (int i =0; i < 3; ++i)
      rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];
  }

}

void WallAdiabatic::addBoundaryFlux(IdealGasRhs* rhs) const {

  assert( solver != NULL);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv           = zone_ptr->cvobf[ibf];
    const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];

    double visc_work = 0.0;
    for (int i =0; i < 3; ++i) {

      double tauijnj    = visc_coeff*(solver->u[icv][i] - u_bc[ibf][i]);
      rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i] + tauijnj;
      visc_work        += u_bc[ibf][i]*tauijnj;

    }

    rhs[icv].rhoE      -= visc_work;

  }

}

void WallAdiabatic::query(const string& param_str) const {

  double my_buf[3] = {0.0,0.0,0.0};
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

      const double y1         = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
      const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];

      double du_mag = 0.0;
      for (int i =0; i < 3; ++i)
        du_mag  += (solver->u[icv][i] - u_bc[ibf][i])*(solver->u[icv][i] - u_bc[ibf][i]);

      const double tau_wall   = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
      const double y_plus      = y1*sqrt(tau_wall*solver->rho[icv])/solver->mu_lam[icv];

      my_buf[0]              += zone_ptr->area_bf[ibf];
      my_buf[1]              += zone_ptr->area_bf[ibf]*tau_wall;
      my_buf[2]              += zone_ptr->area_bf[ibf]*y_plus;

    }
  }

  double buf[3];
  MPI_Reduce(my_buf,buf,3,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if ( mpi_rank == 0 ) {

    cout << zone_ptr->getName() << " time, int(tau_wall)dA, avg y_plus = "
         << solver->time << "     "
         << buf[1] << "     "
         << buf[2]/buf[0] << endl;

  }
}

void computeForceValues(CtiRegister::CtiData& v,BfZone* zone_ptr, double (*f_bf)[9],const int component,const string type) {

  double wgt[3] = {1.0,1.0,1.0};
  if (type == "conv") {
    wgt[1] = wgt[2] = 0.0;
  }
  else if (type == "press") {
    wgt[0] = wgt[2] = 0.0;
  }
  else if (type == "visc") {
    wgt[1] = wgt[0] = 0.0;
  }

  if (component == 0) {
    // sum all
    double (*tmp)[3] = zone_ptr->createBfD2Data(v);

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      // NOTE: cannot put into dp form without access to p_inf, which is passed in WRITE_FORCE case
      // FOR_I3 f_bf[ibf][3+i] -= solver->p_inf*zone_ptr->n_bf[ibf][i];

      // sum all components
      for (int i=0; i<3; ++i) {
        tmp[ibf][i] = wgt[0]*f_bf[ibf][i] + wgt[1]*f_bf[ibf][3+i] + wgt[2]*f_bf[ibf][6+i];
        tmp[ibf][i] /= zone_ptr->area_bf[ibf];
      }
    }
  }
  else {
    // component access; index already vetted
    assert((component > 0) && (component < 4));

    double *tmp = zone_ptr->createBfD1Data(v);
    int index = component-1;
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      // NOTE: cannot put into dp form without access to p_inf, which is passed in WRITE_FORCE case
      // FOR_I3 f_bf[ibf][3+i] -= solver->p_inf*zone_ptr->n_bf[ibf][i];

      tmp[ibf] = wgt[0]*f_bf[ibf][index] + wgt[1]*f_bf[ibf][3+index] + wgt[2]*f_bf[ibf][6+index];
      tmp[ibf] /= zone_ptr->area_bf[ibf];
    }
  }
}

void computeMomentValues(CtiRegister::CtiData& v,BfZone* zone_ptr, double (*f_bf)[9],const double (&origin)[3],const int component,const string type) {

  double wgt[3] = {1.0,1.0,1.0};
  if (type == "conv") {
    wgt[1] = wgt[2] = 0.0;
  }
  else if (type == "press") {
    wgt[0] = wgt[2] = 0.0;
  }
  else if (type == "visc") {
    wgt[1] = wgt[0] = 0.0;
  }

  if (component == 0) {
    // sum all
    double (*tmp)[3] = zone_ptr->createBfD2Data(v);

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      // NOTE: cannot put into dp form without access to p_inf, which is passed in WRITE_FORCE case
      // FOR_I3 f_bf[ibf][3+i] -= solver->p_inf*zone_ptr->n_bf[ibf][i];

      // sum all components
      double my_force[3];
      FOR_I3 my_force[i] = wgt[0]*f_bf[ibf][i] + wgt[1]*f_bf[ibf][3+i] + wgt[2]*f_bf[ibf][6+i];

      const double my_bf[3] = CROSS_PRODUCT(zone_ptr->x_bf[ibf],my_force);
      const double my_bf0[3] = CROSS_PRODUCT(origin,my_force);

      FOR_I3 tmp[ibf][i] = (my_bf[i]-my_bf0[i])/zone_ptr->area_bf[ibf];
    }
  }
  else {
    // component access; index already vetted
    assert((component > 0) && (component < 4));

    double *tmp = zone_ptr->createBfD1Data(v);
    int index = component-1;
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      // NOTE: cannot put into dp form without access to p_inf, which is passed in WRITE_FORCE case
      // FOR_I3 f_bf[ibf][3+i] -= solver->p_inf*zone_ptr->n_bf[ibf][i];

      // sum all components
      double my_force[3];
      FOR_I3 my_force[i] = wgt[0]*f_bf[ibf][i] + wgt[1]*f_bf[ibf][3+i] + wgt[2]*f_bf[ibf][6+i];

      const double my_bf[3] = CROSS_PRODUCT(zone_ptr->x_bf[ibf],my_force);
      const double my_bf0[3] = CROSS_PRODUCT(origin,my_force);

      tmp[ibf] = (my_bf[index]-my_bf0[index])/zone_ptr->area_bf[ibf];
    }
  }
}

CtiRegister::CtiDataError WallAdiabatic::funcEvalCtiData(CtiRegister::CtiData& v,
                                                         const string& name, list<CtiRegister::CtiData>& args,
                                                         const bool b_eval_func) {

  using namespace CtiRegister;

  if ( zone_ptr->isLeadingStrMatched(name) ) {

    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix.  we have to check for the actual data now..

    const string tau_wall_str = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str  = zone_ptr->getName() + ":" + "y_plus";
    const string force_str  = zone_ptr->getName() + ":" + "force";
    const string moment_str  = zone_ptr->getName() + ":" + "moment";
    // all of the following functions do not take any arguments -- check first

    if ( name == tau_wall_str) {

      double * tmp = zone_ptr->createBfD1Data(v);

      // if stats have been requested on this quantity, this function will be called
      // but its possible that the mu_lam, rho necessary have not been populated yet.
      // check the data flags for these values -- revisit if this impacts performance
      // during the solution run-time.

      if (b_eval_func) {
        if (args.size() == 0) {
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            const int icv           = zone_ptr->cvobf[ibf];
            const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];

            double du_mag = 0.0;
            for (int i =0; i < 3; ++i)
              du_mag  += (solver->u[icv][i] - u_bc[ibf][i])*(solver->u[icv][i] - u_bc[ibf][i]);

            tmp[ibf]                = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
          }
        }
        else if (args.size() == 1) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          double dir[3] = {0.0,0.0,0.0};

          if (arg->getType() == D3_DATA) {
            FOR_I3 dir[i] = arg->d3(i);
            NORMALIZE(dir);
          }
          else if (arg->getType() == D_DATA) {
            const int index = int(arg->d());
            if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;
            dir[index] = 1.0;
          }
          else return CTI_DATA_NOT_VALID;

          // loop faces and compute
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            const int icv           = zone_ptr->cvobf[ibf];
            const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];

            double local_u[3] = DIFF(solver->u[icv],u_bc[ibf]);
            tmp[ibf]          = visc_coeff*DOT_PRODUCT(dir,local_u)/zone_ptr->area_bf[ibf];
          }
        }
        else return CTI_DATA_ARG_COUNT;
      }

      return CTI_DATA_OK;

    }
    else if ( name == yp_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          const int icv           = zone_ptr->cvobf[ibf];
          const double y1         = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
          const double visc_coeff = solver->mu_lam[icv]/y1;

          double du_mag = 0.0;
          for (int i =0; i < 3; ++i)
            du_mag  += (solver->u[icv][i] - u_bc[ibf][i])*(solver->u[icv][i] - u_bc[ibf][i]);

          const double tau_wall   = visc_coeff*sqrt(du_mag);
          tmp[ibf]                = y1* sqrt(tau_wall*solver->rho[icv])/solver->mu_lam[icv];

        }

      }

      return CTI_DATA_OK;

    }
    else if ( name == force_str) {

          if (b_eval_func) {
            if (args.size() == 0) {

              // pressure, viscous, convective force components
              double (*f_bf)[9] = new double[zone_ptr->nbf][9];
              force_bf(f_bf);

              computeForceValues(v,zone_ptr,f_bf,0,"total");

              delete[] f_bf;
            }
            else if (args.size() == 1) {
              list<CtiRegister::CtiData>::iterator arg = args.begin();

              int index = -1;
              if (arg->getType() == D_DATA) {
                index = int(arg->d());
                if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;

                double (*f_bf)[9] = new double[zone_ptr->nbf][9];
                force_bf(f_bf);

                computeForceValues(v,zone_ptr,f_bf,(index+1),"total");

                delete[] f_bf;
              }
              else if (arg->getType() == S_DATA) {

                // convective, pressure, viscous force components
                if ((arg->s() != "conv") && (arg->s() != "visc") && (arg->s() != "press")) return CTI_DATA_NOT_VALID;

                double (*f_bf)[9] = new double[zone_ptr->nbf][9];
                force_bf(f_bf);

                computeForceValues(v,zone_ptr,f_bf,0,arg->s());

                delete[] f_bf;
              }
              else return CTI_DATA_NOT_VALID;
            }
            else if (args.size() == 2) {
              list<CtiRegister::CtiData>::iterator arg = args.begin();

              if (arg->getType() != S_DATA) return CTI_DATA_NOT_VALID;
              const string component = arg->s();
              arg++;

              int index = -1;
              if (arg->getType() != D_DATA) return CTI_DATA_NOT_VALID;
              index = int(arg->d());
              if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;

              if ((component != "conv") && (component != "visc") && (component != "press")) return CTI_DATA_NOT_VALID;

              double (*f_bf)[9] = new double[zone_ptr->nbf][9];
              force_bf(f_bf);

              computeForceValues(v,zone_ptr,f_bf,(index+1),component);

              delete[] f_bf;
            }
            else return CTI_DATA_ARG_COUNT;
          }

          return CTI_DATA_OK;

    }
    else if ( name == moment_str) {

          if (b_eval_func) {
            double origin[3] = {0.0,0.0,0.0};
            if (args.size() == 0) {

              // convective, pressure, viscous moment components
              double (*f_bf)[9] = new double[zone_ptr->nbf][9];
              force_bf(f_bf);

              computeMomentValues(v,zone_ptr,f_bf,origin,0,"total");

              delete[] f_bf;
            }
            else if (args.size() == 1) {
              list<CtiRegister::CtiData>::iterator arg = args.begin();

              if (arg->getType() == D_DATA) {
                int index = int(arg->d());
                if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;

                double (*f_bf)[9] = new double[zone_ptr->nbf][9];
                force_bf(f_bf);

                computeMomentValues(v,zone_ptr,f_bf,origin,(index+1),"total");

                delete[] f_bf;
              }
              else if (arg->getType() == S_DATA) {

                // convective, pressure, viscous force components
                if ((arg->s() != "conv") && (arg->s() != "visc") && (arg->s() != "press")) return CTI_DATA_NOT_VALID;

                double (*f_bf)[9] = new double[zone_ptr->nbf][9];
                force_bf(f_bf);

                computeMomentValues(v,zone_ptr,f_bf,origin,0,arg->s());

                delete[] f_bf;
              }
              else return CTI_DATA_NOT_VALID;
            }
            else if (args.size() == 3) {
              list<CtiRegister::CtiData>::iterator arg = args.begin();

              FOR_I3 {
                if (arg->getType() == D_DATA) {
                  origin[i] = arg->d();
                  arg++;
                }
                else return CTI_DATA_NOT_VALID;
              }

              double (*f_bf)[9] = new double[zone_ptr->nbf][9];
              force_bf(f_bf);

              computeMomentValues(v,zone_ptr,f_bf,origin,0,"total");

              delete[] f_bf;
            }
            else if (args.size() == 4) {
              list<CtiRegister::CtiData>::iterator arg = args.begin();

              if (arg->getType() != S_DATA) return CTI_DATA_NOT_VALID;
              const string component = arg->s();
              arg++;

              FOR_I3 {
                if (arg->getType() == D_DATA) {
                  origin[i] = arg->d();
                  arg++;
                }
                else return CTI_DATA_NOT_VALID;
              }

              if ((component != "conv") && (component != "visc") && (component != "press")) return CTI_DATA_NOT_VALID;

              double (*f_bf)[9] = new double[zone_ptr->nbf][9];
              force_bf(f_bf);

              computeMomentValues(v,zone_ptr,f_bf,origin,0,component);

              delete[] f_bf;
            }
            else return CTI_DATA_ARG_COUNT;
          }

          return CTI_DATA_OK;

    }
  }

  return IdealGasBc::funcEvalCtiData(v,name,args,b_eval_func);
}

void WallAdiabatic::initData() {

  assert( u_bc == NULL); u_bc = new double[zone_ptr->nbf][3];

  // default behavior (which is overridden below..) is u=0 wall

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    for (int i = 0; i < 3; ++i)
      u_bc[ibf][i] = 0.0;
  }

  Param * param = getParam(getName());
  assert( param->getString(0) == "WALL_ADIABATIC");

  int iarg = 1;
  while ( iarg < param->size()) {

    string token = param->getString(iarg++);
    if ( token == "TRANSLATE") {

      double u_wall[3]; FOR_I3 u_wall[i] = param->getDouble(iarg++);
      for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
        for (int i = 0; i  < 3; ++i)
          u_bc[ibf][i] = u_wall[i];
      }

    } else if ( (token == "STATIONARY") || (token == "STATIONARY_FRAME")) {

      if ( solver->frame_rotation != NULL) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          for (int i = 0; i < 3; ++i) {

            u_bc[ibf][0] = -solver->frame_rotation[1]*zone_ptr->x_bf[ibf][2] +
                            solver->frame_rotation[2]*zone_ptr->x_bf[ibf][1];

            u_bc[ibf][1] = -solver->frame_rotation[2]*zone_ptr->x_bf[ibf][0] +
                            solver->frame_rotation[0]*zone_ptr->x_bf[ibf][2];

            u_bc[ibf][2] = -solver->frame_rotation[0]*zone_ptr->x_bf[ibf][1] +
                            solver->frame_rotation[1]*zone_ptr->x_bf[ibf][0];
          }
        }

      }
    }
  }
}


void WmExchange::initData() {

  assert( tau_wall == NULL); tau_wall = new double[zone_ptr->nbf];
  assert( q_wall  == NULL);  q_wall   = new double[zone_ptr->nbf];
  assert( T_wall  == NULL);  T_wall   = new double[zone_ptr->nbf];
  assert( u1    == NULL);    u1       = new double[zone_ptr->nbf];
  assert( T1    == NULL);    T1       = new double[zone_ptr->nbf];

  Param * param = getParam(getName());

  del_exch_ratio = 3.0;

  if ( param->getString(0) == "WM_EXCHANGE_ADIABATIC") {

    if ( mpi_rank == 0)
      cout << " > wm exchange adiabatic for zone: " << this->getName() << endl;

    T_wall_ = -1.0; // to denote adiabatic conditions.

    if (param->size() > 1)
      del_exch_ratio = param->getDouble(1);

  } else if ( param->getString(0) == "WM_EXCHANGE_ISOTHERMAL") {

    if ( mpi_rank == 0)
      cout << " > wm exchange isothermal for zone: " << this->getName() << endl;

    T_wall_ = param->getDouble(1);
    assert( T_wall_ > 0.0);

    if (param->size() > 2)
      del_exch_ratio = param->getDouble(2);

  }

  double (*x_ex)[3] = new double[zone_ptr->nbf][3];
  unit_n            = new double[zone_ptr->nbf][3];

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    double nmag = MAG(zone_ptr->n_bf[ibf]);

    if ( nmag > 0.0) {

      for (int i =0; i <3; ++i)
        unit_n[ibf][i] = zone_ptr->n_bf[ibf][i]/nmag;

      const double delta = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];

      for (int i = 0; i < 3; ++i)
        x_ex[ibf][i] = zone_ptr->x_bf[ibf][i] - del_exch_ratio*unit_n[ibf][i]*delta;

    } else {

      const int icv = zone_ptr->cvobf[ibf];

      for (int i = 0; i < 3; ++i) {

        unit_n[ibf][i] = 0.0;
        x_ex[ibf][i]   = solver->x_cv[icv][i];

      }
    }

    tau_wall[ibf] = 0.0;
    q_wall[ibf]   = 0.0;
    T_wall[ibf]   = 0.0;
    u1[ibf]       = 0.0;
    T1[ibf]       = 0.0;
  }

  exchanger = new DataExchanger(solver);
  exchanger->init(x_ex,zone_ptr->nbf);

  delete[] x_ex;

}

void WmExchange::initialHook() {

  if ( ( !zone_ptr->checkDataFlag("u1") ) ||
       ( !zone_ptr->checkDataFlag("tau_wall")) ||
       ( !zone_ptr->checkDataFlag("T1")) ||
       ( !zone_ptr->checkDataFlag("q_wall")) ) {

    if ( mpi_rank == 0 )
      cout << " > " << zone_ptr->getName() << "  resetting everything... " << endl;


    const double cp = solver->gamma*solver->R_gas/(solver->gamma - 1.0);

    if ( T_wall_ < 0.0) {

      // adiabatic ..

      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

        const int icv = zone_ptr->cvobf[ibf];
        u1[ibf]       = MAG(solver->u[icv]);
        T1[ibf]       = solver->T[icv];
        T_wall[ibf]   = solver->T[icv] + 0.5*u1[ibf]*u1[ibf]/cp;

        if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

          const double y1 = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
          tau_wall[ibf]   = solver->mu_lam[icv] * MAG(solver->u[icv])/y1;

        } else {

          tau_wall[ibf]   = 0.0;

        }

        q_wall[ibf]     = 0.0;

      }

    } else {

      // isothermal..


      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

        const int icv = zone_ptr->cvobf[ibf];
        u1[ibf]       = MAG(solver->u[icv]);
        T1[ibf]       = solver->T[icv];
        T_wall[ibf]   = T_wall_;

        if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

          const double y1 = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
          tau_wall[ibf]   = solver->mu_lam[icv] * MAG(solver->u[icv])/y1;
          q_wall[ibf]     = cp*solver->loc_lam[icv]*(T1[ibf] - T_wall_)/y1;

        } else {

          tau_wall[ibf]   = 0.0;
          q_wall[ibf]     = 0.0;

        }


      }

    }
  }
}

int WmExchange::storeData(double * buf) const {

  // the wall model needs to store tau_wall
  // everything else gets recomputed from nearby fluid state

  if (buf) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      buf[ibf] = tau_wall[ibf];
    }
  }

  return zone_ptr->nbf;

}

int WmExchange::restoreData(double * buf) {

  assert(buf);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    tau_wall[ibf] = buf[ibf];
  }

  return zone_ptr->nbf;

}

WmExchange::~WmExchange() {

  DELETE(u1);
  DELETE(T1);
  DELETE(tau_wall);
  DELETE(q_wall);
  DELETE(T_wall);
  DELETE(unit_n);

  if ( exchanger) delete exchanger;

}

CtiRegister::CtiDataError WmExchange::funcEvalCtiData(CtiRegister::CtiData& v,
    const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if ( zone_ptr->isLeadingStrMatched(name) ) {

    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix.  we have to check for the actual data now..

    // for consistency, we are supporting a tau_wall(), y_plus()
    // function to have the same syntax as all of the other wall zones..

    const string tau_wall_str = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str       = zone_ptr->getName() + ":" + "y_plus";
    const string qw_str       = zone_ptr->getName() + ":" + "q_wall";
    const string Tw_str       = zone_ptr->getName() + ":" + "T_wall";

    // all of the following functions do not take any arguments -- check first

    if ( name == tau_wall_str) {

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {

        if (args.size() == 0) {
          // use what was computed by wall model
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
            tmp[ibf] = tau_wall[ibf];
          }
        }
        else if (args.size() == 1) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          double dir[3] = {0.0,0.0,0.0};

          if (arg->getType() == D3_DATA) {
            FOR_I3 dir[i] = arg->d3(i);
            NORMALIZE(dir);
          }
          else if (arg->getType() == D_DATA) {
            const int index = int(arg->d());
            if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;
            dir[index] = 1.0;
          }
          else return CTI_DATA_NOT_VALID;

          // loop faces and compute
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            double unit_u[3];
            const int icv      = zone_ptr->cvobf[ibf];
            const double mag_u = MAG(solver->u[icv]);

            if ( mag_u > 0.0) {

              for (int i = 0; i < 3; ++i)
                unit_u[i] = solver->u[icv][i]/mag_u;

              tmp[ibf] = tau_wall[ibf]*DOT_PRODUCT(dir,unit_u);

            } else {

              tmp[ibf] = 0.0;

            }

          }
        }
        else
          return CTI_DATA_ARG_COUNT;
      }

      return CTI_DATA_OK;

    } else if ( name == yp_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          const int icv      = zone_ptr->cvobf[ibf];
          const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
          tmp[ibf]           = y1*sqrt(solver->rho[icv]*tau_wall[ibf])/solver->mu_lam[icv];
        }

      }

      return CTI_DATA_OK;

    } else if ( name == qw_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
          tmp[ibf] = q_wall[ibf];
      }

      return CTI_DATA_OK;

    } else if ( name == Tw_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
          tmp[ibf] = T_wall[ibf];
      }

      return CTI_DATA_OK;

    }
  }

  return IdealGasBc::funcEvalCtiData(v,name,args,b_eval_func);
}

void WmExchange::calcRhs(const double time, const int rk_stage) {

  if ( rk_stage == 1 ) {

    const double cp = solver->R_gas*solver->gamma/(solver->gamma-1.0);

    // XXX we are presently not exchanging the wall temperature
    // and density -- these are all low Mach number approximations...

    double (*u_tmp)[3] = new double[zone_ptr->nbf][3];

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
      for (int i = 0; i < 3; ++i)
        u_tmp[ibf][i] = HUGE_VAL;

    exchanger->exchangeCvDataWithGhosts(u_tmp,solver->u,solver->dudx);

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      const int icv = zone_ptr->cvobf[ibf];

      if ( zone_ptr->area_over_delta_bf[ibf] < 0.0 ) {
        tau_wall[ibf]  = 0.0;
        q_wall[ibf]    = 0.0;
        T_wall[ibf]    = solver->T[icv];
        continue;
      }

      double y1;

      if ( u_tmp[ibf][0] == HUGE_VAL) {

        assert( u_tmp[ibf][1] == HUGE_VAL);
        assert( u_tmp[ibf][2] == HUGE_VAL);

        //
        // if we are unable ot locate this point (because its potenitally
        // outside of the domain.. we'll revert to the nearby cell)
        //

        for (int i = 0; i < 3; ++i)
          u_tmp[ibf][i] = solver->u[icv][i];

        y1  = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];

      } else {

        y1 = del_exch_ratio*zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];

      }

      double u_par[3];
      const double un = DOT_PRODUCT(u_tmp[ibf],unit_n[ibf]);

      for (int i = 0; i < 3; ++i)
        u_par[i] = u_tmp[ibf][i] - un*unit_n[ibf][i];

      const double upar_mag = MAG(u_par);

      // if you wanted to implement the time filtering, you would do so here ...

      u1[ibf]    = upar_mag;


      // there are low mach number approximations buried in the following.  we need
      // to exchagne T,rho otherwise -- see above.

      T1[ibf] = solver->T[icv];

      if ( T_wall_ < 0.0)  { // adiabatic case ..

        const double u_tau = sqrt( tau_wall[ibf] / solver->rho[icv]);
        tau_wall[ibf]      = AlgebraicWM::solve_tau(u1[ibf], y1, solver->rho[icv], solver->mu_lam[icv], u_tau);
        q_wall[ibf]        = 0.0;

        T_wall[ibf]        = solver->T[icv] + 0.5*u1[ibf]*u1[ibf]/cp;

      } else {

        const double u_tau = sqrt(tau_wall[ibf] / solver->rho[icv]);
        tau_wall[ibf]      = AlgebraicWM::solve_tau(u1[ibf], y1, solver->rho[icv],solver->mu_lam[icv], u_tau);

        q_wall[ibf]        = AlgebraicWM::compute_q_wall_approx(solver->T[icv], T_wall_, tau_wall[ibf],
                                                                solver->rho[icv], solver->mu_lam[icv],
                                                                solver->loc_lam[icv], cp, y1);
        T_wall[ibf]        = T_wall_;

      }

    }


    delete[] u_tmp;
  }
}

void WmExchange::addBoundaryFlux(IdealGasRhs* rhs) const {

  assert( solver != NULL);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    // tau_wall has already been computed at this point, but we'll assume
    // that the wall stress is instantaneously aligned with the velocity vector
    // for the component that is parallel to the wall...

    double u_par[3] = {0.0,0.0,0.0};
    const double un = DOT_PRODUCT(solver->u[icv],unit_n[ibf]);

    for (int i = 0; i < 3; ++i)
      u_par[i]  = solver->u[icv][i] - un*unit_n[ibf][i];

    const double u_mag = MAG(u_par);

    if ( u_mag > 0.0) {

      for (int i =0; i < 3;++i) {
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];
        rhs[icv].rhou[i] -= tau_wall[ibf]*u_par[i]/u_mag*zone_ptr->area_bf[ibf];
      }


    } else {

      for (int i =0; i < 3; ++i)
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];

    }

    // XXX if u_bc \neq 0, there should be a u \dot \tau term added here..
    rhs[icv].rhoE -= q_wall[ibf]*zone_ptr->area_bf[ibf];


  }
}

void WmExchange::query(const string& param_str) const {

  // report an estimate of the local y^+ on the boundary and the integrated tau_wall

  double my_buf[4] = {0.0,0.0,0.0,0.0};
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

      const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
      const double u_tau = sqrt( tau_wall[ibf]/solver->rho[icv]);
      const double nu    = solver->mu_lam[icv] / solver->rho[icv];
      const double y_plus = y1*u_tau/nu;

      my_buf[0] += zone_ptr->area_bf[ibf];
      my_buf[1] += zone_ptr->area_bf[ibf]*tau_wall[ibf];
      my_buf[2] += zone_ptr->area_bf[ibf]*y_plus;
      my_buf[3] += zone_ptr->area_bf[ibf]*q_wall[ibf];

    }
  }

  double buf[4];
  MPI_Reduce(my_buf,buf,4,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if ( mpi_rank == 0 ) {

    cout << zone_ptr->getName() << " time, int tau_wall dA, int q_w dA, avg y_plus = "
      << solver->time << "    "
      << buf[1] << "    "
      << buf[3] << "    "
      << buf[2]/buf[0] << endl;

  }

}


void WmAdiabatic::initData() {

  // allocate and set u_bc
  assert( u_bc == NULL); u_bc = new double[zone_ptr->nbf][3];

  // default behavior (which is overridden below..) is u=0 wall

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    for (int i = 0; i < 3; ++i)
      u_bc[ibf][i] = 0.0;
  }

  Param * param = getParam(getName());
  assert( (param->getString(0) == "WM_ALG_ADIABATIC") || (param->getString(0) == "WM_EQ_ADIABATIC"));

  int iarg = 1;
  while ( iarg < param->size()) {

    string token = param->getString(iarg++);
    if ( token == "TRANSLATE") {

      double u_wall[3]; FOR_I3 u_wall[i] = param->getDouble(iarg++);
      for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
        for (int i = 0; i  < 3; ++i)
          u_bc[ibf][i] = u_wall[i];
      }

    } else if ( (token == "STATIONARY") || (token == "STATIONARY_FRAME")) {

      if ( solver->frame_rotation != NULL) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          for (int i = 0; i < 3; ++i) {

            u_bc[ibf][0] = -solver->frame_rotation[1]*zone_ptr->x_bf[ibf][2] +
                            solver->frame_rotation[2]*zone_ptr->x_bf[ibf][1];

            u_bc[ibf][1] = -solver->frame_rotation[2]*zone_ptr->x_bf[ibf][0] +
                            solver->frame_rotation[0]*zone_ptr->x_bf[ibf][2];

            u_bc[ibf][2] = -solver->frame_rotation[0]*zone_ptr->x_bf[ibf][1] +
                            solver->frame_rotation[1]*zone_ptr->x_bf[ibf][0];
          }
        }

      }
    }
    else if ( token == "RELAX") {
      relax = param->getDouble(iarg++);
    }
    else if ( token == "Y_PLUS0") {
      y_plus0 = param->getDouble(iarg++);
    }
    else if ( token == "TOL") {
      tol = param->getDouble(iarg++);
    }

  }


  assert( tau_wall == NULL); tau_wall = new double[zone_ptr->nbf];
  assert( u1    == NULL);    u1       = new double[zone_ptr->nbf];

  if ( wm_type == WM_EQUILIBRIUM_ADIABATIC) {

    assert( T_wall == NULL); T_wall   = new double[zone_ptr->nbf];

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      T_wall[ibf] = -1.0; // set to invalid
    }

  }

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    tau_wall[ibf] = 0.0;
    u1[ibf]       = 0.0;
  }

}

void WmAdiabatic::initialHook() {

  if ( ( !zone_ptr->checkDataFlag("u1") ) || (!zone_ptr->checkDataFlag("tau_wall")) ) {

    if ( mpi_rank == 0 )
      cout << " > " << zone_ptr->getName() << "  resetting tau_wall " << endl;

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      const int icv = zone_ptr->cvobf[ibf];
      double u_minus_ubc[3];
      FOR_I3 u_minus_ubc[i] = solver->u[icv][i] - u_bc[ibf][i];
      u1[ibf] = MAG(u_minus_ubc);

      // set an initial value for tau_wall, this will constitute the
      // initial guess for the solution when it is finally constructed.

      if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

        const double y1 = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        tau_wall[ibf] = solver->mu_lam[icv] * u1[ibf]/y1;

      }
      else {
        tau_wall[ibf] = 0.0;
        u1[ibf] = 0.0;
      }
    }

  }

  if ( T_wall && !zone_ptr->checkDataFlag("T_wall") ) {

    if ( mpi_rank == 0)
      cout << " > " << zone_ptr->getName() << " resetting T_wall " << endl;

    const double cp = solver->R_gas*solver->gamma/(solver->gamma - 1.0);

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      const int icv = zone_ptr->cvobf[ibf];
      //assert( solver->frame_rotation == NULL);  // need to modify the following equation if there's frame rotation...
      T_wall[ibf]   = solver->T[icv] + 0.5*u1[ibf]*u1[ibf]/cp;

    }

  }
}

int WmAdiabatic::storeData(double * buf) const {

  // the wall model needs to store tau_wall because it is used as
  // an initial guess. No need to store u1. It is not used...

  if (buf) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      buf[ibf] = tau_wall[ibf];
    }
    if (T_wall) {
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
        buf[ibf+zone_ptr->nbf] = T_wall[ibf];
      }
    }
  }

  if (T_wall) {
    // if allocated then this is EQ model and needs T_wall as well
    return 2*zone_ptr->nbf;
  } else {
    return zone_ptr->nbf;  // just tau_wall
  }

}

int WmAdiabatic::restoreData(double * buf) {

  assert(buf);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    tau_wall[ibf] = buf[ibf];
  }

  if (T_wall) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      T_wall[ibf] = buf[ibf+zone_ptr->nbf];
    }
    return 2*zone_ptr->nbf;
  }

  return zone_ptr->nbf;  // just tau_wall
}

WmAdiabatic::~WmAdiabatic() {

  DELETE(tau_wall);
  DELETE(u1);
  DELETE(T_wall);
  DELETE(u_bc);

}

CtiRegister::CtiDataError WmAdiabatic::funcEvalCtiData(CtiRegister::CtiData& v,
                                                       const string& name, list<CtiRegister::CtiData>& args,
                                                       const bool b_eval_func) {

  using namespace CtiRegister;

  if ( zone_ptr->isLeadingStrMatched(name) ) {

    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix.  we have to check for the actual data now..

    // for consistency, we are supporting a tau_wall(), y_plus()
    // function to have the same syntax as all of the other wall zones..

    const string tau_wall_str = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str  = zone_ptr->getName() + ":" + "y_plus";

    // all of the following functions do not take any arguments -- check first


    if ( name == tau_wall_str) {

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {

        if (args.size() == 0) {
          // use what was computed by wall model
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
            tmp[ibf] = tau_wall[ibf];
          }
        }
        else if (args.size() == 1) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          double dir[3] = {0.0,0.0,0.0};

          if (arg->getType() == D3_DATA) {
            FOR_I3 dir[i] = arg->d3(i);
            NORMALIZE(dir);
          }
          else if (arg->getType() == D_DATA) {
            const int index = int(arg->d());
            if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;
            dir[index] = 1.0;
          }
          else return CTI_DATA_NOT_VALID;

          // loop faces and compute
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            double unit_u[3];
            const int icv      = zone_ptr->cvobf[ibf];
	    double u_minus_ubc[3];
	    FOR_I3 u_minus_ubc[i] = solver->u[icv][i] - u_bc[ibf][i];
	    const double mag_u = MAG(u_minus_ubc);

            if ( mag_u > 0.0) {

              for (int i = 0; i < 3; ++i)
                unit_u[i] = (solver->u[icv][i] - u_bc[ibf][i])/mag_u;

              tmp[ibf] = tau_wall[ibf]*DOT_PRODUCT(dir,unit_u);

            } else {

              tmp[ibf] = 0.0;

            }

          }
        }
        else
          return CTI_DATA_ARG_COUNT;
      }

      return CTI_DATA_OK;

    } else if ( name == yp_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          const int icv      = zone_ptr->cvobf[ibf];
          const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
          tmp[ibf]           = y1*sqrt(solver->rho[icv]*tau_wall[ibf])/solver->mu_lam[icv];
        }

      }

      return CTI_DATA_OK;

    }
  }

  return IdealGasBc::funcEvalCtiData(v,name,args,b_eval_func);
}

void WmAdiabatic::calcRhs(const double time, const int rk_stage) {

  // in order to avoid artificial correlation of the velocity
  // signal with the wall stress, the time filtering approach of
  // Yang et al is adopted here.  the filtered time scale is
  // relaxed based on an estimate of the outer units.  viscous units
  // could be an acceptable scaling here as well, but that would lead
  // to singular limits near separation points.
  // this results in the wall stress being computed once per step
  // and is performed inside of the calcRhs call for the first rk stage.

  if ( rk_stage == 1 ) {

    const double cp = solver->R_gas*solver->gamma/(solver->gamma-1.0);

    if ( T_wall) {

      // expecting that the wall temperature is the negative of its actual
      // value based on the convention for the adibatic solution...

      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
        T_wall[ibf] *= -1.0;

    }

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      const int icv = zone_ptr->cvobf[ibf];

      if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

        const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        //const double mag_u = MAG(solver->u[icv]);

        // need to project the normal component out of the velocity at this location ...

        double u_par[3];
        const double n_mag = MAG(zone_ptr->n_bf[ibf]);

        if ( n_mag >  0.0) {

          double unit_n[3];
          for (int i = 0; i < 3; ++i)
            unit_n[i] = zone_ptr->n_bf[ibf][i]/n_mag;

          const double un = DOT_PRODUCT(solver->u[icv],unit_n);

          for (int i = 0; i < 3; ++i)
            u_par[i] = solver->u[icv][i] - un*unit_n[i] - u_bc[ibf][i];

        } else {

          for (int i =0; i < 3; ++i)
            u_par[i] = solver->u[icv][i] - u_bc[ibf][i];

        }

        const double mag_u = MAG(u_par);


        // time filtering of the input velocity ...
        // \phi_new = eps*u + (1.0 -eps)*\phi_old
        // \phi_new - \phi_old = eps*(u - \phi_old)

        //const double alpha = 2.0;
        //const double eps   = mag_u*solver->dt / (alpha*y1 + solver->dt*mag_u);  // two times the local time scale
        //u1[ibf]           += eps* ( mag_u - u1[ibf]);

        // no time filtering ..

        u1[ibf]   = mag_u;

        // use the existing tau_wall to set a guess for the u_tau

        if ( wm_type == WM_ALGEBRAIC_ADIABATIC) {

          const double u_tau = sqrt(tau_wall[ibf] / solver->rho[icv]);
          tau_wall[ibf] = AlgebraicWM::solve_tau(u1[ibf], y1, solver->rho[icv],
                                                 solver->mu_lam[icv], u_tau);

        } else {

          assert( wm_type == WM_EQUILIBRIUM_ADIABATIC);

          double q_w = 0.0;
          EquilibriumWM::solve_qt(tau_wall[ibf],q_w,u1[ibf],solver->T[icv],solver->rho[icv],
                                  y1, solver->mu_lam[icv], solver->Pr_lam, cp, T_wall[ibf],
                                  relax, y_plus0, tol, true);

        }

      } else {

        //assert ( zone_ptr->area_bf[ibf] < 1.0e-16); // this should be a collapsed face.
        tau_wall[ibf] = 0.0;
	double u_minus_ubc[3];
	FOR_I3 u_minus_ubc[i] = solver->u[icv][i] - u_bc[ibf][i];
	u1[ibf] = MAG(u_minus_ubc);

      }
    }

    // on return if the wall temperature is available, we need to flip the sign so that it
    // can be queried ..

    if ( T_wall) {

      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
        T_wall[ibf] *= -1.0;

    }

  }
}

void WmAdiabatic::addBoundaryFlux(IdealGasRhs* rhs) const {

  assert( solver != NULL);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    // tau_wall has already been computed at this point, but we'll assume
    // that the wall stress is instantaneously aligned with the velocity vector
    // for the component that is parallel to the wall...

    double u_par_hat[3] = {0.0,0.0,0.0};
    double unit_n[3]    = {0.0,0.0,0.0};
    const double n_mag  = MAG(zone_ptr->n_bf[ibf]);

    if ( n_mag > 0.0) {

      for (int i = 0; i < 3; ++i)
        unit_n[i] = zone_ptr->n_bf[ibf][i]/n_mag;

    }

    const double un = DOT_PRODUCT(solver->u[icv],unit_n);

    for (int i = 0; i < 3; ++i)
      u_par_hat[i]  = solver->u[icv][i] - un*unit_n[i] - u_bc[ibf][i];

    const double u_mag = MAG(u_par_hat);

    double force[3];
    if ( u_mag > 0.0) {

      FOR_I3 force[i] = tau_wall[ibf]*u_par_hat[i]/u_mag*zone_ptr->area_bf[ibf];

      for (int i =0; i < 3;++i) {
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];
        rhs[icv].rhou[i] -= force[i];
      }
      rhs[icv].rhoE -= DOT_PRODUCT(force, u_bc[ibf]);

    } else {

      for (int i =0; i < 3; ++i)
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];

    }
  }
}

void WmAdiabatic::query(const string& param_str) const {

  // report an estimate of the local y^+ on the boundary and the integrated tau_wall

  double my_buf[3] = {0.0,0.0,0.0};
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

      const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
      const double u_tau = sqrt( tau_wall[ibf]/solver->rho[icv]);
      const double nu    = solver->mu_lam[icv] / solver->rho[icv];
      const double y_plus = y1*u_tau/nu;

      my_buf[0] += zone_ptr->area_bf[ibf];
      my_buf[1] += zone_ptr->area_bf[ibf]*tau_wall[ibf];
      my_buf[2] += zone_ptr->area_bf[ibf]*y_plus;

    }
  }

  double buf[3];
  MPI_Reduce(my_buf,buf,3,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if ( mpi_rank == 0 ) {

    cout << zone_ptr->getName() << " time, int tau_wall dA, avg y_plus = "
         << solver->time << "    "
         << buf[1] << "    "
         << buf[2]/buf[0] << endl;

  }
}

void WallIsothermal::initData() {
  Param* param = getParam(getName());
  T_bc         = param->getDouble(1);
}

void WallIsothermal::addBoundaryFlux(IdealGasRhs* rhs) const {

  assert( solver != NULL);

  const double cp = solver->R_gas*solver->gamma/(solver->gamma-1.0);
  assert( cp > 0.0);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv = zone_ptr->cvobf[ibf];

    //const double visc_coeff = (solver->mu_lam[icv] + solver->mu_sgs[icv]) *zone_ptr->area_over_delta_bf[ibf];
    //const double k_coeff    = (solver->loc_lam[icv]+ solver->loc_sgs[icv])*zone_ptr->area_over_delta_bf[ibf]*cp;

    const double visc_coeff = (solver->mu_lam[icv] ) *zone_ptr->area_over_delta_bf[ibf];
    const double k_coeff    = (solver->loc_lam[icv])*zone_ptr->area_over_delta_bf[ibf]*cp;
    assert( k_coeff >= 0.0);


    for (int i =0; i < 3; ++i)
      rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i] + visc_coeff*(solver->u[icv][i] - u_bc[i]);

    // XXX if u_bc \neq 0, there should be a u \dot \tau term added here..
    rhs[icv].rhoE -= k_coeff*(solver->T[icv] - T_bc);
  }

}

void WallIsothermal::query(const string& param_str) const {

  double my_buf[4] = {0.0,0.0,0.0,0.0};
  const double cp = solver->R_gas*solver->gamma/(solver->gamma-1.0);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

      const double y1         = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
      const double visc_coeff = solver->mu_lam[icv]* zone_ptr->area_over_delta_bf[ibf];
      double du_mag = 0.0;
      for (int i =0; i < 3; ++i)
        du_mag  += (solver->u[icv][i] - u_bc[i])*(solver->u[icv][i] - u_bc[i]);

      const double tau_wall   = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
      const double y_plus      = y1*sqrt(tau_wall*solver->rho[icv])/solver->mu_lam[icv];

      const double k_coeff    = solver->loc_lam[icv]*zone_ptr->area_over_delta_bf[ibf]*cp;
      const double q_wall      = k_coeff*abs(solver->T[icv] - T_bc)/zone_ptr->area_bf[ibf];

      my_buf[0]              += zone_ptr->area_bf[ibf];
      my_buf[1]              += zone_ptr->area_bf[ibf]*tau_wall;
      my_buf[2]              += zone_ptr->area_bf[ibf]*y_plus;
      my_buf[3]              += zone_ptr->area_bf[ibf]*q_wall;

    }
  }

  double buf[4];
  MPI_Reduce(my_buf,buf,4,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if ( mpi_rank == 0 ) {

    cout << " > " << zone_ptr->getName() << " time, int(tau_wall)dA, int(q_wall)dA, avg y_plus = "
         << solver->time << "     "
         << buf[1] << "     "
         << buf[3] << "     "
         << buf[2]/buf[0] << endl;

  }
}

CtiRegister::CtiDataError WallIsothermal::funcEvalCtiData(CtiRegister::CtiData& v,
                                                          const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if ( zone_ptr->isLeadingStrMatched(name) ) {

    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix.  we have to check for the actual data now..

    const string tau_wall_str = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str  = zone_ptr->getName() + ":" + "y_plus";
    const string qw_str  = zone_ptr->getName() + ":" + "q_wall";

    // all of the following functions do not take any arguments -- check first


    if ( name == tau_wall_str) {

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {
        if ( args.size() == 0) {
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            const int icv           = zone_ptr->cvobf[ibf];
            const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
            double du_mag = 0.0;
            for (int i =0; i < 3; ++i)
              du_mag  += (solver->u[icv][i] - u_bc[i])*(solver->u[icv][i] - u_bc[i]);

            tmp[ibf]                = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
          }
        }
        else if (args.size() == 1) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          double dir[3] = {0.0,0.0,0.0};

          if (arg->getType() == D3_DATA) {
            FOR_I3 dir[i] = arg->d3(i);
            NORMALIZE(dir);
          }
          else if (arg->getType() == D_DATA) {
            const int index = int(arg->d());
            if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;
            dir[index] = 1.0;
          }
          else return CTI_DATA_NOT_VALID;

          // loop faces and compute
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            const int icv           = zone_ptr->cvobf[ibf];
            const double visc_coeff = solver->mu_lam[icv] *zone_ptr->area_over_delta_bf[ibf];

            double local_u[3] = DIFF(solver->u[icv],u_bc);
            tmp[ibf]          = visc_coeff*DOT_PRODUCT(dir,local_u)/zone_ptr->area_bf[ibf];
          }

        }
        else return CTI_DATA_ARG_COUNT;
      }

      return CTI_DATA_OK;

    } else if ( name == yp_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          const int icv           = zone_ptr->cvobf[ibf];
          const double y1         = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
          const double visc_coeff = solver->mu_lam[icv]/y1;

          double du_mag = 0.0;
          for (int i =0; i < 3; ++i)
            du_mag  += (solver->u[icv][i] - u_bc[i])*(solver->u[icv][i] - u_bc[i]);

          const double tau_wall      = visc_coeff*sqrt(du_mag);
          tmp[ibf]                = y1* sqrt(tau_wall*solver->rho[icv])/solver->mu_lam[icv];

        }

      }

      return CTI_DATA_OK;

    } else if ( name == qw_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {

        const double cp = solver->R_gas*solver->gamma/(solver->gamma-1.0);

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          const int icv        = zone_ptr->cvobf[ibf];
          const double k_coeff = solver->loc_lam[icv]*zone_ptr->area_over_delta_bf[ibf]*cp;

          assert( k_coeff >= 0.0);

          tmp[ibf]             = k_coeff*(solver->T[icv] - T_bc)/zone_ptr->area_bf[ibf];
        }

      }

      return CTI_DATA_OK;

    }
  }

  return IdealGasBc::funcEvalCtiData(v,name,args,b_eval_func);
}

void WallCht::initData() {

  if (solver->cht == NULL) {
    CERR("WALL_CHT bcs require CHT");
  }

  int * ssp_flag = new int[solver->subSurface->nsp];
  for (int issp = 0; issp < solver->subSurface->nsp; ++issp)
    ssp_flag[issp] = -1;

  for (int ibf = zone_ptr->ibf_f; ibf <= zone_ptr->ibf_l; ++ibf) {
    for (int sob = solver->sspobf_i[ibf]; sob != solver->sspobf_i[ibf+1]; ++sob) {
      const int issp = solver->sspobf_v[sob];
      ssp_flag[issp] = 0;
    }
  }

  np = 0;
  for (int issp = 0; issp < solver->subSurface->nsp; ++issp)
    if (ssp_flag[issp] == 0)
      ssp_flag[issp] = np++;

  assert(Tp == NULL);
  Tp = new double[np];
  for (int ip = 0; ip < np; ++ip) Tp[ip] = HUGE_VAL; // put garbage here to ensure it is set before being used

  assert(qp == NULL);
  qp = new double[np];
  for (int ip = 0; ip < np; ++ip) qp[ip] = HUGE_VAL; // put garbage here to ensure it is set before being used

  double (*xp)[3] = new double[np][3];
  int ip = 0;
  for (int issp = 0; issp < solver->subSurface->nsp; ++issp) {
    if (ssp_flag[issp] >= 0) {
      FOR_I3 xp[ip][i] = solver->subSurface->xp[issp][i];
      const int bits = int(solver->subSurface->isp_global_and_bits[issp]>>52);
      if (bits)
        PeriodicData::periodicTranslate(xp[ip],1,BitUtils::flipPeriodicBits(bits));
      ++ip;
    }
  }
  assert(ip == np);

  if ( mpi_rank == 0)
    cout << " > ... " << zone_ptr->getName() << "  about to register cht data " << endl;

  //solver->cht->registerBoundaryData(xp,Tp,qp,np,zone_ptr->index);
  solver->cht->registerBoundaryData(xp,Tp,qp,ssp_flag,np,solver->subSurface->nsp,zone_ptr->index);
  delete[] xp; // does not need to be persistent.

  // build the ipobf_i/v/wgt needed to move information to/from the Tp/qp arrays...
  assert(ipobf_i == NULL);
  assert(zone_ptr->ibf_l-zone_ptr->ibf_f+1 == zone_ptr->nbf);
  ipobf_i = new int[zone_ptr->nbf+1];
  ipobf_i[0] = 0;
  const int ipobf_s = solver->sspobf_i[zone_ptr->ibf_l+1] - solver->sspobf_i[zone_ptr->ibf_f];
  assert(ipobf_v == NULL);
  ipobf_v = new int[ipobf_s];
  assert(ipobf_wgt == NULL);
  ipobf_wgt = new double[ipobf_s];
  int iob = 0;
  for (int ibf = zone_ptr->ibf_f; ibf <= zone_ptr->ibf_l; ++ibf) {
    double wgt_sum = 0.0;
    const int iob0 = iob;
    bool b_unmatched = false;
    for (int sob = solver->sspobf_i[ibf]; sob != solver->sspobf_i[ibf+1]; ++sob) {
      const int issp = solver->sspobf_v[sob];
      //assert((ssp_flag[issp] >= 0)&&(ssp_flag[issp] < np));
      if (ssp_flag[issp] >= 0) {
        assert(ssp_flag[issp] < np);
        ipobf_v[iob] = ssp_flag[issp];
        ipobf_wgt[iob] = solver->sspobf_wgt[sob];
        wgt_sum += ipobf_wgt[iob];
        ++iob;
      }
      else {
        b_unmatched = true;
      }
    }
    //assert(fabs(1.0-wgt_sum) < 1.0E-12);
    ipobf_i[ibf-zone_ptr->ibf_f+1] = iob;
    if (b_unmatched) {
      // redistribute unmatched ip weights...
      assert(wgt_sum > 1.0E-12);
      const double inv_wgt_sum = 1.0/wgt_sum;
      for (int iob_ = iob0; iob_ < iob; ++iob_)
        ipobf_wgt[iob_] *= inv_wgt_sum;
    }
  }
  //assert(iob == ipobf_s);
  assert(iob <= ipobf_s);
  delete[] ssp_flag;

}

void WallCht::addBoundaryFlux(IdealGasRhs* rhs) const {

  //if (mpi_rank == 0) cout << "in WallCht::addBoundaryFlux for zone: " << getName() << endl;

  // confirm that the Tp data is properly updated from the fem solver...

  /*
  int * ssp_flag = new int[solver->subSurface->nsp];
  for (int issp = 0; issp < solver->subSurface->nsp; ++issp)
    ssp_flag[issp] = -1;
  for (int ibf = zone_ptr->ibf_f; ibf <= zone_ptr->ibf_l; ++ibf) {
    for (int sob = solver->sspobf_i[ibf]; sob != solver->sspobf_i[ibf+1]; ++sob) {
      const int issp = solver->sspobf_v[sob];
      ssp_flag[issp] = 0;
    }
  }
  double (*xp)[3] = new double[np][3];
  int ip = 0;
  for (int issp = 0; issp < solver->subSurface->nsp; ++issp) {
    if (ssp_flag[issp] == 0) {
      FOR_I3 xp[ip][i] = solver->subSurface->xp[issp][i];
      ++ip;
    }
  }
  assert(ip == np);
  for (ip = 0; ip < np; ++ip) {
    const double Tp_check = xp[ip][0] + 1.1*xp[ip][1] + 1.2*xp[ip][2];
    assert(fabs(Tp[ip]-Tp_check) < 1.0E-12);
  }
  MPI_Pause("YOU ROCK!");
  */

  assert( solver != NULL);
  const double cp = solver->R_gas*solver->gamma/(solver->gamma-1.0);

  // zero the qp so it can be summed below...
  for (int ip = 0; ip < np; ++ip) qp[ip] = 0.0;

  // loop on bf's...
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    // reconstruct the boundary temperature from the Tp data that should
    // be current in the Tp[ip] array. (updated with an earlier call to
    // cht->updateChtTemperatureToFlowBcs())...
    double T_bc = 0.0;
    for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
      const int ip = ipobf_v[iob];
      const double wgt = ipobf_wgt[iob];
      T_bc += wgt*Tp[ip];
    }

    const int icv = zone_ptr->cvobf[ibf];
    const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
    const double k_coeff    = solver->loc_lam[icv]*zone_ptr->area_over_delta_bf[ibf]*cp;

    for (int i =0; i < 3; ++i)
      rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i] + visc_coeff*(solver->u[icv][i] - 0.0); // u_bc == 0 for CHT

    // XXX if u_bc \neq 0, there should be a u \dot \tau term added here..
    rhs[icv].rhoE -= k_coeff*(solver->T[icv] - T_bc);

    // scatter/sum q_bf to the qp[ip] data points for exchange to the FEM solve this sub-step...
    for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
      const int ip = ipobf_v[iob];
      const double wgt = ipobf_wgt[iob];
      qp[ip] += wgt*k_coeff*(solver->T[icv] - Tp[ip]); // should match -q_bf when summed over ip
    }

  }

}

void WallCht::query(const string& param_str) const {

  double my_buf[5] = {0.0,0.0,0.0,0.0,0.0};
  const double cp = solver->R_gas*solver->gamma/(solver->gamma-1.0);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    // compute this...
    double T_bc = 0.0;
    for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
      const int ip = ipobf_v[iob];
      const double wgt = ipobf_wgt[iob];
      T_bc += wgt*Tp[ip];
    }

    if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

      const double y1         = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
      const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
      double du_mag = 0.0;
      for (int i =0; i < 3; ++i)
        du_mag  += (solver->u[icv][i] - 0.0)*(solver->u[icv][i] - 0.0); // u_bc == 0.0

      const double tau_wall   = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
      const double y_plus      = y1*sqrt(tau_wall*solver->rho[icv])/solver->mu_lam[icv];

      const double k_coeff    = solver->loc_lam[icv]*zone_ptr->area_over_delta_bf[ibf]*cp;
      const double q_wall      = k_coeff*abs(solver->T[icv] - T_bc)/zone_ptr->area_bf[ibf];

      my_buf[0]              += zone_ptr->area_bf[ibf];
      my_buf[1]              += zone_ptr->area_bf[ibf]*tau_wall;
      my_buf[2]              += zone_ptr->area_bf[ibf]*y_plus;
      my_buf[3]              += zone_ptr->area_bf[ibf]*q_wall;
      my_buf[4]              += zone_ptr->area_bf[ibf]*T_bc;

    }
  }

  double buf[5];
  MPI_Reduce(my_buf,buf,5,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if ( mpi_rank == 0 ) {
    cout << zone_ptr->getName() << " time, int(tau_wall)dA, int(q_wall)dA, avg y_plus, avg T_bc = "
         << solver->time << " "
         << buf[1] << " "
         << buf[3] << " "
         << buf[2]/buf[0] << " "
         << buf[4]/buf[0] << endl;
  }
}

CtiRegister::CtiDataError WallCht::funcEvalCtiData(CtiRegister::CtiData& v,
                                                   const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if ( zone_ptr->isLeadingStrMatched(name) ) {

    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix.  we have to check for the actual data now..

    const string tau_wall_str = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str  = zone_ptr->getName() + ":" + "y_plus";
    const string qw_str  = zone_ptr->getName() + ":" + "q_wall";
    const string t_wall_str = zone_ptr->getName() + ":" + "T_wall";

    // all of the following functions do not take any arguments -- check first

    if ( name == tau_wall_str) {

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {
        if ( args.size() == 0) {
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            const int icv           = zone_ptr->cvobf[ibf];
            const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
            double du_mag = 0.0;
            for (int i =0; i < 3; ++i)
              du_mag  += (solver->u[icv][i] - 0.0)*(solver->u[icv][i] - 0.0);

            tmp[ibf]                = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
          }
        }
        else if (args.size() == 1) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          double dir[3] = {0.0,0.0,0.0};

          if (arg->getType() == D3_DATA) {
            FOR_I3 dir[i] = arg->d3(i);
            NORMALIZE(dir);
          }
          else if (arg->getType() == D_DATA) {
            const int index = int(arg->d());
            if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;
            dir[index] = 1.0;
          }
          else return CTI_DATA_NOT_VALID;

          // loop faces and compute
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            const int icv           = zone_ptr->cvobf[ibf];
            const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];

            tmp[ibf]          = visc_coeff*DOT_PRODUCT(dir,solver->u[icv])/zone_ptr->area_bf[ibf];
          }

        }
        else return CTI_DATA_ARG_COUNT;
      }

      return CTI_DATA_OK;

    } else if ( name == yp_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          const int icv           = zone_ptr->cvobf[ibf];
          const double y1         = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
          const double visc_coeff = solver->mu_lam[icv]/y1;

          double du_mag = 0.0;
          for (int i =0; i < 3; ++i)
            du_mag  += (solver->u[icv][i] - 0.0)*(solver->u[icv][i] - 0.0);

          const double tau_wall      = visc_coeff*sqrt(du_mag);
          tmp[ibf]                = y1* sqrt(tau_wall*solver->rho[icv])/solver->mu_lam[icv];

        }

      }

      return CTI_DATA_OK;

    } else if ( name == qw_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {

        const double cp = solver->R_gas*solver->gamma/(solver->gamma-1.0);

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          double T_bc = 0.0;
          for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
            const int ip = ipobf_v[iob];
            const double wgt = ipobf_wgt[iob];
            T_bc += wgt*Tp[ip];
          }

          const int icv        = zone_ptr->cvobf[ibf];
          const double k_coeff = solver->loc_lam[icv]*zone_ptr->area_over_delta_bf[ibf]*cp;

          tmp[ibf]             = k_coeff*abs(solver->T[icv] - T_bc)/zone_ptr->area_bf[ibf];
        }

      }

      return CTI_DATA_OK;

    }
    else if (name == t_wall_str) {

      if (args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          // reconstruct wall temperature from Tp...
          double T_wall = 0.0;
          for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
            const int ip = ipobf_v[iob];
            const double wgt = ipobf_wgt[iob];
            T_wall += wgt*Tp[ip];
          }
          tmp[ibf] = T_wall;
        }
      }

      return CTI_DATA_OK;

    }
  }

  return IdealGasBc::funcEvalCtiData(v,name,args,b_eval_func);
}

void WmIsothermal::initData() {

  assert( tau_wall  == NULL); tau_wall  = new double[zone_ptr->nbf];
  assert( u1     == NULL); u1     = new double[zone_ptr->nbf];
  assert( q_wall == NULL); q_wall = new double[zone_ptr->nbf];

  /*
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    tau_wall[ibf]  = 0.0;
    q_wall[ibf] = 0.0;
  }
  */

  Param * param = getParam(zone_ptr->getName());
  T_bc          = param->getDouble(1);

  int iarg = 2;
  while ( iarg < param->size()) {

    string token = param->getString(iarg++);
    if ( token == "RELAX") {
      relax = param->getDouble(iarg++);
    }
    else if ( token == "Y_PLUS0") {
      y_plus0 = param->getDouble(iarg++);
    }
    else if ( token == "TOL") {
      tol = param->getDouble(iarg++);
    }

  }

}

void WmIsothermal::initialHook() {

  if ( !zone_ptr->checkDataFlag("u1")) {

    const double cp = solver->R_gas*solver->gamma/(solver->gamma-1.0);

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      const int icv = zone_ptr->cvobf[ibf];
      u1[ibf] = MAG(solver->u[icv]);

      // set an initial value for tau_wall, this will constitute the
      // initial guess for the solution when it is finally constructed.

      if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

        const double y1 = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        tau_wall[ibf]      = solver->mu_lam[icv] * MAG(solver->u[icv])/y1;
        q_wall[ibf]        = solver->mu_lam[icv] * (cp*solver->T[icv] + 0.5*DOT_PRODUCT(solver->u[icv],solver->u[icv]) - cp*T_bc)/y1;
      }
      else {
        tau_wall[ibf] = 0.0;
        q_wall[ibf]   = 0.0;
      }
    }
  }
}

int WmIsothermal::storeData(double * buf) const {

  // the wall model needs to store tau_wall because it is used as
  // an initial guess. No need to store u1. It is not used...
  if ( wm_type == WM_ALGEBRAIC_ISOTHERMAL) {
    if (buf) {
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
        buf[ibf] = tau_wall[ibf];
      }
    }

    return zone_ptr->nbf;
  }

  assert(wm_type == WM_EQUILIBRIUM_ISOTHERMAL);
  return 0;  // EQ model computes tau_wall
}

int WmIsothermal::restoreData(double * buf) {

  if ( wm_type == WM_ALGEBRAIC_ISOTHERMAL) {
    assert(buf);
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      tau_wall[ibf] = buf[ibf];
    }

    return zone_ptr->nbf;
  }

  return 0;
}

WmIsothermal::~WmIsothermal() {

  DELETE(u1);
  DELETE(tau_wall);
  DELETE(q_wall);

}

void WmIsothermal::calcRhs(const double time, const int rk_stage) {

  // see notes left in IdealGasSolverBcs.cpp for the time
  // filtering approach

  if ( rk_stage == 1 ) {

    const double cp = solver->R_gas*solver->gamma/(solver->gamma-1.0);

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      const int icv = zone_ptr->cvobf[ibf];

      if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

        const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        const double mag_u = MAG(solver->u[icv]);
        //const double eps   = mag_u*solver->dt / (y1 + solver->dt*mag_u);

        // \phi_new = eps*u + (1.0 -eps)*\phi_old
        // \phi_new - \phi_old = eps*(u - \phi_old)

        //u1[ibf]   += eps* ( mag_u - u1[ibf]);
        u1[ibf]   = mag_u;

        // use the existing tau_wall to set a guess for the u_tau

        if ( wm_type == WM_ALGEBRAIC_ISOTHERMAL) {

          const double u_tau = sqrt(tau_wall[ibf] / solver->rho[icv]);
          tau_wall[ibf]      = AlgebraicWM::solve_tau(u1[ibf], y1, solver->rho[icv],
                                                      solver->mu_lam[icv], u_tau);

          // approximate material properties from the interior -- these are low Ma approximations
          // that ignore property changes in the near wall cell due to compressibility or strong heating

          q_wall[ibf] = AlgebraicWM::compute_q_wall_approx(solver->T[icv], T_bc, tau_wall[ibf],
                                                           solver->rho[icv], solver->mu_lam[icv],
                                                           solver->loc_lam[icv], cp, y1);

        } else {

          assert( wm_type == WM_EQUILIBRIUM_ISOTHERMAL);

          EquilibriumWM::solve_qt(tau_wall[ibf], q_wall[ibf], u1[ibf], solver->T[icv], solver->rho[icv],
                                  y1, solver->mu_lam[icv], solver->Pr_lam,cp, T_bc,
                                  relax, y_plus0, tol, true);

        }

      } else {

        //assert ( zone_ptr->area_bf[ibf] < 1.0e-16); // this should be a collapsed face.
        tau_wall[ibf]  = 0.0;
        q_wall[ibf] = 0.0;
        u1[ibf]     = MAG(solver->u[icv]);

      }
    }
  }
}

void WmIsothermal::addBoundaryFlux(IdealGasRhs* rhs) const {

  assert( solver != NULL);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    // tau_wall has already been computed at this point, but we'll assume
    // that the wall stress is instantaneously aligned with the velocity vector

    double u_mag = 0.0;
    for (int i =0; i < 3; ++i)
      u_mag += solver->u[icv][i] * solver->u[icv][i];
    u_mag = sqrt( max(0.0, u_mag));

    if ( u_mag > 0.0) {

      for (int i =0; i < 3;++i) {
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];
        rhs[icv].rhou[i] -= tau_wall[ibf]*solver->u[icv][i]/u_mag*zone_ptr->area_bf[ibf];
      }
      // XXX if u_bc \neq 0, there should be a u \dot \tau term added here..
      rhs[icv].rhoE -= q_wall[ibf]*zone_ptr->area_bf[ibf];

    } else {

      for (int i =0; i < 3; ++i)
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];

      // XXX if u_bc \neq 0, there should be a u \dot \tau term added here..
      rhs[icv].rhoE -= q_wall[ibf]*zone_ptr->area_bf[ibf];

    }
  }
}

CtiRegister::CtiDataError WmIsothermal::funcEvalCtiData(CtiRegister::CtiData& v,
                                                        const string& name,
                                                        list<CtiRegister::CtiData>& args,
                                                        const bool b_eval_func) {
  using namespace CtiRegister;

  if ( zone_ptr->isLeadingStrMatched(name) ) {

    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix.  we have to check for the actual data now..

    // for consistency, we are supporting a tau_wall(), y_plus()
    // function to have the same syntax as all of the other wall zones..

    const string tau_wall_str = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str  = zone_ptr->getName() + ":" + "y_plus";
    const string qw_str  = zone_ptr->getName() + ":" + "q_wall";

    // all of the following functions do not take any arguments -- check first


    if ( name == tau_wall_str) {

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {
        if (args.size() == 0) {
          // use what was computed by wall model
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
            tmp[ibf] = tau_wall[ibf];
          }
        }
        else if (args.size() == 1) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          double dir[3] = {0.0,0.0,0.0};

          if (arg->getType() == D3_DATA) {
            FOR_I3 dir[i] = arg->d3(i);
            NORMALIZE(dir);
          }
          else if (arg->getType() == D_DATA) {
            const int index = int(arg->d());
            if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;
            dir[index] = 1.0;
          }
          else return CTI_DATA_NOT_VALID;

          // loop faces and compute
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            double unit_u[3];
            const int icv      = zone_ptr->cvobf[ibf];
            const double mag_u = MAG(solver->u[icv]);

            if ( mag_u > 0.0) {

              for (int i = 0; i < 3; ++i)
                unit_u[i] = solver->u[icv][i]/mag_u;

              tmp[ibf] = tau_wall[ibf]*DOT_PRODUCT(dir,unit_u);

            } else {

              tmp[ibf] = 0.0;

            }

          }

        }
        else
          return CTI_DATA_ARG_COUNT;
      }

      return CTI_DATA_OK;

    } else if ( name == yp_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

          const int icv      = zone_ptr->cvobf[ibf];
          const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
          tmp[ibf]           = y1*sqrt(solver->rho[icv]*tau_wall[ibf])/solver->mu_lam[icv];
        }

      }

      return CTI_DATA_OK;

    } else if ( name == qw_str) {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {

        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
          tmp[ibf] = q_wall[ibf];

      }

      return CTI_DATA_OK;

    }
  }

  return IdealGasBc::funcEvalCtiData(v,name,args,b_eval_func);
}

void WmIsothermal::query(const string& param_str) const {

  // report an estimate of the local y^+ on the boundary and the integrated tau_wall

  double my_buf[4] = {0.0,0.0,0.0,0.0};
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

      const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
      const double u_tau = sqrt( tau_wall[ibf]/solver->rho[icv]);
      const double nu    = solver->mu_lam[icv] / solver->rho[icv];
      const double y_plus = y1*u_tau/nu;

      my_buf[0] += zone_ptr->area_bf[ibf];
      my_buf[1] += zone_ptr->area_bf[ibf]*tau_wall[ibf];
      my_buf[2] += zone_ptr->area_bf[ibf]*y_plus;
      my_buf[3] += zone_ptr->area_bf[ibf]*q_wall[ibf];

    }
  }

  double buf[4];
  MPI_Reduce(my_buf,buf,4,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if ( mpi_rank == 0 ) {

    cout << zone_ptr->getName() << " time, int tau_wall dA, int q_w dA, avg y_plus = "
         << solver->time << "    "
         << buf[1] << "    "
         << buf[3] << "    "
         << buf[2]/buf[0] << endl;

  }
}

void WmAlgConductive::initData() {
  Param * param  = getParam(getName());
  int i_param     = 1;
  T_end_bc        = param->getDouble(i_param++);

  // Read layer lengths and conductivities
  thermal_resistance_external_times_area = 0.0;
  while (param->getString(i_param) != "END") {
    layer_length_bc.push_back(param->getDouble(i_param++)); // [m]
    thermal_conductivity_bc.push_back(param->getDouble(i_param++)); // [W/m*K]
    thermal_resistance_external_times_area += layer_length_bc.back() / thermal_conductivity_bc.back(); // [m] / [W/m*K] = [W*m^2/K]
  }

  // Layer length and conductivity vectors are the same size!
  assert(layer_length_bc.size() == thermal_conductivity_bc.size());

  assert(tau_wall == NULL); tau_wall  = new double[zone_ptr->nbf];
  assert(u1       == NULL); u1        = new double[zone_ptr->nbf];
  assert(q_wall   == NULL); q_wall    = new double[zone_ptr->nbf];
  assert(T_bc     == NULL); T_bc      = new double[zone_ptr->nbf];
}

void WmAlgConductive::initialHook() {
  if (!zone_ptr->checkDataFlag("u1")) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ ibf) {
      const int icv = zone_ptr->cvobf[ibf];
      u1[ibf]       = MAG(solver->u[icv]);

      // set an initial value for tau_wall, this will constitute the
      // initial guess for the solution when it is finally constructed

      if (zone_ptr->area_over_delta_bf[ibf] > 0.0) {
        const double y1 = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        tau_wall[ibf]   = solver->mu_lam[icv] * MAG(solver->u[icv])/y1;
        q_wall[ibf]     = 0.0;
        T_bc[ibf]       = 0.0;
      }
      else {
        tau_wall[ibf]   = 0.0;
        q_wall[ibf]     = 0.0;
        T_bc[ibf]       = 0.0;
      }
    }
  }
}

int WmAlgConductive::storeData(double * buf) const {

  // the wall model needs to store tau_wall because it is used as
  // an initial guess. No need to store u1. It is not used...

  if (buf) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      buf[ibf] = tau_wall[ibf];
    }
  }

  return zone_ptr->nbf;

}

int WmAlgConductive::restoreData(double * buf) {

  assert(buf);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    tau_wall[ibf] = buf[ibf];
  }

  return zone_ptr->nbf;

}

WmAlgConductive::~WmAlgConductive() {
  DELETE(u1);
  DELETE(tau_wall);
  DELETE(q_wall);
  DELETE(T_bc);
}

void WmAlgConductive::calcRhs(const double time, const int rk_stage) {
  // see notes left in IdealGasSolverBcs.cpp for the time filtering approach
  bool laminar_q_wall = false;
  if (rk_stage == 1) {
    const double cp = solver->R_gas*solver->gamma/(solver->gamma - 1.0);

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv = zone_ptr->cvobf[ibf];
      // Ensure that the cell actually has area
      if (zone_ptr->area_over_delta_bf[ibf] > 0.0) {
        const double y1     = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        const double mag_u  = MAG(solver->u[icv]);

        u1[ibf]             = mag_u;

        // use the existing tau_wall to set a guess for the u_tau
        const double u_tau = sqrt(tau_wall[ibf] / solver->rho[icv]);
        tau_wall[ibf] = AlgebraicWM::solve_tau(u1[ibf], y1, solver->rho[icv], solver->mu_lam[icv], u_tau);

        if (laminar_q_wall) {
          const double L_cv = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf]; // [m]
          const double k_cv = solver->loc_lam[icv]*cp; // gas thermal conductivity [W/m*K]
          const double thermal_resistance_cv = L_cv/k_cv/zone_ptr->area_bf[ibf]; // [m] / [W/m*K] / [m^2] = [K/W]
          const double thermal_resistance_overall = thermal_resistance_external_times_area/zone_ptr->area_bf[ibf] + thermal_resistance_cv; // [K/W]
          q_wall[ibf] = (solver->T[icv] - T_end_bc)/thermal_resistance_overall/zone_ptr->area_bf[ibf]; // [K]/[K/W]/[m^2] = [W/m^2]
          T_bc[ibf]   = solver->T[icv] - q_wall[ibf]*thermal_resistance_cv*zone_ptr->area_bf[ibf]; // [K]  - [W/m^2]*[K/W]*[m^2] = [K] - [K]
        }
        else {
          const double thermal_resistance_cv = AlgebraicWM::compute_R_approx(tau_wall[ibf], solver->rho[icv], solver->mu_lam[icv],
          solver->loc_lam[icv], cp, y1);
          const double thermal_resistance_overall = thermal_resistance_external_times_area/zone_ptr->area_bf[ibf] + thermal_resistance_cv; // [K/W]
          q_wall[ibf] = (solver->T[icv] - T_end_bc)/thermal_resistance_overall/zone_ptr->area_bf[ibf]; // [K]/[K/W]/[m^2] = [W/m^2]
          T_bc[ibf]   = solver->T[icv] - q_wall[ibf]*thermal_resistance_cv*zone_ptr->area_bf[ibf];// [K] - [W/m^2]*[K/W]*[m^2] = [K] - [K]
        }
      }
    }
  }
}

void WmAlgConductive::addBoundaryFlux(IdealGasRhs* rhs) const {
  assert( solver != NULL);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv = zone_ptr->cvobf[ibf];
    // tau_wall has already been computed at this point, but we'll assume
    // that the wall stress is instantaneously aligned with the velocity vector
    double u_mag = 0.0;
    for (int i =0; i < 3; ++i)
      u_mag += solver->u[icv][i] * solver->u[icv][i];
    u_mag = sqrt( max(0.0, u_mag));

    if ( u_mag > 0.0) {
      for (int i =0; i < 3;++i) {
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];
        rhs[icv].rhou[i] -= tau_wall[ibf]*solver->u[icv][i]/u_mag*zone_ptr->area_bf[ibf];
      }
      // XXX if u_bc \neq 0, there should be a u \dot \tau term added here..
      rhs[icv].rhoE -= q_wall[ibf]*zone_ptr->area_bf[ibf]; // [W]
    }
    else {
      for (int i =0; i < 3; ++i)
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];
      // XXX if u_bc \neq 0, there should be a u \dot \tau term added here..
      rhs[icv].rhoE -= q_wall[ibf]*zone_ptr->area_bf[ibf]; // [W/m^2]*[m^2] = [W]
    }
  }
}

CtiRegister::CtiDataError WmAlgConductive::funcEvalCtiData(CtiRegister::CtiData& v,
    const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if (zone_ptr->isLeadingStrMatched(name)) {
    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix.  we have to check for the actual data now..

    // for consistency, we are supporting a tau_wall(), y_plus()
    // function to have the same syntax as all of the other wall zones..
    const string tau_wall_str = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str  = zone_ptr->getName() + ":" + "y_plus";
    const string qw_str  = zone_ptr->getName() + ":" + "q_wall";
    const string Tbc_str = zone_ptr->getName() + ":" + "T_bc";

    // all of the following functions do not take any arguments -- check first
    if (name == tau_wall_str) {
      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {
        if (args.size() == 0) {
          // use what was computed by wall model
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
            tmp[ibf] = tau_wall[ibf];
          }
        }
        else if (args.size() == 1) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          double dir[3] = {0.0,0.0,0.0};

          if (arg->getType() == D3_DATA) {
            FOR_I3 dir[i] = arg->d3(i);
            NORMALIZE(dir);
          }
          else if (arg->getType() == D_DATA) {
            const int index = int(arg->d());
            if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;
            dir[index] = 1.0;
          }
          else return CTI_DATA_NOT_VALID;

          // loop faces and compute
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
            const int icv      = zone_ptr->cvobf[ibf];
            const double mag_u = MAG(solver->u[icv]);

            double unit_u[3];
            if ( mag_u > 0.0) {
              for (int i = 0; i < 3; ++i) unit_u[i] = solver->u[icv][i]/mag_u;

              tmp[ibf] = tau_wall[ibf]*DOT_PRODUCT(dir,unit_u);
            }
            else {
              tmp[ibf] = 0.0;
            }
          }
        }
        else return CTI_DATA_ARG_COUNT;
      }

      return CTI_DATA_OK;
    }

    else if ( name == yp_str) {
      if ( args.size() != 0)  return CTI_DATA_ARG_COUNT;
      double * tmp = zone_ptr->createBfD1Data(v);
      if ( b_eval_func) {
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          const int icv      = zone_ptr->cvobf[ibf];
          const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
          tmp[ibf]           = y1*sqrt(solver->rho[icv]*tau_wall[ibf])/solver->mu_lam[icv];
        }
      }
      return CTI_DATA_OK;
    }

    else if ( name == qw_str) {
      if ( args.size() != 0) return CTI_DATA_ARG_COUNT;
      double * tmp = zone_ptr->createBfD1Data(v);
      if ( b_eval_func) {
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
          tmp[ibf] = q_wall[ibf];
      }
      return CTI_DATA_OK;
    }

    else if (name == Tbc_str) {
      if (args.size() != 0) return CTI_DATA_ARG_COUNT;
      double * tmp = zone_ptr->createBfD1Data(v);
      if (b_eval_func) {
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
          tmp[ibf] = T_bc[ibf];
      }
      return CTI_DATA_OK;
    }
  }

  return IdealGasBc::funcEvalCtiData(v,name,args,b_eval_func);
}

void WmAlgConductive::query(const string& param_str) const {

  // report an estimate of the local y^+ on the boundary and the integrated tau_wall

  double my_buf[5] = {0.0,0.0,0.0,0.0,0.0};
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {

      const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
      const double u_tau = sqrt( tau_wall[ibf]/solver->rho[icv]);
      const double nu    = solver->mu_lam[icv] / solver->rho[icv];
      const double y_plus = y1*u_tau/nu;

      my_buf[0] += zone_ptr->area_bf[ibf];
      my_buf[1] += zone_ptr->area_bf[ibf]*tau_wall[ibf];
      my_buf[2] += zone_ptr->area_bf[ibf]*y_plus;
      my_buf[3] += zone_ptr->area_bf[ibf]*q_wall[ibf];
      my_buf[4] += zone_ptr->area_bf[ibf]*T_bc[ibf];

    }
  }

  double buf[5];
  MPI_Reduce(my_buf,buf,5,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if ( mpi_rank == 0 ) {

    cout << zone_ptr->getName() << " time, int tau_wall dA, int q_w dA, avg y_plus, avg T_bc = "
         << solver->time  << "    "
         << buf[1]        << "    "
         << buf[3]        << "    "
         << buf[2]/buf[0] << "    "
         << buf[4]/buf[0] << endl;
  }
}


void WallConductive::initData() {
  Param* param = getParam(getName());

  int i_param = 1;
  T_end_bc = param->getDouble(i_param++);

  while (param->getString(i_param) != "END") {
    layer_length_bc.push_back(param->getDouble(i_param++));
    thermal_conductivity_bc.push_back(param->getDouble(i_param++));
    // XXX Assuming layers are not parallel.
    thermal_resistance_external_times_area += layer_length_bc.back() / thermal_conductivity_bc.back();
  }

  assert(layer_length_bc.size() == thermal_conductivity_bc.size());
}

void WallConductive::addBoundaryFlux(IdealGasRhs* rhs) const {
  assert( solver != NULL);
  const double cp = solver->R_gas * solver->gamma / (solver->gamma - 1.0);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv = zone_ptr->cvobf[ibf];

    // Thermal resistance is denominated in [K/W]
    const double L_cv = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];  // This is an improved estimate of the length scale that balances robustness w accuracy
    const double k_cv = solver->loc_lam[icv]*cp; //gas thermal conductivity, [W/(m*K)]
    const double thermal_resistance_cv = L_cv/k_cv/zone_ptr->area_bf[ibf]; // [m]/[W/(m*K)]/[m^2] = [K/W]
    double thermal_resistance_overall = thermal_resistance_external_times_area/zone_ptr->area_bf[ibf] + thermal_resistance_cv; //[m^2*K/W]/[m^2] + [K/W] = [K/W]

    // heat transfer through layers is temperature difference between both ends of the thermal circuit divided by the overall thermal resitance
    const double q_wall = (solver->T[icv] - T_end_bc) / thermal_resistance_overall / zone_ptr->area_bf[ibf]; // [K]/[K/W]/[m^2] = [W/m^2]
    // heat transfer through each layer is constant
    //const double T_wall = solver->T[icv] - q_wall * thermal_resistance_cv*zone_ptr->area_bf[ibf]; // [K] - [W/m^2]*[K/W][m^2] = [K]

    const double visc_coeff = solver->mu_lam[icv]* zone_ptr->area_over_delta_bf[ibf];
    //Don't need to do this anymore b/c q has already been computed

    for (int i = 0; i < 3; ++i)
      rhs[icv].rhou[i] -= solver->p[icv] * zone_ptr->n_bf[ibf][i] + visc_coeff * (solver->u[icv][i] - u_bc[i]);

    // XXX if u_bc \neq 0, there should be a u \dot \tau term added here..
    //rhs[icv].rhoE -= k_coeff * (solver->T[icv] - T_wall); // [W/K]*[K] = [W]
    rhs[icv].rhoE -= q_wall*zone_ptr->area_bf[ibf]; // [W/m^2] * [m^2] = [W]
  }
}

void WallConductive::query(const string& param_str) const {

  // report an estimate of the local y^+ on the boundary and the integrated tau_wall
  double my_buf[5] = {0.0,0.0,0.0,0.0,0.0};

  const double cp = solver->R_gas * solver->gamma / (solver->gamma - 1.0);
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    if ( zone_ptr->area_over_delta_bf[ibf] > 0.0) {
      const int icv = zone_ptr->cvobf[ibf];

      // Shear stress stuff...
      const double y1         = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
      const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];

      double du_mag = 0.0;
      for (int i =0; i < 3; ++i)
        du_mag  += (solver->u[icv][i] - u_bc[i])*(solver->u[icv][i] - u_bc[i]);

      const double tau_wall   = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
      const double y_plus      = y1*sqrt(tau_wall*solver->rho[icv])/solver->mu_lam[icv];

      // Heat flux stuff...
      //const double L_cv                         = pow(solver->vol_cv[icv], 1.0 / 3.0);
      const double L_cv                         = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
      const double k_cv                         = solver->loc_lam[icv]* cp;
      const double thermal_resistance_cv        = L_cv/k_cv/zone_ptr->area_bf[ibf];
      const double thermal_resistance_overall   = thermal_resistance_external_times_area/zone_ptr->area_bf[ibf] + thermal_resistance_cv;
      const double q_wall                       = (solver->T[icv] - T_end_bc) / thermal_resistance_overall / zone_ptr->area_bf[ibf];
      const double T_wall                       = solver->T[icv] - q_wall * thermal_resistance_cv*zone_ptr->area_bf[ibf];
      //const double k_coeff                      = (solver->loc_lam[icv] + solver->loc_sgs[icv]) * zone_ptr->area_over_delta_bf[ibf] * cp; //[kg/m*s]*[m]*[W*s/kg*K] = [W/K]
      //const double q_wall                       = k_coeff*abs(solver->T[icv] - T_wall)/zone_ptr->area_bf[ibf];

      //NOTE: We are reporting the actual q_wall value (based on rhs.rhoE) rather than the laminar q wall value (which defines q)
      //There seems to be some inconsistency here...

      // Sum for reduction...
      my_buf[0] += zone_ptr->area_bf[ibf];
      my_buf[1] += zone_ptr->area_bf[ibf]*tau_wall;
      my_buf[2] += zone_ptr->area_bf[ibf]*y_plus;
      my_buf[3] += zone_ptr->area_bf[ibf]*q_wall;
      my_buf[4] += zone_ptr->area_bf[ibf]*T_wall;
    }
  }

  double buf[5];
  MPI_Reduce(my_buf,buf,5,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if ( mpi_rank == 0 ) {

    cout << zone_ptr->getName() << " time, int tau_wall dA, int q_w dA, avg y_plus, avg T_bc = "
         << solver->time  << "    "
         << buf[1]        << "    "
         << buf[3]        << "    "
         << buf[2]/buf[0] << "    "
         << buf[4]/buf[0] << endl;
  }
}

CtiRegister::CtiDataError WallConductive::funcEvalCtiData(CtiRegister::CtiData& v,
    const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  if (zone_ptr->isLeadingStrMatched(name)) {
    // we checked if the requested function is a potential match for the data
    // based on the zone name prefix.  we have to check for the actual data now..

    const string tau_wall_str = zone_ptr->getName() + ":" + "tau_wall";
    const string yp_str  = zone_ptr->getName() + ":" + "y_plus";
    const string qw_str  = zone_ptr->getName() + ":" + "q_wall";
    const string Tbc_str  = zone_ptr->getName() + ":" + "T_bc";

    // all of the following functions do not take any arguments -- check first
    if (name == tau_wall_str) {
      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {
        if (args.size() == 0) {
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

            const int icv           = zone_ptr->cvobf[ibf];
            const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];
            double du_mag = 0.0;
            for (int i =0; i < 3; ++i)
              du_mag  += (solver->u[icv][i] - u_bc[i])*(solver->u[icv][i] - u_bc[i]);

            tmp[ibf]                = visc_coeff*sqrt(du_mag)/zone_ptr->area_bf[ibf];
          }
        }
        else if (args.size() == 1) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          double dir[3] = {0.0,0.0,0.0};
          if (arg->getType() == D3_DATA) {
            FOR_I3 dir[i] = arg->d3(i);
            NORMALIZE(dir);
          }
          else if (arg->getType() == D_DATA) {
            const int index = int(arg->d());
            if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;
            dir[index] = 1.0;
          }
          else return CTI_DATA_NOT_VALID;

          // loop faces and compute
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
            const int icv           = zone_ptr->cvobf[ibf];
            const double visc_coeff = solver->mu_lam[icv]*zone_ptr->area_over_delta_bf[ibf];

            double local_u[3] = DIFF(solver->u[icv],u_bc);
            tmp[ibf]          = visc_coeff*DOT_PRODUCT(dir,local_u)/zone_ptr->area_bf[ibf];
          }
        }
        else return CTI_DATA_ARG_COUNT;
      }
      return CTI_DATA_OK;
    }

    else if (name == yp_str) {
      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          const int icv           = zone_ptr->cvobf[ibf];
          const double y1         = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
          const double visc_coeff = solver->mu_lam[icv]/y1;

          double du_mag = 0.0;
          for (int i =0; i < 3; ++i)
            du_mag  += (solver->u[icv][i] - u_bc[i])*(solver->u[icv][i] - u_bc[i]);

          const double tau_wall      = visc_coeff*sqrt(du_mag);
          tmp[ibf]                = y1* sqrt(tau_wall*solver->rho[icv])/solver->mu_lam[icv];
        }
      }
      return CTI_DATA_OK;
    }

    else if (name == qw_str) {

      if (args.size() != 0) return CTI_DATA_ARG_COUNT;
      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {
        const double cp = solver->R_gas*solver->gamma/(solver->gamma-1.0);
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          const int icv        = zone_ptr->cvobf[ibf];
          //const double L_cv = pow(solver->vol_cv[icv], 1.0 / 3.0);  // cv length scale
          const double L_cv = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
          const double k_cv = solver->loc_lam[icv]*cp;  // gas thermal conductivity, [W/(m*K)]
          const double thermal_resistance_cv = L_cv/k_cv/zone_ptr->area_bf[ibf]; // [m]/[W/(m*K)]/[m^2] = [K/W]
          double thermal_resistance_overall = thermal_resistance_external_times_area/zone_ptr->area_bf[ibf] + thermal_resistance_cv; //[m^2*K/W]/[m^2] + [K/W] = [K/W]

          tmp[ibf] = (solver->T[icv] - T_end_bc) / thermal_resistance_overall / zone_ptr->area_bf[ibf]; // [K]/[K/W]/[m^2] = [W/m^2]
          //const double T_wall = solver->T[icv] - q * thermal_resistance_cv*zone_ptr->area_bf[ibf]; // [K] - [W/m^2]*[K/W][m^2] = [K]

          //const double k_coeff = (solver->loc_lam[icv] + solver->loc_sgs[icv]) * zone_ptr->area_over_delta_bf[ibf] * cp; //[kg/m*s]*[m]*[W*s/kg*K] = [W/K]
          //tmp[ibf]             = k_coeff*abs(solver->T[icv] - T_wall)/zone_ptr->area_bf[ibf];
        }
      }
      return CTI_DATA_OK;
    }

    else if (name == Tbc_str) {
      if (args.size() != 0) return CTI_DATA_ARG_COUNT;
      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {
        const double cp = solver->R_gas * solver->gamma / (solver->gamma - 1.0);
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ ibf) {
          const int icv = zone_ptr->cvobf[ibf];
          //const double L_cv = pow(solver->vol_cv[icv], 1.0/3.0);
          const double L_cv = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
          const double k_cv = solver->loc_lam[icv]*cp;
          const double thermal_resistance_cv = L_cv/k_cv/zone_ptr->area_bf[ibf];
          const double thermal_resistance_overall = thermal_resistance_external_times_area/zone_ptr->area_bf[ibf] + thermal_resistance_cv;
          const double q_wall = (solver->T[icv] - T_end_bc)/thermal_resistance_overall / zone_ptr->area_bf[ibf];
          tmp[ibf] = solver->T[icv] = q_wall * thermal_resistance_cv*zone_ptr->area_bf[ibf];
        }
      }
      return CTI_DATA_OK;
    }
  }

  return IdealGasBc::funcEvalCtiData(v,name,args,b_eval_func);
}


void Cbc::query(const string& param_str) const {

  // report integrated values for the p_bc, rho_bc, u_bc
  double my_buf[4] = {0.0,0.0,0.0,0.0};
  double buf[4];

  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv         = zone_ptr->cvobf[ibf];
    const double area     = MAG(zone_ptr->n_bf[ibf]);
    const double undA_bc  = DOT_PRODUCT(zone_ptr->n_bf[ibf],solver->u[icv]);
    my_buf[0]            += area;
    my_buf[1]            += area*solver->rho[icv];
    my_buf[2]            += undA_bc*solver->rho[icv];
    my_buf[3]            += area*solver->p[icv];
  }

  MPI_Reduce(my_buf,buf,4,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  if ( mpi_rank == 0 ) {
    cout << " > " << zone_ptr->getName() << " :  time, rho_bc, mdot, p_bc = "
         << solver->time << "   " << buf[1]/buf[0] << "    "
         << buf[2] << "    " << buf[3]/buf[0] << "    " << endl;
  }
}

void Cbc::addBoundaryFlux(IdealGasRhs* rhs) const {
  assert( solver != NULL);
  const double gamma_ = solver->gamma;
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv = zone_ptr->cvobf[ibf];
    IdealGasRhs flux;
    calcRiemannFlux(flux,zone_ptr->n_bf[ibf],solver->cv_light[icv],
                    bf_state[ibf],gamma_, gamma_);

    mf[ibf]       = flux.rho;

    rhs[icv].rho -= flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv].rhou[i] -= flux.rhou[i];
    rhs[icv].rhoE -= flux.rhoE;
  }
}

void Cbc::addScalarBoundaryFlux(double * rhs_sc) const {
  BcTemplate<IdealGasSolver,IdealGasRhs>::_addScalarBoundaryFlux(rhs_sc);
}

void computeBfWallDistance(double * wall_dist, const BfZone * zone_ptr, const IdealGasSolver* solver) {

  const int ncv = solver->ncv;
  double * cv_flag = new double[ncv];
  FOR_ICV cv_flag[icv] = -2.0;

  const int nbf = zone_ptr->nbf;
  // flag cells adjacent to faces of interest
  FOR_IBF {
    cv_flag[zone_ptr->cvobf[ibf]] = -1.0;
  }

  // We have flagged the cvs adjacent to this zone and we can loop over other boundary conditions to find cells at the corner of this zone and wall BCs. In only these cells we will compute the wall distance
  for(vector<IdealGasBc*>::const_iterator it = solver->bcs.begin(); it != solver->bcs.end(); ++it) {
    string s=(*it)->getName();
    Param* p=getParam(s);
    if ((p->getString()=="WALL_ADIABATIC") || (p->getString()=="WALL_ISOTHERMAL") || (p->getString()=="WM_ALG_ADIABATIC") || (p->getString()=="WM_ALG_ISOTHERMAL")) {
      // TODO This must recognize all wall bcs. How?
      BfZone* zone_ptr_it=(*it)->zone_ptr;
      for (int ibf=0; ibf < zone_ptr_it->nbf; ++ibf) {
        const int icv = zone_ptr_it->cvobf[ibf];

        if (cv_flag[icv] == -1.0) {
          // first time hit, identify as a corner cell
          cv_flag[icv] = HUGE_VAL;
        }

        if (cv_flag[icv] > 0.0) {
          // we are at a zone corner cv, so compute wall-distance and store minimum in cv_flag
          const double dist = DIST(zone_ptr_it->x_bf[ibf],solver->x_cv[icv]);
          if (dist < cv_flag[icv]) {
            cv_flag[icv] = dist;
          }
        }
      }
    }
  }  // bcs iterator

  solver->computeBfDistanceFromFlaggedCvs(wall_dist,zone_ptr,cv_flag);
  DELETE(cv_flag);
}

void CbcProfile::initData() {
  Param * param  = getParam(getName());
  assert( param != NULL);

  int type = -1;
  double type_double[2] = {0.0,0.0};

  FluentBPReader profile;
  int iarg = 1;
  while (iarg < param->size()) {
    string token = param->getString(iarg++);
    if (token == "FILE") {
      profile.init(param->getString(iarg++));
    }
    else if (token == "TYPE") {
      const string _type = param->getString(iarg++);
      if (_type == "UPT") {
        type = 0;
      }
      else if (_type == "RUP") {
        type = 1;
      }
      else if (_type == "MPT") {
        type = 2;
        type_double[0] = param->getDouble(iarg++);
      }
      else if (_type == "U_CONSTANT_TP") {
        type = 3;
        type_double[0] = param->getDouble(iarg++);
        type_double[1] = param->getDouble(iarg++);
      }
      else {
        CERR("unsupported CBC_PROFILE TYPE " << _type << "; please choose from UPT, RUP, MPT, U_CONSTANT_TP");
      }
    }
    else if (token == "FROM_MEAN") {
      profile.setUseMean(true);
    }
  }

  if (!profile.isInitialized()) {
    CERR("boundary profile was not properly initialized");
  }
  // profile.printVarsInFile();

  if (profile.getType() == LINE_WD) {
    // need to pass wall-distance for each ibf
    double * wall_dist = new double[zone_ptr->nbf];

    // compute wall distance
    computeBfWallDistance(wall_dist,zone_ptr,solver);

    profile.setPoints(wall_dist,zone_ptr->nbf);
    DELETE(wall_dist);
  }
  else {
    // pass in actual ibf locations
    double (*x_fa_zone)[3] = new double[zone_ptr->nbf][3];

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      FOR_I3 x_fa_zone[ibf][i] = zone_ptr->x_bf[ibf][i];
    }
    profile.setPoints(x_fa_zone,zone_ptr->x_global,zone_ptr->nbf);
    DELETE(x_fa_zone);
  }

  if (type == -1) {
    CERR("CBC_PROFILE \"TYPE\" must be specified");
  }
  else if (type == 0) {
    profile.ensureVar("x-velocity");
    profile.ensureVar("y-velocity");
    profile.ensureVar("z-velocity");
    profile.ensureVar("absolute-pressure");
    profile.ensureVar("temperature");

    assert( bf_state == NULL);
    bf_state = new IdealGasState[zone_ptr->nbf];

    for (int ibf = 0, nbf=zone_ptr->nbf; ibf < nbf; ibf++) {
      bf_state[ibf].u[0] = profile.getData(ibf,"x-velocity");
      bf_state[ibf].u[1] = profile.getData(ibf,"y-velocity");
      bf_state[ibf].u[2] = profile.getData(ibf,"z-velocity");
      bf_state[ibf].p    = profile.getData(ibf,"absolute-pressure");
      const double T_ibf = profile.getData(ibf,"temperature");
      bf_state[ibf].sp_vol = solver->R_gas*T_ibf/bf_state[ibf].p ;
      bf_state[ibf].h      = solver->R_gas*T_ibf*solver->gamma/(solver->gamma - 1.0);
    }
  }
  else if (type == 1) {
    profile.ensureVar("x-velocity");
    profile.ensureVar("y-velocity");
    profile.ensureVar("z-velocity");
    profile.ensureVar("absolute-pressure");
    profile.ensureVar("density");

    assert( bf_state == NULL);
    bf_state = new IdealGasState[zone_ptr->nbf];

    for (int ibf = 0, nbf=zone_ptr->nbf; ibf < nbf; ibf++) {
      bf_state[ibf].u[0] = profile.getData(ibf,"x-velocity");
      bf_state[ibf].u[1] = profile.getData(ibf,"y-velocity");
      bf_state[ibf].u[2] = profile.getData(ibf,"z-velocity");
      bf_state[ibf].p    = profile.getData(ibf,"absolute-pressure");
      bf_state[ibf].sp_vol = 1.0/profile.getData(ibf,"density");
      const double T_ibf = bf_state[ibf].p*bf_state[ibf].sp_vol/solver->R_gas;
      bf_state[ibf].h      = solver->R_gas*T_ibf*solver->gamma/(solver->gamma - 1.0);
    }
  }
  else if (type == 2) {
    profile.ensureVar("absolute-pressure");
    profile.ensureVar("temperature");
    if (type_double[0] == 0.0) CERR("Mass flux cannot be zero");

    assert( bf_state == NULL);
    bf_state = new IdealGasState[zone_ptr->nbf];

    for (int ibf = 0, nbf=zone_ptr->nbf; ibf < nbf; ibf++) {
      double T_ibf;
      bf_state[ibf].p    = profile.getData(ibf,"absolute-pressure");
      T_ibf = profile.getData(ibf,"temperature");

      const double rho_ibf = bf_state[ibf].p/T_ibf/solver->R_gas;
      const double un_ibf  = -1.0*type_double[0]/(rho_ibf*zone_ptr->area_global);
      const double bf_area = MAG(zone_ptr->n_bf[ibf]);
      FOR_I3 bf_state[ibf].u[i] = un_ibf*zone_ptr->n_bf[ibf][i]/bf_area;
      bf_state[ibf].sp_vol = 1.0/rho_ibf;
      bf_state[ibf].h      = solver->R_gas*T_ibf*solver->gamma/(solver->gamma - 1.0);
    }
  }
  else if (type == 3) {
    profile.ensureVar("x-velocity");
    profile.ensureVar("y-velocity");
    profile.ensureVar("z-velocity");
    if (type_double[0] == 0.0) CERR("temperature cannot be zero");
    if (type_double[1] == 0.0) CERR("pressure cannot be zero");

    assert( bf_state == NULL);
    bf_state = new IdealGasState[zone_ptr->nbf];

    for (int ibf = 0, nbf=zone_ptr->nbf; ibf < nbf; ibf++) {
      bf_state[ibf].p        = type_double[1];
      const double T_ibf     = type_double[0];
      bf_state[ibf].sp_vol   = (solver->R_gas*T_ibf)/bf_state[ibf].p;  // 1.0/rho
      bf_state[ibf].h        = solver->R_gas*T_ibf*solver->gamma/(solver->gamma - 1.0);
      bf_state[ibf].u[0]     = profile.getData(ibf, "x-velocity");
      bf_state[ibf].u[1]     = profile.getData(ibf, "y-velocity");
      bf_state[ibf].u[2]     = profile.getData(ibf, "z-velocity");
    }
  }
  else CERR("unrecognized type somehow set...");

  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<IdealGasSolver,IdealGasRhs>::parseScalarBc(param);

  // profile.dumpTecplot("profile_test_dump.dat");
}

void CbcUpt::initData() {
  Param * param  = getParam(getName());
  assert( param != NULL);

  double u_bc[3], p_bc, T_bc;

  u_bc[0]        = param->getDouble(1);
  u_bc[1]        = param->getDouble(2);
  u_bc[2]        = param->getDouble(3);
  p_bc           = param->getDouble(4);
  T_bc           = param->getDouble(5);

  assert( bf_state == NULL);
  bf_state = new IdealGasState[zone_ptr->nbf];

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    FOR_I3 bf_state[ibf].u[i] = u_bc[i];
    bf_state[ibf].p      = p_bc;
    bf_state[ibf].sp_vol = solver->R_gas*T_bc/p_bc;
    bf_state[ibf].h      = solver->R_gas*T_bc*solver->gamma/(solver->gamma - 1.0);
  }

  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<IdealGasSolver,IdealGasRhs>::parseScalarBc(param);

}

void CbcMpt::initData() {
  Param* param   = getParam(getName());
  assert( param != NULL);
  double mdot_bc = param->getDouble(1);
  double p_bc    = param->getDouble(2);
  double T_bc    = param->getDouble(3);

  //double my_area = 0.0, area = 0.0;
  //for (int ibf =0; ibf < zone_ptr->nbf; ++ibf)
  //  my_area += MAG(zone_ptr->n_bf[ibf]);

  // if the assert passes, then we can avoid the reduction here..
  //MPI_Allreduce(&my_area,&area,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
  //assert( abs(area-zone_ptr->area_global) < 1.0e-15);

  double rho_bc = p_bc/T_bc/solver->R_gas;
  double un_bc  = -mdot_bc/(rho_bc*zone_ptr->area_global);

  assert( bf_state == NULL);
  bf_state = new IdealGasState[zone_ptr->nbf];
  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    const double bf_area = MAG(zone_ptr->n_bf[ibf]);
    FOR_I3 bf_state[ibf].u[i] = un_bc*zone_ptr->n_bf[ibf][i]/bf_area;
    bf_state[ibf].sp_vol = 1.0/rho_bc;
    bf_state[ibf].p      = p_bc;
    bf_state[ibf].h      = solver->R_gas*T_bc*solver->gamma/(solver->gamma - 1.0);
  }

  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<IdealGasSolver,IdealGasRhs>::parseScalarBc(param);

}

void CbcRunp::initData() {
  Param * param = getParam(getName());
  assert( param != NULL);

  double rho_bc  = param->getDouble(1);
  double un_bc   = param->getDouble(2);
  double p_bc    = param->getDouble(3);

  double T_bc    = p_bc/rho_bc/solver->R_gas;

  assert( bf_state == NULL);
  bf_state = new IdealGasState[zone_ptr->nbf];
  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    const double bf_area = MAG(zone_ptr->n_bf[ibf]);
    for (int i =0; i < 3; ++i)
      bf_state[ibf].u[i] = -un_bc*zone_ptr->n_bf[ibf][i]/bf_area;
    bf_state[ibf].sp_vol = 1.0/rho_bc;
    bf_state[ibf].p      = p_bc;
    bf_state[ibf].h      = solver->R_gas*T_bc*solver->gamma/(solver->gamma - 1.0);
  }

  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<IdealGasSolver,IdealGasRhs>::parseScalarBc(param);

}

void CbcRup::initData() {

  Param * param = getParam(getName());
  assert( param != NULL);

  double rho_bc, p_bc, T_bc;
  rho_bc  = param->getDouble(1);
  u_bc[0] = param->getDouble(2);
  u_bc[1] = param->getDouble(3);
  u_bc[2] = param->getDouble(4);
  p_bc    = param->getDouble(5);

  T_bc    = p_bc/rho_bc/solver->R_gas;

  assert( bf_state == NULL);
  bf_state = new IdealGasState[zone_ptr->nbf];

  u_fr          = new double[zone_ptr->nbf][3];
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
    for (int i = 0; i < 3; ++i)
      u_fr[ibf][i] = 0.0;

  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    FOR_I3 bf_state[ibf].u[i] = u_bc[i];
    bf_state[ibf].sp_vol = 1.0/rho_bc;
    bf_state[ibf].p      = p_bc;
    bf_state[ibf].h      = solver->R_gas*T_bc*solver->gamma/(solver->gamma - 1.0);
  }

  if ( (param->size() == 7) && (param->getString(6) == "STATIONARY_FRAME")) {

    if ( solver->frame_rotation != NULL) {

      omega = MAG(solver->frame_rotation);
      FOR_I3 rot_axis[i] = solver->frame_rotation[i]/omega;

      double ubc_x_axis[3] = CROSS_PRODUCT(u_bc, rot_axis);
      FOR_I3 ubc_x_axis[i] /= MAG(u_bc);
      if ( MAG(ubc_x_axis)>1.0e-12 ) {
	b_aligned = false;
	if (mpi_rank == 0) cout << getName() << " - CBC_RUP: u_bc and the rotation axis is not aligned." << endl;
      }

      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

        u_fr[ibf][0] = -solver->frame_rotation[1]*zone_ptr->x_bf[ibf][2] +
                        solver->frame_rotation[2]*zone_ptr->x_bf[ibf][1];

        u_fr[ibf][1] = -solver->frame_rotation[2]*zone_ptr->x_bf[ibf][0] +
                        solver->frame_rotation[0]*zone_ptr->x_bf[ibf][2];

        u_fr[ibf][2] = -solver->frame_rotation[0]*zone_ptr->x_bf[ibf][1] +
                        solver->frame_rotation[1]*zone_ptr->x_bf[ibf][0];

	if (!b_aligned)
	  calcBfState(ibf, solver->time);
	else
	  FOR_I3 bf_state[ibf].u[i] += u_fr[ibf][i];
      }
    }
  }

  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<IdealGasSolver,IdealGasRhs>::parseScalarBc(param);

}

inline void CbcRup::calcBfState(const int ibf, const double time) {

  assert(!b_aligned);
  // note: the bf states sp_vol, p and h are initialized in initData and would not change
  // the only thing that needs an update is bf_state.u
  // rotate the u_bc(in stationary frame) in the opposite direction
  const double theta = -omega*time;
  // use the Rodrigues formula
  // U_rot = U cos(theta) + (K x U) sin(theta) + K (K dot U) (1-cos(theta))
  const double axis_x_ubc[3] = CROSS_PRODUCT(rot_axis,u_bc);
  const double axis_dot_ubc  = DOT_PRODUCT(rot_axis,u_bc);
  const double cos_theta = cos(theta);
  const double sin_theta = sin(theta);
  for (int i = 0; i < 3; ++i)
    bf_state[ibf].u[i] = u_bc[i]*cos_theta + axis_x_ubc[i]*sin_theta +
                         rot_axis[i]*axis_dot_ubc*(1.0-cos_theta) + u_fr[ibf][i];
}

void CbcRup::calcRhs(const double time, const int rk_stage) {
  if (!b_aligned) {
    // b_aligned can be false iff STATIONARY_FRAME is specified
    // and solver->frame_rotation is not NULL
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
	calcBfState(ibf, time);
    }
  }
}

void CbcExtrapolate::initData() {

  assert( bf_state == NULL);
  bf_state      = new IdealGasState[zone_ptr->nbf];
  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  //BcTemplate<IdealGasSolver,IdealGasRhs>::parseScalarBc(param);

}

void CbcTotalPt::initData() {

  Param * param = getParam(getName());
  total_p       = param->getDouble(1);
  total_t       = param->getDouble(2);
  t_relax       = param->getDouble(3);

  assert( bf_state == NULL);
  bf_state      = new IdealGasState[zone_ptr->nbf];
  un_bc         = new double[zone_ptr->nbf];

  u_fr          = new double[zone_ptr->nbf][3];
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
    for (int i = 0; i < 3; ++i)
      u_fr[ibf][i] = 0.0;

  if ( (param->size() >= 5) && (param->getString(4) == "USE_AVG"))
    b_use_avg = true;

  if ( ((param->size() == 5) && (param->getString(4) == "STATIONARY_FRAME")) ||
       ((param->size() == 6) && (param->getString(5) == "STATIONARY_FRAME")) ) {

    if ( solver->frame_rotation != NULL) {

      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

        u_fr[ibf][0] = -solver->frame_rotation[1]*zone_ptr->x_bf[ibf][2] +
                        solver->frame_rotation[2]*zone_ptr->x_bf[ibf][1];

        u_fr[ibf][1] = -solver->frame_rotation[2]*zone_ptr->x_bf[ibf][0] +
                        solver->frame_rotation[0]*zone_ptr->x_bf[ibf][2];

        u_fr[ibf][2] = -solver->frame_rotation[0]*zone_ptr->x_bf[ibf][1] +
                        solver->frame_rotation[1]*zone_ptr->x_bf[ibf][0];
      }
    }
  }

  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<IdealGasSolver,IdealGasRhs>::parseScalarBc(param);
}

inline void CbcTotalPt::calcBfState(const int ibf, const double unit_n[3]) {

  double ut[3];
  double ufr_n = DOT_PRODUCT(u_fr[ibf],unit_n);
  for (int i = 0; i < 3; ++i)
    ut[i] = u_fr[ibf][i] - ufr_n*unit_n[i];

  for (int i = 0; i < 3; ++i)
    bf_state[ibf].u[i] = un_bc[ibf]*unit_n[i] + ut[i];

  const double total_rh = total_p/solver->R_gas/total_t;
  const double invgm1   = 1.0/( solver->gamma - 1.0);
  const double cp       = solver->R_gas * solver->gamma*invgm1;
  const double T_bc     = total_t - 0.5*un_bc[ibf]*un_bc[ibf]/cp;
  const double rho_bc   = total_rh* pow(T_bc/total_t, invgm1);
  const double p_bc     = total_p * pow(T_bc/total_t, solver->gamma*invgm1);

  bf_state[ibf].sp_vol  = 1.0/rho_bc;
  bf_state[ibf].p       = p_bc;
  bf_state[ibf].h       = cp*T_bc;
}

inline void CbcExtrapolate::calcBfState(const int ibf, const double unit_n[3]) {
  const double cp = solver->R_gas*solver->gamma/(solver->gamma-1.0);
  const int icv = solver->cvobf[ibf];
  bf_state[ibf].sp_vol = 1.0/solver->rho[icv];
  for (int i = 0; i < 3; ++i )
    bf_state[ibf].u[i] = solver->u[icv][i];
  bf_state[ibf].p = solver->p[icv];
  bf_state[ibf].h = cp*solver->T[icv];

}

void CbcTotalPt::initialHook() {
  if ( !zone_ptr->checkDataFlag("un_bc")) {
    if ( mpi_rank == 0 )
      cout << " > zone : " << getName() << " : reinitializing un_bc for cbc total pt .. " << endl;

    if (b_use_avg) {
      double my_buf[2] = {0.0, 0.0};    // rhoA, massflux
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
	const int icv        = zone_ptr->cvobf[ibf];
	my_buf[0] += solver->rho[icv] * MAG(zone_ptr->n_bf[ibf]);
	my_buf[1] += solver->rho[icv] * DOT_PRODUCT(solver->u[icv],zone_ptr->n_bf[ibf]);
      }
      double buf[2];
      MPI_Allreduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,mpi_comm);
      un_bc_avg = buf[1]/buf[0];
    }

    for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv        = zone_ptr->cvobf[ibf];
      const double bf_area = MAG(zone_ptr->n_bf[ibf]);
      double unit_n[3]     = {zone_ptr->n_bf[ibf][0]/bf_area,
                              zone_ptr->n_bf[ibf][1]/bf_area,
                              zone_ptr->n_bf[ibf][2]/bf_area};
      if (b_use_avg) {
	un_bc[ibf]         = un_bc_avg;
      }
      else {
	un_bc[ibf]         = DOT_PRODUCT(unit_n,solver->u[icv]);
      }
      calcBfState(ibf,unit_n);
    }
  }
  else {
    if (b_use_avg) {
      double my_buf[2] = {0.0, 0.0};    // rhoA, massflux
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
	const int icv        = zone_ptr->cvobf[ibf];
	const double bf_area = MAG(zone_ptr->n_bf[ibf]);
	my_buf[0] += solver->rho[icv] * bf_area;
	my_buf[1] += solver->rho[icv] * un_bc[ibf] * bf_area;
      }
      double buf[2];
      MPI_Allreduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,mpi_comm);
      un_bc_avg = buf[1]/buf[0];
    }

    for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
      const double bf_area = MAG(zone_ptr->n_bf[ibf]);
      double unit_n[3]     = {zone_ptr->n_bf[ibf][0]/bf_area,
                              zone_ptr->n_bf[ibf][1]/bf_area,
                              zone_ptr->n_bf[ibf][2]/bf_area};
      if (b_use_avg) un_bc[ibf] = un_bc_avg;
      calcBfState(ibf,unit_n);
    }
  }
}

void CbcExtrapolate::initialHook() {

  if ( mpi_rank == 0) {
    cout << " im here.. " << endl;
    cout.flush();
  }

  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    const double bf_area = MAG(zone_ptr->n_bf[ibf]);
    double unit_n[3]     = {zone_ptr->n_bf[ibf][0]/bf_area,
      zone_ptr->n_bf[ibf][1]/bf_area,
      zone_ptr->n_bf[ibf][2]/bf_area};
    calcBfState(ibf,unit_n);
  }

  if ( mpi_rank == 0) {
    cout << " im out. " << endl;
    cout.flush();
  }

}

void CbcTotalPt::rk3Step(const double dt, const int rk_stage) {
  if ( rk_stage == 1) { // updates the value at the beginning of the step...
    
    const double eps   = 1.0 / (t_relax/dt + 1.0);
    
    if (b_use_avg) {
      double my_buf[2] = {0.0, 0.0};    // rhoA, massflux
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
	const int icv        = zone_ptr->cvobf[ibf];
	my_buf[0] += solver->rho[icv] * MAG(zone_ptr->n_bf[ibf]);
	my_buf[1] += solver->rho[icv] * DOT_PRODUCT(solver->u[icv],zone_ptr->n_bf[ibf]);
      }
      double buf[2];
      MPI_Allreduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,mpi_comm);
      const double un_bc_avg0 = buf[1]/buf[0];
      un_bc_avg  = (1.0-eps)*un_bc_avg + eps*un_bc_avg0;
    }

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv        = zone_ptr->cvobf[ibf];
      const double bf_area = MAG(zone_ptr->n_bf[ibf]);
      double unit_n[3]     = {zone_ptr->n_bf[ibf][0]/bf_area,
                              zone_ptr->n_bf[ibf][1]/bf_area,
                              zone_ptr->n_bf[ibf][2]/bf_area};

      if (b_use_avg) {
	un_bc[ibf] = un_bc_avg;
      }
      else {
	const double un0 = DOT_PRODUCT(solver->u[icv],unit_n);
	un_bc[ibf]         = (1.0-eps)*un_bc[ibf] + eps*un0;
      }

      // with the new un_bc, we are re-compute the bf_state.  this
      // is held frozen for the rest of the step...
      calcBfState(ibf,unit_n);
    }
  }
}

void CbcExtrapolate::rk3Step(const double dt, const int rk_stage) {

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const double bf_area = MAG(zone_ptr->n_bf[ibf]);
    double unit_n[3]     = {zone_ptr->n_bf[ibf][0]/bf_area,
      zone_ptr->n_bf[ibf][1]/bf_area,
      zone_ptr->n_bf[ibf][2]/bf_area};
    calcBfState(ibf,unit_n);
  }
}

int Nscbc::storeData(double * buf) const {

  if (buf) {
    int count = 0;
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      buf[count++] = rho_bc[ibf];
      buf[count++] = p_bc[ibf];
      FOR_I3 buf[count++] = u_bc[ibf][i];
    }
    assert(count == zone_ptr->nbf*5);
  }

  return zone_ptr->nbf*5; // rho,p,u(3)

}

int Nscbc::restoreData(double * buf) {

  assert(buf);
  int count = 0;
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    rho_bc[ibf] = buf[count++];
    p_bc[ibf] = buf[count++];
    FOR_I3 u_bc[ibf][i] = buf[count++];
  }
  assert(count == zone_ptr->nbf*5);

  return zone_ptr->nbf*5;
}

Nscbc::~Nscbc() {
  DELETE(p_bc);
  DELETE(rho_bc);
  DELETE(u_bc);
  DELETE(rhs_bc);
}

void Nscbc::initData() {
  p_bc   = new double[zone_ptr->nbf];
  rho_bc = new double[zone_ptr->nbf];
  u_bc   = new double[zone_ptr->nbf][3];
  rhs_bc = new double[zone_ptr->nbf][3][5];

  // TODO: do we need to do this?
  // yes we do. STB,CBI,JOB
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    for (int i =0; i < 3; ++i)
      for (int j =0; j < 5; ++j)
        rhs_bc[ibf][i][j] = 0.0;
  }
}

void Nscbc::initialHook() {

  if ( (!zone_ptr->checkDataFlag("p_bc")) ||
       (!zone_ptr->checkDataFlag("rho_bc")) ||
       (!zone_ptr->checkDataFlag("u_bc"))     ) {

    // if any of the variables are not instantiated, then
    // it is considered that all of the variables are not instantiated..
    if ( mpi_rank == 0 )
      cout << " > zone : " << getName() << " : reinitializing nscbc vars .. " << endl;

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv    = zone_ptr->cvobf[ibf];
      rho_bc[ibf]      = solver->rho[icv];
      p_bc[ibf]        = solver->p[icv];
      for (int i =0; i < 3; ++i)
        u_bc[ibf][i]   = solver->u[icv][i];
    }
  }
}

void Nscbc::rk3Step(const double dt, const int rk_stage) {
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    for (int i =0; i < rk_stage; ++i) {
      rho_bc[ibf] += erk_wgt3[rk_stage-1][i]*dt*rhs_bc[ibf][i][0];
      p_bc[ibf]   += erk_wgt3[rk_stage-1][i]*dt*rhs_bc[ibf][i][1];
      for (int j =0; j < 3; ++j)
        u_bc[ibf][j] += erk_wgt3[rk_stage-1][i]*dt*rhs_bc[ibf][i][2+j];
    }
  }
}

void Nscbc::addBoundaryFlux(IdealGasRhs* rhs) const {

  assert( solver != NULL);
  assert( mf != NULL);

  const double gamma_ = solver->gamma;
  const double gogm1  = solver->gamma/(solver->gamma-1.0);

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    const int icv = zone_ptr->cvobf[ibf];
    IdealGasRhs flux;

    // reconstruct the packed state for the riemann flux..
    // if this becomes performance critical, we can define
    // an alternative riemann flux function signature ..
    IdealGasState bf_state;
    for (int i =0; i < 3; ++i)
      bf_state.u[i] = u_bc[ibf][i];
    bf_state.sp_vol = 1.0/rho_bc[ibf];
    bf_state.p      = p_bc[ibf];
    bf_state.h      = gogm1*p_bc[ibf]*(1.0/rho_bc[ibf]);

    calcRiemannFlux(flux,zone_ptr->n_bf[ibf],solver->cv_light[icv],bf_state,gamma_, gamma_);

    mf[ibf]       = flux.rho;

    rhs[icv].rho -= flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv].rhou[i] -= flux.rhou[i];
    rhs[icv].rhoE -= flux.rhoE;
  }
}

void Nscbc::addScalarBoundaryFlux(double* rhs_sc) const {
  BcTemplate<IdealGasSolver,IdealGasRhs>::_addScalarBoundaryFlux(rhs_sc);
}

void Nscbc::query(const string& param_str) const {

  // report integrated values for the p_bc, rho_bc, u_bc
  double my_buf[4] = {0.0,0.0,0.0,0.0};
  double buf[4];

  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    const double area    = MAG(zone_ptr->n_bf[ibf]);
    const double undA_bc = DOT_PRODUCT(zone_ptr->n_bf[ibf],u_bc[ibf]);
    my_buf[0]        += area;
    my_buf[1]        += area*rho_bc[ibf];
    my_buf[2]        += undA_bc*rho_bc[ibf];
    my_buf[3]        += area*p_bc[ibf];
  }

  MPI_Reduce(my_buf,buf,4,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  if ( mpi_rank == 0 ) {
    cout << " > " << zone_ptr->getName() << " :  time, rho_bc, mdot, p_bc = "
         << solver->time << "   " << buf[1]/buf[0] << "    "
         << buf[2] << "    " << buf[3]/buf[0] << "    " << endl;
  }
}

void NscbcMt::initData() {

  Nscbc::initData(); // takes care of u_bc,rho_bc,p_bc registration..

  Param * param        = getParam(getName());
  const double mdot_bc = param->getDouble(1);

  if (param->getString(2) == "T_TOTAL") {

    T_bc = -param->getDouble(3);
    if ( mpi_rank == 0)
      cout << " > zone : " << zone_ptr->getName() << " using T_total: " << T_bc << endl;

  }
  else {

    T_bc                 = param->getDouble(2);
    if ( mpi_rank == 0 )
      cout << " > zone : " << zone_ptr->getName() << " using T: " << T_bc << endl;

  }

  // the following can be removed shortly... use the area_global instead..
  double my_area = 0.0, area = 0.0;
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
    my_area += MAG(zone_ptr->n_bf[ibf]);
  MPI_Allreduce(&my_area,&area,1,MPI_DOUBLE,MPI_SUM,mpi_comm);

  rhoun_bc = -mdot_bc/area;

  u_fr          = new double[zone_ptr->nbf][3];
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
    for (int i = 0; i < 3; ++i)
      u_fr[ibf][i] = 0.0;

  if ( solver->frame_rotation != NULL) {

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      u_fr[ibf][0] = -solver->frame_rotation[1]*zone_ptr->x_bf[ibf][2] +
                      solver->frame_rotation[2]*zone_ptr->x_bf[ibf][1];

      u_fr[ibf][1] = -solver->frame_rotation[2]*zone_ptr->x_bf[ibf][0] +
                      solver->frame_rotation[0]*zone_ptr->x_bf[ibf][2];

      u_fr[ibf][2] = -solver->frame_rotation[0]*zone_ptr->x_bf[ibf][1] +
                      solver->frame_rotation[1]*zone_ptr->x_bf[ibf][0];
    }
  }

  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<IdealGasSolver,IdealGasRhs>::parseScalarBc(param);

}

void NscbcMt::calcRhs(const double time, const int rk_stage) {

  double static_T_bc;

  if ( T_bc > 0.0) {

    static_T_bc = T_bc;

  } else {

    double my_buf[2] = {0.0,0.0};
    double buf[2];

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      const int icv      = zone_ptr->cvobf[ibf];
      const double un_bc = rhoun_bc/solver->rho[icv];
      const double Ma    = un_bc/solver->sos[icv];
      my_buf[0]         += zone_ptr->area_bf[ibf] * Ma*Ma;
      my_buf[1]         += zone_ptr->area_bf[ibf];

    }

    MPI_Allreduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,mpi_comm);

    const double Ma2_avg = buf[0]/buf[1];
    assert( Ma2_avg > 0.0);

    static_T_bc = -T_bc / ( 1.0 + 0.5*(solver->gamma - 1.0)*Ma2_avg);

  }

  const double RT_bc = static_T_bc*solver->R_gas;

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const double bf_area = MAG(zone_ptr->n_bf[ibf]);
    double unit_n[3]     = {zone_ptr->n_bf[ibf][0]/bf_area,
                            zone_ptr->n_bf[ibf][1]/bf_area,
                            zone_ptr->n_bf[ibf][2]/bf_area};

    // ensure consistency of the u_bc, p_bc with the specified bc values..
    // this assumes semantically that the calcRhs of the boundary conditions
    // is called before the calcRhs of the volumetric cells... (sweep order)

    double ut[3];
    double ufr_n = DOT_PRODUCT(u_fr[ibf],unit_n);
    for (int i = 0; i < 3; ++i)
      ut[i] = u_fr[ibf][i] - ufr_n*unit_n[i];


    for (int i =0; i < 3; ++i)
      u_bc[ibf][i] = (rhoun_bc*unit_n[i] + rho_bc[ibf]*ut[i])/rho_bc[ibf];

    p_bc[ibf]      = rho_bc[ibf] * RT_bc;

    // and compute the rhs.. which is really only for rho...
    const int icv      = zone_ptr->cvobf[ibf];
    const double delta = bf_area/zone_ptr->area_over_delta_bf[ibf];
    calcNscbcIgSubsonicInletRhs(rhs_bc[ibf][rk_stage-1],unit_n,delta,
                                rho_bc[ibf],u_bc[ibf],p_bc[ibf],
                                solver->rho[icv],solver->u[icv],solver->p[icv],solver->gamma);

    // the coriolis force here is not relevant because there is only
    // 1 characteristic equation being solve (which is for rho)

  }
}

void NscbcOutletP::initData() {

  Nscbc::initData();

  Param * param = getParam(getName());
  if ( param->getString(0) == "NSCBC_OUTLET_P") {

    p_ref         = param->getDouble(1);
    sigma         = param->getDouble(2);
    L_ref         = param->getDouble(3);

    type          = OUTLET_PRESSURE;

  } else if ( param->getString(0) == "NSCBC_OUTLET_PRESSURE") {

    int iarg = 1;
    sigma    = 0.2;
    p_ref    = -1.0; // set to a negative value for error checking below..
    L_ref    = 1.0;  // XXX estimate a better value from the geometry here..

    while ( iarg < param->size()) {

      string token = param->getString(iarg++);
      if ( token == "SIGMA") {
        sigma = param->getDouble(iarg++);
      } else if ( token == "L") {
        L_ref = param->getDouble(iarg++);
      } else if ( token == "P_REF") {
        p_ref = param->getDouble(iarg++);
      } else if ( token == "CONSTANT_PRESSURE") {
        sigma = -1.0; // negative sigma value used to denote a constant pressure condition...
      } else {
        CWARN( " > skipping unrecognized token in NSCBC_OUTLET_PRESSURE: " << token);
      }
    }

    if ( p_ref < 0.0) {
      CERR( " > did not find P_REF specified for bc: " << zone_ptr->getName() );
    }

    type = OUTLET_PRESSURE;

  } else if ( param->getString(0) == "NSCBC_OUTLET_MDOT") {

    int iarg    = 1;
    rhou_target = -1.0;
    T_relax     = -1.0;

    while ( iarg < param->size()) {

      string token = param->getString(iarg++);
      if ( token == "T_RELAX") {
        T_relax = param->getDouble(iarg++);
      } else if ( token == "MDOT") {
        rhou_target = param->getDouble(iarg++);
      } else {
        CWARN( " > skipping unrecognized token in NSCBC_OUTLET_PRESSURE: " << token);
      }
    }

    if ( rhou_target < 0.0) {
      CERR( " > did not find MDOT specified for bc: " << zone_ptr->getName() );
    } else {
      rhou_target = rhou_target / zone_ptr->area_global;
    }

    if ( T_relax < 0.0) {
      CERR( " > did not find T_RELAX specified for bc: " << zone_ptr->getName());
    }

    type = OUTLET_MDOT_WEAK;

  } else {
    // if you get here there was a boundary condition parsing problem
    assert(0);
  }

  assert( mf == NULL); mf = new double[zone_ptr->nbf];

  // if reverse flow is encountered at the outlet, the passive scalars
  // will be reintroduced with 0 value as a default.

  assert( scbc_vals == NULL); scbc_vals = new double[solver->nsc];
  for (int isc = 0; isc < solver->nsc; ++isc)
    scbc_vals[isc] = 0.0;

}

void NscbcOutletP::calcRhsMdot(const double time, const int rk_stage) {

  // compute the present mdot ..

  double my_buf[3] = {0.0,0.0,0.0};

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const double bf_area = MAG(zone_ptr->n_bf[ibf]);
    const double un_dA   = DOT_PRODUCT(u_bc[ibf],zone_ptr->n_bf[ibf]);
    my_buf[0]           += bf_area;
    my_buf[1]           += rho_bc[ibf]*un_dA;
    my_buf[2]           += un_dA;

  }

  double buf[3];
  MPI_Allreduce(my_buf,buf,3,MPI_DOUBLE,MPI_SUM,mpi_comm);

  // by convention, the mdot that we are penalizing against is
  // positive for outgoing flow; same as the outward pointing normal here.

  double rhou_avg = buf[1]/buf[0];
  double u_avg    = buf[2]/buf[0];

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const double bf_area = MAG(zone_ptr->n_bf[ibf]);
    double unit_n[3]     = {zone_ptr->n_bf[ibf][0]/bf_area,
                            zone_ptr->n_bf[ibf][1]/bf_area,
                            zone_ptr->n_bf[ibf][2]/bf_area};

    const int icv      = zone_ptr->cvobf[ibf];
    const double delta = bf_area/zone_ptr->area_over_delta_bf[ibf];
    calcNscbcIgOutletMdotWeakRhs(rhs_bc[ibf][rk_stage-1],unit_n,delta,
                                 rho_bc[ibf],u_bc[ibf],p_bc[ibf],
                                 solver->rho[icv],solver->u[icv],solver->p[icv],solver->gamma,
                                 rhou_avg,u_avg,rhou_target, T_relax);
  }


  if ( solver->frame_rotation) {

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      // term 1: -omega x omega x r
      // centrifugal acceleration

      double coeff[3], rot_accel[3];
      coeff[0] = solver->frame_rotation[1]*zone_ptr->x_bf[ibf][2] - solver->frame_rotation[2]*zone_ptr->x_bf[ibf][1];
      coeff[1] = solver->frame_rotation[2]*zone_ptr->x_bf[ibf][0] - solver->frame_rotation[0]*zone_ptr->x_bf[ibf][2];
      coeff[2] = solver->frame_rotation[0]*zone_ptr->x_bf[ibf][1] - solver->frame_rotation[1]*zone_ptr->x_bf[ibf][0];

      rot_accel[0] = -(solver->frame_rotation[1]*coeff[2] - solver->frame_rotation[2]*coeff[1]);
      rot_accel[1] = -(solver->frame_rotation[2]*coeff[0] - solver->frame_rotation[0]*coeff[2]);
      rot_accel[2] = -(solver->frame_rotation[0]*coeff[1] - solver->frame_rotation[1]*coeff[0]);

      // term 2: -2 omega x v
      // Coriolis acceleration

      rot_accel[0] -= 2.0*(solver->frame_rotation[1]*u_bc[ibf][2] - solver->frame_rotation[2]*u_bc[ibf][1]);
      rot_accel[1] -= 2.0*(solver->frame_rotation[2]*u_bc[ibf][0] - solver->frame_rotation[0]*u_bc[ibf][2]);
      rot_accel[2] -= 2.0*(solver->frame_rotation[0]*u_bc[ibf][1] - solver->frame_rotation[1]*u_bc[ibf][0]);

      for (int i = 0; i <3; ++i)
        rhs_bc[ibf][rk_stage-1][2+i] += rot_accel[i];

     }

  }


}


void NscbcOutletP::calcRhsPressure(const double time, const int rk_stage) {

  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    const double bf_area = MAG(zone_ptr->n_bf[ibf]);
    double unit_n[3]     = {zone_ptr->n_bf[ibf][0]/bf_area,
                            zone_ptr->n_bf[ibf][1]/bf_area,
                            zone_ptr->n_bf[ibf][2]/bf_area};

    const int icv      = zone_ptr->cvobf[ibf];
    const double delta = bf_area/zone_ptr->area_over_delta_bf[ibf];
    calcNscbcIgOutletRhs(rhs_bc[ibf][rk_stage-1],unit_n,delta,
                         rho_bc[ibf],u_bc[ibf],p_bc[ibf],
                         solver->rho[icv],solver->u[icv],solver->p[icv],solver->gamma,
                         L_ref,sigma,p_ref);

  }


  if ( solver->frame_rotation) {

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      // term 1: -omega x omega x r
      // centrifugal acceleration

      double coeff[3], rot_accel[3];
      coeff[0] = solver->frame_rotation[1]*zone_ptr->x_bf[ibf][2] - solver->frame_rotation[2]*zone_ptr->x_bf[ibf][1];
      coeff[1] = solver->frame_rotation[2]*zone_ptr->x_bf[ibf][0] - solver->frame_rotation[0]*zone_ptr->x_bf[ibf][2];
      coeff[2] = solver->frame_rotation[0]*zone_ptr->x_bf[ibf][1] - solver->frame_rotation[1]*zone_ptr->x_bf[ibf][0];

      rot_accel[0] = -(solver->frame_rotation[1]*coeff[2] - solver->frame_rotation[2]*coeff[1]);
      rot_accel[1] = -(solver->frame_rotation[2]*coeff[0] - solver->frame_rotation[0]*coeff[2]);
      rot_accel[2] = -(solver->frame_rotation[0]*coeff[1] - solver->frame_rotation[1]*coeff[0]);

      // term 2: -2 omega x v
      // Coriolis acceleration

      rot_accel[0] -= 2.0*(solver->frame_rotation[1]*u_bc[ibf][2] - solver->frame_rotation[2]*u_bc[ibf][1]);
      rot_accel[1] -= 2.0*(solver->frame_rotation[2]*u_bc[ibf][0] - solver->frame_rotation[0]*u_bc[ibf][2]);
      rot_accel[2] -= 2.0*(solver->frame_rotation[0]*u_bc[ibf][1] - solver->frame_rotation[1]*u_bc[ibf][0]);

      for (int i = 0; i <3; ++i)
        rhs_bc[ibf][rk_stage-1][2+i] += rot_accel[i];

     }

  }
}

void NscbcOutletP::calcRhs(const double time, const int rk_stage) {

  if ( type == OUTLET_PRESSURE) {

    this->calcRhsPressure(time,rk_stage);

  } else {

    assert( type == OUTLET_MDOT_WEAK);

    this->calcRhsMdot(time, rk_stage);

  }
}

//
// sponge boundary condition implementation below.
// XXX the following implementation assumes that there is only one active sponge.
// otherwise the registration that the bc forces in the interior will run into
// name conflicts.  the bc registration is postponed until the initData
//

Sponge::Sponge(BfZone* p, IdealGasSolver* s) : IdealGasBc(p,s) {
  for (int i =0; i < 5; ++i)
    sponge_data[i] = 0.0;
  eps_p_sponge     = 0.0;
  t_relax          = 0.0;
  sponge_strength  = 0.0;
  sponge_length    = 0.0;
  sponge_type      = SPONGE_UNDEFINED;
  sponge_pressure  = -1.0;

  T_sponge  = NULL; s->registerCvData(T_sponge,"T_sponge",READWRITE_DATA);
  u_sponge = NULL; s->registerCvData(u_sponge,"u_sponge",READWRITE_DATA);
}

void Sponge::initData() {

  Param* param = getParam(getName());

  // default params..
  eps_p_sponge = 0.01; // not used any more!
  sponge_pressure = solver->p_ref;

  parseSpongeParam(param,sponge_type,sponge_strength,eps_p_sponge,
                   t_relax,sponge_length,sponge_pressure);

  T_sponge  = new double[solver->ncv];
  u_sponge = new double[solver->ncv][3];

  // use nodes to get the domain limits...
  double my_buf[5] = { -1.0E-20, -1.0E-20, -1.0E-20, -1.0E-20, -1.0E-20 };
  for (int ino = 0; ino < solver->nno; ++ino) {
    my_buf[0] = max(my_buf[0],solver->x_no[ino][0]);
    my_buf[1] = max(my_buf[1],solver->x_no[ino][1]);
    my_buf[2] = max(my_buf[2],solver->x_no[ino][2]);
    my_buf[3] = max(my_buf[3],sqrt(solver->x_no[ino][1]*solver->x_no[ino][1] + solver->x_no[ino][2]*solver->x_no[ino][2]));
    my_buf[4] = max(my_buf[4],sqrt(solver->x_no[ino][0]*solver->x_no[ino][0] + solver->x_no[ino][1]*solver->x_no[ino][1]));
  }
  MPI_Allreduce(my_buf,sponge_data,5,MPI_DOUBLE,MPI_MAX,mpi_comm);

  assert( mf == NULL); mf = new double[zone_ptr->nbf];

  // if reverse flow is encountered at the outlet, the passive scalars
  // will be reintroduced with 0 value as a default.

  assert( scbc_vals == NULL); scbc_vals = new double[solver->nsc];
  for (int isc = 0; isc < solver->nsc; ++isc)
    scbc_vals[isc] = 0.0;

}

int Sponge::storeData(double * buf) const {

  // sponge stores cell values for T and u

  if (buf) {
    int count = 0;
    for (int icv = 0; icv < solver->ncv; ++icv) {
      buf[count++] = T_sponge[icv];
      FOR_I3 buf[count++] = u_sponge[icv][i];
    }
    assert(count = solver->ncv*4);
  }

  return 4*solver->ncv;

}

int Sponge::restoreData(double * buf) {

  assert(buf);
  int count = 0;
  for (int icv = 0; icv < solver->ncv; ++icv) {
    T_sponge[icv] = buf[count++];
    FOR_I3 u_sponge[icv][i] = buf[count++];
  }
  assert(count = solver->ncv*4);

  return 4*solver->ncv;

}

Sponge::~Sponge() {
  DELETE(T_sponge);
  DELETE(u_sponge);
}

void Sponge::addBoundaryFlux(IdealGasRhs* rhs) const {

  assert( mf != NULL);

  const double gm1    = solver->gamma - 1.0;
  const double invgm1 = 1.0/gm1;
  const double Rgogm1 = solver->R_gas*solver->gamma*invgm1;

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv       = zone_ptr->cvobf[ibf];

    // calcRiemannFlux requires we define u,p,sp_vol, and h....
    IdealGasState bf_state;

    //recall T = p/(rho*R_gas)
    bf_state.sp_vol = T_sponge[icv]*solver->R_gas/sponge_pressure;
    for (int i =0; i < 3; ++i)
      bf_state.u[i] = u_sponge[icv][i];
    bf_state.p      = sponge_pressure;
    bf_state.h      = Rgogm1*T_sponge[icv];

    IdealGasRhs flux;
    calcRiemannFlux(flux,zone_ptr->n_bf[ibf],solver->cv_light[icv],
                    bf_state,solver->gamma,solver->gamma);

    mf[ibf]       = flux.rho;

    rhs[icv].rho -= flux.rho;
    for (int i =0; i < 3; ++i)
      rhs[icv].rhou[i] -= flux.rhou[i];
    rhs[icv].rhoE -= flux.rhoE;

  }

  switch( sponge_type) {
  case SPONGE_LX:
  case SPONGE_LY:
  case SPONGE_LZ:
    {

      int idir;
      if (sponge_type == SPONGE_LX) {
        idir = 0;
      } else if ( sponge_type == SPONGE_LY) {
        idir = 1;
      } else {
        assert( sponge_type == SPONGE_LZ);
        idir = 2;
      }

      const double Lmax = sponge_data[idir];

      for (int icv = 0; icv < solver->ncv; ++icv) {
        if ( solver->x_cv[icv][idir] > Lmax-sponge_length) {
          const double x_sp = (solver->x_cv[icv][idir] - (Lmax-sponge_length))/sponge_length;
          // sponge profile = a*x^2 + b*x^8..
          // const double sponge_coeff = sponge_strength*(0.068*x_sp*x_sp + 0.845*pow(x_sp,8));
          // quadratic sponge profile (see JCP 212 681-702 (2006) & JCP 231 704-716 (2012))
	  const double sponge_coeff = sponge_strength*x_sp*x_sp;
	  const double rho_sponge = sponge_pressure/(T_sponge[icv]*solver->R_gas);
          rhs[icv].rho             -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv] - rho_sponge);
          for (int i =0; i < 3; ++i)
            rhs[icv].rhou[i]       -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv]*solver->u[icv][i] - rho_sponge*u_sponge[icv][i]);
	  const double rhoE_sponge = sponge_pressure/gm1 + 0.5*rho_sponge*DOT_PRODUCT(u_sponge[icv],u_sponge[icv]);
          rhs[icv].rhoE            -= solver->vol_cv[icv]*sponge_coeff*(solver->rhoE[icv] - rhoE_sponge);
        }
      }

    }
    break;

  case SPONGE_LRX:
    {

      const double Lmax = sponge_data[3];
      for (int icv = 0; icv < solver->ncv; ++icv) {
        const double rx = sqrt( solver->x_cv[icv][1]*solver->x_cv[icv][1] +
                                solver->x_cv[icv][2]*solver->x_cv[icv][2]   );
        if ( rx > Lmax - sponge_length) {
          const double x_sp = (rx - (Lmax-sponge_length))/sponge_length;
          // const double sponge_coeff = sponge_strength*(0.068*x_sp*x_sp + 0.845*pow(x_sp,8));
          // quadratic sponge profile (see JCP 212 681-702 (2006) & JCP 231 704-716 (2012))
          const double sponge_coeff = sponge_strength*x_sp*x_sp;
	  const double rho_sponge = sponge_pressure/(T_sponge[icv]*solver->R_gas);
	  rhs[icv].rho             -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv] - rho_sponge);
          for (int i =0; i < 3; ++i)
            rhs[icv].rhou[i]       -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv]*solver->u[icv][i] - rho_sponge*u_sponge[icv][i]);
	  const double rhoE_sponge = sponge_pressure/gm1 + 0.5*rho_sponge*DOT_PRODUCT(u_sponge[icv],u_sponge[icv]);
	  rhs[icv].rhoE            -= solver->vol_cv[icv]*sponge_coeff*(solver->rhoE[icv] - rhoE_sponge);
        }
      }
    }
    break;

  case SPONGE_LRZ:
    {

      const double Lmax = sponge_data[4];
      for (int icv = 0; icv < solver->ncv; ++icv) {
        const double rz = sqrt ( solver->x_cv[icv][0]*solver->x_cv[icv][0] +
                                 solver->x_cv[icv][1]*solver->x_cv[icv][1]  );
        if ( rz > Lmax - sponge_length) {
          const double x_sp = (rz-(Lmax-sponge_length))/sponge_length;
          // const double sponge_coeff = sponge_strength*(0.068*x_sp*x_sp + 0.845*pow(x_sp,8));
    	  // quadratic sponge profile (see JCP 212 681-702 (2006) & JCP 231 704-716 (2012))
          const double sponge_coeff = sponge_strength*x_sp*x_sp;
	  const double rho_sponge = sponge_pressure/(T_sponge[icv]*solver->R_gas);
          rhs[icv].rho             -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv] - rho_sponge);
          for (int i =0; i < 3; ++i)
            rhs[icv].rhou[i]       -= solver->vol_cv[icv]*sponge_coeff*(solver->rho[icv]*solver->u[icv][i] - rho_sponge*u_sponge[icv][i]);
	  const double rhoE_sponge = sponge_pressure/gm1 + 0.5*rho_sponge*DOT_PRODUCT(u_sponge[icv],u_sponge[icv]);
          rhs[icv].rhoE            -= solver->vol_cv[icv]*sponge_coeff*(solver->rhoE[icv] - rhoE_sponge);
        }
      }
    }
    break;

  default:
    assert(0);  // these errors have been checked earlier..
  }
}

void Sponge::rk3Step(const double dt, const int rk_stage) {

  if ( rk_stage == 1) {

    // relax the state at the beginning of the step (solver
    // data is consistent at time level n.. )

    const double eps = 1.0/(t_relax/dt + 1.0);
    assert(dt > 0.0);
    assert(eps > 0.0);

    for (int icv = 0; icv < solver->ncv; ++icv) {
      // do we have T available? should have...
      T_sponge[icv]   += eps*(solver->T[icv] - T_sponge[icv]);
      for (int i =0; i < 3; ++i)
        u_sponge[icv][i] += eps*(solver->u[icv][i] - u_sponge[icv][i]);
    }

  }

}

void Sponge::initialHook() {

  if ( !(solver->checkDataFlag("T_sponge")) ||
       !(solver->checkDataFlag("u_sponge")) ||
       checkParam("RESET_SPONGE") ) {

    // get these from somewhere...
    // TODO: GB: how do you want to set this?
    /*
    double T_sponge_init = solver->T_ref;
    double u_sponge_init[3];
    u_sponge_init[0] = 0.009;
    u_sponge_init[1] = 0.0;
    u_sponge_init[2] = 0.0;

    if ( mpi_rank == 0)
      cout << " > Setting values for sponge to u,T: " << COUT_VEC(u_sponge_init) << " " << T_sponge_init << endl;

    for (int icv = 0; icv < solver->ncv; ++icv) {
      T_sponge[icv] = T_sponge_init;
      for (int i =0 ; i < 3; ++i)
        u_sponge[icv][i] = u_sponge_init[i];
    }
    */

    if ( mpi_rank == 0)
      cout << " > Setting sponge to current u and T" << endl;

    for (int icv = 0; icv < solver->ncv; ++icv) {
      T_sponge[icv] = solver->T[icv];
      for (int i =0 ; i < 3; ++i)
	u_sponge[icv][i] = solver->u[icv][i];
    }

  }
}

void Sponge::query(const string& param_str) const {

  // report integrated values for the p_bc, rho_bc, u_bc, T_bc
  double my_buf[5] = {0.0,0.0,0.0,0.0,0.0};
  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv           = zone_ptr->cvobf[ibf];
    const double area       = MAG(zone_ptr->n_bf[ibf]);
    const double rho_sponge = sponge_pressure/(T_sponge[icv]*solver->R_gas);
    const double rhoundA_bc = rho_sponge*DOT_PRODUCT(zone_ptr->n_bf[ibf],u_sponge[icv]);

    my_buf[0]        += area;
    my_buf[1]        += area*rho_sponge;
    my_buf[2]        += rhoundA_bc;
    my_buf[3]        += area*sponge_pressure;
    my_buf[4]        += area*T_sponge[icv];
  }

  double buf[5];
  MPI_Reduce(my_buf,buf,5,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  if ( mpi_rank == 0 ) {
    cout << " > " << zone_ptr->getName() << ": time, rho_bc, mdot, p_bc, T_bc: "
         << solver->time << "   " << buf[1]/buf[0] << "    "
         << buf[2] << "    " << buf[3]/buf[0] << "    " << buf[4]/buf[0] << endl;
  }

  MiscUtils::dumpRange(T_sponge,solver->ncv,"T_sponge");
  MiscUtils::dumpRange(u_sponge,solver->ncv,"u_sponge");

}

void Sponge::addScalarBoundaryFlux(double* rhs_sc) const {
  BcTemplate<IdealGasSolver,IdealGasRhs>::_addScalarBoundaryFlux(rhs_sc);
}

// This is the beginning of inflow turbulence routines
void InflowTurbulenceIG::initData() {

  um_it   = new double[zone_ptr->nbf][3];
  Rd_it   = new double[zone_ptr->nbf][3];
  Rod_it   = new double[zone_ptr->nbf][3];
  Tm_it= new double[zone_ptr->nbf];
  pm_it  = new double[zone_ptr->nbf];
  length_scale  = new double[zone_ptr->nbf];
  u_bc   = new double[zone_ptr->nbf][3];
  T_bc= new double[zone_ptr->nbf];
  p_bc  = new double[zone_ptr->nbf];

  avg   = new double[zone_ptr->nbf][3];
  rms   = new double[zone_ptr->nbf][3];

  //  This store the state at the boundary
  assert( bf_state == NULL);
  bf_state = new IdealGasState[zone_ptr->nbf];

  IturbWgt=0.0;
  max_length=1.0E+200;
  in_plane_length=1.0E+200;
  Un=1.0E+200;
  inflow_zero=0.001;
  inflow_maxiter=30;

  Param * param = getParam(getName());
  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<IdealGasSolver,IdealGasRhs>::parseScalarBc(param);

}

void InflowTurbulenceIG::initialHook() {

  // set the statistics of the turbulence...SHould this only happen at step 0 ? Or when the user asks?
  setStatistics();
  //assert(1==0);
  // we need the reduced connectivity of the boundary attached cvs. These are computed here
  buildITBoundaryStructures();

  //now we compute and store the lengthscale based on the distance to the wall and/or whatever other input the user has specified
  // same problem as set statistics
  setLengthscale();

  //create the operators
  buildFilterOperators();

  // This should only be done when we are at step 0....or if the user would like to reset the inflow turbulence...
  initiateStats();
}

int InflowTurbulenceIG::storeData(double * buf) const {

  // should only need to store avg, rms
  // other values are static or updated from solver state

  if (buf) {
    int count = 0;
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      FOR_I3 buf[count++] = avg[ibf][i];
      FOR_I3 buf[count++] = rms[ibf][i];
    }
    assert(count == zone_ptr->nbf*6);
  }

  return 6*zone_ptr->nbf;

}

int InflowTurbulenceIG::restoreData(double * buf) {

  assert(buf);
  int count = 0;
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    FOR_I3 avg[ibf][i] = buf[count++];
    FOR_I3 rms[ibf][i] = buf[count++];
  }
  assert(count == zone_ptr->nbf*6);

  return 6*zone_ptr->nbf;

}

void InflowTurbulenceIG::calcRhs(const double time, const int rk_stage) {
  if(rk_stage==1){
    updateInflowTurbulence(u_bc);
  }
}


WmAlgCht::WmAlgCht(BfZone* p, IdealGasSolver* s) : IdealGasBc(p,s) {

  if (mpi_rank == 0)
    cout << "WmAlgCht() for zone: " << p->getName() << endl;

  u1       = NULL; zone_ptr->registerBfData(u1,"u1",READWRITE_DATA); // saved for filtering (if used)
  tau_wall = NULL; zone_ptr->registerBfData(tau_wall,"tau_wall",READWRITE_DATA); // saved because it is an initial guess
  q_wall   = NULL; zone_ptr->registerBfData(q_wall,"q_wall",READWRITE_DATA); // saved for no reason
  k_eff    = NULL;

  // allow up to 3 layers.
  nlayer = 0;
  FOR_I3 T_layer[i] = NULL; // double * for the layer fluid-side temperatures: T_layer1, T_layer2, T_layer3
  FOR_I3 k_layer[i] = NULL; // SimpleFunc pointers for the k's...

  np = 0;
  Tp = NULL;
  qp = NULL;
  ipobf_i = NULL;
  ipobf_v = NULL;
  ipobf_wgt = NULL;

  // add to funcEvalList
  funcEvalList.push_back(p->getName()+":tau_wall()");
  funcEvalList.push_back(p->getName()+":y_plus()");
  funcEvalList.push_back(p->getName()+":q_wall()");
  funcEvalList.push_back(p->getName()+":T_wall()");
  funcEvalList.push_back(p->getName()+":moment()");
  funcEvalList.push_back(p->getName()+":force()");
  // set the load balance cost. Default for bc is 1
  // see note above
  /*
    if (mpi_rank == 0) cout << "XXXXX setting lb_cost to 5" << endl;
    assert(zone_ptr->lb_cost == 1);
    zone_ptr->lb_cost = 5;
    assert(p->lb_cost == 5);
  */

  // look for optional thin layer model associated with this bc...

  Param * param = getParam(getName());
  assert(param);

  vector<pair<string,double> > layerVec;
  int iarg = 1;
  while (iarg < param->size()) {
    string token = param->getString(iarg++);
    if (token == "LAYER") {
      token = param->getString(iarg++);
      const double delta = param->getDouble(iarg++);
      if (mpi_rank == 0) cout << " > adding LAYER " << token << " with thickness: " << delta << endl;
      layerVec.push_back(pair<string,double>(token,delta));
    }
    else {
      CERR("unrecognized WM_ALG_CHT param: " << token);
    }
  }

  if (!layerVec.empty()) {
    nlayer = layerVec.size();
    if (nlayer > 3) {
      CERR("number of LAYERs in WM_ALG_CHT limited to 3 for now");
    }
    // register a temperature associated with the fluid side of each layer...
    char vname[32];
    for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
      // register the data with read/write...
      // T_layer1,T_layer2...
      sprintf(vname,"T_layer%d",ilayer+1);
      zone_ptr->registerBfData(T_layer[ilayer],vname,READWRITE_DATA);
      // set the layer delta...
      delta_layer[ilayer] = layerVec[ilayer].second;
      // and k+layer is a SimpleFunc...
      assert(k_layer[ilayer] == NULL);
      Param * k_param = getParam(layerVec[ilayer].first+".K");
      if (k_param == NULL) {
        CERR("LAYER " << layerVec[ilayer].first << " requires param for thermal conductivity " << layerVec[ilayer].first << ".K");
      }
      k_layer[ilayer] = processSimpleFunc(k_param);
      assert(k_layer[ilayer] != NULL);
    }
  }

}

void WmAlgCht::initData() {

  if (mpi_rank == 0)
    cout << "WmAlgCht::initData() for zone " << zone_ptr->getName() << endl;

  // we should have CHT by now...
  if (solver->cht == NULL) {
    CERR("WM_ALG_CHT bcs require CHT");
  }

  assert(u1 == NULL);        u1 = new double[zone_ptr->nbf];
  assert(tau_wall == NULL);  tau_wall = new double[zone_ptr->nbf];
  assert(q_wall == NULL);    q_wall = new double[zone_ptr->nbf];
  assert(k_eff == NULL);     k_eff = new double[zone_ptr->nbf]; // not registered

  // set everything to HUGE_VAL...

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    u1[ibf] = HUGE_VAL;
    tau_wall[ibf] = HUGE_VAL;
    q_wall[ibf] = HUGE_VAL;
    k_eff[ibf] = HUGE_VAL;
  }

  // and allocate and set any layer temperatures to HUGE_VAL...
  assert(nlayer <= 3);
  for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
    assert(T_layer[ilayer] == NULL);
    T_layer[ilayer] = new double[zone_ptr->nbf];
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
      T_layer[ilayer][ibf] = HUGE_VAL;
  }

  int * ssp_flag = new int[solver->subSurface->nsp];
  for (int issp = 0; issp < solver->subSurface->nsp; ++issp)
    ssp_flag[issp] = -1;

  for (int ibf = zone_ptr->ibf_f; ibf <= zone_ptr->ibf_l; ++ibf) {
    for (int sob = solver->sspobf_i[ibf]; sob != solver->sspobf_i[ibf+1]; ++sob) {
      const int issp = solver->sspobf_v[sob];
      ssp_flag[issp] = 0;
    }
  }

  np = 0;
  for (int issp = 0; issp < solver->subSurface->nsp; ++issp)
    if (ssp_flag[issp] == 0)
      ssp_flag[issp] = np++;

  assert(Tp == NULL);
  Tp = new double[np];
  for (int ip = 0; ip < np; ++ip) Tp[ip] = HUGE_VAL; // put garbage here to ensure it is set before being used

  assert(qp == NULL);
  qp = new double[np];
  for (int ip = 0; ip < np; ++ip) qp[ip] = HUGE_VAL; // put garbage here to ensure it is set before being used

  double (*xp)[3] = new double[np][3];
  int ip = 0;
  for (int issp = 0; issp < solver->subSurface->nsp; ++issp) {
    if (ssp_flag[issp] >= 0) {
      FOR_I3 xp[ip][i] = solver->subSurface->xp[issp][i];
      const int bits = int(solver->subSurface->isp_global_and_bits[issp]>>52);
      if (bits)
        PeriodicData::periodicTranslate(xp[ip],1,BitUtils::flipPeriodicBits(bits));
      ++ip;
    }
  }
  assert(ip == np);

  // register this bc with the CHT registration...
  //solver->cht->registerBoundaryData(xp,Tp,qp,np,zone_ptr->index);
  solver->cht->registerBoundaryData(xp,Tp,qp,ssp_flag,np,solver->subSurface->nsp,zone_ptr->index);
  delete[] xp;

  // build the ipobf_i/v/wgt needed to move information to/from the Tp/qp arrays...
  assert(ipobf_i == NULL);
  assert(zone_ptr->ibf_l-zone_ptr->ibf_f+1 == zone_ptr->nbf);
  ipobf_i = new int[zone_ptr->nbf+1];
  ipobf_i[0] = 0;
  const int ipobf_s = solver->sspobf_i[zone_ptr->ibf_l+1] - solver->sspobf_i[zone_ptr->ibf_f];
  assert(ipobf_v == NULL);
  ipobf_v = new int[ipobf_s];
  assert(ipobf_wgt == NULL);
  ipobf_wgt = new double[ipobf_s];
  int iob = 0;
  for (int ibf = zone_ptr->ibf_f; ibf <= zone_ptr->ibf_l; ++ibf) {
    double wgt_sum = 0.0;
    const int iob0 = iob;
    bool b_unmatched = false;
    for (int sob = solver->sspobf_i[ibf]; sob != solver->sspobf_i[ibf+1]; ++sob) {
      const int issp = solver->sspobf_v[sob];
      //assert((ssp_flag[issp] >= 0)&&(ssp_flag[issp] < np));
      if (ssp_flag[issp] >= 0) {
        assert(ssp_flag[issp] < np);
        ipobf_v[iob] = ssp_flag[issp];
        ipobf_wgt[iob] = solver->sspobf_wgt[sob];
        wgt_sum += ipobf_wgt[iob];
        ++iob;
      }
      else {
        b_unmatched = true;
      }
    }
    //assert(fabs(1.0-wgt_sum) < 1.0E-12);
    ipobf_i[ibf-zone_ptr->ibf_f+1] = iob;
    if (b_unmatched) {
      // redistribute unmatched ip weights...
      assert(wgt_sum > 1.0E-12);
      const double inv_wgt_sum = 1.0/wgt_sum;
      for (int iob_ = iob0; iob_ < iob; ++iob_)
        ipobf_wgt[iob_] *= inv_wgt_sum;
    }
  }
  //assert(iob == ipobf_s);
  assert(iob <= ipobf_s);
  delete[] ssp_flag;

}

void WmAlgCht::initialHook() {

  // This is a better way to set boundary data. If we find any
  // single HUGE_VAL in any of the boundary data that should have been
  // read from the restart, we reset everything...

  bool b_init = false;
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    if (u1[ibf] == HUGE_VAL) {
      b_init = true;
      break;
    }
    if (tau_wall[ibf] == HUGE_VAL) {
      b_init = true;
      break;
    }
    if (q_wall[ibf] == HUGE_VAL) {
      b_init = true;
      break;
    }
    for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
      if (T_layer[ilayer][ibf] == HUGE_VAL) {
        b_init = true;
        break;
      }
    }
    if (b_init)
      break;
  }

  // --------------------------------------------------------------------
  // if even one stored value was HUGE_VAL, initialize the whole bc...
  // --------------------------------------------------------------------

  if (b_init) {

    if (mpi_rank == 0)
      cout << "WM_ALG_CHT " << zone_ptr->getName() << ": initializing all data" << endl;

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      const int icv = zone_ptr->cvobf[ibf];
      u1[ibf] = MAG(solver->u[icv]);
      // set an initial value for tau_wall, this will constitute the
      // initial guess for the solution when it is finally constructed.
      if (zone_ptr->area_over_delta_bf[ibf] > 0.0) {
        const double y1 = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        tau_wall[ibf] = solver->mu_lam[icv]*MAG(solver->u[icv])/y1;
      }
      else {
        tau_wall[ibf] = 0.0;
      }
      q_wall[ibf] = 0.0;
      k_eff[ibf] = 0.0;
      for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
        T_layer[ilayer][ibf] = solver->T[icv]; // initialize any layer temps to the fluid T
      }
    }

  }

}

int WmAlgCht::storeData(double * buf) const {

  // the wall model needs to store tau_wall because it is used as
  // an initial guess. No need to store u1. It is not used...
  // q_wall is computed. If there are any layer temperatures, these
  // need to be stored too...

  // TODO: if you reintroduce filtering, we need u1 here...

  if (buf) {
    int count = 0;
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      buf[count++] = tau_wall[ibf];
    }
    for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
      for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
        buf[count++] = T_layer[ilayer][ibf];
      }
    }
    assert(count == zone_ptr->nbf*(1+nlayer));
  }

  return zone_ptr->nbf*(1+nlayer);

}

int WmAlgCht::restoreData(double * buf) {

  assert(buf);
  int count = 0;
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    tau_wall[ibf] = buf[count++];
  }
  for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
      T_layer[ilayer][ibf] = buf[count++];
    }
  }
  assert(count == zone_ptr->nbf*(1+nlayer));

  return zone_ptr->nbf*(1+nlayer);

}

WmAlgCht::~WmAlgCht() {

  DELETE(u1);
  DELETE(tau_wall);
  DELETE(q_wall);
  DELETE(k_eff);

  FOR_I3 DELETE(T_layer[i]);
  FOR_I3 if (k_layer[i]) delete k_layer[i];

  DELETE(Tp);
  DELETE(qp);
  DELETE(ipobf_i);
  DELETE(ipobf_v);
  DELETE(ipobf_wgt);

}

void WmAlgCht::calcRhs(const double time, const int rk_stage) {

  // see notes left in IdealGasSolverBcs.cpp for the time
  // filtering approach

  assert((rk_stage >= 1)&&(rk_stage <= 3));
  if (rk_stage == 1) {

    // on the first stage, do the wall model calcs...

    double thermal_resistance[3];
    const int nip = 7; // number of interpolation points for T-dependent k in each layer
    double Tip[nip];
    double kip[nip];

    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

      const int icv = zone_ptr->cvobf[ibf];

      if (zone_ptr->area_over_delta_bf[ibf] > 0.0) {

        const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
        const double mag_u = MAG(solver->u[icv]);
        //const double eps   = mag_u*solver->dt / (y1 + solver->dt*mag_u);

        // \phi_new = eps*u + (1.0 -eps)*\phi_old
        // \phi_new - \phi_old = eps*(u - \phi_old)

        //u1[ibf]   += eps* ( mag_u - u1[ibf]);
        u1[ibf]   = mag_u;

        // use the existing tau_wall to set a guess for the u_tau

        const double u_tau = sqrt(tau_wall[ibf] / solver->rho[icv]);
        tau_wall[ibf] = AlgebraicWM::solve_tau(u1[ibf], y1, solver->rho[icv],
					       solver->mu_lam[icv], u_tau);

        // approximate material properties from the interior -- these are low Ma approximations
        // that ignore property changes in the near wall cell due to compressibility

        const double cp = solver->R_gas*solver->gamma/(solver->gamma-1.0);

        // reconstruct the wall temperature from the Tp data that should
        // be current in the Tp[ip] array. (updated with an earlier call to
        // cht->updateChtTemperatureToFlowBcs())...
        double T_wall = 0.0;
        for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
          const int ip = ipobf_v[iob];
          const double wgt = ipobf_wgt[iob];
          T_wall += wgt*Tp[ip];
        }

        // recall that flux and thermal resistance are related by:
        //
        // q = deltaT / sum(R)
        //
        // where sum(R) = dx1/k1 + dx2/k2...

        // determine the thermal resistance of each layer and sum it...
        double thermal_resistance_sum = 0.0;
        double T_layer_prev = T_wall; // start at the cht wall
        for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
          // interpolate the temperatures to the nip points inside the layer...
          for (int ip = 0; ip < nip; ++ip) {
            const double w1 = (double(ip)+0.5)/double(nip);
            Tip[ip] = (1.0-w1)*T_layer_prev + w1*T_layer[ilayer][ibf];
          }
          k_layer[ilayer]->eval(kip,Tip,nip);
          thermal_resistance[ilayer] = 0.0;
          for (int ip = 0; ip < nip; ++ip) {
            thermal_resistance[ilayer] += delta_layer[ilayer]/(kip[ip]*double(nip));
          }
          thermal_resistance_sum += thermal_resistance[ilayer];
          T_layer_prev = T_layer[ilayer][ibf];
        }

        // now call the wall model routine that includes a thermal resistance...

        q_wall[ibf] = AlgebraicWM::compute_q_wall_approx_r(solver->T[icv], T_wall, thermal_resistance_sum, tau_wall[ibf],
                                                           solver->rho[icv], solver->mu_lam[icv], solver->loc_lam[icv], cp, y1, false);

        // reset the layer temperatures based on thie new flux and the current T_wall...
        T_layer_prev = T_wall; // start at the cht wall
        for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
          T_layer[ilayer][ibf] = T_layer_prev + q_wall[ibf]*thermal_resistance[ilayer];
          T_layer_prev = T_layer[ilayer][ibf];
        }

        // and compute the effective resistance of the wall model + thermal resistance of the layers...
        // and save as k_eff (1/R) including the bf area...

        const double dT = solver->T[icv] - T_wall;

        if (dT != 0.0) {
          k_eff[ibf] = zone_ptr->area_bf[ibf]*q_wall[ibf]/dT;
        }
	else {
	  // just use the laminar value...
	  k_eff[ibf] = zone_ptr->area_bf[ibf]/(y1/(solver->loc_lam[icv]*cp) + thermal_resistance_sum);
	}

      }
      else {

        // this should be a collapsed face.
        u1[ibf] = MAG(solver->u[icv]);
        tau_wall[ibf] = 0.0;
        q_wall[ibf] = 0.0;
	k_eff[ibf] = 0.0;

      }
    }

  }

}

void WmAlgCht::addBoundaryFlux(IdealGasRhs* rhs) const {

  assert( solver != NULL);

  // zero the qp so it can be summed below...
  for (int ip = 0; ip < np; ++ip) qp[ip] = 0.0;

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    const int icv = zone_ptr->cvobf[ibf];

    // tau_wall has already been computed at this point, but we'll assume
    // that the wall stress is instantaneously aligned with the velocity vector

    double u_mag = 0.0;
    for (int i =0; i < 3; ++i)
      u_mag += solver->u[icv][i] * solver->u[icv][i];
    u_mag = sqrt(max(0.0,u_mag));

    if ( u_mag > 0.0) {
      for (int i =0; i < 3;++i) {
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];
        rhs[icv].rhou[i] -= tau_wall[ibf]*solver->u[icv][i]/u_mag*zone_ptr->area_bf[ibf];
      }
    }
    else {
      for (int i =0; i < 3; ++i)
        rhs[icv].rhou[i] -= solver->p[icv]*zone_ptr->n_bf[ibf][i];
    }

    double T_wall = 0.0;
    for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
      const int ip = ipobf_v[iob];
      const double wgt = ipobf_wgt[iob];
      T_wall += wgt*Tp[ip];
    }

    // XXX if u_bc \neq 0, there should be a u \dot \tau term added here..
    //rhs[icv].rhoE -= q_wall[ibf]*zone_ptr->area_bf[ibf];
    rhs[icv].rhoE -= k_eff[ibf]*(solver->T[icv] - T_wall);

    // scatter/sum q_bf to the qp[ip] data points for exchange to the FEM solve this sub-step...
    for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
      const int ip = ipobf_v[iob];
      const double wgt = ipobf_wgt[iob];
      qp[ip] += wgt*k_eff[ibf]*(solver->T[icv] - Tp[ip]); // should match -q_bf when summed over ip
    }

  }

}

CtiRegister::CtiDataError WmAlgCht::funcEvalCtiData(CtiRegister::CtiData& v,
    const string& name, list<CtiRegister::CtiData>& args, const bool b_eval_func) {

  using namespace CtiRegister;

  // all of the following functions do not take any arguments -- check first
  if ( zone_ptr->isLeadingStrMatched(name) ) {
    if (name == zone_ptr->getName()+":tau_wall") {

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {
        if (args.size() == 0) {
          // use what was computed by wall model
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
            tmp[ibf] = tau_wall[ibf];
          }
        }
        else if (args.size() == 1) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          double dir[3] = {0.0,0.0,0.0};

          if (arg->getType() == D3_DATA) {
            FOR_I3 dir[i] = arg->d3(i);
            NORMALIZE(dir);
          }
          else if (arg->getType() == D_DATA) {
            const int index = int(arg->d());
            if ((index < 0)||(index > 2)) return CTI_DATA_NOT_VALID;
            dir[index] = 1.0;
          }
          else {
            return CTI_DATA_NOT_VALID;
          }

          // loop faces and compute
          for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
            const int icv      = zone_ptr->cvobf[ibf];
            const double mag_u = MAG(solver->u[icv]);
            double unit_u[3];
            if ( mag_u > 0.0) {
              for (int i = 0; i < 3; ++i) unit_u[i] = solver->u[icv][i]/mag_u;
              tmp[ibf] = tau_wall[ibf]*DOT_PRODUCT(dir,unit_u);
            }
            else {
              tmp[ibf] = 0.0;
            }
          }
        }
        else {
          return CTI_DATA_ARG_COUNT;
        }
      }

      return CTI_DATA_OK;

    }
    else if (name == zone_ptr->getName()+":y_plus") {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          const int icv      = zone_ptr->cvobf[ibf];
          const double y1    = zone_ptr->area_bf[ibf] / zone_ptr->area_over_delta_bf[ibf];
          tmp[ibf]           = y1*sqrt(solver->rho[icv]*tau_wall[ibf])/solver->mu_lam[icv];
        }
      }

      return CTI_DATA_OK;

    }
    else if (name == zone_ptr->getName()+":q_wall") {

      if ( args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if ( b_eval_func) {
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          if (zone_ptr->area_bf[ibf] > 0.0) {
            double T_wall = 0.0;
            for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
              const int ip = ipobf_v[iob];
              const double wgt = ipobf_wgt[iob];
              T_wall += wgt*Tp[ip];
            }
            const int icv = zone_ptr->cvobf[ibf];
            // recall the definition of k_eff:
            //k_eff[ibf] = q_wall[ibf]*zone_ptr->area_bf[ibf]/dT;
            tmp[ibf] = k_eff[ibf]*(solver->T[icv] - T_wall)/zone_ptr->area_bf[ibf];
          }
          else {
            tmp[ibf] = 0.0;
          }
        }
      }

      return CTI_DATA_OK;

    }
    else if (name == zone_ptr->getName()+":T_wall") {

      if (args.size() != 0)
        return CTI_DATA_ARG_COUNT;

      double * tmp = zone_ptr->createBfD1Data(v);

      if (b_eval_func) {
        for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
          // reconstruct wall temperature from Tp...
          double T_wall = 0.0;
          for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
            const int ip = ipobf_v[iob];
            const double wgt = ipobf_wgt[iob];
            T_wall += wgt*Tp[ip];
          }
          tmp[ibf] = T_wall;
        }
      }

      return CTI_DATA_OK;

    }
    else if ( name == zone_ptr->getName()+":force") {

          if (b_eval_func) {
            if (args.size() == 0) {

              // pressure, viscous, convective force components
              double (*f_bf)[9] = new double[zone_ptr->nbf][9];
              force_bf(f_bf);

              computeForceValues(v,zone_ptr,f_bf,0,"total");

              delete[] f_bf;
            }
            else if (args.size() == 1) {
              list<CtiRegister::CtiData>::iterator arg = args.begin();

              int index = -1;
              if (arg->getType() == D_DATA) {
                index = int(arg->d());
                if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;

                double (*f_bf)[9] = new double[zone_ptr->nbf][9];
                force_bf(f_bf);

                computeForceValues(v,zone_ptr,f_bf,(index+1),"total");

                delete[] f_bf;
              }
              else if (arg->getType() == S_DATA) {

                // convective, pressure, viscous force components
                if ((arg->s() != "conv") && (arg->s() != "visc") && (arg->s() != "press")) return CTI_DATA_NOT_VALID;

                double (*f_bf)[9] = new double[zone_ptr->nbf][9];
                force_bf(f_bf);

                computeForceValues(v,zone_ptr,f_bf,0,arg->s());

                delete[] f_bf;
              }
              else return CTI_DATA_NOT_VALID;
            }
            else if (args.size() == 2) {
              list<CtiRegister::CtiData>::iterator arg = args.begin();

              if (arg->getType() != S_DATA) return CTI_DATA_NOT_VALID;
              const string component = arg->s();
              arg++;

              int index = -1;
              if (arg->getType() != D_DATA) return CTI_DATA_NOT_VALID;
              index = int(arg->d());
              if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;

              if ((component != "conv") && (component != "visc") && (component != "press")) return CTI_DATA_NOT_VALID;

              double (*f_bf)[9] = new double[zone_ptr->nbf][9];
              force_bf(f_bf);

              computeForceValues(v,zone_ptr,f_bf,(index+1),component);

              delete[] f_bf;
            }
            else return CTI_DATA_ARG_COUNT;
          }

          return CTI_DATA_OK;

    }
    else if ( name == zone_ptr->getName()+":moment") {
      if (b_eval_func) {
        double origin[3] = {0.0,0.0,0.0};
        if (args.size() == 0) {

          // convective, pressure, viscous moment components
          double (*f_bf)[9] = new double[zone_ptr->nbf][9];
          force_bf(f_bf);

          computeMomentValues(v,zone_ptr,f_bf,origin,0,"total");

          delete[] f_bf;
        }
        else if (args.size() == 1) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          if (arg->getType() == D_DATA) {
            int index = int(arg->d());
            if (index < 0 || index > 2) return CTI_DATA_NOT_VALID;

            double (*f_bf)[9] = new double[zone_ptr->nbf][9];
            force_bf(f_bf);

            computeMomentValues(v,zone_ptr,f_bf,origin,(index+1),"total");

            delete[] f_bf;
          }
          else if (arg->getType() == S_DATA) {

            // convective, pressure, viscous force components
            if ((arg->s() != "conv") && (arg->s() != "visc") && (arg->s() != "press")) return CTI_DATA_NOT_VALID;

            double (*f_bf)[9] = new double[zone_ptr->nbf][9];
            force_bf(f_bf);

            computeMomentValues(v,zone_ptr,f_bf,origin,0,arg->s());

            delete[] f_bf;
          }
          else return CTI_DATA_NOT_VALID;
        }
        else if (args.size() == 3) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          for (int i=0; i<3; ++i) {
            if (arg->getType() == D_DATA) {
              origin[i] = arg->d();
              arg++;
            }
            else return CTI_DATA_NOT_VALID;
          }

          double (*f_bf)[9] = new double[zone_ptr->nbf][9];
          force_bf(f_bf);

          computeMomentValues(v,zone_ptr,f_bf,origin,0,"total");

          delete[] f_bf;
        }
        else if (args.size() == 4) {
          list<CtiRegister::CtiData>::iterator arg = args.begin();

          if (arg->getType() != S_DATA) return CTI_DATA_NOT_VALID;
          const string component = arg->s();
          arg++;

          for (int i=0; i<3; ++i) {
            if (arg->getType() == D_DATA) {
              origin[i] = arg->d();
              arg++;
            }
            else return CTI_DATA_NOT_VALID;
          }

          if ((component != "conv") && (component != "visc") && (component != "press")) return CTI_DATA_NOT_VALID;

          double (*f_bf)[9] = new double[zone_ptr->nbf][9];
          force_bf(f_bf);

          computeMomentValues(v,zone_ptr,f_bf,origin,0,component);

          delete[] f_bf;
        }
        else return CTI_DATA_ARG_COUNT;
      }
      return CTI_DATA_OK;
    }
  }

  return IdealGasBc::funcEvalCtiData(v,name,args,b_eval_func);

}

void WmAlgCht::query(const string& param_str) const {

  // these statically sized buffers have room for up to 3 layers of temperature
  // data. Expand them if you implement more...

  double my_buf[9] = { 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0 };
  double my_minmax[10] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL,
                           HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL,
                           HUGE_VAL, HUGE_VAL };

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    if (zone_ptr->area_over_delta_bf[ibf] > 0.0) {

      const int icv = zone_ptr->cvobf[ibf];
      const double y1    = zone_ptr->area_bf[ibf]/zone_ptr->area_over_delta_bf[ibf];
      const double u_tau = sqrt( tau_wall[ibf]/solver->rho[icv]);
      const double nu    = solver->mu_lam[icv]/solver->rho[icv];
      const double y_plus = y1*u_tau/nu;

      my_buf[0] += zone_ptr->area_bf[ibf];
      my_buf[1] += zone_ptr->area_bf[ibf]*tau_wall[ibf];
      my_buf[2] += zone_ptr->area_bf[ibf]*y_plus;
      my_buf[3] += zone_ptr->area_bf[ibf]*q_wall[ibf];

      double T_wall = 0.0;
      for (int iob = ipobf_i[ibf]; iob != ipobf_i[ibf+1]; ++iob) {
        const int ip = ipobf_v[iob];
        const double wgt = ipobf_wgt[iob];
        T_wall += wgt*Tp[ip];
      }
      my_buf[4] += zone_ptr->area_bf[ibf]*T_wall;
      my_minmax[0] = min(my_minmax[0],T_wall);
      my_minmax[1] = min(my_minmax[1],-T_wall);

      my_buf[5] += zone_ptr->area_bf[ibf]*solver->T[icv];
      my_minmax[2] = min(my_minmax[2],solver->T[icv]);
      my_minmax[3] = min(my_minmax[3],-solver->T[icv]);

      for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
        my_buf[6+ilayer] += zone_ptr->area_bf[ibf]*T_layer[ilayer][ibf];
        my_minmax[4+ilayer*2] = min(my_minmax[4+ilayer*2],T_layer[ilayer][ibf]);
        my_minmax[4+ilayer*2+1] = min(my_minmax[4+ilayer*2+1],-T_layer[ilayer][ibf]);
      }

    }
  }

  double buf[9];
  MPI_Reduce(my_buf,buf,9,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  double minmax[10];
  MPI_Reduce(my_minmax,minmax,10,MPI_DOUBLE,MPI_MIN,0,mpi_comm);

  if ( mpi_rank == 0 ) {

    cout << " > " << zone_ptr->getName() << ": time, int(tau_wall)dA, int(q_w)dA, avg(q_w), avg(y_plus): "
         << solver->time << " "
         << buf[1] << " "
         << buf[3] << " "
         << buf[3]/buf[0] << " "
         << buf[2]/buf[0] << " T_wall(min,avg,max): "
         << minmax[0] << " " << buf[4]/buf[0] << " " << -minmax[1];

    for (int ilayer = 0; ilayer < nlayer; ++ilayer)
      cout << " T_layer" << ilayer+1 << ": " << minmax[4+ilayer*2] << " "
           << buf[6+ilayer]/buf[0] << " " << -minmax[4+ilayer*2+1];

    cout << " T[cvobf]: " << minmax[2] << " "
         << buf[5]/buf[0] << " " << -minmax[3] << endl;

  }

}

void CbcProfileSimple::initData() { 

  Param * param  = getParam(getName());
  assert( param != NULL);

  assert( param->getString(0) == "CBC_PROFILE_SIMPLE");
  assert( param->getString(1) == "FILE_NAME");

  oneDProfile * prof = new oneDProfile(); 
  prof->read(param->getString(2)); 

  assert( param->getString(3) == "P");
  double p_bc = param->getDouble(4);
         
  assert( param->getString(5) == "T");
  double T_bc = param->getDouble(6);

  assert( bf_state == NULL);
  bf_state = new IdealGasState[zone_ptr->nbf];

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {

    bf_state[ibf].u[0]   = prof->getData(zone_ptr->x_bf[ibf][1]);
    bf_state[ibf].u[1]   = 0.0;
    bf_state[ibf].u[2]   = 0.0;

    bf_state[ibf].p      = p_bc;
    bf_state[ibf].sp_vol = solver->R_gas*T_bc/p_bc ;
    bf_state[ibf].h      = solver->R_gas*T_bc*solver->gamma/(solver->gamma - 1.0);
  
  }


  assert( mf == NULL); mf = new double[zone_ptr->nbf];
  BcTemplate<IdealGasSolver,IdealGasRhs>::parseScalarBc(param);


  delete prof;

} 


