
#include "FlowSolver.hpp"
#include "HelmholtzSolver.hpp"

//===============================
// HelmholtzSolver
//===============================

class MyHelmholtzSolver : public HelmholtzSolver {
public:

  MyHelmholtzSolver() {}

  void initData() {

    HelmholtzSolver::initData();

  }

  ~MyHelmholtzSolver() {}

  void calcInviscidVortexLocal(double &rho,double* u,double &p,const double* x,
      const double t,const double rho_inf,const double u_inf,
      const bool b_periodic,const double Lx,const double Ly) {

    //const double Lx = 10.0;
    //const double Ly = 10.0*sqrt(3.0)/2.0;

    // for periodic grids, the xmin, xmax needs to be set... 
    const double xmin = -Lx;
    const double xmax =  Lx;
    const double ymin = -Ly;
    const double ymax =  Ly;

    // except the prism grids, which have (-5..5)*cos(30)...
    //const double xmin = -4.330127019;
    //const double xmax =  4.330127019;

    // x,y position of the vortex
    // centered...
    //const double x0 = 0.0; //-3.5;
    //const double y0 = 0.0; //-3.5;
    const double x0 = -2.5; 
    const double y0 = -2.5; 
    //const double x0 = -Lx/2.0; 
    //const double y0 = -Ly/2.0;

    // direction of propagation
    const double theta = M_PI/3.0; // default
    //const double theta = M_PI/4.0;
    //const double theta = 0.0;
    const double cos_theta = cos(theta);
    const double sin_theta = sin(theta);


    const double Ma_inf = u_inf/helmholtz_sos;
    const double rc = 1.0;

    // hack: uniform freestream...
    /*
       {
       rho = rho_inf;
       u[0] = u_inf*cos_theta;
       u[1] = u_inf*sin_theta;
       u[2] = 0.0;
       p = 0.0;
       return;
       }
       */

    // circulation parameter...
    //const double e_twopi = 0.001; // very weak
    //const double e_twopi = 0.005;
    const double e_twopi = 0.08; // normal
    //double e_twopi = 0.1;
    //double e_twopi = 0.4; //very strong 
    //double e_twopi = 1.0; //very very strong 

    // setup...
    const double coeff = 0.5 * e_twopi*e_twopi * Ma_inf*Ma_inf;

    double dx = x[0] - x0 - u_inf*cos_theta*t;
    double dy = x[1] - y0 - u_inf*sin_theta*t;

    // if periodic, shift the exact solution so it is aligned with
    // the charles calculation...
    if (b_periodic) {
      while(dx > xmax) dx -= (xmax-xmin);
      while(dx < xmin) dx += (xmax-xmin);
      while(dy > ymax) dy -= (ymax-ymin);
      while(dy < ymin) dy += (ymax-ymin);
    }

    const double f0 = - (( dx*dx ) + ( dy*dy ))/( rc*rc );
    rho = rho_inf*exp(-coeff*exp(f0));
    u[0] = u_inf*( cos_theta - e_twopi * ( dy )/rc * exp( f0 / 2.0 ) );
    u[1] = u_inf*( sin_theta + e_twopi * ( dx )/rc * exp( f0 / 2.0 ) );
    u[2] = 0.0;
    p = helmholtz_sos*helmholtz_sos*(rho-rho_inf);

  }

  void initialHook() {

    const string test_case = getStringParam("TEST_CASE","NONE");
    const double rho_ref = getDoubleParam("RHO");
    if (test_case == "TAYLOR_GREEN") {
      COUT1(" > running inviscid taylor green vortex...");
      const double sos2 = helmholtz_sos*helmholtz_sos;
      FOR_ICV {
        FOR_I3 u[icv][i] = 0.0;
        u[icv][0] =  cos(x_cv[icv][0])*sin(x_cv[icv][1]);
        u[icv][1] = -sin(x_cv[icv][0])*cos(x_cv[icv][1]);
        u[icv][2] = 0.0;
        rho[icv]  = rho_ref*exp(-0.5/sos2*(cos(x_cv[icv][0])*cos(x_cv[icv][0])+cos(x_cv[icv][1])*cos(x_cv[icv][1])));
        p[icv]    = sos2*(rho[icv]-rho_ref);
      }
    }
    else if (test_case == "ACOUSTIC_WAVE") {
      COUT1(" > running 1d linear acoustic wave...");
      const double sos2 = helmholtz_sos*helmholtz_sos;
      setDt(dt,0); // const dt only
      FOR_ICV {
        const double eps = 1.0E-4*sin(2.0*M_PI*x_vv[icv][0]);
        u[icv][0] = helmholtz_sos*eps;
        u[icv][1] = 0.0;
        u[icv][2] = 0.0;
        rho[icv]  = rho_ref*(1.0+eps);
        p[icv]    = sos2*(rho[icv]-rho_ref);
        // need previous time step's state for BDF2
        const double eps0 = 1.0E-4*sin(2.0*M_PI*(x_cv[icv][0]-helmholtz_sos*(-dt)));
        u0[icv][0] = helmholtz_sos*eps0;
        u0[icv][1] = 0.0;
        u0[icv][2] = 0.0;
        rho0[icv]  = rho_ref*(1.0+eps0);
      }
      setDataFlag("rho0",1);
      setDataFlag("u0",1);
    }
    else if (test_case == "EULER_VORTEX") {
      COUT1(" > running euler vortex...");
      const double Lx = getDoubleParam("EULER_LX",10.0);
      const double Ly = getDoubleParam("EULER_LY",10.0*sqrt(3.0)/2.0);
      const bool b_periodic = getBoolParam("EULER_PERIODIC",true);
      const double u_inf = getDoubleParam("EULER_UINF",0.5);
      setDt(dt,0); // const dt only
      double p0;
      FOR_ICV calcInviscidVortexLocal(rho0[icv],u0[icv],p0,x_cv[icv],-dt,
          rho_ref,u_inf,b_periodic,Lx,Ly);
      FOR_ICV calcInviscidVortexLocal(rho[icv],u[icv],p[icv],x_cv[icv],0.0,
          rho_ref,u_inf,b_periodic,Lx,Ly);
      setDataFlag("rho0",1);
      setDataFlag("u0",1);
    }
    setDataFlag("rho",1);
    setDataFlag("u",1);
    setDataFlag("p",1);

  }

  void temporalHook() {

    if (step%check_interval == 0) {
      const string test_case = getStringParam("TEST_CASE","NONE");
      const double rho_ref = getDoubleParam("RHO");
      if (test_case == "TAYLOR_GREEN") {
        double my_buf[7] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
        const double sos2 = helmholtz_sos*helmholtz_sos;
        FOR_ICV {
          double u_exact[3];
          u_exact[0] =  cos(x_cv[icv][0])*sin(x_cv[icv][1]);
          u_exact[1] = -sin(x_cv[icv][0])*cos(x_cv[icv][1]);
          u_exact[2] = 0.0;
          const double rho_exact = rho_ref*exp(-0.5/sos2*(cos(x_cv[icv][0])*cos(x_cv[icv][0])+cos(x_cv[icv][1])*cos(x_cv[icv][1])));
          const double p_exact = sos2*(rho_exact-rho_ref);
          FOR_I3 my_buf[i] += fabs(u[icv][i]-u_exact[i])*vol_cv[icv];
          my_buf[3] += fabs(rho[icv]-rho_exact)*vol_cv[icv];
          my_buf[4] += fabs(p[icv]-p_exact)*vol_cv[icv];
          my_buf[5] += 0.5*DOT_PRODUCT(u[icv],u[icv])*vol_cv[icv];
          my_buf[6] += vol_cv[icv];
        }
        double buf[7];
        MPI_Reduce(my_buf,buf,7,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
        if (mpi_rank == 0) {
          cout << " > time,E_u,E_v,E_w,E_rho,E_p,KE: " << time;
          FOR_I6 cout << " " << buf[i]/buf[6];
          cout << endl;
        }
      }
      else if (test_case == "ACOUSTIC_WAVE") {
        double my_buf[7] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
        const double sos2 = helmholtz_sos*helmholtz_sos;
        FOR_ICV {
          const double eps = 1.0E-4*sin(2.0*M_PI*(x_vv[icv][0]-helmholtz_sos*time));
          double u_exact[3];
          u_exact[0] = helmholtz_sos*eps;
          u_exact[1] = 0.0;
          u_exact[2] = 0.0;
          const double rho_exact = rho_ref*(1.0+eps);
          const double p_exact = sos2*(rho_exact-rho_ref);
          FOR_I3 my_buf[i] += fabs(u[icv][i]-u_exact[i])*vol_cv[icv];
          my_buf[3] += fabs(rho[icv]-rho_exact)*vol_cv[icv];
          my_buf[4] += fabs(p[icv]-p_exact)*vol_cv[icv];
          my_buf[5] += 0.5*DOT_PRODUCT(u[icv],u[icv])*vol_cv[icv];
          my_buf[6] += vol_cv[icv];
        }
        double buf[7];
        MPI_Reduce(my_buf,buf,7,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
        if (mpi_rank == 0) {
          cout << " > time,E_u,E_v,E_w,E_rho,E_p,KE: " << time;
          FOR_I6 cout << " " << buf[i]/buf[6];
          cout << endl;
        }
      }
      else if (test_case == "EULER_VORTEX") {
        COUT1(" > running euler vortex...");

        const double Lx = getDoubleParam("EULER_LX",10.0);
        const double Ly = getDoubleParam("EULER_LY",10.0*sqrt(3.0)/2.0);
        const bool b_periodic = getBoolParam("EULER_PERIODIC",true);
        const double u_inf = getDoubleParam("EULER_UINF",0.5);

        double rho_exact;
        double u_exact[3];
        double p_exact;

        double my_l2[6],my_linf[5];
        for (int i = 0; i < 6; i++)
          my_l2[i] = 0.0;
        for (int i = 0; i < 5; i++)
          my_linf[i] = 0.0;

        dumpRange(rho_uncorrected,ncv,"rho_uncorrected");
        FOR_ICV {

          calcInviscidVortexLocal(rho_exact,u_exact,p_exact,x_cv[icv],time,
              rho_ref,u_inf,b_periodic,Lx,Ly);

          double delta = rho[icv] - rho_exact;
          my_l2[0]    += vol_cv[icv]*delta*delta;
          my_linf[0]   = max( fabs(delta), my_linf[0] );

          for (int i = 0; i < 3; i++) {
            delta        = u[icv][i] - u_exact[i];
            my_l2[i+1]  += vol_cv[icv]*delta*delta;
            my_linf[i+1] = max( fabs(delta), my_linf[i+1] );
          }

          delta      = p[icv] - p_exact;
          my_l2[4]  += vol_cv[icv]*delta*delta;
          my_linf[4] = max( fabs(delta), my_linf[4] );
          my_l2[5]  += vol_cv[icv];

        }

        double l2[6],linf[5];
        MPI_Reduce(my_l2,l2,6,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
        MPI_Reduce(my_linf,linf,5,MPI_DOUBLE,MPI_MAX,0,mpi_comm);

        if (mpi_rank == 0) {

          const double l2s[5] = {
            sqrt(l2[0]/l2[5]), sqrt(l2[1]/l2[5]),
            sqrt(l2[2]/l2[5]), sqrt(l2[3]/l2[5]),
            sqrt(l2[4]/l2[5])
          };

          cout << " > time, Ltwo: " << time << " " << l2s[0] << " " <<
          l2s[1] << " " <<
          l2s[2] << " " <<
          l2s[3] << " " <<
          l2s[4] << endl;

          cout << " > time, Linf: " << time << " " <<
          linf[0] << " " <<
          linf[1] << " " <<
          linf[2] << " " <<
          linf[3] << " " <<
          linf[4] << endl;
        }

      }
    }

  }

  void finalHook() {}

  HelmholtzBc* initHookBc(BfZone* p) {
    CERR(" > user must supply a hook boundary condition for zone: " << p->getName());
  }


  // the Helmholtz solver has implicit time advancement in a fractional
  // step setting; as a result, the hooks for add source hooks are slightly
  // different.

  void momentumSourceHook(double * A,double (*rhs)[3]) {}
  void massSourceHook(double * rhs) {}

};

int main(int argc, char* argv[]) {

  try {

    CTI_Init(argc,argv,"charles.in");

    {

      const bool b_post = checkParam("POST");

      if (Param * param = getParam("EOS")) {
        const string eos = param->getString();
        if ( eos == "HELMHOLTZ") {

          MyHelmholtzSolver solver;

          if ( b_post) {

            solver.initMin();
            solver.runPost();

          } else {

            solver.init();
            solver.run();
          }
        }
        else {
          CERR("unrecognized EOS: " << eos << ", possible choices are \n" << "EOS HELMHOLTZ\n");
        }
      }
    }

    CTI_Finalize();
  }
  catch (int e) {
    if (e >= 0) {
      CTI_Finalize();
    }
    else {
      CTI_Abort();
    }
  }
  catch(...) {
    CTI_Abort();
  }

  return 0;

}
