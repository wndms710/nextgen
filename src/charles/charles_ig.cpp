
#include "FlowSolver.hpp"
#include "IdealGasSolver.hpp"

//===============================
// IdealGasSolver
//===============================

class MyIdealGasSolver : public IdealGasSolver {
public:

  /*
  double (* vort)[3];         // = curl(u);
  double * vort_mag;          // = |curl(u)|
  */

  MyIdealGasSolver() {

    /*
    vort     = NULL;  registerCvData( vort, "vort", READWRITE_DATA);
    vort_mag = NULL;  registerCvData( vort_mag, "vort_mag", NO_READWRITE_DATA);
    */
  }

  void initData() {

    IdealGasSolver::initData();

    /*
    assert( vort == NULL);     vort     = new double[ncv_g2][3];  // include vort calculation in ghosts...
    assert( vort_mag == NULL); vort_mag = new double[ncv];        // no ghosts ...
    */

  }

  ~MyIdealGasSolver() {

    /*
    DELETE(vort);
    DELETE(vort_mag);
    */

  }

  void initialHook() {}

  void temporalHook() {

    /*
    if ( step % check_interval == 0) {

      if ( mpi_rank == 0 )
        cout << " > computing vorticity ... " << endl;

      for (int icv = 0; icv < ncv; ++icv) {

        vort[icv][0]  = dudx[icv][2][1] - dudx[icv][1][2];
        vort[icv][1]  = dudx[icv][0][2] - dudx[icv][2][0];
        vort[icv][2]  = dudx[icv][1][0] - dudx[icv][0][1];

        vort_mag[icv] = MAG(vort[icv]);

      }

      updateCv2Data(vort, REPLACE_ROTATE_DATA);

    }
    */

  }

  void finalHook() {}

  IdealGasBc* initHookBc(BfZone* p) {
    CERR(" > user must supply a hook boundary condition for zone: " << p->getName());
  }

  void addSourceHook(IdealGasRhs* rhs, const double time, const int rk_stage) {}

};

// =====================================================
// Euler Vortex solver
//
// activate from the input file by EOS = IDEAL_GAS_EULER
// =====================================================

void calcInviscidVortexLocal(double &rho,double* u,double &rhoE,const double* x,const int n,
                             const double t,const double gamma,const double rho_inf,
                             const double p_inf,const double u_inf,const bool b_periodic,const double Lx,const double Ly) {
  
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
  const double theta = M_PI/3.0;
  //const double theta = M_PI/4.0;
  //const double theta = 0.0;
  const double cos_theta = cos(theta);
  const double sin_theta = sin(theta);
  const double Ma_inf = u_inf/sqrt(gamma*p_inf/rho_inf);
  const double rc = 1.0;

  // hack: uniform freestream...
  /*
  {
    rho = rho_inf;
    u[0] = u_inf*cos_theta;
    u[1] = u_inf*sin_theta;
    u[2] = 0.0;
    rhoE = p_inf/(gamma-1.0) + 0.5*rho*(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    return;
  }
  */
  
  // circulation parameter...
  //const double e_twopi = 0.001; // very weak
  //const double e_twopi = 0.005;
  //const double e_twopi = 0.08; // normal
  const double e_twopi = 0.32; // strong: matching cuda version
  //double e_twopi = 0.1;
  //double e_twopi = 0.4; //very strong 
  //double e_twopi = 1.0; //very very strong 

  // setup...
  const double coeff = 0.5 * e_twopi*e_twopi * (gamma-1.0) * Ma_inf*Ma_inf;
    
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

  const double f0 = 1.0 - (( dx*dx ) + ( dy*dy ))/( rc*rc );
  rho = rho_inf*pow( 1.0 - coeff * exp( f0 ) , 1.0/(gamma-1.0) );
  u[0] = u_inf*( cos_theta - e_twopi * ( dy )/rc * exp( f0 / 2.0 ) );
  u[1] = u_inf*( sin_theta + e_twopi * ( dx )/rc * exp( f0 / 2.0 ) );
  u[2] = 0.0;

  const double p = p_inf*pow( 1.0 - coeff * exp( f0 ) , gamma/(gamma-1.0) );
  rhoE = p/(gamma-1.0) + 0.5*rho*(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
  
}

class IdealGasEulerSolver : public IdealGasSolver {
public:
  
  // declare your own data here...
  
  double * rho_exact;
  double Lx,Ly,u_inf;
  bool b_periodic;
  
  IdealGasEulerSolver() {
    
    if (mpi_rank == 0)
      cout << "IdealGasEulerSolver()" << endl;
    
    rho_exact = NULL;
    registerCvData(rho_exact,"rho_exact",READWRITE_DATA);
    
    Lx = getDoubleParam("EULER_LX",10.0);
    Ly = getDoubleParam("EULER_LY",10.0*sqrt(3.0)/2.0);
    b_periodic = getBoolParam("EULER_PERIODIC",true);
    u_inf = getDoubleParam("EULER_UINF",0.5);

  }

  virtual ~IdealGasEulerSolver() {
  
    DELETE(rho_exact);

  }
  
  void initData() {

    IdealGasSolver::initData();

    assert(rho_exact == NULL);
    rho_exact = new double[ncv];
    
  }
  
  void initialHook() {

    if (step == 0) {
    
      if (mpi_rank == 0)
        cout << " > IdealGasEulerSolver: initialHook(): " << endl;
      
      FOR_ICV calcInviscidVortexLocal(rho[icv],u[icv],rhoE[icv],x_cv[icv],1,0.0,gamma,rho_ref,p_ref,u_inf,b_periodic,Lx,Ly);
      
      temporalHook();

    }

  }

  void temporalHook() {

    if (step%check_interval == 0) {

      if (mpi_rank == 0)
	cout << "Errors:" << endl;
      
      //double rho_exact;
      double u_exact[3];
      double rhoE_exact;
      
      double my_l2[6],my_linf[5];
      for (int i = 0; i < 6; i++)
	my_l2[i] = 0.0;
      for (int i = 0; i < 5; i++)
	my_linf[i] = 0.0;
      
      FOR_ICV {

        calcInviscidVortexLocal(rho_exact[icv],u_exact,rhoE_exact,x_cv[icv],1,time,gamma,rho_ref,p_ref,u_inf,b_periodic,Lx,Ly);
        
	double delta = rho[icv] - rho_exact[icv];
	my_l2[0]    += vol_cv[icv]*delta*delta;
	my_linf[0]   = max( fabs(delta), my_linf[0] );
	
	for (int i = 0; i < 3; i++) {
	  delta        = u[icv][i] - u_exact[i];
	  my_l2[i+1]  += vol_cv[icv]*delta*delta;
	  my_linf[i+1] = max( fabs(delta), my_linf[i+1] );
	}
	
	delta      = rhoE[icv] - rhoE_exact;
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
          sqrt(l2[4]/l2[5])};
        
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

};

int main(int argc, char* argv[]) {

  try {

    CTI_Init(argc,argv,"charles.in");

    {

      const bool b_post = checkParam("POST");

      if (Param * param = getParam("EOS")) {
        const string eos = param->getString();
        if ( eos == "IDEAL_GAS") {

          MyIdealGasSolver solver;

          if (b_post) {

            solver.initMin();
            solver.runPost();

          } else {

            solver.init();
            solver.run();
          }
        }
        else if ( eos == "IDEAL_GAS_EULER") {
	  
	  assert(!b_post);
	  IdealGasEulerSolver solver;
	  solver.init();
	  solver.run();
	  
	}
        else {
          CERR("unrecognized EOS: " << eos << ", possible choices are \n" <<
               "EOS IDEAL_GAS\n" <<
               "EOS IDEAL_GAS_LSP\n" << 
	       "EOS IDEAL_GAS_EULER");
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
