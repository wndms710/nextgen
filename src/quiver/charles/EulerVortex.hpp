#ifndef EULERVORTEX_HPP
#define EULERVORTEX_HPP

#include "Double3.hpp"

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
  const double e_twopi = 0.08; // normal
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

class CbcEulerVortex : public Cbc { 
public: 
  double Lx,Ly,u_inf;
  CbcEulerVortex(BfZone* p, IdealGasSolver* s) : Cbc(p,s) {
    Lx = getDoubleParam("EULER_LX",10.0);
    Ly = getDoubleParam("EULER_LY",10.0*sqrt(3.0)/2.0);
    u_inf = getDoubleParam("EULER_UINF",0.5);
  }
  ~CbcEulerVortex() {}
  void initData() {
    assert( mf == NULL); mf = new double[zone_ptr->nbf];
  }
  void addBoundaryFlux(IdealGasRhs* rhs) const { 
    assert( solver != NULL);
    const double gamma_ = solver->gamma;
    for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) { 
      const int icv = zone_ptr->cvobf[ibf];
      IdealGasRhs flux;
      double rho,u[3],rhoE;
      calcInviscidVortexLocal(rho,u,rhoE,zone_ptr->x_bf[ibf],1,solver->time,gamma_,solver->rho_ref,solver->p_ref,u_inf,false,Lx,Ly);
      IdealGasState bf_state;
      FOR_I3 bf_state.u[i] = u[i];
      bf_state.sp_vol  = 1.0/rho;
      bf_state.p       = (rhoE-0.5*rho*DOT_PRODUCT(u,u))*(gamma_-1.0);
      bf_state.h       = (rhoE*bf_state.sp_vol-0.5*DOT_PRODUCT(u,u))*gamma_;
      calcRiemannFlux(flux,zone_ptr->n_bf[ibf],solver->cv_light[icv],bf_state,gamma_, gamma_);

      mf[ibf]       = flux.rho;
      
      rhs[icv].rho -= flux.rho;
      for (int i =0; i < 3; ++i) 
        rhs[icv].rhou[i] -= flux.rhou[i];
      rhs[icv].rhoE -= flux.rhoE;

    }
  }

};

class EulerVortex : public IdealGasSolver { 
public:

  // declare your own data here...

  double * rho_exact;
  double Lx,Ly,u_inf;
  bool b_periodic;
 
  EulerVortex() {

    if (mpi_rank == 0)
      cout << "EulerVortex()" << endl;
    
    // register your own data here in the constructor if you 
    // want it read in from the sles file. Otherwise, you can 
    // register elsewhere (initData, for example). The bits at the 
    // end of the call can be:
    //   READWRITE_DATA    - for read/write access to sles file
    //   READ_DATA         - for read access to sles file
    //   WRITE_DATA        - for write acces to sles file
    //   NO_READWRITE_DATA - for no read/write access to sles file
    // Note the name specified in the register call is its identifier
    // (e.g., in the charles.in file or the Cascade app)...
    
    rho_exact = NULL;
    registerCvData(rho_exact,"rho_exact",READWRITE_DATA);
    
    Lx = getDoubleParam("EULER_LX",10.0);
    Ly = getDoubleParam("EULER_LY",10.0*sqrt(3.0)/2.0);
    b_periodic = getBoolParam("EULER_PERIODIC",true);
    u_inf = getDoubleParam("EULER_UINF",0.5);
  }

  virtual ~EulerVortex() {
  
    // clear your own data here in the destructor...
  
    DELETE(rho_exact);

  }
  
  void initData() {

    IdealGasSolver::initData();

    // allocate your own data here, if it needs be read
    // from the sles file. Otherwise, you can allocate it
    // later (like initialHook).
   
    rho_exact = new double[ncv];
    
    // HACK: if requested, override the current nef stuff...

    if (checkParam("NEW_OPS")) {

      if (mpi_rank == 0) cout << "XXXXXXXXXXXXXXXXXXXX NEW_OPS XXXXXXXXXXXXXXXXXXXXX" << endl;
      
      assert(nef > 0);
      assert(nef_i > 0);
      assert(n_ef); delete[] n_ef;
      assert(c_ef); delete[] c_ef;
      assert(cvoef); delete[] cvoef;
      assert(group_ef); delete[] group_ef; 
      
      //for (int ief = 0; ief < nef; ++ief) cout << "group_ef: " << group_ef[ief] << endl;

      int * no_flag = new int[nno];
      FOR_INO no_flag[ino] = -1;

      double * no_dist = new double[nno];
      FOR_INO no_dist[ino] = HUGE_VAL;
    
      int (*cvono)[4] = new int[nno][4];
      FOR_INO FOR_I4 cvono[ino][i] = -1;
      
      int my_ierr = 0;
      FOR_ICV {
        for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
          const int ifa = faocv_v[foc];
          for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
            const int ino = noofa_v[nof];
            if (no_flag[ino] != icv) {
              // advance to the first -1...
              int i;
              for (i = 0; i < 4; ++i) 
                if (cvono[ino][i] == -1)
                  break;
              if (i < 4) {
                cvono[ino][i] = icv;
              }
              else {
                cout << "WARNING: got degenerate node at: " << COUT_VEC(x_no[ino]) << endl;
                my_ierr = -1;
              }
              no_flag[ino] = icv;
              const double dist = DIST(x_vv[icv],x_no[ino]);
              if (no_dist[ino] == HUGE_VAL) {
                no_dist[ino] = dist;
              }
              else {
                const double ddist = no_dist[ino] - dist;
                assert(fabs(ddist) < 1.0E-14);
              }
            }
          }
        }
      }
      
      // now all nodes should have a meaningful length scale...
      
      FOR_INO assert(no_dist[ino] < HUGE_VAL);
      
      // complete the cvono groups with ghost data...
      
      // put all ghost cvs into a point adt...
      
      Adt<double> adt(ncv_g2-ncv,x_vv+ncv,x_vv+ncv);
      
      // now search adt for nearby cvs to the nodes...
    
      vector<int> intVec;
      FOR_INO { 
        assert(intVec.empty());
        adt.buildListForSphere(intVec,x_no[ino],no_dist[ino]+1.0E-12);
        for (int ii = 0; ii < intVec.size(); ++ii) {
          const int icv = intVec[ii]+ncv;
          const double dist = DIST(x_vv[icv],x_no[ino]);
          const double ddist = no_dist[ino] - dist;
          assert(fabs(ddist) < 1.0E-14);
          int i;
          for (i = 0; i < 4; ++i) 
            if (cvono[ino][i] == -1)
              break;
          if (i < 4) {
            cvono[ino][i] = icv;
          }
          else {
            cout << "ERROR: got degenerate node at: " << COUT_VEC(x_no[ino]) << endl;
            my_ierr = -1;
          }
        }
        intVec.clear();
      }

      int ierr;
      MPI_Allreduce(&my_ierr,&ierr,1,MPI_INT,MPI_MIN,mpi_comm);
      if (ierr == -1) {
        CERR("fix degenerate nodes");
      }
      
      FOR_INO FOR_I4 assert(cvono[ino][i] != -1);
      
      FOR_INO {
        const int icv0 = cvono[ino][0];
        const int icv1 = cvono[ino][1];
        const int icv2 = cvono[ino][2];
        const int icv3 = cvono[ino][3];
        double signed_tet_vol_6 = SIGNED_TET_VOLUME_6(x_vv[icv0],x_vv[icv1],x_vv[icv2],x_vv[icv3]);
        if (signed_tet_vol_6 < 0.0) {
          if ((icv0 < ncv)&&(icv1 < ncv)) {
            // if the first 2 are both internal, switch them to change volume sign...
            cvono[ino][1] = icv0;
            cvono[ino][0] = icv1;
            signed_tet_vol_6 = SIGNED_TET_VOLUME_6(x_vv[icv1],x_vv[icv0],x_vv[icv2],x_vv[icv3]);
          }
          else {
            // else switch the last 2 to change sign...
            cvono[ino][3] = icv2;
            cvono[ino][2] = icv3;
            signed_tet_vol_6 = SIGNED_TET_VOLUME_6(x_vv[icv0],x_vv[icv1],x_vv[icv3],x_vv[icv2]);
          }
        }
        assert(signed_tet_vol_6 > 0.0);
      }
    
      // now build the ops from the properly-oriented tets in cvono...
      
      double exact_grad[3] = { 1.2234123, 3.235235, -2.2463 };
      
      map<const pair<int,int>,Double3> n_ef_map;
      
      double *vol_vv = new double[ncv];
      FOR_ICV vol_vv[icv] = 0.0;
    
      double my_grad_err_max = 0.0;
      FOR_INO {
        const int icv0 = cvono[ino][0];
        const int icv1 = cvono[ino][1];
        const int icv2 = cvono[ino][2];
        const int icv3 = cvono[ino][3];
        const double signed_tet_vol_6 = SIGNED_TET_VOLUME_6(x_vv[icv0],x_vv[icv1],x_vv[icv2],x_vv[icv3]);
        assert(signed_tet_vol_6 > 0.0);
        // grad weights are as follows...
        double gl[4][3];
        {
          const double n[3] = TRI_NORMAL_2(x_vv[icv1],x_vv[icv2],x_vv[icv3]);
          FOR_I3 gl[0][i] = -n[i]/signed_tet_vol_6;
        }
        {
          const double n[3] = TRI_NORMAL_2(x_vv[icv0],x_vv[icv2],x_vv[icv3]);
          FOR_I3 gl[1][i] = n[i]/signed_tet_vol_6;
        }
        {
          const double n[3] = TRI_NORMAL_2(x_vv[icv0],x_vv[icv1],x_vv[icv3]);
          FOR_I3 gl[2][i] = -n[i]/signed_tet_vol_6;
        }
        {
          const double n[3] = TRI_NORMAL_2(x_vv[icv0],x_vv[icv1],x_vv[icv2]);
          FOR_I3 gl[3][i] = n[i]/signed_tet_vol_6;
        }
        // test...
        const double phi0 = DOT_PRODUCT(exact_grad,x_vv[icv0]);
        const double phi1 = DOT_PRODUCT(exact_grad,x_vv[icv1]);
        const double phi2 = DOT_PRODUCT(exact_grad,x_vv[icv2]);
        const double phi3 = DOT_PRODUCT(exact_grad,x_vv[icv3]);
        // check...
        double grad_check[3];
        FOR_I3 grad_check[i] = gl[0][i]*phi0 + gl[1][i]*phi1 + gl[2][i]*phi2 + gl[3][i]*phi3;
        FOR_I3 grad_check[i] = grad_check[i]/exact_grad[i] - 1.0;
        FOR_I3 my_grad_err_max = max(my_grad_err_max,fabs(grad_check[i]));
        // now assemble operators...
        // sum(1/2*(int(v[i]*grad[j])+int(v[j]*grad[i])),j)...
        FOR_I4 {
          FOR_J4 {
            // "i" determines the row. Only consider local row...
            if (cvono[ino][i] < ncv) {
              // lumped mass...
              vol_vv[cvono[ino][i]] += signed_tet_vol_6/24.0/4.0;
            }
            if ((cvono[ino][i] < ncv)&&(cvono[ino][j] > cvono[ino][i])) {
              double coeff[3]; FOR_K3 coeff[k] = signed_tet_vol_6/24.0*gl[j][k];
              // add coeff to ef normal in (i,j) order...
              map<const pair<int,int>,Double3>::iterator it = n_ef_map.find(pair<int,int>(cvono[ino][i],cvono[ino][j]));
              if (it == n_ef_map.end()) {
                n_ef_map[pair<int,int>(cvono[ino][i],cvono[ino][j])] = Double3(coeff);
              }
              else {
                FOR_K3 it->second.data[k] += coeff[k];
              }
            }
            else if ((cvono[ino][j] < ncv)&&(cvono[ino][j] < cvono[ino][i])) {
              // SUBTRACT coeff to ef normal in (j,i) order...
              double coeff[3]; FOR_K3 coeff[k] = -signed_tet_vol_6/24.0*gl[j][k];
              map<const pair<int,int>,Double3>::iterator it = n_ef_map.find(pair<int,int>(cvono[ino][j],cvono[ino][i]));
              if (it == n_ef_map.end()) {
                n_ef_map[pair<int,int>(cvono[ino][j],cvono[ino][i])] = Double3(coeff);
              }
              else {
                FOR_K3 it->second.data[k] += coeff[k];
              }
            }
          }
        }
      }

      {
        double my_buf[2] = { 0.0, 0.0 };
        assert(vol_cv);
        FOR_ICV {
          my_buf[0] += vol_cv[icv];
          my_buf[1] += vol_vv[icv];
        }
        double buf[2];
        MPI_Reduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
        if (mpi_rank == 0) cout << " compare vols: " << buf[0] << " " << buf[1] << endl;
      }
      
      double grad_err_max;
      MPI_Reduce(&my_grad_err_max,&grad_err_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      if (mpi_rank == 0)
        cout << " > tet grad_err_max: " << grad_err_max << endl;
      
      // now try the normals...
      
      double * phi = new double[ncv_g];
      FOR_ICV_G {
        phi[icv] = DOT_PRODUCT(exact_grad,x_vv[icv]);
      }
      
      double (*grad_phi)[3] = new double[ncv][3];
      FOR_ICV FOR_I3 grad_phi[icv][i] = 0.0;
      
      double (*gcl)[3] = new double[ncv][3];
      FOR_ICV FOR_I3 gcl[icv][i] = 0.0;
      
      nef_i = 0;
      nef = 0;
      for (map<const pair<int,int>,Double3>::iterator it = n_ef_map.begin(); it != n_ef_map.end(); ++it) {
        const int icv0 = it->first.first;
        assert((icv0 >= 0)&&(icv0 < ncv));
        const int icv1 = it->first.second;
        assert((icv1 >= 0)&&(icv1 < ncv_g));
        assert(icv0 < icv1);
        const double phi_ef = 0.5*(phi[icv0]+phi[icv1]);
        FOR_I3 grad_phi[icv0][i] += it->second.data[i]*phi_ef;
        FOR_I3 gcl[icv0][i] += it->second.data[i];
        if (icv1 < ncv) { 
          ++nef_i;
          FOR_I3 grad_phi[icv1][i] -= it->second.data[i]*phi_ef;
          FOR_I3 gcl[icv1][i] -= it->second.data[i];
        }
        else {
          ++nef;
        }
      }

      nef += nef_i;
      
      dumpRange(gcl,ncv,"GCL");
      
      delete[] gcl;

      n_ef = new double[nef][3];
      c_ef = new double[nef][3];
      cvoef = new int[nef][2];
      group_ef = new int[nef];
      
      int ief_i = 0;
      int ief = nef_i;
      for (map<const pair<int,int>,Double3>::iterator it = n_ef_map.begin(); it != n_ef_map.end(); ++it) {
        const int icv0 = it->first.first;
        assert((icv0 >= 0)&&(icv0 < ncv));
        const int icv1 = it->first.second;
        assert((icv1 >= 0)&&(icv1 < ncv_g));
        assert(icv0 < icv1);
        if (icv1 < ncv) { 
          FOR_I3 n_ef[ief_i][i] = it->second.data[i];
          FOR_I3 c_ef[ief_i][i] = 0.0;
          cvoef[ief_i][0] = icv0;
          cvoef[ief_i][1] = icv1;
          group_ef[ief_i] = -1;
          ++ief_i;
        }
        else {
          FOR_I3 n_ef[ief][i] = it->second.data[i];
          FOR_I3 c_ef[ief][i] = 0.0;
          cvoef[ief][0] = icv0;
          cvoef[ief][1] = icv1;
          group_ef[ief] = -1;
          ++ief;
        }
      }
      assert(ief_i == nef_i);
      assert(ief == nef);

      // overwrite vol_cv with vol_vv and update into ghosts...
      
      FOR_ICV {
        vol_cv[icv] = vol_vv[icv];
      }
      updateCvData(vol_cv);
      
      delete[] vol_vv;
      
      my_grad_err_max = 0.0;
      FOR_ICV {
        FOR_I3 grad_phi[icv][i] /= vol_cv[icv];
        // check relative to exact...
        //cout << " icv: " << icv << " grad_phi[icv]/exact_grad: " << 
        //  grad_phi[icv][0]/exact_grad[0]-1.0 << " " << 
        //  grad_phi[icv][1]/exact_grad[1]-1.0 << " " << 
        //  grad_phi[icv][2]/exact_grad[2]-1.0 << endl;
        //getchar();
        FOR_I3 my_grad_err_max = max(my_grad_err_max,fabs(grad_phi[icv][i]/exact_grad[i] - 1.0));
      }
      
      delete[] grad_phi;

      MPI_Reduce(&my_grad_err_max,&grad_err_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      if (mpi_rank == 0)
        cout << " > operator grad_err_max: " << grad_err_max << endl;
      
      
      
      MPI_Pause("ALL looks good");
      
    }
    
  }
  
  void initialHook() {

    // initialize base solver conservative data here (e.g., rho, u & rhoE for 
    // the IdealGasSolver). Note, that this is not necessary if you initialize 
    // from paramaters in the in file (INIT_RUP), or if you are restarting a
    // calculation from an sles file...

    if (step == 0) {
    
      if (mpi_rank == 0)
        cout << " > EulerVortex: initialHook(): " << endl;
      
      FOR_ICV calcInviscidVortexLocal(rho[icv],u[icv],rhoE[icv],x_cv[icv],1,0.0,gamma,rho_ref,p_ref,u_inf,b_periodic,Lx,Ly);
      
      temporalHook();
    }

  }

  void temporalHook() {

    // the temporal hook is called at the end of time step. it can be used to 
    // overwrite variables, calculate global metrics and error norms. it can also
    // be used to issue diagnostic output (tecplot, images, etc.)...
  
    //double * rho_error = new double[ncv];
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

        //const double r = MAG(x_cv[icv]);
        //if ( r > 30.0) continue;

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

    /*
      if ( step%500 == 0) { 
      // dump the rho error.. trying to see where we are going wrong...
      assert(mpi_size == 1);
      char filename[128]; sprintf(filename,"work.%04d.dat",step);
      FILE* fp = fopen(filename,"w");
      fprintf(fp,"#VARIABLES X Y Z RHO_ERROR");
      for (int icv = 0; icv < ncv; ++icv) {
      fprintf(fp,"%12.8g  %12.8g  %12.8g  %12.8g\n", 
      x_cv[icv][0], x_cv[icv][1], x_cv[icv][2], rho_error[icv]);
      }
      fclose(fp);
      }
    */
    //delete[] rho_error;
  }

  void finalHook() {

    // finalHook is called after a simulation run has completed (e.g., step == NSTEPS). It can be used
    // like the temporalHook. It is called prior to the final result is written, so any changes to 
    // written data will take effect...

  }

  IdealGasBc* initHookBc(BfZone* p) { 
    return new CbcEulerVortex(p,this);
  }

};

#endif
