#include "VofSolver.hpp"
#include "BasicPostpro.hpp"
 
class SurfaceWaveVofSolver: public VofSolver {

  void initialHook(){

    if (step == 0) {

      COUT1("SurfaceWaveVofSolver:initialHook()");

      const double init_center[2] = { 0.0, M_PI};
      const double init_amplitude = 0.01*2.0*M_PI;
      const double init_wavelength =  2.0*M_PI;

      FOR_ICV_G {
        g[icv] = -x_cv[icv][1] +  init_center[1] + init_amplitude*cos(2.0*M_PI*(x_cv[icv][0]-init_center[0])/init_wavelength);
        u[icv][0] = 0.0;
        u[icv][1] = 0.0;
        u[icv][2] = 0.0;
      }

      buildVofFromG();
      calcNormal();
      calcGfromVof();

    }
      
  }
 
  void temporalHook(){
    const double eps = 1.0E-7;
    if (mpi_rank == 0)
      cout << "temporalHook()" << endl;

    double myAmp=1E10;
    double myX=1E10;
    double X=1E10;
    double Amp=1E10;
    int my_icv = -1;
    FOR_ICV{
      if(abs(x_cv[icv][0]-M_PI)< myX ){
        myX=fabs(x_cv[icv][0]-M_PI);
      }
    }
    MPI_Allreduce(&myX,&X,1,MPI_DOUBLE,MPI_MIN,mpi_comm);

    double myliquid = 0.0;
    double myvol = 0.0;
    double vol = 0.0;
    double liquid = 0.0;
    
    FOR_ICV {
      if (fabs(x_cv[icv][0]-M_PI) <= X+ eps ) { 
        myliquid += vof[icv]*vol_cv[icv];
        myvol  += vol_cv[icv];
      }
    }
    
      MPI_Allreduce(&myliquid, &liquid, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
      MPI_Allreduce(&myvol, &vol, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

      //if (mpi_rank == 0 ) cout << "my volume = " << vol << endl;
      //if (mpi_rank == 0 ) cout << "my liquid = " << liquid << " " << liquid/vol << endl;
      
      //double delta = 2.0*M_PI/32.0;
      double delta = pow(vol_cv[0],1.0/3.0); // TODO may be wrong
      
      Amp = liquid;
      
      if (mpi_rank==0){
        cout << "amplitude: "<<" "<< time << " " << Amp <<endl;
      }
  }

  /*
  void temporalHook() {
    const double eps = 1.0E-8;
    if (mpi_rank == 0)
      cout << "Prosperetti::temporalHook()" << endl;
    
    ofstream myfile;
    myfile.open("amplitude.dat",ios::app);
    double myAmp=1E10;
    double myX=1E10;
    double X=1E10;
    double Amp=1E10;
    int my_icv = -1;
    FOR_ICV{
      if(cv_flag[icv] >= 1){
        if(abs(x_cv[icv][0])<myX && abs(g[icv])<0.196){
          myX=abs(x_cv[icv][0]);
          myAmp=x_cv[icv][1]+g[icv];
          my_icv = icv;
        }
      }
    }
    MPI_Allreduce(&myX,&X,1,MPI_DOUBLE,MPI_MIN,mpi_comm);
    if(myX!=X){			
      myAmp=1E10;
    }
    MPI_Reduce(&myAmp,&Amp,1,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
    if (mpi_rank==0){
      myfile<<time<<" "<< dt << " " << (Amp-3.141592)/2.0/3.141592 <<endl;
    }
    myfile.close();

  }
  */
};

class StaticDropVofSolver : public VofSolver {
public:

  void initialHook() {

    COUT1("StaticDropSolver::initialHook()");
    
    if (!checkDataFlag("vof")) {
      assert(!checkDataFlag("u"));

      // drop of vof = 1 at (0,1), r = 0.5
      // plus a liquid at the bottom up to y = -1...

      double x[4][3];
      FOR_ICV {
        FOR_I3 u[icv][i] = 0.0;
        vof[icv] = 0.0;
      }
   
    }
  }

  void temporalHook() {

    // velocity measures...
    double umax = 0.0;
    double utot = 0.0;
  
    // pressure measures...

    double pmin = +HUGE_VAL;
    double pmax = -HUGE_VAL;
    double pin_tot  = 0.0;
    double pout_tot = 0.0;
    double pin_part  = 0.0;
    double pout_part = 0.0;

    // volume measures...

    double vol = 0.0;
    double vin_tot  = 0.0;
    double vout_tot = 0.0;
    double vin_part  = 0.0;
    double vout_part = 0.0;

    // drop stuff...

    const double x_drop = 0.0;
    const double y_drop = 0.0;
    //const double z_drop = 0.0;
    const double r_drop = 2.0;
    double dp_ex = sigma/r_drop;
    //double dp_ex = 2.0*sigma/r_drop;

    FOR_ICV {
      
      // circle...
      const double r = sqrt( pow(x_cv[icv][0]-x_drop,2) + pow(x_cv[icv][1]-y_drop,2) );

      // sphere...
      //const double r = sqrt( pow(x_cv[icv][0]-x_drop,2) + pow(x_cv[icv][1]-y_drop,2) + pow(x_cv[icv][2]-z_drop,2) );

      const double magu = MAG(u[icv]);
      umax = max(umax,magu);
      utot += magu*vol_cv[icv];
      
      pmin = min(pmin,p[icv]);
      pmax = max(pmax,p[icv]);

      if (r <= r_drop) {
        vin_tot += vol_cv[icv];
        pin_tot += p[icv]*vol_cv[icv];
        if (r <= 0.5*r_drop) {
          vin_part += vol_cv[icv];
          pin_part += p[icv]*vol_cv[icv];
        }
      }
      else {
        vout_tot += vol_cv[icv];
        pout_tot += p[icv]*vol_cv[icv];
        if (r >= 1.5*r_drop) {
          vout_part += vol_cv[icv];
          pout_part += p[icv]*vol_cv[icv];
        }
      }
      vol += vol_cv[icv];

    }

    double my_buf_max[3] = {umax,pmax,-pmin};
    double my_buf_sum[10] = {utot,vol,pin_tot,vin_tot,pout_tot,vout_tot,pin_part,vin_part,pout_part,vout_part};
    double buf_max[3], buf_sum[10];
    MPI_Reduce(my_buf_max,buf_max,3,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    MPI_Reduce(my_buf_sum,buf_sum,10,MPI_DOUBLE,MPI_SUM,0,mpi_comm); 

//    if (mpi_rank == 0) {
//      cout << "|pin_tot-pout_tot-dp|/dp: " << fabs(buf_sum[2]/buf_sum[3]-buf_sum[4]/buf_sum[5]-dp_ex)/dp_ex << endl;
//      cout << "|pin_part-pout_part-dp|/dp: " << fabs(buf_sum[6]/buf_sum[7]-buf_sum[8]/buf_sum[9]-dp_ex)/dp_ex << endl;
//      cout << "|pmax-pmin-dp|/dp: " << fabs(buf_max[1]+buf_max[2]-dp_ex)/dp_ex << endl;
//      cout << "Linf_u: " << buf_max[0] << endl;
//      cout << "L1_u: " << buf_sum[0]/buf_sum[1] << endl; 
//    }

  }

};

class PoolVofSolver : public VofSolver {
public:

  void initialHook() {

    COUT1("PoolVofSolver::initialHook()");
    
    if (!checkDataFlag("vof")) {
      assert(!checkDataFlag("u"));

      // liquid at the bottom up to y = -1...

      const double y_pool = -0.5;
      FOR_ICV_G {
        
        FOR_I3 u[icv][i] = 0.0;
        g[icv] = y_pool-x_cv[icv][1];

      }

      buildVofFromG();
      
      /*
      updateCvData(vof);
      FOR_ICV {
        const double r = sqrt( pow(x_cv[icv][0]-x_drop,2) + pow(x_cv[icv][1]-y_drop,2) );
        if (r <= r_drop+r_vv[icv]) {
          u[icv][1] = -0.1*vof[icv];
        }
      }
      updateCvData(u);
      */
      
    }
    
  }

};

class DropPoolVofSolver : public VofSolver {
public:

  void initialHook() {

    COUT1("DropPoolVofSolver::initialHook()");
    
    if (!checkDataFlag("vof")) {
      assert(!checkDataFlag("u"));

      // drop of vof = 1 at (0,0.5), r = 0.25
      // plus a liquid at the bottom up to y = -0.5...

      const double x_drop = 0.0;
      const double y_drop = 0.5;
      const double r_drop = 0.25;
      const double y_pool = -0.5;
      FOR_ICV_G {
        
        const double r = sqrt( pow(x_cv[icv][0]-x_drop,2) + pow(x_cv[icv][1]-y_drop,2) );
        
        FOR_I3 u[icv][i] = 0.0;
        g[icv] = max(r_drop-r,y_pool-x_cv[icv][1]);

      }

      buildVofFromG();
      
      /*
      updateCvData(vof);
      FOR_ICV {
        const double r = sqrt( pow(x_cv[icv][0]-x_drop,2) + pow(x_cv[icv][1]-y_drop,2) );
        if (r <= r_drop+r_vv[icv]) {
          u[icv][1] = -0.1*vof[icv];
        }
      }
      updateCvData(u);
      */
      
    }
    
  }

};
class ConvectionVofSolver : public VofSolver {
public:

  double* vof0;
  ConvectionVofSolver() {
    vof0 = NULL; 
  }

  void initData() {
    VofSolver::initData();
    vof0 = new double[ncv];
  }

  ~ConvectionVofSolver() {
    DELETE(vof0);
  }

  void initialHook() {

    COUT1("ConvectionVofSolver::initialHook()");
    
    assert(!checkDataFlag("u"));
    assert(!checkDataFlag("vof"));
      
     
    const double radius = getDoubleParam("RADIUS",0.1);
    const double xc[2] = {0.0, 0.0};
    
    FOR_ICV_G {
      double dx[2] = { x_cv[icv][0]-xc[0], x_cv[icv][1]-xc[1] };
      g[icv] = radius - sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
    }

  
    buildVofFromG();
    calcNormal();
    calcGfromVof();
    FOR_ICV_G {
      u[icv][0] = 1.0;
      u[icv][1] = 0.0;
      u[icv][2] = 0.0;
    }
  
    FOR_ICV vof0[icv] = vof[icv];

  }

  void temporalHook() {

   if (step%check_interval == 0) { 
    double shape_err = 0.0;
    double vol_err = 0.0;
    double vol = 0.0;
    FOR_ICV {
      vol += vof0[icv]*vol_cv[icv];
      shape_err += fabs(vof[icv]-vof0[icv])*vol_cv[icv];
      vol_err += (vof[icv]-vof0[icv])*vol_cv[icv];
    }
    
    // print total L2/L1/Linf errs
    double my_buf_sum[3] = {shape_err,vol_err,vol};
    double buf_sum[3]; 
    MPI_Reduce(my_buf_sum,buf_sum,3,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    
    if (mpi_rank == 0) {
      
      cout << " > Shape error: " << buf_sum[0]/buf_sum[2] << endl;
      cout << " > Volume error: " << buf_sum[1] << endl;
      cout << " > Total volume: " << buf_sum[2] << endl;
      
    }
   } 
    
  }
  
};

class LiquidjetVofSolver : public VofSolver {
public:

  LiquidjetVofSolver() {
  }

  void initData() {
    VofSolver::initData();
  }

  ~LiquidjetVofSolver() {
  }

  void initialHook() {

    COUT1("LiquidjetVofSolver::initialHook()");
    
    assert(!checkDataFlag("u"));
    assert(!checkDataFlag("vof"));
    FOR_ICV_G {
      if (x_cv[icv][0] < -3.0e-4 ) vof[icv] = 1.0;
      else vof[icv] = 0.0;
    }
     
    calcGfromVof();
    calcNormal();
    FOR_ICV_G {
      u[icv][0] = 0.0;
      u[icv][1] = 0.0;
      u[icv][2] = 0.0;
    }
  
  }

  void temporalHook() {
  }
  
};


class DropBreakSolver : public VofSolver {
public:

  DropBreakSolver() {
  }

  void initData() {
    VofSolver::initData();
  }

  ~DropBreakSolver() {
  }

  void initialHook() {

    COUT1("DropBreakSolver::initialHook()");
     
    if (!checkDataFlag("vof")) {
      assert(!checkDataFlag("u"));
      const double radius = getDoubleParam("DROP_RADIUS",0.5);
      const double xc[3] = {0.0, 0.0, 0.0};
      
      FOR_ICV_G {
        double dx[3] = { x_cv[icv][0]-xc[0], x_cv[icv][1]-xc[1], x_cv[icv][2]-xc[2] };
        g[icv] = radius - sqrt(dx[0]*dx[0] + dx[1]*dx[1]+ dx[2]*dx[2]);
      }
      
      
      buildVofFromG();
      calcNormal();
      calcGfromVof();
      FOR_ICV_G {
        const double R = x_cv[icv][0]*x_cv[icv][0] + x_cv[icv][1]*x_cv[icv][1];
        u[icv][0] = 0.0;
        u[icv][1] = 0.0;
        u[icv][2] = 0.0;
        //u[icv][2] = 0.5*x_cv[icv][2];
      }
    } 
  }

  void calcEulerVortex(double* u, double* x, const double t) {

     // for periodic grids, the xmin, xmax needs to be set...
    const double xmin = -5.0;
    const double xmax =  5.0;
    const double ymin = -5.0*sqrt(3.0)/2.0;
    const double ymax =  5.0*sqrt(3.0)/2.0;

    const double u_inf = 0.5;

    // x,y position of the vortex
    // centered...
    const double x0 = 0.0; //-3.5;
    const double y0 = 0.0; //-3.5;

    // direction of propagation
    const double theta = M_PI/3.0;
    //const double theta = M_PI/4.0;
    //const double theta = 0.0;
    const double cos_theta = cos(theta);
    const double sin_theta = sin(theta);
    //const double Ma_inf = u_inf/sqrt(gamma*p_inf/rho_inf);
    const double rc = 1.0;

    // circulation parameter...
    //const double e_twopi = 0.001; // very weak
    //const double e_twopi = 0.005;
    const double e_twopi = 0.08; // normal
    //double e_twopi = 0.1;
    //double e_twopi = 0.4; //very strong
    //double e_twopi = 1.0; //very very strong

    // setup...
   // const double coeff = 0.5 * e_twopi*e_twopi * (gamma-1.0) * Ma_inf*Ma_inf;

    double dx = x[0] - x0 - u_inf*cos_theta*t;
    double dy = x[1] - y0 - u_inf*sin_theta*t;

    // if periodic, shift the exact solution so it is aligned with
    // the charles calculation...
  //  if (periodic_flag == 1) {
  //    while(dx > xmax) dx -= (xmax-xmin);
  //    while(dx < xmin) dx += (xmax-xmin);
  //    while(dy > ymax) dy -= (ymax-ymin);
  //    while(dy < ymin) dy += (ymax-ymin);
  //  }

    const double f0 = 1.0 - (( dx*dx ) + ( dy*dy ))/( rc*rc );
    u[0] = u_inf*( cos_theta - e_twopi * ( dy )/rc * exp( f0 / 2.0 ) );
    u[1] = u_inf*( sin_theta + e_twopi * ( dx )/rc * exp( f0 / 2.0 ) );
    u[2] = 0.0;


  }
  void temporalHook() {

   if (step%check_interval == 0) { 
   
   }
  }
  
};


class RTIVofSolver: public VofSolver{
   
  void initialHook(){

    if (step == 0) {

      COUT1("RTIVofSolver:initialHook()");

      const double init_center[2] = { 0.0, 0.0};
      const double init_amplitude = 0.05;
      const double init_wavelength =  1.0;
      
      FOR_ICV_G {
	g[icv] = x_cv[icv][1] - init_center[1] - init_amplitude*cos(2.0*M_PI*(x_cv[icv][0]-init_center[0])/init_wavelength);
	FOR_I3 u[icv][i] = 0.0; 
      }
      
      buildVofFromG();
      calcNormal();
      calcGfromVof();
    }
      
  }

};

class HsquareSolver : public VofSolver {
public:

  double* vof_ini;
  double* vof_ex;
  HsquareSolver() {
    vof_ini = NULL; 
    vof_ex = NULL; registerCvData(vof_ex,"vof_ex",CAN_WRITE_DATA);
  }

  void initData() {
    VofSolver::initData();
    vof_ini = new double[ncv_g];
    vof_ex = new double[ncv_g];
  }

  ~HsquareSolver() {
    DELETE(vof_ex);
  }

  void initialHook() {

    COUT1("HSquareSolver::initialHook()");
   
   
    if (!checkDataFlag("u") && !checkDataFlag("vof") ) { 
      
      FOR_ICV_G {
        u[icv][0] = 2.0;
        u[icv][1] = 1.0;
        u[icv][2] = 0.0;
      } 


      const double x_disk[2] = { 0.8, 0.8};
      const double width1 = 0.2;
      const double width2 =0.4;

      FOR_ICV_G {
        double dx[2] = { x_cv[icv][0]-x_disk[0], x_cv[icv][1]-x_disk[1] };
        double G_hollow_square = -1.0E20;
        double do1=1.0E20, do2=1.0E20, do3=1.0E20, do4=1.0E20, do5=1.0E20, do6=1.0E20, do7=1.0E20, do8 = 1.0E20;
        double di1=1.0E20, di2=1.0E20, di3=1.0E20, di4=1.0E20, di5=1.0E20, di6=1.0E20, di7=1.0E20, di8 = 1.0E20;
        
        if (fabs(dx[1]) <= width2) {
          do1 = fabs(dx[0] - width2); //right side
          do2 = fabs(dx[0] + width2); //left side
        }
        if (fabs(dx[1]) <= width1) {
          di1 = fabs(dx[0] - width1); //right side
          di2 = fabs(dx[0] + width1); //left side
        }
        if (fabs(dx[0]) <= width2) {
          do3 = fabs(dx[1] - width2); // top
          do4 = fabs(dx[1] + width2); // bottom
        }
        if (fabs(dx[0]) <= width1) {
          di3 = fabs(dx[1] - width1); // top
          di4 = fabs(dx[1] + width1); // bottom
        }
        
        do5 = sqrt( (dx[0]-width2)*(dx[0]-width2) + (dx[1]-width2)*(dx[1]-width2)); //upper right corner
        do6 = sqrt( (dx[0]+width2)*(dx[0]+width2) + (dx[1]-width2)*(dx[1]-width2)); //upper left corner
        do7 = sqrt( (dx[0]+width2)*(dx[0]+width2) + (dx[1]+width2)*(dx[1]+width2)); //lower left corner
        do8 = sqrt( (dx[0]-width2)*(dx[0]-width2) + (dx[1]+width2)*(dx[1]+width2)); //lower right corner
        di5 = sqrt( (dx[0]-width1)*(dx[0]-width1) + (dx[1]-width1)*(dx[1]-width1)); //upper right corner
        di6 = sqrt( (dx[0]+width1)*(dx[0]+width1) + (dx[1]-width1)*(dx[1]-width1)); //upper left corner
        di7 = sqrt( (dx[0]+width1)*(dx[0]+width1) + (dx[1]+width1)*(dx[1]+width1)); //lower left corner
        di8 = sqrt( (dx[0]-width1)*(dx[0]-width1) + (dx[1]+width1)*(dx[1]+width1)); //lower right corner
        
        double d_min = 1.0E20;
        
        d_min = min(d_min,do1);
        d_min = min(d_min,do2);
        d_min = min(d_min,do3);
        d_min = min(d_min,do4);
        d_min = min(d_min,do5);
        d_min = min(d_min,do6);
        d_min = min(d_min,do7);
        d_min = min(d_min,do8);
        
        
        d_min = min(d_min,di1);
        d_min = min(d_min,di2);
        d_min = min(d_min,di3);
        d_min = min(d_min,di4);
        d_min = min(d_min,di5);
        d_min = min(d_min,di6);
        d_min = min(d_min,di7);
        d_min = min(d_min,di8);
        
        // inside square ...
        if (fabs(dx[0]) >= width1 && fabs(dx[0]) <= width2 && fabs(dx[1]) <= width2) g[icv] = d_min;
        else if (fabs(dx[1]) >= width1 && fabs(dx[1]) <= width2 && fabs(dx[0]) <= width2) g[icv] = d_min;
        else  g[icv]  = -d_min;
        
      } 
      /*  
      FOR_ICV_G {
        vof[icv] = 0.0;
        double x[2] = {x_cv[icv][0] - x_disk[0], x_cv[icv][1] - x_disk[1]};
        if (fabs(x[0]) >= width1 && fabs(x[0]) <= width2 && fabs(x[1]) <= width2) vof[icv] = 1.0;
        if (fabs(x[1]) >= width1 && fabs(x[1]) <= width2 && fabs(x[0]) <= width2) vof[icv] = 1.0;
      }
     */ 

      //updateInterface();
      buildVofFromG();
      calcNormal();
      calcGfromVof();
      FOR_ICV_G vof_ini[icv] = vof[icv];
    }


  }

  void temporalHook() {

     FOR_ICV_G {
      u[icv][0] = 2.0;
      u[icv][1] = 1.0;
      u[icv][2] = 0.0;
     }

     // exact solution ...
    
     const double x_disk[2] = { 0.8 + 2.0*time, 0.8+1.0*time};
     const double width1 = 0.2;
     const double width2 =0.4;
     double g_ex[ncv_g];
     
     FOR_ICV_G {
       double dx[2] = { x_cv[icv][0]-x_disk[0], x_cv[icv][1]-x_disk[1] };
       double G_hollow_square = -1.0E20;
       double do1=1.0E20, do2=1.0E20, do3=1.0E20, do4=1.0E20, do5=1.0E20, do6=1.0E20, do7=1.0E20, do8 = 1.0E20;
       double di1=1.0E20, di2=1.0E20, di3=1.0E20, di4=1.0E20, di5=1.0E20, di6=1.0E20, di7=1.0E20, di8 = 1.0E20;
       
       if (fabs(dx[1]) <= width2) {
         do1 = fabs(dx[0] - width2); //right side
         do2 = fabs(dx[0] + width2); //left side
       }
       if (fabs(dx[1]) <= width1) {
         di1 = fabs(dx[0] - width1); //right side
         di2 = fabs(dx[0] + width1); //left side
       }
       if (fabs(dx[0]) <= width2) {
         do3 = fabs(dx[1] - width2); // top
         do4 = fabs(dx[1] + width2); // bottom
       }
       if (fabs(dx[0]) <= width1) {
         di3 = fabs(dx[1] - width1); // top
         di4 = fabs(dx[1] + width1); // bottom
       }
       
       do5 = sqrt( (dx[0]-width2)*(dx[0]-width2) + (dx[1]-width2)*(dx[1]-width2)); //upper right corner
       do6 = sqrt( (dx[0]+width2)*(dx[0]+width2) + (dx[1]-width2)*(dx[1]-width2)); //upper left corner
       do7 = sqrt( (dx[0]+width2)*(dx[0]+width2) + (dx[1]+width2)*(dx[1]+width2)); //lower left corner
       do8 = sqrt( (dx[0]-width2)*(dx[0]-width2) + (dx[1]+width2)*(dx[1]+width2)); //lower right corner
       di5 = sqrt( (dx[0]-width1)*(dx[0]-width1) + (dx[1]-width1)*(dx[1]-width1)); //upper right corner
       di6 = sqrt( (dx[0]+width1)*(dx[0]+width1) + (dx[1]-width1)*(dx[1]-width1)); //upper left corner
       di7 = sqrt( (dx[0]+width1)*(dx[0]+width1) + (dx[1]+width1)*(dx[1]+width1)); //lower left corner
       di8 = sqrt( (dx[0]-width1)*(dx[0]-width1) + (dx[1]+width1)*(dx[1]+width1)); //lower right corner
       
       double d_min = 1.0E20;
       
       d_min = min(d_min,do1);
       d_min = min(d_min,do2);
       d_min = min(d_min,do3);
       d_min = min(d_min,do4);
       d_min = min(d_min,do5);
       d_min = min(d_min,do6);
       d_min = min(d_min,do7);
       d_min = min(d_min,do8);
       
       
       d_min = min(d_min,di1);
       d_min = min(d_min,di2);
       d_min = min(d_min,di3);
       d_min = min(d_min,di4);
       d_min = min(d_min,di5);
       d_min = min(d_min,di6);
       d_min = min(d_min,di7);
       d_min = min(d_min,di8);
       
       // inside square ...
       if (fabs(dx[0]) >= width1 && fabs(dx[0]) <= width2 && fabs(dx[1]) <= width2) g_ex[icv] = d_min;
       else if (fabs(dx[1]) >= width1 && fabs(dx[1]) <= width2 && fabs(dx[0]) <= width2) g_ex[icv] = d_min;
       else  g_ex[icv]  = -d_min;
       
     } 
     
     FOR_ICV_G {
        vof_ex[icv] = 0.5*(1.0+tanh(beta[icv]*(g_ex[icv])));
        if (vof_ex[icv] < vof_zero) vof_ex[icv] = 0.0;
        if (vof_ex[icv] > 1.0-vof_zero) vof_ex[icv] = 1.0;
      }
    
      double shape_err = 0.0;
      double denom = 0.0;
      FOR_ICV {
        shape_err += fabs(vof[icv]-vof_ex[icv]);
        denom += vof_ini[icv];

      }
      
      double my_buf_sum[2] = {shape_err,denom};
      double buf_sum[2]; 
      MPI_Reduce(my_buf_sum,buf_sum,2,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      
      if (mpi_rank == 0) {
        
        cout << " > Relative Shape error: " << buf_sum[0]/buf_sum[1] << endl;
        
      }
      
      
  }

};


class SquareSolver : public VofSolver {
public:

  double* vof0;
  SquareSolver() {
    vof0 = NULL; 
  }

  void initData() {
    VofSolver::initData();
    vof0 = new double[ncv];
  }

  ~SquareSolver() {
    DELETE(vof0);
  }

  void initialHook() {

    COUT1("SquareSolver::initialHook()");
    
    assert(!checkDataFlag("u"));
    assert(!checkDataFlag("vof"));
      
    FOR_ICV_G {
      u[icv][0] = 0.5-x_cv[icv][1];
      u[icv][1] = x_cv[icv][0]-0.5;
      u[icv][2] = 0.0;
    }
      
    const double x_disk[2] = { 0.5, 0.75};
    const double init_width = 0.05;
    const double init_height = 0.25;
    const double init_radius = 0.15;
    double G_notched_circle = 0.0;
    
    
    FOR_ICV_G {
      double dx[2] = { x_cv[icv][0]-x_disk[0], x_cv[icv][1]-x_disk[1] };
      double c = init_radius - sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
      double b1 = x_disk[0]-0.5*init_width;
      double b2 = x_disk[0]+0.5*init_width;
      double h1 = x_disk[1]-init_radius*cos(asin(0.5*init_width/init_radius));
      double h2 = x_disk[1]-init_radius+init_height;
      double xyz[2] = {x_cv[icv][0], x_cv[icv][1]};
      double b = 0.0;
      
      if  (c >= 0.0 && xyz[0] <= b1 && xyz[1] <= h2 ) {
        b = b1-xyz[0];
        G_notched_circle = min(c,b);
      }

      else if (c >= 0.0 && xyz[0] >= b2 && xyz[1] <= h2) {
        b = xyz[0]-b2;
        G_notched_circle = min(c,b);
      }
      
      else if (c >= 0.0 && xyz[0] >= b1 && xyz[0] <= b2 && xyz[1] >= h2) {
        b = xyz[1]-h2;
        G_notched_circle = min(c,b);
      }
      else if (c >= 0.0 && xyz[0] <= b1 && xyz[1] >= h2) {
        b = sqrt( (xyz[0]-b1)*(xyz[0]-b1)+(xyz[1]-h2)*(xyz[1]-h2));
        G_notched_circle = min(c,b);
      }
      else if (c >= 0.0 && xyz[0] >= b2 && xyz[1] >= h2) {
        b = sqrt( (xyz[0]-b2)*(xyz[0]-b2)+(xyz[1]-h2)*(xyz[1]-h2));
        G_notched_circle = min(c,b);
      }
      else if (xyz[0] >= b1 && xyz[0] <= b2 && xyz[1] <= h2 && xyz[1] >= h1) {
        G_notched_circle = min(fabs(xyz[0]-b1),fabs(xyz[0]-b2));
        G_notched_circle = -min(G_notched_circle,fabs(xyz[1]-h2));
      }
      else if (xyz[0] >= b1 && xyz[0] <= b2 && xyz[1] <= h1) {
        double a1 =sqrt((xyz[0]-b1)*(xyz[0]-b1)+(xyz[1]-h1)*(xyz[1]-h1));
        double a2 =sqrt((xyz[0]-b2)*(xyz[0]-b2)+(xyz[1]-h1)*(xyz[1]-h1));
        G_notched_circle = -min(a1,a2);
      }
      else
        G_notched_circle = c;
      
      g[icv] = G_notched_circle;
    }

    buildVofFromG();
    calcNormal();
    calcGfromVof();

    FOR_ICV vof0[icv] = vof[icv];

  }

  void temporalHook() {

     FOR_ICV_G {
      u[icv][0] = 0.5-x_cv[icv][1];
      u[icv][1] = x_cv[icv][0]-0.5;
      u[icv][2] = 0.0;
    }
    
    // solvePAndCorrectU();
    // FOR_IFA q0_fa[ifa] = q_fa[ifa];
    
    
      double shape_err = 0.0;
      double vol_err = 0.0;
      double vol = 0.0;
      FOR_ICV {
        vol += vof0[icv]*vol_cv[icv];
        shape_err += fabs(vof[icv]-vof0[icv])*vol_cv[icv];
        vol_err += (vof[icv]-vof0[icv])*vol_cv[icv];
      }
      
      // print total L2/L1/Linf errs
      double my_buf_sum[3] = {shape_err,vol_err,vol};
      double buf_sum[3]; 
      MPI_Reduce(my_buf_sum,buf_sum,3,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      
      if (mpi_rank == 0) {
        
        cout << " > Shape error: " << buf_sum[0]/buf_sum[2] << endl;
        cout << " > Volume error: " << buf_sum[1] << endl;
        cout << " > Total volume: " << buf_sum[2] << endl;
        
      }
      
      
  }

};


class ZalesaksSolver : public VofSolver {
public:

  double* vof0;
  ZalesaksSolver() {
    vof0 = NULL; 
  }

  void initData() {
    VofSolver::initData();
    vof0 = new double[ncv];
  }

  ~ZalesaksSolver() {
    DELETE(vof0);
  }

  void initialHook() {

    COUT1("ZalesaksSolver::initialHook()");
    
    assert(!checkDataFlag("u"));
    assert(!checkDataFlag("vof"));
      
    FOR_ICV_G {
      u[icv][0] = 0.5-x_cv[icv][1];
      u[icv][1] = x_cv[icv][0]-0.5;
      u[icv][2] = 0.0;
    }
      
    const double x_disk[2] = { 0.5, 0.75};
    const double init_width = 0.05;
    const double init_height = 0.25;
    const double init_radius = 0.15;
    double G_notched_circle = 0.0;
    
    
    FOR_ICV_G {
      double dx[2] = { x_cv[icv][0]-x_disk[0], x_cv[icv][1]-x_disk[1] };
      double c = init_radius - sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
      double b1 = x_disk[0]-0.5*init_width;
      double b2 = x_disk[0]+0.5*init_width;
      double h1 = x_disk[1]-init_radius*cos(asin(0.5*init_width/init_radius));
      double h2 = x_disk[1]-init_radius+init_height;
      double xyz[2] = {x_cv[icv][0], x_cv[icv][1]};
      double b = 0.0;
      
      if  (c >= 0.0 && xyz[0] <= b1 && xyz[1] <= h2 ) {
        b = b1-xyz[0];
        G_notched_circle = min(c,b);
      }

      else if (c >= 0.0 && xyz[0] >= b2 && xyz[1] <= h2) {
        b = xyz[0]-b2;
        G_notched_circle = min(c,b);
      }
      
      else if (c >= 0.0 && xyz[0] >= b1 && xyz[0] <= b2 && xyz[1] >= h2) {
        b = xyz[1]-h2;
        G_notched_circle = min(c,b);
      }
      else if (c >= 0.0 && xyz[0] <= b1 && xyz[1] >= h2) {
        b = sqrt( (xyz[0]-b1)*(xyz[0]-b1)+(xyz[1]-h2)*(xyz[1]-h2));
        G_notched_circle = min(c,b);
      }
      else if (c >= 0.0 && xyz[0] >= b2 && xyz[1] >= h2) {
        b = sqrt( (xyz[0]-b2)*(xyz[0]-b2)+(xyz[1]-h2)*(xyz[1]-h2));
        G_notched_circle = min(c,b);
      }
      else if (xyz[0] >= b1 && xyz[0] <= b2 && xyz[1] <= h2 && xyz[1] >= h1) {
        G_notched_circle = min(fabs(xyz[0]-b1),fabs(xyz[0]-b2));
        G_notched_circle = -min(G_notched_circle,fabs(xyz[1]-h2));
      }
      else if (xyz[0] >= b1 && xyz[0] <= b2 && xyz[1] <= h1) {
        double a1 =sqrt((xyz[0]-b1)*(xyz[0]-b1)+(xyz[1]-h1)*(xyz[1]-h1));
        double a2 =sqrt((xyz[0]-b2)*(xyz[0]-b2)+(xyz[1]-h1)*(xyz[1]-h1));
        G_notched_circle = -min(a1,a2);
      }
      else
        G_notched_circle = c;
      
      g[icv] = G_notched_circle;
    }

    buildVofFromG();
    calcNormal();
    calcGfromVof();

    FOR_ICV vof0[icv] = vof[icv];

  }

  void temporalHook() {

     FOR_ICV_G {
      u[icv][0] = 0.5-x_cv[icv][1];
      u[icv][1] = x_cv[icv][0]-0.5;
      u[icv][2] = 0.0;
    }
    
    // solvePAndCorrectU();
    // FOR_IFA q0_fa[ifa] = q_fa[ifa];
    
    
      double shape_err = 0.0;
      double vol_err = 0.0;
      double vol = 0.0;
      FOR_ICV {
        vol += vof0[icv]*vol_cv[icv];
        shape_err += fabs(vof[icv]-vof0[icv])*vol_cv[icv];
        vol_err += (vof[icv]-vof0[icv])*vol_cv[icv];
      }
      
      // print total L2/L1/Linf errs
      double my_buf_sum[3] = {shape_err,vol_err,vol};
      double buf_sum[3]; 
      MPI_Reduce(my_buf_sum,buf_sum,3,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      
      if (mpi_rank == 0) {
        
        cout << " > Shape error: " << buf_sum[0]/buf_sum[2] << endl;
        cout << " > Volume error: " << buf_sum[1] << endl;
        cout << " > Total volume: " << buf_sum[2] << endl;
        
      }
      
      
  }

};

class MZalesaksSolver : public VofSolver {
public:

  double* vof0;
  MZalesaksSolver() {
    vof0 = NULL; 
  }

  void initData() {
    VofSolver::initData();
    vof0 = new double[ncv];
  }

  ~MZalesaksSolver() {
    DELETE(vof0);
  }

  void initialHook() {

    COUT1("ZalesaksSolver::initialHook()");
    
    assert(!checkDataFlag("u"));
    assert(!checkDataFlag("vof"));
      
    FOR_ICV_G {
      u[icv][0] = 0.5*(2.0-x_cv[icv][1]);
      u[icv][1] = 0.5*(x_cv[icv][0]-2.0);
      u[icv][2] = 0.0;
    }
      
    const double x_disk[2] = { 2.0, 2.75};
    const double init_width = 0.12;
    const double init_height = 0.55;
    const double init_radius = 0.5;
    double G_notched_circle = 0.0;

    
    FOR_ICV_G {
      double dx[2] = { x_cv[icv][0]-x_disk[0], x_cv[icv][1]-x_disk[1] };
      double c = init_radius - sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
      double b1 = x_disk[0]-0.5*init_width;
      double b2 = x_disk[0]+0.5*init_width;
      double h1 = x_disk[1]-init_radius*cos(asin(0.5*init_width/init_radius));
      double h2 = x_disk[1]-init_radius+init_height;
      double xyz[2] = {x_cv[icv][0], x_cv[icv][1]};
      double b = 0.0;
      
      if  (c >= 0.0 && xyz[0] <= b1 && xyz[1] <= h2 ) {
        b = b1-xyz[0];
        G_notched_circle = min(c,b);
      }

      else if (c >= 0.0 && xyz[0] >= b2 && xyz[1] <= h2) {
        b = xyz[0]-b2;
        G_notched_circle = min(c,b);
      }
      
      else if (c >= 0.0 && xyz[0] >= b1 && xyz[0] <= b2 && xyz[1] >= h2) {
        b = xyz[1]-h2;
        G_notched_circle = min(c,b);
      }
      else if (c >= 0.0 && xyz[0] <= b1 && xyz[1] >= h2) {
        b = sqrt( (xyz[0]-b1)*(xyz[0]-b1)+(xyz[1]-h2)*(xyz[1]-h2));
        G_notched_circle = min(c,b);
      }
      else if (c >= 0.0 && xyz[0] >= b2 && xyz[1] >= h2) {
        b = sqrt( (xyz[0]-b2)*(xyz[0]-b2)+(xyz[1]-h2)*(xyz[1]-h2));
        G_notched_circle = min(c,b);
      }
      else if (xyz[0] >= b1 && xyz[0] <= b2 && xyz[1] <= h2 && xyz[1] >= h1) {
        G_notched_circle = min(fabs(xyz[0]-b1),fabs(xyz[0]-b2));
        G_notched_circle = -min(G_notched_circle,fabs(xyz[1]-h2));
      }
      else if (xyz[0] >= b1 && xyz[0] <= b2 && xyz[1] <= h1) {
        double a1 =sqrt((xyz[0]-b1)*(xyz[0]-b1)+(xyz[1]-h1)*(xyz[1]-h1));
        double a2 =sqrt((xyz[0]-b2)*(xyz[0]-b2)+(xyz[1]-h1)*(xyz[1]-h1));
        G_notched_circle = -min(a1,a2);
      }
      else
        G_notched_circle = c;
      
      g[icv] = G_notched_circle;
    }

    buildVofFromG();
    calcNormal();
    calcGfromVof();
    FOR_ICV vof0[icv] = vof[icv];

  }

  void temporalHook() {

     FOR_ICV_G {
      u[icv][0] = 0.5*(2.0-x_cv[icv][1]);
      u[icv][1] = 0.5*(x_cv[icv][0]-2.0);
      u[icv][2] = 0.0;
    }
    
    //solvePAndCorrectU();
    //FOR_IFA q0_fa[ifa] = q_fa[ifa];
    
    
    double shape_err = 0.0;
    double shape_err2 = 0.0;
    double vol_err = 0.0;
    double vol = 0.0;
    double vof_sum = 0.0;
    double vof_sum2 = 0.0;
    FOR_ICV {
      vol += vol_cv[icv];
      shape_err += fabs(vof[icv]-vof0[icv])*vol_cv[icv];
      shape_err2 += fabs(vof[icv]-vof0[icv]);
      vof_sum   += vof0[icv];
      vof_sum2   += vof0[icv]*vol_cv[icv];
      vol_err += (vof[icv]-vof0[icv])*vol_cv[icv];
    }
      
    // print total L2/L1/Linf errs
    double my_buf_sum[6] = {shape_err,vol_err,vol,shape_err2, vof_sum,vof_sum2};
    double buf_sum[6]; 
    MPI_Reduce(my_buf_sum,buf_sum,6,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    
    if (mpi_rank == 0) {
      
      cout << "> Shape error1: " << buf_sum[0]/buf_sum[5] << endl;
      cout << "> Volume error: " << buf_sum[1] << endl;
      cout << "> Total volume: " << buf_sum[2] << endl;
      cout << "> Shape error2: " << buf_sum[3]/buf_sum[4] << endl;
      
    }
    
    FOR_ICV {
      if (fabs(x_cv[icv][1]) < 6.26e-06 && fabs(x_cv[icv][2])< 6.26e-06 ) {
        cout << x_cv[icv][0] << " " << x_cv[icv][1] << " " << x_cv[icv][2] <<  "  " << vof[icv] <<  endl;
      }
    }
  }
    
  
};



class TestVofSolver : public VofSolver {
public:

  double* vof0;
  TestVofSolver() {
    vof0 = NULL; 
  }

  void initData() {
    VofSolver::initData();
    vof0 = new double[ncv];
  }

  ~TestVofSolver() {
    DELETE(vof0);
  }

  void initialHook() {

    COUT1("TestVofSolver::initialHook()");
    
    if (!checkDataFlag("vof")) {
      assert(!checkDataFlag("u"));

      // translation...
      //const double x_drop = 0.0;
      //const double y_drop = 0.0;
      //const double r_drop = 0.5;
      
      // rotation...
      //const double x_drop = 0.0;
      //const double y_drop = 0.0;
      //const double r_drop = 0.25;

      // deformation...
      const double x_drop = 0.0;
      const double y_drop = 0.0;
      const double r_drop = getDoubleParam("RADIUS",0.1);
      const double ax = 1.0;
      const double ay = 1.0;
      //const double ax = 3.0;
      //const double ay = 0.5;

      FOR_ICV_G {
        const double r = sqrt( pow(x_cv[icv][0]-x_drop,2)/ax/ax + pow(x_cv[icv][1]-y_drop,2)/ay/ay );
        g[icv] = r_drop-r;
      }
      buildVofFromG();
      calcNormal();
      calcGfromVof();

      FOR_ICV_G {
        // translation...
        //u[icv][0] = 1.0/sqrt(2.0);
        //u[icv][1] = 1.0/sqrt(2.0);
        //u[icv][2] = 0.0; 

        // rotation...
        //u[icv][0] = -2.0*M_PI*x_cv[icv][1];
        //u[icv][1] = +2.0*M_PI*x_cv[icv][0];
        //u[icv][2] = 0.0; 

        // deformation...
        u[icv][0] = 2.0;
        u[icv][1] = 1.0;
        u[icv][2] = 0.0; 
      }

      FOR_ICV vof0[icv] = vof[icv];

    }
    
  }

  void temporalHook() {

    
    //  FOR_ICV_G {
    //    u[icv][0] = 1.0;
    //    u[icv][1] = 0.5;
    //    u[icv][2] = 0.0; 
    //  }
    
   // solvePAndCorrectU();
   // FOR_IFA q0_fa[ifa] = q_fa[ifa];
    
    const double T = 4.0;


    double shape_err = 0.0;
      double vol_err = 0.0;
      double vol = 0.0;
      FOR_ICV {
        vol += vof0[icv]*vol_cv[icv];;
        shape_err += fabs(vof[icv]-vof0[icv])*vol_cv[icv];
        vol_err += (vof[icv]-vof0[icv])*vol_cv[icv];
      }

      // print total L2/L1/Linf errs
      double my_buf_sum[3] = {shape_err,vol_err,vol};
      double buf_sum[3];
      MPI_Reduce(my_buf_sum,buf_sum,3,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

      if (mpi_rank == 0) {

        cout << "> Shape error: " << buf_sum[0]/buf_sum[2] << endl;

      }


  } 
  
};

class KinematicVofSolver : public VofSolver {
public:

  double* vof0;
  KinematicVofSolver() {
    vof0 = NULL; 
  }

  void initData() {
    VofSolver::initData();
    vof0 = new double[ncv];
  }

  ~KinematicVofSolver() {
    DELETE(vof0);
  }

  void initialHook() {

    COUT1("KinematicVofSolver::initialHook()");
    
    if (!checkDataFlag("vof")) {
      assert(!checkDataFlag("u"));

      // translation...
      //const double x_drop = 0.0;
      //const double y_drop = 0.0;
      //const double r_drop = 0.5;
      
      // rotation...
      //const double x_drop = 0.0;
      //const double y_drop = 0.0;
      //const double r_drop = 0.25;

      // deformation...
      const double x_drop = 0.5;
      const double y_drop = 0.75;
      const double r_drop = 0.15;

      FOR_ICV_G {
        const double r = sqrt( pow(x_cv[icv][0]-x_drop,2) + pow(x_cv[icv][1]-y_drop,2) );
        g[icv] = r_drop-r;
      }
      buildVofFromG();
      calcNormal();
      calcGfromVof();
      buildIband();
      
      FOR_ICV_G {
        // translation...
        //u[icv][0] = 1.0/sqrt(2.0);
        //u[icv][1] = 1.0/sqrt(2.0);
        //u[icv][2] = 0.0; 

        // rotation...
        //u[icv][0] = -2.0*M_PI*x_cv[icv][1];
        //u[icv][1] = +2.0*M_PI*x_cv[icv][0];
        //u[icv][2] = 0.0; 

        // deformation...
        u[icv][0] = -2.0*pow(sin(M_PI*x_cv[icv][0]),2)*sin(M_PI*x_cv[icv][1])*cos(M_PI*x_cv[icv][1]);
        u[icv][1] = 2.0*pow(sin(M_PI*x_cv[icv][1]),2)*sin(M_PI*x_cv[icv][0])*cos(M_PI*x_cv[icv][0]);
        u[icv][2] = 0.0; 
      }

      FOR_ICV vof0[icv] = vof[icv];

    }
    
  }

  void temporalHook() {

    
    // translation...
    /*
    const double T = 1.0;
    if (time <= 0.5*T) {
      FOR_ICV_G {
        u[icv][0] = 1.0/sqrt(2.0);
        u[icv][1] = 1.0/sqrt(2.0);
        u[icv][2] = 0.0;
      }
    }
    else {
      FOR_ICV_G {
        u[icv][0] = -1.0/sqrt(2.0);
        u[icv][1] = -1.0/sqrt(2.0);
        u[icv][2] = 0.0;
      }
    }
    */

    /*
    // rotation...
    const double T = 1.0;
    FOR_ICV_G {
      u[icv][0] = -2.0*M_PI*x_cv[icv][1];
      u[icv][1] = +2.0*M_PI*x_cv[icv][0];
      u[icv][2] = 0.0;
    }
    */

    const double T = 2.0;
    FOR_ICV_G {
      u[icv][0] = -2.0*pow(sin(M_PI*x_cv[icv][0]),2)*sin(M_PI*x_cv[icv][1])*cos(M_PI*x_cv[icv][1])*cos(M_PI*(time+0.5*dt)/T);
      u[icv][1] =  2.0*pow(sin(M_PI*x_cv[icv][1]),2)*sin(M_PI*x_cv[icv][0])*cos(M_PI*x_cv[icv][0])*cos(M_PI*(time+0.5*dt)/T);
      u[icv][2] = 0.0; 
    }

    solvePAndCorrectU();
    FOR_IFA q0_fa[ifa] = q_fa[ifa];
    
      // calculate shape error...

      double shape_err = 0.0;
      double vol_err = 0.0;
      double vol = 0.0;
      FOR_ICV {
        vol += vol_cv[icv];
        shape_err += fabs(vof[icv]-vof0[icv])*pow(vol_cv[icv],2.0/3.0);
        vol_err += (vof[icv]-vof0[icv])*vol_cv[icv];
      }

      // print total L2/L1/Linf errs
      double my_buf_sum[3] = {shape_err,vol_err,vol};
      double buf_sum[3]; 
      MPI_Reduce(my_buf_sum,buf_sum,3,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

      if (mpi_rank == 0) {

        cout << "> Shape error: " << buf_sum[0]*0.25 << endl;
        cout << "> Volume error: " << buf_sum[1]/buf_sum[2] << endl;
        cout << "> Total volume: " << buf_sum[2] << endl;

      }

      flushImages();

    //  throw(0);
 
  }

};

class Sphere3DVofSolver : public VofSolver {
public:

  double* vof0;
  Sphere3DVofSolver() {
    vof0 = NULL; 
  }

  void initData() {
    VofSolver::initData();
    vof0 = new double[ncv];
  }

  ~Sphere3DVofSolver() {
    DELETE(vof0);
  }

  void initialHook() {

    COUT1("Sphere3DVofSolver::initialHook()");
    
    assert(!checkDataFlag("u"));
    assert(!checkDataFlag("vof"));
      
     
    const double radius = 0.15;
    const double xc[3] = {0.35, 0.35,0.35};
    
    FOR_ICV_G {
      double dx[3] = { x_cv[icv][0]-xc[0], x_cv[icv][1]-xc[1], x_cv[icv][2] - xc[2] };
      g[icv] = radius - sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    }

  
    buildVofFromG();
    calcNormal();
    calcGfromVof();
    FOR_ICV_G {
      const double x = x_cv[icv][0];
      const double y = x_cv[icv][1];
      const double z = x_cv[icv][2];
      const double pi = M_PI;
      const double T = 3.0;
      u[icv][0] = 2.0*sin(pi*x)*sin(pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)*cos(pi*time*T);
      u[icv][1] =   - sin(2.0*pi*x)*sin(pi*y)*sin(pi*y)*sin(2.0*pi*z)*cos(pi*time/T);
      u[icv][2] =    -sin(2.0*pi*x)*sin(2.0*pi*y)*sin(pi*z)*sin(pi*z)*cos(pi*time/T);
    }
  
    FOR_ICV vof0[icv] = vof[icv];

  }

  void temporalHook() {

   if (step%check_interval == 0) { 
    double shape_err = 0.0;
    double vol_err = 0.0;
    double vol = 0.0;
    FOR_ICV {
      vol += vof0[icv]*vol_cv[icv];
      shape_err += fabs(vof[icv]-vof0[icv])*vol_cv[icv];
      vol_err += (vof[icv]-vof0[icv])*vol_cv[icv];
    }
    
    // print total L2/L1/Linf errs
    double my_buf_sum[3] = {shape_err,vol_err,vol};
    double buf_sum[3]; 
    MPI_Reduce(my_buf_sum,buf_sum,3,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    
    if (mpi_rank == 0) {
      
      cout << " > L1 Shape error: " << buf_sum[0] << endl;
      cout << " > Volume error: " << buf_sum[1] << endl;
      cout << " > Total volume: " << buf_sum[2] << endl;
      
    }
   } 
    
  }
  
};


class StaticVofSolver : public VofSolver {
public:

  double* kappa_err;
  double* kappa_ex;
  double* n_err;

  StaticVofSolver() {
    kappa_err = NULL; registerCvData(kappa_err,"kappa_err",0);
    kappa_ex = NULL; registerCvData(kappa_ex,"kappa_ex",0);
    n_err = NULL; registerCvData(n_err,"n_err",0);
  }

  void initData() {
    VofSolver::initData();
    kappa_err = new double[ncv];
    kappa_ex = new double[ncv];
    n_err = new double[ncv];
  }

  ~StaticVofSolver() {
    DELETE(kappa_err);
    DELETE(kappa_ex);
    DELETE(n_err);
  }

  void initialHook() {

    COUT1("StaticVofSolver::initialHook()");
    
    if (!checkDataFlag("vof")) {
      assert(!checkDataFlag("u"));
      const double r_drop = getDoubleParam("DROP_RADIUS",0.325);
      const double x_drop = 0.525;
      const double y_drop = 0.464;
      const double z_drop = 0.516;

      FOR_ICV_G {
        const double r = sqrt( pow(x_cv[icv][0]-x_drop,2) + pow(x_cv[icv][1]-y_drop,2) + pow(x_cv[icv][2]-z_drop,2) );
        g[icv] = r_drop-r;
      }
      buildVofFromG();
      calcNormal();
      calcGfromVof();

      FOR_ICV_G {
        u[icv][0] = 0.0;
        u[icv][1] = 0.0;
        u[icv][2] = 0.0;
      }
    
    }
    
  }

  void temporalHook() {

  //  const double eps = 1.0E-4;
    //check icv
  //  FOR_ICV {
  //    if (x_cv[icv][0] > -1.12499 -eps && x_cv[icv][0] < -1.12499 +eps) {
  //      if (x_cv[icv][1] > -1.67792 -eps && x_cv[icv][1] < -1.67792 +eps) {
  //        cout << "x_cv = " << icv << " " << COUT_VEC(x_cv[icv]) << " " << vof[icv] << " " << kappa[icv] << endl;
  //      }
  //    }
  //  }

    const double r_drop = getDoubleParam("DROP_RADIUS",0.325);
    const double x_drop = 0.525;
    const double y_drop = 0.464;
    const double z_drop = 0.516;
    
    {double L2_kappa = 0.0;
    double L1_kappa = 0.0;
    double Linf_kappa = 0.0;
    double L2_n = 0.0;
    double L1_n = 0.0;
    double Linf_n = 0.0;
    double cnt_int = 0.0;
    double vol_1 = 0.0;

    FOR_ICV {

      vol_1 += vof[icv]*vol_cv[icv];

      if (cv_flag[icv] >= 1) {

      
        // 45 deg line...

        //const double n_true[3] = {-1.0/sqrt(2.0),1.0/sqrt(2.0),0.0};
        //const double kappa_true = 0.0;

        // sin wave...

        //const double h = y0_wave + a0_wave*sin(2.0*M_PI*(x_cv[icv][0]-x0_wave)/l_wave) - x_cv[icv][1];
        //const double hp = 2.0*a0_wave/l_wave*M_PI*cos(2.0*M_PI*(x_cv[icv][0]-x0_wave)/l_wave);
        //const double hpp = -4.0*a0_wave*M_PI*M_PI/l_wave/l_wave*sin(2.0*M_PI*(x_cv[icv][0]-x0_wave)/l_wave);
        //const double n_true[3] = {-hp/sqrt(1.0+hp*hp),1.0/sqrt(1.0+hp*hp),0.0};
        //const double den = 1.0+hp*hp;
        //const double kappa_true = -hpp/sqrt(den*den*den);

        // circle...
        //const double r = sqrt( pow(x_cv[icv][0]-x_drop,2) + pow(x_cv[icv][1]-y_drop,2) );
        //const double n_true[3] = {(x_cv[icv][0]-x_drop)/r,(x_cv[icv][1]-y_drop)/r,0.0};
        //const double kappa_true = 1.0/r_drop;

        // sphere...
        const double r = sqrt( pow(x_cv[icv][0]-x_drop,2) + pow(x_cv[icv][1]-y_drop,2) + pow(x_cv[icv][2]-z_drop,2) );
        const double n_true[3] = {(x_cv[icv][0]-x_drop)/r,(x_cv[icv][1]-y_drop)/r,(x_cv[icv][2]-z_drop)/r};
        const double kappa_true = 2.0/r_drop;

        kappa_err[icv] = fabs(kappa_true-kappa[icv]);
        //if (kappa_err[icv] > 0.0015) cout << "ocv = " << icv << " " << vof[icv] << " " << cv_flag_real[icv] << " " << COUT_VEC(x_cv[icv]) << endl;
        kappa_ex[icv] = kappa_true;
        n_err[icv] = 1.0-fabs(DOT_PRODUCT(n[icv],n_true)); 
        //cout << "E(x): " << x_cv[icv][0] << " " << vof[icv] << " " << acos(MAX3(abs(n_true[0]),abs(n_true[1]),abs(n_true[2]))) << " " << kappa[icv] << " " << kappa_true << endl;
        
        cnt_int += 1.0;
        L2_kappa += kappa_err[icv]*kappa_err[icv];
        L1_kappa += fabs(kappa_err[icv]);
        Linf_kappa = max(Linf_kappa,fabs(kappa_err[icv]));
        L2_n += n_err[icv]*n_err[icv];
        L1_n += n_err[icv];
        Linf_n = max(Linf_n,n_err[icv]);
      }
      else {
        kappa_err[icv] = 0.0;
        kappa_ex[icv] = 0.0;
        n_err[icv] = 0.0;
      }
    }

    // print total L1/Linf kappa and n error and interface count...
    double my_buf_sum[6] = {L2_kappa,L1_kappa,L2_n,L1_n,cnt_int,vol_1};
    double my_buf_max[2] = {Linf_kappa,Linf_n};
    double buf_sum[6], buf_max[2];
    MPI_Reduce(my_buf_sum,buf_sum,6,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    MPI_Reduce(my_buf_max,buf_max,2,MPI_DOUBLE,MPI_MAX,0,mpi_comm);

    if (mpi_rank == 0) {

      const double kappa_ref = 1.0/r_drop;
      //const double kappa_ref = 2.0/r_drop;
      //const double kappa_ref = 4.0*a0_wave*M_PI*M_PI/l_wave/l_wave;

      const double vol = 4.0*M_PI/3.0*r_drop*r_drop*r_drop;
      cout << "> vol err: " << fabs(buf_sum[5]-vol)/vol << endl;
      cout << "> L2 curvature error: " << sqrt(buf_sum[0]/buf_sum[4]) << endl;
      cout << "> L1 curvature error: " << buf_sum[1]/buf_sum[4]<< endl;
      cout << "> Linf curvature error: " << buf_max[0] << endl;
      cout << "> L2 normal error: " << sqrt(buf_sum[2]/buf_sum[4]) << endl;;
      cout << "> L1 normal error: " << buf_sum[3]/buf_sum[4] << endl;;
      cout << "> Linf normal error: " << buf_max[1] << endl;

    }
    }

    /*
    // velocity measures...
    double umax = 0.0;
    double kinetic = 0.0;
  
    // pressure measures...

    double pmin = +HUGE_VAL;
    double pmax = -HUGE_VAL;
    double pin_tot  = 0.0;
    double pout_tot = 0.0;
    double pin_part  = 0.0;
    double pout_part = 0.0;

    // volume measures...

    double vol = 0.0;
    double vin_tot  = 0.0;
    double vout_tot = 0.0;
    double vin_part  = 0.0;
    double vout_part = 0.0;

    // drop stuff...

    //const double z_drop = 0.0;
    double dp_ex = sigma/r_drop;
    //double dp_ex = 2.0*sigma/r_drop;

    FOR_ICV {
      
      // circle...
      const double r = sqrt( pow(x_cv[icv][0]-x_drop,2) + pow(x_cv[icv][1]-y_drop,2) );

      // sphere...
      //const double r = sqrt( pow(x_cv[icv][0]-x_drop,2) + pow(x_cv[icv][1]-y_drop,2) + pow(x_cv[icv][2]-z_drop,2) );

      const double magu = MAG(u[icv]);
      umax = max(umax,magu);
      kinetic += 0.5*magu*magu*vol_cv[icv];
      
      pmin = min(pmin,p[icv]);
      pmax = max(pmax,p[icv]);

      if (r <= r_drop) {
        vin_tot += vol_cv[icv];
        pin_tot += p[icv]*vol_cv[icv];
        if (r <= 0.5*r_drop) {
          vin_part += vol_cv[icv];
          pin_part += p[icv]*vol_cv[icv];
        }
      }
      else {
        vout_tot += vol_cv[icv];
        pout_tot += p[icv]*vol_cv[icv];
        if (r >= 1.5*r_drop) {
          vout_part += vol_cv[icv];
          pout_part += p[icv]*vol_cv[icv];
        }
      }
      vol += vol_cv[icv];

    }

    double my_buf_max[3] = {umax,pmax,-pmin};
    double my_buf_sum[10] = {kinetic,vol,pin_tot,vin_tot,pout_tot,vout_tot,pin_part,vin_part,pout_part,vout_part};
    double buf_max[3], buf_sum[10];
    MPI_Reduce(my_buf_max,buf_max,3,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    MPI_Reduce(my_buf_sum,buf_sum,10,MPI_DOUBLE,MPI_SUM,0,mpi_comm); 

    if (mpi_rank == 0) {
      cout << "|pin_tot-pout_tot-dp|/dp: " << fabs(buf_sum[2]/buf_sum[3]-buf_sum[4]/buf_sum[5]-dp_ex)/dp_ex << endl;
      cout << "|pin_part-pout_part-dp|/dp: " << fabs(buf_sum[6]/buf_sum[7]-buf_sum[8]/buf_sum[9]-dp_ex)/dp_ex << endl;
      cout << "|pmax-pmin-dp|/dp: " << fabs(buf_max[1]+buf_max[2]-dp_ex)/dp_ex << endl;
      cout << "Lmax_u: " << time << " " <<  buf_max[0] << endl;
      cout << "kinetic energy: " << time << " " <<  buf_sum[0]/buf_sum[1]<< endl; 
    }
*/

  }

};

class DropTransferSolver : public VofSolver {
public:


  DropTransferSolver() {
  }

  void initData() {
    VofSolver::initData();
  }

  ~DropTransferSolver() {
  }

  void initialHook() {

    COUT1("DropTransferSolver::initialHook()");

    int nDrop;
    int nd;
    string line;
    ifstream dropFile;
    double dx[3];

    assert(!checkDataFlag("vof")); 
    assert(!checkDataFlag("u"));
    
    // read random drop placements from file
    if (mpi_rank == 0) {
      dropFile.open("random_drops.txt");
      if (dropFile.is_open()) {
        getline(dropFile, line);
        stringstream ss(line);
        ss >> nDrop;
        cout << "read integer=" << nDrop << endl;
      }
    }
    MPI_Bcast(&nDrop,1,MPI_INT,0,mpi_comm);
    double drop_xyz[nDrop][3];
    double dropR[nDrop];
    if (mpi_rank == 0) {
      for (nd=0; nd<nDrop; nd++) {
        FOR_I3 {
          getline(dropFile,line);
          stringstream ss(line);
          ss >> drop_xyz[nd][i];
        }
        getline(dropFile,line);
        stringstream ss(line);
        ss >> dropR[nd];
      }
    }
    MPI_Bcast(drop_xyz,3*nDrop,MPI_DOUBLE,0,mpi_comm);
    MPI_Bcast(dropR,nDrop,MPI_DOUBLE,0,mpi_comm);
    
    FOR_ICV_G {
      g[icv] = -1.0E10;
      for (nd=0; nd<nDrop; nd++) {
	FOR_I3 dx[i] =  x_cv[icv][i]-drop_xyz[nd][i];
        double gg = dropR[nd] - sqrt(DOT_PRODUCT(dx,dx));
        if (fabs(gg) < fabs(g[icv])) {
          if (!(g[icv] > 0.0 && gg < 0.0)) {
	    g[icv] = gg;
	  }
	} else {
	  if (g[icv] < 0.0 && gg > 0.0) {
	    g[icv] = gg;
	  }
	}
      }
    }
    buildVofFromG();
    calcNormal();
    calcGfromVof();

    FOR_ICV_G {
      u[icv][0] = vof[icv];
      u[icv][1] = 0.0;
      u[icv][2] = 0.0;
    }
  }
  
  void temporalHook() {
    
  }
  

};
class PjetSolver : public VofSolver {

  public:
    /*
  double erf(const double &eta) {
    
    double z = fabs(eta);
    double t = 1.0/(1.0+0.5*z);
    double func = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+
            t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+
                    t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
    if (eta < 0.0) func = 2.0-func;
    func = 1.0 - func;
    return func;
  }
*/
  void initialHook() {

    COUT1("PjetSolver::initialHook()");
    
    assert(!checkDataFlag("vof"));
    assert(!checkDataFlag("u"));

    double Uliq, Ugas, aspect;

    Uliq = getDoubleParam("UL");
    Ugas = getDoubleParam("UG");
    aspect = getDoubleParam("ASPECT");

    COUT1("UL is set to " << Uliq);
    COUT1("UG is set to " << Ugas);
    COUT1("Aspect ratio is set to " << aspect);
    const double r_cyl = 50E-6;
    const double L = 2.0*r_cyl*aspect;
    const double tt = fabs(L/Uliq);
    const double delta_l = 3.64*sqrt(mu_l/rho_l*tt);
    const double delta_g = sqrt(mu_g*rho_l/mu_l/rho_g)*delta_l;

    FOR_ICV_G {
      g[icv] = r_cyl - sqrt(x_cv[icv][2]*x_cv[icv][2]+x_cv[icv][1]*x_cv[icv][1]);
      if (g[icv] >= 0.0) u[icv][0] = Uliq*erf(0.5*3.64*g[icv]/delta_l);
      else               u[icv][0] = Ugas*erf(0.5*3.64*fabs(g[icv])/delta_g);
      u[icv][1] = 0.0;
      u[icv][2] = 0.0;
    }

    buildVofFromG();
    calcNormal();
    calcGfromVof();
    
  }
  
  void temporalHook() {
  } 
  
};


class RayleighSolver : public VofSolver {

public:

  RayleighSolver() {
  }

  void initData() {
    VofSolver::initData();
  }

  ~RayleighSolver() {
  }

  void initialHook() {

    COUT1("RayleighSolver::initialHook()");
    
    if (!checkDataFlag("vof")) {
      assert(!checkDataFlag("u"));

      // translation...
      //const double x_drop = 0.0;
      //const double y_drop = 0.0;
      //const double r_drop = 0.5;
      
      // rotation...
      //const double x_drop = 0.0;
      //const double y_drop = 0.0;
      //const double r_drop = 0.25;
      const double r_cyl = 0.14;
      const double ap = getDoubleParam("AMP");
      const double L = 4.0*M_PI*r_cyl;
      const double lambda = 9.0*r_cyl;
      const double k = 2.0*M_PI/lambda;
      FOR_ICV_G {
        g[icv] = r_cyl - sqrt(x_cv[icv][0]*x_cv[icv][0]+x_cv[icv][1]*x_cv[icv][1]) -  ap*sin(k*x_cv[icv][2]);
        u[icv][0] = 0.0;
        u[icv][1] = 0.0;
        u[icv][2] = 0.0;
      }
 
      buildVofFromG();
      calcNormal();
      calcGfromVof();

    }
    
  }

  void temporalHook() {
    const double eps = getDoubleParam("RESOL");
    if (step%check_interval == 0) {
      
      if (mpi_rank == 0)
        cout << "Rayleigh::temporalHook()" << endl;
      
      // calculate the amplitude at the 3/4 Lambdafirst cell
      const double r_cyl = 0.14;
      const double lambda = 9.0*r_cyl;
      const double target = 0.75*lambda;
      double my_l2[2];
      double l2[2];
      my_l2[0] = my_l2[1] = 0.0;
      l2[0] = l2[1] = 0.0;
      double sum_dist = 0.0;
      FOR_ICV {
        //if (fabs(x_cv[icv][0] - 6.283184/16.0)< 6.283184/16.0/2.0) {
        if (fabs(x_cv[icv][2] - target ) < eps ) {
          if ( iband[icv] == 2 ) {
            double x[3];
            assert(vof[icv] <= 0.5);
            FOR_I3 x[i] = x_cv[icv][i] + g[icv]*n[icv][i];
            my_l2[0] += 2.0*vof[icv]*sqrt(x[0]*x[0] + x[1]*x[1]);
            my_l2[1] += 2.0*vof[icv];
          }
        }
      }
      
      MPI_Reduce(my_l2,l2,2,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      
      if (mpi_rank == 0) {
        assert(l2[1] > 0.0);
        l2[0] = l2[0]/l2[1];
        cout << "Time, Amax " << time << " "  <<  (l2[0]-r_cyl)/r_cyl  << " " << endl;
      }
      
      } 
    } 
  };




class FilamentSolver : public VofSolver {

public:

  FilamentSolver() {
  }

  void initData() {
    VofSolver::initData();
  }

  ~FilamentSolver() {
  }

  void initialHook() {

    COUT1("FilamentSolver::initialHook()");
    
    if (!checkDataFlag("vof")) {
      assert(!checkDataFlag("u"));

      // translation...
      //const double x_drop = 0.0;
      //const double y_drop = 0.0;
      //const double r_drop = 0.5;
      
      // rotation...
      //const double x_drop = 0.0;
      //const double y_drop = 0.0;
      //const double r_drop = 0.25;
      const double r_cyl = 0.4;
      const double ap = r_cyl*0.2;
      const double L = 4.0*M_PI*r_cyl;
      FOR_ICV_G {
        g[icv] = r_cyl - sqrt(x_cv[icv][0]*x_cv[icv][0]+x_cv[icv][1]*x_cv[icv][1])+ap*cos(2.0*M_PI*x_cv[icv][2]/L);
        u[icv][0] = 0.0;
        u[icv][1] = 0.0;
        u[icv][2] = 0.0;
      }
 
      buildVofFromG();
      calcNormal();
      calcGfromVof();


    }
    
  }

  void temporalHook() {
    const double eps = 1.0E-3;
    if (step%check_interval == 0) {
      
      if (mpi_rank == 0)
        cout << "Filament::temporalHook()" << endl;
      
      // calculate the amplitude at the first cell
      const double r_cyl = 0.4;
      const double L = 4.0*M_PI*r_cyl;
      double my_l2[2];
      double l2[2];
      my_l2[0] = my_l2[1] = 1.0E10;
      l2[0] = l2[1] = 0.0;
      FOR_ICV {
        //if (fabs(x_cv[icv][0] - 6.283184/16.0)< 6.283184/16.0/2.0) {
        if (fabs(x_cv[icv][2] )< eps && fabs(x_cv[icv][0])<eps && x_cv[icv][1] > 0.0) {
          if ( vof[icv] > eps && vof[icv] <= 1.0-eps) {
            my_l2[0] = x_cv[icv][1] + g[icv]*n[icv][1];
          }
        }
      }
      
      MPI_Reduce(my_l2,l2,2,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
      
      if (mpi_rank == 0) {
        cout << "Time, Amax " << time << " "  <<
          (l2[0]-r_cyl)/r_cyl  << " " << endl;
      }
      
      
      my_l2[0] = my_l2[1] = 1.0E10;
      l2[0] = l2[1] = 0.0;
      FOR_ICV {
        //if (fabs(x_cv[icv][0] - 6.283184/16.0)< 6.283184/16.0/2.0) {
        if (fabs(x_cv[icv][2]-0.5*L)< eps && fabs(x_cv[icv][0])<eps && x_cv[icv][1] > 0.0) {
          if ( vof[icv] > eps && vof[icv] <= 1.0-eps) {
            my_l2[0] = x_cv[icv][1] + g[icv]*n[icv][1];
          }
        }
      }
      
      
      MPI_Reduce(my_l2,l2,2,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
      
      if (mpi_rank == 0) {
        cout << "Time, Amin " << time << " " <<
          (l2[0]-r_cyl)/r_cyl  << " " << endl;
      }
      
      }
      } 
      
    };



class OscillatingDropVofSolver : public VofSolver {
public:

  // declare data here...

  OscillatingDropVofSolver() {
    // nullify and register data here...
  }

  void initData() {
    VofSolver::initData();
    // allocate data here...
  }

  ~OscillatingDropVofSolver() {
    // delete data here..
  }

  void initialHook() {
    // init vof and u here...

    COUT1("OscillatingDropVofSolver::initialHook()");
    
    if (!checkDataFlag("vof")) {
      assert(!checkDataFlag("u"));

      /*
      const double x_drop = 0.0;
      const double y_drop = 0.0;
      const double z_drop = 0.0;
      const double r0 = 2.0;
      const double a0 = 0.01*r0;
      //const double a0 = 0.0;
      FOR_ICV_G {
        
        const double r = sqrt( pow(x_cv[icv][0]-x_drop,2) + pow(x_cv[icv][1]-y_drop,2));
        const double cos2theta = (x_cv[icv][0]*x_cv[icv][0] - x_cv[icv][1]*x_cv[icv][1])/(x_cv[icv][0]*x_cv[icv][0] + x_cv[icv][1]*x_cv[icv][1]);
        const double r_drop = r0 + a0*cos2theta;
        
        FOR_I3 u[icv][i] = 0.0;
        g[icv] = r_drop-r;
      }
      buildVofFromG();
      */
      
      }
    
  }

  void temporalHook() {
   
    double my_tke = 0.0;
    FOR_ICV my_tke += 0.5*DOT_PRODUCT(u[icv],u[icv]);
    double tke;
    MPI_Reduce(&my_tke,&tke,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm); 
    if (mpi_rank == 0) {
      cout << "TKE: " << tke << endl;
    }
  
  };
  void finalHook() {};

};

class MyVofSolver : public VofSolver {
public:

  // declare data here...

  MyVofSolver() {
    // nullify and register data here...
  }

  void initData() {
    VofSolver::initData();
    // allocate data here...
  }

  ~MyVofSolver() {
    // delete data here..
  }

  void initialHook() {
    // init vof and u here...

    COUT1("MyVofSolver::initialHook()");
     
    if (!checkDataFlag("vof")) {
      assert(!checkDataFlag("u"));
      FOR_ICV vof[icv] = 0.0;
      FOR_ICV FOR_I3 u[icv][i] = 0.0;
    }
    
  }

  void temporalHook() {};
  void finalHook() {};

};

int main(int argc, char* argv[]) {
  
  try {
   
    CTI_Init(argc,argv,"vof.in");

    const bool b_post = checkParam("POST");

    { 
      
      if (Param * param = getParam("TEST_CASE")) { 
        const string test = param->getString();
        if ( test == "POOL") {
          PoolVofSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "DROP_POOL") {
          DropPoolVofSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "STATIC_DROP") {
          StaticDropVofSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "RTI" ) {
          RTIVofSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "KINEMATIC" ) {
          KinematicVofSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "OSCILLATING_DROP") {
          OscillatingDropVofSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "STATIC") {
          StaticVofSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "CONVECTION") {
          ConvectionVofSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "LIQUID_JET") {
          LiquidjetVofSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "SURFACE_WAVE") {
          SurfaceWaveVofSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "DROP_TRANSFER") {
          DropTransferSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "HSQUARE") {
          HsquareSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "TEST_VOF") {
          TestVofSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "ZALESAKS") {
          ZalesaksSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "M_ZALESAKS") {
          MZalesaksSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "FILAMENT") {
          FilamentSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "PJET") {
          PjetSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "DROP_BREAKUP") {
          DropBreakSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "RAYLEIGH") {
          RayleighSolver solver;
          solver.init();
          solver.run();
        }
        else if ( test == "Sphere3D") {
          Sphere3DVofSolver solver;
          solver.init();
          solver.run();
        } 
        else { 
          CERR("unrecognized TEST_CASE: " << test << ", possible choices are \n \
                  TEST_CASE POOL\n					\
                  TEST_CASE DROP_POOL\n					\
                  TEST_CASE STATIC_DROP\n				\
                  TEST_CASE RTI\n					\
                  TEST_CASE KINEMATIC\n					\
                  TEST_CASE OSCILLATING_DROP\n                          \
                  TEST_CASE STATIC\n");
	}
      }
      else {
        VofSolver solver;
        //MyVofSolver solver;
        //
        if (b_post) {
          solver.initMin();
          solver.runPost();
        }
        else {
          solver.init();
          solver.run();
        }
      }
    
    } 
    
    CTI_Finalize();

  } 
  catch (int e) {
    if (e == 0)
      CTI_Finalize();
    else
      CTI_Abort();
  } 
  catch(...) {
    CTI_Abort();
  }
  
  return 0;

} 
