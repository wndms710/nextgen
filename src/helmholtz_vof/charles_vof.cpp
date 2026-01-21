#include "HelmholtzVofSolver.hpp"
#include "BasicPostpro.hpp"
 

class HsquareSolver : public HelmholtzVofSolver {
public:

  double* vof_ini;
  double* vof_ex;
  HsquareSolver() {
    vof_ini = NULL; 
    vof_ex = NULL; registerCvData(vof_ex,"vof_ex",CAN_WRITE_DATA);
  }

  void initData() {
    HelmholtzVofSolver::initData();
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


class SquareSolver : public HelmholtzVofSolver {
public:

  double* vof0;
  SquareSolver() {
    vof0 = NULL; 
  }

  void initData() {
    HelmholtzVofSolver::initData();
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

class CircleSolver : public HelmholtzVofSolver {
public:

  double* vof0;
  CircleSolver() {
    vof0 = NULL; 
  }

  void initData() {
    HelmholtzVofSolver::initData();
    vof0 = new double[ncv_g];
  }

  ~CircleSolver() {
    DELETE(vof0);
  }

  void initialHook() {

    COUT1("CircleSolver::initialHook()");
    
    if (!checkDataFlag("vof")) {  
      FOR_ICV_G {
        u[icv][0] = 0.5-x_cv[icv][1];
        u[icv][1] = x_cv[icv][0]-0.5;
        u[icv][2] = 0.0;
      }
      
      
      FOR_IFA {
        const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
        const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
        double u_fa[3]; FOR_I3 u_fa[i] = 0.5*(u[icv0][i]+u[icv1][i]);
        q_fa[ifa] = DOT_PRODUCT(u_fa,n_fa[ifa]);
      }
      
      const double x_disk[2] = { 0.5, 0.75};
      const double init_radius = 0.15;
      
      
      FOR_ICV_G {
        double dx[2] = { x_cv[icv][0]-x_disk[0], x_cv[icv][1]-x_disk[1] };
        g[icv] = init_radius - sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
      }
      
      //buildVofFromG();
      //calcNormal();
      //calcGfromVof();
      
      // vof subdivision...
      double x[4][3];


      for (int icv = 0; icv < ncv; ++icv) {
        
        vof[icv] = 0.0;

        for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
          const int ifa = faocv_v[foc];
          
          int ino0 = noofa_v[noofa_i[ifa+1]-1];
          for (int nof1 = noofa_i[ifa]; nof1 != noofa_i[ifa+1]; ++nof1) {
            const int ino1 = noofa_v[nof1];
            
            // we need to order the nodes according to dx dot n...
            
            FOR_J3 x[0][j] = x_cv[icv][j];
            FOR_J3 x[1][j] = x_no[ino0][j];
            FOR_J3 x[2][j] = x_no[ino1][j];
            FOR_J3 x[3][j] = x_fa[ifa][j];
            vof[icv] += addSubVols(x,0,init_radius);

            ino0 = ino1;
          }
          
        }
        vof[icv] *= inv_vol[icv];
        vof[icv] = min(1.0,max(0.0,vof[icv])); // get rid of eps's
      }

    }

    calcGfromVof();
    calcNormal();
  
    FOR_ICV vof0[icv] = vof[icv];
    
  }

  void temporalHook() {

    int vof_start = getDoubleParam("VOF_START");
    if (step == vof_start ) FOR_ICV_G  vof0[icv] = vof[icv];

     FOR_ICV_G {
      u[icv][0] = 0.5-x_cv[icv][1];
      u[icv][1] = x_cv[icv][0]-0.5;
      u[icv][2] = 0.0;
    }
    
    // solvePAndCorrectU();
    
    
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
        
        cout << " > Shape error: " << buf_sum[0] << endl;
        cout << " > Volume error: " << buf_sum[1] << endl;
        cout << " > Total volume: " << buf_sum[2] << endl;
        
      }
      
      
  }

  void reorderTetByS(double s[4],double dx[4][3]) {

    double dxtemp[4][3];
    double stemp[4];

    // copy the info for reordering ...

    FOR_I4 FOR_J3 dxtemp[i][j] = dx[i][j];
    FOR_I4 stemp[i] = s[i];

    // start with first point in dx[0][i] ...

    s[0] = stemp[0];
    FOR_J3 dx[0][j] = dxtemp[0][j];

    // second point of subtet ...

    double this_s = stemp[1];
    int i;
    for (i = 1; i > 0; --i) {
      if (this_s >= s[i-1]) break;
      FOR_J3 dx[i][j] = dx[i-1][j];
      s[i] = s[i-1];
    }
    FOR_J3 dx[i][j] = dxtemp[1][j];
    s[i] = this_s;

    // third point...

    this_s = stemp[2];
    for (i = 2; i > 0; --i) {
      if (this_s >= s[i-1]) break;
      FOR_J3 dx[i][j] = dx[i-1][j];
      s[i] = s[i-1];
    }
    FOR_J3 dx[i][j] = dxtemp[2][j];
    s[i] = this_s;

    // fourth point...

    this_s = stemp[3];
    for (i = 3; i > 0; --i) {
      if (this_s>= s[i-1]) break;
      FOR_J3 dx[i][j] = dx[i-1][j];
      s[i] = s[i-1];
    }
    FOR_J3 dx[i][j] = dxtemp[3][j];
    s[i] = this_s;

  }

  double calcIntervalVolume(const double this_ds,const double ds,const double A0,const double A1,const double dV) {

    //assert(this_ds >= 0.0);
    //assert(ds > 0.0);

    const double c2 = (3.0*dV/ds - 2.0*A0 - A1)/ds; // quadratic coeff
    const double c3 = (A0 + A1 - 2.0*dV/ds)/(ds*ds); // cubic coeff

    return( (A0 + (c2 + c3*this_ds)*this_ds)*this_ds );

  }

  double calcOrderedTetVofVolume(const double s[4],const double dx[4][3],const double normal[3]) {

    /*
    FOR_I4 cout << s[i] << " ";
    FOR_I4 {
      cout << endl;
      FOR_J3 cout << dx[i][j] << " ";
    }
    cout << endl;
    */

    // double-check ordering...
    //assert(s[1] >= s[0]);
    //assert(s[2] >= s[1]);
    //assert(s[3] >= s[2]);

    // we can either work from the minus g side or from the plus g side towards g = 0.
    // Efficiency depends on where the zero-crossing occurs...

    if (s[3] <= 0.0) {
      // all g's are negative, so volume above g = 0 is formally 0...
      return(0.0);
    }
    else {
      const double tet_volume = fabs(SIGNED_TET_VOLUME_6(dx[0],dx[1],dx[2],dx[3]))/6.0;
      if (s[0] >= 0.0) {
        // all g's are positive, so volume above g = 0 is the full volume...
        return(tet_volume);
      }
      else if (s[1] >= 0.0) {
        //assert(s[0] < 0.0);
        // most of the g's are positive, so calculate the minus area, then reverse...
        double v1_02[3]; FOR_I3 v1_02[i] = ((s[2]-s[1])*dx[0][i] + (s[1]-s[0])*dx[2][i])/(s[2]-s[0]) - dx[1][i];
        double v1_03[3]; FOR_I3 v1_03[i] = ((s[3]-s[1])*dx[0][i] + (s[1]-s[0])*dx[3][i])/(s[3]-s[0]) - dx[1][i];
        double A = 0.5*fabs(CROSS_DOT(v1_02,v1_03,normal));
        double V = A*(s[1]-s[0])/3.0;
        return(tet_volume - calcIntervalVolume(-s[0],s[1]-s[0],0.0,A,V));
      }
      else if (s[2] <= 0.0) {
        //assert(s[3] > 0.0);
        // most of the g's are negative, so calculate the little bit of positive area directly...
        double v2_03[3]; FOR_I3 v2_03[i] = ((s[3]-s[2])*dx[0][i] + (s[2]-s[0])*dx[3][i])/(s[3]-s[0]) - dx[2][i];
        double v2_13[3]; FOR_I3 v2_13[i] = ((s[3]-s[2])*dx[1][i] + (s[2]-s[1])*dx[3][i])/(s[3]-s[1]) - dx[2][i];
        double A = 0.5*fabs(CROSS_DOT(v2_03,v2_13,normal));
        double V = A*(s[3]-s[2])/3.0;
        return(calcIntervalVolume(s[3],s[3]-s[2],0.0,A,V));
      }
      else {
        // intersection is between s[1] and s[2]...
        //assert(s[1] < 0.0);
        //assert(s[2] > 0.0);
        double v1_02[3]; FOR_I3 v1_02[i] = ((s[2]-s[1])*dx[0][i] + (s[1]-s[0])*dx[2][i])/(s[2]-s[0]) - dx[1][i];
        double v1_03[3]; FOR_I3 v1_03[i] = ((s[3]-s[1])*dx[0][i] + (s[1]-s[0])*dx[3][i])/(s[3]-s[0]) - dx[1][i];
        double A1 = 0.5*fabs(CROSS_DOT(v1_02,v1_03,normal));
        double V1 = A1*(s[1]-s[0])/3.0;
        double v2_03[3]; FOR_I3 v2_03[i] = ((s[3]-s[2])*dx[0][i] + (s[2]-s[0])*dx[3][i])/(s[3]-s[0]) - dx[2][i];
        double v2_13[3]; FOR_I3 v2_13[i] = ((s[3]-s[2])*dx[1][i] + (s[2]-s[1])*dx[3][i])/(s[3]-s[1]) - dx[2][i];
        double A2 = 0.5*fabs(CROSS_DOT(v2_03,v2_13,normal));
        double V2 = A2*(s[3]-s[2])/3.0;
        return(V2 + calcIntervalVolume(s[2],s[2]-s[1],A2,A1,tet_volume-V1-V2));
      }
    }
  }

  double addSubVols(const double x_p[4][3],const int level_p, const double r_drop) {
    const int max_level = 0; // 4

    // get signed distances...
    const double x_drop = 0.5;
    const double y_drop = 0.75;
    const double z_drop = 0.5;
    double s[4];
    bool all_pos = true;
    bool all_neg = true;

    FOR_I4 {

      // circle...
       s[i] = r_drop - sqrt( pow(x_p[i][0]-x_drop,2) + pow(x_p[i][1]-y_drop,2) );

      // sphere...
      //s[i] = r_drop - sqrt( pow(x_p[i][0]-x_drop,2) + pow(x_p[i][1]-y_drop,2) + pow(x_p[i][2]-z_drop,2) );

      if (all_pos && s[i] < 0.0) all_pos = false;
      else if (all_neg && s[i] > 0.0) all_neg = false;
    }
    if (all_pos) {
      return fabs(SIGNED_TET_VOLUME_6(x_p[0],x_p[1],x_p[2],x_p[3]))/6.0;
    }
    else if (all_neg) {
      return 0.0;
    }
    else if (level_p < max_level) {
      double vol = 0.0;
      double x[4][3];
      FOR_I3 x[0][i] = 0.25*(x_p[0][i]+x_p[1][i]+x_p[2][i]+x_p[3][i]);
      FOR_J4 {
        FOR_I3 x[1][i] = x_p[(j  )%4][i];
        FOR_I3 x[2][i] = x_p[(j+1)%4][i];
        FOR_I3 x[3][i] = x_p[(j+2)%4][i];
        vol += addSubVols(x,level_p+1,r_drop);
      }
      return vol;
    }
    else {
      double xc[3]; FOR_I3 xc[i] = 0.25*(x_p[0][i]+x_p[1][i]+x_p[2][i]+x_p[3][i]);

      double x0[3]; FOR_I3 x0[i] = (x_p[1][i]+x_p[2][i]+x_p[3][i])/3.0;
      const double dx0[3] = DIFF(x0,xc);
      double n0[3] = TRI_NORMAL_2(x_p[1],x_p[2],x_p[3]);
      double sgn0 = SGN(DOT_PRODUCT(n0,dx0));
      FOR_I3 n0[i] *= sgn0;

      double x1[3]; FOR_I3 x1[i] = (x_p[0][i]+x_p[2][i]+x_p[3][i])/3.0;
      const double dx1[3] = DIFF(x1,xc);
      double n1[3] = TRI_NORMAL_2(x_p[2],x_p[3],x_p[0]);
      double sgn1 = SGN(DOT_PRODUCT(n1,dx1));
      FOR_I3 n1[i] *= sgn1;

      double x2[3]; FOR_I3 x2[i] = (x_p[0][i]+x_p[1][i]+x_p[3][i])/3.0;
      const double dx2[3] = DIFF(x2,xc);
      double n2[3] = TRI_NORMAL_2(x_p[3],x_p[0],x_p[1]);
      double sgn2 = SGN(DOT_PRODUCT(n2,dx2));
      FOR_I3 n2[i] *= sgn2;

      double x3[3]; FOR_I3 x3[i] = (x_p[0][i]+x_p[1][i]+x_p[2][i])/3.0;
      const double dx3[3] = DIFF(x3,xc);
      double n3[3] = TRI_NORMAL_2(x_p[0],x_p[1],x_p[2]);
      double sgn3 = SGN(DOT_PRODUCT(n3,dx3));
      FOR_I3 n3[i] *= sgn3;

      const double vol = fabs(SIGNED_TET_VOLUME_6(x_p[0],x_p[1],x_p[2],x_p[3]));
      double grad_s[3];
      FOR_I3 grad_s[i] = -(n0[i]*s[0]+n1[i]*s[1]+n2[i]*s[2]+n3[i]*s[3])/vol;
      const double mag = MAG(grad_s);
      FOR_I3 grad_s[i] /= -mag;

      double my_x_p[4][3]; FOR_J4 FOR_I3 my_x_p[j][i] = x_p[j][i]; // need local copy
      reorderTetByS(s,my_x_p);
      return calcOrderedTetVofVolume(s,my_x_p,grad_s);
    }

  }
};


class ZalesaksSolver : public HelmholtzVofSolver {
public:

  double* vof0;
  ZalesaksSolver() {
    vof0 = NULL; 
  }

  void initData() {
    HelmholtzVofSolver::initData();
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


    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      double u_fa[3]; FOR_I3 u_fa[i] = 0.5*(u[icv0][i]+u[icv1][i]);
      q_fa[ifa] = DOT_PRODUCT(u_fa,n_fa[ifa]);
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

class MZalesaksSolver : public HelmholtzVofSolver {
public:

  double* vof0;
  MZalesaksSolver() {
    vof0 = NULL; 
  }

  void initData() {
    HelmholtzVofSolver::initData();
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

    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      double u_fa[3]; FOR_I3 u_fa[i] = 0.5*(u[icv0][i]+u[icv1][i]);
      q_fa[ifa] = DOT_PRODUCT(u_fa,n_fa[ifa]);
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
    
     FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      double u_fa[3]; FOR_I3 u_fa[i] = 0.5*(u[icv0][i]+u[icv1][i]);
      q_fa[ifa] = DOT_PRODUCT(u_fa,n_fa[ifa]);
    }
    
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



class KinematicHelmholtzVofSolver : public HelmholtzVofSolver {
public:

  double* vof0;
  KinematicHelmholtzVofSolver() {
    vof0 = NULL; 
  }

  void initData() {
    HelmholtzVofSolver::initData();
    vof0 = new double[ncv];
  }

  ~KinematicHelmholtzVofSolver() {
    DELETE(vof0);
  }

  void initialHook() {

    COUT1("KinematicHelmholtzVofSolver::initialHook()");
    
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


class Sphere3DHelmholtzVofSolver : public HelmholtzVofSolver {
public:

  double* vof0;
  Sphere3DHelmholtzVofSolver() {
    vof0 = NULL; 
  }

  void initData() {
    HelmholtzVofSolver::initData();
    vof0 = new double[ncv];
  }

  ~Sphere3DHelmholtzVofSolver() {
    DELETE(vof0);
  }

  void initialHook() {

    COUT1("Sphere3DHelmholtzVofSolver::initialHook()");
    
    assert(!checkDataFlag("u"));
    assert(!checkDataFlag("vof"));
      
     
    const double radius = 0.15;
    const double xc[3] = {0.35, 0.35, 0.35};
    
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
      u[icv][0] = 2.0*sin(pi*x)*sin(pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)*cos(pi*time/T);
      u[icv][1] =   - sin(2.0*pi*x)*sin(pi*y)*sin(pi*y)*sin(2.0*pi*z)*cos(pi*time/T);
      u[icv][2] =    -sin(2.0*pi*x)*sin(2.0*pi*y)*sin(pi*z)*sin(pi*z)*cos(pi*time/T);
    }

    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      double u_fa[3]; FOR_I3 u_fa[i] = 0.5*(u[icv0][i]+u[icv1][i]);
      q_fa[ifa] = DOT_PRODUCT(u_fa,n_fa[ifa]);
    }
   
    FOR_ICV vof0[icv] = vof[icv];

  }

  void temporalHook() {

    FOR_ICV_G {
      const double x = x_cv[icv][0];
      const double y = x_cv[icv][1];
      const double z = x_cv[icv][2];
      const double pi = M_PI;
      const double T = 3.0;
      u[icv][0] = 2.0*sin(pi*x)*sin(pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)*cos(pi*time/3.0);
      u[icv][1] =   - sin(2.0*pi*x)*sin(pi*y)*sin(pi*y)*sin(2.0*pi*z)*cos(pi*time/3.0);
      u[icv][2] =    -sin(2.0*pi*x)*sin(2.0*pi*y)*sin(pi*z)*sin(pi*z)*cos(pi*time/3.0);
    }
    
    FOR_IFA {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv_g));
      double u_fa[3]; FOR_I3 u_fa[i] = 0.5*(u[icv0][i]+u[icv1][i]);
      q_fa[ifa] = DOT_PRODUCT(u_fa,n_fa[ifa]);
    }
   
 //   solvePAndCorrectU();

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


class StaticHelmholtzVofSolver : public HelmholtzVofSolver {
public:

  double* kappa_err;
  double* kappa_ex;
  double* n_err;

  StaticHelmholtzVofSolver() {
    kappa_err = NULL; registerCvData(kappa_err,"kappa_err",0);
    kappa_ex = NULL; registerCvData(kappa_ex,"kappa_ex",0);
    n_err = NULL; registerCvData(n_err,"n_err",0);
  }

  void initData() {
    HelmholtzVofSolver::initData();
    kappa_err = new double[ncv];
    kappa_ex = new double[ncv];
    n_err = new double[ncv];
  }

  ~StaticHelmholtzVofSolver() {
    DELETE(kappa_err);
    DELETE(kappa_ex);
    DELETE(n_err);
  }

  void initialHook() {

    COUT1("StaticHelmholtzVofSolver::initialHook()");
    
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


class OscillatingDropHelmholtzVofSolver : public HelmholtzVofSolver {
public:

  // declare data here...

  OscillatingDropHelmholtzVofSolver() {
    // nullify and register data here...
  }

  void initData() {
    HelmholtzVofSolver::initData();
    // allocate data here...
  }

  ~OscillatingDropHelmholtzVofSolver() {
    // delete data here..
  }

  void initialHook() {
    // init vof and u here...

    COUT1("OscillatingDropHelmholtzVofSolver::initialHook()");
    
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

class MyHelmholtzVofSolver : public HelmholtzVofSolver {
public:

  // declare data here...

  MyHelmholtzVofSolver() {
    // nullify and register data here...
  }

//  void initData() {
//    HelmholtzVofSolver::initData();
//  }

  ~MyHelmholtzVofSolver() {
    // delete data here..
  }

  void initialHook() {

    string test_case;

    if ( Param * param  = getParam("TEST_CASE")) { 

      test_case = param->getString(0);
      if ( mpi_rank == 0 ) 
        cout << " > starting test case: " << test_case << endl;
    }

    if ( step == 0) { 
     
      if ( test_case == "DAMPED_SURFACE_WAVE") { 
        initialHookSurfaceWave();
      } else if (test_case == "SPRAY_A") {
        initialHookLiquidJet();
      } else if (test_case == "RTI") {
        initialHookRTI();
      } else if (test_case == "DROP") {
        initialHookDrop();
      } else { 
        CWARN(" > initial hook not setting .. ");
      }
    }

  }

  void temporalHook() {

    string test_case;

    if ( Param * param  = getParam("TEST_CASE")) { 
      test_case = param->getString(0);
      
      if ( test_case == "DAMPED_SURFACE_WAVE") { 
        temporalHookSurfaceWave();
      } 
      else if (test_case == "RTI") {
        temporalHookRTI();
      }
    }

  }

  void initialHookSurfaceWave(){

    if (step == 0) {

      COUT1("SurfaceWaveHelmholtzVofSolver:initialHook()");

      // check min X position ...
      double myXmin = 1.0E10;
      FOR_ICV {
        myXmin = min(x_cv[icv][0],myXmin);
      }
      double Xmin;
      
      MPI_Allreduce(&myXmin,&Xmin,1,MPI_DOUBLE,MPI_MIN,mpi_comm);

      if (mpi_rank == 0 ) cout << "Xmin global = " << Xmin << endl;


      const double init_center[2] = { Xmin, M_PI};
      const double init_amplitude = 0.01*2.0*M_PI;
      const double init_wavelength =  2.0*M_PI;

      FOR_ICV_G {
        g[icv] = x_cv[icv][1] - init_center[1] + init_amplitude*cos(x_cv[icv][0]-init_center[0]);
        u[icv][0] = 0.0;
        u[icv][1] = 0.0;
        u[icv][2] = 0.0;
      }

      buildVofFromG();
      calcNormal();
      calcGfromVof();

      double ymin = 1.0E10;
      FOR_IFA {
        const int icv0 = cvofa[ifa][0];
        const int icv1 = cvofa[ifa][1];
        if ( (fabs( x_cv[icv0][0] - Xmin) < 1.0E-10 ) && (fabs(x_cv[icv1][0] -Xmin) < 1.0E-10) ) {
          if ( g[icv0]*g[icv1] <= 0.0) {
            ymin = (fabs(g[icv0])*x_cv[icv1][1] + fabs(g[icv1])*x_cv[icv0][1] ) / (fabs(g[icv0]) + fabs(g[icv1]));
          }
        }
      }

      double ymin_global;
      MPI_Allreduce(&ymin,&ymin_global,1,MPI_DOUBLE,MPI_MIN,mpi_comm);
      if (mpi_rank == 0 ) cout << " time, amplitude = " << time << " " << -(ymin_global-M_PI)/2.0/M_PI << endl;
    }
      
  }
 
  void temporalHookSurfaceWave(){

    if (step%check_interval == 0 ) {

      COUT1(" > temporalHookSurfaceWave()");
      // check min X position ...
      double myXmin = 1.0E10;
      FOR_ICV {
        myXmin = min(x_cv[icv][0],myXmin);
      }
      double Xmin;
      
      MPI_Allreduce(&myXmin,&Xmin,1,MPI_DOUBLE,MPI_MIN,mpi_comm);


      double ymin = 1.0E10;
      FOR_IFA {
        const int icv0 = cvofa[ifa][0];
        const int icv1 = cvofa[ifa][1];
        if ( (fabs( x_cv[icv0][0] - Xmin) < 1.0E-10 ) && (fabs(x_cv[icv1][0] -Xmin) < 1.0E-10) ) {
          if ( g[icv0]*g[icv1] <= 0.0) {
            ymin = (fabs(g[icv0])*x_cv[icv1][1] + fabs(g[icv1])*x_cv[icv0][1] ) / (fabs(g[icv0]) + fabs(g[icv1]));
          }
        }
      }

      double ymin_global;
      MPI_Allreduce(&ymin,&ymin_global,1,MPI_DOUBLE,MPI_MIN,mpi_comm);
      if (mpi_rank == 0 ) cout << " time, amplitude = " << time << " " << -(ymin_global-M_PI)/2.0/M_PI << endl;
    }
     
  }

  void initialHookLiquidJet() {

    COUT1("LiquidjetHelmholtzVofSolver::initialHook()");
    
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
    FOR_ICV {
      rho0[icv] = rho0_ref;
      rho1[icv] = rho1_ref;
      rhoY0[icv] = rho0_ref*vof[icv];
      rho[icv] = rho0_ref*vof[icv] + rho1_ref*(1.0-vof[icv]);
      p[icv] = p_ref;
    }

    updateCvData(rho0);
    updateCvData(rho1);
    updateCvData(p);
    updateCvData(rhoY0);
    updateCvData(rho);

  }

  void initialHookRTI(){

    if (step == 0) {

      COUT1("RTIHelmholtzVofSolver:initialHook()");

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

      FOR_ICV_G {
        rho0[icv] = rho0_ref;
        rho1[icv] = rho1_ref;
        rhoY0[icv] = rho0[icv]*vof[icv];
        rho[icv] = rho0[icv]*vof[icv] +  rho1[icv]*(1.0-vof[icv]);
        p[icv] = p_ref;
      }


    }
      
  }

  void temporalHookRTI(){


  }

  void initialHookDrop(){

    if (step == 0) {

      COUT1("DropHelmholtzVofSolver:initialHook()");

      const double radius = getDoubleParam("RADIUS", 0.1);
      const double xc[3] = {0.0,0.0,0.0};

      if (mpi_rank == 0 ) cout << " >> RADIUS = " << radius << endl; 

      FOR_ICV_G {
        double dx[3] = { x_cv[icv][0]-xc[0], x_cv[icv][1]-xc[1], x_cv[icv][2] - xc[2] };
        g[icv] = radius - sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
        u[icv][0] = 0.0;
        u[icv][1] = 0.0;
        u[icv][2] = 0.0;
      }

      buildVofFromG();
      calcNormal();
      calcGfromVof();

     
    } 

  }

};

int main(int argc, char* argv[]) {
  
  try {
   
    CTI_Init(argc,argv,"charles_vof.in");

    const bool b_post = checkParam("POST");
    const bool test_case = checkParam("TEST_CASE");

    if (b_post) {
      HelmholtzVofSolver solver;
      solver.initMin();
      solver.runPost();
    }
    else if (test_case) {
      MyHelmholtzVofSolver solver;
      solver.init();
      solver.run();
    }
    else {
      HelmholtzVofSolver solver;
      solver.init();
      solver.run();
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
