
#include "CTI.hpp"
using namespace CTI;

#include "StaticSolver.hpp"
#include <complex>

typedef complex<double> dcmplx;

class ComplexComm { 
public: 
  dcmplx* var;
  ComplexComm() : var(NULL) {}
  void setVar(dcmplx* _var) { var = _var; }
  static int data_size() { return 2;}
};

inline void pack_class(double* buf, const ComplexComm* s, const int icv) { 
  buf[0] = s->var[icv].real();
  buf[1] = s->var[icv].imag();
}

inline void pack_class(double* buf, const ComplexComm* s, const int icv, 
                       const double* R) { 

  buf[0] = s->var[icv].real();
  buf[1] = s->var[icv].imag();
}

inline void unpack_class(double* buf, const ComplexComm* s, const int icv) { 

  s->var[icv] = dcmplx(buf[0],buf[1]);
}

class Resonance : public StaticSolver { 
public: 

  std::fstream sweepFile;
  dcmplx * A_op;
  double * sos;

  // the following are used to try to determine a speed of sound..
  string T_type;
  double * T_avg;
  string sos_type;
  double * sos_avg;
  
  dcmplx * psi;

  double * mode_real;
  double * mode_imag;

  Resonance() { 

    if ( mpi_rank == 0) 
      cout << "Resonance()" << endl;

    A_op = NULL;
    sos  = NULL;
    psi  = NULL;

    if ( checkParam("USE_INSTANTANEOUS_SOS"))
      sos_type = "sos";
    else
      sos_type = "sos_avg";
    
    if ( checkParam("USE_INSTANTANEOUS_T"))
      T_type = "T";
    else 
      T_type = "T_avg";

    sos_avg = NULL; registerCvData(sos_avg,sos_type,READWRITE_DATA);
    T_avg = NULL; registerCvData(T_avg,T_type, READWRITE_DATA);

    mode_real = NULL; registerCvData(mode_real,"mode_real", READWRITE_DATA);
    mode_imag = NULL; registerCvData(mode_imag,"mode_imag", READWRITE_DATA);

  } 

  void init() { 

    if ( mpi_rank == 0) 
      logger->setKillFilename("killcharles");

    StaticSolver::init(INIT_COMPACT_FACES|INIT_CV_GRAD);

  }

  void initData() {
  
    assert( T_avg == NULL); T_avg = new double[ncv];
    assert( sos_avg == NULL); sos_avg = new double[ncv];

    assert( mode_real == NULL); mode_real = new double[ncv_g];
    assert( mode_imag == NULL); mode_imag = new double[ncv_g];
  }

  virtual void set_speed_of_sound(double* &sos) { 

    if ( sos == NULL)  
      sos = new double[ncv_g];

    if ( CtiRegister::checkDataFlag(sos_type)) { 

      if ( mpi_rank == 0) { 
        cout << " > using the sos from the result file: " << sos_type << endl;
      }

      for (int icv = 0; icv < ncv; ++icv) 
        sos[icv] = sos_avg[icv];
    
    } 
    else if ( CtiRegister::checkDataFlag(T_type)) { 

      if ( mpi_rank == 0) { 
        cout << " > approximating sos from T_avg " << endl;
        cout << " > sos_avg = sqrt(gamma*R*T_avg) " << endl;
        cout << " > checking for params: " << endl;
        cout << "     - GAMMA = [1.4] " << endl;
        cout << "     - R_GAS = [287 #J/kg/K]" << endl;
      }

      const double gamma = getDoubleParam("GAMMA", 1.4);
      const double R_gas = getDoubleParam("R_GAS", 287.05);

      for (int icv = 0; icv < ncv; ++icv) 
        sos[icv] = sqrt(gamma*R_gas*T_avg[icv]);

    } 
    else { 

      if ( mpi_rank == 0) { 
        cout << " > unable to find sos_avg or T_avg from the result file" << endl;
        cout << " > speed of sound now set from a constant value : " << endl;
        cout << "    - CONSTANT_SOS = [1.0] " << endl;
      }

      const double sos_nom = getDoubleParam("CONSTANT_SOS", 1.0);

      for (int icv = 0; icv < ncv; ++icv) 
        sos[icv] = sos_nom;

    }
    
    updateCvData(sos);
    MiscUtils::dumpRange(sos,ncv,"sos");

  } 

  void run() { 

    if ( mpi_rank == 0) 
      cout << " Resonance::run() " << endl;

    set_speed_of_sound(sos);
    
    if ( Param* param = getParam("SWEEP_FREQUENCY")) { 

      double freq_start, freq_end, df;
      assert( param->getString(0) == "START");
      freq_start = param->getDouble(1);
      assert( param->getString(2) == "DELTA");
      df         = param->getDouble(3);
      assert( param->getString(4) == "END");
      freq_end   = param->getDouble(5);
  
      if (mpi_rank == 0) {
        sweepFile.open("sweep.modes", fstream::out);
        sweepFile << "#Frequency [Hz] Gain [-]" << endl;
        cout << "tsup" << endl;
      }

      double freq = freq_start;
      while ( freq <= freq_end) { 
        solve_freq(freq);
        freq += df;
      }
      if (mpi_rank == 0) sweepFile.close();

    } else if ( Param * param = getParam("FIND_MODE")) { 
    
      assert( param->getString(0) == "AT_FREQUENCY");
      const double freq = param->getDouble(1); 
      compute_mode(freq);

      if ( mpi_rank == 0) 
        cout << " > finished computing mode -- processStep/interactive" << endl;

      if (mpi_rank==0) logger->insertSolverAction(SIMULATION_RUN);
      processStep(0);

      if (mpi_rank == 0) 
        cout << " > writing a result and finishing .. " << endl;
       
      writeResult(0);

    } else if (checkParam("INSPECT")) { 
      
      if (mpi_rank==0) logger->insertSolverAction(SIMULATION_RUN);
      processStep(0);
      
    } else { 

      CERR(" > Syntax error in run. Looking for: \n \
             > For sweeping through frequency space \n \n \
                SWEEP_FREQUENCY START <freq-start> DELTA <delta-f> END <freq-end> \n \
                X_SEED <x0> <y0> <z0> \n \n \
             > or to identify an eigenfunction: \n \
                FIND_MODE AT_FREQUENCY <freq0> \n \
                DF_RANGE <freq_int> # within a certain Delta(f) of freq0 \n \
               ");

    }

  }

  void compute_mode(const double freq0) { 

    // for the following -- we are assuming that the boundary condition 
    // is simple resulting in the ability to write the helmholtz operator 
    // as a real, symmetric matrix.  as a result, these bcs result in a 
    // multiplicity of one for the approximate frequency given 
    
    double freq = freq0;

    assert ( mode_real != NULL);
    assert ( mode_imag != NULL);

    double my_sum = 0.0;

    for (int icv = 0; icv < ncv; ++icv) { 

      mode_imag[icv] = 0.0; // unsued in the following.

      // set an initial vector... 

      mode_real[icv] = 2.0*double(rand())/double(RAND_MAX) -1.0;
      my_sum        += mode_real[icv]*vol_cv[icv]*mode_real[icv];

    } 
    
    double sum;
    MPI_Allreduce(&my_sum,&sum,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
    sum = sqrt(sum);

    for (int icv = 0; icv < ncv; ++icv) 
      mode_real[icv] /= sum;
   
    build_internal_operator(freq,sos);

    addBoundaryConditions(freq);

    // extract the real part of the matirx (imaginary part assumed zero below)..

    double * A_rl = new double[cvocv_i[ncv]];
    for (int ii = 0; ii < cvocv_i[ncv]; ++ii) { 

      assert( abs(A_op[ii].imag()) < 1.0e-15);
      A_rl[ii] = A_op[ii].real();

    }

    double * v = new double[ncv_g];
    double * w = new double[ncv_g];
    double * y = new double[ncv_g];

    updateCvData(mode_real);

    for (int icv = 0; icv < ncv_g; ++icv) { 
      v[icv] = mode_real[icv];
      y[icv] = mode_real[icv];
    }

    int done = 0;
    int iter = 0;

    const double zero    = getDoubleParam("SOLVER_ZERO", 1.0e-6);
    const int maxiter    = getIntParam("MAX_ITER", 10000);
    const int maxiter_lambda    = getIntParam("MAX_ITER_LAMBDA", 7);
    const double eig_tol = getDoubleParam("EIG_TOL", 1.0e-6); 
    const double df_range = getDoubleParam("DF_RANGE", 10.0);
    const bool b_verbose  = getBoolParam("VERBOSE_EIG_SOLVER", false);

    double lambda    = 0.0;

    while ( done == 0) { 

      ++iter;

      if ( mpi_rank == 0) 
        cout << " ****** starting iter : " << iter << endl;

      for (int icv = 0; icv < ncv_g; ++icv) 
        w[icv] = vol_cv[icv]*y[icv];

      solveCvCg(y,A_rl,w,zero,maxiter,b_verbose);

      double my_sum = 0.0;
      double sum;

      for (int icv = 0; icv < ncv; ++icv) 
        my_sum += y[icv]*vol_cv[icv]*y[icv];

      MPI_Allreduce(&my_sum,&sum,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
      sum = sqrt(sum);

      for (int icv = 0; icv < ncv; ++icv)
        y[icv] /= sum;

      updateCvData(y);

      my_sum = 0.0;
      for (int icv = 0; icv < ncv; ++icv) { 

        double Ay = 0.0;
        for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) { 
          const int icv_nbr = cvocv_v[coc];
          Ay += (A_op[coc].real())*y[icv_nbr];
        }

        my_sum += y[icv]*Ay;
      }

      MPI_Allreduce(&my_sum,&sum,1,MPI_DOUBLE,MPI_SUM,mpi_comm);

      //if ( mpi_rank == 0) 
      //  cout << " sum : " << sum << endl;

      const double lambda_new = sum;

      if ( abs(lambda_new) < 8.0*M_PI*M_PI*df_range*freq0) { 
      
        const double relax      = 1.0;
        
        const double lambda_old = lambda;
        lambda                  = (1.0-relax)*lambda + relax*lambda_new;
        
        if ( mpi_rank == 0) 
          cout << " >> iter, new perturb [delta f*4*pi^2], rel diff [prev iter] : " << iter << "   " 
            << lambda << "   " << abs((lambda-lambda_old)/lambda) << endl;
     
        if ( (abs(lambda-lambda_old) <= eig_tol*abs(lambda))) { 
          
          done = 1;
        
        } else { 

          // udpate the operator  ... 
          
          for (int icv = 0; icv < ncv; ++icv) { 
            
            for (int coc = cvocv_i[icv]; coc < cvocv_i[icv+1]; ++coc) 
              A_rl[coc] = A_op[coc].real();
            
            A_rl[cvocv_i[icv]] -= lambda*vol_cv[icv];
            
          }
          
          done = 0;

        }

      } else if ( iter >= maxiter_lambda) { 

        done = -1;

      } else { 

        if ( mpi_rank == 0) 
          cout << " >> iter, new perturb (no update): " << iter << "   "  << lambda << endl;

        done = 0;
      }

    }

    if ( done == 1) { 

      if ( mpi_rank == 0) 
        cout << " converged frequency = " << sqrt(4.0*M_PI*M_PI*freq*freq + lambda)/(2.0*M_PI) << endl;

    }

    for (int icv = 0; icv < ncv_g; ++icv) 
      mode_real[icv] = y[icv];

    delete[] w;
    delete[] v;
    delete[] y;
    delete[] A_rl;

  } 

  void solve_freq(const double freq) { 

    static bool first = true;

    if ( mpi_rank == 0) 
      cout << " > starting freq : " << freq << endl;

    const double zero = getDoubleParam("SOLVER_ZERO", 1.0e-06);

    build_internal_operator(freq,sos);
    
    addBoundaryConditions(freq);

    dcmplx * rhs = new dcmplx[ncv_g];

    if ( psi == NULL)  
      psi = new dcmplx[ncv_g];

    double x_seed[3];
    if ( Param * param = getParam("X_SEED")) { 
      for (int i = 0; i < 3; ++i)
        x_seed[i] = param->getDouble(i);
    } else { 
      CERR(" > must specify an x_seed ");
    }

    int icv_seed = findClosestCv(x_seed);
   
    for (int icv = 0; icv < ncv; ++icv) {

      if ( icv_seed == icv) {  
        rhs[icv] = dcmplx(1.0,0.0); 
      }
      else  { 
        rhs[icv] = dcmplx(0.0,0.0);
      }

    }

    if ( first) { 

      for (int icv = 0; icv <ncv; ++icv) { 
        if ( icv_seed == icv) 
          psi[icv] = dcmplx(1.0,0.0);
        else 
          psi[icv] = dcmplx(0.0,0.0);
      }
    }

    // now, let's try to get a measure of resonance by trying to solve 
    // this system.  the solutions will explode in norm if we find an 
    // eigenvalue .. 
    
    updateCvData_(psi);
    updateCvData_(rhs);

    const int maxiter = getIntParam("MAX_ITER", 100000);

    solve(psi,A_op,rhs,zero,maxiter,false);

    // report the magnitude of the solution .. 

    double my_sum[2] = {0.0,0.0};
    double sum[2];

    for (int icv = 0; icv < ncv; ++icv) {  
      my_sum[0] += vol_cv[icv]*(psi[icv].real()*psi[icv].real() + psi[icv].imag()*psi[icv].imag());
      my_sum[1] += vol_cv[icv];
    }

    MPI_Reduce(my_sum,sum,2,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

    if ( mpi_rank == 0) {
      cout << " XXXXX  " << freq << "   " << sqrt(sum[0]/sum[1]) << endl;
      sweepFile << freq << "\t" << sqrt(sum[0]/sum[1]) << endl;
    }


    delete[] rhs;

    first = false;
  }

  void updateCvData_(dcmplx* phi) { 

    // this may not be the most efficient .. 

    ComplexComm* cc = new ComplexComm();
    cc->setVar(phi);

    updateCvClassStart<ComplexComm,double>(cc);
    updateCvClassFinish<ComplexComm,double>(cc);

    delete cc;
  }

  int solve(dcmplx * phi,  const dcmplx * A, const dcmplx* rhs, 
            const double zero, const int maxiter, const bool verbose) { 

    // in the initial case, where the A_op is real, then we can split 
    // the real and imaginary parts of the sol.. 

    double * phi_rl = new double[ncv_g];
    double * phi_im = new double[ncv_g];
    double * rhs_rl = new double[ncv_g];
    double * rhs_im = new double[ncv_g];

    for (int icv = 0; icv < ncv_g; ++icv) { 

      phi_rl[icv] = phi[icv].real();
      phi_im[icv] = phi[icv].imag();
      rhs_rl[icv] = rhs[icv].real();
      rhs_im[icv] = rhs[icv].imag();

    }

    double * A_rl = new double[cvocv_i[ncv]];
    for (int ii = 0; ii < cvocv_i[ncv]; ++ii) { 

      assert( abs(A[ii].imag()) < 1.0e-15);
      A_rl[ii] = A[ii].real();

    }


    // solve via cg.. 

    int done0 = solveCvCg(phi_rl,A_rl,rhs_rl,zero,maxiter,verbose);
    int done1 = solveCvCg(phi_im,A_rl,rhs_im,zero,maxiter,verbose);

    updateCvData(phi_rl);
    updateCvData(phi_im);

    for (int icv = 0; icv < ncv; ++icv) { 

      phi[icv] = dcmplx(phi_rl[icv],phi_im[icv]);

    }

    delete[] A_rl;
    delete[] phi_rl;
    delete[] phi_im;
    delete[] rhs_rl;
    delete[] rhs_im;

    return ( (done0 == 1)&&(done1 == 1));

  }

  int solve_(dcmplx * phi,  const dcmplx * A, const dcmplx* rhs, 
            const double zero, const int maxiter, const bool verbose) { 

    const double relax = getDoubleParam("JACOBI_RELAX", 0.6);

    dcmplx *res      = new dcmplx[ncv];

    int iter = 0;
    int done = 0;
    while (done == 0) {

      iter++;

      // calculate the residual...
      FOR_ICV {
        const int coc_f = cvocv_i[icv];
        res[icv] = (rhs[icv] - A[coc_f]*phi[icv])/vol_cv[icv];
        for (int coc = coc_f+1; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          res[icv] -= A[coc]*phi[icv_nbr]/vol_cv[icv];
        }
        res[icv] /= A[coc_f]/vol_cv[icv];
      }

      // update the active u's...
      FOR_ICV phi[icv]   += relax*res[icv];

      // and ghosts...
      updateCvData_(phi);

      // check...
      double my_res_max = 0.0;
      FOR_ICV my_res_max = max(my_res_max,max(abs(res[icv].real()),abs(res[icv].imag())));
      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank == 0) {
        if ((verbose)||(iter > maxiter/2))
          cout << " > solveCvJacobi iter, res_max: " << iter << " " << res_max << endl;
        if (res_max < zero) {
          done = 1;
        }
        else if (iter > maxiter) {
          cout << " > Warning: solveCvJacobi did not converge after " << maxiter << " iters, res_max: " << res_max << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

    }

    delete[] res;

    return( done == 1 );
  } 

  int solve__(dcmplx * phi,
      const dcmplx * const A,
      const dcmplx * const rhs,
      const double zero,
      const int maxiter,
      const bool verbose) {
    
    assert(0);

    // assume we come in with a consistent initial condition...

    // we need the following work arrays...

    dcmplx* inv_diag   = new dcmplx[ncv_g];
    //for (int icv = 0; icv < ncv_g; ++icv) 
    //  inv_diag[icv] = dcmplx(1.0/vol_cv[icv]/vol_cv[icv],0.0);

    for (int icv = 0; icv < ncv_g; ++icv) 
      inv_diag[icv] = dcmplx(1.0,0.0);

    dcmplx* p         = new dcmplx[ncv_g];
    dcmplx* v         = new dcmplx[ncv_g]; // Apk
    dcmplx* res       = new dcmplx[ncv_g];
    dcmplx* tmp       = new dcmplx[ncv_g];

    int iter = 0; 
    int done = 0;
    while ( done == 0) { 

      ++iter;

      // compute residual .. 

      for (int icv = 0; icv < ncv; ++icv) {
        res[icv] = rhs[icv] - A[cvocv_i[icv]]*phi[icv];
        for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          res[icv] -= A[coc]*phi[icv_nbr];
        }

        tmp[icv]    = res[icv]*inv_diag[icv]; 
      }

      updateCvData_(res);
      updateCvData_(tmp);

      // check the residual function .. 

      double my_sum = 0.0;

      for (int icv = 0; icv < ncv; ++icv) { 
        my_sum += inv_diag[icv].real()*(res[icv].real()*res[icv].real() + 
                                        res[icv].imag()*res[icv].imag()  );
      }

      double sum;
      MPI_Allreduce(&my_sum,&sum,1,MPI_DOUBLE,MPI_SUM,mpi_comm);

      if ( verbose) { 
        cout << " iter, res : " << iter << "    " << sum << endl;
      }

      if ( sum < zero) { 
        done = 1;
      }

      // build pk .. 

      for (int icv = 0; icv < ncv; ++icv) { 

        p[icv] = 0.0;
        for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) { 
          const int icv_nbr = cvocv_v[coc];
          p[icv]           += A[coc]*tmp[icv_nbr]; // A= A^H
        }
      }

      updateCvData_(p);

      // now build v ... = A pk.. 

      for (int icv = 0; icv < ncv; ++icv) { 
        
        v[icv] = 0.0;

        for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) { 
          const int icv_nbr = cvocv_v[coc];
          v[icv]           += A[coc]*p[icv_nbr];
        }

      }

      updateCvData_(v);

      // now compute the dot products.. 
         
      double my_buf[2] = {0.0,0.0};

      for (int icv = 0; icv < ncv; ++icv) { 

        my_buf[0] += p[icv].real()*p[icv].real() + p[icv].imag()*p[icv].imag();
        my_buf[1] += inv_diag[icv].real()*(v[icv].real()*v[icv].real() + 
                                           v[icv].imag()*v[icv].imag()  );

      }

      double buf[2];
      MPI_Allreduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,mpi_comm);

      const double alpha = buf[0]/(buf[1] + 1.0e-20);

      if ( mpi_rank == 0 ) 
        cout << "   alpha : " << alpha << endl;

      for (int icv = 0; icv < ncv; ++icv) 
        phi[icv] += dcmplx(alpha)*p[icv];

      updateCvData_(phi);

      if ( iter == maxiter) 
        done = -1;
    }

    // XXX debug .. 

    for (int icv = 0; icv < ncv; ++icv) 
      cout << " YYYY  " << x_cv[icv][0] << "    " << res[icv].real() << "    " << res[icv].imag() << endl;


    delete[] tmp;
    delete[] v;
    delete[] p;
    delete[] res;
    delete[] inv_diag;

    return (done == 1);
  }


  int findClosestCv(const double x[3]) { 

    if ( cvAdt == NULL) buildCvAdt();

    double my_d2    = HUGE_VAL;
    int icv_closest = -1;

    vector<int> bbox_vec;
    assert( bbox_vec.empty());
    cvAdt->buildListForPoint(bbox_vec,x);

    for (int ii = 0, lim = bbox_vec.size(); ii < lim; ++ii) { 
     
      const int icv = bbox_vec[ii];
      const double d2 = DIST2(x_vv[icv],x);
      if ( (icv_closest == -1) || (d2 < my_d2)) { 
        icv_closest = icv;
        my_d2       = d2;
      }
    }

    double d2;
    MPI_Allreduce(&my_d2,&d2,1,MPI_DOUBLE,MPI_MIN,mpi_comm);

    int my_count = 0;
    
    if ( my_d2 == d2 )  
      my_count++;
    else 
      icv_closest = -1; // not yours.

    int count;
    MPI_Allreduce(&my_count,&count,1,MPI_INT,MPI_SUM,mpi_comm);

    if ( count != 1) { 
      CERR( " > oops .. need to implement tie break: " << count );
    }

    return icv_closest;
  }

  void addBoundaryConditions(const double freq) { 

    // for now, we presume that d\hat{p}/dn = 0 at all boundaries..
    // which leads to no modification to the operator... 
 
    FOR_PARAM_MATCHING("ADD_BC") { 

      const string zone_name = param->getString(0);
      const string bc_type   = param->getString(1); 

      if ( bc_type == "DIRICHLET") { 

        const double bc_val = param->getDouble(2);
        assert( bc_val == 0.0);

        map<const string,int>::const_iterator iter = bfZoneNameMap.find(zone_name);
        if ( iter == bfZoneNameMap.end()) { 

          CWARN( " > skipping bc for zone : " << zone_name << " ; zone not found " ); 
        
        } else { 

          const int iz = iter->second;
          assert( bfZoneVec[iz].getName() == zone_name);

          for (int ibf = bfZoneVec[iz].ibf_f; ibf <= bfZoneVec[iz].ibf_l; ++ibf) { 

            const int icv = cvobf[ibf];
            const double sos2 = sos[icv] * sos[icv];
            const dcmplx fax  = dcmplx(sos2*area_over_delta_bf[ibf]);


            const int coc00 = cvocv_i[icv];
            A_op[coc00]    += fax;

          } 

        }


      } else { 

        CERR( " > unknown bc type for acoustic solver : " << bc_type );

      } 


    } 

 
  }

  void build_internal_operator(const double freq, const double* sos) { 

    // the following builds a complex helmholtz operator ... 

    assert( cvocv_i != NULL);
    assert( cvocv_v != NULL);

    // i*omega
    const dcmplx freq_factor(-4.0*M_PI*M_PI*freq*freq,0.0);

    if ( A_op == NULL) 
      A_op = new dcmplx[cvocv_i[ncv]];

    for (int ii = 0; ii < cvocv_i[ncv]; ++ii) 
      A_op[ii] = 0.0;

    for (int icv = 0; icv < ncv; ++icv) 
      A_op[cvocv_i[icv]] = freq_factor*dcmplx(vol_cv[icv]);

    // use compact faces to discretize the laplacian 
  
    for (int ifa = 0; ifa < nfa; ++ifa) { 

      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];

      const double sos_avg = 0.5*(sos[icv0] + sos[icv1]);
      const double sos2    = sos_avg*sos_avg;
      const dcmplx fax     = dcmplx(sos2*area_over_delta_fa[ifa]);

      const int coc00   = cvocv_i[icv0];
      int coc01         = coc00;
      while ( coc01 != cvocv_i[icv0+1]) { 
        if ( cvocv_v[++coc01] == icv1) 
          break;
      }
      assert( coc01 != cvocv_i[icv0+1]);

      A_op[coc00]      += fax;
      A_op[coc01]      -= fax;

      if ( icv1 < ncv) { 

        const int coc11 = cvocv_i[icv1];
        int coc10       = coc11;
        while ( coc10 != cvocv_i[icv1+1]) { 
          if ( cvocv_v[++coc10] == icv0)
            break;
        }
        assert( coc10 != cvocv_i[icv1+1]);

        A_op[coc10]    -= fax;
        A_op[coc11]    += fax;

      }
    }
  }


};


int main(int argc, char* argv[]) { 

  try { 

    CTI_Init(argc,argv,"resonance.in");

    { 

      Resonance solver;
      solver.init();
      solver.run();

    }

    CTI_Finalize();

  }
  catch( int e) { 

    if ( e == 0) 
      CTI_Finalize();
    else 
      CTI_Abort();

  }
  catch(...) { 

    CTI_Abort();

  }

  return 0;
} 
