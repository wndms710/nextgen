
void vrem_with_comments() {


  for (int icv = 0; icv < ncv; ++icv) { 
    
    const double dx2 = pow( vol_cv[icv], 2.0/3.0 ); // XXX precompute..
    
    const double alpha11   = dudx[icv][0][0];
    const double alpha11_2 = alpha11*alpha11;
    const double alpha22   = dudx[icv][1][1];
    const double alpha22_2 = alpha22*alpha22;
    const double alpha33   = dudx[icv][2][2];
    const double alpha33_2 = alpha33*alpha33;
    
    const double alpha12   = dudx[icv][0][1];
    const double alpha12_2 = alpha12*alpha12;
    const double alpha13   = dudx[icv][0][2];
    const double alpha13_2 = alpha13*alpha13;
    const double alpha23   = dudx[icv][1][2];
    const double alpha23_2 = alpha23*alpha23;
    
    const double alpha21   = dudx[icv][1][0];
    const double alpha21_2 = alpha21*alpha21;
    const double alpha31   = dudx[icv][2][0];
    const double alpha31_2 = alpha31*alpha31;
    const double alpha32   = dudx[icv][2][1];
    const double alpha32_2 = alpha32*alpha32;

    const double beta11  = dx2*(alpha11_2+alpha21_2+alpha31_2);
    const double beta22  = dx2*(alpha12_2+alpha22_2+alpha32_2);
    const double beta33  = dx2*(alpha13_2+alpha23_2+alpha33_2);
    
    const double beta12  = dx2*(alpha11*alpha12+alpha21*alpha22+alpha31*alpha32);
    const double beta13  = dx2*(alpha11*alpha13+alpha21*alpha23+alpha31*alpha33);
    const double beta23  = dx2*(alpha12*alpha13+alpha22*alpha23+alpha32*alpha33);
    
    double B       = (beta11*beta22-beta12*beta12)+(beta11*beta33-beta13*beta13)+(beta22*beta33-beta23*beta23);
    B              = (B + abs(B))*0.5;
    
    const double alal    =
      ((alpha11_2+alpha22_2)+(alpha33_2+alpha12_2))+
      ((alpha13_2+alpha23_2)+(alpha21_2+alpha31_2))+alpha32_2;
    
    
    const double s_mag = sqrt(B/(alal+1.0E-20)); // includes lengthscale squared too
    
    // mu_sgs...
    mu_sgs[icv] = rho[icv]*vreman_coeff*s_mag;
    
    // clip at the constant smagorinsky level with c2 = 0.5^2...
    //mu_sgs[icv] = min(mu_sgs[icv],rho[icv]*0.05*dx2*sqrt(2.0*alal));
    
    // acoustic dissipation (to treat unresolved density fluctuations; see notes)
    const double sos2      = sos[icv]*sos[icv]; 
    //a_sgs[icv]             = a_sgs_coeff*dx2*abs(lap_smag[icv]);
    //a_sgs[icv]             = a_sgs_coeff*abs(lap_smag[icv]);
    
    a_sgs[icv]               = a_sgs_coeff*s_mag; // using inv time scale from the vreman kernel... 
    
    const double Masq      = (DOT_PRODUCT(u[icv],u[icv]))/(sos[icv]*sos[icv]);
    
    const double dil_      = dudx[icv][0][0] + dudx[icv][1][1] + dudx[icv][2][2];
    const double tmp_      = abs(dil_)*sqrt(dx2)/sos[icv];
    if ( tmp_ > 0.01) 
      a_sgs[icv]          += a_sgs_coeff2*sa_to_vol_ratio[icv]*Masq*sos[icv];
    
    a_sgs[icv]            /= sos2; // for use with the new flux.. 
    
  }
}  

#ifdef NEW_
  void setFgrCvs() { 

    const double gr_nrm = getDoubleParam("GR_NRM", 0.03);

    assert( vv2 != NULL);
    //double * vv2    = new double[ncv_g2];
    double * T_stag = new double[ncv_g2];  
    double * tmp      = new double[ncv_g2];

    const double cp  = R_gas*gamma/(gamma-1.0);
    const double gm1 = gamma - 1.0;
    
    for (int icv = 0; icv < ncv; ++icv)  { 

      vv2[icv]    = 0.0;
      tmp[icv]    = 0.0;
      T_stag[icv] = (1.0 + 0.5*gm1*(DOT_PRODUCT(u[icv],u[icv]))/(sos[icv]*sos[icv]))*T[icv];

    }

    updateCv2Data(T_stag);

    for (int ief = 0; ief < nef; ++ief) { 

      const int icv0 = cvoef[ief][0];
      const int icv1 = cvoef[ief][1];

      const double dp       = p[icv1] - p[icv0];
      const double dh       = cp*(T[icv1] - T[icv0]);
      const double v_avg    = 0.5*(1.0/rho[icv0] + 1.0/rho[icv1]);
      const double invT_avg = 0.5*(1.0/(T[icv0]) + 1.0/(T[icv1]));
      const double ds       = R_gas*(ent[icv1] - ent[icv0]);  // the entropy already has an R factor built in...

      //const double gr       = -ds + invT_avg*(dh - v_avg*dp);
      const double gr       = R_gas*log(p[icv1]/p[icv0]) - v_avg*dp*invT_avg;

      double u_avg[3];
      for (int i = 0; i < 3; ++i) 
        u_avg[i] = 0.5*(u[icv0][i] + u[icv1][i]);


      const double flux     = gr*DOT_PRODUCT(u_avg,n_ef[ief])/v_avg;
      tmp[icv0]            += flux;
      if ( icv1 < ncv) 
        tmp[icv1]          -= flux;

    }

    for (int icv = 0; icv < ncv; ++icv) { 
      
      // normalize .. 
      const double delta  = pow(vol_cv[icv], 1.0/3.0);
      tmp[icv]           *= delta/vol_cv[icv]/R_gas/rho[icv]/sos[icv];
      tmp[icv]            = abs(tmp[icv]);
    }

    updateCv2Data(tmp);


    for (int ief = 0; ief < nef; ++ief) { 

      const int icv0 = cvoef[ief][0];
      const int icv1 = cvoef[ief][1];
      const double tt = max(tmp[icv0],tmp[icv1]);

      vv2[icv0] = max(tt,vv2[icv0]);
      if ( icv1 < ncv) 
        vv2[icv1] = max(tt,vv2[icv1]);
    }

    updateCv2Data(vv2);

    MiscUtils::dumpRange(vv2,ncv_g2, "vv2");

    for (int icv = 0; icv < ncv_g2; ++icv)  
      cv_compact[icv].fgr = gr_func(vv2[icv],gr_nrm);

    if ( checkParam("NO_DISS") ) { 
      for (int icv = 0; icv < ncv_g2; ++icv) 
        cv_compact[icv].fgr = 0.0;
    }

    delete[] T_stag;
    //delete[] vv2;

  } 
#endif

void ducros_stuff() { 

  double * tmp = new double[ncv_g2];

    for (int icv = 0; icv < ncv; ++icv) { 

      // the code seems stable but is not activating dissipation
      // near shocks.  we are activating a modified ducros senor 
      // to try to provide some level of dissipation near the shocks.
      
      const double dil = dudx[icv][0][0] + dudx[icv][1][1] + dudx[icv][2][2];
      double vort[3];
      vort[0]  = dudx[icv][2][1] - dudx[icv][1][2];
      vort[1]  = dudx[icv][0][2] - dudx[icv][2][0];
      vort[2]  = dudx[icv][1][0] - dudx[icv][0][1];

      const double eps   = 1.0e-10;
      const double v_mag = MAG(vort);
      double ds          = dil*dil/(dil*dil + v_mag*v_mag + eps);
      ds                *= abs(dil)*pow(vol_cv[icv],1.0/3.0)/sos[icv];
      tmp[icv]           = ds;

    }

    MiscUtils::dumpRange(tmp,ncv,"ds_tmp");

    updateCv2Data(tmp);

    const double scal2 = getDoubleParam("DS_COEFF", 0.5);

    for (int icv = 0; icv < ncv_g2; ++icv) { 
      //cv_compact[icv].fgr += 2.0*tmp[icv]*tmp[icv]; 
      //cv_compact[icv].fgr  += scal2*tmp[icv]; 
      cv_compact[icv].fgr += gr_func(scal2*tmp[icv],gr_nrm);
      cv_compact[icv].fgr  = min(cv_compact[icv].fgr,0.5);
    }

    delete[] tmp;
}

    // conservation check...
    if (cht) {
      double my_buf[2] = { 0.0, 0.0 };
      FOR_ICV my_buf[0] += rhoE[icv]*vol_cv[icv];
      my_buf[1] = cht->calcEnergyLocal();
      double buf[2];
      MPI_Reduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      if (mpi_rank == 0) cout << "ENERGY_BALANCE:   rk3 " << buf[0] << " " << buf[1] << " " << buf[0]+buf[1] << " change: " << buf[0]+buf[1]-esum << endl;
    }



      double q_sens = 0.0;
      double q_rhs = 0.0;
      double T_min = HUGE_VAL;
      double T_max = -HUGE_VAL;
      FOR_INO {
        q_sens += rho*cp*T_no[ino]*vol_no[ino];
        q_rhs += q_no[ino];
        T_min = min(T_min,T_no[ino]);
        T_max = max(T_max,T_no[ino]);
      }
      cout << "q_sens at start: " << q_sens << " q_rhs: " << q_rhs << " T_min, T_max: " << T_min << " " << T_max << endl;
      
      cout << "about to cg solve" << endl;



      cout << "back from gg solve" << endl;
      
      MiscUtils::solveGaussSidelSerial(T_no,nno,A,rhs,noono_i,noono_v,1.0E-12,1000,0.5,true);

      double q_sens1 = 0.0;
      T_min = HUGE_VAL;
      T_max = -HUGE_VAL;
      FOR_INO {
        q_sens1 += rho*cp*T_no[ino]*vol_no[ino];
        T_min = min(T_min,T_no[ino]);
        T_max = max(T_max,T_no[ino]);
      }
      cout << "q_sens1 at end: " << q_sens1 << " T_min, T_max: " << T_min << " " << T_max << " check: " << q_sens1-q_sens - q_rhs*dt << endl;

      
      





// fast fake build of sspobf_i/v/wgt...
      

      // build a fake sspobf_i/v...
      vector<int> sspobf_v_vec;
      vector<double> sspobf_wgt_vec;
      assert(sspobf_i == NULL);
      sspobf_i = new int[nbf+1];
      sspobf_i[0] = 0;
      set<int> sspSet;
      for (int ibf = 0; ibf < nbf; ++ibf) {
        assert(sspSet.empty());
        for (int sob = sstobf_i[ibf]; sob != sstobf_i[ibf+1]; ++sob) {
          const int isst = sstobf_v[sob];
          //int ist,bits;
          //subSurface->getIstGlobalAndBits(ist,bits,isst);
          //assert(bits == 0); // no transforms for now
          FOR_I3 {
            const int issp = subSurface->spost[isst][i];
            sspSet.insert(issp);
          }
        }
        // sspSet contains the subsurface points associated with this bf. Assume
        // a uniform weighting for now...
        assert(!sspSet.empty());
        const double wgt = 1.0/double(sspSet.size());
        for (set<int>::iterator iter = sspSet.begin(); iter != sspSet.end(); ++iter) {
          sspobf_v_vec.push_back(*iter);
          sspobf_wgt_vec.push_back(wgt);
        }
        sspSet.clear();
        sspobf_i[ibf+1] = sspobf_v_vec.size();
      }
      // copy the local vec's into the class sspobf structure...
      assert(sspobf_v == NULL);
      sspobf_v = new int[sspobf_i[nbf]];
      assert(sspobf_wgt == NULL);
      sspobf_wgt = new double[sspobf_i[nbf]];
      for (int sob = 0; sob < sspobf_i[nbf]; ++sob) {
        sspobf_v[sob] = sspobf_v_vec[sob];
        sspobf_wgt[sob] = sspobf_wgt_vec[sob];
      }
      sspobf_v_vec.clear();
      sspobf_wgt_vec.clear();

      






void dumpGhostCounts() const { 
    for (int rank = 0; rank < mpi_size; ++rank) { 
      if ( mpi_rank == rank) { 
        cout << " mpi_rank, cv ghost size: " << mpi_rank << "   " << ncv << "   " << ncv_g2-ncv 
          << "   " << double(ncv_g2-ncv)/double(ncv) << endl;
        cout.flush();
      }
      MPI_Barrier(mpi_comm);
    }
  }

void calcCvGrad2(double (*__restrict__ dudx)[3][3], double (*__restrict__ dpdx)[3], 
                  const double (*__restrict__ u)[3], const double *__restrict__ p) { 
    
    for (int icv = 0; icv < ncv; ++icv) { 
      for (int i =0; i < 3; ++i) 
        dpdx[icv][i] = 0.0;
      for (int i =0; i < 3; ++i) 
        for (int j = 0; j < 3; ++j) 
          dudx[icv][i][j] = 0.0;
    }

    for (int ifa = 0; ifa < nfa_i; ++ifa) { 
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      
      const double p_avg = 0.5*(p[icv0] + p[icv1]);
      double u_avg[3];
      for (int i =0; i < 3; ++i) 
        u_avg[i] = 0.5*(u[icv0][i] + u[icv1][i]);
     
      for (int i =0; i < 3; ++i) { 
        const double flux = n_fa[ifa][i]*p_avg;
        dpdx[icv0][i]    += flux;
        dpdx[icv1][i]    -= flux;
      }

      for (int i =0; i < 3; ++i) { 
        for (int j =0; j < 3; ++j) { 
          const double flux = n_fa[ifa][j]*u_avg[i];
          dudx[icv0][i][j] += flux;
          dudx[icv1][i][j] -= flux;
        }
      }
    }

    for (int ifa = nfa_i; ifa < nfa; ++ifa) { 
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      
      const double p_avg = 0.5*(p[icv0] + p[icv1]);
      double u_avg[3];
      for (int i =0; i < 3; ++i) 
        u_avg[i] = 0.5*(u[icv0][i] + u[icv1][i]);
     
      for (int i =0; i < 3; ++i) { 
        const double flux = n_fa[ifa][i]*p_avg;
        dpdx[icv0][i]    += flux;
      }

      for (int i =0; i < 3; ++i) { 
        for (int j =0; j < 3; ++j) { 
          const double flux = n_fa[ifa][j]*u_avg[i];
          dudx[icv0][i][j] += flux;
        }
      }
    }

    // neumann closure for the gradient at the boundary...
    for (int ibf = 0; ibf < nbf; ++ibf) { 
      const int icv = cvobf[ibf];
      for (int i =0; i < 3; ++i) { 
        dpdx[icv][i] += n_bf[ibf][i]*p[icv];
      }

      for (int i =0; i < 3; ++i) { 
        for (int j =0; j < 3; ++j) { 
          dudx[icv][i][j] += n_bf[ibf][j]*u[icv][i];
        }
      }
    }

    for (int icv =0; icv < ncv; ++icv) { 
      
      for (int i =0; i < 3; ++i) 
        dpdx[icv][i] *= inv_vol[icv];

      for (int i = 0; i < 3; ++i) 
        for (int j = 0; j <3 ; ++j) 
          dudx[icv][i][j] *= inv_vol[icv];
    }
  }
  void addCrappyInternalRhs(IdealGasRhs* rhs, const double time, const int rk_stage) {

    // u, rho, rhoE, p, T 

    const double gm1    = gamma - 1.0;
    const double invgm1 = 1.0/gm1;
    const double Rgogm1 = R_gas*gamma*invgm1;
    
    FOR_IEF {
      
      const int icv0 = cvoef[ief][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvoef[ief][1]; assert((icv1 >= 0)&&(icv1 < ncv_g2));
      
      // normal-flux...
      IdealGasRhs n_flux[2];
      IdealGasRhs c_flux[2];

      FOR_I2 {
	const int icv = cvoef[ief][i];
	// n...
	const double h = Rgogm1*T[icv];
	n_flux[i].rho            = rho[icv]*DOT_PRODUCT(u[icv],n_ef[ief]);
	FOR_J3 n_flux[i].rhou[j] = n_flux[i].rho*u[icv][j] + p[icv]*n_ef[ief][j];
	n_flux[i].rhoE           = n_flux[i].rho*(h + 0.5*DOT_PRODUCT(u[icv],u[icv]));
	// c flux...
	c_flux[i].rho            = rho[icv]*DOT_PRODUCT(u[icv],c_ef[ief]);
	FOR_J3 c_flux[i].rhou[j] = c_flux[i].rho*u[icv][j] + p[icv]*c_ef[ief][j];
	c_flux[i].rhoE           = c_flux[i].rho*(h + 0.5*DOT_PRODUCT(u[icv],u[icv]));
      }
      
      //double flux[3];	FOR_I3 flux[i] = 0.5*n_ef[ief][i]*(phi[icv0]+phi[icv1]) + 0.5*c_ef[ief][i]*(phi[icv1]-phi[icv0]);
      
      rhs[icv0].rho            -= 0.5*(n_flux[0].rho     + n_flux[1].rho)     + 0.5*(c_flux[1].rho     - c_flux[0].rho);
      FOR_I3 rhs[icv0].rhou[i] -= 0.5*(n_flux[0].rhou[i] + n_flux[1].rhou[i]) + 0.5*(c_flux[1].rhou[i] - c_flux[0].rhou[i]);
      rhs[icv0].rhoE           -= 0.5*(n_flux[0].rhoE    + n_flux[1].rhoE)    + 0.5*(c_flux[1].rhoE    - c_flux[0].rhoE);
      
      if (icv1 < ncv) {
	rhs[icv1].rho            += 0.5*(n_flux[0].rho     + n_flux[1].rho)     + 0.5*(c_flux[1].rho     - c_flux[0].rho);
	FOR_I3 rhs[icv1].rhou[i] += 0.5*(n_flux[0].rhou[i] + n_flux[1].rhou[i]) + 0.5*(c_flux[1].rhou[i] - c_flux[0].rhou[i]);
	rhs[icv1].rhoE           += 0.5*(n_flux[0].rhoE    + n_flux[1].rhoE)    + 0.5*(c_flux[1].rhoE    - c_flux[0].rhoE);
      }
    }

  }
class BasicSolver : public VgpWithTools { 

private:

  KillfileReader kfr;

public: 

  vector<Stats*> stats_vec;

  int step; 
  double time, dt;
  ImageManager * imageManager;

  BasicSolver() : VgpWithTools(), step(0), time(0.0), imageManager(NULL) {
    dt = 1.0;  // stats should behave as if it was an ensemble
    stats_vec.clear();
    //TODO move to VgpWithTools?
    b_boundingBox = false;
  }

  ~BasicSolver() { 
    if (imageManager  != NULL) delete imageManager;
    FOR_STATS delete *it;
  }

  void init() { 
    readRestart(this); 
    registerData(step,"step");
    registerData(time,"time"); 
    registerFromParams(); 

    initStats();
    imageManager = new ImageManager(this);

    if (checkParam("REPORT_ZONE_AREAS"))
      reportZoneAreas();

    
    //setup flags for imaging
    assert(bf_flag.isNull());
    bf_flag.setLength(nbf);

    assert(zone_flag.isNull());
    zone_flag.setLength(zone_map.size());

    for (int ibf = 0; ibf < nbf; ++ibf) {
      bf_flag[ibf] = znobf[ibf];
    }


  }

  void run() { 
    const string filename = getStringParam("RESTART", "restart.les");
    readData(static_cast<Vgp*>(this),filename);
    checkDataInitialization();
    dt = 1.0; // dt tracks an snapshot ensemble count 
    updateConservativeAndPrimitiveData();

    cout << "Hold mode on..." << endl;
    if (mpi_rank==0) logger->insertSolverAction(SIMULATION_RUN);
    kfr.setHold(true);
    while (Param * param = kfr.getParam("killsurfer")) {
      if (mpi_rank == 0) cout << "just got param \"" << param->str() << "\"" << endl;

      if (param->getName() == "WRITE_IMAGE") {
        processWriteImage(param);
      }
      else if (param->getName() == "WRITE_JSON") {
        processWriteJson(param);
      }
      else {
        if (mpi_rank == 0) cout << "Warning: skipping unrecognized param \"" << param->getName() << "\"" << endl;
      }
      cout << "done" << endl;

    }

  } 

  void checkDataInitialization() const {}
  void updateConservativeAndPrimitiveData() { syncData(); }  

  void initStats() { _initStatsCommon(stats_vec,this); }
  void resetStats() { FOR_STATS (*it)->reset(); }
  void updateStats(const cti_real _dt) { 
    FOR_STATS (*it)->update(step,_dt,true);
  }

  void syncData() { 
    for (vector<DataContainer*>::iterator it = vio_vec.begin(); it != vio_vec.end(); ++it) { 
      if (CvDataContainer* dc = dynamic_cast<CvDataContainer*>(*it)) {
        if (dc->name.find("stats-") == string::npos) { // stats variables don't have ghost vars, so not mp sync
          for (list<CtiRealScalar>::iterator sit = dc->sdata_vec.begin(); sit != dc->sdata_vec.end(); ++sit)  
            cv_comm->updateCvData(sit->r1);
	
          for (list<CtiRealVector>::iterator vit = dc->vdata_vec.begin(); vit != dc->vdata_vec.end(); ++vit) 
            cv_comm->updateCvData(vit->r2);
        }
      }
    }
  }

  int getCvR1(cti_real * r1, const string& name, const bool& doEval=false) { 
    return _getCvR1(r1,name);
  }
  int getCvR2(cti_real (*r2)[3], const string& name, const bool& doEval=false) { 
    return _getCvR2(r2,name);
  }
  int setCvR1(cti_real* r1, const string& name) { 
    return _setCvR1(r1,name);
  }
  int setCvR2(cti_real (*r2)[3], const string& name) { 
    return _setCvR2(r2,name);
  }
  
  virtual void temporalHook() {}
  virtual void finalHook() {}

  void runKillFileParams(KillFileManager * _kfm){
  }
 
  void processWriteJson(Param * param) {

    assert(param->getString(0) == "NAME");
    const string prefix = param->getString(1);

    cout << "JUST GOT WRITE_JSON command: " << param->getString() << endl;

    // write the JSON file...

    FILE * fp = fopen((prefix+".json").c_str(),"w");
    assert(fp);

    double xbb[3];
    getBoundingBoxCenter(xbb);
    fprintf(fp,"{\n\"BoundingBoxCentroid\":[%f,%f,%f]",xbb[0],xbb[1],xbb[2]);

    double diag = 2.0*getBoundingBoxRmax();
    fprintf(fp,",\n\"BoundingBoxDiagonal\":%f",diag);

    fprintf(fp,",\n\"LightLocation\":[0.039503,0,1.01189,-0.530497,-0.57,0.581893]");

    fprintf(fp,",\n\"webUIOutput\":[]");

    fprintf(fp,",\n\"zoneAreas\":[\n");
    for (int i = 0; i < zone_map.size(); ++i) {
      if (i == 0) fprintf(fp," %f",1.1234+i);
      else fprintf(fp,",\n %f",1.1234+i);
    }
    fprintf(fp,"\n]");

    fprintf(fp,",\n\"zoneCentroids\":[\n");
    for (int i = 0; i < zone_map.size(); ++i) {
      if (i == 0) fprintf(fp," [%f,%f,%f]",1.1234+i,2.1,345.1+i);
      else fprintf(fp,",\n [%f,%f,%f]",1.1234+i,2.1,345.1+i);
    }
    fprintf(fp,"\n]");

    fprintf(fp,",\n\"zoneIds\":[\n");
    int i=0;
    for(map<int,ZoneDescriptor>::const_iterator it = zone_map.begin(); it != zone_map.end(); ++it) {
      if (i == 0) fprintf(fp," %d",it->second.index);
      else fprintf(fp,",\n %d",it->second.index);
      ++i;
    }
    fprintf(fp,"\n]");

    fprintf(fp,",\n\"zoneNames\":[\n");
    i=0;
    for(map<int,ZoneDescriptor>::const_iterator it = zone_map.begin(); it != zone_map.end(); ++it) {
      if (i == 0) fprintf(fp," \"%s\"",it->second.name.c_str());
      else fprintf(fp,",\n \"%s\"",it->second.name.c_str());
      ++i;
    }
    fprintf(fp,"\n]");

    fprintf(fp,",\n\"zonePeriodicPairs\":[]");

    fprintf(fp,"\n}\n");

    fclose(fp);

    //rename(("."+prefix+".json").c_str(),(prefix+".json").c_str());

  }

  void processWriteImage(Param * param) {

    // note - any INTERVAL param in the WRITE_IMAGE command will be ignored here...

    CtiScene scene(param);
    if (!scene.skip(step)){

      scene.setRmax(getBoundingBoxRmax());
      double center[3];
      getBoundingBoxCenter(center);
      scene.setCenter(center[0],center[1],center[2]);
       
      //TODO explicitly add zone id (second.index, may not be sequential?)
      for(map<int,ZoneDescriptor>::const_iterator it = zone_map.begin(); it != zone_map.end(); ++it) {
        scene.addZoneName(it->second.name);
      }
 
      scene.initCanvas();

      //Approach #1 get the name of the data requested from the scene
      string surfaceVarStr;
      cti_real* no_data = NULL;
      if (scene.getSurfaceVar(surfaceVarStr)){
        cti_real* cv_data = new cti_real[ncv];
        if (getCvR1(cv_data,surfaceVarStr)>=0){
          no_data = new cti_real[nno];
          interpCvToNo(no_data,cv_data,this);
        }
        delete[] cv_data;
      }
      //if the scene's pointer is not null it will write data
      //what should be done if the user requested bad data?
      scene.setSurfaceData(no_data); 
      
      // Approach #2 somehow pass the data access function to the scene (or the full solver...)
      //scene.addData(this->BasicSolver::getCvR1);
      
      scene.addFaces(x_no, x_bf, noobf_i, noobf_v, znobf, nbf, bf_flag, zone_flag);
      //scene.addEdges(ss.xsp, ss.spost, ss.teost, ss.nst);
      scene.writeImage();

      if (!no_data){
        delete[] no_data;
        no_data = NULL;
      }
    }
  }


  private:

  //TODO where should these flags go...
  IntFlag bf_flag;
  IntFlag zone_flag;

  
  //methods that could be moved to VgpWithTools.hpp
  bool b_boundingBox; 
  double boundingBox[6];

  void calcBoundingBox() {
    assert(!b_boundingBox);

    double my_buf[6] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };
  
    FOR_INO{
      my_buf[0] = min(my_buf[0],x_no[ino][0]);
      my_buf[1] = min(my_buf[1],-x_no[ino][0]);
      my_buf[2] = min(my_buf[2],x_no[ino][1]);
      my_buf[3] = min(my_buf[3],-x_no[ino][1]);
      my_buf[4] = min(my_buf[4],x_no[ino][2]);
      my_buf[5] = min(my_buf[5],-x_no[ino][2]);
    }
    my_buf[1]*=-1.0;
    my_buf[3]*=-1.0;
    my_buf[5]*=-1.0;
    MPI_Allreduce(my_buf, boundingBox, 6, MPI_DOUBLE, MPI_SUM,mpi_comm);
    

    cout << "calcBoundingBox(): X " << boundingBox[0] << ":" << boundingBox[1] <<
                             ", Y " << boundingBox[2] << ":" << boundingBox[3] <<
                             ", Z " << boundingBox[4] << ":" << boundingBox[5] << endl;
    b_boundingBox = true;
  }
  
  void getBoundingBox(double _bBox[6]) {
    if (!b_boundingBox) calcBoundingBox();
    assert(b_boundingBox);
    _bBox[0] = boundingBox[0];
    _bBox[1] = boundingBox[1];
    _bBox[2] = boundingBox[2];
    _bBox[3] = boundingBox[3];
    _bBox[4] = boundingBox[4];
    _bBox[5] = boundingBox[5];
  }
  
  double getBoundingBoxRmax() {
    if (!b_boundingBox) calcBoundingBox();
    assert(b_boundingBox);
  
    return 0.5*sqrt((boundingBox[0]-boundingBox[1])*(boundingBox[0]-boundingBox[1])+
                    (boundingBox[2]-boundingBox[3])*(boundingBox[2]-boundingBox[3])+
                    (boundingBox[4]-boundingBox[5])*(boundingBox[4]-boundingBox[5]));
  }
  
  void getBoundingBoxCenter(double _bBoxCenter[3]) {
    if (!b_boundingBox) calcBoundingBox();
    assert(b_boundingBox);
    _bBoxCenter[0] = 0.5*(boundingBox[0]+boundingBox[1]);
    _bBoxCenter[1] = 0.5*(boundingBox[2]+boundingBox[3]);
    _bBoxCenter[2] = 0.5*(boundingBox[4]+boundingBox[5]);
  }




};
if ( Param * param = getParam("IGNITE")) { 

      // populate the current p, T..
      updateConservativeAndPrimitiveData(); 

      // XXX revisit the error parsing and forms of ignition..
      // IGNITE GEOM SPHERE POINT <x0> <x1> <x2> RADIUS <r0>
      assert( param->getString(0) == "GEOM");
      assert( param->getString(1) == "SPHERE");
      assert( param->getString(2) == "POINT");

      double x_center[3];
      x_center[0] = param->getDouble(3);
      x_center[1] = param->getDouble(4);
      x_center[2] = param->getDouble(5);

      assert( param->getString(6) == "RADIUS");
      const double radius = param->getDouble(7);

      // lookup Cmax everywhere... 
      double * Cmax_cv = new double[ncv];
      chemtable->lookupReduced(Cmax_cv,"upper",Z,ncv);
      chemtable->lookup(T_p,"T",e,"e",Z,Cmax_cv,ncv);

      for (int icv = 0; icv < ncv; ++icv) { 
        const double dist = DIST(x_cv[icv],x_center);
        if ( dist <= radius) {
          rho[icv]          = p[icv]/(R[icv]*T_p[icv]);
          C[icv]            = Cmax_cv[icv]; // advance the reaction..
          rhoE[icv]         = rho[icv]*(e[icv] + 0.5*DOT_PRODUCT(u[icv],u[icv]));
        }
      }

      delete[] Cmax_cv;
      updateConservativeAndPrimitiveData();
    }

 // non-compact viscous operators... 

void addViscousFlux(IdealGasRhs& flux, const FaNC& fa, const double x_vv0[3], const double x_vv1[3], 
                      const IdealGasState& cv0, const IdealGasState& cv1, 
                      const AuxGasProp& cv0_compact, const AuxGasProp& cv1_compact) { 

    const double mu_total_c             =   0.5*(cv0_compact.mu_total+cv1_compact.mu_total);
    const double loc_coeff              =   0.5*(cv0_compact.loc_total + cv1_compact.loc_total);

    double duds[3]; 
    double ds[3] = DIFF(x_vv1,x_vv0);
    double s_mag = 0.0;
    for (int i =0; i < 3; ++i) 
      s_mag += ds[i]*ds[i];

    s_mag = sqrt(s_mag);
    // ds now houses a unit vector
    for (int i =0; i < 3; ++i) 
      ds[i] /= s_mag;

    for (int i =0; i < 3; ++i) 
      duds[i] = (cv1.u[i] - cv0.u[i])/s_mag;

    const double dhds = (cv1.h - cv0.h)/s_mag;

    const double skk = duds[0]*ds[0] + duds[1]*ds[1] + duds[2]*ds[2];

    double tauij[3][3];
    for (int i =0; i < 3; ++i) { 
      
      for (int j =0; j < 3; ++j) { 
        tauij[i][j] = mu_total_c*(duds[i]*ds[j] + duds[j]*ds[i]);
      }
    
      tauij[i][i] -= 2.0/3.0*mu_total_c*skk;

    }

    double tauijnj[3];
    for (int i =0; i < 3; ++i) 
      tauijnj[i] = DOT_PRODUCT(tauij[i],fa.n);

    for (int i = 0; i < 3; ++i) { 
      flux.rhou[i] -= tauijnj[i];
    }


    double u_avg[3]; 
    for (int i =0; i < 3; ++i) 
      u_avg[i] = 0.5*(cv0.u[i] + cv1.u[i]);

    flux.rhoE -= loc_coeff*dhds*DOT_PRODUCT(ds,fa.n);
    for (int i =0; i < 3; ++i)
      flux.rhoE -= u_avg[i]*tauijnj[i];
  }

  void calcRhs(IdealGasRhs* rhs, const double time, const int rk_stage) {

    double (*x_vv_tmp)[3] = new double[ncv_g2][3];
    for (int icv = 0;icv < ncv; ++icv) 
      for (int i =0; i < 3; ++i) 
        x_vv_tmp[icv][i] = x_vv[icv][i];

    updateCv2Data(x_vv_tmp);

    for(vector<IdealGasBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)  
      (*it)->addBoundaryFlux(rhs);
    
    for (int ief = 0; ief < n_e__i; ++ief) { 
      const int icv0           = fa_hashed[ief].cvofa[0];
      const int icv1           = fa_hashed[ief].cvofa[1];
      const int ief_compressed = fa_hashed[ief].ief_compressed;
     
      cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
      cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],1,0);
      cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
      cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);
      
      IdealGasRhs flux;
      double u_avg[3], inv_v_avg, p_avg, h_stag;
      calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
                          c_diss,u_avg,inv_v_avg,h_stag,p_avg);

      addViscousFlux(flux,ef_geom[ief_compressed],x_vv_tmp[icv0],x_vv_tmp[icv1],cv_light[icv0],cv_light[icv1],
                     cv_compact[icv0],cv_compact[icv1]);

      rhs[icv0].rho -= flux.rho;
      rhs[icv1].rho += flux.rho;
      for (int i =0; i <3 ; ++i) { 
        rhs[icv0].rhou[i] -= flux.rhou[i];
        rhs[icv1].rhou[i] += flux.rhou[i];
      }
      rhs[icv0].rhoE -= flux.rhoE;
      rhs[icv1].rhoE += flux.rhoE;
    }

    // these faces have compact contributions to the viscous & modeling terms...
    for (int ief = n_e__i; ief < n_c__i; ++ief) { 
      const int icv0           = fa_hashed[ief].cvofa[0];
      const int icv1           = fa_hashed[ief].cvofa[1];
      const int ief_compressed = fa_hashed[ief].ief_compressed;
      const int icf_compressed = fa_hashed[ief].icf_compressed;
      
      cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
      cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],1,0);
      cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
      cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);
      
      IdealGasRhs flux;
      double u_avg[3], inv_v_avg, p_avg, h_stag;
      calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
                          c_diss,u_avg,inv_v_avg,h_stag,p_avg);

      cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
      cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);
 
      addViscousFlux(flux,ef_geom[ief_compressed],x_vv_tmp[icv0],x_vv_tmp[icv1],cv_light[icv0],cv_light[icv1],
                     cv_compact[icv0],cv_compact[icv1]);


      // aggregate into the flux...
      calcInternalCompactFluxNew(flux,cf_geom[icf_compressed],cv_light[icv0],cv_light[icv1],
                                 cv_compact[icv0],cv_compact[icv1],u_avg,inv_v_avg,h_stag,c_diss);

      rhs[icv0].rho -= flux.rho;
      rhs[icv1].rho += flux.rho;
      for (int i =0; i <3 ; ++i) { 
        rhs[icv0].rhou[i] -= flux.rhou[i];
        rhs[icv1].rhou[i] += flux.rhou[i];
      }
      rhs[icv0].rhoE -= flux.rhoE;
      rhs[icv1].rhoE += flux.rhoE;
    }

    // faces have a c component but are purely extended (associated with the convection)..
    for (int ief = n_c__i; ief < n_ec_i; ++ief) { 
      const int icv0           = fa_hashed[ief].cvofa[0];
      const int icv1           = fa_hashed[ief].cvofa[1];
      const int ief_compressed = fa_hashed[ief].ief_compressed;

      cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
      cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],1,0);
      cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
      cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);
      
      IdealGasRhs flux;
      double u_avg[3], inv_v_avg, p_avg, h_stag;
      calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
                          c_diss,u_avg,inv_v_avg,h_stag,p_avg);
  
      addViscousFlux(flux,ef_geom[ief_compressed],x_vv_tmp[icv0],x_vv_tmp[icv1],cv_light[icv0],cv_light[icv1],
                     cv_compact[icv0],cv_compact[icv1]);


      IdealGasRhs c_flux;
      calcInternalFluxSymmetricNew(c_flux,ef_geom[ief_compressed],u_avg,inv_v_avg,
                                   h_stag,p_avg,cv_light[icv0],cv_light[icv1]);

      rhs[icv0].rho -= flux.rho + c_flux.rho;
      for (int i =0; i < 3; ++i)
        rhs[icv0].rhou[i] -= flux.rhou[i] + c_flux.rhou[i];
      rhs[icv0].rhoE -= flux.rhoE + c_flux.rhoE;

      rhs[icv1].rho += flux.rho + c_flux.rho;
      for (int i =0; i < 3; ++i)
        rhs[icv1].rhou[i] += flux.rhou[i] + c_flux.rhou[i];
      rhs[icv1].rhoE += flux.rhoE + c_flux.rhoE;
    }

    // have c, and are also compact faces ... 
    for (int ief = n_ec_i; ief < n_cc_i; ++ief) { 
      const int icv0           = fa_hashed[ief].cvofa[0];
      const int icv1           = fa_hashed[ief].cvofa[1];
      const int ief_compressed = fa_hashed[ief].ief_compressed;
      const int icf_compressed = fa_hashed[ief].icf_compressed;

      cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
      cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],1,0);
      cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
      cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);
      
      IdealGasRhs flux;
      double u_avg[3], inv_v_avg, p_avg, h_stag;
      calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
                          c_diss,u_avg,inv_v_avg,h_stag,p_avg);
      
      addViscousFlux(flux,ef_geom[ief_compressed],x_vv_tmp[icv0],x_vv_tmp[icv1],cv_light[icv0],cv_light[icv1],
                     cv_compact[icv0],cv_compact[icv1]);


      IdealGasRhs c_flux;
      calcInternalFluxSymmetricNew(c_flux,ef_geom[ief_compressed],u_avg,inv_v_avg,
                                   h_stag,p_avg,cv_light[icv0],cv_light[icv1]); 

      cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
      cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);
 
      // aggregate into the flux...
      calcInternalCompactFluxNew(flux,cf_geom[icf_compressed],cv_light[icv0],cv_light[icv1],
                                 cv_compact[icv0],cv_compact[icv1],u_avg,inv_v_avg,h_stag,c_diss);

      rhs[icv0].rho -= flux.rho + c_flux.rho;
      for (int i =0; i < 3; ++i)
        rhs[icv0].rhou[i] -= flux.rhou[i] + c_flux.rhou[i];
      rhs[icv0].rhoE -= flux.rhoE + c_flux.rhoE;

      rhs[icv1].rho += flux.rho + c_flux.rho;
      for (int i =0; i < 3; ++i)
        rhs[icv1].rhou[i] += flux.rhou[i] + c_flux.rhou[i];
      rhs[icv1].rhoE += flux.rhoE + c_flux.rhoE;
    }

    // boundary extended faces with no c, no cmpact 
    for (int ief = n_cc_i; ief < n_e__b; ++ief) { 
      const int icv0           = fa_hashed[ief].cvofa[0];
      const int icv1           = fa_hashed[ief].cvofa[1];
      const int ief_compressed = fa_hashed[ief].ief_compressed;

      cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
      cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
      cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);
 
      IdealGasRhs flux;
      double u_avg[3], inv_v_avg, p_avg, h_stag;
      calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
                          c_diss,u_avg,inv_v_avg,h_stag,p_avg);

      addViscousFlux(flux,ef_geom[ief_compressed],x_vv_tmp[icv0],x_vv_tmp[icv1],cv_light[icv0],cv_light[icv1],
                     cv_compact[icv0],cv_compact[icv1]);


      rhs[icv0].rho -= flux.rho;
      for (int i =0; i < 3; ++i) 
        rhs[icv0].rhou[i] -= flux.rhou[i];
      rhs[icv0].rhoE -= flux.rhoE;
    }

    // boundary extended faces with no c, but has compact closure terms
    for (int ief = n_e__b; ief < n_c__b; ++ief) { 
      const int icv0           = fa_hashed[ief].cvofa[0];
      const int icv1           = fa_hashed[ief].cvofa[1];
      const int ief_compressed = fa_hashed[ief].ief_compressed;
      const int icf_compressed = fa_hashed[ief].icf_compressed;

      cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
      cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
      cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);
 
      IdealGasRhs flux;
      double u_avg[3], inv_v_avg,p_avg,h_stag;
      calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
                          c_diss,u_avg,inv_v_avg,h_stag,p_avg);

      addViscousFlux(flux,ef_geom[ief_compressed],x_vv_tmp[icv0],x_vv_tmp[icv1],cv_light[icv0],cv_light[icv1],
                     cv_compact[icv0],cv_compact[icv1]);


      cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
      cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);
 
      // aggregate into the flux...
      calcInternalCompactFluxNew(flux,cf_geom[icf_compressed],cv_light[icv0],cv_light[icv1],
                                 cv_compact[icv0],cv_compact[icv1],u_avg,inv_v_avg,h_stag,c_diss);

      rhs[icv0].rho -= flux.rho;
      for (int i =0; i < 3; ++i) 
        rhs[icv0].rhou[i] -= flux.rhou[i];
      rhs[icv0].rhoE -= flux.rhoE;
    }

    // proc boundary with c, but no compact closures...
    for (int ief = n_c__b; ief < n_ec_b; ++ief) { 
      const int icv0           = fa_hashed[ief].cvofa[0];
      const int icv1           = fa_hashed[ief].cvofa[1];
      const int ief_compressed = fa_hashed[ief].ief_compressed;

      cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
      cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
      cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);
 
      IdealGasRhs flux;
      double u_avg[3],inv_v_avg,p_avg,h_stag;
      calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
                          c_diss,u_avg,inv_v_avg,h_stag,p_avg);

      addViscousFlux(flux,ef_geom[ief_compressed],x_vv_tmp[icv0],x_vv_tmp[icv1],cv_light[icv0],cv_light[icv1],
                     cv_compact[icv0],cv_compact[icv1]);


      IdealGasRhs c_flux;
      calcInternalFluxSymmetricNew(c_flux,ef_geom[ief_compressed],u_avg,inv_v_avg,
                                   h_stag,p_avg,cv_light[icv0],cv_light[icv1]); 

      rhs[icv0].rho -= flux.rho + c_flux.rho;
      for (int i =0; i < 3; ++i)
        rhs[icv0].rhou[i] -= flux.rhou[i] + c_flux.rhou[i];
      rhs[icv0].rhoE -= flux.rhoE + c_flux.rhoE;
    }

    // proc boundary with c and compact closures
    for (int ief = n_ec_b; ief < n_cc_b; ++ief) { 
      const int icv0           = fa_hashed[ief].cvofa[0];
      const int icv1           = fa_hashed[ief].cvofa[1];
      const int ief_compressed = fa_hashed[ief].ief_compressed;
      const int icf_compressed = fa_hashed[ief].icf_compressed;

      cti_prefetch(&rhs[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],1,0);
      cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
      cti_prefetch(&cv_light[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);
 
      IdealGasRhs flux;
      double u_avg[3],inv_v_avg,p_avg,h_stag;
      calcInternalFluxNew(flux,ef_geom[ief_compressed],cv_light[icv0],cv_light[icv1],
                          c_diss,u_avg,inv_v_avg,h_stag,p_avg);

      addViscousFlux(flux,ef_geom[ief_compressed],x_vv_tmp[icv0],x_vv_tmp[icv1],cv_light[icv0],cv_light[icv1],
                     cv_compact[icv0],cv_compact[icv1]);


      IdealGasRhs c_flux;
      calcInternalFluxSymmetricNew(c_flux,ef_geom[ief_compressed],u_avg,inv_v_avg,
                                   h_stag,p_avg,cv_light[icv0],cv_light[icv1]); 

      cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[0]],0,0);
      cti_prefetch(&cv_compact[fa_hashed[ief+PFD_DISTANCE].cvofa[1]],0,0);
 
      // aggregate into the flux...
      calcInternalCompactFluxNew(flux,cf_geom[icf_compressed],cv_light[icv0],cv_light[icv1],
                                 cv_compact[icv0],cv_compact[icv1],u_avg,inv_v_avg,h_stag,c_diss);

      rhs[icv0].rho -= flux.rho + c_flux.rho;
      for (int i =0; i < 3; ++i)
        rhs[icv0].rhou[i] -= flux.rhou[i] + c_flux.rhou[i];
      rhs[icv0].rhoE -= flux.rhoE + c_flux.rhoE;
    }

    addSourceHook(rhs,time,rk_stage);

    delete[] x_vv_tmp;
  }
 /*
  void calcRhs(IdealGasRhs* rhs, double * rhs_sc, const double time, const int rk_stage) { 

    // HACK HACK HACK ... build another gradient.

    // overwrite the velocity field.. 

    for (int icv = 0; icv < ncv_g2; ++icv) { 
      u[icv][0] = 0.5*(1.0 - x_cv[icv][1]*x_cv[icv][1]);
      u[icv][1] = 0.0;
      u[icv][2] = 0.0;
    }

    double (*dudx_g)[3][3] = new double[ncv_g][3][3];
    StaticSolver::calcCvGrad(dudx_g, u);

    MiscUtils::dumpRange(u,ncv,"u");
    MiscUtils::dumpRange(dudx_g,ncv,"dudx");

    // need to update into the ghosts .. 
    
    double (*tmp)[3]  = new double[ncv_g][3];
    for (int i =0; i < 3; ++i) { 
      
      for (int icv = 0; icv < ncv; ++icv) {
        for (int j = 0; j < 3; ++j) 
          tmp[icv][j] = dudx_g[icv][i][j];
      }
      updateCvData(tmp, REPLACE_DATA);

      for (int icv = ncv; icv < ncv_g; ++icv) {
        for (int j = 0; j < 3; ++j)
          dudx_g[icv][i][j] = tmp[icv][j];
      }
    }
    delete[] tmp;
    


    for(vector<IdealGasBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)  
      (*it)->addBoundaryFlux(rhs);

    // we're only going to compute the viscous flux to avoid any ambiguity.. 
    // and this only occurs on the compact faces.

    for (int ifa = 0; ifa < nfa; ++ifa) { 

      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];

      const double area = MAG(n_fa[ifa]);
      const double delta = area/area_over_delta_fa[ifa];
 
      IdealGasRhs flux;
      flux.rho = 0.0;
      for (int i =0 ; i < 3; ++i) 
        flux.rhou[i] = 0.0;
      flux.rhoE = 0.0;

      double unit_n[3]; 
      for (int i =0; i < 3; ++i) 
        unit_n[i] = n_fa[ifa][i]/area;

      cout << " ifa compact_grad " << (u[icv1][0] - u[icv0][0])/delta*unit_n[0] << "    " 
                                   << (u[icv1][0] - u[icv0][0])/delta*unit_n[1] << "    " 
                                   << (u[icv1][0] - u[icv0][0])/delta*unit_n[2] << endl;

      double duidxj_fa[3][3];
      for (int i =0; i < 3; ++i) {
        double duds = (u[icv1][i] - u[icv0][i])/delta;
        double duidxj_cv[3]; FOR_J3 duidxj_cv[j] = 0.5*(dudx_g[icv0][i][j] + dudx_g[icv1][i][j]);
        duidxj_fa[i][0] = (1.0 - unit_n[0]*unit_n[0])*duidxj_cv[0] + unit_n[0]*(duds - unit_n[1]*duidxj_cv[1] - unit_n[2]*duidxj_cv[2]);
        duidxj_fa[i][1] = (1.0 - unit_n[1]*unit_n[1])*duidxj_cv[1] + unit_n[1]*(duds - unit_n[2]*duidxj_cv[2] - unit_n[0]*duidxj_cv[0]);
        duidxj_fa[i][2] = (1.0 - unit_n[2]*unit_n[2])*duidxj_cv[2] + unit_n[2]*(duds - unit_n[0]*duidxj_cv[0] - unit_n[1]*duidxj_cv[1]);
      }
      
      double dukdxk_fa = duidxj_fa[0][0] + duidxj_fa[1][1] + duidxj_fa[2][2];
      const double mu_avg = 0.5*( mu_lam[icv0] + mu_lam[icv1]); 
      
      // inspect .. 
      cout << "ifa, n_fa, dudx, dvdx, dudy, divu: " << ifa << "   " 
            << COUT_VEC(n_fa[ifa]) << "   " << duidxj_fa[0][0] << "    " << duidxj_fa[1][0] << "    " << duidxj_fa[0][1]  << "    " 
            << dukdxk_fa << endl;


      double tauij[3][3];
      for (int i =0; i < 3; ++i) {
        for (int j =0; j < 3; ++j) tauij[i][j] = ( duidxj_fa[i][j] + duidxj_fa[j][i] );
        //for (int j =0; j < 3; ++j) tauij[i][j] = duidxj_fa[i][j];
        tauij[i][i]       -= 2.0/3.0*dukdxk_fa;
      }

      cout << " ifa , tau0j dot nj = " << ifa << "    " << DOT_PRODUCT(tauij[0],n_fa[ifa]);

      double tauijnj[3];
      for (int i =0; i < 3; ++i) { 
        
        // note tauij is currently storing 2(Sij - 1/3*Skk*delta_ij)
        tauijnj[i] = mu_avg*DOT_PRODUCT(tauij[i],unit_n)*area;
      }
      
      double u_avg[3];
      for (int i =0; i < 3; ++i) 
        u_avg[i] = 0.5*(u[icv0][i] + u[icv1][i]);

      for (int i = 0; i < 3; ++i) 
        flux.rhou[i] -= tauijnj[i];
      flux.rhoE -= DOT_PRODUCT(tauijnj, u_avg);
     
      cout << " ==============" << endl;
      getchar();
    }


    delete[] dudx_g;

    throw(0);
  }
  */
 /*
  void calcRhs(IdealGasRhs* rhs, double * rhs_sc, const double time, const int rk_stage) { 

    // HACK HACK HACK ... build another gradient.

    // overwrite the velocity field.. 

    for (int icv = 0; icv < ncv_g2; ++icv) { 
      u[icv][0] = 0.5*(1.0 - x_cv[icv][1]*x_cv[icv][1]);
      u[icv][1] = 0.0;
      u[icv][2] = 0.0;
    }

    double (*dudx_g)[3][3] = new double[ncv_g][3][3];
    StaticSolver::calcCvGrad(dudx_g, u);

    MiscUtils::dumpRange(u,ncv,"u");
    MiscUtils::dumpRange(dudx_g,ncv,"dudx");

    // need to update into the ghosts .. 
    
    double (*tmp)[3]  = new double[ncv_g][3];
    for (int i =0; i < 3; ++i) { 
      
      for (int icv = 0; icv < ncv; ++icv) {
        for (int j = 0; j < 3; ++j) 
          tmp[icv][j] = dudx_g[icv][i][j];
      }
      updateCvData(tmp, REPLACE_DATA);

      for (int icv = ncv; icv < ncv_g; ++icv) {
        for (int j = 0; j < 3; ++j)
          dudx_g[icv][i][j] = tmp[icv][j];
      }
    }
    delete[] tmp;
    


    for(vector<IdealGasBc*>::iterator it = bcs.begin(); it != bcs.end(); ++it)  
      (*it)->addBoundaryFlux(rhs);

    // we're only going to compute the viscous flux to avoid any ambiguity.. 
    // and this only occurs on the compact faces.

    for (int ifa = 0; ifa < nfa; ++ifa) { 

      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];

      const double area = MAG(n_fa[ifa]);
      const double delta = area/area_over_delta_fa[ifa];
 
      IdealGasRhs flux;
      flux.rho = 0.0;
      for (int i =0 ; i < 3; ++i) 
        flux.rhou[i] = 0.0;
      flux.rhoE = 0.0;

      double unit_n[3]; 
      for (int i =0; i < 3; ++i) 
        unit_n[i] = n_fa[ifa][i]/area;

      cout << " ifa compact_grad " << (u[icv1][0] - u[icv0][0])/delta*unit_n[0] << "    " 
                                   << (u[icv1][0] - u[icv0][0])/delta*unit_n[1] << "    " 
                                   << (u[icv1][0] - u[icv0][0])/delta*unit_n[2] << endl;

      double duidxj_fa[3][3];
      for (int i =0; i < 3; ++i) {
        double duds = (u[icv1][i] - u[icv0][i])/delta;
        double duidxj_cv[3]; FOR_J3 duidxj_cv[j] = 0.5*(dudx_g[icv0][i][j] + dudx_g[icv1][i][j]);
        duidxj_fa[i][0] = (1.0 - unit_n[0]*unit_n[0])*duidxj_cv[0] + unit_n[0]*(duds - unit_n[1]*duidxj_cv[1] - unit_n[2]*duidxj_cv[2]);
        duidxj_fa[i][1] = (1.0 - unit_n[1]*unit_n[1])*duidxj_cv[1] + unit_n[1]*(duds - unit_n[2]*duidxj_cv[2] - unit_n[0]*duidxj_cv[0]);
        duidxj_fa[i][2] = (1.0 - unit_n[2]*unit_n[2])*duidxj_cv[2] + unit_n[2]*(duds - unit_n[0]*duidxj_cv[0] - unit_n[1]*duidxj_cv[1]);
      }
      
      double dukdxk_fa = duidxj_fa[0][0] + duidxj_fa[1][1] + duidxj_fa[2][2];
      const double mu_avg = 0.5*( mu_lam[icv0] + mu_lam[icv1]); 
      
      // inspect .. 
      cout << "ifa, n_fa, dudx, dvdx, dudy, divu: " << ifa << "   " 
            << COUT_VEC(n_fa[ifa]) << "   " << duidxj_fa[0][0] << "    " << duidxj_fa[1][0] << "    " << duidxj_fa[0][1]  << "    " 
            << dukdxk_fa << endl;


      double tauij[3][3];
      for (int i =0; i < 3; ++i) {
        for (int j =0; j < 3; ++j) tauij[i][j] = ( duidxj_fa[i][j] + duidxj_fa[j][i] );
        //for (int j =0; j < 3; ++j) tauij[i][j] = duidxj_fa[i][j];
        tauij[i][i]       -= 2.0/3.0*dukdxk_fa;
      }

      cout << " ifa , tau0j dot nj = " << ifa << "    " << DOT_PRODUCT(tauij[0],n_fa[ifa]);

      double tauijnj[3];
      for (int i =0; i < 3; ++i) { 
        
        // note tauij is currently storing 2(Sij - 1/3*Skk*delta_ij)
        tauijnj[i] = mu_avg*DOT_PRODUCT(tauij[i],unit_n)*area;
      }
      
      double u_avg[3];
      for (int i =0; i < 3; ++i) 
        u_avg[i] = 0.5*(u[icv0][i] + u[icv1][i]);

      for (int i = 0; i < 3; ++i) 
        flux.rhou[i] -= tauijnj[i];
      flux.rhoE -= DOT_PRODUCT(tauijnj, u_avg);
     
      cout << " ==============" << endl;
      getchar();
    }


    delete[] dudx_g;

    throw(0);
  }
  */
/*
  void calcRhs(IdealGasRhs* rhs, double * rhs_sc, const double time, const int rk_stage) { 

    // just computing the simple diffusive closure.. without the tranpose of grad(u)

    double (*dudx)[3][3] = new double[ncv_g][3][3];

    for (int icv = 0; icv < ncv_g; ++icv) { 
      
      // specify the exact solution for the viscous transpose operator ... 

      FOR_I3 FOR_J3 dudx[icv][i][j] = 0.0;
      dudx[icv][0][1] = -(1.0/mu_ref)*x_cv[icv][1]; 

    } 

    const double cp = R_gas*gamma/(gamma-1.0);
    for (int ifa = 0; ifa < nfa; ++ifa) { 

      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];

      IdealGasRhs flux;
      flux.rho  = 0.0;
      flux.rhoE = 0.0;

      //for (int i = 0; i < 3; ++i) { 
      //  flux.rhou[i] = -mu_ref * (u[icv1][i] - u[icv0][i]) * area_over_delta_fa[ifa];
      //}

      // apply the generalized tau_ij = 2mu(sij - 1/3Skk*delta_ij)

      const double fa_area = MAG(n_fa[ifa]);
      const double delta   = fa_area/area_over_delta_fa[ifa];
      
      double unit_n[3];
      for (int i =0; i <3; ++i)
        unit_n[i] = n_fa[ifa][i]/fa_area;
      
      double duidxj_fa[3][3];
      for (int i =0; i < 3; ++i) {
        double duds = (u[icv1][i] - u[icv0][i])/delta;
        double duidxj_cv[3]; FOR_J3 duidxj_cv[j] = 0.5*(dudx[icv0][i][j] + dudx[icv1][i][j]); 
        duidxj_fa[i][0] = (1.0 - unit_n[0]*unit_n[0])*duidxj_cv[0] + unit_n[0]*(duds - unit_n[1]*duidxj_cv[1] - unit_n[2]*duidxj_cv[2]);
        duidxj_fa[i][1] = (1.0 - unit_n[1]*unit_n[1])*duidxj_cv[1] + unit_n[1]*(duds - unit_n[2]*duidxj_cv[2] - unit_n[0]*duidxj_cv[0]);
        duidxj_fa[i][2] = (1.0 - unit_n[2]*unit_n[2])*duidxj_cv[2] + unit_n[2]*(duds - unit_n[0]*duidxj_cv[0] - unit_n[1]*duidxj_cv[1]);
      }
  
      double dukdxk_fa = duidxj_fa[0][0] + duidxj_fa[1][1] + duidxj_fa[2][2];

      double tauij[3][3];
      for (int i =0; i < 3; ++i) {
        for (int j =0; j < 3; ++j) tauij[i][j] = ( duidxj_fa[i][j] + duidxj_fa[j][i] );
        tauij[i][i]       -= 2.0/3.0*dukdxk_fa;
      }

      double tauijnj[3];
      for (int i =0; i < 3; ++i) { 
        // note tauij is currently storing 2(Sij - 1/3*Skk*delta_ij)
        tauijnj[i] = mu_ref*DOT_PRODUCT(tauij[i],n_fa[ifa]);
        flux.rhou[i] = -tauijnj[i];
      }

      // total energy equation terms...

      double u_avg[3];
      for (int i =0; i < 3; ++i) 
        u_avg[i] = 0.5*(u[icv0][i] + u[icv1][i]);

      flux.rhoE += DOT_PRODUCT(u_avg,flux.rhou);
      flux.rhoE -= mu_ref*cp/Pr_lam*(T[icv1] - T[icv0]) * area_over_delta_fa[ifa];

      rhs[icv0].rho -= flux.rho;
      for (int i = 0; i < 3; ++i) 
        rhs[icv0].rhou[i] -= flux.rhou[i];
      rhs[icv0].rhoE -= flux.rhoE;

      if ( icv1 < ncv) { 
        rhs[icv1].rho += flux.rho;
        for (int i = 0; i < 3; ++i) 
          rhs[icv1].rhou[i] += flux.rhou[i];
        rhs[icv1].rhoE += flux.rhoE;
      }
    }

    // add the pressure flux on the extended faces
    if ( false) { 

      for (int ief = 0; ief < nef; ++ief) { 

        const int icv0 = cvoef[ief][0];
        const int icv1 = cvoef[ief][1];

        IdealGasRhs flux;
        const double p_avg = 0.5*(p[icv0] + p[icv1]);
        double u_avg[3];
        for (int i =0; i < 3; ++i) 
          u_avg[i] = 0.5*(u[icv0][i] + u[icv1][i]);

        flux.rho = 0.0;
        for (int i =0; i < 3; ++i) 
          flux.rhou[i] = p_avg * n_ef[ief][i];

        flux.rhoE = DOT_PRODUCT(u_avg,flux.rhou);

        rhs[icv0].rho -= flux.rho;
        for (int i = 0; i < 3; ++i) 
          rhs[icv0].rhou[i] -= flux.rhou[i];
        rhs[icv0].rhoE -= flux.rhoE;
        
        if ( icv1 < ncv) { 
          rhs[icv1].rho += flux.rho;
          for (int i = 0; i < 3; ++i) 
            rhs[icv1].rhou[i] += flux.rhou[i];
          rhs[icv1].rhoE += flux.rhoE;
        }
      }

    } 

    // boundary conditions... everyone is a cold isothermal wall...
    for (int ibf = 0; ibf < nbf; ++ibf) { 
      const int icv = cvobf[ibf];
      for (int i =0 ; i < 3; ++i) 
        rhs[icv].rhou[i] -= mu_ref*(u[icv][i] - 0.0)*area_over_delta_bf[ibf];
      rhs[icv].rhoE -= mu_ref/Pr_lam*cp*(T[icv] - 1.0)*area_over_delta_bf[ibf];
    }

    // add the forcing term.. 
    for (int icv = 0; icv < ncv; ++icv) { 
      rhs[icv].rhou[0] += 1.0 * vol_cv[icv];
      rhs[icv].rhoE += u[icv][0]* vol_cv[icv];
    }

    delete[] dudx;
  }
  */
struct MM{ 
  double val;
  int ind;
} ;

class MyIdealGasSolver : public IdealGasSolver { 
public: 

  double * dil;

  MyIdealGasSolver() { 

    dil = NULL; registerCvData(dil,"dil",0);
  }
  
  void initialHook() { 
    
    if ( checkParam("DEBUG_IC")) { 
      const double rho_ref = getDoubleParam("RHO_REF");
      const double p_ref   = getDoubleParam("P_REF");
      const double gamma   = getDoubleParam("GAMMA", 1.4);
      
      const double eps = 1.0E-12;

      for (int icv = 0; icv < ncv; ++icv) { 
        const double sos0      = sqrt(gamma*p_ref/rho_ref);
        const double rho_prime = eps*double(rand())/double(RAND_MAX);
        rho[icv] = rho_ref + rho_prime;
        FOR_I3 { 
          u[icv][i] = eps*2.0*(2.0*double(rand())/double(RAND_MAX) - 1.0);
        }
        rhoE[icv]  = (p_ref + sos0*sos0*rho_prime)/(gamma-1.0);
        rhoE[icv] += 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);
      }

      dumpRange(rho,ncv,"rho_init");
      dumpRange(u,ncv,"u_init");
      dumpRange(rhoE,ncv,"rhoE_init");
    }
    
  }

  void temporalHook() { 

    if ( !dil) { 
      dil = new double[ncv];
    }

    FOR_ICV { 
      dil[icv] = dudx[icv][0][0] + dudx[icv][1][1] + dudx[icv][2][2];
    }

    if ( step%check_interval == 0) { 

      double * vv = new double[ncv];
      double * vv2 = new double[ncv];
      double (*x_vv_tmp)[3] = new double[ncv_g2][3];
      
      if (true) { 

        for (int icv = 0; icv < ncv; ++icv) 
          vv[icv] = 0.0;

        for (int icv = 0; icv < ncv; ++icv) 
          vv2[icv] = 0.0;

        for (int icv = 0; icv < ncv; ++icv) 
          for (int i = 0; i < 3; ++i) 
            x_vv_tmp[icv][i] = x_vv[icv][i];

        updateCv2Data(x_vv_tmp,REPLACE_TRANSLATE_DATA);

        const double cp = R_gas*gamma/(gamma-1.0);

        for (int ifa = 0; ifa < nfa; ++ifa) { 

          const int icv0 = cvofa[ifa][0];
          const int icv1 = cvofa[ifa][1];

          const double area     = MAG(n_fa[ifa]);
          const double delta    = area/area_over_delta_fa[ifa];
          const double dp       = p[icv1] - p[icv0];
          const double dh       = cp*(T[icv1] - T[icv0]);
          const double v_avg    = 0.5*(1.0/rho[icv0] + 1.0/rho[icv1]);
          const double beta_avg = 0.5*(1.0/(R_gas*T[icv0]) + 1.0/(R_gas*T[icv1]));
          const double dsor     = cv_compact[icv1].sor - cv_compact[icv0].sor;
          const double gr       = -dsor + beta_avg*(dh - v_avg*dp);
         
          const double c_avg    = fabs(gr)*area*delta/6.0;
          vv[icv0]           += c_avg;
          if ( icv1 < ncv) 
            vv[icv1]         += c_avg;

        }

        for (int icv = 0; icv < ncv; ++icv) 
          vv[icv] /= vol_cv[icv];


        for (int ief = 0; ief < nef; ++ief) { 

          const int icv0 = cvoef[ief][0];
          const int icv1 = cvoef[ief][1];

          // not sure how to average this exactly ... 

          const double area     = MAG(n_ef[ief]);
          const double delta    = DIST(x_vv_tmp[icv1],x_vv_tmp[icv0]); 
          const double dp       = p[icv1] - p[icv0];
          const double dh       = cp*(T[icv1] - T[icv0]);
          const double v_avg    = 0.5*(1.0/rho[icv0] + 1.0/rho[icv1]);
          const double beta_avg = 0.5*(1.0/(R_gas*T[icv0]) + 1.0/(R_gas*T[icv1]));
          const double dsor     = cv_compact[icv1].sor - cv_compact[icv0].sor;
          const double gr       = -dsor + beta_avg*(dh - v_avg*dp);
         
          const double c_avg    = fabs(gr)*area*delta/6.0;
          vv2[icv0]           += c_avg;
          if ( icv1 < ncv) 
            vv2[icv1]         += c_avg;
        }


        for (int icv = 0; icv < ncv; ++icv)
          vv2[icv] /= vol_cv[icv];

      }


      MM tmp; 
      tmp.val = 1.0e+16;
      tmp.ind = -1;
      int icv_min = -1;

      for (int icv = 0; icv < ncv; ++icv) { 
        if ( u[icv][0] < tmp.val) { 
          tmp.val = u[icv][0];
          tmp.ind = mpi_rank;
          icv_min = icv;

        }
      }

      MM tmp_global;

      MPI_Allreduce(&tmp,&tmp_global,1,MPI_DOUBLE_INT,MPI_MINLOC,mpi_comm);

      if ( tmp_global.ind == mpi_rank ) { 

        const double divu = dudx[icv_min][0][0] + dudx[icv_min][1][1] + dudx[icv_min][2][2];

        cout << " >> min vel, gr : " << COUT_VEC(x_cv[icv_min]) << "    " << COUT_VEC(u[icv_min]) << "   " << vv[icv_min] << "    " << vv2[icv_min] << "   " << pow(vol_cv[icv_min],1.0/3.0)*divu/sos[icv_min] << endl;
        cout << " >> rho, Ma, p : " << rho[icv_min] << "   " << MAG(u[icv_min])/sos[icv_min] << "     " << p[icv_min] << "    " << T[icv_min] << endl;
        //cout << " >> div_norm   : " << pow(vol_cv[icv_min],1.0/3.0)*divu/sos[icv_min] << endl;
        cout.flush();
      }

      delete[] vv;
      delete[] vv2;
      delete[] x_vv_tmp;

      MPI_Barrier(mpi_comm);

    } 


  } 

};
void calcInternalFluxNew(IdealGasRhs& flux, const FaNC& fa, const IdealGasState& cv0, 
                         const IdealGasState& cv1, const double& c_diss, double u_avg[3],
                         double& inv_v_avg, double& h_stag, double& p_avg) { 
  
  inv_v_avg =   2.0/(cv0.sp_vol + cv1.sp_vol); // specific volume 
  p_avg     =   0.5*(cv1.p + cv0.p); 

  const double rho_avg = 0.5*(1.0/cv0.sp_vol + 1.0/cv1.sp_vol);
  const double beta_avg = 0.5*(1.0/(cv0.sp_vol*cv0.p) + 1.0/(cv1.sp_vol*cv1.p));
  const double p_tilde  = rho_avg/beta_avg;
  const double v_avg    = 0.5*(cv0.sp_vol + cv1.sp_vol);

  double un_avg = 0.0;
  double u0u1   = 0.0;
  for (int i =0; i < 3; ++i) {
    double uu  = 0.5*(cv0.u[i] + cv1.u[i]); 
    un_avg      += uu * fa.n[i];
    u0u1        += cv0.u[i] * cv1.u[i]; 
    u_avg[i]     = uu;
    flux.rhou[i] = p_tilde*fa.n[i];
  }
  
  double Frho =     un_avg*inv_v_avg;
  h_stag      =   cv1.h + cv0.h; 
  flux.rho    =   Frho; 
  
  for (int i =0; i < 3; ++i)   
    flux.rhou[i]     += Frho*u_avg[i];

  h_stag     = 0.5*(h_stag + u0u1);
  flux.rhoE  =  Frho*h_stag;
  flux.rhoE +=  Frho*v_avg*(p_tilde - p_avg);
}

