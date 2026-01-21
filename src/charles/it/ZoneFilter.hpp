#ifndef ZONEFILTER_HPP
#define ZONEFILTER_HPP

class BfFilter { 

public: 

  StaticSolver* solver;
  int *faofi; // face-of-filter
  double * lfa_fi;
  int nfi;    // nfilter-internal

  int *bfofi; // boundary-face of filter
  int nfb;    // nfilter-boundary

  BfFilter(StaticSolver* _s, const int* bf_flag) : solver(_s) {
  
    int* cv_flag = new int[solver->ncv_g];
    for (int icv = 0; icv < solver->ncv; ++icv) 
      cv_flag[icv] = -1;
  
    faofi = NULL;
    bfofi = NULL;

    for (int iter = 0; iter < 2; ++iter) { 

      nfb = 0;

      for (int ibf = 0; ibf < solver->nbf; ++ibf) {
        if ( bf_flag[ibf] > 0) { 
          if ( iter == 0 ) { 
            ++nfb;
          } else { 
            cv_flag[solver->cvobf[ibf]] = 1;
            bfofi[nfb++] = ibf;
          }
        }
      }

      if ( iter == 0) { 
        assert( bfofi == NULL); 
        bfofi = new int[nfb];
      }
    }
    
    solver->updateCvData(cv_flag);

    //
    // build the list of internal faces that touch 2 flagged cvs.. 
    // 
   
    for (int iter = 0; iter < 2; ++iter) { 

      nfi = 0;

      for (int ifa = 0; ifa < solver->nfa; ++ifa) {

        const int icv0 = solver->cvofa[ifa][0];
        const int icv1 = solver->cvofa[ifa][1];

        if ( (cv_flag[icv0] > 0) && (cv_flag[icv1] > 0)) { 

          if ( iter == 0) { 
            ++nfi;
          } else { 
            faofi[nfi++] = ifa;
          }
        }
      }

      if ( iter == 0) { 

        assert( faofi == NULL);
        faofi = new int[nfi];
      
      }
    }

    delete[] cv_flag;

    // set an initial guess for the lfa_fi .. should be exposed ... 
    
    lfa_fi = new double[nfi];

    for (int ifi = 0; ifi < nfi; ++ifi) { 

      const int ifa     = faofi[ifi];
      const int icv0    = solver->cvofa[ifa][0];
      const int icv1    = solver->cvofa[ifa][1];

      const double del2 = 0.5*( pow(solver->vol_cv[icv0],2.0/3.0) + 
                                pow(solver->vol_cv[icv1],2.0/3.0));

      lfa_fi[ifi]       = 0.1*del2;
    }

  }

  ~BfFilter() {

    DELETE(faofi);
    DELETE(bfofi);
    DELETE(lfa_fi);

  }


  void doFilter(double* phi_f, const double* phi) const {

    double * phi_cv = new double[solver->ncv_g];

    //MiscUtils::dumpRange(phi,solver->nbf,"BEFORE");

    for (int icv = 0; icv < solver->ncv; ++icv) 
      phi_cv[icv] = 0.0;

    for (int ibf = 0; ibf < solver->nbf; ++ibf) { 
      phi_cv[solver->cvobf[ibf]] = phi[ibf];
      phi_f[ibf]                 = 0.0;
    }

    solver->updateCvData(phi_cv);
    
    double* tmp_cv = new double[solver->ncv_g];

    for (int icv = 0; icv < solver->ncv_g; ++icv) 
      tmp_cv[icv] = 0.0;

    for (int ifi = 0; ifi < nfi; ++ifi) { 

      const int ifa  = faofi[ifi];
      const int icv0 = solver->cvofa[ifa][0];
      const int icv1 = solver->cvofa[ifa][1];
      
      // l_fi is of length^2 units... this is a diffusion coefficient.
      const double flux = lfa_fi[ifi]*(phi_cv[icv1] - phi_cv[icv0])*solver->area_over_delta_fa[ifa];
      tmp_cv[icv0]     += flux;
      tmp_cv[icv1]     -= flux;
    
    }

    for (int ifb = 0; ifb < nfb; ++ifb) { 

      const int ibf    = bfofi[ifb];
      const int icv    = solver->cvobf[ibf];
      phi_f[ibf]       = phi[ibf] + tmp_cv[icv]/solver->vol_cv[icv];

    }

    //MiscUtils::dumpRange(phi_f,solver->nbf,"AFTER");

    delete[] tmp_cv;
    delete[] phi_cv;
  }

};

#endif
