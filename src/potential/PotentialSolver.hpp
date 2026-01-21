#ifndef _POTENTIAL_SOLVER_HPP_
#define _POTENTIAL_SOLVER_HPP_

#include "StaticSolver.hpp"


class PotentialSolver : public StaticSolver {
  
private:

protected:

  double * phi;
  double * rhoun_bf;
  double * c;
  double * c_bf;

public:
  
  PotentialSolver() {
    
    COUT1("PotentialSolver()");

    // variables...
    
    phi = NULL; registerCvData(phi,"phi",READWRITE_DATA);
    c = NULL; registerCvData(c,"c",READWRITE_DATA);

    // no need to register rhoun_bf -- it gets set before solve based on bc...
    
    rhoun_bf = NULL; 
    c_bf = NULL;

  }
  
  virtual ~PotentialSolver() {
    
    COUT1("~PotentialSolver()");
    
    // variable cleanup...
    
    DELETE(phi);
    DELETE(rhoun_bf);
    DELETE(c);
    DELETE(c_bf);

  }
  
  void initData() {

    assert(c == NULL); c = new double[ncv_g];
    assert(phi == NULL); phi = new double[ncv_g];
    assert(rhoun_bf == NULL); rhoun_bf = new double[nbf];
    assert(c_bf == NULL); c_bf = new double[nbf];

  }

  void initBoundaryConditions() {
    StaticSolver::initBoundaryConditions();
  }

  void initialHook() {

    // if the data was not read in, initialize...
    
    if (!CtiRegister::checkDataFlag("phi")) {
      COUT1(" > setting phi and c to zero...");
      FOR_ICV_G phi[icv] = 0.0;
      FOR_ICV_G c[icv] = 0.0;
    }
    else {
      updateCvData(phi);
      updateCvData(c);
    }
    
    // bcs...

    if (mpi_rank == 0)
      cout << "Applying bcs (MS first)..." << endl;
    
    // zero rhoun_bf and c_bf first...
    FOR_IBF rhoun_bf[ibf] = 0.0;
    FOR_IBF c_bf[ibf] = 0.0;
    
    // cycle through non-outlets first...
    double inlet_mdot_sum = 0.0;
    double inlet_mdotc_sum = 0.0;
    double outlet_rhoun = 0.0;
    int ierr = 0;

    bool found = false;
    FOR_IZONE(bfZoneVec) {
      if (Param * param = getParam(bfZoneVec[izone].getName())) {
        found = true;
	const string bc_name = param->getString(0);
	if (bc_name == "MS") {
	  // MS = mass, scalar (2 value required)
	  const double mdot = param->getDouble(1);
	  inlet_mdot_sum += mdot;
	  const double c = param->getDouble(2);
	  inlet_mdotc_sum += mdot*c;
	  const double rhoun = mdot/bfZoneVec[izone].area_global;
	  if (mpi_rank == 0)
	    cout << " > zone \"" << bfZoneVec[izone].getName() << "\" MS=" << mdot << " " << c << " (rhoun=" << rhoun << ")" << endl;
	  for (int ibf = bfZoneVec[izone].ibf_f; ibf <= bfZoneVec[izone].ibf_l; ++ibf) {
	    rhoun_bf[ibf] = -rhoun*area_bf[ibf]; // inlet negative by convention
	    c_bf[ibf] = c;
	  }
	}
	else if (bc_name == "OUTLET") {
	  // accumulate areas in outlet_rhoun...
	  outlet_rhoun += bfZoneVec[izone].area_global;
	}
	else {
	  if (mpi_rank == 0)
	    cout << "Error: unrecognized bc: " << bc_name << endl;
	  ierr = -1;
	}
      }
      /*
      else {
	if (mpi_rank == 0)
	  cout << " > zone \"" << bfZoneVec[izone].getName() << "\" SLIP (default)" << endl;
      }
      */
    }
    if (!found) {
      CWARN("No bc zones we found.");
    }

    if (ierr != 0) {
      CERR("bc problem");
    }

    // note that no reduction is necessary, because the zone area info is available 
    // already reduced...
    
    if (inlet_mdot_sum != 0.0) {
      if (outlet_rhoun == 0.0) { // recall this is holding the area right now -- should be positive
	CERR("Inlets do not balance. You must specify one or more zones as OUTLET");
      }
      outlet_rhoun = inlet_mdot_sum/outlet_rhoun; // recall outlet_rhoun contained the outlet area
    }
    
    if (mpi_rank == 0)
      cout << "Outlet(s) second (total outlet mdot=" << inlet_mdot_sum << ", outlet mdot c estimate=" << inlet_mdotc_sum << ")..." << endl;
    
    FOR_IZONE(bfZoneVec) {
      if (Param * param = getParam(bfZoneVec[izone].getName())) {
	const string bc_name = param->getString();
	if (bc_name == "OUTLET") {
	  if (mpi_rank == 0)
	    cout << " > zone \"" << bfZoneVec[izone].getName() << "\" OUTLET (computed MDOT=" << outlet_rhoun*bfZoneVec[izone].area_global << ", rhoun=" << outlet_rhoun << ")" << endl;
	  for (int ibf = bfZoneVec[izone].ibf_f; ibf <= bfZoneVec[izone].ibf_l; ++ibf)
	    rhoun_bf[ibf] = outlet_rhoun*area_bf[ibf]; // outlet positive by convention
	}
      }
    }

  }

  void init() {

    COUT1("PotentialSolver::init()");
    StaticSolver::init(INIT_COMPACT_FACES|INIT_CV_GRAD);

    if (mpi_rank == 0) 
      logger->setKillFilename("killcharles");

  }
  
  void run() {

    COUT1("PotentialSolver::run()");

    initialHook();

    // initialHook could of changed data...
    CtiRegister::clearCurrentData();

    // tell the central log a solver is running...
    if (mpi_rank==0) logger->insertSolverAction(SIMULATION_RUN);

    // build the rhs...
    
    double * rhs = new double[ncv];
    for (int icv = 0; icv < ncv; ++icv)
      rhs[icv] = 0.0;

    for (int ibf = 0; ibf < nbf; ++ibf) {
      const int icv = cvobf[ibf];
      rhs[icv] += rhoun_bf[ibf];
    }

    // check rhs...

    double my_sum = 0.0;
    for (int icv = 0; icv < ncv; ++icv)
      my_sum += rhs[icv];
    double sum;
    MPI_Reduce(&my_sum,&sum,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    if (mpi_rank == 0)
      cout << " > sum(rhs) (should be zero): " << sum << endl;

    // build a crs grid description of the mesh using the partition... 

#ifdef CRS_STUFF

    int nfa_crs = 0;
    for (int ii = 0; ii < cvPrcommVec.size(); ++ii)
      if (cvPrcommVec[ii].getRank() != mpi_rank)
	++nfa_crs;

    double * coeff_buf = new double[nfa_crs];
    int * rank_buf = new int[nfa_crs];
    for (int ifa_crs = 0; ifa_crs < nfa_crs; ++ifa_crs) {
      coeff_buf[ifa_crs] = 0.0;
      rank_buf[ifa_crs] = -1;
    }

    {
      map<const int,int> rankMap;
      for (int ifa = nfa_i; ifa < nfa; ++ifa) {
	const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
	const int icv1 = cvofa[ifa][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g));
	int rank,bits,index;
	BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv1-ncv]);
	assert(rank != mpi_rank);
	assert(bits == 0);
	map<const int,int>::const_iterator iter = rankMap.find(rank);
	int ifa_crs;
	if (iter == rankMap.end()) {
	  ifa_crs = rankMap.size();
	  assert(rank_buf[ifa_crs] == -1);
	  rank_buf[ifa_crs] = rank;
	  rankMap[rank] = ifa_crs;
	}
	else {
	  ifa_crs = iter->second;
	}
	const double coeff = fa[ifa].area/DIST(x_cv[icv0],x_cv[icv1]);
	coeff_buf[ifa_crs] += coeff;
      }
      assert(rankMap.size() == nfa_crs);
    }

    // build the rhs...

    double rhs_crs = 0;

    for (int ibf = 0; ibf < nbf; ++ibf) {
      const int icv = cvobf[ibf];
      rhs_crs += rhoun_bf[ibf];
    }
    
    for (int ifa = nfa_i; ifa < nfa; ++ifa) {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g));
      const double coeff = fa[ifa].area/DIST(x_cv[icv0],x_cv[icv1]);
      rhs_crs -= coeff*(phi[icv1]-phi[icv0]);
    }

    //cout << "rank: " << mpi_rank << " rhs_crs first: " << rhs_crs << endl;
    
    /*
      FOR_RANK {
      if (mpi_rank == rank) {
      cout << "rank: " << mpi_rank << " nfa_crs: " << nfa_crs << endl;
      for (int ifa_crs = 0; ifa_crs < nfa_crs; ++ifa_crs) {
      cout << " > rank_buf: " << rank_buf[ifa_crs] << " coeff_buf: " << coeff_buf[ifa_crs] << endl;
      }
      }
      MPI_Pause("OK");
      }
    */

    // reduce to one processor...
    // in the future might want to use a subset of the processors, but for now...
    
    int * nfa_crs_full = NULL;
    if (mpi_rank == 0)
      nfa_crs_full = new int[mpi_size];
    MPI_Gather(&nfa_crs,1,MPI_INT,nfa_crs_full,1,MPI_INT,0,mpi_comm);
    
    int * disp_crs_full = NULL;
    int * rank_buf_full = NULL;
    double * coeff_buf_full = NULL;
    
    int full_size;
    if (mpi_rank == 0) {
      disp_crs_full = new int[mpi_size];
      disp_crs_full[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
	disp_crs_full[rank] = disp_crs_full[rank-1] + nfa_crs_full[rank-1];
      full_size = disp_crs_full[mpi_size-1] + nfa_crs_full[mpi_size-1];
      rank_buf_full = new int[full_size];
      coeff_buf_full = new double[full_size];
    }
    
    MPI_Gatherv(rank_buf,nfa_crs,MPI_INT,rank_buf_full,nfa_crs_full,disp_crs_full,MPI_INT,0,mpi_comm);
    MPI_Gatherv(coeff_buf,nfa_crs,MPI_DOUBLE,coeff_buf_full,nfa_crs_full,disp_crs_full,MPI_DOUBLE,0,mpi_comm);
    
    int * cvocv_i_crs = NULL;
    int * cvocv_v_crs = NULL;
    double * A_crs_full = NULL;
    
    if (mpi_rank == 0) {

      cvocv_i_crs = new int[mpi_size+1];
      cvocv_i_crs[0] = 0;

      cvocv_v_crs = new int[mpi_size+full_size]; // add space for the diagonal
      A_crs_full = new double[mpi_size+full_size];
      
      int coc = 0;
      FOR_RANK {
	const int coc0 = coc++;
	cvocv_v_crs[coc0] = rank;
	A_crs_full[coc0] = 0.0;
	for (int ii = disp_crs_full[rank]; ii < disp_crs_full[rank]+nfa_crs_full[rank]; ++ii) {
	  cvocv_v_crs[coc] = rank_buf_full[ii];
	  assert(rank_buf_full[ii] != rank);
	  A_crs_full[coc] = coeff_buf_full[ii];
	  A_crs_full[coc0] -= A_crs_full[coc];
	  ++coc;
	}
	cvocv_i_crs[rank+1] = coc;
      }
      assert(coc == mpi_size+full_size);
      
    }
    
    double * rhs_crs_full = NULL;
    if (mpi_rank == 0)
      rhs_crs_full = new double[mpi_size];
    MPI_Gather(&rhs_crs,1,MPI_DOUBLE,rhs_crs_full,1,MPI_DOUBLE,0,mpi_comm);


    double * phi_crs_full = NULL;

    if (mpi_rank == 0) {

      phi_crs_full = new double[mpi_size];
      FOR_RANK phi_crs_full[rank] = 0.0;
      
      FOR_RANK {
	cout << " rank: " << rank << " rhs_crs_full: " << rhs_crs_full[rank] << " coeff: ";
	for (int coc = cvocv_i_crs[rank]; coc != cvocv_i_crs[rank+1]; ++coc)
	  cout << cvocv_v_crs[coc] << ":" << A_crs_full[coc] << " ";
	cout << endl;
      }
      
      solveCgLocal(phi_crs_full,mpi_size,A_crs_full,rhs_crs_full,cvocv_i_crs,cvocv_v_crs);
      
    }
      
    MPI_Scatter(phi_crs_full,1,MPI_DOUBLE,&rhs_crs,1,MPI_DOUBLE,0,mpi_comm);
    
    for (int icv = 0; icv < ncv; ++icv)
      phi[icv] += rhs_crs;
    updateCvData(phi);

    // now check that the partitions EXACTLY satisfy mass by recomputing the rhs...

    {
    
      double rhs_crs = 0;
      for (int ibf = 0; ibf < nbf; ++ibf) {
	const int icv = cvobf[ibf];
	rhs_crs += rhoun_bf[ibf];
      }
      
      for (int ifa = nfa_i; ifa < nfa; ++ifa) {
	const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
	const int icv1 = cvofa[ifa][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g));
	const double coeff = fa[ifa].area/DIST(x_cv[icv0],x_cv[icv1]);
	rhs_crs -= coeff*(phi[icv1]-phi[icv0]);
      }
      
      //cout << "rank: " << mpi_rank << " rhs_crs after: " << rhs_crs << endl;

    }

#endif

    /*
    if (checkParam("READ_PHI")) {
      readCvDataCheckpoint(phi,ncv,"phi");
      updateCvData(phi);
    }
    */

    // build the Laplacian operator...
    
    double * A = new double[cvocv_i[ncv]];
    buildCvLaplacian(A);
    
    solveCvCg(phi,A,rhs);

    delete[] A;
    delete[] rhs;
    
    //if (checkParam("WRITE_PHI"))
    //writeCvDataCheckpoint(phi,ncv,"phi");

    MiscUtils::dumpRange(phi,ncv,"phi");
    
    // now solve the c scalar equation assuming first order upwind...
    
    if (mpi_rank == 0)
      cout << "solving c scalar..." << endl;

    // step 1: build a rank-based order that sweeps from largest potential to
    // smallest potential...

    const double dc_zero = 1.0E-12;
      
    {
      
      int * sweep_order = new int[ncv];
      {
	vector<pair<double,int> > diPairVec(ncv);
	for (int icv = 0; icv < ncv; ++icv) {
	  diPairVec[icv].first = -phi[icv]; // we want largest to smallest in phi, so sort based on -phi
	  diPairVec[icv].second = icv;
	}
	sort(diPairVec.begin(),diPairVec.end());
	for (int icv = 0; icv < ncv; ++icv) {
	  sweep_order[icv] = diPairVec[icv].second;
	}
      }
      
      double * mdot_fa = new double[nfa];
      for (int ifa = 0; ifa < nfa; ++ifa) {
	const int icv0 = cvofa[ifa][0];
	const int icv1 = cvofa[ifa][1];
	mdot_fa[ifa] = area_over_delta_fa[ifa]*(phi[icv0]-phi[icv1]); // note sign here: decreasing potential = positive mass
      }
      
      int iter = 0;
      int done = 0;
      while (done == 0) {

	++iter;
	
	double my_dc_max = 0.0;
	for (int ii = 0; ii < ncv; ++ii) {
	  
	  const int icv = sweep_order[ii];
	  double sum_mdot[2] = { 0.0, 0.0 };
	  
	  for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
	    const int ibf = bfocv_v[boc];
	    if (rhoun_bf[ibf] <= 0.0) {
	      // negative rhoun_bf means INFLOW to icv...
	      sum_mdot[0] += rhoun_bf[ibf]*c_bf[ibf];
	      sum_mdot[1] += rhoun_bf[ibf];
	    }
	  }
	  
	  for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc) {
	    const int ifa = faocv_v[foc];
	    if ((mdot_fa[ifa] < 0.0)&&(cvofa[ifa][0] == icv)) {
	      // this flux is leaving icv1 and going to icv0...
	      sum_mdot[0] += mdot_fa[ifa]*c[cvofa[ifa][1]];
	      sum_mdot[1] += mdot_fa[ifa];
	    }
	    else if ((mdot_fa[ifa] > 0.0)&&(cvofa[ifa][1] == icv)) {
	      // this flux is leaving icv0 and going to icv1...
	      sum_mdot[0] -= mdot_fa[ifa]*c[cvofa[ifa][0]];
	      sum_mdot[1] -= mdot_fa[ifa];
	    }
	  }
	  
	  if (sum_mdot[1] != 0.0) {
	    const double c_new = sum_mdot[0]/sum_mdot[1];
	    my_dc_max = max(my_dc_max,fabs(c[icv]-c_new));
	    c[icv] = c_new;
	  }
	  
	}

	updateCvData(c);

	double dc_max;
	MPI_Reduce(&my_dc_max,&dc_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
	if (mpi_rank == 0) {
	  cout << " > solveC iter: " << iter << " dc_max: " << dc_max << endl;
	  if (dc_max <= dc_zero)
	    done = 1;
	}
	MPI_Bcast(&done,1,MPI_INT,0,mpi_comm);

      }
      
      delete[] sweep_order;
      delete[] mdot_fa;

    }

    // old jacobi solver...

    {
      

      double (*sum_mdot)[2] = new double[ncv][2];
      
      int iter = 0;
      int done = 0;
      while (done == 0) {

	++iter;
	
	for (int icv = 0; icv < ncv; ++icv) {
	  sum_mdot[icv][0] = 0.0;
	  sum_mdot[icv][1] = 0.0;
	}

	for (int ibf = 0; ibf < nbf; ++ibf) {
	  const int icv = cvobf[ibf];
	  if (rhoun_bf[ibf] <= 0.0) {
	    // negative rhoun_bf means INFLOW to icv...
	    sum_mdot[icv][0] += rhoun_bf[ibf]*c_bf[ibf];
	    sum_mdot[icv][1] += rhoun_bf[ibf];
	  }
	}
	
	for (int ifa = 0; ifa < nfa; ++ifa) {
	  const int icv0 = cvofa[ifa][0];
	  const int icv1 = cvofa[ifa][1];
	  //const double delta = DIST(x_cv[icv0],x_cv[icv1]); 
	  const double mdot = area_over_delta_fa[ifa]*(phi[icv0]-phi[icv1]); // note sign here: decreasing potential = positive mass
	  if (mdot >= 0.0) {
	    // this flux is leaving icv0 and going to icv1...
	    if (icv1 < ncv) {
	      sum_mdot[icv1][0] -= mdot*c[icv0];
	      sum_mdot[icv1][1] -= mdot;
	    }
	  }
	  else {
	    // this flux is leaving icv1 and going to icv0...
	    sum_mdot[icv0][0] += mdot*c[icv1];
	    sum_mdot[icv0][1] += mdot;
	  }
	}
	
	double my_dc_max = 0.0;
	for (int icv = 0; icv < ncv; ++icv) {
	  if (sum_mdot[icv][1] != 0.0) { 
	    const double c_new = sum_mdot[icv][0]/sum_mdot[icv][1];
	    my_dc_max = max(my_dc_max,fabs(c[icv]-c_new));
	    c[icv] = c_new;
	  }
	}
	updateCvData(c);

	double dc_max;
	MPI_Reduce(&my_dc_max,&dc_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
	if (mpi_rank == 0) {
	  cout << " > solveC iter: " << iter << " dc_max: " << dc_max << endl;
	  if (dc_max <= dc_zero)
	    done = 1;
	}
	MPI_Bcast(&done,1,MPI_INT,0,mpi_comm);
	
      }

      delete[] sum_mdot;
      
    }

    MiscUtils::dumpRange(c,ncv,"c");
    
    StaticSolver::processStep(0);

    const int nsteps = getIntParam("NSTEPS",-1);

    for (int step = 0; step < nsteps; ++step) {
      
      COUT1("\n============================\nstep " << step << "\n============================");

      //FOR_ICV_G c[icv] *= 0.9;
      //FOR_ICV_G phi[icv] *= 0.9;
      
      doProbes(step);

      if (step%48 == 0)
	flushProbes();
    
      CtiRegister::clearCurrentData();

    }

    flushProbes();
    
    writeResult();

    if ( Param* param = getParam("GET_STREAMLINE")) { 

      // initial point defines streamline...
      const string name = param->getString(0);
      const double dx = param->getDouble(1);
      if (dx <= 0.0) {
        CERR("GET_STREAMLINE must have a postive displacement: GET_STREAMLINE <zone_name> <dx>");
      }

      bool found = false;
      double xp[3];
      FOR_IZONE(bfZoneVec) {
        if (name == bfZoneVec[izone].getName()) {
          FOR_I3 xp[i] = bfZoneVec[izone].x_global[i];
          found = true;
        }
      }
      if (!found) {
        CERR("GET_STREAMLINE must have a zone name: GET_STREAMLINE <zone_name> <dx>");
      }

      if (mpi_rank == 0) cout << "Starting position: " << COUT_VEC(xp) << ", dx: " << dx << endl;

      // need nodal velocity for interpolation to point...

      double (*u)[3] = new double[ncv][3];
      calcCvGrad(u,phi);
      double (*u_no)[3] = new double[nno][3];
      averageCvToNo(u_no,u);
      delete[] u;
      
      // build the adt so we can find point...

      if (cvAdt == NULL)
        buildCvAdt();
      
      // locate point... 
      
      vector<int> intVec;
      assert(intVec.empty());
      cvBboxAdt->buildListForPoint(intVec,xp);

      int icv_p = -1; // when > 0, lies on this rank
      if (intVec.size() > 0) {
        intVec.clear();
        cvAdt->buildListForPoint(intVec,xp);
        double dist = HUGE_VAL;
        for (int ii = 0; ii < intVec.size(); ++ii) {
          const int this_icv = intVec[ii]; assert((this_icv >= 0)&&(this_icv < ncv));
          const double this_dist = DIST(xp,x_cv[this_icv]); // x_vv?
          if (this_dist <= dist) {
            icv_p = this_icv;
            dist = this_dist;
          }
        }
      }
      
      // tell everyone your icv_p, so we can determine the owner rank...
      int *recv_buf = new int[mpi_size];
      MPI_Allgather(&icv_p,1,MPI_INT,recv_buf,1,MPI_INT,mpi_comm);
      int owner_rank = -1;
      FOR_RANK {
        // cout << recv_buf[rank] << " " << rank << endl;
        if (recv_buf[rank] >= 0) {
          owner_rank = rank;
          break;
        }
      }
      delete[] recv_buf;

      char filename[128];
      sprintf(filename,"streamline.dat");
      if (mpi_rank == owner_rank) {
        FILE * fp = fopen(filename,"w");
        fprintf(fp,"%18.16e %18.16e %18.16e\n",xp[0],xp[1],xp[2]);
        fclose(fp);
      }

      double bbox[6];
      getBoundingBox(bbox);

      int done = 0;
      while (done == 0) {

        if (mpi_rank == owner_rank) {
          // interp velocity to point...

          set<int> no_list; 
          double up[3] = {0.0,0.0,0.0};
          for (int foc = faocv_i[icv_p]; foc != faocv_i[icv_p+1]; ++foc) {
            const int ifa = faocv_v[foc];
            for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
              const int ino = noofa_v[nof];
              if (no_list.find(ino) == no_list.end()) {
                no_list.insert(ino);
                const double dist = DIST(x_no[ino],xp); 
                assert(dist > 0.0);
                const double wgt = 1.0/dist; // inverse distance weighted (normalized below)
                FOR_I3 up[i] += u_no[ino][i]*wgt;
              }
            }
          }
          const double mag = MAG(up);

          FOR_I3 xp[i] -= dx*up[i]/mag;

          // write out point...
          FILE * fp = fopen(filename,"a");
          fprintf(fp,"%18.16e %18.16e %18.16e\n",xp[0],xp[1],xp[2]);
          fclose(fp);

          // locate cell...
          
          FOR_I3 {
            if ((xp[i] < bbox[2*i]) || (xp[i] > bbox[2*i+1])) {
              done = 1;
              break;
            }
          }
          
          //cout << mpi_rank << " " << icv_p << "/" << ncv << " " << COUT_VEC(xp) << " " << done << " " <<  endl;
          if (done == 0) {
            assert((icv_p >= 0)&&(icv_p < ncv));
            double dist = HUGE_VAL;
            int new_icv = -1;
            for (int coc = cvocv_i[icv_p]; coc != cvocv_i[icv_p+1]; ++coc) {
              const int this_icv = cvocv_v[coc]; 
              const double this_dist = DIST(xp,x_cv[this_icv]); // x_vv?
              if (this_dist <= dist) {
                new_icv = this_icv;
                dist = this_dist;
              }
            }
            assert(new_icv >= 0);
            icv_p = new_icv;
          }
        }

        // if the cell is a ghost then send to nbr...
        
        int rank = -1,bits = -1,index = -1;
        if (mpi_rank == owner_rank) {
          if (icv_p > ncv) {
            BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_p-ncv]);
            assert(bits == 0); // can add later
            icv_p = -1;
          }
          else {
            rank = mpi_rank;
            index = icv_p;
          }
        }

        int buf[3] = {done,rank,index};
        MPI_Bcast(buf,3,MPI_INT,owner_rank,mpi_comm);

        done = buf[0];
        owner_rank = buf[1];
        if (mpi_rank == owner_rank)
          icv_p = buf[2];
      }
      
      delete[] u_no;

      if (mpi_rank == owner_rank) cout << "Final position: " << COUT_VEC(xp) << endl;

    }

    /*
      {
    
      double rhs_crs = 0;
      for (int ibf = 0; ibf < nbf; ++ibf) {
	const int icv = cvobf[ibf];
	rhs_crs += rhoun_bf[ibf];
      }
      
      for (int ifa = nfa_i; ifa < nfa; ++ifa) {
	const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
	const int icv1 = cvofa[ifa][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g));
	const double coeff = fa[ifa].area/DIST(x_cv[icv0],x_cv[icv1]);
	rhs_crs -= coeff*(phi[icv1]-phi[icv0]);
      }
      
      //cout << "rank: " << mpi_rank << " rhs_crs after full: " << rhs_crs << endl;

    }
    */

    /*
    Histogram hist(phi,ncv,200);
    hist.write("hist.dat");
    
    writeImages(0);
    
    doProbes(0);
    */

  }

  /*
  void writeCvDataCheckpoint(const double * const var,const int n,const string& prefix) const {
    
    char filename[256];
    sprintf(filename,"%s-%04d.bin",prefix.c_str(),mpi_rank);
    FILE * fp = fopen(filename,"wb");
    assert(fp);
    fwrite(&n,sizeof(int),1,fp);
    fwrite(var,sizeof(double),n,fp);
    fclose(fp);

  }
  
  void readCvDataCheckpoint(double * var,const int n,const string& prefix) const {
    
    char filename[256];
    sprintf(filename,"%s-%04d.bin",prefix.c_str(),mpi_rank);
    FILE * fp = fopen(filename,"rb");
    assert(fp);
    int n_check;
    fread(&n_check,sizeof(int),1,fp);
    assert(n_check == n);
    fread(var,sizeof(double),n,fp);
    fclose(fp);
    
  }
  */

  /*
  void writeImages(const int step) {

    FOR_PARAM_MATCHING("WRITE_IMAGE") {

      processWriteImage(&(*param),step);
      
    }
    
  }
  */
  
  void processWriteImage(Param * param,const int step) {

    // once again, needs re-plumbing for imaging...
    
    assert(0);
    
#ifdef DKDASLKDLKAS

    WriteImageData wid(param);
    if ((step%wid.interval == 0)&&(wid.got_var)) {
	
      int ierr = 0;
      if (wid.got_frac) {
	if (mpi_rank == 0)
	  cout << " > fractional geometry not supported yet" << endl;
	ierr = -1;
      }
      if (!wid.got_width) {
	if (mpi_rank == 0)
	  cout << " > you mut specify WIDTH param for now" << endl;
	ierr = -1;
      }
      if (!wid.got_geom) {
	if (mpi_rank == 0)
	  cout << " > you mut specify valid GEOM param for now" << endl;
	ierr = -1;
      }
      if (ierr != 0)
	return;
      
      CtiCanvas canvas(wid.xp,wid.np,wid.width,wid.nx,wid.ny);
	
      if (wid.got_range) 
	canvas.setRange(wid.rmin,wid.rmax);

      // ---------------------------------------------------------------------------------
      // figure out the fraction of tris we should render. In the future, we will have
      // a subset of the tris anyways as part of SubSurface...
      // ---------------------------------------------------------------------------------
      
      int * stora = NULL;
      buildUniformXora(stora,surface->nst);
      canvas.addSurfaceTris(stora[mpi_rank+1]-stora[mpi_rank],surface->spost+stora[mpi_rank],surface->xp);
      delete[] stora;
      
      canvas.activateScalarPlane();

      if (wid.var == "PHI")
	canvas.addScalarPlane(phi,ncv,x_cv,r_cv);
      else if (wid.var == "C")
	canvas.addScalarPlane(c,ncv,x_cv,r_cv);
      else {
	if (mpi_rank == 0)
	  cout << "Skipping unrecognized VAR: " << wid.var << endl;
      }
      
      canvas.writeImage(wid.name,step);
      
    }

#endif

  }
  
  void buildCvLaplacian(double * A) {

    assert(cvocv_i);
    assert(cvocv_v);
    
    for (int coc = 0; coc < cvocv_i[ncv]; ++coc)
      A[coc] = 0.0;

    for (int ifa = 0; ifa < nfa_i; ++ifa) {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
      const double coeff = area_over_delta_fa[ifa];
      {
	int coc = cvocv_i[icv0];
	assert(cvocv_v[coc] == icv0);
	A[coc] -= coeff;
	++coc;
	while (cvocv_v[coc] != icv1)
	  ++coc;
	assert(coc < cvocv_i[icv0+1]);
	A[coc] += coeff;
      }
      {
	int coc = cvocv_i[icv1];
	assert(cvocv_v[coc] == icv1);
	A[coc] -= coeff;
	++coc;
	while (cvocv_v[coc] != icv0)
	  ++coc;
	assert(coc < cvocv_i[icv1+1]);
	A[coc] += coeff;
      }
    }

    for (int ifa = nfa_i; ifa < nfa; ++ifa) {
      const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
      const int icv1 = cvofa[ifa][1]; assert((icv1 >= ncv)&&(icv1 < ncv_g));
      const double coeff = area_over_delta_fa[ifa];
      {
	int coc = cvocv_i[icv0];
	assert(cvocv_v[coc] == icv0);
	A[coc] -= coeff;
	++coc;
	while (cvocv_v[coc] != icv1)
	  ++coc;
	assert(coc < cvocv_i[icv0+1]);
	A[coc] += coeff;
      }
    }
  
  }
  
  int solveCvCg(double * phi,const double * const A,const double * const rhs) {

    bool verbose = true;
    int maxiter = 10000;
    double zero = getDoubleParam("PHI_ZERO",1.0E-06);
    
    // assume we come in with a consistent initial condition...

    // we need the following work arrays...
    
    double * res      = new double[ncv];
    double * v        = new double[ncv];
    double * p        = new double[ncv_g];
    double * inv_diag = new double[ncv];
    
    // initialize...
    for (int icv = 0; icv < ncv; ++icv) 
      inv_diag[icv] = 1.0/A[cvocv_i[icv]];

    for (int icv = 0; icv < ncv; ++icv) 
      p[icv] = 0.0;
    double rho = 1.0;
    
    // calculate the residual in rhs format...   
    for (int icv = 0; icv < ncv; ++icv) {
      res[icv] = rhs[icv] - A[cvocv_i[icv]]*phi[icv];
      for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
	const int icv_nbr = cvocv_v[coc];
	res[icv] -= A[coc]*phi[icv_nbr];
      }
    }
    
    // diagonal precon/compute normalized residual...
    for (int icv = 0; icv < ncv; ++icv)
      v[icv] = res[icv]*inv_diag[icv];

    int iter = 0;
    int done = 0;
    while (done == 0) {
    
      ++iter;
      
      double rho_prev = rho;
      if (fabs(rho_prev) < 1.0E-20)
	rho_prev = -1.0E-20; // -1.0E-20? seems to help
      
      double my_rho = 0.0; 
      for (int icv = 0; icv < ncv; ++icv)
	my_rho += res[icv]*v[icv];
      MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
      
      double beta = rho/rho_prev;
      for (int icv = 0; icv < ncv; ++icv)
	p[icv] = v[icv] + beta*p[icv];
      updateCvData(p);
      
      // v = [Ap]{p}...
      for (int icv = 0; icv < ncv; ++icv) {
	v[icv] = A[cvocv_i[icv]]*p[icv];
	for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
	  const int icv_nbr = cvocv_v[coc];
	  v[icv] += A[coc]*p[icv_nbr];
	}
      }
      
      double my_gamma = 0.0; 
      for (int icv = 0; icv < ncv; ++icv)
	my_gamma += p[icv]*v[icv];
      double gamma; 
      MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
      if (fabs(gamma) < 1.0E-20)
	gamma = 1.0E-20;
      
      const double alpha = rho/gamma;
      
      // check if we are done...
      if (iter%3 == 0) {
	
	for (int icv = 0; icv < ncv; ++icv) 
	  phi[icv] += alpha*p[icv];
	updateCvData(phi);
	
	// recompute the residual...
	for (int icv = 0; icv < ncv; ++icv) {
	  res[icv] = rhs[icv] - A[cvocv_i[icv]]*phi[icv];
	  for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
	    const int icv_nbr = cvocv_v[coc];
	    res[icv] -= A[coc]*phi[icv_nbr];
	  }
	}

	for (int icv = 0; icv < ncv; ++icv) 
	  v[icv] = res[icv]*inv_diag[icv];
	
	// compute the max (L-infinity) normalized residual...
	double  my_res_max = 0.0;        
	for (int icv = 0; icv < ncv; ++icv) 
	  my_res_max = max( my_res_max, fabs(v[icv]) );
	double res_max;
	MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
	if (mpi_rank == 0) {
	  // only share the last half of the convergence behaviour...
	  if (verbose || (iter > maxiter/2))
	    cout << " > solveCvCg iter " << iter << " res_max " << res_max << endl;
	  if (res_max <= zero) {
	    cout << "-> Successfully converged error to " << res_max << endl;
	    done = 1;
	  }
	  else if (iter > maxiter) {
	    cout << "Warning: solveCvCg did not converge after " << maxiter << 
	      " iters, res_max: " << res_max << endl;
	    done = 2;
	  }
	}
	MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
        
      }
      else {
	
	// update full phi including ghosts...
	for (int icv = 0; icv < ncv_g; ++icv)
	  phi[icv] += alpha*p[icv];
	
	for (int icv = 0; icv < ncv; ++icv) {
	  // on the other iterations, use this approximation to update
	  // the unreduced residual...
	  res[icv] -= alpha*v[icv];
	  // still need to compute v, diag precon for next iteration...
	  v[icv] = res[icv]*inv_diag[icv];
	}
	
      }
      
    }
    
    delete[] res;
    delete[] v;
    delete[] p;
    delete[] inv_diag;
    
    // let the calling routine know if we were successful...
    return( done == 1 );
    
  }
  
  int solveCgLocal(double * phi,const int ncv,const double * const A,const double * const rhs,
		   const int * const cvocv_i,const int * const cvocv_v) {
    
    bool verbose = true;
    int maxiter = 10000;
    double zero = 1.0E-08;
    
    // assume we come in with a consistent initial condition...

    // we need the following work arrays...
    
    double * res      = new double[ncv];
    double * v        = new double[ncv];
    double * p        = new double[ncv];
    double * inv_diag = new double[ncv];
    
    // initialize...
    for (int icv = 0; icv < ncv; ++icv)  {
      if (A[cvocv_i[icv]] == 0.0) {
	// must be just one detached unknown...
	assert(cvocv_i[icv+1]-cvocv_i[icv] == 1);
	assert(rhs[icv] == 0.0);
	inv_diag[icv] = 0.0;
      }
      else {
	inv_diag[icv] = 1.0/A[cvocv_i[icv]];
      }
    }

    for (int icv = 0; icv < ncv; ++icv) 
      p[icv] = 0.0;
    double rho = 1.0;
    
    // calculate the residual in rhs format...   
    for (int icv = 0; icv < ncv; ++icv) {
      res[icv] = rhs[icv] - A[cvocv_i[icv]]*phi[icv];
      for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
	const int icv_nbr = cvocv_v[coc];
	res[icv] -= A[coc]*phi[icv_nbr];
      }
    }
    
    // diagonal precon/compute normalized residual...
    for (int icv = 0; icv < ncv; ++icv)
      v[icv] = res[icv]*inv_diag[icv];

    int iter = 0;
    int done = 0;
    while (done == 0) {
      
      ++iter;
      
      double rho_prev = rho;
      if (fabs(rho_prev) < 1.0E-20)
	rho_prev = -1.0E-20; // -1.0E-20? seems to help
      
      double rho = 0.0; 
      for (int icv = 0; icv < ncv; ++icv)
	rho += res[icv]*v[icv];
      
      double beta = rho/rho_prev;
      for (int icv = 0; icv < ncv; ++icv)
	p[icv] = v[icv] + beta*p[icv];
      
      // v = [Ap]{p}...
      for (int icv = 0; icv < ncv; ++icv) {
	v[icv] = A[cvocv_i[icv]]*p[icv];
	for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
	  const int icv_nbr = cvocv_v[coc];
	  v[icv] += A[coc]*p[icv_nbr];
	}
      }
      
      double gamma = 0.0; 
      for (int icv = 0; icv < ncv; ++icv)
	gamma += p[icv]*v[icv];
      if (fabs(gamma) < 1.0E-20)
	gamma = 1.0E-20;
      
      const double alpha = rho/gamma;
      
      // check if we are done...
      if (iter%3 == 0) {
	
	for (int icv = 0; icv < ncv; ++icv) 
	  phi[icv] += alpha*p[icv];
	
	// recompute the residual...
	for (int icv = 0; icv < ncv; ++icv) {
	  res[icv] = rhs[icv] - A[cvocv_i[icv]]*phi[icv];
	  for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
	    const int icv_nbr = cvocv_v[coc];
	    res[icv] -= A[coc]*phi[icv_nbr];
	  }
	}
	
	for (int icv = 0; icv < ncv; ++icv) 
	  v[icv] = res[icv]*inv_diag[icv];
	
	// compute the max (L-infinity) normalized residual...
	double  res_max = 0.0;        
	for (int icv = 0; icv < ncv; ++icv) 
	  res_max = max( res_max, fabs(v[icv]) );
	// only share the last half of the convergence behaviour...
	if (verbose || (iter > maxiter/2))
	  cout << " > solveCgLocal iter " << iter << " res_max " << res_max << endl;
	if (res_max <= zero) {
	  cout << "-> Successfully converged error to " << res_max << endl;
	  done = 1;
	}
	else if (iter > maxiter) {
	  cout << "Warning: solveCgLocal did not converge after " << maxiter << 
	    " iters, res_max: " << res_max << endl;
	  done = 2;
	}
        
      }
      else {
	
	// update full phi...
	for (int icv = 0; icv < ncv; ++icv)
	  phi[icv] += alpha*p[icv];
	
	for (int icv = 0; icv < ncv; ++icv) {
	  // on the other iterations, use this approximation to update
	  // the unreduced residual...
	  res[icv] -= alpha*v[icv];
	  // still need to compute v, diag precon for next iteration...
	  v[icv] = res[icv]*inv_diag[icv];
	}
	
      }
      
    }
    
    delete[] res;
    delete[] v;
    delete[] p;
    delete[] inv_diag;
    
    // let the calling routine know if we were successful...
    return( done == 1 );
    
  }

  virtual CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,const string& name,list<CtiRegister::CtiData>& args,const bool b_eval_func) {
    
    //if (mpi_rank == 0) 
    //  cout << "PotentialSolver::funcEvalCtiData() being asked to evaluate: \"" << name << "\"" << endl;
    
    if (name == "mdot_fa") {
      
      if (!args.empty()) 
	return(CtiRegister::CTI_DATA_ARG_COUNT);

      double * v_ptr = createSignedFaD1Data(v);
      
      if (b_eval_func) {
        FOR_IFA {
          v_ptr[ifa] = area_over_delta_fa[ifa]*(phi[cvofa[ifa][0]]-phi[cvofa[ifa][1]]);
        }
      }
      
      return CtiRegister::CTI_DATA_OK;
      
    }
    else if (name == "mdotc_fa") {
      
      if (!args.empty()) 
	return(CtiRegister::CTI_DATA_ARG_COUNT);

      double * v_ptr = createSignedFaD1Data(v);
      
      if (b_eval_func) {
        FOR_IFA {
          const double this_mdot = area_over_delta_fa[ifa]*(phi[cvofa[ifa][0]]-phi[cvofa[ifa][1]]);
          if (this_mdot >= 0.0) {
            v_ptr[ifa] = this_mdot*c[cvofa[ifa][0]];
          }
          else {
            v_ptr[ifa] = this_mdot*c[cvofa[ifa][1]];
          }
        }
      }
      
      return CtiRegister::CTI_DATA_OK;
      
    }
    
    return StaticSolver::funcEvalCtiData(v,name,args,b_eval_func);
    
  }




};

#endif
