#include "StaticSolver.hpp"
#include "LpInjector.hpp"
#include "MiscUtils.hpp"
#include "SubSurface.hpp"

double error_function2(const double eta);

LpInjectorTracer::LpInjectorTracer() {
  IF_RANK0 cout << "LpInjectorTracer()" << endl;

  InjectorKind = -1;
  name = "";
  npdot = 0.0;
  np = 0;
  dp = -1;

  GeomKind = -1;
  for (int i = 0; i < 11; ++i) GeomData[i] = 0.0;

  ResidualCount = 0.0;

  time0 = -1.0;
  injectDt = -1.0;

  nst = 0;
  nsp = 0;
  spost = NULL;
  xp = NULL;
  cdf_area = NULL;

  FOR_I3 {
    e1[i] = 0;
    e2[i] = 0;
  }

  triVec.clear();
  cdf_area_tri = NULL;
  cdf_area_rank = NULL;

  solver = NULL;

  injector_id = -1;

}

LpInjectorTracer::~LpInjectorTracer() {
  IF_RANK0 cout << "~LpInjectorTracer()" << endl;
  
  if (spost != NULL) delete[] spost;
  if (xp != NULL) delete[] xp;
  if (cdf_area != NULL) delete[] cdf_area;
  if (cdf_area_tri != NULL) delete[] cdf_area_tri;
  if (cdf_area_rank != NULL) delete[] cdf_area_rank;
}

void LpInjectorTracer::initFromParams(Param * param, StaticSolver * solver, double time0, int injector_id) {

  this->solver = solver;
  this->time0 = time0;
  this->injector_id = injector_id;
  assert(injector_id >= 0);

  string token = param->getString();
  if (token != "NAME") CERR("LP_TRACER.INJECTOR needs NAME. Unrecogniced token "<<token<<"\n");
  string name = param->getString(1);
  this->name = name;

  if (mpi_rank == 0)
    cout << " > initializing LP_TRACER.INJECTOR: " << name << endl;

  int iarg = 2;
  while (iarg < param->size()) {
    string token = param->getString(iarg++);
    if (token == "TYPE") {
      // ==========================================
      // TYPE of injector:
      // ==========================================
      token = param->getString(iarg++);
      if (token == "COUNT") {
        assert( InjectorKind == -1 );
        InjectorKind = LpInjectorNS::COUNT_INJ;
        np = param->getInt(iarg++);
        if (mpi_rank == 0)
          cout << "   > injector TYPE = COUNT " << np << endl;
      }
      else if ( token == "NDOT" ) {
        assert( InjectorKind == -1 );
        InjectorKind = LpInjectorNS::COUNTDOT_INJ;
        npdot =  param->getDouble(iarg++);
        if (mpi_rank == 0)
          cout << "   > injector TYPE = NDOT " << npdot << endl;
      }
      else if ( token == "STRUCTURED" ) {
        assert( InjectorKind == -1 );
        InjectorKind = LpInjectorNS::STRUCTURED_INJ;
        this->npx = param->getInt(iarg++);
        this->npy = param->getInt(iarg++);
        this->npz = param->getInt(iarg++);
        this->np = npx * npy * npz;
        if (mpi_rank == 0)
          cout << "   > injector TYPE = STRUCTURED npx: " << npx << " npy: " << npy << " npz: " << npz << endl;
      }
      else if (token == "EACH_CELL") {
        assert( InjectorKind == -1 );
        InjectorKind = LpInjectorNS::EACHCELL_INJ;
        this->np_in_cell = param->getInt(iarg++); // number of particles in each cell
        int np_local = solver->ncv * np_in_cell;
        MPI_Reduce(&np_local,&np,1,MPI_INT,MPI_SUM,0,mpi_comm);
      }
      else {
        CERR( "\n\n\n*********************************************************\n" <<
            "Error: unrecognized TYPE in LP.INJECTOR: " << token << ". Valid types are: \n" <<
            "TYPE COUNT <n>                     -- injects n particles once\n" <<
            "TYPE NDOT <npdot>                  -- injects at rate npdot\n" <<
            "TYPE STRUCTURED <npx> <npy> <npz>  -- injects structured points in axis-aligned manner\n" <<
            "TYPE EACH_CELL <np_in_cell>        -- injects at each cell np_in_cell tracers\n" <<
            "*********************************************************\n\n\n");
      }
    }
    else if (token == "GEOM") {
      // ==========================================
      // GEOM of particles:
      // ==========================================
      token = param->getString(iarg++);
      if (token == "POINT") {
        assert( GeomKind == -1 );
        GeomKind = LpInjectorNS::POINT_GEOMETRY;
        GeomData[0] = param->getDouble(iarg++); // <x>
        GeomData[1] = param->getDouble(iarg++); // <y>
        GeomData[2] = param->getDouble(iarg++); // <z>
        if (mpi_rank == 0)
          cout << "   > injector GEOM = POINT  X=" << GeomData[0] << " " << GeomData[1] << " " << GeomData[2] << endl;
      }
      else if (token == "BOX") {
        assert( GeomKind == -1 );
        GeomKind = LpInjectorNS::BOX_GEOMETRY;
        GeomData[0] = param->getDouble(iarg++); // <xmin>
        GeomData[1] = param->getDouble(iarg++); // <xmax>
        GeomData[2] = param->getDouble(iarg++); // <ymin>
        GeomData[3] = param->getDouble(iarg++); // <ymax>
        GeomData[4] = param->getDouble(iarg++); // <zmin>
        GeomData[5] = param->getDouble(iarg++); // <zmax>
        if (mpi_rank == 0)
          cout << "   > injector GEOM = BOX X=" << GeomData[0] << ":" << GeomData[1] << 
            " Y=" << GeomData[2] << ":" << GeomData[3] <<
            " Z=" << GeomData[4] << ":" << GeomData[5] << endl; 
        if ( (GeomData[0] > GeomData[1]) or (GeomData[2] > GeomData[3]) or (GeomData[4] > GeomData[5]) ) 
          CERR("BOX is defined incorrectly in LSP.STATIC_INJECTOR. Check box dimensions...");
      } 
      else if (token == "CIRCLE") {
        assert( GeomKind == -1 );
        GeomKind = LpInjectorNS::CIRCLE_GEOMETRY;
        GeomData[0] = param->getDouble(iarg++); // <x>
        GeomData[1] = param->getDouble(iarg++); // <y>
        GeomData[2] = param->getDouble(iarg++); // <z>

        GeomData[3] = param->getDouble(iarg++);
        GeomData[4] = param->getDouble(iarg++);
        GeomData[5] = param->getDouble(iarg++);
        NORMALIZE((GeomData+3));

        GeomData[6] = param->getDouble(iarg++);
        if (mpi_rank == 0)
          cout << "   > injector GEOM = CIRCLE X=" << GeomData[0] << " " << GeomData[1] << " " << GeomData[2] << " ,N= " << GeomData[3] << " " << GeomData[4] << " " << GeomData[5] << " ,R= " << GeomData[6] << endl; 
      }
      else if (token == "ZONE") {
        assert( GeomKind == -1 );
        GeomKind = LpInjectorNS::ZONE_GEOMETRY;

        string zone_name = param->getString(iarg++);

        int ierr = solver->initZoneInjector(this->cdf_area, this->spost, this->xp, this->nst, this->nsp, zone_name);
        if (ierr == 1)
          CERR("Problem in ZONE injector initiation, check the zone name")

      }
      else if (token == "PLANE") {
        assert( GeomKind == -1 );
        GeomKind = LpInjectorNS::PLANE_GEOMETRY;

        GeomData[0] = param->getDouble(iarg++);
        GeomData[1] = param->getDouble(iarg++);
        GeomData[2] = param->getDouble(iarg++);
        
        GeomData[3] = param->getDouble(iarg++);
        GeomData[4] = param->getDouble(iarg++);
        GeomData[5] = param->getDouble(iarg++);
        NORMALIZE((GeomData+3));

        int ierr = solver->initPlaneInjector(cdf_area_rank, cdf_area_tri, triVec, GeomData, (GeomData+3));
        if (ierr == 1) 
          CERR("Problem in PLANE injector initiation, check the plane location and normal");

      }
      else {
        CERR( "Error: unrecognized GEOM in LP.INJECTOR: " << token << ". Valid geoms are: \n" <<
            "GEOM POINT <x> <y> <z>\n" <<
            "GEOM BOX <xmin> <xmax> <ymin> <ymax> <zmin> <zmax>\n" <<
            "GEOM CIRCLE <x> <y> <z> <nx> <ny> <nz> <r>\n" <<
            "GEOM ZONE <zone-name>\n" <<
            "GEOM PLANE <x> <y> <z> <nx> <ny> <nz>");
      }
    }
    else if (token == "TIME_SPAN") {
      injectDt = param->getDouble(iarg++);
      if (injectDt <= 0) CERR("TIME_SPAN should be positive.");
      if (mpi_rank == 0)
        cout << "   > injector time span:  " << injectDt << endl;
    } else if (token == "DP") {
      dp = param->getDouble(iarg++);
      if (dp <=0 ) CERR("DP must be positive.");
    }
    else {
      CERR( "Error: unrecognized keyword in LP.INJECTOR: " << token);
    }
  }

  // now check...
  if ( InjectorKind == -1 ) {
     CERR( "Error: unrecognized TYPE in LP.INJECTOR: " << token << ". Valid types are: \n" <<
         "TYPE COUNT <n>                     -- injects n particles once\n" <<
         "TYPE NDOT <npdot>                  -- injects at rate npdot\n" <<
         "TYPE STRUCTURED <npx> <npy> <npz>  -- injects structured points in axis-aligned manner\n" <<
         "TYPE EACH_CELL <np_in_cell>        -- injects at each cell np_in_cell tracers");
  }


  if ( GeomKind == -1 ) {
    CERR( "Error: LP.INJECTOR: " << name << " did not specified GEOM: Valid geoms are: \n" <<
        "GEOM POINT <x> <y> <z>\n" <<
        "GEOM BOX <xmin> <xmax> <ymin> <ymax> <zmin> <zmax>\n" <<
        "GEOM CIRCLE <x> <y> <z> <nx> <ny> <nz> <r>\n" <<
        "GEOM ZONE <zone-name>\n" <<
        "GEOM PLANE <x> <y> <z> <nx> <ny> <nz>");
  }

  if ( dp == -1)
    CERR("Missing DP from LP_TRACER.INJECTOR");

}

void LpInjectorTracer::addNewLp(vector<LpDataBase>& lpDataVec, const double time, const double dt, const bool verbose) {

  // time: the relative time from starting the simulation
  // Don't add particles if larger than the injection time
  // for one-time injection we set injectDt to zero after injection
  if ((injectDt >= 0) && ((time-time0) > injectDt)) return;

  //int from_ind = lpDataVec.size();

  // for STRUCTURED_INJ/EACHCELL_INJ one time injection is intended
  if ( (InjectorKind == LpInjectorNS::STRUCTURED_INJ) || (InjectorKind == LpInjectorNS::EACHCELL_INJ) ) {

    if (GeomKind == LpInjectorNS::BOX_GEOMETRY) {

      if ( InjectorKind == LpInjectorNS::STRUCTURED_INJ ) {

        //int my_np_count = solver->addStructuredParticlesInBox(lpDataVec, GeomData, npx, npy, npz, dp, name);
        solver->addStructuredParticlesInBox(lpDataVec, GeomData, npx, npy, npz, dp, name);
      
      } else if ( InjectorKind == LpInjectorNS::EACHCELL_INJ ) {

        //int my_np_count = solver->addEachCellParticlesInBox(lpDataVec, GeomData, np_in_cell, dp, name);
        solver->addEachCellParticlesInBox(lpDataVec, GeomData, np_in_cell, dp, name);

      }
  
    }
    else if (GeomKind == LpInjectorNS::CIRCLE_GEOMETRY) {
      CERR("CIRCLE GEOMETRY STILL NOT IMPLEMENTED");
    }
    else {
      CERR("STRUCTURED or EACHCELL injectors are not designed for this geometry...");
    }

    injectDt = 0;
    return;
  }


  int npnew = 0;                         // number of particles injected at each time step
  double (*xpnew)[3] = NULL;             // particle position
  int * np_count_rank = NULL;            // np for each rank for plane injectors

  // ============================================
  // on rank 0 build the injected particle list:
  // This is necessary because random numbers are involved
  // in some of the particle distributions, and we want
  // all ranks to start with the same distribution.
  // ============================================

  IF_RANK0 {

    // if the injector is NDOT modify the residual count

    if ( InjectorKind == LpInjectorNS::COUNTDOT_INJ )
      ResidualCount += npdot*dt;

    // -----------------------------
    // number of particles...
    // -----------------------------

    if ( InjectorKind == LpInjectorNS::COUNTDOT_INJ ){
      double npcheck = ResidualCount;
      if (npcheck > 1000000.0) 
        CERR( "Error: LP Tracer Injector " << getName() << " is trying to inject " << npcheck << " particles this time step. Check LP_TRACER.INJECT settings.");
      
      npnew = (int)npcheck;
      ResidualCount -= double(npnew);
    }
    else if (InjectorKind == LpInjectorNS::COUNT_INJ) {
      // for these injectors, the np stores the count. Use it once, then
      // reset to zero...
      npnew = np;
      np = 0;
    }
    else{
      assert(0);
    }

    // -----------------------
    // then geometry...
    // -----------------------

    if ( npnew > 0 ) {

      xpnew = new double[npnew][3];

      switch (GeomKind) {

      case LpInjectorNS::POINT_GEOMETRY:
      {
        for (int j = 0; j < npnew; ++j) {
          FOR_I3 xpnew[j][i] = GeomData[i];
        }
        break;
      }
      case LpInjectorNS::BOX_GEOMETRY:
      {
        for (int j = 0; j < npnew; ++j) {
          double rand_x[3];
          FOR_I3 rand_x[i] = double(rand())/double(RAND_MAX);
          FOR_I3 xpnew[j][i] = GeomData[2*i] + (GeomData[2*i+1] - GeomData[2*i])*rand_x[i];
        }
      }
      break;

      case LpInjectorNS::CIRCLE_GEOMETRY:
      {
        double X[3], n[3], R;
        FOR_I3 X[i] = GeomData[i]; 
        FOR_I3 n[i] = GeomData[i+3];
        R = GeomData[6];
        int i_min = -1;
        double n_min = 1e20;
        FOR_I3 {
          if (abs(n[i]) < n_min) {
            i_min = i;
            n_min = n[i];
          }
        }
        double dir[3];
        dir[i_min] = 0.0;
        if (i_min == 0) {
          dir[1] = n[2];
          dir[2] = -n[1];
        } else if (i_min == 1) {
          dir[0] = n[2];
          dir[2] = -n[0];
        } else if (i_min == 2) {
          dir[0] = n[1];
          dir[1] = -n[0];
        } else {assert(0);}
        double mag_dir = MAG(dir);
        assert(mag_dir>0);
        FOR_I3 dir[i] /= mag_dir;
        assert(DOT_PRODUCT(dir,n) < 1e-12);
        for (int j = 0; j < npnew; ++j) {
          double rand_r_ = rand()/double(RAND_MAX);
          double rand_t = rand()/double(RAND_MAX);   // random variables for r and theta
          double r = R * sqrt(rand_r_);
          double theta = rand_t * 2. * M_PI;
          double Rot[3][3];
          // build three dimensional rotation matrix
          FOR_I3 Rot[i][i] = cos(theta) + n[i]*n[i]*(1.-cos(theta));
          Rot[0][1] = n[0]*n[1]*(1.-cos(theta)) - n[2]*sin(theta);
          Rot[1][0] = n[0]*n[1]*(1.-cos(theta)) + n[2]*sin(theta);
          Rot[0][2] = n[0]*n[2]*(1.-cos(theta)) + n[1]*sin(theta);
          Rot[2][0] = n[0]*n[2]*(1.-cos(theta)) - n[1]*sin(theta);
          Rot[2][1] = n[1]*n[2]*(1.-cos(theta)) + n[0]*sin(theta);
          Rot[1][2] = n[1]*n[2]*(1.-cos(theta)) - n[0]*sin(theta);
          double rand_dir[3] = MATRIX_PRODUCT(Rot,dir);
          assert(MAG(rand_dir) - 1. < 1e-12);
          assert(DOT_PRODUCT(rand_dir,n) < 1e-12);
          FOR_I3 xpnew[j][i] = X[i] + r * rand_dir[i];
        }
      }
      break;

      case LpInjectorNS::ZONE_GEOMETRY:
      {
        assert(this->nst>0);
        assert(this->nsp>0);
        assert(this->spost);
        assert(this->xp);
        assert(this->cdf_area);
        for (int j = 0; j < npnew; ++j) {
          // first randomly choose a tri then choose the point in the tri randomly
          double rand_val[3];
          FOR_I3 rand_val[i] = rand()/double(RAND_MAX);
          int ist_rand = 0;
          for (int ist = 1; ist < this->nst+1; ist++) {
            if (rand_val[0] < cdf_area[ist]) {
              ist_rand = ist-1;
              break;
            }
          }
          //IF_RANK0 cout << "rand_val: " << COUT_VEC(rand_val) << " ist_rand: " << ist_rand << endl;
          assert((ist_rand>=0)&&(ist_rand<this->nst));

          const double * const x0 = xp[spost[ist_rand][0]];
          const double * const x1 = xp[spost[ist_rand][1]];
          const double * const x2 = xp[spost[ist_rand][2]];

          // set the normal to be inward
          double normal[3] = TRI_NORMAL_2(x0,x2,x1);
          double mag = MAG(normal);
          FOR_I3 normal[i] /= mag;
          
          // locate a random point inside the tri by partitioning 1 into three segments
          double rand_min = min(rand_val[1],rand_val[2]);
          double rand_max = max(rand_val[1],rand_val[2]);
          // move the particle slightly inside
          double small = 1e-7;
          FOR_I3 xpnew[j][i] = rand_min*x0[i] + (rand_max-rand_min)*x1[i] + (1. - rand_max)*x2[i] + small*normal[i];
        }
      }
      break;

      case LpInjectorNS::PLANE_GEOMETRY:
      {
        assert(cdf_area_rank);
        assert(np_count_rank==NULL);
        np_count_rank = new int [mpi_size];
        FOR_RANK np_count_rank[rank] = 0;

        for (int ip = 0; ip < npnew; ++ip) {
          double rand_rank = double(rand())/(double(RAND_MAX) + 1.0);
          FOR_RANK {
            if (rand_rank < cdf_area_rank[rank]) {
              assert(rank < mpi_size);
              np_count_rank[rank]++;
              break;
            }
          }
        }

        int np_sum = 0;
        FOR_RANK np_sum += np_count_rank[rank];
        assert(np_sum == npnew);
      }
      break;

      default:
        assert(0);
      }

    }

  }

  // ======================================
  // rank 0 bcast to everyone...
  // every rank should deside which particles
  // should be injected in their rank.
  // for plane injector each rank should 
  // determine the location of injection as well.
  // ======================================

  MPI_Bcast(&npnew,1,MPI_INT,0,mpi_comm);
  if (npnew == 0) return;

  switch (GeomKind) {

    case LpInjectorNS::PLANE_GEOMETRY:    // in case of a plane injector the location of the points are not still determined
                          // so each rank prepares particles for itself
    {

    //const double * const xp_plane = GeomData;
    //const double * const np_plane = GeomData+3;

    assert(cdf_area_tri);
    IF_RANK0 {
      assert(np_count_rank);
      assert(cdf_area_rank);
    }

    int my_npnew;
    MPI_Scatter(np_count_rank, 1, MPI_INT, &my_npnew, 1, MPI_INT, 0, mpi_comm);
    
    int npold = lpDataVec.size();
    lpDataVec.resize(npold+my_npnew);

    int ipnew = 0;
    while (ipnew < my_npnew) {
      
      int itri = -1;
      double rand_tri = double(rand()) / (double(RAND_MAX) + 1.0);
      for (int ii = 0; ii < triVec.size(); ++ii) {
        if (rand_tri < cdf_area_tri[ii]) {
          itri = ii;
          break;
        }
      }
      assert((itri>=0)&&(itri<triVec.size()));

      // locate a random point inside the tri by partitioning 1 into three segments
      double rand_val[2]; 
      FOR_I2 rand_val[i] = double(rand()) / (double(RAND_MAX) + 1.0);
      double rand_min = min(rand_val[0],rand_val[1]);
      double rand_max = max(rand_val[0],rand_val[1]);
  
      const double * const x0 = triVec[itri].first.x0;
      const double * const x1 = triVec[itri].first.x1;
      const double * const x2 = triVec[itri].first.x2;

      double xp_seed[3];
      FOR_I3 xp_seed[i] = rand_min*x0[i] + (rand_max-rand_min)*x1[i] + (1. - rand_max)*x2[i];

      int icv_closest;
      const int isInside = solver->pointIsInside(xp_seed, icv_closest);

      if ( isInside ) {
        assert( (icv_closest >= 0) && (icv_closest < solver->ncv));
        lpDataVec[npold+ipnew].icv = icv_closest;
        short int * material_injector_flag = (short int *)&lpDataVec[npold+ipnew].flag;
        material_injector_flag[0] = 0;
        material_injector_flag[1] = injector_id;
        FOR_I3 lpDataVec[npold+ipnew].xp[i] = xp_seed[i];
        lpDataVec[npold+ipnew].dp = dp;
        ipnew++;
      }
    }

    int np_count;
    MPI_Allreduce(&ipnew, &np_count, 1, MPI_INT, MPI_SUM, mpi_comm);
    assert(np_count == npnew);

    delete[] np_count_rank;
    }
    break;

    default:      // in the default case rank 0 sets all particles' xp and we just decide which rank owns which particles
    {
    if (mpi_rank != 0) {
      xpnew    = new double[npnew][3];
    }

    MPI_Bcast((double*)xpnew, 3*npnew,MPI_DOUBLE,0,mpi_comm);

    // add to the lpDataVec...

    int my_np_count = 0;
    for (int ip = 0; ip < npnew; ++ip) {
      const double xp_seed[3] = {xpnew[ip][0], xpnew[ip][1], xpnew[ip][2]};
      int icv_closest;
      const int isInside = solver->pointIsInside(xp_seed,icv_closest);
      if ( isInside ) {
        assert( (icv_closest >= 0) && (icv_closest < solver->ncv));
        lpDataVec.resize(lpDataVec.size()+1);
        int iback = lpDataVec.size()-1;
        lpDataVec[iback].icv = icv_closest;
        short int * material_injector_flag = (short int *)&lpDataVec[iback].flag;
        material_injector_flag[0] = 0;
        material_injector_flag[1] = injector_id;
        FOR_I3 lpDataVec[iback].xp[i] = xp_seed[i];
        lpDataVec[iback].dp = dp;
        my_np_count++;
      }
    }

    int np_count;
    MPI_Allreduce(&my_np_count, &np_count, 1, MPI_INT, MPI_SUM, mpi_comm);
    assert(np_count == npnew);

  }
  break;
}

  delete[] xpnew;

}

LpInjector::LpInjector() {
  IF_RANK0 cout << "LpInjector()" << endl;

  InjectorKind = -1;
  name = "";
  Tp = 0.0;
  rhop = 0.0;
  mdot = 0.0;
  np = 0;

  DistKind = -1;
  for (int i = 0; i < 5; ++i) DistData[i] = 0.0;

  cdfKind = LpInjectorNS::NA;
  countCdf_x = NULL;
  countCdf_y = NULL;
  countCdf_n = 500;
  c0 = c1 = -1;

  GeomKind = -1;
  for (int i = 0; i < 11; ++i) GeomData[i] = 0.0;

  ResidualMass = 0.0;

  injectDt = -1.0;

  nst = 0;
  nsp = 0;
  spost = NULL;
  xp = NULL;
  cdf_area = NULL;

  FOR_I3 {
    e1[i] = 0;
    e2[i] = 0;
  }

  triVec.clear();
  cdf_area_tri = NULL;
  cdf_area_rank = NULL;

  solver = NULL;
  u = NULL;

  overWriteVelocity_b = false;

  injector_id = -1;
  material_id = -1;
  fuel_id = -1;
  dust_id = -1;

}

LpInjector::~LpInjector() {
  IF_RANK0 cout << "~LpInjector()" << endl;
  
  if (spost != NULL) delete[] spost;
  if (xp != NULL) delete[] xp;
  if (cdf_area != NULL) delete[] cdf_area;
  if (countCdf_x != NULL) delete[] countCdf_x;
  if (countCdf_y != NULL) delete[] countCdf_y;
  if (cdf_area_tri != NULL) delete[] cdf_area_tri;
  if (cdf_area_rank != NULL) delete[] cdf_area_rank;
}

void LpInjector::initFromParams(Param * param, StaticSolver * solver, double time0, vector<CtiLiquid*> fuelVec, vector<CtiSolid*> dustVec, double (*u) [3], int injector_id)  {

  this->solver = solver;
  this->u = u;
  this->time0 = time0;
  this->injector_id = injector_id;
  assert(this->injector_id >= 0);
  CtiMaterial * material = NULL;

  string token = param->getString();
  if (token != "NAME") CERR("LP_SOLID/LIQUID.INJECTOR needs NAME. Unrecogniced token "<<token<<"\n");
  string name = param->getString(1);
  this->name = name;

  if (mpi_rank == 0)
    cout << " > initializing LP_SOLID/LIQUID.INJECTOR: " << name << endl;

  bool TP_flag = false;

  int iarg = 2;
  while (iarg < param->size()) {
    string token = param->getString(iarg++);
    if (token == "TYPE") {
      // ==========================================
      // TYPE of injector:
      // ==========================================
      token = param->getString(iarg++);
      if (token == "COUNT") {
        assert( InjectorKind == -1 );
        InjectorKind = LpInjectorNS::COUNT_INJ;
        np = param->getInt(iarg++);
        if (mpi_rank == 0)
          cout << "   > injector TYPE = COUNT " << np << endl;
      }
      else if ( token == "MDOT" ) {
        assert( InjectorKind == -1 );
        InjectorKind = LpInjectorNS::MDOT_INJ;
        mdot =  param->getDouble(iarg++);
        if (mpi_rank == 0)
          cout << "   > injector TYPE = MDOT " << mdot << endl;
      }
      else {
        CERR( "\n\n\n*********************************************************\n" <<
            "Error: unrecognized TYPE in LP_SOLID/LIQUID.INJECTOR: " << token << ". Valid types are: \n" <<
            "TYPE COUNT <n>   -- injects n particles once\n" <<
            "TYPE MDOT <mdot> -- injects at rate mdot\n" <<
            "*********************************************************\n\n\n");
      }
    }
    else if (token == "DIST") {
      // ==========================================
      // DIST of particles:
      // ==========================================
      token = param->getString(iarg++);
      if (token == "UNIFORM") {
        assert( DistKind == -1 );
        DistKind = LpInjectorNS::UNIFORM_DIST;
        DistData[0] = param->getDouble(iarg++); // <diam>
        if (mpi_rank == 0)
          cout << "   > injector DIST = UNIFORM " << " Particle diam: " << DistData[0] << endl;
      }
      else if (token == "RR") {
        assert( DistKind == -1 );
        DistKind = LpInjectorNS::RR_DIST;
        string token = param->getString(iarg++);
        if (token=="COUNT_DIST") {
          cdfKind = LpInjectorNS::CDF_COUNT_DIST;
        }
        else if (token=="MASS_DIST") {
          cdfKind = LpInjectorNS::CDF_MASS_DIST;
        } else {
          CERR("RR cdf distribution type is missing. Options: COUNT_DIST, MASS_DIST");
        }
        DistData[0] = param->getDouble(iarg++); // <dmin>
        DistData[1] = param->getDouble(iarg++); // <dmax>
        DistData[2] = param->getDouble(iarg++); // <dbar>
        if (DistData[2] <= DistData[0]) CERR("Wrong diameter values in LP_SOLID/LIQUID.INJECTOR->DIST");
        if (DistData[2] >= DistData[1]) CERR("Wrong diameter values in LP_SOLID/LIQUID.INJECTOR->DIST");
        DistData[3] = param->getDouble(iarg++); // <nrr>
        if (mpi_rank == 0)
          cout << "   > injector DIST = RR, type: " << cdfKind << 
            " (count_dist: " << LpInjectorNS::CDF_COUNT_DIST << ", mass_dist: " << LpInjectorNS::CDF_MASS_DIST << ")" << 
            " dmin: " << DistData[0] << " dmax: " << DistData[1] << " davg: " <<  DistData[2] <<
            " n: " << DistData[3] << endl;
      }
      else if (token == "GAUSSIAN") {
        assert( DistKind == -1);
        DistKind = LpInjectorNS::GAUSSIAN_DIST;
        string token = param->getString(iarg++);
        if (token=="COUNT_DIST") {
          cdfKind = LpInjectorNS::CDF_COUNT_DIST;
        }
        else if (token=="MASS_DIST") {
          cdfKind = LpInjectorNS::CDF_MASS_DIST;
        } else {
          CERR("Gaussian cdf distribution type is missing. Options: COUNT_DIST, MASS_DIST");
        }
        DistData[0] = param->getDouble(iarg++); // min diameter
        DistData[1] = param->getDouble(iarg++); // max diameter
        DistData[2] = param->getDouble(iarg++); // mean diameter
        DistData[3] = param->getDouble(iarg++); // standard deviation
        if (DistData[2] <= DistData[0]) CERR("Wrong diameter values in LP_SOLID/LIQUID.INJECTOR->DIST");
        if (DistData[2] >= DistData[1]) CERR("Wrong diameter values in LP_SOLID/LIQUID.INJECTOR->DIST");
        if (mpi_rank == 0)
          cout << "   > injector DIST = GAUSSIAN type: " << cdfKind << 
            " (count_dist: " << LpInjectorNS::CDF_COUNT_DIST << ", mass_dist: " << LpInjectorNS::CDF_MASS_DIST << ")" << 
            " dmin: " << DistData[0] << " dmax: " << DistData[1] << " davg: " << DistData[2] <<
            " dstd: " << DistData[3] << endl;
      }
      else {
        CERR( "\n\n\n*********************************************************\n" <<
            "Error: unrecognized DIST in LP_SOLID/LIQUID.INJECTOR: " << token << ". Valid dists are: \n" <<
            "DIST UNIFORM <diam>                                            -- injects constant diameter particles\n" <<
            "DIST RR COUNT_DIST/MASS_DIST <dmin> <dmax> <dbar> <nrr>        -- Rosin-Rammler distribution\n" <<
            "DIST GAUSSIAN COUNT_DIST/MASS_DIST <dmin> <dmax> <davg> <std>  -- Gaussian distribution\n" <<
            "*********************************************************\n\n\n");
      }
      build_countCdf_stuff();
    }
    else if (token == "GEOM") {
      // ==========================================
      // GEOM of particles:
      // ==========================================
      token = param->getString(iarg++);
      if (token == "POINT") {
        assert( GeomKind == -1 );
        GeomKind = LpInjectorNS::POINT_GEOMETRY;
        token = param->getString(iarg++);
        if (token != "X")
          CERR("LP_SOLID/LIQUID.INJECTOR expects X");
        GeomData[0] = param->getDouble(iarg++); // <x>
        GeomData[1] = param->getDouble(iarg++); // <y>
        GeomData[2] = param->getDouble(iarg++); // <z>
        token = param->getString(iarg++);
        if (token != "N")
          CERR("LP_SOLID/LIQUID.INJECTOR expects N");
        GeomData[3] = param->getDouble(iarg++);
        GeomData[4] = param->getDouble(iarg++);
        GeomData[5] = param->getDouble(iarg++);
        NORMALIZE((GeomData+3));
        token = param->getString(iarg++);
        if (token != "VMAG")
          CERR("LP_SOLID/LIQUID.INJECTOR expects VMAG");
        GeomData[6] = param->getDouble(iarg++);
        if (mpi_rank == 0)
          cout << "   > injector GEOM = POINT  X=" << GeomData[0] << " " << GeomData[1] << " " << GeomData[2] <<
            " N=" << GeomData[3] << " " << GeomData[4] << " " << GeomData[5] << " VMAG=" << GeomData[6] << endl;
      }
      else if (token == "RING") {
        assert( GeomKind == -1 );
        GeomKind = LpInjectorNS::RING_GEOMETRY;
        token = param->getString(iarg++);
        if (token != "X")
          CERR("LP_SOLID/LIQUID.INJECTOR expects X");
        GeomData[0] = param->getDouble(iarg++); // <x>
        GeomData[1] = param->getDouble(iarg++); // <y>
        GeomData[2] = param->getDouble(iarg++); // <z>
        token = param->getString(iarg++);
        if (token != "N")
          CERR("LP_SOLID/LIQUID.INJECTOR expects N");
        GeomData[3] = param->getDouble(iarg++); // <nx>
        GeomData[4] = param->getDouble(iarg++); // <ny>
        GeomData[5] = param->getDouble(iarg++); // <nz>
        token = param->getString(iarg++);
        if (token != "R")
          CERR("LP_SOLID/LIQUID.INJECTOR expects R");
        GeomData[6] = param->getDouble(iarg++); // <r>
        token = param->getString(iarg++);
        if (token != "ANGLE")
          CERR("LP_SOLID/LIQUID.INJECTOR expects ANGLE");
        GeomData[7] = param->getDouble(iarg++); // <cone angle>
        token = param->getString(iarg++);
        if (token != "V")
          CERR("LP_SOLID/LIQUID.INJECTOR expects V");
        GeomData[8] = param->getDouble(iarg++); // <injection speed>
        token = param->getString(iarg++);
        if (token != "SF")
          CERR("LP_SOLID/LIQUID.INJECTOR expects SF");
        GeomData[9] = param->getDouble(iarg++); // <swirl fraction: 1 == 45 degrees>
        token = param->getString(iarg++);
        if (token != "SPREAD")
          CERR("LP_SOLID/LIQUID.INJECTOR expects SPREAD");
        GeomData[10] = param->getDouble(iarg++); // <cone angle spread fraction>
        // make sure that the normal has a unity length
        double mag = 0; FOR_I3 mag += GeomData[i+3]*GeomData[i+3];
        FOR_I3 GeomData[i+3] /= sqrt(mag);
        if (mpi_rank == 0)
          cout << "   > injector GEOM = RING X=" << GeomData[0] << " " << GeomData[1] << " " << GeomData[2] <<
            " V=" << GeomData[3] << " " << GeomData[4] << " " << GeomData[5] << " R=" << GeomData[6] <<
            " ANGLE=" << GeomData[7] << " V=" << GeomData[8] << " SF=" << GeomData[9] << " SPREAD=" << GeomData[10] << endl;
      }
      else if (token == "BOX") {
        assert( GeomKind == -1 );
        GeomKind = LpInjectorNS::BOX_GEOMETRY;
        token = param->getString(iarg++);
        if (token != "X")
          CERR("LP_SOLID/LIQUID.INJECTOR expects X");
        GeomData[0] = param->getDouble(iarg++); // <xmin>
        GeomData[1] = param->getDouble(iarg++); // <xmax>
        GeomData[2] = param->getDouble(iarg++); // <ymin>
        GeomData[3] = param->getDouble(iarg++); // <ymax>
        GeomData[4] = param->getDouble(iarg++); // <zmin>
        GeomData[5] = param->getDouble(iarg++); // <zmax>
        token = param->getString(iarg++);
        if (token != "V")
          CERR("LP_SOLID/LIQUID.INJECTOR expects V");
        GeomData[6] = param->getDouble(iarg++); // <vx>
        GeomData[7] = param->getDouble(iarg++); // <vy>
        GeomData[8] = param->getDouble(iarg++); // <vz>
        if (mpi_rank == 0)
          cout << "   > injector GEOM = BOX X=" << GeomData[0] << ":" << GeomData[1] << 
            " Y=" << GeomData[2] << ":" << GeomData[3] <<
            " Z=" << GeomData[4] << ":" << GeomData[5] << 
            " V=" << GeomData[6] << " " << GeomData[7] << " " << GeomData[8] << endl;
        if ( (GeomData[0] > GeomData[1]) or (GeomData[2] > GeomData[3]) or (GeomData[4] > GeomData[5]) ) 
          CERR("BOX is defined incorrectly in LSP.STATIC_INJECTOR. Check box dimensions...");
      } 
      else if (token == "CIRCLE") {
        assert( GeomKind == -1 );
        GeomKind = LpInjectorNS::CIRCLE_GEOMETRY;
        token = param->getString(iarg++);
        if (token != "X")
          CERR("LP_SOLID/LIQUID.INJECTOR expects X");
        GeomData[0] = param->getDouble(iarg++); // <x>
        GeomData[1] = param->getDouble(iarg++); // <y>
        GeomData[2] = param->getDouble(iarg++); // <z>
        token = param->getString(iarg++);
        if (token != "N")
          CERR("LP_SOLID/LIQUID.INJECTOR expects N");
        GeomData[3] = param->getDouble(iarg++);
        GeomData[4] = param->getDouble(iarg++);
        GeomData[5] = param->getDouble(iarg++);
        NORMALIZE((GeomData+3));
        token = param->getString(iarg++);
        if (token != "R")
          CERR("LP_SOLID/LIQUID.INJECTOR expects R");
        GeomData[6] = param->getDouble(iarg++);
        token = param->getString(iarg++);
        if (token != "V")
          CERR("LP_SOLID/LIQUID.INJECTOR expects V");
        GeomData[7] = param->getDouble(iarg++);
        if (mpi_rank == 0)
          cout << "   > injector GEOM = CIRCLE X=" << GeomData[0] << " " << GeomData[1] << " " << GeomData[2] << " ,N= " << GeomData[3] << " " << GeomData[4] << " " << GeomData[5] << " ,R= " << GeomData[6] << 
          " ,V= " << GeomData[7] << endl;
      }
      else if (token == "ZONE") {
        assert( GeomKind == -1 );
        GeomKind = LpInjectorNS::ZONE_GEOMETRY;

        string zone_name = param->getString(iarg++);

        token = param->getString(iarg++);
        if (token != "VMAG")
          CERR("LP_SOLID/LIQUID.INJECTOR expects VMAG");
        GeomData[0] = param->getDouble(iarg++);

        int ierr = solver->initZoneInjector(this->cdf_area, this->spost, this->xp, this->nst, this->nsp, zone_name);
        if (ierr == 1)
          CERR("Problem in ZONE injector initiation, check the zone name")

      }
      else if (token == "PLANE") {
        assert( GeomKind == -1 );
        GeomKind = LpInjectorNS::PLANE_GEOMETRY;

        token = param->getString(iarg++);
        if (token != "X")
          CERR("LP_SOLID/LIQUID.INJECTOR expects X");
        GeomData[0] = param->getDouble(iarg++);
        GeomData[1] = param->getDouble(iarg++);
        GeomData[2] = param->getDouble(iarg++);
        
        token = param->getString(iarg++);
        if (token != "N")
          CERR("LP_SOLID/LIQUID.INJECTOR expects N");
        GeomData[3] = param->getDouble(iarg++);
        GeomData[4] = param->getDouble(iarg++);
        GeomData[5] = param->getDouble(iarg++);
        NORMALIZE((GeomData+3));

        token = param->getString(iarg++);
        if (token != "VMAG")
          CERR("LP_SOLID/LIQUID.INJECTOR expects VMAG");
        GeomData[6] = param->getDouble(iarg++);

        int ierr = solver->initPlaneInjector(cdf_area_rank, cdf_area_tri, triVec, GeomData, (GeomData+3));
        if (ierr == 1) 
          CERR("Problem in PLANE injector initiation, check the plane location and normal");

      }
      else if (token == "CONE") {
        assert( GeomKind == -1 );
        GeomKind = LpInjectorNS::CONE_GEOMETRY;
        
        bool got_x = false;
        bool got_n = false;
        bool got_a = false;
        bool got_v = false;
        int nParam = 4;
        try {
          for (int iter = 0; iter < nParam; iter++) {
            string token = param->getString(iarg++);
            if (token == "X") {
              got_x = true;
              GeomData[0] = param->getDouble(iarg++);
              GeomData[1] = param->getDouble(iarg++);
              GeomData[2] = param->getDouble(iarg++);
            }
            else if (token == "N") {
              got_n = true;
              GeomData[3] = param->getDouble(iarg++);
              GeomData[4] = param->getDouble(iarg++);
              GeomData[5] = param->getDouble(iarg++);
              NORMALIZE((GeomData+3));

              // form e1 and e2
              MiscUtils::getBestE1E2FromE0(e1, e2, (GeomData+3));
            }
            else if (token == "CONE_ANGLE") {
              got_a = true;
              GeomData[6] = param->getDouble(iarg++);
            }
            else if (token == "VMAG") {
              got_v = true;
              GeomData[7] = param->getDouble(iarg++);
            }
            else {
              throw token;
            }
          }
          if (not got_x) throw 1;
          if (not got_n) throw 2;
          if (not got_a) throw 3;
          if (not got_v) throw 4;
        }
        catch (int ierr) {
          if (ierr == 1) {
            CERR("couldn't find X for CONE");
          } else if (ierr == 2) {
            CERR("couldn't find N for CONE");
          } else if (ierr == 3) {
            CERR("couldn't find CONE_ANGLE for CONE");
          } else if (ierr == 4) {
            CERR("couldn't find VMAG for CONE");
          } else {
            CERR("unrecognized parameter for LP_SOLID/LIQUID.INJECTOR GEOM CONE");
          }
        }
        catch (string token) {
          CERR("expected parameter for CONE: X, N, CONE_ANGLE, VMAG. Make sure you have them all. reached to unrecognized token: "+token);
        }
        catch(...) {
          CERR("unknown error in GEOM=CONE");
        }

        IF_RANK0 {
          cout << "   > injector GEOM = CONE X=" << GeomData[0] << " " << GeomData[1] << " " << GeomData[2] << " ,N= " << GeomData[3] << " " << GeomData[4] << " " << GeomData[5] << " ,CONE_ANGLE= " << GeomData[6] << " ,VMAG= " << GeomData[7] << " e1: " << COUT_VEC(e1) << " e2: " << COUT_VEC(e2) << endl;
        }
      }
      else {
        CERR( "\n\n\n*********************************************************\n" <<
            "Error: unrecognized GEOM in LP_SOLID/LIQUID.INJECTOR: " << token << ". Valid geoms are: \n" <<
            "GEOM POINT X <x> <y> <z> N <nx> <ny> <nz> VMAG <vmag>\n" <<
            "GEOM RING <x> <y> <z> <nx> <ny> <nz> <r> <cone angle> <injection speed> <swirl ratio> <spread>\n" <<
            "GEOM BOX X <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> V <vx> <vy> <vz> \n" <<
            "GEOM CIRCLE X <x> <y> <z> N <nx> <ny> <nz> R <r> V <vmag> \n" <<
            "GEOM ZONE <zone-name> VMAG <vmag> \n" <<
            "GEOM CONE X <x> <y> <z> N <nx> <ny> <nz> CONE_ANGLE <angle> VMAG <vmag> \n" <<
            "GEOM PLANE X <x> <y> <z> N <nx> <ny> <nz> VMAG <vmag> \n" <<
            "*********************************************************\n\n\n");
      }
    }
    else if ( token == "TP" ) {
      // ==========================================
      // TP <temperature>
      // ==========================================
      assert( !TP_flag );
      TP_flag = true;
      Tp = param->getDouble(iarg++);
      if (mpi_rank == 0) cout << "   > injector temperature: " << Tp << endl;
    }
    else if (token == "TIME_SPAN") {
      injectDt = param->getDouble(iarg++);
      if (injectDt <= 0) CERR("TIME_SPAN should be positive.");
      if (mpi_rank == 0) cout << "   > injector time span:  " << injectDt << endl;
    }
    else if ( token == "OVER_WRITE_VELOCITY" ) {
      overWriteVelocity_b = true;
      IF_RANK0 cout << "   > WARNING: over writing the injection velocity to the fluid velocity" << endl;
    }
    else if ( token == "MATERIAL") {
      assert(material == NULL);
      token = param->getString(iarg++);
      if (param->getName() == "LP_LIQUID.INJECTOR") {
        for (int ivec = 0; ivec < fuelVec.size(); ++ivec) {
          if ( token == fuelVec[ivec]->getName()) {
            material = fuelVec[ivec];
            fuel_id = ivec;
            break;
          }
        }
      }
      else if (param->getName() == "LP_SOLID.INJECTOR") {
        for (int ivec = 0; ivec < dustVec.size(); ++ivec) {
          if ( token == dustVec[ivec]->getName()) {
            material = dustVec[ivec];
            dust_id = ivec;
            break;
          }
        }
      }
      else {
        assert(0);
      }

      material_id = (fuel_id>=0) ? fuel_id:dust_id;

      if (material == NULL) CERR("Error: LP_SOLID/LIQUID.INJECTOR could not recognize MATERIAL: "<< token << ". Make sure this material is defined for the SOLID or LIQUID phase...");
      if (mpi_rank == 0) cout << "   > injector material: " << material->getName() << endl;
    }
    else {
      if (mpi_rank == 0)
      CERR( "\n\n\n*********************************************************\n" <<
          "Error: unrecognized keyword in LP_SOLID/LIQUID.INJECTOR: " << token <<
          "\n*********************************************************\n\n\n");
    }
  }

  // now check...
  if ( InjectorKind == -1 ) {
    CERR( "\n\n\n*********************************************************\n" <<
        "Error: LP_SOLID/LIQUID.INJECTOR " << name << ": TYPE not specified. Valid types are: \n" <<
        "TYPE COUNT <n>   -- injects n particles per time step\n" <<
        "TYPE MDOT <mdot> -- injects at rate mdot\n" <<
          "*********************************************************\n\n\n");
  }

  if ( DistKind == -1 ) {
    CERR( "\n\n\n*********************************************************\n" <<
        "Error: LP_SOLID/LIQUID.INJECTOR: " << name << ": DIST not specified. Valid dists are: \n" <<
        "DIST UNIFORM <diam>                                            -- injects constant diameter particles\n" <<
        "DIST RR COUNT_DIST/MASS_DIST <dmin> <dmax> <dbar> <nrr>        -- Rosin-Rammler distribution\n" <<
        "DIST GAUSSIAN COUNT_DIST/MASS_DIST <dmin> <dmax> <davg> <std>  -- Gaussian distribution\n" <<
        "*********************************************************\n\n\n");
  }

  if ( GeomKind == -1 ) {
    CERR( "\n\n\n*********************************************************\n" <<
        "Error: LP_SOLID/LIQUID.INJECTOR: " << name << " did not specified GEOM: Valid geoms are: \n" <<
        "GEOM POINT X <x> <y> <z> N <nx> <ny> <nz> VMAG <vmag>\n" <<
        "GEOM RING <x> <y> <z> <nx> <ny> <nz> <r> <cone angle> <injection speed> <swirl ratio> \n" <<
        "GEOM BOX X <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> V <vx> <vy> <vz> \n" <<
        "GEOM CIRCLE X <x> <y> <z> N <nx> <ny> <nz> R <r> V <vmag> \n" <<
        "GEOM ZONE <zone-name> VMAG <vmag> \n" <<
        "GEOM CONE X <x> <y> <z> N <nx> <ny> <nz> CONE_ANGLE <angle> VMAG <vmag> \n" <<
        "GEOM PLANE X <x> <y> <z> N <nx> <ny> <nz> VMAG <vmag> \n" <<
        "*********************************************************\n\n\n");
  }

  if ( !TP_flag ) {
    CERR( "\n\n\n*********************************************************\n" <<
        "Error: LP_SOLID/LIQUID.INJECTOR: " << name << ": TP (temperature) not specified.\n" <<
        "*********************************************************\n\n\n");
  }

  if ( material==NULL ) {
    CERR("Error: LP_SOLID/LIQUID.INJECTOR: " << name << ": MATERIAL was not determined. Please put MATERIAL <name> in the param line...");
  }

  assert(material);
  rhop = material->calcRho(Tp);
  if (mpi_rank == 0) cout << "   > injector Lp density: " << rhop << endl;

}

void LpInjector::addNewLp(vector<LpData>& lpDataVec, const double time,const double dt,const bool verbose) {

  // time: the relative time from starting the simulation
  // Don't add particles if larger than the injection time
  // for one-time injection we set injectDt to zero after injection
  if ((injectDt >= 0) && ((time-time0) > injectDt)) return;

  int from_ind = lpDataVec.size();

  int npnew = 0;                         // number of particles injected at each time step
  double * mpnew = NULL;                 // particle mass
  double * Tpnew = NULL;                 // particle temperature
  double * my_mpnew = NULL;              // particle mass for each rank
  double * my_Tpnew = NULL;              // particle temperature for each rank
  double (*xpnew)[3] = NULL;             // particle position
  double (*upnew)[3] = NULL;             // particle velocity
  int * np_count_rank = NULL;            // np for each rank for plane injectors

  // ============================================
  // on rank 0 build the injected particle list:
  // This is necessary because random numbers are involved
  // in some of the particle distributions, and we want
  // all ranks to start with the same distribution.
  // ============================================

  IF_RANK0 {

    // if the injector is MDOT modify the residual mass

    if ( InjectorKind == LpInjectorNS::MDOT_INJ )
      ResidualMass += mdot*dt;

    // -----------------------------
    // distribution...
    // -----------------------------

    if ( DistKind == LpInjectorNS::UNIFORM_DIST ) {

      double dp   = DistData[0];
      double mp   = rhop*M_PI/6.0*pow(dp,3);

      if ( InjectorKind == LpInjectorNS::MDOT_INJ ){
        double npcheck = ResidualMass/mp;
        if (npcheck > 1000000.0) {
          cerr << "Error: LpInjector " << getName() << " is trying to inject " << npcheck << " particles this time step. Check LSP.INJECT settings." << endl;
          throw(-1);
        }
        npnew = (int)npcheck;
        ResidualMass -= npnew*mp;
      }
      else if ( InjectorKind == LpInjectorNS::COUNT_INJ ) {
        // for COUNT injectors, the np stores the count. Use it once, then
        // reset to zero...
        npnew = np;
        np = 0;
      }
      else{
        assert(0);
      }

      // allocate data 
      if ( npnew > 0 ) {
        mpnew    = new double[npnew];
        Tpnew    = new double[npnew];
        for (int j = 0; j < npnew; ++j) {
          mpnew[j]   = mp;
          Tpnew[j]   = Tp;
        }
      }
    }

    else if (DistKind == LpInjectorNS::RR_DIST ) {

      /// for RR injection, the dist_data is the
      // user-specified dmin, dmax, dbar, nrr, so add particles 
      // until the stored mass is negative...

      vector<double> dpVec; // manage memory with a vector
      bool done = (InjectorKind==LpInjectorNS::MDOT_INJ) ? (ResidualMass<0) : (np<=0);
      while (!done) {

        // get the next particle diameter...
        double dp = 0;
        if (cdfKind == LpInjectorNS::CDF_COUNT_DIST) {
          dp = diam_rr_countCdf(DistData[0],DistData[1],DistData[2],DistData[3]);
        } else if (cdfKind == LpInjectorNS::CDF_MASS_DIST) {
          dp = diam_rr_massCdf(DistData[0],DistData[1],DistData[2],DistData[3]);
        } else {
          assert(0);
        }
        const double mp = rhop*M_PI/6.0*pow(dp,3);

        dpVec.push_back(dp);

        switch(InjectorKind) {
        case LpInjectorNS::MDOT_INJ:
          ResidualMass -= mp;
          done = (ResidualMass<0);
          break;
        case LpInjectorNS::COUNT_INJ:
          // for COUNT injectors, we inject a count of particles...
          np--;
          done = (np<=0);
          break;
        default:
          assert(0);
        }

      }

      npnew = dpVec.size();

      if ( npnew > 0 ) {
        mpnew    = new double[npnew];
        Tpnew    = new double[npnew];
        for (int j = 0; j < npnew; ++j) {
          double dp = dpVec[j];
          double mp = rhop*M_PI/6.0*pow(dp,3);
          mpnew[j]   = mp;
          Tpnew[j]   = Tp;
        }
      }

    }
    else if (DistKind == LpInjectorNS::GAUSSIAN_DIST ) {

      // for Gaussian injection, the dist_data is the
      // user-specified dmin, dmax, davg, dstd, so add particles 
      // until the stored mass is smaller than one particle mass...

      vector<double> dpVec; // manage memory with a vector
      bool done = (InjectorKind==LpInjectorNS::MDOT_INJ) ? (ResidualMass<0) : (np<=0);
      while (!done) {

        // get the next particle diameter...
        const double dp = diam_gaussian(DistData[0],DistData[1],DistData[2],DistData[3]);
        const double mp = rhop*M_PI/6.0*pow(dp,3);

        dpVec.push_back(dp);

        switch(InjectorKind) {
        case LpInjectorNS::MDOT_INJ:
          ResidualMass -= mp;
          done = (ResidualMass<0);
          break;
        case LpInjectorNS::COUNT_INJ:
          // for COUNT injectors, we inject a count of particles...
          np--;
          done = (np<=0);
          break;
        default:
          assert(0);
        }

      }

      npnew = dpVec.size();

      if ( npnew > 0 ) {
        mpnew    = new double[npnew];
        Tpnew    = new double[npnew];
        for (int j = 0; j < npnew; ++j) {
          double dp = dpVec[j];
          double mp = rhop*M_PI/6.0*pow(dp,3);
          mpnew[j]   = mp;
          Tpnew[j]   = Tp;
        }
      }
    }

    else {
      assert(0);
    }

    // -----------------------
    // then geometry...
    // -----------------------

    if ( npnew > 0 ) {

      xpnew = new double[npnew][3];
      upnew = new double[npnew][3];

      switch (GeomKind) {

      case LpInjectorNS::POINT_GEOMETRY:
      {
        for (int j = 0; j < npnew; ++j) {
          FOR_I3 xpnew[j][i] = GeomData[i];
          FOR_I3 upnew[j][i] = GeomData[i+3]*GeomData[6];
        }
        break;
      }

      case LpInjectorNS::RING_GEOMETRY:
        {
          // copy the data
          double XInjector[3];
          int i = 0;
          XInjector[0] = GeomData[i++];
          XInjector[1] = GeomData[i++];
          XInjector[2] = GeomData[i++];
          double Vinjector[3];
          Vinjector[0] = GeomData[i++];
          Vinjector[1] = GeomData[i++];
          Vinjector[2] = GeomData[i++];
          double RInjector        = GeomData[i++];
          double ConeAngleDegrees = GeomData[i++];
          double InjectSpeed      = GeomData[i++];
          double SwirlRatio       = GeomData[i++];
          double spreadFrac       = GeomData[i++];
          assert( i==11 );

          // generate random numbers for both location and injection angle
          double * RandX = new double[npnew];
          for ( int j =0; j < npnew; ++j) {
            RandX[j] = 2.0*M_PI*(double)(rand())/(double)(RAND_MAX); // 0 to 2pi
          }

          // Box-Muller
          double * RandA = new double[npnew];
          for ( int j =0; j < npnew; ++j) {
            // not sure id rand() can return 0, so make sure it doesn't as follows...
            double u1 = (double(rand())+1.0)/(double(RAND_MAX)+1.0);
            double u2 = (double(rand())+1.0)/(double(RAND_MAX)+1.0);
            double normal_rand_number = sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
            // a second normal (not used here) would be 
            // = sqrt(-2.0*log(u1))*sin(2.0*M_PI*u2);
            // scatter the points about the user-defined cone half-angle according to spreadFrac
            RandA[j] = ConeAngleDegrees*M_PI/180.0*(1.0 + spreadFrac*normal_rand_number);
            // should be very rare, but limit between -90 and 90...
            RandA[j] = max(-M_PI/2.0,min(M_PI/2.0,RandA[j]));
          }

          //dumpHist0(RandA,npnew,"RandA");
          //dumpHist0(RandX,npnew,"RandX");
          //getchar();

          // build a set of orthogonal vectors
          // first vector is just the injection direction (already normalized)...
          double e1[3]; FOR_I3 e1[i] = Vinjector[i];

          // for e2, first choose the min component of e1...
          double e3[3] = { 0.0, 0.0, 0.0 };
                if      (fabs(e1[0]) <= min(fabs(e1[1]),fabs(e1[2])))
                  e3[0] = 1.0;
                else if (fabs(e1[1]) <= min(fabs(e1[0]),fabs(e1[2])))
                  e3[1] = 1.0;
                else
                  e3[2] = 1.0;

          // e2 is the cross product 
          double e2[3];
          e2[0] = e1[1]*e3[2] - e1[2]*e3[1];
          e2[1] = e1[2]*e3[0] - e1[0]*e3[2];
          e2[2] = e1[0]*e3[1] - e1[1]*e3[0];
          double inv_mag = 1.0/sqrt(DOT_PRODUCT(e2,e2));
          FOR_I3 e2[i] *= inv_mag;

          // e3 is e1 x e2... 
          e3[0] = e1[1]*e2[2] - e1[2]*e2[1];
          e3[1] = e1[2]*e2[0] - e1[0]*e2[2];
          e3[2] = e1[0]*e2[1] - e1[1]*e2[0];

          // should already be normal...
          assert( fabs(1.0-sqrt(DOT_PRODUCT(e3,e3))) < 1.0E-10);

          // break the velocity into a tangential part and a "normal" part. The normal
          // part is the projection in the "e1-R" plane...
          const double u_n = InjectSpeed/pow((1.0+pow(SwirlRatio,2)),0.5);
          const double u_t = u_n*SwirlRatio;

          // for all the new particles
          for (int j = 0; j < npnew; ++j) {

            // first randomly generate particles on a 360 ring
            double Rdirection[3]; FOR_I3 Rdirection[i] = cos(RandX[j])*e2[i] + sin(RandX[j])*e3[i];
            FOR_I3 xpnew[j][i] = XInjector[i] + RInjector*Rdirection[i];

            // swirl direction T = e1 x R...
            double Tdirection[3];
            Tdirection[0] = e1[1]*Rdirection[2] - e1[2]*Rdirection[1];
            Tdirection[1] = e1[2]*Rdirection[0] - e1[0]*Rdirection[2];
            Tdirection[2] = e1[0]*Rdirection[1] - e1[1]*Rdirection[0];

            FOR_I3 upnew[j][i] = u_t*Tdirection[i] + u_n*(cos(RandA[j])*e1[i] + sin(RandA[j])*Rdirection[i]);

          }

          delete[] RandX;
          delete[] RandA;

        }
        break;
    
        case LpInjectorNS::BOX_GEOMETRY:
        {
          for (int j = 0; j < npnew; ++j) {
            double rand_x[3];
            FOR_I3 rand_x[i] = rand()/double(RAND_MAX);
            FOR_I3 xpnew[j][i] = GeomData[2*i] + (GeomData[2*i+1] - GeomData[2*i])*rand_x[i];
            FOR_I3 upnew[j][i] = GeomData[i+6];
          }
        }
        break;

        case LpInjectorNS::CIRCLE_GEOMETRY:
        {
          double X[3], n[3], R;
          FOR_I3 X[i] = GeomData[i]; 
          FOR_I3 n[i] = GeomData[i+3];
          R = GeomData[6];
          int i_min = -1;
          double n_min = 1e20;
          FOR_I3 {
            if (abs(n[i]) < n_min) {
              i_min = i;
              n_min = n[i];
            }
          }
          double dir[3];
          dir[i_min] = 0.0;
          if (i_min == 0) {
            dir[1] = n[2];
            dir[2] = -n[1];
          } else if (i_min == 1) {
            dir[0] = n[2];
            dir[2] = -n[0];
          } else if (i_min == 2) {
            dir[0] = n[1];
            dir[1] = -n[0];
          } else {assert(0);}
          double mag_dir = MAG(dir);
          assert(mag_dir>0);
          FOR_I3 dir[i] /= mag_dir;
          assert(DOT_PRODUCT(dir,n) < 1e-12);
          for (int j = 0; j < npnew; ++j) {
            double rand_r_ = rand()/double(RAND_MAX);
            double rand_t = rand()/double(RAND_MAX);   // random variables for r and theta
            double r = R * sqrt(rand_r_);
            double theta = rand_t * 2. * M_PI;
            double Rot[3][3];
            // build three dimensional rotation matrix
            FOR_I3 Rot[i][i] = cos(theta) + n[i]*n[i]*(1.-cos(theta));
            Rot[0][1] = n[0]*n[1]*(1.-cos(theta)) - n[2]*sin(theta);
            Rot[1][0] = n[0]*n[1]*(1.-cos(theta)) + n[2]*sin(theta);
            Rot[0][2] = n[0]*n[2]*(1.-cos(theta)) + n[1]*sin(theta);
            Rot[2][0] = n[0]*n[2]*(1.-cos(theta)) - n[1]*sin(theta);
            Rot[2][1] = n[1]*n[2]*(1.-cos(theta)) + n[0]*sin(theta);
            Rot[1][2] = n[1]*n[2]*(1.-cos(theta)) - n[0]*sin(theta);
            double rand_dir[3] = MATRIX_PRODUCT(Rot,dir);
            assert(MAG(rand_dir) - 1. < 1e-12);
            assert(DOT_PRODUCT(rand_dir,n) < 1e-12);
            FOR_I3 xpnew[j][i] = X[i] + r * rand_dir[i];
            FOR_I3 upnew[j][i] = n[i] * GeomData[7];
          }
        }
        break;

        case LpInjectorNS::ZONE_GEOMETRY:
        {
          assert(this->nst>0);
          assert(this->nsp>0);
          assert(spost);
          assert(xp);
          assert(cdf_area);
          for (int j = 0; j < npnew; ++j) {
            // first randomly choose a tri then choose the point in the tri randomly
            double rand_val[3];
            FOR_I3 rand_val[i] = rand()/double(RAND_MAX);
            int ist_rand = 0;
            for (int ist = 1; ist < this->nst+1; ist++) {
              if (rand_val[0] < cdf_area[ist]) {
                ist_rand = ist-1;
                break;
              }
            }
            //IF_RANK0 cout << "rand_val: " << COUT_VEC(rand_val) << " ist_rand: " << ist_rand << endl;
            assert((ist_rand>=0)&&(ist_rand<this->nst));

            const double * const x0 = xp[spost[ist_rand][0]];
            const double * const x1 = xp[spost[ist_rand][1]];
            const double * const x2 = xp[spost[ist_rand][2]];

            // set the normal to be inward
            double normal[3] = TRI_NORMAL_2(x0,x2,x1);
            double mag = MAG(normal);
            FOR_I3 normal[i] /= mag;
            
            // locate a random point inside the tri by partitioning 1 into three segments
            double rand_min = min(rand_val[1],rand_val[2]);
            double rand_max = max(rand_val[1],rand_val[2]);
            // move the particle slightly inside
            double small = 1e-7;
            FOR_I3 xpnew[j][i] = rand_min*x0[i] + (rand_max-rand_min)*x1[i] + (1. - rand_max)*x2[i] + small*normal[i];
            FOR_I3 upnew[j][i] = GeomData[0]*normal[i];
          }
        }
        break;

        case LpInjectorNS::PLANE_GEOMETRY:
        {
          assert(cdf_area_rank);
          assert(np_count_rank==NULL);
          np_count_rank = new int [mpi_size];
          FOR_RANK np_count_rank[rank] = 0;
  
          for (int ip = 0; ip < npnew; ++ip) {
            double rand_rank = double(rand())/(double(RAND_MAX) + 1.0);
            FOR_RANK {
              if (rand_rank < cdf_area_rank[rank]) {
                assert(rank < mpi_size);
                np_count_rank[rank]++;
                break;
              }
            }
          }

          int np_sum = 0;
          FOR_RANK np_sum += np_count_rank[rank];
          assert(np_sum == npnew);

        }
        break;

        case LpInjectorNS::CONE_GEOMETRY:
        {
          for (int j = 0; j < npnew; ++j) {
            double rand_val[2];
            FOR_I2 rand_val[i] = rand()/double(RAND_MAX);

            double theta = (rand_val[0] - 0.5) * 2.0 * (GeomData[6]/180.*M_PI);
            double phi = rand_val[1] * 2.0 * M_PI;
            double sin_theta = sin(theta);
            double cos_theta = cos(theta);
            double sin_phi = sin(phi);
            double cos_phi = cos(phi);
            // if we take e1, e2, and n as the main axes, 
            // we rotate in the e1-n plane with theta, then rotate in e1-e2 plane with phi
            // so n = (0,0,1) goes to...
            double normal[3] = {0,0,0};
            FOR_I3 normal[i] += -cos_phi*sin_theta*e1[i];
            FOR_I3 normal[i] += -sin_phi*sin_theta*e2[i];
            FOR_I3 normal[i] += cos_theta*GeomData[3+i];
            assert((DOT_PRODUCT(normal,normal)-1.)<1e-6);

            FOR_I3 xpnew[j][i] = GeomData[i];
            FOR_I3 upnew[j][i] = GeomData[7]*normal[i];
          }
          
        }
        break;

      default:
        assert(0);
      }

    }

  }

  // ======================================
  // rank 0 bcast to everyone...
  // every rank should deside which particles
  // should be injected in their rank.
  // for plane injector each rank should 
  // determine the location of injection as well.
  // ======================================

  MPI_Bcast(&npnew,1,MPI_INT,0,mpi_comm);
  if (npnew == 0) return;

  switch (GeomKind) {

    case LpInjectorNS::PLANE_GEOMETRY:  // in case of a plane injector the location of the points are not still determined
                          // so each rank prepares particles for itself
    {

    //const double * const xp_plane = GeomData;
    const double * const np_plane = GeomData+3;

    assert(cdf_area_tri);
    IF_RANK0 {
      assert(np_count_rank);
      assert(cdf_area_rank);
    }

    int my_npnew;
    MPI_Scatter(np_count_rank, 1, MPI_INT, &my_npnew, 1, MPI_INT, 0, mpi_comm);
    
    //if (mpi_rank != 0) {
    my_mpnew    = new double[my_npnew];
    my_Tpnew    = new double[my_npnew];
    //}

    int * send_disp;
    IF_RANK0 {
      send_disp = new int [mpi_size];
      send_disp[0] = 0;
      for (int irank = 1; irank < mpi_size; irank++) {
        send_disp[irank] = send_disp[irank-1] + np_count_rank[irank-1];
      }
      assert((send_disp[mpi_size-1]+np_count_rank[mpi_size-1])==npnew);
    }
    MPI_Scatterv(mpnew, np_count_rank, send_disp, MPI_DOUBLE, my_mpnew, my_npnew, MPI_DOUBLE, 0, mpi_comm);
    MPI_Scatterv(Tpnew, np_count_rank, send_disp, MPI_DOUBLE, my_Tpnew, my_npnew, MPI_DOUBLE, 0, mpi_comm);
    
    IF_RANK0 delete[] send_disp;

    int npold = lpDataVec.size();
    lpDataVec.resize(npold+my_npnew);

    int ipnew = 0;
    while (ipnew < my_npnew) {
      
      int itri = -1;
      double rand_tri = double(rand()) / (double(RAND_MAX) + 1.0);
      for (int ii = 0; ii < triVec.size(); ++ii) {
        if (rand_tri < cdf_area_tri[ii]) {
          itri = ii;
          break;
        }
      }
      assert((itri>=0)&&(itri<triVec.size()));

      // locate a random point inside the tri by partitioning 1 into three segments
      double rand_val[2]; 
      FOR_I2 rand_val[i] = double(rand()) / (double(RAND_MAX) + 1.0);
      double rand_min = min(rand_val[0],rand_val[1]);
      double rand_max = max(rand_val[0],rand_val[1]);
  
      const double * const x0 = triVec[itri].first.x0;
      const double * const x1 = triVec[itri].first.x1;
      const double * const x2 = triVec[itri].first.x2;

      double xp_seed[3];
      FOR_I3 xp_seed[i] = rand_min*x0[i] + (rand_max-rand_min)*x1[i] + (1. - rand_max)*x2[i];

      int icv_closest;
      const int isInside = solver->pointIsInside(xp_seed, icv_closest);

      if ( isInside ) {
        assert( (icv_closest >= 0) && (icv_closest < solver->ncv));
        lpDataVec[npold+ipnew].icv = icv_closest;
        char * material_keep_flag = (char *)&lpDataVec[npold+ipnew].flag;
        material_keep_flag[0] = material_id;
        material_keep_flag[1] = KEEP;
        short int * injector_flag = (short int *)&lpDataVec[npold+ipnew].flag;
        injector_flag[1] = injector_id;
        lpDataVec[npold+ipnew].mp = my_mpnew[ipnew];
        lpDataVec[npold+ipnew].dp = pow( my_mpnew[ipnew] / (M_PI/6.0*rhop) , 1./3.);
        lpDataVec[npold+ipnew].Tp = my_Tpnew[ipnew];
        FOR_I3 lpDataVec[npold+ipnew].xp[i] = xp_seed[i];
        FOR_I3 lpDataVec[npold+ipnew].up[i] = GeomData[6]*np_plane[i];
        ipnew++;
      }
    }

    int np_count;
    MPI_Allreduce(&ipnew, &np_count, 1, MPI_INT, MPI_SUM, mpi_comm);
    assert(np_count == npnew);

    delete[] np_count_rank;
    }
    break;

    default:      // in the default case rank 0 sets all mp, Tp, xp, up and we just decide which rank owns which particles
    {
    if (mpi_rank != 0) {
      mpnew    = new double[npnew];
      Tpnew    = new double[npnew];
      xpnew    = new double[npnew][3];
      upnew    = new double[npnew][3];
    }

    MPI_Bcast(mpnew,            npnew,MPI_DOUBLE,0,mpi_comm);
    MPI_Bcast(Tpnew,            npnew,MPI_DOUBLE,0,mpi_comm);
    MPI_Bcast((double*)xpnew, 3*npnew,MPI_DOUBLE,0,mpi_comm);
    MPI_Bcast((double*)upnew, 3*npnew,MPI_DOUBLE,0,mpi_comm);

    // add to the lpDataVec...

    int my_np_count = 0;
    for (int ip = 0; ip < npnew; ++ip) {
      const double xp_seed[3] = {xpnew[ip][0], xpnew[ip][1], xpnew[ip][2]};
      int icv_closest;
      const int isInside = solver->pointIsInside(xp_seed,icv_closest);
      if ( isInside ) {
        assert( (icv_closest >= 0) && (icv_closest < solver->ncv));
        lpDataVec.resize(lpDataVec.size()+1);
        int iback = lpDataVec.size()-1;
        lpDataVec[iback].icv = icv_closest;
        char * material_keep_flag = (char *)&lpDataVec[iback].flag;
        material_keep_flag[0] = material_id;
        material_keep_flag[1] = KEEP;
        short int * injector_flag = (short int *)&lpDataVec[iback].flag;
        injector_flag[1] = injector_id;
        //short int * material_injector_flag = (short int *)&lpDataVec[iback].flag;
        //material_injector_flag[0] = material_id;
        //material_injector_flag[1] = injector_id;
        lpDataVec[iback].mp = mpnew[ip];
        lpDataVec[iback].dp = pow( mpnew[ip] / (M_PI/6.0*rhop) , 1./3.);
        lpDataVec[iback].Tp = Tpnew[ip];
        FOR_I3 lpDataVec[iback].xp[i] = xp_seed[i];
        FOR_I3 lpDataVec[iback].up[i] = upnew[ip][i];
        my_np_count++;
      }
    }

    int np_count;
    MPI_Allreduce(&my_np_count, &np_count, 1, MPI_INT, MPI_SUM, mpi_comm);
    assert(np_count == npnew);

    }
    break;
  }

  if (overWriteVelocity_b) overWriteVelocity(lpDataVec,from_ind);

  delete[] mpnew;
  delete[] Tpnew;
  delete[] xpnew;
  delete[] upnew;
  delete[] my_mpnew;
  delete[] my_Tpnew;

}

double LpInjector::diam_rr_countCdf(const double &dmin,const double &dmax,const double &davg,const double &q) {

  // clipped Rosin-Rammler (Weibull) distribution
  //
  // dmin --> minimum diameter
  // dmax --> maximum diameter
  // davg --> scale parameter (not a true mean diameter)
  // q    --> shape parameter (exponent in distribution)

  assert((c0>=0)&&(c0<1));
  assert((c1>0)&&(c1<=1));
  const double r = double(rand())/double(RAND_MAX); // 0 < r < 1
  const double rand_val = c0 + r*(c1-c0);

  double d_sel = davg*pow(-log(1.-rand_val),1.0/q);

  d_sel = min(d_sel,dmax);
  d_sel = max(d_sel,dmin);

  return d_sel;
}

double LpInjector::diam_rr_massCdf(const double &dmin,const double &dmax,const double &davg,const double &q) {

  // clipped Rosin-Rammler (Weibull) distribution
  //
  // dmin --> minimum diameter
  // dmax --> maximum diameter
  // davg --> scale parameter (not a true mean diameter)
  // q    --> shape parameter (exponent in distribution)

  assert((c0>=0)&&(c0<1));
  assert((c1>0)&&(c1<=1));
  const double r = double(rand())/double(RAND_MAX); // 0 < r < 1
  const double rand_val = c0 + r*(c1-c0);

  int bin_sel = -1;

  for (int ibin = 1; ibin < countCdf_n; ibin++) {

    if (countCdf_y[ibin] >= rand_val) {
      bin_sel = ibin - 1;
      break;
    }
      
  }
  assert(bin_sel >= 0);

  // add randomness inside each bin
  const double dbin = countCdf_x[2] - countCdf_x[1];
  const double r_inbin = double(rand())/double(RAND_MAX); // 0 < r_inbin < 1
  double d_sel = countCdf_x[bin_sel] + r_inbin*dbin;
  d_sel = min(d_sel,dmax);
  d_sel = max(d_sel,dmin);

  return d_sel;
  
}

double LpInjector::diam_gaussian(const double &dmin, const double &dmax, const double &davg, const double &d_std) {
  
  assert((c0>=0)&&(c0<1));
  assert((c1>0)&&(c1<=1));
  const double r = double(rand())/double(RAND_MAX); // 0 < r < 1
  const double rand_val = c0 + (c1 - c0)*r;

  int bin_sel = -1;

  for (int ibin = 1;ibin < countCdf_n; ibin++) { 

    if (countCdf_y[ibin] >= rand_val) {
      bin_sel = ibin - 1;
      break;
    }

  }
  assert(bin_sel >= 0);

  // add randomness inside each bin
  const double dbin = countCdf_x[2] - countCdf_x[1];
  const double r_inbin = double(rand())/double(RAND_MAX); // 0 < r_inbin < 1
  double d_sel = countCdf_x[bin_sel] + r_inbin*dbin;
  d_sel = min(d_sel,dmax);
  d_sel = max(d_sel,dmin);

  return d_sel;
}

void LpInjector::build_countCdf_stuff() {
  
  IF_RANK0 cout << "   > build_countCdf_stuff()" << endl;

  assert(countCdf_x==NULL);
  assert(countCdf_y==NULL);
  assert(countCdf_n>0);

  switch (DistKind) {
  case LpInjectorNS::UNIFORM_DIST:
    break;
  case LpInjectorNS::RR_DIST:
    {
      const double dmin = DistData[0];
      const double dmax = DistData[1];
      const double davg = DistData[2];
      const double n = DistData[3];

      switch (cdfKind) {
      case LpInjectorNS::CDF_COUNT_DIST:
        // count cdf is analytical and only need the beginning and the end of y axis
        c0 = 1. - exp(-pow((dmin/davg),n));
        c1 = 1. - exp(-pow((dmax/davg),n));
        break;
      case LpInjectorNS::CDF_MASS_DIST:
        // count cdf should be build from analytical mass cdf
        {
        countCdf_x = new double [countCdf_n];
        countCdf_y = new double [countCdf_n];
        for (int ibin = 0; ibin < countCdf_n; ibin++) {
          countCdf_x[ibin] = 0.;
          countCdf_y[ibin] = 0.;
        }

        const double dmax_ext = dmax*5.;
        const double dmin_ext = 0.;
        const double dbin = (dmax_ext-dmin_ext)/double(countCdf_n);

        const int dmin_ind = int(floor(dmin / dbin)); 
        const int dmax_ind = int(floor(dmax / dbin)); 
        
        for (int ibin = 0; ibin < countCdf_n; ibin++) {
          countCdf_x[ibin] = dmin_ext + (double(ibin))*dbin;
          assert(countCdf_y[ibin]==0.);
          if (ibin > 0) {
            countCdf_y[ibin] = countCdf_y[ibin-1] + n*pow(countCdf_x[ibin]/davg,n)/pow(countCdf_x[ibin],4)*exp(-pow(countCdf_x[ibin]/davg,n));
          }
        }
        
        // normalize the cdf with the last element of cdf
        for (int ibin = 0; ibin < countCdf_n; ibin++)
          countCdf_y[ibin] /= countCdf_y[countCdf_n-1];

        c0 = countCdf_y[dmin_ind];
        c1 = countCdf_y[dmax_ind];
        }
        break;
      default:
        assert(0);
      }
    }
    break;
  case LpInjectorNS::GAUSSIAN_DIST:
    {
      switch (cdfKind) {
      case LpInjectorNS::CDF_COUNT_DIST:
        {
        const double dmin = DistData[0];
        const double dmax = DistData[1];
        const double davg = DistData[2];
        const double std  = DistData[3];

        // only count cdf is supported for now
        countCdf_x = new double [countCdf_n];
        countCdf_y = new double [countCdf_n];
        for (int ibin = 0; ibin < countCdf_n; ibin++) {
          countCdf_x[ibin] = 0.;
          countCdf_y[ibin] = 0.;
        }

        const double dmax_ext = dmax*5.;
        const double dmin_ext = 0.;
        const double dbin = (dmax_ext-dmin_ext)/double(countCdf_n);

        const int dmin_ind = int(floor(dmin / dbin)); 
        const int dmax_ind = int(floor(dmax / dbin)); 
        
        for (int ibin = 0; ibin < countCdf_n; ibin++) {
          countCdf_x[ibin] = dmin_ext + (double(ibin))*dbin;
          countCdf_y[ibin] = 0.5*(1.0 + error_function2((countCdf_x[ibin]-davg)/(sqrt(2)*std)));
        }
        
        // normalize the cdf with the last element of cdf
        for (int ibin = 0; ibin < countCdf_n; ibin++)
          countCdf_y[ibin] /= countCdf_y[countCdf_n-1];

        c0 = countCdf_y[dmin_ind];
        c1 = countCdf_y[dmax_ind];
        }
        break;
      case LpInjectorNS::CDF_MASS_DIST:
        {
        const double dmin = DistData[0];
        const double dmax = DistData[1];
        const double davg = DistData[2];
        const double std  = DistData[3];

        // only count cdf is supported for now
        countCdf_x = new double [countCdf_n];
        countCdf_y = new double [countCdf_n];
        for (int ibin = 0; ibin < countCdf_n; ibin++) {
          countCdf_x[ibin] = 0.;
          countCdf_y[ibin] = 0.;
        }

        const double dmax_ext = dmax*5.;
        const double dmin_ext = 0.;
        const double dbin = (dmax_ext-dmin_ext)/double(countCdf_n);

        const int dmin_ind = int(floor(dmin / dbin)); 
        const int dmax_ind = int(floor(dmax / dbin)); 
        
        for (int ibin = 0; ibin < countCdf_n; ibin++) {
          countCdf_x[ibin] = dmin_ext + (double(ibin))*dbin;
          assert(countCdf_y[ibin]==0.);
          if (ibin > 0) {
            countCdf_y[ibin] = countCdf_y[ibin-1] + exp(-pow(countCdf_x[ibin]-davg,2)/(2.*pow(std,2))) / pow(countCdf_x[ibin],3);
          }
        }
        
        // normalize the cdf with the last element of cdf
        for (int ibin = 0; ibin < countCdf_n; ibin++)
          countCdf_y[ibin] /= countCdf_y[countCdf_n-1];

        c0 = countCdf_y[dmin_ind];
        c1 = countCdf_y[dmax_ind];
        }
        break;
      default:
        assert(0);
      }
    }
    break;
  default:
    assert(0);
  }

}

void LpInjector::overWriteVelocity(vector<LpData>& lpDataVec,int from_ind, bool verbose) {
    
  assert(u != NULL);
  IF_RANK0 
    if (verbose)
      cout << "   > over writing injector velocity..." << endl;

  for (int ivec = from_ind; ivec < lpDataVec.size(); ivec++) {
    FOR_I3 lpDataVec[ivec].up[i] = u[lpDataVec[ivec].icv][i];
  }

}

double error_function2(const double eta) {

  double z = fabs(eta);
  double t = 1.0/(1.0+0.5*z);
  double erfcc = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+
              t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+
                        t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
  if (eta < 0.0) erfcc = 2.0-erfcc;

  double err_func = 1.0 - erfcc;
  return err_func;

}


