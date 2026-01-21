#include "CTI.hpp"
using namespace CTI;

#include "Octree.hpp"
#include "StaticSolver.hpp"

// reasonable abstract base-class for solver bcs...

class SampleRhs {
public:
  double a;
  double b;
};

class SampleState {
public:
  int it;
  double a;
  double b;  
  double c;
  double d;
  double x[3];

  // since the sample state is of a mixed 
  // type, we;re going to need to do the 
  // pack and unpack through char.. 
  static int data_size() { return 4*8 + 1*4 + 8*3; }

  // pack and unpack routines moved outside of the 
  // class definition to support a more general interface..
};

inline void pack_class(char* buf, const SampleState* s, 
                       const int icv) { 
  assert( sizeof(int) == 4);
  assert( sizeof(double) == 8);
  memcpy(&buf[0] ,&s[icv].it,4);
  memcpy(&buf[4] ,&s[icv].a, 8);
  memcpy(&buf[12],&s[icv].b, 8);
  memcpy(&buf[20],&s[icv].c, 8);
  memcpy(&buf[28],&s[icv].d, 8);

  double * xx = (double*)(&buf[36]);
  for (int i =0 ; i < 3; ++i)
    xx[i] = s[icv].x[i];
}

inline void pack_class(char* buf, const SampleState* s, 
                       const int icv, const double* R) { 
  assert( sizeof(int) == 4);
  assert( sizeof(double) == 8);
  memcpy(&buf[0] ,&s[icv].it,4);
  memcpy(&buf[4] ,&s[icv].a, 8);
  memcpy(&buf[12],&s[icv].b, 8);
  memcpy(&buf[20],&s[icv].c, 8);
  memcpy(&buf[28],&s[icv].d, 8);

  double * xx = (double*)(&buf[36]);
  MiscUtils::matVecMult(xx,R,s[icv].x);
}

inline void unpack_class(const char* buf, SampleState* s, const int icv) { 
  assert( sizeof(int) == 4);
  assert( sizeof(double) == 8);
  memcpy(&s[icv].it  , &buf[0] ,4);
  memcpy(&s[icv].a   , &buf[4] ,8);
  memcpy(&s[icv].b   , &buf[12],8);
  memcpy(&s[icv].c   , &buf[20],8);
  memcpy(&s[icv].d   , &buf[28],8);

  double * xx = (double*)(&buf[36]);
  for (int i =0 ; i < 3; ++i)
    s[icv].x[i] = xx[i]; 

}

class SampleSolverBaseBc {
public:
  BfZone * bfZonePtr; // a pointer to the corresponding bfZone data for this bc
  SampleSolverBaseBc(BfZone& bfZone) {
    bfZonePtr = &bfZone; // just store a ptr
  }
  virtual ~SampleSolverBaseBc() {
  }
  // pull the name from bfZone...
  string getName() const {
    assert(bfZonePtr);
    return bfZonePtr->getName();
  }
  // some virtual methods that each instance must provide...
  virtual void calcRhs() = 0;
};

// lagrangian particles can be managed in a class 
// that starts with members int icv, double xp[3], double xp0[3]...

class MyLsp {
public:
  int icv;
  int flag; // put the flag here if you need it or not: it is required
  double xp[3];
  
  double xp0[3];
  int c;
  double a;

  double b;
  double e[3];
  double f[3];
};

class MyState {
public:
  int i;
  double d;
  double d3[3];
};

class SampleSolver : public StaticSolver {
public:

  vector<SampleSolverBaseBc*> bcVec;
  SampleState * state;

  int np,np_max;
  MyLsp * lsp;

  MyState * my_state;
  double * my_d;
  double (*my_d3)[3];
  int * my_i;

  double (*my_n_fa)[3];

  SampleSolver() {

    COUT1("SampleSolver()");
    
    state = NULL;
  
    // particles...
    np = np_max = 0; // np_max for memory management
    lsp = NULL;
    //StaticSolver::registerLp<MyLsp>(lsp,"lsp",np);  // start : must do this (assumes icv, flag, xp[3] included)
    //StaticSolver::registerLpData<MyLsp>(lsp,lsp->xp0,"xp0"); // double 3 : optional
    //StaticSolver::registerLpData<MyLsp>(lsp,lsp->c,"c"); // int 
    //StaticSolver::registerLpData<MyLsp>(lsp,lsp->a,"a"); // double 

    // register using ncv 
    assert(ncv == 0);
    my_state = NULL;
    my_d = NULL;
    my_d3 = NULL;
    my_i = NULL;
    StaticSolver::registerCvData(my_state,my_state->i,"my_state::i",READ_DATA|WRITE_DATA);
    StaticSolver::registerCvData(my_state,my_state->d,"my_state::d",READ_DATA|WRITE_DATA);
    StaticSolver::registerCvData(my_state,my_state->d3,"my_state::d3",READ_DATA|WRITE_DATA);
    StaticSolver::registerCvData(my_i,"my_i",READ_DATA|WRITE_DATA);
    StaticSolver::registerCvData(my_d,"my_d",READ_DATA|WRITE_DATA);
    StaticSolver::registerCvData(my_d3,"my_d3",READ_DATA|WRITE_DATA);
    StaticSolver::registerSignedFaData(my_n_fa,"my_n_fa",READ_DATA|WRITE_DATA);

    my_n_fa = NULL;
  }

  SampleSolver(const int icg) : StaticSolver(icg) {

    state = NULL;
    np = np_max = 0; // np_max for memory management
    lsp = NULL;
    assert(ncv == 0);
    my_state = NULL;
    my_d = NULL;
    my_d3 = NULL;
    my_i = NULL;
    my_n_fa = NULL;

  }

  void initCtiRegisterManagedData() {

    // also try to register data where CtiRegister manages memory...
    
    registerI("my_int",READWRITE_DATA);
    CtiRegister::CtiData * data0 = CtiRegister::getRegisteredCtiData("my_int"); // we know it is registered
    assert(data0->getType() == I_DATA);
    FOR_ICV data0->i() = 4321;
    //cout << "FOUND my_int: " << *data0 << endl;

    /*
    registerCvIN("index_cv",READWRITE_DATA);
    CtiRegister::CtiData * data1 = CtiRegister::getRegisteredCtiData("index_cv");
    assert(data1->getType() == IN_DATA);
    assert(data1->getUnindexedTopology() == CV_DATA);
    FOR_ICV data1->in(icv) = icv;
    //cout << "FOUND index_cv: " << *data1 << endl;
    //
    */

    registerD("my_double",READWRITE_DATA);
    CtiRegister::CtiData * data2 = CtiRegister::getRegisteredCtiData("my_double");
    assert(data2->getType() == D_DATA);
    FOR_ICV data2->d() = double(1234);
    //cout << "FOUND my_double: " << *data2 << endl;

    /*
    registerD3("my_double3",READWRITE_DATA);
    CtiRegister::CtiData * data3 = CtiRegister::getRegisteredCtiData("my_double3");
    assert(data3->getType() == D3_DATA);
    FOR_I3 data3->d3(i) = double(i);
    //cout << "FOUND my_double3: " << *data3 << endl;
    */

    //cout << "FOUND my_double: " << *data2 << endl;
    registerCvDN("rho_cv",READWRITE_DATA);
    CtiRegister::CtiData * data4 = CtiRegister::getRegisteredCtiData("rho_cv");
    assert(data4->getType() == DN_DATA);
    assert(data4->getUnindexedTopology() == CV_DATA);
    FOR_ICV data4->dn(icv) = double(icv);
    //cout << "FOUND rho_cv: " << *data4 << endl;
    
    registerCvDN3("rhou_cv",READWRITE_DATA);
    CtiRegister::CtiData * data5 = CtiRegister::getRegisteredCtiData("rhou_cv");
    assert(data5->getType() == DN3_DATA);
    assert(data5->getUnindexedTopology() == CV_DATA);
    FOR_ICV FOR_I3 data5->dn3(icv,i) = double(icv)*double(3)+double(i);
    //cout << "FOUND rhou_cv: " << *data5 << endl;

    FOR_IZONE(bfZoneVec) {
         
      bfZoneVec[izone].registerBfDN(bfZoneVec[izone].getName()+":rho_bf",READWRITE_DATA);
      CtiRegister::CtiData * data6 = CtiRegister::getRegisteredCtiData(bfZoneVec[izone].getName()+":rho_bf");
      assert(data6->getType() == DN_DATA);
      assert(data6->getUnindexedTopology() == BF_DATA);
      for (int ibf = 0; ibf < bfZoneVec[izone].nbf; ++ibf)  data6->dn(ibf) = double(ibf);
      //cout << "FOUND rho_bf[" << izone << "]: " << *data6 << endl;
      
      bfZoneVec[izone].registerBfDN3(bfZoneVec[izone].getName()+":rhou_bf",READWRITE_DATA);
      CtiRegister::CtiData * data7 = CtiRegister::getRegisteredCtiData(bfZoneVec[izone].getName()+":rhou_bf");
      assert(data7->getType() == DN3_DATA);
      assert(data7->getUnindexedTopology() == BF_DATA);
      for (int ibf = 0; ibf < bfZoneVec[izone].nbf; ++ibf) FOR_I3 data7->dn3(ibf,i) = double(ibf)*double(3)+double(i);
      //cout << "FOUND rhou_bf[" <<izone << "]: " << *data7 << endl;

    }

  }
  
  void initData() {

    //cout << "Lsp: np: " << np << endl;
    assert(lsp == NULL);
    
    //cout << "np: " << np << " ncv*2: " << ncv*2 << " lp_flag: " << lp_flag << endl;
    //MPI_Pause("is this the same");
    
    lsp = new MyLsp[np];
    my_state = new MyState[ncv];
    my_i = new int[ncv];
    my_d = new double[ncv];
    my_d3 = new double[ncv_g2][3];
    my_n_fa = new double[nfa][3];

  }

  void testNaturalNbrInterp() {

    cout << "testNaturalNbrInterp()" << endl;

    double * phi = new double[ncv];
    const double exact_grad[3] = { 1.1234, 2.1243, 1.5433 };
    FOR_ICV phi[icv] = DOT_PRODUCT(x_vv[icv],exact_grad);
    //FOR_ICV phi[icv] = DOT_PRODUCT(x_cv[icv],exact_grad);
    
    // build gradient...
    
    double (*phi_grad)[3] = new double[ncv][3];
    
    FOR_ICV FOR_I3 phi_grad[icv][i] = 0.0;
    FOR_IFA {
      const int icv0 = cvofa[ifa][0];
      const int icv1 = cvofa[ifa][1];
      
      double x_mid[3]; FOR_I3 x_mid[i] = 0.5*(x_vv[icv0][i]+x_vv[icv1][i]);
      const double dx[3] = DIFF(x_fa[ifa],x_mid);
      if (MAG(dx) > 1.0E-10) {
        cout << "GOT LARGISH MAG(dx): " << MAG(dx) << endl;
        cout << "icv: x_vv: " << icv0 << " " << COUT_VEC(x_vv[icv0]) << " " << icv1 << " " << COUT_VEC(x_vv[icv1]) << endl;
        getchar();
      }
      


      FOR_I3 phi_grad[icv0][i] += n_fa[ifa][i]*phi[icv1];
      FOR_I3 phi_grad[icv1][i] -= n_fa[ifa][i]*phi[icv0];
    }
    
    FOR_ICV {
      FOR_I3 phi_grad[icv][i] /= 2.0*vol_cv[icv];
      // convert to an error:
      FOR_I3 phi_grad[icv][i] -= exact_grad[i];
      cout << "got phi_grad: " << COUT_VEC(phi_grad[icv]) << endl;
    }

    FILE * fp = fopen("check.dat","w");
    FOR_ICV {
      if ((fabs(x_vv[icv][0])<0.9)&&(fabs(x_vv[icv][1])<0.9))
        fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le\n",
                x_vv[icv][0],x_vv[icv][1],x_vv[icv][2],
                phi[icv],
                phi_grad[icv][0],phi_grad[icv][1],phi_grad[icv][2]);
    }
    fclose(fp);

    cout << "take a look at check.dat" << endl; 

    getchar();


  }
  
  void run() {

    //buildDelaunayDual();
    //testNaturalNbrInterp();

    int8 ncolor = getIntParam("COLOR_CVS",-1);
    double* color_cv_double = NULL;
    if (ncolor >= 0) {
      int8* color_cv = new int8[ncv_g];
      Multigrid::colorCvsPadt(color_cv,ncolor,x_cv,icv_global,ncv,ncv_global);
      updateCvData(color_cv);
      Multigrid::splitOrphanedColors(color_cv,ncolor,cvocv_i,cvocv_v,icv_global,
          rbi_g,NULL,ncv,ncv_g,ncv_g,ncv_global); // NULL because we are using compact grid only
      color_cv_double = new double[ncv];
      FOR_ICV {
        // randomize it to spread out color for viz
        const double x = sin(double(color_cv[icv]))*43758.5453123;
        color_cv_double[icv] = x-floor(x);
        //color_cv_double[icv] = (double)color_cv[icv];
      }
      delete[] color_cv;
      registerCvData(color_cv_double,"color",NO_READWRITE_DATA);
    }



    double* A_cv = new double[cvocv_i[ncv]];
    buildCvLaplacian(A_cv);
    
    //Multigrid::AlgebraicCoarseGrid acg;
    //acg.init(vol_cv,x_cv,A_cv,cvocv_i,cvocv_v,icv_global,rbi_g,NULL,ncv,ncv_g,ncv_g,ncv_global,1,8,true);

    if (cvocv2_i == NULL) {
      assert(cvocv2_v == NULL);
      buildCvocv2();
    }
    double* A_cv2 = new double[cvocv2_i[ncv]];
    buildCv2Laplacian(A_cv2);

    //Multigrid::AlgebraicCoarseGrid acg2;
    //acg2.init(vol_cv,x_cv,A_cv2,cvocv2_i,cvocv2_v,icv_global,rbi_g,rbi_g2,ncv,ncv_g,ncv_g2,ncv_global,1,8,true);

    double* phi = new double[ncv_g2];
    FOR_ICV phi[icv] = 0.0; //double(rand())/double(RAND_MAX);
    updateCv2Data(phi);
    double* rhs = new double[ncv];
    //FOR_ICV rhs[icv] = sin(2.0*M_PI*x_cv[icv][0])*sin(2.0*M_PI*x_cv[icv][1]);
    FOR_ICV rhs[icv] = sin(2.0*M_PI*x_cv[icv][0]);
    const double zero = getDoubleParam("ZERO",1.0E-6);
    const double maxcycle = getDoubleParam("MAXCYCLE",100);
    const int nsmooth = getIntParam("NSMOOTH",10);
    const string smoother = getStringParam("SMOOTHER","JACOBI");
    const int maxiter_coarsest = getIntParam("MAXITER_COARSEST",-1);
    const string solver_coarsest = getStringParam("SOLVER_COARSEST",smoother);
    const double relax = getDoubleParam("RELAX",0.7); // jacobi needs something around this
    const double agglomeration_factor = getDoubleParam("AGGLOMERATION_FACTOR",8.0); // 2^3 = 8
    const bool split_orphaned_colors = getBoolParam("SPLIT_ORPHANED_COLORS",true); 
    const int ncg = getIntParam("NCG",1); 
    Multigrid::AlgebraicMultigrid am;
    //am.init(vol_cv,x_cv,A_cv2,cvocv2_i,cvocv2_v,icv_global,rbi_g,rbi_g2,ncv,ncv_g,ncv_g2,ncv_global,ncg,
    //    agglomeration_factor,split_orphaned_colors);
    //am.update(A_cv2);
    am.init(vol_cv,x_cv,A_cv,cvocv_i,cvocv_v,icv_global,rbi_g,NULL,ncv,ncv_g,ncv_g,ncv_global,ncg,
        agglomeration_factor,split_orphaned_colors);
    am.update(A_cv);
    am.solve(phi,rhs,zero,maxcycle,nsmooth,smoother,relax,maxiter_coarsest,solver_coarsest,true);

    delete[] rhs;
    delete[] phi;
    delete[] A_cv;
    delete[] A_cv2;

    /*
    CoarseGrid cg0(this);
    cg0.initCoarseGrid(1,8,true);
    SampleSolver cs0(0); 
    cs0.initFromCoarseGrid(&cg0);
    CoarseGrid cg1(&cs0); 
    cg1.initCoarseGrid(1,8,true);
    SampleSolver cs1(1); 
    cs1.initFromCoarseGrid(&cg1);
    */
    
    if (checkParam("INTERACTIVE"))
      kfr.setHold(true);
    
    for (int i = 1; i < 5; ++i)
      processStep(i);
    
    DELETE(color_cv_double);
    
    flushProbes();
    flushImages();
    flushFwhSurfaces();

    return;

    cout << "particle location routine" << endl;
	
    {
      
    ensureCvAdt();

    double bbox[6];
    getBoundingBox(bbox);
    
    const int np = 1000;
    double (*xp)[3] = new double[np][3];

    if (mpi_rank == 0) {
      for (int ip = 0; ip < np; ++ip) {
	FOR_I3 {
	  double wgt = double(rand())/double(RAND_MAX);
	  xp[ip][i] = wgt*bbox[2*i] + (1.0-wgt)*bbox[2*i+1];
	}
      }
    }
    MPI_Bcast(xp,np*3,MPI_DOUBLE,0,mpi_comm);
    
    vector<int> candidateVec;
    
    int * flag = new int[np];
    for (int ip = 0; ip < np; ++ip) {
      assert(candidateVec.empty());
      cvAdt->buildListForPoint(candidateVec,xp[ip]);
      cout << "RANK: " << mpi_rank << " got candidateVec.size(): " << candidateVec.size() << endl;
      flag[ip] = 0;
      for (int ii = 0; ii < candidateVec.size(); ++ii) {
	const int icv = candidateVec[ii];
	if (bfocv_i[icv+1]-bfocv_i[icv] > 0) {
	  // this is a cv with boundaries, so we check its distance to its
	  // closest face...
	  int ibf_closest = -1;
	  double d2_closest;
	  for (int boc = bfocv_i[icv]; boc != bfocv_i[icv+1]; ++boc) {
	    const int ibf = bfocv_v[boc];
	    // this bf should have a node loop...
	    int ino1 = noobf_v[noobf_i[ibf+1]-1];
	    for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
	      const int ino0 = ino1;
	      ino1 = noobf_v[nob];
	      // join this edge to the bf centroid and compute distance...
	      const double d2 = MiscUtils::getPointToTriDist2(xp[ip],x_bf[ibf],x_no[ino0],x_no[ino1]);
	      if ((ibf_closest == -1)||(d2 < d2_closest)) {
		ibf_closest = ibf;
		d2_closest = d2;
	      }
	    }
	  }
	  assert(ibf_closest != -1);
	}
	else {
	  // totally internal cv...
	  //const double d2 = DIST2(xp[ip],x_vv[icv]);
	  //...
	}
      }
      candidateVec.clear();
      MPI_Pause("OK");
    }


    MPI_Pause("done");
    throw(0);
    
    }

    if (checkParam("TEST_CYL_X")) {
      
      double (*x_cv_tmp)[3] = new double[ncv_g2][3];
      FOR_ICV FOR_I3 x_cv_tmp[icv][i] = x_cv[icv][i];
      for (int icv = ncv; icv < ncv_g2; ++icv) FOR_I3 x_cv_tmp[icv][i] = 1.123E+20;
      updateCv2Data(x_cv_tmp,REPLACE_TRANSLATE_DATA);
      for (int icv = ncv; icv < ncv_g2; ++icv) FOR_I3 assert(x_cv_tmp[icv][i] != 1.123E+20);
      for (int icv = ncv; icv < ncv_g; ++icv) FOR_I3 assert(x_cv_tmp[icv][i] == x_cv[icv][i]);
      
      MPI_Pause("step 1 done");
      
      // check that the extended face proximity is resonable...
      
      double *dn = new double[max(nef,nfa)];
      FOR_IFA {
	const int icv0 = cvofa[ifa][0];
	const int icv1 = cvofa[ifa][1];
	dn[ifa] = DIST(x_cv_tmp[icv0],x_cv_tmp[icv1]);
      }
      MiscUtils::dumpRange(dn,nfa,"dn - compact");	
      
      FOR_IEF {
	const int icv0 = cvoef[ief][0];
	const int icv1 = cvoef[ief][1];
	dn[ief] = DIST(x_cv_tmp[icv0],x_cv_tmp[icv1]);
      }
      MiscUtils::dumpRange(dn,nef,"dn - extended");
      delete[] dn;
      
      MiscUtils::dumpRange(r_vv,ncv,"r_vv");

      MPI_Pause("step 2 done");
      
      // make sure the x points have been rotated properly...
      
      double (*v_cv_tmp)[3] = new double[ncv_g2][3];
      FOR_ICV {
	const double x_cv_mag = MAG(x_cv[icv]);
	assert(x_cv_mag > 0.0);
	FOR_I3 v_cv_tmp[icv][i] = 0.1*x_cv[icv][i]/x_cv_mag; // put a unit outward-pointing normal 
      }
      for (int icv = ncv; icv < ncv_g2; ++icv) FOR_I3 v_cv_tmp[icv][i] = 1.1E+20;
      updateCv2Data(v_cv_tmp,REPLACE_ROTATE_DATA);
      for (int icv = ncv; icv < ncv_g2; ++icv) FOR_I3 assert(v_cv_tmp[icv][i] != 1.1E+20);
      
      double * tmp = new double[ncv_g2];
      for (int icv = 0; icv < ncv_g2; ++icv) {
	const double x_cv_mag = MAG(x_cv_tmp[icv]);
	tmp[icv] = DOT_PRODUCT(x_cv_tmp[icv],v_cv_tmp[icv])/x_cv_mag;
      }
      MiscUtils::dumpRange(tmp,ncv_g2,"x dot v");
      delete[] tmp;
      
      MPI_Pause("HOW WAS THAT?");

    }
    
    int lp_flag = 0;
    if (lp_flag == 0) {

      assert(np == 0);
      FOR_ICV {
	np += 2 + getIcvGlobal(icv)%3;
      }
      
      assert(lsp);
      delete[] lsp;
      lsp = new MyLsp[np];
      
      int ip1 = 0;
      FOR_ICV {
        const int ip0 = ip1;
        ip1 += 2 + getIcvGlobal(icv)%3;
        for (int ip = ip0; ip < ip1; ++ip) {
          lsp[ip].icv = icv;
          lsp[ip].a = double(getIcvGlobal(icv));
          lsp[ip].c = int(getIcvGlobal(icv));
          FOR_I3 lsp[ip].xp[i] = x_cv[icv][i];
          FOR_I3 lsp[ip].xp0[i] = x_vv[icv][i];
        }
      }
      assert(ip1 == np);

      FOR_ICV {
        my_state[icv].i = getIcvGlobal(icv);
        my_state[icv].d = (double)getIcvGlobal(icv);
        FOR_I3 my_state[icv].d3[i] = (double)(3*getIcvGlobal(icv)+i);
        my_i[icv] = getIcvGlobal(icv);
        my_d[icv] = (double)getIcvGlobal(icv);
        FOR_I3 my_d3[icv][i] = (double)(3*getIcvGlobal(icv)+i);
      }
      
    }
    else {

      /*
      {
        char filename[128];
        sprintf(filename,"x.%01d.dat",mpi_rank);
        FILE * fp = fopen(filename,"w");
        int ip1 = 0;
        FOR_ICV {
          const int ip0 = ip1;
          ip1 += 2 + getIcvGlobal(icv)%3;
          for (int ip = ip0; ip < ip1; ++ip) {
            //fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %18.15le %18.15le\n",
            //        x_vv[icv][0],x_vv[icv][1],x_vv[icv][2],
            //        lsp[ip].xp0[0],lsp[ip].xp0[1],lsp[ip].xp0[2]);
            fprintf(fp,"%18.15le %18.15le\n",lsp[ip].a,double(getIcvGlobal(icv)));
          }
        }
        fclose(fp);
      }
      */

      int ip1 = 0;
      FOR_ICV {
        const int ip0 = ip1;
        ip1 += 2 + getIcvGlobal(icv)%3;
        for (int ip = ip0; ip < ip1; ++ip) {
          assert(lsp[ip].icv == icv);
	  assert(lsp[ip].a == double(getIcvGlobal(icv)));
	  assert(lsp[ip].c == int(getIcvGlobal(icv)));
          FOR_I3 assert(lsp[ip].xp[i] == x_cv[icv][i]);
          FOR_I3 assert(lsp[ip].xp0[i] == x_vv[icv][i]);
        }
      }
      assert(ip1 == np);

      //MPI_Pause("DONE!");
      FOR_ICV {
        assert(my_state[icv].i == getIcvGlobal(icv));
        assert(my_state[icv].d == (double)getIcvGlobal(icv));
        FOR_I3 assert(my_state[icv].d3[i] == (double)(3*getIcvGlobal(icv)+i));
        assert(my_i[icv] == getIcvGlobal(icv));
        assert(my_d[icv] == (double)getIcvGlobal(icv));
        FOR_I3 assert(my_d3[icv][i] == (double)(3*getIcvGlobal(icv)+i));
      }


    }

    if (checkDataFlag("my_n_fa")) {
      FOR_IFA {
        FOR_I3 {
          if (fabs(my_n_fa[ifa][i]-n_fa[ifa][i]) > 1.0E-14)
            cout << my_n_fa[ifa][i] << " " << n_fa[ifa][i] << endl;
        }
      }
    }
    else {
      FOR_IFA FOR_I3 my_n_fa[ifa][i] = n_fa[ifa][i];
    }

    // test evaluation of registered data...
    if (checkParam("TEST_EVAL")) {
      while(1) {
        try {
          char c_expression[100]; 
          int size_expression = 0;
          if (mpi_rank == 0) {
            cout << "\n\nEnter expression to evaluate (exit ends loop): " << endl;
            string expression0;
            cin >> expression0;
            cout << "Evaluating expression \"" << expression0 << "\"..." << endl;
            assert(expression0.size()+1 < 100); 
            strcpy(c_expression,expression0.c_str());
            size_expression = expression0.size()+1;
          }
	  //         MPI_Barrier(mpi_comm); // need this???
          MPI_Bcast(&size_expression,1,MPI_INT,0,mpi_comm);
          MPI_Bcast(c_expression,size_expression,MPI_CHAR,0,mpi_comm);
          string expression(c_expression);
          if (expression == "exit")
            break;
          CtiRegister::CtiData * data = CtiRegister::getCtiData(expression);
          assert(data != NULL);
          if (mpi_rank == 0) {
            cout << "Got data: " << *data << "\n" << endl;
            CtiRegister::dumpCurrentData();
          }
        }
        catch(int ierr) {
          if (mpi_rank == 0) cout << "caught error: " << ierr << endl;
        }
        catch(...) {
          if (mpi_rank == 0) cout << "caught unknown error" << endl;
        }
      }
    }
      

    // now write lsp to the restart file...
    
    //StaticSolver::writeLpResult(1);
    
    if (checkParam("TEST_CYL_X")) {

      if ( mpi_rank == 0 ) { 
        cout << " ===================================  " << endl;
        cout << "  === starting tests ===              " << endl;
        cout << " ===================================  " << endl;
      }
      
      assert( state == NULL);
      state = new SampleState[ncv_g2];
      
      // populate the state with some known values..
      for (int icv = 0; icv < ncv; ++icv) { 
        state[icv].it = 1;
        state[icv].a  = 2.0;
        state[icv].b  = 3.1;
        state[icv].c  = 4.2;
        state[icv].d =  5.3;
        
        for (int i =0; i < 3; ++i) 
          state[icv].x[i] = x_cv[icv][i];
      }
      
      // set the ghosts to invalid values.. 
      for (int icv = ncv; icv < ncv_g2; ++icv) { 
        state[icv].it = -1;
        state[icv].a  =  1.0e+20;
        state[icv].b  =  1.0e+20;
        state[icv].c  =  1.0e+20;
        state[icv].d  =  1.0e+20;
        
        for (int i = 0; i < 3; ++i) 
          state[icv].x[i] = 1.0e+20; 
      } 
      
      // check the level 1 ghost packings.. 
      if ( mpi_rank == 0 ) 
        cout << " > checking level 1 class update ... " ;
      updateCvClass<SampleState,char>(state,REPLACE_ROTATE_DATA);
      
      const double TOL = 1.0e-16;
      for (int icv = ncv; icv < ncv_g; ++icv) { 
        assert( state[icv].it == 1); 
        assert( abs(state[icv].a - 2.0) < TOL);
        assert( abs(state[icv].b - 3.1) < TOL);
        assert( abs(state[icv].c - 4.2) < TOL);
        assert( abs(state[icv].d - 5.3) < TOL);
      }
      
      if ( mpi_rank == 0 ) 
        cout << " OK. " << endl;

      // and check some facts of the x replace_rotate_data

      double * dn = new double[max(nef,nfa)];
      for (int ifa = 0; ifa < nfa; ++ifa) { 
        const int icv0 = cvofa[ifa][0];
        const int icv1 = cvofa[ifa][1];
        dn[ifa]        = DIST(state[icv1].x,state[icv0].x);
      }
      
      MiscUtils::dumpRange(dn,nfa, "dn -- compact check via state");

      // populate the state with some different values..
      for (int icv = 0; icv < ncv; ++icv) { 
        state[icv].it = 2;
        state[icv].a  = 6.7;
        state[icv].b  = 8.1;
        state[icv].c  = 9.3;
        state[icv].d =  0.4;
        
        for (int i =0; i < 3; ++i) 
          state[icv].x[i] = x_cv[icv][i];
      }
      
      // set the ghosts to invalid values.. 
      for (int icv = ncv; icv < ncv_g2; ++icv) { 
        state[icv].it = -1;
        state[icv].a  =  1.0e+20;
        state[icv].b  =  1.0e+20;
        state[icv].c  =  1.0e+20;
        state[icv].d  =  1.0e+20;
        
        for (int i =0; i < 3; ++i)
          state[icv].x[i] = 1.0e+20;
      } 
      
      // check the level 1 ghost packings.. 
      if ( mpi_rank == 0 ) 
        cout << " > checking level 1+2 class update ... " ;
      updateCv2Class<SampleState,char>(state, REPLACE_ROTATE_DATA);
      
      for (int icv = ncv; icv < ncv_g2; ++icv) { 
        assert( state[icv].it == 2); 
        assert( abs(state[icv].a - 6.7) < TOL);
        assert( abs(state[icv].b - 8.1) < TOL);
        assert( abs(state[icv].c - 9.3) < TOL);
        assert( abs(state[icv].d - 0.4) < TOL);
      }
      
      if ( mpi_rank == 0 ) 
        cout << " OK. " << endl;

      // check the spacing on the extended faces
      for (int ief = 0; ief < nef; ++ief) { 
        const int icv0 = cvoef[ief][0];
        const int icv1 = cvoef[ief][1];
        dn[ief]        = DIST(state[icv1].x,state[icv0].x);
      }
 
      MiscUtils::dumpRange(dn,nef, "dn -- extended check via state");

      delete[] dn;
      
      double (*v_cv_tmp)[3] = new double[ncv_g2][3];
      FOR_ICV {
	const double x_cv_mag = MAG(x_cv[icv]);
	assert(x_cv_mag > 0.0);
	FOR_I3 v_cv_tmp[icv][i] = 0.1*x_cv[icv][i]/x_cv_mag; // put a unit outward-pointing normal 
      }
      for (int icv = ncv; icv < ncv_g2; ++icv) FOR_I3 v_cv_tmp[icv][i] = 1.1E+20;
      updateCv2Data(v_cv_tmp,REPLACE_ROTATE_DATA);
      for (int icv = ncv; icv < ncv_g2; ++icv) FOR_I3 assert(v_cv_tmp[icv][i] != 1.1E+20);
      
      double * tmp = new double[ncv_g2];
      for (int icv = 0; icv < ncv_g2; ++icv) {
	const double x_cv_mag = MAG(state[icv].x);
	tmp[icv] = DOT_PRODUCT(state[icv].x,v_cv_tmp[icv])/x_cv_mag;
      }
      MiscUtils::dumpRange(tmp,ncv_g2,"x dot v");
      delete[] tmp;
 
    }

    for (int i = 0; i < 3; ++i)
      processStep(i);
    
    if ( checkParam("TEST_DN33_UPDATE")) { 

      double (*T)[3][3] = new double[ncv_g][3][3];
      int *icv_local    = new int[ncv_g];

      // pack the data .. 

      for (int icv = 0; icv < ncv; ++icv) { 
        icv_local[icv] = icv;
        for (int i =0; i < 3; ++i) { 
          for (int j = 0; j < 3; ++j) { 
            T[icv][i][j] = double(9*icv + 3*i + j);
          }
        }
      }

      for (int icv = ncv; icv < ncv_g; ++icv) { 
        icv_local[icv] = -1;
        for (int i =0; i < 3; ++i) 
          for (int j = 0; j < 3; ++j) 
            T[icv][i][j] = 1.0e+20;
      }

      updateCvData(icv_local);
      //updateCvData(T,REPLACE_ROTATE_DATA); // test only true if there is no cyl periodicity.. 
      updateCvData(T,REPLACE_DATA);

      for (int icv = 0; icv < ncv_g; ++icv) { 
        const int ii = icv_local[icv];
        for (int i =0; i < 3; ++i) { 
          for (int j = 0; j < 3; ++j) { 
            assert( abs(T[icv][i][j] - double(9*ii+3*i+j)) < 1.0e-12);
          }
        }
      }

      delete[] T;
      delete[] icv_local;

      if ( mpi_rank == 0 ) 
        cout << " passed test dn33 update " << endl;
    }

    //writeResult(0);
    
  } 
  
  virtual ~SampleSolver() {
    FOR_IZONE(bcVec) {
      delete bcVec[izone];
    }
    DELETE(state);
    DELETE(lsp);
    DELETE(my_i);
    DELETE(my_d);
    DELETE(my_d3);
    DELETE(my_n_fa);
    DELETE(my_state);
  }
  void initBoundaryConditions();

  void init() { 

    StaticSolver::init(); // extended faces
    //StaticSolver::init(INIT_COMPACT_FACES|INIT_CV_GRAD);

    if (mpi_rank == 0) 
      logger->setKillFilename("killcharles");
    
  }

};

// specific Bc implementations need to know about the solver...

class SampleWallBc : public SampleSolverBaseBc {
private:
  SampleSolver * solver;
public:
  SampleWallBc(Param * param,BfZone& bfZone,SampleSolver * s) : SampleSolverBaseBc(bfZone), solver(s) {
    COUT1(" > SampleWallBc: " << getName());
    // set the bc cost here like this...
    bfZone.lb_cost = 100;
  }
  void calcRhs() {
    COUT1("SampleWallBc::calcRhs()");
  }
};

class SampleInletBc : public SampleSolverBaseBc {
private:
  SampleSolver * solver;
public:
  SampleInletBc(Param * param,BfZone& bfZone,SampleSolver * s) : SampleSolverBaseBc(bfZone), solver(s) {
    COUT1(" > SampleInletBc: " << getName());
    // you can set the bc cost...
    bfZone.lb_cost = 1000;
  }
  void calcRhs() {
    COUT1("SampleInletBc::calcRhs()");
  }
};

void SampleSolver::initBoundaryConditions() {

  StaticSolver::initBoundaryConditions();

  COUT1("SampleSolver::initBoundaryConditions()");

  int ierr = 0;
  FOR_IZONE(bfZoneVec) {
    
    Param * param = getParam(bfZoneVec[izone].getName());
    if (param == NULL) {
      if (mpi_rank == 0) cout << "Error: BC not found for zone \"" << bfZoneVec[izone].getName() << "\"" << endl;
      ierr = -1;
    }
    else {
      
      const string token = param->getString();
      if (token == "WALL") {
	bcVec.push_back(new SampleWallBc(param,bfZoneVec[izone],this));
      }
      else if (token == "INLET") {
	bcVec.push_back(new SampleInletBc(param,bfZoneVec[izone],this));
      }
      else {
	if (mpi_rank == 0) cout << "Error: BC token \"" << token << "\" not recognized for zone \"" << bfZoneVec[izone].getName() << "\"" << endl;
	ierr = -1;
      }
	
    }

  }

  if (ierr != 0) {
    CWARN("bc problem in initBoundaryConditions()");
  }
  
}

#include "FluentReader.hpp"

void processMsh() {
  
  // cap transformations, surface and tet writing...
  // Feb 2020
  
  FluentMsh * msh = new FluentMsh("metal_2.msh");
  
  // now transform...
  for (int ino = 0; ino < msh->nno; ++ino) {
    msh->x_no[ino][0] = (msh->x_no[ino][0]+132.81981002)*0.0254;
    msh->x_no[ino][1] = (msh->x_no[ino][1]+0.0)*0.0254;
    msh->x_no[ino][2] = (msh->x_no[ino][2]+0.0)*0.0254;
  }
  
  // write the full surface mesh out. Delete zones in surfer to get the mesh to 
  // replace the old cap with... 
  //msh->writeSurfaceSbin("cht.sbin");
  //MPI_Pause("how was that?");
  
  // try to build tets...
  cout << "trying to build tets..." << endl;
  int nte;
  int (*noote)[4] = NULL;
  msh->buildTets(nte,noote);
  cout << " > got nte: " << nte << endl;
  
  double vol_sum = 0.0;
  double vol_min = HUGE_VAL;
  double vol_max = -HUGE_VAL;
  for (int ite = 0; ite < nte; ++ite) {
    // flip first 2 nodes...
    const int tmp = noote[ite][0];
    noote[ite][0] = noote[ite][1];
    noote[ite][1] = tmp;
    // and check volume...
    const int ino0 = noote[ite][0];
    const int ino1 = noote[ite][1];
    const int ino2 = noote[ite][2];
    const int ino3 = noote[ite][3];
    const double vol = SIGNED_TET_VOLUME_6(msh->x_no[ino0],msh->x_no[ino1],msh->x_no[ino2],msh->x_no[ino3])/6.0;
    vol_sum += vol;
    vol_min = min(vol_min,vol);
    vol_max = max(vol_max,vol);
  }
  
  // we require positive volumes...
  cout << " > total volume: " << vol_sum << ", individual tet volume (min,avg,max): " << 
    vol_min << " " << vol_sum/double(nte) << " " << vol_max << endl;
  
  // write the tets to a simple binary file...
  /*
  cout << "writing cap.tet_bin with nno: " << msh->nno << ", nte: " << nte << endl;
  FILE * fp = fopen("cap.tet_bin","wb");
  assert(fp != NULL);
  fwrite(&msh->nno,sizeof(int),1,fp);
  fwrite(&nte,sizeof(int),1,fp);
  fwrite(msh->x_no,sizeof(double),msh->nno*3,fp);
  fwrite(noote,sizeof(int),nte*4,fp);
  fclose(fp);
  */
  
  MPI_Pause("and how was that?");
  
}


#include "SimpleFunc.hpp"

void parseSimpleFunc() {
  
  Param * param = getParam("CHT.RHO");
  assert(param);
  SimpleFunc * rho_func = processSimpleFunc(param);
  assert(rho_func != NULL);
  
  for (int x = 200.0; x < 3000.0; x += 7.1) {
    const double y = rho_func->eval(x);
    cout << "XXXXX " << x << " " << y << endl;
  }
    
  MPI_Pause("how was that?");

}

int main(int argc, char * argv[]) {

  try {
    
    CTI_Init(argc,argv,"testStaticSolver.in");

    {
      
      //parseSimpleFunc();
      //processMsh();
      //testOctree();
      //MPI_Run_check();
      
      SampleSolver solver;
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

  return(0);

}

