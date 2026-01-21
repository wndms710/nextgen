#ifndef __SINGLEDROP_HPP__
#define __SINGLEDROP_HPP__

#include "CTI.hpp"
#include "Lsp.hpp"
#include "CtiLiquid.hpp"
#include "CtiSolid.hpp"
#include "LspInjector.hpp"

class SingleLspSolver {

protected:
  
  LpClass<LspState> * lpCls;

  CtiMaterial * material;
  CtiLiquid * fuel;
  CtiSolid  * dust;
  bool SOL_LP, LIQ_LP;

  double u_g[3];
  double T_g;
  double rho_g;
  double mu_g;
  double p_g;
  double R_UNIVERSAL; 
  double MW_g;
  double R_g;
  double gamma;
  double Pr_g;

  int check_interval;
  int step;
  int nsteps;
  double time;
  double dt;
  // time integration method EXPLICIT/IMPLICIT
  string time_int;

  // monitoring file
  ofstream fid;

public:

  SingleLspSolver() {

    COUT1("SingleLspSolver()");

    if (mpi_size != 1) {
      CERR("SingleLspSolver should be run with single cpu. Make sure the you are not running in parallel.");
    }

    // create the lp (not allocating at this stage)
    lpCls   = new LpClass<LspState>; // sets the lp to NULL, np and NP_MAX to 0
    material= NULL;
    fuel    = NULL;
    dust    = NULL;
    SOL_LP  = LIQ_LP = false;

    step = 0;
    time = 0;
    nsteps = getIntParam("NSTEPS");
    dt = getDoubleParam("DT");
    check_interval = getIntParam("CHECK_INTERVAL");

    {
      Param * param = getParam("LSP.MATERIAL");
      if (param == NULL) CERR("could not find param LSP.MATERIAL");
      string token = param->getString();
      if (token == "LIQUID") {
        fuel = newCtiLiquid(param->getString(1),101325);
        material = fuel;
        LIQ_LP = true;
        IF_RANK0 fuel->info(300,350,10);
        //MPI_Pause("check it out");
      } else if (token == "SOLID") {
        dust = newCtiSolid(param->getString(1));
        IF_RANK0 dust->info(300,350,10);
        //MPI_Pause("check it out");
        material = dust;
        SOL_LP = true;
      } else {
        CERR("token for LSP.MATERIAL was not recognized!");
      }
    }

    {
      bool got_u, got_rho, got_T, got_mu, got_gamma;
      got_u = got_rho = got_T = got_mu = got_gamma = false;
      Param * param = getParam("GAS_PROP");
      int iarg = 0;
      while (iarg < param->size()) {
        string token = param->getString(iarg++);
        if (token == "U") {
          FOR_I3 u_g[i] = param->getDouble(iarg++);
          got_u = true;
        } else if (token == "RHO") {
          rho_g = param->getDouble(iarg++);
          got_rho = true;
        } else if (token == "T") {
          T_g = param->getDouble(iarg++);
          got_T = true;
        } else if (token == "MU") {
          mu_g = param->getDouble(iarg++);
          got_mu = true;
        } else if (token == "GAMMA") {
          gamma = param->getDouble(iarg++);
          got_gamma = true;
        } else {
          CERR("Unknown parameter for GAS_PROP");
        }
      }

      if (!got_u) 
        CERR("U missing in GAS_PROP");
      if (!got_rho) 
        CERR("RHO missing in GAS_PROP");
      if (!got_T) 
        CERR("T missing in GAS_PROP");
      if (!got_mu) 
        CERR("MU missing in GAS_PROP");
      if (!got_gamma) 
        CERR("GAMMA missing in GAS_PROP");
      
      R_UNIVERSAL = 8314.472; // J/kmol/K
      MW_g = 29.; // kg/kmol
      Pr_g = 0.7; // Prandtl number
      R_g = R_UNIVERSAL / MW_g;
      p_g = rho_g * R_g * T_g;
    }

    {
      Param * param = getParam("TIME_INT");
      if (!param) {
        CERR("SingleLspSolver requires param TIME_INT <stirng>");
      } else {
        string token = param->getString();
        if ((token == "EXPLICIT") || (token == "IMPLICIT")) {
          time_int = param->getString();
        } else {
          CERR("SingleLspSolver requires param TIME_INT <stirng>");
        }
      }
    }

    fid.open("regression.fit");    
  }

  virtual ~SingleLspSolver() {

    COUT1("~SingleLspSolver()");

    if (material != NULL) {material = NULL;}
    if (fuel != NULL) {delete fuel; fuel = NULL;}
    if (dust != NULL) {delete dust; dust = NULL;}
    if (lpCls != NULL) {delete lpCls; lpCls = NULL;}

    if (fid.is_open()) fid.close();
  }

  void init() {

    initLsp();

  }

  void initLsp() {

    lpCls->resize(1);

    FOR_I3 lpCls->lp[0].xp[i] = 0.0;
    FOR_I3 lpCls->lp[0].xp0[i] = 0.0;
    lpCls->lp[0].flag = KEEP;
    lpCls->lp[0].icv = 0;
    lpCls->lp[0].npar = 1;
    lpCls->lp[0].tbu = 0.0;
    lpCls->lp[0].k = 0.0;
    lpCls->lp[0].kt = 0.0;
    lpCls->lp[0].kb = 0.0;
    lpCls->lp[0].ktb = 0.0;
    lpCls->lp[0].mp0 = 0.0;

    bool got_U, got_T, got_D;
    got_U = got_T = got_D = false;

    // setting the particles in each cv 
    IF_RANK0 cout << "initLsp()" << endl;
    Param * param = getParam("LSP.INIT");
    if (param == NULL) CERR("input file needs LSP.INIT");
    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "U") {
        FOR_I3 lpCls->lp[0].up[i] = param->getDouble(iarg++);
        got_U = true;
      } else if (token == "T") {
        lpCls->lp[0].Tp = param->getDouble(iarg++);
        got_T = true;
      } else if (token == "D") { 
        lpCls->lp[0].dp = param->getDouble(iarg++);
        got_D = true;
      } else {
        CERR("Unknown parameter for LSP_INIT");
      }
    }

    if (!got_U) 
      CERR("U missing in LSP.INIT");
    if (!got_T)
      CERR("T missing in LSP.INIT");
    if (!got_D)
      CERR("D missing in LSP.INIT");

    lpCls->lp[0].mp = M_PI/6.*pow(lpCls->lp[0].dp,3)*material->calcRho(lpCls->lp[0].Tp);

  }

  void run() {

    // put a header for the solid particle test case
    if (SOL_LP) {
      double tau_p = dust->calcRho(lpCls->lp[0].Tp) * pow(lpCls->lp[0].dp,2) / (18. * mu_g);
      double up0 = MAG(lpCls->lp[0].up);
      fid << "# " << "tau: " << tau_p << " up0: " << up0 << " ug: " << MAG(u_g) << endl;
    }
    report();

    int done = 0;
    while (done == 0) { 

      ++step;
      
      if (step%check_interval == 0) {
        cout <<
          "\n----------------------------------------------------------\n" <<
          " starting step: " << step << " time: " << time+dt << " dt: " << dt <<
          "\n----------------------------------------------------------" << endl;
      }

      advanceSolution(); // time advanced in here with solution to support time-dependent evaluations

      report();

      // additional diagonistics potentially can go in here... 
      // write snapshots, images, tecplot, probes, fwh, boundary surfaces... 
      if ((nsteps >= 0) && (step >= nsteps)) { 
        done = 1;
        COUT1(" > reach NSTEPS " << nsteps << ". Stopping.");
      }
    }

  }

  void report () {
    if (step%check_interval == 0) {
      // some reports
      LspState* lsp = &lpCls->lp[0];
      cout << "time: " << time << " xp: " << COUT_VEC(lsp->xp) << endl;
      cout << "time: " << time << " up: " << COUT_VEC(lsp->up) << endl;
      cout << "time: " << time << " Tp: " << (lsp->Tp) << endl;
      cout << "time: " << time << " dp: " << (lsp->dp) << endl;
      cout << "time: " << time << " dp^2: " << (pow(lsp->dp,2)) << endl;
      cout << "time: " << time << " mp: " << (lsp->mp) << endl;

      if (SOL_LP)
        fid << time << "   " << lsp->up[0] << endl;
      else if (LIQ_LP)
        fid << time << "   " << (pow(lsp->dp,2) * 1e6) << endl;
        
    }
  }
  void advanceSolution() {

    LspRhs * rhslsp0 = new LspRhs [lpCls->size()];
    LspRhs * rhslsp1 = new LspRhs [lpCls->size()];
    LspRhs * rhslsp2 = new LspRhs [lpCls->size()];
    LspRhs * rhslsp_arr[3] = {rhslsp0, rhslsp1, rhslsp2};

    LspRhs * rhslsp_d = new LspRhs [lpCls->size()]; // each substep sets to zero first (so one variable is enough)

    if (time_int == "EXPLICIT") {
      // rk3 update...
      calcRhsLsp_ex(rhslsp0);
      //addLspSrc(rhs0);
      rk3StepLsp(erk_wgt3[0],rhslsp_arr,1);

      time += dt;
      updateLspPrimitiveData();

      // stage 2...
      calcRhsLsp_ex(rhslsp1);
      //addLspSrc(rhs1);
      rk3StepLsp(erk_wgt3[1],rhslsp_arr,2);

      time -= 0.5*dt;
      updateLspPrimitiveData();

      // stage 3...
      calcRhsLsp_ex(rhslsp2);
      //addLspSrc(rhs1);
      rk3StepLsp(erk_wgt3[2],rhslsp_arr,3);

      time += 0.5*dt;
      updateLspPrimitiveData();
    }
    else if (time_int == "IMPLICIT") {
      // rk3 update...
      calcRhsLsp_im(rhslsp0,rhslsp_d);
      //addLspSrc(rhs0);
      rk3StepLsp_im(erk_wgt3[0], rhslsp_arr, rhslsp_d, 1);

      time += dt;
      updateLspPrimitiveData();

      // stage 2...
      calcRhsLsp_im(rhslsp1,rhslsp_d);
      //addLspSrc(rhs1);
      rk3StepLsp_im(erk_wgt3[1], rhslsp_arr, rhslsp_d, 2);

      time -= 0.5*dt;
      updateLspPrimitiveData();

      // stage 3... 
      calcRhsLsp_im(rhslsp2,rhslsp_d);
      //addLspSrc(rhs2);
      rk3StepLsp_im(erk_wgt3[2], rhslsp_arr, rhslsp_d, 3);

      time += 0.5*dt;
      updateLspPrimitiveData();

    }

    delete [] rhslsp0;
    delete [] rhslsp1;
    delete [] rhslsp2;
    delete [] rhslsp_d;
  }

  void calcRhsLsp_ex(LspRhs* rhs) {
    if (SOL_LP)
      calcRhsLsp_ex_sol(rhs);
    else if (LIQ_LP)
      calcRhsLsp_ex_liq(rhs);
    else
      CERR("LP type is not specified");
  }

  void calcRhsLsp_im(LspRhs* rhs, LspRhs* rhs_d) {
    if (SOL_LP)
      calcRhsLsp_im_sol(rhs, rhs_d);
    else if (LIQ_LP)
      calcRhsLsp_im_liq(rhs, rhs_d);
    else
      CERR("LP type is not specified");
  }

  // rk3 explicit time advancement
  void rk3StepLsp(const double * rk_wgt, LspRhs ** rhslsp_arr,const int rkstep) {
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      LspRhs rhs_agg;
      rhs_agg.zero();
      for (int irk = 0; irk < rkstep; ++irk) {
        double wgt = dt*rk_wgt[irk];
        rhs_agg.axpy(wgt,rhslsp_arr[irk][ip]);
      }
      lpCls->lp[ip].add(rhs_agg);
    }
  }
    
  // rk3 implicit time advancement
  void rk3StepLsp_im(const double * rk_wgt, LspRhs ** rhslsp_arr, LspRhs* rhs_d, const int& rkstep) {
    for (int ip=0; ip < lpCls->size(); ++ip) {
      // build the new state
      LspRhs rhs_agg;
      rhs_agg.zero();
      for (int irk=0; irk < rkstep; ++irk) {
        double wgt = dt*rk_wgt[irk];
        rhs_agg.axpy(wgt,rhslsp_arr[irk][ip]);
      }   
      const double temp = dt*rk_wgt[rkstep-1];
      lpCls->lp[ip].add_im(rhs_agg,rhs_d[ip],temp);
      // correct the rhs of previous state
      FOR_I3 rhslsp_arr[rkstep-1][ip].xp[i] += rhs_d[ip].xp[i]*lpCls->lp[ip].xp[i];
      FOR_I3 rhslsp_arr[rkstep-1][ip].up[i] += rhs_d[ip].up[i]*lpCls->lp[ip].up[i];
      rhslsp_arr[rkstep-1][ip].Tp += rhs_d[ip].Tp*lpCls->lp[ip].Tp;
      rhslsp_arr[rkstep-1][ip].mp += rhs_d[ip].mp*lpCls->lp[ip].mp;
    }
  }

  void updateLspPrimitiveData() {
    // should not set xp0 here. it should be set only in the begining of the 
    // time step not after each substep, otherwise the collision detection fails if
    // we set it equal to xp at the end of the time step
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      updateDp(ip);
    }
  }

  void updateDp(const int &ip) {
    lpCls->lp[ip].dp = pow(lpCls->lp[ip].mp/(M_PI/6.0*lpCls->lp[ip].npar*material->calcRho(lpCls->lp[ip].Tp)) , 1./3.);
  }

  void calcRhsLsp_ex_sol(LspRhs* rhs) {

    for (int ip = 0; ip < lpCls->size(); ++ip) {
      // initialize rhs with zero
      rhs[ip].zero();

      //#############################################
      // npar: we compute the Rhs of one droplet, 
      //       the effect of npar is taken into acount
      //       in energy transfer and mass transfer to 
      //       the gas only
      //#############################################
      // make sure mass is not zero otherwise:
      // rhs and rhs_d would remain zero for ip -> no change in xp, mp
      // set k and kt to zero -> no change in up, Tp (coupled to gas)
      // continue
      //#############################################
      if ((lpCls->lp[ip].flag == TRASH) || (lpCls->lp[ip].mp <= 0.0)) {
        lpCls->lp[ip].flag = TRASH; // make sure it will be trashed
        continue;
      }

      // gas properties
      const double rhog = rho_g;
      const double mug  = mu_g;
      const double cpg  = R_g*gamma/(gamma-1.0);
      const double kg   = cpg * mug / Pr_g; // Pr: cp * mu / k
      const double ug[3]= {u_g[0], u_g[1], u_g[2]};
      const double Tg   = T_g;

      // compute liquid properties at parcel temperature
      const double rhop = dust->calcRho(lpCls->lp[ip].Tp);
      const double Dp   = lpCls->lp[ip].dp;
      assert(Dp > 0.0);

      double du[3];
      FOR_I3 du[i] = ug[i] - lpCls->lp[ip].up[i];
      double dus = MAG(du);
      // compute non-dimensional and model parameters
      const double tau = rhop*Dp*Dp/(18.*mug);
      const double Rep = rhog*Dp*dus/mug;
      const double Prg = mug*cpg/kg;
      const double Nu  = 2. + 0.552*sqrt(Rep)*pow(Prg,1./3.);

      // drag law from  Crowe et al, Multiphase Flows w/ Droplets & Particles, 1998
      //const double f1 = 1.0+0.15*pow(Rep,0.687) + 0.0175*Rep/(1.0+4.25e4*pow(Rep,-1.16));
      // drag law from White, eqn 3-225
      // CD = 24.0/Rep + 6.0/(1.0+sqrt(Rep)) + 0.4 and f1 = Rep/24.0*CD but for numerical reasons when Rep=0 we add a correction to f1
      //const double CD_c = 6.0/(1.0+sqrt(Rep)) + 0.4;
      //const double f1 =  1.0 + Rep/24.0*CD_c;
      const double f1 = 1.0; // only for verificatiion step where exact solution needs to be compared

      const double f2 = 1.0; // for no evaporative case

      ////////////////
      // setting rhs
      ////////////////
      // set xp rhs
      FOR_I3 rhs[ip].xp[i] = lpCls->lp[ip].up[i];
      //FOR_I3 rhs[ip].xp[i] = 0.0;

      // set up rhs 
      FOR_I3 rhs[ip].up[i] = f1/tau*(ug[i] - lpCls->lp[ip].up[i]);;

      // energy equation
      double const coef = f2*Nu*(cpg/dust->calcCp(lpCls->lp[ip].Tp))/(3.0*Prg*tau);
      rhs[ip].Tp = coef * (Tg - lpCls->lp[ip].Tp);

      // mass transfer equation
      rhs[ip].mp = 0.0;

      // monitor if needed
      if (step%check_interval==0) {
      //  cout << " > ip: "    << ip   << endl <<
      //          " > mug = "  << mug  << endl <<
      //          " > rhog = " << rhog << endl << 
      //          " > ug = " << COUT_VEC(ug) << endl;

      //  cout << " > rhop = " << rhop << endl <<
      //          " > Dp = " << Dp << endl;

        cout << " > tau = " << tau << endl;
  
      //  cout << " > Rep = " << Rep << endl;

      //  cout << " > f1 = " << f1 << endl;

      }

    }
  }

  void calcRhsLsp_ex_liq(LspRhs* rhs) {

    //FOR_ICV lspSrc[icv].zero(); // only used for adding evaporated mass to Y
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      // initialize rhs with zero
      rhs[ip].zero();

      //#############################################
      // npar: we compute the Rhs of one droplet, 
      //       the effect of npar is taken into acount
      //       in energy transfer and mass transfer to 
      //       the gas only
      //#############################################
      // make sure mass is not zero otherwise:
      // rhs and rhs_d would remain zero for ip -> no change in xp, mp
      // set k and kt to zero -> no change in up, Tp (coupled to gas)
      // continue
      //#############################################
      if ((lpCls->lp[ip].flag == TRASH) || (lpCls->lp[ip].mp <= 0.0)) {
        lpCls->lp[ip].k = 0.0;
        lpCls->lp[ip].kt = 0.0;
        lpCls->lp[ip].kb = 0.0;
        lpCls->lp[ip].ktb = 0.0;
        lpCls->lp[ip].flag = TRASH; // make sure it will be trashed
        continue;
      }

      // gas properties
      //const int icv     = lpCls->lp[ip].icv;
      const double Tr   = fuel->Tboil;;
      //const double Tg   = T[icv];
      const double Pg   = p_g;
      const double rhog = rho_g;
      //const double mug  = mu_g;

      const double mug  = calcGasMu(Tr);
      const double kg   = calcGasK(Tr);
      const double cpg  = calcGasCp(Tr, mug, kg);

      //const double cpg  = R_g*gamma/(gamma-1.0);
      //const double kg   = cpg * mug / Pr_g; // Pr: cp * mu / k
      //const double Yfg  = max(0.0,cv[icv].Y); 
      const double Yfg  = 0.0; 
      const double ug[3]= {u_g[0], u_g[1], u_g[2]};
      const double MWg  = MW_g; // kg/m^3 for air

      // compute liquid properties at parcel temperature
      const double rhof = fuel->calcRho(lpCls->lp[ip].Tp);  // liquid density
      const double cpf  = fuel->calcCp(lpCls->lp[ip].Tp);   // liquid heat capacity
      const double Lv   = fuel->calcHvap(lpCls->lp[ip].Tp); // latent heat of evaporation
      const double Dp   = lpCls->lp[ip].dp;
      assert(Dp > 0.0);

      // equilibrium fuel composition at droplet surface
      const double nearOne = 0.9999999999;
      const double Xfseq   = min(nearOne,fuel->calcPv(lpCls->lp[ip].Tp)/Pg);
      const double Yfseq   = Xfseq/(Xfseq+(1.0-Xfseq)*MWg/fuel->MW);
      const double BM      = (Yfseq - Yfg)/(1.0 - Yfseq);            // Spalding mass transfer number

      // compute reference prop
      //const double Yfr = (2.*Yfseq + Yfg)/3.;  // 1/3 rule
      const double Yfr = Yfg;
      const double mur = Yfr*fuel->calcMuv(Tr) + (1.-Yfr)*mug;
      const double rhor= Yfr*fuel->calcRhov(Tr) + (1.-Yfr)*rhog;
      const double kr  = Yfr*fuel->calcKv(Tr) + (1.-Yfr)*kg;
      const double cpr = Yfr*fuel->calcCpv(Tr) + (1.-Yfr)*cpg;

      double du[3];
      FOR_I3 du[i] = ug[i] - lpCls->lp[ip].up[i];
      double dus = MAG(du);
      // compute non-dimensional and model parameters
      const double tau = rhof*Dp*Dp/(18.*mur);
      const double Rep = rhor*Dp*dus/mur;
      const double Prg = mur*cpr/kr;
      const double Scg = mur/(rhor*fuel->calcDv(Tr));
      const double Nu  = 2. + 0.552*sqrt(Rep)*pow(Prg,1./3.);
      const double Sh  = 2. + 0.552*sqrt(Rep)*pow(Scg,1./3.);
      const double Lk  = mur*sqrt(2.0*M_PI*lpCls->lp[ip].Tp*fuel->R_UNIVERSAL/fuel->MW)/(Scg*Pg);

      // Newton-iteration for evaporation rate mpdot
      const double F1 = (Sh/Scg/3.0)*(lpCls->lp[ip].mp/tau); // mdot + F1*ln(1+BM) = 0, see eq (4)
      const double F2 = -1.5*Prg*tau/lpCls->lp[ip].mp;        // beta = F2*mdot, see eq (17)
      double beta, Xfsneq, Yfsneq, BMneq, Fm;
      double mdot  = -1.0E-08;
      double mdot0 =  0.0;
      double Fm0    =  mdot0 + F1*log(1.0+BM); // initialize with equilib value
      double eps    =  1.0E-15;
      double resid  =  eps+1.0;
      int    iter   =  0;
      while (resid > eps) {
        iter += 1;
        beta   = F2*mdot;
        Xfsneq = Xfseq - (Lk/Dp*2.0)*beta;
        Yfsneq = Xfsneq/(Xfsneq+(1.0-Xfsneq)*MWg/fuel->MW);
        BMneq  = max(-nearOne,(Yfsneq-Yfg)/(1.0-Yfsneq));
        Fm     = mdot + F1*log(1.0+BMneq);
        double tmp = mdot;
        if (fabs(Fm-Fm0) < eps ) break;
        mdot = mdot - Fm/(Fm-Fm0)*(mdot-mdot0);
        resid = sqrt(pow(Fm-Fm0,2));
        Fm0 = Fm;
        mdot0 = tmp;
        if (iter > 20) cout << "ERROR: LSP Newton iteration did not converge !!!" << endl;
      }
      //cout << " > iter = " << iter << endl;
      beta = F2*mdot;

      // Dont let possitive mdot
      mdot = min(mdot, 0.0);
      assert(mdot<=0);

      const double mdot_over_mp = mdot/lpCls->lp[ip].mp;
      // drag law from  Crowe et al, Multiphase Flows w/ Droplets & Particles, 1998
      //const double f1 = 1.0+0.15*pow(Rep,0.687) + 0.0175*Rep/(1.0+4.25e4*pow(Rep,-1.16));
      // drag law from White, eqn 3-225
      // CD = 24.0/Rep + 6.0/(1.0+sqrt(Rep)) + 0.4 and f1 = Rep/24.0*CD but for numerical reasons when Rep=0 we add a correction to f1
      //const double CD_c = 6.0/(1.0+sqrt(Rep)) + 0.4;
      //const double f1 =  1.0 + Rep/24.0*CD_c;
      double f2 = beta/(exp(beta)-1.0);
      //double f2 = 1.0; // only for non-evaporating case
      if (beta < 1e-16)
        f2 = 1.0;

      // set the xp
      FOR_I3 rhs[ip].xp[i] = lpCls->lp[ip].up[i];

      // set the up
      //FOR_I3 rhs[ip].up[i] = f1/tau * (u_g[i] - lpCls->lp[ip].up[i]);

      // for evaporation the droplet should be still
      FOR_I3 rhs[ip].up[i] = 0.0;

      // set the Tp 
      const double theta1 = cpr/cpf;
      rhs[ip].Tp = f2*Nu*theta1/(3.*Prg*tau)*(T_g - lpCls->lp[ip].Tp) + Lv/cpf*mdot_over_mp;

      // mass transfer equation
      rhs[ip].mp = mdot;

      // Y is the mass fraction of vapor in the gas (vapor mass/ carrier mass)
      // the lspSrc.rhoY is the evaporated vapor mass per time (integrated in vol)
      //
      // YOU SHOULD ADD THIS WHEN YOU HAVE RHOY XXX
      //
      //lspSrc[icv].rhoY += -mdot*lpCls->lp[ip].npar;

      // monitor if needed:
      //if (step%check_interval=0 && monitor_bool) {
      //  cout << " > ip: "    << ip   << endl <<
      //          " > Tp: "    << lpCls->lp[ip].Tp << endl;;
      //  cout << " > Yfg = "  << Yfg  << endl <<
      //          " > mug = "  << mug  << endl <<
      //          " > rhog = " << rhog << endl <<
      //          " > kg = "   << kg   << endl <<
      //          " > Pg = "   << Pg   << endl <<
      //          " > Tg = "   << Tg   << endl <<
      //          " > cpg = "  << cpg  << endl;

      //  cout << " > rhof = " << rhof << endl <<
      //          " > cpf = " << cpf << endl <<
      //          " > Lv = " << Lv << endl <<
      //          " > Dp = " << Dp << endl;

      //  cout << " > Xfseq = " << Xfseq << endl <<
      //          " > Yfseq = " << Yfseq << endl <<
      //          " > BM = " << BM << endl <<
      //          " > Tr " << Tr << endl;

      //  cout << " > Yfr = " << Yfr << endl <<
      //          " > mur = " << mur << endl <<
      //          " > rhor = " << rhor << endl <<
      //          " > kr = " << kr << endl <<
      //          " > cpr = " << cpr << endl;

      //  cout << " > tau = " << tau << endl <<
      //          " > Rep = " << Rep << endl <<
      //          " > Prg = " << Prg << endl <<
      //          " > Scg = " << Scg << endl <<
      //          " > Nu = " << Nu << endl <<
      //          " > Sh = " << Sh << endl <<
      //          " > Lk = " << Lk << endl;

      //  cout << "mdot: " << mdot << endl;

      //  cout << " > f1 = " << f1 << endl <<
      //          " > f2 = " << f2 << endl <<
      //          " > beta = " << beta << endl;
      //}
    }
  }

  void calcRhsLsp_im_sol(LspRhs* rhs, LspRhs* rhs_d) {

    for (int ip = 0; ip < lpCls->size(); ++ip) {
      // initialize rhs with zero
      rhs[ip].zero();

      //#############################################
      // npar: we compute the Rhs of one droplet, 
      //       the effect of npar is taken into acount
      //       in energy transfer and mass transfer to 
      //       the gas only
      //#############################################
      // make sure mass is not zero otherwise:
      // rhs and rhs_d would remain zero for ip -> no change in xp, mp
      // set k and kt to zero -> no change in up, Tp (coupled to gas)
      // continue
      //#############################################
      if ((lpCls->lp[ip].flag == TRASH) || (lpCls->lp[ip].mp <= 0.0)) {
        lpCls->lp[ip].flag = TRASH; // make sure it will be trashed
        continue;
      }

      // gas properties
      const double rhog = rho_g;
      const double mug  = mu_g;
      const double cpg  = R_g*gamma/(gamma-1.0);
      const double kg   = cpg * mug / Pr_g; // Pr: cp * mu / k
      const double ug[3]= {u_g[0], u_g[1], u_g[2]};
      const double Tg   = T_g;

      // compute liquid properties at parcel temperature
      const double rhop = dust->calcRho(lpCls->lp[ip].Tp);
      const double Dp   = lpCls->lp[ip].dp;
      assert(Dp > 0.0);

      double du[3];
      FOR_I3 du[i] = ug[i] - lpCls->lp[ip].up[i];
      double dus = MAG(du);
      // compute non-dimensional and model parameters
      const double tau = rhop*Dp*Dp/(18.*mug);
      const double Rep = rhog*Dp*dus/mug;
      const double Prg = mug*cpg/kg;
      const double Nu  = 2. + 0.552*sqrt(Rep)*pow(Prg,1./3.);

      // drag law from  Crowe et al, Multiphase Flows w/ Droplets & Particles, 1998
      //const double f1 = 1.0+0.15*pow(Rep,0.687) + 0.0175*Rep/(1.0+4.25e4*pow(Rep,-1.16));
      // drag law from White, eqn 3-225
      // CD = 24.0/Rep + 6.0/(1.0+sqrt(Rep)) + 0.4 and f1 = Rep/24.0*CD but for numerical reasons when Rep=0 we add a correction to f1
      //const double CD_c = 6.0/(1.0+sqrt(Rep)) + 0.4;
      //const double f1 =  1.0 + Rep/24.0*CD_c;
      const double f1 = 1.0; // only for verificatiion step where exact solution needs to be compared

      const double f2 = 1.0;

      ///////////////////////////////
      // setting rhs and rhs_d
      // rhs_total = rhs + rhs_d*var
      ///////////////////////////////
      // set xp rhs
      FOR_I3 rhs[ip].xp[i] = lpCls->lp[ip].up[i];
      FOR_I3 rhs_d[ip].xp[i] = 0.0;

      // set up rhs 
      FOR_I3 rhs[ip].up[i] = f1/tau*(ug[i]);
      FOR_I3 rhs_d[ip].up[i] = -f1/tau;

      // energy equation
      //rhs[ip].Tp = f2*Nu*(cpg/dust->calcCp(lpCls->lp[ip].Tp))/(3.0*Prg*tau) * (Tg - lpCls->lp[ip].Tp);
      double const coef = f2*Nu*(cpg/dust->calcCp(lpCls->lp[ip].Tp))/(3.0*Prg*tau);
      rhs[ip].Tp = coef * (Tg);
      rhs_d[ip].Tp = -coef;

      // mass transfer equation
      // mass does not change
      rhs[ip].mp = 0.0;
      rhs_d[ip].mp = 0.0;

      // monitor if needed
      if (step%check_interval==0) {
      //  cout << " > ip: "    << ip   << endl <<
      //          " > mug = "  << mug  << endl <<
      //          " > rhog = " << rhog << endl << 
      //          " > ug = " << COUT_VEC(ug) << endl;

      //  cout << " > rhop = " << rhop << endl <<
      //          " > Dp = " << Dp << endl;

        cout << " > tau = " << tau << endl;

      //  cout << " > Rep = " << Rep << endl;

      //  cout << " > f1 = " << f1 << endl;

      }

    }

  }

  void calcRhsLsp_im_liq(LspRhs* rhs, LspRhs* rhs_d) {

    // IMPORTANT NOTE: the velocity is not advanced since the initial droplet velocity in the 
    // validation case is 0 and not changed.

    //FOR_ICV lspSrc[icv].zero(); // only used for adding evaporated mass to Y
    for (int ip = 0; ip < lpCls->size(); ++ip) {
      // initialize rhs with zero
      rhs[ip].zero();
      rhs_d[ip].zero();

      //#############################################
      // npar: we compute the Rhs of one droplet, 
      //       the effect of npar is taken into acount
      //       in energy transfer and mass transfer to 
      //       the gas only
      //#############################################
      // make sure mass is not zero otherwise:
      // rhs and rhs_d would remain zero for ip -> no change in xp, mp
      // set k and kt to zero -> no change in up, Tp (coupled to gas)
      // continue
      //#############################################
      if ((lpCls->lp[ip].flag == TRASH) || (lpCls->lp[ip].mp <= 0.0)) {
        lpCls->lp[ip].k = 0.0;
        lpCls->lp[ip].kt = 0.0;
        lpCls->lp[ip].kb = 0.0;
        lpCls->lp[ip].ktb = 0.0;
        lpCls->lp[ip].flag = TRASH; // make sure it will be trashed
        continue;
      }

      // gas properties
      //const int icv     = lpCls->lp[ip].icv;
      const double Tr   = fuel->Tboil;
      //const double Tg   = T[icv];
      const double Pg   = p_g;
      const double rhog = rho_g;
      //const double mug  = mu_g;

      const double mug  = calcGasMu(Tr);
      const double kg   = calcGasK(Tr);
      const double cpg  = calcGasCp(Tr, mug, kg);

      //const double cpg  = R_g*gamma/(gamma-1.0);
      //const double kg   = cpg * mug / Pr_g; // Pr: cp * mu / k
      //const double Yfg  = max(0.0,cv[icv].Y); 
      const double Yfg  = 0.0; 
      const double ug[3]= {u_g[0], u_g[1], u_g[2]};
      const double MWg  = MW_g; // kg/m^3 for air

      // compute liquid properties at parcel temperature
      const double rhof = fuel->calcRho(lpCls->lp[ip].Tp);  // liquid density
      const double cpf  = fuel->calcCp(lpCls->lp[ip].Tp);   // liquid heat capacity
      const double Lv   = fuel->calcHvap(lpCls->lp[ip].Tp); // latent heat of evaporation
      const double Dp   = lpCls->lp[ip].dp;
      assert(Dp > 0.0);

      // equilibrium fuel composition at droplet surface
      const double nearOne = 0.9999999999;
      const double Xfseq   = min(nearOne,fuel->calcPv(lpCls->lp[ip].Tp)/Pg);
      const double Yfseq   = Xfseq/(Xfseq+(1.0-Xfseq)*MWg/fuel->MW);
      const double BM      = (Yfseq - Yfg)/(1.0 - Yfseq);            // Spalding mass transfer number

      // compute reference prop
      //const double Yfr = (2.*Yfseq + Yfg)/3.;  // 1/3 rule
      const double Yfr = Yfg;
      const double mur = Yfr*fuel->calcMuv(Tr) + (1.-Yfr)*mug;
      const double rhor= Yfr*fuel->calcRhov(Tr) + (1.-Yfr)*rhog;
      const double kr  = Yfr*fuel->calcKv(Tr) + (1.-Yfr)*kg;
      const double cpr = Yfr*fuel->calcCpv(Tr) + (1.-Yfr)*cpg;

      double du[3];
      FOR_I3 du[i] = ug[i] - lpCls->lp[ip].up[i];
      double dus = MAG(du);
      // compute non-dimensional and model parameters
      const double tau = rhof*Dp*Dp/(18.*mur);
      const double Rep = rhor*Dp*dus/mur;
      const double Prg = mur*cpr/kr;
      const double Scg = mur/(rhor*fuel->calcDv(Tr));
      const double Nu  = 2. + 0.552*sqrt(Rep)*pow(Prg,1./3.);
      const double Sh  = 2. + 0.552*sqrt(Rep)*pow(Scg,1./3.);
      const double Lk  = mur*sqrt(2.0*M_PI*lpCls->lp[ip].Tp*fuel->R_UNIVERSAL/fuel->MW)/(Scg*Pg);

      // Newton-iteration for evaporation rate mpdot
      const double F1 = (Sh/Scg/3.0)*(lpCls->lp[ip].mp/tau); // mdot + F1*ln(1+BM) = 0, see eq (4)
      const double F2 = -1.5*Prg*tau/lpCls->lp[ip].mp;        // beta = F2*mdot, see eq (17)
      double beta, Xfsneq, Yfsneq, BMneq, Fm;
      double mdot  = -1.0E-08;
      double mdot0 =  0.0;
      double Fm0    =  mdot0 + F1*log(1.0+BM); // initialize with equilib value
      double eps    =  1.0E-15;
      double resid  =  eps+1.0;
      int    iter   =  0;
      while (resid > eps) {
        iter += 1;
        beta   = F2*mdot;
        Xfsneq = Xfseq - (Lk/Dp*2.0)*beta;
        Yfsneq = Xfsneq/(Xfsneq+(1.0-Xfsneq)*MWg/fuel->MW);
        BMneq  = max(-nearOne,(Yfsneq-Yfg)/(1.0-Yfsneq));
        Fm     = mdot + F1*log(1.0+BMneq);
        double tmp = mdot;
        if (fabs(Fm-Fm0) < eps ) break;
        mdot = mdot - Fm/(Fm-Fm0)*(mdot-mdot0);
        resid = sqrt(pow(Fm-Fm0,2));
        Fm0 = Fm;
        mdot0 = tmp;
        if (iter > 20) cout << "ERROR: LSP Newton iteration did not converge !!!" << endl;
      }
      //cout << " > iter = " << iter << endl;
      beta = F2*mdot;

      // Dont let possitive mdot
      mdot = min(mdot, 0.0);
      assert(mdot<=0);

      const double mdot_over_mp = mdot/lpCls->lp[ip].mp;
      // drag law from  Crowe et al, Multiphase Flows w/ Droplets & Particles, 1998
      //const double f1 = 1.0+0.15*pow(Rep,0.687) + 0.0175*Rep/(1.0+4.25e4*pow(Rep,-1.16));
      // drag law from White, eqn 3-225
      // CD = 24.0/Rep + 6.0/(1.0+sqrt(Rep)) + 0.4 and f1 = Rep/24.0*CD but for numerical reasons when Rep=0 we add a correction to f1
      const double CD_c = 6.0/(1.0+sqrt(Rep)) + 0.4;
      const double f1 =  1.0 + Rep/24.0*CD_c;
      double f2 = beta/(exp(beta)-1.0);
      //double f2 = 1.0; // only for non-evaporating case
      if (beta < 1e-16)
        f2 = 1.0;

      ///////////////////////////////
      // setting rhs and rhs_d
      // rhs_total = rhs + rhs_d*var
      ///////////////////////////////
      // set the xp
      FOR_I3 rhs[ip].xp[i] = lpCls->lp[ip].up[i];
      FOR_I3 rhs_d[ip].xp[i] = 0.0;

      lpCls->lp[ip].k = f1/tau;
      lpCls->lp[ip].kt = f2*Nu*(cpr/cpf)/(3.0*Prg*tau);
      lpCls->lp[ip].kb =  0.0;
      lpCls->lp[ip].ktb = 0.0;

      assert(lpCls->lp[ip].k>0);
      assert(lpCls->lp[ip].kt>0);
      // up would be solved implicitly coupled with ug later
      // In this test case we cannot couple to the gas so we use implicit formulation independent of the gas
      // set the up
      //FOR_I3 rhs[ip].up[i] = lpCls->lp[ip].k * u_g[i];
      //FOR_I3 rhs_d[ip].up[i] = -lpCls->lp[ip].k;

      // for evaporation the droplet should be still
      FOR_I3 rhs[ip].up[i] = 0.0;
      FOR_I3 rhs_d[ip].up[i] = 0.0;

      // energy equation
      //const double theta1 = cpr/cpf;
      //rhs[ip].Tp = f2*Nu*theta1/(3.*Prg*tau)*(Tg - P[ip].Tp) + Lv/cpf*mdot/P[ip].mp;
      //rhs[ip].Tp = f2*Nu*theta1/(3.*Prg*tau)*(Tg) + Lv/cpf*mdot_over_mp;
      //rhs_d[ip].Tp = - f2*Nu*theta1/(3.*Prg*tau);

      //#############################################
      // only evaporation latent heat is treated explicitly
      // coupling with the gas is treated implicitly and will
      // be solved coupled to gas later
      //#############################################
      //rhs[ip].Tp = 0.0;
      // In this test case we cannot couple to the gas so we use implicit formulation independent of the gas
      rhs[ip].Tp = lpCls->lp[ip].kt * T_g + Lv/cpf*mdot_over_mp;
      rhs_d[ip].Tp = -lpCls->lp[ip].kt;

      //#############################################
      // if you like to test a zero evaporation case, make sure
      // in addition to zeroing mass and temperature source terms
      // you make f2=1. This is the limit for non-evaporative case
      // otherwise your kt coeff is not correct
      //#############################################

      // mass transfer equation
      rhs[ip].mp = 0.0;
      rhs_d[ip].mp = mdot_over_mp;
      //rhs[ip].mp = mdot;
      //rhs_d[ip].mp = 0.0;

      // Y is the mass fraction of vapor in the gas (vapor mass/ carrier mass)
      // the lspSrc.rhoY is the evaporated vapor mass per time (integrated in vol)
      //
      // YOU SHOULD ADD THIS WHEN YOU HAVE RHOY XXX
      //
      //lspSrc[icv].rhoY += -mdot*lpCls->lp[ip].npar;

      // monitor if needed:
      //if (step%check_interval=0 && monitor_bool) {
      //  cout << " > ip: "    << ip   << endl <<
      //          " > Tp: "    << lpCls->lp[ip].Tp << endl;;
      //  cout << " > Yfg = "  << Yfg  << endl <<
      //          " > mug = "  << mug  << endl <<
      //          " > rhog = " << rhog << endl <<
      //          " > kg = "   << kg   << endl <<
      //          " > Pg = "   << Pg   << endl <<
      //          " > Tg = "   << Tg   << endl <<
      //          " > cpg = "  << cpg  << endl;

      //  cout << " > rhof = " << rhof << endl <<
      //          " > cpf = " << cpf << endl <<
      //          " > Lv = " << Lv << endl <<
      //          " > Dp = " << Dp << endl;

      //  cout << " > Xfseq = " << Xfseq << endl <<
      //          " > Yfseq = " << Yfseq << endl <<
      //          " > BM = " << BM << endl <<
      //          " > Tr " << Tr << endl;

      //  cout << " > Yfr = " << Yfr << endl <<
      //          " > mur = " << mur << endl <<
      //          " > rhor = " << rhor << endl <<
      //          " > kr = " << kr << endl <<
      //          " > cpr = " << cpr << endl;

      //  cout << " > tau = " << tau << endl <<
      //          " > Rep = " << Rep << endl <<
      //          " > Prg = " << Prg << endl <<
      //          " > Scg = " << Scg << endl <<
      //          " > Nu = " << Nu << endl <<
      //          " > Sh = " << Sh << endl <<
      //          " > Lk = " << Lk << endl;

      //  cout << "mdot: " << mdot << endl;

      //  cout << " > f1 = " << f1 << endl <<
      //          " > f2 = " << f2 << endl <<
      //          " > beta = " << beta << endl;
      //}
    }
  }

  // compute gas properties dynamically
  double calcGasMu(const double& Tr) {
    double mug = 6.109e-06 +  4.604e-08*Tr -  1.051e-11*Tr*Tr; // kg/m/s
    return mug;
  }

  double calcGasK(const double& Tr) {
    double kg = 3.227e-03 + 8.3894e-05*Tr - 1.9858e-08*Tr*Tr; // W/m/K
    return kg;
  }

  double calcGasCp(const double& Tr, const double& mug, const double& kg) {
    double cpg;
    if (Tr > 600) {
      cpg = (0.647 + 5.5e-05*Tr)*kg/mug;
    }
    else {
      cpg = (0.815 - 4.958e-04*Tr + 4.514e-07*Tr*Tr)*kg/mug; // J/kg/K
    }
    return cpg;
  }

};

#endif
