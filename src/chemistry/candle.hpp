/*!
 * flamelet generation using cantera
 * Author: Amirreza Saghafian 1/11/16
 */

#ifndef CANDLE_HPP
#define CANDLE_HPP

#define EPSC 1.0e-6

//#include "Domain1DLe.h"
//#include "cantera/oneD/Sim1D.h"
#include "Sim1DLe.h"
//#include "cantera/oneD/Inlet1D.h"
#include "Inlet1DLe.h"
//#include "cantera/oneD/StFlow.h"
#include "StFlowLe.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"
#include "cantera/kinetics.h"
#include <ctime>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include "myMem.h"

using namespace CanteraLe;
//using CanteraLe::Domain1D;
using Cantera::IdealGasMix;
using std::cout;
using std::endl;
using std::string;

class Flame {

  protected:
    IdealGasMix *mix;
//    FreeFlame *f;
    Sim1D *flame;
    Transport *trans;

    double *Le;

    int flowdomain;

    doublereal p;
    doublereal Tu;

    string mechanism;
    string transport;
    bool SORET;
    IdealGasMix *mix_i;
    IdealGasMix *mix_j;
    size_t nsp;

    string fuel_compo;
    vector_fp Xf;
    double Tf;
    string oxidizer_compo;
    vector_fp Xox;
    double Tox;

    // grid refinement parameters
    double ratio;
    double slope;
    double curve;
    double prune;

    size_t i_O2;

    int loglevel;

  public:

    Flame(const string fuel, string ox, const double pu, const double Tfuel, const double Toxidizer, const string mech, const string tr,
          const double Lewis = 1.0, const int lglvl = 0, const bool soret = false, const double rat = 5.0, const double slo = 0.05, const double cur = 0.1, const double pru = 0.03) {

      mechanism = mech;
      transport = tr;

      if (mechanism=="GRI3.0") mechanism = "gri30.xml";
      if (ox=="AIR") ox = "O2:1 N2:3.76";
      if (transport=="MIXTURE_AVERAGED")    transport = "Mix";
      else if (transport=="MULTICOMPONENT") transport = "Multi";
      SORET = soret;//false;

      mix   = new IdealGasMix(mechanism);
      mix_i = new IdealGasMix(mechanism);
      mix_j = new IdealGasMix(mechanism);
      trans = newTransportMgr(transport, mix_i); // XXX changed 1 to 0
      nsp   = mix->nSpecies();

      p   = pu;
      Tf  = Tfuel;
      Tox = Toxidizer;

      fuel_compo = fuel;
//      printf("========= fuel compo ============\n");
      Xf = getMixVecFromString(fuel_compo);

      oxidizer_compo = ox;
//      printf("========= oxidizer compo ============\n");
      Xox = getMixVecFromString(oxidizer_compo);

      // grid refinement parameters
      ratio = rat;
      slope = slo;
      curve = cur;
      prune = pru;

      i_O2 = mix_i->speciesIndex("O2");

      loglevel = lglvl;

      Le = new double[nsp];
      setLewis(Lewis);
      if (Lewis > 0 ) cout << " >>>>>>>>>> Constant Lewis Number Assumption: Le = " << Lewis << " <<<<<<<<<<" << endl;
    }

//    Flame() {
//      mechanism = "gri30.xml";
//      transport = "Mix";
//
//      mix   = new IdealGasMix(mechanism);
//      mix_i = new IdealGasMix(mechanism);
//      trans = newTransportMgr(transport, mix_i, 0);
//      nsp = mix->nSpecies();
//
//      p   = 1.0e5;             // constant background pressure
//      Tu = 298.0;             // inlet temperature
//      fuel_compo     = "CH4";
//      oxidizer_compo = "Air";
//
//      ratio = 10.0;
//      slope = 0.40;//0.2;
//      curve = 0.04;//0.02;
//
//      i_O2 = mix_i->speciesIndex("O2");
//
//      loglevel = 1;
//    }

    ~Flame() {
      delete flame;
      delete mix;
      delete mix_i;
      delete mix_j;
      delete[] Le;
    }

    doublereal calcO2Stoich() {

      bool CinMech = 0, HinMech = 0, OinMech = 0;
      for (size_t i = 0; i < mix_i->nElements(); ++i) { 
        if (mix_i->elementName(i) == "C") {
          CinMech = true;
        } else if (mix_i->elementName(i) == "H") {
          HinMech = true;
        } else if (mix_i->elementName(i) == "O") {
          OinMech = true;
        }
      }

      double N_O2 = 0.0;
      for (int isp=0; isp<nsp; isp++) {
        double nC = 0.0, nH = 0.0, nO = 0.0; 
        if (CinMech) nC = mix_i->nAtoms(isp, mix_i->elementIndex("C"));
        if (HinMech) nH = mix_i->nAtoms(isp, mix_i->elementIndex("H"));
        if (OinMech) nO = mix_i->nAtoms(isp, mix_i->elementIndex("O"));
        N_O2 += Xf[isp] * ( nC + 0.25*nH - 0.5*nO ); 

        //N_O2 += Xf[isp] * ( 1.00 * mix_i->nAtoms(isp, mix_i->elementIndex("C")) +
        //                    0.25 * mix_i->nAtoms(isp, mix_i->elementIndex("H")) -
        //                    0.50 * mix_i->nAtoms(isp, mix_i->elementIndex("O")) );
      }
      return N_O2;
    }

    double calcZStoich() {
      mix_i->setState_TPX(Tf, p, &Xf[0]);
      const double mf_st  = mix_i->meanMolecularWeight();

      mix_i->setState_TPX(Tox, p, &Xox[0]);
      const double mox_st = mix_i->meanMolecularWeight() * calcO2Stoich()/Xox[i_O2];

      return mf_st / (mf_st+mox_st);
    }

    double calcf_aStoich() {
      mix_i->setState_TPX(Tf, p, &Xf[0]);
      const double mf_st  = mix_i->meanMolecularWeight();

      mix_i->setState_TPX(Tox, p, &Xox[0]);
      const double mox_st = mix_i->meanMolecularWeight() * calcO2Stoich()/Xox[i_O2];

      return mf_st/mox_st;
    }

    vector_fp getMixVecFromString(const string compo) {
      vector_fp X(nsp);
      mix_i->setState_TPX(300.0, 1.0e5, compo);
      for (int isp=0; isp<nsp; isp++) {
        X[isp] = mix_i->moleFraction(isp);
//        printf("species %s has X = %f\n", mix->speciesName(isp).c_str(), X[isp]);
      }

      return X;
    }

    template <class flowType>
    std::vector<doublereal> getVariableVec(flowType *flow, string varName) {
      int np = flow->nPoints();
      std::vector<doublereal> varVec;

      for (int n=0; n<np; n++)
        varVec.push_back( flame->value(flowdomain, flow->componentIndex(varName), n) );

      return varVec;
    }

    template <class flowType>
    vector_fp getMassFracVec(flowType *flow, const int n) {
      vector_fp y(nsp);

      for (size_t k=0; k<nsp; k++) {
        y[k] = flame->value( flowdomain, flow->componentIndex(mix->speciesName(k)), n );
      }

      return y;
    }

    template <class flowType>
    void calcThermoProp(flowType *flow, doublereal *rho, doublereal *mw, doublereal *cp, doublereal *mu, doublereal *lambda, doublereal *h, std::vector<doublereal> T) {
      int np = flow->nPoints();

      for (int n=0; n<np; n++) {
        mix_i->setState_TPY( T[n], p, &(getMassFracVec(flow, n))[0] );
        rho[n]    = mix_i->density();
        mw[n]     = mix_i->meanMolecularWeight();
        cp[n]     = mix_i->cp_mass();
        mu[n]     = trans->viscosity();
        lambda[n] = trans->thermalConductivity();
        h[n]      = mix_i->enthalpy_mass();
      }
    }

//    std::vector<doublereal> calcDensity(std::vector<doublereal> T) {
//      int np = f->nPoints();
//      std::vector<doublereal> rho;
//
//      for (int n=0; n<np; n++) {
//        mix_i->setState_TPY( T[n], p, &(getMassFracVec(f, n))[0] );
//        rho.push_back(mix_i->density());
////        rho[n] = mix_i->density();
//      }
//
//      return rho;
//    }

    double calcSpeciesEnthalpy(const doublereal T, const int isp) {
      vector_fp y(nsp);
      for (int i = 0; i < nsp; ++i) {
        y[i] = 0.0;
      }
      y[isp] = 1.0;
      mix_j->setState_TPY(T, p, &y[0]);
      return mix_j->enthalpy_mass();
    }

    // Computes the heat release rate given temperature and mass fraction vector
    double calcHeatRelease(const doublereal T, const vector_fp y) {
      vector_fp prodRate(nsp);
      mix_i->getNetProductionRates(&prodRate[0]);
      double HeatRelease = 0.0;
      mix_i->setState_TPY(T, p, &y[0]);
      for (int isp = 0; isp < nsp; ++isp) {
        HeatRelease -= prodRate[isp] * mix_i->molecularWeight(isp) * calcSpeciesEnthalpy(T, isp);
      }
      return HeatRelease;
    }

    // Compute entropy given temperature and mass fraction vector
    double calcEntropy(const doublereal T, const vector_fp y) {
      mix_i->setState_TPY(T, p, &y[0]);
      return mix_i->entropy_mass();
    }

    vector_fp calcProdRate_i(doublereal T, vector_fp y) {
      vector_fp prodRate(nsp);
      mix_i->setState_TPY(T, p, &y[0] );
      mix_i->getNetProductionRates(&prodRate[0]);
      return prodRate;
    }

//    void calcProdRate() {
//      T = getVariableVec("T");
//      int np = f->nPoints();
//
//      for (int n = 0; n < np; n++) {
//        vector_fp this_ProdRate = calcProdRate_i(T[n], getMassFracVec(f, n) );
//        for (size_t k=0; k<nsp; k++)
//          ProdRate[n][k] = this_ProdRate[k];
//      }
//
//    }

    template <class flowType>
    double calcFlameThickness(flowType *flow, std::vector<doublereal> T ) {
      int np = flow->nPoints();
      double dTdy_max = 0.0;

      for (int n=1; n<np; n++) {
        double dTdy = (T[n-1]-T[n]) / (flow->grid(n-1)-flow->grid(n));
        if (dTdy>dTdy_max) dTdy_max = dTdy;
      }

      double delta_f = 384.4e6;
      if (fabs(dTdy_max)>1.0e-10) delta_f = (T[np-1]-T[0]) / dTdy_max;

      return delta_f;
    }

    void writeFlameletList(std::vector<string> flLi) {
      string fileName = "FlameletList.txt";
      FILE* FP = fopen(fileName.c_str(), "w");
      if (!FP) {
          printf("Failure to open file\n");
          exit(-1);
      }

      for (int i=0; i<flLi.size(); i++)
        fprintf(FP, "%s\n", flLi[i].c_str());

      fclose(FP);

    }

    void setLewis(const double this_Le) {
//      Le.resize(nsp);
      for (int isp=0; isp<nsp; isp++) {
        Le[isp] = this_Le;
      }
    }

    template <class flowType>
    std::vector<double> calcMixtureFraction(flowType *f, const string element) {
      int np = f->nPoints();
      vector_fp Z(np);

      std::vector<doublereal> T = getVariableVec(f, "T");

      mix_i->setState_TPX(Tox, p, &Xox[0]);
      const double Z_elem_ox = mix_i->elementalMassFraction( mix_i->elementIndex(element) );
      mix_i->setState_TPX(Tf, p, &Xf[0]);
      const double Z_elem_f  = mix_i->elementalMassFraction( mix_i->elementIndex(element) );

      for (int ip=0; ip<np; ip++) {
        mix_i->setState_TPY( T[ip], p, &(getMassFracVec(f, ip))[0] );
        double Z_elem = mix_i->elementalMassFraction( mix_i->elementIndex(element) );
        Z[ip] = (Z_elem - Z_elem_ox) / (Z_elem_f - Z_elem_ox);
      }

      return Z;
    }

    template <class flowType>
    void writeFlameletCSV(flowType *f, string fileName) {
      int np = f->nPoints();
      std::vector<doublereal> Tvec, Uvec, Zvec;

      for (int n=0; n<np; n++) {
          Tvec.push_back(flame->value( flowdomain, f->componentIndex("T"), n) );
          Uvec.push_back(flame->value( flowdomain, f->componentIndex("u"), n) );
      }

      Zvec = calcMixtureFraction(f, "H");

      FILE* FP = fopen(fileName.c_str(), "w");
      if (!FP) {
          printf("Failure to open file\n");
          exit(-1);
      }

      fprintf( FP, "VARIABLES = y, Z, prog, T, u, " );
      for (int isp=0; isp< nsp; isp++)
        fprintf( FP, " %s,", mix_i->speciesName(isp).c_str() );
      fprintf(FP, "\n");

      for (int n = 0; n < np; n++) {

        const double C = flame->value( flowdomain, f->componentIndex("CO2"), n) +
                         flame->value( flowdomain, f->componentIndex("CO"), n)  +
                         flame->value( flowdomain, f->componentIndex("H2O"), n) +
                         flame->value( flowdomain, f->componentIndex("H2"), n);

        fprintf(FP," %11.3e, %11.3e, %11.3e, %11.3e, %11.3e",
                f->grid(n), Zvec[n], C, Tvec[n], Uvec[n]);
        for (int isp=0; isp< nsp; isp++)
          fprintf(FP," %11.3e,",
                  flame->value( flowdomain, f->componentIndex(mix_i->speciesName(isp)), n) );
        fprintf(FP, "\n");

      }
      fclose(FP);
    }



};

class Premixed: public Flame {

  private:
    Inlet1D inlet;
    Outlet1D outlet;
    FreeFlameLe *f;

    bool writeCSV;
    bool createInit;
    string initFlamelet;
    bool generateOnlyFlameletRange;

  public:
    Premixed(const string fuel, string ox, const double pu, const double Tfuel, const double Toxidizer, const string mech, const string tr,
             const double Lewis = 1.0, const int lglvl = 0,
             const bool init = true, const string initFile = "init.xml", const bool onlyRange = false,
             const bool csv = false, const bool soret = false,
//             const double rat = 5.0, const double slo = 0.05, const double cur = 0.1, const double pru = 0.03) :
             const double rat = 10.0, const double slo = 0.075, const double cur = 0.15, const double pru = 0.05) :
               Flame(fuel, ox, pu, Tfuel, Toxidizer, mech, tr, Lewis, lglvl, soret, rat, slo, cur, pru),
               writeCSV(csv), createInit(init), initFlamelet(initFile), generateOnlyFlameletRange(onlyRange) {
    }

    ~Premixed() {
      delete f;
    }

    int createFlamelets(const double phi1, const double phi2, const int nfl) {

      try {
        doublereal uin = 0.0;//0.3;  // inlet velocity
        const double phi_init = phi1;

        //mix->setState_TPX(Tu, p, "CH4:1.0, O2:2.0, N2:7.52");
        vector_fp x(nsp);
        doublereal N_O2_st = calcO2Stoich();
        cout << " >> moles of O2 at st for 1 mole of fuel: " << N_O2_st << endl;
        cout << " >> Zst = " << calcZStoich() << "; (f/a)st = " << calcf_aStoich() << endl;

        doublereal rho_in;
        doublereal rho_out;
        doublereal Tad;
        vector_fp yin(nsp);
        vector_fp yout(nsp);
        doublereal mdot;
        doublereal lz = 0.15;

        if (createInit) {
          for (size_t k=0; k<nsp; k++) {
            x[k] = Xf[k]*phi_init + Xox[k]/Xox[i_O2]*N_O2_st;
          }

          setInlet(phi_init);
  //        mix->setState_TPX(Tu, p, &x[0]);
          rho_in = mix->density();
        }

        const double rtol    = 1.0e-8;
        const double atol    = 1.0e-12;
        const double rtol_tr = 1.0e-8;
        const double atol_tr = 1.0e-10;
        doublereal dt = 1.0e-8;  // time step
        size_t nStepsInit = 7;
        size_t nSteps     = 2;
        integer steps[] = {2, 5, 10, 20, 40, 80, 120};

        flowdomain = 1;

        double PHI_INF;
        (p>10e5) ? PHI_INF = 0.96 : PHI_INF = 0.94;

        if (createInit) {
          mix->getMassFractions(&yin[0]);

          try {
            mix->equilibrate("HP");
          }
          catch (CanteraError& err) {
            std::cout << err.what() << std::endl;
          }

          mix->getMassFractions(&yout[0]);
          rho_out = mix->density();
          Tad = mix->temperature();
        }

        clock_t tb = clock();

        //=============  build each domain ========================

        //-------- step 1: create the flow -------------
        f = new FreeFlameLe(mix);

        if (createInit) {
          // create an initial grid
          int nz = 21;
          vector_fp z(nz+1);
          doublereal dz = lz/((doublereal)(nz-1));
          for (int iz=0; iz<nz; iz++) {
              z[iz] = ((doublereal)iz)*dz;
          }
          //add one node onto end of domain to help with zero gradient at outlet
          z[nz] = lz*1.05;
          nz++;

          f->setupGrid(nz, &z[0]);
        }

        // specify the objects to use to compute kinetic rates and
        // transport properties

        std::auto_ptr<Transport> trmix(newTransportMgr("Mix", mix));
        std::auto_ptr<Transport> trmulti(newTransportMgr("Multi", mix));

//        setLewis(1.0);
        f->setTransport(*trmix, Le);
        f->setKinetics(*mix);
        f->setPressure(p);

//        f->setTransientTolerances(rtol_tr, atol_tr);
//        f->setSteadyTolerances(rtol, atol);

        if (createInit) {
          //------- step 2: create the inlet  -----------------------
          inlet.setMoleFractions(&x[0]);
          mdot = uin*rho_in;
          inlet.setMdot(mdot);
          inlet.setTemperature(Tu);

          //------- step 3: create the outlet  ---------------------
  //        Outlet1D outlet;
        }

        //=================== create the container and insert the domains =====

        std::vector<Domain1D*> domains;
        domains.push_back(&inlet);
        domains.push_back(f);
        domains.push_back(&outlet);

        flame = new Sim1D(domains);

        std::vector<string> flameletList;
        // ============================
        // Extreme cases
        // phi=0
        if (!generateOnlyFlameletRange) {
          writeExtremeFlamelets(0.0);
          writeExtremeFlamelets(384.4e6);
          flameletList.push_back("flamelet_oxidizer");
          flameletList.push_back("flamelet_fuel");
        }

        flame->setRefineCriteria(flowdomain, ratio, slope, curve, prune);
        bool refine_grid = true;


        if (createInit) {
          //----------- Supply initial guess----------------------

          vector_fp locs;
          vector_fp value;

          locs.resize(3);
          value.resize(3);

          //ramp values from inlet to adiabatic flame conditions
          //  over 70% of domain and then level off at equilibrium
          double z1 = 0.7;

          flame->setTimeStep(dt, nStepsInit, steps);

          locs[0] = 0.0;
          locs[1] = z1;
          locs[2] = 1.0;
  /*
          double uout;
          uout = 0.0;//inlet.mdot()/rho_out;
          uin = 0.0;//inlet.mdot()/rho_in;
          value[0] = uin;
          value[1] = uout;
          value[2] = uout;
          flame->setInitialGuess("u",locs,value);
  */
          value[0] = Tu;
          value[1] = Tad;
          value[2] = Tad;
          flame->setInitialGuess("T", locs, value);

          for (size_t i=0; i<nsp; i++) {
            value[0] = yin[i];
            value[1] = yout[i];
            value[2] = yout[i];
            flame->setInitialGuess(mix->speciesName(i), locs, value);
          }

          inlet.setMoleFractions(&x[0]);
          inlet.setMdot(mdot);
          inlet.setTemperature(Tu);

  //        flame->showSolution();
          // Solve freely propagating flame

          refine_grid = false;
          f->fixTemperature();
          flame->setFixedTemperature(900.0);
          flame->solve(loglevel, refine_grid);
//          flame->save(initFlamelet,"fl", "init", loglevel);
        }
//        else {
//          flame->restore(initFlamelet, "fl", loglevel);
//        }

        const double lz_delta_max = 2500.0;
        const double lz_delta_min = 1900.0;

        const double PHI1  = phi1/(1.0+phi1);
        const double PHI2  = phi2/(1.0+phi2);
        const double nfl12 = 0.7*nfl;            // 70% of flamelets are in the phi1<phi<phi2 range, and 30% in the phi>phi2 range
        double this_phi = phi1;
        double dphi;
        double dphi12   = (phi2-phi1)/(nfl12-1.0);
        double dPHI12   = (PHI2-PHI1)/(nfl12-1.0);
        double dPHI2inf = (PHI_INF-PHI2) /(double(nfl)-nfl12-1.0);
        int ifl = 1;

        flame->setTimeStep(dt, nSteps, steps);
        flame->setRefineCriteria(flowdomain, ratio, slope, curve, prune);

        if (createInit) {
          refine_grid = true;
          f->solveEnergyEqn();
          flame->solve(loglevel, refine_grid);
          flame->save(initFlamelet,"fl", "init", loglevel);
        }
        else {
          flame->restore(initFlamelet, "fl", loglevel);
        }

        while ( ifl<=nfl ) {
          if (this_phi > 11.0 && p<10.0e5) break;   // XXX should find a better solution
          setInlet(this_phi);
          refine_grid = true;
          f->solveEnergyEqn();
          flame->solve(loglevel, refine_grid);

          std::vector<doublereal> Tvec = getVariableVec(f, "T");

          while ( f->grid(f->nPoints()-1) < lz_delta_min*calcFlameThickness(f, Tvec) || f->grid(f->nPoints()-1) > lz_delta_max*calcFlameThickness(f, Tvec) ) {
            int np = f->nPoints();
            if (lz < lz_delta_min*calcFlameThickness(f, Tvec))      lz = 1.05*lz_delta_min * calcFlameThickness(f, Tvec);
            else if (lz > lz_delta_max*calcFlameThickness(f, Tvec)) lz = 0.95*lz_delta_max * calcFlameThickness(f, Tvec);

            double scf = fmax( fmin(lz / f->grid(np-1), 1.5), 0.75 );
            cout << "=========== Domain length changed by a factor of " << scf <<  " ===========" << endl;
            double *y = new double [np];
            for (int ip=0; ip<np; ip++)
              y[ip] = f->grid(ip) * scf;
            f->setupGrid(np, &y[0]);
            f->solveEnergyEqn();
            flame->solve(loglevel, refine_grid);
            Tvec = getVariableVec(f, "T");
            delete[] y;
          }

          if (transport == "Multi") {
//            setLewis(-1);
            f->setTransport(*trmulti, Le);
            if (ifl == 1) {
              cout << "========= Activating Multicomponent transport ==========" << endl;
              if (SORET) {
                cout << "========= Activating Soret effect ==========" << endl;
                f->enableSoret(true);
              }
              flame->solve(loglevel, refine_grid);
              Tvec = getVariableVec(f, "T");
            }
          }

          if (writeCSV) writeFlameletCSV(f, flameletName(this_phi)+".dat");
          flame->save(flameletName(this_phi)+".xml", "fl", "flm", loglevel);
          writeFlamelet(f, flameletName(this_phi), this_phi, Tvec);
          flameletList.push_back(flameletName(this_phi));
          cout << " ((((( " << ifl++ << ": solve premixed flamelet: phi = " << this_phi << "; Tu = " << Tu << "; Tmax = " << *(std::max_element(Tvec.begin(), Tvec.end())) << "; np = " << f->nPoints() << " ))))) " << endl;

          if (this_phi < phi2)
            dphi = dphi12;//(1.0 + this_phi*this_phi) * dPHI12;
          else
            dphi = (1.0 + this_phi)*(1.0 + this_phi) * dPHI2inf;

          this_phi += dphi;

          if (generateOnlyFlameletRange && this_phi>phi2)
            break;
        }

        writeFlameletList(flameletList);
        clock_t te = clock();
        cout << ".... Flamelet generation time: " << double(te-tb)/CLOCKS_PER_SEC << endl;

      }
      catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << "program terminating." << endl;
        return -1;
      }

    return 0;
    }

    void setInlet(double phi) {
      doublereal N_O2_st = calcO2Stoich();
      if (phi>=100.0) N_O2_st = 0.0;
      vector_fp x(nsp);

      for (size_t k=0; k<nsp; k++) {
        x[k]   = Xf[k]*phi + Xox[k]/Xox[i_O2]*N_O2_st;
      }

      mix_i->setState_TPX(Tf, p, &Xf[0]);
      double h_mix = mix_i->enthalpy_mass() * phi*mix_i->meanMolecularWeight();
      double m_tot = phi*mix_i->meanMolecularWeight();
      mix_i->setState_TPX(Tox, p, &Xox[0]);
      h_mix += mix_i->enthalpy_mass() *N_O2_st/Xox[i_O2]*mix_i->meanMolecularWeight();
      m_tot += N_O2_st/Xox[i_O2]*mix_i->meanMolecularWeight();

      mix->setMoleFractions(&x[0]);  // set the mixture mole fractions
      mix->setState_HP(h_mix/m_tot, p);    // set mixture enthalpy [J/Kg] and pressure

      Tu = mix->temperature();

      inlet.setMoleFractions(&x[0]);
      inlet.setTemperature(Tu);
    }

    string flameletName(double phi) {
      string prefix = "";
      int is = 0;
      while (fuel_compo[is]!=':')
        prefix += fuel_compo[is++];
//      for (int is=0; is<fuel_compo.size(); is++)
//        if (fuel_compo[is]!=':' && fuel_compo[is]!=' ') prefix += fuel_compo[is];

      std::stringstream strm;
            strm << std::fixed << "_p" << std::setprecision(2) << p/1.0e5 << "_Tu" << Tu << "_phi" << std::setprecision(5) << phi;

      string name = prefix + strm.str();
      return name;
    }

    void writeExtremeFlamelets(const double phi) {
      flame->setFlatProfile(flowdomain, f->componentIndex("u"), 0.0);

      if (phi==0.0) {
        flame->setFlatProfile(flowdomain, f->componentIndex("T"), Tox);
        mix_i->setState_TPX(Tox, p, &Xox[0]);
      }
      else {
        flame->setFlatProfile(flowdomain, f->componentIndex("T"), Tf);
        mix_i->setState_TPX(Tf, p, &Xf[0]);
      }

      for (int isp=0; isp<nsp; isp++)
        if (phi==0.0)
          flame->setFlatProfile(flowdomain, f->componentIndex(mix_i->speciesName(isp)), Xox[isp]*mix_i->molecularWeight(isp)/mix_i->meanMolecularWeight());
        else
          flame->setFlatProfile(flowdomain, f->componentIndex(mix_i->speciesName(isp)), Xf[isp]*mix_i->molecularWeight(isp)/mix_i->meanMolecularWeight());

      std::vector<doublereal> Tvec = getVariableVec(f, "T");
      if (phi==0.0)
        writeFlamelet(f, "flamelet_oxidizer", 0.0, Tvec);
      else {
        writeFlamelet(f, "flamelet_fuel", 384.4e6, Tvec);
      }
//      flameletList.push_back(flameletName(0.0));

    }

    template <class flowType>
    void writeFlamelet(flowType *flow, string fileName, double phi, std::vector<doublereal> Tvec) {
//      void *flow = f;
//      if (f != Null) flow = (FreeFlame*)f;
//      if (f != Null && fstgn == Null)      flow = f;
//      else if (fstgn != Null && f == Null) flow = fstgn;
      int np = flow->nPoints();

      FILE* FP = fopen(fileName.c_str(), "w");
      if (!FP) {
          printf("Failure to open file\n");
          exit(-1);
      }

      std::vector<doublereal> Uvec = getVariableVec(flow, "u");
//      std::vector<doublereal> Tvec = getVariableVec(f, "T");
      doublereal *rhovec = new doublereal[np];
      doublereal *mw     = new doublereal[np];
      doublereal *cp     = new doublereal[np];
      doublereal *mu     = new doublereal[np];
      doublereal *lambda = new doublereal[np];
      doublereal *h      = new doublereal[np];
      calcThermoProp(flow, rhovec, mw, cp, mu, lambda, h, Tvec);

      fprintf(FP, "header\n\n");
      fprintf(FP, "title = \"free premixed flame\"\n");
      fprintf(FP, "mechanism = \"%s\"\n", mechanism.c_str());
      fprintf(FP, "author = \"Cascade Technologies\"\n");

      fprintf(FP, "\nfuel = \"%s\"\n", fuel_compo.c_str());
      fprintf(FP, "oxidizer = \"%s\"\n", oxidizer_compo.c_str());
      fprintf(FP, "pressure = %7.3f [bar]\n", p/1.0e5);
      fprintf(FP, "fuel-air-equivalence-ratio = %8.6f\n", phi);
      fprintf(FP, "mixture-fraction-stoichiometric = %8.6f\n", calcZStoich());
      fprintf(FP, "Tmax = %7.2f [K]\n", *(std::max_element(Tvec.begin(), Tvec.end())));

      fprintf(FP, "\nunburnt\n");
      fprintf(FP, "begin\n");
      fprintf(FP, "\tTemperature = %7.3f [K]\n", Tu);
      vector_fp Yu = getMassFracVec(flow, 0);
      for (size_t k=0; k<nsp; k++) {
        if (Yu[k]>1.0e-8) {
//          fprintf(FP, "\tMolefraction-%s = %7.6f\n", mix->speciesName(k).c_str(), Xu[k]);
          fprintf(FP, "\tMassfraction-%s = %7.6f\n", mix->speciesName(k).c_str(), Yu[k]);
        }
      }
      fprintf(FP, "end\n");

      fprintf(FP, "\nburningVelocity = %8.4f [cm/sec]\n", Uvec[0]*100.0);
      fprintf(FP, "FlameThickness = %11.5e [m]\n", calcFlameThickness(flow, Tvec));
      fprintf(FP, "numOfSpecies = %d\n", int(nsp));
      fprintf(FP, "gridPoints = %d\n", int(flow->nPoints()));

      fprintf(FP,"\nbody\n");

      //////////////////////////////////////////////////////////////////
      //////////////////// write grid, rho, mdot, and T /////////////////////
      fprintf(FP,"y [m]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", flow->grid(n));
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"density [kg/m^3]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", rhovec[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"massflowrate [kg/m^2s]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", rhovec[n]*Uvec[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"temperature [K]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", Tvec[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"W [kg/kmole]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", mw[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"cp [J/kgK]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", cp[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"mu [kg/ms]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", mu[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"lambda [J/mKs]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", lambda[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"TotalEnthalpy [J/kg]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", h[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }

      /////////////////////////////////////////////////////////////////
      ///////////////// write species mass fraction ///////////////////

      for (size_t k=0; k<nsp; k++) {
        const char *spName = mix->speciesName(k).c_str();
        std::vector<doublereal> Xsp = getVariableVec(flow, spName);
        string headerName = "massfraction-"+ mix->speciesName(k);
        fprintf(FP,"%s\n", headerName.c_str());

        for (int n = 0; n < np; n++) {
//          double this_y = Xsp[n] * mix->molecularWeight(spName)/
          const double tmp = ((abs(Xsp[n])<1e-99) ? 0.0 : Xsp[n]);
          fprintf(FP,"\t%11.6e",tmp);
//          fprintf( FP,"\t%11.6e", flame->value( flowdomain,f->componentIndex(spName),n ) );
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
        }
      }

      /////////////////////////////////////////////////////////////////
      //////////////// write species production rate //////////////////
      double **ProdRate;
      getMem2D(&ProdRate, 0, np-1, 0, nsp-1, "ProdRate", true);

      for (int n = 0; n < np; ++n) {
        vector_fp this_ProdRate = calcProdRate_i(Tvec[n], getMassFracVec(flow, n));
        for (size_t k = 0; k < nsp; ++k)
          ProdRate[n][k] = this_ProdRate[k];
      }

      for (size_t k = 0; k < nsp; ++k) {
        const char *spName = mix->speciesName(k).c_str();
        string headerName = "ProdRate-"+ mix->speciesName(k)+" [kg/m^3s]";
        fprintf(FP,"%s\n", headerName.c_str());

        for (int n = 0; n < np; ++n) {
          double this_pr = ProdRate[n][k] * mix->molecularWeight(k);
          const double tmp = ((abs(this_pr)<1e-99) ? 0.0 : this_pr);
          fprintf(FP,"\t%11.6e",tmp);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
        }
      }

      /////////////////////////////////////////////////////////////////
      /////////////////// write heat release rate /////////////////////
      fprintf(FP,"HeatRelease [J/m^3s]\n");
      for (int n = 0; n < np; ++n) {
        fprintf( FP,"\t%11.6e", calcHeatRelease(Tvec[n], getMassFracVec(flow, n)) );
        if ( n%5 == 4 || n == np-1 ) fprintf(FP,"\n");
      }

      ///////////////////////////////////////////////////////
      /////////////////// write entropy /////////////////////
      fprintf(FP,"entropy [J/kgK]\n");
      for (int n = 0; n < np; ++n) {
        fprintf( FP,"\t%11.6e", calcEntropy(Tvec[n], getMassFracVec(flow, n)) );
        if ( n%5 == 4 || n == np-1 ) fprintf(FP,"\n");
      }

      fclose(FP);
      freeMem2D(ProdRate, 0, np-1, 0, nsp-1);
      delete[] rhovec;
      delete[] mw;
      delete[] mu;
      delete[] cp;
      delete[] lambda;
      delete[] h;
    }


};

class Diffusion: public Flame {
  private:
    CanteraLe::AxiStagnFlowLe *f;
//    StFlow *f;
    Inlet1D inletL;
    Inlet1D inletR;
    doublereal lz;
    int nz;
    static constexpr doublereal Lmax = 10.0e5;
    static constexpr doublereal Lmin = 2.0e-6;
    doublereal rho_ox;
    doublereal rho_f;
    bool createInit;
    string initFlamelet;
    bool writeCSV;
    double Terror;

    doublereal aInit;

    std::vector<double> TmaxVec;
    std::vector<double> aVec;

  public:
    Diffusion(const string fuel, string ox, const double pu, const double Tfuel, const double Toxidizer, const string mech, const string tr,
              const double Lewis = 1.0, const int lglvl = 0,
              const bool init = true, const int ng = 31, const string initFile = "init.xml", doublereal a1 = 0.5,
              const bool csv = false, const bool soret = false,
//              const double rat = 10.0, const double slo = 0.1, const double cur = 0.2, const double pru = 0.075) :
              const double rat = 15.0, const double slo = 0.25, const double cur = 0.5, const double pru = 0.15) :
              Flame(fuel, ox, pu, Tfuel, Toxidizer, mech, tr, Lewis, lglvl, soret, rat, slo, cur, pru),
              nz(ng), createInit(init), initFlamelet(initFile), Terror(1.0e-1), writeCSV(csv), aInit(a1) {}

    int createFlamelets(const int nfl = 1) {

      try {
        doublereal a          = aInit;//1.0;//0.05;//1.0;//0.5;  // initial strain rate
        const double phi_init = 1.0;

        lz = 0.25;
//        int nz = 31;
        doublereal dt = 1.0e-8;  // time step

        double ds          = 0.1125;//0.25;//0.35;//0.05;
        const double dsMax = 0.5;//0.5;//0.35;//0.25;
        const double dsMin = 0.01;//0.01;//0.025;
        const double a_ratio_max   = 5.0;
        double a_T_min             = 1.0e-6;  // should be a small number
        const double a_T_min_scale = 2.0;
        const double aFactor2nd    = 10.0;
        const double verySmallA    = 1.0e-2;

        const double rtol    = 1.0e-8;
        const double atol    = 1.0e-12;
        const double rtol_tr = 1.0e-5;
        const double atol_tr = 1.0e-8;
        const double dx_min  = 1.0e-6;
        const double T_2ndError = 25.0;//50.0;  // Kelvin
        const double dT_max     = 10.0;  // Kelvin
        const double Textinction = 350.0;

        vector_fp yox(nsp);
        mix->setState_TPX(Tox, p, &Xox[0]);
        mix->getMassFractions(&yox[0]);
        rho_ox = mix->density();

        vector_fp yf(nsp);
        mix->setState_TPX(Tf, p, &Xf[0]);
        mix->getMassFractions(&yf[0]);
        rho_f = mix->density();

        setInlet(a);

        setPremixedMix(phi_init);

        try {
          mix->equilibrate("HP");
        }
        catch (CanteraError& err) {
          std::cout << err.what() << std::endl;
        }

        vector_fp yeq(nsp);
        mix->getMassFractions(&yeq[0]);
        doublereal rho_eq = mix->density();
        doublereal Tad = mix->temperature();
//        doublereal ueq = mdot/mix->density();
        cout << phi_init << ' ' << Tad << endl;

        clock_t tb = clock();

        //=============  build each domain ========================
        size_t nStepsInit = 2;
        size_t nSteps     = 4;       // should find the optimized combination
        integer steps[]   = {2, 5, 10, 20, 40, 80, 120};
        flowdomain = 1;

        //-------- step 1: create the flow -------------
        f = new CanteraLe::AxiStagnFlowLe(mix);

        if (createInit) {
          // create an initial grid
          vector_fp z(nz);
          doublereal dz = lz/((doublereal)(nz-1));
          for (int iz=0; iz<nz; iz++) {
            z[iz] = ((doublereal)iz)*dz;
          }
          f->setupGrid(nz, &z[0]);
        }

        // specify the objects to use to compute kinetic rates and
        // transport properties
        std::auto_ptr<Transport> trmix(newTransportMgr("Mix", mix));
        std::auto_ptr<Transport> trmulti(newTransportMgr("Multi", mix));

        f->setTransport(*trmix, Le);
        f->setKinetics(*mix);
        f->setPressure(p);

        //=================== create the container and insert the domains ====================//

        std::vector<Domain1D*> domains;
        domains.push_back(&inletL);
        domains.push_back(f);
        domains.push_back(&inletR);

        flame = new Sim1D(domains);
        flame->setTimeStep(dt, nStepsInit, steps);
        //        flame->showSolution();
        flame->setRefineCriteria(flowdomain, ratio, slope, curve, prune);
        flame->setGridMin(flowdomain, dx_min);
        f->setTransientTolerances(rtol_tr, atol_tr);
        f->setSteadyTolerances(rtol, atol);
        flame->setTimeStep(dt, nSteps, steps);


        std::vector<string> flameletList;

        //----------- Supply initial guess----------------------

        if (createInit) {
          vector_fp locs;
          vector_fp value;
          locs.resize(6);
          value.resize(6);

          locs[0] = 0.0;
          locs[1] = 0.2;
          locs[2] = 0.3;
          locs[3] = 0.7;
          locs[4] = 0.8;
          locs[5] = 1.0;

          value[0] = Tox;
          value[1] = Tox;
          value[2] = Tad;
          value[3] = Tad;
          value[4] = Tf;
          value[5] = Tf;
          flame->setProfile(flowdomain, f->componentIndex("T"), locs, value);

          for (size_t i=0; i<nsp; i++) {
            value[0] = yox[i];
            value[1] = yox[i];
            value[2] = yeq[i];
            value[3] = yeq[i];
            value[4] = yf[i];
            value[5] = yf[i];
            flame->setProfile(flowdomain, f->componentIndex(mix->speciesName(i)), locs, value);
          }
        }

        setInlet(a);

        bool refine_grid = true;

        int ifl = 1;
        int iflExt;
        double TmaxExt;
        double a1, a2, T1, T2;
        double Tmax = 2.0*fmax(Tf, Tox);
        double a_old;
        double dsa_dsT, dsa_dsT_old;

        if (createInit) {
          cout << " ======================================================================== " << endl;
          cout << " >>>>>>>>>>>>>>>>>>> generate initial flamelet <<<<<<<<<<<<<<<<<<<<<<<<<< " << endl;
          cout << " ======================================================================== " << endl;
          setInlet(a);

          f->solveEnergyEqn();
          setInlet(a);
          flame->solve(loglevel, false);
          writeFlameletCSV(f, "CH4_noRefine_couple.dat");

          setInlet(a);
          flame->solve(loglevel, true);
          writeFlameletCSV(f, "CH4_refine_couple.dat");

          flame->save(initFlamelet,"fl", "init", loglevel);
        }
        else {
          flame->restore(initFlamelet, "fl", loglevel);
        }

        refine_grid = true;

//        bool isExtinct               = false;
        bool inUnstableBranch        = false;
        bool sCurveDone              = false;
        bool reached1stTurning       = false;
        bool firstFlameletInUnstable = false;

        scaleGrid(1.0);

        cout << " =================================================================== " << endl;
        cout << " >>>>>>>>>>>>>>>>>>> starting the s-curve <<<<<<<<<<<<<<<<<<<<<<<<<< " << endl;
        cout << " =================================================================== " << endl;
        flame->setRefineCriteria(flowdomain, ratio, slope, curve, prune);
        cout << "ratio = " << ratio << "; slope = " <<  slope << "; curve = " <<  curve << "; prune = " << prune << endl;

        while ( ifl<=nfl && !sCurveDone) {

          setInlet(a);

          f->solveEnergyEqn();
//          cout << " ((((( " << ifl << ": solve diffusion flamelet: a = " << a << "; Tmax = " << calcMax("T") << "; np = " << f->nPoints() << " ))))) " << endl;
          if (!firstFlameletInUnstable) flame->solve(loglevel, refine_grid);

          std::vector<doublereal> Tvec = getVariableVec(f, "T");

          writeFlamelet(f, flameletName(ifl, a), a, Tvec);
          if (writeCSV) writeFlameletCSV(f, flameletName(ifl, a)+".dat");
          flameletList.push_back(flameletName(ifl, a));
          Tmax = calcMax("T");
//          isExtinct = (fabs(Tmax - fmax(Tf, Tox)) < Terror);

          if (feelTheBurn()) {
            a_old       = a;
            dsa_dsT_old = dsa_dsT;

            flame->save(flameletName(ifl, a)+".xml", "fl", "flm", loglevel);
            TmaxVec.push_back(Tmax);
            aVec.push_back(a);
            cout << " ((((( " << ifl++ << ": solve diffusion flamelet: a = " << a << "; Tmax = " << Tmax << "; np = " << f->nPoints() << " ))))) " << endl;
            firstFlameletInUnstable = false;
          }
          else {
            cout << " ((((( " << " --- extinct flamelet: a = " << a << "; Tmax = " << Tmax << "; np = " << f->nPoints() << " ))))) " << endl;
          }

          if ( (Tmax - fmax(Tf, Tox) <= T_2ndError || Tmax < Textinction) && Tmax - fmax(Tf, Tox) > Terror) {
            cout << "S-curve is complete. Going to make 5 more flamelets." << endl;
            rescaleT(0.0);
            for (int isp=0; isp<nsp; isp++) {
              rescaleVar(mix_i->speciesName(isp), 0.0);
            }

            double this_a = a;
            for (int i=0; i<5; i++) {
              this_a *= aFactor2nd / (5.0-double(i));
              setInlet(this_a);

//              double delta_ratio = 1.0 / sqrt(aFactor2nd);
//              scaleGrid(delta_ratio);
//              scaleSelfsimilarSolution(aFactor2nd);

              setReactionsMultiplier(0.0);

              refine_grid = false;
              f->solveEnergyEqn();
              flame->solve(loglevel, refine_grid);
              Tmax = calcMax("T");
              sCurveDone = (fabs(Tmax - fmax(Tf, Tox)) < Terror);
              if (sCurveDone) {
                cout << " ((((( " << ifl++ << ": solve diffusion flamelet: a = " << a << "; Tmax = " << Tmax << "; np = " << f->nPoints() << " ))))) " << endl;
                writeFlamelet(f, flameletName(ifl-1, this_a), this_a, Tvec);
//                flameletList.push_back(flameletName(ifl-1, this_a));
                TmaxVec.push_back(Tmax);
                aVec.push_back(this_a);

                setReactionsMultiplier(1.0);
              }
            }
            refine_grid = true;
          }
          else if (ifl>2) {
            firstFlameletInUnstable = false;

            if (!feelTheBurn()) {
              flame->save("mixing.xml", "mixing", "mixing", loglevel);
              flame->restore(flameletName(ifl-1, a_old)+".xml", "fl", loglevel);
              a = a_old;
              a_T_min = 1.5*dsa_dsT_old;

              std::vector<doublereal> Tvec = getVariableVec(f, "T");
              Tmax = calcMax("T");
//              isExtinct = (fabs(Tmax - fmax(Tf, Tox)) < Terror);
              reached1stTurning = true;
              iflExt  = ifl;
              TmaxExt = Tmax;
            }

            a1 = aVec[ifl-3];
            a2 = aVec[ifl-2];
            T1 = TmaxVec[ifl-3];
            T2 = TmaxVec[ifl-2];

            double ds_T = T2-T1;
            double ds_a = 1/a2 - 1/a1;
            double c = sqrt( ((1-T1/T2)*(1-T1/T2)) / ((1-a2/a1)*(1-a2/a1)) + 1.0 );
            dsa_dsT = ds_a/ds_T;

            double this_ds = ds / ( fabs(ds_T/dT_max));//*fabs(ds_T/dT_max));//*fabs(ds_T/dT_max) );
            if (dsa_dsT >= 0) {
              ds = fmax( fmin(this_ds, dsMax), dsMin );
            }
            else {
//              if (a>verySmallA)
              if (TmaxExt/Tmax-1<0.1)//(ifl-iflExt>10)
                ds = fmax(fmin(this_ds, dsMax/2.0), dsMin*5);//dsMax;
              else
                ds = fmax(fmin(this_ds, dsMax/4.0), dsMin*5);
//              else
//                ds = fmin(this_ds, dsMax/10.0);//dsMax;
            }

//            assert(ds_T<0);
//            cout << "ds = " << ds << "; ds_T = " << ds_T << "; ds_a = " << ds_a << "; ds_a/ds_T = " << dsa_dsT << endl << endl;

            if ( reached1stTurning ) {
              if ( !inUnstableBranch ) {
                flame->save("extinct.xml","id2", "ext", loglevel);
                double this_Tmax = Tmax;
                double sc = 0.9;//       = 0.8;//0.8;
                double sc_min   = 0.05;
                double a_ratio  = 0.9;
                double sc_ratio = 0.9;//0.9;
//                a  *= a_ratio;
                a = a1;
//                a_ratio  = a_ratio;//fmin(a1/a2, a_ratio);
//                sc_ratio = sc_ratio;//fmin(T2/T1, sc_ratio);
                const double dtemp = T1 - T2;
                const double T3    = T2 - dtemp;
                sc = fmin(T3/T2, sc);
                cout << "+++++++++ a1/a2 = " << a_ratio << "; T1/T2 = " << sc_ratio << "; sc = " << sc << endl;
                double delta_ratio = 1.0 / sqrt(a_ratio);
                scaleGrid(delta_ratio);
                scaleSelfsimilarSolution(a_ratio);

                int k = 0;
                while (this_Tmax - Tmax >= 0 || !feelTheBurn() ) {
                  flame->restore("extinct.xml", "id2", loglevel);
//                  flame->restore("extinct_rescaled.xml", "id2", loglevel);

                  if (this_Tmax >= Tmax) rescaleT(sc);
                  setInlet(a);

                  f->fixTemperature();
                  flame->solve(loglevel, refine_grid);
                  f->solveEnergyEqn();
                  flame->solve(loglevel, refine_grid);

//                  if (k++ > 6) a_ratio = 1.0;

                  this_Tmax = calcMax("T");
                  cout << ">>>>>>>>>> Reached 1st turning point. This Tmax = " << this_Tmax << " should be smaller than " << Tmax
                       << "; a = " << a << "; next scaling factor = " << sc << "; Np = " << f->nPoints() << " <<<<<<<<<<" << endl;

                  if ( this_Tmax < Tmax ) {
                    firstFlameletInUnstable = true;
                  }

                  if (this_Tmax >= Tmax || !feelTheBurn()) {
                    a *= a_ratio;
                    double this_a_ratio = a/a1;
                    double delta_ratio = 1.0 / sqrt(this_a_ratio);
                    scaleGrid(delta_ratio);
                    scaleSelfsimilarSolution(this_a_ratio);
                    sc = fmax(sc*sc_ratio, sc_min);
                  }

                }

              }
              else a = fmin( 1.0/(ds/c+1.0), a_ratio_max ) * a2;

              c *= -1;
            }
            else
              a = fmin( 1.0/(-ds/c+1.0), a_ratio_max ) * a2;

            if ( reached1stTurning ) inUnstableBranch = true;

            if (!firstFlameletInUnstable) {
              double a_ratio = a/a2;
              double delta_ratio = 1.0 / sqrt(a_ratio);
              scaleGrid(delta_ratio);
              scaleSelfsimilarSolution(a_ratio);
            }
          }
          else {
            a2 = a;
            T2 = *(std::max_element(Tvec.begin(), Tvec.end()));

            double a_ratio = 1.1;
            a  *= a_ratio;
            double delta_ratio = 1.0 / sqrt(a_ratio);
            scaleGrid(delta_ratio);
            scaleSelfsimilarSolution(a_ratio);
          }

        }

        writeFlameletList(flameletList);
        writeScurve();
        clock_t te = clock();
        cout << ".... Flamelet generation time: " << double(te-tb)/CLOCKS_PER_SEC << endl;

      }
      catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << "program terminating." << endl;
        return -1;
      }

    return 0;
    }

    void lowerBranch(double a, const string firstFlamelet) {

      cout << "lowerBranch: a = " << a << endl;
      flame->restore(firstFlamelet, "mixing", loglevel);
      bool refine_grid = true;
      int ifl = 1;
      while ( ifl<=1000 ) {

        double a_ratio = 0.9;
        a *= a_ratio;
        double delta_ratio = 1.0 / sqrt(a_ratio);
        scaleGrid(delta_ratio);
        scaleSelfsimilarSolution(a_ratio);

        setInlet(a);
        f->solveEnergyEqn();
        flame->solve(loglevel, refine_grid);

        std::vector<doublereal> Tvec = getVariableVec(f, "T");

        writeFlamelet(f, flameletName(ifl, a), a, Tvec);
        if (writeCSV) writeFlameletCSV(f, flameletName(ifl, a)+".dat");
//        flameletList.push_back(flameletName(ifl, a));
        double Tmax = calcMax("T");
//        isExtinct = (fabs(Tmax - fmax(Tf, Tox)) < Terror);

//        if (!isExtinct) {
//        a_old       = a;
//        dsa_dsT_old = dsa_dsT;

        flame->save(flameletName(ifl, a)+".xml", "fl", "flm", loglevel);
//        TmaxVec.push_back(Tmax);
//        aVec.push_back(a);
        cout << " ((((( " << ifl++ << ": solve diffusion flamelet: a = " << a << "; Tmax = " << Tmax << "; np = " << f->nPoints() << " ))))) " << endl;
//        }

      }

    }

    void setExoReactionsMultiplier(const double sc) {
      int Nr = mix->nReactions();
      vector_fp dh(Nr);
      mix->getDeltaEnthalpy(&dh[0]);

      for (int ir=0; ir<Nr; ir++) {
//        cout << "reaction " << ir << " isRiversible = " << mix->isReversible(ir) << "; delta_enthalpy = " << dh[ir] << endl;
        if (dh[ir] < 0)
          mix->setMultiplier(ir, sc);
      }
    }

    void setReactionsMultiplier(const double sc) {
      int Nr = mix->nReactions();
      for (int ir=0; ir<Nr; ir++) {
          mix->setMultiplier(ir, sc);
      }
    }

    double calcMax(const string var) {
      std::vector<doublereal> varVec = getVariableVec(f, var);
      return *(std::max_element(varVec.begin(), varVec.end()));
    }

    void scaleGrid(const double scf) {
      int np = f->nPoints();
      double *y = new double [np];

      for (int ip=0; ip<np; ip++) {
        double this_l = (f->grid(np-1) - f->grid(0)) * scf;
        double this_scf = fmax(fmin(this_l, Lmax), Lmin ) / (f->grid(np-1) - f->grid(0));
        y[ip] = f->grid(ip) * this_scf;
      }

      f->setupGrid(np, &y[0]);
      lz = f->grid(np-1) - f->grid(0);
      cout << "************** changed domain size. new L = " << lz << endl;
      delete[] y;
    }

    bool feelTheBurn() {
      return ( fabs( calcMax("T") - fmax(Tf,Tox) ) > Terror );
    }

    void scaleSelfsimilarSolution(const double a_ratio, const double p_ratio = 1.0) {
      int np = f->nPoints();
      vector_fp locs(np);
      vector_fp values(np);

      // scale axial velocity
      for (int i=0; i<np; i++) {
        locs[i]   = f->grid(i) / lz;
        values[i] = flame->value( flowdomain, f->componentIndex("u"), i) * sqrt(a_ratio/p_ratio);
      }
      flame->setProfile(flowdomain, f->componentIndex("u"), locs, values);

      // scale radial valocity; note: V = v/r
      for (int i=0; i<np; i++) {
        values[i] = flame->value( flowdomain, f->componentIndex("V"), i) * a_ratio;
      }
      flame->setProfile(flowdomain, f->componentIndex("V"), locs, values);

      // scale pressure curvature
      for (int i=0; i<np; i++) {
        values[i] = flame->value( flowdomain, f->componentIndex("lambda"), i) * a_ratio*a_ratio * p_ratio;
      }
      flame->setProfile(flowdomain, f->componentIndex("lambda"), locs, values);
    }

    void rescaleT(const double scf) {
      int np = f->nPoints();
      vector_fp locs(np);
      vector_fp values(np);

      // scale axial velocity
      for (int i=0; i<np; i++) {
        locs[i]   = f->grid(i) / lz;
        const double T  = flame->value( flowdomain, f->componentIndex("T"), i);
        const double Ts = Tox + (Tf-Tox)*locs[i];
        values[i] = (T - Ts)*scf + Ts;
      }
      flame->setProfile(flowdomain, f->componentIndex("T"), locs, values);
    }

    void rescaleVar(const string varName, const double scf) {
      int np = f->nPoints();
      vector_fp locs(np);
      vector_fp values(np);
      std::vector<doublereal> varOld = getVariableVec(f, varName);

      // scale axial velocity
      for (int i=0; i<np; i++) {
        locs[i]   = f->grid(i) / lz;
        const double var  = flame->value( flowdomain, f->componentIndex(varName), i);
        const double varS = varOld[0] + (varOld[np-1]-varOld[0])*locs[i];
        values[i] = (var - varS)*scf + varS;
      }
      flame->setProfile(flowdomain, f->componentIndex(varName), locs, values);
    }

    void setPremixedMix(double phi) {
      doublereal N_O2_st = calcO2Stoich();
      if (phi>=100.0) N_O2_st = 0.0;
      vector_fp x(nsp);

      for (size_t k=0; k<nsp; k++) {
        x[k]   = Xf[k]*phi + Xox[k]/Xox[i_O2]*N_O2_st;
      }

      mix->setState_TPX(Tf, p, &Xf[0]);
      double h_mix = mix->enthalpy_mass() * phi*mix->meanMolecularWeight();
      double m_tot = phi*mix->meanMolecularWeight();
      mix->setState_TPX(Tox, p, &Xox[0]);
      h_mix += mix->enthalpy_mass() *N_O2_st/Xox[i_O2]*mix->meanMolecularWeight();
      m_tot += N_O2_st/Xox[i_O2]*mix->meanMolecularWeight();

      mix->setMoleFractions(&x[0]);  // set the mixture mole fractions
      mix->setState_HP(h_mix/m_tot, p);    // set mixture enthalpy [J/Kg] and pressure
    }

    void setInlet(const double a) {

      doublereal uL = a * lz;
      double mdot_ox = rho_ox*uL;
      double mdot_f  = mdot_ox * sqrt(rho_f/rho_ox);

      inletL.setMoleFractions(&Xox[0]);
      inletL.setTemperature(Tox);
      inletL.setMdot(mdot_ox);

      inletR.setMoleFractions(&Xf[0]);
      inletR.setTemperature(Tf);
      inletR.setMdot(mdot_f);

//      inletL.showSolution(0);
//      inletR.showSolution(0);
    }

    string flameletName(const int ifl, const double a) {
      string prefix = "";
      int is = 0;
      while (fuel_compo[is]!=':')
        prefix += fuel_compo[is++];

      std::stringstream strm;
      strm << std::fixed << "_p" << std::setprecision(2) << p/1.0e5 << "_Tf" << Tf << "_Tox" << Tox << "_fl" << std::setw(3) << std::setfill('0') << ifl << std::setprecision(5) << "_a" << a;

      string name = prefix + strm.str();
      return name;
    }

    template <class flowType>
    void writeFlamelet(flowType *flow, string fileName, const double a, std::vector<doublereal> Tvec) {

      int np = flow->nPoints();

      FILE* FP = fopen(fileName.c_str(), "w");
      if (!FP) {
          printf("Failure to open file\n");
          exit(-1);
      }

      std::vector<doublereal> Uvec = getVariableVec(flow, "u");
      doublereal *rhovec = new doublereal[np];
      doublereal *mw     = new doublereal[np];
      doublereal *cp     = new doublereal[np];
      doublereal *mu     = new doublereal[np];
      doublereal *lambda = new doublereal[np];
      doublereal *h      = new doublereal[np];
      calcThermoProp(flow, rhovec, mw, cp, mu, lambda, h, Tvec);

      std::vector<doublereal> Zvec = calcMixtureFraction(f, "H");
      for (int n = 0; n < np; n++) {
        if ( std::abs(Zvec[n]-1.0) < EPSC) {
          np = n;
          break;
        }
      }

      fprintf(FP, "header\n\n");
      fprintf(FP, "title = \"counterflow diffusion flame\"\n");
      fprintf(FP, "mechanism = \"%s\"\n", mechanism.c_str());
      fprintf(FP, "author = \"Cascade Technologies\"\n");

      fprintf(FP, "\nfuel = \"%s\"\n", fuel_compo.c_str());
      fprintf(FP, "pressure = %7.3f [bar]\n", p/1.0e5);
      fprintf(FP, "strain-rate = %8.6f [1/s]\n", a);
      fprintf(FP, "Tmax = %7.2f [K]\n", *(std::max_element(Tvec.begin(), Tvec.end())));

      fprintf(FP, "\nFuelSide\n");
      fprintf(FP, "begin\n");
      fprintf(FP, "\tTemperature = %7.3f [K]\n", Tf);
      vector_fp Yf = getMassFracVec(flow, np-1);
      for (size_t k=0; k<nsp; k++) {
        if (Yf[k]>1.0e-8) {
          fprintf(FP, "\tMassfraction-%s = %7.6f\n", mix->speciesName(k).c_str(), Yf[k]);
        }
      }
      fprintf(FP, "end\n");

      fprintf(FP, "\nOxidizerSide\n");
      fprintf(FP, "begin\n");
      fprintf(FP, "\tTemperature = %7.3f [K]\n", Tox);
      vector_fp Yox = getMassFracVec(flow, 0);
      for (size_t k=0; k<nsp; k++) {
        if (Yox[k]>1.0e-8) {
          fprintf(FP, "\tMassfraction-%s = %7.6f\n", mix->speciesName(k).c_str(), Yox[k]);
        }
      }
      fprintf(FP, "end\n\n");

//      fprintf(FP, "FlameThickness = %11.5e [m]\n", calcFlameThickness(flow, Tvec));
      fprintf(FP, "numOfSpecies = %d\n", int(nsp));
      fprintf(FP, "gridPoints = %d\n", np);

      fprintf(FP,"\nbody\n");

      //////////////////////////////////////////////////////////////////////////
      //////////////////// write Z, grid, rho, mdot, and T /////////////////////
      fprintf(FP,"Z\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", Zvec[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"y [m]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", flow->grid(n));
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"density [kg/m^3]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", rhovec[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"massflowrate [kg/m^2s]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", rhovec[n]*Uvec[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"temperature [K]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", Tvec[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"W [kg/kmole]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", mw[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"cp [J/kgK]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", cp[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"mu [kg/ms]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", mu[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"lambda [J/mKs]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", lambda[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"lambdaOverCp [kg/ms]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", lambda[n]/cp[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }
      fprintf(FP,"TotalEnthalpy [J/kg]\n");
      for (int n = 0; n < np; n++) {
          fprintf(FP,"\t%11.6e", h[n]);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
      }

      /////////////////////////////////////////////////////////////////
      ///////////////// write species mass fraction ///////////////////

      for (size_t k=0; k<nsp; k++) {
        const char *spName = mix->speciesName(k).c_str();
        std::vector<doublereal> Xsp = getVariableVec(flow, spName);
        string headerName = "massfraction-"+ mix->speciesName(k);
        fprintf(FP,"%s\n", headerName.c_str());

        for (int n = 0; n < np; n++) {
//          double this_y = Xsp[n] * mix->molecularWeight(spName)/
          const double tmp = ((abs(Xsp[n])<1e-99) ? 0.0 : Xsp[n]);
          fprintf(FP,"\t%11.6e",tmp);
//          fprintf( FP,"\t%11.6e", flame->value( flowdomain,f->componentIndex(spName),n ) );
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
        }
      }

      /////////////////////////////////////////////////////////////////
      //////////////// write species production rate //////////////////
      double **ProdRate;
      getMem2D(&ProdRate, 0, np-1, 0, nsp-1, "ProdRate", true);

      for (int n = 0; n < np; n++) {
        vector_fp this_ProdRate = calcProdRate_i(Tvec[n], getMassFracVec(flow, n) );
        for (size_t k=0; k<nsp; k++)
          ProdRate[n][k] = this_ProdRate[k];
      }

      for (size_t k=0; k<nsp; k++) {
        const char *spName = mix->speciesName(k).c_str();
        string headerName = "ProdRate-"+ mix->speciesName(k)+" [kg/m^3s]";
        fprintf(FP,"%s\n", headerName.c_str());

        for (int n = 0; n < np; n++) {
          double this_pr = ProdRate[n][k] * mix->molecularWeight(k);
          const double tmp = ((abs(this_pr)<1e-99) ? 0.0 : this_pr);
          fprintf(FP,"\t%11.6e",tmp);
          if ( n%5==4 || n==np-1 ) fprintf(FP,"\n");
        }
      }

      fclose(FP);
      freeMem2D(ProdRate, 0, np-1, 0, nsp-1);
      delete[] rhovec;
      delete[] mw;
      delete[] mu;
      delete[] cp;
      delete[] lambda;
      delete[] h;
    }

    void writeScurve() {
      int n = aVec.size();
      FILE* FP = fopen("sCurve.dat", "w");
      fprintf(FP, "VARIABLES = a_inverse[s],  a[1/s],  Tmax[K]\n");

      for (int i=0; i<n; i++) {
        fprintf(FP, "%11.3e, %11.3e, %11.3e\n", 1.0/aVec[i], aVec[i], TmaxVec[i]);
      }
    }

};

#endif
