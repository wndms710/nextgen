#include "CtiLiquid.hpp"
#include "MiscUtils.hpp"

CtiLiquid::CtiLiquid(const string& name, double Pref) : CtiMaterial() {
  R_UNIVERSAL = 8314.472; // J/kmol/K
  this->name = name;
  this->Pref = Pref;
  info_prefix = "";
}

void CtiLiquid::setInfoPrefix(const string& s) {
  info_prefix = s+" ";
}

void CtiLiquid::info() {
  cout << "======================= liquid info ========================" << endl;
  cout << "CtiLiquid: " << name << endl;
  cout << "Pref [Pa]: " << Pref << endl;
  cout << "Tboil [K]: " << Tboil << endl;
  cout << "MW [kg/kmol]: " << MW << endl;
  cout << "density at 300K [kg/m3]: " << calcRho(300) << endl;
  cout << "=================== end of liquid info ======================" << endl;
}

void CtiLiquid::info(const double T0,const double T1,const int n = 25) {
  const int nset = 15;
  cout << "# CtiLiquid " << name << ": Pref=" << Pref << "[Pa], Tboil=" << Tboil << "[K], MW=" << MW << "[kg/kmol]" << endl;
  cout << setw(nset) << left << 
          "# name" << setw(nset) << left <<
          "T" << setw(nset) << left <<
          "Hvap" << setw(nset) << left <<   
          "rho" << setw(nset) << left <<     
          "rhov" << setw(nset) << left <<    
          "Cp" << setw(nset) << left <<         
          "Cpv" << setw(nset) << left <<        
          "kv" << setw(nset) << left <<       
          "muv" << setw(nset) << left <<      
          "Dv" << setw(nset) << left <<
          "Pv" << setw(nset) << left <<
          "sigma" << endl;
  cout << setw(nset) << left <<
          "# [none]" << setw(nset) << left <<
          "[K]" << setw(nset) << left <<
          "[J/kg]" << setw(nset) << left <<
          "[kg/m3]" << setw(nset) << left <<
          "[kg/m3]" << setw(nset) << left <<
          "[J/(kg*K)]" << setw(nset) << left <<
          "[J/(kg*K)]" << setw(nset) << left <<
          "[W/(m*K)]" << setw(nset) << left <<
          "[kg/(m*s)]" << setw(nset) << left <<
          "[m^2/s]" << setw(nset) << left <<
          "[Pa]" << setw(nset) << left <<
          "[N/m]" << endl;
  for (int i = 0; i < n; ++i) {
    const double T = T0 + double(i)/double(n-1)*(T1-T0);
    cout << setw(nset) << left <<
            info_prefix << setw(nset) << left <<
            T << setw(nset) << left << 
            calcHvap(T) << setw(nset) << left <<
            calcRho(T) << setw(nset) << left <<
            calcRhov(T) << setw(nset) << left <<
            calcCp(T) << setw(nset) << left <<
            calcCpv(T) << setw(nset) << left <<
            calcKv(T) << setw(nset) << left <<
            calcMuv(T) << setw(nset) << left <<
            calcDv(T) << setw(nset) << left <<
            calcPv(T) << setw(nset) << left <<
            calcSigma(T) << endl;
  }
}

YawsLiquid::YawsLiquid(const string& name, double Pref) : CtiLiquid(name, Pref) {

  if (mpi_rank == 0)
    cout << "YawsLiquid(): " << name << endl;

  if (Pref != 101325.0)
    cout << " > WARNING!!! Tboil is for Pref = 101325 Pa (not " << Pref << ")" << endl;

  if (name == "WATER_YAWS") {

    Tboil = 373.15;
    Tcrit = 647.3;
    MW = 18.02; // kg/kmol

    // density...
    Ar = 0.3471;
    Br = 0.274;
    nr = 0.2857;

    // latent heat of vaporization...
    AL = 52.023;
    nL = 0.321;

    // liquid heat capacity at constant pressure...
    AC = 92.053 ;
    BC = -0.039953 ;
    CC = -0.00021103 ;
    DC = 0.0000005347 ;

    // liquid vapor heat capacity at constant pressure...
    AV = 33.933 ;
    BV = -0.0084186 ;
    CV = 0.000029906 ;
    DV = -0.000000017825 ;
    EV = 3.6934E-12 ;

    // vapor thermal conductivity...
    Ak = 0.00053 ;
    Bk = 0.000047093 ;
    Ck = 0.000000049551 ;

    // vapor pressure...
    AP = 29.8605 ;
    BP = -3152.2 ;
    CP = -7.3037 ;
    DP = 0.00000024247 ;
    EP = 0.000001809 ;

    // sigma...
    Asigma = 6.26E-2; // N/m
    Tsigma = 353.15;

  }
  else if (name == "JET-A_YAWS") {

    Tboil = 443.83;
    Tcrit = 737.0;
    MW = 181.321; // kg/kmol

    // density...
    Ar = 0.29292;
    Br = 0.26661;
    nr = 0.298;

    // latent heat of vaporization...
    AL = 73.509;
    nL = 0.347;

    // liquid heat capacity at constant pressure...
    AC = 142.238;
    BC = 1.5261;
    CC = -0.0034477;
    DC = 0.0000032968;

    // liquid vapor heat capacity at constant pressure...
    AV = -128.032 ;
    BV = 1.4622 ;
    CV = -0.00086193 ;
    DV = 0.00000018462 ;
    EV = 0.00000000000036227 ;

    // vapor thermal conductivity...
    Ak = -0.01184 ;
    Bk = 0.000061839 ;
    Ck = 0.000000025082 ;

    // vapor pressure...
    AP = -50.5512 ;
    BP = -2705.3 ;
    CP = 28.273 ;
    DP = -0.045702 ;
    EP = 0.000020443 ;

    // sigma: NIST data for C12H26 dodecane
    Asigma = 0.024752;
    Tsigma = 300.0;

  }
  else {

    CERR("unrecognized Yaws liquid: " << name);

  }

}

double YawsLiquid::calcRho(const double T) {
  // density in kg/m3...
  return( 1000.0*Ar*pow(Br,-pow(1.0-T/Tcrit,nr)) );
}

double YawsLiquid::calcHvap(const double T) {
  // latent heat of vaporization in kJ/mol...
  const double t = min(1.0,T/Tcrit);
  const double molar_value = AL*pow(1.0-t,nL);
  return( molar_value*1.0E+6/MW ); // convert to J/kg
}

double YawsLiquid::calcCp(const double T) {
  // heat capacity at constant pressure in J / (mol*K)
  const double molar_value = AC + BC*T + CC*T*T + DC*T*T*T;
  return( molar_value*1000.0/MW ); // convert to J/(kg*K)
}

double YawsLiquid::calcCpv(const double T) {
  // heat capacity at constant pressure in J / (mol*K)
  const double molar_value =  AV + BV*T + CV*T*T + DV*T*T*T + EV*T*T*T*T;
  return( molar_value*1000.0/MW ); // convert to J/(kg*K)
}

double YawsLiquid::calcKv(const double T) {
  // thermal conductivity in W/(m*K)
  return( Ak + Bk*T + Ck*T*T );
}

double YawsLiquid::calcPv(const double T) {
  // vapor pressure in Pa (original fit in mmHg (or Torr)
  const double log10Pv = AP + BP/T + CP*log10(T) + DP*T + EP*T*T;
  return( 101325.0/760.0*pow(10.0,log10Pv) ); // converted to Pa by Pa = mmHg*101325/760
}

double YawsLiquid::calcSigma(const double T) {
  return( max(0.0,Asigma*(Tcrit - T)/(Tcrit - Tsigma)) );
}
double YawsLiquid::calcDv(const double T) {
  return( calcKv(T)/calcCpv(T)/calcRhov(T) );
  // NOTE: assuming Le = 1 -> Pr = Sc
}
double YawsLiquid::calcMuv(const double T) {
  return( 1.0e-5 ); // Warning!!
}
double YawsLiquid::calcRhov(const double T) {
  return( Pref*MW/R_UNIVERSAL/T );
}

MHB98Liquid::MHB98Liquid(const string& name, double Pref) : CtiLiquid(name, Pref) {

  if (mpi_rank == 0)
    cout << "MHB98Liquid(): " << name << endl;

  if (Pref != 101325.0)
    cout << " > WARNING!!! Tboil is for Pref = 101325 Pa (not " << Pref << ")" << endl;

  if (name == "WATER_MHB98") {

    Tboil = 373.15;
    MW = 18.015; // kg/kmol

    // Liquid density...
    rho = 997.0;

    // latent heat of vaporization...
    AL = 2.257E6;
    BL = 2.595E3;
    CL = 373.15;
    DL = 1.0;
    nL = 1;

    // liquid heat capacity at constant pressure...
    cp = 4184.0 ;

    // liquid vapor heat capacity at constant pressure...
    AV = 8137.0 ;
    BV = -37.34;
    CV = 0.07482 ;
    DV = -4.956E-5 ;

    // vapor thermal conductivity...
    Ak = 1.024E-2 ;
    Bk = -8.21E-6 ;
    Ck = 1.41E-7 ;
    Dk = -4.51E-11;

    // vapor viscosity
    Am = -3.077E-6;
    Bm = 4.07E-8;

  }
  else if (name == "HEXANE_MHB98") {

    Tboil = 344.6;
    MW = 86.178; // kg/kmol

    // density...
    rho = 664.0;

    // latent heat of vaporization...
    AL = 0.0;
    BL = 5.1478E5;
    CL = 1.0;
    DL = 512.0;
    nL = 0.3861;

    // liquid heat capacity at constant pressure...
    cp = 2302.0 ;

    // liquid vapor heat capacity at constant pressure...
    AV = -51.31 ;
    BV = 6.767;
    CV = -3.626E-3 ;
    DV = 0.0 ;

    // vapor thermal conductivity...
    Ak = 1.112E-2 ;
    Bk = 3.837E-5;
    Ck = 3.778E-8 ;
    Dk = 0.0;

    // vapor viscosity
    Am = 5.592E-6;
    Bm = 5.622E-9;

 }
  else {
    CERR("unrecognized MHB98 liquid: " << name);
  }

}

double MHB98Liquid::calcRho(const double T) {
  // density in kg/m3... 
  return( rho );
}

double MHB98Liquid::calcHvap(const double T) {
  // latent heat of vaporization in J/Kg...
  double T_mod = min(T,Tboil);
  Hvap = AL+BL*pow((CL-T_mod/DL),nL);
  return( Hvap ); //  J/kg
}

double MHB98Liquid::calcCp(const double T) {
  // heat capacity at constant pressure 
  return( cp ); // convert to J/(kg*K)
}

double MHB98Liquid::calcCpv(const double T) {
  // heat capacity of liquid vapor at constant pressure in J / (mol*K)
  cpv = AV + BV*T+CV*T*T + DV*T*T*T;
  return( cpv ); // convert to J/(kg*K) ????
}

double MHB98Liquid::calcKv(const double T) {
  // thermal conductivity of liquid vapor  in W/(m*K)
  kv = Ak + Bk*T + Ck*T*T + Dk*T*T*T;
  return( kv );
}

double MHB98Liquid::calcDv(const double T ){
  // For water and hexane, Lewis number is assumed 1.0 and Diffusive coefficient is evaluated
  // based on Prandtl number of air 
  double mua =  6.109E-6 + 4.604E-8*T-1.051E-11*T*T;
  double Prg;
  if (T < 600 )
    Prg = 0.815-4.958E-4*T+4.514E-7*T*T;
  else
    Prg = 0.647+5.5E-5*T;

  double rhoa = ( 355.91 * pow(T,-1.0032) ); // kg/m^3

  return(mua/(rhoa*Prg)); // Lewis number = 1 assumed for water and hexane
}

double MHB98Liquid::calcPv(const double T) {
  return(101325.0*exp(calcHvap(T)*MW/R_UNIVERSAL*(1.0/Tboil-1.0/T)));
}

double MHB98Liquid::calcSigma(const double T) {
  return( 0.02 ); // Warning !!
}

double MHB98Liquid::calcMuv(const double T) {
  return( Am + Bm*T );
}

double MHB98Liquid::calcRhov(const double T) {
  return( Pref*MW/R_UNIVERSAL/T );
}

MHB98heptane::MHB98heptane(const string& name, double Pref) : CtiLiquid(name, Pref) { 

  if (mpi_rank == 0)
    cout << "MHB98heptane(): " << name << endl;

  //if (Pref != 101325.0)
  //  cout << " > WARNING!!! Tboil is for Pref = 101325 Pa (not " << Pref << ")" << endl;

  if (name == "HEPTANE_MHB98") {

    //Tboil = 371.6; //ok
    Tboil = 437.0; // for 0.5 MPa
    MW = 100.0; // kg/kmol //ok

    // Liquid density...
    rho = 649.38; //ok

    // latent heat of vaporization...
    AL = 0.0; //ok
    BL = 3.163E5;
    CL = 3.204;
    DL = 168.6;
    nL = 0.38;

    // liquid heat capacity at constant pressure...
    cp = 2383.89 ; //ok

    // liquid vapor heat capacity at constant pressure...
    AV = -51.56 ; //ok
    BV = 6.776;
    CV = -3.658E-3 ;
    DV = -7.673E-7 ;

    // vapor thermal conductivity...
    Ak = -4.401E-2 ; //ok
    Bk = 2.514E-4 ;
    Ck = -3.173E-7 ;
    Dk = 2.487E-10;

    // vapor viscosity
    Am = 3.83E-6; //ok
    Bm = -3.613E-9;
    Cm = 4.911E-11;
    Dm = -3.577E-14;

    // vapor mass diffusivity
    Ad = 5.94E-6;
    Bd = 1.0;
    Cd = 273.0;
    Dd = 1.6;
  } else {
    CERR("unrecognized MHB98 liquid: " << name);
  }
}

double MHB98heptane::calcRho(const double T) {
  // density in kg/m3... 
  return( rho );
}

double MHB98heptane::calcHvap(const double T) {
  // latent heat of vaporization in J/Kg...
  const double T_mod = min(T,Tboil);
  Hvap = AL+BL*pow((CL-T_mod/DL),nL);
  return( Hvap ); //  J/kg
}
double MHB98heptane::calcCp(const double T) {
  // heat capacity at constant pressure 
  return( cp ); // convert to J/(kg*K)
}

double MHB98heptane::calcCpv(const double T) {
  // heat capacity of liquid vapor at constant pressure in J / (kg*K)
  cpv = AV + BV*T+CV*T*T + DV*T*T*T;
  return( cpv );
}

double MHB98heptane::calcKv(const double T) {
  // thermal conductivity of liquid vapor  in W/(m*K)
  kv = Ak + Bk*T + Ck*T*T + Dk*T*T*T;
  return( kv );
}

double MHB98heptane::calcDv(const double T ){
  // mass diffusivity of liquid vapor in m^2/s
  Dv = Ad*pow(Bd*T/Cd,Dd)/(Pref/101325.0);
  return( Dv ); // Lewis number = 1 assumed for water and hexane
}

double MHB98heptane::calcPv(const double T) {
  //return(101325.0*exp(calcHvap(T)*MW/R_UNIVERSAL*(1.0/Tboil-1.0/T)));
  return(Pref*exp(calcHvap(T)*MW/R_UNIVERSAL*(1.0/Tboil-1.0/T)));
}

double MHB98heptane::calcSigma(const double T) {
  return( 0.02 ); // Warning !!
}

double MHB98heptane::calcMuv(const double T) {
  return( Am + Bm*T + Cm*T*T + Dm*T*T*T);
}

double MHB98heptane::calcRhov(const double T) {
  return( Pref*MW/R_UNIVERSAL/T );
}

MHB98decane::MHB98decane(const string& name, double Pref) : CtiLiquid(name, Pref) {

  if (mpi_rank == 0)
    cout << "MHB98decane(): " << name << endl;

  if (Pref != 101325.0)
    cout << " > WARNING!!! Tboil is for Pref = 101325 Pa (not " << Pref << ")" << endl;

  Tboil = 447.7;
  MW = 142.0; // kg/kmol

}

double MHB98decane::calcRho(const double T) {
  // density in kg/m3... 
  return(642.0);
}

double MHB98decane::calcHvap(const double T) {
  // latent heat of vaporization in J/Kg...
  double T_mod = min(T,Tboil);
  return( 3.958E4*pow((619.0-T_mod),0.38) ); //  J/kg
}

double MHB98decane::calcCp(const double T) {
  // heat capacity at constant pressure 
  return(2520.5 ); // convert to J/(kg*K)
}

double MHB98decane::calcCpv(const double T) {
  // heat capacity of liquid vapor at constant pressure in J / (mol*K)
  double AV,BV,CV,DV;
  if (T/1000.0 < 0.8 ) {
    AV = 106.6 ;
    BV = 5765.0;
    CV = -1675.0 ;
    DV = 473.1; }
  else {
    AV = 411.1;
    BV = 5460.0;
    CV = -2483.0;
    DV = 422.9; }

  double Tstar = T/1000.0;

  return(AV + BV*Tstar+CV*Tstar*Tstar + DV*Tstar*Tstar*Tstar ); // convert to J/(kg*K)
}

double MHB98decane::calcKv(const double T) {
  // thermal conductivity of liquid vapor in W/(m*K)
  return(1.124E-2*pow((T/300.0),1.8));
}

double MHB98decane::calcDv(const double T){
  return(5.46E-6*pow((T/300.0),1.583));
}

double MHB98decane::calcPv(const double T) {
  return(101325.0*exp(calcHvap(T)*MW/R_UNIVERSAL*(1.0/Tboil-1.0/T)));
}

double MHB98decane::calcSigma(const double T) {
  return(23.83E-3); // Warning !!
}

double MHB98decane::calcMuv(const double T) {
  return(5.64E-6+1.75E-8*(T-300));
}

double MHB98decane::calcRhov(const double T) {
  return( Pref*MW/R_UNIVERSAL/T );
}

VLSATLiquid::VLSATLiquid(const string& name, double Pref) : CtiLiquid(name, Pref) {

  if (mpi_rank == 0)
    cout << "VLSATLiquid(): " << name << endl;

  if (name == "WATER_VLSAT") {

    // saturated liquid/vapor data from NIST (REFPROP)
    // valid for 275K < T < Tcrit , 0.007 bar < P < Pcrit

    MW = 18.0153; // molecular weight [kg/kmol]

    // critical properties
    Pcrit = 2.2064e+07;
    Tcrit = 647.096;
    rho_c = 322;

    // vapor pressure
    PV0 = -6.65491;
    PV1 = -0.0549419;

    // heat of evaporation
    HV0 = 2.89961e+06;
    HV1 = 0.168367;
    HV2 = 0.402826;

    // surface tension
    ST0 = 0.102838;
    ST1 = 0.977705;
    ST2 = 1.30648;

    // liquid density
    DL0 = 799.064;
    DL1 = 0.187675;
    DL2 = 0.401723;

    // liquid heat capacity
    CL0 = 6.38117;
    CL1 = -0.0486798;
    CL2 = -0.13008;
    CL3 = -0.262305;
    CL4 = 0.3628;

    // vapor density
    DV0 = 0.089767;
    DV1 = 0.324282;

    // vapor heat capacity
    CV0 = 7.42043;
    CV1 = 0.125143;
    CV2 = -0.0343777;
    CV3 = 0.0566652;
    CV4 = 0.0290291;

    // vapor conductivity
    KV0 = 0.0247523;
    KV1 = 0.718914;
    KV2 = -0.516384;

    // vapor viscosity
    VV0 = 4.66541e-05;
    VV1 = 0.814516;
    VV2 = -0.0775211;
    VV3 = 0.164386;

  } else if (name == "DECANE_VLSAT") {

    // saturated liquid/vapor data from NIST (REFPROP)
    // valid for 290K < T < Tcrit , 0.001 bar < P < Pcrit

    MW = 142.2817; // molecular weight [kg/kmol]

    // critical properties
    Pcrit = 2.103e+06;
    Tcrit = 617.7;
    rho_c = 233.34;

    // vapor pressure
    PV0 = -6.8082;
    PV1 = -0.10566;

    // heat of evaporation
    HV0 = 458105;
    HV1 = 0.0816545;
    HV2 = 0.434169;

    // surface tension
    ST0 = 0.0534539;
    ST1 = 0.000509676;
    ST2 = 1.25993;

    // liquid density
    DL0 = 666.436;
    DL1 = -0.0846286;
    DL2 = 0.404108;

    // liquid heat capacity
    CL0 = 5.74942;
    CL1 = 0.247807;
    CL2 = -0.145615;
    CL3 = -0.617016;
    CL4 = 0.707891;

    // vapor density
    DV0 = 0.117082;
    DV1 = 0.369648;

    // vapor heat capacity
    CV0 = 5.6758;
    CV1 = 0.0680354;
    CV2 = -0.0752261;
    CV3 = -0.492396;
    CV4 = 0.593234;

    // vapor conductivity
    KV0 = 0.0279638;
    KV1 = 1.63483;
    KV2 = -0.170038;

    // vapor viscosity
    VV0 = 2.932e-05;
    VV1 = 0.871425;
    VV2 = -0.0440559;
    VV3 = 0.139206;

  } else if (name == "DODECANE_VLSAT") {

    // saturated liquid/vapor data from NIST (REFPROP)
    // valid for 265K < T < Tcrit , 7e-6 bar < P < Pcrit

    MW = 170.3348; // molecular weight [kg/kmol]

    // critical properties
    Pcrit = 1.8176e+06;
    Tcrit = 658.1;
    rho_c = 226.55;

    // vapor pressure
    PV0 = -7.09559;
    PV1 = -0.115611;

    // heat of evaporation
    HV0 = 503310;
    HV1 = -0.236108;
    HV2 = 0.371665;

    // surface tension
    ST0 = 0.0524489;
    ST1 = 0.385223;
    ST2 = 1.52454;

    // liquid density
    DL0 = 737.627;
    DL1 = -0.308591;
    DL2 = 0.347158;

    // liquid heat capacity
    CL0 = 6.47165;
    CL1 = 0.170849;
    CL2 = -0.0793192;
    CL3 = -0.40694;
    CL4 = 0.469055;

    // vapor density
    DV0 = 0.0956629;
    DV1 = 0.35027;

    // vapor heat capacity
    CV0 = 6.11236;
    CV1 = 0.0668121;
    CV2 = -0.0364225;
    CV3 = -0.420555;
    CV4 = 0.49785;

    // vapor conductivity
    KV0 = 0.0348778;
    KV1 = 1.83907;
    KV2 = -0.140195;

    // vapor viscosity
    VV0 = 1.8929e-05;
    VV1 = 0.76997;
    VV2 = -0.13408;
    VV3 = 0.189985;

  } else if (name == "JETFUEL_VLSAT") {

    // jet-A/jet-A1/JP8 surrogate
    // https://doi.org/10.1016/j.combustflame.2015.12.013
    // 
    // species            name       mol%
    // ----------         --------   ----
    // n-dodecane         N-C12H26   30.3
    // m-xylene           A1CH3CH3   21.2
    // methylcyclohexane  MCH-C7H14  48.5
    //
    // saturated liquid/vapor data from PR-EOS (Kij = 0) 
    // valid for 218K < T < Tcrit and 1e-5 bar < P < Pcrit

    MW = 121.741; // molecular weight [kg/kmol]

    // critical properties
    Pcrit = 2.67887e+06;
    Tcrit = 606.359;
    rho_c = 201.283;

    // vapor pressure
    PV0 = -6.60763;
    PV1 = -0.0732803;

    // heat of evaporation
    HV0 = 423450;
    HV1 = 0.269276;
    HV2 = 0.451625;

    // surface tension
    ST0 = 0.0518332;
    ST1 = 0.40444;
    ST2 = 1.5221;

    // liquid density
    DL0 = 631.632;
    DL1 = 0.354879;
    DL2 = 0.507017;

    // liquid heat capacity
    CL0 = 6.86592;
    CL1 = -0.0482263;
    CL2 = 0.0300829;
    CL3 = -0.171197;
    CL4 = 0.226898;

    // vapor density
    DV0 = 0.160749;
    DV1 = 0.484891;

    // vapor heat capacity
    CV0 = 6.21319;
    CV1 = -0.106451;
    CV2 = 0.0446096;
    CV3 = -0.30314;
    CV4 = 0.37452;

    // vapor conductivity
    KV0 = 0.0330895;
    KV1 = 1.9179;
    KV2 = -0.126413;

    // vapor viscosity
    VV0 = 1.14239e-05;
    VV1 = 0.812368;
    VV2 = -0.0932924;
    VV3 = 0.311746;

  }

  // invert calcPv for Tboil
  Tboil = Tcrit*(PV0-PV1*log(Pref/Pcrit))/(PV0+log(Pref/Pcrit));
}

double VLSATLiquid::calcPv(const double T) {
  // vapor pressure [Pa]
  double t = min(1.0,T/Tcrit);
  return( Pcrit*exp(PV0*(1.0-t)/(t+PV1)) );
  // NOTE: Pv = Pcrit when T/Tcrit = 1
}

double VLSATLiquid::calcHvap(const double T) {
  // heat of evaporation [J/kg]
  double t = min(1.0,T/Tcrit);
  return( HV0*exp(HV1*t)*pow(1.0-t,HV2) );
  // NOTE: dHv = 0 when T/Tcrit = 1
}

double VLSATLiquid::calcSigma(const double T) {
  // liquid surface tension [N/m]
  double t = min(1.0,T/Tcrit);
  return( ST0*exp(ST1*t)*pow(1.0-t,ST2) );
  // NOTE: sigma = 0 when T/Tcrit = 1
}

double VLSATLiquid::calcRho(const double T) {
  // liquid density [kg/m3]
  double t = min(1.0,T/Tcrit);
  return( DL0*exp(DL1*t)*pow(1.0-t,DL2)+rho_c );
  // NOTE: rho = rho_c when T/Tcrit = 1
}

double VLSATLiquid::calcRhov(const double T) {
  // vapor density [kg/m3]
  double t = min(1.0,T/Tcrit);
  double p = min(1.0,calcPv(T)/Pcrit);
  return( rho_c*(1.0-exp(DV0*p/t)*pow(1.0-p/t,DV1)) ); 
  // NOTE: rhov = rho_c when T/Tcrit = 1

  //return( Pref*MW/R_UNIVERSAL/T ); // ideal gas
}

double VLSATLiquid::calcCp(const double T) {
  // liquid heat capacity [J/kg/K]
  double t = min(.999999,T/Tcrit);
  return( exp(CL0*pow(t,CL1*t+CL2)/pow(1.0-t,CL3*t+CL4)) );
  // NOTE: cp = Inf when T/Tcrit = 1
}

double VLSATLiquid::calcCpv(const double T) {
  // vapor heat capacity [J/kg/K]
  double t = min(.999999,T/Tcrit);
  return( exp(CV0*pow(t,CV1*t+CV2)/pow(1.0-t,CV3*t+CV4)) );
  // NOTE: cpv = Inf when T/Tcrit = 1
}

double VLSATLiquid::calcKv(const double T ){
  // vapor conductivity [W/m/K]
  double t = min(.999999,T/Tcrit);
  return( KV0*pow(t,KV1)*pow(1.0-t,KV2) );
  // NOTE: kv = Inf when T/Tcrit = 1
}

double VLSATLiquid::calcMuv(const double T ){
  // vapor viscosity [kg/m/s]
  double t = min(1.0,T/Tcrit);
  return( VV0*(1.0-VV1*pow(t,VV2)*pow(1.0-t,VV3)) );
}

double VLSATLiquid::calcDv(const double T ){
  // vapor mass diffusivity [m2/s]
  return( calcKv(T)/calcCpv(T)/calcRhov(T) );
  // NOTE: assuming Le = 1 -> Pr = Sc (Dv -> 0 at Tcrit)
}

UserLiquid::UserLiquid(const string& name, double Pref) : CtiLiquid(name, Pref) {

  if (mpi_rank == 0)
    cout << "UserLiquid(): " << name << endl;

  MW    = getDoubleParam(name+".MW");
  rho   = getDoubleParam(name+".RHO");
  Hvap  = getDoubleParam(name+".HVAP");
  cp    = getDoubleParam(name+".CP");
  cpv   = getDoubleParam(name+".CPV");
  kv    = getDoubleParam(name+".KV");
  pv    = getDoubleParam(name+".PV");
  muv   = getDoubleParam(name+".MUV");
  sigma = getDoubleParam(name+".SIGMA");
  Tboil = getDoubleParam(name+".TBOIL");

}

double UserLiquid::calcRho(const double T) {
  return(rho);
}

double UserLiquid::calcHvap(const double T) {
  return(Hvap);
}

double UserLiquid::calcCp(const double T) {
  return(cp);
}

double UserLiquid::calcCpv(const double T) {
  return(cpv);
}

double UserLiquid::calcKv(const double T) {
  return(kv);
}

double UserLiquid::calcPv(const double T) {
  return(pv);
}

double UserLiquid::calcSigma(const double T) {
  return(sigma);
}

double UserLiquid::calcDv(const double T) {
  return(kv/cpv/calcRhov(T));
}

double UserLiquid::calcMuv(const double T) {
  return(muv);
}

double UserLiquid::calcRhov(const double T) {
  return( Pref*MW/R_UNIVERSAL/T );
}

CtiLiquid * newCtiLiquid(Param * param,int &iarg) {

  // first token should be a recognized name...

  if (iarg >= param->size()) {
    CERR("LP.MATERIAL LIQUID missing id. Possible choices include:\n" << 
         "USER, WATER_YAWS, JET-A_YAWS, WATER_MHB98, HEXANE_MHB98, HEPTANE_MHB98, DECANE_MHB98, WATER_VLSAT, DECANE_VLSAT, DODECANE_VLSAT, JETFUEL_VLSAT");
  }
  const string id = param->getString(iarg++);
  string name = id;

  // parse anything else...
  
  //bool b_Pref = false;
  double Pref = 101325;
  while (iarg < param->size()) {
    string token = MiscUtils::toUpperCase(param->getString(iarg++));
    if (token == "P_REF") {
      //b_Pref = true;
      Pref = param->getDouble(iarg++);
      if (mpi_rank == 0) cout << " > P_REF " << Pref << endl;
    }
    else if (token == "NAME") {
      name = param->getString(iarg++);
    }
    else {
      if (mpi_rank == 0) cout << " > skipping unrecognized LP.MATERIAL token: " << token << endl;
    }
  }
  
  if (id == "USER")
    return(new UserLiquid(name,Pref));
  else if (id == "WATER_YAWS" || id == "JET-A_YAWS")
    return(new YawsLiquid(name,Pref) );
  else if (id == "WATER_MHB98" || id=="HEXANE_MHB98")
    return(new MHB98Liquid(name,Pref) );
  else if (id == "HEPTANE_MHB98")
    return(new MHB98heptane(name,Pref) );
  else if (id == "DECANE_MHB98")
    return(new MHB98decane(name,Pref) );
  else if (id == "WATER_VLSAT" || "DECANE_VLSAT" || id == "DODECANE_VLSAT" || id == "JETFUEL_VLSAT")
    return(new VLSATLiquid(name,Pref) );
  else {
    CERR("unrecognized LP.MATERIAL id: " << id << ". Possible choices include:\n" << 
         "USER, WATER_YAWS, JET-A_YAWS, WATER_MHB98, HEXANE_MHB98, HEPTANE_MHB98, DECANE_MHB98, WATER_VLSAT, DECANE_VLSAT, DODECANE_VLSAT, JETFUEL_VLSAT");
  }

  // should never get here...
  return NULL;
}
