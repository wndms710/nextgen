#include "candle.hpp"

int premixed1() {
  const string fuel = "CH4:.95799, C2H6:.02683, C3H8:.00289, CO2:.00779, N2:.00450";
  const string oxid = "O2:0.21 N2:0.79";
  const string mech = "gri30.xml";
  const string tran = "Mix";

  const double Lewis = 1.0;
  const double p     = 20.739e5; // 300.8 psia
  const double tf    = 483.15;   // 410 F
  const double to    = 729.26;   // 853 F

  int nfl;
  int loglev = 0;
  double phi1, phi2;
  std::cout << "enter phi1, phi2, nfl: ";
  std::cin >> phi1;
  std::cin >> phi2;
  std::cin >> nfl;
  std::cout << endl;
  const bool writeTecplot = false;
  const bool createInitFlamelet = true;

  Premixed f(fuel, oxid, p, tf, to, mech, tran, Lewis, loglev, createInitFlamelet, "init_p23_0", writeTecplot);
  return f.createFlamelets(phi1, phi2, nfl);
}

int diffusion() {
  const double Lewis      = 1.0;
  const int loglevel      = 0;//1;
  const bool createInitFlamelet = false;
  const bool writeTecplot = true;
  const string fCompo     = "CH4:0.25 O2:0.1575 N2:0.5925";
  const double Tf         = 294.0;
  const string oxCompo    = "O2:0.21 N2:0.79";
  const double Tox        = 294.0;
  const double pc         = 0.993 * 101325.0;

//    Premixed f("CH4:1", "O2:1 N2:3.76", 1.0e5, 300.0, 300.0, "gri30.xml", "Mix", Lewis, 0);
//    return f.createFlamelets(0.65, 0.85, 9);

  Diffusion f(fCompo, oxCompo, pc, Tf, Tox, "gri30.xml", "Mix", Lewis, loglevel, createInitFlamelet, 31, "init_D_a200.62669", 200.62669, writeTecplot);
  return f.createFlamelets(400);

}

int premixed2() {
  const string fuel = "CH4:0.936, C2H6:0.043, C3H8:0.003, CO2:0.011, N2:0.007";
//  const string fuel = "H2:0.1, CH4:0.8424, C2H6:0.0387, C3H8:0.0027, CO2:0.0099, N2:0.0063"
//  const string fuel = "H2:0.2, CH4:0.67392, C2H6:0.03096, C3H8:0.00216, CO2:0.00792, N2:0.00504";
//  const string fuel = "H2:0.3, CH4:0.471744, C2H6:0.021672, C3H8:0.001512, CO2:0.005544, N2:0.003528";
//  const string fuel = "H2:0.4, CH4:0.2830464, C2H6:0.0130032, C3H8:0.0009072, CO2:0.0033264, N2:0.0021168";
//  const string fuel = "H2:0.5, CH4:0.1415232, C2H6:0.0065016, C3H8:0.0004536, CO2:0.0016632, N2:0.0010584";
//  const string fuel = "H2:0.6, CH4:0.05660928, C2H6:0.00260064, C3H8:0.00018144, CO2:0.00066528, N2:0.00042336";

  const string oxid = "O2:0.21 N2:0.79";
  const string mech = "gri30.xml";
  const string tran = "Mix";

  const double Lewis = 1.0;
  const double p     = 23.40e5;  // 339.4 psia
  const double tf    = 499.82;   // 440.0 F
  const double to    = 754.15;   // 897.8 F

  int nfl = 200;
  const bool writeTecplot = true;
  const bool createInitFlamelet = true;//false;
  int loglev = 1;
  double phi1, phi2;
  phi1 = 0.4;
  phi2 = 1.5;

  Premixed f(fuel, oxid, p, tf, to, mech, tran, Lewis, loglev, createInitFlamelet, "init_p23_0", writeTecplot);
  return f.createFlamelets(phi1, phi2, nfl);
}

int premixed3() {
  const string fuel = "CH4:1.0";

  const string oxid = "O2:0.21 N2:0.79";
  const string mech = "gri30.xml";
  const string tran = "Mix";

  const double Lewis = 1.0;
  const double p     = 24.54e5;  // 355.9 psia
  const double tf    = 499.82;   // 440.0 F
  const double to    = 754.15;   // 897.8 F

  int nfl = 200;
  const bool writeTecplot = true;
  const bool createInitFlamelet = true;//false;
  int loglev = 0;
  double phi1, phi2;
  phi1 = 0.4;
  phi2 = 1.5;

  Premixed f(fuel, oxid, p, tf, to, mech, tran, Lewis, loglev, createInitFlamelet, "init_p23_0", writeTecplot);
  return f.createFlamelets(phi1, phi2, nfl);
}

int premixed4() {
//  const string fuel = "CH4:0.75, H2:0.25";
  const string fuel = "CH4:0.5, H2:0.5";

  const string oxid = "O2:0.21 N2:0.79";
  const string mech = "gri30.xml";
  const string tran = "Mix";

  const double Lewis = 1.0;
  const double p     = 16.20e5;  // 235 psia
  const double tf    = 458.15;   // 365 F
  const double to    = 688.71;   // 780 F

  int nfl = 200;
  const bool writeTecplot = true;
  const bool createInitFlamelet = true;//false;
  int loglev = 0;
  double phi1, phi2;
  phi1 = 0.4;
  phi2 = 1.5;

  Premixed f(fuel, oxid, p, tf, to, mech, tran, Lewis, loglev, createInitFlamelet, "init_p16_25", writeTecplot);
  return f.createFlamelets(phi1, phi2, nfl);
}

int diffusion4() {
  const string fuel = "CH4:0.75, H2:0.25";
//  const string fuel = "CH4:0.5, H2:0.5";

  const string oxid = "O2:0.21 N2:0.79";
  const string mech = "gri30.xml";
  const string tran = "Mix";

  const double Lewis = 1.0;
  const double p     = 16.20e5;  // 235 psia
  const double tf    = 458.15;   // 365 F
  const double to    = 688.71;   // 780 F

  int nfl = 200;
  const bool writeTecplot = true;
  const bool createInitFlamelet = true;//false;
  int loglev = 0;
  double phi1, phi2;
  phi1 = 0.4;
  phi2 = 1.5;

  Diffusion f(fuel, oxid, p, tf, to, mech, tran, Lewis, loglev, createInitFlamelet, 31, "init_p16_0.5_diff", 0.5, writeTecplot);
  return f.createFlamelets(400);

}

int diffusion2() {
  const string fuel = "CH4:0.936, C2H6:0.043, C3H8:0.003, CO2:0.011, N2:0.007";
//  const string fuel = "H2:0.1, CH4:0.8424, C2H6:0.0387, C3H8:0.0027, CO2:0.0099, N2:0.0063"
//  const string fuel = "H2:0.2, CH4:0.67392, C2H6:0.03096, C3H8:0.00216, CO2:0.00792, N2:0.00504";
//  const string fuel = "H2:0.3, CH4:0.471744, C2H6:0.021672, C3H8:0.001512, CO2:0.005544, N2:0.003528";
//  const string fuel = "H2:0.4, CH4:0.2830464, C2H6:0.0130032, C3H8:0.0009072, CO2:0.0033264, N2:0.0021168";
//  const string fuel = "H2:0.5, CH4:0.1415232, C2H6:0.0065016, C3H8:0.0004536, CO2:0.0016632, N2:0.0010584";
//  const string fuel = "H2:0.6, CH4:0.05660928, C2H6:0.00260064, C3H8:0.00018144, CO2:0.00066528, N2:0.00042336";

  const string oxCompo = "O2:0.21 N2:0.79";
  const string mech    = "gri30.xml";
  const string tran    = "Mix";

  const double Lewis = 1.0;
  const double p     = 23.164e5; // 335.9696 psia
  const double Tf    = 499.82;   // 440.0 F
  const double Tox   = 754.15;   // 897.8 F


  const int loglevel            = 0;
  const bool createInitFlamelet = false;
  const bool writeTecplot       = true;

  Diffusion f(fuel, oxCompo, p, Tf, Tox, "gri30.xml", "Mix", Lewis, loglevel, createInitFlamelet, 31, "init_p23_0_diff", 0.5, writeTecplot);
  return f.createFlamelets(400);

}

int diffusion5() {
  const string fuel = "CH4:1.0";

  const string oxCompo = "O2:1 N2:3.76";
  const string mech    = "gri30.xml";
  const string tran    = "Mix";

  const double Lewis = 1.0;
  const double p     = 17.24e5;  // 250 psia
  const double Tf    = 483.0;    // 410.0 F
  const double Tox   = 729.0;    // 853 F


  const int loglevel            = 1;
  const bool createInitFlamelet = true;
  const bool writeTecplot       = true;

  Diffusion f(fuel, oxCompo, p, Tf, Tox, "gri30.xml", "Mix", Lewis, loglevel, createInitFlamelet, 31, "init_p17_diff", 0.5, writeTecplot);
  return f.createFlamelets(400);

}

int main() {
//  premixed1();
  premixed2();
////  diffusion5();
//  diffusion2();
}
