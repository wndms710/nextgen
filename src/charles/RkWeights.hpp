#ifndef RKWEIGHTS_HPP
#define RKWEIGHTS_HPP

namespace RkWeights {

  const double erk_wgt3[3][3] = {{1.0,0.0,0.0},
                               {-0.75,0.25,0.0},
                               {-1.0/12.0,-1.0/12.0,8.0/12.0}};

  const double erk_wgt_sum3[3] = {
    1.0,
    -0.5,
    0.5
  };

  const double erk_c3[3] = {
    0.0,
    1.0,
    0.5
  };

#define ALPHA2 0.06700195705857155
#define ALPHA3 0.4504266495780916
#define ALPHA4 0.5166954721912165
#define ALPHA5 0.3961936038613332
#define ALPHA6 0.04623353107167757
#define ALPHA7 0.1981318996397357
  
  const double erk_wgt_stb53[5][5] = {
    { 1.0, 0.0, 0.0, 0.0, 0.0 },
    { ALPHA2-1.0, ALPHA2, 0.0, 0.0, 0.0 },
    { ALPHA2*(ALPHA3-1.0), ALPHA2*(ALPHA3-1.0), ALPHA3, 0.0, 0.0 },
    { ALPHA2*ALPHA3*(ALPHA4-1.0), ALPHA2*ALPHA3*(ALPHA4-1.0), ALPHA3*(ALPHA4-1.0), ALPHA4, 0.0 },
    { ALPHA2*(ALPHA7+ALPHA3*(ALPHA6+ALPHA4*(ALPHA5-1.0))), ALPHA2*(ALPHA7+ALPHA3*(ALPHA6+ALPHA4*(ALPHA5-1.0))), ALPHA7+ALPHA3*(ALPHA6+ALPHA4*(ALPHA5-1.0)), ALPHA6+ALPHA4*(ALPHA5-1.0), ALPHA5 }
  };

  const double erk_wgt_sum_stb53[5] = {
    1.0,
    erk_wgt_stb53[1][0]+erk_wgt_stb53[1][1],
    erk_wgt_stb53[2][0]+erk_wgt_stb53[2][1]+erk_wgt_stb53[2][2],
    erk_wgt_stb53[3][0]+erk_wgt_stb53[3][1]+erk_wgt_stb53[3][2]+erk_wgt_stb53[3][3],
    erk_wgt_stb53[4][0]+erk_wgt_stb53[4][1]+erk_wgt_stb53[4][2]+erk_wgt_stb53[4][3]+erk_wgt_stb53[4][4]
  };

  const double erk_c_stb53[5] = {
    0.0,
    1.0,
    2.0*ALPHA2,
    ALPHA3*(2.0*ALPHA2+1.0),
    ALPHA4*(2.0*ALPHA2*ALPHA3+ALPHA3+1.0),
  };

#undef ALPHA2
#undef ALPHA3
#undef ALPHA4
#undef ALPHA5
#undef ALPHA6
#undef ALPHA7
      
  // from stb:
  //const double c2 =   1.0000000000000000;
  //const double c3 =   0.1340039141171431;
  //const double c4 =   0.5107855836442267;
  //const double c5 =   0.7806160705207363;

  //extern int nerk;
  //extern double erk_wgt[5][5];
  //extern double erk_c[5];
  //void processErk(Param * param,const bool b_help);

};

namespace Ascher343 { 

  const double c[3] = {0.43586652150845899942, 0.71793326075422949975, 1.0};

  // Butcher Tab. DIRK

  const double A[3][3] = { {0.43586652150845899942, 0.0, 0.0},
                           {0.28206673924577050029, 0.43586652150845899942, 0.0}, 
                           {1.2084966491760100703, -0.6443631706844690698, 0.43586652150845899942} };

  const double bt[3] =     {1.2084966491760100703, -0.6443631706844690698, 0.43586652150845899942};

  // Butcher Tab. ERK

  const double A_hat[4][4] = { {0.0,0.0,0.0,0.0},
                               {0.43586652150845899942, 0.0, 0.0, 0.0},
                               {0.3212788862720422546, 0.3966543744821872450, 0.0, 0.0}, 
                               {-0.1058582958, 0.5529291479, 0.5529291479, 0.0} };

  const double bt_hat[4] = {0.0,1.2084966491760100703, -0.6443631706844690698, 0.43586652150845899942};

};



#endif
