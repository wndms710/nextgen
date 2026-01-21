// two point minus...

{

  assert(PHIP0 >= 0.0);
  assert(PHIP1 >= 0.0);
  assert(PHIM0 < 0.0);
  assert(PHIM1 < 0.0);

  double x02[3]; FOR_I3 x02[i] = (-PHIM0*XP0[i]+PHIP0*XM0[i])/(PHIP0-PHIM0);
  double x12[3]; FOR_I3 x12[i] = (-PHIM1*XP0[i]+PHIP0*XM1[i])/(PHIP0-PHIM1);
  double x03[3]; FOR_I3 x03[i] = (-PHIM0*XP1[i]+PHIP1*XM0[i])/(PHIP1-PHIM0);
  double x13[3]; FOR_I3 x13[i] = (-PHIM1*XP1[i]+PHIP1*XM1[i])/(PHIP1-PHIM1);
  
  triVec.push_back(SimpleTriWithWgts(x12,x02,x13,
                                     IP0,IM1,-PHIM1/(PHIP0-PHIM1),
                                     IP0,IM0,-PHIM0/(PHIP0-PHIM0),
                                     IP1,IM1,-PHIM1/(PHIP1-PHIM1)));
  
  triVec.push_back(SimpleTriWithWgts(x03,x13,x02,
                                     IP1,IM0,-PHIM0/(PHIP1-PHIM0),
                                     IP1,IM1,-PHIM1/(PHIP1-PHIM1),
                                     IP0,IM0,-PHIM0/(PHIP0-PHIM0)));
  
  return;

  //return( SIGNED_TET_VOLUME_6(XP0,XP1,x02,x03)*(PHIP0+PHIP1) + 
  //	  SIGNED_TET_VOLUME_6(x02,XP1,x12,x03)*PHIP1 + 
  //	  SIGNED_TET_VOLUME_6(x03,x13,XP1,x12)*PHIP1 -
  //	  SIGNED_TET_VOLUME_6(x02,x12,XM0,XM1)*(PHIM0+PHIM1) -
  //	  SIGNED_TET_VOLUME_6(x12,x02,x13,XM1)*PHIM1 -
  //	  SIGNED_TET_VOLUME_6(x03,x13,x02,XM1)*PHIM1 );
//  return( SIGNED_TET_VOLUME_6(XP0,XP1,x02,x03)*(PHIP0+PHIP1) + 
//	  SIGNED_TET_VOLUME_6(x02,XP1,x12,x03)*PHIP1 + 
//	  SIGNED_TET_VOLUME_6(x03,x13,XP1,x12)*PHIP1 +
//	  SIGNED_TET_VOLUME_6(x02,x12,XM0,XM1)*(PHIM0+PHIM1) +
//	  SIGNED_TET_VOLUME_6(x12,x02,x13,XM1)*PHIM1 +
//	  SIGNED_TET_VOLUME_6(x03,x13,x02,XM1)*PHIM1 );
}
