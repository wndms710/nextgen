// three point minus...
{
  
  assert(PHIP >= 0.0);
  assert(PHIM0 < 0.0);
  assert(PHIM1 < 0.0);
  assert(PHIM2 < 0.0);

  double x01[3]; FOR_I3 x01[i] = (-PHIM0*XP[i]+PHIP*XM0[i])/(PHIP-PHIM0);
  double x02[3]; FOR_I3 x02[i] = (-PHIM1*XP[i]+PHIP*XM1[i])/(PHIP-PHIM1);
  double x03[3]; FOR_I3 x03[i] = (-PHIM2*XP[i]+PHIP*XM2[i])/(PHIP-PHIM2);

  const double d01 = (-PHIM0*DP+PHIP*DM0)/(PHIP-PHIM0);
  const double d02 = (-PHIM1*DP+PHIP*DM1)/(PHIP-PHIM1);
  const double d03 = (-PHIM2*DP+PHIP*DM2)/(PHIP-PHIM2);

  triVec.push_back(SimpleTriWithData(x01,x02,x03,d01,d02,d03));
  return;

  //  return( SIGNED_TET_VOLUME_6(XP,x01,x02,x03)*PHIP -
  //	  SIGNED_TET_VOLUME_6(x01,x02,x03,XM1)*PHIM1 - 
  //	  SIGNED_TET_VOLUME_6(x01,x03,XM0,XM1)*(PHIM0+PHIM1) -
  //	  SIGNED_TET_VOLUME_6(x03,XM0,XM1,XM2)*(PHIM0+PHIM1+PHIM2) );
//  return( SIGNED_TET_VOLUME_6(XP,x01,x02,x03)*PHIP +
//	  SIGNED_TET_VOLUME_6(x01,x02,x03,XM1)*PHIM1 + 
//	  SIGNED_TET_VOLUME_6(x01,x03,XM0,XM1)*(PHIM0+PHIM1) +
//	  SIGNED_TET_VOLUME_6(x03,XM0,XM1,XM2)*(PHIM0+PHIM1+PHIM2) );
}


