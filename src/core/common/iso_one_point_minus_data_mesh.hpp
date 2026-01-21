
// one point minus...
{

  assert(PHIP0 >= 0.0);
  assert(PHIP1 >= 0.0);
  assert(PHIP2 >= 0.0);
  assert(PHIM < 0.0);
  double x30[3]; FOR_I3 x30[i] = (-PHIM*XP0[i]+PHIP0*XM[i])/(PHIP0-PHIM);
  double x31[3]; FOR_I3 x31[i] = (-PHIM*XP1[i]+PHIP1*XM[i])/(PHIP1-PHIM);
  double x32[3]; FOR_I3 x32[i] = (-PHIM*XP2[i]+PHIP2*XM[i])/(PHIP2-PHIM);

  const double d30 = (-PHIM*DP0+PHIP0*DM)/(PHIP0-PHIM);
  const double d31 = (-PHIM*DP1+PHIP1*DM)/(PHIP1-PHIM);
  const double d32 = (-PHIM*DP2+PHIP2*DM)/(PHIP2-PHIM);

  triVec.push_back(pair<SimpleTriWithData,int> (SimpleTriWithData(x30,x31,x32,d30,d31,d32),MESH_EDGES));
  return;

  //return( -SIGNED_TET_VOLUME_6(x30,x31,x32,XM)*PHIM +
  //	   SIGNED_TET_VOLUME_6(x32,x31,x30,XP2)*PHIP2 +
  //	   SIGNED_TET_VOLUME_6(x30,x31,XP1,XP2)*(PHIP1+PHIP2) +
  //	   SIGNED_TET_VOLUME_6(XP0,XP2,x30,XP1)*(PHIP0+PHIP1+PHIP2) );
  //  return( SIGNED_TET_VOLUME_6(x30,x31,x32,XM)*PHIM +
  //	  SIGNED_TET_VOLUME_6(x32,x31,x30,XP2)*PHIP2 +
  //	  SIGNED_TET_VOLUME_6(x30,x31,XP1,XP2)*(PHIP1+PHIP2) +
  //	  SIGNED_TET_VOLUME_6(XP0,XP2,x30,XP1)*(PHIP0+PHIP1+PHIP2) );

}
