
 case 90364472:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13941632:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 354451650:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 115873256:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373557750:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13976448:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 117999512:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13979544:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 271546746:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 115873580:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373444022:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13976460:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 90364148:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13941620:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 354451974:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 32968352:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373557762:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13862720:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 39346148:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13871972:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 288060158:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 348074340:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 14213130:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373548516:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 288060644:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 14213148:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 39345662:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 348074664:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13871954:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373548528:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 39345824:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13871960:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 288060482:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 99359844:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 14213142:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373207340:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 271547232:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373444040:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 117999026:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 297056016:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13979526:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373478856:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 354452136:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373557768:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 90363986:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 297056340:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13941614:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373478868:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 271546908:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373444028:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 117999350:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 269420976:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13979538:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373440944:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 115873742:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13976466:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 32968838:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13862738:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 32968514:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13862726:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 348074826:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373548534:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 99360330:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373207358:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 99360006:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373207346:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 297056502:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373478874:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 269421462:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373440962:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 269421138:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373440950:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 290338026:
   return 0;

 case 124755414:
   return 0;

 case 41623044:
   return 0;

 case 345797282:
   return 0;

 case 124414220:
   return 0;

 case 263006262:
   return 0;

 case 97120518:
   return 0;

 case 290299634:
   return 0;

 case 345835176:
   return 0;

 case 41509334:
   return 0;

 case 290640816:
   return 0;

 case 41509802:
   return 0;

 case 97196658:
   return 0;

 case 345569510:
   return 0;

 case 345910668:
   return 0;

 case 345797426:
   return 0;

 case 263005782:
   return 0;

 case 96892920:
   return 0;

 case 290527406:
   return 0;

 case 124490360:
   return 0;

 case 262930122:
   return 0;

 case 41585294:
   return 0;

 case 96779186:
   return 0;

 case 290565324:
   return 0;

 case 124528446:
   return 0;

 case 124528290:
   return 0;

 case 262892036:
   return 0;

 case 290299644:
   return 0;

 case 97120788:
   return 0;

 case 41850974:
   return 0;

 case 41509824:
   return 0;

 case 290223722:
   return 0;

 case 41509928:
   return 0;

 case 96779196:
   return 0;

 case 41585564:
   return 0;

 case 345911150:
   return 0;

 case 41585316:
   return 0;

 case 345910578:
   return 0;

 case 96779780:
   return 0;

 case 290527416:
   return 0;

 case 262930392:
   return 0;

 case 41623202:
   return 0;

 case 345797448:
   return 0;

 case 124414118:
   return 0;

 case 263006376:
   return 0;

 case 345569496:
   return 0;

 case 97196280:
   return 0;

 case 96855486:
   return 0;

 case 97196672:
   return 0;

 case 96855122:
   return 0;

 case 345569888:
   return 0;

 case 41509320:
   return 0;

 case 41509424:
   return 0;

 case 290224302:
   return 0;

 case 97120532:
   return 0;

 case 41850618:
   return 0;

 case 290300012:
   return 0;

 case 345797268:
   return 0;

 case 263005884:
   return 0;

 case 97082610:
   return 0;

 case 290338040:
   return 0;

 case 262664702:
   return 0;

 case 124755792:
   return 0;

 case 124755392:
   return 0;

 case 290337432:
   return 0;

 case 262665090:
   return 0;

 case 263006252:
   return 0;

 case 97082894:
   return 0;

 case 345797012:
   return 0;

 case 263005760:
   return 0;

 case 345796832:
   return 0;

 case 124414722:
   return 0;

 case 262930112:
   return 0;

 case 41623494:
   return 0;

 case 290527136:
   return 0;

 case 124528268:
   return 0;

 case 124527852:
   return 0;

 case 262892214:
   return 0;

 case 124528436:
   return 0;

 case 262892474:
   return 0;

 case 124528020:
   return 0;

 case 97082624:
   return 0;

 case 262665080:
   return 0;

 case 96855500:
   return 0;

 case 290224316:
   return 0;

 case 41850996:
   return 0;

 case 262892204:
   return 0;

 case 124414712:
   return 0;

 case 41623224:
   return 0;

 case 345911172:
   return 0;

 case 290564988:
   return 0;

 case 41850978:
   return 0;

 case 41509820:
   return 0;

 case 290223830:
   return 0;

 case 97196172:
   return 0;

 case 345569492:
   return 0;

 case 345911154:
   return 0;

 case 41585312:
   return 0;

 case 345910686:
   return 0;

 case 96779672:
   return 0;

 case 290337864:
   return 0;

 case 124755408:
   return 0;

 case 41623206:
   return 0;

 case 345797444:
   return 0;

 case 124414226:
   return 0;

 case 263006268:
   return 0;

 case 97120854:
   return 0;

 case 290299970:
   return 0;

 case 41509316:
   return 0;

 case 290641302:
   return 0;

 case 96855164:
   return 0;

 case 345835194:
   return 0;

 case 345797264:
   return 0;

 case 263005776:
   return 0;

 case 96893082:
   return 0;

 case 290527568:
   return 0;

 case 124490366:
   return 0;

 case 262930128:
   return 0;

 case 262665074:
   return 0;

 case 97082462:
   return 0;

 case 124414706:
   return 0;

 case 41623062:
   return 0;

 case 124528284:
   return 0;

 case 262892198:
   return 0;

 case 124528452:
   return 0;

 case 262892042:
   return 0;

 case 100427262:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 376397450:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 344767920:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 42573674:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 287450712:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 42574142:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 100424670:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 374271698:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 344846328:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 376785762:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 383153838:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 12912108:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 374508218:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 114808916:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 272611566:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 42652550:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 99969290:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 287372304:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 374508362:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 272611086:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 10786356:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 376633970:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 114806324:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 272614158:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 383163066:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 383162910:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 4257416:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 286996304:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13148808:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 10643960:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 10634888:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 4266656:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 114809408:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 12912288:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 4257584:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 344846832:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 11023038:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 286993226:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 374271680:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 344846346:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 100424184:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 344846814:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 287451198:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 344767938:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 42573656:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 10786518:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 114806330:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 374508200:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 12912126:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 272611080:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 114809402:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13148790:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 286995818:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 376776528:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 12912270:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 114808922:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 376785600:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 10634726:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 383153832:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 4266650:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 4257578:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 4257422:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 383162904:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 100048184:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 99969776:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 42652568:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 376634132:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 272614164:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 42574160:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 374508380:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 272611572:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 383163072:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 41623058:
   return 0;

 case 124414598:
   return 0;

 case 345835190:
   return 0;

 case 290641194:
   return 0;

 case 345910682:
   return 0;

 case 345911046:
   return 0;

 case 96892910:
   return 0;

 case 124490090:
   return 0;

 case 290565314:
   return 0;

 case 290565054:
   return 0;

 case 290640806:
   return 0;

 case 345834906:
   return 0;

 case 262892058:
   return 0;

 case 262892630:
   return 0;

 case 124490382:
   return 0;

 case 96893514:
   return 0;

 case 124414242:
   return 0;

 case 41623638:
   return 0;

 case 13858854:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 29780318:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373561640:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373194234:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 357640332:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 89794338:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373439670:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 268358694:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13980824:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373439502:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 119061956:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 268358526:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13940178:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 89301698:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373480316:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13858362:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 298118952:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 29780138:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13858838:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 29779886:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373561656:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373194218:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 357640764:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 89793906:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373439654:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 268358262:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13980840:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373439486:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 119062388:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 268358094:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13940162:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 89301266:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373480332:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13858346:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 298119384:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 29779706:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13858850:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 29780210:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373561652:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373194222:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 357640656:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 89794014:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373439666:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 268358586:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13980836:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373439490:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 119062280:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 268358202:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13940174:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 89301590:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373480328:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13858350:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 298119276:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 29779814:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 1068788:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 4253528:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 384074318:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 3346332:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 272624694:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 114795800:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 28712900:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 10668744:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 300932894:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 86566488:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 100324070:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 344492052:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 28715492:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 12794496:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 300854486:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 3582852:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 42928454:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 286982696:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 9723480:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 114805028:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 377697014:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 9723636:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 272615622:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 114804872:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 258674496:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 287333100:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 128745998:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 258753384:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 100087550:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 344728572:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 258752904:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 344728716:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 128667590:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 9960156:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 42691934:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 286991768:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 3194540:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 4256120:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 386351718:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 1068932:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 383167446:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 4253048:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 86108516:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 10747152:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 358707606:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 28715960:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 376752230:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 12794028:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 86111108:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 12872904:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 358705014:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 1078004:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 374626478:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 10630352:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 86566002:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 344492034:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 358704996:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 1078022:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 374625992:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 10630838:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 3346170:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 114795794:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 386351700:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 1068950:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 383166960:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 4253534:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 86487594:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 287096418:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 358707588:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 28715978:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 376751744:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 12794514:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 28715474:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 12794010:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 301309380:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 3203774:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 374547584:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 10633430:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 1068770:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 4253042:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 384225948:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 3194702:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 383164368:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 4256126:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 28712882:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 10668258:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 301311972:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 86111594:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 376673336:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 12872922:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 258752898:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 344728554:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 128667584:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 9960162:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 42691772:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 286991930:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 9723474:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 114804866:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 377697008:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 9723642:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 272615460:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 114805034:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 258674490:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 287332938:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 128745992:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 258753390:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 100087388:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 344728734:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373562126:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 357640350:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13858368:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373480310:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 29780300:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 298118790:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13980986:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 119061962:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373439508:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13980818:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 268358688:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 119061794:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 14226254:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 297626150:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373194240:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373561634:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 89794500:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 357640170:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373562138:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 357640674:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13858364:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373480314:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 29780192:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 298118898:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13980998:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 119062286:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373439504:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13980822:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 268358580:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 119061902:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 14226266:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 297626474:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373194236:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373561638:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 89794392:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 357640278:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373562142:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 357640782:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13858352:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373480326:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 29779868:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 298119222:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13981002:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 119062394:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373439492:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13980834:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 268358256:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 119062226:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 14226270:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 297626582:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373194224:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373561650:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 89794068:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 357640602:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 345797430:
   return 0;

 case 263005890:
   return 0;

 case 97082448:
   return 0;

 case 290337878:
   return 0;

 case 262664696:
   return 0;

 case 124755786:
   return 0;

 case 41585298:
   return 0;

 case 96779294:
   return 0;

 case 290299956:
   return 0;

 case 97196186:
   return 0;

 case 97120476:
   return 0;

 case 345569870:
   return 0;

 case 41509806:
   return 0;

 case 41509442:
   return 0;

 case 290223816:
   return 0;

 case 290565002:
   return 0;

 case 41850600:
   return 0;

 case 290565366:
   return 0;

 case 124528430:
   return 0;

 case 124527858:
   return 0;

 case 262892052:
   return 0;

 case 124528274:
   return 0;

 case 262892468:
   return 0;

 case 124528014:
   return 0;

 case 262930106:
   return 0;

 case 290526974:
   return 0;

 case 124490376:
   return 0;

 case 263005766:
   return 0;

 case 96893352:
   return 0;

 case 345796994:
   return 0;

 case 263006246:
   return 0;

 case 345796850:
   return 0;

 case 124414236:
   return 0;

 case 124755398:
   return 0;

 case 41623476:
   return 0;

 case 290337594:
   return 0;

 case 290527578:
   return 0;

 case 262930398:
   return 0;

 case 41623040:
   return 0;

 case 345797286:
   return 0;

 case 124414112:
   return 0;

 case 263006370:
   return 0;

 case 96855174:
   return 0;

 case 96855434:
   return 0;

 case 345835172:
   return 0;

 case 41509338:
   return 0;

 case 290640708:
   return 0;

 case 41509910:
   return 0;

 case 96779682:
   return 0;

 case 41585582:
   return 0;

 case 345910664:
   return 0;

 case 345569514:
   return 0;

 case 345910560:
   return 0;

 case 97196766:
   return 0;

 case 373562136:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 357640620:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13980996:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 119062232:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 14226264:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 297626420:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373562124:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 357640296:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13980984:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 119061908:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 14226252:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 297626096:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373562120:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 357640188:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13980980:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 119061800:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 14226248:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 297625988:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 355476984:
   return 0;

 case 345923804:
   return 0;

 case 87630342:
   return 0;

 case 355135448:
   return 0;

 case 345556374:
   return 0;

 case 97209780:
   return 0;

 case 299903964:
   return 0;

 case 124768536:
   return 0;

 case 32057106:
   return 0;

 case 355363220:
   return 0;

 case 124401098:
   return 0;

 case 263019384:
   return 0;

 case 355400844:
   return 0;

 case 290653928:
   return 0;

 case 87554850:
   return 0;

 case 299789660:
   return 0;

 case 290286522:
   return 0;

 case 41864096:
   return 0;

 case 127603244:
   return 0;

 case 41627432:
   return 0;

 case 261943002:
   return 0;

 case 291362516:
   return 0;

 case 345795810:
   return 0;

 case 97121976:
   return 0;

 case 348986072:
   return 0;

 case 263010156:
   return 0;

 case 95830038:
   return 0;

 case 291590288:
   return 0;

 case 124488902:
   return 0;

 case 262931580:
   return 0;

 case 127678736:
   return 0;

 case 96897284:
   return 0;

 case 261867510:
   return 0;

 case 125477000:
   return 0;

 case 290525958:
   return 0;

 case 41624660:
   return 0;

 case 291704076:
   return 0;

 case 345836648:
   return 0;

 case 38320778:
   return 0;

 case 127944060:
   return 0;

 case 41504946:
   return 0;

 case 290342400:
   return 0;

 case 125591328:
   return 0;

 case 124529748:
   return 0;

 case 259703390:
   return 0;

 case 127716936:
   return 0;

 case 262887662:
   return 0;

 case 124532820:
   return 0;

 case 291627936:
   return 0;

 case 290566772:
   return 0;

 case 38396918:
   return 0;

 case 349099224:
   return 0;

 case 96774822:
   return 0;

 case 345915524:
   return 0;

 case 383167440:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 386351556:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 272615616:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 377696852:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 272624688:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 384074156:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 374626460:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 358704528:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 42691916:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 128667104:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 42928436:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 300854000:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 376790136:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 386342484:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 100428720:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 377460332:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 100437792:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 383837636:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 259817730:
   return 0;

 case 345793074:
   return 0;

 case 38662350:
   return 0;

 case 290219942:
   return 0;

 case 38434578:
   return 0;

 case 124410338:
   return 0;

 case 31943990:
   return 0;

 case 41496702:
   return 0;

 case 87289562:
   return 0;

 case 96842378:
   return 0;

 case 87516686:
   return 0;

 case 262651958:
   return 0;

 case 95716898:
   return 0;

 case 41583858:
   return 0;

 case 261602198:
   return 0;

 case 97081166:
   return 0;

 case 261829322:
   return 0;

 case 262890746:
   return 0;

 case 95792552:
   return 0;

 case 96853716:
   return 0;

 case 349023570:
   return 0;

 case 38321264:
   return 0;

 case 290645666:
   return 0;

 case 41504964:
   return 0;

 case 95716412:
   return 0;

 case 41583840:
   return 0;

 case 349099710:
   return 0;

 case 259476428:
   return 0;

 case 345915542:
   return 0;

 case 97078088:
   return 0;

 case 261829160:
   return 0;

 case 262890740:
   return 0;

 case 127717098:
   return 0;

 case 259703552:
   return 0;

 case 124532826:
   return 0;

 case 262887668:
   return 0;

 case 32019644:
   return 0;

 case 96766560:
   return 0;

 case 299865638:
   return 0;

 case 87630828:
   return 0;

 case 97133966:
   return 0;

 case 345556392:
   return 0;

 case 31943504:
   return 0;

 case 41496684:
   return 0;

 case 299790146:
   return 0;

 case 32285040:
   return 0;

 case 41864114:
   return 0;

 case 290210708:
   return 0;

 case 87516524:
   return 0;

 case 262651952:
   return 0;

 case 355363382:
   return 0;

 case 32057268:
   return 0;

 case 263019390:
   return 0;

 case 124401104:
   return 0;

 case 259741752:
   return 0;

 case 290523204:
   return 0;

 case 125552978:
   return 0;

 case 261943488:
   return 0;

 case 96894530:
   return 0;

 case 345795828:
   return 0;

 case 259817244:
   return 0;

 case 345793056:
   return 0;

 case 125477486:
   return 0;

 case 96057972:
   return 0;

 case 41624678:
   return 0;

 case 290298512:
   return 0;

 case 38434416:
   return 0;

 case 124410332:
   return 0;

 case 291590450:
   return 0;

 case 95830200:
   return 0;

 case 262931586:
   return 0;

 case 124488908:
   return 0;

 case 377460326:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 100428558:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 377696846:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 272615454:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 128667098:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 42691754:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 386342466:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 376789650:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 386351538:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 383166954:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 358704510:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 374625974:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 384216714:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 376787058:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 384225786:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 383164362:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 301308894:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 374547566:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 348985910:
   return 0;

 case 263010150:
   return 0;

 case 127602758:
   return 0;

 case 41627414:
   return 0;

 case 348758138:
   return 0;

 case 97200546:
   return 0;

 case 299903802:
   return 0;

 case 124768530:
   return 0;

 case 355476498:
   return 0;

 case 345923786:
   return 0;

 case 300130926:
   return 0;

 case 290578110:
   return 0;

 case 125591166:
   return 0;

 case 124529742:
   return 0;

 case 291703590:
   return 0;

 case 345836630:
   return 0;

 case 125818290:
   return 0;

 case 290339322:
   return 0;

 case 124414604:
   return 0;

 case 41623220:
   return 0;

 case 124490096:
   return 0;

 case 96893072:
   return 0;

 case 262892636:
   return 0;

 case 262892220:
   return 0;

 case 345911064:
   return 0;

 case 345911168:
   return 0;

 case 345834924:
   return 0;

 case 290641292:
   return 0;

 case 41623656:
   return 0;

 case 124414728:
   return 0;

 case 290224208:
   return 0;

 case 41850992:
   return 0;

 case 290299700:
   return 0;

 case 97120844:
   return 0;

 case 97083056:
   return 0;

 case 262665096:
   return 0;

 case 89301104:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13940156:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 89301536:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13940172:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 89301212:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13940160:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 268358100:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373439648:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 268358532:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373439664:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 268358208:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373439652:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 29779724:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13858832:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 29780156:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13858848:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 29779832:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13858836:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 349025276:
   return 0;

 case 291707964:
   return 0;

 case 95790834:
   return 0;

 case 291708548:
   return 0;

 case 95791094:
   return 0;

 case 349025028:
   return 0;

 case 349112432:
   return 0;

 case 355480872:
   return 0;

 case 95703678:
   return 0;

 case 125831780:
   return 0;

 case 32018186:
   return 0;

 case 299905004:
   return 0;

 case 127720856:
   return 0;

 case 261825390:
   return 0;

 case 125592344:
   return 0;

 case 259702386:
   return 0;

 case 125592624:
   return 0;

 case 259702094:
   return 0;

 case 127721472:
   return 0;

 case 261824786:
   return 0;

 case 125831412:
   return 0;

 case 299905260:
   return 0;

 case 259463306:
   return 0;

 case 349112832:
   return 0;

 case 87512150:
   return 0;

 case 355481480:
   return 0;

 case 38316890:
   return 0;

 case 38316318:
   return 0;

 case 32017902:
   return 0;

 case 95703290:
   return 0;

 case 300144048:
   return 0;

 case 31930746:
   return 0;

 case 300144440:
   return 0;

 case 31930382:
   return 0;

 case 87512730:
   return 0;

 case 259462950:
   return 0;

 case 345924186:
   return 0;

 case 355477106:
   return 0;

 case 96766172:
   return 0;

 case 97133598:
   return 0;

 case 32019360:
   return 0;

 case 299865894:
   return 0;

 case 41628030:
   return 0;

 case 127603374:
   return 0;

 case 290522600:
   return 0;

 case 96894810:
   return 0;

 case 259741460:
   return 0;

 case 125553258:
   return 0;

 case 345836382:
   return 0;

 case 291704174:
   return 0;

 case 96853976:
   return 0;

 case 290645082:
   return 0;

 case 95792292:
   return 0;

 case 349023818:
   return 0;

 case 124768898:
   return 0;

 case 299903546:
   return 0;

 case 262651596:
   return 0;

 case 263018990:
   return 0;

 case 87517104:
   return 0;

 case 355362774:
   return 0;

 case 263010734:
   return 0;

 case 348985662:
   return 0;

 case 124409760:
   return 0;

 case 262931834:
   return 0;

 case 38434988:
   return 0;

 case 291589866:
   return 0;

 case 124529462:
   return 0;

 case 125590886:
   return 0;

 case 262891032:
   return 0;

 case 124532210:
   return 0;

 case 261829764:
   return 0;

 case 127716482:
   return 0;

 case 290578502:
   return 0;

 case 300131318:
   return 0;

 case 41496320:
   return 0;

 case 41863722:
   return 0;

 case 31943868:
   return 0;

 case 299789754:
   return 0;

 case 97201154:
   return 0;

 case 348758538:
   return 0;

 case 345792476:
   return 0;

 case 41624934:
   return 0;

 case 259817600:
   return 0;

 case 125477118:
   return 0;

 case 290339066:
   return 0;

 case 125818658:
   return 0;

 case 41584124:
   return 0;

 case 345914934:
   return 0;

 case 95716800:
   return 0;

 case 349099310:
   return 0;

 case 290337446:
   return 0;

 case 124755770:
   return 0;

 case 41623652:
   return 0;

 case 345796998:
   return 0;

 case 124414620:
   return 0;

 case 263005874:
   return 0;

 case 97120802:
   return 0;

 case 290300022:
   return 0;

 case 345834920:
   return 0;

 case 41509914:
   return 0;

 case 290641184:
   return 0;

 case 41509446:
   return 0;

 case 97196294:
   return 0;

 case 345569874:
   return 0;

 case 345911060:
   return 0;

 case 345796854:
   return 0;

 case 263006354:
   return 0;

 case 96893504:
   return 0;

 case 290527146:
   return 0;

 case 124490112:
   return 0;

 case 262930382:
   return 0;

 case 41585586:
   return 0;

 case 96779790:
   return 0;

 case 290565044:
   return 0;

 case 124527842:
   return 0;

 case 124527998:
   return 0;

 case 262892652:
   return 0;

 case 262664712:
   return 0;

 case 97082880:
   return 0;

 case 124414128:
   return 0;

 case 41623472:
   return 0;

 case 262892484:
   return 0;

 case 41850596:
   return 0;

 case 290223708:
   return 0;

 case 345910556:
   return 0;

 case 96855108:
   return 0;

 case 262888272:
   return 0;

 case 259703844:
   return 0;

 case 124488648:
   return 0;

 case 95830460:
   return 0;

 case 124401492:
   return 0;

 case 32057552:
   return 0;

 case 97078668:
   return 0;

 case 259476072:
   return 0;

 case 290298228:
   return 0;

 case 96057584:
   return 0;

 case 290211072:
   return 0;

 case 32284676:
   return 0;

 case 41505536:
   return 0;

 case 38320692:
   return 0;

 case 345795536:
   return 0;

 case 261942884:
   return 0;

 case 345556748:
   return 0;

 case 87630248:
   return 0;

 case 124490106:
   return 0;

 case 96893342:
   return 0;

 case 262930376:
   return 0;

 case 124414134:
   return 0;

 case 290526984:
   return 0;

 case 41623634:
   return 0;

 case 124414614:
   return 0;

 case 41623490:
   return 0;

 case 263005868:
   return 0;

 case 262664718:
   return 0;

 case 345796836:
   return 0;

 case 97083042:
   return 0;

 case 262892646:
   return 0;

 case 262892490:
   return 0;

 case 124527836:
   return 0;

 case 290299686:
   return 0;

 case 97120466:
   return 0;

 case 41585568:
   return 0;

 case 345910574:
   return 0;

 case 96779304:
   return 0;

 case 345911042:
   return 0;

 case 290224194:
   return 0;

 case 41850614:
   return 0;

 case 41509428:
   return 0;

 case 345834902:
   return 0;

 case 290640698:
   return 0;

 case 96855444:
   return 0;

 case 28711424:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 9605376:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 301313430:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 86107220:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 377736218:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 9684276:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 28702352:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 3228072:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 301322502:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 3190652:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 384113522:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 1067492:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 1064396:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 384230322:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 3193244:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 386353014:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 258635124:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 128785358:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 258635292:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 128785202:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 258398604:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 86448228:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 129021878:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 9605868:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 300972098:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 28711604:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 377815106:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 358708902:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 358709046:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 377814626:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 86211708:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 358718118:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 1064900:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 384191930:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 386356074:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 386355606:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 358708884:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 377814620:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 86448390:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 300972260:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 258398610:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 129021884:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 386352996:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 384229836:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 3306966:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 384192416:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 86097986:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 358718136:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 386355588:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 3228558:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 301208780:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 28702370:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 128785196:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 258635286:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 128785364:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 258635130:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 377736212:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 301313268:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 9684270:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 377815112:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 86107058:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 358709064:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 9605862:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 28711442:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 28711586:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 9605382:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 384227244:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 1067474:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 386356092:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 3190166:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 1064882:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 1064414:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 124528004:
   return 0;

 case 124755776:
   return 0;

 case 290337608:
   return 0;

 case 263006360:
   return 0;

 case 345797016:
   return 0;

 case 345569892:
   return 0;

 case 97196780:
   return 0;

 case 41509932:
   return 0;

 case 290565380:
   return 0;

 case 38395212:
   return 0;

 case 95712524:
   return 0;

 case 291629654:
   return 0;

 case 95711940:
   return 0;

 case 291629394:
   return 0;

 case 38395460:
   return 0;

 case 38308056:
   return 0;

 case 31939616:
   return 0;

 case 291716810:
   return 0;

 case 261588708:
   return 0;

 case 355402302:
   return 0;

 case 87515484:
   return 0;

 case 259699632:
   return 0;

 case 125595098:
   return 0;

 case 261828144:
   return 0;

 case 127718102:
   return 0;

 case 261827864:
   return 0;

 case 127718394:
   return 0;

 case 259699016:
   return 0;

 case 125595702:
   return 0;

 case 261589076:
   return 0;

 case 87515228:
   return 0;

 case 127957182:
   return 0;

 case 38307656:
   return 0;

 case 299908338:
   return 0;

 case 31939008:
   return 0;

 case 349103598:
   return 0;

 case 349104170:
   return 0;

 case 355402586:
   return 0;

 case 291717198:
   return 0;

 case 87276440:
   return 0;

 case 355489742:
   return 0;

 case 87276048:
   return 0;

 case 355490106:
   return 0;

 case 299907758:
   return 0;

 case 127957538:
   return 0;

 case 124532216:
   return 0;

 case 127716644:
   return 0;

 case 262888278:
   return 0;

 case 124529456:
   return 0;

 case 259704006:
   return 0;

 case 125590724:
   return 0;

 case 262931840:
   return 0;

 case 291590028:
   return 0;

 case 124488654:
   return 0;

 case 263010728:
   return 0;

 case 95830622:
   return 0;

 case 348985500:
   return 0;

 case 263018996:
   return 0;

 case 355362936:
   return 0;

 case 124401498:
   return 0;

 case 124768892:
   return 0;

 case 32057714:
   return 0;

 case 299903384:
   return 0;

 case 290341820:
   return 0;

 case 127944416:
   return 0;

 case 41505554:
   return 0;

 case 345836364:
   return 0;

 case 38321178:
   return 0;

 case 291703688:
   return 0;

 case 97122260:
   return 0;

 case 291362904:
   return 0;

 case 345795554:
   return 0;

 case 41628012:
   return 0;

 case 261943370:
   return 0;

 case 127602888:
   return 0;

 case 97209416:
   return 0;

 case 355135812:
   return 0;

 case 345556766:
   return 0;

 case 345924168:
   return 0;

 case 87630734:
   return 0;

 case 355476620:
   return 0;

 case 345914952:
   return 0;

 case 349099796:
   return 0;

 case 96775406:
   return 0;

 case 290566512:
   return 0;

 case 38396670:
   return 0;

 case 291628196:
   return 0;

 case 41624952:
   return 0;

 case 125477604:
   return 0;

 case 290525678:
   return 0;

 case 96897888:
   return 0;

 case 261867230:
   return 0;

 case 127679028:
   return 0;

 case 41863740:
   return 0;

 case 299790240:
   return 0;

 case 290286890:
   return 0;

 case 290654316:
   return 0;

 case 87554594:
   return 0;

 case 355401128:
   return 0;

 case 290219334:
   return 0;

 case 38661950:
   return 0;

 case 97081422:
   return 0;

 case 261601830:
   return 0;

 case 96841986:
   return 0;

 case 87289170:
   return 0;

 case 124409754:
   return 0;

 case 38434826:
   return 0;

 case 262891026:
   return 0;

 case 261829602:
   return 0;

 case 262651590:
   return 0;

 case 87516942:
   return 0;

 case 345792458:
   return 0;

 case 259817114:
   return 0;

 case 41584106:
   return 0;

 case 95716314:
   return 0;

 case 41496302:
   return 0;

 case 31943382:
   return 0;
