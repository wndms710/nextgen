
 case 359202300:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[3] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   return 2;

 case 349649484:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 56916008:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 330425586:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[3] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   return 2;

 case 37810532:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[3] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   return 2;

 case 292214322:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13445916:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[3] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   return 2;

 case 116404460:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 375037460:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 12381408:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[3] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   return 2;

 case 271017648:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[3] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   return 2;

 case 115339952:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13680978:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[3] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   return 2;

 case 287528474:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 374802398:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 330499620:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[3] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   return 2;

 case 99893634:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[3] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   return 2;

 case 346421292:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 244126718:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[3] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   return 2;

 case 142229598:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 143293776:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 244123634:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[3] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   return 2;

 case 245191052:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[3] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   return 2;

 case 140103666:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 9099602:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[3] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   return 2;

 case 47310866:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 378320892:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 9060068:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[3] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   return 2;

 case 340109784:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[3] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   return 2;

 case 18612884:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 9069308:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[3] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   return 2;

 case 24990356:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 378351186:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 244008128:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[3] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   return 2;

 case 362430294:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[3] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   return 2;

 case 56135976:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 330220172:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[3] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   return 2;

 case 142348020:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 28218206:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 359123388:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[3] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   return 2;

 case 37771490:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[3] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   return 2;

 case 292253364:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 12288260:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[3] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   return 2;

 case 47315564:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373974590:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 13444278:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 2; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[3] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   return 2;

 case 271016514:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[3] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   return 2;

 case 115341086:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 2; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 12257966:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 1; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[3] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   return 2;

 case 24995054:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 1; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 373739528:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 359197422:
   // intersection 0...
   idata[0] = TRI_EDGE_INT;
   idata[1] = 0; // edge index on tri1
   idata[2] = 0; // unused
   ddata[0] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[3] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   return 2;

 case 99892500:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[3] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   return 2;

 case 346460334:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   // intersection 1...
   idata[3] = EDGE_TRI_INT;
   idata[4] = 0; // edge index on tri0
   idata[5] = 0; // unused
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 37067645:
   return 0;

 case 290217633:
   return 0;

 case 295348349:
   return 0;

 case 161084031:
   return 0;

 case 290571941:
   return 0;

 case 41673467:
   return 0;

 case 216770411:
   return 0;

 case 345733895:
   return 0;

 case 226336943:
   return 0;

 case 161159523:
   return 0;

 case 345747039:
   return 0;

 case 96943319:
   return 0;

 case 36840521:
   return 0;

 case 124408053:
   return 0;

 case 295120577:
   return 0;

 case 92300073:
   return 0;

 case 124762337:
   return 0;

 case 262658157:
   return 0;

 case 261070395:
   return 0;

 case 97080855:
   return 0;

 case 125287217:
   return 0;

 case 110065707:
   return 0;

 case 290338337:
   return 0;

 case 41603807:
   return 0;

 case 81367413:
   return 0;

 case 41564585:
   return 0;

 case 277355267:
   return 0;

 case 110141199:
   return 0;

 case 345816699:
   return 0;

 case 96873659:
   return 0;

 case 261298167:
   return 0;

 case 262890459:
   return 0;

 case 125059445:
   return 0;

 case 262361205:
   return 0;

 case 124528733:
   return 0;

 case 262891761:
   return 0;

 case 96589025:
   return 0;

 case 290298957:
   return 0;

 case 347164215:
   return 0;

 case 302864321:
   return 0;

 case 97198967:
   return 0;

 case 345851525:
   return 0;

 case 276291791:
   return 0;

 case 345815219:
   return 0;

 case 84556653:
   return 0;

 case 302788181:
   return 0;

 case 41568981:
   return 0;

 case 290581649:
   return 0;

 case 96361901:
   return 0;

 case 124489377:
   return 0;

 case 347391339:
   return 0;

 case 40029311:
   return 0;

 case 263008547:
   return 0;

 case 124411947:
   return 0;

 case 98810178:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359920586:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 345378984:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 42118940:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 345378516:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 98374724:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 333561732:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 253954620:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27267738:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 360152912:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 119611574:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 267808920:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 333670758:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27499884:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 42041990:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 288609824:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 42041486:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 289046250:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 98451674:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 345301566:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 360152732:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 267808428:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 53839800:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 333580850:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 119648348:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 267772146:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 54071946:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 288646598:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 119612054:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27267594:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 253845594:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 133574888:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 253845762:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 133574732:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 133465862:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 53858594:
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 124490117:
   return 0;

 case 96892911:
   return 0;

 case 262930365:
   return 0;

 case 124414739:
   return 0;

 case 290527415:
   return 0;

 case 41623225:
   return 0;

 case 124414139:
   return 0;

 case 41623041:
   return 0;

 case 263006343:
   return 0;

 case 262665161:
   return 0;

 case 345797285:
   return 0;

 case 97082627:
   return 0;

 case 262892009:
   return 0;

 case 262892035:
   return 0;

 case 124528473:
   return 0;

 case 262892447:
   return 0;

 case 124528291:
   return 0;

 case 262892213:
   return 0;

 case 290565081:
   return 0;

 case 290565315:
   return 0;

 case 41585537:
   return 0;

 case 345911199:
   return 0;

 case 96779195:
   return 0;

 case 345911173:
   return 0;

 case 290640735:
   return 0;

 case 345835173:
   return 0;

 case 41509883:
   return 0;

 case 41851077:
   return 0;

 case 41509337:
   return 0;

 case 290224319:
   return 0;

 case 96892893:
   return 0;

 case 124490359:
   return 0;

 case 345797453:
   return 0;

 case 41623467:
   return 0;

 case 263005783:
   return 0;

 case 124414721:
   return 0;

 case 345834933:
   return 0;

 case 290640807:
   return 0;

 case 97120761:
   return 0;

 case 290224343:
   return 0;

 case 290299643:
   return 0;

 case 41850997:
   return 0;

 case 345910587:
   return 0;

 case 345910665:
   return 0;

 case 97196739:
   return 0;

 case 96855581:
   return 0;

 case 345569513:
   return 0;

 case 96855503:
   return 0;

 case 41623017:
   return 0;

 case 124414219:
   return 0;

 case 290338053:
   return 0;

 case 97082867:
   return 0;

 case 124755415:
   return 0;

 case 262665089:
   return 0;

 case 3361241615:
   return 0;

 case 2622300633:
   return 0;

 case 2804144291:
   return 0;

 case 3359197661:
   return 0;

 case 3543085273:
   return 0;

 case 1132140075:
   return 0;

 case 3359197175:
   return 0;

 case 1132140057:
   return 0;

 case 2806188731:
   return 0;

 case 2796990193:
   return 0;

 case 738278553:
   return 0;

 case 2622490175:
   return 0;

 case 2803122055:
   return 0;

 case 2798004553:
   return 0;

 case 3362263851:
   return 0;

 case 2803122217:
   return 0;

 case 3367381353:
   return 0;

 case 2798004559:
   return 0;

 case 3550297929:
   return 0;

 case 3550563311:
   return 0;

 case 1122536297:
   return 0;

 case 750084055:
   return 0;

 case 2421302363:
   return 0;

 case 1054068073:
   return 0;

 case 3552347553:
   return 0;

 case 750008095:
   return 0;

 case 1120335041:
   return 0;

 case 1130242155:
   return 0;

 case 816350555:
   return 0;

 case 3743666557:
   return 0;

 case 2616112801:
   return 0;

 case 3362028959:
   return 0;

 case 746781681:
   return 0;

 case 1123637091:
   return 0;

 case 2941607531:
   return 0;

 case 3223778381:
   return 0;

 case 747882313:
   return 0;

 case 3744083543:
   return 0;

 case 2615087977:
   return 0;

 case 3552348039:
   return 0;

 case 2614822595:
   return 0;

 case 750008113:
   return 0;

 case 750083569:
   return 0;

 case 1054068055:
   return 0;

 case 2613038353:
   return 0;

 case 2622263131:
   return 0;

 case 1120410515:
   return 0;

 case 3550297957:
   return 0;

 case 1123636929:
   return 0;

 case 3223778375:
   return 0;

 case 3549273105:
   return 0;

 case 2616112963:
   return 0;

 case 2803356947:
   return 0;

 case 3362028965:
   return 0;

 case 89832707:
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

 case 13940891:
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

 case 356046297:
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

 case 31374029:
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

 case 373559949:
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

 case 13860533:
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

 case 89833031:
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

 case 13940903:
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

 case 356045973:
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

 case 72826535:
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

 case 373559937:
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

 case 13917399:
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

 case 103650605:
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

 case 13959861:
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

 case 314593467:
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

 case 72826859:
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

 case 373503071:
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

 case 13917411:
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

 case 34562855:
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

 case 13865399:
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

 case 292843451:
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

 case 94576875:
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

 case 14219703:
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

 case 373200779:
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

 case 34563179:
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

 case 13865411:
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

 case 292843127:
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

 case 218934177:
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

 case 14219691:
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

 case 373371369:
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

 case 158920481:
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

 case 14036001:
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

 case 168485825:
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

 case 218934501:
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

 case 14049101:
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

 case 373371381:
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

 case 269952585:
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

 case 373441841:
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

 case 118530791:
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

 case 268889535:
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

 case 13980267:
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

 case 373440215:
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

 case 269952909:
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

 case 373441853:
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

 case 118530467:
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

 case 282707109:
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

 case 13980255:
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

 case 373459173:
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

 case 311405415:
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

 case 373498719:
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

 case 104712893:
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

 case 282707433:
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

 case 13961297:
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

 case 373459185:
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

 case 5851757:
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

