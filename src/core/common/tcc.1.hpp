
 case 216655391:
   return 0;

 case 261766223:
   return 0;

 case 216416603:
   return 0;

 case 87453587:
   return 0;

 case 36722585:
   return 0;

 case 38314617:
   return 0;

 case 81441933:
   return 0;

 case 95771573:
   return 0;

 case 81354777:
   return 0;

 case 31998665:
   return 0;

 case 261293955:
   return 0;

 case 259701819:
   return 0;

 case 276176771:
   return 0;

 case 261847547:
   return 0;

 case 275937983:
   return 0;

 case 87534911:
   return 0;

 case 96243965:
   return 0;

 case 38395941:
   return 0;

 case 347518509:
   return 0;

 case 355479293:
   return 0;

 case 84674751:
   return 0;

 case 127662423:
   return 0;

 case 347430705:
   return 0;

 case 291706361:
   return 0;

 case 125300339:
   return 0;

 case 299904275:
   return 0;

 case 277359641:
   return 0;

 case 349005345:
   return 0;

 case 125060903:
   return 0;

 case 125591615:
   return 0;

 case 295361471:
   return 0;

 case 300137879:
   return 0;

 case 226341317:
   return 0;

 case 348935685:
   return 0;

 case 295122035:
   return 0;

 case 125825219:
   return 0;

 case 277241381:
   return 0;

 case 262911897:
   return 0;

 case 277013609:
   return 0;

 case 97102293:
   return 0;

 case 124945559:
   return 0;

 case 41623931:
   return 0;

 case 226223057:
   return 0;

 case 262842237:
   return 0;

 case 225995285:
   return 0;

 case 97032633:
   return 0;

 case 295006691:
   return 0;

 case 41857535:
   return 0;

 case 84670215:
   return 0;

 case 124473771:
   return 0;

 case 84897339:
   return 0;

 case 290283351:
   return 0;

 case 347504901:
   return 0;

 case 345913337:
   return 0;

 case 36726473:
   return 0;

 case 41503245:
   return 0;

 case 216656687:
   return 0;

 case 262829099:
   return 0;

 case 36802613:
   return 0;

 case 96773121:
   return 0;

 case 261411567:
   return 0;

 case 345795243:
   return 0;

 case 81481137:
   return 0;

 case 124469381:
   return 0;

 case 261336075:
   return 0;

 case 290525391:
   return 0;

 case 96247853:
   return 0;

 case 41584569:
   return 0;

 case 276178067:
   return 0;

 case 262910423:
   return 0;

 case 96323993:
   return 0;

 case 96854445:
   return 0;

 case 330217094:
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

 case 140222250:
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

 case 330101588:
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

 case 56254560:
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

 case 12248732:
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

 case 18617744:
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

 case 359123406:
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

 case 292253850:
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

 case 359197440:
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

 case 346460820:
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

 case 13444296:
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

 case 115341572:
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

 case 244123640:
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

 case 140103828:
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

 case 244008134:
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

 case 56136138:
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

 case 9060074:
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

 case 18613046:
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

 case 57200316:
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

 case 245072468:
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

 case 375132228:
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

 case 340104924:
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

 case 375162522:
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

 case 362425434:
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

 case 28218188:
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

 case 37771004:
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

 case 373974572:
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

 case 271016028:
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

 case 373739510:
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

 case 99892014:
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

 case 143293770:
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

 case 245190890:
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

 case 378320886:
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

 case 340109622:
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

 case 378351180:
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

 case 362430132:
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

 case 27267756:
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

 case 119612060:
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

 case 333580688:
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

 case 267772140:
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

 case 360152894:
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

 case 267808434:
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

 case 53839638:
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

 case 119648342:
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

 case 27499902:
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

 case 288610310:
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

 case 333348542:
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

 case 98773890:
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

 case 42041972:
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

 case 42041504:
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

 case 289045764:
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

 case 345301548:
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

 case 345379002:
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

 case 98374238:
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

 case 42118922:
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

 case 288968814:
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

 case 133574894:
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

 case 253845600:
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

 case 253845756:
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

 case 133574726:
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

 case 133465868:
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

 case 53858756:
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

 case 253954626:
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

 case 333561894:
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

 case 267808914:
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

 case 360152750:
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

 case 119611568:
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

 case 27267576:
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

 case 53749730:
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

 case 359920604:
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

 case 98810664:
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

 case 345378498:
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

 case 247198238:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[3] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   return 2;

 case 57203394:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 95166638:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[3] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   return 2;

 case 28297082:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 247316660:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[3] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   return 2;

 case 143296848:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 331165928:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[3] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   return 2;

 case 57318900:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 40959668:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[3] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   return 2;

 case 28223048:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 331284350:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[3] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   return 2;

 case 143412354:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 368802744:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[3] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   return 2;

 case 375171756:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 272078916:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[3] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   return 2;

 case 373976192:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 368807442:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[3] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   return 2;

 case 378360414:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 116403974:
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

 case 13445898:
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

 case 272079402:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[3] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   return 2;

 case 373976210:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 47310704:
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

 case 9099596:
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

 case 368807604:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[3] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   return 2;

 case 378360420:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 116402840:
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

 case 12383028:
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

 case 272080536:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 2; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[3] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   return 2;

 case 375039080:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 2; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 287527988:
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

 case 13680960:
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

 case 40960154:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[3] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   return 2;

 case 28223066:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 24990194:
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

 case 9069302:
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

 case 331284512:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[3] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   return 2;

 case 143412360:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 287526854:
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

 case 12618090:
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

 case 40999196:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 0; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[3] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   return 2;

 case 56920868:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 0; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 349648998:
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

 case 359202282:
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

 case 95167124:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[3] = (vp1[0][0]+vp1[0][1]+vp1[0][2])/(vp1[0][0]+vp1[0][1]+vp1[0][2]-vp1[1][0]-vp1[1][1]-vp1[1][2]); // edge1 wgt
   return 2;

 case 28297100:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 142229436:
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

 case 244126712:
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

 case 247316822:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[3] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   return 2;

 case 143296854:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 349609956:
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

 case 330504480:
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

 case 95206166:
   // intersection 0...
   idata[0] = EDGE_TRI_INT;
   idata[1] = 1; // edge index on tri0
   idata[2] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_EDGE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[3] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge1 wgt
   return 2;

 case 56994902:
   // intersection 0...
   idata[0] = EDGE_EDGE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = (vp1[2][0]+vp1[2][1]+vp1[2][2])/(vp1[2][0]+vp1[2][1]+vp1[2][2]-vp1[0][0]-vp1[0][1]-vp1[0][2]); // edge1 wgt
   // intersection 1...
   idata[3] = TRI_EDGE_INT;
   idata[4] = 1; // edge index on tri1
   idata[5] = 0; // unused
   ddata[2] = (vp1[1][0]+vp1[1][1]+vp1[1][2])/(vp1[1][0]+vp1[1][1]+vp1[1][2]-vp1[2][0]-vp1[2][1]-vp1[2][2]); // edge wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 130384784:
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

 case 132507476:
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

 case 254909934:
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

 case 257038782:
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

 case 130264094:
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

 case 44288282:
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

 case 255030624:
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

 case 343129128:
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

 case 353915484:
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

 case 369836948:
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

 case 90900638:
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

 case 17662434:
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

 case 33514238:
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

 case 23960850:
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

 case 296510616:
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

 case 363380744:
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

 case 353867046:
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

 case 334761830:
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

 case 90949076:
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

 case 52737552:
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

 case 52654770:
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

 case 30364814:
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

 case 363463526:
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

 case 357094878:
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

 case 127426832:
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

 case 300438632:
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

 case 258930780:
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

 case 86980560:
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

 case 260230338:
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

 case 259168758:
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

 case 128253026:
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

 case 127188854:
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

 case 17579652:
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

 case 30316376:
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

 case 343249818:
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

 case 300676610:
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

 case 363109232:
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

 case 98814552:
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

 case 53009064:
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

 case 288645140:
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

 case 259285074:
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

 case 345260886:
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

 case 127072538:
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

 case 42158306:
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

 case 128134928:
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

 case 42159584:
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

 case 260348436:
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

 case 345262200:
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

 case 259208124:
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

 case 288928152:
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

 case 127149488:
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

 case 98491040:
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

 case 363341378:
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

 case 267812802:
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

 case 52776918:
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

 case 119646890:
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

 case 24078948:
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

 case 119607680:
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

 case 334643732:
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

 case 267773604:
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

 case 255017502:
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

 case 333563190:
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

 case 130277216:
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

 case 53854220:
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

 case 90913760:
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

 case 27228372:
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

 case 353902362:
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

 case 360271010:
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

 case 296506242:
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

 case 360192098:
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

 case 33518612:
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

 case 27149496:
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

 case 90950534:
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

 case 53800434:
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

 case 353865588:
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

 case 333698948:
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

 case 254908476:
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

 case 253847052:
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

 case 130386242:
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

 case 133570358:
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

 case 132511850:
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

 case 133573430:
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

 case 257034408:
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

 case 253850136:
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

 case 343236696:
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

 case 333683880:
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

 case 17933946:
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

 case 288596702:
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

 case 30329498:
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

 case 27145590:
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

 case 300322316:
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

 case 42396284:
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

 case 357090504:
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

 case 360274880:
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

 case 87098658:
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

 case 345024222:
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

 case 30366272:
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

 case 53717652:
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

 case 300399266:
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

 case 98729018:
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

 case 343127670:
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

 case 253967742:
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

 case 17701800:
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

 case 119598452:
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

 case 44292656:
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

 case 133452740:
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

 case 369718850:
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

 case 267822042:
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

 case 128489708:
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

 case 300439928:
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

 case 259993656:
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

 case 86981856:
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

 case 23956962:
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

 case 30325610:
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

 case 334765718:
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

 case 357055674:
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

 case 363459638:
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

 case 353906250:
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

 case 52658658:
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

 case 33553442:
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

 case 24039744:
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

 case 90909872:
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

 case 334682936:
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

 case 296471412:
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

 case 128251730:
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

 case 127190150:
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

 case 260231634:
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

 case 259167462:
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

 case 296519850:
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

 case 369758054:
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

 case 33505004:
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

 case 17583540:
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

 case 132389864:
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

 case 44291360:
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

 case 257156394:
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

 case 343132206:
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

 case 254913012:
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

 case 257035704:
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

 case 130381706:
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

 case 132510554:
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

 case 357104112:
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

 case 369840836:
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

 case 86743878:
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

 case 44170670:
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

 case 301799625:
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

