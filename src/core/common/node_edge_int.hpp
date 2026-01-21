
 case 2806188839:
   return 0;

 case 738278557:
   return 0;

 case 3362263419:
   return 0;

 case 3367381337:
   return 0;

 case 3368396091:
   return 0;

 case 3542895745:
   return 0;

 case 2613038461:
   return 0;

 case 1120410519:
   return 0;

 case 3549272673:
   return 0;

 case 2803356931:
   return 0;

 case 3543123153:
   return 0;

 case 2615087963:
   return 0;

 case 1120335149:
   return 0;

 case 816350559:
   return 0;

 case 746781249:
   return 0;

 case 2941607515:
   return 0;

 case 740176833:
   return 0;

 case 2421719363:
   return 0;

 case 202690310:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[1][0]/(e1d[1][2]+e1d[1][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 153864598:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 157095256:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[1][2]/(e1d[1][1]+e1d[1][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 230365786:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[0][2]/(e1d[0][1]+e1d[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 233518252:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 183662922:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 202726598:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[0][0]/(e1d[0][2]+e1d[0][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 180436642:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 157058968:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[0][2]/(e1d[0][1]+e1d[0][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 23057086:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[1][2]/(e1d[1][1]+e1d[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 206946208:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 149366386:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 50569664:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[2][0]/(e1d[2][2]+e1d[2][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 60122480:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[2][0]/(e0d[2][2]+e0d[2][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 364485766:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[2][2]/(e1d[2][1]+e1d[2][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 22935910:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[2][2]/(e1d[2][1]+e1d[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 327336186:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[2][2]/(e0d[2][1]+e0d[2][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 61147174:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[2][2]/(e0d[2][1]+e0d[2][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 316589278:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[1][2]/(e1d[1][1]+e1d[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 278378014:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 69768338:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[1][0]/(e1d[1][2]+e1d[1][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 317692704:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[0][0]/(e1d[0][2]+e1d[0][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 109041286:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 308139888:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 316511842:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[0][2]/(e1d[0][1]+e1d[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 222045262:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 69845774:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[0][0]/(e1d[0][2]+e1d[0][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 54962508:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[1][0]/(e1d[1][2]+e1d[1][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 165374038:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 163227808:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 26142340:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[2][2]/(e1d[2][1]+e1d[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 73906628:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[2][2]/(e0d[2][1]+e0d[2][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 332580344:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[2][0]/(e1d[2][2]+e1d[2][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 54879240:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[2][0]/(e1d[2][2]+e1d[2][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 313474764:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[2][0]/(e0d[2][2]+e0d[2][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 102643528:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[2][0]/(e0d[2][2]+e0d[2][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 234747358:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[1][1]/(e1d[1][0]+e1d[1][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 278265586:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 181541536:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[1][1]/(e1d[1][0]+e1d[1][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 205919506:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[0][1]/(e1d[0][0]+e1d[0][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 233551624:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 183629550:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 234669922:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[0][1]/(e1d[0][0]+e1d[0][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 221932834:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 181505248:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[0][1]/(e1d[0][0]+e1d[0][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 136785466:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[1][1]/(e1d[1][0]+e1d[1][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 206979580:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 149522554:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 137744944:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[2][1]/(e1d[2][0]+e1d[2][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 74059556:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[2][1]/(e0d[2][0]+e0d[2][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 250757386:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[2][1]/(e1d[2][0]+e1d[2][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 136664290:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[2][1]/(e1d[2][0]+e1d[2][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 327180018:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[2][1]/(e0d[2][0]+e0d[2][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 61303342:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[2][1]/(e0d[2][0]+e0d[2][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 331985825:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 267238519:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27799737:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359618861:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 120144223:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 266213367:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 55398635:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 93609665:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359656791:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27685991:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 293848893:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 37238833:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 55325897:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 40465571:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359729529:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359388335:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 346992987:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 98278005:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359620769:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 267276751:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 55472573:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 331908117:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 147816955:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 210905703:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27763715:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 93572081:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 331870223:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 55586283:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 183271149:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 230721373:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27690977:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 40427987:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 331715513:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 332056671:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 70605651:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 319319661:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 249043137:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 239490295:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 138377357:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 249003333:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 147930355:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 210792303:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 138455079:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 204263549:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 248965415:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 138491091:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 183157101:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 230835421:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 138609789:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 316929047:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 248810705:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 249151887:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 70491603:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 319206261:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 331925810:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[0][0]/(e1d[0][2]+e1d[0][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[1][0]/(e1d[1][2]+e1d[1][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 223660510:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 27859756:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[1][2]/(e1d[1][1]+e1d[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[0][2]/(e1d[0][1]+e1d[0][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 27859738:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[1][2]/(e1d[1][1]+e1d[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[0][2]/(e1d[0][1]+e1d[0][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 163722340:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 163721854:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 55372316:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[2][0]/(e1d[2][2]+e1d[2][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[0][0]/(e1d[0][2]+e1d[0][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 74477948:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[2][0]/(e0d[2][2]+e0d[2][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 359683114:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[0][2]/(e1d[0][1]+e1d[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[2][2]/(e1d[2][1]+e1d[2][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 27777928:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[2][2]/(e1d[2][1]+e1d[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[1][2]/(e1d[1][1]+e1d[1][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 312980718:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[2][2]/(e0d[2][1]+e0d[2][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 104200456:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[2][2]/(e0d[2][1]+e0d[2][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 332007620:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[1][0]/(e1d[1][2]+e1d[1][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[2][0]/(e1d[2][2]+e1d[2][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 283181908:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[2][0]/(e0d[2][2]+e0d[2][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 27777946:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[2][2]/(e1d[2][1]+e1d[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[1][2]/(e1d[1][1]+e1d[1][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 359683096:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[0][2]/(e1d[0][1]+e1d[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[2][2]/(e1d[2][1]+e1d[2][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 104200942:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[2][2]/(e0d[2][1]+e0d[2][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 312980232:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[2][2]/(e0d[2][1]+e0d[2][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 359560750:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[0][2]/(e1d[0][1]+e1d[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[1][2]/(e1d[1][1]+e1d[1][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 223698634:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 55494680:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[1][0]/(e1d[1][2]+e1d[1][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[0][0]/(e1d[0][2]+e1d[0][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 55494678:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[1][0]/(e1d[1][2]+e1d[1][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[0][0]/(e1d[0][2]+e1d[0][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 163760032:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 163759978:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 27737392:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[2][2]/(e1d[2][1]+e1d[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[0][2]/(e1d[0][1]+e1d[0][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 74440256:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[2][2]/(e0d[2][1]+e0d[2][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 332048174:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[0][0]/(e1d[0][2]+e1d[0][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[2][0]/(e1d[2][2]+e1d[2][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 55412868:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[2][0]/(e1d[2][2]+e1d[2][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[1][0]/(e1d[1][2]+e1d[1][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 312942594:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[2][0]/(e0d[2][2]+e0d[2][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 104238580:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[2][0]/(e0d[2][2]+e0d[2][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 359642560:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[1][2]/(e1d[1][1]+e1d[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[2][2]/(e1d[2][1]+e1d[2][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 283220032:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[2][2]/(e0d[2][1]+e0d[2][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 55412870:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[2][0]/(e1d[2][2]+e1d[2][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[1][0]/(e1d[1][2]+e1d[1][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 332048172:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[0][0]/(e1d[0][2]+e1d[0][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[2][0]/(e1d[2][2]+e1d[2][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 104238634:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[2][0]/(e0d[2][2]+e0d[2][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 312942540:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[2][0]/(e0d[2][2]+e0d[2][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 249021016:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[0][1]/(e1d[0][0]+e1d[0][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[1][1]/(e1d[1][0]+e1d[1][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 223546840:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 138399478:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[1][1]/(e1d[1][0]+e1d[1][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[0][1]/(e1d[0][0]+e1d[0][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 138399472:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[1][1]/(e1d[1][0]+e1d[1][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[0][1]/(e1d[0][0]+e1d[0][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 163873810:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 163873648:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 138277114:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[2][1]/(e1d[2][0]+e1d[2][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[0][1]/(e1d[0][0]+e1d[0][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 74591726:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[2][1]/(e0d[2][0]+e0d[2][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 249143380:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[0][1]/(e1d[0][0]+e1d[0][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[2][1]/(e1d[2][0]+e1d[2][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 138317662:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[2][1]/(e1d[2][0]+e1d[2][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[1][1]/(e1d[1][0]+e1d[1][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 312828924:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[2][1]/(e0d[2][0]+e0d[2][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 104352250:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[2][1]/(e0d[2][0]+e0d[2][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 249102826:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[1][1]/(e1d[1][0]+e1d[1][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[2][1]/(e1d[2][0]+e1d[2][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 283068238:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[2][1]/(e0d[2][0]+e0d[2][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 138317668:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[2][1]/(e1d[2][0]+e1d[2][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[1][1]/(e1d[1][0]+e1d[1][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 249143374:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[0][1]/(e1d[0][0]+e1d[0][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[2][1]/(e1d[2][0]+e1d[2][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 104352412:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[2][1]/(e0d[2][0]+e0d[2][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 312828762:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[2][1]/(e0d[2][0]+e0d[2][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 331834205:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 156699115:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27686009:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359656773:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 37239319:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 293848407:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 55512371:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 176514785:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359618879:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27799719:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 266213853:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 120143737:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 55363817:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 68100827:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359388353:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359729511:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 98278491:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 346992501:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359734497:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 350181655:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 55320929:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 332021853:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 37276903:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 293810823:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27801627:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 121207121:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 331983935:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 55434663:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 266175621:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 120181969:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 28032153:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 289142483:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 331753409:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 332094591:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 98240259:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 346954917:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 248929397:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 156585067:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 138491097:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 248965409:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 230835583:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 183156939:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 138417155:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 176628185:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 249003339:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 138377351:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 210792465:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 147930193:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 138268601:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 68214227:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 249151893:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 248810699:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 319206423:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 70491441:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 331829831:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 153510469:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27804107:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359617407:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 123332761:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 265150593:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 55510913:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 175451903:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359658245:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27681621:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 294911667:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 34050295:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 55350695:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 58534889:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359742647:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359375217:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 356558817:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 88712175:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359616399:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 264088213:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 55439027:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 331982487:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 123370345:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 265113009:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27762261:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 92509307:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 332023301:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 55316565:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 294873435:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 34088527:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27677859:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 30862157:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 332107703:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 331740297:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 356520585:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 88674591:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 248925023:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 153396421:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 138495471:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 248963951:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 234024229:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 182094057:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 138415697:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 175565303:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 249004797:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 138372977:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 211855347:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 144741547:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 138255479:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 58648289:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 249165015:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 248797577:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 328772361:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 60925503:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 332103923:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 353331961:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27681639:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359658227:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 34050781:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 294911181:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 55438001:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 122307479:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359617425:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27804089:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 265151079:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 123332275:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 55680191:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 298745897:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359375235:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359742629:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 88712661:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 356558331:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 359738867:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 353370193:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 55468199:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 331909575:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 144628309:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 211968585:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 27803081:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 122269895:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 331868765:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 55590657:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 182208267:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 233910019:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 28045271:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 298708313:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 331702391:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 332069793:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 61039713:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 328885599:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 249047511:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 242678941:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 138372983:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 249004791:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 144741709:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 211855185:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 138456537:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 205326431:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 248963957:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 138495465:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 182094219:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 234024067:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 138622911:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 326494985:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 248797583:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 249165009:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 60925665:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 328772199:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 332457980:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[1][0]/(e1d[1][2]+e1d[1][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 224192680:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 26264704:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[1][2]/(e1d[1][1]+e1d[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 70908646:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[0][2]/(e1d[0][1]+e1d[0][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 163188712:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 165375226:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 69727784:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[0][0]/(e1d[0][2]+e1d[0][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 79280600:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 316629832:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[0][2]/(e1d[0][1]+e1d[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 70831210:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[1][2]/(e1d[1][1]+e1d[1][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 308138700:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 109042474:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 332541248:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[2][0]/(e1d[2][2]+e1d[2][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 284776960:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[2][0]/(e0d[2][2]+e0d[2][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 26181436:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[2][2]/(e1d[2][1]+e1d[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 361278148:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[2][2]/(e1d[2][1]+e1d[2][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 102604432:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 0; // edge1 index on tri1
   ddata[2] = e0d[2][2]/(e0d[2][1]+e0d[2][2]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 313513860:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[2][2]/(e0d[2][1]+e0d[2][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 364363402:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[1][2]/(e1d[1][1]+e1d[1][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 238054102:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 50692028:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[1][0]/(e1d[1][2]+e1d[1][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 184693890:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[0][0]/(e1d[0][2]+e1d[0][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 149404564:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 206983846:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 157054702:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[0][2]/(e1d[0][1]+e1d[0][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 203757566:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 202730864:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[0][0]/(e1d[0][2]+e1d[0][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 184730178:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[1][0]/(e1d[1][2]+e1d[1][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 183625284:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 233555890:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 364484578:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 0; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[2][2]/(e1d[2][1]+e1d[2][2]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 326273314:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   ddata[0] = e0d[2][2]/(e0d[2][1]+e0d[2][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 50570852:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[2][0]/(e1d[2][2]+e1d[2][0]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 336850824:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 1; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[2][0]/(e1d[2][2]+e1d[2][0]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 61185352:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 1; // edge1 index on tri1
   ddata[2] = e0d[2][0]/(e0d[2][2]+e0d[2][0]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 327298008:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   ddata[0] = e0d[2][0]/(e0d[2][2]+e0d[2][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 250635022:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[1][1]/(e1d[1][0]+e1d[1][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 237897934:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 137867308:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = e1d[1][1]/(e1d[1][0]+e1d[1][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 152750566:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[0][1]/(e1d[0][0]+e1d[0][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 163341640:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 165487654:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 181500982:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[0][1]/(e1d[0][0]+e1d[0][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 203790938:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 2; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 234787912:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = e1d[0][1]/(e1d[0][0]+e1d[0][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 152673130:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = e1d[1][1]/(e1d[1][0]+e1d[1][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 308026272:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 109154902:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 250756198:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = e1d[2][1]/(e1d[2][0]+e1d[2][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 326117146:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[2][1]/(e0d[2][0]+e0d[2][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 137784040:
   // intersection 0...
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = e1d[2][1]/(e1d[2][0]+e1d[2][1]); // edge0 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 1; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 249675544:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = EDGE_NODE_INT;
   idata[4] = 2; // edge0 index on tri0
   idata[5] = 2; // node1 index on tri1
   ddata[2] = e1d[2][1]/(e1d[2][0]+e1d[2][1]); // edge0 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 102757360:
   // intersection 0...
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_EDGE_INT;
   idata[4] = 0; // node0 index on tri0
   idata[5] = 2; // edge1 index on tri1
   ddata[2] = e0d[2][1]/(e0d[2][0]+e0d[2][1]); // edge1 wgt
   ddata[3] = 0.0; // unused
   return 2;

 case 313360932:
   // intersection 0...
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   ddata[0] = e0d[2][1]/(e0d[2][0]+e0d[2][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   // intersection 1...
   idata[3] = NODE_NODE_INT;
   idata[4] = 1; // node0 index on tri0
   idata[5] = 0; // node1 index on tri1
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 2;

 case 345913338:
   return 0;

 case 347504928:
   return 0;

 case 97137120:
   return 0;

 case 290319656:
   return 0;

 case 302522690:
   return 0;

 case 111469842:
   return 0;

 case 41623932:
   return 0;

 case 124945586:
   return 0;

 case 290318178:
   return 0;

 case 97141514:
   return 0;

 case 110406420:
   return 0;

 case 305711876:
   return 0;

 case 41857536:
   return 0;

 case 295006718:
   return 0;

 case 290387838:
   return 0;

 case 290400980:
   return 0;

 case 161424744:
   return 0;

 case 170991222:
   return 0;

 case 124473754:
   return 0;

 case 84669756:
   return 0;

 case 262946718:
   return 0;

 case 124510058:
   return 0;

 case 302750300:
   return 0;

 case 111242232:
   return 0;

 case 262911880:
   return 0;

 case 277240922:
   return 0;

 case 124508592:
   return 0;

 case 262951100:
   return 0;

 case 110179134:
   return 0;

 case 305939162:
   return 0;

 case 262842220:
   return 0;

 case 226222598:
   return 0;

 case 124578252:
   return 0;

 case 124591382:
   return 0;

 case 161197458:
   return 0;

 case 170763612:
   return 0;

 case 290643480:
   return 0;

 case 347429274:
   return 0;

 case 41507128:
   return 0;

 case 345835924:
   return 0;

 case 39914966:
   return 0;

 case 291172770:
   return 0;

 case 96893802:
   return 0;

 case 125021564:
   return 0;

 case 345796534:
   return 0;

 case 41625250:
   return 0;

 case 262474308:
   return 0;

 case 126009056:
   return 0;

 case 97127406:
   return 0;

 case 295082696:
   return 0;

 case 345562930:
   return 0;

 case 345917248:
   return 0;

 case 92413176:
   return 0;

 case 350694150:
   return 0;

 case 290584430:
   return 0;

 case 304382526:
   return 0;

 case 41566200:
   return 0;

 case 345816218:
   return 0;

 case 82962308:
   return 0;

 case 276823242:
   return 0;

 case 96874118:
   return 0;

 case 110672630:
   return 0;

 case 345816240:
   return 0;

 case 41566178:
   return 0;

 case 276823836:
   return 0;

 case 82961714:
   return 0;

 case 96950258:
   return 0;

 case 165942506:
   return 0;

 case 345740100:
   return 0;

 case 345740078:
   return 0;

 case 221553960:
   return 0;

 case 221553366:
   return 0;

 case 345854306:
   return 0;

 case 304458666:
   return 0;

 case 97139324:
   return 0;

 case 290318910:
   return 0;

 case 304117472:
   return 0;

 case 110937942:
   return 0;

 case 41604266:
   return 0;

 case 110597138:
   return 0;

 case 290318924:
   return 0;

 case 97139310:
   return 0;

 case 110938320:
   return 0;

 case 304117094:
   return 0;

 case 41680406:
   return 0;

 case 165867014:
   return 0;

 case 290394416:
   return 0;

 case 290394402:
   return 0;

 case 166208172:
   return 0;

 case 166207794:
   return 0;

 case 124471590:
   return 0;

 case 83076054:
   return 0;

 case 262948904:
   return 0;

 case 124509330:
   return 0;

 case 304344596:
   return 0;

 case 110710818:
   return 0;

 case 262911174:
   return 0;

 case 276710102:
   return 0;

 case 124509320:
   return 0;

 case 262948914:
   return 0;

 case 110710548:
   return 0;

 case 304344866:
   return 0;

 case 262835682:
   return 0;

 case 221440250:
   return 0;

 case 124584812:
   return 0;

 case 124584822:
   return 0;

 case 165980400:
   return 0;

 case 165980670:
   return 0;

 case 419459085:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 3695383239:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 732451299:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 1455846195:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 2469017221:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 922821469:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 654433173:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 921722299:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 414572415:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 1137967311:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 947597141:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 3696368685:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 2666944583:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 3418698209:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 723395293:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 1147023317:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 2748772799:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 3416613107:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 1137967329:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 3696369171:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 1451035343:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 654243631:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 3962289049:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 1485973721:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 1455846213:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 922821955:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][2]/(e0d[1][1]+e0d[1][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 1216174979:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 419383267:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 384444889:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][0]/(e0d[1][2]+e0d[1][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 2203096857:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 1147023335:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 0; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[0][0]+vp0[0][1]+vp0[0][2])/(vp0[0][0]+vp0[0][1]+vp0[0][2]-vp0[1][0]-vp0[1][1]-vp0[1][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 3416613593:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][2]/(e0d[0][1]+e0d[0][2]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 3498175969:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 2667209937:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 1; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[1][0]+vp0[1][1]+vp0[1][2])/(vp0[1][0]+vp0[1][1]+vp0[1][2]-vp0[2][0]-vp0[2][1]-vp0[2][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 1818652899:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 51765711:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][0]/(e0d[0][2]+e0d[0][0]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 2558825885:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 2200140519:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 3606560027:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 2544971855:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 3965245549:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 1482676049:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 2544971861:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 1482676211:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 3620414051:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 2558825879:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 387742561:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[1][1]/(e0d[1][0]+e0d[1][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 2200140357:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 2932006439:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 58021017:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 3233379473:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 2932006433:
   idata[0] = EDGE_NODE_INT;
   idata[1] = 2; // edge0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = (vp0[2][0]+vp0[2][1]+vp0[2][2])/(vp0[2][0]+vp0[2][1]+vp0[2][2]-vp0[0][0]-vp0[0][1]-vp0[0][2]); // edge0 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 1812397755:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 58020855:
   idata[0] = NODE_EDGE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // edge1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = e0d[0][1]/(e0d[0][0]+e0d[0][1]); // edge1 wgt
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 746781519:
   return 0;

 case 2941607525:
   return 0;

 case 740176455:
   return 0;

 case 2421719349:
   return 0;

 case 1120334555:
   return 0;

 case 816350537:
   return 0;

 case 3362263689:
   return 0;

 case 3367381347:
   return 0;

 case 3368395713:
   return 0;

 case 3542895731:
   return 0;

 case 2806188245:
   return 0;

 case 738278535:
   return 0;

 case 3549272943:
   return 0;

 case 2803356941:
   return 0;

 case 3543122775:
   return 0;

 case 2615087949:
   return 0;

 case 2613037867:
   return 0;

 case 1120410497:
   return 0;

 case 97121247:
   return 0;

 case 290831075:
   return 0;

 case 97203219:
   return 0;

 case 350352479:
   return 0;

 case 290340213:
   return 0;

 case 126349737:
   return 0;

 case 41604977:
   return 0;

 case 111128093:
   return 0;

 case 41686949:
   return 0;

 case 170649497:
   return 0;

 case 345856475:
   return 0;

 case 306052503:
   return 0;

 case 262930851:
   return 0;

 case 291058847:
   return 0;

 case 263012823:
   return 0;

 case 350580251:
   return 0;

 case 124530633:
   return 0;

 case 126122613:
   return 0;

 case 41584835:
   return 0;

 case 96247755:
   return 0;

 case 290565783:
   return 0;

 case 96777593:
   return 0;

 case 291096755:
   return 0;

 case 39990993:
   return 0;

 case 41502863:
   return 0;

 case 36726351:
   return 0;

 case 290647755:
   return 0;

 case 290293451:
   return 0;

 case 350618159:
   return 0;

 case 92337563:
   return 0;

 case 345794645:
   return 0;

 case 261411437:
   return 0;

 case 96895701:
   return 0;

 case 290526407:
   return 0;

 case 126084705:
   return 0;

 case 262398671:
   return 0;

 case 262910709:
   return 0;

 case 276178509:
   return 0;

 case 124509773:
   return 0;

 case 262947327:
   return 0;

 case 111241817:
   return 0;

 case 302750727:
   return 0;

 case 262828737:
   return 0;

 case 216657105:
   return 0;

 case 124591745:
   return 0;

 case 124578645:
   return 0;

 case 170763221:
   return 0;

 case 161197877:
   return 0;

 case 124468803:
   return 0;

 case 81481547:
   return 0;

 case 262951679:
   return 0;

 case 124508337:
   return 0;

 case 305938779:
   return 0;

 case 110179529:
   return 0;

 case 97101105:
   return 0;

 case 275950737:
   return 0;

 case 345835635:
   return 0;

 case 41507741:
   return 0;

 case 291172247:
   return 0;

 case 39915501:
   return 0;

 case 97019133:
   return 0;

 case 216429333:
   return 0;

 case 345917607:
   return 0;

 case 345563327:
   return 0;

 case 350693651:
   return 0;

 case 92413703:
   return 0;

 case 290278383:
   return 0;

 case 81708671:
   return 0;

 case 41625825:
   return 0;

 case 345796283:
   return 0;

 case 126008565:
   return 0;

 case 262474811:
   return 0;

 case 276766344:
   return 0;

 case 304363576:
   return 0;

 case 83019196:
   return 0;

 case 304363816:
   return 0;

 case 83018572:
   return 0;

 case 276767000:
   return 0;

 case 110691612:
   return 0;

 case 304363792:
   return 0;

 case 83019220:
   return 0;

 case 276766352:
   return 0;

 case 110767266:
   return 0;

 case 165961470:
   return 0;

 case 304288138:
   return 0;

 case 304288162:
   return 0;

 case 221496494:
   return 0;

 case 221497142:
   return 0;

 case 221496474:
   return 0;

 case 304287598:
   return 0;

 case 221497122:
   return 0;

 case 165961706:
   return 0;

 case 304287622:
   return 0;

 case 110767814:
   return 0;

 case 165961686:
   return 0;

 case 110767274:
   return 0;

 case 166037124:
   return 0;

 case 166037664:
   return 0;

 case 166037684:
   return 0;

 case 166037144:
   return 0;

 case 83018548:
   return 0;

 case 276766992:
   return 0;

 case 110691836:
   return 0;

 case 304363600:
   return 0;

 case 110691828:
   return 0;

 case 110691620:
   return 0;

 case 110767806:
   return 0;

 case 165961490:
   return 0;

 case 304117025:
   return 0;

 case 97139631:
   return 0;

 case 82961699:
   return 0;

 case 41566501:
   return 0;

 case 304344635:
   return 0;

 case 262949229:
   return 0;

 case 110937865:
   return 0;

 case 290319015:
   return 0;

 case 276823219:
   return 0;

 case 345816325:
   return 0;

 case 110710579:
   return 0;

 case 124509429:
   return 0;

 case 166207723:
   return 0;

 case 290394669:
   return 0;

 case 221553349:
   return 0;

 case 345740347:
   return 0;

 case 165980437:
   return 0;

 case 124585083:
   return 0;

 case 165924014:
   return 0;

 case 83132890:
   return 0;

 case 165923366:
   return 0;

 case 221458782:
   return 0;

 case 83132866:
   return 0;

 case 276652674:
   return 0;

 case 221459018:
   return 0;

 case 276653222:
   return 0;

 case 221458802:
   return 0;

 case 165923346:
   return 0;

 case 276653214:
   return 0;

 case 83132326:
   return 0;

 case 221383364:
   return 0;

 case 221382824:
   return 0;

 case 221382804:
   return 0;

 case 221383344:
   return 0;

 case 304401940:
   return 0;

 case 110653496:
   return 0;

 case 276728652:
   return 0;

 case 83056888:
   return 0;

 case 83056912:
   return 0;

 case 110654144:
   return 0;

 case 276728660:
   return 0;

 case 110653488:
   return 0;

 case 276728868:
   return 0;

 case 83056672:
   return 0;

 case 276652682:
   return 0;

 case 221458998:
   return 0;

 case 304401292:
   return 0;

 case 304401916:
   return 0;

 case 276728876:
   return 0;

 case 83056696:
   return 0;

 case 304401268:
   return 0;

 case 110654136:
   return 0;

 case 83132350:
   return 0;

 case 165923994:
   return 0;

 case 2613038434:
   return 0;

 case 1120410518:
   return 0;

 case 750083488:
   return 0;

 case 740176374:
   return 0;

 case 1054068052:
   return 0;

 case 2421719346:
   return 0;

 case 3549273132:
   return 0;

 case 2803356948:
   return 0;

 case 1123636902:
   return 0;

 case 746781276:
   return 0;

 case 3223778374:
   return 0;

 case 2941607516:
   return 0;

 case 2615088220:
   return 0;

 case 2614822604:
   return 0;

 case 747882070:
   return 0;

 case 1120334528:
   return 0;

 case 3744083534:
   return 0;

 case 816350536:
   return 0;

 case 1120335122:
   return 0;

 case 816350558:
   return 0;

 case 3552347472:
   return 0;

 case 3543122694:
   return 0;

 case 750008092:
   return 0;

 case 2615087946:
   return 0;

 case 746781708:
   return 0;

 case 2941607532:
   return 0;

 case 2616112774:
   return 0;

 case 3549272700:
   return 0;

 case 3362028958:
   return 0;

 case 2803356932:
   return 0;

 case 1122536540:
   return 0;

 case 2421302372:
   return 0;

 case 3550297686:
   return 0;

 case 2613037840:
   return 0;

 case 3550563302:
   return 0;

 case 1120410496:
   return 0;

 case 2806188812:
   return 0;

 case 738278556:
   return 0;

 case 3359197094:
   return 0;

 case 3368395632:
   return 0;

 case 1132140054:
   return 0;

 case 3542895728:
   return 0;

 case 3362263878:
   return 0;

 case 3367381354:
   return 0;

 case 2803122028:
   return 0;

 case 3362263446:
   return 0;

 case 2798004552:
   return 0;

 case 3367381338:
   return 0;

 case 2804144534:
   return 0;

 case 3543085282:
   return 0;

 case 3361241372:
   return 0;

 case 2806188218:
   return 0;

 case 2622300624:
   return 0;

 case 738278534:
   return 0;

 case 750008114:
   return 0;

 case 3552348066:
   return 0;

 case 1132140076:
   return 0;

 case 3359197688:
   return 0;

 case 1054068074:
   return 0;

 case 750084082:
   return 0;

 case 3550297960:
   return 0;

 case 2622263212:
   return 0;

 case 2622490178:
   return 0;

 case 2796990274:
   return 0;

 case 3743666560:
   return 0;

 case 1130242236:
   return 0;

 case 3362028974:
   return 0;

 case 2616113206:
   return 0;

 case 2798004568:
   return 0;

 case 2803122460:
   return 0;

 case 3223778390:
   return 0;

 case 1123637334:
   return 0;

 case 107977418:
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 69235780:
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 108091304:
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 152140582:
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 204851528:
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 179908194:
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 169624580:
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 69320506:
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 169738466:
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 152225308:
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 321768554:
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 180068736:
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 310266178:
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 1; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 318227404:
   idata[0] = NODE_NODE_INT;
   idata[1] = 1; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 310152616:
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 235322614:
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 241027432:
   idata[0] = NODE_NODE_INT;
   idata[1] = 0; // node0 index on tri0
   idata[2] = 2; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 207592914:
   idata[0] = NODE_NODE_INT;
   idata[1] = 2; // node0 index on tri0
   idata[2] = 0; // node1 index on tri1
   idata[3] = 0; // unused
   idata[4] = 0; // unused
   ddata[0] = 0.0; // unused
   ddata[1] = 0.0; // unused
   ddata[2] = 0.0; // unused
   ddata[3] = 0.0; // unused
   return 1;

 case 348928723:
   return 0;

 case 221557713:
   return 0;

 case 95887377:
   return 0;

 case 125648081:
   return 0;

 case 165941075:
   return 0;

 case 165982115:
   return 0;

 case 127659619:
   return 0;

 case 83079785:
   return 0;

 case 261886617:
   return 0;

 case 125571941:
   return 0;

 case 304343187:
   return 0;

 case 110712239:
   return 0;

 case 349004863:
   return 0;

 case 276827589:
   return 0;

 case 95811237:
   return 0;

 case 291647321:
   return 0;

 case 110671199:
   return 0;

 case 304384227:
   return 0;

 case 125647695:
   return 0;

 case 165981885:
   return 0;

 case 259647013:
   return 0;

 case 348929125:
   return 0;

 case 221435255:
   return 0;

 case 221558375:
   return 0;

 case 291646719:
   return 0;

 case 304383989:
   return 0;

 case 38378125:
   return 0;

 case 349004617:
   return 0;

 case 82957335:
   return 0;

 case 276828227:
   return 0;

 case 125572203:
   return 0;

 case 110712033:
   return 0;

 case 259722505:
   return 0;

 case 127660237:
   return 0;

 case 276705107:
   return 0;

 case 83080455:
   return 0;

 case 299960337:
   return 0;

 case 166220835:
   return 0;

 case 32114451:
   return 0;

 case 299960735:
   return 0;

 case 165853433:
   return 0;

 case 166221389:
   return 0;

 case 355419633:
   return 0;

 case 304471307:
   return 0;

 case 87573963:
   return 0;

 case 299884595:
   return 0;

 case 304103913:
   return 0;

 case 110951513:
   return 0;

 case 299884845:
   return 0;

 case 110950983:
   return 0;

 case 32038311:
   return 0;

 case 355420247:
   return 0;

 case 110583557:
   return 0;

 case 304471869:
   return 0;

 case 276705103:
   return 0;

 case 259722397:
   return 0;

 case 83080443:
   return 0;

 case 304343393:
   return 0;

 case 127659913:
   return 0;

 case 261886355:
   return 0;

 case 110671195:
   return 0;

 case 95811129:
   return 0;

 case 304384215:
   return 0;

 case 82957973:
   return 0;

 case 291646997:
   return 0;

 case 38377879:
   return 0;

 case 110583553:
   return 0;

 case 32038203:
   return 0;

 case 304471857:
   return 0;

 case 304104443:
   return 0;

 case 355419923:
   return 0;

 case 87573713:
   return 0;

 case 221435245:
   return 0;

 case 259646743:
   return 0;

 case 221558365:
   return 0;

 case 165941287:
   return 0;

 case 348928855:
   return 0;

 case 95887277:
   return 0;

 case 165941065:
   return 0;

 case 95887107:
   return 0;

 case 165982105:
   return 0;

 case 221435899:
   return 0;

 case 125647811:
   return 0;

 case 259646929:
   return 0;

 case 165853423:
   return 0;

 case 32114181:
   return 0;

 case 166221379:
   return 0;

 case 165853969:
   return 0;

 case 299960465:
   return 0;

 case 32114363:
   return 0;

 case 82957323:
   return 0;

 case 38377801:
   return 0;

 case 276828223:
   return 0;

 case 110671429:
   return 0;

 case 349004509:
   return 0;

 case 95811623:
   return 0;

 case 304343175:
   return 0;

 case 261886293:
   return 0;

 case 110712235:
   return 0;

 case 276705769:
   return 0;

 case 125571833:
   return 0;

 case 259722907:
   return 0;

 case 304103901:
   return 0;

 case 87573639:
   return 0;

 case 110951509:
   return 0;

 case 110584111:
   return 0;

 case 299884487:
   return 0;

 case 32038709:
   return 0;

 case 263012435:
   return 0;

 case 350579967:
   return 0;

 case 124408059:
   return 0;

 case 124762331:
   return 0;

 case 36840683:
   return 0;

 case 295120415:
   return 0;

 case 124530029:
   return 0;

 case 126122321:
   return 0;

 case 262890465:
   return 0;

 case 124528727:
   return 0;

 case 261298329:
   return 0;

 case 125059283:
   return 0;

 case 262931111:
   return 0;

 case 291058587:
   return 0;

 case 124489383:
   return 0;

 case 263008541:
   return 0;

 case 96362063:
   return 0;

 case 347391177:
   return 0;

 case 97202855:
   return 0;

 case 350352843:
   return 0;

 case 345733913:
   return 0;

 case 345747021:
   return 0;

 case 216770897:
   return 0;

 case 226336457:
   return 0;

 case 290339633:
   return 0;

 case 126350093:
   return 0;

 case 41564603:
   return 0;

 case 345816681:
   return 0;

 case 81367899:
   return 0;

 case 277354781:
   return 0;

 case 97121531:
   return 0;

 case 290831463:
   return 0;

 case 345815237:
   return 0;

 case 41568963:
   return 0;

 case 276292277:
   return 0;

 case 84556167:
   return 0;

 case 41686593:
   return 0;

 case 170650077:
   return 0;

 case 290464037:
   return 0;

 case 290477169:
   return 0;

 case 216694757:
   return 0;

 case 226260965:
   return 0;

 case 345855903:
   return 0;

 case 306053075:
   return 0;

 case 96834455:
   return 0;

 case 290546829:
   return 0;

 case 81443391:
   return 0;

 case 277279289:
   return 0;

 case 41605269:
   return 0;

 case 111128697:
   return 0;

 case 290545361:
   return 0;

 case 96838839:
   return 0;

 case 276216137:
   return 0;

 case 84632307:
   return 0;

 case 41604517:
   return 0;

 case 110596635:
   return 0;

 case 41680009:
   return 0;

 case 165866487:
   return 0;

 case 345853693:
   return 0;

 case 304458131:
   return 0;

 case 262911429:
   return 0;

 case 276709707:
   return 0;

 case 262835289:
   return 0;

 case 221439831:
   return 0;

 case 124470981:
   return 0;

 case 83075627:
   return 0;

 case 97101831:
   return 0;

 case 276482097:
   return 0;

 case 97025691:
   return 0;

 case 221212221:
   return 0;

 case 290280567:
   return 0;

 case 83302913:
   return 0;
