if (cutCubeNodeVec[IED].empty()) {
  // if there are no intersections along this edge, we may still have to
  // add a new edge from corner to corner...
  if (corner_flag[INO0] == 0) {
    assert(corner_flag[INO1] == 0);
    const int ied = new_edge();
    nooed[ied][0] = nno_orig+INO0; // recall INO0,INO1 are in corner reference frame: 0,1,2,3,4,5,6,7 
    nooed[ied][1] = nno_orig+INO1;
    
    // these 2 nodes are not out for sure (they may still be out, but not for sure)...
    no_flag[nooed[ied][0]] = 0;
    no_flag[nooed[ied][1]] = 0;
    
    faoed[ied][0] = IFA0;
    faoed[ied][1] = IFA1;
  }
 }
 else {
   // start ino_prev off with the corner at INO0...
   int ino_prev = -1;
   if (corner_flag[INO0] == 0) 
     ino_prev = nno_orig+INO0;
   for (vector<pair<double,int> >::const_iterator iter = cutCubeNodeVec[IED].begin(); iter != cutCubeNodeVec[IED].end(); ++iter) {
     if ((ino_prev >= 0)&&(iter->second < 0)) {
       // this is an edge that should be added (orientation of the nodes is
       // correct)...
       const int ied = new_edge();
       nooed[ied][0] = ino_prev;
       nooed[ied][1] = -iter->second-1;

       // see note above and longer explanation in calling routine...
       no_flag[nooed[ied][0]] = 0;
       no_flag[nooed[ied][1]] = 0;
       
       faoed[ied][0] = IFA0;
       faoed[ied][1] = IFA1;
     }
     ino_prev = iter->second;
   }
   // if we ended positive and the second corner is in, then
   // add a last edge...
   if ((ino_prev >= 0)&&(corner_flag[INO1] == 0)) {
     const int ied = new_edge();
     nooed[ied][0] = ino_prev;
     nooed[ied][1] = nno_orig+INO1;
     
     // see above...
     no_flag[nooed[ied][0]] = 0;
     no_flag[nooed[ied][1]] = 0;

     faoed[ied][0] = IFA0;
     faoed[ied][1] = IFA1;
   }
 }
