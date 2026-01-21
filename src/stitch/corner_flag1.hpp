if (!cutCubeNodeVec[IED].empty()) {
  // sort if there are 2 or more entries along this edge...
  if (cutCubeNodeVec[IED].size() > 1) 
    sort(cutCubeNodeVec[IED].begin(),cutCubeNodeVec[IED].end());
  if (cutCubeNodeVec[IED].front().second < 0) {
    corner_flag[INO0] = min(corner_flag[INO0],0); // in
  }
  else {
    corner_flag[INO0] = -1; // out
  }
  if (cutCubeNodeVec[IED].back().second < 0) {
    corner_flag[INO1] = -1; // out
  }
  else {
    corner_flag[INO1] = min(corner_flag[INO1],0); // in
  }
 }
