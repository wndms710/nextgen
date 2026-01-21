if (cutCubeNodeVec[IED].empty()) {
  if (corner_flag[INO0] == -1) {
    if (corner_flag[INO1] != -1) {
      done = 0;
      corner_flag[INO1] = -1;
    }
  }
  else if (corner_flag[INO1] == -1) {
    if (corner_flag[INO0] != -1) {
      done = 0;
      corner_flag[INO0] = -1;
    }
  }
  else if (corner_flag[INO0] == 0) {
    if (corner_flag[INO1] != 0) {
      done = 0;
      corner_flag[INO1] = 0;
    }
  }
  else if (corner_flag[INO1] == 0) {
    if (corner_flag[INO0] != 0) {
      done = 0;
      corner_flag[INO0] = 0;
    }
  }
 }
