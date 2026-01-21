

if (!edge_flag[IED]) {
  if (corner_flag[INO0].first >= 0) {
    // building edge so flag corners
    ++corner_flag[INO0].second;
    ++corner_flag[INO1].second;

    // one of the corners is set for this edge, so set the other corner (if necessary), and build the edge...
    if (corner_flag[INO1].first == -1) {
      const int ino_new = new_node();
      x_no[ino_new][0] = XNO1; x_no[ino_new][1] = YNO1; x_no[ino_new][2] = ZNO1;
      corner_flag[INO1].first = ino_new;
      // adding a new node means we are not done...
      done = false;
    }
    const int ied_new = new_edge();
    nooed[ied_new][0] = corner_flag[INO0].first;
    nooed[ied_new][1] = corner_flag[INO1].first;
    faoed[ied_new][0] = IFA0;
    faoed[ied_new][1] = IFA1;
    edge_flag[IED] = true;
  }
  else if (corner_flag[INO1].first >= 0) {
    // building edge so flag corners
    ++corner_flag[INO0].second;
    ++corner_flag[INO1].second;

    // add INO0...
    assert(corner_flag[INO0].first == -1);
    const int ino_new = new_node();
    x_no[ino_new][0] = XNO0; x_no[ino_new][1] = YNO0; x_no[ino_new][2] = ZNO0;
    corner_flag[INO0].first = ino_new;
    // adding a new node means we are not done...
    done = false;
    const int ied_new = new_edge();
    nooed[ied_new][0] = corner_flag[INO0].first;
    nooed[ied_new][1] = corner_flag[INO1].first;
    faoed[ied_new][0] = IFA0;
    faoed[ied_new][1] = IFA1;
    edge_flag[IED] = true;
  }
 }
