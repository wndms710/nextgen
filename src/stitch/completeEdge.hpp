
switch(current_id) {
  // =========================================
  // z-edges...
  // =========================================
 case 0:
   // complete z-edge to node 4...
   if (corner_flag[4].first == -1) {
     const int ino_new = new_node();
     x_no[ino_new][0] = -L; x_no[ino_new][1] = -L; x_no[ino_new][2] = L;
     corner_flag[4].first = ino_new;
   }
   {
     const int ied_new = new_edge();
     nooed[ied_new][0] = current_ino;
     nooed[ied_new][1] = corner_flag[4].first;
     faoed[ied_new][0] = -2;
     faoed[ied_new][1] = -4;
   }
   ++corner_flag[4].second;
   break;
 case 1:
   // complete z-edge to node 5...
   if (corner_flag[5].first == -1) {
     const int ino_new = new_node();
     x_no[ino_new][0] = L; x_no[ino_new][1] = -L; x_no[ino_new][2] = L;
     corner_flag[5].first = ino_new;
   }
   {
     const int ied_new = new_edge();
     nooed[ied_new][0] = current_ino;
     nooed[ied_new][1] = corner_flag[5].first;
     faoed[ied_new][0] = -4;
     faoed[ied_new][1] = -3;
   }
   ++corner_flag[5].second;
   break;
 case 2:
   // complete z-edge to node 6...
   if (corner_flag[6].first == -1) {
     const int ino_new = new_node();
     x_no[ino_new][0] = L; x_no[ino_new][1] = L; x_no[ino_new][2] = L;
     corner_flag[6].first = ino_new;
   }
   {
     const int ied_new = new_edge();
     nooed[ied_new][0] = current_ino;
     nooed[ied_new][1] = corner_flag[6].first;
     faoed[ied_new][0] = -3;
     faoed[ied_new][1] = -5;
   }
   ++corner_flag[6].second;
   break;
 case 3:
   // complete z-edge to node 7...
   if (corner_flag[7].first == -1) {
     const int ino_new = new_node();
     x_no[ino_new][0] = -L; x_no[ino_new][1] = L; x_no[ino_new][2] = L;
     corner_flag[7].first = ino_new;
   }
   {
     const int ied_new = new_edge();
     nooed[ied_new][0] = current_ino;
     nooed[ied_new][1] = corner_flag[7].first;
     faoed[ied_new][0] = -5;
     faoed[ied_new][1] = -2;
   }
   ++corner_flag[7].second;
   break;
   // =========================================
   // x-edges...
   // =========================================
 case 4:
   // complete x-edge to node 1...
   if (corner_flag[1].first == -1) {
     const int ino_new = new_node();
     x_no[ino_new][0] = L; x_no[ino_new][1] = -L; x_no[ino_new][2] = -L;
     corner_flag[1].first = ino_new;
   }
   {
     const int ied_new = new_edge();
     nooed[ied_new][0] = current_ino;
     nooed[ied_new][1] = corner_flag[1].first;
     faoed[ied_new][0] = -4;
     faoed[ied_new][1] = -6;
   }
   ++corner_flag[1].second;
   break;
 case 5:
   // complete x-edge to node 2...
   if (corner_flag[2].first == -1) {
     const int ino_new = new_node();
     x_no[ino_new][0] = L; x_no[ino_new][1] = L; x_no[ino_new][2] = -L;
     corner_flag[2].first = ino_new;
   }
   {
     const int ied_new = new_edge();
     nooed[ied_new][0] = current_ino;
     nooed[ied_new][1] = corner_flag[2].first;
     faoed[ied_new][0] = -6;
     faoed[ied_new][1] = -5;
   }
   ++corner_flag[2].second;
   break;
 case 6:
   // complete x-edge to node 6...
   if (corner_flag[6].first == -1) {
     const int ino_new = new_node();
     x_no[ino_new][0] = L; x_no[ino_new][1] = L; x_no[ino_new][2] = L;
     corner_flag[6].first = ino_new;
   }
   {
     const int ied_new = new_edge();
     nooed[ied_new][0] = current_ino;
     nooed[ied_new][1] = corner_flag[6].first;
     faoed[ied_new][0] = -5;
     faoed[ied_new][1] = -7;
   }
   ++corner_flag[6].second;
   break;
 case 7:
   // complete x-edge to node 5...
   if (corner_flag[5].first == -1) {
     const int ino_new = new_node();
     x_no[ino_new][0] = L; x_no[ino_new][1] = -L; x_no[ino_new][2] = L;
     corner_flag[5].first = ino_new;
   }
   {
     const int ied_new = new_edge();
     nooed[ied_new][0] = current_ino;
     nooed[ied_new][1] = corner_flag[5].first;
     faoed[ied_new][0] = -7;
     faoed[ied_new][1] = -4;
   }
   ++corner_flag[5].second;
   break;
   // =========================================
   // y-edges...
   // =========================================
 case 8:
   // complete y-edge to node 3...
   if (corner_flag[3].first == -1) {
     const int ino_new = new_node();
     x_no[ino_new][0] = -L; x_no[ino_new][1] = L; x_no[ino_new][2] = -L;
     corner_flag[3].first = ino_new;
   }
   {
     const int ied_new = new_edge();
     nooed[ied_new][0] = current_ino;
     nooed[ied_new][1] = corner_flag[3].first;
     faoed[ied_new][0] = -6;
     faoed[ied_new][1] = -2;
   }
   ++corner_flag[3].second;
   break;
 case 9:
   // complete y-edge to node 7...
   if (corner_flag[7].first == -1) {
     const int ino_new = new_node();
     x_no[ino_new][0] = -L; x_no[ino_new][1] = L; x_no[ino_new][2] = L;
     corner_flag[7].first = ino_new;
   }
   {
     const int ied_new = new_edge();
     nooed[ied_new][0] = current_ino;
     nooed[ied_new][1] = corner_flag[7].first;
     faoed[ied_new][0] = -2;
     faoed[ied_new][1] = -7;
   }
   ++corner_flag[7].second;
   break;
 case 10:
   // complete y-edge to node 6...
   if (corner_flag[6].first == -1) {
     const int ino_new = new_node();
     x_no[ino_new][0] = L; x_no[ino_new][1] = L; x_no[ino_new][2] = L;
     corner_flag[6].first = ino_new;
   }
   {
     const int ied_new = new_edge();
     nooed[ied_new][0] = current_ino;
     nooed[ied_new][1] = corner_flag[6].first;
     faoed[ied_new][0] = -7;
     faoed[ied_new][1] = -3;
   }
   ++corner_flag[6].second;
   break;
 case 11:
   // complete y-edge to node 2...
   if (corner_flag[2].first == -1) {
     const int ino_new = new_node();
     x_no[ino_new][0] = L; x_no[ino_new][1] = L; x_no[ino_new][2] = -L;
     corner_flag[2].first = ino_new;
   }
   {
     const int ied_new = new_edge();
     nooed[ied_new][0] = current_ino;
     nooed[ied_new][1] = corner_flag[2].first;
     faoed[ied_new][0] = -3;
     faoed[ied_new][1] = -6;
   }
   ++corner_flag[2].second;
   break;
 default:
   assert(0);
 }
