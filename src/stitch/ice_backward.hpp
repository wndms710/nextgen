
// ice_backward.hpp
// define ICE 

assert(nooce[ICE][1] == icorner);
assert(cdoce[ICE] == icd);
if (cutCubeNodeVec[ICE].empty()) {
  if (corner_flag[nooce[ICE][0]][icd] == 2) {
    // someone else already added this cube edge.
    assert(corner_flag[nooce[ICE][1]][icd] == 2);
  }
  else {
    assert((corner_flag[nooce[ICE][0]][icd] == 0)||(corner_flag[nooce[ICE][0]][icd] == 1));
    assert(corner_flag[nooce[ICE][1]][icd] == 1);
    corner_flag[nooce[ICE][1]][icd] = 2;
    corner_flag[nooce[ICE][0]][icd] = 2; 
    const int ied_new = new_edge();
    faoed[ied_new][0] = faoce[ICE][0];
    faoed[ied_new][1] = faoce[ICE][1];
    nooed[ied_new][0] = nooce[ICE][0]+nno_cube;
    nooed[ied_new][1] = nooce[ICE][1]+nno_cube;
    // we need to travel from this corner in the 2 other directions...
    if (corner_flag[nooce[ICE][0]][(icd+1)%3] == 0) {
      corner_flag[nooce[ICE][0]][(icd+1)%3] = 1; 
      cornerStack.push(pair<int,int>(nooce[ICE][0],(icd+1)%3));
    }
    if (corner_flag[nooce[ICE][0]][(icd+2)%3] == 0) {
      corner_flag[nooce[ICE][0]][(icd+2)%3] = 1; 
      cornerStack.push(pair<int,int>(nooce[ICE][0],(icd+2)%3));
    }
  }
 }
 else if (cutCubeNodeVec[ICE].back().second.second == 1) {
   // add the edge...
   const int ied_new = new_edge();
   faoed[ied_new][0] = faoce[ICE][0];
   faoed[ied_new][1] = faoce[ICE][1];
   nooed[ied_new][0] = cutCubeNodeVec[ICE].back().second.first;
   nooed[ied_new][1] = nooce[ICE][1]+nno_cube;
   cutCubeNodeVec[ICE].back().second.second = 2; // no need to ever start marching from this one...
 }
 else {
   // TODO: get rid of this...
   //if (!(cutCubeNodeVec[ICE].back().second.second <= -1))
   //  cout << "ice_backward: cutCubeNodeVec[ICE].back().second.second: " << cutCubeNodeVec[ICE].back().second.second << endl;
   assert(cutCubeNodeVec[ICE].back().second.second <= -1);
   //cout << "setting ierr = -1" << endl;
   ierr = -1;
 }
