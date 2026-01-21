#include "SimpleSurface.hpp"

void SimpleSurface::buildTeost() {

  COUT1("SimpleSurface::buildTeost()");

  assert(!b_teost);
  b_teost = teost.build(nsp,nst,spost);
}

void SimpleSurface::clearTeost() {

  COUT1("SimpleSurface::clearTeost()");

  b_teost = false;
  // for now, delete and clear everything. Could do memory management later...
  teost.clear();

  // open edge groups depend on teost, so they must be cleared as well
  clearOpenEdgeGroups();

  // sposp is constructed from teost, so it should be updated as well
  clearSposp();
}

void SimpleSurface::ensureTeost() {

  if (!b_teost) {
    buildTeost();
    assert(b_teost);
  }

}

bool SimpleSurface::gotTeost() const {

  return b_teost;

}

bool SimpleSurface::getTriNbrData(int& ist_nbr,int& i_nbr,int& orient_nbr,const int ist,const int i) const {
  return teost.getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i);
}

int SimpleSurface::getTriNbrDataFull(int& ist_nbr,int& i_nbr,int& orient_nbr,const int ist,const int i) const {
  return teost.getTriNbrDataFull(ist_nbr,i_nbr,orient_nbr,ist,i);
}

void SimpleSurface::setTriNbrData(const int ist_nbr,const int i_nbr,const int orient_nbr,const int ist,const int i) {
  teost.setTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i);
}

void SimpleSurface::setTriNbrOpen(const int ist,const int i) {
  teost.setTriNbrOpen(ist,i);
}

void SimpleSurface::setTriNbrMulti(const int ime,const int orient_ime,const int ist,const int i) {
  teost.setTriNbrMulti(ime,orient_ime,ist,i);
}

void SimpleSurface::setTriNbrToSameAs(const int ist_ref,const int i_ref,const int ist,const int i,const bool b_flip_orient) {
  teost.setTriNbrToSameAs(ist_ref,i_ref,ist,i,b_flip_orient);
}

void SimpleSurface::getTriNbrs(int (&ist_nbr)[3],const int ist) const {
  assert(b_teost);
  teost.getTriNbrs(ist_nbr,ist);
}

bool SimpleSurface::getAlignedTriNbr(int& ist_nbr,const int ist,const int i) const {
  return teost.getAlignedTriNbr(ist_nbr,ist,i);
}

bool SimpleSurface::isNbrAligned(const int ist,const int i) const {
  return teost.isNbrAligned(ist,i);
}

bool SimpleSurface::isEdgeOpen(const int ist,const int i) const {
  return teost.isEdgeOpen(ist,i);
}

bool SimpleSurface::isEdgeMulti(const int ist,const int i) const {
  int ime,orient_ime;
  return teost.isEdgeMulti(ime,orient_ime,ist,i);
}

bool SimpleSurface::isEdgeMulti(int& ime,int& orient_ime,const int ist,const int i) const {
  return teost.isEdgeMulti(ime,orient_ime,ist,i);
}

void SimpleSurface::setOpenEdgeGroup(const int igr,const int ist,const int i) {
  teost.setOpenEdgeGroup(igr,ist,i);
}

bool SimpleSurface::getOpenEdgeGroup(int& igr,const int ist,const int i) const {
  return teost.getOpenEdgeGroup(igr,ist,i);
}

bool SimpleSurface::isEdgeZoneBoundary(const int ist, const int i) const {
  return teost.isEdgeZoneBoundary(ist,i,znost);
}

bool SimpleSurface::isEdgeSubzoneBoundary(const int ist, const int i) const {
  return teost.isEdgeZoneBoundary(ist,i,szost);  // use szost as boundary flag instead of znost
}

bool SimpleSurface::isEdgeSelectionBoundary(const int ist, const int i) const {
  return teost.isEdgeZoneBoundary(ist,i,st_flag);  // use st_flag
}

bool SimpleSurface::isEdgeMisaligned(const int ist, const int i) const {
  if (teost.isEdgeOpen(ist,i) || teost.isEdgeMulti(ist,i)) return false;

  return (teost.isNbrAligned(ist,i)) ? false:true;
}

// bool SimpleSurface::isEdgeMetalBoundary(const int ist, const int i) const {
//   if (!zoneVec[znost[ist]].isMetal()) return false;  // only report half edge from metal zone
//
//   int ist_nbr = -1;
//   if (teost.getTriNbr(ist_nbr,ist,i)) {
//     return (zoneVec[znost[ist_nbr]].isMetal()) ? false:true;
//   }
//   else if (teost.isEdgeMulti(ist,i)) {
//     vector<pair<int,int> > ist_nbrs;  // <ist,orient> for ime nbrs
//     teost.getMultiNbrs(ist_nbrs,ist,i);
//
//     int count = 0;
//     for (int ii=0,ii_end=ist_nbrs.size(); ii < ii_end; ++ii) {
//       int ist_nbr = ist_nbrs[ii].first;
//       if ((ist_nbr != ist) && zoneVec[znost[ist_nbr]].isMetal()) ++count;
//     }
//
//     return (count == 1) ? false:true;  // even if multi, if there one metal nbr then edge is not boundary
//   }
//
//   assert(isEdgeOpen(ist,i));
//   return false;  // must be a metal open edge, not a metal boundary...?
// }
//
// bool SimpleSurface::isEdgeFluidBoundary(const int ist, const int i) const {
//   if (!zoneVec[znost[ist]].isFluid()) return false;  // only report half edge from fluid zone
//
//   int ist_nbr = -1;
//   if (teost.getTriNbr(ist_nbr,ist,i)) {
//     return (zoneVec[znost[ist_nbr]].isFluid()) ? false:true;
//   }
//   else if (teost.isEdgeMulti(ist,i)) {
//     vector<pair<int,int> > ist_nbrs;  // <ist,orient> for ime nbrs
//     teost.getMultiNbrs(ist_nbrs,ist,i);
//
//     int count = 0;
//     for (int ii=0,ii_end=ist_nbrs.size(); ii < ii_end; ++ii) {
//       int ist_nbr = ist_nbrs[ii].first;
//       if ((ist_nbr != ist) && zoneVec[znost[ist_nbr]].isFluid()) ++count;
//     }
//
//     return (count == 1) ? false:true;  // even if multi, if there one fluid nbr then edge is not boundary
//   }
//
//   assert(isEdgeOpen(ist,i));
//   return false;  // must be a fluid open edge...
// }

bool SimpleSurface::isEdgeNmeBoundary(const int ist,const int i) const {
  if (teost.isEdgeOpen(ist,i) || teost.isEdgeMulti(ist,i)) return false;

  int ist_nbr,i_nbr,orient_nbr;
  getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i);

  return (nmTriData.isClean(ist) != nmTriData.isClean(ist_nbr)) ? true:false;

}

bool SimpleSurface::isEdgeCrease(const int ist,const int i) const {
  if (teost.isEdgeOpen(ist,i) || teost.isEdgeMulti(ist,i)) return false;

  int ist_nbr,i_nbr,orient_nbr;
  getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i);
  return isEdgeCrease(ist,ist_nbr,orient_nbr);

}

bool SimpleSurface::isEdgeCrease(const int ist0,const int ist1, const bool misaligned) const {
  double unit_n_ist[3];
  double ist_n[3] = TRI_NORMAL_2(xsp[spost[ist0][0]],xsp[spost[ist0][1]],xsp[spost[ist0][2]]);
  if (misaligned) {
    FOR_I3 ist_n[i] *= -1.0;
  }
  double mag = MAG(ist_n);
  if (mag <= 0.0) {
    return false;  // linear tri, cannot compute
  }
  else {
    FOR_I3 unit_n_ist[i] = ist_n[i]/mag;
  }

  double unit_n_nbr[3];
  double nbr_n[3] = TRI_NORMAL_2(xsp[spost[ist1][0]],xsp[spost[ist1][1]],xsp[spost[ist1][2]]);
  if (misaligned) {
    FOR_I3 nbr_n[i] *= -1.0;
  }
  mag = MAG(nbr_n);
  if (mag <= 0.0) {
    return false;  // linear tri, cannot compute
  }
  else {
    FOR_I3 unit_n_nbr[i] = nbr_n[i]/mag;
  }

  return (DOT_PRODUCT(unit_n_ist,unit_n_nbr) >= feature_cos) ? false:true;
}

double SimpleSurface::edgeCreaseCos(const int ist,const int i) const {

  int ist_nbr,i_nbr,orient_nbr;
  getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i);
  return edgeCreaseCos(ist,ist_nbr,orient_nbr);

}

double SimpleSurface::edgeCreaseCos(const int ist0,const int ist1, const bool misaligned) const {

  double unit_n_ist[3];
  double ist_n[3] = TRI_NORMAL_2(xsp[spost[ist0][0]],xsp[spost[ist0][1]],xsp[spost[ist0][2]]);
  if (misaligned) FOR_I3 ist_n[i] *= -1.0;
  double mag = MAG(ist_n);
  FOR_I3 unit_n_ist[i] = ist_n[i]/mag;

  double unit_n_nbr[3];
  double nbr_n[3] = TRI_NORMAL_2(xsp[spost[ist1][0]],xsp[spost[ist1][1]],xsp[spost[ist1][2]]);
  if (misaligned) FOR_I3 nbr_n[i] *= -1.0;
  mag = MAG(nbr_n);
  FOR_I3 unit_n_nbr[i] = nbr_n[i]/mag;

  return DOT_PRODUCT(unit_n_ist,unit_n_nbr);

}


bool SimpleSurface::isEdgeCrease(const double unit_n_ist[3],const int ist1, const bool misaligned) const {
  double unit_n_nbr[3];
  double nbr_n[3] = TRI_NORMAL_2(xsp[spost[ist1][0]],xsp[spost[ist1][1]],xsp[spost[ist1][2]]);
  if (misaligned) {
    FOR_I3 nbr_n[i] *= -1.0;
  }
  double mag = MAG(nbr_n);
  if (mag <= 0.0) {
    return false;  // linear tri, cannot compute
  }
  else {
    FOR_I3 unit_n_nbr[i] = nbr_n[i]/mag;
  }

  return (DOT_PRODUCT(unit_n_ist,unit_n_nbr) >= feature_cos) ? false:true;
}

bool SimpleSurface::isEdgeFeature(const int ist, const int i) const {
  double unit_n_ist[3];
  const double ist_n[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
  double mag = MAG(ist_n);
  if (mag <= 0.0) {
    return false;  // linear tri, treat as not a feature
  }
  else {
    FOR_I3 unit_n_ist[i] = ist_n[i]/mag;
  }

  int ist_nbr, i_nbr, orient_nbr;
  if (teost.getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
    // single valid neighbor
    return isEdgeCrease(unit_n_ist,ist_nbr,bool(orient_nbr));
  }

  // if reached here then is either multi or open edge
  int ime, orient_ime;
  if (teost.isEdgeMulti(ime,orient_ime,ist,i)) {
    vector<pair<int,int> > ist_nbrs;  // <ist,orient> for ime nbrs
    teost.getMultiNbrs(ist_nbrs,ist,i);

    for (int _ist=0,_nst=ist_nbrs.size(); _ist < _nst; ++_ist) {
      int ist_nbr = ist_nbrs[_ist].first;
      if ((ist_nbr != ist)) {
        int orient_nbr = ist_nbrs[_ist].second;
        const bool misaligned = (orient_nbr != orient_ime) ? true:false;
        if (isEdgeCrease(unit_n_ist,ist_nbr,misaligned)) return true;
      }
    }
  }

  // no mutli-edge crease or this is an open edge
  return false;
}

void SimpleSurface::multiEdgesToDynamicEdges() {
  ensureTeost();
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      int ime, orient_ime;
      if (isEdgeMulti(ime,orient_ime,ist,i)) {
        eoi_to_group.insert(pair<uint,int> (packEdge(ist,i),ime));
      }
    }
  }
}

void SimpleSurface::flagTrisTouchingMultiEdges() {
  cout << "SimpleSurface::flagTrisTouchingMultiEdges()" << endl;
  ensureTeost();
  st_flag.resize(nst);
  st_flag.setAll(0);
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      int ime,orient_ime;
      if (isEdgeMulti(ime,orient_ime,ist,i)) {
        cout << " > ist: " << ist << " touching ime: " << ime << " nodes: " << spost[ist][i] << " " << spost[ist][(i+1)%3] <<
          " coords: " << COUT_VEC(xsp[spost[ist][i]]) << " " << COUT_VEC(xsp[spost[ist][(i+1)%3]]) << endl;
        st_flag[ist] = 1;
      }
    }
  }
}

void SimpleSurface::flagTrisTouchingMultiEdges(set<int>& specificMultiEdges) {
  cout << "SimpleSurface::flagTrisTouchingMultiEdges()" << endl;
  ensureTeost();
  st_flag.resize(nst);
  st_flag.setAll(0);
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      int ime,orient_ime;
      if (isEdgeMulti(ime,orient_ime,ist,i)) {
        if (specificMultiEdges.find(ime) != specificMultiEdges.end()) {
          cout << " > ist: " << ist << " touching ime: " << ime << " nodes: " << spost[ist][i] << " " << spost[ist][(i+1)%3] <<
            " coords: " << COUT_VEC(xsp[spost[ist][i]]) << " " << COUT_VEC(xsp[spost[ist][(i+1)%3]]) << endl;
          st_flag[ist] = 1;
        }
      }
    }
  }
}

void SimpleSurface::flagTrisTouchingOpenEdges() {
  cout << "SimpleSurface::flagTrisTouchingOpenEdges()" << endl;
  ensureTeost();
  st_flag.resize(nst);
  st_flag.setAll(0);
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      if (isEdgeOpen(ist,i)) {
        const double n[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
        cout << " > ist: " << ist << " edge: " << i << " is open edge, tris area: " << 0.5*MAG(n) << endl;
        st_flag[ist] = 1;
      }
    }
  }
}

void SimpleSurface::flagTrisTouchingFlaggedTris() {
  cout << "SimpleSurface::flagTrisTouchingFlaggedTris()" << endl;
  ensureTeost();
  int st_flag_max = 0;
  for (int ist = 0; ist < nst; ++ist)
    st_flag_max = max(st_flag_max,st_flag[ist]);
  if (st_flag_max == 0) {
    cout << "no tris are currently flagged" << endl;
    return;
  }
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] == 0) {
      FOR_I3 {
        int ist_nbr;
        if (getAlignedTriNbr(ist_nbr,ist,i)) {
          if ((st_flag[ist_nbr] > 0)&&(st_flag[ist_nbr] <= st_flag_max)) {
            st_flag[ist] = st_flag_max+1;
            break;
          }
        }
      }
    }
  }
}

void SimpleSurface::flagTrisTouchingFlaggedTrisInSameZone() {
  cout << "SimpleSurface::flagTrisTouchingFlaggedTrisInSameZone()" << endl;
  ensureTeost();
  int st_flag_max = 0;
  for (int ist = 0; ist < nst; ++ist) st_flag_max = max(st_flag_max,st_flag[ist]);
  if (st_flag_max == 0) {
    cout << "no tris are currently flagged" << endl;
    return;
  }
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] == 0) {
      FOR_I3 {
        int ist_nbr;
        if (getAlignedTriNbr(ist_nbr,ist,i)) {
          if ((st_flag[ist_nbr] > 0)&&(st_flag[ist_nbr] <= st_flag_max)&&(znost[ist]==znost[ist_nbr])) {
            st_flag[ist] = st_flag_max+1;
            break;
          }
        }
      }
    }
  }
}

void SimpleSurface::flagTrisWithAreaLessThan(const double area) {
  cout << "SimpleSurface::flagTrisWithAreaLessThan(): " << area << endl;
  st_flag.resize(nst);
  st_flag.setAll(0);
  for (int ist = 0; ist < nst; ++ist) {
    const double n[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    if (0.5*MAG(n) < area) {
      st_flag[ist] = 1;
    }
  }
}

void SimpleSurface::flagTrisWithSubzoneAreaLessThan(const double area) {
  cout << "SimpleSurface::flagTrisWithSubzoneAreaLessThan(): " << area << endl;
  ensureSubzoneData();
  st_flag.resize(nst);
  st_flag.setAll(0);
  for (int ist = 0; ist < nst; ++ist) {
    const int isz = szost[ist];
    assert((isz >= 0)&&(isz < nsz));
    if (subzoneDataVec[isz].area < area) {
      st_flag[ist] = 1;
    }
  }
}

void SimpleSurface::flagNmeTris(const vector<int>& zone_indices,const int nme_bits,const bool b_add) {
  // always allow user to downselect which zones
  cout << "SimpleSurface::flagNmeTrisInZones()" << endl;

  // don't simply ensure since you don't know how the user may want to populate the nme data
  if (!hasNonManifoldData()) {
    WUI(WARN,"please execute RUN_DIAGNOSTICS to identify non-manifold elements before execution; skipping");
    return;
  }

  if (b_add) {
    if (st_flag.size() != nst) {
      WUI(WARN,"ADD cannot be used on the first FLAG_TRIS* command; skipping");
      return;
    }
  }
  else {
    st_flag.resize(nst);
    st_flag.setAll(0);
  }



}

void SimpleSurface::flagTrisInZones(const vector<int>& zone_indices,const bool b_add) {

  cout << "SimpleSurface::flagTrisInZones()" << endl;

  if (b_add) {
    if (st_flag.size() != nst) {
      WUI(WARN,"ADD cannot be used on the first FLAG_TRIS* command; skipping");
      return;
    }
  }
  else {
    st_flag.resize(nst);
    st_flag.setAll(0);
  }

  int nzones = zoneVec.size();
  int * zn_flag = new int[nzones];
  for (int izn = 0; izn < nzones; ++izn)
    zn_flag[izn] = 0;

  for (int ii = 0; ii < zone_indices.size(); ++ii) {
    const int izn = zone_indices[ii];
    assert((izn >= 0)&&(izn < nzones));
    zn_flag[izn] = 1;
  }

  int nflagged_new = 0;
  int nflagged_old = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] == 0) {
      const int izn = znost[ist];
      assert((izn >= 0)&&(izn < nzones));
      if (zn_flag[izn]) {
        st_flag[ist] = 1;
        ++nflagged_new;
      }
    }
    else {
      ++nflagged_old;
    }
  }

  delete[] zn_flag;

  if (b_add) {
    WUI(INFO,"flagged " << nflagged_new << " additional tris. Total flagged tris: " << nflagged_old+nflagged_new);
  }
  else {
    WUI(INFO,"flagged " << nflagged_new << " tris.");
  }

}

void SimpleSurface::flagTrisTouchingZones(const vector<int>& zone_indices,const bool b_add) {

  cout << "SimpleSurface::flagTrisTouchingZones()" << endl;

  if (b_add) {
    if (st_flag.size() != nst) {
      WUI(WARN,"ADD cannot be used on the first FLAG_TRIS* command");
      return;
    }
  }
  else {
    st_flag.resize(nst);
    st_flag.setAll(0);
  }

  int nzones = zoneVec.size();
  int * zn_flag = new int[nzones];
  for (int izn = 0; izn < nzones; ++izn)
    zn_flag[izn] = 0;

  for (int ii = 0; ii < zone_indices.size(); ++ii) {
    const int izn = zone_indices[ii];
    assert((izn >= 0)&&(izn < nzones));
    zn_flag[izn] = 1;
  }

  sp_flag.resize(nsp);
  sp_flag.setAll(0);

  for (int ist = 0; ist < nst; ++ist) {
    const int izn = znost[ist];
    assert((izn >= 0)&&(izn < nzones));
    if (zn_flag[izn]) {
      FOR_I3 {
        const int isp = spost[ist][i];
        sp_flag[isp] = 1;
      }
    }
  }

  delete[] zn_flag;

  int nflagged_new = 0;
  int nflagged_old = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] == 0) {
      if (sp_flag[spost[ist][0]] || sp_flag[spost[ist][1]] || sp_flag[spost[ist][2]]) {
        st_flag[ist] = 1;
        ++nflagged_new;
      }
    }
    else {
      ++nflagged_old;
    }
  }

  if (b_add) {
    WUI(INFO,"flagged " << nflagged_new << " additional tris. Total flagged tris: " << nflagged_old+nflagged_new);
  }
  else {
    WUI(INFO,"flagged " << nflagged_new << " tris.");
  }

}

void SimpleSurface::countFlaggedTris() {

  if (st_flag.size() != nst) {
    WUI(WARN,"tri flag is not set");
    return;
  }

  int n = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist]) {
      ++n;
    }
  }

  WUI(INFO,n << " out of " << nst << " tris are flagged");

}

void SimpleSurface::invertFlaggedTris() {

  if (st_flag.size() != nst) {
    WUI(WARN,"tri flag is not set");
    return;
  }

  int n = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist]) {
      st_flag[ist] = 0;
    }
    else {
      st_flag[ist] = 1;
      ++n;
    }
  }

  WUI(INFO,"tri flagging inverted. Now " << n << " out of " << nst << " tris are flagged");

}
