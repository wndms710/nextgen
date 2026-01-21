#include "SimpleSurface.hpp"
#include "../core/common/GeomUtils.hpp"
#include "WebUI.hpp"

uint SimpleSurface::packEdge(const uint ist,const uchar i) const {
  assert( (ist<(uint(1)<<30)) && (i<3) );
  // ist (30 bits) | i (2 bits)
  uint edge = ((ist << 2) | uint(i));
  return edge;
}

uint SimpleSurface::packEdge(const int ist,const int i) const {
  return packEdge(uint(ist),uchar(i));
}

void SimpleSurface::unpackEdge(uint& ist,uchar& i,const uint edge) const {
  // ist (30 bits) | i (2 bits)
  ist = (edge >> 2) & MASK_30BITS;
  i = edge & 3;  // mask
}

void SimpleSurface::buildOpenEdgeGroups(const double crease_angle) {
  cout << "SimpleSurface::buildOpenEdgeGroups()" << endl;

  n_open_edge_groups = buildHalfEdgeGroups(openEdgeGroupDataVec, oe_to_group,&SimpleSurface::isEdgeOpen);

  // open edges store index data in teost, put there as well for now
  for (map<uint,int>::const_iterator it=oe_to_group.begin(); it!=oe_to_group.end(); ++it) {
    uint ist; uchar i;
    unpackEdge(ist,i,it->first);
    setOpenEdgeGroup(it->second,ist,i);
  }

  if (n_open_edge_groups) {
    COUT1(" > number of open edge groups: " << n_open_edge_groups << " (" << oe_to_group.size() << " half-edges)");
  }
  else {
    COUT1(" > no open edges were found");
  }
  b_open_edge_groups = true;

}

int SimpleSurface::getNOpenEdgeGroups() const {
  return n_open_edge_groups;
}

void SimpleSurface::ensureOpenEdgeGroups(const double crease_angle) {

  if (!b_open_edge_groups) {

    buildOpenEdgeGroups(crease_angle);
    assert(b_open_edge_groups);

  }

}

void SimpleSurface::clearOpenEdgeGroups() {

  b_open_edge_groups = false;

  n_open_edge_groups = 0; // not really necessary, but prevents possible confusion

  // openEdgeGroupDataVec must be cleared as well
  openEdgeGroupDataVec.clear();
  oe_to_group.clear();

}

int SimpleSurface::reportOpenEdgeGroupData(const bool detailed) {

  ensureOpenEdgeGroups();

  // TODO: replace this with WUI(MESSAGE... or similar construction when ready...
  UIMessage_ msg(MESSAGE);
  msg.addLine("");
  msg.addLine(stringstream().flush() << " ----- open edge group information ----- ");
  int n_open_edges = 0;
  if (openEdgeGroupDataVec.size()) {

    // for now only open edge loops go in here, so assume reporting open edge loops.
    // when edgeloop vec gets used for other edges we will need to revisit this reporting
    int group = 0;
    int n_loops = 0;
    for (vector<HalfEdgeGroupData>::iterator it=openEdgeGroupDataVec.begin(); it!=openEdgeGroupDataVec.end(); ++it) {
      // we a are a loop if either we have no end points or they match
      if (!it->end_ed[0].first || !it->end_ed[1].first) {
        ++n_loops;
      }
      else {
        uint ist0; uchar i0;
        unpackEdge(ist0,i0,it->end_ed[0].second);
        const int isp0 = spost[ist0][i0];
        uint ist1; uchar i1;
        unpackEdge(ist1,i1,it->end_ed[1].second);
        const int isp1 = spost[ist1][i1];
        if (isp0 == isp1) ++n_loops;
      }

      if (detailed) {
        msg.addLine("");
        msg.addLine(stringstream().flush() << "-- openEdgeGroup[" << group++ << "]:");
        msg.addLine(stringstream().flush() << " > number of edges: " << it->ned);
        msg.addLine(stringstream().flush() << " > length of group: " << it->length);
        msg.addLine(stringstream().flush() << " > centroid: " << COUT_VEC(it->xc));
        msg.addLine(stringstream().flush() << " > group end half-edges: " << ((it->end_ed[0].first) ? it->end_ed[0].second:-1) << " " << ((it->end_ed[1].first) ? it->end_ed[1].second:-1));
        for (set<int>::const_iterator it2 = it->zones.begin(); it2 != it->zones.end(); ++it2)
          msg.addLine(stringstream().flush() << " > zone: " << *it2);
      }
      n_open_edges += it->ned;

    }

    msg.addLine("");
    msg.addLine(stringstream().flush() << " > " << n_open_edges << " open edges were found in this surface");
    msg.addLine(stringstream().flush() << " > " << n_open_edge_groups << " open edge groups (connected-links) were constructed, comprising " << n_loops << " loops");
  }
  else {
    msg.addLine("");
    msg.addLine(stringstream().flush() << "There are no open edges on this surface");
  }
  WebUI::webUIOutput.add_(msg);
  return n_open_edges;
}

bool groupsD2SortCrit(pair<int,double> g0, pair<int,double> g1) {
  return (g0.second < g1.second);
}

void SimpleSurface::findNestedHalfEdgeLoops(vector<pair<int,int> >& gr_pairs,const double tol,const vector<int>& group_indices,const vector<HalfEdgeGroupData>& halfEdgeGroupDataVec) const {
  // use ADT to find pairs
  const double tol2 = tol*tol;
  const int n_groups = group_indices.size();

  // create ADT
  const double start_time = MPI_Wtime();

  // build Adts for all group centroids
  cout << " > building half-edge groups Adt: ";
  cout.flush();
  double (*bbmin)[3] = new double[n_groups][3];
  double (*bbmax)[3] = new double[n_groups][3];
  for (vector<int>::const_iterator it=group_indices.begin(); it!=group_indices.end(); ++it) {
    const int igr = it-group_indices.begin();
    FOR_J3 bbmin[igr][j] = halfEdgeGroupDataVec[*it].xc[j] - tol;
    FOR_J3 bbmax[igr][j] = halfEdgeGroupDataVec[*it].xc[j] + tol;
  }
  Adt<double> grAdt(n_groups,bbmin,bbmax);
  DELETE(bbmin);
  DELETE(bbmax);
  const double adt_time = MPI_Wtime()-start_time;
  cout << " finished in " << adt_time << "s" << endl;
  cout.flush();

  // loop through groups and find matches; only processing un-set candidates
  vector<int> candidates;
  vector<int>::iterator c_it;
  IntFlag gr_flag(halfEdgeGroupDataVec.size());
  gr_flag.setAll(0);

  for (vector<int>::const_iterator grit=group_indices.begin(); grit!=group_indices.end(); ++grit) {
    const int igroup = *grit;
    if (gr_flag[igroup] == 0) {
      // hasn't been matched yet

      // create list of centroids within tolerance
      // then match candidates with each other based on a sorted mean_d
      // leverage the fact that a set is ordered already
      vector<pair<double,int> > d2_to_gr;
      d2_to_gr.push_back(pair<double,int> (halfEdgeGroupDataVec[igroup].mean_d,igroup));

      candidates.clear();
      grAdt.buildListForPoint(candidates,halfEdgeGroupDataVec[igroup].xc);
      for (c_it=candidates.begin(); c_it!=candidates.end(); ++c_it) {
        const int c_igroup = group_indices[*c_it];
        if ((c_igroup!=igroup) && (gr_flag[c_igroup]==0)) {
          // for all remaining candidates who are not me
          const double dist2 = DIST2(halfEdgeGroupDataVec[igroup].xc,halfEdgeGroupDataVec[c_igroup].xc);
          if (dist2 <= tol2) {
            // valid candidate found; add to list of potential nesters
            d2_to_gr.push_back(pair<double,int> (halfEdgeGroupDataVec[c_igroup].mean_d,c_igroup));
          }
        }
      }

      if (int(d2_to_gr.size()) > 1) {
        sort(d2_to_gr.begin(),d2_to_gr.end());  // default sorts ascending on first element

        // match pairs based on mean_d, starting with smallest
        //TODO loop properly, can't do with iterator easily...?
        for (int igr=0,ngr=d2_to_gr.size(); igr < (ngr-1); igr+=2) {
          const int igr0 = d2_to_gr[igr].second;
          const int igr1 = d2_to_gr[igr+1].second;
          gr_pairs.push_back(pair<int,int>(igr0,igr1));
          gr_flag[igr0] = 1;
          gr_flag[igr1] = 1;
        }
      }

      // even if no valid candidates were found, remove from search
      gr_flag[igroup] = 1;
    }
  }
  candidates.clear();

  const int n_nested_groups = gr_pairs.size();

  COUT1(" > found " << n_nested_groups << " pairs of nested loops");
}

void SimpleSurface::buildDynamicEdgeGroups() {
  COUT2("SimpleSurface::buildDynamicEdgeGroups()");
  assert(!b_eoi);
  string eType = "unknown";
  switch (eoi_type) {
    // case METAL_BOUNDARY:
    //   buildMetalZoneEdgeGroups();
    //   eType = "metal boundaries";
    //   break;
    // case FLUID_BOUNDARY:
    //   buildFluidZoneEdgeGroups();
    //   eType = "fluid boundaries";
    //   break;
    case ZONE_BOUNDARY:
      buildZoneEdgeGroups();
      eType = "zone boundaries";
      break;
    case FEATURE:
      buildFeatureEdgeGroups();
      eType = "feature edges";
      break;
    case SUBZONE_BOUNDARY:
      buildSubzoneEdgeGroups();
      eType = "subzone boundaries";
      break;
    case MULTI:
      multiEdgesToDynamicEdges();
      eType = "multi-tri (3+) edges";
      break;
    case OPEN:
      buildOpenEdgeGroups2();
      eType = "open edges";
      break;
    case NME:
      ensureNonManifoldData();
      buildNmeEdgeGroups();
      eType = "non-manifold tri edges";
      break;
    case SELECTED_BOUNDARY:
      buildSelectedEdgeGroups();
      eType = "selection boundary edges";
      break;
    default:
      clearDynamicEdgeGroups();
  }
  COUT2(" > edge group type: " << eType);
}

void SimpleSurface::ensureDynamicEdgeGroups() {
  // simply rebuild with current settings
  if (!b_eoi) {
    buildDynamicEdgeGroups();
  }
  else {
    clearDynamicEdgeGroups();
    buildDynamicEdgeGroups();
  }
}

void SimpleSurface::ensureDynamicEdgeGroups(const DYNAMIC_EDGE_TYPE requested_type) {

  if (!b_eoi) {
    eoi_type = requested_type;
    buildDynamicEdgeGroups();
  }
  else if (requested_type != eoi_type) {
    eoi_type = requested_type;
    clearDynamicEdgeGroups();
    buildDynamicEdgeGroups();
  }
}

void SimpleSurface::clearDynamicEdgeGroups() {
  COUT2("SimpleSurface::clearDynamicEdgeGroups()");
  eoi_to_group.clear();
  eoiGroupDataVec.clear();
  b_eoi = false;
}

void SimpleSurface::buildZoneEdgeGroups() {
  int n_edge_groups = buildHalfEdgeGroups(eoiGroupDataVec, eoi_to_group,&SimpleSurface::isEdgeZoneBoundary);

  if (n_edge_groups) {
    COUT1(" > found " << n_edge_groups << " half-edge groups on zone boundaries (" << eoi_to_group.size() << " half-edges)");
    b_eoi = true;
  }
  else {
    COUT1("no zone boundary edge groups were found");
  }
}

void SimpleSurface::buildSubzoneEdgeGroups() {
  int n_edge_groups = buildHalfEdgeGroups(eoiGroupDataVec, eoi_to_group,&SimpleSurface::isEdgeSubzoneBoundary);

  if (n_edge_groups) {
    COUT1(" > found " << n_edge_groups << " half-edge groups on subzone boundaries (" << eoi_to_group.size() << " half-edges)");
    b_eoi = true;
  }
  else {
    COUT1("no subzone boundary edge groups were found");
  }
}

void SimpleSurface::buildSelectedEdgeGroups() {
  sz_flag.resize(nsz);
  sz_flag.setAll(0);
  for (vector<int>::iterator it=selectedSubzoneVec.begin(); it!= selectedSubzoneVec.end(); ++it) {
    if (*it>=0 && *it < nsz) {
      sz_flag[*it] = 1;
    }
  }
  st_flag.resize(nst);
  st_flag.setAll(0);
  for (int ist=0; ist<nst; ++ist) {
    if (sz_flag[szost[ist]]) st_flag[ist] = 1;
  }
  int n_edge_groups = buildHalfEdgeGroups(eoiGroupDataVec, eoi_to_group,&SimpleSurface::isEdgeSelectionBoundary);

  if (n_edge_groups) {
    COUT1(" > found " << n_edge_groups << " half-edge groups on subzone boundaries (" << eoi_to_group.size() << " half-edges)");
    b_eoi = true;
  }
  else {
    COUT1("no subzone boundary edge groups were found");
  }
}

void SimpleSurface::buildFeatureEdgeGroups() {
  int n_edge_groups = buildHalfEdgeGroups(eoiGroupDataVec, eoi_to_group,&SimpleSurface::isEdgeFeature);

  if (n_edge_groups) {
    COUT1(" > found " << n_edge_groups << " feature half-edge groups (" << eoi_to_group.size() << " half-edges)");
    b_eoi = true;
  }
  else {
    COUT1("no feature edge groups were found");
  }
}

void SimpleSurface::buildNmeEdgeGroups() {
  nmTriData.dump();

  int n_edge_groups = buildHalfEdgeGroups(eoiGroupDataVec, eoi_to_group,&SimpleSurface::isEdgeNmeBoundary);

  if (n_edge_groups) {
    COUT1(" > found " << n_edge_groups << " non-manifold half-edge groups (" << eoi_to_group.size() << " half-edges)");
    b_eoi = true;
  }
  else {
    COUT1("no feature edge groups were found");
  }
}

void SimpleSurface::buildOpenEdgeGroups2() {
  int n_edge_groups = buildHalfEdgeGroups(eoiGroupDataVec, eoi_to_group,&SimpleSurface::isEdgeOpen);

  if (n_edge_groups) {
    COUT1(" > found " << n_edge_groups << " open-edge groups (" << eoi_to_group.size() << " half-edges)");
    b_eoi = true;
  }
  else {
    COUT1("no open edge groups were found");
  }
}

void SimpleSurface::selectEdgesInWindow(const double window[4][3],const bool b_strictly_inside,const bool open) {

  const double dx01[3] = DIFF(window[1],window[0]);
  const double dx12[3] = DIFF(window[2],window[1]);
  const double n[3] = CROSS_PRODUCT(dx01,dx12);

  bool* b_inside_sp = new bool[nsp];
  FOR_ISP b_inside_sp[isp] = true;
  int ii0 = 3;
  for (int ii = 0; ii < 4; ++ii) {

    // get plane orthogonal to this edge and n
    const double dx_ii[3] = DIFF(window[ii],window[ii0]);
    const double n_ii[3] = CROSS_PRODUCT(dx_ii,n);

    for (int isp = 0; isp < nsp; ++isp) {
      if (b_inside_sp[isp]) {
        const double dxp[3] = DIFF(xsp[isp],window[ii]);
        const double sd = DOT_PRODUCT(n_ii,dxp);
        if (sd > 0.0)
          b_inside_sp[isp] = false;
      }
    }

    ii0 = ii;
  }

  bool* b_hide_sz = new bool[nsz];
  for (int isz = 0; isz < nsz; ++isz) b_hide_sz[isz] = false;
  for (int ii = 0, lim = hiddenSubzoneVec.size(); ii < lim; ++ii) {
    const int isz = hiddenSubzoneVec[ii];
    if ((isz >= 0)&&(isz < nsz)) {
      b_hide_sz[isz] = true;
    }
    else {
      CWARN(" > hidden subzone: " << isz << " not in subzone range, [0," << nsz << "); skipping..." )
    }
  }

  bool* b_inside_gr = NULL;
  if (open) {
    const int ngr = openEdgeGroupDataVec.size();
    //cout << ngr << " " << oe_to_group.size() << endl;
    b_inside_gr = new bool[ngr];
    if (b_strictly_inside) {
      for (int igr = 0; igr < ngr; ++igr) b_inside_gr[igr] = true;
      for (map<uint,int>::const_iterator it=oe_to_group.begin(); it!=oe_to_group.end(); ++it) {
        const int igr = it->second;
        if (b_inside_gr[igr]) {
          uint ist; uchar i;
          unpackEdge(ist,i,it->first);
          if (b_hide_sz[szost[ist]]||(!(b_inside_sp[spost[ist][i]]&&(b_inside_sp[spost[ist][(i+1)%3]]))))
            b_inside_gr[igr] = false;
        }
      }
    }
    else {
      for (int igr = 0; igr < ngr; ++igr) b_inside_gr[igr] = false;
      for (map<uint,int>::const_iterator it=oe_to_group.begin(); it!=oe_to_group.end(); ++it) {
        const int igr = it->second;
        if (!b_inside_gr[igr]) {
          uint ist; uchar i;
          unpackEdge(ist,i,it->first);
          if ((!b_hide_sz[szost[ist]])&&(b_inside_sp[spost[ist][i]]||b_inside_sp[spost[ist][(i+1)%3]]))
            b_inside_gr[igr] = true;
        }
      }
    }
    for (int igr = 0; igr < ngr; ++igr) {
      if (b_inside_gr[igr]) {
        selectedSubzoneVec.push_back(igr+OPEN_E_SZ_MIN);
        //cout << igr+OPEN_E_SZ_MIN << endl;
      }
    }
  }
  else {
    const int ngr = eoiGroupDataVec.size();
    //cout << ngr << " " << eoi_to_group.size() << endl;
    b_inside_gr = new bool[ngr];
    if (b_strictly_inside) {
      for (int igr = 0; igr < ngr; ++igr) b_inside_gr[igr] = true;
      for (map<uint,int>::const_iterator it=eoi_to_group.begin(); it!=eoi_to_group.end(); ++it) {
        const int igr = it->second;
        if (b_inside_gr[igr]) {
          uint ist; uchar i;
          unpackEdge(ist,i,it->first);
          if (b_hide_sz[szost[ist]]||(!(b_inside_sp[spost[ist][i]]&&(b_inside_sp[spost[ist][(i+1)%3]]))))
            b_inside_gr[igr] = false;
        }
      }
    }
    else {
      for (int igr = 0; igr < ngr; ++igr) b_inside_gr[igr] = false;
      for (map<uint,int>::const_iterator it=eoi_to_group.begin(); it!=eoi_to_group.end(); ++it) {
        const int igr = it->second;
        if (!b_inside_gr[igr]) {
          uint ist; uchar i;
          unpackEdge(ist,i,it->first);
          if ((!b_hide_sz[szost[ist]])&&(b_inside_sp[spost[ist][i]]||b_inside_sp[spost[ist][(i+1)%3]]))
            b_inside_gr[igr] = true;
        }
      }
    }
    for (int igr = 0; igr < ngr; ++igr) {
      if (b_inside_gr[igr]) {
        selectedSubzoneVec.push_back(igr+DYN_E_SZ_MIN);
        //cout << igr+DYN_E_SZ_MIN << endl;
      }
    }
  }
  delete[] b_inside_sp;
  delete[] b_inside_gr;
  delete[] b_hide_sz;

}

void SimpleSurface::orderEdgesInGroup(vector<uint>& edgeVec,vector<int>& group_indices,const map<uint,int>& edge_to_group,const vector<HalfEdgeGroupData>& heGroupDataVec) {

  if (!edge_to_group.size()) {
    CWARN("dynamic edge groups have not been computed; skipping...");
    return;
  }

  map<int,int> sp_flag_m;  // valence count for group start/end nodes only
  multimap<int,int> head2group;

  // get the ordered set of edges for a single group
  // this will provide start/end edges, and use similar grouping logic to determine the ordering of groups tpass into the edgeVec
  vector<vector<uint> > group_edges;
  for (int igr=0, ngr=group_indices.size(); igr<ngr; ++igr) {

    vector<uint> _edgeVec;
    orderEdgesInGroup(_edgeVec,group_indices[igr],edge_to_group,heGroupDataVec);

    // process head/tail for the group (don't currently trust the data stored in eoiEdgeDataVec.end)
    uint ist0; uchar i0;
    unpackEdge(ist0,i0,_edgeVec[0]);
    const int isp0 = spost[ist0][i0];
    uint ist1; uchar i1;
    unpackEdge(ist1,i1,_edgeVec.back());
    const int isp1 = spost[ist1][(i1+1)%3];
    head2group.insert(pair<int,int>(isp0,igr));

    if (isp0 != isp1) {
      // chain, not a loop
      //add start/end nodes for chain-stacking
      map<int,int>::iterator sp_it=sp_flag_m.find(isp0);
      if (sp_it!=sp_flag_m.end()) sp_it->second += 1;
      else sp_flag_m.insert(pair<int,int> (isp0,1));
      sp_it=sp_flag_m.find(isp1);
      if (sp_it!=sp_flag_m.end()) sp_it->second += 1;
      else sp_flag_m.insert(pair<int,int> (isp1,1));
    }

    group_edges.push_back(_edgeVec);
  }

  // now order based on edge group connectivity
  // use sp_flag_m to ensure we always start marching from valence = 1 or >2 nodes
  vector<int> starting_isps;
  for (map<int,int>::iterator sp_it=sp_flag_m.begin(); sp_it!=sp_flag_m.end(); ++sp_it) {
    if (sp_it->second != 2) {
      // this is a node of interest
      // see if any edges leave from here (head is here...)
      if (head2group.count(sp_it->first)) starting_isps.push_back(sp_it->first);
    }
  }
  sp_flag_m.clear();

  int isp1;  // group tail
  multimap<int,int>::iterator h2g_it;
  vector<int> groups_ordered;
  while (!head2group.empty()) {

    // identify starting node and order groups; set isp1 for the search
    if (starting_isps.size()) {
      isp1 = starting_isps.back();
      starting_isps.pop_back();
    }
    else {
      // prune first group
      h2g_it = head2group.begin();
      const int igr = h2g_it->second;
      uint ist; uchar i;
      unpackEdge(ist,i,group_edges[igr].back());
      isp1 = spost[ist][(i+1)%3];
      head2group.erase(h2g_it);
      groups_ordered.push_back(igr);
    }

    // search groups in chain; random new ordering at valence > 2
    while (true) {
      // find all edges in this chain
      if (head2group.count(isp1) == 1) {
        h2g_it = head2group.find(isp1);
        const int igr = h2g_it->second;
        uint ist; uchar i;
        unpackEdge(ist,i,group_edges[igr].back());
        isp1 = spost[ist][(i+1)%3];
        head2group.erase(h2g_it);
        groups_ordered.push_back(igr);
      }
      else {
        // either 0 (terminal) or multiple candidates (intersection) found; either way we start a new group
        break;
      }

    }
  }
  assert(starting_isps.size() == 0);
  assert(groups_ordered.size() == group_indices.size());

  // push into single edge vector for all groups
  for (vector<int>::const_iterator it=groups_ordered.begin(); it!=groups_ordered.end(); ++it) {
    edgeVec.insert(edgeVec.end(), group_edges[*it].begin(), group_edges[*it].end());
  }
}

void SimpleSurface::orderHalfEdgesInVec(vector<uint>& edgeVec) {
  // assumes all walkable loops for now... (otherwise need directed graph stuff)
  // pass in half edge indices, will get them sorted in a walkable chain...
  const int n_original_edges = edgeVec.size();

  // this counts edges and sets sp_flag for walking
  sp_flag.setLength(nsp);
  sp_flag.setAll(-1);
  map<int,int> sp_flag_m;  // valence count for touched nodes only
  multimap<int,uint> head2edge;

  for (vector<uint>::const_iterator it=edgeVec.begin(); it!=edgeVec.end(); ++it) {
    uint ist; uchar i;
    unpackEdge(ist,i,*it);

    // construct a map from head node to the edge
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    head2edge.insert(pair<int,uint>(isp0,*it));

    // update sp_flag for group-local valence
    map<int,int>::iterator sp_it=sp_flag_m.find(isp0);
    if (sp_it!=sp_flag_m.end()) {
      // already created, so increase valence count
      sp_it->second += 1;
    }
    else {
      // first time this isp is hit, so create a count of 1
      sp_flag_m.insert(pair<int,int> (isp0,1));
    }

    sp_it=sp_flag_m.find(isp1);
    if (sp_it!=sp_flag_m.end()) sp_it->second += 1;  // add to valence
    else sp_flag_m.insert(pair<int,int> (isp1,1));  // first time, create a count of 1
  }

  // use sp_flag_m to ensure we always start marching from valence = 1 or >2 nodes
  vector<int> starting_isps;
  for (map<int,int>::const_iterator sp_it=sp_flag_m.begin(); sp_it!=sp_flag_m.end(); ++sp_it) {
    if (sp_it->second != 2) {
      // this is a node of interest
      // see if any edges leave from here (head is here...)
      if (head2edge.count(sp_it->first)) starting_isps.push_back(sp_it->first);
    }
  }
  sp_flag_m.clear();

  // edgemap in place, now connect loops and repopulate
  edgeVec.clear();
  int ned = 0;
  int isp1;
  uint ied;
  multimap<int,uint>::iterator h2t_it;
  while (!head2edge.empty()) {

    // identify starting node and march edges; effectively set isp1 for the search
    if (starting_isps.size()) {
      // if terminal nodes existed, start from here
      isp1 = starting_isps.back();
      starting_isps.pop_back();
    }
    else {
      // otherwise these are loops, so any edge is sufficient to start
      // prune first edge
      h2t_it = head2edge.begin();
      ied = h2t_it->second;
      uint ist; uchar i;
      unpackEdge(ist,i,ied);
      isp1 = spost[ist][(i+1)%3];
      head2edge.erase(h2t_it);
      edgeVec.push_back(ied);
      ++ned;
    }

    // search edges in chain; new group at valence > 2
    while (true) {
      // find all edges connected to isp1
      if (head2edge.count(isp1) == 1) {
        h2t_it = head2edge.find(isp1);
        ied = h2t_it->second;
        uint ist; uchar i;
        unpackEdge(ist,i,ied);
        isp1 = spost[ist][(i+1)%3];
        head2edge.erase(h2t_it);
        edgeVec.push_back(ied);
        ++ned;
      }
      else {
        // either 0 (terminal) or multiple candidates (intersection) found; either way we start a new group
        break;
      }

    }
  }
  assert(starting_isps.size() == 0);
  assert(n_original_edges == int(edgeVec.size()));
}

void SimpleSurface::orderEdgesInGroup(vector<uint>& edgeVec,const int igroup,const map<uint,int>& edge_to_group,const vector<HalfEdgeGroupData>& heGroupDataVec) {
  // assumes all walkable loops for now... (otherwise need directed graph stuff)
  if (edge_to_group.empty()) {
    CWARN("Edges-of-interest have not been grouped; skipping");
    return;
  }

  if (igroup < 0 || igroup >= int(heGroupDataVec.size())) {
    CWARN("invalid edges-of-interest group " << igroup << "; skipping");
    return;
  }

  // this counts edges and sets sp_flag for walking
  sp_flag.setLength(nsp);
  sp_flag.setAll(-1);
  map<int,int> sp_flag_m;  // valence count for touched nodes only
  multimap<int,uint> head2edge;

  for (map<uint,int>::const_iterator it=edge_to_group.begin(); it!=edge_to_group.end(); ++it) {
    if (it->second == igroup) {
      uint ist; uchar i;
      unpackEdge(ist,i,it->first);

      // construct a map from head node to the edge
      const int isp0 = spost[ist][i];
      const int isp1 = spost[ist][(i+1)%3];
      head2edge.insert(pair<int,uint>(isp0,it->first));

      // update sp_flag for group-local valence
      map<int,int>::iterator sp_it=sp_flag_m.find(isp0);
      if (sp_it!=sp_flag_m.end()) {
        // already created, so increase valence count
        sp_it->second += 1;
      }
      else {
        // first time this isp is hit, so create a count of 1
        sp_flag_m.insert(pair<int,int> (isp0,1));
      }

      sp_it=sp_flag_m.find(isp1);
      if (sp_it!=sp_flag_m.end()) sp_it->second += 1;  // add to valence
      else sp_flag_m.insert(pair<int,int> (isp1,1));  // first time, create a count of 1
    }
  }

  // use sp_flag_m to ensure we always start marching from valence = 1 or >2 nodes
  vector<int> starting_isps;
  for (map<int,int>::const_iterator sp_it=sp_flag_m.begin(); sp_it!=sp_flag_m.end(); ++sp_it) {
    if (sp_it->second != 2) {
      // this is a node of interest
      // see if any edges leave from here (head is here...)
      if (head2edge.count(sp_it->first)) starting_isps.push_back(sp_it->first);
    }
  }
  sp_flag_m.clear();

  // edgemap in place, now connect loops
  int ned = 0;
  uint ied;
  int isp1;
  multimap<int,uint>::iterator h2t_it;
  while (!head2edge.empty()) {

    // identify starting node and march edges; effectively set isp1 for the search
    if (starting_isps.size()) {
      // if terminal nodes existed, start from here
      isp1 = starting_isps.back();
      starting_isps.pop_back();
    }
    else {
      // otherwise these are loops, so any edge is sufficient to start
      // prune first edge
      h2t_it = head2edge.begin();
      ied = h2t_it->second;
      uint ist; uchar i;
      unpackEdge(ist,i,ied);
      isp1 = spost[ist][(i+1)%3];
      head2edge.erase(h2t_it);
      edgeVec.push_back(ied);
      ++ned;
    }

    // search edges in chain; new group at valence > 2
    while (true) {
      // find all edges connected to isp1
      if (head2edge.count(isp1) == 1) {
        h2t_it = head2edge.find(isp1);
        ied = h2t_it->second;
        uint ist; uchar i;
        unpackEdge(ist,i,ied);
        isp1 = spost[ist][(i+1)%3];
        head2edge.erase(h2t_it);
        edgeVec.push_back(ied);
        ++ned;
      }
      else {
        // either 0 (terminal) or multiple candidates (intersection) found; either way we start a new group
        break;
      }

    }
  }
  assert(starting_isps.size() == 0);

  assert(heGroupDataVec[igroup].ned == ned);
}

bool SimpleSurface::getHalfEdgeGroup(int& igr,const uint ist,const uchar i) {
  map<uint,int>::iterator it;
  it = eoi_to_group.find(packEdge(ist,i));
  if (it != eoi_to_group.end()) {
    // element found
    igr = it->second;
    return true;
  }
  else {
    return false;
  }
}

int SimpleSurface::buildHalfEdgeGroups(vector<HalfEdgeGroupData>& halfEdgeGroupDataVec, map<uint,int>& he_to_gr, bool (SimpleSurface::*halfEdgeCriterion)(const int,const int) const) {
  // pass half edge criteria based on ist,i (half edge index)

  cout << "SimpleSurface::buildHalfEdgeGroups()" << endl;

  ensureTeost();

  sp_flag.setLength(nsp); // -1 unvisited, -2 intersection/crease, < -2 seed, >= 0 link
  sp_flag.setAll(-1);

  int *noe_sp = new int[nsp]; // counts number of half edges that touch this node (if negative the node is "special")
  FOR_ISP noe_sp[isp] = 0;
  FOR_IST {
    FOR_I3 {
      if ((this->*halfEdgeCriterion)(ist,i) || isEdgeOpen(ist,i)) {  // open edge required for correct valence...
        const int isp0 = spost[ist][i];
        const int isp1 = spost[ist][(i+1)%3];
        assert(isp0 != isp1);

        // isp0
        if (noe_sp[isp0] == 0) {
          noe_sp[isp0] = 1;
        }
        else if (noe_sp[isp0] == 1) {
          noe_sp[isp0] = 2;
        }
        // flag intersection
        else {
          assert(noe_sp[isp0] >= 2);
          ++noe_sp[isp0];
          sp_flag[isp0] = -2;
        }

        // isp1
        if (noe_sp[isp1] == 0) {
          noe_sp[isp1] = 1;
        }
        else if (noe_sp[isp1] == 1) {
          noe_sp[isp1] = 2;
        }
        // flag intersection
        else {
          assert(noe_sp[isp1] >= 2);
          ++noe_sp[isp1];
          sp_flag[isp1] = -2;
        }
      }
    }
  }

  int terminal_count = 0;
  int high_valence_count = 0;
  FOR_ISP {
    if (noe_sp[isp] == 1) ++terminal_count;
    if (sp_flag[isp] == -2) ++high_valence_count;
  }
  if (terminal_count) cout << "Warning: " << terminal_count << " terminal points were detected" << endl;
  if (high_valence_count) cout << "Warning: " << high_valence_count << " points with valence > 2 were detected" << endl;

  // at this pt noe_sp and sp_flag are set up correctly to indicate creases/intersections
  // now go through again and do internal nodes
  int count = 0;
  FOR_IST {
    FOR_I3 {
      if ((this->*halfEdgeCriterion)(ist,i)) {
        const int isp0 = spost[ist][i];
        const int isp1 = spost[ist][(i+1)%3];
        assert(isp0 != isp1);
        if ((sp_flag[isp0] == -1)&&(sp_flag[isp1] == -1)) {
          // first time we've visted either of these nodes...
          ++count;
          sp_flag[isp0] = -count-2; // -3, -4, etc...
          sp_flag[isp1] = isp0; // positive indicates link to value in isp0
        }
        else if (sp_flag[isp0] == -1) {
          sp_flag[isp0] = isp1;
        }
        else if (sp_flag[isp1] == -1) {
          sp_flag[isp1] = isp0;
        }
        else if ((sp_flag[isp0] != -2)&&(sp_flag[isp1] != -2)) {
          // both have been set, so traverse both stacks to their negative index...
          int isp0_ = isp0;
          int isp0p_ = isp0_;
          while (sp_flag[isp0_] >= 0) {
            isp0p_ = isp0_;
            isp0_ = sp_flag[isp0_];
          }
          assert(sp_flag[isp0_] < -1);
          int isp1_ = isp1;
          int isp1p_ = isp1_;
          while (sp_flag[isp1_] >= 0) {
            isp1p_ = isp1_;
            isp1_ = sp_flag[isp1_];
          }
          assert(sp_flag[isp1_] < -1);
          // now create a link between the loser and the winner (the lowest index)...
          if (sp_flag[isp0_] > sp_flag[isp1_]) {
            sp_flag[isp1_] = isp0p_;
          }
          else if (sp_flag[isp0_] < sp_flag[isp1_]) {
            sp_flag[isp0_] = isp1p_;
          }
          else if (sp_flag[isp0_] == sp_flag[isp1_]) {
            if (sp_flag[isp0_] == -2) sp_flag[isp1p_] = isp0p_;
          }
        }
      }
    }
  }
  int min_internal_flag = -count-2;

  map<const pair<int,int>, int> ed2grp; // first pt is crease/intersection, second is open edge nbr, it points to a grp id

  if (high_valence_count) {
    cout << " > processing high-valence edge groups...." << endl;
    // only traverse this IST loop when high-valence nodes are present
    //ed2grp sits outside because used later for standard edges too

    const int count0 = -min_internal_flag-2;

    map<pair<int,int>, int> szsz2grp;  // subzones across half edge : group
    multimap<int,pair<int,int> > grp2ed_mmap;  // for each group list the edges
    FOR_IST {
      FOR_I3 {
        if ((this->*halfEdgeCriterion)(ist,i)) {
          const int isp0 = spost[ist][i];
          const int isp1 = spost[ist][(i+1)%3];
          if (sp_flag[isp0] == -2) {
            int group;
            int ist_nbr,i_nbr,orient_nbr;
            const bool valid = getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i);
            if (valid && (ist_nbr > -1)) {
              const int sz_ist = szost[ist];
              const int sz_ist_nbr = szost[ist_nbr];
              map<pair<int,int>,int>::iterator it = szsz2grp.find(pair<int,int>(sz_ist,sz_ist_nbr));
              if (it != szsz2grp.end()) {
                // sz-pair already registered; use this group index
                group = -it->second-2;
              }
              else {
                ++count;
                group = -count-2;
                szsz2grp.insert(pair<pair<int,int>,int>(pair<int,int>(sz_ist,sz_ist_nbr),count));
              }
            }
            else {
              // open or multi_edge; push to new group
              ++count;
              group = -count-2;
            }
            grp2ed_mmap.insert(pair<int,pair<int,int> >((group-count0),pair<int,int>(isp0,isp1)));
          }
          else if (sp_flag[isp1] == -2) {
            int group;
            int ist_nbr,i_nbr,orient_nbr;
            const bool valid = getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i);
            if (valid && (ist_nbr > -1)) {
              const int sz_ist = szost[ist];
              const int sz_ist_nbr = szost[ist_nbr];
              map<pair<int,int>,int>::iterator it = szsz2grp.find(pair<int,int>(sz_ist_nbr,sz_ist));  // switch order because endge orientation backwards
              if (it != szsz2grp.end()) {
                // sz-pair already registered; use this group index
                group = -it->second-2;
              }
              else {
                ++count;
                group = -count-2;
                szsz2grp.insert(pair<pair<int,int>,int>(pair<int,int>(sz_ist_nbr,sz_ist),count));
              }
            }
            else {
              ++count;
              group = -count-2;
            }
            grp2ed_mmap.insert(pair<int,pair<int,int> >((group-count0),pair<int,int>(isp1,isp0)));
          }
        }
      }
    }
    szsz2grp.clear();

    // now process these groups into disjoint loops...
    const int n_szsz_groups0 = count - count0;
    COUT2(" > number of sz-pair-based groups: " << n_szsz_groups0);

    // within each group determine disjoint subgroups...
    count = count0;
    for (multimap<int,pair<int,int> >::iterator it=grp2ed_mmap.begin(); it!=grp2ed_mmap.end(); ) {
      // starting new group
      const int my_group = it->first;  // 0 indexed; not the -3 indexed group number
      map<int,int> sp_flag_m;  // valence count for touched nodes only
      multimap<int,int> head2tail;
      int isp0,isp1;

      do {
        // for all members of this group, perform following
        // create map to search edges based on head node
        isp0 = it->second.first;
        map<int,int>::iterator sp_it=sp_flag_m.find(isp0);
        if (sp_it!=sp_flag_m.end()) {
          // already created, so increase valence count
          sp_it->second += 1;
        }
        else {
          // first time this isp is hit, so create a count of 1
          sp_flag_m.insert(pair<int,int> (isp0,1));
        }

        isp1 = it->second.second;
        sp_it=sp_flag_m.find(isp1);
        if (sp_it!=sp_flag_m.end()) {
          // already created, so increase valence count
          sp_it->second += 1;
        }
        else {
          // first time this isp is hit, so create a count of 1
          sp_flag_m.insert(pair<int,int> (isp1,1));
        }

        // create edge map
        head2tail.insert(pair<int,int>(isp0,isp1));
        ++it;
      } while (it->first == my_group);

      // use sp_flag_m to ensure we always start marching from valence = 1 or >2 nodes
      vector<int> starting_isps;
      for (map<int,int>::iterator sp_it=sp_flag_m.begin(); sp_it!=sp_flag_m.end(); ++sp_it) {
        if (sp_it->second != 2) {
          // this is a node of interest
          // see if any edges leave from here (head is here...)
          if (head2tail.count(sp_it->first)) starting_isps.push_back(sp_it->first);
        }
      }
      sp_flag_m.clear();

      // edgemap in place, now connect loops
      multimap<int,int>::iterator h2t_it;
      while (!head2tail.empty()) {
        // starting new group
        ++count;
        const int _group = -count-2;

        // identify starting node and march edges; effectively set isp1 for the search
        if (starting_isps.size()) {
          isp1 = starting_isps.back();
          starting_isps.pop_back();
        }
        else {
          // prune first edge
          h2t_it = head2tail.begin();
          isp0 = h2t_it->first;
          isp1 = h2t_it->second;
          head2tail.erase(h2t_it);
          ed2grp.insert(pair< pair<int,int>, int>(pair<int,int>(isp0,isp1),_group));
        }

        // search edges in chain; new group at valence > 2
        while (true) {
          // find all edges in this group
          if (head2tail.count(isp1) == 1) {
            h2t_it = head2tail.find(isp1);
            isp0 = h2t_it->first;
            isp1 = h2t_it->second;
            head2tail.erase(h2t_it);
            ed2grp.insert(pair< pair<int,int>, int>(pair<int,int>(isp0,isp1),_group));
          }
          else {
            // either 0 (terminal) or multiple candidates (intersection) found; either way we start a new group
            break;
          }

        }
      }
      assert(starting_isps.size() == 0);
    }
    COUT2("    > decimated into " << (count - count0) << " szsz-disjoint edge groups");
  }

  // traverse the sp_flag setting everyone to unique values...

  FOR_ISP noe_sp[isp] = sp_flag[isp];
  FOR_ISP {
    if (sp_flag[isp] >= 0) {
      // connected to a root...
      int isp_ = isp;
      int ispp_ = isp_;
      while (sp_flag[isp_] >= 0) {
        assert(isp_ != sp_flag[isp_]);
        ispp_ = isp_;
        isp_ = sp_flag[isp_];
      }
      // regular
      if (sp_flag[isp_] < -2) {
        sp_flag[isp] = sp_flag[isp_];
      }
      // crease or intersection
      else {
        assert(sp_flag[isp_] == -2);
        map<const pair<int,int>,int>::iterator it = ed2grp.find(pair<int,int>(isp_,ispp_));
        assert(it != ed2grp.end());
        sp_flag[isp] = it->second;
      }
    }
  }

  // get unique half-edge group index

  // count the internal ones that remained (pure loops)
  int * heg_index = new int[count]; // know unique count <= count
  int n_half_edge_groups = 0;
  FOR_ISP {
    if ((noe_sp[isp] == sp_flag[isp])&&(sp_flag[isp] < -2)) {
      assert(sp_flag[isp] >= min_internal_flag);
      heg_index[-sp_flag[isp]-3] = n_half_edge_groups++; // map to unique index
    }
  }

  for (int igr=(-min_internal_flag-2); igr<count; ++igr) {
    heg_index[igr] = n_half_edge_groups++;  // grouping already handled by ed2grp above; so simple index increase per group
  }

  // one pass through the teost_data to store this index.
  // also build edge group data while we are here...

  halfEdgeGroupDataVec.resize(n_half_edge_groups);
  FOR_IST {
    FOR_I3 {
      if ((this->*halfEdgeCriterion)(ist,i)) {
        const int isp0 = spost[ist][i];
        const int isp1 = spost[ist][(i+1)%3];
        assert((sp_flag[isp0] <= -2)&&(sp_flag[isp1] <= -2));
        int igr = -1;
        if ((sp_flag[isp0] < -2)&&(sp_flag[isp1] < -2)) {
          if (sp_flag[isp0] != sp_flag[isp1]) {
            cout << "Warning: internal edge: " << isp0 << " " << isp1;
            cout << " doesn't have unique value: " << sp_flag[isp0] << " " << sp_flag[isp1] << endl;
            //problems = true;
          }
          igr = heg_index[-sp_flag[isp0]-3];
          he_to_gr.insert(pair<uint,int> (packEdge(ist,i),igr));
        }
        else if (sp_flag[isp0] < -2) {
          igr = heg_index[-sp_flag[isp0]-3];
          he_to_gr.insert(pair<uint,int> (packEdge(ist,i),igr));
          // add end point to open edge group data
          //assert(halfEdgeGroupDataVec[igr].end_ed[1] == -1);
          // const uint edge = packEdge(ist,(i+1)%3);
          const uint edge = packEdge(ist,i);  // last edge should just be this edge...
          if (halfEdgeGroupDataVec[igr].end_ed[1].first) {
            cout << "Warning: half-edge group: " << igr << " is overwriting end_ed[1] value ";
            cout << halfEdgeGroupDataVec[igr].end_ed[1].second << " with " << edge << endl;
            cout << "Ensure that your normals are aligned for your disjoint parts and try again." << endl;
            //problems = true;
          };
          halfEdgeGroupDataVec[igr].end_ed[1].first = true;
          halfEdgeGroupDataVec[igr].end_ed[1].second = edge;
        }
        else if (sp_flag[isp1] < -2) {
          igr = heg_index[-sp_flag[isp1]-3];
          he_to_gr.insert(pair<uint,int> (packEdge(ist,i),igr));

          // add end point to start edge group data
          const uint edge = packEdge(ist,i);
          if (halfEdgeGroupDataVec[igr].end_ed[0].first) {
            cout << "Warning: half-edge group: " << igr << " is overwriting end_ed[0] value ";
            cout << halfEdgeGroupDataVec[igr].end_ed[0].second << " with " << edge << endl;
            cout << "Ensure that your normals are aligned for your disjoint parts and try again." << endl;
            //problems = true;
          };
          halfEdgeGroupDataVec[igr].end_ed[0].first = true;
          halfEdgeGroupDataVec[igr].end_ed[0].second = edge;
        }
        else {
          // both are ">2 edge" : >-< or |_| or combination
          map<const pair<int,int>,int>::iterator it = ed2grp.find(pair<int,int>(isp0,isp1));
          assert(it != ed2grp.end());
          igr = heg_index[-it->second-3];
          he_to_gr.insert(pair<uint,int> (packEdge(ist,i),igr));

          // add end points to edge group data
          // cannot do edge based because can have multiple intersections...
          // but previous grouping should address this...so use nodes as sufficient check of loop or not
          if ((!halfEdgeGroupDataVec[igr].end_ed[0].first)&&(!halfEdgeGroupDataVec[igr].end_ed[1].first)) {
            halfEdgeGroupDataVec[igr].end_ed[0] = pair<bool,uint> (true,packEdge(ist,i) );
            halfEdgeGroupDataVec[igr].end_ed[1] = pair<bool,uint> (true,packEdge(ist,i) );
          }
          else {
            // after being set always take min edge at each position; should at least help indicate if connected loops...
            halfEdgeGroupDataVec[igr].end_ed[0].second = min(packEdge(ist,i),halfEdgeGroupDataVec[igr].end_ed[0].second);
            halfEdgeGroupDataVec[igr].end_ed[1].second = min(packEdge(ist,i),halfEdgeGroupDataVec[igr].end_ed[1].second);
          }
        }
        assert(igr >= 0);
        // fill half edge group data
        const double dl = DIST(xsp[isp0],xsp[isp1]);
        FOR_J3 halfEdgeGroupDataVec[igr].xc[j] += (xsp[isp0][j]+xsp[isp1][j])*dl;
        halfEdgeGroupDataVec[igr].length += dl;
        halfEdgeGroupDataVec[igr].ned++;
      }
    }
  }
  for (int igr = 0; igr < n_half_edge_groups; ++igr) {
    FOR_I3 halfEdgeGroupDataVec[igr].xc[i] /= (2.0*halfEdgeGroupDataVec[igr].length);
  }

  // second pass for mean distance
  map<uint,int>::iterator it;
  FOR_IST {
    FOR_I3 {
      it = he_to_gr.find(packEdge(ist,i));

      if (it!=he_to_gr.end()) {
        int igr = it->second;
        halfEdgeGroupDataVec[igr].mean_d += DIST(xsp[spost[ist][i]],halfEdgeGroupDataVec[igr].xc);
      }
    }
  }
  for (int igr = 0; igr < n_half_edge_groups; ++igr) {
    halfEdgeGroupDataVec[igr].mean_d /= halfEdgeGroupDataVec[igr].ned;
    // isoperimetric radii quotient: ratio of circle radius with same perimiter to mean_d
    //TODO
    // double
    // halfEdgeGroupDataVec[igr].iso_per_q
  }

  // finished building half-edge groups
  delete[] heg_index;
  delete[] noe_sp;

  return n_half_edge_groups;
}


void SimpleSurface::closeHalfEdgeLoop(vector<NewNode>& xsp_new,vector<NewTri>& spoed_new,const int igroup, const bool ear_clipping,const map<uint,int>& edge_to_group,const vector<HalfEdgeGroupData>& heGroupDataVec) {

  // assumes half-edge groups ahve already been created
  if (edge_to_group.empty()) {
    CWARN("Edges-of-interest have not been grouped; skipping");
    return;
  }

  // get vertex count...
  set<int> dg_v;
  vector<int> sp0_sp1_index;
  bool b_self_arc = false;

  const uint ied0 = heGroupDataVec[igroup].end_ed[0].second;
  const uint ied1 = heGroupDataVec[igroup].end_ed[1].second;

  // bool debug = true;
  // if (debug) {
  //   uint ist0; uchar i0;
  //   unpackEdge(ist0,i0,ied0);
  //   const int isp0 = spost[ist0][i0];
  //   cout << "end_ed[0]: " << (heGroupDataVec[igroup].end_ed[0].first ? "true":"false") << " ist,i: " << ist0 << "," << int(i0) << ", isp: " << isp0 << endl;
  //   uint ist1; uchar i1;
  //   unpackEdge(ist1,i1,ied1);
  //   const int isp1 = spost[ist1][i1];
  //   cout << "end_ed[1]: " << (heGroupDataVec[igroup].end_ed[1].first ? "true":"false") << " ist,i: " << ist1 << "," << int(i1) << ", isp: " << isp1 << endl;
  // }

  if ((ied0 == ied1)&&(heGroupDataVec[igroup].end_ed[0].first)) {
    b_self_arc = true;
  }

  if (heGroupDataVec[igroup].end_ed[0].first == heGroupDataVec[igroup].end_ed[1].first) {
    // if both false, should be a loop whose start/end edges weren't set during grouping b/c it was a loop (like open edges where single chain forms a loop)
    // if both true, start/end edges have been set; should verify that they indicate identical nodes (i.e., subzone boundary loops where redundant edge loops exist)

    // either case, let's ensure the endpoints are the same node...
    uint ist0; uchar i0;
    unpackEdge(ist0,i0,ied0);
    const int isp0 = spost[ist0][i0];
    uint ist1; uchar i1;
    unpackEdge(ist1,i1,ied1);
    const int isp1 = spost[ist1][i1];
    if (isp0 == isp1) b_self_arc = true;
  }

  if (!b_self_arc) {
    CWARN("cannot close edge group " << igroup << " because it is not a closed loop; skipping");
    heGroupDataVec[igroup].dump();
    return;
  }

  vector<uint> halfEdgeVec;
  orderEdgesInGroup(halfEdgeVec,igroup,edge_to_group,heGroupDataVec);

  // have a simple walkable loop describe by halfEdgeVec
  if (halfEdgeVec.size() == 3) {
    closeHalfEdgeLoopWith3Nodes(spoed_new,halfEdgeVec);
  }
  else if (halfEdgeVec.size() == 4) {
    closeHalfEdgeLoopWith4Nodes(spoed_new,halfEdgeVec);
  }
  else if (ear_clipping) {
    // use ear-clipping directly
    closeHalfEdgeLoopEarClipping(spoed_new,halfEdgeVec);
  }
  else {
    // try assuming the mean is visible and, if that fails, use ear-clipping (remove some of the user choices)
    closeHalfEdgeLoopMeanVisible(xsp_new,spoed_new,halfEdgeVec);
  }

}

void SimpleSurface::closeHalfEdgeLoops(vector<int>& group_indices, const bool ear_clipping,const bool open) {
  // closing individual holes
  // assumes half-edge groups have already been created

  if (open) {
    if (n_open_edge_groups == 0) {
      CWARN("open edges have not been grouped; skipping");
      return;
    }
  }
  else {
    if (eoi_to_group.empty()) {
      CWARN("Edges-of-interest have not been grouped; skipping");
      return;
    }
  }

  vector<NewNode> xsp_new;
  vector<NewTri> spoed_new;

  if (open) {
    for (int igr=0, limit=group_indices.size(); igr<limit; igr++) {
      const int igroup = group_indices[igr];
      COUT2("    > closing group: " << igroup);
      closeHalfEdgeLoop(xsp_new,spoed_new,igroup,ear_clipping,oe_to_group,openEdgeGroupDataVec);
    }
  }
  else {
    for (int igr=0, limit=group_indices.size(); igr<limit; igr++) {
      const int igroup = group_indices[igr];
      COUT2("    > closing group: " << igroup);
      closeHalfEdgeLoop(xsp_new,spoed_new,igroup,ear_clipping,eoi_to_group,eoiGroupDataVec);
    }
  }

  // update data
  const int ss_nst0 = nst;
  const int ss_nsp0 = nsp;

  // update nsp dependent data
  if (xsp_new.size() > 0) {

    // update xsp
    nsp += xsp_new.size();
    growNspData(nsp,ss_nsp0);
    for (int isp = ss_nsp0; isp < nsp; ++isp) {
      FOR_I3 xsp[isp][i] = xsp_new[isp-ss_nsp0].xsp[i];
    }
    xsp_new.clear();
  }

  // update nst dependent data
  if (spoed_new.size() > 0) {
    nst += spoed_new.size();

    growNstData(nst,ss_nst0);
    for (int ist = ss_nst0; ist < nst; ++ist) {
      FOR_I3 spost[ist][i] = spoed_new[ist-ss_nst0].spost[i];
      znost[ist] = spoed_new[ist-ss_nst0].znost;
    }
    spoed_new.clear();

    // update zoneVec
    stringstream oss; oss << n_closed_edge_groups++;
    string zone_name = "Closed_Loops_" + oss.str();
    addNewZone(zone_name); // also resizes szozn_i

    // update szost
    szost.resize(nst);

    for (int ist = ss_nst0; ist < nst; ++ist) szost[ist] = spoed_new[ist-ss_nst0].szost; // everything thrown into one sz
    assert(szozn_i.getLength() == int(zoneVec.size())+1);
    szozn_i[zoneVec.size()] += 1; // add sub-zone in new zone
    nsz++;

    // teost and (implicitly) open edge groups need to be updated
    clearTeost();

    if (open) clearOpenEdgeGroups();
    clearDynamicEdgeGroups();  // could be open, so just clear

    // addition of new zones requires rebuilding zone/subzone info
    clearSubzoneData();
    clearZoneData();
  }
}

void SimpleSurface::closeHalfEdgeDonutLoops(vector<pair<int,int> >& donutLoopsVec,const bool open) {
  if (donutLoopsVec.empty()) return;

  // Actually close the pairs and return the newTris formed by the closure
  vector<NewTri> newTris;
  for (vector<pair<int,int> >::iterator it=donutLoopsVec.begin(); it!=donutLoopsVec.end(); ++it) {
    closeHalfEdgeDonutLoop(newTris,it->first,it->second,zoneVec.size(),nsz,open);
  }

  // Add teh new tries to the zones of surface tris and subzones of surface tris structure
  if (!newTris.empty()) {
    // set zone and znost to a new zone
    const int nzn = zoneVec.size();
    for (vector<NewTri>::iterator it=newTris.begin(); it!=newTris.end(); ++it) {
      it->znost = nzn;
      it->szost = nsz;
    }
  }

  // update data
  const int ss_nst0 = nst;
  nst += int(newTris.size());

  // update spost
  assert(ss_nst0 > 0);
  assert(znost != NULL);
  growNstData(nst,ss_nst0);

  for (int ist = ss_nst0; ist < nst; ++ist) {
    FOR_I3 spost[ist][i] = newTris[ist-ss_nst0].spost[i];
    znost[ist] = newTris[ist-ss_nst0].znost;
  }

  stringstream oss; oss << n_closed_edge_groups++;
  string zone_name = "Closed_Loops_" + oss.str();
  addNewZone(zone_name); // also resizes szozn_i

  // update szost
  szost.resize(nst);
  for (int ist = ss_nst0; ist < nst; ++ist) szost[ist] = newTris[ist-ss_nst0].szost; // everything thrown into one sz
  assert(szozn_i.getLength() == int(zoneVec.size())+1);
  szozn_i[zoneVec.size()] += 1; // add sub-zone from new zone
  nsz++;

  // teost and (implicitly) edge groups need to be updated
  clearTeost();
  if (open) clearOpenEdgeGroups();
  clearDynamicEdgeGroups();  // could be open, so just clear

  // addition of new zones requires rebuilding zone/subzone info
  clearSubzoneData();
  clearZoneData();
}

void SimpleSurface::closeHalfEdgeDonutLoop(vector<NewTri>& newTris,const int group0, const int group1, const int izone, const int isz,const bool open) {

  ensureTeost();

  //check validity of groups passed in
  int n_groups;
  if (open) n_groups = n_open_edge_groups;
  else n_groups = eoiGroupDataVec.size();

  if (group0<0 || group0>n_groups) {
    CWARN("group index " << group0 << " is out-of-bounds; skipping donut hole closing");
    return;
  }
  if (group1<0 || group1>n_groups) {
    CWARN("group index " << group1 << " is out-of-bounds; skipping donut hole closing");
    return;
  }

  // Build the half edge vector associated with each open edge to be closed
  vector<uint> halfEdgeVec1;
  vector<uint> halfEdgeVec0;

  if (open) {
  // assumes open-edge groups have already been created
    if (oe_to_group.empty()) {
      CWARN("open edges have not been grouped; skipping");
      return;
    }
    orderEdgesInGroup(halfEdgeVec0,group0,oe_to_group,openEdgeGroupDataVec);
    orderEdgesInGroup(halfEdgeVec1,group1,oe_to_group,openEdgeGroupDataVec);
  }
  else {
    // assumes half-edge groups have already been created
    if (eoi_to_group.empty()) {
      CWARN("Edges-of-interest have not been grouped; skipping");
      return;
    }
    orderEdgesInGroup(halfEdgeVec0,group0,eoi_to_group,eoiGroupDataVec);
    orderEdgesInGroup(halfEdgeVec1,group1,eoi_to_group,eoiGroupDataVec);
  }

  //Get the node indices on each halfEdgeVec
  vector<int> isp0Vec(halfEdgeVec0.size());
  for (vector<uint>::iterator it=halfEdgeVec0.begin(); it!=halfEdgeVec0.end(); ++it) {
    uint ist; uchar i;
    unpackEdge(ist,i,*it);
    const int isp = spost[ist][i];
    isp0Vec[it-halfEdgeVec0.begin()] = isp;
  }
  cout << " > nodes on loop0: " << isp0Vec.size() << endl;

  vector<int> isp1Vec(halfEdgeVec1.size());
  for (vector<uint>::iterator it=halfEdgeVec1.begin(); it!=halfEdgeVec1.end(); ++it) {
    uint ist; uchar i;
    unpackEdge(ist,i,*it);
    const int isp = spost[ist][i];
    isp1Vec[it-halfEdgeVec1.begin()] = isp;
  }
  cout << " > nodes on loop1: " << isp1Vec.size() << endl;

  // Close the damn thing
  facetGap(newTris,isp0Vec,isp1Vec,true);

  // Associate the new tris with an additional zone and subzone (izone and isz are equal to the N_zone and N_sz)
  for (vector<NewTri>::iterator it=newTris.begin(); it!=newTris.end(); ++it) {
    it->znost = izone;
    it->szost = isz;
  }
}//closeHalfEdgeDonutLoop()

void writeXpToTecplotFile(const string& filename,const vector<NewNode>& nodes) {
  FILE * fp = fopen(filename.c_str(),"w");
  // if ((data==NULL) && (flag==NULL)) {
    fprintf(fp,"VARIABLES=\"X\" \"Y\" \"Z\"\n");
    fprintf(fp,"ZONE T=\"nodes\"\n");
    for (vector<NewNode>::const_iterator it=nodes.begin(); it!=nodes.end(); ++it) {
      fprintf(fp,"%18.15e %18.15e %18.15e\n",it->xsp[0],it->xsp[1],it->xsp[2]);
    }
  // }
  // else if (data==NULL) {
  //   fprintf(fp,"VARIABLES=\"X\" \"Y\" \"Z\" \"flag\"\n");
  //   for (int ip = 0; ip < np; ++ip) {
  //     fprintf(fp,"%18.15e %18.15e %18.15e %d\n",xp[ip][0],xp[ip][1],xp[ip][2],flag[ip]);
  //   }
  // }
  // else {
  //   assert(data!=NULL && flag !=NULL);
  //   fprintf(fp,"VARIABLES=\"X\" \"Y\" \"Z\" \"flag\" \"data\"\n");
  //   for (int ip = 0; ip < np; ++ip) {
  //     fprintf(fp,"%18.15e %18.15e %18.15e %d %18.15e\n",xp[ip][0],xp[ip][1],xp[ip][2],flag[ip],data[ip]);
  //   }
  // }
  fclose(fp);
}

void writeTrisToTecplotFile(const string& filename,const vector<NewTri>& tris,const vector<NewNode>& nodes) {
  FILE * fp = fopen(filename.c_str(),"w");
  // if ((data==NULL) && (flag==NULL)) {
    fprintf(fp,"VARIABLES=\"X\" \"Y\" \"Z\"\n");
    fprintf(fp,"ZONE T=\"tris\"\n");
    const int nst_new = tris.size();
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nst_new*3,nst_new);
    for (vector<NewTri>::const_iterator it=tris.begin(); it!=tris.end(); ++it) {
      FOR_I3 fprintf(fp,"%18.15e %18.15e %18.15e\n",nodes[it->spost[i]].xsp[0],nodes[it->spost[i]].xsp[1],nodes[it->spost[i]].xsp[2]);
    }
    for (vector<NewTri>::const_iterator it=tris.begin(); it!=tris.end(); ++it) {
      FOR_I3 fprintf(fp,"%d %d %d\n",it->spost[0],it->spost[1],it->spost[2]);
    }

  fclose(fp);
}

void writeTrisToTecplotFile(const string& filename,const vector<NewTri>& tris,const double (*xsp)[3]) {
  FILE * fp = fopen(filename.c_str(),"w");
  fprintf(fp,"VARIABLES=\"X\" \"Y\" \"Z\"\n");
  fprintf(fp,"ZONE T=\"tris\"\n");
  const int nst_new = tris.size();
  fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nst_new*3,nst_new);
  for (vector<NewTri>::const_iterator it=tris.begin(); it!=tris.end(); ++it) {
    FOR_I3 fprintf(fp,"%18.15e %18.15e %18.15e\n",xsp[it->spost[i]][0],xsp[it->spost[i]][1],xsp[it->spost[i]][2]);
  }
  int count = 1;  // tecplot is 1-indexed
  for (vector<NewTri>::const_iterator it=tris.begin(); it!=tris.end(); ++it) {
    fprintf(fp,"%d %d %d\n",count,count+1,count+2);
    count += 3;
  }

  fclose(fp);
}

void SimpleSurface::getNodesFromHalfEdges(vector<int>& nodes,const vector<uint>& edges,const bool loop) const {
  uint ist; uchar i;
  for (vector<uint>::const_iterator it=edges.begin(); it!=edges.end(); ++it) {
    unpackEdge(ist,i,*it);
    nodes.push_back(spost[ist][i]);
  }
  if (!loop) {
    // need to process end node if not a loop
    unpackEdge(ist,i,edges.back());
    nodes.push_back(spost[ist][i]);
  }
}

double SimpleSurface::getNodesAndLengthFromHalfEdges(vector<pair<int,double> >& nodes,const vector<uint>& edges,const bool loop) const {
  double length = 0.0;
  uint ist; uchar i;
  for (vector<uint>::const_iterator it=edges.begin(); it!=edges.end(); ++it) {
    unpackEdge(ist,i,*it);
    nodes.push_back(pair<int,double> (spost[ist][i],length));
    length += DIST(xsp[spost[ist][i]],xsp[spost[ist][(i+1)%3]]);
  }

  if (!loop) {
    unpackEdge(ist,i,edges.back());
    nodes.push_back(pair<int,double> (spost[ist][(i+1)%3],length));  // need to process end node if not a loop
  }
  return length;  // total length if a loop
}

void SimpleSurface::facetGap(vector<NewTri>& newTris,const vector<int>& isp0Vec,const vector<int>& isp1Vec,const bool loop) const {
  // expected new tri count for loops
  int n_expected = isp0Vec.size() + isp1Vec.size();
  if (!loop) n_expected -= 2;

  int (*tmp_spost)[3] = new int[n_expected][3];
  const int n_new = GeomUtils::facetGap(tmp_spost,isp0Vec,isp1Vec,xsp,loop);
  assert(n_new == n_expected);

  // add to new tris vec
  for (int ist=0; ist<n_new; ++ist) {
    newTris.push_back(NewTri(tmp_spost[ist],0,0));
  }

  DELETE(tmp_spost);
}

void SimpleSurface::facetGap(vector<NewTri>& newTris,const vector<int>& isp0Vec,const vector<int>& isp1Vec,const double (*my_xsp)[3],const bool loop) const {
  // expected new tri count for loops
  int n_expected = isp0Vec.size() + isp1Vec.size();
  if (!loop) n_expected -= 2;

  int (*tmp_spost)[3] = new int[n_expected][3];
  const int n_new = GeomUtils::facetGap(tmp_spost,isp0Vec,isp1Vec,my_xsp,loop);
  assert(n_new == n_expected);

  // add to new tris vec
  for (int ist=0; ist<n_new; ++ist) {
    newTris.push_back(NewTri(tmp_spost[ist],0,0));
  }

  DELETE(tmp_spost);
}

void SimpleSurface::facetGapFromHalfEdges(vector<NewTri>& newTris,const vector<uint>& ise0Vec,const vector<uint>& ise1Vec,const bool loop) const {

  // convert half-edges to nodes
  vector<int> isp0Vec;
  getNodesFromHalfEdges(isp0Vec,ise0Vec,loop);
  vector<int> isp1Vec;
  getNodesFromHalfEdges(isp1Vec,ise1Vec,loop);

  facetGap(newTris,isp0Vec,isp1Vec,loop);
}

bool sortBySecondFirst(const pair<int,double> &a,const pair<int,double> &b) {
    if (a.second == b.second) {
      return (a.first < b.first);  // secondary sort if length is the same
    }
    else {
      return (a.second < b.second);
    }
}

void SimpleSurface::reDiscretizeEdges(const vector<int>& edge_indices,const double delta,const int type,const bool b_keep_edge_nodes) {
  COUT2("SimpleSurface::reDiscretizeEdges()");

  // need dynamic edge groups
  if (eoi_type == NONE || eoi_type == MULTI) {
    CWARN("cannot discretize dynamic edges when type is NONE or MULTI; skipping");
    return;
  }
  ensureDynamicEdgeGroups();

  // writeTecplot("pre_discretize.dat");

  COUT2(" > user specified delta: " << delta);
  COUT2(" > type: " << type);
  bool no_adj_feature = false;
  if (type == 1) no_adj_feature = true;

  // assume for now that user is only passing dynamic edges (can use to process open edges if needed)
  // assuming selected number of edges won't be too large, we can build unique edge map based on
  // (isp0,isp1). Order is min,max isp value, where a -index indicates order in this chain direction is reversed
  // we should be able to discretize the chains, keeping track of new nodes and old, and re-faceting adjacent tris
  // as necessary

  vector<NewNode> fixedNodes;  // new fixed nodes

  vector<pair<int,int> > featureEdgeVec;
  vector<pair<int,int> > adj_tris;
  vector<pair<int,int> > isp_to_nearest_fixed;
  identifyFeatureXp(fixedNodes,featureEdgeVec,adj_tris,isp_to_nearest_fixed,edge_indices,delta,no_adj_feature,b_keep_edge_nodes);

  multimap <pair<int,int>,int> edge_to_fixedNodes;
  map<int,pair<int,bool> > isp_to_fixedNodes;  // pair<int,bool> indicates fixed new node index, if protected node associated
  // create multimap from (isp_min,isp_max) to protected nodes
  // create map from isp to protected nodes
  for (vector<NewNode>::const_iterator it=fixedNodes.begin(); it!=fixedNodes.end(); ++it) {
    if (it->isp < 0) {
      // registered to an old edge
      const int ied = -it->isp-1;
      const int ino = it-fixedNodes.begin();
      const int isp0 = featureEdgeVec[ied].first;
      const int isp1 = featureEdgeVec[ied].second;
      edge_to_fixedNodes.insert(pair<pair<int,int>,int> (pair<int,int>(isp0,isp1),ino));
    }
    else {
      assert(it->isp >=0 && it->isp < nsp);
      const int ino = it-fixedNodes.begin();
      isp_to_fixedNodes.insert(pair<int,pair<int,bool> >(it->isp,pair<int,bool> (ino,true)));
    }
  }
  for (vector<pair<int,int> >::iterator it=isp_to_nearest_fixed.begin(); it!=isp_to_nearest_fixed.end(); ++it) {
    isp_to_fixedNodes.insert(pair<int,pair<int,bool> >(it->first,pair<int,bool> (it->second,false)));
  }

  // containers to hold all new surface information
  vector<NewTri> new_adj_tris;

  // process adjacent tris and add replacement tris (some new nodes potentially as well)
  zipAdjacentTris(new_adj_tris,adj_tris,edge_to_fixedNodes,isp_to_fixedNodes,fixedNodes);
  const int np_fixed = fixedNodes.size();

  // positive values will be deleted
  st_flag.resize(nst);
  st_flag.setAll(-1);
  for (vector<pair<int,int> >::iterator it=adj_tris.begin(); it!=adj_tris.end(); ++it) {
    st_flag[it->first] = 1;
  }

  // cull old surface and replace with new
  vector<NewTri> new_tris;  // empty placeholder
  updateSurfaceWithRetessellation(fixedNodes,np_fixed,new_adj_tris,new_tris);
}
