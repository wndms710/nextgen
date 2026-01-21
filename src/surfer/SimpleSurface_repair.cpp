#include "SimpleSurface.hpp"
#include "Adt.hpp"
#include "WebUI.hpp"

void SimpleSurface::deleteTrisWithCollocatedNodes() {
  COUT1("SimpleSurface::deleteTrisWithCollocatedNodes()");
  int nst_new = 0;
  FOR_IST {
    if ( (spost[ist][0] != spost[ist][1]) && (spost[ist][1] != spost[ist][2]) && (spost[ist][2] != spost[ist][0]) ) {
      const int ist_new = nst_new++;
      if (ist_new != ist) {
        assert(ist_new < ist);
        FOR_I3 spost[ist_new][i] = spost[ist][i];
        znost[ist_new] = znost[ist];
        szost[ist_new] = szost[ist];
      }
    }
  }
  cout << "nst: " << nst << " nst_new: " << nst_new << endl;
  if (nst_new != nst) {
    clearTeost();
    clearNonManifoldData();
  }
  nst = nst_new;
}

void SimpleSurface::deleteTrisWithIdenticalNodes() {
  cout << "SimpleSurface::deleteTrisWithIdenticalNodes()" << endl;
  // we may not be able to "teost" across the nodes, so build stosp_i/v...
  ensureTeost();
  st_flag.resize(nst);
  st_flag.setAll(0);
  sp_flag.resize(nst);
  sp_flag.setAll(-1);
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 sp_flag[spost[ist][i]] = ist;
    FOR_I3 {
      int ime,orient_ime;
      if (isEdgeMulti(ime,orient_ime,ist,i)) {
        vector<pair<int,int> > ist_nbrs;
        teost.getMultiNbrs(ist_nbrs,ist,i);
        for (vector<pair<int,int> >::const_iterator it=ist_nbrs.begin(); it!=ist_nbrs.end(); ++it) {
          const int ist_nbr = it->first;
          if (ist_nbr != ist) {
            if ((sp_flag[spost[ist_nbr][0]] == ist)&&
                (sp_flag[spost[ist_nbr][1]] == ist)&&
                (sp_flag[spost[ist_nbr][2]] == ist)) {
              // this is a match!
              if (st_flag[ist] == 0) cout << " > ist " << ist << " and " << ist_nbr << " have identical nodes" << endl;
              st_flag[ist] = 1;
              st_flag[ist_nbr] = 1;
            }
          }
        }
      }
    }
  }
  deleteFlaggedTris();
}

void SimpleSurface::mergeCollocatedNodes() {

  COUT1("SimpleSurface::mergeCollocatedNodes()");

  // surface may change slightly, recompute geoms
  // ME: is this right?...
  b_centroid = false;
  b_rmax = false;
  b_boundingBox = false;
  clearTeost();
  clearNonManifoldData();
  b_update_hidden_subzones = true;  // will update hidden sz in the prune method

  // use position to sort and group collocated nodes
  vector<pair<DoubleVertex,int> > xspIspVec(nsp);
  FOR_ISP xspIspVec[isp] = pair<DoubleVertex,int>(DoubleVertex(xsp[isp]),isp);
  sort(xspIspVec.begin(),xspIspVec.end());

  // make sp_flag point to unique isp
  sp_flag.setLength(nsp);
  FOR_ISP sp_flag[isp] = isp;
  for (int isp = 1; isp < nsp; ++isp) {
    if (xspIspVec[isp].first == xspIspVec[isp-1].first) {
      const int min_isp = min(xspIspVec[isp].second,xspIspVec[isp-1].second);
      sp_flag[xspIspVec[isp].second] = sp_flag[min_isp];
      sp_flag[xspIspVec[isp-1].second] = sp_flag[min_isp];
    }
  }

  // now compress sp_flag & xsp
  int nsp_new = 0;
  FOR_ISP {
    if (sp_flag[isp] == isp) {
      if (nsp_new != isp) {
        FOR_I3 xsp[nsp_new][i] = xsp[isp][i];
      }
      sp_flag[isp] = nsp_new;
      ++nsp_new;
    }
    else {
      sp_flag[isp] = sp_flag[sp_flag[isp]];
    }
    assert(sp_flag[isp] < nsp_new);
  }

  // update spost
  FOR_IST {
    FOR_I3 {
      spost[ist][i] = sp_flag[spost[ist][i]];
    }
  }
  cout << " > nsp: " << nsp << " nsp_new: " << nsp_new << endl;
  xspIspVec.clear();
  if (nsp_new != nsp) clearTeost();
  nsp = nsp_new;
}

void SimpleSurface::mergeCollocatedNodes(const double eps) {

  COUT1("SimpleSurface::mergeCollocatedNodes: eps: " << eps);

  if (eps <= 0.0) {
    mergeCollocatedNodes();
    return;
  }

  // surface may change slightly, recompute geoms
  // ME: is this right?...
  b_centroid = false;
  b_rmax = false;
  b_boundingBox = false;
  clearTeost();
  clearNonManifoldData();
  b_update_hidden_subzones = true;  // will update hidden sz in the prune method

  // use the bbox adt to store the zero-sized bbox for each point...
  Adt<double> xsp_adt(nsp,xsp,xsp);

  // use sp_flag to store the index, starting with the current isp...
  sp_flag.resize(nsp);
  for (int isp = 0; isp < nsp; ++isp)
    sp_flag[isp] = isp;

  vector<int> intVec;
  for (int isp = 0; isp < nsp; ++isp) {
    assert(intVec.empty());
    xsp_adt.buildListForSphere(intVec,xsp[isp],eps);
    // for (int ii = 0; ii < intVec.size(); ++ii) {
    for (vector<int>::const_iterator it=intVec.begin(); it!=intVec.end(); ++it) {
      const int isp_nbr = *it;
      if (isp != isp_nbr) {
        // these guys need to be connected...
        int isp_ = sp_flag[isp];
        while (isp_ != sp_flag[isp_])
          isp_ = sp_flag[isp_];
        int isp_nbr_ = sp_flag[isp_nbr];
        while (isp_nbr_ != sp_flag[isp_nbr_])
          isp_nbr_ = sp_flag[isp_nbr_];
        // set both to the minimum...
        sp_flag[isp_] = sp_flag[isp_nbr_] = min(isp_,isp_nbr_);
      }
    }
    intVec.clear();
  }

  // now loop through the nodes and set nodes that are going to stay using -1 indexing, and nodes that are going to be deleted
  // to their correct -1 index...

  double d2_max = 0.0;
  int nsp_new = 0;
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] == isp) {
      // this is a keeper!...
      // for now just take the coords of the minimum index node. We could average
      // also, and then normalize at the end?...
      FOR_I3 xsp[nsp_new][i] = xsp[isp][i];
      sp_flag[isp] = -nsp_new-1; // -1,-2,-3...
      ++nsp_new;
    }
    else {
      // traverse down to the bottom of this condesed set of nodes...
      int isp_ = sp_flag[isp];
      while (sp_flag[isp_] >= 0)
        isp_ = sp_flag[isp_];
      sp_flag[isp] = sp_flag[isp_];
      // we know node isp is being replace by the node at new location -sp_flag[isp_]-1:
      const int isp_new = -sp_flag[isp_]-1;
      d2_max = max(d2_max,DIST2(xsp[isp_new],xsp[isp]));
    }
  }

  cout << " > nsp: " << nsp << " nsp_new: " << nsp_new << ", max dist: " << sqrt(d2_max) << endl;
  nsp = nsp_new;

  // now sp_flag contains the new node numbers, but in -1 indexing...
  // use this as an opportunity to remove collapsed tris as well...
  int nst_new = 0;
  for (int ist = 0; ist < nst; ++ist) {
    const int isp0 = -sp_flag[spost[ist][0]]-1;
    const int isp1 = -sp_flag[spost[ist][1]]-1;
    const int isp2 = -sp_flag[spost[ist][2]]-1;
    if ((isp0 != isp1)&&(isp0 != isp2)&&(isp1 != isp2)) {
      spost[nst_new][0] = isp0;
      spost[nst_new][1] = isp1;
      spost[nst_new][2] = isp2;
      znost[nst_new] = znost[ist];
      szost[nst_new] = szost[ist];
      ++nst_new;
    }
  }

  if (nst_new != nst) {
    cout << " > " << nst-nst_new << " tris with collocated nodes were removed." << endl;
    nst = nst_new;
  }

}

string SimpleSurface::replaceSpacesWithUnderscores(const string& newname) {
  string new_name = newname;
  // automatically replace spaces with underscores
  const string space = " ";
  std::size_t found = new_name.find(space);
  while (found != string::npos) {
    new_name.replace(found,space.length(),"_");
    found = new_name.find(space);
  }
  if (new_name != newname) {
    WUI(WARN,"Zone names cannot contain spaces; replaced spaces with _");
  }
  return new_name;

}


// decided to build one loop at a time due to isp's begin part of multiple loops
// (i.e. can't use sp_flag to walk everyone at once). Maybe we can do better
void SimpleSurface::buildHalfEdgeVec(vector<uint>& halfEdgeVec,const map<const int,int> groupToLoop,const int iloop) {

  halfEdgeVec.clear();

  // need open edge groups (and teost)
  if (eoi_type == NONE) {
    CWARN("cannot build half-edge structures when dynamic edge type is NONE; skipping");
    return;
  }
  ensureDynamicEdgeGroups();

  // this counts edges and sets sp_flag for walking
  sp_flag.setLength(nsp);
  sp_flag.setAll(-1);
  int ned = 0;
  uint ied = 0;
  map<const int,int>::const_iterator cit;
  FOR_IST {
    FOR_I3 {
      int igr;
      if (getHalfEdgeGroup(igr,ist,i)) {
        // put our edge index in one of our nodes...
        if (sp_flag[spost[ist][i]] == -1) {
          cit = groupToLoop.find(igr);
          // is this group in a loop?
          if (cit != groupToLoop.end()) {
            // is this group in our loop?
            if (cit->second == iloop) {
              ied = packEdge(ist,i);
              sp_flag[spost[ist][i]] = int(ied);
              ned++;
            }
          }
        }
        else {
          assert(0); // unwalkable loop. should not get here
        }
      }
    }
  }
  assert(ied != 0);

  // walk the loop
  halfEdgeVec.resize(ned);
  uint ist; uchar i;
  unpackEdge(ist,i,ied);
  ned = 0;
  while (sp_flag[spost[ist][i]] != -1) {
    sp_flag[spost[ist][i]] = -1; // set to -1 indicating its been visited
    halfEdgeVec[ned++] = ied; // put in current open edge
    ied = uint(sp_flag[spost[ist][(i+1)%3]]); // go to next adjacent edge
    unpackEdge(ist,i,ied);
  }
  //cout << ned << " " << halfEdgeVec.size() << endl;
  assert(ned == int(halfEdgeVec.size())); // walked
}

void SimpleSurface::closeHalfEdgeLoopMeanVisible(vector<NewNode>& xsp_new,vector<NewTri>& spoed_new,vector<uint>& halfEdgeVec) {

  // populate spoed by walking loop
  int ned_loop = halfEdgeVec.size();
  double loop_centroid[3] = {0.0,0.0,0.0};
  double circ = 0.0;
  for (vector<uint>::const_iterator it=halfEdgeVec.begin(); it!=halfEdgeVec.end(); ++it) {
    uint ist; uchar i;
    unpackEdge(ist,i,*it);
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    assert(isp0 != isp1);
    const double dl = DIST(xsp[isp0],xsp[isp1]);
    FOR_I3 loop_centroid[i] += (xsp[isp0][i]+xsp[isp1][i])*dl;
    circ += dl;
  }
  FOR_I3 loop_centroid[i] /= (2.0*circ);

  // the mean is visible if all tris created from centroid and open edges are not flipped wrt to each other
  bool mean_visible = true;
  {
    uint ist; uchar i;
    unpackEdge(ist,i,halfEdgeVec[0]);
    const double n0[3] = TRI_NORMAL_2(xsp[spost[ist][(i+1)%3]],xsp[spost[ist][i]],loop_centroid);
    for (vector<uint>::const_iterator it=halfEdgeVec.begin()+1; it!=halfEdgeVec.end(); ++it) {
      unpackEdge(ist,i,*it);
      const double n1[3] = TRI_NORMAL_2(xsp[spost[ist][(i+1)%3]],xsp[spost[ist][i]],loop_centroid);
      if (DOT_PRODUCT(n0,n1) < 0.0) {
        mean_visible = false;
        break;
      }
    }
  }

  if (!mean_visible) {
    cout << "The loop mean is not visible to all the vertices. Attempting to use ear-clipping algorithm..." << endl;
    closeHalfEdgeLoopEarClipping(spoed_new,halfEdgeVec);
    return;
  }

  // ========================================
  // perform triangulation...
  // ========================================

  // populate xsp with mean
  const int ss_nst0 = nst;
  const int xsp_new_size0 = xsp_new.size();
  xsp_new.push_back(NewNode(loop_centroid));
  cout << "Added vertex: " << xsp_new[xsp_new_size0].xsp[0] << " " << xsp_new[xsp_new_size0].xsp[1] << " " << xsp_new[xsp_new_size0].xsp[2] << endl;

  // populate spoed by walking loop
  const int spoed_new_size0 = spoed_new.size();
  spoed_new.resize(spoed_new_size0 + ned_loop);
  int ist_nbr = ss_nst0;
  for (vector<uint>::const_iterator it=halfEdgeVec.begin(); it!=halfEdgeVec.end(); ++it) {
    uint ist; uchar i;
    unpackEdge(ist,i,*it);
    spoed_new[(ist_nbr-ss_nst0)+spoed_new_size0].spost[0] = spost[ist][(i+1)%3];
    spoed_new[(ist_nbr-ss_nst0)+spoed_new_size0].spost[1] = spost[ist][i];
    spoed_new[(ist_nbr-ss_nst0)+spoed_new_size0].spost[2] = xsp_new_size0+nsp;
    spoed_new[(ist_nbr-ss_nst0)+spoed_new_size0].znost = zoneVec.size();
    spoed_new[(ist_nbr-ss_nst0)+spoed_new_size0].szost = nsz;
    ist_nbr++;
  }

  cout << "Finished closing loop" << endl;

}

void SimpleSurface::closeHalfEdgeLoopWith3Nodes(vector<NewTri>& spoed_new,vector<uint>& halfEdgeVec) {

  // Close by adding a new tri...

  const int spoed_new_size0 = spoed_new.size();
  spoed_new.resize(spoed_new_size0 + 1);

  FOR_J3 {
    uint ist; uchar i;
    unpackEdge(ist,i,halfEdgeVec[j]);
    spoed_new[spoed_new_size0].spost[2-j] = spost[ist][(i+1)%3];
  }
  spoed_new[spoed_new_size0].znost = zoneVec.size();
  spoed_new[spoed_new_size0].szost = nsz;

}

void SimpleSurface::closeHalfEdgeLoopWith4Nodes(vector<NewTri>& spoed_new,vector<uint>& halfEdgeVec) {

  // Close by adding a single edge b/w two new tris. Edge is chosen to minimize crease
  // angles between adjacent tris....

  const int spoed_new_size0 = spoed_new.size();
  spoed_new.resize(spoed_new_size0 + 2);

  double cos_sum[2] = {0.0,0.0};
  FOR_K2 {
    FOR_J2 {
      uint ied0 = halfEdgeVec[2*j+k];
      uint ist0; uchar i0;
      unpackEdge(ist0,i0,ied0);
      const double n0[3] = TRI_NORMAL_2(xsp[spost[ist0][i0]],xsp[spost[ist0][(i0+1)%3]],xsp[spost[ist0][(i0+2)%3]]);
      const double n0_mag = MAG(n0);
      if (n0_mag > 0.0) {
        uint ist; uchar i;
        unpackEdge(ist,i,halfEdgeVec[(2*j+k+1)%4]);
        const double n[3] = TRI_NORMAL_2(xsp[spost[ist][i]],xsp[spost[ist][(i+1)%3]],xsp[spost[ist][(i+2)%3]]);
        const double n_mag = MAG(n);
        if (n_mag > 0.0)
          cos_sum[k] += DOT_PRODUCT(n0,n)/(n0_mag*n_mag);
      }
    }
  }

  int edge_start;
  if (cos_sum[0] >= cos_sum[1])
    edge_start = 3; // form tris from 3,0 and 1,2
  else
    edge_start = 0; // form tris from 0,1 and 2,3

  FOR_J2 {
    uint ist0; uchar i0;
    unpackEdge(ist0,i0,halfEdgeVec[(edge_start+2*j)%4]);
    spoed_new[spoed_new_size0+j].spost[1] = spost[ist0][i0];
    spoed_new[spoed_new_size0+j].spost[0] = spost[ist0][(i0+1)%3];
    uint ist; uchar i;
    unpackEdge(ist,i,halfEdgeVec[(edge_start+2*j+1)%4]);
    assert(spoed_new[spoed_new_size0+j].spost[0] == spost[ist][i]);
    spoed_new[spoed_new_size0+j].spost[2] = spost[ist][(i+1)%3];
  }

  FOR_I2 {
    spoed_new[spoed_new_size0+i].znost = zoneVec.size();
    spoed_new[spoed_new_size0+i].szost = nsz;
  }

}

void SimpleSurface::closeHalfEdgeLoopEarClipping(vector<NewTri>& spoed_new,vector<uint>& halfEdgeVec) {

  // extract vertices bounding the open loop...
  int nv = halfEdgeVec.size();
  double (*xloop)[3] = new double[nv][3];
  int *isp_loop = new int[nv];
  double loop_centroid[3] = {0.0,0.0,0.0};
  double circ = 0;
  for (vector<uint>::const_iterator it=halfEdgeVec.begin(); it!=halfEdgeVec.end(); ++it) {
    const int ii = it-halfEdgeVec.begin();
    uint ist; uchar i;
    unpackEdge(ist,i,*it);
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    assert(isp0 != isp1);
    const double dl = DIST(xsp[isp0],xsp[isp1]);
    FOR_I3 loop_centroid[i] += (xsp[isp0][i]+xsp[isp1][i])*dl;
    FOR_J3 xloop[ii][j] = xsp[isp0][j];
    isp_loop[ii] = isp0;
    circ += dl;
  }
  FOR_I3 loop_centroid[i] /= (2.0*circ);
  // cout << "number of vertices: " << nv << " centroid: " << COUT_VEC(loop_centroid) << endl;

  // shift loop points by mean
  //for (int iv = 0; iv < nv; ++iv) FOR_I3 xloop[iv][i] -= loop_centroid[i];
  //FOR_I3 loop_centroid[i] = 0.0;

  // scale each direction to the unit cube
  //for (int iv = 0; iv < nv; ++iv) FOR_I3 xloop[iv][i] /= max(1E-4,bb_dx[i]);

  double area_vec_2[3] = {0.0,0.0,0.0};
  for (int iv = 0; iv < nv; ++iv) {
    int iv1 = iv + 1;
    if (iv1 == nv) iv1 = 0;
    double p0[3] = DIFF(xloop[iv1],loop_centroid);
    double p1[3] = DIFF(xloop[iv],loop_centroid);
    area_vec_2[0] += CROSS_PRODUCT_0(p0,p1);
    area_vec_2[1] += CROSS_PRODUCT_1(p0,p1);
    area_vec_2[2] += CROSS_PRODUCT_2(p0,p1);
  }
  // double mag_area_vec_2 = MAG(area_vec_2);
  double loop_normal[3];
  FOR_I3 loop_normal[i] = area_vec_2[i];
  // cout << "loop normal vector: " << COUT_VEC(loop_normal) << endl;
  // double area = DOT_PRODUCT(area_vec_2,area_vec_2)/mag_area_vec_2/2.0;
  // cout << "area of loop: " << area << endl;
  double loop_tangent0[3] = DIFF(xloop[0],loop_centroid);
  // double mag_tan0 = MAG(loop_tangent0);
  // cout << "loop tangent vector: " << COUT_VEC(loop_tangent0) << endl;
  double loop_tangent1[3] = CROSS_PRODUCT(loop_normal,loop_tangent0);
  // double mag_tan1 = MAG(loop_tangent1);
  // cout << "loop tangent vector: " << COUT_VEC(loop_tangent1) << endl;

  // project onto tangent
  // x = t_0*e_0 + t_1*e_1 + s*n
  // x = xloop[3]
  // e_0 = loop_tangent0[3]
  // e_1 = loop_tangent1[3]
  // n = loop_normal[3]
  double *t0 = new double[nv];
  double *t1 = new double[nv];
  double *s = new double[nv];
  double max_dist = -1.0;
  double mag_area_vec_2_sq = DOT_PRODUCT(area_vec_2,area_vec_2);
  double mag_tan0_sq = DOT_PRODUCT(loop_tangent0,loop_tangent0);
  double mag_tan1_sq = DOT_PRODUCT(loop_tangent1,loop_tangent1);
  for (int iv = 0; iv < nv; ++iv) {
    t0[iv] = DOT_PRODUCT(xloop[iv],loop_tangent0)/mag_tan0_sq;
    t1[iv] = DOT_PRODUCT(xloop[iv],loop_tangent1)/mag_tan1_sq;
    s[iv] = DOT_PRODUCT(xloop[iv],loop_normal)/mag_area_vec_2_sq;
    double xloop_check[3];
    FOR_I3 xloop_check[i] = t0[iv]*loop_tangent0[i]+
      t1[iv]*loop_tangent1[i]+
      s[iv]*loop_normal[i];
    max_dist = max(max_dist,DIST(xloop[iv],xloop_check));
  }
  // cout << "max projection error : " << max_dist << endl;
  double signed_area = polygonSignedArea(nv, t0, t1);

  // reverse loop to make it ccw
  int *loop_ind;
  if (signed_area < 0.0) {
    loop_ind = new int[nv];
    double* t0_ = new double[nv];
    double* t1_ = new double[nv];
    for (int iv = 0; iv < nv; ++iv) {
      t0_[nv-1-iv] = t0[iv];
      t1_[nv-1-iv] = t1[iv];
      loop_ind[iv] = nv-1-iv;
    }
    for (int iv = 0; iv < nv; ++iv) {
      t0[iv] = t0_[iv];
      t1[iv] = t1_[iv];
    }
    delete[] t0_;
    delete[] t1_;
  }

  // triangulate projected points (ccw ordering)
  // tris[3*(n-2)]
  int* tris = new int[3*(nv-2)];
  const bool triangulated = polygonTriangulate(tris, nv, t0, t1);
  if (!triangulated) {
    cout << " > could not triangulate loop" << endl;
    return;
    delete[] isp_loop;
    delete[] xloop;
    delete[] t0;
    delete[] t1;
    delete[] s;
    delete[] tris;
  }
  int nt = (nv-2);
  // cout << "nv and nt: " << nv << " and " << nt << endl;

  // put indices back in original order
  if (signed_area < 0.0) {
    for (int it = 0; it < nt; ++it) {
      FOR_I3 tris[3*it+i] = loop_ind[tris[3*it+i]];
    }
    delete[] loop_ind;
  }

  double tri_area = 0.0;
  area_vec_2[0] = 0.0;
  area_vec_2[1] = 0.0;
  area_vec_2[2] = 0.0;
  for (int it = 0; it < nt; ++it) {
    double v0[3],v1[3],v2[3];
    FOR_I3 v0[i] = xloop[tris[3*it+0]][i];
    FOR_I3 v1[i] = xloop[tris[3*it+1]][i];
    FOR_I3 v2[i] = xloop[tris[3*it+2]][i];
    double p0[3] = DIFF(v2,v0);
    double p1[3] = DIFF(v1,v0);
    area_vec_2[0] = CROSS_PRODUCT_0(p0,p1);
    area_vec_2[1] = CROSS_PRODUCT_1(p0,p1);
    area_vec_2[2] = CROSS_PRODUCT_2(p0,p1);
    tri_area += 0.5*MAG(area_vec_2);
  }
  //cout << "triangulated area : " << area << endl;
  // cout << "area error : " << tri_area-area << endl;

  // write triangles in OFF file to view

  // ========================================
  // set the simple surface stuff...
  // ========================================

  const int ss_nst0 = nst;
  map<const pair<int,int>,int> edgeMap;
  vector< int > edgeVec;
  edgeVec.resize(nt*3); // index represents my edge id, storage contains my edge neighbor
  const int spoed_new_size0 = spoed_new.size();
  // cout << "ear clipping n new tris (pre): " << spoed_new_size0 << endl;
  spoed_new.resize(spoed_new_size0+nt);
  int outside = 0, inside = 0;
  for (int ist = ss_nst0, nst_new = (nst+nt); ist < nst_new; ++ist) {
    // cout << " > adding info to new tri index " << (ist-ss_nst0)+spoed_new_size0 << endl;
    FOR_I3 {
      int iv1 = tris[(ist-ss_nst0)*3+i]; // local node index
      int iv = tris[(ist-ss_nst0)*3+(i+1)%3]; // local node index + 1 (cw)
      int diff_iv = min(abs(iv1-iv),abs(nv+iv1-iv));
      int ie = (ist-ss_nst0)*3+i; // local edge index (add ss_nst0*3 to make ied)
      spoed_new[(ist-ss_nst0)+spoed_new_size0].spost[i] = isp_loop[iv1];
      // external edges are consecutive
      // if tris[ist*3+i] == j then tris[ist*3+(i+1)%3] == j+1
      if (diff_iv == 1) {
        edgeVec[ie] = -2; // just a flag to skip
        // int ist_nbr = ied_nbr/3;
        // int i_nbr = ied_nbr-ist_nbr*3;
        outside++;
      }
      else {
        map<const pair<int,int>,int>::iterator iter = edgeMap.find(pair<int,int>(iv,iv1));
        // find internal edge in map
        if (iter != edgeMap.end()) {
          assert(edgeVec[iter->second] == -1);
          edgeVec[iter->second] = ie;
          edgeVec[ie] = iter->second;
          // found a match...
          edgeMap.erase(iter);
        }
        // add reverse edge to map
        else {
          edgeMap[pair<int,int>(iv1,iv)] = ie;
          edgeVec[ie] = -1;
        }
      }
      if (edgeVec[ie] == -1) inside++;
    }
    // znost,szost for new tri
    spoed_new[(ist-ss_nst0)+spoed_new_size0].znost = zoneVec.size();
    spoed_new[(ist-ss_nst0)+spoed_new_size0].szost = nsz;
  }

  // internal edges should be matched and we should have nv/ned outside edges
  if (edgeMap.size() != 0 || outside != nv) {
    // remove the new edges
    COUT2(" > triangulation was topologically inconsistent; reverting to original open edge group");
    spoed_new.resize(spoed_new_size0);
  }
  else {
    assert(outside == nv);
    COUT2(" > finished closing loop");
  }

  // cleanup
  edgeVec.clear();
  delete[] isp_loop;
  delete[] xloop;
  delete[] t0;
  delete[] t1;
  delete[] s;
  delete[] tris;

}

bool SimpleSurface::triangulateMarchingFront(vector<NewTri>& newTrisVec,set< std::pair<int,int> >& frontSet,const double (*x)[2],const int n) {

  // zone, subzone for new tris
  const int izone = zoneVec.size();
  // expected number of new tris is nEdges-2 if all are used
  const int n_expected = frontSet.size()-2;
  const int percent_10 = int(floor(double(n_expected)/10.0));
  int next_output = percent_10;
  if (cti_verbose) {
    cout << " > marching front progress: 0%";
    cout.flush();
  }
  // int count[n];
  // for (int i = 0; i < n; ++i)
  // count[i] = 0;

  IntFlag count(n);
  count.setAll(0);

  // int failed[n];
  // for (int i = 0; i < n; ++i)
  // failed[i] = 0;

  IntFlag failed(n);
  failed.setAll(0);

  int index = 0;

  for (set< std::pair<int,int> >::iterator iter = frontSet.begin(); iter != frontSet.end(); ++iter) {
    int i0 = iter->first; assert((i0 >= 0)&&(i0 < n));
    int i1 = iter->second; assert((i1 >= 0)&&(i1 < n));

    count[i0] += 1;
    count[i1] += 1;
  }
  for (int i = 0; i < n; ++i) {
    if (count[i] != 2) {
      // cout << "count[i] < 2: " << count[i] << " debug_count=" << debug_count << endl;
      // assert(0);
      CWARN("some nodes on the passed loops have a valence != 2. Cannot cap these edges; skipping.");
      return false;
    }
  }

  while (!frontSet.empty()) {

    // take the first element from the front...
    set< std::pair<int,int> >::iterator iter = frontSet.begin();
    const int i0 = iter->first; assert((i0 >= 0)&&(i0 < n)); assert(count[i0] > 0);
    const int i1 = iter->second; assert((i1 >= 0)&&(i1 < n)); assert(count[i1] > 0);

    // and remove the element...

    frontSet.erase(iter);
    --count[i0];
    --count[i1];

    ++index;

    // look for a matching edge pair...

    set< std::pair<int,int> >::iterator iter0 = frontSet.find( std::pair<int,int>(i1,i0) );
    if ( (iter0 != frontSet.end()) && ((count[i0] == 1)||(count[i1] == 1)) ) {

      frontSet.erase(iter0);
      --count[i0];
      --count[i1];

      if (percent_10 && (int(newTrisVec.size()) == next_output)) {
        if (cti_verbose) {
          cout << ".." << ceil(double(next_output)/double(n_expected)*100.0) << "\"";
          cout.flush();
        }

        next_output += percent_10;
      }
      newTrisVec.push_back(NewTri(i0,i1,i0,izone,nsz));

      CWARN("Adding a colapsed tri: " << i0+1 << " " << i1+1 << " " << i0+1);
    }
    else {

      int i2;
      for (i2 = 0; i2 < n; ++i2) if ((i2 != i0)&&(i2 != i1)&&(count[i2] > 0)&&(failed[i2] < index)) {

          // check orientation of node...
          const double dx02[2] = { x[i2][0]-x[i0][0], x[i2][1]-x[i0][1] };
          const double dx21[2] = { x[i1][0]-x[i2][0], x[i1][1]-x[i2][1] };
          if ((dx02[0]*dx21[1] - dx02[1]*dx21[0]) < 0.0) {

            // check that no other valid node is in-circle...
            int i;
            for (i = 0; i < n; ++i) if ((i != i0)&&(i != i1)&&(i != i2)&&(count[i] > 0)&&(failed[i] < index)) {
                const double dx00 = x[i0][0] - x[i][0];
                const double dx01 = x[i0][1] - x[i][1];
                const double dx10 = x[i1][0] - x[i][0];
                const double dx11 = x[i1][1] - x[i][1];
                const double dx20 = x[i2][0] - x[i][0];
                const double dx21 = x[i2][1] - x[i][1];
                const double det  =
                  dx00* ( dx11*(dx20*dx20 + dx21*dx21) - dx21*( dx10*dx10 + dx11*dx11)) -
                  dx01* ( dx10*(dx20*dx20 + dx21*dx21) - dx20*( dx10*dx10 + dx11*dx11)) +
                  (dx00*dx00 + dx01*dx01)*( dx10*dx21 - dx11*dx20) ;

                if (det > 0.0) {

                  // cout << " > > potentially failed in-circle test, i: " << i+1 << " det: " << det << " i0: " << i0+1 << " i1: " << i1+1 << endl;

                  // point i is in-circle, but if tri i0->i1->i is invalid, then this is not
                  // a problem...

                  const double dx0i[2] = { x[i][0]-x[i0][0], x[i][1]-x[i0][1] };
                  const double dxi1[2] = { x[i1][0]-x[i][0], x[i1][1]-x[i][1] };
                  if ((dx0i[0]*dxi1[1] - dx0i[1]*dxi1[0]) >= 0.0) {
                    // cout << " > > potential check: failed orientation test, i: " << i+1 << endl;
                    failed[i] = index;
                    continue;
                  }

                  // cout << "1 looking for front: " << i+1 << " " << i0+1 << endl;

                  // new edge i0->i first...
                  set< std::pair<int,int> >::iterator iter0 = frontSet.find( std::pair<int,int>(i,i0) );
                  if (iter0 == frontSet.end()) {

                    // cout << "1 did not find: " << i+1 << " " << i0+1 << endl;

                    // we did not find the front edge i->i0, so we need to ensure that new
                    // front edge i0->i does NOT intersect any other front edges...
                    set< std::pair<int,int> >::iterator iter_check;
                    for (iter_check = frontSet.begin(); iter_check != frontSet.end(); ++iter_check) {
                      const int i0_check = iter_check->first; assert((i0_check >= 0)&&(i0_check < n)); assert(count[i0_check] > 0);
                      const int i1_check = iter_check->second; assert((i1_check >= 0)&&(i1_check < n)); assert(count[i1_check] > 0);

                      // cout << "1 doing edge check: " << i0 << " " << i << " against " << i0_check << " " << i1_check << endl;

                      if ((i0_check != i0)&&(i0_check != i)&&(i1_check != i0)&&(i1_check != i)) {
                        const double dx0_check[2] = { x[i0_check][0]-x[i0][0], x[i0_check][1]-x[i0][1] };
                        const double dx01_check[2] = { x[i1_check][0]-x[i0_check][0], x[i1_check][1]-x[i0_check][1] };
                        const double denom = dx01_check[0]*dx0i[1] - dx01_check[1]*dx0i[0];
                        if (denom != 0.0) {
                          const double s = (dx01_check[0]*dx0_check[1] - dx01_check[1]*dx0_check[0])/denom;
                          if ((s >= 0.0)&&(s <= 1.0)) {
                            // the first line intersected. Now check the second...
                            const double t = (dx0i[0]*dx0_check[1] - dx0i[1]*dx0_check[0])/denom;
                            if ((t >= 0.0)&&(t <= 1.0)) {
                              // cout << " > > potential check: got intersection with edge: " << i0_check+1 << " " << i1_check+1 << endl;
                              failed[i] = index;
                              break;
                            }
                          }
                        }
                      }
                    }
                    if (iter_check != frontSet.end()) {
                      assert(failed[i] == index);
                      // cout << " > > failed in-circle ignored, i: " << i+1 << endl;
                      continue;
                    }

                  }

                  // cout << "2 looking for front: " << i1+1 << " " << i+1 << endl;

                  // new edge i->i1 next...
                  set< std::pair<int,int> >::iterator iter1 = frontSet.find( std::pair<int,int>(i1,i) );
                  if (iter1 == frontSet.end()) {

                    // cout << "2 did not find: " << i1+1 << " " << i+1 << endl;

                    // we did not find the front edge i1->i, so we need to ensure that new
                    // front edge i->i1 does NOT intersect any other front edges...
                    set< std::pair<int,int> >::iterator iter_check;
                    for (iter_check = frontSet.begin(); iter_check != frontSet.end(); ++iter_check) {
                      const int i0_check = iter_check->first; assert((i0_check >= 0)&&(i0_check < n)); assert(count[i0_check] > 0);
                      const int i1_check = iter_check->second; assert((i1_check >= 0)&&(i1_check < n)); assert(count[i1_check] > 0);

                      // cout << "2 doing edge check: " << i << " " << i1 << " against " << i0_check << " " << i1_check << endl;

                      if ((i0_check != i1)&&(i0_check != i)&&(i1_check != i1)&&(i1_check != i)) {
                        const double dx0_check[2] = { x[i0_check][0]-x[i][0], x[i0_check][1]-x[i][1] };
                        const double dx01_check[2] = { x[i1_check][0]-x[i0_check][0], x[i1_check][1]-x[i0_check][1] };
                        const double denom = dx01_check[0]*dxi1[1] - dx01_check[1]*dxi1[0];
                        if (denom != 0.0) {
                          const double s = (dx01_check[0]*dx0_check[1] - dx01_check[1]*dx0_check[0])/denom;
                          if ((s >= 0.0)&&(s <= 1.0)) {
                            // the first line intersected. Now check the second...
                            const double t = (dxi1[0]*dx0_check[1] - dxi1[1]*dx0_check[0])/denom;
                            if ((t >= 0.0)&&(t <= 1.0)) {
                              // if (debug) cout << " > > potential check: got intersection with edge: " << i0_check+1 << " " << i1_check+1 << endl;
                              failed[i] = index;
                              break;
                            }
                          }
                        }
                      }
                    }
                    if (iter_check != frontSet.end()) {
                      assert(failed[i] == index);
                      // if (debug) cout << " > > failed in-circle ignored, i: " << i+1 << endl;
                      continue;
                    }
                  }

                  break;

                }
              }

            // if we did not get through all the points, then i2 will not work, so continue...
            if (i < n) {
              // cout << " > > failed in-circle test, i: " << i+1 << endl;
              failed[i2] = index;
              continue;
            }

            // we are considering forming the tri i0->i1->i2...

            // new edge i0->i2 first...
            set< std::pair<int,int> >::iterator iter0 = frontSet.find( std::pair<int,int>(i2,i0) );
            if (iter0 == frontSet.end()) {
              // we did not find the front edge i2->i0, so we need to ensure that new
              // front edge i0->i2 does NOT intersect any other front edges...
              // we already have dx02...
              set< std::pair<int,int> >::iterator iter_check;
              for (iter_check = frontSet.begin(); iter_check != frontSet.end(); ++iter_check) {
                const int i0_check = iter_check->first; assert((i0_check >= 0)&&(i0_check < n)); assert(count[i0_check] > 0);
                const int i1_check = iter_check->second; assert((i1_check >= 0)&&(i1_check < n)); assert(count[i1_check] > 0);
                if ((i0_check != i0)&&(i0_check != i2)&&(i1_check != i0)&&(i1_check != i2)) {
                  const double dx0_check[2] = { x[i0_check][0]-x[i0][0], x[i0_check][1]-x[i0][1] };
                  const double dx01_check[2] = { x[i1_check][0]-x[i0_check][0], x[i1_check][1]-x[i0_check][1] };
                  const double denom = dx01_check[0]*dx02[1] - dx01_check[1]*dx02[0];
                  if (denom != 0.0) {
                    const double s = (dx01_check[0]*dx0_check[1] - dx01_check[1]*dx0_check[0])/denom;
                    if ((s >= 0.0)&&(s <= 1.0)) {
                      // the first line intersected. Now check the second...
                      const double t = (dx02[0]*dx0_check[1] - dx02[1]*dx0_check[0])/denom;
                      if ((t >= 0.0)&&(t <= 1.0)) {
                        // cout << " > > got intersection with edge: " << i0_check+1 << " " << i1_check+1 << endl;
                        failed[i2] = index;
                        break;
                      }
                    }
                  }
                }
              }
              if (iter_check != frontSet.end())
                continue;
            }

            // new edge i2->i1 next...
            set< std::pair<int,int> >::iterator iter1 = frontSet.find( std::pair<int,int>(i1,i2) );
            if (iter1 == frontSet.end()) {
              // we did not find the front edge i1->i2, so we need to ensure that new
              // front edge i2->i1 does NOT intersect any other front edges...
              set< std::pair<int,int> >::iterator iter_check;
              for (iter_check = frontSet.begin(); iter_check != frontSet.end(); ++iter_check) {
                const int i0_check = iter_check->first; assert((i0_check >= 0)&&(i0_check < n)); assert(count[i0_check] > 0);
                const int i1_check = iter_check->second; assert((i1_check >= 0)&&(i1_check < n)); assert(count[i1_check] > 0);
                if ((i0_check != i1)&&(i0_check != i2)&&(i1_check != i1)&&(i1_check != i2)) {
                  const double dx0_check[2] = { x[i0_check][0]-x[i2][0], x[i0_check][1]-x[i2][1] };
                  const double dx01_check[2] = { x[i1_check][0]-x[i0_check][0], x[i1_check][1]-x[i0_check][1] };
                  const double denom = dx01_check[0]*dx21[1] - dx01_check[1]*dx21[0];
                  if (denom != 0.0) {
                    const double s = (dx01_check[0]*dx0_check[1] - dx01_check[1]*dx0_check[0])/denom;
                    if ((s >= 0.0)&&(s <= 1.0)) {
                      // the first line intersected. Now check the second...
                      const double t = (dx21[0]*dx0_check[1] - dx21[1]*dx0_check[0])/denom;
                      if ((t >= 0.0)&&(t <= 1.0)) {
                        // cout << " > > got intersection with edge: " << i0_check+1 << " " << i1_check+1 << endl;
                        failed[i2] = index;
                        break;
                      }
                    }
                  }
                }
              }
              if (iter_check != frontSet.end())
                continue;

            }

            // if we made it here, then i2 is our guy!...
            if (percent_10 && (int(newTrisVec.size()) == next_output)) {
              if (cti_verbose) {
                cout << ".." << ceil(double(next_output)/double(n_expected)*100.0) << "%";
                cout.flush();
              }

              next_output += percent_10;
            }
            newTrisVec.push_back(NewTri(i0,i1,i2,izone,nsz));

            // delete the current front element from the front...
            if (iter0 == frontSet.end()) {
              frontSet.insert( std::pair<int,int>(i0,i2) );
              ++count[i0];
              ++count[i2];
            }
            else {
              frontSet.erase(iter0);
              --count[i0];
              --count[i2];
            }

            if (iter1 == frontSet.end()) {
              frontSet.insert( std::pair<int,int>(i2,i1) );
              ++count[i1];
              ++count[i2];
            }
            else {
              frontSet.erase(iter1);
              --count[i1];
              --count[i2];
            }

            break;

          }
        }

      if (i2 == n) {
        cout << endl << "WARNING triangulate marching front: no valid i2" << endl;
        cout.flush();
        // assert(0);
        return false;
      }
    }
  }

  // front is gone. Check that all nodes are cleared...
  if (cti_verbose) cout << endl;

  // int ierr = 0;
  for (int i = 0; i < n; ++i) {
    if (count[i] != 0) {
      CWARN("front is finished but some nodes were not processed");
      return false;
    }
  }
  return true;
  // assert(ierr == 0);

}



void SimpleSurface::createPlanarOpenEdgeLoopTris(vector<NewTri>& newTrisVec,const vector<int>& group_indices) {
  // grab relevant open edges...

  vector<pair<int8,char> > closeTeostVec;
  int igr = -1;
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      if (getOpenEdgeGroup(igr,ist,i)) {
        std::vector<int>::const_iterator it = find(group_indices.begin(), group_indices.end(), igr);
        if (it != group_indices.end()) closeTeostVec.push_back(pair<int8,char>(ist,i)); // save this edge + triangle
      }
    }
  }

  assert(!closeTeostVec.empty());

  createPlanarEdgeLoopTris(newTrisVec,closeTeostVec);
}

void SimpleSurface::createPlanarHalfEdgeLoopTris(vector<NewTri>& newTrisVec, const vector<int>& group_indices) {
  // grab relevant open edges...

  vector<pair<int8,char> > closeTeostVec;
  int igr = -1;
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      if (getHalfEdgeGroup(igr,ist,i)) {
        std::vector<int>::const_iterator it = find(group_indices.begin(), group_indices.end(), igr);
        if (it != group_indices.end()) closeTeostVec.push_back(pair<int8,char>(ist,i)); // save this edge + triangle
      }
    }
  }

  assert(!closeTeostVec.empty());

  createPlanarEdgeLoopTris(newTrisVec,closeTeostVec);
}

void SimpleSurface::createPlanarEdgeLoopTris(vector<NewTri>& newTrisVec, vector<pair<int8,char> >& closeTeostVec) {
  // start by confirming that the nodes associated with these edges
  // all have one edge in and one edge out...

  sp_flag.setLength(nsp);
  sp_flag.setAll(0);

  for (int ii = 0, ii_end=closeTeostVec.size(); ii < ii_end; ++ii) {
    const int ist = closeTeostVec[ii].first;
    const int i  = closeTeostVec[ii].second;
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    assert((sp_flag[isp0]&1) == 0);
    sp_flag[isp0] |= 1;
    assert((sp_flag[isp1]&2) == 0);
    sp_flag[isp1] |= 2;
  }

  int nsp_ = 0;
  double xp[3] = { 0.0, 0.0, 0.0 }; // centroid
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] == 3) {
      sp_flag[isp] = nsp_++;
      FOR_I3 xp[i] += xsp[isp][i];
    }
    else {
      assert(sp_flag[isp] == 0);
      sp_flag[isp] = -1;
    }
  }

  FOR_I3 xp[i] /= double(nsp_);

  // good -- all nodes are touched by two edges, and we have the centroid.
  // Now build the plane for projection...

  double xx = 0.0;
  double xy = 0.0;
  double xz = 0.0;
  double yy = 0.0;
  double yz = 0.0;
  double zz = 0.0;

  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] >= 0) {
      double dx[3] = DIFF(xsp[isp],xp);
      xx += dx[0]*dx[0];
      xy += dx[0]*dx[1];
      xz += dx[0]*dx[2];
      yy += dx[1]*dx[1];
      yz += dx[1]*dx[2];
      zz += dx[2]*dx[2];
    }
  }

  double det_x = yy*zz - yz*yz;
  double det_y = xx*zz - xz*xz;
  double det_z = xx*yy - xy*xy;

  double np[3];
  if (fabs(det_x) > max(fabs(det_y),fabs(det_z))) {
    np[0] = det_x;
    np[1] = xz*yz - xy*zz;
    np[2] = xy*yz - xz*yy;
  }
  else if (fabs(det_y) > fabs(det_z)) {
    np[0] = xz*yz - xy*zz;
    np[1] = det_y;
    np[2] = xy*xz - yz*xx;
  }
  else {
    np[0] = xy*yz - xz*yy;
    np[1] = xy*xz - yz*xx;
    np[2] = det_z;
  }

  NORMALIZE(np);

  COUT2("    > planar closing info: nsp " << nsp_ << "; centroid " << COUT_VEC(xp) << "; normal " << COUT_VEC(np) );

  // build 2 vectors in this plane and set the coordinated for each point...

  double up[3];
  if (fabs(np[0]) < min(fabs(np[1]),fabs(np[2]))) {
    up[0] = 1.0;
    up[1] = 0.0;
    up[2] = 0.0;
  }
  else if (fabs(np[1]) < min(fabs(np[0]),fabs(np[2]))) {
    up[0] = 0.0;
    up[1] = 1.0;
    up[2] = 0.0;
  }
  else {
    up[0] = 0.0;
    up[1] = 0.0;
    up[2] = 1.0;
  }
  double e0[3] = CROSS_PRODUCT(up,np);
  NORMALIZE(e0);
  double e1[3] = CROSS_PRODUCT(np,e0);
  NORMALIZE(e1);

  double (*xsp_)[2] = new double[nsp_][2];
  IntFlag isp_local_to_global(nsp_);
  isp_local_to_global.setAll(-1);
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] >= 0) {
      int isp_ = sp_flag[isp];
      isp_local_to_global[isp_] = isp;
      xsp_[isp_][0] =
        (xsp[isp][0]-xp[0])*e0[0] +
        (xsp[isp][1]-xp[1])*e0[1] +
        (xsp[isp][2]-xp[2])*e0[2];
      xsp_[isp_][1] =
        (xsp[isp][0]-xp[0])*e1[0] +
        (xsp[isp][1]-xp[1])*e1[1] +
        (xsp[isp][2]-xp[2])*e1[2];
      //cout << " XXXXX " << xsp_[isp_][0] << " " << xsp_[isp_][1] << endl;
    }
  }
  assert(isp_local_to_global.countNegative() == 0);

  // now get imputs ready for triangulation algo...


  set<pair<int,int> > frontSet;
  for (int ii = 0, ii_end=closeTeostVec.size(); ii < ii_end; ++ii) {
    const int ist = closeTeostVec[ii].first;
    const int i  = closeTeostVec[ii].second;
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    frontSet.insert(pair<int,int>(sp_flag[isp0],sp_flag[isp1]));
  }


  vector<NewTri> _newTrisVec;
  bool success = triangulateMarchingFront(_newTrisVec,frontSet,xsp_,nsp_);
  DELETE(xsp_);
  if (!success || _newTrisVec.empty()) {
    CWARN("there was a problem with the marching front; skipping closing of these holes");
    _newTrisVec.clear();
    return;
  }
  else {
    for (int ist = 0, nst_new = _newTrisVec.size(); ist < nst_new; ++ist) {
      // convert back to global isps while mapping is available
      _newTrisVec[ist].spost[0] = isp_local_to_global[_newTrisVec[ist].spost[0]];
      _newTrisVec[ist].spost[1] = isp_local_to_global[_newTrisVec[ist].spost[1]];
      _newTrisVec[ist].spost[2] = isp_local_to_global[_newTrisVec[ist].spost[2]];
    }
    newTrisVec.insert(newTrisVec.end(),_newTrisVec.begin(),_newTrisVec.end());
  }
}

void SimpleSurface::closePlanarHalfEdgeLoops(vector<int>& group_indices,vector<HalfEdgeGroupData>& halfEdgeGroupDataVec,const bool nested, const double nested_tol, const string zonename,const bool b_open) {

  ensureTeost();

  vector<NewTri> newTrisVec;

  if (nested) {
    vector<pair<int,int> > gr_pairs;
    findNestedHalfEdgeLoops(gr_pairs,nested_tol,group_indices,halfEdgeGroupDataVec);

    const int n_nested_groups = gr_pairs.size();
    COUT1(" > found " << n_nested_groups << " pairs of nested loops to close");

    // close pairs of holes using fmm planar closing
    vector<int> nestedIndices(2);
    int count = 0;
    int success = 0;
    int nst_new = newTrisVec.size();
    for (vector<pair<int,int> >::iterator it=gr_pairs.begin(); it!=gr_pairs.end(); ++it) {
      nestedIndices[0] = it->first;
      nestedIndices[1] = it->second;
      COUT2("    > closing loop pair " << (count+1) << "/" << n_nested_groups);
      if (b_open) createPlanarOpenEdgeLoopTris(newTrisVec,nestedIndices);
      else createPlanarHalfEdgeLoopTris(newTrisVec,nestedIndices);
      if (nst_new < int(newTrisVec.size())) {
        ++success;
        nst_new = newTrisVec.size();
      }
      ++count;
    }

    COUT2(" > successfully closed " << success << "/"<< count << " pairs of concentric loops");
  }
  else {
    if (b_open) createPlanarOpenEdgeLoopTris(newTrisVec,group_indices);
    else createPlanarHalfEdgeLoopTris(newTrisVec,group_indices);
  }

  if (newTrisVec.empty()) {
    COUT2(" > no new tris were created");
    return;
  }
  COUT2(" > new tris being added: " << newTrisVec.size());

  // dump as tecplot...
  // {
  //
  //   FILE * fp = fopen("tris.dat","w");
  //   assert(fp != NULL);
  //   fprintf(fp,"TITLE = \"tris\"\n");
  //   fprintf(fp,"VARIABLES = \"X\"\n");
  //   fprintf(fp,"\"Y\"\n");
  //   fprintf(fp,"\"Z\"\n");
  //   // zone header
  //   fprintf(fp,"ZONE T=\"blah\"\n");
  //   fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nsp_,triVec.size());
  //
  //   // order should be good...
  //   for (int isp = 0; isp < nsp; ++isp) {
  //     if (sp_flag[isp] >= 0) {
  //       fprintf(fp,"%lf %lf %lf\n",xsp[isp][0],xsp[isp][1],xsp[isp][2]);
  //     }
  //   }
  //
  //   for (int ii = 0; ii < triVec.size(); ++ii) {
  //     fprintf(fp,"%d %d %d\n",triVec[ii].ip[0]+1,triVec[ii].ip[1]+1,triVec[ii].ip[2]+1);
  //   }
  //
  //   fclose(fp);
  //
  // }

  // add tris to surface
  const int ss_nst0 = nst;
  nst += newTrisVec.size();
  szost.resize(nst);

  growNstData(nst,ss_nst0);

  int count = 0;
  for (int ist = 0, nst_new = newTrisVec.size(); ist < nst_new; ++ist) {
    // marching front flips the normal alignment; do the flipping here when adding tris
    spost[ist+ss_nst0][0] =newTrisVec[ist].spost[0];
    spost[ist+ss_nst0][2] =newTrisVec[ist].spost[1];
    spost[ist+ss_nst0][1] =newTrisVec[ist].spost[2];
    znost[ist+ss_nst0] = newTrisVec[ist].znost;  // assumes going to a new zone
    szost[ist+ss_nst0] = newTrisVec[ist].szost;
    ++count;
  }
  assert(count+ss_nst0 == nst);

  const int zone_index = getZoneIndex(zonename);
  if (zone_index == -1) {
    stringstream oss; oss << n_closed_edge_groups++;
    string zone_name = "Closed_Loops_" + oss.str();
    if (zonename != "") zone_name = zonename;
    addNewZone(zone_name); // also resizes szozn_i
    szozn_i[zoneVec.size()] += 1; // add sub-zone in new zone
    ++nsz;  // add a new sz
  }
  else {
    // moving new tris to existing zone
    const int new_sz = szozn_i[zone_index+1] - szozn_i[zone_index];  // number of subzones in this zone
    globalSzostToLocal();
    for (int ist = 0, nst_new = newTrisVec.size(); ist < nst_new; ++ist) {
      znost[ist+ss_nst0] = zone_index;
      szost[ist+ss_nst0] = new_sz;
      ++count;
    }
    buildSzoznFromLocalSzost();  // computes nsz here, so don't need to increment manually
    localSzostToGlobal();
  }

  // teost and (implicitly) open edge groups need to be updated
  clearTeost();
  if (b_open) clearOpenEdgeGroups();
  else clearDynamicEdgeGroups();

  // addition of new zones requires rebuilding zone/subzone info
  clearSubzoneData();
  clearZoneData();
}

void SimpleSurface::flipFlaggedTris() {

  bool changed = false;
  FOR_IST {
    // st_flag is boolean which determines whether tri gets flipped or not
    if (st_flag[ist]) {
      changed = true;
      // leave 0 the same and flip nodes associated with spost 1,2...
      int tmp = spost[ist][1];
      spost[ist][1] = spost[ist][2];
      spost[ist][2] = tmp;
    }
  }

  if (changed) {
    // need to rebuild edge structures
    clearTeost();
    clearDynamicEdgeGroups();

    // metadata for zones needs to update
    clearSubzoneData();
    clearZoneData();
  }
}

void SimpleSurface::separateOverlappingTris(const double dn) {

  ensureTeost();

  // flag nodes touching overlapping tris and unflag ones that are along crease or at border...
  sp_flag.resize(nsp);
  sp_flag.setAll(0);
  FOR_IST {
    if (nmTriData.isOverlapping(ist)) {
      FOR_I3 {
        int ist_nbr,i_nbr,orient;
        if (getTriNbrData(ist_nbr,i_nbr,orient,ist,i)) {
          if (nmTriData.isOverlapping(ist_nbr)) {
            const double dp = MiscUtils::getCosAngleBetweenTris(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],
                                                                xsp[spost[ist_nbr][0]],xsp[spost[ist_nbr][1]],xsp[spost[ist_nbr][2]]);
            if (dp >= 0.0) {
              sp_flag[spost[ist][i]] = 1;
              sp_flag[spost[ist][(i+1)%3]] = 1;
            }
          }
        }
      }
    }
  }

  // move flagged nodes in local node normal direction...
  double (*n_sp)[3] = new double[nsp][3];
  FOR_ISP FOR_I3 n_sp[isp][i] = 0.0;
  FOR_IST {
    const double n[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    FOR_I3 {
      const int isp = spost[ist][i];
      FOR_J3 n_sp[isp][j] += n[j];
    }
  }
  FOR_ISP {
    const double mag = MAG(n_sp[isp]);
    assert(mag > 0.0);
    FOR_I3 xsp[isp][i] -= n_sp[isp][i]*dn/mag;
  }
  delete[] n_sp;

}

void SimpleSurface::buildNonManifoldData() {
  COUT2("SimpleSurface::buildNonManifoldData()");
  diagnoseManifold();
}

void SimpleSurface::ensureNonManifoldData() {
  if (nmTriData.isNull()) buildNonManifoldData();
}

bool SimpleSurface::hasNonManifoldData() {
  return (!nmTriData.isNull());
}

void SimpleSurface::clearNonManifoldData() {
  COUT2("SimpleSurface::clearNonManifoldData()");
  if ( !nmTriData.isNull() ) nmTriData.clear();
}

void SimpleSurface::setNonManifoldCheckProperties(const double dist_tol,const double angle_tol_degrees,const bool check_self_intersections,const string output_type,const int n_samples,const int nbr_layers) {
  // overwrites default values from nmTriData class instanciation
  nmTriData.setDistTol(dist_tol);
  nmTriData.setAngleTol(angle_tol_degrees);
  nmTriData.setCheckSelfIntersections(check_self_intersections);
  nmTriData.setFormat(output_type);
  nmTriData.setNSeedSamples(n_samples);
  nmTriData.setNbrLayers(nbr_layers);
}

void SimpleSurface::diagnoseManifold() {

  // no need to ensure NonManifoldData b/c this is what builds it. so regardless we will overwrite
  nmTriData.resize(nst);
  nmTriData.setAll(0);

  ensureTeost();
  // check for linear tris; both:
  // 1: node-collapsed tris | 2: zero-area tris
  int ime, ime_orientation;
  int ist_nbr[3];
  for (int ist=0; ist<nst; ++ist) {
    const int isp0 = spost[ist][0];
    const int isp1 = spost[ist][1];
    const int isp2 = spost[ist][2];
    if ( (isp0 == isp1) || (isp0 == isp2) || (isp1 == isp2) ) nmTriData.setCollapsed(ist);  // node collapsed

    const double n2[3] = TRI_NORMAL_2(xsp[isp0],xsp[isp1],xsp[isp2]);
    const double mag_n2 = MAG(n2);
    if ( !(mag_n2 > 0.0) ) nmTriData.setLinear(ist);  // zero surface area

    getTriNbrs(ist_nbr,ist);
    if ( (ist_nbr[0] == ist_nbr[1]) || (ist_nbr[0] == ist_nbr[2]) ) {
      if (ist_nbr[0] >= 0) nmTriData.setMultiNeighbor(ist);
    }
    else if ( (ist_nbr[1] == ist_nbr[2]) && (ist_nbr[1] >= 0) ) nmTriData.setMultiNeighbor(ist);

    FOR_I3 {
      if (isEdgeMulti(ime,ime_orientation,ist,i)) {
        nmTriData.setNmeAdjacent(ist);
      }
    }
  }

  nmTriData.getCounts();  // compute counts of each type

  UIMessage_ msg(MESSAGE); // TODO: replace
  msg.addLine("");
  msg.addLine("----- Non-manifold Element Report -----");
  msg.addLine("");
  msg.addLine(stringstream().flush() << " > detected " << nmTriData.countCollapsed() << " node-collapsed tris");
  if (nmTriData.countCollapsed()) {
    int zone_counts[zoneVec.size()];
    for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) zone_counts[izn] = 0;
    for (int ist=0; ist<nst; ++ist) {
      if (nmTriData.isCollapsed(ist)) ++zone_counts[znost[ist]];  // increment zone counter
    }
    for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
      if (zone_counts[izn]) msg.addLine(stringstream().flush() << "    > detected " << zone_counts[izn] << " problems in " << zoneVec[izn].getName());
    }
  }
  msg.addLine(stringstream().flush() << " > detected " << nmTriData.countLinear() << " linear, i.e., zero surface area, tris (NOTE: will include collapsed tris in this count)");
  if (nmTriData.countLinear()) {
    int zone_counts[zoneVec.size()];
    for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) zone_counts[izn] = 0;
    for (int ist=0; ist<nst; ++ist) {
      if (nmTriData.isLinear(ist)) ++zone_counts[znost[ist]];  // increment zone counter
    }
    for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
      if (zone_counts[izn]) msg.addLine(stringstream().flush() << "    > detected " << zone_counts[izn] << " problems in " << zoneVec[izn].getName());
    }
  }
  msg.addLine(stringstream().flush() << " > detected " << nmTriData.countNmeAdjacent() << " non-manifold multi-edge adjacent tris");
  if (nmTriData.countNmeAdjacent()) {
    int zone_counts[zoneVec.size()];
    for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) zone_counts[izn] = 0;
    for (int ist=0; ist<nst; ++ist) {
      if (nmTriData.isNmeAdjacent(ist)) ++zone_counts[znost[ist]];  // increment zone counter
    }
    for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
      if (zone_counts[izn]) msg.addLine(stringstream().flush() << "    > detected " << zone_counts[izn] << " problems in " << zoneVec[izn].getName());
    }
  }
  msg.addLine(stringstream().flush() << " > detected " << nmTriData.countMultiNeighbor() << " tris with a repeated neighbor");
  if (nmTriData.countMultiNeighbor()) {
    int zone_counts[zoneVec.size()];
    for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) zone_counts[izn] = 0;
    for (int ist=0; ist<nst; ++ist) {
      if (nmTriData.isMultiNeighbor(ist)) ++zone_counts[znost[ist]];  // increment zone counter
    }
    for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
      if (zone_counts[izn]) msg.addLine(stringstream().flush() << "    > detected " << zone_counts[izn] << " problems in " << zoneVec[izn].getName());
    }
  }

  if (nmTriData.checkSelfIntersections()) {
    searchForAndCategorizeIntersectingTris(nmTriData.distTol(),nmTriData.angleTol());

    nmTriData.getCounts();  // re-compute counts of each type since new types added

    //FOR_IST if (nmTriData.isOverlapping(ist)) cout << ist << ",";
    //cout << endl;
    msg.addLine(stringstream().flush() << " > detected " << nmTriData.countOverlapping() << " overlapping (coplanar) tris");
    if (nmTriData.countOverlapping()) {
      int zone_counts[zoneVec.size()];
      for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) zone_counts[izn] = 0;
      for (int ist=0; ist<nst; ++ist) {
        if (nmTriData.isOverlapping(ist)) ++zone_counts[znost[ist]];  // increment zone counter
      }
      for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
        if (zone_counts[izn]) msg.addLine(stringstream().flush() << "    > detected " << zone_counts[izn] << " problems in " << zoneVec[izn].getName());
      }
    }
    msg.addLine(stringstream().flush() << " > detected " << nmTriData.countImprinting() << " imprinting (node or edge adjacent) tris");
    if (nmTriData.countImprinting()) {
      int zone_counts[zoneVec.size()];
      for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) zone_counts[izn] = 0;
      for (int ist=0; ist<nst; ++ist) {
        if (nmTriData.isImprinting(ist)) ++zone_counts[znost[ist]];  // increment zone counter
      }
      for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
        if (zone_counts[izn]) msg.addLine(stringstream().flush() << "    > detected " << zone_counts[izn] << " problems in " << zoneVec[izn].getName());
      }
    }
    msg.addLine(stringstream().flush() << " > detected " << nmTriData.countIntersecting() << " intersecting (penetrating) tris");
    if (nmTriData.countIntersecting()) {
      int zone_counts[zoneVec.size()];
      for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) zone_counts[izn] = 0;
      for (int ist=0; ist<nst; ++ist) {
        if (nmTriData.isIntersecting(ist)) ++zone_counts[znost[ist]];  // increment zone counter
      }
      for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
        if (zone_counts[izn]) msg.addLine(stringstream().flush() << "    > detected " << zone_counts[izn] << " problems in " << zoneVec[izn].getName());
      }
    }
  }
  msg.dumpLines();
  WebUI::webUIOutput.add_(msg);

  // produce desired output files
  const string output_type = nmTriData.getFormat();
  if (output_type != "none") {
    IntFlag st_show(nst);

    // collapsed and linear need to be edge-drawn...

    // following zones will overwrite each other (steal already moved tris)
    // maybe dump selectively to sbin/tecplot...

    //also manually hack szost for now - don't add anything to szdata
    if (nmTriData.countNmeAdjacent()) {
      flagNonManifoldTrisAndNeighborhood(st_show,&SimpleSurface::isNmeAdjacent,nmTriData.nNbrLayers());
      if (output_type == "tecplot") writeSelectedFacesByZoneToTecplot("nmeAdjacent.dat",st_show,true);
      else if (output_type == "sbin") writeSelectedTrisToBinary("nmeAdjacent.sbin",st_show);
    }

    if (nmTriData.countLinear()) {
      flagNonManifoldTrisAndNeighborhood(st_show,&SimpleSurface::isLinear,nmTriData.nNbrLayers());
      if (output_type == "tecplot") writeSelectedFacesByZoneToTecplot("linear.dat",st_show,true);
      else if (output_type == "sbin") writeSelectedTrisToBinary("linear.sbin",st_show);
    }

    if (nmTriData.countIntersecting()) {
      flagNonManifoldTrisAndNeighborhood(st_show,&SimpleSurface::isIntersecting,nmTriData.nNbrLayers());
      if (output_type == "tecplot") writeSelectedFacesByZoneToTecplot("intersecting.dat",st_show,true);
      else if (output_type == "sbin") writeSelectedTrisToBinary("intersecting.sbin",st_show);
    }

    if (nmTriData.countMultiNeighbor()) {
      flagNonManifoldTrisAndNeighborhood(st_show,&SimpleSurface::isMultiNeighbor,nmTriData.nNbrLayers());
      if (output_type == "tecplot") writeSelectedFacesByZoneToTecplot("multi-neighbor.dat",st_show,true);
      else if (output_type == "sbin") writeSelectedTrisToBinary("multi-neighbor.sbin",st_show);
    }

    if (nmTriData.countOverlapping()) {
      flagNonManifoldTrisAndNeighborhood(st_show,&SimpleSurface::isOverlapping,nmTriData.nNbrLayers());
      if (output_type == "tecplot") writeSelectedFacesByZoneToTecplot("overlapping.dat",st_show,true);
      else if (output_type == "sbin") writeSelectedTrisToBinary("overlapping.sbin",st_show);
    }

    if (nmTriData.countImprinting()) {
      flagNonManifoldTrisAndNeighborhood(st_show,&SimpleSurface::isImprinting,nmTriData.nNbrLayers());
      if (output_type == "tecplot") writeSelectedFacesByZoneToTecplot("imprinting.dat",st_show,true);
      else if (output_type == "sbin") writeSelectedTrisToBinary("imprinting.sbin",st_show);
    }
  }

}

bool SimpleSurface::isMultiNeighbor(const int ist) const {
  return nmTriData.isMultiNeighbor(ist);
}

bool SimpleSurface::isNmeAdjacent(const int ist) const {
  return nmTriData.isNmeAdjacent(ist);
}

bool SimpleSurface::isLinear(const int ist) const {
  return nmTriData.isLinear(ist);
}

bool SimpleSurface::isIntersecting(const int ist) const {
  return nmTriData.isIntersecting(ist);
}

bool SimpleSurface::isOverlapping(const int ist) const {
  return nmTriData.isOverlapping(ist);
}

bool SimpleSurface::isImprinting(const int ist) const {
  return nmTriData.isImprinting(ist);
}

void SimpleSurface::flagNonManifoldTrisAndNeighborhood(IntFlag& st_show,bool (SimpleSurface::*nmtCriterion)(const int) const,const int nbr_layers) const {
  for (int ist=0; ist<nst; ++ist) {
    if ((this->*nmtCriterion)(ist)) st_show[ist] = 1;
    else st_show[ist] = 0;
  }

  if (nbr_layers == 0) return;

  int ist_nbr,i_nbr,orient;
  for (int layer=0; layer<nbr_layers; ++layer) {
    for (int ist=0; ist<nst; ++ist) {
      if (st_show[ist] == (layer+1)) {
        for (int i=0; i<3; ++i) {
          if (getTriNbrData(ist_nbr,i_nbr,orient,ist,i)) {
            if (st_show[ist_nbr] == 0) st_show[ist_nbr] = st_show[ist]+1;
          }
        }
      }
    }
  }
}

void SimpleSurface::flagNonManifoldTris() {
  cout << "SimpleSurface::flagNonManifoldTris()" << endl;
  assert(!nmTriData.isNull());
  st_flag.resize(nst);
  st_flag.setAll(0);
  FOR_IST if (!nmTriData.isClean(ist)) st_flag[ist] = 1;
}

void SimpleSurface::flagNonManifoldTris(const int nme_bits) {
  cout << "SimpleSurface::flagNonManifoldTris()" << endl;
  assert(!nmTriData.isNull());
  st_flag.resize(nst);
  st_flag.setAll(0);
  for (int ist=0; ist<nst; ++ist) {
    if (nmTriData.hasBits(ist,nme_bits)) st_flag[ist] = 1;
  }
}

double SimpleSurface::findMinNonzeroEdgeLength() {
  ensureTeost();
  double minD2=HUGE_VAL;
  bool * tri_edge_flag = new bool[nst*3];
  for (int ist = 0, end=nst*3; ist < end; ++ist) tri_edge_flag[ist] = false;

  for (int ist = 0; ist < nst; ++ist) {
    for (int edge=0; edge<3; ++edge) {
      // only process this edge if it hasn't been set already
      if (!tri_edge_flag[ist*3 + edge]) {
        const double dx2 = DIST2( xsp[spost[ist][edge]],xsp[spost[ist][(edge+1)%3]] );
        if ((dx2<minD2) && (dx2 != 0.0)) minD2 = dx2;

        tri_edge_flag[ist*3 + edge] = true;

        if (!(isEdgeOpen(ist,edge) || isEdgeMulti(ist,edge))) {
          int ist_nbr, edge_nbr, orient;
          getTriNbrData(ist_nbr,edge_nbr,orient,ist,edge);
          tri_edge_flag[ist_nbr*3 + edge_nbr] = true;  // flag neighboring tri edge as computed
        }
      }
    }
  }
  delete[] tri_edge_flag;

  if (minD2 != HUGE_VAL) return sqrt(minD2);
  else return 0.0;  // condition which means an error is present - no non-zero edge lengths
}

void SimpleSurface::searchForAndCategorizeIntersectingTris(const double dist_tol,const double angle_tol_degrees) {
  //TODO anything needed here to restrict memory footprint?

  // tolerances
  // should use a global min-edge length value as distance tolerance
  const double min_edge_l = findMinNonzeroEdgeLength();

  const double angle_tol = cos(M_PI * (1.0-(angle_tol_degrees/180.0)));  // radians ~= 0.57 degrees
  COUT1(" > searching for self-intersections");
  COUT1("    > angle tolerance: " << angle_tol_degrees << " degrees");
  COUT1("    > distance tolerance: " << dist_tol);
  COUT1("    > minimum detected edge length: " << min_edge_l);

  const double start_time = MPI_Wtime();

  // build Adts for all tris
  COUT1("    > building surface-tri Adt: ");
  double (*bbmin)[3] = new double[nst][3];
  double (*bbmax)[3] = new double[nst][3];
  for (int ist = 0; ist < nst; ++ist) {
    // create buffer around tri based on min allowable edge length
    FOR_J3 bbmin[ist][j] = min(xsp[spost[ist][0]][j],min(xsp[spost[ist][1]][j],xsp[spost[ist][2]][j])) - dist_tol;
    FOR_J3 bbmax[ist][j] = max(xsp[spost[ist][0]][j],max(xsp[spost[ist][1]][j],xsp[spost[ist][2]][j])) + dist_tol;
  }
  Adt<double> stAdt(nst,bbmin,bbmax);
  DELETE(bbmin);
  DELETE(bbmax);

  const double adt_time = MPI_Wtime()-start_time;
  COUT1("       > finished in " << adt_time << "s");

  // just want to see which tris intersect/overlap, if already set to yes don't do any further processing
  // use method of Moller 1997b
  vector<int> candidates;
  vector<int>::iterator c_it;
  int ist_nbr[3];
  int shareNode[3][2];  // index of (c_ist,ist) node that is shared
  double signedDistToIst[3];
  IntFlag signToIst(3);
  double signedDistToCIst[3];
  IntFlag signToCIst(3);

  // int _ist_nbr,_i_nbr,_orient;
  const int interval = int(floor(0.1*(double)nst));
  cout << "    > comparing tris in proximity to one another" << endl;
  cout << "      0%" << std::flush;
  int percent_done = 10;
  for (int ist = 0; ist < nst; ++ist) {
    if ((ist > 0) && (ist%interval == 0)) {
      if (mpi_rank == 0) cout << ".." << percent_done << "%" << std::flush;
      percent_done += 10;
    }
    // if already set as intersection/overlapping, no need to check this tri...?
    // if (nmTriData.isIntersecting(ist) || nmTriData.isOverlapping(ist)) continue;

    // store walkable neighbors
    //TODO currently ignores multi-edge-based nbrs
    getTriNbrs(ist_nbr,ist);

    // compute plane equation coeffs (N,d) where N*x + d = 0 for points x on the plane
    double N_ist[3];
    double d_ist;
    double N_c_ist[3];
    double d_c_ist;
    MiscUtils::computePlaneCoeffsFromPoints(N_ist,d_ist,xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);

    // find bbox intersecting tris
    candidates.clear();
    double ist_bbmin[3];
    double ist_bbmax[3];
    FOR_J3 ist_bbmin[j] = min(xsp[spost[ist][0]][j],min(xsp[spost[ist][1]][j],xsp[spost[ist][2]][j])) - dist_tol;
    FOR_J3 ist_bbmax[j] = max(xsp[spost[ist][0]][j],max(xsp[spost[ist][1]][j],xsp[spost[ist][2]][j])) + dist_tol;
    stAdt.buildListForBBox(candidates,ist_bbmin,ist_bbmax);

    for (c_it=candidates.begin(); c_it!=candidates.end(); ++c_it) {
      const int c_ist = *c_it;

      if (c_ist == ist) continue;  // skip if returned myself

      // distance from c_ist nodes to ist plane; augment with shareNode info
      FOR_I3 signedDistToIst[i] = DOT_PRODUCT(N_ist,xsp[spost[c_ist][i]]) + d_ist;
      signToIst.setAll(0);
      FOR_I3 {
        if (signedDistToIst[i] < -dist_tol) signToIst[i] = -1;
        else if (signedDistToIst[i] > dist_tol) signToIst[i] = 1;
      }

      // if candidate nodes are all to one side of this tri, ignore immediately
      if ((signToIst.countPositive() == 3) || (signToIst.countNegative() == 3)) continue;

      // compute dist of ist nodes to c_ist plane; augment with shareNode info
      MiscUtils::computePlaneCoeffsFromPoints(N_c_ist,d_c_ist,xsp[spost[c_ist][0]],xsp[spost[c_ist][1]],xsp[spost[c_ist][2]]);
      FOR_I3 signedDistToCIst[i] = DOT_PRODUCT(N_c_ist,xsp[spost[ist][i]]) + d_c_ist;

      // MiscUtils::getTriUnitNormal(N_c_ist,xsp[spost[c_ist][0]],xsp[spost[c_ist][1]],xsp[spost[c_ist][2]]);
      // for (int i=0; i<3; ++i) {
      //   const double dx[3] = DIFF(xsp[spost[ist][i]],xsp[spost[c_ist][0]]);
      //   signedDistToCIst[i] = DOT_PRODUCT(N_c_ist,dx);
      // }
      signToCIst.setAll(0);
      FOR_I3 {
        if (signedDistToCIst[i] < -dist_tol)  signToCIst[i] = -1;
        else if (signedDistToCIst[i] > dist_tol) signToCIst[i] = 1;
      }
      // if my nodes are all to one side of candidate tri, ignore immediately
      if ((signToCIst.countPositive() == 3) || (signToCIst.countNegative() == 3)) continue;

      // if all nodes within other tri's plane (within tol) then assume co-planar
      bool isCoplanar = ((signToIst.countZero() == 3) || (signToCIst.countZero() == 3)) ? true:false;

      // tris further than tolerance distance from each other have been rejected
      // now look at nodal connectivity of ist and c_ist to organize categorization based on whether
      // these tris share a node or not
      int isEdgeNbr = -1;
      FOR_I3 if ((ist_nbr[i]>=0) && (c_ist == ist_nbr[i])) isEdgeNbr = i;  // based on Teost we are edge nbrs

      // edge neighbors who have already processed their neighbors should be skipped (assumes checks are symmetric)
      if ((isEdgeNbr != -1) && nmTriData.isVisited(c_ist)) continue;

      // determine if formally a node neighbor
      // store which nodes of c_ist are shared with ist and vice_versa
      int shareNodeCount = 0;
      FOR_I3 FOR_J2 shareNode[i][j] = -1;
      FOR_J3 {
        if (spost[ist][0] == spost[c_ist][j]) {
          shareNode[0][0] = j;
          shareNode[j][1] = 0;
          ++shareNodeCount;
        }
        else if (spost[ist][1] == spost[c_ist][j]) {
          shareNode[1][0] = j;
          shareNode[j][1] = 1;
          ++shareNodeCount;
        }
        else if (spost[ist][2] == spost[c_ist][j]) {
          shareNode[2][0] = j;
          shareNode[j][1] = 2;
          ++shareNodeCount;
        }
      }
      // node neighbors who have already processed their neighbors should be skipped
      if (shareNodeCount && nmTriData.isVisited(c_ist)) continue;

      // special case of identically overlapping
      if (shareNodeCount == 3) {
        nmTriData.setMultiNeighbor(ist);
        nmTriData.setMultiNeighbor(c_ist);
        continue;
      }

      // check that data is consistent between teost and spost
      if (isEdgeNbr != -1) assert(shareNodeCount == 2);
      if (shareNodeCount == 1) assert(isEdgeNbr == -1);

      if (isEdgeNbr != -1) {
        assert(isEdgeNbr < 3);
        const bool aligned = isNbrAligned(ist,isEdgeNbr);
        categorizeEdgeNeighbor(ist,c_ist,aligned,angle_tol);
        continue;  // must have been categorized; otherwise skip
      }

      // can have a case where multiply node connected but NOT edge neighbor: flipped adjacent tris
      if (shareNodeCount == 2) {
        // got here is an indication that spost shows these tris are neighbors but teost did not build properly...
        categorizeMultiNodeButNotEdgeNeighbor(ist,c_ist,angle_tol);
        continue;
      }

      // if you share a node based on spost
      if (shareNodeCount == 1) {
        categorizeNodeNeighbor(ist,c_ist,shareNode,signToIst,signToCIst,isCoplanar,nmTriData.distTol());
        continue;  // must have been categorized; otherwise skip
      }

      // -- Vicinity based on Adt --
      // if you reach here then you are not topologically connected

      // check if coplanar based on tri normals and angle_tol
      const double dp = MiscUtils::getCosAngleBetweenTris(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],xsp[spost[c_ist][0]],xsp[spost[c_ist][1]],xsp[spost[c_ist][2]]);
      if (abs(dp) > -angle_tol) isCoplanar = true;

      if (isCoplanar) {
        // could have situations where dist tol and angle tol disagree on whether coplanar...consider
        // how to report difference
        categorizeCoplanarVicinityTri(ist,c_ist,nmTriData.distTol());
        continue;  // should have been categorized
      }

      // determine if intersecting/imprinting
      // note that many have already been rejected by this point
      categorizeVicinityTri(ist,c_ist,signToIst,signToCIst,nmTriData.distTol());

    }  // done looping candidates

    // set this tri as having been checked against its node & edge neighbors
    nmTriData.setVisited(ist);
  }
  if (mpi_rank == 0) cout << "..100%" << endl;

  const double diagnose_time = MPI_Wtime()-adt_time-start_time;
  COUT1("    > completed in: " << diagnose_time <<  "s");

  if (nmTriData.nSeedSamples()) {
    nmTriData.getCounts();  // re-compute counts of each type

    UIMessage_ msg(MESSAGE);
    msg.addLine(stringstream().flush() << "    > detected " << nmTriData.countOverlapping() << " overlapping (coplanar) tris");
    if (nmTriData.countOverlapping()) {
      int zone_counts[zoneVec.size()];
      for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) zone_counts[izn] = 0;
      for (int ist=0; ist<nst; ++ist) {
        if (nmTriData.isOverlapping(ist)) ++zone_counts[znost[ist]];  // increment zone counter
      }
      for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
        if (zone_counts[izn]) msg.addLine(stringstream().flush() << "       > detected " << zone_counts[izn] << " problems in " << zoneVec[izn].getName());
      }
    }
    msg.addLine(stringstream().flush() << "    > detected " << nmTriData.countImprinting() << " imprinting (node or edge adjacent) tris");
    if (nmTriData.countImprinting()) {
      int zone_counts[zoneVec.size()];
      for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) zone_counts[izn] = 0;
      for (int ist=0; ist<nst; ++ist) {
        if (nmTriData.isImprinting(ist)) ++zone_counts[znost[ist]];  // increment zone counter
      }
      for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
        if (zone_counts[izn]) msg.addLine(stringstream().flush() << "       > detected " << zone_counts[izn] << " problems in " << zoneVec[izn].getName());
      }
    }
    msg.addLine(stringstream().flush() << "    > detected " << nmTriData.countIntersecting() << " intersecting (penetrating) tris");
    if (nmTriData.countIntersecting()) {
      int zone_counts[zoneVec.size()];
      for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) zone_counts[izn] = 0;
      for (int ist=0; ist<nst; ++ist) {
        if (nmTriData.isIntersecting(ist)) ++zone_counts[znost[ist]];  // increment zone counter
      }
      for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
        if (zone_counts[izn]) msg.addLine(stringstream().flush() << "       > detected " << zone_counts[izn] << " problems in " << zoneVec[izn].getName());
      }
    }
    msg.dumpLines();
    WebUI::webUIOutput.add_(msg);

    testNonmanifoldSeeds(&stAdt,nmTriData.nSeedSamples());
  }
}

bool SimpleSurface::categorizeTwoNodesInPlaneVicinityTri(const int ist,const int c_ist,IntFlag& signToIst,const double dist_tol) {
  const double d2_tol = dist_tol*dist_tol;
  // only possibility is to imprint or not
  // there is a c_ist edge on the plane of ist;
  // see whether either node is inside ist
  FOR_I3 {
    if (signToIst[i] == 0) {
      const int isp = spost[c_ist][i];
      const bool inTri = pointInTriangle(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],xsp[isp]);
      if (inTri) return true;  // only one needs to be in the tri for this to imprint
    }
  }

  // if got here we did not find imprint
  // see if the edge intersect any of the ist edges
  int ispa = -1;
  int ispb = -1;
  FOR_I3 {
    if (ispa == -1 && signToIst[i] == 0) ispa = spost[c_ist][i];
    else if (signToIst[i] == 0) ispb = spost[c_ist][i];
  }
  assert(ispa != -1 && ispb != -1 && ispa != ispb);

  FOR_I3 {
    const double d2 = MiscUtils::getEdgeToEdgeDist2(xsp[ispa],xsp[ispb],xsp[spost[ist][i]],xsp[spost[ist][(i+1)%3]]);

    if (d2 < d2_tol) {
      return true;
    }
  }

  return false;  // did not imprint
}

int SimpleSurface::categorizeOneNodeInPlaneVicinityTri(const int ist,const int c_ist,IntFlag& signToIst,const double dist_tol) {
  // return value by bit:
  // 1: imprint
  // 2: intersect
  int ret_val = 0;

  // node may sit in ist
  bool b_inTri = false;
  FOR_I3 {
    if (signToIst[i] == 0) {
      const int isp = spost[c_ist][i];
      b_inTri = pointInTriangle(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],xsp[isp]);
      if (b_inTri) break;
    }
  }

  // cout << " > b_inTri: " << b_inTri << endl;
  const int npos = signToIst.countPositive();
  const int nneg = signToIst.countNegative();

  if (b_inTri) {
    if (npos == 2 || nneg == 2) ret_val |= (1<<0);
    else ret_val |= (1<<1);
  }
  else {
    // in-plane point is not on ist; check for intersect
    if (npos == 2 || nneg == 2) return ret_val;  // cannot intersect

    assert(npos == 1 && nneg == 1);
    // check whether edge between pos/neg c_ist nodes pierces ist
    int ispa,ispb;
    FOR_I3 {
      if (signToIst[i] == 1) ispa = spost[c_ist][i];
      else if (signToIst[i] == -1) ispb = spost[c_ist][i];
    }
    int intersection = MiscUtils::howDoesLinePierceTri(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],xsp[ispa],xsp[ispb],dist_tol);

    if (intersection == 1) {
      // cout << " > intersection = 1" << endl;
      ret_val |= (1<<1);
    }
    else if (intersection != 0) {
      ret_val |= (1<<0);
    }
    else if (intersection == 0) {
      // ist may pierce c_ist; should be categorized when c_ist processed. But to be safe do it here too
      //NOTE if can get rid of asymmetric categorizations then can ignore previously flagged tris...

      const int pierce = checkForValidTriPiercing(ist,c_ist);
      if ((pierce == 2) || (pierce == 3)) ret_val |= (1<<1);
      if ((pierce == 1) || (pierce == 3)) ret_val |= (1<<0);

      // check if any edge of ist pierces c_ist
      // FOR_I3 {
      //   ispa = spost[ist][i];
      //   ispb = spost[ist][(i+1)%3];
      //   intersection = MiscUtils::howDoesLinePierceTri(xsp[spost[c_ist][0]],xsp[spost[c_ist][1]],xsp[spost[c_ist][2]],xsp[ispa],xsp[ispb],dist_tol);
      //
      //   // only check for interior intersection; edge/node matches would have registered when processing c_ist:ist
      //   if (intersection == 1) {
      //     cout << " > intersection = 0; ist pierce c_ist" << endl;
      //     ret_val |= (1<<1);
      //     break;
      //   }
      // }
    }
  }
  return ret_val;
}

void SimpleSurface::categorizeEdgeNeighbor(const int ist,const int c_ist,const bool aligned,const double angle_tol) {

  // only two nodes shared
  // check angle at this edge and whether it passes tolerance
  if (nmTriData.isLinear(ist) || nmTriData.isLinear(c_ist) || nmTriData.isCollapsed(ist) || nmTriData.isCollapsed(c_ist)) {
    return;  // normal is undefined for a linear tri, so skip
  }

  // check of overlapping co-planar tris:
  // use dot product of unit normals - if within tolerance then consider overlapping
  const double dp = MiscUtils::getCosAngleBetweenTris(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],xsp[spost[c_ist][0]],xsp[spost[c_ist][1]],xsp[spost[c_ist][2]]);

  if ( (aligned && (dp < angle_tol)) || ( (!aligned) && (dp > -angle_tol)) ) {
    // this is an overlapping pair of neighboring tris
    nmTriData.setOverlapping(ist);
    nmTriData.setOverlapping(c_ist);
  }
}

void SimpleSurface::categorizeMultiNodeButNotEdgeNeighbor(const int ist,const int c_ist,const double angle_tol) {

  // check angle at this edge and whether it passes tolerance
  if (nmTriData.isLinear(ist) || nmTriData.isLinear(c_ist) || nmTriData.isCollapsed(ist) || nmTriData.isCollapsed(c_ist)) {
    return;  // normal is undefined for a linear tri, so skip
  }

  // check of overlapping co-planar tris:
  // use dot product of unit normals - if within tolerance then consider overlapping
  const double dp = MiscUtils::getCosAngleBetweenTris(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],xsp[spost[c_ist][0]],xsp[spost[c_ist][1]],xsp[spost[c_ist][2]]);

  if (dp > -angle_tol) {
    // use negative tol b/c normals are in the same direction if flipped
    nmTriData.setOverlapping(ist);
    nmTriData.setOverlapping(c_ist);
  }
}

void SimpleSurface::categorizeNodeNeighbor(const int ist,const int c_ist,const int (&shareNode)[3][2],const IntFlag& signToIst,const IntFlag& signToCIst,const bool isCoplanar,const double dist_tol) {
  // when called should already have checked that only a single node is shared
  // determine shared node index for each tri

  int sn_ist = -1;
  int sn_c_ist = -1;
  FOR_I3 {
    if (shareNode[i][0] != -1) {
      sn_ist = i;
      break;
    }
  }
  FOR_I3 {
    if (shareNode[i][1] != -1) {
      sn_c_ist = i;
      break;
    }
  }
  assert(sn_ist != -1);
  assert(sn_c_ist != -1);

  // b/c of shared node at least one dist should be zero
  const int nzero = signToIst.countZero();
  const int npos = signToIst.countPositive();
  const int nneg = signToIst.countNegative();
  if (nzero == 0) {
    WUI(WARN,"Shared node being seeing as non-coplanar. RUN_DIAGNOSTICS DIST_TOL (default 1E-10) not set appropiately and/or ist's: " << ist << " " << c_ist << " are of poor quality.");
  }
  assert(npos+nneg+nzero == 3);

  const int c_nzero = signToCIst.countZero();
  const int c_npos = signToCIst.countPositive();
  const int c_nneg = signToCIst.countNegative();
  if (c_nzero == 0) {
    WUI(WARN,"Shared node being seeing as non-coplanar. RUN_DIAGNOSTICS DIST_TOL (default 1E-10) not set appropiately and/or ist's: " << ist << " " << c_ist << " are of poor quality.");
  }
  assert(c_npos+c_nneg+c_nzero == 3);

  if (isCoplanar || nzero == 3 || c_nzero == 3) {
    bool b_overlap = false;

    // any edge-edge intersection indicates overlap
    // b/c of shared node, the edge between the free nodes of one of the tris is
    // guaranteed to intersect an edge on the other
    // use sign of tet_vol to determine edge intersection

    // check free nodes of c_ist against all edges of ist
    int ispa = -1;
    int ispb = -1;
    FOR_I3 {
      if (i != sn_c_ist) {
        if (ispa == -1) ispa = spost[c_ist][i];
        else ispb = spost[c_ist][i];
      }
    }
    assert(ispa != -1 && ispb != -1);
    FOR_I3 {
      const double d2 = MiscUtils::getEdgeToEdgeDist2(xsp[ispa],xsp[ispb],xsp[spost[ist][i]],xsp[spost[ist][(i+1)%3]]);
      if (d2 < 0.0) continue;  // parallel edge; should only be possible for other free edge...

      if (d2 < dist_tol*dist_tol) {
        // edges are coplanar and intersect within tolerance
        // COUT2("ZONE T=\"edges\" N=4, E=2, F=FEPOINT, ET=TRIANGLE");
        // COUT2(xsp[ispa][0] << " " << xsp[ispa][1] << " " << xsp[ispa][2]);
        // COUT2(xsp[ispb][0] << " " << xsp[ispb][1] << " " << xsp[ispb][2]);
        // COUT2(xsp[spost[ist][i]][0] << " " << xsp[spost[ist][i]][1] << " " << xsp[spost[ist][i]][2]);
        // COUT2(xsp[spost[ist][(i+1)%3]][0] << " " << xsp[spost[ist][(i+1)%3]][1] << " " << xsp[spost[ist][(i+1)%3]][2]);
        // COUT2("1 2 2\n3 4 4");
        b_overlap = true;
        break;
      }
    }

    if (!b_overlap) {
      // check free nodes of ist against all edges of c_ist
      ispa = -1;
      ispb = -1;
      FOR_I3 {
        if (i != sn_ist) {
          if (ispa == -1) ispa = spost[ist][i];
          else ispb = spost[ist][i];
        }
      }
      assert(ispa != -1 && ispb != -1);
      FOR_I3 {
        const double d2 = MiscUtils::getEdgeToEdgeDist2(xsp[ispa],xsp[ispb],xsp[spost[c_ist][i]],xsp[spost[c_ist][(i+1)%3]]);
        if (d2 < 0.0) continue;  // doesn't intersect...

        if (d2 < dist_tol*dist_tol) {
          // edges are coplanar and intersect within tolerance
          b_overlap = true;
          break;
        }
      }
    }

    if (!b_overlap) {
      // check if c_ist free nodes are inside ist
      FOR_I3 {
        if (i == sn_c_ist) continue;  // skip shared node

        // check non-shared nodes and test if in tri
        const bool inTri = pointInTriangle(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],xsp[spost[c_ist][i]]);

        if (inTri) {
          b_overlap = true;
          break;  // only one needs to be in tri to overlap
        }
      }
    }

    if (!b_overlap) {
      // check if ist free nodes are inside c_ist
      FOR_I3 {
        if (i != sn_ist) {
          const bool inTri = pointInTriangle(xsp[spost[c_ist][0]],xsp[spost[c_ist][1]],xsp[spost[c_ist][2]],xsp[spost[ist][i]]);

          if (inTri) {
            b_overlap = true;
            break;  // only one needs to be in tri to overlap
          }
        }
      }
    }

    if (b_overlap) {
      // this is an overlapping pair of neighboring tris
      nmTriData.setOverlapping(ist);
      nmTriData.setOverlapping(c_ist);
    }
  }
  else {
    // shares a single node with ist and NOT coplanar

    // check for intersection of other nodes or imprinting
    if (npos == 2 || nneg == 2 || c_npos == 2 || c_nneg == 2) return;  // no intersection

    bool b_intersect = false;
    bool b_imprint = false;

    if (npos == 1 && nneg == 1) {
      // one node of c_ist on either side of ist plane
      // check if free c_ist edge intersects ist or not
      int ispa = -1;
      int ispb = -1;
      FOR_I3 {
        if (i != sn_c_ist) {
          if (ispa == -1) ispa = spost[c_ist][i];
          else ispb = spost[c_ist][i];
        }
      }
      assert(ispa != -1 && ispb != -1 && ispa != ispb);

      int intersection = MiscUtils::howDoesLinePierceTri(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],xsp[ispa],xsp[ispb],dist_tol);

      if (intersection == 0) {
        // check if free ist edge intersects c_ist or not
        ispa = -1;
        ispb = -1;
        FOR_I3 {
          if (i != sn_ist) {
            if (ispa == -1) ispa = spost[ist][i];
            else ispb = spost[ist][i];
          }
        }
        assert(ispa != -1 && ispb != -1 && ispa != ispb);

        intersection = MiscUtils::howDoesLinePierceTri(xsp[spost[c_ist][0]],xsp[spost[c_ist][1]],xsp[spost[c_ist][2]],xsp[ispa],xsp[ispb],dist_tol);
      }

      // catergorize based on response
      if (intersection == 1) b_intersect = true;
      else if (intersection != 0) b_imprint = true;  // node or edge hits
    }
    else {
      assert(nzero == 2);
      // two (one unshared) zero distance nodes
      // check if other zero dist point is in ist
      FOR_I3 {
        if (signToIst[i] == 0 && i != sn_c_ist) {
          const int isp = spost[c_ist][i];
          // const bool inTri = pointInTriangle(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],xsp[isp]);
          const double d2 = MiscUtils::getPointToTriDist2(xsp[isp],xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);

          if (d2 < dist_tol*dist_tol) {
            b_imprint = true;
            break;  // only one can be on tri since not coplanar
          }
        }
      }

      // check if edge from shared node to this point intersects ist free edge
      if (!b_imprint) {
        // get free edge nodes of ist
        int ispa = -1;
        int ispb = -1;
        FOR_I3 {
          if (i != sn_ist) {
            if (ispa == -1) ispa = spost[ist][i];
            else ispb = spost[ist][i];
          }
        }
        assert(ispa != -1 && ispb != -1 && ispa != ispb);

        int count = 0;
        FOR_I3 {
          if (signToIst[i] == 0 && i != sn_c_ist) {
            ++count;
            const int sn_isp = spost[c_ist][sn_c_ist];  // shared node isp
            const int isp = spost[c_ist][i];

            const double d2 = MiscUtils::getEdgeToEdgeDist2(xsp[ispa],xsp[ispb],xsp[isp],xsp[sn_isp]);
            if (d2 < 0.0) continue;  // doesn't imprint...

            if (d2 < dist_tol*dist_tol) {
              // edges are close enough that considered an edge-intersection
              b_imprint = true;
              break;
            }
          }
        }
        if (count != 1) {
          WUI(WARN,"Shared node being seeing as not imprinting. RUN_DIAGNOSTICS DIST_TOL (default 1E-10) not set appropiately and/or ist's     : " << ist << " " << c_ist << " are of poor quality.");
        }
      }
    }

    if (b_intersect) {
      nmTriData.setIntersecting(ist);
      nmTriData.setIntersecting(c_ist);
    }
    else if (b_imprint) {
      nmTriData.setImprinting(ist);
      nmTriData.setImprinting(c_ist);
    }
  }
}

void SimpleSurface::categorizeCoplanarVicinityTri(const int ist,const int c_ist,const double dist_tol) {
  bool b_overlap = false;
  bool b_imprint = false;

  FOR_I3 {
    // check if any c_ist nodes inside ist
    const int isp = spost[c_ist][i];

    const int onTri = howIsPointOnTriangle(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],xsp[isp],dist_tol);

    if (onTri == 1) {
      b_overlap = true;
      break;  // only one needs to be in tri to overlap
    }
    else if (onTri < 0 || onTri > 1) {
      // node on node or edge of c_ist
      b_imprint = true;
    }
  }

  if (!b_overlap) {
    FOR_I3 {
      // check if any ist nodes inside c_ist
      const int isp = spost[ist][i];
      const int onTri = howIsPointOnTriangle(xsp[spost[c_ist][0]],xsp[spost[c_ist][1]],xsp[spost[c_ist][2]],xsp[isp],dist_tol);

      if (onTri == 1) {
        b_overlap = true;
        break;  // only one needs to be in tri to overlap
      }
      else if (onTri < 0 || onTri > 1) {
        b_imprint = true;
      }
    }
  }

  if (!b_overlap) {
    // check for edge intersections (area overlap)
    FOR_I3 {
      const int ispa = spost[ist][i];
      const int ispb = spost[ist][(i+1)%3];
      FOR_J3 {
        const int c_ispa = spost[c_ist][j];
        const int c_ispb = spost[c_ist][(j+1)%3];
        const double d2 = MiscUtils::getEdgeToEdgeDist2(xsp[ispa],xsp[ispb],xsp[c_ispa],xsp[c_ispb]);

        if (d2 >= 0.0 && d2 < dist_tol*dist_tol) {
          b_overlap = true;
          break;  // only one needs to be in tri to overlap
        }
      }
      if (b_overlap) break;
    }
  }

  if (b_overlap) {
    // this is an overlapping pair of neighboring tris;
    nmTriData.setOverlapping(ist);
    nmTriData.setOverlapping(c_ist);
    return;
  }
  if (b_imprint) {
    nmTriData.setImprinting(ist);
    nmTriData.setImprinting(c_ist);
  }
}

void SimpleSurface::categorizeVicinityTri(const int ist,const int c_ist,IntFlag& signToIst,IntFlag& signToCIst,const double dist_tol) {
  UNUSED(dist_tol);
  bool b_imprint = false;
  bool b_intersect = false;

  // below commented code is now handled earlier
  // // checks that reach here have no shared nodes & are not coplanar
  const int nzero = signToIst.countZero();
  // const int npos = signToIst.countPositive();
  // const int nneg = signToIst.countNegative();
  // assert(npos+nneg+nzero == 3);
  //
  // const int c_npos = signToCIst.countPositive();
  // const int c_nneg = signToCIst.countNegative();
  const int c_nzero = signToCIst.countZero();
  // assert(c_npos+c_nneg+c_nzero == 3);
  //
  // // check if nodes of each tri are to one side of the others' plane
  // if (npos==3 || nneg==3 || c_npos==3 || c_nneg==3) return;
  //
  // // cout << "----------------------" << endl;
  //
  // if (nzero == 3) {
  //   // cout << "nzero == 3" << endl;
  //   categorizeCoplanarVicinityTri(ist,c_ist,nmTriData.distTol());
  //   return;
  // }
  // if (c_nzero == 3) {
  //   // cout << "c_nzero == 3" << endl;
  //   categorizeCoplanarVicinityTri(c_ist,ist,nmTriData.distTol());
  //   return;
  // }

  if (nzero == 2) {
    // cout << "nzero == 2" << endl;
    b_imprint = categorizeTwoNodesInPlaneVicinityTri(ist,c_ist,signToIst,nmTriData.distTol());
  }
  else if (nzero == 1) {
    // cout << "nzero == 1" << endl;
    const int type = categorizeOneNodeInPlaneVicinityTri(ist,c_ist,signToIst,nmTriData.distTol());
    b_imprint = (type & (1<<0));  // if imprint bit set
    b_intersect = (type & (1<<1));  // if intersect bit set
    // if (b_intersect) cout << " > intersect" << endl;
  }
  else {
    assert(nzero == 0);
    // cout << "nzero == 0" << endl;
    if (c_nzero == 1) {
      // cout << "c_nzero == 1" << endl;
      const int type = categorizeOneNodeInPlaneVicinityTri(c_ist,ist,signToCIst,nmTriData.distTol());
      b_imprint = (type & (1<<0));  // if imprint bit set
      b_intersect = (type & (1<<1));  // if intersect bit set

    }
    else if (c_nzero == 2) {
      // cout << "c_nzero == 2" << endl;
      b_imprint = categorizeTwoNodesInPlaneVicinityTri(c_ist,ist,signToCIst,nmTriData.distTol());
    }
    else {
      assert(c_nzero == 0);

      // cout << "c_nzero == 0" << endl;
      const int ret_val = checkForValidTriPiercing(ist,c_ist);
      if ((ret_val == 2) || (ret_val == 3)) b_intersect = true;
      if ((ret_val == 1) || (ret_val == 3)) b_imprint = true;
    }
  }

  if (b_intersect) {
    // takes precedence over cases where both intersect & imprint found
    nmTriData.setIntersecting(ist);
    nmTriData.setIntersecting(c_ist);
  }
  else if (b_imprint) {
    nmTriData.setImprinting(ist);
    nmTriData.setImprinting(c_ist);
  }
}

int SimpleSurface::checkForValidTriPiercing(const int ist,const int c_ist) const {

  bool b_imprint = false;
  bool b_intersect = false;
  // check if any ist edges pierce c_ist
  int c_ist_intersection = 0;
  int ist_intersection = 0;
  int counts_ist[2] = {0,0}; // pierces, boundary intersections
  int counts_c_ist[2] = {0,0}; // pierces, boundary intersections
  FOR_I3 {
    const int ispa = spost[ist][i];
    const int ispb = spost[ist][(i+1)%3];
    c_ist_intersection = MiscUtils::howDoesLinePierceTri(xsp[spost[c_ist][0]],xsp[spost[c_ist][1]],xsp[spost[c_ist][2]],xsp[ispa],xsp[ispb],nmTriData.distTol());

    if (c_ist_intersection == 1) {
      ++counts_c_ist[0];
    }
    else if (c_ist_intersection != 0) {
      ++counts_c_ist[1];  // boundary (edge/node) intersection
      b_imprint = true;
    }
  }


  // check if any c_ist edges pierce ist
  FOR_I3 {
    const int ispa = spost[c_ist][i];
    const int ispb = spost[c_ist][(i+1)%3];
    ist_intersection = MiscUtils::howDoesLinePierceTri(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],xsp[ispa],xsp[ispb],nmTriData.distTol());
    if (ist_intersection == 1) {
      ++counts_ist[0];
    }
    else if (ist_intersection != 0) {
      ++counts_ist[1];  // boundary (edge/node) intersection
      b_imprint = true;
    }
  }

  // piercings must be understandably symmetric or assymetric, this protects against issues with tet-vol based computation in nearly parallel tris
  if (
      ((counts_ist[0] == counts_c_ist[0]) && (counts_ist[0] == 1)) ||             // symmetric piercing
      ((counts_ist[0] == 2) && (counts_c_ist[0] == 0)) ||                         // one tri pierces interior of another
      ((counts_ist[0] == 0) && (counts_c_ist[0] == 2)) ||
      ((counts_ist[0] == 1) && (counts_ist[1] == 1) && (counts_c_ist[1] == 1)) || // edge and interior piercing
      ((counts_c_ist[0] == 1) && (counts_c_ist[1] == 1) && (counts_ist[1] == 1))
      ) {
    // cout << "nzero == pierce" << endl;
    b_intersect = true;
  }

  // set return value
  int ret_val = 0;
  if (b_intersect) ret_val |= (1<<1);
  if (b_imprint) ret_val |= (1<<0);
  return ret_val;
}

void SimpleSurface::mergeNodesOfOpenEdges(const double merge_tol) {

  COUT2("SimpleSurface::mergeNodesOfOpenEdges");
  COUT2(" > merging nodes within tol: "<< merge_tol);

  // needed to determine which edges are open
  ensureTeost();

  sp_flag.setLength(nsp);
  sp_flag.setAll(-1);

  // flag nodes on open edges
  FOR_IST {
    FOR_I3 {
      if (isEdgeOpen(ist,i)) {
        sp_flag[spost[ist][i]] = 1;
        sp_flag[spost[ist][(i+1)%3]] = 1;
      }
    }
  }
  const int n_oe_nodes = sp_flag.countPositive();

  // now index and build bbox
  double (*bbmin)[3] = new double[n_oe_nodes][3];
  double (*bbmax)[3] = new double[n_oe_nodes][3];

  int isp_index = 0;
  int *oen_to_isp = new int[n_oe_nodes];
  FOR_ISP {
    if (sp_flag[isp] != -1) {
      oen_to_isp[isp_index] = isp;
      FOR_I3 bbmin[isp_index][i] = xsp[isp][i] - merge_tol;
      FOR_I3 bbmax[isp_index][i] = xsp[isp][i] + merge_tol;
      ++isp_index;
    }
    sp_flag[isp] = isp;  // reset sp_flag to own index - used next for compress
  }
  assert(n_oe_nodes==isp_index);
  Adt<double> nodesAdt = Adt<double>(n_oe_nodes,bbmin,bbmax);
  delete[] bbmin;
  delete[] bbmax;

  // match candidates
  vector<int> candidateVec;
  for (int ino = 0; ino < n_oe_nodes; ++ino) {
    candidateVec.clear();
    const int isp = oen_to_isp[ino];
    nodesAdt.buildListForPoint(candidateVec,xsp[isp]);

    if (candidateVec.empty()) continue;  // no matches found; skip to next node

    for (vector<int>::iterator it=candidateVec.begin(); it!=candidateVec.end(); ++it) {
      const int c_oen_index = *it;
      if (c_oen_index == ino) continue;  // don't process myself

      const int c_isp = oen_to_isp[c_oen_index];
      if (c_isp < isp) sp_flag[isp] = c_isp;
      if (c_isp > isp) sp_flag[c_isp] = isp;
    }
  }
  delete[] oen_to_isp;

  // sp_flag should now be indexed such that nodes who are being compressed index their
  // lower partner

  // push sp_flag to unique value
  FOR_ISP {
    if (sp_flag[isp] != isp) {
      int isp_t = sp_flag[isp];
      while (sp_flag[isp_t] != isp_t) isp_t = sp_flag[isp_t];
      sp_flag[isp] = isp_t;
    }
  }

  // now compress sp_flag & xsp
  int nsp_new = 0;
  double * nsp_merged = new double[nsp]; // nsp_new <= nsp
  FOR_ISP {
    nsp_merged[isp] = 0;
    if (sp_flag[isp] == isp) {
      if (nsp_new != isp) {
        FOR_I3 xsp[nsp_new][i] = xsp[isp][i];
      }
      assert(nsp_merged[nsp_new] == 0);
      ++nsp_merged[nsp_new];
      sp_flag[isp] = nsp_new;
      ++nsp_new;
    }
    else {
      sp_flag[isp] = sp_flag[sp_flag[isp]];  // point to appropriate node
      assert(nsp_merged[sp_flag[isp]] >= 1);
      FOR_I3 xsp[sp_flag[isp]][i] += xsp[isp][i];
      ++nsp_merged[sp_flag[isp]];
    }
  }
  for (int isp = 0; isp < nsp_new; ++isp) {
    assert(nsp_merged[isp] >= 1);
    if (nsp_merged[isp] > 1) {
      double inv_n = 1.0/(double)nsp_merged[isp];
      FOR_I3 xsp[isp][i] *= inv_n;
    }
  }
  delete[] nsp_merged;

  // update spost
  FOR_IST {
    FOR_I3 {
      int isp = spost[ist][i];
      spost[ist][i] = sp_flag[isp];
    }
  }

  // if we compressed any nodes, teost is outdated
  if (nsp != nsp_new) {
    clearTeost();
    clearNonManifoldData();
  }

  // resize nsp
  nsp = nsp_new;

  // znost, szost, and szozn_i don't change

}

void SimpleSurface::mergeNodesOfOpenEdgeGroups(const int group0, const int group1, const double merge_tol) {

  COUT2("SimpleSurface::mergeNodesOfOpenEdgeGroups");
  COUT2(" > merging nodes from groups: " << group0 << "," << group1);
  ensureOpenEdgeGroups();

  COUT2(" > building local edge information");
  vector<oegEdge> edges0;
  vector<oegEdge> edges1;
  for (int ist=0; ist < nst; ++ist) {
    FOR_I3 {
      int igr;
      if (getOpenEdgeGroup(igr,ist,i)) {
        if (igr == group0 || igr == group1) {
          const int isp0 = spost[ist][i];
          const int isp1 = spost[ist][(i+1)%3];
          const double d2 = DIST2(xsp[isp0],xsp[isp1]);

          if (igr == group0) edges0.push_back(oegEdge(igr,isp0,isp1,d2));
          else if (igr == group1) edges1.push_back(oegEdge(igr,isp0,isp1,d2));
        }
      }
    }
  }

  mergeNodesOfOpenEdgeVectors(edges0,edges1,merge_tol);

}

void SimpleSurface::mergeNodesOfOpenEdgeGroupsOnSubzones(const int zone0, const int zone1, const double merge_tol) {

  COUT2("SimpleSurface::mergeNodesOfOpenEdgeGroupsOnZones");
  COUT2(" > merging nodes from subzones: " << zone0 << "," << zone1);
  ensureOpenEdgeGroups();

  COUT2(" > building local edge information");
  vector<oegEdge> edges0;
  vector<oegEdge> edges1;
  for (int ist=0; ist < nst; ++ist) {
    if (szost[ist] == zone0 || szost[ist] == zone1) {
      FOR_I3 {
        int igr;
        if (getOpenEdgeGroup(igr,ist,i)) {
          const int isp0 = spost[ist][i];
          const int isp1 = spost[ist][(i+1)%3];
          const double d2 = DIST2(xsp[isp0],xsp[isp1]);

          if (szost[ist] == zone0) edges0.push_back(oegEdge(igr,isp0,isp1,d2));
          else if (szost[ist] == zone1) edges1.push_back(oegEdge(igr,isp0,isp1,d2));
        }
      }
    }
  }

  mergeNodesOfOpenEdgeVectors(edges0,edges1,merge_tol);

}

void SimpleSurface::mergeNodesOfOpenEdgeVectors(vector<oegEdge>& edges0, vector<oegEdge>& edges1, const double merge_tol) {
  COUT2(" > determining edge lengthscale");
  // compute min edge distance for each node in group
  IntFlag _sp_flag(nsp);
  _sp_flag.setAll(-1);
  vector<oegNode> nodes0;
  for (int ied=0, limit=edges0.size(); ied<limit; ++ied) {
    // update start node
    if (_sp_flag[edges0[ied].isp0] < 0) {
      // node not processed; add to vector
      _sp_flag[edges0[ied].isp0] = nodes0.size();
      nodes0.push_back(oegNode(edges0[ied].isp0,edges0[ied].d2));
    }
    else {
      // update min edge distance
      const int index = _sp_flag[edges0[ied].isp0];
      nodes0[index].setMinD2(edges0[ied].d2);  // this automatically takes min
    }

    // update end node
    if (_sp_flag[edges0[ied].isp1] < 0) {
      // node not processed; add to vector
      _sp_flag[edges0[ied].isp1] = nodes0.size();
      nodes0.push_back(oegNode(edges0[ied].isp1,edges0[ied].d2));
    }
    else {
      // update min edge distance
      const int index = _sp_flag[edges0[ied].isp1];
      nodes0[index].setMinD2(edges0[ied].d2);  // this automatically takes min
    }
  }


  // compute min edge distance for each node in group
  _sp_flag.setAll(-1);
  vector<oegNode> nodes1;
  for (int ied=0, limit=edges1.size(); ied<limit; ++ied) {
    // update start node
    if (_sp_flag[edges1[ied].isp0] < 0) {
      // node not processed; add to vector
      _sp_flag[edges1[ied].isp0] = nodes1.size();
      nodes1.push_back(oegNode(edges1[ied].isp0,edges1[ied].d2));
    }
    else {
      // update min edge distance
      const int index = _sp_flag[edges1[ied].isp0];
      nodes1[index].setMinD2(edges1[ied].d2);  // this automatically takes min
    }

    // update end node
    if (_sp_flag[edges1[ied].isp1] < 0) {
      // node not processed; add to vector
      _sp_flag[edges1[ied].isp1] = nodes1.size();
      nodes1.push_back(oegNode(edges1[ied].isp1,edges1[ied].d2));
    }
    else {
      // update min edge distance
      const int index = _sp_flag[edges1[ied].isp1];
      nodes1[index].setMinD2(edges1[ied].d2);  // this automatically takes min
    }
  }

  const int nno0 = nodes0.size();
  const int nno1 = nodes1.size();
  COUT2(" > nodes on each group: " << nno0 << "," << nno1);
  if (nno0 == 0 || nno1 == 0) {
    COUT2(" > no nodes available for matching on at least one group; skipping");
    return;
  }

  COUT2(" > building node Adts:");
  // now build Adts for each set of nodes
  double (*bbmin)[3] = new double[nno0][3];
  double (*bbmax)[3] = new double[nno0][3];
  for (int ino0=0; ino0 < nno0; ++ino0) {
    const int isp = nodes0[ino0].isp;
    FOR_I3 bbmin[ino0][i] = xsp[isp][i] - 0.5*sqrt(nodes0[ino0].minD2);
    FOR_I3 bbmax[ino0][i] = xsp[isp][i] + 0.5*sqrt(nodes0[ino0].minD2);
  }
  Adt<double> nodes0Adt = Adt<double>(nno0,bbmin,bbmax);


  delete[] bbmin;
  delete[] bbmax;
  bbmin = new double[nno1][3];
  bbmax = new double[nno1][3];
  for (int ino1=0; ino1 < nno1; ++ino1) {
    const int isp = nodes1[ino1].isp;
    FOR_I3 bbmin[ino1][i] = xsp[isp][i] - 0.5*sqrt(nodes1[ino1].minD2);
    FOR_I3 bbmax[ino1][i] = xsp[isp][i] + 0.5*sqrt(nodes1[ino1].minD2);
  }
  Adt<double> nodes1Adt = Adt<double>(nno1,bbmin,bbmax);
  DELETE(bbmin);
  DELETE(bbmax);

  // set node merge tolerances
  //NOTE must be set before can apply matches
  if (merge_tol != -1) {
    COUT2(" > merge tolerance: " << merge_tol);
    for (int ino = 0; ino < nno0; ++ino) nodes0[ino].setTol(merge_tol);
    for (int ino = 0; ino < nno1; ++ino) nodes1[ino].setTol(merge_tol);
  }
  else {
    // default is to use 1% of smallest adjacent edge length
    COUT2(" > merge tolerance: " << 0.05 << "% of min(adjacent edge length)");
    for (int ino = 0; ino < nno0; ++ino) nodes0[ino].setTol(0.05*0.05*nodes0[ino].minD2);
    for (int ino = 0; ino < nno1; ++ino) nodes1[ino].setTol(0.05*0.05*nodes1[ino].minD2);
  }


  // reuse minD2 to hold matching node distance
  for (int ino = 0; ino < nno0; ++ino) nodes0[ino].clearMinD2();
  for (int ino = 0; ino < nno1; ++ino) nodes1[ino].clearMinD2();


  COUT2(" > searching for node-matches");
  // loop and find candidates for each node
  vector<int> candidateVec;

  for (int ino = 0; ino < nno1; ++ino) {
    candidateVec.clear();
    const int isp1 = nodes1[ino].isp;
    nodes0Adt.buildListForPoint(candidateVec,xsp[isp1]);

    if (candidateVec.empty()) continue;  // no matches found; skip to next node

    for (int cno=0, limit=candidateVec.size(); cno < limit; ++cno) {
      const int ino0 = candidateVec[cno];
      const int isp0 = nodes0[ino0].isp;
      const double d2 = DIST2(xsp[isp0],xsp[isp1]);

      // store the nodes1-index of the match, NOT the actual isp
      nodes0[ino0].setMatch(ino,d2);  // in setMatch this takes the closest node only
    }
  }

  for (int ino = 0; ino < nno0; ++ino) {
    candidateVec.clear();
    const int isp0 = nodes0[ino].isp;
    nodes1Adt.buildListForPoint(candidateVec,xsp[isp0]);

    if (candidateVec.empty()) continue;  // no matches found; skip to next node

    for (int cno=0, limit=candidateVec.size(); cno < limit; ++cno) {
      const int ino1 = candidateVec[cno];
      const int isp1 = nodes1[ino1].isp;
      const double d2 = DIST2(xsp[isp1],xsp[isp0]);

      // store the nodes0-index of the match, NOT the actual isp
      nodes1[ino1].setMatch(ino,d2);  // in setMatch this takes the closest node only
    }
  }

  COUT2(" > checking symmetry of matches");
  // check for match symmetry; if not symmetric then remove match on both
  for (int ino = 0; ino < nno0; ++ino) {
    const int ino1 = nodes0[ino].match;
    if (ino1 == -1) {
      nodes0[ino].match = -1;
      continue;
    }

    const int their_match = nodes1[ino1].match;
    if (their_match != ino) {
      nodes1[ino1].match = -1;
      nodes0[ino].match = -1;
    }
  }

  for (int isp=0; isp < nsp; ++isp) {
    _sp_flag[isp] = isp;
  }

  int symm_count = 0;
  for (int ino = 0; ino < nno1; ++ino) {
    const int ino0 = nodes1[ino].match;
    if (ino0 == -1) continue;  // not symmetric, so skip processing

    const int their_match = nodes0[ino0].match;
    if (their_match != ino) {
      nodes0[ino0].match = -1;
      nodes1[ino].match = -1;
    }
    else {
      // check match is within tolerance
      assert(nodes0[ino0].minD2 == nodes1[ino].minD2);
      if (nodes0[ino0].minD2 <= nodes0[ino0].tol && nodes1[ino].minD2 <= nodes1[ino].tol) {
        // we match after checking from both sides
        ++symm_count;
        const int isp_min = min(nodes0[ino0].isp,nodes1[ino].isp);
        const int isp_max = max(nodes0[ino0].isp,nodes1[ino].isp);
        _sp_flag[isp_max] = -isp_min-1;  // negative 1-indexed of actual matching node

        // modify xsp[isp_min] to be average of two nodes
        // should not clobber elsewhere b/c nodes can only have a single match
        FOR_I3 xsp[isp_min][i] = 0.5*(xsp[isp_min][i] + xsp[isp_max][i]);
      }
    }
  }

  if (symm_count == 0) {
    COUT2(" > no symmetric matches found");
    return;
  }
  COUT2(" > found " << symm_count << " symmetric node matches");

  // simply reduce nsp and over-write values
  const int nsp_new = nsp - symm_count;
  double (*xsp0)[3] = xsp;
  xsp = new double[nsp_new][3];
  int _nsp = 0;
  for (int isp=0; isp < nsp; ++isp) {
    if (_sp_flag[isp] >= 0) {
      _sp_flag[isp] = _nsp;  // for values > 0 this is the map to new indexing
      FOR_I3 xsp[_nsp][i] = xsp0[isp][i];
      ++_nsp;
    }
  }
  assert(_nsp == nsp_new);

  if (nsp != _nsp) {
    clearTeost();
    clearNonManifoldData();
  }
  nsp = _nsp;

  // spost has same number of tris
  // place current node-index in place
  for (int ist=0; ist < nst; ++ist) {
    FOR_I3 {
      int isp = _sp_flag[spost[ist][i]];
      if (isp < 0) isp = -(isp+1);  // matched a previous node
      spost[ist][i] = isp;
    }
  }
  DELETE(xsp0);

}

void SimpleSurface::openLinearTris(IntFlag& st_flag) {
  UNUSED(st_flag);
  assert(0);
#ifdef RETHINK

  ensureNonManifoldData();

  if (nmTriData.countLinear() == 0) {
    COUT1(" > currently no linear tris detected; skipping");
    return;
  }

  COUT1(" > attempting to \"open\" " << nmTriData.countLinear() << " linear surface tris");

  int n_linear = 0;
  int n_collapsed = 0;

  for (int ist=0; ist<nst; ++ist) {
    if (nmTriData.isLinear(ist)) {
      // can't use this "open up" approach to fix collapsed node tris
      if (nmTriData.isCollapsed(ist)) {
        ++n_collapsed;
        continue;  // skip this tri
      }

      double d2_no[3] = { 0.0, 0.0, 0.0 };
      FOR_I3 {
        const double dist = DIST2(xsp[spost[ist][i]],xsp[spost[ist][(i+1)%3]]);
        d2_no[i] += dist;
        d2_no[(i+1)%3] += dist;
      }
      // now the node with the smallest sum of edge distances is the one we need to open...
      int isp_index = -1;
      if (d2_no[0] < min(d2_no[1],d2_no[2])) {
        isp_index = 0;
      }
      else if (d2_no[1] < d2_no[2]) {
        isp_index = 1;
      }
      else {
        assert(d2_no[2] <= d2_no[1]);
        isp_index = 2;
      }
      const int isp = spost[ist][isp_index];
      assert(isp >= 0 && isp < nsp);

      COUT2("      [" << n_linear << "] moving mid-node originally at " << COUT_VEC(xsp[isp]) << " towards nbrs");
      ++n_linear;

      // we are going to move isp. Figure out where to move by looking at nbrs...
      double dx[3] = { 0.0, 0.0, 0.0 };
      int count = 0;

      set<int> st_processed;
      set<int>::iterator set_it;
      stack<int> nbr_edge;
      // push both edges that touch this node onto the tri
      st_processed.insert(ist);
      nbr_edge.push(teost[ist][isp_index]);
      nbr_edge.push(teost[ist][(isp_index+2)%3]);

      while (!nbr_edge.empty()) {
        const int ied_nbr = nbr_edge.top();
        nbr_edge.pop();

        if (ied_nbr < 0) continue;  // skip if open edge

        const int ist_nbr = ied_nbr/3;
        set_it = st_processed.find(ist_nbr);
        if (set_it != st_processed.end()) continue;  // tri already processed; skip
        else st_processed.insert(ist_nbr);

        const int old_edge = ied_nbr-ist_nbr*3;
        assert(old_edge >= 0 && old_edge < 3);

        // find and process the other edge that touches the node of interest
        FOR_I3 {
          if (i != old_edge) {
            if (spost[ist_nbr][i] == isp) {
              // found the second edge
              FOR_J3 dx[j] += xsp[spost[ist_nbr][(i+1)%3]][j] - xsp[isp][j];  // add displacement
              ++count;
              nbr_edge.push(teost[ist_nbr][i]);  // add edge to stack
              break;  // only one other edge, so can exit if found
            }
            else if (spost[ist_nbr][(i+1)%3] == isp) {
              // found the second edge
              FOR_J3 dx[j] += xsp[spost[ist_nbr][i]][j] - xsp[isp][j];  // add displacement
              ++count;
              nbr_edge.push(teost[ist_nbr][i]);  // add edge to stack
              break;  // only one other edge, so can exit if found
            }
          }
        }
      }

      if (count) {
        // open the collapsed tri up in the direction of the nbrs...
        FOR_I3 xsp[isp][i] += 1.0E-6*dx[i]/double(count);
        COUT2("       > node displaced by 1.0e-6*" << COUT_VEC(dx));
        const double this_n_mod[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
        const double mag_mod = MAG(this_n_mod);
        if (mag_mod > 0.0) {
          COUT2("       > successfully opened");
          nmTriData.unsetLinear(ist);
        }
        else {
          COUT2("       > unsuccessful; translation was not enough. Run again to open more");
        }
      }
      else {
        COUT2("       > could not open up; no valid edges available to define direction");
      }
    }
  }
  if (n_collapsed) COUT2(" > " << n_collapsed << " collapsed tris were skipped");

  COUT1("Nonmanifold diagnoses after openLinearTris():");
  nmTriData.dump();

  clearNonManifoldData();

#endif

}

void SimpleSurface::splitLinearTriNeighbor() {

  ensureNonManifoldData();

  if (nmTriData.countLinear() == 0) {
    COUT1(" > currently no linear tris detected; skipping");
    return;
  }

  COUT1(" > attempting to split linear-tri neighbors for " << nmTriData.countLinear() << " linear tris");

  //int n_linear = 0;
  int n_collapsed = 0;
  int n_split = 0;

  st_flag.resize(nst);
  st_flag.setAll(0);

  for (int ist=0; ist<nst; ++ist) {
    if (nmTriData.isLinear(ist)) {
      // don't use this approach to fix collapsed node tris
      if (nmTriData.isCollapsed(ist)) {
        ++n_collapsed;
        continue;  // skip this tri
      }

      // determine who is the long edge
      int longest_edge = -1;
      double edge_d2 = 0.0;
      FOR_I3 {
        const double this_edge_d2 = DIST2(xsp[spost[ist][i]],xsp[spost[ist][(i+1)%3]]);
        if (this_edge_d2 > edge_d2) {
          edge_d2 = this_edge_d2;
          longest_edge = i;
        }
      }
      assert(longest_edge != -1);

      // determine whether to split or delete this tri
      // case 0: long edge is open -> delete
      // case 1: 2+ open edges -> delete
      // case 2: single short edge is open -> split long edge nbr
      // case 3: fully connected -> split long edge nbr
      bool b_del = false;
      bool b_split = false;
      int ist_split;
      int ist_split_edge;
      int ist_split_orient;

      if (isEdgeOpen(ist,longest_edge)) b_del = true;
      else getTriNbrData(ist_split,ist_split_edge,ist_split_orient,ist,longest_edge);  // longest edge nbr info
      //TODO may get into trouble if nbr is multi-edge...not accounted for yet

      if (!b_del) {
        // check how many shorter edges are open
        int oe_count = 0;
        FOR_I2 {
          if (isEdgeOpen(ist,(longest_edge+i)%3)) ++oe_count;
        }

        if (oe_count == 2) b_del = true;
        else b_split = true;
      }

      if (b_del) {
        assert (!b_split);
        st_flag[ist] = 1;
        nmTriData.unsetLinear(ist);
        continue;
      }

      if (b_split) {
        const int splitting_node = spost[ist][(longest_edge+2)%3];

        // get nbr info for short edges of linear tri (used to set split tri teost later)
        int nbr_ist1,nbr_i1,nbr_o1;  // nbr info for longest_edge+1
        const int nbr_ist1_bits = getTriNbrDataFull(nbr_ist1,nbr_i1,nbr_o1,ist,(longest_edge+1)%3);
        int nbr_ist2,nbr_i2,nbr_o2;  // nbr info for longest_edge+2
        const int nbr_ist2_bits = getTriNbrDataFull(nbr_ist2,nbr_i2,nbr_o2,ist,(longest_edge+2)%3);

        // new portion of split tri replaces current linear tri
        znost[ist] = znost[ist_split];
        szost[ist] = szost[ist_split];
        spost[ist][0] = splitting_node;
        spost[ist][1] = spost[ist_split][(ist_split_edge+1)%3];
        spost[ist][2] = spost[ist_split][(ist_split_edge+2)%3];
        // update my teost since I know all the local nbrs
        setTriNbrData(ist_split,(ist_split_edge+1)%3,0,ist,2);  // new edge splitting ist_split
        setTriNbrToSameAs(ist_split,(ist_split_edge+1)%3,ist,1);  // existing +1 edge of ist_split

        // split portion of ist,longest_edge
        if (ist_split_orient) {
          // long edge is misaligned, short nbr for ist is +1
          if (nbr_ist1_bits == 0) setTriNbrOpen(ist,0);
          else if ((nbr_ist1_bits == 1) || (nbr_ist1_bits == 2)) setTriNbrData(nbr_ist1,nbr_i1,nbr_o1^ist_split_orient,ist,0);
          else if (nbr_ist1_bits == 3) setTriNbrMulti(nbr_ist1,nbr_o1^ist_split_orient,ist,0);
        }
        else {
          // long edge is aligned, short nbr for ist is +2
          if (nbr_ist2_bits == 0) setTriNbrOpen(ist,0);
          else if ((nbr_ist2_bits == 1) || (nbr_ist2_bits == 2)) setTriNbrData(nbr_ist2,nbr_i2,nbr_o2,ist,0);
          else if (nbr_ist2_bits == 3) setTriNbrMulti(nbr_ist2,nbr_o2,ist,0);
        }

        // update teost for original ist_split +1 neighbor in case this tri gets hit again
        int edge_nbr,edge_i,edge_o;
        // only set info if it is a valid nbr
        if (getTriNbrData(edge_nbr,edge_i,edge_o,ist_split,(ist_split_edge+1)%3)) setTriNbrData(ist,1,edge_o,edge_nbr,edge_i);

        // adjust the tri being split
        // split the long edge by the splitting node
        spost[ist_split][(ist_split_edge+1)%3] = splitting_node;
        // update teost
        setTriNbrData(ist,2,0,ist_split,(ist_split_edge+1)%3);  // new edge splitting ist_split

        // match the short edge to the linear edge's nbr if appropriate
        if (ist_split_orient) {
          // long edge is misaligned, short nbr for ist_split is +2
          if (nbr_ist2_bits == 0) setTriNbrOpen(ist_split,ist_split_edge);
          else if ((nbr_ist2_bits == 1) || (nbr_ist2_bits == 2)) setTriNbrData(nbr_ist2,nbr_i2,nbr_o2^ist_split_orient,ist_split,ist_split_edge);
          else if (nbr_ist2_bits == 3) setTriNbrMulti(nbr_ist2,nbr_o2^ist_split_orient,ist_split,ist_split_edge);
        }
        else {
          // long edge is aligned, short nbr for ist_split is +1
          if (nbr_ist1_bits == 0) setTriNbrOpen(ist_split,ist_split_edge);
          else if ((nbr_ist1_bits == 1) || (nbr_ist1_bits == 2)) setTriNbrData(nbr_ist1,nbr_i1,nbr_o1,ist_split,ist_split_edge);
          else if (nbr_ist1_bits == 3) setTriNbrMulti(nbr_ist1,nbr_o1,ist_split,ist_split_edge);
        }

        // update teost for original le+1 neighbor
        // only set info if it is a valid nbr
        if ((nbr_ist1_bits==1) || (nbr_ist1_bits==2)) setTriNbrData(ist_split,ist_split_edge,nbr_o1^ist_split_orient,nbr_ist1,nbr_i1);

        nmTriData.unsetLinear(ist);
        ++n_split;
      }
    }
  }

  if (n_collapsed) COUT2(" > " << n_collapsed << " collapsed tris were skipped");

  if (n_split) {
    COUT2(" > " << n_split << " tris were split");
    clearNonManifoldData();
  }

  if (st_flag.count()) {
    COUT2(" > " << st_flag.count() << " tris were deleted");
    deleteFlaggedTris();  // already clears nonManifold and teost
  }

  pruneEmptyZonesAndSubzones();

}

void SimpleSurface::flipSelectedSubzones(const vector<int>& subzonesVec) {

  flagTrisFromSubzoneVec(subzonesVec);  // flagged tris in st_flag
  flipFlaggedTris();
}


void SimpleSurface::flipSelectedGroups(const vector<int>& groupsVec) {

  double volume, area;
  double gcl[3];

  vector<double> groupVolVec;
  const int ngroups = countAndIndexDisjointSurfaceGroups(volume,area,gcl,groupVolVec,false);
  // st_flag holds group index now

  IntFlag gr_flag(ngroups);
  gr_flag.setAll(0);
  for (vector<int>::const_iterator it=groupsVec.begin(); it!=groupsVec.end(); ++it) {
    gr_flag[*it] = 1;
  }
  FOR_IST {
    if (gr_flag[st_flag[ist]]) st_flag[ist] = 1;
    else st_flag[ist] = 0;
  }
  flipFlaggedTris();
}

void SimpleSurface::flipSurface() {

  st_flag.setLength(nst);
  st_flag.setAll(1);
  flipFlaggedTris();

}

void SimpleSurface::alignNormalsAuto() {

  double volume, area;
  double gcl[3];

  vector<double> groupVolVec;
  const int ngroups = countAndIndexDisjointSurfaceGroups(volume,area,gcl,groupVolVec,false);
  // st_flag holds group index now

  // find largest group and set this to positive, all else negative
  int max_group = -1;
  double maxAbsVol = -1.0;
  for (vector<double>::iterator it=groupVolVec.begin(); it!=groupVolVec.end(); ++it) {
    if (fabs(*it) > maxAbsVol) {
      maxAbsVol = fabs(*it);
      max_group = int(it-groupVolVec.begin());
    }
  }
  assert(max_group != -1);
  COUT1(" > largest abs(volume) disjoint group: " << max_group);


  IntFlag gr_flag(ngroups);
  gr_flag.setAll(0);
  for (int i=0; i<ngroups; ++i) {
    if (i == max_group) {
      gr_flag[i] = (groupVolVec[i] < 0.0) ? 1:0;
    }
    else {
      gr_flag[i] = (groupVolVec[i] < 0.0) ? 0:1;
    }
  }
  FOR_IST {
    if (gr_flag[st_flag[ist]]) st_flag[ist] = 1;
    else st_flag[ist] = 0;
  }
  flipFlaggedTris();
}

void SimpleSurface::alignNormals(const int ist_seed) {

  assert((ist_seed >= 0)&&(ist_seed < nst));
  ensureTeost();

  // put a flag in each tri indicating if it has been visited yet, and whether
  // it should be flipped...

  int * st_flag = new int[nst];
  for (int ist = 0; ist < nst; ++ist)
    st_flag[ist] = -1; // -1 means not visited yet, 0 means visited and orientation correct, 1 means visited and orientation needs flipping

  stack<int> stStack;
  st_flag[ist_seed] = 0; // 0 means leave in current orientation
  stStack.push(ist_seed);

  int flip_count = 0;
  int disjoint_count = 0;
  int ist0 = 0;
  bool done = false;
  while (!done) {

    while (!stStack.empty()) {
      const int ist = stStack.top(); stStack.pop();
      assert(st_flag[ist] != -1); // should be set
      // visit our nbrs...
      FOR_I3 {
        int ist_nbr,i_nbr,orient_nbr;
        if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
          if (st_flag[ist_nbr] == -1) {
            // this nbr has not been visited.
            if (orient_nbr == 1) {
              // nbr is flipped wrt me...
              st_flag[ist_nbr] = 1 - st_flag[ist];
            }
            else {
              assert(orient_nbr == 0);
              st_flag[ist_nbr] = st_flag[ist];
            }
            //st_flag[ist_nbr] = orient_nbr + st_flag[ist] - 2*orient_nbr*st_flag[ist];
            assert((st_flag[ist_nbr] == 0)||(st_flag[ist_nbr] == 1));
            stStack.push(ist_nbr);
          }
          else {
            // should check if multiedge here - the below assert could be caused by that...?
            //assert(st_flag[ist_nbr] == orient_nbr + st_flag[ist] - 2*orient_nbr*st_flag[ist]); // surface is mobius / folded
          }
        }
      }
    }

    ++disjoint_count;
    done = true;
    for (int ist = ist0; ist < nst; ++ist) {
      if (st_flag[ist] == -1) {
        // surface must be multiply-connected. Insert another seed and continue...
        st_flag[ist] = 0; // 0 means leave in current orientation
        stStack.push(ist);
        ist0 = ist;
        done = false;
        break;
      }
      if (st_flag[ist] == 1) {
        ++flip_count;
        // flip this tri by exchanging first 2 nodes...
        const int tmp = spost[ist][0];
        spost[ist][0] = spost[ist][1];
        spost[ist][1] = tmp;
      }
    }

  }

  delete[] st_flag;

  cout << " > flipped " << flip_count << " of " << nst << " tris in " << disjoint_count << " walkable regions." << endl;
  if (flip_count > 0) clearTeost();

  clearDynamicEdgeGroups();
  // zone metadata changed
  clearSubzoneData();
  clearZoneData();
}

void SimpleSurface::deleteLinearTris() {

  ensureNonManifoldData();

  if (nmTriData.countLinear() == 0) {
    COUT1(" > currently no linear tris detected; skipping");
    return;
  }

  clearTeost();

  st_flag.resize(nst);
  st_flag.setAll(0);

  for (int ist=0; ist<nst; ++ist) {
    if (nmTriData.isLinear(ist)) st_flag[ist] = 1;  // don't need to unset b/c mnTriData gets recomputed
  }

  deleteFlaggedTris();
  pruneEmptyZonesAndSubzones();

  clearNonManifoldData();

}

void SimpleSurface::deleteMultiNeighborTris() {

  assert(0);
#ifdef RETHINK

  ensureNonManifoldData();

  if (nmTriData.countMultiNeighbor() == 0) {
    COUT1(" > currently no multi-neighbor tris detected; skipping");
    return;
  }

  IntFlag st_del(nst);
  st_del.setAll(0);

  for (int ist=0; ist<nst; ++ist) {
    if (nmTriData.isMultiNeighbor(ist)) st_del[ist] = 1;  // don't need to unset b/c nmTriData gets recomputed in buildTeost() below
  }

  SimpleSurface::deleteFlaggedTris(st_del);

  pruneEmptyZonesAndSubzones();

#endif

}

void SimpleSurface::refineFlaggedTris(const double delta) {

  cout << "Made it! refineFlaggedTris: delta: " << delta << endl;

  ensureTeost();

  // change the meaning of st_flag to indicate which edges need to be considered for splitting
  // in the tri...

  // make sure it is 1 or zero -- this is important because we use bits...

  for (int ist = 0; ist < nst; ++ist) {
    assert((st_flag[ist] == 1)||(st_flag[ist] == 0));
  }

  for (int ist = 0; ist < nst; ++ist) {

    if (st_flag[ist] & 1) {

      // here we change st_flag to be 14 (all 3 edge bits set)...

      st_flag[ist] = 14;

      // and loop through and see if any of our nbrs are not flagged. Their
      // edges will require refinement...

      FOR_I3 {
        int ist_nbr,i_nbr,orient_nbr;
        if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
          // nbr may already be flagged, so this does nothing, but does
          // set the edges of unflagged tris for refinement...
          st_flag[ist_nbr] |= (2<<i_nbr); // note we shift 2 to start at second bit
        }
      }

    }

  }

  const int nst_old = nst;

  // used for memory management...
  int nsp_max = nsp;
  int nst_max = nst;

  map<const pair<int,int>,int> newNodeMap;

  // we change (grow) nst, so store nst_old...
  for (int ist = 0; ist < nst_old; ++ist) {

    // selected tris only...

    if (st_flag[ist] != 0) {

      // loop on the edges of this tri...

      recursivelyRefineTri(ist,delta,newNodeMap,nsp_max,nst_max);

    }

  }

  // we added new tris, so clear some stuff that might be set...

  clearTeost();
  clearStoszSzosz();

}

void SimpleSurface::recursivelyRefineTri(const int ist,const double delta,map<const pair<int,int>,int>& newNodeMap,int& nsp_max,int& nst_max) {

  // find its longest edge...
  int i_longest = -1;
  double delta_longest;
  FOR_I3 {
    // only consider edges flagged as candidates...
    if (st_flag[ist] & (2<<i)) {
      const int isp0 = spost[ist][i];
      const int isp1 = spost[ist][(i+1)%3];
      const double this_delta = DIST(xsp[isp0],xsp[isp1]);
      // this is an edge we need to consider...
      if ((i_longest == -1)||(this_delta > delta_longest)) {
        i_longest = i;
        delta_longest = this_delta;
      }
    }
  }
  assert(i_longest != -1);
  if (delta_longest > delta) {
    // split this tri at this edge...
    const int isp0 = spost[ist][i_longest];
    const int isp1 = spost[ist][(i_longest+1)%3];
    map<const pair<int,int>,int>::iterator iter = newNodeMap.find(pair<int,int>(min(isp0,isp1),max(isp0,isp1)));
    int isp_new;
    if (iter == newNodeMap.end()) {
      // no tri refinement has added this node yet, so add it...
      if (nsp == nsp_max) {
        // expand the number of points...
        nsp_max = max(nsp/4*5,nsp+100);
        // resize xsp...
        double (*xsp0)[3] = xsp;
        xsp = new double[nsp_max][3];
        for (int isp = 0; isp < nsp; ++isp) FOR_I3 xsp[isp][i] = xsp0[isp][i];
        delete[] xsp0;
      }
      assert(nsp_max > nsp);
      isp_new =  nsp++;
      FOR_I3 xsp[isp_new][i] = 0.5*(xsp[isp0][i] + xsp[isp1][i]);
      // and add to map, so we don't make it twice...
      newNodeMap[pair<int,int>(min(isp0,isp1),max(isp0,isp1))] = isp_new;
    }
    else {
      // someone already added this node...
      isp_new = iter->second;
      // we can delete, because every entry should only be used once...
      newNodeMap.erase(iter);
    }
    // we are also going to need a new tri...
    if (nst == nst_max) {
      // expand the number of tris...
      nst_max = max(nst/4*5,nst+100);
      // resize spost...
      int (*spost0)[3] = spost;
      spost = new int[nst_max][3];
      for (int ist = 0; ist < nst; ++ist) FOR_I3 spost[ist][i] = spost0[ist][i];
      delete[] spost0;
      int *znost0 = znost;
      znost = new int[nst_max];
      for (int ist = 0; ist < nst; ++ist) znost[ist] = znost0[ist];
      delete[] znost0;
      // and the IntFlag's...
      st_flag.ensureSize(nst_max); // does not change the size, just reserves and preserves current values...
      szost.ensureSize(nst_max); // does not change the size, just reserves and preserves current values...
    }
    assert(nst_max > nst);
    assert(st_flag.getMaxLength() >= nst_max);
    assert(st_flag.getLength() == nst);
    assert(szost.getMaxLength() >= nst_max);
    assert(szost.getLength() == nst);
    const int ist_new = nst++;
    st_flag.setLength(nst);
    szost.setLength(nst);
    // for the new tri...
    spost[ist_new][i_longest]       = spost[ist][i_longest];
    spost[ist_new][(i_longest+1)%3] = isp_new;
    spost[ist_new][(i_longest+2)%3] = spost[ist][(i_longest+2)%3];
    znost[ist_new] = znost[ist];
    szost[ist_new] = szost[ist];
    // reuse the original tri -- only need to change one spost...
    spost[ist][i_longest] = isp_new;
    // set the bits of the new tri...
    st_flag[ist_new] = st_flag[ist];
    // call this routine recursively for the new tris...
    recursivelyRefineTri(ist,delta,newNodeMap,nsp_max,nst_max);
    recursivelyRefineTri(ist_new,delta,newNodeMap,nsp_max,nst_max);
  }
}

void SimpleSurface::splitMultiEdges() {

  ensureTeost();

  st_flag.resize(nst);
  st_flag.setAll(0);
  sp_flag.resize(nsp);
  sp_flag.setAll(0);
  int nsp_new = nsp;
  FOR_IST {
    FOR_I3 {
      if (isEdgeMulti(ist,i)) {
        st_flag[ist] = 1;
        const int isp = spost[ist][i];
        if (sp_flag[isp] > 0)
          ++nsp_new;
        ++sp_flag[spost[ist][i]];
      }
    }
  }
  if (nsp_new == nsp) {
    cout << " no multi edges found. returning..." << endl;
    return;
  }
  sp_flag.setAll(0);
  cout << nsp_new << endl;

  growNspData(nsp_new,nsp);
  nsp_new = nsp;
  FOR_IST {
    if (st_flag[ist] == 1) {
      FOR_I3 {
        if (isEdgeMulti(ist,i)) {
          st_flag[ist] = 1;
          const int isp = spost[ist][i];
          if (sp_flag[isp] > 0) {
            FOR_J3 xsp[nsp_new][j] = xsp[isp][j];
            spost[ist][i] = nsp_new;
            ++nsp_new;
          }
          ++sp_flag[isp];
        }
      }
    }
  }

  cout << " > nsp: " << nsp << " nsp_new: " << nsp_new << endl;

  clearTeost();
  clearNonManifoldData();
  nsp = nsp_new;

}

void SimpleSurface::checkSelfIntersectionsNew() {

  // double-precision version of self-intersection checking...

  cout << "SimpleSurface::checkSelfIntersectionsNew()" << endl;

  cout << " > building adt for " << nst << " tris..." << endl;
  double (*bbmin)[3] = new double[nst][3];
  double (*bbmax)[3] = new double[nst][3];
  FOR_IST {
    FOR_I3 bbmin[ist][i] = xsp[spost[ist][0]][i];
    FOR_I3 bbmax[ist][i] = xsp[spost[ist][0]][i];
    for (int j = 1; j < 3; ++j) {
      FOR_I3 bbmin[ist][i] = min(bbmin[ist][i],xsp[spost[ist][j]][i]);
      FOR_I3 bbmax[ist][i] = max(bbmax[ist][i],xsp[spost[ist][j]][i]);
    }
  }
  Adt<double> adt(nst,bbmin,bbmax);

  int ierr = 0;
  FILE * fp = NULL;
  vector<int> intVec;
  double * x0[3];
  double * x1[3];
  int idata[6];
  double ddata[4];
  FOR_IST {
    if (ist%1000 == 0) cout << " > finished ist " << ist << " out of " << nst << endl;
    assert(intVec.empty());
    adt.buildListForBBox(intVec,bbmin[ist],bbmax[ist]);
    for (int ii = 0; ii < intVec.size(); ++ii) {
      const int ist2 = intVec[ii];
      // the tri-tri intersection algo should be consistent, so only consider a given
      // pair of tris once...
      if (ist2 > ist) {
        // look for edge-edge connections first. For these, just check the orientation and angle...
        bool b_edge = false;
        FOR_I3 {
          FOR_J3 {
            if ((spost[ist][i] == spost[ist2][(j+1)%3])&&(spost[ist][(i+1)%3] == spost[ist2][j])) {
              // this is a regular edge connection. Make sure the angle is OK...
              assert(spost[ist][(i+2)%3] != spost[ist2][(j+2)%3]);
              // TODO: angle check...
              b_edge = true;
              break;
            }
          }
          if (b_edge)
            break;
        }
        if (!b_edge) {
          // look for an flipped edge connection...
          FOR_I3 {
            FOR_J3 {
              if ((spost[ist][i] == spost[ist2][j])&&(spost[ist][(i+1)%3] == spost[ist2][(j+1)%3])) {
                // this is a flipped edge connection. Make sure the angle is OK...
                assert(spost[ist][(i+2)%3] != spost[ist2][(j+2)%3]);
                // TODO: angle check...
                cout << "Warning: got flipped edge connection" << endl;
                b_edge = true;
                break;
              }
            }
            if (b_edge)
              break;
          }
          if (!b_edge) {
            // we did not find an edge, so call the double intersection routine...
            FOR_I3 x0[i] = xsp[spost[ist][i]];
            FOR_I3 x1[i] = xsp[spost[ist2][i]];
            if (const int n = MiscUtils::calcTriTriIntersection(idata,ddata,x0,x1)) {
              if (n == -2) {
                MiscUtils::mkdir_for_file("debug/tmp.dat");
                MiscUtils::writeTriTriBinDebug(ierr,x0,x1);
                ++ierr;
              }
              else {
                if (n == 1) {
                  if ((idata[0] != MiscUtils::NODE_NODE_INT)||(spost[ist][idata[1]] != spost[ist2][idata[2]])) {
                    cout << "Warning got n: " << n << " problem a" << endl;
                  }
                  // make sure the other nodes do not match...
                  assert(spost[ist][(idata[1]+1)%3] != spost[ist2][(idata[2]+1)%3]);
                  assert(spost[ist][(idata[1]+2)%3] != spost[ist2][(idata[2]+1)%3]);
                  assert(spost[ist][(idata[1]+1)%3] != spost[ist2][(idata[2]+2)%3]);
                  assert(spost[ist][(idata[1]+2)%3] != spost[ist2][(idata[2]+2)%3]);
                }
                else if (n == 2) {
                  // these tris intersect at 2 locations...
                  double xi0[3]; MiscUtils::getTriTriX(xi0,x0,x1,idata[0],idata+1,ddata);
                  double xi1[3]; MiscUtils::getTriTriX(xi1,x0,x1,idata[3],idata+4,ddata+2);
                  if (fp == NULL) fp = fopen("edges.dat","w");
                  fprintf(fp,"%18.15le %18.15le %18.15le\n",xi0[0],xi0[1],xi0[2]);
                  fprintf(fp,"%18.15le %18.15le %18.15le\n",xi1[0],xi1[1],xi1[2]);
                }
                else {
                  assert(n == -1); // coplanar tris
                  //MiscUtils::writeTriTri(0,x0,x1);
                }
              }
            }
          }
        }
      }
    }
    intVec.clear();
  }

  if (fp) {
    cout << "Warning: surface had self-intersections. See edges.dat" << endl;
    fclose(fp);
  }

  if (ierr > 0) {
    CERR("This case not handled properly in calcTriTriIntersection. There should be one or more\n" <<
         "failing tritri.*.bin files in the debug directory. Run surfer with DEBUG_TRITRI <filename> to\n" <<
         "generate the required code, which should then be appended to ../core/common/node_node_int.hpp");
  }

  delete[] bbmin;
  delete[] bbmax;

  MPI_Pause("made it here!");

}
