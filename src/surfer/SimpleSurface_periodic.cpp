#include "SimpleSurface.hpp"
#include "Adt.hpp"
#include "Utils.hpp"
#include "WebUI.hpp"

void SimpleSurface::clearPeriodicity() {
  COUT2("clearPeriodicity()");
  //clearTeost();  // necessary?

  for (int izn=0,nzn=zoneVec.size(); izn<nzn; ++izn) {
    if (zoneVec[izn].getPeriodicBits() != 0) zoneVec[izn].setPeriodicBits(0);
  }
  PeriodicData::periodicTransformVec.clear();
  DELETE(pbi);
  b_pbi_bits = false;

  WUI(INFO,"Periodicity has been removed from this surface");
}

bool SimpleSurface::setTeostVecs(vector<uint8>& p0TeostVec,vector<uint8>& p1TeostVec) const {
  assert(p0TeostVec.empty());
  assert(p1TeostVec.empty());

  bool got_open = false;
  bool got_multi = false;
  bool got_opposite_orient = false;
  bool got_common_axis = false;

  // loop through open edges and flag nodes that are on periodic edge boundaries
  IntFlag sp_periodic(nsp);
  sp_periodic.setAll(0);
  if (pbi != NULL) {
    for (int isp = 0; isp < nsp; ++isp) {
      const int bits = int(pbi[isp]>>52);
      const int isp_ghost = pbi[isp]&MASK_52BITS;
      if (bits) {
        sp_periodic[isp] = 1;
        sp_periodic[isp_ghost] = 1;
      }
    }
  }

  for (int ist = 0; ist < nst; ++ist) {
    const int izn = znost[ist];
    if (zone_flag[izn] == 1) {
      // this is a periodic boundary tri on the main side...
      FOR_I3 {
        int ist_nbr,i_nbr,orient_nbr;
        const bool ierr = getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i);
        if (!ierr) {
          if (isEdgeOpen(ist,i)) {
            const bool b_per0 = sp_periodic[spost[ist][i]];
            const bool b_per1 = sp_periodic[spost[ist][(i+1)%3]];
            if (b_per0 && b_per1) {
              // open edge that has been set as periodic (must come first), so include
              p0TeostVec.push_back( (uint8(i)<<60) | ist );
            }
            else {
              got_open = true;
            }
          }
          else got_multi = true;
          continue;
        }
        else if (orient_nbr != 0) {
          got_opposite_orient = true;
          continue;
        }

        const int izn_nbr = znost[ist_nbr];
        if (zone_flag[izn_nbr] == 1) {
          // this is a nbr tri on the same periodic boundary -- nothing to do...
        }
        else if (zone_flag[izn_nbr] == 2) {
          // this is a nbr tri on the OTHER periodic boundary -- we allow this
          // and do not reconnect these edges (along a common axis for example)
          // but should make a note of it...
          got_common_axis = true;
        }
        else {
          assert(zone_flag[izn_nbr] == 0);
          p0TeostVec.push_back( (uint8(i_nbr)<<60) | ist_nbr ); // edge + triangle
          // no longer formally split. Leave everything connected...
        }
      }
    }
    else if (zone_flag[izn] == 2) {
      // this is a periodic boundary tri on the shadow side...
      FOR_I3 {
        int ist_nbr,i_nbr,orient_nbr;
        const bool ierr = getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i);
        if (!ierr) {
          if (isEdgeOpen(ist,i)) {
            const bool b_per0 = sp_periodic[spost[ist][i]];
            const bool b_per1 = sp_periodic[spost[ist][(i+1)%3]];
            if (b_per0 && b_per1) {
              // open edge that has been set as periodic (must come first), so include
              p1TeostVec.push_back( (uint8(i)<<60) | ist );
            }
            else {
              got_open = true;
            }
          }
          else got_multi = true;
          continue;
        }
        else if (orient_nbr != 0) {
          got_opposite_orient = true;
          continue;
        }

        const int izn_nbr = znost[ist_nbr];
        if (zone_flag[izn_nbr] == 2) {
          // this is a nbr tri on the same periodic boundary -- nothing to do...
        }
        else if (zone_flag[izn_nbr] == 1) {
          // this is a nbr tri on the OTHER periodic boundary -- we allow this
          // and do not reconnect these edges (along a common axis for example)
          // but should make a note of it...
          got_common_axis = true;
        }
        else {
          assert(zone_flag[izn_nbr] == 0);
          p1TeostVec.push_back( (uint8(i_nbr)<<60) | ist_nbr ); // edge + triangle
        }
      }
    }
    else {
      assert(zone_flag[izn] == 0);
    }
  }

  COUT2(" > p0 edges: " << p0TeostVec.size() << " p1 shadow edges: " << p1TeostVec.size() );

  if ( p0TeostVec.empty() || p1TeostVec.empty() || got_open || got_multi || got_opposite_orient || got_common_axis ) {
    // CWARN("could not set periodic: open_or_multi: " << got_open_or_multi << " opposite_orient: " << got_opposite_orient << " common_axis: " << got_common_axis );
    WUI(WARN,"unable to set periodicity; issue(s) present in surface \n" <<
        " > non-periodic open edges: " << got_open << "\n" <<
        " > multi edges: " << got_multi << "\n" <<
        " > opposite_orient: " << got_opposite_orient << "\n" <<
        " > common_axis: " << got_common_axis);
    return false;
  }

  return true;
}

void SimpleSurface::setTeostVec(vector<uint8>& pTeostVec,const vector<int>& open_edge_groups) const {
  assert(pTeostVec.empty());

  IntFlag oe_group_flag(n_open_edge_groups);
  oe_group_flag.setAll(0);
  for (int i=0, lim=open_edge_groups.size(); i<lim; ++i) {
    const int index = open_edge_groups[i];
    if (index >= 0 && (index < n_open_edge_groups)) {
      oe_group_flag[index] = 1;
    }
    else {
      WUI(WARN,"invalid loop index " << index << " specified for periodicity; ignoring");
    }
  }

  for (map<uint,int>::const_iterator it=oe_to_group.begin(); it!=oe_to_group.end(); ++it) {
    if (oe_group_flag[it->second]) {
      uint ist; uchar i;
      unpackEdge(ist,i,it->first);
      pTeostVec.push_back( (uint8(i)<<60) | ist ); // edge + triangle
    }
  }
}

bool SimpleSurface::setTeostVecs(vector<uint8>& p0TeostVec,vector<uint8>& p1TeostVec,const vector<int>& edges0,const vector<int>& edges1) const {
  assert(p0TeostVec.empty());
  assert(p1TeostVec.empty());

  setTeostVec(p0TeostVec,edges0);
  setTeostVec(p1TeostVec,edges1);

  COUT2(" > p0 edges: " << p0TeostVec.size() << " p1 shadow edges: " << p1TeostVec.size() );

  if ( p0TeostVec.empty() || p1TeostVec.empty()) {
    WUI(WARN,"unable to set periodicity; empty edge loop(s) on one or both sides");
    return false;
  }

  return true;
}

int SimpleSurface::resetPreviousPeriodic() {
  cout << "resetPreviousPeriodic()" << endl;

  // should only be in here if the bits get reset...
  assert(!b_pbi_bits);

  ensureTeost();

  vector<uint8> p0TeostVec;
  vector<uint8> p1TeostVec;
  zone_flag.setLength(zoneVec.size());

  // note that you don't do it on last transform...
  for (int ipt = 0, npt = PeriodicData::periodicTransformVec.size(); ipt < npt-1; ++ipt) {

    // first flag the zones related to this transform...
    zone_flag.setAll(0);
    for (int izn = 0, nzn = zoneVec.size(); izn < nzn; ++izn) {
      const int bits = zoneVec[izn].getPeriodicBits();
      if (bits & (1<<(2*ipt)))
        zone_flag[izn] = 1;
      else if (bits & (1<<(2*ipt+1)))
        zone_flag[izn] = 2;
    }

    // already did this once, so it should work here...
    assert(setTeostVecs(p0TeostVec,p1TeostVec));
    assert(setPeriodic(p0TeostVec,p1TeostVec,ipt,false,0.0) == 0);

    p0TeostVec.clear();
    p1TeostVec.clear();
  }

  return 0;

}

int SimpleSurface::setPeriodicOpenEdgeGroups(const vector<int>& edges0,const vector<int>& edges1,PeriodicTransform& pt) {

  cout << "setPeriodicOpenEdgesGroups()" << endl;

  // need to first ensure that the transform is either cylx or a translation without any x movement
  switch (pt.getKind()) {
  case PERIODIC_TRANSFORM_NULL:
  case PERIODIC_TRANSFORM_CYL_X:
  case PERIODIC_TRANSFORM_CART:
    break;
  default:
    WUI(WARN,"only CART and CYL_X transforms are currently supported for edge-loop periodicity");
    return 2;
  }

  // there are some flagged open edge loops that
  // need to be set as periodic

  ensureTeost();
  ensureOpenEdgeGroups();

  vector<uint8> p0TeostVec;
  vector<uint8> p1TeostVec;
  if (!setTeostVecs(p0TeostVec,p1TeostVec,edges0,edges1)) return 1;

  if (!checkLoopOrthogonalToX(p0TeostVec) || !checkLoopOrthogonalToX(p1TeostVec)) {
    WUI(WARN,"edge loops must be exactly orthogonal to X to use for periodicity");
    return 2;
  }

  const int npt0 = PeriodicData::periodicTransformVec.size();
  const int ipt = addPeriodicTransform(pt);
  assert((ipt >= 0)&&(ipt < int(PeriodicData::periodicTransformVec.size())));

  int ierr = setPeriodic(p0TeostVec,p1TeostVec,ipt,false,0.0);

  if (ierr == -1) {
    WUI(WARN,"unable to set periodicity; use SUGGEST to estimate parameters");
    // clear the transform if it was new...
    if (ipt == npt0)
      PeriodicData::periodicTransformVec.resize(npt0);
    return 1;
  }
  else if (ierr == -2) {
    // for a suggest, it must have been the top transform...
    assert(ipt == npt0);
    PeriodicData::periodicTransformVec.resize(npt0);
    return 2;  // user requested suggestion; return without clearing
  }

  if (ierr == -3) {
    return 1;
  }
  assert(ierr == 0);

  b_pbi_bits = true;  //HACK no zone bits set; do we set this state or not?
  COUT1(" > OH YEAH! periodicity set successfully (light applause)");
  return 0;  // success

}

int SimpleSurface::setPeriodicFlaggedZones(PeriodicTransform& pt,const bool b_force,const double crease_angle) {

  cout << "setPeriodicFlaggedZones()" << endl;

  // there are some flagged subzones that are touching edge loops that
  // need to be set as periodic

  ensureTeost();

  vector<uint8> p0TeostVec;
  vector<uint8> p1TeostVec;
  if (!setTeostVecs(p0TeostVec,p1TeostVec))
    return 1;

  const int npt0 = PeriodicData::periodicTransformVec.size();
  const int ipt = addPeriodicTransform(pt);
  assert((ipt >= 0)&&(ipt < int(PeriodicData::periodicTransformVec.size())));

  int ierr = setPeriodic(p0TeostVec,p1TeostVec,ipt,b_force,crease_angle);

  if (ierr == -1) {
    // CWARN("could not set periodicity; use SUGGEST to estimate required parameters");
    WUI(WARN,"unable to set periodicity; use SUGGEST to estimate parameters");
    // clear the transform if it was new...
    if (ipt == npt0)
      PeriodicData::periodicTransformVec.resize(npt0);
    return 1;
  }
  else if (ierr == -2) {
    // for a suggest, it must have been the top transform...
    assert(ipt == npt0);
    PeriodicData::periodicTransformVec.resize(npt0);
    return 2;  // user requested suggestion; return without clearing
  }

  // everything worked, so set the zoneVec's periodic_bits...
  for (int izn = 0,nzn=zoneVec.size(); izn < nzn; ++izn) {
    if (zone_flag[izn] == 1) {
      assert(zoneVec[izn].getPeriodicBits() == 0);
      zoneVec[izn].setPeriodicBits(1<<(2*ipt));
    }
    else if (zone_flag[izn] == 2) {
      assert(zoneVec[izn].getPeriodicBits() == 0);
      zoneVec[izn].setPeriodicBits(1<<(2*ipt+1));
    }
    else {
      assert(zone_flag[izn] == 0);
    }
  }

  // this gets done AFTER because is uses zone_flag...
  if (ierr == -3) {
    // this algorithm fails with multiple periodicities when nodes are added because the parents of these nodes
    // do not have their pbi updated. -3 return indicates that the periodicity was set properly, but other
    // periodicity perviously set was cleared. This call resets those periodicities (and asserts that
    // no new nodes are added in the process!)...
    ierr = resetPreviousPeriodic();
  }
  assert(ierr == 0);

  b_pbi_bits = true;
  COUT1(" > OH YEAH! periodicity set successfully (light applause)");
  return 0;  // success

}

int SimpleSurface::setPeriodic(vector<uint8>& p0TeostVec,vector<uint8>& p1TeostVec,const int ipt,const bool b_force,const double crease_angle) {

  //static int debug_counter = 0;
  //++debug_counter;

  cout << "setPeriodic()" << endl;

  // p0TeostVec and p1TeostVec are going to be split. They are on the primary, non-periodic zones.

  // build node-node and edge-node merges...

  vector<NodeNodeMerge> node0Node1MergeVec;
  vector<EdgeNodeMerge> edge0Node1MergeVec;
  vector<EdgeNodeMerge> edge1Node0MergeVec;

  int ierr;
  if (b_force) {
    ierr = setPeriodicMergeVecsForce(node0Node1MergeVec,edge0Node1MergeVec,edge1Node0MergeVec,p0TeostVec,p1TeostVec,ipt,crease_angle);
  }
  else {
    ierr = setPeriodicMergeVecs(node0Node1MergeVec,edge0Node1MergeVec,edge1Node0MergeVec,p0TeostVec,p1TeostVec,ipt);
  }

  if (ierr == -2) {
    COUT1(" > user requested a mesh suggestion");
    return -2;
  }
  else if (ierr != 0) {
    cout << "setPeriodic() Error: could not setPeriodicMergeVecs: " << ierr << endl;
    return -1;
  }

  COUT2(" > node0Node1MergeVec.size(): " << node0Node1MergeVec.size() << " edge0Node1MergeVec.size(): " << edge0Node1MergeVec.size() << " edge1Node0MergeVec.size(): " << edge1Node0MergeVec.size() );

  // now turn any edge-node merges into node-node merges by adding nodes to edges...
  // note that, unlike earlier implementations, this only splits tris with 2 or more edges
  // when those edges require splitting...

  if ((!edge0Node1MergeVec.empty())||(!edge1Node0MergeVec.empty())) {
    if (b_pbi_bits) {
      // this algorithm fails with multiple periodicities when nodes are added because the parents of these nodes
      // do not have their pbi updated. A simple fix is to redo the process starting from the first transform with
      // the updated nodes.
      ierr = -3;
      FOR_ISP pbi[isp] = uint8(isp);
    }
    splitEdgesAndAddNodeNodeMerges(node0Node1MergeVec,edge0Node1MergeVec,edge1Node0MergeVec,p0TeostVec,p1TeostVec,ipt);
  }

  /*
    cout << "debug_counter: " << debug_counter << endl;
    if (debug_counter == 4) {
    writeTecplot("debug.dat");
    }
    getchar();
  */

  // and now we have all the node-node merges, so finally put the periodic transform into
  // the periodicTransformVec. This gives us the index 0,1,or 2. This also allows the user to set part of
  // the periodicity at one time, and then more of the same periodicity later. Note that this does not
  // look to see if this transform is the inverse of another, which should be considered and combined in a
  // single transform...

  // and set the Pbi stuff...

  mergeNodesAndSetPeriodicPbi(node0Node1MergeVec,ipt);

  /*
    static int index  = 0;
    ++index;
    {
    char filename[128];
    sprintf(filename,"pbi.%08d.dat",index);
    FILE * fp = fopen(filename,"w");
    assert(pbi);
    for (int isp = 0; isp < nsp; ++isp) {
    const int bits = int(pbi[isp]>>52);
    if (bits != 0) {
    fprintf(fp,"%18.15le %18.15le %18.15le %d\n",xsp[isp][0],xsp[isp][1],xsp[isp][2],bits);
    const int isp_parent = int(pbi[isp]&MASK_52BITS);
    fprintf(fp,"%18.15le %18.15le %18.15le %d\n",xsp[isp_parent][0],xsp[isp_parent][1],xsp[isp_parent][2],bits);
    }
    }
    fclose(fp);
    cout << "take a look at pbi" << endl;
    getchar();
    }
  */

  if (ierr == -3) {
    assert(b_pbi_bits);
    b_pbi_bits = false;
    return ierr; // return of -3 means we must redo all periodicities that came before this one
  }
  else
    return 0;

}


int SimpleSurface::setPeriodicMergeVecsForce(vector<NodeNodeMerge>& node0Node1MergeVec,vector<EdgeNodeMerge>& edge0Node1MergeVec,vector<EdgeNodeMerge>& edge1Node0MergeVec,const vector<uint8>& p0TeostVec,const vector<uint8>& p1TeostVec,const int ipt,const double crease_angle) {

  // this FORCE routine is more robust for certain geometries, but does not currently support
  // suggestion.

  // TODO: confirm the weighting of the new nodes in the edge-node pairings works and
  // has selected the proper weight. if you choose a fract_threshold that is very low, then
  // only matched nodes will be connected, and unmatched with split tris along the
  // periodic boundary. It should be confirmed that these nodes match closely...

  // for grep'ing for paired new/delete
  // $> grep -e "= new" -e "delete" ~/codes/nextgen/src/surfer/SimpleSurface_periodic.cpp
  // = new delete[] START ********** setPeriodicMergeVecsForce **********

  cout << "setPeriodicForce() with crease_angle: " << crease_angle << endl;

  // For the math in this function, fract_threshold must be less than 0.25. If you make it very
  // small, then nodes are only combined when they are essentially coincident, and the maximum
  // number of edges will be split...

  const double fract_threshold = 0.2;

  // here we use the concept of corners to divide the periodic edge loops into segments
  // that start and end at corners. Once all such segments are processed, then we
  // look for any remaining segments that must be loops of some sort (i.e. no corners)...

  int * sp_flag = new int[nsp];
  for (int isp = 0; isp < nsp; ++isp)
    sp_flag[isp] = -1;

  int nlp0 = p0TeostVec.size(); // nlp == number of local pts
  int * p0Teost_flag = new int[nlp0];
  int (*spolp0)[3] = new int[nlp0][3];
  // we need a way to lookup the edge index (in the teostVec in this case) from a
  // particular local point. Normally 2 edges are associated with a particular local
  // point, but since we always walk the edges in one direction only, we can populate it
  // with just one...
  int *edolp0 = new int[nlp0];
  for (int ilp = 0; ilp < nlp0; ++ilp) {
    p0Teost_flag[ilp] = 0; // gets set to 1 when used
    FOR_I3 spolp0[ilp][i] = -1; // i-1, i, i+1
    edolp0[ilp] = -1;
  }

  int nlp_check = 0;
  double dx_edge_min = HUGE_VAL;
  for (int ii = 0; ii < p0TeostVec.size(); ++ii) {
    const int ist = int( p0TeostVec[ii] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
    const int i = int( p0TeostVec[ii]>>60 ); assert((i >= 0)&&(i < 3));
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    const double dx = DIST(xsp[isp0],xsp[isp1]);
    dx_edge_min = min(dx_edge_min,dx);
    // isp0...
    int ilp;
    if (sp_flag[isp0] == -1) {
      ilp = nlp_check++;
      assert(nlp_check <= nlp0);
      sp_flag[isp0] = ilp;
      assert(spolp0[ilp][1] == -1);
      spolp0[ilp][1] = isp0;
    }
    else {
      ilp = sp_flag[isp0];
      assert(spolp0[ilp][1] == isp0);
    }
    assert(spolp0[ilp][2] == -1);
    spolp0[ilp][2] = isp1;
    // isp1...
    if (sp_flag[isp1] == -1) {
      ilp = nlp_check++;
      assert(nlp_check <= nlp0);
      sp_flag[isp1] = ilp;
      assert(spolp0[ilp][1] == -1);
      spolp0[ilp][1] = isp1;
    }
    else {
      ilp = sp_flag[isp1];
      assert(spolp0[ilp][1] == isp1);
    }
    assert(spolp0[ilp][0] == -1);
    spolp0[ilp][0] = isp0;
    // provide a way to access the edge from ilp...
    assert(edolp0[ilp] == -1);
    edolp0[ilp] = ii;
  }
  assert(nlp_check == nlp0);

  int nlp1 = p1TeostVec.size();
  int * p1Teost_flag = new int[nlp1];
  int (*spolp1)[3] = new int[nlp1][3];
  int *edolp1 = new int[nlp1];
  for (int ilp = 0; ilp < nlp1; ++ilp) {
    p1Teost_flag[ilp] = 0;
    FOR_I3 spolp1[ilp][i] = -1; // i+1, i, i-1 NOTE: reversed
    edolp1[ilp] = -1;
  }

  // same as above, but use -2 indexing to ensure there are no common nodes...

  nlp_check = 0;
  for (int ii = 0; ii < p1TeostVec.size(); ++ii) {
    const int ist = int( p1TeostVec[ii] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
    const int i = int( p1TeostVec[ii]>>60 ); assert((i >= 0)&&(i < 3));
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    const double dx = DIST(xsp[isp0],xsp[isp1]);
    dx_edge_min = min(dx_edge_min,dx);
    // isp0...
    int ilp;
    if (sp_flag[isp0] == -1) {
      ilp = nlp_check++;
      assert(nlp_check <= nlp1);
      sp_flag[isp0] = -ilp-2;
      assert(spolp1[ilp][1] == -1);
      spolp1[ilp][1] = isp0;
    }
    else {
      ilp = -sp_flag[isp0]-2;
      assert((ilp >= 0)&&(ilp < nlp_check));
      assert(spolp1[ilp][1] == isp0);
    }
    assert(spolp1[ilp][0] == -1);
    spolp1[ilp][0] = isp1;
    // treat this edge index using the opposite sp as for 0 above
    // (opposite direction of looping)...
    assert(edolp1[ilp] == -1);
    edolp1[ilp] = ii;
    // isp1...
    if (sp_flag[isp1] == -1) {
      ilp = nlp_check++;
      assert(nlp_check <= nlp1);
      sp_flag[isp1] = -ilp-2;
      assert(spolp1[ilp][1] == -1);
      spolp1[ilp][1] = isp1;
    }
    else {
      ilp = -sp_flag[isp1]-2;
      assert((ilp >= 0)&&(ilp < nlp_check));
      assert(spolp1[ilp][1] == isp1);
    }
    assert(spolp1[ilp][2] == -1);
    spolp1[ilp][2] = isp0;
  }
  assert(nlp_check == nlp1);

  cout << " > links built successfully: no bow-ties or common points. Smallest edge length: " <<
    dx_edge_min << endl;

  // now identify "corners" on both sizes...

  const double dp_tol = cos((180.0-crease_angle)*M_PI/180.0);

  vector<pair<int,double> > corner0_vec;
  double dp_max_in = -1.0;
  double dp_min_in = 1.0;
  double dp_min_out = 1.0;
  for (int ilp = 0; ilp < nlp0; ++ilp) {
    const double dx0[3] = DIFF(xsp[spolp0[ilp][1]],xsp[spolp0[ilp][0]]);
    const double dx0_mag = MAG(dx0);
    assert(dx0_mag > 0.0);
    const double dx1[3] = DIFF(xsp[spolp0[ilp][2]],xsp[spolp0[ilp][1]]);
    const double dx1_mag = MAG(dx1);
    assert(dx1_mag > 0.0);
    const double dp = DOT_PRODUCT(dx0,dx1)/(dx0_mag*dx1_mag);
    if (dp <= dp_tol) {
      dp_max_in = max(dp_max_in,dp);
      dp_min_in = min(dp_min_in,dp);
      corner0_vec.push_back(pair<int,double>(ilp,180.0/M_PI*(1.0-acos(dp))));
    }
    else {
      dp_min_out = min(dp_min_out,dp);
    }
  }

  cout << " > side0 has " << corner0_vec.size() << " corners" << endl;
  if (!corner0_vec.empty()) {
    cout << "   > most acute corner: " << 180.0-acos(dp_min_in)*180.0/M_PI << " degrees" << endl;
    cout << "   > closest corner to crease angle: " << 180.0-acos(dp_max_in)*180.0/M_PI << " degrees" << endl;
    cout << "   > closest non-corner to crease angle: " << 180.0-acos(dp_min_out)*180.0/M_PI << " degrees" << endl;
  }

  vector<pair<int,double> > corner1_vec;
  dp_max_in = -1.0;
  dp_min_in = 1.0;
  dp_min_out = 1.0;
  for (int ilp = 0; ilp < nlp1; ++ilp) {
    const double dx0[3] = DIFF(xsp[spolp1[ilp][1]],xsp[spolp1[ilp][0]]);
    const double dx0_mag = MAG(dx0);
    assert(dx0_mag > 0.0);
    const double dx1[3] = DIFF(xsp[spolp1[ilp][2]],xsp[spolp1[ilp][1]]);
    const double dx1_mag = MAG(dx1);
    assert(dx1_mag > 0.0);
    const double dp = DOT_PRODUCT(dx0,dx1)/(dx0_mag*dx1_mag);
    if (dp <= dp_tol) {
      dp_max_in = max(dp_max_in,dp);
      dp_min_in = min(dp_min_in,dp);
      corner1_vec.push_back(pair<int,double>(ilp,180.0/M_PI*(1.0-acos(dp))));
    }
    else {
      dp_min_out = min(dp_min_out,dp);
    }
  }

  cout << " > side1 has " << corner1_vec.size() << " corners" << endl;
  if (!corner1_vec.empty()) {
    cout << "   > most acute corner: " << 180.0-acos(dp_min_in)*180.0/M_PI << " degrees" << endl;
    cout << "   > closest corner to crease angle: " << 180.0-acos(dp_max_in)*180.0/M_PI << " degrees" << endl;
    cout << "   > closest non-corner to crease angle: " << 180.0-acos(dp_min_out)*180.0/M_PI << " degrees" << endl;
  }

  if (corner0_vec.size() != corner1_vec.size()) {
    WUI(WARN,"corner counts differ. Consider changing CREASE_ANGLE if you know these surfaces should be periodic");
    delete[] edolp1;
    delete[] spolp1;
    delete[] p1Teost_flag;
    delete[] edolp0;
    delete[] spolp0;
    delete[] p0Teost_flag;
    delete[] sp_flag;
    return -1;
  }

  // now try and match the corners using BOTH angle and proximity
  // after translation...

  const int nc = corner0_vec.size();
  int * corner_flag = new int[nc];
  for (int ic = 0; ic < nc; ++ic)
    corner_flag[ic] = 0;

  vector<pair<int,double> > corner1_matched_vec;
  double d2_max = 0.0;
  double dangle_max = 0.0;
  for (int i0 = 0; i0 < nc; ++i0) {
    const int ilp0 = corner0_vec[i0].first;
    const int isp0 = spolp0[ilp0][1];
    double xsp0_t[3];
    FOR_I3 xsp0_t[i] = xsp[isp0][i];
    PeriodicData::periodicTransformVec[ipt].translate(xsp0_t);
    int i1_closest = -1;
    double d2_closest;
    for (int i1 = 0; i1 < nc; ++i1) {
      const int ilp1 = corner1_vec[i1].first;
      const int isp1 = spolp1[ilp1][1];
      const double d2 = DIST2(xsp0_t,xsp[isp1]);
      if ((i1_closest == -1)||(d2 < d2_closest)) {
	i1_closest = i1;
	d2_closest = d2;
      }
    }
    assert(i1_closest != -1);
    d2_max = max(d2_max,d2_closest);
    // make sure i1_closest is not already in the corner1_vec
    assert(corner_flag[i1_closest] == 0);
    corner_flag[i1_closest] = 1;
    corner1_matched_vec.push_back(corner1_vec[i1_closest]);
    const double dangle = fabs(corner0_vec[i0].second-corner1_vec[i1_closest].second);
    dangle_max = max(dangle_max,dangle);
  }

  for (int ic = 0; ic < nc; ++ic)
    assert(corner_flag[ic] == 1);
  delete[] corner_flag;

  if (nc > 0) {
    cout << " > all corners matched uniquely: max distance tol: " << sqrt(d2_max) <<
      ", max angle difference: " << dangle_max << " degrees." << endl;
  }

  // now loop through the links. Start by labeling the corners
  // using -1 indexing...

  for (int i = 0; i < nc; ++i) {
    const int ilp0 = corner0_vec[i].first;
    const int isp0 = spolp0[ilp0][1];
    assert((isp0 >= 0)&&(isp0 < nsp));
    spolp0[ilp0][1] = -isp0-1;
    const int ilp1 = corner1_matched_vec[i].first;
    const int isp1 = spolp1[ilp1][1];
    assert((isp1 >= 0)&&(isp1 < nsp));
    spolp1[ilp1][1] = -isp1-1;
    // merge these corners once here...
    node0Node1MergeVec.push_back(NodeNodeMerge(isp0,isp1,0.0));
  }

  for (int i = 0; i < nc; ++i) {
    double length0 = 0.0;
    int ilp0 = corner0_vec[i].first;
    int isp0 = -spolp0[ilp0][1]-1;
    assert(isp0 >= 0);
    assert(sp_flag[isp0] == ilp0);
    while (isp0 >= 0) {
      const int isp0_prev = isp0;
      const int isp0_next = spolp0[ilp0][2];
      length0 += DIST(xsp[isp0_next],xsp[isp0_prev]);
      ilp0 = sp_flag[isp0_next];
      isp0 = spolp0[ilp0][1];
      assert(isp0_next == max(isp0,-isp0-1));
    }
    // ilp0 should be at another corner (or the same if a loop)...
    int iend;
    for (iend = 0; iend < nc; ++iend) {
      if (ilp0 == corner0_vec[iend].first)
	break;
    }
    assert(iend < nc);
    cout << " > link0 from corner " << i << " to corner " << iend << " has length: " << length0 << endl;

    // and on the other side...
    double length1 = 0.0;
    int ilp1 = corner1_matched_vec[i].first;
    int isp1 = -spolp1[ilp1][1]-1;
    assert(isp1 >= 0);
    assert(sp_flag[isp1] == -ilp1-2); // recall -2 indexing
    while (isp1 >= 0) {
      const int isp1_prev = isp1;
      const int isp1_next = spolp1[ilp1][2];
      length1 += DIST(xsp[isp1_next],xsp[isp1_prev]);
      ilp1 = -sp_flag[isp1_next]-2;
      isp1 = spolp1[ilp1][1];
      assert(isp1_next == max(isp1,-isp1-1));
    }
    // ilp1 should be at another corner (or the same if a loop)...
    int iend_check;
    for (iend_check = 0; iend_check < nc; ++iend_check) {
      if (ilp1 == corner1_matched_vec[iend_check].first)
	break;
    }
    assert(iend_check < nc);
    cout << " > link1 from corner " << i << " to corner " << iend_check << " has length: " << length1 << endl;
    assert(iend_check == iend);
    cout << " > difference between link lengths: " << length1 - length0 << endl;
    // ensure lengths are within 10%...
    assert(fabs(length1-length0) < 0.1*(length0+length1));

    // -----------------------------------------------------------------------------------------------
    // now advance both cases together determining which edges get cut and which nodes get matched...
    // -----------------------------------------------------------------------------------------------

    double l0 = 0.0;
    double l1 = 0.0;
    ilp0 = corner0_vec[i].first;
    isp0 = -spolp0[ilp0][1]-1;
    assert(isp0 >= 0);
    assert(sp_flag[isp0] == ilp0);
    ilp1 = corner1_matched_vec[i].first;
    isp1 = -spolp1[ilp1][1]-1;
    assert(isp1 >= 0);
    assert(sp_flag[isp1] == -ilp1-2); // recall -2 indexing
    bool b_advance0 = true;
    bool b_advance1 = true;
    int isp0_prev,isp0_next;
    int isp1_prev,isp1_next;
    double l0_prev,l1_prev;
    while (1) {
      // advance isp0...
      if (b_advance0) {
	assert(isp0 >= 0);
	isp0_prev = isp0;
	isp0_next = spolp0[ilp0][2];
	l0_prev = l0;
	l0 += DIST(xsp[isp0_next],xsp[isp0_prev]);
	ilp0 = sp_flag[isp0_next];
	isp0 = spolp0[ilp0][1];
	assert(isp0_next == max(isp0,-isp0-1));
        // flag the edge in the teost flag as touched...
        const int ied0 = edolp0[ilp0];
        p0Teost_flag[ied0] = 1;
      }
      // advance isp1...
      if (b_advance1) {
	assert(isp1 >= 0);
	isp1_prev = isp1;
	isp1_next = spolp1[ilp1][2];
	l1_prev = l1;
	l1 += DIST(xsp[isp1_next],xsp[isp1_prev]);
	ilp1 = -sp_flag[isp1_next]-2;
	isp1 = spolp1[ilp1][1];
	assert(isp1_next == max(isp1,-isp1-1));
        // flag the edge in the teost flag as touched...
        const int ied1 = edolp1[ilp1];
        p1Teost_flag[ied1] = 1;
      }
      // we are done when we have advanced to the next corner...
      if ((isp0 < 0)&&(isp1 < 0))
        break;
      // otherwise...
      if (l0/length0 <= l1/length1) {
	// l0 is not as far along as l1. We should be able to get the
        // next distance for l0 as follows:
        assert(isp0 >= 0);
        const int isp0_next_next = spolp0[ilp0][2];
        const double dx0_next = DIST(xsp[isp0_next_next],xsp[isp0_next]);
        // we can join the nodes together if the new location of the
        // combined node results in a motion that is less than some
        // fraction of the relevant edge lengths less than 50%...
        const double fract0 = (l1*length0-l0*length1)/(2.0*length1*dx0_next);
        const double fract1 = (l1*length0-l0*length1)/(2.0*length0*(l1-l1_prev));
	if ((fract0 < fract_threshold)&&(fract1 < fract_threshold)) {
          // neither of the nodes should be corner. We have already checked
          // isp0 above...
          assert(isp1 >= 0);
          node0Node1MergeVec.push_back(NodeNodeMerge(isp0,isp1,0.0));
          b_advance0 = b_advance1 = true;
	}
	else {
          const int ied1 = edolp1[ilp1];
          {
            const int ist1 = int( p1TeostVec[ied1] & MASK_60BITS ); assert((ist1 >= 0)&&(ist1 < nst));
            const int i = int( p1TeostVec[ied1]>>60 ); assert((i >= 0)&&(i < 3));
            const int isp0_ = spost[ist1][i];
            const int isp1_ = spost[ist1][(i+1)%3];
            assert(isp1_ == isp1_prev);
            assert(isp0_ == isp1_next);
          }
          double w1 = l0/length0-l1_prev/length1;
	  double w0 = l1/length1-l0/length0;
	  double inv_sum = 1.0/(w0+w1);
	  w0 *= inv_sum;
	  w1 *= inv_sum;
	  assert((w0 > 0.0)&&(w0 < 1.0));
	  assert((w1 > 0.0)&&(w1 < 1.0));
          edge1Node0MergeVec.push_back(EdgeNodeMerge(ied1,isp0,w0,0.0)); // HEREHERE
	  b_advance0 = true;
	  b_advance1 = false;
	}
      }
      else {
	// l1 is not as far along as l0. We should be able to get the
        // next distance for l1 as follows:
        assert(isp1 >= 0);
        const int isp1_next_next = spolp1[ilp1][2];
        const double dx1_next = DIST(xsp[isp1_next_next],xsp[isp1_next]);
        // we can join the nodes together if the new location of the
        // combined node results in a motion that is less than some
        // fraction of the relevant edge lengths less than 50%...
        const double fract1 = (l0*length1-l1*length0)/(2.0*length0*dx1_next);
        const double fract0 = (l0*length1-l1*length0)/(2.0*length1*(l0-l0_prev));
	if ((fract0 < fract_threshold)&&(fract1 < fract_threshold)) {
          // neither of the nodes should be corner. We have already checked
          // isp0 above...
          assert(isp0 >= 0);
          node0Node1MergeVec.push_back(NodeNodeMerge(isp0,isp1,0.0));
	  b_advance0 = b_advance1 = true;
	}
	else {
          const int ied0 = edolp0[ilp0];
          {
            const int ist0 = int( p0TeostVec[ied0] & MASK_60BITS ); assert((ist0 >= 0)&&(ist0 < nst));
            const int i = int( p0TeostVec[ied0]>>60 ); assert((i >= 0)&&(i < 3));
            const int isp0_ = spost[ist0][i];
            const int isp1_ = spost[ist0][(i+1)%3];
            assert(isp0_ == isp0_prev);
            assert(isp1_ == isp0_next);
          }

	  double w1 = l1/length1-l0_prev/length0;
	  double w0 = l0/length0-l1/length1;
	  double inv_sum = 1.0/(w0+w1);
	  w0 *= inv_sum;
	  w1 *= inv_sum;
	  assert((w0 > 0.0)&&(w0 < 1.0));
	  assert((w1 > 0.0)&&(w1 < 1.0));
          edge0Node1MergeVec.push_back(EdgeNodeMerge(ied0,isp1,w1,0.0)); // HEREHERE
	  b_advance0 = false;
	  b_advance1 = true;
	}
      }
      //cout << "isp0: " << isp0 << " l0/length0: " << l0/length0 << " isp1: " << isp1 << " l1/length1: " << l1/length1 << endl;
    }

  }

  // ==================================================================
  // part 2: there may be some edges in p0Teost,p1Teost that had no
  // corners...
  // ==================================================================

  int n0 = 0;
  for (int ii = 0; ii < p0TeostVec.size(); ++ii) {
    if (p0Teost_flag[ii] == 0)
      ++n0;
  }
  int n1 = 0;
  for (int ii = 0; ii < p1TeostVec.size(); ++ii) {
    if (p1Teost_flag[ii] == 0)
      ++n1;
  }

  if ((n0 != 0)||(n1 != 0)) {
    cout << " > some edges were not handled in corner graph. Must be loops: " << n0 << " " << n1 << endl;

    vector<int> lp0Vec;
    vector<double> cc0Vec;
    for (int ii = 0; ii < p0TeostVec.size(); ++ii) {
      if (p0Teost_flag[ii] == 0) {
        p0Teost_flag[ii] = 2;
        // start a new loop from this untouched edge...
        const int iloop0 = lp0Vec.size();
        const int ist = int( p0TeostVec[ii] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
        const int i = int( p0TeostVec[ii]>>60 ); assert((i >= 0)&&(i < 3));
        int isp0 = spost[ist][i];
        int isp1 = spost[ist][(i+1)%3];
        const int ilp0 = sp_flag[isp0];
        assert((ilp0 >= 0)&&(ilp0 < nlp0));
        lp0Vec.push_back(ilp0);
        assert(spolp0[ilp0][1] == isp0);
        assert(spolp0[ilp0][2] == isp1);
        int ilp1 = sp_flag[isp1];
        assert((ilp1 >= 0)&&(ilp1 < nlp0));
        const int ied = edolp0[ilp1];
        assert(ied == ii);
        const double dx = DIST(xsp[isp0],xsp[isp1]);
        assert(iloop0*4 == cc0Vec.size());
        cc0Vec.push_back(dx);
        FOR_I3 cc0Vec.push_back(dx*(xsp[isp0][i]+xsp[isp1][i]));
        while (ilp1 != ilp0) {
          isp0 = isp1;
          assert(spolp0[ilp1][1] == isp1);
          isp1 = spolp0[ilp1][2];
          ilp1 = sp_flag[isp1];
          const double dx = DIST(xsp[isp0],xsp[isp1]);
          cc0Vec[iloop0*4  ] += dx;
          FOR_I3 cc0Vec[iloop0*4+1+i] += dx*(xsp[isp0][i]+xsp[isp1][i]);
          const int ied = edolp0[ilp1];
          assert(p0Teost_flag[ied] == 0);
          p0Teost_flag[ied] = 2;
        }
      }
    }

    cout << " > side0 has " << lp0Vec.size() << " loops without corners" << endl;

    vector<int> lp1Vec;
    vector<double> cc1Vec;
    for (int ii = 0; ii < p1TeostVec.size(); ++ii) {
      if (p1Teost_flag[ii] == 0) {
        p1Teost_flag[ii] = 2;
        // start a new loop from this untouched edge...
        const int iloop1 = lp1Vec.size();
        const int ist = int( p1TeostVec[ii] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
        const int i = int( p1TeostVec[ii]>>60 ); assert((i >= 0)&&(i < 3));
        int isp1 = spost[ist][i]; // note flip here
        int isp0 = spost[ist][(i+1)%3];
        const int ilp0 = -sp_flag[isp0]-2;
        assert((ilp0 >= 0)&&(ilp0 < nlp1));
        lp1Vec.push_back(ilp0);
        assert(spolp1[ilp0][1] == isp0);
        assert(spolp1[ilp0][2] == isp1);
        int ilp1 = -sp_flag[isp1]-2;
        assert((ilp1 >= 0)&&(ilp1 < nlp1));
        const int ied = edolp1[ilp1];
        assert(ied == ii);
        const double dx = DIST(xsp[isp0],xsp[isp1]);
        assert(iloop1*4 == cc1Vec.size());
        cc1Vec.push_back(dx);
        FOR_I3 cc1Vec.push_back(dx*(xsp[isp0][i]+xsp[isp1][i]));
        while (ilp1 != ilp0) {
          isp0 = isp1;
          assert(spolp1[ilp1][1] == isp1);
          isp1 = spolp1[ilp1][2];
          ilp1 = -sp_flag[isp1]-2;
          assert((ilp1 >= 0)&&(ilp1 < nlp1));
          const double dx = DIST(xsp[isp0],xsp[isp1]);
          cc1Vec[iloop1*4  ] += dx;
          FOR_I3 cc1Vec[iloop1*4+1+i] += dx*(xsp[isp0][i]+xsp[isp1][i]);
          const int ied = edolp1[ilp1];
          assert(p1Teost_flag[ied] == 0);
          p1Teost_flag[ied] = 2;
        }
      }
    }

    cout << " > side1 has " << lp1Vec.size() << " loops without corners" << endl;

    if (lp0Vec.size() != lp1Vec.size()) {
      // if the loop counts don't match, then periodicity is not possible...
      cout << " > loop counts differ. Geometry cannot be periodic" << endl;
      delete[] edolp1;
      delete[] spolp1;
      delete[] p1Teost_flag;
      delete[] edolp0;
      delete[] spolp0;
      delete[] p0Teost_flag;
      delete[] sp_flag;
      return -1;
    }

    // build a new vec called lp1MatchedVec that corresponds to lp0Vec...

    int * lp1_flag = new int[lp1Vec.size()];
    for (int ilp1 = 0; ilp1 < lp1Vec.size(); ++ilp1)
      lp1_flag[ilp1] = 0;
    vector<int> lp1MatchedVec;
    double d2_closest_max = 0.0;
    for (int iloop0 = 0; iloop0 < lp0Vec.size(); ++iloop0) {
      // grab our length...
      const double length0 = cc0Vec[iloop0*4];
      double x0_t[3]; FOR_I3 x0_t[i] = cc0Vec[iloop0*4+1+i]/(2.0*length0);
      PeriodicData::periodicTransformVec[ipt].translate(x0_t);
      // do n^2 loop on lp1Vec...
      int iloop1_closest = -1;
      double d2_closest;
      for (int iloop1 = 0; iloop1 < lp1Vec.size(); ++iloop1) {
        // grab our length...
        const double length1 = cc1Vec[iloop1*4];
        double x1[3]; FOR_I3 x1[i] = cc1Vec[iloop1*4+1+i]/(2.0*length1);
        // need a single consistent measure of "closeness"...
        const double d2 = DIST2(x1,x0_t) + (length1-length0)*(length1-length0);
        if ((iloop1_closest == -1)||(d2 < d2_closest)) {
          iloop1_closest = iloop1;
          d2_closest = d2;
        }
      }
      assert(iloop1_closest != -1);
      d2_closest_max = max(d2_closest_max,d2_closest);
      assert(lp1_flag[iloop1_closest] == 0);
      lp1_flag[iloop1_closest] = 1; // if you hit this, d2 must be an issue
      lp1MatchedVec.push_back(lp1Vec[iloop1_closest]);
    }
    delete[] lp1_flag;
    lp1Vec.clear();

    cout << " > d2_closest_max: " << d2_closest_max << endl;

    // now loop through paired ilp's and connect...

    for (int iloop = 0; iloop < lp0Vec.size(); ++iloop) {

      cout << " > working on iloop: " << iloop << endl;

      int ilp0 = lp0Vec[iloop];
      int isp0 = spolp0[ilp0][1];
      assert((isp0 >= 0)&&(isp0 < nsp));
      assert(sp_flag[isp0] == ilp0);
      double x0_t[3];
      FOR_I3 x0_t[i] = xsp[isp0][i];
      PeriodicData::periodicTransformVec[ipt].translate(x0_t);
      const double length0 = cc0Vec[iloop*4];
      int isp0_prev = spolp0[ilp0][0];
      double l0 = 0.0;
      double l0_prev = -DIST(xsp[isp0],xsp[isp0_prev]);
      int isp0_next = isp0; // rename to isp0_current

      /*
      {
        FILE * fp = fopen("p0.dat","w");
        double x_t[3];
        FOR_I3 x_t[i] = xsp[isp0][i];
        PeriodicData::periodicTransformVec[ipt].translate(x_t);
        fprintf(fp,"%18.15le %18.15le %18.15le\n",x_t[0],x_t[1],x_t[2]);
        //FOR_I3 x_t[i] = xsp[isp0_prev][i];
        //PeriodicData::periodicTransformVec[ipt].translate(x_t);
        //fprintf(fp,"%18.15le %18.15le %18.15le\n",x_t[0],x_t[1],x_t[2]);
        fclose(fp);
      }
      */

      // loop lp1 and find closest...
      int ilp1 = lp1MatchedVec[iloop];
      int isp1 = spolp1[ilp1][1];
      assert(-sp_flag[isp1]-2 == ilp1);
      int isp1_next = spolp1[ilp1][2];
      int ilp1_next = -sp_flag[isp1_next]-2;
      assert((ilp1_next >= 0)&&(ilp1_next < nlp1));
      int ilp1_closest = ilp1;
      double d2_closest = MiscUtils::getPointToEdgeDist2(x0_t,xsp[isp1],xsp[isp1_next]);
      assert(d2_closest == d2_closest);
      // recompute length1...
      double length1 = DIST(xsp[isp1],xsp[isp1_next]);
      // loop through ilp1...
      while (ilp1_next != ilp1) {
        isp1 = isp1_next;
        const int ilp1_prev = ilp1_next;
        isp1_next = spolp1[ilp1_next][2];
        ilp1_next = -sp_flag[isp1_next]-2;
        assert((ilp1_next >= 0)&&(ilp1_next < nlp1));
        const double d2 = MiscUtils::getPointToEdgeDist2(x0_t,xsp[isp1],xsp[isp1_next]);
        assert(d2 == d2);
        if (d2 < d2_closest) {
          ilp1_closest = ilp1_prev;
          d2_closest = d2;
        }
        length1 += DIST(xsp[isp1],xsp[isp1_next]);
      }

      // we now have length0 and length1, and can start the loop matching
      // at the matched point along the edge...

      ilp1 = ilp1_closest;
      isp1 = spolp1[ilp1][1];
      assert((isp1 >= 0)&&(isp1 < nsp));
      isp1_next = spolp1[ilp1][2];
      ilp1_next = -sp_flag[isp1_next]-2;
      double xc1[3];
      MiscUtils::getClosestPointOnEdgeRobust(xc1,x0_t,xsp[isp1],xsp[isp1_next]);
      // xc1 corresponds to the zero in l1, so...
      double l1 = -DIST(xsp[isp1],xc1);

      /*
      {
        FILE * fp = fopen("p1.dat","w");
        fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp1][0],xsp[isp1][1],xsp[isp1][2]);
        fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp1_next][0],xsp[isp1_next][1],xsp[isp1_next][2]);
        fclose(fp);
      }

      {
        FILE * fp = fopen("xc1.dat","w");
        fprintf(fp,"%18.15le %18.15le %18.15le\n",xc1[0],xc1[1],xc1[2]);
        fclose(fp);
      }

      cout << "TAKE A LOOK" << endl;
      getchar();
      */

      // unless we find a node-node match right away, we will advance ONLY
      // loop1 (it's l1 is currently -ve to l0's zero value)...
      bool b_advance0 = false;
      bool b_advance1 = true;

      // consider the case when we have a node-node match...
      const double dx1_next = DIST(xsp[isp1_next],xsp[isp1]);
      // we can join the nodes together if the new location of the
      // combined node results in a motion that is less than some
      // fraction of the relevant edge lengths less than 50%...
      const double fract1 = (l0*length1-l1*length0)/(2.0*length0*dx1_next);
      const double fract0 = (l0*length1-l1*length0)/(2.0*length1*(l0-l0_prev));
      if ((fract0 < fract_threshold)&&(fract1 < fract_threshold)) {
        assert(isp0 >= 0);
        node0Node1MergeVec.push_back(NodeNodeMerge(isp0,isp1,0.0));
        b_advance0 = b_advance1 = true;
        // use the corner-marking trick to stop loop...
        assert(spolp0[ilp0][1] == isp0);
        spolp0[ilp0][1] = -isp0-1;
        assert(spolp1[ilp1][1] == isp1);
        spolp1[ilp1][1] = -isp1-1;
      }

      // this code is similar to the loop code between corners...
      assert(isp0 >= 0);
      assert(sp_flag[isp0] == ilp0);
      assert(isp1 >= 0);
      assert(sp_flag[isp1] == -ilp1-2); // recall -2 indexing
      assert(b_advance1);
      int isp1_prev;
      double l1_prev;
      while (1) {

	//cout << "top of while: b_advance0: " << b_advance0 << " b_advance1: " << b_advance1 << endl;

        // advance isp0...
        if (b_advance0) {
          assert(isp0 >= 0);
          isp0_prev = isp0;
          isp0_next = spolp0[ilp0][2];
          l0_prev = l0;
          l0 += DIST(xsp[isp0_next],xsp[isp0_prev]);
          ilp0 = sp_flag[isp0_next];
          isp0 = spolp0[ilp0][1];
          assert(isp0_next == max(isp0,-isp0-1));
          // flag the edge in the teost flag as touched...
          const int ied0 = edolp0[ilp0];
          assert(p0Teost_flag[ied0] == 2);
        }
        // advance isp1...
        if (b_advance1) {
          assert(isp1 >= 0);
          isp1_prev = isp1;
          isp1_next = spolp1[ilp1][2];
          l1_prev = l1;
          l1 += DIST(xsp[isp1_next],xsp[isp1_prev]);
          ilp1 = -sp_flag[isp1_next]-2;
          isp1 = spolp1[ilp1][1];
          assert(isp1_next == max(isp1,-isp1-1));
          // flag the edge in the teost flag as touched...
          const int ied1 = edolp1[ilp1];
          assert(p1Teost_flag[ied1] == 2);
        }
        // we are done when we have advanced to the next corner...
        if ((isp0 < 0)&&(isp1 < 0))
          break;
        // otherwise...
        if (l0/length0 <= l1/length1) {
          // l0 is not as far along as l1. We should be able to get the
          // next distance for l0 as follows:
          assert(isp0 >= 0);
          const int isp0_next_next = spolp0[ilp0][2];
          const double dx0_next = DIST(xsp[isp0_next_next],xsp[isp0_next]);
          // we can join the nodes together if the new location of the
          // combined node results in a motion that is less than some
          // fraction of the relevant edge lengths less than 50%...
          const double fract0 = (l1*length0-l0*length1)/(2.0*length1*dx0_next);
          const double fract1 = (l1*length0-l0*length1)/(2.0*length0*(l1-l1_prev));
          if ((fract0 < fract_threshold)&&(fract1 < fract_threshold)) {
            // neither of the nodes should be corner. We have already checked
            // isp0 above...
            assert(isp1 >= 0);
            node0Node1MergeVec.push_back(NodeNodeMerge(isp0,isp1,0.0));
            b_advance0 = b_advance1 = true;
            // mark the nodes as handled...
            assert(spolp0[ilp0][1] == isp0);
            spolp0[ilp0][1] = -isp0-1;
            assert(spolp1[ilp1][1] == isp1);
            spolp1[ilp1][1] = -isp1-1;
          }
          else {
            const int ied1 = edolp1[ilp1];
            {
              const int ist1 = int( p1TeostVec[ied1] & MASK_60BITS ); assert((ist1 >= 0)&&(ist1 < nst));
              const int i = int( p1TeostVec[ied1]>>60 ); assert((i >= 0)&&(i < 3));
              const int isp0_ = spost[ist1][i];
              const int isp1_ = spost[ist1][(i+1)%3];
              assert(isp1_ == isp1_prev);
              assert(isp0_ == isp1_next);
            }
            double w1 = l0/length0-l1_prev/length1;
            double w0 = l1/length1-l0/length0;
            double inv_sum = 1.0/(w0+w1);
            w0 *= inv_sum;
            w1 *= inv_sum;
            assert((w0 > 0.0)&&(w0 < 1.0));
            assert((w1 > 0.0)&&(w1 < 1.0));
            edge1Node0MergeVec.push_back(EdgeNodeMerge(ied1,isp0,w0,0.0)); // HEREHERE
            b_advance0 = true;
            b_advance1 = false;
            // mark isp0 as handled...
            assert(spolp0[ilp0][1] == isp0);
            spolp0[ilp0][1] = -isp0-1;
          }
        }
        else {
          // l1 is not as far along as l0. We should be able to get the
          // next distance for l1 as follows:
          assert(isp1 >= 0);
          const int isp1_next_next = spolp1[ilp1][2];
          const double dx1_next = DIST(xsp[isp1_next_next],xsp[isp1_next]);
          // we can join the nodes together if the new location of the
          // combined node results in a motion that is less than some
          // fraction of the relevant edge lengths less than 50%...
          const double fract1 = (l0*length1-l1*length0)/(2.0*length0*dx1_next);
          const double fract0 = (l0*length1-l1*length0)/(2.0*length1*(l0-l0_prev));
          if ((fract0 < fract_threshold)&&(fract1 < fract_threshold)) {
            // neither of the nodes should be corner. We have already checked
            // isp0 above...
            assert(isp0 >= 0);
            node0Node1MergeVec.push_back(NodeNodeMerge(isp0,isp1,0.0));
            b_advance0 = b_advance1 = true;
            // mark the nodes as handled...
            assert(spolp0[ilp0][1] == isp0);
            spolp0[ilp0][1] = -isp0-1;
            assert(spolp1[ilp1][1] == isp1);
            spolp1[ilp1][1] = -isp1-1;
          }
          else {
            const int ied0 = edolp0[ilp0];
            {
              const int ist0 = int( p0TeostVec[ied0] & MASK_60BITS ); assert((ist0 >= 0)&&(ist0 < nst));
              const int i = int( p0TeostVec[ied0]>>60 ); assert((i >= 0)&&(i < 3));
              const int isp0_ = spost[ist0][i];
              const int isp1_ = spost[ist0][(i+1)%3];
              assert(isp0_ == isp0_prev);
              assert(isp1_ == isp0_next);
            }

            double w1 = l1/length1-l0_prev/length0;
            double w0 = l0/length0-l1/length1;
            double inv_sum = 1.0/(w0+w1);
            w0 *= inv_sum;
            w1 *= inv_sum;
            assert((w0 > 0.0)&&(w0 < 1.0));
            assert((w1 > 0.0)&&(w1 < 1.0));
            edge0Node1MergeVec.push_back(EdgeNodeMerge(ied0,isp1,w1,0.0)); // HEREHERE
            b_advance0 = false;
            b_advance1 = true;
            // mark isp1 as handled...
            assert(spolp1[ilp1][1] == isp1);
            spolp1[ilp1][1] = -isp1-1;
          }
        }
        //cout << "isp0: " << isp0 << " l0/length0: " << l0/length0 << " isp1: " << isp1 << " l1/length1: " << l1/length1 << endl;
        //getchar();

      }

      //cout << "!!!!!!!!!!!!!!!! COMPLETED loop: " << iloop << " !!!!!!!!!!!!!!!!!!!!" << endl;

    }

    //cout << "DONE FIRST" << endl;
    //getchar();

  }

  delete[] edolp1;
  delete[] spolp1;
  delete[] p1Teost_flag;
  delete[] edolp0;
  delete[] spolp0;
  delete[] p0Teost_flag;
  delete[] sp_flag;

  return 0;

  // = new delete[] END ********** setPeriodicMergeVecsForce **********

}

bool SimpleSurface::checkLoopOrthogonalToX(const vector<uint8>& edgeVec) const {
  double wgt0;
  double x0[3],n0[3];
  computeEdgeLoopWeightCentroidNormal(wgt0,x0,n0,edgeVec);

  if (n0[0] != 0.0) {
    WUI(WARN," > loop has component in X-direction; normal: " << COUT_VEC(n0));
    return false;
  }
  return true;
}

void SimpleSurface::computeEdgeLoopWeightCentroidNormal(double& wgt,double (&x0)[3],double (&n0)[3],const vector<uint8>& edgeVec) const {
  // compute edge-length-weighted centers of both edge loops as a check...
  FOR_I3 x0[i] = 0.0;
  FOR_I3 n0[i] = 0.0;
  wgt = 0.0;
  double * xfirst = NULL;
  const int ne0 = edgeVec.size();
  for (int ie0 = 0; ie0 < ne0; ++ie0) {
    const int ist = edgeVec[ie0] & MASK_60BITS;
    const int i = edgeVec[ie0] >> 60;
    const double this_wgt = DIST(xsp[spost[ist][i]],xsp[spost[ist][(i+1)%3]]);
    wgt += this_wgt;
    FOR_J3 x0[j] += this_wgt*(xsp[spost[ist][i]][j] + xsp[spost[ist][(i+1)%3]][j]);
    if (xfirst == NULL) {
      xfirst = xsp[spost[ist][i]];
    }
    else {
      const double this_n[3] = TRI_NORMAL_2(xfirst,xsp[spost[ist][i]],xsp[spost[ist][(i+1)%3]]);
      if (isEdgeOpen(ist,i)) {
        // open edges need to be reversed, so just reverse the normal direction...
        FOR_J3 n0[j] += 0.5*this_n[j]; // add to get outward-facing normal
      }
      else {
        FOR_J3 n0[j] -= 0.5*this_n[j]; // subtract to get outward-facing normal
      }
    }
  }
  assert(wgt > 0.0);
  FOR_I3 x0[i] /= 2.0*wgt;
}

void SimpleSurface::computeEdgeLoopWeightCentroidNormal(double& wgt,double (&x0)[3],double (&n0)[3],const vector<uint>& edgeVec) const {
  // compute edge-length-weighted centers of both edge loops as a check...
  FOR_I3 x0[i] = 0.0;
  FOR_I3 n0[i] = 0.0;
  wgt = 0.0;
  double * xfirst = NULL;
  const int ne0 = edgeVec.size();
  for (int ie0 = 0; ie0 < ne0; ++ie0) {
    uint ist; uchar i;
    unpackEdge(ist,i,edgeVec[ie0]);
    assert(isEdgeOpen(ist,i));
    const double this_wgt = DIST(xsp[spost[ist][i]],xsp[spost[ist][(i+1)%3]]);
    wgt += this_wgt;
    FOR_J3 x0[j] += this_wgt*(xsp[spost[ist][i]][j] + xsp[spost[ist][(i+1)%3]][j]);
    if (xfirst == NULL) {
      xfirst = xsp[spost[ist][i]];
    }
    else {
      const double this_n[3] = TRI_NORMAL_2(xfirst,xsp[spost[ist][i]],xsp[spost[ist][(i+1)%3]]);
      FOR_J3 n0[j] -= 0.5*this_n[j]; // subtract to get outward-facing normal
    }
  }
  assert(wgt > 0.0);
  FOR_I3 x0[i] /= 2.0*wgt;
}

void SimpleSurface::computeOpenEdgeLoopNormal(double (&n0)[3],const int igroup) const {
  // build normals for each "face"...
  FOR_I3 n0[i] = 0.0;
  double * xfirst = NULL;
  for (map<uint,int>::const_iterator it=oe_to_group.begin(); it!=oe_to_group.end(); ++it) {
    if (it->second == igroup) {
      uint ist; uchar i;
      unpackEdge(ist,i,it->first);
      assert(isEdgeOpen(ist,i));
      if (xfirst == NULL) {
        xfirst = xsp[spost[ist][i]];
      }
      else {
        const double this_n[3] = TRI_NORMAL_2(xfirst,xsp[spost[ist][i]],xsp[spost[ist][(i+1)%3]]);
        FOR_J3 n0[j] -= 0.5*this_n[j]; // subtract to get outward-facing normal
      }
    }
  }
}

int SimpleSurface::setPeriodicMergeVecs(vector<NodeNodeMerge>& node0Node1MergeVec,vector<EdgeNodeMerge>& edge0Node1MergeVec,vector<EdgeNodeMerge>& edge1Node0MergeVec,const vector<uint8>& p0TeostVec,const vector<uint8>& p1TeostVec,const int ipt) {

  // returns true if this worked, false otherwise

  double wgt0;
  double x0[3],n0[3];
  computeEdgeLoopWeightCentroidNormal(wgt0,x0,n0,p0TeostVec);

  double wgt1;
  double x1[3],n1[3];
  computeEdgeLoopWeightCentroidNormal(wgt1,x1,n1,p1TeostVec);

  COUT2(" > got x0: " << COUT_VEC(x0) << " and x1: " << COUT_VEC(x1));

  if (fabs(wgt0-wgt1) > 0.5E-2*(wgt0+wgt1)) {
    CWARN("setPeriodic edge loop length differs by more than 1%: " << wgt0 << " " << wgt1);
    WUI(WARN,"periodic edge-loops differ in length by more than 1 percent");
    return -1;
  }

  COUT2(" > average length of periodic edges: " << 0.5*(wgt0+wgt1) << ", normalized diff: " << fabs(wgt0-wgt1)/(wgt0+wgt1)*2.0);
  COUT2(" > got n0: " << COUT_VEC(n0) << " and n1: " << COUT_VEC(n1) );

  const double n0_mag = MAG(n0);
  const double n1_mag = MAG(n1);
  if (fabs(n0_mag-n1_mag) > 0.5E-2*(n0_mag+n1_mag)) {
    CWARN("setPeriodic edge loop area differs by more than 1%: " << n0_mag << " " << n1_mag);
    WUI(WARN,"periodic edge-loops differ in area by more than 1 percent");
    return -1;
  }

  COUT2(" > average area of periodic faces: " << 0.5*(n0_mag+n1_mag) << ", normalized diff: " << fabs(n0_mag-n1_mag)/(n0_mag+n1_mag)*2.0 );

  // check the transform of the edge loop centroid...
  const double rmax_tol = 1.0E-4; // normalized distance tol
  const double angle_tol = 1.0E-4; // angle tol in radians: (sin(x) = x...)

  int suggest = 0;
  if (PeriodicData::periodicTransformVec[ipt].getKind() == PERIODIC_TRANSFORM_NULL) {
    suggest = 2; // the
  }
  else {
    // we have a valid periodic transform from the user. Ensure it transforms the x and n accurately...
    double n0_t[3]; FOR_I3 n0_t[i] = n0[i];
    PeriodicData::periodicTransformVec[ipt].rotate(n0_t);
    const double dpp1 = DOT_PRODUCT(n0_t,n1)/n0_mag/n1_mag + 1.0;
    if (dpp1 > angle_tol) {
      CWARN("setPeriodic transform seems innacurate: dpp1: " << dpp1);
      suggest = 1;
    }
    // ...
    double x0_t[3]; FOR_I3 x0_t[i] = x0[i];
    PeriodicData::periodicTransformVec[ipt].translate(x0_t); // bad name. This is just rotate:
    const double dist = DIST(x0_t,x1);
    if (dist > rmax_tol*getBoundingBoxRmax()) {
      CWARN("setPeriodic transform seems innacurate: dx: [ " << x1[0]-x0_t[0] << " " << x1[1]-x0_t[1] << " " << x1[2]-x0_t[2] << " ]");
      suggest = 1;
    }
  }

  if (suggest != 0) {
    // if the edge loop normals point in opposite directions, then
    // this is cartesian periodicity...
    const double dpp1 = DOT_PRODUCT(n0,n1)/n0_mag/n1_mag + 1.0;
    COUT2(setprecision(8) << " > unit normal alignment dpp1: " << dpp1);
    if (dpp1 < angle_tol) {
      // looks like Cartesian periodicity. Tolerance the small differences to
      // zero if near zero relative to the bounding box...
      double dx[3];
      FOR_I3 {
        dx[i] = 0.0;
        if (fabs(x1[i]-x0[i]) > 1.0E-8*getBoundingBoxRmax()) dx[i] = x1[i]-x0[i];
      }
      COUT1(setprecision(8) << " > Suggested periodicity: CART " << dx[0] << " " << dx[1] << " " << dx[2]);
      UIMessage_ msg(INFO); // TODO: replace with WUI(INFO...
      msg.addLine(stringstream().flush() << "Suggested periodicity: CART " << dx[0] << " " << dx[1] << " " << dx[2]);
      msg.setClosable(true);
      WebUI::webUIOutput.add_(msg);
    }
    else {
      // otherwise this is rotational periodicity...
      double thetax=  0.0,thetay = 0.0,thetaz = 0.0;
      int count = 0;
      const double n0x_mag = sqrt(n0[1]*n0[1]+n0[2]*n0[2]);
      const double n1x_mag = sqrt(n1[1]*n1[1]+n1[2]*n1[2]);
      if ((n0x_mag > 0.5*n0_mag)&&(n1x_mag > 0.5*n1_mag)&&(fabs(n0x_mag-n1x_mag) < 0.5E-2*(n0x_mag+n1x_mag))) {
        thetax = atan2(-n0[1]*n1[2]+n0[2]*n1[1],-n0[1]*n1[1]-n0[2]*n1[2]);
        COUT1(" > thetax [degrees]: " << thetax*180.0/M_PI << " (360/" << fabs(thetax)*180.0/M_PI << "=" << 360.0/(fabs(thetax)*180.0/M_PI) << " units in a full revolution)");
        ++count;
      }
      const double n0y_mag = sqrt(n0[2]*n0[2]+n0[0]*n0[0]);
      const double n1y_mag = sqrt(n1[2]*n1[2]+n1[0]*n1[0]);
      if ((n0y_mag > 0.5*n0_mag)&&(n1y_mag > 0.5*n1_mag)&&(fabs(n0y_mag-n1y_mag) < 0.5E-2*(n0y_mag+n1y_mag))) {
        thetay = atan2(-n0[2]*n1[0]+n0[0]*n1[2],-n0[2]*n1[2]-n0[0]*n1[0]);
        COUT1(" > thetay [degrees]: " << thetay*180.0/M_PI << " (360/" << fabs(thetay)*180.0/M_PI << "=" << 360.0/(fabs(thetay)*180.0/M_PI) << " units in a full revolution)");
        ++count;
      }
      const double n0z_mag = sqrt(n0[0]*n0[0]+n0[1]*n0[1]);
      const double n1z_mag = sqrt(n1[0]*n1[0]+n1[1]*n1[1]);
      if ((n0z_mag > 0.5*n0_mag)&&(n1z_mag > 0.5*n1_mag)&&(fabs(n0z_mag-n1z_mag) < 0.5E-2*(n0z_mag+n1z_mag))) {
        thetaz = atan2(-n0[0]*n1[1]+n0[1]*n1[0],-n0[0]*n1[0]-n0[1]*n1[1]);
        COUT1(" > thetaz [degrees]: " << thetaz*180.0/M_PI << " (360/" << fabs(thetaz)*180.0/M_PI << "=" << 360.0/(fabs(thetaz)*180.0/M_PI) << " units in a full revolution)");
        ++count;
      }
      if (count != 1) {
        CWARN("Cannot determine periodicity; skipping");
        return -2;
      }
      // at this point, one of thetax, thetay, thetaz should be non-zero...
      if (thetax != 0.0) {
        assert(cos(thetax)-1.0 != 0.0);
        const double degrees = thetax*180.0/M_PI;
        // calculate the center of rotation...
        const double yc = 0.5*((x0[1]+x1[1])*(cos(thetax)-1.0)+(x1[2]-x0[2])*sin(thetax))/(cos(thetax)-1.0);
        const double zc = 0.5*((x0[2]+x1[2])*(cos(thetax)-1.0)+(x0[1]-x1[1])*sin(thetax))/(cos(thetax)-1.0);
        UIMessage_ msg(INFO); // TODO: replace with WUI(INFO,...
        msg.setClosable(true);
        if ((fabs(yc) > rmax_tol*getBoundingBoxRmax())||(fabs(zc) > rmax_tol*getBoundingBoxRmax())) {
          msg.addLine(stringstream().flush() << "Current center of rotation (0," << yc << "," << zc << ") does not support periodicity. May we suggest first applying the global translation:");
          msg.addLine(stringstream().flush() << "   TRANSLATE ALL DX 0 " << -yc << " 0 " << -zc);
          COUT1(" > NOTE: Current center of rotation (0," << yc << "," << zc << ") does not support periodicity");
          COUT1("    > first translate via: TRANSLATE ALL DX 0 " << -yc << " " << -zc << "; then");
        }
        msg.addLine("");
        msg.addLine(stringstream().flush() << "Suggested periodicity: CYL_X " << degrees);
        WebUI::webUIOutput.add_(msg);
        COUT1(setprecision(8) << " > Suggest periodicity: CYL_X " << degrees);
      }
      else if (thetay != 0.0) {
        assert(cos(thetay)-1.0 != 0.0);
        const double degrees = thetay*180.0/M_PI;
        // calculate the center of rotation...
        const double zc = 0.5*((x0[2]+x1[2])*(cos(thetay)-1.0)+(x1[0]-x0[0])*sin(thetay))/(cos(thetay)-1.0);
        const double xc = 0.5*((x0[0]+x1[0])*(cos(thetay)-1.0)+(x0[2]-x1[2])*sin(thetay))/(cos(thetay)-1.0);
        UIMessage_ msg(INFO); // TODO: replace with WUI(INFO,...
        msg.setClosable(true);
        if ((fabs(zc) > rmax_tol*getBoundingBoxRmax())||(fabs(xc) > rmax_tol*getBoundingBoxRmax())) {
          msg.addLine(stringstream().flush() << "Current center of rotation (" << -xc << " 0 " << -zc << ") does not support periodicity. May we suggest first applying the global translation:");
          msg.addLine(stringstream().flush() << "   TRANSLATE ALL DX " << -xc << " 0 " << -zc);
          COUT1(" > NOTE: Current center of rotation (" << xc << ",0," << zc << ") does not support periodicity");
          COUT1("    > first translate via: TRANSLATE ALL DX " << -xc << " 0 " << -zc << "; then");
        }
        msg.addLine("");
        msg.addLine(stringstream().flush() << "Suggested periodicity: CYL_Y " << degrees);
        WebUI::webUIOutput.add_(msg);
        COUT1(setprecision(8) << " > Suggest periodicity: CYL_Y " << degrees);
      }
      else {
        assert(thetaz != 0.0);
        assert(cos(thetaz)-1.0 != 0.0);
        const double degrees = thetaz*180.0/M_PI;
        // calculate the center of rotation...
        const double xc = 0.5*((x0[0]+x1[0])*(cos(thetaz)-1.0)+(x1[1]-x0[1])*sin(thetaz))/(cos(thetaz)-1.0);
        const double yc = 0.5*((x0[1]+x1[1])*(cos(thetaz)-1.0)+(x0[0]-x1[0])*sin(thetaz))/(cos(thetaz)-1.0);
        UIMessage_ msg(INFO); // TODO: replace
        msg.setClosable(true);
        if ((fabs(xc) > rmax_tol*getBoundingBoxRmax())||(fabs(yc) > rmax_tol*getBoundingBoxRmax())) {
          msg.addLine(stringstream().flush() << "Current center of rotation (" << xc << "," << yc << ",0) does not support periodicity. May we suggest first applying the global translation:");
          msg.addLine(stringstream().flush() << "   TRANSLATE ALL DX " << xc << "," << yc << ",0)");
          COUT1(" > NOTE: Current center of rotation (" << xc << "," << yc << ",0) does not support periodicity");
          COUT1("    > first translate via: TRANSLATE ALL DX " << -xc << " " << -yc << " 0; then");
        }
        msg.addLine("");
        msg.addLine(stringstream().flush() << "Suggested periodicity: CYL_Z " << degrees );
        WebUI::webUIOutput.add_(msg);
        COUT1(setprecision(8) << " > Suggest periodicity: CYL_Z " << degrees);
      }

    }

    // TODO: ME consider flipping this so it is consistent with UI error index 1:suggestion ,2:failure
    return -suggest; // either -2: suggestion, or -1: failure+suggestion
  }

  // now build the mergeVec's...

  // use bits to count edge touches at each node...
  sp_flag.resize(nsp);
  sp_flag.setAll(0);

  IntFlag sp_periodic(nsp);
  sp_periodic.setAll(0);
  if (pbi != NULL) {
    // if periodicity was set via edges only and not zones, these edges are oriented opposite
    // of the zone-based edges. For this reason we have to flag their nodes with
    // the opposite values, so that the sum is still the same.
    // Here we simply identify which nodes live on
    for (int isp = 0; isp < nsp; ++isp) {
      const int bits = int(pbi[isp]>>52);
      const int isp_ghost = pbi[isp]&MASK_52BITS;
      if (bits) {
        sp_periodic[isp] = 1;
        sp_periodic[isp_ghost] = 1;
      }
    }
  }

  // loop through edges and flag nodes that are on periodic edge boundaries
  const int ne0 = p0TeostVec.size();
  for (int ie0 = 0; ie0 < ne0; ++ie0) {
    const int ist = int( p0TeostVec[ie0] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
    const int i = int( p0TeostVec[ie0]>>60 ); assert((i >= 0)&&(i < 3));
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    // assert((sp_flag[spost[ist][i]]&1) == 0);
    // assert((sp_flag[spost[ist][(i+1)%3]]&2) == 0);
    if (isEdgeOpen(ist,i) && sp_periodic[isp0] && sp_periodic[isp1]) {
      sp_flag[isp0] |= 2;
      sp_flag[isp1] |= 1;
    }
    else {
      sp_flag[isp0] |= 1;
      sp_flag[isp1] |= 2;
    }

  }

  const int ne1 = p1TeostVec.size();
  for (int ie1 = 0; ie1 < ne1; ++ie1) {
    const int ist = int( p1TeostVec[ie1] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
    const int i = int( p1TeostVec[ie1]>>60 ); assert((i >= 0)&&(i < 3));
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    // assert((sp_flag[spost[ist][i]]&4) == 0);
    // assert((sp_flag[spost[ist][(i+1)%3]]&8) == 0);
    if (isEdgeOpen(ist,i) && sp_periodic[isp0] && sp_periodic[isp1]) {
      sp_flag[isp0] |= 8;
      sp_flag[isp1] |= 4;
    }
    else {
      sp_flag[isp0] |= 4;
      sp_flag[isp1] |= 8;
    }
  }

  int nsp0_ = 0;
  int nsp1_ = 0;
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] == 3) {
      sp_flag[isp] = nsp0_; // positive indexing for p0 nodes
      ++nsp0_;
    }
    else if (sp_flag[isp] == 12) {
      sp_flag[isp] = -nsp1_-2; // use -2 indexing for p1 nodes
      ++nsp1_;
    }
    else if (sp_flag[isp] == 0) {
      sp_flag[isp] = -1;
    }
    else {
      CERR("Unexpected sp_flag value: " << sp_flag[isp]); // could be an axis, figure 8, or non-loop request
      assert(0);
    }
  }

  // store the 2 nbrs associated with each isp_...

  int (*nbosp0_)[2] = new int[nsp0_][2];
  int *sposp0_ = new int[nsp0_];
  double *delta0_ = new double[nsp0_];
  for (int isp_ = 0; isp_ < nsp0_; ++isp_) {
    FOR_I2 nbosp0_[isp_][i] = -1;
    sposp0_[isp_] = -1;
    delta0_[isp_] = 0.0;
  }

  for (int ie0 = 0; ie0 < ne0; ++ie0) {
    const int ist = int( p0TeostVec[ie0] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
    const int i = int( p0TeostVec[ie0]>>60 ); assert((i >= 0)&&(i < 3));

    int i0 = i;
    int i1 = (i+1)%3;
    if (isEdgeOpen(ist,i)) {
      swap(i0,i1);
    }

    const int isp0 = spost[ist][i0];
    const int isp0_ = sp_flag[isp0]; assert((isp0_ >= 0)&&(isp0_ < nsp0_));
    const int isp1 = spost[ist][i1];
    const int isp1_ = sp_flag[isp1]; assert((isp1_ >= 0)&&(isp1_ < nsp0_));
    assert(nbosp0_[isp0_][1] == -1);
    nbosp0_[isp0_][1] = isp1;
    assert(nbosp0_[isp1_][0] == -1);
    nbosp0_[isp1_][0] = isp0;
    if (sposp0_[isp0_] == -1) sposp0_[isp0_] = isp0;
    assert(sposp0_[isp0_] == isp0);
    if (sposp0_[isp1_] == -1) sposp0_[isp1_] = isp1;
    assert(sposp0_[isp1_] == isp1);
    const double this_delta = DIST(xsp[isp0],xsp[isp1]);
    if ((delta0_[isp0_] == 0.0)||(this_delta < delta0_[isp0_])) delta0_[isp0_] = this_delta;
    if ((delta0_[isp1_] == 0.0)||(this_delta < delta0_[isp1_])) delta0_[isp1_] = this_delta;
  }

  int (*nbosp1_)[2] = new int[nsp1_][2];
  int *sposp1_ = new int[nsp1_];
  double *delta1_ = new double[nsp1_];
  for (int isp_ = 0; isp_ < nsp1_; ++isp_) {
    FOR_I2 nbosp1_[isp_][i] = -1;
    sposp1_[isp_] = -1;
    delta1_[isp_] = 0.0;
  }

  for (int ie1 = 0; ie1 < ne1; ++ie1) {
    const int ist = int( p1TeostVec[ie1] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
    const int i = int( p1TeostVec[ie1]>>60 ); assert((i >= 0)&&(i < 3));

    int i0 = i;
    int i1 = (i+1)%3;
    if (isEdgeOpen(ist,i)) {
      swap(i0,i1);
    }

    const int isp0 = spost[ist][i0];
    const int isp0_ = -sp_flag[isp0]-2; assert((isp0_ >= 0)&&(isp0_ < nsp1_));
    const int isp1 = spost[ist][i1];
    const int isp1_ = -sp_flag[isp1]-2; assert((isp1_ >= 0)&&(isp1_ < nsp1_));
    assert(nbosp1_[isp0_][1] == -1);
    nbosp1_[isp0_][1] = isp1;
    assert(nbosp1_[isp1_][0] == -1);
    nbosp1_[isp1_][0] = isp0;
    if (sposp1_[isp0_] == -1) sposp1_[isp0_] = isp0;
    assert(sposp1_[isp0_] == isp0);
    if (sposp1_[isp1_] == -1) sposp1_[isp1_] = isp1;
    assert(sposp1_[isp1_] == isp1);
    const double this_delta = DIST(xsp[isp0],xsp[isp1]);
    if ((delta1_[isp0_] == 0.0)||(this_delta < delta1_[isp0_])) delta1_[isp0_] = this_delta;
    if ((delta1_[isp1_] == 0.0)||(this_delta < delta1_[isp1_])) delta1_[isp1_] = this_delta;
  }

  COUT2(" > zip node counts: nsp0_: " << nsp0_ << ", nsp1_: " << nsp1_);

  double (*bbmin)[3] = new double[max(ne0,ne1)][3];
  double (*bbmax)[3] = new double[max(ne0,ne1)][3];
  vector<int> candidateVec;

  // ---------------------------------
  // connect p0 to 1 edges...
  // ---------------------------------

  {

    // build the edge adt for e1 edges...

    for (int ie1 = 0; ie1 < ne1; ++ie1) {
      const int ist = int( p1TeostVec[ie1] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
      const int i = int( p1TeostVec[ie1]>>60 ); assert((i >= 0)&&(i < 3));
      const int isp0 = spost[ist][i];
      const int isp1 = spost[ist][(i+1)%3];
      // set the bounding box for the edge...
      FOR_J3 bbmin[ie1][j] = min(xsp[isp0][j],xsp[isp1][j]);
      FOR_J3 bbmax[ie1][j] = max(xsp[isp0][j],xsp[isp1][j]);
    }

    Adt<double> e1Adt(ne1,bbmin,bbmax);

    // now loop on p0 points and find p1 edges...

    for (int isp_ = 0; isp_ < nsp0_; ++isp_) {

      const int isp = sposp0_[isp_]; assert((isp >= 0)&&(isp < nsp));
      double xsp_t[3]; FOR_J3 xsp_t[j] = xsp[isp][j];
      PeriodicData::periodicTransformVec[ipt].translate(xsp_t);

      // cout << "working on point: " << COUT_VEC(xsp[isp]) << " transformed to: " << COUT_VEC(xsp_t) << endl;
      // getchar();
      // {
      // FILE * fp = fopen("p.dat","w");
      // fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp][0],xsp[isp][1],xsp[isp][2]);
      // fclose(fp);
      // fp = fopen("p_t.dat","w");
      // fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_t[0],xsp_t[1],xsp_t[2]);
      // fclose(fp);
      // }

      assert(candidateVec.empty());
      assert(delta0_[isp_] > 0.0);
      e1Adt.buildListForSphere(candidateVec,xsp_t,delta0_[isp_]);
      // if we are empty, multiply delta by a factor until we find something...
      if (candidateVec.empty()) {
        double factor = 1.0;
        while (candidateVec.empty()) {
          factor *= 2.0;
          CWARN("doubling search radius to: " << delta0_[isp_]*factor << " for point at: " << COUT_VEC(xsp_t));
          e1Adt.buildListForSphere(candidateVec,xsp_t,delta0_[isp_]*factor);
        }
      }

      int isp0_min,isp1_min,ie1_min;
      double t_min;
      double d2_min = 1.0E+20;
      for (int ii = 0,ii_end=candidateVec.size(); ii < ii_end; ++ii) {
        const int ie1 = candidateVec[ii];
        // for now, allow basically anything else...
        const int ist = int( p1TeostVec[ie1] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
        const int i = int( p1TeostVec[ie1]>>60 ); assert((i >= 0)&&(i < 3));

        int i0 = i;
        int i1 = (i+1)%3;
        if (isEdgeOpen(ist,i)) {
          swap(i0,i1);
        }

        const int isp0 = spost[ist][i0];       assert(isp0 != isp);
        const int isp1 = spost[ist][i1]; assert(isp1 != isp);

        // {
        //   FILE * fp = fopen("edge.dat","w");
        //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp0][0],xsp[isp0][1],xsp[isp0][2]);
        //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp1][0],xsp[isp1][1],xsp[isp1][2]);
        //   fclose(fp);
        // }

        const double dx[3] = DIFF(xsp[isp1],xsp[isp0]);
        const double dxp[3] = DIFF(xsp_t,xsp[isp0]);
        const double dp = DOT_PRODUCT(dxp,dx);
        if (dp <= 0.0) {
          // we are closest to the first point. Note this includes the case when dx is zero, and/or dxp is zero...
          const double d2 = dxp[0]*dxp[0] + dxp[1]*dxp[1] + dxp[2]*dxp[2];
          //cout << " > case 1: got d2: " << d2 << endl;
          if (d2 < d2_min) {
            d2_min   = d2;
            isp0_min = isp0;
            isp1_min = isp1;
            ie1_min  = ie1;
            t_min    = 0.0;
          }
        }
        else {
          const double dx2 = DOT_PRODUCT(dx,dx);
          assert(dx2 > 0.0);
          if (dp >= dx2) {
            // we are closest to the second point...
            const double d2 =
              (xsp_t[0]-xsp[isp1][0])*(xsp_t[0]-xsp[isp1][0]) +
              (xsp_t[1]-xsp[isp1][1])*(xsp_t[1]-xsp[isp1][1]) +
              (xsp_t[2]-xsp[isp1][2])*(xsp_t[2]-xsp[isp1][2]);
            //cout << " > case 2: got d2: " << d2 << endl;
            if (d2 < d2_min) {
              d2_min   = d2;
              isp0_min = isp0;
              isp1_min = isp1;
              ie1_min  = ie1;
              t_min    = 1.0;
            }
          }
          else {
            // we are closest to a point along the edge...
            const double this_t = dp/dx2;
            const double d2 = max(0.0,dxp[0]*dxp[0] + dxp[1]*dxp[1] + dxp[2]*dxp[2] - dp*dp/dx2);
            if (d2 < d2_min) {
              d2_min   = d2;
              isp0_min = isp0;
              isp1_min = isp1;
              ie1_min  = ie1;
              t_min    = this_t;
            }
          }
        }
      }
      candidateVec.clear();

      assert(d2_min < 1.0E+20);

      // record this intersection...

      if (t_min == 0.0) {
        // this is a node-node merge...
        node0Node1MergeVec.push_back(NodeNodeMerge(isp,isp0_min,d2_min)); // note: put the p0 points first and p1 points second
      }
      else if (t_min == 1.0) {
        // node-node...
        node0Node1MergeVec.push_back(NodeNodeMerge(isp,isp1_min,d2_min));
      }
      else {
        // --------------------------------------------------------------------
        // this is a node-edge intersection. Use the nearby geometry to decide
        // if we can make it a node-node...
        // --------------------------------------------------------------------
        bool b_node_node = false;
        if (t_min >= 0.5) {
          // consider the consequences on the geometry of connecting isp (transformed) and isp1_min...
          const int isp1_ = -sp_flag[isp1_min]-2; assert((isp1_ >= 0)&&(isp1_ < nsp1_));
          const double d2 = DIST2(xsp_t,xsp[isp1_min]);
          if (sqrt(d2) < min(delta0_[isp_],delta1_[isp1_])) { // movement should not be larger than delta: this is dangerous...
            assert(nbosp1_[isp1_][0] == isp0_min);
            const int isp1_next = nbosp1_[isp1_][1]; assert((isp1_next >= 0)&&(isp1_next < nsp));
            // also bring over the other edge points. This requires transformation...
            double xsp_nbr_t[2][3];
            const int isp_prev = nbosp0_[isp_][0]; assert((isp_prev >= 0)&&(isp_prev < nsp));
            const int isp_next = nbosp0_[isp_][1]; assert((isp_next >= 0)&&(isp_next < nsp));
            FOR_I3 xsp_nbr_t[0][i] = xsp[isp_prev][i];
            FOR_I3 xsp_nbr_t[1][i] = xsp[isp_next][i];
            PeriodicData::periodicTransformVec[ipt].translate(xsp_nbr_t,2);
            const double angle1 = MiscUtils::calcAngle(xsp[isp0_min],xsp[isp1_min],xsp[isp1_next]);
            const double angle0 = MiscUtils::calcAngle(xsp_nbr_t[0],xsp_t,xsp_nbr_t[1]);
            // now try the same thing with the xsp merged...
            double xsp_avg[3]; FOR_I3 xsp_avg[i] = 0.5*(xsp[isp1_min][i] + xsp_t[i]);
            const double angle1_avg = MiscUtils::calcAngle(xsp[isp0_min],xsp_avg,xsp[isp1_next]);
            const double angle0_avg = MiscUtils::calcAngle(xsp_nbr_t[0],xsp_avg,xsp_nbr_t[1]);
            // if by using the average point, the segment angles only change by
            // one degree, then connect the points...
            if ((fabs(angle1-angle1_avg) < 3.0/180.0*M_PI)&&(fabs(angle0-angle0_avg) < 3.0/180.0*M_PI)) {
              node0Node1MergeVec.push_back(NodeNodeMerge(isp,isp1_min,d2));
              b_node_node = true;
            }

            // cout << "angle1: " << angle1 << endl;
            // cout << "angle0: " << angle0 << endl;
            // cout << "angle1_avg: " << angle1_avg << endl;
            // cout << "angle0_avg: " << angle0_avg << endl;
            // {
            //   FILE * fp = fopen("sp.dat","w");
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp0_min][0],xsp[isp0_min][1],xsp[isp0_min][2]);
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp1_min][0],xsp[isp1_min][1],xsp[isp1_min][2]);
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp1_next][0],xsp[isp1_next][1],xsp[isp1_next][2]);
            //   fclose(fp);
            // }
            // {
            //   FILE * fp = fopen("sp_t.dat","w");
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_nbr_t[0][0],xsp_nbr_t[0][1],xsp_nbr_t[0][2]);
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_t[0],xsp_t[1],xsp_t[2]);
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_nbr_t[1][0],xsp_nbr_t[1][1],xsp_nbr_t[1][2]);
            //   fclose(fp);
            // }

          }
        }
        else {
          // consider the consequences on the geometry of connecting isp and isp0_min...
          const int isp0_ = -sp_flag[isp0_min]-2; assert((isp0_ >= 0)&&(isp0_ < nsp1_));
          const double d2 = DIST2(xsp_t,xsp[isp0_min]);
          if (sqrt(d2) < min(delta0_[isp_],delta1_[isp0_])) {
            // movement should not be larger than delta: this is dangerous...
            assert(nbosp1_[isp0_][1] == isp1_min);
            const int isp0_prev = nbosp1_[isp0_][0]; assert((isp0_prev >= 0)&&(isp0_prev < nsp));
            // also bring over the other edge points. This requires transformation...
            double xsp_nbr_t[2][3];
            const int isp_prev = nbosp0_[isp_][0]; assert((isp_prev >= 0)&&(isp_prev < nsp));
            const int isp_next = nbosp0_[isp_][1]; assert((isp_next >= 0)&&(isp_next < nsp));
            FOR_I3 xsp_nbr_t[0][i] = xsp[isp_prev][i];
            FOR_I3 xsp_nbr_t[1][i] = xsp[isp_next][i];
            PeriodicData::periodicTransformVec[ipt].translate(xsp_nbr_t,2);
            const double angle1 = MiscUtils::calcAngle(xsp[isp0_prev],xsp[isp0_min],xsp[isp1_min]);
            const double angle0 = MiscUtils::calcAngle(xsp_nbr_t[0],xsp_t,xsp_nbr_t[1]);
            // now try the same thing with the xsp merged...
            double xsp_avg[3]; FOR_I3 xsp_avg[i] = 0.5*(xsp[isp0_min][i] + xsp_t[i]);
            const double angle1_avg = MiscUtils::calcAngle(xsp[isp0_prev],xsp_avg,xsp[isp1_min]);
            const double angle0_avg = MiscUtils::calcAngle(xsp_nbr_t[0],xsp_avg,xsp_nbr_t[1]);
            if ((fabs(angle1-angle1_avg) < 3.0/180.0*M_PI)&&(fabs(angle0-angle0_avg) < 3.0/180.0*M_PI)) {
              node0Node1MergeVec.push_back(NodeNodeMerge(isp,isp0_min,d2));
              b_node_node = true;
            }

            // cout << "angle1: " << angle1 << endl;
            // cout << "angle0: " << angle0 << endl;
            // cout << "angle1_avg: " << angle1_avg << endl;
            // cout << "angle0_avg: " << angle0_avg << endl;
            // {
            //   FILE * fp = fopen("sp.dat","w");
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp0_prev][0],xsp[isp0_prev][1],xsp[isp0_prev][2]);
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp0_min][0],xsp[isp0_min][1],xsp[isp0_min][2]);
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp1_min][0],xsp[isp1_min][1],xsp[isp1_min][2]);
            //   fclose(fp);
            // }
            // {
            //   FILE * fp = fopen("sp_t.dat","w");
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_nbr_t[0][0],xsp_nbr_t[0][1],xsp_nbr_t[0][2]);
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_t[0],xsp_t[1],xsp_t[2]);
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_nbr_t[1][0],xsp_nbr_t[1][1],xsp_nbr_t[1][2]);
            //   fclose(fp);
            // }

          }
        }
        if (!b_node_node) {
          edge1Node0MergeVec.push_back(EdgeNodeMerge(ie1_min,isp,t_min,d2_min));

          // {
          //   FILE * fp = fopen("sp.dat","w");
          //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp0_min][0],xsp[isp0_min][1],xsp[isp0_min][2]);
          //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp1_min][0],xsp[isp1_min][1],xsp[isp1_min][2]);
          //   fclose(fp);
          // }
          // {
          //   FILE * fp = fopen("sp_t.dat","w");
          //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_t[0],xsp_t[1],xsp_t[2]);
          //   fclose(fp);
          // }

        }
      }
    }

  }

  // ---------------------------------
  // connect p1 to 0 edges...
  // ---------------------------------

  {

    // build the edge adt for e0 edges...

    for (int ie0 = 0; ie0 < ne0; ++ie0) {
      const int ist = int( p0TeostVec[ie0] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
      const int i = int( p0TeostVec[ie0]>>60 ); assert((i >= 0)&&(i < 3));
      const int isp0 = spost[ist][i];
      const int isp1 = spost[ist][(i+1)%3];
      // set the bounding box for the edge...
      FOR_J3 bbmin[ie0][j] = min(xsp[isp0][j],xsp[isp1][j]);
      FOR_J3 bbmax[ie0][j] = max(xsp[isp0][j],xsp[isp1][j]);
    }

    Adt<double> e0Adt(ne0,bbmin,bbmax);

    // now loop on p1 points and find p0 edges...

    for (int isp_ = 0; isp_ < nsp1_; ++isp_) {

      const int isp = sposp1_[isp_]; assert((isp >= 0)&&(isp < nsp));
      double xsp_t[3]; FOR_J3 xsp_t[j] = xsp[isp][j];
      PeriodicData::periodicTransformVec[ipt].inv_translate(xsp_t);

      // cout << "working on point: " << COUT_VEC(xsp[isp]) << " transformed to: " << COUT_VEC(xsp_t) << endl;
      // getchar();
      // {
      //   FILE * fp = fopen("p.dat","w");
      //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp][0],xsp[isp][1],xsp[isp][2]);
      //   fclose(fp);
      //   fp = fopen("p_t.dat","w");
      //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_t[0],xsp_t[1],xsp_t[2]);
      //   fclose(fp);
      // }

      assert(candidateVec.empty());
      assert(delta1_[isp_] > 0.0);
      e0Adt.buildListForSphere(candidateVec,xsp_t,delta1_[isp_]);
      // if we are empty, multiply delta by a factor until we find something...
      if (candidateVec.empty()) {
        double factor = 1.0;
        while (candidateVec.empty()) {
          factor *= 2.0;
          CWARN("doubling search radius to: " << delta1_[isp_]*factor << " for point at: " << COUT_VEC(xsp_t));
          e0Adt.buildListForSphere(candidateVec,xsp_t,delta1_[isp_]*factor);
        }
      }

      int isp0_min,isp1_min,ie0_min;
      double t_min;
      double d2_min = 1.0E+20;
      for (int ii = 0,ii_end=candidateVec.size(); ii < ii_end; ++ii) {
        const int ie0 = candidateVec[ii];
        // for now, allow basically anything else...
        const int ist = int( p0TeostVec[ie0] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
        const int i = int( p0TeostVec[ie0]>>60 ); assert((i >= 0)&&(i < 3));
        int i0 = i;
        int i1 = (i+1)%3;
        if (isEdgeOpen(ist,i)) {
          swap(i0,i1);
        }
        const int isp0 = spost[ist][i0];       assert(isp0 != isp);
        const int isp1 = spost[ist][i1]; assert(isp1 != isp);

        // {
        //   FILE * fp = fopen("edge.dat","w");
        //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp0][0],xsp[isp0][1],xsp[isp0][2]);
        //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp1][0],xsp[isp1][1],xsp[isp1][2]);
        //   fclose(fp);
        // }

        const double dx[3] = DIFF(xsp[isp1],xsp[isp0]);
        const double dxp[3] = DIFF(xsp_t,xsp[isp0]);
        const double dp = DOT_PRODUCT(dxp,dx);
        if (dp <= 0.0) {
          // we are closest to the first point. Note this includes the case when dx is zero, and/or dxp is zero...
          const double d2 = dxp[0]*dxp[0] + dxp[1]*dxp[1] + dxp[2]*dxp[2];
          //cout << " > case 1: got d2: " << d2 << endl;
          if (d2 < d2_min) {
            d2_min   = d2;
            isp0_min = isp0;
            isp1_min = isp1;
            ie0_min  = ie0;
            t_min    = 0.0;
          }
        }
        else {
          const double dx2 = DOT_PRODUCT(dx,dx);
          assert(dx2 > 0.0);
          if (dp >= dx2) {
            // we are closest to the second point...
            const double d2 =
              (xsp_t[0]-xsp[isp1][0])*(xsp_t[0]-xsp[isp1][0]) +
              (xsp_t[1]-xsp[isp1][1])*(xsp_t[1]-xsp[isp1][1]) +
              (xsp_t[2]-xsp[isp1][2])*(xsp_t[2]-xsp[isp1][2]);
            //cout << " > case 2: got d2: " << d2 << endl;
            if (d2 < d2_min) {
              d2_min   = d2;
              isp0_min = isp0;
              isp1_min = isp1;
              ie0_min  = ie0;
              t_min    = 1.0;
            }
          }
          else {
            // we are closest to a point along the edge...
            const double this_t = dp/dx2;
            const double d2 = max(0.0,dxp[0]*dxp[0] + dxp[1]*dxp[1] + dxp[2]*dxp[2] - dp*dp/dx2);
            if (d2 < d2_min) {
              d2_min   = d2;
              isp0_min = isp0;
              isp1_min = isp1;
              ie0_min  = ie0;
              t_min    = this_t;
            }
          }
        }
      }
      candidateVec.clear();

      assert(d2_min < 1.0E+20);

      // record this intersection...

      if (t_min == 0.0) {
        // this is a node-node merge...
        node0Node1MergeVec.push_back(NodeNodeMerge(isp0_min,isp,d2_min));
      }
      else if (t_min == 1.0) {
        // node-node...
        node0Node1MergeVec.push_back(NodeNodeMerge(isp1_min,isp,d2_min));
      }
      else {
        // --------------------------------------------------------------------
        // this is a node-edge intersection. Use the nearby geometry to decide
        // if we can make it a node-node...
        // --------------------------------------------------------------------
        bool b_node_node = false;
        if (t_min >= 0.5) {
          // consider the consequences on the geometry of connecting isp (transformed) and isp1_min...
          const int isp1_ = sp_flag[isp1_min]; assert((isp1_ >= 0)&&(isp1_ < nsp0_));
          const double d2 = DIST2(xsp_t,xsp[isp1_min]);
          if (sqrt(d2) < min(delta1_[isp_],delta0_[isp1_])) { // movement should not be larger than delta: this is dangerous...
            assert(nbosp0_[isp1_][0] == isp0_min);
            const int isp1_next = nbosp0_[isp1_][1]; assert((isp1_next >= 0)&&(isp1_next < nsp));
            // also bring over the other edge points. This requires transformation...
            double xsp_nbr_t[2][3];
            const int isp_prev = nbosp1_[isp_][0]; assert((isp_prev >= 0)&&(isp_prev < nsp));
            const int isp_next = nbosp1_[isp_][1]; assert((isp_next >= 0)&&(isp_next < nsp));
            FOR_I3 xsp_nbr_t[0][i] = xsp[isp_prev][i];
            FOR_I3 xsp_nbr_t[1][i] = xsp[isp_next][i];
            PeriodicData::periodicTransformVec[ipt].inv_translate(xsp_nbr_t,2);
            const double angle1 = MiscUtils::calcAngle(xsp[isp0_min],xsp[isp1_min],xsp[isp1_next]);
            const double angle0 = MiscUtils::calcAngle(xsp_nbr_t[0],xsp_t,xsp_nbr_t[1]);
            // now try the same thing with the xsp merged...
            double xsp_avg[3]; FOR_I3 xsp_avg[i] = 0.5*(xsp[isp1_min][i] + xsp_t[i]);
            const double angle1_avg = MiscUtils::calcAngle(xsp[isp0_min],xsp_avg,xsp[isp1_next]);
            const double angle0_avg = MiscUtils::calcAngle(xsp_nbr_t[0],xsp_avg,xsp_nbr_t[1]);
            // if by using the average point, the segment angles only change by
            // one degree, then connect the points...
            if ((fabs(angle1-angle1_avg) < 3.0/180.0*M_PI)&&(fabs(angle0-angle0_avg) < 3.0/180.0*M_PI)) {
              node0Node1MergeVec.push_back(NodeNodeMerge(isp1_min,isp,d2));
              b_node_node = true;
            }
            // cout << "angle1: " << angle1 << endl;
            // cout << "angle0: " << angle0 << endl;
            // cout << "angle1_avg: " << angle1_avg << endl;
            // cout << "angle0_avg: " << angle0_avg << endl;
            // {
            //   FILE * fp = fopen("sp.dat","w");
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp0_min][0],xsp[isp0_min][1],xsp[isp0_min][2]);
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp1_min][0],xsp[isp1_min][1],xsp[isp1_min][2]);
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp1_next][0],xsp[isp1_next][1],xsp[isp1_next][2]);
            //   fclose(fp);
            // }
            // {
            //   FILE * fp = fopen("sp_t.dat","w");
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_nbr_t[0][0],xsp_nbr_t[0][1],xsp_nbr_t[0][2]);
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_t[0],xsp_t[1],xsp_t[2]);
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_nbr_t[1][0],xsp_nbr_t[1][1],xsp_nbr_t[1][2]);
            //   fclose(fp);
            // }

          }
        }
        else {
          // consider the consequences on the geometry of connecting isp and isp0_min...
          const int isp0_ = sp_flag[isp0_min]; assert((isp0_ >= 0)&&(isp0_ < nsp0_));
          const double d2 = DIST2(xsp_t,xsp[isp0_min]);
          if (sqrt(d2) < min(delta1_[isp_],delta0_[isp0_])) {
            // movement should not be larger than delta: this is dangerous...
            assert(nbosp0_[isp0_][1] == isp1_min);
            const int isp0_prev = nbosp0_[isp0_][0]; assert((isp0_prev >= 0)&&(isp0_prev < nsp));
            // also bring over the other edge points. This requires transformation...
            double xsp_nbr_t[2][3];
            const int isp_prev = nbosp1_[isp_][0]; assert((isp_prev >= 0)&&(isp_prev < nsp));
            const int isp_next = nbosp1_[isp_][1]; assert((isp_next >= 0)&&(isp_next < nsp));
            FOR_I3 xsp_nbr_t[0][i] = xsp[isp_prev][i];
            FOR_I3 xsp_nbr_t[1][i] = xsp[isp_next][i];
            PeriodicData::periodicTransformVec[ipt].inv_translate(xsp_nbr_t,2);
            const double angle1 = MiscUtils::calcAngle(xsp[isp0_prev],xsp[isp0_min],xsp[isp1_min]);
            const double angle0 = MiscUtils::calcAngle(xsp_nbr_t[0],xsp_t,xsp_nbr_t[1]);
            // now try the same thing with the xsp merged...
            double xsp_avg[3]; FOR_I3 xsp_avg[i] = 0.5*(xsp[isp0_min][i] + xsp_t[i]);
            const double angle1_avg = MiscUtils::calcAngle(xsp[isp0_prev],xsp_avg,xsp[isp1_min]);
            const double angle0_avg = MiscUtils::calcAngle(xsp_nbr_t[0],xsp_avg,xsp_nbr_t[1]);
            if ((fabs(angle1-angle1_avg) < 3.0/180.0*M_PI)&&(fabs(angle0-angle0_avg) < 3.0/180.0*M_PI)) {
              node0Node1MergeVec.push_back(NodeNodeMerge(isp0_min,isp,d2));
              b_node_node = true;
            }

            // cout << "angle1: " << angle1 << endl;
            // cout << "angle0: " << angle0 << endl;
            // cout << "angle1_avg: " << angle1_avg << endl;
            // cout << "angle0_avg: " << angle0_avg << endl;
            // {
            //   FILE * fp = fopen("sp.dat","w");
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp0_prev][0],xsp[isp0_prev][1],xsp[isp0_prev][2]);
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp0_min][0],xsp[isp0_min][1],xsp[isp0_min][2]);
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp1_min][0],xsp[isp1_min][1],xsp[isp1_min][2]);
            //   fclose(fp);
            // }
            // {
            //   FILE * fp = fopen("sp_t.dat","w");
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_nbr_t[0][0],xsp_nbr_t[0][1],xsp_nbr_t[0][2]);
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_t[0],xsp_t[1],xsp_t[2]);
            //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_nbr_t[1][0],xsp_nbr_t[1][1],xsp_nbr_t[1][2]);
            //   fclose(fp);
            // }

          }
        }
        if (!b_node_node) {
          edge0Node1MergeVec.push_back(EdgeNodeMerge(ie0_min,isp,t_min,d2_min));

          // {
          //   FILE * fp = fopen("sp.dat","w");
          //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp0_min][0],xsp[isp0_min][1],xsp[isp0_min][2]);
          //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp1_min][0],xsp[isp1_min][1],xsp[isp1_min][2]);
          //   fclose(fp);
          // }
          // {
          //   FILE * fp = fopen("sp_t.dat","w");
          //   fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_t[0],xsp_t[1],xsp_t[2]);
          //   fclose(fp);
          // }
        }
      }
    }

  }

  // = new
  delete[] bbmin;
  delete[] bbmax;

  delete[] nbosp0_;
  delete[] nbosp1_;
  delete[] sposp0_;
  delete[] sposp1_;
  delete[] delta0_;
  delete[] delta1_;

  return 0;
}

void SimpleSurface::splitEdgesAndAddNodeNodeMerges(vector<NodeNodeMerge>& node0Node1MergeVec,vector<EdgeNodeMerge>& edge0Node1MergeVec,vector<EdgeNodeMerge>& edge1Node0MergeVec,vector<uint8>& p0TeostVec,vector<uint8>& p1TeostVec,const int ipt) {

  assert( (!edge0Node1MergeVec.empty()) || (!edge1Node0MergeVec.empty()) );

  st_flag.setLength(nst);
  st_flag.setAll(0);

  sort(edge0Node1MergeVec.begin(),edge0Node1MergeVec.end());
  int ie_current = -1;
  int nst0_new = 0;
  int this_nst_new = 0;
  for (int ii = 0, limit = edge0Node1MergeVec.size(); ii < limit; ++ii) {
    if (edge0Node1MergeVec[ii].ioe != ie_current) {
      ie_current = edge0Node1MergeVec[ii].ioe;
      const int ist = int( p0TeostVec[ie_current] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
      st_flag[ist] += 1;
      this_nst_new = 1;
      const int i = int( p0TeostVec[ie_current]>>60 ); assert((i >= 0)&&(i < 3));
      int ist_nbr;
      if (getAlignedTriNbr(ist_nbr,ist,i)) {
        st_flag[ist_nbr] += 1;
        this_nst_new = 2;
      }
    }
    nst0_new += this_nst_new; // 1 or 2
  }

  sort(edge1Node0MergeVec.begin(),edge1Node0MergeVec.end());
  ie_current = -1;
  int nst1_new = 0;
  this_nst_new = 0;
  for (int ii = 0, limit = edge1Node0MergeVec.size(); ii < limit; ++ii) {
    if (edge1Node0MergeVec[ii].ioe != ie_current) {
      ie_current = edge1Node0MergeVec[ii].ioe;
      const int ist = int( p1TeostVec[ie_current] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
      st_flag[ist] += 1;
      this_nst_new = 1;
      const int i = int( p1TeostVec[ie_current]>>60 ); assert((i >= 0)&&(i < 3));
      int ist_nbr;
      if (getAlignedTriNbr(ist_nbr,ist,i)) {
        st_flag[ist_nbr] += 1;
        this_nst_new = 2;
      }
    }
    nst1_new += this_nst_new; // 1 or 2
  }

  // check if there are any tris that have 2 (or even 3?!) edges participating
  // in the edges requiring splitting. These tris need to be split...

  int nst_split_tri_new = 0;
  int nsp_split_tri_new = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] > 1) {
      st_flag[ist] = nsp_split_tri_new++; // used to index new tris and new nodes...
      nst_split_tri_new += 2; // two new tris
    }
    else {
      st_flag[ist] = -1;
    }
  }

  // we now know exactly how big everything has to be, so grow memory...

  const int nst_old = nst;
  nst = nst_old + nst_split_tri_new + nst0_new + nst1_new;
  const int nsp_old = nsp;
  nsp = nsp_old + nsp_split_tri_new + edge0Node1MergeVec.size() + edge1Node0MergeVec.size();

  // resize st and sp stuff...
  // note we will be clearing teost later...
  int (*spost_old)[3] = spost;
  spost = new int[nst][3];
  for (int ist = 0; ist < nst_old; ++ist) FOR_I3 spost[ist][i] = spost_old[ist][i];
  delete[] spost_old;

  int *znost_old = znost;
  znost = new int[nst];
  for (int ist = 0; ist < nst_old; ++ist) znost[ist] = znost_old[ist];
  delete[] znost_old;

  // also the subzone indexing...
  szost.resize(nst);

  double (*xsp_old)[3] = xsp;
  xsp = new double[nsp][3];
  for (int isp = 0; isp < nsp_old; ++isp) FOR_I3 xsp[isp][i] = xsp_old[isp][i];
  delete[] xsp_old;

  // and pbi, if it exists...
  if (pbi) {
    uint8 * pbi_old = pbi;
    pbi = new uint8[nsp];
    for (int isp = 0; isp < nsp_old; ++isp) pbi[isp] = pbi_old[isp];
    for (int isp = nsp_old; isp < nsp; ++isp) pbi[isp] = uint8(isp); // put isp in new nodes
    delete[] pbi_old;
  }

  // build new tris and nodes...

  int ist_new = nst_old;
  int isp_new = nsp_old;

  // start by splitting the tris with double edges...

  for (int ist = 0; ist < nst_old; ++ist) {
    if (st_flag[ist] >= 0) {
      assert(ist_new == nst_old+2*st_flag[ist]);
      assert(isp_new == nsp_old+st_flag[ist]);
      const int isp0 = spost[ist][0];
      const int isp1 = spost[ist][1];
      const int isp2 = spost[ist][2];
      FOR_I3 xsp[isp_new][i] = (xsp[isp0][i] + xsp[isp1][i] + xsp[isp2][i])/3.0;
      // split the original tri into three tris such that the outer edges of
      // these tris are consistent with the original edge ordering...
      spost[ist][2] = isp_new;
      // the first new tri...
      spost[ist_new][0] = isp_new;
      spost[ist_new][1] = isp1;
      spost[ist_new][2] = isp2;
      znost[ist_new] = znost[ist];
      szost[ist_new] = szost[ist];
      // the second new tri...
      spost[ist_new+1][0] = isp0;
      spost[ist_new+1][1] = isp_new;
      spost[ist_new+1][2] = isp2;
      znost[ist_new+1] = znost[ist];
      szost[ist_new+1] = szost[ist];
      // increment the new counts...
      isp_new += 1;
      ist_new += 2;
    }
  }

  assert(isp_new == nsp_old + nsp_split_tri_new);
  assert(ist_new == nst_old + nst_split_tri_new);

  // also a new node and 1 or 2 new tris are required for each edge...

  ie_current = -1;
  double x0[3],x1[3];
  for (int ii = 0, limit = edge0Node1MergeVec.size(); ii < limit; ++ii) {
    const int ie = edge0Node1MergeVec[ii].ioe;
    const int ist_orig = int( p0TeostVec[ie] & MASK_60BITS ); assert((ist_orig >= 0)&&(ist_orig < nst));
    const int i = int( p0TeostVec[ie]>>60 ); assert((i >= 0)&&(i < 3));
    // if ist_orig is a tri that has been duplicated, then adjust ist to the duplicate. Note
    // that because of the way we oriented the new tris, there is no need to change the
    // edge access pattern. Just modify the ist to the appropriate new tri...
    int ist = ist_orig;
    if ((st_flag[ist] >= 0)&&(i > 0)) {
      ist = nst_old+2*st_flag[ist]+i-1; // i = 1,2
    }
    if (ie != ie_current) {
      // on the first time we touch an edge, grap the end coordinates. These are used
      // to compute all the new nodes...
      ie_current = ie;
      const int isp0 = spost[ist][i];
      FOR_J3 x0[j] = xsp[isp0][j];
      const int isp1 = spost[ist][(i+1)%3];
      FOR_J3 x1[j] = xsp[isp1][j];

    }
    // the new node is at position "t"...
    FOR_J3 xsp[isp_new][j] = (1.0-edge0Node1MergeVec[ii].t)*x0[j] + edge0Node1MergeVec[ii].t*x1[j];
    node0Node1MergeVec.push_back(NodeNodeMerge(isp_new,edge0Node1MergeVec[ii].isp,edge0Node1MergeVec[ii].d2));

    // {
    //   // check proximity to the transformed point...
    //   double xsp_t[3]; FOR_J3 xsp_t[j] = xsp[edge0Node1MergeVec[ii].isp][j];
    //   periodicTransformVec[ipt].inv_translate(xsp_t);
    //   cout << " SECOND d2: " << DIST2(xsp_t,xsp[isp_new]) << " edge0Node1MergeVec[ii].d2: " << edge0Node1MergeVec[ii].d2 << endl;
    //   getchar();
    // }

    // the new tri goes in where the old tri was...
    spost[ist_new][i]       = spost[ist][i];
    spost[ist_new][(i+1)%3] = isp_new;
    spost[ist_new][(i+2)%3] = spost[ist][(i+2)%3]; // same internal node
    // the new node becomes isp0 of the original tri...
    spost[ist][i] = isp_new; // do this AFTER above
    // new zone is the same zone...
    znost[ist_new] = znost[ist];
    // new subzone is the same...
    szost[ist_new] = szost[ist];
    ++ist_new;
    // now check the edge neighbor if it exists...
    int ist_nbr,i_nbr,orient_nbr;
    if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist_orig,i)) {
      assert(orient_nbr == 0); // oriented properly for now
      // similar to above, grab the correct split tri if split...
      if ((st_flag[ist_nbr] >= 0)&&(i_nbr > 0)) {
	ist_nbr = nst_old+2*st_flag[ist_nbr]+i_nbr-1; // i = 1,2
      }
      assert(spost[ist_nbr][i_nbr] == spost[ist][(i+1)%3]);
      // the 2nd new tri...
      spost[ist_new][(i_nbr+1)%3] = spost[ist_nbr][(i_nbr+1)%3];
      spost[ist_new][i_nbr]       = isp_new;
      spost[ist_new][(i_nbr+2)%3] = spost[ist_nbr][(i_nbr+2)%3];
      // and readjust ist_nbr...
      spost[ist_nbr][(i_nbr+1)%3] = isp_new; // do this AFTER above
      // new zone is the same zone...
      znost[ist_new] = znost[ist_nbr];
      // new subzone is the same...
      szost[ist_new] = szost[ist_nbr];
      ++ist_new;
    }
    // and finally increment the node...
    ++isp_new;
  }
  edge0Node1MergeVec.clear();

  ie_current = -1;
  for (int ii = 0, limit = edge1Node0MergeVec.size(); ii < limit; ++ii) {
    const int ie = edge1Node0MergeVec[ii].ioe;
    const int ist_orig = int( p1TeostVec[ie] & MASK_60BITS ); assert((ist_orig >= 0)&&(ist_orig < nst));
    const int i = int( p1TeostVec[ie]>>60 ); assert((i >= 0)&&(i < 3));
    // if ist_orig is a tri that has been duplicated, then adjust ist to the duplicate. Note
    // that because of the way we oriented the new tris, there is no need to change the
    // edge access pattern. Just modify the ist to the appropriate new tri...
    int ist = ist_orig;
    if ((st_flag[ist] >= 0)&&(i > 0)) {
      ist = nst_old+2*st_flag[ist]+i-1; // i = 1,2
    }
    if (ie != ie_current) {
      // on the first time we touch an edge, grap the end coordinates. These are used
      // to compute all the new nodes...
      ie_current = ie;
      const int isp0 = spost[ist][i];
      FOR_J3 x0[j] = xsp[isp0][j];
      const int isp1 = spost[ist][(i+1)%3];
      FOR_J3 x1[j] = xsp[isp1][j];
    }
    // the new node is at position "t"...
    FOR_J3 xsp[isp_new][j] = (1.0-edge1Node0MergeVec[ii].t)*x0[j] + edge1Node0MergeVec[ii].t*x1[j];
    node0Node1MergeVec.push_back(NodeNodeMerge(edge1Node0MergeVec[ii].isp,isp_new,edge1Node0MergeVec[ii].d2));

    // {
    //   // check proximity to the transformed point...
    //   double xsp_t[3]; FOR_J3 xsp_t[j] = xsp[edge1Node0MergeVec[ii].isp][j];
    //   periodicTransformVec[ipt].translate(xsp_t);
    //   cout << " SECOND d2: " << DIST2(xsp_t,xsp[isp_new]) << " edge1Node0MergeVec[ii].d2: " << edge1Node0MergeVec[ii].d2 << endl;
    //   getchar();
    // }

    // the new tri goes in where the old tri was...
    spost[ist_new][i]       = spost[ist][i];
    spost[ist_new][(i+1)%3] = isp_new;
    spost[ist_new][(i+2)%3] = spost[ist][(i+2)%3]; // same internal node
    // the new node becomes isp0 of the original tri...
    spost[ist][i] = isp_new;
    // new zone is the same zone...
    znost[ist_new] = znost[ist];
    // new subzone is the same...
    szost[ist_new] = szost[ist];
    ++ist_new;
    // now check the edge neighbor if it exists...
    int ist_nbr,i_nbr,orient_nbr;
    if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist_orig,i)) {
      assert(orient_nbr == 0); // oriented properly for now
      // similar to above, grab the correct split tri if split...
      if ((st_flag[ist_nbr] >= 0)&&(i_nbr > 0)) {
	ist_nbr = nst_old+2*st_flag[ist_nbr]+i_nbr-1; // i = 1,2
      }
      assert(spost[ist_nbr][i_nbr] == spost[ist][(i+1)%3]);
      // the 2nd new tri...
      spost[ist_new][(i_nbr+1)%3] = spost[ist_nbr][(i_nbr+1)%3];
      spost[ist_new][i_nbr]       = isp_new;
      spost[ist_new][(i_nbr+2)%3] = spost[ist_nbr][(i_nbr+2)%3];
      // and readjust ist_nbr...
      spost[ist_nbr][(i_nbr+1)%3] = isp_new;
      // new zone is the same zone...
      znost[ist_new] = znost[ist_nbr];
      // new subzone is the same...
      szost[ist_new] = szost[ist_nbr];
      ++ist_new;
    }
    // and finally increment the node...
    ++isp_new;
  }
  edge1Node0MergeVec.clear();

  assert(isp_new == nsp);
  assert(ist_new == nst);

  // even though we used teost above to get to edge-connected tris of both
  // p0TeostVec and p1TeostVec, it is now broken/incomplete because of
  // new nodes and tris...

  clearTeost();

}

int SimpleSurface::addPeriodicTransform(PeriodicTransform& pt) {

  // returns the index of the periodic transform.

  int ipt=0;
  for (int npt=PeriodicData::periodicTransformVec.size(); ipt < npt; ++ipt) {
    if (PeriodicData::periodicTransformVec[ipt].isEqual(pt)) return ipt;
  }
  PeriodicData::periodicTransformVec.push_back(pt);
  return ipt;

}

void SimpleSurface::mergeNodesAndSetPeriodicPbi(vector<NodeNodeMerge>& node0Node1MergeVec,const int ipt) {

  // some notes: when setting multiple periodicity, it is not necessary to set
  // pbi on split edges where both tris are periodic. In general periodicity
  // should work even if splitting is required to match nodes.

  COUT1(" > mergeNodesAndSetPeriodicPbi: ipt: " << ipt );

  // if we don't have a pbi, build it for the first time...
  if (pbi == NULL) {
    pbi = new uint8[nsp];
    for (int isp = 0; isp < nsp; ++isp)
      pbi[isp] = uint8(isp);
  }

  // convert existing pbi into sp_bits and sp_flag...
  // pbi should ONLY contain odd bits?
  unsigned char * sp_bits = new unsigned char[nsp];
  sp_flag.setLength(nsp);
  for (int isp = 0; isp < nsp; ++isp) {
    sp_bits[isp] = 0;
    const int bits = int(pbi[isp]>>52); // use the top 12 places for bits -- allows for double bits x 3 pairs
    FOR_I3 {
      if (bits&(1<<(2*i))) sp_bits[isp] |= (1<<i); // should be just 0,2,4 bits
      assert((bits&(1<<(2*i+1))) == 0); // never any odd bits in pbi
    }
    sp_flag[isp] = (pbi[isp]&MASK_52BITS);
  }

  // now push down connections...
  for (int ii = 0, limit = node0Node1MergeVec.size(); ii < limit; ++ii) {
    const int isp0 = node0Node1MergeVec[ii].isp0;
    const int isp1 = node0Node1MergeVec[ii].isp1;
    // set the bits on the isp1 side...
    sp_bits[isp1] |= (1<<ipt);
    // now pound down the heirarchy...
    int isp0_ = sp_flag[isp0];
    while (sp_flag[isp0_] != isp0_)
      isp0_ = sp_flag[isp0_];
    int isp1_ = sp_flag[isp1];
    while (sp_flag[isp1_] != isp1_)
      isp1_ = sp_flag[isp1_];
    sp_flag[isp0_] = sp_flag[isp1_] = min(isp0_,isp1_);
  }

  // now set the new points...
  double (*xsp2)[3] = new double[nsp][3];

  int nsp_new = 0;
  int * nsp_merged = new int[nsp]; // nsp_new <= nsp
  for (int isp = 0; isp < nsp; ++isp) {
    nsp_merged[isp] = 0; // for checking
    if (sp_flag[isp] == isp) {
      const int isp_new = nsp_new++;
      nsp_merged[isp_new] = 1;
      if (sp_bits[isp] == 0) {
        FOR_I3 xsp2[isp_new][i] = xsp[isp][i]*xsp[isp][i];
        FOR_I3 xsp[isp_new][i] = xsp[isp][i];
      }
      else {
        double xsp_t[3];
        FOR_I3 xsp_t[i] = xsp[isp][i];
        FOR_I3 if (sp_bits[isp] & (1<<i)) PeriodicData::periodicTransformVec[i].inv_translate(xsp_t);
        FOR_I3 xsp2[isp_new][i] = xsp_t[i]*xsp_t[i];
        FOR_I3 xsp[isp_new][i] = xsp_t[i];
      }
      sp_flag[isp] = -isp_new-1; // -1 indexing
    }
    else {
      int isp_ = sp_flag[isp]; assert(isp_ >= 0);
      while (sp_flag[isp_] >= 0)
        isp_ = sp_flag[isp_];
      sp_flag[isp] = sp_flag[isp_];
      const int isp_new = -sp_flag[isp_]-1;
      assert(nsp_merged[isp_new] >= 1);
      assert(isp_new < isp);
      if (sp_bits[isp] == 0) {
        FOR_I3 xsp2[isp_new][i] += xsp[isp][i]*xsp[isp][i];
        FOR_I3 xsp[isp_new][i] += xsp[isp][i];
      }
      else {
        double xsp_t[3];
        FOR_I3 xsp_t[i] = xsp[isp][i];
        FOR_I3 if (sp_bits[isp] & (1<<i)) {
          assert(i < int(PeriodicData::periodicTransformVec.size()));
          PeriodicData::periodicTransformVec[i].inv_translate(xsp_t);
        }
        FOR_I3 xsp2[isp_new][i] += xsp_t[i]*xsp_t[i];
        FOR_I3 xsp[isp_new][i] += xsp_t[i];
      }
      ++nsp_merged[isp_new];
    }
  }

  // normalize xsp...
  double rms_max = 0.0;
  for (int isp_new = 0; isp_new < nsp_new; ++isp_new) {
    assert(nsp_merged[isp_new] >= 1);
    if (nsp_merged[isp_new] > 1) {
      const double inv_n = 1.0/(double)nsp_merged[isp_new];
      FOR_I3 xsp[isp_new][i] *= inv_n;
      FOR_I3 {
        const double rms = sqrt(max(0.0,xsp2[isp_new][i]*inv_n - xsp[isp_new][i]*xsp[isp_new][i]));
        rms_max = max(rms_max,rms);
      }
    }
    pbi[isp_new] = isp_new; // these first nodes are not on odd periodic zones
  }

  COUT2(" > max rms point motion: absolute: " << rms_max << ", normalized: " << rms_max/getBoundingBoxRmax() );

  delete[] xsp2;
  delete[] nsp_merged;

  map<const pair<unsigned char,int>,int> pMap;
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_bits[isp]) {
      const int isp_new = -sp_flag[isp]-1;
      assert((isp_new >= 0)&&(isp_new < nsp_new));
      map<const pair<unsigned char,int>,int>::iterator iter = pMap.find(pair<unsigned char,int>(sp_bits[isp],isp_new));
      if (iter == pMap.end()) {
        const int isp_new_p = nsp_new++;
        pMap[pair<unsigned char,int>(sp_bits[isp],isp_new)] = isp_new_p;
        FOR_I3 xsp[isp_new_p][i] = xsp[isp_new][i];
        pbi[isp_new_p] = isp_new;
        FOR_I3 if (sp_bits[isp]&(1<<i)) {
          PeriodicData::periodicTransformVec[i].translate(xsp[isp_new_p]);
          pbi[isp_new_p] |= (uint8(1)<<(52+2*i)); // set the bit associated with the tranformation to this point: 0,2,4
        }
        sp_flag[isp] = -isp_new_p-1;
      }
      else {
        sp_flag[isp] = -iter->second-1;
      }
    }
  }

  delete[] sp_bits;

  // the sp_flag contains a -1 indexing of the new sp's...
  int collapsed_tri_count = 0;
  int nst_new = 0;
  FOR_IST {
    FOR_I3 {
      const int isp = spost[ist][i];
      spost[ist][i] = -sp_flag[isp]-1;
    }
    if ((spost[ist][0] != spost[ist][1])&&(spost[ist][0] != spost[ist][2])&&(spost[ist][1] != spost[ist][2])) {
      const int ist_new = nst_new++;
      if (ist_new != ist) {
        assert(ist_new < ist);
        FOR_I3 spost[ist_new][i] = spost[ist][i];
        znost[ist_new] = znost[ist];
        szost[ist_new] = szost[ist];
      }
    }
    else {
      ++collapsed_tri_count;
    }
  }

  // resize...
  assert(nsp_new <= nsp); nsp = nsp_new;
  assert(nst_new <= nst); nst = nst_new;

  // leave szozn_i

  // clear teost if tris were collapsed...
  if (collapsed_tri_count > 0) {
    COUT1(" > node merging created " << collapsed_tri_count << " collapsed tris");
    clearTeost();
  }


}

// ================================================================================
// this zip routine shares similarity with the periodic force routine, so
// it is put here...
// ================================================================================

int SimpleSurface::zipOpenEdgesChtFlaggedZones(const double crease_angle) {

  // ======START====== new delete[] SimpleSurface::zipOpenEdgesChtFlaggedZones ======START=======

  cout << "SimpleSurface::zipOpenEdgesChtFlaggedZones, crease_angle: " << crease_angle << endl;

  ensureTeost();

  vector<uint8> p0TeostVec;
  vector<uint8> p1TeostVec;
  for (int ist = 0; ist < nst; ++ist) {
    const int izn = znost[ist];
    if (zone_flag[izn] == 1) {
      FOR_I3 {
        if (isEdgeOpen(ist,i)) {
          p0TeostVec.push_back( (uint8(i)<<60) | ist );
        }
      }
    }
    else if (zone_flag[izn] == 2) {
      FOR_I3 {
        if (isEdgeOpen(ist,i)) {
          p1TeostVec.push_back( (uint8(i)<<60) | ist );
        }
      }
    }
  }

  cout << " > p0TeostVec.size(): " << p0TeostVec.size() << endl;
  cout << " > p1TeostVec.size(): " << p1TeostVec.size() << endl;

  const double fract_threshold = 0.2;

  // here we use the concept of corners to divide the periodic edge loops into segments
  // that start and end at corners. Once all such segments are processed, then we
  // look for any remaining segments that must be loops of some sort (i.e. no corners)...

  int * sp_flag = new int[nsp];
  for (int isp = 0; isp < nsp; ++isp)
    sp_flag[isp] = -1;

  int nlp0 = p0TeostVec.size(); // nlp == number of local pts
  int * p0Teost_flag = new int[nlp0];
  int (*spolp0)[3] = new int[nlp0][3];
  // we need a way to lookup the edge index (in the teostVec in this case) from a
  // particular local point. Normally 2 edges are associated with a particular local
  // point, but since we always walk the edges in one direction only, we can populate it
  // with just one...
  int *edolp0 = new int[nlp0];
  for (int ilp = 0; ilp < nlp0; ++ilp) {
    p0Teost_flag[ilp] = 0; // gets set to 1 when used
    FOR_I3 spolp0[ilp][i] = -1; // i-1, i, i+1
    edolp0[ilp] = -1;
  }

  int nlp_check = 0;
  double dx_edge_min = HUGE_VAL;
  for (int ii = 0; ii < p0TeostVec.size(); ++ii) {
    const int ist = int( p0TeostVec[ii] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
    const int i = int( p0TeostVec[ii]>>60 ); assert((i >= 0)&&(i < 3));
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    const double dx = DIST(xsp[isp0],xsp[isp1]);
    dx_edge_min = min(dx_edge_min,dx);
    // isp0...
    int ilp;
    if (sp_flag[isp0] == -1) {
      ilp = nlp_check++;
      assert(nlp_check <= nlp0);
      sp_flag[isp0] = ilp;
      assert(spolp0[ilp][1] == -1);
      spolp0[ilp][1] = isp0;
    }
    else {
      ilp = sp_flag[isp0];
      assert(spolp0[ilp][1] == isp0);
    }
    assert(spolp0[ilp][2] == -1);
    spolp0[ilp][2] = isp1;
    // isp1...
    if (sp_flag[isp1] == -1) {
      ilp = nlp_check++;
      assert(nlp_check <= nlp0);
      sp_flag[isp1] = ilp;
      assert(spolp0[ilp][1] == -1);
      spolp0[ilp][1] = isp1;
    }
    else {
      ilp = sp_flag[isp1];
      assert(spolp0[ilp][1] == isp1);
    }
    assert(spolp0[ilp][0] == -1);
    spolp0[ilp][0] = isp0;
    // provide a way to access the edge from ilp...
    assert(edolp0[ilp] == -1);
    edolp0[ilp] = ii;
  }
  assert(nlp_check == nlp0);

  int nlp1 = p1TeostVec.size();
  int * p1Teost_flag = new int[nlp1];
  int (*spolp1)[3] = new int[nlp1][3];
  int *edolp1 = new int[nlp1];
  for (int ilp = 0; ilp < nlp1; ++ilp) {
    p1Teost_flag[ilp] = 0;
    FOR_I3 spolp1[ilp][i] = -1; // i+1, i, i-1 NOTE: reversed
    edolp1[ilp] = -1;
  }

  // same as above, but use -2 indexing to ensure there are no common nodes...

  nlp_check = 0;
  for (int ii = 0; ii < p1TeostVec.size(); ++ii) {
    const int ist = int( p1TeostVec[ii] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
    const int i = int( p1TeostVec[ii]>>60 ); assert((i >= 0)&&(i < 3));
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    const double dx = DIST(xsp[isp0],xsp[isp1]);
    dx_edge_min = min(dx_edge_min,dx);
    // isp0...
    int ilp;
    if (sp_flag[isp0] == -1) {
      ilp = nlp_check++;
      assert(nlp_check <= nlp1);
      sp_flag[isp0] = -ilp-2;
      assert(spolp1[ilp][1] == -1);
      spolp1[ilp][1] = isp0;
    }
    else {
      ilp = -sp_flag[isp0]-2;
      assert((ilp >= 0)&&(ilp < nlp_check));
      assert(spolp1[ilp][1] == isp0);
    }
    assert(spolp1[ilp][0] == -1);
    spolp1[ilp][0] = isp1;
    // treat this edge index using the opposite sp as for 0 above
    // (opposite direction of looping)...
    assert(edolp1[ilp] == -1);
    edolp1[ilp] = ii;
    // isp1...
    if (sp_flag[isp1] == -1) {
      ilp = nlp_check++;
      assert(nlp_check <= nlp1);
      sp_flag[isp1] = -ilp-2;
      assert(spolp1[ilp][1] == -1);
      spolp1[ilp][1] = isp1;
    }
    else {
      ilp = -sp_flag[isp1]-2;
      assert((ilp >= 0)&&(ilp < nlp_check));
      assert(spolp1[ilp][1] == isp1);
    }
    assert(spolp1[ilp][2] == -1);
    spolp1[ilp][2] = isp0;
  }
  assert(nlp_check == nlp1);

  cout << " > links built successfully: no bow-ties or common points. Smallest edge length: " <<
    dx_edge_min << endl;

  // node merges are going into these types...

  vector<pair<int,int> > node0Node1MergeVec;
  vector<pair<pair<int,double>,int> > edge1Node0MergeVec; // note that for sorting we put edge 1 first

  // now identify "corners" on both sizes...

  const double dp_tol = cos((180.0-crease_angle)*M_PI/180.0);

  vector<pair<int,double> > corner0_vec;
  double dp_max_in = -1.0;
  double dp_min_in = 1.0;
  double dp_min_out = 1.0;
  for (int ilp = 0; ilp < nlp0; ++ilp) {
    const double dx0[3] = DIFF(xsp[spolp0[ilp][1]],xsp[spolp0[ilp][0]]);
    const double dx0_mag = MAG(dx0);
    assert(dx0_mag > 0.0);
    const double dx1[3] = DIFF(xsp[spolp0[ilp][2]],xsp[spolp0[ilp][1]]);
    const double dx1_mag = MAG(dx1);
    assert(dx1_mag > 0.0);
    const double dp = DOT_PRODUCT(dx0,dx1)/(dx0_mag*dx1_mag);
    if (dp <= dp_tol) {
      dp_max_in = max(dp_max_in,dp);
      dp_min_in = min(dp_min_in,dp);
      corner0_vec.push_back(pair<int,double>(ilp,180.0/M_PI*(1.0-acos(dp))));
    }
    else {
      dp_min_out = min(dp_min_out,dp);
    }
  }

  cout << " > side0 has " << corner0_vec.size() << " corners" << endl;
  if (!corner0_vec.empty()) {
    cout << "   > most acute corner: " << 180.0-acos(dp_min_in)*180.0/M_PI << " degrees" << endl;
    cout << "   > closest corner to crease angle: " << 180.0-acos(dp_max_in)*180.0/M_PI << " degrees" << endl;
    cout << "   > closest non-corner to crease angle: " << 180.0-acos(dp_min_out)*180.0/M_PI << " degrees" << endl;
  }

  vector<pair<int,double> > corner1_vec;
  dp_max_in = -1.0;
  dp_min_in = 1.0;
  dp_min_out = 1.0;
  for (int ilp = 0; ilp < nlp1; ++ilp) {
    const double dx0[3] = DIFF(xsp[spolp1[ilp][1]],xsp[spolp1[ilp][0]]);
    const double dx0_mag = MAG(dx0);
    assert(dx0_mag > 0.0);
    const double dx1[3] = DIFF(xsp[spolp1[ilp][2]],xsp[spolp1[ilp][1]]);
    const double dx1_mag = MAG(dx1);
    assert(dx1_mag > 0.0);
    const double dp = DOT_PRODUCT(dx0,dx1)/(dx0_mag*dx1_mag);
    if (dp <= dp_tol) {
      dp_max_in = max(dp_max_in,dp);
      dp_min_in = min(dp_min_in,dp);
      corner1_vec.push_back(pair<int,double>(ilp,180.0/M_PI*(1.0-acos(dp))));
    }
    else {
      dp_min_out = min(dp_min_out,dp);
    }
  }

  cout << " > side1 has " << corner1_vec.size() << " corners" << endl;
  if (!corner1_vec.empty()) {
    cout << "   > most acute corner: " << 180.0-acos(dp_min_in)*180.0/M_PI << " degrees" << endl;
    cout << "   > closest corner to crease angle: " << 180.0-acos(dp_max_in)*180.0/M_PI << " degrees" << endl;
    cout << "   > closest non-corner to crease angle: " << 180.0-acos(dp_min_out)*180.0/M_PI << " degrees" << endl;
  }

  if (corner0_vec.size() != corner1_vec.size()) {
    assert(0);  // TODO: actually, corner counts can be different in this code.
    WUI(WARN,"corner counts differ. Consider changing CREASE_ANGLE if you know these surfaces should be periodic");
    delete[] edolp1;
    delete[] spolp1;
    delete[] p1Teost_flag;
    delete[] edolp0;
    delete[] spolp0;
    delete[] p0Teost_flag;
    delete[] sp_flag;
    return -1;
  }

  // now try and match the corners using BOTH angle and proximity
  // after translation...

  const int nc = corner0_vec.size();
  int * corner_flag = new int[nc];
  for (int ic = 0; ic < nc; ++ic)
    corner_flag[ic] = 0;

  vector<pair<int,double> > corner1_matched_vec;
  double d2_max = 0.0;
  double dangle_max = 0.0;
  for (int i0 = 0; i0 < nc; ++i0) {
    const int ilp0 = corner0_vec[i0].first;
    const int isp0 = spolp0[ilp0][1];
    int i1_closest = -1;
    double d2_closest;
    for (int i1 = 0; i1 < nc; ++i1) {
      const int ilp1 = corner1_vec[i1].first;
      const int isp1 = spolp1[ilp1][1];
      const double d2 = DIST2(xsp[isp0],xsp[isp1]);
      if ((i1_closest == -1)||(d2 < d2_closest)) {
        i1_closest = i1;
        d2_closest = d2;
      }
    }
    assert(i1_closest != -1);
    d2_max = max(d2_max,d2_closest);
    // make sure i1_closest is not already in the corner1_vec
    assert(corner_flag[i1_closest] == 0);
    corner_flag[i1_closest] = 1;
    corner1_matched_vec.push_back(corner1_vec[i1_closest]);
    const double dangle = fabs(corner0_vec[i0].second-corner1_vec[i1_closest].second);
    dangle_max = max(dangle_max,dangle);
  }

  for (int ic = 0; ic < nc; ++ic)
    assert(corner_flag[ic] == 1);
  delete[] corner_flag;

  if (nc > 0) {
    cout << " > all corners matched uniquely: max distance tol: " << sqrt(d2_max) <<
      ", max angle difference: " << dangle_max << " degrees." << endl;
  }

  // now loop through the links. Start by labeling the corners
  // using -1 indexing...


  // ###### link logic from periodic connectivity
  // previously assumes nc == 0, so skips directly to loops logic

  for (int i = 0; i < nc; ++i) {
    const int ilp0 = corner0_vec[i].first;
    const int isp0 = spolp0[ilp0][1];
    assert((isp0 >= 0)&&(isp0 < nsp));
    spolp0[ilp0][1] = -isp0-1;
    const int ilp1 = corner1_matched_vec[i].first;
    const int isp1 = spolp1[ilp1][1];
    assert((isp1 >= 0)&&(isp1 < nsp));
    spolp1[ilp1][1] = -isp1-1;
    // merge these corners once here...
    node0Node1MergeVec.push_back(pair<int,int>(isp0,isp1));
  }

  for (int i = 0; i < nc; ++i) {
    double length0 = 0.0;
    int ilp0 = corner0_vec[i].first;
    int isp0 = -spolp0[ilp0][1]-1;
    assert(isp0 >= 0);
    assert(sp_flag[isp0] == ilp0);
    while (isp0 >= 0) {
      const int isp0_prev = isp0;
      const int isp0_next = spolp0[ilp0][2];
      length0 += DIST(xsp[isp0_next],xsp[isp0_prev]);
      ilp0 = sp_flag[isp0_next];
      isp0 = spolp0[ilp0][1];
      assert(isp0_next == max(isp0,-isp0-1));
    }
    // ilp0 should be at another corner (or the same if a loop)...
    int iend;
    for (iend = 0; iend < nc; ++iend) {
      if (ilp0 == corner0_vec[iend].first)
      break;
    }
    assert(iend < nc);
    cout << " > link0 from corner " << i << " to corner " << iend << " has length: " << length0 << endl;

    // and on the other side...
    double length1 = 0.0;
    int ilp1 = corner1_matched_vec[i].first;
    int isp1 = -spolp1[ilp1][1]-1;
    assert(isp1 >= 0);
    assert(sp_flag[isp1] == -ilp1-2); // recall -2 indexing
    while (isp1 >= 0) {
      const int isp1_prev = isp1;
      const int isp1_next = spolp1[ilp1][2];
      length1 += DIST(xsp[isp1_next],xsp[isp1_prev]);
      ilp1 = -sp_flag[isp1_next]-2;
      isp1 = spolp1[ilp1][1];
      assert(isp1_next == max(isp1,-isp1-1));
    }
    // ilp1 should be at another corner (or the same if a loop)...
    int iend_check;
    for (iend_check = 0; iend_check < nc; ++iend_check) {
      if (ilp1 == corner1_matched_vec[iend_check].first)
        break;
    }
    assert(iend_check < nc);
    cout << " > link1 from corner " << i << " to corner " << iend_check << " has length: " << length1 << endl;
    assert(iend_check == iend);
    cout << " > difference between link lengths: " << length1 - length0 << endl;
    // ensure lengths are within 10%...
    assert(fabs(length1-length0) < 0.1*(length0+length1));

    // -----------------------------------------------------------------------------------------------
    // now advance both cases together determining which edges get cut and which nodes get matched...
    // -----------------------------------------------------------------------------------------------

    double l0 = 0.0;
    double l1 = 0.0;
    ilp0 = corner0_vec[i].first;
    isp0 = -spolp0[ilp0][1]-1;
    assert(isp0 >= 0);
    assert(sp_flag[isp0] == ilp0);
    ilp1 = corner1_matched_vec[i].first;
    isp1 = -spolp1[ilp1][1]-1;
    assert(isp1 >= 0);
    assert(sp_flag[isp1] == -ilp1-2); // recall -2 indexing
    bool b_advance0 = true;
    bool b_advance1 = true;
    int isp0_prev,isp0_next;
    int isp1_prev,isp1_next;
    double l0_prev,l1_prev;
    while (1) {
      // advance isp0...
      if (b_advance0) {
        assert(isp0 >= 0);
        isp0_prev = isp0;
        isp0_next = spolp0[ilp0][2];
        l0_prev = l0;
        l0 += DIST(xsp[isp0_next],xsp[isp0_prev]);
        ilp0 = sp_flag[isp0_next];
        isp0 = spolp0[ilp0][1];
        assert(isp0_next == max(isp0,-isp0-1));
        // flag the edge in the teost flag as touched...
        const int ied0 = edolp0[ilp0];
        p0Teost_flag[ied0] = 1;
      }
      // advance isp1...
      if (b_advance1) {
        assert(isp1 >= 0);
        isp1_prev = isp1;
        isp1_next = spolp1[ilp1][2];
        l1_prev = l1;
        l1 += DIST(xsp[isp1_next],xsp[isp1_prev]);
        ilp1 = -sp_flag[isp1_next]-2;
        isp1 = spolp1[ilp1][1];
        assert(isp1_next == max(isp1,-isp1-1));
        // flag the edge in the teost flag as touched...
        const int ied1 = edolp1[ilp1];
        p1Teost_flag[ied1] = 1;
      }
      // we are done when we have advanced to the next corner...
      if ((isp0 < 0)&&(isp1 < 0))
        break;
      // otherwise...
      if (l0/length0 <= l1/length1) {
        // l0 is not as far along as l1. We should be able to get the
        // next distance for l0 as follows:
        assert(isp0 >= 0);
        const int isp0_next_next = spolp0[ilp0][2];
        const double dx0_next = DIST(xsp[isp0_next_next],xsp[isp0_next]);
        // we can join the nodes together if the new location of the
        // combined node results in a motion that is less than some
        // fraction of the relevant edge lengths less than 50%...
        const double fract0 = (l1*length0-l0*length1)/(2.0*length1*dx0_next);
        const double fract1 = (l1*length0-l0*length1)/(2.0*length0*(l1-l1_prev));
        if ((fract0 < fract_threshold)&&(fract1 < fract_threshold)) {
          // neither of the nodes should be corner. We have already checked
          // isp0 above...
          assert(isp1 >= 0);
          node0Node1MergeVec.push_back(pair<int,int>(isp0,isp1));
          b_advance0 = b_advance1 = true;
        }
        else {
          const int ied1 = edolp1[ilp1];
          {
            const int ist1 = int( p1TeostVec[ied1] & MASK_60BITS ); assert((ist1 >= 0)&&(ist1 < nst));
            const int i = int( p1TeostVec[ied1]>>60 ); assert((i >= 0)&&(i < 3));
            const int isp0_ = spost[ist1][i];
            const int isp1_ = spost[ist1][(i+1)%3];
            assert(isp1_ == isp1_prev);
            assert(isp0_ == isp1_next);
          }
          double w1 = l0/length0-l1_prev/length1;
          double w0 = l1/length1-l0/length0;
          double inv_sum = 1.0/(w0+w1);
          w0 *= inv_sum;
          w1 *= inv_sum;
          assert((w0 > 0.0)&&(w0 < 1.0));
          assert((w1 > 0.0)&&(w1 < 1.0));
          // on this side we still split edges
          edge1Node0MergeVec.push_back(pair<pair<int,double>,int>(pair<int,double>(ied1,w0),isp0));
          b_advance0 = true;
          b_advance1 = false;
        }
      }
      else {
        // l1 is not as far along as l0. We should be able to get the
        // next distance for l1 as follows:
        assert(isp1 >= 0);
        const int isp1_next_next = spolp1[ilp1][2];
        const double dx1_next = DIST(xsp[isp1_next_next],xsp[isp1_next]);
        // we can join the nodes together if the new location of the
        // combined node results in a motion that is less than some
        // fraction of the relevant edge lengths less than 50%...
        const double fract1 = (l0*length1-l1*length0)/(2.0*length0*dx1_next);
        const double fract0 = (l0*length1-l1*length0)/(2.0*length1*(l0-l0_prev));
        if ((fract0 < fract_threshold)&&(fract1 < fract_threshold)) {
          // neither of the nodes should be corner. We have already checked
          // isp0 above...
          assert(isp0 >= 0);
          node0Node1MergeVec.push_back(pair<int,int>(isp0,isp1));
          b_advance0 = b_advance1 = true;
        }
        else {
          const int ied0 = edolp0[ilp0];
          {
            const int ist0 = int( p0TeostVec[ied0] & MASK_60BITS ); assert((ist0 >= 0)&&(ist0 < nst));
            const int i = int( p0TeostVec[ied0]>>60 ); assert((i >= 0)&&(i < 3));
            const int isp0_ = spost[ist0][i];
            const int isp1_ = spost[ist0][(i+1)%3];
            assert(isp0_ == isp0_prev);
            assert(isp1_ == isp0_next);
          }

          double w1 = l1/length1-l0_prev/length0;
          double w0 = l0/length0-l1/length1;
          double inv_sum = 1.0/(w0+w1);
          w0 *= inv_sum;
          w1 *= inv_sum;
          assert((w0 > 0.0)&&(w0 < 1.0));
          assert((w1 > 0.0)&&(w1 < 1.0));
          // this is the part that differs from the periodic version above. Here we
          // cannot allow isp0 to move, so we link isp1 with either isp0_prev or isp0_next,
          // depending on which is closer...
          if (w0 > w1) {
            // we are closer to the start...
            node0Node1MergeVec.push_back(pair<int,int>(isp0_prev,isp1));
          }
          else {
            // we are closer to the end...
            node0Node1MergeVec.push_back(pair<int,int>(isp0_next,isp1));
          }
          b_advance0 = false;
          b_advance1 = true;
        }
      }
      //cout << "isp0: " << isp0 << " l0/length0: " << l0/length0 << " isp1: " << isp1 << " l1/length1: " << l1/length1 << endl;
    }

  }



  // loops...

  // ==================================================================
  // part 2: there may be some edges in p0Teost,p1Teost that had no
  // corners...
  // ==================================================================

  int n0 = 0;
  for (int ii = 0; ii < p0TeostVec.size(); ++ii) {
    if (p0Teost_flag[ii] == 0)
      ++n0;
  }
  int n1 = 0;
  for (int ii = 0; ii < p1TeostVec.size(); ++ii) {
    if (p1Teost_flag[ii] == 0)
      ++n1;
  }

  if ((n0 != 0)||(n1 != 0)) {
    cout << " > some edges were not handled in corner graph. Must be loops: " << n0 << " " << n1 << endl;

    vector<int> lp0Vec;
    vector<double> cc0Vec;
    for (int ii = 0; ii < p0TeostVec.size(); ++ii) {
      if (p0Teost_flag[ii] == 0) {
        p0Teost_flag[ii] = 2;
        // start a new loop from this untouched edge...
        const int iloop0 = lp0Vec.size();
        const int ist = int( p0TeostVec[ii] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
        const int i = int( p0TeostVec[ii]>>60 ); assert((i >= 0)&&(i < 3));
        int isp0 = spost[ist][i];
        int isp1 = spost[ist][(i+1)%3];
        const int ilp0 = sp_flag[isp0];
        assert((ilp0 >= 0)&&(ilp0 < nlp0));
        lp0Vec.push_back(ilp0);
        assert(spolp0[ilp0][1] == isp0);
        assert(spolp0[ilp0][2] == isp1);
        int ilp1 = sp_flag[isp1];
        assert((ilp1 >= 0)&&(ilp1 < nlp0));
        const int ied = edolp0[ilp1];
        assert(ied == ii);
        const double dx = DIST(xsp[isp0],xsp[isp1]);
        assert(iloop0*4 == cc0Vec.size());
        cc0Vec.push_back(dx);
        FOR_I3 cc0Vec.push_back(dx*(xsp[isp0][i]+xsp[isp1][i]));
        while (ilp1 != ilp0) {
          isp0 = isp1;
          assert(spolp0[ilp1][1] == isp1);
          isp1 = spolp0[ilp1][2];
          ilp1 = sp_flag[isp1];
          const double dx = DIST(xsp[isp0],xsp[isp1]);
          cc0Vec[iloop0*4  ] += dx;
          FOR_I3 cc0Vec[iloop0*4+1+i] += dx*(xsp[isp0][i]+xsp[isp1][i]);
          const int ied = edolp0[ilp1];
          assert(p0Teost_flag[ied] == 0);
          p0Teost_flag[ied] = 2;
        }
      }
    }

    cout << " > side0 has " << lp0Vec.size() << " loops without corners" << endl;

    vector<int> lp1Vec;
    vector<double> cc1Vec;
    for (int ii = 0; ii < p1TeostVec.size(); ++ii) {
      if (p1Teost_flag[ii] == 0) {
        p1Teost_flag[ii] = 2;
        // start a new loop from this untouched edge...
        const int iloop1 = lp1Vec.size();
        const int ist = int( p1TeostVec[ii] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
        const int i = int( p1TeostVec[ii]>>60 ); assert((i >= 0)&&(i < 3));
        int isp1 = spost[ist][i]; // note flip here
        int isp0 = spost[ist][(i+1)%3];
        const int ilp0 = -sp_flag[isp0]-2;
        assert((ilp0 >= 0)&&(ilp0 < nlp1));
        lp1Vec.push_back(ilp0);
        assert(spolp1[ilp0][1] == isp0);
        assert(spolp1[ilp0][2] == isp1);
        int ilp1 = -sp_flag[isp1]-2;
        assert((ilp1 >= 0)&&(ilp1 < nlp1));
        const int ied = edolp1[ilp1];
        assert(ied == ii);
        const double dx = DIST(xsp[isp0],xsp[isp1]);
        assert(iloop1*4 == cc1Vec.size());
        cc1Vec.push_back(dx);
        FOR_I3 cc1Vec.push_back(dx*(xsp[isp0][i]+xsp[isp1][i]));
        while (ilp1 != ilp0) {
          isp0 = isp1;
          assert(spolp1[ilp1][1] == isp1);
          isp1 = spolp1[ilp1][2];
          ilp1 = -sp_flag[isp1]-2;
          assert((ilp1 >= 0)&&(ilp1 < nlp1));
          const double dx = DIST(xsp[isp0],xsp[isp1]);
          cc1Vec[iloop1*4  ] += dx;
          FOR_I3 cc1Vec[iloop1*4+1+i] += dx*(xsp[isp0][i]+xsp[isp1][i]);
          const int ied = edolp1[ilp1];
          assert(p1Teost_flag[ied] == 0);
          p1Teost_flag[ied] = 2;
        }
      }
    }

    cout << " > side1 has " << lp1Vec.size() << " loops without corners" << endl;

    // pair loops that are close enough to be zipped...

    vector<pair<int,int> > loopPairVec;
    const double d2_threshold = 1.0E-8;
    int * lp1_flag = new int[lp1Vec.size()];
    for (int ilp1 = 0; ilp1 < lp1Vec.size(); ++ilp1)
      lp1_flag[ilp1] = 0;
    for (int iloop0 = 0; iloop0 < lp0Vec.size(); ++iloop0) {
      // grab our length...
      const double length0 = cc0Vec[iloop0*4];
      double x0[3]; FOR_I3 x0[i] = cc0Vec[iloop0*4+1+i]/(2.0*length0);
      // do n^2 loop on lp1Vec...
      int iloop1_closest = -1;
      double d2_closest;
      for (int iloop1 = 0; iloop1 < lp1Vec.size(); ++iloop1) {
        // grab our length...
        const double length1 = cc1Vec[iloop1*4];
        double x1[3]; FOR_I3 x1[i] = cc1Vec[iloop1*4+1+i]/(2.0*length1);
        // need a single consistent measure of "closeness"...
        const double d2 = DIST2(x1,x0) + (length1-length0)*(length1-length0);
        if ((iloop1_closest == -1)||(d2 < d2_closest)) {
          iloop1_closest = iloop1;
          d2_closest = d2;
        }
      }
      cout << " > iloop0 " << iloop0 << " is closest to iloop1 " << iloop1_closest << " with d2_closest: " << d2_closest << endl;
      if (d2_closest < d2_threshold) {
        cout << "   > this is closer than the threshold " << d2_threshold << " so creating a loop pair" << endl;
        assert(lp1_flag[iloop1_closest] == 0); // if you hit this, d2 must be an issue
        lp1_flag[iloop1_closest] = 1;
        loopPairVec.push_back(pair<int,int>(iloop0,iloop1_closest));
      }
    }
    delete[] lp1_flag;

    cout << "finished pairing loops: loopPairVec.size(): " << loopPairVec.size() << endl;

    // now loop through paired ilp's and connect...

    for (int ii = 0; ii < loopPairVec.size(); ++ii) {

      const int iloop0 = loopPairVec[ii].first;
      const int iloop1 = loopPairVec[ii].second;

      cout << " > working on pair iloop0,iloop1: " << iloop0 << " " << iloop1 << endl;

      int ilp0 = lp0Vec[iloop0];
      int isp0 = spolp0[ilp0][1];
      assert((isp0 >= 0)&&(isp0 < nsp));
      assert(sp_flag[isp0] == ilp0);
      double x0[3];
      FOR_I3 x0[i] = xsp[isp0][i];
      const double length0 = cc0Vec[iloop0*4];
      int isp0_prev = spolp0[ilp0][0];
      double l0 = 0.0;
      double l0_prev = -DIST(xsp[isp0],xsp[isp0_prev]);
      int isp0_next = isp0; // rename to isp0_current

      /*
      {
        FILE * fp = fopen("p0.dat","w");
        double x_t[3];
        FOR_I3 x_t[i] = xsp[isp0][i];
        fprintf(fp,"%18.15le %18.15le %18.15le\n",x_t[0],x_t[1],x_t[2]);
        FOR_I3 x_t[i] = xsp[isp0_prev][i];
        fprintf(fp,"%18.15le %18.15le %18.15le\n",x_t[0],x_t[1],x_t[2]);
        fclose(fp);
      }
      */

      // loop lp1 and find closest...
      int ilp1 = lp1Vec[iloop1];
      int isp1 = spolp1[ilp1][1];
      assert(-sp_flag[isp1]-2 == ilp1);
      int isp1_next = spolp1[ilp1][2];
      int ilp1_next = -sp_flag[isp1_next]-2;
      assert((ilp1_next >= 0)&&(ilp1_next < nlp1));
      int ilp1_closest = ilp1;
      double d2_closest = MiscUtils::getPointToEdgeDist2(x0,xsp[isp1],xsp[isp1_next]);
      assert(d2_closest == d2_closest);
      // recompute length1...
      double length1 = DIST(xsp[isp1],xsp[isp1_next]);
      // loop through ilp1...
      while (ilp1_next != ilp1) {
        isp1 = isp1_next;
        const int ilp1_prev = ilp1_next;
        isp1_next = spolp1[ilp1_next][2];
        ilp1_next = -sp_flag[isp1_next]-2;
        assert((ilp1_next >= 0)&&(ilp1_next < nlp1));
        const double d2 = MiscUtils::getPointToEdgeDist2(x0,xsp[isp1],xsp[isp1_next]);
        assert(d2 == d2);
        if (d2 < d2_closest) {
          ilp1_closest = ilp1_prev;
          d2_closest = d2;
        }
        length1 += DIST(xsp[isp1],xsp[isp1_next]);
      }

      // we now have length0 and length1, and can start the loop matching
      // at the matched point along the edge...

      ilp1 = ilp1_closest;
      isp1 = spolp1[ilp1][1];
      assert((isp1 >= 0)&&(isp1 < nsp));
      isp1_next = spolp1[ilp1][2];
      ilp1_next = -sp_flag[isp1_next]-2;
      double xc1[3];
      MiscUtils::getClosestPointOnEdgeRobust(xc1,x0,xsp[isp1],xsp[isp1_next]);
      // xc1 corresponds to the zero in l1, so...
      double l1 = -DIST(xsp[isp1],xc1);

      /*
      {
        FILE * fp = fopen("p1.dat","w");
        fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp1][0],xsp[isp1][1],xsp[isp1][2]);
        fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp1_next][0],xsp[isp1_next][1],xsp[isp1_next][2]);
        fclose(fp);
      }

      {
        FILE * fp = fopen("xc1.dat","w");
        fprintf(fp,"%18.15le %18.15le %18.15le\n",xc1[0],xc1[1],xc1[2]);
        fclose(fp);
      }

      cout << "length0: " << length0 << " length1: " << length1 << endl;
      cout << "TAKE A LOOK" << endl;
      getchar();
      */

      // unless we find a node-node match right away, we will advance ONLY
      // loop1 (it's l1 is currently -ve to l0's zero value)...
      bool b_advance0 = false;
      bool b_advance1 = true;

      // consider the case when we have a node-node match...
      const double dx1_next = DIST(xsp[isp1_next],xsp[isp1]);
      // we can join the nodes together if the new location of the
      // combined node results in a motion that is less than some
      // fraction of the relevant edge lengths less than 50%...
      const double fract1 = (l0*length1-l1*length0)/(2.0*length0*dx1_next);
      const double fract0 = (l0*length1-l1*length0)/(2.0*length1*(l0-l0_prev));
      if ((fract0 < fract_threshold)&&(fract1 < fract_threshold)) {
        assert(isp0 >= 0);
        node0Node1MergeVec.push_back(pair<int,int>(isp0,isp1));
        b_advance0 = b_advance1 = true;
        // use the corner-marking trick to stop loop...
        assert(spolp0[ilp0][1] == isp0);
        spolp0[ilp0][1] = -isp0-1;
        assert(spolp1[ilp1][1] == isp1);
        spolp1[ilp1][1] = -isp1-1;
      }

      // this code is similar to the loop code between corners...
      assert(isp0 >= 0);
      assert(sp_flag[isp0] == ilp0);
      assert(isp1 >= 0);
      assert(sp_flag[isp1] == -ilp1-2); // recall -2 indexing
      assert(b_advance1);
      int isp1_prev;
      double l1_prev;
      while (1) {

	//cout << "top of while: b_advance0: " << b_advance0 << " b_advance1: " << b_advance1 << endl;

        // advance isp0...
        if (b_advance0) {
          assert(isp0 >= 0);
          isp0_prev = isp0;
          isp0_next = spolp0[ilp0][2];
          l0_prev = l0;
          l0 += DIST(xsp[isp0_next],xsp[isp0_prev]);
          ilp0 = sp_flag[isp0_next];
          isp0 = spolp0[ilp0][1];
          assert(isp0_next == max(isp0,-isp0-1));
          // flag the edge in the teost flag as touched...
          const int ied0 = edolp0[ilp0];
          assert(p0Teost_flag[ied0] == 2);
        }
        // advance isp1...
        if (b_advance1) {
          assert(isp1 >= 0);
          isp1_prev = isp1;
          isp1_next = spolp1[ilp1][2];
          l1_prev = l1;
          l1 += DIST(xsp[isp1_next],xsp[isp1_prev]);
          ilp1 = -sp_flag[isp1_next]-2;
          isp1 = spolp1[ilp1][1];
          assert(isp1_next == max(isp1,-isp1-1));
          // flag the edge in the teost flag as touched...
          const int ied1 = edolp1[ilp1];
          assert(p1Teost_flag[ied1] == 2);
        }
        // we are done when we have advanced to the next corner...
        if ((isp0 < 0)&&(isp1 < 0))
          break;
        // otherwise...
        if (l0/length0 <= l1/length1) {
          // l0 is not as far along as l1. We should be able to get the
          // next distance for l0 as follows:
          assert(isp0 >= 0);
          const int isp0_next_next = spolp0[ilp0][2];
          const double dx0_next = DIST(xsp[isp0_next_next],xsp[isp0_next]);
          // we can join the nodes together if the new location of the
          // combined node results in a motion that is less than some
          // fraction of the relevant edge lengths less than 50%...
          const double fract0 = (l1*length0-l0*length1)/(2.0*length1*dx0_next);
          const double fract1 = (l1*length0-l0*length1)/(2.0*length0*(l1-l1_prev));
          if ((fract0 < fract_threshold)&&(fract1 < fract_threshold)) {
            // neither of the nodes should be corner. We have already checked
            // isp0 above...
            assert(isp1 >= 0);
            node0Node1MergeVec.push_back(pair<int,int>(isp0,isp1));
            b_advance0 = b_advance1 = true;
            // mark the nodes as handled...
            assert(spolp0[ilp0][1] == isp0);
            spolp0[ilp0][1] = -isp0-1;
            assert(spolp1[ilp1][1] == isp1);
            spolp1[ilp1][1] = -isp1-1;
          }
          else {
            const int ied1 = edolp1[ilp1];
            {
              const int ist1 = int( p1TeostVec[ied1] & MASK_60BITS ); assert((ist1 >= 0)&&(ist1 < nst));
              const int i = int( p1TeostVec[ied1]>>60 ); assert((i >= 0)&&(i < 3));
              const int isp0_ = spost[ist1][i];
              const int isp1_ = spost[ist1][(i+1)%3];
              assert(isp1_ == isp1_prev);
              assert(isp0_ == isp1_next);
            }
            double w1 = l0/length0-l1_prev/length1;
            double w0 = l1/length1-l0/length0;
            double inv_sum = 1.0/(w0+w1);
            w0 *= inv_sum;
            w1 *= inv_sum;
            assert((w0 > 0.0)&&(w0 < 1.0));
            assert((w1 > 0.0)&&(w1 < 1.0));
	    edge1Node0MergeVec.push_back(pair<pair<int,double>,int>(pair<int,double>(ied1,w0),isp0));
            b_advance0 = true;
            b_advance1 = false;
            // mark isp0 as handled...
            assert(spolp0[ilp0][1] == isp0);
            spolp0[ilp0][1] = -isp0-1;
          }
        }
        else {
          // l1 is not as far along as l0. We should be able to get the
          // next distance for l1 as follows:
          assert(isp1 >= 0);
          const int isp1_next_next = spolp1[ilp1][2];
          const double dx1_next = DIST(xsp[isp1_next_next],xsp[isp1_next]);
          // we can join the nodes together if the new location of the
          // combined node results in a motion that is less than some
          // fraction of the relevant edge lengths less than 50%...
          const double fract1 = (l0*length1-l1*length0)/(2.0*length0*dx1_next);
          const double fract0 = (l0*length1-l1*length0)/(2.0*length1*(l0-l0_prev));
          if ((fract0 < fract_threshold)&&(fract1 < fract_threshold)) {
            // neither of the nodes should be corner. We have already checked
            // isp0 above...
            assert(isp0 >= 0);
            node0Node1MergeVec.push_back(pair<int,int>(isp0,isp1));
            b_advance0 = b_advance1 = true;
            // mark the nodes as handled...
            assert(spolp0[ilp0][1] == isp0);
            spolp0[ilp0][1] = -isp0-1;
            assert(spolp1[ilp1][1] == isp1);
            spolp1[ilp1][1] = -isp1-1;
          }
          else {
            const int ied0 = edolp0[ilp0];
            {
              const int ist0 = int( p0TeostVec[ied0] & MASK_60BITS ); assert((ist0 >= 0)&&(ist0 < nst));
              const int i = int( p0TeostVec[ied0]>>60 ); assert((i >= 0)&&(i < 3));
              const int isp0_ = spost[ist0][i];
              const int isp1_ = spost[ist0][(i+1)%3];
              assert(isp0_ == isp0_prev);
              assert(isp1_ == isp0_next);
            }
	    double w1 = l1/length1-l0_prev/length0;
            double w0 = l0/length0-l1/length1;
            double inv_sum = 1.0/(w0+w1);
            w0 *= inv_sum;
            w1 *= inv_sum;
            assert((w0 > 0.0)&&(w0 < 1.0));
            assert((w1 > 0.0)&&(w1 < 1.0));
	    // this is the part that differs from the periodic version above. Here we
	    // cannot allow isp0 to move, so we link isp1 with either isp0_prev or isp0_next,
	    // depending on which is closer...
	    if (w0 > w1) {
	      // we are closer to the start...
	      node0Node1MergeVec.push_back(pair<int,int>(isp0_prev,isp1));
	    }
	    else {
	      // we are closer to the end...
	      node0Node1MergeVec.push_back(pair<int,int>(isp0_next,isp1));
	    }
	    b_advance0 = false;
	    b_advance1 = true;
            // mark isp1 as handled...
            assert(spolp1[ilp1][1] == isp1);
            spolp1[ilp1][1] = -isp1-1;
          }
        }

        //cout << "isp0: " << isp0 << " l0/length0: " << l0/length0 << " isp1: " << isp1 << " l1/length1: " << l1/length1 << endl;
        //getchar();

      }

      //cout << "!!!!!!!!!!!!!!!! COMPLETED loop: " << iloop << " !!!!!!!!!!!!!!!!!!!!" << endl;

    }

    //cout << "DONE FIRST" << endl;
    //getchar();

  }

  delete[] edolp1;
  delete[] spolp1;
  delete[] p1Teost_flag;
  delete[] edolp0;
  delete[] spolp0;
  delete[] p0Teost_flag;

  /*
    {
    FILE * fp = fopen("p0.dat","w");
    for (int ii = 0; ii < node0Node1MergeVec.size(); ++ii) {
    const int isp0 = node0Node1MergeVec[ii].first;
    fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp0][0],xsp[isp0][1],xsp[isp0][2]);
    }
    fclose(fp);
    fp = fopen("p1.dat","w");
    for (int ii = 0; ii < node0Node1MergeVec.size(); ++ii) {
    const int isp1 = node0Node1MergeVec[ii].second;
    fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp1][0],xsp[isp1][1],xsp[isp1][2]);
    }
    fclose(fp);
    }
  */

  // use sp_flag to store the new node indices...
  for (int isp = 0; isp < nsp; ++isp)
    sp_flag[isp] = isp;

  // node-node merges...
  d2_max = 0.0;
  for (int ii = 0; ii < node0Node1MergeVec.size(); ++ii) {
    const int isp0 = node0Node1MergeVec[ii].first;
    const int isp1 = node0Node1MergeVec[ii].second;
    const double d2 = DIST2(xsp[isp0],xsp[isp1]);
    d2_max = max(d2_max,d2);
    sp_flag[isp1] = isp0;
  }
  cout << " > node-node merge dist max: " << sqrt(d2_max) << endl;

  // then node-edge...
  st_flag.setLength(nst);
  st_flag.setAll(0);

  // sort in order of weight...
  sort(edge1Node0MergeVec.begin(),edge1Node0MergeVec.end());

  int ie_current = -1;
  int nst_new = 0;
  for (int ii = 0, limit = edge1Node0MergeVec.size(); ii < limit; ++ii) {
    if (edge1Node0MergeVec[ii].first.first != ie_current) {
      ie_current = edge1Node0MergeVec[ii].first.first;
      const int ist = int( p1TeostVec[ie_current] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
      st_flag[ist] += 1;
      const int i = int( p1TeostVec[ie_current]>>60 ); assert((i >= 0)&&(i < 3));
      assert(isEdgeOpen(ist,i));
    }
    nst_new += 1;
  }

  // check if there are any tris that have 2 (or even 3?!) edges participating
  // in the edges requiring splitting. These tris need to be split...

  int nst_split_tri_new = 0;
  int nsp_split_tri_new = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] > 1) {
      st_flag[ist] = nsp_split_tri_new++; // used to index new tris and new nodes...
      nst_split_tri_new += 2; // two new tris
    }
    else {
      st_flag[ist] = -1;
    }
  }

  // we now know exactly how big everything has to be, so grow memory...

  const int nst_old = nst;
  nst = nst_old + nst_split_tri_new + nst_new;
  const int nsp_old = nsp;
  nsp = nsp_old + nsp_split_tri_new;

  // resize st and sp stuff...
  // note we will be clearing teost later...
  if (nst > nst_old) {
    int (*spost_old)[3] = spost;
    spost = new int[nst][3];
    for (int ist = 0; ist < nst_old; ++ist) FOR_I3 spost[ist][i] = spost_old[ist][i];
    delete[] spost_old;
    int *znost_old = znost;
    znost = new int[nst];
    for (int ist = 0; ist < nst_old; ++ist) znost[ist] = znost_old[ist];
    delete[] znost_old;
    // also the subzone indexing...
    szost.resize(nst);
  }

  if (nsp > nsp_old) {
    double (*xsp_old)[3] = xsp;
    xsp = new double[nsp][3];
    for (int isp = 0; isp < nsp_old; ++isp) FOR_I3 xsp[isp][i] = xsp_old[isp][i];
    delete[] xsp_old;
  }

  // no pbi: clear periodicity as part of this routine if you ever get this...
  assert(pbi == NULL);

  // build new tris and nodes...

  int ist_new = nst_old;
  int isp_new = nsp_old;

  // start by splitting the tris with double edges...

  for (int ist = 0; ist < nst_old; ++ist) {
    if (st_flag[ist] >= 0) {
      assert(ist_new == nst_old+2*st_flag[ist]);
      assert(isp_new == nsp_old+st_flag[ist]);
      const int isp0 = spost[ist][0];
      const int isp1 = spost[ist][1];
      const int isp2 = spost[ist][2];
      FOR_I3 xsp[isp_new][i] = (xsp[isp0][i] + xsp[isp1][i] + xsp[isp2][i])/3.0;
      // split the original tri into three tris such that the outer edges of
      // these tris are consistent with the original edge ordering...
      spost[ist][2] = isp_new;
      // the first new tri...
      spost[ist_new][0] = isp_new;
      spost[ist_new][1] = isp1;
      spost[ist_new][2] = isp2;
      znost[ist_new] = znost[ist];
      szost[ist_new] = szost[ist];
      // the second new tri...
      spost[ist_new+1][0] = isp0;
      spost[ist_new+1][1] = isp_new;
      spost[ist_new+1][2] = isp2;
      znost[ist_new+1] = znost[ist];
      szost[ist_new+1] = szost[ist];
      // increment the new counts...
      isp_new += 1;
      ist_new += 2;
    }
  }

  assert(isp_new == nsp_old + nsp_split_tri_new);
  assert(ist_new == nst_old + nst_split_tri_new);

  // also a new node and a new tri is required for each edge...

  ie_current = -1;
  double x0[3],x1[3]; // for checking
  for (int ii = 0, limit = edge1Node0MergeVec.size(); ii < limit; ++ii) {
    const int ie = edge1Node0MergeVec[ii].first.first;
    const int ist_orig = int( p1TeostVec[ie] & MASK_60BITS ); assert((ist_orig >= 0)&&(ist_orig < nst));
    const int i = int( p1TeostVec[ie]>>60 ); assert((i >= 0)&&(i < 3));
    assert(isEdgeOpen(ist_orig,i));
    int isp_new = edge1Node0MergeVec[ii].second;
    assert(sp_flag[isp_new] == isp_new); // should be an unchanging node on the cht side
    // if ist_orig is a tri that has been duplicated, then adjust ist to the duplicate. Note
    // that because of the way we oriented the new tris, there is no need to change the
    // edge access pattern. Just modify the ist to the appropriate new tri...
    int ist = ist_orig;
    if ((st_flag[ist] >= 0)&&(i > 0)) {
      ist = nst_old+2*st_flag[ist]+i-1; // i = 1,2
    }
    if (ie != ie_current) {
      // on the first time we touch an edge, grap the end coordinates. These are used
      // to compute all the new nodes...
      ie_current = ie;
      const int isp0 = spost[ist][i];
      FOR_J3 x0[j] = xsp[isp0][j];
      const int isp1 = spost[ist][(i+1)%3];
      FOR_J3 x1[j] = xsp[isp1][j];
    }
    // the new node is at position "first.second"...
    /*
    {
      double x_check[3];
      FOR_J3 x_check[j] = (1.0-edge1Node0MergeVec[ii].first.second)*x0[j] + edge1Node0MergeVec[ii].first.second*x1[j];
      // check proximity to the transformed point...
      const double dist = DIST(x_check,xsp[isp_new]);
      // try opposite orientation...
      FOR_J3 x_check[j] = (1.0-edge1Node0MergeVec[ii].first.second)*x1[j] + edge1Node0MergeVec[ii].first.second*x0[j];
      // check proximity to the transformed point...
      const double dist_opp = DIST(x_check,xsp[isp_new]);
      cout << "NODE dist during edge split: " << dist << " dist_opp: " << dist_opp << endl;
      getchar();
    }
    */
    // the new tri goes in where the old tri was...
    spost[ist_new][i]       = spost[ist][i];
    spost[ist_new][(i+1)%3] = isp_new;
    spost[ist_new][(i+2)%3] = spost[ist][(i+2)%3]; // same internal node
    // the new node becomes isp0 of the original tri...
    // for the case of multiple splits, the oriinal tri ust keeps getting
    // pushed to the end...
    spost[ist][i] = isp_new; // do this AFTER above
    // new zone is the same zone...
    znost[ist_new] = znost[ist];
    // new subzone is the same...
    szost[ist_new] = szost[ist];
    ++ist_new;
  }
  assert(isp_new == nsp);
  assert(ist_new == nst);

  // teost is now broken...
  clearTeost();

  // now loop through tris updating nodes and removing slivers (2 or 3 nodes the same)...
  ist_new = 0;
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      const int isp = spost[ist][i];
      spost[ist][i] = max(sp_flag[isp],-sp_flag[isp]-1);
      assert(spost[ist][i] == max(sp_flag[spost[ist][i]],-sp_flag[spost[ist][i]]-1)); // should never be multiple links
    }
    if ((spost[ist][0] != spost[ist][1])&&(spost[ist][0] != spost[ist][2])&&(spost[ist][1] != spost[ist][2])) {
      if (ist_new != ist) {
	FOR_I3 spost[ist_new][i] = spost[ist][i];
	znost[ist_new] = znost[ist];
	szost[ist_new] = szost[ist];
      }
      FOR_I3 {
	// flip sp_flag to -1 indexing in touched nodes...
	if (sp_flag[spost[ist_new][i]] >= 0) {
	  assert(sp_flag[spost[ist_new][i]] == spost[ist_new][i]);
	  sp_flag[spost[ist_new][i]] = -spost[ist_new][i]-1;
	}
      }
      ++ist_new;
    }
  }
  nst = ist_new;

  // now nodes being kept are flagged with sp_flag < 0...
  isp_new = 0;
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] < 0) {
      if (isp_new != isp) {
	FOR_I3 xsp[isp_new][i] = xsp[isp][i];
      }
      sp_flag[isp] = isp_new++;
    }
    else {
      sp_flag[isp] = -1;
    }
  }
  nsp = isp_new;

  // finally one last loop through the tris to update the point index...
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      const int isp_new = sp_flag[spost[ist][i]];
      assert((isp_new >= 0)&&(isp_new < nsp));
      spost[ist][i] = isp_new;
    }
  }

  delete[] sp_flag;

  // =======END======= new delete[] SimpleSurface::zipOpenEdgesChtFlaggedZones =======END========

}
