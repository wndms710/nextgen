#include "SimpleSurface.hpp"
#include "WebUI.hpp"

void SimpleSurface::translateSurface(const double dx[3]) {
  sp_flag.resize(nsp);
  sp_flag.setAll(1);

  translateFlaggedNodes(dx);
}

void SimpleSurface::translateFlaggedTris(const double dx[3]) {
  flagNodesOfFlaggedTris();
  translateFlaggedNodes(dx);
}

void SimpleSurface::translateFlaggedNodes(const double dx[3]) {

  if (sp_flag.count() == 0) {
    COUT2(" > no surface nodes were selected for translation; skipping");
    return;
  }

  COUT2(" > dx: " << dx[0] << " " << dx[1] << " " << dx[2]);

  // surface will change, recompute geoms
  clearDynamicEdgeGroups();
  clearOpenEdgeGroups();
  clearSubzoneData();
  clearZoneData();
  b_centroid = false;
  b_rmax = false;
  b_boundingBox = false;
  clearNonManifoldData();

  int nsp_move = 0;
  FOR_ISP {
    if (sp_flag[isp]) {
      FOR_J3 xsp[isp][j] += dx[j];
      ++nsp_move;
    }
  }
  COUT2(" > " << nsp_move << " surface nodes were translated");

  // if we had periodic data, we need to reset it here
  if (pbi) {
    clearPeriodicity();
    WUI(WARN,"Good news: we've moved some surfaces around...Bad news: stuff is different now, so you need to reset periodicity.");
  }
}

void SimpleSurface::mirrorFlaggedTris(const double xc[3],const double np[3]) {
  flagNodesOfFlaggedTris();
  mirrorFlaggedNodes(xc,np);
}

void SimpleSurface::mirrorFlaggedNodes(const double xc[3],const double np[3]) {

  if (sp_flag.count() == 0) {
    COUT2(" > no surface nodes were selected for mirroring; skipping");
    return;
  }

  COUT2(" > mirror plane (pt,normal): " << COUT_VEC(xc) << " " << COUT_VEC(np));

  // surface will change, recompute geoms
  clearDynamicEdgeGroups();
  clearOpenEdgeGroups();
  clearSubzoneData();
  clearZoneData();
  b_centroid = false;
  b_rmax = false;
  b_boundingBox = false;
  clearNonManifoldData();

  int nsp_move = 0;
  FOR_ISP {
    if (sp_flag[isp]) {
      double x_plane[3];
      int ierr = MiscUtils::computeLinePlaneIntersection(x_plane,xsp[isp],np,xc,np);
      if (ierr == 0) {
        FOR_J3 xsp[isp][j] = x_plane[j] + (x_plane[j]-xsp[isp][j]);
      }
      ++nsp_move;
    }
  }
  COUT2(" > " << nsp_move << " surface nodes were translated");

  // if we had periodic data, we need to reset it here
  if (pbi) {
    clearPeriodicity();
    WUI(WARN,"Good news: we've moved some surfaces around...Bad news: stuff is different now, so you need to reset periodicity.");
  }
}

// st_flag needs to be sized and set prior to calling this
void SimpleSurface::translateFlaggedTrisNormal(const double dn, bool b_flagged_only) {

  if (st_flag.count() == 0) {
    COUT2(" > no surface tris were selected for translation; skipping");
    return;
  }

  COUT2(" > dn: " << dn);

  /*
  if ( pbi) {
    clearPeriodicity();
  }
  */

  // going to need node normal to do the perturbation -- building it now...

  double (*sp_normal)[3] = new double[nsp][3];

  for (int isp=0; isp<nsp; ++isp) {
    FOR_I3 sp_normal[isp][i] = 0.0;
  }

  for (int ist=0; ist<nst; ++ist) {

    if ( (b_flagged_only && (st_flag[ist] == 1)) || !b_flagged_only) {
      const double this_n[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
      FOR_I3 FOR_J3 sp_normal[spost[ist][i]][j] += this_n[j];
    }
  }

  for (int isp=0; isp<nsp; ++isp) {
    const double mag = MAG(sp_normal[isp]);
    if ( mag > 0.0) {
      FOR_I3 sp_normal[isp][i] /= mag;
    }
    else {
      FOR_I3 sp_normal[isp][i] = 0.0;
    }
  }


  // surface will change, recompute geoms
  clearSubzoneData();
  clearZoneData();
  b_centroid = false;
  b_rmax = false;
  b_boundingBox = false;
  clearNonManifoldData();

  // keep track of which nodes will  be moved
  sp_flag.setLength(nsp);
  sp_flag.setAll(0);  // if stays 0 then don't modify

  double sp_normal_avg[3] = { 0.0, 0.0, 0.0 };
  int nsp_move = 0;
  FOR_IST {
    if (st_flag[ist]) {
      // move this tri
      FOR_I3 {
        const int isp = spost[ist][i];
        if (!sp_flag[isp]) {
          // only update if first time touched
          sp_flag[isp] = 1;
          //cout << " > got one: " <<  isp << "   " << COUT_VEC(sp_normal[isp]) << "    " << COUT_VEC(xsp[isp]) << endl;
          FOR_J3 xsp[isp][j] += dn*sp_normal[isp][j];
          FOR_J3 sp_normal_avg[j] += sp_normal[isp][j];
          ++nsp_move;
        }
      }
    }
  }

  // report the average normal direction...
  const double sp_normal_avg_mag = MAG(sp_normal_avg);
  if (sp_normal_avg_mag > 0.0) {
    FOR_I3 sp_normal_avg[i] /= sp_normal_avg_mag;
    COUT2(" > average normal direction: " << COUT_VEC(sp_normal_avg));
  }

  delete[] sp_normal;
  //TODO output of nodes moved?
}

void SimpleSurface::rotateFlaggedTris(const double _axis[3],const double point[3],const double angle_deg) {
  flagNodesOfFlaggedTris();
  rotateFlaggedNodes(_axis,point,angle_deg);
}

void SimpleSurface::rotateFlaggedNodes(const double _axis[3],const double point[3],const double angle_deg) {

  if (sp_flag.count() == 0) {
    COUT2(" > no surface nodes were selected for rotation; skipping");
    return;
  }

  COUT2(" > rotation axis (pt,normal) and angle: " << COUT_VEC(point) << " " << COUT_VEC(_axis) << " " << angle_deg << " degrees");

  // surface will change, recompute geoms
  clearDynamicEdgeGroups();
  clearOpenEdgeGroups();
  clearSubzoneData();
  clearZoneData();
  b_centroid = false;
  b_rmax = false;
  b_boundingBox = false;
  clearNonManifoldData();

  double axis[3] = {_axis[0],_axis[1],_axis[2]};
  NORMALIZE(axis);
  const double angle = angle_deg/180.0*M_PI;  // radians

  const double cp_mat[3][3] = {
    {0.0,-axis[2],axis[1]},
    {axis[2],0.0,-axis[0]},
    {-axis[1],axis[0],0.0}
  };

  // create rotation matrix
  // R = cos0*I + sin0*[axis]_x + (1-cos0)*(u tensor-product u)
  double R[3][3];
  FOR_I3 {
    R[i][i] = cos(angle);
    FOR_J3 {
      if (i != j) R[i][j] = sin(angle)*cp_mat[i][j];
      R[i][j] += (1-cos(angle))*axis[i]*axis[j];
    }
  }


  int nsp_move = 0;
  FOR_ISP {
    if (sp_flag[isp]) {
      double _xsp[3];
      FOR_J3 _xsp[j] = xsp[isp][j] - point[j];
      const double temp[3] = MATRIX_PRODUCT(R,_xsp);
      FOR_J3 xsp[isp][j] = temp[j] + point[j];
      ++nsp_move;
    }
  }
  COUT2(" > " << nsp_move << " surface nodes were rotated");

  // if we had periodic data, we need to reset it here
  if (pbi) {
    clearPeriodicity();
    WUI(WARN,"Good news: we've moved some surfaces around...Bad news: stuff is different now, so you need to reset periodicity.");
  }
}

void SimpleSurface::scaleFlaggedTris(const double sx[3],const double x0[3],bool b_norm) {
  flagNodesOfFlaggedTris();
  scaleFlaggedNodes(sx,x0,b_norm);
}

void SimpleSurface::scaleFlaggedTrisRadially(const double factor,const double x0[3],const double axis[3],bool b_norm) {
  flagNodesOfFlaggedTris();
  scaleFlaggedNodesRadially(factor,x0,axis,b_norm);
}

void SimpleSurface::scaleFlaggedNodes(const double sx[3],const double x0[3],const bool b_norm) {

  if (sp_flag.count() == 0) {
    COUT2(" > no surface nodes were selected for scaling; skipping");
    return;
  }

  COUT2(" > scaling about point: " << COUT_VEC(x0) << " with factors (x,y,z): " << COUT_VEC(sx));

  // surface will change, recompute geoms
  clearDynamicEdgeGroups();
  clearOpenEdgeGroups();
  clearSubzoneData();
  clearZoneData();
  b_centroid = false;
  b_rmax = false;
  b_boundingBox = false;
  clearNonManifoldData();

  int nsp_move = 0;
  FOR_ISP {
    if (sp_flag[isp]) {
      // FOR_J3 xsp[isp][j] = x0[j] + sx[j]*(xsp[isp][j] - x0[j]);
      const double dx[3] = DIFF(xsp[isp],x0);
      double vmag = 0.0;
      FOR_I3 vmag += dx[i]*dx[i];
      vmag = sqrt(vmag);

      FOR_I3 {
        if (sx[i] == 0.0) continue;  // zero is treated as no-op instead of multiple by zero...

        if (b_norm) xsp[isp][i] += (sx[i]-1.0)*dx[i];
        else xsp[isp][i] += sx[i]*dx[i]/vmag;
      }
      ++nsp_move;
    }
  }
  COUT2(" > " << nsp_move << " surface nodes were scaled");

  // if we had periodic data, we need to reset it here
  if (pbi && (sp_flag.count() == nsp)) {
    const int npt = PeriodicData::periodicTransformVec.size();
    for (int ipt = 0; ipt < npt; ++ipt) PeriodicData::periodicTransformVec[ipt].scalePeriodicTransform(sx);

    // this sbin has periodic information...
    cout << " > sbin has " << npt << " periodic transforms..." << endl;
    for (int ipt = 0; ipt < npt; ++ipt) {
      cout << "    > " << ipt << ": "; PeriodicData::periodicTransformVec[ipt].dump();
    }

    WUI(WARN,"The CART periodicities have also been scaled.");
  }
  else if (pbi) {
    clearPeriodicity();
    WUI(WARN,"Good news: we've moved some surfaces around...Bad news: stuff is different now, so you need to reset periodicity.");
  }
}

void SimpleSurface::scaleFlaggedNodesRadially(const double factor,const double x0[3],const double axis[3],const bool b_norm) {

  if (sp_flag.count() == 0) {
    COUT2(" > no surface nodes were selected for scaling; skipping");
    return;
  }

  COUT2(" > scaling about axis (x0,normal): " << COUT_VEC(x0) << ", " << COUT_VEC(axis));
  if (b_norm) {
    COUT2("    > by a factor of: " << factor);
  }
  else {
    COUT2("    > by a distance of: " << factor);
  }

  // surface will change, recompute geoms
  clearDynamicEdgeGroups();
  clearOpenEdgeGroups();
  clearSubzoneData();
  clearZoneData();
  b_centroid = false;
  b_rmax = false;
  b_boundingBox = false;
  clearNonManifoldData();

  double e0[3] = {axis[0],axis[1],axis[2]};
  NORMALIZE(e0);
  double e1[3] = {0.0,0.0,0.0};
  double e2[3] = {0.0,0.0,0.0};
  MiscUtils::getOrthogonalVectors(e1,e2,x0,axis,e1);

  int nsp_move = 0;
  FOR_ISP {
    if (sp_flag[isp]) {
      const double dx[3] = DIFF(xsp[isp],x0);
      double imag = DOT_PRODUCT(dx,e1);
      double jmag = DOT_PRODUCT(dx,e2);
      const double kmag = DOT_PRODUCT(dx,axis);  // constant

      // compute new distances to get desired r
      if (b_norm) {
        imag *= factor;
        jmag *= factor;
      }
      else {
        const double rmag = sqrt(imag*imag + jmag*jmag);
        const double fac = (rmag-factor)/rmag;
        imag *= fac;
        jmag *= fac;
      }
      // set new point location
      FOR_I3 xsp[isp][i] = x0[i] + kmag*axis[i] + imag*e1[i] + jmag*e2[i];

      ++nsp_move;
    }
  }
  COUT2(" > " << nsp_move << " surface nodes were radially scaled");

  // if we had periodic data, we need to reset it here
  if (pbi) {
    clearPeriodicity();
    WUI(WARN,"Good news: we've moved some surfaces around...Bad news: stuff is different now, so you need to reset periodicity.");
  }
}

void SimpleSurface::mirrorSelectedSubzones(const vector<int>& subzonesVec,const double xp[3],const double np[3]) {
  flagTrisFromSubzoneVec(subzonesVec);  // flagged tris in st_flag
  mirrorFlaggedTris(xp,np);
}

void SimpleSurface::translateSelectedSubzones(const vector<int>& subzonesVec,const double dx[3]) {
  flagTrisFromSubzoneVec(subzonesVec);  // flagged tris in st_flag
  translateFlaggedTris(dx);
}

void SimpleSurface::translateOpenEdges(const vector<int>& edge_indices, const double dx[3],const int n_duplicates) {
  ensureOpenEdgeGroups();

  if (!n_open_edge_groups) {
    WUI(WARN,"no open edge loops to process; skipping");
    return;
  }

  sp_flag.resize(nsp);

  IntFlag oeg_flag(n_open_edge_groups);

  if (edge_indices.empty()) {
    // if not specified then assume all open edge nodes should be scaled
    oeg_flag.setAll(1);
  }
  else {
    oeg_flag.setAll(0);
    for (vector<int>::const_iterator it=edge_indices.begin(); it!=edge_indices.end(); ++it) {
      oeg_flag[*it] = 1;
    }
  }

  // populate gr_to_edge structure
  vector<pair<int,uint> > group_to_edges;
  for (map<uint,int>::const_iterator it=oe_to_group.begin(); it!= oe_to_group.end(); ++it) {
    if (oeg_flag[it->second]) {
      group_to_edges.push_back(pair<int,uint> (it->second,it->first));
    }
  }
  COUT2(" > translating " << group_to_edges.size() << " edges in total");

  // sort by group index; should only contain the flagged edges
  sort(group_to_edges.begin(),group_to_edges.end());

  if (n_duplicates) {
    // first populate duplicates in nsp, then use scale to modify appropriate cloned points
    const int ss_nsp0 = nsp;
    nsp += n_duplicates*group_to_edges.size();  // because open edges, should only be applying this to loops anyways, so for now assume 1 node per edge
    growNspData(nsp,ss_nsp0);
    sp_flag.resize(nsp);

    double dx_dup[3];
    FOR_I3 dx_dup[i] = dx[i]/double(n_duplicates);

    vector<NewTri> newTris;
    int new_isp = ss_nsp0;
    for (vector<pair<int,uint> >::const_iterator it=group_to_edges.begin(); it!=group_to_edges.end(); ) {

      const int group = it->first;
      vector<uint> heVec;
      while ((it->first == group)&&(it!=group_to_edges.end())) {
        heVec.push_back(it->second);
        ++it;
      }
      orderHalfEdgesInVec(heVec);  // uses sp_flag, so make sure no clobbering

      vector<int> isp0Vec;
      vector<int> isp1Vec;

      // peel first loop b/c isp0 indices not ordered; others are in xsp
      sp_flag.setAll(0);
      for (vector<uint>::const_iterator hit=heVec.begin(); hit!=heVec.end(); ++hit) {
        uint ist; uchar i;
        unpackEdge(ist,i,*hit);
        const int isp0 = spost[ist][i];
        isp0Vec.push_back(isp0);
        isp1Vec.insert(isp1Vec.begin(),new_isp);  // needs to be in reverse order...
        FOR_I3 xsp[new_isp][i] = xsp[isp0][i];
        sp_flag[new_isp] = 1;
        ++new_isp;
      }
      // sp_flag should have flagged nodes just on this group
      // apply translation for the selected nodes
      translateFlaggedNodes(dx_dup);
      facetGap(newTris,isp0Vec,isp1Vec,true);  // HACK: revisit loop boolean

      if (n_duplicates > 1) {
        // additional loops if present
        const int n_loop_nodes = isp0Vec.size();
        for (int idup=1; idup<n_duplicates; ++idup) {
          sp_flag.setAll(0);
          for (int ino=0; ino<n_loop_nodes; ++ino) {
            isp0Vec[n_loop_nodes-ino-1] = isp1Vec[ino]; // needs to be reverse order
            isp1Vec[ino] = new_isp;
            FOR_I3 xsp[new_isp][i] = xsp[isp0Vec[n_loop_nodes-ino-1]][i];
            sp_flag[new_isp] = 1;
            ++new_isp;
          }
          translateFlaggedNodes(dx_dup);
          facetGap(newTris,isp0Vec,isp1Vec,true);  // HACK: revisit loop boolean
        }
      }
    }

    if (!newTris.empty()) {
      // all new tris get put into a new zone, so we don't have to keep track of which tri they were originally attached to (and for easy deletion if desired)
      const int ss_nst0 = nst;
      nst += newTris.size();
      growNstData(nst,ss_nst0);
      const int nzn0 = zoneVec.size();
      const int nsz0 = nsz++;
      zoneVec.push_back(SurfaceZone("translated_edges"));
      szozn_i.resize(zoneVec.size()+1);
      szozn_i[zoneVec.size()] = nsz;

      for (vector<NewTri>::const_iterator it=newTris.begin(); it!=newTris.end(); ++it) {
        const int ist = ss_nst0 + (it-newTris.begin());
        FOR_I3 spost[ist][i] = it->spost[i];
        znost[ist] = nzn0;
        szost[ist] = nsz0;
      }
      clearTeost();
    }
  }
  else {
    for (vector<pair<int,uint> >::const_iterator it=group_to_edges.begin(); it!=group_to_edges.end(); ) {
      const int group = it->first;
      sp_flag.setAll(0);
      while ((it->first == group)&&(it!=group_to_edges.end())) {
        uint ist; uchar i;
        unpackEdge(ist,i,it->second);
        const int isp0 = spost[ist][i];
        const int isp1 = spost[ist][(i+1)%3];
        sp_flag[isp0]  = sp_flag[isp1] = 1;
        ++it;
      }

      // sp_flag should have flagged nodes just on this group
      // apply translation for the selected nodes
      translateFlaggedNodes(dx);
    }
  }
}

void SimpleSurface::translateSelectedSubzonesNormal(const vector<int>& subzonesVec,const double dn) {

  flagTrisFromSubzoneVec(subzonesVec);  // flagged tris in st_flag

  // restrict the normal curvature to the flagged subzones...
  translateFlaggedTrisNormal(dn,true);

  // if we had periodic data, we need to reset it here
  if (pbi) {
    clearPeriodicity();
    WUI(WARN,"Bad news, friend. You need to reset periodicity.");
  }
}


void SimpleSurface::rotateSelectedSubzones(const vector<int>& subzonesVec,const double axis[3], const double point[3],const double angle_deg) {
  flagTrisFromSubzoneVec(subzonesVec);  // flagged tris in st_flag
  rotateFlaggedTris(axis,point,angle_deg);
}

void SimpleSurface::rotateSurface(const double axis[3],const double point[3], const double angle_deg) {
  sp_flag.setLength(nsp);
  sp_flag.setAll(1);
  rotateFlaggedNodes(axis,point,angle_deg);
}

void SimpleSurface::rotateOpenEdges(const vector<int>& edge_indices, const double axis[3],const double point[3],const double angle_deg,const int n_duplicates) {
  ensureOpenEdgeGroups();

  if (!n_open_edge_groups) {
    WUI(WARN,"no open edge loops to process; skipping");
    return;
  }

  sp_flag.resize(nsp);

  IntFlag oeg_flag(n_open_edge_groups);

  if (edge_indices.empty()) {
    // if not specified then assume all open edge nodes should be scaled
    oeg_flag.setAll(1);
  }
  else {
    oeg_flag.setAll(0);
    for (vector<int>::const_iterator it=edge_indices.begin(); it!=edge_indices.end(); ++it) {
      oeg_flag[*it] = 1;
    }
  }

  // populate gr_to_edge structure
  vector<pair<int,uint> > group_to_edges;
  for (map<uint,int>::const_iterator it=oe_to_group.begin(); it!= oe_to_group.end(); ++it) {
    if (oeg_flag[it->second]) {
      group_to_edges.push_back(pair<int,uint> (it->second,it->first));
    }
  }
  COUT2(" > rotating " << group_to_edges.size() << " edges in total");

  // sort by group index; should only contain the flagged edges
  sort(group_to_edges.begin(),group_to_edges.end());

  if (n_duplicates) {
    // first populate duplicates in nsp, then use scale to modify appropriate cloned points
    const int ss_nsp0 = nsp;
    nsp += n_duplicates*group_to_edges.size();  // because open edges, should only be applying this to loops anyways, so for now assume 1 node per edge
    growNspData(nsp,ss_nsp0);
    sp_flag.resize(nsp);

    const double angle_dup = angle_deg/double(n_duplicates);

    vector<NewTri> newTris;
    int new_isp = ss_nsp0;
    for (vector<pair<int,uint> >::const_iterator it=group_to_edges.begin(); it!=group_to_edges.end(); ) {

      const int group = it->first;
      vector<uint> heVec;
      while ((it->first == group)&&(it!=group_to_edges.end())) {
        heVec.push_back(it->second);
        ++it;
      }
      orderHalfEdgesInVec(heVec);  // uses sp_flag, so make sure no clobbering

      vector<int> isp0Vec;
      vector<int> isp1Vec;

      // peel first loop b/c isp0 indices not ordered; others are in xsp
      sp_flag.setAll(0);
      for (vector<uint>::const_iterator hit=heVec.begin(); hit!=heVec.end(); ++hit) {
        uint ist; uchar i;
        unpackEdge(ist,i,*hit);
        const int isp0 = spost[ist][i];
        isp0Vec.push_back(isp0);
        isp1Vec.insert(isp1Vec.begin(),new_isp);  // needs to be in reverse order...
        FOR_I3 xsp[new_isp][i] = xsp[isp0][i];
        sp_flag[new_isp] = 1;
        ++new_isp;
      }
      // sp_flag should have flagged nodes just on this group
      // apply translation for the selected nodes
      rotateFlaggedNodes(axis,point,angle_dup);
      facetGap(newTris,isp0Vec,isp1Vec,true);  // HACK: revisit loop boolean

      if (n_duplicates > 1) {
        // additional loops if present
        const int n_loop_nodes = isp0Vec.size();
        for (int idup=1; idup<n_duplicates; ++idup) {
          sp_flag.setAll(0);
          for (int ino=0; ino<n_loop_nodes; ++ino) {
            isp0Vec[n_loop_nodes-ino-1] = isp1Vec[ino]; // needs to be reverse order
            isp1Vec[ino] = new_isp;
            FOR_I3 xsp[new_isp][i] = xsp[isp0Vec[n_loop_nodes-ino-1]][i];
            sp_flag[new_isp] = 1;
            ++new_isp;
          }
          rotateFlaggedNodes(axis,point,angle_dup);
          facetGap(newTris,isp0Vec,isp1Vec,true);  // HACK: revisit loop boolean
        }
      }
    }

    if (!newTris.empty()) {
      // all new tris get put into a new zone, so we don't have to keep track of which tri they were originally attached to (and for easy deletion if desired)
      const int ss_nst0 = nst;
      nst += newTris.size();
      growNstData(nst,ss_nst0);
      const int nzn0 = zoneVec.size();
      const int nsz0 = nsz++;
      zoneVec.push_back(SurfaceZone("rotated_edges"));
      szozn_i.resize(zoneVec.size()+1);
      szozn_i[zoneVec.size()] = nsz;

      for (vector<NewTri>::const_iterator it=newTris.begin(); it!=newTris.end(); ++it) {
        const int ist = ss_nst0 + (it-newTris.begin());
        FOR_I3 spost[ist][i] = it->spost[i];
        znost[ist] = nzn0;
        szost[ist] = nsz0;
      }
      clearTeost();
    }
  }
  else {
    for (vector<pair<int,uint> >::const_iterator it=group_to_edges.begin(); it!=group_to_edges.end(); ) {
      const int group = it->first;
      sp_flag.setAll(0);
      while ((it->first == group)&&(it!=group_to_edges.end())) {
        uint ist; uchar i;
        unpackEdge(ist,i,it->second);
        const int isp0 = spost[ist][i];
        const int isp1 = spost[ist][(i+1)%3];
        sp_flag[isp0]  = sp_flag[isp1] = 1;
        ++it;
      }

      // sp_flag should have flagged nodes just on this group
      // apply rotation for the selected nodes
      rotateFlaggedNodes(axis,point,angle_deg);
    }
  }
}

void SimpleSurface::scaleSelectedSubzones(const vector<int>& subzonesVec,const double sx[3],const double x0[3],const bool b_scale_centroid,const bool b_norm) {
  flagTrisFromSubzoneVec(subzonesVec);  // flagged tris in st_flag
  double ref_x0[3] = {x0[0],x0[1],x0[2]};
  if (b_scale_centroid) calcFlaggedTrisCentroid(ref_x0);  // calculate centroid of tris being scaled
  scaleFlaggedTris(sx,ref_x0,b_norm);
}

void SimpleSurface::scaleSelectedSubzonesRadially(const vector<int>& subzonesVec,const double factor,const double x0[3],const double axis[3],const bool b_scale_centroid,const bool b_norm) {
  flagTrisFromSubzoneVec(subzonesVec);  // flagged tris in st_flag
  double ref_x0[3] = {x0[0],x0[1],x0[2]};
  if (b_scale_centroid) calcFlaggedTrisCentroid(ref_x0);  // calculate centroid of tris being scaled
  scaleFlaggedTrisRadially(factor,ref_x0,axis,b_norm);
}

void SimpleSurface::scaleSurface(const double sx[3],const double x0[3],const bool b_scale_centroid,const bool b_norm) {
  sp_flag.setLength(nsp);
  sp_flag.setAll(1);
  double ref_x0[3] = {x0[0],x0[1],x0[2]};
  if (b_scale_centroid) getCentroid(ref_x0);  // entire surface's centroid already computed
  scaleFlaggedNodes(sx,ref_x0,b_norm);
}

void SimpleSurface::scaleOpenEdges(const vector<int>& edge_indices, const double sx[3], const double x0[3], const bool b_scale_norm,const bool b_centroids,const int n_duplicates) {
  ensureOpenEdgeGroups();

  if (!n_open_edge_groups) {
    WUI(WARN,"no open edge loops to process; skipping");
    return;
  }

  sp_flag.resize(nsp);

  IntFlag oeg_flag(n_open_edge_groups);

  if (edge_indices.empty()) {
    // if not specified then assume all open edge nodes should be scaled
    oeg_flag.setAll(1);
  }
  else {
    oeg_flag.setAll(0);
    for (vector<int>::const_iterator it=edge_indices.begin(); it!=edge_indices.end(); ++it) {
      oeg_flag[*it] = 1;
    }
  }

  double (*centroids)[3] = new double[n_open_edge_groups][3];
  if (!b_centroids) {
    for (int ig=0; ig<n_open_edge_groups; ++ig) {
      FOR_I3 centroids[ig][i] = x0[i];
    }
  }
  else {
    for (int ig=0; ig<n_open_edge_groups; ++ig) {
      FOR_I3 centroids[ig][i] = openEdgeGroupDataVec[ig].xc[i];
    }
  }

  // populate gr_to_edge structure
  vector<pair<int,uint> > group_to_edges;
  for (map<uint,int>::const_iterator it=oe_to_group.begin(); it!= oe_to_group.end(); ++it) {
    if (oeg_flag[it->second]) {
      group_to_edges.push_back(pair<int,uint> (it->second,it->first));
    }
  }
  COUT2(" > scaling " << group_to_edges.size() << " edges in total");

  // sort by group index; should only contain the flagged edges
  sort(group_to_edges.begin(),group_to_edges.end());

  if (n_duplicates) {
    // first populate duplicates in nsp, then use scale to modify appropriate cloned points
    const int ss_nsp0 = nsp;
    nsp += n_duplicates*group_to_edges.size();  // because open edges, should only be applying this to loops anyways, so for now assume 1 node per edge
    growNspData(nsp,ss_nsp0);
    sp_flag.resize(nsp);

    double sx_dup[3];
    if (b_scale_norm) {
      // FOR_I3 sx_dup[i] = 1.0 + (sx[i]-1.0)/double(n_duplicates);
      FOR_I3 sx_dup[i] = pow(sx[i],1.0/double(n_duplicates));
    }
    else {
      // distance based description
      FOR_I3 sx_dup[i] = sx[i]/double(n_duplicates);
    }

    vector<NewTri> newTris;
    int new_isp = ss_nsp0;
    for (vector<pair<int,uint> >::const_iterator it=group_to_edges.begin(); it!=group_to_edges.end(); ) {

      const int group = it->first;
      vector<uint> heVec;
      while ((it->first == group)&&(it!=group_to_edges.end())) {
        heVec.push_back(it->second);
        ++it;
      }
      orderHalfEdgesInVec(heVec);  // uses sp_flag, so make sure no clobbering

      vector<int> isp0Vec;
      vector<int> isp1Vec;

      // peel first loop b/c isp0 indices not ordered; others are in xsp
      sp_flag.setAll(0);
      for (vector<uint>::const_iterator hit=heVec.begin(); hit!=heVec.end(); ++hit) {
        uint ist; uchar i;
        unpackEdge(ist,i,*hit);
        const int isp0 = spost[ist][i];
        isp0Vec.push_back(isp0);
        isp1Vec.insert(isp1Vec.begin(),new_isp);  // needs to be in reverse order...
        FOR_I3 xsp[new_isp][i] = xsp[isp0][i];
        sp_flag[new_isp] = 1;
        ++new_isp;
      }
      // sp_flag should have flagged nodes just on this group
      // apply scaling for the selected nodes
      scaleFlaggedNodes(sx_dup,centroids[group],b_scale_norm);
      facetGap(newTris,isp0Vec,isp1Vec,true);  // HACK: revisit loop boolean

      if (n_duplicates > 1) {
        // additional loops if present
        const int n_loop_nodes = isp0Vec.size();
        for (int idup=1; idup<n_duplicates; ++idup) {
          sp_flag.setAll(0);
          for (int ino=0; ino<n_loop_nodes; ++ino) {
            isp0Vec[n_loop_nodes-ino-1] = isp1Vec[ino]; // needs to be reverse order
            isp1Vec[ino] = new_isp;
            FOR_I3 xsp[new_isp][i] = xsp[isp0Vec[n_loop_nodes-ino-1]][i];
            sp_flag[new_isp] = 1;
            ++new_isp;
          }
          scaleFlaggedNodes(sx_dup,centroids[group],b_scale_norm);
          facetGap(newTris,isp0Vec,isp1Vec,true);  // HACK: revisit loop boolean
        }
      }
    }

    if (!newTris.empty()) {
      const int ss_nst0 = nst;
      nst += newTris.size();
      growNstData(nst,ss_nst0);
      const int nzn0 = zoneVec.size();
      const int nsz0 = nsz++;
      zoneVec.push_back(SurfaceZone("scaled_edges"));
      szozn_i.resize(zoneVec.size()+1);
      szozn_i[zoneVec.size()] = nsz;

      for (vector<NewTri>::const_iterator it=newTris.begin(); it!=newTris.end(); ++it) {
        const int ist = ss_nst0 + (it-newTris.begin());
        FOR_I3 spost[ist][i] = it->spost[i];
        znost[ist] = nzn0;
        szost[ist] = nsz0;
      }
      clearTeost();
    }
  }
  else {
    for (vector<pair<int,uint> >::const_iterator it=group_to_edges.begin(); it!=group_to_edges.end(); ) {
      const int group = it->first;
      sp_flag.setAll(0);
      while ((it->first == group)&&(it!=group_to_edges.end())) {
        uint ist; uchar i;
        unpackEdge(ist,i,it->second);
        const int isp0 = spost[ist][i];
        const int isp1 = spost[ist][(i+1)%3];
        sp_flag[isp0]  = sp_flag[isp1] = 1;
        ++it;
      }

      // sp_flag should have flagged nodes just on this group
      // apply scaling for the selected nodes
      scaleFlaggedNodes(sx,centroids[group],b_scale_norm);
    }
  }
  DELETE(centroids);
}

void SimpleSurface::forceOrthogonalityOpenEdges(const vector<int>& open_edge_groups,const double dir[3]) {
  COUT2("forceOrthogonalityOpenEdges()");

  double dir_[3] = {dir[0],dir[1],dir[2]};
  NORMALIZE(dir_);

  sp_flag.resize(nsp);

  for (int oeg=0, lim=open_edge_groups.size(); oeg<lim; ++oeg) {
    const int index = open_edge_groups[oeg];
    if (index < 0 || (index>n_open_edge_groups)) {
      continue;
    }
    else {
      COUT2(" > working on loop: " << index);
      double n0[3];
      computeOpenEdgeLoopNormal(n0,index);
      NORMALIZE(n0);
      const double angle = fabs(acos(DOT_PRODUCT(n0,dir_)) - 0.5*M_PI);
      if (angle > 0.5*M_PI/180.0) {
        COUT2("    > deviation from orthogonality is large (greater than 0.5 degrees); this could cause problems with the node translations (i.e., folding)");
      }
      COUT2("    > angle dev between n,dir:" << COUT_VEC(n0) << "," << COUT_VEC(dir) << "," << angle);

      // flag nodes on group
      sp_flag.setAll(0);
      for (map<uint,int>::iterator it=oe_to_group.begin(); it!=oe_to_group.end(); ++it) {
        if (it->second == index) {
          uint ist; uchar i;
          unpackEdge(ist,i,it->first);
          const int isp0 = spost[ist][i];
          const int isp1 = spost[ist][(i+1)%3];
          sp_flag[isp0] = 1;
          sp_flag[isp1] = 1;
        }
      }

      // now project onto new plane
      double n_new[3];
      FOR_I3 n_new[i] = n0[i] - n0[i]*dir_[i];
      NORMALIZE(n_new);
      double dist_sum = 0.0;
      for (int isp=0; isp<nsp; ++isp) {
        if (sp_flag[isp]) {
          // translate onto properly orthogonal plane
          const double dist = MiscUtils::getPointToPlaneDist(xsp[isp],openEdgeGroupDataVec[index].xc,n_new);  // signed
          dist_sum += dist;
          FOR_I3 xsp[isp][i] -= dist*n_new[i];
        }
      }
      COUT2("    > total nodes translated: " << sp_flag.count());
      COUT2("    > average distance translated: " << dist_sum/double(sp_flag.count()));
    }

  }

  // nodes change so metadata needsd to be recomputed
  clearDynamicEdgeGroups();
  clearOpenEdgeGroups();
  clearSubzoneData();
  clearZoneData();
  b_centroid = false;
  b_rmax = false;
  b_boundingBox = false;
}

void SimpleSurface::copyTranslateSelectedSubzones(const vector<int>& subzonesVec,const int n_copies,const double dx[3]) {

  //const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;
  copySelectedSubzones(subzonesVec,n_copies);
  const int tris_per_copy = (nst-ss_nst0)/n_copies;
  assert((nst-ss_nst0)%tris_per_copy==0);

  st_flag.resize(nst);
  for (int icopy=0; icopy<n_copies; ++icopy) {
    st_flag.setAll(0);
    for (int ist=0; ist < tris_per_copy; ++ist) {
      st_flag[ss_nst0 + icopy*tris_per_copy + ist] = 1;
    }

    double this_dx[3];
    FOR_I3 this_dx[i] = (icopy+1)*dx[i];
    translateFlaggedTris(this_dx);
  }
}

void SimpleSurface::copyRotateSelectedSubzones(const vector<int>& subzonesVec,const int n_copies,const double axis[3],const double point[3],const double angle_deg) {

  //const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;
  copySelectedSubzones(subzonesVec,n_copies);
  const int tris_per_copy = (nst-ss_nst0)/n_copies;
  assert((nst-ss_nst0)%tris_per_copy==0);

  st_flag.resize(nst);
  for (int icopy=0; icopy<n_copies; ++icopy) {
    st_flag.setAll(0);
    for (int ist=0; ist < tris_per_copy; ++ist) {
      st_flag[ss_nst0 + icopy*tris_per_copy + ist] = 1;
    }

    const double this_angle = (icopy+1)*angle_deg;
    rotateFlaggedTris(axis,point,this_angle);
  }
}

void SimpleSurface::copySelectedSubzones(const vector<int>& subzonesVec,const int n_copies) {

  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;
  const int ss_nsz0 = nsz;
  if (ss_nsp0==0 || ss_nst0==0 || subzonesVec.empty()) {
    WUI(WARN,"no surface currently available or selected; skipping");
    return;
  }

  // count new tris and nodes to be added
  int sp_new,st_new;
  indexNodesAndTrisFromSubzoneVec(sp_new,st_new,subzonesVec);
  nsp = ss_nsp0 + sp_new*n_copies;
  nst = ss_nst0 + st_new*n_copies;

  // need to keep track of which subzones are hidden as well as which are selected
  multimap<int,int> track_szn_local;  // positive is selected, negative are hidden (+/-1 indexed)
  sz_flag.resize(nsz);
  sz_flag.setAll(0);
  for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
    for (int isz=szozn_i[izn]; isz<szozn_i[izn+1]; ++isz) sz_flag[isz] = izn;  // sz_flag holds zone now
  }

  // hidden subzones
  if (!hiddenSubzoneVec.empty()) {
    b_update_hidden_subzones = true;
    for (vector<int>::iterator it=hiddenSubzoneVec.begin(); it!=hiddenSubzoneVec.end(); ++it) {
      const int isz = *it;
      const int local_isz = isz-szozn_i[sz_flag[isz]];
      track_szn_local.insert(pair<int,int> (sz_flag[isz],-(local_isz+1)) );
    }
  }

  // selected subzones
  for (vector<int>::const_iterator it=subzonesVec.begin(); it!=subzonesVec.end(); ++it) {
    const int isz = *it;
    const int local_isz = isz-szozn_i[sz_flag[isz]];
    track_szn_local.insert(pair<int,int> (sz_flag[isz],(local_isz+1)) );
  }

  // convert szost to local indices
  sz_flag.setAll(0);
  for (vector<int>::const_iterator it=subzonesVec.begin(); it!=subzonesVec.end(); ++it) sz_flag[*it] = 1;
  int sz_new = sz_flag.count();
  nsz = ss_nsz0 + sz_new*n_copies;
  // convert global to local subzones for old tris
  for (int ist=0; ist<ss_nst0; ++ist) {
    szost[ist] -= szozn_i[znost[ist]];
  }

  growNspData(nsp,ss_nsp0);
  growNstData(nst,ss_nst0);

  // populate new data
  for (int icopy=0; icopy<n_copies; ++icopy) {
    // nodes
    int sp_count=0;
    for (int isp=0; isp<ss_nsp0; ++isp) {
      if (sp_flag[isp] != -1) {
        FOR_I3 xsp[ss_nsp0 + sp_new*icopy + sp_flag[isp]][i] = xsp[isp][i];
        ++sp_count;
      }
    }
    assert(sp_count == sp_new);

    // tri structures
    int st_count=0;
    for (int ist=0; ist<ss_nst0; ++ist) {
      if (st_flag[ist] != -1) {
        assert(st_count==st_flag[ist]);
        const int ist_new = ss_nst0 + st_new*icopy + st_count++;
        FOR_I3 {
          spost[ist_new][i] = ss_nsp0 + sp_new*icopy + sp_flag[spost[ist][i]];
        }
        znost[ist_new] = znost[ist];
        szost[ist_new] = szozn_i[znost[ist]+1]-szozn_i[znost[ist]] + icopy;  // each copy to own subzone
      }
    }
    assert(st_count == st_new);
  }

  // if we had periodic data, we need to reset it here
  if (pbi) {
    clearPeriodicity();
    WUI(WARN,"We've added copies of the surface, but things changed enough that you'll need to reset periodicity.");
  }

  // is any of this overkill?
  clearStoszSzosz();
  clearTeost();
  clearNonManifoldData();

  // rebuild sz stuff
  buildSzoznFromLocalSzost();
  localSzostToGlobal();

  // set hidden & selected subzones
  if (!track_szn_local.empty()) {
    hiddenSubzoneVec.clear();
    selectedSubzoneVec.clear();
    for (multimap<int,int>::iterator it=track_szn_local.begin(); it!=track_szn_local.end();) {
      const int izn = it->first;
      const int offset = szozn_i[izn];
      const int isz_i = it->second;

      do {
        if (isz_i < 0) {
          hiddenSubzoneVec.push_back((-isz_i-1) + offset);
        }
        else if (isz_i > 0) {
          selectedSubzoneVec.push_back((isz_i-1) + offset);
        }
        else { assert(0);}  // shouldn't have inserted

        ++it;
      } while (it->first == izn && it!=track_szn_local.end());
    }
  }

  pruneEmptyZonesAndSubzones();
}

void SimpleSurface::roughenFlaggedSubzones(const double range[2],const int n_filter) {

  // flag nodes on the interior of the flagged subzones...

  sp_flag.resize(nsp);
  sp_flag.setAll(0);
  for (int ist = 0; ist < nst; ++ist) {
    if (sz_flag[szost[ist]]) {
      // set the 2 bit in all nodes of this tri...
      FOR_I3 {
        const int isp = spost[ist][i];
        sp_flag[isp] |= 2;
      }
    }
    else {
      // set the one bit in all nodes of this tri...
      FOR_I3 {
        const int isp = spost[ist][i];
        sp_flag[isp] |= 1;
      }
    }
  }

  // now the 2's can be moved. 1's and 3's should be left alone...

  int nsp_roughen = 0;
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] == 2) {
      sp_flag[isp] = nsp_roughen++;
    }
    else {
      sp_flag[isp] = -1;
    }
  }

  if (nsp_roughen == 0) {
    WUI(WARN," > no surface nodes available to perturb; skipping");
    return;
  }

  // reduce the tri normals to the flagged nodes...

  double (*n_sp)[3] = new double[nsp_roughen][3];
  double *A_sp = new double[nsp_roughen];
  for (int isp = 0; isp < nsp_roughen; ++isp) {
    FOR_I3 n_sp[isp][i] = 0.0;
    A_sp[isp] = 0.0;
  }
  FOR_IST {
    const int isp0 = spost[ist][0];
    const int isp1 = spost[ist][1];
    const int isp2 = spost[ist][2];
    if ((sp_flag[isp0] >= 0)||(sp_flag[isp1] >= 0)||(sp_flag[isp2] >= 0)) {
      const double n[3] = TRI_NORMAL_2(xsp[isp0],xsp[isp1],xsp[isp2]);
      const double n_mag = MAG(n);
      if (n_mag > 0.0) {
        if (sp_flag[isp0] >= 0) {
          FOR_I3 n_sp[sp_flag[isp0]][i] += n[i];
          A_sp[sp_flag[isp0]] += n_mag;
        }
        if (sp_flag[isp1] >= 0) {
          FOR_I3 n_sp[sp_flag[isp1]][i] += n[i];
          A_sp[sp_flag[isp1]] += n_mag;
        }
        if (sp_flag[isp2] >= 0) {
          FOR_I3 n_sp[sp_flag[isp2]][i] += n[i];
          A_sp[sp_flag[isp2]] += n_mag;
        }
      }
    }
  }

  // init pertubation to uniform random distribution...

  double *dn_sp = new double[nsp_roughen];
  for (int isp = 0; isp < nsp_roughen; ++isp)
    dn_sp[isp] = MiscUtils::uniformRand(range[0],range[1]);

  // smooth and rescale perturbation...

  double *dn_tmp = new double[nsp_roughen];
  for (int ii = 0; ii < n_filter; ++ii) {

    // distribute node perturbations to other nodes of tri, weight by face area...

    for (int isp = 0; isp < nsp_roughen; ++isp)
      dn_tmp[isp] = 0.0;
    FOR_IST {
      const int isp0 = spost[ist][0];
      const int isp1 = spost[ist][1];
      const int isp2 = spost[ist][2];
      if ((sp_flag[isp0] >= 0)||(sp_flag[isp1] >= 0)||(sp_flag[isp2] >= 0)) {
        const double n[3] = TRI_NORMAL_2(xsp[isp0],xsp[isp1],xsp[isp2]);
        const double n_mag = MAG(n);
        if (n_mag > 0.0) {
          double tmp = 0.0;
          if (sp_flag[isp0] >= 0) tmp += dn_sp[sp_flag[isp0]];
          if (sp_flag[isp1] >= 0) tmp += dn_sp[sp_flag[isp1]];
          if (sp_flag[isp2] >= 0) tmp += dn_sp[sp_flag[isp2]];
          tmp /= 3.0;

          if (sp_flag[isp0] >= 0) dn_tmp[sp_flag[isp0]] += n_mag*tmp;
          if (sp_flag[isp1] >= 0) dn_tmp[sp_flag[isp1]] += n_mag*tmp;
          if (sp_flag[isp2] >= 0) dn_tmp[sp_flag[isp2]] += n_mag*tmp;
        }
      }
    }
    for (int isp = 0; isp < nsp_roughen; ++isp) {
      if (A_sp[isp] > 0.0)
        dn_sp[isp] = dn_tmp[isp]/A_sp[isp];
    }

    // compute stats...

    double _mean = 0.0;
    double _min = HUGE_VAL;
    double _max = -HUGE_VAL;
    double _M2 = 0.0;
    for (int isp = 0; isp < nsp_roughen; ++isp) {
      double delta = dn_sp[isp] - _mean;
      _mean += delta/nsp_roughen;
      _M2 += delta*(dn_sp[isp] - _mean);
      _min = min(dn_sp[isp],_min);
      _max = max(dn_sp[isp],_max);
    }

    double _std = sqrt(_M2/(nsp_roughen-1));
    COUT2("       > pre-scale stats (mean,std,min,max): " << _mean << " " << _std << " " << _min << " " << _max);

    const double width = range[1]-range[0];
    const double midpoint = 0.5*(range[0]+range[1]);
    const double _width = _max - _min;
    const double _midpoint = 0.5*(_max + _min);
    assert(_width > 0.0);
    const double scale_fac = width/_width;

    // this implementation works to preserve the min and max of the distribution above all else, so mean may float a bit...

    for (int isp = 0; isp < nsp_roughen; ++isp) {
      dn_sp[isp] -= _midpoint;
      dn_sp[isp] *= scale_fac;
      dn_sp[isp] += midpoint;
    }

    // post statistics...
    _mean = 0.0;
    _min = HUGE_VAL;
    _max = -HUGE_VAL;
    _M2 = 0.0;

    for (int isp = 0; isp < nsp_roughen; ++isp) {
      double delta = dn_sp[isp] - _mean;
      _mean += delta/nsp_roughen;
      _M2 += delta*(dn_sp[isp] - _mean);
      _min = min(dn_sp[isp],_min);
      _max = max(dn_sp[isp],_max);
    }
    assert(nsp_roughen > 1);
    _std = sqrt(_M2/(nsp_roughen-1));
    COUT2("       > post-scale stats (mean,std,min,max): " << _mean << " " << _std << " " << _min << " " << _max);

  }
  delete[] dn_tmp;

  // apply perturbation to flagged points...

  FOR_ISP {
    if (sp_flag[isp] >= 0) {
      const int isp_roughen = sp_flag[isp];
      if (A_sp[isp_roughen] > 0.0)
        FOR_I3 xsp[isp][i] += dn_sp[isp_roughen]*n_sp[isp_roughen][i]/A_sp[isp_roughen];
    }
  }
  delete[] A_sp;
  delete[] n_sp;
  delete[] dn_sp;

}

