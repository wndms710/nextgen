#include "SimpleSurface.hpp"
#include "GeomUtils.hpp"

void SimpleSurface::intersectFlaggedSubzones(const string& mode) {

  cout << "SimpleSurface::intersectFlaggedSubzones: mode: " << mode << endl;

  int mode_int = 0; // default subtract mode...
  if ((mode == "SUBTRACT")||(mode == "subtract")) {
    mode_int = 0;
  }
  else if ((mode == "ADD")||(mode == "add")) {
    mode_int = 1;
  }
  else {
    WUI(WARN,"INTERSECT failed: unrecognized mode: " << mode);
    return;
  }

  // coming into this routine we expect 2 sets of flagged subzones:
  // sz_flag == 0 in one set, and
  // sz_flag == 1 in the other
  // and
  // sz_flag == -1 in tris not being considered for intersection

  //writeTecplot("debug.dat");
  //bool debug = false;

  // start by building the bbox of the 2 regions, and figure out
  // how many tris are participating...

  double bbmin0[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
  double bbmax0[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
  double bbmin1[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
  double bbmax1[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };

  int nst0 = 0;
  int nst1 = 0;
  FOR_IST {
    const int isz = szost[ist];
    if (sz_flag[isz] == 0) {
      ++nst0;
      FOR_I3 {
        const int isp = spost[ist][i];
        FOR_J3 bbmin0[j] = min(bbmin0[j],xsp[isp][j]);
        FOR_J3 bbmax0[j] = max(bbmax0[j],xsp[isp][j]);
      }
    }
    else if (sz_flag[isz] == 1) {
      ++nst1;
      FOR_I3 {
        const int isp = spost[ist][i];
        FOR_J3 bbmin1[j] = min(bbmin1[j],xsp[isp][j]);
        FOR_J3 bbmax1[j] = max(bbmax1[j],xsp[isp][j]);
      }
    }
    else {
      assert(sz_flag[isz] == -1);
    }
  }

  cout << " > flagged tri count:" << endl;
  cout << "   > nst0: " << nst0 << " bbox: " << COUT_VEC(bbmin0) << " : " << COUT_VEC(bbmax0) << endl;
  cout << "   > nst1: " << nst1 << " bbox: " << COUT_VEC(bbmin1) << " : " << COUT_VEC(bbmax1) << endl;

  // now, the only tris that can possibly intersect are the ones touching the 
  // common part of the bbox...
  
  double bbmin[3],bbmax[3];
  FOR_I3 {
    bbmin[i] = max(bbmin0[i],bbmin1[i]);
    bbmax[i] = min(bbmax0[i],bbmax1[i]);
  }
  
  // check that there is some common volume...

  if ((bbmax[0] < bbmin[0])||(bbmax[1] < bbmin[1])||(bbmax[2] < bbmin[2])) {
    WUI(INFO,"passed subzones do not intersect - no change to surface");
    return;
  }
  
  // now store tris in 2 vectors...
  vector<int> fa0Vec;
  vector<int> fa1Vec;
  
  FOR_IST {
    const int isz = szost[ist];
    if (sz_flag[isz] == 0) {
      if ((max(xsp[spost[ist][0]][0],max(xsp[spost[ist][1]][0],xsp[spost[ist][2]][0])) < bbmin[0])||
          (min(xsp[spost[ist][0]][0],min(xsp[spost[ist][1]][0],xsp[spost[ist][2]][0])) > bbmax[0]))
        continue;
      if ((max(xsp[spost[ist][0]][1],max(xsp[spost[ist][1]][1],xsp[spost[ist][2]][1])) < bbmin[1])||
          (min(xsp[spost[ist][0]][1],min(xsp[spost[ist][1]][1],xsp[spost[ist][2]][1])) > bbmax[1]))
        continue;
      if ((max(xsp[spost[ist][0]][2],max(xsp[spost[ist][1]][2],xsp[spost[ist][2]][2])) < bbmin[2])||
          (min(xsp[spost[ist][0]][2],min(xsp[spost[ist][1]][2],xsp[spost[ist][2]][2])) > bbmax[2]))
        continue;
      // if we made it here, this one is part of nfa0...
      fa0Vec.push_back(ist);
    }
    else if (sz_flag[isz] == 1) {
      if ((max(xsp[spost[ist][0]][0],max(xsp[spost[ist][1]][0],xsp[spost[ist][2]][0])) < bbmin[0])||
          (min(xsp[spost[ist][0]][0],min(xsp[spost[ist][1]][0],xsp[spost[ist][2]][0])) > bbmax[0]))
        continue;
      if ((max(xsp[spost[ist][0]][1],max(xsp[spost[ist][1]][1],xsp[spost[ist][2]][1])) < bbmin[1])||
          (min(xsp[spost[ist][0]][1],min(xsp[spost[ist][1]][1],xsp[spost[ist][2]][1])) > bbmax[1]))
        continue;
      if ((max(xsp[spost[ist][0]][2],max(xsp[spost[ist][1]][2],xsp[spost[ist][2]][2])) < bbmin[2])||
          (min(xsp[spost[ist][0]][2],min(xsp[spost[ist][1]][2],xsp[spost[ist][2]][2])) > bbmax[2]))
        continue;
      // if we made it here, this one is part of nfa0...
      fa1Vec.push_back(ist);
    }
  }
  
  const int nfa0 = fa0Vec.size();
  const int nfa1 = fa1Vec.size();
  
  // put nfa0 in an adt...
  double (*bbmin_fa0)[3] = new double[nfa0][3];
  double (*bbmax_fa0)[3] = new double[nfa0][3];
  for (int ifa0 = 0; ifa0 < nfa0; ++ifa0) {
    FOR_I3 bbmin_fa0[ifa0][i] = HUGE_VAL;
    FOR_I3 bbmax_fa0[ifa0][i] = -HUGE_VAL;
    const int ist = fa0Vec[ifa0]; assert((ist >= 0)&&(ist < nst));
    FOR_I3 {
      const int isp = spost[ist][i];
      FOR_J3 {
	bbmin_fa0[ifa0][j] = min(bbmin_fa0[ifa0][j],xsp[isp][j]);
	bbmax_fa0[ifa0][j] = max(bbmax_fa0[ifa0][j],xsp[isp][j]);
      }
    }
  }
  Adt<double> * adt = new Adt<double>(nfa0,bbmin_fa0,bbmax_fa0);
  delete[] bbmin_fa0;
  delete[] bbmax_fa0;
  
  // figure out what faces intersect based on bbox check...
  vector<pair<int,int> > facePairVec;
  vector<int> intVec;
  for (int ifa1 = 0; ifa1 < nfa1; ++ifa1) {
    // we need the bbox of this face to query the adt...
    double bbmin_fa1[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
    double bbmax_fa1[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
    const int ist = fa1Vec[ifa1]; assert((ist >= 0)&&(ist < nst));
    FOR_I3 {
      const int isp = spost[ist][i];
      FOR_J3 {
	bbmin_fa1[j] = min(bbmin_fa1[j],xsp[isp][j]);
	bbmax_fa1[j] = max(bbmax_fa1[j],xsp[isp][j]);
      }
    }
    // now grab the faces may intersect based on the bbox overlap...
    assert(intVec.empty());
    adt->buildListForBBox(intVec,bbmin_fa1,bbmax_fa1);
    if (!intVec.empty()) {
      for (vector<int>::const_iterator it=intVec.begin(); it!=intVec.end(); ++it) {
        facePairVec.push_back(pair<int,int>(*it,ifa1));
      }
      intVec.clear();
    }
  }

  // done with adt...
  delete adt; adt = NULL;
  
  cout << " > got nfa0,nfa1: " << nfa0 << " " << nfa1 << " facePairVec.size(): " << facePairVec.size() << endl;
  
  map<const pair<int,int>,int> edgeMap;
  vector<pair<int,int> > stoedVec;
  vector<IntersectionData> intersectionVec;
  vector<pair<int,int> > linkVec; // new edges introduced by tri-tri intersections
  vector<pair<int,int> > stolkVec; // tris on either link...
  double * x0[3];
  double * x1[3];
  int idata[6];
  double ddata[4];
  for (vector<pair<int,int> >::const_iterator fpit=facePairVec.begin(); fpit!=facePairVec.end(); ++fpit) {
    
    const int ist0 = fa0Vec[fpit->first];
    const int ist1 = fa1Vec[fpit->second];
    assert(ist1 != ist0); // ist's should always be unique

    FOR_I3 x0[i] = xsp[spost[ist0][i]];
    FOR_I3 x1[i] = xsp[spost[ist1][i]];
    if (const int n = MiscUtils::calcTriTriIntersection(idata,ddata,x0,x1)) {
      switch (n) {
      case -2:
        assert(0);
      case 1:
        processTriTriIntersection(edgeMap,stoedVec,intersectionVec,idata,ddata,ist0,ist1);
        break;
      case 2:
        // add a link using the index of the 2 intersections...
        linkVec.push_back(pair<int,int>(intersectionVec.size(),intersectionVec.size()+1));
        stolkVec.push_back(pair<int,int>(ist0,ist1));
        processTriTriIntersection(edgeMap,stoedVec,intersectionVec,idata,ddata,ist0,ist1);
        processTriTriIntersection(edgeMap,stoedVec,intersectionVec,idata+3,ddata+2,ist0,ist1);
        break;
      default:
        assert(n == -1);
      }
    }

  }
  
  if (intersectionVec.empty()) {
    WUI(INFO,"passed subzones do not intersect - no change to surface");
    return;
  }
  
  cout << " > intersectionVec.size(): " << intersectionVec.size() << " linkVec.size(): " << linkVec.size() << endl;
  
  // a sufficient check for closed loops is that stoedVec.second is set...
  for (int ii = 0; ii < stoedVec.size(); ++ii) {
    if (stoedVec[ii].second == -1) {
      cout << " > intersections do not form closed loop" << endl;
      assert(0);
    }
  }

  cout << " > intersections appear to form one or more closed loops" << endl;
  
  // use the map to build the spoedVec, AND set bits in a global st_flag to 
  // indicate which edges of the tris have been replaced...
  st_flag.resize(nst);
  st_flag.setAll(0);
  vector<pair<int,int> > spoedVec(stoedVec.size());
  for (map<const pair<int,int>,int>::iterator it = edgeMap.begin(); it != edgeMap.end(); ++it) {
    const int ied = it->second;
    spoedVec[ied].first = it->first.first;
    spoedVec[ied].second = it->first.second;
    // on the left side of this edge, the spost should be aligned with one of the edges...
    {
      const int ist0 = stoedVec[ied].first; assert((ist0 >= 0)&&(ist0 < nst));
      if ((spoedVec[ied].first == spost[ist0][0])&&(spoedVec[ied].second == spost[ist0][1])) {
        st_flag[ist0] |= (1<<0);
      }
      else if ((spoedVec[ied].first == spost[ist0][1])&&(spoedVec[ied].second == spost[ist0][2])) {
        st_flag[ist0] |= (1<<1);
      }
      else {
        assert((spoedVec[ied].first == spost[ist0][2])&&(spoedVec[ied].second == spost[ist0][0]));
        st_flag[ist0] |= (1<<2);
      }
    }
    // on the right side of this edge, the spost should be reverse-aligned with one of the edges...
    {
      const int ist1 = stoedVec[ied].second; assert((ist1 >= 0)&&(ist1 < nst));
      if ((spoedVec[ied].first == spost[ist1][1])&&(spoedVec[ied].second == spost[ist1][0])) {
        st_flag[ist1] |= (1<<0);
      }
      else if ((spoedVec[ied].first == spost[ist1][2])&&(spoedVec[ied].second == spost[ist1][1])) {
        st_flag[ist1] |= (1<<1);
      }
      else {
        assert((spoedVec[ied].first == spost[ist1][0])&&(spoedVec[ied].second == spost[ist1][2]));
        st_flag[ist1] |= (1<<2);
      }
    }
  }
  edgeMap.clear();
  
  /*
  // write intersections...
  {
  FILE * fp = fopen("xi.dat","w");
  for (int ii = 0; ii < intersectionVec.size(); ++ii) {
  double xi[3];
  intersectionVec[ii].calcXi(xi,xsp,spoedVec);
  fprintf(fp,"%18.15le %18.15le %18.15le\n",xi[0],xi[1],xi[2]);
  }
  fclose(fp);
  }
  */
  
  // sort the intersections so duplicates are contiguous...
  for (int ii = 0; ii < intersectionVec.size(); ++ii) 
    intersectionVec[ii].ino = ii;
  
  sort(intersectionVec.begin(),intersectionVec.end());
  
  int * order = new int[intersectionVec.size()];
  for (int ii = 0; ii < intersectionVec.size(); ++ii) 
    order[intersectionVec[ii].ino] = ii;
  
  // adjust the linkVec references based on the new intersectionVec order. These
  // get converted to new isp indices later...
  for (int ii = 0; ii < linkVec.size(); ++ii) {
    linkVec[ii].first = order[linkVec[ii].first];
    linkVec[ii].second = order[linkVec[ii].second];
  }
  delete[] order;
  
  // add new nodes...
  vector<pair<pair<int,double>,int> > cutEdgeDataVec; // ((edge,wgt),node)
  vector<double> xspNewVec;
  int nsp_new = nsp;
  int ii_prev = -1;
  for (int ii = 0; ii < intersectionVec.size(); ++ii) {
    // skip duplicates...
    if ((ii_prev == -1)||
        (intersectionVec[ii].kind != intersectionVec[ii_prev].kind)||
        (intersectionVec[ii].idata[0] != intersectionVec[ii_prev].idata[0])||
        (intersectionVec[ii].idata[1] != intersectionVec[ii_prev].idata[1])) {
      ii_prev = ii;
      switch (intersectionVec[ii].kind) {
      case NODE_NODE_INTERSECTION:
        {
          assert(0);
          // TODO: when you hit this, use -ve indexing here for the node maybe, plus a node-node pair vector
          /*
          // combine nodes in idata[0] and idata[1]...
          int8 ino0 = no_flag[intersectionVec[ii].idata[0]];
          while (ino0 != no_flag[ino0])
          ino0 = no_flag[ino0];
          int8 ino1 = no_flag[intersectionVec[ii].idata[1]];
          while (ino1 != no_flag[ino1])
          ino1 = no_flag[ino1];
          no_flag[intersectionVec[ii].ino] = no_flag[ino0] = no_flag[ino1] = min(ino0,ino1);
          // and add this to a special vec of node-node intersections...
          intersectionNodeNode.push_back(min(ino0,ino1));
          */
        }
        break;
      case NODE_EDGE_INTERSECTION:
        {
          // use the node to define the x...
          intersectionVec[ii].ino = intersectionVec[ii].idata[0]; // the isp
          const int ied = intersectionVec[ii].idata[1];
          const double wgt = intersectionVec[ii].ddata[1];
          cutEdgeDataVec.push_back(pair<pair<int,double>,int>(pair<int,double>(ied,wgt),intersectionVec[ii].idata[0]));
        }
        break;
      case NODE_FACE_INTERSECTION:
        {
          // use the node to define the x...
          intersectionVec[ii].ino = intersectionVec[ii].idata[0];
        }
        break;
      case EDGE_EDGE_INTERSECTION:
        {
          // this requires a new node. Use the simple average of the two edge intersections...
          const int ied0 = intersectionVec[ii].idata[0];
          const double wgt0 = intersectionVec[ii].ddata[0];
          const int ied1 = intersectionVec[ii].idata[1];
          const double wgt1 = intersectionVec[ii].ddata[1];
          const int isp_new = nsp_new++;
          intersectionVec[ii].ino = isp_new;
          FOR_I3 {
            const double x = 0.5*( wgt0*xsp[spoedVec[ied0].second][i] + (1.0-wgt0)*xsp[spoedVec[ied0].first][i] +
                                   wgt1*xsp[spoedVec[ied1].second][i] + (1.0-wgt1)*xsp[spoedVec[ied1].first][i] );
            xspNewVec.push_back(x);
          }
          cutEdgeDataVec.push_back(pair<pair<int,double>,int>(pair<int,double>(ied0,wgt0),isp_new));
          cutEdgeDataVec.push_back(pair<pair<int,double>,int>(pair<int,double>(ied1,wgt1),isp_new));
        }
        break;
      case EDGE_FACE_INTERSECTION:
        {
          const int ied = intersectionVec[ii].idata[0];
          const double wgt = intersectionVec[ii].ddata[0];
          const int isp_new = nsp_new++;
          intersectionVec[ii].ino = isp_new;
          FOR_I3 {
            const double x = wgt*xsp[spoedVec[ied].second][i] + (1.0-wgt)*xsp[spoedVec[ied].first][i];
            xspNewVec.push_back(x);
          }
          cutEdgeDataVec.push_back(pair<pair<int,double>,int>(pair<int,double>(ied,wgt),isp_new));
        }
        break;
      default:
        assert(0);
      }
    }
    else {
      assert(ii_prev != -1);
      intersectionVec[ii].ino = intersectionVec[ii_prev].ino;
    }
  }

  // now update nooed for ALL edges to the no_flag value...
  // TODO: why have these changed? they have not unless node-node above, so skip
  // and consider later...
  /*
    for (int ied = 0,ned=edgeVec.size(); ied < ned; ++ied) {
    FOR_I2 {
    const int ino = edgeVec[ied].nooed[i];
    assert((no_flag[ino] >= 0)&&(no_flag[ino] < int8(nno_new)));
    edgeVec[ied].nooed[i] = int(no_flag[ino]);
    }
    }
    for (int ied = 0,ned=intersectionEdgeVec.size(); ied < ned; ++ied) {
    FOR_I2 {
    const int ino = intersectionEdgeVec[ied].nooed[i];
    assert((no_flag[ino] >= 0)&&(no_flag[ino] < int8(nno_new)));
    intersectionEdgeVec[ied].nooed[i] = int(no_flag[ino]);
    }
    }
  */

  // and cut the edges that have intersections...
  sort(cutEdgeDataVec.begin(),cutEdgeDataVec.end());
  int ied_current = -1;
  int isp_current = -1;
  for (int ii = 0,nii=cutEdgeDataVec.size(); ii < nii; ++ii) {
    const int ied = cutEdgeDataVec[ii].first.first;
    if (ied != ied_current) {
      // this is the first of a new edge. Complete the old edge first...
      if (ied_current >= 0) {
        assert(isp_current >= 0);
        spoedVec[ied_current].first = isp_current;
      }
      ied_current = ied;
      isp_current = spoedVec[ied].first;
    }
    const int isp_new = cutEdgeDataVec[ii].second;
    assert(isp_new != isp_current);
    spoedVec.push_back(pair<int,int>(isp_current,isp_new));
    stoedVec.push_back(pair<int,int>(stoedVec[ied_current].first,stoedVec[ied_current].second));
    isp_current = isp_new;
  }
  cutEdgeDataVec.clear();
  // complete last edge (the original edge)...
  if (ied_current >= 0) {
    assert(isp_current >= 0);
    spoedVec[ied_current].first = isp_current;
  }
  
  // now we can walk...
  
  sp_flag.resize(nsp_new);
  sp_flag.setAll(0);
  set<pair<int,int> > forwardFaNoSet;
  set<pair<int,int> > backwardFaNoSet;
  for (int ii = 0,nii=linkVec.size(); ii < nii; ++ii) {
    // the link is currently referencing the index of the intersectionVec. Convert
    // to a isp index in intersectionVec[].ino...
    linkVec[ii].first = intersectionVec[linkVec[ii].first].ino;
    linkVec[ii].second = intersectionVec[linkVec[ii].second].ino;
    assert(linkVec[ii].first != linkVec[ii].second);
    // these nodes stay...
    sp_flag[linkVec[ii].first] = -2;
    sp_flag[linkVec[ii].second] = -2;
    bool f0 = (forwardFaNoSet.find(pair<int,int>(stolkVec[ii].first,linkVec[ii].first)) == forwardFaNoSet.end());
    bool f1 = (forwardFaNoSet.find(pair<int,int>(stolkVec[ii].second,linkVec[ii].second)) == forwardFaNoSet.end());
    bool b0 = (backwardFaNoSet.find(pair<int,int>(stolkVec[ii].first,linkVec[ii].second)) == backwardFaNoSet.end());
    bool b1 = (backwardFaNoSet.find(pair<int,int>(stolkVec[ii].second,linkVec[ii].first)) == backwardFaNoSet.end());
    if (f0 && f1 && b0 && b1) {
      // looks good...
      forwardFaNoSet.insert(pair<int,int>(stolkVec[ii].first,linkVec[ii].first));
      forwardFaNoSet.insert(pair<int,int>(stolkVec[ii].second,linkVec[ii].second));
      backwardFaNoSet.insert(pair<int,int>(stolkVec[ii].first,linkVec[ii].second));
      backwardFaNoSet.insert(pair<int,int>(stolkVec[ii].second,linkVec[ii].first));
    }
    else {
      // one or more failures. We must have introduced the edge in the wrong
      // orientation above (during the dp check). This should not happen unless the loop
      // logic above if flawed...
      cout << "Error: intersection edge orientation FAILURE ii: " << ii << " f0: " << f0 << " f1: " << f1 << " b0: " << b0 << " b1: " << b1 << endl;
      throw(-1);
    }
  }
  
  intersectionVec.clear();

  cout << " > links appear properly ordered" << endl;
  
  // now loop on the edges and flag them in or out...
  
  assert(stoedVec.size() == spoedVec.size());
  const int ned = stoedVec.size();
  int * ed_flag = new int[ned];
  for (int ied = 0; ied < ned; ++ied) {
    ed_flag[ied] = 0;
    if (sp_flag[spoedVec[ied].first] == -2) {
      if ( (forwardFaNoSet.find(pair<int,int>(stoedVec[ied].first,spoedVec[ied].first)) != forwardFaNoSet.end()) ||
           (backwardFaNoSet.find(pair<int,int>(stoedVec[ied].second,spoedVec[ied].first)) != backwardFaNoSet.end()) ) {
        ed_flag[ied] = -1;
        if (sp_flag[spoedVec[ied].second] == 0)
          sp_flag[spoedVec[ied].second] = -1;
        assert((sp_flag[spoedVec[ied].second] == -1)||(sp_flag[spoedVec[ied].second] == -2)); // -1 is out: could be -2 as well
      }
      else {
        ed_flag[ied] = 1;
        if (sp_flag[spoedVec[ied].second] == 0)
          sp_flag[spoedVec[ied].second] = 1;
        assert((sp_flag[spoedVec[ied].second] == 1)||(sp_flag[spoedVec[ied].second] == -2)); // 1 is in: could be -2 as well
      }
    }
    if (sp_flag[spoedVec[ied].second] == -2) {
      if ( (forwardFaNoSet.find(pair<int,int>(stoedVec[ied].second,spoedVec[ied].second)) != forwardFaNoSet.end()) ||
           (backwardFaNoSet.find(pair<int,int>(stoedVec[ied].first,spoedVec[ied].second)) != backwardFaNoSet.end()) ) {
        assert((ed_flag[ied] == 0)||(ed_flag[ied] == -1));
        ed_flag[ied] = -1;
        if (sp_flag[spoedVec[ied].first] == 0)
          sp_flag[spoedVec[ied].first] = -1;
        assert((sp_flag[spoedVec[ied].first] == -1)||(sp_flag[spoedVec[ied].first] == -2)); // -1 is out: could be -2 as well
      }
      else {
        assert((ed_flag[ied] == 0)||(ed_flag[ied] == 1));
        ed_flag[ied] = 1;
        if (sp_flag[spoedVec[ied].first] == 0)
          sp_flag[spoedVec[ied].first] = 1;
        assert((sp_flag[spoedVec[ied].first] == 1)||(sp_flag[spoedVec[ied].first] == -2)); // 1 is in: could be -2 as well
      }
    }
    // maek sure the edge was set...
    assert(ed_flag[ied] != 0);
  }

  cout << " > surface points around links flagged successfully" << endl;
  
  // now walk it out...
  
  /*
    {
    FILE * fp = fopen("selected_sp.dat","w");
    for (int ied = 0; ied < ned; ++ied) {
    if (ed_flag[ied] == 1) {
    {
    const int isp0 = spoedVec[ied].first;
    if (isp0 < nsp) {
    fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp0][0],xsp[isp0][1],xsp[isp0][2]);
    }
    else {
    fprintf(fp,"%18.15le %18.15le %18.15le\n",xspNewVec[(isp0-nsp)*3],xspNewVec[(isp0-nsp)*3+1],xspNewVec[(isp0-nsp)*3+2]);
    }
    }
    {
    const int isp1 = spoedVec[ied].second;
    if (isp1 < nsp) {
    fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp1][0],xsp[isp1][1],xsp[isp1][2]);
    }
    else {
    fprintf(fp,"%18.15le %18.15le %18.15le\n",xspNewVec[(isp1-nsp)*3],xspNewVec[(isp1-nsp)*3+1],xspNewVec[(isp1-nsp)*3+2]);
    }
    }
    }
    }
    fclose(fp);
    }
  */
  
  // now build half-edges (ist,(isp0,isp1))...
  
  vector<pair<int,pair<int,int> > > halfEdgeVec;
  
  // edges...

  for (int ied = 0; ied < ned; ++ied) {
    if (ed_flag[ied] == 1) {
      assert(spoedVec[ied].first != spoedVec[ied].second);
      halfEdgeVec.push_back(pair<int,pair<int,int> >(stoedVec[ied].first,pair<int,int>(spoedVec[ied].first,spoedVec[ied].second)));
      halfEdgeVec.push_back(pair<int,pair<int,int> >(stoedVec[ied].second,pair<int,int>(spoedVec[ied].second,spoedVec[ied].first)));
    }
  }
  
  delete[] ed_flag;
  
  // and the links...
  
  for (int ii = 0,nii=linkVec.size(); ii < nii; ++ii) {
    assert(linkVec[ii].first != linkVec[ii].second);
    halfEdgeVec.push_back(pair<int,pair<int,int> >(stolkVec[ii].first,pair<int,int>(linkVec[ii].first,linkVec[ii].second)));
    halfEdgeVec.push_back(pair<int,pair<int,int> >(stolkVec[ii].second,pair<int,int>(linkVec[ii].second,linkVec[ii].first)));
  }

  // now sort...
  
  sort(halfEdgeVec.begin(),halfEdgeVec.end());
  
  // and walk the faces: these are portions of tris...
  
  map<const int,int> spMap;
  set<pair<int,int> > frontSet;
  vector<NewTri> newTriVec;
  int ii = 0;
  while (ii != halfEdgeVec.size()) {
    assert(spMap.empty());
    assert(frontSet.empty());
    const int ist = halfEdgeVec[ii].first;
    int nsp_local = 0;
    spMap[halfEdgeVec[ii].second.first] = nsp_local++;
    spMap[halfEdgeVec[ii].second.second] = nsp_local++;
    frontSet.insert(pair<int,int>(0,1));
    ++ii;
    while ((ii != halfEdgeVec.size())&&(halfEdgeVec[ii].first == ist)) {
      int isp_local0;
      map<const int,int>::iterator it = spMap.find(halfEdgeVec[ii].second.first);
      if (it == spMap.end()) {
        spMap[halfEdgeVec[ii].second.first] = isp_local0 = nsp_local++;
      }
      else {
        isp_local0 = it->second;
      }
      int isp_local1;
      it = spMap.find(halfEdgeVec[ii].second.second);
      if (it == spMap.end()) {
        spMap[halfEdgeVec[ii].second.second] = isp_local1 = nsp_local++;
      }
      else {
        isp_local1 = it->second;
      }
      assert(frontSet.find(pair<int,int>(isp_local0,isp_local1)) == frontSet.end());
      frontSet.insert(pair<int,int>(isp_local0,isp_local1));
      ++ii;
    }
    // finally, check the edges of the tri to see if there is more to add...
    assert(st_flag[ist] != 0);
    FOR_I3 {
      if ((st_flag[ist]&(1<<i)) == 0) {
        // this edge was not added in the half-edges, so we need to add it 
        // if its nodes are flagged in. This means either a 1 in one or both...
        if ((sp_flag[spost[ist][i]] == 1)||(sp_flag[spost[ist][(i+1)%3]] == 1)) {
          // make sure neither point is out...
          assert((sp_flag[spost[ist][i]] != -1)&&(sp_flag[spost[ist][(i+1)%3]] != -1));
          // need to add this segment...
          int isp_local0;
          map<const int,int>::iterator it = spMap.find(spost[ist][i]);
          if (it == spMap.end()) {
            spMap[spost[ist][i]] = isp_local0 = nsp_local++;
          }
          else {
            isp_local0 = it->second;
          }
          int isp_local1;
          it = spMap.find(spost[ist][(i+1)%3]);
          if (it == spMap.end()) {
            spMap[spost[ist][(i+1)%3]] = isp_local1 = nsp_local++;
          }
          else {
            isp_local1 = it->second;
          }
          assert(frontSet.find(pair<int,int>(isp_local0,isp_local1)) == frontSet.end());
          frontSet.insert(pair<int,int>(isp_local0,isp_local1));
        }
      }
    }
    // we now have a front that should be continuous. We also have a local set of nodes that need to 
    // be projected into 2d...
    int * sposp_local = new int[nsp_local];
    // first check that the front touches each node twice...
    for (int isp_local = 0; isp_local < nsp_local; ++isp_local) 
      sposp_local[isp_local] = 0;
    for (set<pair<int,int> >::iterator it = frontSet.begin(); it != frontSet.end(); ++it) {
      --sposp_local[it->first];
      --sposp_local[it->second];
    }
    for (int isp_local = 0; isp_local < nsp_local; ++isp_local) 
      assert(sposp_local[isp_local] == -2);
    const double normal[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    double e1[3],e2[3];
    MiscUtils::getBestE1E2FromE0(e1,e2,normal);
    // and scale e1 and e2 by 1/sqrt(area) to give O(1) 2d coordinates...
    const double area = MAG(normal); assert(area > 0.0);
    const double factor = 1.0/sqrt(area); assert(factor > 0.0);
    FOR_I3 {
      e1[i] *= factor;
      e2[i] *= factor;
    }
    double (*xsp_local_2d)[2] = new double[nsp_local][2];
    for (map<const int,int>::iterator it = spMap.begin(); it != spMap.end(); ++it) {
      const int isp = it->first;
      const int isp_local = it->second; assert((isp_local >= 0)&&(isp_local < nsp_local)); 
      assert(sposp_local[isp_local] == -2); // from check above
      sposp_local[isp_local] = isp;
      double dx[3];
      if (isp < nsp) {
        FOR_I3 dx[i] = xsp[isp][i] - xsp[spost[ist][0]][i]; // use spost[ist][0] to zero everything
      }
      else {
        FOR_I3 dx[i] = xspNewVec[(isp-nsp)*3+i] - xsp[spost[ist][0]][i];
      }
      xsp_local_2d[isp_local][0] = DOT_PRODUCT(dx,e1);
      xsp_local_2d[isp_local][1] = DOT_PRODUCT(dx,e2);
    }
    const int newTriVec_size0 = newTriVec.size();
    const int ierr = GeomUtils::addTrisMarchingFront2d(newTriVec,frontSet,xsp_local_2d,nsp_local);
    if (ierr != 0) {
      cout << "addTrisMarchingFront2d failed." << endl;
      throw(0);
    }
    // for the added tris, we need to set the spost back to the 3d node indexing,
    // and set the izn,isz...
    for (int itri = newTriVec_size0,ntri=newTriVec.size(); itri < ntri; ++itri) {
      FOR_I3 {
        const int isp_local = newTriVec[itri].spost[i];
        assert((isp_local >= 0)&&(isp_local < nsp_local));
        newTriVec[itri].spost[i] = sposp_local[isp_local];
      }
      newTriVec[itri].znost = znost[ist];
      newTriVec[itri].szost = szost[ist];
    }
    // cleanup...
    delete[] xsp_local_2d;
    delete[] sposp_local;
    spMap.clear();
    frontSet.clear();
  }

  // ok - we have new tris...
  
  cout << " > done creating new tris by marching front: newTriVec.size(): " << newTriVec.size() << endl;

  int8 * stosp_i = new int8[nsp+1];
  for (int isp = 0; isp < nsp; ++isp)
    stosp_i[isp+1] = 0;
  
  stack<int> stStack;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] == 0) {
      if ((sp_flag[spost[ist][0]] == 1)||(sp_flag[spost[ist][1]] == 1)||(sp_flag[spost[ist][2]] == 1)) {
        assert((sp_flag[spost[ist][0]] != -1)&&(sp_flag[spost[ist][1]] != -1)&&(sp_flag[spost[ist][2]] != -1));
        // this tri is in, so mark as 1...
        st_flag[ist] = 1;
      }
      else if ((sp_flag[spost[ist][0]] == -1)||(sp_flag[spost[ist][1]] == -1)||(sp_flag[spost[ist][2]] == -1)) {
        assert((sp_flag[spost[ist][0]] != 1)&&(sp_flag[spost[ist][1]] != 1)&&(sp_flag[spost[ist][2]] != 1));
        // this tri is out, so add to stack...
        st_flag[ist] = -1;
        stStack.push(ist);
      }
      else {
        // this tri may be in or out, so leave as zero...
        FOR_I3 {
          const int isp = spost[ist][i];
          ++stosp_i[isp+1];
        }
      }        
    }
    else {
      // this tri is out. Its edges have been replaced with new triangles...
      st_flag[ist] = -1;
    }
  }
  
  // build csr structure...
  
  stosp_i[0] = 0;
  for (int isp = 0; isp < nsp; ++isp)
    stosp_i[isp+1] += stosp_i[isp];
  const int8 stosp_s = stosp_i[nsp];
  int * stosp_v = new int[stosp_s];
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] == 0) {
      FOR_I3 {
        const int isp = spost[ist][i];
        stosp_v[stosp_i[isp]++] = ist;
      }
    }
  }
  
  // and rewind...

  for (int isp = nsp-1; isp > 0; --isp)
    stosp_i[isp] = stosp_i[isp-1];
  stosp_i[0] = 0;
  
  // and march through the tris...
  
  while (!stStack.empty()) {
    const int ist = stStack.top(); stStack.pop();
    assert(st_flag[ist] == -1);
    // use stosp_i/v to grab our node-based nbrs...
    FOR_I3 {
      const int isp = spost[ist][i];
      for (int sos = stosp_i[isp]; sos != stosp_i[isp+1]; ++sos) {
        const int ist_nbr = stosp_v[sos];
        if (st_flag[ist_nbr] == 0) {
          st_flag[ist_nbr] = -1;
          stStack.push(ist_nbr);
        }
      }
    }
  }

  delete[] stosp_i;
  delete[] stosp_v;
  
  // now all tris are marked with -1 are out, 0 and 1 are in...

  sp_flag.setAll(-1);
  int nst_final = 0;
  int nsp_final = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] != -1) {
      assert((st_flag[ist] == 0)||(st_flag[ist] == 1));
      ++nst_final;
      FOR_I3 {
        const int isp = spost[ist][i];
        assert((isp >= 0)&&(isp < nsp_new));
        if (sp_flag[isp] == -1) {
          ++nsp_final;
          sp_flag[isp] = 0;
        }
      }
    }
  }
  nst_final += newTriVec.size();
  for (int ii = 0; ii < newTriVec.size(); ++ii) {
    FOR_I3 {
      const int isp = newTriVec[ii].spost[i];
      assert((isp >= 0)&&(isp < nsp_new));
      if (sp_flag[isp] == -1) {
        ++nsp_final;
        sp_flag[isp] = 0;
      }
    }
  }

  cout << " > final nst,nsp: " << nst_final << " " << nsp_final << endl;
  
  // index sp_flag and copy down xsp...
  double (*xsp_final)[3] = new double[nsp_final][3];
  int isp_final = 0;
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] == 0) {
      sp_flag[isp] = isp_final;
      FOR_I3 xsp_final[isp_final][i] = xsp[isp][i];
      ++isp_final;
    }
  }
  for (int isp = nsp; isp < nsp_new; ++isp) {
    assert(sp_flag[isp] == 0);
    sp_flag[isp] = isp_final;
    FOR_I3 xsp_final[isp_final][i] = xspNewVec[(isp-nsp)*3+i];
    ++isp_final;
  }
  assert(isp_final == nsp_final);
  nsp = nsp_final;
  delete[] xsp;
  xsp = xsp_final;

  // and now the tris...
  zone_flag.resize(zoneVec.size());
  zone_flag.setAll(0);
  sz_flag.resize(nsz);
  sz_flag.setAll(0);

  int (*spost_final)[3] = new int[nst_final][3];
  int *znost_final = new int[nst_final];
  int *szost_final = new int[nst_final];
  int ist_final = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] != -1) {
      assert((st_flag[ist] == 0)||(st_flag[ist] == 1));
      FOR_I3 {
        const int isp = spost[ist][i];
        spost_final[ist_final][i] = sp_flag[isp];
      }
      const int izn = znost[ist]; assert((izn >= 0)&&(izn < int(zoneVec.size())));
      const int isz = szost[ist]; assert((isz >= 0)&&(isz < nsz));
      znost_final[ist_final] = izn;
      szost_final[ist_final] = isz;
      // count zones and subzones...
      ++zone_flag[izn];
      ++sz_flag[isz];
      ++ist_final;
    }
  }
  for (int ii = 0; ii < newTriVec.size(); ++ii) {
    FOR_I3 {
      const int isp = newTriVec[ii].spost[i];
      assert((isp >= 0)&&(isp < nsp_new));
      spost_final[ist_final][i] = sp_flag[isp];
    }
    const int izn = newTriVec[ii].znost; assert((izn >= 0)&&(izn < int(zoneVec.size())));
    const int isz = newTriVec[ii].szost; assert((isz >= 0)&&(isz < nsz));
    znost_final[ist_final] = izn;
    szost_final[ist_final] = isz;
    // count zones and subzones...
    ++zone_flag[izn];
    ++sz_flag[isz];
    ++ist_final;
  }
  assert(ist_final == nst_final);

  delete[] spost; spost = spost_final;
  delete[] znost; znost = znost_final;
  szost.resize(nst_final);
  for (int ist = 0; ist < nst_final; ++ist)
    szost[ist] = szost_final[ist];
  delete[] szost_final;
  nst = nst_final;

  // at this point we have counted the zone_flag and sz_flag tris. We should
  // clear out any zones that have changed...
  // ME: what is the way to do this?

  for (int izn = 0,nzn=zoneVec.size(); izn < nzn; ++izn)
    if (zone_flag[izn] == 0)
      cout << "Warning: zone " << izn << " reduced to 0 tris" << endl;

  for (int isz = 0; isz < nsz; ++isz)
    if (sz_flag[isz] == 0)
      cout << "Warning: subzone " << isz << " reduced to 0 tris" << endl;
  
  // I think this is what has to be run to cleanup...
  // taken from SimpleSurface_imprint...
  
  // new zones have local-indexed subzones from source zones
  // to rebuild properly pruned szozn_i and szost we need to run these
  buildSzoznFromLocalSzost();
  localSzostToGlobal();
  pruneEmptyZonesAndSubzones();
  clearSubzoneData();
  clearZoneData();
  clearTeost();
  clearDynamicEdgeGroups();
  clearStoszSzosz();
  clearOpenEdgeGroups();
  clearSposp();
  clearNonManifoldData();  // can remove once multiedges are handled

  WUI(INFO,"INTERSECT successful");

}

void SimpleSurface::processTriTriIntersection(map<const pair<int,int>,int>& edgeMap,vector<pair<int,int> >& stoedVec,vector<IntersectionData>& intersectionVec,
                                              const int * const idata,const double * const ddata,const int ist0,const int ist1) {
  
  switch (idata[0]) {
  case MiscUtils::EDGE_TRI_INT:
    {
      // edge on tri0 is in idata[1]...
      // look for the reverse edge first: most common...
      int ied0;
      double wgt0;
      map<const pair<int,int>,int>::iterator it = edgeMap.find(pair<int,int>(spost[ist0][(idata[1]+1)%3],spost[ist0][idata[1]]));
      if (it != edgeMap.end()) {
        ied0 = it->second;
        wgt0 = 1.0-ddata[0];
        if (stoedVec[ied0].second == -1)
          stoedVec[ied0].second = ist0;
        assert(stoedVec[ied0].second == ist0);
      }
      else {
        // look in forward direction...
        it = edgeMap.find(pair<int,int>(spost[ist0][idata[1]],spost[ist0][(idata[1]+1)%3]));
        if (it != edgeMap.end()) {
          ied0 = it->second;
          wgt0 = ddata[0];
          assert(stoedVec[ied0].first == ist0);
        }
        else {
          // add in forward direction...
          edgeMap[pair<int,int>(spost[ist0][idata[1]],spost[ist0][(idata[1]+1)%3])] = ied0 = stoedVec.size();
          stoedVec.push_back(pair<int,int>(ist0,-1));
          wgt0 = ddata[0];
        }
      }
      intersectionVec.push_back(IntersectionData(EDGE_FACE_INTERSECTION,ied0,ist1,wgt0));
    }
    break;
  case MiscUtils::TRI_EDGE_INT:
    {
      // edge on tri1 is in idata[1]...
      // look for the reverse edge first: most common...
      int ied1;
      double wgt1;
      map<const pair<int,int>,int>::iterator it = edgeMap.find(pair<int,int>(spost[ist1][(idata[1]+1)%3],spost[ist1][idata[1]]));
      if (it != edgeMap.end()) {
        ied1 = it->second;
        wgt1 = 1.0-ddata[0];
        if (stoedVec[ied1].second == -1)
          stoedVec[ied1].second = ist1;
        assert(stoedVec[ied1].second == ist1);
      }
      else {
        // look in forward direction...
        it = edgeMap.find(pair<int,int>(spost[ist1][idata[1]],spost[ist1][(idata[1]+1)%3]));
        if (it != edgeMap.end()) {
          ied1 = it->second;
          wgt1 = ddata[0];
          assert(stoedVec[ied1].first == ist1);
        }
        else {
          // add in forward direction...
          edgeMap[pair<int,int>(spost[ist1][idata[1]],spost[ist1][(idata[1]+1)%3])] = ied1 = stoedVec.size();
          stoedVec.push_back(pair<int,int>(ist1,-1));
          wgt1 = ddata[0];
        }
      }
      intersectionVec.push_back(IntersectionData(EDGE_FACE_INTERSECTION,ied1,ist0,wgt1));
    }
    break;
  default:
    // we will hit this eventually: just add new intersections...
    assert(0);
  }
  
}

