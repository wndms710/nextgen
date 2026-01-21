#include "SimpleSurface.hpp"
#include "Adt.hpp"
#include "Histogram.hpp"

void SimpleSurface::zipOpenEdges(const vector<int>& group_indices,const double edge_factor,const double delta_max) {

  // need open edge groups and (implicitly) teost
  ensureOpenEdgeGroups();

  // zip routines require a selection of open edges for zipping...

  vector<uint8> zipTeostVec;
  const int ngr = getNOpenEdgeGroups();
  bool *gr_flag = new bool[ngr];
  for (int igr = 0; igr < ngr; ++igr) gr_flag[igr] = false;

  for (vector<int>::const_iterator it=group_indices.begin(); it!=group_indices.end(); ++it) {
    if ((*it < 0) || (*it >= ngr)) {
      WUI(WARN,"invalid open edge group " << *it << " selected; skipping");
    }
    else {
      gr_flag[*it] = true;
    }
  }

  for (map<uint,int>::const_iterator it=oe_to_group.begin(); it!=oe_to_group.end(); ++it) {
    if (gr_flag[it->second]) {
      uint ist; uchar i;
      unpackEdge(ist,i,it->first);
      zipTeostVec.push_back( uint8(i)<<60 | ist ); // edge + triangle
    }
  }

  delete[] gr_flag;

  // then zip...
  zipOpenEdges(zipTeostVec,edge_factor,delta_max);
}

void SimpleSurface::zipOpenEdges(const double edge_factor,const double delta_max) {
  // tries to zip everything...

  ensureTeost();
  // grab ALL open edges...
  vector<uint8> zipTeostVec;

  for (map<uint,int>::const_iterator it=oe_to_group.begin(); it!=oe_to_group.end(); ++it) {
    uint ist; uchar i;
    unpackEdge(ist,i,it->first);
    zipTeostVec.push_back( (uint8(i)<<60) | ist ); // edge + triangle
  }

  if (zipTeostVec.empty())
    return;

  // pass these edges into the zip routine...
  zipOpenEdges(zipTeostVec,edge_factor,delta_max);
}

void SimpleSurface::zipOpenEdges(vector<uint8>& zipTeostVec,const double edge_factor,const double delta_max) {

  vector<NodeNodeMerge> nodeNodeMergeVec;
  vector<EdgeNodeMerge> edgeNodeMergeVec;
  setZipMergeVecs(nodeNodeMergeVec,edgeNodeMergeVec,zipTeostVec,edge_factor,delta_max);

  cout << " > full nodeNodeMergeVec.size(): " << nodeNodeMergeVec.size() << " edgeNodeMergeVec.size(): " << edgeNodeMergeVec.size() << endl;

  /*

  // build the histogram of tolerances...

  {

    Histogram hist(1.0E-8*getBoundingBoxRmax(),0.01*getBoundingBoxRmax(),MPI_COMM_NULL);
    hist.setLog(true);

    const int nnn = nodeNodeMergeVec.size();
    const int nen = edgeNodeMergeVec.size();
    double * tol = new double[nnn+nen];
    for (int ii = 0, limit = nodeNodeMergeVec.size(); ii < limit; ++ii) {
      tol[ii] = sqrt(nodeNodeMergeVec[ii].d2);
    }
    for (int ii = 0, limit = edgeNodeMergeVec.size(); ii < limit; ++ii) {
      tol[nnn+ii] = sqrt(edgeNodeMergeVec[ii].d2);
    }
    MiscUtils::dumpRange(tol,nnn+nen,"d");
    hist.add(tol,nnn+nen);
    hist.write("hist.dat");
    delete[] tol;

    cout << "take a look at hist" << endl;
    getchar();

  }

  // remove everything above tol...

  const double tol = 1.0E-4;

  cout << " > removing merges above tol: " << tol << endl;

  int new_size = 0;
  for (int ii = 0, limit = nodeNodeMergeVec.size(); ii < limit; ++ii) {
    if (nodeNodeMergeVec[ii].d2 < tol*tol) {
      if (new_size != ii)
	nodeNodeMergeVec[new_size] = nodeNodeMergeVec[ii];
      ++new_size;
    }
  }
  nodeNodeMergeVec.resize(new_size);

  new_size = 0;
  for (int ii = 0, limit = edgeNodeMergeVec.size(); ii < limit; ++ii) {
    if (edgeNodeMergeVec[ii].d2 < tol*tol) {
      if (new_size != ii)
	edgeNodeMergeVec[new_size] = edgeNodeMergeVec[ii];
      ++new_size;
    }
  }
  edgeNodeMergeVec.resize(new_size);

  cout << " > after tol nodeNodeMergeVec.size(): " << nodeNodeMergeVec.size() << " edgeNodeMergeVec.size(): " << edgeNodeMergeVec.size() << endl;

  */

  // convert any edgeNodeMergeVec's to nodeNodeMergeVecs...

  if (!edgeNodeMergeVec.empty())
    splitEdgesAndAddNodeNodeMerges(nodeNodeMergeVec,edgeNodeMergeVec,zipTeostVec);

  // finally merge the nodes...

  if (!nodeNodeMergeVec.empty())
    mergeNodes(nodeNodeMergeVec);

}

int SimpleSurface::setZipMergeVecs(vector<NodeNodeMerge>& nodeNodeMergeVec,vector<EdgeNodeMerge>& edgeNodeMergeVec,vector<uint8>& zipTeostVec,const double edge_factor,const double delta_max) {

  int isp_debug = -1;

  /*
  double xsp_debug[3] = { 18.138, 0.846643, 10.6514 };
  double d2_debug;
  for (int isp = 0; isp < nsp; ++isp) {
    const double d2 = DIST2(xsp_debug,xsp[isp]);
    if ((isp_debug == -1)||(d2 < d2_debug)) {
      isp_debug = isp;
      d2_debug = d2;
    }
  }
  cout << "got isp_debug: " << isp_debug << " at: " << COUT_VEC(xsp[isp_debug]) << " d: " << sqrt(d2_debug) << endl;
  getchar();
  */

  // local versions...

  //vector<NodeNodeMerge> nodeNodeMergeVec;
  //vector<EdgeNodeMerge> edgeNodeMergeVec;

  // no need to return int for this routine. No error possible?

  // use bits to count edge touches at each node...
  sp_flag.resize(nsp);
  sp_flag.setAll(0);

  const int ne = zipTeostVec.size();
  for (int ie = 0; ie < ne; ++ie) {
    const int ist = int( zipTeostVec[ie] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
    const int i = int( zipTeostVec[ie]>>60 ); assert((i >= 0)&&(i < 3));
    sp_flag[spost[ist][i]] += 1;
    sp_flag[spost[ist][(i+1)%3]] += 1;
  }

  // find the participating nodes...

  int nsp_ = 0;
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] > 0) {

      if (isp == isp_debug) cout << "DEBUG: sp_flag[isp]: " << sp_flag[isp] << endl;

      nsp_ += 1;
    }
  }

  cout << " > nsp_: " << nsp_ << endl;

  // build a CSR structure at each node that includes its associated edges...

  int *sposp_ = new int[nsp_];
  nsp_ = 0;
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] > 0) {

      sposp_[nsp_] = isp;
      sp_flag[isp] = nsp_++;
    }
    else {
      sp_flag[isp] = -1;
    }
  }

  double *delta_ = new double[nsp_];
  for (int isp_ = 0; isp_ < nsp_; ++isp_) {
    delta_[isp_] = 1.0E+20;
  }

  // delta is a fraction of the minimum edge length for ALL
  // edges that touch a node. This attempts to prevent inadvertently grabbing
  // a node that we shouldn't...

  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      const int isp_ = sp_flag[spost[ist][i]];
      if (isp_ >= 0) {
	// store d2 for now...
	delta_[isp_] = min(delta_[isp_],
			   min(DIST2(xsp[spost[ist][i]],xsp[spost[ist][(i+1)%3]]),
			       DIST2(xsp[spost[ist][i]],xsp[spost[ist][(i+2)%3]])));
      }
    }
  }

  // then take sqrt and multiply by factor...
  for (int isp_ = 0; isp_ < nsp_; ++isp_) {
    delta_[isp_] = min(delta_max,edge_factor*sqrt(delta_[isp_])); // TODO: this factor is important -- there is a value that will prevent folding...
  }
  MiscUtils::dumpRange(delta_,nsp_,"delta_ range");

  double (*bbmin)[3] = new double[ne][3];
  double (*bbmax)[3] = new double[ne][3];

  for (int ie = 0; ie < ne; ++ie) {
    const int ist = int( zipTeostVec[ie] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
    const int i = int( zipTeostVec[ie]>>60 ); assert((i >= 0)&&(i < 3));
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    FOR_J3 bbmin[ie][j] = min(xsp[isp0][j],xsp[isp1][j]);
    FOR_J3 bbmax[ie][j] = max(xsp[isp0][j],xsp[isp1][j]);
  }

  Adt<double> eAdt(ne,bbmin,bbmax);

  delete[] bbmin;
  delete[] bbmax;



  // now loop on points...

  vector<int> candidateVec;
  for (int isp_ = 0; isp_ < nsp_; ++isp_) {

    const int isp = sposp_[isp_]; assert((isp >= 0)&&(isp < nsp));

    assert(candidateVec.empty());
    eAdt.buildListForSphere(candidateVec,xsp[isp],delta_[isp_]);
    if (!candidateVec.empty()) {

      int isp0_min,isp1_min,ie_min = -1;
      double t_min;
      double d2_min = delta_[isp_]*delta_[isp_]; // start at delta^2
      for (int ii = 0,ii_end=candidateVec.size(); ii < ii_end; ++ii) {
	const int ie = candidateVec[ii];
	const int ist = int( zipTeostVec[ie] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
	const int i = int( zipTeostVec[ie]>>60 ); assert((i >= 0)&&(i < 3));
	const int isp0 = spost[ist][i];
	const int isp1 = spost[ist][(i+1)%3];
	if ((isp0 == isp)||(isp1 == isp))
	  continue;
	const double dx[3] = DIFF(xsp[isp1],xsp[isp0]);
	const double dxp[3] = DIFF(xsp[isp],xsp[isp0]);
	const double dp = DOT_PRODUCT(dxp,dx);
	if (dp <= 0.0) {
	  // we are closest to the first point. Note this includes the case when dx is zero, and/or dxp is zero...
	  const double d2 = dxp[0]*dxp[0] + dxp[1]*dxp[1] + dxp[2]*dxp[2];
	  //cout << " > case 1: got d2: " << d2 << endl;
	  if (d2 < d2_min) {
	    d2_min   = d2;
	    isp0_min = isp0;
	    isp1_min = isp1;
	    ie_min   = ie;
	    t_min    = 0.0;
	  }
	}
	else {
	  const double dx2 = DOT_PRODUCT(dx,dx);
	  assert(dx2 > 0.0);
	  if (dp >= dx2) {
	    // we are closest to the second point...
	    const double d2 =
	      (xsp[isp][0]-xsp[isp1][0])*(xsp[isp][0]-xsp[isp1][0]) +
	      (xsp[isp][1]-xsp[isp1][1])*(xsp[isp][1]-xsp[isp1][1]) +
	      (xsp[isp][2]-xsp[isp1][2])*(xsp[isp][2]-xsp[isp1][2]);
	    if (d2 < d2_min) {
	      d2_min   = d2;
	      isp0_min = isp0;
	      isp1_min = isp1;
	      ie_min   = ie;
	      t_min    = 1.0;
	    }
	  }
	  else {
	    // we are closest to a point along the edge...
	    const double this_t = max(0.0,min(1.0,dp/dx2));
	    const double d2 = max(0.0,dxp[0]*dxp[0] + dxp[1]*dxp[1] + dxp[2]*dxp[2] - dp*dp/dx2);
	    if (d2 < d2_min) {
	      d2_min   = d2;
	      isp0_min = isp0;
	      isp1_min = isp1;
	      ie_min   = ie;
	      t_min    = this_t;
	    }
	  }
	}
      }
      candidateVec.clear();

      if (ie_min != -1) {

	// got one!

	if (t_min == 0.0) {
	  nodeNodeMergeVec.push_back(NodeNodeMerge(isp,isp0_min,d2_min));

	  if (isp == isp_debug) cout << "DEBUG Adding node-node merge 0 for isp_debug" << endl;
	  if (isp0_min == isp_debug) cout << "DEBUG Adding node-node merge 1 for isp_debug" << endl;

	}
	else if (t_min == 1.0) {
	  nodeNodeMergeVec.push_back(NodeNodeMerge(isp,isp1_min,d2_min));

	  if (isp == isp_debug) cout << "DEBUG Adding node-node merge 0 for isp_debug" << endl;
	  if (isp1_min == isp_debug) cout << "DEBUG Adding node-node merge 1 for isp_debug" << endl;

	}
	else {

	  // --------------------------------------------------------------------
	  // this is a node-edge intersection. Use the nearby geometry to decide
	  // if we can make it a node-node...
	  // --------------------------------------------------------------------
	  bool b_node_node = false;
	  if (t_min >= 0.5) {
	    const double d2_isp1 = DIST2(xsp[isp],xsp[isp1_min]);
	    if (d2_isp1 < 2.0*d2_min) {
	      // this is the 45 degree triangle rule. We are about to move the point
	      // a minimum distance of sqrt(d2_min) -- really half that because of the
	      // node averaging that will occur -- so if a corner is within sqrt(2)
	      // of that distance, we take the corner. This is the "gelato" rule ;)
	      nodeNodeMergeVec.push_back(NodeNodeMerge(isp,isp1_min,d2_isp1));
	      b_node_node = true;
	    }
	  }
	  else {
	    const double d2_isp0 = DIST2(xsp[isp],xsp[isp0_min]);
	    if (d2_isp0 < 2.0*d2_min) {
	      // see note above...
	      nodeNodeMergeVec.push_back(NodeNodeMerge(isp,isp0_min,d2_isp0));
	      b_node_node = true;
	    }
	  }

	  // if we didn't create a node-node, then do the original edge-node...
	  if (!b_node_node) {
	    edgeNodeMergeVec.push_back(EdgeNodeMerge(ie_min,isp,t_min,d2_min));
	  }
	}
      }
    }
  }

  // = new

  delete[] sposp_;
  delete[] delta_;

  return 0;

}

void SimpleSurface::splitEdgesAndAddNodeNodeMerges(vector<NodeNodeMerge>& nodeNodeMergeVec,vector<EdgeNodeMerge>& edgeNodeMergeVec,vector<uint8>& zipTeostVec) {

  assert(!edgeNodeMergeVec.empty());

  st_flag.setLength(nst);
  st_flag.setAll(0);

  sort(edgeNodeMergeVec.begin(),edgeNodeMergeVec.end());
  int ie_current = -1;
  for (int ii = 0, limit = edgeNodeMergeVec.size(); ii < limit; ++ii) {
    if (edgeNodeMergeVec[ii].ioe != ie_current) {
      ie_current = edgeNodeMergeVec[ii].ioe;
      const int ist = int( zipTeostVec[ie_current] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
      st_flag[ist] += 1;
    }
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
  nst = nst_old + nst_split_tri_new + edgeNodeMergeVec.size();
  const int nsp_old = nsp;
  nsp = nsp_old + nsp_split_tri_new + edgeNodeMergeVec.size();

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

  // set periodicity AFTER zippering...
  assert(pbi == NULL);

  // build new tris and nodes...

  int ist_new = nst_old;
  int isp_new = nsp_old;

  // start by splitting the tris with double edges...

  for (int ist = 0; ist < nst_old; ++ist) {
    //cout << "WORKING ON ist: " << ist << endl;
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

  // also a new node and 1 new tris are required for each edge...

  ie_current = -1;
  double x0[3],x1[3];
  for (int ii = 0, limit = edgeNodeMergeVec.size(); ii < limit; ++ii) {
    const int ie = edgeNodeMergeVec[ii].ioe;
    int ist = int( zipTeostVec[ie] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
    const int i = int( zipTeostVec[ie]>>60 ); assert((i >= 0)&&(i < 3));
    // if ist_orig is a tri that has been duplicated, then adjust ist to the duplicate. Note
    // that because of the way we oriented the new tris, there is no need to change the
    // edge access pattern. Just modify the ist to the appropriate new tri...
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
    FOR_J3 xsp[isp_new][j] = (1.0-edgeNodeMergeVec[ii].t)*x0[j] + edgeNodeMergeVec[ii].t*x1[j];
    nodeNodeMergeVec.push_back(NodeNodeMerge(isp_new,edgeNodeMergeVec[ii].isp,edgeNodeMergeVec[ii].d2));
    /*
    {
      // check proximity to the point...
      double xsp_t[3]; FOR_J3 xsp_t[j] = xsp[edgeNodeMergeVec[ii].isp][j];
      cout << " SECOND d2: " << DIST2(xsp_t,xsp[isp_new]) << " edgeNodeMergeVec[ii].d2: " << edgeNodeMergeVec[ii].d2 << endl;
      getchar();
    }
    */
    // the new tri goes in where the old tri was...
    spost[ist_new][i]       = spost[ist][i];
    spost[ist_new][(i+1)%3] = isp_new;
    spost[ist_new][(i+2)%3] = spost[ist][(i+2)%3]; // same internal node -- mod if tri split?
    // the new node becomes isp0 of the original tri...
    spost[ist][i] = isp_new; // do this AFTER above
    // new zone is the same zone...
    znost[ist_new] = znost[ist];
    // new subzone is the same...
    szost[ist_new] = szost[ist];
    // increment...
    ++ist_new;
    ++isp_new;
  }
  edgeNodeMergeVec.clear();

  assert(isp_new == nsp);
  assert(ist_new == nst);

  // because of new nodes and tris...

  clearTeost();
  clearDynamicEdgeGroups();
  clearNonManifoldData();

}

void SimpleSurface::mergeNodes(vector<NodeNodeMerge>& nodeNodeMergeVec) {

  cout << " > mergeNodes" << endl;

  // node merges should clear Teost
  clearTeost();
  clearDynamicEdgeGroups();
  clearNonManifoldData();

  // index into sp_flag...
  sp_flag.setLength(nsp);
  for (int isp = 0; isp < nsp; ++isp)
    sp_flag[isp] = isp;

  // now push down connections...
  for (int ii = 0, limit = nodeNodeMergeVec.size(); ii < limit; ++ii) {
    const int isp0 = nodeNodeMergeVec[ii].isp0;
    const int isp1 = nodeNodeMergeVec[ii].isp1;
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
      FOR_I3 xsp2[isp_new][i] = xsp[isp][i]*xsp[isp][i];
      FOR_I3 xsp[isp_new][i] = xsp[isp][i];
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
      FOR_I3 xsp2[isp_new][i] += xsp[isp][i]*xsp[isp][i];
      FOR_I3 xsp[isp_new][i] += xsp[isp][i];
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
  }

  cout << " > max rms point motion: absolute: " << rms_max << ", normalized: " << rms_max/getBoundingBoxRmax() << endl;

  delete[] xsp2;
  delete[] nsp_merged;

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
  if (collapsed_tri_count > 0)
    cout << " > node merging created " << collapsed_tri_count << " collapsed tris" << endl;

}



















class PossibleMerge {
public:
  double d2,t;
  int ioe,isp0,isp1;
  PossibleMerge(const double d2_,const double t_,const int ioe_,const int isp0_,const int isp1_) {
    d2 = d2_;
    t = t_;
    ioe = ioe_;
    isp0 = isp0_;
    isp1 = isp1_;
  }
};

void SimpleSurface::buildZipMergeVecs(vector<pair<int,int> >& nodeNodeMergeVec,vector<pair<pair<int,double>,int> >& edgeNodeMergeVec,
				   const vector<uint8>& zipTeostVec) {

  assert(0);

  // count edge touches at each node...
  sp_flag.resize(nsp);
  sp_flag.setAll(0);
  const int noe = zipTeostVec.size();
  for (int ioe = 0; ioe < noe; ++ioe) {
    const int ist = int( zipTeostVec[ioe] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
    const int i = int( zipTeostVec[ioe]>>60 ); assert((i >= 0)&&(i < 3));
    assert(isEdgeOpen(ist,i));
    sp_flag[spost[ist][i]] += 1;
    sp_flag[spost[ist][(i+1)%3]] += 1;
  }

  // only nodes with 1 or 2 touches should be zip'd...
  // here we use an underscore suffix to indicate these reduced counts...
  int nsp_ = 0;
  int nsp_many_ = 0;
  for (int isp = 0; isp < nsp; ++isp) {
    if ((sp_flag[isp] >= 1)&&(sp_flag[isp] <= 2)) {
      sp_flag[isp] = nsp_++;
    }
    else if (sp_flag[isp] == 0) {
      sp_flag[isp] = -1; // no touches -- don't touch this node...
    }
    else {
      ++nsp_many_;
      sp_flag[isp] = -2; // this is a node touching many open edges -- we do not explicitly reconnect, but may want to know...
    }
  }

  cout << " > zip node counts: nsp_: " << nsp_ << ", nsp_many_: " << nsp_many_ << endl;

  int (*sp_nbr_of_sp_)[2] = new int[nsp_][2];
  for (int isp_ = 0; isp_ < nsp_; ++isp_) FOR_I2 sp_nbr_of_sp_[isp_][i] = -1;

  int *sposp_ = new int[nsp_];
  for (int isp_ = 0; isp_ < nsp_; ++isp_) sposp_[isp_] = -1;

  double (*bbmin)[3] = new double[noe][3];
  double (*bbmax)[3] = new double[noe][3];

  for (int ioe = 0; ioe < noe; ++ioe) {
    const int ist = int( zipTeostVec[ioe] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
    const int i = int( zipTeostVec[ioe]>>60 ); assert((i >= 0)&&(i < 3));
    assert(isEdgeOpen(ist,i));
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    // set the bounding box for the edge...
    FOR_I3 bbmin[ioe][i] = min(xsp[isp0][i],xsp[isp1][i]);
    FOR_I3 bbmax[ioe][i] = max(xsp[isp0][i],xsp[isp1][i]);
    const int isp0_ = sp_flag[isp0]; assert((isp0_ >= 0)||(isp0_ == -2));
    if (isp0_ >= 0) {
      assert(sp_nbr_of_sp_[isp0_][1] == -1);
      sp_nbr_of_sp_[isp0_][1] = isp1;
      if (sposp_[isp0_] == -1)
	sposp_[isp0_] = isp0;
      assert(sposp_[isp0_] == isp0);
    }
    const int isp1_ = sp_flag[isp1]; assert((isp1_ >= 0)||(isp1_ == -2));
    if (isp1_ >= 0) {
      assert(sp_nbr_of_sp_[isp1_][0] == -1);
      sp_nbr_of_sp_[isp1_][0] = isp0;
      if (sposp_[isp1_] == -1)
	sposp_[isp1_] = isp1;
      assert(sposp_[isp1_] == isp1);
    }
  }

  // check...
  for (int isp_ = 0; isp_ < nsp_; ++isp_) assert(sposp_[isp_] >= 0);

  // now sp_nbr_of_sp_ should be set with atleast one isp nbr for each isp_...
  // and we have the bbox min and max for each open edge oe...

  Adt<double> oeAdt(noe,bbmin,bbmax);
  delete[] bbmin; bbmin = NULL;
  delete[] bbmax; bbmax = NULL;

  // and finally we need a length scale at each node. Use the MAX of ANY touching
  // edge (not just open edges). This will result in more hits and take longer, but
  // more stuff will get zipped...

  double *delta_sp_ = new double[nsp_];
  for (int isp_ = 0; isp_ < nsp_; ++isp_) delta_sp_[isp_] = 0.0;
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] >= 0) {
	const int isp_ = sp_flag[isp];
	delta_sp_[isp_] = max( delta_sp_[isp_], max( DIST2(xsp[spost[ist][(i+2)%3]],xsp[isp]), DIST2(xsp[isp],xsp[spost[ist][(i+1)%3]]) ) );
      }
    }
  }

  // now take the sqrt...
  for (int isp_ = 0; isp_ < nsp_; ++isp_) {
    assert(delta_sp_[isp_] > 0.0);
    delta_sp_[isp_] = sqrt( delta_sp_[isp_] );
  }

  // -----------------------------------------------------
  // debugging...
  // -----------------------------------------------------

  int isp_debug_ = -1;

  // select isp_debug_ to match a coordinate...
  /*
  //const double xsp_debug[3] = { -0.008333, 1.00833, 0.008333 };
  const double xsp_debug[3] = { 0.0, -3.0, 0.0 };
  double d2_closest;
  for (int isp_ = 0; isp_ < nsp_; ++isp_) {
    const int isp = sposp_[isp_];
    const double this_d2 = DIST2(xsp_debug,xsp[isp]);
    if ((isp_debug_ == -1)||(this_d2 < d2_closest)) {
      isp_debug_ = isp_;
      d2_closest = this_d2;
    }
  }
  cout << "got isp_debug_: " << isp_debug_ << " dist: " << sqrt(d2_closest) << endl;
  */

  // select isp_debug_ to match an isp...
  /*
  for (int isp_ = 0; isp_ < nsp_; ++isp_) {
    const int isp = sposp_[isp_];
    if (isp == 114) {
      isp_debug_ = isp_;
    }
  }
  cout << "got isp_debug_: " << isp_debug_ << endl;
  */

  // or directly set...
  /*
  isp_debug_ = 78;
  cout << "got isp_debug_: " << isp_debug_ << endl;
  */

  // -----------------------------------------------------
  // now find the closest edge to each node...
  // -----------------------------------------------------

  vector<int> candidateVec;
  vector<PossibleMerge> possibleMergeVec;
  set<int> ispMergeSet;

  for (int isp_ = 0; isp_ < nsp_; ++isp_) {
    const int isp = sposp_[isp_];

    const bool debug = isp_ == isp_debug_;

    if (debug) {

      cout << "working on isp_: " << isp_ << ", isp: " << isp << ", sp_nbr_of_sp_: " << sp_nbr_of_sp_[isp_][0] << " " << sp_nbr_of_sp_[isp_][1] << endl;

      FILE * fp = fopen("p.dat","w");
      if (sp_nbr_of_sp_[isp_][0] >= 0)
	fprintf(fp,"%lf %lf %lf\n",xsp[sp_nbr_of_sp_[isp_][0]][0],xsp[sp_nbr_of_sp_[isp_][0]][1],xsp[sp_nbr_of_sp_[isp_][0]][2]);
      fprintf(fp,"%lf %lf %lf\n",xsp[isp][0],xsp[isp][1],xsp[isp][2]);
      if (sp_nbr_of_sp_[isp_][1] >= 0)
	fprintf(fp,"%lf %lf %lf\n",xsp[sp_nbr_of_sp_[isp_][1]][0],xsp[sp_nbr_of_sp_[isp_][1]][1],xsp[sp_nbr_of_sp_[isp_][1]][2]);
      fclose(fp);
    }


    assert(candidateVec.empty());
    oeAdt.buildListForSphere(candidateVec,xsp[isp],delta_sp_[isp_]);

    double d2_min = 1.0E+20;
    assert(possibleMergeVec.empty());

    for (int ii = 0,ii_end=candidateVec.size(); ii < ii_end; ++ii) {
      const int ioe = candidateVec[ii];
      // for this edge to be eligible, its dp must be the opposite of either one of ours...
      const int ist = int( zipTeostVec[ioe] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
      const int i = int( zipTeostVec[ioe]>>60 ); assert((i >= 0)&&(i < 3));
      assert(isEdgeOpen(ist,i));
      const int isp0 = spost[ist][i];
      const int isp1 = spost[ist][(i+1)%3];
      // make sure it is neither of our open edges -- we can't zip ourselves...
      if ( ((sp_nbr_of_sp_[isp_][0] == isp0)&&(isp == isp1)) || ((isp == isp0)&&(sp_nbr_of_sp_[isp_][1] == isp1)) )
	continue;
      // also make sure that atleast ONE of our touching edges dots negatively with this edge...
      const double dx[3] = DIFF(xsp[isp1],xsp[isp0]);
      if ( ( (sp_nbr_of_sp_[isp_][0] >= 0)&&( (xsp[isp][0]-xsp[sp_nbr_of_sp_[isp_][0]][0])*dx[0] +
					      (xsp[isp][1]-xsp[sp_nbr_of_sp_[isp_][0]][1])*dx[1] +
					      (xsp[isp][2]-xsp[sp_nbr_of_sp_[isp_][0]][2])*dx[2] < 0.0 ) ) ||
	   ( (sp_nbr_of_sp_[isp_][1] >= 0)&&( (xsp[sp_nbr_of_sp_[isp_][1]][0]-xsp[isp][0])*dx[0] +
					      (xsp[sp_nbr_of_sp_[isp_][1]][1]-xsp[isp][1])*dx[1] +
					      (xsp[sp_nbr_of_sp_[isp_][1]][2]-xsp[isp][2])*dx[2] < 0.0 ) ) ) {
	// we got one! now figure out how close it is...
	const double dxp[3] = DIFF(xsp[isp],xsp[isp0]);
	const double dp = DOT_PRODUCT(dxp,dx);
	double this_d2,this_t;
	if (dp <= 0.0) {
	  // we are closest to the first point. Note this includes the case when dx is zero, and/or dxp is zero...
	  this_t = 0.0;
	  this_d2 = dxp[0]*dxp[0] + dxp[1]*dxp[1] + dxp[2]*dxp[2];
	}
	else {
	  const double dx2 = DOT_PRODUCT(dx,dx);
	  if (dp >= dx2) {
	    // we are closest to the second point...
	    this_t = 1.0;
	    this_d2 =
	      (xsp[isp][0]-xsp[isp1][0])*(xsp[isp][0]-xsp[isp1][0]) +
	      (xsp[isp][1]-xsp[isp1][1])*(xsp[isp][1]-xsp[isp1][1]) +
	      (xsp[isp][2]-xsp[isp1][2])*(xsp[isp][2]-xsp[isp1][2]);
	  }
	  else {
	    // we are closest to a point along the edge...
	    this_t = dp/dx2;
	    this_d2 = dxp[0]*dxp[0] + dxp[1]*dxp[1] + dxp[2]*dxp[2] - dp*dp/dx2;
	  }
	}

	d2_min = min(d2_min,this_d2);
	possibleMergeVec.push_back(PossibleMerge(this_d2,this_t,ioe,isp0,isp1));

	if (debug) {
	  cout << " got potentially matching edge: " << isp0 << " " << isp1 << " this_t: " << this_t << " this_d: " << sqrt(this_d2) << endl;
	  FILE * fp = fopen("pe.dat","w");
	  fprintf(fp,"%lf %lf %lf\n",xsp[isp0][0],xsp[isp0][1],xsp[isp0][2]);
	  fprintf(fp,"%lf %lf %lf\n",xsp[isp1][0],xsp[isp1][1],xsp[isp1][2]);
	  fclose(fp);
	  fp = fopen("pm.dat","w");
	  fprintf(fp,"%lf %lf %lf\n",
		  (1.0-this_t)*xsp[isp0][0]+this_t*xsp[isp1][0],
		  (1.0-this_t)*xsp[isp0][1]+this_t*xsp[isp1][1],
		  (1.0-this_t)*xsp[isp0][2]+this_t*xsp[isp1][2]);
	  fclose(fp);
	  getchar();
	}

      }
    }
    candidateVec.clear();

    // now consider every merge close to d2_min...

    assert(ispMergeSet.empty());
    for (int ii = 0,ii_end=possibleMergeVec.size(); ii < ii_end; ++ii) {

      // any merge even slightly larger than the minimum should be considered. here we use
      // a factor of 2 on the d2, or 1.414 on dist.

      if (possibleMergeVec[ii].d2 <= 2.0*d2_min) {

	// apply this merge. Here we

	const int isp0 = possibleMergeVec[ii].isp0;
	const int isp1 = possibleMergeVec[ii].isp1;
	const double delta_ed = DIST(xsp[isp0],xsp[isp1]);
	const double d_min = sqrt(d2_min);
	if ((possibleMergeVec[ii].t <= 0.5)&&(possibleMergeVec[ii].t*delta_ed <= 0.5*d_min)) {
	  // push into a local set to avoid duplication in the nodeNodeMergeVec...
	  ispMergeSet.insert(isp0);
	}
	else if ((possibleMergeVec[ii].t > 0.5)&&((1.0-possibleMergeVec[ii].t)*delta_ed <= 0.5*d_min)) {
	  // push into a local set to avoid duplication in the nodeNodeMergeVec...
	  ispMergeSet.insert(isp1);
	}
	else {
	  //push directly into the global edgeNodeMergeVec...
	  edgeNodeMergeVec.push_back(pair<pair<int,double>,int>(pair<int,double>(possibleMergeVec[ii].ioe,possibleMergeVec[ii].t),isp));
	  /*
	  cout << "XXXXXXXXXXXXXXXXXXXXXXXXX ADDING A few HACK edgeNodeMerges" << endl;
	  for (int i = 0; i < 7; ++i) {
	    const double this_t = max(0.05,min(0.95,double(rand())/double(RAND_MAX)));
	    if (fabs(possibleMergeVec[ii].t - this_t) > 0.05)
	      edgeNodeMergeVec.push_back(pair<pair<int,double>,int>(pair<int,double>(possibleMergeVec[ii].ioe,this_t),isp));
	  }
	  */
	}

      }

    }
    possibleMergeVec.clear();

    // now push any nodes in the set into the global nodeNodeMergeVec with isp as the match...
    for (set<int>::iterator iter = ispMergeSet.begin(); iter != ispMergeSet.end(); ++iter) {
      nodeNodeMergeVec.push_back(pair<int,int>(*iter,isp));
      if (debug) cout << " DEBUG adding node match: " << *iter << " " << isp << endl;
    }
    ispMergeSet.clear();

  }

}

void SimpleSurface::splitZipEdgesAndCreateNodeNodeMerges(vector<pair<int,int> >& nodeNodeMergeVec,vector<pair<pair<int,double>,int> >& edgeNodeMergeVec,
							 const vector<uint8>& zipTeostVec) {

  cout << " > splitZipEdgesAndCreateNodeNodeMerges()..." << endl;

  // this routine processes the edgeNodeMergeVec, splits the isdentified edges, and intropduces new nodes, AND creates new
  // node/node pairs in nodeNodeMergeVec...
  if (edgeNodeMergeVec.empty())
    return;

  // we are not going to keep teost constant through this process...
  clearTeost();
  clearDynamicEdgeGroups();
  clearNonManifoldData();

  // we will need to add (at most) one new st and one new sp for each entry in
  // edgeNodeMergeVec. Some of the entries may be so close together that we simply
  // combine them into 1, so the actual size may be less...
  const int nst_new = nst + edgeNodeMergeVec.size();
  const int nsp_new = nsp + edgeNodeMergeVec.size();

  // resize spost and xsp...
  int (*spost0)[3] = spost;
  spost = new int[nst_new][3];
  for (int ist = 0; ist < nst; ++ist) FOR_I3 spost[ist][i] = spost0[ist][i];
  delete[] spost0;

  int *znost0 = znost;
  znost = new int[nst_new];
  for (int ist = 0; ist < nst; ++ist) znost[ist] = znost0[ist];
  delete[] znost0;

  // also the subzone indexing...
  szost.resize(nst_new);

  double (*xsp0)[3] = xsp;
  xsp = new double[nsp_new][3];
  for (int isp = 0; isp < nsp; ++isp) FOR_I3 xsp[isp][i] = xsp0[isp][i];
  delete[] xsp0;

  // sort nodes along each edge. This puts everything in edge-t-order where
  // t is the linear weight from ist0 to ist1 along the edge...

  sort(edgeNodeMergeVec.begin(),edgeNodeMergeVec.end());

  vector<pair<pair<int,double>,int> >::const_iterator iter_end = edgeNodeMergeVec.begin();
  while (iter_end != edgeNodeMergeVec.end()) {
    vector<pair<pair<int,double>,int> >::const_iterator iter_begin = iter_end;
    const int ioe = iter_begin->first.first;
    const int ist = int( zipTeostVec[ioe] & MASK_60BITS ); assert((ist >= 0)&&(ist < nst));
    const int i = int( zipTeostVec[ioe]>>60 ); assert((i >= 0)&&(i < 3));
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    const int isp2 = spost[ist][(i+2)%3];
    ++iter_end;
    while ((iter_end != edgeNodeMergeVec.end())&&(iter_end->first.first == ioe))
      ++iter_end;
    // we now have a sorted range of edge splits associated with this edge...
    for (vector<pair<pair<int,double>,int> >::const_iterator iter = iter_begin; iter != iter_end; ++iter) {
      // add the node...
      const int isp_new = nsp++;
      FOR_J3 xsp[isp_new][j] = (1.0-iter->first.second)*xsp[isp0][j] + iter->first.second*xsp[isp1][j];
      nodeNodeMergeVec.push_back(pair<int,int>(isp_new,iter->second)); // link the new point just created with the node that produced the edge split
      // and the new tri...
      const int ist_new = nst++;
      spost[ist_new][i]       = isp_new;
      spost[ist_new][(i+1)%3] = isp1; // attach terminal end of edge to isp1 for now (gets changed below if not the last tri)...
      spost[ist_new][(i+2)%3] = isp2; // all tris attach to isp2
      znost[ist_new]          = znost[ist];
      szost[ist_new]          = szost[ist];
      // for the tri, it depends on if we are the first or not...
      if (iter == iter_begin) {
	// modify the connections of the first tri...
	spost[ist][(i+1)%3] = isp_new;
      }
      else {
	// the previous trui was a new one, so adjust its isp1...
	spost[ist_new-1][(i+1)%3] = isp_new; // attach to isp1 for now...
      }
    }
  }
  assert(nsp == nsp_new);
  assert(nst == nst_new);

}

void SimpleSurface::compressZipEdgeNodes(const vector<pair<int,int> >& nodeNodeMergeVec) {

  sp_flag.setLength(nsp);
  for (int isp = 0; isp < nsp; ++isp)
    sp_flag[isp] = isp;

  for (int ii = 0, limit = nodeNodeMergeVec.size(); ii < limit; ++ii) {
    int isp0_ = sp_flag[nodeNodeMergeVec[ii].first];
    while (sp_flag[isp0_] != isp0_)
      isp0_ = sp_flag[isp0_];
    int isp1_ = sp_flag[nodeNodeMergeVec[ii].second];
    while (sp_flag[isp1_] != isp1_)
      isp1_ = sp_flag[isp1_];
    sp_flag[isp0_] = sp_flag[isp1_] = min(isp0_,isp1_);
  }

  int nsp_new = 0;
  double * nsp_merged = new double[nsp]; // nsp_new <= nsp
  for (int isp = 0; isp < nsp; ++isp) {
    nsp_merged[isp] = 0; // for checking
    if (sp_flag[isp] == isp) {
      const int isp_new = nsp_new++;
      nsp_merged[isp_new] = 1;
      if (isp_new != isp) {
	assert(isp_new < isp);
	FOR_I3 xsp[isp_new][i] = xsp[isp][i];
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
      FOR_I3 xsp[isp_new][i] += xsp[isp][i];
      ++nsp_merged[isp_new];
    }
  }

  // normalize xsp...
  for (int isp_new = 0; isp_new < nsp_new; ++isp_new) {
    assert(nsp_merged[isp_new] >= 1);
    if (nsp_merged[isp_new] > 1) {
      const double inv_n = 1.0/(double)nsp_merged[isp_new];
      FOR_I3 xsp[isp_new][i] *= inv_n;
    }
  }

  delete[] nsp_merged;

  // the sp_flag contains a -1 indexing of the new sp's...
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
  }

  if ((nsp != nsp_new)||(nst != nst_new)) {
    clearTeost();
    clearDynamicEdgeGroups();
    clearNonManifoldData();
  }

  // resize...
  nsp = nsp_new;
  nst = nst_new;

  // leave szozn_i

}
