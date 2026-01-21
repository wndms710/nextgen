#include "SimpleSurface.hpp"
#include "DoubleVertex.hpp"
#include "Adt.hpp"
#include "Histogram.hpp"
#include "WebUI.hpp"

void SimpleSurface::dumpZones() const {
  COUT1("SimpleSurface::dumpZones(): zoneVec.size: " << zoneVec.size());
  for (int izone = 0,nzn=zoneVec.size(); izone < nzn; ++izone) {
    if (mpi_rank == 0) {
      cout << " > " << izone << " ";
      zoneVec[izone].dump();
    }
  }
}

int SimpleSurface::isManifold() {
  return teost.isManifold();
}

void SimpleSurface::calcRmax() {
  if (!b_centroid) calcCentroid();
  rmax = 0.f;
  for (int isp = 0; isp < nsp; ++isp) {
    double this_r2 = DIST2(centroid,xsp[isp]);
    rmax = max(rmax,this_r2);
  }
  rmax = sqrt(rmax);
  cout << "Surface::calcRmax(): " << rmax << endl;
  b_rmax = true;
}

double SimpleSurface::getRmax() {
  if (!b_rmax) calcRmax();
  assert(b_rmax);
  return rmax;
}

void SimpleSurface::ensureCentroid() {
  if (!b_centroid) {
    calcCentroid();
    assert(b_centroid);
  }
}

void SimpleSurface::calcCentroid() {
  assert(!b_centroid);
  double my_dbuf[4] = { 0.0, 0.0, 0.0, 0.0 };
  FOR_IST {
    const double n[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    const double area = MAG(n);
    FOR_I3 my_dbuf[i] += area*(xsp[spost[ist][0]][i]+xsp[spost[ist][1]][i]+xsp[spost[ist][2]][i]);
    my_dbuf[3] += area;
  }
  assert(my_dbuf[3] > 0.0);
  FOR_I3 centroid[i] = my_dbuf[i]/(3.0*my_dbuf[3]);
  cout << "Surface::calcCentroid(): " << COUT_VEC(centroid) << endl;
  b_centroid = true;
}

void SimpleSurface::calcCenterOfMass(double xcom[3]) {
  double my_buf[4] = { 0.0, 0.0, 0.0, 0.0 };
  FOR_IST {
    const double this_vol = SIGNED_TET_VOLUME_6(xsp[0],xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    FOR_I3 my_buf[i] += this_vol*(xsp[0][i]+xsp[spost[ist][0]][i]+xsp[spost[ist][1]][i]+xsp[spost[ist][2]][i]);
    my_buf[3] += this_vol;
  }
  FOR_I3 xcom[i] = my_buf[i]/(4.0*my_buf[3]);
  cout << "Center of mass [xcom]: " << COUT_VEC(xcom) << endl;
}

void SimpleSurface::calcFlaggedTrisCentroid(double (&_centroid)[3]) {
  double my_dbuf[4] = { 0.0, 0.0, 0.0, 0.0 };
  FOR_IST {
    if (st_flag[ist]) {
      const double n[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
      const double area = MAG(n);
      FOR_I3 my_dbuf[i] += area*(xsp[spost[ist][0]][i]+xsp[spost[ist][1]][i]+xsp[spost[ist][2]][i]);
      my_dbuf[3] += area;
    }
  }
  assert(my_dbuf[3] > 0.0);
  FOR_I3 _centroid[i] = my_dbuf[i]/(3.0*my_dbuf[3]);
  cout << "Surface::calcFlaggedTrisCentroid(): " << COUT_VEC(_centroid) << endl;
}

void SimpleSurface::getCentroid(double x[3]) {
  if (!b_centroid) calcCentroid();
  assert(b_centroid);
  x[0] = centroid[0];
  x[1] = centroid[1];
  x[2] = centroid[2];
}

void SimpleSurface::ensureBoundingBox() {
  if (!b_boundingBox) {
    calcBoundingBox();
    assert(b_boundingBox);
  }
}

void SimpleSurface::calcBoundingBox() {
  assert(!b_boundingBox);
  double buf[6] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };

  FOR_ISP {
    buf[0] = min(buf[0],xsp[isp][0]);
    buf[1] = min(buf[1],-xsp[isp][0]);
    buf[2] = min(buf[2],xsp[isp][1]);
    buf[3] = min(buf[3],-xsp[isp][1]);
    buf[4] = min(buf[4],xsp[isp][2]);
    buf[5] = min(buf[5],-xsp[isp][2]);
  }

  boundingBox[0] =  buf[0];
  boundingBox[1] =  -buf[1];
  boundingBox[2] =  buf[2];
  boundingBox[3] =  -buf[3];
  boundingBox[4] =  buf[4];
  boundingBox[5] =  -buf[5];
  cout << "Surface::calcBoundingBox(): X " << boundingBox[0] << ":" << boundingBox[1] <<
  ", Y " << boundingBox[2] << ":" << boundingBox[3] <<
  ", Z " << boundingBox[4] << ":" << boundingBox[5] << endl;
  b_boundingBox = true;
}

void SimpleSurface::getBoundingBox(double _bBox[6]) {
  if (!b_boundingBox) calcBoundingBox();
  assert(b_boundingBox);
  _bBox[0] = boundingBox[0];
  _bBox[1] = boundingBox[1];
  _bBox[2] = boundingBox[2];
  _bBox[3] = boundingBox[3];
  _bBox[4] = boundingBox[4];
  _bBox[5] = boundingBox[5];
}

double SimpleSurface::getBoundingBoxRmax() {
  if (!b_boundingBox) calcBoundingBox();
  assert(b_boundingBox);

  return 0.5*sqrt((boundingBox[0]-boundingBox[1])*(boundingBox[0]-boundingBox[1])+
  (boundingBox[2]-boundingBox[3])*(boundingBox[2]-boundingBox[3])+
  (boundingBox[4]-boundingBox[5])*(boundingBox[4]-boundingBox[5]));
}

void SimpleSurface::getBoundingBoxCenter(double _bBoxCenter[3]) {
  if (!b_boundingBox) calcBoundingBox();
  assert(b_boundingBox);
  _bBoxCenter[0] = 0.5*(boundingBox[0]+boundingBox[1]);
  _bBoxCenter[1] = 0.5*(boundingBox[2]+boundingBox[3]);
  _bBoxCenter[2] = 0.5*(boundingBox[4]+boundingBox[5]);
}

void SimpleSurface::getAreaAndVolume(double& volume,double& area) {
  volume = 0.0;
  area = 0.0;

  // use first node as reference location
  double x0[3];
  FOR_I3 x0[i] = xsp[spost[0][0]][i];

  FOR_IST {
    int isp0 = spost[ist][0];
    int isp1 = spost[ist][1];
    int isp2 = spost[ist][2];

    const double n2[3] = TRI_NORMAL_2(xsp[isp0],xsp[isp1],xsp[isp2]);
    area += 0.5*MAG(n2);
    const double dx[3] = DIFF(xsp[isp0],x0);
    volume += DOT_PRODUCT(dx,n2);  // division by 6.0 after loop
  }
  volume /= 6.0;
}

void SimpleSurface::flipTrisRandom() {

  cout << "SimpleSurface::flipTrisRandom()" << endl;

  FOR_IST {
    const double rand_number = double(rand())/double(RAND_MAX);
    if (rand_number > 0.5) {
      const int tmp = spost[ist][0];
      spost[ist][0] = spost[ist][1];
      spost[ist][1] = tmp;
    }
  }

  clearTeost();

}

void SimpleSurface::flagNodesOfFlaggedTris() {
  if (st_flag.countPositive() == 0) {
    COUT2(" > no tris currently selected; skipping node flagging");
    return;
  }

  sp_flag.resize(nsp);
  sp_flag.setAll(0);

  FOR_IST {
    if (st_flag[ist]) {
      FOR_I3 sp_flag[spost[ist][i]] |= 1;
    }
  }
}

void SimpleSurface::deleteIsolatedNodes() {

  if (nsp == 0) return;

  sp_flag.setLength(nsp);
  sp_flag.setAll(-1);  // if stays 0 then don't modify

  int nsp_new = 0;
  FOR_IST {
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == -1) sp_flag[isp] = nsp_new++;
    }
  }

  FOR_IST {
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] >= 0) spost[ist][i] = sp_flag[isp];
    }
  }

  double (*xsp0)[3] = xsp;
  xsp = new double[nsp_new][3];
  int _nsp = 0;
  FOR_ISP {
    if (sp_flag[isp] >= 0) {
      FOR_I3 xsp[sp_flag[isp]][i] = xsp0[isp][i];
      ++_nsp;
    }
  }
  assert(_nsp == nsp_new);
  delete[] xsp0;
  nsp = nsp_new;

}

void SimpleSurface::deleteFlaggedTris() {

  // surface will change, recompute geoms
  b_centroid = false;
  b_rmax = false;
  b_boundingBox = false;
  clearTeost();
  clearNonManifoldData();
  clearDynamicEdgeGroups();
  clearOpenEdgeGroups();
  clearSubzoneData();
  clearZoneData();
  b_update_hidden_subzones = true;  // will update hidden sz in the prune method

  // if we had periodic data, we need to reset it here
  if (pbi) {
    clearPeriodicity();
    WUI(WARN,"Bad news, friend. You need to reset periodicity.");
  }

  if (st_flag.count() == 0) {
    COUT2("no surface tris were selected for deletion; skipping");
    return;
  }
  else {
    COUT2(" > removing " << st_flag.count() << " tris");
  }

  const int nst_new = nst - st_flag.count();
  // resize iterator length, but no need to reallocate spost, znost

  // keep track of which nodes will also be removed
  sp_flag.resize(nsp);
  sp_flag.setAll(-1);  // if stays -1 then throw away node, otherwise put new node index here

  int nsp_new = 0;
  int _nst = 0;
  for (int ist=0; ist<nst; ++ist) {
    if (st_flag[ist] == 0) {
      // keep this tri
      FOR_I3 {
        const int isp = spost[ist][i];
        // only update if first time touched
        if (sp_flag[isp] == -1) sp_flag[isp] = nsp_new++;
        spost[_nst][i] = sp_flag[isp];
      }
      znost[_nst] = znost[ist];
      szost[_nst] = szost[ist];
      ++_nst;
    }
  }
  assert(_nst == nst_new);
  nst = _nst;

  const int n_nodes_deleted = sp_flag.countNegative();
  if (n_nodes_deleted) COUT2(" > removing " << n_nodes_deleted << " associated nodes");

  // resort xsp based on updates to spost
  double (*xsp0)[3] = xsp;
  xsp = new double[nsp_new][3];
  int _nsp = 0;
  for (int isp=0; isp < nsp; ++isp) {
    const int isp_new = sp_flag[isp];
    if (isp_new != -1) {
      FOR_I3 xsp[isp_new][i] = xsp0[isp][i];
      ++_nsp;
    }
  }
  assert(_nsp == nsp_new);
  nsp = _nsp;

  DELETE(xsp0);

  // resize flags as necessary
  st_flag.resize(nst);
  szost.resize(nst);
  sp_flag.resize(nsp);
}

int SimpleSurface::countAndIndexDisjointSurfaceGroups(double& volume,double& area, double (&gcl)[3],vector<double>& groupVolVec,const bool detailed) {

  zone_flag.setLength(zoneVec.size());
  st_flag.setLength(nst);
  st_flag.setAll(-1);

  ensureTeost();

  // TODO: replace this with WUI(MESSAGE,... construction
  UIMessage_ msg(MESSAGE);
  int n_neg_vol = 0;
  volume = 0.0;
  area = 0.0;
  FOR_I3 gcl[i] = 0.0;
  msg.addLine("");
  msg.addLine(" ----- disjoint surface group report ----- ");
  int ngroups = 0;
  while (1) {

    // find a seed face that is still -1...
    int ist_seed = -1;
    FOR_IST {
      if (st_flag[ist] == -1) {
        ist_seed = ist;
        break;
      }
    }

    // we are done if there are no more...
    if (ist_seed == -1) break;

    zone_flag.setAll(0);
    int nOpenEdges = 0;
    int nFaces = 0;
    double group_gcl[3] = { 0.0, 0.0, 0.0 };
    double groupArea = 0.0;
    double groupVolume = 0.0;
    double x0[3] = { 0.0, 0.0, 0.0 };
    bool have_x0 = false;

    const int igroup = ngroups++;

    if (detailed) {
      COUT1("Report for surface group[" << igroup << "]:");
      msg.addLine("");
      msg.addLine(stringstream().flush() << "surface group[" << igroup << "]:");
    }

    st_flag[ist_seed] = igroup;
    stack<int> istStack;
    istStack.push(ist_seed);

    while (!istStack.empty()) {

      // pop the next leaf off the stack...
      const int ist = istStack.top(); istStack.pop();
      assert((ist >= 0)&&(ist < nst));
      assert(st_flag[ist] == igroup);

      // process...
      ++nFaces;
      const int izone = znost[ist];
      if (izone >= 0) zone_flag[izone] = 1;

      int isp0 = spost[ist][0];
      int isp1 = spost[ist][1];
      int isp2 = spost[ist][2];


      const double n2[3] = TRI_NORMAL_2(xsp[isp0],xsp[isp1],xsp[isp2]);
      const double mag_n2 = MAG(n2);

      groupArea += 0.5*mag_n2;
      FOR_I3 group_gcl[i] += 0.5*n2[i];
      // for more accurate volume, store the first tri's
      if (!have_x0) {
        have_x0 = true;
        FOR_I3 x0[i] = xsp[isp0][i];
      }
      else {
        const double dx[3] = DIFF(xsp[isp0],x0);
        groupVolume += DOT_PRODUCT(dx,n2);  // needs additional division by 6, performed after loop
      }

      // check our nbrs...
      FOR_I3 {
        int ist_nbr,i_nbr,orient_nbr;
        if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
          if (st_flag[ist_nbr] == -1) {
            st_flag[ist_nbr] = igroup;
            istStack.push(ist_nbr);
          }
          else {
            assert( st_flag[ist_nbr] == igroup );
          }
        }
        else {
          // either open edge, or a multi-edge, but break anyways...
          ++nOpenEdges;
        }
      }

    }
    groupVolume /= 6.0;
    groupVolVec.push_back(groupVolume);

    volume += groupVolume;
    area += groupArea;
    FOR_I3 gcl[i] += group_gcl[i];

    if (detailed) {
      // group report
      COUT1(" > area: " << groupArea );
      if (groupVolume < 0.0) {
        ++n_neg_vol;
        COUT1(" ! volume: " << groupVolume << " (ensure negative is appropriate)");
      }
      else {
        COUT1(" > volume: " << groupVolume);
      }
      COUT1(" > gcl: " << COUT_VEC(group_gcl) );

      // webUI output
      msg.addLine(stringstream().flush() << "    - surface area: " << groupArea);
      if (groupVolume < 0.0) msg.addLine(stringstream().flush() << "    ! volume: " << groupVolume << " (negative)");
      else msg.addLine(stringstream().flush() << "    - volume: " << groupVolume);
      msg.addLine(stringstream().flush() << "    - gcl: " << COUT_VEC(group_gcl));
      msg.addLine(stringstream().flush() << "    - # of tris: " << nFaces);
      if (nOpenEdges) msg.addLine(stringstream().flush() << "    - # of open edges: " << nOpenEdges);

      stringstream zone_list;
      zone_list << " > participating zones: ";
      for (int izn=0,nzn=zone_flag.size(); izn<nzn; ++izn) {
        if (zone_flag[izn]) zone_list << zoneVec[izn].getName() << ",";
      }
      msg.addLine(zone_list);
    }
  }

  COUT1("Surface group summary:");
  COUT1(ngroups << " disjoint surface groups were detected");
  msg.addLine("");
  msg.addLine(stringstream().flush() << " > " << ngroups << " disjoint surface groups were detected");
  if (n_neg_vol) {
    msg.addLine(stringstream().flush() << " > " << n_neg_vol << " of these groups have negative volume");
    COUT1(n_neg_vol << " of these groups have negative volume");
  }
  if (!detailed) msg.addLine(stringstream().flush() << " > rerun using WITH_DETAILS for additional information");

  WebUI::webUIOutput.add_(msg);
  return ngroups;
}

void SimpleSurface::updateSubzoneData(const int isz) {
  // checks member tri list, adds and checks parent zone id,
  // and updates subzone normal and area

  ensureStoszSzosz(); // maybe its better to rebuild SubzoneData instead of calling this fxn
  ensureSubzoneData();

  assert((isz >= 0)&&(isz < nsz));
  subzoneDataVec[isz].zero();

  for (int toz = stosz_i[isz]; toz != stosz_i[isz+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert((ist >= 0)&&(ist < nst));
    assert(szost[ist] == isz);
    const int izone = znost[ist];

    if (subzoneDataVec[isz].zone == -1) subzoneDataVec[isz].zone = izone;
    else assert(subzoneDataVec[isz].zone == izone);

    const double * const x0 = xsp[spost[ist][0]];
    const double * const x1 = xsp[spost[ist][1]];
    const double * const x2 = xsp[spost[ist][2]];
    const double this_normal2[3] = TRI_NORMAL_2(x0,x1,x2);
    const double this_area2 = MAG(this_normal2);
    FOR_I3 subzoneDataVec[isz].xc[i]     += this_area2*(x0[i]+x1[i]+x2[i]);
    FOR_I3 subzoneDataVec[isz].normal[i] += 0.5*this_normal2[i];
    subzoneDataVec[isz].area             += 0.5*this_area2;
  }

  // normalize centroid
  assert(subzoneDataVec[isz].area > 0.0);
  FOR_I3 subzoneDataVec[isz].xc[i] /=  subzoneDataVec[isz].area*6.0;

}

void SimpleSurface::computeSurfaceCurvatureGaussBonnet(const int n_filter) {
  COUT2("SimpleSurface::computeSurfaceCurvatureGuassBonnet()");

  if (sp_buf == NULL) sp_buf = new double[nsp];

  for (int isp = 0; isp < nsp; ++isp) sp_buf[isp] =2.0*M_PI;

  double dx1[3];
  double dx2[3];
  double dx3[3];
  double length[3];
  double theta[3];
  double thetatot;
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      dx1[i]=xsp[spost[ist][1]][i]-xsp[spost[ist][0]][i];
      dx2[i]=xsp[spost[ist][2]][i]-xsp[spost[ist][1]][i];
      dx3[i]=xsp[spost[ist][0]][i]-xsp[spost[ist][2]][i];
      theta[i]=0.0;
      length[i]=0.0;
    }
    // And normalize the vectors
    FOR_I3 {
      length[0] += dx1[i]*dx1[i];
      length[1] += dx2[i]*dx2[i];
      length[2] += dx3[i]*dx3[i];
    }
    FOR_I3 {
      theta[0] -= dx3[i]*dx1[i];
      theta[1] -= dx1[i]*dx2[i];
      theta[2] -= dx2[i]*dx3[i];
      length[i]=sqrt(length[i]);
    }
    thetatot=0.0;
    theta[0]=acos(theta[0]/(length[0]*length[2]));
    theta[1]=acos(theta[1]/(length[0]*length[1]));
    theta[2]=acos(theta[2]/(length[1]*length[2]));
    FOR_I3 thetatot += theta[i];
    FOR_I3 sp_buf[spost[ist][i]] -= theta[i];
  }

  double theta_sum=0.0;
  for (int isp = 0; isp < nsp; ++isp) theta_sum=theta_sum+sp_buf[isp];

  const int EC=nsp-(3*nst/2)+nst;
  const int genus=1-EC/2;
  COUT2(" > Euler characteristic, genus: "<< EC <<", "<< genus);
  COUT2(" > angle sum :" << theta_sum);

  // tri based filter
  if (st_buf == NULL) st_buf = new double[nst];
  sp_flag.setLength(nsp);

  for (int ifl=0; ifl < n_filter; ++ifl) {

    // push nodal data to tris
    for (int ist = 0; ist < nst; ++ist) {
      st_buf[ist] = 0.0;
      FOR_I3 st_buf[ist] += sp_buf[spost[ist][i]];
      st_buf[ist] /= 3.0;
    }

    // clear nodal data
    for (int isp = 0; isp < nsp; ++isp) sp_buf[isp] = 0.0;
    sp_flag.setAll(0);

    // push tri data to nodes
    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
        sp_buf[spost[ist][i]] += st_buf[ist];
        sp_flag[spost[ist][i]]++;
      }
    }

    // average tris at nodes
    for (int isp=0; isp < nsp; ++isp) sp_buf[isp] /= (double)sp_flag[isp];
  }

}

void SimpleSurface::computeSurfaceCurvatureOld() {
  COUT2("SimpleSurface::computeSurfaceCurvature()");

  ensureStoszSzosz();

  double (*sp_normal)[3] = new double[nsp][3];

  for (int isp=0; isp<nsp; ++isp) {
    FOR_I3 sp_normal[isp][i] = 0.0;
  }

  for (int ist=0; ist<nst; ++ist) {
    const double this_n[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    FOR_I3 FOR_J3 sp_normal[spost[ist][i]][j] += this_n[j];
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

  double * curv_st = new double[nst];  // in stosz_v indexing

  for (int isz=0; isz<nsz; ++isz) {
    for (int toz = stosz_i[isz]; toz != stosz_i[isz+1]; ++toz) {
      const int ist = stosz_v[toz];

      curv_st[toz] = 0.0;
      double n2[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
      const double nn2   = MAG(n2);
      if (nn2 > 0.0) {
        FOR_I3 n2[i] /= nn2;

        FOR_I3 { //edges
          const int ii0      = spost[ist][i];
          const int ii1      = spost[ist][(i+1)%3];
          const double ds[3] = DIFF(xsp[ii1],xsp[ii0]);
          const double dx[3] = CROSS_PRODUCT(n2,ds);
          assert( fabs(MAG(ds)-MAG(dx)) < 1.0e-14);
          const double P0    = DOT_PRODUCT(sp_normal[ii0],dx);
          const double P1    = DOT_PRODUCT(sp_normal[ii1],dx);
          curv_st[toz]       += 0.5*(P0+P1);
        }

        curv_st[toz] *= -1.0/nn2;
        curv_st[toz]  = fabs(curv_st[toz]);
        curv_st[toz]  = min(1.0/(curv_st[toz]+1.0e-12),100.0);
      }
      else {
        // invalid tri normal; set to largest curvature
        curv_st[toz] = 100.0;
      }

    }
  }
  delete[] sp_normal;
  COUT2(" > curvature currently available in images via VAR_ON_SURFACE = G_ST");

  if (st_buf == NULL) st_buf = new double[nst];

  for (int isz = 0; isz < nsz; ++isz) {
    const int ist0 = stosz_i[isz];
    const int isz_nst = stosz_i[isz+1] - stosz_i[isz];
    Histogram h(MPI_COMM_NULL);
    h.setLog(true);
    h.add(&curv_st[ist0],isz_nst);
    h.reduce();

    // values are reduced.. we can perform some diagnostics
    const double r_10 = h.binOfCdf(0.10);
    // const double r_90 = h.binOfCdf(0.90);
    for (int toz = stosz_i[isz]; toz != stosz_i[isz+1]; ++toz) {
      const int ist = stosz_v[toz];
      st_buf[ist] = r_10;
    }
  }

  DELETE(curv_st);
}

void SimpleSurface::computeSubzoneCurvature(const double radius_cdf,const bool dump_histogram) {
  if (st_buf == NULL) st_buf = new double[nst];
  FOR_IST st_buf[ist] = double(znost[ist]);
  return;

  ensureSubzoneData();
  ensureStoszSzosz();

  COUT2("SimpleSurface::computeSubzoneCurvature()");
  double (*sp_normal)[3] = new double[nsp][3];
  double * curv_st = new double[nst];

  const double start_time = MPI_Wtime();
  for (int isz=0; isz<nsz; ++isz) {
    computeSubzoneCurvatureForSubzone(isz,sp_normal,curv_st);
  }
  delete[] sp_normal;
  const double time1 = MPI_Wtime()-start_time;
  COUT2(" > computed curvature on tris: " << time1 << "s");

  // histogram and viz output
  char filename[128];
  if (dump_histogram) COUT2(" > writing subzone histograms to files");
  for (int isz = 0; isz < nsz; ++isz) {
    const int ist0 = stosz_i[isz];
    const int isz_nst = stosz_i[isz+1] - stosz_i[isz];
    Histogram h(MPI_COMM_NULL);
    h.setLog(true);
    h.add(&curv_st[ist0],isz_nst);
    h.reduce();
    if (dump_histogram) {
      sprintf(filename,"kappa_isz_%05d.dat",isz);
      h.write(filename);
    }

    // values are reduced.. we can perform some diagnostics
    const double r_threshold = h.binOfCdf(radius_cdf);

    subzoneDataVec[isz].curvature = r_threshold;
    subzoneDataVec[isz].has_curvature = true;

    for (int toz = stosz_i[isz]; toz != stosz_i[isz+1]; ++toz) {
      const int ist = stosz_v[toz];
      st_buf[ist] = r_threshold;
    }
    COUT2(" > curvature per subzone available in images via VAR_ON_SURFACE = G_ST");
  }
  delete[] curv_st;
  const double time2 = MPI_Wtime()-start_time-time1;
  COUT2(" > processed subzone histograms: " << time2 << "s");
}

void SimpleSurface::computeSubzoneCurvatureForSubzone(const int isz,double (*sp_normal)[3],double * curv_st) {

  ensureStoszSzosz();
  assert((isz >= 0)&&(isz < nsz));

  // put subzone-based normals at the participating nodes
  for (int isp=0; isp<nsp; ++isp) {
    FOR_I3 sp_normal[isp][i] = 0.0;
  }

  // loop over all tris in subzone
  for (int toz = stosz_i[isz]; toz != stosz_i[isz+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert((ist >= 0)&&(ist < nst));
    assert(szost[ist] == isz);

    const double this_n[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    FOR_I3 FOR_J3 sp_normal[spost[ist][i]][j] += this_n[j];
  }

  for (int isp=0; isp<nsp; ++isp) {
    const double mag = MAG(sp_normal[isp]);
    if ( mag > 0.0) FOR_I3 sp_normal[isp][i] /= mag;
    else FOR_I3 sp_normal[isp][i] = 0.0;
  }

  // curv_st is indexed to toz for contiguous access later
  for (int toz = stosz_i[isz]; toz != stosz_i[isz+1]; ++toz) {
    const int ist = stosz_v[toz];

    curv_st[toz] = 0.0;
    double n2[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    const double nn2   = MAG(n2);
    if (nn2 > 0.0) {
      FOR_I3 n2[i] /= nn2;

      FOR_I3 { //edges
        const int ii0      = spost[ist][i];
        const int ii1      = spost[ist][(i+1)%3];
        const double ds[3] = DIFF(xsp[ii1],xsp[ii0]);
        const double dx[3] = CROSS_PRODUCT(n2,ds);
        assert( fabs(MAG(ds)-MAG(dx)) < 1.0e-14);
        const double P0    = DOT_PRODUCT(sp_normal[ii0],dx);
        const double P1    = DOT_PRODUCT(sp_normal[ii1],dx);
        curv_st[toz]         += 0.5*(P0+P1);
      }

      curv_st[toz] *= -1.0/nn2;
      curv_st[toz]  = fabs(curv_st[toz]);
      curv_st[toz]  = min(1.0/(curv_st[toz]+1.0e-12),100.0); // go back to length..
    }
    else {
      // invalid tri normal
      curv_st[toz] = 100.0;
    }
  }
}

void SimpleSurface::buildSzoznFromLocalSzost() {
  // builds szozn_i assuming a zone-local (0-indexed within each zone) indexing szost

  // count subzones in each zone
  zone_flag.setLength(zoneVec.size());
  zone_flag.setAll(-1);
  for (int ist=0; ist < nst; ++ist) {
    const int izone = znost[ist];
    if (szost[ist] > zone_flag[izone]) zone_flag[izone] = szost[ist];
  }
  // assert(zone_flag.countNegative() == 0);  // all zones had a subzone  // will be pruned?

  // increment by one so zone_flag carries nsz within each zone
  nsz=0;
  for (int izone=0, nzn=zone_flag.getLength(); izone < nzn; ++izone) {
    ++zone_flag[izone];
    nsz+=zone_flag[izone];
  }

  // update znost zone-local indices to global indices
  szozn_i.resize(zoneVec.size()+1);
  szozn_i[0] = 0;
  for (int izone=0, nzn=zoneVec.size(); izone < nzn; ++izone) {
    szozn_i[izone+1] = szozn_i[izone] + zone_flag[izone];
    // cout << "szozn["<<izone<<"] " << szozn_i[izone+1] << endl;
  }
  // COUT2("buildSzoznFromLocalSzost nsz: " << nsz);

  assert(szozn_i[zoneVec.size()] == nsz);
}

void SimpleSurface::buildSzoznFromGlobalSzost() {
  // builds szozn_i assuming a global indexing of szost

  // find min sz index in each zone (start index)
  zone_flag.setLength(zoneVec.size());
  zone_flag.setAll(-1);
  nsz=0;
  for (int ist=0; ist < nst; ++ist) {
    const int izone = znost[ist];
    if ((zone_flag[izone] == -1) || (szost[ist] < zone_flag[izone])) zone_flag[izone] = szost[ist];
    if (szost[ist] > nsz) nsz = szost[ist];
  }
  assert(zone_flag.countNegative() == 0);  // all zones had a subzone

  ++nsz;  // increment by one to account for number of subzones

  // update znost zone-local indices to global indices
  szozn_i.resize(zoneVec.size()+1);
  for (int izone=0, nzn=zoneVec.size(); izone < nzn; ++izone) {
    szozn_i[izone] = zone_flag[izone];
    // cout << "szozn["<<izone<<"] " << szozn_i[izone] << endl;
  }
  szozn_i[zoneVec.size()] = nsz;
  // cout << "szozn["<<zoneVec.size()<<"] " << szozn_i[zoneVec.size()] << endl;
  // cout << "nsz: " << nsz << endl;

  assert(szozn_i[zoneVec.size()] == nsz);
}

void SimpleSurface::globalSzostToLocal() {
  // go from globally indexed szost to locally indexed (0-indexed within each zone)
  // assumes szozn_i has been built properly
  for (int ist=0; ist<nst; ++ist) {
    szost[ist] -= szozn_i[znost[ist]];
  }
}

void SimpleSurface::localSzostToGlobal() {
  // go from locally indexed (0-indexed within each zone) szost to globally indexed
  // assumes szozn_i has been built properly
  for (int ist=0; ist<nst; ++ist) {
    szost[ist] += szozn_i[znost[ist]];
  }
}

void SimpleSurface::deleteSelectedSubzones(const vector<int>& subzonesVec) {

  sz_flag.setLength(nsz);
  sz_flag.setAll(0);
  for (int i = 0, limit = subzonesVec.size(); i < limit; ++i) {
    const int isz = subzonesVec[i];
    if (isz < 0 || isz >= nsz) {
      CWARN("subzone index out-of-bounds: " << isz);
    }
    else {
      // flag tris in this subzone
      // if (stosz_i == NULL) {
      //   CWARN("stosz_i/v not available; exiting");
      //   return;
      // }
      // for (int toz = stosz_i[isz]; toz != stosz_i[isz+1]; ++toz) {
      //   const int ist = stosz_v[toz];
      //   st_flag[ist] = 1;
      // }
      sz_flag[isz] = 1;
    }
  }

  if (sz_flag.count() == 0) {
    COUT2("no valid subzones were selected for deletion; skipping");
    return;
  }
  // valid subzones were selected
  st_flag.setLength(nst);
  st_flag.setAll(0);
  for (int ist=0; ist<nst; ++ist) {
    if (sz_flag[szost[ist]]) st_flag[ist] = 1;
  }

  selectedSubzoneVec.clear();
  deleteFlaggedTris();  // b_update_hidden_subzones is set within
  pruneEmptyZonesAndSubzones();
  clearDynamicEdgeGroups();
}

// Eikonal solve Fast Marching Method Loop
void SimpleSurface::fmmLoop(const int * stosp_i, const int * stosp_v, MinHeap& trialHeap) {

  // get trial point from close vector...
  pair<double,int> trial = trialHeap.extractMin();
  const int trial_isp = trial.second;
  trialHeap.deleteMin(); // remove from trial
  assert(sp_flag[trial_isp] == 1);
  sp_flag[trial_isp] = 2; // tag as alive
  //cout << "trial " << trial_isp << " g " << sp_buf[trial_isp] << " g again " << trial.first << endl;

  // mark trial_isp neighbors that aren't alive as close...
  for (int tos = stosp_i[trial_isp]; tos != stosp_i[trial_isp+1]; ++tos) {
    const int ist = stosp_v[tos];
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == 0) {
        sp_flag[isp] = 1; // tag as close
        trialHeap.insert(pair<double,int>(HUGE_VAL,isp)); // push onto heap
      }
    }
  }

  /*
  // XXX Can we do this better using up/downHeap???
  for (int index = 0, limit = trialHeap.size(); index < limit; ++index) {
    // marked as close and maintain back-pointers
    sp_flag[trialHeap.heap[index].second] = -index-1; // maintain back pointers;
  }
  */

  // recompute g's for all close trial_isp neighbors...
  for (int tos = stosp_i[trial_isp]; tos != stosp_i[trial_isp+1]; ++tos) {
    const int ist = stosp_v[tos];
    // count number of alive points
    int count = 0;
    FOR_I3 {
      const int isp = spost[ist][i];
      //cout << "isp " << isp << " flag " << sp_flag[isp] << " g " << sp_buf[isp] << endl;
      if (sp_flag[isp] == 2) ++count;
    }
    if (count == 3) {
      continue; // skip tris where all the pts are alive
    }
    else if (count == 2) {
      // find two alive points
      int isp_a, isp_b;
      double g_a = HUGE_VAL, g_b = HUGE_VAL; // a smallest alive value
      FOR_I3 {
        const int isp = spost[ist][i];
        const double this_g = sp_buf[isp];
        if (sp_flag[isp] == 2) {
          if (this_g < g_a) {
            isp_b = isp_a;
            g_b = g_a;
            isp_a = isp;
            g_a = this_g;
          }
          else if (this_g < g_b) {
            isp_b = isp;
            g_b = this_g;
          }
        }
      }
      // update close points g
      FOR_I3 {
        const int isp = spost[ist][i];
        if (sp_flag[isp] == 1) {
          // calc g using  g_b and g_a
          assert(g_a <= g_b); // make sure our assumption is correct
          const double dx_ca[3] = DIFF(xsp[isp_a],xsp[isp]);
          const double dx_cb[3] = DIFF(xsp[isp_b],xsp[isp]);
          const double a2 = DOT_PRODUCT(dx_cb,dx_cb);
          const double b2 = DOT_PRODUCT(dx_ca,dx_ca);
          // const double a = sqrt(a2);
          const double b = sqrt(b2);
          const double a_cos_theta = DOT_PRODUCT(dx_ca,dx_cb)/b;
          const double a2_sin2_theta = a2-a_cos_theta*a_cos_theta;
          const double u = g_b-g_a;
          // solve quadratic
          const double q2 = (a2+b2-2.0*b*a_cos_theta);
          const double q1 = 2.0*b*u*(a_cos_theta-b);
          const double q0 = b2*(u*u-a2_sin2_theta);
          const double disc = q1*q1-2.0*q2*q0;
          double t;
          if (disc >= 0.0) {
            t = max((-q1+sqrt(disc))/(2.0*q2),(-q1-sqrt(disc))/(2.0*q2));
          }
          // update g with best available estimate
          if ( (disc >= 0.0) && (u < t) && ( (b*(t-u)/t) > a_cos_theta) && ( (b*(t-u)/t) < (a2/a_cos_theta) ) ) {
            sp_buf[isp] = min(sp_buf[isp],t+g_a);
          }
          else {
            // calc g using simple dist
            sp_buf[isp] = min(sp_buf[isp], g_a+MAG(dx_ca));
            sp_buf[isp] = min(sp_buf[isp], g_b+MAG(dx_cb));
          }
          // update trialHeap
          if (trialHeap.index(isp) >= 0)
            trialHeap.updateDouble(trialHeap.index(isp),sp_buf[isp]);
        }
      }
    }
    else {
      assert(count == 1); // we need at least one alive point
      // find alive point
      int isp_a;
      double g_a;
      FOR_I3 {
        const int isp = spost[ist][i];
        const double this_g = sp_buf[isp];
        if (sp_flag[isp] == 2) {
          isp_a = isp;
          g_a = this_g;
        }
      }
      // update close points g
      FOR_I3 {
        const int isp = spost[ist][i];
        if (sp_flag[isp] == 1) {
          // calc g using simple dist
          sp_buf[isp] = min(sp_buf[isp], g_a+DIST(xsp[isp],xsp[isp_a]));
          // update trialHeap
          trialHeap.updateDouble(trialHeap.index(isp),sp_buf[isp]);
        }
      }
    }
  }
}

// solve Eikonal equation starting from a provided point
void SimpleSurface::calcGeoDistFromPoint(const double x[3]) {

  // build surface-tri-of-surface-point to allow use to go from isp to ied

  int * stosp_i = new int[nsp+1]; // node,edge-of-surface-point
  for (int isp = 0; isp < nsp; ++isp) stosp_i[isp+1] = 0;

  // count number of tris for each node
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      stosp_i[spost[ist][i]+1]++;
    }
  }

  // scan to get disp from counts...
  stosp_i[0] = 0;
  for (int isp = 0; isp < nsp; ++isp) stosp_i[isp+1] += stosp_i[isp];
  assert(stosp_i[nsp] == nst*3);  // only true if manifold

  // populate stosp_v....
  int *stosp_v = new int[stosp_i[nsp]];
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      stosp_v[stosp_i[spost[ist][i]]++] = ist;
    }
  }

  // rewind...
  for (int isp = nsp-1; isp > 0; --isp) stosp_i[isp] = stosp_i[isp-1];
  stosp_i[0] = 0;

  // find closest point and put its g/isp into CloseVec
  int isp_closest = -1;
  double dist_closest = HUGE_VAL;
  for (int isp = 0; isp < nsp; ++isp) {
    const double this_dist = DIST(xsp[isp],x);
    if (this_dist < dist_closest) {
      isp_closest = isp;
      dist_closest = this_dist;
    }
  }
  assert(isp_closest >= 0);

  // store g values...
  if (sp_buf == NULL) sp_buf = new double[nsp];
  sp_flag.setLength(nsp);
  for (int isp = 0; isp < nsp; ++isp) {
    sp_buf[isp] = HUGE_VAL; // g
    sp_flag[isp] = 0; // 1: alive, 0: far, < 0: close
  }

  // following notation from Sethian Fast Marching Method, 1998...
  MinHeap trialHeap(nsp); // nbrs to alive pts
  sp_buf[isp_closest] = dist_closest;
  sp_flag[isp_closest] = 1; // first heap element (-1 indexing)
  trialHeap.insert(pair<double,int>(dist_closest,isp_closest)); // negated because heap is descending

  // perform fmm loops until finished...
  while (trialHeap.size() > 0) {
    //cout << trialHeap.size() << endl;
    fmmLoop(stosp_i,stosp_v,trialHeap);
  }

  // negate g values (or add invert to color map)?
  //for (int isp = 0; isp < nsp; ++isp) sp_buf[isp] = -sp_buf[isp];

  delete[] stosp_i;
  delete[] stosp_v;

}

// solve Eikonal equation starting from flagged points
void SimpleSurface::calcGeoDistFromFlaggedNodes() {
  // assumes sp_flag was set outside of here

  // build surface-tri-of-surface-point to allow use to go from isp to ied

  int * stosp_i = new int[nsp+1]; // node,edge-of-surface-point
  for (int isp = 0; isp < nsp; ++isp) stosp_i[isp+1] = 0;

  // count number of tris for each node
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      stosp_i[spost[ist][i]+1]++;
    }
  }

  // scan to get disp from counts...
  stosp_i[0] = 0;
  for (int isp = 0; isp < nsp; ++isp) stosp_i[isp+1] += stosp_i[isp];
  assert(stosp_i[nsp] == nst*3);  // only true if manifold

  // populate stosp_v....
  int *stosp_v = new int[stosp_i[nsp]];
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      stosp_v[stosp_i[spost[ist][i]]++] = ist;
    }
  }

  // rewind...
  for (int isp = nsp-1; isp > 0; --isp) stosp_i[isp] = stosp_i[isp-1];
  stosp_i[0] = 0;

  // store g values...
  if (sp_buf == NULL) sp_buf = new double[nsp]; // g

  // following notation from Sethian Fast Marching Method, 1998...
  MinHeap trialHeap(nsp); // nbrs to alive pts
  for (int isp = 0; isp < nsp; ++isp) {
    // sp_flag -> 1: alive, 0: far, < 0: close
    if (sp_flag[isp] == 1) {
      sp_buf[isp] = 0.0;
      trialHeap.insert(pair<double,int>(0.0,isp)); // negated because heap is descending
    }
    else {
      sp_buf[isp] = HUGE_VAL;
    }
  }


  // perform fmm loops until finished...
  while (trialHeap.size() > 0) {
    //cout << trialHeap.size() << endl;
    fmmLoop(stosp_i,stosp_v,trialHeap);
  }

  // negate g values (or add invert to color map)?
  //for (int isp = 0; isp < nsp; ++isp) sp_buf[isp] = -sp_buf[isp];

  delete[] stosp_i;
  delete[] stosp_v;

  //writeSpDataToTecplot("test.dat",sp_buf);

}

void SimpleSurface::getSubzoneBoundingBoxCenterAndDiagonal(double _bBox[6],double _bBoxCenter[3],double &_diag,IntFlag show_sz_flag) {
  double buf[6] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };

  for (int ist = 0; ist < nst; ++ist) {
    const int index = szost[ist];  // subzone index
    const bool bShow = show_sz_flag[index];
    if (index >= 0 && bShow) {
      for (int i = 0; i < 3; ++i) {
        int isp = spost[ist][i];
        buf[0] = min(buf[0],xsp[isp][0]);
        buf[1] = min(buf[1],-xsp[isp][0]);
        buf[2] = min(buf[2],xsp[isp][1]);
        buf[3] = min(buf[3],-xsp[isp][1]);
        buf[4] = min(buf[4],xsp[isp][2]);
        buf[5] = min(buf[5],-xsp[isp][2]);
      }
    }
  }

  _bBox[0] = buf[0];
  _bBox[1] = -buf[1];
  _bBox[2] = buf[2];
  _bBox[3] = -buf[3];
  _bBox[4] = buf[4];
  _bBox[5] = -buf[5];

  _bBoxCenter[0] = 0.5*(buf[0]-buf[1]);
  _bBoxCenter[1] = 0.5*(buf[2]-buf[3]);
  _bBoxCenter[2] = 0.5*(buf[4]-buf[5]);

  _diag = sqrt((buf[0]+buf[1])*(buf[0]+buf[1])+
               (buf[2]+buf[3])*(buf[2]+buf[3])+
               (buf[4]+buf[5])*(buf[4]+buf[5]));

  cout << "SimpleSurface::getSubzoneBoundingBoxCenterAndDiagonal(): X " << _bBoxCenter[0] <<
                                                                   ", Y " << _bBoxCenter[1] <<
                                                                   ", Z " << _bBoxCenter[2] <<
                                                                   ", D " << _diag << endl;

}

void SimpleSurface::reportZoneData(const bool subzone_info) {

  ensureZoneData();  // esures subzone data within

  // TODO: replace this with TODO:
  UIMessage_ msg(MESSAGE);
  msg.addLine("----- Zone Information -----");
  int parent_id = -1;
  int nsz_parent = 0;
  for (int ii = 0, limit = subzoneDataVec.size(); ii < limit; ++ii) {
    cout << "Subzone: " << ii << " parent zone: " << zoneVec[subzoneDataVec[ii].zone].getName() << endl;
    subzoneDataVec[ii].dump();

    if (parent_id != subzoneDataVec[ii].zone) {
      // switching to a new parent zone; provide current zone summary
      ++parent_id;
      nsz_parent = szozn_i[parent_id+1]-szozn_i[parent_id];
      assert(parent_id == subzoneDataVec[ii].zone);  // subzones should be contiguous by parent, so should only be incremented by 1

      msg.addLine("");
      msg.addLine(stringstream().flush() << "zone[" << parent_id << "]: " << zoneVec[parent_id].getName());
      msg.addLine(stringstream().flush() << "    - number of subzones: " << nsz_parent);
      double area;
      if (zoneVec[parent_id].getArea(area)) msg.addLine(stringstream().flush() << "    - surface area: " << area);
      double vec3[3];
      if (zoneVec[parent_id].getNormal(vec3)) msg.addLine(stringstream().flush() << "    - normal vec: " << COUT_VEC(vec3));
      if (zoneVec[parent_id].getCentroid(vec3)) msg.addLine(stringstream().flush() << "    - centroid: " << COUT_VEC(vec3));
    }

    if (subzone_info) {
      if (nsz_parent != 1) {
        msg.addLine(stringstream().flush() << " > subzone[" << ii << "]:");
        msg.addLine(stringstream().flush() << "    - surface area: " << subzoneDataVec[ii].area);
        msg.addLine(stringstream().flush() << "    - normal vec: " << COUT_VEC(subzoneDataVec[ii].normal));
        msg.addLine(stringstream().flush() << "    - mag(normal): " << MAG(subzoneDataVec[ii].normal));
        msg.addLine(stringstream().flush() << "    - x_centroid: " << COUT_VEC(subzoneDataVec[ii].xc));
      }
      else {
        // single subzone in zone
        msg.addLine(stringstream().flush() << " > contains a single subzone[" << ii << "] (identical properties)");
      }
    }
  }

  WebUI::webUIOutput.add_(msg);
}

void SimpleSurface::showSurfaceDiagnostics(const bool detailed,const bool subzone_info) {

  // TODO: replace with WUI(MESSAGE,... construct
  UIMessage_ msg(MESSAGE);

  reportZoneData(subzone_info);
  reportMaterialInfo();

  const int n_open_edges = reportOpenEdgeGroupData(detailed);  // ensureTeost is buried within here already so no need to repeat...

  double vol;
  double area;
  double gcl[3];
  vector<double> groupVolVec;
  const int n_groups = countAndIndexDisjointSurfaceGroups(vol,area,gcl,groupVolVec,detailed);

  msg.addLine("");
  msg.addLine(stringstream().flush() << " ----- entire surface report ----- ");
  msg.addLine("");
  msg.addLine(stringstream().flush() << " > total volume: " << vol);
  msg.addLine(stringstream().flush() << " > total surface area: " << area);
  msg.addLine(stringstream().flush() << " > gcl: " << COUT_VEC(gcl) << " (should be small)");
  if (n_groups > 1) msg.addLine(stringstream().flush() << " > number of surface groups: " << n_groups);

  msg.addLine("");
  const bool b_vol = (vol <= 0.0);
  if (b_vol || n_open_edges || teost.hasMisaligned() || teost.hasMulti()) {
    msg.addLine(stringstream().flush() << "WARNING: The following criteria indicate the surface may require additional repairs:");
    if (b_vol) msg.addLine(stringstream().flush() << " > negative total volume");
    if (n_open_edges) {
      assert(teost.hasOpen());
      msg.addLine(stringstream().flush() << " ! open edges present: " << n_open_edges);
      calcGcl(gcl,true);
      msg.addLine(stringstream().flush() << " ! gcl (including open edge groups): " << COUT_VEC(gcl));
    }
    if (teost.hasMisaligned()) msg.addLine(stringstream().flush() << " ! misaligned neighbors (based on tri normal) present");
    if (teost.hasMulti()) msg.addLine(stringstream().flush() << " ! multiply connected edges (more than 2 adjacent tris) present");
  }
  else {
    msg.addLine(stringstream().flush() << "No discernable problems based on GCL, volume, and simple edge-based diagnostics. Check the non-manifold report for additional potential problems.");
  }
  msg.dumpLines();
  WebUI::webUIOutput.add_(msg);

}

void SimpleSurface::getGapThicknesses(const double delta, const bool gap_viz) {

  ensureTeost();
  ensureSubzoneData();

  const double start_time = MPI_Wtime();
  COUT1("Subzoned:Surface::getGapThicknesses()");

  /*
  IntFlag st_show(nst);
  for (int ist=0; ist < nst; ++ist) {
    if (szost[ist] == 1579 || szost[ist] == 9669) st_show[ist] = 1;
    else st_show[ist] = 0;
  }
  writeSelectedFacesByZoneToTecplot("min_sz_gap.dat",st_show);
  getchar();
  */

  // build bbs...
  double (*bbmin)[3] = new double[nst][3];
  double (*bbmax)[3] = new double[nst][3];
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 bbmin[ist][i] = HUGE_VAL;
    FOR_I3 bbmax[ist][i] = -HUGE_VAL;
    const double n_st[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    const double nmag = MAG(n_st);
    FOR_I3 {
      const int isp = spost[ist][i];
      // orig
      FOR_J3 bbmin[ist][j] = fmin(bbmin[ist][j],xsp[isp][j]);
      FOR_J3 bbmax[ist][j] = fmax(bbmax[ist][j],xsp[isp][j]);
      // trans
      FOR_J3 bbmin[ist][j] = fmin(bbmin[ist][j],xsp[isp][j]-0.501*delta*n_st[j]/nmag); // normals outward
      FOR_J3 bbmax[ist][j] = fmax(bbmax[ist][j],xsp[isp][j]-0.501*delta*n_st[j]/nmag);
    }
  }
  // shrink a bit
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      assert(bbmax[ist][i] >= bbmin[ist][i]);
      const double dx_i = bbmax[ist][i]-bbmin[ist][i];
      bbmin[ist][i] += 0.001*dx_i;
      bbmax[ist][i] -= 0.001*dx_i;
    }
  }

  // build adt...
  Adt<double> * stAdt = new Adt<double>(nst,bbmin,bbmax);
  const double adt_time = MPI_Wtime()-start_time;
  COUT2(" > built surface group Adt: " << adt_time << "s");

  // get gap thicknesses...
  vector<int> candidates;
  vector<int>::iterator c_it;
  double * h_st = new double[nst]; // gap thickness
  int * pair_st = new int[nst]; // gap neighbor
  for (int ist = 0; ist < nst; ++ist) {
    h_st[ist] = delta;
    pair_st[ist] = -1; // no neighbor
  }
  for (int ist = 0; ist < nst; ++ist) {
    double x_st[3]; FOR_I3 x_st[i] = (xsp[spost[ist][0]][i]+xsp[spost[ist][1]][i]+xsp[spost[ist][2]][i])/3.0;
    double n_st[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    const double nmag = MAG(n_st);
    FOR_I3 n_st[i] /= -nmag; // take inward normal
    const double d_st = -DOT_PRODUCT(n_st,x_st);
    int ist_nbr[3] = {-1,-1,-1};
    FOR_I3 {
      int i_nbr,orient_nbr;
      getTriNbrData(ist_nbr[i],i_nbr,orient_nbr,ist,i);
    }

    // find candidates for this tri
    candidates.clear();
    stAdt->buildListForBBox(candidates,bbmin[ist],bbmax[ist]);
    for (c_it=candidates.begin(); c_it!=candidates.end(); ++c_it) {
      const int ist2 = *c_it;
      if ((szost[ist] != szost[ist2])&&(ist2 != ist_nbr[0])&&(ist2 != ist_nbr[1])&&(ist2 != ist_nbr[2])&&(ist != ist2)) {
        double n_st2[3] = TRI_NORMAL_2(xsp[spost[ist2][0]],xsp[spost[ist2][1]],xsp[spost[ist2][2]]);
        const double nmag2 = MAG(n_st2);
        FOR_I3 n_st2[i] /= -nmag2; // take inward normal
        if (DOT_PRODUCT(n_st,n_st2) < -0.1) { // limit angle
          bool in_front = true;
          FOR_I3 {
            const double dist = DOT_PRODUCT(xsp[spost[ist2][i]],n_st)+d_st;
            if (dist < 0) {
              in_front = false;
              break;
            }
          }
          if (in_front) {
            double x_st2[3]; FOR_I3 x_st2[i] = (xsp[spost[ist2][0]][i]+xsp[spost[ist2][1]][i]+xsp[spost[ist2][2]][i])/3.0;
            // check if the proj is in tri...
            bool in_tri = false;
            FOR_I3 {
              const double dist = DOT_PRODUCT(xsp[spost[ist2][i]],n_st)+d_st;
              double x_proj[3];
              FOR_J3 x_proj[j] = xsp[spost[ist2][i]][j]-dist*n_st[j];
              in_tri = pointInTriangle(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],x_proj,-0.001);
              if (in_tri) break;
            }
            if (!in_tri) {
              // check center point for off-chance that the points are right on top of one another...
              const double dist = DOT_PRODUCT(x_st2,n_st)+d_st;
              double x_proj[3];
              FOR_J3 x_proj[j] = x_st2[j]-dist*n_st[j];
              in_tri = pointInTriangle(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]],x_proj,-0.001);
            }
            if (in_tri) {
              const double this_dx[3] = DIFF(x_st2,x_st);
              const double this_dn[3] = DIFF(n_st,n_st2); // they are of opposite sign
              double this_h = 0.5*DOT_PRODUCT(this_dx,this_dn);
              /*
              if (this_h <= 0.0) {
                cout << "problem st pair: st_pair( " << ist << "," << ist2 << "), sz_pair( "
                     << szost[ist] << "," << szost[ist2] << "), gap/delta = " << this_h/delta << endl;
                cout << xsp[spost[ist][0]][0] << " " << xsp[spost[ist][0]][1] << " " << xsp[spost[ist][0]][2] << endl;
                cout << xsp[spost[ist][1]][0] << " " << xsp[spost[ist][1]][1] << " " << xsp[spost[ist][1]][2] << endl;
                cout << xsp[spost[ist][2]][0] << " " << xsp[spost[ist][2]][1] << " " << xsp[spost[ist][2]][2] << endl;
                cout << xsp[spost[ist2][0]][0] << " " << xsp[spost[ist2][0]][1] << " " << xsp[spost[ist2][0]][2] << endl;
                cout << xsp[spost[ist2][1]][0] << " " << xsp[spost[ist2][1]][1] << " " << xsp[spost[ist2][1]][2] << endl;
                cout << xsp[spost[ist2][2]][0] << " " << xsp[spost[ist2][2]][1] << " " << xsp[spost[ist2][2]][2] << endl;
                cout << endl;
                this_h = fabs(this_h);
              }
              */
              if (this_h > 0.0) {
                if (this_h < h_st[ist]) {
                  h_st[ist] = this_h;
                  pair_st[ist] = ist2;
                }
                if (this_h < h_st[ist2]) {
                  h_st[ist2] = this_h;
                  pair_st[ist2] = ist;
                }
              }
            }
          }
        }
      }
    }
  }
  delete stAdt;
  delete[] bbmin;
  delete[] bbmax;
  const double gap_time = MPI_Wtime()-start_time-adt_time;
  cout << " > finished estimating gaps: " << gap_time << "s" << endl;

  // agglomoreate to subzones...
  vector<pair<double,pair<int,int> > > h_sz(nsz); // < <h,<isz,isz2> >
  for (int isz = 0; isz < nsz; ++isz) h_sz[isz] = pair<double,pair<int,int> >(delta,pair<int,int>(isz,-1));
  for (int ist = 0; ist < nst; ++ist) {
    const int isz = szost[ist];
    if (h_st[ist] < h_sz[isz].first) {
      h_sz[isz].first = h_st[ist];
      h_sz[isz].second.first = isz;
      const int ist2 = pair_st[ist];
      assert(ist2 >= 0);
      h_sz[isz].second.second = szost[ist2];
    }
  }
  delete[] pair_st;

  // update subzoneDataVec
  for (int isz = 0; isz < nsz; ++isz) {
    subzoneDataVec[isz].gap_h = h_sz[isz].first;
    subzoneDataVec[isz].has_gap_h = true;
  }
  // vizualize the gap thicknes using seed values...
  if (gap_viz) {
    if (st_buf == NULL) st_buf = new double[nst];
    for (int ist = 0; ist < nst; ++ist) {
      st_buf[ist] = h_st[ist]/delta;
    }
  }
  delete[] h_st;

  // sort for some print out to the user
  sort(h_sz.begin(),h_sz.end());
  if (gap_viz) {
    UIMessage_ msg(MESSAGE); // TODO: replace with WUI(MESSAGE,...
    msg.addLine(stringstream().flush() << " ----- gap thickness report ----- ");
    msg.addLine(stringstream().flush() << " > total number of zones and subzones: " << zoneVec.size() << " and " << nsz);
    for (int isz = 0; isz < min(nsz,50); ++isz) {
      if (h_sz[isz].first < delta)
        msg.addLine(stringstream().flush() << " > (" << h_sz[isz].second.first << "," << h_sz[isz].second.second << "): gap/delta = " << h_sz[isz].first/delta);
    }
    WebUI::webUIOutput.add_(msg);
  }
  for (int isz = 0; isz < nsz; ++isz) {
    if (h_sz[isz].first < delta)
      COUT2(" > (" << h_sz[isz].second.first<< "," << h_sz[isz].second.second << "): gap/delta = " << h_sz[isz].first/delta);
  }
  h_sz.clear();

  //IntFlag st_show(nst);
  //for (int ist=0; ist < nst; ++ist) {
    //if (szost[ist] == min_isz) st_show[ist] = 1;
    //else st_show[ist] = 0;
  //}
  //writeSelectedFacesByZoneToTecplot("min_sz_gap.dat",st_show);

}

void SimpleSurface::ensureSubzoneData() {

  if (!b_subzone_data) {
    buildSubzoneData();
    assert(b_subzone_data);
  }

}

void SimpleSurface::buildSubzoneData() {

  assert(!b_subzone_data);

  cout << "SimpleSurface::buildSubzoneData()" << endl;

  assert(nsz == szozn_i[zoneVec.size()]);
  subzoneDataVec.resize(nsz);
  for (int isz = 0; isz < nsz; ++isz) subzoneDataVec[isz].zero();

  sz_flag.resize(nsz);
  sz_flag.setAll(0);
  FOR_IST {
    const int isz = szost[ist];
    assert((isz >= 0)&&(isz < nsz));
    sz_flag[isz] = 1; // use 1 to indicate there are tris in this subzone
    const int izone = znost[ist];
    if (subzoneDataVec[isz].zone == -1) subzoneDataVec[isz].zone = izone;
    assert(subzoneDataVec[isz].zone == izone);
    const double * const x0 = xsp[spost[ist][0]];
    const double * const x1 = xsp[spost[ist][1]];
    const double * const x2 = xsp[spost[ist][2]];
    const double this_normal2[3] = TRI_NORMAL_2(x0,x1,x2);
    const double this_area2 = MAG(this_normal2);
    FOR_I3 subzoneDataVec[isz].xc[i]     += this_area2*(x0[i]+x1[i]+x2[i]);
    FOR_I3 subzoneDataVec[isz].normal[i] += 0.5*this_normal2[i];
    subzoneDataVec[isz].area             += 0.5*this_area2;
  }

  // normalize centroid...

  bool got_zero = false;
  for (int isz = 0; isz < nsz; ++isz) {
    if (subzoneDataVec[isz].area > 0.0) {
      FOR_I3 subzoneDataVec[isz].xc[i] /=  subzoneDataVec[isz].area*6.0;
      sz_flag[isz] = 0; // handled
    }
    else {
      // we had zero area...
      FOR_I3 subzoneDataVec[isz].xc[i] = 0.0;
      FOR_I3 subzoneDataVec[isz].normal[i] = 0.0;
      subzoneDataVec[isz].area = 0.0;
      if (sz_flag[isz] == 1) {
        got_zero = true;
      }
    }
  }

  // if we got a zero, then there were zero area subzones. Try and build the centroid
  // using edge length...
  if (got_zero) {
    FOR_IST {
      const int isz = szost[ist];
      assert((isz >= 0)&&(isz < nsz));
      if (sz_flag[isz] == 1) {
        FOR_I3 {
          const double * const x0 = xsp[spost[ist][i]];
          const double * const x1 = xsp[spost[ist][(i+1)%3]];
          double length = DIST(x0,x1);
          FOR_J3 subzoneDataVec[isz].xc[j]     += length*(x0[j]+x1[j]);
          subzoneDataVec[isz].area             += length;
        }
      }
    }

    // normalize centroid...
    got_zero = false;
    for (int isz = 0; isz < nsz; ++isz) {
      if (sz_flag[isz] == 1) {
        if (subzoneDataVec[isz].area > 0.0) {
          FOR_I3 subzoneDataVec[isz].xc[i] /= subzoneDataVec[isz].area*2.0;
          subzoneDataVec[isz].area = 0.0;
          sz_flag[isz] = 0;
        }
        else {
          got_zero = true;
          FOR_I3 subzoneDataVec[isz].xc[i] = 0.0;
          FOR_I3 subzoneDataVec[isz].normal[i] = 0.0;
          subzoneDataVec[isz].area = 0.0;
        }
      }
    }

    if (got_zero) {
      // if we got a zero, then there were zero length tris -- must be just a single node
      FOR_IST {
        const int isz = szost[ist];
        assert((isz >= 0)&&(isz < nsz));
        if (sz_flag[isz] == 1) {
          FOR_I3 {
            FOR_J3 subzoneDataVec[isz].xc[j] += xsp[spost[ist][i]][j];
          }
          subzoneDataVec[isz].area += 3.0;
        }
      }

      // normalize centroid...
      got_zero = false;
      for (int isz = 0; isz < nsz; ++isz) {
        if (sz_flag[isz] == 1) {
          assert(subzoneDataVec[isz].area > 0.0);
          FOR_I3 subzoneDataVec[isz].xc[i] /= subzoneDataVec[isz].area;
          subzoneDataVec[isz].area = 0.0;
          sz_flag[isz] = 0;
        }
      }
    }
  }
  b_subzone_data = true;
}

void SimpleSurface::clearSubzoneData() {

  cout << "SimpleSurface::clearSubzoneData()" << endl;

  b_subzone_data = false;
  subzoneDataVec.clear();

}

void SimpleSurface::ensureZoneData() {
  if (!b_zone_data) {
    buildZoneData();
    assert(b_zone_data);
  }
}

void SimpleSurface::buildZoneData() {
  cout << "SimpleSurface::buildZoneData()" << endl;
  assert(!b_zone_data);

  // built off of subzone data, so make sure this is available
  ensureSubzoneData();

  for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
    double _xc[3] = {0.0,0.0,0.0};
    double _normal[3] = {0.0,0.0,0.0};
    double _area = 0.0;

    int count=0;
    for (int isz=szozn_i[izn]; isz<szozn_i[izn+1]; ++isz) {
      if (subzoneDataVec[isz].zone != izn) {
        cout << "Warning: subzoneDataVec[isz].zone != izn: " << subzoneDataVec[isz].zone << " " << izn << endl;
      }
      assert(subzoneDataVec[isz].zone == izn);
      FOR_I3 _xc[i] += subzoneDataVec[isz].xc[i]*subzoneDataVec[isz].area;  // area weighted
      FOR_I3 _normal[i] += subzoneDataVec[isz].normal[i];
      _area += subzoneDataVec[isz].area;
      ++count;
    }
    if (!count) {
      CWARN("empty zone " << zoneVec[izn].getName() << " detected (no subzones); skipping");
      continue;
    }
    // normalize centroid by zone area
    FOR_I3 _xc[i] /= _area;

    zoneVec[izn].setCentroid(_xc);
    zoneVec[izn].setNormal(_normal);
    zoneVec[izn].setArea(_area);
  }

  b_zone_data = true;
}

void SimpleSurface::clearZoneData() {

  cout << "SimpleSurface::clearZoneData()" << endl;

  for (vector<SurfaceZone>::iterator it=zoneVec.begin(); it!=zoneVec.end(); ++it) {
    it->clearMetadata();
  }
  b_zone_data = false;
}

void SimpleSurface::ensureStoszSzosz() {

  if (!b_stosz_szosz) {
    buildStoszSzosz();
    assert(b_stosz_szosz);
  }

}

void SimpleSurface::buildStoszSzosz() {

  assert(!b_stosz_szosz);
  ensureTeost();

  cout << "SimpleSurface::buildStoszSzosz()" << endl;

  set<pair<int,int> > flagPairSet;

  assert(nsz == szozn_i[zoneVec.size()]);

  assert(stosz_i == NULL);
  stosz_i = new int[nsz+1];
  for (int isz = 0; isz < nsz; ++isz) stosz_i[isz+1] = 0;

  FOR_IST {
    const int isz = szost[ist];
    assert((isz >= 0)&&(isz < nsz));
    ++stosz_i[isz+1];

    // build the subzone-to-subzone graph using the flagPairSet. We
    // need only insert one of the graph directions here...
    FOR_I3 {
      int ist_nbr,i_nbr,orient_nbr;
      if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
        // only add when nbr is greater...
        if (szost[ist_nbr] > isz){
          flagPairSet.insert(pair<int,int>(isz,szost[ist_nbr]));
        }
      }
    }
  }

  // turn stosz into CSR...

  stosz_i[0] = 0;
  for (int isz = 0; isz < nsz; ++isz) stosz_i[isz+1] += stosz_i[isz];
  assert(stosz_i[nsz] == nst);

  assert(stosz_v == NULL);
  stosz_v = new int[nst];

  FOR_IST {
    const int isz = szost[ist];
    stosz_v[stosz_i[isz]++] = ist;
  }

  for (int isz = nsz-1; isz > 0; --isz) stosz_i[isz] = stosz_i[isz-1];
  stosz_i[0] = 0;

  // check...

  for (int isz = 0; isz < nsz; ++isz) {
    for (int toz = stosz_i[isz]; toz != stosz_i[isz+1]; ++toz) {
      const int ist = stosz_v[toz];
      assert((ist >= 0)&&(ist < nst));
      assert(szost[ist] == isz);
    }
  }

  // use the set to build the fast nbr graph...

  assert(szosz_i == NULL);
  szosz_i = new int[nsz+1];
  for (int isz = 0; isz < nsz; ++isz) szosz_i[isz+1] = 0;

  for (set<pair<int,int> >::iterator iter = flagPairSet.begin(); iter != flagPairSet.end(); ++iter) {
    ++szosz_i[iter->first+1];
    ++szosz_i[iter->second+1];
  }

  szosz_i[0] = 0;
  for (int isz = 0; isz < nsz; ++isz) szosz_i[isz+1] += szosz_i[isz];
  const int szosz_s = szosz_i[nsz];

  assert(szosz_v == NULL);
  szosz_v = new int[szosz_s];

  for (set<pair<int,int> >::iterator iter = flagPairSet.begin(); iter != flagPairSet.end(); ++iter) {
    szosz_v[szosz_i[iter->first]] = iter->second;
    szosz_v[szosz_i[iter->second]] = iter->first;
    ++szosz_i[iter->first];
    ++szosz_i[iter->second];
  }

  for (int isz = nsz-1; isz > 0; --isz) szosz_i[isz] = szosz_i[isz-1];
  szosz_i[0] = 0;

  b_stosz_szosz = true;
}

void SimpleSurface::clearStoszSzosz() {

  cout << "SimpleSurface::clearStoszSzosz()" << endl;

  b_stosz_szosz = false;
  DELETE(szosz_i);
  DELETE(szosz_v);
  DELETE(stosz_i);
  DELETE(stosz_v);

}

void SimpleSurface::resizeNspData(const int nsp_new,const int nsp_old) {
  growNspData(nsp_new,nsp_old);
}

void SimpleSurface::growNspData(const int nsp_new,const int nsp_old) {
  // allocate new sizes for nsp surface objects
  if (xsp == NULL) {
    assert(nsp_old == 0);
    xsp = new double[nsp_new][3];
  }
  else {
    assert(nsp_old > 0);
    double (*ss_xsp0)[3] = xsp;
    xsp = new double[nsp_new][3];
    for (int isp = 0,lim=min(nsp_old,nsp_new); isp < lim; ++isp) {
      FOR_I3 xsp[isp][i] = ss_xsp0[isp][i];
    }
    DELETE(ss_xsp0);
  }
}

void SimpleSurface::resizeNstData(const int nst_new,const int nst_old) {
  growNstData(nst_new,nst_old);
}

void SimpleSurface::growNstData(const int nst_new,const int nst_old) {
  if (spost == NULL && znost == NULL) {
    assert(nst_old == 0);
    spost = new int[nst_new][3];
    znost = new int[nst_new];
  }
  else {
    assert(nst_old > 0);
    int (*ss_spost0)[3] = spost;
    int *ss_znost0 = znost;
    spost = new int[nst_new][3];
    znost = new int[nst_new];
    for (int ist = 0,lim=min(nst_old,nst_new); ist < lim; ++ist) {
      FOR_I3 spost[ist][i] = ss_spost0[ist][i];
      znost[ist] = ss_znost0[ist];
    }
    DELETE(ss_spost0);
    DELETE(ss_znost0);
  }
  szost.resize(nst_new);  // automatically resizes and copies
}
