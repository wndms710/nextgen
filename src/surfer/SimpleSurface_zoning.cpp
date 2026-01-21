#include "SimpleSurface.hpp"
#include "WebUI.hpp"

void SimpleSurface::migrateSubzonesToZones(const vector<int> zone_indices) {

  cout << "migrateSubzonesToZones()" << endl;

  int nzn0 = zoneVec.size();
  IntFlag zn_flag(nzn0);
  zn_flag.setAll(0);
  if (zone_indices.size() == 0) {
    // default to all
    zn_flag.setAll(1);
  }
  else {
    for (int ii=0, ni=zone_indices.size(); ii < ni; ++ii) {
      const int index = zone_indices[ii];
      if ((index>0) && (index<nzn0)) zn_flag[index] = 1;
    }
  }

  vector<SurfaceZone> zoneVecOld(zoneVec);
  zoneVec.clear();

  // global subzone indices shouldn't change;
  // but the zone they are associated to will

  // store new zone index for this subzone
  sz_flag.setLength(nsz);
  sz_flag.setAll(0);

  int new_isz = 0;
  for (int izn=0, nzn=zoneVecOld.size(); izn<nzn; ++izn) {
    string old_zone = zoneVecOld[izn].getName();

    const int n_subzones = szozn_i[izn+1] - szozn_i[izn];
    if (zn_flag[izn] && (n_subzones > 1)) {
      for (int isz=szozn_i[izn], iend=szozn_i[izn+1]; isz<iend; ++isz) {
        int local_index = isz-szozn_i[izn];
        std::ostringstream stringStream;
        stringStream << old_zone << "_" << local_index;
        std::string copyOfStr = stringStream.str();

        sz_flag[isz] = zoneVec.size();
        zoneVec.push_back(copyOfStr);
      }
    }
    else {
      // keep as is
      for (int isz=szozn_i[izn], iend=szozn_i[izn+1]; isz<iend; ++isz) {
        sz_flag[isz] = zoneVec.size();
      }
      zoneVec.push_back(old_zone);
    }
  }

  for (int ist=0; ist<nst; ++ist) {
    znost[ist] = sz_flag[szost[ist]];
  }
  buildSzoznFromGlobalSzost();

  assert(nsz == szozn_i[zoneVec.size()]);

}

void SimpleSurface::setSubzonesToZones() {

  cout << "setSubzonesToZones()" << endl;

  clearSubzoneData();
  clearZoneData();
  clearStoszSzosz();

  // szost just takes znost...
  szost.setLength(nst);
  const int nzn = zoneVec.size();
  for (int ist = 0; ist < nst; ++ist) {
    assert((znost[ist] >= 0)&&(znost[ist] < nzn));
    szost[ist] = znost[ist];
  }

  // one subzone per zone...
  nsz = zoneVec.size();
  szozn_i.setLength(zoneVec.size()+1);
  for (int izn = 0; izn <= nzn; ++izn) szozn_i[izn] = izn;

  selectedSubzoneVec.clear();
  hiddenSubzoneVec.clear();
  clearDynamicEdgeGroups();
}

void SimpleSurface::combineZonesWithSameName() {

  cout << "combineZonesWithSameName()" << endl;

  map<const string,int> zoneMap;
  const int nzn0 = zoneVec.size();
  int * znozn0 = new int[nzn0];
  int * zn0ozn = new int[nzn0]; // larger
  int nzn = 0;
  for (int izn0 = 0; izn0 < nzn0; ++izn0) {
    const string zone_name = zoneVec[izn0].getName();
    map<const string,int>::iterator it = zoneMap.find(zone_name);
    if (it == zoneMap.end()) {
      zn0ozn[nzn] = izn0;
      zoneMap[zone_name] = znozn0[izn0] = nzn++;
    }
    else {
      znozn0[izn0] = it->second;
    }
  }
  zoneMap.clear();

  FOR_IST {
    const int izn0 = znost[ist]; assert((izn0 >= 0)&&(izn0 < nzn0));
    const int izn = znozn0[izn0]; assert((izn >= 0)&&(izn < nzn));
    znost[ist] = izn;
  }
  delete[] znozn0;
  vector<SurfaceZone> zoneVec0 = zoneVec;
  zoneVec.resize(nzn);
  for (int izn = 0; izn < nzn; ++izn) {
    const int izn0 = zn0ozn[izn]; assert((izn0 >= 0)&&(izn0 < nzn0));
    zoneVec[izn] = zoneVec0[izn0];
  }
  delete[] zn0ozn;
  setSubzonesToZones();
}

int SimpleSurface::addNewZone(const string& zonename) {

  // add a new zone to the end of the zoneVec, and set szozn_i as empty...

  assert(nsz == szozn_i[zoneVec.size()]);
  const int izone = zoneVec.size(); // put the new zone at the end
  const string zonename_ = replaceSpacesWithUnderscores(zonename);
  zoneVec.push_back(SurfaceZone(zonename_));
  COUT2(" > creating new zone \"" << zonename_ << "\"");
  szozn_i.resize(izone+2);
  szozn_i[izone+1] = szozn_i[izone];

  // addition of new zones requires rebuilding zone/subzone info
  clearSubzoneData();
  clearZoneData();

  return izone;
}

void SimpleSurface::splitFlaggedSubzonesIntoTris() {

  cout << "splitFlaggedSubzonesIntoTris()" << endl;

  // indicate we need to update hidden subzone indices
  b_update_hidden_subzones = true;

  const int nzn = zoneVec.size();
  map<const int,int>* szoznMap = new map<const int,int>[nzn];
  int * nsz_zn = new int[nzn];
  for (int izn = 0; izn < nzn; ++izn) {
    nsz_zn[izn] = 0;
    for (int isz = szozn_i[izn]; isz != szozn_i[izn+1]; ++isz) {
      if (sz_flag[isz] == 0)
        szoznMap[izn][isz] = nsz_zn[izn]++;
    }
  }
  FOR_IST {
    const int isz0 = szost[ist];
    const int izn = znost[ist];
    if (sz_flag[isz0] == 0) {
      map<const int,int>::iterator it = szoznMap[izn].find(isz0);
      assert(it != szoznMap[izn].end());
      szost[ist] = it->second; // local
    }
    else {
      szost[ist] = nsz_zn[izn]++;
    }
  }
  delete[] nsz_zn;

  int* szozn_i0 = new int[nzn+1];
  for (int izn = 0; izn < nzn+1; ++izn)
    szozn_i0[izn] = szozn_i[izn];
  buildSzoznFromLocalSzost();
  FOR_IST szost[ist] += szozn_i[znost[ist]];

  // because there are SO many tri's I am only going to propogate the hidden subzones that are not flagged
  {
    vector<int> hiddenSubzoneVecNew;
    for (int ii = 0, lim = hiddenSubzoneVec.size(); ii < lim; ++ii) {
      const int isz0 = hiddenSubzoneVec[ii];
      if ((isz0 >= 0)&&(isz0 < szozn_i0[nzn])&&(sz_flag[isz0] == 0)) {
        int izn = 0;
        while (szozn_i0[izn+1] <= isz0)
          ++izn;
        map<const int,int>::iterator it = szoznMap[izn].find(isz0);
        assert(it != szoznMap[izn].end()); // do I need to make this an if?
        hiddenSubzoneVecNew.push_back(it->second+szozn_i[izn]);
      }
    }
    hiddenSubzoneVec.swap(hiddenSubzoneVecNew);
  }
  selectedSubzoneVec.clear();

  delete[] szozn_i0;
  delete[] szoznMap;

  clearSubzoneData();
  clearStoszSzosz();

}

void SimpleSurface::splitFlaggedSubzonesInSphere(const double x[3],const double r) {

  cout << "splitFlaggedSubzonesInSphere()" << endl;

  // indicate we need to update hidden subzone indices
  b_update_hidden_subzones = true;

  const int nzn = zoneVec.size();
  int * nsz_zn = new int[nzn];
  for (int izn = 0; izn < nzn; ++izn) {
    nsz_zn[izn] = 0;
    for (int isz = szozn_i[izn]; isz != szozn_i[izn+1]; ++isz)
      ++nsz_zn[izn];
  }
  for (int ii = 0, lim = hiddenSubzoneVec.size(); ii < lim; ++ii)
    sz_flag[hiddenSubzoneVec[ii]] |= (1<<1);
  map<const int,int>* szoznMap = new map<const int,int>[nzn];
  int* szosz0 = new int[szozn_i[nzn]];
  for (int isz = 0; isz < szozn_i[nzn]; ++isz) szosz0[isz] = -1;
  FOR_IST {
    const int isz = szost[ist];
    const int izn = znost[ist];
    if ((sz_flag[isz] & (1<<0)) == 0) {
      szost[ist] -= szozn_i[izn]; // make local
      szosz0[isz] = szost[ist];
    }
    else {
      bool b_in = false;
      FOR_I3 {
        if (DIST2(xsp[spost[ist][i]],x) < r*r) {
          map<const int,int>::iterator it = szoznMap[izn].find(isz);
          if (it == szoznMap[izn].end()) {
            szost[ist] = nsz_zn[izn]++; // move to the new local subzone
            szoznMap[izn][isz] = szost[ist];
          }
          else {
            szost[ist] = it->second;
          }
          b_in = true;
          break;
        }
      }
      if (!b_in) {
        szost[ist] -= szozn_i[izn]; // make local
        szosz0[isz] = szost[ist];
      }
    }
  }
  delete[] nsz_zn;
  int* szozn_i0 = new int[nzn+1];
  for (int izn = 0; izn < nzn+1; ++izn)
    szozn_i0[izn] = szozn_i[izn];
  buildSzoznFromLocalSzost();
  FOR_IST szost[ist] += szozn_i[znost[ist]];

  hiddenSubzoneVec.clear();
  selectedSubzoneVec.clear();
  // add the updated unsplit subzones...
  for (int isz0 = 0; isz0 < szozn_i0[nzn]; ++isz0) {
    if ((szosz0[isz0] >= 0)&&(sz_flag[isz0] >= 1)) {
      int izn = 0;
      while (szozn_i0[izn+1] <= isz0)
        ++izn;
      const int isz = szosz0[isz0]+szozn_i[izn];
      if (sz_flag[isz0] & (1<<0))
        selectedSubzoneVec.push_back(isz);
      if (sz_flag[isz0] & (1<<1))
        hiddenSubzoneVec.push_back(isz);
    }
  }
  delete[] szosz0;
  delete[] szozn_i0;
  // add the split subzones...
  for (int izn = 0; izn < nzn; ++izn) {
    for (map<const int,int>::iterator it = szoznMap[izn].begin(); it != szoznMap[izn].end(); ++it) {
      const int isz0 = it->first;
      const int isz = it->second + szozn_i[izn];
      assert(sz_flag[isz0] >= 1);
      selectedSubzoneVec.push_back(isz);
      if (sz_flag[isz0] & (1<<1))
        hiddenSubzoneVec.push_back(isz);
    }
  }
  delete[] szoznMap;

  clearSubzoneData();
  clearStoszSzosz();
  pruneEmptyZonesAndSubzones(); // some could be empty

}

void SimpleSurface::splitFlaggedSubzonesIntoDisjointSubzones() {

  cout << "splitFlaggedSubzonesIntoDisjointSubzones()" << endl;

  // indicate we need to update hidden subzone indices
  b_update_hidden_subzones = true;


  if (!gotTeost()) {

    cout << " > trying noTeost version..." << endl;

    // if we don't have teost, use a different approach to compute
    // disjoint subzones that does not involve computing teost...

    vector<pair<int,int> > flagVec;
    {
      map<const pair<int,int>,int> edgeMap;
      for (int ist = 0; ist < nst; ++ist) {
	const int isz = szost[ist];
	if (sz_flag[isz]) {
	  int count = flagVec.size();
	  flagVec.push_back(pair<int,int>(ist,count));
	  FOR_I3 {
	    const int isp0 = spost[ist][i];
	    const int isp1 = spost[ist][(i+1)%3];
	    map<const pair<int,int>,int>::iterator iter = edgeMap.find(pair<int,int>(isp1,isp0));
	    if (iter != edgeMap.end()) {
	      int nbr_count = iter->second;
	      edgeMap.erase(iter);
	      while (nbr_count != flagVec[nbr_count].second)
		nbr_count = flagVec[nbr_count].second;
	      if (nbr_count < count) {
		assert(flagVec[count].second == count);
		flagVec[count].second = nbr_count;
		count = nbr_count;
	      }
	      else if (count < nbr_count) {
		assert(flagVec[nbr_count].second == nbr_count);
		flagVec[nbr_count].second = count;
		nbr_count = count;
	      }
	    }
	    else {
	      edgeMap[pair<int,int>(isp0,isp1)] = count;
	    }
	  }
	}
      }
    }

    cout << " > nst flagged: " << flagVec.size() << endl;

    int ngr = 0;
    for (int ii = 0, limit = flagVec.size(); ii < limit; ++ii) {
      int count = flagVec[ii].second;
      while ((count >= 0)&&(count != flagVec[count].second))
	count = flagVec[count].second;
      if (count >= 0) {
	++ngr;
	flagVec[ii].second = -ngr; // -1,-2,...
      }
      else {
	flagVec[ii].second = count;
      }
    }

    cout << " > ngr: " << ngr << endl;

    double (*n_sz)[3] = new double[ngr][3];
    double (*x_sz)[3] = new double[ngr][3];
    double *area_sz   = new double[ngr];

    for (int igr = 0; igr < ngr; ++igr) {
      FOR_I3 n_sz[igr][i] = 0.0;
      FOR_I3 x_sz[igr][i] = 0.0;
      area_sz[igr] = 0.0;
    }

    for (int ii = 0, limit = flagVec.size(); ii < limit; ++ii) {
      const int ist = flagVec[ii].first; assert((ist >= 0)&&(ist < nst));
      const int igr = -flagVec[ii].second-1; assert((igr >= 0)&&(igr < ngr));
      const double n2[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
      FOR_I3 n_sz[igr][i] += n2[i];
      const double area = MAG(n2);
      area_sz[igr] += area;
      FOR_I3 x_sz[igr][i] += area*(xsp[spost[ist][0]][i] + xsp[spost[ist][1]][i] + xsp[spost[ist][2]][i]); // don't forget about the factor of 3
    }

    for (int igr = 0; igr < ngr; ++igr) {
      FOR_I3 n_sz[igr][i] *= 0.5;
      FOR_I3 x_sz[igr][i] /= area_sz[igr]*3.0;
      area_sz[igr] *= 0.5;
      // report...
      cout << " > group: " << igr << " normal: " << COUT_VEC(n_sz[igr]) << " centroid: " << COUT_VEC(x_sz[igr]) << " area: " << area_sz[igr] << endl;
    }

    delete[] n_sz;
    delete[] x_sz;
    delete[] area_sz;

    /*
      cout << "writing test.dat..." << endl;
      st_flag.resize(nst);
      st_flag.setAll(0);
      for (int ii = 0, limit = flagVec.size(); ii < limit; ++ii) {
      const int ist = flagVec[ii].first;
      st_flag[ist] = 1;
      }
      writeSelectedFacesByZoneToTecplot("test.dat",st_flag);
      getchar();
    */

    cout << " > done noTeost version." << endl;

  }


  // walking the surface will require teost...
  ensureTeost();
  clearSubzoneData();
  clearZoneData();
  clearStoszSzosz();

  // the sz_flag contains the zones to be split. We also want to preserve any
  // existing subzoning...

  const int nzn = zoneVec.size();
  assert(nsz == szozn_i[nzn]);
  for (int ist = 0; ist < nst; ++ist) {
    // confirm ranges...
    const int izn = znost[ist];
    assert((izn >= 0)&&(izn < nzn));
    const int isz = szost[ist];
    assert((isz >= szozn_i[izn])&&(isz < szozn_i[izn+1]));
    // change ALL szost to a zero-indexed range in the zone...
    szost[ist] -= szozn_i[izn];
    if (sz_flag[isz]) {
      szost[ist] = -szost[ist]-1; // flip to a -1 indexing, preserving unique subzone index
    }
  }

  stack<int> stStack;
  zone_flag.resize(nzn);
  zone_flag.setAll(0);
  multimap<int,int> hidden_zn_local;
  if (!hiddenSubzoneVec.empty()) {
    for (int izn=0; izn<nzn; ++izn) {
      for (int isz=szozn_i[izn]; isz<szozn_i[izn+1]; ++isz) sz_flag[isz] = izn;  // sz_flag holds zone now
    }

    for (int i=0,end=hiddenSubzoneVec.size(); i<end; ++i) {
      const int isz = hiddenSubzoneVec[i];
      const int local_isz = isz-szozn_i[sz_flag[isz]];
      hidden_zn_local.insert(pair<int,int> (sz_flag[isz],local_isz));
    }
  }


  // switch szozn_i to a count...
  for (int izn = nzn; izn > 0; --izn) szozn_i[izn] -= szozn_i[izn-1];

  multimap<int,int> selected_zn_local;
  for (int ist_seed = 0; ist_seed < nst; ++ist_seed) {
    if (szost[ist_seed] <= -1) {
      // we got a seed that has not been marched yet.
      const int this_izn = znost[ist_seed];
      const int this_isz = szost[ist_seed];
      // const int this_isz_new = zone_flag[this_izn]++; // start new indexing for additionally split zones
      const int this_isz_new = szozn_i[this_izn+1]++; // start new indexing for additionally split zones
      selected_zn_local.insert(pair<int,int> (this_izn,this_isz_new));
      szost[ist_seed] = this_isz_new;
      stStack.push(ist_seed);
      while (!stStack.empty()) {
        const int ist = stStack.top(); stStack.pop();
        // loop through the 3 tri-edge nbrs...
        FOR_I3 {
          int ist_nbr;
          if (getAlignedTriNbr(ist_nbr,ist,i)) {
            // do not cross zone boundaries, or previous subzone boundaries...
            if ( (znost[ist_nbr] == this_izn) && (szost[ist_nbr] == this_isz) ) {
              szost[ist_nbr] = this_isz_new;
              stStack.push(ist_nbr);
            }
          }
        }
      }
    }
  }

  // turn szozn_i back to CSR...
  assert(szozn_i[0] == 0);
  for (int izn = 0; izn < nzn; ++izn) {
    // if (zone_flag[izn]) szozn_i[izn+1] = zone_flag[izn];
    szozn_i[izn+1] += szozn_i[izn];
  }
  nsz = szozn_i[nzn];

  // set which subzones should be highlighted
  selectedSubzoneVec.clear();
  for (multimap<int,int>::iterator it=selected_zn_local.begin(); it!=selected_zn_local.end();) {
    const int izn = it->first;
    const int offset = szozn_i[izn];

    do {
      selectedSubzoneVec.push_back(it->second + offset);
      ++it;
    } while (it->first == izn && it!=selected_zn_local.end());
  }

  // set hidden subzones
  if (hiddenSubzoneVec.size()) {
    hiddenSubzoneVec.clear();
    for (multimap<int,int>::iterator it=hidden_zn_local.begin(); it!=hidden_zn_local.end();) {
      const int izn = it->first;
      const int offset = szozn_i[izn];

      do {
        hiddenSubzoneVec.push_back(it->second + offset);
        ++it;
      } while (it->first == izn && it!=hidden_zn_local.end());
    }
  }

  localSzostToGlobal();

  pruneEmptyZonesAndSubzones();

  clearDynamicEdgeGroups();
}

void SimpleSurface::splitFlaggedSubzonesAtCreaseAngle(const double crease_angle) {

  // NOTE: very similar to above routine. If you make changes here, consider
  // making them above

  cout << "splitFlaggedSubzonesAtCreaseAngle(" << crease_angle << ")" << endl;

  b_update_hidden_subzones = true;

  // walking the surface will require teost...
  ensureTeost();
  clearSubzoneData();
  clearZoneData();
  clearStoszSzosz();

  // the sz_flag contains the zones to be split. We also want to preserve any
  // existing subzoning...

  // unit normals to get crease angle
  double (*unit_n_st)[3] = new double[nst][3];

  const int nzn = zoneVec.size();
  assert(nsz == szozn_i[nzn]);
  for (int ist = 0; ist < nst; ++ist) {
    // confirm ranges...
    const int izn = znost[ist];
    assert((izn >= 0)&&(izn < nzn));
    const int isz = szost[ist];
    assert((isz >= szozn_i[izn])&&(isz < szozn_i[izn+1]));
    // change ALL szost to a zero-indexed range in the zone...
    szost[ist] -= szozn_i[izn];
    if (sz_flag[isz]) {
      szost[ist] = -szost[ist]-1; // flip to a -1 indexing, preserving unique subzone index
      // only need the normal in these ist's...
      const double * const x0 = xsp[spost[ist][0]];
      const double * const x1 = xsp[spost[ist][1]];
      const double * const x2 = xsp[spost[ist][2]];
      const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
      const double mag = MAG(this_n);
      if (mag > 0.0) {
        FOR_I3 unit_n_st[ist][i] = this_n[i]/mag;
      }
      else {
        FOR_I3 unit_n_st[ist][i] = 0.0;
      }
    }
  }

  multimap<int,int> hidden_zn_local;
  if (!hiddenSubzoneVec.empty()) {
    for (int izn=0; izn<nzn; ++izn) {
      for (int isz=szozn_i[izn]; isz<szozn_i[izn+1]; ++isz) sz_flag[isz] = izn;  // sz_flag holds zone now
    }

    for (int i=0,end=hiddenSubzoneVec.size(); i<end; ++i) {
      const int isz = hiddenSubzoneVec[i];
      const int local_isz = isz-szozn_i[sz_flag[isz]];
      hidden_zn_local.insert(pair<int,int> (sz_flag[isz],local_isz));
    }
  }

  // switch szozn_i to a count...
  for (int izn = nzn; izn > 0; --izn) szozn_i[izn] -= szozn_i[izn-1];

  const double dp_tol = cos((180.0-crease_angle)*M_PI/180.0);

  multimap<int,int> selected_zn_local;
  stack<int> stStack;
  for (int ist_seed = 0; ist_seed < nst; ++ist_seed) {
    if (szost[ist_seed] <= -1) {
      // we got a seed that has not been marched yet.
      const int this_izn = znost[ist_seed];
      const int this_isz = szost[ist_seed];
      const int this_isz_new = szozn_i[this_izn+1]++; // note that count is in izn+1
      selected_zn_local.insert(pair<int,int> (this_izn,this_isz_new));
      szost[ist_seed] = this_isz_new;
      stStack.push(ist_seed);
      while (!stStack.empty()) {
        const int ist = stStack.top(); stStack.pop();
        // loop through the 3 tri-edge nbrs...
        FOR_I3 {
          int ist_nbr;
          if (getAlignedTriNbr(ist_nbr,ist,i)) {
            // do not cross zone boundaries, or previous subzone boundaries...
            if ((znost[ist_nbr] == this_izn)&&(szost[ist_nbr] == this_isz)) {
              // check the angle...
              if (DOT_PRODUCT(unit_n_st[ist],unit_n_st[ist_nbr]) >= dp_tol) {
                szost[ist_nbr] = this_isz_new;
                stStack.push(ist_nbr);
              }
            }
          }
        }
      }
    }
  }
  DELETE(unit_n_st);

  // turn szozn_i back to CSR...
  assert(szozn_i[0] == 0);
  for (int izn = 0; izn < nzn; ++izn) szozn_i[izn+1] += szozn_i[izn];
  nsz = szozn_i[nzn];

  // set which subzones should be highlighted
  selectedSubzoneVec.clear();
  for (multimap<int,int>::iterator it=selected_zn_local.begin(); it!=selected_zn_local.end();) {
    const int izn = it->first;
    const int offset = szozn_i[izn];

    do {
      selectedSubzoneVec.push_back(it->second + offset);
      ++it;
    } while (it->first == izn && it!=selected_zn_local.end());
  }

  // set hidden subzones
  if (hiddenSubzoneVec.size()) {
    hiddenSubzoneVec.clear();
    for (multimap<int,int>::iterator it=hidden_zn_local.begin(); it!=hidden_zn_local.end();) {
      const int izn = it->first;
      const int offset = szozn_i[izn];

      do {
        hiddenSubzoneVec.push_back(it->second + offset);

        ++it;
      } while (it->first == izn && it!=hidden_zn_local.end());
    }
  }

  localSzostToGlobal();

  pruneEmptyZonesAndSubzones();

  clearDynamicEdgeGroups();
}

void SimpleSurface::splitFlaggedSubzonesByCoordinate(const int isz_mode,const double* split_data, const string& dest_zone) {
  // split flagged subzones by coordinate value of tri centroid

  clearSubzoneData();
  clearZoneData();
  clearStoszSzosz();
  st_flag.resize(nst);

  vector< pair<int,int> > selected_local;
  selectedSubzoneVec.clear();
  IntFlag znosz(nsz);
  znosz.setAll(-1);
  for (int izone=0, nzn=zoneVec.size(); izone<nzn; ++izone) {
    for (int isz=szozn_i[izone]; isz<szozn_i[izone+1]; ++isz) {
      if (sz_flag[isz]) {
        selected_local.push_back(pair<int,int> (izone,isz-szozn_i[izone]));
      }
      znosz[isz] = izone;  // zone for old globally indexed sz
    }
  }

  b_update_hidden_subzones = true;
  vector< pair<int,int> > hidden_zn_local;
  if (!hiddenSubzoneVec.empty()) {
    // cout << "hidden pre-split: " << endl;
    for (int i=0, end=hiddenSubzoneVec.size(); i<end; ++i) {
      const int global_isz = hiddenSubzoneVec[i];
      if (sz_flag[global_isz]) {
        // need to use global index here for re-indexing, so store with izone = -1 as indication
        hidden_zn_local.push_back(pair<int,int> (-1,global_isz));
      }
      else {
        const int izone = znosz[global_isz];
        const int local_isz = global_isz - szozn_i[izone];
        // cout << "global (local,zone): " << global_isz << " (" << local_isz << "," << izone << ")" << endl;
        hidden_zn_local.push_back(pair<int,int> (izone,local_isz));
      }
    }
    hiddenSubzoneVec.clear(); // will repopulate after isz/zn pair has been updated
  }

  int idir   = -1;
  bool above = false;
  if ( (isz_mode) & (1<<1)) {
    idir = 0;
  } else if ( (isz_mode) & (1<<2)) {
    idir = 1;
  } else {
    assert( (isz_mode) & (1<<3));
    idir = 2;
  }

  if ( (isz_mode) & (1<<4)) {
    above = true;
  } else {
    assert( (isz_mode) & ( 1<<5));
    above = false;
  }

  const double threshold = split_data[0];

  // zone_flag holds original max subzones per zone
  zone_flag.resize(zoneVec.size());
  for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
    zone_flag[izn] = szozn_i[izn+1] - szozn_i[izn];
  }
  assert(zone_flag.countZero() == 0);  // no empty zones
  globalSzostToLocal();
  set<pair<int,int> > toNewZone;  // (zone,local subzone)

  // flag the desired faces in src_zone ..
  const double one_third = 1.0/3.0;
  for (int ist = 0; ist < nst; ++ist) {
    const int global_sz = szozn_i[znost[ist]] + szost[ist];  // local-to-global
    if (sz_flag[global_sz]) {
      double st_centroid = 0.0;
      FOR_I3 st_centroid += xsp[spost[ist][i]][idir];
      st_centroid *= one_third;

      // if in flagged region, then set to new subzone in same zone
      if ( ((st_centroid <= threshold) && !above) || ((st_centroid >= threshold) && above) ) {
        szost[ist] = zone_flag[znost[ist]];  // new largest subzone
        toNewZone.insert(pair<int,int> (znost[ist],szost[ist]));
      }
    }
  }
  buildSzoznFromLocalSzost();
  localSzostToGlobal();

  // at this point szost contains new subzones for tris within criterion
  // toNewZone contains list of (izn,local isz) subzones that should be moved to destination

  // set which subzones should be highlighted, hidden
  vector<pair<int,int> >::iterator it;
  for (it=selected_local.begin(); it!=selected_local.end(); ++it) {
    const int global_isz = szozn_i[it->first] + it->second;
    selectedSubzoneVec.push_back(global_isz);
  }

  if (!hidden_zn_local.empty()) {
    // cout << "post-move hidden: " << endl;
    for (it=hidden_zn_local.begin(); it!=hidden_zn_local.end(); ++it) {
      if (it->first == -1) continue;  // skip hidden zone if part of move - should show it

      const int global_isz = szozn_i[it->first] + it->second;
      // cout << "global (local,zone): " << global_isz << " (" << it->second << "," << it->first << ")" << endl;
      hiddenSubzoneVec.push_back(global_isz);
    }
  }

  // check validity of destination zone
  if (dest_zone == "") {
    // user did not specify; stop at simply subzoning
    COUT2(" > no destination zone specifed; only new subzones were created");
    pruneEmptyZonesAndSubzones();
    clearDynamicEdgeGroups();
    return;
  }

  sz_flag.resize(nsz);
  sz_flag.setAll(0);
  for (set<pair<int,int> >::iterator it=toNewZone.begin(); it!=toNewZone.end(); ++it) {
    const int isz = szozn_i[it->first] + it->second;
    assert(isz < szozn_i[it->first+1]);  // ensure still within same zone
    sz_flag[isz] = 1;
  }

  if (sz_flag.count()) {
    // need to move flagged (new) subzones to destination zone (new or existing)
    int izone = getZoneIndex(dest_zone);
    if (izone == -1) {
      // create a new zone
      izone = addNewZone(dest_zone);
    }
    moveFlaggedSubzonesToZone(izone);
  }
  else {
    COUT2(" > did not find any subzones to move to new zone");
    pruneEmptyZonesAndSubzones();
  }
  clearDynamicEdgeGroups();
}

void SimpleSurface::splitFlaggedSubzonesByMarch(const double* split_data, const string& dest_zone) {
  // march flagged subzones from specified point until crease angle

  const double x_seed[3] = {split_data[0], split_data[1], split_data[2]};
  const double tol       = split_data[3];  // crease criterion
  int itri        = -1;
  double min_dist = 1.0e+15;

  ensureTeost();
  clearSubzoneData();
  clearZoneData();
  clearStoszSzosz();
  st_flag.resize(nst);

  vector< pair<int,int> > selected_local;
  selectedSubzoneVec.clear();
  IntFlag znosz(nsz);
  znosz.setAll(-1);
  for (int izone=0, nzn=zoneVec.size(); izone<nzn; ++izone) {
    for (int isz=szozn_i[izone]; isz<szozn_i[izone+1]; ++isz) {
      if (sz_flag[isz]) {
        selected_local.push_back(pair<int,int> (izone,isz-szozn_i[izone]));
      }
      znosz[isz] = izone;  // zone for old globally indexed sz
    }
  }

  b_update_hidden_subzones = true;
  vector< pair<int,int> > hidden_zn_local;
  if (!hiddenSubzoneVec.empty()) {
    // cout << "hidden pre-split: " << endl;
    for (int i=0, end=hiddenSubzoneVec.size(); i<end; ++i) {
      const int global_isz = hiddenSubzoneVec[i];
      if (sz_flag[global_isz]) {
        // need to use global index here for re-indexing, so store with izone = -1 as indication
        hidden_zn_local.push_back(pair<int,int> (-1,global_isz));
      }
      else {
        const int izone = znosz[global_isz];
        const int local_isz = global_isz - szozn_i[izone];
        // cout << "global (local,zone): " << global_isz << " (" << local_isz << "," << izone << ")" << endl;
        hidden_zn_local.push_back(pair<int,int> (izone,local_isz));
      }
    }
    hiddenSubzoneVec.clear(); // will repopulate after isz/zn pair has been updated
  }

  // find closest tri to seed location
  const double one_third = 1.0/3.0;
  for (int ist = 0; ist < nst; ++ist) {
    st_flag[ist] = 0;

    double x_st[3];
    FOR_I3 x_st[i] = (xsp[spost[ist][0]][i] + xsp[spost[ist][1]][i] + xsp[spost[ist][2]][i])*one_third;
    double this_dist = DIST(x_st,x_seed);
    if ( this_dist < min_dist) {
      min_dist = this_dist;
      itri     = ist;
    }
  }

  assert( (itri >= 0) && (itri < nst));
  // seed tri normal used as the reference for crease criterion
  double tri_normal[3] = TRI_NORMAL_2(xsp[spost[itri][0]],xsp[spost[itri][1]],xsp[spost[itri][2]]);
  double n_mag = MAG(tri_normal);
  if (n_mag > 0.0) FOR_I3 tri_normal[i] /= n_mag;  // otherwise a linear tri

  // zone_flag holds original max subzones per zone
  zone_flag.resize(zoneVec.size());
  for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
    zone_flag[izn] = szozn_i[izn+1] - szozn_i[izn];
  }
  assert(zone_flag.countZero() == 0);  // no empty zones
  globalSzostToLocal();

  // stack-based search until crease angle hit
  stack<int> tris_;
  tris_.push(itri);

  set<pair<int,int> > toNewZone;  // (zone,local subzone)


  while (!tris_.empty()) {

    const int ist = tris_.top();
    tris_.pop();

    if (st_flag[ist] == 0) {
      double n[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
      n_mag = MAG(n);
      if (n_mag > 0.0) FOR_I3 n[i] /= n_mag;

      if ( DOT_PRODUCT(tri_normal,n) > tol) {
        const int global_sz = szozn_i[znost[ist]] + szost[ist];
        if (sz_flag[global_sz]) {
          szost[ist] = zone_flag[znost[ist]];  // push to new subzone
          toNewZone.insert(pair<int,int> (znost[ist],szost[ist]));  // zone, local subzone to move to new zone
        }

        // add unprocessed neighbors to stack
        for (int ii = 0; ii < 3; ++ii) {
          int ist_nbr,i_nbr,orient_nbr;
          if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,ii)) {
            if (st_flag[ist_nbr] == 0) tris_.push(ist_nbr);
          }
        }
      }
      st_flag[ist] = 1;  // flag as processed
    }
  }
  buildSzoznFromLocalSzost();
  localSzostToGlobal();

  // at this point szost contains new subzones for tris within criterion
  // toNewZone contains list of (izn,local isz) subzones that should be moved to destination

  // set which subzones should be highlighted, hidden
  vector<pair<int,int> >::iterator it;
  for (it=selected_local.begin(); it!=selected_local.end(); ++it) {
    const int global_isz = szozn_i[it->first] + it->second;
    selectedSubzoneVec.push_back(global_isz);
  }

  if (!hidden_zn_local.empty()) {
    // cout << "post-move hidden: " << endl;
    for (it=hidden_zn_local.begin(); it!=hidden_zn_local.end(); ++it) {
      if (it->first == -1) continue;  // skip hidden zone if part of move - should show it

      const int global_isz = szozn_i[it->first] + it->second;
      // cout << "global (local,zone): " << global_isz << " (" << it->second << "," << it->first << ")" << endl;
      hiddenSubzoneVec.push_back(global_isz);
    }
  }

  // check validity of destination zone
  if (dest_zone == "") {
    // user did not specify; stop at simply subzoning
    COUT2(" > no destination zone specifed; only new subzones were created");
    pruneEmptyZonesAndSubzones();
    clearDynamicEdgeGroups();
    return;
  }

  sz_flag.resize(nsz);
  sz_flag.setAll(0);
  for (set<pair<int,int> >::iterator it=toNewZone.begin(); it!=toNewZone.end(); ++it) {
    const int isz = szozn_i[it->first] + it->second;
    assert(isz < szozn_i[it->first+1]);  // ensure still within same zone
    sz_flag[isz] = 1;
  }

  if (sz_flag.count()) {
    // need to move flagged (new) subzones to destination zone (new or existing)
    int izone = getZoneIndex(dest_zone);
    if (izone == -1) {
      // create a new zone
      izone = addNewZone(dest_zone);
    }
    moveFlaggedSubzonesToZone(izone);
  }
  else {
    COUT2(" > did not find any subzones to move to new zone");
    pruneEmptyZonesAndSubzones();
  }
  clearDynamicEdgeGroups();
}

void SimpleSurface::splitFlaggedSubzonesByPlane(const double x_plane[3],const double n_plane[3]) {
  // split flagged subzones by coordinate value of tri centroid

  clearSubzoneData();
  clearZoneData();
  clearStoszSzosz();
  st_flag.resize(nst);

  vector< pair<int,int> > selected_local;
  selectedSubzoneVec.clear();
  IntFlag znosz(nsz);
  znosz.setAll(-1);
  for (int izone=0, nzn=zoneVec.size(); izone<nzn; ++izone) {
    for (int isz=szozn_i[izone]; isz<szozn_i[izone+1]; ++isz) {
      if (sz_flag[isz]) {
        selected_local.push_back(pair<int,int> (izone,isz-szozn_i[izone]));
      }
      znosz[isz] = izone;  // zone for old globally indexed sz
    }
  }

  b_update_hidden_subzones = true;
  vector< pair<int,int> > hidden_zn_local;
  if (!hiddenSubzoneVec.empty()) {
    // cout << "hidden pre-split: " << endl;
    for (int i=0, end=hiddenSubzoneVec.size(); i<end; ++i) {
      const int global_isz = hiddenSubzoneVec[i];
      if (sz_flag[global_isz]) {
        // need to use global index here for re-indexing, so store with izone = -1 as indication
        hidden_zn_local.push_back(pair<int,int> (-1,global_isz));
      }
      else {
        const int izone = znosz[global_isz];
        const int local_isz = global_isz - szozn_i[izone];
        // cout << "global (local,zone): " << global_isz << " (" << local_isz << "," << izone << ")" << endl;
        hidden_zn_local.push_back(pair<int,int> (izone,local_isz));
      }
    }
    hiddenSubzoneVec.clear(); // will repopulate after isz/zn pair has been updated
  }

  // zone_flag holds original max subzones per zone
  zone_flag.resize(zoneVec.size());
  for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
    zone_flag[izn] = szozn_i[izn+1] - szozn_i[izn];
  }
  assert(zone_flag.countZero() == 0);  // no empty zones
  globalSzostToLocal();
  set<pair<int,int> > toNewZone;  // (zone,local subzone)

  // flag the desired faces in src_zone ..
  const double one_third = 1.0/3.0;
  for (int ist = 0; ist < nst; ++ist) {
    const int global_sz = szozn_i[znost[ist]] + szost[ist];  // local-to-global
    if (sz_flag[global_sz]) {
      double x_st[3] = {0.0,0.0,0.0};
      FOR_I3 FOR_J3 x_st[j] += xsp[spost[ist][i]][j];
      FOR_I3 x_st[i] *= one_third;
      const double dx[3] = DIFF(x_st,x_plane);

      // if in flagged region, then set to new subzone in same zone
      if (DOT_PRODUCT(n_plane,dx) >= 0.0) {
        szost[ist] = zone_flag[znost[ist]];  // new largest subzone
        toNewZone.insert(pair<int,int> (znost[ist],szost[ist]));
      }
    }
  }
  buildSzoznFromLocalSzost();
  localSzostToGlobal();

  // at this point szost contains new subzones for tris within criterion
  // toNewZone contains list of (izn,local isz) subzones that should be moved to destination

  // set which subzones should be highlighted, hidden
  vector<pair<int,int> >::iterator it;
  for (it=selected_local.begin(); it!=selected_local.end(); ++it) {
    const int global_isz = szozn_i[it->first] + it->second;
    selectedSubzoneVec.push_back(global_isz);
  }

  if (!hidden_zn_local.empty()) {
    // cout << "post-move hidden: " << endl;
    for (it=hidden_zn_local.begin(); it!=hidden_zn_local.end(); ++it) {
      if (it->first == -1) continue;  // skip hidden zone if part of move - should show it

      const int global_isz = szozn_i[it->first] + it->second;
      // cout << "global (local,zone): " << global_isz << " (" << it->second << "," << it->first << ")" << endl;
      hiddenSubzoneVec.push_back(global_isz);
    }
  }

  pruneEmptyZonesAndSubzones();
  clearDynamicEdgeGroups();
}

int SimpleSurface::getZoneIndex(const string& zonename) const {
  int izone, nzones;
  for (izone = 0, nzones=zoneVec.size(); izone < nzones; ++izone) {
    if (zoneVec[izone].getName() == zonename) return izone;
  }
  return -1;  // zone doesn't currently exist
}

void SimpleSurface::renameZone(const string& name, const string& newname) {
  const int izone = getZoneIndex(name);
  if (izone == -1) {
    WUI(WARN,"could not find zone: \\\"" + name + "\\\" to rename; skipping");
  }
  else {
    const string new_name = replaceSpacesWithUnderscores(newname);
    zoneVec[izone].setName(new_name);
  }
}

void SimpleSurface::renameZone(const int index, const string& newname) {
  if (index < 0 || index > int(zoneVec.size())) {
    WUI(WARN,"could not find zone with index: \\\"" << index << "\\\" to rename; skipping");
  }
  else {
    const string new_name = replaceSpacesWithUnderscores(newname);
    zoneVec[index].setName(new_name);
  }
}

void SimpleSurface::moveFlaggedTrisToZone(const int izn_dest) {

  cout << "moveFlaggedTrisToZone()" << endl;

  // for persistent highlighting of selected subzones
  // push subzones to keep highlighted here
  selectedSubzoneVec.clear();

  clearSubzoneData();
  clearZoneData();
  clearStoszSzosz();

  // if we had periodic data, we need to reset it here
  if (pbi) {
    clearPeriodicity();
    WUI(WARN,"Bad news, friend. You need to reset periodicity.");
  }
  // anything else to clear?...

  const int nzn = zoneVec.size();
  assert((izn_dest >= 0)&&(izn_dest < nzn));
  assert(nsz == szozn_i[nzn]);

  // use zone_flag to store the subzone count...
  zone_flag.setLength(nzn);
  for (int izn = 0; izn < nzn; ++izn) {
    zone_flag[izn] = szozn_i[izn+1] - szozn_i[izn];
  }

  sz_flag.setLength(nsz);
  sz_flag.setAll(-1);
  for (int ist = 0; ist < nst; ++ist) {
    // confirm ranges...
    const int izn = znost[ist];
    assert((izn >= 0)&&(izn < nzn));
    const int isz = szost[ist];
    assert((isz >= szozn_i[izn])&&(isz < szozn_i[izn+1]));
    if (st_flag[ist] && (izn != izn_dest)) {
      znost[ist] = izn_dest;
      // check if a flagged member of this subzone has already been indexed in the new zone...
      if (sz_flag[isz] == -1) {
        sz_flag[isz] = zone_flag[izn_dest];
        szost[ist] = zone_flag[izn_dest]++;
      }
      else {
        szost[ist] = sz_flag[isz];
      }
    }
    else {
      szost[ist] -= szozn_i[izn];
    }
  }

  // turn szozn_i back to CSR...
  assert(szozn_i[0] == 0);
  for (int izn = 0; izn < nzn; ++izn) szozn_i[izn+1] = szozn_i[izn] + zone_flag[izn];
  nsz = szozn_i[nzn];

  localSzostToGlobal();

  clearDynamicEdgeGroups();
}

void SimpleSurface::moveFlaggedSubzonesToZone(const int izn_dest) {

  cout << "moveFlaggedSubzonesToZone()" << endl;

  // for persistent highlighting of selected subzones
  // push subzones to keep highlighted here
  selectedSubzoneVec.clear();
  b_update_hidden_subzones = true;

  clearSubzoneData();
  clearZoneData();
  clearStoszSzosz();

  // if we had periodic data, we need to reset it here
  if (pbi) {
    clearPeriodicity();
    WUI(WARN,"Bad news, friend. You need to reset periodicity.");
  }
  // anything else to clear?...

  const int nzn = zoneVec.size();
  assert((izn_dest >= 0)&&(izn_dest < nzn));
  assert(nsz == szozn_i[nzn]);

  vector<int> selected_global;
  IntFlag znosz(nsz);
  znosz.setAll(-1);
  for (int izone=0; izone<nzn; ++izone) {
    for (int isz=szozn_i[izone]; isz<szozn_i[izone+1]; ++isz) {
      if (sz_flag[isz]) selected_global.push_back(isz);
      znosz[isz] = izone;
    }
  }

  vector< pair<int,int> > hidden_zn_local;
  if (!hiddenSubzoneVec.empty()) {
    // cout << "hidden pre-move: " << endl;
    for (int i=0, end=hiddenSubzoneVec.size(); i<end; ++i) {
      const int global_isz = hiddenSubzoneVec[i];
      if (sz_flag[global_isz]) {
        // need to use global index here for re-indexing, so store with izone = -1 as indication
        hidden_zn_local.push_back(pair<int,int> (-1,global_isz));
      }
      else {
        const int izone = znosz[global_isz];
        const int local_isz = global_isz - szozn_i[izone];
        // cout << "global (local,zone): " << global_isz << " (" << local_isz << "," << izone << ")" << endl;
        hidden_zn_local.push_back(pair<int,int> (izone,local_isz));
      }
    }
    hiddenSubzoneVec.clear(); // will repopulate after isz/zn pair has been updated
  }

  // use zone_flag to store the subzone count...
  zone_flag.setLength(nzn);
  for (int izn = 0; izn < nzn; ++izn) {
    zone_flag[izn] = szozn_i[izn+1] - szozn_i[izn];
  }

  // use new_sz as a map from old global isz to new local isz; only for moved subzones
  // if new_sz == -1 then use old local index
  IntFlag new_sz(nsz);
  new_sz.setAll(-1);
  for (int ist = 0; ist < nst; ++ist) {
    // confirm ranges...
    const int izn = znost[ist];
    assert((izn >= 0)&&(izn < nzn));
    const int isz = szost[ist];
    assert((isz >= szozn_i[izn])&&(isz < szozn_i[izn+1]));
    if (sz_flag[isz]) {
      if (izn != izn_dest) {
        znost[ist] = izn_dest;
        // check if a flagged member of this subzone has already been indexed in the new zone...
        if (new_sz[isz] == -1) {
          new_sz[isz] = zone_flag[izn_dest];
          szost[ist] = zone_flag[izn_dest]++;
        }
        else {
          szost[ist] = new_sz[isz];
        }
      }
      else {
        if (new_sz[isz] == -1) new_sz[isz] = isz - szozn_i[izn];  // current local index b/c already in zone
        szost[ist] = new_sz[isz];  // use a local index
      }
    }
    else {
      szost[ist] -= szozn_i[izn];  // old local index
    }
  }

  // turn szozn_i back to CSR...
  assert(szozn_i[0] == 0);
  for (int izn = 0; izn < nzn; ++izn) szozn_i[izn+1] = szozn_i[izn] + zone_flag[izn];
  nsz = szozn_i[nzn];

  // set which subzones should be highlighted
  const int offset = szozn_i[izn_dest];
  for (vector<int>::iterator it=selected_global.begin(); it!=selected_global.end(); ++it) {
    const int isz = new_sz[*it] + offset;
    selectedSubzoneVec.push_back(isz);
  }

  if (!hidden_zn_local.empty()) {
    // cout << "post-move hidden: " << endl;
    for (vector< pair<int,int> >::iterator it=hidden_zn_local.begin(); it!=hidden_zn_local.end(); ++it) {
      if (it->first == -1) continue;  // skip hidden zone if part of move - should show it

      const int global_isz = it->second + szozn_i[it->first];
      // cout << "global (local,zone): " << global_isz << " (" << it->second << "," << it->first << ")" << endl;
      hiddenSubzoneVec.push_back(global_isz);
    }
  }

  localSzostToGlobal();

  pruneEmptyZonesAndSubzones();

  clearDynamicEdgeGroups();
}

bool SimpleSurface::pruneEmptyZonesAndSubzones() {
  COUT2("SimpleSurface::pruneEmptyZonesAndSubzones()");
  // assumes szost is globally indexed
  // and that szozn_i has been built...could call here i suppose...

  // determine which subzones have tris; will by association determine empty zones
  // put tri counts into sz_flag
  sz_flag.resize(nsz);
  sz_flag.setAll(0);
  for (int ist=0; ist < nst; ++ist) {
    const int isz = szost[ist];
    if (sz_flag[isz] == 0) sz_flag[isz] = 1;
  }
  const int n_empty_subzones = sz_flag.countZero();

  // determine if any zones are inappropriately without subzones
  int n_empty_zones = 0;
  for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
    if (szozn_i[izn+1]-szozn_i[izn] == 0) {
      ++n_empty_zones;
    }
  }

  if (n_empty_subzones || n_empty_zones) {
    COUT2(" > found " << n_empty_subzones << " empty subzones; pruning");
    COUT2(" > found " << n_empty_zones << " empty zones in szozn_i; pruning");

    // to keep track of changing indices place selected zones into local-subzone & zone pair; adjust as necessary
    vector< pair<int,int> > tracked_sz_local;
    vector< pair<int,int> >::iterator it;

    // znosz holds zn for each sz
    IntFlag znosz(nsz);
    znosz.setAll(-1);

    for (int izone=0, limit=zoneVec.size(); izone < limit; ++izone) {
      for (int isz=szozn_i[izone]; isz<szozn_i[izone+1]; ++isz) znosz[isz] = izone;
    }

    if (!selectedSubzoneVec.empty()) {
      for (int i=0, end=selectedSubzoneVec.size(); i<end; ++i) {
        const int global_isz = selectedSubzoneVec[i];
        const int izone = znosz[global_isz];
        if ((izone<0) || (izone>int(zoneVec.size()))) continue;  // don't process subzones that seem to be improperly indexed

        const int local_isz = global_isz - szozn_i[izone];
        tracked_sz_local.push_back(pair<int,int> (izone,local_isz+1));  // selected are 1-indexed
      }

      selectedSubzoneVec.clear(); // will repopulate after isz/zn pair has been updated
    }

    if (b_update_hidden_subzones && !hiddenSubzoneVec.empty()) {
      for (int i=0, end=hiddenSubzoneVec.size(); i<end; ++i) {
        const int global_isz = hiddenSubzoneVec[i];
        const int izone = znosz[global_isz];
        if ((izone<0) || (izone>int(zoneVec.size()))) continue;  // don't process subzones that seem to be improperly indexed

        const int local_isz = global_isz - szozn_i[izone];
        tracked_sz_local.push_back(pair<int,int> (izone,-(local_isz+1)));  // hidden are -1-indexed
      }
      hiddenSubzoneVec.clear(); // will repopulate after isz/zn pair has been updated
    }

    nsz = sz_flag.count();  // set new nsz
    IntFlag new_sz(sz_flag.getLength());  // map to new local subzone index
    new_sz.setAll(-1);
    // use zone_flag to store counts of subzones
    zone_flag.setLength(zoneVec.size());
    zone_flag.setAll(0);
    for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
      for (int isz=szozn_i[izn],end=szozn_i[izn+1]; isz<end; ++isz) {
        if (sz_flag[isz]) {
          new_sz[isz] = zone_flag[izn];  // stores new local sz index
          ++zone_flag[izn];
        }
      }
      // COUT2("szones in zone ["<<izn<<"]: " << zone_flag[izn]);
    }
    const int n_empty_zones = zone_flag.countZero();

    // update szost with local indexing
    for (int ist=0; ist<nst; ++ist) {
      szost[ist] = new_sz[szost[ist]];
      assert(szost[ist] != -1);
    }

    // update selected & highlighted local subzone indices
    for (it=tracked_sz_local.begin(); it!=tracked_sz_local.end(); ++it) {
      const int isz_i = it->second;
      const int isz_old = (isz_i < 0) ? (-(isz_i+1)):(isz_i-1);
      const int isz_new = new_sz[isz_old + szozn_i[it->first]];  // new_sz operates on old global indices

      if (isz_i < 0) it->second = -(isz_new+1);
      else if (isz_i > 0) it->second = isz_new+1;
      else {
        CWARN("unexpected subzone present (isz_i): " << isz_i);
      }
    }

    // build szozn_i from zone_flag (subzone counts)
    const int nzn = zoneVec.size()-n_empty_zones;
    szozn_i.resize(nzn+1);
    szozn_i[0] = 0;

    // and prune zoneVec if necessary
    if (n_empty_zones) {
      COUT2(" > found " << n_empty_zones << " empty zones; pruning");
      vector<SurfaceZone> zoneVec_old(zoneVec);
      zoneVec.clear();

      IntFlag new_zone(zone_flag.getLength());
      new_zone.setAll(-1);
      int zone_index = 0;
      for (int izone=0, limit=zone_flag.getLength(); izone < limit; ++izone) {
        if (zone_flag[izone] != 0) {
          szozn_i[zone_index+1] = szozn_i[zone_index] + zone_flag[izone];
          zoneVec.push_back(zoneVec_old[izone]);
          new_zone[izone] = zone_index;
          ++zone_index;
        }
        else {
          COUT2("    removing zone: " << zoneVec_old[izone].getName());
        }
      }
      assert(new_zone.countNegative() == n_empty_zones);
      assert(szozn_i[nzn] == nsz);

      // update znost
      for (int ist=0; ist<nst; ++ist) znost[ist] = new_zone[znost[ist]];

      // update zone info for selected/highlighted subzones
      for (it=tracked_sz_local.begin(); it!=tracked_sz_local.end(); ++it) it->first = new_zone[it->first];
    }
    else {
      // no zones pruned, so just adjust szozn_i
      for (int izone=0, limit=zone_flag.getLength(); izone < limit; ++izone) {
        szozn_i[izone+1] = szozn_i[izone] + zone_flag[izone];
      }
      assert(szozn_i[nzn] == nsz);
    }

    // adjust selected/highlighted subzones based on pruning
    for (it=tracked_sz_local.begin(); it!=tracked_sz_local.end(); ++it) {
      if ( (it->first >= 0) && (it->first<int(zoneVec.size())) ) {
        const int isz_i = it->second;
        const int isz = (isz_i < 0) ? (-(isz_i+1)):(isz_i-1);
        const int global_isz = isz + szozn_i[it->first];

        if (isz_i > 0) selectedSubzoneVec.push_back(global_isz);
        else if (isz_i < 0) hiddenSubzoneVec.push_back(global_isz);
        else {
          CWARN("unexpected subzone present (global_isz): " << global_isz);
        }
      }
      else {
        CWARN("pruning hidden/selected subzone (local sz,zn), no longer exists: " << it->second << "," << it->first);
      }
    }

    // update to global szost
    localSzostToGlobal();

    return true;
  }
  else return false;  // no pruning required

}

void SimpleSurface::flagNodesFromSubzoneVec(const vector<int>& subzonesVec,const bool b_exclude_boundaries) {
  assert(nsz == szozn_i[zoneVec.size()]);
  sz_flag.setLength(nsz);
  sz_flag.setAll(0);
  for (int i = 0, limit = subzonesVec.size(); i < limit; ++i) {
    const int isz = subzonesVec[i];
    if (isz < 0 || isz >= nsz) {
      CWARN("subzone index out-of-bounds: " << isz);
    }
    else {
      sz_flag[isz] = 1;
    }
  }

  sp_flag.setLength(nsp);
  sp_flag.setAll(0);
  for (int ist = 0; ist < nst; ++ist) {
    const int izn = znost[ist];
    assert((izn >= 0)&&(izn < int(zoneVec.size())));
    const int isz = szost[ist];
    assert((isz >= 0)&&(isz < nsz));
    assert((isz >= szozn_i[izn])&&(isz < szozn_i[izn+1]));
    if (sz_flag[isz]) {
      FOR_I3 sp_flag[spost[ist][i]] = 1;
    }
  }

  if (b_exclude_boundaries) {
    FOR_IST {
      const int izn = znost[ist];
      assert((izn >= 0)&&(izn < int(zoneVec.size())));
      const int isz = szost[ist];
      assert((isz >= 0)&&(isz < nsz));
      assert((isz >= szozn_i[izn])&&(isz < szozn_i[izn+1]));
      if (!sz_flag[isz])
        FOR_I3 sp_flag[spost[ist][i]] = 0;
    }
  }
}

void SimpleSurface::flagTrisFromSubzoneVec(const vector<int>& subzonesVec) {
  assert(nsz == szozn_i[zoneVec.size()]);
  sz_flag.setLength(nsz);
  sz_flag.setAll(0);
  for (int i = 0, limit = subzonesVec.size(); i < limit; ++i) {
    const int isz = subzonesVec[i];
    if (isz < 0 || isz >= nsz) {
      CWARN("subzone index out-of-bounds: " << isz);
    }
    else {
      sz_flag[isz] = 1;
    }
  }
  st_flag.setLength(nst);
  st_flag.setAll(0);
  for (int ist = 0; ist < nst; ++ist) {
    const int izn = znost[ist];
    assert((izn >= 0)&&(izn < int(zoneVec.size())));
    const int isz = szost[ist];
    assert((isz >= 0)&&(isz < nsz));
    assert((isz >= szozn_i[izn])&&(isz < szozn_i[izn+1]));
    if (sz_flag[isz])
      st_flag[ist] = 1;
  }
}

void SimpleSurface::flagNodesAndTrisFromSubzoneVec(const vector<int>& subzonesVec) {
  assert(nsz == szozn_i[zoneVec.size()]);
  sz_flag.setLength(nsz);
  sz_flag.setAll(0);
  for (int i = 0, limit = subzonesVec.size(); i < limit; ++i) {
    const int isz = subzonesVec[i];
    if (isz < 0 || isz >= nsz) {
      CWARN("subzone index out-of-bounds: " << isz);
    }
    else {
      sz_flag[isz] = 1;
    }
  }
  st_flag.setLength(nst);
  st_flag.setAll(0);
  sp_flag.setLength(nsp);
  sp_flag.setAll(0);
  for (int ist = 0; ist < nst; ++ist) {
    const int izn = znost[ist];
    assert((izn >= 0)&&(izn < int(zoneVec.size())));
    const int isz = szost[ist];
    assert((isz >= 0)&&(isz < nsz));
    assert((isz >= szozn_i[izn])&&(isz < szozn_i[izn+1]));
    if (sz_flag[isz]) {
      FOR_I3 sp_flag[spost[ist][i]] = 1;
      st_flag[ist] = 1;
    }
  }
}

void SimpleSurface::indexNodesAndTrisFromSubzoneVec(int& sp_count, int& st_count,const vector<int>& subzonesVec) {
  assert(nsz == szozn_i[zoneVec.size()]);
  sz_flag.setLength(nsz);
  sz_flag.setAll(0);
  for (vector<int>::const_iterator it=subzonesVec.begin(); it!=subzonesVec.end(); ++it) {
    const int isz = *it;
    if (isz < 0 || isz >= nsz) {
      CWARN("subzone index out-of-bounds: " << isz);
    }
    else {
      sz_flag[isz] = 1;
    }
  }

  sp_count = 0;
  sp_flag.setLength(nsp);
  sp_flag.setAll(-1);
  st_flag.setLength(nst);
  st_flag.setAll(-1);
  st_count = 0;
  for (int ist = 0; ist < nst; ++ist) {
    const int izn = znost[ist];
    assert((izn >= 0)&&(izn < int(zoneVec.size())));
    const int isz = szost[ist];
    assert((isz >= 0)&&(isz < nsz));
    assert((isz >= szozn_i[izn])&&(isz < szozn_i[izn+1]));
    if (sz_flag[isz]) {
      FOR_I3 {
        if (sp_flag[spost[ist][i]] == -1) sp_flag[spost[ist][i]] = sp_count++;
      }
      if (st_flag[ist] == -1) st_flag[ist] = st_count++;
    }
  }
}

void SimpleSurface::flagZonesFromZoneNameVec(const vector<string>& zoneNamesVec) {
  const int nzn=zoneVec.size();
  zone_flag.setLength(nzn);
  zone_flag.setAll(0);

  for (int i = 0, limit = zoneNamesVec.size(); i < limit; ++i) {
    int zone_index;
    for (zone_index=0; zone_index < nzn; ++zone_index) {
      if (zoneNamesVec[i] == zoneVec[zone_index].getName()) break;
    }
    if (zone_index != nzn) {
      zone_flag[zone_index] = 1;
    }
  }
}

void SimpleSurface::flagZonesFromZoneIndexVec(const vector<int>& zone_indices) {
  const int nzn=zoneVec.size();
  zone_flag.setLength(nzn);
  zone_flag.setAll(0);

  for (int i = 0, limit = zone_indices.size(); i < limit; ++i) {
    int zone_index = zone_indices[i];
    if (zone_index < nzn && zone_index >= 0) {
      zone_flag[zone_index] = 1;
    }
    else {
      WUI(WARN,"invalid zone index " << zone_index << " being skipped");
    }
  }
}

void SimpleSurface::selectSimilarSubzone(const int sz_index,const bool bSameZone,const bool byArea,
                                         const bool byNormal, const double diffTol, const double normalTol) {

  selectedSubzoneVec.clear();
  ensureSubzoneData();

  // reference zone data
  const int _zone = subzoneDataVec[sz_index].zone;
  assert(_zone != -1);  // subzone doesn't have specified parent
  const double _area = subzoneDataVec[sz_index].area;
  assert(_area > 0.0);
  const double _normal[3] = {subzoneDataVec[sz_index].normal[0],subzoneDataVec[sz_index].normal[1],subzoneDataVec[sz_index].normal[2]};

  // search through subzones for similarities
  for (int isz = 0, limit = subzoneDataVec.size(); isz < limit; ++isz) {
    if (isz == sz_index) {
      selectedSubzoneVec.push_back(isz);
    } else {
      // check similarity to reference subzone
      const bool sameParent = (subzoneDataVec[isz].zone == _zone) ? true:false;
      if (bSameZone && !sameParent) continue;

      const bool area_match = (fabs(_area - subzoneDataVec[isz].area) < diffTol*_area) ? true:false;
      if (byArea && !area_match) continue;

      const double dp = DOT_PRODUCT(_normal,subzoneDataVec[isz].normal);
      const double mag1 = MAG(_normal);
      const double mag2 = MAG(subzoneDataVec[isz].normal);
      //TODO user pass in this tolerance
      const bool normal_match = (dp > normalTol*mag1*mag2) ? true:false;
      if (byNormal && !normal_match) continue;

      selectedSubzoneVec.push_back(isz);
    }
  }
  const int new_zones = selectedSubzoneVec.size()-1;
  COUT1(" > similar subzones found: " << new_zones);
  if (!new_zones) {
    WUI(INFO,"No similar subzones were found");
  }
  else {
    WUI(INFO,new_zones << " additional subzones were selected");
  }
}

void SimpleSurface::selectSimilarOpenEdgeGroup(const int group_index,const bool byLength,const bool byAngleSum,const bool byCentroidMean, const double diffTol, const double angleSumTol, const double centroidMeanTol) {

  UNUSED(angleSumTol);
  UNUSED(byAngleSum);
  selectedSubzoneVec.clear();

  ensureOpenEdgeGroups();

  // reference group data
  const double _length = openEdgeGroupDataVec[group_index].length;
  assert(_length > 0.0);
  const double _mean_d = openEdgeGroupDataVec[group_index].mean_d;
  assert(_mean_d > 0.0);

  // search through subzones for similarities
  for (int igroup = 0, limit = openEdgeGroupDataVec.size(); igroup < limit; ++igroup) {
    if (igroup == group_index) {
      selectedSubzoneVec.push_back(igroup+32769);
    }
    else {
      // check similarity to reference group

      const bool length_match = (fabs(_length - openEdgeGroupDataVec[igroup].length) < diffTol*_length) ? true:false;
      if (byLength && !length_match) continue;

      const bool mean_d_match = (fabs(_mean_d - openEdgeGroupDataVec[igroup].mean_d) < centroidMeanTol*_mean_d) ? true:false;
      if (byCentroidMean && !mean_d_match) continue;

      selectedSubzoneVec.push_back(igroup+32769);
    }
  }

  const int new_groups = selectedSubzoneVec.size()-1;
  COUT1(" > similar edge groups found: " << new_groups);
  if (!new_groups) {
    WUI(INFO,"No similar edge groups were found");
  }
  else {
    WUI(INFO,new_groups << " additional edge groups were selected");
  }

}

void SimpleSurface::makeInjectorForFlaggedSubzones() {

  cout << "SimpleSurface::makeInjectorForFlaggedSubzones()" << endl;

  cout << " > trying noTeost version..." << endl;

  // if we don't have teost, use a different approach to compute
  // disjoint subzones that does not involve computing teost...

  vector<pair<int,int> > flagVec;
  {
    map<const pair<int,int>,int> edgeMap;
    for (int ist = 0; ist < nst; ++ist) {
      const int isz = szost[ist];
      if (sz_flag[isz]) {
        int count = flagVec.size();
        flagVec.push_back(pair<int,int>(ist,count));
        FOR_I3 {
          const int isp0 = spost[ist][i];
          const int isp1 = spost[ist][(i+1)%3];
          map<const pair<int,int>,int>::iterator iter = edgeMap.find(pair<int,int>(isp1,isp0));
          if (iter != edgeMap.end()) {
            int nbr_count = iter->second;
            edgeMap.erase(iter);
            while (nbr_count != flagVec[nbr_count].second)
              nbr_count = flagVec[nbr_count].second;
            if (nbr_count < count) {
              assert(flagVec[count].second == count);
              flagVec[count].second = nbr_count;
              count = nbr_count;
            }
            else if (count < nbr_count) {
              assert(flagVec[nbr_count].second == nbr_count);
              flagVec[nbr_count].second = count;
              nbr_count = count;
            }
          }
          else {
            edgeMap[pair<int,int>(isp0,isp1)] = count;
          }
        }
      }
    }
  }

  cout << " > nst flagged: " << flagVec.size() << endl;

  int ngr = 0;
  for (int ii = 0, limit = flagVec.size(); ii < limit; ++ii) {
    int count = flagVec[ii].second;
    while ((count >= 0)&&(count != flagVec[count].second))
      count = flagVec[count].second;
    if (count >= 0) {
      ++ngr;
      flagVec[ii].second = -ngr; // -1,-2,...
    }
    else {
      flagVec[ii].second = count;
    }
  }

  cout << " > ngr: " << ngr << endl;

  double (*n_sz)[3] = new double[ngr][3];
  double (*x_sz)[3] = new double[ngr][3];
  double *area_sz   = new double[ngr];

  for (int igr = 0; igr < ngr; ++igr) {
    FOR_I3 n_sz[igr][i] = 0.0;
    FOR_I3 x_sz[igr][i] = 0.0;
    area_sz[igr] = 0.0;
  }

  for (int ii = 0, limit = flagVec.size(); ii < limit; ++ii) {
    const int ist = flagVec[ii].first; assert((ist >= 0)&&(ist < nst));
    const int igr = -flagVec[ii].second-1; assert((igr >= 0)&&(igr < ngr));
    const double n2[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    FOR_I3 n_sz[igr][i] += n2[i];
    const double area = MAG(n2);
    area_sz[igr] += area;
    FOR_I3 x_sz[igr][i] += area*(xsp[spost[ist][0]][i] + xsp[spost[ist][1]][i] + xsp[spost[ist][2]][i]); // don't forget about the factor of 3
  }

  for (int igr = 0; igr < ngr; ++igr) {
    FOR_I3 n_sz[igr][i] *= 0.5;
    FOR_I3 x_sz[igr][i] /= area_sz[igr]*3.0;
    area_sz[igr] *= 0.5;
    const double n_mag = MAG(n_sz[igr]);
    if (n_mag/area_sz[igr] > 0.95) {
      // this is a planar injector, so add injection based on computed area...
      cout << "this is planar: " << endl;
    }
    else {
      // for now, we assume this group represents a cylinder. Try to figure out its
      // orientation and radius...
      cout << "this is a cylinder?: " << endl;
      st_flag.resize(nst);
      st_flag.setAll(0);
      for (int ii = 0, limit = flagVec.size(); ii < limit; ++ii) {
        const int ist = flagVec[ii].first; assert((ist >= 0)&&(ist < nst));
        const int this_igr = -flagVec[ii].second-1; assert((igr >= 0)&&(igr < ngr));
        if (this_igr == igr)
          st_flag[ist] = 1;
      }
      double xc[3],nc[3],r;
      calcClosestCylinderForFlaggedTris(xc,nc,r);
    }
    // report...
    cout << " > group: " << igr << " normal: " << COUT_VEC(n_sz[igr]) << " centroid: " << COUT_VEC(x_sz[igr]) << " area: " << area_sz[igr] << " nmag/area: " << n_mag/area_sz[igr] << endl;

  }

  delete[] n_sz;
  delete[] x_sz;
  delete[] area_sz;

  cout << " > done noTeost version." << endl;

}

void calcHessianGrad(double hessian[3][3], double grad[3], double xi[2], const double& wi, double xc[2], double& r) {
  double dx[2] = DIFF_2D(xi,xc);
  double temp = dx[0]*dx[0] + dx[1]*dx[1] - r*r;
  FOR_I2 {
    grad[i] = -4.0 * temp * dx[i];// * wi;
  }
  grad[2] = -4.0 * temp * r;// * wi;

  hessian[0][0] = 8.0 * dx[0]*dx[0] + 4.0 * temp;
  hessian[1][1] = 8.0 * dx[1]*dx[1] + 4.0 * temp;
  hessian[2][2] = 8.0 * r*r - 4.0 * temp;
  hessian[1][0] = hessian[0][1] = 8.0 * dx[0] * dx[1];
  hessian[2][0] = hessian[0][2] = 8.0 * r * dx[0];
  hessian[2][1] = hessian[1][2] = 8.0 * r * dx[1];
}

double calcCostGrad(double grad[3], const double xi[2], const double& wi, double xc[2], double& r) {
  double dx[2] = DIFF_2D(xi,xc);
  double temp = dx[0]*dx[0] + dx[1]*dx[1] - r*r;
  double cost = temp*temp;// * wi;
  FOR_I2 {
    grad[i] = -4.0 * temp * dx[i];// * wi;
  }
  grad[2] = -4.0 * temp * r;// * wi;
  return cost;
}

void SimpleSurface::characterizePipe(double &l_on_d,double theta[2],const int isz) {

  vector<pair<int,int> > flagVec;
  map<const pair<int,int>,pair<int,int> > edgeMap;
  for (int toz = stosz_i[isz]; toz != stosz_i[isz+1]; ++toz) {
    const int ist = stosz_v[toz];
    int count = flagVec.size();
    flagVec.push_back(pair<int,int>(ist,count));
    FOR_I3 {
      const int isp0 = spost[ist][i];
      const int isp1 = spost[ist][(i+1)%3];
      map<const pair<int,int>,pair<int,int> >::iterator iter = edgeMap.find(pair<int,int>(isp1,isp0));
      if (iter != edgeMap.end()) {
        int nbr_count = iter->second.first;
        edgeMap.erase(iter);
        while (nbr_count != flagVec[nbr_count].second)
          nbr_count = flagVec[nbr_count].second;
        if (nbr_count < count) {
          assert(flagVec[count].second == count);
          flagVec[count].second = nbr_count;
          count = nbr_count;
        }
        else if (count < nbr_count) {
          assert(flagVec[nbr_count].second == nbr_count);
          flagVec[nbr_count].second = count;
          nbr_count = count;
        }
      }
      else {
        edgeMap[pair<int,int>(isp0,isp1)] = pair<int,int>(count,ist);
      }
    }
  }

  cout << " > nst flagged: " << flagVec.size() << endl;

  int ngr = 0;
  for (int ii = 0, limit = flagVec.size(); ii < limit; ++ii) {
    int count = flagVec[ii].second;
    while ((count >= 0)&&(count != flagVec[count].second))
      count = flagVec[count].second;
    if (count >= 0) {
      ++ngr;
      flagVec[ii].second = -ngr; // -1,-2,...
    }
    else {
      flagVec[ii].second = count;
    }
  }

  cout << " > ngr: " << ngr << endl;

  // for the following code, we need ngr to be 1 for now...
  assert(ngr == 1);

  for (int igr = 0; igr < ngr; ++igr) {

    cout << " > might be a cylinder..." << endl;
    // flag everybody...
    st_flag.resize(nst);
    st_flag.setAll(0);
    for (int ii = 0, limit = flagVec.size(); ii < limit; ++ii) {
      const int ist = flagVec[ii].first; assert((ist >= 0)&&(ist < nst));
      const int this_igr = -flagVec[ii].second-1; assert((igr >= 0)&&(igr < ngr));
      if (this_igr == igr)
        st_flag[ist] = 1;
    }
    // now find the edge loops associated with the ends...
    map<const int,int> linkMap;
    for (map<const pair<int,int>,pair<int,int> >::iterator iter = edgeMap.begin(); iter != edgeMap.end(); ++iter) {
      const int ist = iter->second.second;
      if (st_flag[ist] == 1) {
        const int isp0 = iter->first.first;
        const int isp1 = iter->first.second;
        linkMap.insert(pair<int,int>(isp0,isp1));
      }
    }

    // how many loops...
    int nlink = linkMap.size();
    int * spoli_v = new int[nlink];

    vector<double> doubleVec;
    vector<int> spoli_i;
    int ilink = 0;
    while (!linkMap.empty()) {

      double xc[3] = { 0.0, 0.0, 0.0 };
      double wgt = 0.0;

      map<const int,int>::iterator iter = linkMap.begin();
      assert(iter != linkMap.end());
      const int start = iter->first;
      int next = start;
      do {

        iter = linkMap.find(next);
        assert(iter != linkMap.end());
        const int prev = next;
        next = iter->second;
        linkMap.erase(iter);

        spoli_v[ilink++] = prev;

        // add segment prev->next to xc,wgt...
        const double this_wgt = DIST(xsp[next],xsp[prev]);
        FOR_I3 xc[i] += this_wgt*(xsp[prev][i]+xsp[next][i]);
        wgt += this_wgt; // don't forget factor of 2 when normalizing!

      } while (next != start);

      FOR_I3 {
        xc[i] /= 2.0*wgt;
        doubleVec.push_back(xc[i]);
      }

      cout << " > link: " << spoli_i.size() << " xc: " << COUT_VEC(xc) << endl;

      spoli_i.push_back(ilink);

    }
    assert(ilink == nlink);

    // we expect 2 links. If we don't then do something else eventually...
    assert(spoli_i.size() == 2);
    const double x0[3] = { doubleVec[0], doubleVec[1], doubleVec[2] };
    const double x1[3] = { doubleVec[3], doubleVec[4], doubleVec[5] };

    double e0[3] = DIFF(x1,x0);
    double mag_e0 = MAG(e0);
    FOR_I3 e0[i] /= mag_e0;

    double e1[3],e2[3];
    MiscUtils::getBestE1E2FromE0(e1,e2,e0);

    const double L =
    (x1[0]-x0[0])*e0[0] +
    (x1[1]-x0[1])*e0[1] +
    (x1[2]-x0[2])*e0[2];
    cout << " > L: " << L << endl;

    // and figure out the radius from the tris...
    double r_avg = 0.0;
    double r_min = HUGE_VAL;
    double r_max = 0.0;
    double wgt = 0.0;
    for (int toz = stosz_i[isz]; toz != stosz_i[isz+1]; ++toz) {
      const int ist = stosz_v[toz];
      const double normal[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
      const double this_wgt = MAG(normal);
      FOR_I3 {
        const int isp = spost[ist][i];
        const double dx[3] = DIFF(xsp[isp],x0);
        const double yp = DOT_PRODUCT(dx,e1);
        const double zp = DOT_PRODUCT(dx,e2);
        const double rp = sqrt(yp*yp + zp*zp);
        wgt += this_wgt;
        r_avg += this_wgt*rp;
        r_min = min(r_min,rp);
        r_max = max(r_max,rp);
      }
      // the minimum radius for reasonable tris will probably occur along an edge...
      int isp1 = spost[ist][2];
      FOR_I3 {
        const int isp0 = isp1;
        isp1 = spost[ist][i];
        const double dx[3] = {
          0.5*(xsp[isp0][0]+xsp[isp1][0]) - x0[0],
          0.5*(xsp[isp0][1]+xsp[isp1][1]) - x0[1],
          0.5*(xsp[isp0][2]+xsp[isp1][2]) - x0[2] };
        const double yp = DOT_PRODUCT(dx,e1);
        const double zp = DOT_PRODUCT(dx,e2);
        const double rp = sqrt(yp*yp + zp*zp);
        r_min = min(r_min,rp);
        r_max = max(r_max,rp);
      }
    }
    r_avg /= wgt;
    cout << " > r_avg: " << r_avg << " plus/minus: " << r_min-r_avg << " " << r_max-r_avg << endl;

    l_on_d = 0.5*L/r_avg;

    // now check the length by projecting both links in e0...
    double xp_min = 0.0;
    double xp_max = 0.0;
    for (int ilink = 0; ilink < spoli_i[0]; ++ilink) {
      const int isp = spoli_v[ilink];
      double dx[3] = DIFF(xsp[isp],x0);
      double xp = DOT_PRODUCT(dx,e0);
      xp_min = min(xp_min,xp);
      xp_max = max(xp_max,xp);
    }
    cout << " > cylinder starts at 0 plus/minus: " << xp_min << " " << xp_max << endl;

    theta[0] = 0.5*M_PI-atan2(2.0*r_avg,xp_max-xp_min);

    // second loop...
    xp_min = L;
    xp_max = L;
    for (int ilink = spoli_i[0]; ilink < spoli_i[1]; ++ilink) {
      const int isp = spoli_v[ilink];
      double dx[3] = DIFF(xsp[isp],x0);
      double xp = DOT_PRODUCT(dx,e0);
      xp_min = min(xp_min,xp);
      xp_max = max(xp_max,xp);
    }
    cout << " > cylinder ends at L plus/minus: " << xp_min-L << " " << xp_max-L << endl;

    theta[1] = 0.5*M_PI-atan2(2.0*r_avg,xp_max-xp_min);

  }

}

void SimpleSurface::calcClosestCylinderForFlaggedTris(double xc[3],double nc[3],double& r) {

  //writeSelectedFacesByZoneToTecplot("test.dat",st_flag);
  //cout << "TAKE A LOOK" << endl;

  double LS_mat[3][3] = { { 0.0, 0.0, 0.0 },
			  { 0.0, 0.0, 0.0 },
			  { 0.0, 0.0, 0.0 } };

  for (int ist = 0; ist < nst; ist++) {
    if (st_flag[ist] == 1) {
      double unit_n_st[3];
      const double * const x0 = xsp[spost[ist][0]];
      const double * const x1 = xsp[spost[ist][1]];
      const double * const x2 = xsp[spost[ist][2]];
      const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
      const double mag = MAG(this_n);
      assert(mag > 0.0);
      FOR_I3 unit_n_st[i] = this_n[i]/mag;

      for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
          LS_mat[i][j] += unit_n_st[i] * unit_n_st[j];
        }
      }

    }
  }

  double eig_vec[3][3];
  double eig_val[3];
  MiscUtils::eigenDecomposition(LS_mat, eig_vec, eig_val);
  FOR_I3 cout << "eig_val i: " << i << " " << eig_val[i] << endl;
  FOR_I3
    FOR_J3
    cout << "eig_vec[" << i << "][" << j << "]: "<< eig_vec[i][j] << endl;

  // cylinder normal
  FOR_I3 nc[i] = eig_vec[i][0];

  // cylinder other axes
  double uc[3]; FOR_I3 uc[i] = eig_vec[i][1];
  double vc[3] = CROSS_PRODUCT(nc,uc);
  assert(MAG(nc) - 1.0 < 1e-6);
  assert(MAG(uc) - 1.0 < 1e-6);
  assert(MAG(vc) - 1.0 < 1e-6);
  assert(DOT_PRODUCT(nc,uc) < 1e-6);
  assert(DOT_PRODUCT(uc,vc) < 1e-6);
  assert(DOT_PRODUCT(vc,nc) < 1e-6);
  cout << "nc: " << COUT_VEC(nc) << " ,uc: " << COUT_VEC(uc) << " ,vc: " << COUT_VEC(vc) << endl;

  // use gradient descent to find the center and radius of cylinder
  double xc_cur[2] = { 0.0, 0.0 };
  double sum_wgt = 0.0;

  for (int ist = 0; ist < nst; ist++) {
    if (st_flag[ist] == 1){
      const double * const x0 = xsp[spost[ist][0]];
      const double * const x1 = xsp[spost[ist][1]];
      const double * const x2 = xsp[spost[ist][2]];
      const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
      const double area = MAG(this_n);
      assert(area > 0.0);
      double x0_2d[2] = {DOT_PRODUCT(x0,uc), DOT_PRODUCT(x0,vc)};
      FOR_I2 xc_cur[i] += area*x0_2d[i];
      double x1_2d[2] = {DOT_PRODUCT(x1,uc), DOT_PRODUCT(x1,vc)};
      FOR_I2 xc_cur[i] += area*x1_2d[i];
      double x2_2d[2] = {DOT_PRODUCT(x2,uc), DOT_PRODUCT(x2,vc)};
      FOR_I2 xc_cur[i] += area*x2_2d[i];
      sum_wgt += area;
    }
  }
  FOR_I2 xc_cur[i] /= sum_wgt*3.0;

  // for the radius, just take the mean radius fromthe guess...
  double r_cur = 0.0;
  sum_wgt = 0.0;
  for (int ist = 0; ist < nst; ist++) {
    if (st_flag[ist] == 1){
      const double * const x0 = xsp[spost[ist][0]];
      const double * const x1 = xsp[spost[ist][1]];
      const double * const x2 = xsp[spost[ist][2]];
      const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
      const double area = MAG(this_n);
      assert(area > 0.0);
      double x0_2d[2] = {DOT_PRODUCT(x0,uc), DOT_PRODUCT(x0,vc)};
      r_cur += area*DIST_2D(xc_cur,x0_2d);
      double x1_2d[2] = {DOT_PRODUCT(x1,uc), DOT_PRODUCT(x1,vc)};
      r_cur += area*DIST_2D(xc_cur,x1_2d);
      double x2_2d[2] = {DOT_PRODUCT(x2,uc), DOT_PRODUCT(x2,vc)};
      r_cur += area*DIST_2D(xc_cur,x2_2d);
      sum_wgt += area;
    }
  }
  r_cur /= sum_wgt*3.0;
  assert(r_cur > 0.0);

  cout << "initial guess: xc_cur: " << COUT_VEC_2D(xc_cur) << " r_cur: " << r_cur << endl;

  enum {
    NEWTON,
    GRAD_DES
  } OPT_METHOD;

  OPT_METHOD = NEWTON;
  //OPT_METHOD = GRAD_DES;

  switch (OPT_METHOD) {
  case NEWTON:
    // Newton's method optimization
    {
      int iter_max = 1000;
      int iter = 1;
      bool done = false;
      double tol = 5e-16;
      while ((!done) && (iter<iter_max)) {
        int num_st = 0;
        double Hessian_sum[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        double Grad_sum[3] = {0.0, 0.0, 0.0};
        for (int ist = 0; ist < nst; ist++) {
          if (st_flag[ist] == 1){
            const double * const x0 = xsp[spost[ist][0]];
            const double * const x1 = xsp[spost[ist][1]];
            const double * const x2 = xsp[spost[ist][2]];
            const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
            const double area = MAG(this_n);
            assert(area > 0.0);
            double Hessian_cur[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
            double Grad_cur[3] = {0.0, 0.0, 0.0};
            double x0_2d [2] = {DOT_PRODUCT(x0,uc), DOT_PRODUCT(x0,vc)};
            calcHessianGrad(Hessian_cur, Grad_cur, x0_2d, area, xc_cur, r_cur);
            FOR_I3 Grad_sum[i] += Grad_cur[i];
            FOR_I3 {
              FOR_J3 {
                Hessian_sum[i][j] += Hessian_cur[i][j];
              }
            }
            double x1_2d [2] = {DOT_PRODUCT(x1,uc), DOT_PRODUCT(x1,vc)};
            calcHessianGrad(Hessian_cur, Grad_cur, x1_2d, area, xc_cur, r_cur);
            FOR_I3 Grad_sum[i] += Grad_cur[i];
            FOR_I3 {
              FOR_J3 {
                Hessian_sum[i][j] += Hessian_cur[i][j];
              }
            }
            double x2_2d [2] = {DOT_PRODUCT(x2,uc), DOT_PRODUCT(x2,vc)};
            calcHessianGrad(Hessian_cur, Grad_cur, x2_2d, area, xc_cur, r_cur);
            FOR_I3 Grad_sum[i] += Grad_cur[i];
            FOR_I3 {
              FOR_J3 {
                Hessian_sum[i][j] += Hessian_cur[i][j];
              }
            }
            num_st++;
          }
        }
        assert(num_st > 0);
        FOR_I3 Grad_sum[i] /= num_st*3.0;
        FOR_I3 {
          FOR_J3 {
            Hessian_sum[i][j] /= num_st*3.0;
          }
        }

        double inv_Hessian[3][3];
        MiscUtils::invertMat(inv_Hessian, Hessian_sum);

        double dx[3] = MATRIX_PRODUCT(inv_Hessian, Grad_sum);

        FOR_I2 {
          xc_cur[i] -= dx[i];
        }
        double r_old = r_cur;
        r_cur -= dx[2];

        double r_update = abs(r_cur - r_old) / r_old;
        if (r_update < tol) {
          done = true;
        }

        iter++;
        cout << "iter: " << iter <<  " error: " << r_update << " xc_cur: " << COUT_VEC_2D(xc_cur) << " r_cur: " << r_cur << endl;

      }
    }
    break;

  case GRAD_DES:
    // Gradient descent optimization
    {
      int iter_max = 100000;
      int iter = 1;
      bool done = false;
      double tol = 5e-16;
      double lambda = 1e-1;
      int lam_inc_frq = 100;
      double lam_inc_fac = 1.05;
      // momentum parameters
      double mu = 0.9;
      double v[3] = {0.0, 0.0, 0.0};
      // Adagrad parameters
      //double cache[3] = {0.0, 0.0, 0.0};
      //double eps = 1e-8;
      while ((!done) && (iter<iter_max)) {
        cout << "here" << endl;
        // gradients for xc_cur and r_cur
        double grad_sum[3] = {0.0, 0.0, 0.0};
        double cost_sum = 0.0;
        int num_st = 0;
        // optional: can increase lambda for larger iterations
        if (iter%lam_inc_frq == 0) lambda *= lam_inc_fac;
        for (int ist = 0; ist < nst; ist++) {
          if (st_flag[ist] == 1){
            const double * const x0 = xsp[spost[ist][0]];
            const double * const x1 = xsp[spost[ist][1]];
            const double * const x2 = xsp[spost[ist][2]];
            const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
            const double area = MAG(this_n);
            assert(area > 0.0);
            double grad_cur[3] = {0.0, 0.0, 0.0};
            double cost_cur = 0.0;
            double x0_2d [2] = {DOT_PRODUCT(x0,uc), DOT_PRODUCT(x0,vc)};
            cost_cur = calcCostGrad(grad_cur, x0_2d, area, xc_cur, r_cur);
            cost_sum += cost_cur;
            FOR_I3 grad_sum[i] += grad_cur[i];
            double x1_2d [2] = {DOT_PRODUCT(x1,uc), DOT_PRODUCT(x1,vc)};
            cost_cur = calcCostGrad(grad_cur, x1_2d, area, xc_cur, r_cur);
            cost_sum += cost_cur;
            FOR_I3 grad_sum[i] += grad_cur[i];
            double x2_2d [2] = {DOT_PRODUCT(x2,uc), DOT_PRODUCT(x2,vc)};
            cost_cur = calcCostGrad(grad_cur, x2_2d, area, xc_cur, r_cur);
            cost_sum += cost_cur;
            FOR_I3 grad_sum[i] += grad_cur[i];
            num_st++;
          }
        }
        assert(num_st > 0);
        FOR_I3 grad_sum[i] /= num_st*3.0;
        cost_sum /= num_st*3.0;

        // naive update
        //FOR_I2 xc_cur[i] -= lambda * grad_sum[i];
        //r_cur -= lambda * grad_sum[2];
        // momentum update
        FOR_I3 v[i] = mu * v[i] - lambda * grad_sum[i];
        FOR_I2 xc_cur[i] += v[i];
        r_cur += v[2];
        // Adagrad update
        //FOR_I3 cache[i] += grad_sum[i] * grad_sum[i];
        //FOR_I2 xc_cur[i] -= lambda * grad_sum[i] / (sqrt(cache[i]) + eps);
        //r_cur -= lambda * grad_sum[2] / (sqrt(cache[2]) + eps);
        cout << "iter: " << iter << " cost_sum: " << cost_sum << " xc_cur: " << COUT_VEC_2D(xc_cur) << " r_cur: " << r_cur << " lambda: " << lambda << endl;
        if (cost_sum < tol) {
          done = true;
        }
        iter++;
      }
    }
    break;
  } // switch

  // compute the nc component for xc
  double area_sum = 0.0;
  double xc_n = 0.0;
  for (int ist = 0; ist < nst; ist++) {
    if (st_flag[ist] == 1){
      const double * const x0 = xsp[spost[ist][0]];
      const double * const x1 = xsp[spost[ist][1]];
      const double * const x2 = xsp[spost[ist][2]];
      const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
      const double area = MAG(this_n);
      assert(area > 0.0);
      double xi[3];
      FOR_I3 xi[i] = (x0[i] + x1[i] + x2[i])/3.0;
      xc_n += DOT_PRODUCT(xi,nc) * area;
      area_sum += area;
    }
  }
  xc_n /= area_sum;
  r = r_cur;
  FOR_I3 xc[i] = xc_cur[0] * uc[i] + xc_cur[1] * vc[i] + xc_n * nc[i];
  cout << "injection location xc: " << COUT_VEC(xc) << endl;

  // find the neighboring zone's normal vectors
  // The injection normal should be oposite to this direction
  vector<int> zoneVec;
  map<int, double*> normalMap;
  for (int ist = 0; ist < nst; ist++) {
    if (st_flag[ist] == 1) {
      int ist_nbr[3];
      getTriNbrs(ist_nbr, ist);
      FOR_I3 {
        int sz_nbr = szost[ist_nbr[i]];
        int sz_self = szost[ist];
        if (sz_nbr != sz_self) {
          // compute this triangle which should be added to the map
          const double * const x0 = xsp[spost[ist_nbr[i]][0]];
          const double * const x1 = xsp[spost[ist_nbr[i]][1]];
          const double * const x2 = xsp[spost[ist_nbr[i]][2]];
          double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
          const double mag = MAG(this_n);
          assert(mag > 0.0);
          FOR_J3 this_n[j] = this_n[j]/mag;

          map<int, double*>::iterator iter = normalMap.find(sz_nbr);
          if (iter != normalMap.end()) {
            // add the normal to the existing vecotor in the map
            FOR_J3 normalMap[sz_nbr][j] += this_n[j];
          }
          else {
            // create a new vector in the map and add the normal to it
            double* st_norm = new double[3];
            FOR_J3 st_norm[j] = this_n[j];
            normalMap[sz_nbr] = st_norm;
            zoneVec.push_back(sz_nbr);
          }
        }
      }
    }
  }

  for (int i = 0, nzn=zoneVec.size(); i < nzn; i++) {
    int isz = zoneVec[i];
    assert(normalMap.find(isz) != normalMap.end());
    double mag = MAG(normalMap[isz]);
    FOR_J3 normalMap[isz][j] /= mag;
  }

  // average neighbor's normal
  double normal_avg[3] = {0.0, 0.0, 0.0};
  for (int i = 0,nzn=zoneVec.size(); i < nzn; i++) {
    int isz = zoneVec[i];
    assert(normalMap.find(isz) != normalMap.end());
    FOR_J3 normal_avg[j] += normalMap[isz][j];
  }
  double mag = MAG(normal_avg);
  if (mag > 0) {
    FOR_J3 normal_avg[j] /= mag;
  }
  else {
    cout << "WARNING: there is no preferencial direction for injection... " << endl;
    cout << "The best injection direction proposed is nc: " << COUT_VEC(nc) << endl;
    for (int i = 0,nzn=zoneVec.size(); i < nzn; i++) {
      delete [] normalMap[zoneVec[i]];
    }

    return;
  }

  double sign = DOT_PRODUCT(normal_avg, nc);
  sign /= abs(sign);
  //cout << "normal_avg: " << COUT_VEC(normal_avg) << " sign: " << sign << endl;

  FOR_I3 nc[i] *= -sign;

  cout << "injection normal: " << COUT_VEC(nc) << endl;

  for (int i = 0,nzn=zoneVec.size(); i < nzn; i++) {
    delete [] normalMap[zoneVec[i]];
  }

}

void SimpleSurface::calcClosestCylinderForFlaggedTris(double nc[3], double xc[3]) {

  //writeSelectedFacesByZoneToTecplot("test.dat",st_flag);
  //cout << "TAKE A LOOK" << endl;

  FOR_I3 xc[i] = 0.0;
  double area_sum = 0.0;
  double LS_mat[3][3] = { { 0.0, 0.0, 0.0 },
                          { 0.0, 0.0, 0.0 },
                          { 0.0, 0.0, 0.0 } };

  for (int ist = 0; ist < nst; ist++) {
    if (st_flag[ist] == 1) {
      double unit_n_st[3];
      const double * const x0 = xsp[spost[ist][0]];
      const double * const x1 = xsp[spost[ist][1]];
      const double * const x2 = xsp[spost[ist][2]];
      const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
      const double mag = MAG(this_n);
      assert(mag > 0.0);
      FOR_I3 unit_n_st[i] = this_n[i]/mag;

      FOR_I3 xc[i] += (x0[i]+x1[i]+x2[i])*mag;
      area_sum += mag;

      for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
          LS_mat[i][j] += unit_n_st[i] * unit_n_st[j];
        }
      }

    }
  }

  FOR_I3 xc[i] /= (area_sum*3.);

  double eig_vec[3][3];
  double eig_val[3];
  MiscUtils::eigenDecomposition(LS_mat, eig_vec, eig_val);
  FOR_I3 {
    COUT2(" > eig_val i: " << i << " " << eig_val[i]);
  }

  // cylinder normal
  FOR_I3 nc[i] = eig_vec[i][0];

  // cylinder other axes
  double uc[3]; FOR_I3 uc[i] = eig_vec[i][1];
  double vc[3] = CROSS_PRODUCT(nc,uc);
  assert(MAG(nc) - 1.0 < 1e-6);
  assert(MAG(uc) - 1.0 < 1e-6);
  assert(MAG(vc) - 1.0 < 1e-6);
  assert(DOT_PRODUCT(nc,uc) < 1e-6);
  assert(DOT_PRODUCT(uc,vc) < 1e-6);
  assert(DOT_PRODUCT(vc,nc) < 1e-6);
  cout << " > nc: " << COUT_VEC(nc) << " ,uc: " << COUT_VEC(uc) << " ,vc: " << COUT_VEC(vc) << endl;

}

void SimpleSurface::makePartBlForFlaggedSubzones(const string& name,const double dn,const double dt,const int n) {

  cout << "SimpleSurface::makePartBlForFlaggedSubzones: name: " << name << " dn: " << dn << " dt: " << dt << " n: " << n << endl;

  vector<pair<int,int> > flagVec;
  map<const pair<int,int>,pair<int,int> > edgeMap;
  for (int ist = 0; ist < nst; ++ist) {
    const int isz = szost[ist];
    if (sz_flag[isz]) {
      int count = flagVec.size();
      flagVec.push_back(pair<int,int>(ist,count));
      FOR_I3 {
        const int isp0 = spost[ist][i];
        const int isp1 = spost[ist][(i+1)%3];
        map<const pair<int,int>,pair<int,int> >::iterator iter = edgeMap.find(pair<int,int>(isp1,isp0));
        if (iter != edgeMap.end()) {
          int nbr_count = iter->second.first;
          edgeMap.erase(iter);
          while (nbr_count != flagVec[nbr_count].second)
            nbr_count = flagVec[nbr_count].second;
          if (nbr_count < count) {
            assert(flagVec[count].second == count);
            flagVec[count].second = nbr_count;
            count = nbr_count;
          }
          else if (count < nbr_count) {
            assert(flagVec[nbr_count].second == nbr_count);
            flagVec[nbr_count].second = count;
            nbr_count = count;
          }
        }
        else {
          edgeMap[pair<int,int>(isp0,isp1)] = pair<int,int>(count,ist);
        }
      }
    }
  }

  cout << " > nst flagged: " << flagVec.size() << endl;

  int ngr = 0;
  for (int ii = 0, limit = flagVec.size(); ii < limit; ++ii) {
    int count = flagVec[ii].second;
    while ((count >= 0)&&(count != flagVec[count].second))
      count = flagVec[count].second;
    if (count >= 0) {
      ++ngr;
      flagVec[ii].second = -ngr; // -1,-2,...
    }
    else {
      flagVec[ii].second = count;
    }
  }

  cout << " > ngr: " << ngr << endl;

  double (*n_sz)[3] = new double[ngr][3];
  double (*x_sz)[3] = new double[ngr][3];
  double *area_sz   = new double[ngr];

  for (int igr = 0; igr < ngr; ++igr) {
    FOR_I3 n_sz[igr][i] = 0.0;
    FOR_I3 x_sz[igr][i] = 0.0;
    area_sz[igr] = 0.0;
  }

  for (int ii = 0, limit = flagVec.size(); ii < limit; ++ii) {
    const int ist = flagVec[ii].first; assert((ist >= 0)&&(ist < nst));
    const int igr = -flagVec[ii].second-1; assert((igr >= 0)&&(igr < ngr));
    const double n2[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    FOR_I3 n_sz[igr][i] += n2[i];
    const double area = MAG(n2);
    area_sz[igr] += area;
    FOR_I3 x_sz[igr][i] += area*(xsp[spost[ist][0]][i] + xsp[spost[ist][1]][i] + xsp[spost[ist][2]][i]); // don't forget about the factor of 3
  }

  // for the following code, we need ngr to be 1 for now...
  assert(ngr == 1);

  for (int igr = 0; igr < ngr; ++igr) {

    cout << "working on group: " << igr << endl;

    FOR_I3 n_sz[igr][i] *= 0.5;
    FOR_I3 x_sz[igr][i] /= area_sz[igr]*3.0;
    area_sz[igr] *= 0.5;
    const double n_mag = MAG(n_sz[igr]);
    if (n_mag/area_sz[igr] > 0.95) {
      cout << " > seems planar: skipping." << endl;
    }
    else {
      cout << " > might be a cylinder..." << endl;
      // flag everybody...
      st_flag.resize(nst);
      st_flag.setAll(0);
      for (int ii = 0, limit = flagVec.size(); ii < limit; ++ii) {
        const int ist = flagVec[ii].first; assert((ist >= 0)&&(ist < nst));
        const int this_igr = -flagVec[ii].second-1; assert((igr >= 0)&&(igr < ngr));
        if (this_igr == igr)
          st_flag[ist] = 1;
      }
      // now find the edge loops associated with the ends...
      map<const int,int> linkMap;
      for (map<const pair<int,int>,pair<int,int> >::iterator iter = edgeMap.begin(); iter != edgeMap.end(); ++iter) {
        const int ist = iter->second.second;
        if (st_flag[ist] == 1) {
          const int isp0 = iter->first.first;
          const int isp1 = iter->first.second;
          linkMap.insert(pair<int,int>(isp0,isp1));
        }
      }

      // how many loops...
      int nlink = linkMap.size();
      int * spoli_v = new int[nlink];

      vector<double> doubleVec;
      vector<int> spoli_i;
      int ilink = 0;
      while (!linkMap.empty()) {

        double xc[3] = { 0.0, 0.0, 0.0 };
        double wgt = 0.0;

        map<const int,int>::iterator iter = linkMap.begin();
        assert(iter != linkMap.end());
        const int start = iter->first;
        int next = start;
        do {

          iter = linkMap.find(next);
          assert(iter != linkMap.end());
          const int prev = next;
          next = iter->second;
          linkMap.erase(iter);

          spoli_v[ilink++] = prev;

          // add segment prev->next to xc,wgt...
          const double this_wgt = DIST(xsp[next],xsp[prev]);
          FOR_I3 xc[i] += this_wgt*(xsp[prev][i]+xsp[next][i]);
          wgt += this_wgt; // don't forget factor of 2 when normalizing!

        } while (next != start);

        FOR_I3 {
          xc[i] /= 2.0*wgt;
          doubleVec.push_back(xc[i]);
        }

        cout << " > link: " << spoli_i.size() << " xc: " << COUT_VEC(xc) << endl;

        spoli_i.push_back(ilink);

      }
      assert(ilink == nlink);

      // we expect 2 links. If we don't then do something else eventually...
      assert(spoli_i.size() == 2);
      const double x0[3] = { doubleVec[0], doubleVec[1], doubleVec[2] };
      const double x1[3] = { doubleVec[3], doubleVec[4], doubleVec[5] };

      double e0[3] = DIFF(x1,x0);
      double mag_e0 = MAG(e0);
      FOR_I3 e0[i] /= mag_e0;

      double e1[3],e2[3];
      MiscUtils::getBestE1E2FromE0(e1,e2,e0);

      const double L =
        (x1[0]-x0[0])*e0[0] +
        (x1[1]-x0[1])*e0[1] +
        (x1[2]-x0[2])*e0[2];
      cout << " > L: " << L << endl;

      // and figure out the radius from the tris...
      double r_avg = 0.0;
      double r_min = HUGE_VAL;
      double r_max = 0.0;
      double wgt = 0.0;
      for (int ist = 0; ist < nst; ++ist) {
        if (st_flag[ist] == 1) {
          const double normal[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
          const double this_wgt = MAG(normal);
          FOR_I3 {
            const int isp = spost[ist][i];
            const double dx[3] = DIFF(xsp[isp],x0);
            const double yp = DOT_PRODUCT(dx,e1);
            const double zp = DOT_PRODUCT(dx,e2);
            const double rp = sqrt(yp*yp + zp*zp);
            wgt += this_wgt;
            r_avg += this_wgt*rp;
            r_min = min(r_min,rp);
            r_max = max(r_max,rp);
          }
          // the minimum radius for reasonable tris will probably occur along an edge...
          int isp1 = spost[ist][2];
          FOR_I3 {
            const int isp0 = isp1;
            isp1 = spost[ist][i];
            const double dx[3] = {
              0.5*(xsp[isp0][0]+xsp[isp1][0]) - x0[0],
              0.5*(xsp[isp0][1]+xsp[isp1][1]) - x0[1],
              0.5*(xsp[isp0][2]+xsp[isp1][2]) - x0[2] };
            const double yp = DOT_PRODUCT(dx,e1);
            const double zp = DOT_PRODUCT(dx,e2);
            const double rp = sqrt(yp*yp + zp*zp);
            r_min = min(r_min,rp);
            r_max = max(r_max,rp);
          }
        }
      }
      r_avg /= wgt;
      cout << " > r_avg: " << r_avg << " plus/minus: " << r_min-r_avg << " " << r_max-r_avg << endl;
      if ((r_avg-r_min > 0.05*dn)||(r_max-r_avg > 0.05*dn)) {
        cout << "Warning: cylinder may not be discretized accurately enough: dr > 5% of dn: " << max((r_avg-r_min)/dn,(r_max-r_avg)/dn) << endl;
      }

      //cout << "HACK: setting r_avg to 0.5" << endl;
      //r_avg = 0.5;

      // now check the length by projecting both links in e0...
      double xp_min = 0.0;
      double xp_max = 0.0;
      for (int ilink = 0; ilink < spoli_i[0]; ++ilink) {
        const int isp = spoli_v[ilink];
        double dx[3] = DIFF(xsp[isp],x0);
        double xp = DOT_PRODUCT(dx,e0);
        xp_min = min(xp_min,xp);
        xp_max = max(xp_max,xp);
      }
      cout << " > cylinder starts at 0 plus/minus: " << xp_min << " " << xp_max << endl;
      if ((-xp_min > 0.05*dn)||(xp_max > 0.05*dn)) {
        cout << "Error: cylinder start not orthogonal: " << xp_min << " " << xp_max << endl;
        assert(0);
      }

      // second loop...
      xp_min = L;
      xp_max = L;
      for (int ilink = spoli_i[0]; ilink < spoli_i[1]; ++ilink) {
        const int isp = spoli_v[ilink];
        double dx[3] = DIFF(xsp[isp],x0);
        double xp = DOT_PRODUCT(dx,e0);
        xp_min = min(xp_min,xp);
        xp_max = max(xp_max,xp);
      }
      cout << " > cylinder ends at L plus/minus: " << xp_min-L << " " << xp_max-L << endl;
      if ((L-xp_min > 0.05*dn)||(xp_max-L > 0.05*dn)) {
        cout << "Error: cylinder end not orthogonal: " << xp_min << " " << xp_max << endl;
        assert(0);
      }

      // spacing function depends on dn,dt,n...
      assert(dt>dn);
      const double sf = pow(dt/dn,1.0/double(n-1));
      const double delta = (dt*sf-dn)/(sf-1.0);
      cout << " > BL sf: " << sf << " delta for " << n << " layers: " << delta << endl;

      assert(L>2.0*delta);
      const int nL = int((L-2.0*delta)/dt); // + 2*n
      const int nTheta = int(2.0*M_PI*r_avg/dt);

      cout << " > total point estimate: nL+2*n: " << nL+2*n << " x nTheta: " << nTheta << " x n: " << n << " = " << (nL+2*n)*nTheta*n << endl;

      // ===================================================
      // surface points: xsp_ff...
      // ===================================================

      const int nsp_ff = spoli_i[1] + nTheta*(4*n-2);
      double (*xsp_ff)[3] = new double[nsp_ff][3];

      // HACK: check everything gets set...
      for (int isp = 0; isp < nsp_ff; ++isp)
        FOR_I3 xsp_ff[isp][i] = HUGE_VAL;

      // copy the links into xsp_ff first...
      for (ilink = 0; ilink < spoli_i[1]; ++ilink) {
        const int isp = spoli_v[ilink];
        FOR_I3 xsp_ff[ilink][i] = xsp[isp][i];
      }

      // then the layers start at spoli_i[1]...
      for (int k = 0; k <= 2*(n-1); ++k) {
        for (int j = 0; j < nTheta; ++j) {
          double xp;
          double rp;
          if (k < n) {
            xp = dn*(pow(sf,k)-1.0)/(sf-1.0);
            rp = r_avg - dn*(pow(sf,k+1)-1.0)/(sf-1.0);
          }
          else {
            const int kp = 2*(n-1)-k;
            xp = L - dn*(pow(sf,kp+1)-1.0)/(sf-1.0);
            rp = r_avg - dn*(pow(sf,kp+1)-1.0)/(sf-1.0);
          }
          //if (j == 0) cout << "p0 k: " << k << " xp: " << xp << " L-xp: " << L-xp << endl;
          const double yp = rp*cos(2.0*M_PI*double(j)/double(nTheta));
          const double zp = rp*sin(2.0*M_PI*double(j)/double(nTheta));
          int isp_ff = spoli_i[1]+nTheta*2*k+j;
          assert(isp_ff < nsp_ff);
          FOR_I3 xsp_ff[isp_ff][i] = xp*e0[i] + yp*e1[i] + zp*e2[i];
          if (k < n-1) {
            xp = dn*(pow(sf,k+1)-1.0)/(sf-1.0);
          }
          else {
            const int kp = 2*(n-1)-k;
            xp = L-dn*(pow(sf,kp)-1.0)/(sf-1.0);
          }
          //if (j == 0) cout << "p1 k: " << k << " xp: " << xp << " L-xp: " << L-xp << endl;
          isp_ff = spoli_i[1]+nTheta*(2*k+1)+j;
          assert(isp_ff < nsp_ff);
          FOR_I3 xsp_ff[isp_ff][i] = xp*e0[i] + yp*e1[i] + zp*e2[i];
        }
      }

      // HACK: check everything gets set...
      for (int isp = 0; isp < nsp_ff; ++isp)
        FOR_I3 assert(xsp_ff[isp][i] != HUGE_VAL);

      // ===================================================
      // surface tris: xst_ff...
      // ===================================================

      const int nst_ff = spoli_i[1] + nTheta*4*(2*n-1);
      int (*spost_ff)[3] = new int[nst_ff][3];
      double *dxost_ff = new double[nst_ff];

      // 1. the lip tris are special because they have to respect both discretizations...
      int ist_ff = 0;
      {
        int isp_ff_prev = spoli_i[0]-2;
        int isp_ff = spoli_i[0]-1;
        // compute the theta of the middle...
        const double dx[3] = {
          0.5*(xsp_ff[isp_ff_prev][0]+xsp_ff[isp_ff][0]) - x0[0],
          0.5*(xsp_ff[isp_ff_prev][1]+xsp_ff[isp_ff][1]) - x0[1],
          0.5*(xsp_ff[isp_ff_prev][2]+xsp_ff[isp_ff][2]) - x0[2] };
        double yp = DOT_PRODUCT(dx,e1);
        double zp = DOT_PRODUCT(dx,e2);
        double theta = atan2(zp,yp);
        int j_prev = (int)floor(theta/(2.0*M_PI)*double(nTheta)+0.5);
        if (j_prev < 0) j_prev += nTheta;
        assert((j_prev >= 0)&&(j_prev < nTheta));
        isp_ff_prev = isp_ff;
        for (isp_ff = 0; isp_ff < spoli_i[0]; ++isp_ff) {
          // compute the theta...
          const double dx[3] = {
            0.5*(xsp_ff[isp_ff_prev][0]+xsp_ff[isp_ff][0]) - x0[0],
            0.5*(xsp_ff[isp_ff_prev][1]+xsp_ff[isp_ff][1]) - x0[1],
            0.5*(xsp_ff[isp_ff_prev][2]+xsp_ff[isp_ff][2]) - x0[2] };
          double yp = DOT_PRODUCT(dx,e1);
          double zp = DOT_PRODUCT(dx,e2);
          double theta = atan2(zp,yp);
          int j = (int)floor(theta/(2.0*M_PI)*double(nTheta)+0.5);
          if (j < 0) j += nTheta;
          assert((j >= 0)&&(j < nTheta));
          // tris from the nTheta side...
          while (j_prev != j) {
            int j_next = j_prev + 1;
            if (j_next == nTheta) j_next = 0;
            spost_ff[ist_ff][0] = spoli_i[1]+j_next;
            spost_ff[ist_ff][1] = isp_ff_prev;
            spost_ff[ist_ff][2] = spoli_i[1]+j_prev;
            dxost_ff[ist_ff] = dn;
            ++ist_ff;
            j_prev = j_next;
          }
          // tris from the lipline side...
          spost_ff[ist_ff][0] = isp_ff_prev;
          spost_ff[ist_ff][1] = spoli_i[1]+j;
          spost_ff[ist_ff][2] = isp_ff;
          dxost_ff[ist_ff] = dn;
          ++ist_ff;
          j_prev = j;
          isp_ff_prev = isp_ff;
        }
        assert(ist_ff == spoli_i[0]+nTheta);
      }

      // now the step tris...
      for (int k = 0; k <= 2*(n-1); ++k) {
        double dxr;
        double dxx;
        if (k < (n-1)) {
          dxr = dn*pow(sf,k);
          dxx = dn*pow(sf,k+1);
        }
        else {
          const int kp = 2*(n-1)-k;
          dxr = dn*pow(sf,kp);
          dxx = dn*pow(sf,kp);
        }

        //cout << "k: " << k << " dxr: " << dxr << " dxx: " << dxx << endl;

        int j_prev = nTheta-1;
        for (int j = 0; j < nTheta; ++j) {
          // 2 tris at the same r...
          assert(ist_ff < nst_ff);
          //assert(ist_ff == spoli_i[1]+2*k*nTheta+j);
          spost_ff[ist_ff][0] = spoli_i[1]+nTheta*2*k+j_prev;
          spost_ff[ist_ff][1] = spoli_i[1]+nTheta*(2*k+1)+j;
          spost_ff[ist_ff][2] = spoli_i[1]+nTheta*2*k+j;
          dxost_ff[ist_ff] = dxr;
          ++ist_ff;
          assert(ist_ff < nst_ff);
          spost_ff[ist_ff][0] = spoli_i[1]+nTheta*2*k+j_prev;
          spost_ff[ist_ff][1] = spoli_i[1]+nTheta*(2*k+1)+j_prev;
          spost_ff[ist_ff][2] = spoli_i[1]+nTheta*(2*k+1)+j;
          dxost_ff[ist_ff] = dxr;
          ++ist_ff;
          // 2 tris at the same x...
          if (k < 2*(n-1)) {
            assert(ist_ff < nst_ff);
            spost_ff[ist_ff][0] = spoli_i[1]+nTheta*(2*k+1)+j_prev;
            spost_ff[ist_ff][1] = spoli_i[1]+nTheta*(2*k+2)+j;
            spost_ff[ist_ff][2] = spoli_i[1]+nTheta*(2*k+1)+j;
            dxost_ff[ist_ff] = dxx;
            ++ist_ff;
            assert(ist_ff < nst_ff);
            spost_ff[ist_ff][0] = spoli_i[1]+nTheta*(2*k+1)+j_prev;
            spost_ff[ist_ff][1] = spoli_i[1]+nTheta*(2*k+2)+j_prev;
            spost_ff[ist_ff][2] = spoli_i[1]+nTheta*(2*k+2)+j;
            dxost_ff[ist_ff] = dxx;
            ++ist_ff;
          }
          // copy j into j_prev for next one...
          j_prev = j;
        }
      }

      // and the second ring of lip tris...
      {
        int isp_ff_prev = spoli_i[1]-2;
        int isp_ff = spoli_i[1]-1;
        // compute the theta of the middle...
        const double dx[3] = {
          0.5*(xsp_ff[isp_ff_prev][0]+xsp_ff[isp_ff][0]) - x0[0],
          0.5*(xsp_ff[isp_ff_prev][1]+xsp_ff[isp_ff][1]) - x0[1],
          0.5*(xsp_ff[isp_ff_prev][2]+xsp_ff[isp_ff][2]) - x0[2] };
        double yp = DOT_PRODUCT(dx,e1);
        double zp = DOT_PRODUCT(dx,e2);
        double theta = atan2(zp,yp);
        int j_prev = (int)floor(theta/(2.0*M_PI)*double(nTheta)+0.5);
        if (j_prev < 0) j_prev += nTheta;
        assert((j_prev >= 0)&&(j_prev < nTheta));
        isp_ff_prev = isp_ff;
        for (isp_ff = spoli_i[0]; isp_ff < spoli_i[1]; ++isp_ff) {
          // compute the theta...
          const double dx[3] = {
            0.5*(xsp_ff[isp_ff_prev][0]+xsp_ff[isp_ff][0]) - x0[0],
            0.5*(xsp_ff[isp_ff_prev][1]+xsp_ff[isp_ff][1]) - x0[1],
            0.5*(xsp_ff[isp_ff_prev][2]+xsp_ff[isp_ff][2]) - x0[2] };
          double yp = DOT_PRODUCT(dx,e1);
          double zp = DOT_PRODUCT(dx,e2);
          double theta = atan2(zp,yp);
          int j = (int)floor(theta/(2.0*M_PI)*double(nTheta)+0.5);
          if (j < 0) j += nTheta;
          assert((j >= 0)&&(j < nTheta));
          // tris from the nTheta side...
          while (j_prev != j) {
            int j_next = j_prev - 1;
            if (j_next == -1) j_next = nTheta-1;
            assert(ist_ff < nst_ff);
            spost_ff[ist_ff][0] = spoli_i[1]+(4*n-3)*nTheta+j_next;
            spost_ff[ist_ff][1] = isp_ff_prev;
            spost_ff[ist_ff][2] = spoli_i[1]+(4*n-3)*nTheta+j_prev;
            dxost_ff[ist_ff] = dn;
            ++ist_ff;
            j_prev = j_next;
          }
          // tris from the lipline side...
          assert(ist_ff < nst_ff);
          spost_ff[ist_ff][0] = isp_ff_prev;
          spost_ff[ist_ff][1] = spoli_i[1]+(4*n-3)*nTheta+j;
          spost_ff[ist_ff][2] = isp_ff;
          dxost_ff[ist_ff] = dn;
          ++ist_ff;
          j_prev = j;
          isp_ff_prev = isp_ff;
        }
      }
      assert(ist_ff = nst_ff);

      writeTecplotFF("part_ff.dat",xsp_ff,nsp_ff,spost_ff,nst_ff);
      writeSbinFF("part_ff.sbin",xsp_ff,nsp_ff,spost_ff,nst_ff);

      // ==================================================================
      // now the points...
      // ==================================================================

      FILE * fp = fopen("part_pts.dat","w");
      vector<double> ptsVec;

      for (int i = 0; i < n; ++i) {
        const double xp = dn*(pow(sf,0.5+double(i))-1.0)/(sf-1.0);
        // the number of points in the azimuthal direction depends on the spacing...
        const double dx = dn*pow(sf,i);
        const double this_nTheta = int(2.0*M_PI*r_avg/dx);
        if (i == n-1) assert(this_nTheta == nTheta);
        cout << "i: " << i << " dx: " << dx << " this_nTheta: " << this_nTheta << " nTheta: " << nTheta << endl;
        for (int k = 0; k < this_nTheta; ++k) {
          double theta = (double(k)+0.5)/double(this_nTheta)*2.0*M_PI;
          for (int j = 0; j <= i; ++j) {
            double rp = r_avg - dn*(pow(sf,0.5+double(j))-1.0)/(sf-1.0);
            double yp = rp*cos(theta);
            double zp = rp*sin(theta);
            double x[3];
            FOR_I3 x[i] = xp*e0[i] + yp*e1[i] + zp*e2[i];
            ptsVec.push_back(x[0]);
            ptsVec.push_back(x[1]);
            ptsVec.push_back(x[2]);
            ptsVec.push_back(dx*1.75);
            fprintf(fp,"%18.15le %18.15le %18.15le %18.15le\n",x[0],x[1],x[2],dx*1.75); // i.e. a little more than sqrt(3)
          }
        }
      }

      for (int i = 0; i < nL; ++i) {
        double xp = delta + (double(i)+0.5)/double(nL)*(L-2.0*delta);
        for (int k = 0; k < nTheta; ++k) {
          double theta = (double(k)+0.5)/double(nTheta)*2.0*M_PI;
          for (int j = 0; j < n; ++j) {
            double rp = r_avg - dn*(pow(sf,0.5+double(j))-1.0)/(sf-1.0);
            double yp = rp*cos(theta);
            double zp = rp*sin(theta);
            double x[3];
            FOR_I3 x[i] = xp*e0[i] + yp*e1[i] + zp*e2[i];
            ptsVec.push_back(x[0]);
            ptsVec.push_back(x[1]);
            ptsVec.push_back(x[2]);
            ptsVec.push_back((L-2.0*delta)/nL*1.75);
            fprintf(fp,"%18.15le %18.15le %18.15le %18.15le\n",x[0],x[1],x[2],(L-2.0*delta)/nL*1.75);
          }
        }
      }

      for (int i = n-1; i >= 0; --i) {
        double xp = L - dn*(pow(sf,0.5+double(i))-1.0)/(sf-1.0);
        const double dx = dn*pow(sf,i);
        const double this_nTheta = int(2.0*M_PI*r_avg/dx);
        if (i == n-1) assert(this_nTheta == nTheta);
        for (int k = 0; k < this_nTheta; ++k) {
          double theta = (double(k)+0.5)/double(this_nTheta)*2.0*M_PI;
          for (int j = 0; j <= i; ++j) {
            double rp = r_avg - dn*(pow(sf,0.5+double(j))-1.0)/(sf-1.0);
            double yp = rp*cos(theta);
            double zp = rp*sin(theta);
            double x[3];
            FOR_I3 x[i] = xp*e0[i] + yp*e1[i] + zp*e2[i];
            ptsVec.push_back(x[0]);
            ptsVec.push_back(x[1]);
            ptsVec.push_back(x[2]);
            ptsVec.push_back(dx*1.75);
            fprintf(fp,"%18.15le %18.15le %18.15le %18.15le\n",x[0],x[1],x[2],dx*1.75);
          }
        }
      }

      fclose(fp);

      // now write a "part" file...

      {

        string name = "bl.part";
        cout << "writing part file: " << name << "..." << endl;

        FILE * fp = fopen(name.c_str(),"wb");
        assert(fp != NULL);

        int ibuf[5] = { PART_IO_MAGIC_NUMBER, PART_IO_VERSION, 0, 1, 1 }; // surface, ff_surface, points
        fwrite(ibuf,sizeof(int),5,fp);

        // ff surface...

        ibuf[0] = 1;
        ibuf[1] = nsp_ff;
        ibuf[2] = nst_ff;
        fwrite(ibuf,sizeof(int),3,fp);

        string ff_name = "FF";
        const int length = ff_name.length();
        fwrite(&length,sizeof(int),1,fp);
        fwrite(ff_name.c_str(),sizeof(char),length,fp);
        cout << " > ff_surface zone " << ff_name << endl;

        // surface points...

        fwrite(xsp_ff,sizeof(double),nsp_ff*3,fp);

        // surface tris...

        fwrite(spost_ff,sizeof(int),nst_ff*3,fp);

        int * znost_ff = new int[nst_ff];
        for (int ist = 0; ist < nst_ff; ++ist)
          znost_ff[ist] = 0;
        fwrite(znost_ff,sizeof(int),nst_ff,fp);
        delete[] znost_ff;

        // tri length scale...

        fwrite(dxost_ff,sizeof(double),nst_ff,fp);

        // ----------------------------
        // finally the points...
        // ----------------------------

        const int np = ptsVec.size()/4;
        cout << " > writing np points: " << np << endl;

        assert(ptsVec.size()%4 == 0);
        fwrite(&np,sizeof(int),1,fp);

        for (int ip = 0; ip < np; ++ip) {
          fwrite(&ptsVec[ip*4],sizeof(double),3,fp);
        }

        for (int ip = 0; ip < np; ++ip) {
          fwrite(&ptsVec[ip*4+3],sizeof(double),1,fp);
        }

        fclose(fp);

      }

      cout << "Looking good!" << endl;
      //getchar();

    }

    // report...
    cout << " > group: " << igr << " normal: " << COUT_VEC(n_sz[igr]) << " centroid: " << COUT_VEC(x_sz[igr]) << " area: " << area_sz[igr] << " nmag/area: " << n_mag/area_sz[igr] << endl;

  }

  delete[] n_sz;
  delete[] x_sz;
  delete[] area_sz;

}

void SimpleSurface::writeTecplotFF(const string& name,const double (* const xsp_ff)[3],const int nsp_ff,const int (* const spost_ff)[3],const int nst_ff) {

  cout << "SimpleSurface::writeTecplotFF: " << name << " nsp_ff: " << nsp_ff << " nst_ff: " << nst_ff << endl;

  FILE * fp = fopen(name.c_str(),"w");
  assert(fp != NULL);

  fprintf(fp,"TITLE = \"%s\"\n",name.c_str());
  fprintf(fp,"VARIABLES = \"X\"\n");
  fprintf(fp,"\"Y\"\n");
  fprintf(fp,"\"Z\"\n");

  fprintf(fp,"ZONE T=\"FF\"\n");
  fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nsp_ff,nst_ff);

  for (int isp = 0; isp < nsp_ff; ++isp) {
    fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_ff[isp][0],xsp_ff[isp][1],xsp_ff[isp][2]);
  }

  for (int ist = 0; ist < nst_ff; ++ist) {
    fprintf(fp,"%d %d %d\n",spost_ff[ist][0]+1,spost_ff[ist][1]+1,spost_ff[ist][2]+1);  //tecplot file is 1-indexed
  }

  fclose(fp);

}

void SimpleSurface::writeSbinFF(const string& name,const double (* const xsp_ff)[3],const int nsp_ff,const int (* const spost_ff)[3],const int nst_ff) {

  cout << "SimpleSurface::writeSbinFF: " << name << " nsp_ff: " << nsp_ff << " nst_ff: " << nst_ff << endl;

  FILE * fp = fopen(name.c_str(),"wb");
  assert(fp != NULL);

  const int version = 1;
  fwrite(&version,sizeof(int),1,fp);

  vector<string> zoneVec;
  zoneVec.push_back("FF");

  const int count = zoneVec.size();
  fwrite(&count,sizeof(int),1,fp);
  cout << " > writing " << count << " zones" << endl;
  for (int izone = 0, limit = zoneVec.size(); izone < limit; ++izone) {
    const int length = zoneVec[izone].length();
    fwrite(&length,sizeof(int),1,fp);
    fwrite(zoneVec[izone].c_str(),sizeof(char),length,fp);
    cout << " > zone " << zoneVec[izone] << endl;
  }

  fwrite(&nsp_ff,sizeof(int),1,fp);
  fwrite(xsp_ff,sizeof(double),nsp_ff*3,fp);

  fwrite(&nst_ff,sizeof(int),1,fp);
  fwrite(spost_ff,sizeof(int),nst_ff*3,fp);

  int * znost_ff = new int[nst_ff];
  for (int ist = 0; ist < nst_ff; ++ist)
    znost_ff[ist] = 0;
  fwrite(znost_ff,sizeof(int),nst_ff,fp);
  delete[] znost_ff;

  fclose(fp);

}


void SimpleSurface::queryFlaggedSubzones() {

  cout << "SimpleSurface::queryFlaggedSubzones()" << endl;

  vector<pair<int,int> > flagVec;
  map<const pair<int,int>,pair<int,int> > edgeMap;
  for (int ist = 0; ist < nst; ++ist) {
    const int isz = szost[ist];
    if (sz_flag[isz]) {
      int count = flagVec.size();
      flagVec.push_back(pair<int,int>(ist,count));
      FOR_I3 {
        const int isp0 = spost[ist][i];
        const int isp1 = spost[ist][(i+1)%3];
        map<const pair<int,int>,pair<int,int> >::iterator iter = edgeMap.find(pair<int,int>(isp1,isp0));
        if (iter != edgeMap.end()) {
          int nbr_count = iter->second.first;
          edgeMap.erase(iter);
          while (nbr_count != flagVec[nbr_count].second)
            nbr_count = flagVec[nbr_count].second;
          if (nbr_count < count) {
            assert(flagVec[count].second == count);
            flagVec[count].second = nbr_count;
            count = nbr_count;
          }
          else if (count < nbr_count) {
            assert(flagVec[nbr_count].second == nbr_count);
            flagVec[nbr_count].second = count;
            nbr_count = count;
          }
        }
        else {
          edgeMap[pair<int,int>(isp0,isp1)] = pair<int,int>(count,ist);
        }
      }
    }
  }

  cout << " > total tris involved: " << flagVec.size() << " out of " << nst << endl;

  int ngr = 0;
  for (int ii = 0, limit = flagVec.size(); ii < limit; ++ii) {
    int count = flagVec[ii].second;
    while ((count >= 0)&&(count != flagVec[count].second))
      count = flagVec[count].second;
    if (count >= 0) {
      ++ngr;
      flagVec[ii].second = -ngr; // -1,-2,...
    }
    else {
      flagVec[ii].second = count;
    }
  }

  cout << " > disjoint group count: " << ngr << endl;

  double (*n_sz)[3] = new double[ngr][3];
  double (*x_sz)[3] = new double[ngr][3];
  double *area_sz   = new double[ngr];

  for (int igr = 0; igr < ngr; ++igr) {
    FOR_I3 n_sz[igr][i] = 0.0;
    FOR_I3 x_sz[igr][i] = 0.0;
    area_sz[igr] = 0.0;
  }

  for (int ii = 0, limit = flagVec.size(); ii < limit; ++ii) {
    const int ist = flagVec[ii].first; assert((ist >= 0)&&(ist < nst));
    const int igr = -flagVec[ii].second-1; assert((igr >= 0)&&(igr < ngr));
    const double n2[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    FOR_I3 n_sz[igr][i] += n2[i];
    const double area = MAG(n2);
    area_sz[igr] += area;
    FOR_I3 x_sz[igr][i] += area*(xsp[spost[ist][0]][i] + xsp[spost[ist][1]][i] + xsp[spost[ist][2]][i]); // don't forget about the factor of 3
  }

  for (int igr = 0; igr < ngr; ++igr) {

    cout << " > ======================= group " << igr << " ==========================" << endl;

    FOR_I3 n_sz[igr][i] *= 0.5;
    FOR_I3 x_sz[igr][i] /= area_sz[igr]*3.0;
    area_sz[igr] *= 0.5;
    const double n_mag = MAG(n_sz[igr]);
    FOR_I3 n_sz[igr][i] /= n_mag;
    cout << " > area-weighted x: " << COUT_VEC(x_sz[igr]) << " unit n: " << COUT_VEC(n_sz[igr]) << " area: " << area_sz[igr] << " n_mag/area: " << n_mag/area_sz[igr] << endl;
    // flag everybody...
    st_flag.resize(nst);
    st_flag.setAll(0);
    for (int ii = 0, limit = flagVec.size(); ii < limit; ++ii) {
      const int ist = flagVec[ii].first; assert((ist >= 0)&&(ist < nst));
      const int this_igr = -flagVec[ii].second-1; assert((igr >= 0)&&(igr < ngr));
      if (this_igr == igr)
        st_flag[ist] = 1;
    }
    // now find the edge loops associated with the ends...
    map<const int,int> linkMap;
    for (map<const pair<int,int>,pair<int,int> >::iterator iter = edgeMap.begin(); iter != edgeMap.end(); ++iter) {
      const int ist = iter->second.second;
      if (st_flag[ist] == 1) {
        const int isp0 = iter->first.first;
        const int isp1 = iter->first.second;
        linkMap.insert(pair<int,int>(isp0,isp1));
      }
    }
    // how many loops...
    int nlink = 0;
    vector<double> doubleVec;
    while (!linkMap.empty()) {
      double xc[3] = { 0.0, 0.0, 0.0 };
      double wgt = 0.0;
      int count = 0;
      map<const int,int>::iterator iter = linkMap.begin();
      assert(iter != linkMap.end());
      const int start = iter->first;
      int next = start;
      do {
        iter = linkMap.find(next);
        assert(iter != linkMap.end());
        const int prev = next;
        next = iter->second;
        linkMap.erase(iter);
        // add segment prev->next to xc,wgt...
        const double this_wgt = DIST(xsp[next],xsp[prev]);
        FOR_I3 xc[i] += this_wgt*(xsp[prev][i]+xsp[next][i]);
        wgt += this_wgt; // don't forget factor of 2 when normalizing!
        ++count;
      } while (next != start);
      FOR_I3 {
        xc[i] /= 2.0*wgt;
        doubleVec.push_back(xc[i]);
      }
      ++nlink;
      cout << " > found open edge loop " << nlink << " of " << count << " edges, length: " << wgt << " with center: " << COUT_VEC(xc) << endl;
    }
    // --------------------------------------------------------------------------------------
    // depending on the number of links and other geometric checks, explain wtf this is...
    // --------------------------------------------------------------------------------------
    if ((nlink == 1)&&(n_mag/area_sz[igr] > 0.8)) {
      cout << " > seems pretty planar, and one open edge loop: might be a disk..." << endl;
      const double xc[3] = { doubleVec[0], doubleVec[1], doubleVec[2] };
      double r_min = HUGE_VAL;
      double r_max = 0.0;
      // go through edges again and record min/max radius...
      for (map<const pair<int,int>,pair<int,int> >::iterator iter = edgeMap.begin(); iter != edgeMap.end(); ++iter) {
        const int ist = iter->second.second;
        if (st_flag[ist] == 1) {
          const int isp0 = iter->first.first;
          const int isp1 = iter->first.second;
          const double r0 = DIST(xsp[isp0],xc);
          const double r1 = DIST(xsp[isp0],xc);
          r_min = min(r_min,min(r0,r1));
          r_max = max(r_max,max(r0,r1));
          double x_mid[3]; FOR_I3 x_mid[i] = 0.5*(xsp[isp0][i]+xsp[isp1][i]);
          r_min = min(r_min,DIST(x_mid,xc));
        }
      }
      cout << " > disk rmin: " << r_min << " rmax: " << r_max << " rmax-rmin: " << r_max-r_min << endl;
    }
    else if ((nlink == 2)&&(n_mag/area_sz[igr] < 0.8)) {
      cout << " > looking more and more like a cylinder: 2 open edge loops!" << endl;
      const double x0[3] = { doubleVec[0], doubleVec[1], doubleVec[2] };
      const double x1[3] = { doubleVec[3], doubleVec[4], doubleVec[5] };
      double e0[3] = DIFF(x1,x0);
      double mag_e0 = MAG(e0);
      FOR_I3 e0[i] /= mag_e0;
      double e1[3],e2[3];
      MiscUtils::getBestE1E2FromE0(e1,e2,e0);
      const double L =
        (x1[0]-x0[0])*e0[0] +
        (x1[1]-x0[1])*e0[1] +
        (x1[2]-x0[2])*e0[2];
      double x_mid[3]; FOR_I3 x_mid[i] = 0.5*(x0[i]+x1[i]);
      double dx[3] = DIFF(x1,x0);
      const double mag_dx = MAG(dx);
      FOR_I3 dx[i] /= mag_dx;
      cout << " > cylinder length: " << L << " center: " << COUT_VEC(x_mid) << " axis: " << COUT_VEC(dx) << endl;
      // and figure out the radius from the tris...
      double r_min = HUGE_VAL;
      double r_max = 0.0;
      for (int ist = 0; ist < nst; ++ist) {
        if (st_flag[ist] == 1) {
          FOR_I3 {
            const int isp = spost[ist][i];
            const double dx[3] = DIFF(xsp[isp],x0);
            const double yp = DOT_PRODUCT(dx,e1);
            const double zp = DOT_PRODUCT(dx,e2);
            const double rp = sqrt(yp*yp + zp*zp);
            r_min = min(r_min,rp);
            r_max = max(r_max,rp);
          }
          // the minimum radius for reasonable tris will probably occur along an edge...
          int isp1 = spost[ist][2];
          FOR_I3 {
            const int isp0 = isp1;
            isp1 = spost[ist][i];
            const double dx[3] = {
              0.5*(xsp[isp0][0]+xsp[isp1][0]) - x0[0],
              0.5*(xsp[isp0][1]+xsp[isp1][1]) - x0[1],
              0.5*(xsp[isp0][2]+xsp[isp1][2]) - x0[2] };
            const double yp = DOT_PRODUCT(dx,e1);
            const double zp = DOT_PRODUCT(dx,e2);
            const double rp = sqrt(yp*yp + zp*zp);
            r_min = min(r_min,rp);
            r_max = max(r_max,rp);
          }
        }
      }
      cout << " > cylinder rmin: " << r_min << " rmax: " << r_max << " rmax-rmin: " << r_max-r_min << endl;
    }
    else {
      cout << " > found " << nlink << " open edge loops. Not sure what this is (yet)." << endl;
    }
  }
  cout << " > ==========================================================" << endl;


}

void SimpleSurface::selectZonesInWindow(const double window[4][3],const bool b_strictly_inside) {

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

  bool* b_inside_sz = new bool[nsz];
  if (b_strictly_inside) {
    for (int isz = 0; isz < nsz; ++isz) b_inside_sz[isz] = true;
    for (int ist = 0; ist < nst; ++ist) {
      if (b_inside_sz[szost[ist]]) {
        // count points inside window...
        int n_inside = 0;
        FOR_I3 if (b_inside_sp[spost[ist][i]]) n_inside++;
        if (n_inside < 3)
          b_inside_sz[szost[ist]] = false;
      }
    }
  }
  else {
    for (int isz = 0; isz < nsz; ++isz) b_inside_sz[isz] = false;
    for (int ist = 0; ist < nst; ++ist) {
      if (!b_inside_sz[szost[ist]]) {
        FOR_I3 {
          if (b_inside_sp[spost[ist][i]]) {
            b_inside_sz[szost[ist]] = true;
            break;
          }
        }
      }
    }
  }
  delete[] b_inside_sp;

  // remove hidden zones...
  for (int ii = 0, lim = hiddenSubzoneVec.size(); ii < lim; ++ii) {
    const int isz = hiddenSubzoneVec[ii];
    if ((isz >= 0)&&(isz < nsz)) {
      b_inside_sz[isz] = false;
    }
    else {
      CWARN(" > hidden subzone: " << isz << " not in subzone range, [0," << nsz << "); skipping..." )
    }
  }

  for (int isz = 0; isz < nsz; ++isz) {
    if (b_inside_sz[isz]) {
      selectedSubzoneVec.push_back(isz);
      //cout << isz << endl;
    }
  }
  delete[] b_inside_sz;

}

void SimpleSurface::buildStostForFlaggedZones(int (**stost)[3]) {
  ensureTeost();

  // build stost structure where:
  // -2: indicates an invalid edge, should only exist for tris that are not participating, if found in valid tri then an error occurred
  // -1: an open or boundary edge, no nbr tri but this edge needs protecting
  // 0+: ist index of ist_nbr, who is part of the flagged region

  *stost = new int[nst][3];
  for (int ist=0; ist<nst; ist++) {
    FOR_I3 (*stost)[ist][i] = -2;  // default, should only remain for tris that aren't valid
  }

  for (int ist=0; ist<nst; ist++) {
    if (zone_flag[znost[ist]] == 1) {
      FOR_I3 {
        if (isEdgeOpen(ist,i)) {
          (*stost)[ist][i] = -1;
        }
        else if (isEdgeMulti(ist,i)) {
          vector<pair<int,int> > ist_nbrs;  // <ist,orient> for ime nbrs
          teost.getMultiNbrs(ist_nbrs,ist,i);

          bool has_flagged_nbr = false;
          for (vector<pair<int,int> >::iterator it=ist_nbrs.begin(); it!=ist_nbrs.end(); ++it) {
            int ist_nbr = it->first;
            if (ist_nbr != ist) {
              if (zone_flag[znost[ist_nbr]] == 1) {
                has_flagged_nbr = true;
                if ((*stost)[ist][i] != -2) {
                  WUI(WARN,"cannot currently retessellate across multi-edges, please select a manifold set of zones; exiting");
                  throw(-1);
                }
                (*stost)[ist][i] = ist_nbr;
              }
            }
          }
          if (!has_flagged_nbr) {
            (*stost)[ist][i] = -1;  // like an open edge, stopped at a multi edge
          }
        }
        else {
          // should have valid tri nbr
          if (!isNbrAligned(ist,i)) {
            WUI(WARN,"problem with tri alignement, please fix before attempting retessellation; skipping");
            throw(-1);
          }
          int ist_nbr;
          getAlignedTriNbr(ist_nbr,ist,i);
          if (zone_flag[znost[ist_nbr]] == 1) (*stost)[ist][i] = ist_nbr;
          else (*stost)[ist][i] = -1;
        }
      }
    }
  }

  // stost check
  for (int ist=0; ist<nst; ++ist) {
    if (zone_flag[znost[ist]] == 1) {
      //cout << "* tri " << ist << " " << COUT_VEC((*stost)[ist]) << endl;
      FOR_I3 assert((*stost)[ist][i] != -2);
    }
    else {
      //cout << "  tri " << ist << " " << COUT_VEC((*stost)[ist]) << endl;
      FOR_I3 assert((*stost)[ist][i] == -2);
    }
  }
}

double SimpleSurface::identifyFeatureXp(vector<NewNode>& fixedXp,vector<pair<int,int> >& featureEdgeVec,vector<pair<int,int> >& adj_tris,vector<pair<int,int> >& isp_to_nearest_fixed,const vector<int>& edge_indices,const double delta,const bool no_adj_feature,const bool b_keep_edge_nodes) {
  // given a collection of edge loop indices, retesselate
  //const double dp_tol = feature_cos;

  // need a node to feature edge search structure too...
  multimap<int,int> ispToEdgeMap;

  set<pair<int,int> > featureEdgeSet;  // to avoid unique entries use set, then convert

  IntFlag egrFlag(eoiGroupDataVec.size());
  egrFlag.setAll(0);
  for (vector<int>::const_iterator eit=edge_indices.begin(); eit!=edge_indices.end(); ++eit) {
    egrFlag[*eit] = 1;
  }

  // count feature edge touches for valence
  sp_flag.resize(nsp);
  sp_flag.setAll(0);

  // loop half edges and add to feature edges without duplication
  for (map<uint,int>::const_iterator eit=eoi_to_group.begin(); eit!=eoi_to_group.end(); ++eit) {
    if (egrFlag[eit->second]) {
      uint ist; uchar i;
      unpackEdge(ist,i,eit->first);
      const int isp0 = spost[ist][i];
      const int isp1 = spost[ist][(i+1)%3];
      pair<set<pair<int,int> >::iterator,bool> ret = featureEdgeSet.insert( pair<int,int>(min(isp0,isp1),max(isp0,isp1)) );
      if (ret.second == true) {
        // new element inserted, so add to vector
        sp_flag[isp0] += 2;  // increment by 2 b/c we are uniquely adding this edge
        sp_flag[isp1] += 2;
        const int ise = featureEdgeVec.size();
        featureEdgeVec.push_back(pair<int,int>(min(isp0,isp1),max(isp0,isp1)));
        ispToEdgeMap.insert(pair<int,int> (isp0,ise));
        ispToEdgeMap.insert(pair<int,int> (isp1,ise));
      }
    }
  }

  featureEdgeSet.clear();

  int f_nodes = sp_flag.countPositive();
  cout << " > number of unique feature edges,nodes: " << featureEdgeVec.size() << " " << f_nodes << endl;

  if (false) {
    // feature diagnostics
    FILE * fp = fopen("feature_edges.dat","w");
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");
    fprintf(fp,"\"valence\"\n");
    fprintf(fp,"ZONE T=\"feature edges\"\n");
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",int(featureEdgeVec.size()*3),int(featureEdgeVec.size()));
    for (vector<pair<int,int> >::iterator it=featureEdgeVec.begin(); it!=featureEdgeVec.end(); ++it) {
      fprintf(fp,"%lf %lf %lf %d\n",xsp[it->first][0],xsp[it->first][1],xsp[it->first][2],sp_flag[it->first]);
      fprintf(fp,"%lf %lf %lf %d\n",xsp[it->second][0],xsp[it->second][1],xsp[it->second][2],sp_flag[it->second]);
      fprintf(fp,"%lf %lf %lf %d\n",xsp[it->first][0],xsp[it->first][1],xsp[it->first][2],sp_flag[it->first]);
    }
    for (vector<pair<int,int> >::iterator it=featureEdgeVec.begin(); it!=featureEdgeVec.end(); ++it) {
      const int offset = 3*(it-featureEdgeVec.begin())+1;  // tecplot is 1-indexed
      fprintf(fp,"%d %d %d\n",offset,offset+1,offset+2);
    }
    fclose(fp);
  }

  zone_flag.resize(zoneVec.size());
  zone_flag.setAll(0); // adjacent tris can only live on unflagged surface

  return rediscretizeEdges(fixedXp,adj_tris,isp_to_nearest_fixed,ispToEdgeMap,featureEdgeVec,delta,no_adj_feature,true,b_keep_edge_nodes);
}

double SimpleSurface::identifyFeatureXp(vector<NewNode>& fixedXp,vector<pair<int,int> >& featureEdgeVec,vector<pair<int,int> >& adj_tris,vector<pair<int,int> >& isp_to_nearest_fixed,const int (*stost)[3],const double delta,const bool no_adj_feature,const bool b_keep_edge_nodes) {
  const double dp_tol = feature_cos;

  // need a node to feature edge search structure too...
  multimap<int,int> ispToEdgeMap;

  // count feature edge touches for valence
  sp_flag.resize(nsp);
  sp_flag.setAll(0);

  // use bits to indicate which half-edges are boundaries or features AND whether this tri has been visited or not:
  st_flag.resize(nst);
  st_flag.setAll(0);
  // & 1: has been visited
  // & 2: 0-1 edge is a feature
  // & 4: 1-2 edge is a feature
  // & 8: 2-0 edge is a feature
  for (int ist = 0; ist < nst; ++ist) {
    if (zone_flag[znost[ist]] == 1) {
      st_flag[ist] |= 1;
      double * n_ist = NULL;
      FOR_I3 {
        const int ist_nbr = stost[ist][i];
        if (ist_nbr < 0) {
          // open or boundary edge should be protected
          assert(ist_nbr == -1);

          const int isp0 = spost[ist][i];
          const int isp1 = spost[ist][(i+1)%3];
          // special treatment of boundary/open edges so gegt processed as chains
          // if we don't do this their nodes only have half-edge valence 2, which looks like a terminus
          sp_flag[isp0] += 2;
          sp_flag[isp1] += 2;
          // ++sp_flag[isp0];
          // ++sp_flag[isp1];

          // flag this tri's edge
          switch (i) {
            case 0:
              st_flag[ist] |= 2;
              break;
            case 1:
              st_flag[ist] |= 4;
              break;
            case 2:
              st_flag[ist] |= 8;
              break;
            default:
              assert(0);  // how did you get here?
          }

          // add edge to feature-edge search structures
          const int ise = featureEdgeVec.size();
          featureEdgeVec.push_back(pair<int,int>(min(isp0,isp1),max(isp0,isp1)));
          ispToEdgeMap.insert(pair<int,int> (isp0,ise));
          ispToEdgeMap.insert(pair<int,int> (isp1,ise));
        }
        else if (ist_nbr >= 0) {
          assert(zone_flag[znost[ist_nbr]] == 1);

          // check if this edge is a crease
          if (n_ist == NULL) {
            n_ist = new double[3];
            const double n_tmp[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
            FOR_I3 n_ist[i] = n_tmp[i];
          }

          const double n_ist_nbr[3] = TRI_NORMAL_2(xsp[spost[ist_nbr][0]],xsp[spost[ist_nbr][1]],xsp[spost[ist_nbr][2]]);

          const double n_ist_mag = sqrt(DOT_PRODUCT(n_ist,n_ist)); assert(n_ist_mag > 0.0);
          const double n_ist_nbr_mag = sqrt(DOT_PRODUCT(n_ist_nbr,n_ist_nbr)); assert(n_ist_nbr_mag > 0.0);

          // criteria for edge protection within flagged zones...
          const bool is_crease = DOT_PRODUCT(n_ist,n_ist_nbr)/(n_ist_mag*n_ist_nbr_mag) < dp_tol;
          const bool is_zone_boundary = (znost[ist] == znost[ist_nbr]) ? false:true;
          const bool is_subzone_boundary = (szost[ist] == szost[ist_nbr]) ? false:true;

          if ( is_crease || is_zone_boundary || is_subzone_boundary ) {
            // crease angle criterion hit
            const int isp0 = spost[ist][i];
            const int isp1 = spost[ist][(i+1)%3];
            ++sp_flag[isp0];
            ++sp_flag[isp1];

            // flag this tri's feature edge
            switch (i) {
              case 0:
                st_flag[ist] |= 2;
                break;
              case 1:
                st_flag[ist] |= 4;
                break;
              case 2:
                st_flag[ist] |= 8;
                break;
              default:
                assert(0);  // how did you get here?
            }


            // add edge to feature-edge search structures if hasn't been added already
            if (st_flag[ist_nbr] == 0) {
              const int ise = featureEdgeVec.size();
              featureEdgeVec.push_back(pair<int,int>(min(isp0,isp1),max(isp0,isp1)));
              ispToEdgeMap.insert(pair<int,int> (isp0,ise));
              ispToEdgeMap.insert(pair<int,int> (isp1,ise));
            }
          }
        }
      }
      DELETE(n_ist);
    }
  }

  int f_nodes = sp_flag.countPositive();
  cout << " > number of unique feature edges,nodes: " << featureEdgeVec.size() << " " << f_nodes << endl;

  if (false) {
    // feature diagnostics
    FILE * fp = fopen("feature_edges.dat","w");
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");
    fprintf(fp,"\"valence\"\n");
    fprintf(fp,"ZONE T=\"feature edges\"\n");
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",int(featureEdgeVec.size())*3,int(featureEdgeVec.size()));
    for (vector<pair<int,int> >::iterator it=featureEdgeVec.begin(); it!=featureEdgeVec.end(); ++it) {
      fprintf(fp,"%lf %lf %lf %d\n",xsp[it->first][0],xsp[it->first][1],xsp[it->first][2],sp_flag[it->first]);
      fprintf(fp,"%lf %lf %lf %d\n",xsp[it->second][0],xsp[it->second][1],xsp[it->second][2],sp_flag[it->second]);
      fprintf(fp,"%lf %lf %lf %d\n",xsp[it->first][0],xsp[it->first][1],xsp[it->first][2],sp_flag[it->first]);
    }
    for (vector<pair<int,int> >::iterator it=featureEdgeVec.begin(); it!=featureEdgeVec.end(); ++it) {
      const int offset = 3*(it-featureEdgeVec.begin())+1;  // tecplot is 1-indexed
      fprintf(fp,"%d %d %d\n",offset,offset+1,offset+2);
    }
    fclose(fp);
  }

  return rediscretizeEdges(fixedXp,adj_tris,isp_to_nearest_fixed,ispToEdgeMap,featureEdgeVec,delta,no_adj_feature,false,b_keep_edge_nodes);
}

double SimpleSurface::rediscretizeEdges(vector<NewNode>& fixedXp,vector<pair<int,int> >& adj_tris,vector<pair<int,int> >& isp_to_nearest_fixed,const multimap<int,int>& ispToEdgeMap,const vector<pair<int,int> >& featureEdgeVec,const double delta,const bool no_adj_feature,const bool edge_only,const bool b_keep_edge_nodes) {

  // Notes/assumptions:
  // 1. half-edge touches per node are stored in sp_flag
  // 2. zone_flag stores which zones we should look for adj_tris within


  map<int,int> ispToFixedXp;  // for collocated nodes only...
  double min_effective_delta = -1.0;

  // go through and register points of interest (endpoints and chain intersections)
  // loop through to assess nodes with non-zero valence:
  // 0: no valence, i.e., no protected edges touch this node
  // 1+: half edge valence, positive indicates this node hasn't been processed yet
  // 1-: single-edge valence, negative indicates edge has been visited

  // also keep track of adjacent tris on the unflagged region in adj_tris
  for (int ist = 0; ist < nst; ++ist) {
    if (edge_only || zone_flag[znost[ist]] == 1) {
      FOR_I3 {
        const int isp=spost[ist][i];
        if (sp_flag[isp] > 0) {
          assert(sp_flag[isp]%2 == 0);  // valence counts half edge touches on the node. If selection is manifold, then should always be even...?
          const int valence = sp_flag[isp]/2;
          sp_flag[isp] = -valence;  // -indexed valence; negative so we don't visit again in this loop

          if (valence == 2) {
            // not an endpoint, but sits between feature edges
          }
          else {
            if (valence == 2)
              --sp_flag[isp]; // treat as a terminus and not a chain node

            // store this point, is a terminus for a feature edge segment
            fixedXp.push_back(NewNode());
            FOR_I3 fixedXp.back().xsp[i] = xsp[isp][i];
            fixedXp.back().isp = isp;
            ispToFixedXp.insert(pair<int,int> (isp,fixedXp.size()-1));
          }
        }
      }
    }
    if (edge_only || zone_flag[znost[ist]] == 0) {
      int touches = 0;
      FOR_I3 {
        const int isp=spost[ist][i];
        if (sp_flag[isp] != 0) touches |= (1 << i);
      }
      // bits of touches indicate which edges are on adjacent
      if (touches) adj_tris.push_back(pair<int,int>(ist,touches));
    }
  }
  COUT2(" > counted " << fixedXp.size() << " terminal nodes based on valence");
  COUT2(" > counted " << adj_tris.size() << " adjacent tris");
  const int n_terminal = fixedXp.size();

  if (!no_adj_feature) {
    // there may be features from the adjacent tris side that require protecting that are not detected on the flagged zones
    // by looping these tris and checking the edges that touch nodes on the protected loop(s) can determine if requires protecting
    for (vector<pair<int,int> >::iterator it=adj_tris.begin(); it!=adj_tris.end(); ++it) {

      vector<int> check_edge;  // tells which edges to test for crease, and which tri node to update
      switch (it->second) {
        case 1:
        // node 0 touches, so check edges 0 & 2
          check_edge.push_back(0);
          check_edge.push_back(2);
          break;
        case 2:
          // node 1 touches, so check edges 1 & 0
          check_edge.push_back(1);
          check_edge.push_back(0);
          break;
        case 4:
          // node 2 touches, so check edges 2 & 1
          check_edge.push_back(2);
          check_edge.push_back(1);
          break;
        case 3:
          // edge 0 touches, so check edges 2 & 1
          check_edge.push_back(2);
          check_edge.push_back(1);
          break;
        case 5:
          // edge 2 touches, so check edges 0 & 1
          check_edge.push_back(0);
          check_edge.push_back(1);
          break;
        case 6:
          // edge 1 touches, so check edges 0 & 2
          check_edge.push_back(0);
          check_edge.push_back(2);
          break;
        case 7:
          // all nodes touch
          check_edge.push_back(0);
          check_edge.push_back(1);
          check_edge.push_back(2);
          break;
        default:
          CWARN("got a bad touches value of: " << it->second);
          assert(0);  // should not arrive here b/c bits must be set
      }

      if (!check_edge.empty() && (it->second!=7)) {  // case of all edges touching, do nothing....
        const int ist=it->first;
        for (vector<int>::const_iterator ceit=check_edge.begin(); ceit!=check_edge.end(); ++ceit) {
          if (isEdgeFeature(ist,*ceit)) {
            // determine the node to protect based on sp_flag
            int isp_to_protect;
            if (sp_flag[spost[ist][*ceit]] < 0) isp_to_protect = spost[ist][*ceit];
            else isp_to_protect = spost[ist][(*ceit+1)%3];

            if (sp_flag[isp_to_protect] == -2) {
              // only update if a chain node that hasn't been flagged as a terminus
              --sp_flag[isp_to_protect];  // treat as a terminus and not a chain node
              // store this point

              fixedXp.push_back(NewNode());
              FOR_I3 fixedXp.back().xsp[i] = xsp[isp_to_protect][i];
              fixedXp.back().isp = isp_to_protect;
              ispToFixedXp.insert(pair<int,int> (isp_to_protect,fixedXp.size()-1));
            }
          }
        }
      }

    }
  }

  const int n_endpoints = fixedXp.size();
  COUT2(" > adjacent feature termini: " << n_endpoints-n_terminal);
  COUT2(" > total high valence/feature nodes detected (endpoints):  " << n_endpoints);

  if (false) {
    // high valence diagnostics
    FILE * fp = fopen("valence_nodes.dat","w");
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");
    fprintf(fp,"\"valence\"\n");
    fprintf(fp,"ZONE T=\"feature nodes\"\n");
    if (n_endpoints) {
      for (vector<NewNode>::const_iterator it=fixedXp.begin(); it!=fixedXp.end(); ++it) {
        fprintf(fp,"%lf %lf %lf %d\n",it->xsp[0],it->xsp[1],it->xsp[2],sp_flag[it->isp]);
      }
    }
    else {
      // populate dummy data
      fprintf(fp,"0.0 0.0 0.0 0\n");
    }
    fclose(fp);
  }

  if (false) {
    // feature diagnostics
    FILE * fp = fopen("adjacent_tris.dat","w");
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");

    int ntris=0;
    for (vector<pair<int,int> >::iterator it=adj_tris.begin(); it!=adj_tris.end(); ++it) {
      if ((it->second == 1) || (it->second == 2) || (it->second == 4)) ++ntris;
    }

    fprintf(fp,"ZONE T=\"1-touch\"\n");
    if (ntris) {
      fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",ntris*3,ntris);
      for (vector<pair<int,int> >::iterator it=adj_tris.begin(); it!=adj_tris.end(); ++it) {
        if ((it->second == 1) || (it->second == 2) || (it->second == 4)) {
          FOR_I3 fprintf(fp,"%lf %lf %lf\n",xsp[spost[it->first][i]][0],xsp[spost[it->first][i]][1],xsp[spost[it->first][i]][2]);
        }
      }
      for (int itri=0; itri<ntris; ++itri) {
        const int offset = 3*itri + 1;  // tecplot is 1-indexed
        fprintf(fp,"%d %d %d\n",offset,offset+1,offset+2);
      }
    }
    else {
      // dummy value
      fprintf(fp,"N=3, E=1, F=FEPOINT, ET=TRIANGLE\n");
      FOR_I3 fprintf(fp,"0.0 0.0 0.0\n");
      fprintf(fp,"1 2 3\n");
    }

    ntris=0;
    for (vector<pair<int,int> >::iterator it=adj_tris.begin(); it!=adj_tris.end(); ++it) {
      if ((it->second == 3) || (it->second == 5) || (it->second == 6)) ++ntris;
    }

    fprintf(fp,"ZONE T=\"2-touch\"\n");
    if (ntris) {
      fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",ntris*3,ntris);
      for (vector<pair<int,int> >::iterator it=adj_tris.begin(); it!=adj_tris.end(); ++it) {
        if ((it->second == 3) || (it->second == 5) || (it->second == 6)) {
          FOR_I3 fprintf(fp,"%lf %lf %lf\n",xsp[spost[it->first][i]][0],xsp[spost[it->first][i]][1],xsp[spost[it->first][i]][2]);
        }
      }
      for (int itri=0; itri<ntris; ++itri) {
        const int offset = 3*itri + 1;  // tecplot is 1-indexed
        fprintf(fp,"%d %d %d\n",offset,offset+1,offset+2);
      }
    }
    else {
      // dummy value
      fprintf(fp,"N=3, E=1, F=FEPOINT, ET=TRIANGLE\n");
      FOR_I3 fprintf(fp,"0.0 0.0 0.0\n");
      fprintf(fp,"1 2 3\n");
    }

    ntris=0;
    for (vector<pair<int,int> >::iterator it=adj_tris.begin(); it!=adj_tris.end(); ++it) {
      if (it->second == 7) ++ntris;
    }

    fprintf(fp,"ZONE T=\"3-touch\"\n");
    if (ntris) {
      fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",ntris*3,ntris);
      for (vector<pair<int,int> >::iterator it=adj_tris.begin(); it!=adj_tris.end(); ++it) {
        if (it->second == 7) {
          FOR_I3 fprintf(fp,"%lf %lf %lf\n",xsp[spost[it->first][i]][0],xsp[spost[it->first][i]][1],xsp[spost[it->first][i]][2]);
        }
      }
      for (int itri=0; itri<ntris; ++itri) {
        const int offset = 3*itri + 1;  // tecplot is 1-indexed
        fprintf(fp,"%d %d %d\n",offset,offset+1,offset+2);
      }
    }
    else {
      // dummy value
      fprintf(fp,"N=3, E=1, F=FEPOINT, ET=TRIANGLE\n");
      FOR_I3 fprintf(fp,"0.0 0.0 0.0\n");
      fprintf(fp,"1 2 3\n");
    }

    fclose(fp);
  }

  // now we can group feature edges into chains by marching from the endpoints
  // if feature edges remain after all this marching, then loop-chains exist and
  // we must process these as well
  const int f_edges = featureEdgeVec.size();
  IntFlag se_flag(f_edges);
  se_flag.setAll(1);  // track which edges have already been visited

  // knowing edges are unique, we should be able to march from endpoint nodes
  // along edge-segments until all non-loops have been processed
  stack<pair<int,int> > isp_starts;  // old isp, newnode isp
  for (vector<NewNode>::const_iterator it=fixedXp.begin(); it!=fixedXp.end(); ++it) {
    isp_starts.push(pair<int,int> (it->isp,it-fixedXp.begin()));
  }

  // currently all nodes with any valence are negative
  // all midchain nodes (val == 2) should be positive,
  // we keep terminal nodes negative since we use the negative
  // criterion to stop chain marching
  for (int isp=0; isp<nsp; ++isp) {
    if (sp_flag[isp] == -2) sp_flag[isp] *= -1;
  }

  int nse_processed = 0;
  int isp0,isp1;
  vector<int> chain_edges;  // index is edge, +/- indicates ordering of nodes (isp_min,max or vice versa)
  int n_chains = 0;
  int isp_start_loop = 0;
  bool b_loop = false;
  while (nse_processed != f_edges) {

    // starting nodes are either identified endpoints or random node in a loop
    int isp_start;
    int new_node_start;
    if (!isp_starts.empty()) {
      // grab edges emanating from a terminal node
      isp_start = isp_starts.top().first;
      new_node_start = isp_starts.top().second;
      isp_starts.pop();
    }
    else {
      // all remaining edges are in loops, so take a random point and march
      b_loop = true;
      for ( ; isp_start_loop<nsp; ++isp_start_loop) {
        if (sp_flag[isp_start_loop] == 2) {
          --sp_flag[isp_start_loop];  // don't loop around on ourselves again, so 1 is a protected value
          break;
        }
      }
      isp_start = isp_start_loop;
      new_node_start = -1;
    }
    // COUT2(" > starting edge chain at (isp,flag): " << isp_start << " " << sp_flag[isp_start]);

    // for each touching edge, see if processed. if not, start marching
    pair <multimap<int,int>::const_iterator, multimap<int,int>::const_iterator> ret,ret2;
    ret = ispToEdgeMap.equal_range(isp_start);  // grab edges touching this node
    for (multimap<int,int>::const_iterator it=ret.first; it!=ret.second; ++it) {
      int ise = it->second;

      if (se_flag[ise]) {  // not visited yet

        isp0 = isp_start;

        // start of a chain
        double length = 0.0;
        chain_edges.clear();
        bool done = false;

        while (!done) {
          se_flag[ise] = 0;  // set edge as visited
          const int isp_min = featureEdgeVec[ise].first;
          const int isp_max = featureEdgeVec[ise].second;
          //bool reverse = false;
          if (isp0 == isp_min) {
            isp1 = isp_max;
            chain_edges.push_back(ise);
          }
          else {
            assert(isp0 == isp_max);
            isp1 = isp_min;
            //reverse = true;
            chain_edges.push_back(-ise-1);  // -1 indexed
          }
          // COUT2("   - link (edge,isp0,isp1): " << ise << " " << isp0 << " " << isp1);
          length += DIST(xsp[isp0],xsp[isp1]);
          isp0 = isp1;  // set new search head

          if ((sp_flag[isp1] < 0) || (sp_flag[isp1] == 1)) {
            // end of this chain
            done = true;
          }
          else {
            assert(sp_flag[isp1] == 2);
            --sp_flag[isp1];  // internal chain node, set as visited so don't loop around on ourselves

            // should have a valid next link in chain
            assert(ispToEdgeMap.count(isp1) == 2);
            ret2 = ispToEdgeMap.equal_range(isp1);
            for (multimap<int,int>::const_iterator it2=ret2.first; it2!=ret2.second; ++it2) {
              if (it2->second != ise) {
                ise = it2->second;  // update next edge
                break;
              }
            }
          }
        }
        // now all edges and total length should have been found
        ++n_chains;
        COUT2(" > feature edge chain[" << n_chains << "]:");
        COUT2("    > n_edges: " << chain_edges.size());
        COUT2("    > length: " << length);
        nse_processed += int(chain_edges.size());
        // COUT2(" > nse_processed/f_edges: " << nse_processed << "/" << f_edges);


        // the chain links are stored in chain_edges
        // we can now insert our fixed new points along this chain
        const double edge_delta = 1.0*delta;
        double effective_delta = length;
        if (edge_delta < length) {
          const int internal_pts = (int)floor(length/edge_delta);
          effective_delta = length/internal_pts;
        }
        COUT2("    > effective delta (uniform integer divisible): " << effective_delta);


        // store min edge discretization length to ensure interior points are seeded with appropriate lengthscale
        if (min_effective_delta < 0) min_effective_delta = effective_delta;
        else min_effective_delta = min(min_effective_delta,effective_delta);

        // want to keep track of nearest new node for all non-protected old-nodes
        // create a vector of isp (-1-indexed for old nodes, 0+ indexed for new nodes) and arc-length to do a sort and then determine nearest new nbr
        vector<pair<double,int> > len_to_isp;
        double len = 0.0;

        const double snap_tol = 5.0e-3;
        double to_next = effective_delta;
        for (vector<int>::iterator link=chain_edges.begin(); link!=chain_edges.end(); ++link) {
          int ise = *link;
          int isp0,isp1;
          bool reverse = false;
          if (ise < 0) {
            ise = -ise-1;
            isp1 = featureEdgeVec[ise].first;
            isp0 = featureEdgeVec[ise].second;
            reverse = true;
          }
          else {
            isp0 = featureEdgeVec[ise].first;
            isp1 = featureEdgeVec[ise].second;
          }

          // for loops we need to add the starting node
          if (b_loop && (link==chain_edges.begin())) {
            fixedXp.push_back(NewNode());
            FOR_I3 fixedXp.back().xsp[i] = xsp[isp0][i];
            fixedXp.back().isp = isp0;
            ispToFixedXp.insert(pair<int,int> (isp0,fixedXp.size()-1));
            assert(sp_flag[isp0] == 1);
            --sp_flag[isp0];  // indicate this point already has a node on it
          }

          if (link==chain_edges.begin()) {
            // whether loop or chain this node is protected new node
            if (new_node_start == -1) len_to_isp.push_back(pair<double,int> (len,int(fixedXp.size()-1)));
            else len_to_isp.push_back(pair<double,int> (len,new_node_start));
          }

          const double dx[3] = DIFF(xsp[isp1],xsp[isp0]);
          const double my_length = MAG(dx);

          if (b_keep_edge_nodes) {
            if (edge_delta < my_length) {
              const int internal_pts = (int)floor(my_length/edge_delta);
              effective_delta = my_length/internal_pts;
            }
            else {
              effective_delta = my_length;
            }
            to_next = 0.0;  // reset for each edge segment so original edge nodes are preserved
          }

          double my_traversed = 0.0;
          while (my_traversed<my_length) {
            if (to_next <= (my_length-my_traversed)) {
              my_traversed += to_next;
              len += to_next;
              to_next = effective_delta;

              const double frac = my_traversed/my_length;

              bool snap_to_isp1 = (((1.0-frac)*my_length) < (snap_tol*effective_delta)) ? true:false;
              bool snap_to_isp0 = (((frac)*my_length) < (snap_tol*effective_delta)) ? true:false;
              // if (snap_to_isp1&&snap_to_isp0) cout << "very small edge relative to lengthscale detected" << endl;  // if both true indicates nearly collocated nodes...

              if (!snap_to_isp1 && !snap_to_isp0) {
                // insert a new node, no snapping
                fixedXp.push_back(NewNode());
                FOR_I3 fixedXp.back().xsp[i] = xsp[isp0][i] + frac*dx[i];
                fixedXp.back().isp = -ise-1;  // negative index indicates node on an edge as opposed to on an existing node
                if (reverse) fixedXp.back().data_d = (1.0-frac);
                else fixedXp.back().data_d = frac;

                len_to_isp.push_back(pair<double,int> (len,int(fixedXp.size()-1)));
              }
              else if (snap_to_isp1) {

                if (snap_to_isp0 && ((sp_flag[isp1] == 0) || (sp_flag[isp0] == 0))) {
                  // both nodes are snapped to, which indicates the edge length is much smaller than the effective delta
                  // if either one of these has already been snapped to, we should not do anything
                  // if neither snapped, then fall into below loop which snaps us to isp1 (chosen convention)

                  if (sp_flag[isp1] == 0) {
                    // isp1 is an already protected node, so add its new node to the len search
                    map<int,int>::iterator mit=ispToFixedXp.find(isp1);
                    assert(mit!=ispToFixedXp.end());
                    len_to_isp.push_back(pair<double,int> (len,mit->second));
                  }
                }
                else {
                  // pretty much on isp1
                  if (sp_flag[isp1] == 1) {
                    // only add if a link-link node that hasn't already been added
                    fixedXp.push_back(NewNode());
                    FOR_I3 fixedXp.back().xsp[i] = xsp[isp1][i];
                    fixedXp.back().isp = isp1;
                    ispToFixedXp.insert(pair<int,int> (isp1,fixedXp.size()-1));
                    --sp_flag[isp1];  // set as visited
                    to_next -= (my_length-my_traversed);  // compensate for extra movement
                    len += (my_length-my_traversed);
                    len_to_isp.push_back(pair<double,int> (len,int(fixedXp.size()-1)));
                  }
                  else {
                    // is an already protected node, add its new node
                    map<int,int>::iterator mit=ispToFixedXp.find(isp1);
                    assert(mit!=ispToFixedXp.end());
                    len_to_isp.push_back(pair<double,int> (len,mit->second));
                  }
                }
              }
              else {
                assert(snap_to_isp0);
                // pretty much on isp0
                if (sp_flag[isp0] == 1) {
                  // only add if a link-link node that hasn't already been added
                  fixedXp.push_back(NewNode());
                  FOR_I3 fixedXp.back().xsp[i] = xsp[isp0][i];
                  fixedXp.back().isp = isp0;
                  ispToFixedXp.insert(pair<int,int> (isp0,fixedXp.size()-1));
                  --sp_flag[isp0];  // set as visited
                  to_next += my_traversed;  // compensate for extra movement
                  len -= my_traversed;
                  len_to_isp.push_back(pair<double,int> (len,int(fixedXp.size()-1)));
                }
                else {
                  // is an already protected node, add its new node
                  map<int,int>::iterator mit=ispToFixedXp.find(isp0);
                  assert(mit!=ispToFixedXp.end());
                  len_to_isp.push_back(pair<double,int> (len,mit->second));
                }
              }
            }
            else {
              to_next -= (my_length-my_traversed);
              len += (my_length-my_traversed);
              my_traversed = my_length;
            }
          }

          if (chain_edges.end()-link == 1) {
            // sometimes is not caught in the isp1 catch above, so manually register
            // if the last map entry isn't this node (based only on len for now...)
            if (len_to_isp.back().first != len) {
              map<int,int>::iterator mit=ispToFixedXp.find(isp1);
              assert(mit!=ispToFixedXp.end());
              len_to_isp.push_back(pair<double,int> (len,mit->second));
            }
          }
        }

        // second pass to determine isp_to_fixed for non-protected old nodes
        // protected state is present in sp_flag; if 0 then protected
        // have to do in second pass to sensure flagged state of nodes is known
        // due to isp0 snapping that may occur after each link is processed
        len = 0.0;
        for (vector<int>::iterator link=chain_edges.begin(); link!=chain_edges.end(); ++link) {
          int ise = *link;
          int isp0,isp1;
          //bool reverse = false;
          if (ise < 0) {
            ise = -ise-1;
            isp1 = featureEdgeVec[ise].first;
            isp0 = featureEdgeVec[ise].second;
            //reverse = true;
          }
          else {
            isp0 = featureEdgeVec[ise].first;
            isp1 = featureEdgeVec[ise].second;
          }

          const double dx[3] = DIFF(xsp[isp1],xsp[isp0]);
          len += MAG(dx);

          // only add isp1 for each edge if not protected
          // know that chain start isp0 is protected by construction, so only isp1's required
          if (sp_flag[isp1] == 1) len_to_isp.push_back(pair<double,int> (len,-isp1-1));
        }
        chain_edges.clear();


        // all links processed
        sort(len_to_isp.begin(),len_to_isp.end());  // sorts based on arc-length

        if (false) {

          // diagnsotic output for fixed/old nodes
          char filename[128];
          sprintf(filename,"chain_sorted.%03d.dat",n_chains-1);
          FILE * fp = fopen(filename,"w");
          fprintf(fp,"VARIABLES=\"X\" \"Y\" \"Z\" \"perc_len\"\n");
          fprintf(fp,"ZONE T=\"old_not_protected_%d\"\n",n_chains-1);
          for (vector<pair<double,int> >::iterator it=len_to_isp.begin(); it!=len_to_isp.end(); ++it) {
            if (it->second < 0) {
              // old surface node
              const int isp = -it->second-1;
              fprintf(fp,"%18.15e %18.15e %18.15e %18.15e\n",xsp[isp][0],xsp[isp][1],xsp[isp][2],it->first/len);
            }
          }
          fprintf(fp,"ZONE T=\"new_nodes_%d\"\n",n_chains-1);
          for (vector<pair<double,int> >::iterator it=len_to_isp.begin(); it!=len_to_isp.end(); ++it) {
            if (it->second >= 0) {
              // new  node
              fprintf(fp,"%18.15e %18.15e %18.15e %18.15e\n",fixedXp[it->second].xsp[0],fixedXp[it->second].xsp[1],fixedXp[it->second].xsp[2],it->first/len);

            }
          }
          fclose(fp);
        }

        // for each old surface node, determine the nearest new node
        const int n_len_to_isp = len_to_isp.size();
        for (vector<pair<double,int> >::const_iterator it0=len_to_isp.begin(); it0!=len_to_isp.end(); ++it0) {
          if (it0->second < 0) {
            // add entry for this node
            const int isp = -it0->second-1;
            isp_to_nearest_fixed.push_back(pair<int,int> (isp,-1));
            // look forward and backward for nearest new node based on arc length
            int index = it0-len_to_isp.begin();
            while ((len_to_isp[index].second < 0) && (index<n_len_to_isp)) ++index;
            // if (index>=n_len_to_isp) {
            //   for (vector<pair<double,int> >::const_iterator iterr=len_to_isp.begin(); iterr!=len_to_isp.end(); ++iterr) cout << "len,isp: " << iterr->first << " " << iterr->second << endl;
            //   cout << "bad index: " << index << " " << n_len_to_isp << endl;
            // }
            // assert(index>=0 && index<n_len_to_isp);
            // sometimes len computation puts new node at end fo sort, so ignore and take the lower value
            double dist;
            if (index==n_len_to_isp) dist = HUGE_VAL;
            else dist = len_to_isp[index].first - it0->first;

            isp_to_nearest_fixed.back().second = len_to_isp[index].second;

            index = it0-len_to_isp.begin();
            while ((len_to_isp[index].second < 0) && (index>=0)) --index;
            // assert(index>=0 && index<n_len_to_isp);
            if (index>=0 && ((it0->first - len_to_isp[index].first) < dist) ) {  // index criterion ensurethat if sort put new node at bottom, we took higher one
              isp_to_nearest_fixed.back().second = len_to_isp[index].second;
            }

            // cout << "isp to fixed new node for non-protected (isp,new): " << isp_to_nearest_fixed.back().first << " " << isp_to_nearest_fixed.back().second << endl;
          }
        }
      }
    }
  }
  COUT2(" > feature edge nodes:  " << int(fixedXp.size())-n_endpoints);
  COUT2(" > minimum edge delta used:  " << min_effective_delta);
  return min_effective_delta;
}

int SimpleSurface::groupFlaggedTrisByFeatures(const int (*stost)[3]) {

  // currently st_flag stores bits that indicate which edges are features
  // in terms of passing to group processing, this is redundant with stost
  // for the edges that matter (boundaries)

  // here we will over-write st_flag with the disparate flagged group these
  // tris belong to

  int n_groups = 0;
  stack<int> next_tris;
  int * st_bits = new int[nst];
  for (int ist = 0; ist < nst; ++ist) {
    st_bits[ist] = st_flag[ist];
    st_flag[ist] = -1;
  }

  // use 5th bit to set as visited (1-4 in use)
  for (int ist = 0; ist < nst; ++ist) {
    if ((st_bits[ist] >> 0) & 1) {
      assert(zone_flag[znost[ist]] == 1);
      // non zero means part of flagged region

      if (!((st_bits[ist] >> 4) & 1)) {
        // hasn't been processed yet so add to stack
        next_tris.push(ist);

        // loop through neighbors and add those who are not separated by a feature boundary
        while (!next_tris.empty()) {
          const int my_ist = next_tris.top(); next_tris.pop();
          st_flag[my_ist] = n_groups;
          st_bits[my_ist] |= (1 << 4);  // set as visited

          FOR_I3 {
            const int ist_nbr = stost[my_ist][i];
            if (ist_nbr >= 0) {
              const bool edge_is_feature = ((st_bits[my_ist] >> (i+1)) & 1) ? true:false;
              const bool nbr_valid = ((st_bits[ist_nbr] >> 0) & 1) ? true:false;
              const bool nbr_visited = ((st_bits[ist_nbr] >> 4) & 1) ? true:false;
              if (nbr_valid && !nbr_visited && !edge_is_feature) {
                next_tris.push(ist_nbr);
              }
            }
          }
        }
        ++n_groups;
      }

    }
  }

  DELETE(st_bits);

  return n_groups;
}

void SimpleSurface::sampleTriUniformRandom(vector<NewNode>& v_samples,const int n_samples,const int ist,const double delta) const {
  // uniform random distribution on a tri...
  for (int is=0; is < n_samples; ++is) {
    const double r0 = 0.999*double(rand())/double(RAND_MAX)+0.0005;
    const double r1 = 0.999*double(rand())/double(RAND_MAX)+0.0005;
    v_samples.push_back(NewNode());
    FOR_I3 v_samples.back().xsp[i] = (1.0-sqrt(r0))*xsp[spost[ist][0]][i] +
    sqrt(r0)*(1.0-r1)*xsp[spost[ist][1]][i] +
    sqrt(r0)*r1*xsp[spost[ist][2]][i];

    v_samples.back().isp = ist;
    v_samples.back().data_d = delta*1.25;
  }
}

bool sortDescendD2 (const pair<double,int>& i,const pair<double,int>& j) {
  return (i.first>j.first);
}

void SimpleSurface::sampleTriR2(vector<NewNode>& v_samples,const int n_samples,const int ist,const double delta,const int useEdge) const {
  // R^2 method for quasi-random sampling

  // want to find the largest angle node and work from there, so first
  // determine largest edge
  vector<pair<double,int> > edge_d2;
  for (int e=0; e<3; ++e) {
    const double dist2 = DIST2(xsp[spost[ist][e]],xsp[spost[ist][(e+1)%3]]);
    edge_d2.push_back(pair<double,int>(dist2,e));
  }
  sort(edge_d2.begin(), edge_d2.end(), sortDescendD2);  // descending order so largest edge is first

  // largest edge opposite node
  const int inode = (edge_d2[useEdge].second+2)%3;
  // cout << " > using node: " << inode << " at " << COUT_VEC(xsp[spost[ist][inode]]) << " opposite edge with length: " << sqrt(edge_d2[useEdge].first) << endl;

  // basis vectors for the parallelogram
  const double r1[3] = DIFF(xsp[spost[ist][(inode+1)%3]],xsp[spost[ist][inode]]);
  const double r2[3] = DIFF(xsp[spost[ist][(inode+2)%3]],xsp[spost[ist][inode]]);
  // cout << " > r1: " << COUT_VEC(r1) << ", r2: " << COUT_VEC(r2) << endl;

  // could make these defines....
  const double one_o_g = 1.0/1.32471795724474602596;
  const double one_o_g2 = one_o_g*one_o_g;

  for (int is=1; is <= n_samples; ++is) {
    // use quasi-random R^2 method
    double r2_coeff[2];
    r2_coeff[0] = fmod(0.5 + double(is)*one_o_g,1.0);
    r2_coeff[1] = fmod(0.5 + double(is)*one_o_g2,1.0);

    // cout << "point[" << ip_sa << "] r coeffs: " << r2_coeff[0] << ", " << r2_coeff[1] << endl;

    // reflection technique
     if ((r2_coeff[0]+r2_coeff[1]) > 1.0) {
      FOR_I2 r2_coeff[i] = 1.0 - r2_coeff[i];
    }

    v_samples.push_back(NewNode());
    FOR_I3 v_samples.back().xsp[i] = xsp[spost[ist][inode]][i] + r2_coeff[0]*r1[i] + r2_coeff[1]*r2[i];

    v_samples.back().isp = ist;
    v_samples.back().data_d = delta*1.25;
  }
}

void SimpleSurface::sampleTriR2AR(vector<NewNode>& v_samples,const int n_samples,const int ist,const double delta,const int useEdge) const {
  // R^2 method for quasi-random sampling

  // want to find the largest angle node and work from there, so first
  // determine largest edge
  vector<pair<double,int> > edge_d2;
  for (int e=0; e<3; ++e) {
    const double dist2 = DIST2(xsp[spost[ist][e]],xsp[spost[ist][(e+1)%3]]);
    edge_d2.push_back(pair<double,int>(dist2,e));
  }
  sort(edge_d2.begin(), edge_d2.end(), sortDescendD2);  // descending order so largest edge is first

  // largest edge opposite node
  const int inode = (edge_d2[useEdge].second+2)%3;
  // cout << " > largest angle node: " << inode << " opposite edge[" << l_edge << "] with length: " << sqrt(edge_d2) << endl;

  // basis vectors for the parallelogram
  const double r1[3] = DIFF(xsp[spost[ist][(inode+1)%3]],xsp[spost[ist][inode]]);
  const double r2[3] = DIFF(xsp[spost[ist][(inode+2)%3]],xsp[spost[ist][inode]]);
  const double ar=edge_d2[(useEdge+1)%3].first/edge_d2[(useEdge-1)%3].first;  // ratio of these edge lengths

  // cout << " > r1: " << COUT_VEC(r1) << ", r2: " << COUT_VEC(r2) << endl;

  // could make these defines....
  const double one_o_g = 1.0/1.32471795724474602596;
  const double one_o_g2 = one_o_g*one_o_g;

  for (int is=1; is <= n_samples; ++is) {
    // use quasi-random R^2 method
    double r2_coeff[2];
    r2_coeff[0] = fmod(0.5 + double(is)*one_o_g,1.0);
    r2_coeff[1] = fmod(0.5 + double(is)*one_o_g2,1.0);

    // cout << "point[" << ip_sa << "] r coeffs: " << r2_coeff[0] << ", " << r2_coeff[1] << endl;

    // reflection technique
     if ((r2_coeff[0]+r2_coeff[1]) > 1.0) {
      FOR_I2 r2_coeff[i] = 1.0 - r2_coeff[i];
    }

    v_samples.push_back(NewNode());
    FOR_I3 v_samples.back().xsp[i] = xsp[spost[ist][inode]][i] + r2_coeff[0]*r1[i] + r2_coeff[1]*r2[i];

    v_samples.back().isp = ist;
    v_samples.back().data_d = delta*1.25;
  }
}

double SimpleSurface::computeTriAspectRatio(const int ist) const {
  // want to find the largest angle node and work from there, so first
  // determine largest edge
  double edge_d2[2] = {HUGE_VAL,-HUGE_VAL};
  for (int e=0; e<3; ++e) {
    const double tmp_d2 = DIST2(xsp[spost[ist][e]],xsp[spost[ist][(e+1)%3]]);
    edge_d2[0] = min(tmp_d2,edge_d2[0]);
    edge_d2[1] = max(tmp_d2,edge_d2[1]);
  }
  return sqrt(edge_d2[1])/sqrt(edge_d2[0]);
}

int SimpleSurface::init2DVorPointsOnGroup(vector<NewNode>& internal_nodes,const double delta,const int np_fixed,const double factor,const int igr, const bool b_keep_edge_nodes,const int init_type) const {

  // store xsp, ist (in isp), and delta as data_d
  // in internal_nodes

  // based on count of fixed nodes we may want to insert much fewer nodes than based
  // purely on the tri surface areas
  // assuming each boundary node covers ~0.5 an internal, we remove this area from the calculation
  const double area_vor = factor*delta*delta;

  double fixed_area = HUGE_VAL;
  if (b_keep_edge_nodes) fixed_area = 0.0;
  else fixed_area = double(np_fixed)*0.5*area_vor;

  double area_st = 0.0;  // running flagged surface area
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] == igr) {
      const double normal2[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
      const double this_area = 0.5*MAG(normal2);
      area_st += this_area;
    }
  }
  const double group_area = area_st;
  COUT2(" > group surface area: " << group_area);


  if (fixed_area >= group_area) {
    COUT2(" > introducing 0 surface vor points on this group");
    return 0;  // no internal nodes added
  }

  // use a factor to keep processing all tris but seed more sparsely
  const double area_ratio = (group_area-fixed_area)/group_area;

  area_st = 0.0;
  int type_count[2] = {0,0};
  int ip_begin = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] == igr) {
      const double normal2[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
      const double this_area = 0.5*MAG(normal2);
      area_st += this_area;

      const int ip_end = int(area_ratio*area_st/area_vor);
      if (ip_end > ip_begin) {
        // this tri gets some points...
        const int n_samples = ip_end - ip_begin;
        if (init_type == 0) {
          sampleTriUniformRandom(internal_nodes,n_samples,ist,delta);
          ++type_count[0];
        }
        else if (init_type == 1) {
          sampleTriR2(internal_nodes,n_samples,ist,delta);
          ++type_count[1];
        }
        else if (init_type == 2) {
          // mixed depending on n_samples and aspect ratio
          const double ar = computeTriAspectRatio(ist);
          const double n_o_ar = double(n_samples)/ar;
          if ((n_samples < 5) || (n_o_ar < 2.0)) {
            sampleTriUniformRandom(internal_nodes,n_samples,ist,delta);
            ++type_count[0];
          }
          else {
            // even if skewed it seems as though this approach isn't bad; simply a density problem if too few points
            sampleTriR2(internal_nodes,n_samples,ist,delta,2);
            ++type_count[1];
          }
        }
        // else if (init_type == 3) {
        //   sampleTriR2(internal_nodes,n_samples,ist,delta,2);
        //   ++type_count[1];
        // }
        // else if (init_type == 4) {
        //   sampleTriR2AR(internal_nodes,n_samples,ist,delta,0);
        //   ++type_count[1];
        // }
        // else if (init_type == 5) {
        //   sampleTriR2AR(internal_nodes,n_samples,ist,delta,2);
        //   ++type_count[1];
        // }
      }
      ip_begin = ip_end;
    }
  }
  assert(ip_begin == int(internal_nodes.size()));

  const int n_internal = internal_nodes.size();
  COUT2(" > introducing " << n_internal << " surface vor points on this group");
  COUT2(" > initial condition type tri count; random: " << type_count[0] << ", R2: " << type_count[1]);

  return n_internal;
}

void edgeFlipPreserveTeost(Teost& my_teost, int (*spost_new)[3], const int ist, const int i) {
  int ist_nbr;
  int i_nbr;
  int orient;
  if(my_teost.getTriNbrData(ist_nbr,i_nbr,orient, ist,i)){

    // TmpTri tri1;
    // TmpTri tri2;
    NewTri tri1;
    NewTri tri2;
    // tri1.ist=ist;
    // tri2.ist=ist_nbr;
    tri1.spost[0]=spost_new[ist][(i+2)%3];
    tri1.spost[1]=spost_new[ist_nbr][(i_nbr+2)%3];
    tri1.spost[2]=spost_new[ist_nbr][i_nbr];
    tri2.spost[0]=spost_new[ist_nbr][(i_nbr+2)%3];
    tri2.spost[1]=spost_new[ist][(i+2)%3];
    tri2.spost[2]=spost_new[ist][i];

    spost_new[ist][0]=tri1.spost[0];
    spost_new[ist][1]=tri1.spost[1];
    spost_new[ist][2]=tri1.spost[2];
    spost_new[ist_nbr][0]=tri2.spost[0];
    spost_new[ist_nbr][1]=tri2.spost[1];
    spost_new[ist_nbr][2]=tri2.spost[2];

    int ist_nbr2[4],i_nbr2[4], orient2[4]; bool valid[4];
    // get external neighbour data
    valid[0]=my_teost.getTriNbrData(ist_nbr2[0],i_nbr2[0],orient2[0], ist_nbr,(i_nbr+2)%3);
    valid[1]=my_teost.getTriNbrData(ist_nbr2[1],i_nbr2[1],orient2[1], ist,(i+1)%3);
    valid[2]=my_teost.getTriNbrData(ist_nbr2[2],i_nbr2[2],orient2[2], ist,(i+2)%3);
    valid[3]=my_teost.getTriNbrData(ist_nbr2[3],i_nbr2[3],orient2[3], ist_nbr,(i_nbr+1)%3);

    // set neigbours for  the new ist
    if(valid[0]){
      my_teost.setTriNbrData(ist,1,orient2[0],ist_nbr2[0],i_nbr2[0]);
      my_teost.setTriNbrData(ist_nbr2[0],i_nbr2[0],orient2[0],ist,1);
    }else{
      my_teost.setTriNbrOpen(ist, 1);
    }
    if(valid[1]){
      my_teost.setTriNbrData(ist,2,orient2[1],ist_nbr2[1],i_nbr2[1]);
      my_teost.setTriNbrData(ist_nbr2[1],i_nbr2[1],orient2[1],ist,2);
    }else{
      my_teost.setTriNbrOpen(ist, 2);
    }
    // 	 // set neigbours for new ist_nbr
    if(valid[2]){
      my_teost.setTriNbrData(ist_nbr,1,orient2[2],ist_nbr2[2],i_nbr2[2]);
      my_teost.setTriNbrData(ist_nbr2[2],i_nbr2[2],orient2[2],ist_nbr,1);
    }else{
      my_teost.setTriNbrOpen(ist_nbr, 1);
    }

    if(valid[3]){
      my_teost.setTriNbrData(ist_nbr,2,orient2[3],ist_nbr2[3],i_nbr2[3]);
      my_teost.setTriNbrData(ist_nbr2[3],i_nbr2[3],orient2[3],ist_nbr,2);
    }else{
      my_teost.setTriNbrOpen(ist_nbr, 2);
    }

    // 	 //the external neihbours have been set correctly
    my_teost.setTriNbrData(ist,0,orient,ist_nbr,0);
    my_teost.setTriNbrData(ist_nbr,0,orient,ist,0);
  }else{
    assert(1==2); // Can't flip if there is no valid neighbour
  }
}


void makeTrisDelaunay(Teost& my_teost,int (*spost_new)[3],const int nst_new,const double (*xsp_new)[3]) {

  double (*angle)[3];
  angle=new double[nst_new][3]; // each edge needs an angle value

  for (int ist=0; ist<nst_new; ++ist) {
    FOR_I3 {
      int isp0=spost_new[ist][i];
      int isp1=spost_new[ist][(i+1)%3];
      int isp2=spost_new[ist][(i+2)%3];
      double x10[3];
      double x12[3];
      FOR_J3 {
        x10[j]=xsp_new[isp0][j]-xsp_new[isp1][j];
        x12[j]=xsp_new[isp2][j]-xsp_new[isp1][j];
      }
      double d0=DIST(xsp_new[isp0],xsp_new[isp1]);
      double d1=DIST(xsp_new[isp1],xsp_new[isp2]);
      double dp=DOT_PRODUCT(x10,x12);
      dp=acos(dp/(d0*d1));
      angle[ist][(i+2)%3]=dp;
    }
  }

  // remove this later
  int counter=0;
  int t_counter=0;
  for (int ist=0; ist<nst_new; ++ist) {
    FOR_I3 {
      t_counter++;
      int ist_nbr;
      int i_nbr;
      int orient_nbr;
      if (my_teost.getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
        //cout << angle<<endl;
        if ((angle[ist][i]+angle[ist_nbr][i_nbr])>M_PI) {
          counter++;
        }
      }
    }
  }
  cout << " > non Delaunay edges found: "<<counter << " of "<< t_counter << " total edges processed" << endl;

  // this implementation is likely far from optimal but it doesn't seem to slow the process down by much.
  bool done=false;
  bool converged=false;
  int iter=0;
  cout << " > starting Delaunay edge flipping iterations" <<endl;
  while(!done){
    int n_flips=0;
    iter++;
    for (int ist=0; ist<nst_new; ++ist) {
      FOR_I3{
        int ist_nbr;
        int i_nbr;
        int orient;
        if(my_teost.getTriNbrData(ist_nbr,i_nbr,orient, ist,i)){
          if((angle[ist][i]+angle[ist_nbr][i_nbr])>M_PI){
            // we need to flip
            edgeFlipPreserveTeost(my_teost,spost_new,ist,i);
            n_flips=n_flips+1;

            //recompute angles
            FOR_K3 {
              int isp0=spost_new[ist][k];
              int isp1=spost_new[ist][(k+1)%3];
              int isp2=spost_new[ist][(k+2)%3];
              double x10[3];
              double x12[3];
              FOR_J3 {
                x10[j]=xsp_new[isp0][j]-xsp_new[isp1][j];
                x12[j]=xsp_new[isp2][j]-xsp_new[isp1][j];
              }
              double d0=DIST(xsp_new[isp0],xsp_new[isp1]);
              double d1=DIST(xsp_new[isp1],xsp_new[isp2]);
              double dp=DOT_PRODUCT(x10,x12);
              dp=acos(dp/(d0*d1));
              angle[ist][(k+2)%3]=dp;

            }
            FOR_K3 {
              int isp0=spost_new[ist_nbr][k];
              int isp1=spost_new[ist_nbr][(k+1)%3];
              int isp2=spost_new[ist_nbr][(k+2)%3];
              double x10[3];
              double x12[3];
              FOR_J3 {
                x10[j]=xsp_new[isp0][j]-xsp_new[isp1][j];
                x12[j]=xsp_new[isp2][j]-xsp_new[isp1][j];
              }
              double d0=DIST(xsp_new[isp0],xsp_new[isp1]);
              double d1=DIST(xsp_new[isp1],xsp_new[isp2]);
              double dp=DOT_PRODUCT(x10,x12);
              dp=acos(dp/(d0*d1));
              angle[ist_nbr][(k+2)%3]=dp;

            }
          }
        }
      }
    }
    cout << "    > iteration [" << iter << "]: edges flipped: "<< n_flips <<endl;
    if(n_flips==0){
      done=true;
      converged=true;
    }else if(iter==20){
      done=true;
      converged=false;
    }

  }

  DELETE(angle);

  if(converged){
    cout<< " > surface mesh is Delaunay"<<endl;
  }
  else {
    cout << " ! Warning: Delauney surface mesh not achieved in "<<iter<< " iterations" << endl;
  }
}

void SimpleSurface::makeSelectedTrisDelaunay(const int iter_limit) {
  ensureTeost();

  // each edge needs an angle value...

  double (*angle)[3];
  angle = new double[nst][3];
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      int isp0 = spost[ist][i];
      int isp1 = spost[ist][(i+1)%3];
      int isp2 = spost[ist][(i+2)%3];
      double x10[3];
      double x12[3];
      FOR_J3 {
        x10[j] = xsp[isp0][j]-xsp[isp1][j];
        x12[j] = xsp[isp2][j]-xsp[isp1][j];
      }
      double d0 = DIST(xsp[isp0],xsp[isp1]);
      double d1 = DIST(xsp[isp1],xsp[isp2]);
      double dp = DOT_PRODUCT(x10,x12);
      dp = acos(dp/(d0*d1));
      angle[ist][(i+2)%3] = dp;
    }
  }

  /*
  // Loop through and see if any of our nbrs are not flagged.
  // Their edges might require splitting...

  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] != 0) {
      FOR_I3 {
        int ist_nbr,i_nbr,orient_nbr;
        if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i))
          st_flag[ist_nbr] = st_flag[ist];
      }
    }
  }

  // used for memory management...

  int nsp_max = nsp;
  int nst_max = nst;

  map<const pair<int,int>,int> newNodeMap;

  // we change (grow) nst, so store nst_old...
  const int nst_old = nst;

  int nsplit = 0;
  const double split_cos = cos((180.0-60.0)*M_PI/180.0);
  for (int ist = 0; ist < nst_old; ++ist) {
    if (st_flag[ist] != 0) {
      FOR_I3 {
        if (isEdgeCrease(ist,i)&&(edgeCreaseCos(ist,i) < split_cos)) {
          int ist_nbr;
          int i_nbr;
          int orient;
          if (teost.getTriNbrData(ist_nbr,i_nbr,orient,ist,i)) {
            cout << edgeCreaseCos(ist,i) << " " << split_cos << " " << angle[ist][i]+angle[ist_nbr][i_nbr] << endl;
            if ( (angle[ist][i]+angle[ist_nbr][i_nbr]) > M_PI ) {
              ++nsplit;

              // split this tri at this edge...
              const int isp0 = spost[ist][i];
              const int isp1 = spost[ist][(i+1)%3];
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
                szost.ensureSize(nst_max);
                double (*angle0)[3] = angle;
                angle = new double[nst_max][3];
                for (int ist = 0; ist < nst; ++ist) FOR_I3 angle[ist][i] = angle0[ist][i];
                delete[] angle0;
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
              spost[ist_new][i]       = spost[ist][i];
              spost[ist_new][(i+1)%3] = isp_new;
              spost[ist_new][(i+2)%3] = spost[ist][(i+2)%3];
              znost[ist_new] = znost[ist];
              szost[ist_new] = szost[ist];
              // reuse the original tri -- only need to change one spost...
              spost[ist][i] = isp_new;
              // set the bits of the new tri...
              st_flag[ist_new] = st_flag[ist];
              FOR_K3 {
                const int isp0 = spost[ist_new][k];
                const int isp1 = spost[ist_new][(k+1)%3];
                const int isp2 = spost[ist_new][(k+2)%3];
                double x10[3];
                double x12[3];
                FOR_J3 {
                  x10[j] = xsp[isp0][j]-xsp[isp1][j];
                  x12[j] = xsp[isp2][j]-xsp[isp1][j];
                }
                double d0 = DIST(xsp[isp0],xsp[isp1]);
                double d1 = DIST(xsp[isp1],xsp[isp2]);
                double dp = DOT_PRODUCT(x10,x12);
                dp = acos(dp/(d0*d1));
                angle[ist_new][(k+2)%3] = dp;
              }
              FOR_K3 {
                const int isp0 = spost[ist][k];
                const int isp1 = spost[ist][(k+1)%3];
                const int isp2 = spost[ist][(k+2)%3];
                double x10[3];
                double x12[3];
                FOR_J3 {
                  x10[j] = xsp[isp0][j]-xsp[isp1][j];
                  x12[j] = xsp[isp2][j]-xsp[isp1][j];
                }
                double d0 = DIST(xsp[isp0],xsp[isp1]);
                double d1 = DIST(xsp[isp1],xsp[isp2]);
                double dp = DOT_PRODUCT(x10,x12);
                dp = acos(dp/(d0*d1));
                angle[ist][(k+2)%3] = dp;
              }
            }
          }
        }
      }
    }
  }
  cout << " > split feature edges: " << nsplit/2 << endl;

  // we added new tris, so clear some stuff that might be set...

  clearTeost();
  clearStoszSzosz();
  buildTeost();
  */

  // counting loop; remove this later
  int counter=0;
  int t_counter=0;
  for (int ist=0; ist<nst; ++ist) {
    if (st_flag[ist] != 0) {
      FOR_I3 {
        t_counter++;
        int ist_nbr;
        int i_nbr;
        int orient_nbr;
        if (teost.getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
          if ((angle[ist][i]+angle[ist_nbr][i_nbr])>M_PI) {
            counter++;
          }
        }
      }
    }
  }
  cout << " > non Delaunay half-edges found: "<<counter << " of "<< t_counter << " total half-edges" << endl;

  // this implementation is likely far from optimal but it doesn't seem to slow the process down by much.
  bool done=false;
  bool converged=false;
  int iter=0;
  cout << " > starting Delaunay edge flipping iterations" <<endl;
  while(!done){
    int n_flips=0;
    iter++;
    for (int ist=0; ist<nst; ++ist) {
      if (st_flag[ist]) {
        FOR_I3{
          if (isEdgeCrease(ist,i) || isEdgeSubzoneBoundary(ist,i)) {
            continue;  // skip processing of this edge so we don't break a feature to preserve
          }

          int ist_nbr;
          int i_nbr;
          int orient;
          if(teost.getTriNbrData(ist_nbr,i_nbr,orient, ist,i)){
            if((angle[ist][i]+angle[ist_nbr][i_nbr])>M_PI){
              // we need to flip
              edgeFlipPreserveTeost(teost,spost,ist,i);
              n_flips=n_flips+1;

              //recompute angles
              FOR_K3 {
                int isp0=spost[ist][k];
                int isp1=spost[ist][(k+1)%3];
                int isp2=spost[ist][(k+2)%3];
                double x10[3];
                double x12[3];
                FOR_J3 {
                  x10[j]=xsp[isp0][j]-xsp[isp1][j];
                  x12[j]=xsp[isp2][j]-xsp[isp1][j];
                }
                double d0=DIST(xsp[isp0],xsp[isp1]);
                double d1=DIST(xsp[isp1],xsp[isp2]);
                double dp=DOT_PRODUCT(x10,x12);
                dp=acos(dp/(d0*d1));
                angle[ist][(k+2)%3]=dp;

              }
              FOR_K3 {
                int isp0=spost[ist_nbr][k];
                int isp1=spost[ist_nbr][(k+1)%3];
                int isp2=spost[ist_nbr][(k+2)%3];
                double x10[3];
                double x12[3];
                FOR_J3 {
                  x10[j]=xsp[isp0][j]-xsp[isp1][j];
                  x12[j]=xsp[isp2][j]-xsp[isp1][j];
                }
                double d0=DIST(xsp[isp0],xsp[isp1]);
                double d1=DIST(xsp[isp1],xsp[isp2]);
                double dp=DOT_PRODUCT(x10,x12);
                dp=acos(dp/(d0*d1));
                angle[ist_nbr][(k+2)%3]=dp;

              }
            }
          }
        }
      }
    }
    cout << "    > iteration [" << iter << "]: edges flipped: "<< n_flips <<endl;
    if(n_flips==0){
      done=true;
      converged=true;
    }else if(iter==iter_limit){
      done=true;
      converged=false;
    }

  }

  DELETE(angle);

  if(converged){
    cout<< " > surface mesh is Delaunay"<<endl;
  }
  else {
    cout << " ! Warning: Delauney surface mesh not achieved in "<<iter<< " iterations" << endl;
  }
}

int SimpleSurface::triangulateVorPoints(vector<NewTri>& triVec,const double (*xp)[3],const int * stoxp,const int * nboxp_i,const vector<int>& nboxp_v,const int np) const {
  int n_tris = 0;
  for (int ip0 = 0; ip0 < np; ++ip0) {
    for (int nox0 = nboxp_i[ip0]; nox0 != nboxp_i[ip0+1]; ++nox0) {
      const int ip1 = nboxp_v[nox0];
      if (ip1 > ip0) {
        int nox1;
        for (nox1 = nboxp_i[ip1]; nox1 != nboxp_i[ip1+1]; ++nox1) {
          if (nboxp_v[nox1] == ip0)
            break;
        }
        if (nox1 != nboxp_i[ip1+1]) {
          for (nox1 = nboxp_i[ip1]; nox1 != nboxp_i[ip1+1]; ++nox1) {
            const int ip2 = nboxp_v[nox1];
            if (ip2 > ip1) {
              // ip0 has ip1,
              // ip1 has ip0 and ip2
              // first, confirm that ip0 has ip2...
              int nox0_;
              for (nox0_ = nboxp_i[ip0]; nox0_ != nboxp_i[ip0+1]; ++nox0_) {
                if (nboxp_v[nox0_] == ip2)
                  break;
              }
              if (nox0_ != nboxp_i[ip0+1]) {
                // finally, check for ip0 and ip1 in ip2's nbrs...
                int count = 0;
                for (int nox2 = nboxp_i[ip2]; nox2 != nboxp_i[ip2+1]; ++nox2) {
                  if ((nboxp_v[nox2] == ip0)||(nboxp_v[nox2] == ip1)) {
                    ++count;
                    if (count == 2) {
                      const double this_normal2[3] = TRI_NORMAL_2(xp[ip0],xp[ip1],xp[ip2]);
                      const int ist0 = stoxp[ip0];
                      const double st0_normal2[3] = TRI_NORMAL_2(xsp[spost[ist0][0]],xsp[spost[ist0][1]],xsp[spost[ist0][2]]);
                      const double st0_normal2_mag = MAG(st0_normal2); assert(st0_normal2_mag > 0.0);
                      const int ist1 = stoxp[ip1];
                      const double st1_normal2[3] = TRI_NORMAL_2(xsp[spost[ist1][0]],xsp[spost[ist1][1]],xsp[spost[ist1][2]]);
                      const double st1_normal2_mag = MAG(st1_normal2); assert(st1_normal2_mag > 0.0);
                      const int ist2 = stoxp[ip2];
                      const double st2_normal2[3] = TRI_NORMAL_2(xsp[spost[ist2][0]],xsp[spost[ist2][1]],xsp[spost[ist2][2]]);
                      const double st2_normal2_mag = MAG(st2_normal2); assert(st2_normal2_mag > 0.0);
                      // just sum the dot products here (should we normalize? sts can be any size)...
                      const double dp =
                        DOT_PRODUCT(this_normal2,st0_normal2)/st0_normal2_mag +
                        DOT_PRODUCT(this_normal2,st1_normal2)/st1_normal2_mag +
                        DOT_PRODUCT(this_normal2,st2_normal2)/st2_normal2_mag;
                      if (dp > 0.0) {
                        triVec.push_back(NewTri(ip0,ip1,ip2));
                        ++n_tris;
                      }
                      else {
                        triVec.push_back(NewTri(ip0,ip2,ip1));
                        ++n_tris;
                      }
                      break;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return n_tris;
}

void SimpleSurface::updateSurfaceWithRetessellation(const vector<NewNode>& new_nodes,const int np_fixed,const vector<NewTri>& new_adj_tris,const vector<NewTri>& new_tris) {
  // flagged tris (st_flag) should contain retessellation group information
  // so positive group indices are tris that will be removed

  // need map from current node value to new node value
  const int nsp_all = nsp + int(new_nodes.size());
  sp_flag.resize(nsp_all);
  sp_flag.setAll(-1);

  int sp_count = 0;
  int st_count = 0;
  vector<NewNode> global_nodes;
  vector<NewTri> global_tris;

  // loop through all kept tris and register their nodes
  for (int ist=0; ist<nst; ++ist) {
    if (st_flag[ist] < 0) {
      FOR_I3 {
        const int isp = spost[ist][i];
        if (sp_flag[isp] < 0) {
          sp_flag[isp] = sp_count++;
          global_nodes.push_back(NewNode(xsp[isp]));
        }
      }

      // all node indices have been updated, so can push back this tri with valid indices for new surface
      ++st_count;
      global_tris.push_back(NewTri(sp_flag[spost[ist][0]],sp_flag[spost[ist][1]],sp_flag[spost[ist][2]],znost[ist],szost[ist]));
    }
  }

  // add nodes from adjacent interface tris as well
  for (vector<NewTri>::const_iterator it=new_adj_tris.begin(); it!=new_adj_tris.end(); ++it) {
    ++st_count;

    int new_spost[3] = {-1,-1,-1};
    FOR_I3 {
      if (it->spost[i] < 0) {
        // references an old node
        const int isp = -it->spost[i]-1;
        if (sp_flag[isp] < 0) {
          sp_flag[isp] = sp_count++;
          global_nodes.push_back(NewNode(xsp[isp]));
        }
        new_spost[i] = sp_flag[isp];
      }
      else {
        // references new node
        const int nn_isp = it->spost[i];
        if ((nn_isp < np_fixed) && (new_nodes[nn_isp].isp >= 0)) {
          // collocated with an old surface node, use that isp
          const int isp = new_nodes[nn_isp].isp;
          if (sp_flag[isp] < 0) {
            sp_flag[isp] = sp_count++;
            global_nodes.push_back(NewNode(xsp[isp]));
          }
          new_spost[i] = sp_flag[isp];
        }
        else {
          if (sp_flag[nsp+nn_isp] < 0) {
            sp_flag[nsp+nn_isp] = sp_count++;
            global_nodes.push_back(NewNode(new_nodes[nn_isp].xsp));
          }
          new_spost[i] = sp_flag[nsp+nn_isp];
        }
      }
    }

    // all node indices have been updated, so can push back with valid indices for new surface
    global_tris.push_back(NewTri(new_spost[0],new_spost[1],new_spost[2],it->znost,it->szost));
  }

  // add nodes from new tris as well
  for (vector<NewTri>::const_iterator it=new_tris.begin(); it!=new_tris.end(); ++it) {
    ++st_count;

    int new_spost[3] = {-1,-1,-1};
    FOR_I3 {
      // always references a new_node entry
      const int nn_isp = it->spost[i];

      if ((nn_isp < np_fixed) && (new_nodes[nn_isp].isp >= 0)) {
        // collocated with an old surface node, use that isp
        const int isp = new_nodes[nn_isp].isp;
        if (sp_flag[isp] < 0) {
          sp_flag[isp] = sp_count++;
          global_nodes.push_back(NewNode(xsp[isp]));
        }
        new_spost[i] = sp_flag[isp];
      }
      else {
        if (sp_flag[nsp+nn_isp] < 0) {
          sp_flag[nsp+nn_isp] = sp_count++;
          global_nodes.push_back(NewNode(new_nodes[nn_isp].xsp));
        }
        new_spost[i] = sp_flag[nsp+nn_isp];
      }
    }

    // all node indices have been updated, so can push back with valid indices for new surface
    global_tris.push_back(NewTri(new_spost[0],new_spost[1],new_spost[2],it->znost,it->szost));

  }

  assert(sp_count == int(global_nodes.size()));
  cout << "processing new tris and nodes for updated surface:" << endl;
  cout << " > determined " << sp_count << " total nodes" << endl;
  if (st_count != int(global_tris.size())) {
    WUI(WARN,"not all expected tris were built; some holes may exist where tris are missing at patch edges");
  }
  cout << " > determined " << st_count << " total tris" << endl;

  cout << " > resizing surface...";
  // fill in new nodes
  const int nsp0 = nsp;
  nsp = sp_count;
  resizeNspData(nsp,nsp0);
  for (vector<NewNode>::const_iterator it=global_nodes.begin(); it!=global_nodes.end(); ++it) {
    const int index = it-global_nodes.begin();
    FOR_I3 xsp[index][i] = it->xsp[i];
  }
  global_nodes.clear();

  cout << "...";

  // fill in new tris
  const int nst0 = nst;
  nst = st_count;
  resizeNstData(nst,nst0);
  for (vector<NewTri>::const_iterator it=global_tris.begin(); it!=global_tris.end(); ++it) {
    const int index = it-global_tris.begin();
    FOR_I3 spost[index][i] = it->spost[i];
    znost[index] = it->znost;
    szost[index] = it->szost;
  }
  global_tris.clear();

  cout << "done" << endl;

  // cleanup
  clearStoszSzosz();
  clearTeost();
  clearDynamicEdgeGroups();
  clearNonManifoldData();
  clearSubzoneData();
  clearZoneData();
  pruneEmptyZonesAndSubzones();
}

void SimpleSurface::zipAdjacentTris(vector<NewTri>& new_adj_tris,vector<pair<int,int> >& adj_tris,const multimap <pair<int,int>,int>& edge_to_fixedNodes,const map<int,pair<int,bool> >& isp_to_fixedNodes,vector<NewNode>& fixedNodes) {
  // loop through adjacent tris and determine how to join them to their new fixed nodes
  // old nodes are -1-indexed, new nodes are 0+ indexed
  for (vector<pair<int,int> >::const_iterator it=adj_tris.begin(); it!=adj_tris.end(); ++it) {

    // check feature touching bits to see how to process (node touch, edge touch, all edges)
    const int bits = it->second;
    const int ist = it->first;
    if ((bits == 1) || (bits == 2) || (bits == 4)) {
      // single node is on feature edge
      int new_spost[3];
      FOR_I3 new_spost[i] = -spost[ist][i]-1;  // old nodes -1-indexed
      // update node that sits on features with fixedNode index
      int index = -1;
      switch (bits) {
        case 1:
          index = 0;
          break;
        case 2:
          index = 1;
          break;
        case 4:
          index = 2;
          break;
        default:
          assert(0);
      }
      assert(index != -1);

      // find fixedNode index to use (0-indexed because this tri is in New group)
      map<int,pair<int,bool> >::const_iterator mit = isp_to_fixedNodes.find(spost[ist][index]);
      if (mit!=isp_to_fixedNodes.end()) {
        new_spost[index] = mit->second.first;  // 0-indexed
        new_adj_tris.push_back(NewTri(new_spost[0],new_spost[1],new_spost[2],znost[ist],szost[ist]));
      }
    }
    else if ((bits == 3) || (bits == 5) || (bits == 6)) {
      // single edge requires treatment
      int i0,i1,i2;
      switch (bits) {
        case 3:
          i0 = 0;
          i1 = 1;
          i2 = 2;
          break;
        case 5:
          i0 = 2;
          i1 = 0;
          i2 = 1;
          break;
        case 6:
          i0 = 1;
          i1 = 2;
          i2 = 0;
          break;
        default:
          assert(0);
      }
      // first treat split edges only
      int isp0 = spost[ist][i0];
      int isp1 = spost[ist][i1];
      assert(isp0!=isp1);
      i0 = 0;
      i1 = 1;
      if (isp0 > isp1) {
        swap(isp0,isp1);
        swap(i0,i1);
      }

      // isp0,1 are now ordered from min to max
      map<int,pair<int,bool> >::const_iterator mit = isp_to_fixedNodes.find(isp0);
      assert(mit!=isp_to_fixedNodes.end());
      int isp_min = mit->second.first;

      mit = isp_to_fixedNodes.find(isp1);
      assert(mit!=isp_to_fixedNodes.end());
      int isp_max = mit->second.first;

      vector<pair<double,int> > frac_to_new_isp;
      // add current start/end nodes
      frac_to_new_isp.push_back(pair<double,int> (0.0,isp_min));
      frac_to_new_isp.push_back(pair<double,int> (1.0,isp_max));

      // find splits registered to this edge
      pair<multimap <pair<int,int>,int>::const_iterator,multimap <pair<int,int>,int>::const_iterator> ret;
      ret = edge_to_fixedNodes.equal_range(pair<int,int>(isp0,isp1));
      for (multimap <pair<int,int>,int>::const_iterator it=ret.first; it!=ret.second; ++it) {
        frac_to_new_isp.push_back(pair<double,int> (fixedNodes[it->second].data_d,it->second));
      }
      sort(frac_to_new_isp.begin(),frac_to_new_isp.end());

      // start from isp_min and march along frac
      int isp_last = frac_to_new_isp[0].second;
      const int isp_opposite = -spost[ist][i2]-1;  // old node -1-indexed

      for (int ifrac=1, nfrac=frac_to_new_isp.size(); ifrac<nfrac; ++ifrac) {
        int isp_next = frac_to_new_isp[ifrac].second;

        if (isp_next != isp_last) {
          // only add if not collapsed
          new_adj_tris.push_back(NewTri());
          new_adj_tris.back().spost[i0] = isp_opposite;
          new_adj_tris.back().spost[i1] = isp_last;
          new_adj_tris.back().spost[2]  = isp_next;
          new_adj_tris.back().znost = znost[ist];
          new_adj_tris.back().szost = szost[ist];
          // cout << " + added tri spost: " << COUT_VEC(new_adj_tris.back().spost) << " frac: " << it->first << endl;
        }
        isp_last = isp_next;
      }
    }
    else {
      assert(bits==7);
      // either 2 or 3 edges have splits registered
      // introduce fanning that doesn't require new node...

      // first treat split edges only
      int isp0 = spost[ist][0];
      int isp1 = spost[ist][1];
      int isp2 = spost[ist][2];

      // for each edge find the registered splits
      vector<pair<double,int> > e0_frac_to_new_isp;
      vector<pair<double,int> > e1_frac_to_new_isp;
      vector<pair<double,int> > e2_frac_to_new_isp;
      bool e_invert[3] = {false,false,false};


      if (isp0 > isp1) e_invert[0] = true;
      if (isp1 > isp2) e_invert[1] = true;
      if (isp2 > isp0) e_invert[2] = true;

      map<int,pair<int,bool> >::const_iterator mit = isp_to_fixedNodes.find(isp0);
      assert(mit!=isp_to_fixedNodes.end());
      const int new_isp0 = mit->second.first;

      mit = isp_to_fixedNodes.find(isp1);
      assert(mit!=isp_to_fixedNodes.end());
      const int new_isp1 = mit->second.first;

      mit = isp_to_fixedNodes.find(isp2);
      assert(mit!=isp_to_fixedNodes.end());
      const int new_isp2 = mit->second.first;

      if ( (new_isp0 == new_isp1) || (new_isp1 == new_isp2) || (new_isp0 == new_isp2) ) continue;  // avoid introducing linear tris

      // edge0
      {
        // add current start/end nodes
        int isp_min = (e_invert[0]) ? new_isp1:new_isp0;
        int isp_max = (e_invert[0]) ? new_isp0:new_isp1;
        e0_frac_to_new_isp.push_back(pair<double,int> (0.0,isp_min));
        e0_frac_to_new_isp.push_back(pair<double,int> (1.0,isp_max));

        // find splits registered to this edge
        pair<multimap <pair<int,int>,int>::const_iterator,multimap <pair<int,int>,int>::const_iterator> ret;
        if (e_invert[0]) ret = edge_to_fixedNodes.equal_range(pair<int,int>(isp1,isp0));
        else ret = edge_to_fixedNodes.equal_range(pair<int,int>(isp0,isp1));

        for (multimap <pair<int,int>,int>::const_iterator it=ret.first; it!=ret.second; ++it) {
          e0_frac_to_new_isp.push_back(pair<double,int> (fixedNodes[it->second].data_d,it->second));
        }
        sort(e0_frac_to_new_isp.begin(),e0_frac_to_new_isp.end());
      }

      // edge1
      {
        // add current start/end nodes
        int isp_min = (e_invert[1]) ? new_isp2:new_isp1;
        int isp_max = (e_invert[1]) ? new_isp1:new_isp2;
        e1_frac_to_new_isp.push_back(pair<double,int> (0.0,isp_min));
        e1_frac_to_new_isp.push_back(pair<double,int> (1.0,isp_max));

        // find splits registered to this edge
        pair<multimap <pair<int,int>,int>::const_iterator,multimap <pair<int,int>,int>::const_iterator> ret;
        if (e_invert[1]) ret = edge_to_fixedNodes.equal_range(pair<int,int>(isp2,isp1));
        else ret = edge_to_fixedNodes.equal_range(pair<int,int>(isp1,isp2));

        for (multimap <pair<int,int>,int>::const_iterator it=ret.first; it!=ret.second; ++it) {
          e1_frac_to_new_isp.push_back(pair<double,int> (fixedNodes[it->second].data_d,it->second));
        }
        sort(e1_frac_to_new_isp.begin(),e1_frac_to_new_isp.end());
      }

      // edge2
      {
        // add current start/end nodes
        int isp_min = (e_invert[2]) ? new_isp0:new_isp2;
        int isp_max = (e_invert[2]) ? new_isp2:new_isp0;
        e2_frac_to_new_isp.push_back(pair<double,int> (0.0,isp_min));
        e2_frac_to_new_isp.push_back(pair<double,int> (1.0,isp_max));

        // find splits registered to this edge
        pair<multimap <pair<int,int>,int>::const_iterator,multimap <pair<int,int>,int>::const_iterator> ret;
        if (e_invert[2]) ret = edge_to_fixedNodes.equal_range(pair<int,int>(isp0,isp2));
        else ret = edge_to_fixedNodes.equal_range(pair<int,int>(isp2,isp0));

        for (multimap <pair<int,int>,int>::const_iterator it=ret.first; it!=ret.second; ++it) {
          e2_frac_to_new_isp.push_back(pair<double,int> (fixedNodes[it->second].data_d,it->second));
        }
        sort(e2_frac_to_new_isp.begin(),e2_frac_to_new_isp.end());
      }

      // add new centerpoint node
      const int isp_middle = fixedNodes.size();
      fixedNodes.push_back(NewNode(0.0,0.0,0.0));
      FOR_I3 fixedNodes.back().xsp[i] = (xsp[isp0][i] + xsp[isp1][i] + xsp[isp2][i])/3.0;

      // march and fill
      int index = 0;
      int index_max = e0_frac_to_new_isp.size();
      // cout << "e0 has " << index_max << " splits to process" << endl;
      int isp_right = e0_frac_to_new_isp[index++].second;
      int isp_left;
      int ir = (e_invert[0]) ? 2:1;
      int il = (e_invert[0]) ? 1:2;


      while (index < index_max) {
        isp_left = e0_frac_to_new_isp[index++].second;

        new_adj_tris.push_back(NewTri());
        new_adj_tris.back().spost[0] = isp_middle;
        new_adj_tris.back().spost[ir] = isp_right;
        new_adj_tris.back().spost[il] = isp_left;
        new_adj_tris.back().znost = znost[ist];
        new_adj_tris.back().szost = szost[ist];

        isp_right = isp_left;
      }

      index = 0;
      index_max = e1_frac_to_new_isp.size();
      // cout << "e1 has " << index_max << " splits to process" << endl;
      isp_right = e1_frac_to_new_isp[index++].second;
      ir = (e_invert[1]) ? 2:1;
      il = (e_invert[1]) ? 1:2;

      while (index < index_max) {
        isp_left = e1_frac_to_new_isp[index++].second;

        new_adj_tris.push_back(NewTri());
        new_adj_tris.back().spost[0] = isp_middle;
        new_adj_tris.back().spost[ir] = isp_right;
        new_adj_tris.back().spost[il] = isp_left;
        new_adj_tris.back().znost = znost[ist];
        new_adj_tris.back().szost = szost[ist];

        isp_right = isp_left;
      }

      index = 0;
      index_max = e2_frac_to_new_isp.size();
      // cout << "e2 has " << index_max << " splits to process" << endl;
      isp_right = e2_frac_to_new_isp[index++].second;
      ir = (e_invert[2]) ? 2:1;
      il = (e_invert[2]) ? 1:2;

      while (index < index_max) {
        isp_left = e2_frac_to_new_isp[index++].second;

        new_adj_tris.push_back(NewTri());
        new_adj_tris.back().spost[0] = isp_middle;
        new_adj_tris.back().spost[ir] = isp_right;
        new_adj_tris.back().spost[il] = isp_left;
        new_adj_tris.back().znost = znost[ist];
        new_adj_tris.back().szost = szost[ist];

        isp_right = isp_left;
      }
    }
  }
}

void SimpleSurface::reDiscretizeZones(const vector<int>& zone_indices,const double delta,const double edge_factor,const double seed_factor,const int n_lloyd,const int type,const bool b_delaunay,const bool b_keep_edge_nodes,const int init_type,const bool b_power_diagram,const double growth_factor,const double growth_power) {
  bool debug = false;
  bool no_adj_feature = false;
  if (type == 1) no_adj_feature = true;

  // currently the  regions to refine are determined by flagging zones
  zone_flag.resize(zoneVec.size());
  zone_flag.setAll(0);

  for (vector<int>::const_iterator it=zone_indices.begin(); it<zone_indices.end(); ++it) {
    if (*it >= 0 && *it < int(zoneVec.size())) {  // valid zone index
      zone_flag[*it] = 1;
    }
  }

  if (b_power_diagram) {
    sp_flag.resize(nsp);
    sp_flag.setAll(0);
    FOR_IST {
      const int izn = znost[ist];
      assert((izn >= 0)&&(izn < int(zoneVec.size())));
      FOR_I3 {
        const int isp = spost[ist][i];
        if (zone_flag[izn])
          sp_flag[isp] |= 2;
        else
          sp_flag[isp] |= 1;
      }
    }

    // 1 exterior, 2 interior, 3 zone boundary...

    FOR_ISP {
      if (sp_flag[isp] == 3)
        sp_flag[isp] = 1;
      else
        sp_flag[isp] = 0;
    }
    calcGeoDistFromFlaggedNodes();

  }

  // build temporary stost structure
  int (*stost)[3];
  buildStostForFlaggedZones(&stost);

  // identify chains and protected nodes along chains
  vector<NewNode> fixedNodes;
  // also collect edges into a vector by unique node-pair (isp0<isp1, direction of edge doesn't matter)
  vector<pair<int,int> > featureEdgeVec;
  vector<pair<int,int> > adj_tris;
  vector<pair<int,int> > isp_to_nearest_fixed;  // for non-protected old nodes, store their nearest new node for zipping later
  identifyFeatureXp(fixedNodes,featureEdgeVec,adj_tris,isp_to_nearest_fixed,stost,edge_factor*delta,no_adj_feature,b_keep_edge_nodes);  // scale edge discretization smaller than internal

  // fixedNodes.isp: stores where this node sits
  // 0+: collocated with an old node, where isp is the isp of the old node
  // -1-: on a feature edge, whih is -1-indexed

  if (debug) {
    // protected node diagnostics
    FILE * fp = fopen("protected_nodes.dat","w");
    fprintf(fp,"VARIABLES=\"X\" \"Y\" \"Z\"\n");
    fprintf(fp,"ZONE T=\"protected_nodes_isps\"\n");
    int count = 0;
    for (vector<NewNode>::const_iterator it=fixedNodes.begin(); it!=fixedNodes.end(); ++it) {
      if (it->isp >= 0) {
        fprintf(fp,"%18.15e %18.15e %18.15e\n",it->xsp[0],it->xsp[1],it->xsp[2]);
        ++count;
      }
    }
    if (count == 0) fprintf(fp,"0.0 0.0 0.0\n");  // dummy

    count = 0;
    fprintf(fp,"ZONE T=\"protected_nodes_edges\"\n");
    for (vector<NewNode>::const_iterator it=fixedNodes.begin(); it!=fixedNodes.end(); ++it) {
      if (it->isp < 0) {
        fprintf(fp,"%18.15e %18.15e %18.15e\n",it->xsp[0],it->xsp[1],it->xsp[2]);
        ++count;
      }
    }
    if (count == 0) fprintf(fp,"0.0 0.0 0.0\n");  // dummy

    fclose(fp);
  }

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
  vector<NewTri> new_tris;

  // process adjacent tris
  // add replacement tris in new_tris (some new nodes potentially as well)
  zipAdjacentTris(new_adj_tris,adj_tris,edge_to_fixedNodes,isp_to_fixedNodes,fixedNodes);
  const int np_fixed = fixedNodes.size();

  vector<NewNode> new_nodes(fixedNodes);

  // in st_flag now write which retesselation group a tri belongs to (patch that will be faceted independently):
  // -1: not a participating tri
  // 0+: index of group the tri belongs to
  const int n_groups = groupFlaggedTrisByFeatures(stost);
  COUT2(" > found " << n_groups << " feature-based patches to rediscretize");

  // const double internal_delta = 1.25*min_edge_delta;
  const double internal_delta = 1.0*delta;
  COUT2(" > internal patches using effective delta of " << internal_delta);

  for (int igr=0; igr<n_groups; ++igr) {
    IntFlag stosp(np_fixed);
    stosp.setAll(-1);
    int gr_izn = -1;
    int gr_isz = -1;
    for (int ist=0; ist<nst; ++ist) {
      if (st_flag[ist] == igr) {
        // set zone,subzone for this group
        if (gr_izn == -1) {
          assert(gr_isz == -1);
          gr_izn = znost[ist];
          gr_isz = szost[ist];
        }

        FOR_I3 {
          // const int ist_nbr = stost[ist][i];
          // bool is_boundary = false;
          // if (ist_nbr < 0) is_boundary = true;
          // else if (st_flag[ist_nbr] != igr) is_boundary = true;

          // criteria was boundary (see commented above) but this ignores internal
          // feature nodes not on a boundary, thus always do a search....
          if (true) {
            const int isp0 = min(spost[ist][i],spost[ist][(i+1)%3]);
            const int isp1 = max(spost[ist][i],spost[ist][(i+1)%3]);
            pair<multimap <pair<int,int>,int>::iterator,multimap <pair<int,int>,int>::iterator> ret;
            ret = edge_to_fixedNodes.equal_range(pair<int,int>(isp0,isp1));
            for (multimap <pair<int,int>,int>::iterator it=ret.first; it!=ret.second; ++it) {
              const int ino = it->second;
              assert(ino>=0 && ino<np_fixed);
              stosp[ino] = ist;
            }

            // check for protected collocated nodes
            map<int,pair<int,bool> >::iterator mit;
            mit = isp_to_fixedNodes.find(isp0);
            if (mit!=isp_to_fixedNodes.end()) {
              if (mit->second.second) {
                const int ino = mit->second.first;
                assert(ino>=0 && ino<np_fixed);
                stosp[ino] = ist;  // ok to overwrite if hit again, just want to register to one of these tris
              }
            }
            mit = isp_to_fixedNodes.find(isp1);
            if (mit!=isp_to_fixedNodes.end()) {
              if (mit->second.second) {
                const int ino = mit->second.first;
                assert(ino>=0 && ino<np_fixed);
                stosp[ino] = ist;  // ok to overwrite if hit again, just want to register to one of these tris
              }
            }
          }
        }
      }
    }

    // want to only pass patch-local fixed nodes to 2D vor, but
    // we need a map from their local index to the global fixed index
    const int np_fixed_local = np_fixed - stosp.countNegative();
    cout << " > passing " << np_fixed_local << " fixed nodes to local node set" << endl;
    IntFlag local_to_fixed(np_fixed_local);
    int count=0;
    for (int ino=0; ino<np_fixed; ++ino) {
      if (stosp[ino] >= 0) {
        local_to_fixed[count++] = ino;
      }
    }
    assert(count == np_fixed_local);

    // initialize internal points on surface
    vector<NewNode> internal_nodes;
    const int np_internal = init2DVorPointsOnGroup(internal_nodes,internal_delta,np_fixed_local,seed_factor,igr,b_keep_edge_nodes,init_type);

    // populate local arrays
    const int np_local = np_fixed_local + np_internal;
    double (*local_xp)[3] = new double[np_local][3];
    int * local_stoxp = new int[np_local];
    double * local_delta = new double[np_local];

    for (int ino=0; ino<np_fixed_local; ++ino) {
      const int isp_fixed = local_to_fixed[ino];
      FOR_I3 local_xp[ino][i] = fixedNodes[isp_fixed].xsp[i];
      // local_stoxp[ino] = fixedNodes[isp_fixed].isp;
      // local_delta[ino] = fixedNodes[isp_fixed].data_d;
      local_stoxp[ino] = stosp[isp_fixed];  // earelier handled isp vs ist ambiguity in fixedNodes
      local_delta[ino] = 1.25*delta;  // data_d holds frac, so don't use here
    }

    for (vector<NewNode>::const_iterator it=internal_nodes.begin(); it!=internal_nodes.end(); ++it) {
      const int ino = np_fixed_local + (it-internal_nodes.begin());
      FOR_I3 local_xp[ino][i] = it->xsp[i];
      local_stoxp[ino] = it->isp;
      local_delta[ino] = it->data_d;
    }


    // xp neighbor information for eventual triangulation
    int * nboxp_i = new int[np_local+1];
    nboxp_i[0] = 0;
    vector<int> nboxp_v;
    nboxp_v.reserve(np_local*6);

    lloydIterateInternalVorPoints(local_xp,local_stoxp,local_delta,np_fixed_local,np_local,nboxp_i,nboxp_v,stost,n_lloyd,delta,0.00001,igr,b_power_diagram,growth_factor,growth_power);

    // populate global new nodes
    const int igr_no_offset = new_nodes.size();
    new_nodes.resize(igr_no_offset + np_internal);
    count = 0;
    for (int ino=np_fixed_local,index=igr_no_offset; ino < np_local; ++ino,++index) {  // don't add fixed nodes
      FOR_I3 new_nodes[index].xsp[i] = local_xp[ino][i];
      ++count;
    }
    assert(count == np_internal);

    vector<NewTri> group_tris;
    const int n_gr_tris = triangulateVorPoints(group_tris,local_xp,local_stoxp,nboxp_i,nboxp_v,np_local);

    // populate global new tris
    const int igr_tri_offset = new_tris.size();
    new_tris.resize(igr_tri_offset + n_gr_tris);
    count = 0;
    for (vector<NewTri>::iterator it=group_tris.begin(); it!=group_tris.end(); ++it) {
      // update spost for global node indices
      const int index = it-group_tris.begin();
      FOR_I3 {
        if (it->spost[i] < np_fixed_local) {
          // this references a fixed node; determine the global index
          new_tris[igr_tri_offset+index].spost[i] = local_to_fixed[it->spost[i]];
        }
        else {
          // new node, use it's global index
          new_tris[igr_tri_offset+index].spost[i] = igr_no_offset + (it->spost[i]-np_fixed_local);
        }
      }
      // for now each patch into zone/subzone of the group (should be subzone patches)
      // this means zn/sz structures shouldn't need to be updated...
      new_tris[igr_tri_offset+index].znost = gr_izn;
      new_tris[igr_tri_offset+index].szost = gr_isz;
      ++count;
    }
    assert(count == n_gr_tris);

    if (debug) {
      // diagnostic output
      char filename[128];
      sprintf(filename,"group_tris_and_nodes.%03d.dat",igr);
      FILE * fp = fopen(filename,"w");
      fprintf(fp,"VARIABLES=\"X\" \"Y\" \"Z\"\n");
      fprintf(fp,"ZONE T=\"tris_group_%d\"\n",igr);
      int count = st_flag.countEqualTo(igr);
      fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",3*count,count);
      for (int ist=0; ist<nst; ++ist) {
        if (st_flag[ist] == igr) {
          FOR_I3 fprintf(fp,"%18.15e %18.15e %18.15e\n",xsp[spost[ist][i]][0],xsp[spost[ist][i]][1],xsp[spost[ist][i]][2]);
        }
      }
      count = 0;
      for (int ist=0; ist<nst; ++ist) {
        if (st_flag[ist] == igr) {
          const int offset = count*3 + 1;
          fprintf(fp,"%d %d %d\n",offset,offset+1,offset+2);
          ++count;
        }
      }

      fprintf(fp,"ZONE T=\"protect_nodes_group_%d\"\n",igr);
      for (int ino=0; ino<np_fixed_local; ++ino) {
        fprintf(fp,"%18.15e %18.15e %18.15e\n",local_xp[ino][0],local_xp[ino][1],local_xp[ino][2]);
      }
      fprintf(fp,"ZONE T=\"internal_nodes_group_%d\"\n",igr);
      for (int ino=np_fixed_local; ino<np_local; ++ino) {
        fprintf(fp,"%18.15e %18.15e %18.15e\n",local_xp[ino][0],local_xp[ino][1],local_xp[ino][2]);
      }

      count = group_tris.size();
      if (count) {
        fprintf(fp,"ZONE T=\"newtris_group_%d\"\n",igr);
        fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",3*count,count);
        for (vector<NewTri>::iterator it=group_tris.begin(); it!=group_tris.end(); ++it) {
          FOR_I3 fprintf(fp,"%18.15e %18.15e %18.15e\n",local_xp[it->spost[i]][0],local_xp[it->spost[i]][1],local_xp[it->spost[i]][2]);
        }
        count = 0;
        for (int ist=0, lim=group_tris.size(); ist<lim; ++ist) {
          const int offset = count*3 + 1;
          fprintf(fp,"%d %d %d\n",offset,offset+1,offset+2);
          ++count;
        }
      }
      else {
        fprintf(fp,"ZONE T=\"newtris_group_%d\"\n",igr);
        fprintf(fp,"0.0 0.0 0.0\n");
      }
      fclose(fp);
    }

    DELETE(local_xp);
    DELETE(local_stoxp);
    DELETE(local_delta);
    DELETE(nboxp_i);
    nboxp_v.clear();
  }

  if (debug) {
    // diagnostic for new nodes, tris
    char filename[128];
    sprintf(filename,"new_tris.dat");
    FILE * fp = fopen(filename,"w");
    fprintf(fp,"VARIABLES=\"X\" \"Y\" \"Z\"\n");
    fprintf(fp,"ZONE T=\"retessellated\"\n");
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",int(new_nodes.size()),int(new_tris.size()));
    for (vector<NewNode>::iterator it=new_nodes.begin(); it!=new_nodes.end(); ++it) {
      fprintf(fp,"%18.15e %18.15e %18.15e\n",it->xsp[0],it->xsp[1],it->xsp[2]);
    }
    for (vector<NewTri>::iterator it=new_tris.begin(); it!=new_tris.end(); ++it) {
      fprintf(fp,"%d %d %d\n",it->spost[0]+1,it->spost[1]+1,it->spost[2]+1);  // tecplot is 1-indexed
    }
    fclose(fp);
  }

  DELETE(stost);

  // build a local teost to do Delaunay edge flipping with
  if (b_delaunay) {
    double (*local_xsp)[3] = new double[new_nodes.size()][3];
    for (int isp=0, lim=new_nodes.size(); isp<lim; ++isp) {
      FOR_I3 local_xsp[isp][i] = new_nodes[isp].xsp[i];
    }
    int (*local_spost)[3] = new int[new_tris.size()][3];
    for (int ist=0, lim=new_tris.size(); ist<lim; ++ist) {
      FOR_I3 local_spost[ist][i] = new_tris[ist].spost[i];
    }
    Teost local_teost;
    bool b_teost = local_teost.build(new_nodes.size(),int(new_tris.size()),local_spost);
    COUT2(" > retiangulation teost build status: " << (b_teost ? "success":"failed"));

    makeTrisDelaunay(local_teost,local_spost,int(new_tris.size()),local_xsp);

    // now update
    for (int ist=0, lim=new_tris.size(); ist<lim; ++ist) {
      FOR_I3 new_tris[ist].spost[i] = local_spost[ist][i];
    }
    DELETE(local_xsp);
    DELETE(local_spost);
  }

  // process adjacent tris and make them conformal with new tessellation at feature boundaries
  // remove adjacent tris from old surface by faking them as a separate group in st_flag
  for (vector<pair<int,int> >::iterator it=adj_tris.begin(); it!=adj_tris.end(); ++it) {
    st_flag[it->first] = n_groups;
  }

  // cull old surface and replace with new
  updateSurfaceWithRetessellation(new_nodes,np_fixed,new_adj_tris,new_tris);
}

void SimpleSurface::splitZonesAfterCopyAll(const int n) {

  cout << "in splitZonesAfterCopyAll: n=" << n << endl;
  int nzn = zoneVec.size();

  // blow up the zoneVec...
  for (int ii = 1; ii < n; ++ii) {
    for (int izn = 0; izn < nzn; ++izn) {
      stringstream ss;
      ss << zoneVec[izn].getName() << "_" << ii;
      int izn_check = addNewZone(ss.str());
      assert(izn_check == nzn*ii+izn);
    }
  }
  // and add "_0" to the first ones...
  for (int izn = 0; izn < nzn; ++izn) {
    stringstream ss;
    ss << zoneVec[izn].getName() << "_0";
    zoneVec[izn].setName(ss.str());
  }

  for (int izn = 0; izn < nzn; ++izn) {
    cout << "working on zone: " << zoneVec[izn].getName() << endl;
    bool b_on = false;
    int ist_on;
    int count = 0;
    for (int ist = 0; ist < nst; ++ist) {
      const int izn_ = znost[ist];
      if (izn_ == izn) {
	++count;
	if (!b_on) {
	  cout << " > starts at ist: " << ist << endl;
	  ist_on = ist;
	  b_on = true;
	}
      }
      else if (b_on) {
	cout << " > stops at ist: " << ist-1 << " ntri: " << ist-ist_on << endl;
	b_on = false;
      }
    }
    if (b_on) {
      cout << " > stops at ist: " << nst-1 << " ntri: " << nst-ist_on << endl;
      b_on = false;
    }
    assert(count%n == 0);
    int count_flip = count/n;
    int ii = 0;
    count = 0;
    for (int ist = 0; ist < nst; ++ist) {
      const int izn_ = znost[ist];
      if (izn_ == izn) {
	znost[ist] = izn + nzn*ii;
	++count;
	if (count == count_flip) {
	  count = 0;
	  ++ii;
	}
      }
    }

  }

  setSubzonesToZones();

  double d2_min = HUGE_VAL;
  for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
      const int isp0 = spost[ist][i];
      const int isp1 = spost[ist][(i+1)%3];
      double d2 = DIST2(xsp[isp0],xsp[isp1]);
      d2_min = min(d2_min,d2);
    }
  }
  cout << " XXXXXXXXXXXXXXXXX d_min (for MERGE_COLLOCATED_NODES): " << sqrt(d2_min) << endl;

}
