#include "SimpleSurface.hpp"

void SimpleSurface::imprintPlane(const double x_plane[3],const double _n_plane[3],const double node_tol,const bool splitSubZones) {

  // TODO: figure out how to condense this code into and use/modify the
  // general imprint(sp_dist,... routine below, like imprintCyl below

  const int nsz0 = nsz;
  ensureTeost();

  COUT1(" > imprinting edges from plane:");
  COUT1("    > point: " << COUT_VEC(x_plane));
  // unit normal important?...
  double n_plane[3] = {_n_plane[0],_n_plane[1],_n_plane[2]};
  NORMALIZE(n_plane);
  COUT1("    > normal: " << COUT_VEC(n_plane) << " (normalized)");

  sp_flag.resize(nsp);
  sp_flag.setAll(-2);  // -1,0,1 will be used for sign of signed distance

  double * sp_dist = new double[nsp];

  // compute signed dist to plane for each point on surface
  for (int isp=0; isp<nsp; ++isp) {
    sp_dist[isp] =
      (xsp[isp][0]-x_plane[0])*n_plane[0] +
      (xsp[isp][1]-x_plane[1])*n_plane[1] +
      (xsp[isp][2]-x_plane[2])*n_plane[2];

    // node_tol determines distance to plane within which we say we are on the plane
    // by default this is zero but the user can specify a specific value
    //TODO: we may need to add in checks that ensurethis is an appropriate length scale...
    if (sp_dist[isp] < (0.0-node_tol)) sp_flag[isp] = -1;
    else if (sp_dist[isp] > (0.0+node_tol)) sp_flag[isp] = 1;
    else sp_flag[isp] = 0;

    sp_dist[isp] = fabs(sp_dist[isp]);  // store abs value so weights later easier. Sign is stored in sp_flag
  }

  // loop through tris to see who gets split; also perform split
  stack<int> next_tri_stack;
  typedef map<pair<int,int>,int > SplitEdgeMap;  // ist,ege : new_isp
  SplitEdgeMap split_edge_node;
  map<int,pair<int,int> > forced_split_edge_node;
  IntFlag n_split(3);
  IntFlag isp_sign(3);
  bool b_multi = false;  // did we hit multiedges

  vector<NewNode> newNodesVec;
  vector<NewTri> newTrisVec;
  map<pair<string,int>,pair<int,int> > newZonesMap;
  map<pair<string,int>,pair<int,int> >::iterator zm_it;

  st_flag.resize(nst);
  st_flag.setAll(0);  // indicator of visited yet or not

  IntFlag szost_global(szost);  // globally indexed subzones
  globalSzostToLocal();

  zone_flag.resize(zoneVec.size());
  zone_flag.setAll(0);
  // store local max sz index in zone_flag
  for (int izn=0,nzn=zoneVec.size(); izn<nzn; ++izn) {
    zone_flag[izn] = szozn_i[izn+1] - szozn_i[izn];
  }

  for (int ist=0; ist<nst; ++ist) {
    if (st_flag[ist] == 0) {
      next_tri_stack.push(ist);

      while (!next_tri_stack.empty()) {
        const int this_ist = next_tri_stack.top(); next_tri_stack.pop();

        // don't process if subzone isn't flagged
        if (sz_flag[szost_global[this_ist]] == 0) st_flag[this_ist] = 1;

        if (st_flag[this_ist]) continue;  // don't process
        else st_flag[this_ist] = 1;  // set as visited

        n_split.setAll(-1);
        isp_sign.setAll(-1);
        int isp0,isp1 = -1;
        FOR_I3 {
          isp0 = spost[this_ist][i];
          isp_sign[i] = sp_flag[isp0];
          isp1 = spost[this_ist][(i+1)%3];
          if ( (sp_flag[isp0] == -1 && sp_flag[isp1] == 1) || (sp_flag[isp0] == 1 && sp_flag[isp1] == -1) ) {
            n_split[i] = 0;
          }
          // else if (sp_flag[isp0] == 0 && sp_flag[isp1] == 0) {
          //   // edge is in plane
          // }
        }

        if (splitSubZones && (isp_sign.countNegative() == 0) ) {
          // to one side (or touching) plane, so move to new zone
          // add map from current zone name to new zone,sz index
          const string this_zonename = zoneVec[znost[this_ist]].getName()+"_split";

          int zone_exists = -1;
          for (int izn=0,nzn=zoneVec.size(); izn<nzn; ++izn) {
            if (this_zonename == zoneVec[izn].getName() ) {
              zone_exists = izn;
              break;
            }
          }

          zm_it = newZonesMap.find(pair<string,int> (this_zonename,szost[this_ist]));
          if (zm_it != newZonesMap.end()) {
            // already in map, so use this value
            znost[this_ist] = zm_it->second.first;
            szost[this_ist] = zm_it->second.second;
          }
          else {
            // hasn't been added yet, so add to map and zone list
            int nzn;
            if (zone_exists == -1) {
              nzn = zoneVec.size();
              zoneVec.push_back(SurfaceZone(this_zonename));
            }
            else {
              nzn = zone_exists;
            }

            newZonesMap.insert( pair<pair<string,int>,pair<int,int> > (pair<string,int> (this_zonename,szost[this_ist]),pair<int,int>(nzn,nsz)));
            znost[this_ist] = nzn;  // change parent zone
            szost[this_ist] = nsz;
            ++nsz;
          }
        }
        else if (isp_sign.countNegative() == 0) {
          const string this_zonename = zoneVec[znost[this_ist]].getName();
          zm_it = newZonesMap.find(pair<string,int> (this_zonename,szost[this_ist]));
          if (zm_it != newZonesMap.end()) {
            // already in map, so use this value
            znost[this_ist] = zm_it->second.first;
            szost[this_ist] = zm_it->second.second;
          }
          else {
            // hasn't been added yet, so add to map and zone list
            newZonesMap.insert( pair<pair<string,int>,pair<int,int> > (pair<string,int> (this_zonename,szost[this_ist]),pair<int,int>(znost[this_ist],nsz)));
            szost[this_ist] = nsz;
            ++nsz;
          }
        }

        const int n_splits = n_split.countZero();
        // only process further if edge(s) are split
        if (n_splits) {
          // for each split edge determine/insert the new splitting node
          FOR_I3 {
            if (n_split[i] == 0) {
              // this edge is split, so process

              // store new node's isp in n_split
              // first search in map to see if this edge already has a node defined
              SplitEdgeMap::iterator iter=split_edge_node.find(pair<int,int>(this_ist,i));
              if (iter != split_edge_node.end()) {
                // found a node to use already
                n_split[i] = iter->second;
                split_edge_node.erase(iter);
              }
              else {
                // need to introduce the new node
                n_split[i] = nsp + newNodesVec.size();  // new global index

                //TODO account for multiedges too
                int ist_nbr,i_nbr,orient_nbr;
                const int nbr_type = getTriNbrDataFull(ist_nbr,i_nbr,orient_nbr,this_ist,i);
                if (nbr_type == 1 || nbr_type == 2) {
                  if (sz_flag[szost_global[ist_nbr]] == 0) {
                    forced_split_edge_node.insert(pair<int,pair<int,int> > (ist_nbr,pair<int,int>(i_nbr,n_split[i])));
                  }
                  else {
                    split_edge_node.insert(pair<pair<int,int>,int > (pair<int,int>(ist_nbr,i_nbr),n_split[i]));
                    next_tri_stack.push(ist_nbr);
                  }
                }
                else if (nbr_type == 3) b_multi = true;

                // compute and push the new node
                isp0 = spost[this_ist][i];
                isp1 = spost[this_ist][(i+1)%3];
                const double dx[3] = DIFF(xsp[isp1],xsp[isp0]);
                const double frac = sp_dist[isp0]/(sp_dist[isp0]+sp_dist[isp1]);
                const double new_x = xsp[isp0][0] + frac*dx[0];
                const double new_y = xsp[isp0][1] + frac*dx[1];
                const double new_z = xsp[isp0][2] + frac*dx[2];
                newNodesVec.push_back(NewNode(new_x,new_y,new_z));
              }
            }
          }

          int new_zone = znost[this_ist];
          int new_sz = szost[this_ist];

          // if tris to one side of plane will be given to a new zone, we need to keep track of the new zone(s)
          if (splitSubZones) {
            // add map from current zone name to new zone index
            const string this_zonename = zoneVec[new_zone].getName()+"_split";

            int zone_exists = -1;
            for (int izn=0,nzn=zoneVec.size(); izn<nzn; ++izn) {
              if (this_zonename == zoneVec[izn].getName() ) {
                zone_exists = izn;
                break;
              }
            }

            zm_it = newZonesMap.find(pair<string,int> (this_zonename,szost[this_ist]));
            if (zm_it != newZonesMap.end()) {
              // already in map, so use this value
              new_zone = zm_it->second.first;
              new_sz = zm_it->second.second;
            }
            else {
              // hasn't been added yet, so add to map and zone list
              // only add new zone if this zone wasn't found yet

              if (zone_exists == -1) {
                new_zone = zoneVec.size();
                zoneVec.push_back(SurfaceZone(this_zonename));
              }
              else {
                new_zone = zone_exists;
              }

              newZonesMap.insert( pair<pair<string,int>,pair<int,int> > (pair<string,int> (this_zonename,szost[this_ist]),pair<int,int> (new_zone,nsz) ));
              new_sz = nsz;
              ++nsz;
            }
          }
          else {
            const string this_zonename = zoneVec[new_zone].getName();
            zm_it = newZonesMap.find(pair<string,int> (this_zonename,szost[this_ist]));
            if (zm_it != newZonesMap.end()) {
              // already in map, so use this value
              new_zone = zm_it->second.first;
              new_sz = zm_it->second.second;
            }
            else {
              // hasn't been added yet, so add to map and zone list
              newZonesMap.insert( pair<pair<string,int>,pair<int,int> > (pair<string,int> (this_zonename,szost[this_ist]),pair<int,int> (new_zone,nsz) ));
              new_sz = nsz;
              ++nsz;
            }
          }

          // insert the new tris and update the old tris in spost
          if (n_splits == 1) {
            // split tri down the middle
            FOR_I3 {
              if (n_split[i] > -1) {

                if (sp_flag[spost[this_ist][i]] > 0) {
                  newTrisVec.push_back(NewTri(spost[this_ist][i],n_split[i],spost[this_ist][(i+2)%3],new_zone,new_sz));
                  NewTri old_tri(n_split[i],spost[this_ist][(i+1)%3],spost[this_ist][(i+2)%3],0,0);  // temp placeholder so we don't overwrite
                  FOR_J3 spost[this_ist][j] = old_tri.spost[j];

                }
                else {
                  newTrisVec.push_back(NewTri(spost[this_ist][i],n_split[i],spost[this_ist][(i+2)%3],znost[this_ist],szost[this_ist]));
                  NewTri old_tri(n_split[i],spost[this_ist][(i+1)%3],spost[this_ist][(i+2)%3],0,0);  // temp placeholder so we don't overwrite
                  FOR_J3 spost[this_ist][j] = old_tri.spost[j];
                  znost[this_ist] = new_zone;
                  szost[this_ist] = new_sz;
                }

                break;  // only one, so done
              }
            }
          }
          else if (n_splits == 2) {
            // split tri between nodes
            FOR_I3 {
              if (n_split[i] == -1) {
                // find the edge that isn't split
                if (sp_flag[spost[this_ist][(i+2)%3]] > 0) {
                  newTrisVec.push_back(NewTri(spost[this_ist][i],spost[this_ist][(i+1)%3],n_split[(i+1)%3],znost[this_ist],szost[this_ist]));
                  newTrisVec.push_back(NewTri(spost[this_ist][i],n_split[(i+1)%3],n_split[(i+2)%3],znost[this_ist],szost[this_ist]));

                  NewTri old_tri(n_split[(i+1)%3],spost[this_ist][(i+2)%3],n_split[(i+2)%3],0,0);  // temp placeholder so we don't overwrite
                  FOR_J3 spost[this_ist][j] = old_tri.spost[j];
                  znost[this_ist] = new_zone;
                  szost[this_ist] = new_sz;
                }
                else {
                  newTrisVec.push_back(NewTri(spost[this_ist][i],spost[this_ist][(i+1)%3],n_split[(i+1)%3],new_zone,new_sz));
                  newTrisVec.push_back(NewTri(spost[this_ist][i],n_split[(i+1)%3],n_split[(i+2)%3],new_zone,new_sz));

                  NewTri old_tri(n_split[(i+1)%3],spost[this_ist][(i+2)%3],n_split[(i+2)%3],0,0);  // temp placeholder so we don't overwrite
                  FOR_J3 spost[this_ist][j] = old_tri.spost[j];
                }

                break;  // only one, so done
              }
            }
          }
          assert(n_splits != 3);
        }
      }

      assert(split_edge_node.empty());
    }
  }
  DELETE(sp_dist);

  // process tris forced to be split in neighboring zones
  if (!forced_split_edge_node.empty()) {
    for (map<int,pair<int,int> >::iterator it=forced_split_edge_node.begin(); it!=forced_split_edge_node.end(); ) {
      const int this_ist = it->first;
      IntFlag n_split(3);
      n_split.setAll(-1);
      while (it->first == this_ist) {
        n_split[it->second.first] = it->second.second;  // set edge splitting node as value
        ++it;
      }

      // insert the new tris and update the old tris in spost
      const int n_splits = 3-n_split.countNegative();
      if (n_splits == 1) {
        // split tri down the middle
        FOR_I3 {
          if (n_split[i] > -1) {
            newTrisVec.push_back(NewTri(spost[this_ist][i],n_split[i],spost[this_ist][(i+2)%3],znost[this_ist],szost[this_ist]));
            NewTri old_tri(n_split[i],spost[this_ist][(i+1)%3],spost[this_ist][(i+2)%3],0,0);  // temp placeholder so we don't overwrite
            FOR_J3 spost[this_ist][j] = old_tri.spost[j];

            break;  // only one, so done
          }
        }
      }
      else if (n_splits == 2) {
        // split tri between nodes
        FOR_I3 {
          if (n_split[i] == -1) {
            // find the edge that isn't split
            newTrisVec.push_back(NewTri(spost[this_ist][i],spost[this_ist][(i+1)%3],n_split[(i+1)%3],znost[this_ist],szost[this_ist]));
            newTrisVec.push_back(NewTri(spost[this_ist][i],n_split[(i+1)%3],n_split[(i+2)%3],znost[this_ist],szost[this_ist]));

            NewTri old_tri(n_split[(i+1)%3],spost[this_ist][(i+2)%3],n_split[(i+2)%3],0,0);  // temp placeholder so we don't overwrite
            FOR_J3 spost[this_ist][j] = old_tri.spost[j];

            break;  // done
          }
        }
      }
      assert(n_splits != 3);
    }

  }

  if (newTrisVec.size()) {
    if (b_multi) CWARN("multi-edges were intersected but double counted; you may need to collapse collocated nodes");
    COUT2(" > new nodes introduced: " << newNodesVec.size());
    COUT2(" > new tris introduced: " << newTrisVec.size());
    if (splitSubZones) COUT2(" > new zones introduced: " << newZonesMap.size() << " (with suffix '_split')");

    // append new nodes
    double (*xsp0)[3] = xsp;
    xsp = new double[nsp + newNodesVec.size()][3];
    for (int isp = 0; isp < nsp; ++isp) {
      FOR_I3 xsp[isp][i] = xsp0[isp][i];
    }
    DELETE(xsp0);
    for (int isp=0, limit=newNodesVec.size(); isp<limit; ++isp) {
      FOR_I3 xsp[nsp+isp][i] = newNodesVec[isp].xsp[i];
    }
    nsp += newNodesVec.size();

    // append new tris
    int (*spost0)[3] = spost;
    spost = new int[nst + newTrisVec.size()][3];
    int * znost0 = znost;
    znost = new int[nst + newTrisVec.size()];
    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 spost[ist][i] = spost0[ist][i];
      znost[ist] = znost0[ist];
    }
    DELETE(spost0);
    DELETE(znost0);
    szost.resize(nst+newTrisVec.size());
    for (int ist=0, limit=newTrisVec.size(); ist<limit; ++ist) {
      FOR_I3 spost[nst+ist][i] = newTrisVec[ist].spost[i];
      znost[nst+ist] = newTrisVec[ist].znost;
      szost[nst+ist] = newTrisVec[ist].szost;
    }
    nst += newTrisVec.size();

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
  }
  else if (splitSubZones || (nsz != nsz0)) {
    // new zones have local-indexed subzones from source zones
    // to rebuild properly pruned szozn_i and szost we need to run these
    buildSzoznFromLocalSzost();
    localSzostToGlobal();
    pruneEmptyZonesAndSubzones();

    clearSubzoneData();
    clearZoneData();
    clearStoszSzosz();
    clearDynamicEdgeGroups();
  }
  else {
    localSzostToGlobal();  // convert back
    COUT1(" > no intersections detected; surface wasn't modified");
  }
}


void SimpleSurface::imprintCyl(const double x_cyl[3],const double n_cyl[3],const double r_cyl,const bool splitZones) {

  // coming into this routine, sz_flag contains 1 in the relevant subzones...

  // flag the sp's involved...

  sp_flag.resize(nsp);
  sp_flag.setAll(0);

  for (int ist = 0; ist < nst; ++ist) {
    if (sz_flag[szost[ist]]) {
      FOR_I3 sp_flag[spost[ist][i]] = 1;
    }
  }

  double * sp_dist = new double[nsp];

  // compute signed dist to cyl for each point on surface
  const double n_cyl_mag2 = DOT_PRODUCT(n_cyl,n_cyl);
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] == 1) {
      const double dx[3] = { xsp[isp][0]-x_cyl[0], xsp[isp][1]-x_cyl[1], xsp[isp][2]-x_cyl[2] };
      const double dx_mag2 = DOT_PRODUCT(dx,dx);
      const double dp = DOT_PRODUCT(dx,n_cyl);
      const double r2 = dx_mag2 - dp*dp/n_cyl_mag2;
      if (r2 > 0.0) {
        sp_dist[isp] = sqrt(r2) - r_cyl;
      }
      else {
        // this can happen because of roundoff in certain rare cases...
        sp_dist[isp] = - r_cyl;
      }
    }
  }

  // now we have the signed distance in every sp, so split...

  imprint(sp_dist,splitZones);

  delete[] sp_dist;

}

void SimpleSurface::imprint(const double * const sp_dist,const bool splitZones) {

  map<const pair<int,int>,int> edMap;
  vector<double> xspVec; // new points associated with split edges
  int nsp_new = nsp;
  int edost[3]; // temporary edge index storage
  vector<NewTri> newTriVec; // new tris spost and znost

  buildSzoznFromGlobalSzost();
  IntFlag nsz0_zone(zoneVec.size());  // store count of subzones per zone
  for (int izn=0, nzn=zoneVec.size(); izn<nzn; ++izn) {
    nsz0_zone[izn] = szozn_i[izn+1]-szozn_i[izn];  // index for new subzone in this zone
    // cout << "nsz for zone " << izn << ": " << nsz0_zone[izn] << endl;
  }
  globalSzostToLocal();  // make szost a local index

  for (int ist = 0; ist < nst; ++ist) {
    const int global_isz = szozn_i[znost[ist]] + szost[ist];
    if (sz_flag[global_isz]) {

      // loop through the edges...
      int tri_bits = 0;
      if (sp_dist[spost[ist][0]] > 0.0) tri_bits |= 1;
      else if (sp_dist[spost[ist][0]] == 0.0) tri_bits |= 2;
      if (sp_dist[spost[ist][1]] > 0.0) tri_bits |= 4;
      else if (sp_dist[spost[ist][1]] == 0.0) tri_bits |= 8;
      if (sp_dist[spost[ist][2]] > 0.0) tri_bits |= 16;
      else if (sp_dist[spost[ist][2]] == 0.0) tri_bits |= 32;

      // set edost. This int[3] stores the new node index along each edge,
      // if present. This mirrors the vertlist[3] structure in the marching
      // cubes implementation...
      // assumes only a single new node is inserted...
      FOR_I3 {
        if (sp_dist[spost[ist][i]]*sp_dist[spost[ist][(i+1)%3]] < 0.0) {
          map<const pair<int,int>,int>::iterator iter = edMap.find(pair<int,int>(spost[ist][(i+1)%3],spost[ist][i]));
          if (iter == edMap.end()) {
            edMap[pair<int,int>(spost[ist][i],spost[ist][(i+1)%3])] = edost[i] = nsp_new++;
            const double w1 = sp_dist[spost[ist][i]]/(sp_dist[spost[ist][i]]-sp_dist[spost[ist][(i+1)%3]]);
            FOR_J3 xspVec.push_back(w1*xsp[spost[ist][(i+1)%3]][j]+(1.0-w1)*xsp[spost[ist][i]][j]);
          }
          else {
            edost[i] = iter->second;
            edMap.erase(iter);
          }
        }
        else {
          // not strictly necessary, because it should not be called, but
          // useful for checking...
          edost[i] = -1;
        }
      }

      switch (tri_bits) {
      case 2:
      case 8:
      case 32:
      case 40:
      case 10:
      case 34:
      case 0:
        // entire tri is negative
        szost[ist] += nsz0_zone[znost[ist]];
        break;
      case 21:
      case 25:
      case 22:
      case 37:
      case 41:
      case 26:
      case 38:
        // entire tri is positive
        break;
      case 18:
        // new tri...
        newTriVec.push_back(NewTri(spost[ist][0],spost[ist][1],edost[1],znost[ist],szost[ist]+nsz0_zone[znost[ist]]));
        // and reset the original tri...
        spost[ist][1] = edost[1];
        break;
      case 6:
        // new tri...
        newTriVec.push_back(NewTri(spost[ist][0],spost[ist][1],edost[1],znost[ist],szost[ist]));
        // and reset the original tri...
        spost[ist][1] = edost[1];
        szost[ist] += nsz0_zone[znost[ist]];
        break;
      case 24:
        // new tri...
        newTriVec.push_back(NewTri(spost[ist][0],spost[ist][1],edost[2],znost[ist],szost[ist]+nsz0_zone[znost[ist]]));
        // and reset the original tri...
        spost[ist][0] = edost[2];
        break;
      case 9:
        // new tri...
        newTriVec.push_back(NewTri(spost[ist][0],spost[ist][1],edost[2],znost[ist],szost[ist]));
        // and reset the original tri...
        spost[ist][0] = edost[2];
        szost[ist] += nsz0_zone[znost[ist]];
        break;
      case 36:
        // new tri...
        newTriVec.push_back(NewTri(spost[ist][0],edost[0],spost[ist][2],znost[ist],szost[ist]+nsz0_zone[znost[ist]]));
        // and reset the original tri...
        spost[ist][0] = edost[0];
        break;
      case 33:
        // new tri...
        newTriVec.push_back(NewTri(spost[ist][0],edost[0],spost[ist][2],znost[ist],szost[ist]));
        // and reset the original tri...
        spost[ist][0] = edost[0];
        szost[ist] += nsz0_zone[znost[ist]];
        break;
      case 1:
        // new tri 0...
        newTriVec.push_back(NewTri(spost[ist][0],edost[0],edost[2],znost[ist],szost[ist]));
        // new tri 1...
        newTriVec.push_back(NewTri(edost[0],spost[ist][2],edost[2],znost[ist],szost[ist]+nsz0_zone[znost[ist]]));
        // and reset the original tri...
        spost[ist][0] = edost[0];
        szost[ist] += nsz0_zone[znost[ist]];
        break;
      case 4:
        // new tri 0...
        newTriVec.push_back(NewTri(spost[ist][0],edost[0],edost[1],znost[ist],szost[ist]+nsz0_zone[znost[ist]]));
        // new tri 1...
        newTriVec.push_back(NewTri(edost[0],spost[ist][1],edost[1],znost[ist],szost[ist]));
        // and reset the original tri...
        spost[ist][1] = edost[1];
        szost[ist] += nsz0_zone[znost[ist]];
        break;
      case 5:
        // new tri 0...
        newTriVec.push_back(NewTri(edost[2],spost[ist][1],edost[1],znost[ist],szost[ist]));
        // new tri 1...
        newTriVec.push_back(NewTri(edost[2],edost[1],spost[ist][2],znost[ist],szost[ist]+nsz0_zone[znost[ist]]));
        // and reset the original tri...
        spost[ist][2] = edost[2];
        break;
      case 16:
        // new tri 0...
        newTriVec.push_back(NewTri(edost[2],spost[ist][1],edost[1],znost[ist],szost[ist]+nsz0_zone[znost[ist]]));
        // new tri 1...
        newTriVec.push_back(NewTri(edost[2],edost[1],spost[ist][2],znost[ist],szost[ist]));
        // and reset the original tri...
        spost[ist][2] = edost[2];
        szost[ist] += nsz0_zone[znost[ist]];
        break;
      case 17:
        // new tri 0...
        newTriVec.push_back(NewTri(spost[ist][0],edost[0],edost[1],znost[ist],szost[ist]));
        // new tri 1...
        newTriVec.push_back(NewTri(edost[0],spost[ist][1],edost[1],znost[ist],szost[ist]+nsz0_zone[znost[ist]]));
        // and reset the original tri...
        spost[ist][1] = edost[1];
        break;
      case 20:
        // new tri 0...
        newTriVec.push_back(NewTri(spost[ist][0],edost[0],edost[2],znost[ist],szost[ist]+nsz0_zone[znost[ist]]));
        // new tri 1...
        newTriVec.push_back(NewTri(edost[0],spost[ist][2],edost[2],znost[ist],szost[ist]));
        // and reset the original tri...
        spost[ist][0] = edost[0];
        break;
      default:
        cout << "default tri_bits: " << tri_bits << endl;
        cout << "sp_dist[spost[ist][0]]: " << sp_dist[spost[ist][0]] << " sp_dist[spost[ist][1]]: " << sp_dist[spost[ist][1]] << " sp_dist[spost[ist][2]]: " << sp_dist[spost[ist][2]] << endl;
        assert(0);
      }
    }
  }

  // add the new points...
  int nsp_orig = nsp; nsp = nsp_new;
  growNspData(nsp,nsp_orig);

  for (int isp = nsp_orig; isp < nsp; ++isp) {
    FOR_I3 xsp[isp][i] = xspVec[(isp-nsp_orig)*3+i];
  }
  xspVec.clear();

  // and add the new tris...
  const int nst_inc = newTriVec.size();
  int nst_orig = nst; nst += nst_inc;
  growNstData(nst,nst_orig);

  for (int ist = 0; ist < nst_inc; ++ist) {
    FOR_I3 spost[nst_orig+ist][i] = newTriVec[ist].spost[i];
    FOR_I3 assert(spost[nst_orig+ist][i] >= 0);
    znost[nst_orig+ist] = newTriVec[ist].znost;
    szost[nst_orig+ist] = newTriVec[ist].szost;
  }

  if (splitZones) {
    // move tris in new subzones to new zones, create those zones if not created yet
    // should only be required for new zones
    IntFlag zone_split_index(zoneVec.size());
    zone_split_index.setAll(-1);

    FOR_IST {
      if (szost[ist] == nsz0_zone[znost[ist]]) {
        // this tri is in a new subzone
        if (zone_split_index[znost[ist]] == -1) {
          // new zone hasn't been created/assigned yet
          zone_split_index[znost[ist]] = zoneVec.size();
          assert(zone_split_index[znost[ist]] >= 0);
          zoneVec.push_back(SurfaceZone(zoneVec[znost[ist]].getName()+"_split"));
        }
        szost[ist] = 0;  // only subzone in this zone
        znost[ist] = zone_split_index[znost[ist]];
      }
    }
  }
  buildSzoznFromLocalSzost();
  localSzostToGlobal();

  // lazy; just clearing these for now
  selectedSubzoneVec.clear();
  hiddenSubzoneVec.clear();

  pruneEmptyZonesAndSubzones();

  clearSubzoneData();
  clearZoneData();
  clearStoszSzosz();
  clearDynamicEdgeGroups();
}
