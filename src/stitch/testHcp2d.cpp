#include "CTI.hpp"
using namespace CTI;

#include "SurfaceShm.hpp"
#include "CuttableVoronoiData.hpp"
#include <stack>
#include "GeomUtils.hpp"

class TripleInt {
public:
  int first,second,third;
  TripleInt(const int first,const int second,const int third) {
    this->first = first;
    this->second = second;
    this->third = third;
  }
};


struct NewXp {
public:
  double xyz[3];
  double delta;
  int isp;
  int st;

  NewXp() {
    FOR_I3 xyz[i] = -1.0;
    delta = -1.0;
    st = -1;
    isp = -1;
  }

  NewXp(const NewXp& old) {
    FOR_I3 xyz[i] = old.xyz[i];
    delta = old.delta;
    st = old.st;
    isp = old.isp;
  }
};

void writeTriToFp(FILE * fp,SurfaceShm * surface,const int ist) {
  FOR_I3 fprintf(fp,"%18.15e %18.15e %18.15e\n",surface->xsp[surface->spost[ist][i]][0],surface->xsp[surface->spost[ist][i]][1],surface->xsp[surface->spost[ist][i]][2]);
  int i = 0;
  fprintf(fp,"%18.15e %18.15e %18.15e\n",surface->xsp[surface->spost[ist][i]][0],surface->xsp[surface->spost[ist][i]][1],surface->xsp[surface->spost[ist][i]][2]);
}

void writeXpToFp(const int step,const double (*xp)[3],const int np,const int* flag = NULL,const double* data = NULL) {
  char filename[128];
  sprintf(filename,"xp.%06d.dat",step);
  FILE * fp = fopen(filename,"w");
  if ((data==NULL) && (flag==NULL)) {
    fprintf(fp,"VARIABLES=\"X\" \"Y\" \"Z\"\n");
    for (int ip = 0; ip < np; ++ip) {
      fprintf(fp,"%18.15e %18.15e %18.15e\n",xp[ip][0],xp[ip][1],xp[ip][2]);
    }
  }
  else if (data==NULL) {
    fprintf(fp,"VARIABLES=\"X\" \"Y\" \"Z\" \"flag\"\n");
    for (int ip = 0; ip < np; ++ip) {
      fprintf(fp,"%18.15e %18.15e %18.15e %d\n",xp[ip][0],xp[ip][1],xp[ip][2],flag[ip]);
    }
  }
  else {
    assert(data!=NULL && flag !=NULL);
    fprintf(fp,"VARIABLES=\"X\" \"Y\" \"Z\" \"flag\" \"data\"\n");
    for (int ip = 0; ip < np; ++ip) {
      fprintf(fp,"%18.15e %18.15e %18.15e %d %18.15e\n",xp[ip][0],xp[ip][1],xp[ip][2],flag[ip],data[ip]);
    }
  }
  fclose(fp);
}

void writeDebugFile(const string& name,const vector<int>& vec,const double (* const xp)[3]) {
  FILE * fp = fopen(name.c_str(),"w");
  for (int ii = 0; ii < vec.size(); ++ii) {
    const int ip = vec[ii];
    fprintf(fp,"%18.15e %18.15e %18.15e\n",xp[ip][0],xp[ip][1],xp[ip][2]);
  }
  fclose(fp);
}

double getFlaggedZoneArea(SurfaceShm * surface) {
  if (mpi_rank == 0) cout << " > flagged zones:" << endl;
  for (int izone = 0,nzone=surface->zoneVec.size(); izone < nzone; ++izone) {
    if (surface->zoneVec[izone].flag) {
      if (mpi_rank == 0) cout << "    > " << surface->zoneVec[izone].getName() << endl;
    }
  }
  surface->ensureStost();

  // get the total flagged area:
  double area_st = 0.0;
  for (int ist = 0; ist < surface->nst; ++ist) {
    const int izone = surface->znost[ist];
    if (surface->zoneVec[izone].flag) {
      const double normal2[3] = TRI_NORMAL_2(surface->xsp[surface->spost[ist][0]],surface->xsp[surface->spost[ist][1]],surface->xsp[surface->spost[ist][2]]);
      area_st += 0.5*MAG(normal2);
    }
  }
  if (mpi_rank == 0) cout << " > flagged surface area: " << area_st << endl;
  return area_st;
}

int groupTrisByFeatures(int * groups,SurfaceShm * const surface) {
  surface->ensureStost();

  // use 5th bit to set as visited (1-4 in use)

  int n_groups = 0;
  stack<int> next_tris;
  int * st_flag = new int[surface->nst];
  for (int ist = 0; ist < surface->nst; ++ist) {
    st_flag[ist] = groups[ist];
    groups[ist] = -1;
  }

  for (int ist = 0; ist < surface->nst; ++ist) {
    if ((st_flag[ist] >> 0) & 1) {
      // non zero means part of flagged region

      if (!((st_flag[ist] >> 4) & 1)) {
        // hasn't been processed yet so add to stack
        next_tris.push(ist);

        // loop through neighbors and add those who are not separated by a feature boundary
        while (!next_tris.empty()) {
          const int my_ist = next_tris.top(); next_tris.pop();
          groups[my_ist] = n_groups;
          st_flag[my_ist] |= (1 << 4);  // set as visited

          FOR_I3 {
            const int ist_nbr = surface->stost[my_ist][i];
            const bool edge_is_feature = ((st_flag[my_ist] >> (i+1)) & 1) ? true:false;
            const bool nbr_valid = ((st_flag[ist_nbr] >> 0) & 1) ? true:false;
            const bool nbr_visited = ((st_flag[ist_nbr] >> 4) & 1) ? true:false;
            if (nbr_valid && !nbr_visited && !edge_is_feature) {
              next_tris.push(ist_nbr);
            }
          }
        }
        ++n_groups;
      }

    }
  }

  DELETE(st_flag);

  return n_groups;
}

void identifyFeatureXp(vector<NewXp>& fixedXp,int * const st_flag,SurfaceShm * const surface,const double delta,const double feature_angle) {
  const double dp_tol = cos((180.0-feature_angle)*M_PI/180.0);

  int * sp_flag = new int[surface->nsp];  // count feature edge touches for valence
  for (int isp=0; isp<surface->nsp; ++isp) sp_flag[isp] = 0;

  surface->ensureStost();

  // compute valence at surface nodes (to determine feature endpoints/intersections) based on half-edge touches
  // considerations:
  //  - ignore points on selection boundaries (handled by ist_nbr >= 0 check)
  //  - valence of 2 is an endpoint, 6+ is an intersection
  //  - must be divisible by two because using half-edges for counting

  // also collect edges into a vector by unique node-pair (isp0<isp1, direction of edge doesn't matter)
  // need a node to edge search structure too...
  vector<TripleInt> featureEdgeVec;  // nodes of edge, ist of tri it lives on (min ist)
  multimap<int,int> ispToEdgeMap;

  // flag processed tris
  // use bits to indicate which half-edges are features AND whether this tri has been visited or not:
  // & 1: has been visited
  // & 2: 0-1 edge is a feature
  // & 4: 1-2 edge is a feature
  // & 8: 2-0 edge is a feature
  for (int ist = 0; ist < surface->nst; ++ist) st_flag[ist] = 0;

  for (int ist = 0; ist < surface->nst; ++ist) {
    if (surface->zoneVec[surface->znost[ist]].flag) {
      st_flag[ist] |= 1;
      double * n_ist = NULL;
      FOR_I3 {
        const int ist_nbr = surface->stost[ist][i];
        if (ist_nbr < 0) {
          // open edge should be protected
          const int isp0 = surface->spost[ist][i];
          const int isp1 = surface->spost[ist][(i+1)%3];
          sp_flag[isp0] += 1;
          sp_flag[isp1] += 1;

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

          // add edge to feature-edge search structures
          const int ise = featureEdgeVec.size();
          featureEdgeVec.push_back(TripleInt(min(isp0,isp1),max(isp0,isp1),ist));
          ispToEdgeMap.insert(pair<int,int> (isp0,ise));
          ispToEdgeMap.insert(pair<int,int> (isp1,ise));
        }
        else if ((ist_nbr >= 0) && (surface->zoneVec[surface->znost[ist_nbr]].flag)) {
          // determine crease angle
          if (n_ist == NULL) {
            n_ist = new double[3];
            const double n_tmp[3] = TRI_NORMAL_2(surface->xsp[surface->spost[ist][0]],surface->xsp[surface->spost[ist][1]],surface->xsp[surface->spost[ist][2]]);
            FOR_I3 n_ist[i] = n_tmp[i];
          }

          const double n_ist_nbr[3] = TRI_NORMAL_2(surface->xsp[surface->spost[ist_nbr][0]],surface->xsp[surface->spost[ist_nbr][1]],surface->xsp[surface->spost[ist_nbr][2]]);

          const double n_ist_mag = sqrt(DOT_PRODUCT(n_ist,n_ist)); assert(n_ist_mag > 0.0);
          const double n_ist_nbr_mag = sqrt(DOT_PRODUCT(n_ist_nbr,n_ist_nbr)); assert(n_ist_nbr_mag > 0.0);

          if ( DOT_PRODUCT(n_ist,n_ist_nbr)/(n_ist_mag*n_ist_nbr_mag) < dp_tol ) {
            // crease angle criterion hit
            const int isp0 = surface->spost[ist][i];
            const int isp1 = surface->spost[ist][(i+1)%3];
            sp_flag[isp0] += 1;
            sp_flag[isp1] += 1;

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
              // cout << "feature edge[" << ise << "]: " << min(isp0,isp1) << " " << max(isp0,isp1) << endl;
              featureEdgeVec.push_back(TripleInt(min(isp0,isp1),max(isp0,isp1),ist));
              ispToEdgeMap.insert(pair<int,int> (isp0,ise));
              ispToEdgeMap.insert(pair<int,int> (isp1,ise));
            }
          }
        }
      }
      delete[] n_ist;
    }
  }

  int f_nodes = 0;
  for (int isp=0; isp<surface->nsp; ++isp) {
    if (sp_flag[isp]) ++f_nodes;
  }
  cout << " > number of unique feature edges,nodes: " << featureEdgeVec.size() << " " << f_nodes << endl;

  // valence is stored in sp_flag
  // go through and register points of interest (endpoints and chain intersections)
  // loop tris again to register a single valid stoxp (even though an edge...)
  for (int ist = 0; ist < surface->nst; ++ist) {
    if (surface->zoneVec[surface->znost[ist]].flag) {
      FOR_I3 {
        const int isp=surface->spost[ist][i];
        if (sp_flag[isp] <= 0) continue;  // not of interest, so skip

        assert(sp_flag[isp]%2 == 0);
        const int valence = sp_flag[isp]/2;
        sp_flag[isp] = -valence;  // -indexed valence; negative so we don't visit again in this loop

        if (valence == 2) {
          // not an endpoint, but sits between feature edges
        }
        else {
          // store this point of interest
          fixedXp.push_back(NewXp());
          FOR_I3 fixedXp.back().xyz[i] = surface->xsp[isp][i];
          fixedXp.back().delta = delta;
          fixedXp.back().isp = isp;
          fixedXp.back().st = ist;
        }
      }
    }
  }

  const int n_endpoints = fixedXp.size();
  cout << " > high valence nodes (endpoints):  " << n_endpoints << endl;

  // now we can group our feature edges into chains by marching from these points.
  // if feature edges remain after all this marching, then loop-chains exist and
  // we must process these as well
  const int f_edges = featureEdgeVec.size();
  int * se_flag = new int[f_edges];
  for (int ise=0; ise<f_edges; ++ise) se_flag[ise] = 1;  // track which edges have already been visited

  // knowing edges are unique, we should be able to march from endpoint nodes
  // along edge-chains until all non-loops have been processed
  stack<int> isp_starts;
  for (vector<NewXp>::iterator it=fixedXp.begin(); it!=fixedXp.end(); ++it) {
    isp_starts.push(it->isp);
  }

  for (int isp=0; isp<surface->nsp; ++isp) {
    if (sp_flag[isp] == -2) sp_flag[isp] *= -1;  // non-endpoints are positive indexed, so negative are valence nodes
  }

  int nse_processed = 0;
  int isp0,isp1;
  vector<int> chain_edges;  // index is edge, +/- indicates ordering of nodes (isp_min,max or vice versa)
  int chains = 0;
  int isp_start_loop = 0;
  bool b_loop = false;
  while (nse_processed != f_edges) {

    // starting nodes are either identified endpoints or rando mnode in a loop
    int isp_start;
    if (!isp_starts.empty()) {
      // grab edges emanating from a start node
      isp_start = isp_starts.top(); isp_starts.pop();
    }
    else {
      // all remaining edges are in loops, so take a random point and march
      b_loop = true;
      for ( ; isp_start_loop<surface->nsp; ++isp_start_loop) {
        if (sp_flag[isp_start_loop] == 2) {
          --sp_flag[isp_start_loop];  // don't loop around on ourselves again
          break;
        }
      }
      isp_start = isp_start_loop;
    }
    // cout << "isp_start (isp,flag): " << isp_start << " " << sp_flag[isp_start] << endl;
    // cout << "nse_processed/f_edges: " << nse_processed << "/" << f_edges << endl;
    // getchar();

    // for each touching edge, see if processed. if not, start marching
    pair <multimap<int,int>::iterator, multimap<int,int>::iterator> ret,ret2;
    ret = ispToEdgeMap.equal_range(isp_start);
    for (multimap<int,int>::iterator it=ret.first; it!=ret.second; ++it) {
      int ise = it->second;

      if (se_flag[ise]) {

        isp0 = isp_start;

        // start of a chain
        double length = 0.0;
        chain_edges.clear();
        bool done = false;

        while (!done) {
          se_flag[ise] = 0;  // set as visited
          const int isp_min = featureEdgeVec[ise].first;
          const int isp_max = featureEdgeVec[ise].second;
          bool reverse = false;
          if (isp0 == isp_min) {
            isp1 = isp_max;
            chain_edges.push_back(ise);
          }
          else {
            assert(isp0 == isp_max);
            isp1 = isp_min;
            reverse = true;
            chain_edges.push_back(-ise-1);  // -1 indexed
          }
          // cout << "link (edgeIndex,isp0,isp1): " << ise << " " << isp0 << " " << isp1 << endl;
          length += DIST(surface->xsp[isp0],surface->xsp[isp1]);
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
            for (multimap<int,int>::iterator it2=ret2.first; it2!=ret2.second; ++it2) {
              if (it2->second != ise) {
                ise = it2->second;  // update next edge
                break;
              }
            }
          }
        }
        // now all edges and total length should have been found
        ++chains;
        cout << " > feature edge chain " << chains << " (n_edges,length): " << chain_edges.size() << " " << length << endl;
        nse_processed += int(chain_edges.size());


        // the chain links are stored in chain_edges
        // we can now insert our fixed new points along this chain
        const double edge_delta = 0.8*delta;
        if (edge_delta < length) {
          const int internal_pts = floor(length/edge_delta);
          const double effective_delta = length/internal_pts;
          double to_next = effective_delta;
          for (vector<int>::iterator link=chain_edges.begin(); link!=chain_edges.end(); ++link) {
            const int ise = *link;
            int isp0,isp1,ist;
            if (ise < 0) {
              isp1 = featureEdgeVec[-ise-1].first;
              isp0 = featureEdgeVec[-ise-1].second;
              ist  = featureEdgeVec[-ise-1].third;
            }
            else {
              isp0 = featureEdgeVec[ise].first;
              isp1 = featureEdgeVec[ise].second;
              ist  = featureEdgeVec[ise].third;
            }

            // for loops we need to add the starting node
            if (b_loop && (link==chain_edges.begin())) {
              fixedXp.push_back(NewXp());
              FOR_I3 fixedXp.back().xyz[i] = surface->xsp[isp0][i];
              fixedXp.back().delta = delta;
              fixedXp.back().st = ist;
              assert(sp_flag[isp0] == 1);
              --sp_flag[isp0];  // indicate this point already has a node on it
            }

            const double dx[3] = DIFF(surface->xsp[isp1],surface->xsp[isp0]);
            const double my_length = MAG(dx);

            double my_traversed = 0.0;
            while (my_traversed<my_length) {
              if (to_next < (my_length-my_traversed)) {
                my_traversed += to_next;
                to_next = effective_delta;

                // insert point
                const double tol = 5.0e-3;
                const double remainder = (my_length-my_traversed);  // diff between new node and edge end
                if (remainder > tol*effective_delta) {
                  fixedXp.push_back(NewXp());
                  const double frac = my_traversed/my_length;
                  FOR_I3 fixedXp.back().xyz[i] = surface->xsp[isp0][i] + frac*dx[i];
                  fixedXp.back().delta = delta;
                  fixedXp.back().st = ist;
                }
                else {
                  // pretty much on isp1
                  if (sp_flag[isp1] == 1) {
                    // only add if a link-link node that hasn't already been added
                    fixedXp.push_back(NewXp());
                    FOR_I3 fixedXp.back().xyz[i] = surface->xsp[isp1][i];
                    fixedXp.back().delta = delta;
                    fixedXp.back().st = ist;
                    --sp_flag[isp1];  // set as visited
                    to_next -= remainder;  // compensate for extra movement
                  }
                }
              }
              else {
                to_next -= (my_length-my_traversed);
                my_traversed = my_length;
              }
            }
          }
        }

        chain_edges.clear();
      }
    }
  }
  cout << " > feature edge nodes:  " << int(fixedXp.size())-n_endpoints << endl;

  delete[] se_flag;
  delete[] sp_flag;

  // if (mpi_rank == 0) {
  //   cout << " > " << fixedXp.size() << " fixed points identified and added (xyz,delta,st):" << endl;
  //   for (vector<NewXp>::iterator it=fixedXp.begin(); it!=fixedXp.end(); ++it) {
  //     cout << "    > " << COUT_VEC(it->xyz) << " " << it->delta << " " << it->st << endl;
  //   }
  // }
}

void initSurfaceVorPoints(double (*xp)[3],int* stoxp,double* deltap,int np_fixed,const int np,SurfaceShm * surface,const double delta,const double c) {

  double area_st = 0.0;  // running flagged surface area
  int ip_begin = 0;
  for (int ist = 0; ist < surface->nst; ++ist) {
    if (surface->zoneVec[surface->znost[ist]].flag) {
      const double normal2[3] = TRI_NORMAL_2(surface->xsp[surface->spost[ist][0]],surface->xsp[surface->spost[ist][1]],surface->xsp[surface->spost[ist][2]]);
      const double this_area = 0.5*MAG(normal2);
      area_st += this_area;

      const int ip_end = int(area_st/(c*delta*delta));
      if (ip_end > ip_begin) {
        // this tri gets some points...
        for (int ip_sa = ip_begin; ip_sa < ip_end; ++ip_sa) {
          const int ip = ip_sa + np_fixed;
          // uniform distribution on a tri...
          const double r0 = 0.999*double(rand())/double(RAND_MAX)+0.0005;
          const double r1 = 0.999*double(rand())/double(RAND_MAX)+0.0005;
          FOR_I3 xp[ip][i] =
            (1.0-sqrt(r0))*surface->xsp[surface->spost[ist][0]][i] +
            sqrt(r0)*(1.0-r1)*surface->xsp[surface->spost[ist][1]][i] +
            sqrt(r0)*r1*surface->xsp[surface->spost[ist][2]][i];
          stoxp[ip] = ist;
          deltap[ip] = delta*1.25;
        }
      }
      ip_begin = ip_end;
    }
  }
  assert(ip_begin == (np-np_fixed));
}

void lloydIterateSurfaceVorPoints(int& iter,double (*xp)[3],double * deltap,double (*surface_area_xp)[2],int * nboxp_i,vector<int>& nboxp_v,int *stoxp,int * xp_flag,const int maxiter,const double zero,const int np,SurfaceShm * surface, const double delta) {
  int * st_flag = new int[surface->nst];

  double (*dxp)[3] = new double[np][3]; // lloyd motion
  int *stoxp_new = new int[np];
  int * xpost_i = new int[surface->nst+1];
  int * xpost_v = new int[np];

  CuttableVoronoiData cvd;
  stack<int> stack;
  vector<pair<double,int> > nbrVec;
  map<const pair<int,int>,int> edgeMap;
  map<const int,int> nodeMap;
  int noost[3],edost[3];
  map<const int,double*> x0Map;
  set<int> nbrSet;

  // incoming flag value determines if point is allowed to move or not
  int * xp_fixed = new int[np];
  int count = 0;
  for (int ip=0; ip<np; ++ip) {
    xp_fixed[ip] = xp_flag[ip];
    if (xp_fixed[ip]) ++count;
  }
  const int np_fixed = count;

  int debug_iter = 0;
  int done = 0;
  while (done != 2) {
    ++iter;

    // build reverse-lookup xpost_i/v...
    for (int ist = 0; ist < surface->nst; ++ist) xpost_i[ist+1] = 0;
    for (int ip = 0; ip < np; ++ip) {
      const int ist = stoxp[ip];
      ++xpost_i[ist+1];
    }
    xpost_i[0] = 0;
    for (int ist = 0; ist < surface->nst; ++ist) xpost_i[ist+1] += xpost_i[ist];
    assert(xpost_i[surface->nst] == np);
    for (int ip = 0; ip < np; ++ip) {
      const int ist = stoxp[ip];
      xpost_v[xpost_i[ist]++] = ip;
    }
    for (int ist = surface->nst-1; ist > 0; --ist) xpost_i[ist] = xpost_i[ist-1];
    xpost_i[0] = 0;

    // check...

    for (int ist = 0; ist < surface->nst; ++ist) {
      for (int xos = xpost_i[ist]; xos != xpost_i[ist+1]; ++xos) {
        const int ip = xpost_v[xos];
        assert(stoxp[ip] == ist);
      }
    }

    // st_flag is used to build the surface patches. Start it off
    // with -1...

    for (int ist = 0; ist < surface->nst; ++ist) st_flag[ist] = -1;

    // loop through our ip's...
    int my_count[3] = { 0, 0, 0 };
    for (int ip = 0; ip < np; ++ip) {
      // NOTE: if you start at np_fixed, then only internal tris are built....
      //cout << "working on ip: " << ip << " stoxp: " << stoxp[ip] << endl;

      ++debug_iter;

      assert(cvd.nno == 0);
      assert(cvd.ned == 0);
      assert(nbrVec.empty());
      assert(edgeMap.empty());
      assert(nodeMap.empty());
      assert(st_flag[stoxp[ip]] != ip);
      st_flag[stoxp[ip]] = ip;
      assert(stack.empty());
      stack.push(-stoxp[ip]-1);
      while (!stack.empty()) {
        int ist = stack.top(); stack.pop();
        if (ist < 0) {
          // -----------------------------------------------------------------
          // a negative ist means both tris AND nbrs, so add the tris here...
          // -----------------------------------------------------------------
          ist = -ist-1;
          // use edge-matching to set as much of the noost and edost that we can...
          FOR_I3 noost[i] = -1;
          FOR_I3 {
            edost[i] = -1;
            // here the edgeMap is used to find edges in terms of surface node pairs...
            map<const pair<int,int>,int>::iterator it = edgeMap.find(pair<int,int>(surface->spost[ist][(i+1)%3],surface->spost[ist][i]));
            if (it != edgeMap.end()) {
              edost[i] = it->second;
              edgeMap.erase(it);
              if (noost[i] == -1) noost[i] = cvd.nooed[edost[i]][1];
              else assert(noost[i] == cvd.nooed[edost[i]][1]);

              if (noost[(i+1)%3] == -1) noost[(i+1)%3] = cvd.nooed[edost[i]][0];
              else assert(noost[(i+1)%3] == cvd.nooed[edost[i]][0]);
            }
          }
          // step 3: now noost and edost are set from edges where possible. At this point, complete any nodes not set...
          FOR_I3 {
            if (noost[i] == -1) {
              map<const int,int>::iterator it = nodeMap.find(surface->spost[ist][i]);
              if (it == nodeMap.end()) {
                const int ino = cvd.new_node();
                noost[i] = ino;
                nodeMap[surface->spost[ist][i]] = ino;
                FOR_J3 cvd.x_no[ino][j] = surface->xsp[surface->spost[ist][i]][j] - xp[ip][j];
              }
              else {
                noost[i] = it->second;
              }
            }
          }
          // step 4: finally complete the edges...
          FOR_I3 {
            if (edost[i] == -1) {
              const int ied = cvd.new_edge();
              edgeMap[pair<int,int>(surface->spost[ist][i],surface->spost[ist][(i+1)%3])] = ied;
              cvd.nooed[ied][0] = noost[i];
              cvd.nooed[ied][1] = noost[(i+1)%3];
              cvd.faoed[ied][0] = ist;
              const int ist_nbr = surface->stost[ist][i];
              if ((ist_nbr >= 0)&&(surface->zoneVec[surface->znost[ist_nbr]].flag)) cvd.faoed[ied][1] = -1;
              else cvd.faoed[ied][1] = -2; // used to mark external boundary
            }
            else {
              const int ied = edost[i];
              assert(cvd.nooed[ied][1] == noost[i]);
              assert(cvd.nooed[ied][0] == noost[(i+1)%3]);
              assert(cvd.faoed[ied][0] != ist); // no edge can have the same tri on both sides
              assert(cvd.faoed[ied][1] == -1);
              cvd.faoed[ied][1] = ist;
            }
          }
        }
        // -----------------------------------------------------------------
        // always add nbrs...
        // -----------------------------------------------------------------
        assert(st_flag[ist] == ip);
        //cout << "adding nbrs from tri: " << ist << endl;
        for (int xos = xpost_i[ist]; xos != xpost_i[ist+1]; ++xos) {
          const int ip_nbr = xpost_v[xos];
          if (ip_nbr != ip) {
            const double d2 = DIST2(xp[ip_nbr],xp[ip]);
            if (d2 == 0.0) {
              cout << "got zero: " << ip << " " << ip_nbr << endl;
              FILE * fp = fopen("xp.dat","w");
              fprintf(fp,"VARIABLES=\"X\" \"Y\" \"Z\"\n");
              fprintf(fp,"%18.15e %18.15e %18.15e\n",xp[ip][0],xp[ip][1],xp[ip][2]);
              fclose(fp);
            }
            assert(d2 > 0.0); // possible to seed 2 points at the same location?
            if (d2 < deltap[ip]*deltap[ip]) nbrVec.push_back(pair<double,int>(d2,ip_nbr));
          }
        }
        // if the neighboring tri is closer than 0.5*deltap[ip], then we need to add it AND
        // its nbrs. If it is closer than deltap[ip], then we just need to add the nbrs.
        FOR_I3 {
          const int ist_nbr = surface->stost[ist][i];
          assert(ist_nbr != ist);
          if ((ist_nbr >= 0)&&(surface->zoneVec[surface->znost[ist_nbr]].flag)&&(st_flag[ist_nbr] != ip)) {
            st_flag[ist_nbr] = ip; // do not ever revisit this tri
            // the 3 nodes of ist_nbr...
            const int isp0 = surface->spost[ist_nbr][0];
            const int isp1 = surface->spost[ist_nbr][1];
            const int isp2 = surface->spost[ist_nbr][2];
            const double d2 = MiscUtils::getPointToTriDist2(xp[ip],surface->xsp[isp0],surface->xsp[isp1],surface->xsp[isp2]);
            // if this tri is closer than 0.5*deltap[ip], we add it to the surface, and if less than
            // deltap[ip], we need to potentially add the nbrs (but not the surface)...
            if (d2 < 0.2501*deltap[ip]*deltap[ip])  {
              // both tri and nbrs...
              stack.push(-ist_nbr-1);
            }
            else if (d2 < deltap[ip]*deltap[ip])  {
              // just nbrs...
              stack.push(ist_nbr);
            }
          }
        }
      }
      nodeMap.clear();
      edgeMap.clear();

      cvd.setD2Max();
      sort(nbrVec.begin(),nbrVec.end());
      // we now have a sorted list of nbrs. Cut the cvd against the nbrs until it cannot possibly
      // be cut anymore...
      for (int in = 0; in < nbrVec.size(); ++in) {
        // recall that the .first contains the d2 nbr distance. and cvd.d2_max contains
        // the maximum node radius. as soon as d2 >= 4*cvd.d2_max, there cannot be any
        // more nbrs that cut this cv...
        if (nbrVec[in].first > 4.0*cvd.d2_max) break; // done: all remaining nbrs are too far away to possibly cut the current vd...
        const int ip_nbr = nbrVec[in].second;
        double dn[3]; FOR_I3 dn[i] = 0.5*(xp[ip_nbr][i] - xp[ip][i]);
        // only cut if we are inside the 6 paraboloids...
        if ((2.0*dn[0] > cvd.Lmax[0])||
        (2.0*dn[1] > cvd.Lmax[1])||
        (2.0*dn[2] > cvd.Lmax[2])||
        (2.0*dn[0] < cvd.Lmin[0])||
        (2.0*dn[1] < cvd.Lmin[1])||
        (2.0*dn[2] < cvd.Lmin[2]) ) {
          continue;
        }
        // if we got here, we are cutting...
        cvd.cut_surf(dn,-8-ip_nbr); // note that local nbrs use a -8 convention. cvd.d2_max updated automatically
      }
      nbrVec.clear();

      // compute the centroid of the remaining parts...
      bool has_seed_boundary = false;
      bool has_boundary = false;
      double dxc[3] = { 0.0, 0.0, 0.0 };
      double area2_sum = 0.0;
      assert(x0Map.empty());
      for (int ied = 0; ied < cvd.ned; ++ied) {
        if (cvd.faoed[ied][0] >= 0) {
          map<const int,double*>::const_iterator it = x0Map.find(cvd.faoed[ied][0]);
          if (it == x0Map.end()) {
            // on the first time we visit any face, store a node to use in all
            // the tris with other edges (the face is planar and convex, so this is fine)...
            x0Map[cvd.faoed[ied][0]] = cvd.x_no[cvd.nooed[ied][0]];
          }
          else if (cvd.x_no[cvd.nooed[ied][1]] != it->second) {
            const double normal2[3] = TRI_NORMAL_2(it->second,cvd.x_no[cvd.nooed[ied][0]],cvd.x_no[cvd.nooed[ied][1]]);
            const double area2 = MAG(normal2);
            FOR_I3 dxc[i] += area2*(it->second[i]+cvd.x_no[cvd.nooed[ied][0]][i]+cvd.x_no[cvd.nooed[ied][1]][i]);
            area2_sum += area2;
          }
        }
        if (cvd.faoed[ied][1] >= 0) {
          map<const int,double*>::const_iterator it = x0Map.find(cvd.faoed[ied][1]);
          if (it == x0Map.end()) {
            x0Map[cvd.faoed[ied][1]] = cvd.x_no[cvd.nooed[ied][1]];
          }
          else if (cvd.x_no[cvd.nooed[ied][0]] != it->second) {
            const double normal2[3] = TRI_NORMAL_2(it->second,cvd.x_no[cvd.nooed[ied][1]],cvd.x_no[cvd.nooed[ied][0]]);
            const double area2 = MAG(normal2);
            FOR_I3 dxc[i] += area2*(it->second[i]+cvd.x_no[cvd.nooed[ied][0]][i]+cvd.x_no[cvd.nooed[ied][1]][i]);
            area2_sum += area2;
          }
        }
        else if (cvd.faoed[ied][1] == -1) {
          has_seed_boundary = true;
          break;
        }
        else {
          // must be a -2: this is ok (flagged region boundary)...
          assert(cvd.faoed[ied][1] == -2);
          has_boundary = true;
        }
      }

      if (has_seed_boundary) {
        // we are too small. Increase deltap[ip], and do NOT apply any motion to the point. The remaining
        // patch is arbitrarily large...
        ++my_count[0]; xp_flag[ip] = 0;
        assert(done == 0);
        deltap[ip] *= 1.5;
        stoxp_new[ip] = stoxp[ip];
        FOR_I3 dxp[ip][i] = 0.0;
      }
      else {
        // Even though we have cut away all the original seed boundary of this patch, there still
        // may be nbrs out there that can cut some of our corners off if deltap was not large
        // enough...
        const double delta_half = sqrt(cvd.d2_max); // cvd.d2_max stores the furthest distance of any node
        if (deltap[ip] <= 2.0*delta_half) {
          // other neighbors might exist that could still cut this volume. Find them in
          // the next iteration...
          ++my_count[1]; xp_flag[ip] = 1;
          assert(done == 0);
          deltap[ip] = 2.5*delta_half; // could be just 2.0 here? - use a little larger so face checks pass for sure
          stoxp_new[ip] = stoxp[ip];
          FOR_I3 dxp[ip][i] = 0.0;
        }
        else {
          ++my_count[2]; xp_flag[ip] = 2;
          // this is good. Do some more work to decide where to move this point...
          if (done == 0) {
            // regular lloyd iteration pass
            assert(area2_sum > 0.0);
            surface_area_xp[ip][0] = sqrt(area2_sum);
            FOR_I3 dxc[i] /= area2_sum*3.0;
            deltap[ip] = 2.5*delta_half;
            // in addition to lloyd, add random perturbation...
            //FOR_I3 dxc[i] += (double(rand())/double(RAND_MAX)-0.5)*deltap[ip]*factor;
            // now we need to figure out which patch we are closest to...
            int ist_closest = -1;
            double d2_closest;
            for (int ied = 0; ied < cvd.ned; ++ied) {
              if (cvd.faoed[ied][0] >= 0) {
                map<const int,double*>::const_iterator it = x0Map.find(cvd.faoed[ied][0]);
                assert(it != x0Map.end());
                if ((cvd.x_no[cvd.nooed[ied][0]] != it->second)&&(cvd.x_no[cvd.nooed[ied][1]] != it->second)) {
                  double this_dxp[3];
                  MiscUtils::getClosestPointOnTriRobust(this_dxp,dxc,it->second,cvd.x_no[cvd.nooed[ied][0]],cvd.x_no[cvd.nooed[ied][1]]);
                  const double this_d2 = DIST2(this_dxp,dxc);
                  if ((ist_closest == -1)||(this_d2 < d2_closest)) {
                    FOR_I3 dxp[ip][i] = this_dxp[i];
                    d2_closest = this_d2;
                    ist_closest = cvd.faoed[ied][0];
                  }
                }
              }
              if (cvd.faoed[ied][1] >= 0) {
                map<const int,double*>::const_iterator it = x0Map.find(cvd.faoed[ied][1]);
                assert(it != x0Map.end());
                if ((cvd.x_no[cvd.nooed[ied][0]] != it->second)&&(cvd.x_no[cvd.nooed[ied][1]] != it->second)) {
                  double this_dxp[3];
                  MiscUtils::getClosestPointOnTriRobust(this_dxp,dxc,it->second,cvd.x_no[cvd.nooed[ied][0]],cvd.x_no[cvd.nooed[ied][1]]);
                  const double this_d2 = DIST2(this_dxp,dxc);
                  if ((ist_closest == -1)||(this_d2 < d2_closest)) {
                    FOR_I3 dxp[ip][i] = this_dxp[i];
                    d2_closest = this_d2;
                    ist_closest = cvd.faoed[ied][1];
                  }
                }
              }
            }
            // make sure we found one!...
            assert(ist_closest != -1);
            stoxp_new[ip] = ist_closest;
          }
          else {
            assert(done == 1);
            // this is a recompute of the last time...
            if (has_boundary) xp_flag[ip] = 1;
            assert(nbrSet.empty());
            for (int ied = 0; ied < cvd.ned; ++ied) {
              if (cvd.faoed[ied][0] < 0) {
                assert(cvd.faoed[ied][0] <= -8);
                nbrSet.insert(-cvd.faoed[ied][0]-8);
              }
            }
            assert(!nbrSet.empty());
            for (set<int>::const_iterator it = nbrSet.begin(); it != nbrSet.end(); ++it) nboxp_v.push_back(*it);
            nboxp_i[ip+1] = nboxp_v.size();
            nbrSet.clear();
          }
        }
      }
      // clear cvd for next time...
      x0Map.clear();
      cvd.clear();
    } // ip loop

    double dx2_max = 0.0;
    double dx2_sum = 0.0;
    for (int ip = np_fixed; ip < np; ++ip) {
      const double dx2 = DOT_PRODUCT(dxp[ip],dxp[ip]);
      dx2_max = max(dx2_max,dx2);
      dx2_sum += dx2;
    }
    cout << "iter: " << iter << " counts: " << my_count[0] << " " << my_count[1] << " " << my_count[2] << " l2, linf: " << sqrt(dx2_sum/double(np))/delta << " " << sqrt(dx2_max)/delta << endl;

    if ((done == 0)&&(my_count[0]+my_count[1]==0)&&((iter >= maxiter)||(sqrt(dx2_sum/double(np))/delta < zero))) {
      // we need to go through again and compute the nbrs...
      done = 1;
    }
    else if (done == 1) {
      // we are done...
      done = 2;
    }

    // only modify xp by dxp when done == 0. If done has been switched
    // to 1, we leave xp unmodified and recompute, storing nbr csr structures...

    if (done == 0) {
      for (int ip = 0; ip < np; ++ip)  {
        if (xp_fixed[ip] == 0) {
          stoxp[ip] = stoxp_new[ip];
          FOR_I3 xp[ip][i] += dxp[ip][i];
        }
      }

      if (iter%10 == 0) {
        writeXpToFp(iter,xp,np,xp_flag,NULL);
      }
    }
  }

  // reset xp_flag to who is fixed
  for (int ip=0; ip<np; ++ip) xp_flag[ip]=xp_fixed[ip];

  delete[] xp_fixed;
  delete[] st_flag;
  delete[] dxp;
  delete[] stoxp_new;
  delete[] xpost_i;
  delete[] xpost_v;

  {
    // output final points and area-ratio
    // compute area ratio of neighbors per ip
    double * area_ratio = new double[np];
    for (int ip = 0; ip < np; ++ip) {
      double areas[2] = {HUGE_VAL,-HUGE_VAL};
      for (int noxp = nboxp_i[ip]; noxp != nboxp_i[ip+1]; ++noxp) {
        const int ip_nbr = nboxp_v[noxp];
        areas[0] = min(areas[0],surface_area_xp[ip_nbr][0]);
        areas[1] = max(areas[1],surface_area_xp[ip_nbr][0]);
      }
      area_ratio[ip] = areas[1]/areas[0];
      surface_area_xp[ip][1] = area_ratio[ip];
    }
    writeXpToFp(iter,xp,np,xp_flag,area_ratio);
    DELETE(area_ratio);
  }
}

void jiggerXp(double (*xp)[3],int * xp_flag,const int np,const double jigger_tol,const double (*surface_area_xp)[2],const SurfaceShm * surface,const double * deltap) {
  //HACK assumes will stay on same tri... should re-project at this point but don't have a good way to constrain tri search...
  const double jigger_factor = 0.25;
  for (int ip = 0; ip < np; ++ip)  {
    if ((xp_flag[ip] == 0)&&(surface_area_xp[ip][1]>jigger_tol)) {
      FOR_I3 xp[ip][i] += (0.5-(double(rand())/double(RAND_MAX)))*deltap[ip]*jigger_factor;

    }
    else if (xp_flag[ip] == 0) {
      // should be fixed based on our tolerance, so update "fixed" flag
      xp_flag[ip] = 1;
    }
  }
}

void triangulateSurfaceVorPoints(vector<TripleInt>& triVec,int (**spost)[3],int& nst,int& np_internal,vector<int>& xpotl_i,vector<int>& xpotl_v,int& ntl,int * nboxp_i,vector<int>& nboxp_v,double (**xp)[3],double (*normalp)[3],int& np,SurfaceShm * surface,const bool skip_loops) {

  // building nboxp structure
  for (int ip0 = 0; ip0 < np; ++ip0) {
    for (int nox0 = nboxp_i[ip0]; nox0 != nboxp_i[ip0+1]; ++nox0) {
      const int ip1 = nboxp_v[nox0];
      if (ip1 > ip0) {
        int nox1;
        for (nox1 = nboxp_i[ip1]; nox1 != nboxp_i[ip1+1]; ++nox1) {
          const int ip2 = nboxp_v[nox1];
          if (ip2 == ip0) break;
        }
        if (nox1 != nboxp_i[ip1+1]) {
          for (nox1 = nboxp_i[ip1]; nox1 != nboxp_i[ip1+1]; ++nox1) {
            const int ip2 = nboxp_v[nox1];
            if (ip2 > ip1) {
              // check for ip0 and ip1 in ip2's nbrs...
              int count = 0;
              for (int nox2 = nboxp_i[ip2]; nox2 != nboxp_i[ip2+1]; ++nox2) {
                if ((nboxp_v[nox2] == ip0)||(nboxp_v[nox2] == ip1)) {
                  ++count;
                  if (count == 2) {
                    const double this_normal2[3] = TRI_NORMAL_2((*xp)[ip0],(*xp)[ip1],(*xp)[ip2]);
                    double dp = DOT_PRODUCT(this_normal2,normalp[ip0]);
                    dp += DOT_PRODUCT(this_normal2,normalp[ip1]);
                    dp += DOT_PRODUCT(this_normal2,normalp[ip2]);
                    if (dp > 0.0) {
                      triVec.push_back(TripleInt(ip0,ip1,ip2));
                    }
                    else {
                      triVec.push_back(TripleInt(ip0,ip2,ip1));
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
  cout << " > done building nboxp and new Tris" << endl;

  const int ntr = triVec.size();
  int * troxp_i = new int[np+1];
  for (int ip = 0; ip < np; ++ip) troxp_i[ip+1] = 0;
  for (int it = 0; it < ntr; ++it) {
    ++troxp_i[triVec[it].first+1];
    ++troxp_i[triVec[it].second+1];
    ++troxp_i[triVec[it].third+1];
  }
  troxp_i[0] = 0;
  for (int ip = 0; ip < np; ++ip) troxp_i[ip+1] += troxp_i[ip];
  int * troxp_v = new int[troxp_i[np]];
  for (int it = 0; it < ntr; ++it) {
    troxp_v[troxp_i[triVec[it].first]++] = it;
    troxp_v[troxp_i[triVec[it].second]++] = it;
    troxp_v[troxp_i[triVec[it].third]++] = it;
  }
  for (int ip = np-1; ip > 0; --ip) troxp_i[ip] = troxp_i[ip-1];
  troxp_i[0] = 0;

  // now find the edges in the tri surface...

  // int ntl = 0; // number of tri loops...
  // vector<int> xpotl_i;
  // xpotl_i.push_back(0);
  // vector<int> xpotl_v;

  cout << "starting tri loop computation" << endl;

  int * xp_flag = new int[np];
  for (int ip = 0; ip < np; ++ip) xp_flag[ip] = -1;

  for (int ip = 0; ip < np; ++ip) {
    if (xp_flag[ip] == -1) {
      int ip1 = ip;
      int done = 0;
      while (done == 0) {
        done = 1;
        const int ip0 = ip1;
        for (int tox = troxp_i[ip0]; tox != troxp_i[ip0+1]; ++tox) {
          const int it = troxp_v[tox];
          if (triVec[it].first == ip0) {
            ip1 = triVec[it].second;
          }
          else if (triVec[it].second == ip0) {
            ip1 = triVec[it].third;
          }
          else {
            assert(triVec[it].third == ip0);
            ip1 = triVec[it].first;
          }
          // if this edge ip0->ip1 has another tri that has ip1->ip0, then
          // we skip it...
          int tox2;
          for (tox2 = troxp_i[ip0]; tox2 != troxp_i[ip0+1]; ++tox2) {
            if (tox2 != tox) {
              const int it2 = troxp_v[tox2];
              if ( ((triVec[it2].first == ip1)&&(triVec[it2].second == ip0)) ||
                   ((triVec[it2].second == ip1)&&(triVec[it2].third == ip0)) ||
                   ((triVec[it2].third == ip1)&&(triVec[it2].first == ip0)) )
                break;
            }
          }
          if (tox2 == troxp_i[ip0+1]) {
            // there was no other tri in the opposite dirction, so this must be
            // an edge. For now, do not consider bow-ties...
            if (ip1 >= np_internal) assert(xp_flag[ip1] == -1);  //HACK

            xp_flag[ip1] = ntl;
            xpotl_v.push_back(ip1);
            if (ip1 == ip) {
              xpotl_i.push_back(xpotl_v.size());
              ++ntl;
            }
            else {
              done = 0;
            }
            break;
          }
        }
      }
    }
  }

  delete[] xp_flag;
  cout << "got ntl (number of tri loops): " << ntl << endl;

  {
    FILE * fp = fopen("tri_edges.dat","w");
    fprintf(fp,"VARIABLES=\"X\" \"Y\" \"Z\"\n");
    for (int itl = 0; itl < ntl; ++itl) {
      for (int xol = xpotl_i[itl]; xol != xpotl_i[itl+1]; ++xol) {
        const int ip = xpotl_v[xol];
        fprintf(fp,"%18.15e %18.15e %18.15e\n",(*xp)[ip][0],(*xp)[ip][1],(*xp)[ip][2]);
      }
    }
    fclose(fp);
  }

  // ===================================================================
  // and get the surface loops, if any...
  // ===================================================================

  int * sp_flag = new int[surface->nsp];
  // -1: un-touched
  // -2 and below: -2-indexed next node along edge
  for (int isp = 0; isp < surface->nsp; ++isp) sp_flag[isp] = -1;

  // sum boundary edge touches at each node: expect 0 or 2...
  for (int ist = 0; ist < surface->nst; ++ist) {
    if (surface->zoneVec[surface->znost[ist]].flag) {
      FOR_I3 {
        const int ist_nbr = surface->stost[ist][i];
        if ((ist_nbr < 0)||(surface->zoneVec[surface->znost[ist_nbr]].flag == 0)) {
          // this is a boundary edge (no nbr or unflagged zone)
          assert(sp_flag[surface->spost[ist][i]] == -1);
          sp_flag[surface->spost[ist][i]] = -surface->spost[ist][(i+1)%3]-2; // use -2 indexing
        }
      }
    }
  }

  int nsl = 0; // n surface loops
  vector<int> sposl_i; // surface-points of surface-loop
  sposl_i.push_back(0);
  vector<int> sposl_v;
  for (int isp = 0; isp < surface->nsp; ++isp) {
    if (sp_flag[isp] < -1) {
      int isp1 = isp;
      int done = 0;
      while (done == 0) {
        const int isp0 = isp1;
        isp1 = -sp_flag[isp0]-2;
        sp_flag[isp0] = nsl;
        sposl_v.push_back(isp0);
        if (isp1 == isp) {
          sposl_i.push_back(sposl_v.size());
          ++nsl;
          done = 1;
        }
      }
    }
  }
  delete[] sp_flag;
  cout << "nsl (number of surface loops): " << nsl << endl;


  {
    FILE * fp = fopen("surface_edges.dat","w");
    fprintf(fp,"VARIABLES=\"X\" \"Z\" \"Y\"\n");
    for (int isl = 0; isl < nsl; ++isl) {
      for (int sos = sposl_i[isl]; sos != sposl_i[isl+1]; ++sos) {
        const int isp = sposl_v[sos];
        fprintf(fp,"%18.15e %18.15e %18.15e\n",surface->xsp[isp][0],surface->xsp[isp][1],surface->xsp[isp][2]);
      }
    }
    fclose(fp);
  }

  assert(nsl == ntl);

  // now decide who goes with who...

  vector<pair<int,int> > loopPairVec;
  if ((ntl == 0)&&(nsl == 0)) {
    // leave loopPairVec empty...
  }
  else if ((ntl == 1)&&(nsl == 1)) {
    // one loop goes with one loop...
    loopPairVec.push_back(pair<int,int>(0,0));
  }
  else {
    // figure out who goes with who...
    assert(0);
  }

  // expand xp to include all edge points...
  if (!skip_loops) {

    np_internal = np;
    np += sposl_i[nsl];
    if (np > np_internal) {
      double (*xp_old)[3] = *xp;
      *xp = new double[np][3];
      for (int ip = 0; ip < np_internal; ++ip) {
        FOR_I3 (*xp)[ip][i] = xp_old[ip][i];
      }
      delete[] xp_old;
    }
    for (int sos = 0; sos < sposl_i[nsl]; ++sos) {
      const int isp = sposl_v[sos];
      FOR_I3 (*xp)[np_internal+sos][i] = surface->xsp[isp][i];
    }
    nst = ntr + xpotl_i[ntl] + sposl_i[nsl];
    *spost = new int[nst][3];
    for (int itr = 0; itr < ntr; ++itr) {
      (*spost)[itr][0] = triVec[itr].first;
      (*spost)[itr][1] = triVec[itr].second;
      (*spost)[itr][2] = triVec[itr].third;
    }
    vector<int> vec1;
    vector<int> vec2;
    int offset = ntr;

    for (int ii = 0; ii < loopPairVec.size(); ++ii) {
      const int itl = loopPairVec[ii].first;
      const int isl = loopPairVec[ii].second;
      vec1.resize(xpotl_i[itl+1]-xpotl_i[itl]);
      for (int xot = xpotl_i[itl]; xot != xpotl_i[itl+1]; ++xot) vec1[xot-xpotl_i[itl]] = xpotl_v[xot];
      vec2.resize(sposl_i[isl+1]-sposl_i[isl]);
      for (int sos = sposl_i[isl]; sos != sposl_i[isl+1]; ++sos) vec2[sposl_i[isl+1]-sos-1] = np_internal + sos;
      const int count = GeomUtils::facetGap(*spost+offset,vec1,vec2,*xp,true); // true: loop
      assert(count == vec1.size()+vec2.size());
      offset += count;
      // take a look...
      writeDebugFile("vec1.dat",vec1,*xp);
      writeDebugFile("vec2.dat",vec2,*xp);
      cout << "checkout vec1 and vec2 files" << endl;
      getchar();
    }
    assert(offset == nst);
  }
  else {
    // don't do edge processing
    nst = ntr;
      *spost = new int[nst][3];
      for (int itr = 0; itr < ntr; ++itr) {
        (*spost)[itr][0] = triVec[itr].first;
        (*spost)[itr][1] = triVec[itr].second;
        (*spost)[itr][2] = triVec[itr].third;
      }
  }

  GeomUtils::writeTecplot("tris.dat",*spost,nst,*xp);
  GeomUtils::writeSbin("tris.sbin",*spost,nst,*xp);

  cout << " > done internal triangulation" << endl;
}

void buildHcp2d(SurfaceShm * surface,const double delta,const int iters,const double feature_angle,const int jiggers,const bool skip_loops) {

  srand(0);

  if (mpi_rank == 0) cout << "buildHcp2d: delta: " << delta << endl;

  const double area_st = getFlaggedZoneArea(surface);
  if (!(area_st > 0.0)) {
    if (mpi_rank == 0) cerr << "no flagged area; skipping" << endl;
    return;
  }

  // identify points that should be fixed, i.e., the constrained features
  vector<NewXp> fixedXp;
  int * st_flag = new int[surface->nst];
  identifyFeatureXp(fixedXp,st_flag,surface,delta,feature_angle);
  const int np_fixed = fixedXp.size();
  if (mpi_rank == 0) cout << " > np_features: " << np_fixed << endl;

  // fixed feature nodes have been inserted
  // st_flag contains bits which tell the edges of tris with features; can be used to color tris
  // by disjoint feature-patch
  const int n_groups = groupTrisByFeatures(st_flag,surface);
  // st_flag now contains disjoint feature group index of tris


  FILE * fp = fopen("init_tris_colored.dat","w");
  fprintf(fp,"VARIABLES = \"X\"\n");
  fprintf(fp,"\"Y\"\n");
  fprintf(fp,"\"Z\"\n");
  for (int ig = 0; ig<n_groups; ++ig) {
    GeomUtils::writeTecplotZone(fp,surface->spost,surface->nst,surface->xsp,st_flag,ig);
  }
  fclose(fp);


  // the number of points to add on the surface is equal to area_st/(c*delta*delta)
  const double c = 0.67;
  int np_surface = int(area_st/(c*delta*delta));
  // TODO: eventually split up np_surface over ranks. For now, one rank...
  if (mpi_rank == 0) cout << " > np_surface: " << np_surface << endl;


  int np = np_fixed + np_surface;
  double (*xp)[3] = new double[np][3]; // location of the points
  int *stoxp = new int[np]; // the tri they live on
  double *deltap = new double[np]; // their max distance to nbrs

  int * xp_flag = new int[np];
  for (int ip=0; ip<np; ++ip) xp_flag[ip] = 0;

  int count=0;
  for (vector<NewXp>::iterator it=fixedXp.begin(); it!=fixedXp.end(); ++it,++count) {
    FOR_I3 xp[count][i] = it->xyz[i];
    deltap[count] = it->delta;
    stoxp[count] = it->st;
    xp_flag[count] = 1;  // indicates point is fixed
  }
  assert(count == np_fixed);

  initSurfaceVorPoints(xp,stoxp,deltap,np_fixed,np,surface,delta,c);
  writeXpToFp(0,xp,np,xp_flag,NULL);  // output initial strand locations

  // xp neighbor information for eventual triangulation
  int * nboxp_i = new int[np+1];
  nboxp_i[0] = 0;
  vector<int> nboxp_v;
  nboxp_v.reserve(np*6);

  //double factor = 0.25;
  double (*surface_area_xp)[2] = new double[np][2];  // potentially use to jigger points selectively and re-lloyd

  double zero = 0.001;
  int iters_performed = 0;
  lloydIterateSurfaceVorPoints(iters_performed,xp,deltap,surface_area_xp,nboxp_i,nboxp_v,stoxp,xp_flag,(iters_performed+iters),zero,np,surface,delta);

  const double jigger_tol = 1.025;
  for (int ijig=0; ijig < jiggers; ++ijig) {
    jiggerXp(xp,xp_flag,np,jigger_tol,surface_area_xp,surface,deltap);
    lloydIterateSurfaceVorPoints(iters_performed,xp,deltap,surface_area_xp,nboxp_i,nboxp_v,stoxp,xp_flag,(iters_performed+iters),zero,np,surface,delta);
  }
  DELETE(surface_area_xp);

  cout << "triangulating..." << endl;

  // tri-based...

  // add tris from the connectivity that are in increasing node order
  // (so they are only added once)...

  // to get orientation correct, we need the normal of the surface
  // patch we are living on...
  //TODO how does this work for points forced to crease edges?
  double (*normalp)[3] = new double[np][3];
  for (int ip = 0; ip < np; ++ip) {
    const int ist = stoxp[ip];
    const double normal2[3] = TRI_NORMAL_2(surface->xsp[surface->spost[ist][0]],surface->xsp[surface->spost[ist][1]],surface->xsp[surface->spost[ist][2]]);
    const double nmag = MAG(normal2);
    assert(nmag > 0.0);
    FOR_I3 normalp[ip][i] = normal2[i]/nmag;
  }

  int (*spost)[3] = NULL;
  int nst,np_internal;

  int ntl = 0; // number of tri loops...
  vector<int> xpotl_i;
  xpotl_i.push_back(0);
  vector<int> xpotl_v;

  vector<TripleInt> triVec;

  triangulateSurfaceVorPoints(triVec,&spost,nst,np_internal,xpotl_i,xpotl_v,ntl,nboxp_i,nboxp_v,&xp,normalp,np,surface,skip_loops);

  if (skip_loops) return;

  if (np>np_internal) {
    // need to resize xp_flag after potential addition of boundary
    delete[] xp_flag;
    xp_flag = new int[np];
  }

  // recall we have the normalp allocated for all internal ip's...

  for (int ip = 0; ip < np_internal; ++ip) {
    FOR_I3 normalp[ip][i] = 0.0;
  }
  for (int ist = 0; ist < nst; ++ist) {
    const double this_normal[3] = TRI_NORMAL_2(xp[spost[ist][0]],xp[spost[ist][1]],xp[spost[ist][2]]);
    FOR_I3 {
      const int ip = spost[ist][i];
      if (ip < np_internal) {
        FOR_J3 normalp[ip][j] += this_normal[j];
      }
    }
  }
  for (int ip = 0; ip < np_internal; ++ip) {
    const double mag = MAG(normalp[ip]);
    assert(mag > 0.0);
    FOR_I3 normalp[ip][i] /= mag;
  }

  // now use the xp_flag to store the local level...

  // everyone gets nlayers
  int nlayers = 3;
  for (int ip = 0; ip < np_internal; ++ip) xp_flag[ip] = nlayers;

  // set nlayers to 1 in points surrounding edges...

  for (int xot = 0; xot < xpotl_i[ntl]; ++xot) {
    const int ip = xpotl_v[xot];
    assert((ip >= 0)&&(ip < np_internal));
    xp_flag[ip] = 1;
  }

  // smooth this into surrounding layers...

  int done = 0;
  while (done == 0) {
    done = 1;
    for (int ist = 0,ntr=triVec.size(); ist < ntr; ++ist) {
      int layers_min = nlayers;
      int layers_max = 0;
      FOR_I3 {
        const int isp = spost[ist][i];
        assert((isp >= 0)&&(isp < np_internal));
        layers_min = min(layers_min,xp_flag[isp]);
        layers_max = max(layers_max,xp_flag[isp]);
      }
      if (layers_max > layers_min+1) {
        layers_max = layers_min+1;
        done = 0;
        FOR_I3 {
          const int isp = spost[ist][i];
          xp_flag[isp] = min(xp_flag[isp],layers_min+1);
        }
      }
      // also consider how folded this tri has become at layers_max,
      // and reduce if necessary...
      if (layers_max >= 1) {
        bool failed_checks = false;
        // loop through tri edges...
        FOR_I3 {
          const int isp0 = spost[ist][i];
          const int isp1 = spost[ist][(i+1)%3];
          const double dx[3] = DIFF(xp[isp1],xp[isp0]);
          double dx_ff[3]; FOR_I3 dx_ff[i] = dx[i] - double(layers_max)*delta*(normalp[isp1][i]-normalp[isp0][i]);
          const double dx_mag = MAG(dx); assert(dx_mag > 0.0);
          const double dx_ff_mag = MAG(dx_ff);
          if ((dx_ff_mag < dx_mag*0.5)||(dx_ff_mag > dx_mag*1.5)) {
            failed_checks = true;
            break;
          }
          /*
          const double dp = DOT_PRODUCT(dx_ff,dx)/(dx_mag*dx_ff_mag);
          if ((dp < 0.8)||(dp > 1.2)) {
            failed_checks = true;
            break;
          }
          */
        }
        if (failed_checks) {
          done = 0;
          // set everyone to layers_max-1...
          FOR_I3 {
            const int isp = spost[ist][i];
            xp_flag[isp] = layers_max-1;
          }
        }
      }
    }
  }

  double (*xp_ff)[3] = new double[np][3];
  for (int ip = 0; ip < np_internal; ++ip) {
    FOR_I3 xp_ff[ip][i] = xp[ip][i] - double(xp_flag[ip])*delta*normalp[ip][i];
  }
  for (int ip = np_internal; ip < np; ++ip) {
    FOR_I3 xp_ff[ip][i] = xp[ip][i];
  }

  GeomUtils::writeSbin("tris_ff.sbin",spost,nst,xp_ff);
  GeomUtils::writeTecplot("tris_ff.dat",spost,nst,xp_ff);

  cout << "check out sbin: SO FAR SO GOOD" << endl;
  getchar();


  // write out the mesh...
  fp = fopen("surf.dat","w");
  fprintf(fp,"TITLE = \"%s\"\n","triangles");
  fprintf(fp,"VARIABLES = \"X\"\n");
  fprintf(fp,"\"Y\"\n");
  fprintf(fp,"\"Z\"\n");
  fprintf(fp,"ZONE T=\"%s\"\n","main");
  fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",np,triVec.size());
  for (int ip = 0; ip < np; ++ip) {
    fprintf(fp,"%18.15e %18.15e %18.15e\n",xp[ip][0],xp[ip][1],xp[ip][2]);
  }
  for (int it = 0; it < triVec.size(); ++it) {
    fprintf(fp,"%d %d %d\n",triVec[it].first+1,triVec[it].second+1,triVec[it].third+1);
  }
  fclose(fp);

  cout << "checkout surf.dat" << endl;
  getchar();



  // look for edge nodes...



  for (int ip = 0; ip < np; ++ip) {
    if (xp_flag[ip] == 1) {
      // how many edges do we touch?...
      int boundary_nbr_count = 0;
      for (int nox = nboxp_i[ip]; nox != nboxp_i[ip+1]; ++nox) {
        const int ip_nbr = nboxp_v[nox];
        assert(ip_nbr != ip);
        if (xp_flag[ip_nbr] == 1)
          ++boundary_nbr_count;
      }
      assert(boundary_nbr_count == 2);
    }
  }







  MPI_Pause("WE ARE DONE");


}




int main(int argc, char * argv[]) {

  try {

    CTI_Init(argc,argv,"testHcp2d.in");

    {

      SurfaceShm * surface = NULL;

      Param * param = getParam("SURF");
      assert(param);
      int iarg = 0;
      while (iarg < param->size()) {
        string token = param->getString(iarg++);
        if ((token.size() > 5)&&(token.compare(token.size()-5,5,".sbin") == 0)) {
          surface = new SurfaceShm(token);
        }
        else {
          if (mpi_rank == 0) cout << "WARNING: skipping unrecognized token: " << token << endl;
        }
      }

      param = getParam("TEST_HCP2D");
      assert(param);
      double delta = -1.0;
      int iters = 50;
      int jiggers = 0;
      double feature_angle = 150.0;
      bool skip_loops = false;
      iarg = 0;
      while (iarg < param->size()) {
        string token = param->getString(iarg++);
        if ((token == "ZONE")||(token == "FAZONE")) {
          token = param->getString(iarg++);
          surface->flagZonesMatching(token);
        }
        else if (token == "DELTA") {
          delta = param->getDouble(iarg++);
          cout << "GOT delta: " << delta << endl;
        }
        else if (token == "ITERS") {
          iters = param->getInt(iarg++);
          cout << "GOT iters: " << iters << endl;
        }
        else if (token == "JIGGERS") {
          jiggers = param->getInt(iarg++);
          cout << "GOT jiggers: " << jiggers << endl;
        }
        else if (token == "FEATURE_ANGLE") {
          feature_angle = param->getDouble(iarg++);
          cout << "GOT feature angle: " << feature_angle << endl;
        }
        else if (token == "SKIP_LOOPS") {
          skip_loops = true;
          cout << "GOT skip surface/new-tri loops" << endl;
        }
        else {
          if (mpi_rank == 0) cout << "WARNING: skipping unrecognized token: " << token << endl;
        }
      }
      assert(delta != -1);
      buildHcp2d(surface,delta,iters,feature_angle,jiggers,skip_loops);







      MPI_Pause("SOOOO FARR SO GOOOOOOOOOOD");






      // cleanup...
      if (surface) delete surface;

    }

    CTI_Finalize();

  }
  catch (int e) {
    if (e == 0)
      CTI_Finalize();
    else
      CTI_Abort();
  }
  catch(...) {
    CTI_Abort();
  }

  return(0);

}
