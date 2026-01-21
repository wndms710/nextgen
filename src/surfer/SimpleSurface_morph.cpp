#include "SimpleSurface.hpp"
#include "WebUI.hpp"

void SimpleSurface::morphTranslateFlaggedZones(const double dx[3]) {
  if (sp_buf == NULL) sp_buf = new double[nsp];
  double dx_mag = MAG(dx);
  for (int isp = 0; isp < nsp; ++isp) sp_buf[isp] = 3.0*dx_mag;
  double bbmin[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
  double bbmax[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
  for (int ist = 0; ist < nst; ++ist) {
    if (zoneVec[znost[ist]].flag) {
      FOR_I3 {
	const int isp = spost[ist][i];
	sp_buf[isp] = 0.0;
	FOR_J3 {
	  bbmin[j] = min(bbmin[j],xsp[isp][j]);
	  bbmax[j] = max(bbmax[j],xsp[isp][j]);
	}
      }
    }
  }
  // we now have a bounding box. enlarge to include a set of tris we are going to work on...
  FOR_I3 bbmin[i] -= 2.0*dx_mag;
  FOR_I3 bbmax[i] += 2.0*dx_mag;
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_buf[isp] == 3.0*dx_mag) {
      if ((xsp[isp][0] >= bbmin[0])&&
	  (xsp[isp][0] <= bbmax[0])&&
	  (xsp[isp][1] >= bbmin[1])&&
	  (xsp[isp][1] <= bbmax[1])&&
	  (xsp[isp][2] >= bbmin[2])&&
	  (xsp[isp][2] <= bbmax[2])) {
	sp_buf[isp] = 2.0*dx_mag;
      }
    }
  }
  // now the tris we are going to consider must have ay least 2 
  // nodes less than or equal to 2.0*dx_mag and atleast one greater
  // than zero...
  vector<int> stVec;
  for (int ist = 0; ist < nst; ++ist) {
    int count = 0;
    FOR_I3 {
      int isp = spost[ist][i];
      if (sp_buf[isp] <= 2.0*dx_mag)
	++count;
    }
    if (count >= 2) {
      if ((sp_buf[spost[ist][0]] > 0.0)||
	  (sp_buf[spost[ist][1]] > 0.0)||
	  (sp_buf[spost[ist][2]] > 0.0)) {
	stVec.push_back(ist);
      }
    }
  }
  bool done = false;
  while (!done) {
    done = true;
    for (int ii = 0,size = stVec.size(); ii < size; ++ii) {
      const int ist = stVec[ii];
      FOR_I3 {
	int isp = spost[ist][i];
	// 1D distance first...
	FOR_J2 {
	  const int isp_nbr = spost[ist][(i+j)%3];
	  if (sp_buf[isp] > sp_buf[isp_nbr]) {
	    const double sp_check = sp_buf[isp_nbr]+DIST(xsp[isp],xsp[isp_nbr]);
	    if (sp_check < sp_buf[isp]) {
	      done = false;
	      sp_buf[isp] = sp_check;
	    }
	  }
	}
      }
      FOR_I3 {
	int isp = spost[ist][i];
	// 2D distance next...
	if ((sp_buf[isp] > sp_buf[spost[ist][(i+1)%3]])&&(sp_buf[isp] > sp_buf[spost[ist][(i+2)%3]])) {
	  const int isp_nbr1 = spost[ist][(i+1)%3];
	  const int isp_nbr2 = spost[ist][(i+2)%3];
	  const double dx_nbr_mag = DIST(xsp[isp_nbr1],xsp[isp_nbr2]);
	  const double dudx = (sp_buf[isp_nbr2]-sp_buf[isp_nbr1])/dx_nbr_mag;
	  const double dudy2 = sqrt(1.0 - dudx*dudx);
	  if (dudy2 > 0.0) {
	    const double dudy = sqrt(dudy2);
	    const double dx[3] = DIFF(xsp[isp],xsp[isp_nbr1]);
	    const double dx_nbr[3] = DIFF(xsp[isp_nbr2],xsp[isp_nbr1]);
	    const double dxp = DOT_PRODUCT(dx,dx_nbr)/dx_nbr_mag;
	    if (dxp > 0.0) {
	      const double dyp = sqrt(DOT_PRODUCT(dx,dx)-dxp*dxp);
	      const double sp_check = sp_buf[isp_nbr1]+dudx*dxp+dudy*dyp;
	      if ((sp_check < sp_buf[isp])&&(sp_check > sp_buf[isp_nbr1])&&(sp_check > sp_buf[isp_nbr2])) {
		done = false;
		sp_buf[isp] = sp_check;
	      }
	    }
	  }
	}
      }
    }
  }
  // now compute the morph factor in sp_buf and morph the points...
  for (int isp = 0; isp < nsp; ++isp) {
    sp_buf[isp] = max(0.0,(2.0*dx_mag-sp_buf[isp])/(2.0*dx_mag));
    if (sp_buf[isp] > 0.0) {
      FOR_I3 xsp[isp][i] += sp_buf[isp]*dx[i];
    }
  }
  return;
}

bool SimpleSurface::checkMorphedZone(const int zone) {

  // need teost to check for local folding
  ensureTeost();

  // check to see if neighboring points are inside my tri
  bool bad_zone = false;
  // double sum_n[3] = { 0.0, 0.0, 0.0 };
  // double sum_n0[3] = { 0.0, 0.0, 0.0 };
  for (int toz = stosz_i[zone]; toz != stosz_i[zone+1]; ++toz) {
    const int ist = stosz_v[toz];
    // TODO replace with zero area check and crease-angle check
    assert(szost[ist] == zone);
    FOR_I3 {
      const double* a = xsp[spost[ist][i]];
      const double* b = xsp[spost[ist][(i+1)%3]];
      const double* c = xsp[spost[ist][(i+2)%3]];
      double per = 0.0;
      FOR_J3 per += DIST(xsp[spost[ist][(i+j)%3]],xsp[spost[ist][(i+1+j)%3]]);
      const double n_st[3] = TRI_NORMAL_2(a,b,c);
      if (MAG(n_st) > 1.0E-6*per*per) { // skip zero area tris
        int ist_nbr, i_nbr, orient_nbr;
        if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
          // is this correct?
          if (orient_nbr)
            assert(spost[ist][i] == spost[ist_nbr][i_nbr] && spost[ist][(i+1)%3] == spost[ist_nbr][(i_nbr+1)%3]);
          else
            assert(spost[ist][i] == spost[ist_nbr][(i_nbr+1)%3] && spost[ist][(i+1)%3] == spost[ist_nbr][i_nbr]);
          const double* d = xsp[spost[ist_nbr][(i_nbr+2)%3]];
          if (pointInTriangle(a,b,c,d)) {
            double x_st[3];
            FOR_I3 x_st[i] = (a[i]+b[i]+c[i])/3.0;
            cout << "Warning: neighbor triangle point is surface triangle at x = "
                 << COUT_VEC(x_st) << " n " << COUT_VEC(n_st) << " on zone " << zone
                 << ". Please try again with different parameters." << endl;
            bad_zone = true;
            break;
          }
        }
      }
    }
    if (bad_zone) break;
  }

  return bad_zone;

}

void SimpleSurface::morphTranslate(const set<int>& moved, const set<int>& freed,const double dx[3], const int nrelax, const bool use_bump) {

  ensureStoszSzosz();

  sp_flag.setLength(nsp);

  // init points to freed (-1)...
  FOR_ISP sp_flag[isp] = -1;

  // store
  double (*xsp0)[3] = new double[nsp][3];
  FOR_ISP FOR_I3 xsp0[isp][i] = xsp[isp][i];

  // set points on moved to 1...
  for (set<int>::iterator it = moved.begin(); it != moved.end(); ++it) {
    const int isz = *it;
    for (int toz = stosz_i[isz]; toz != stosz_i[isz+1]; ++toz) {
      const int ist = stosz_v[toz];
      FOR_I3 {
        const int isp = spost[ist][i];
        if (sp_flag[isp] == -1) {
          sp_flag[isp] = 1;
        }
      }
    }
  }

  // set points on frozen (not moved/free) to 0.
  // this overwrite shared free pts (boundaries have precedent)...
  for (int isz = 0; isz < nsz; ++isz) {
    if (moved.find(isz) != moved.end()) continue;
    if (freed.find(isz) != freed.end()) continue;
    for (int toz = stosz_i[isz]; toz != stosz_i[isz+1]; ++toz) {
      const int ist = stosz_v[toz];
      FOR_I3 {
        const int isp = spost[ist][i];
        if (sp_flag[isp] == -1) {
          sp_flag[isp] = 0;
        }
      }
    }
  }

  // build sposp...
  ensureSposp();

  // solve Laplace's equation using relaxation method...
  if (sp_buf == NULL) sp_buf = new double[nsp];
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] != 1)
      sp_buf[isp] = 0.0;
    else
      sp_buf[isp] = 1.0;
  }
  for (int ii = 0; ii < nrelax; ++ii) {
    for (int isp = 0; isp < nsp; ++isp) {
      if (sp_flag[isp] == -ii-1) {
        sp_flag[isp] -= 1;
        sp_buf[isp] = 0.0;
        // note that I don't include myself
        for (int pop = sposp_i[isp]+1; pop != sposp_i[isp+1]; ++pop) {
          const int isp2 = sposp_v[pop];
          //sp_buf[isp] += sp_buf[isp2]*sp_buf[isp2];
          sp_buf[isp] += sp_buf[isp2];
        }
        if ((sposp_i[isp+1]-(sposp_i[isp]+1)) > 1) {
          sp_buf[isp] /= sposp_i[isp+1]-(sposp_i[isp]+1);
        }
        //sp_buf[isp] = sqrt(sp_buf[isp]);
        //sp_buf[isp] = exp(1.0-sp_buf[isp]);
      }
    }
  }

  //writeSpDataToTecplot("laplacian_solve.plt",sp_buf);

  // apply tranformation...
  if (use_bump) {
    for (int isp = 0; isp < nsp; ++isp) {
      FOR_I3 xsp[isp][i] += dx[i]*exp(1.0-1.0/(1.0-(1.0-sp_buf[isp])*(1.0-sp_buf[isp])));
    }
  }
  else {
    for (int isp = 0; isp < nsp; ++isp) {
      FOR_I3 xsp[isp][i] += dx[i]*sp_buf[isp];
    }
  }

  // check...
  bool bad_morph = false;
  for (set<int>::iterator it = moved.begin(); it != moved.end(); ++it) {
    const int isz = *it;
    bad_morph = checkMorphedZone(isz);
    if (bad_morph) break;
  }
  if (!bad_morph) {
    for (set<int>::iterator it = freed.begin(); it != freed.end(); ++it) {
      const int isz = *it;
      bad_morph = checkMorphedZone(isz);
      if (bad_morph) break;
    }
  }

  if (bad_morph) {
    // revert back so you can redo the morph
    for (int isp = 0; isp < nsp; ++isp) FOR_I3 xsp[isp][i] = xsp0[isp][i];
  }
  else if (b_subzone_data) {
    // surface will change, recompute geoms
    clearSubzoneData();
    clearZoneData();
    b_centroid = false;
    b_rmax = false;
    b_boundingBox = false;
    clearNonManifoldData();

    // if we had periodic data, we need to reset it here
    if (pbi) {
      clearPeriodicity();
      WUI(WARN,"Bad news, friend. You need to reset periodicity.");
    }
  }

  delete[] xsp0;
}

void SimpleSurface::morphRotateHalfCylinder(const int cyl_id, const int top_id, const int base_id, const int side_id, const double phi,const int nrelax) {

  // need stosz and teost (implicit)
  ensureStoszSzosz();

  assert((cyl_id >= 0)&&(cyl_id < nsz));
  assert((top_id >= 0)&&(top_id < nsz));
  assert((base_id >= 0)&&(base_id < nsz));
  assert((side_id >= 0)&&(side_id < nsz));

  // make sure cylinder/top, cylinder/bottom, cylinder/side are connected in the graph...
  {
    int sos;
    for (sos = szosz_i[cyl_id]; sos != szosz_i[cyl_id+1]; ++sos) {
      if (szosz_v[sos] == top_id) break;
    }
    if (sos == szosz_i[cyl_id+1]) {
      cout << "Warning: CYLINDER AND TOP are not connected subsurfaces. Skipping." << endl;
      return;
    }
    for (sos = szosz_i[cyl_id]; sos != szosz_i[cyl_id+1]; ++sos) {
      if (szosz_v[sos] == base_id) break;
    }
    if (sos == szosz_i[cyl_id+1]) {
      cout << "Warning: CYLINDER AND BASE are not connected subsurfaces. Skipping." << endl;
      return;
    }
    for (sos = szosz_i[cyl_id]; sos != szosz_i[cyl_id+1]; ++sos) {
      if (szosz_v[sos] == side_id) break;
    }
    if (sos == szosz_i[cyl_id+1]) {
      cout << "Warning: CYLINDER AND SIDE are not connected subsurfaces. Skipping." << endl;
      return;
    }
  }

  // will use flag for marking subzone boundaries and to demark sides of square prism
  sp_flag.setLength(nsp);
  sp_flag.setAll(-1);

  // start by building vectors of top and bottom edges shared with cyl
  set<int> topSpSet;
  set<int> botSpSet;
  for (int toz = stosz_i[cyl_id]; toz != stosz_i[cyl_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == cyl_id);
    FOR_I3 {
      int ist_nbr, i_nbr, orient_nbr;
      if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
        if (szost[ist_nbr] == top_id) {
          topSpSet.insert(spost[ist][i]);
          topSpSet.insert(spost[ist][(i+1)%3]);
        }
        else if (szost[ist_nbr] == base_id) {
          botSpSet.insert(spost[ist][i]);
          botSpSet.insert(spost[ist][(i+1)%3]);
        }
      }
    }
  }
  for (int toz = stosz_i[top_id]; toz != stosz_i[top_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == top_id);
    FOR_I3 {
      int ist_nbr, i_nbr, orient_nbr;
      if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
        if (szost[ist_nbr] == cyl_id) {
          topSpSet.insert(spost[ist][i]);
          topSpSet.insert(spost[ist][(i+1)%3]);
        }
      }
    }
  }
  for (int toz = stosz_i[base_id]; toz != stosz_i[base_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == base_id);
    FOR_I3 {
      int ist_nbr, i_nbr, orient_nbr;
      if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
        if (szost[ist_nbr] == cyl_id) {
          botSpSet.insert(spost[ist][i]);
          botSpSet.insert(spost[ist][(i+1)%3]);
        }
      }
    }
  }
  vector<int> topSpVec(topSpSet.size());
  vector<int> botSpVec(botSpSet.size());
  int index = 0;
  for (set<int>::iterator iter = topSpSet.begin(); iter != topSpSet.end(); ++iter) {
    topSpVec[index++] = (*iter);
  }
  index = 0;
  for (set<int>::iterator iter = botSpSet.begin(); iter != botSpSet.end(); ++iter) {
    botSpVec[index++] = (*iter);
  }
  botSpSet.clear();
  topSpSet.clear();
  // TODO just stick with sets instead of copying into vector

  const int n_points_top_loop = topSpVec.size();
  double x_top[3] = {0.0,0.0,0.0};
  int min_x_top_id = -1;
  int max_x_top_id = -1;
  double min_x_top = HUGE_VAL, max_x_top = -HUGE_VAL;
  for (int ii = 0; ii < n_points_top_loop; ++ii) {
    const int isp = topSpVec[ii];
    FOR_J3 x_top[j] += xsp[isp][j];
    if (xsp[isp][0] > max_x_top) {
      max_x_top_id = isp;
      max_x_top    = xsp[isp][0];
    }
    if (xsp[isp][0] < min_x_top) {
      min_x_top_id = isp;
      min_x_top    = xsp[isp][0];
    }
  }
  FOR_J3 x_top[j] /= double(n_points_top_loop);
  const double r_top = 0.5*DIST(xsp[max_x_top_id],xsp[min_x_top_id]);
  x_top[2] = 0.5*(xsp[max_x_top_id][2]+xsp[min_x_top_id][2]);

  const int n_points_bot_loop = botSpVec.size();
  double x_bot[3] = {0.0,0.0,0.0};
  int min_x_bot_id = -1;
  int max_x_bot_id = -1;
  double min_x_bot = HUGE_VAL, max_x_bot = -HUGE_VAL;
  for (int ii = 0; ii < n_points_bot_loop; ++ii) {
    const int isp = botSpVec[ii];
    FOR_J3 x_bot[j] += xsp[isp][j];
    if (xsp[isp][0] > max_x_bot) {
      max_x_bot_id = isp;
      max_x_bot    = xsp[isp][0];
    }
    if (xsp[isp][0] < min_x_bot) {
      min_x_bot_id = isp;
      min_x_bot    = xsp[isp][0];
    }
  }
  FOR_J3 x_bot[j] /= double(n_points_bot_loop);
  const double r_bot = 0.5*DIST(xsp[max_x_bot_id],xsp[min_x_bot_id]);
  x_bot[2] = 0.5*(xsp[max_x_bot_id][2]+xsp[min_x_bot_id][2]);

  // calculate the central axis
  double dx_bot_to_top[3] = DIFF(x_top,x_bot);
  // const double dx_bot_to_top_mag2 = DOT_PRODUCT(dx_bot_to_top,dx_bot_to_top);

  const double H0 = MAG(dx_bot_to_top);
  const double R0 = 0.5*(r_top+r_bot);

  // report the cylinder characterization
  cout << "top: x_top[3] = " << COUT_VEC(x_top) << "bot: x_bot[3] = " << COUT_VEC(x_bot) << endl;
  cout << "central axis: x_top[3]-x_bot[3] = " << COUT_VEC(dx_bot_to_top) << " H0 " << H0 << endl;
  cout << "radii: R0 = " << R0 << ", r_top = " << r_top << ", r_bot = " << r_bot << endl;

  // make sure our spread is less than 5%
  assert(fabs(r_top-r_bot)/R0 < 0.5);

  if (phi == 0.0) return;

  // right now this doesn't do anything after checking we have a cylinder
  // XXXX add morphing here
  double (*xsp0)[3] = new double[nsp][3];
  for (int isp = 0; isp < nsp; ++isp) FOR_I3 xsp0[isp][i] = xsp[isp][i];

  // rotate us such that the cylinder is y and down the plane is x
  // const double cos_angle = dx_bot_to_top[0]/sqrt(dx_bot_to_top[0]*dx_bot_to_top[0]+dx_bot_to_top[1]*dx_bot_to_top[1]);
  // const double sin_angle = dx_bot_to_top[1]/sqrt(dx_bot_to_top[0]*dx_bot_to_top[0]+dx_bot_to_top[1]*dx_bot_to_top[1]);

  // now translate and rotate cylinder...

  for (int isp = 0; isp < nsp; ++isp) sp_flag[isp] = -1;
  for (int toz = stosz_i[cyl_id]; toz != stosz_i[cyl_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == cyl_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == -1) {
        sp_flag[isp] = 0;
        // cylinder frame...
        FOR_J3 xsp[isp][j] -= x_bot[j];
        // const double x0[2] = {xsp[isp][0],xsp[isp][1]};
//        xsp[isp][0] = x0[0]*cos_angle-x0[1]*sin_angle;
//        xsp[isp][1] = x0[0]*sin_angle+x0[1]*cos_angle;
      }
    }
  }
  for (int toz = stosz_i[top_id]; toz != stosz_i[top_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == top_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == -1) {
        sp_flag[isp] = 0;
        // cylinder frame...
        FOR_J3 xsp[isp][j] -= x_bot[j];
        // const double x0[2] = {xsp[isp][0],xsp[isp][1]};
//        xsp[isp][0] = x0[0]*cos_angle-x0[1]*sin_angle;
//        xsp[isp][1] = x0[0]*sin_angle+x0[1]*cos_angle;
      }
    }
  }

  // get h's ...

  if (sp_buf3 == NULL) sp_buf3 = new double[nsp][3];
  for (int isp = 0; isp < nsp; ++isp) {
    sp_flag[isp] = -1;
    FOR_J3 sp_buf3[isp][j] = 0.0;
  }
  for (int toz = stosz_i[cyl_id]; toz != stosz_i[cyl_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == cyl_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == -1) {
        // find closest point on top edge...
        double min_dist = HUGE_VAL;
        int ii_top = -1;
        for (int ii = 0; ii < n_points_top_loop; ++ii) {
          const int isp_top = topSpVec[ii];
          const double this_dist = DIST(xsp[isp_top],xsp[isp]);
          if (this_dist < min_dist) {
            min_dist = this_dist;
            ii_top = ii;
          }
        }
        assert(ii_top >= 0);
        const int isp_top = topSpVec[ii_top];
        FOR_J3 sp_buf3[isp][j] = xsp[isp][j]-xsp[isp_top][j];
        sp_flag[isp] = ii_top; // store closest top node index here
      }
    }
  }

  double x_top_cyl[3] = {0.0,0.0,0.0};
  double x_bot_cyl[3] = {0.0,0.0,0.0};
  for (int ii = 0; ii < n_points_top_loop; ++ii)  {
    const int isp_top = topSpVec[ii];
    FOR_I3 x_top_cyl[i] += xsp[isp_top][i];
  }
  FOR_I3 x_top_cyl[i] /= (double)n_points_top_loop;
  for (int ii = 0; ii < n_points_bot_loop; ++ii)  {
    const int isp_bot = botSpVec[ii];
    FOR_I3 x_bot_cyl[i] += xsp[isp_bot][i];
  }
  FOR_I3 x_bot_cyl[i] /= (double)n_points_bot_loop;
  double dx_top_to_bot_cyl[3] = DIFF(x_bot_cyl,x_top_cyl);
  cout << "orientation: " << COUT_VEC(dx_top_to_bot_cyl)  << endl;

  for (int isp = 0; isp < nsp; ++isp) {
    const double hdx = DOT_PRODUCT(sp_buf3[isp],dx_top_to_bot_cyl);
    const double dxdx = DOT_PRODUCT(dx_top_to_bot_cyl,dx_top_to_bot_cyl);
    FOR_J3 sp_buf3[isp][j] -= hdx/dxdx*dx_top_to_bot_cyl[j];
  }

  // now rescale cylinder...

  for (int toz = stosz_i[cyl_id]; toz != stosz_i[cyl_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == cyl_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] >= 0) {
        const double h_over_h0 = 1.0+SGN(dx_bot_to_top[1])*xsp[isp][0]*tan(phi)/H0;
        const int isp_top = topSpVec[sp_flag[isp]];
        const double h = DIST(xsp[isp_top],xsp[isp]);
        FOR_J2 xsp[isp][j] = xsp[isp_top][j]+h_over_h0*h*dx_top_to_bot_cyl[j]/MAG(dx_top_to_bot_cyl)+sp_buf3[isp][j];
        //FOR_J3 xsp[isp][j] = xsp[isp_top][j]+h_over_h0*h*dx_top_to_bot_cyl[j]/MAG(dx_top_to_bot_cyl)+sp_buf3[isp][j];
        sp_flag[isp] = -1;
      }
    }
  }

  // back to lab frame
  for (int isp = 0; isp < nsp; ++isp) sp_flag[isp] = -1;
  for (int toz = stosz_i[cyl_id]; toz != stosz_i[cyl_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == cyl_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == -1) {
        sp_flag[isp] = 0;
        // const double x0[2] = {xsp[isp][0],xsp[isp][1]};
//        xsp[isp][0] = x0[0]*cos_angle+x0[1]*sin_angle;
//        xsp[isp][1] =-x0[0]*sin_angle+x0[1]*cos_angle;
        // phi
        const double x1[2] = {xsp[isp][0],xsp[isp][1]};
        xsp[isp][0] = x1[0]*cos(phi)-x1[1]*sin(phi);
        xsp[isp][1] = x1[0]*sin(phi)+x1[1]*cos(phi);
        FOR_J3 xsp[isp][j] += x_bot[j];
      }
    }
  }
  for (int toz = stosz_i[top_id]; toz != stosz_i[top_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == top_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == -1) {
        sp_flag[isp] = 0;
        // const double x0[2] = {xsp[isp][0],xsp[isp][1]};
//        xsp[isp][0] = x0[0]*cos_angle+x0[1]*sin_angle;
//        xsp[isp][1] =-x0[0]*sin_angle+x0[1]*cos_angle;
        // phi
        const double x1[2] = {xsp[isp][0],xsp[isp][1]};
        xsp[isp][0] = x1[0]*cos(phi)-x1[1]*sin(phi);
        xsp[isp][1] = x1[0]*sin(phi)+x1[1]*cos(phi);
        FOR_J3 xsp[isp][j] += x_bot[j];
      }
    }
  }

  // build sposp...
  ensureSposp();

  // now we just need to redistribute base to prevent folding...

  for (int isp = 0; isp < nsp; ++isp) sp_flag[isp] = 0; // imovable
  for (int toz = stosz_i[base_id]; toz != stosz_i[base_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == base_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == 0) sp_flag[isp] = -1; // interior
    }
  }
  for (int toz = stosz_i[side_id]; toz != stosz_i[side_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == side_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == 0) sp_flag[isp] = -1; // interior
    }
  }
  for (int ii = 0; ii < n_points_bot_loop; ++ii)  {
    const int isp = botSpVec[ii];
    if (sp_flag[isp] == -1)
      sp_flag[isp] = 1;
  }
  for (int toz = stosz_i[side_id]; toz != stosz_i[side_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == side_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      const int isp2 = spost[ist][i+1];
      int ist_nbr, i_nbr, orient_nbr;
      if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
        if (szost[ist_nbr] == cyl_id || szost[ist_nbr] == top_id) {
          if (sp_flag[isp] == -1) sp_flag[isp] = 1;
          if (sp_flag[isp] == -1) sp_flag[isp2] = 1;
        }
      }
    }
  }

  // solve Laplace's equation using relaxation method...
  for (int isp = 0; isp < nsp; ++isp) FOR_I3 sp_buf3[isp][i] = 0.0;
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] != 1)
      FOR_I3 sp_buf3[isp][i] = 0.0;
    else
      FOR_I3 sp_buf3[isp][i] = xsp[isp][i]-xsp0[isp][i];
  }
  for (int ii = 0; ii < nrelax; ++ii) {
    for (int toz = stosz_i[base_id]; toz != stosz_i[base_id+1]; ++toz) {
      const int ist = stosz_v[toz];
      assert(szost[ist] == base_id);
      FOR_I3 {
        const int isp = spost[ist][i];
        if (sp_flag[isp] == -ii-1) {
          sp_flag[isp] -= 1;
          FOR_J3 sp_buf3[isp][j] = 0.0;
          for (int pop = sposp_i[isp]+1; pop != sposp_i[isp+1]; ++pop) {
            const int isp2 = sposp_v[pop];
            FOR_J3 sp_buf3[isp][j] += sp_buf3[isp2][j];
          }
          if ((sposp_i[isp+1]-(sposp_i[isp]+1)) > 1) {
            FOR_J3 sp_buf3[isp][j] /= sposp_i[isp+1]-(sposp_i[isp]+1);
          }
        }
      }
    }
    for (int toz = stosz_i[side_id]; toz != stosz_i[side_id+1]; ++toz) {
      const int ist = stosz_v[toz];
      assert(szost[ist] == side_id);
      FOR_I3 {
        const int isp = spost[ist][i];
        if (sp_flag[isp] == -ii-1) {
          sp_flag[isp] -= 1;
          FOR_J3 sp_buf3[isp][j] = 0.0;
          for (int pop = sposp_i[isp]+1; pop != sposp_i[isp+1]; ++pop) {
            const int isp2 = sposp_v[pop];
            FOR_J3 sp_buf3[isp][j] += sp_buf3[isp2][j];
          }
          if ((sposp_i[isp+1]-(sposp_i[isp]+1)) > 1) {
            FOR_J3 sp_buf3[isp][j] /= sposp_i[isp+1]-(sposp_i[isp]+1);
          }
        }
      }
    }
  }

  /*
  cout << "write out laplacian_solve.plt" << endl;
  {
    for (int isp = 0; isp < nsp; ++isp) {
      sp_buf[isp] = (double)sp_flag[isp];
      sp_buf[isp] = SGN(sp_buf[isp]);
    }
   writeSpDataToTecplot("laplacian_solve.plt",sp_buf);
  }
  */
//  writeSpDataToTecplot("laplacian_solve.plt",sp_buf3);

  // apply tranformation...
  for (int isp = 0; isp < nsp; ++isp) {
 //   if (sp_flag[isp] < 0) FOR_I3 xsp[isp][i] += sp_buf3[isp][i];
    if (sp_flag[isp] < 0) FOR_I2 xsp[isp][i] += sp_buf3[isp][i];
  }

  // check the morph...
  bool bad_top = checkMorphedZone(top_id);
  bool bad_bot = checkMorphedZone(base_id);
  bool bad_cyl = checkMorphedZone(cyl_id);
  bool bad_side = checkMorphedZone(side_id);

  if (bad_top || bad_bot || bad_cyl || bad_side) {
     // revert back so you can redo the morph
     FOR_ISP FOR_I3 xsp[isp][i] = xsp0[isp][i];
  }
  else {
    // surface will change, recompute geoms
    clearSubzoneData();
    clearZoneData();
    b_centroid = false;
    b_rmax = false;
    b_boundingBox = false;
    clearNonManifoldData();

    // if we had periodic data, we need to reset it here
    if (pbi) {
      clearPeriodicity();
      WUI(WARN,"Bad news, friend. You need to reset periodicity.");
    }
  }
  delete[] xsp0;

}

void SimpleSurface::morphRotateCylinder(const int cyl_id, const int top_id, const int base_id, const double phi,const int nrelax) {

  // need stosz and teost (implicit)
  ensureStoszSzosz();

  assert((cyl_id >= 0)&&(cyl_id < nsz));
  assert((top_id >= 0)&&(top_id < nsz));
  assert((base_id >= 0)&&(base_id < nsz));

  // make sure cylinder/top and cylinder/bottom are connected in the graph...
  {
    int sos;
    for (sos = szosz_i[cyl_id]; sos != szosz_i[cyl_id+1]; ++sos) {
      if (szosz_v[sos] == top_id) break;
    }
    if (sos == szosz_i[cyl_id+1]) {
      cout << "Warning: CYLINDER AND TOP are not connected subsurfaces. Skipping." << endl;
      return;
    }
    for (sos = szosz_i[cyl_id]; sos != szosz_i[cyl_id+1]; ++sos) {
      if (szosz_v[sos] == base_id) break;
    }
    if (sos == szosz_i[cyl_id+1]) {
      cout << "Warning: CYLINDER AND BASE are not connected subsurfaces. Skipping." << endl;
      return;
    }
  }

  // will use flag for marking subzone boundaries and to demark sides of square prism
  sp_flag.setLength(nsp);
  sp_flag.setAll(-1);
  vector<int> topEdgeVec;
  vector<int> botEdgeVec;

  // start by building vectors of top and bottom edges shared with cyl
  for (int toz = stosz_i[cyl_id]; toz != stosz_i[cyl_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == cyl_id);
    FOR_I3 {
      int ist_nbr, i_nbr, orient_nbr;
      if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
        if (szost[ist_nbr] == top_id) {
          topEdgeVec.push_back(ist*3+i);
          // put our edge index in one of our nodes
          if (sp_flag[spost[ist][i]] == -1) {
            sp_flag[spost[ist][i]] = ist*3+i;
          }
          else {
            // more than 2 open edges must touch this node
            cout << "Warning: top edge loop algorithm problem. Skipping." << endl;
            return;
          }
        }
        else if (szost[ist_nbr] == base_id) {
          botEdgeVec.push_back(ist*3+i);
          // put our edge index in one of our nodes
          if (sp_flag[spost[ist][i]] == -1) {
            sp_flag[spost[ist][i]] = ist*3+i;
          }
          else {
            // more than 2 bottom edges must touch this node (or top/bot loops touch)
            cout << "Warning: bottom edge loop algorithm problem. Skipping." << endl;
            return;
          }
        }
      }
    }
  }

  // Check that the top and bottom loops are closed by walking them. also get the loop centers.
  // TODO replace simple average with centroid
  int n_edges_top_loop = 0;
  double x_top[3] = { 0.0, 0.0, 0.0 };
  vector<int> auxEdgeVec(topEdgeVec.size());
  for (int ii = 0,ii_end=topEdgeVec.size(); ii < ii_end; ++ii) {
    int ied = topEdgeVec[ii];
    int ist = ied/3;
    int i = ied-ist*3;
    int ist_nbr,i_nbr,orient_nbr;
    if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
      if (szost[ist_nbr] == top_id) {
        FOR_J3 x_top[j] += xsp[spost[ist][i]][j];
        szost[ist_nbr] = -top_id-1;
        auxEdgeVec[n_edges_top_loop] = ied;
        n_edges_top_loop++;
        while (sp_flag[spost[ist][(i+1)%3]] != -1) {
          ied = sp_flag[spost[ist][(i+1)%3]];
          ist = ied/3;
          i = ied-ist*3;
          if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
            if (szost[ist_nbr] == top_id) {
              FOR_J3 x_top[j] += xsp[spost[ist][i]][j];
              szost[ist_nbr] = -top_id-1;
              auxEdgeVec[n_edges_top_loop] = ied;
              n_edges_top_loop++;
            }
            else {
              assert(szost[ist_nbr] == -top_id-1);
              break;
            }
          }
        }
      }
    }
  }
  assert(n_edges_top_loop == int(topEdgeVec.size()));
  FOR_I3 x_top[i] /= n_edges_top_loop;
  for (int ii = 0; ii < n_edges_top_loop; ++ii) {
    topEdgeVec[ii] = auxEdgeVec[ii];
    const int ied = topEdgeVec[ii];
    const int ist = ied/3;
    const int i = ied-ist*3;
    int ist_nbr,i_nbr,orient_nbr;
    if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
      assert(szost[ist_nbr] == -top_id-1);
      szost[ist_nbr] = top_id;
    }
  }

  // now repeat for bottom loop
  int n_edges_bot_loop = 0;
  double x_bot[3] = { 0.0, 0.0, 0.0 };
  auxEdgeVec.resize(botEdgeVec.size());
  for (int ii = 0,ii_end=botEdgeVec.size(); ii < ii_end; ++ii) {
    int ied = botEdgeVec[ii];
    int ist = ied/3;
    int i = ied-ist*3;
    int ist_nbr,i_nbr,orient_nbr;
    if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
      if (szost[ist_nbr] == base_id) {
        FOR_J3 x_bot[j] += xsp[spost[ist][i]][j];
        szost[ist_nbr] = -base_id-1;
        auxEdgeVec[n_edges_bot_loop] = ied;
        n_edges_bot_loop++;
        while (sp_flag[spost[ist][(i+1)%3]] != -1) {
          ied = sp_flag[spost[ist][(i+1)%3]];
          ist = ied/3;
          i = ied-ist*3;
          if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
            if (szost[ist_nbr] == base_id) {
              FOR_J3 x_bot[j] += xsp[spost[ist][i]][j];
              szost[ist_nbr] = -base_id-1;
              auxEdgeVec[n_edges_bot_loop] = ied;
              n_edges_bot_loop++;
            }
            else {
              assert(szost[ist_nbr] == -base_id-1);
              break;
            }
          }
        }
      }
    }
  }

  assert(n_edges_bot_loop == int(botEdgeVec.size()));
  FOR_I3 x_bot[i] /= n_edges_bot_loop;
  for (int ii = 0; ii < n_edges_bot_loop; ++ii) {
    botEdgeVec[ii] = auxEdgeVec[ii];
    const int ied = botEdgeVec[ii];
    const int ist = ied/3;
    const int i = ied-ist*3;
    int ist_nbr,i_nbr,orient_nbr;
    if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
      assert(szost[ist_nbr] == -base_id-1);
      szost[ist_nbr] = base_id;
    }
  }

  // calculate the central axis
  double dx_bot_to_top[3] = DIFF(x_top,x_bot);
  const double dx_bot_to_top_mag2 = DOT_PRODUCT(dx_bot_to_top,dx_bot_to_top);

  // calculate the min and max extensions between top and bot
  double min_extent = +1E+20;
  double max_extent = -1E+20;
  for (int itop = 0; itop < n_edges_top_loop; ++itop) {
    const int ist = topEdgeVec[itop]/3;
    const int i = topEdgeVec[itop]-ist*3;
    // const double * const x1 = xsp[spost[ist][i]];
    for (int ibot = 0; ibot < n_edges_bot_loop; ++ibot) {
      const int ist_ = botEdgeVec[ibot]/3;
      const int i_ = botEdgeVec[ibot]-ist_*3;
      const double dx[3] = DIFF(xsp[spost[ist][i]],xsp[spost[ist_][i_]]);
      const double extent = DOT_PRODUCT(dx,dx_bot_to_top)/sqrt(dx_bot_to_top_mag2);
      min_extent = min(min_extent,extent);
      max_extent = max(max_extent,extent);
    }
  }

  // now lets walk them again using ordered top/botEdgeVec to calculate the loop normals and areas
  double n_top[3] = { 0.0, 0.0, 0.0 }; // 2x area weighted
  double area_top_2 = 0.0; // 2x the area
  for (int ii = 0; ii < n_edges_top_loop; ++ii) {
    int iip1 = ii+1;
    if (iip1 == n_edges_top_loop) iip1 = 0;
    const int ist = topEdgeVec[ii]/3;
    const int i = topEdgeVec[ii]-ist*3;
    const int istp1 = topEdgeVec[iip1]/3;
    const int ip1 = topEdgeVec[iip1]-istp1*3;
    const double * const x0 = x_top;
    const double * const x1 = xsp[spost[ist][i]];
    const double * const x2 = xsp[spost[istp1][ip1]];
    const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
    FOR_I3 n_top[i] += this_n[i];
    const double this_area = MAG(this_n);
    area_top_2 += this_area;
  }
  // make sure walks and normals are pointing out of cylinder
  if (DOT_PRODUCT(n_top,dx_bot_to_top) < 0.0) {
    reverse(topEdgeVec.begin(),topEdgeVec.end());
    FOR_I3 n_top[i] *= -1;
  }

  // now repeat for bottom loop
  double n_bot[3] = { 0.0, 0.0, 0.0 }; // 2x area weighted
  double area_bot_2 = 0.0; // 2x the area
  for (int ii = 0; ii < n_edges_bot_loop; ++ii) {
    int iip1 = ii+1;
    if (iip1 == n_edges_bot_loop) iip1 = 0;
    const int ist = botEdgeVec[ii]/3;
    const int i = botEdgeVec[ii]-ist*3;
    const int istp1 = botEdgeVec[iip1]/3;
    const int ip1 = botEdgeVec[iip1]-istp1*3;
    const double * const x0 = x_bot;
    const double * const x1 = xsp[spost[ist][i]];
    const double * const x2 = xsp[spost[istp1][ip1]];
    const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
    FOR_I3 n_bot[i] += this_n[i];
    const double this_area = MAG(this_n);
    area_bot_2 += this_area;
  }
  // make sure walks and normals are pointing out of cylinder
  if (DOT_PRODUCT(n_bot,dx_bot_to_top) > 0.0) {
    reverse(botEdgeVec.begin(),botEdgeVec.end());
    FOR_I3 n_bot[i] *= -1;
  }

  // check how close top_cyl surfaces approximate a cylinder by getting min/max/avg r_mag
  double r_max_2 = -1E+20;
  double r_min_2 = +1E+20;
  double r_avg_2 = 0.0;
 // for (int isp = 0; isp < nsp; ++isp) sp_flag[isp] = -1;
  for (int toz = stosz_i[cyl_id]; toz != stosz_i[cyl_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == cyl_id);
    FOR_I3 {
      const int isp = spost[ist][i];
//      if (sp_flag[isp] == -1) {
//        sp_flag[isp] = 0;
        const double * const this_x =xsp[isp];
        const double dx_bot[3] = DIFF(this_x,x_bot);
        const double dx_top[3] = DIFF(this_x,x_top);
        double r_vec[3];
        // calc r using the closest of x_bot and x_top
        if (DOT_PRODUCT(dx_bot,dx_bot) <= DOT_PRODUCT(dx_top,dx_top)) {
          const double dx_frac = DOT_PRODUCT(dx_bot,dx_bot_to_top)/dx_bot_to_top_mag2;
          FOR_J3 r_vec[j] = dx_bot[j] - dx_frac * dx_bot_to_top[j];
        }
        else {
          const double dx_frac = DOT_PRODUCT(dx_top,dx_bot_to_top)/dx_bot_to_top_mag2;
          FOR_J3 r_vec[j] = dx_top[j] - dx_frac * dx_bot_to_top[j];
        }
        const double r_mag_2 = DOT_PRODUCT(r_vec,r_vec);
        r_max_2 = max(r_max_2,r_mag_2);
        r_min_2 = min(r_min_2,r_mag_2);
        r_avg_2 += r_mag_2;
  //    }
    }
  }
  r_avg_2 /= (3*(stosz_i[cyl_id+1]-stosz_i[cyl_id]));
  double r_max = sqrt(r_max_2);
  double r_min = sqrt(r_min_2);
  double r_avg = sqrt(r_avg_2);

  // report the cylinder characterization
  cout << "top: x_top[3] = " << COUT_VEC(x_top) << ", 2 A n_top[3] = " << COUT_VEC(n_top) << ", A_top = " << 0.5*area_top_2 << endl;
  cout << "bot: x_bot[3] = " << COUT_VEC(x_bot) << ", 2 A n_bot[3] = " << COUT_VEC(n_bot) << ", A_bot = " << 0.5*area_bot_2 << endl;
  cout << "central axis: x_top[3]-x_bot[3] = " << COUT_VEC(dx_bot_to_top) << endl;
  cout << "radii: mean(r) = " << r_avg << ", min(r) = " << r_min << ", max(r) = " << r_max << endl;
  cout << "extent: h1 = " << min_extent << ", h2 = " << max_extent << endl;
  cout << "volume: 1/2*(h1+h2)*pi*r^2 = " << 0.5*(min_extent+max_extent)*M_PI*r_avg*r_avg << endl;

  // make sure our spread is less than 5%
  assert(fabs(r_max-r_min)/r_avg < 0.05);

  if (phi == 0.0) return;

  // right now this doesn't do anything after checking we have a cylinder
  // XXXX add morphing here
  double (*xsp0)[3] = new double[nsp][3];
  for (int isp = 0; isp < nsp; ++isp) FOR_I3 xsp0[isp][i] = xsp[isp][i];

  // now we are going to perform the transform...

  const double H0 = MAG(dx_bot_to_top);
  // const double R0 = r_avg;

  // rotate us such that the cylinder is y and down the plane is x
  // const double cos_angle = dx_bot_to_top[0]/sqrt(dx_bot_to_top[0]*dx_bot_to_top[0]+dx_bot_to_top[1]*dx_bot_to_top[1]);
  // const double sin_angle = dx_bot_to_top[1]/sqrt(dx_bot_to_top[0]*dx_bot_to_top[0]+dx_bot_to_top[1]*dx_bot_to_top[1]);

  // now translate and rotate cylinder...

  for (int isp = 0; isp < nsp; ++isp) sp_flag[isp] = -1;
  for (int toz = stosz_i[cyl_id]; toz != stosz_i[cyl_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == cyl_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == -1) {
        sp_flag[isp] = 0;
        // cylinder frame...
        FOR_J3 xsp[isp][j] -= x_bot[j];
        // const double x0[2] = {xsp[isp][0],xsp[isp][1]};
//        xsp[isp][0] = x0[0]*cos_angle-x0[1]*sin_angle;
//        xsp[isp][1] = x0[0]*sin_angle+x0[1]*cos_angle;
      }
    }
  }
  for (int toz = stosz_i[top_id]; toz != stosz_i[top_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == top_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == -1) {
        sp_flag[isp] = 0;
        // cylinder frame...
        FOR_J3 xsp[isp][j] -= x_bot[j];
        // const double x0[2] = {xsp[isp][0],xsp[isp][1]};
//        xsp[isp][0] = x0[0]*cos_angle-x0[1]*sin_angle;
//        xsp[isp][1] = x0[0]*sin_angle+x0[1]*cos_angle;
      }
    }
  }

  // get h's ...

  if (sp_buf3 == NULL) sp_buf3 = new double[nsp][3];
  for (int isp = 0; isp < nsp; ++isp) {
    sp_flag[isp] = -1;
    FOR_J3 sp_buf3[isp][j] = 0.0;
  }
  for (int toz = stosz_i[cyl_id]; toz != stosz_i[cyl_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == cyl_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == -1) {
        // find closest point on top edge...
        double min_dist = HUGE_VAL;
        int ii_top = -1;
        for (int ii = 0; ii < n_edges_top_loop; ++ii) {
          const int ist_top = topEdgeVec[ii]/3;
          const int i_top = topEdgeVec[ii]-ist_top*3;
          const double this_dist = DIST(xsp[spost[ist_top][i_top]],xsp[isp]);
          if (this_dist < min_dist) {
            min_dist = this_dist;
            ii_top = ii;
          }
        }
        assert(ii_top >= 0);
        const int ist_top = topEdgeVec[ii_top]/3;
        const int i_top = topEdgeVec[ii_top]-ist_top*3;
        FOR_J3 sp_buf3[isp][j] = xsp[isp][j]-xsp[spost[ist_top][i_top]][j];
        sp_flag[isp] = ii_top; // store closest top node index here
      }
    }
  }

  double x_top_cyl[3] = {0.0,0.0,0.0};
  double x_bot_cyl[3] = {0.0,0.0,0.0};
  for (int ii = 0; ii < n_edges_top_loop; ++ii)  {
    const int ist_top = topEdgeVec[ii]/3;
    const int i_top = topEdgeVec[ii]-ist_top*3;
    FOR_I3 x_top_cyl[i] += xsp[spost[ist_top][i_top]][i];
  }
  FOR_I3 x_top_cyl[i] /= (double)n_edges_top_loop;
  for (int ii = 0; ii < n_edges_bot_loop; ++ii)  {
    const int ist_bot = botEdgeVec[ii]/3;
    const int i_bot = botEdgeVec[ii]-ist_bot*3;
    FOR_I3 x_bot_cyl[i] += xsp[spost[ist_bot][i_bot]][i];
  }
  FOR_I3 x_bot_cyl[i] /= (double)n_edges_bot_loop;
  double dx_top_to_bot_cyl[3] = DIFF(x_bot_cyl,x_top_cyl);
  cout << "orientation: " << COUT_VEC(dx_top_to_bot_cyl)  << endl;

  for (int isp = 0; isp < nsp; ++isp) {
    const double hdx = DOT_PRODUCT(sp_buf3[isp],dx_top_to_bot_cyl);
    const double dxdx = DOT_PRODUCT(dx_top_to_bot_cyl,dx_top_to_bot_cyl);
    FOR_J3 sp_buf3[isp][j] -= hdx/dxdx*dx_top_to_bot_cyl[j];
  }

  // now rescale cylinder...

  for (int toz = stosz_i[cyl_id]; toz != stosz_i[cyl_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == cyl_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] >= 0) {
        const double h_over_h0 = 1.0+SGN(dx_bot_to_top[1])*xsp[isp][0]*tan(phi)/H0;
        const int ist_top = topEdgeVec[sp_flag[isp]]/3;
        const int i_top = topEdgeVec[sp_flag[isp]]-ist_top*3;
        const double h = DIST(xsp[spost[ist_top][i_top]],xsp[isp]);
        FOR_J3 xsp[isp][j] = xsp[spost[ist_top][i_top]][j]+h_over_h0*h*dx_top_to_bot_cyl[j]/MAG(dx_top_to_bot_cyl)+sp_buf3[isp][j];
        sp_flag[isp] = -1;
      }
    }
  }

  // back to lab frame
  for (int isp = 0; isp < nsp; ++isp) sp_flag[isp] = -1;
  for (int toz = stosz_i[cyl_id]; toz != stosz_i[cyl_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == cyl_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == -1) {
        sp_flag[isp] = 0;
        // const double x0[2] = {xsp[isp][0],xsp[isp][1]};
//        xsp[isp][0] = x0[0]*cos_angle+x0[1]*sin_angle;
//        xsp[isp][1] =-x0[0]*sin_angle+x0[1]*cos_angle;
        // phi
        const double x1[2] = {xsp[isp][0],xsp[isp][1]};
        xsp[isp][0] = x1[0]*cos(phi)-x1[1]*sin(phi);
        xsp[isp][1] = x1[0]*sin(phi)+x1[1]*cos(phi);
        FOR_J3 xsp[isp][j] += x_bot[j];
      }
    }
  }
  for (int toz = stosz_i[top_id]; toz != stosz_i[top_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == top_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == -1) {
        sp_flag[isp] = 0;
        // const double x0[2] = {xsp[isp][0],xsp[isp][1]};
//        xsp[isp][0] = x0[0]*cos_angle+x0[1]*sin_angle;
//        xsp[isp][1] =-x0[0]*sin_angle+x0[1]*cos_angle;
        // phi
        const double x1[2] = {xsp[isp][0],xsp[isp][1]};
        xsp[isp][0] = x1[0]*cos(phi)-x1[1]*sin(phi);
        xsp[isp][1] = x1[0]*sin(phi)+x1[1]*cos(phi);
        FOR_J3 xsp[isp][j] += x_bot[j];
      }
    }
  }

  // now we just need to redistribute base to prevent folding...

  for (int isp = 0; isp < nsp; ++isp) sp_flag[isp] = 0; // imovable
  for (int toz = stosz_i[base_id]; toz != stosz_i[base_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == base_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == 0) sp_flag[isp] = -1; // interior
    }
  }
  for (int ii = 0; ii < n_edges_bot_loop; ++ii)  {
    const int ist_bot = botEdgeVec[ii]/3;
    const int i_bot = botEdgeVec[ii]-ist_bot*3;
    const int isp = spost[ist_bot][i_bot];
    if (sp_flag[isp] == -1) {
      sp_flag[isp] = 1;
    }
  }

  // build sposp...
  ensureSposp();

  // solve Laplace's equation using relaxation method...
  for (int isp = 0; isp < nsp; ++isp) FOR_I3 sp_buf3[isp][i] = 0.0;
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] != 1)
      FOR_I3 sp_buf3[isp][i] = 0.0;
    else
      FOR_I3 sp_buf3[isp][i] = xsp[isp][i]-xsp0[isp][i];
  }
  for (int ii = 0; ii < nrelax; ++ii) {
    for (int toz = stosz_i[base_id]; toz != stosz_i[base_id+1]; ++toz) {
      const int ist = stosz_v[toz];
      assert(szost[ist] == base_id);
      FOR_I3 {
        const int isp = spost[ist][i];
        if (sp_flag[isp] == -ii-1) {
          sp_flag[isp] -= 1;
          FOR_J3 sp_buf3[isp][j] = 0.0;
          for (int pop = sposp_i[isp]+1; pop != sposp_i[isp+1]; ++pop) {
            const int isp2 = sposp_v[pop];
            FOR_J3 sp_buf3[isp][j] += sp_buf3[isp2][j];
          }
          if ((sposp_i[isp+1]-(sposp_i[isp]+1)) > 1) {
            FOR_J3 sp_buf3[isp][j] /= sposp_i[isp+1]-(sposp_i[isp]+1);
          }
        }
      }
    }
  }
  cout << "write out laplacian_solve.plt" << endl;

 // writeSpDataToTecplot("laplacian_solve.plt",sp_buf3);

  // apply tranformation...
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] < 0) FOR_I3 xsp[isp][i] += sp_buf3[isp][i];
  }

  // check the morph...
  bool bad_top = checkMorphedZone(top_id);
  bool bad_bot = checkMorphedZone(base_id);
  bool bad_cyl = checkMorphedZone(cyl_id);

  if (bad_top || bad_bot || bad_cyl) {
    // revert back so you can redo the morph
    for (int isp = 0; isp < nsp; ++isp) FOR_I3 xsp[isp][i] = xsp0[isp][i];
  }
  else {
    // surface will change, recompute geoms
    clearSubzoneData();
    clearZoneData();
    b_centroid = false;
    b_rmax = false;
    b_boundingBox = false;
    clearNonManifoldData();

    // if we had periodic data, we need to reset it here
    if (pbi) {
      clearPeriodicity();
      WUI(WARN,"Bad news, friend. You need to reset periodicity.");
    }
  }
  delete[] xsp0;

}

void SimpleSurface::morphRotateSelectedSubzones(const vector<int>& subzonesVec,const double _axis[3], const double point[3],const double angle_deg) {

  // need stosz and teost (implicit)
  ensureStoszSzosz();

  // untouched subzones = 0
  // selected  subzones = 1
  // neighbor  subzones = 2

  sz_flag.resize(nsz);
  for (int isz = 0; isz < nsz; ++isz) 
    sz_flag[isz] = 0;
  for (int ii = 0, lim = subzonesVec.size(); ii < lim; ++ii) {
    const int isz = subzonesVec[ii];
    if ((isz >= 0)&&(isz < nsz)) {
      sz_flag[isz] = 1;
      for (int sos = szosz_i[isz]; sos < szosz_i[isz+1]; ++sos) {
        const int isz_nbr = szosz_v[sos]; assert((isz_nbr >= 0)&&(isz_nbr < nsz));
        if (sz_flag[isz_nbr] != 1)
          sz_flag[isz_nbr] = 2;
      }
    }
    else {
      cout << " > subzone " << isz << " outside of [0," << nsz << "). Skipping..." << endl;
    }
  }

  vector<int> nbrSubzoneVec;
  for (int isz = 0; isz < nsz; ++isz) {
    if (sz_flag[isz] == 2) {
      nbrSubzoneVec.push_back(isz);
    }
  }

  // untouched points = 0
  // selected  points = 1
  // neighbor  points = 2
  // border    points = 3
  
  sp_flag.setLength(nsp);
  sp_flag.setAll(2);
  FOR_IST {
    const int isz = szost[ist];
    if (sz_flag[isz] == 1) 
      FOR_I3 sp_flag[spost[ist][i]] |= sz_flag[isz];
    else if (sz_flag[isz] == 0)
      FOR_I3 sp_flag[spost[ist][i]] = 0;
  }

  COUT2(" > rotation axis (pt,normal) and angle: " << COUT_VEC(point) << " " << COUT_VEC(_axis) << " " << angle_deg << " degrees");

  // surface will change, recompute geoms
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

  double (*dx)[3] = new double[nsp][3];
  FOR_ISP {
    if (sp_flag[isp] & 1) {
      double _xsp[3];
      FOR_I3 _xsp[i] = xsp[isp][i] - point[i];
      const double temp[3] = MATRIX_PRODUCT(R,_xsp);
      FOR_I3 dx[isp][i] = temp[i] + point[i] - xsp[isp][i];
    }
    else {
      FOR_I3 dx[isp][i] = 0.0;
    }
  }

  // build sposp...
  ensureSposp();

  // solve Laplace's equation using relaxation method...
  for (int ii = 0; ii < 200; ++ii) {
    for (int jj = 0, lim = nbrSubzoneVec.size(); jj < lim; ++jj) {
      const int isz = nbrSubzoneVec[jj];
      for (int toz = stosz_i[isz]; toz != stosz_i[isz+1]; ++toz) {
        const int ist = stosz_v[toz];
        assert(szost[ist] == isz);
        FOR_I3 {
          const int isp = spost[ist][i];
          if (sp_flag[isp] == 2) {
            FOR_J3 dx[isp][j] = 0.0;
            for (int pop = sposp_i[isp]+1; pop != sposp_i[isp+1]; ++pop) {
              const int isp2 = sposp_v[pop];
              FOR_J3 dx[isp][j] += dx[isp2][j];
            }
            if ((sposp_i[isp+1]-(sposp_i[isp]+1)) > 1) {
              FOR_J3 dx[isp][j] /= sposp_i[isp+1]-(sposp_i[isp]+1);
            }
          }
        }
      }
    }
  }

  //cout << "write out laplacian_solve.plt" << endl;
  //writeSpDataToTecplot("laplacian_solve.plt",dx);

  // apply tranformation...
  for (int isp = 0; isp < nsp; ++isp) 
    if (sp_flag[isp] > 0) FOR_I3 xsp[isp][i] += dx[isp][i];
  delete[] dx;

  // surface will change, recompute geoms
  clearDynamicEdgeGroups();
  clearOpenEdgeGroups();
  clearSubzoneData();
  clearZoneData();
  b_centroid = false;
  b_rmax = false;
  b_boundingBox = false;
  clearNonManifoldData();

  // if we had periodic data, we need to reset it here
  if (pbi) {
    clearPeriodicity();
    WUI(WARN,"Bad news, friend. You need to reset periodicity.");
  }

}

void SimpleSurface::morphCorner(const int top, const int side, const double dx, const double dl,
    const bool use_bump, const double dw) {

  // need stosz, szosz and teost (implicit)
  ensureStoszSzosz();

  assert((top >= 0)&&(top < nsz));
  assert((side >= 0)&&(side < nsz));

  // make sure side and top are connected in the graph...

  {
    int sos;
    for (sos = szosz_i[top]; sos != szosz_i[top+1]; ++sos) {
      if (szosz_v[sos] == side) break;
    }
    if (sos == szosz_i[top+1]) {
      cout << "Warning: TOP and SIDE are not connected subsurfaces. Skipping." << endl;
      return;
    }
  }

  // and do the morph...

  // if (sp_flag == NULL) sp_flag = new int[nsp];
  sp_flag.setLength(nsp);

  double x_top[3] = { 0.0, 0.0, 0.0 };
  double n_top[3] = { 0.0, 0.0, 0.0 };
  double area_top = 0.0;

  assert((top >= 0)&&(top < nsz));
  for (int toz = stosz_i[top]; toz != stosz_i[top+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == top);
    FOR_I3 sp_flag[spost[ist][i]] = 1;
    // these tris contribute to the top normal...
    const double * const x0 = xsp[spost[ist][0]];
    const double * const x1 = xsp[spost[ist][1]];
    const double * const x2 = xsp[spost[ist][2]];
    const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
    FOR_I3 n_top[i] += this_n[i];
    const double this_area = MAG(this_n);
    FOR_I3 x_top[i] += this_area*(x0[i]+x1[i]+x2[i]);
    area_top += this_area; // actually 2x the area
  }

  FOR_I3 x_top[i] /= area_top*3.0;
  const double mag_top = MAG(n_top);
  FOR_I3 n_top[i] /= mag_top;

  // check planarity of top...

  double dn_max = 0.0;
  for (int toz = stosz_i[top]; toz != stosz_i[top+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == top);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == 1) {
        sp_flag[isp] = 0;
        // this is a top node. Check the distance...
        const double dx_top[3] = DIFF(x_top,xsp[isp]);
        const double this_dn = DOT_PRODUCT(dx_top,n_top);
        dn_max = max(dn_max,fabs(this_dn));
      }
    }
  }

  cout << "[ACTION] MORPH_CORNER x_top: " << COUT_VEC(x_top) << " n_top: " << COUT_VEC(n_top) << " planarity: " << dn_max << endl;

  // now the side -- this is what gets morphed...

  if (sp_buf3 == NULL) sp_buf3 = new double[nsp][3];
  double (*sp_normal)[3] = sp_buf3; // for clarity

  //for (int isp = 0; isp < nsp; ++isp)
  //  FOR_I3 sp_normal[isp][i] = 0.f;

  assert((side >= 0)&&(side < nsz));
  sp_flag.setAll(-1);
  for (int toz = stosz_i[side]; toz != stosz_i[side+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == side);
    FOR_I3 {
      const int isp = spost[ist][i];
      const double dx_top[3] = DIFF(xsp[isp],x_top);
      const double this_dl = DOT_PRODUCT(dx_top,n_top);
      if (this_dl < dl) {
        sp_flag[isp] = 1;
        sp_normal[isp][0] = 0.f;
        sp_normal[isp][1] = 0.f;
        sp_normal[isp][2] = 0.f;
      }
      else {
        sp_flag[isp] = 0;
      }
    }
  }

  // now build the normal in sp_flag == 1 surface points...

  for (int toz = stosz_i[side]; toz != stosz_i[side+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == side);
    if ((sp_flag[spost[ist][0]] == 1)||(sp_flag[spost[ist][1]] == 1)||(sp_flag[spost[ist][2]] == 1)) {
      // these tris contribute to the individual side normals...
      const double * const x0 = xsp[spost[ist][0]];
      const double * const x1 = xsp[spost[ist][1]];
      const double * const x2 = xsp[spost[ist][2]];
      const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
      FOR_I3 {
        if (sp_flag[spost[ist][i]] == 1)
          FOR_J3 sp_normal[spost[ist][i]][j] += this_n[j];
      }
    }
  }

  // now go through and morph...

  int nsp_morph = 0;
  double (*xsp0)[3] = new double[nsp][3];
  for (int isp = 0; isp < nsp; ++isp) FOR_I3 xsp0[isp][i] = xsp[isp][i];
  for (int toz = stosz_i[side]; toz != stosz_i[side+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == side);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == 1) {
        sp_flag[isp] = 0;
        // this is a side node: move in based on distance from top...
        const double dx_top[3] = DIFF(xsp[isp],x_top);
        const double this_dl = DOT_PRODUCT(dx_top,n_top);
        assert(this_dl < dl);
        // this guy gets moved in the sp_normal direction...
        // first remove any projection of sp_normal on n_top...
        const double dp = DOT_PRODUCT(sp_normal[isp],n_top);
        FOR_J3 sp_normal[isp][j] -= dp*n_top[j];
        const double mag = MAG(sp_normal[isp]);
        if (use_bump) {
          if (this_dl/dl < 1.0) {
            FOR_J3 xsp[isp][j] += dx*exp(1.0-1.0/(1.0-this_dl*this_dl/dl/dl))*sp_normal[isp][j]/mag;
          }
        }
        else {
          FOR_J3 xsp[isp][j] += dx*(1.0-this_dl/dl)*(1.0-this_dl/dl)*sp_normal[isp][j]/mag;
        }
        ++nsp_morph;
        int ist_nbr, i_nbr, orient_nbr;
        if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
          if (szost[ist_nbr] == top) {
            sp_flag[isp] = 2;
          }
        }
      }
    }
  }

  // morph internal points of top subsurface to help prevent folding

  for (int toz = stosz_i[top]; toz != stosz_i[top+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == top);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == -1) {
        sp_flag[isp] = 1;
      }
    }
  }
  for (int toz = stosz_i[top]; toz != stosz_i[top+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == top);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == 1) {
        sp_flag[isp] = 0;

        // find closest distance to top border
        double min_dist = 1.0E+20;
        double this_w[3];
        for (int toz2 = stosz_i[top]; toz2 != stosz_i[top+1]; ++toz2) {
          const int ist2 = stosz_v[toz2];
          assert(szost[ist2] == top);
          FOR_J3 {
            const int isp2 = spost[ist2][j];
            if (sp_flag[isp2] == 2) {
              if (DIST(xsp[isp],xsp[isp2]) <= min_dist) {
                min_dist = DIST(xsp[isp],xsp[isp2]);
                FOR_K3 this_w[k] = xsp[isp][k]-xsp[isp2][k];
              }
            }
          }
        }

        // push points inward using bump
        const double sp_dir[3] = DIFF(x_top,xsp[isp]);
        const double mag = MAG(sp_dir);
        const double this_dw = DOT_PRODUCT(this_w,sp_dir)/mag;
        // always use bump function to locally redistribute points
        //if (use_bump) {
          if (this_dw/dw < 1.0) {
            FOR_J3 xsp[isp][j] += dx*exp(1.0-1.0/(1.0-this_dw*this_dw/dw/dw))*sp_dir[j]/mag;
          }
        //}
        ++nsp_morph;
      }
    }
  }

  cout << "[ACTION] MORPH_CORNER nsp_morph: " << nsp_morph << endl;

  // check the orientation
  bool bad_top = checkMorphedZone(top);
  bool bad_side = checkMorphedZone(side);
  if (bad_top || bad_side) {
    // revert back so you can redo the morph
    for (int isp = 0; isp < nsp; ++isp) FOR_I3 xsp[isp][i] = xsp0[isp][i];
  }
  else {
    // surface will change, recompute geoms
    clearSubzoneData();
    clearZoneData();
    b_centroid = false;
    b_rmax = false;
    b_boundingBox = false;
    clearNonManifoldData();

    // if we had periodic data, we need to reset it here
    if (pbi) {
      WUI(WARN,"Bad news, friend. You need to reset periodicity.");
    }
  }
  delete[] xsp0;

}

void SimpleSurface::morphStretch(const int top, const double dn, const bool use_bump, const double dl) {

  // need stosz and teost (implicit)
  ensureStoszSzosz();

  assert((top >= 0)&&(top < nsz));

  sp_flag.setLength(nsp);

  double x_top[3] = { 0.0, 0.0, 0.0 };
  double n_top[3] = { 0.0, 0.0, 0.0 };
  double area_top = 0.0;

  for (int toz = stosz_i[top]; toz != stosz_i[top+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == top);
    FOR_I3 sp_flag[spost[ist][i]] = 1;
    // these tris contribute to the top normal...
    const double * const x0 = xsp[spost[ist][0]];
    const double * const x1 = xsp[spost[ist][1]];
    const double * const x2 = xsp[spost[ist][2]];
    const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
    FOR_I3 n_top[i] += this_n[i];
    const double this_area = MAG(this_n);
    FOR_I3 x_top[i] += this_area*(x0[i]+x1[i]+x2[i]);
    area_top += this_area; // actually 2x the area
  }

  FOR_I3 x_top[i] /= area_top*3.0;
  const double mag_top = MAG(n_top);
  FOR_I3 n_top[i] /= mag_top;

  // check planarity of top...

  stack<int> sideTriStack;

  double dn_max = 0.0;
  for (int toz = stosz_i[top]; toz != stosz_i[top+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == top);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == 1) {
        sp_flag[isp] = 0;
        // this is a top node. Check the distance...
        const double dx_top[3] = DIFF(x_top,xsp[isp]);
        const double this_dn = DOT_PRODUCT(dx_top,n_top);
        dn_max = max(dn_max,fabs(this_dn));
      }
    }
    // also start to build the set of tris connected to this surface...
    FOR_I3 {
      int ist_nbr,i_nbr,orient_nbr;
      if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
        if ((szost[ist_nbr] >= 0)&&(szost[ist_nbr] != top)) {
          szost[ist_nbr] = -szost[ist_nbr]-1;
          sideTriStack.push(ist_nbr);
        }
      }
    }

  }
  for (int toz = stosz_i[top]; toz != stosz_i[top+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == top);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == 0) {
        sp_flag[isp] = 1;
      }
    }
  }


  cout << "[ACTION] MORPH_STRETCH x_top: " << COUT_VEC(x_top) << " n_top: " << COUT_VEC(n_top) << " planarity: " << dn_max << endl;

  vector<int> sideTriVec;

  while (!sideTriStack.empty()) {
    // off the stack and into the vec...
    const int ist = sideTriStack.top(); sideTriStack.pop();
    assert(szost[ist] < 0);
    sideTriVec.push_back(ist);
    // figure out which nodes are beyond 2*dn...
    double this_dn[3];
    FOR_I3 {
      // this is a top node. Check the distance...
      const int isp = spost[ist][i];
      sp_flag[isp] = 1;
      const double dx_top[3] = DIFF(x_top,xsp[isp]);
      this_dn[i] = fabs(DOT_PRODUCT(dx_top,n_top));
    }
    FOR_I3 {
      // only grab a tri beyond an edge if either one of its nodes is within the dn range...
      int ist_nbr,i_nbr,orient_nbr;
      if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
        if ( (this_dn[i] < dn+dl)||(this_dn[(i+1)%3] < dn+dl) ) {
          if ((szost[ist_nbr] >= 0)&&(szost[ist_nbr] != top)) {
            szost[ist_nbr] = -szost[ist_nbr]-1;
            sideTriStack.push(ist_nbr);
          }
        }
      }
    }
  }

  cout << "sideTriVec.size(): " << sideTriVec.size() << endl;

  // now everything went through the stack and into the vec, so we can
  // morph...

  set<int> subzoneSet;
  double (*xsp0)[3] = new double[nsp][3];
  for (int isp = 0; isp < nsp; ++isp) FOR_I3 xsp0[isp][i] = xsp[isp][i];
  for (int ii = 0,ii_end=sideTriVec.size(); ii < ii_end; ++ii) {
    const int ist = sideTriVec[ii];
    assert(szost[ist] < 0);
    szost[ist] = -szost[ist]-1;
    subzoneSet.insert(szost[ist]);
    // now morph...
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == 1) {
        sp_flag[isp] = 0;
        const double dx_top[3] = DIFF(x_top,xsp[isp]);
        const double this_dn = fabs(DOT_PRODUCT(dx_top,n_top));
        if (use_bump) {
          if (this_dn/dl < 1.0) {
            FOR_J3 xsp[isp][j] += dn*exp(1.0-1.0/(1.0-this_dn*this_dn/dl/dl))*n_top[j];
          }
        }
        else {
          if (this_dn/dl < 1.0) {
            FOR_J3 xsp[isp][j] += dn*(1.0-this_dn/dl)*n_top[j];
          }
        }
      }
    }
  }

  // stretch rest of top
  for (int toz = stosz_i[top]; toz != stosz_i[top+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == top);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == 1) {
        sp_flag[isp] = 0;
        FOR_J3 xsp[isp][j] += dn*n_top[j];
      }
    }
  }

  // check the morph
  bool bad_top = checkMorphedZone(top);
  bool bad_sz;
  for (set<int>::iterator iter = subzoneSet.begin(); iter != subzoneSet.end(); ++iter) {
    bad_sz = checkMorphedZone(*iter);
    if (bad_sz) break;
  }

  if (bad_sz || bad_top) {
    // revert back so you can redo the morph
    for (int isp = 0; isp < nsp; ++isp) FOR_I3 xsp[isp][i] = xsp0[isp][i];
  }
  else {
    // surface will change, recompute geoms
    clearSubzoneData();
    clearZoneData();
    b_centroid = false;
    b_rmax = false;
    b_boundingBox = false;
    clearNonManifoldData();

    // if we had periodic data, we need to reset it here
    if (pbi) {
      clearPeriodicity();
      WUI(WARN,"Bad news, friend. You need to reset periodicity.");
    }
  }
  delete[] xsp0;

}

void SimpleSurface::morphCylinder(const int cyl_id, const int top_id, const int bot_id, const double L,const double theta0, const double a, const double b, const int n) {
  UNUSED(L);
  UNUSED(theta0);
  // need szosz
  ensureStoszSzosz();

  assert((cyl_id >= 0)&&(cyl_id < nsz));
  assert((top_id >= 0)&&(top_id < nsz));
  assert((bot_id >= 0)&&(bot_id < nsz));

  // make sure cylinder/top and cylinder/bottom are connected in the graph...
  {
    int sos;
    for (sos = szosz_i[cyl_id]; sos != szosz_i[cyl_id+1]; ++sos) {
      if (szosz_v[sos] == top_id) break;
    }
    if (sos == szosz_i[cyl_id+1]) {
      cout << "Warning: CYLINDER AND TOP are not connected subsurfaces. Skipping." << endl;
      return;
    }
    for (sos = szosz_i[cyl_id]; sos != szosz_i[cyl_id+1]; ++sos) {
      if (szosz_v[sos] == bot_id) break;
    }
    if (sos == szosz_i[cyl_id+1]) {
      cout << "Warning: CYLINDER AND BOTTOM are not connected subsurfaces. Skipping." << endl;
      return;
    }
  }

  // will use flag for marking subzone boundaries and to demark sides of square prism
  sp_flag.setLength(nsp);
  sp_flag.setAll(-1);
  vector<int> topEdgeVec;
  vector<int> botEdgeVec;

  // start by building vectors of top and bottom edges shared with cyl
  for (int toz = stosz_i[cyl_id]; toz != stosz_i[cyl_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == cyl_id);
    FOR_I3 {
      int ist_nbr,i_nbr,orient_nbr;
      if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
        if (szost[ist_nbr] == top_id) {
          topEdgeVec.push_back(ist*3+i);
          // put our edge index in one of our nodes
          if (sp_flag[spost[ist][i]] == -1) {
            sp_flag[spost[ist][i]] = ist*3+i;
          }
          else {
            // more than 2 open edges must touch this node
            cout << "Warning: top edge loop algorithm problem. Skipping." << endl;
            return;
          }
        }
        else if (szost[ist_nbr] == bot_id) {
          botEdgeVec.push_back(ist*3+i);
          // put our edge index in one of our nodes
          if (sp_flag[spost[ist][i]] == -1) {
            sp_flag[spost[ist][i]] = ist*3+i;
          }
          else {
            // more than 2 bottom edges must touch this node (or top/bot loops touch)
            cout << "Warning: bottom edge loop algorithm problem. Skipping." << endl;
            return;
          }
        }
      }
    }
  }

  // Check that the top and bottom loops are closed by walking them. also get the loop centers.
  int n_edges_top_loop = 0;
  double x_top[3] = { 0.0, 0.0, 0.0 };
  vector<int> auxEdgeVec(topEdgeVec.size());
  for (int ii = 0,ii_end=topEdgeVec.size(); ii < ii_end; ++ii) {
    int ied = topEdgeVec[ii];
    int ist = ied/3;
    int i = ied-ist*3;
    int ist_nbr,i_nbr,orient_nbr;
    if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
      if (szost[ist_nbr] == top_id) {
        FOR_J3 x_top[j] += xsp[spost[ist][i]][j];
        szost[ist_nbr] = -top_id-1;
        auxEdgeVec[n_edges_top_loop] = ied;
        n_edges_top_loop++;
        while (sp_flag[spost[ist][(i+1)%3]] != -1) {
          ied = sp_flag[spost[ist][(i+1)%3]];
          ist = ied/3;
          i = ied-ist*3;
          if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
            if (szost[ist_nbr] == top_id) {
              FOR_J3 x_top[j] += xsp[spost[ist][i]][j];
              szost[ist_nbr] = -top_id-1;
              auxEdgeVec[n_edges_top_loop] = ied;
              n_edges_top_loop++;
            }
            else {
              assert(szost[ist_nbr] == -top_id-1);
              break;
            }
          }
        }
      }
    }
  }
  assert(n_edges_top_loop == int(topEdgeVec.size()));
  FOR_I3 x_top[i] /= n_edges_top_loop;
  for (int ii = 0; ii < n_edges_top_loop; ++ii) {
    topEdgeVec[ii] = auxEdgeVec[ii];
    const int ied = topEdgeVec[ii];
    const int ist = ied/3;
    const int i = ied-ist*3;
    int ist_nbr,i_nbr,orient_nbr;
    if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
      assert(szost[ist_nbr] == -top_id-1);
      szost[ist_nbr] = top_id;
    }
  }

  // now repeat for bottom loop
  int n_edges_bot_loop = 0;
  double x_bot[3] = { 0.0, 0.0, 0.0 };
  auxEdgeVec.resize(botEdgeVec.size());
  for (int ii = 0,ii_end=botEdgeVec.size(); ii < ii_end; ++ii) {
    int ied = botEdgeVec[ii];
    int ist = ied/3;
    int i = ied-ist*3;
    int ist_nbr,i_nbr,orient_nbr;
    if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
      if (szost[ist_nbr] == bot_id) {
        FOR_J3 x_bot[j] += xsp[spost[ist][i]][j];
        szost[ist_nbr] = -bot_id-1;
        auxEdgeVec[n_edges_bot_loop] = ied;
        n_edges_bot_loop++;
        while (sp_flag[spost[ist][(i+1)%3]] != -1) {
          ied = sp_flag[spost[ist][(i+1)%3]];
          ist = ied/3;
          i = ied-ist*3;
          if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
            if (szost[ist_nbr] == bot_id) {
              FOR_J3 x_bot[j] += xsp[spost[ist][i]][j];
              szost[ist_nbr] = -bot_id-1;
              auxEdgeVec[n_edges_bot_loop] = ied;
              n_edges_bot_loop++;
            }
            else {
              assert(szost[ist_nbr] == -bot_id-1);
              break;
            }
          }
        }
      }
    }
  }

  assert(n_edges_bot_loop == int(botEdgeVec.size()));
  FOR_I3 x_bot[i] /= n_edges_bot_loop;
  for (int ii = 0; ii < n_edges_bot_loop; ++ii) {
    botEdgeVec[ii] = auxEdgeVec[ii];
    const int ied = botEdgeVec[ii];
    const int ist = ied/3;
    const int i = ied-ist*3;
    int ist_nbr,i_nbr,orient_nbr;
    if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
      assert(szost[ist_nbr] == -bot_id-1);
      szost[ist_nbr] = bot_id;
    }
  }

  // calculate the central axis
  double dx_bot_to_top[3] = DIFF(x_top,x_bot);
  const double dx_bot_to_top_mag2 = DOT_PRODUCT(dx_bot_to_top,dx_bot_to_top);

  // calculate the min and max extensions between top and bot
  double min_extent = +1E+20;
  double max_extent = -1E+20;
  for (int itop = 0; itop < n_edges_top_loop; ++itop) {
    const int ist = topEdgeVec[itop]/3;
    const int i = topEdgeVec[itop]-ist*3;
    // const double * const x1 = xsp[spost[ist][i]];
    for (int ibot = 0; ibot < n_edges_bot_loop; ++ibot) {
      const int ist_ = botEdgeVec[ibot]/3;
      const int i_ = botEdgeVec[ibot]-ist_*3;
      const double dx[3] = DIFF(xsp[spost[ist][i]],xsp[spost[ist_][i_]]);
      const double extent = DOT_PRODUCT(dx,dx_bot_to_top)/sqrt(dx_bot_to_top_mag2);
      min_extent = min(min_extent,extent);
      max_extent = max(max_extent,extent);
    }
  }

  // now lets walk them again using ordered top/botEdgeVec to calculate the loop normals and areas
  double n_top[3] = { 0.0, 0.0, 0.0 }; // 2x area weighted
  double area_top_2 = 0.0; // 2x the area
  for (int ii = 0; ii < n_edges_top_loop; ++ii) {
    int iip1 = ii+1;
    if (iip1 == n_edges_top_loop) iip1 = 0;
    const int ist = topEdgeVec[ii]/3;
    const int i = topEdgeVec[ii]-ist*3;
    const int istp1 = topEdgeVec[iip1]/3;
    const int ip1 = topEdgeVec[iip1]-istp1*3;
    const double * const x0 = x_top;
    const double * const x1 = xsp[spost[ist][i]];
    const double * const x2 = xsp[spost[istp1][ip1]];
    const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
    FOR_I3 n_top[i] += this_n[i];
    const double this_area = MAG(this_n);
    area_top_2 += this_area;
  }
  // make sure walks and normals are pointing out of cylinder
  if (DOT_PRODUCT(n_top,dx_bot_to_top) < 0.0) {
    reverse(topEdgeVec.begin(),topEdgeVec.end());
    FOR_I3 n_top[i] *= -1;
  }

  // now repeat for bottom loop
  double n_bot[3] = { 0.0, 0.0, 0.0 }; // 2x area weighted
  double area_bot_2 = 0.0; // 2x the area
  for (int ii = 0; ii < n_edges_bot_loop; ++ii) {
    int iip1 = ii+1;
    if (iip1 == n_edges_bot_loop) iip1 = 0;
    const int ist = botEdgeVec[ii]/3;
    const int i = botEdgeVec[ii]-ist*3;
    const int istp1 = botEdgeVec[iip1]/3;
    const int ip1 = botEdgeVec[iip1]-istp1*3;
    const double * const x0 = x_bot;
    const double * const x1 = xsp[spost[ist][i]];
    const double * const x2 = xsp[spost[istp1][ip1]];
    const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
    FOR_I3 n_bot[i] += this_n[i];
    const double this_area = MAG(this_n);
    area_bot_2 += this_area;
  }
  // make sure walks and normals are pointing out of cylinder
  if (DOT_PRODUCT(n_bot,dx_bot_to_top) > 0.0) {
    reverse(botEdgeVec.begin(),botEdgeVec.end());
    FOR_I3 n_bot[i] *= -1;
  }

  // check how close top_cyl surfaces approximate a cylinder by getting min/max/avg r_mag
  double r_max_2 = -1E+20;
  double r_min_2 = +1E+20;
  double r_avg_2 = 0.0;
  for (int toz = stosz_i[cyl_id]; toz != stosz_i[cyl_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == cyl_id);
    FOR_I3 {
      const double * const this_x =xsp[spost[ist][i]];
      const double dx_bot[3] = DIFF(this_x,x_bot);
      const double dx_top[3] = DIFF(this_x,x_top);
      double r_vec[3];
      // calc r using the closest of x_bot and x_top
      if (DOT_PRODUCT(dx_bot,dx_bot) <= DOT_PRODUCT(dx_top,dx_top)) {
        const double dx_frac = DOT_PRODUCT(dx_bot,dx_bot_to_top)/dx_bot_to_top_mag2;
        FOR_J3 r_vec[j] = dx_bot[j] - dx_frac * dx_bot_to_top[j];
      }
      else {
        const double dx_frac = DOT_PRODUCT(dx_top,dx_bot_to_top)/dx_bot_to_top_mag2;
        FOR_J3 r_vec[j] = dx_top[j] - dx_frac * dx_bot_to_top[j];
      }
      const double r_mag_2 = DOT_PRODUCT(r_vec,r_vec);
      r_max_2 = max(r_max_2,r_mag_2);
      r_min_2 = min(r_min_2,r_mag_2);
      r_avg_2 += r_mag_2;
    }
  }
  r_avg_2 /= (3*(stosz_i[cyl_id+1]-stosz_i[cyl_id]));
  double r_max = sqrt(r_max_2);
  double r_min = sqrt(r_min_2);
  double r_avg = sqrt(r_avg_2);

  // report the cylinder characterization
  cout << "top: x_top[3] = " << COUT_VEC(x_top) << ", 2 A n_top[3] = " << COUT_VEC(n_top) << ", A_top = " << 0.5*area_top_2 << endl;
  cout << "bot: x_bot[3] = " << COUT_VEC(x_bot) << ", 2 A n_bot[3] = " << COUT_VEC(n_bot) << ", A_bot = " << 0.5*area_bot_2 << endl;
  cout << "central axis: x_top[3]-x_bot[3] = " << COUT_VEC(dx_bot_to_top) << endl;
  cout << "radii: mean(r) = " << r_avg << ", min(r) = " << r_min << ", max(r) = " << r_max << endl;
  cout << "extent: h1 = " << min_extent << ", h2 = " << max_extent << endl;
  cout << "volume: 1/2*(h1+h2)*pi*r^2 = " << 0.5*(min_extent+max_extent)*M_PI*r_avg*r_avg << endl;

  // make sure our spread is less than 5%
  assert(fabs(r_max-r_min)/r_avg < 0.05);

  // leave if they gave circle of the same radius
  const double EPS = 2.220446049250313E-16;
  if ( fabs(a-r_avg) <= EPS && fabs(b-r_avg) <= EPS && n == 2) return;

  // right now this doesn't do anything after checking we have a cylinder
  // XXXX add morphing here
  double (*xsp0)[3] = new double[nsp][3];
  for (int isp = 0; isp < nsp; ++isp) FOR_I3 xsp0[isp][i] = xsp[isp][i];

  // check the morph...
  bool bad_top = checkMorphedZone(top_id);
  bool bad_bot = checkMorphedZone(bot_id);
  bool bad_cyl = checkMorphedZone(cyl_id);

  if (bad_top || bad_bot || bad_cyl) {
    // revert back so you can redo the morph
    for (int isp = 0; isp < nsp; ++isp) FOR_I3 xsp[isp][i] = xsp0[isp][i];
  }
  else {
    // surface will change, recompute geoms
    clearSubzoneData();
    clearZoneData();
    b_centroid = false;
    b_rmax = false;
    b_boundingBox = false;
    clearNonManifoldData();

    // if we had periodic data, we need to reset it here
    if (pbi) {
      clearPeriodicity();
      WUI(WARN,"Bad news, friend. You need to reset periodicity.");
    }
  }
  delete[] xsp0;

}

void SimpleSurface::morphCylinderZ(const int cyl_id, const int top_id, const int bot_id, const double x0, const double y0,
        const double r0, const double L, const double theta0, const double a, const double b, const int n) {

  // need stosz
  ensureStoszSzosz();

  assert((cyl_id >= 0)&&(cyl_id < nsz));
  assert((top_id >= 0)&&(top_id < nsz));
  assert((bot_id >= 0)&&(bot_id < nsz));

  // leave if they gave circle of the same radius
  const double EPS = 2.220446049250313E-16;
  if ( fabs(a-r0) <= EPS && fabs(b-r0) <= EPS && n == 2) return;

  sp_flag.setLength(nsp); sp_flag.setAll(1);
  double (*xsp0)[3] = new double[nsp][3];
  for (int isp = 0; isp < nsp; ++isp) FOR_I3 xsp0[isp][i] = xsp[isp][i];
  for (int toz = stosz_i[cyl_id]; toz != stosz_i[cyl_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == cyl_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == 1) {
        sp_flag[isp] = 0; // operate once on point
        const double my_theta = atan2(xsp[isp][1]-y0,xsp[isp][0]-x0);
        xsp[isp][0] = x0 + pow(fabs(cos(my_theta)),2.0/double(n))*a*SGN(cos(my_theta));
        xsp[isp][1] = y0 + pow(fabs(sin(my_theta)),2.0/double(n))*b*SGN(sin(my_theta));
        const double my_dx = xsp[isp][0]-x0;
        const double my_dy = xsp[isp][1]-y0;
        xsp[isp][0] = my_dx*cos(-theta0)-my_dy*sin(-theta0)+x0;
        xsp[isp][1] = my_dx*sin(-theta0)+my_dy*cos(-theta0)+y0;
      }
    }
  }

  // now update top and bottom
  for (int toz = stosz_i[top_id]; toz != stosz_i[top_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == top_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == 1) {
        sp_flag[isp] = 0; // operate once on point
        // corresponding cyl_id point
        const double cyl_theta = atan2(xsp[isp][1]-y0,xsp[isp][0]-x0);
        const double cyl_x = x0 + pow(fabs(cos(cyl_theta)),2.0/double(n))*a*SGN(cos(cyl_theta));
        const double cyl_y = y0 + pow(fabs(sin(cyl_theta)),2.0/double(n))*b*SGN(sin(cyl_theta));
        // const double cyl_r = sqrt( pow(cyl_x-x0,2) + pow(cyl_y-y0,2) );
        const double dr = sqrt( pow(xsp[isp][0]-x0,2) + pow(xsp[isp][1]-y0,2) ) - r0;
        // bump-based displacement psi(a) = exp(-1/(1-a^2)) for |a| < 1
        if (dr < L && dr > 0.0) {
          xsp[isp][0] += (cyl_x-x0)*exp(1.0-1.0/(1.0-dr*dr/L/L));
          xsp[isp][1] += (cyl_y-y0)*exp(1.0-1.0/(1.0-dr*dr/L/L));
          const double my_dx = xsp[isp][0]-x0;
          const double my_dy = xsp[isp][1]-y0;
          const double my_theta = -theta0*exp(1.0-1.0/(1.0-dr*dr/L/L));
          xsp[isp][0] = my_dx*cos(my_theta)-my_dy*sin(my_theta)+x0;
          xsp[isp][1] = my_dx*sin(my_theta)+my_dy*cos(my_theta)+y0;
          //cout << "theta " << cyl_theta << " dr/L " << dr/L << endl;
          //cout << " > cyl_x " << cyl_x << " cyl_y " << cyl_y << " cyl_r " << cyl_r << endl;
          //cout << " > my_x " << xsp[isp][0] << " my_y " << xsp[isp][1] << endl;
        }
      }
    }
  }
  for (int toz = stosz_i[bot_id]; toz != stosz_i[bot_id+1]; ++toz) {
    const int ist = stosz_v[toz];
    assert(szost[ist] == bot_id);
    FOR_I3 {
      const int isp = spost[ist][i];
      if (sp_flag[isp] == 1) {
        sp_flag[isp] = 0; // operate once on point
        // corresponding cyl_id point
        const double cyl_theta = atan2(xsp[isp][1]-y0,xsp[isp][0]-x0)-theta0;
        const double cyl_x = x0 + pow(fabs(cos(cyl_theta)),2/n)*a*SGN(cos(cyl_theta));
        const double cyl_y = y0 + pow(fabs(sin(cyl_theta)),2/n)*b*SGN(sin(cyl_theta));
        // const double cyl_r = sqrt( pow(cyl_x-x0,2) + pow(cyl_y-y0,2) );
        const double dr = sqrt( pow(xsp[isp][0]-x0,2) + pow(xsp[isp][1]-y0,2) ) - r0;
        // bump-based displacement psi(a) = exp(-1/(1-a^2)) for |a| < 1
        if (dr < L && dr > 0.0) {

          xsp[isp][0] += (cyl_x-x0)*exp(1.0-1.0/(1.0-dr*dr/L/L));
          xsp[isp][1] += (cyl_y-y0)*exp(1.0-1.0/(1.0-dr*dr/L/L));
          const double my_dx = xsp[isp][0]-x0;
          const double my_dy = xsp[isp][1]-y0;
          const double my_theta = -theta0*exp(1.0-1.0/(1.0-dr*dr/L/L));
          xsp[isp][0] = my_dx*cos(my_theta)-my_dy*sin(my_theta)+x0;
          xsp[isp][1] = my_dx*sin(my_theta)+my_dy*cos(my_theta)+y0;
          //cout << "theta " << cyl_theta << " dr/L " << dr/L << endl;
          //cout << " > cyl_x " << cyl_x << " cyl_y " << cyl_y << " cyl_r " << cyl_r << endl;
          //cout << " > my_x " << xsp[isp][0] << " my_y " << xsp[isp][1] << endl;
        }
      }
    }
  }

  // check the morph...
  bool bad_top = checkMorphedZone(top_id);
  bool bad_bot = checkMorphedZone(bot_id);
  bool bad_cyl = checkMorphedZone(cyl_id);

  if (bad_top || bad_bot || bad_cyl) {
    // revert back so you can redo the morph
    for (int isp = 0; isp < nsp; ++isp) FOR_I3 xsp[isp][i] = xsp0[isp][i];
  }
  else {
    // surface will change, recompute geoms
    clearSubzoneData();
    clearZoneData();
    b_centroid = false;
    b_rmax = false;
    b_boundingBox = false;
    clearNonManifoldData();

    // if we had periodic data, we need to reset it here
    if (pbi) {
      clearPeriodicity();
      WUI(WARN,"Bad news, friend. You need to reset periodicity.");
    }
  }
  delete[] xsp0;

}

void SimpleSurface::morphCylinderR(const int cyl_id,const double factor) {

  cout << "morphCylinderR: cyl_id: " << cyl_id << " factor: " << factor << endl;

  assert((cyl_id >= 0)&&(cyl_id < nsz));

  // we need teost to do edge walking...

  ensureTeost();

  st_flag.setLength(nst);
  st_flag.setAll(0);

  sp_flag.setLength(nsp);
  sp_flag.setAll(-1);

  double (*dx)[3] = new double[nsp][3];
  for (int isp = 0; isp < nsp; ++isp)
    FOR_I3 dx[isp][i] = 0.0;

  vector<int> spVec;
  stack<int> stStack;
  for (int ist_seed = 0; ist_seed < nst; ++ist_seed) {
    if ((szost[ist_seed] == cyl_id)&&(st_flag[ist_seed] == 0)) {
      // this is a seed -- march out to get the cylinder...
      double xc0[3] = { 0.0, 0.0, 0.0 };
      double wgt0 = 0.0;
      int n0 = 0;
      int isz0 = -1;
      double xc1[3] = { 0.0, 0.0, 0.0 };
      double wgt1 = 0.0;
      int n1 = 0;
      int isz1 = -1;
      st_flag[ist_seed] = 1;
      stStack.push(ist_seed);
      spVec.clear();
      while (!stStack.empty()) {
	const int ist = stStack.top(); stStack.pop();
	FOR_I3 {
	  const int isp = spost[ist][i];
	  if (sp_flag[isp] == -1) {
	    spVec.push_back(isp);
	    sp_flag[isp] = ist_seed;
	  }
	  else {
	    assert(sp_flag[isp] == ist_seed);
	  }
	}
	FOR_I3 {
	  // use teost to check edge nbr...
	  int ist_nbr,i_nbr,orient_nbr;
	  if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
	    if (szost[ist_nbr] != cyl_id) {
	      // this is an edge along another subzone...
	      const int isp0 = spost[ist][i];
	      const int isp1 = spost[ist][(i+1)%3];
	      const double wgt = DIST(xsp[isp0],xsp[isp1]);
	      if ((isz0 == -1)||(isz0 == szost[ist_nbr])) {
		isz0 = szost[ist_nbr];
		FOR_I3 xc0[i] += wgt*(xsp[isp0][i]+xsp[isp1][i]);
		wgt0 += wgt;
		++n0; // number of segments...
	      }
	      else if ((isz1 == -1)||(isz1 == szost[ist_nbr])) {
		isz1 = szost[ist_nbr];
		FOR_I3 xc1[i] += wgt*(xsp[isp0][i]+xsp[isp1][i]);
		wgt1 +=wgt;
		++n1;
	      }
	      else {
		cout << "Error: cylinder can only touch two other subzones" << endl;
		assert(0);
	      }
	    }
	    else if (st_flag[ist_nbr] == 0) {
	      st_flag[ist_nbr] = 1;
	      stStack.push(ist_nbr);
	    }
	  }
	  else {
	    cout << "Error: open edges not supported in cylinder selection" << endl;
	    assert(0);
	  }
	}
      }
      // make sure we got the two defining edges...
      assert(isz0 != -1);
      assert(isz1 != -1);
      // note that wgt0, wgt1 contains the circumference...
      cout << "isz0: " << isz0 << " d0: " << wgt0/M_PI << " n0: " << n0 << " isz1: " << isz1 << " d1: " << wgt1/M_PI << " n1: " << n1 << endl;
      FOR_I3 xc0[i] /= 2.0*wgt0;
      FOR_I3 xc1[i] /= 2.0*wgt1;
      cout << "xc0: " << COUT_VEC(xc0) << " xc1: " << COUT_VEC(xc1) << endl;
      // now go through and set dx...
      // start by turning xc1 into the cylinder axis unit normal...
      FOR_I3 xc1[i] -= xc0[i];
      const double mag = MAG(xc1);
      assert(mag > 0.0);
      FOR_I3 xc1[i] /= mag;

      // now loop on nodes and perturb...
      for (int ii = 0, limit = spVec.size(); ii < limit; ++ii) {
        const int isp = spVec[ii];
        assert(sp_flag[isp] == ist_seed);
        const double dxsp[3] = DIFF(xsp[isp],xc0);
        const double dp = DOT_PRODUCT(dxsp,xc1);
        double dr[3]; FOR_I3 dr[i] = dxsp[i] - dp*xc1[i];
        // const double r = MAG(dr);
        FOR_I3 dx[isp][i] = (factor-1.0)*dr[i];
      }

    }
  }

  // at this point, the morph is in dx, and nodes with dx have both non-zero dx and non-negative-one sp_flag...

  //smoothMorph(dx);

  MiscUtils::dumpRange(dx,nsp,"morph dx");
  for (int isp = 0; isp < nsp; ++isp)
    FOR_I3 xsp[isp][i] += dx[isp][i];

  delete[] dx;

}

void SimpleSurface::morphGaussianBump() {

  sp_flag.resize(nsp);
  sp_flag.setAll(0);
  for (int ist = 0; ist < nst; ++ist) {
    if (zoneVec[znost[ist]].flag) {
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

  const double L = 1.0;
  const double h = 0.085*L;
  const double x0 = 0.195*L;
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] >= 2) {
      assert(xsp[isp][1] == 0.0);
      xsp[isp][1] = h*exp(-xsp[isp][0]*xsp[isp][0]/(x0*x0));
    }
  }

}

void SimpleSurface::smoothMorph(double (*dx)[3]) {
  (void)(dx);
  // take some fraction of the motion...

  double (*dx0)[3] = new double[nsp][3];
  double *wgt0 = new double[nsp];

  int iter = 0;
  while (1) {

    ++iter;
    cout << "smoothing iter: " << iter << endl;

    // zero the correction...

    for (int isp = 0; isp < nsp; ++isp) {
      FOR_I3 dx0[isp][i] = 0.0;
      wgt0[isp] = 0.0;
    }

    // loop on tris to apply correction...

    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
        // consider this edge...
        const int isp0 = spost[ist][i];
        const int isp1 = spost[ist][(i+1)%3];
        if ((sp_flag[isp0] == -1)||(sp_flag[isp1] == -1)) {
          // const double delta = DIST(xsp[isp0],xsp[isp1]);
          if (sp_flag[isp0] == -1) {
            // const double dx_mag = MAG(dx[isp1]);
            // weight the nbr dx's with inverse distance, and
            // also only consider a fraction of these values...


          }
        }
      }
    }



  }

  assert(0);

  delete[] dx0;
  delete[] wgt0;

}

void SimpleSurface::nsmoothFlaggedSubzones(const int niter,const double relax,const bool b_poisson) {
  
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

  int nsp_smooth = 0;
  for (int isp = 0; isp < nsp; ++isp) {
    if (sp_flag[isp] == 2) {
      sp_flag[isp] = nsp_smooth++;
    }
    else {
      sp_flag[isp] = -1;
    }
  }
    
  if (b_poisson) {

    ensureTeost();
    st_flag.resize(nst);
    int nst_smooth = 0;
    for (int ist = 0; ist < nst; ++ist) {
      if (sz_flag[szost[ist]]) 
        st_flag[ist] = nst_smooth++;
      else
        st_flag[ist] = -1;
    }

    // build target normals for smooth tris via simple smoothing...

    double (*n_st_smooth)[3] = new double[nst_smooth][3]; 
    for (int ist = 0; ist < nst; ++ist) {
      if (st_flag[ist] >= 0) {
        const int ist_smooth = st_flag[ist];
        const double * const x0 = xsp[spost[ist][0]];
        const double * const x1 = xsp[spost[ist][1]];
        const double * const x2 = xsp[spost[ist][2]];
        const double n[3] = TRI_NORMAL_2(x0,x1,x2);
        FOR_I3 n_st_smooth[ist_smooth][i] = n[i];
      }
    }
    for (int iter = 0; iter < 1000; ++iter) {
      for (int ist = 0; ist < nst; ++ist) {
        if (st_flag[ist] >= 0) {
          const int ist_smooth = st_flag[ist];
          const double area_2 = MAG(n_st_smooth[ist_smooth]);
          FOR_I3 n_st_smooth[ist_smooth][i] = 0.0;
          FOR_I3 {
            int ist_nbr, i_nbr, orient_nbr;
            if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
              if (st_flag[ist_nbr] >= 0) {
                const int ist_nbr_smooth = st_flag[ist_nbr];
                FOR_J3 n_st_smooth[ist_smooth][j] += n_st_smooth[ist_nbr_smooth][j];
              }
              else {
                const double * const x0 = xsp[spost[ist_nbr][0]];
                const double * const x1 = xsp[spost[ist_nbr][1]];
                const double * const x2 = xsp[spost[ist_nbr][2]];
                const double n[3] = TRI_NORMAL_2(x0,x1,x2);
                FOR_J3 n_st_smooth[ist_smooth][j] += n[j];
              }
            }
          }
          const double n_mag = MAG(n_st_smooth[ist_smooth]);
          if (n_mag > 0.0) 
            FOR_I3 n_st_smooth[ist_smooth][i] *= area_2/n_mag;
        }
      }
    }

    // build surface-tri-of-surface-point to allow use to go from isp_smooth to ist

    int * stosp_i = new int[nsp_smooth+1]; 
    for (int isp = 0; isp < nsp_smooth; ++isp) stosp_i[isp+1] = 0;

    // count number of tris for each node
    for (int ist = 0; ist < nst; ++ist) {
      if (st_flag[ist] >= 0) {
        FOR_I3 {
          if (sp_flag[spost[ist][i]] >= 0) 
            stosp_i[sp_flag[spost[ist][i]]+1]++;
        }
      }
    }

    // scan to get disp from counts...
    stosp_i[0] = 0;
    for (int isp = 0; isp < nsp_smooth; ++isp) stosp_i[isp+1] += stosp_i[isp];

    // populate stosp_v....
    int *stosp_v = new int[stosp_i[nsp_smooth]];
    for (int ist = 0; ist < nst; ++ist) {
      if (st_flag[ist] >= 0) {
        FOR_I3 {
          if (sp_flag[spost[ist][i]] >= 0) 
            stosp_v[stosp_i[sp_flag[spost[ist][i]]]++] = ist;
        }
      }
    }

    // rewind...
    for (int isp = nsp_smooth-1; isp > 0; --isp) stosp_i[isp] = stosp_i[isp-1];
    stosp_i[0] = 0;

    double (*kappa_n_sp_smooth)[3] = new double[nsp_smooth][3];
    for (int iter = 0; iter < niter; ++iter) {

      // build target nodal normal from target normals...
     
      for (int isp = 0; isp < nsp; ++isp) {
        if (sp_flag[isp] >= 0) {
          const int isp_smooth = sp_flag[isp];
          FOR_I3 kappa_n_sp_smooth[isp_smooth][i] = 0.0;
          for (int top = stosp_i[isp_smooth]; top != stosp_i[isp_smooth+1]; ++top) {
            const int ist = stosp_v[top];
            assert(st_flag[ist] >= 0); 
            const int ist_smooth = st_flag[ist];
            FOR_J3 kappa_n_sp_smooth[isp_smooth][j] += n_st_smooth[ist_smooth][j];
          }
          const double n_mag = MAG(kappa_n_sp_smooth[isp_smooth]);
          if (n_mag > 0.0)
            FOR_I3 kappa_n_sp_smooth[isp_smooth][i] /= n_mag;
        }
      }

      // build target nodal curvature from div of target normals...

      for (int isp = 0; isp < nsp; ++isp) {
        if (sp_flag[isp] >= 0) {
          const int isp_smooth = sp_flag[isp];
          double kappa = 0.0;
          for (int top = stosp_i[isp_smooth]; top != stosp_i[isp_smooth+1]; ++top) {
            const int ist = stosp_v[top];
            assert(st_flag[ist] >= 0); 
            const int ist_smooth = st_flag[ist];
            int i = -1;
            for (i = 0; i < 3; ++i) {
              if (spost[ist][i] == isp)
                break;
            }
            const double * x0 = xsp[spost[ist][i]];
            const double * x1 = xsp[spost[ist][(i+1)%3]];
            const double * x2 = xsp[spost[ist][(i+2)%3]];
            const double dx[3] = DIFF(x2,x1);
            const double dx_mag = MAG(dx);
            const double n[3] = TRI_NORMAL_2(x0,x1,x2);
            const double n_mag = MAG(n);
            const double n_st_mag = MAG(n_st_smooth[ist_smooth]);
            if ((dx_mag > 0.0)&&(n_mag > 0.0)&&(n_st_mag > 0.0)) {
              const double wgt[3] = CROSS_PRODUCT(n,dx);
              kappa -= DOT_PRODUCT(n_st_smooth[ist_smooth],wgt)/(dx_mag*n_mag*n_st_mag);
            }
          }
          FOR_I3 kappa_n_sp_smooth[isp_smooth][i] *= kappa;
        }
      }

      // estimate eigenvalue

      double d2_min = HUGE_VAL;
      for (int ist = 0; ist < nst; ++ist) {
        if (st_flag[ist] >= 0) {
          FOR_I3 {
            const double * x0 = xsp[spost[ist][i]];
            const double * x1 = xsp[spost[ist][(i+1)%3]];
            d2_min = min(d2_min,DIST2(x1,x0));
          }
        }
      }
      //cout << " > iter, d2_min: " << iter << " " << d2_min << endl; 

      // update nodal locations by solving poisson equation with target curvature on rhs...

      for (int isp = 0; isp < nsp; ++isp) {
        if (sp_flag[isp] >= 0) {
          const int isp_smooth = sp_flag[isp];
          double kappa_n[3] = {0.0,0.0,0.0};
          double n_mag = 0.0;
          for (int top = stosp_i[isp_smooth]; top != stosp_i[isp_smooth+1]; ++top) {
            const int ist = stosp_v[top];
            assert(st_flag[ist] >= 0); 
            int i = -1;
            for (i = 0; i < 3; ++i) {
              if (spost[ist][i] == isp)
                break;
            }
            const double * x0 = xsp[spost[ist][i]];
            const double * x1 = xsp[spost[ist][(i+1)%3]];
            const double * x2 = xsp[spost[ist][(i+2)%3]];
            const double dx12[3] = DIFF(x2,x1);
            const double dx10[3] = DIFF(x0,x1);
            const double dx20[3] = DIFF(x0,x2);
            const double num1 = DOT_PRODUCT(dx12,dx10);
            const double den1 = CROSS_PRODUCT_MAG(dx12,dx10);
            const double num2 = -DOT_PRODUCT(dx12,dx20); // dx21 = -dx12
            const double den2 = CROSS_PRODUCT_MAG(dx12,dx20);
            if ((den1 > 0.0)&&(den2 > 0.0)) {
              const double cot1 = num1/den1;
              const double cot2 = num2/den2;
              FOR_J3 kappa_n[j] -= 0.5*(cot1*dx20[j]+cot2*dx10[j]); // 2*area
              const double n[3] = TRI_NORMAL_2(x0,x1,x2);
              n_mag += MAG(n);
            }
          }
          if (n_mag > 0.0) {
            FOR_I3 kappa_n[i] /= n_mag;
            //cout << COUT_VEC(kappa_n) << " " << COUT_VEC(kappa_n_sp_smooth[isp_smooth]) << endl; 
            FOR_I3 xsp[isp][i] += relax*0.5*d2_min*(kappa_n_sp_smooth[isp_smooth][i] - kappa_n[i]); 
          }
        }
      }
    
    }

    delete[] n_st_smooth;
    delete[] kappa_n_sp_smooth;
    delete[] stosp_i;
    delete[] stosp_v;

  }
  else {

    double (*xsp_smooth)[3] = new double[nsp_smooth][3];
    double *count_smooth = new double[nsp_smooth];
    for (int iter = 0; iter < niter; ++iter) {
      // zero xsp_smooth...
      for (int isp_smooth = 0; isp_smooth < nsp_smooth; ++isp_smooth) {
        FOR_I3 xsp_smooth[isp_smooth][i] = 0.0;
        count_smooth[isp_smooth] = 0.0;
      }
      for (int ist = 0; ist < nst; ++ist) {
        if (sz_flag[szost[ist]]) {
          FOR_I3 {
            const int isp = spost[ist][i];
            if (sp_flag[isp] >= 0) {
              const int isp_smooth = sp_flag[isp];
              FOR_J3 xsp_smooth[isp_smooth][j] += xsp[spost[ist][(i+1)%3]][j];
              FOR_J3 xsp_smooth[isp_smooth][j] += xsp[spost[ist][(i+2)%3]][j];
              count_smooth[isp_smooth] += 2.0;
            }
          }
        }
      }
      // and set the xsp to the smoothed values...
      for (int isp = 0; isp < nsp; ++isp) {
        if (sp_flag[isp] >= 0) {
          const int isp_smooth = sp_flag[isp];
          FOR_I3 xsp[isp][i] = relax*(xsp_smooth[isp_smooth][i]/count_smooth[isp_smooth])+(1.0-relax)*xsp[isp][i];
        }
      }
    }
  }
  
}

void SimpleSurface::nsmoothNearPoint(const int niter,const double x[3],const double r) {

  sp_flag.resize(nsp);
  for (int iter = 0; iter < niter; ++iter) {
    sp_flag.setAll(-1);
    int nsp_smooth = 0;
    for (int isp = 0; isp < nsp; ++isp) {
      if (DIST2(xsp[isp],x) < r*r) {
	sp_flag[isp] = nsp_smooth++;
      }
    }
    cout << " > smoothing iter " << iter << " of " << niter << ": " << nsp_smooth << " surface points..." << endl;
    double (*xsp_smooth)[3] = new double[nsp_smooth][3];
    int *count_smooth = new int[nsp_smooth];
    for (int isp_smooth = 0; isp_smooth < nsp_smooth; ++isp_smooth) {
      FOR_I3 xsp_smooth[isp_smooth][i] = 0.0;
      count_smooth[isp_smooth] = 0;
    }
    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
	const int isp = spost[ist][i]; 
	if (sp_flag[isp] >= 0) {
	  const int isp_smooth = sp_flag[isp];
	  FOR_J3 xsp_smooth[isp_smooth][j] += xsp[spost[ist][(i+1)%3]][j];
	  FOR_J3 xsp_smooth[isp_smooth][j] += xsp[spost[ist][(i+2)%3]][j];
	  count_smooth[isp_smooth] += 2;
	}
      }
    }
    for (int isp = 0; isp < nsp; ++isp) {
      if (sp_flag[isp] >= 0) {
	const int isp_smooth = sp_flag[isp];
	assert(count_smooth[isp_smooth] > 0);
	FOR_I3 xsp[isp][i] = xsp_smooth[isp_smooth][i]/double(count_smooth[isp_smooth]);
      }
    }
    delete[] xsp_smooth;
    delete[] count_smooth;
  }
  
}
					  
    


