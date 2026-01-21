void CUT(const double n[3],const int ifa_value) {
  CUT(n,n,ifa_value);
}

void CUT(const double x[3],const double n[3],const int ifa_value) {

  //static int debug = 0;
  //debug++;

  // volume cutting and surface cutting are VERY similar, except surface cutting
  // has incomplete edges -- i.e. the faoed on one side is cut.

  // cut the current volume with normal n, and
  // every new edge introduced gets faoed[ied][1] set to ifa_value

  assert(ifa_value != -1);

  /*
  double phi[nno];
  bool phi_pos = false;
  FOR_INO {
    phi[ino] = (x_no[ino][0]-n[0])*n[0] + (x_no[ino][1]-n[1])*n[1] + (x_no[ino][2]-n[2])*n[2];
    phi_pos = phi_pos | (phi[ino] > 0.0);
  }
  */

  double phi[nno];
  bool phi_pos = false;
  
  // this one tries a tolerance...
  // I did some experiments with the tolerance. In theory it should be something like 1.0E-32
  // (for 16 digits of precision). 1.0E-30 seems to be good and has significantly fewer unmatched
  // faces early on, while still meeting strict tolerances. 1.0E-28 gave some (small) tolerance
  // checking errors.
  
  const double tol2 = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  FOR_INO {
    phi[ino] = (x_no[ino][0]-x[0])*n[0] + (x_no[ino][1]-x[1])*n[1] + (x_no[ino][2]-x[2])*n[2];
    if (phi[ino]*phi[ino] < 1.0E-30*tol2*tol2)
      phi[ino] = 0.0;
    else
      phi_pos = (phi_pos || (phi[ino] > 0.0));
  }

  // if no positive, there is nothing to cut away...
  // could check for negative too, and quickly cut away everything, but this will only
  // happen when dealing with multi-volumes...

  if (!phi_pos)
    return;

  vector< std::pair<int,int> > fano0Vec;
  vector< std::pair<int,int> > fano1Vec;

  int nno_orig = nno;
  int no_flag[nno];
  FOR_INO no_flag[ino] = -1;

  int ned_new = 0;
  FOR_IED {
    const int ino0 = nooed[ied][0];
    const int ino1 = nooed[ied][1];
    if (phi[ino0] < 0.0) {
      if (phi[ino1] > 0.0) {
	// connect to a new node...
	const int ino_new = new_node();
	FOR_I3 {
	  // for the new intersection point, we can simplify the math under certain
	  // circumstances...
	  if ((n[(i+1)%3]==0.0)&&(n[(i+2)%3]==0.0)) {
	    x_no[ino_new][i] = x[i];
	  }
	  else if (x_no[ino0][i] == x_no[ino1][i]) {
	    x_no[ino_new][i] = x_no[ino0][i];
	  }
	  else {
	    x_no[ino_new][i] = (phi[ino1]*x_no[ino0][i] - phi[ino0]*x_no[ino1][i])/(phi[ino1]-phi[ino0]);
	  }
	}
	// set new edge start and end points...
	IF_FAOED_0 fano0Vec.push_back( std::pair<int,int>(faoed[ied][0],ino_new) );
	IF_FAOED_1 fano1Vec.push_back( std::pair<int,int>(faoed[ied][1],ino_new) );
	// keep ino0...
	no_flag[ino0] = 1;
	nooed[ned_new][0] = ino0;
	nooed[ned_new][1] = ino_new;
	faoed[ned_new][0] = faoed[ied][0];
	faoed[ned_new][1] = faoed[ied][1];
	++ned_new;
      }
      else if (phi[ino1] == 0.0) {
	IF_FAOED_0 fano0Vec.push_back( std::pair<int,int>(faoed[ied][0],ino1) );
	IF_FAOED_1 fano1Vec.push_back( std::pair<int,int>(faoed[ied][1],ino1) );
	// keep ino0,ino1...
	no_flag[ino0] = 1;
	no_flag[ino1] = 1;
	nooed[ned_new][0] = ino0;
	nooed[ned_new][1] = ino1;
	faoed[ned_new][0] = faoed[ied][0];
	faoed[ned_new][1] = faoed[ied][1];
	++ned_new;
      }
      else {
	// keep ino0,ino1...
	no_flag[ino0] = 1;
	no_flag[ino1] = 1;
	nooed[ned_new][0] = ino0;
	nooed[ned_new][1] = ino1;
	faoed[ned_new][0] = faoed[ied][0];
	faoed[ned_new][1] = faoed[ied][1];
	++ned_new;
      }
    }
    else if (phi[ino1] < 0.0) {
      if (phi[ino0] > 0.0) {
	// connect to a new node...
	const int ino_new = new_node();
	FOR_I3 {
	  // for the new intersection point, we can simplify the math under certain
	  // circumstances...
	  if ((n[(i+1)%3]==0.0)&&(n[(i+2)%3]==0.0)) {
	    x_no[ino_new][i] = x[i];
	  }
	  else if (x_no[ino0][i] == x_no[ino1][i]) {
	    x_no[ino_new][i] = x_no[ino0][i];
	  }
	  else {
	    x_no[ino_new][i] = (phi[ino0]*x_no[ino1][i] - phi[ino1]*x_no[ino0][i])/(phi[ino0]-phi[ino1]);
	  }
	}
	// set new edge start and end points...
	IF_FAOED_1 fano0Vec.push_back( std::pair<int,int>(faoed[ied][1],ino_new) );
	IF_FAOED_0 fano1Vec.push_back( std::pair<int,int>(faoed[ied][0],ino_new) );
	// keep ino1...
	no_flag[ino1] = 1;
	nooed[ned_new][0] = ino_new;
	nooed[ned_new][1] = ino1;
	faoed[ned_new][0] = faoed[ied][0];
	faoed[ned_new][1] = faoed[ied][1];
	++ned_new;
      }
      else if (phi[ino0] == 0.0) {
	// set new edge start and end points...
	IF_FAOED_1 fano0Vec.push_back( std::pair<int,int>(faoed[ied][1],ino0) );
	IF_FAOED_0 fano1Vec.push_back( std::pair<int,int>(faoed[ied][0],ino0) );
	// keep ino0,ino1...
	no_flag[ino0] = 1;
	no_flag[ino1] = 1;
	nooed[ned_new][0] = ino0;
	nooed[ned_new][1] = ino1;
	faoed[ned_new][0] = faoed[ied][0];
	faoed[ned_new][1] = faoed[ied][1];
	++ned_new;
      }
      else {
	// keep ino0,ino1...
	no_flag[ino0] = 1;
	no_flag[ino1] = 1;
	nooed[ned_new][0] = ino0;
	nooed[ned_new][1] = ino1;
	faoed[ned_new][0] = faoed[ied][0];
	faoed[ned_new][1] = faoed[ied][1];
	++ned_new;
      }
    }
  }

  ned = ned_new;

  // now add new edges from the fano0 and fano1 vecs...
  // start by sorting by face...

  sort(fano0Vec.begin(),fano0Vec.end());
  sort(fano1Vec.begin(),fano1Vec.end());
  assert(fano0Vec.size() == fano1Vec.size());

  /*
    for (int ii = 0; ii < fano0Vec.size(); ++ii)
    cout << "fano0Vec: " << ii << " " << fano0Vec[ii].first << " " << fano0Vec[ii].second << endl;
    for (int ii = 0; ii < fano1Vec.size(); ++ii)
    cout << "fano1Vec: " << ii << " " << fano1Vec[ii].first << " " << fano1Vec[ii].second << endl;
  */

  int i1 = 0;
  const int fano0Vec_size = fano0Vec.size();
  while (i1 != fano0Vec_size) {
    const int i0 = i1;
    ++i1;
    while ((i1 < fano0Vec_size)&&(fano0Vec[i1].first == fano0Vec[i0].first))
      ++i1;
    // the i0 through, but not including i1 range contains start and finish nodes...
    // check that the face is the same in the sorted lists...
    for (int i = i0; i < i1; ++i) {
      if (!(fano0Vec[i].first == fano1Vec[i].first)) {

	cout << "mpi_rank: " << mpi_rank << endl;
	for (int ii = 0; ii < fano0Vec_size; ++ii)
	  cout << "fano0Vec: " << ii << " " << fano0Vec[ii].first << " " << fano0Vec[ii].second << endl;
	for (int ii = 0; ii < fano0Vec_size; ++ii)
	  cout << "fano1Vec: " << ii << " " << fano1Vec[ii].first << " " << fano1Vec[ii].second << endl;

      }
      assert(fano0Vec[i].first == fano1Vec[i].first);
    }
    if (i1 == i0+1) {
      // this is typical case -- one start and one finish per face...
      // it is possible that the nodes are the same, in which
      // case no edge is required...
      if (fano0Vec[i0].second != fano1Vec[i0].second) {
	// this one is simple. Just join...
	const int ied_new = new_edge();
	// add the edge in the reverse direction to correspond with the way it was done earlier --
	// i.e. the faoed[][0] is the one that get ifa_value
	// I guess this was because the [0] reference eventually gets the internal face, rather
	// than the boundary face...
	nooed[ied_new][1] = fano0Vec[i0].second;
	nooed[ied_new][0] = fano1Vec[i0].second;
	faoed[ied_new][1] = fano0Vec[i0].first;
	faoed[ied_new][0] = ifa_value;
      }
    }
    else {
      // find the vector between the furthest points...
      double dx[3] = { 0.0, 0.0, 0.0 };
      double d2_max = 0.0;
      for (int i = i0; i < i1; ++i) {
	for (int j = i0; j < i1; ++j) {
	  const int ino0 = fano0Vec[i].second;
	  const int ino1 = fano1Vec[j].second;
	  const double this_dx[3] = DIFF(x_no[ino1],x_no[ino0]);
	  const double this_d2 = DOT_PRODUCT(this_dx,this_dx);
	  if (this_d2 > d2_max) {
	    FOR_K3 dx[k] = this_dx[k];
	    d2_max = this_d2;
	  }
	}
      }
      vector< std::pair<int,double> > no0Vec;
      vector< std::pair<int,double> > no1Vec;
      for (int i = i0; i < i1; ++i) {
	const int ino0 = fano0Vec[i].second;
	//cout << "[" << mpi_rank << "] x_no0 : " << ino0 << " " << COUT_VEC(x_no[ino0]) << endl;
	no0Vec.push_back( std::pair<int,double>(ino0,DOT_PRODUCT(x_no[ino0],dx)) );
	const int ino1 = fano1Vec[i].second;
	//cout << "[" << mpi_rank << "] x_no1 : " << ino1 << " " << COUT_VEC(x_no[ino1]) << endl;
	no1Vec.push_back( std::pair<int,double>(ino1,DOT_PRODUCT(x_no[ino1],dx)) );
      }
      std::sort(no0Vec.begin(),no0Vec.end(),MiscUtils::compareSecondIntDoublePair);
      std::sort(no1Vec.begin(),no1Vec.end(),MiscUtils::compareSecondIntDoublePair);
      for (int i = i0; i < i1; ++i) {
	//cout << "[" << mpi_rank << "] after sort no0Vec[i-i0].first: " << no0Vec[i-i0].first << " no0Vec[i-i0].second: " << no0Vec[i-i0].second << endl;
	//cout << "[" << mpi_rank << "] after sort no1Vec[i-i0].first: " << no1Vec[i-i0].first << " no1Vec[i-i0].second: " << no1Vec[i-i0].second << endl;
	if (no0Vec[i-i0].first != no1Vec[i-i0].first) {
	  const int ied_new = new_edge();
	  nooed[ied_new][1] = no0Vec[i-i0].first;
	  nooed[ied_new][0] = no1Vec[i-i0].first;
	  faoed[ied_new][1] = fano0Vec[i0].first;
	  faoed[ied_new][0] = ifa_value;
	}
      }
    }
  }

  // now condense nodes, and reset d2_max...
  d2_max = 0.0;
  FOR_I3 Lmax[i] = 0.0;
  FOR_I3 Lmin[i] = 0.0;
  int nno_new = 0;
  for (int ino = 0; ino < nno_orig; ++ino) {
    if (no_flag[ino] == 1) {
      // keep...
      FOR_I3 x_no[nno_new][i] = x_no[ino][i];
      const double d2 = DOT_PRODUCT(x_no[nno_new],x_no[nno_new]);
      d2_max = max(d2_max,d2);
      const double r = sqrt(d2);
      FOR_I3 {
	Lmax[i] = max(Lmax[i],r+x_no[nno_new][i]);
	Lmin[i] = min(Lmin[i],-r+x_no[nno_new][i]);
      }
      no_flag[ino] = nno_new++;
    }
    else {
      assert(no_flag[ino] == -1);
    }
  }
  for (int ino = nno_orig; ino < nno; ++ino) {
    FOR_I3 x_no[nno_new][i] = x_no[ino][i];
    const double d2 = DOT_PRODUCT(x_no[nno_new],x_no[nno_new]);
    d2_max = max(d2_max,d2);
    const double r = sqrt(d2);
    FOR_I3 {
      Lmax[i] = max(Lmax[i],r+x_no[nno_new][i]);
      Lmin[i] = min(Lmin[i],-r+x_no[nno_new][i]);
    }
    ++nno_new;
  }

  // and blow through edges and change node numbering...
  FOR_IED {
    FOR_I2 {
      const int ino_orig = nooed[ied][i];
      if (ino_orig < nno_orig) {
	nooed[ied][i] = no_flag[ino_orig];
	//if (!(nooed[ied][i] >= 0))
	//  cout << "ERROR rank: " << mpi_rank << " debug: " << debug << endl;
	assert(nooed[ied][i] >= 0);
      }
      else {
	nooed[ied][i] -= (nno-nno_new);
      }
    }
  }

  nno = nno_new;

}
