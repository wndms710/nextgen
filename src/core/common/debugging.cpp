

// ==============================
// Cht Helper Stuff

  // Conjugate heat transfer...
  
  class ChtHelper {
  public:
    string name;
    int index;
    int (**noote_ptr)[4];
    int (**noost_ptr)[3];
    int **znost_ptr;
    int* nte_ptr;
    int* nst_ptr;
    int* nno_ptr;

    ChtHelper() {
      name = "";
      index = -1;
      noote_ptr = NULL;
      noost_ptr = NULL;
      znost_ptr = NULL;
      nte_ptr = NULL;
      nst_ptr = NULL;
      nno_ptr = NULL;
    }

  };
  vector<ChtHelper> chtHelperVec;
  map<const string,int> chtHelperNameMap;

  void registerCht(const string &name,int (*&noote)[4],int (*&noost)[3],int *&znost,int &nte,int &nst,int &nno) {

    size_t colon_pos = name.find(":");
    assert(colon_pos == string::npos); // should be name of particle class
    map<const string,int>::iterator iter = chtHelperNameMap.find(name);
    assert(iter == chtHelperNameMap.end()); // no name repeats

    // populate helper vec

    chtHelperVec.push_back(ChtHelper());
    const int cht_index = chtHelperVec.size()-1;
    chtHelperNameMap[name] = cht_index;
    chtHelperVec.back().name = name;
    chtHelperVec.back().index = cht_index;
    chtHelperVec.back().noote_ptr = &noote;
    chtHelperVec.back().noost_ptr = &noost;
    chtHelperVec.back().znost_ptr = &znost;
    chtHelperVec.back().nte_ptr = &nte;
    chtHelperVec.back().nst_ptr = &nst;
    chtHelperVec.back().nno_ptr = &nno;

  }

  void registerChtData(double *&data,const string& name) {

    size_t colon_pos = name.find(":");
    assert(colon_pos != string::npos); // names must be prepended with particle class name
    const string cht_name  = name.substr(0,colon_pos);

    // get index into helper vec

    map<const string,int>::iterator iter = chtHelperNameMap.find(cht_name);
    assert(iter != chtHelperNameMap.end());
    const int cht_index = iter->second;
    assert((cht_index >= 0)&&(cht_index < chtHelperVec.size()));
    assert(chtHelperVec[cht_index].name  == cht_name);
    assert(chtHelperVec[cht_index].index == cht_index);

    // now register the data

    const int topo = (((cht_index+1)<<INDEX_SHIFT)|NO_DATA);
    CtiRegister::_registerData(data,name,topo,NO_READWRITE_DATA,*chtHelperVec[cht_index].nno_ptr);

  }

  void registerChtData(double (*&data)[3],const string& name) {

    size_t colon_pos = name.find(":");
    assert(colon_pos != string::npos); // names must be prepended with particle class name
    const string cht_name  = name.substr(0,colon_pos);

    // get index into helper vec

    map<const string,int>::iterator iter = chtHelperNameMap.find(cht_name);
    assert(iter != chtHelperNameMap.end());
    const int cht_index = iter->second;
    assert((cht_index >= 0)&&(cht_index < chtHelperVec.size()));
    assert(chtHelperVec[cht_index].name  == cht_name);
    assert(chtHelperVec[cht_index].index == cht_index);

    // now register the data

    const int topo = (((cht_index+1)<<INDEX_SHIFT)|NO_DATA);
    CtiRegister::_registerData(data,name,topo,NO_READWRITE_DATA,*chtHelperVec[cht_index].nno_ptr);

  }

// end of Cht Helper Stuff
// ==============================


// used to be in StaticSolver bbox calc

    // add in cht nodes?...
    for (int ii = 0, lim = chtHelperVec.size(); ii < lim; ++ii) {
      CtiRegister::CtiData * data = CtiRegister::getRegisteredCtiData(chtHelperVec[ii].name+":x_no",false);
      assert((data)&&(data->getType() == DN3_DATA)&&(data->getUnindexedTopology() == NO_DATA));
      const int nno_= *chtHelperVec[ii].nno_ptr;
      double (*x_no_)[3] = data->getDN3ptr();
      for (int ino = 0; ino < nno_; ++ino) {
        buf[0] = min(buf[0],x_no_[ino][0]);
        buf[1] = min(buf[1],-x_no_[ino][0]);
        buf[2] = min(buf[2],x_no_[ino][1]);
        buf[3] = min(buf[3],-x_no_[ino][1]);
        buf[4] = min(buf[4],x_no_[ino][2]);
        buf[5] = min(buf[5],-x_no_[ino][2]);
      }
    }












    int my_nst_flagged = 0; 
    for (int ist = ist_begin; ist < ist_end; ++ist) {
      const int isz = znost[ist];
      assert((isz >= 0)&&(isz < nsz));
      if (sz_flag[isz])
	++my_nst_flagged;
    }

    int * nst_flagged = new int[mpi_size];
    MPI_Allgather(&my_nst_flagged,1,MPI_INT,nst_flagged,1,MPI_INT,mpi_comm);

    int nst_flagged_sum = 0;
    for (int rank = 0; rank < mpi_size; ++rank) {
      nst_flagged_sum += nst_flagged[rank];
    }

    const int nst_flagged_local_target = (nst_flagged_sum + mpi_size - 1)/mpi_size;
    const int ist_flagged_begin = min(nst_flagged_sum,mpi_rank*nst_flagged_local_target);
    const int ist_flagged_end   = min(nst_flagged_sum,(mpi_rank+1)*nst_flagged_local_target);
    

    my_nst_flagged = 0;
    ist_end = nst;
    for (int ist = 0; ist < nst; ++ist) {
      const int isz = znost[ist];
      assert((isz >= 0)&&(isz < nsz));
      if (sz_flag[isz]) {
	if (ist_flagged_begin == my_nst_flagged)
	  ist_begin = ist;
	else if (ist_flagged_end == my_nst_flagged)
	  ist_end = ist;
	my_nst_flagged++;
      }
    }
    assert(my_nst_flagged == nst_flagged_sum);

    if (mpi_rank == 0) cout << "nst_flagged_sum = " << nst_flagged_sum << endl;

    int check_nst_flagged = 0;
    for (int ist = ist_begin; ist < ist_end; ++ist) {
      const int isz = znost[ist];
      assert((isz >= 0)&&(isz < nsz));
      if (sz_flag[isz])
	++check_nst_flagged;
    }

    assert(check_nst_flagged == ist_flagged_end - ist_flagged_begin);
    /*
    FOR_RANK {
      if (rank == mpi_rank)
	cout << "check_nst_flagged = " << check_nst_flagged << ", ist_flagged_end-ist_flagged_begin = " << ist_flagged_end - ist_flagged_begin 
	     << ", ist_begin = " << ist_begin << ", ist_end = " << ist_end << endl;
      MPI_Barrier(mpi_comm);
    }
    MPI_Pause("looks good!");
    */
    
    nst_local = ist_end-ist_begin;
    int * st_flag = new int[nst_local];
    for (int i = 0; i < nst_local; ++i) st_flag[i] = -1;

    //cout << "serial implementation of grouping" << endl;


#ifdef GATHERV_VERSION

    MPI_Sync("about to gather and unpack");
    
    // reduce to rank0 and build a tri-nbr structure...
    // recall we are sending 4 ints per edge. For huge surfaces (nst ~ 300M), this 
    // will eventually fail, so check for it...
    
    int my_ned_times_4 = intVec.size();
    MPI_Gather(&my_ned_times_4,1,MPI_INT,send_count,1,MPI_INT,0,mpi_comm);
    
    int * int_buf = NULL;
    int ned_times_4;
    if (mpi_rank == 0) {
      // we may already have send_disp...
      if (send_disp == NULL) send_disp = new int[mpi_size];
      int8 ned_times_4_check = 0;
      FOR_RANK ned_times_4_check += send_count[rank];
      assert(ned_times_4_check < TWO_BILLION);
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
      ned_times_4 = send_disp[mpi_size-1] + send_count[mpi_size-1];
      assert(ned_times_4 == ned_times_4_check);
      int_buf = new int[ned_times_4];
    }

    int * my_int_buf = (intVec.empty() ? NULL : &intVec.front());
    MPI_Gatherv(my_int_buf,my_ned_times_4,MPI_INT,int_buf,send_count,send_disp,MPI_INT,0,mpi_comm);
    intVec.clear();
    delete[] send_count;
    if (send_disp) delete[] send_disp;
    
    assert(teost == NULL);
    CTI_Mmap(teost,nst);
    if (mpi_rank == 0) {

      // as an initial condition, just set ist. At the end of this routine,
      // everyone should be set EXCEPT teost in periodic tris. We should not 
      // ever need these...

      for (int ist = 0; ist < nst; ++ist) 
        FOR_I3 teost[ist][i] = ist;
      
      cout << "I have an int buf with size: " << ned_times_4 << endl;
      
      assert(ned_times_4%4 == 0);
      for (int ii = 0; ii < ned_times_4; ii += 4) {
        const int ist =       int_buf[ii]; assert((ist >= 0)&&(ist < nst));
        const int i_and_bit = int_buf[ii+1];
        const int i = (i_and_bit&3); assert((i >= 0)&&(i < 3));
        const int bit = (i_and_bit>>2); assert((bit >= 0)&&(bit <= 1)); // for now 0 or 1
        const int ist2 =      int_buf[ii+2]; assert((ist2 >= 0)&&(ist2 < nst));
        const int j =         int_buf[ii+3]; assert((j >= 0)&&(j < 3));
        assert(teost[ist][i] == ist);
        teost[ist][i] = ist2;
        assert(teost[ist2][j] == ist2);
        teost[ist2][j] = ist;
      }

      delete[] int_buf;

    }
   
#else





  void CubicSpline::init(const double (*xp)[3],const int np) {
    // if we have only two points, then the bessier 
    /*
      for (int i = 0; i < n; ++i)
      cout << "XXXXX " << COUT_VEC(x[i]) << endl;
      getchar();
    */
    
    /*
      if (n == 2) {
      }
      else {
      }
    */
    

    double (*rhs)[3] = new double[np][3];
    FOR_I3 rhs[0][i] = 3.0*(xp[1][i]-xp[0][i]);
    for (int ip = 1; ip < np-1; ++ip)
      FOR_I3 rhs[ip][i] = 3.0*(xp[ip+1][i]-xp[ip-1][i]);
    FOR_I3 rhs[np-1][i] = 3.0*(xp[np-1][i]-xp[np-2][i]);
    
    const double a = 1.0, b = 4.0, c = 1.0;   
    double *cp = new double[np-1];
    assert(dp == NULL); 
    dp = new double[np][3];
    
    // forward sweep to get cp and dp...
    cp[0] = 2.0*c/b;
    FOR_I3 dp[0][i] = 2.0*rhs[0][i]/b;
    for (int ip = 1; ip < np-1; ++ip) {
      cp[ip] = c/(b-a*cp[ip-1]);
      FOR_I3 dp[ip][i] = (rhs[ip][i]-a*dp[ip-1][i])/(b-a*cp[ip-1]);
    }
    FOR_I3 dp[np-1][i] = (rhs[np-1][i]-a*dp[np-2][i])/(0.5*b-a*cp[np-2]);






  /*
  bool checkCurrentParam() {

    using namespace ParamsData;

    // simply validates the currentParam...
    return(currentParam != paramList.end());
  }

  bool checkCurrentParam(const int i) {

    using namespace ParamsData;

    // simply validates that the currentParam has a token at position i...
    return( (currentParam != paramList.end())&&(int(currentParam->tokens.size()) > i) );
  }

  bool setCurrentParam(const string& name) {

    using namespace ParamsData;

    // use the currentParam to locate this param's name in the list...
    for (currentParam = paramList.begin(); currentParam != paramList.end(); ++currentParam) {
      if (currentParam->name == name) {
        currentParam->touch();
        return(true);
      }
    }

    // if it was not found, add it to the unknownParamList
    for (list<string>::iterator up = unknownParamList.begin(); up != unknownParamList.end(); ++up)
      if (name == *up)
        return(false);
    unknownParamList.push_back(name);
    return(false);
  }
  */


  void addParamsFromArgs(int argc,char * argv[],const bool verbose) {

    assert(0);

    //if (verbose)
    //cout << " > addParamsFromArgs" << endl;

    string line;
    assert( line.length() == 0 );

    int iargc = 1;
    while (iargc < argc) {

      // look for args that start with --, and skip anything before this...
      if ( (strlen(argv[iargc]) >= 3) && (argv[iargc][0] == '-') && (argv[iargc][1] == '-') ) {

	if (line.length() > 0)
	  addParamFromString(line);

	// start the next...
	line = &(argv[iargc][2]);

      }
      else {

	line.append(" ");
	line.append(argv[iargc]);

      }

      iargc++;

    }

    // last one could be out of the loop...
    if (line.length() > 0)
      addParamFromString(line);

  }

  /*
  void clearCurrentParamFlag() {

    using namespace ParamsData;
    using namespace MpiStuff;

    if (!checkCurrentParam()) {
      if (mpi_rank == 0)
        cerr << "\n\n\n********************************************************************\n" <<
          "Error: currentParam is not valid." <<
          "\n********************************************************************\n\n\n" << endl;
      throw(0);
    }
    currentParam->clearFlag();

  }

  void setCurrentParamFlag() {

    using namespace ParamsData;
    using namespace MpiStuff;

    if (!checkCurrentParam()) {
      if (mpi_rank == 0)
        cerr << "\n\n\n********************************************************************\n" <<
          "Error: currentParam is not valid." <<
          "\n********************************************************************\n\n\n" << endl;
      throw(0);
    }
    currentParam->setFlag();

  }
  */



  /*
  void dumpCurrentParam(ostream& ofile) {

    using namespace ParamsData;
    using namespace MpiStuff;

    if (!checkCurrentParam()) {
      if (mpi_rank == 0)
        cerr << "\n\n\n********************************************************************\n" <<
          "Error: currentParam is not valid." <<
          "\n********************************************************************\n\n\n" << endl;
      throw(0);
    }
    currentParam->dump(ofile);

  }

  string getCurrentParamName() {

    using namespace ParamsData;
    using namespace MpiStuff;

    if (!checkCurrentParam()) {
      if (mpi_rank == 0)
        cerr << "\n\n\n********************************************************************\n" <<
          "Error: currentParam is not valid." <<
          "\n********************************************************************\n\n\n" << endl;
      throw(0);
    }
    return(currentParam->name);
  }

  int getCurrentParamSize() {

    using namespace ParamsData;
    using namespace MpiStuff;

    if (!checkCurrentParam()) {
      if (mpi_rank == 0)
        cerr << "\n\n\n********************************************************************\n" <<
          "Error: currentParam is not valid." <<
          "\n********************************************************************\n\n\n" << endl;
      throw(0);
    }
    return(currentParam->tokens.size());

  }

  string getCurrentParamString(const int i) {

    using namespace ParamsData;
    using namespace MpiStuff;

    if (!checkCurrentParam()) {
      if (mpi_rank == 0)
        cerr << "\n\n\n********************************************************************\n" <<
          "Error: currentParam is not valid." <<
          "\n********************************************************************\n\n\n" << endl;
      throw(0);
    }
    return(currentParam->getString(i));
  }

  float getCurrentParamFloat(const int i) {

    using namespace ParamsData;
    using namespace MpiStuff;

    if (!checkCurrentParam()) {
      if (mpi_rank == 0)
        cerr << "\n\n\n********************************************************************\n" <<
          "Error: currentParam is not valid." <<
          "\n********************************************************************\n\n\n" << endl;
      throw(0);
    }
    return(currentParam->getFloat(i));
  }

  double getCurrentParamDouble(const int i) {

    using namespace ParamsData;
    using namespace MpiStuff;

    if (!checkCurrentParam()) {
      if (mpi_rank == 0)
        cerr << "\n\n\n********************************************************************\n" <<
          "Error: currentParam is not valid." <<
          "\n********************************************************************\n\n\n" << endl;
      throw(0);
    }
    return(currentParam->getDouble(i));
  }

  int getCurrentParamInt(const int i) {

    using namespace ParamsData;
    using namespace MpiStuff;

    if (!checkCurrentParam()) {
      if (mpi_rank == 0)
        cerr << "\n\n\n********************************************************************\n" <<
          "Error: currentParam is not valid." <<
          "\n********************************************************************\n\n\n" << endl;
      throw(0);
    }
    return(currentParam->getInt(i));

  }

  int8 getCurrentParamInt8(const int i) {

    using namespace ParamsData;
    using namespace MpiStuff;

    if (!checkCurrentParam()) {
      if (mpi_rank == 0)
        cerr << "\n\n\n********************************************************************\n" <<
          "Error: currentParam is not valid." <<
          "\n********************************************************************\n\n\n" << endl;
      throw(0);
    }
    return(currentParam->getInt8(i));

  }

  bool getCurrentParamBool(const int i) {

    using namespace ParamsData;
    using namespace MpiStuff;

    if (!checkCurrentParam()) {
      if (mpi_rank == 0)
        cerr << "\n\n\n********************************************************************\n" <<
          "Error: currentParam is not valid." <<
          "\n********************************************************************\n\n\n" << endl;
      throw(0);
    }
    return(currentParam->getBool(i));

  }
  */




  void expandLoops() {
    
    //vector<list<Param>::iterator> forVec;
    //vector<list<Param>::iterator> endforVec;
    
    vector<pair<list<Param>::const_iterator,int> > forVec;
    stack<list<Param>::const_iterator> forStack;
    
    list<Param>::const_iterator iter = ParamsData::paramList.begin();
    while (iter != ParamsData::paramList.end()) {
      cout << "got param: " << iter->str() << endl;
      const string token = iter->getName();
      if (token == "FOR") {
	forStack.push(iter);
	int ii;
	for (ii = 0; ii < forVec.size(); ++ii) {
	  if (forVec[ii].first == iter)
	    break;
	}
	if (ii == forVec.size()) {
	  forVec.push_back(pair<list<Param>::const_iterator,int>(iter,5));
	  cout << " > this is the first time at this for: 5 " << endl;
	}
	else {
	  cout << " > found in forVec: " << forVec[ii].second << endl;
	}
      }
      else if (token == "ENDFOR") {
	// the latest "for" should be on the stack...
	if (forStack.empty()) {
	  CERR("Params::expandLoops: FOR/ENDFOR pairs not matched. Check input file.");
	}
	list<Param>::const_iterator iter_for = forStack.top(); forStack.pop();
	// check the counter in the forMap...
	int ii;
	for (ii = 0; ii < forVec.size(); ++ii) {
	  if (forVec[ii].first == iter_for)
	    break;
	}
	assert(ii < forVec.size()); // must find it
	cout << " > this is end for with loop index: " << forVec[ii].second << endl;
	--forVec[ii].second;
	if (forVec[ii].second > 0) {
	  iter = forVec[ii].first;
	  continue;
	}
      }
      // just proceed to the next param...
      ++iter;
      getchar();
    }


    cout << "LOOKS GOOD" << endl;
    getchar();

  }







    else {

      // on the second and later calls, we need to check and potentially modify the bins

      int r_bins = 0;
      if (-data_minmax[1] > range[1]) {
	// decide how many bins we need to add on the plus side to capture the largest data...
	r_bins = (int)ceil((-data_minmax[1]-range[1])/(range[1]-range[0])*double(nbin));
      }

      int l_bins = 0;
      if (data_minmax[0] < range[0]) {
	// decide how many bins we need to add on the minus side to capture the smallest data...
	l_bins = (int)ceil((range[0]-data_minmax[0])/(range[1]-range[0])*double(nbin));
      }

      // consider whether we need to condense bins...

      //if (mpi_rank == rank0) cout << "l_bins: " << l_bins << " r_bins: " << r_bins << endl;
      
      if ( (l_bins+r_bins+nbin) > 2*nbin_target) {
	
	// we can combine bins by some factor...
	int nbin_collapse = int((l_bins+r_bins+nbin)/nbin_target);
	//if (mpi_rank == rank0) cout << "l_bins + r_bins + nbin: " << l_bins + r_bins + nbin << " nbin_collapse: " << nbin_collapse << endl;

	// builds new bins and collapse existing data at the same time...
	
	// first collapse the old bins
	// range stays the same
	const int nbin_old = nbin;
	int nbin_temp = nbin_old + l_bins + r_bins;
	int add_bins = 0;
	int extra_l_bins = 0;
	int extra_r_bins = 0;
	bool bypass = false;
	bool inner_first = false;	

	if ( (nbin_temp%nbin_collapse) == 0 ) {
	  //if (mpi_rank == rank0) cout << "Resizing all bins - l_, r_ and old bins ALL AT ONCE!" << "\n";
	  //if (mpi_rank == rank0) cout << "nbin_collapse stays at " << nbin_collapse << "\n";
	  bypass = true;
	}

	if (((nbin_old%nbin_collapse)==0) && bypass==false) {
	  //if (mpi_rank == rank0) cout << "Resizing old bins first, then adding l_ and r_ bins!" << "\n";
	  //if (mpi_rank == rank0) cout << "nbin_collapse stays at " << nbin_collapse << "\n";
	  bypass = true;
	  inner_first = true;	
	}

	if (((nbin_old%nbin_collapse)!=0) && ((nbin_temp%nbin_collapse)!=0)) {
	  //if (mpi_rank == rank0) cout << "nbin_collapse not acceptable, must add extra empty bins..." << "\n";
	  add_bins = (nbin_collapse) - (nbin_temp % nbin_collapse);
	  //if (mpi_rank == rank0) cout << "add bins is " << add_bins << "\n";
          for (int abin = 1; abin <= add_bins; ++abin) {
            if (abin % 2) {
              // Add to left side
              extra_l_bins += 1;
            }
            else {
              extra_r_bins += 1;
            }
          }
          inner_first = false;
          nbin_temp += (extra_l_bins + extra_r_bins);
          l_bins += extra_l_bins;
          r_bins += extra_r_bins;
          //if (mpi_rank == rank0) cout << "After adding " << extra_l_bins << " extra l_bins and " << extra_r_bins << " extra r_bins, nbin_temp is now " << nbin_temp << " which is divisible by " << nbin_collapse << "\n"; 
	}
	assert(((nbin_old % nbin_collapse)==0) || ((nbin_temp % nbin_collapse)==0));
	// End of analysis

	double * count_old = count;	
	count = new double[nbin_temp];
	// recalculate l_bins and r_bins	
        double *count_inner = new double[nbin_old];
        for (int ibin = 0; ibin < nbin_old; ++ibin) 
          count_inner[ibin] = 0.0;
	int nbin_new;
	int nbin_inner_new;
	double old_bin_width = ((range[1]-range[0])/double(nbin_old));
	double new_bin_width;

	for (int cbin = 0; cbin < nbin_temp; ++cbin) {
	  count[cbin] = 0.0;
	}

        // default: resize l_ r_ and inner bins all at once...
	if (!inner_first) {
          //if (mpi_rank == rank0) cout << "Starting resizing now for all bins!" << "\n";
          nbin_new = nbin_temp / nbin_collapse;
          //if (mpi_rank == rank0) cout << "Taking " << nbin_temp << " bins down to " << nbin_new << "\n\n";
  
          for (int ibin = 0; ibin < nbin_new; ++ibin) {
            for (int jbin = 0; jbin < nbin_collapse; ++jbin) {
              if ((((ibin*nbin_collapse)+jbin-l_bins) >= 0) && (((ibin*nbin_collapse)+jbin-l_bins) < (nbin_old))) {
                count[ibin] += count_old[(ibin*nbin_collapse)+jbin-l_bins];
                //if (mpi_rank == rank0) cout << "New bin count " << ibin << ": " << count[ibin] << "\n";
              }
            }	
          }

          // Account for l_ and r_bins
          new_bin_width = old_bin_width*double(nbin_collapse);
          range[0] -= (old_bin_width*l_bins);
          range[1] += (old_bin_width*r_bins);	
          //if (mpi_rank == rank0) cout << " > setting range to: " << range[0] << " " << range[1] << endl;
          //if (mpi_rank == rank0) cout << "Old bin width was " << old_bin_width << ", new bin width will be " << new_bin_width << "\n";
	}

        // default: resize inner (old) bins first
	if (inner_first) {
          nbin_inner_new = nbin_old / nbin_collapse;
          new_bin_width = old_bin_width*double(nbin_collapse);
          //if (mpi_rank == rank0) cout << "Old bin width was " << old_bin_width << ", new bin width will be " << new_bin_width << "\n";
          //if (mpi_rank == rank0) cout << "nbin_inner_new is " << nbin_inner_new << "\n";

          for (int ibin = 0; ibin < nbin_inner_new; ++ibin) {
            for (int jbin = 0; jbin < nbin_collapse; ++jbin) {
              count_inner[ibin] += count_old[(ibin*nbin_collapse)+jbin];
              //if (mpi_rank == rank0) cout << "Inner ind + " << count_old[(ibin*nbin_collapse)+jbin] << " from " << (ibin*nbin_collapse)+jbin << "\n";
            }
            //if (mpi_rank == rank0) cout << "Inner indice " << ibin << " has count " << count_inner[ibin] << "\n\n";
          }

          l_bins = (int)ceil(double(l_bins)/double(nbin_collapse));
          r_bins = (int)ceil(double(r_bins)/double(nbin_collapse));
          
          // now resize l_bins and r_bins...
          nbin_new = nbin_inner_new + l_bins + r_bins;
          //if (mpi_rank == rank0) cout << "New l_bins and r_bins are " << l_bins << " and " << r_bins << ", respectively." << "\n";
          //if (mpi_rank == rank0) cout << "Total new nbin is " << nbin_new << "\n";

          for (int kbin = 0; kbin < nbin_new-r_bins; ++kbin) {
            count[kbin+l_bins] = count_inner[kbin];
          }

          //if (mpi_rank == rank0) {
          //  for (int lbin = 0; lbin < nbin_new; ++lbin) {	
          //    cout << "New count " << lbin << " has value " << count[lbin] << "\n";
          //  }
          //}
          
          range[0] -= (new_bin_width*l_bins);
          range[1] += (new_bin_width*r_bins);
          //if (mpi_rank == rank0) cout << " > setting range to: " << range[0] << " " << range[1] << endl;
	}	
	nbin = nbin_new;
	
	delete[] count_old;
	delete[] count_inner;
      }
      else if ((l_bins+r_bins) > 0) {
	
	// just add the new bins...
	const int nbin_old = nbin;
	nbin = nbin_old + l_bins + r_bins;

	// adjust the range...
	const double bin_width = (range[1]-range[0])/double(nbin_old);
	//if (mpi_rank == rank0) cout << "Check: bin_width is " << bin_width << "\n";
	//if (mpi_rank == rank0) cout << "Since r1 and r0 are " << range[1] << " and " << range[0] << " and nbin_old is " << nbin_old << "\n";
	range[0] -= bin_width*double(l_bins);
	range[1] += bin_width*double(r_bins);
	//if (mpi_rank == rank0) cout << " > setting range to: " << range[0] << " " << range[1] << endl;
	//if (mpi_rank == rank0)cout << "Check: r_ and l_bins are : " << r_bins << " and " << l_bins << "\n";		
	// allocate new count, storing a ptr to the old count...
	double * count_old = count;
	count = new double[nbin];
	
	// zero the new count...
	for (int ibin = 0; ibin < nbin; ++ibin)
	  count[ibin] = 0.0;
	
	// copy old data into the new count in the correct bin...
	for (int ibin = 0; ibin < nbin_old; ++ibin)
	  count[ibin+l_bins] = count_old[ibin];
	delete[] count_old;

      }
      
    }






  void initProbes() {

    FOR_PARAM_MATCHING("PROBE") {
      addPointProbe(&(*param));
    }
    
    FOR_PARAM_MATCHING("MULTIFLUX_PROBE") {
      addMultiFluxProbe(&(*param));
    }
    
  }

void processMultiFluxProbe(Param * param) {
    
    MultiFluxProbe mfp;
    const int ierr = initMultiFluxProbe(mfp,param);
    if (ierr != 0) {
      CWARN("MultiFluxProbe could not be initialized");
    }
    else {
      mfp.doProbe();
      mfp.flush();
    }
    
    // mfp is destroyed on exit...
    
  }
  
  void addMultiFluxProbe(Param * param) {
    
    // put an empty MultiFluxProbe() into the mfpMap based on the full param as a string...
    
    pair<map<const string,MultiFluxProbe>::iterator,bool> return_pair = mfpMap.insert(pair<const string,MultiFluxProbe>(param->str(),MultiFluxProbe()));
    assert(return_pair.second); // should insert successfully

    // now initialize...

    const int ierr = initMultiFluxProbe(return_pair.first->second,&(*param));

    // if init failed, remove from the map...
    if (ierr != 0) {
      CWARN("MultiFluxProbe could not be added");
      mfpMap.erase(return_pair.first);
    }

  }

void addPointProbe(Param * param) {
    
    // put an empty PointProbe() into the ppMap based on the full param as a string...
    
    pair<map<const string,PointProbe>::iterator,bool> return_pair = ppMap.insert(pair<const string,PointProbe>(param->str(),PointProbe()));
    assert(return_pair.second); // should insert successfully

    // now initialize...

    const int ierr = initPointProbe(return_pair.first->second,&(*param));

    // if init failed, remove from the map...
    if (ierr != 0) {
      CWARN("PointProbe could not be added");
      ppMap.erase(return_pair.first);
    }

  }
  /*
  CtiDataError _sqrtCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {
  
    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);
    
    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();
    
    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == D_DATA) {
      v.new_d();
      if ( b_eval_func) v.d() = sqrt(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      if ( b_eval_func) FOR_I3 v.d3(i) = sqrt(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      if ( b_eval_func) { 
        for (int i = 0; i < v.size(); ++i) {
          v.dn(i) = sqrt(arg->dn(i));
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      if ( b_eval_func) { 
        for (int i = 0; i < v.size(); ++i) {
          FOR_J3 v.dn3(i,j) = sqrt(arg->dn3(i,j));
        }
      }
      return(CTI_DATA_OK);
    }
    
    return(CTI_DATA_NOT_VALID);
    
  }
  
  CtiDataError _expCtiData(CtiData& v,list<CtiData>& args, const bool b_eval_func) {
  
    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);
    
    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();
    
    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == D_DATA) {
      v.new_d();
      if ( b_eval_func) v.d() = exp(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      if ( b_eval_func) FOR_I3 v.d3(i) = exp(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      if ( b_eval_func) { 
        for (int i = 0; i < v.size(); ++i) {
          v.dn(i) = exp(arg->dn(i));
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      if ( b_eval_func) { 
        for (int i = 0; i < v.size(); ++i) {
          FOR_J3 v.dn3(i,j) = exp(arg->dn3(i,j));
        }
      }
      return(CTI_DATA_OK);
    }
    
    return(CTI_DATA_NOT_VALID);
    
  }

  CtiDataError _logCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {
  
    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);
    
    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();
    
    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == D_DATA) {
      v.new_d();
      if ( b_eval_func) v.d() = log(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      if ( b_eval_func) FOR_I3 v.d3(i) = log(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      if ( b_eval_func) { 
        for (int i = 0; i < v.size(); ++i) {
          v.dn(i) = log(arg->dn(i));
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      if ( b_eval_func) { 
        for (int i = 0; i < v.size(); ++i) {
          FOR_J3 v.dn3(i,j) = log(arg->dn3(i,j));
        }
      }
      return(CTI_DATA_OK);
    }
    
    return(CTI_DATA_NOT_VALID);
    
  }

  CtiDataError _log10CtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {
  
    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);
    
    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();
    
    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == D_DATA) {
      v.new_d();
      if ( b_eval_func) v.d() = log10(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      if ( b_eval_func) FOR_I3 v.d3(i) = log10(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      if ( b_eval_func) { 
        for (int i = 0; i < v.size(); ++i) {
          v.dn(i) = log10(arg->dn(i));
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      if ( b_eval_func) { 
        for (int i = 0; i < v.size(); ++i) {
          FOR_J3 v.dn3(i,j) = log10(arg->dn3(i,j));
        }
      }
      return(CTI_DATA_OK);
    }
    
    return(CTI_DATA_NOT_VALID);
    
  }

  CtiDataError _cosCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {
  
    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);
    
    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();
    
    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == D_DATA) {
      v.new_d();
      if ( b_eval_func) v.d() = cos(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      if ( b_eval_func) FOR_I3 v.d3(i) = cos(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      if ( b_eval_func) { 
        for (int i = 0; i < v.size(); ++i) {
          v.dn(i) = cos(arg->dn(i));
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      if ( b_eval_func) { 
        for (int i = 0; i < v.size(); ++i) {
          FOR_J3 v.dn3(i,j) = cos(arg->dn3(i,j));
        }
      }
      return(CTI_DATA_OK);
    }
    
    return(CTI_DATA_NOT_VALID);
    
  }

  CtiDataError _sinCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {
  
    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);
    
    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();
    
    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == D_DATA) {
      v.new_d();
      if ( b_eval_func) v.d() = sin(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      if ( b_eval_func) FOR_I3 v.d3(i) = sin(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      if ( b_eval_func) { 
        for (int i = 0; i < v.size(); ++i) {
          v.dn(i) = sin(arg->dn(i));
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      if ( b_eval_func) { 
        for (int i = 0; i < v.size(); ++i) {
          FOR_J3 v.dn3(i,j) = sin(arg->dn3(i,j));
        }
      }
      return(CTI_DATA_OK);
    }
    
    return(CTI_DATA_NOT_VALID);
    
  }

  CtiDataError _tanCtiData(CtiData& v,list<CtiData>& args, const bool b_eval_func) {
  
    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);
    
    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();
    
    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == D_DATA) {
      v.new_d();
      if ( b_eval_func) v.d() = tan(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      if ( b_eval_func) FOR_I3 v.d3(i) = tan(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      if ( b_eval_func) { 
        for (int i = 0; i < v.size(); ++i) {
          v.dn(i) = tan(arg->dn(i));
        }
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      if ( b_eval_func) { 
        for (int i = 0; i < v.size(); ++i) {
          FOR_J3 v.dn3(i,j) = tan(arg->dn3(i,j));
        }
      }
      return(CTI_DATA_OK);
    }
    
    return(CTI_DATA_NOT_VALID);
    
  }

  CtiDataError _coshCtiData(CtiData& v,list<CtiData>& args,const bool b_eval_func) {
  
    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);
    
    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();
    
    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == D_DATA) {
      v.new_d();
      if ( b_eval_func) v.d() = cosh(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      if ( b_eval_func) FOR_I3 v.d3(i) = cosh(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      for (int i = 0; i < v.size(); ++i) {
	v.dn(i) = cosh(arg->dn(i));
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      for (int i = 0; i < v.size(); ++i) {
	FOR_J3 v.dn3(i,j) = cosh(arg->dn3(i,j));
      }
      return(CTI_DATA_OK);
    }
    
    return(CTI_DATA_NOT_VALID);
    
  }

  CtiDataError _sinhCtiData(CtiData& v,list<CtiData>& args) {
  
    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);
    
    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();
    
    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == D_DATA) {
      v.new_d();
      v.d() = sinh(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      FOR_I3 v.d3(i) = sinh(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      for (int i = 0; i < v.size(); ++i) {
	v.dn(i) = sinh(arg->dn(i));
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      for (int i = 0; i < v.size(); ++i) {
	FOR_J3 v.dn3(i,j) = sinh(arg->dn3(i,j));
      }
      return(CTI_DATA_OK);
    }
    
    return(CTI_DATA_NOT_VALID);
    
  }

  CtiDataError _tanhCtiData(CtiData& v,list<CtiData>& args) {
  
    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);
    
    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();
    
    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == D_DATA) {
      v.new_d();
      v.d() = tanh(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      FOR_I3 v.d3(i) = tanh(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      for (int i = 0; i < v.size(); ++i) {
	v.dn(i) = tanh(arg->dn(i));
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      for (int i = 0; i < v.size(); ++i) {
	FOR_J3 v.dn3(i,j) = tanh(arg->dn3(i,j));
      }
      return(CTI_DATA_OK);
    }
    
    return(CTI_DATA_NOT_VALID);
    
  }

  CtiDataError _acosCtiData(CtiData& v,list<CtiData>& args) {
  
    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);
    
    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();
    
    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == D_DATA) {
      v.new_d();
      v.d() = acos(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      FOR_I3 v.d3(i) = acos(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      for (int i = 0; i < v.size(); ++i) {
	v.dn(i) = acos(arg->dn(i));
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      for (int i = 0; i < v.size(); ++i) {
	FOR_J3 v.dn3(i,j) = acos(arg->dn3(i,j));
      }
      return(CTI_DATA_OK);
    }
    
    return(CTI_DATA_NOT_VALID);
    
  }

  CtiDataError _asinCtiData(CtiData& v,list<CtiData>& args) {
  
    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);
    
    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();
    
    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == D_DATA) {
      v.new_d();
      v.d() = asin(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      FOR_I3 v.d3(i) = asin(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      for (int i = 0; i < v.size(); ++i) {
	v.dn(i) = asin(arg->dn(i));
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      for (int i = 0; i < v.size(); ++i) {
	FOR_J3 v.dn3(i,j) = asin(arg->dn3(i,j));
      }
      return(CTI_DATA_OK);
    }
    
    return(CTI_DATA_NOT_VALID);
    
  }

  CtiDataError _atanCtiData(CtiData& v,list<CtiData>& args) {
  
    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);
    
    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();
    
    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == D_DATA) {
      v.new_d();
      v.d() = atan(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      FOR_I3 v.d3(i) = atan(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      for (int i = 0; i < v.size(); ++i) {
	v.dn(i) = atan(arg->dn(i));
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      for (int i = 0; i < v.size(); ++i) {
	FOR_J3 v.dn3(i,j) = atan(arg->dn3(i,j));
      }
      return(CTI_DATA_OK);
    }
    
    return(CTI_DATA_NOT_VALID);
    
  }

  CtiDataError _ceilCtiData(CtiData& v,list<CtiData>& args) {
  
    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);
    
    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();
    
    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == D_DATA) {
      v.new_d();
      v.d() = ceil(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      FOR_I3 v.d3(i) = ceil(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      for (int i = 0; i < v.size(); ++i) {
	v.dn(i) = ceil(arg->dn(i));
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      for (int i = 0; i < v.size(); ++i) {
	FOR_J3 v.dn3(i,j) = ceil(arg->dn3(i,j));
      }
      return(CTI_DATA_OK);
    }
    
    return(CTI_DATA_NOT_VALID);
    
  }

  CtiDataError _floorCtiData(CtiData& v,list<CtiData>& args) {
  
    // 1. check arg count...
    if (args.size() != 1)
      return(CTI_DATA_ARG_COUNT);
    
    // 2. grab the arg...
    list<CtiData>::iterator arg = args.begin();
    
    // 3. and apply the function...
    const int datatype = arg->getType();
    if (datatype == D_DATA) {
      v.new_d();
      v.d() = floor(arg->d());
      return(CTI_DATA_OK);
    }
    else if (datatype == D3_DATA) {
      v.new_d3();
      FOR_I3 v.d3(i) = floor(arg->d3(i));
      return(CTI_DATA_OK);
    }
    else if (datatype == DN_DATA) {
      v.new_dn(*arg);
      for (int i = 0; i < v.size(); ++i) {
	v.dn(i) = floor(arg->dn(i));
      }
      return(CTI_DATA_OK);
    }
    else if (datatype == DN3_DATA) {
      v.new_dn3(*arg);
      for (int i = 0; i < v.size(); ++i) {
	FOR_J3 v.dn3(i,j) = floor(arg->dn3(i,j));
      }
      return(CTI_DATA_OK);
    }
    
    return(CTI_DATA_NOT_VALID);
    
  }
  */


void CtiCanvas::calcVolvisIntensityAndDepth(float &intensity,float &depth,vector<VolVisData>::iterator& iter_begin,vector<VolVisData>::iterator& iter_end) const {

  // recall flag can be 0,1,2 for internal face, or
  // flag = 3 +
  // if (show) flag |= 4;
  // if (dir) flag |= 8;

  const double intensity_eps = 1.0E-32;

  // integrate in double, then return float
  double intensity_dble = 0.0;
  double fabs_intensity_dble = 0.0;
  double depth_dble = 0.0;
  double intensity_current = 0.0;
  float depth_current = iter_begin->depth;
  float delta_sum = 0.f;
  bool b_v0 = false;
  bool b_v1 = false;
  float v0,v1;

  if (volvis_type == "X-RAY") {

    //if (debug) cout << "starting depth: " << depth_current << endl;
    bool b_closed = true;
    for (vector<VolVisData>::iterator iter = iter_begin; iter != iter_end; ++iter) {

      // surface
      if ((iter->flag & 3) == 3) {
        // the 8 bit determines the direction...
        if (iter->flag & 8) {
          // this is a start.
          if (!b_closed)
            delta_sum += depth_current-iter->depth; // opposite sign due to direction of depth
          depth_current = iter->depth;
          b_closed = false; // we use this just in case pixel round off switches order b/w open and close
        }
        else {
          // this is an end.
          delta_sum += depth_current-iter->depth; // opposite sign due to direction of depth
          intensity_dble += delta_sum;
          depth_dble += 0.5*(depth_current+iter->depth)*(intensity_dble-intensity_current);
          intensity_current = intensity_dble;
          depth_current = iter->depth;
          delta_sum = 0.f;
          b_closed = true;
        }
        // the 4 bit determines show or no-show...
        if (iter->flag & 4) {
          // if the surface is "show", then we clear out the integrated intensity and depth
          intensity_dble = 0.0;
          depth_dble = 0.0;
        }
        else {
          //if (debug) cout << "This is no-show" << endl;
        }
      }
      // internal/interproc face (assumed empty!)
      else {
        assert(0);
      }

    }

  }
  else if (volvis_type == "LINEAR") {

    //if (debug) cout << "starting depth: " << depth_current << endl;
    for (vector<VolVisData>::iterator iter = iter_begin; iter != iter_end; ++iter) {

    // surface
    if ((iter->flag & 3) == 3) {
    // the 8 bit determines the direction...
    if (iter->flag & 8) {
    // this is a start. There should be no delta in sum (should have already been delt with)...
    //assert(delta_sum == 0.f); // bet this fails some day
    depth_current = iter->depth;
    b_v0 = b_v1 = false;
    //if (debug) cout << "This is start" << endl;
    }
    else {
    // this is an end. process anything in delta_sum and zero it...
    delta_sum += depth_current - iter->depth;
    if (b_v0) {
    // TODO v0 should be replaced with some sort of kernel
    //intensity_dble += v0*delta_sum;
    intensity_dble += max(min(v0,volvis_aux_data[1]),volvis_aux_data[0])*delta_sum;
    depth_dble += 0.5*(depth_current+iter->depth)*(intensity_dble-intensity_current);
    intensity_current = intensity_dble;
    b_v0 = false;
    delta_sum = 0.f;
    }
    depth_current = iter->depth;
    //if (debug) cout << "This is end" << endl;
    }
    // the 4 bit determines show or no-show...
    if (iter->flag & 4) {
    //if (debug) cout << "This is show" << endl;
    // if the surface is "show", then we clear out the integrated intensity and depth
    intensity_dble = 0.0;
    depth_dble = 0.0;
    }
    else {
    //if (debug) cout << "This is no-show" << endl;
    }
    }
    // internal/interproc face
    else {
    if (iter->flag == 0) {
    //if (debug) cout << "0 set" << endl;
    b_v0 = true;
    v0 = iter->v0;
    }
    else if (iter->flag == 1) {
    //if (debug) cout << "1 set" << endl;
    b_v1 = true;
    v1 = iter->v1;
    }
    else if (iter->flag == 2) {
    //if (debug) cout << "both set" << endl;
    b_v0 = true;
    v0 = iter->v0;
    b_v1 = true;
    v1 = iter->v1;
    }
    else {
    assert(0);
    }
    // process any depth and reset...
    delta_sum += depth_current - iter->depth;
    if (b_v1) {
    //intensity_dble += v1*delta_sum;
    intensity_dble += max(min(v1,volvis_aux_data[1]),volvis_aux_data[0])*delta_sum;
    depth_dble += 0.5*(depth_current+iter->depth)*(intensity_dble-intensity_current);
    intensity_current = intensity_dble;
    b_v1 = false;
    delta_sum = 0.f;
    }
    depth_current = iter->depth;
    }

    }

  }
  else if (volvis_type == "GAUSSIAN") {

    //if (debug) cout << "starting depth: " << depth_current << endl;
    for (vector<VolVisData>::iterator iter = iter_begin; iter != iter_end; ++iter) {

      // surface
      if ((iter->flag & 3) == 3) {
        // the 8 bit determines the direction...
        if (iter->flag & 8) {
          // this is a start. There should be no delta in sum (should have already been delt with)...
          //assert(delta_sum == 0.f); // bet this fails some day
          depth_current = iter->depth;
          b_v0 = b_v1 = false;
          //if (debug) cout << "This is start" << endl;
        }
        else {
          // this is an end. process anything in delta_sum and zero it...
          delta_sum += depth_current - iter->depth;

          if (b_v0) {
            // TODO v0 should be replaced with some sort of kernel
            //intensity_dble += v0*delta_sum;
            intensity_dble += exp(-0.5*(v0-volvis_aux_data[0])*(v0-volvis_aux_data[0])/
                                  (volvis_aux_data[1]*volvis_aux_data[1]))/(volvis_aux_data[1]*sqrt(2.0*M_PI))*delta_sum;
            depth_dble += 0.5*(depth_current+iter->depth)*(intensity_dble-intensity_current);
            intensity_current = intensity_dble;
            b_v0 = false;
            delta_sum = 0.f;
          }
          depth_current = iter->depth;
          //if (debug) cout << "This is end" << endl;
        }
        // the 4 bit determines show or no-show...
        if (iter->flag & 4) {
          //if (debug) cout << "This is show" << endl;
          // if the surface is "show", then we clear out the integrated intensity and depth
          intensity_dble = 0.0;
          depth_dble = 0.0;
        }
        else {
          //if (debug) cout << "This is no-show" << endl;
        }
      }
      // internal/interproc face
      else {
        if (iter->flag == 0) {
          //if (debug) cout << "0 set" << endl;
          b_v0 = true;
          v0 = iter->v0;
        }
        else if (iter->flag == 1) {
          //if (debug) cout << "1 set" << endl;
          b_v1 = true;
          v1 = iter->v1;
        }
        else if (iter->flag == 2) {
          //if (debug) cout << "both set" << endl;
          b_v0 = true;
          v0 = iter->v0;
          b_v1 = true;
          v1 = iter->v1;
        }
        else {
          assert(0);
        }
        // process any depth and reset...
        delta_sum += depth_current - iter->depth;

        if (b_v1) {
          //intensity_dble += v1*delta_sum;
          intensity_dble += exp(-0.5*(v1-volvis_aux_data[0])*(v1-volvis_aux_data[0])/
                                (volvis_aux_data[1]*volvis_aux_data[1]))/(volvis_aux_data[1]*sqrt(2.0*M_PI))*delta_sum;
          depth_dble += 0.5*(depth_current+iter->depth)*(intensity_dble-intensity_current);
          intensity_current = intensity_dble;
          b_v1 = false;
          delta_sum = 0.f;
        }
        depth_current = iter->depth;
      }

    }
  
  }
  else {
    assert(volvis_type == "SURFACE");

    //if (debug) cout << "starting depth: " << depth_current << endl;

    for (vector<VolVisData>::iterator iter = iter_begin; iter != iter_end; ++iter) {

      // surface
      if ((iter->flag & 3) == 3) {
        // the 8 bit determines the direction...
        if ((depth_current-iter->depth > 1.0E-6)||iter == iter_begin) {
          intensity_dble = intensity_dble*volvis_aux_data[0]+(1.0-volvis_aux_data[0])*iter->v0; // fresnel effect
          depth_dble += 0.5*(depth_current+iter->depth)*(intensity_dble-intensity_current);
          intensity_current = intensity_dble;
          depth_current = iter->depth;
          //depth_dble = depth_current;
        }
        if (iter->flag & 8) {
          // this is a start.
        }
        else {
          // this is an end.
        }
        // the 4 bit determines show or no-show...
        if (iter->flag & 4) {
          // if the surface is "show", then we clear out the integrated intensity and depth
          intensity_dble = 0.0;
          depth_dble = 0.0;
        }
        else {
          //if (debug) cout << "This is no-show" << endl;
        }
      }
      // internal/interproc face (assumed empty!)
      else {
        assert(0);
      }

    }
  }

  if (intensity_dble != 0.0) {
    intensity = float(intensity_dble);
    // TODO: hack: if fabs_intensity_dble not worked through above, set it to intensity_dble here...
    // at present, it is only worked through LINEAR...
    if (fabs_intensity_dble == 0.0)
      fabs_intensity_dble = intensity_dble;
  }
  else {
    intensity = 0.f;
  }

  if (fabs_intensity_dble != 0.0) {
    depth = float(depth_dble/fabs_intensity_dble); // center of intensity based depth
  }
  else {
    depth = 0.f;
  }
  
  return;
}

void CtiCanvas::addEdge(const double xe0[3],const double xe1[3],const uint2 zone) {

  const float d0 = CALC_PIXEL_DEPTH(xe0);
  const float d1 = CALC_PIXEL_DEPTH(xe1);

  if (max(d0,d1) <= 0.f) return;

  const float i0 = CALC_PIXEL_I(xe0);
  const float i1 = CALC_PIXEL_I(xe1);

  if ( (max(i0,i1) < 0.f) || (min(i0,i1) > float(ni-1)) ) return;

  const float j0 = CALC_PIXEL_J(xe0);
  const float j1 = CALC_PIXEL_J(xe1);

  if ( (max(j0,j1) < 0.f) || (min(j0,j1) > float(nj-1)) ) return;

  // render...

  float normal[3];
  float r = 3.f;

  const float denom = (i1-i0)*(i1-i0) + (j1-j0)*(j1-j0);
  if (denom == 0.f) {

    // render point at i0,j0...
    const int i_f = max(0,(int)ceil(i0-r));
    const int i_l = min(ni-1,(int)floor(i0+r));
    const int j_f = max(0,(int)ceil(j0-r));
    const int j_l = min(nj-1,(int)floor(j0+r));
    for (int i = i_f; i <= i_l; ++i) {
      for (int j = j_f; j <= j_l; ++j) {
        normal[0] = float(i)-i0;
        normal[1] = float(j)-j0;
        normal[2] = r*r - normal[0]*normal[0] - normal[1]*normal[1];
        if (normal[2] > 0.f) {
          normal[2] = sqrt(normal[2])/r;
          const float d = d0-normal[2]*r;
          if ((d > 0.f)&&(!isBlanked(float(i),float(j),d))) {
            // finish normalization...
            normal[0] /= r;
            normal[1] /= r;
            setPixelData(i,j,normal,d,0.f,zone); // i,j,normal,depth,data,zone
          }
        }
      }
    }

  }
  else {

    // render edge...
    // NOT quite right yet...
    const int i_f = max(0,(int)ceil(min(i0,i1)-r));
    const int i_l = min(ni-1,(int)floor(max(i0,i1)+r));
    const int j_f = max(0,(int)ceil(min(j0,j1)-r));
    const int j_l = min(nj-1,(int)floor(max(j0,j1)+r));
    for (int i = i_f; i <= i_l; ++i) {
      for (int j = j_f; j <= j_l; ++j) {
        const float t = ((i1-i0)*(float(i)-i0)+(j1-j0)*(float(j)-j0))/denom;
        if (t <= 0.f) {
          // sphere at i0,j0...
          normal[0] = float(i)-i0;
          normal[1] = float(j)-j0;
          normal[2] = r*r - normal[0]*normal[0] - normal[1]*normal[1];
          if (normal[2] > 0.f) {
            normal[2] = sqrt(normal[2])/r;
            const float d = d0-normal[2]*r;
            if ((d > 0.f)&&(!isBlanked(float(i),float(j),d))) {
              normal[0] /= r;
              normal[1] /= r;
              setPixelData(i,j,normal,d,0.f,zone); // i,j,normal,depth,data,zone
            }
          }
        }
        else if (t >= 1.f) {
          // sphere at i1,j1...
          normal[0] = float(i)-i1;
          normal[1] = float(j)-j1;
          normal[2] = r*r - normal[0]*normal[0] - normal[1]*normal[1];
          if (normal[2] > 0.f) {
            normal[2] = sqrt(normal[2])/r;
            const float d = d1-normal[2]*r;
            if ((d > 0.f)&&(!isBlanked(float(i),float(j),d))) {
              normal[0] /= r;
              normal[1] /= r;
              setPixelData(i,j,normal,d,0.f,zone); // i,j,normal,depth,data,zone
            }
          }
        }
        else {
          // along cylinder...
          normal[0] = float(i)-((1.f-t)*i0 + t*i1);
          normal[1] = float(j)-((1.f-t)*j0 + t*j1);
          normal[2] = r*r - normal[0]*normal[0] - normal[1]*normal[1];
          if (normal[2] > 0.f) {
            normal[2] = sqrt(normal[2])/r;
            const float d = (1.f-t)*d0+t*d1-normal[2]*r;
            if ((d > 0.f)&&(!isBlanked(float(i),float(j),d))) {
              normal[0] /= r;
              normal[1] /= r;
              setPixelData(i,j,normal,d,0.f,zone); // i,j,normal,depth,data,zone
            }
          }
        }
      }
    }
  }
}

void CtiCanvas::writeImageVolVis(const string& filename) {

  unsigned char (*rgb)[3] = new unsigned char[ni*nj][3];

  for (int ij = 0; ij < ni*nj; ++ij) {
    rgb[ij][0] = 73;
    rgb[ij][1] = 175;
    rgb[ij][2] = 205;
  }

  float * intensity = new float[ni*nj];
  for (int ij = 0; ij < ni*nj; ++ij) intensity[ij] = 0.f;

  /*
    class VolVisData {
    public:
    int i,j;
    float depth;
    bool in;
    VolVisData(const int i,const int j,const float depth,const bool in) {
    this->i = i;
    this->j = j;
    this->sdepth = depth;
    this->in = in;
    }
    };
  */

  sort(volVisDataVec.begin(),volVisDataVec.end());

  int ij_end = 0;
  while (ij_end < int(volVisDataVec.size())) {

    const int ij_begin = ij_end;
    const int i = volVisDataVec[ij_begin].i;
    const int j = volVisDataVec[ij_begin].j;

    ij_end = ij_begin+1;
    while ((ij_end < int(volVisDataVec.size()))&&(volVisDataVec[ij_end].i == i)&&(volVisDataVec[ij_end].j == j))
      ++ij_end;

    // take a look...
    int state[2] = { 0, 0 };
    float depth_prev[2];

    int image_ij = j*ni+i;
    rgb[image_ij][0] = 0;

    for (int ij = ij_begin; ij < ij_end; ++ij) {
      //cout << "ij: " << ij << " i: " << i << " j: " << j << " depth: " << volVisDataVec[ij].depth << " flag: " << volVisDataVec[ij].flag << endl;
      //getchar();
      switch(volVisDataVec[ij].flag) {
      case 2:
	// into zone 1...
	state[0] = 1;
	depth_prev[0] = volVisDataVec[ij].depth;
	break;
      case 1:
	if (state[0] == 1) {
	  state[0] = 0;
	  intensity[image_ij] += volVisDataVec[ij].depth-depth_prev[0];
	}
	depth_prev[0] = volVisDataVec[ij].depth;
	break;
      case 4:
	// into zone 2...
	state[1] = 1;
	depth_prev[1] = volVisDataVec[ij].depth;
	break;
      case 3:
	if (state[1] == 1) {
	  state[1] = 0;
	  intensity[image_ij] -= volVisDataVec[ij].depth-depth_prev[1];
	}
	depth_prev[1] = volVisDataVec[ij].depth;
	break;
      default:
	assert(0);
      }
    }

  }

  // now use the intensity to set rgb...

  float abs_intensity_max = -1.0E+20;
  for (int ij = 0; ij < ni*nj; ++ij) {
    if (abs(intensity[ij]) > abs_intensity_max) abs_intensity_max = abs(intensity[ij]);
  }

  cout << "abs_intensity_max: " << abs_intensity_max << endl;

  abs_intensity_max = 10.f;

  for (int ij = 0; ij < ni*nj; ++ij) {
    if (rgb[ij][0] == 0) {
      rgb[ij][0] = max(0,min(255,(int)(128.f + 128.f*intensity[ij]/abs_intensity_max)));
      rgb[ij][1] = rgb[ij][0];
      rgb[ij][2] = rgb[ij][0];
    }
  }

  PngImage png(ni,nj,rgb);

  string tmp_filename = MiscUtils::makeTmpPrefix(filename);
  png.write(tmp_filename.c_str());
  rename(tmp_filename.c_str(),filename.c_str());

  delete[] rgb; rgb = NULL;

  delete[] intensity;

}

void CtiCanvas::addInternalTriVolVis(const double xt0[3],const double xt1[3],const double xt2[3],const float v0,const float v1,const int flag) {

  // this routine adds VolVisData to the volVisDataVec for the internal face intersections that occur
  // along each pixel line-of-sight. 1 or 2 data values are passed into v0 and v1, depending on flag...

  // v0 should be the closest to the camera (i.e. smallest depth), and v1 further (i.e. larger depth)...

  // flag should be one of 3 types:
  // 0: only v0 is valid
  // 1: only v1 is valid
  // 2: v0 and v1 both valid

  assert((flag >= 0)&&(flag <= 2));

  const float d0 = CALC_PIXEL_DEPTH(xt0);
  const float d1 = CALC_PIXEL_DEPTH(xt1);
  const float d2 = CALC_PIXEL_DEPTH(xt2);

  //if (max(d0,max(d1,d2)) <= 0.f) return; // allow negative depth here for now

  const float i0 = CALC_PIXEL_I(xt0);
  const float i1 = CALC_PIXEL_I(xt1);
  const float i2 = CALC_PIXEL_I(xt2);

  if ( (max(i0,max(i1,i2)) < 0.f) || (min(i0,min(i1,i2)) > float(ni-1)) ) return;

  const float j0 = CALC_PIXEL_J(xt0);
  const float j1 = CALC_PIXEL_J(xt1);
  const float j2 = CALC_PIXEL_J(xt2);

  if ( (max(j0,max(j1,j2)) < 0.f) || (min(j0,min(j1,j2)) > float(nj-1)) ) return;

  // blanked applied later...

  // this routine doesn't use the normal, so I think we can skip this, unless it contributes
  // to robustness in some way...
  //  const float normal[3] = { (j1-j0)*(d2-d0)-(d1-d0)*(j2-j0), (d1-d0)*(i2-i0)-(i1-i0)*(d2-d0), (i1-i0)*(j2-j0)-(j1-j0)*(i2-i0) };
  //  const double mag = sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );
  //  if (mag == 0.f)
  //  return;
  // FOR_I3 normal[i] /= mag;

  // sort the data by y...
  float x0_s[3],x1_s[3],x2_s[3];
  if (j0 <= j1) {
    if (j1 <= j2) {
      x0_s[0] = i0; x0_s[1] = j0; x0_s[2] = d0;
      x1_s[0] = i1; x1_s[1] = j1; x1_s[2] = d1;
      x2_s[0] = i2; x2_s[1] = j2; x2_s[2] = d2;
    }
    else if (j0 <= j2) {
      x0_s[0] = i0; x0_s[1] = j0; x0_s[2] = d0;
      x1_s[0] = i2; x1_s[1] = j2; x1_s[2] = d2;
      x2_s[0] = i1; x2_s[1] = j1; x2_s[2] = d1;
    }
    else {
      x0_s[0] = i2; x0_s[1] = j2; x0_s[2] = d2;
      x1_s[0] = i0; x1_s[1] = j0; x1_s[2] = d0;
      x2_s[0] = i1; x2_s[1] = j1; x2_s[2] = d1;
    }
  }
  else if (j0 <= j2) {
    x0_s[0] = i1; x0_s[1] = j1; x0_s[2] = d1;
    x1_s[0] = i0; x1_s[1] = j0; x1_s[2] = d0;
    x2_s[0] = i2; x2_s[1] = j2; x2_s[2] = d2;
  }
  else if (j1 <= j2) {
    x0_s[0] = i1; x0_s[1] = j1; x0_s[2] = d1;
    x1_s[0] = i2; x1_s[1] = j2; x1_s[2] = d2;
    x2_s[0] = i0; x2_s[1] = j0; x2_s[2] = d0;
  }
  else {
    x0_s[0] = i2; x0_s[1] = j2; x0_s[2] = d2;
    x1_s[0] = i1; x1_s[1] = j1; x1_s[2] = d1;
    x2_s[0] = i0; x2_s[1] = j0; x2_s[2] = d0;
  }

  // check y-sort...

  assert(x0_s[1] <= x1_s[1]);
  assert(x1_s[1] <= x2_s[1]);

  // the sorting has divided the tri into 2 vertical pieces...

  const float pixel_eps = 1.0E-6f;

  // these can be equal, but we skip this case...

  if (x0_s[1] < x1_s[1]) {
    const int j0_plus = max(0,(int)ceil(x0_s[1]-pixel_eps));
    const int j1_minus = min(nj-1,(int)floor(x1_s[1]+pixel_eps));
    for (int j = j0_plus; j <= j1_minus; ++j) {
      // get the x01 as a linear weight on the x0_s->x1_s line...
      const float w01 = max(0.f,min(1.f,(x1_s[1]-float(j))/(x1_s[1] - x0_s[1])));
      const float x01 = w01*x0_s[0] + (1.0f-w01)*x1_s[0];
      const float w02 = max(0.f,min(1.f,(x2_s[1]-float(j))/(x2_s[1] - x0_s[1])));
      const float x02 = w02*x0_s[0] + (1.0f-w02)*x2_s[0];
      // the x values can be in either order...
      if (x01 < x02) {
        const int i_f = max(0,(int)ceil(x01-pixel_eps));
        const int i_l = min(ni-1,(int)floor(x02+pixel_eps));
        if (i_f <= i_l) {
          const float d01 = w01*x0_s[2] + (1.0f-w01)*x1_s[2];
          const float d02 = w02*x0_s[2] + (1.0f-w02)*x2_s[2];
          for (int i = i_f; i <= i_l; ++i) {
            const float w = max(0.f,min(1.f,(x02 - float(i))/(x02 - x01)));
            const float d = w*d01 + (1.0f-w)*d02;
	    volVisDataVec.push_back(VolVisData(j*ni+i,d,v0,v1,flag));
          }
        }
      }
      else if (x02 < x01) {
        const int i_f = max(0,(int)ceil(x02-pixel_eps));
        const int i_l = min(ni-1,(int)floor(x01+pixel_eps));
        if (i_f <= i_l) {
          const float d01 = w01*x0_s[2] + (1.0f-w01)*x1_s[2];
          const float d02 = w02*x0_s[2] + (1.0f-w02)*x2_s[2];
          for (int i = i_f; i <= i_l; ++i) {
            const float w = max(0.f,min(1.f,(x01 - float(i))/(x01 - x02)));
            const float d = w*d02 + (1.0f-w)*d01;
	    volVisDataVec.push_back(VolVisData(j*ni+i,d,v0,v1,flag));
          }
        }
      }
    }
  }

  if (x1_s[1] < x2_s[1]) {
    const int j1_plus = max(0,(int)ceil(x1_s[1]-pixel_eps));
    const int j2_minus = min(nj-1,(int)floor(x2_s[1]+pixel_eps));
    for (int j = j1_plus; j <= j2_minus; ++j) {
      // get the x12 as a linear weight on the x1_s->x2_s line...
      const float w12 = max(0.f,min(1.f,(x2_s[1]-float(j))/(x2_s[1] - x1_s[1])));
      const float x12 = w12*x1_s[0] + (1.0f-w12)*x2_s[0];
      const float w02 = max(0.f,min(1.f,(x2_s[1]-float(j))/(x2_s[1] - x0_s[1])));
      const float x02 = w02*x0_s[0] + (1.0f-w02)*x2_s[0];
      // the x values can be in either order...
      if (x12 < x02) {
        const int i_f = max(0,(int)ceil(x12-pixel_eps));
        const int i_l = min(ni-1,(int)floor(x02+pixel_eps));
        if (i_f <= i_l) {
          const float d12 = w12*x1_s[2] + (1.0f-w12)*x2_s[2];
          const float d02 = w02*x0_s[2] + (1.0f-w02)*x2_s[2];
          for (int i = i_f; i <= i_l; ++i) {
            const float w = max(0.f,min(1.f,(x02 - float(i))/(x02 - x12)));
            const float d = w*d12 + (1.0f-w)*d02;
	    volVisDataVec.push_back(VolVisData(j*ni+i,d,v0,v1,flag));
          }
        }
      }
      else if (x02 < x12) {
        const int i_f = max(0,(int)ceil(x02-pixel_eps));
        const int i_l = min(ni-1,(int)floor(x12+pixel_eps));
        if (i_f <= i_l) {
          const float d12 = w12*x1_s[2] + (1.0f-w12)*x2_s[2];
          const float d02 = w02*x0_s[2] + (1.0f-w02)*x2_s[2];
          for (int i = i_f; i <= i_l; ++i) {
            const float w = max(0.f,min(1.f,(x12 - float(i))/(x12 - x02)));
            const float d = w*d02 + (1.0f-w)*d12;
	    volVisDataVec.push_back(VolVisData(j*ni+i,d,v0,v1,flag));
	  }
        }
      }
    }
  }

}

void CtiCanvas::addPlaneMesh(const double (*x_vv)[3],const double * r_vv,const int *i_vv,const int n,const double factor) {

  double ijk[3];
  double nijk[3];
  FOR_I3 ijk[i] = dataPlaneData->center[i]; 
  FOR_I3 nijk[i] = dataPlaneData->normal[i];
  if (nijk[2] == 0.f) return;

  const double nijk_mag = MAG(nijk); assert(nijk_mag > 0.0);
  const double unit_nijk[3] = {nijk[0]/nijk_mag,nijk[1]/nijk_mag,nijk[2]/nijk_mag};
  const float normal[3] = {(float)unit_nijk[0],(float)unit_nijk[1],(float)unit_nijk[2]};

  for (int ivv = 0; ivv < n; ++ivv) {

    // compute formulas below with i,j,k in delta mode relative to the plane ijk.
    // this must be corrected when
    double i = CALC_PIXEL_I(x_vv[ivv]);
    double j = CALC_PIXEL_J(x_vv[ivv]);
    double k = CALC_PIXEL_DEPTH(x_vv[ivv]);

    // compute distance from the plane...

    const double dp = (i-ijk[0])*unit_nijk[0] + (j-ijk[1])*unit_nijk[1] + (k-ijk[2])*unit_nijk[2];
    // the radius of the point on the plane...
    const double r2 = r_vv[ivv]*r_vv[ivv]*factor*factor*double(ni*ni)/(width*width) - dp*dp;
    if (r2 > 0.0) {
      // move the i,j,k point onto the plane...
      i -= dp*unit_nijk[0];
      j -= dp*unit_nijk[1];
      k -= dp*unit_nijk[2];
      // now i,j,k is in the imaging data plane. Figure out the limits of the projected sphere r2 in ip,jp...
      // for these details see maple worksheet voronoi plane.mw...
      const double tmp1 = (nijk[0]*nijk[0]+nijk[2]*nijk[2]); assert(tmp1 > 0.0);
      const double djp = sqrt(r2*tmp1)/nijk_mag;
      const int jp_min = max((int)ceil(j-djp),0);
      const int jp_max = min((int)floor(j+djp),nj-1);
      const double d2_min = dp*dp;
      const double d2_max = r2+d2_min;
      for (int jp = jp_min; jp <= jp_max; ++jp) {
	const double tmp2 = sqrt(r2*tmp1 - (double(jp)-j)*(double(jp)-j)*nijk_mag*nijk_mag)*nijk[2]*(nijk[0] >= 0.0 ? 1.0 : -1.0);
	const double ip0 = -(tmp2 + (double(jp)-j)*nijk[0]*nijk[1])/tmp1;
	const double ip1 = (tmp2 - (double(jp)-j)*nijk[0]*nijk[1])/tmp1;
	const int ip_min = max((int)ceil(i+min(ip0,ip1)),0);
	const int ip_max = min((int)floor(i+max(ip0,ip1)),ni-1);
	for (int ip = ip_min; ip <= ip_max; ++ip) {
          const float depth = dataPlaneData->center[2] - (dataPlaneData->normal[0]*(ip-dataPlaneData->center[0]) +
                                                          dataPlaneData->normal[1]*(jp-dataPlaneData->center[1]))/dataPlaneData->normal[2]; 
          if (!isBlankedSkipDataPlane(ip,jp,depth,true)) {
            int IJ,ij; 
            getIJij(IJ,ij,ip,jp);
            float paux; 
            uint2 bits;
            if (getPixelPauxBits(paux,bits,IJ,ij)) {
              if (bits&INTERNAL_DATA_BIT) {
                assert(bits&MESH_BIT);
                assert(bits&MASK_BIT);
                assert(bits&AUX_BIT); // should not have other volume data concurrently
                const double dkp = ((double(ip)-i)*nijk[0]+(double(jp)-j)*nijk[1])/nijk[2];
                const double my_d2 = dp*dp + (double(ip)-i)*(double(ip)-i) + (double(jp)-j)*(double(jp)-j) + dkp*dkp;
                if (my_d2 < paux) 
                  setPixelPnormalPdepthPdataPauxPzone(IJ,ij,normal,depth,1.0-(my_d2-d2_min)/(d2_max-d2_min),my_d2,i_vv[ivv]); // leave bits unchanged
                  //setPixelPnormalPdepthPdataPauxPzone(IJ,ij,normal,depth,sqrt(my_d2),my_d2,i_vv[ivv]); // leave bits unchanged
              }
              else {
                // use the paux to store the SQUARE of the Euclidean distance of the ORIGINAL point to the rendered pixel...
                const double dkp = ((double(ip)-i)*nijk[0]+(double(jp)-j)*nijk[1])/nijk[2];
                paux = dp*dp + (double(ip)-i)*(double(ip)-i) + (double(jp)-j)*(double(jp)-j) + dkp*dkp;
                // TODO throw unique cell qualifier into data and do reduction on flush
                setPixelPnormalPdepthPdataPauxPzoneBits(IJ,ij,normal,depth,1.0-(paux-d2_min)/(d2_max-d2_min),paux,i_vv[ivv],bits|INTERNAL_DATA_BIT|MESH_BIT|AUX_BIT|MASK_BIT);
                //setPixelPnormalPdepthPdataPauxPzoneBits(IJ,ij,normal,depth,sqrt(paux),paux,i_vv[ivv],bits|INTERNAL_DATA_BIT|MESH_BIT|AUX_BIT|MASK_BIT);
              }
            }
            else {
              pbd[IJ] = new PixelBlockData();
              bits = 0;
              // use the paux to store the SQUARE of the Euclidean distance of the ORIGINAL point to the rendered pixel...
              const double dkp = ((double(ip)-i)*nijk[0]+(double(jp)-j)*nijk[1])/nijk[2];
              paux = dp*dp + (double(ip)-i)*(double(ip)-i) + (double(jp)-j)*(double(jp)-j) + dkp*dkp;
              setPixelPnormalPdepthPdataPauxPzoneBits(IJ,ij,normal,depth,1.0-(paux-d2_min)/(d2_max-d2_min),paux,i_vv[ivv],bits|INTERNAL_DATA_BIT|MESH_BIT|AUX_BIT|MASK_BIT);
              //setPixelPnormalPdepthPdataPauxPzoneBits(IJ,ij,normal,depth,sqrt(paux),paux,i_vv[ivv],bits|INTERNAL_DATA_BIT|MESH_BIT|AUX_BIT|MASK_BIT);
            }
          }
	}
      }
    }
  }

}

void CtiCanvas::setPixelData(const int i,const int j,const float normal[3],const float depth,const float data,const uint2 zone) {
  assert((i >= 0)&&(i < ni));
  assert((j >= 0)&&(j < nj));
  // assert(zone < 65535); // last zone is reserved
  // determine the index of the PixelBlock...
  const int IJ = (j>>3)*nI + (i>>3);
  assert((IJ >= 0)&&(IJ < nI*nJ));

  if (pbd[IJ] == NULL) pbd[IJ] = new PixelBlockData();

  const int ij = (j&7)*8+(i&7);
  assert((ij >= 0)&&(ij < 64));
  if ((pbd[IJ]->szone[ij] == 65535) || (depth < pbd[IJ]->sdepth[ij])) { // let depth always win
    //if ((pbd[IJ]->szone[ij] == 65535) || ((depth < pbd[IJ]->sdepth[ij])&&(pbd[IJ]->szone[ij] < 32768)) || (zone >= 32768) ) { // let edge always win
    pbd[IJ]->snormal[ij][0] = normal[0];
    pbd[IJ]->snormal[ij][1] = normal[1];
    pbd[IJ]->snormal[ij][2] = normal[2];
    pbd[IJ]->sdepth[ij] = depth;
    pbd[IJ]->pdata[ij] = data;
    pbd[IJ]->szone[ij] = zone;
  }
}

bool CtiCanvas::getPixelData(float normal[3],float& depth,float& data,uint2& zone,const int i,const int j) const {
  assert((i >= 0)&&(i < ni));
  assert((j >= 0)&&(j < nj));
  // assert(zone < 65535); // last zone is reserved
  // determine the index of the PixelBlock...
  const int IJ = (j>>3)*nI + (i>>3);
  assert((IJ >= 0)&&(IJ < nI*nJ));

  if (pbd[IJ] == NULL) return false;

  const int ij = (j&7)*8+(i&7);
  assert((ij >= 0)&&(ij < 64));

  normal[0] = pbd[IJ]->snormal[ij][0];
  normal[1] = pbd[IJ]->snormal[ij][1];
  normal[2] = pbd[IJ]->snormal[ij][2];

  depth = pbd[IJ]->sdepth[ij];

  data = pbd[IJ]->pdata[ij];

  zone = pbd[IJ]->szone[ij];

  return true;

}

  void interpFromData() {

    if (checkParam("INTERP_FROM_RESTART")) {

      // preallocate
      int * send_count = new int[mpi_size];
      int * send_disp = new int[mpi_size];
      int * recv_count = new int[mpi_size];
      int * recv_disp = new int[mpi_size];
      double * cv_wgt = new double[ncv];
      int * fa_flag = new int[nfa];
      double * fa_dn_data = new double[nfa];
      double (*fa_dn3_data)[3] = new double[nfa][3];
      bool * cv_interp = new bool[ncv]; FOR_ICV cv_interp[icv] = false; // flag if interpolated

      stringstream ss_hash; //build a string to identify the interpolated data in a hash

      int ndc = 0;
      FOR_PARAM_MATCHING("INTERP_FROM_RESTART") ndc++;

      int idc = 0;
      bool b_reset = false; // reset time and step to 0 (true if any INTERP_FROM_RESTART requests it)
      FOR_PARAM_MATCHING("INTERP_FROM_RESTART") {

        // put the data into the data container...

        int iarg = 0;
        const string filename = param->getString(iarg++);

        // ===============================================================
        // we can pass hints to the data container. Use the currently
        // registered data...
        // ===============================================================

        vector<string> nameVec;

        CtiRegister::setRegisteredCvDnAndDn3Names(nameVec);
        set<string> varSet;
        for (int ii = 0, ii_end=nameVec.size(); ii < ii_end; ++ii)  {
          // skip registered coordinate data: it should never be interpolated!...
          if ((nameVec[ii] != "x_cv")&&(nameVec[ii] != "x_vv"))
	    varSet.insert(nameVec[ii]);
        }

        nameVec.clear();
        CtiRegister::setRegisteredDAndINames(nameVec);
        for (int ii = 0, ii_end=nameVec.size(); ii < ii_end; ++ii)  {
          //if (mpi_rank == 0) cout << "XXXXXXXXXX d and i: " << nameVec[ii] << endl;
          varSet.insert(nameVec[ii]);
        }

        // now build the data container...

        DataContainer dc;

        if (filename.find(".les") != string::npos) {
          //dc.initFromOldRestart(filename,"u,T,step,time,dt"); // another example of passing hints
          int ierr = dc.initFromOldRestart(filename,varSet);
          if (ierr != 0) {
            CERR("initFromOldRestart failed");
          }
        }
	else if (filename.find(".sles") != string::npos) {
          int ierr = dc.initFromSles(filename,varSet);
          if (ierr != 0) {
            CERR("DataContainer::initFromSles() failed");
          }
	  // ss_hash stuff
          ss_hash << RestartHashUtilities::slesInterpHash;
          RestartHashUtilities::slesInterpHash.clear();
	}
        else if (filename.find(".mles") != string::npos) {
	  // this version may not be necessary going forward. x_vv is now included in sles
          const string sles_filename = param->getString(iarg++);
	  assert(sles_filename.find(".sles") != string::npos);
          int ierr = dc.initFromRestart(filename,sles_filename,varSet);
          if (ierr != 0) {
            CERR("initFromRestart failed");
          }
          //add interpolated file hashes to ss_hash;
          ss_hash << RestartHashUtilities::mlesInterpHash << RestartHashUtilities::slesInterpHash;
          RestartHashUtilities::mlesInterpHash.clear();
          RestartHashUtilities::slesInterpHash.clear();
        }
        else if ((filename.find(".cas") != string::npos) || (filename.find(".msh") != string::npos)) {
          const string dat_filename = param->getString(iarg++);
          int ierr = dc.initFromFluent(filename,dat_filename,varSet);
          if (ierr != 0) {
            CERR("initFromFluent failed");
          }
        }
        else if (filename.find(".pbin") != string::npos) {
          const string filename2 = param->getString(iarg++);
          if (filename2.find(".dat") != string::npos) { // coult eventually include others
            int ierr = dc.initFromPbinFluent(filename,filename2,varSet);
            if (ierr != 0) {
              CERR("initFromPbinFluent failed");
            }
          }
          else {
            CERR("Unknown data file extension provided after pbin in INTERP_FROM_RESTART.");
          }
        }
        else {
          CERR("Unrecognized file extension in INTERP_FROM_RESTART");
        }

        // TODO pbin, fluent, ASCII

        // store hash in RestartHashUtilities::slesHash
        RestartHashUtilities::slesHash.clear();
        RestartHashUtilities::slesHash.init(ss_hash,RestartHashUtilities::sha1hashlength);

        int double_count = prepareInterp(dc);
        double_count += 3; // for the x_vv's

        // if the dc requires transform, process these in order...

        while (iarg < param->size()) {
          const string token = param->getString(iarg++);
          if (token == "ROTATE_X") {
            const double degrees = param->getDouble(iarg++);
            if (mpi_rank == 0) cout << " > ROTATE_X " << degrees << " degrees" << endl;
            dc.rotate_x(degrees);

            MiscUtils::dumpRange(dc.x_cv,dc.ncv,"dc.x_cv");

          }
          else if (token == "TRANSLATE") {
            const double dx = param->getDouble(iarg++);
            const double dy = param->getDouble(iarg++);
            const double dz = param->getDouble(iarg++);
            if (mpi_rank == 0) cout << " > TRANSLATE " << dx << " " << dy << " " << dz << endl;
            dc.translate(dx,dy,dz);

            MiscUtils::dumpRange(dc.x_cv,dc.ncv,"dc.x_cv");

          }
          else if (token == "RESET") {
            b_reset = true;
          }
          else {
            if (mpi_rank == 0) cout << "Warning: skipping unrecognized INTERP_FROM_RESTART token: " << token << endl;
          }
        }

        // build the bbox...

        // make sure we have the stuff we need...

        if (cvAdt == NULL) buildCvAdt();

        // ========================================
        // we now have everything we need...
        // ========================================

        FOR_RANK send_count[rank] = 0;

        vector<int> intVec;
        for (int icv_ = 0; icv_ < dc.ncv; ++icv_) {
          assert(intVec.empty());
          cvBboxAdt->buildListForPoint(intVec,dc.x_cv[icv_]);
          for (int ii = 0, ii_end=intVec.size(); ii < ii_end; ++ii) {
            const int rank = intVec[ii]; assert((rank >= 0)&&(rank < mpi_size));
            send_count[rank] += double_count;
          }
          intVec.clear();
        }

        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
        const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

        double * send_buf = new double[send_count_sum];
        for (int icv_ = 0; icv_ < dc.ncv; ++icv_) {
          assert(intVec.empty());
          cvBboxAdt->buildListForPoint(intVec,dc.x_cv[icv_]);
          for (int ii = 0, ii_end=intVec.size(); ii < ii_end; ++ii) {
            const int rank = intVec[ii]; assert((rank >= 0)&&(rank < mpi_size));
            send_buf[send_disp[rank]  ] = dc.x_cv[icv_][0];
            send_buf[send_disp[rank]+1] = dc.x_cv[icv_][1];
            send_buf[send_disp[rank]+2] = dc.x_cv[icv_][2];
            send_disp[rank] += 3;
            // then the doubles...
            for (list<DnData>::iterator iter = dc.dnList.begin(); iter != dc.dnList.end(); ++iter) {
              if (iter->ctiData != NULL) {
                send_buf[send_disp[rank]] = iter->data[icv_];
                send_disp[rank] += 1;
              }
            }
            for (list<Dn3Data>::iterator iter = dc.dn3List.begin(); iter != dc.dn3List.end(); ++iter) {
              if (iter->ctiData != NULL) {
                send_buf[send_disp[rank]  ] = iter->data[icv_][0];
                send_buf[send_disp[rank]+1] = iter->data[icv_][1];
                send_buf[send_disp[rank]+2] = iter->data[icv_][2];
                send_disp[rank] += 3;
              }
            }
          }
          intVec.clear();
        }

        // rewind...

        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

        // now send...

        MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

        recv_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
        const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

        double * recv_buf = new double[recv_count_sum];
        MPI_Alltoallv(send_buf,send_count,send_disp,MPI_DOUBLE,
                      recv_buf,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
        delete[] send_buf;

        // now unpack and set...

        FOR_ICV cv_wgt[icv] = 0.0;

        double my_dist_max = 0.0;
        for (int irecv = 0; irecv < recv_count_sum; irecv += double_count) {
          double xp[3]; FOR_I3 xp[i] = recv_buf[irecv+i];
          assert(intVec.empty());
          cvAdt->buildListForPoint(intVec,xp);
          for (int ii = 0, ii_end=intVec.size(); ii < ii_end; ++ii) {
            const int icv = intVec[ii]; assert((icv >= 0)&&(icv < ncv));
            // here we need to check if the point is inside the cv. For now, we use
            // point proximity...
            const double dist = DIST(xp,x_vv[icv]);
            if (dist <= r_vv[icv]) { // HACK -- adjust tol here for exact point match
              my_dist_max = max(my_dist_max,dist/r_vv[icv]); // remember the normalized distance
              const double wgt = 1.0; // TODO: need something like r_vv[icv]/dist here
              cv_wgt[icv] += wgt;
              int i = 3;
              for (list<DnData>::iterator iter = dc.dnList.begin(); iter != dc.dnList.end(); ++iter) {
                if (iter->ctiData != NULL) {
                  iter->ctiData->dn(icv) += wgt*recv_buf[irecv+i];
                  i += 1;
                }
              }
              for (list<Dn3Data>::iterator iter = dc.dn3List.begin(); iter != dc.dn3List.end(); ++iter) {
                if (iter->ctiData != NULL) {
                  iter->ctiData->dn3(icv,0) += wgt*recv_buf[irecv+i  ];
                  iter->ctiData->dn3(icv,1) += wgt*recv_buf[irecv+i+1];
                  iter->ctiData->dn3(icv,2) += wgt*recv_buf[irecv+i+2];
                  i += 3;
                }
              }
              assert(i == double_count);
            }
          }
          intVec.clear();
        }

        delete[] recv_buf;

        // check normalization and prox...

        double dist_max;
        MPI_Reduce(&my_dist_max,&dist_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
        if (mpi_rank == 0)
          cout << " > dist_max: " << dist_max << endl;

        int8 my_count[2] = { ncv, 0 };
        FOR_ICV if (cv_wgt[icv] != 0.0) ++my_count[1];

        int8 count[2];
        MPI_Reduce(my_count,count,2,MPI_INT8,MPI_SUM,0,mpi_comm);
        if (mpi_rank == 0)
          cout << " > set " << count[1] << " out of " << count[0] << " cvs." << endl;

        // normalize...

        for (list<DnData>::iterator iter = dc.dnList.begin(); iter != dc.dnList.end(); ++iter) {
          if (iter->ctiData != NULL) {
            FOR_ICV {
              if (cv_wgt[icv] != 0.0) {
                iter->ctiData->dn(icv) /= cv_wgt[icv];
              }
              else {
                assert(cv_wgt[icv] == 0.0);
                assert(iter->ctiData->dn(icv) == 0.0);
              }
            }
          }
        }
        for (list<Dn3Data>::iterator iter = dc.dn3List.begin(); iter != dc.dn3List.end(); ++iter) {
          if (iter->ctiData != NULL) {
            FOR_ICV {
              if (cv_wgt[icv] != 0.0) {
                FOR_I3 iter->ctiData->dn3(icv,i) /= cv_wgt[icv];
              }
              else {
                assert(cv_wgt[icv] == 0.0);
                FOR_I3 assert(iter->ctiData->dn3(icv,i) == 0.0);
              }
            }
          }
        }

        // keep track of interpolated cells from all data sources
        FOR_ICV if (cv_wgt[icv] != 0.0) cv_interp[icv] = true;

        // only march after finishing all interpolations
        if (idc == (ndc-1) ) {

          // now use a marching algorithm to set any vars still unset...
          // here we use a face-based apporach to handle inter-processor data to allow
          // the registered CtiData to have arbitrary stride...

          if (mpi_rank == 0) cout << " > marching data:" << endl;

          // data set by the above process gets -1...

          FOR_ICV if (cv_interp[icv]) cv_wgt[icv] = -1.0;
          delete[] cv_interp;

          {

            for (list<DnData>::iterator iter = dc.dnList.begin(); iter != dc.dnList.end(); ++iter) {
              if (iter->ctiData != NULL) {

                if (mpi_rank == 0) cout << " > " << iter->name << "...";

                // reset cv_wgt...

                FOR_ICV if (cv_wgt[icv] != -1.0) cv_wgt[icv] = 0.0;

                int done = 0;
                while (done == 0) {

                  if (mpi_rank == 0) {
                    cout << ".";
                    cout.flush();
                  }

                  int my_done = 1;

                  FOR_INTERPROC_IFA {
                    const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
                    if (cv_wgt[icv0] <= -1.0) {
                      fa_flag[ifa] = 1;
                      fa_dn_data[ifa] = iter->ctiData->dn(icv0);
                    }
                    else {
                      assert(cv_wgt[icv0] == 0.0);
                      fa_flag[ifa] = 0;
                      fa_dn_data[ifa] = 0.0;
                    }
                  }

                  updateFaData(fa_flag,REPLACE_DATA);
                  updateFaData(fa_dn_data,REPLACE_DATA);

                  FOR_INTERPROC_IFA {
                    const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
                    if (cv_wgt[icv0] >= 0.0) {
                      if (fa_flag[ifa] == 1) {
                        // this nbr is good...
                        const double wgt = MAG(n_fa[ifa]); assert(wgt > 0.0);
                        cv_wgt[icv0] += wgt;
                        iter->ctiData->dn(icv0) += wgt*fa_dn_data[ifa];
                      }
                    }
                  }

                  FOR_INTERNAL_IFA {
                    const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
                    const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
                    if ((cv_wgt[icv0] >= 0.0)&&(cv_wgt[icv1] <= -1.0)) {
                      const double wgt = MAG(n_fa[ifa]); assert(wgt > 0.0);
                      cv_wgt[icv0] += wgt;
                      iter->ctiData->dn(icv0) += wgt*iter->ctiData->dn(icv1);
                    }
                    else if ((cv_wgt[icv0] <= -1.0)&&(cv_wgt[icv1] >= 0.0)) {
                      const double wgt = MAG(n_fa[ifa]); assert(wgt > 0.0);
                      cv_wgt[icv1] += wgt;
                      iter->ctiData->dn(icv1) += wgt*iter->ctiData->dn(icv0);
                    }
                  }

                  FOR_ICV {
                    if (cv_wgt[icv] > 0.0) {
                      iter->ctiData->dn(icv) /= cv_wgt[icv];
                      cv_wgt[icv] = -2.0; // use -2 to indicate a value that has been set by the marching
                    }
                    else if (cv_wgt[icv] == 0.0) {
                      // any zeros left indicate we are not done yet...
                      my_done = 0;
                    }
                    else {
                      assert(cv_wgt[icv] <= -1.0);
                    }
                  }

                  MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);

                }

                if (mpi_rank == 0) cout << "OK" << endl;

              }
            }

          }

          {

            for (list<Dn3Data>::iterator iter = dc.dn3List.begin(); iter != dc.dn3List.end(); ++iter) {
              if (iter->ctiData != NULL) {

                if (mpi_rank == 0) cout << " > " << iter->name << "...";

                // reset cv_wgt...

                FOR_ICV if (cv_wgt[icv] != -1.0) cv_wgt[icv] = 0.0;

                int done = 0;
                while (done == 0) {

                  if (mpi_rank == 0) {
                    cout << ".";
                    cout.flush();
                  }

                  int my_done = 1;

                  FOR_INTERPROC_IFA {
                    const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
                    if (cv_wgt[icv0] <= -1.0) {
                      fa_flag[ifa] = 1;
                      FOR_I3 fa_dn3_data[ifa][i] = iter->ctiData->dn3(icv0,i);
                    }
                    else {
                      fa_flag[ifa] = 0;
                      FOR_I3 fa_dn3_data[ifa][i] = 0.0;
                    }
                  }

                  updateFaData(fa_flag,REPLACE_DATA);
                  updateFaData(fa_dn3_data,REPLACE_ROTATE_DATA);

                  FOR_INTERPROC_IFA {
                    const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
                    if (cv_wgt[icv0] >= 0.0) {
                      if (fa_flag[ifa] == 1) {
                        // this nbr is good...
                        const double wgt = MAG(n_fa[ifa]); assert(wgt > 0.0);
                        cv_wgt[icv0] += wgt;
                        FOR_I3 iter->ctiData->dn3(icv0,i) += wgt*fa_dn3_data[ifa][i];
                      }
                    }
                  }

                  FOR_INTERNAL_IFA {
                    const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
                    const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
                    if ((cv_wgt[icv0] >= 0.0)&&(cv_wgt[icv1] <= -1.0)) {
                      const double wgt = MAG(n_fa[ifa]); assert(wgt > 0.0);
                      cv_wgt[icv0] += wgt;
                      FOR_I3 iter->ctiData->dn3(icv0,i) += wgt*iter->ctiData->dn3(icv1,i);
                    }
                    else if ((cv_wgt[icv0] <= -1.0)&&(cv_wgt[icv1] >= 0.0)) {
                      const double wgt = MAG(n_fa[ifa]); assert(wgt > 0.0);
                      cv_wgt[icv1] += wgt;
                      FOR_I3 iter->ctiData->dn3(icv1,i) += wgt*iter->ctiData->dn3(icv0,i);
                    }
                  }

                  FOR_ICV {
                    if (cv_wgt[icv] > 0.0) {
                      FOR_I3 iter->ctiData->dn3(icv,i) /= cv_wgt[icv];
                      cv_wgt[icv] = -2.0; // use -2 to indicate a value that has been set by the marching
                    }
                    else if (cv_wgt[icv] == 0.0) {
                      my_done = 0;
                    }
                    else {
                      assert(cv_wgt[icv] <= -1.0);
                    }
                  }

                  MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);

                }

                if (mpi_rank == 0) cout << "OK" << endl;

              }
            }
          }
        }
        ++idc;
      }

      // reset time and step if requested. we put this here because if the user only
      // reset the time on one of the INTERP_FROM_RESTART calls, we still want to
      // reset the time and step.
      if (b_reset) {
        COUT1(" > reseting time and step to 0");
        if ( CtiRegister::CtiData* data = CtiRegister::getRegisteredCtiData("step") ) {
          data->i() = 0.0;
        }
        if ( CtiRegister::CtiData* data = CtiRegister::getRegisteredCtiData("time") ) {
          data->d() = 0.0;
        }
      }

      delete[] fa_dn3_data;
      delete[] fa_dn_data;
      delete[] fa_flag;
      delete[] cv_wgt;
      delete[] recv_disp;
      delete[] recv_count;
      delete[] send_disp;
      delete[] send_count;
    }
  }
