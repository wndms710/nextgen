#include "CTI.hpp"
using namespace CTI;

#include "SurfaceShm.hpp"

// i,j,k can be stored using successive shifts of 21 bits...
#define TWENTYONE_BIT_MASK 2097151

void splitAndRedist(vector<pair<uint8,int> >& cubeTriPairVec,MPI_Comm &mpi_split_comm) {

  // get the rank and size in the mpi_split_comm...
  int mpi_split_rank,mpi_split_size;
  MPI_Comm_rank(mpi_split_comm, &mpi_split_rank);
  MPI_Comm_size(mpi_split_comm, &mpi_split_size);
  assert(mpi_split_size > 1);

  // get the range in each direction...
  int my_bbox[6] = { TWENTYONE_BIT_MASK, TWENTYONE_BIT_MASK, TWENTYONE_BIT_MASK, TWENTYONE_BIT_MASK, TWENTYONE_BIT_MASK, TWENTYONE_BIT_MASK };
  const int size = cubeTriPairVec.size();
  for (int ii = 0; ii < size; ++ii) {
    const uint8 ijk = cubeTriPairVec[ii].first;
    const int i = (ijk&TWENTYONE_BIT_MASK);
    my_bbox[0] = min(my_bbox[0],i);
    my_bbox[3] = min(my_bbox[3],-i);
    const int j = ((ijk>>21)&TWENTYONE_BIT_MASK);
    my_bbox[1] = min(my_bbox[1],j);
    my_bbox[4] = min(my_bbox[4],-j);
    const int k = ((ijk>>42)&TWENTYONE_BIT_MASK);
    my_bbox[2] = min(my_bbox[2],k);
    my_bbox[5] = min(my_bbox[5],-k);
  }
  int bbox[6];
  MPI_Allreduce(my_bbox,bbox,6,MPI_INT,MPI_MIN,mpi_split_comm);

  // which direction is the largest?...
  int id;
  if ( (-bbox[3]-bbox[0]) >= max(-bbox[4]-bbox[1],-bbox[5]-bbox[2]) ) {
    id = 0;
  }
  else if ( (-bbox[4]-bbox[1]) >= max(-bbox[3]-bbox[0],-bbox[5]-bbox[2]) ) {
    id = 1;
  }
  else {
    id = 2;
  }

  // collectively determine the best split point xmid. Note that everyone
  // equal to xmid gets included in the lower bin... 
  const int mpi_split_size_half = mpi_split_size/2;
  int8 my_count[2],count[2];
  vector<pair<int,int> > diVec(size);
  for (int ii = 0; ii < size; ++ii) {
    const uint8 ijk = cubeTriPairVec[ii].first;
    diVec[ii].first = ((ijk>>(21*id))&TWENTYONE_BIT_MASK);
    diVec[ii].second = ii;
  }
  sort(diVec.begin(),diVec.end());
  int x0 = bbox[id];
  int x1 = -bbox[id+3];
  int i0 = 0;
  int i1 = size-1;
  while (x0+1 < x1) {
    const int xmid = (x0+x1)/2;
    if (mpi_rank == 0) cout << "x0,x1,xmid: " << x0 << " " << x1 << " " << xmid << endl;
    int i0_ = i0;
    int i1_ = i1;
    if (size == 0) {
      my_count[0] = 0;
      my_count[1] = 0;
    }
    else if (diVec[size-1].first <= xmid) {
      my_count[0] = size;
      my_count[1] = 0;
    }
    else if (!(diVec[0].first <= xmid)) {
      my_count[0] = 0;
      my_count[1] = size;
    }
    else {
      while (i0_ < i1_-1) {
	const int imid = (i0_+i1_)/2;
	if (diVec[imid].first <= xmid)
	  i0_ = imid;
	else
	  i1_ = imid;
      }
      // advance to the xmid boundary if necessary...
      while ((i0_+1 < size)&&(diVec[i0_+1].first == xmid)) {
	++i0_;
      }
      my_count[0] = i0_+1;
      my_count[1] = size-i0_-1;
    }
    MPI_Allreduce(my_count,count,2,MPI_INT8,MPI_SUM,mpi_split_comm);
    if (count[1]*mpi_split_size_half > count[0]*(mpi_split_size-mpi_split_size_half)) {
      x0 = xmid;
      i0 = i0_;
    }
    else {
      x1 = xmid;
      i1 = i1_;
    }
  }
  cout << "got x0: " << x0 << endl;

  MPI_Pause("OKOKKOKOKO");

  // TODO: complete this some day 



#ifdef JFDUISFDSFD
  
    // we now know how to split the data...
    // we will redistribute with individual messages rather than all2allv. start by
    // exchanging everyone's counts and we all compute what messages will be going where...

    int8 (*my_count_of_rank)[2] = new int8[mpi_split_size][2];
    MPI_Allgather(my_count,2,MPI_INT8,(int8*)my_count_of_rank,2,MPI_INT8,mpi_split_comm);

    // the target number of cvs per core is...

    int8 my_count_target[2];
    // allowing for load imbalance, we need different targets for 0 and 1 processor groups...
    if (count[0]%mpi_split_size_half == 0)
      my_count_target[0] = count[0]/mpi_split_size_half;
    else
      my_count_target[0] = count[0]/mpi_split_size_half + 1;

    if (count[1]%(mpi_split_size-mpi_split_size_half) == 0)
      my_count_target[1] = count[1]/(mpi_split_size-mpi_split_size_half);
    else
      my_count_target[1] = count[1]/(mpi_split_size-mpi_split_size_half) + 1;

    assert(int(diVec.size()) == ncv);
    double (*x_cv_new)[3] = NULL;
    int8 *icv_global_new = NULL;
    int icv_new = 0;
    int icv_old = 0;
    if (mpi_split_rank < mpi_split_size_half) {
      // this rank is part of the first "half"...
      const int n_new = max(my_count[0],my_count_target[0]);
      x_cv_new = new double[n_new][3];
      icv_global_new = new int8[n_new];
      // simultaneously copy values that are staying into _new arrays, and compress old...
      for (int ii = 0, ii_max=diVec.size(); ii < ii_max; ++ii) {
	const int icv = diVec[ii].second;
	assert(icv_global[icv] >= 0);
	if (diVec[ii].first <= xmid) {
	  x_cv_new[icv_new][0] = x_cv[icv][0];
	  x_cv_new[icv_new][1] = x_cv[icv][1];
	  x_cv_new[icv_new][2] = x_cv[icv][2];
	  icv_global_new[icv_new] = icv_global[icv];
	  icv_global[icv] = -1;
	  ++icv_new;
	}
      }
      for (int icv = 0; icv < ncv; ++icv) {
	if (icv_global[icv] >= 0) {
	  if (icv_old < icv) {
	    x_cv[icv_old][0] = x_cv[icv][0];
	    x_cv[icv_old][1] = x_cv[icv][1];
	    x_cv[icv_old][2] = x_cv[icv][2];
	    icv_global[icv_old] = icv_global[icv];
	  }
	  ++icv_old;
	}
      }
      assert(icv_new == my_count[0]);
      assert(icv_old == my_count[1]);
    }
    else {
      // this rank is part of the second "half"...
      const int n_new = max(my_count[1],my_count_target[1]);
      x_cv_new = new double[n_new][3];
      icv_global_new = new int8[n_new];
      // same as above, but flipped...
      for (int ii = 0,ii_max = diVec.size(); ii < ii_max; ++ii) {
	const int icv = diVec[ii].second;
	assert(icv_global[icv] >= 0);
	if (!(diVec[ii].first <= xmid)) {
	  x_cv_new[icv_new][0] = x_cv[icv][0];
	  x_cv_new[icv_new][1] = x_cv[icv][1];
	  x_cv_new[icv_new][2] = x_cv[icv][2];
	  icv_global_new[icv_new] = icv_global[icv];
	  icv_global[icv] = -1;
	  ++icv_new;
	}
      }
      for (int icv = 0; icv < ncv; ++icv) {
	if (icv_global[icv] >= 0) {
	  if (icv_old < icv) {
	    x_cv[icv_old][0] = x_cv[icv][0];
	    x_cv[icv_old][1] = x_cv[icv][1];
	    x_cv[icv_old][2] = x_cv[icv][2];
	    icv_global[icv_old] = icv_global[icv];
	  }
	  ++icv_old;
	}
      }
      assert(icv_new == my_count[1]);
      assert(icv_old == my_count[0]);
    }

    // reset the reference to the old buffers...
    icv_old = 0;

    // we now have the new array with new data in it up to icv_new, and the
    // old array with old data that needs to be sent somewhere in old...

    // turn the new data counter into a space indicator...

    for (int rank0 = 0; rank0 < mpi_split_size_half; ++rank0)
      my_count_of_rank[rank0][0] = max(0,int(my_count_target[0]-my_count_of_rank[rank0][0]));
    for (int rank1 = mpi_split_size_half; rank1 < mpi_split_size; ++rank1)
      my_count_of_rank[rank1][1] = max(0,int(my_count_target[1]-my_count_of_rank[rank1][1]));

    MPI_Request * sendRequestArray = NULL;
    MPI_Request * recvRequestArray = NULL;

    // loop twice because we need to count and then post the send/recv pairs...

    int send_count,recv_count;
    for (int iter = 0; iter < 2; ++iter) {

      int8 (*my_count_of_rank_copy)[2] = NULL;
      if (iter == 0) {
	my_count_of_rank_copy = new int8[mpi_split_size][2];
	for (int rank = 0; rank < mpi_split_size; ++rank) {
	  my_count_of_rank_copy[rank][0] = my_count_of_rank[rank][0];
	  my_count_of_rank_copy[rank][1] = my_count_of_rank[rank][1];
	}
      }

      send_count = 0;
      recv_count = 0;

      // start by sending the count[1] stuff from the lower half to the upper half...
      {
	int rank1 = mpi_split_size_half;
	for (int rank0 = 0; rank0 < mpi_split_size_half; ++rank0) {
	  while (my_count_of_rank[rank0][1] > 0) {
	    // we have some data that has to be sent to rank1...
	    //cout << mpi_rank << " rank1: " << rank1 << " mpi_split_size: " << mpi_split_size << " mpi_split_size_half: " << mpi_split_size_half << endl;
	    assert(rank1 < mpi_split_size);
	    if (my_count_of_rank[rank1][1] > 0) {
	      const int n_send_recv = min(my_count_of_rank[rank0][1],my_count_of_rank[rank1][1]);
	      // rank0 needs to send rank1 n_send_recv...
	      if (rank0 == mpi_split_rank) {
		if (iter == 1) {
		  // post sends...
		  MPI_Issend(x_cv+icv_old,n_send_recv*3,MPI_DOUBLE,rank1,
			     12345,mpi_split_comm,&(sendRequestArray[send_count]));
		  MPI_Issend(icv_global+icv_old,n_send_recv,MPI_INT8,rank1,
			     12346,mpi_split_comm,&(sendRequestArray[send_count+1]));
		  icv_old += n_send_recv;
		}
		send_count += 2;
	      }
	      else if (rank1 == mpi_split_rank) {
		if (iter == 1) {
		  // post recvs...
		  MPI_Irecv(x_cv_new+icv_new,n_send_recv*3,MPI_DOUBLE,rank0,
			    12345,mpi_split_comm,&(recvRequestArray[recv_count]));
		  MPI_Irecv(icv_global_new+icv_new,n_send_recv,MPI_INT8,rank0,
			    12346,mpi_split_comm,&(recvRequestArray[recv_count+1]));
		  icv_new += n_send_recv;
		}
		recv_count += 2;
	      }
	      my_count_of_rank[rank0][1] -= n_send_recv;
	      my_count_of_rank[rank1][1] -= n_send_recv;
	    }
	    if (my_count_of_rank[rank1][1] == 0)
	      ++rank1;
	  }
	}
      }

      // and then the count[0] stuff from the upper half to the lower half...
      {
	int rank0 = 0;
	for (int rank1 = mpi_split_size_half; rank1 < mpi_split_size; ++rank1) {
	  while (my_count_of_rank[rank1][0] > 0) {
	    // we have some data that has to be sent to rank0...
	    //cout << mpi_rank << " rank0: " << rank0 << " mpi_split_size: " << mpi_split_size << " mpi_split_size_half: " << mpi_split_size_half << endl;
	    assert(rank0 < mpi_split_size);
	    if (my_count_of_rank[rank0][0] > 0) {
	      const int n_send_recv = min(my_count_of_rank[rank1][0],my_count_of_rank[rank0][0]);
	      // rank1 needs to send rank0 n_send_recv...
	      if (rank1 == mpi_split_rank) {
		if (iter == 1) {
		  // post sends...
		  MPI_Issend(x_cv+icv_old,n_send_recv*3,MPI_DOUBLE,rank0,
			     12347,mpi_split_comm,&(sendRequestArray[send_count]));
		  MPI_Issend(icv_global+icv_old,n_send_recv,MPI_INT8,rank0,
			     12348,mpi_split_comm,&(sendRequestArray[send_count+1]));
		  icv_old += n_send_recv;
		}
		send_count += 2;
	      }
	      else if (rank0 == mpi_split_rank) {
		if (iter == 1) {
		  // post recv...
		  MPI_Irecv(x_cv_new+icv_new,n_send_recv*3,MPI_DOUBLE,rank1,
			    12347,mpi_split_comm,&(recvRequestArray[recv_count]));
		  MPI_Irecv(icv_global_new+icv_new,n_send_recv,MPI_INT8,rank1,
			    12348,mpi_split_comm,&(recvRequestArray[recv_count+1]));
		  icv_new += n_send_recv;
		}
		recv_count += 2;
	      }
	      my_count_of_rank[rank0][0] -= n_send_recv;
	      my_count_of_rank[rank1][0] -= n_send_recv;
	    }
	    if (my_count_of_rank[rank0][0] == 0)
	      ++rank0;
	  }
	}
      }

      // on the first time through reset my_count_of_rank and allocate the MPI_Request arrays...

      if (iter == 0) {

	for (int rank = 0; rank < mpi_split_size; ++rank) {
	  my_count_of_rank[rank][0] = my_count_of_rank_copy[rank][0];
	  my_count_of_rank[rank][1] = my_count_of_rank_copy[rank][1];
	}
	delete[] my_count_of_rank_copy;

	if (send_count > 0) sendRequestArray = new MPI_Request[send_count];
	if (recv_count > 0) recvRequestArray = new MPI_Request[recv_count];

      }

    }

    delete[] my_count_of_rank;

    // --------------------------------------------------------------------------------
    // ok -- everything is posted, so now just allow them to complete. Unpacking
    // is done directly into the relevant buffers...
    // --------------------------------------------------------------------------------

    if ((send_count > 0)||(recv_count > 0)) {

      MPI_Status * statusArray = new MPI_Status[max(send_count,recv_count)];

      if (recv_count > 0) {
	MPI_Waitall(recv_count,recvRequestArray,statusArray);
	delete[] recvRequestArray;
      }

      if (send_count > 0) {
	MPI_Waitall(send_count,sendRequestArray,statusArray);
	delete[] sendRequestArray;
      }

      delete[] statusArray;

    }

    if (x_cv != NULL) delete[] x_cv;
    x_cv = x_cv_new;

    if (icv_global != NULL) delete[] icv_global;
    icv_global = icv_global_new;

    ncv = icv_new;

    // wait for everyont to complete...

    MPI_Barrier(mpi_split_comm);

    // --------------------------------------------------------------------------------
    // now split the communicator before returning...
    // --------------------------------------------------------------------------------

    int mpi_key;
    if (mpi_split_rank < mpi_split_size_half)
      mpi_key = 0;
    else
      mpi_key = 1;

    MPI_Comm mpi_split_comm_copy;
    MPI_Comm_dup(mpi_split_comm,&mpi_split_comm_copy);
    MPI_Comm_free(&mpi_split_comm);
    MPI_Comm_split(mpi_split_comm_copy, mpi_key, mpi_split_rank, &mpi_split_comm);
    MPI_Comm_free(&mpi_split_comm_copy);

    // take a look...

    /*
      int mpi_split_rank_new,mpi_split_size_new;
      MPI_Comm_rank(mpi_split_comm, &mpi_split_rank_new);
      MPI_Comm_size(mpi_split_comm, &mpi_split_size_new);
      cout << "mpi_split_rank: " << mpi_split_rank << " mpi_split_size_new: " << mpi_split_size_new << " mpi_split_rank_new: " << mpi_split_rank_new << endl;
    */

#endif
    
}

class MySurfaceShm : public SurfaceShm {
  
public:
  
  MySurfaceShm(const string& name) : SurfaceShm(name) {

  }

  void buildLiftedSurfaceFlaggedTris(const double delta) {

    if (mpi_rank == 0) cout << "buildLiftedSurfaceFlaggedTris: delta=" << delta << endl;
    
    assert(flag_st);
    assert((1<<21) == TWENTYONE_BIT_MASK+1);

    // figure out the bbox...
    double bbminmax[6]; // stored as xmin,ymin,zmin,xmax,ymax,zmax
    {
      double my_bbminmax[6] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };
      for (int ist = mpi_rank; ist < nst; ist += mpi_size) {
	if (flag_st[ist]) {
	  FOR_I3 {
	    const int isp = spost[ist][i];
	    FOR_J3 {
	      my_bbminmax[j] = min(my_bbminmax[j],xsp[isp][j]);
	      my_bbminmax[3+j] = min(my_bbminmax[3+j],-xsp[isp][j]);
	    }
	  }
	}
      }
      MPI_Allreduce(my_bbminmax,bbminmax,6,MPI_DOUBLE,MPI_MIN,mpi_comm);
      FOR_I3 bbminmax[3+i] = -bbminmax[3+i];
    }

    // and the spacing dx...
    const double dx = delta*0.2; // this factor is important - it should be as large as possible

    if (mpi_rank == 0) {
      cout << " > bbox: x: " << bbminmax[0] << " " << bbminmax[3] << 
	" y: " << bbminmax[1] << " " << bbminmax[4] << 
	" z: " << bbminmax[2] << " " << bbminmax[5] << endl;
      cout << " > cartesian grid size: " << 
	(bbminmax[3]-bbminmax[0]+2.0*delta)/dx << " " << 
	(bbminmax[4]-bbminmax[1]+2.0*delta)/dx << " " << 
	(bbminmax[5]-bbminmax[2]+2.0*delta)/dx << endl;
    }

    // grid WSB corner: corresponds to a node [0,0,0]...
    double x0[3];
    FOR_I3 x0[i] = bbminmax[i]-delta-dx;

    // now build the cube-tri pairs that link a cube and 
    // a tri bounding box...
    vector<pair<uint8,int> > cubeTriPairVec;
    for (int ist = mpi_rank; ist < nst; ist += mpi_size) {
      if (flag_st[ist]) {
	int ijk_min[3],ijk_max[3];
	FOR_I3 {
	  ijk_min[i] = int((min(xsp[spost[ist][0]][i],min(xsp[spost[ist][1]][i],xsp[spost[ist][2]][i]))-x0[i])/dx);
	  assert(ijk_min[i] > 0);
	  ijk_max[i] = int((max(xsp[spost[ist][0]][i],max(xsp[spost[ist][1]][i],xsp[spost[ist][2]][i]))-x0[i])/dx);
	  assert(ijk_min[i] < (1<<21)); // use a 21 bit shift to store the i,j,k as a single uint8. If we ever hit this, we could modify 
	}
	for (int i = ijk_min[0]; i < ijk_max[0]; ++i) {
	  for (int j = ijk_min[1]; j < ijk_max[1]; ++j) {
	    for (int k = ijk_min[2]; k < ijk_max[2]; ++k) {
	      cubeTriPairVec.push_back(pair<uint8,int>(uint8(i)|(uint8(j)<<21)|(uint8(k)<<42),ist));
	    }
	  }
	}
      }
    }

    // now load balance these cube/tri interactions...
    int8 my_count = cubeTriPairVec.size();
    int8 count;
    MPI_Reduce(&my_count,&count,1,MPI_INT8,MPI_SUM,0,mpi_comm);
    int8 count_max;
    MPI_Reduce(&my_count,&count_max,1,MPI_INT8,MPI_MAX,0,mpi_comm);
    if (mpi_rank == 0) {
      cout << " > cube-tri pair load balance: avg: " << double(count)/double(mpi_size) << " max: " << count_max << endl;
    }
 
    int mpi_split_size;
    MPI_Comm_size(mpi_comm, &mpi_split_size);
    MPI_Comm mpi_split_comm; MPI_Comm_dup(mpi_comm,&mpi_split_comm);
    while (mpi_split_size > 1) {
      // this routine returns a split cubeTriPairVec and mpi_split_comm modified to continue...
      splitAndRedist(cubeTriPairVec,mpi_split_comm); // TODO: not done yet: see above
      // and get the size again...
      MPI_Comm_size(mpi_split_comm,&mpi_split_size);
    }
    MPI_Comm_free(&mpi_split_comm);
    
    // TODO: now we need to determine the corners associated with the initial condition, and
    // also add a suitable number of surrounding unknowns. This is then load balanced 
    // AGAIN to produce the final grid where the fast-sweeping method is performed and the
    // final surface produced by marching cubes. How can this be fast?...

    // could build a 3d grid with sufficient ghost right now based on the load-balanced cubes,
    // with an appropriate amount of ghost data...
    
    // now convert to initial condition points by evaluating closest distance to tri...
    for (int ii = 0, size =  cubeTriPairVec.size(); ii < size; ++ii) {
      const uint8 ijk = cubeTriPairVec[ii].first;
      const int ist = cubeTriPairVec[ii].second;
      const int i = (ijk&TWENTYONE_BIT_MASK);
      const int j = ((ijk>>21)&TWENTYONE_BIT_MASK);
      const int k = ((ijk>>42)&TWENTYONE_BIT_MASK);
      double x[3];
      x[0] = x0[0] + dx*i;
      x[1] = x0[1] + dx*j;
      x[2] = x0[2] + dx*k;
      
      
    }
    
    
    MPI_Pause("OKOK");


  }
  
};

int main(int argc, char * argv[]) {

  try {
      
    CTI_Init(argc,argv,"testSurfaceShm.in");
    
    {
      
      Param * param = getParam("SURF");
      if (!param) {
	CERR("missing SURF");
      }
      assert(param->size() >= 1);
      assert(param->getString(0) == "SBIN");
      string filename = param->getString(1);
      MySurfaceShm * surface = new MySurfaceShm(filename);

      
      // NAME=rr_rh ./rotating_rr_rh_180.mles MOTION ROTATION XC 2.65895327 0.89150566 0.35316522 AXIS 0.00168939 -0.99969830 -0.02450393 RPT 11.942114

      double xc[3] = { 2.65895327, 0.89150566, 0.35316522 };
      double axis[3] = { 0.00168939, -0.99969830, -0.02450393 }; 
      surface->applyInterface("farfield_rr_rh",axis,xc);
      surface->writeSbin("test.sbin");
      MPI_Pause("have a look!");
      assert(0);
      
      // try the znozn graph...
      //surface->ensureZnozn();

      // is flag_st available...
      surface->ensureFlagSt();
      
      // set the flag in certain tris...
      if (mpi_rank == 0) {
	int count = 0;
	for (int ist = 0; ist < surface->nst; ++ist) {
	  surface->flag_st[ist] = 0;
	  if (rand() < RAND_MAX/2) {
	    surface->flag_st[ist] = 1;
	    ++count;
	  }
	}
	cout << " > flag set in " << count << " tris out of " << surface->nst << endl;
      }
      if (mpi_rank_shared == 0) {
	MPI_Bcast(surface->flag_st,surface->nst,MPI_INT,0,mpi_comm_internode);
      }
      MPI_Barrier(mpi_comm_shared);

      surface->buildLiftedSurfaceFlaggedTris(0.01);
      
      MPI_Pause("XXXXXXXXXXXXXXXXXXXx HOW WAS THA? XXXXXXXXXXXXXXXXXXXXXXXXXX");
	
      delete surface;
      
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

