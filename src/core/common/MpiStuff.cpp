#include "MpiStuff.hpp"
#include "Macros.hpp" // for MPI_TYPE_FREE
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <assert.h>
#include <cstring>

using std::cout;
using std::cerr;
using std::endl;
using std::max;

#include <sys/stat.h> // for stat fileExists

#if defined(BGP_MEMCHECK) || defined(BGQ_MEMCHECK)
#include "bluegene/memcheck2.hpp"
#elif defined(LINUX_MEMCHECK)
#include "linux/memcheck.hpp"
#endif

namespace MpiStuffHelper {

  int part_kind = 0; // 0 == regular PADT, 1 == BAD

  // these routines are used in the padt partitioning...

  void calcUniformDistNew(int8 * &xod, const int8 n, const int ndist) {

    assert(xod == NULL);
    xod = new int8[ndist+1];
    xod[0] = 0;

    int8 remainder = n%ndist;
    if (remainder == 0) {
      for (int id = 1; id <= ndist; id++) {
	xod[id] = xod[id-1] + n/ndist;
      }
    }
    else {
      for (int id = 1; id <= ndist; id++) {
	xod[id] = min(n,xod[id-1] + n/ndist + 1);
      }
    }
    assert(xod[ndist] == n);

    /*
      int8 n_max = 0;
      for (int id = 0; id < ndist; ++id)
      n_max = max(n_max,xod[id+1]-xod[id]);
      if (n_max > TWO_BILLION) {
      if (mpi_rank == 0)
      cout << "Error: max n per core exceeds 2 billion: " << n_max << ". This size problem requires ndist (np?) >= " << ceil(double(ndist)*double(n_max)/2.0E+9) << "." << endl;
      throw(0);
      }
    */

  }

  int getRankInXoraNew(const int8 ix,const int8 * xora,const int size) {

    assert((ix >= 0)&&(ix < xora[size]));

    // bracket...
    int left_rank = 0;
    int right_rank = size;
    while ((right_rank - left_rank) > 1) {
      const int middle_rank = (left_rank + right_rank)/2;   // equivalent to floor..
      if (ix >= xora[middle_rank])
	left_rank = middle_rank;
      else
	right_rank = middle_rank;
    }
    return(left_rank);

  }

  void splitAndRedist(double (* &x_cv)[3],int8 * &icv_global,int &ncv,MPI_Comm &mpi_split_comm,const int bits) {

    //int mpi_rank;
    //MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    // this routine geometrically splits the passed arrays between the processors in the passed
    // mpi_split_comm. Assume the x_cv, icv_global, ncv, and mpi_split_comm are all changed on return...

    // get the rank and size in the mpi_split_comm...

    int mpi_split_rank,mpi_split_size;
    MPI_Comm_rank(mpi_split_comm, &mpi_split_rank);
    MPI_Comm_size(mpi_split_comm, &mpi_split_size);
    assert(mpi_split_size > 1);

    // get the range in each direction...

    double my_bbox[6] = {HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL};
    FOR_ICV {
      FOR_I3 my_bbox[i] = min(my_bbox[i],x_cv[icv][i]);
      FOR_I3 my_bbox[i+3] = min(my_bbox[i+3],-x_cv[icv][i]);
    }
    double bbox[6];
    MPI_Allreduce(my_bbox,bbox,6,MPI_DOUBLE,MPI_MIN,mpi_split_comm);

    // which direction is the largest?...

    int i;
    if (bits == 1) {
      i = 0;
    }
    else if (bits == 2) {
      i = 1;
    }
    else if (bits == 3) {
      if ( (-bbox[3]-bbox[0]) >= (-bbox[4]-bbox[1]) ) {
        i = 0;
      }
      else {
        i = 1;
      }
    }
    else if (bits == 4) {
      i = 2;
    }
    else if (bits == 5) {
      if ( (-bbox[3]-bbox[0]) >= (-bbox[5]-bbox[2]) ) {
        i = 0;
      }
      else {
        i = 2;
      }
    }
    else if (bits == 6) {
      if ( (-bbox[4]-bbox[1]) >= (-bbox[5]-bbox[2]) ) {
        i = 0;
      }
      else {
        i = 2;
      }
    }
    else { 
      assert(bits == 7);
      if ( (-bbox[3]-bbox[0]) >= max(-bbox[4]-bbox[1],-bbox[5]-bbox[2]) ) {
        i = 0;
      }
      else if ( (-bbox[4]-bbox[1]) >= max(-bbox[3]-bbox[0],-bbox[5]-bbox[2]) ) {
        i = 1;
      }
      else {
        i = 2;
      }
    }

    //i = 0; // HACK -- fix split direction for testing...

    // sort the data in this direction...

    const double eps1 = 1.1234E-6; // include a little y and z to break sort ties (and hopefully not cause any)...
    const double eps2 = 2.7531E-6; // include a little y and z to break sort ties (and hopefully not cause any)...

    // the same values caused problems with y = -z points...
    //const double eps1 = 1.0E-6;
    //const double eps2 = 1.0E-6;

    // reduce scope for diVec...

    const int mpi_split_size_half = mpi_split_size/2;
    double xmid;
    int8 my_count[2],count[2];

    vector<pair<double,int> > diVec(ncv);
    double x0,x1;

    switch (part_kind) {
    case 0: // PADT

      FOR_ICV {
	diVec[icv].first = x_cv[icv][i] + eps1*x_cv[icv][(i+1)%3] + eps2*x_cv[icv][(i+2)%3];
	diVec[icv].second = icv;
      }
      sort(diVec.begin(),diVec.end());
      x0 = bbox[i] + eps1*bbox[(i+1)%3] + eps2*bbox[(i+2)%3];
      x1 = -bbox[i+3] - eps1*bbox[(i+1)%3+3] - eps2*bbox[(i+2)%3+3]; // the bbox max is negative...
      break;

    case 1: // BAD

      FOR_ICV {
	diVec[icv].first = double(rand())/double(RAND_MAX);
	diVec[icv].second = icv;
      }
      sort(diVec.begin(),diVec.end());
      x0 = 0.0;
      x1 = 1.0;
      break;

    default:

      assert(0);

    }

    int i0 = 0;
    int i1 = ncv-1;
    int iter = 0;
    while (1) {

      xmid = 0.5*(x0+x1);

      int i0_ = i0;
      int i1_ = i1;
      if (ncv == 0) {
	my_count[0] = 0;
	my_count[1] = 0;
      }
      else if (diVec[ncv-1].first <= xmid) {
	my_count[0] = ncv;
	my_count[1] = 0;
      }
      else if (!(diVec[0].first <= xmid)) {
	my_count[0] = 0;
	my_count[1] = ncv;
      }
      else {
	while (i0_ < i1_-1) {
	  const int imid = (i0_+i1_)/2;
	  if (diVec[imid].first <= xmid)
	    i0_ = imid;
	  else
	    i1_ = imid;
	}
	my_count[0] = i0_+1;
	my_count[1] = ncv-i0_-1;
      }

      /*
	{
	int my_count_check[2] = { 0, 0 };
	FOR_ICV {
	if (x_cv[icv][i] + eps1*x_cv[icv][(i+1)%3] + eps2*x_cv[icv][(i+2)%3] <= xmid)
	my_count_check[0] += 1;
	else
	my_count_check[1] += 1;
	}
	//cout << mpi_rank << " my_count: " << my_count[0] << " " << my_count[1] << " my_count_check: " << my_count_check[0] << " " << my_count_check[1] << endl;
	//MPI_Pause("OKOK");
	assert(my_count[0] == my_count_check[0]);
	assert(my_count[1] == my_count_check[1]);
	//my_count[0] = my_count_check[0];
	//my_count[1] = my_count_check[1];
	}
      */

      MPI_Allreduce(my_count,count,2,MPI_INT8,MPI_SUM,mpi_split_comm); // should be MPI_INT8

      if ( (count[1]*mpi_split_size_half > count[0]*(mpi_split_size-mpi_split_size_half)) &&
	   ((count[1]-1)*mpi_split_size_half >= (count[0]+1)*(mpi_split_size-mpi_split_size_half)) ) {
	x0 = xmid;
	i0 = i0_;
      }
      else if ( (count[1]*mpi_split_size_half < count[0]*(mpi_split_size-mpi_split_size_half)) &&
		((count[1]+1)*mpi_split_size_half <= (count[0]-1)*(mpi_split_size-mpi_split_size_half)) ) {
	x1 = xmid;
	i1 = i1_;
      }
      else
	break;

      ++iter;
      if (iter > 50) {
	if (mpi_split_rank == 0)
	  cout << " Warning: 50 iters exceed in bisection routine: count[0]: " << count[0] << " count[1]: " << count[1] << " x0: " << x0 << " xmid: " << xmid << " x1: " << x1 << endl;
	if (iter > 75)
	  break;  // can't seem to balance. Just move ahead anyways...
      }

    }

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

  }

  void splitAndRedistWeighted(double (* &x_cv)[3],int8 (* &cv_global_wgt)[2],int &ncv,MPI_Comm &mpi_split_comm,const int bits) {

    assert(bits != 13); // 13 is used for a "bad" partion

    //int mpi_rank;
    //MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    // this routine geometrically splits the passed arrays between the processors in the passed
    // mpi_split_comm. Assume the x_cv, cv_global_wgt, ncv, and mpi_split_comm are all changed on return...

    // get the rank and size in the mpi_split_comm...

    int mpi_split_rank,mpi_split_size;
    MPI_Comm_rank(mpi_split_comm, &mpi_split_rank);
    MPI_Comm_size(mpi_split_comm, &mpi_split_size);
    assert(mpi_split_size > 1);

    // get the range in each direction...

    double my_bbox[6] = {HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL};
    FOR_ICV {
      FOR_I3 my_bbox[i] = min(my_bbox[i],x_cv[icv][i]);
      FOR_I3 my_bbox[i+3] = min(my_bbox[i+3],-x_cv[icv][i]);
    }
    double bbox[6];
    MPI_Allreduce(my_bbox,bbox,6,MPI_DOUBLE,MPI_MIN,mpi_split_comm);

    // which direction is the largest?...

    int idir;
    if (bits == 1) {
      idir = 0;
    }
    else if (bits == 2) {
      idir = 1;
    }
    else if (bits == 3) {
      if ( (-bbox[3]-bbox[0]) >= (-bbox[4]-bbox[1]) ) {
        idir = 0;
      }
      else {
        idir = 1;
      }
    }
    else if (bits == 4) {
      idir = 2;
    }
    else if (bits == 5) {
      if ( (-bbox[3]-bbox[0]) >= (-bbox[5]-bbox[2]) ) {
        idir = 0;
      }
      else {
        idir = 2;
      }
    }
    else if (bits == 6) {
      if ( (-bbox[4]-bbox[1]) >= (-bbox[5]-bbox[2]) ) {
        idir = 0;
      }
      else {
        idir = 2;
      }
    }
    else { 
      assert(bits == 7);
      if ( (-bbox[3]-bbox[0]) >= max(-bbox[4]-bbox[1],-bbox[5]-bbox[2]) ) {
        idir = 0;
      }
      else if ( (-bbox[4]-bbox[1]) >= max(-bbox[3]-bbox[0],-bbox[5]-bbox[2]) ) {
        idir = 1;
      }
      else {
        idir = 2;
      }
    }

    // sort the data in this direction...

    const double eps1 = 1.1234E-6; // include a little y and z to break sort ties (and hopefully not cause any)...
    const double eps2 = 2.7531E-6; // include a little y and z to break sort ties (and hopefully not cause any)...

    // the same values caused problems with y = -z points...
    //const double eps1 = 1.0E-6;
    //const double eps2 = 1.0E-6;

    // reduce scope for diVec...

    const int mpi_split_size_half = mpi_split_size/2;
    double xmid = 0.0; // init to avoid warning
    int8 my_count[4],count[4];

    {

      vector< pair<double,int8> > diVec(ncv);
      FOR_ICV {
	diVec[icv].first = x_cv[icv][idir] + eps1*x_cv[icv][(idir+1)%3] + eps2*x_cv[icv][(idir+2)%3];
	diVec[icv].second = cv_global_wgt[icv][1]; // the weight
      }
      sort(diVec.begin(),diVec.end());
      for (int icv = 1; icv < ncv; ++icv)
	diVec[icv].second += diVec[icv-1].second;

      double x0 = bbox[idir] + eps1*bbox[(idir+1)%3] + eps2*bbox[(idir+2)%3];
      double x1 = -bbox[idir+3] - eps1*bbox[(idir+1)%3+3] - eps2*bbox[(idir+2)%3+3]; // the bbox max is negative...

      int i0 = 0;
      int i1 = ncv-1;
      int iter = 0;
      int8 count0_x0 = -1;
      int8 count0_x1 = -1;
      while ((count0_x0 == -1)||(count0_x1 == -1)||(count0_x1 > count0_x0+1)) {

	xmid = 0.5*(x0+x1);

	int i0_ = i0;
	int i1_ = i1;
	if (ncv == 0) {
	  my_count[0] = 0;
	  my_count[1] = 0;
	  my_count[2] = 0; // the cv-count[0] -- used for testing convergence
	  my_count[3] = 0;
	}
	else if (diVec[ncv-1].first <= xmid) {
	  my_count[0] = diVec[ncv-1].second;
	  my_count[1] = 0;
	  // cv counts...
	  my_count[2] = ncv;
	  my_count[3] = 0;
	}
	else if (!(diVec[0].first <= xmid)) {
	  my_count[0] = 0;
	  my_count[1] = diVec[ncv-1].second;
	  // cv counts...
	  my_count[2] = 0;
	  my_count[3] = ncv;
	}
	else {
	  while (i0_ < i1_-1) {
	    const int imid = (i0_+i1_)/2;
	    if (diVec[imid].first <= xmid)
	      i0_ = imid;
	    else
	      i1_ = imid;
	  }
	  my_count[0] = diVec[i0_].second;
	  my_count[1] = diVec[ncv-1].second-diVec[i0_].second;
	  // cv counts...
	  my_count[2] = i0_+1;
	  my_count[3] = ncv-i0_-1;
	}

	MPI_Allreduce(my_count,count,4,MPI_INT8,MPI_SUM,mpi_split_comm);

	// it should never take more than 53ish iterations to saturate double-precision (53 bit mantissa)...

	++iter;
	if ((xmid == x0)||(xmid == x1)) {
	  break;
	}
	else if (iter > 64) {
	  if (mpi_split_rank == 0)
	    cout << " > Warning: 64 iters exceed in bisection routine: count[0]: " <<
	      count[0] << " count[1]: " << count[1] << " x0: " << x0 << " xmid: " << xmid << " x1: " << x1 << endl;
	  break;  // can't seem to balance. Warn then move ahead anyways...
	}

	if (count[1]*mpi_split_size_half > count[0]*(mpi_split_size-mpi_split_size_half)) {
	  count0_x0 = count[2]; // for convergence
	  x0 = xmid;
	  i0 = i0_;
	}
	else if (count[1]*mpi_split_size_half < count[0]*(mpi_split_size-mpi_split_size_half)) {
	  count0_x1 = count[2]; // for convergence
	  x1 = xmid;
	  i1 = i1_;
	}
	else {
	  // perfectly balanced - done!
	  break;
	}

      }

    }

    /*
      if (mpi_split_rank == 0) {
      cout << "SPLIT size: " << mpi_split_size << " count[0]: " << count[0] << " count[1]: " << count[1] << " sum: " << count[0]+count[1] << endl;
      }
    */

    /*
      cout << "rank: " << mpi_split_rank << " of " << mpi_split_size << " made it" << endl;
      if (mpi_split_rank == 0) getchar();
      MPI_Barrier(mpi_split_comm);
    */

    // we now know how to split the data...
    // we will redistribute with individual messages rather than all2allv. start by
    // exchanging everyone's counts and we all compute what messages will be going where...

    int8 (*my_count_of_rank)[2] = new int8[mpi_split_size][2];
    MPI_Allgather(my_count+2,2,MPI_INT8,(int8*)my_count_of_rank,2,MPI_INT8,mpi_split_comm);

    // the target number of cvs per core is...

    //cout << "XXXXX count[0]: " << count[0] << " count[1]: " << count[1] << endl;
    //cout << "XXXXX count[2]: " << count[2] << " count[3]: " << count[3] << endl;

    int8 my_count_target[2];
    // allowing for load imbalance, we need different targets for 0 and 1 processor groups...
    if (count[2]%mpi_split_size_half == 0)
      my_count_target[0] = count[2]/mpi_split_size_half;
    else
      my_count_target[0] = count[2]/mpi_split_size_half + 1;

    if (count[3]%(mpi_split_size-mpi_split_size_half) == 0)
      my_count_target[1] = count[3]/(mpi_split_size-mpi_split_size_half);
    else
      my_count_target[1] = count[3]/(mpi_split_size-mpi_split_size_half) + 1;

    double (*x_cv_new)[3] = NULL;
    int8 (*cv_global_wgt_new)[2] = NULL;
    int icv_new = 0;
    int icv_old = 0;
    if (mpi_split_rank < mpi_split_size_half) {
      // this rank is part of the first "half"...
      const int n_new = max(my_count[2],my_count_target[0]);
      x_cv_new          = new double[n_new][3];
      cv_global_wgt_new = new int8[n_new][2];
      // simultaneously copy values that are staying into _new arrays, and compress old...
      FOR_ICV {
	if (x_cv[icv][idir] + eps1*x_cv[icv][(idir+1)%3] + eps2*x_cv[icv][(idir+2)%3] <= xmid) {
	  x_cv_new[icv_new][0] = x_cv[icv][0];
	  x_cv_new[icv_new][1] = x_cv[icv][1];
	  x_cv_new[icv_new][2] = x_cv[icv][2];
	  cv_global_wgt_new[icv_new][0] = cv_global_wgt[icv][0];
	  cv_global_wgt_new[icv_new][1] = cv_global_wgt[icv][1];
	  ++icv_new;
	}
	else {
	  if (icv_old < icv) {
	    x_cv[icv_old][0] = x_cv[icv][0];
	    x_cv[icv_old][1] = x_cv[icv][1];
	    x_cv[icv_old][2] = x_cv[icv][2];
	    cv_global_wgt[icv_old][0] = cv_global_wgt[icv][0];
	    cv_global_wgt[icv_old][1] = cv_global_wgt[icv][1];
	  }
	  ++icv_old;
	}
      }
      assert(icv_new == my_count[2]);
      assert(icv_old == my_count[3]);
    }
    else {
      // this rank is part of the second "half"...
      const int n_new = max(my_count[3],my_count_target[1]);
      x_cv_new          = new double[n_new][3];
      cv_global_wgt_new = new int8[n_new][2];
      // same as above, but flipped...
      FOR_ICV {
	if (x_cv[icv][idir] + eps1*x_cv[icv][(idir+1)%3] + eps2*x_cv[icv][(idir+2)%3] <= xmid) {
	  if (icv_old < icv) {
	    x_cv[icv_old][0] = x_cv[icv][0];
	    x_cv[icv_old][1] = x_cv[icv][1];
	    x_cv[icv_old][2] = x_cv[icv][2];
	    cv_global_wgt[icv_old][0] = cv_global_wgt[icv][0];
	    cv_global_wgt[icv_old][1] = cv_global_wgt[icv][1];
	  }
	  ++icv_old;
	}
	else {
	  x_cv_new[icv_new][0] = x_cv[icv][0];
	  x_cv_new[icv_new][1] = x_cv[icv][1];
	  x_cv_new[icv_new][2] = x_cv[icv][2];
	  cv_global_wgt_new[icv_new][0] = cv_global_wgt[icv][0];
	  cv_global_wgt_new[icv_new][1] = cv_global_wgt[icv][1];
	  ++icv_new;
	}
      }
      assert(icv_new == my_count[3]);
      assert(icv_old == my_count[2]);
    }

    /*
      {
      int my_count = 0;
      for (int icv = 0; icv < icv_new; ++icv)
      my_count += cv_global_wgt_new[icv][1];
      int my_count2 = 0;
      for (int icv = 0; icv < icv_old; ++icv)
      my_count2 += cv_global_wgt[icv][1];
      cout << "CCCCCC mpi_split_rank: " << mpi_split_rank << " my_count: " << my_count << " " << my_count2 << endl;
      }
    */

    /*
      cout << "2 rank: " << mpi_split_rank << " of " << mpi_split_size << " made it" << endl;
      if (mpi_split_rank == 0) getchar();
      MPI_Barrier(mpi_split_comm);
    */

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
		  MPI_Issend(cv_global_wgt+icv_old,n_send_recv*2,MPI_INT8,rank1,
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
		  MPI_Irecv(cv_global_wgt_new+icv_new,n_send_recv*2,MPI_INT8,rank0,
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
		  MPI_Issend(cv_global_wgt+icv_old,n_send_recv*2,MPI_INT8,rank0,
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
		  MPI_Irecv(cv_global_wgt_new+icv_new,n_send_recv*2,MPI_INT8,rank1,
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

    if (cv_global_wgt != NULL) delete[] cv_global_wgt;
    cv_global_wgt = cv_global_wgt_new;

    ncv = icv_new;

    // check and report cv_global_wgt[1]...
    /*
      {
      int my_wgt = 0;
      for (int icv = 0; icv < ncv; ++icv)
      my_wgt += cv_global_wgt[icv][1];
      cout << "mpi_split_rank: " << mpi_split_rank << " key: " << (mpi_split_rank < mpi_split_size_half) << " wgt: " << my_wgt << endl;
      if (mpi_split_rank == 0) getchar();
      MPI_Barrier(mpi_split_comm);
      }
    */

    // wait for everyone to complete...

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
  }

  void reorderRecursive(vector< pair<double,int> >& pVec,const int i_begin,const int i_end,const int id,const double (*x_cv)[3]) {

    const double eps1 = 1.1234E-6; // include a little y and z to break sort ties (and hopefully not cause any)...
    const double eps2 = 2.7531E-6;

    assert(i_end-1 > i_begin);

    // pVec is currently sorted in id_sorted, so there is no need to figure out the range of this one

    // check...

    double xmin[3] = { 1.0E+20, 1.0E+20, 1.0E+20 };
    double xmax[3] = { -1.0E+20, -1.0E+20, -1.0E+20 };
    xmin[id] = pVec[i_begin].first;
    xmax[id] = pVec[i_end-1].first;
    for (int i = i_begin; i < i_end; ++i) {
      const int icv = pVec[i].second;
      FOR_I3 if (i != id) {
	const double x = x_cv[icv][i] + eps1*x_cv[icv][(i+1)%3] + eps2*x_cv[icv][(i+2)%3];
	xmin[i] = min(xmin[i],x);
	xmax[i] = max(xmax[i],x);
      }
    }

    int id_next;
    if (xmax[0]-xmin[0] >= max(xmax[1]-xmin[1],xmax[2]-xmin[2]))
      id_next = 0;
    else if (xmax[1]-xmin[1] >= max(xmax[0]-xmin[0],xmax[2]-xmin[2]))
      id_next = 1;
    else
      id_next = 2;

    if (id_next != id) {
      // id changed, so replace the "first" in pVec, sort, then pass...
      for (int i = i_begin; i < i_end; ++i) {
	const int icv = pVec[i].second;
	pVec[i].first = x_cv[icv][id_next] + eps1*x_cv[icv][(id_next+1)%3] + eps2*x_cv[icv][(id_next+2)%3];
      }
      vector< pair<double,int> >::iterator iter_f = pVec.begin(); iter_f += i_begin;
      vector< pair<double,int> >::iterator iter_l = pVec.begin(); iter_l += i_end;
      sort(iter_f,iter_l);
    }

    // no need to change or sort, just continue...
    const int i_mid = (i_begin+i_end)/2;
    if (i_mid-1 > i_begin) reorderRecursive(pVec,i_begin,i_mid,id_next,x_cv);
    if (i_end-1 > i_mid) reorderRecursive(pVec,i_mid,i_end,id_next,x_cv);

  }

}

// =======================================================================================
// =======================================================================================
// =======================================================================================
// =======================================================================================
// MpiStuff namespace routines...
// =======================================================================================
// =======================================================================================
// =======================================================================================
// =======================================================================================

namespace MpiStuff {

  // give single processor initial values to these
  // in practice, the user should call initMpiStuff
  // near the start of there code, after MPI_Init(*,*)...

  int mpi_rank = 0;
  int mpi_size = 1;
  MPI_Comm mpi_comm = MPI_COMM_WORLD;

  double mpi_wtime_0 = 0.0;

  int mpi_rank_shared = 0;
  int mpi_size_shared = 1;
  MPI_Comm mpi_comm_shared = MPI_COMM_NULL;

  int mpi_rank_internode = 0;
  int mpi_size_internode = 1;
  MPI_Comm mpi_comm_internode = MPI_COMM_NULL;

  int * rank_of_rank_internode = NULL;

  int mpi_rank_internode_actual = 0; 

  MPI_Offset int_size,
    float_size,
    int8_size,
    uint8_size,
    double_size,
    header_size;

  MPI_Datatype MPI_Header = MPI_DATATYPE_NULL;

  int ishm = 0; // index to help increment unique filenames 
  map<void*,SharedMemDescriptor> shm_map;


  //int phys_cores_per_node = 1;

  void initMpiStuff() {
    MPI_Comm_dup(MPI_COMM_WORLD, &mpi_comm);
    initMpiStuffCommon();
  }

  void initMpiStuff(MPI_Comm& comm) {
    MPI_Comm_dup(comm, &mpi_comm);
    initMpiStuffCommon();
  }

  void initMpiStuffCommon() {

    // set world rank and size...

    MPI_Comm_rank(mpi_comm, &mpi_rank);
    MPI_Comm_size(mpi_comm, &mpi_size);

    // record the start time...

    mpi_wtime_0 = MPI_Wtime();

    // check expected sizes...

    if ((sizeof(int) != 4) || (sizeof(int8) != 8) || (sizeof(MPI_Offset) != 8)|| (sizeof(float) != 4) || (sizeof(double) != 8)) {
      if (mpi_rank == 0) {
	cerr << "Error: unexpected sizes:" << endl;
	cerr << "int       : " << sizeof(int) << endl;
	cerr << "int8      : " << sizeof(int8) << endl;
	cerr << "MPI_Offset: " << sizeof(MPI_Offset) << endl;
	cerr << "float     : " << sizeof(float) << endl;
	cerr << "double    : " << sizeof(double) << endl;
      }
      throw(0);
    }

    // set type sizes...

    int_size    = 4ll;
    float_size  = 4ll;
    int8_size   = 8ll;
    uint8_size  = 8ll;
    double_size = 8ll;
    header_size = 256ll; // checked below

    MPI_Datatype oldtypes[6];
    int blockcounts[6];
    MPI_Aint offsets[6];

    // the name...
    offsets[0] = 0;
    oldtypes[0] = MPI_CHAR;
    blockcounts[0] = UGP_IO_HEADER_NAME_LEN;

    // the id...
    offsets[1] = offsets[0] + 1*UGP_IO_HEADER_NAME_LEN;
    oldtypes[1] = MPI_INT;
    blockcounts[1] = 1;

    // the offset skip...
    offsets[2] = offsets[1] + 4;
    oldtypes[2] = MPI_OFFSET_DATATYPE;
    blockcounts[2] = 1;

    // the 16 ints...
    offsets[3] = offsets[2] + 8;
    oldtypes[3] = MPI_INT;
    blockcounts[3] = 16;

    // the 12 (was 16) reals...
    offsets[4] = offsets[3] + 4*16;
    oldtypes[4] = MPI_DOUBLE;
    blockcounts[4] = 12;

    // the 4 uint8s...
    offsets[5] = offsets[4] + 8*12;
    oldtypes[5] = MPI_UINT8;
    blockcounts[5] = 4;

    // build the header type...
    MPI_Type_create_struct(6, blockcounts, offsets, oldtypes, &MPI_Header);
    MPI_Type_commit(&MPI_Header);

    int header_size_int;
    MPI_Aint header_lb,header_extent;
    MPI_Type_get_extent(MPI_Header, &header_lb, &header_extent);
    MPI_Type_size(MPI_Header, &header_size_int);

    // check that is is 256
    if ((header_extent != header_size_int) || (header_size_int != 256)) {
      if (mpi_rank == 0) {
	cerr << "Error: header size/extent problem: "<< endl;
	cerr << "header_extent : "<< header_extent << endl;
	cerr << "header_size  : "<< header_size_int << endl;
	int long_size;
	MPI_Type_size(MPI_LONG,&long_size);
	cerr << "long_size: " << long_size << endl;
      }
      throw(0);
    }

    // ===================================================================================================
    // currently the shared mem stuff is implemented using shm_open / mmap / etc. In the future, we could
    // use MPI3 routines to hide this stuff, but MPI3 is not uniformly supported for now...
    //
    // We recommend you add -DWITH_SHM to your flags, and you may have to add -lrt to CXXFLAGS
    // ===================================================================================================

#ifdef WITH_MPI3

    /*
      MPI_Comm_split_type(mpi_comm,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,&mpi_comm_shared);
      if (mpi_comm_shared == MPI_COMM_NULL)
      CERR("could not build mpi_comm_shared. Perhaps you are compiling with -D WITH_MPI3 but do not have MPI3");

      MPI_Comm_rank(mpi_comm_shared,&mpi_rank_shared);
      MPI_Comm_size(mpi_comm_shared,&mpi_size_shared);
    */

    CERR("WITH_MPI3 is not implemented properly. Use WITH_SHM for now.");

#endif

#ifdef WITH_SHM

    // first remove file in /dev/shm if present (possibly from a prior unclean exit)...
    shm_unlink(SHM_NAME);
    MPI_Barrier(mpi_comm);

    // everyone tries to create file in /dev/shm, but operation is atomic...
    int fd_shm = shm_open(SHM_NAME, O_CREAT | O_EXCL | O_RDWR, 0600);
    int one_core = 0;
    if (fd_shm != -1) {
      one_core = 1;
    }
    else {
      // just open file descriptor since already created...
      fd_shm = shm_open(SHM_NAME, O_RDWR, 0600);
      if (fd_shm == -1) {
	//cerr << " mpi_rank " << mpi_rank << " had shm_open error: " << strerror(errno) << endl;
	one_core = mpi_size+1;
      }
    }

    int nnodes;
    MPI_Allreduce(&one_core, &nnodes, 1, MPI_INT, MPI_SUM, mpi_comm);

    // error checking...

    if (nnodes > mpi_size) {

      if (mpi_rank == 0)
	cerr << "\n\n========================================\n" <<
	  "ERROR: shm_open error on these ranks:" << endl;

      int ierr = 0;
      if (one_core == mpi_size+1)
	ierr = 1;

      int * ierr_gather = NULL;
      if (mpi_rank == 0)
	ierr_gather = new int[mpi_size];
      MPI_Gather(&ierr,1,MPI_INT,ierr_gather,1,MPI_INT,0,mpi_comm);

      // reduce the offending names...

      char name[MPI_MAX_PROCESSOR_NAME];
      int len;
      MPI_Get_processor_name(name,&len);
      assert(len < MPI_MAX_PROCESSOR_NAME);
      assert(name[len] == '\0');

      char * all_names = NULL;
      if (mpi_rank == 0)
	all_names = new char[MPI_MAX_PROCESSOR_NAME*mpi_size];
      MPI_Gather(name,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,
		 all_names,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,0,mpi_comm);

      if (mpi_rank == 0) {

	std::set<std::string> badNameSet;
	FOR_RANK {
	  if (ierr_gather[rank] != 0) {
	    cerr << " > rank " << rank << ", ierr " << ierr_gather[rank] << endl;
	    badNameSet.insert(all_names+MPI_MAX_PROCESSOR_NAME*rank);
	  }
	}

	cerr << "Check these processors for problems: " << endl;
	for (std::set<std::string>::const_iterator iter = badNameSet.begin(); iter != badNameSet.end(); ++iter)
	  cerr << *iter << endl;

	delete[] all_names;
	delete[] ierr_gather;

	cerr << "========================================\n\n" << endl;

      }

      throw(0);

    }

    // now build intranode communicator...
    if (one_core) {
      // set file size..
      if (ftruncate(fd_shm, sizeof(int)) == -1) {
	cerr << " mpi_rank " << mpi_rank << " had ftrucate error: " << strerror(errno) << endl;
	throw(-1);
      }
    }

    MPI_Barrier(mpi_comm);

    int * rankptr;
    if (one_core) {
      // only one_core ranks are able to write...
      rankptr = (int *)mmap(NULL, sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, fd_shm, 0);
      *rankptr = mpi_rank;
    }
    else {
      rankptr = (int *)mmap(NULL, sizeof(int), PROT_READ, MAP_SHARED, fd_shm, 0);
    }

    MPI_Barrier(mpi_comm);

    MPI_Comm_split(mpi_comm, *rankptr, mpi_rank, &mpi_comm_shared);

    if (munmap(rankptr, sizeof(int)) == -1) {
      cerr << "rank " << mpi_rank << " had munmap() error: " << strerror(errno) << endl;
      throw(-1);
    }

    // don't need dev/shm file anymore, unlink removes it...
    if (one_core) {
      if (shm_unlink(SHM_NAME) == -1) {
	cerr << " mpi_rank " << mpi_rank << " had shm_unlink error: " << strerror(errno) << endl;
	throw(-1);
      }
    }

    // also don't need file descriptor..
    if (close(fd_shm) == -1) {
      cerr << " mpi_rank " << mpi_rank << " had close error: " << strerror(errno) << endl;
      throw(-1);
    }

    MPI_Comm_rank(mpi_comm_shared,&mpi_rank_shared);
    MPI_Comm_size(mpi_comm_shared,&mpi_size_shared);

    // now build internode communicator...
    // switch one_core to the first core...
    if (mpi_rank_shared == 0)
      one_core = 0;
    else
      one_core = 1;

    MPI_Comm_split(mpi_comm, one_core, mpi_rank, &mpi_comm_internode);

    if (one_core == 0) {
      MPI_Comm_rank(mpi_comm_internode, &mpi_rank_internode);
      MPI_Comm_size(mpi_comm_internode, &mpi_size_internode);
      assert(nnodes == mpi_size_internode);
    }
    else {
      // complement of communicator is not usefull..
      mpi_rank_internode = -1;
      //mpi_size_internode = -1; // now stored everywhere...
      MPI_Comm_free(&mpi_comm_internode); // XXXXXX should we set like this, or free?...
      mpi_comm_internode = MPI_COMM_NULL;
    }

    mpi_rank_internode_actual = mpi_rank_internode;
    MPI_Bcast(&mpi_rank_internode_actual,1,MPI_INT,0,mpi_comm_shared);

    // and share the internode size with all ranks: this is the number of nodes
    // and everyone needs to know this...
    MPI_Bcast(&mpi_size_internode,1,MPI_INT,0,mpi_comm);

    // at this point, all nodes in the shared communicator know the internode size, but
    // do NOT have a valid internode communicator, and their rank in the internode
    // comm is -1.

    // check that rank 0 is zero in all communicators...
    if (mpi_rank == 0) {
      assert(mpi_rank_shared == 0);
      assert(mpi_rank_internode == 0);
    }

    // use the inter-node communicator to gather the rank of the first
    // rank on each node...

    assert(rank_of_rank_internode == NULL);
    rank_of_rank_internode = new int[mpi_size_internode];

    if (mpi_rank_shared == 0) {
      MPI_Allgather(&mpi_rank,1,MPI_INT,rank_of_rank_internode,1,MPI_INT,mpi_comm_internode);
    }

    // now the list of ranks at the start of each node are in rank_of_rank_internode on
    // the first rank in each shared memory comm. Distribute to all members of each node...

    MPI_Bcast(rank_of_rank_internode,mpi_size_internode,MPI_INT,0,mpi_comm_shared);

#else // WITH_SHM

    // if the user has not specified to use shared memory, then setup the
    // communicators as if each core was its own independent node. In this
    // mode, CTI_Mmap and CTI_Munmap default to new and delete...

    MPI_Comm_split(mpi_comm, mpi_rank, mpi_rank, &mpi_comm_shared);
    MPI_Comm_rank(mpi_comm_shared,&mpi_rank_shared); assert(mpi_rank_shared == 0);
    MPI_Comm_size(mpi_comm_shared,&mpi_size_shared); assert(mpi_size_shared == 1);

    // and the internode communicator is the same as the standard communicator...

    MPI_Comm_dup(mpi_comm, &mpi_comm_internode);
    mpi_rank_internode = mpi_rank;
    mpi_size_internode = mpi_size;

    mpi_rank_internode_actual = mpi_rank;

#endif

  }

  void destroyMpiStuff() {

    MPI_TYPE_FREE(MPI_Header);
    MPI_Comm_free(&mpi_comm) ;

    if (mpi_comm_shared != MPI_COMM_NULL)
      MPI_Comm_free(&mpi_comm_shared) ;

    if (mpi_comm_internode != MPI_COMM_NULL)
      MPI_Comm_free(&mpi_comm_internode) ;

    if (rank_of_rank_internode) delete[] rank_of_rank_internode;

  }

  void MPI_Check(const int ierr,const char * message) {
    if (ierr != 0) {
      char str[MPI_MAX_ERROR_STRING];
      int length;
      MPI_Error_string(ierr,str,&length);
      cerr << "ERROR: " << message << ", " << str << endl;
      throw(-1);
    }
  }

  void MPI_Pause(const char * message) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank == 0) {
      cout << message << endl;
      getchar();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  void MPI_File_pause(const char * message) {
    if (mpi_rank == 0)
      cout << " > MPI_File_pause: " << message << ". Touch file \"go\" to continue..." << endl;
    int done = 0;
    while (done != 1) {
      if (mpi_rank == 0) {
	struct stat stFileInfo;
	//bool blnReturn;
	int intStat;
	// Attempt to get the file attributes
	intStat = stat("go",&stFileInfo);
	if(intStat == 0) {
	  // We were able to get the file attributes
	  // so the file obviously exists.
	  MPI_File_delete("go",MPI_INFO_NULL);
	  done = 1;
	}
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
      // do something costly...
      if (done != 1) {
	usleep(500000);
      }
    }
  }

  void MPI_Sync(const char * message,const bool reset) {
    // reset takes the default value of false in the header...
    // pass reset == true to reset the cumulative time.
    int my_data = 0;
    int data;
    MPI_Allreduce(&my_data,&data,1,MPI_INT,MPI_SUM,mpi_comm);
    if (data != 0) throw(0);
    if (mpi_rank == 0) {
      static double wtime0 = 0.0;
      static double wtime_prev = 0.0;
      double wtime = MPI_Wtime();
      if ( reset || (wtime0 == 0.0) ) {
	wtime0 = wtime_prev = wtime;
	cout << "MPI_Sync: " << message << string(max(3,int(30-strlen(message))),'.') << "reset." << endl;
      }
      else {
	cout << "MPI_Sync: " << message << string(max(3,int(30-strlen(message))),'.') << "time since last sync: " <<
	  setw(12) << wtime-wtime_prev << ", cumulative: " << setw(12) << wtime-wtime0 << endl;
	wtime_prev = wtime;
      }
    }
    my_data = 1;
    MPI_Allreduce(&my_data,&data,1,MPI_INT,MPI_SUM,mpi_comm);
    if (data != mpi_size) throw(0);
  }

  void MPI_Run_check() {
    
    if (mpi_rank == 0)
      cout << "======================= MPI_Run_check() ============================" << endl;

    char name[MPI_MAX_PROCESSOR_NAME];
    int namelen;
    MPI_Get_processor_name(name,&namelen);
    assert(namelen < MPI_MAX_PROCESSOR_NAME);
    name[namelen] = '\0';

    char * all_names = NULL;
    if (mpi_rank == 0) all_names = new char[MPI_MAX_PROCESSOR_NAME*mpi_size];
    MPI_Gather(name,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,all_names,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,0,mpi_comm);

    // now try and run something with some cost...

    double wtime,work;
    int8 loop_size = 1000000;
    int done = 0;
    while (done == 0) {

      loop_size *= 2;
      if (mpi_rank == 0) cout << " > trying work loop size: " << loop_size << endl;

      double wtime0 = MPI_Wtime();

      work = 0.0;
      for (int i = 0; i < loop_size; ++i) {
	work += double(rand());
      }
      work /= double(RAND_MAX)*loop_size;

      wtime = MPI_Wtime() - wtime0;

      double wtime_max;
      MPI_Reduce(&wtime,&wtime_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      if (mpi_rank == 0) {
	if (wtime_max > 1.0)
	  done = 1;
      }
      MPI_Bcast(&done,1,MPI_INT,0,mpi_comm);

    }

    double * dbuf = NULL;
    if (mpi_rank == 0) dbuf = new double[2*mpi_size];
    double my_dbuf[2] = { wtime, work };
    MPI_Gather(my_dbuf,2,MPI_DOUBLE,dbuf,2,MPI_DOUBLE,0,mpi_comm);

    if (mpi_rank == 0) {

      // do something with the work...
      double work_sum = 0.0;
      FOR_RANK {
	work_sum += dbuf[2*rank+1];
      }
      cout << " > work sum (ignore this, but must be reported so compiler does the work): " << work_sum/double(mpi_size) << endl;

      // sort the timing...
      vector<pair<double,int> > costRankPairVec(mpi_size);
      FOR_RANK {
	costRankPairVec[rank].first = dbuf[2*rank];
	costRankPairVec[rank].second = rank;
      }

      sort(costRankPairVec.begin(),costRankPairVec.end());

      const double fastest = costRankPairVec[0].first;
      cout << " > fastest time to do work: " << fastest << " [s]" << endl;

      bool b_slower = false;
      FOR_RANK {
	int sorted_rank = costRankPairVec[rank].second;
	if (costRankPairVec[rank].first >= 1.05*fastest) { // HACK 1.05
	  if (!b_slower)
	    cout << " > the following ranks are atleast 5% slower than the fastest\n" <<
	      " > rank, name, relative performance" << endl;
	  b_slower = true;
	  cout << " > " << sorted_rank << " " << all_names+MPI_MAX_PROCESSOR_NAME*sorted_rank << " " << dbuf[2*sorted_rank]/fastest << endl;
	}
      }

      if (!b_slower) cout << " > ALL ranks are within 5% of the fastest (this is good!)" << endl;

      delete[] dbuf;
      delete[] all_names;

      cout << "===================== MPI_Run_check() done =========================" << endl;

    }

  }

  int MPI_Type_indexed_clean(int n_zones, int *my_zone_count, int *my_zone_disp,
			     MPI_Datatype old_type, MPI_Datatype *new_type) {

    // on some systems - e.g. Argonne Intrepid BGL, JWN's mississippi account?
    // the zero counts present in this indexed type can cause problems. This
    // "clean" type removes these zeros, but leaves the non-zero counts and displacements
    // unchanged.

    int new_n_zones = 0;
    int * new_zone_count = new int[n_zones];
    int * new_zone_disp = new int[n_zones];

    // remove zero-length zones (blocks)
    for (int i=0; i<n_zones; ++i) {
      if (my_zone_count[i] > 0) {
	new_zone_count[new_n_zones] = my_zone_count[i];
	new_zone_disp[new_n_zones] = my_zone_disp[i];
	new_n_zones++;
      }
    }

    // do not allow a zero in new_n_zones. is this a problem?
    if (new_n_zones == 0) {
      new_n_zones = 1;
      new_zone_count[0] = my_zone_count[0];
      new_zone_disp[0] = my_zone_disp[0];
    }

    // now look for huge displacements...
    int my_flag_large = 0;
    for (int i = 0; i < new_n_zones; ++i)
      if (new_zone_disp[i] < 0)
	my_flag_large = 1;

    // make sure we all build the same type...
    int flag_large;
    MPI_Allreduce(&my_flag_large,&flag_large,1,MPI_INT,MPI_MAX,mpi_comm);
    if (flag_large == 1) {
      if (mpi_rank == 0)
	cerr << "Error: negatives detected in MPI_Type_indexed_clean. Global count must exceed 2^31" << endl;
      throw(0);
      //new_zone_disp2[i] += 4294967296ll;
    }

    // call the original routine
    int ierr = MPI_Type_indexed(new_n_zones,new_zone_count,new_zone_disp,
				old_type,new_type);

    delete[] new_zone_count;
    delete[] new_zone_disp;

    return(ierr);

  }

  int MPI_Type_indexed_clean(int n_zones, int *my_zone_count, int8 *my_zone_disp,
			     MPI_Datatype old_type, MPI_Datatype *new_type) {

    // on some systems - e.g. Argonne Intrepid BGL, JWN's mississippi account?
    // the zero counts present in this indexed type can cause problems. This
    // "clean" type removes these zeros, but leaves the non-zero counts and displacements
    // unchanged.

    int8 typesize;
    if (old_type == MPI_INT)
      typesize = 4ll;
    else if (old_type == MPI_DOUBLE)
      typesize = 8ll;
    else
      CERR("unsupported type in MPI_Type_indexed_clean: " << old_type);

    // this particular routine also builds an hindexed type because it is the
    // only one available with 64-bit int addressing...

    int new_n_zones = 0;
    int * new_zone_count = new int[n_zones];
    MPI_Aint * new_zone_disp_bytes = new MPI_Aint[n_zones];
    assert( sizeof(MPI_Aint) == 8 ); // assume the arcitectural int is 8 bytes
    assert( sizeof(int) == 4 ); // int better be 4

    // remove zero-length zones (blocks)
    for (int i=0; i<n_zones; ++i) {
      if (my_zone_count[i] > 0) {
	new_zone_count[new_n_zones] = my_zone_count[i];
	new_zone_disp_bytes[new_n_zones] = MPI_Aint(my_zone_disp[i])*typesize; // byte displacements
	new_n_zones++;
      }
    }

    // do not allow a zero in new_n_zones. is this a problem?
    if (new_n_zones == 0) {
      new_n_zones = 1;
      new_zone_count[0] = my_zone_count[0];
      new_zone_disp_bytes[0] = my_zone_disp[0]*typesize; // byte displacements
    }

    /*
      cout << "about to create hindexed type: " << new_n_zones << endl;
      for (int i = 0; i < new_n_zones; ++i)
      cout << "i=" << i << ", new_zone_count=" << new_zone_count[i] << ", new_zone_disp_bytes=" << new_zone_disp_bytes[i] << endl;
    */

    // call the original routine
    int ierr = MPI_Type_create_hindexed(new_n_zones,new_zone_count,new_zone_disp_bytes,old_type,new_type);

    delete[] new_zone_count;
    delete[] new_zone_disp_bytes;

    return(ierr);

  }

  int MPI_Bcast_string(std::string& str,int rank,MPI_Comm& comm) {
    (void)(comm);
    int nchar;

    if (rank == mpi_rank)
      nchar = int(str.length())+1; // add one for trailing \0

    // send size and max length to everyone...

    MPI_Bcast(&nchar, 1, MPI_INT, rank, mpi_comm);

    // exchange char buffer...

    char * cbuf = new char[nchar];

    if (rank == mpi_rank) {
      for (int i = 0; i < int(str.length()); ++i) {
	cbuf[i] = str[i];
      }
      cbuf[nchar-1] = '\0';
    }

    MPI_Bcast(cbuf, nchar, MPI_CHAR, rank, mpi_comm);

    if (rank != mpi_rank) {
      str = cbuf;
    }

    delete[] cbuf;

    return(0);

  }

  int MPI_Bcast_stringVec(std::vector<std::string>& stringVec,int rank,MPI_Comm& comm) {
    (void)(comm);
    // this routine copies the stringVec on rank to all other ranks...

    int ibuf[2]; // stores the size of the stringVec, and the max string length + 1

    if (rank == mpi_rank) {
      ibuf[0] = stringVec.size();
      ibuf[1] = 0;
      for (unsigned int i = 0; i < stringVec.size(); ++i)
	ibuf[1] = max(ibuf[1],int(stringVec[i].length())+1); // add one for trailing \0
    }

    // send size and max length to everyone...

    MPI_Bcast(ibuf, 2, MPI_INT, rank, mpi_comm);

    // exchange char buffer...

    const int nchar = ibuf[0]*ibuf[1];
    char * cbuf = new char[nchar];

    if (rank == mpi_rank) {
      for (unsigned int i = 0; i < stringVec.size(); ++i) {
	for (int j = 0; j < int(stringVec[i].length()); ++j) {
	  cbuf[ibuf[1]*i+j] = stringVec[i][j];
	}
	for (int j = int(stringVec[i].length()); j < ibuf[1]; ++j) {
	  cbuf[ibuf[1]*i+j] = '\0';
	}
      }
    }

    MPI_Bcast(cbuf, nchar, MPI_CHAR, rank, mpi_comm);

    if (rank != mpi_rank) {
      //assert(stringVec.size() == 0);
      stringVec.resize(ibuf[0]);
      for (unsigned int i = 0; i < stringVec.size(); ++i) {
	stringVec[i] = cbuf + ibuf[1]*i;
	//cout << "stringVec[i] = \"" << stringVec[i] << "\"" << endl;
      }
    }

    delete[] cbuf;
    return 0;
  }

  void checkMemoryUsage() {

    // report minimum memory available.  the type of
    // hardware must be defined at compile time bc this
    // check is architecture dependent...

#if defined(BGP_MEMCHECK) || defined(BGQ_MEMCHECK)

    bg_long free_mem, min_free_mem ;
    BlueGeneMemCheck::getFreeMem(&free_mem) ;

    MPI_Reduce(&free_mem, &min_free_mem, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, mpi_comm);
    if ( mpi_rank == 0 )
      cout << " >> Min mem available = " << min_free_mem <<endl ;

#elif defined(LINUX_MEMCHECK)

    long free_mem, min_free_mem ;
    LinuxMemcheck::getFreeMem(&free_mem);


    MPI_Reduce(&free_mem,&min_free_mem, 1, MPI_LONG, MPI_MIN, 0, mpi_comm) ;
    if ( mpi_rank == 0 )
      cout << " >> Min mem available (MB) = " << min_free_mem <<endl  ;

#else
    if ( mpi_rank == 0 )
      cout << " Mem check not implemented yet ... " << endl ;
#endif
  }

  template<>
  MPI_Datatype getMpiDatatype<char>() {
    return MPI_CHAR ;
  }

  template<>
  MPI_Datatype getMpiDatatype<int>() {
    return MPI_INT ;
  }

  template<>
  MPI_Datatype getMpiDatatype<int8>() {
    return MPI_INT8 ;
  }

  template<>
  MPI_Datatype getMpiDatatype<uint8>() {
    return MPI_UINT8 ;
  }

  template<>
  MPI_Datatype getMpiDatatype<float>() {
    return MPI_FLOAT ;
  }

  template<>
  MPI_Datatype getMpiDatatype<double>() {
    return MPI_DOUBLE ;
  }

  template<>
  MPI_Datatype getMpiDatatype<unsigned short>() {
    return MPI_UNSIGNED_SHORT ;
  }

  size_t sizeOfMpiDatatype(const MPI_Datatype datatype) {

    if (datatype == MPI_CHAR)
      return sizeof(char);
    else if (datatype == MPI_UNSIGNED_CHAR)
      return sizeof(unsigned char);
    else if (datatype == MPI_UNSIGNED_SHORT)
      return sizeof(unsigned short);
    else if (datatype == MPI_INT) 
      return sizeof(int);
    else if (datatype == MPI_INT8) 
      return sizeof(int8);
    else if (datatype == MPI_UINT8) 
      return sizeof(uint8);
    else if (datatype == MPI_FLOAT)
      return sizeof(float);
    else if (datatype == MPI_DOUBLE) 
      return sizeof(double);
    else if (datatype == MPI_Header) 
      return sizeof(Header);
    else 
      assert(0);
    return 0;

  }

#define T unsigned char
# include "Mmap_Munmap.hpp"
#undef T

#define T int
# include "Mmap_Munmap.hpp"
#undef T

#define T uint
# include "Mmap_Munmap.hpp"
#undef T

#define T float
# include "Mmap_Munmap.hpp"
#undef T

#define T double
# include "Mmap_Munmap.hpp"
#undef T

#define T int8
# include "Mmap_Munmap.hpp"
#undef T

#define T uint8
# include "Mmap_Munmap.hpp"
#undef T

#define T uint2
# include "Mmap_Munmap.hpp"
#undef T

  // ==============================================================
  // padt partitioning routines...
  // ==============================================================

  void setPartKind(const string& name) {

    if (name == "PADT") {
      MpiStuffHelper::part_kind = 0;
    }
    else if (name == "BAD") {
      MpiStuffHelper::part_kind = 1;
    }
    else {
      if (mpi_rank == 0) cout << "Warning: unrecognized SET_PART_KIND \"" << name << "\". Skipping." << endl;
    }

  }

  void repartXcvPadt(double (*&x_cv)[3],int8 *&icv_global,int &ncv,const MPI_Comm& mpi_comm,const int bits) {

    // ========================================================================================
    // warning: this changes, potentially reallocates and reorders x_cv,icv_global,and ncv...
    // ========================================================================================

    int mpi_split_size;
    MPI_Comm_size(mpi_comm, &mpi_split_size);
    MPI_Comm mpi_split_comm; MPI_Comm_dup(mpi_comm,&mpi_split_comm);

    while (mpi_split_size > 1) {

      // this routine returns x_cv, icv_global, ncv, and mpi_split_comm modified to continue...
      MpiStuffHelper::splitAndRedist(x_cv,icv_global,ncv,mpi_split_comm,bits);

      // and get the size again...
      MPI_Comm_size(mpi_split_comm,&mpi_split_size);

    }

    MPI_Comm_free(&mpi_split_comm);

  }

  void repartXcvRandom(double (*&x_cv)[3],int8 (*&cv_global_wgt)[2],int &ncv,const MPI_Comm& mpi_comm) {

    int mpi_size;
    MPI_Comm_size(mpi_comm,&mpi_size);

    int * send_count = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;
    
    int * raocv = new int[ncv];
    FOR_ICV {
      raocv[icv] = rand()%mpi_size;
      assert((raocv[icv] >= 0)&&(raocv[icv] < mpi_size));
      ++send_count[raocv[icv]];
    }

    int * send_disp = new int[mpi_size];
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
    assert(send_count_sum == ncv);
    
    int8 * send_buf_int8 = new int8[send_count_sum*2];
    double * send_buf_double = new double[send_count_sum*3];

    FOR_ICV {
      const int rank = raocv[icv];
      FOR_I2 send_buf_int8[send_disp[rank]*2+i] = cv_global_wgt[icv][i];
      FOR_I3 send_buf_double[send_disp[rank]*3+i] = x_cv[icv][i];
      ++send_disp[rank];
    }
    delete[] raocv;
    delete[] cv_global_wgt;
    delete[] x_cv;

    // rewind...
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    
    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    ncv = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    FOR_RANK {
      send_count[rank] *= 2;
      send_disp[rank] *= 2;
      recv_count[rank] *= 2;
      recv_disp[rank] *= 2;
    }
    
    cv_global_wgt = new int8[ncv][2];
    MPI_Alltoallv(send_buf_int8,send_count,send_disp,MPI_INT8,
                  (int8*)cv_global_wgt,recv_count,recv_disp,MPI_INT8,
                  mpi_comm);
    
    delete[] send_buf_int8;

    FOR_RANK {
      send_count[rank] = (send_count[rank]/2)*3;
      send_disp[rank] = (send_disp[rank]/2)*3;
      recv_count[rank] = (recv_count[rank]/2)*3;
      recv_disp[rank] = (recv_disp[rank]/2)*3;
    }

    x_cv = new double[ncv][3];
    MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
                  (double*)x_cv,recv_count,recv_disp,MPI_DOUBLE,
                  mpi_comm);
    
    delete[] send_buf_double;
    delete[] send_count;
    delete[] send_disp;
    delete[] recv_count;
    delete[] recv_disp;
    
  }

  void repartXcvWeightedPadt(double (*&x_cv)[3],int8 (*&cv_global_wgt)[2],int &ncv,const MPI_Comm& mpi_comm,const int bits) {

    // ========================================================================================
    // warning: this changes, potentially reallocates and reorders x_cv,icv_global,and ncv...
    // ========================================================================================

    if (bits == 13) {
      // use "lucky number" 13 to indicate a "bad" partition...
      repartXcvRandom(x_cv,cv_global_wgt,ncv,mpi_comm);
      return;
    }
    
    int mpi_split_size;
    MPI_Comm_size(mpi_comm, &mpi_split_size);
    MPI_Comm mpi_split_comm; MPI_Comm_dup(mpi_comm,&mpi_split_comm);

    while (mpi_split_size > 1) {

      // this routine returns x_cv, icv_global, ncv, and mpi_split_comm modified to continue...
      MpiStuffHelper::splitAndRedistWeighted(x_cv,cv_global_wgt,ncv,mpi_split_comm,bits);

      // and get the size again...
      MPI_Comm_size(mpi_split_comm,&mpi_split_size);

    }

    MPI_Comm_free(&mpi_split_comm);

  }

  void calcCvPartPadt(int * cv_part,const double (*x_cv_orig)[3],const int ncv,const MPI_Comm& mpi_comm,const int bits) {

    // -------------------------------------------------------------------------------
    // this routine uses parallel recursive bisection to split the passed points into
    // logically blocked partitions with near-perfect load balance...
    // -------------------------------------------------------------------------------

    // here we allow for a different passed comm -- for example the nodal comm...

    int mpi_rank,mpi_size;
    MPI_Comm_rank(mpi_comm, &mpi_rank);
    MPI_Comm_size(mpi_comm, &mpi_size);

    if (mpi_rank == 0)
      cout << "calcCvPartPadt(), size: " << mpi_size << endl;


    //dumpRange(x_cv_orig,ncv,"X_CV_ORIG");


    if (mpi_size == 0) {
      FOR_ICV cv_part[icv] = 0;
      return;
    }

    // set cv_part to -1 everywhere for checking...

    FOR_ICV cv_part[icv] = -1;

    // and we need cvora...

    int8 * cvora = new int8[mpi_size+1];
    int8 ncv_local = (int8)ncv;
    MPI_Allgather(&ncv_local,1,MPI_INT8,cvora+1,1,MPI_INT8,mpi_comm);
    cvora[0] = 0;
    FOR_RANK cvora[rank+1] += cvora[rank];
    const int8 ncv_global = cvora[mpi_size];

    if (mpi_rank == 0) {
      if (ncv_global < 1000)
	cout << " > ncv_global: " << ncv_global << endl;
      else if (ncv_global < 1000000)
	cout << " > ncv_global: " << ncv_global << " (" << double(ncv_global/100)*0.1 << "K)" << endl;
      else if (ncv_global < 1000000000)
	cout << " > ncv_global: " << ncv_global << " (" << double(ncv_global/100000)*0.1 << "M)" << endl;
      else
	cout << " > ncv_global: " << ncv_global << " (" << double(ncv_global/100000000)*0.1 << "B)" << endl;
    }

    assert( ncv == cvora[mpi_rank+1]-cvora[mpi_rank] );

    // copy the coordinate data, because it gets rearranged during the partitioning...

    double (*x_cv)[3] = new double[ncv][3];
    int8 * icv_global = new int8[ncv]; // gets changed too, but this is local so it does not matter
    FOR_ICV {
      x_cv[icv][0]    = x_cv_orig[icv][0];
      x_cv[icv][1]    = x_cv_orig[icv][1];
      x_cv[icv][2]    = x_cv_orig[icv][2];
      icv_global[icv] = cvora[mpi_rank] + icv;
    }
    int ncv_new = ncv;

    // recursive subdivision and reditribution of x_cv and icv_global. Note:
    // this (potentially) changes ncv_new...

    double wtime0 = MPI_Wtime();
    const double wtime00 = wtime0;
    vector<double> wtimeVec;
    vector<int> ncvVec;

    {

      int mpi_split_size = mpi_size;
      MPI_Comm mpi_split_comm; MPI_Comm_dup(mpi_comm,&mpi_split_comm);

      while (mpi_split_size > 1) {

	// this routine returns x_cv, ncv_new, and mpi_split_comm modified to continue...
	MpiStuffHelper::splitAndRedist(x_cv,icv_global,ncv_new,mpi_split_comm,bits);

	// and get the size again...
	MPI_Comm_size(mpi_split_comm, &mpi_split_size);

	const double new_wtime = MPI_Wtime();
	wtimeVec.push_back(new_wtime-wtime0);
	wtime0 = new_wtime;

	ncvVec.push_back(ncv_new);

      }

      MPI_Comm_free(&mpi_split_comm);

    }

    delete[] x_cv;

    {

      // now the x_cv and icv_global are partitioned. We need to return the
      // partition index to cv_part in the original ncv-striped data...

      int * send_count = new int[mpi_size];
      FOR_RANK send_count[rank] = 0;

      int * send_disp = new int[mpi_size];
      int * send_buf = NULL;

      for (int iter = 0; iter < 2; ++iter) {
	for (int icv_new = 0; icv_new < ncv_new; ++icv_new) {
	  const int rank = MpiStuffHelper::getRankInXoraNew(icv_global[icv_new],cvora,mpi_size);
	  const int icv = int(icv_global[icv_new] - cvora[rank]);
	  assert((icv >= 0)&&(icv < int(cvora[rank+1]-cvora[rank])));
	  if (iter == 0) {
	    // only pack non-local...
	    if (rank != mpi_rank) ++send_count[rank];
	  }
	  else {
	    if (rank != mpi_rank) {
	      send_buf[send_disp[rank]++] = icv;
	    }
	    else {
	      // for local icv's, set cv_part directly...
	      assert(cv_part[icv] == -1);
	      cv_part[icv] = mpi_rank;
	    }
	  }
	}
	// always compute send_disp...
	send_disp[0] = 0;
	for (int rank = 1; rank < mpi_size; ++rank)
	  send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
	// on the first iter, allocate the buffer...
	if (iter == 0) {
	  const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
	  send_buf = new int[send_count_sum];
	}
      }

      int * recv_count = new int[mpi_size];
      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      int * recv_disp = new int[mpi_size];
      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
	recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

      const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
      int * recv_buf = new int[recv_count_sum];

      MPI_Alltoallv(send_buf,send_count,send_disp,MPI_INT,
		    recv_buf,recv_count,recv_disp,MPI_INT,
		    mpi_comm);
      delete[] send_buf;
      delete[] send_count;
      delete[] send_disp;

      FOR_RANK {
	for (int irecv = recv_disp[rank]; irecv != recv_disp[rank]+recv_count[rank]; ++irecv) {
	  const int icv = recv_buf[irecv];
	  assert((icv >= 0)&&(icv < ncv));
	  assert(cv_part[icv] == -1);
	  cv_part[icv] = rank;
	}
      }
      delete[] recv_buf;
      delete[] recv_count;
      delete[] recv_disp;

      // check...
      FOR_ICV {
	assert((cv_part[icv] >= 0)&&(cv_part[icv] < mpi_size));
      }

    }

    delete[] icv_global;
    delete[] cvora;

    // report timing...

    int my_n_levels = wtimeVec.size();
    int n_levels;
    MPI_Allreduce(&my_n_levels,&n_levels,1,MPI_INT,MPI_MAX,mpi_comm);

    // make sure everyone has the same number of levels. This will not
    // be exactly true in general...

    for (int i = my_n_levels; i < n_levels; ++i) {
      wtimeVec.push_back(0.0);
      ncvVec.push_back(0);
    }

    double * wtime_max = NULL;
    if (mpi_rank == 0) wtime_max = new double[n_levels];
    MPI_Reduce(&wtimeVec.front(),wtime_max,n_levels,MPI_DOUBLE,MPI_MAX,0,mpi_comm);

    double * wtime_min = NULL;
    if (mpi_rank == 0) wtime_min = new double[n_levels];
    MPI_Reduce(&wtimeVec.front(),wtime_min,n_levels,MPI_DOUBLE,MPI_MIN,0,mpi_comm);

    int * ncv_max = NULL;
    if (mpi_rank == 0) ncv_max = new int[n_levels];
    MPI_Reduce(&ncvVec.front(),ncv_max,n_levels,MPI_INT,MPI_MAX,0,mpi_comm);

    if (mpi_rank == 0) {
      cout << " > partition time: " << MPI_Wtime()-wtime00 << " [s]. Partition timing per level: " << endl;
      for (int i = 0; i < n_levels; ++i)
	cout << " > level: " << i << " max wtime: " << wtime_max[i] << " [s], min wtime: " << wtime_min[i] <<
	  " imbalance: " << double(ncv_max[i])/(double(ncv_global)/double(mpi_size))-1.0 << endl;
      delete[] wtime_max;
      delete[] wtime_min;
      delete[] ncv_max;
    }

    //dumpRange(x_cv,ncv,"X_CV FINAL");

    /*
      FOR_RANK {
      if (mpi_rank == rank) {
      FOR_ICV {
      cout << "XXXXX " << x_cv[icv][0] << " " << x_cv[icv][1] << " " << x_cv[icv][2] << " " << mpi_rank << endl;
      }
      }
      MPI_Pause("X");
      }
    */

    // check the partition quality...

    {
      int ncv_new_max;
      MPI_Reduce(&ncv_new,&ncv_new_max,1,MPI_INT,MPI_MAX,0,mpi_comm);
      if (mpi_rank == 0)
	cout << " > final partition imbalance: " << double(ncv_new_max)/(double(ncv_global)/double(mpi_size))-1.0 << endl;
    }

  }

  void reorderXcv(double (*x_cv)[3],int8 *icv_global,const int ncv) {

    // go through the data and figure out which direction has the largest variation...

    const double eps1 = 1.1234E-6; // include a little y and z to break sort ties (and hopefully not cause any)...
    const double eps2 = 2.7531E-6; // include a little y and z to break sort ties (and hopefully not cause any)...

    double xmin[3] = { 1.0E+20, 1.0E+20, 1.0E+20 };
    double xmax[3] = { -1.0E+20, -1.0E+20, -1.0E+20 };
    FOR_ICV {
      FOR_I3 {
	const double x = x_cv[icv][i] + eps1*x_cv[icv][(i+1)%3] + eps2*x_cv[icv][(i+2)%3];
	xmin[i] = min(xmin[i],x);
	xmax[i] = max(xmax[i],x);
      }
    }

    int id;
    if (xmax[0]-xmin[0] >= max(xmax[1]-xmin[1],xmax[2]-xmin[2]))
      id = 0;
    else if (xmax[1]-xmin[1] >= max(xmax[0]-xmin[0],xmax[2]-xmin[2]))
      id = 1;
    else
      id = 2;

    vector< pair<double,int> > pVec(ncv);
    FOR_ICV {
      pVec[icv].first = x_cv[icv][id] + eps1*x_cv[icv][(id+1)%3] + eps2*x_cv[icv][(id+2)%3];
      pVec[icv].second = icv;
    }
    sort(pVec.begin(),pVec.end());

    if (ncv/2-1 > 0) MpiStuffHelper::reorderRecursive(pVec,0,ncv/2,id,x_cv);
    if (ncv-1 > ncv/2) MpiStuffHelper::reorderRecursive(pVec,ncv/2,ncv,id,x_cv);

    double (*x_cv_tmp)[3] = new double[ncv][3];
    int8 *icv_global_tmp = new int8[ncv];
    FOR_ICV {
      FOR_I3 x_cv_tmp[icv][i] = x_cv[icv][i];
      icv_global_tmp[icv] = icv_global[icv];
    }
    FOR_ICV {
      const int icv_new = pVec[icv].second;
      FOR_I3 x_cv[icv][i] = x_cv_tmp[icv_new][i];
      icv_global[icv]     = icv_global_tmp[icv_new];
    }
    delete[] x_cv_tmp;
    delete[] icv_global_tmp;

    // reordering done, now print...

    /*
      char filename[128];
      sprintf(filename,"xcv.%08d.dat",mpi_rank);
      FILE * fp = fopen(filename,"w");
      FOR_ICV {
      fprintf(fp,"%lf %lf %lf %d\n",x_cv[icv][0],x_cv[icv][1],x_cv[icv][2],int(icv/(ncv/12)));
      }
      fclose(fp);
    */

  }

#ifdef SERIAL_IO

// in the mpich implementation, MPI_File is an int...  

int ifile_desc = 0;
const int max_file_desc = 32;
int file_desc[max_file_desc];
map<MPI_File,pair<FILE*,MPI_Comm> > file_map;

// functions must be outside the namespace to overwrite MPI...

#endif

}

#ifdef SERIAL_IO

int MPI_File_open(MPI_Comm comm, const char *filename, int amode,
    MPI_Info info, MPI_File * fh) { 
  using namespace MpiStuff;

  //if ( mpi_rank == 0) 
  //  cout << " intercepted an mpi file open call " << endl;

  // create the file handle for everyone ; this is a 
  // collective operation... 

  
  (*fh) = (MPI_File)(&file_desc[ifile_desc]);
  ++ifile_desc;
  assert(ifile_desc <= max_file_desc);

  FILE * fp = NULL;

  if ( mpi_rank == 0) { 

    if ( (amode & MPI_MODE_WRONLY)) { 

      if (amode & MPI_MODE_CREATE)
        fp = fopen(filename, "wb");
      else
        fp = fopen(filename, "ab"); // maybe we should just use this

    } else if ( (amode & MPI_MODE_RDONLY)) { 

      fp = fopen(filename, "rb");

    } else { 

      assert(0);

    }

  } 

  file_map[*fh] = pair<FILE*,MPI_Comm>(fp,comm);

  return MPI_SUCCESS;

} 

int MPI_File_read(MPI_File fh, void *buf, int count,
    MPI_Datatype datatype, MPI_Status *status) { 
  using namespace MpiStuff;
  
  //if ( mpi_rank == 0) 
  //  cout << " intercepted a read call... " << endl;

  FILE * fp; 
  map<MPI_File,pair<FILE*,MPI_Comm> >::iterator it = file_map.find(fh);
  assert( it != file_map.end());

  fp = it->second.first;
  assert( fp != NULL); 

  size_t nbytes = count*sizeOfMpiDatatype(datatype);

  fread(buf,sizeof(char),nbytes,fp);

  return MPI_SUCCESS;

} 

int MPI_File_read_at(MPI_File fh, MPI_Offset offset, void *buf,
    int count, MPI_Datatype datatype, MPI_Status * status) { 
  using namespace MpiStuff;

  //if ( mpi_rank == 0) 
  //  cout << " intercepted a read at call... " << endl;

  FILE * fp; 
  map<MPI_File,pair<FILE*,MPI_Comm> >::iterator it = file_map.find(fh);
  assert( it != file_map.end());

  fp = it->second.first;
  assert( fp != NULL); 

  size_t nbytes = count*sizeOfMpiDatatype(datatype);

  fseek(fp,offset,SEEK_SET);
  fread(buf,sizeof(char),nbytes,fp);

  return MPI_SUCCESS;

} 

int MPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void *buf,
    int count, MPI_Datatype datatype, MPI_Status * status) { 
  using namespace MpiStuff;

  //if ( mpi_rank == 0) 
  //  cout << " intercepted a read at all call... " << endl;

  FILE * fp; 
  map<MPI_File,pair<FILE*,MPI_Comm> >::iterator it = file_map.find(fh);
  assert( it != file_map.end());

  fp = it->second.first;

  MPI_Comm mpi_file_comm = it->second.second;
  int mpi_file_rank,mpi_file_size;
  MPI_Comm_rank(mpi_file_comm, &mpi_file_rank);
  MPI_Comm_size(mpi_file_comm, &mpi_file_size);

  int *counts = NULL;
  if (mpi_file_rank == 0) counts = new int[mpi_file_size];
  MPI_Gather(&count,1,MPI_INT,counts,1,MPI_INT,0,mpi_file_comm);
  MPI_Offset *offsets = NULL;
  if (mpi_file_rank == 0) offsets = new MPI_Offset[mpi_file_size];
  MPI_Gather(&offset,1,MPI_OFFSET,offsets,1,MPI_OFFSET,0,mpi_file_comm);

  void* tmp_buf = NULL;
  if (mpi_file_rank == 0) {
    assert( fp != NULL); 
    int max_count = 0;
    for (int rank = 0; rank < mpi_file_size; ++rank) 
      max_count = max(max_count,counts[rank]);
    size_t max_nbytes = max_count*sizeOfMpiDatatype(datatype);
    tmp_buf = (void*)(new char[max_nbytes]);
    for (int rank = 0; rank < mpi_file_size; ++rank) {
      size_t nbytes = counts[rank]*sizeOfMpiDatatype(datatype);
      fseek(fp,offsets[rank],SEEK_SET);
      if (rank == 0) {
        fread(buf,sizeof(char),nbytes,fp);
      }
      else {
        fread(tmp_buf,sizeof(char),nbytes,fp);
        MPI_Send(tmp_buf,counts[rank],datatype,rank,rank,mpi_file_comm);
      }
    }
    delete[] counts;
    delete[] offsets;
  }
  else {
    MPI_Recv(buf,count,datatype,0,mpi_file_rank,mpi_file_comm,MPI_STATUS_IGNORE);
  }

  if (tmp_buf) delete[] (char*)(tmp_buf);

  return MPI_SUCCESS;

} 

int MPI_File_write(MPI_File fh, const void *buf, int count,
    MPI_Datatype datatype, MPI_Status * status) { 
  using namespace MpiStuff;

  //if ( mpi_rank == 0) 
  //  cout << " intercepted a write call... " << endl;

  FILE * fp; 
  map<MPI_File,pair<FILE*,MPI_Comm> >::iterator it = file_map.find(fh);
  assert( it != file_map.end());

  fp = it->second.first;
  assert( fp != NULL); 

  size_t nbytes = count*sizeOfMpiDatatype(datatype);

  fwrite((void*)buf,sizeof(char),nbytes,fp);

  return MPI_SUCCESS;

} 

int MPI_File_write_at(MPI_File fh, MPI_Offset offset, const void *buf,
    int count, MPI_Datatype datatype, MPI_Status * status) { 
  using namespace MpiStuff;

  //if ( mpi_rank == 0) 
  //  cout << " intercepted a write at call... " << endl;

  FILE * fp; 
  map<MPI_File,pair<FILE*,MPI_Comm> >::iterator it = file_map.find(fh);
  assert( it != file_map.end());

  fp = it->second.first;
  assert(fp != NULL);

  size_t nbytes = count*sizeOfMpiDatatype(datatype);

  fseek(fp,offset,SEEK_SET);
  fwrite((void*)buf,sizeof(char),nbytes,fp);

  return MPI_SUCCESS;

} 

int MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, const void *buf,
    int count, MPI_Datatype datatype, MPI_Status * status) { 
  using namespace MpiStuff;

  //if ( mpi_rank == 0) 
  //  cout << " intercepted a write at all call... " << endl;

  FILE * fp; 
  map<MPI_File,pair<FILE*,MPI_Comm> >::iterator it = file_map.find(fh);
  assert( it != file_map.end());

  fp = it->second.first;

  MPI_Comm mpi_file_comm = it->second.second;
  int mpi_file_rank,mpi_file_size;
  MPI_Comm_rank(mpi_file_comm, &mpi_file_rank);
  MPI_Comm_size(mpi_file_comm, &mpi_file_size);

  int *counts = NULL;
  if (mpi_file_rank == 0) counts = new int[mpi_file_size];
  MPI_Gather(&count,1,MPI_INT,counts,1,MPI_INT,0,mpi_file_comm);
  MPI_Offset *offsets = NULL;
  if (mpi_file_rank == 0) offsets = new MPI_Offset[mpi_file_size];
  MPI_Gather(&offset,1,MPI_OFFSET,offsets,1,MPI_OFFSET,0,mpi_file_comm);

  if (mpi_file_rank == 0) {
    assert( fp != NULL); 
    int max_count = 0;
    for (int rank = 0; rank < mpi_file_size; ++rank) 
      max_count = max(max_count,counts[rank]);
    size_t max_nbytes = max_count*sizeOfMpiDatatype(datatype);
    void* tmp_buf = (void*)(new char[max_nbytes]);
    for (int rank = 0; rank < mpi_file_size; ++rank) {
      size_t nbytes = counts[rank]*sizeOfMpiDatatype(datatype);
      fseek(fp,offsets[rank],SEEK_SET);
      if (rank == 0) {
        fwrite(buf,sizeof(char),nbytes,fp);
      }
      else {
        MPI_Recv(tmp_buf,counts[rank],datatype,rank,rank,mpi_file_comm,MPI_STATUS_IGNORE);
        fwrite(tmp_buf,sizeof(char),nbytes,fp);
      }
    }
    delete[] counts;
    delete[] offsets;
    delete[] (char*)(tmp_buf);
  }
  else {
    MPI_Send(buf,count,datatype,0,mpi_file_rank,mpi_file_comm);
  }

  return MPI_SUCCESS;

} 

int MPI_File_close(MPI_File * fh) { 
  using namespace MpiStuff;

  //if ( mpi_rank == 0) 
  //  cout << " intercepted a close call... " << endl;

  FILE * fp; 
  map<MPI_File,pair<FILE*,MPI_Comm> >::iterator it = file_map.find(*fh);
  assert( it != file_map.end());

  fp = it->second.first;

  if ( fp != NULL) { 
    fclose(fp);
  }

  file_map.erase(it);
  if (file_map.empty())
    ifile_desc = 0;

  return MPI_SUCCESS;

} 

int MPI_File_delete(const char *filename, MPI_Info info) {
  using namespace MpiStuff;

  //if (mpi_rank == 0)
  //  cout << " intercepted a delete call..." << endl;

  if (mpi_rank == 0) 
    remove(filename);

  return MPI_SUCCESS;

}

int MPI_File_get_size(MPI_File fh,MPI_Offset *size) {
  using namespace MpiStuff;

  //if ( mpi_rank == 0) 
  //  cout << " intercepted a get size call... " << endl;

  FILE * fp; 
  map<MPI_File,pair<FILE*,MPI_Comm> >::iterator it = file_map.find(fh);
  assert( it != file_map.end());

  fp = it->second.first;

  MPI_Comm mpi_file_comm = it->second.second;
  int mpi_file_rank,mpi_file_size;
  MPI_Comm_rank(mpi_file_comm, &mpi_file_rank);
  MPI_Comm_size(mpi_file_comm, &mpi_file_size);

  if (mpi_file_rank == 0) {
    assert( fp != NULL); 
    fseek(fp, 0, SEEK_END); 
    *size = ftell(fp); 
    fseek(fp, 0, SEEK_SET); 
  }
  MPI_Bcast(size,1,MPI_OFFSET,0,mpi_file_comm);

  return MPI_SUCCESS;

}

int MPI_File_set_size(MPI_File fh,MPI_Offset size) {
  using namespace MpiStuff;

  //if (mpi_rank == 0)
  //  cout << " intercepted a set size call..." << endl;

  return MPI_SUCCESS;
}

#endif
