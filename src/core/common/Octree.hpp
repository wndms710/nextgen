#ifndef _OCTREE_HPP_
#define _OCTREE_HPP_

class Octree {

private:
  
  double bbmin_global[3];
  double bbmax_global[3];
  double (*x_cv)[3];
  int * leaf_buf;
  size_t leaf_buf_size;
  bool b_copy;
  
public:

  Octree() {
    FOR_I3 {
      bbmin_global[i] = HUGE_VAL;
      bbmax_global[i] = -HUGE_VAL;
    }
    x_cv = NULL;
    leaf_buf = NULL;
    leaf_buf_size = 0;
    b_copy = false;
  }

  ~Octree() {
    if (mpi_rank == 0)
      cout << " < cleanup octree..." << endl;
    if (b_copy) 
      DELETE(x_cv); // see below
    DELETE(leaf_buf);
  }

  void getBbox(double bbmin[3],double bbmax[3]) {
    FOR_I3 {
      bbmin[i] = bbmin_global[i];
      bbmax[i] = bbmax_global[i];
    }
  }
  
  void build(double (*x_cv)[3],const int ncv,const bool b_copy = false) {
    
    if (b_copy) {
      // the moving solver dynamically resizes x_vv for ghosts, so you need a copy
      this->x_cv = new double[ncv][3];
      FOR_ICV FOR_I3 this->x_cv[icv][i] = x_cv[icv][i];
      this->b_copy = true;
    }
    else {
      // just get pointer to data
      this->x_cv = x_cv;
    }

    // compute the bbox...
    FOR_I3 {
      assert(bbmin_global[i] == HUGE_VAL);
      assert(bbmax_global[i] == -HUGE_VAL);
    }
    for (int icv = 0; icv < ncv; ++icv) {
      FOR_I3 bbmin_global[i] = min(bbmin_global[i],x_cv[icv][i]);
      FOR_I3 bbmax_global[i] = max(bbmax_global[i],x_cv[icv][i]);
    }
    
    if (mpi_rank == 0)
      cout << " < build octree..." << endl;
    //cout << "got points bbox: " << COUT_VEC(bbmin_global) << " " << COUT_VEC(bbmax_global) << endl;
    
    // and build...
    assert(leaf_buf == NULL);
    //leaf_buf_size = ncv*5+2; // +2 is because ncv could be zero, and we use leaf_buf[0,1] always
    leaf_buf_size = ncv*7+2; // moving rotor-stator case needed this
    leaf_buf = new int[leaf_buf_size];
    memset (leaf_buf,0,sizeof(int)*leaf_buf_size);
    leaf_buf[0] = -8-1;
    leaf_buf[1] = 16; // next leaf location
    
    double bbmin[3],bbmax[3];
    FOR_I3 bbmin[i] = bbmin_global[i]; 
    FOR_I3 bbmax[i] = bbmax_global[i]; 

    // leaf_buf is divided into sets of 8 ints representing the 
    // 2x2x2 octants of the leaf. Values are:
    //  0: empty leaf octant (ready to take a value if one lands here)
    // -1: locked: someone is working on me or my leaf
    // +1,+2,+3...: icv+1
    // -2,-3,-4: -2 index into the child  
    
    int icv = 0;
    int stride = 1;
    int leaf = 0;
    while (icv < ncv) {

      // we use a flag to indicate when a thread needs to add a child leaf...

      bool flag = false; 

      // we should be able to safely perform this traversal of the top of the 
      // tree because once a value is set to a value less than -1, it represents
      // an -1=indexed offset into leaf_buf, and wil never change...
    
      assert(leaf >= 0);
      while (leaf_buf[leaf] < -1) {
        leaf = -leaf_buf[leaf]-1; // update leaf to the next child location: anything less than -1 should never change
        assert((leaf > 0)&&(leaf < leaf_buf_size));
        FOR_I3 {
          const double mid = 0.5*(bbmin[i]+bbmax[i]);
          if (x_cv[icv][i] >= mid) {
            leaf += (1<<i);
            bbmin[i] = mid;
          }
          else {
            bbmax[i] = mid;
          }
        }
        assert((leaf > 0)&&(leaf < leaf_buf_size));
      }
      
      const int value = leaf_buf[leaf];
      if (value > -1) { // skip if locked, AND skip if it happens to have just changed
        if (value == 0) {
          // 0 means empty, so just try and switch to our icv value (1-indexed)...
          if (atomicCAS(leaf_buf+leaf,0,icv+1) == 0) {
            // this leaf was empty and now contains us (as icv+1). so move on...
            FOR_I3 bbmin[i] = bbmin_global[i]; 
            FOR_I3 bbmax[i] = bbmax_global[i]; 
            icv += stride;
            leaf = 0;
          }
        }
        else if (atomicCAS(leaf_buf+leaf,value,-1) == value) {
          // this path means we found an icv+1 (== value) in the octant already
          // we have successfully locked this leaf. Do the work below...
          assert(value >= 1);
          flag = true;
        }
      }

      if (flag) {

        // Put the existing "value" into the correct
        // bin of the child. The new "icv" is then NOT incremented and gets handled 
        // on the next time through. This saves us having to write a special nested/recursive
        // routine here to ensure the the existing and new values go into different child bins...

        // if flag is set, value should contain a +1-indexed icv that has been previously inserted...
      
        int icv0 = value-1; assert((icv0 >= 0)&&(icv0 < ncv));
        int leaf0 = atomicAdd(leaf_buf+1,8); // get the offset to the next free leaf...
        int unlock = -leaf0-1; // this is the value that needs to go at leaf_buf[leaf]...
        FOR_I3 {
          const double mid0 = 0.5*(bbmin[i]+bbmax[i]);
          if (x_cv[icv0][i] >= mid0) {
            leaf0 += (1<<i); // use comparisons to increment the leaf0 to the correct octant
          }
        }
        assert((leaf0 > 0)&&(leaf0 < leaf_buf_size));
        
        // put the previous value back in the child at the right spot. Note that we are not advancing
        // icv; it will get handled the next time through (potentially, unless it hits the child 
        // again because they are in the same bin. Note that perfectly coincident points will fail
        // because they will never be in different bins at any child level)..
        //int old = atomicExch(leaf_buf+leaf0,value);
        //assert(old == 0); // should have been "empty"
      
        // and release the lock, leaving the proper -1-indexed offset to the 
        // first leaf of the child...
        //old = atomicExch(leaf_buf+leaf,unlock);
        //assert(old == -1); // should have been a lock

        assert(leaf_buf[leaf0] == 0);
        leaf_buf[leaf0] = value;

        //__threadfence(); // make sure value is flushed out of cache before unlocking

        // do not increment icv, just continue, because we may have to split again

        assert(leaf_buf[leaf] == -1);
        leaf_buf[leaf] = unlock; // unlock with the new -1 indexed offset to the new child
        
      }
      
    }
    
  }
  
  void buildListForSphere(vector<int>& nbrVec,const double x[3],const double r) {
    
    // segfaults w/o this when ncv == 0
    if (leaf_buf_size == 2)
      return;

    const int nstack_max = 64; // should be big enough
    int * stack_int = new int[nstack_max];
    double * stack_double = new double[nstack_max*6];
  
    int nstack = 1;
    stack_int[0] = 0;
    FOR_I3 stack_double[i  ] = bbmin_global[i];
    FOR_I3 stack_double[i+3] = bbmax_global[i];
    
    //int nstack_max_actual = 1;
    const double r2 = r*r;
    int leaf;
    double bbox[9];
    while (nstack > 0) {
      // pop the top off the stack...
      --nstack;
      leaf = stack_int[nstack];
      FOR_I3 bbox[i]   = stack_double[nstack*6+i];
      FOR_I3 bbox[i+6] = stack_double[nstack*6+3+i];
      FOR_I3 bbox[i+3] = 0.5*(bbox[i]+bbox[i+6]);
      assert(leaf_buf[leaf] < -1);
      leaf = -leaf_buf[leaf]-1;
      // child leaves are in leaf+0...leaf+7, with possible values of
      // -1: empty
      // -2,-3,-4...: -1-indexed link to other children
      // 0,1,2,3...: value to be returned
      // 000...
      if (leaf_buf[leaf] > 0) {
        // terminal leaf with +1 indexed value...
        const int icv = leaf_buf[leaf]-1;
        if (DIST2(x,x_cv[icv]) <= r2)
          nbrVec.push_back(icv);
      }
      else if ((leaf_buf[leaf] < -1)&&(calcMinD2PointToBbox(x,bbox[0],bbox[1],bbox[2],bbox[3],bbox[4],bbox[5]) <= r2)) { 
        stack_int[nstack] = leaf;
        stack_double[nstack*6  ] = bbox[0];
        stack_double[nstack*6+1] = bbox[1];
        stack_double[nstack*6+2] = bbox[2];
        stack_double[nstack*6+3] = bbox[3];
        stack_double[nstack*6+4] = bbox[4];
        stack_double[nstack*6+5] = bbox[5];
        ++nstack;
        //nstack_max_actual = max(nstack_max_actual,nstack);
      }
      // 001
      if (leaf_buf[leaf+1] > 0) {
        // terminal leaf with value...
        const int icv = leaf_buf[leaf+1]-1;
        if (DIST2(x,x_cv[icv]) <= r2)
          nbrVec.push_back(icv);
      }
      else if ((leaf_buf[leaf+1] < -1)&&(calcMinD2PointToBbox(x,bbox[0+3],bbox[1],bbox[2],bbox[3+3],bbox[4],bbox[5]) <= r2)) { 
        assert(nstack < nstack_max);
        stack_int[nstack] = leaf+1;
        stack_double[nstack*6  ] = bbox[0+3];
        stack_double[nstack*6+1] = bbox[1];
        stack_double[nstack*6+2] = bbox[2];
        stack_double[nstack*6+3] = bbox[3+3];
        stack_double[nstack*6+4] = bbox[4];
        stack_double[nstack*6+5] = bbox[5];
        ++nstack;
        //nstack_max_actual = max(nstack_max_actual,nstack);
      }
      // 010...
      if (leaf_buf[leaf+2] > 0) {
        // terminal leaf with value...
        const int icv = leaf_buf[leaf+2]-1;
        if (DIST2(x,x_cv[icv]) <= r2)
          nbrVec.push_back(icv);
      }
      else if ((leaf_buf[leaf+2] < -1)&&(calcMinD2PointToBbox(x,bbox[0],bbox[1+3],bbox[2],bbox[3],bbox[4+3],bbox[5]) <= r2)) { 
        assert(nstack < nstack_max);
        stack_int[nstack] = leaf+2;
        stack_double[nstack*6  ] = bbox[0];
        stack_double[nstack*6+1] = bbox[1+3];
        stack_double[nstack*6+2] = bbox[2];
        stack_double[nstack*6+3] = bbox[3];
        stack_double[nstack*6+4] = bbox[4+3];
        stack_double[nstack*6+5] = bbox[5];
        ++nstack;
        //nstack_max_actual = max(nstack_max_actual,nstack);
      }
      // 011
      if (leaf_buf[leaf+3] > 0) {
        // terminal leaf with value...
        const int icv = leaf_buf[leaf+3]-1;
        if (DIST2(x,x_cv[icv]) <= r2)
          nbrVec.push_back(icv);
      }
      else if ((leaf_buf[leaf+3] < -1)&&(calcMinD2PointToBbox(x,bbox[0+3],bbox[1+3],bbox[2],bbox[3+3],bbox[4+3],bbox[5]) <= r2)) { 
        assert(nstack < nstack_max);
        stack_int[nstack] = leaf+3;
        stack_double[nstack*6  ] = bbox[0+3];
        stack_double[nstack*6+1] = bbox[1+3];
        stack_double[nstack*6+2] = bbox[2];
        stack_double[nstack*6+3] = bbox[3+3];
        stack_double[nstack*6+4] = bbox[4+3];
        stack_double[nstack*6+5] = bbox[5];
        ++nstack;
        //nstack_max_actual = max(nstack_max_actual,nstack);
      }
      // 100...
      if (leaf_buf[leaf+4] > 0) {
        // terminal leaf with value...
        const int icv = leaf_buf[leaf+4]-1;
        if (DIST2(x,x_cv[icv]) <= r2)
          nbrVec.push_back(icv);
      }
      else if ((leaf_buf[leaf+4] < -1)&&(calcMinD2PointToBbox(x,bbox[0],bbox[1],bbox[2+3],bbox[3],bbox[4],bbox[5+3]) <= r2)) { 
        assert(nstack < nstack_max);
        stack_int[nstack] = leaf+4;
        stack_double[nstack*6  ] = bbox[0];
        stack_double[nstack*6+1] = bbox[1];
        stack_double[nstack*6+2] = bbox[2+3];
        stack_double[nstack*6+3] = bbox[3];
        stack_double[nstack*6+4] = bbox[4];
        stack_double[nstack*6+5] = bbox[5+3];
        ++nstack;
        //nstack_max_actual = max(nstack_max_actual,nstack);
      }
      // 101
      if (leaf_buf[leaf+5] > 0) {
        // terminal leaf with value...
        const int icv = leaf_buf[leaf+5]-1;
        if (DIST2(x,x_cv[icv]) <= r2)
          nbrVec.push_back(icv);
      }
      else if ((leaf_buf[leaf+5] < -1)&&(calcMinD2PointToBbox(x,bbox[0+3],bbox[1],bbox[2+3],bbox[3+3],bbox[4],bbox[5+3]) <= r2)) { 
        assert(nstack < nstack_max);
        stack_int[nstack] = leaf+5;
        stack_double[nstack*6  ] = bbox[0+3];
        stack_double[nstack*6+1] = bbox[1];
        stack_double[nstack*6+2] = bbox[2+3];
        stack_double[nstack*6+3] = bbox[3+3];
        stack_double[nstack*6+4] = bbox[4];
        stack_double[nstack*6+5] = bbox[5+3];
        ++nstack;
        //nstack_max_actual = max(nstack_max_actual,nstack);
      }
      // 110...
      if (leaf_buf[leaf+6] > 0) {
        // terminal leaf with value...
        const int icv = leaf_buf[leaf+6]-1;
        if (DIST2(x,x_cv[icv]) <= r2)
          nbrVec.push_back(icv);
      }
      else if ((leaf_buf[leaf+6] < -1)&&(calcMinD2PointToBbox(x,bbox[0],bbox[1+3],bbox[2+3],bbox[3],bbox[4+3],bbox[5+3]) <= r2)) { 
        assert(nstack < nstack_max);
        stack_int[nstack] = leaf+6;
        stack_double[nstack*6  ] = bbox[0];
        stack_double[nstack*6+1] = bbox[1+3];
        stack_double[nstack*6+2] = bbox[2+3];
        stack_double[nstack*6+3] = bbox[3];
        stack_double[nstack*6+4] = bbox[4+3];
        stack_double[nstack*6+5] = bbox[5+3];
        ++nstack;
        //nstack_max_actual = max(nstack_max_actual,nstack);
      }
      // 011
      if (leaf_buf[leaf+7] > 0) {
        // terminal leaf with value...
        const int icv = leaf_buf[leaf+7]-1;
        if (DIST2(x,x_cv[icv]) <= r2)
          nbrVec.push_back(icv);
      }
      else if ((leaf_buf[leaf+7] < -1)&&(calcMinD2PointToBbox(x,bbox[0+3],bbox[1+3],bbox[2+3],bbox[3+3],bbox[4+3],bbox[5+3]) <= r2)) { 
        assert(nstack < nstack_max);
        stack_int[nstack] = leaf+7;
        stack_double[nstack*6  ] = bbox[0+3];
        stack_double[nstack*6+1] = bbox[1+3];
        stack_double[nstack*6+2] = bbox[2+3];
        stack_double[nstack*6+3] = bbox[3+3];
        stack_double[nstack*6+4] = bbox[4+3];
        stack_double[nstack*6+5] = bbox[5+3];
        ++nstack;
        //nstack_max_actual = max(nstack_max_actual,nstack);
      }
    }
  
    //cout << " > got nstack_max_actual: " << nstack_max_actual << endl;
    
    delete[] stack_int;
    delete[] stack_double;
    
  }
  
  // supporting functions: put these in MiscUtils at global scope eventually... 

  int atomicCAS(int * buf,const int compare,const int val) {

    // compare and swap...

    if (*buf == compare) {
      *buf = val;
      return compare;
    }
    else {
      return *buf;
    }

  }

  int atomicAdd(int * buf,const int val) {
  
    const int old = *buf;
    *buf += val;
    return old;

  }

  double calcMinD2PointToBbox(const double x[3],
                              const double bbmin0,const double bbmin1,const double bbmin2,
                              const double bbmax0,const double bbmax1,const double bbmax2) {
  
    // unrolled version of above...
  
    double d2;

    // x...
    if (x[0] < bbmin0) {
      d2 = (bbmin0 - x[0])*(bbmin0 - x[0]);
    }
    else if (x[0] > bbmax0) {
      d2 = (x[0] - bbmax0)*(x[0] - bbmax0);
    }
    else {
      d2 = 0.0;
    }

    // y...
    if (x[1] < bbmin1) {
      d2 += (bbmin1 - x[1])*(bbmin1 - x[1]);
    }
    else if (x[1] > bbmax1) {
      d2 += (x[1] - bbmax1)*(x[1] - bbmax1);
    }
  
    // z...
    if (x[2] < bbmin2) {
      d2 += (bbmin2 - x[2])*(bbmin2 - x[2]);
    }
    else if (x[2] > bbmax2) {
      d2 += (x[2] - bbmax2)*(x[2] - bbmax2);
    }

    return(d2);

  }
  
};

#ifdef TEST_OCTREE

#include "Adt.hpp"

void testOctree() {
  
  const int ncv = getIntParam("N",10000);
  double (*x_cv)[3] = new double[ncv][3];
  
  for (int icv = 0; icv < ncv; ++icv) {
    FOR_I3 x_cv[icv][i] = double(rand())/double(RAND_MAX)-0.5 + double(i);
  }

  double wtime0 = MPI_Wtime();
  
  Octree oct;
  oct.build(x_cv,ncv);
  
  double wtime1 = MPI_Wtime();
  cout << "octree build time: " << (wtime1-wtime0) << " seconds, " << (wtime1-wtime0)*1E6/double(ncv) << " us/cv " << endl;

  // now try to find the points within delta = 2.0*pow(1.0/double(n),1.0/3.0);
  const double r = 2.0*pow(1.0/double(ncv),1.0/3.0);
  const double r2 = r*r;
  
  wtime0 = MPI_Wtime();
  
  vector<int> nbrVec;
  for (int icv = 0; icv < ncv; ++icv) {

    assert(nbrVec.empty());
    oct.buildListForSphere(nbrVec,x_cv[icv],r);

    /*
    // confirm we found ourselves...
    bool found = false;
    for (int ii = 0; ii < nbrVec.size(); ++ii) {
      const int icv_nbr = nbrVec[ii];
      if (icv_nbr == icv) {
        found = true;
        break;
      }
    }
    assert(found);
    
    // now do hard loop through ALL cv's and find the ones in the sphere
    // and confirm they are in the int vect returned by the octree...
    int count = 0;
    for (int icv_nbr = 0; icv_nbr < ncv; ++icv_nbr) {
      const double dist2 = DIST2(x_cv[icv],x_cv[icv_nbr]);
      if (dist2 <= r2) {
        // make sure icv_nbr is in the nbrVec...
        int ii;
        for (ii = 0; ii < nbrVec.size(); ++ii)
          if (nbrVec[ii] == icv_nbr)
            break;
        assert(ii < nbrVec.size());
        ++count;
      }
    }
    
    cout << "icv " << icv << " out of " << ncv << " looks good with count: " << count << " out of nbrVec.size(): " << nbrVec.size() << endl;
    assert(count == nbrVec.size());
    */

    nbrVec.clear();
    
  }
  
  wtime1 = MPI_Wtime();
  cout << "octree lookup time: " << (wtime1-wtime0) << " seconds, " << (wtime1-wtime0)*1E6/double(ncv) << " us/cv " << endl;

  // try the adt...

  wtime0 = MPI_Wtime();
  
  Adt<double> adt(ncv,x_cv,x_cv);

  wtime1 = MPI_Wtime();
  cout << "adt build time: " << (wtime1-wtime0) << " seconds, " << (wtime1-wtime0)*1E6/double(ncv) << " us/cv " << endl;

  wtime0 = MPI_Wtime();

  for (int icv = 0; icv < ncv; ++icv) {
  
    assert(nbrVec.empty());
    adt.buildListForSphere(nbrVec,x_cv[icv],r);

    /*
    // confirm we found ourselves...
    bool found = false;
    for (int ii = 0; ii < nbrVec.size(); ++ii) {
      const int icv_nbr = nbrVec[ii];
      if (icv_nbr == icv) {
        found = true;
        break;
      }
    }
    assert(found);
    
    // now do hard loop through ALL cv's and find the ones in the sphere
    // and confirm they are in the int vect returned by the octree...
    int count = 0;
    for (int icv_nbr = 0; icv_nbr < ncv; ++icv_nbr) {
      const double dist2 = DIST2(x_cv[icv],x_cv[icv_nbr]);
      if (dist2 <= r2) {
        // make sure icv_nbr is in the nbrVec...
        int ii;
        for (ii = 0; ii < nbrVec.size(); ++ii)
          if (nbrVec[ii] == icv_nbr)
            break;
        assert(ii < nbrVec.size());
        ++count;
      }
    }
    
    cout << "icv " << icv << " out of " << ncv << " looks good with count: " << count << " out of nbrVec.size(): " << nbrVec.size() << endl;
    assert(count == nbrVec.size());

    */

    nbrVec.clear();
    
  }
  
  wtime1 = MPI_Wtime();
  cout << "adt lookup time: " << (wtime1-wtime0) << " seconds, " << (wtime1-wtime0)*1E6/double(ncv) << " us/cv " << endl;

  MPI_Pause("done");
  
}

#endif

#endif
